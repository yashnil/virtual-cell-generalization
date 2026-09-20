"""Zero-shot recoverability diagnostic (v1).

Leave-one-context-out study of how much held-out perturbation response simple,
transparent, source-only baselines recover — and how much reproducible
context-specific signal is left over.

Diagnostic only. **No learned model of any kind is built here.** Every baseline
is a closed-form function of source-context responses and basal (control)
expression. The held-out context's perturbation responses are used only to
evaluate.

Usage::

    uv run python scripts/run_zero_shot_recoverability.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from virtual_cell.analysis import loco, robustness
from virtual_cell.data import scperteval
from virtual_cell.decomposition import anova

REPO_ROOT = Path(__file__).resolve().parents[1]
DESIGN_DIR = REPO_ROOT / "data" / "splits" / "four_context_v1"
CANONICAL_DIR = REPO_ROOT / "data" / "processed" / "four_context_v1"
SPLITS_DIR = REPO_ROOT / "data" / "splits" / "loco_v1"
WORK_DIR = REPO_ROOT / "data" / "processed" / "zero_shot_v1"
OUT_DIR = REPO_ROOT / "outputs" / "zero_shot_v1"

SCRIPT_VERSION = "1.0.0"


def sha256_text(t: str) -> str:
    return hashlib.sha256(t.encode()).hexdigest()


# --------------------------------------------------------------------------
# 1. freeze the LOCO splits
# --------------------------------------------------------------------------


def freeze_splits(perts: list[str], genes: list[str]) -> list[loco.LocoFold]:
    print("=" * 70)
    print("1. FREEZE LEAVE-ONE-CONTEXT-OUT SPLITS")
    print("=" * 70)
    folds = loco.make_folds(scperteval.CONTEXTS)
    SPLITS_DIR.mkdir(parents=True, exist_ok=True)
    files = {}
    for f in folds:
        text = f"# leave-one-context-out fold\ntarget\t{f.target}\n" + "".join(
            f"source\t{s}\n" for s in f.sources
        )
        name = f"fold_{f.target}.tsv"
        (SPLITS_DIR / name).write_text(text)
        files[name] = sha256_text(text)
        print(
            f"  {'+'.join(scperteval.dataset(s).cell_line for s in f.sources):>28s}"
            f"  ->  {scperteval.dataset(f.target).cell_line}"
        )

    manifest = {
        "split_version": "loco_v1",
        "generated": date.today().isoformat(),
        "generation_script": "scripts/run_zero_shot_recoverability.py",
        "script_version": SCRIPT_VERSION,
        "design": "four_context_v1",
        "n_perturbations": len(perts),
        "n_genes": len(genes),
        "folds": [{"target": f.target, "sources": list(f.sources)} for f in folds],
        "allowed_inputs": [
            "basal/control cells of the held-out context",
            "perturbation identity",
            "perturbation responses of the three source contexts",
            "external priors free of target perturbation outcomes",
        ],
        "forbidden_inputs": [
            "any perturbation-response measurement from the held-out context",
            "feature selection using held-out perturbation responses",
            "normalisation fitted to held-out perturbation responses",
            "hyperparameter selection using held-out perturbation outcomes",
            "four-context beta or gamma (both contain the held-out response)",
        ],
        "evaluation_only_quantities": [
            "held-out delta",
            "four-context gamma",
            "target response template (mean over perturbations)",
            "target split halves",
        ],
        "file_sha256": files,
    }
    (SPLITS_DIR / "loco_manifest.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False))
    print(f"\n  frozen to {SPLITS_DIR}")
    return folds


# --------------------------------------------------------------------------
# 2. target split halves (evaluation-only reliability machinery)
# --------------------------------------------------------------------------


def build_half_deltas(
    data_dir: Path, perts: list[str], genes: list[str], *, n_repeats: int, seed: int
) -> Path:
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    path = WORK_DIR / "half_deltas.npy"
    if path.exists():
        print(f"\n  reusing cached half deltas at {path.name}")
        return path
    n_c, n_p, n_g = len(scperteval.CONTEXTS), len(perts), len(genes)
    arr = np.lib.format.open_memmap(
        path, mode="w+", dtype=np.float32, shape=(n_repeats, 2, n_c, n_p, n_g)
    )
    print(f"\n  building half deltas: {arr.nbytes / 1e9:.1f} GB, {n_repeats} repeats")
    for ci, ds in enumerate(scperteval.DATASETS):
        t0 = time.time()
        X, group, X_ctrl = robustness.load_cells(
            ds.path(data_dir), genes=genes, perturbations=perts
        )
        for r in range(n_repeats):
            a, b = robustness.half_deltas(
                X,
                group,
                X_ctrl,
                n_perturbations=n_p,
                rng=np.random.default_rng([seed, ci, r]),
                control_scheme="shared",
                aggregation="mean_log",
            )
            arr[r, 0, ci] = a.astype(np.float32)
            arr[r, 1, ci] = b.astype(np.float32)
        del X, X_ctrl
        arr.flush()
        print(f"    {ds.name:16s} [{time.time() - t0:.0f}s]", flush=True)
    del arr
    return path


# --------------------------------------------------------------------------
# baselines
# --------------------------------------------------------------------------


def build_predictions(D: np.ndarray, control: np.ndarray, fold: loco.LocoFold) -> dict:
    """Every zero-shot baseline for one fold. Source responses + basal only."""
    src = fold.source_indices
    preds = {"source_mean": loco.source_mean(D, src)}
    for s in src:
        preds[f"single_{scperteval.dataset(scperteval.CONTEXTS[s]).cell_line}"] = (
            loco.single_source(D, s)
        )
    nearest, sim = loco.nearest_context(control, fold.target_index, src)
    preds["nearest_basal"] = loco.single_source(D, nearest)
    w_affine = loco.basal_affine_weights(control, fold.target_index, src)
    preds["basal_affine"] = loco.weighted_source(D, src, w_affine)
    w_simplex = loco.basal_simplex_weights(control, fold.target_index, src)
    preds["basal_simplex"] = loco.weighted_source(D, src, w_simplex)
    preds["null_zero"] = np.zeros_like(preds["source_mean"])
    meta = {
        "nearest_source": scperteval.CONTEXTS[nearest],
        "basal_similarity": {scperteval.CONTEXTS[s]: float(sim[s]) for s in src},
        "affine_weights": {
            scperteval.CONTEXTS[s]: float(w) for s, w in zip(src, w_affine, strict=True)
        },
        "simplex_weights": {
            scperteval.CONTEXTS[s]: float(w) for s, w in zip(src, w_simplex, strict=True)
        },
    }
    return preds, meta


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", type=Path, default=REPO_ROOT / scperteval.DATA_SUBDIR)
    ap.add_argument("--n-repeats", type=int, default=10)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--n-boot", type=int, default=2000)
    args = ap.parse_args()

    t_start = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    perts = (DESIGN_DIR / "shared_perturbations.txt").read_text().split()
    genes = (DESIGN_DIR / "shared_genes.txt").read_text().split()
    D = np.load(CANONICAL_DIR / "delta_tensor.npy").astype(np.float64)
    control = np.load(CANONICAL_DIR / "control_means.npy").astype(np.float64)
    counts = np.load(CANONICAL_DIR / "cell_counts.npy")
    folds = freeze_splits(perts, genes)

    # evaluation-only four-context decomposition
    dec = anova.decompose(D, cell_lines=scperteval.CONTEXTS, perturbations=perts)

    half_path = build_half_deltas(
        args.data_dir, perts, genes, n_repeats=args.n_repeats, seed=args.seed
    )
    halves = np.load(half_path, mmap_mode="r")

    sim_matrix = loco.basal_similarity(control)
    pd.DataFrame(
        sim_matrix, index=list(scperteval.CONTEXTS), columns=list(scperteval.CONTEXTS)
    ).to_csv(OUT_DIR / "basal_similarity.csv")
    print("\n  basal similarity (control profiles, Pearson):")
    print(
        pd.DataFrame(
            sim_matrix,
            index=[scperteval.dataset(c).cell_line for c in scperteval.CONTEXTS],
            columns=[scperteval.dataset(c).cell_line for c in scperteval.CONTEXTS],
        ).round(4)
    )

    metric_rows, reliability_rows, gamma_rows, fold_meta = [], [], [], {}

    print("\n" + "=" * 70)
    print("2/5. SOURCE-ONLY TRANSFER AND SIMPLE ZERO-SHOT BASELINES")
    print("=" * 70)

    for fold in folds:
        c = fold.target_index
        cl = scperteval.dataset(fold.target).cell_line
        Y = D[c]
        T_eval = loco.target_template_evaluation_only(Y)
        preds, meta = build_predictions(D, control, fold)
        fold_meta[fold.target] = meta
        print(
            f"\n--- held out {cl} "
            f"(nearest basal source: {scperteval.dataset(meta['nearest_source']).cell_line}) ---"
        )

        for name, P in preds.items():
            m = loco.per_perturbation_metrics(Y, P)
            m_tr = loco.per_perturbation_metrics(
                loco.template_removed(Y, T_eval),
                loco.template_removed(P, T_eval),
                with_spearman=False,
            )
            m["pearson_template_removed"] = m_tr["pearson"]
            m["energy_explained_template_removed"] = m_tr["energy_explained"]
            m["fold"] = fold.target
            m["cell_line"] = cl
            m["baseline"] = name
            m["perturbation"] = perts
            m["target_cells"] = counts[c]
            metric_rows.append(m)
            b = loco.bootstrap_median(m["pearson"].to_numpy(), n_boot=args.n_boot, seed=1)
            bt = loco.bootstrap_median(
                m["pearson_template_removed"].to_numpy(), n_boot=args.n_boot, seed=1
            )
            print(
                f"  {name:22s} r={b['median']:+.4f} [{b['lo']:+.4f},{b['hi']:+.4f}]  "
                f"r_templ-rm={bt['median']:+.4f}  "
                f"energy={np.nanmedian(m['energy_explained']):+.4f}  "
                f"cos={np.nanmedian(m['cosine']):+.4f}"
            )

        # --- 4. reliability-normalised evaluation, averaged over repeats ---
        for name, P in preds.items():
            if name == "null_zero":
                continue
            acc = []
            for r in range(halves.shape[0]):
                acc.append(
                    loco.reliability_normalised_evaluation(
                        P,
                        np.asarray(halves[r, 0, c], dtype=np.float64),
                        np.asarray(halves[r, 1, c], dtype=np.float64),
                    )
                )
            mean_tab = sum(acc) / len(acc)
            mean_tab["fold"] = fold.target
            mean_tab["cell_line"] = cl
            mean_tab["baseline"] = name
            mean_tab["perturbation"] = perts
            reliability_rows.append(mean_tab)

        # --- 6. gamma recoverability ---------------------------------------
        A = preds["source_mean"]
        Ac = A - A.mean(axis=0, keepdims=True)
        gamma_true = dec.gamma[c]  # evaluation only
        residual = Y - A
        residual_centred = residual - residual.mean(axis=0, keepdims=True)

        # gamma reliability.
        #
        # gamma[c,p] is a function of ALL FOUR contexts, so replacing only the
        # target row with a half would leave the three source contexts
        # contributing an identical component to both "halves" and inflate the
        # correlation badly. Every context is therefore split, gamma is
        # decomposed separately per repeat, and the per-repeat correlations are
        # averaged -- averaging gamma across repeats first would average the
        # noise away and inflate the estimate a second time.
        per_repeat = []
        for r in range(halves.shape[0]):
            g1 = anova.decompose(np.asarray(halves[r, 0], dtype=np.float64)).gamma[c]
            g2 = anova.decompose(np.asarray(halves[r, 1], dtype=np.float64)).gamma[c]
            per_repeat.append(loco.split_half_reliability(g1, g2))
        gamma_rho = np.nanmean(np.vstack(per_repeat), axis=0)
        gamma_rho_full = np.asarray(loco.spearman_brown(gamma_rho))

        for name, P in preds.items():
            if name in ("null_zero", "source_mean"):
                continue
            Pc = P - P.mean(axis=0, keepdims=True)
            gamma_pred = Pc - Ac
            r_g = np.array(
                [loco._safe_pearson(gamma_true[i], gamma_pred[i]) for i in range(len(perts))]
            )
            gamma_rows.append(
                pd.DataFrame(
                    {
                        "fold": fold.target,
                        "cell_line": cl,
                        "baseline": name,
                        "perturbation": perts,
                        "r_gamma": r_g,
                        "gamma_rho_half": gamma_rho,
                        "gamma_rho_full": gamma_rho_full,
                        "gamma_ceiling_half": loco.ceiling_from_reliability(gamma_rho),
                        "gamma_ceiling_full": loco.ceiling_from_reliability(gamma_rho_full),
                        "r_gamma_disattenuated": loco.disattenuate(r_g, gamma_rho_full),
                        "gamma_norm": np.linalg.norm(gamma_true, axis=1),
                        "residual_norm": np.linalg.norm(residual_centred, axis=1),
                        "gamma_pred_norm": np.linalg.norm(gamma_pred, axis=1),
                    }
                )
            )

        bg = loco.bootstrap_median(gamma_rho, n_boot=args.n_boot, seed=2)
        print(
            f"  [gamma] reliability median={bg['median']:.4f} "
            f"[{bg['lo']:.4f},{bg['hi']:.4f}]  ceiling={np.sqrt(max(bg['median'], 0)):.4f}"
        )

    metrics = pd.concat(metric_rows, ignore_index=True)
    reliability = pd.concat(reliability_rows, ignore_index=True)
    gammas = pd.concat(gamma_rows, ignore_index=True)
    metrics.to_csv(OUT_DIR / "loco_metrics.csv", index=False)
    reliability.to_csv(OUT_DIR / "loco_reliability.csv", index=False)
    gammas.to_csv(OUT_DIR / "loco_gamma.csv", index=False)

    # --- cross-context gamma structure ------------------------------------
    print("\n" + "=" * 70)
    print("6. CROSS-CONTEXT GAMMA STRUCTURE")
    print("=" * 70)
    n_c = len(scperteval.CONTEXTS)
    pair_rows = []
    for i in range(n_c):
        for j in range(i + 1, n_c):
            rs = np.array(
                [loco._safe_pearson(dec.gamma[i, p], dec.gamma[j, p]) for p in range(len(perts))]
            )
            pair_rows.append(
                {
                    "context_a": scperteval.CONTEXTS[i],
                    "context_b": scperteval.CONTEXTS[j],
                    "cell_line_a": scperteval.dataset(scperteval.CONTEXTS[i]).cell_line,
                    "cell_line_b": scperteval.dataset(scperteval.CONTEXTS[j]).cell_line,
                    "median_gamma_r": float(np.nanmedian(rs)),
                    "basal_similarity": float(sim_matrix[i, j]),
                }
            )
    pairs = pd.DataFrame(pair_rows)
    pairs["null_expectation"] = -1.0 / (n_c - 1)
    pairs["excess_over_null"] = pairs["median_gamma_r"] - pairs["null_expectation"]
    pairs.to_csv(OUT_DIR / "gamma_cross_context.csv", index=False)
    print(pairs.round(4).to_string(index=False))
    print(
        f"\n  note: sum_c gamma[c,p] = 0 forces a null expectation of "
        f"-1/(C-1) = {-1.0 / (n_c - 1):.4f}, not 0."
    )
    if len(pairs) > 2:
        rho = pairs["basal_similarity"].corr(pairs["excess_over_null"], method="spearman")
        print(
            f"  Spearman(basal similarity, gamma-correlation excess) = {rho:+.4f} "
            f"over {len(pairs)} context pairs"
        )

    # --- 7. transferability descriptives ----------------------------------
    print("\n" + "=" * 70)
    print("7. TRANSFERABILITY DESCRIPTIVES (no target defined)")
    print("=" * 70)
    sm = metrics[metrics.baseline == "source_mean"].reset_index(drop=True)
    rel_sm = reliability[reliability.baseline == "source_mean"].reset_index(drop=True)
    desc = sm[
        [
            "fold",
            "cell_line",
            "perturbation",
            "pearson",
            "pearson_template_removed",
            "energy_explained",
            "obs_norm",
            "target_cells",
        ]
    ].copy()
    desc["target_rho_half"] = rel_sm["rho_half"]
    desc["ceiling"] = rel_sm["ceiling"]
    desc["r_disattenuated"] = rel_sm["r_disattenuated"]
    # source agreement: mean pairwise correlation between the three source responses
    agree = []
    for _, row in desc.iterrows():
        f = [x for x in folds if x.target == row["fold"]][0]
        pi = perts.index(row["perturbation"])
        rs = [
            loco._safe_pearson(D[a, pi], D[b, pi])
            for ii, a in enumerate(f.source_indices)
            for b in f.source_indices[ii + 1 :]
        ]
        agree.append(float(np.nanmean(rs)))
    desc["source_agreement"] = agree
    gsub = gammas[gammas.baseline == "nearest_basal"].reset_index(drop=True)
    desc["gamma_rho_half"] = gsub["gamma_rho_half"]
    desc["gamma_norm"] = gsub["gamma_norm"]
    # target-gene basal expression, where the target gene is in the panel
    gene_index = {g: i for i, g in enumerate(genes)}
    desc["target_gene_basal"] = [
        float(
            control[
                [x for x in folds if x.target == row["fold"]][0].target_index,
                gene_index[row["perturbation"]],
            ]
        )
        if row["perturbation"] in gene_index
        else np.nan
        for _, row in desc.iterrows()
    ]
    desc["basal_context_distance"] = [
        1.0
        - float(
            np.mean(
                [
                    sim_matrix[[x for x in folds if x.target == row["fold"]][0].target_index, s]
                    for s in [x for x in folds if x.target == row["fold"]][0].source_indices
                ]
            )
        )
        for _, row in desc.iterrows()
    ]
    desc.to_csv(OUT_DIR / "transferability_descriptives.csv", index=False)

    corr_targets = [
        "source_agreement",
        "target_rho_half",
        "obs_norm",
        "target_cells",
        "target_gene_basal",
        "gamma_norm",
        "basal_context_distance",
    ]
    print("\n  Spearman vs conserved-transfer success (r, pooled over folds):")
    for col in corr_targets:
        rho = desc["pearson"].corr(desc[col], method="spearman")
        rho_d = desc["r_disattenuated"].corr(desc[col], method="spearman")
        print(f"    {col:26s} raw={rho:+.4f}   disattenuated={rho_d:+.4f}")

    # --- summary ----------------------------------------------------------
    summary = {
        "script_version": SCRIPT_VERSION,
        "generated": date.today().isoformat(),
        "runtime_seconds": time.time() - t_start,
        "n_repeats": args.n_repeats,
        "seed": args.seed,
        "folds": {f.target: fold_meta[f.target] for f in folds},
        "basal_similarity": pd.DataFrame(
            sim_matrix, index=list(scperteval.CONTEXTS), columns=list(scperteval.CONTEXTS)
        ).to_dict(),
        "by_fold_baseline": {},
    }
    for (fold_name, baseline), g in metrics.groupby(["fold", "baseline"]):
        summary["by_fold_baseline"][f"{fold_name}|{baseline}"] = {
            "pearson": loco.bootstrap_median(g["pearson"].to_numpy(), n_boot=args.n_boot, seed=1),
            "pearson_template_removed": loco.bootstrap_median(
                g["pearson_template_removed"].to_numpy(), n_boot=args.n_boot, seed=1
            ),
            "energy_explained": float(np.nanmedian(g["energy_explained"])),
            "cosine": float(np.nanmedian(g["cosine"])),
            "spearman": float(np.nanmedian(g["spearman"])),
            "mse": float(np.nanmedian(g["mse"])),
        }
    for (fold_name, baseline), g in reliability.groupby(["fold", "baseline"]):
        key = f"{fold_name}|{baseline}"
        summary["by_fold_baseline"].setdefault(key, {})
        summary["by_fold_baseline"][key]["reliability"] = {
            "rho_half": float(np.nanmedian(g["rho_half"])),
            "ceiling": float(np.nanmedian(g["ceiling"])),
            "r_half_mean": float(np.nanmedian(g["r_half_mean"])),
            "r_disattenuated": float(np.nanmedian(g["r_disattenuated"])),
        }
    summary["gamma"] = {}
    for (fold_name, baseline), g in gammas.groupby(["fold", "baseline"]):
        summary["gamma"][f"{fold_name}|{baseline}"] = {
            "r_gamma": loco.bootstrap_median(g["r_gamma"].to_numpy(), n_boot=args.n_boot, seed=3),
            "gamma_rho_half": float(np.nanmedian(g["gamma_rho_half"])),
            "gamma_rho_full": float(np.nanmedian(g["gamma_rho_full"])),
            "gamma_ceiling_half": float(np.nanmedian(g["gamma_ceiling_half"])),
            "gamma_ceiling_full": float(np.nanmedian(g["gamma_ceiling_full"])),
            "r_gamma_disattenuated": float(np.nanmedian(g["r_gamma_disattenuated"])),
        }
    summary["gamma_cross_context"] = pairs.to_dict("records")
    summary["transferability_spearman"] = {
        col: {
            "vs_pearson": float(desc["pearson"].corr(desc[col], method="spearman")),
            "vs_disattenuated": float(desc["r_disattenuated"].corr(desc[col], method="spearman")),
        }
        for col in corr_targets
    }
    (OUT_DIR / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"\nTotal runtime: {(time.time() - t_start) / 60:.1f} min")
    print(f"Wrote {OUT_DIR}")


if __name__ == "__main__":
    main()
