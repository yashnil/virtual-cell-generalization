"""Pathway residual model v2 — CLEAN (beta-free) gamma target.

The single permitted follow-up to v1. **Exactly one thing changes**: the
training target is now the unshrunk centred residual, which is beta-free by
construction (``centre_p(Y - A) = gamma - mean_S gamma``). Splits, features,
baseline, leakage contract, families, grids and metrics are all reused verbatim
from v1.

Also reports ``lambda_theory = (3 + s)/4`` as a diagnostic and an oracle sweep as
a post-hoc ceiling only.

No deep model, no M3, no new features, no Arc predictions, no gene-level
lift-back.

Usage::

    uv run python scripts/run_pathway_residual_model_v2.py
"""

from __future__ import annotations

import argparse
import json
import time
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps

from virtual_cell.analysis import falsification as fal
from virtual_cell.analysis import foundations as fx
from virtual_cell.analysis import loco
from virtual_cell.data import scperteval
from virtual_cell.modelling import pathway_gamma_v2 as v2
from virtual_cell.modelling import pathway_residual as pr

REPO = Path(__file__).resolve().parents[1]
DESIGN = REPO / "data" / "splits" / "four_context_v1"
CANON = REPO / "data" / "processed" / "four_context_v1"
HALVES = REPO / "data" / "processed" / "zero_shot_v1" / "half_deltas.npy"
GMT = REPO / "data" / "raw" / "msigdb"
OUT = REPO / "outputs" / "pathway_gamma_v2"
CL = {c: scperteval.dataset(c).cell_line for c in scperteval.CONTEXTS}
ORDER = ["K562", "RPE1", "HepG2", "Jurkat"]


def gene_level_stats(D, halves, counts, control, perts, genes):
    n_c, n_p = D.shape[0], D.shape[1]
    pair = np.full((n_c, n_c, n_p), np.nan)
    for a in range(n_c):
        for b in range(a + 1, n_c):
            r = fal.rowwise_pearson(D[a], D[b])
            pair[a, b] = pair[b, a] = r
    rel = np.vstack(
        [
            np.nanmean(
                np.vstack(
                    [
                        fal.rowwise_pearson(halves[r, 0, c], halves[r, 1, c])
                        for r in range(halves.shape[0])
                    ]
                ),
                axis=0,
            )
            for c in range(n_c)
        ]
    )
    gi = {g: i for i, g in enumerate(genes)}
    basal = np.vstack(
        [np.array([control[c, gi[p]] if p in gi else np.nan for p in perts]) for c in range(n_c)]
    )
    return pr.GeneLevelStats(
        pair_corr=pair,
        reliability=rel,
        cells=counts.astype(float),
        magnitude=np.linalg.norm(D, axis=2),
        gene_basal=np.nan_to_num(basal),
    )


def representation(M, D, halves, control):
    W = fal.membership_to_weights(M)
    return fal.project(D, W), fal.project(halves, W), fal.project(control, W)


def run_all_folds(
    Y,
    Hp,
    basal_pw,
    in_pw,
    size_pw,
    stats,
    control,
    folds,
    *,
    families=("M0", "M1", "M2"),
    verbose=True,
):
    """Nested-LOCO selection + single outer evaluation for every fold."""

    def weight_fn(target, sources):
        return loco.basal_affine_weights(control, target, list(sources))

    kw = dict(
        stats=stats,
        weight_fn=weight_fn,
        basal_pathway=basal_pw,
        in_pathway=in_pw,
        pathway_size=size_pw,
    )
    rows, details = [], {}
    for fold in folds:
        t, src = fold.target_index, list(fold.source_indices)
        sel = v2.select_by_inner_loco(Y, src, families=families, **kw)
        out = v2.fit_and_predict_outer(Y, t, src, sel, **kw)
        B, R_hat = out["baseline"], out["gamma_hat"]
        obs = Y[t]

        # every family at its own best lambda, for the ladder table
        per_family = {}
        for fam in families:
            sub = {k: v for k, v in sel.inner_scores.items() if k.startswith(fam + "|")}
            best_key = max(sub, key=sub.get)
            a = float(best_key.split("alpha=")[1].split("|")[0])
            lam = float(best_key.split("lam=")[1])
            o = v2.fit_and_predict_outer(Y, t, src, v2.Selection(fam, a, lam), **kw)
            per_family[fam] = {
                "alpha": a,
                "lam": lam,
                **pr.prediction_metrics(obs, o["prediction"]),
                "inner_score": sub[best_key],
            }

        base_m = pr.prediction_metrics(obs, B)
        corr_m = pr.prediction_metrics(obs, out["prediction"])
        r_base = pr.per_perturbation_pearson(obs, B)
        r_corr = pr.per_perturbation_pearson(obs, out["prediction"])
        theory_m = pr.prediction_metrics(obs, out["prediction_theory"])
        boot = v2.bootstrap_delta(obs, B, out["prediction"])
        rows.append(
            {
                "fold": fold.target,
                "cell_line": CL[fold.target],
                "selected_family": sel.family,
                "selected_alpha": sel.alpha,
                "selected_lambda": sel.lam,
                "scale": out["scale"],
                "lambda_theory": out["lambda_theory"],
                **{f"theory_{k}": v for k, v in theory_m.items()},
                "boot_delta": boot["delta_median_pearson"],
                "boot_lo": boot["lo"],
                "boot_hi": boot["hi"],
                **{f"base_{k}": v for k, v in base_m.items()},
                **{f"corr_{k}": v for k, v in corr_m.items()},
                "delta_pearson": corr_m["pearson"] - base_m["pearson"],
                "delta_energy": corr_m["energy_explained"] - base_m["energy_explained"],
                "frac_improved": float(np.nanmean(r_corr > r_base)),
            }
        )
        details[fold.target] = {
            "selection": sel,
            "out": out,
            "per_family": per_family,
            "r_base": r_base,
            "r_corr": r_corr,
            "obs": obs,
            "B": B,
            "R_hat": R_hat,
        }
        if verbose:
            print(
                f"  {CL[fold.target]:8s} sel={sel.family} a={sel.alpha:g} "
                f"lam={sel.lam:g} s={out['scale']:.3f} | "
                f"r {base_m['pearson']:+.4f}->{corr_m['pearson']:+.4f} "
                f"({corr_m['pearson'] - base_m['pearson']:+.4f}) | "
                f"energy {base_m['energy_explained']:+.4f}->"
                f"{corr_m['energy_explained']:+.4f} | "
                f"improved {rows[-1]['frac_improved']:.3f} | "
                f"boot [{boot['lo']:+.4f},{boot['hi']:+.4f}] | "
                f"theory(lam={out['lambda_theory']:.2f}) r={theory_m['pearson']:+.4f}",
                flush=True,
            )
    return pd.DataFrame(rows), details


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-null", type=int, default=20)
    ap.add_argument("--seed", type=int, default=20260920)
    args = ap.parse_args()
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)

    perts = (DESIGN / "shared_perturbations.txt").read_text().split()
    genes = (DESIGN / "shared_genes.txt").read_text().split()
    D = np.load(CANON / "delta_tensor.npy").astype(np.float64)
    control = np.load(CANON / "control_means.npy").astype(np.float64)
    counts = np.load(CANON / "cell_counts.npy")
    halves = np.asarray(np.load(HALVES, mmap_mode="r"), dtype=np.float64)
    folds = loco.make_folds(scperteval.CONTEXTS)
    stats = gene_level_stats(D, halves, counts, control, perts, genes)

    gi = {g: i for i, g in enumerate(genes)}
    results = {}
    for label, fname in (
        ("hallmark", "h.all.v2024.1.Hs.symbols.gmt"),
        ("reactome", "c2.cp.reactome.v2024.1.Hs.symbols.gmt"),
    ):
        sets = fx.read_gmt(GMT / fname)
        names, M = fx.pathway_membership(sets, genes, min_genes=10)
        Y, Hp, basal_pw = representation(M, D, halves, control)
        in_pw = np.vstack([M[:, gi[p]] if p in gi else np.zeros(M.shape[0]) for p in perts])
        size_pw = M.sum(axis=1)
        print("=" * 70)
        print(f"{label.upper()}  ({len(names)} pathways)")
        print("=" * 70)
        tbl, det = run_all_folds(Y, Hp, basal_pw, in_pw, size_pw, stats, control, folds)
        results[label] = {
            "table": tbl,
            "details": det,
            "Y": Y,
            "Hp": Hp,
            "names": names,
            "M": M,
            "in_pw": in_pw,
            "size_pw": size_pw,
            "basal_pw": basal_pw,
        }
        tbl.to_csv(OUT / f"{label}_folds.csv", index=False)

    # ---------------- capacity ladder, Hallmark -------------------------
    print("\n" + "=" * 70)
    print("CAPACITY LADDER (Hallmark, each family at its own inner-best lambda)")
    print("=" * 70)
    ladder = []
    for fold in folds:
        d = results["hallmark"]["details"][fold.target]
        for fam, m in d["per_family"].items():
            ladder.append({"cell_line": CL[fold.target], "family": fam, **m})
    ladder = pd.DataFrame(ladder)
    ladder.to_csv(OUT / "capacity_ladder.csv", index=False)
    print(
        ladder.pivot_table(
            index="cell_line", columns="family", values=["pearson", "energy_explained", "lam"]
        )
        .round(4)
        .to_string()
    )

    # ---------------- diagnostics: what did R_hat learn? ----------------
    print("\n" + "=" * 70)
    print("WHAT DID THE CORRECTION LEARN?  (beta contamination diagnostic)")
    print("=" * 70)
    diag = []
    Yh = results["hallmark"]["Y"]
    for fold in folds:
        t, src = fold.target_index, list(fold.source_indices)
        d = results["hallmark"]["details"][fold.target]
        gamma_true = fal.gamma_of(Yh)[t]
        A = Yh[src].mean(axis=0)
        beta_like = A - A.mean(axis=0, keepdims=True)
        R_true = pr.centre_predictions(d["obs"] - d["B"])
        row = {
            "cell_line": CL[fold.target],
            "r_Rhat_Rtrue": float(np.nanmedian(pr.rowwise_pearson(R_true, d["R_hat"]))),
            "r_Rhat_gamma": float(np.nanmedian(pr.rowwise_pearson(gamma_true, d["R_hat"]))),
            "r_Rhat_betalike": float(np.nanmedian(pr.rowwise_pearson(beta_like, d["R_hat"]))),
            "norm_Rhat": float(np.linalg.norm(d["R_hat"])),
            "norm_Rtrue": float(np.linalg.norm(R_true)),
        }
        diag.append(row)
        print(
            f"  {row['cell_line']:8s} r(Rhat,Rtrue)={row['r_Rhat_Rtrue']:+.4f}  "
            f"r(Rhat,gamma)={row['r_Rhat_gamma']:+.4f}  "
            f"r(Rhat,beta-like)={row['r_Rhat_betalike']:+.4f}  "
            f"||Rhat||/||Rtrue||={row['norm_Rhat'] / row['norm_Rtrue']:.3f}"
        )
    pd.DataFrame(diag).to_csv(OUT / "residual_diagnostics.csv", index=False)

    # ---------------- confidence stratification --------------------------
    print("\n" + "=" * 70)
    print("SOURCE-AGREEMENT CONFIDENCE (raw, inference-available)")
    print("=" * 70)
    conf_rows = []
    for fold in folds:
        t, src = fold.target_index, list(fold.source_indices)
        d = results["hallmark"]["details"][fold.target]
        pairs = [(a, b) for i, a in enumerate(src) for b in src[i + 1 :]]
        agree = np.nanmean(np.vstack([stats.pair_corr[a, b] for a, b in pairs]), axis=0)
        q = pd.qcut(agree, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
        for lab in ["Q1", "Q2", "Q3", "Q4"]:
            m = np.asarray(q == lab)
            if m.sum() == 0:
                continue
            conf_rows.append(
                {
                    "cell_line": CL[fold.target],
                    "quartile": lab,
                    "n": int(m.sum()),
                    "source_agreement": float(np.nanmedian(agree[m])),
                    "base_r": float(np.nanmedian(d["r_base"][m])),
                    "corr_r": float(np.nanmedian(d["r_corr"][m])),
                    "delta_r": float(np.nanmedian(d["r_corr"][m] - d["r_base"][m])),
                }
            )
        sp = sps.spearmanr(agree, d["r_corr"], nan_policy="omit").statistic
        print(f"  {CL[fold.target]:8s} Spearman(source agreement, corrected r) = {sp:+.4f}")
    conf = pd.DataFrame(conf_rows)
    conf.to_csv(OUT / "confidence_strata.csv", index=False)
    print()
    print(
        conf.pivot_table(index="cell_line", columns="quartile", values=["base_r", "corr_r"])
        .round(4)
        .to_string()
    )

    # ---------------- matched-random null --------------------------------
    print("\n" + "=" * 70)
    print(f"MATCHED-RANDOM CONTROL ({args.n_null} replicates, full nested selection each)")
    print("=" * 70)
    M_h = results["hallmark"]["M"]
    null_rows = []
    for i in range(args.n_null):
        rng = np.random.default_rng([args.seed, 4242, i])
        Mn = fal.permuted_membership(M_h, rng)
        Yn, Hn, basal_n = representation(Mn, D, halves, control)
        in_n = np.vstack([Mn[:, gi[p]] if p in gi else np.zeros(Mn.shape[0]) for p in perts])
        tbl, _ = run_all_folds(
            Yn, Hn, basal_n, in_n, Mn.sum(axis=1), stats, control, folds, verbose=False
        )
        tbl["replicate"] = i
        null_rows.append(tbl)
        if (i + 1) % 5 == 0:
            print(f"  {i + 1}/{args.n_null} [{time.time() - t0:.0f}s]", flush=True)
    nulls = pd.concat(null_rows, ignore_index=True)
    nulls.to_csv(OUT / "null_replicates.csv", index=False)

    print("\n  Hallmark vs matched-random null:")
    null_cmp = []
    for metric in ("delta_pearson", "delta_energy", "corr_pearson"):
        for cl in ORDER:
            obs = (
                results["hallmark"]["table"]
                .loc[results["hallmark"]["table"].cell_line == cl, metric]
                .iloc[0]
            )
            st = fal.empirical_percentile(obs, nulls[nulls.cell_line == cl][metric].to_numpy())
            null_cmp.append({"metric": metric, "cell_line": cl, "observed": obs, **st})
            if metric == "delta_pearson":
                print(
                    f"    {cl:8s} delta_r obs={obs:+.4f}  null={st['null_mean']:+.4f}"
                    f"+-{st['null_sd']:.4f}  p={st['p_value']:.3f}  z={st['z']:+.2f}"
                )
    pd.DataFrame(null_cmp).to_csv(OUT / "null_comparison.csv", index=False)

    # ---------------- coefficients ---------------------------------------
    print("\n" + "=" * 70)
    print("M2 STANDARDISED COEFFICIENTS (Hallmark)")
    print("=" * 70)
    coef_rows = []
    for fold in folds:
        d = results["hallmark"]["details"][fold.target]
        model = d["out"]["model"]
        if d["selection"].family != "M2" or not isinstance(model, np.ndarray):
            m2 = d["per_family"].get("M2")
            if m2 is None:
                continue
            o = pr.fit_and_predict_outer(
                results["hallmark"]["Y"],
                fold.target_index,
                list(fold.source_indices),
                pr.Selection("M2", m2["alpha"], m2["lam"]),
                stats=stats,
                weight_fn=lambda t, s: loco.basal_affine_weights(control, t, list(s)),
                basal_pathway=results["hallmark"]["basal_pw"],
                in_pathway=results["hallmark"]["in_pw"],
                pathway_size=results["hallmark"]["size_pw"],
            )
            model = o["model"]
        for name, c in zip(pr.FEATURE_NAMES, model[1:], strict=True):
            coef_rows.append(
                {"cell_line": CL[fold.target], "feature": name, "coefficient": float(c)}
            )
    coefs = pd.DataFrame(coef_rows)
    coefs.to_csv(OUT / "coefficients.csv", index=False)
    piv = coefs.pivot_table(index="feature", columns="cell_line", values="coefficient")
    piv["mean_abs"] = piv.abs().mean(axis=1)
    piv["sign_consistent"] = np.sign(piv[ORDER]).abs().sum(axis=1) == np.abs(
        np.sign(piv[ORDER]).sum(axis=1)
    )
    print(piv.sort_values("mean_abs", ascending=False).round(4).to_string())

    (OUT / "summary.json").write_text(
        json.dumps(
            {
                "generated": date.today().isoformat(),
                "runtime_seconds": time.time() - t0,
                "n_null": args.n_null,
                "lambda_grid": list(pr.LAMBDA_GRID),
                "alpha_grid": list(pr.ALPHA_GRID),
                "features": list(pr.FEATURE_NAMES),
                "hallmark_folds": results["hallmark"]["table"].to_dict("records"),
                "reactome_folds": results["reactome"]["table"].to_dict("records"),
                "capacity_ladder": ladder.to_dict("records"),
                "residual_diagnostics": diag,
                "null_comparison": null_cmp,
                "confidence": conf.to_dict("records"),
                "coefficients": coefs.to_dict("records"),
            },
            indent=2,
            default=str,
        )
    )
    print(f"\nTotal runtime: {(time.time() - t0) / 60:.1f} min")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
