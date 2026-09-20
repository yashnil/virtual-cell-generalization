"""Transferability foundations (v1): template recovery, source agreement, pathway gamma.

Three questions, all leave-one-context-out on the frozen design:

  A  can basal control expression recover the context template alpha?
  B  does source agreement predict transferability within each fold, and does it
     survive inference-available confounds?
  C  does gamma become more recoverable at pathway resolution?
  D  candidate D targets are computed and characterised, but NOT selected.

No model is trained beyond a single global scalar and a rank-<=2 ridge map.

Usage::

    uv run python scripts/run_transferability_foundations.py
"""

from __future__ import annotations

import argparse
import json
import time
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from virtual_cell.analysis import foundations as fx
from virtual_cell.analysis import loco
from virtual_cell.data import scperteval
from virtual_cell.decomposition import anova

REPO_ROOT = Path(__file__).resolve().parents[1]
DESIGN_DIR = REPO_ROOT / "data" / "splits" / "four_context_v1"
CANONICAL_DIR = REPO_ROOT / "data" / "processed" / "four_context_v1"
HALVES = REPO_ROOT / "data" / "processed" / "zero_shot_v1" / "half_deltas.npy"
GMT_DIR = REPO_ROOT / "data" / "raw" / "msigdb"
OUT = REPO_ROOT / "outputs" / "transferability_foundations_v1"

CL = {c: scperteval.dataset(c).cell_line for c in scperteval.CONTEXTS}
MIN_PATHWAY_GENES = 10


def med(x) -> float:
    return float(np.nanmedian(np.asarray(x, dtype=np.float64)))


# --------------------------------------------------------------------------


def part_a(D, control, folds, halves) -> tuple[pd.DataFrame, pd.DataFrame]:
    print("=" * 70)
    print("A. BASAL TEMPLATE / SCALE RECOVERABILITY")
    print("=" * 70)

    diag, perf = [], []
    print("\n--- descriptive: basal deviation vs true alpha (evaluation only) ---")
    for fold in folds:
        c = fold.target_index
        alpha_true = fx.oracle_alpha_evaluation_only(D, c)
        dev_all = fx.basal_deviation(control, c, range(len(scperteval.CONTEXTS)))
        r = float(np.corrcoef(alpha_true, dev_all)[0, 1])
        cos = float(alpha_true @ dev_all / (np.linalg.norm(alpha_true) * np.linalg.norm(dev_all)))
        r_abs = float(stats.spearmanr(np.abs(alpha_true), np.abs(dev_all)).statistic)
        # subspace alignment: how much of alpha lies in the span of source basal deviations
        devs_S = np.vstack(
            [fx.basal_deviation(control, s, fold.source_indices) for s in fold.source_indices]
        )
        q, _ = np.linalg.qr(devs_S.T)
        proj = q @ (q.T @ alpha_true)
        frac_in_span = float((proj @ proj) / (alpha_true @ alpha_true))
        alphas_S = fx.source_alphas(D, fold.source_indices)
        k = fx.fit_global_scale(alphas_S, devs_S)
        diag.append(
            {
                "fold": fold.target,
                "cell_line": CL[fold.target],
                "r_alpha_basaldev": r,
                "cosine": cos,
                "spearman_abs": r_abs,
                "alpha_frac_in_source_basal_span": frac_in_span,
                "global_scalar_k": k,
                "alpha_norm": float(np.linalg.norm(alpha_true)),
                "basal_dev_norm": float(np.linalg.norm(dev_all)),
            }
        )
        print(
            f"  {CL[fold.target]:8s} r={r:+.4f} cos={cos:+.4f} "
            f"spearman(|a|,|dev|)={r_abs:+.4f} frac_in_span={frac_in_span:.4f} k={k:+.4f}"
        )

    ks = [d["global_scalar_k"] for d in diag]
    print(
        f"\n  global scalar k across folds: {[round(x, 4) for x in ks]}  "
        f"(sign consistent: {all(np.sign(x) == np.sign(ks[0]) for x in ks)})"
    )

    print("\n--- LOCO performance with template correction ---")
    for fold in folds:
        c = fold.target_index
        Y = D[c]
        A = loco.source_mean(D, fold.source_indices)
        ests = fx.template_estimators(D, control, c, fold.source_indices)
        ests.append(
            fx.TemplateEstimate("oracle_alpha_EVAL_ONLY", fx.oracle_alpha_evaluation_only(D, c), {})
        )
        shrink = fx.fit_response_shrinkage(D, fold.source_indices)
        extra = [
            ("scale_only", fx.apply_scale_correction(A, shrink)),
            (
                "scale_plus_oracle_alpha_EVAL_ONLY",
                fx.apply_scale_correction(A, shrink, fx.oracle_alpha_evaluation_only(D, c)),
            ),
        ]
        print(f"\n  held out {CL[fold.target]}:  (fitted shrinkage s={shrink:.4f})")
        named = [(e.name, fx.apply_template_correction(A, e.alpha_hat)) for e in ests] + extra
        for est_name, P in named:
            m = loco.per_perturbation_metrics(Y, P, with_spearman=False)
            pooled = fx.pooled_unexplained_fraction(
                np.asarray(halves[:, 0, c], dtype=np.float64).mean(axis=0),
                np.asarray(halves[:, 1, c], dtype=np.float64).mean(axis=0),
                P,
            )
            perf.append(
                {
                    "fold": fold.target,
                    "cell_line": CL[fold.target],
                    "estimator": est_name,
                    "shrinkage": shrink,
                    "pearson": med(m["pearson"]),
                    "cosine": med(m["cosine"]),
                    "mse": med(m["mse"]),
                    "energy_explained": med(m["energy_explained"]),
                    "pooled_unexplained_reliable": pooled,
                }
            )
            print(
                f"    {est_name:34s} r={med(m['pearson']):+.4f} cos={med(m['cosine']):+.4f} "
                f"mse={med(m['mse']):.4f} energy={med(m['energy_explained']):+.4f} "
                f"unexpl_reliable={pooled:+.4f}"
            )
    return pd.DataFrame(diag), pd.DataFrame(perf)


def part_b(D, control, folds, halves, perts, genes) -> tuple[pd.DataFrame, pd.DataFrame]:
    print("\n" + "=" * 70)
    print("B. SOURCE AGREEMENT VALIDATION")
    print("=" * 70)
    counts = np.load(CANONICAL_DIR / "cell_counts.npy")
    gene_index = {g: i for i, g in enumerate(genes)}
    rows, summary = [], []

    for fold in folds:
        c = fold.target_index
        S = list(fold.source_indices)
        Y = D[c]
        A = loco.source_mean(D, S)
        m = loco.per_perturbation_metrics(Y, A, with_spearman=False)

        # target reliability and normalised success (evaluation-only reliability)
        rel = (
            sum(
                loco.reliability_normalised_evaluation(
                    A,
                    np.asarray(halves[r, 0, c], dtype=np.float64),
                    np.asarray(halves[r, 1, c], dtype=np.float64),
                )
                for r in range(halves.shape[0])
            )
            / halves.shape[0]
        )

        # ---- inference-available features (source contexts + target basal only)
        pair_r, src_mag, src_rel, min_pair = [], [], [], []
        for pi in range(len(perts)):
            rs = [
                float(np.corrcoef(D[a, pi], D[b, pi])[0, 1])
                for ii, a in enumerate(S)
                for b in S[ii + 1 :]
            ]
            pair_r.append(float(np.nanmean(rs)))
            min_pair.append(float(np.nanmin(rs)))
            src_mag.append(float(np.mean([np.linalg.norm(D[s, pi]) for s in S])))
        # source-side reliability of the source mean, from source cells only
        src_rel = np.zeros(len(perts))
        for r in range(halves.shape[0]):
            a = np.asarray(halves[r, 0][S], dtype=np.float64).mean(axis=0)
            b = np.asarray(halves[r, 1][S], dtype=np.float64).mean(axis=0)
            src_rel += loco.split_half_reliability(a, b) / halves.shape[0]
        src_cells = counts[S].mean(axis=0)
        tgt_gene_basal = np.array(
            [control[c, gene_index[p]] if p in gene_index else np.nan for p in perts]
        )

        frame = pd.DataFrame(
            {
                "fold": fold.target,
                "cell_line": CL[fold.target],
                "perturbation": perts,
                "success_raw": m["pearson"].to_numpy(),
                "success_norm": rel["r_disattenuated"].to_numpy(),
                "target_rho_half": rel["rho_half"].to_numpy(),
                "source_agreement": pair_r,
                "source_min_pair_agreement": min_pair,
                "source_magnitude": src_mag,
                "source_reliability": src_rel,
                "source_cells": src_cells,
                "target_gene_basal": tgt_gene_basal,
            }
        )
        rows.append(frame)

        def boot(x, y, method):
            ok = np.isfinite(x) & np.isfinite(y)
            xv, yv = np.asarray(x)[ok], np.asarray(y)[ok]
            f = (
                (lambda u, v: float(stats.spearmanr(u, v).statistic))
                if method == "spearman"
                else (lambda u, v: float(np.corrcoef(u, v)[0, 1]))
            )
            point = f(xv, yv)
            rng = np.random.default_rng(0)
            bs = [f(xv[i], yv[i]) for i in (rng.integers(0, len(xv), size=(400, len(xv))))]
            return (
                point,
                float(np.percentile(bs, 2.5)),
                float(np.percentile(bs, 97.5)),
                int(ok.sum()),
            )

        print(f"\n  held out {CL[fold.target]}:")
        for target_col in ("success_raw", "success_norm"):
            sp, lo, hi, n = boot(frame["source_agreement"], frame[target_col], "spearman")
            pe, plo, phi, _ = boot(frame["source_agreement"], frame[target_col], "pearson")
            print(
                f"    vs {target_col:12s} Spearman={sp:+.4f} [{lo:+.4f},{hi:+.4f}]  "
                f"Pearson={pe:+.4f} [{plo:+.4f},{phi:+.4f}]  n={n}"
            )
            summary.append(
                {
                    "fold": fold.target,
                    "cell_line": CL[fold.target],
                    "target": target_col,
                    "spearman": sp,
                    "sp_lo": lo,
                    "sp_hi": hi,
                    "pearson": pe,
                    "pe_lo": plo,
                    "pe_hi": phi,
                    "n": n,
                }
            )

        # NOTE: source_min_pair_agreement is a *component* of source_agreement
        # (the min of the same pairwise correlations whose mean is the feature).
        # Adjusting for it is over-adjustment, so it is reported individually but
        # excluded from the joint confound set.
        confounds = [
            "source_magnitude",
            "source_reliability",
            "source_cells",
            "target_gene_basal",
        ]
        for extra_conf in ["source_min_pair_agreement"]:
            pc = fx.partial_spearman(
                frame["source_agreement"],
                frame["success_norm"],
                frame[[extra_conf]].to_numpy(),
            )
            print(f"      [over-adjustment, not a confound] {extra_conf} = {pc:+.4f}")
        for conf in confounds:
            pc = fx.partial_spearman(
                frame["source_agreement"], frame["success_norm"], frame[[conf]].to_numpy()
            )
            print(f"      partial (control {conf:26s}) = {pc:+.4f}")
        allc = fx.partial_spearman(
            frame["source_agreement"], frame["success_norm"], frame[confounds].to_numpy()
        )
        print(f"      partial (control ALL confounds)          = {allc:+.4f}")
        summary.append(
            {
                "fold": fold.target,
                "cell_line": CL[fold.target],
                "target": "partial_all_confounds",
                "spearman": allc,
            }
        )
    return pd.concat(rows, ignore_index=True), pd.DataFrame(summary)


def part_c(D, folds, halves, genes) -> pd.DataFrame:
    print("\n" + "=" * 70)
    print("C. PATHWAY-LEVEL GAMMA RECOVERABILITY")
    print("=" * 70)
    out = []
    for label, fname in (
        ("hallmark", "h.all.v2024.1.Hs.symbols.gmt"),
        ("reactome", "c2.cp.reactome.v2024.1.Hs.symbols.gmt"),
    ):
        sets = fx.read_gmt(GMT_DIR / fname)
        names, M = fx.pathway_membership(sets, genes, min_genes=MIN_PATHWAY_GENES)
        print(
            f"\n  {label}: {len(sets)} sets released, {len(names)} with "
            f">={MIN_PATHWAY_GENES} genes in the frozen space"
        )
        Dp = fx.pathway_scores(D, M)
        dec_p = anova.decompose(Dp)
        Hp = np.stack(
            [fx.pathway_scores(np.asarray(halves[:, h], dtype=np.float64), M) for h in (0, 1)],
            axis=1,
        )
        for fold in folds:
            c = fold.target_index
            Ap = loco.source_mean(Dp, fold.source_indices)
            Apc = Ap - Ap.mean(axis=0, keepdims=True)
            gp_true = dec_p.gamma[c]
            # reliability of pathway gamma: split every context, per repeat
            rho = np.nanmean(
                np.vstack(
                    [
                        loco.split_half_reliability(
                            anova.decompose(Hp[r, 0]).gamma[c], anova.decompose(Hp[r, 1]).gamma[c]
                        )
                        for r in range(Hp.shape[0])
                    ]
                ),
                axis=0,
            )
            rho_full = np.asarray(loco.spearman_brown(rho))
            control_means = np.load(CANONICAL_DIR / "control_means.npy").astype(np.float64)
            nearest, _ = loco.nearest_context(control_means, c, fold.source_indices)
            w_affine = loco.basal_affine_weights(control_means, c, fold.source_indices)
            baselines = (
                ("nearest_basal", loco.single_source(Dp, nearest)),
                ("basal_affine", loco.weighted_source(Dp, fold.source_indices, w_affine)),
            )
            for bname, P in baselines:
                Pc = P - P.mean(axis=0, keepdims=True)
                g_hat = Pc - Apc
                rg = np.array(
                    [loco._safe_pearson(gp_true[i], g_hat[i]) for i in range(gp_true.shape[0])]
                )
                out.append(
                    {
                        "collection": label,
                        "n_pathways": len(names),
                        "fold": fold.target,
                        "cell_line": CL[fold.target],
                        "baseline": bname,
                        "r_gamma": med(rg),
                        "gamma_rho_half": med(rho),
                        "gamma_rho_full": med(rho_full),
                        "gamma_ceiling_full": float(np.sqrt(max(med(rho_full), 0.0))),
                        "r_gamma_norm": med(rg) / max(np.sqrt(max(med(rho_full), 1e-9)), 1e-9),
                    }
                )
                print(
                    f"    {CL[fold.target]:8s} {bname:14s} r_gamma={med(rg):+.4f} "
                    f"rho_full={med(rho_full):.4f} ceiling={np.sqrt(max(med(rho_full), 0)):.4f} "
                    f"normalised={out[-1]['r_gamma_norm']:+.4f}"
                )
    return pd.DataFrame(out)


def part_d(D, folds, halves) -> pd.DataFrame:
    print("\n" + "=" * 70)
    print("D. CANDIDATE D TARGETS (characterised, not selected)")
    print("=" * 70)
    rows = []
    for fold in folds:
        c = fold.target_index
        A = loco.source_mean(D, fold.source_indices)
        h1 = np.asarray(halves[:, 0, c], dtype=np.float64).mean(axis=0)
        h2 = np.asarray(halves[:, 1, c], dtype=np.float64).mean(axis=0)
        tab = fx.candidate_targets(h1, h2, A)
        tab["fold"] = fold.target
        tab["cell_line"] = CL[fold.target]
        rows.append(tab)
        print(
            f"  {CL[fold.target]:8s} pooled_unexplained="
            f"{fx.pooled_unexplained_fraction(h1, h2, A):+.4f}  "
            f"per-pair median={med(tab['D_unexplained_fraction']):+.4f}  "
            f"usable={float(tab['usable'].mean()):.3f}  "
            f"frac>1 (worse than zero)={float((tab['D_unexplained_fraction'] > 1).mean()):.3f}"
        )
    return pd.concat(rows, ignore_index=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.parse_args()
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    perts = (DESIGN_DIR / "shared_perturbations.txt").read_text().split()
    genes = (DESIGN_DIR / "shared_genes.txt").read_text().split()
    D = np.load(CANONICAL_DIR / "delta_tensor.npy").astype(np.float64)
    control = np.load(CANONICAL_DIR / "control_means.npy").astype(np.float64)
    halves = np.load(HALVES, mmap_mode="r")
    folds = loco.make_folds(scperteval.CONTEXTS)

    diag, perf = part_a(D, control, folds, halves)
    diag.to_csv(OUT / "template_diagnostics.csv", index=False)
    perf.to_csv(OUT / "template_performance.csv", index=False)

    feats, bsum = part_b(D, control, folds, halves, perts, genes)
    feats.to_csv(OUT / "source_agreement_features.csv", index=False)
    bsum.to_csv(OUT / "source_agreement_summary.csv", index=False)

    pathc = part_c(D, folds, halves, genes)
    pathc.to_csv(OUT / "pathway_gamma.csv", index=False)

    dcand = part_d(D, folds, halves)
    dcand.to_csv(OUT / "candidate_d_targets.csv", index=False)

    (OUT / "summary.json").write_text(
        json.dumps(
            {
                "generated": date.today().isoformat(),
                "runtime_seconds": time.time() - t0,
                "template_diagnostics": diag.to_dict("records"),
                "template_performance": perf.to_dict("records"),
                "source_agreement": bsum.to_dict("records"),
                "pathway_gamma": pathc.to_dict("records"),
            },
            indent=2,
            default=str,
        )
    )
    print(f"\nTotal runtime: {(time.time() - t0) / 60:.1f} min")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
