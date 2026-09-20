"""Reliability-aware zero-shot transferability / confidence model v1.

The point predictor is **frozen**: scale-calibrated source-only conserved
transfer. This study only predicts how trustworthy that prediction will be,
before the target perturbation is measured.

Nested LOCO throughout: every standardisation, model choice, penalty and
calibration is fitted on the three source contexts only. The outer target's
responses are used once, for evaluation.

Usage::

    uv run python scripts/run_transferability_confidence.py
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
from virtual_cell.analysis import loco
from virtual_cell.data import scperteval
from virtual_cell.modelling import pathway_residual as pr
from virtual_cell.modelling import transferability as tf

REPO = Path(__file__).resolve().parents[1]
DESIGN = REPO / "data" / "splits" / "four_context_v1"
CANON = REPO / "data" / "processed" / "four_context_v1"
HALVES = REPO / "data" / "processed" / "zero_shot_v1" / "half_deltas.npy"
OUT = REPO / "outputs" / "transferability_v1"
CL = {c: scperteval.dataset(c).cell_line for c in scperteval.CONTEXTS}
ORDER = ["K562", "RPE1", "HepG2", "Jurkat"]


def build_stats(D, halves, counts, control, perts, genes):
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
    return tf.SourceStats(
        pair_corr=pair,
        reliability=rel,
        cells=counts.astype(float),
        magnitude=np.linalg.norm(D, axis=2),
        gene_basal=basal,
        basal_sim=loco.basal_similarity(control),
    )


def fold_targets(D, halves, target, sources, min_signal):
    """Frozen predictor B plus every reliability-aware evaluation target."""
    scale = pr.fit_scale(D, sources)
    B = pr.baseline(D, sources, scale)
    h1 = np.asarray(halves[:, 0, target], dtype=np.float64).mean(axis=0)
    h2 = np.asarray(halves[:, 1, target], dtype=np.float64).mean(axis=0)
    sig = tf.signal_energy(h1, h2)
    res = tf.residual_energy(h1, h2, B)
    d, stable = tf.normalised_d(h1, h2, B, min_signal=min_signal)
    rel_tab = (
        sum(
            loco.reliability_normalised_evaluation(
                B,
                np.asarray(halves[r, 0, target], dtype=np.float64),
                np.asarray(halves[r, 1, target], dtype=np.float64),
            )
            for r in range(halves.shape[0])
        )
        / halves.shape[0]
    )
    return {
        "scale": scale,
        "B": B,
        "signal": sig,
        "residual": res,
        "D": d,
        "stable": stable,
        "r_disatt": rel_tab["r_disattenuated"].to_numpy(),
        "pearson": fal.rowwise_pearson(D[target], B),
        "cosine": np.einsum("pg,pg->p", D[target], B)
        / (np.linalg.norm(D[target], axis=1) * np.linalg.norm(B, axis=1) + 1e-12),
    }


def coverage_report(confidence, t, coverages=tf.COVERAGE_LEVELS):
    c = np.asarray(confidence, dtype=np.float64)
    ok = np.isfinite(c)
    order = np.argsort(-c[ok], kind="stable")
    idx = np.flatnonzero(ok)[order]
    rows = []
    for cov in coverages:
        k = max(1, int(round(cov * idx.size)))
        s = idx[:k]
        rows.append(
            {
                "coverage": cov,
                "n": k,
                "mean_residual_energy": float(np.mean(t["residual"][s])),
                "median_residual_energy": float(np.median(t["residual"][s])),
                "median_D": float(np.nanmedian(t["D"][s])),
                "median_r_disatt": float(np.nanmedian(t["r_disatt"][s])),
                "median_pearson": float(np.nanmedian(t["pearson"][s])),
                "median_cosine": float(np.nanmedian(t["cosine"][s])),
            }
        )
    return rows


def boot_ci(values, n_boot=2000, seed=0):
    v = np.asarray(values, dtype=np.float64)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    b = np.mean(rng.choice(v, size=(n_boot, v.size), replace=True), axis=1)
    return float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seed", type=int, default=20260920)
    args = ap.parse_args()
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)

    rule = json.loads((OUT / "stability_rule.json").read_text())
    min_signal = float(rule["min_signal_energy"])
    print(f"PREDECLARED stability rule: min signal energy = {min_signal:.4f}")

    perts = (DESIGN / "shared_perturbations.txt").read_text().split()
    genes = (DESIGN / "shared_genes.txt").read_text().split()
    D = np.load(CANON / "delta_tensor.npy").astype(np.float64)
    control = np.load(CANON / "control_means.npy").astype(np.float64)
    counts = np.load(CANON / "cell_counts.npy")
    halves = np.asarray(np.load(HALVES, mmap_mode="r"), dtype=np.float64)
    folds = loco.make_folds(scperteval.CONTEXTS)
    stats = build_stats(D, halves, counts, control, perts, genes)

    # ---- per-fold targets and features ------------------------------------
    T, F = {}, {}
    for fold in folds:
        t = fold.target_index
        src = list(fold.source_indices)
        T[t] = fold_targets(D, halves, t, src, min_signal)
        F[t] = tf.build_confidence_features(D, halves, t, src, stats=stats)
    print(
        "\nstable fraction by context: "
        + ", ".join(f"{CL[scperteval.CONTEXTS[t]]}={T[t]['stable'].mean():.3f}" for t in T)
    )

    # ---- M0 / M1 / M2 with nested inner LOCO ------------------------------
    print("\n" + "=" * 70)
    print("MODEL LADDER (inner pseudo-LOCO selection, outer evaluated once)")
    print("=" * 70)
    ag_idx = tf.FEATURE_NAMES.index("source_agreement")
    scores, chosen, coefs = {}, {}, []
    for fold in folds:
        t, src = fold.target_index, list(fold.source_indices)
        # inner rows: each source context as a pseudo-target, its own sources
        inner = {}
        for c in src:
            isrc = [s for s in src if s != c]
            inner[c] = {
                "X": tf.build_confidence_features(D, halves, c, isrc, stats=stats),
                **fold_targets(D, halves, c, isrc, min_signal),
            }

        def inner_score(make_conf, src=src, inner=inner):
            """Mean negative AURC over the inner pseudo-LOCO folds.

            ``src`` and ``inner`` are bound as defaults so the closure cannot
            pick up a later loop iteration's values.
            """
            vals = []
            for held in src:
                train = [c for c in src if c != held]
                Xtr = np.vstack([inner[c]["X"] for c in train])
                ytr = np.concatenate([inner[c]["D"] for c in train])
                mtr = np.concatenate([inner[c]["stable"] for c in train])
                conf = make_conf(Xtr[mtr], ytr[mtr], inner[held]["X"])
                vals.append(-tf.area_under_risk_coverage(conf, inner[held]["residual"]))
            return float(np.mean(vals))

        def m0(_Xtr, _ytr, Xte):
            return Xte[:, ag_idx]

        def m1(Xtr, ytr, Xte):
            iso = tf.fit_monotone(Xtr[:, ag_idx], ytr)
            return -iso.predict(Xte[:, ag_idx])

        cands = {"M0": (m0, None), "M1": (m1, None)}
        for name, feats in tf.ABLATIONS.items():
            cols = [tf.FEATURE_NAMES.index(f) for f in feats]
            for a in tf.CONFIDENCE_ALPHAS:

                def m2(Xtr, ytr, Xte, cols=cols, a=a):
                    sc = tf.Standardiser().fit(Xtr[:, cols])
                    coef = tf.fit_ridge(sc.transform(Xtr[:, cols]), ytr, a)
                    return -tf.predict_ridge(sc.transform(Xte[:, cols]), coef)

                cands[f"M2:{name}:a={a:g}"] = (m2, (cols, a))

        results = {k: inner_score(fn) for k, (fn, _) in cands.items()}
        best = max(results, key=results.get)
        chosen[t] = {"best": best, "inner_scores": results}

        # refit on all three source pseudo-targets, apply to the outer target
        Xtr = np.vstack([inner[c]["X"] for c in src])
        ytr = np.concatenate([inner[c]["D"] for c in src])
        mtr = np.concatenate([inner[c]["stable"] for c in src])
        fn = cands[best][0]
        learned = fn(Xtr[mtr], ytr[mtr], F[t])
        scores[t] = {
            "M0": F[t][:, ag_idx],
            "M1": cands["M1"][0](Xtr[mtr], ytr[mtr], F[t]),
            "learned": learned,
            "random": np.random.default_rng([args.seed, t]).normal(size=len(perts)),
            "source_reliability": F[t][:, tf.FEATURE_NAMES.index("source_reliability")],
            "source_magnitude": -F[t][:, tf.FEATURE_NAMES.index("source_magnitude")],
        }
        # coefficients of the full ridge, for interpretation
        cols = list(range(len(tf.FEATURE_NAMES)))
        sc = tf.Standardiser().fit(Xtr[mtr][:, cols])
        cf = tf.fit_ridge(sc.transform(Xtr[mtr][:, cols]), ytr[mtr], 10.0)
        for nm, v in zip(tf.FEATURE_NAMES, cf[1:], strict=True):
            coefs.append({"cell_line": CL[fold.target], "feature": nm, "coefficient": float(v)})
        print(f"  {CL[fold.target]:8s} selected: {best}")

    pd.DataFrame(coefs).to_csv(OUT / "coefficients.csv", index=False)
    (OUT / "selection.json").write_text(
        json.dumps({CL[scperteval.CONTEXTS[t]]: chosen[t] for t in chosen}, indent=2, default=str)
    )

    # ---- E/H: how well does each score rank? ------------------------------
    print("\n" + "=" * 70)
    print("SELECTIVE PREDICTION — Spearman(score, -D) and AURC on residual energy")
    print("=" * 70)
    rank_rows, rc_rows, prio_rows, calib_rows = [], [], [], []
    for fold in folds:
        t = fold.target_index
        cl = CL[fold.target]
        tt = T[t]
        for name, sc in scores[t].items():
            m = np.isfinite(sc) & tt["stable"]
            rho_d = sps.spearmanr(sc[m], -tt["D"][m]).statistic if m.sum() > 10 else np.nan
            rho_r = sps.spearmanr(sc[np.isfinite(sc)], -tt["residual"][np.isfinite(sc)]).statistic
            rho_sim = sps.spearmanr(sc, tt["r_disatt"], nan_policy="omit").statistic
            aurc = tf.area_under_risk_coverage(sc, tt["residual"])
            rank_rows.append(
                {
                    "cell_line": cl,
                    "score": name,
                    "spearman_vs_negD": rho_d,
                    "spearman_vs_neg_residual": rho_r,
                    "spearman_vs_similarity": rho_sim,
                    "aurc": aurc,
                }
            )
            for r in coverage_report(sc, tt):
                rc_rows.append({"cell_line": cl, "score": name, **r})
            for r in tf.prioritisation_capture(sc, tt["residual"]):
                prio_rows.append({"cell_line": cl, "score": name, **r})
        for b in tf.calibration_curve(scores[t]["M0"], tt["r_disatt"], n_bins=5):
            calib_rows.append({"cell_line": cl, "score": "M0", **b})
        for b in tf.calibration_curve(scores[t]["learned"], tt["r_disatt"], n_bins=5):
            calib_rows.append({"cell_line": cl, "score": "learned", **b})

    ranks = pd.DataFrame(rank_rows)
    rc = pd.DataFrame(rc_rows)
    prio = pd.DataFrame(prio_rows)
    calib = pd.DataFrame(calib_rows)
    for df, nm in (
        (ranks, "ranking"),
        (rc, "risk_coverage"),
        (prio, "prioritisation"),
        (calib, "calibration"),
    ):
        df.to_csv(OUT / f"{nm}.csv", index=False)

    print(
        ranks.pivot_table(index="cell_line", columns="score", values="spearman_vs_negD")
        .round(4)
        .to_string()
    )
    print("\nAURC on reproducible residual energy (lower is better):")
    print(ranks.pivot_table(index="cell_line", columns="score", values="aurc").round(4).to_string())

    print("\n" + "=" * 70)
    print("RISK-COVERAGE (median reliability-normalised similarity)")
    print("=" * 70)
    for name in ("M0", "learned", "random"):
        print(f"\n  {name}:")
        print(
            rc[rc.score == name]
            .pivot_table(index="cell_line", columns="coverage", values="median_r_disatt")
            .round(4)
            .to_string()
        )

    print("\n" + "=" * 70)
    print("EXPERIMENT PRIORITISATION — fraction of reproducible error captured")
    print("=" * 70)
    for name in ("M0", "learned", "random"):
        print(f"\n  {name}:")
        print(
            prio[prio.score == name]
            .pivot_table(index="cell_line", columns="budget", values="captured_fraction")
            .round(4)
            .to_string()
        )

    summary = {
        "generated": date.today().isoformat(),
        "runtime_seconds": time.time() - t0,
        "stability_rule": {k: v for k, v in rule.items() if k != "curve"},
        "stable_fraction": {CL[scperteval.CONTEXTS[t]]: float(T[t]["stable"].mean()) for t in T},
        "scale": {CL[scperteval.CONTEXTS[t]]: float(T[t]["scale"]) for t in T},
        "selected": {CL[scperteval.CONTEXTS[t]]: chosen[t]["best"] for t in chosen},
        "ranking": ranks.to_dict("records"),
        "risk_coverage": rc.to_dict("records"),
        "prioritisation": prio.to_dict("records"),
        "calibration": calib.to_dict("records"),
        "coefficients": coefs,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"\nTotal runtime: {(time.time() - t0) / 60:.1f} min")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
