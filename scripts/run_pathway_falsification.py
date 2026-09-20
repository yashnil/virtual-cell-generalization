"""Representation falsification battery for pathway-level gamma recoverability.

Asks whether the pathway-level gamma gain is *biology* or merely aggregation:

  B  Hallmark vs size/geometry-matched random gene sets
  C  Hallmark vs dense random projections of equal dimensionality
  D  an independent ontology (Reactome) under the same rule
  E  K562 without Jurkat, and Jurkat without K562

All representations are linear maps applied identically; only the map differs.
Null constructions are fixed in advance and are never adjusted after a result.

No model is trained.

Usage::

    uv run python scripts/run_pathway_falsification.py
"""

from __future__ import annotations

import argparse
import json
import time
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

from virtual_cell.analysis import falsification as fal
from virtual_cell.analysis import foundations as fx
from virtual_cell.analysis import loco
from virtual_cell.data import scperteval

REPO_ROOT = Path(__file__).resolve().parents[1]
DESIGN_DIR = REPO_ROOT / "data" / "splits" / "four_context_v1"
CANONICAL = REPO_ROOT / "data" / "processed" / "four_context_v1"
HALVES = REPO_ROOT / "data" / "processed" / "zero_shot_v1" / "half_deltas.npy"
GMT_DIR = REPO_ROOT / "data" / "raw" / "msigdb"
OUT = REPO_ROOT / "outputs" / "pathway_falsification_v1"

CL = {c: scperteval.dataset(c).cell_line for c in scperteval.CONTEXTS}
ORDER = ["K562", "RPE1", "HepG2", "Jurkat"]
MIN_PATHWAY_GENES = 10
BASELINE = "basal_affine"  # predeclared: the better of the two in foundations


def evaluate(W, D, halves, control, folds, *, with_reliability=True, source_override=None):
    """Recoverability for one representation across all folds."""
    Dp = fal.project(D, W)
    Hp = fal.project(halves, W) if with_reliability else None
    rows = []
    for fold in folds:
        c = fold.target_index
        sources = (
            list(source_override[fold.target]) if source_override else list(fold.source_indices)
        )
        w = loco.basal_affine_weights(control, c, sources)
        r_gamma = fal.gamma_recoverability(Dp, c, sources, weights=w)
        row = {
            "fold": fold.target,
            "cell_line": CL[fold.target],
            "n_sources": len(sources),
            "sources": "+".join(CL[scperteval.CONTEXTS[s]] for s in sources),
            "r_gamma": float(np.nanmedian(r_gamma)),
        }
        if with_reliability:
            rho = fal.gamma_reliability(Hp, c)
            rho_full = fal.spearman_brown(rho)
            ceil = float(np.sqrt(max(float(np.nanmedian(rho_full)), 0.0)))
            row.update(
                {
                    "rho_half": float(np.nanmedian(rho)),
                    "rho_full": float(np.nanmedian(rho_full)),
                    "ceiling": ceil,
                    "r_gamma_normalised": row["r_gamma"] / ceil if ceil > 1e-9 else np.nan,
                }
            )
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-null", type=int, default=100)
    ap.add_argument("--seed", type=int, default=20260920)
    args = ap.parse_args()

    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    genes = (DESIGN_DIR / "shared_genes.txt").read_text().split()
    D = np.load(CANONICAL / "delta_tensor.npy").astype(np.float64)
    control = np.load(CANONICAL / "control_means.npy").astype(np.float64)
    halves = np.asarray(np.load(HALVES, mmap_mode="r"), dtype=np.float64)
    folds = loco.make_folds(scperteval.CONTEXTS)
    n_genes = len(genes)

    print("=" * 70)
    print("OBSERVED BIOLOGICAL REPRESENTATIONS")
    print("=" * 70)
    observed, memberships = {}, {}
    for label, fname in (
        ("hallmark", "h.all.v2024.1.Hs.symbols.gmt"),
        ("reactome", "c2.cp.reactome.v2024.1.Hs.symbols.gmt"),
    ):
        sets = fx.read_gmt(GMT_DIR / fname)
        names, M = fx.pathway_membership(sets, genes, min_genes=MIN_PATHWAY_GENES)
        memberships[label] = M
        res = evaluate(fal.membership_to_weights(M), D, halves, control, folds)
        res["representation"] = label
        observed[label] = res
        print(f"\n  {label}: {len(names)} sets (of {len(sets)} released)")
        print(
            res[["cell_line", "r_gamma", "rho_full", "ceiling", "r_gamma_normalised"]]
            .round(4)
            .to_string(index=False)
        )

    gene_level = evaluate(np.eye(n_genes), D, halves, control, folds)
    gene_level["representation"] = "genes"
    print("\n  gene level (identity map):")
    print(
        gene_level[["cell_line", "r_gamma", "rho_full", "ceiling", "r_gamma_normalised"]]
        .round(4)
        .to_string(index=False)
    )

    # ---------------------------------------------------------------- nulls
    print("\n" + "=" * 70)
    print(f"NULL REPRESENTATIONS ({args.n_null} replicates each)")
    print("=" * 70)
    M_h = memberships["hallmark"]
    sizes = M_h.sum(axis=1).astype(int)
    k = M_h.shape[0]
    null_rows = []
    for kind in ("permuted", "resampled", "gaussian"):
        t1 = time.time()
        for i in range(args.n_null):
            rng = np.random.default_rng([args.seed, hash(kind) % (2**31), i])
            if kind == "permuted":
                W = fal.membership_to_weights(fal.permuted_membership(M_h, rng))
            elif kind == "resampled":
                W = fal.membership_to_weights(fal.resampled_membership(sizes, n_genes, rng))
            else:
                W = fal.gaussian_projection(k, n_genes, rng)
            res = evaluate(W, D, halves, control, folds)
            res["representation"] = kind
            res["replicate"] = i
            null_rows.append(res)
        print(f"  {kind:10s} done [{time.time() - t1:.0f}s]", flush=True)
    nulls = pd.concat(null_rows, ignore_index=True)
    nulls.to_csv(OUT / "null_replicates.csv", index=False)

    # ------------------------------------------------------------ summaries
    print("\n" + "=" * 70)
    print("HALLMARK VS NULLS")
    print("=" * 70)
    summary = []
    for metric in ("r_gamma", "r_gamma_normalised"):
        print(f"\n--- {metric} ---")
        for cl in ORDER:
            obs = observed["hallmark"].loc[observed["hallmark"].cell_line == cl, metric].iloc[0]
            line = f"  {cl:8s} Hallmark={obs:+.4f}"
            for kind in ("permuted", "resampled", "gaussian"):
                vals = nulls[(nulls.cell_line == cl) & (nulls.representation == kind)][metric]
                st = fal.empirical_percentile(obs, vals.to_numpy())
                summary.append(
                    {"metric": metric, "cell_line": cl, "null": kind, "observed": obs, **st}
                )
                line += (
                    f" | {kind[:4]}: {st['null_mean']:+.4f}+-{st['null_sd']:.4f}"
                    f" p={st['p_value']:.3f} z={st['z']:+.1f}"
                )
            print(line)
    pd.DataFrame(summary).to_csv(OUT / "null_comparison.csv", index=False)

    # Reactome against its own matched permuted null
    print("\n" + "=" * 70)
    print("REACTOME VS ITS OWN PERMUTED NULL")
    print("=" * 70)
    M_r = memberships["reactome"]
    react_null = []
    n_react_null = max(20, args.n_null // 5)
    for i in range(n_react_null):
        rng = np.random.default_rng([args.seed, 777, i])
        res = evaluate(
            fal.membership_to_weights(fal.permuted_membership(M_r, rng)), D, halves, control, folds
        )
        res["replicate"] = i
        react_null.append(res)
    react_null = pd.concat(react_null, ignore_index=True)
    react_null.to_csv(OUT / "reactome_null.csv", index=False)
    react_summary = []
    for metric in ("r_gamma", "r_gamma_normalised"):
        for cl in ORDER:
            obs = observed["reactome"].loc[observed["reactome"].cell_line == cl, metric].iloc[0]
            st = fal.empirical_percentile(
                obs, react_null[react_null.cell_line == cl][metric].to_numpy()
            )
            react_summary.append({"metric": metric, "cell_line": cl, "observed": obs, **st})
            if metric == "r_gamma_normalised":
                print(
                    f"  {cl:8s} Reactome={obs:+.4f}  null={st['null_mean']:+.4f}"
                    f"+-{st['null_sd']:.4f}  p={st['p_value']:.3f}  z={st['z']:+.1f}"
                )
    pd.DataFrame(react_summary).to_csv(OUT / "reactome_comparison.csv", index=False)

    # ------------------------------------------------- E. dependence test
    print("\n" + "=" * 70)
    print("E. K562 / JURKAT DEPENDENCE (sensitivity diagnostic, not a benchmark)")
    print("=" * 70)
    idx = {c: i for i, c in enumerate(scperteval.CONTEXTS)}
    restricted = {
        "replogle22k562": [idx["replogle22rpe1"], idx["nadig25hepg2"]],
        "nadig25jurkat": [idx["replogle22rpe1"], idx["nadig25hepg2"]],
        "replogle22rpe1": [idx["replogle22k562"], idx["nadig25hepg2"]],
        "nadig25hepg2": [idx["replogle22k562"], idx["replogle22rpe1"]],
    }
    dep_rows = []
    for label in ("hallmark", "reactome"):
        W = fal.membership_to_weights(memberships[label])
        res = evaluate(W, D, halves, control, folds, source_override=restricted)
        res["representation"] = label
        dep_rows.append(res)
        print(f"\n  {label}:")
        for _, r in res.iterrows():
            full = observed[label].loc[observed[label].cell_line == r.cell_line]
            drop = r.r_gamma - float(full.r_gamma.iloc[0])
            print(
                f"    {r.cell_line:8s} sources={r.sources:22s} "
                f"r_gamma={r.r_gamma:+.4f} (3-source={float(full.r_gamma.iloc[0]):+.4f}, "
                f"change={drop:+.4f})"
            )
    dependence = pd.concat(dep_rows, ignore_index=True)
    dependence.to_csv(OUT / "dependence_test.csv", index=False)

    allobs = pd.concat([*observed.values(), gene_level], ignore_index=True)
    allobs.to_csv(OUT / "observed_representations.csv", index=False)
    (OUT / "summary.json").write_text(
        json.dumps(
            {
                "generated": date.today().isoformat(),
                "runtime_seconds": time.time() - t0,
                "n_null": args.n_null,
                "n_reactome_null": n_react_null,
                "seed": args.seed,
                "baseline": BASELINE,
                "observed": allobs.to_dict("records"),
                "null_comparison": summary,
                "reactome_comparison": react_summary,
                "dependence": dependence.to_dict("records"),
            },
            indent=2,
            default=str,
        )
    )
    print(f"\nTotal runtime: {(time.time() - t0) / 60:.1f} min")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
