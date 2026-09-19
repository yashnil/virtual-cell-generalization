"""Predeclared robustness battery for the independent four-context decomposition.

Runs variants A-E of the frozen sensitivity plan against the frozen
1,264 x 6,640 design. Canonical v1 outputs are never touched; everything lands in
``outputs/four_context_sensitivity/``.

  A  independent control split  --control-scheme split
  B  controlled cell-depth      --depth-experiment
  C  feature space              gene subsets, control-derived HVG ranking
  D  aggregation order          --aggregation log_mean
  E  seed stability             --seed

No option is tuned toward any published percentage, and no predictive model is
built.

Usage::

    uv run python scripts/run_four_context_sensitivity.py --all
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from virtual_cell.analysis import robustness
from virtual_cell.data import scperteval
from virtual_cell.decomposition import anova

REPO_ROOT = Path(__file__).resolve().parents[1]
DESIGN_DIR = REPO_ROOT / "data" / "splits" / "four_context_v1"
CANONICAL_DIR = REPO_ROOT / "data" / "processed" / "four_context_v1"
OUT_DIR = REPO_ROOT / "outputs" / "four_context_sensitivity"
WORK_DIR = REPO_ROOT / "data" / "processed" / "four_context_sensitivity"

GENE_SUBSETS = {"all": None, "hvg4000": 4000, "hvg2000": 2000}
DEPTHS = (15, 30, 50, 100)


def load_design() -> tuple[list[str], list[str]]:
    perts = (DESIGN_DIR / "shared_perturbations.txt").read_text().split()
    genes = (DESIGN_DIR / "shared_genes.txt").read_text().split()
    return perts, genes


def hvg_ranking(data_dir: Path, genes: list[str]) -> pd.Series:
    """Control-derived global HVG ranking, cached (basal expression only)."""
    cache = WORK_DIR / "control_gene_variance.csv"
    if cache.exists():
        frame = pd.read_csv(cache, index_col=0)
        if list(frame.index) == genes:
            return robustness.global_hvg_ranking({c: frame[c].to_numpy() for c in frame}, genes)
    variances = {}
    for ds in scperteval.DATASETS:
        t0 = time.time()
        variances[ds.name] = robustness.control_gene_variance(ds.path(data_dir), genes=genes)
        print(f"    control variance {ds.name:16s} [{time.time() - t0:.1f}s]", flush=True)
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(variances, index=pd.Index(genes, name="gene")).to_csv(cache)
    return robustness.global_hvg_ranking(variances, genes)


def subset_columns(genes: list[str], ranking: pd.Series, k: int | None) -> np.ndarray:
    """Column indices of a gene subset, kept in the frozen gene order."""
    if k is None:
        return np.arange(len(genes))
    chosen = set(ranking.index[:k])
    return np.array([i for i, g in enumerate(genes) if g in chosen])


def decompose_subset(D: np.ndarray, cols: np.ndarray, perts: list[str]) -> dict:
    sub = D[:, :, cols]
    dec = anova.decompose(sub, cell_lines=scperteval.CONTEXTS, perturbations=perts)
    total = anova.total_sum_of_squares(sub)
    ss = anova.sums_of_squares(dec)
    fractions = {k: v / total for k, v in ss.items()}
    fractions["template"] = fractions["mu"] + fractions["alpha"]
    return {
        "n_genes": int(len(cols)),
        "total_ss": total,
        "sums_of_squares": ss,
        "uncorrected_fractions": fractions,
        "reconstruction_max_abs_err": float(np.abs(dec.reconstruct() - sub).max()),
        "ss_partition_rel_err": float(abs(sum(ss.values()) - total) / total),
        "zero_sum_max": float(
            max(
                np.abs(dec.alpha.sum(axis=0)).max(),
                np.abs(dec.beta.sum(axis=0)).max(),
                np.abs(dec.gamma.sum(axis=0)).max(),
                np.abs(dec.gamma.sum(axis=1)).max(),
            )
        ),
    }


def run_variant(
    data_dir: Path,
    perts: list[str],
    genes: list[str],
    ranking: pd.Series,
    *,
    control_scheme: str,
    aggregation: str,
    seed: int,
    n_splits: int,
) -> dict:
    """One decomposition + split-half configuration, evaluated on all gene subsets."""
    tag = f"{control_scheme}_{aggregation}_seed{seed}_n{n_splits}"
    print(f"\n--- variant {tag} ---", flush=True)
    n_c, n_p, n_g = len(scperteval.CONTEXTS), len(perts), len(genes)
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    mm_path = WORK_DIR / f"_half_{tag}.npy"
    half = np.lib.format.open_memmap(
        mm_path, mode="w+", dtype=np.float32, shape=(n_splits, 2, n_c, n_p, n_g)
    )

    full_delta = np.zeros((n_c, n_p, n_g))
    for ci, ds in enumerate(scperteval.DATASETS):
        t0 = time.time()
        X, group, X_ctrl = robustness.load_cells(
            ds.path(data_dir), genes=genes, perturbations=perts
        )
        ctrl_profile = robustness.aggregate(X_ctrl, aggregation)
        full_delta[ci] = (
            robustness.group_aggregate(X, group, n_p, aggregation) - ctrl_profile[None, :]
        )
        print(
            f"  {ds.name:16s} cells={X.shape[0]:,} ctrl={X_ctrl.shape[0]:,} "
            f"[{time.time() - t0:.0f}s load]",
            flush=True,
        )
        t0 = time.time()
        for s in range(n_splits):
            rng = np.random.default_rng([seed, ci, s])
            a, b = robustness.half_deltas(
                X,
                group,
                X_ctrl,
                n_perturbations=n_p,
                rng=rng,
                control_scheme=control_scheme,
                aggregation=aggregation,
            )
            half[s, 0, ci] = a.astype(np.float32)
            half[s, 1, ci] = b.astype(np.float32)
        print(f"    {n_splits} splits [{time.time() - t0:.0f}s]", flush=True)
        del X, X_ctrl
        half.flush()

    result = {
        "tag": tag,
        "control_scheme": control_scheme,
        "aggregation": aggregation,
        "seed": seed,
        "n_splits": n_splits,
        "subsets": {},
    }

    for name, k in GENE_SUBSETS.items():
        cols = subset_columns(genes, ranking, k)
        dec_info = decompose_subset(full_delta, cols, perts)
        per_split = []
        for s in range(n_splits):
            D1 = np.asarray(half[s, 0][:, :, cols], dtype=np.float64)
            D2 = np.asarray(half[s, 1][:, :, cols], dtype=np.float64)
            per_split.append(anova.cross_half_signal(D1, D2))
        signal = {k2: float(np.mean([p[k2] for p in per_split])) for k2 in per_split[0]}
        corrected = anova.noise_corrected_fractions(signal, dec_info["total_ss"])
        per_split_fracs = pd.DataFrame(
            [anova.noise_corrected_fractions(p, dec_info["total_ss"]) for p in per_split]
        )
        reproducibility = {
            k2: signal[k2] / dec_info["sums_of_squares"][k2]
            for k2 in ("mu", "alpha", "beta", "gamma")
        }
        result["subsets"][name] = {
            **dec_info,
            "signal": signal,
            "corrected_fractions": corrected,
            "corrected_sd": {c: float(per_split_fracs[c].std()) for c in per_split_fracs},
            "corrected_min": {c: float(per_split_fracs[c].min()) for c in per_split_fracs},
            "corrected_max": {c: float(per_split_fracs[c].max()) for c in per_split_fracs},
            "reproducibility": reproducibility,
        }
        print(
            f"  [{name:8s} n_genes={dec_info['n_genes']:>5}] "
            f"template={corrected['template'] * 100:5.2f}% beta={corrected['beta'] * 100:5.2f}% "
            f"gamma={corrected['gamma'] * 100:5.2f}% noise={corrected['noise'] * 100:5.2f}%  "
            f"| repro beta={reproducibility['beta'] * 100:.1f}% "
            f"gamma={reproducibility['gamma'] * 100:.1f}%",
            flush=True,
        )

    del half
    mm_path.unlink()
    return result


def run_depth_experiment(
    data_dir: Path, perts: list[str], genes: list[str], *, seed: int, n_repeats: int
) -> pd.DataFrame:
    print("\n--- B: controlled cell-depth experiment ---", flush=True)
    counts = np.load(CANONICAL_DIR / "cell_counts.npy")
    need = 2 * max(DEPTHS)
    frames = []
    for ci, ds in enumerate(scperteval.DATASETS):
        rows = np.flatnonzero(counts[ci] >= need)
        print(f"  {ds.name:16s} pairs with >= {need} cells: {len(rows)}", flush=True)
        if not len(rows):
            continue
        t0 = time.time()
        X, group, X_ctrl = robustness.load_cells(
            ds.path(data_dir), genes=genes, perturbations=perts
        )
        frame = robustness.depth_reliability(
            X,
            group,
            X_ctrl,
            pair_rows=rows,
            depths=DEPTHS,
            n_repeats=n_repeats,
            rng=np.random.default_rng([seed, ci]),
        )
        frame["context"] = ds.name
        frame["cell_line"] = ds.cell_line
        frame["perturbation"] = [perts[r] for r in frame["pert_row"]]
        frames.append(frame)
        del X, X_ctrl
        print(f"    {len(frame)} (pair, depth) estimates [{time.time() - t0:.0f}s]", flush=True)
    return pd.concat(frames, ignore_index=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=REPO_ROOT / scperteval.DATA_SUBDIR)
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--n-splits", type=int, default=20)
    parser.add_argument("--seeds", type=int, nargs="+", default=[42, 7, 101, 2024, 31337])
    parser.add_argument("--depth-repeats", type=int, default=25)
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--skip-depth", action="store_true")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    perts, genes = load_design()
    print(f"frozen design: {len(perts):,} perturbations x {len(genes):,} genes")

    print("\n--- control-derived HVG ranking (basal expression only) ---")
    ranking = hvg_ranking(args.data_dir, genes)
    ranking.to_csv(args.out_dir / "hvg_ranking.csv")
    print(f"  ranked {len(ranking):,} genes; top 5: {list(ranking.index[:5])}")

    t_start = time.time()
    variants = []

    # E (seed stability) + C (feature space), canonical scheme
    for seed in args.seeds:
        variants.append(
            run_variant(
                args.data_dir,
                perts,
                genes,
                ranking,
                control_scheme="shared",
                aggregation="mean_log",
                seed=seed,
                n_splits=args.n_splits,
            )
        )
    # A: independent control split
    variants.append(
        run_variant(
            args.data_dir,
            perts,
            genes,
            ranking,
            control_scheme="split",
            aggregation="mean_log",
            seed=args.seeds[0],
            n_splits=args.n_splits,
        )
    )
    # D: aggregation order
    variants.append(
        run_variant(
            args.data_dir,
            perts,
            genes,
            ranking,
            control_scheme="shared",
            aggregation="log_mean",
            seed=args.seeds[0],
            n_splits=args.n_splits,
        )
    )

    (args.out_dir / "variants.json").write_text(json.dumps(variants, indent=2, default=str))

    if not args.skip_depth:
        depth = run_depth_experiment(
            args.data_dir, perts, genes, seed=args.seeds[0], n_repeats=args.depth_repeats
        )
        depth.to_csv(args.out_dir / "depth_experiment.csv", index=False)

    print(f"\nTotal runtime: {(time.time() - t_start) / 60:.1f} min")
    print(f"Wrote {args.out_dir}")


if __name__ == "__main__":
    main()
