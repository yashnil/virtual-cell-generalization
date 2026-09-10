"""Basic exploration of the synthetic control contexts.

Prints the AnnData structure, per-context summary statistics, library-size
distributions, and a comparison of basal mean expression between contexts.
Writes tables and figures to ``--out-dir`` (default ``outputs/exploration``,
which is git-ignored).

Usage (from the repository root)::

    uv run python scripts/explore_synthetic_contexts.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from virtual_cell.data.io import load_contexts  # noqa: E402
from virtual_cell.data.summary import library_sizes, summarize_contexts  # noqa: E402
from virtual_cell.preprocessing.pseudobulk import (  # noqa: E402
    basal_log_fold_change,
    basal_mean_correlation,
    basal_mean_table,
)

REPO_ROOT = Path(__file__).resolve().parents[1]


def describe_structure(label: str, adata) -> None:
    X = adata.X
    print(f"\n=== Context {label} ===")
    print(adata)
    print(f"rows (cells): {adata.n_obs}   columns (genes): {adata.n_vars}")
    print(f"X type: {type(X).__name__}, dtype: {X.dtype}, nnz: {X.nnz}")
    print(f"first obs_names: {adata.obs_names[:3].tolist()}")
    print(f"first var_names: {adata.var_names[:3].tolist()}")
    print("obs head:")
    print(adata.obs.head(3))
    dense_corner = X[:5, :8].toarray()
    print("count matrix corner (5 cells x 8 genes):")
    print(dense_corner)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=REPO_ROOT / "data" / "raw" / "synthetic")
    parser.add_argument("--out-dir", type=Path, default=REPO_ROOT / "outputs" / "exploration")
    args = parser.parse_args()

    paths = sorted(args.data_dir.glob("context_*.h5ad"))
    if not paths:
        raise SystemExit(
            f"No context_*.h5ad files in {args.data_dir}. "
            "Run scripts/make_synthetic_controls.py first."
        )
    contexts = load_contexts(paths)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    for label, adata in contexts.items():
        describe_structure(label, adata)

    summary = summarize_contexts(contexts)
    print("\n=== Per-context summary ===")
    with pd.option_context("display.width", 200, "display.float_format", "{:.3f}".format):
        print(summary)
    summary.to_csv(args.out_dir / "context_summary.csv")

    means = basal_mean_table(contexts)
    corr = basal_mean_correlation(means)
    print("\n=== Basal mean expression (normalised log1p), first genes ===")
    print(means.head())
    print("\n=== Pearson correlation of basal mean profiles between contexts ===")
    print(corr.round(4))
    means.to_csv(args.out_dir / "basal_means.csv")
    corr.to_csv(args.out_dir / "basal_mean_correlation.csv")

    labels = list(contexts)
    if len(labels) >= 2:
        lfc = basal_log_fold_change(means, labels[0], labels[1])
        print(
            f"\nlog2 fold change {labels[0]} vs {labels[1]}: "
            f"mean={lfc.mean():.3f}, sd={lfc.std():.3f}, "
            f"max |lfc| gene={lfc.abs().idxmax()} ({lfc.abs().max():.3f})"
        )

    # --- figures -------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for label, adata in contexts.items():
        axes[0].hist(library_sizes(adata), bins=30, alpha=0.5, label=label)
    axes[0].set_xlabel("library size (total counts per cell)")
    axes[0].set_ylabel("cells")
    axes[0].set_title("Library size per context")
    axes[0].legend()

    if len(labels) >= 2:
        a, b = labels[0], labels[1]
        axes[1].scatter(means[a], means[b], s=8, alpha=0.7)
        lim = [0, float(np.nanmax(means[[a, b]].to_numpy())) * 1.05]
        axes[1].plot(lim, lim, "k--", lw=0.8)
        axes[1].set_xlabel(f"basal mean expression, context {a}")
        axes[1].set_ylabel(f"basal mean expression, context {b}")
        axes[1].set_title(f"Basal means: {a} vs {b} (r={corr.loc[a, b]:.3f})")
    fig.tight_layout()
    fig_path = args.out_dir / "synthetic_context_overview.png"
    fig.savefig(fig_path, dpi=120)
    print(f"\nWrote tables and figure to {args.out_dir}")


if __name__ == "__main__":
    main()
