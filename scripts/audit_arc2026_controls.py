"""Read-only audit of the official Arc Virtual Cell Challenge 2026 control bundle.

Streams the official ``context_{A,B,C}.h5ad`` control matrices, verifies the
manifest and side-car tables, checks every cross-context invariant, and writes
descriptive tables and figures to ``--out-dir``.

The script never writes to ``data/`` and never uses A/B/C as perturbation-response
training data: these contexts contain unperturbed control cells only and are
evaluation inputs.

Usage (from the repository root)::

    uv run python scripts/audit_arc2026_controls.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import sparse  # noqa: E402
from sklearn.decomposition import PCA  # noqa: E402

from virtual_cell.data import arc2026  # noqa: E402
from virtual_cell.preprocessing.pseudobulk import normalized_matrix  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
CONTEXT_COLORS = {"A": "#1f77b4", "B": "#d62728", "C": "#2ca02c"}


def _color(label: str) -> str:
    return CONTEXT_COLORS.get(label, "#777777")


def report_metadata(controls_dir: Path) -> tuple[arc2026.Arc2026Manifest, pd.Index, pd.DataFrame]:
    manifest = arc2026.load_manifest(controls_dir)
    gene_names = arc2026.load_gene_names(controls_dir)
    pert_counts = arc2026.load_pert_counts(controls_dir)

    print("=== manifest.json ===")
    print(json.dumps(dict(manifest.raw), indent=2))

    print("\n=== gene_names.csv ===")
    print(f"rows: {len(gene_names)}   column: {gene_names.name!r}   dtype: {gene_names.dtype}")
    print(f"duplicates: {int(gene_names.duplicated().sum())}")
    print(f"first: {gene_names[:5].tolist()}")
    print(f"last:  {gene_names[-5:].tolist()}")

    print("\n=== pert_counts.csv ===")
    print(f"rows: {len(pert_counts)}   columns: {list(pert_counts.columns)}")
    print(f"dtypes: {pert_counts.dtypes.to_dict()}")
    duplicated = int(pert_counts.duplicated().sum())
    print(f"duplicate rows: {duplicated}")
    for column in pert_counts.columns:
        print(f"  {column}: {pert_counts[column].nunique()} unique values")
    print(f"first: {pert_counts.head(3).to_dict('records')}")

    print("\n=== panel composition (gene families) ===")
    print(arc2026.panel_composition(gene_names).to_string())

    overlap = sorted(set(pert_counts[manifest.pert_col]) & set(gene_names))
    print(
        f"\ntarget genes also present in the 18,533-gene panel: {len(overlap)} / {len(pert_counts)}"
    )
    return manifest, gene_names, pert_counts


def report_context(audit: arc2026.ContextAudit) -> None:
    lib, det = audit.library_sizes, audit.genes_detected
    print(f"\n=== context {audit.context} ({audit.path.name}) ===")
    print(f"shape: ({audit.n_cells}, {audit.n_genes})  cells x genes")
    print(f"X: {audit.x_encoding} (sparse={audit.is_sparse}), dtype={audit.x_dtype}")
    print(
        f"stored entries: {audit.n_stored:,}  non-zero: {audit.n_nonzero:,}  "
        f"explicit zeros: {audit.n_explicit_zeros:,}  sparsity: {audit.sparsity:.4f}"
    )
    print(f"counts: min={audit.min_value:g} max={audit.max_value:g}")
    print(
        f"finite={audit.all_finite} nonnegative={audit.all_nonnegative} "
        f"integer-valued={audit.all_integer}"
    )
    print(f"obs columns: {list(audit.obs_columns)}   var columns: {list(audit.var_columns)}")
    print(
        f"obs_names unique: {audit.obs_names_unique} "
        f"(duplicates: {audit.n_duplicate_obs_names})   "
        f"var_names unique: {audit.var_names_unique} "
        f"(duplicates: {audit.n_duplicate_var_names})"
    )
    print(
        f"library size: min={lib.min():.0f} median={np.median(lib):.0f} "
        f"mean={lib.mean():.1f} max={lib.max():.0f}"
    )
    print(
        f"genes detected: min={det.min()} median={np.median(det):.0f} "
        f"mean={det.mean():.1f} max={det.max()}"
    )
    for column, counts in audit.obs_value_counts.items():
        head = counts.head(3).to_dict()
        print(
            f"obs[{column!r}]: {len(counts)} distinct values, "
            f"counts min={counts.min()} max={counts.max()}, e.g. {head}"
        )


def plot_distributions(audits: dict[str, arc2026.ContextAudit], out_dir: Path) -> list[Path]:
    paths = []

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for label, audit in audits.items():
        axes[0].hist(
            audit.library_sizes, bins=80, histtype="step", lw=1.6, label=label, color=_color(label)
        )
        axes[1].hist(
            np.log10(audit.library_sizes),
            bins=80,
            histtype="step",
            lw=1.6,
            label=label,
            color=_color(label),
        )
    axes[0].set_xlabel("library size (total UMI counts per cell)")
    axes[1].set_xlabel("log10 library size")
    for ax in axes:
        ax.set_ylabel("control cells")
        ax.legend(title="context")
    axes[0].set_title("Library size, official control cells")
    axes[1].set_title("Library size (log10)")
    fig.tight_layout()
    path = out_dir / "library_size_distributions.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    paths.append(path)

    fig, ax = plt.subplots(figsize=(6, 4))
    for label, audit in audits.items():
        ax.hist(
            audit.genes_detected,
            bins=80,
            histtype="step",
            lw=1.6,
            label=label,
            color=_color(label),
        )
    ax.set_xlabel("genes detected per cell (non-zero counts)")
    ax.set_ylabel("control cells")
    ax.set_title("Genes detected, official control cells")
    ax.legend(title="context")
    fig.tight_layout()
    path = out_dir / "genes_detected_distributions.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    paths.append(path)
    return paths


def plot_pseudobulk_scatter(table: pd.DataFrame, correlations: pd.DataFrame, out_dir: Path) -> Path:
    labels = list(table.columns)
    pairs = [(a, b) for i, a in enumerate(labels) for b in labels[i + 1 :]]
    fig, axes = plt.subplots(1, len(pairs), figsize=(4.2 * len(pairs), 4.2))
    axes = np.atleast_1d(axes)
    for ax, (a, b) in zip(axes, pairs, strict=True):
        ax.scatter(table[a], table[b], s=3, alpha=0.2, color="#333333", edgecolors="none")
        lim = [0, float(table[[a, b]].to_numpy().max()) * 1.05]
        ax.plot(lim, lim, "--", lw=0.9, color="#d62728")
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_xlabel(f"basal pseudobulk, context {a}")
        ax.set_ylabel(f"basal pseudobulk, context {b}")
        ax.set_title(f"{a} vs {b}  (Pearson r = {correlations.loc[a, b]:.4f})")
    fig.suptitle("Basal pseudobulk mean log1p(CP10K) expression, 18,533 genes")
    fig.tight_layout()
    path = out_dir / "pseudobulk_scatter.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def plot_control_pca(
    controls_dir: Path,
    contexts: list[str],
    table: pd.DataFrame,
    out_dir: Path,
    *,
    n_cells: int,
    n_genes: int,
    seed: int,
) -> tuple[Path, pd.DataFrame]:
    """PCA of a random subsample of control cells on the most expressed genes."""
    genes = table.mean(axis=1).sort_values(ascending=False).head(n_genes).index
    blocks, labels = [], []
    for context in contexts:
        adata = arc2026.subsample_context(
            arc2026.context_path(controls_dir, context),
            n_cells=n_cells,
            seed=seed,
            genes=genes,
        )
        blocks.append(normalized_matrix(adata))
        labels.extend([context] * adata.n_obs)

    X = sparse.vstack(blocks, format="csr").toarray()
    pca = PCA(n_components=10, svd_solver="randomized", random_state=seed)
    scores = pca.fit_transform(X - X.mean(axis=0, keepdims=True))
    frame = pd.DataFrame(scores[:, :4], columns=[f"PC{i + 1}" for i in range(4)])
    frame.insert(0, "context", labels)

    fig, ax = plt.subplots(figsize=(6, 5))
    for context in contexts:
        mask = frame["context"] == context
        ax.scatter(
            frame.loc[mask, "PC1"],
            frame.loc[mask, "PC2"],
            s=4,
            alpha=0.4,
            label=context,
            color=_color(context),
            edgecolors="none",
        )
    var = pca.explained_variance_ratio_
    ax.set_xlabel(f"PC1 ({var[0] * 100:.1f}% var)")
    ax.set_ylabel(f"PC2 ({var[1] * 100:.1f}% var)")
    ax.set_title(f"PCA of {n_cells} control cells per context\n({n_genes} most expressed genes)")
    ax.legend(title="context", markerscale=3)
    fig.tight_layout()
    path = out_dir / "control_cell_pca.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)

    frame.attrs["explained_variance_ratio"] = var.tolist()
    return path, frame


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--controls-dir", type=Path, default=REPO_ROOT / arc2026.CONTROLS_SUBDIR)
    parser.add_argument(
        "--out-dir", type=Path, default=REPO_ROOT / "outputs" / "arc2026_controls_audit"
    )
    parser.add_argument("--chunk-size", type=int, default=arc2026.DEFAULT_CHUNK_SIZE)
    parser.add_argument("--pca-cells", type=int, default=2000)
    parser.add_argument("--pca-genes", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--no-pca", action="store_true", help="skip the PCA figure")
    parser.add_argument("--top-genes", type=int, default=50)
    args = parser.parse_args()

    manifest, gene_names, pert_counts = report_metadata(args.controls_dir)

    missing = arc2026.missing_bundle_files(args.controls_dir, manifest.contexts)
    if missing:
        raise SystemExit(f"Missing official files: {[str(p) for p in missing]}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    arc2026.panel_composition(gene_names).to_csv(args.out_dir / "panel_composition.csv")

    audits = {}
    for context in manifest.contexts:
        audits[context] = arc2026.audit_context_file(
            arc2026.context_path(args.controls_dir, context),
            expected_context=context,
            context_col=manifest.context_col,
            chunk_size=args.chunk_size,
        )
        report_context(audits[context])

    stats = arc2026.audit_summary_table(audits)
    stats.to_csv(args.out_dir / "context_stats.csv")
    print("\n=== per-context statistics ===")
    with pd.option_context("display.width", 220, "display.max_columns", 40):
        print(stats.T)

    depth = arc2026.depth_quantile_table(audits)
    depth.to_csv(args.out_dir / "depth_quantiles.csv")
    print("\n=== per-cell depth quantiles and low-depth cell counts ===")
    with pd.option_context("display.width", 220, "display.max_columns", 40):
        print(depth.T.round(4))

    checks = arc2026.check_cross_context_invariants(
        audits, manifest=manifest, gene_names=gene_names
    )
    check_frame = pd.DataFrame([vars(c) for c in checks])
    check_frame.to_csv(args.out_dir / "invariant_checks.csv", index=False)
    failures = [c for c in checks if not c.passed]
    print("\n=== cross-context invariants ===")
    for check in checks:
        print(f"[{'PASS' if check.passed else 'FAIL'}] {check.name}  {check.detail}")

    overlaps = arc2026.cell_id_overlaps(audits)
    overlaps.to_csv(args.out_dir / "cell_id_overlap.csv")
    print("\n=== pairwise shared cell identifiers ===")
    print(overlaps)

    if failures:
        raise SystemExit(
            f"\nSTOP: {len(failures)} invariant(s) failed: {[c.name for c in failures]}"
        )

    table = arc2026.basal_mean_table(audits, normalized=True)
    table.to_csv(args.out_dir / "pseudobulk_normalized_mean.csv")
    arc2026.basal_mean_table(audits, normalized=False).to_csv(
        args.out_dir / "pseudobulk_raw_mean.csv"
    )
    pearson = table.corr(method="pearson")
    spearman = table.corr(method="spearman")
    pearson.to_csv(args.out_dir / "pseudobulk_pearson.csv")
    spearman.to_csv(args.out_dir / "pseudobulk_spearman.csv")
    print("\n=== basal pseudobulk correlations (mean log1p CP10K) ===")
    print("Pearson:\n", pearson.round(4))
    print("Spearman:\n", spearman.round(4))

    top = arc2026.top_basal_difference_genes(table, n=args.top_genes)
    top.to_csv(args.out_dir / "top_basal_difference_genes.csv")
    print(f"\n=== top {min(15, args.top_genes)} genes by cross-context spread (descriptive) ===")
    print(top.head(15).round(3))

    figures = plot_distributions(audits, args.out_dir)
    figures.append(plot_pseudobulk_scatter(table, pearson, args.out_dir))
    pca_variance = None
    if not args.no_pca:
        pca_path, pca_frame = plot_control_pca(
            args.controls_dir,
            list(manifest.contexts),
            table,
            args.out_dir,
            n_cells=args.pca_cells,
            n_genes=args.pca_genes,
            seed=args.seed,
        )
        figures.append(pca_path)
        pca_frame.to_csv(args.out_dir / "control_cell_pca_scores.csv", index=False)
        pca_variance = pca_frame.attrs["explained_variance_ratio"]
        print("\n=== control-cell PCA explained variance ratio (first 5) ===")
        print([round(v, 4) for v in pca_variance[:5]])

    summary = {
        "controls_dir": str(args.controls_dir),
        "manifest": dict(manifest.raw),
        "n_genes_csv": int(len(gene_names)),
        "n_perturbations_csv": int(len(pert_counts)),
        "panel_composition": arc2026.panel_composition(gene_names).to_dict(),
        "contexts": {label: audit.stats_row() for label, audit in audits.items()},
        "invariants_passed": len(checks) - len(failures),
        "invariants_total": len(checks),
        "pearson": pearson.round(6).to_dict(),
        "spearman": spearman.round(6).to_dict(),
        "cell_id_overlap": overlaps.to_dict(),
        "depth_quantiles": depth.round(4).to_dict(),
        "pca_explained_variance_ratio": pca_variance,
        "figures": [str(p.relative_to(REPO_ROOT)) for p in figures],
    }
    (args.out_dir / "audit_summary.json").write_text(json.dumps(summary, indent=2, default=str))

    print(f"\nAll {len(checks)} invariants passed.")
    print(f"Wrote tables and {len(figures)} figures to {args.out_dir}")


if __name__ == "__main__":
    main()
