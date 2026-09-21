"""Building gene-level feature matrices from the catalogued prior sources.

Every builder takes the gene symbols it must describe and returns a dense
``(n_genes, n_features)`` matrix together with a boolean ``covered`` vector
saying which of those genes the source actually knows about. Coverage is
returned rather than silently imputed: a gene absent from STRING and a gene
present in STRING with no strong partner are different situations, and a
benchmark that cannot tell them apart will misattribute its own failures.

Dimensionality reduction is deliberately **not** done here. The unseen-
perturbation protocol requires any projection to be fitted on training
perturbations alone, so builders emit the raw representation and
:func:`fit_pca` is applied inside the cross-validation loop.
"""

from __future__ import annotations

import gzip
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from virtual_cell.analysis.foundations import read_gmt

__all__ = [
    "FeatureBlock",
    "pathway_features",
    "string_features",
    "depmap_features",
    "basal_features",
    "coexpression_features",
    "fit_pca",
    "apply_pca",
]


@dataclass(frozen=True)
class FeatureBlock:
    """A gene-level feature matrix and the coverage that produced it."""

    name: str
    genes: pd.Index
    matrix: np.ndarray
    covered: np.ndarray
    feature_names: tuple[str, ...] | None = None

    def __post_init__(self) -> None:
        n = len(self.genes)
        if self.matrix.shape[0] != n or self.covered.shape != (n,):
            raise ValueError("matrix and covered must have one row per gene")

    @property
    def n_features(self) -> int:
        return int(self.matrix.shape[1])

    @property
    def coverage(self) -> float:
        return float(self.covered.mean()) if len(self.genes) else 0.0

    def subset(self, genes: Sequence[str]) -> FeatureBlock:
        idx = self.genes.get_indexer(pd.Index(genes).astype(object))
        if (idx < 0).any():
            missing = [g for g, i in zip(genes, idx, strict=True) if i < 0]
            raise KeyError(f"{len(missing)} genes are not in this block, e.g. {missing[:5]}")
        return FeatureBlock(
            name=self.name,
            genes=pd.Index(genes).astype(object),
            matrix=self.matrix[idx],
            covered=self.covered[idx],
            feature_names=self.feature_names,
        )


def _index(genes: Sequence[str]) -> pd.Index:
    return pd.Index(genes).astype(object)


# --------------------------------------------------------------------------
# pathway membership
# --------------------------------------------------------------------------


def pathway_features(
    genes: Sequence[str], gmt_path: str | Path, *, name: str, min_size: int = 5
) -> FeatureBlock:
    """Binary pathway membership, one column per gene set.

    Sets smaller than ``min_size`` *within the supplied gene list* are dropped:
    a set contributing one or zero genes is a constant column that carries no
    information and inflates the feature count.
    """
    sets = read_gmt(gmt_path)
    index = _index(genes)
    pos = pd.Series(np.arange(len(index)), index=index)

    columns: list[np.ndarray] = []
    labels: list[str] = []
    for set_name, members in sorted(sets.items()):
        hit = pos.reindex(_index(members)).dropna().astype(int).to_numpy()
        if hit.size < min_size:
            continue
        col = np.zeros(len(index), dtype=np.float64)
        col[hit] = 1.0
        columns.append(col)
        labels.append(set_name)

    matrix = np.column_stack(columns) if columns else np.zeros((len(index), 0))
    covered = matrix.sum(axis=1) > 0
    return FeatureBlock(name, index, matrix, covered, tuple(labels))


# --------------------------------------------------------------------------
# STRING association network
# --------------------------------------------------------------------------


def _string_alias(info_path: str | Path) -> pd.Series:
    """Ensembl protein id -> preferred gene name."""
    with gzip.open(info_path, "rt") as fh:
        info = pd.read_csv(fh, sep="\t", usecols=[0, 1])
    info.columns = ["protein", "gene"]
    return pd.Series(info["gene"].astype(str).to_numpy(), index=info["protein"].astype(str))


def string_features(
    genes: Sequence[str],
    links_path: str | Path,
    info_path: str | Path,
    *,
    min_score: int = 700,
    name: str = "string",
) -> FeatureBlock:
    """A symmetric association profile over the supplied genes.

    Row ``i`` holds the STRING combined score between gene ``i`` and every gene
    in the list, scaled to ``[0, 1]``. ``min_score`` is STRING's own high-
    confidence threshold; lower-scoring edges are dropped rather than
    down-weighted, because the score is a confidence, not a strength.

    The profile is over the supplied genes only, so it describes a gene by the
    company it keeps *within the panel under study* -- which is what a
    neighbour-based predictor needs.
    """
    index = _index(genes)
    alias = _string_alias(info_path)
    pos = pd.Series(np.arange(len(index)), index=index)

    matrix = np.zeros((len(index), len(index)), dtype=np.float32)
    with gzip.open(links_path, "rt") as fh:
        reader = pd.read_csv(fh, sep=" ", chunksize=2_000_000)
        for chunk in reader:
            chunk = chunk[chunk["combined_score"] >= min_score]
            if chunk.empty:
                continue
            g1 = alias.reindex(chunk["protein1"].to_numpy())
            g2 = alias.reindex(chunk["protein2"].to_numpy())
            i = pos.reindex(g1.to_numpy()).to_numpy()
            j = pos.reindex(g2.to_numpy()).to_numpy()
            keep = ~(np.isnan(i) | np.isnan(j))
            if not keep.any():
                continue
            score = chunk["combined_score"].to_numpy()[keep] / 1000.0
            matrix[i[keep].astype(int), j[keep].astype(int)] = score

    matrix = np.maximum(matrix, matrix.T)
    np.fill_diagonal(matrix, 0.0)
    covered = matrix.sum(axis=1) > 0
    return FeatureBlock(name, index, matrix.astype(np.float64), covered, tuple(index))


# --------------------------------------------------------------------------
# DepMap gene effect
# --------------------------------------------------------------------------


def depmap_features(
    genes: Sequence[str],
    gene_effect_path: str | Path,
    *,
    name: str = "depmap",
    summaries_only: bool = True,
) -> FeatureBlock:
    """Gene-effect features for each gene, pooled across all DepMap lines.

    ``summaries_only`` is the default and is a **leakage control, not a
    convenience**: it emits distributional summaries over cell lines (mean,
    standard deviation, quantiles, fraction of lines where the gene is
    selectively essential) instead of the per-line profile. A per-line column
    would let a model key on a specific cell line, and since Arc's contexts are
    unidentified, selecting lines by similarity to them would be an identity
    inference. See ``data/provenance/depmap/depmap.md``.
    """
    if not summaries_only:
        raise NotImplementedError(
            "per-line DepMap profiles are withheld deliberately; see the leakage "
            "note in data/provenance/depmap/depmap.md"
        )

    header = pd.read_csv(gene_effect_path, nrows=0)
    raw = list(header.columns[1:])
    # DepMap column labels are "SYMBOL (entrez)".
    symbol = {col.split(" (")[0]: col for col in raw}

    index = _index(genes)
    wanted = [symbol[g] for g in index if g in symbol]
    frame = pd.read_csv(gene_effect_path, usecols=[header.columns[0], *wanted], low_memory=False)
    frame = frame.set_index(frame.columns[0])
    frame.columns = [c.split(" (")[0] for c in frame.columns]

    labels = ("mean", "std", "q10", "q50", "q90", "frac_dependent", "n_lines")
    matrix = np.full((len(index), len(labels)), np.nan)
    covered = np.zeros(len(index), dtype=bool)
    for i, gene in enumerate(index):
        if gene not in frame.columns:
            continue
        col = frame[gene].to_numpy(dtype=np.float64)
        col = col[np.isfinite(col)]
        if col.size == 0:
            continue
        covered[i] = True
        matrix[i] = (
            col.mean(),
            col.std(),
            np.quantile(col, 0.10),
            np.quantile(col, 0.50),
            np.quantile(col, 0.90),
            float((col < -0.5).mean()),
            float(col.size),
        )
    return FeatureBlock(name, index, matrix, covered, labels)


# --------------------------------------------------------------------------
# features read from the control cells themselves
# --------------------------------------------------------------------------


def basal_features(
    genes: Sequence[str],
    control_means: np.ndarray,
    gene_axis: Sequence[str],
    *,
    context_names: Sequence[str] | None = None,
    name: str = "basal",
) -> FeatureBlock:
    """Each gene's own basal expression across contexts, plus simple summaries.

    ``control_means`` is ``(n_contexts, n_genes_on_axis)``. Reads control cells
    only, so it carries no perturbation response.
    """
    index = _index(genes)
    axis = _index(gene_axis)
    control_means = np.asarray(control_means, dtype=np.float64)
    if control_means.shape[1] != len(axis):
        raise ValueError("control_means must have one column per gene on gene_axis")

    n_ctx = control_means.shape[0]
    names = tuple(f"basal_{c}" for c in (context_names or range(n_ctx))) + (
        "basal_mean",
        "basal_std",
        "basal_range",
    )
    loc = axis.get_indexer(index)
    matrix = np.full((len(index), len(names)), np.nan)
    covered = loc >= 0
    if covered.any():
        block = control_means[:, loc[covered]].T
        matrix[covered] = np.column_stack(
            [block, block.mean(axis=1), block.std(axis=1), block.max(axis=1) - block.min(axis=1)]
        )
    return FeatureBlock(name, index, matrix, covered, names)


def coexpression_features(
    genes: Sequence[str],
    control_profiles: np.ndarray,
    gene_axis: Sequence[str],
    *,
    name: str = "coexpression",
) -> FeatureBlock:
    """Correlation of each gene's basal profile with every gene on the axis.

    ``control_profiles`` is ``(n_samples, n_genes_on_axis)`` of control-cell
    expression -- pseudobulk per context, or per replicate. This is a prior over
    gene relatedness computed entirely from unperturbed cells.
    """
    index = _index(genes)
    axis = _index(gene_axis)
    x = np.asarray(control_profiles, dtype=np.float64)
    if x.shape[1] != len(axis):
        raise ValueError("control_profiles must have one column per gene on gene_axis")

    centred = x - x.mean(axis=0, keepdims=True)
    norm = np.linalg.norm(centred, axis=0)
    safe = np.where(norm > 0, norm, 1.0)
    unit = centred / safe

    loc = axis.get_indexer(index)
    covered = (loc >= 0) & np.isin(loc, np.flatnonzero(norm > 0))
    matrix = np.zeros((len(index), len(axis)))
    if covered.any():
        matrix[covered] = unit[:, loc[covered]].T @ unit
    return FeatureBlock(name, index, matrix, covered, tuple(axis))


# --------------------------------------------------------------------------
# projection, fitted on training rows only
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class PCA:
    """A projection fitted on a specific set of rows."""

    mean: np.ndarray
    components: np.ndarray
    explained_variance_ratio: np.ndarray


def fit_pca(matrix: np.ndarray, n_components: int) -> PCA:
    """Fit a PCA on ``matrix`` alone. Callers pass *training* rows only."""
    x = np.asarray(matrix, dtype=np.float64)
    mean = x.mean(axis=0)
    centred = x - mean
    n_components = int(min(n_components, *centred.shape))
    if n_components < 1:
        raise ValueError("n_components must be at least 1")
    _, s, vt = np.linalg.svd(centred, full_matrices=False)
    var = s**2
    total = var.sum()
    ratio = var[:n_components] / total if total > 0 else np.zeros(n_components)
    return PCA(mean=mean, components=vt[:n_components], explained_variance_ratio=ratio)


def apply_pca(matrix: np.ndarray, pca: PCA) -> np.ndarray:
    return (np.asarray(matrix, dtype=np.float64) - pca.mean) @ pca.components.T
