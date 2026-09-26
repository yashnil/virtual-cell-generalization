r"""Split-half measurement quality of a public perturbation source.

Machinery for the Kaden source-reliability diagnostic
(``reports/kaden_source_reliability_diagnostic_v1.md``). Nothing here fits a
model or touches the frozen Arc predictor; it measures how reproducible a
dataset's own responses are, per perturbation and averaged into a main effect.

Block split-halves
------------------
A cell-level split half needs one pass over the count matrix per repeat, and
the largest source (``kaden25rpe1``, 2.45e9 nonzeros) does not fit in memory.
So every group's cells are randomly dealt, once, into ``n_blocks`` blocks of
near-equal size (cell ``k`` of a random permutation goes to block
``k mod n_blocks``), block sums are accumulated in **one** streaming pass, and
each repeat draws a uniformly random half of each group's blocks. Each half is
still a uniformly random subset of the group's cells at block granularity, so
expectations match a cell-level split; only the repeat-to-repeat variation is
slightly reduced. Controls are blocked the same way and split independently in
every repeat, so the two halves carry independent control noise -- the
protocol of :func:`virtual_cell.modelling.external_benchmark.target_reliability`.

Reliability of a half versus the full estimate
----------------------------------------------
``corr(half_a, half_b)`` is the reliability of one half; the full estimate's
reliability is its Spearman-Brown projection ``2 r / (1 + r)``
(:func:`virtual_cell.analysis.loco.spearman_brown`). Two independent datasets
measuring the same latent response can correlate at most at
``sqrt(R_x R_y)`` (their *noise ceiling*); an observed correlation far below a
high ceiling is a genuine disagreement, not noise.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import sparse

from virtual_cell.analysis.loco import spearman_brown
from virtual_cell.data import arc2026, scperteval

__all__ = [
    "BlockLayout",
    "assign_blocks",
    "accumulate_block_sums",
    "draw_half_masks",
    "half_means",
    "control_half_means",
    "row_pearson",
    "row_cosine",
    "pearson",
    "cosine",
    "norm_ratio",
    "sign_agreement_top",
    "noise_ceiling",
    "disattenuated",
    "quality_band",
    "spearman_brown",
]


@dataclass(frozen=True)
class BlockLayout:
    """Where every cell's contribution lands in the block-sum matrix.

    Rows ``g * n_blocks + b`` hold block ``b`` of group ``g`` for the
    ``n_groups`` perturbation (and pseudo-perturbation) groups; the final
    ``n_control_blocks`` rows hold the control blocks.
    """

    n_groups: int
    n_blocks: int
    n_control_blocks: int
    cell_block: np.ndarray  # (n_cells,) block row per cell, -1 when excluded

    @property
    def n_rows(self) -> int:
        return self.n_groups * self.n_blocks + self.n_control_blocks

    @property
    def control_offset(self) -> int:
        return self.n_groups * self.n_blocks


def assign_blocks(
    group: np.ndarray,
    *,
    n_groups: int,
    control_code: int,
    n_blocks: int,
    n_control_blocks: int,
    rng: np.random.Generator,
) -> BlockLayout:
    """Deal each group's cells into near-equal random blocks.

    ``group[i]`` is ``0..n_groups-1`` for a perturbation group, ``control_code``
    for a control cell, and any negative value for a cell to ignore.
    """
    group = np.asarray(group, dtype=np.int64)
    cell_block = np.full(group.shape, -1, dtype=np.int64)
    for g in range(n_groups):
        rows = np.flatnonzero(group == g)
        if len(rows):
            perm = rng.permutation(rows)
            cell_block[perm] = g * n_blocks + (np.arange(len(perm)) % n_blocks)
    rows = np.flatnonzero(group == control_code)
    if len(rows):
        perm = rng.permutation(rows)
        cell_block[perm] = n_groups * n_blocks + (np.arange(len(perm)) % n_control_blocks)
    return BlockLayout(n_groups, n_blocks, n_control_blocks, cell_block)


def accumulate_block_sums(
    path: str | Path,
    *,
    genes: Sequence[str],
    layout: BlockLayout,
    chunk_size: int = 20_000,
) -> tuple[np.ndarray, np.ndarray]:
    """One streaming pass: per-block expression sums and cell counts.

    Returns ``(sums, counts)`` with shapes ``(layout.n_rows, len(genes))`` and
    ``(layout.n_rows,)``. The count matrix is never materialised.
    """
    path = Path(path)
    selector = scperteval._column_selector(scperteval.read_var_names(path), genes)
    sums = np.zeros((layout.n_rows, len(genes)), dtype=np.float64)
    counts = np.zeros(layout.n_rows, dtype=np.int64)
    for start, chunk in arc2026.stream_row_chunks(path, chunk_size=chunk_size):
        blocks = layout.cell_block[start : start + chunk.shape[0]]
        keep = np.flatnonzero(blocks >= 0)
        if not len(keep):
            continue
        x = (chunk[keep] @ selector).astype(np.float64)
        indicator = sparse.csr_matrix(
            (np.ones(len(keep)), (blocks[keep], np.arange(len(keep)))),
            shape=(layout.n_rows, len(keep)),
        )
        prod = (indicator @ x).tocoo()
        prod.sum_duplicates()
        sums[prod.row, prod.col] += prod.data
        counts += np.bincount(blocks[keep], minlength=layout.n_rows)
    return sums, counts


def draw_half_masks(n_groups: int, n_blocks: int, rng: np.random.Generator) -> np.ndarray:
    """``(n_groups, n_blocks)`` boolean masks, each row exactly half True."""
    if n_blocks % 2:
        raise ValueError("n_blocks must be even")
    keys = rng.random((n_groups, n_blocks))
    order = np.argsort(keys, axis=1)
    masks = np.zeros((n_groups, n_blocks), dtype=bool)
    np.put_along_axis(masks, order[:, : n_blocks // 2], True, axis=1)
    return masks


def half_means(
    block_sums: np.ndarray, block_counts: np.ndarray, masks: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Mean profile of each half from ``(n_groups, n_blocks, G)`` block sums.

    Returns ``(mean_a, mean_b, n_a, n_b)``; a half with no cells gives NaN.
    """
    m = masks.astype(np.float64)
    sum_a = np.matmul(m[:, None, :], block_sums)[:, 0, :]
    sum_b = block_sums.sum(axis=1) - sum_a
    n_a = (m * block_counts).sum(axis=1)
    n_b = block_counts.sum(axis=1) - n_a
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_a = sum_a / n_a[:, None]
        mean_b = sum_b / n_b[:, None]
    return mean_a, mean_b, n_a, n_b


def control_half_means(
    control_sums: np.ndarray, control_counts: np.ndarray, rng: np.random.Generator
) -> tuple[np.ndarray, np.ndarray]:
    """Independent control halves: a random half of the control blocks each."""
    mask = draw_half_masks(1, control_sums.shape[0], rng)[0]
    a = control_sums[mask].sum(axis=0) / control_counts[mask].sum()
    b = control_sums[~mask].sum(axis=0) / control_counts[~mask].sum()
    return a, b


# --------------------------------------------------------------------------
# agreement statistics
# --------------------------------------------------------------------------


def row_pearson(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Row-wise Pearson correlation across columns; NaN where undefined."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    ac = a - a.mean(axis=1, keepdims=True)
    bc = b - b.mean(axis=1, keepdims=True)
    denom = np.sqrt((ac * ac).sum(axis=1) * (bc * bc).sum(axis=1))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(denom > 0, (ac * bc).sum(axis=1) / denom, np.nan)


def row_cosine(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Row-wise uncentred cosine; NaN where either row is zero."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    denom = np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(denom > 0, (a * b).sum(axis=1) / denom, np.nan)


def pearson(a: np.ndarray, b: np.ndarray) -> float:
    return float(row_pearson(np.atleast_2d(a), np.atleast_2d(b))[0])


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    return float(row_cosine(np.atleast_2d(a), np.atleast_2d(b))[0])


def norm_ratio(a: np.ndarray, b: np.ndarray) -> float:
    """``||a|| / ||b||``."""
    nb = float(np.linalg.norm(b))
    return float(np.linalg.norm(a)) / nb if nb > 0 else float("nan")


def sign_agreement_top(a: np.ndarray, b: np.ndarray, k: int) -> np.ndarray:
    """Per row: fraction of genes with matching sign, over the union of each
    row's top-``k`` genes by absolute value in ``a`` and in ``b``."""
    a = np.atleast_2d(np.asarray(a, dtype=np.float64))
    b = np.atleast_2d(np.asarray(b, dtype=np.float64))
    k = min(k, a.shape[1])
    out = np.full(a.shape[0], np.nan)
    for i in range(a.shape[0]):
        if not (np.isfinite(a[i]).all() and np.isfinite(b[i]).all()):
            continue
        top = np.union1d(
            np.argpartition(-np.abs(a[i]), k - 1)[:k], np.argpartition(-np.abs(b[i]), k - 1)[:k]
        )
        out[i] = float(np.mean(np.sign(a[i, top]) == np.sign(b[i, top])))
    return out


def noise_ceiling(rel_x: np.ndarray | float, rel_y: np.ndarray | float) -> np.ndarray | float:
    """Largest correlation two measurements of one latent response can reach:
    ``sqrt(R_x R_y)`` with negative reliabilities clipped to 0."""
    rx = np.clip(np.asarray(rel_x, dtype=np.float64), 0.0, None)
    ry = np.clip(np.asarray(rel_y, dtype=np.float64), 0.0, None)
    return np.sqrt(rx * ry)


def disattenuated(
    r_obs: np.ndarray, rel_x: np.ndarray, rel_y: np.ndarray, *, min_reliability: float
) -> np.ndarray:
    """``r_obs / sqrt(R_x R_y)``, NaN unless both reliabilities reach the floor.

    A model-based correction, not ground truth: it assumes both datasets
    measure the same latent response plus independent noise.
    """
    r = np.asarray(r_obs, dtype=np.float64)
    rx = np.asarray(rel_x, dtype=np.float64)
    ry = np.asarray(rel_y, dtype=np.float64)
    ok = (rx >= min_reliability) & (ry >= min_reliability) & np.isfinite(r)
    out = np.full(r.shape, np.nan)
    out[ok] = r[ok] / np.sqrt(rx[ok] * ry[ok])
    return out


def quality_band(reliability: float, *, high: float, moderate: float) -> str:
    """Descriptive band for a Spearman-Brown reliability; never a filter."""
    if not np.isfinite(reliability):
        return "not estimable"
    if reliability >= high:
        return "high"
    if reliability >= moderate:
        return "moderate"
    return "low / uninformative"
