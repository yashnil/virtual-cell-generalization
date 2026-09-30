"""Source-centred direct-response effects: gene harmonisation, lookup, centering.

Attribution: the effect spaces and the per-source centering are AtlasShift's
(``model.aligned_effect`` / ``model.aligned_cd4`` @ ``d24ce4f``, MIT). This is our own
implementation of the same definitions; dtypes follow the upstream arithmetic so that
``tests/test_competition_v2_c1.py`` can require machine equality.

Two effect spaces are used:

* ``log2fc``     — ``log2((mean_cpm_pert + 1) / (mean_cpm_ctrl + 1))``, the per-cell
  mean-CPM moment;
* ``bulk_delta`` — ``log1p(5e4 p_pert) - log1p(5e4 p_ctrl)``, the pseudobulk moment in the
  exact space PDS and MSE are computed in.

**Centering** (the frozen definition): for each source, the mean effect over that
source's own usable retained targets, with each target's own gene excluded from its
row, is subtracted from every row. That removes the source's generic perturbation
response and keeps the target-specific part. The centering set never contains a
held-out context's responses, because a held-out source is simply not loaded.
"""

from __future__ import annotations

import numpy as np

from virtual_cell.competition_v2.sources import MINIMUM_CELLS, SourceStats

BULK_SCALE = 50_000.0
SPACES = ("log2fc", "bulk_delta")


def _index(labels: np.ndarray, wanted: np.ndarray) -> np.ndarray:
    lookup = {str(x): i for i, x in enumerate(labels)}
    return np.asarray([lookup.get(str(x), -1) for x in wanted], dtype=np.int64)


def source_effect(
    source: SourceStats,
    targets: np.ndarray,
    genes: np.ndarray,
    *,
    space: str,
    minimum_cells: int = MINIMUM_CELLS,
    pseudocount_cpm: float = 1.0,
    center: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """One source's (optionally centred) effect on a (targets x genes) panel.

    Returns ``(effect, mask)``: ``effect`` is float32 with zeros where the source has no
    usable measurement, ``mask`` marks usable (target, gene) cells. The centering mean is
    taken over the rows of ``targets`` that are usable in this source.
    """
    if space not in SPACES:
        raise ValueError(space)
    targets = np.asarray(targets).astype(str)
    genes = np.asarray(genes).astype(str)
    rows = _index(source.targets, targets)
    cols = _index(source.genes, genes)
    row_ok = rows >= 0
    row_ok[row_ok] &= source.n_cells[rows[row_ok]] >= minimum_cells
    col_ok = cols >= 0
    if source.measured is not None:
        col_ok[col_ok] &= source.measured[cols[col_ok]]
    result = np.zeros((len(targets), len(genes)), dtype=np.float32)
    mask = row_ok[:, None] & col_ok[None, :]
    ri, ci = rows[row_ok], cols[col_ok]
    if not len(ri) or not len(ci):
        return result, mask

    ctrl = source.control_probability
    ctrl = ctrl[ci][None, :] if ctrl.ndim == 1 else ctrl[np.ix_(ri, ci)]
    pert = source.probability[np.ix_(ri, ci)]
    if space == "bulk_delta":
        effect = np.log1p(BULK_SCALE * pert) - np.log1p(BULK_SCALE * ctrl)
    else:
        if source.mean_cpm is None:
            pmean, cmean = 1e6 * pert, 1e6 * ctrl
        else:
            pmean = source.mean_cpm[np.ix_(ri, ci)]
            c = source.control_mean_cpm
            cmean = c[ci][None, :] if c.ndim == 1 else c[np.ix_(ri, ci)]
        effect = np.log2((pmean + pseudocount_cpm) / (cmean + pseudocount_cpm))

    if center:
        # exclude each row's own knocked-down gene from the centering mean
        values = effect.copy()
        kept_genes = genes[col_ok]
        own = kept_genes[None, :] == targets[row_ok][:, None]
        values[own] = np.nan
        common = np.nan_to_num(np.nanmean(values, axis=0))
        effect -= common[None, :]
    result[np.ix_(np.flatnonzero(row_ok), np.flatnonzero(col_ok))] = effect
    return result, mask


def cd4_effect(
    stats: dict,
    targets: np.ndarray,
    genes: np.ndarray,
    *,
    minimum_cells: int = MINIMUM_CELLS,
    center: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """CD4 DE log2FC, centred per culture condition on the source-wide mean, then
    averaged over the conditions in which the target passes the publisher quality flags.

    Returns ``(effect float32, available bool)`` on the (targets x genes) panel.
    """
    src_targets = np.asarray(stats["targets"]).astype(str)
    src_genes = np.asarray(stats["genes"]).astype(str)
    values = np.asarray(stats["log2fc"], dtype=float)
    usable = (
        np.asarray(stats["available"])
        & np.asarray(stats["quality_pass"])
        & (np.asarray(stats["n_cells"]) >= minimum_cells)
    )
    targets = np.asarray(targets).astype(str)
    genes = np.asarray(genes).astype(str)
    ri = _index(src_targets, targets)
    ci = _index(src_genes, genes)
    rows, cols = np.flatnonzero(ri >= 0), np.flatnonzero(ci >= 0)
    num = np.zeros((len(targets), len(genes)))
    den = np.zeros_like(num)
    gene_pos = {g: i for i, g in enumerate(src_genes)}
    own = [(k, gene_pos[t]) for k, t in enumerate(src_targets) if t in gene_pos]
    for c in range(len(values)):
        effect = values[c][np.ix_(ri[rows], ci[cols])].copy()
        mask = usable[c, ri[rows], None] & np.isfinite(effect)
        if center:
            centering = np.where(usable[c, :, None], values[c], np.nan)
            for k, g in own:
                centering[k, g] = np.nan
            count = np.isfinite(centering).sum(axis=0)
            full = np.divide(
                np.nansum(centering, axis=0), count, out=np.zeros(len(src_genes)), where=count > 0
            )
            effect -= full[ci[cols]][None, :]
        num[np.ix_(rows, cols)] += np.where(mask, effect, 0)
        den[np.ix_(rows, cols)] += mask
    effect = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    return effect.astype(np.float32), den > 0


def cd4_to_bulk_delta(cd4: np.ndarray, control_probability: np.ndarray) -> np.ndarray:
    """Carry a log2FC into ``bulk_delta`` through the destination control composition."""
    base = np.asarray(control_probability) * BULK_SCALE
    return np.log1p(base[None, :] * np.exp2(np.clip(cd4, -10, 10))) - np.log1p(base[None, :])
