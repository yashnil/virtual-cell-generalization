"""Turning a predicted response into raw single-cell counts.

The challenge scores raw integer counts (``input_type: counts``,
``allow_discrete: false``, ``allow_fractional_counts: false``), so a predictor
that stops at a continuous profile is not yet a submission. Rounding a
log-normalised matrix is not a substitute: it does not reproduce the per-cell
library-size distribution or the sparsity the metrics read, and the
sampling-noise correction in ``expr_mse_unbiased_capped_norm`` is estimated
from exactly that within-group spread.

Three generators are provided, in increasing order of what they are willing to
assume. All three take their per-cell library sizes from the context's own
control cells, so the depth distribution is the measured one rather than an
invention.

``G0`` resampling
    Draw control cells with replacement. The predicted effect is zero. This is
    the honest floor, and the reference point every other generator is judged
    against.

``G1`` transport
    Take a control cell and re-express it under a multiplicative per-gene
    effect, redrawing its counts at its own depth. The cell's identity and
    depth survive; only the composition moves. Redrawing from a single cell's
    own empirical composition would destroy its singletons and cost roughly a
    sixth of its detected genes, so the composition is first blended with the
    pooled control composition -- see ``smoothing``.

``G2`` count model
    Draw from a negative binomial whose mean is the predicted profile and
    whose dispersion is fitted on the controls. This can express an effect no
    control cell exhibits, at the cost of assuming a parametric shape.
"""

from __future__ import annotations

import numpy as np

__all__ = ["MAX_COUNTS_PER_CELL", "resample_controls", "transport_controls", "count_model"]

#: ``max_counts_per_cell`` from ``configs/vcc2026.yaml``.
MAX_COUNTS_PER_CELL = 1_000_000


def _library_sizes(controls: np.ndarray, n_cells: int, rng: np.random.Generator) -> np.ndarray:
    lib = np.asarray(controls.sum(axis=1)).ravel().astype(np.int64)
    lib = lib[lib > 0]
    if lib.size == 0:
        raise ValueError("control cells carry no counts")
    return rng.choice(lib, size=n_cells, replace=True)


def _as_dense(x: np.ndarray) -> np.ndarray:
    return np.asarray(x.toarray() if hasattr(x, "toarray") else x, dtype=np.float64)


def _check(out: np.ndarray) -> np.ndarray:
    totals = out.sum(axis=1)
    if totals.max(initial=0) > MAX_COUNTS_PER_CELL:
        raise ValueError(
            f"generated a cell with {totals.max()} counts, over the {MAX_COUNTS_PER_CELL} cap"
        )
    return out


def resample_controls(
    controls: np.ndarray, n_cells: int, *, rng: np.random.Generator
) -> np.ndarray:
    """``G0``: draw ``n_cells`` control cells with replacement, unchanged."""
    ctrl = _as_dense(controls)
    idx = rng.integers(0, ctrl.shape[0], size=n_cells)
    return _check(ctrl[idx].astype(np.int64))


def transport_controls(
    controls: np.ndarray,
    log2_fold_change: np.ndarray,
    n_cells: int,
    *,
    rng: np.random.Generator,
    smoothing: float = 0.5,
) -> np.ndarray:
    """``G1``: redraw control cells under a multiplicative per-gene effect.

    Each drawn cell keeps its own library size; its composition is multiplied
    by ``2 ** log2_fold_change``, renormalised, and recounted by a multinomial
    draw. A cell's zeros are therefore not forced to stay zero, which is what
    lets an upregulated gene appear.

    ``smoothing`` is the weight given to the pooled control composition when
    forming the per-cell draw probabilities, the rest going to the cell's own.
    It is not cosmetic. A cell's empirical composition assigns probability zero
    to every gene it happens not to have caught, so redrawing straight from it
    is a second round of sampling loss on top of the one the measurement
    already made: at ``smoothing=0`` this generator emits roughly 17% fewer
    detected genes than the control cells it was built from, which would bias
    every differential-expression member of the score. Blending restores the
    detection rate. ``0.5`` is the default because it matches the measured
    control density closely; the value is a modelling choice, not a fitted one.
    """
    if not 0.0 <= smoothing <= 1.0:
        raise ValueError("smoothing must lie in [0, 1]")
    ctrl = _as_dense(controls)
    lfc = np.asarray(log2_fold_change, dtype=np.float64)
    if lfc.shape != (ctrl.shape[1],):
        raise ValueError(f"log2_fold_change must have {ctrl.shape[1]} entries")
    mult = np.exp2(np.clip(lfc, -30.0, 30.0))

    pooled = ctrl.sum(axis=0)
    pooled = pooled / pooled.sum() if pooled.sum() > 0 else pooled

    idx = rng.integers(0, ctrl.shape[0], size=n_cells)
    drawn = ctrl[idx]
    lib = drawn.sum(axis=1).astype(np.int64)

    own = drawn.sum(axis=1, keepdims=True)
    own = np.divide(drawn, own, out=np.zeros_like(drawn), where=own > 0)
    blended = (1.0 - smoothing) * own + smoothing * pooled[None, :]
    weighted = blended * mult[None, :]
    totals = weighted.sum(axis=1, keepdims=True)
    probs = np.divide(weighted, totals, out=np.zeros_like(weighted), where=totals > 0)

    out = np.zeros_like(drawn, dtype=np.int64)
    for i in range(n_cells):
        if lib[i] > 0 and probs[i].sum() > 0:
            out[i] = rng.multinomial(int(lib[i]), probs[i] / probs[i].sum())
    return _check(out)


def count_model(
    controls: np.ndarray,
    target_cpm: np.ndarray,
    n_cells: int,
    *,
    rng: np.random.Generator,
    min_dispersion: float = 0.01,
) -> np.ndarray:
    """``G2``: negative-binomial draws at the predicted composition.

    ``target_cpm`` is the predicted per-gene composition (any positive scale;
    it is renormalised). Per-gene dispersion is fitted on the control cells by
    method of moments, ``phi = (var - mean) / mean**2``, floored at
    ``min_dispersion`` so an under-dispersed gene cannot produce a degenerate
    draw.
    """
    ctrl = _as_dense(controls)
    n_genes = ctrl.shape[1]
    comp = np.asarray(target_cpm, dtype=np.float64)
    if comp.shape != (n_genes,):
        raise ValueError(f"target_cpm must have {n_genes} entries")
    comp = np.clip(comp, 0.0, None)
    if comp.sum() <= 0:
        raise ValueError("target_cpm sums to zero")
    comp = comp / comp.sum()

    ctrl_lib = ctrl.sum(axis=1, keepdims=True)
    frac = np.divide(ctrl, ctrl_lib, out=np.zeros_like(ctrl), where=ctrl_lib > 0)
    mean_f = frac.mean(axis=0)
    var_f = frac.var(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        phi = np.where(mean_f > 0, (var_f - mean_f) / mean_f**2, min_dispersion)
    phi = np.clip(np.nan_to_num(phi, nan=min_dispersion), min_dispersion, None)

    lib = _library_sizes(ctrl, n_cells, rng)
    mu = lib[:, None] * comp[None, :]
    # NB as a gamma-Poisson mixture: shape 1/phi, scale mu*phi.
    rate = rng.gamma(shape=1.0 / phi[None, :], scale=mu * phi[None, :])
    return _check(rng.poisson(rate).astype(np.int64))
