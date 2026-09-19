"""Robustness variants for the independent four-context decomposition.

Everything here is a *sensitivity analysis*. The canonical v1 result and the
module that produced it (``virtual_cell.data.scperteval``) are frozen — see
``data/provenance/scperteval/canonical_v1_freeze.txt`` — so every variant lives
in this separate module and writes to its own output directory.

Variants implemented:

* **control split scheme** — ``"shared"`` reuses one full-context control mean in
  both halves (canonical); ``"split"`` divides the control cells into two
  disjoint halves and references each perturbation half against its own.
* **aggregation order** — ``"mean_log"`` averages ``log1p(CP10K)`` over cells
  (canonical); ``"log_mean"`` recovers CP10K with ``expm1``, averages, then takes
  ``log1p``. This probes Jensen / aggregation-order sensitivity.
* **feature space** — a control-derived global HVG ranking, so a decomposition
  can be recomputed on gene subsets without ever consulting a perturbation
  response.
* **controlled depth** — repeated equal-size subsampling of the *same*
  (context, perturbation) pairs at several cell depths.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from virtual_cell.data.arc2026 import stream_row_chunks
from virtual_cell.data.io import DataIntegrityError
from virtual_cell.data.scperteval import (
    CONTROL_LABEL,
    DEFAULT_CHUNK_SIZE,
    cell_labels,
    read_var_names,
)

CONTROL_SCHEMES = ("shared", "split")
AGGREGATIONS = ("mean_log", "log_mean")


def _column_selector(var_names: pd.Index, genes: Sequence[str]) -> sparse.csr_matrix:
    idx = var_names.get_indexer(pd.Index(genes))
    if (idx < 0).any():
        missing = [g for g, i in zip(genes, idx, strict=True) if i < 0]
        raise DataIntegrityError(f"{len(missing)} requested genes absent, e.g. {missing[:5]}")
    n_sel = len(genes)
    return sparse.csr_matrix(
        (np.ones(n_sel, dtype=np.float64), (idx, np.arange(n_sel))),
        shape=(len(var_names), n_sel),
    )


# --------------------------------------------------------------------------
# aggregation
# --------------------------------------------------------------------------


def aggregate(X: sparse.spmatrix, aggregation: str) -> np.ndarray:
    """Collapse a cell x gene block to one profile under the chosen order.

    ``mean_log`` returns ``mean_cells log1p(CP10K)`` — the values are already on
    that scale, so it is a plain column mean.

    ``log_mean`` returns ``log1p(mean_cells CP10K)``: ``expm1`` is applied to the
    stored values first (``expm1(0) == 0``, so sparsity is preserved), the mean is
    taken, and ``log1p`` applied to the resulting dense profile.
    """
    if aggregation not in AGGREGATIONS:
        raise ValueError(f"Unknown aggregation {aggregation!r}; expected one of {AGGREGATIONS}")
    if X.shape[0] == 0:
        raise DataIntegrityError("Cannot aggregate an empty cell block.")
    if aggregation == "mean_log":
        return np.asarray(X.mean(axis=0)).ravel()
    linear = X.copy().astype(np.float64)
    linear.data = np.expm1(linear.data)
    return np.log1p(np.asarray(linear.mean(axis=0)).ravel())


def group_aggregate(
    X: sparse.csr_matrix, group: np.ndarray, n_groups: int, aggregation: str
) -> np.ndarray:
    """Aggregate rows of ``X`` by ``group`` into an ``(n_groups, G)`` profile matrix."""
    counts = np.bincount(group, minlength=n_groups)
    if (counts == 0).any():
        raise DataIntegrityError(f"{int((counts == 0).sum())} group(s) have no cells.")
    source = X
    if aggregation == "log_mean":
        source = X.copy().astype(np.float64)
        source.data = np.expm1(source.data)
    indicator = sparse.csr_matrix(
        (np.ones(len(group), dtype=np.float64), (group, np.arange(len(group)))),
        shape=(n_groups, X.shape[0]),
    )
    means = np.asarray((indicator @ source).todense()) / counts[:, None]
    return np.log1p(means) if aggregation == "log_mean" else means


# --------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------


def load_cells(
    path: str | Path,
    *,
    genes: Sequence[str],
    perturbations: Sequence[str],
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> tuple[sparse.csr_matrix, np.ndarray, sparse.csr_matrix]:
    """Load one context's perturbation cells **and** its control cells.

    Returns ``(X_pert, group, X_control)`` restricted to ``genes``. Only one
    context should be held at a time.
    """
    path = Path(path)
    selector = _column_selector(read_var_names(path), genes)
    labels = cell_labels(path)
    row_of = {p: i for i, p in enumerate(perturbations)}
    group_all = np.array(
        [n_groups_lookup(lab, row_of) for lab in labels],
        dtype=np.int64,
    )

    pert_blocks, groups, ctrl_blocks = [], [], []
    for start, chunk in stream_row_chunks(path, chunk_size=chunk_size):
        g_chunk = group_all[start : start + chunk.shape[0]]
        pert_rows = np.flatnonzero(g_chunk >= 0)
        ctrl_rows = np.flatnonzero(g_chunk == -2)
        if len(pert_rows):
            pert_blocks.append((chunk[pert_rows] @ selector).astype(np.float32))
            groups.append(g_chunk[pert_rows])
        if len(ctrl_rows):
            ctrl_blocks.append((chunk[ctrl_rows] @ selector).astype(np.float32))
    if not pert_blocks:
        raise DataIntegrityError(f"{path.name}: no cells for the frozen perturbations.")
    if not ctrl_blocks:
        raise DataIntegrityError(f"{path.name}: no control cells.")
    return (
        sparse.vstack(pert_blocks, format="csr"),
        np.concatenate(groups),
        sparse.vstack(ctrl_blocks, format="csr"),
    )


def n_groups_lookup(label: str, row_of: dict[str, int]) -> int:
    """-2 for a control cell, the perturbation row index, or -1 to ignore."""
    if label == CONTROL_LABEL:
        return -2
    return row_of.get(label, -1)


def control_gene_variance(
    path: str | Path, *, genes: Sequence[str], chunk_size: int = DEFAULT_CHUNK_SIZE
) -> np.ndarray:
    """Per-gene variance of ``X`` across this context's **control cells only**.

    Basal quantity: no perturbed cell contributes, so a ranking built from it
    cannot leak a perturbation response.
    """
    path = Path(path)
    selector = _column_selector(read_var_names(path), genes)
    labels = cell_labels(path)
    is_ctrl = labels == CONTROL_LABEL
    n_g = len(genes)
    s1 = np.zeros(n_g)
    s2 = np.zeros(n_g)
    n = 0
    for start, chunk in stream_row_chunks(path, chunk_size=chunk_size):
        rows = np.flatnonzero(is_ctrl[start : start + chunk.shape[0]])
        if not len(rows):
            continue
        sub = (chunk[rows] @ selector).astype(np.float64)
        s1 += np.asarray(sub.sum(axis=0)).ravel()
        sq = sub.copy()
        sq.data = sq.data**2
        s2 += np.asarray(sq.sum(axis=0)).ravel()
        n += sub.shape[0]
    if n == 0:
        raise DataIntegrityError(f"{path.name}: no control cells.")
    mean = s1 / n
    return np.maximum(s2 / n - mean**2, 0.0)


def global_hvg_ranking(variances: dict[str, np.ndarray], genes: Sequence[str]) -> pd.Series:
    """One global gene ranking from per-context control variances.

    **Documented rule**, applied identically to all four contexts: score each gene
    by the *mean of the per-context control-cell variances*, then rank descending.
    A single ranking is used for every context, so the feature space never varies
    by context and never depends on a perturbation response.
    """
    frame = pd.DataFrame(variances, index=pd.Index(genes, name="gene"))
    score = frame.mean(axis=1)
    return score.sort_values(ascending=False, kind="stable")


# --------------------------------------------------------------------------
# split-half with a selectable control scheme
# --------------------------------------------------------------------------


def half_deltas(
    X_pert: sparse.csr_matrix,
    group: np.ndarray,
    X_ctrl: sparse.csr_matrix,
    *,
    n_perturbations: int,
    rng: np.random.Generator,
    control_scheme: str,
    aggregation: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Two half-response matrices ``(P, G)`` for one context and one resample.

    Perturbation cells are always split into two disjoint equal halves (an odd
    cell is dropped). The control reference depends on ``control_scheme``:

    * ``"shared"`` — both halves use the full-context control profile (canonical).
    * ``"split"`` — the control cells are themselves split into two disjoint
      equal halves and each perturbation half is referenced against its own.
    """
    if control_scheme not in CONTROL_SCHEMES:
        raise ValueError(f"Unknown control_scheme {control_scheme!r}")

    order = np.argsort(group, kind="stable")
    Xs = X_pert[order]
    gs = group[order]
    bounds = np.searchsorted(gs, np.arange(n_perturbations + 1))

    idx_a, idx_b, rows_a, rows_b = [], [], [], []
    for p in range(n_perturbations):
        idx = np.arange(bounds[p], bounds[p + 1])
        rng.shuffle(idx)
        h = len(idx) // 2
        if h == 0:
            raise DataIntegrityError(f"Perturbation row {p} has {len(idx)} cell(s); need >= 2.")
        idx_a.append(idx[:h])
        idx_b.append(idx[h : 2 * h])
        rows_a.append(np.full(h, p))
        rows_b.append(np.full(h, p))

    out = []
    if control_scheme == "shared":
        ctrl_a = ctrl_b = aggregate(X_ctrl, aggregation)
    else:
        perm = rng.permutation(X_ctrl.shape[0])
        hc = len(perm) // 2
        if hc == 0:
            raise DataIntegrityError("Too few control cells to split.")
        ctrl_a = aggregate(X_ctrl[np.sort(perm[:hc])], aggregation)
        ctrl_b = aggregate(X_ctrl[np.sort(perm[hc : 2 * hc])], aggregation)

    for idx_list, row_list, ctrl in ((idx_a, rows_a, ctrl_a), (idx_b, rows_b, ctrl_b)):
        sel = np.concatenate(idx_list)
        grp = np.concatenate(row_list)
        means = group_aggregate(Xs[sel], grp, n_perturbations, aggregation)
        out.append(means - ctrl[None, :])
    return out[0], out[1]


# --------------------------------------------------------------------------
# controlled depth experiment
# --------------------------------------------------------------------------


def depth_reliability(
    X_pert: sparse.csr_matrix,
    group: np.ndarray,
    X_ctrl: sparse.csr_matrix,
    *,
    pair_rows: Sequence[int],
    depths: Sequence[int],
    n_repeats: int,
    rng: np.random.Generator,
    aggregation: str = "mean_log",
) -> pd.DataFrame:
    """Split-half reliability of the *same* pairs re-estimated at several depths.

    For each perturbation row and each depth ``n``, ``2n`` cells are drawn without
    replacement and split into two disjoint halves of ``n``; the two half
    responses are correlated. Perturbation identity and context are held fixed, so
    depth is the only thing that varies — this is the controlled version of the
    reliability-versus-depth question, and it never conditions on the observed
    response magnitude.
    """
    order = np.argsort(group, kind="stable")
    Xs = X_pert[order]
    gs = group[order]
    n_p = int(gs.max()) + 1
    bounds = np.searchsorted(gs, np.arange(n_p + 1))
    ctrl = aggregate(X_ctrl, aggregation)

    records = []
    for p in pair_rows:
        available = bounds[p + 1] - bounds[p]
        base = np.arange(bounds[p], bounds[p + 1])
        for n in depths:
            if available < 2 * n:
                continue
            correlations = []
            for _ in range(n_repeats):
                pick = rng.choice(base, size=2 * n, replace=False)
                a = aggregate(Xs[np.sort(pick[:n])], aggregation) - ctrl
                b = aggregate(Xs[np.sort(pick[n:])], aggregation) - ctrl
                sa, sb = a.std(), b.std()
                correlations.append(
                    float(np.corrcoef(a, b)[0, 1]) if sa > 1e-12 and sb > 1e-12 else 0.0
                )
            records.append(
                {
                    "pert_row": int(p),
                    "n_cells": int(n),
                    "available_cells": int(available),
                    "reliability": float(np.mean(correlations)),
                    "reliability_sd": float(np.std(correlations)),
                    "n_repeats": int(n_repeats),
                }
            )
    return pd.DataFrame.from_records(records)
