"""Per-source pseudobulk statistics for the atlas backbone (our implementation).

Attribution: the statistics, their names and the pseudocount shrinkage are those of
AtlasShift (https://github.com/kaipengm2/Virtual-Cell-Challenge-2026 @ ``d24ce4f``,
MIT, ``prepare.py`` / ``model.load_xatlas``). This module re-derives them from the raw
public files in our own code so that the C1 pipeline never reads a vendored file or a
fused vector, and so a source can be dropped (licensing) or held out (public folds)
without touching ``third_party/``. ``tests/test_competition_v2_c1.py`` pins numerical
equality with the vendored functions.

A source reduces to five arrays over (retained target, gene):

* ``target_count_sums`` — raw counts summed over the target's cells;
* ``target_mean_cpm`` — the mean of per-cell CPM (cell-weighted);
* ``matched_control_probability`` — control composition matched to the target's batches,
  weighted by the target's library mass per batch;
* ``matched_control_mean_cpm`` — the matched control's mean per-cell CPM, weighted by the
  target's cell count per batch;
* ``n_cells``.

Only GREEN sources are prepared here (K562 GWPS, VCC 2025 H1, CD4 DE); X-Atlas is
deliberately not implemented in this namespace. See
``reports/competition_v2/data_license_register.md``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

#: Label of control cells in the Perturb-seq sources we read.
CONTROL_LABEL = "non-targeting"
#: AtlasShift ``load_xatlas(prior_counts=...)`` default used by ``predict.py``.
PRIOR_COUNTS = 100_000.0
#: AtlasShift ``minimum_cells``.
MINIMUM_CELLS = 20


@dataclass
class SourceStats:
    """One source's shrunk per-target statistics on its own gene axis."""

    name: str
    targets: np.ndarray
    genes: np.ndarray
    probability: np.ndarray  # shrunk pseudobulk composition, float32
    control_probability: np.ndarray  # matched control composition, float32
    n_cells: np.ndarray
    mean_cpm: np.ndarray | None = None  # shrunk mean per-cell CPM, float32
    control_mean_cpm: np.ndarray | None = None
    measured: np.ndarray | None = None

    def usable(self, minimum_cells: int = MINIMUM_CELLS) -> np.ndarray:
        """Retained targets with enough cells to be used as direct evidence."""
        return self.targets[self.n_cells >= minimum_cells]


def shrink_statistics(raw: dict, *, prior_counts: float = PRIOR_COUNTS, name: str | None = None):
    """Apply the sequencing-mass pseudocount shrinkage to raw per-target sums.

    ``prior_counts`` pseudo-counts of the matched control composition are added to each
    target's count sums; mean CPM is blended toward the control mean by
    ``N / (N + prior_counts)`` with ``N`` the target's total counts. dtypes mirror the
    upstream loader (float64 arithmetic, float32 storage) so equality is exact.
    """
    counts = np.asarray(raw["target_count_sums"], dtype=np.float64)
    ctrl = np.asarray(raw["matched_control_probability"], dtype=np.float64)
    prior = counts + prior_counts * ctrl
    row_mass = prior.sum(axis=1, keepdims=True)
    probability = np.divide(prior, row_mass, out=np.zeros_like(prior), where=row_mass > 0)
    total = counts.sum(axis=1)
    fraction = np.divide(
        total, total + prior_counts, out=np.zeros(len(counts)), where=total + prior_counts > 0
    )
    control_cpm = np.asarray(raw["matched_control_mean_cpm"], dtype=np.float64)
    mean_cpm = fraction[:, None] * raw["target_mean_cpm"] + (1 - fraction[:, None]) * control_cpm
    return SourceStats(
        name=str(raw["source"]) if name is None else name,
        targets=np.asarray(raw["targets"]).astype(str),
        genes=np.asarray(raw["genes"]).astype(str),
        probability=probability.astype(np.float32),
        control_probability=ctrl.astype(np.float32),
        n_cells=np.asarray(raw["n_cells"]),
        mean_cpm=mean_cpm.astype(np.float32),
        control_mean_cpm=control_cpm.astype(np.float32),
        measured=np.asarray(raw["measured_genes"]),
    )


def load_source(path: Path, *, prior_counts: float = PRIOR_COUNTS) -> SourceStats:
    """Load a prepared ``*_statistics.npz`` and shrink it."""
    with np.load(path, allow_pickle=False) as d:
        return shrink_statistics({k: d[k] for k in d.files}, prior_counts=prior_counts)


def retained_targets(official_targets, h1_train_targets) -> list[str]:
    """The frozen retained-target set: 300 official targets, then H1 2025 training targets.

    Every single-cell source keeps exactly these rows, and each source's centering mean
    runs over them, so this ordering is part of the frozen backbone definition.
    """
    extra = [t for t in h1_train_targets if t != CONTROL_LABEL]
    return list(dict.fromkeys(list(official_targets) + extra))


def _check_counts(x: sparse.csr_matrix) -> np.ndarray:
    if not np.isfinite(x.data).all() or (x.data < 0).any() or (x.data != np.floor(x.data)).any():
        raise ValueError("source counts must be finite nonnegative integers")
    depth = np.asarray(x.sum(axis=1)).ravel()
    if (depth <= 0).any():
        raise ValueError("zero-depth cell")
    return depth


def symbol_projection(symbols: np.ndarray) -> tuple[np.ndarray, sparse.csr_matrix]:
    """Collapse duplicated gene symbols by summation (first-occurrence order)."""
    genes = np.asarray(list(dict.fromkeys(symbols)), dtype=str)
    lookup = {g: i for i, g in enumerate(genes)}
    projection = sparse.csr_matrix(
        (np.ones(len(symbols)), (np.arange(len(symbols)), [lookup[g] for g in symbols])),
        shape=(len(symbols), len(genes)),
    )
    return genes, projection


def prepare_k562(
    path: Path,
    targets: list[str],
    *,
    eval_controls: int = 8000,
    eval_cap: int = 1000,
    seed: int = 20260910,
    chunk: int = 1024,
    log=print,
):
    """K562 GWPS raw single cells -> statistics, plus held-out evaluation cells.

    Batches are ``gem_group``. The statistics match upstream ``prepare_k562``. The same
    pass also keeps (for the public K562 fold only) every retained target's cells, capped
    at ``eval_cap`` by a seeded draw, and ``eval_controls`` random non-targeting cells.
    Evaluation cells are on the collapsed-symbol axis.
    """
    import anndata as ad

    a = ad.read_h5ad(path, backed="r")
    labels = a.obs["gene"].astype(str).to_numpy()
    batch_codes, batches = pd.factorize(a.obs["gem_group"].astype(str))
    genes, projection = symbol_projection(a.var["gene_name"].astype(str).to_numpy())
    n_t, n_b = len(targets), len(batches)
    lookup = {t: i for i, t in enumerate(targets)}
    group = np.asarray([lookup.get(t, -1) for t in labels], dtype=np.int64)
    control = labels == CONTROL_LABEL
    group[control] = n_t + batch_codes[control]
    rows = np.flatnonzero(group >= 0)

    rng = np.random.default_rng(seed)
    keep_eval = np.zeros(len(labels), dtype=bool)
    for t in range(n_t):
        members = np.flatnonzero(group == t)
        if len(members) > eval_cap:
            members = rng.choice(members, eval_cap, replace=False)
        keep_eval[members] = True
    keep_eval[rng.choice(np.flatnonzero(control), eval_controls, replace=False)] = True

    shape = (n_t + n_b, len(genes))
    count_sums = np.zeros(shape)
    cpm_sums = np.zeros(shape)
    group_n = np.bincount(group[rows], minlength=shape[0])
    tb_n = np.zeros((n_t, n_b), dtype=np.int64)
    tb_lib = np.zeros((n_t, n_b))
    library = 0.0
    eval_rows, eval_blocks = [], []
    for k, left in enumerate(range(0, len(rows), chunk)):
        sel = rows[left : left + chunk]
        x = sparse.csr_matrix(a.X[sel]).astype(np.float64)
        depth = _check_counts(x)
        library += depth.sum()
        x = x @ projection
        g = group[sel]
        assign = sparse.csr_matrix(
            (np.ones(len(sel)), (g, np.arange(len(sel)))), shape=(shape[0], len(sel))
        )
        count_sums += (assign @ x).toarray()
        cpm_sums += (assign @ x.multiply((1e6 / depth)[:, None])).toarray()
        pert = g < n_t
        np.add.at(tb_n, (g[pert], batch_codes[sel][pert]), 1)
        np.add.at(tb_lib, (g[pert], batch_codes[sel][pert]), depth[pert])
        ev = keep_eval[sel]
        if ev.any():
            eval_rows.append(sel[ev])
            eval_blocks.append(sparse.csr_matrix(x[ev]))
        if k % 100 == 0:
            log(f"K562: {left + len(sel)}/{len(rows)} rows")
    a.file.close()
    np.testing.assert_allclose(count_sums.sum(), library, rtol=1e-12)
    stats = _matched_statistics(
        "K562_GWPS_CPM", targets, genes, count_sums, cpm_sums, group_n, tb_n, tb_lib
    )
    order = np.concatenate(eval_rows)
    cells = {
        "counts": sparse.vstack(eval_blocks, format="csr"),
        "labels": labels[order],
        "genes": genes,
    }
    return stats, cells


def _matched_statistics(name, targets, genes, count_sums, cpm_sums, group_n, tb_n, tb_lib):
    n_t = len(targets)
    control_n = group_n[n_t:]
    if (control_n == 0).any():
        raise ValueError("a batch lacks control cells")
    controls = count_sums[n_t:]
    control_probability = controls / controls.sum(axis=1, keepdims=True)
    control_cpm = cpm_sums[n_t:] / control_n[:, None]
    n = group_n[:n_t]
    counts = count_sums[:n_t]
    means = np.divide(cpm_sums[:n_t], n[:, None], out=np.zeros_like(counts), where=n[:, None] > 0)
    cell_w = np.divide(tb_n, n[:, None], out=np.zeros_like(tb_lib), where=n[:, None] > 0)
    lib = tb_lib.sum(axis=1, keepdims=True)
    lib_w = np.divide(tb_lib, lib, out=np.zeros_like(tb_lib), where=lib > 0)
    return {
        "source": np.asarray(name),
        "targets": np.asarray(targets),
        "genes": genes,
        "n_cells": n,
        "target_count_sums": counts.astype(np.float32),
        "target_mean_cpm": means.astype(np.float32),
        "matched_control_probability": (lib_w @ control_probability).astype(np.float32),
        "matched_control_mean_cpm": (cell_w @ control_cpm).astype(np.float32),
        "global_control_probability": (controls.sum(axis=0) / controls.sum()).astype(np.float32),
        "global_control_mean_cpm": (cpm_sums[n_t:].sum(axis=0) / control_n.sum()).astype(
            np.float32
        ),
        "measured_genes": np.ones(len(genes), dtype=bool),
    }


def summarise_h1_split(path: Path, *, chunk: int = 4096, log=print) -> dict:
    """Per-label count sums and mean CPM over every cell of one VCC 2025 H1 split."""
    import anndata as ad

    a = ad.read_h5ad(path, backed="r")
    try:
        labels = a.obs["target_gene"].astype(str).to_numpy()
        groups = np.asarray(sorted(set(labels)), dtype=str)
        ids = pd.Categorical(labels, categories=groups).codes
        n = np.bincount(ids, minlength=len(groups))
        genes = a.var_names.astype(str).to_numpy(dtype=str)
        counts = np.zeros((len(groups), len(genes)))
        cpm = np.zeros_like(counts)
        for left in range(0, len(labels), chunk):
            right = min(left + chunk, len(labels))
            x = sparse.csr_matrix(a.X[left:right], dtype=np.float64)
            depth = _check_counts(x)
            member = sparse.csr_matrix(
                (np.ones(right - left), (ids[left:right], np.arange(right - left))),
                shape=(len(groups), right - left),
            )
            counts += (member @ x).toarray()
            cpm += (member @ x.multiply(1e6 / depth[:, None])).toarray()
            if left % (chunk * 16) == 0:
                log(f"{Path(path).name}: {right}/{len(labels)} cells")
        return {
            "targets": groups,
            "genes": genes,
            "n_cells": n,
            "count_sums": counts,
            "mean_cpm": cpm / n[:, None],
        }
    finally:
        a.file.close()


def assemble_h1(summaries: dict[str, dict]) -> dict:
    """Join the three 2025 splits; each target is matched to its own split's controls."""
    out = {k: [] for k in ["targets", "counts", "means", "n", "cp", "cm", "split"]}
    genes = summaries["train"]["genes"]
    global_cp = global_cm = None
    for split in ["train", "validation", "test"]:
        d = summaries[split]
        if not np.array_equal(genes, d["genes"]):
            raise ValueError("H1 gene axes differ between splits")
        c = np.flatnonzero(d["targets"] == CONTROL_LABEL)
        if len(c) != 1:
            raise ValueError("expected one control group per split")
        c = int(c[0])
        ctrl = d["count_sums"][c] / d["count_sums"][c].sum()
        cmean = d["mean_cpm"][c]
        if split == "train":
            global_cp, global_cm = ctrl.copy(), cmean.copy()
        for i, t in enumerate(d["targets"]):
            if t == CONTROL_LABEL:
                continue
            if str(t) in out["targets"]:
                raise ValueError("a target appears in two H1 splits")
            out["targets"].append(str(t))
            out["counts"].append(d["count_sums"][i])
            out["means"].append(d["mean_cpm"][i])
            out["n"].append(int(d["n_cells"][i]))
            out["cp"].append(ctrl)
            out["cm"].append(cmean)
            out["split"].append(split)
    return {
        "source": np.asarray("H1_2025_public"),
        "targets": np.asarray(out["targets"]),
        "genes": genes,
        "n_cells": np.asarray(out["n"]),
        "target_count_sums": np.asarray(out["counts"], dtype=np.float32),
        "target_mean_cpm": np.asarray(out["means"], dtype=np.float32),
        "matched_control_probability": np.asarray(out["cp"], dtype=np.float32),
        "matched_control_mean_cpm": np.asarray(out["cm"], dtype=np.float32),
        "global_control_probability": global_cp.astype(np.float32),
        "global_control_mean_cpm": global_cm.astype(np.float32),
        "measured_genes": np.ones(len(genes), dtype=bool),
        "public_2025_split": np.asarray(out["split"]),
    }


#: CD4 publisher quality flags; a (condition, target) is usable only if every row passes.
CD4_QUALITY = "n_guides>=2 & !single_guide & ontarget_significant & !distal_offtarget & !low_gex"


def prepare_cd4(path: Path, targets: list[str]) -> dict:
    """CD4 genome-wide DE table -> per (condition, target) mean log2FC and quality flags."""
    import h5py
    from anndata.io import read_elem

    with h5py.File(path, "r") as f:
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])
        genes = (
            var["gene_name"].astype(str).to_numpy(dtype=str)
            if "gene_name" in var
            else var.index.astype(str).to_numpy(dtype=str)
        )
        if len(set(genes)) != len(genes):
            raise ValueError("duplicate CD4 gene symbols")
        selected = np.flatnonzero(obs["target_contrast_gene_name"].astype(str).isin(targets))
        chosen = obs.iloc[selected]
        layers = {
            name: np.asarray(f["layers"][name][selected, :], dtype=np.float32)
            for name in ["log_fc", "adj_p_value", "lfcSE"]
        }
    conditions = sorted(chosen["culture_condition"].astype(str).unique())
    shape = (len(conditions), len(targets), len(genes))
    effects = np.zeros(shape, np.float32)
    q = np.ones(shape, np.float32)
    se = np.zeros(shape, np.float32)
    available = np.zeros(shape[:2], bool)
    quality = np.zeros(shape[:2], bool)
    n_cells = np.zeros(shape[:2], np.int64)
    n_rows = np.zeros(shape[:2], np.int32)
    cond = chosen["culture_condition"].astype(str).to_numpy()
    tgt = chosen["target_contrast_gene_name"].astype(str).to_numpy()
    for c, condition in enumerate(conditions):
        for t, target in enumerate(targets):
            rows = np.flatnonzero((cond == condition) & (tgt == target))
            if not len(rows):
                continue
            meta = chosen.iloc[rows]
            values = layers["log_fc"][rows]
            if not np.isfinite(values).all():
                raise ValueError("nonfinite publisher log2FC")
            effects[c, t] = values.mean(axis=0)
            q[c, t] = layers["adj_p_value"][rows].max(axis=0)
            se[c, t] = layers["lfcSE"][rows].mean(axis=0)
            available[c, t] = True
            n_rows[c, t] = len(rows)
            n_cells[c, t] = int(meta["n_cells_target"].min())
            quality[c, t] = bool(
                (
                    (meta["n_guides"] >= 2)
                    & ~meta["single_guide_estimate"].astype(bool)
                    & meta["ontarget_significant"].astype(bool)
                    & ~meta["distal_offtarget_flag"].astype(bool)
                    & ~meta["low_target_gex"].astype(bool)
                ).all()
            )
    return {
        "targets": np.asarray(targets),
        "genes": genes,
        "conditions": np.asarray(conditions),
        "log2fc": effects,
        "adjusted_p": q,
        "lfcSE": se,
        "available": available,
        "quality_pass": quality,
        "n_cells": n_cells,
        "rows_per_result": n_rows,
    }
