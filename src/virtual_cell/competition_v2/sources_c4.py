"""C4: KOLF2.1J genome-scale CRISPRi atlas -> C1-format source statistics (ours).

Predeclared in ``reports/competition_v2/c4_predeclaration.md``. The statistics use the
**same definition** as the K562 GWPS source (:func:`sources._matched_statistics`):

* retained targets as passed (the frozen 437);
* NTC controls, matched per KOLF ``batch``;
* duplicated symbols summed.

The KOLF file stores raw counts column-compressed (CSC, 2.66 M cells x 37,567 genes), so
everything is computed in **one streaming pass over gene blocks**, never holding cells in
memory. The same pass accumulates the split-half means that the C3 reliability protocol
needs, with the **exact** half assignment :func:`fusion_c3.split_half` would draw (same
seed, same per-target permutation order), and an NTC pseudo-target null.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy import sparse

from virtual_cell.competition_v2 import sources
from virtual_cell.competition_v2.fusion_c3 import MIN_CPM, MIN_RELIABILITY_CELLS

KOLF_NAME = "KOLF2.1J_iPSC"
KOLF_CONTROL = "NTC"


def _strings(ds) -> np.ndarray:
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in ds[()]])


def split_half_assignment(labels, targets, *, seed: int, min_cells: int = MIN_RELIABILITY_CELLS):
    """Per-cell half (0/1, -1 elsewhere) and the kept targets, drawn as ``split_half`` does."""
    rng = np.random.default_rng(seed)
    half = np.full(len(labels), -1, dtype=np.int8)
    owner = np.full(len(labels), -1, dtype=np.int64)
    kept = []
    for t in targets:
        rows = np.flatnonzero(labels == t)
        if len(rows) < min_cells:
            continue
        rows = rng.permutation(rows)
        h = len(rows) // 2
        half[rows[:h]], half[rows[h:]] = 0, 1
        owner[rows] = len(kept)
        kept.append(str(t))
    return half, owner, kept


def reliability_from_halves(mean_a, mean_b, control_mean_cpm, kept) -> dict:
    """The :func:`fusion_c3.split_half` statistics, from per-half mean-CPM rows."""
    ctrl = np.asarray(control_mean_cpm, dtype=np.float64)
    keep = ctrl >= MIN_CPM
    a = np.log2((np.asarray(mean_a)[:, keep] + 1) / (ctrl[keep] + 1))
    b = np.log2((np.asarray(mean_b)[:, keep] + 1) / (ctrl[keep] + 1))
    a -= a.mean(axis=0, keepdims=True)
    b -= b.mean(axis=0, keepdims=True)
    rel, noise, energy, pearson = {}, {}, {}, {}
    for i, t in enumerate(kept):
        r = np.corrcoef(a[i], b[i])[0, 1]
        pearson[t] = float(r)
        rel[t] = float(max(0.0, 2 * r / (1 + r))) if np.isfinite(r) and r > -1 else 0.0
        noise[t] = float(((a[i] - b[i]) ** 2).sum() / 4)
        energy[t] = float((((a[i] + b[i]) / 2) ** 2).sum())
    return {"genes": keep, "rel": rel, "noise": noise, "energy": energy, "pearson": pearson}


def kolf_statistics(
    path: Path,
    targets: list[str],
    *,
    seed: int,
    half_seed: int,
    null_n: int = 100,
    null_size: int = 218,
    block_genes: int = 512,
    log=print,
) -> tuple[dict, dict]:
    """Stream the KOLF h5ad once -> ``(stats, halves)``.

    ``stats`` has exactly the keys of the other single-cell sources. ``halves`` holds the
    per-half mean-CPM rows for targets and null pseudo-targets, plus the control mean.
    """
    import h5py

    f = h5py.File(path, "r")
    obs = f["obs"]
    labels = _strings(obs["gene_target"]["categories"])[obs["gene_target"]["codes"][()]]
    batch_names = _strings(obs["batch"]["categories"])
    batch_codes = obs["batch"]["codes"][()].astype(np.int64)
    depth = obs["total_counts"][()].astype(np.float64)
    raw_genes = _strings(f["var"][f["var"].attrs.get("_index", "_index")])
    genes, projection = sources.symbol_projection(raw_genes)
    n_cells, n_raw = len(labels), len(raw_genes)
    x = f["X"]
    if x.attrs["encoding-type"] != "csc_matrix" or tuple(x.attrs["shape"]) != (n_cells, n_raw):
        raise ValueError("unexpected KOLF X encoding")
    if (depth <= 0).any():
        raise ValueError("zero-depth cell")

    n_t, n_b = len(targets), len(batch_names)
    lookup = {t: i for i, t in enumerate(targets)}
    control = labels == KOLF_CONTROL
    rng_null = np.random.default_rng(seed + 1)
    null_cells = rng_null.choice(np.flatnonzero(control), null_n * null_size, replace=False)
    null_labels = np.full(n_cells, "", dtype=object)
    if len(null_cells):
        null_labels[null_cells] = [f"__null_{i // null_size}" for i in range(len(null_cells))]
    in_null = np.zeros(n_cells, dtype=bool)
    in_null[null_cells] = True
    group = np.asarray([lookup.get(t, -1) for t in labels], dtype=np.int64)
    main_control = control & ~in_null
    group[main_control] = n_t + batch_codes[main_control]
    group[in_null] = -1

    present = sorted(t for t in targets if (labels == t).sum() >= MIN_RELIABILITY_CELLS)
    half, owner, kept = split_half_assignment(labels, present, seed=half_seed)
    null_names = [f"__null_{i}" for i in range(null_n)]
    nhalf, nowner, nkept = split_half_assignment(null_labels.astype(str), null_names, seed=seed + 2)
    n_h = 2 * len(kept) + 2 * len(nkept)
    hgroup = np.full(n_cells, -1, dtype=np.int64)
    sel = owner >= 0
    hgroup[sel] = 2 * owner[sel] + half[sel]
    sel = nowner >= 0
    hgroup[sel] = 2 * len(kept) + 2 * nowner[sel] + nhalf[sel]

    def assign(g, n, w):
        ok = g >= 0
        return sparse.csr_matrix((w[ok], (g[ok], np.flatnonzero(ok))), shape=(n, n_cells))

    inv = 1e6 / depth
    a_main = assign(group, n_t + n_b, np.ones(n_cells))
    a_main_cpm = assign(group, n_t + n_b, inv)
    a_half_cpm = assign(hgroup, n_h, inv)
    count_sums = np.zeros((n_t + n_b, n_raw))
    cpm_sums = np.zeros_like(count_sums)
    half_sums = np.zeros((n_h, n_raw))
    depth_seen = np.zeros(n_cells)
    indptr = x["indptr"][()]
    for left in range(0, n_raw, block_genes):
        right = min(left + block_genes, n_raw)
        lo, hi = int(indptr[left]), int(indptr[right])
        data = x["data"][lo:hi].astype(np.float64)
        idx = x["indices"][lo:hi]
        if not np.isfinite(data).all() or (data < 0).any() or (data != np.floor(data)).any():
            raise ValueError("KOLF counts must be finite nonnegative integers")
        block = sparse.csc_matrix(
            (data, idx, indptr[left : right + 1] - lo), shape=(n_cells, right - left)
        )
        count_sums[:, left:right] = (a_main @ block).toarray()
        cpm_sums[:, left:right] = (a_main_cpm @ block).toarray()
        half_sums[:, left:right] = (a_half_cpm @ block).toarray()
        depth_seen += np.bincount(idx, weights=data, minlength=n_cells)
        log(f"KOLF genes {right}/{n_raw}")
    f.close()
    if not np.allclose(depth_seen, depth, rtol=1e-6, atol=0.5):
        raise ValueError("obs total_counts differs from the X row sums")

    count_sums = count_sums @ projection
    cpm_sums = cpm_sums @ projection
    half_sums = half_sums @ projection
    group_n = np.bincount(group[group >= 0], minlength=n_t + n_b)
    pert = (group >= 0) & (group < n_t)
    tb_n = np.zeros((n_t, n_b), dtype=np.int64)
    tb_lib = np.zeros((n_t, n_b))
    np.add.at(tb_n, (group[pert], batch_codes[pert]), 1)
    np.add.at(tb_lib, (group[pert], batch_codes[pert]), depth[pert])
    stats = sources._matched_statistics(
        KOLF_NAME, targets, genes, count_sums, cpm_sums, group_n, tb_n, tb_lib
    )
    stats["all_target_cells"] = np.asarray(
        [int((labels == t).sum()) for t in targets], dtype=np.int64
    )
    h_n = np.bincount(hgroup[hgroup >= 0], minlength=n_h).astype(np.float64)
    means = half_sums / np.maximum(h_n, 1)[:, None]
    k = len(kept)
    halves = {
        "kept": np.asarray(kept),
        "mean_a": means[0 : 2 * k : 2],
        "mean_b": means[1 : 2 * k : 2],
        "null_kept": np.asarray(nkept),
        "null_a": means[2 * k :: 2],
        "null_b": means[2 * k + 1 :: 2],
        "control_mean_cpm": stats["global_control_mean_cpm"].astype(np.float64),
        "genes": genes,
    }
    return stats, halves
