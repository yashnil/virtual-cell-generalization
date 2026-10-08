"""Build the N5 K562 GWPS data object (reports/n5_protocol.md §1).

scPertEval recipe on GWPS raw counts: drop cells with < 200 non-zero genes on the
8,248-gene axis; per-cell CP10K over that axis; log1p; restrict to the N5 genes.
Canonical deltas, disjoint F/E1/E2 parts (5 repeats) and a depth-matched variant
(n = RPE1 canonical cell count per perturbation). Gates G2–G4 stop the build.

Reproduce: ``uv run python scripts/research_v3/build_n5_gwps.py``
"""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

REPO = Path(__file__).resolve().parents[2]
GWPS = REPO / "data" / "raw" / "competition_v2" / "K562_gwps_raw_singlecell_01.h5ad"
AXES = REPO / "data" / "splits" / "n5_k562"
N3 = REPO / "outputs" / "n3" / "data"
N3AX = REPO / "data" / "splits" / "six_context_n3"
OUT = REPO / "outputs" / "n5" / "gwps"
SEED = 20261006
N_REPEATS = 5
BLOCK = 4096
MIN_GENES = 200


def three_way(idx):
    n = len(idx)
    h, q = n // 2, n // 4
    return idx[:h], idx[h : h + q], idx[h + q : h + 2 * q]


def read_obs(f):
    o, cats = f["obs"], f["obs"]["__categories"]
    gene = pd.Categorical.from_codes(o["gene"][:], cats["gene"][:].astype(str)).astype(str)
    v = f["var"]
    names = pd.Categorical.from_codes(
        v["gene_name"][:], v["__categories"]["gene_name"][:].astype(str)
    )
    return np.asarray(gene), np.asarray(names.astype(str))


def assign_parts(groups: np.ndarray, n_groups: int, seed_parts, rows: np.ndarray) -> np.ndarray:
    """rowid[r, i] = part * n_groups + group (or -1). ``rows`` orders cells within a group."""
    out = np.full((N_REPEATS, len(groups)), -1, dtype=np.int64)
    order = np.lexsort((rows, groups))
    bounds = np.searchsorted(groups[order], np.arange(n_groups + 1))
    for r in range(N_REPEATS):
        rng = np.random.default_rng([*seed_parts, r])
        for g in range(n_groups):
            idx = order[bounds[g] : bounds[g + 1]].copy()
            rng.shuffle(idx)
            for j, part in enumerate(three_way(idx)):
                out[r, part] = j * n_groups + g
    return out


def main() -> None:
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    perts = (AXES / "shared_perturbations.txt").read_text().split()
    genes = (AXES / "shared_genes.txt").read_text().split()
    P, G = len(perts), len(genes)
    f = h5py.File(GWPS, "r")
    gene, var_names = read_obs(f)
    n_cells, n_var = f["X"].shape
    X = np.memmap(
        GWPS, dtype=np.float32, mode="r", offset=f["X"].id.get_offset(), shape=(n_cells, n_var)
    )
    col = {g: i for i, g in enumerate(var_names)}
    gcols = np.array([col[g] for g in genes])
    row_of = {p: i for i, p in enumerate(perts)}
    grp_all = np.array([row_of.get(g, P if g == "non-targeting" else -1) for g in gene])
    sel = np.flatnonzero(grp_all >= 0)

    # pass 1: cell filter (non-zero genes on the full axis)
    nnz = np.empty(len(sel), dtype=np.int64)
    for s in range(0, len(sel), BLOCK):
        nnz[s : s + BLOCK] = np.count_nonzero(X[sel[s : s + BLOCK]], axis=1)
    keep = nnz >= MIN_GENES
    rows = sel[keep]
    grp = grp_all[rows]
    print(
        f"pass 1: {len(sel)} selected, {int(keep.sum())} kept [{time.time() - t0:.0f}s]", flush=True
    )

    # split assignments: main parts (perturbations and controls together, group P = control)
    rowid = assign_parts(grp, P + 1, (SEED, 6), rows)
    # depth-matched subset: n = RPE1 canonical count per perturbation
    p3 = (N3AX / "shared_perturbations.txt").read_text().split()
    n_rpe1 = np.load(N3 / "cell_counts6.npy")[1][[p3.index(p) for p in perts]].astype(int)
    in_depth = np.zeros(len(rows), dtype=bool)
    order = np.lexsort((rows, grp))
    bounds = np.searchsorted(grp[order], np.arange(P + 2))
    capped = 0
    for p in range(P):
        idx = order[bounds[p] : bounds[p + 1]]
        n = min(len(idx), n_rpe1[p])
        capped += int(len(idx) < n_rpe1[p])
        pick = np.random.default_rng([SEED, 7, p]).choice(idx, n, replace=False)
        in_depth[pick] = True
    in_depth[grp == P] = True  # controls: full set, split as in the main object
    sub = np.flatnonzero(in_depth)
    rowid_depth = np.full((N_REPEATS, len(rows)), -1, dtype=np.int64)
    pert_sub = sub[grp[sub] < P]
    rowid_depth[:, pert_sub] = assign_parts(grp[pert_sub], P, (SEED, 8), rows[pert_sub])
    ctrl_rows = np.flatnonzero(grp == P)
    for r in range(N_REPEATS):  # controls keep their main-object part
        part = rowid[r, ctrl_rows] // (P + 1)
        rowid_depth[r, ctrl_rows] = np.where(rowid[r, ctrl_rows] >= 0, part * (P + 1) + P, -1)
        # shift perturbation codes from base P to base P + 1
        m = rowid_depth[r, pert_sub] >= 0
        code = rowid_depth[r, pert_sub][m]
        rowid_depth[r, pert_sub[m]] = (code // P) * (P + 1) + code % P

    # pass 2: accumulate log1p CP10K sums
    full = np.zeros((P + 1, G))
    full_depth = np.zeros((P + 1, G))
    parts = np.zeros((N_REPEATS, 3 * (P + 1), G))
    parts_depth = np.zeros((N_REPEATS, 3 * (P + 1), G))
    for s in range(0, len(rows), BLOCK):
        r_idx = np.arange(s, min(s + BLOCK, len(rows)))
        block = np.asarray(X[rows[r_idx]], dtype=np.float64)
        lib = block.sum(axis=1, keepdims=True)
        vals = np.log1p(1e4 * block[:, gcols] / lib)

        def add(target, codes, n_rows, vals=vals):
            m = codes >= 0
            ind = sparse.csr_matrix(
                (np.ones(m.sum()), (codes[m], np.flatnonzero(m))), shape=(n_rows, len(codes))
            )
            target += ind @ vals

        add(full, grp[r_idx], P + 1)
        add(full_depth, np.where(in_depth[r_idx], grp[r_idx], -1), P + 1)
        for r in range(N_REPEATS):
            add(parts[r], rowid[r, r_idx], 3 * (P + 1))
            add(parts_depth[r], rowid_depth[r, r_idx], 3 * (P + 1))
    print(f"pass 2 done [{time.time() - t0:.0f}s]", flush=True)

    def means(sums, codes, n_rows):
        n = np.bincount(codes[codes >= 0], minlength=n_rows).astype(float)
        return sums / np.maximum(n, 1)[:, None], n

    full_m, n_full = means(full, grp, P + 1)
    depth_m, n_depth = means(full_depth, np.where(in_depth, grp, -1), P + 1)
    pm = np.zeros((N_REPEATS, 3, P + 1, G))
    pmd = np.zeros_like(pm)
    n_parts = np.zeros((N_REPEATS, 3, P + 1))
    n_parts_d = np.zeros_like(n_parts)
    for r in range(N_REPEATS):
        m, n = means(parts[r], rowid[r], 3 * (P + 1))
        pm[r], n_parts[r] = m.reshape(3, P + 1, G), n.reshape(3, P + 1)
        m, n = means(parts_depth[r], rowid_depth[r], 3 * (P + 1))
        pmd[r], n_parts_d[r] = m.reshape(3, P + 1, G), n.reshape(3, P + 1)

    # gates
    exact = n_parts[0, :, :P].sum(axis=0) == n_full[:P]
    recon = (parts[0].reshape(3, P + 1, G)[:, :P].sum(axis=0))[exact] / n_full[:P][exact, None]
    g2 = float(np.max(np.abs(recon - full_m[:P][exact])) / max(1e-12, np.max(np.abs(full_m[:P]))))
    if g2 > 1e-6:
        raise RuntimeError(f"G2 failed: {g2}")
    delta = full_m[:P] - full_m[P]
    gi = {g: i for i, g in enumerate(genes)}
    on = [delta[i, gi[p]] for i, p in enumerate(perts) if p in gi]
    g3 = float(np.median(on))
    if not g3 < 0:
        raise RuntimeError(f"G3 failed: {g3}")
    sizes = np.bincount(grp, minlength=P + 1)
    for r in range(N_REPEATS):
        if not np.array_equal(n_parts[r].sum(axis=0), sizes // 2 + 2 * (sizes // 4)):
            raise AssertionError("G4 failed")

    np.save(OUT / "delta.npy", delta.astype(np.float32))
    np.save(OUT / "delta_depth.npy", (depth_m[:P] - depth_m[P]).astype(np.float32))
    np.save(OUT / "pert_part_means.npy", pm[:, :, :P].astype(np.float32))
    np.save(OUT / "ctrl_part_means.npy", pm[:, :, P].astype(np.float32))
    np.save(OUT / "pert_part_means_depth.npy", pmd[:, :, :P].astype(np.float32))
    np.save(OUT / "ctrl_part_means_depth.npy", pmd[:, :, P].astype(np.float32))
    np.save(OUT / "cell_counts.npy", n_full[:P])
    np.save(OUT / "cell_counts_depth.npy", n_depth[:P])
    gates = {
        "cells_selected": int(len(sel)),
        "cells_kept": int(len(rows)),
        "controls_kept": int(n_full[P]),
        "G2_max_rel_err": g2,
        "G2_n_checked": int(exact.sum()),
        "G3_median_on_target": g3,
        "G3_n": len(on),
        "G4": "asserted",
        "depth_matched_capped_perturbations": capped,
        "median_cells_full": float(np.median(n_full[:P])),
        "median_cells_depth": float(np.median(n_depth[:P])),
        "min_part_cells": n_parts[:, :, :P].min(axis=(0, 2)).tolist(),
    }
    (OUT / "build_gates.json").write_text(json.dumps(gates, indent=2))
    with open(OUT / "sha256.txt", "w") as fh:
        for p in sorted(OUT.glob("*.npy")):
            h = hashlib.sha256()
            with open(p, "rb") as src:
                for blk in iter(lambda: src.read(1 << 24), b""):
                    h.update(blk)
            fh.write(f"{h.hexdigest()}  {p.name}\n")
    print(json.dumps(gates, indent=2))


if __name__ == "__main__":
    main()
