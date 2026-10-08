"""Build the N6 data object (reports/n6_protocol.md §1).

1. VIPerturb-seq K562 (converted h5ad): scPertEval recipe (cells with ≥ 200 genes,
   CP10K over the full probe axis, log1p), restricted to the N6 axes; canonical
   deltas against all NO-TARGET cells; disjoint F/E1/E2 parts (5 repeats).
2. GWPS subsampled per perturbation to the VIPerturb cell count (positive control
   at VIPerturb depth), same recipe as N5.
3. Target, GWPS (full), RPE1 and Jurkat subset exactly from the frozen N3/N5 objects.

Reproduce: ``uv run python scripts/research_v3/build_n6.py``
"""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import h5py
import numpy as np
from build_n5_gwps import BLOCK, GWPS, MIN_GENES, N_REPEATS, assign_parts, read_obs
from scipy import sparse

REPO = Path(__file__).resolve().parents[2]
VIP = REPO / "data" / "raw" / "viperturb" / "viperturb_k562_genome_wide.h5ad"
AX = REPO / "data" / "splits" / "n6_k562"
OUT = REPO / "outputs" / "n6" / "data"
SEED = 20261006
CONTROL_VIP = "NO-TARGET"


def means_from_sums(full, parts, grp, rowid, P, G):
    n_full = np.bincount(grp, minlength=P + 1).astype(float)
    fm = full / np.maximum(n_full, 1)[:, None]
    pm = np.zeros((N_REPEATS, 3, P + 1, G))
    n_parts = np.zeros((N_REPEATS, 3, P + 1))
    for r in range(N_REPEATS):
        n = np.bincount(rowid[r][rowid[r] >= 0], minlength=3 * (P + 1)).astype(float)
        pm[r] = (parts[r] / np.maximum(n, 1)[:, None]).reshape(3, P + 1, G)
        n_parts[r] = n.reshape(3, P + 1)
    return fm, pm, n_full, n_parts


def accumulate(blocks, grp, rowid, P, G):
    """``blocks`` yields (row_positions, values (n, G) log1p CP10K)."""
    full = np.zeros((P + 1, G))
    parts = np.zeros((N_REPEATS, 3 * (P + 1), G))
    for pos, vals in blocks:

        def add(target, codes, n_rows, vals=vals):
            m = codes >= 0
            ind = sparse.csr_matrix(
                (np.ones(m.sum()), (codes[m], np.flatnonzero(m))), shape=(n_rows, len(codes))
            )
            target += ind @ vals

        add(full, grp[pos], P + 1)
        for r in range(N_REPEATS):
            add(parts[r], rowid[r, pos], 3 * (P + 1))
    return full, parts


def gates(fm, pm, n_full, n_parts, P, genes, perts, label):
    delta = fm[:P] - fm[P]
    gi = {g: i for i, g in enumerate(genes)}
    on = [delta[i, gi[p]] for i, p in enumerate(perts) if p in gi]
    med = float(np.median(on)) if on else float("nan")
    if not med < 0:
        raise RuntimeError(f"{label}: G3 failed (median on-target {med})")
    sizes = n_full.astype(int)
    for r in range(N_REPEATS):
        if not np.array_equal(n_parts[r].sum(axis=0).astype(int), sizes // 2 + 2 * (sizes // 4)):
            raise AssertionError(f"{label}: G4 failed")
    return {"G3_median_on_target": med, "G3_n": len(on)}


def build_viperturb(perts, vip_labels, genes, vip_genes):
    P, G = len(perts), len(genes)
    f = h5py.File(VIP, "r")
    obs = f["obs"]
    lab = obs["gene"]
    if isinstance(lab, h5py.Group):  # categorical
        cats = lab["categories"][:].astype(str)
        gene = cats[lab["codes"][:]]
    else:
        gene = lab[:].astype(str)
    var = f["var"]
    var_names = var[var.attrs["_index"]][:].astype(str)
    gcol = {g: i for i, g in enumerate(var_names)}
    gcols = np.array([gcol[g] for g in vip_genes])  # alias-mapped VIPerturb symbols
    row_of = {v: i for i, v in enumerate(vip_labels)}
    grp_all = np.array([row_of.get(g, P if g == CONTROL_VIP else -1) for g in gene])
    sel = np.flatnonzero(grp_all >= 0)
    X = f["X"]
    indptr = X["indptr"][:]
    data, indices = X["data"], X["indices"]

    def row_block(rows):
        mats = []
        for r0 in rows:
            a, b = indptr[r0], indptr[r0 + 1]
            mats.append((indices[a:b], data[a:b]))
        ip = np.concatenate([[0], np.cumsum([len(m[0]) for m in mats])])
        return sparse.csr_matrix(
            (
                np.concatenate([m[1] for m in mats]).astype(np.float64),
                np.concatenate([m[0] for m in mats]),
                ip,
            ),
            shape=(len(rows), len(var_names)),
        )

    nnz = np.diff(indptr)[sel]
    keep = nnz >= MIN_GENES
    rows = sel[keep]
    grp = grp_all[rows]
    rowid = assign_parts(grp, P + 1, (SEED, 30), rows)

    def blocks():
        for s in range(0, len(rows), BLOCK):
            pos = np.arange(s, min(s + BLOCK, len(rows)))
            m = row_block(rows[pos])
            lib = np.asarray(m.sum(axis=1)).ravel()
            sub = m[:, gcols].toarray()
            yield pos, np.log1p(1e4 * sub / lib[:, None])

    full, parts = accumulate(blocks(), grp, rowid, P, G)
    fm, pm, n_full, n_parts = means_from_sums(full, parts, grp, rowid, P, G)
    g = gates(fm, pm, n_full, n_parts, P, genes, perts, "VIPerturb")
    g.update(
        {
            "cells_selected": int(len(sel)),
            "cells_kept": int(len(rows)),
            "controls_kept": int(n_full[P]),
        }
    )
    return fm, pm, n_full, n_parts, g


def build_gwps_at_depth(perts, genes, n_target):
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
    nnz = np.empty(len(sel), dtype=np.int64)
    for s in range(0, len(sel), BLOCK):
        nnz[s : s + BLOCK] = np.count_nonzero(X[sel[s : s + BLOCK]], axis=1)
    rows_all = sel[nnz >= MIN_GENES]
    grp_k = grp_all[rows_all]
    keep = grp_k == P  # all controls
    order = np.lexsort((rows_all, grp_k))
    bounds = np.searchsorted(grp_k[order], np.arange(P + 2))
    for p in range(P):
        idx = order[bounds[p] : bounds[p + 1]]
        n = int(min(len(idx), n_target[p]))
        keep[np.random.default_rng([SEED, 31, p]).choice(idx, n, replace=False)] = True
    rows = rows_all[keep]
    grp = grp_all[rows]
    rowid = assign_parts(grp, P + 1, (SEED, 32), rows)

    def blocks():
        for s in range(0, len(rows), BLOCK):
            pos = np.arange(s, min(s + BLOCK, len(rows)))
            blk = np.asarray(X[rows[pos]], dtype=np.float64)
            lib = blk.sum(axis=1, keepdims=True)
            yield pos, np.log1p(1e4 * blk[:, gcols] / lib)

    full, parts = accumulate(blocks(), grp, rowid, P, G)
    fm, pm, n_full, n_parts = means_from_sums(full, parts, grp, rowid, P, G)
    g = gates(fm, pm, n_full, n_parts, P, genes, perts, "GWPS_at_VIP_depth")
    return fm, pm, n_full, n_parts, g


def save(prefix, fm, pm, n_full, n_parts, P):
    np.save(OUT / f"{prefix}_delta.npy", (fm[:P] - fm[P]).astype(np.float32))
    np.save(OUT / f"{prefix}_pert_part_means.npy", pm[:, :, :P].astype(np.float32))
    np.save(OUT / f"{prefix}_ctrl_part_means.npy", pm[:, :, P].astype(np.float32))
    np.save(OUT / f"{prefix}_cell_counts.npy", n_full[:P])


def main() -> None:
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    perts = (AX / "shared_perturbations.txt").read_text().split()
    vip_labels = (AX / "viperturb_labels.txt").read_text().split()
    genes = (AX / "shared_genes.txt").read_text().split()
    vip_genes = (AX / "viperturb_gene_symbols.txt").read_text().split()
    P = len(perts)
    vip = build_viperturb(perts, vip_labels, genes, vip_genes)
    save("viperturb", *vip[:4], P)
    print(f"VIPerturb built [{time.time() - t0:.0f}s]", flush=True)
    gw = build_gwps_at_depth(perts, genes, vip[2][:P])
    save("gwps_vipdepth", *gw[:4], P)
    print(f"GWPS at VIPerturb depth built [{time.time() - t0:.0f}s]", flush=True)
    report = {
        "viperturb": vip[4],
        "gwps_vipdepth": gw[4],
        "median_cells_viperturb": float(np.median(vip[2][:P])),
        "median_cells_gwps_vipdepth": float(np.median(gw[2][:P])),
    }
    (OUT / "build_gates.json").write_text(json.dumps(report, indent=2))
    with open(OUT / "sha256.txt", "w") as fh:
        for p in sorted(OUT.glob("*.npy")):
            fh.write(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
