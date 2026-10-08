"""Build the N3 six-context data object (reports/n3_protocol.md §0).

Contexts 0–3 are subset from frozen scPertEval objects (canonical tensor, N1
split-part means). Contexts 4–5 (X-Atlas HCT116, HEK293T) are streamed from the
Lance expression tables once. Each cell's log1p(CP10K) values on the shared genes
are accumulated into canonical, full-control and disjoint F/E1/E2 part sums for 5
repeats, using the N1 split algorithm. Gates G1–G4 stop the build on failure.

Research use only (reports/xatlas_license_memo.md). Outputs stay in git-ignored
``outputs/n3/data``.

Reproduce: ``uv run python scripts/research_v3/build_n3_six_context.py``
"""

from __future__ import annotations

import hashlib
import json
import time
from multiprocessing import Pool
from pathlib import Path

import lance
import numpy as np

REPO = Path(__file__).resolve().parents[2]
XATLAS = REPO / "data" / "raw" / "competition_v2" / "xatlas_orion" / "data"
DESIGN4 = REPO / "data" / "splits" / "four_context_v1"
DESIGN6 = REPO / "data" / "splits" / "six_context_n3"
CANON = REPO / "data" / "processed" / "four_context_v1"
N1SPLITS = REPO / "outputs" / "n1_n4" / "splits"
OUT = REPO / "outputs" / "n3" / "data"
TMP = REPO / "outputs" / "n3" / "tmp"
SEED = 20261006
N_REPEATS = 5
CONTROL = "Non-Targeting"
LINES = ("HCT116", "HEK293T")
CONTEXTS = ("K562", "RPE1", "HepG2", "Jurkat") + LINES
N_WORKERS = 8


def three_way(idx: np.ndarray):
    n = len(idx)
    h, q = n // 2, n // 4
    return idx[:h], idx[h : h + q], idx[h + q : h + 2 * q]


def axes() -> tuple[list[str], list[str]]:
    p4 = (DESIGN4 / "shared_perturbations.txt").read_text().split()
    g4 = (DESIGN4 / "shared_genes.txt").read_text().split()
    keep_p, keep_g = set(p4), set(g4)
    for line in LINES:
        cells = lance.dataset(f"{XATLAS}/{line}/cells.lance").to_table(columns=["gene_target"])
        vc = cells.to_pandas().gene_target.astype(str).value_counts()
        keep_p &= set(vc[vc >= 30].index)
        gid = lance.dataset(f"{XATLAS}/{line}/genes.lance").to_table(columns=["gene_id"])
        gs = gid.to_pandas().gene_id.astype(str)
        keep_g &= set(gs[~gs.duplicated(keep=False)])
    return [p for p in p4 if p in keep_p], [g for g in g4 if g in keep_g]


# --------------------------------------------------------------------------
# X-Atlas streaming
# --------------------------------------------------------------------------


def line_layout(line: str, perts: list[str], genes: list[str], ci: int) -> dict:
    """Cell grouping and split assignments; reads identifiers and library sizes only."""
    cells = (
        lance.dataset(f"{XATLAS}/{line}/cells.lance")
        .to_table(columns=["cell_integer_id", "gene_target", "total_counts", "n_genes_by_counts"])
        .to_pandas()
    )
    cells["gene_target"] = cells.gene_target.astype(str)
    row_of = {p: i for i, p in enumerate(perts)}
    P = len(perts)
    grp = cells.gene_target.map(lambda x: row_of.get(x, P if x == CONTROL else -1)).to_numpy()
    sel = grp >= 0
    cells = cells[sel].sort_values("cell_integer_id").reset_index(drop=True)
    grp = cells.gene_target.map(lambda x: row_of.get(x, P)).to_numpy().astype(np.int64)
    ids = cells.cell_integer_id.to_numpy().astype(np.int64)
    n = len(ids)
    # split assignment: rowid[r, i] = part * (P + 1) + group, or -1 (dropped remainder)
    rowid = np.full((N_REPEATS, n), -1, dtype=np.int64)
    order = np.argsort(grp, kind="stable")  # within group, ascending cell_integer_id
    bounds = np.searchsorted(grp[order], np.arange(P + 2))
    for r in range(N_REPEATS):
        rng = np.random.default_rng([SEED, ci, r])
        for p in range(P):
            idx = order[bounds[p] : bounds[p + 1]].copy()
            rng.shuffle(idx)
            for j, part in enumerate(three_way(idx)):
                rowid[r, part] = j * (P + 1) + p
        ctrl = order[bounds[P] : bounds[P + 1]]
        perm = ctrl[rng.permutation(len(ctrl))]
        for j, part in enumerate(three_way(perm)):
            rowid[r, part] = j * (P + 1) + P
    # G4: every selected cell gets at most one part per repeat (single array slot), and
    # the number assigned per group equals the three-way sizes
    for r in range(N_REPEATS):
        assigned = np.bincount(grp[rowid[r] >= 0], minlength=P + 1)
        sizes = np.bincount(grp, minlength=P + 1)
        if not np.array_equal(assigned, sizes // 2 + 2 * (sizes // 4)):
            raise AssertionError("G4: part sizes inconsistent")
    gtab = lance.dataset(f"{XATLAS}/{line}/genes.lance").to_table(
        columns=["gene_id", "gene_integer_id"]
    )
    gdf = gtab.to_pandas()
    gene_of = {g: i for i, g in enumerate(genes)}
    gene_map = np.full(int(gdf.gene_integer_id.max()) + 1, -1, dtype=np.int64)
    for g, gi in zip(gdf.gene_id.astype(str), gdf.gene_integer_id, strict=True):
        gene_map[int(gi)] = gene_of.get(g, -1)
    lookup = np.full(int(ids.max()) + 1, -1, dtype=np.int64)
    lookup[ids] = np.arange(n)
    return {
        "ids": ids,
        "grp": grp,
        "total": cells.total_counts.to_numpy().astype(np.float64),
        "ngenes": cells.n_genes_by_counts.to_numpy().astype(np.int64),
        "rowid": rowid,
        "gene_map": gene_map,
        "lookup": lookup,
        "P": P,
        "G": len(genes),
    }


def _scan(args) -> str:
    line, frag_ids, layout_path, wid = args
    lay = dict(np.load(layout_path, allow_pickle=False))
    P, G = int(lay["P"]), int(lay["G"])
    full = np.zeros((P + 1) * G)
    parts = np.zeros((N_REPEATS, 3 * (P + 1) * G))
    entries = 0
    ds = lance.dataset(f"{XATLAS}/{line}/expression.lance")
    frags = {f.fragment_id: f for f in ds.get_fragments()}
    lookup, gene_map = lay["lookup"], lay["gene_map"]
    for fid in frag_ids:
        for attempt in range(6):
            try:
                tab = frags[fid].to_table(columns=["cell_integer_id", "gene_integer_id", "value"])
                break
            except OSError:
                if attempt == 5:
                    raise
                time.sleep(5 * (attempt + 1))
        ci = tab.column("cell_integer_id").to_numpy().astype(np.int64)
        ok = ci < len(lookup)
        loc = np.full(len(ci), -1, dtype=np.int64)
        loc[ok] = lookup[ci[ok]]
        sel = loc >= 0
        entries += int(sel.sum())
        if not sel.any():
            continue
        loc = loc[sel]
        gcol = gene_map[tab.column("gene_integer_id").to_numpy()[sel].astype(np.int64)]
        keep = gcol >= 0
        loc, gcol = loc[keep], gcol[keep]
        val = tab.column("value").to_numpy()[sel][keep].astype(np.float64)
        val = np.log1p(1e4 * val / lay["total"][loc])
        full += np.bincount(lay["grp"][loc] * G + gcol, weights=val, minlength=(P + 1) * G)
        for r in range(N_REPEATS):
            rid = lay["rowid"][r][loc]
            m = rid >= 0
            parts[r] += np.bincount(rid[m] * G + gcol[m], weights=val[m], minlength=3 * (P + 1) * G)
    out = TMP / f"{line}_w{wid}.npz"
    np.savez(out, full=full, parts=parts, entries=entries)
    return str(out)


def build_line(line: str, ci: int, perts, genes) -> dict:
    t0 = time.time()
    lay = line_layout(line, perts, genes, ci)
    TMP.mkdir(parents=True, exist_ok=True)
    layout_path = TMP / f"{line}_layout.npz"
    np.savez(layout_path, **{k: np.asarray(v) for k, v in lay.items()})
    frag_ids = [
        f.fragment_id for f in lance.dataset(f"{XATLAS}/{line}/expression.lance").get_fragments()
    ]
    jobs = [(line, frag_ids[w::N_WORKERS], str(layout_path), w) for w in range(N_WORKERS)]
    with Pool(N_WORKERS) as pool:
        paths = pool.map(_scan, jobs)
    P, G = lay["P"], lay["G"]
    full = np.zeros((P + 1) * G)
    parts = np.zeros((N_REPEATS, 3 * (P + 1) * G))
    entries = 0
    for p in paths:
        z = np.load(p)
        full += z["full"]
        parts += z["parts"]
        entries += int(z["entries"])
        Path(p).unlink()
    coverage = entries / lay["ngenes"].sum()
    if coverage < 0.99:  # G1
        raise RuntimeError(f"G1 failed for {line}: coverage {coverage:.4f}")
    n_full = np.bincount(lay["grp"], minlength=P + 1).astype(np.float64)
    full_means = full.reshape(P + 1, G) / n_full[:, None]
    n_parts = np.zeros((N_REPEATS, 3 * (P + 1)))
    for r in range(N_REPEATS):
        rid = lay["rowid"][r]
        n_parts[r] = np.bincount(rid[rid >= 0], minlength=3 * (P + 1))
    part_means = (
        parts.reshape(N_REPEATS, 3, P + 1, G) / n_parts.reshape(N_REPEATS, 3, P + 1)[..., None]
    )
    print(f"  {line}: scan done, coverage {coverage:.5f} [{time.time() - t0:.0f}s]", flush=True)
    return {
        "delta": full_means[:P] - full_means[P],
        "ctrl_full": full_means[P],
        "pert_parts": part_means[:, :, :P],
        "ctrl_parts": part_means[:, :, P],
        "n_pert": n_parts.reshape(N_REPEATS, 3, P + 1)[:, :, :P],
        "n_ctrl": n_parts.reshape(N_REPEATS, 3, P + 1)[:, :, P],
        "n_full": n_full[:P],
        "coverage": coverage,
        "full_sums": full.reshape(P + 1, G),
        "part_sums": parts.reshape(N_REPEATS, 3, P + 1, G),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    DESIGN6.mkdir(parents=True, exist_ok=True)
    perts, genes = axes()
    (DESIGN6 / "shared_perturbations.txt").write_text("\n".join(perts) + "\n")
    (DESIGN6 / "shared_genes.txt").write_text("\n".join(genes) + "\n")
    (DESIGN6 / "contexts.txt").write_text("\n".join(CONTEXTS) + "\n")
    print(f"axes: {len(perts)} perturbations x {len(genes)} genes", flush=True)
    p4 = (DESIGN4 / "shared_perturbations.txt").read_text().split()
    g4 = (DESIGN4 / "shared_genes.txt").read_text().split()
    pi = np.array([p4.index(p) for p in perts])
    gi = np.array([g4.index(g) for g in genes])
    P, G = len(perts), len(genes)

    delta = np.zeros((6, P, G), dtype=np.float32)
    ctrl_full = np.zeros((6, G), dtype=np.float32)
    counts = np.zeros((6, P))
    D4 = np.load(CANON / "delta_tensor.npy", mmap_mode="r")
    C4 = np.load(CANON / "control_means.npy")
    N4 = np.load(CANON / "cell_counts.npy")
    for c in range(4):
        delta[c] = np.asarray(D4[c])[np.ix_(pi, gi)]
        ctrl_full[c] = C4[c, gi]
        counts[c] = N4[c, pi]
    pert = np.lib.format.open_memmap(
        OUT / "pert_part_means.npy", mode="w+", dtype=np.float32, shape=(N_REPEATS, 3, 6, P, G)
    )
    ctrl = np.zeros((N_REPEATS, 3, 6, G), dtype=np.float32)
    n_pert = np.zeros((N_REPEATS, 3, 6, P), dtype=np.int64)
    n_ctrl = np.zeros((N_REPEATS, 3, 6), dtype=np.int64)
    pm4 = np.load(N1SPLITS / "pert_part_means.npy", mmap_mode="r")
    cm4 = np.load(N1SPLITS / "ctrl_part_means.npy")
    np4 = np.load(N1SPLITS / "n_pert_cells.npy")
    nc4 = np.load(N1SPLITS / "n_ctrl_cells.npy")
    for c in range(4):
        for r in range(N_REPEATS):
            for j in range(3):
                pert[r, j, c] = np.asarray(pm4[r, j, c])[np.ix_(pi, gi)]
        ctrl[:, :, c] = cm4[:, :, c][..., gi]
        n_pert[:, :, c] = np4[:, :, c][..., pi]
        n_ctrl[:, :, c] = nc4[:, :, c]

    gates = {"G1_coverage": {}, "G2_max_rel_err": {}, "G3_median_on_target": {}, "G4": "asserted"}
    for j, line in enumerate(LINES):
        c = 4 + j
        res = build_line(line, c, perts, genes)
        delta[c] = res["delta"].astype(np.float32)
        ctrl_full[c] = res["ctrl_full"].astype(np.float32)
        counts[c] = res["n_full"]
        pert[:, :, c] = res["pert_parts"].astype(np.float32)
        ctrl[:, :, c] = res["ctrl_parts"].astype(np.float32)
        n_pert[:, :, c] = res["n_pert"]
        n_ctrl[:, :, c] = res["n_ctrl"]
        gates["G1_coverage"][line] = res["coverage"]
        # G2: canonical = cell-weighted combination of part sums + dropped remainder is
        # excluded, so compare sums over cells that are in some part against full sums
        # minus nothing: check the identity on repeat 0 for groups with no remainder.
        n_parts_tot = res["n_pert"][0].sum(axis=0)
        exact = n_parts_tot == res["n_full"]
        recon = res["part_sums"][0, :, :P].sum(axis=0)[exact] / n_parts_tot[exact, None]
        target = res["full_sums"][:P][exact] / res["n_full"][exact, None]
        err = float(np.max(np.abs(recon - target)) / max(1e-12, np.max(np.abs(target))))
        gates["G2_max_rel_err"][line] = {"max_rel_err": err, "n_checked": int(exact.sum())}
        if err > 1e-6:
            raise RuntimeError(f"G2 failed for {line}: {err}")
        on = [res["delta"][i, genes.index(p)] for i, p in enumerate(perts) if p in genes]
        med = float(np.median(on))
        gates["G3_median_on_target"][line] = {"median": med, "n": len(on)}
        if not med < 0:
            raise RuntimeError(f"G3 failed for {line}: median on-target delta {med}")
        pert.flush()
    del pert
    np.save(OUT / "delta6.npy", delta)
    np.save(OUT / "control_means6.npy", ctrl_full)
    np.save(OUT / "cell_counts6.npy", counts)
    np.save(OUT / "ctrl_part_means.npy", ctrl)
    np.save(OUT / "n_pert_cells.npy", n_pert)
    np.save(OUT / "n_ctrl_cells.npy", n_ctrl)
    gates["axes"] = {"n_perturbations": P, "n_genes": G, "contexts": list(CONTEXTS)}
    gates["min_cells_per_part"] = {
        CONTEXTS[c]: n_pert[:, :, c].min(axis=(0, 2)).tolist() for c in range(6)
    }
    gates["control_cells_per_part"] = {CONTEXTS[c]: n_ctrl[0, :, c].tolist() for c in range(6)}
    (OUT / "build_gates.json").write_text(json.dumps(gates, indent=2))
    with open(DESIGN6 / "sha256.txt", "w") as fh:
        for f in sorted(DESIGN6.glob("*.txt")):
            if f.name != "sha256.txt":
                fh.write(f"{hashlib.sha256(f.read_bytes()).hexdigest()}  {f.name}\n")
    with open(OUT / "sha256.txt", "w") as fh:
        for f in sorted(OUT.glob("*.npy")):
            h = hashlib.sha256()
            with open(f, "rb") as src:
                for block in iter(lambda: src.read(1 << 24), b""):
                    h.update(block)
            fh.write(f"{h.hexdigest()}  {f.name}\n")
    print(json.dumps(gates, indent=2))


if __name__ == "__main__":
    main()
