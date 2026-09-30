"""Exact comparison of two prediction .h5ad files: obs, var and every CSR array.

    uv run python scripts/competition_v2/compare_predictions_exact.py A.h5ad B.h5ad [--out J]

It reports whether the files are identical array for array. If they are not, it lists
the (context, target) blocks whose rows differ, the number of differing cells, and the
largest per-cell depth and count differences. Everything is streamed by row block, so no
full matrix is loaded.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse


def _col(f, name):
    node = f["obs"][name]
    if isinstance(node, h5py.Group):
        cats = np.array(
            [c.decode() if isinstance(c, bytes) else str(c) for c in node["categories"][()]]
        )
        return cats[node["codes"][()]]
    return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in node[()]])


def _var(f):
    var = f["var"]
    node = var[var.attrs.get("_index", "_index")]
    return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in node[()]])


def compare(a_path: Path, b_path: Path, block_rows: int = 400) -> dict:
    with h5py.File(a_path, "r") as a, h5py.File(b_path, "r") as b:
        out = {
            "a": str(a_path),
            "b": str(b_path),
            "shape_equal": tuple(a["X"].attrs["shape"]) == tuple(b["X"].attrs["shape"]),
            "var_equal": bool(np.array_equal(_var(a), _var(b))),
        }
        ta, tb = _col(a, "target_gene"), _col(b, "target_gene")
        ca, cb = _col(a, "context"), _col(b, "context")
        out["obs_equal"] = bool(np.array_equal(ta, tb) and np.array_equal(ca, cb))
        out["nnz"] = [int(a["X/data"].shape[0]), int(b["X/data"].shape[0])]
        ipa, ipb = a["X/indptr"][()], b["X/indptr"][()]
        n = len(ipa) - 1
        diff_blocks, diff_cells, max_depth, max_count, diff_entries = [], 0, 0.0, 0.0, 0
        for lo in range(0, n, block_rows):
            hi = min(lo + block_rows, n)
            sa = slice(int(ipa[lo]), int(ipa[hi]))
            sb = slice(int(ipb[lo]), int(ipb[hi]))
            same_ptr = np.array_equal(ipa[lo : hi + 1] - ipa[lo], ipb[lo : hi + 1] - ipb[lo])
            da, db = a["X/data"][sa], b["X/data"][sb]
            ia, ib = a["X/indices"][sa], b["X/indices"][sb]
            if same_ptr and np.array_equal(da, db) and np.array_equal(ia, ib):
                continue
            ncol = int(a["X"].attrs["shape"][1])
            ma = sparse.csr_matrix((da, ia, ipa[lo : hi + 1] - ipa[lo]), shape=(hi - lo, ncol))
            mb = sparse.csr_matrix((db, ib, ipb[lo : hi + 1] - ipb[lo]), shape=(hi - lo, ncol))
            d = (ma.astype(np.float64) - mb.astype(np.float64)).tocsr()
            d.eliminate_zeros()
            rows = np.diff(d.indptr) > 0
            if not rows.any():
                continue
            diff_cells += int(rows.sum())
            depth = np.abs(np.asarray(d.sum(axis=1)).ravel())
            max_depth = max(max_depth, float(depth.max()))
            max_count = max(max_count, float(np.abs(d.data).max()))
            diff_entries += int(d.nnz)
            diff_blocks.append({"context": str(ca[lo]), "target": str(ta[lo]), "rows": [lo, hi]})
        out.update(
            identical=out["shape_equal"]
            and out["var_equal"]
            and out["obs_equal"]
            and not diff_blocks,
            differing_blocks=diff_blocks,
            differing_cells=diff_cells,
            differing_entries=diff_entries,
            max_abs_depth_difference=max_depth,
            max_abs_count_difference=max_count,
        )
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--out")
    args = ap.parse_args()
    res = compare(Path(args.a), Path(args.b))
    if args.out:
        Path(args.out).write_text(json.dumps(res, indent=2) + "\n")
    summary = {k: v for k, v in res.items() if k != "differing_blocks"}
    summary["n_differing_blocks"] = len(res["differing_blocks"])
    summary["first_blocks"] = res["differing_blocks"][:10]
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
