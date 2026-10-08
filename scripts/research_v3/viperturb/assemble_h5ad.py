"""Assemble and validate the VIPerturb-seq genome-wide K562 h5ad from R exports.

Each bin's export (export_seurat.R) is rebuilt as a cells × genes CSR matrix of
integer counts and validated against the statistics R computed on the source
object: dimensions, nnz, total counts, sampled column/row sums and sampled
entries. Bins must share the gene axis. Cells present in several bins (the
non-targeting controls, by design) must be identical across bins; they are kept
once. Writes ``viperturb_k562_genome_wide.h5ad`` + SHA-256 + a validation JSON.

Usage: uv run python scripts/research_v3/viperturb/assemble_h5ad.py
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

REPO = Path(__file__).resolve().parents[3]
EXP = REPO / "data" / "raw" / "viperturb" / "export"
OUT = REPO / "data" / "raw" / "viperturb" / "viperturb_k562_genome_wide.h5ad"
VAL = REPO / "data" / "provenance" / "research_v3" / "viperturb_conversion_validation.json"
BINS = ("genome_wide_binA", "genome_wide_binB", "genome_wide_binC")


def load_bin(name: str) -> tuple[sparse.csr_matrix, list[str], list[str], pd.DataFrame, dict, dict]:
    d = EXP / name
    info = json.loads((d / "object.json").read_text())
    i = np.fromfile(d / "i.int32", dtype="<i4")
    p = np.fromfile(d / "p.int32", dtype="<i4")
    x = np.fromfile(d / "x.int32", dtype="<i4")
    n_genes, n_cells = info["dims_genes_cells"]
    csc = sparse.csc_matrix((x, i, p), shape=(n_genes, n_cells))
    genes = (d / "genes.txt").read_text().splitlines()
    cells = (d / "barcodes.txt").read_text().splitlines()
    meta = pd.read_csv(d / "meta.csv", index_col=0)
    checks = {
        "dims": [n_genes, n_cells] == [len(genes), len(cells)] == list(csc.shape),
        "nnz": int(csc.nnz) == int(info["nnz"]),
        "total_counts": int(x.astype(np.int64).sum()) == int(round(info["total_counts"])),
        "col_sums": bool(
            np.array_equal(
                np.asarray(csc[:, np.array(info["check_cols"]) - 1].sum(axis=0)).ravel(),
                np.array(info["check_col_sums"]),
            )
        ),
        "row_sums": bool(
            np.array_equal(
                np.asarray(csc[np.array(info["check_rows"]) - 1].sum(axis=1)).ravel(),
                np.array(info["check_row_sums"]),
            )
        ),
        "entries": bool(
            np.array_equal(
                np.asarray(
                    csc[
                        np.array(info["check_entries"])[:, 0] - 1,
                        np.array(info["check_entries"])[:, 1] - 1,
                    ]
                ).ravel(),
                np.array(info["check_entry_values"]),
            )
        ),
        "meta_rows_match_cells": list(meta.index.astype(str)) == cells,
        "integer_nonnegative": bool((x >= 0).all()),
    }
    return csc.T.tocsr(), genes, cells, meta, info, checks


def main() -> None:
    mats, metas, report = [], [], {"bins": {}}
    genes0 = None
    seen: dict[str, tuple[int, int]] = {}
    dup_identical, dup_total = 0, 0
    keep_rows = []
    for b, name in enumerate(BINS):
        X, genes, cells, meta, info, checks = load_bin(name)
        if genes0 is None:
            genes0 = genes
        checks["gene_axis_identical_to_binA"] = genes == genes0
        report["bins"][name] = {
            "checks": checks,
            "dims_genes_cells": info["dims_genes_cells"],
            "nnz": info["nnz"],
            "total_counts": info["total_counts"],
            "md5_input": info["md5_input"],
            "class": info["class"],
            "assays": info["assays"],
            "layers": info["layers"],
            "meta_columns": info["meta_columns"],
            "r_version": info["r_version"],
            "SeuratObject": info["SeuratObject"],
        }
        if not all(checks.values()):
            raise RuntimeError(f"{name}: validation failed {checks}")
        keep = np.ones(len(cells), dtype=bool)
        for r, c in enumerate(cells):
            if c in seen:
                dup_total += 1
                pb, pr = seen[c]
                same = (X[r] != mats[pb][pr]).nnz == 0
                dup_identical += int(same)
                keep[r] = False
            else:
                seen[c] = (b, r)
        mats.append(X)
        meta = meta.assign(source_bin=name)
        metas.append(meta[keep])
        keep_rows.append(np.flatnonzero(keep))
    if dup_identical != dup_total:
        raise RuntimeError(f"{dup_total - dup_identical} duplicated barcodes differ across bins")
    X = sparse.vstack([m[k] for m, k in zip(mats, keep_rows, strict=True)], format="csr")
    obs = pd.concat(metas)
    obs.index = obs.index.astype(str)
    if obs.index.duplicated().any():
        raise RuntimeError("duplicate barcodes after de-duplication")
    adata = ad.AnnData(
        X=X.astype(np.int32), obs=obs, var=pd.DataFrame(index=pd.Index(genes0, name="gene"))
    )
    adata.uns["provenance"] = {
        "zenodo_record": "10.5281/zenodo.18460279",
        "licence": "CC BY 4.0",
        "publication": "Bradu et al., bioRxiv 2026.02.12.705613",
        "bins": list(BINS),
    }
    adata.write_h5ad(OUT, compression="gzip")
    h = hashlib.sha256()
    with open(OUT, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 24), b""):
            h.update(blk)
    report.update(
        {
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "nnz": int(X.nnz),
            "total_counts": int(X.sum()),
            "duplicated_barcodes": dup_total,
            "duplicated_identical": dup_identical,
            "h5ad": str(OUT.relative_to(REPO)),
            "h5ad_sha256": h.hexdigest(),
        }
    )
    VAL.write_text(json.dumps(report, indent=2))
    print(json.dumps({k: v for k, v in report.items() if k != "bins"}, indent=2))


if __name__ == "__main__":
    main()
