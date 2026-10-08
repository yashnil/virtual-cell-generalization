"""Freeze the N6 axes from identifiers and QC counts only (no expression values).

Inclusion (reports/n6_protocol.md §1): N5 panel perturbations whose VIPerturb label
(frozen alias map) has >= 30 cells passing the scPertEval cell filter (>= 200
detected genes), pooled over guides. Genes: N5 genes present (via the same alias
map) on the VIPerturb probe axis.

Reproduce: ``uv run python scripts/research_v3/make_n6_axes.py``
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
VIP = REPO / "data" / "raw" / "viperturb" / "viperturb_k562_genome_wide.h5ad"
AX5 = REPO / "data" / "splits" / "n5_k562"
OUT = REPO / "data" / "splits" / "n6_k562"
ALIAS = {"GET1": "WRB", "GET3": "ASNA1", "MICOS10": "MINOS1"}  # HGNC current -> VIPerturb symbol
MIN_CELLS = 30
MIN_GENES = 200
CONTROL = "NO-TARGET"


def col(obs, name):
    o = obs[name]
    if isinstance(o, h5py.Group):
        return o["categories"][:].astype(str)[o["codes"][:]]
    return o[:].astype(str)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    f = h5py.File(VIP, "r")
    gene = col(f["obs"], "gene")
    guide = col(f["obs"], "guide")
    nnz = np.diff(f["X"]["indptr"][:])
    var = f["var"][f["var"].attrs["_index"]][:].astype(str)
    qc = nnz >= MIN_GENES
    vc = pd.Series(gene[qc]).value_counts()
    p5 = (AX5 / "shared_perturbations.txt").read_text().split()
    g5 = (AX5 / "shared_genes.txt").read_text().split()
    perts, labels, rows = [], [], []
    for p in p5:
        lab = ALIAS.get(p, p)
        n = int(vc.get(lab, 0))
        guides = sorted(set(guide[qc & (gene == lab)]))
        rows.append(
            {
                "perturbation": p,
                "viperturb_label": lab,
                "cells_qc": n,
                "n_guides": len(guides),
                "alias_used": p in ALIAS,
                "included": n >= MIN_CELLS,
            }
        )
        if n >= MIN_CELLS:
            perts.append(p)
            labels.append(lab)
    vset = set(var)

    # gene axis: current symbol if present, else the frozen alias
    # (VIPerturb perturbation labels and its probe axis use different symbol vintages)
    def vip_symbol(g):
        if g in vset:
            return g
        return ALIAS[g] if g in ALIAS and ALIAS[g] in vset else None

    genes = [g for g in g5 if vip_symbol(g) is not None]
    vip_genes = [vip_symbol(g) for g in genes]
    pd.DataFrame(rows).to_csv(OUT / "panel_audit.csv", index=False)
    (OUT / "shared_perturbations.txt").write_text("\n".join(perts) + "\n")
    (OUT / "viperturb_labels.txt").write_text("\n".join(labels) + "\n")
    (OUT / "shared_genes.txt").write_text("\n".join(genes) + "\n")
    (OUT / "viperturb_gene_symbols.txt").write_text("\n".join(vip_genes) + "\n")
    summary = {
        "n_cells_total": int(len(gene)),
        "n_cells_qc": int(qc.sum()),
        "controls_qc": int(vc.get(CONTROL, 0)),
        "n_vip_labels": int(len(set(gene)) - 1),
        "n5_panel": len(p5),
        "n6_panel": len(perts),
        "n5_panel_present_any": int(sum(r["cells_qc"] > 0 for r in rows)),
        "median_cells_included": float(np.median([r["cells_qc"] for r in rows if r["included"]])),
        "guides_per_included_perturbation": pd.Series(
            [r["n_guides"] for r in rows if r["included"]]
        )
        .value_counts()
        .to_dict(),
        "n6_genes": len(genes),
        "n5_genes": len(g5),
        "vip_probe_genes": len(var),
        "alias_map": ALIAS,
    }
    (OUT / "axes_summary.json").write_text(json.dumps(summary, indent=2, default=int))
    with open(OUT / "sha256.txt", "w") as fh:
        for p in sorted(OUT.glob("*")):
            if p.name != "sha256.txt":
                fh.write(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n")
    print(json.dumps(summary, indent=2, default=int))


if __name__ == "__main__":
    main()
