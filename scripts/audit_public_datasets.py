"""Local audit of every public scPertEval dataset, and its Arc coverage.

Reads the downloaded ``.h5ad`` files directly -- no network, no cached
summaries -- and rebuilds ``public_label_inventory.json`` from what is on disk.
Prints the seven-dataset Arc coverage table so it can be compared against the
previously reported values rather than assumed to match.

Usage::

    uv run python scripts/audit_public_datasets.py
"""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

from virtual_cell.data import arc2026

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "scperteval"
CONTROLS = ROOT / "data" / "raw" / "arc2026" / "controls"
# The v1 inventory is frozen (arc_bridge_v1_freeze.txt) and read remotely over
# HTTP range requests. This audit reads the downloaded files instead, so it
# writes a separate v2 rather than overwriting a frozen artifact. Every Arc
# coverage figure the two produce is identical; v2 adds the fields that only a
# local read can supply.
INVENTORY = ROOT / "data" / "provenance" / "scperteval" / "public_label_inventory_v2_local.json"

DATASETS = (
    "replogle22k562",
    "replogle22rpe1",
    "nadig25hepg2",
    "nadig25jurkat",
    "arch1",
    "kaden25rpe1",
    "wessels23",
)
CONTROL_LABELS = {"control", "non-targeting"}


def read_dataset(path: Path) -> dict:
    with h5py.File(path, "r") as handle:
        x = handle["X"]
        shape = tuple(int(n) for n in x.attrs["shape"]) if "shape" in x.attrs else x.shape
        encoding = str(x.attrs.get("encoding-type", "dense"))
        sample = x["data"][:300_000] if isinstance(x, h5py.Group) else np.asarray(x[:30]).ravel()
        obs = read_elem(handle["obs"])
        var_group = handle["var"]
        node = var_group[var_group.attrs.get("_index", "_index")]
        values = node["values"] if isinstance(node, h5py.Group) else node
        genes = [g.decode() if isinstance(g, bytes) else str(g) for g in values[:]]

    column = "perturbation" if "perturbation" in obs.columns else "gene"
    labels = obs[column].astype(str)
    counts = labels.value_counts()
    perturbations = sorted(set(labels) - CONTROL_LABELS)
    per_pert = counts.drop(list(CONTROL_LABELS), errors="ignore")

    return {
        "label_col": column,
        "cells": int(shape[0]),
        "genes": genes,
        "perturbations": perturbations,
        "control_cells": int(sum(counts.get(c, 0) for c in CONTROL_LABELS)),
        "x_encoding": encoding,
        "integer_valued": bool(np.all(sample == np.floor(sample))),
        "value_min": float(sample.min()),
        "value_max": float(sample.max()),
        "cells_per_pert_median": float(per_pert.median()),
        "cells_per_pert_min": int(per_pert.min()),
        "cells_per_pert_max": int(per_pert.max()),
        "n_perts_ge_400_cells": int((per_pert >= 400).sum()),
        "n_perts_ge_100_cells": int((per_pert >= 100).sum()),
    }


def main() -> None:
    arc_genes = set(arc2026.load_gene_names(CONTROLS))
    arc_targets = set(arc2026.load_pert_counts(CONTROLS)["target_gene"])

    inventory: dict[str, dict] = {}
    rows = []
    for name in DATASETS:
        path = RAW / f"{name}_processed_complete.h5ad"
        if not path.is_file():
            raise FileNotFoundError(f"{path} is missing; the audit must be complete")
        print(f"  reading {name} ...", flush=True)
        record = read_dataset(path)
        inventory[name] = record

        genes = set(record["genes"])
        perts = set(record["perturbations"])
        rows.append(
            {
                "dataset": name,
                "cells": record["cells"],
                "genes": len(record["genes"]),
                "perturbations": len(perts),
                "control_cells": record["control_cells"],
                "integer_valued": record["integer_valued"],
                "median_cells_per_pert": record["cells_per_pert_median"],
                "perts_ge_400_cells": record["n_perts_ge_400_cells"],
                "arc_genes_measured": len(arc_genes & genes),
                "pct_arc_genes": 100 * len(arc_genes & genes) / len(arc_genes),
                "arc_targets_perturbed": len(arc_targets & perts),
                "pct_arc_targets": 100 * len(arc_targets & perts) / len(arc_targets),
            }
        )

    table = pd.DataFrame(rows)
    print("\n" + "=" * 78)
    print("SEVEN-DATASET AUDIT, RECOMPUTED FROM LOCAL FILES")
    print("=" * 78)
    print(table.to_string(index=False, float_format=lambda v: f"{v:.2f}"))

    covered_g = {g for rec in inventory.values() for g in rec["genes"]} & arc_genes
    covered_t = {p for rec in inventory.values() for p in rec["perturbations"]} & arc_targets
    n_ctx = {
        t: sum(t in set(rec["perturbations"]) for rec in inventory.values()) for t in arc_targets
    }
    print(
        f"\n  union Arc genes measured      {len(covered_g):>6}/18533 "
        f"({100 * len(covered_g) / 18533:.1f}%)"
    )
    print(
        f"  union Arc targets perturbed   {len(covered_t):>6}/300   "
        f"({100 * len(covered_t) / 300:.1f}%)"
    )
    print(f"  Arc targets in >=2 datasets   {sum(v >= 2 for v in n_ctx.values()):>6}/300")

    table.to_csv(ROOT / "outputs" / "unseen_perturbation_v1" / "dataset_audit.csv", index=False)
    with INVENTORY.open("w") as fh:
        json.dump(inventory, fh, indent=1, sort_keys=True)
    print(f"\n  rebuilt {INVENTORY}")


if __name__ == "__main__":
    main()
