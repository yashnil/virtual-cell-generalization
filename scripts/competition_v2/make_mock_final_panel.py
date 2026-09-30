"""Build a MOCK final panel from the validation controls. LOCAL TESTING ONLY; never submit.

    uv run python scripts/competition_v2/make_mock_final_panel.py \
        --out outputs/final_mock/controls [--n 111] [--seed 20261022]

* Copies of ``context_{A,B,C}.h5ad`` are relabelled D, E, F. Only ``obs['context']``
  changes; X, genes and the other obs columns are byte-for-byte the same values.
* ``pert_counts.csv`` is a stratified random subset of the validation targets, in a
  shuffled order. It includes unsupported, single-source, two-source and 3+-source
  targets (by the C1 license-clean coverage).
* ``manifest.json`` is marked ``partition: MOCK-LOCAL-TEST-ONLY`` and ``do_not_submit``.
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
VAL = ROOT / "data" / "raw" / "arc2026" / "controls"
SUPPORT = ROOT / "data" / "splits" / "arc_target_support_license_clean_v1.csv"
RELABEL = {"A": "D", "B": "E", "C": "F"}
STRATA = {0: 4, 1: 25, 2: 70, 3: 12}  # sources per target -> how many (3 means 3+)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(ROOT / "outputs" / "final_mock" / "controls"))
    ap.add_argument("--seed", type=int, default=20261022)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    sup = pd.read_csv(SUPPORT)
    n_src = sup.n_usable_contexts_c1.clip(upper=3)
    chosen = []
    for k, n in STRATA.items():
        pool = sup.arc_target[n_src == k].to_numpy()
        chosen += list(rng.choice(pool, min(n, len(pool)), replace=False))
    targets = rng.permutation(np.array(chosen, dtype=object))
    pd.DataFrame({"target_gene": targets}).to_csv(out / "pert_counts.csv", index=False)
    shutil.copyfile(VAL / "gene_names.csv", out / "gene_names.csv")
    manifest = json.loads((VAL / "manifest.json").read_text())
    per = {}
    for old, new in RELABEL.items():
        a = ad.read_h5ad(VAL / f"context_{old}.h5ad")
        a.obs["context"] = pd.Categorical([new] * a.n_obs)
        a.obs_names = pd.Index([f"{new}_{x}" for x in a.obs_names], dtype=object)
        a.var_names = pd.Index(a.var_names.astype(str), dtype=object)
        for c in a.obs.columns:
            if c != "context":
                a.obs[c] = pd.Categorical(a.obs[c].astype(str))
        a.write_h5ad(out / f"context_{new}.h5ad")
        per[new] = {
            "n_perturbations": len(targets),
            "control_cells": int(a.n_obs),
            "ground_truth_cells": None,
            "n_ntc_ids": manifest["per_context"][old]["n_ntc_ids"],
        }
    mock = {
        **{
            k: manifest[k]
            for k in ("season", "pert_col", "context_col", "control_label", "n_genes")
        },
        "partition": "MOCK-LOCAL-TEST-ONLY",
        "panel_id": "mock-final-DEF",
        "do_not_submit": True,
        "contexts": list(RELABEL.values()),
        "n_constructs": len(targets),
        "per_context": per,
        "cells_per_pert": manifest["cells_per_pert"],
        "built_from": {
            "validation_panel": manifest["panel_id"],
            "seed": args.seed,
            "strata": STRATA,
        },
    }
    (out / "manifest.json").write_text(json.dumps(mock, indent=2) + "\n")
    (out / "README.md").write_text(
        "MOCK final panel built from the validation controls, for local integration tests "
        "only. **Never submit anything built from it.**\n"
    )
    counts = n_src.set_axis(sup.arc_target).loc[targets].value_counts().sort_index().to_dict()
    print(json.dumps({"out": str(out), "n_targets": len(targets), "by_c1_sources": counts}))


if __name__ == "__main__":
    main()
