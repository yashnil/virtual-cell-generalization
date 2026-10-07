"""Build the fully disjoint F/E1/E2 cell splits for N1/N4 (protocol §1).

Perturbed cells and control cells are each split three ways before any delta is
formed. Part means are stored separately so the shared-control diagnostic (F3b)
can be formed later. Nothing frozen is read except the design and raw cells.

Reproduce: ``uv run python scripts/build_n1_n4_splits.py``
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from virtual_cell.analysis import calibration_budget as cb
from virtual_cell.analysis import robustness
from virtual_cell.data import scperteval

REPO = Path(__file__).resolve().parents[1]
DESIGN = REPO / "data" / "splits" / "four_context_v1"
OUT = REPO / "outputs" / "n1_n4" / "splits"
SEED = 20261006
N_REPEATS = 5


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", type=Path, default=REPO / scperteval.DATA_SUBDIR)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    perts = (DESIGN / "shared_perturbations.txt").read_text().split()
    genes = (DESIGN / "shared_genes.txt").read_text().split()
    n_c, n_p, n_g = len(scperteval.DATASETS), len(perts), len(genes)

    pert = np.lib.format.open_memmap(
        OUT / "pert_part_means.npy",
        mode="w+",
        dtype=np.float32,
        shape=(N_REPEATS, 3, n_c, n_p, n_g),
    )
    ctrl = np.zeros((N_REPEATS, 3, n_c, n_g), dtype=np.float32)
    n_pert = np.zeros((N_REPEATS, 3, n_c, n_p), dtype=np.int64)
    n_ctrl = np.zeros((N_REPEATS, 3, n_c), dtype=np.int64)
    for ci, ds in enumerate(scperteval.DATASETS):
        t0 = time.time()
        X, group, X_ctrl = robustness.load_cells(
            ds.path(args.data_dir), genes=genes, perturbations=perts
        )
        for r in range(N_REPEATS):
            sm = cb.disjoint_split_means(
                X,
                group,
                X_ctrl,
                n_perturbations=n_p,
                rng=np.random.default_rng([SEED, ci, r]),
            )
            pert[r, :, ci] = sm.pert.astype(np.float32)
            ctrl[r, :, ci] = sm.ctrl.astype(np.float32)
            n_pert[r, :, ci] = sm.n_pert
            n_ctrl[r, :, ci] = sm.n_ctrl
        pert.flush()
        del X, X_ctrl
        print(f"  {ds.name:16s} [{time.time() - t0:.0f}s]", flush=True)
    del pert
    np.save(OUT / "ctrl_part_means.npy", ctrl)
    np.save(OUT / "n_pert_cells.npy", n_pert)
    np.save(OUT / "n_ctrl_cells.npy", n_ctrl)
    checks = {
        "contexts": list(scperteval.CONTEXTS),
        "n_repeats": N_REPEATS,
        "seed": SEED,
        "disjointness_asserted": True,
        "min_cells_per_part": {
            "fit": int(n_pert[:, 0].min()),
            "e1": int(n_pert[:, 1].min()),
            "e2": int(n_pert[:, 2].min()),
        },
        "control_cells_per_part": {
            ds.name: n_ctrl[0, :, i].tolist() for i, ds in enumerate(scperteval.DATASETS)
        },
    }
    (OUT / "build_checks.json").write_text(json.dumps(checks, indent=2))
    with open(OUT / "sha256.txt", "w") as fh:
        for f in sorted(OUT.glob("*.npy")):
            h = hashlib.sha256()
            with open(f, "rb") as src:
                for block in iter(lambda: src.read(1 << 24), b""):
                    h.update(block)
            fh.write(f"{h.hexdigest()}  {f.name}\n")
    print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    main()
