"""EXPLORATORY (not preregistered): does averaging halves over repeats inflate agreement?

Recomputes N4's primary quality two ways on the same disjoint splits:
(a) per-repeat unbiased energies averaged over repeats (the N4 implementation);
(b) halves averaged over repeats first, then energies (the construction used in
transferability_confidence_model_v1). Reports Spearman(agreement, −D) and the
stable fraction for each. No pass/fail rule depends on this.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from virtual_cell.analysis import agreement_null as an
from virtual_cell.analysis import loco
from virtual_cell.data import scperteval
from virtual_cell.modelling import pathway_residual as pr

REPO = Path(__file__).resolve().parents[1]
D = np.load(REPO / "data/processed/four_context_v1/delta_tensor.npy").astype(np.float64)
pert = np.load(REPO / "outputs/n1_n4/splits/pert_part_means.npy", mmap_mode="r")
ctrl = np.load(REPO / "outputs/n1_n4/splits/ctrl_part_means.npy")
ms = float(
    json.loads((REPO / "outputs/transferability_v1/stability_rule.json").read_text())[
        "min_signal_energy"
    ]
)
rows = []
for fold in loco.make_folds(scperteval.CONTEXTS):
    t, src = fold.target_index, list(fold.source_indices)
    B = pr.baseline(D, src, pr.fit_scale(D, src))
    a = an.agreement(D[src])
    P = np.asarray(pert[:, :, t], dtype=np.float64)
    C = ctrl[:, :, t].astype(np.float64)
    h1 = [P[r, 0] - C[r, 0] for r in range(5)]
    h2 = [((P[r, 1] - C[r, 1]) + (P[r, 2] - C[r, 2])) / 2 for r in range(5)]
    q_new, _ = an.quality_d(h1, h2, B, min_signal=ms)
    q_old, _ = an.quality_d([np.mean(h1, axis=0)], [np.mean(h2, axis=0)], B, min_signal=ms)
    rows.append(
        {
            "context": scperteval.dataset(scperteval.CONTEXTS[t]).cell_line,
            "rho_per_repeat": an.spearman(a, q_new),
            "stable_per_repeat": float(np.isfinite(q_new).mean()),
            "rho_averaged_halves": an.spearman(a, q_old),
            "stable_averaged_halves": float(np.isfinite(q_old).mean()),
        }
    )
df = pd.DataFrame(rows)
df.to_csv(REPO / "outputs/n1_n4/n4/exploratory_halves_bias.csv", index=False)
print(df.round(3).to_string(index=False))
