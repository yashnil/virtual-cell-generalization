"""EXPLORATORY (not preregistered): rank structure of E4's learned γ⊥.

Question: is E4's γ⊥ gain carried by a few target-specific response axes
(a context-specific program scaled per perturbation) or by many
perturbation-specific patterns? E4's γ⊥ prediction on the test set is truncated
to its top-r singular components (r = 1, 3, 10) and M3 is recomputed.
No pass/fail rule depends on this.

Reproduce: ``uv run python scripts/explore_n1_rank.py``
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from virtual_cell.analysis import calibration_budget as cb
from virtual_cell.analysis import loco
from virtual_cell.data import scperteval
from virtual_cell.modelling import pathway_residual as pr

REPO = Path(__file__).resolve().parents[1]
CANON = REPO / "data" / "processed" / "four_context_v1"
SPLITS = REPO / "outputs" / "n1_n4" / "splits"
OUT = REPO / "outputs" / "n1_n4" / "n1"
SEED = 20261006
CL = {c: scperteval.dataset(c).cell_line for c in scperteval.CONTEXTS}


def orth(x: np.ndarray, u: np.ndarray) -> np.ndarray:
    return x - np.einsum("pg,pg->p", x, u)[:, None] * u


def m3(ev: cb.EvalContext, P_orth: np.ndarray) -> float:
    a1, a2 = orth(ev.e1c, ev.u), orth(ev.e2c, ev.u)
    num = np.einsum("pg,pg->p", a1 - P_orth, a2 - P_orth).sum()
    den = np.einsum("pg,pg->p", a1, a2).sum()
    return 1 - num / den


def main() -> None:
    D = np.load(CANON / "delta_tensor.npy").astype(np.float64)
    pert = np.load(SPLITS / "pert_part_means.npy", mmap_mode="r")
    ctrl = np.load(SPLITS / "ctrl_part_means.npy")
    rows = []
    for fold in loco.make_folds(scperteval.CONTEXTS):
        t, src = fold.target_index, list(fold.source_indices)
        sv = cb.make_source_view(D, src, pr.fit_scale(D, src))
        perm = np.random.default_rng([SEED, 1, t]).permutation(D.shape[1])
        T, pool = np.sort(perm[:379]), np.sort(perm[379:])
        r = 0
        P = np.asarray(pert[r, :, t], dtype=np.float64)
        C = ctrl[r, :, t].astype(np.float64)
        fit = P[0] - C[0]
        ev = cb.make_eval_context((P[1] - C[1])[T], (P[2] - C[2])[T], sv.A[T])
        for k in (20, 50, 200, len(pool)):
            n_draw = 1 if k == len(pool) else 10
            for d in range(n_draw):
                K = (
                    pool
                    if k == len(pool)
                    else np.sort(
                        np.random.default_rng([SEED, 2, t, k, r, d]).choice(pool, k, replace=False)
                    )
                )
                Pe4, _ = cb.e4(sv, T, K, fit[K])
                Po = orth(Pe4 - Pe4.mean(axis=0), ev.u)
                U, S, Vt = np.linalg.svd(Po, full_matrices=False)
                energy = S**2 / np.sum(S**2)
                row = {
                    "context": CL[scperteval.CONTEXTS[t]],
                    "k": k,
                    "draw": d,
                    "M3_full": m3(ev, Po),
                }
                for rank in (1, 3, 10):
                    Pr = orth((U[:, :rank] * S[:rank]) @ Vt[:rank], ev.u)
                    row[f"M3_rank{rank}"] = m3(ev, Pr)
                    row[f"pred_energy_top{rank}"] = float(energy[:rank].sum())
                rows.append(row)
        print(CL[scperteval.CONTEXTS[t]], flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "exploratory_rank.csv", index=False)
    print(df.groupby(["context", "k"]).median(numeric_only=True).drop(columns="draw").round(3))


if __name__ == "__main__":
    main()
