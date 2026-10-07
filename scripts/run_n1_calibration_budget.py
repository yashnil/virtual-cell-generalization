"""N1: calibration-budget learning curve (reports/n1_n4_protocol.md §2).

Per outer LOCO fold: fixed 30 % test set, k anchors from the remaining pool,
200 draws per k (5 cell-split repeats x 40), estimators E0–E4, metrics M0–M4.
Also runs the F3a anchor-permutation null and the F3b shared-control diagnostic.

Reproduce: ``uv run python scripts/run_n1_calibration_budget.py``
(requires ``scripts/build_n1_n4_splits.py`` first).
"""

from __future__ import annotations

import json
import time
from multiprocessing import Pool
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
KS = (0, 1, 2, 5, 10, 20, 50, 100, 200)
DRAWS_PER_REPEAT = 40
N_TEST = 379
N_BOOT = 2000
CL = {c: scperteval.dataset(c).cell_line for c in scperteval.CONTEXTS}


def target_parts(t: int):
    pert = np.load(SPLITS / "pert_part_means.npy", mmap_mode="r")
    ctrl = np.load(SPLITS / "ctrl_part_means.npy")
    n_ctrl = np.load(SPLITS / "n_ctrl_cells.npy")
    P = np.asarray(pert[:, :, t], dtype=np.float64)  # (R, 3, P, G)
    C = ctrl[:, :, t].astype(np.float64)  # (R, 3, G)
    w = n_ctrl[:, :, t].astype(np.float64)  # (R, 3)
    C_all = (C * w[:, :, None]).sum(axis=1) / w.sum(axis=1)[:, None]  # (R, G)
    return P, C, C_all


def run_fold(t: int) -> dict:
    t0 = time.time()
    D = np.load(CANON / "delta_tensor.npy").astype(np.float64)
    fold = next(f for f in loco.make_folds(scperteval.CONTEXTS) if f.target_index == t)
    src = list(fold.source_indices)
    sv = cb.make_source_view(D, src, pr.fit_scale(D, src))
    del D
    n_p = sv.A.shape[0]
    perm = np.random.default_rng([SEED, 1, t]).permutation(n_p)
    T, pool = np.sort(perm[:N_TEST]), np.sort(perm[N_TEST:])
    Pparts, Cparts, C_all = target_parts(t)
    R = Pparts.shape[0]

    rows, lam_rows, f3a_rows, f3b_rows = [], [], [], []
    acc: dict = {}  # (k, est, metric) -> [num_sum, den_sum, n]

    def record(k, r, d, preds, lams, ev, tag_rows, extra=None):
        for est, P_ in preds.items():
            terms = cb.metric_terms(ev, P_)
            row = {
                "context": CL[scperteval.CONTEXTS[t]],
                "k": k,
                "repeat": r,
                "draw": d,
                "estimator": est,
            }
            row.update({m: cb.ratio(*terms[m]) for m in cb.METRIC_TERMS})
            row.update(cb.correlation_terms(ev, P_))
            if extra:
                row.update(extra)
            tag_rows.append(row)
            if tag_rows is rows:
                for m in cb.METRIC_TERMS:
                    key = (k, est, m)
                    if key not in acc:
                        acc[key] = [np.zeros(len(T)), np.zeros(len(T)), 0]
                    acc[key][0] += terms[m][0]
                    acc[key][1] += terms[m][1]
                    acc[key][2] += 1
        for est, lam in lams.items():
            lam_rows.append(
                {
                    "context": CL[scperteval.CONTEXTS[t]],
                    "k": k,
                    "repeat": r,
                    "draw": d,
                    "estimator": est,
                    "lambda": lam,
                }
            )

    for r in range(R):
        fit = Pparts[r, 0] - Cparts[r, 0]
        e1 = Pparts[r, 1] - Cparts[r, 1]
        e2 = Pparts[r, 2] - Cparts[r, 2]
        ev = cb.make_eval_context(e1[T], e2[T], sv.A[T])
        # shared-control versions (F3b only)
        fit_s = Pparts[r, 0] - C_all[r]
        ev_s = cb.make_eval_context(Pparts[r, 1][T] - C_all[r], Pparts[r, 2][T] - C_all[r], sv.A[T])

        for k in KS + ("ref",):
            if k == 0:
                preds, lams = cb.predict_all(sv, T, None, None)
                record(0, r, 0, preds, lams, ev, rows)
                continue
            if k == "ref":
                preds, lams = cb.predict_all(sv, T, pool, fit[pool])
                record(len(pool), r, 0, preds, lams, ev, rows)
                continue
            for d in range(DRAWS_PER_REPEAT):
                K = np.sort(
                    np.random.default_rng([SEED, 2, t, k, r, d]).choice(pool, k, replace=False)
                )
                preds, lams = cb.predict_all(sv, T, K, fit[K])
                record(k, r, d, preds, lams, ev, rows)
                if d < 10 and k in (20, 100):  # F3a anchor-permutation null
                    pi = np.random.default_rng([SEED, 6, t, k, r, d]).permutation(k)
                    pp, _ = cb.predict_all(sv, T, K, fit[K][pi])
                    record(
                        k,
                        r,
                        d,
                        {e: pp[e] for e in ("E3", "E4")},
                        {},
                        ev,
                        f3a_rows,
                        {"variant": "permuted"},
                    )
                if d < 10 and k in (5, 20, 100):  # F3b shared-control diagnostic
                    ps, _ = cb.predict_all(sv, T, K, fit_s[K])
                    record(
                        k,
                        r,
                        d,
                        {e: ps[e] for e in ("E1", "E3", "E4")},
                        {},
                        ev_s,
                        f3b_rows,
                        {"variant": "shared_control"},
                    )
                    record(
                        k,
                        r,
                        d,
                        {e: preds[e] for e in ("E1", "E3", "E4")},
                        {},
                        ev,
                        f3b_rows,
                        {"variant": "split_control"},
                    )
        print(f"  {CL[scperteval.CONTEXTS[t]]} repeat {r} [{time.time() - t0:.0f}s]", flush=True)

    boot = []
    for (k, est, m), (num, den, n) in acc.items():
        lo, hi = cb.bootstrap_ratio(
            num / n, den / n, n_boot=N_BOOT, rng=np.random.default_rng([SEED, 3, t, int(k)])
        )
        boot.append(
            {
                "context": CL[scperteval.CONTEXTS[t]],
                "k": k,
                "estimator": est,
                "metric": m,
                "point": cb.ratio(num, den),
                "boot_lo": lo,
                "boot_hi": hi,
            }
        )
    return {
        "draws": rows,
        "lambdas": lam_rows,
        "f3a": f3a_rows,
        "f3b": f3b_rows,
        "boot": boot,
        "test_set": T.tolist(),
        "pool_size": int(len(pool)),
        "context": CL[scperteval.CONTEXTS[t]],
        "scale": sv.scale,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with Pool(4) as pool:
        results = pool.map(run_fold, range(len(scperteval.CONTEXTS)))
    for name in ("draws", "lambdas", "f3a", "f3b", "boot"):
        pd.DataFrame([row for res in results for row in res[name]]).to_csv(
            OUT / f"n1_{name}.csv", index=False
        )
    meta = {
        res["context"]: {
            "test_set": res["test_set"],
            "pool_size": res["pool_size"],
            "frozen_scale": res["scale"],
        }
        for res in results
    }
    (OUT / "n1_meta.json").write_text(json.dumps(meta))
    print(f"done [{time.time() - t0:.0f}s]")


if __name__ == "__main__":
    main()
