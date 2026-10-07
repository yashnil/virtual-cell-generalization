"""N3-A: six-context replication of the N1 calibration-budget curve (n3_protocol.md §1).

N1 protocol §2 with 6 LOCO folds (5 sources each), test set 319, pool 743,
200 draws per k, estimators E0–E4, metrics M0–M4, F3a/F3b diagnostics.

Reproduce: ``uv run python scripts/research_v3/run_n3a.py`` (after build_n3_six_context.py).
"""

from __future__ import annotations

import json
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

from virtual_cell.analysis import calibration_budget as cb
from virtual_cell.modelling import pathway_residual as pr

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "outputs" / "n3" / "data"
OUT = REPO / "outputs" / "n3" / "n3a"
SEED = 20261006
CONTEXTS = ("K562", "RPE1", "HepG2", "Jurkat", "HCT116", "HEK293T")
KS = (0, 1, 2, 5, 10, 20, 50, 100, 200)
DRAWS_PER_REPEAT = 40
N_BOOT = 2000


def test_split(t: int, n_p: int) -> tuple[np.ndarray, np.ndarray]:
    n_test = round(0.3 * n_p)
    perm = np.random.default_rng([SEED, 11, t]).permutation(n_p)
    return np.sort(perm[:n_test]), np.sort(perm[n_test:])


def anchor_draw(t: int, k: int, r: int, d: int, pool: np.ndarray) -> np.ndarray:
    return np.sort(np.random.default_rng([SEED, 12, t, k, r, d]).choice(pool, k, replace=False))


def target_parts(t: int):
    pert = np.load(DATA / "pert_part_means.npy", mmap_mode="r")
    ctrl = np.load(DATA / "ctrl_part_means.npy")
    n_ctrl = np.load(DATA / "n_ctrl_cells.npy")
    P = np.asarray(pert[:, :, t], dtype=np.float64)
    C = ctrl[:, :, t].astype(np.float64)
    w = n_ctrl[:, :, t].astype(np.float64)
    C_all = (C * w[:, :, None]).sum(axis=1) / w.sum(axis=1)[:, None]
    return P, C, C_all


def source_view(D: np.ndarray, sources) -> cb.SourceView:
    return cb.make_source_view(D, list(sources), pr.fit_scale(D, list(sources)))


def run_fold(t: int) -> dict:
    t0 = time.time()
    D = np.load(DATA / "delta6.npy").astype(np.float64)
    src = [c for c in range(len(CONTEXTS)) if c != t]
    sv = source_view(D, src)
    del D
    T, pool = test_split(t, sv.A.shape[0])
    Pp, Cp, C_all = target_parts(t)
    cl = CONTEXTS[t]
    rows, lam_rows, f3a_rows, f3b_rows = [], [], [], []
    acc: dict = {}

    def record(k, r, d, preds, lams, ev, sink, extra=None):
        for est, P_ in preds.items():
            terms = cb.metric_terms(ev, P_)
            row = {"context": cl, "k": k, "repeat": r, "draw": d, "estimator": est}
            row.update({m: cb.ratio(*terms[m]) for m in cb.METRIC_TERMS})
            row.update(cb.correlation_terms(ev, P_))
            if extra:
                row.update(extra)
            sink.append(row)
            if sink is rows:
                for m in cb.METRIC_TERMS:
                    a = acc.setdefault((k, est, m), [np.zeros(len(T)), np.zeros(len(T)), 0])
                    a[0] += terms[m][0]
                    a[1] += terms[m][1]
                    a[2] += 1
        for est, lam in lams.items():
            lam_rows.append(
                {"context": cl, "k": k, "repeat": r, "draw": d, "estimator": est, "lambda": lam}
            )

    for r in range(Pp.shape[0]):
        fit = Pp[r, 0] - Cp[r, 0]
        ev = cb.make_eval_context((Pp[r, 1] - Cp[r, 1])[T], (Pp[r, 2] - Cp[r, 2])[T], sv.A[T])
        fit_s = Pp[r, 0] - C_all[r]
        ev_s = cb.make_eval_context(Pp[r, 1][T] - C_all[r], Pp[r, 2][T] - C_all[r], sv.A[T])
        for k in (*KS, "ref"):
            if k == 0:
                record(0, r, 0, *cb.predict_all(sv, T, None, None), ev, rows)
                continue
            if k == "ref":
                record(len(pool), r, 0, *cb.predict_all(sv, T, pool, fit[pool]), ev, rows)
                continue
            for d in range(DRAWS_PER_REPEAT):
                K = anchor_draw(t, k, r, d, pool)
                preds, lams = cb.predict_all(sv, T, K, fit[K])
                record(k, r, d, preds, lams, ev, rows)
                if d < 10 and k in (20, 100):
                    pi = np.random.default_rng([SEED, 16, t, k, r, d]).permutation(k)
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
                if d < 10 and k in (5, 20, 100):
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
        print(f"  {cl} repeat {r} [{time.time() - t0:.0f}s]", flush=True)

    boot = []
    for (k, est, m), (num, den, n) in acc.items():
        lo, hi = cb.bootstrap_ratio(
            num / n, den / n, n_boot=N_BOOT, rng=np.random.default_rng([SEED, 13, t, int(k)])
        )
        boot.append(
            {
                "context": cl,
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
        "meta": {
            "context": cl,
            "test_set": T.tolist(),
            "pool_size": int(len(pool)),
            "frozen_scale": sv.scale,
        },
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with Pool(len(CONTEXTS)) as pool:
        results = pool.map(run_fold, range(len(CONTEXTS)))
    for name in ("draws", "lambdas", "f3a", "f3b", "boot"):
        pd.DataFrame([row for res in results for row in res[name]]).to_csv(
            OUT / f"n3a_{name}.csv", index=False
        )
    (OUT / "n3a_meta.json").write_text(
        json.dumps({r["meta"]["context"]: r["meta"] for r in results})
    )
    print(f"done [{time.time() - t0:.0f}s]")


if __name__ == "__main__":
    main()
