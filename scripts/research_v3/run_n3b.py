"""N3-B: source-count ablation R_gamma(k, m) (n3_protocol.md §3).

For each target, every non-empty subset of its 5 sources (31) is evaluated on the
same test set and the same anchor draws (the first 8 per repeat of N3-A). E4 is
the calibration estimator; E0s (k = 0) and E2 are references. γ⊥ is scored
against the fixed 5-source consensus (M3_ref, primary) and against the subset's
own consensus (M3_own, secondary). No subset is selected.

Reproduce: ``uv run python scripts/research_v3/run_n3b.py``
"""

from __future__ import annotations

import itertools
import time
from multiprocessing import Pool
from pathlib import Path

import n3b_fast as fast
import numpy as np
import pandas as pd
from run_n3a import CONTEXTS, DATA, anchor_draw, source_view, target_parts, test_split

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "outputs" / "n3" / "n3b"
KS = (5, 10, 20, 50, 100, 200)
DRAWS_PER_REPEAT = 8


def run_target(t: int) -> list[dict]:
    t0 = time.time()
    D = np.load(DATA / "delta6.npy").astype(np.float64)
    src_all = [c for c in range(len(CONTEXTS)) if c != t]
    A_all = source_view(D, src_all).A
    T, pool = test_split(t, D.shape[1])
    subsets = [s for m in range(1, 6) for s in itertools.combinations(src_all, m)]
    views = {s: fast.light_view(source_view(D, s)) for s in subsets}
    del D
    Pp, Cp, _ = target_parts(t)
    rows = []
    for r in range(Pp.shape[0]):
        fit = Pp[r, 0] - Cp[r, 0]
        e1, e2 = (Pp[r, 1] - Cp[r, 1])[T], (Pp[r, 2] - Cp[r, 2])[T]
        fe = {s: fast.fast_eval(e1, e2, {"ref": A_all[T], "own": views[s].A[T]}) for s in subsets}
        plan = [(0, None)] + [(k, d) for k in KS for d in range(DRAWS_PER_REPEAT)] + [("ref", 0)]
        for k, d in plan:
            if k == 0:
                K = None
            elif k == "ref":
                K = pool
            else:
                K = anchor_draw(t, k, r, d, pool)
            kk = 0 if k == 0 else len(K)
            for s in subsets:
                lv = views[s]
                base = {
                    "context": CONTEXTS[t],
                    "k": kk,
                    "repeat": r,
                    "draw": d or 0,
                    "subset": "+".join(CONTEXTS[c] for c in s),
                    "m": len(s),
                }
                if K is None:
                    mt = fast.fast_metrics(fe[s], fast.e0s(lv, T))
                    rows.append(
                        {
                            **base,
                            "estimator": "E0s",
                            "M0": mt["M0"],
                            "M1": mt["M1"],
                            "M3_ref": mt["M3_ref"],
                            "M3_own": mt["M3_own"],
                        }
                    )
                    continue
                Y = fit[K]
                P4, lam = fast.e4(lv, T, K, Y)
                for est, P_ in (("E2", fast.e2(lv, T, K, Y)), ("E4", P4)):
                    mt = fast.fast_metrics(fe[s], P_)
                    row = {
                        **base,
                        "estimator": est,
                        "M0": mt["M0"],
                        "M1": mt["M1"],
                        "M3_ref": mt["M3_ref"],
                    }
                    if est == "E4":
                        row["M3_own"] = mt["M3_own"]
                        row["lambda"] = lam
                    rows.append(row)
        print(f"  {CONTEXTS[t]} repeat {r} [{time.time() - t0:.0f}s]", flush=True)
    return rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    todo = [t for t in range(len(CONTEXTS)) if not (OUT / f"n3b_draws_{CONTEXTS[t]}.csv").exists()]
    with Pool(len(todo) or 1) as pool:
        for t, rows in zip(todo, pool.imap(run_target, todo), strict=True):
            pd.DataFrame(rows).to_csv(OUT / f"n3b_draws_{CONTEXTS[t]}.csv", index=False)
    pd.concat([pd.read_csv(OUT / f"n3b_draws_{c}.csv") for c in CONTEXTS]).to_csv(
        OUT / "n3b_draws.csv", index=False
    )
    print(f"done [{time.time() - t0:.0f}s]")


if __name__ == "__main__":
    main()
