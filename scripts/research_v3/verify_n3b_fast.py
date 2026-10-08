"""Check that n3b_fast reproduces calibration_budget E0s/E2/E4 and metrics on real N3 data."""

from __future__ import annotations

import json
import sys

import n3b_fast as fast
import numpy as np
from run_n3a import DATA, anchor_draw, source_view, target_parts, test_split

from virtual_cell.analysis import calibration_budget as cb


def main() -> None:
    D = np.load(DATA / "delta6.npy").astype(np.float64)
    worst = {}
    for t, subset in ((4, (0, 1, 2, 3, 5)), (5, (1, 4)), (0, (2,))):
        src_all = [c for c in range(6) if c != t]
        sv_all = source_view(D, src_all)
        sv = source_view(D, list(subset))
        lv = fast.light_view(sv)
        T, pool = test_split(t, D.shape[1])
        Pp, Cp, _ = target_parts(t)
        fit = Pp[0, 0] - Cp[0, 0]
        e1, e2 = (Pp[0, 1] - Cp[0, 1])[T], (Pp[0, 2] - Cp[0, 2])[T]
        ev_ref = cb.make_eval_context(e1, e2, sv_all.A[T])
        ev_own = cb.make_eval_context(e1, e2, sv.A[T])
        fe = fast.fast_eval(e1, e2, {"ref": sv_all.A[T], "own": sv.A[T]})
        for K in (anchor_draw(t, 5, 0, 0, pool), anchor_draw(t, 50, 0, 1, pool), pool):
            Y = fit[K]
            P_ref, lam_ref = cb.e4(sv, T, K, Y)
            P_new, lam_new = fast.e4(lv, T, K, Y)
            assert lam_ref == lam_new, (lam_ref, lam_new)
            pairs = {
                "E4_pred": np.max(np.abs(P_ref - P_new)) / np.max(np.abs(P_ref)),
                "E2_pred": np.max(np.abs(cb.e2(sv, T, K, Y) - fast.e2(lv, T, K, Y))),
            }
            for name, Pq in (("E4", P_ref), ("E2", cb.e2(sv, T, K, Y)), ("E0s", cb.e0s(sv, T))):
                mt = fast.fast_metrics(fe, Pq)
                ref = cb.metric_terms(ev_ref, Pq)
                own = cb.metric_terms(ev_own, Pq)
                pairs[f"{name}_M0"] = abs(mt["M0"] - cb.ratio(*ref["M0"]))
                pairs[f"{name}_M1"] = abs(mt["M1"] - cb.ratio(*ref["M1"]))
                pairs[f"{name}_M3ref"] = abs(mt["M3_ref"] - cb.ratio(*ref["M3"]))
                pairs[f"{name}_M3own"] = abs(mt["M3_own"] - cb.ratio(*own["M3"]))
            for key, v in pairs.items():
                worst[key] = max(worst.get(key, 0.0), float(v))
    print(json.dumps(worst, indent=2))
    ok = all(v < 1e-8 for v in worst.values())
    print("VERIFIED" if ok else "MISMATCH")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
