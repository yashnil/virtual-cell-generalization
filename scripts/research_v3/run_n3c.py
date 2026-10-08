"""N3-C: rank structure of learnable γ⊥ (n3_protocol.md §4, procedure fixed in advance).

Estimator side: SVD of E4's centred, consensus-orthogonal test prediction; gain
fraction f(r) of the top-r reconstruction; r_eff90; participation ratio; axis
stability across draws against a permuted-anchor null. Data side: unbiased
fraction o(r) of reliable test γ⊥ energy inside the top-r subspace of the pool
perturbations' fit-part γ⊥ (independent of E4).

Reproduce: ``uv run python scripts/research_v3/run_n3c.py``
"""

from __future__ import annotations

import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from run_n3a import CONTEXTS, DATA, SEED, anchor_draw, source_view, target_parts, test_split

from virtual_cell.analysis import calibration_budget as cb

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "outputs" / "n3" / "n3c"
KS = (20, 50, 200)
DRAWS_PER_REPEAT = 4
RANKS = (1, 2, 3, 5, 10)
N_PERM = 10


def orth(x: np.ndarray, u: np.ndarray) -> np.ndarray:
    return x - np.einsum("pg,pg->p", x, u)[:, None] * u


def m3(ev: cb.EvalContext, Po: np.ndarray) -> float:
    a1, a2 = orth(ev.e1c, ev.u), orth(ev.e2c, ev.u)
    num = np.einsum("pg,pg->p", a1 - Po, a2 - Po).sum()
    return float(1 - num / np.einsum("pg,pg->p", a1, a2).sum())


def decompose(ev: cb.EvalContext, P: np.ndarray, *, cumulative: bool) -> dict:
    Po = orth(P - P.mean(axis=0), ev.u)
    U, S, Vt = np.linalg.svd(Po, full_matrices=False)
    out = {
        "M3_full": m3(ev, Po),
        "pr": float(np.sum(S**2) ** 2 / np.sum(S**4)),
        "v1": Vt[0],
        "V3": Vt[:3],
    }
    rmax = 50 if cumulative else max(RANKS)
    for r in range(1, rmax + 1):
        if r in RANKS or cumulative:
            out[f"M3_rank{r}"] = m3(ev, orth((U[:, :r] * S[:r]) @ Vt[:r], ev.u))
    return out


def run_target(t: int) -> dict:
    t0 = time.time()
    D = np.load(DATA / "delta6.npy").astype(np.float64)
    src = [c for c in range(len(CONTEXTS)) if c != t]
    sv = source_view(D, src)
    del D
    T, pool = test_split(t, sv.A.shape[0])
    Pp, Cp, _ = target_parts(t)
    est_rows, axes, null_axes = [], {}, {}
    for r in range(Pp.shape[0]):
        fit = Pp[r, 0] - Cp[r, 0]
        ev = cb.make_eval_context((Pp[r, 1] - Cp[r, 1])[T], (Pp[r, 2] - Cp[r, 2])[T], sv.A[T])
        plan = [(k, d) for k in KS for d in range(DRAWS_PER_REPEAT)] + [("ref", 0)]
        for k, d in plan:
            K = pool if k == "ref" else anchor_draw(t, k, r, d, pool)
            kk = len(K)
            P4, _ = cb.e4(sv, T, K, fit[K])
            res = decompose(ev, P4, cumulative=True)
            axes.setdefault(kk, []).append((res.pop("v1"), res.pop("V3")))
            est_rows.append({"context": CONTEXTS[t], "k": kk, "repeat": r, "draw": d, **res})
        if r < 2:  # permuted-anchor null: 10 fits per k over repeats 0-1
            for k in KS:
                for d in range(N_PERM // 2):
                    K = anchor_draw(t, k, r, d, pool)
                    pi = np.random.default_rng([SEED, 18, t, k, r * 100 + d]).permutation(k)
                    P4, _ = cb.e4(sv, T, K, fit[K][pi])
                    Po = orth(P4 - P4.mean(axis=0), ev.u)
                    _, _, Vt = np.linalg.svd(Po, full_matrices=False)
                    null_axes.setdefault(k, []).append((Vt[0], Vt[:3]))
        print(f"  {CONTEXTS[t]} repeat {r} [{time.time() - t0:.0f}s]", flush=True)

    def stability(pairs):
        cos, sub = [], []
        for i in range(len(pairs)):
            for j in range(i + 1, len(pairs)):
                cos.append(abs(float(pairs[i][0] @ pairs[j][0])))
                sub.append(float(np.sum((pairs[i][1] @ pairs[j][1].T) ** 2)) / 3)
        return np.array(cos), np.array(sub)

    stab_rows = []
    for k, pairs in axes.items():
        c, s = stability(pairs)
        row = {
            "context": CONTEXTS[t],
            "k": k,
            "n_fits": len(pairs),
            "v1_abscos_median": float(np.median(c)),
            "top3_overlap_median": float(np.median(s)),
        }
        if k in null_axes:
            cn, sn = stability(null_axes[k])
            row.update(
                {
                    "null_v1_abscos_median": float(np.median(cn)),
                    "null_v1_abscos_p975": float(np.percentile(cn, 97.5)),
                    "null_top3_overlap_median": float(np.median(sn)),
                    "null_top3_overlap_p975": float(np.percentile(sn, 97.5)),
                }
            )
        stab_rows.append(row)

    # data side: pool fit-part γ⊥ subspace, repeat 0
    fit0 = Pp[0, 0] - Cp[0, 0]
    Fp = fit0[pool] - fit0[pool].mean(axis=0)
    Ap = sv.A[pool] - sv.A[pool].mean(axis=0)
    up = Ap / np.maximum(np.linalg.norm(Ap, axis=1, keepdims=True), 1e-12)
    _, _, Vt_pool = np.linalg.svd(orth(Fp, up), full_matrices=False)
    ev0 = cb.make_eval_context((Pp[0, 1] - Cp[0, 1])[T], (Pp[0, 2] - Cp[0, 2])[T], sv.A[T])
    a1, a2 = orth(ev0.e1c, ev0.u), orth(ev0.e2c, ev0.u)
    den = np.einsum("pg,pg->p", a1, a2).sum()
    data_rows = []
    for r in (1, 3, 10, 30):
        V = Vt_pool[:r]
        num = np.einsum("pg,pg->p", (a1 @ V.T) @ V, (a2 @ V.T) @ V).sum()
        data_rows.append({"context": CONTEXTS[t], "rank": r, "o_r": float(num / den)})
    return {"est": est_rows, "stab": stab_rows, "data": data_rows}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with Pool(len(CONTEXTS)) as pool:
        res = pool.map(run_target, range(len(CONTEXTS)))
    for name in ("est", "stab", "data"):
        pd.DataFrame([r for rr in res for r in rr[name]]).to_csv(
            OUT / f"n3c_{name}.csv", index=False
        )
    print(f"done [{time.time() - t0:.0f}s]")


if __name__ == "__main__":
    main()
