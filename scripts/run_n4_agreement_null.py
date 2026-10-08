"""N4: does source agreement predict held-out transfer beyond simple statistics?

Preregistered in reports/n1_n4_protocol.md §3. Per held-out context: T1 marginal
comparison against b1–b6, T2 partial Spearman, T3 exchangeable random-effects null.

Declared deviation (recorded before any N4 output existed): target and source
half energies are computed **per repeat** and then averaged, rather than averaging
the halves over repeats first. Averaged halves share cells across repeats, which
correlates their noise and biases ⟨h1, h2⟩.

Reproduce: ``uv run python scripts/run_n4_agreement_null.py``
"""

from __future__ import annotations

import json
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

from virtual_cell.analysis import agreement_null as an
from virtual_cell.analysis import foundations, loco
from virtual_cell.data import scperteval
from virtual_cell.modelling import pathway_residual as pr

REPO = Path(__file__).resolve().parents[1]
CANON = REPO / "data" / "processed" / "four_context_v1"
DESIGN = REPO / "data" / "splits" / "four_context_v1"
SPLITS = REPO / "outputs" / "n1_n4" / "splits"
OUT = REPO / "outputs" / "n1_n4" / "n4"
STABILITY = REPO / "outputs" / "transferability_v1" / "stability_rule.json"
SEED = 20261006
N_BOOT = 2000
N_SIM = 200
CL = {c: scperteval.dataset(c).cell_line for c in scperteval.CONTEXTS}
BASELINES = (
    "b1_magnitude",
    "b2_source_reliability",
    "b3_source_reliable_energy",
    "b4_min_source_cells",
    "b5_target_gene_basal",
    "b6_noise_agreement_ceiling",
)
REQUIRED = (
    "b1_magnitude",
    "b2_source_reliability",
    "b3_source_reliable_energy",
    "b6_noise_agreement_ceiling",
)


def halves(ctx: int, centre: bool = False):
    pert = np.load(SPLITS / "pert_part_means.npy", mmap_mode="r")
    ctrl = np.load(SPLITS / "ctrl_part_means.npy")
    P = np.asarray(pert[:, :, ctx], dtype=np.float64)
    C = ctrl[:, :, ctx].astype(np.float64)
    h1 = [P[r, 0] - C[r, 0] for r in range(P.shape[0])]
    h2 = [((P[r, 1] - C[r, 1]) + (P[r, 2] - C[r, 2])) / 2 for r in range(P.shape[0])]
    if centre:
        h1 = [h - h.mean(axis=0) for h in h1]
        h2 = [h - h.mean(axis=0) for h in h2]
    return h1, h2


def run_fold(t: int) -> dict:
    t0 = time.time()
    min_signal = float(json.loads(STABILITY.read_text())["min_signal_energy"])
    D = np.load(CANON / "delta_tensor.npy").astype(np.float64)
    control = np.load(CANON / "control_means.npy").astype(np.float64)
    counts = np.load(CANON / "cell_counts.npy").astype(float)
    perts = (DESIGN / "shared_perturbations.txt").read_text().split()
    genes = (DESIGN / "shared_genes.txt").read_text().split()
    fold = next(f for f in loco.make_folds(scperteval.CONTEXTS) if f.target_index == t)
    src = list(fold.source_indices)
    cl = CL[scperteval.CONTEXTS[t]]
    scale = pr.fit_scale(D, src)
    B = pr.baseline(D, src, scale)
    a = an.agreement(D[src])

    # ---- quality (target halves, evaluation only) ---------------------------
    th1, th2 = halves(t)
    q, sig = an.quality_d(th1, th2, B, min_signal=min_signal)
    full_t = np.mean([(x + y) / 2 for x, y in zip(th1, th2, strict=True)], axis=0)
    r_sec = an.rowwise_pearson(B, full_t)

    # ---- simple statistics (source-only or target basal) ---------------------
    rel, ren, sig2_src = [], [], []
    for s in src:
        h1, h2 = halves(s)
        rel.append(np.mean([an.rowwise_pearson(x, y) for x, y in zip(h1, h2, strict=True)], axis=0))
        ren.append(
            np.mean([np.einsum("pg,pg->p", x, y) for x, y in zip(h1, h2, strict=True)], axis=0)
        )
        sig2_src.append(an.noise_variance(h1, h2))
    rel = np.vstack(rel)
    rho_full = np.clip(2 * rel / (1 + rel), 0, 1)
    pairs = [(i, j) for i in range(len(src)) for j in range(i + 1, len(src))]
    gi = {g: i for i, g in enumerate(genes)}
    stats_tab = {
        "agreement": a,
        "b1_magnitude": np.linalg.norm(D[src].mean(axis=0), axis=1),
        "b2_source_reliability": rel.mean(axis=0),
        "b3_source_reliable_energy": np.vstack(ren).mean(axis=0),
        "b4_min_source_cells": counts[src].min(axis=0),
        "b5_target_gene_basal": np.array([control[t, gi[p]] if p in gi else np.nan for p in perts]),
        "b6_noise_agreement_ceiling": np.mean(
            [np.sqrt(rho_full[i] * rho_full[j]) for i, j in pairs], axis=0
        ),
    }

    # ---- T1 -----------------------------------------------------------------
    t1 = []
    for name in ("agreement",) + BASELINES:
        row = {
            "context": cl,
            "statistic": name,
            "spearman_q": an.spearman(stats_tab[name], q),
            "spearman_r_secondary": an.spearman(stats_tab[name], r_sec),
        }
        if name != "agreement":
            pt, lo, hi = an.paired_bootstrap_spearman_diff(
                a,
                stats_tab[name],
                q,
                n_boot=N_BOOT,
                rng=np.random.default_rng([SEED, 4, t, BASELINES.index(name)]),
            )
            row.update({"diff_vs_agreement": pt, "diff_lo": lo, "diff_hi": hi})
        t1.append(row)

    # ---- T2 -----------------------------------------------------------------
    cov = np.column_stack([stats_tab[b] for b in BASELINES[:5]])
    ok = np.isfinite(q) & np.all(np.isfinite(cov), axis=1) & np.isfinite(a)
    part = foundations.partial_spearman(a[ok], q[ok], cov[ok])
    rng = np.random.default_rng([SEED, 4, t, 99])
    idx_ok = np.flatnonzero(ok)
    boots = []
    for _ in range(N_BOOT):
        i = rng.choice(idx_ok, len(idx_ok))
        boots.append(foundations.partial_spearman(a[i], q[i], cov[i]))
    t2 = {
        "context": cl,
        "partial_spearman": part,
        "lo": float(np.nanpercentile(boots, 2.5)),
        "hi": float(np.nanpercentile(boots, 97.5)),
        "n": int(ok.sum()),
    }
    print(f"  {cl}: T1/T2 done [{time.time() - t0:.0f}s]", flush=True)

    # ---- T3: exchangeable random-effects null --------------------------------
    S_c = D[src] - D[src].mean(axis=1, keepdims=True)
    sig2_src = np.vstack(sig2_src)
    m, tau2 = an.fit_exchangeable(S_c, sig2_src)
    ch1, ch2 = halves(t, centre=True)
    sig2_tgt = an.noise_variance(ch1, ch2)
    a_c = an.agreement(S_c)
    q_c, _ = an.quality_d(ch1, ch2, scale * S_c.mean(axis=0), min_signal=min_signal)
    observed = an.spearman(a_c, q_c)
    null = [
        an.simulate_null(
            m,
            tau2,
            sig2_src,
            sig2_tgt,
            scale=scale,
            min_signal=min_signal,
            rng=np.random.default_rng([SEED, 5, t, i]),
        )
        for i in range(N_SIM)
    ]
    lo, hi = np.percentile(null, [2.5, 97.5])
    cls = "below null" if observed < lo else ("above null" if observed > hi else "within null")
    t3 = {
        "context": cl,
        "observed_centred": observed,
        "null_median": float(np.median(null)),
        "null_lo": float(lo),
        "null_hi": float(hi),
        "class": cls,
        "tau2_median": float(np.median(tau2)),
        "sigma2_src_median": float(np.median(sig2_src)),
        "sigma2_tgt_median": float(np.median(sig2_tgt)),
        "n_stable_observed": int(np.isfinite(q_c).sum()),
    }
    print(f"  {cl}: T3 done [{time.time() - t0:.0f}s]", flush=True)

    per_pert = pd.DataFrame(
        {
            "context": cl,
            "perturbation": perts,
            "q": q,
            "signal": sig,
            "r_secondary": r_sec,
            **stats_tab,
        }
    )
    return {"t1": t1, "t2": t2, "t3": t3, "null": null, "per_pert": per_pert}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with Pool(4) as pool:
        res = pool.map(run_fold, range(len(scperteval.CONTEXTS)))
    t1 = pd.DataFrame([row for r in res for row in r["t1"]])
    t2 = pd.DataFrame([r["t2"] for r in res])
    t3 = pd.DataFrame([r["t3"] for r in res])
    t1.to_csv(OUT / "n4_t1.csv", index=False)
    t2.to_csv(OUT / "n4_t2.csv", index=False)
    t3.to_csv(OUT / "n4_t3.csv", index=False)
    pd.DataFrame({r["t3"]["context"]: r["null"] for r in res}).to_csv(
        OUT / "n4_t3_null.csv", index=False
    )
    pd.concat([r["per_pert"] for r in res]).to_csv(OUT / "n4_per_perturbation.csv", index=False)

    verdict = {}
    for r in res:
        cl = r["t2"]["context"]
        sub = t1[(t1.context == cl) & t1.statistic.isin(REQUIRED)]
        verdict[cl] = {
            "partial_ok": bool(r["t2"]["lo"] > 0),
            "beats_required": bool((sub.diff_lo > 0).all()),
            "t3_class": r["t3"]["class"],
        }
    n_strong = sum(v["partial_ok"] and v["beats_required"] for v in verdict.values())
    n_weak = sum(v["partial_ok"] for v in verdict.values())
    decision = {
        "per_context": verdict,
        "n_pass_strong": n_strong,
        "n_partial_ok": n_weak,
        "N4": "PASS" if n_strong >= 3 else ("WEAK" if n_weak >= 3 else "FAIL"),
    }
    (OUT / "n4_decision.json").write_text(json.dumps(decision, indent=2))
    print(json.dumps(decision, indent=2))
    print(f"done [{time.time() - t0:.0f}s]")


if __name__ == "__main__":
    main()
