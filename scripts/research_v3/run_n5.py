"""N5: source compatibility endpoints (reports/n5_protocol.md §2–4).

Order of operations enforces the protocol: source-only quantities and the
matched comparator are computed and written **before** the target is loaded.

Reproduce: ``uv run python scripts/research_v3/run_n5.py`` (after build_n5_gwps.py).
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import n3b_fast as fast
import numpy as np
import pandas as pd

from virtual_cell.analysis import calibration_budget as cb
from virtual_cell.analysis import source_compatibility as sc

REPO = Path(__file__).resolve().parents[2]
AXES = REPO / "data" / "splits" / "n5_k562"
N3AX = REPO / "data" / "splits" / "six_context_n3"
N3 = REPO / "outputs" / "n3" / "data"
GW = REPO / "outputs" / "n5" / "gwps"
OUT = REPO / "outputs" / "n5"
SEED = 20261006
N_BOOT = 2000
N3_INDEX = {"K562_essential": 0, "RPE1": 1, "HepG2": 2, "Jurkat": 3, "HCT116": 4, "HEK293T": 5}
SOURCES = ["K562_GWPS", "RPE1", "HepG2", "Jurkat", "HCT116", "HEK293T"]


def halves_from(pert, ctrl):
    """F half and (E1+E2)/2 half per repeat."""
    a = [pert[r, 0] - ctrl[r, 0] for r in range(pert.shape[0])]
    b = [((pert[r, 1] - ctrl[r, 1]) + (pert[r, 2] - ctrl[r, 2])) / 2 for r in range(pert.shape[0])]
    return a, b


def load_n3_context(c: int, pi, gi):
    D = np.load(N3 / "delta6.npy", mmap_mode="r")
    pert = np.load(N3 / "pert_part_means.npy", mmap_mode="r")
    ctrl = np.load(N3 / "ctrl_part_means.npy")
    full = np.asarray(D[c], dtype=np.float64)[np.ix_(pi, gi)]
    pm = np.stack(
        [
            np.stack(
                [np.asarray(pert[r, j, c], dtype=np.float64)[np.ix_(pi, gi)] for j in range(3)]
            )
            for r in range(pert.shape[0])
        ]
    )
    cm = ctrl[:, :, c][..., gi].astype(np.float64)
    a, b = halves_from(pm, cm)
    cells = np.load(N3 / "cell_counts6.npy")[c][pi]
    return full, a, b, cells


def load_gwps(depth: bool):
    tag = "_depth" if depth else ""
    full = np.load(GW / f"delta{tag}.npy").astype(np.float64)
    pm = np.load(GW / f"pert_part_means{tag}.npy").astype(np.float64)
    cm = np.load(GW / f"ctrl_part_means{tag}.npy").astype(np.float64)
    a, b = halves_from(pm, cm)
    return full, a, b, np.load(GW / f"cell_counts{tag}.npy")


def main() -> None:
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    perts = (AXES / "shared_perturbations.txt").read_text().split()
    genes = (AXES / "shared_genes.txt").read_text().split()
    p3 = (N3AX / "shared_perturbations.txt").read_text().split()
    g3 = (N3AX / "shared_genes.txt").read_text().split()
    pi = np.array([p3.index(p) for p in perts])
    gi = np.array([g3.index(g) for g in genes])

    # ---- 1. source-only quantities (no target read) -------------------------
    src = {"K562_GWPS": load_gwps(False), "K562_GWPS_depth": load_gwps(True)}
    for s in SOURCES[1:]:
        src[s] = load_n3_context(N3_INDEX[s], pi, gi)
    rel = {s: sc.split_half_reliability(v[1], v[2]) for s, v in src.items()}
    mag = {s: np.linalg.norm(sc.centre(v[0]), axis=1) for s, v in src.items()}
    cells = {s: np.asarray(v[3], dtype=float) for s, v in src.items()}
    comp = sc.matched_comparator(
        {s: rel[s] for s in SOURCES}, {s: cells[s] for s in SOURCES}, reference="K562_GWPS"
    )
    (OUT / "matched_comparator.json").write_text(json.dumps(comp, indent=2))
    pd.DataFrame({s: rel[s] for s in src}, index=perts).to_csv(OUT / "source_reliability.csv")
    print(
        f"source-only stage done; matched comparator = {comp['matched']} [{time.time() - t0:.0f}s]",
        flush=True,
    )

    # ---- 2. target ----------------------------------------------------------
    T_full, T_a, T_b, _ = load_n3_context(N3_INDEX["K562_essential"], pi, gi)
    per = {}
    for s, (S_full, S_a, S_b, _) in src.items():
        lt = sc.latent_cosine_terms(S_full, T_full, S_a, S_b, T_a, T_b)
        rc = sc.pooled_raw_cosine(S_full, T_full)
        per[s] = {
            "num": lt["num"],
            "den_S": lt["den_S"],
            "den_T": lt["den_T"],
            "raw_num": rc["num"],
            "raw_ss": rc["ss"],
            "raw_tt": rc["tt"],
            "r": sc.rowwise_pearson(sc.centre(S_full), sc.centre(T_full)),
            "dir_acc": sc.directional_accuracy(S_full, T_a, T_b),
            "rel": rel[s],
            "log_mag": np.log(np.maximum(mag[s], 1e-12)),
        }
    terms = pd.concat(
        {s: pd.DataFrame(v, index=perts) for s, v in per.items()}, names=["source", "perturbation"]
    )
    terms.to_csv(OUT / "per_perturbation_terms.csv")

    # paired perturbation bootstrap of C_S
    rng = np.random.default_rng([SEED, 9])
    idx = rng.integers(0, len(perts), size=(N_BOOT, len(perts)))
    boot = {}
    for s, v in per.items():
        num, dS, dT = v["num"][idx].sum(1), v["den_S"][idx].sum(1), v["den_T"][idx].sum(1)
        boot[s] = num / np.sqrt(dS * dT)
    pd.DataFrame(boot).to_csv(OUT / "bootstrap_C.csv", index=False)
    point = {s: sc.pooled_latent_cosine(v["num"], v["den_S"], v["den_T"]) for s, v in per.items()}
    raw = {
        s: float(v["raw_num"].sum() / np.sqrt(v["raw_ss"].sum() * v["raw_tt"].sum()))
        for s, v in per.items()
    }
    summary = {
        "C": point,
        "R2_raw_cosine": raw,
        "R1_median_r": {s: float(np.nanmedian(v["r"])) for s, v in per.items()},
        "R3_dir_acc": {s: float(np.nanmedian(v["dir_acc"])) for s, v in per.items()},
        "median_source_reliability": {s: float(np.nanmedian(rel[s])) for s in src},
        "median_source_cells": {s: float(np.median(cells[s])) for s in src},
    }

    # Adj-2: fixed-effect regression (two versions: full GWPS, depth-matched GWPS)
    for tag, gw in (("full", "K562_GWPS"), ("depth", "K562_GWPS_depth")):
        names = [gw] + SOURCES[1:]
        y = np.concatenate([per[s]["r"] for s in names])
        pert_id = np.tile(np.arange(len(perts)), len(names))
        source = np.repeat(np.array(["K562_GWPS" if s == gw else s for s in names]), len(perts))
        cov = np.column_stack(
            [
                np.concatenate([per[s]["rel"] for s in names]),
                np.concatenate([per[s]["rel"] for s in names]) ** 2,
                np.concatenate([per[s]["log_mag"] for s in names]),
            ]
        )
        est = sc.fe_regression(y, pert_id, source, cov, SOURCES, "RPE1")
        ci = sc.cluster_bootstrap_fe(
            y,
            pert_id,
            source,
            cov,
            SOURCES,
            "RPE1",
            n_boot=N_BOOT,
            rng=np.random.default_rng([SEED, 10]),
        )
        summary[f"Adj2_{tag}"] = {
            s: {"beta_minus_RPE1": est[s], "lo": ci[s][0], "hi": ci[s][1]} for s in est
        }

    # Adj-3: reliability-matched perturbations, GWPS vs RPE1
    m = np.abs(rel["K562_GWPS"] - rel["RPE1"]) <= 0.05
    d = (per["K562_GWPS"]["r"] - per["RPE1"]["r"])[m]
    d = d[np.isfinite(d)]
    rb = np.random.default_rng([SEED, 11])
    meds = [np.median(d[rb.integers(0, len(d), len(d))]) for _ in range(N_BOOT)]
    summary["Adj3"] = {
        "n": int(len(d)),
        "median_diff": float(np.median(d)),
        "lo": float(np.percentile(meds, 2.5)),
        "hi": float(np.percentile(meds, 97.5)),
    }
    # raw R1 difference (Outcome B test)
    rawd = per["K562_GWPS"]["r"] - per["RPE1"]["r"]
    rawd = rawd[np.isfinite(rawd)]
    rr = np.random.default_rng([SEED, 12])
    meds = [np.median(rawd[rr.integers(0, len(rawd), len(rawd))]) for _ in range(N_BOOT)]
    summary["R1_diff_GWPS_RPE1"] = {
        "median": float(np.median(rawd)),
        "lo": float(np.percentile(meds, 2.5)),
        "hi": float(np.percentile(meds, 97.5)),
    }
    print(f"endpoints done [{time.time() - t0:.0f}s]", flush=True)

    # ---- 3. R4: single-source calibration (N3 machinery) -------------------
    n_test = round(0.3 * len(perts))
    perm = np.random.default_rng([SEED, 21]).permutation(len(perts))
    T_idx, pool = np.sort(perm[:n_test]), np.sort(perm[n_test:])
    pert = np.load(N3 / "pert_part_means.npy", mmap_mode="r")
    ctrl = np.load(N3 / "ctrl_part_means.npy")
    tgt = N3_INDEX["K562_essential"]
    parts = [
        [
            np.asarray(pert[r, j, tgt], dtype=np.float64)[np.ix_(pi, gi)] - ctrl[r, j, tgt][gi]
            for j in range(3)
        ]
        for r in range(pert.shape[0])
    ]
    r4 = []
    for s in SOURCES + ["K562_GWPS_depth"]:
        sv = cb.make_source_view(src[s][0][None], [0], 1.0)
        lv = fast.light_view(sv)
        for r, (fit, e1, e2) in enumerate(parts):
            fe = fast.fast_eval(e1[T_idx], e2[T_idx], {"own": sv.A[T_idx]})
            for k in (20, 50):
                for d_ in range(8):
                    K = np.sort(
                        np.random.default_rng([SEED, 22, k, r, d_]).choice(pool, k, replace=False)
                    )
                    P4, _ = fast.e4(lv, T_idx, K, fit[K])
                    mt = fast.fast_metrics(fe, P4)
                    r4.append(
                        {
                            "source": s,
                            "k": k,
                            "repeat": r,
                            "draw": d_,
                            "M3_own": mt["M3_own"],
                            "M1": mt["M1"],
                        }
                    )
    pd.DataFrame(r4).to_csv(OUT / "r4_calibration.csv", index=False)
    (OUT / "n5_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"done [{time.time() - t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
