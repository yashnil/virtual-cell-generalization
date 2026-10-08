"""N6: independent-lab K562 source compatibility (reports/n6_protocol.md §2–4).

Frozen N5 endpoint (pooled both-sided noise-corrected latent cosine), on the N6
panel. Source-only quantities and the reliability gate are written before the
target is read.

Reproduce: ``uv run python scripts/research_v3/run_n6.py`` (after build_n6.py).
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from virtual_cell.analysis import source_compatibility as sc

REPO = Path(__file__).resolve().parents[2]
AX6 = REPO / "data" / "splits" / "n6_k562"
AX5 = REPO / "data" / "splits" / "n5_k562"
AX3 = REPO / "data" / "splits" / "six_context_n3"
N3 = REPO / "outputs" / "n3" / "data"
N5GW = REPO / "outputs" / "n5" / "gwps"
D6 = REPO / "outputs" / "n6" / "data"
OUT = REPO / "outputs" / "n6"
SEED = 20261006
N_BOOT = 2000
N3_INDEX = {"K562_essential": 0, "RPE1": 1, "Jurkat": 3}
SOURCES = ["VIPerturb_K562", "K562_GWPS", "K562_GWPS_vipdepth", "RPE1", "Jurkat"]


def halves(pm, cm):
    a = [pm[r, 0] - cm[r, 0] for r in range(pm.shape[0])]
    b = [((pm[r, 1] - cm[r, 1]) + (pm[r, 2] - cm[r, 2])) / 2 for r in range(pm.shape[0])]
    return a, b


def from_arrays(delta, pm, cm, cells):
    a, b = halves(pm.astype(np.float64), cm.astype(np.float64))
    return delta.astype(np.float64), a, b, np.asarray(cells, dtype=float)


def main() -> None:
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    perts = (AX6 / "shared_perturbations.txt").read_text().split()
    genes = (AX6 / "shared_genes.txt").read_text().split()
    p5, g5 = (
        (AX5 / "shared_perturbations.txt").read_text().split(),
        (AX5 / "shared_genes.txt").read_text().split(),
    )
    p3, g3 = (
        (AX3 / "shared_perturbations.txt").read_text().split(),
        (AX3 / "shared_genes.txt").read_text().split(),
    )
    i5 = np.array([p5.index(p) for p in perts])
    j5 = np.array([g5.index(g) for g in genes])
    i3 = np.array([p3.index(p) for p in perts])
    j3 = np.array([g3.index(g) for g in genes])

    def n3_context(c):
        D = np.load(N3 / "delta6.npy", mmap_mode="r")
        pert = np.load(N3 / "pert_part_means.npy", mmap_mode="r")
        ctrl = np.load(N3 / "ctrl_part_means.npy")
        pm = np.stack(
            [
                np.stack([np.asarray(pert[r, j, c])[np.ix_(i3, j3)] for j in range(3)])
                for r in range(pert.shape[0])
            ]
        )
        return from_arrays(
            np.asarray(D[c])[np.ix_(i3, j3)],
            pm,
            ctrl[:, :, c][..., j3],
            np.load(N3 / "cell_counts6.npy")[c][i3],
        )

    # ---- 1. source-only -----------------------------------------------------
    src = {
        "VIPerturb_K562": from_arrays(
            np.load(D6 / "viperturb_delta.npy"),
            np.load(D6 / "viperturb_pert_part_means.npy"),
            np.load(D6 / "viperturb_ctrl_part_means.npy"),
            np.load(D6 / "viperturb_cell_counts.npy"),
        ),
        "K562_GWPS": from_arrays(
            np.load(N5GW / "delta.npy")[np.ix_(i5, j5)],
            np.load(N5GW / "pert_part_means.npy")[:, :, i5][..., j5],
            np.load(N5GW / "ctrl_part_means.npy")[..., j5],
            np.load(N5GW / "cell_counts.npy")[i5],
        ),
        "K562_GWPS_vipdepth": from_arrays(
            np.load(D6 / "gwps_vipdepth_delta.npy"),
            np.load(D6 / "gwps_vipdepth_pert_part_means.npy"),
            np.load(D6 / "gwps_vipdepth_ctrl_part_means.npy"),
            np.load(D6 / "gwps_vipdepth_cell_counts.npy"),
        ),
        "RPE1": n3_context(N3_INDEX["RPE1"]),
        "Jurkat": n3_context(N3_INDEX["Jurkat"]),
    }
    rel = {s: sc.split_half_reliability(v[1], v[2]) for s, v in src.items()}
    mag = {s: np.linalg.norm(sc.centre(v[0]), axis=1) for s, v in src.items()}
    va, vb = src["VIPerturb_K562"][1], src["VIPerturb_K562"][2]
    pooled_rel = float(
        np.mean(
            [
                np.sum(sc.centre(a) * sc.centre(b))
                / np.sqrt(np.sum(sc.centre(a) ** 2) * np.sum(sc.centre(b) ** 2))
                for a, b in zip(va, vb, strict=True)
            ]
        )
    )
    den_v = sc.reliable_energy(va, vb)
    rng = np.random.default_rng([SEED, 40])
    idx = rng.integers(0, len(perts), size=(N_BOOT, len(perts)))
    gate = {
        "n_panel": len(perts),
        "viperturb_pooled_split_half_reliability": pooled_rel,
        "viperturb_reliable_energy_sum": float(den_v.sum()),
        "viperturb_reliable_energy_boot_lo": float(np.percentile(den_v[idx].sum(1), 2.5)),
        "median_source_reliability": {s: float(np.nanmedian(rel[s])) for s in src},
        "median_source_cells": {s: float(np.median(src[s][3])) for s in src},
        "panel_B_n": int(np.sum(rel["VIPerturb_K562"] >= 0.10)),
    }
    (OUT / "source_gate.json").write_text(json.dumps(gate, indent=2))
    print(f"source-only stage written [{time.time() - t0:.0f}s]", flush=True)

    # ---- 2. target ----------------------------------------------------------
    T_full, T_a, T_b, _ = n3_context(N3_INDEX["K562_essential"])
    per = {}
    for s, (S_full, S_a, S_b, _) in src.items():
        lt = sc.latent_cosine_terms(S_full, T_full, S_a, S_b, T_a, T_b)
        rc = sc.pooled_raw_cosine(S_full, T_full)
        per[s] = {
            **lt,
            "raw_num": rc["num"],
            "raw_ss": rc["ss"],
            "raw_tt": rc["tt"],
            "r": sc.rowwise_pearson(sc.centre(S_full), sc.centre(T_full)),
            "dir_acc": sc.directional_accuracy(S_full, T_a, T_b),
            "rel": rel[s],
            "log_mag": np.log(np.maximum(mag[s], 1e-12)),
        }
    pd.concat(
        {s: pd.DataFrame(v, index=perts) for s, v in per.items()}, names=["source", "perturbation"]
    ).to_csv(OUT / "per_perturbation_terms.csv")

    def c_of(v, ix):
        return v["num"][ix].sum(-1) / np.sqrt(v["den_S"][ix].sum(-1) * v["den_T"][ix].sum(-1))

    bidx = np.random.default_rng([SEED, 41]).integers(0, len(perts), size=(N_BOOT, len(perts)))
    pd.DataFrame({s: c_of(v, bidx) for s, v in per.items()}).to_csv(
        OUT / "bootstrap_C.csv", index=False
    )
    allix = np.arange(len(perts))
    summary = {
        "C": {s: float(c_of(v, allix)) for s, v in per.items()},
        "R1_median_r": {s: float(np.nanmedian(v["r"])) for s, v in per.items()},
        "R2_raw_cosine": {
            s: float(v["raw_num"].sum() / np.sqrt(v["raw_ss"].sum() * v["raw_tt"].sum()))
            for s, v in per.items()
        },
        "R3_dir_acc": {s: float(np.nanmedian(v["dir_acc"])) for s, v in per.items()},
    }

    # panel B (VIPerturb per-perturbation reliability >= 0.10; same perturbations for all sources)
    pb = np.flatnonzero(rel["VIPerturb_K562"] >= 0.10)
    bb = np.random.default_rng([SEED, 42]).integers(0, len(pb), size=(N_BOOT, len(pb)))
    summary["panel_B"] = {"n": int(len(pb)), "C": {s: float(c_of(v, pb)) for s, v in per.items()}}
    pd.DataFrame({s: c_of(v, pb[bb]) for s, v in per.items()}).to_csv(
        OUT / "bootstrap_C_panelB.csv", index=False
    )

    # Adj-2: FE regression, reference Jurkat
    fe_sources = ["VIPerturb_K562", "K562_GWPS", "RPE1", "Jurkat"]
    y = np.concatenate([per[s]["r"] for s in fe_sources])
    pid = np.tile(np.arange(len(perts)), len(fe_sources))
    lab = np.repeat(np.array(fe_sources), len(perts))
    cov = np.column_stack(
        [
            np.concatenate([per[s]["rel"] for s in fe_sources]),
            np.concatenate([per[s]["rel"] for s in fe_sources]) ** 2,
            np.concatenate([per[s]["log_mag"] for s in fe_sources]),
        ]
    )
    est = sc.fe_regression(y, pid, lab, cov, fe_sources, "Jurkat")
    ci = sc.cluster_bootstrap_fe(
        y, pid, lab, cov, fe_sources, "Jurkat", n_boot=N_BOOT, rng=np.random.default_rng([SEED, 43])
    )
    summary["Adj2_vs_Jurkat"] = {
        s: {"beta_minus_Jurkat": est[s], "lo": ci[s][0], "hi": ci[s][1]} for s in est
    }

    # Adj-3: reliability-matched perturbations, VIPerturb vs Jurkat
    m = np.abs(rel["VIPerturb_K562"] - rel["Jurkat"]) <= 0.05
    d = (per["VIPerturb_K562"]["r"] - per["Jurkat"]["r"])[m]
    d = d[np.isfinite(d)]
    r3 = np.random.default_rng([SEED, 44])
    meds = (
        [np.median(d[r3.integers(0, len(d), len(d))]) for _ in range(N_BOOT)]
        if len(d)
        else [np.nan]
    )
    summary["Adj3_VIP_vs_Jurkat"] = {
        "n": int(len(d)),
        "median_diff": float(np.median(d)) if len(d) else None,
        "lo": float(np.percentile(meds, 2.5)),
        "hi": float(np.percentile(meds, 97.5)),
    }
    (OUT / "n6_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"done [{time.time() - t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
