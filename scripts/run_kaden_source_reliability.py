"""Kaden source-reliability diagnostic v1 (predeclared; changes no model).

Question: is the arch1 / kaden25rpe1 disagreement about the mean perturbation
response caused primarily by Kaden measurement unreliability, and is Kaden
reliable enough to serve as (A) a perturbation-specific beta source and (B) a
context-main-effect (m_hat) source? The two are answered separately.

Everything the decision depends on -- seeds, repeat counts, gene axes,
perturbation panels, quality bands and the CASE A-E rules -- is fixed in the
constants below and written to ``predeclaration.json`` (with this script's
SHA-256) BEFORE any expression value is read. No Arc outcome, leaderboard value,
or frozen model parameter is read or changed; the frozen Arc model is not
evaluated here.

Reproduce: ``uv run python scripts/run_kaden_source_reliability.py``
Outputs:   ``outputs/kaden_source_reliability_v1/``
"""

from __future__ import annotations

import hashlib
import json
import time
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

from virtual_cell.analysis import source_reliability as sr
from virtual_cell.data import arc2026, scperteval

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "scperteval"
CONTROLS = ROOT / "data" / "raw" / "arc2026" / "controls"
SPLITS = ROOT / "data" / "splits"
FROZEN_RESPONSE_GENES = ROOT / "outputs" / "arc_dry_run_v1" / "response_genes.txt"
FROZEN_DRY_RUN_SUMMARY = ROOT / "outputs" / "arc_dry_run_v1" / "summary.json"
OUTDIR = ROOT / "outputs" / "kaden_source_reliability_v1"

# ---------------------------------------------------------------- predeclared
SEED = 20260925
N_REPEATS = 50
N_BLOCKS = 20  # per perturbation group; each half = 10 random blocks
N_CONTROL_BLOCKS = 200
MIN_CELLS = 20  # a perturbation needs >= 10 cells per half to get a reliability
N_NULL = 40  # control pseudo-perturbations per dataset
NULL_CONTROL_FRACTION = 0.5  # at most this share of controls goes into the null
NULL_QUANTILE = 0.95
TOP_K_SIGN = 100
DISATTENUATION_FLOOR = 0.1
N_BOOT = 1000
CHUNK = 20_000

BAND_HIGH = 0.5  # Spearman-Brown reliability of the full estimate
BAND_MODERATE = 0.2

DECISION_RULES = {
    "beta_source": {
        "statistic": "Spearman-Brown reliability of the perturbation-CENTRED response "
        "(the beta_hat quantity) for Kaden-supported Arc targets on the frozen "
        "16,494-gene response axis",
        "unsuitable": "median < 0.2 AND fraction with signal above the null 95th percentile < 0.5",
        "suitable": "median >= 0.5 AND fraction with detected signal >= 0.75",
        "partial": "otherwise",
    },
    "m_hat_source": {
        "statistic": "Spearman-Brown reliability of the main effect (mean over the "
        "natural perturbation panel) on the frozen response axis",
        "reliable": ">= 0.8",
        "unreliable": "< 0.5",
        "moderate": "otherwise",
    },
    "noise_explanation": {
        "statistic": "observed cross-dataset correlation vs the noise ceiling sqrt(R_x R_y)",
        "explained_by_noise": "ceiling < 0.3",
        "not_explained_by_noise": "ceiling >= 0.5 AND |observed| < 0.5 * ceiling",
        "partially_explained": "otherwise",
    },
    "cases": {
        "A": "beta unsuitable AND m_hat unreliable",
        "B": "beta unsuitable AND m_hat reliable",
        "C": "m_hat unreliable AND beta not unsuitable (partial or suitable) with "
        ">= 25% of Kaden Arc targets in the high/moderate band",
        "D": "beta suitable AND m_hat reliable AND the arch1-Kaden main-effect "
        "disagreement not explained by noise",
        "E": "anything else (mixed / inconclusive)",
    },
}

DATASETS = {
    "arch1": ["G_ARC", "G_3"],
    "kaden25rpe1": ["G_ARC", "G_RPE", "G_3", "G_4CTX"],
    "replogle22rpe1": ["G_RPE", "G_3", "G_4CTX"],
    "replogle22k562": ["G_4CTX"],
    "nadig25hepg2": ["G_4CTX"],
    "nadig25jurkat": ["G_4CTX"],
}
AXIS_DOC = {
    "G_ARC": "frozen Arc response space: panel n arch1 n kaden25rpe1 "
    "(outputs/arc_dry_run_v1/response_genes.txt)",
    "G_RPE": "kaden25rpe1 n replogle22rpe1 genes",
    "G_3": "arch1 n kaden25rpe1 n replogle22rpe1 genes (identical axis for the "
    "three-source comparison)",
    "G_4CTX": "frozen four-context 6,640-gene axis (data/splits/four_context_v1)",
}
PANEL_DOC = {
    "all": "every perturbation with >= MIN_CELLS cells",
    "arc_targets": "the dataset's Arc validation targets",
    "shared_arch1_kaden": "perturbations measured in both arch1 and kaden25rpe1",
    "arc_matched": "the Arc targets perturbed in BOTH arch1 and kaden25rpe1",
    "shared_kaden_rpe1": "perturbations measured in both kaden25rpe1 and replogle22rpe1",
    "shared_four_context": "the frozen 1,264 four-context perturbations",
}


def rule(title: str) -> None:
    print("\n" + "=" * 78 + f"\n{title}\n" + "=" * 78, flush=True)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def boot_median_ci(values: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    v = np.asarray(values, dtype=np.float64)
    v = v[np.isfinite(v)]
    if len(v) < 3:
        return float("nan"), float("nan")
    idx = rng.integers(0, len(v), size=(N_BOOT, len(v)))
    meds = np.median(v[idx], axis=1)
    return float(np.quantile(meds, 0.025)), float(np.quantile(meds, 0.975))


def labels_of(name: str) -> list[str]:
    path = RAW / f"{name}_processed_complete.h5ad"
    return sorted(set(scperteval.cell_labels(path).tolist()) - {scperteval.CONTROL_LABEL})


# --------------------------------------------------------------------------
# one dataset
# --------------------------------------------------------------------------


def analyse_dataset(name, axes, panels, rng):
    """Block split-halves for one dataset on each requested axis."""
    path = RAW / f"{name}_processed_complete.h5ad"
    labels = scperteval.cell_labels(path)
    perts = sorted(set(labels.tolist()) - {scperteval.CONTROL_LABEL})
    code = {p: i for i, p in enumerate(perts)}
    group = np.array([code.get(lab, -1) for lab in labels], dtype=np.int64)
    control_code = -2
    group[labels == scperteval.CONTROL_LABEL] = control_code
    n_cells = np.bincount(group[group >= 0], minlength=len(perts))
    eligible = n_cells >= MIN_CELLS

    # pseudo-perturbations carved from the controls (null reference only)
    ctrl_rows = np.flatnonzero(group == control_code)
    median_cells = int(np.median(n_cells[eligible]))
    null_size = max(
        MIN_CELLS, min(median_cells, int(NULL_CONTROL_FRACTION * len(ctrl_rows) / N_NULL))
    )
    chosen = rng.choice(ctrl_rows, size=N_NULL * null_size, replace=False)
    for j in range(N_NULL):
        group[chosen[j * null_size : (j + 1) * null_size]] = len(perts) + j
    n_groups = len(perts) + N_NULL

    layout = sr.assign_blocks(
        group,
        n_groups=n_groups,
        control_code=control_code,
        n_blocks=N_BLOCKS,
        n_control_blocks=N_CONTROL_BLOCKS,
        rng=rng,
    )
    union = sorted(set().union(*(axes[a] for a in DATASETS[name])))
    t0 = time.time()
    sums, counts = sr.accumulate_block_sums(path, genes=union, layout=layout, chunk_size=CHUNK)
    print(
        f"  {name}: {len(perts)} perturbations ({int(eligible.sum())} with >= {MIN_CELLS} "
        f"cells), {len(ctrl_rows):,} controls, null {N_NULL} x {null_size} cells, "
        f"{len(union):,} genes, pass {time.time() - t0:.0f}s",
        flush=True,
    )

    off = layout.control_offset
    blocks = sums[:off].reshape(n_groups, N_BLOCKS, len(union))
    bcounts = counts[:off].reshape(n_groups, N_BLOCKS)
    ctrl_rest = sums[off:]
    ctrl_rest_n = counts[off:]
    # real analysis uses ALL controls: the remaining blocks plus the null blocks
    ctrl_all = np.vstack([ctrl_rest, blocks[len(perts) :].reshape(-1, len(union))])
    ctrl_all_n = np.concatenate([ctrl_rest_n, bcounts[len(perts) :].ravel()])

    # full-data responses (every perturbation, every control: the frozen definition)
    ctrl_mean = ctrl_all.sum(axis=0) / ctrl_all_n.sum()
    with np.errstate(invalid="ignore", divide="ignore"):
        full = blocks[: len(perts)].sum(axis=1) / bcounts[: len(perts)].sum(axis=1)[:, None]
    full_delta = full - ctrl_mean[None, :]

    col = {g: i for i, g in enumerate(union)}
    axis_cols = {a: np.array([col[g] for g in axes[a]]) for a in DATASETS[name]}
    panel_rows = {
        (a, pn): np.array([code[p] for p in plist if p in code and eligible[code[p]]], dtype=int)
        for (a, pn), plist in panels.items()
    }
    all_rows = np.flatnonzero(eligible)

    rho = {a: np.full((N_REPEATS, len(perts)), np.nan) for a in axis_cols}
    rho_c = {a: np.full((N_REPEATS, len(perts)), np.nan) for a in axis_cols}
    energy = {a: np.full((N_REPEATS, len(perts)), np.nan) for a in axis_cols}
    null_rho = {a: np.full((N_REPEATS, N_NULL), np.nan) for a in axis_cols}
    me_stats = {k: [] for k in panel_rows}
    me_halves = {k: [] for k in panel_rows}

    for r in range(N_REPEATS):
        masks = sr.draw_half_masks(n_groups, N_BLOCKS, rng)
        ma, mb, _, _ = sr.half_means(blocks, bcounts, masks)
        ca, cb = sr.control_half_means(ctrl_all, ctrl_all_n, rng)
        na, nb = sr.control_half_means(ctrl_rest, ctrl_rest_n, rng)
        for a, cols in axis_cols.items():
            da = ma[: len(perts), cols] - ca[cols]
            db = mb[: len(perts), cols] - cb[cols]
            ok = eligible
            rho[a][r, ok] = sr.row_pearson(da[ok], db[ok])
            energy[a][r, ok] = (da[ok] * db[ok]).sum(axis=1)
            cen_a = da - da[all_rows].mean(axis=0)
            cen_b = db - db[all_rows].mean(axis=0)
            rho_c[a][r, ok] = sr.row_pearson(cen_a[ok], cen_b[ok])
            pa = ma[len(perts) :, cols] - na[cols]
            pb = mb[len(perts) :, cols] - nb[cols]
            null_rho[a][r] = sr.row_pearson(pa, pb)
            for (pa_axis, pn), rows in panel_rows.items():
                if pa_axis != a or not len(rows):
                    continue
                m1 = da[rows].mean(axis=0)
                m2 = db[rows].mean(axis=0)
                me_stats[(a, pn)].append(
                    (sr.pearson(m1, m2), sr.cosine(m1, m2), sr.norm_ratio(m1, m2))
                )
                me_halves[(a, pn)].append((m1, m2))

    # ------------------------------------------------------------- summaries
    per_pert = []
    null_rows = []
    for a, cols in axis_cols.items():
        mean_rho = np.nanmean(rho[a], axis=0)
        mean_rho_c = np.nanmean(rho_c[a], axis=0)
        sb = sr.spearman_brown(mean_rho)
        sb_c = sr.spearman_brown(mean_rho_c)
        null_sb = sr.spearman_brown(np.nanmean(null_rho[a], axis=0))
        null_q = float(np.quantile(null_sb, NULL_QUANTILE))
        for j, v in enumerate(null_sb):
            null_rows.append(
                {
                    "dataset": name,
                    "axis": a,
                    "pseudo_perturbation": j,
                    "cells": null_size,
                    "spearman_brown": float(v),
                }
            )
        norms = np.linalg.norm(full_delta[:, cols], axis=1)
        for i, p in enumerate(perts):
            per_pert.append(
                {
                    "dataset": name,
                    "axis": a,
                    "perturbation": p,
                    "n_cells": int(n_cells[i]),
                    "eligible": bool(eligible[i]),
                    "rho_half_mean": mean_rho[i],
                    "rho_half_sd": float(np.nanstd(rho[a][:, i])) if eligible[i] else np.nan,
                    "spearman_brown": sb[i],
                    "rho_half_centred_mean": mean_rho_c[i],
                    "spearman_brown_centred": sb_c[i],
                    "signal_energy": float(np.nanmean(energy[a][:, i])) if eligible[i] else np.nan,
                    "delta_norm": float(norms[i]),
                    "null_q95": null_q,
                    "signal_detected": bool(eligible[i] and sb[i] > null_q),
                    "signal_detected_centred": bool(eligible[i] and sb_c[i] > null_q),
                }
            )
    me_rows = []
    for (a, pn), vals in me_stats.items():
        if not vals:
            continue
        arr = np.array(vals)
        mp = float(np.nanmean(arr[:, 0]))
        me_rows.append(
            {
                "dataset": name,
                "axis": a,
                "panel": pn,
                "n_perturbations": len(panel_rows[(a, pn)]),
                "half_pearson_mean": mp,
                "half_pearson_q05": float(np.quantile(arr[:, 0], 0.05)),
                "half_pearson_q95": float(np.quantile(arr[:, 0], 0.95)),
                "half_cosine_mean": float(np.nanmean(arr[:, 1])),
                "half_norm_ratio_mean": float(np.nanmean(arr[:, 2])),
                "half_norm_ratio_q05": float(np.quantile(arr[:, 2], 0.05)),
                "half_norm_ratio_q95": float(np.quantile(arr[:, 2], 0.95)),
                "spearman_brown": float(sr.spearman_brown(mp)),
            }
        )
    return {
        "perts": perts,
        "code": code,
        "n_cells": n_cells,
        "eligible": eligible,
        "union": union,
        "axis_cols": axis_cols,
        "full_delta": full_delta,
        "per_pert": pd.DataFrame(per_pert),
        "null": pd.DataFrame(null_rows),
        "main_effect": pd.DataFrame(me_rows),
        "me_halves": {k: v for k, v in me_halves.items() if v},
        "null_size": null_size,
        "n_controls": int(len(ctrl_rows)),
    }


# --------------------------------------------------------------------------
# cross-dataset comparisons
# --------------------------------------------------------------------------


def main_effect_full(res, axis, plist):
    rows = [res["code"][p] for p in plist if p in res["code"]]
    return res["full_delta"][np.ix_(rows, res["axis_cols"][axis])].mean(axis=0), len(rows)


def centred(res, axis):
    """Perturbation-centred responses over the whole natural panel (beta_hat's
    quantity: frozen beta_hat centres within each source over all perturbations)."""
    d = res["full_delta"][:, res["axis_cols"][axis]]
    return d - d.mean(axis=0, keepdims=True)


def reliability_lookup(res, axis, column):
    df = res["per_pert"]
    df = df[df.axis == axis]
    return dict(zip(df.perturbation, df[column], strict=True))


def per_pert_agreement(rx, ry, axis_x, axis_y, plist, label, rng, arc_set, arc_matched):
    """Per-perturbation agreement between two datasets on the same gene axis."""
    rows = []
    xs = [rx["code"][p] for p in plist]
    ys = [ry["code"][p] for p in plist]
    for kind in ("raw", "centred"):
        if kind == "raw":
            X = rx["full_delta"][np.ix_(xs, rx["axis_cols"][axis_x])]
            Y = ry["full_delta"][np.ix_(ys, ry["axis_cols"][axis_y])]
            col = "spearman_brown"
        else:
            X = centred(rx, axis_x)[xs]
            Y = centred(ry, axis_y)[ys]
            col = "spearman_brown_centred"
        relx = np.array([reliability_lookup(rx, axis_x, col)[p] for p in plist])
        rely = np.array([reliability_lookup(ry, axis_y, col)[p] for p in plist])
        r = sr.row_pearson(X, Y)
        cos = sr.row_cosine(X, Y)
        ratio = np.linalg.norm(X, axis=1) / np.linalg.norm(Y, axis=1)
        sign = sr.sign_agreement_top(X, Y, TOP_K_SIGN)
        ceil = sr.noise_ceiling(relx, rely)
        dis = sr.disattenuated(r, relx, rely, min_reliability=DISATTENUATION_FLOOR)
        for i, p in enumerate(plist):
            rows.append(
                {
                    "comparison": label,
                    "kind": kind,
                    "perturbation": p,
                    "arc_target": p in arc_set,
                    "arc_matched": p in arc_matched,
                    "pearson": r[i],
                    "cosine": cos[i],
                    "norm_ratio_x_over_y": ratio[i],
                    "sign_agreement_top100": sign[i],
                    "reliability_x": relx[i],
                    "reliability_y": rely[i],
                    "noise_ceiling": ceil[i],
                    "disattenuated": dis[i],
                }
            )
    return pd.DataFrame(rows)


def summarise_agreement(df, rng):
    out = []
    for (comp, kind), g in df.groupby(["comparison", "kind"], sort=False):
        for subset, h in (("all shared", g), ("Arc targets", g[g.arc_target])):
            if not len(h):
                continue
            lo, hi = boot_median_ci(h.pearson.values, rng)
            dlo, dhi = boot_median_ci(h.disattenuated.values, rng)
            out.append(
                {
                    "comparison": comp,
                    "kind": kind,
                    "subset": subset,
                    "n": len(h),
                    "median_pearson": float(np.nanmedian(h.pearson)),
                    "pearson_ci_lo": lo,
                    "pearson_ci_hi": hi,
                    "median_cosine": float(np.nanmedian(h.cosine)),
                    "median_norm_ratio": float(np.nanmedian(h.norm_ratio_x_over_y)),
                    "median_sign_agreement": float(np.nanmedian(h.sign_agreement_top100)),
                    "median_reliability_x": float(np.nanmedian(h.reliability_x)),
                    "median_reliability_y": float(np.nanmedian(h.reliability_y)),
                    "median_noise_ceiling": float(np.nanmedian(h.noise_ceiling)),
                    "n_disattenuable": int(np.isfinite(h.disattenuated).sum()),
                    "median_disattenuated": float(np.nanmedian(h.disattenuated))
                    if np.isfinite(h.disattenuated).any()
                    else float("nan"),
                    "disattenuated_ci_lo": dlo,
                    "disattenuated_ci_hi": dhi,
                    "fraction_below_half_ceiling": float(
                        np.mean(h.pearson < 0.5 * h.noise_ceiling)
                    ),
                }
            )
    return pd.DataFrame(out)


def main_effect_agreement(rx, ry, nx, ny, axis, px, py, pn_x, pn_y, label):
    """Full-data main-effect agreement plus split-half noise ceilings."""
    mx, kx = main_effect_full(rx, axis, px)
    my, ky = main_effect_full(ry, axis, py)
    me = pd.concat([rx["main_effect"], ry["main_effect"]])
    relx = me[(me.dataset == nx) & (me.axis == axis) & (me.panel == pn_x)]
    rely = me[(me.dataset == ny) & (me.axis == axis) & (me.panel == pn_y)]
    Rx = float(relx.spearman_brown.iloc[0]) if len(relx) else float("nan")
    Ry = float(rely.spearman_brown.iloc[0]) if len(rely) else float("nan")
    # half-level check: cross-dataset half correlation vs within-dataset halves
    hx = rx["me_halves"].get((axis, pn_x), [])
    hy = ry["me_halves"].get((axis, pn_y), [])
    cross, within_x, within_y = [], [], []
    for (a1, a2), (b1, b2) in zip(hx, hy, strict=False):
        cross.append(np.mean([sr.pearson(a1, b1), sr.pearson(a2, b2)]))
        within_x.append(sr.pearson(a1, a2))
        within_y.append(sr.pearson(b1, b2))
    ceiling = float(sr.noise_ceiling(Rx, Ry))
    obs = sr.pearson(mx, my)
    half_ceiling = (
        float(sr.noise_ceiling(np.mean(within_x), np.mean(within_y))) if cross else np.nan
    )
    return {
        "comparison": label,
        "axis": axis,
        "panel_x": pn_x,
        "panel_y": pn_y,
        "n_x": kx,
        "n_y": ky,
        "cosine": sr.cosine(mx, my),
        "pearson": obs,
        "norm_x": float(np.linalg.norm(mx)),
        "norm_y": float(np.linalg.norm(my)),
        "main_effect_reliability_x": Rx,
        "main_effect_reliability_y": Ry,
        "noise_ceiling": ceiling,
        "pearson_over_ceiling": obs / ceiling if ceiling > 0 else np.nan,
        "half_cross_pearson": float(np.mean(cross)) if cross else np.nan,
        "half_noise_ceiling": half_ceiling,
        "noise_explanation": classify_noise(obs, ceiling),
    }


def classify_noise(observed: float, ceiling: float) -> str:
    if not np.isfinite(ceiling):
        return "not estimable"
    if ceiling < 0.3:
        return "explained_by_noise"
    if ceiling >= 0.5 and abs(observed) < 0.5 * ceiling:
        return "not_explained_by_noise"
    return "partially_explained"


# --------------------------------------------------------------------------


def main() -> None:
    started = time.time()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    rule("0. PREDECLARATION (written before any expression value is read)")
    predeclared = {
        "written": date.today().isoformat(),
        "script": "scripts/run_kaden_source_reliability.py",
        "script_sha256": sha256(Path(__file__)),
        "module_sha256": sha256(ROOT / "src/virtual_cell/analysis/source_reliability.py"),
        "seed": SEED,
        "n_repeats": N_REPEATS,
        "n_blocks": N_BLOCKS,
        "n_control_blocks": N_CONTROL_BLOCKS,
        "min_cells": MIN_CELLS,
        "n_null": N_NULL,
        "null_control_fraction": NULL_CONTROL_FRACTION,
        "null_quantile": NULL_QUANTILE,
        "top_k_sign": TOP_K_SIGN,
        "disattenuation_floor": DISATTENUATION_FLOOR,
        "quality_bands": {"high": BAND_HIGH, "moderate": BAND_MODERATE},
        "decision_rules": DECISION_RULES,
        "axes": AXIS_DOC,
        "panels": PANEL_DOC,
        "datasets": DATASETS,
        "not_done_in_this_phase": [
            "no Arc model evaluated, refitted or modified",
            "no alternative m_hat, tier weight, magnitude calibration or G1 setting tested",
            "no Arc hidden outcome or leaderboard value read",
        ],
    }
    (OUTDIR / "predeclaration.json").write_text(json.dumps(predeclared, indent=2))
    print(f"  written: {OUTDIR / 'predeclaration.json'}")

    rule("1. AXES AND PANELS (identifier presence only)")
    arc_genes = set(arc2026.load_gene_names(CONTROLS))
    arc_targets = list(arc2026.load_pert_counts(CONTROLS)["target_gene"])
    support = pd.read_csv(SPLITS / "arc_target_support_v1.csv")
    var = {
        n: set(scperteval.read_var_names(RAW / f"{n}_processed_complete.h5ad")) for n in DATASETS
    }
    frozen_axis = FROZEN_RESPONSE_GENES.read_text().split()
    recomputed = sorted(arc_genes & var["arch1"] & var["kaden25rpe1"])
    if frozen_axis != recomputed:
        raise SystemExit("frozen response axis does not match panel n arch1 n kaden")
    axes = {
        "G_ARC": frozen_axis,
        "G_RPE": sorted(var["kaden25rpe1"] & var["replogle22rpe1"]),
        "G_3": sorted(var["arch1"] & var["kaden25rpe1"] & var["replogle22rpe1"]),
        "G_4CTX": (SPLITS / "four_context_v1" / "shared_genes.txt").read_text().split(),
    }
    for a, g in axes.items():
        print(f"  {a:7s} {len(g):>6,} genes   {AXIS_DOC[a]}")

    labs = {n: labels_of(n) for n in ("arch1", "kaden25rpe1", "replogle22rpe1")}
    arc_set = set(arc_targets)
    shared_ak = sorted(set(labs["arch1"]) & set(labs["kaden25rpe1"]))
    shared_kr = sorted(set(labs["kaden25rpe1"]) & set(labs["replogle22rpe1"]))
    arc_matched = sorted(set(shared_ak) & arc_set)
    four = (SPLITS / "four_context_v1" / "shared_perturbations.txt").read_text().split()
    print(
        f"  arch1 n kaden perturbations: {len(shared_ak)} (Arc targets among them: "
        f"{len(arc_matched)})"
    )
    print(f"  kaden n replogle22rpe1 perturbations: {len(shared_kr)}")

    def panels_for(name):
        own = labs.get(name) or labels_of(name)
        p = {}
        for a in DATASETS[name]:
            p[(a, "all")] = own
        if name in ("arch1", "kaden25rpe1"):
            p[("G_ARC", "arc_targets")] = sorted(set(own) & arc_set)
            p[("G_ARC", "shared_arch1_kaden")] = shared_ak
            p[("G_ARC", "arc_matched")] = arc_matched
        if name in ("kaden25rpe1", "replogle22rpe1"):
            p[("G_RPE", "shared_kaden_rpe1")] = shared_kr
        if "G_4CTX" in DATASETS[name] and name != "kaden25rpe1":
            p[("G_4CTX", "shared_four_context")] = four
        return p

    rule("2. PER-DATASET SPLIT-HALF RELIABILITY (50 block split-halves)")
    results = {}
    for name in DATASETS:
        results[name] = analyse_dataset(name, axes, panels_for(name), rng)

    per_pert = pd.concat([r["per_pert"] for r in results.values()], ignore_index=True)
    null = pd.concat([r["null"] for r in results.values()], ignore_index=True)
    main_eff = pd.concat([r["main_effect"] for r in results.values()], ignore_index=True)
    arc_tier = dict(zip(support.arc_target, support.support_tier, strict=True))
    per_pert["arc_target"] = per_pert.perturbation.isin(arc_set)
    per_pert["arc_tier"] = per_pert.perturbation.map(arc_tier)
    per_pert.to_csv(OUTDIR / "per_perturbation_reliability.csv", index=False)
    null.to_csv(OUTDIR / "null_reliability.csv", index=False)
    main_eff.to_csv(OUTDIR / "main_effect_reliability.csv", index=False)

    dist = []
    for (d, a), g in per_pert[per_pert.eligible].groupby(["dataset", "axis"], sort=False):
        q = g.spearman_brown.quantile([0.1, 0.25, 0.5, 0.75, 0.9]).values
        nq = float(g.null_q95.iloc[0])
        dist.append(
            {
                "dataset": d,
                "axis": a,
                "n_perturbations": len(g),
                "median_cells": float(g.n_cells.median()),
                "sb_q10": q[0],
                "sb_q25": q[1],
                "sb_median": q[2],
                "sb_q75": q[3],
                "sb_q90": q[4],
                "sb_centred_median": float(g.spearman_brown_centred.median()),
                "median_delta_norm": float(g.delta_norm.median()),
                "median_signal_energy": float(g.signal_energy.median()),
                "null_q95": nq,
                "null_median": float(
                    null[(null.dataset == d) & (null.axis == a)].spearman_brown.median()
                ),
                "fraction_detected": float(g.signal_detected.mean()),
                "fraction_high": float((g.spearman_brown >= BAND_HIGH).mean()),
                "fraction_low": float((g.spearman_brown < BAND_MODERATE).mean()),
            }
        )
    dist = pd.DataFrame(dist)
    dist.to_csv(OUTDIR / "reliability_distribution.csv", index=False)
    print(dist.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print()
    print(main_eff.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    rule("3. MAIN-EFFECT AGREEMENT ACROSS SOURCES (full data + noise ceilings)")
    A, K, R = results["arch1"], results["kaden25rpe1"], results["replogle22rpe1"]
    arc_a = sorted(set(labs["arch1"]) & arc_set)
    arc_k = sorted(set(labs["kaden25rpe1"]) & arc_set)
    me_rows = [
        main_effect_agreement(
            A,
            K,
            "arch1",
            "kaden25rpe1",
            "G_ARC",
            labs["arch1"],
            labs["kaden25rpe1"],
            "all",
            "all",
            "arch1 vs Kaden, natural panels (frozen m_hat inputs)",
        ),
        main_effect_agreement(
            A,
            K,
            "arch1",
            "kaden25rpe1",
            "G_ARC",
            arc_a,
            arc_k,
            "arc_targets",
            "arc_targets",
            "arch1 vs Kaden, Arc-target subsets",
        ),
        main_effect_agreement(
            A,
            K,
            "arch1",
            "kaden25rpe1",
            "G_ARC",
            shared_ak,
            shared_ak,
            "shared_arch1_kaden",
            "shared_arch1_kaden",
            "arch1 vs Kaden, all shared perturbations",
        ),
        main_effect_agreement(
            A,
            K,
            "arch1",
            "kaden25rpe1",
            "G_ARC",
            arc_matched,
            arc_matched,
            "arc_matched",
            "arc_matched",
            "arch1 vs Kaden, matched Arc targets",
        ),
        main_effect_agreement(
            K,
            R,
            "kaden25rpe1",
            "replogle22rpe1",
            "G_RPE",
            labs["kaden25rpe1"],
            labs["replogle22rpe1"],
            "all",
            "all",
            "Kaden vs Replogle RPE1, natural panels",
        ),
        main_effect_agreement(
            K,
            R,
            "kaden25rpe1",
            "replogle22rpe1",
            "G_RPE",
            shared_kr,
            shared_kr,
            "shared_kaden_rpe1",
            "shared_kaden_rpe1",
            "Kaden vs Replogle RPE1, shared perturbations",
        ),
    ]
    four_ctx = ["replogle22k562", "replogle22rpe1", "nadig25hepg2", "nadig25jurkat"]
    for i, x in enumerate(four_ctx):
        for y in four_ctx[i + 1 :]:
            me_rows.append(
                main_effect_agreement(
                    results[x],
                    results[y],
                    x,
                    y,
                    "G_4CTX",
                    four,
                    four,
                    "shared_four_context",
                    "shared_four_context",
                    f"{x} vs {y}, four-context reference",
                )
            )
    me_agree = pd.DataFrame(me_rows)
    me_agree.to_csv(OUTDIR / "main_effect_agreement.csv", index=False)
    print(me_agree.drop(columns=["axis"]).to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    frozen = json.loads(FROZEN_DRY_RUN_SUMMARY.read_text())
    checks = {
        "cosine_natural": (me_rows[0]["cosine"], frozen["source_main_effect_cosine"]),
        "cosine_arc_subset": (me_rows[1]["cosine"], frozen["source_main_effect_cosine_arc_subset"]),
        "cosine_matched": (me_rows[3]["cosine"], frozen["source_main_effect_cosine_matched"]),
    }
    print("\n  reproduction of the frozen dry-run main-effect cosines:")
    for k, (mine, theirs) in checks.items():
        print(
            f"    {k:18s} recomputed {mine:+.6f}   frozen {theirs:+.6f}   "
            f"|diff| {abs(mine - theirs):.2e}"
        )

    rule("4. PER-PERTURBATION AGREEMENT ON SHARED PERTURBATIONS")
    agree_kr = per_pert_agreement(
        K,
        R,
        "G_RPE",
        "G_RPE",
        shared_kr,
        "Kaden vs Replogle RPE1 (same cell line)",
        rng,
        arc_set,
        set(arc_matched),
    )
    agree_ak = per_pert_agreement(
        A, K, "G_ARC", "G_ARC", shared_ak, "arch1 vs Kaden", rng, arc_set, set(arc_matched)
    )
    agree = pd.concat([agree_kr, agree_ak], ignore_index=True)
    agree.to_csv(OUTDIR / "per_perturbation_agreement.csv", index=False)
    agree_summary = summarise_agreement(agree, rng)
    agree_summary.to_csv(OUTDIR / "agreement_summary.csv", index=False)
    print(agree_summary.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    rule("5. ARC-SUPPORTED TARGETS: SOURCE RELIABILITY (descriptive; nothing removed)")
    supported = support[support.support_tier > 0]
    rows = []
    for t in supported.itertuples():
        for src in str(t.datasets_perturbed).split("|"):
            if src not in ("arch1", "kaden25rpe1"):
                continue
            rec = per_pert[
                (per_pert.dataset == src)
                & (per_pert.axis == "G_ARC")
                & (per_pert.perturbation == t.arc_target)
            ]
            if not len(rec):
                continue
            rec = rec.iloc[0]
            rows.append(
                {
                    "target": t.arc_target,
                    "tier": int(t.support_tier),
                    "source_datasets": t.datasets_perturbed,
                    "dataset": src,
                    "cells": int(rec.n_cells),
                    "response_norm": rec.delta_norm,
                    "signal_energy": rec.signal_energy,
                    "reliability": rec.spearman_brown,
                    "reliability_centred": rec.spearman_brown_centred,
                    "null_q95": rec.null_q95,
                    "signal_detected": bool(rec.signal_detected),
                    "signal_detected_centred": bool(rec.signal_detected_centred),
                    "quality_flag": sr.quality_band(
                        rec.spearman_brown_centred, high=BAND_HIGH, moderate=BAND_MODERATE
                    ),
                }
            )
    arc_table = pd.DataFrame(rows).sort_values(["dataset", "tier", "target"])
    arc_table.to_csv(OUTDIR / "arc_supported_target_reliability.csv", index=False)
    print(arc_table.groupby(["dataset", "quality_flag"]).size().to_string())

    rule("6. PREDECLARED CLASSIFICATION")
    kaden_arc = arc_table[arc_table.dataset == "kaden25rpe1"]
    beta_med = float(kaden_arc.reliability_centred.median())
    beta_det = float(kaden_arc.signal_detected_centred.mean())
    if beta_med < 0.2 and beta_det < 0.5:
        beta_verdict = "unsuitable"
    elif beta_med >= 0.5 and beta_det >= 0.75:
        beta_verdict = "suitable"
    else:
        beta_verdict = "partial"
    k_me = main_eff[
        (main_eff.dataset == "kaden25rpe1") & (main_eff.axis == "G_ARC") & (main_eff.panel == "all")
    ]
    m_rel = float(k_me.spearman_brown.iloc[0])
    m_verdict = "reliable" if m_rel >= 0.8 else ("unreliable" if m_rel < 0.5 else "moderate")
    noise = me_rows[0]["noise_explanation"]
    qualified = float(kaden_arc.quality_flag.isin(["high", "moderate"]).mean())
    if beta_verdict == "unsuitable" and m_verdict == "unreliable":
        case = "A"
    elif beta_verdict == "unsuitable" and m_verdict == "reliable":
        case = "B"
    elif m_verdict == "unreliable" and beta_verdict != "unsuitable" and qualified >= 0.25:
        case = "C"
    elif (
        beta_verdict == "suitable" and m_verdict == "reliable" and noise == "not_explained_by_noise"
    ):
        case = "D"
    else:
        case = "E"
    verdict = {
        "beta_source": {
            "median_centred_reliability_kaden_arc_targets": beta_med,
            "fraction_detected_centred": beta_det,
            "fraction_high_or_moderate": qualified,
            "verdict": beta_verdict,
        },
        "m_hat_source": {"kaden_main_effect_reliability": m_rel, "verdict": m_verdict},
        "arch1_kaden_main_effect_noise_explanation": noise,
        "case": case,
    }
    print(json.dumps(verdict, indent=2))

    summary = {
        "generated": date.today().isoformat(),
        "runtime_minutes": (time.time() - started) / 60,
        "axes": {a: len(g) for a, g in axes.items()},
        "n_shared_arch1_kaden": len(shared_ak),
        "n_shared_kaden_rpe1": len(shared_kr),
        "n_arc_matched": len(arc_matched),
        "arc_matched_targets": arc_matched,
        "null_sizes": {n: r["null_size"] for n, r in results.items()},
        "n_controls": {n: r["n_controls"] for n, r in results.items()},
        "frozen_reproduction": {k: {"recomputed": a, "frozen": b} for k, (a, b) in checks.items()},
        "verdict": verdict,
    }
    (OUTDIR / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    print(f"\n  done in {summary['runtime_minutes']:.1f} min")


if __name__ == "__main__":
    main()
