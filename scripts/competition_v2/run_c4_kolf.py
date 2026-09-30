"""C4 KOLF2.1J qualification: statistics, reliability audit, held-out transfer, decision.

Predeclaration: ``reports/competition_v2/c4_predeclaration.md`` (SHA-256 checked first).

    uv run python scripts/competition_v2/run_c4_kolf.py --stage stats
    uv run python scripts/competition_v2/run_c4_kolf.py --stage reliability
    uv run python scripts/competition_v2/run_c4_kolf.py --stage transfer [--jobs 5]
    uv run python scripts/competition_v2/run_c4_kolf.py --stage decide

Nothing here changes C1: KOLF is only ever added as one more equal-weight source in the
public held-out folds. Outputs: ``outputs/competition_v2/c4_kolf/``.
"""

from __future__ import annotations

import os

os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import argparse  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import analyse_c2 as a2  # noqa: E402
import run_c1_public_folds as c1f  # noqa: E402
import run_c3_fusion as c3  # noqa: E402

from virtual_cell.arc import metrics as M  # noqa: E402
from virtual_cell.competition_v2 import (  # noqa: E402
    evaluation,
    fusion_c3,
    generator,
    licensing,
    sources,
    sources_c4,
)
from virtual_cell.competition_v2.atlas import source_effect  # noqa: E402
from virtual_cell.visualization import style  # noqa: E402

OUT = ROOT / "outputs" / "competition_v2" / "c4_kolf"
FIG = ROOT / "reports" / "competition_v2" / "figures"
PREDECL = ROOT / "reports" / "competition_v2" / "c4_predeclaration.md"
PREDECL_SHA = "bcedc1d07d7b5f6f7fc1dde35223011aef1c15deea35bdf2644dc6166f00b397"
KOLF_FILE = ROOT / "data" / "raw" / "competition_v2" / "kolf" / "KOLF_Pan_Genome_QC_Filtered.h5ad"
KOLF_MD5 = "afd30fde1e6ad32969c29868394385d1"
STATS = OUT / "KOLF_full_statistics.npz"
HALVES = OUT / "KOLF_split_halves.npz"
C2_FOLDS = ROOT / "outputs" / "competition_v2" / "c2_calibration" / "folds"
C3_OUT = ROOT / "outputs" / "competition_v2" / "c3_fusion"
ATLASES = ["H1", "K562", "CD4"]
KOLF = "KOLF"


def log(msg: str) -> None:
    print(msg, flush=True)


def retained() -> list[str]:
    """The frozen retained-target set: the K562 GWPS source's rows (every C1 cell source
    keeps exactly these rows, built by ``prepare_c1_sources.retained``)."""
    k562 = sources.load_source(c1f.SRC / "K562_GWPS_CPM_full_statistics.npz")
    official = pd.read_csv(ROOT / "data/raw/arc2026/controls/pert_counts.csv").target_gene
    ret = [str(t) for t in k562.targets]
    if ret[:300] != official.astype(str).tolist() or len(ret) != 437:
        raise SystemExit("unexpected frozen retained-target set")
    return ret


def load_kolf() -> sources.SourceStats:
    return sources.load_source(STATS)


# ----------------------------------------------------------------------------- stats
def stage_stats() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    log_rows = [
        json.loads(x)
        for x in (ROOT / "data/provenance/competition_v2/download_log.jsonl")
        .read_text()
        .splitlines()
        if x.strip()
    ]
    ok = [r for r in log_rows if r.get("md5") == KOLF_MD5 and r.get("md5_ok")]
    if not ok:
        raise SystemExit("KOLF download is not MD5-verified in the download log")
    t0 = time.time()
    stats, halves = sources_c4.kolf_statistics(
        KOLF_FILE, retained(), seed=generator.SEED, half_seed=generator.SEED + 31, log=log
    )
    np.savez_compressed(STATS, **stats)
    np.savez_compressed(HALVES, **halves)
    log(f"KOLF statistics saved ({time.time() - t0:.0f}s)")


# ----------------------------------------------------------------------------- reliability
def _on_target(src: sources.SourceStats) -> dict:
    out = {}
    gi = {g: i for i, g in enumerate(src.genes)}
    for k, t in enumerate(src.targets):
        if t in gi and src.n_cells[k] >= 20 and src.mean_cpm is not None:
            c = (
                src.control_mean_cpm[k, gi[t]]
                if src.control_mean_cpm.ndim == 2
                else src.control_mean_cpm[gi[t]]
            )
            if c >= 1:
                out[t] = float(np.log2((src.mean_cpm[k, gi[t]] + 1) / (c + 1)))
    return out


def _null_cells(atlas: str, n_pseudo: int, size: int) -> dict:
    x, labels, _ = evaluation.load_cells(c1f.SRC / f"fold_cells_{atlas}.npz")
    ctrl = np.flatnonzero(labels == sources.CONTROL_LABEL)
    rng = np.random.default_rng(generator.SEED + 1)
    pick = rng.choice(ctrl, n_pseudo * size, replace=False)
    lab = labels.copy().astype(object)
    lab[pick] = [f"__null_{i // size}" for i in range(len(pick))]
    keep = np.ones(len(labels), dtype=bool)
    names = [f"__null_{i}" for i in range(n_pseudo)]
    res = fusion_c3.split_half(
        x[keep], lab.astype(str), names, sources.CONTROL_LABEL, seed=generator.SEED + 2
    )
    return res["rel"]


def stage_reliability() -> None:
    halves = dict(np.load(HALVES, allow_pickle=True))
    kolf = sources_c4.reliability_from_halves(
        halves["mean_a"], halves["mean_b"], halves["control_mean_cpm"], list(halves["kept"])
    )
    null = sources_c4.reliability_from_halves(
        halves["null_a"], halves["null_b"], halves["control_mean_cpm"], list(halves["null_kept"])
    )
    c3rel = pd.read_csv(C3_OUT / "source_reliability.csv")
    src_k = load_kolf()
    rows = [{"source": "KOLF", "target": t, "reliability": r} for t, r in kolf["rel"].items()]
    rows += c3rel.to_dict("records")
    rel = pd.DataFrame(rows)
    rel.to_csv(OUT / "reliability_per_target.csv", index=False)
    nulls = {"KOLF": list(null["rel"].values())}
    for atlas in ("H1", "K562"):
        nulls[atlas] = list(_null_cells(atlas, 40, 150).values())
    srcs = {
        "KOLF": src_k,
        "H1": sources.load_source(c1f.SRC / "H1_2025_full_statistics.npz"),
        "K562": sources.load_source(c1f.SRC / "K562_GWPS_CPM_full_statistics.npz"),
    }
    official = set(pd.read_csv(ROOT / "data/raw/arc2026/controls/pert_counts.csv").target_gene)
    summary = []
    for name, grp in rel.groupby("source"):
        r = grp.reliability.to_numpy()
        row = {
            "source": name,
            "protocol": "C3 split-half SB (centred log2fc, ctrl CPM>=5)"
            if name != "CD4"
            else "C3 SE-based 1-SE^2/var(lfc)",
            "n_targets": len(r),
            "median_reliability": float(np.median(r)),
            "q25": float(np.quantile(r, 0.25)),
            "q75": float(np.quantile(r, 0.75)),
            "arc_targets_median": float(np.median(grp[grp.target.isin(official)].reliability))
            if grp.target.isin(official).any()
            else np.nan,
        }
        if name in nulls:
            thr = float(np.quantile(nulls[name], 0.95))
            row.update(
                null_p95=thr,
                null_median=float(np.median(nulls[name])),
                detectable_fraction=float(np.mean(r > thr)),
            )
        if name in srcs:
            s = srcs[name]
            usable = s.n_cells >= 20
            row["cells_per_target_median"] = float(np.median(s.n_cells[usable]))
            full, m = source_effect(s, s.targets, s.genes, space="log2fc")
            norms = np.linalg.norm(np.where(m, full, 0)[usable], axis=1)
            row["response_norm_median"] = float(np.median(norms))
            ot = _on_target(s)
            row["on_target_log2fc_median"] = float(np.median(list(ot.values()))) if ot else np.nan
            row["on_target_n"] = len(ot)
        summary.append(row)
    # main-effect reliability: SB of the half-level Pearson over all KOLF targets
    a, b = np.asarray(halves["mean_a"]), np.asarray(halves["mean_b"])
    ctrl = np.asarray(halves["control_mean_cpm"])
    keep = ctrl >= fusion_c3.MIN_CPM
    la = np.log2((a[:, keep] + 1) / (ctrl[keep] + 1))
    lb = np.log2((b[:, keep] + 1) / (ctrl[keep] + 1))
    ma, mb = la.mean(0), lb.mean(0)
    r_main = float(np.corrcoef(ma, mb)[0, 1])
    quoted = [
        {
            "source": "Kaden RPE1 (quoted)",
            "protocol": "research-track split-half, G_ARC",
            "median_reliability": 0.167,
            "cells_per_target_median": 400,
        },
        {
            "source": "arch1 H1 (quoted)",
            "protocol": "research-track split-half, G_ARC",
            "median_reliability": 0.906,
            "cells_per_target_median": 1045,
        },
    ]
    table = pd.DataFrame(summary + quoted)
    table.to_csv(OUT / "reliability_summary.csv", index=False)
    (OUT / "reliability_extra.json").write_text(
        json.dumps(
            {
                "kolf_main_effect_half_pearson": r_main,
                "kolf_main_effect_sb": 2 * r_main / (1 + r_main),
                "nulls": {
                    k: {"n": len(v), "p95": float(np.quantile(v, 0.95))} for k, v in nulls.items()
                },
            },
            indent=2,
        )
    )
    log(table.round(3).to_string())


# ----------------------------------------------------------------------------- transfer
def predictor_set(atlas: str, green: dict, cd4: dict, kolf: sources.SourceStats | None):
    srcs = {k: v for k, v in green.items() if k != atlas}
    if kolf is not None:
        srcs[KOLF] = kolf
    use_cd4 = cd4 if atlas != "CD4" else None
    names = [n if n != KOLF else sources_c4.KOLF_NAME for n in srcs] + (
        ["CD4"] if use_cd4 is not None else []
    )
    licensing.assert_sources_allowed(names)
    return srcs, use_cd4


def stage_transfer(jobs: int) -> None:
    green, cd4, kolf = c1f.load_green(), c1f.load_cd4(), load_kolf()
    mean_rows, pairs_rows = [], []
    for atlas in ATLASES:
        targets, genes = c3.panel(atlas)
        t_eff, t_mask = c3.truth(atlas, targets, genes, green, cd4)
        valid = t_mask & ~np.isin(genes, targets)[None, :]
        preds = {}
        for arm, extra in (("C1a", None), ("C1a_KOLF", kolf)):
            srcs, use_cd4 = predictor_set(atlas, green, cd4, extra)
            comp = fusion_c3.components(srcs, use_cd4, targets, genes)
            pred = fusion_c3.combine(comp, "log2fc", None)
            anym = np.zeros_like(valid)
            for _, _, m in fusion_c3.source_list(comp):
                anym |= m
            preds[arm] = (pred, anym)
            if arm == "C1a_KOLF":
                singles = {n: (e, m) for n, e, m in fusion_c3.source_list(comp)}
        # both arms on the same cells: those the C1 set predicts, plus KOLF-only cells separately
        base_valid = valid & preds["C1a"][1]
        for arm, (pred, anym) in preds.items():
            for scope, v in (("c1_cells", base_valid), ("all_cells", valid & anym)):
                per, sse, sst = fusion_c3.row_metrics(pred, t_eff, v)
                row = {
                    "heldout": atlas,
                    "arm": arm,
                    "scope": scope,
                    **fusion_c3.summarise(per, sse, sst),
                }
                if atlas == "CD4":
                    excl = np.isin(genes, targets)
                    row["effect_pds"] = float(M.pds_cosine(pred, t_eff, exclude=excl).mean())
                mean_rows.append(row)
        # section M: single-source transfer and agreement with KOLF
        for name, (e, m) in singles.items():
            per, sse, sst = fusion_c3.row_metrics(np.where(m, e, 0), t_eff, valid & m)
            pairs_rows.append(
                {
                    "heldout": atlas,
                    "source": name,
                    "kind": "single_source_transfer",
                    **fusion_c3.summarise(per, sse, sst),
                }
            )
        for name, (e, m) in singles.items():
            if name == KOLF:
                continue
            ek, mk = singles[KOLF]
            v = mk & m & ~np.isin(genes, targets)[None, :]
            per, _, _ = fusion_c3.row_metrics(np.where(m, e, 0), np.where(mk, ek, 0), v)
            pairs_rows.append(
                {
                    "heldout": atlas,
                    "source": f"{name}~KOLF",
                    "kind": "source_agreement",
                    "cosine": float(np.nanmean(per["cosine"])),
                    "n_targets": int(np.isfinite(per["cosine"]).sum()),
                }
            )
        log(f"[mean] {atlas} done")
    pd.DataFrame(mean_rows).to_csv(OUT / "transfer_mean_level.csv", index=False)
    pd.DataFrame(pairs_rows).to_csv(OUT / "kolf_single_source_and_agreement.csv", index=False)

    for name in ("H1", "K562"):
        t0 = time.time()
        fold = c1f.fold_h1(False) if name == "H1" else c1f.fold_k562(False)
        scorer = evaluation.FoldScorer(fold, jobs=jobs)
        pairs = generator.promoter_pairs(
            generator.gencode_tss(c1f.RAW / "gencode.v47.annotation.gtf.gz"),
            fold.targets,
            fold.genes,
        )
        rows = {}
        (r, _), _ = scorer.split_half(generator.SEED + 2)
        rows["ANCHOR_split_half"] = r
        real_cpm = np.stack([(b / b.sum(1, keepdims=True)).mean(0) for b in fold.real])
        real_bulk = np.stack([b.sum(0) / b.sum() for b in fold.real])
        n = len(fold.targets)
        moments = {
            "ANCHOR_mean_response": (
                np.tile(real_cpm.mean(0), (n, 1)),
                np.tile(real_bulk.mean(0), (n, 1)),
            )
        }
        for arm, extra in (("C1a", None), ("C1a_KOLF", kolf)):
            srcs, use_cd4 = predictor_set(name, green, cd4, extra)
            comp = fusion_c3.components(srcs, use_cd4, fold.targets, fold.genes)
            eff = fusion_c3.fused(comp, fold.controls)
            moments[arm] = generator.expected_moments(
                eff, fold.controls, fold.targets, fold.genes, pairs=pairs
            )
        for arm, (p_cpm, p_bulk) in moments.items():
            st = scorer.emit(
                arm, p_cpm / p_cpm.sum(1, keepdims=True), p_bulk / p_bulk.sum(1, keepdims=True)
            )
            rows[arm], _ = scorer.members(st)
            log(
                f"[{name}] {arm:22s} PDS={rows[arm]['pds_cosine']:.4f} "
                f"FID={rows[arm]['de_wilcoxon_direction_fidelity_yield_raw']:.4f} "
                f"({time.time() - t0:.0f}s)"
            )
        raw = pd.DataFrame(rows).T
        raw.index.name = "arm"
        out = OUT / "folds" / name
        out.mkdir(parents=True, exist_ok=True)
        raw.to_csv(out / "scores_raw.csv")
        c2 = pd.read_csv(C2_FOLDS / name / "phase1_scores_raw.csv", index_col=0)
        diff = {
            m: abs(float(raw.loc["C1a", m]) - float(c2.loc["G0_a1.00", m]))
            for m in evaluation.MEMBERS
        }
        ok = max(diff.values()) <= 1e-12
        (out / "reproduction.json").write_text(json.dumps({"diffs": diff, "pass": ok}, indent=2))
        log(f"[{name}] C1a reproduces C2 G0_a1.00: {ok}")
        if not ok:
            raise SystemExit("C1a reproduction failed")


# ----------------------------------------------------------------------------- decide
def stage_decide() -> None:
    vcc = {
        f: a2.scale(
            pd.read_csv(OUT / "folds" / f / "scores_raw.csv", index_col=0), "ANCHOR_mean_response"
        )
        for f in ("H1", "K562")
    }
    mean = pd.read_csv(OUT / "transfer_mean_level.csv")
    ml = mean[mean.scope == "all_cells"].set_index(["heldout", "arm"])
    base, cand = "C1a", "C1a_KOLF"
    gains = {f: float(vcc[f].loc[cand, "Overall"] - vcc[f].loc[base, "Overall"]) for f in vcc}
    total = sum(gains.values())
    pds_folds = {f: float(vcc[f].loc[cand, "PDS"] / vcc[f].loc[base, "PDS"]) for f in vcc}
    pds_folds["CD4_effect"] = float(
        ml.loc[("CD4", cand), "effect_pds"] / ml.loc[("CD4", base), "effect_pds"]
    )
    pds_mean_c = np.mean(
        [
            vcc["H1"].loc[cand, "PDS"],
            vcc["K562"].loc[cand, "PDS"],
            ml.loc[("CD4", cand), "effect_pds"],
        ]
    )
    pds_mean_b = np.mean(
        [
            vcc["H1"].loc[base, "PDS"],
            vcc["K562"].loc[base, "PDS"],
            ml.loc[("CD4", base), "effect_pds"],
        ]
    )
    fm = lambda col, arm: float(np.mean([ml.loc[(h, arm), col] for h in ATLASES]))  # noqa: E731
    overall = {a: float(np.mean([vcc[f].loc[a, "Overall"] for f in vcc])) for a in (base, cand)}
    rel = pd.read_csv(OUT / "reliability_summary.csv").set_index("source")
    crit = {
        "1_mean_overall_up": overall[cand] > overall[base],
        "2_pds_preserved": bool(
            pds_mean_c >= pds_mean_b and all(v >= 0.99 for v in pds_folds.values())
        ),
        "3_cosine_up": fm("cosine", cand) > fm("cosine", base),
        "4_sign_acc_up": fm("sign_acc", cand) > fm("sign_acc", base),
        "5_both_folds_up": all(v > 0 for v in gains.values()),
        "6_no_anomalous_fold": bool(total > 0 and all(v / total <= 0.75 for v in gains.values())),
        "7_green": licensing.STATUS[sources_c4.KOLF_NAME] == licensing.GREEN,
        "I_signal_above_null": bool(
            rel.loc["KOLF", "median_reliability"] > rel.loc["KOLF", "null_p95"]
        ),
    }
    decision = {
        "kolf_pass": all(crit.values()),
        "criteria": {k: bool(v) for k, v in crit.items()},
        "overall": overall,
        "fold_gain": gains,
        "pds_retention": pds_folds,
        "cosine": {a: fm("cosine", a) for a in (base, cand)},
        "sign_acc": {a: fm("sign_acc", a) for a in (base, cand)},
    }
    (OUT / "c4_decision.json").write_text(json.dumps(decision, indent=2, default=float))
    rows = [{"fold": f, "arm": a, **vcc[f].loc[a].to_dict()} for f in vcc for a in (base, cand)]
    pd.DataFrame(rows).to_csv(OUT / "vcc_scaled.csv", index=False)
    figures(vcc, mean)
    print(json.dumps(decision, indent=1, default=float))


def figures(vcc, mean) -> None:
    style.apply()
    src = FIG / "sources"
    src.mkdir(parents=True, exist_ok=True)

    def save(fig, name):
        fig.tight_layout()
        for ext in ("png", "svg"):
            fig.savefig(FIG / f"{name}.{ext}", dpi=200)
        plt.close(fig)

    # P1 coverage heatmap
    cov = pd.read_csv(ROOT / "data/provenance/competition_v2/c4/kolf_arc_coverage.csv").set_index(
        "arc_target"
    )
    sup = pd.read_csv(ROOT / "data/splits/arc_target_support_license_clean_v1.csv").set_index(
        "arc_target"
    )
    mat = pd.DataFrame(
        {
            "K562 GWPS": sup.K562_GWPS__usable,
            "H1 2025": sup.H1_2025_full__usable,
            "CD4 DE": sup.CD4_DE__usable,
            "KOLF iPSC (GREEN)": cov.kolf_usable,
            "VIPerturb K562 (GREEN, not dl)": cov.viperturb_cells >= 20,
            "Jurkat library (UNKNOWN)": cov.jurkat_in_library,
            "X-Atlas (BLOCKED)": sup.n_usable_contexts_c0 > sup.n_usable_contexts_c1,
        }
    ).astype(int)
    order = mat.sum(1).sort_values(ascending=False).index
    mat = mat.loc[order]
    mat.to_csv(src / "c4_P1_coverage.csv")
    fig, ax = plt.subplots(figsize=(style.FULL_WIDTH, 2.8))
    ax.imshow(mat.T.to_numpy(), aspect="auto", cmap="Greys", interpolation="nearest")
    ax.set_yticks(range(mat.shape[1]), mat.columns, fontsize=7)
    ax.set_xlabel("300 Arc targets (sorted by number of sources)", fontsize=7)
    ax.set_xticks([])
    ax.set_title("Direct-evidence coverage of the Arc panel (usable = ≥ 20 cells)", fontsize=9)
    save(fig, "c4_P1_coverage")
    # P2 reliability comparison
    rel = pd.read_csv(OUT / "reliability_per_target.csv")
    rel.to_csv(src / "c4_P2_reliability.csv", index=False)
    summ = pd.read_csv(OUT / "reliability_summary.csv")
    fig, ax = plt.subplots(figsize=(style.HALF_WIDTH * 1.6, 2.8))
    names = ["H1", "KOLF", "K562", "CD4"]
    ax.boxplot(
        [rel[rel.source == n].reliability for n in names], tick_labels=names, showfliers=False
    )
    for i, q in enumerate(["arch1 H1 (quoted)", "Kaden RPE1 (quoted)"]):
        v = summ.set_index("source").loc[q, "median_reliability"]
        ax.axhline(
            v, ls="--", lw=0.8, color=[style.GREEN, style.VERMILION][i], label=f"{q}: {v:.2f}"
        )
    ax.set_ylabel("per-target split-half reliability")
    ax.legend(fontsize=6)
    save(fig, "c4_P2_reliability")
    # P3 transfer delta
    cols = ["PDS", "MSE", "NMAE", "FID", "REACH", "JAC", "Overall"]
    rows = [
        {"fold": f, **(vcc[f].loc["C1a_KOLF", cols] - vcc[f].loc["C1a", cols]).to_dict()}
        for f in vcc
    ]
    ml = mean[mean.scope == "all_cells"]
    for h in ATLASES:
        a = ml[(ml.heldout == h) & (ml.arm == "C1a_KOLF")].iloc[0]
        b = ml[(ml.heldout == h) & (ml.arm == "C1a")].iloc[0]
        rows.append(
            {
                "fold": f"{h} (mean level)",
                "cosine": a.cosine - b.cosine,
                "sign_acc": a.sign_acc - b.sign_acc,
            }
        )
    delta = pd.DataFrame(rows)
    delta.to_csv(src / "c4_P3_transfer_delta.csv", index=False)
    fig, axes = plt.subplots(1, 2, figsize=(style.FULL_WIDTH, 2.6))
    d = delta[delta.fold.isin(["H1", "K562"])].set_index("fold")[cols]
    d.T.plot.bar(ax=axes[0])
    style.reference_line(axes[0], 0)
    axes[0].set_title("VCC Δ (C1 set + KOLF − C1 set)", fontsize=8)
    m = delta[delta.fold.str.contains("mean")].set_index("fold")[["cosine", "sign_acc"]]
    m.plot.bar(ax=axes[1])
    style.reference_line(axes[1], 0)
    axes[1].set_title("mean-level Δ", fontsize=8)
    axes[1].tick_params(labelsize=6)
    save(fig, "c4_P3_transfer_delta")
    # P4 context similarity (basal controls) + source agreement
    kolf = load_kolf()
    h1 = sources.load_source(c1f.SRC / "H1_2025_full_statistics.npz")
    k562 = sources.load_source(c1f.SRC / "K562_GWPS_CPM_full_statistics.npz")
    import anndata as ad

    profiles = {}
    for name, s in (("H1", h1), ("K562", k562), ("KOLF", kolf)):
        c = s.control_mean_cpm
        profiles[name] = pd.Series(c.mean(0) if c.ndim == 2 else c, index=s.genes)
    for ctx in "ABC":
        a = ad.read_h5ad(ROOT / f"data/raw/arc2026/controls/context_{ctx}.h5ad")
        x = a.X
        lib = np.asarray(x.sum(1)).ravel()
        cpm = np.asarray((x.multiply(1e6 / lib[:, None])).mean(0)).ravel()
        profiles[f"Arc {ctx}"] = pd.Series(cpm, index=a.var_names.astype(str))
    common = sorted(set.intersection(*[set(p.index) for p in profiles.values()]))
    mat = np.log1p(np.stack([profiles[k].loc[common].to_numpy(dtype=float) for k in profiles]))
    mat = mat - mat.mean(0)
    nrm = mat / np.linalg.norm(mat, axis=1, keepdims=True)
    sim = pd.DataFrame(nrm @ nrm.T, index=list(profiles), columns=list(profiles))
    sim.to_csv(src / "c4_P4_context_similarity.csv")
    fig, ax = plt.subplots(figsize=(style.HALF_WIDTH * 1.4, 3.0))
    ax.imshow(sim.to_numpy(), cmap="RdBu_r", vmin=-1, vmax=1)
    for i in range(len(sim)):
        for j in range(len(sim)):
            ax.text(j, i, f"{sim.iloc[i, j]:.2f}", ha="center", va="center", fontsize=6)
    ax.set_xticks(range(len(sim)), sim.columns, fontsize=6, rotation=45)
    ax.set_yticks(range(len(sim)), sim.index, fontsize=6)
    ax.set_title(f"basal control profiles (centred log1p CPM, {len(common)} genes)", fontsize=7)
    save(fig, "c4_P4_context_similarity")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage", choices=["stats", "reliability", "transfer", "decide"], required=True
    )
    parser.add_argument("--jobs", type=int, default=5)
    args = parser.parse_args()
    if hashlib.sha256(PREDECL.read_bytes()).hexdigest() != PREDECL_SHA:
        sys.exit("C4 predeclaration changed after it was frozen")
    OUT.mkdir(parents=True, exist_ok=True)
    {"stats": stage_stats, "reliability": stage_reliability, "decide": stage_decide}.get(
        args.stage, lambda: stage_transfer(args.jobs)
    )()


if __name__ == "__main__":
    main()
