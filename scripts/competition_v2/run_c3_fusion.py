"""C3 public folds: mean-response diagnostics, source fusion candidates, VCC scoring.

Predeclaration: ``reports/competition_v2/c3_predeclaration.md`` (SHA-256 checked first).

    uv run python scripts/competition_v2/run_c3_fusion.py --stage mean
    uv run python scripts/competition_v2/run_c3_fusion.py --stage vcc --arms C1a R1 R2 Q
    uv run python scripts/competition_v2/run_c3_fusion.py --stage decide

* ``mean``: for H1 / K562 / CD4, the diagnostic object (per-source, fused and true
  responses, compact), source quality, magnitude decomposition, gene-wise calibration,
  sign consensus, reliability, nested inner-holdout weights and scales, and mean-level
  metrics plus effect PDS for every candidate. No cells are generated.
* ``vcc``: named candidates through the unchanged C1 emitter on H1 / K562, scored on the
  frozen local anchors.
* ``decide``: rule O over the capacity ladder, tables and figures.

Outputs: ``outputs/competition_v2/c3_fusion/``.
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

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats as sstats  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import run_c1_public_folds as c1f  # noqa: E402

from virtual_cell.arc import metrics as M  # noqa: E402
from virtual_cell.competition_v2 import (  # noqa: E402
    evaluation,
    fusion_c3,
    generator,
    licensing,
    sources,
)
from virtual_cell.competition_v2.atlas import cd4_effect  # noqa: E402

OUT = ROOT / "outputs" / "competition_v2" / "c3_fusion"
PREDECL = ROOT / "reports" / "competition_v2" / "c3_predeclaration.md"
PREDECL_SHA = "2aa25c728b3ffbd38b70057c74f4204d8b7f5da67601497a23ca46882c3ac9c8"
C2_FOLDS = ROOT / "outputs" / "competition_v2" / "c2_calibration" / "folds"
SRC = c1f.SRC
ATLASES = ["H1", "K562", "CD4"]
SOURCE_OF = {"H1": "H1", "K562": "K562", "CD4": "CD4"}
CONSENSUS = {"C3d_m": 0.75, "C3d_s": 0.5}
LADDER = {2: ["R1", "R2", "Q"], 3: ["S1", "S2"], 5: ["P", "C3d_m", "C3d_s"], 6: ["K"]}


def log(msg: str) -> None:
    print(msg, flush=True)


# ----------------------------------------------------------------------------- panels
def panel(atlas: str) -> tuple[np.ndarray, np.ndarray]:
    """The fold panel exactly as the C1 folds define it (targets, genes)."""
    official = c1f.official_genes()
    if atlas == "H1":
        with np.load(SRC / "fold_cells_H1.npz") as d:
            labels = d["labels"].astype(str)
        h1_genes = pd.read_csv(c1f.RAW / "gene_names.csv").iloc[:, 0].astype(str).to_numpy()
        return np.array(sorted(set(labels) - {sources.CONTROL_LABEL})), official[
            np.isin(official, h1_genes)
        ]
    if atlas == "K562":
        with np.load(SRC / "fold_cells_K562.npz") as d:
            labels, genes = d["labels"].astype(str), d["genes"].astype(str)
        stats = sources.load_source(SRC / "K562_GWPS_CPM_full_statistics.npz")
        present = set(labels)
        return np.array([t for t in stats.usable() if t in present]), official[
            np.isin(official, genes)
        ]
    cd4 = c1f.load_cd4()
    usable = cd4["available"] & cd4["quality_pass"] & (cd4["n_cells"] >= 20)
    return cd4["targets"].astype(str)[usable.any(axis=0)], official[
        np.isin(official, cd4["genes"].astype(str))
    ]


def truth(atlas, targets, genes, green, cd4):
    if atlas == "CD4":
        return cd4_effect(cd4, targets, genes)
    return c1f.heldout_truth_effect(atlas, targets, genes, green)


def predictors(atlas, green, cd4, exclude=()):
    srcs = {k: v for k, v in green.items() if k != atlas and k not in exclude}
    use_cd4 = cd4 if atlas != "CD4" and "CD4" not in exclude else None
    names = list(srcs) + (["CD4"] if use_cd4 is not None else [])
    licensing.assert_sources_allowed(names)
    return srcs, use_cd4


# ----------------------------------------------------------------------------- noise
def noise_fraction(atlas: str, targets, genes, cd4, reliab) -> np.ndarray:
    """Per target, the fraction of the truth's energy that is measurement noise."""
    out = np.full(len(targets), np.nan)
    if atlas in ("H1", "K562"):
        sh = reliab[atlas]
        for i, t in enumerate(targets):
            if t in sh["noise"] and sh["energy"][t] > 0:
                out[i] = sh["noise"][t] / sh["energy"][t]
        return np.clip(out, 0, 0.95)
    tg = cd4["targets"].astype(str)
    gi = pd.Index(cd4["genes"].astype(str)).get_indexer(genes)
    usable = cd4["available"] & cd4["quality_pass"] & (cd4["n_cells"] >= 20)
    eff, _ = cd4_effect(cd4, targets, genes)
    for i, t in enumerate(targets):
        k = int(np.flatnonzero(tg == t)[0])
        conds = np.flatnonzero(usable[:, k])
        if not len(conds):
            continue
        se = cd4["lfcSE"][conds][:, k][:, gi[gi >= 0]].astype(float)
        var = np.nansum(se**2, axis=0) / len(conds) ** 2
        energy = float((eff[i][gi >= 0].astype(float) ** 2).sum())
        if energy > 0:
            out[i] = float(np.nansum(var)) / energy
    return np.clip(out, 0, 0.95)


# ----------------------------------------------------------------------------- mean stage
def reliabilities() -> dict:
    out = {}
    for atlas in ("H1", "K562"):
        x, labels, _ = evaluation.load_cells(SRC / f"fold_cells_{atlas}.npz")
        targets = sorted(set(labels) - {sources.CONTROL_LABEL})
        out[atlas] = fusion_c3.split_half(
            x, labels, targets, sources.CONTROL_LABEL, seed=generator.SEED + 31
        )
        log(f"[reliability] {atlas}: {len(out[atlas]['rel'])} targets")
    cd4 = c1f.load_cd4()
    out["CD4"] = {"rel": fusion_c3.cd4_reliability(cd4)}
    return out


def single_metrics(eff, msk, t, valid):
    per, sse, sst = fusion_c3.row_metrics(np.where(msk, eff, 0), t, valid & msk)
    return per, fusion_c3.summarise(per, sse, sst)


def inner_quantities(outer, green, cd4, reliab) -> dict:
    """For each predictor source of ``outer``: transfer cosine and scale to the third atlas."""
    res = {}
    for s in [a for a in ATLASES if a != outer]:
        h = next(a for a in ATLASES if a not in (outer, s))
        targets, genes = panel(h)
        t_eff, t_mask = truth(h, targets, genes, green, cd4)
        srcs = {s: green[s]} if s != "CD4" else {}
        comp = fusion_c3.components(srcs, cd4 if s == "CD4" else None, targets, genes)
        e = fusion_c3.combine(comp, "log2fc", None)
        m = comp["log2fc"][s][1] if s != "CD4" else comp["CD4"][1]
        valid = t_mask & ~np.isin(genes, targets)[None, :] & m
        per, _, _ = fusion_c3.row_metrics(e, t_eff, valid)
        nf = noise_fraction(h, targets, genes, cd4, reliab)
        corrected = per["norm_ratio"] / np.sqrt(1 - np.nan_to_num(nf, nan=0.0))
        res[s] = {
            "inner_heldout": h,
            "q_transfer_cosine": float(np.nanmean(per["cosine"])),
            "a_scale": float(1 / np.nanmedian(corrected)),
            "n_targets": int(np.isfinite(per["cosine"]).sum()),
        }
    return res


def candidate_specs(outer, inner, reliab, targets, comp) -> dict:
    """Weights / scales / factors for every candidate on this outer fold."""
    names = [n for n, _, _ in fusion_c3.source_list(comp)]
    rel = {n: float(np.median(list(reliab[n]["rel"].values()))) for n in names}
    norm = lambda d: {k: v / sum(d.values()) * len(d) for k, v in d.items()}  # noqa: E731
    q = {n: max(inner[n]["q_transfer_cosine"], 0.0) for n in names}
    specs = {
        "C1a": {},
        "R1": {"weights": norm(rel)},
        "R2": {"weights": norm({k: v**0.5 for k, v in rel.items()})},
        "Q": {"weights": norm(q) if sum(q.values()) > 0 else {}},
        "S1": {"scales": {n: inner[n]["a_scale"] for n in names}},
        "S2": {"scales": {n: inner[n]["a_scale"] for n in names}, "weights": norm(rel)},
        "P": {
            "weights": {
                n: np.array([reliab[n]["rel"].get(str(t), rel[n]) for t in targets]) for n in names
            }
        },
    }
    classes = fusion_c3.consensus_classes(comp)
    for k, f in CONSENSUS.items():
        specs[k] = {"factor": fusion_c3.consensus_factor(classes, f)}
    meta = {
        "reliability_global": rel,
        "inner": inner,
        "weights": {k: v.get("weights") for k, v in specs.items() if k in ("R1", "R2", "Q", "S2")},
        "scales": {k: v.get("scales") for k, v in specs.items() if k in ("S1", "S2")},
    }
    return specs, meta


def stage_mean() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    green, cd4 = c1f.load_green(), c1f.load_cd4()
    reliab = reliabilities()
    rel_rows = [
        {"source": s, "target": t, "reliability": r}
        for s in ATLASES
        for t, r in reliab[s]["rel"].items()
    ]
    pd.DataFrame(rel_rows).to_csv(OUT / "source_reliability.csv", index=False)

    transfer, magnitude, cand_rows, consensus_rows, pert_rows = [], [], [], [], []
    gene_tables = {}
    for outer in ATLASES:
        t0 = time.time()
        targets, genes = panel(outer)
        t_eff, t_mask = truth(outer, targets, genes, green, cd4)
        srcs, use_cd4 = predictors(outer, green, cd4)
        comp = fusion_c3.components(srcs, use_cd4, targets, genes)
        exclude = np.isin(genes, targets)
        valid = t_mask & ~exclude[None, :]
        nf = noise_fraction(outer, targets, genes, cd4, reliab)
        fused = fusion_c3.combine(comp, "log2fc", None)
        pred_any = np.zeros_like(valid)
        single = {}
        # per-source quality vs the held-out truth (section D / E)
        for name, e, m in fusion_c3.source_list(comp):
            per, summ = single_metrics(e, m, t_eff, valid)
            single[name] = per
            pred_any |= m
            covered = np.flatnonzero(m[:, ~exclude].any(axis=1))
            pds = (
                float(
                    M.pds_cosine(np.where(m, e, 0)[covered], t_eff[covered], exclude=exclude).mean()
                )
                if len(covered) > 1
                else np.nan
            )
            gene_bias = np.nanmean(np.where(valid & m, e - t_eff, np.nan), axis=0)
            pert_bias = np.nanmean(np.where(valid & m, e - t_eff, np.nan), axis=1)
            top = []
            for i in covered:
                v = valid[i] & m[i]
                idx = np.flatnonzero(v)[np.argsort(-np.abs(t_eff[i][v]))[:200]]
                if len(idx):
                    top.append(
                        np.abs(e[i, idx] - t_eff[i, idx]).sum() / np.abs(t_eff[i, idx]).sum()
                    )
            transfer.append(
                {
                    "source": name,
                    "heldout": outer,
                    **summ,
                    "norm_ratio_noise_corrected": float(
                        np.nanmedian(per["norm_ratio"] / np.sqrt(1 - np.nan_to_num(nf)))
                    ),
                    "effect_pds": pds,
                    "nmae_top200": float(np.mean(top)) if top else np.nan,
                    "gene_bias_abs_mean": float(np.nanmean(np.abs(gene_bias))),
                    "pert_bias_mean": float(np.nanmean(pert_bias)),
                    "pert_bias_abs_mean": float(np.nanmean(np.abs(pert_bias))),
                    "reliability_global": float(np.median(list(reliab[name]["rel"].values()))),
                }
            )
        # per-target winner (perturbation dependence)
        names = list(single)
        both = np.isfinite(single[names[0]]["cosine"]) & np.isfinite(single[names[1]]["cosine"])
        diff = single[names[0]]["cosine"] - single[names[1]]["cosine"]
        for i in np.flatnonzero(both):
            pert_rows.append(
                {
                    "heldout": outer,
                    "target": targets[i],
                    f"cos_{names[0]}": single[names[0]]["cosine"][i],
                    f"cos_{names[1]}": single[names[1]]["cosine"][i],
                    "cos_diff": diff[i],
                    "rel_diff": reliab[names[0]]["rel"].get(str(targets[i]), np.nan)
                    - reliab[names[1]]["rel"].get(str(targets[i]), np.nan),
                }
            )
        valid_pred = valid & pred_any
        # magnitude decomposition (section S.1)
        per_f, sse, sst = fusion_c3.row_metrics(fused, t_eff, valid_pred)
        single_norms = np.nanmean(np.stack([single[n]["norm_ratio"] for n in names]), axis=0)
        magnitude.append(
            {
                "heldout": outer,
                "single_source_norm_ratio_median": float(np.nanmedian(single_norms)),
                "fused_norm_ratio_median": float(np.nanmedian(per_f["norm_ratio"])),
                "fused_over_single_mean_median": float(
                    np.nanmedian(per_f["norm_ratio"][both] / single_norms[both])
                ),
                "fused_norm_ratio_noise_corrected": float(
                    np.nanmedian(per_f["norm_ratio"] / np.sqrt(1 - np.nan_to_num(nf)))
                ),
                "after_emission_amplitude_0.6": float(0.6 * np.nanmedian(per_f["norm_ratio"])),
                "truth_noise_fraction_median": float(np.nanmedian(nf)),
                "two_source_targets": int(both.sum()),
                "cosine_between_sources_median": float(
                    np.nanmedian(
                        [
                            fusion_c3.row_metrics(comp_row_a[None], comp_row_b[None], vm[None])[0][
                                "cosine"
                            ][0]
                            for comp_row_a, comp_row_b, vm in _pairs(comp, valid)
                        ]
                    )
                ),
            }
        )
        # gene-wise calibration (section J)
        with np.errstate(invalid="ignore", divide="ignore"):
            vp = np.where(valid_pred, fused, np.nan)
            vt = np.where(valid_pred, t_eff, np.nan)
            gene_tables[outer] = pd.DataFrame(
                {
                    "gene": genes,
                    "mean_pred": np.nanmean(vp, axis=0),
                    "mean_true": np.nanmean(vt, axis=0),
                    "bias": np.nanmean(vt - vp, axis=0),
                    "sign_acc": np.nanmean(
                        np.where(
                            valid_pred & (fused != 0), np.sign(fused) == np.sign(t_eff), np.nan
                        ),
                        axis=0,
                    ),
                    "scale_ratio": np.nansum(vp * vt, axis=0) / np.nansum(vp * vp, axis=0),
                    "var_true": np.nanvar(vt, axis=0),
                    "n_obs": valid_pred.sum(axis=0),
                }
            )
        # sign consensus (section L)
        classes = fusion_c3.consensus_classes(comp)
        hits = {c: [0, 0] for c in (1, 2, 3)}
        for i in range(len(targets)):
            v = np.flatnonzero(valid_pred[i] & (fused[i] != 0))
            if not len(v):
                continue
            top = v[np.argsort(-np.abs(t_eff[i][v]))[:200]]
            for c in (1, 2, 3):
                sel = top[classes[i, top] == c]
                hits[c][0] += int((np.sign(fused[i, sel]) == np.sign(t_eff[i, sel])).sum())
                hits[c][1] += len(sel)
        for c, lab in [(1, "single"), (2, "agree"), (3, "conflict")]:
            consensus_rows.append(
                {
                    "heldout": outer,
                    "class": lab,
                    "n_top200_cells": hits[c][1],
                    "sign_acc": hits[c][0] / hits[c][1] if hits[c][1] else np.nan,
                    "fraction_of_top200": hits[c][1] / max(sum(h[1] for h in hits.values()), 1),
                }
            )
        # candidates (sections F-M, N)
        inner = inner_quantities(outer, green, cd4, reliab)
        specs, meta = candidate_specs(outer, inner, reliab, targets, comp)
        (OUT / f"specs_{outer}.json").write_text(json.dumps(meta, indent=2, default=float))
        for name, spec in specs.items():
            pred = fusion_c3.combine(comp, "log2fc", None, **spec)
            per, sse, sst = fusion_c3.row_metrics(pred, t_eff, valid_pred)
            row = {"heldout": outer, "candidate": name, **fusion_c3.summarise(per, sse, sst)}
            row["norm_ratio_noise_corrected"] = float(
                np.nanmedian(per["norm_ratio"] / np.sqrt(1 - np.nan_to_num(nf)))
            )
            if outer == "CD4":
                covered = np.flatnonzero(np.abs(pred[:, ~exclude]).sum(1) > 0)
                row["effect_pds"] = float(M.pds_cosine(pred, t_eff, exclude=exclude).mean())
                row["effect_pds_covered"] = int(len(covered))
            cand_rows.append(row)
        np.savez_compressed(
            OUT / f"diagnostic_{outer}.npz",
            targets=targets,
            genes=genes,
            truth=t_eff.astype(np.float16),
            truth_mask=t_mask,
            fused_c1a=fused.astype(np.float16),
            classes=classes,
            noise_fraction=nf,
            **{f"source_{n}": e.astype(np.float16) for n, e, _ in fusion_c3.source_list(comp)},
            **{f"mask_{n}": m for n, _, m in fusion_c3.source_list(comp)},
        )
        log(f"[mean] {outer}: {len(targets)} targets, {len(genes)} genes ({time.time() - t0:.0f}s)")

    pd.DataFrame(transfer).to_csv(OUT / "source_transfer_matrix.csv", index=False)
    pd.DataFrame(magnitude).to_csv(OUT / "magnitude_decomposition.csv", index=False)
    pd.DataFrame(consensus_rows).to_csv(OUT / "sign_consensus.csv", index=False)
    pd.DataFrame(pert_rows).to_csv(OUT / "per_target_source_cosine.csv", index=False)
    pd.DataFrame(cand_rows).to_csv(OUT / "candidates_mean_level.csv", index=False)
    for k, g in gene_tables.items():
        g.to_csv(OUT / f"gene_calibration_{k}.csv", index=False)
    triggers(gene_tables, pd.DataFrame(consensus_rows))


def _pairs(comp, valid):
    srcs = fusion_c3.source_list(comp)
    (_, ea, ma), (_, eb, mb) = srcs[0], srcs[1]
    for i in range(len(ea)):
        vm = valid[i] & ma[i] & mb[i]
        if vm.any():
            yield ea[i], eb[i], vm


def triggers(gene_tables: dict, consensus: pd.DataFrame) -> dict:
    """The two predeclared conditional triggers (§J → K, §L → C3d)."""
    rhos = {}
    for a, b in [("H1", "K562"), ("H1", "CD4"), ("K562", "CD4")]:
        ga, gb = gene_tables[a].set_index("gene"), gene_tables[b].set_index("gene")
        common = ga.index.intersection(gb.index)
        x, y = ga.loc[common, "bias"], gb.loc[common, "bias"]
        ok = x.notna() & y.notna()
        rhos[f"{a}~{b}"] = {
            "spearman": float(sstats.spearmanr(x[ok], y[ok]).statistic),
            "n_genes": int(ok.sum()),
        }
    j = all(v["spearman"] >= 0.3 for v in rhos.values())
    piv = consensus.pivot(index="heldout", columns="class", values="sign_acc")
    gap = (piv["agree"] - piv["conflict"]).to_dict()
    ell = all(np.isfinite(v) and v >= 0.05 for v in gap.values())
    out = {
        "J_gene_bias_spearman": rhos,
        "J_triggers_K": j,
        "L_agree_minus_conflict": gap,
        "L_triggers_C3d": ell,
    }
    (OUT / "triggers.json").write_text(json.dumps(out, indent=2))
    log(json.dumps(out, indent=1))
    return out


# ----------------------------------------------------------------------------- vcc stage
def stage_vcc(arms: list[str], folds: list[str], jobs: int) -> None:
    trig = json.loads((OUT / "triggers.json").read_text())
    green, cd4 = c1f.load_green(), c1f.load_cd4()
    reliab = None
    for name in folds:
        t0 = time.time()
        fold = c1f.fold_h1(False) if name == "H1" else c1f.fold_k562(False)
        targets, genes = panel(name)
        if not (np.array_equal(targets, fold.targets) and np.array_equal(genes, fold.genes)):
            raise SystemExit(f"[{name}] mean-stage panel differs from the VCC fold")
        scorer = evaluation.FoldScorer(fold, jobs=jobs)
        srcs, use_cd4 = predictors(name, green, cd4)
        comp = fusion_c3.components(srcs, use_cd4, fold.targets, fold.genes)
        if reliab is None:
            reliab = reliabilities()
        inner = json.loads((OUT / f"specs_{name}.json").read_text())["inner"]
        specs, _ = candidate_specs(name, inner, reliab, fold.targets, comp)
        pairs = generator.promoter_pairs(
            generator.gencode_tss(c1f.RAW / "gencode.v47.annotation.gtf.gz"),
            fold.targets,
            fold.genes,
        )
        out = OUT / "folds" / name
        out.mkdir(parents=True, exist_ok=True)
        raw_path = out / "scores_raw.csv"
        raw = pd.read_csv(raw_path, index_col=0) if raw_path.exists() else pd.DataFrame()
        todo = [a for a in arms if a not in raw.index]
        if "ANCHOR_mean_response" not in raw.index:
            todo = ["ANCHOR_mean_response", "ANCHOR_split_half"] + todo
        rows, per_rows = {}, []
        for arm in todo:
            if arm in ("C3d_m", "C3d_s") and not trig["L_triggers_C3d"]:
                log(f"[{name}] {arm} skipped: §L trigger not met")
                continue
            if arm == "K":
                log(f"[{name}] K skipped: {'not ' if not trig['J_triggers_K'] else ''}triggered")
                continue
            if arm == "ANCHOR_split_half":
                (r, _), _ = scorer.split_half(generator.SEED + 2)
                rows[arm] = r
                continue
            if arm == "ANCHOR_mean_response":
                real_cpm = np.stack([(b / b.sum(1, keepdims=True)).mean(0) for b in fold.real])
                real_bulk = np.stack([b.sum(0) / b.sum() for b in fold.real])
                n = len(fold.targets)
                p_cpm = np.tile(real_cpm.mean(0), (n, 1))
                p_bulk = np.tile(real_bulk.mean(0), (n, 1))
            else:
                eff = fusion_c3.fused(comp, fold.controls, **specs[arm])
                p_cpm, p_bulk = generator.expected_moments(
                    eff, fold.controls, fold.targets, fold.genes, pairs=pairs
                )
            p_cpm = p_cpm / p_cpm.sum(1, keepdims=True)
            p_bulk = p_bulk / p_bulk.sum(1, keepdims=True)
            st = scorer.emit(arm, p_cpm, p_bulk)
            r, per = scorer.members(st)
            rows[arm] = r
            frame = pd.DataFrame({k: v for k, v in per.items()})
            frame.insert(0, "target", fold.targets)
            frame.insert(1, "arm", arm)
            per_rows.append(frame)
            log(
                f"[{name}] {arm:22s} PDS={r['pds_cosine']:.4f} "
                f"MSE={r['expr_mse_unbiased_capped_norm']:.4f} "
                f"FID={r['de_wilcoxon_direction_fidelity_yield_raw']:.4f} "
                f"nDE={r['predicted_de_median']:.0f} ({time.time() - t0:.0f}s)"
            )
        new = pd.DataFrame(rows).T
        raw = pd.concat([raw, new]) if len(raw) else new
        raw.index.name = "arm"
        raw.to_csv(raw_path)
        if per_rows:
            per_path = out / "per_perturbation.csv"
            old = pd.read_csv(per_path) if per_path.exists() else None
            pd.concat(([old] if old is not None else []) + per_rows).to_csv(per_path, index=False)
        if "C1a" in raw.index:
            c2 = pd.read_csv(C2_FOLDS / name / "phase1_scores_raw.csv", index_col=0)
            diff = {
                m: abs(float(raw.loc["C1a", m]) - float(c2.loc["G0_a1.00", m]))
                for m in evaluation.MEMBERS
            }
            ok = max(diff.values()) <= 1e-12
            (out / "reproduction.json").write_text(
                json.dumps({"vs": "C2 G0_a1.00", "diffs": diff, "pass": ok}, indent=2)
            )
            log(f"[{name}] C1a reproduces C2 G0_a1.00: {ok}")
            if not ok:
                raise SystemExit("C1a reproduction failed")


# ----------------------------------------------------------------------------- decide
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["mean", "vcc"], required=True)
    parser.add_argument("--arms", nargs="+", default=["C1a"])
    parser.add_argument("--folds", nargs="+", default=["H1", "K562"])
    parser.add_argument("--jobs", type=int, default=5)
    args = parser.parse_args()
    if hashlib.sha256(PREDECL.read_bytes()).hexdigest() != PREDECL_SHA:
        sys.exit("C3 predeclaration changed after it was frozen")
    if args.stage == "mean":
        stage_mean()
    else:
        stage_vcc(args.arms, args.folds, args.jobs)


if __name__ == "__main__":
    main()
