"""C6 uncertainty-aware conserved response estimation (public folds only; no submission).

Predeclaration: ``reports/competition_v2/c6_predeclaration.md`` (SHA-256 checked).

    uv run python scripts/competition_v2/run_c6_uncertainty.py --stage moments
    uv run python scripts/competition_v2/run_c6_uncertainty.py --stage sigma
    uv run python scripts/competition_v2/run_c6_uncertainty.py --stage calibration
    uv run python scripts/competition_v2/run_c6_uncertainty.py --stage diagnose
    uv run python scripts/competition_v2/run_c6_uncertainty.py --stage vcc [--jobs 5]
    uv run python scripts/competition_v2/run_c6_uncertainty.py --stage spectrum
    uv run python scripts/competition_v2/run_c6_uncertainty.py --stage decide

Outputs: ``outputs/competition_v2/c6_uncertainty/``.
"""

from __future__ import annotations

import os

os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

# ruff: noqa: E402, E501

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sstats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import analyse_c2 as a2
import run_c1_public_folds as c1f
import run_c3_fusion as c3

from virtual_cell.arc import metrics as M
from virtual_cell.competition_v2 import evaluation, fusion_c3, generator, sources
from virtual_cell.competition_v2 import uncertainty as U
from virtual_cell.competition_v2.atlas import source_effect

OUT = ROOT / "outputs" / "competition_v2" / "c6_uncertainty"
PREDECL = ROOT / "reports" / "competition_v2" / "c6_predeclaration.md"
PREDECL_SHA = "c0a0ba8206aad13b3b869715672d2f4ab2da56b1f985e47089a32b61b9515cbe"
SRC = c1f.SRC
RAW = c1f.RAW
C2_FOLDS = ROOT / "outputs" / "competition_v2" / "c2_calibration" / "folds"
ATLASES = ["H1", "K562", "CD4"]
CANDIDATES = ["C1", "U1", "U2", "U3"]
HALF_SEED = generator.SEED + 61
CTRL = sources.CONTROL_LABEL
H1_FILES = {
    "train": "adata_Training.h5ad",
    "validation": "adata_Validation.h5ad",
    "test": "adata_Test.h5ad",
}


def log(m):
    print(m, flush=True)


def raw_stats(name: str) -> dict:
    f = {"K562": "K562_GWPS_CPM_full_statistics.npz", "H1": "H1_2025_full_statistics.npz"}[name]
    with np.load(SRC / f) as d:
        return {k: d[k] for k in d.files}


# ----------------------------------------------------------------------------- moments
def stage_moments() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    x, labels, genes = evaluation.load_cells(SRC / "fold_cells_K562.npz")
    groups = sorted(set(labels))
    ids = pd.Categorical(labels, categories=groups).codes.astype(np.int64)
    half = U.half_assignment(labels, [g for g in groups if g != CTRL], seed=HALF_SEED)
    m = U.Moments(groups, x.shape[1])
    for lo in range(0, x.shape[0], 8192):
        m.update(x[lo : lo + 8192], ids[lo : lo + 8192], half[lo : lo + 8192])
    np.savez_compressed(OUT / "moments_K562.npz", genes=genes, **m.result())
    log(f"K562 moments: {len(groups)} groups ({time.time() - t0:.0f}s)")
    import anndata as ad

    for split, fname in H1_FILES.items():
        t0 = time.time()
        a = ad.read_h5ad(RAW / fname, backed="r")
        labels = a.obs["target_gene"].astype(str).to_numpy()
        genes = a.var_names.astype(str).to_numpy(dtype=str)
        groups = sorted(set(labels))
        ids = pd.Categorical(labels, categories=groups).codes.astype(np.int64)
        half = U.half_assignment(labels, [g for g in groups if g != CTRL], seed=HALF_SEED)
        m = U.Moments(groups, len(genes))
        for lo in range(0, len(labels), 4096):
            hi = min(lo + 4096, len(labels))
            m.update(a.X[lo:hi], ids[lo:hi], half[lo:hi])
        a.file.close()
        np.savez_compressed(OUT / f"moments_H1_{split}.npz", genes=genes, **m.result())
        log(
            f"H1 {split} moments: {len(groups)} groups, {len(labels)} cells ({time.time() - t0:.0f}s)"
        )


def _load(p):
    with np.load(p, allow_pickle=False) as d:
        return {k: d[k] for k in d.files}


# ----------------------------------------------------------------------------- sigma2
def stage_sigma() -> None:
    out = {}
    # K562
    s = raw_stats("K562")
    k = sources.load_source(SRC / "K562_GWPS_CPM_full_statistics.npz")
    mo = _load(OUT / "moments_K562.npz")
    if not np.array_equal(mo["genes"].astype(str), k.genes.astype(str)):
        raise SystemExit("K562 cell genes differ from statistics genes")
    gi = {g: i for i, g in enumerate(mo["groups"].astype(str))}
    rows = np.array([gi.get(str(t), -1) for t in k.targets])
    c = gi[CTRL]
    v_t = np.where((rows >= 0)[:, None], mo["var"][np.maximum(rows, 0)], np.nan)
    sig = U.delta_sigma2(
        k.mean_cpm, k.control_mean_cpm, U.shrink_fraction(s["target_count_sums"])[:, None],
        v_t, k.n_cells[:, None], mo["var"][c][None, :], mo["n"][c],
    )  # fmt: skip
    out["K562"] = sig
    log(f"K562 sigma2: {np.isfinite(sig).mean():.3f} finite; median {np.nanmedian(sig):.4f}")
    # H1: each target against its own split's controls
    h = sources.load_source(SRC / "H1_2025_full_statistics.npz")
    hs = raw_stats("H1")
    split_of = dict(
        zip(hs["targets"].astype(str), hs["public_2025_split"].astype(str), strict=True)
    )
    moms = {sp: _load(OUT / f"moments_H1_{sp}.npz") for sp in H1_FILES}
    sig = np.full(h.mean_cpm.shape, np.nan)
    for i, t in enumerate(h.targets.astype(str)):
        sp = split_of.get(t) or next(
            sp for sp, mm in moms.items() if t in set(mm["groups"].astype(str))
        )
        mm = moms[sp]
        if not np.array_equal(mm["genes"].astype(str), h.genes.astype(str)):
            raise SystemExit("H1 cell genes differ from statistics genes")
        g = {x: j for j, x in enumerate(mm["groups"].astype(str))}
        ctrl_i = g[CTRL]
        sig[i] = U.delta_sigma2(
            h.mean_cpm[i], h.control_mean_cpm[i], U.shrink_fraction(hs["target_count_sums"][i : i + 1])[0],
            mm["var"][g[t]], h.n_cells[i], mm["var"][ctrl_i], mm["n"][ctrl_i],
        )  # fmt: skip
    out["H1"] = sig
    log(f"H1 sigma2: {np.isfinite(sig).mean():.3f} finite; median {np.nanmedian(sig):.4f}")
    cd4 = c1f.load_cd4()
    out["CD4"] = U.cd4_sigma2(cd4)
    log(
        f"CD4 sigma2: {np.isfinite(out['CD4']).mean():.3f} finite; median {np.nanmedian(out['CD4']):.4f}"
    )
    np.savez_compressed(OUT / "sigma2.npz", **out)


# ----------------------------------------------------------------------------- calibration
def _decile_table(sig, a, b, source):
    ok = np.isfinite(sig) & np.isfinite(a) & np.isfinite(b) & (a != 0) & (b != 0)
    sig, a, b = sig[ok], a[ok], b[ok]
    edges = np.quantile(sig, np.linspace(0, 1, 11))
    d = np.clip(np.searchsorted(edges, sig, side="right") - 1, 0, 9)
    rows = []
    for q in range(10):
        s = d == q
        rows.append(
            {
                "source": source,
                "decile": q + 1,
                "sigma2_median": float(np.median(sig[s])),
                "n": int(s.sum()),
                "sign_agreement": float(np.mean(np.sign(a[s]) == np.sign(b[s]))),
                "pearson": float(np.corrcoef(a[s], b[s])[0, 1]),
                "median_abs_diff": float(np.median(np.abs(a[s] - b[s]))),
            }
        )
    return rows


def stage_calibration() -> None:
    rows = []
    for name in ("K562", "H1"):
        if name == "K562":
            parts = [(_load(OUT / "moments_K562.npz"), None)]
        else:
            parts = [(_load(OUT / f"moments_H1_{sp}.npz"), sp) for sp in H1_FILES]
        A, B, S = [], [], []
        for mo, _ in parts:
            groups = mo["groups"].astype(str)
            c = int(np.flatnonzero(groups == CTRL)[0])
            ctrl, vc, nc = mo["mean"][c], mo["var"][c], mo["n"][c]
            keep = (mo["half_n"].min(axis=1) >= 20) & (groups != CTRL)
            hm, hv, hn = mo["half_mean"][keep], mo["half_var"][keep], mo["half_n"][keep]
            ea = np.log2((hm[:, 0] + 1) / (ctrl + 1))
            eb = np.log2((hm[:, 1] + 1) / (ctrl + 1))
            A.append(ea - ea.mean(0))
            B.append(eb - eb.mean(0))
            S.append(U.delta_sigma2(hm[:, 0], ctrl, 1.0, hv[:, 0], hn[:, 0:1], vc, nc))
        rows += _decile_table(
            np.concatenate([x.ravel() for x in S]),
            np.concatenate([x.ravel() for x in A]),
            np.concatenate([x.ravel() for x in B]),
            name,
        )
        log(f"calibration {name} done")
    cd4 = c1f.load_cd4()
    usable = cd4["available"] & cd4["quality_pass"] & (cd4["n_cells"] >= 20)
    lfc, se = cd4["log2fc"].astype(float), cd4["lfcSE"].astype(float)
    cent = []
    for ci in range(lfc.shape[0]):
        v = np.where(usable[ci][:, None], lfc[ci], np.nan)
        cent.append(lfc[ci] - np.nan_to_num(np.nanmean(v, axis=0)))
    A, B, S = [], [], []
    for c1i, c2i in ((0, 1), (0, 2), (1, 2)):
        both = usable[c1i] & usable[c2i]
        A.append(cent[c1i][both].ravel())
        B.append(cent[c2i][both].ravel())
        S.append((se[c1i][both] ** 2).ravel())
    rows += _decile_table(
        np.concatenate(S), np.concatenate(A), np.concatenate(B), "CD4 (proxy: condition pairs)"
    )
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "calibration_deciles.csv", index=False)
    verdict = {}
    for src, g in t.groupby("source"):
        rho = float(sstats.spearmanr(g.decile, g.sign_agreement).statistic)
        gap = float(g.sign_agreement.iloc[0] - g.sign_agreement.iloc[-1])
        verdict[src.split(" ")[0]] = {
            "spearman_decile_vs_sign_agreement": rho,
            "lowest_minus_highest": gap,
            "calibrated": bool(rho <= -0.8 and gap >= 0.10),
        }
    (OUT / "calibration_verdict.json").write_text(json.dumps(verdict, indent=2))
    log(json.dumps(verdict, indent=1))


# ----------------------------------------------------------------------------- fold machinery
def fold_inputs(atlas: str, targets, genes, controls=None, raw_sigma: bool = False):
    """Components, panel sigma2 (with calibration fallback), U2 tau2 and U3 multipliers."""
    green, cd4 = c1f.load_green(), c1f.load_cd4()
    srcs = {k: v for k, v in green.items() if k != atlas}
    use_cd4 = cd4 if atlas != "CD4" else None
    comp = fusion_c3.components(srcs, use_cd4, targets, genes)
    sig_src = _load(OUT / "sigma2.npz")
    verdict = json.loads((OUT / "calibration_verdict.json").read_text())
    sig, mult = {}, {}
    names = list(srcs) + (["CD4"] if use_cd4 is not None else [])
    for n in names:
        if n == "CD4":
            panel = U.to_panel(
                sig_src["CD4"], cd4["targets"].astype(str), cd4["genes"].astype(str), targets, genes
            )
            mask = comp["CD4"][1]
        else:
            s = srcs[n]
            panel = U.to_panel(
                sig_src[n],
                s.targets,
                s.genes,
                targets,
                genes,
                row_ok=s.n_cells >= 20,
                col_ok=s.measured,
            )
            mask = comp["log2fc"][n][1]
        panel = np.where(mask, panel, np.nan)
        if not verdict[n]["calibrated"] and not raw_sigma:
            panel = np.where(mask, np.nanmedian(panel[mask]), np.nan)
        sig[n] = panel
        # U3 multipliers from the source's own retained targets / genes; U3 shrinks only
        # sources whose sigma2 passed calibration (the others keep C1's effect, multiplier 1)
        if not verdict[n]["calibrated"]:
            mult[n] = np.ones((len(targets), len(genes)))
            continue
        if n == "CD4":
            eff_full, avail = c3.truth(
                "CD4", cd4["targets"].astype(str), cd4["genes"].astype(str), green, cd4
            )
            usable = avail.any(axis=1)
            mfull, _ = U.eb_multiplier(
                np.where(avail, eff_full, 0), np.where(avail, sig_src["CD4"], np.nan), usable
            )
            mult[n] = np.nan_to_num(
                U.to_panel(
                    mfull, cd4["targets"].astype(str), cd4["genes"].astype(str), targets, genes
                ),
                nan=1.0,
            )
        else:
            s = srcs[n]
            full, m = source_effect(s, s.targets, s.genes, space="log2fc")
            usable = s.n_cells >= 20
            sfull = np.where(m, sig_src[n], np.nan)
            mfull, _ = U.eb_multiplier(np.where(m, full, 0), sfull, usable)
            mult[n] = np.nan_to_num(U.to_panel(mfull, s.targets, s.genes, targets, genes), nan=1.0)
    effects = {
        n: (comp["CD4"][0] if n == "CD4" else comp["log2fc"][n][0]).astype(np.float64)
        for n in names
    }
    masks = {n: (comp["CD4"][1] if n == "CD4" else comp["log2fc"][n][1]) for n in names}
    tau2 = U.pooled_tau2(effects, sig, masks)
    weights = {}
    for n in names:
        med = np.nanmedian(1 / sig[n][masks[n]]) if masks[n].any() else 1.0
        w1 = np.where(np.isfinite(sig[n]) & (sig[n] > 0), 1 / np.where(sig[n] > 0, sig[n], 1), med)
        w2 = np.where(
            np.isfinite(sig[n]),
            1 / (np.nan_to_num(sig[n], nan=0) + tau2[None, :] + 1e-300),
            1 / (np.nanmedian(sig[n][masks[n]]) + np.median(tau2) + 1e-300),
        )
        weights[n] = {"U1": w1, "U2": w2}
    return comp, names, effects, masks, sig, weights, mult, tau2


def candidate(comp, names, weights, mult, which, space, ctrl_prob=None):
    if which == "C1":
        return fusion_c3.combine(comp, space, ctrl_prob)
    ones = {n: 1.0 for n in names}
    if which in ("U1", "U2"):
        return U.fuse(comp, space, ctrl_prob, {n: weights[n][which] for n in names})
    if which == "U3":
        return U.fuse(comp, space, ctrl_prob, ones, mult)
    raise ValueError(which)


# ----------------------------------------------------------------------------- diagnose
def _boot(correct_by_t, total_by_t, reps=2000, seed=0):
    rng = np.random.default_rng(seed)
    c, n = np.asarray(correct_by_t, float), np.asarray(total_by_t, float)
    keep = n > 0
    c, n = c[keep], n[keep]
    if not len(n):
        return np.nan, (np.nan, np.nan), 0
    est = c.sum() / n.sum()
    idx = rng.integers(0, len(n), size=(reps, len(n)))
    bs = c[idx].sum(1) / np.maximum(n[idx].sum(1), 1)
    return float(est), (float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))), int(n.sum())


def stage_diagnose() -> None:
    green, cd4 = c1f.load_green(), c1f.load_cd4()
    diag_rows, mean_rows, conflict_rows, ident_rows = [], [], [], []
    for atlas in ATLASES:
        t0 = time.time()
        targets, genes = c3.panel(atlas)
        t_eff, t_mask = c3.truth(atlas, targets, genes, green, cd4)
        comp, names, effects, masks, sig, weights, mult, tau2 = fold_inputs(atlas, targets, genes)
        excl = np.isin(genes, targets)
        valid = t_mask & ~excl[None, :]
        anym = np.zeros_like(valid)
        for n in names:
            anym |= masks[n]
        preds = {w: candidate(comp, names, weights, mult, w, "log2fc") for w in CANDIDATES}
        a, b = names
        cor_t, tot_t = np.zeros(len(targets)), np.zeros(len(targets))
        src_cor = {a: 0, b: 0}
        sig_right, sig_wrong = [], []
        conf = {w: [0, 0] for w in CANDIDATES}
        cls_counts = {"agree": 0, "conflict": 0, "single": 0}
        for i in range(len(targets)):
            v = np.flatnonzero(valid[i] & anym[i] & (preds["C1"][i] != 0))
            if not len(v):
                continue
            top = v[np.argsort(-np.abs(t_eff[i][v]))[:200]]
            ts = np.sign(t_eff[i, top])
            ea, eb = effects[a][i, top], effects[b][i, top]
            ma, mb = masks[a][i, top] & (ea != 0), masks[b][i, top] & (eb != 0)
            both = ma & mb
            agree = both & (np.sign(ea) == np.sign(eb))
            dis = both & (np.sign(ea) != np.sign(eb))
            cls_counts["agree"] += int(agree.sum())
            cls_counts["conflict"] += int(dis.sum())
            cls_counts["single"] += int((ma ^ mb).sum())
            if dis.any():
                sa, sb = sig[a][i, top][dis], sig[b][i, top][dis]
                pick_a = sa < sb
                lower_sign = np.where(pick_a, np.sign(ea[dis]), np.sign(eb[dis]))
                ok = np.isfinite(sa) & np.isfinite(sb) & (sa != sb)
                cor_t[i] = float((lower_sign[ok] == ts[dis][ok]).sum())
                tot_t[i] = float(ok.sum())
                src_cor[a] += int((np.sign(ea[dis]) == ts[dis]).sum())
                src_cor[b] += int((np.sign(eb[dis]) == ts[dis]).sum())
                right_a = np.sign(ea[dis]) == ts[dis]
                sig_right += list(np.where(right_a, sa, sb))
                sig_wrong += list(np.where(right_a, sb, sa))
                for w in CANDIDATES:
                    conf[w][0] += int((np.sign(preds[w][i, top][dis]) == ts[dis]).sum())
                    conf[w][1] += int(dis.sum())
        est, ci, n_pairs = _boot(cor_t, tot_t)
        diag_rows.append(
            {
                "heldout": atlas,
                "sources": f"{a} vs {b}",
                "p_lower_sigma_correct": est,
                "ci_low": ci[0],
                "ci_high": ci[1],
                "n_conflict_cells": n_pairs,
                "passes": bool(ci[0] > 0.5),
            }
        )
        nconf = conf["C1"][1]
        conflict_rows.append({
            "heldout": atlas, **{f"n_{k}": v for k, v in cls_counts.items()},
            **{f"acc_{w}": conf[w][0] / nconf if nconf else np.nan for w in CANDIDATES},
            f"acc_source_{a}": src_cor[a] / nconf if nconf else np.nan,
            f"acc_source_{b}": src_cor[b] / nconf if nconf else np.nan,
            "median_sigma2_correct_source": float(np.nanmedian(sig_right)) if sig_right else np.nan,
            "median_sigma2_wrong_source": float(np.nanmedian(sig_wrong)) if sig_wrong else np.nan,
        })  # fmt: skip
        vp = valid & anym
        for w in CANDIDATES:
            per, sse, sst = fusion_c3.row_metrics(preds[w], t_eff, vp)
            row = {"heldout": atlas, "candidate": w, **fusion_c3.summarise(per, sse, sst)}
            if atlas == "CD4":
                row["effect_pds"] = float(M.pds_cosine(preds[w], t_eff, exclude=excl).mean())
            mean_rows.append(row)
            cperc, _, _ = fusion_c3.row_metrics(preds[w], preds["C1"], vp)
            changed = np.abs(preds[w] - preds["C1"])[vp] > 0.1
            ident_rows.append({
                "heldout": atlas, "candidate": w,
                "norm_ratio_to_c1_median": float(np.nanmedian(np.linalg.norm(np.where(vp, preds[w], 0), axis=1) / np.maximum(np.linalg.norm(np.where(vp, preds["C1"], 0), axis=1), 1e-12))),
                "cosine_to_c1_mean": float(np.nanmean(cperc["cosine"])),
                "fraction_materially_changed": float(changed.mean()),
                "tau2_median": float(np.median(tau2)) if w == "U2" else np.nan,
            })  # fmt: skip
        log(
            f"[diagnose] {atlas}: P(lower-sigma correct)={est:.3f} CI {ci[0]:.3f}-{ci[1]:.3f} n={n_pairs} ({time.time() - t0:.0f}s)"
        )
    pd.DataFrame(diag_rows).to_csv(OUT / "cross_source_diagnostic.csv", index=False)
    pd.DataFrame(mean_rows).to_csv(OUT / "candidates_mean_level.csv", index=False)
    pd.DataFrame(conflict_rows).to_csv(OUT / "sign_conflicts.csv", index=False)
    pd.DataFrame(ident_rows).to_csv(OUT / "identity_protection.csv", index=False)
    gate = sum(r["passes"] for r in diag_rows) >= 2
    verdict = json.loads((OUT / "calibration_verdict.json").read_text())
    u3_ok = any(verdict[n]["calibrated"] for n in ("K562", "H1", "CD4"))  # shrinks calibrated only
    decision = {
        "weighting_branch_E_gate": gate,
        "u3_calibration_gate": u3_ok,
        "vcc_candidates": (["U1", "U2"] if gate else []) + (["U3"] if u3_ok else []),
    }
    (OUT / "branch_decision.json").write_text(json.dumps(decision, indent=2))
    log(json.dumps(decision))


def stage_rawdiag() -> None:
    """Secondary (non-gating): the §E diagnostic with each source's RAW sigma2."""
    green, cd4 = c1f.load_green(), c1f.load_cd4()
    rows = []
    for atlas in ATLASES:
        targets, genes = c3.panel(atlas)
        t_eff, t_mask = c3.truth(atlas, targets, genes, green, cd4)
        comp, names, effects, masks, sig, _, _, _ = fold_inputs(
            atlas, targets, genes, raw_sigma=True
        )
        valid = t_mask & ~np.isin(genes, targets)[None, :]
        anym = masks[names[0]] | masks[names[1]]
        c1 = fusion_c3.combine(comp, "log2fc", None)
        a, b = names
        cor, tot = np.zeros(len(targets)), np.zeros(len(targets))
        ratio_rows = []
        for i in range(len(targets)):
            v = np.flatnonzero(valid[i] & anym[i] & (c1[i] != 0))
            if not len(v):
                continue
            top = v[np.argsort(-np.abs(t_eff[i][v]))[:200]]
            ea, eb = effects[a][i, top], effects[b][i, top]
            dis = (
                masks[a][i, top]
                & masks[b][i, top]
                & (ea != 0)
                & (eb != 0)
                & (np.sign(ea) != np.sign(eb))
            )
            sa, sb = sig[a][i, top][dis], sig[b][i, top][dis]
            ok = np.isfinite(sa) & np.isfinite(sb) & (sa != sb)
            ts = np.sign(t_eff[i, top][dis])
            lower = np.where(sa < sb, np.sign(ea[dis]), np.sign(eb[dis]))
            cor[i], tot[i] = (lower[ok] == ts[ok]).sum(), ok.sum()
            # within-source: z = |e|/sigma of each source's call, correct vs wrong
            za, zb = np.abs(ea[dis]) / np.sqrt(sa), np.abs(eb[dis]) / np.sqrt(sb)
            higher_z = np.where(za > zb, np.sign(ea[dis]), np.sign(eb[dis]))
            okz = np.isfinite(za) & np.isfinite(zb) & (za != zb)
            ratio_rows.append(((higher_z[okz] == ts[okz]).sum(), okz.sum()))
        est, ci, n = _boot(cor, tot)
        cz = np.array([r[0] for r in ratio_rows], float)
        nz = np.array([r[1] for r in ratio_rows], float)
        estz, ciz, _ = _boot(cz, nz, seed=1)
        rows.append(
            {
                "heldout": atlas,
                "sources": f"{a} vs {b}",
                "p_lower_raw_sigma_correct": est,
                "ci_low": ci[0],
                "ci_high": ci[1],
                "p_higher_z_correct": estz,
                "z_ci_low": ciz[0],
                "z_ci_high": ciz[1],
                "n": n,
            }
        )
        log(
            f"[rawdiag] {atlas}: lower-raw-sigma {est:.3f} ({ci[0]:.3f}-{ci[1]:.3f}); higher |e|/sigma {estz:.3f} ({ciz[0]:.3f}-{ciz[1]:.3f}); n={n}"
        )
    pd.DataFrame(rows).to_csv(OUT / "cross_source_diagnostic_raw_sigma.csv", index=False)


# ----------------------------------------------------------------------------- vcc
def stage_vcc(jobs: int) -> None:
    arms = json.loads((OUT / "branch_decision.json").read_text())["vcc_candidates"]
    for name in ("H1", "K562"):
        t0 = time.time()
        fold = c1f.fold_h1(False) if name == "H1" else c1f.fold_k562(False)
        scorer = evaluation.FoldScorer(fold, jobs=jobs)
        comp, names, _, _, _, weights, mult, _ = fold_inputs(name, fold.targets, fold.genes)
        pairs = generator.promoter_pairs(
            generator.gencode_tss(RAW / "gencode.v47.annotation.gtf.gz"), fold.targets, fold.genes
        )
        rows = {}
        (rows["ANCHOR_split_half"], _), _ = scorer.split_half(generator.SEED + 2)
        real_cpm = np.stack([(b / b.sum(1, keepdims=True)).mean(0) for b in fold.real])
        real_bulk = np.stack([b.sum(0) / b.sum() for b in fold.real])
        n = len(fold.targets)
        moments = {
            "ANCHOR_mean_response": (
                np.tile(real_cpm.mean(0), (n, 1)),
                np.tile(real_bulk.mean(0), (n, 1)),
            )
        }
        for w in ["C1"] + arms:
            eff = {
                sp: candidate(comp, names, weights, mult, w, sp, fold.controls[sp])
                for sp in ("log2fc", "bulk_delta")
            }
            if w == "C1":
                eff = fusion_c3.fused(comp, fold.controls)
            moments[w] = generator.expected_moments(
                eff, fold.controls, fold.targets, fold.genes, pairs=pairs
            )
        for arm, (pc, pb) in moments.items():
            st = scorer.emit(arm, pc / pc.sum(1, keepdims=True), pb / pb.sum(1, keepdims=True))
            rows[arm], _ = scorer.members(st)
            log(
                f"[{name}] {arm:22s} PDS={rows[arm]['pds_cosine']:.4f} FID={rows[arm]['de_wilcoxon_direction_fidelity_yield_raw']:.4f} ({time.time() - t0:.0f}s)"
            )
        raw = pd.DataFrame(rows).T
        raw.index.name = "arm"
        (OUT / "folds" / name).mkdir(parents=True, exist_ok=True)
        raw.to_csv(OUT / "folds" / name / "scores_raw.csv")
        c2 = pd.read_csv(C2_FOLDS / name / "phase1_scores_raw.csv", index_col=0)
        diff = max(
            abs(float(raw.loc["C1", m]) - float(c2.loc["G0_a1.00", m])) for m in evaluation.MEMBERS
        )
        log(f"[{name}] C1 reproduces C2 G0_a1.00: {diff <= 1e-12} (max {diff:.2e})")
        if diff > 1e-12:
            raise SystemExit("C1 reproduction failed")


# ----------------------------------------------------------------------------- spectrum
def stage_spectrum() -> None:
    rows, verdict = [], {}
    for name in ("K562", "H1"):
        parts = (
            [_load(OUT / "moments_K562.npz")]
            if name == "K562"
            else [_load(OUT / f"moments_H1_{sp}.npz") for sp in H1_FILES]
        )
        S, N = [], []
        for mo in parts:
            groups = mo["groups"].astype(str)
            c = int(np.flatnonzero(groups == CTRL)[0])
            ctrl = mo["mean"][c]
            keep = (mo["half_n"].min(axis=1) >= 20) & (groups != CTRL)
            hm = mo["half_mean"][keep]
            ea = np.log2((hm[:, 0] + 1) / (ctrl + 1))
            eb = np.log2((hm[:, 1] + 1) / (ctrl + 1))
            S.append((ea + eb) / 2)
            N.append((ea - eb) / 2)
        Smat, Nmat = np.vstack(S), np.vstack(N)
        mo0 = parts[0]
        c0 = int(np.flatnonzero(mo0["groups"].astype(str) == CTRL)[0])
        ctrl_ok = mo0["mean"][c0] >= 5
        Smat, Nmat = Smat[:, ctrl_ok], Nmat[:, ctrl_ok]
        Smat -= Smat.mean(0)
        Nmat -= Nmat.mean(0)
        s = np.linalg.svd(Smat, compute_uv=False) ** 2
        z = np.linalg.svd(Nmat, compute_uv=False) ** 2
        r_star = int((s > 2 * z.max()).sum())
        signal_energy = max(s.sum() - z.sum(), 0)
        captured = (
            float((s[:r_star] - z[:r_star]).sum() / signal_energy)
            if signal_energy > 0 and r_star
            else 0.0
        )
        verdict[name] = {
            "r_star": r_star,
            "fraction_signal_energy_in_r_star": captured,
            "noise_top": float(z.max()),
            "signal_top": float(s.max()),
            "n_targets": int(Smat.shape[0]),
            "n_genes": int(Smat.shape[1]),
        }
        for i in range(min(60, len(s))):
            rows.append(
                {
                    "source": name,
                    "component": i + 1,
                    "signal_s2": float(s[i]),
                    "noise_s2": float(z[i]) if i < len(z) else np.nan,
                }
            )
    strong = all(
        v["r_star"] >= 3 and v["fraction_signal_energy_in_r_star"] >= 0.20 for v in verdict.values()
    )
    verdict["strong_separation"] = strong
    verdict["run_L1"] = strong
    pd.DataFrame(rows).to_csv(OUT / "spectrum.csv", index=False)
    (OUT / "spectrum_verdict.json").write_text(json.dumps(verdict, indent=2))
    log(json.dumps(verdict, indent=1))


# ----------------------------------------------------------------------------- decide
def stage_decide() -> None:
    mean = pd.read_csv(OUT / "candidates_mean_level.csv").set_index(["heldout", "candidate"])
    arms = json.loads((OUT / "branch_decision.json").read_text())["vcc_candidates"]
    vcc = {
        f: a2.scale(
            pd.read_csv(OUT / "folds" / f / "scores_raw.csv", index_col=0), "ANCHOR_mean_response"
        )
        for f in ("H1", "K562")
    }
    rows = [
        {"fold": f, "arm": a, **vcc[f].loc[a].to_dict()}
        for f in vcc
        for a in vcc[f].index
        if not a.startswith("ANCHOR")
    ]
    pd.DataFrame(rows).to_csv(OUT / "vcc_scaled.csv", index=False)
    verdicts = {}
    for w in arms:
        cos = {h: mean.loc[(h, w), "cosine"] - mean.loc[(h, "C1"), "cosine"] for h in ATLASES}
        sgn = {h: mean.loc[(h, w), "sign_acc"] - mean.loc[(h, "C1"), "sign_acc"] for h in ATLASES}
        gains = {f: float(vcc[f].loc[w, "Overall"] - vcc[f].loc["C1", "Overall"]) for f in vcc}
        total = sum(gains.values())
        pds = {f: float(vcc[f].loc[w, "PDS"] / vcc[f].loc["C1", "PDS"]) for f in vcc}
        pds["CD4"] = float(
            mean.loc[("CD4", w), "effect_pds"] / mean.loc[("CD4", "C1"), "effect_pds"]
        )
        member = {
            m: float(np.mean([vcc[f].loc[w, m] - vcc[f].loc["C1", m] for f in vcc]))
            for m in ("MSE", "NMAE", "FID", "REACH", "JAC")
        }
        crit = {
            "1_cosine_up": float(np.mean(list(cos.values()))) > 0,
            "2_sign_acc_up": float(np.mean(list(sgn.values()))) > 0,
            "3_both_up_in_2_of_3_folds": sum(cos[h] > 0 and sgn[h] > 0 for h in ATLASES) >= 2,
            "4_overall_up_by_0.005": float(np.mean(list(gains.values()))) >= 0.005,
            "5_pds_98pct": all(v >= 0.98 for v in pds.values()),
            "6_no_fold_over_75pct": bool(
                total > 0 and all(g / total <= 0.75 for g in gains.values())
            ),
            "7_no_member_drop_over_0.02": all(v >= -0.02 for v in member.values()),
        }
        verdicts[w] = {
            "pass": all(crit.values()),
            "criteria": crit,
            "cosine_delta": cos,
            "sign_delta": sgn,
            "overall_gain": gains,
            "pds_retention": pds,
            "member_delta": member,
        }
    best = (
        max(arms, key=lambda w: np.mean(list(verdicts[w]["overall_gain"].values())))
        if arms
        else None
    )
    decision = {
        "c6_pass": bool(best and verdicts[best]["pass"]),
        "best": best,
        "verdicts": verdicts,
    }
    (OUT / "c6_decision.json").write_text(json.dumps(decision, indent=2, default=float))
    print(
        json.dumps(
            {
                "c6_pass": decision["c6_pass"],
                "best": best,
                **{w: [k for k, v in verdicts[w]["criteria"].items() if not v] for w in verdicts},
            },
            indent=1,
        )
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--stage",
        required=True,
        choices=[
            "moments",
            "sigma",
            "calibration",
            "diagnose",
            "rawdiag",
            "vcc",
            "spectrum",
            "decide",
        ],
    )
    ap.add_argument("--jobs", type=int, default=5)
    args = ap.parse_args()
    if hashlib.sha256(PREDECL.read_bytes()).hexdigest() != PREDECL_SHA:
        sys.exit("C6 predeclaration changed after it was frozen")
    OUT.mkdir(parents=True, exist_ok=True)
    {
        "moments": stage_moments,
        "sigma": stage_sigma,
        "calibration": stage_calibration,
        "diagnose": stage_diagnose,
        "rawdiag": stage_rawdiag,
        "spectrum": stage_spectrum,
        "decide": stage_decide,
    }.get(args.stage, lambda: stage_vcc(args.jobs))()


if __name__ == "__main__":
    main()
