"""Public leave-one-atlas-out scoring of emitted counts with the local ``vcc2026`` metrics.

A fold holds one public context's real perturbed cells (truth), a disjoint split of its
control cells into a scoring reference and a template pool, and the panel it is scored
on. An arm is a pair of expected compositions per target. The arm is emitted through the
C0-style generator and scored exactly as ``scripts/competition_v2/run_atlasshift_public_h1.py``
scores it: the six members, the DE target index re-expressed on the tested axis, and
structure statistics.

Work is spread over targets with a ``fork`` process pool. Each worker emits one target's
cells and returns only per-target summaries, so no arm is ever resident as dense counts.
"""

from __future__ import annotations

import multiprocessing as mp
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import sparse

from virtual_cell.arc import metrics as M
from virtual_cell.competition_v2 import generator

#: Per-member direction (higher is better) for the six scored metrics.
MEMBERS = list(M.SCORED_METRICS)
DE_MEMBERS = [
    "de_wilcoxon_direction_fidelity_yield_raw",
    "de_wilcoxon_direction_reach_raw",
    "de_wilcoxon_sig_jaccard",
    "de_wilcoxon_lfc_nmae",
]


def load_cells(path) -> tuple[sparse.csr_matrix, np.ndarray, np.ndarray]:
    with np.load(path, allow_pickle=False) as d:
        x = sparse.csr_matrix((d["data"], d["indices"], d["indptr"]), shape=tuple(d["shape"]))
        return x, d["labels"].astype(str), d["genes"].astype(str)


@dataclass
class Fold:
    name: str
    targets: np.ndarray
    genes: np.ndarray
    real: list  # dense int32 blocks, one per target
    reference: np.ndarray  # dense int64 reference controls
    pool: sparse.csr_matrix  # template controls
    seed_prefix: str
    template: np.ndarray = field(init=False)
    depths: np.ndarray = field(init=False)
    ctrl_mean: np.ndarray = field(init=False)
    ctrl_bulk: np.ndarray = field(init=False)

    def __post_init__(self):
        self.template, self.depths, self.ctrl_mean, self.ctrl_bulk = generator.control_template(
            self.pool, 400, generator.POOL_K, generator.SEED
        )
        self.tg = np.array(
            [
                np.flatnonzero(self.genes == t)[0] if t in set(self.genes) else -1
                for t in self.targets
            ]
        )
        self.exclude = np.isin(self.genes, self.targets)
        self.ctrl_profile = M.bulk_profile(self.reference)
        self.ctrl_disp = M.jackknife_dispersion(self.reference)

    @property
    def controls(self) -> dict[str, np.ndarray]:
        return {"log2fc": self.ctrl_mean, "bulk_delta": self.ctrl_bulk}


# --------------------------------------------------------------------------- workers
_G: dict = {}


def _summaries(counts: np.ndarray, ref: np.ndarray, tested: np.ndarray, p_cpm, p_bulk):
    x = counts.astype(np.float64)
    lib = x.sum(axis=1)
    logcpm = np.log1p(1e4 * x / np.where(lib > 0, lib, 1)[:, None])
    de = M.de_table([counts], ref, tested=tested)
    out = {
        "profile": M.bulk_profile(counts),
        "disp": M.jackknife_dispersion(counts),
        "pval": de.pval[0],
        "p_adj": de.p_adj[0],
        "lfc": de.lfc[0],
        "lib": lib,
        "detected": (x > 0).sum(axis=1),
        "spread": np.linalg.norm(logcpm - logcpm.mean(axis=0, keepdims=True), axis=1),
        "nnz": int((x > 0).sum()),
        "size": int(x.size),
    }
    if p_cpm is not None:
        comp = (x / lib[:, None]).mean(axis=0)
        bulk = x.sum(axis=0) / x.sum()
        out["fid_cpm_l1"] = float(np.abs(comp - p_cpm).sum())
        out["fid_bulk_l1"] = float(np.abs(bulk - p_bulk).sum())
    return out


def _emit(i: int):
    g = _G
    t = g["targets"][i]
    counts = generator.dual_moment_counts(
        g["template"],
        g["p_cpm"][i],
        g["p_bulk"][i],
        depths=g["depths"],
        seed=generator.seed_for(g["prefix"] + t),
    )
    return _summaries(counts, g["ref"], g["tested"], g["p_cpm"][i], g["p_bulk"][i])


def _real(i: int):
    g = _G
    return _summaries(g["blocks"][i], g["ref"], g["tested"], None, None)


def _pool_map(fn, n: int, jobs: int, **state):
    _G.clear()
    _G.update(state)
    ctx = mp.get_context("fork")
    with ctx.Pool(jobs) as pool:
        return pool.map(fn, range(n), chunksize=1)


def _stack(results) -> dict:
    return {
        "profiles": np.stack([r["profile"] for r in results]),
        "disp": np.array([r["disp"] for r in results]),
        "de": M.DETable(
            tested=None,
            pval=np.vstack([r["pval"] for r in results]),
            p_adj=np.vstack([r["p_adj"] for r in results]),
            lfc=np.vstack([r["lfc"] for r in results]),
        ),
        "lib": np.concatenate([r["lib"] for r in results]),
        "detected": np.concatenate([r["detected"] for r in results]),
        "spread": np.concatenate([r["spread"] for r in results]),
        "density": sum(r["nnz"] for r in results) / sum(r["size"] for r in results),
        "nnz": sum(r["nnz"] for r in results),
        "fid_cpm_l1": np.array([r.get("fid_cpm_l1", np.nan) for r in results]),
        "fid_bulk_l1": np.array([r.get("fid_bulk_l1", np.nan) for r in results]),
    }


# --------------------------------------------------------------------------- scorer
class FoldScorer:
    """Truth-side tables for one fold, computed once; then any number of arms."""

    def __init__(self, fold: Fold, jobs: int = 9, log=print):
        self.fold, self.jobs, self.log = fold, jobs, log
        real = M.de_table([fold.real[0][:2]], fold.reference)  # tested mask only
        self.tested = real.tested
        res = _pool_map(
            _real, len(fold.real), jobs, blocks=fold.real, ref=fold.reference, tested=self.tested
        )
        st = _stack(res)
        self.real_de = M.DETable(
            tested=self.tested, pval=st["de"].pval, p_adj=st["de"].p_adj, lfc=st["de"].lfc
        )
        self.real_profiles = st["profiles"]
        self.real_disp = st["disp"]
        self.real_structure = st
        pos = np.cumsum(self.tested) - 1
        self.tg_de = np.array([pos[i] if i >= 0 and self.tested[i] else -1 for i in fold.tg])
        self.keep = ~fold.exclude
        self.true_norm = np.linalg.norm(
            (self.real_profiles - fold.ctrl_profile)[:, self.keep], axis=1
        )

    # -- metrics -------------------------------------------------------------------
    def members(self, st, *, real_de=None, cprof=None, cdisp=None, rprof=None, rdisp=None):
        f = self.fold
        real_de = self.real_de if real_de is None else real_de
        cprof = f.ctrl_profile if cprof is None else cprof
        cdisp = f.ctrl_disp if cdisp is None else cdisp
        rprof = self.real_profiles if rprof is None else rprof
        rdisp = self.real_disp if rdisp is None else rdisp
        de = M.DETable(
            tested=real_de.tested, pval=st["de"].pval, p_adj=st["de"].p_adj, lfc=st["de"].lfc
        )
        pds = M.pds_cosine(
            st["profiles"] - cprof[None, :], rprof - cprof[None, :], exclude=f.exclude
        )
        mse = M.expr_mse_unbiased_capped(
            st["profiles"],
            rprof,
            cprof,
            pred_dispersion=st["disp"],
            real_dispersion=rdisp,
            ctrl_dispersion=cdisp,
            target_gene=f.tg,
        )
        tg = self.tg_de
        per = {
            "pds": pds,
            "fid": M.direction_fidelity_yield(de, real_de, target_gene=tg),
            "reach": M.direction_reach(de, real_de, target_gene=tg),
            "jac": M.sig_jaccard(de, real_de, target_gene=tg),
            "nmae": M.lfc_nmae(de, real_de, target_gene=tg),
        }
        n = len(f.targets)
        n_pred, n_real, k = np.zeros(n), np.zeros(n), np.zeros(n)
        for p in range(n):
            rs = M._drop_target(real_de.significant[p], tg[p])
            ps = M._drop_target(de.significant[p] & real_de.adjudicable[p], tg[p])
            n_real[p], n_pred[p] = rs.sum(), ps.sum()
            k[p] = (ps & (np.sign(de.lfc[p]) == np.sign(real_de.lfc[p]))).sum()
        per.update(n_pred=n_pred, n_real=n_real, k_correct=k)
        raw = {
            "pds_cosine": float(np.mean(pds)),
            "expr_mse_unbiased_capped_norm": float(mse.value),
            "de_wilcoxon_direction_fidelity_yield_raw": float(np.nanmean(per["fid"])),
            "de_wilcoxon_direction_reach_raw": float(np.nanmean(per["reach"])),
            "de_wilcoxon_sig_jaccard": float(np.nanmean(per["jac"])),
            "de_wilcoxon_lfc_nmae": float(np.nanmean(per["nmae"])),
        }
        with np.errstate(invalid="ignore", divide="ignore"):
            raw["predicted_de_median"] = float(np.median(n_pred))
            raw["real_de_median"] = float(np.median(n_real))
            raw["directional_precision_pooled"] = float(k.sum() / max(n_pred.sum(), 1))
            raw["de_yield_median"] = float(np.nanmedian(np.minimum(1, n_pred / n_real)))
        return raw, per

    def emit(self, name: str, p_cpm: np.ndarray, p_bulk: np.ndarray):
        f = self.fold
        res = _pool_map(
            _emit,
            len(f.targets),
            self.jobs,
            targets=f.targets,
            template=f.template,
            depths=f.depths,
            p_cpm=p_cpm,
            p_bulk=p_bulk,
            prefix=f.seed_prefix,
            ref=f.reference,
            tested=self.tested,
        )
        return _stack(res)

    def structure(self, name: str, st: dict) -> dict:
        return {
            "group": name,
            "n_cells": int(st["lib"].size),
            "library_median": float(np.median(st["lib"])),
            "library_q10": float(np.percentile(st["lib"], 10)),
            "library_q90": float(np.percentile(st["lib"], 90)),
            "genes_detected_median": float(np.median(st["detected"])),
            "genes_detected_q10": float(np.percentile(st["detected"], 10)),
            "genes_detected_q90": float(np.percentile(st["detected"], 90)),
            "sparsity": float(1 - st["density"]),
            "heterogeneity_median": float(np.median(st["spread"])),
            "nnz": int(st["nnz"]),
            "fidelity_cpm_l1_median": float(np.nanmedian(st["fid_cpm_l1"]))
            if np.isfinite(st["fid_cpm_l1"]).any()
            else np.nan,
            "fidelity_bulk_l1_median": float(np.nanmedian(st["fid_bulk_l1"]))
            if np.isfinite(st["fid_bulk_l1"]).any()
            else np.nan,
        }

    def split_half(self, seed: int):
        """ORACLE local replicate anchor: half of each truth group vs the other half."""
        f = self.fold
        rng = np.random.default_rng(seed)
        rp = rng.permutation(len(f.reference))
        ca, cb = f.reference[rp[: len(rp) // 2]], f.reference[rp[len(rp) // 2 :]]
        ha, hb = [], []
        for b in f.real:
            i = rng.permutation(len(b))
            ha.append(b[i[: len(b) // 2]])
            hb.append(b[i[len(b) // 2 :]])
        truth_b = M.de_table(hb[:1], cb)
        tested_b = truth_b.tested
        res_b = _stack(_pool_map(_real, len(hb), self.jobs, blocks=hb, ref=cb, tested=tested_b))
        res_a = _stack(_pool_map(_real, len(ha), self.jobs, blocks=ha, ref=ca, tested=tested_b))
        real_de = M.DETable(
            tested=tested_b, pval=res_b["de"].pval, p_adj=res_b["de"].p_adj, lfc=res_b["de"].lfc
        )
        pos = np.cumsum(tested_b) - 1
        saved = self.tg_de
        self.tg_de = np.array([pos[i] if i >= 0 and tested_b[i] else -1 for i in f.tg])
        try:
            return self.members(
                res_a,
                real_de=real_de,
                cprof=M.bulk_profile(cb),
                cdisp=M.jackknife_dispersion(cb),
                rprof=res_b["profiles"],
                rdisp=res_b["disp"],
            ), res_a
        finally:
            self.tg_de = saved


def scale_local(raw: pd.DataFrame) -> pd.DataFrame:
    """s = (u - b)/(r - b) per member against the fold's local anchors; MSE clamped."""
    b = raw.loc["ANCHOR_mean_response"]
    r = raw.loc["ANCHOR_split_half"]
    out = pd.DataFrame(index=raw.index)
    for m in MEMBERS:
        span = r[m] - b[m]
        val = (raw[m] - b[m]) / span if span != 0 else np.nan
        if m == "expr_mse_unbiased_capped_norm":
            val = val.clip(0, 1)
        out[m] = val
    out["overall"] = out[MEMBERS].mean(axis=1)
    return out
