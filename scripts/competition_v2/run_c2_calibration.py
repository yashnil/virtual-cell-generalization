"""C2 public folds: mean-response vs generator anatomy, null qualification, amplitude.

Predeclaration: ``reports/competition_v2/c2_predeclaration.md`` (SHA-256 checked first).
The C1 mean predictor, folds, reference/pool splits and seeds are reused unchanged from
``run_c1_public_folds.py``. Nothing here touches the Arc bundle or submits anything.

    uv run python scripts/competition_v2/run_c2_calibration.py --phase 1 [--folds H1 K562]
    uv run python scripts/competition_v2/analyse_c2.py --stage generator
    uv run python scripts/competition_v2/run_c2_calibration.py --phase 2
    uv run python scripts/competition_v2/analyse_c2.py --stage amplitude
    uv run python scripts/competition_v2/run_c2_calibration.py --phase 3   # only if F runs

Phase 1: G0 (C1 reproduction), G1c, G1b, G2, G3 (if justified) at a = 1, their null
arms, both anchors, split half, idealised pseudobulk. Phase 2: the amplitude grid under
G0 and under the selected generator (plus the idealised grid). Phase 3: per-target
amplitude. Outputs: ``outputs/competition_v2/c2_calibration/folds/<fold>/phase<k>_*``.
"""

from __future__ import annotations

import os

os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import argparse  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import run_c1_public_folds as c1f  # noqa: E402

from virtual_cell.arc import metrics as M  # noqa: E402
from virtual_cell.competition_v2 import evaluation, generator, realisation  # noqa: E402

OUT = ROOT / "outputs" / "competition_v2" / "c2_calibration"
PREDECL = ROOT / "reports" / "competition_v2" / "c2_predeclaration.md"
PREDECL_SHA = "dcb0df82849b7d9080002ce369f26e4d2d74c3edea827045b9060734abfd6cd1"
AMENDMENT = ROOT / "reports" / "competition_v2" / "c2_predeclaration_amendment_1.md"
AMENDMENT_SHA = "02c360b96e5e364f665baf70c413a8f52a83619b8f047a6d8c0a77559a5d596b"
AMENDMENT_2 = ROOT / "reports" / "competition_v2" / "c2_predeclaration_amendment_2.md"
AMENDMENT_2_SHA = "2c32a248b5af02b867b064bd05e4d919c6d3b46a2f5f8f2c94b2b0aeb62ddd05"
C1_FOLDS = ROOT / "outputs" / "competition_v2" / "c1_license_clean" / "folds"
AMPLITUDES = [0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00]
GENERATORS = ["G0", "G1c", "G1ci", "G1b", "G2", "G3"]
BETAS = [-0.5, 0.5]
N_DONORS = 400
N_PHI_CELLS = 3000
G3_VARIANCE_GATE = 0.8

_G: dict = {}


# --------------------------------------------------------------------------- workers
def _counts(i: int) -> np.ndarray:
    g = _G
    t = g["targets"][i]
    kind = g["gen"]
    if kind == "G0":
        return generator.dual_moment_counts(
            g["template"],
            g["p_cpm"][i],
            g["p_bulk"][i],
            depths=g["depths"],
            seed=generator.seed_for(g["prefix"] + t),
        )
    rng = np.random.default_rng(generator.seed_for(f"{g['prefix']}{kind}:{t}"))
    if kind == "G1ci":
        drng = np.random.default_rng(generator.seed_for(f"{g['prefix']}donors:{t}"))
        pick = np.sort(drng.choice(g["pool"].shape[0], N_DONORS, replace=False))
        ratio = realisation.effect_ratio(g["p_cpm"][i], g["c_cpm"])
        return realisation.g1_counts(g["pool"][pick].toarray(), ratio, g["c_cpm"], rng=rng)
    if kind == "G1c":
        ratio = realisation.effect_ratio(g["p_cpm"][i], g["c_cpm"])
        return realisation.g1_counts(g["donors"], ratio, g["c_cpm"], rng=rng)
    if kind == "G1b":
        ratio = realisation.effect_ratio(g["p_bulk"][i], g["c_bulk"])
        return realisation.g1_counts(g["donors"], ratio, g["c_bulk"], rng=rng)
    if kind == "G2":
        return realisation.g2_counts(g["donor_depths"], g["p_bulk"][i], rng=rng)
    if kind == "G3":
        return realisation.g3_counts(g["donor_depths"], g["p_bulk"][i], g["phi"], rng=rng)
    raise ValueError(kind)


def _ratio_median(a: np.ndarray, b: np.ndarray) -> float:
    ok = (b > 0) & np.isfinite(a)
    return float(np.median(a[ok] / b[ok])) if ok.any() else np.nan


def _compare(diag: dict, real: dict, tested: np.ndarray) -> dict:
    gv, rv = diag["gene_var"][tested], real["gene_var"][tested]
    gm, rm = diag["gene_mean"][tested], real["gene_mean"][tested]
    with np.errstate(invalid="ignore", divide="ignore"):
        gf = np.where(gm > 0, gv / gm, np.nan)
        rf = np.where(rm > 0, rv / rm, np.nan)
    return {
        "var_ratio": _ratio_median(gv, rv),
        "fano_ratio": _ratio_median(gf, rf),
        "logcpm_var_ratio": _ratio_median(diag["logcpm_var"][tested], real["logcpm_var"][tested]),
        "zero_fraction": diag["zero_fraction"],
        "cosine_distance": diag["cosine_distance_median"],
        "real_zero_fraction": real["zero_fraction"],
        "real_cosine_distance": real["cosine_distance_median"],
    }


def _emit(i: int):
    g = _G
    counts = _counts(i)
    out = evaluation._summaries(counts, g["ref"], g["tested"], g["p_cpm"][i], g["p_bulk"][i])
    diag = realisation.cell_diagnostics(counts, rng=np.random.default_rng(i))
    real = g["compare"][i] if isinstance(g["compare"], list) else g["compare"]
    out["cmp"] = _compare(diag, real, g["tested"])
    return out


def _real_diag(i: int):
    return realisation.cell_diagnostics(_G["blocks"][i], rng=np.random.default_rng(i))


def _pool(fn, n, jobs, **state):
    _G.clear()
    _G.update(state)
    with mp.get_context("fork").Pool(jobs) as pool:
        return pool.map(fn, range(n), chunksize=1)


# --------------------------------------------------------------------------- fold state
class C2Fold:
    def __init__(self, name: str, jobs: int, log, smoke: bool = False):
        self.name, self.jobs, self.log, self.smoke = name, jobs, log, smoke
        t0 = time.time()
        fold = c1f.fold_h1(smoke) if name == "H1" else c1f.fold_k562(smoke)
        self.fold = fold
        self.scorer = evaluation.FoldScorer(fold, jobs=jobs)
        self.tested = self.scorer.tested
        green, cd4 = c1f.load_green(), c1f.load_cd4()
        self.effects, self.names = c1f.arm_effects(
            name, fold.targets, fold.genes, fold.controls, green, cd4, None, "C1a"
        )
        cov = c1f.coverage(self.names, green, cd4, fold.targets, fold.genes)
        self.k_sources = cov.sum(1).to_numpy()
        tss = c1f.generator.gencode_tss(c1f.RAW / "gencode.v47.annotation.gtf.gz")
        self.pairs = generator.promoter_pairs(tss, fold.targets, fold.genes)
        rng = np.random.default_rng(generator.SEED + 7)
        pick = np.sort(rng.choice(fold.pool.shape[0], N_DONORS, replace=False))
        self.donors = fold.pool[pick].toarray().astype(np.int64)
        self.donor_depths = self.donors.sum(axis=1)
        phi_pick = np.sort(
            rng.choice(fold.pool.shape[0], min(N_PHI_CELLS, fold.pool.shape[0]), replace=False)
        )
        self.phi = realisation.fit_dispersion(
            fold.pool[phi_pick].toarray(), fold.ctrl_bulk / fold.ctrl_bulk.sum()
        )
        tested = self.tested
        self.real_diag = [
            {k: v for k, v in d.items()}
            for d in _pool(_real_diag, len(fold.real), jobs, blocks=fold.real)
        ]
        self.ref_diag = realisation.cell_diagnostics(fold.reference, rng=np.random.default_rng(0))
        n = len(fold.targets)
        real_cpm = np.stack([(b / b.sum(1, keepdims=True)).mean(0) for b in fold.real])
        real_bulk = np.stack([b.sum(0) / b.sum() for b in fold.real])
        self.anchor = (
            np.tile(real_cpm.mean(0), (n, 1)),
            np.tile(real_bulk.mean(0), (n, 1)),
        )
        zero = {k: np.zeros_like(v) for k, v in self.effects.items()}
        self.null = generator.expected_moments(zero, fold.controls, fold.targets, fold.genes)
        self.rows, self.per, self.struct, self.nulls = {}, [], [], []
        log(
            f"[{name}] ready: {n} targets, {len(fold.genes)} genes, {int(tested.sum())} tested, "
            f"pool {fold.pool.shape[0]}, predictors {self.names}+CD4 ({time.time() - t0:.0f}s)"
        )

    # ------------------------------------------------------------------ moments
    def moments(self, shrink=None):
        f = self.fold
        return generator.expected_moments(
            self.effects, f.controls, f.targets, f.genes, shrink=shrink, pairs=self.pairs
        )

    # ------------------------------------------------------------------ scoring
    def _record(self, arm: str, gen: str, st: dict, *, is_null=False, ideal=False):
        f, sc = self.fold, self.scorer
        raw, per = sc.members(st)
        tg = sc.tg_de
        sig = st["de"].p_adj < M.P_ADJ_THRESHOLD
        overlap = np.array(
            [
                (
                    M._drop_target(sig[p], tg[p]) & M._drop_target(sc.real_de.significant[p], tg[p])
                ).sum()
                for p in range(len(f.targets))
            ]
        )
        with np.errstate(invalid="ignore", divide="ignore"):
            recall = per["k_correct"] / per["n_real"]
            precision = per["k_correct"] / per["n_pred"]
        pred_norm = np.linalg.norm((st["profiles"] - f.ctrl_profile)[:, sc.keep], axis=1)
        raw.update(
            directional_recall_median=float(np.nanmedian(recall)),
            directional_precision_median=float(np.nanmedian(precision))
            if np.isfinite(precision).any()
            else np.nan,
            de_overlap_median=float(np.median(overlap)),
            effect_norm_ratio_median=float(np.median(pred_norm / sc.true_norm)),
            pred_effect_norm_median=float(np.median(pred_norm)),
            true_effect_norm_median=float(np.median(sc.true_norm)),
        )
        if not ideal:
            real_lib = sc.real_structure["lib"]
            real_det = sc.real_structure["detected"]
            cmp = pd.DataFrame(st["cmp"])
            raw.update(
                library_bias=float(np.median(st["lib"]) / np.median(real_lib) - 1),
                genes_detected_bias=float(np.median(st["detected"]) / np.median(real_det) - 1),
                pseudobulk_dispersion_ratio=float(np.median(st["disp"] / sc.real_disp)),
                **{f"cells_{k}": float(cmp[k].median()) for k in cmp.columns},
            )
        self.rows[arm] = {"generator": gen, **raw}
        frame = pd.DataFrame(per)
        frame.insert(0, "target", f.targets)
        frame.insert(1, "arm", arm)
        frame["pred_norm"] = pred_norm
        frame["true_norm"] = sc.true_norm
        frame["overlap"] = overlap
        self.per.append(frame)
        if not ideal:
            self.struct.append({**sc.structure(arm, st), "generator": gen})
        if is_null:
            n_sig = sig.sum(axis=1)
            up = np.where(sig, st["de"].lfc > 0, False).sum(axis=1)
            self.nulls.append(
                {
                    "arm": arm,
                    "generator": gen,
                    "null_fp_de_median": float(np.median(n_sig)),
                    "null_fp_de_mean": float(np.mean(n_sig)),
                    "null_fp_rate_median": float(np.median(n_sig) / self.tested.sum()),
                    "null_jaccard_mean": float(np.mean(n_sig == 0)),
                    "null_up_fraction": float(up.sum() / max(n_sig.sum(), 1)),
                    "null_expression_error_median": float(
                        np.median(((st["profiles"] - f.ctrl_profile) ** 2).sum(axis=1))
                    ),
                    "tested_genes": int(self.tested.sum()),
                }
            )
        return raw

    def emit(self, arm: str, gen: str, p_cpm, p_bulk, *, is_null=False):
        t0 = time.time()
        f = self.fold
        p_cpm = p_cpm / p_cpm.sum(1, keepdims=True)
        p_bulk = p_bulk / p_bulk.sum(1, keepdims=True)
        res = _pool(
            _emit,
            len(f.targets),
            self.jobs,
            gen=gen,
            targets=f.targets,
            prefix=f.seed_prefix,
            template=f.template,
            depths=f.depths,
            donors=self.donors,
            pool=f.pool,
            donor_depths=self.donor_depths,
            phi=self.phi,
            c_cpm=f.ctrl_mean,
            c_bulk=f.ctrl_bulk,
            p_cpm=p_cpm,
            p_bulk=p_bulk,
            ref=f.reference,
            tested=self.tested,
            compare=self.ref_diag if is_null else self.real_diag,
        )
        st = evaluation._stack(res)
        st["cmp"] = [r["cmp"] for r in res]
        raw = self._record(arm, gen, st, is_null=is_null)
        self.log(
            f"[{self.name}] {arm:28s} PDS={raw['pds_cosine']:.4f} "
            f"MSE={raw['expr_mse_unbiased_capped_norm']:.4f} "
            f"FID={raw['de_wilcoxon_direction_fidelity_yield_raw']:.4f} "
            f"nDE={raw['predicted_de_median']:.0f} var={raw['cells_var_ratio']:.2f} "
            f"({time.time() - t0:.0f}s)"
        )
        return raw

    def ideal(self, arm: str, p_cpm, p_bulk):
        f = self.fold
        p_cpm = p_cpm / p_cpm.sum(1, keepdims=True)
        p_bulk = p_bulk / p_bulk.sum(1, keepdims=True)
        st = {
            "profiles": np.log1p(M.BULK_TARGET_SUM * p_bulk),
            "disp": np.zeros(len(f.targets)),
            "de": realisation.ideal_de_table(p_cpm, f.reference, self.tested),
        }
        raw = self._record(arm, "ideal", st, ideal=True)
        self.log(
            f"[{self.name}] {arm:28s} PDS={raw['pds_cosine']:.4f} "
            f"MSE={raw['expr_mse_unbiased_capped_norm']:.4f} "
            f"FID={raw['de_wilcoxon_direction_fidelity_yield_raw']:.4f} "
            f"nDE={raw['predicted_de_median']:.0f} (ideal)"
        )
        return raw

    def split_half(self):
        (raw, _), st = self.scorer.split_half(generator.SEED + 2)
        self.rows["ANCHOR_split_half"] = {"generator": "real", **raw}
        self.struct.append({**self.scorer.structure("ANCHOR_split_half", st), "generator": "real"})
        self.struct.append(
            {
                **self.scorer.structure("real_perturbed", self.scorer.real_structure),
                "generator": "real",
            }
        )

    def save(self, phase: int):
        out = OUT / ("folds_smoke" if self.smoke else "folds") / self.name
        out.mkdir(parents=True, exist_ok=True)
        raw = pd.DataFrame(self.rows).T
        raw.index.name = "arm"
        raw.to_csv(out / f"phase{phase}_scores_raw.csv")
        pd.concat(self.per).to_csv(out / f"phase{phase}_per_perturbation.csv", index=False)
        if self.struct:
            pd.DataFrame(self.struct).to_csv(out / f"phase{phase}_structure.csv", index=False)
        if self.nulls:
            pd.DataFrame(self.nulls).to_csv(out / f"phase{phase}_null.csv", index=False)
        (out / f"phase{phase}_meta.json").write_text(
            json.dumps(
                {
                    "fold": self.name,
                    "n_targets": int(len(self.fold.targets)),
                    "tested_genes": int(self.tested.sum()),
                    "predictors": self.names + ["CD4"],
                    "n_donors": N_DONORS,
                    "donor_depth_median": float(np.median(self.donor_depths)),
                    "phi_median_tested": float(np.median(self.phi[self.tested])),
                    "k_sources_counts": pd.Series(self.k_sources).value_counts().to_dict(),
                    "predeclaration_sha256": PREDECL_SHA,
                },
                indent=2,
                default=int,
            )
        )


# --------------------------------------------------------------------------- phases
def check_reproduction(c2: C2Fold) -> dict:
    """G0 at a = 1 and both local anchors must equal the frozen C1 fold scores."""
    frozen = pd.read_csv(C1_FOLDS / c2.name / "scores_raw.csv", index_col=0)
    out = {}
    for ours, theirs in [
        ("G0_a1.00", "C1a"),
        ("ANCHOR_mean_response", "ANCHOR_mean_response"),
        ("ANCHOR_split_half", "ANCHOR_split_half"),
        ("G0_null", "G0_control"),
    ]:
        diffs = {
            m: abs(float(c2.rows[ours][m]) - float(frozen.loc[theirs, m]))
            for m in evaluation.MEMBERS
        }
        out[ours] = {"vs": theirs, "max_abs_diff": max(diffs.values()), "diffs": diffs}
    # amendment 2: anchors / null exact; C1a within 1e-6 (post-fold code edits on 09-26)
    tol = {"G0_a1.00": 1e-6}
    out["pass"] = all(
        v["max_abs_diff"] <= tol.get(k, 1e-12) for k, v in out.items() if isinstance(v, dict)
    )
    return out


def phase1(c2: C2Fold) -> None:
    p_cpm, p_bulk = c2.moments()
    c2.emit("G0_a1.00", "G0", p_cpm, p_bulk)
    c2.emit("ANCHOR_mean_response", "G0", *c2.anchor)
    c2.split_half()
    c2.emit("G0_null", "G0", *c2.null, is_null=True)
    repro = {"pass": True, "smoke": True} if c2.smoke else check_reproduction(c2)
    c2.log(f"[{c2.name}] reproduction pass={repro['pass']}")
    rdir = OUT / ("folds_smoke" if c2.smoke else "folds") / c2.name
    rdir.mkdir(parents=True, exist_ok=True)
    (rdir / "reproduction.json").write_text(json.dumps(repro, indent=2))
    if not repro["pass"]:
        c2.save(1)
        raise SystemExit(f"[{c2.name}] C1 reproduction failed; stopping")
    c2.ideal("IDEAL_a1.00", p_cpm, p_bulk)
    c2.ideal("IDEAL_ANCHOR_mean_response", *c2.anchor)
    c2.ideal("IDEAL_null", *c2.null)
    for gen in ["G1c", "G1ci", "G1b", "G2"]:
        c2.emit(f"{gen}_null", gen, *c2.null, is_null=True)
    g2_var = c2.rows["G2_null"]["cells_var_ratio"]
    run_g3 = bool(g2_var < G3_VARIANCE_GATE)
    c2.log(f"[{c2.name}] G2 null variance ratio {g2_var:.3f}; G3 justified here: {run_g3}")
    if run_g3:
        c2.emit("G3_null", "G3", *c2.null, is_null=True)
    for gen in ["G1c", "G1ci", "G1b", "G2"] + (["G3"] if run_g3 else []):
        c2.emit(f"{gen}_a1.00", gen, p_cpm, p_bulk)
    c2.emit("ANCHOR_mean_response_G1c", "G1c", *c2.anchor)
    c2.emit("ANCHOR_mean_response_G1ci", "G1ci", *c2.anchor)
    c2.save(1)


def phase2(c2: C2Fold) -> None:
    decision = json.loads((OUT / "generator_decision.json").read_text())
    best = decision["selected_generator"]
    ones = np.ones(len(c2.fold.targets))
    infeasible = {}
    for a in AMPLITUDES:
        shrink = None if a == 1.0 else a * ones
        p_cpm, p_bulk = c2.moments(shrink)
        c2.ideal(f"IDEAL_a{a:.2f}", p_cpm, p_bulk)
        if a == 1.0:
            continue  # a = 1 arms are in phase 1
        for gen in ["G0"] + ([best] if best != "G0" else []):
            try:
                c2.emit(f"{gen}_a{a:.2f}", gen, p_cpm, p_bulk)
            except ValueError as exc:
                # the frozen emitter cannot realise this amplitude: unbuildable, not a candidate
                infeasible[f"{gen}_a{a:.2f}"] = str(exc)
                c2.log(f"[{c2.name}] {gen}_a{a:.2f} INFEASIBLE: {exc}")
    c2.save(2)
    out = OUT / ("folds_smoke" if c2.smoke else "folds") / c2.name / "phase2_infeasible.json"
    out.write_text(json.dumps(infeasible, indent=2))


def phase3(c2: C2Fold) -> None:
    decision = json.loads((OUT / "amplitude_decision.json").read_text())
    if not decision.get("run_per_target"):
        raise SystemExit("per-target amplitude is not justified by the global study")
    gen, a_star = decision["generator"], decision["a_star"]
    k = c2.k_sources.astype(float)
    kbar = k[k > 0].mean()
    for beta in BETAS:
        shrink = np.where(k > 0, a_star * (k / kbar) ** beta, a_star)
        p_cpm, p_bulk = c2.moments(shrink)
        c2.emit(f"{gen}_a{a_star:.2f}_beta{beta:+.1f}", gen, p_cpm, p_bulk)
        c2.ideal(f"IDEAL_a{a_star:.2f}_beta{beta:+.1f}", p_cpm, p_bulk)
    c2.save(3)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", type=int, choices=[1, 2, 3], required=True)
    parser.add_argument("--folds", nargs="+", default=["H1", "K562"])
    parser.add_argument("--jobs", type=int, default=9)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    if hashlib.sha256(PREDECL.read_bytes()).hexdigest() != PREDECL_SHA:
        sys.exit("C2 predeclaration changed after it was frozen")
    if hashlib.sha256(AMENDMENT.read_bytes()).hexdigest() != AMENDMENT_SHA:
        sys.exit("C2 amendment 1 changed after it was frozen")
    if hashlib.sha256(AMENDMENT_2.read_bytes()).hexdigest() != AMENDMENT_2_SHA:
        sys.exit("C2 amendment 2 changed after it was frozen")
    OUT.mkdir(parents=True, exist_ok=True)

    def log(msg):
        print(msg, flush=True)

    for name in args.folds:
        c2 = C2Fold(name, args.jobs, log, args.smoke)
        {1: phase1, 2: phase2, 3: phase3}[args.phase](c2)


if __name__ == "__main__":
    main()
