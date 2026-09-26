"""C1 public leave-one-atlas-out folds: C1a, the predeclared shrinkage family, the C0
X-Atlas ablation and the promoter ablation. Predeclaration:
``reports/competition_v2/c1_predeclaration.md`` (its SHA-256 is checked before scoring).

    uv run python scripts/competition_v2/run_c1_public_folds.py [--folds H1 K562 CD4]
        [--smoke] [--jobs 9]

Outputs: ``outputs/competition_v2/c1_license_clean/folds/<fold>/``.
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

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from virtual_cell.arc import metrics as M  # noqa: E402
from virtual_cell.competition_v2 import atlas, evaluation, fusion, generator, sources  # noqa: E402
from virtual_cell.competition_v2.licensing import assert_sources_allowed  # noqa: E402

BASE = ROOT / "outputs" / "competition_v2" / "c1_license_clean"
SRC = BASE / "sources"
C0_DATA = ROOT / "outputs" / "competition_v2" / "atlasshift_c0" / "data"
RAW = ROOT / "data" / "raw" / "competition_v2"
OFFICIAL_GENES = ROOT / "data" / "raw" / "arc2026" / "controls" / "gene_names.csv"
PREDECL = ROOT / "reports" / "competition_v2" / "c1_predeclaration.md"
PREDECL_SHA = "fd4efa6f6f69871ec9600ca894465755e726f6f15cc33f45f8d1be4eeda32957"
SEED = generator.SEED
S1_GRID = [0.25, 0.5, 0.75]
S2_GRID = [0.0, 0.25, 0.5]
N_REF = 4000


def official_genes() -> np.ndarray:
    return pd.read_csv(OFFICIAL_GENES).gene_name.astype(str).to_numpy()


def load_green() -> dict[str, sources.SourceStats]:
    return {
        "K562": sources.load_source(SRC / "K562_GWPS_CPM_full_statistics.npz"),
        "H1": sources.load_source(SRC / "H1_2025_full_statistics.npz"),
    }


def load_cd4() -> dict:
    with np.load(SRC / "CD4_DE_statistics.npz") as d:
        return {k: d[k] for k in d.files}


def load_xatlas() -> dict[str, sources.SourceStats]:
    """BLOCKED sources: loaded only for the predeclared X-Atlas ablation arms."""
    return {
        "HCT116": sources.load_source(C0_DATA / "HCT116_full_statistics.npz"),
        "HEK293T": sources.load_source(C0_DATA / "HEK293T_full_statistics.npz"),
    }


# ----------------------------------------------------------------------------- folds
def fold_h1(smoke: bool) -> evaluation.Fold:
    x, labels, all_genes = evaluation.load_cells(SRC / "fold_cells_H1.npz")
    h1_genes = pd.read_csv(RAW / "gene_names.csv").iloc[:, 0].astype(str).to_numpy()
    official = official_genes()
    genes = official[np.isin(official, h1_genes)]
    x = x[:, pd.Index(all_genes).get_indexer(genes)]
    state = json.loads((SRC / "fold_cells_H1_rng_state.json").read_text())
    rng = np.random.default_rng(SEED)
    rng.bit_generator.state = state
    ctrl = x[labels == sources.CONTROL_LABEL]
    perm = rng.permutation(ctrl.shape[0])  # continues the C0 harness's generator
    targets = np.array(sorted(set(labels) - {sources.CONTROL_LABEL}))
    if smoke:
        targets = targets[:12]
    real = [x[labels == t].toarray().astype(np.int32) for t in targets]
    return evaluation.Fold(
        "H1",
        targets,
        genes,
        real,
        reference=ctrl[perm[:N_REF]].toarray().astype(np.int64),
        pool=ctrl[perm[N_REF:]],
        seed_prefix="H1:",
    )


def fold_k562(smoke: bool) -> evaluation.Fold:
    x, labels, all_genes = evaluation.load_cells(SRC / "fold_cells_K562.npz")
    official = official_genes()
    genes = official[np.isin(official, all_genes)]
    x = x[:, pd.Index(all_genes).get_indexer(genes)]
    stats = sources.load_source(SRC / "K562_GWPS_CPM_full_statistics.npz")
    targets = np.array([t for t in stats.usable() if (labels == t).any()])
    if smoke:
        targets = targets[:12]
    ctrl = x[labels == sources.CONTROL_LABEL]
    perm = np.random.default_rng(SEED + 1).permutation(ctrl.shape[0])
    real = [x[labels == t].toarray().astype(np.int32) for t in targets]
    return evaluation.Fold(
        "K562",
        targets,
        genes,
        real,
        reference=ctrl[perm[:N_REF]].toarray().astype(np.int64),
        pool=ctrl[perm[N_REF:]],
        seed_prefix="K562:",
    )


# ----------------------------------------------------------------------------- arms
def arm_effects(fold_name, targets, genes, controls, green, cd4, xatlas, kind):
    """Fused effects for one arm, with the held-out source removed."""
    held = {"H1": "H1", "K562": "K562", "CD4": "CD4"}[fold_name]
    g = {k: v for k, v in green.items() if k != held}
    use_cd4 = cd4 if held != "CD4" else None
    if kind == "C1a":
        names = list(g)
        assert_sources_allowed(names + (["CD4"] if use_cd4 is not None else []))
        return fusion.fused_effects(
            [g[n] for n in names], [1.0] * len(names), use_cd4, 1.0, targets, genes, controls
        ), names
    c0w = {"K562": 2.0, "H1": 2.0, "HCT116": 1.0, "HEK293T": 1.0}
    pool = dict(g)
    if kind == "C0_withX":
        pool.update(xatlas)
    names = list(pool)
    return fusion.fused_effects(
        [pool[n] for n in names], [c0w[n] for n in names], use_cd4, 0.5, targets, genes, controls
    ), names


def coverage(names, srcs, cd4, targets, genes):
    """Usable direct-evidence sources per target (C1 definition: >= 20 cells / CD4 QC)."""
    cols = {}
    for n in names:
        cols[n] = np.isin(targets, srcs[n].usable())
    if cd4 is not None:
        _, avail = atlas.cd4_effect(cd4, targets, genes)
        cols["CD4"] = avail.any(axis=1)
    return pd.DataFrame(cols, index=targets)


def heldout_truth_effect(fold_name, targets, genes, green):
    """Held-out source's own centred log2fc effect: used ONLY by the scorer."""
    src = green[fold_name]
    full, mask = atlas.source_effect(src, src.targets, genes, space="log2fc")
    lookup = {t: i for i, t in enumerate(src.targets)}
    idx = np.array([lookup[t] for t in targets])
    return full[idx], mask[idx]


def row_cosine(a, b, keep):
    a, b = a[:, keep].astype(float), b[:, keep].astype(float)
    na, nb = np.linalg.norm(a, axis=1), np.linalg.norm(b, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where((na > 0) & (nb > 0), (a * b).sum(1) / (na * nb), np.nan)


def run_full_fold(fold: evaluation.Fold, green, cd4, xatlas, pairs_all, jobs, out: Path, log):
    t0 = time.time()
    out.mkdir(parents=True, exist_ok=True)
    scorer = evaluation.FoldScorer(fold, jobs=jobs)
    log(
        f"[{fold.name}] truth tables: {len(fold.targets)} targets, {len(fold.genes)} genes, "
        f"{int(scorer.tested.sum())} tested ({time.time() - t0:.0f}s)"
    )
    pairs = generator.promoter_pairs(pairs_all, fold.targets, fold.genes)
    controls = fold.controls
    effects = {}
    for kind in ["C1a", "C0_noX", "C0_withX"]:
        effects[kind] = arm_effects(
            fold.name, fold.targets, fold.genes, controls, green, cd4, xatlas, kind
        )
    c1_names = effects["C1a"][1]
    held_cd4 = cd4
    green_pred = [green[n] for n in c1_names]
    vectors = fusion.panel_source_vectors(green_pred, held_cd4, fold.targets, fold.genes)
    agree, n_src = fusion.source_agreement(vectors, fold.exclude)

    moments = {}
    base = effects["C1a"][0]
    moments["C1a"] = generator.expected_moments(
        base, controls, fold.targets, fold.genes, pairs=pairs
    )
    moments["C1a_no_promoter"] = generator.expected_moments(
        base, controls, fold.targets, fold.genes
    )
    lambdas = {}
    for s in S1_GRID:
        lam = fusion.shrinkage_factor("S1", len(fold.targets), scalar=s)
        lambdas[f"S1_l{s}"] = lam
    for f in S2_GRID:
        lam = fusion.shrinkage_factor("S2", len(fold.targets), floor=f, agreement=agree)
        lambdas[f"S2_f{f}"] = lam
    for name, lam in lambdas.items():
        moments[name] = generator.expected_moments(
            base, controls, fold.targets, fold.genes, shrink=lam, pairs=pairs
        )
    for kind in ["C0_noX", "C0_withX"]:
        moments[kind] = generator.expected_moments(
            effects[kind][0], controls, fold.targets, fold.genes, pairs=pairs
        )
    zero = {s: np.zeros_like(base[s]) for s in base}
    moments["G0_control"] = generator.expected_moments(zero, controls, fold.targets, fold.genes)
    real_cpm = np.stack([(b / b.sum(1, keepdims=True)).mean(0) for b in fold.real])
    real_bulk = np.stack([b.sum(0) / b.sum() for b in fold.real])
    moments["ANCHOR_mean_response"] = (
        np.tile(real_cpm.mean(0), (len(fold.targets), 1)),
        np.tile(real_bulk.mean(0), (len(fold.targets), 1)),
    )

    raw_rows, per_rows, structure = {}, [], []
    for name, (p_cpm, p_bulk) in moments.items():
        p_cpm = p_cpm / p_cpm.sum(1, keepdims=True)
        p_bulk = p_bulk / p_bulk.sum(1, keepdims=True)
        st = scorer.emit(name, p_cpm, p_bulk)
        raw, per = scorer.members(st)
        raw_rows[name] = raw
        structure.append(scorer.structure(name, st))
        frame = pd.DataFrame({k: v for k, v in per.items()})
        frame.insert(0, "target", fold.targets)
        frame.insert(1, "arm", name)
        frame["pred_norm"] = np.linalg.norm(
            (st["profiles"] - fold.ctrl_profile)[:, scorer.keep], axis=1
        )
        frame["true_norm"] = scorer.true_norm
        per_rows.append(frame)
        log(
            f"[{fold.name}] {name:22s} PDS={raw['pds_cosine']:.4f} "
            f"FID={raw['de_wilcoxon_direction_fidelity_yield_raw']:.4f} "
            f"nDE={raw['predicted_de_median']:.0f} ({time.time() - t0:.0f}s)"
        )
    (sh_raw, _), sh_st = scorer.split_half(SEED + 2)
    raw_rows["ANCHOR_split_half"] = sh_raw
    structure.append(scorer.structure("ANCHOR_split_half", sh_st))
    structure.append(scorer.structure("real_perturbed", scorer.real_structure))

    raw = pd.DataFrame(raw_rows).T
    raw.index.name = "arm"
    scaled = evaluation.scale_local(raw)
    raw.to_csv(out / "scores_raw.csv")
    scaled.to_csv(out / "scores_scaled_local.csv")
    pd.concat(per_rows).to_csv(out / "per_perturbation.csv", index=False)
    pd.DataFrame(structure).to_csv(out / "structure.csv", index=False)

    truth_eff, truth_mask = heldout_truth_effect(fold.name, fold.targets, fold.genes, green)
    keep = scorer.keep
    cov_c1 = coverage(c1_names, green, cd4, fold.targets, fold.genes)
    cov_x = coverage(effects["C0_withX"][1], {**green, **xatlas}, cd4, fold.targets, fold.genes)
    pt = pd.DataFrame(
        {
            "target": fold.targets,
            "agreement": agree,
            "n_sources_c1": cov_c1.sum(1).to_numpy(),
            "n_sources_withX": cov_x.sum(1).to_numpy(),
            "transfer_cosine_c1a": row_cosine(
                np.where(truth_mask, effects["C1a"][0]["log2fc"], 0), truth_eff, keep
            ),
            "transfer_cosine_c0_withX": row_cosine(
                np.where(truth_mask, effects["C0_withX"][0]["log2fc"], 0), truth_eff, keep
            ),
            **{f"lambda_{k}": v for k, v in lambdas.items()},
        }
    )
    pt.to_csv(out / "per_target_agreement.csv", index=False)
    summary = {
        "fold": fold.name,
        "n_targets": int(len(fold.targets)),
        "n_genes": int(len(fold.genes)),
        "tested_genes": int(scorer.tested.sum()),
        "c1_predictors": c1_names + ["CD4"],
        "withX_predictors": effects["C0_withX"][1] + ["CD4"],
        "coverage_c1_any": int((cov_c1.sum(1) > 0).sum()),
        "coverage_withX_any": int((cov_x.sum(1) > 0).sum()),
        "mean_sources_c1": float(cov_c1.sum(1).mean()),
        "mean_sources_withX": float(cov_x.sum(1).mean()),
        "promoter_pairs": int(len(pairs)),
        "agreement_defined": int(np.isfinite(agree).sum()),
        "predeclaration_sha256": PREDECL_SHA,
        "elapsed_seconds": round(time.time() - t0, 1),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    return raw, scaled, pt


def run_cd4_fold(green, cd4, xatlas, out: Path, log, smoke: bool):
    out.mkdir(parents=True, exist_ok=True)
    genes_all = official_genes()
    genes = genes_all[np.isin(genes_all, cd4["genes"].astype(str))]
    usable = cd4["available"] & cd4["quality_pass"] & (cd4["n_cells"] >= 20)
    targets = cd4["targets"].astype(str)[usable.any(axis=0)]
    if smoke:
        targets = targets[:12]
    truth, avail = atlas.cd4_effect(cd4, targets, genes)
    exclude = np.isin(genes, targets)
    dummy = {"log2fc": None, "bulk_delta": None}
    rows, per = {}, {}
    eff = {}
    for kind in ["C1a", "C0_noX", "C0_withX"]:
        e, names = arm_effects("CD4", targets, genes, dummy, green, None, xatlas, kind)
        eff[kind] = (e["log2fc"], names)
    vectors = fusion.panel_source_vectors([green[n] for n in eff["C1a"][1]], None, targets, genes)
    agree, n_src = fusion.source_agreement(vectors, exclude)
    arms = {"C1a": eff["C1a"][0], "C0_noX": eff["C0_noX"][0], "C0_withX": eff["C0_withX"][0]}
    for s in S1_GRID:
        arms[f"S1_l{s}"] = s * eff["C1a"][0]
    for f in S2_GRID:
        lam = fusion.shrinkage_factor("S2", len(targets), floor=f, agreement=agree)
        arms[f"S2_f{f}"] = lam[:, None] * eff["C1a"][0]
    for name, pred in arms.items():
        pds = M.pds_cosine(pred, truth, exclude=exclude)
        per[name] = pds
        rows[name] = {
            "pds_effect": float(pds.mean()),
            "covered": int((np.abs(pred[:, ~exclude]).sum(1) > 0).sum()),
        }
        log(f"[CD4] {name:12s} effect-PDS={rows[name]['pds_effect']:.4f}")
    res = pd.DataFrame(rows).T
    res.index.name = "arm"
    res.to_csv(out / "scores_effect.csv")
    keep = ~exclude
    cov_c1 = coverage(eff["C1a"][1], green, None, targets, genes)
    cov_x = coverage(eff["C0_withX"][1], {**green, **xatlas}, None, targets, genes)
    pd.DataFrame(
        {
            "target": targets,
            "agreement": agree,
            "n_sources_c1": cov_c1.sum(1).to_numpy(),
            "n_sources_withX": cov_x.sum(1).to_numpy(),
            "transfer_cosine_c1a": row_cosine(eff["C1a"][0], truth, keep),
            "transfer_cosine_c0_withX": row_cosine(eff["C0_withX"][0], truth, keep),
            **{f"pds_{k}": v for k, v in per.items()},
        }
    ).to_csv(out / "per_target_agreement.csv", index=False)
    (out / "summary.json").write_text(
        json.dumps(
            {
                "fold": "CD4",
                "n_targets": int(len(targets)),
                "n_genes": int(len(genes)),
                "coverage_c1_any": int((cov_c1.sum(1) > 0).sum()),
                "coverage_withX_any": int((cov_x.sum(1) > 0).sum()),
                "mean_sources_c1": float(cov_c1.sum(1).mean()),
                "mean_sources_withX": float(cov_x.sum(1).mean()),
                "agreement_defined": int(np.isfinite(agree).sum()),
                "predeclaration_sha256": PREDECL_SHA,
            },
            indent=2,
        )
    )
    return res


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--folds", nargs="+", default=["CD4", "H1", "K562"])
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--jobs", type=int, default=9)
    args = parser.parse_args()
    if hashlib.sha256(PREDECL.read_bytes()).hexdigest() != PREDECL_SHA:
        sys.exit("predeclaration changed after it was frozen")
    out_root = BASE / ("folds_smoke" if args.smoke else "folds")

    def log(msg):
        print(msg, flush=True)

    green, cd4, xatlas = load_green(), load_cd4(), load_xatlas()
    tss = generator.gencode_tss(RAW / "gencode.v47.annotation.gtf.gz")
    for name in args.folds:
        if name == "CD4":
            run_cd4_fold(green, cd4, xatlas, out_root / "CD4", log, args.smoke)
            continue
        fold = fold_h1(args.smoke) if name == "H1" else fold_k562(args.smoke)
        run_full_fold(fold, green, cd4, xatlas, tss, args.jobs, out_root / name, log)


if __name__ == "__main__":
    main()
