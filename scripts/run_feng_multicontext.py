"""External multi-context validation of the frozen unseen-perturbation predictor.

Feng et al. screened 444 genes by CRISPRi in 19 iPSC lines under one study and
one protocol. That is the design Kaden was meant to provide and could not: the
study is held fixed while the cellular context varies, so performance can be
read against context without a confound.

**No model is fitted and no constant is changed.** The estimators, ``k``, the
ridge penalty, the STRING version, the PCA width and the scaling all come from
the frozen modules, and section 0 proves the scoring path reproduces the frozen
arch1 result before any Feng line is touched.

Usage::

    uv run python scripts/run_feng_multicontext.py
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from virtual_cell.data import scperteval
from virtual_cell.modelling import multicontext as mc
from virtual_cell.modelling import unseen_perturbation as up
from virtual_cell.modelling.external_benchmark import KNN_K, PRIOR_FAMILY

ROOT = Path(__file__).resolve().parents[1]
DESIGN = ROOT / "data" / "splits" / "four_context_v1"
CANONICAL = ROOT / "data" / "processed" / "four_context_v1"
RAW = ROOT / "data" / "raw" / "scperteval"
FENG = ROOT / "data" / "raw" / "feng" / "TargetedScreen_LFC_byGene-perLine.tsv.gz"
INVENTORY = ROOT / "data" / "provenance" / "scperteval" / "public_label_inventory.json"
FROZEN = ROOT / "outputs" / "unseen_perturbation_v1"
OUTDIR = ROOT / "outputs" / "feng_multicontext_v1"

TOLERANCE = 1e-9
#: Predeclared: a (line, target) cell is usable only with at least this many
#: cells. Set from the dataset's own power distribution before any scoring.
MIN_CELLS_PER_TARGET = 50


def rule(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def load_sources():
    contexts = (DESIGN / "contexts.txt").read_text().split()
    perts = (DESIGN / "shared_perturbations.txt").read_text().split()
    genes = (DESIGN / "shared_genes.txt").read_text().split()
    delta = np.load(CANONICAL / "delta_tensor.npy").astype(np.float64)
    with INVENTORY.open() as fh:
        inventory = json.load(fh)
    seen = set()
    for name in contexts:
        seen |= set(inventory[name]["perturbations"])
    return delta, contexts, perts, genes, seen


# --------------------------------------------------------------------------
# 0. equivalence with the frozen protocol
# --------------------------------------------------------------------------


def check_equivalence(delta, source_genes, source_perts, seen_anywhere) -> None:
    path = RAW / "arch1_processed_complete.h5ad"
    target_genes = set(scperteval.read_var_names(path))
    genes = [g for g in source_genes if g in target_genes]
    labels = set(map(str, scperteval.cell_labels(path)))
    unseen = [p for p in sorted(labels - {"control", "non-targeting"}) if p not in seen_anywhere]

    pseudo = scperteval.pseudobulk(path, genes=genes, perturbations=unseen, context="arch1")
    block = mc.build_source_block(delta, source_genes, source_perts, genes)
    vocabulary = sorted(set(source_perts) | set(unseen))
    matrix, _ = mc.string_matrix(vocabulary, ROOT)
    table, _ = mc.score_unseen(pseudo.delta, unseen, block, matrix, vocabulary)

    frozen = pd.read_csv(FROZEN / "arch1_external.csv").set_index("estimator")
    got = table.set_index("estimator")
    cols = ["beta_pearson", "beta_unexplained"]
    diff = float((got[cols] - frozen[cols]).abs().max().max())
    print(f"  max |difference| against frozen arch1 beta-level result: {diff:.3e}")
    if diff > TOLERANCE:
        raise AssertionError(
            f"the delta-matrix scoring path does not reproduce arch1 ({diff:.3e}); "
            "Feng was not scored"
        )
    print("  REPRODUCED. The scoring path is the frozen one.")


# --------------------------------------------------------------------------
# C. schema audit
# --------------------------------------------------------------------------


def load_feng(genes_wanted: set[str]):
    """Stream the per-line table into everything the phase needs.

    Two passes rather than one. The file holds 53 million rows, and
    accumulating even the filtered subset as DataFrames costs several GB of
    strings; a first pass collects the three label vocabularies and a second
    writes straight into pre-allocated arrays.

    Returns the ``(line, target, gene)`` log-fold-change stack, the per-line
    basal profile taken from ``wt_expr`` (which is exactly constant within
    ``(line, gene)`` across targets, i.e. the unperturbed expression), and a
    ``(line, target)`` count of significantly changed genes.
    """
    peek = pd.read_csv(FENG, sep="\t", nrows=5)
    print(f"  columns: {list(peek.columns)}")
    line_col = "Cell_Line"
    gene_col = "Expressed_Gene_Symbol"
    usecols = [line_col, "Target", gene_col, "wt_expr", "lfc", "pval_adj"]

    lines: set[str] = set()
    targets: set[str] = set()
    genes: set[str] = set()
    n_rows = 0
    for chunk in pd.read_csv(FENG, sep="\t", usecols=usecols, chunksize=8_000_000):
        n_rows += len(chunk)
        keep = chunk[chunk[gene_col].isin(genes_wanted)]
        lines.update(keep[line_col].astype(str).unique())
        targets.update(keep["Target"].astype(str).unique())
        genes.update(keep[gene_col].astype(str).unique())
        print(f"    pass 1: {n_rows / 1e6:.0f}M rows", flush=True)

    line_names, target_names, gene_names = sorted(lines), sorted(targets), sorted(genes)
    li = {n: i for i, n in enumerate(line_names)}
    ti = {n: i for i, n in enumerate(target_names)}
    gi = {n: i for i, n in enumerate(gene_names)}

    shape = (len(line_names), len(target_names), len(gene_names))
    stack = np.full(shape, np.nan, dtype=np.float32)
    basal = np.full((len(line_names), len(gene_names)), np.nan, dtype=np.float32)
    n_sig = np.zeros((len(line_names), len(target_names)), dtype=np.int32)
    written = np.zeros(shape, dtype=bool)
    filled = 0
    duplicated = 0
    for chunk in pd.read_csv(FENG, sep="\t", usecols=usecols, chunksize=8_000_000):
        keep = chunk[chunk[gene_col].isin(genes_wanted)]
        if keep.empty:
            continue
        i = keep[line_col].astype(str).map(li).to_numpy()
        j = keep["Target"].astype(str).map(ti).to_numpy()
        k = keep[gene_col].astype(str).map(gi).to_numpy()
        duplicated += int(written[i, j, k].sum())
        written[i, j, k] = True
        stack[i, j, k] = keep["lfc"].to_numpy(dtype=np.float32)
        basal[i, k] = keep["wt_expr"].to_numpy(dtype=np.float32)
        sig = keep["pval_adj"].to_numpy(dtype=np.float64) < 0.05
        np.add.at(n_sig, (i[sig], j[sig]), 1)
        filled += len(keep)
        print(f"    pass 2: {filled / 1e6:.0f}M rows placed", flush=True)

    meta = {
        "total_rows": n_rows,
        "rows_on_axis": filled,
        "duplicated_rows": duplicated,
        "columns": list(peek.columns),
    }
    return stack, basal, n_sig, line_names, target_names, gene_names, meta


def reliability_bound(stack: np.ndarray, index: list[int]) -> pd.DataFrame:
    """A LOWER bound on each line's reliability, from cross-line agreement.

    There are no split halves in this file and no replicate guides, so true
    split-half reliability cannot be computed. What can be derived is a bound.
    For two lines measuring the same perturbation, the observed correlation
    satisfies ``r_ij <= sqrt(rho_i * rho_j)`` -- equality only if the true
    responses were identical, which across contexts they are not. So

        rho_i * rho_j >= r_ij**2

    gives a floor, and the floor is loose by exactly the amount context
    genuinely changes the response. It is reported as a bound and never as a
    reliability estimate; the distinction matters because the whole Kaden
    lesson was that a benchmark can look fine and measure nothing.
    """
    n_lines = stack.shape[0]
    rows = []
    for a in range(n_lines):
        pairwise = []
        for b in range(n_lines):
            if a == b:
                continue
            values = []
            for p in index:
                x, y = stack[a, p], stack[b, p]
                ok = np.isfinite(x) & np.isfinite(y)
                if ok.sum() < 50:
                    continue
                xc = x[ok] - x[ok].mean()
                yc = y[ok] - y[ok].mean()
                denom = np.linalg.norm(xc) * np.linalg.norm(yc)
                if denom > 0:
                    values.append(float(xc @ yc / denom))
            if values:
                pairwise.append(float(np.median(values)))
        median_r = float(np.median(pairwise)) if pairwise else np.nan
        floor = max(median_r, 0.0) if np.isfinite(median_r) else np.nan
        rows.append(
            {
                "line_index": a,
                "median_cross_line_r": median_r,
                "reliability_lower_bound": floor**2,
                "ceiling_lower_bound": floor,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    started = time.time()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    delta, contexts, source_perts, source_genes, seen_anywhere = load_sources()
    control_means = np.load(CANONICAL / "control_means.npy").astype(np.float64)
    print(f"source tensor {delta.shape} over {contexts}")

    rule("0. PROTOCOL EQUIVALENCE — REPRODUCE THE FROZEN arch1 RESULT")
    check_equivalence(delta, source_genes, source_perts, seen_anywhere)

    rule("C. LOCAL SCHEMA AUDIT")
    stack, basal, n_sig, lines, targets, genes, meta = load_feng(set(source_genes))
    print(f"\n  rows in file             {meta['total_rows']:,}")
    print(f"  rows on our gene axis    {meta['rows_on_axis']:,}")
    print(f"  duplicated (line,target,gene) rows: {meta['duplicated_rows']}")
    print(f"  cell lines               {len(lines)}")
    print(f"  targets                  {len(targets)}")
    print(f"  output genes on our axis {len(genes)} of {len(source_genes)}")
    present = np.isfinite(stack).all(axis=2)
    print(
        f"  (line,target) present    {present.sum()} of {present.size} "
        f"({100 * present.mean():.1f}%)"
    )
    print(f"  missing values overall   {100 * np.isnan(stack).mean():.2f}%")

    unseen = [t for t in targets if t not in seen_anywhere]
    measured = [t for t in targets if t in set(source_perts)]
    ti = {t: i for i, t in enumerate(targets)}
    unseen_idx = [ti[t] for t in unseen]
    measured_idx = [ti[t] for t in measured]
    print(f"\n  globally unseen targets           {len(unseen)}")
    print(f"  directly measured in the sources  {len(measured)}")

    block = mc.build_source_block(delta, source_genes, source_perts, genes)
    vocabulary = sorted(set(source_perts) | set(unseen))
    matrix, covered = mc.string_matrix(vocabulary, ROOT)
    vocab = pd.Index(vocabulary)
    print(
        f"  STRING coverage, unseen targets   "
        f"{covered[vocab.get_indexer(pd.Index(unseen))].mean():.3f}"
    )

    rule("D. RELIABILITY QUALIFICATION")
    print("  No split halves and no replicate guides are in this file, so TRUE")
    print("  split-half reliability cannot be computed. Two things can:")
    print("    (a) a LOWER BOUND from cross-line agreement, r_ij <= sqrt(rho_i rho_j);")
    print("    (b) the dataset's own per-line signal strength, from pval_adj.\n")
    bound = reliability_bound(stack, unseen_idx)
    bound["cell_line"] = [lines[i] for i in bound["line_index"]]
    bound["n_targets_ge10_sig"] = [int((n_sig[i] >= 10).sum()) for i in bound["line_index"]]
    bound["median_sig_genes"] = [float(np.median(n_sig[i])) for i in bound["line_index"]]
    bound["n_unseen_ge10_sig"] = [
        int((n_sig[i, unseen_idx] >= 10).sum()) for i in bound["line_index"]
    ]
    bound = bound[
        [
            "cell_line",
            "median_cross_line_r",
            "ceiling_lower_bound",
            "reliability_lower_bound",
            "n_targets_ge10_sig",
            "n_unseen_ge10_sig",
            "median_sig_genes",
        ]
    ]
    print(bound.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    bound.to_csv(OUTDIR / "reliability.csv", index=False)
    print(
        f"\n  pooled over lines: {len(unseen)} unseen targets, "
        f"{int((n_sig[:, unseen_idx] >= 10).any(axis=0).sum())} with >=10 significant "
        f"genes in at least one line"
    )

    rule("E / F. FROZEN UNSEEN PREDICTION AND MEASURED-TRANSFER CONTROL, PER LINE")
    rows = []
    details = []
    for i, name in enumerate(lines):
        ok_unseen = [p for p in unseen_idx if np.isfinite(stack[i, p]).all()]
        ok_measured = [p for p in measured_idx if np.isfinite(stack[i, p]).all()]
        if len(ok_unseen) < 20 or len(ok_measured) < 20:
            print(f"  {name}: skipped (unseen={len(ok_unseen)}, measured={len(ok_measured)})")
            continue
        table, detail = mc.score_unseen(
            stack[i, ok_unseen].astype(np.float64),
            [targets[p] for p in ok_unseen],
            block,
            matrix,
            vocabulary,
        )
        knn = table.set_index("estimator").loc["U2_knn"]
        ridge = table.set_index("estimator").loc["U3_ridge"]
        control = mc.score_measured_transfer(
            stack[i, ok_measured].astype(np.float64),
            [targets[p] for p in ok_measured],
            block,
        )
        rows.append(
            {
                "cell_line": name,
                "n_unseen": len(ok_unseen),
                "knn_pearson": knn["beta_pearson"],
                "knn_spearman": knn["beta_spearman"],
                "knn_cosine": knn["beta_cosine"],
                "knn_unexplained": knn["beta_unexplained"],
                "ridge_pearson": ridge["beta_pearson"],
                "ridge_unexplained": ridge["beta_unexplained"],
                "n_measured": control["n_seen"],
                "measured_pearson": control["scaled_pearson"],
                "measured_spearman": control["scaled_spearman"],
                "measured_cosine": control["scaled_cosine"],
                "measured_unexplained": control["scaled_unexplained"],
            }
        )
        detail["cell_line"] = name
        details.append(detail)
        print(
            f"  {name:<10} unseen n={len(ok_unseen):<4} knn r={knn['beta_pearson']:+.4f} "
            f"unexp={knn['beta_unexplained']:.4f}   |   measured n={control['n_seen']:<4} "
            f"r={control['scaled_pearson']:+.4f} unexp={control['scaled_unexplained']:.4f}"
        )
    per_line = pd.DataFrame(rows)
    per_line.to_csv(OUTDIR / "per_line.csv", index=False)
    detail_frame = pd.concat(details, ignore_index=True)
    detail_frame.to_csv(OUTDIR / "per_perturbation.csv", index=False)

    print("\n  --- pooled Feng effect (mean over lines) ---")
    pooled = np.nanmean(stack, axis=0)
    pooled_unseen = [p for p in unseen_idx if np.isfinite(pooled[p]).all()]
    pooled_measured = [p for p in measured_idx if np.isfinite(pooled[p]).all()]
    ptable, pdetail = mc.score_unseen(
        pooled[pooled_unseen].astype(np.float64),
        [targets[p] for p in pooled_unseen],
        block,
        matrix,
        vocabulary,
    )
    pcontrol = mc.score_measured_transfer(
        pooled[pooled_measured].astype(np.float64),
        [targets[p] for p in pooled_measured],
        block,
    )
    print(ptable.to_string(index=False, float_format=lambda v: f"{v:+.4f}"))
    print(
        f"  measured transfer (pooled): n={pcontrol['n_seen']} "
        f"r={pcontrol['scaled_pearson']:+.4f} unexp={pcontrol['scaled_unexplained']:.4f}"
    )
    ptable.to_csv(OUTDIR / "pooled_unseen.csv", index=False)
    pd.DataFrame([pcontrol]).to_csv(OUTDIR / "pooled_measured.csv", index=False)
    pdetail.to_csv(OUTDIR / "pooled_per_perturbation.csv", index=False)

    rule("E2. IS THE ARM GAP JUST A DIFFERENCE IN EFFECT STRENGTH?")
    print("  The measured arm is our essential-gene core; the unseen arm is")
    print("  transcription factors. If the measured targets simply have stronger,")
    print("  better-determined responses, the gap is a property of the target sets")
    print("  rather than of the two prediction strategies. Matched here on the")
    print("  dataset's own per-target signal count.\n")
    sig_pooled = n_sig.sum(axis=0)
    pooled_beta = pooled - pooled[pooled_unseen].mean(axis=0)
    match_rows = []
    for label, idx, mode in (
        ("unseen", pooled_unseen, "prior"),
        ("measured", pooled_measured, "transfer"),
    ):
        names = [targets[p] for p in idx]
        truth = pooled[idx] - pooled[idx].mean(axis=0)[None, :]
        if mode == "prior":
            acc = pdetail.set_index("perturbation").loc[names, "knn_pearson"].to_numpy()
        else:
            pos = pd.Index(block.perturbations).get_indexer(pd.Index(names))
            transfer = block.delta[:, pos].mean(axis=0) - block.delta.mean(axis=(0, 1))
            acc = up.evaluate_predictions(block.scale * transfer, truth)["pearson"]
        for name, p, a in zip(names, idx, acc, strict=True):
            match_rows.append(
                {"arm": label, "perturbation": name, "n_sig": int(sig_pooled[p]), "pearson": a}
            )
    match = pd.DataFrame(match_rows)
    edges = np.quantile(match["n_sig"], [0, 0.25, 0.5, 0.75, 1.0])
    match["signal_bin"] = pd.cut(match["n_sig"], bins=np.unique(edges), include_lowest=True)
    table_m = (
        match.groupby(["signal_bin", "arm"], observed=True)["pearson"]
        .agg(["median", "size"])
        .unstack()
    )
    print(table_m.to_string(float_format=lambda v: f"{v:+.4f}"))
    match.to_csv(OUTDIR / "arm_signal_match.csv", index=False)
    del pooled_beta

    rule("G. CONTEXT DISTANCE")
    print("  Basal profiles come from the file's own wt_expr column, which is")
    print("  exactly constant within (line, gene) across targets -- the line's")
    print("  unperturbed expression. It reads no perturbation outcome.\n")
    keep = pd.Index(source_genes).get_indexer(pd.Index(genes))
    source_basal = control_means[:, keep]
    dist_rows = []
    for i, name in enumerate(lines):
        x = basal[i].astype(np.float64)
        sims = []
        for c in range(source_basal.shape[0]):
            y = source_basal[c]
            ok = np.isfinite(x) & np.isfinite(y)
            xc, yc = x[ok] - x[ok].mean(), y[ok] - y[ok].mean()
            sims.append(float(xc @ yc / (np.linalg.norm(xc) * np.linalg.norm(yc))))
        dist_rows.append(
            {
                "cell_line": name,
                "max_similarity_to_sources": max(sims),
                "mean_similarity_to_sources": float(np.mean(sims)),
                **{f"sim_{contexts[c]}": sims[c] for c in range(len(contexts))},
            }
        )
    distance = pd.DataFrame(dist_rows)
    print(distance.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    distance.to_csv(OUTDIR / "context_distance.csv", index=False)

    merged = per_line.merge(distance, on="cell_line").merge(
        bound[["cell_line", "ceiling_lower_bound", "n_unseen_ge10_sig"]], on="cell_line"
    )
    print("\n  Spearman across the 19 lines:")
    tests = []
    for x in ("mean_similarity_to_sources", "ceiling_lower_bound", "n_unseen_ge10_sig"):
        for y in ("knn_pearson", "knn_unexplained", "measured_pearson", "measured_unexplained"):
            ok = np.isfinite(merged[x]) & np.isfinite(merged[y])
            if ok.sum() < 5:
                continue
            rho = float(stats.spearmanr(merged.loc[ok, x], merged.loc[ok, y]).statistic)
            boot = [
                stats.spearmanr(
                    merged.loc[ok, x].to_numpy()[s], merged.loc[ok, y].to_numpy()[s]
                ).statistic
                for s in (
                    np.random.default_rng(k).integers(0, int(ok.sum()), int(ok.sum()))
                    for k in range(500)
                )
            ]
            lo, hi = np.nanpercentile(boot, [2.5, 97.5])
            tests.append({"x": x, "y": y, "spearman": rho, "ci_lo": lo, "ci_hi": hi})
    test_frame = pd.DataFrame(tests)
    print(test_frame.to_string(index=False, float_format=lambda v: f"{v:+.4f}"))
    test_frame.to_csv(OUTDIR / "distance_tests.csv", index=False)
    merged.to_csv(OUTDIR / "per_line_merged.csv", index=False)

    rule("H. WITHIN-FENG CONSISTENCY")
    agree = mc.cross_line_agreement(stack[:, unseen_idx])
    cons = pd.DataFrame(
        {
            "perturbation": unseen,
            "cross_line_agreement": agree["cross_line_agreement"],
            "line_specific_fraction": agree["line_specific_fraction"],
        }
    )
    print(
        f"  unseen targets: median cross-line agreement "
        f"{cons['cross_line_agreement'].median():+.4f}, median line-specific energy "
        f"{cons['line_specific_fraction'].median():.4f}"
    )
    agree_m = mc.cross_line_agreement(stack[:, measured_idx])
    print(
        f"  measured targets: median cross-line agreement "
        f"{np.nanmedian(agree_m['cross_line_agreement']):+.4f}"
    )

    knn_by_pert = detail_frame.groupby("perturbation")["knn_pearson"].median().rename("knn_pearson")
    cons = cons.merge(knn_by_pert, left_on="perturbation", right_index=True, how="left")
    ok = np.isfinite(cons["cross_line_agreement"]) & np.isfinite(cons["knn_pearson"])
    rho = float(
        stats.spearmanr(cons.loc[ok, "cross_line_agreement"], cons.loc[ok, "knn_pearson"]).statistic
    )
    print(
        f"\n  Spearman(cross-line conservation, frozen-predictor accuracy) = {rho:+.4f}"
        f"  over {int(ok.sum())} unseen perturbations"
    )
    print("  Tests whether priors predict the conserved program and fail where")
    print("  context-specific response dominates.")
    cons.to_csv(OUTDIR / "cross_line_consistency.csv", index=False)

    rule("I. NEIGHBOUR-AGREEMENT CONFIDENCE, PER LINE")
    conf_rows = []
    for name, group in detail_frame.groupby("cell_line"):
        for score in ("neighbour_agreement", "support_distance"):
            ok = np.isfinite(group[score]) & np.isfinite(group["knn_pearson"])
            if ok.sum() < 20:
                continue
            conf_rows.append(
                {
                    "cell_line": name,
                    "score": score,
                    "n": int(ok.sum()),
                    "spearman": float(
                        stats.spearmanr(
                            group.loc[ok, score], group.loc[ok, "knn_pearson"]
                        ).statistic
                    ),
                }
            )
    for score in ("neighbour_agreement", "support_distance"):
        ok = np.isfinite(pdetail[score]) & np.isfinite(pdetail["knn_pearson"])
        if ok.sum() >= 20:
            conf_rows.append(
                {
                    "cell_line": "POOLED",
                    "score": score,
                    "n": int(ok.sum()),
                    "spearman": float(
                        stats.spearmanr(
                            pdetail.loc[ok, score], pdetail.loc[ok, "knn_pearson"]
                        ).statistic
                    ),
                }
            )
    conf = pd.DataFrame(conf_rows)
    pivot = conf.pivot(index="cell_line", columns="score", values="spearman")
    print(pivot.to_string(float_format=lambda v: f"{v:+.4f}"))
    per_line_only = pivot.drop(index="POOLED", errors="ignore")
    print(
        f"\n  median across the 19 lines: neighbour_agreement "
        f"{per_line_only['neighbour_agreement'].median():+.4f}, support_distance "
        f"{per_line_only['support_distance'].median():+.4f}"
    )
    if "POOLED" in pivot.index:
        print(
            f"  POOLED (the only adequately powered estimate): neighbour_agreement "
            f"{pivot.loc['POOLED', 'neighbour_agreement']:+.4f}, support_distance "
            f"{pivot.loc['POOLED', 'support_distance']:+.4f}"
        )
    conf.to_csv(OUTDIR / "confidence_per_line.csv", index=False)

    summary = {
        "n_lines": len(lines),
        "n_targets": len(targets),
        "n_genes": len(genes),
        "n_unseen": len(unseen),
        "n_measured": len(measured),
        "prior_family": PRIOR_FAMILY,
        "knn_k": KNN_K,
        "min_cells_per_target": MIN_CELLS_PER_TARGET,
        "median_knn_pearson": float(per_line["knn_pearson"].median()),
        "median_knn_unexplained": float(per_line["knn_unexplained"].median()),
        "median_measured_pearson": float(per_line["measured_pearson"].median()),
        "median_measured_unexplained": float(per_line["measured_unexplained"].median()),
        "median_neighbour_agreement_spearman": float(per_line_only["neighbour_agreement"].median()),
        "pooled_neighbour_agreement_spearman": float(pivot.loc["POOLED", "neighbour_agreement"]),
        "conservation_vs_accuracy_spearman": rho,
        "runtime_minutes": (time.time() - started) / 60.0,
    }
    with (OUTDIR / "summary.json").open("w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)
    print(f"\nWrote {OUTDIR}   [{(time.time() - started) / 60:.1f} min]")


if __name__ == "__main__":
    main()
