"""Arc Challenge Bridge — coverage audit, score accounting, generator validation.

Reads only:
  * the official controls bundle (basal statistics and label lists; never
    perturbation outcomes, never identities),
  * ``data/provenance/scperteval/public_label_inventory.json``, the frozen gene
    and perturbation label lists of the seven public scPertEval datasets,
  * the published ``vcc2026`` metric anchors, transcribed in ``ANCHORS``.

Writes ``outputs/arc_bridge_v1/``. No model is fitted and no submission is made.

Reproduce: ``uv run python scripts/run_arc_bridge.py``
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from virtual_cell.arc import generate, metrics
from virtual_cell.arc.panel import build_panel_map
from virtual_cell.data import arc2026

ROOT = Path(__file__).resolve().parents[1]
CONTROLS = ROOT / "data" / "raw" / "arc2026" / "controls"
INVENTORY = ROOT / "data" / "provenance" / "scperteval" / "public_label_inventory.json"
OUTDIR = ROOT / "outputs" / "arc_bridge_v1"

#: The four contexts the frozen research programme used, then the three that
#: were audited during this phase purely for Arc coverage.
RESEARCH_CONTEXTS = ("replogle22k562", "replogle22rpe1", "nadig25hepg2", "nadig25jurkat")

#: Published reference points, ``vcc2026-metrics-brief.md`` section 8, measured
#: on the three official bundles at ``cell-eval2 0.15.0``, ``rule_version`` 3.
#: Each entry is ``(baseline, replicate)`` as a midpoint of the quoted range.
ANCHORS: dict[str, tuple[float, float]] = {
    "pds_cosine": (0.500, 0.9555),
    "expr_mse_unbiased_capped_norm": (0.989, 0.0365),
    "de_wilcoxon_direction_fidelity_yield_raw": (0.5135, 0.8135),
    "de_wilcoxon_direction_reach_raw": (0.072, 0.968),
    "de_wilcoxon_sig_jaccard": (0.029, 0.399),
    "de_wilcoxon_lfc_nmae": (1.0013, 0.400),
}

#: What each metric reads for a submission that emits the reference control
#: unchanged, i.e. predicts no effect. Every value is stated analytically in
#: the metric reference; see ``reports/arc_bridge_v1.md`` for the derivations.
CONTROL_SUBMISSION: dict[str, float] = {
    "pds_cosine": 0.500,
    "expr_mse_unbiased_capped_norm": 1.0032,
    "de_wilcoxon_direction_fidelity_yield_raw": 0.0,
    "de_wilcoxon_direction_reach_raw": 0.0,
    "de_wilcoxon_sig_jaccard": 0.0,
    "de_wilcoxon_lfc_nmae": 1.0,
}


def _rule(title: str) -> None:
    print("\n" + "=" * 74)
    print(title)
    print("=" * 74)


def load_inventory() -> dict[str, dict]:
    with INVENTORY.open() as fh:
        return json.load(fh)


def coverage_tables(inventory: dict[str, dict]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Which Arc genes are measured, and which Arc targets are perturbed, where."""
    arc_genes = pd.Index(arc2026.load_gene_names(CONTROLS)).astype(object)
    arc_targets = pd.Index(arc2026.load_pert_counts(CONTROLS)["target_gene"]).astype(object)

    genes = pd.DataFrame(index=arc_genes)
    targets = pd.DataFrame(index=arc_targets)
    for name, rec in sorted(inventory.items()):
        gset = set(rec["genes"])
        pset = set(rec["perturbations"]) - {"control", "non-targeting"}
        genes[name] = [g in gset for g in arc_genes]
        targets[name] = [t in pset for t in arc_targets]
    genes["n_datasets"] = genes.sum(axis=1)
    targets["n_datasets"] = targets.sum(axis=1)
    return (
        genes.rename_axis("arc_gene").reset_index(),
        targets.rename_axis("arc_target").reset_index(),
    )


def score_accounting() -> pd.DataFrame:
    """What the two reference submissions score, member by member.

    The mean-response baseline is 0 by construction. The control-emitting
    submission is the one our frozen point predictor degenerates to wherever no
    perturbation-specific evidence exists, so its score is the number that
    matters.
    """
    rows = []
    for member, (baseline, replicate) in ANCHORS.items():
        control = CONTROL_SUBMISSION[member]
        scaled = metrics.scale_score(control, baseline, replicate)
        if member == "expr_mse_unbiased_capped_norm":
            scaled = float(np.clip(scaled, 0.0, 1.0))
        rows.append(
            {
                "member": member,
                "higher_is_better": metrics.SCORED_METRICS[member],
                "baseline_b": baseline,
                "replicate_r": replicate,
                "control_submission_u": control,
                "control_submission_score": scaled,
                "mean_response_score": 0.0,
            }
        )
    return pd.DataFrame(rows)


def validate_generators(rng: np.random.Generator) -> pd.DataFrame:
    """Run each generator on a slice of the real control cells and check the output.

    This uses the official controls as *inference inputs only* — library sizes
    and composition. No perturbation outcome is read.
    """
    rows = []
    for context in ("A", "B", "C"):
        adata = arc2026.subsample_context(
            arc2026.context_path(CONTROLS, context), n_cells=3000, seed=0
        )
        ctrl = np.asarray(adata.X.todense()).astype(np.int64)
        lfc = np.zeros(ctrl.shape[1])
        cpm = ctrl.sum(axis=0) + 1.0
        produced = {
            "G0_resample": generate.resample_controls(ctrl, 400, rng=rng),
            "G1_transport": generate.transport_controls(ctrl, lfc, 400, rng=rng),
            "G2_count_model": generate.count_model(ctrl, cpm, 400, rng=rng),
        }
        for name, out in produced.items():
            lib = out.sum(axis=1)
            rows.append(
                {
                    "context": context,
                    "generator": name,
                    "cells": out.shape[0],
                    "genes": out.shape[1],
                    "integer": bool(np.issubdtype(out.dtype, np.integer)),
                    "non_negative": bool((out >= 0).all()),
                    "median_library": float(np.median(lib)),
                    "max_library": int(lib.max()),
                    "under_cap": bool(lib.max() <= generate.MAX_COUNTS_PER_CELL),
                    "density": float((out > 0).mean()),
                    "control_median_library": float(np.median(ctrl.sum(axis=1))),
                    "control_density": float((ctrl > 0).mean()),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    inventory = load_inventory()

    _rule("1. COVERAGE OF THE ARC PANEL BY PUBLIC PERTURBATION DATA")
    genes, targets = coverage_tables(inventory)
    datasets = [c for c in genes.columns if c not in ("arc_gene", "n_datasets")]
    per_dataset = pd.DataFrame(
        {
            "dataset": datasets,
            "arc_genes_measured": [int(genes[d].sum()) for d in datasets],
            "pct_genes": [100 * genes[d].mean() for d in datasets],
            "arc_targets_perturbed": [int(targets[d].sum()) for d in datasets],
            "pct_targets": [100 * targets[d].mean() for d in datasets],
            "research_context": [d in RESEARCH_CONTEXTS for d in datasets],
        }
    )
    print(per_dataset.to_string(index=False, float_format=lambda v: f"{v:.2f}"))

    covered_g = int((genes["n_datasets"] > 0).sum())
    covered_t = int((targets["n_datasets"] > 0).sum())
    multi_t = int((targets["n_datasets"] > 1).sum())
    research_t = int(targets[list(RESEARCH_CONTEXTS)].any(axis=1).sum())
    print(
        f"\n  union: Arc genes measured anywhere      {covered_g:>6}/18533 "
        f"({100 * covered_g / 18533:.1f}%)"
    )
    print(
        f"  union: Arc targets perturbed anywhere   {covered_t:>6}/300   "
        f"({100 * covered_t / 300:.1f}%)"
    )
    print(
        f"  Arc targets in >=2 datasets             {multi_t:>6}/300   "
        f"({100 * multi_t / 300:.1f}%)  <- what a conserved effect needs"
    )
    print(f"  Arc targets in the 4 research contexts  {research_t:>6}/300")

    genes.to_csv(OUTDIR / "gene_coverage.csv", index=False)
    targets.to_csv(OUTDIR / "target_coverage.csv", index=False)
    per_dataset.to_csv(OUTDIR / "dataset_coverage.csv", index=False)

    _rule("2. PANEL SUPPORT CATEGORIES")
    union_genes = sorted({g for rec in inventory.values() for g in rec["genes"]})
    arc_genes = arc2026.load_gene_names(CONTROLS)
    pmap = build_panel_map(arc_genes, union_genes)
    print("  for a target WITH public response data:")
    print(pmap.counts().to_string())
    print(
        "\n  for a target WITHOUT public response data, every measured gene "
        "falls back to\n  'measured_no_response': no per-perturbation value is "
        "estimable at all."
    )

    _rule("3. WHAT THE TRIVIAL SUBMISSIONS SCORE")
    acct = score_accounting()
    print(acct.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    control_total = acct["control_submission_score"].mean()
    print(f"\n  mean-response baseline, per definition : {0.0:+.4f}")
    print(f"  control-emitting submission            : {control_total:+.4f}")
    print("\n  Predicting no effect is WORSE than predicting the panel's mean")
    print("  response. The scale's zero is the mean response, not the control.")
    acct.to_csv(OUTDIR / "score_accounting.csv", index=False)

    _rule("4. GENERATOR VALIDATION AGAINST THE OFFICIAL CONTROLS")
    gen = validate_generators(np.random.default_rng(20260920))
    print(gen.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    gen.to_csv(OUTDIR / "generator_validation.csv", index=False)
    assert gen["integer"].all() and gen["non_negative"].all() and gen["under_cap"].all()

    summary = {
        "arc_genes_measured_union": covered_g,
        "arc_genes_total": 18533,
        "arc_targets_perturbed_union": covered_t,
        "arc_targets_multi_dataset": multi_t,
        "arc_targets_in_research_contexts": research_t,
        "arc_targets_total": 300,
        "control_submission_score": float(control_total),
        "mean_response_baseline_score": 0.0,
        "panel_support": {k: int(v) for k, v in pmap.counts().items()},
        "datasets": {d: int(targets[d].sum()) for d in datasets},
    }
    with (OUTDIR / "summary.json").open("w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)
    print(f"\nWrote {OUTDIR}")


if __name__ == "__main__":
    main()
