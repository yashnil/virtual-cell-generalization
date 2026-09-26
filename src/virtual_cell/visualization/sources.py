"""Deterministic extraction of figure-source tables from frozen artifacts.

One function per figure. Each reads only frozen output files (or, for the Kaden
diagnostic, that phase's own outputs), computes nothing beyond selection,
relabelling and simple ratios of stored quantities, and returns a
:class:`Extract`. No number is typed in by hand.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from virtual_cell.visualization.common import ROOT

CELL_LINE = {
    "replogle22k562": "K562",
    "replogle22rpe1": "RPE1",
    "nadig25hepg2": "HepG2",
    "nadig25jurkat": "Jurkat",
}


@dataclass
class Extract:
    name: str
    table: pd.DataFrame
    sources: list[str]
    report: str
    figure_script: str
    notes: str = ""
    allow_nan: list[str] = field(default_factory=list)


def _json(rel: str) -> dict:
    return json.loads((ROOT / rel).read_text())


def _csv(rel: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / rel)


# ---------------------------------------------------------------- Figure 1


def decomposition() -> Extract:
    src = "outputs/four_context_v1/summary.json"
    s = _json(src)
    dec, sh = s["decomposition"], s["split_half"]
    frac, ss, sig = dec["fractions"], dec["sums_of_squares"], sh["signal"]
    rows = [
        ("uncorrected", "template", 100 * frac["template(mu+alpha)"], np.nan),
        ("uncorrected", "beta", 100 * frac["beta"], np.nan),
        ("uncorrected", "gamma", 100 * frac["gamma"], np.nan),
    ]
    for comp in ("template", "beta", "gamma", "noise"):
        per = np.asarray(sh["per_split_fractions"][comp]) * 100
        rows.append(("noise-corrected", comp, 100 * sh["fractions"][comp], float(per.std(ddof=1))))
    table = pd.DataFrame(rows, columns=["view", "component", "percent", "sd_over_splits"])
    rep = pd.DataFrame(
        {
            "view": "reproducibility",
            "component": ["mu", "alpha", "beta", "gamma"],
            "percent": [100 * sig[c] / ss[c] for c in ("mu", "alpha", "beta", "gamma")],
            "sd_over_splits": np.nan,
        }
    )
    table = pd.concat([table, rep], ignore_index=True)
    table["n_splits"] = sh["n_splits"]
    return Extract(
        "fig1_decomposition",
        table,
        [src],
        "reports/four_context_decomposition_v1.md",
        "scripts/figures/plot_decomposition.py",
        "Reproducibility = cross-half signal / raw sum of squares per component. "
        "sd_over_splits is the sd over the 50 split-half resamples (NaN where the "
        "quantity has no resampling distribution).",
        allow_nan=["sd_over_splits"],
    )


# ---------------------------------------------------------------- Figure 2


def decomposition_robustness() -> Extract:
    src = "outputs/four_context_sensitivity/variant_table.csv"
    v = _csv(src)
    scheme = {"shared": "shared control mean", "split": "independent control split"}
    agg = {"mean_log": "mean of log1p", "log_mean": "log of mean CP10K"}
    feat = {
        "all": "all 6,640 genes",
        "hvg4000": "4,000 control HVGs",
        "hvg2000": "2,000 control HVGs",
    }
    table = pd.DataFrame(
        {
            "control_scheme": v.scheme.map(scheme),
            "aggregation": v["agg"].map(agg),
            "seed": v.seed,
            "feature_space": v.subset.map(feat),
            "n_genes": v.n_genes,
            "beta_corrected": v.c_beta,
            "gamma_corrected": v.c_gamma,
            "gamma_uncorrected": v.u_gamma,
            "template_corrected": v.c_template,
            "noise": v.c_noise,
        }
    )
    table["canonical"] = (
        (v.scheme == "shared") & (v["agg"] == "mean_log") & (v.seed == 42) & (v.subset == "all")
    ).to_numpy()
    return Extract(
        "fig2_decomposition_robustness",
        table,
        [src],
        "reports/four_context_decomposition_sensitivity.md",
        "scripts/figures/plot_decomposition_robustness.py",
        "All 21 frozen variant x feature-space rows; percentages of response energy.",
    )


# ---------------------------------------------------------------- Figure 3


def gamma_recovery() -> Extract:
    obs_src = "outputs/pathway_falsification_v1/observed_representations.csv"
    null_src = "outputs/pathway_falsification_v1/null_replicates.csv"
    cmp_src = "outputs/pathway_falsification_v1/null_comparison.csv"
    rcmp_src = "outputs/pathway_falsification_v1/reactome_comparison.csv"
    obs = _csv(obs_src)
    rows = []
    for r in obs.itertuples():
        rows.append(
            {
                "record": "observed",
                "cell_line": r.cell_line,
                "representation": r.representation,
                "r_gamma": r.r_gamma,
                "r_gamma_normalised": r.r_gamma_normalised,
                "ceiling": r.ceiling,
                "replicate": -1,
            }
        )
    nulls = _csv(null_src)
    nulls = nulls[nulls.representation == "permuted"]
    for r in nulls.itertuples():
        rows.append(
            {
                "record": "hallmark_permuted_null",
                "cell_line": r.cell_line,
                "representation": "hallmark_null",
                "r_gamma": r.r_gamma,
                "r_gamma_normalised": r.r_gamma_normalised,
                "ceiling": r.ceiling,
                "replicate": int(r.replicate),
            }
        )
    table = pd.DataFrame(rows)
    cmp_ = _csv(cmp_src)
    cmp_ = cmp_[(cmp_.metric == "r_gamma_normalised") & (cmp_.null == "permuted")]
    rc = _csv(rcmp_src)
    rc = rc[rc.metric == "r_gamma_normalised"]
    tests = pd.concat(
        [
            cmp_.assign(record="hallmark_test")[
                ["record", "cell_line", "p_value", "z", "null_mean", "null_sd"]
            ],
            rc.assign(record="reactome_test")[
                ["record", "cell_line", "p_value", "z", "null_mean", "null_sd"]
            ],
        ]
    )
    table = table.merge(tests, on=["record", "cell_line"], how="outer")
    return Extract(
        "fig3_gamma_recovery",
        table,
        [obs_src, null_src, cmp_src, rcmp_src],
        "reports/pathway_gamma_falsification_v1.md",
        "scripts/figures/plot_gamma_recovery.py",
        "Median r(gamma_true, gamma_hat), raw and reliability-normalised. Null = "
        "100 gene-label-permuted Hallmark representations (set sizes and overlaps "
        "preserved); Reactome tested against 20 permutations of its own sets.",
        allow_nan=[
            "r_gamma",
            "r_gamma_normalised",
            "ceiling",
            "replicate",
            "p_value",
            "z",
            "null_mean",
            "null_sd",
        ],
    )


# ---------------------------------------------------------------- Figure 4


def recoverability_vs_utility() -> Extract:
    folds_src = "outputs/pathway_gamma_v2/hallmark_folds.csv"
    diag_src = "outputs/pathway_gamma_v2/residual_diagnostics.csv"
    obs_src = "outputs/pathway_falsification_v1/observed_representations.csv"
    f = _csv(folds_src)
    d = _csv(diag_src)
    o = _csv(obs_src)
    o = o[o.representation == "hallmark"][["cell_line", "r_gamma", "r_gamma_normalised"]]
    t = f.merge(d, on="cell_line").merge(o, on="cell_line")
    table = pd.DataFrame(
        {
            "cell_line": t.cell_line,
            "r_gamma_v2_model": t.r_Rhat_gamma,
            "r_gamma_hallmark_deterministic": t.r_gamma,
            "r_gamma_hallmark_normalised": t.r_gamma_normalised,
            "selected_lambda": t.selected_lambda,
            "delta_pearson_selected": t.boot_delta,
            "delta_pearson_ci_lo": t.boot_lo,
            "delta_pearson_ci_hi": t.boot_hi,
            "lambda_theory": t.lambda_theory,
            "delta_pearson_theory": t.theory_pearson - t.base_pearson,
            "base_pearson": t.base_pearson,
        }
    )
    return Extract(
        "fig4_recoverability_vs_utility",
        table,
        [folds_src, diag_src, obs_src],
        "reports/pathway_residual_model_v2_clean_gamma.md",
        "scripts/figures/plot_recoverability_vs_utility.py",
        "x = correlation of the v2 learned correction with the true interaction "
        "(r(R_hat, gamma)); y = outer change in median per-perturbation response "
        "Pearson after adding the correction at the inner-selected lambda "
        "(bootstrap CI) and at lambda_theory = (3+s)/4.",
    )


# ---------------------------------------------------------------- Figure 5


def external_generalization() -> Extract:
    bench_src = "outputs/unseen_perturbation_v1/benchmark_summary.csv"
    a_src = "outputs/unseen_perturbation_v1/arch1_external.csv"
    aseen_src = "outputs/unseen_perturbation_v1/arch1_seen_perturbation_reference.csv"
    fu_src = "outputs/feng_multicontext_v1/pooled_unseen.csv"
    fm_src = "outputs/feng_multicontext_v1/pooled_measured.csv"
    line_src = "outputs/feng_multicontext_v1/per_line.csv"
    b = _csv(bench_src)
    rows = []
    # Same aggregation as scripts/run_unseen_perturbation.py (median over the
    # summary rows): P1 over its 5 pooled folds, P2 over 5 folds x 4 held-out
    # contexts.
    for regime, stage, interval in (
        ("P1", "internal: held-out perturbation", "min-max over 5 folds"),
        (
            "P2",
            "internal: held-out perturbation + context",
            "min-max over 5 folds x 4 held-out contexts",
        ),
    ):
        sel = b[(b.regime == regime) & (b.prior_family == "string") & (b.estimator == "U2_knn")]
        vals = sel.beta_pearson.to_numpy()
        rows.append(
            {
                "panel": "A",
                "stage": stage,
                "arm": "prior",
                "pearson": float(np.median(vals)),
                "lo": float(vals.min()),
                "hi": float(vals.max()),
                "n": int(sel.n_test.max()),
                "interval": interval,
            }
        )
    a = _csv(a_src).set_index("estimator")
    rows.append(
        {
            "panel": "A",
            "stage": "arch1 (external)",
            "arm": "prior",
            "pearson": float(a.loc["U2_knn", "beta_pearson"]),
            "lo": np.nan,
            "hi": np.nan,
            "n": 100,
            "interval": "none stored",
        }
    )
    s = _csv(aseen_src).iloc[0]
    rows.append(
        {
            "panel": "A",
            "stage": "arch1 (external)",
            "arm": "direct",
            "pearson": float(s.scaled_pearson),
            "lo": np.nan,
            "hi": np.nan,
            "n": int(s.n_seen),
            "interval": "none stored",
        }
    )
    fu = _csv(fu_src).set_index("estimator")
    rows.append(
        {
            "panel": "A",
            "stage": "Feng pooled (external)",
            "arm": "prior",
            "pearson": float(fu.loc["U2_knn", "beta_pearson"]),
            "lo": np.nan,
            "hi": np.nan,
            "n": int(fu.loc["U2_knn", "n_test"]),
            "interval": "none stored",
        }
    )
    fm = _csv(fm_src).iloc[0]
    rows.append(
        {
            "panel": "A",
            "stage": "Feng pooled (external)",
            "arm": "direct",
            "pearson": float(fm.scaled_pearson),
            "lo": np.nan,
            "hi": np.nan,
            "n": int(fm.n_seen),
            "interval": "none stored",
        }
    )
    lines = _csv(line_src)
    for r in lines.itertuples():
        rows.append(
            {
                "panel": "B",
                "stage": r.cell_line,
                "arm": "prior",
                "pearson": r.knn_pearson,
                "lo": np.nan,
                "hi": np.nan,
                "n": int(r.n_unseen),
                "interval": "per line",
            }
        )
        rows.append(
            {
                "panel": "B",
                "stage": r.cell_line,
                "arm": "direct",
                "pearson": r.measured_pearson,
                "lo": np.nan,
                "hi": np.nan,
                "n": int(r.n_measured),
                "interval": "per line",
            }
        )
    return Extract(
        "fig5_external_generalization",
        pd.DataFrame(rows),
        [bench_src, a_src, aseen_src, fu_src, fm_src, line_src],
        "reports/feng_multicontext_external_validation_v1.md",
        "scripts/figures/plot_external_generalization.py",
        "Every value is a median per-perturbation Pearson between the predicted "
        "conserved effect and the held-out measurement centred on its context mean. "
        "Prior = STRING k-NN (k=25) for perturbations measured nowhere; direct = "
        "scale-calibrated transfer of perturbations measured in the sources. "
        "Perturbation sets and gene axes differ between stages.",
        allow_nan=["lo", "hi"],
    )


# ---------------------------------------------------------------- Figure 6


def arc_target_support() -> Extract:
    src = "data/splits/arc_target_support_v1.csv"
    s = _csv(src)
    rows = []
    for tier in (2, 1, 0):
        sub = s[s.support_tier == tier]
        prov = sub.datasets_perturbed.fillna("")
        for label, mask in (
            ("both arch1 and Kaden", prov.str.contains("arch1") & prov.str.contains("kaden25rpe1")),
            ("arch1 only", prov.str.contains("arch1") & ~prov.str.contains("kaden25rpe1")),
            ("Kaden only", ~prov.str.contains("arch1") & prov.str.contains("kaden25rpe1")),
            ("no public perturbation", prov == ""),
        ):
            n = int(mask.sum())
            if n:
                rows.append({"tier": tier, "provenance": label, "n_targets": n})
    table = pd.DataFrame(rows)
    table["percent_of_panel"] = 100 * table.n_targets / len(s)
    table["panel_size"] = len(s)
    return Extract(
        "fig6_arc_target_support",
        table,
        [src],
        "reports/arc_count_space_baseline_v1.md",
        "scripts/figures/plot_arc_target_support.py",
        "Frozen identifier-presence tiers: Tier 2 = perturbed in >= 2 public "
        "datasets, Tier 1 = exactly 1, Tier 0 = none.",
    )


# ---------------------------------------------------------------- Figure 7


KADEN_DIR = "outputs/kaden_source_reliability_v1"


def source_reliability() -> Extract:
    pp_src = f"{KADEN_DIR}/per_perturbation_reliability.csv"
    null_src = f"{KADEN_DIR}/null_reliability.csv"
    me_src = f"{KADEN_DIR}/main_effect_reliability.csv"
    mea_src = f"{KADEN_DIR}/main_effect_agreement.csv"
    agr_src = f"{KADEN_DIR}/agreement_summary.csv"
    pp = _csv(pp_src)
    pp = pp[pp.eligible]
    rows = []
    # A: per-perturbation reliability distribution on the identical three-source axis
    for ds, g in pp[pp.axis == "G_3"].groupby("dataset"):
        q = g.spearman_brown.quantile([0.1, 0.25, 0.5, 0.75, 0.9]).to_numpy()
        e = g.signal_energy.quantile([0.25, 0.5, 0.75]).to_numpy()
        rows.append(
            {
                "panel": "A",
                "dataset": ds,
                "axis": "G_3",
                "label": ds,
                "q10": q[0],
                "q25": q[1],
                "median": q[2],
                "q75": q[3],
                "q90": q[4],
                "value": float(g.n_cells.median()),
                "value_name": "median cells",
                "e25": e[0],
                "e50": e[1],
                "e75": e[2],
            }
        )
    nulls = _csv(null_src)
    for ds, g in nulls[nulls.axis == "G_3"].groupby("dataset"):
        q = g.spearman_brown.quantile([0.1, 0.25, 0.5, 0.75, 0.9]).to_numpy()
        rows.append(
            {
                "panel": "A_null",
                "dataset": ds,
                "axis": "G_3",
                "label": ds,
                "q10": q[0],
                "q25": q[1],
                "median": q[2],
                "q75": q[3],
                "q90": q[4],
                "value": float(g.cells.iloc[0]),
                "value_name": "cells per pseudo-perturbation",
            }
        )
    me = _csv(me_src)
    for r in me.itertuples():
        rows.append(
            {
                "panel": "B",
                "dataset": r.dataset,
                "axis": r.axis,
                "label": r.panel,
                "median": r.spearman_brown,
                "q10": r.half_pearson_q05,
                "q90": r.half_pearson_q95,
                "value": float(r.n_perturbations),
                "value_name": "n perturbations",
            }
        )
    mea = _csv(mea_src)
    for r in mea.itertuples():
        rows.append(
            {
                "panel": "C",
                "dataset": r.comparison,
                "axis": r.axis,
                "label": r.panel_x,
                "median": r.pearson,
                "q10": np.nan,
                "q90": np.nan,
                "value": r.noise_ceiling,
                "value_name": "noise ceiling",
                "e50": r.cosine,
            }
        )
    agr = _csv(agr_src)
    for r in agr.itertuples():
        rows.append(
            {
                "panel": "D",
                "dataset": r.comparison,
                "axis": r.kind,
                "label": r.subset,
                "median": r.median_pearson,
                "q10": r.pearson_ci_lo,
                "q90": r.pearson_ci_hi,
                "value": r.median_noise_ceiling,
                "value_name": "median noise ceiling",
                "e50": r.median_disattenuated,
                "e25": r.n,
            }
        )
    return Extract(
        "fig7_source_reliability",
        pd.DataFrame(rows),
        [pp_src, null_src, me_src, mea_src, agr_src],
        "reports/kaden_source_reliability_diagnostic_v1.md",
        "scripts/figures/plot_source_reliability.py",
        "New diagnostic outputs (not a prior freeze). A: Spearman-Brown reliability "
        "quantiles per dataset on the identical 3-source gene axis, with the "
        "control pseudo-perturbation null. B: main-effect Spearman-Brown reliability "
        "(q10/q90 columns hold the 5th/95th percentile of the half-level Pearson). "
        "C: observed main-effect Pearson (e50 = cosine) vs noise ceiling. D: median "
        "per-perturbation Pearson with bootstrap CI (e50 = median disattenuated, "
        "e25 = n).",
        allow_nan=["q10", "q25", "q75", "q90", "e25", "e50", "e75", "median", "value"],
    )


# ---------------------------------------------------------------- Figure 8


GENERATOR_METRICS = [
    ("pds_cosine", "perturbation discrimination\n(pds_cosine)", "higher"),
    ("de_wilcoxon_direction_fidelity_yield_raw", "DE direction fidelity", "higher"),
    ("de_wilcoxon_direction_reach_raw", "DE direction reach", "higher"),
    ("de_wilcoxon_sig_jaccard", "DE significance overlap\n(Jaccard)", "higher"),
    ("de_wilcoxon_lfc_nmae", "DE log-fold-change error\n(nMAE)", "lower"),
    ("expr_mse_unbiased_capped_norm", "expression error\n(normalised MSE)", "lower"),
]


def count_generator_benchmark() -> Extract:
    src = "outputs/arc_count_space_v1/generator_scores_raw.csv"
    g = _csv(src).set_index("generator")
    rows = []
    for gen in ("G0_control_resample", "G1_transport", "G2_count_model"):
        for key, label, better in GENERATOR_METRICS:
            rows.append(
                {
                    "generator": gen,
                    "metric": key,
                    "metric_label": label,
                    "better": better,
                    "value": float(g.loc[gen, key]),
                }
            )
    return Extract(
        "fig8_count_generator_benchmark",
        pd.DataFrame(rows),
        [src],
        "reports/arc_count_space_baseline_v1.md",
        "scripts/figures/plot_count_generator_benchmark.py",
        "Raw (unnormalised) local vcc2026 metric values; 300 perturbations x 400 "
        "generated cells scored against real held-out K562 cells and a disjoint "
        "2,000-cell control reference.",
    )


# ---------------------------------------------------------------- Figure 9


def transport_fidelity() -> Extract:
    s_src = "outputs/arc_count_space_v1/generator_structure.csv"
    f_src = "outputs/arc_count_space_v1/generator_mean_response_fidelity.csv"
    s = _csv(s_src)
    keep = [
        "real_control_reference",
        "real_perturbed",
        "G0_control_resample",
        "G1_transport",
        "G1_transport_unsmoothed",
        "G2_count_model",
    ]
    s = s[s.group.isin(keep)]
    rows = []
    for r in s.itertuples():
        for q, lib, gd in (
            ("q10", r.library_q10, r.genes_detected_q10),
            ("median", r.library_median, r.genes_detected_median),
            ("q90", r.library_q90, r.genes_detected_q90),
        ):
            rows.append(
                {
                    "panel": "structure",
                    "group": r.group,
                    "quantity": "library_size",
                    "stat": q,
                    "value": lib,
                }
            )
            rows.append(
                {
                    "panel": "structure",
                    "group": r.group,
                    "quantity": "genes_detected",
                    "stat": q,
                    "value": gd,
                }
            )
        rows.append(
            {
                "panel": "structure",
                "group": r.group,
                "quantity": "sparsity",
                "stat": "mean",
                "value": r.sparsity,
            }
        )
        rows.append(
            {
                "panel": "structure",
                "group": r.group,
                "quantity": "heterogeneity",
                "stat": "median",
                "value": r.heterogeneity_median,
            }
        )
    f = _csv(f_src)
    for r in f.itertuples():
        for stat, val in (
            ("realised_norm", r.realised_norm),
            ("pearson_vs_own", r.pearson_vs_own),
            ("slope_vs_own", r.slope_vs_own),
            ("own_intended_norm", r.own_intended_norm),
        ):
            rows.append(
                {
                    "panel": "fidelity",
                    "group": r.generator,
                    "quantity": "mean_response",
                    "stat": stat,
                    "value": val,
                }
            )
    return Extract(
        "fig9_transport_fidelity",
        pd.DataFrame(rows),
        [s_src, f_src],
        "reports/arc_count_space_baseline_v1.md",
        "scripts/figures/plot_transport_fidelity.py",
        "Only summary statistics were frozen (quantiles, sparsity, heterogeneity, "
        "fitted slope and r of realised vs intended pseudobulk response over 50 "
        "perturbations). Per-gene and per-cell values were not retained; the figure "
        "shows the frozen summaries rather than regenerating cells.",
        allow_nan=["value"],
    )


EXTRACTORS: dict[str, Callable[[], Extract]] = {
    "fig1_decomposition": decomposition,
    "fig2_decomposition_robustness": decomposition_robustness,
    "fig3_gamma_recovery": gamma_recovery,
    "fig4_recoverability_vs_utility": recoverability_vs_utility,
    "fig5_external_generalization": external_generalization,
    "fig6_arc_target_support": arc_target_support,
    "fig7_source_reliability": source_reliability,
    "fig8_count_generator_benchmark": count_generator_benchmark,
    "fig9_transport_fidelity": transport_fidelity,
}
