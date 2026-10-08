# ruff: noqa: E501 (figure labels and provenance notes are kept as single strings)
"""Extended Data source tables (ED Figs 1–5). Frozen inputs only; derived quantities are named in each note.

Numbers that exist only in a frozen report (no machine-readable artefact) are transcribed with their section and
recorded under ``transcribed_from_reports`` in the provenance sidecar; where an artefact also exists the value is
asserted against it.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from _common import frozen, write_source

CONTEXTS = ["K562", "RPE1", "HepG2", "Jurkat", "HCT116", "HEK293T"]
B = "scripts/paper_figures/build_ext_sources.py"


# ---------------------------------------------------------------- ED 1: decomposition robustness
def ed1() -> None:
    vt = "outputs/four_context_sensitivity/variant_table.csv"
    de = "outputs/four_context_sensitivity/depth_experiment.csv"
    v = pd.read_csv(frozen(vt))
    feat = {
        "all": "all 6,640 genes",
        "hvg4000": "4,000 control HVGs",
        "hvg2000": "2,000 control HVGs",
    }
    t = pd.DataFrame(
        {
            "control_scheme": v.scheme.map(
                {"shared": "shared control mean", "split": "independent control split"}
            ),
            "aggregation": v["agg"].map(
                {"mean_log": "mean of log1p", "log_mean": "log of mean CP10K"}
            ),
            "seed": v.seed,
            "feature_space": v.subset.map(feat),
            "template_corrected": v.c_template,
            "beta_corrected": v.c_beta,
            "gamma_corrected": v.c_gamma,
            "noise": v.c_noise,
            "beta_uncorrected": v.u_beta,
            "gamma_uncorrected": v.u_gamma,
            "rep_beta": v.rep_beta,
            "rep_gamma": v.rep_gamma,
        }
    )
    t["canonical"] = (
        (v.scheme == "shared") & (v["agg"] == "mean_log") & (v.seed == 42) & (v.subset == "all")
    ).to_numpy()
    write_source(
        "ext1_variants",
        t,
        sources=[vt],
        reports=["reports/four_context_decomposition_sensitivity.md"],
        build_script=B,
        plot_script="scripts/paper_figures/plot_ext_fig1.py",
        notes="All 21 predeclared variant x feature-space rows. Shares are % of total response energy; rep_* = % reproducible (cross-half signal / raw SS).",
    )

    d = pd.read_csv(frozen(de))
    med = d.groupby("n_cells").reliability.median()
    # bootstrap CIs over pairs are stored only in the report table (sensitivity §B)
    ci = {15: (0.0812, 0.1106), 30: (0.1490, 0.2052), 50: (0.2256, 0.2908), 100: (0.3733, 0.4551)}
    rep_med = {15: 0.0966, 30: 0.1762, 50: 0.2626, 100: 0.4188}
    for n, m in rep_med.items():
        assert abs(med.loc[n] - m) < 5e-5, (n, med.loc[n], m)
    dt = pd.DataFrame(
        {
            "cells_per_half": med.index,
            "median_reliability": med.to_numpy(),
            "ci_lo": [ci[n][0] for n in med.index],
            "ci_hi": [ci[n][1] for n in med.index],
            "n_pairs": d.groupby("n_cells").size().to_numpy(),
        }
    )
    write_source(
        "ext1_depth",
        dt,
        sources=[de],
        reports=["reports/four_context_decomposition_sensitivity.md"],
        build_script=B,
        plot_script="scripts/paper_figures/plot_ext_fig1.py",
        notes="Controlled depth experiment: the same 643 (context, perturbation) pairs at n = 15/30/50/100 cells per half; median split-half reliability over pairs.",
        report_values={
            "ci_lo/ci_hi": "four_context_decomposition_sensitivity.md §B, 'Median reliability by depth (bootstrap 95 % CI over pairs)'; medians asserted against depth_experiment.csv"
        },
    )


# ---------------------------------------------------------------- ED 2: N1/N3 validity controls
def _q(x):
    return np.percentile(x, 2.5), np.percentile(x, 97.5)


def ed2() -> None:
    f3a = "outputs/n3/n3a/n3a_f3a.csv"
    f3b = "outputs/n3/n3a/n3a_f3b_summary.csv"
    sm = "outputs/n3/n3a/n3a_summary.csv"
    dr = "outputs/n3/n3a/n3a_draws.csv"
    n1dr = "outputs/n1_n4/n1/n1_draws.csv"
    a = pd.read_csv(frozen(f3a))
    s = pd.read_csv(frozen(sm))
    rows = []
    for c in CONTEXTS:
        for k in (20, 100):
            p = a[(a.context == c) & (a.estimator == "E4") & (a.k == k)].M3
            u = s[(s.context == c) & (s.estimator == "E4") & (s.k == k)].iloc[0]
            lo, hi = _q(p)
            rows.append((c, k, "permuted anchors", p.median(), lo, hi, len(p)))
            rows.append((c, k, "correct anchors", u.M3_median, u.M3_lo, u.M3_hi, int(u.n_draws)))
    perm = pd.DataFrame(
        rows, columns=["context", "k", "variant", "M3_median", "lo", "hi", "n_draws"]
    )

    fb = pd.read_csv(frozen(f3b), header=[0, 1], index_col=[0, 1, 2])
    fb.columns = [f"{m}_{v}" for m, v in fb.columns]
    fb = fb.reset_index()
    fb.columns = ["context", "estimator", "k"] + list(fb.columns[3:])
    fb["dM0_shared_minus_split"] = fb.M0_shared_control - fb.M0_split_control
    fb["dM3_shared_minus_split"] = fb.M3_shared_control - fb.M3_split_control

    d = pd.read_csv(frozen(dr))
    harm = []
    for c in CONTEXTS:
        x = d[d.context == c]
        e0s = x[(x.estimator == "E0s") & (x.k == 0)].set_index("repeat").M0
        for k in (1, 2, 5, 10, 20, 50):
            g = x[(x.estimator == "E1") & (x.k == k)]
            gain = g.M0.to_numpy() - e0s.loc[g.repeat].to_numpy()
            lo, hi = _q(gain)
            harm.append((c, k, np.median(gain), lo, hi))
    harm = pd.DataFrame(harm, columns=["context", "k", "E1_M0_gain_median", "lo", "hi"])

    n1 = pd.read_csv(frozen(n1dr))
    rep = []
    for study, dd, kref, ctxs in (
        ("N1 (3 sources, k_ref 885)", n1, 885, CONTEXTS[:4]),
        ("N3-A (5 sources, k_ref 743)", d, 743, CONTEXTS),
    ):
        for c in ctxs:
            x = dd[(dd.context == c) & (dd.estimator == "E4")]
            ref = x[x.k == kref].set_index("repeat").M3
            g = x[x.k == 20]
            ratio = g.M3.to_numpy() / ref.loc[g.repeat].to_numpy()
            lo, hi = _q(ratio)
            rep.append((study, c, g.M3.median() / ref.median(), lo, hi))
    rep = pd.DataFrame(rep, columns=["study", "context", "ratio_k20_kref", "lo", "hi"])

    common = dict(
        reports=["reports/n1_n4_results.md", "reports/n3_results.md"],
        build_script=B,
        plot_script="scripts/paper_figures/plot_ext_fig2.py",
    )
    write_source(
        "ext2_permutation",
        perm,
        sources=[f3a, sm],
        notes="F3a anchor-permutation null (N3-A, E4): anchors' fit rows permuted among anchors before fitting (50 draws) vs the unpermuted 200 draws; median and 2.5-97.5% over draws.",
        **common,
    )
    write_source(
        "ext2_shared_control",
        fb,
        sources=[f3b],
        notes="F3b (N3-A): medians over 50 paired draws with one shared control mean vs disjoint control parts; differences shared - split.",
        **common,
    )
    write_source(
        "ext2_single_anchor",
        harm,
        sources=[dr],
        notes="E1 (target template from k anchors) M0 gain over the same repeat's zero-shot E0s; median and 2.5-97.5% over 200 draws.",
        **common,
    )
    write_source(
        "ext2_n1_vs_n3",
        rep,
        sources=[n1dr, dr],
        notes="E4 gamma-perp recovery at k = 20 as a fraction of k_ref: ratio of draw medians (preregistered C3 quantity); interval = per-draw ratio with k_ref paired by repeat (descriptive).",
        **common,
    )


# ---------------------------------------------------------------- ED 3: agreement correction
def ed3() -> None:
    hb = "data/figure_sources/n1_n4/exploratory_halves_bias.csv"
    t2a, t2b = "data/figure_sources/n1_n4/n4_t2.csv", "outputs/n3/n4/n4_t2.csv"
    t1b = "outputs/n3/n4/n4_t1.csv"
    h = pd.read_csv(frozen(hb))
    v1 = {
        "K562": 0.554,
        "RPE1": 0.786,
        "HepG2": 0.656,
        "Jurkat": 0.556,
    }  # n1_n4_results §3.4 (frozen v1 report)
    h["rho_frozen_v1_report"] = h.context.map(v1)
    a = pd.read_csv(frozen(t2a)).assign(study="N1/N4 (4 folds, 3 sources)")
    b = pd.read_csv(frozen(t2b)).assign(study="N3 (6 folds, 5 sources)")
    t2 = pd.concat([a, b], ignore_index=True)
    t1 = pd.read_csv(frozen(t1b))
    t1 = t1[
        t1.statistic.isin(
            [
                "b1_magnitude",
                "b2_source_reliability",
                "b3_source_reliable_energy",
                "b6_noise_agreement_ceiling",
            ]
        )
    ]
    common = dict(
        reports=[
            "reports/n1_n4_results.md",
            "reports/n3_results.md",
            "reports/transferability_confidence_model_v1.md",
        ],
        build_script=B,
        plot_script="scripts/paper_figures/plot_ext_fig3.py",
    )
    write_source(
        "ext3_halves_bias",
        h,
        sources=[hb],
        notes="Spearman rho(agreement, transfer quality -D) on the same disjoint N1 splits: per-repeat unbiased energies (N4 estimator) vs halves averaged over repeats first (v1 construction). Exploratory, not preregistered.",
        report_values={
            "rho_frozen_v1_report": "n1_n4_results.md §3.4 table, row 'frozen v1 report (M0, transferability_confidence_model_v1.md §3)'"
        },
        **common,
    )
    write_source(
        "ext3_partial",
        t2,
        sources=[t2a, t2b],
        notes="T2 partial Spearman of agreement with transfer quality, controlling for b1-b5; 2,000-resample bootstrap 95% CI.",
        **common,
    )
    write_source(
        "ext3_t1",
        t1,
        sources=[t1b],
        notes="T1 (N3, six folds): rho(agreement) - |rho(b_j)| with paired bootstrap 95% CI for the four required baselines.",
        **common,
    )


# ---------------------------------------------------------------- ED 4: N6 gate
def ed4() -> None:
    gate = "outputs/n6/source_gate.json"
    dec = "outputs/n6/n6_decision.json"
    terms = "outputs/n6/per_perturbation_terms.csv"
    g = json.loads(frozen(gate).read_text())
    dd = json.loads(frozen(dec).read_text())
    pooled = pd.DataFrame(
        [
            (
                "VIPerturb-seq K562",
                g["viperturb_pooled_split_half_reliability"],
                "frozen gate statistic (d1)",
                g["median_source_cells"]["VIPerturb_K562"],
            ),
            (
                "K562 GWPS at VIPerturb depth",
                0.109,
                "exploratory post-hoc (n6_results §4)",
                g["median_source_cells"]["K562_GWPS_vipdepth"],
            ),
            (
                "K562 GWPS, full depth",
                0.341,
                "exploratory post-hoc (n6_results §4)",
                g["median_source_cells"]["K562_GWPS"],
            ),
        ],
        columns=["source", "pooled_split_half_reliability", "status", "median_cells"],
    )
    pooled["gate"] = 0.10
    t = pd.read_csv(frozen(terms))
    t = t[t.source.isin(["VIPerturb_K562", "K562_GWPS_vipdepth", "K562_GWPS"])][
        ["source", "perturbation", "rel"]
    ]
    pc = pd.DataFrame(
        [
            (
                "K562 GWPS, full depth",
                dd["C_S"]["K562_GWPS"]["point"],
                dd["C_S"]["K562_GWPS"]["lo"],
                dd["C_S"]["K562_GWPS"]["hi"],
            ),
            (
                "K562 GWPS at VIPerturb depth",
                dd["C_S"]["K562_GWPS_vipdepth"]["point"],
                dd["C_S"]["K562_GWPS_vipdepth"]["lo"],
                dd["C_S"]["K562_GWPS_vipdepth"]["hi"],
            ),
        ],
        columns=["source", "C_S", "lo", "hi"],
    )
    gl = pd.DataFrame(
        [
            (
                "d1",
                "Pooled split-half reliability ≥ 0.10",
                f"{g['viperturb_pooled_split_half_reliability']:.4f}",
                dd["gate_failures"]["d1_low_pooled_reliability"],
            ),
            (
                "d2",
                "Reliable energy 2.5th pct > 0",
                f"{g['viperturb_reliable_energy_boot_lo']:,.0f}",
                dd["gate_failures"]["d2_reliable_energy_not_positive"],
            ),
            (
                "d3",
                "Positive control stable at VIPerturb depth (|ΔC| ≤ 0.10)",
                f"{abs(dd['C_S']['K562_GWPS_vipdepth']['point'] - dd['C_S']['K562_GWPS']['point']):.3f}",
                dd["gate_failures"]["d3_positive_control_unstable_at_vip_depth"],
            ),
            (
                "d4",
                "Panel ≥ 300 perturbations",
                str(g["n_panel"]),
                dd["gate_failures"]["d4_panel_too_small"],
            ),
        ],
        columns=["item", "condition_to_pass", "observed", "failed"],
    )
    assert dd["outcome"] == "D"
    common = dict(
        reports=["reports/n6_results.md", "reports/n6_protocol.md"],
        build_script=B,
        plot_script="scripts/paper_figures/plot_ext_fig4.py",
    )
    write_source(
        "ext4_pooled_reliability",
        pooled,
        sources=[gate],
        notes="Pooled split-half reliability (mean over repeats of the pooled Pearson of centred split halves). Only the VIPerturb value is a frozen gate statistic.",
        report_values={
            "0.109 and 0.341": "n6_results.md §4 table (exploratory post-hoc, read-only, not in outputs/n6)"
        },
        **common,
    )
    write_source(
        "ext4_per_perturbation_rel",
        t,
        sources=[terms],
        notes="Per-perturbation split-half reliability (637 panel perturbations) for the K562 sources.",
        **common,
    )
    write_source(
        "ext4_positive_control",
        pc,
        sources=[dec],
        notes="C_S vs the K562-essential target on the N6 panel (637 x 6,083); 95% paired bootstrap. Only the positive-control rows are interpretable under Outcome D.",
        **common,
    )
    write_source(
        "ext4_gate",
        gl,
        sources=[gate, dec],
        notes="Preregistered Outcome-D gate items; one failing item yields Outcome D.",
        **common,
    )


# ---------------------------------------------------------------- ED 5: N3 reference frame
def ed5() -> None:
    rg = "outputs/n3/n3b/n3b_R_gamma.csv"
    dr = "outputs/n3/n3b/n3b_draws.csv"
    dec = "outputs/n3/n3_decision.json"
    r = pd.read_csv(frozen(rg))
    d = pd.read_csv(frozen(dr))
    nb = json.loads(frozen(dec).read_text())["N3B"]["per_target"]
    e4 = d[(d.estimator == "E4") & (d.k == 20)]
    pdm = e4.groupby(["context", "m", "repeat", "draw"]).M3_own.mean().unstack("m")
    rows = []
    for c in CONTEXTS:
        p = nb[c]
        rows.append(
            (
                c,
                "preregistered fixed frame",
                p["delta20_5v3_median"],
                p["delta20_5v3_lo"],
                p["delta20_5v3_hi"],
            )
        )
        diff = (pdm.loc[c][5] - pdm.loc[c][3]).to_numpy()
        lo, hi = _q(diff)
        rows.append((c, "own frame (exploratory)", float(np.median(diff)), lo, hi))
    ch = pd.DataFrame(rows, columns=["context", "frame", "delta20_m5_minus_m3", "lo", "hi"])
    write_source(
        "ext5_delta20",
        ch,
        sources=[dec, dr],
        reports=["reports/n3_results.md", "reports/n3_protocol.md"],
        build_script=B,
        plot_script="scripts/paper_figures/plot_ext_fig5.py",
        notes=(
            "Delta(20) = R(20, m=5) - R(20, m=3), paired per draw (40 draws), 2.5-97.5% draw interval. "
            "Fixed frame = the preregistered Q1 quantity, copied from n3_decision.json. Own frame = the same "
            "pairing on M3_own from n3b_draws.csv (exploratory)."
        ),
    )
    write_source(
        "ext5_reference_frame",
        r,
        sources=[rg],
        reports=["reports/n3_results.md", "reports/n3_protocol.md"],
        build_script=B,
        plot_script="scripts/paper_figures/plot_ext_fig5.py",
        notes="Frozen N3-B R(k, m) table: R_gamma (preregistered, fixed 5-source reference) and R_gamma_own (subset's own consensus).",
    )


if __name__ == "__main__":
    ed1()
    ed2()
    ed3()
    ed4()
    ed5()
