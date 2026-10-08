# ruff: noqa: E501 (figure labels and provenance notes are kept as single strings)
"""Figure 4 source tables: within-study source compatibility for the K562-essential target (N5).

Every estimate and interval is a stored N5 value. The per-perturbation differences in panel C are differences of
the stored per-perturbation Pearson r; the summary fractions are asserted to match the frozen report.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from _common import frozen, write_source

DEC = "outputs/n5/n5_decision.json"
SUM = "outputs/n5/n5_summary.json"
TERMS = "outputs/n5/per_perturbation_terms.csv"
REPORTS = ["reports/n5_results.md", "reports/n5_protocol.md"]

RELATION = {
    "K562_GWPS": ("K562 GWPS", "same cell line · same study · different screen"),
    "K562_GWPS_depth": ("K562 GWPS, depth-matched", "subsampled to RPE1 cell counts"),
    "Jurkat": ("Jurkat", "different line · companion study (Nadig 2025)"),
    "HepG2": ("HepG2", "different line · companion study (Nadig 2025)"),
    "HCT116": ("HCT116", "different line · other lab (X-Atlas)"),
    "HEK293T": ("HEK293T", "different line · other lab (X-Atlas)"),
    "RPE1": ("RPE1", "different line · same study & library"),
}


def main() -> None:
    dec = json.loads(frozen(DEC).read_text())
    summ = json.loads(frozen(SUM).read_text())
    rel = summ["median_source_reliability"]

    comp = pd.DataFrame(
        [
            {
                "source": s,
                "label": RELATION[s][0],
                "relation": RELATION[s][1],
                "C_S": v["point"],
                "lo": v["lo"],
                "hi": v["hi"],
                "median_source_reliability": rel[s],
                "median_source_cells": summ["median_source_cells"][s],
            }
            for s, v in dec["C_S"].items()
        ]
    )

    c = dec["contrasts"]
    a2 = dec["Adj2_GWPS_minus_RPE1"]
    rob = pd.DataFrame(
        [
            (
                "cosine",
                "Primary: pooled latent cosine, Δ = C_GWPS − C_RPE1",
                c["delta_BC"]["point"],
                c["delta_BC"]["lo"],
                c["delta_BC"]["hi"],
                "paired perturbation bootstrap, 2,000",
            ),
            (
                "cosine",
                "GWPS subsampled to RPE1 cell counts (Adj-4)",
                c["delta_BC_depth"]["point"],
                c["delta_BC_depth"]["lo"],
                c["delta_BC_depth"]["hi"],
                "paired perturbation bootstrap, 2,000",
            ),
            (
                "pearson",
                "Raw median per-perturbation r difference",
                dec["R1_diff_GWPS_RPE1"]["median"],
                dec["R1_diff_GWPS_RPE1"]["lo"],
                dec["R1_diff_GWPS_RPE1"]["hi"],
                "perturbation bootstrap, 2,000",
            ),
            (
                "pearson",
                "Reliability + magnitude adjusted, perturbation FE (Adj-2)",
                a2["full"]["beta_minus_RPE1"],
                a2["full"]["lo"],
                a2["full"]["hi"],
                "perturbation-cluster bootstrap, 2,000",
            ),
            (
                "pearson",
                "Adj-2 with depth-matched GWPS",
                a2["depth"]["beta_minus_RPE1"],
                a2["depth"]["lo"],
                a2["depth"]["hi"],
                "perturbation-cluster bootstrap, 2,000",
            ),
            (
                "pearson",
                f"Reliability-matched perturbations, n = {dec['Adj3']['n']} (Adj-3)",
                dec["Adj3"]["median_diff"],
                dec["Adj3"]["lo"],
                dec["Adj3"]["hi"],
                "perturbation bootstrap",
            ),
        ],
        columns=["estimand", "variant", "estimate", "lo", "hi", "interval"],
    )

    terms = pd.read_csv(frozen(TERMS))
    r = terms.pivot(index="perturbation", columns="source", values="r")
    pp = pd.DataFrame(
        {
            "perturbation": r.index,
            "diff_full": (r["K562_GWPS"] - r["RPE1"]).to_numpy(),
            "diff_depth_matched": (r["K562_GWPS_depth"] - r["RPE1"]).to_numpy(),
        }
    )
    f_full = float((pp.diff_full > 0).mean())
    f_depth = float((pp.diff_depth_matched > 0).mean())
    assert len(pp) == 1054 and abs(f_full - 0.858) < 0.0005 and abs(f_depth - 0.705) < 0.0005, (
        f_full,
        f_depth,
    )
    assert abs(np.median(pp.diff_full) - 0.138) < 0.0005

    common = dict(
        reports=REPORTS,
        build_script="scripts/paper_figures/build_fig4_sources.py",
        plot_script="scripts/paper_figures/plot_fig4.py",
    )
    write_source(
        "fig4_compatibility",
        comp,
        sources=[DEC, SUM],
        notes=(
            "C_S = pooled noise-corrected (both-sided) latent cosine between template-removed source and target "
            "(K562 essential) response fields, 1,054 perturbations x 6,408 genes; 95% paired perturbation "
            "bootstrap (2,000)."
        ),
        **common,
    )
    write_source(
        "fig4_robustness",
        rob,
        sources=[DEC],
        notes=(
            "GWPS minus RPE1 advantage under each preregistered control. Two estimands: pooled latent cosine "
            "difference, and per-perturbation template-removed Pearson difference (raw median; Adj-2 coefficient "
            "difference; Adj-3 matched-pair median)."
        ),
        **common,
    )
    write_source(
        "fig4_per_perturbation",
        pp,
        sources=[TERMS],
        notes="Per-perturbation r_GWPS - r_RPE1 (template-removed Pearson across 6,408 genes), n = 1,054.",
        **common,
    )


if __name__ == "__main__":
    main()
