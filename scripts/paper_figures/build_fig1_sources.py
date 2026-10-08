# ruff: noqa: E501 (figure labels and provenance notes are kept as single strings)
"""Figure 1 source table: component energy shares and reproducibility (four-context decomposition v1).

Reads only frozen artefacts. Template reproducibility (μ+α combined) is the SS-weighted combination of the
stored μ and α values; every other number is a stored value.
"""

from __future__ import annotations

import json

import pandas as pd
from _common import frozen, write_source

SUMMARY = "outputs/four_context_v1/summary.json"
VARIANTS = "outputs/four_context_sensitivity/variant_table.csv"


def main() -> None:
    s = json.loads(frozen(SUMMARY).read_text())
    ss = s["decomposition"]["sums_of_squares"]
    sig = s["split_half"]["signal"]
    frac = s["split_half"]["fractions"]
    raw = s["decomposition"]["fractions"]
    v = pd.read_csv(frozen(VARIANTS))

    v_rep_template = (v.rep_mu * v.u_mu + v.rep_alpha * v.u_alpha) / (v.u_mu + v.u_alpha)
    rows = []
    for comp, label in [
        ("template", "Template (μ + α)"),
        ("beta", "Conserved effect β"),
        ("gamma", "Context interaction γ"),
        ("noise", "Measurement noise"),
    ]:
        row = {"component": comp, "label": label, "energy_pct_corrected": 100 * frac[comp]}
        if comp == "template":
            row["energy_pct_uncorrected"] = 100 * raw["template(mu+alpha)"]
            row["reproducible_pct"] = 100 * (sig["mu"] + sig["alpha"]) / (ss["mu"] + ss["alpha"])
            row["reproducible_pct_variant_min"] = float(v_rep_template.min())
            row["reproducible_pct_variant_max"] = float(v_rep_template.max())
            row["energy_pct_variant_min"], row["energy_pct_variant_max"] = (
                v.c_template.min(),
                v.c_template.max(),
            )
        elif comp in ("beta", "gamma"):
            row["energy_pct_uncorrected"] = 100 * raw[comp]
            row["reproducible_pct"] = 100 * sig[comp] / ss[comp]
            row["reproducible_pct_variant_min"] = float(v[f"rep_{comp}"].min())
            row["reproducible_pct_variant_max"] = float(v[f"rep_{comp}"].max())
            row["energy_pct_variant_min"] = float(v[f"c_{comp}"].min())
            row["energy_pct_variant_max"] = float(v[f"c_{comp}"].max())
        else:
            row["energy_pct_uncorrected"] = float("nan")
            row["reproducible_pct"] = float("nan")
            row["reproducible_pct_variant_min"] = row["reproducible_pct_variant_max"] = float("nan")
            row["energy_pct_variant_min"], row["energy_pct_variant_max"] = (
                v.c_noise.min(),
                v.c_noise.max(),
            )
        rows.append(row)
    t = pd.DataFrame(rows)
    # sd over 50 split-half resamples, from the frozen extraction (percentage points)
    f1 = pd.read_csv(frozen("data/figure_sources/fig1_decomposition.csv"))
    sd = f1[f1.view == "noise-corrected"].set_index("component").sd_over_splits
    t["energy_sd_over_50_splits_pp"] = t.component.map(sd)
    t["n_variants"] = len(v)

    write_source(
        "fig1_components",
        t,
        sources=[SUMMARY, VARIANTS, "data/figure_sources/fig1_decomposition.csv"],
        reports=[
            "reports/four_context_decomposition_v1.md",
            "reports/four_context_decomposition_sensitivity.md",
        ],
        build_script="scripts/paper_figures/build_fig1_sources.py",
        plot_script="scripts/paper_figures/plot_fig1.py",
        notes=(
            "Energy shares are % of total response energy (4 contexts x 1,264 perturbations x 6,640 genes). "
            "Corrected = split-half noise correction, mean of 50 resamples. Reproducible % = cross-half signal SS "
            "/ raw SS. Template reproducibility combines mu and alpha weighted by SS. Variant min/max = range over "
            "the 21 predeclared preprocessing variants (not a confidence interval)."
        ),
    )


if __name__ == "__main__":
    main()
