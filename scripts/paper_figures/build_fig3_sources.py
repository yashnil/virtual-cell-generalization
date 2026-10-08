# ruff: noqa: E501 (figure labels and provenance notes are kept as single strings)
"""Figure 3 source tables: source-count ablation in each subset's own frame (N3-B, exploratory reading).

Medians are the frozen ``n3b_R_gamma.csv`` values; the per-(target, k, m) draw intervals for the own-frame
γ⊥ are recomputed with exactly the frozen aggregation (per draw: mean over the subsets of size m; then 2.5 /
97.5 percentiles over the 40 draws) from the frozen per-draw table. The medians recomputed this way are asserted
equal to the frozen table.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from _common import frozen, write_source

RG = "outputs/n3/n3b/n3b_R_gamma.csv"
DRAWS = "outputs/n3/n3b/n3b_draws.csv"
REPORTS = ["reports/n3_results.md", "reports/n3_protocol.md"]


def main() -> None:
    r = pd.read_csv(frozen(RG))
    d = pd.read_csv(frozen(DRAWS))
    e4 = d[d.estimator == "E4"]
    per_draw = e4.groupby(["context", "k", "m", "repeat", "draw"])[
        ["M3_own", "M3_ref", "M1"]
    ].mean()
    agg = (
        per_draw.groupby(["context", "k", "m"])
        .agg(
            own_median=("M3_own", "median"),
            own_lo=("M3_own", lambda x: np.percentile(x, 2.5)),
            own_hi=("M3_own", lambda x: np.percentile(x, 97.5)),
        )
        .reset_index()
    )
    t = r.merge(agg, on=["context", "k", "m"], how="left")
    chk = t[t.k > 0]
    assert np.allclose(chk.own_median, chk.R_gamma_own, atol=1e-12), (
        "own-frame medians differ from frozen table"
    )
    t.loc[t.k == 0, ["R_gamma_own", "own_median", "own_lo", "own_hi"]] = (
        0.0  # own-frame zero-shot ≡ 0
    )
    # zero-shot (k = 0) has one value per cell-split repeat: record the range over the 5 repeats
    z = (
        d[(d.estimator == "E0s") & (d.k == 0)]
        .groupby(["context", "m", "repeat"])[["M1", "M3_ref"]]
        .mean()
    )
    zr = (
        z.groupby(["context", "m"])
        .agg(
            R_full_k0_rep_min=("M1", "min"),
            R_full_k0_rep_max=("M1", "max"),
            R_gamma_k0_rep_min=("M3_ref", "min"),
            R_gamma_k0_rep_max=("M3_ref", "max"),
        )
        .reset_index()
    )
    zr["k"] = 0
    t = t.merge(zr, on=["context", "k", "m"], how="left")
    t = t.drop(columns=["own_median"]).sort_values(["context", "k", "m"]).reset_index(drop=True)
    write_source(
        "fig3_source_count",
        t,
        sources=[RG, DRAWS],
        reports=REPORTS,
        build_script="scripts/paper_figures/build_fig3_sources.py",
        plot_script="scripts/paper_figures/plot_fig3.py",
        notes=(
            "N3-B: E4 with every non-empty subset of the 5 sources (31 per target), 40 anchor draws per k, the "
            "N3-A test set (319). R_gamma = preregistered fixed-reference gamma-perp (5-source consensus); "
            "R_gamma_own = gamma-perp relative to the subset's own consensus (exploratory; 0 at k = 0 by "
            "construction); R_full = template-removed full response M1 (k = 0: zero-shot E0s). Values: per draw "
            "mean over subsets of size m, then median over draws; own_lo/own_hi and R_*_lo/hi = 2.5/97.5% over "
            "draws; for k = 0 (one value per cell-split repeat) R_*_k0_rep_min/max = range over the 5 repeats. m = 1 zero-shot uses an unshrunk single source (fit_scale = 1.0)."
        ),
    )


if __name__ == "__main__":
    main()
