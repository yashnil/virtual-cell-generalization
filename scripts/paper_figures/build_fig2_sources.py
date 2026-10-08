# ruff: noqa: E501 (figure labels and provenance notes are kept as single strings)
"""Figure 2 source tables: component-resolved target-budget learning curves (N3-A, six contexts).

Reads the frozen N3-A summaries, perturbation bootstrap and per-draw metrics. Derived quantities:

* template/scale recovery = [M0_E2(k) − M0_E0s] / [M0_E2(k_ref) − M0_E0s] from draw medians (the preregistered
  k_T50 construct, N3 protocol §1.2);
* G(k) = median M3_E4(k) / median M3_E4(k_ref) (the preregistered C3 ratio at k = 20) and the same ratio for
  template/scale; its interval is the 2.5–97.5 % range over draws of the per-draw ratio, the k_ref value
  being paired by cell-split repeat. This interval is descriptive and was not part of any decision.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from _common import frozen, write_source

SUMMARY = "outputs/n3/n3a/n3a_summary.csv"
BOOT = "outputs/n3/n3a/n3a_boot.csv"
DRAWS = "outputs/n3/n3a/n3a_draws.csv"
DECISION = "outputs/n3/n3_decision.json"
KREF = 743
CONTEXTS = ["K562", "RPE1", "HepG2", "Jurkat", "HCT116", "HEK293T"]
REPORTS = ["reports/n3_results.md", "reports/n3_protocol.md", "reports/n1_n4_protocol.md"]


def curves() -> pd.DataFrame:
    s = pd.read_csv(frozen(SUMMARY))
    b = pd.read_csv(frozen(BOOT))
    dec = json.loads(frozen(DECISION).read_text())["N3A"]
    rows = []
    for c in CONTEXTS:
        x = s[s.context == c]
        e0s_m0 = float(x[(x.estimator == "E0s") & (x.k == 0)].M0_median.iloc[0])
        e2 = x[x.estimator == "E2"].set_index("k").M0_median
        gain_ref = float(e2.loc[KREF] - e0s_m0)
        e4 = x[x.estimator == "E4"].set_index("k")
        bb = b[(b.context == c) & (b.estimator == "E4") & (b.metric == "M3")].set_index("k")
        for k in [0, 1, 2, 5, 10, 20, 50, 100, 200, KREF]:
            if k == 0:
                tmpl_gain, m3, m3lo, m3hi = 0.0, 0.0, 0.0, 0.0
                bp = bl = bh = 0.0
            else:
                tmpl_gain = float(e2.loc[k] - e0s_m0)
                m3, m3lo, m3hi = (float(e4.loc[k, f"M3_{q}"]) for q in ("median", "lo", "hi"))
                bp, bl, bh = (float(bb.loc[k, q]) for q in ("point", "boot_lo", "boot_hi"))
            rows.append(
                {
                    "context": c,
                    "k": k,
                    "template_gain_M0": tmpl_gain,
                    "template_gain_ref_M0": gain_ref,
                    "template_frac": tmpl_gain / gain_ref,
                    "gamma_M3_median": m3,
                    "gamma_M3_draw_lo": m3lo,
                    "gamma_M3_draw_hi": m3hi,
                    "gamma_M3_boot_point": bp,
                    "gamma_M3_boot_lo": bl,
                    "gamma_M3_boot_hi": bh,
                    "k_T50": dec[c]["k_T50"],
                    "k_gamma50": dec[c]["k_gamma50"],
                    "class_N1_rule": dec[c]["class"],
                    "E0s_M0": e0s_m0,
                }
            )
    t = pd.DataFrame(rows)
    # integrity: our template_frac reproduces the frozen k_T50
    for c in CONTEXTS:
        x = t[(t.context == c) & (t.k > 0)]
        kt = int(x[x.template_frac >= 0.5].k.min())
        assert kt == dec[c]["k_T50"], (c, kt, dec[c]["k_T50"])
        r = float(x[x.k == 20].gamma_M3_median.iloc[0] / x[x.k == KREF].gamma_M3_median.iloc[0])
        assert abs(r - dec[c]["E4_M3_ratio_k20_kref"]) < 1e-9, (c, r)
    return t


def gain_fraction() -> pd.DataFrame:
    d = pd.read_csv(frozen(DRAWS))
    rows = []
    for c in CONTEXTS:
        x = d[d.context == c]
        e0s = x[(x.estimator == "E0s") & (x.k == 0)].set_index("repeat").M0
        ref_e4 = x[(x.estimator == "E4") & (x.k == KREF)].set_index("repeat").M3
        ref_e2 = x[(x.estimator == "E2") & (x.k == KREF)].set_index("repeat").M0
        for k in [20, 50, 100, 200]:
            g = x[(x.estimator == "E4") & (x.k == k)]
            ratio = g.M3.to_numpy() / ref_e4.loc[g.repeat].to_numpy()
            point = g.M3.median() / ref_e4.median()
            rows.append(
                {
                    "context": c,
                    "k": k,
                    "component": "gamma_perp",
                    "point": point,
                    "lo": np.percentile(ratio, 2.5),
                    "hi": np.percentile(ratio, 97.5),
                    "n_draws": len(g),
                }
            )
            t = x[(x.estimator == "E2") & (x.k == k)]
            num = t.M0.to_numpy() - e0s.loc[t.repeat].to_numpy()
            den = ref_e2.loc[t.repeat].to_numpy() - e0s.loc[t.repeat].to_numpy()
            point_t = (t.M0.median() - e0s.median()) / (ref_e2.median() - e0s.median())
            rt = num / den
            rows.append(
                {
                    "context": c,
                    "k": k,
                    "component": "template_scale",
                    "point": point_t,
                    "lo": np.percentile(rt, 2.5),
                    "hi": np.percentile(rt, 97.5),
                    "n_draws": len(t),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    write_source(
        "fig2_curves",
        curves(),
        sources=[SUMMARY, BOOT, DECISION],
        reports=REPORTS,
        build_script="scripts/paper_figures/build_fig2_sources.py",
        plot_script="scripts/paper_figures/plot_fig2.py",
        notes=(
            "N3-A, E4 (gamma) and E2 (template/scale), 319 held-out test perturbations per context, k_ref = 743. "
            "gamma_M3_* = unbiased fraction of reliable held-out gamma-perp energy explained; *_median/_draw_* = "
            "median and 2.5-97.5% over 200 draws; *_boot_* = perturbation bootstrap (2,000) on draw-averaged "
            "terms. template_frac = (M0_E2(k) - M0_E0s) / (M0_E2(k_ref) - M0_E0s) from draw medians; reproduces "
            "the frozen k_T50 (asserted). k = 0 is the zero-shot E0s predictor (gamma_M3 = 0 by construction)."
        ),
    )
    write_source(
        "fig2_gain_fraction",
        gain_fraction(),
        sources=[DRAWS],
        reports=REPORTS,
        build_script="scripts/paper_figures/build_fig2_sources.py",
        plot_script="scripts/paper_figures/plot_fig2.py",
        notes=(
            "G(k) = fraction of the k_ref (743) value reached at k. point = ratio of draw medians (for gamma_perp "
            "at k = 20 this is the preregistered C3 quantity). lo/hi = 2.5/97.5 percentiles over the 200 draws "
            "of the per-draw ratio, with the k_ref value taken from the same cell-split repeat (descriptive). "
            "template_scale uses M0 gain of E2 over the repeat's zero-shot E0s."
        ),
    )


if __name__ == "__main__":
    main()
