"""Figure 2 — corrected beta and gamma shares across all 21 frozen robustness
variants."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from virtual_cell.visualization import common, style

SHORT = {
    "shared control mean": "shared ctrl",
    "independent control split": "split ctrl",
    "mean of log1p": "mean log1p",
    "log of mean CP10K": "log mean",
    "all 6,640 genes": "6,640 genes",
    "4,000 control HVGs": "4,000 HVGs",
    "2,000 control HVGs": "2,000 HVGs",
}
MARKERS = {"all 6,640 genes": "o", "4,000 control HVGs": "s", "2,000 control HVGs": "^"}


def main() -> None:
    style.apply()
    t = common.load_source("fig2_decomposition_robustness")
    t["condition"] = t.control_scheme + "|" + t.aggregation
    conds = list(dict.fromkeys(t.condition))
    feats = list(MARKERS)
    fig, axes = plt.subplots(1, 2, figsize=(style.FULL_WIDTH, 3.4), sharey=True)
    ticks, labels = [], []
    for ax, col, comp, letter in (
        (axes[0], "beta_corrected", "beta", "A"),
        (axes[1], "gamma_corrected", "gamma", "B"),
    ):
        c = style.COMPONENTS[comp]
        lo, hi = t[col].min(), t[col].max()
        ax.axvspan(lo, hi, color=c["color"], alpha=0.13, lw=0)
        y, ticks, labels = 0.0, [], []
        for cond in conds:
            scheme, agg = cond.split("|")
            for f in feats:
                g = t[(t.condition == cond) & (t.feature_space == f)]
                if not len(g):
                    continue
                vals = g[col].to_numpy()
                canon = bool(g.canonical.any())
                ax.plot([vals.min(), vals.max()], [y, y], color=c["color"], lw=2.5)
                ax.plot(
                    np.median(vals),
                    y,
                    marker=MARKERS[f],
                    ms=6,
                    ls="",
                    mfc="black" if canon else "white",
                    mec="black" if canon else c["color"],
                )
                ticks.append(y)
                seeds = f" ({len(g)} seeds)" if len(g) > 1 else ""
                labels.append(f"{SHORT[scheme]} · {SHORT[agg]} · {SHORT[f]}{seeds}")
                y += 1
            y += 0.6
        canon_v = float(t.loc[t.canonical, col].iloc[0])
        style.reference_line(ax, canon_v, axis="x", lw=0.8)
        ax.set_title(
            f"{c['label']}\n{lo:.2f}–{hi:.2f}% over 21 variants", fontsize=8.5, loc="center"
        )
        ax.set_xlabel("% of response energy\n(noise-corrected)")
        ax.set_xlim(lo - 2.0, hi + 2.0)
        style.panel_label(ax, letter)
    axes[0].set_yticks(ticks, labels, fontsize=7)
    axes[0].invert_yaxis()
    fig.text(
        0.01,
        0.005,
        "Filled marker and dashed line = canonical v1. Shaded band = full range. Seed "
        "ranges (<0.02 pp) are narrower than the markers. x-axes span ±2 pp around the range.",
        fontsize=6.5,
        color=style.GREY,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    common.save_figure(fig, "fig2_decomposition_robustness")


if __name__ == "__main__":
    main()
