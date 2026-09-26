"""Figure 5 — internal vs external unseen-perturbation generalization, with the
direct-transfer positive control."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from virtual_cell.visualization import common, style


def main() -> None:
    style.apply()
    t = common.load_source("fig5_external_generalization")
    a = t[t.panel == "A"]
    b = t[t.panel == "B"]
    fig, axes = plt.subplots(
        1, 2, figsize=(style.FULL_WIDTH, 3.5), gridspec_kw={"width_ratios": [1.45, 1]}
    )
    ax = axes[0]
    stages = list(dict.fromkeys(a.stage))
    for i, stage in enumerate(stages):
        for arm, meta, dx in (("prior", style.PRIOR, -0.13), ("direct", style.DIRECT, 0.13)):
            row = a[(a.stage == stage) & (a.arm == arm)]
            if not len(row):
                continue
            r = row.iloc[0]
            yerr = None
            if np.isfinite(r.lo):
                yerr = [[r.pearson - r.lo], [r.hi - r.pearson]]
            ax.errorbar(
                i + dx,
                r.pearson,
                yerr=yerr,
                fmt=meta["marker"],
                color=meta["color"],
                mec="black",
                ms=8,
                capsize=3,
                lw=1,
            )
            ax.annotate(
                f"{r.pearson:+.3f}\nn={int(r.n)}",
                (i + dx, r.pearson),
                xytext=(0, 10),
                textcoords="offset points",
                ha="center",
                fontsize=6.5,
            )
    style.reference_line(ax, 0)
    ax.axvspan(1.5, len(stages) - 0.5, color=style.LIGHT_GREY, alpha=0.35, lw=0)
    ax.text(len(stages) - 0.55, -0.08, "external", ha="right", fontsize=7.5, color=style.GREY)
    ax.text(-0.45, -0.08, "internal (4 screens)", ha="left", fontsize=7.5, color=style.GREY)
    labels = ["held-out\nperturbation", "held-out pert.\n+ context", "arch1", "Feng\n(pooled)"]
    ax.set_xticks(range(len(stages)), labels)
    ax.set_ylim(-0.1, 0.66)
    ax.set_ylabel("median per-perturbation Pearson\n(predicted vs held-out conserved effect)")
    ax.set_title("Priors collapse externally; direct transfer does not", fontsize=9)
    ax.plot(
        [],
        [],
        style.PRIOR["marker"],
        color=style.PRIOR["color"],
        mec="black",
        label="unseen-perturbation prior (STRING k-NN)",
    )
    ax.plot(
        [],
        [],
        style.DIRECT["marker"],
        color=style.DIRECT["color"],
        mec="black",
        label="direct transfer (measured in sources)",
    )
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.93), fontsize=7)
    style.panel_label(ax, "A")

    ax = axes[1]
    pr = b[b.arm == "prior"].set_index("stage").pearson
    di = b[b.arm == "direct"].set_index("stage").pearson
    order = di.sort_values().index
    y = np.arange(len(order))
    for yi, line in zip(y, order, strict=True):
        ax.plot([pr[line], di[line]], [yi, yi], color=style.LIGHT_GREY, lw=1, zorder=0)
    ax.plot(pr[order], y, style.PRIOR["marker"], color=style.PRIOR["color"], mec="black", ms=5)
    ax.plot(di[order], y, style.DIRECT["marker"], color=style.DIRECT["color"], mec="black", ms=5)
    style.reference_line(ax, 0, axis="x")
    ax.set_yticks(y, order, fontsize=6.5)
    ax.set_xlabel("median per-perturbation Pearson")
    n_pos = int((di > 0).sum())
    ax.set_title(f"Feng, per iPSC line: direct > 0 in {n_pos}/{len(di)}", fontsize=9)
    style.panel_label(ax, "B")
    fig.text(
        0.01,
        0.005,
        "Same statistic throughout, but perturbation sets and gene axes differ between stages. "
        "Internal bars: min–max over folds (× held-out contexts for the second).",
        fontsize=6.3,
        color=style.GREY,
    )
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    common.save_figure(fig, "fig5_external_generalization")


if __name__ == "__main__":
    main()
