"""Figure 1 — response decomposition (uncorrected vs noise-corrected) and
per-component reproducibility, from the frozen four-context analysis."""

from __future__ import annotations

import matplotlib.pyplot as plt

from virtual_cell.visualization import common, style


def main() -> None:
    style.apply()
    t = common.load_source("fig1_decomposition")
    fig, (ax, ax2) = plt.subplots(
        1, 2, figsize=(style.FULL_WIDTH, 3.0), gridspec_kw={"width_ratios": [1.5, 1]}
    )

    views = ["uncorrected", "noise-corrected"]
    order = ["template", "beta", "gamma", "noise"]
    for y, view in enumerate(views):
        left = 0.0
        sub = t[t.view == view].set_index("component")
        for comp in order:
            if comp not in sub.index:
                continue
            v = float(sub.loc[comp, "percent"])
            c = style.COMPONENTS[comp]
            ax.barh(
                y,
                v,
                left=left,
                color=c["color"],
                hatch=c["hatch"],
                edgecolor="black",
                lw=0.6,
                height=0.6,
                label=c["label"] if y == 1 else None,
            )
            ax.text(
                left + v / 2,
                y,
                f"{v:.1f}%",
                ha="center",
                va="center",
                fontsize=8,
                bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none", "alpha": 0.85},
            )
            left += v
    ax.set_yticks([0, 1], ["uncorrected", "noise-corrected\n(50 split halves)"])
    ax.set_xlim(0, 100)
    ax.set_xlabel("share of response energy (%)")
    ax.invert_yaxis()
    ax.set_title(
        "Share of response energy\n(4 contexts × 1,264 perturbations × 6,640 genes)", fontsize=9
    )
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=2, fontsize=7.5)
    style.panel_label(ax, "A")

    rep = t[t.view == "reproducibility"].set_index("component")
    comps = [
        ("mu", "global template μ", style.COMPONENTS["template"]),
        ("alpha", "context template α", style.COMPONENTS["template"]),
        ("beta", "conserved effect β", style.COMPONENTS["beta"]),
        ("gamma", "context-specific\ninteraction γ", style.COMPONENTS["gamma"]),
    ]
    for i, (key, _label, c) in enumerate(comps):
        v = float(rep.loc[key, "percent"])
        ax2.barh(i, v, color=c["color"], hatch=c["hatch"], edgecolor="black", lw=0.6, height=0.6)
        ax2.text(v + 1.5, i, f"{v:.1f}%", va="center", fontsize=8)
    ax2.set_yticks(range(len(comps)), [c[1] for c in comps])
    ax2.invert_yaxis()
    ax2.set_xlim(0, 115)
    ax2.set_xticks([0, 25, 50, 75, 100])
    style.reference_line(ax2, 100, axis="x", color=style.GREY)
    ax2.set_xlabel(
        "reproducible share of the component (%)\n(cross-half signal / raw sum of squares)"
    )
    ax2.set_title("Reproducibility of each component", fontsize=9)
    style.panel_label(ax2, "B")

    fig.tight_layout()
    common.save_figure(fig, "fig1_decomposition")


if __name__ == "__main__":
    main()
