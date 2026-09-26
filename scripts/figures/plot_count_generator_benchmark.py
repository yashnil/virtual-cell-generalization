"""Figure 8 — count generators scored with the six local vcc2026 metrics
against real held-out K562 cells."""

from __future__ import annotations

import matplotlib.pyplot as plt

from virtual_cell.visualization import common, style


def main() -> None:
    style.apply()
    t = common.load_source("fig8_count_generator_benchmark")
    metrics = list(dict.fromkeys(t.metric))
    gens = ["G0_control_resample", "G1_transport", "G2_count_model"]
    fig, axes = plt.subplots(2, 3, figsize=(style.FULL_WIDTH, 3.9))
    for ax, m in zip(axes.ravel(), metrics, strict=True):
        sub = t[t.metric == m].set_index("generator")
        better = sub.better.iloc[0]
        for i, g in enumerate(gens):
            meta = style.GENERATORS[g]
            v = float(sub.loc[g, "value"])
            ax.barh(
                i,
                v,
                color=meta["color"],
                hatch=meta["hatch"],
                edgecolor="black",
                lw=0.6,
                height=0.65,
            )
            ax.text(v, i, f" {v:.3f}", va="center", fontsize=7)
        ax.set_yticks(range(3), ["G0", "G1", "G2"])
        ax.invert_yaxis()
        arrow = "→ higher is better" if better == "higher" else "← lower is better"
        ax.set_title(f"{sub.metric_label.iloc[0]}", fontsize=8)
        ax.set_xlabel(arrow, fontsize=7, color=style.GREY)
        vmax = float(sub.value.max())
        ax.set_xlim(0, vmax * 1.35)
        if m == "pds_cosine":
            style.reference_line(ax, 0.5, axis="x", color=style.GREY)
        if m == "expr_mse_unbiased_capped_norm":
            style.reference_line(ax, 1.0, axis="x", color=style.VERMILION)
    handles = [
        plt.Rectangle(
            (0, 0),
            1,
            1,
            fc=style.GENERATORS[g]["color"],
            hatch=style.GENERATORS[g]["hatch"],
            ec="black",
        )
        for g in gens
    ]
    fig.legend(
        handles,
        [style.GENERATORS[g]["label"] for g in gens],
        loc="lower center",
        ncol=3,
        fontsize=7.5,
        bbox_to_anchor=(0.5, 0.0),
    )
    fig.suptitle(
        "Count generators vs real held-out K562 cells (raw local metrics; 300 perts × 400 cells)",
        fontsize=9,
        fontweight="bold",
        x=0.01,
        ha="left",
    )
    fig.text(
        0.01,
        0.075,
        "Dashed lines: pds_cosine = 0.5, the value any constant prediction "
        "gets; expression error = 1.0 (above it, squared error exceeds the effect energy; "
        "no generator gets below it).",
        fontsize=6.3,
        color=style.GREY,
    )
    fig.tight_layout(rect=(0, 0.11, 1, 0.97))
    common.save_figure(fig, "fig8_count_generator_benchmark")


if __name__ == "__main__":
    main()
