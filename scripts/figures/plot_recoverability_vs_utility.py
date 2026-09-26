"""Figure 4 — recovering the interaction does not make the response prediction
better (pathway residual model v2, clean beta-free target)."""

from __future__ import annotations

import matplotlib.pyplot as plt

from virtual_cell.visualization import common, style

OFFSETS = {"K562": (0, -13), "Jurkat": (0, -13), "RPE1": (9, 0), "HepG2": (0, -13)}


def main() -> None:
    style.apply()
    t = common.load_source("fig4_recoverability_vs_utility")
    fig, ax = plt.subplots(figsize=(style.FULL_WIDTH * 0.72, 3.6))
    style.reference_line(ax, 0)
    for r in t.itertuples():
        c = style.CONTEXTS[r.cell_line]
        ax.errorbar(
            r.r_gamma_v2_model,
            r.delta_pearson_selected,
            yerr=[
                [r.delta_pearson_selected - r.delta_pearson_ci_lo],
                [r.delta_pearson_ci_hi - r.delta_pearson_selected],
            ],
            fmt=c["marker"],
            color=c["color"],
            mec="black",
            ms=8,
            capsize=3,
            lw=1,
        )
        ax.plot(
            r.r_gamma_v2_model,
            r.delta_pearson_theory,
            marker=c["marker"],
            ms=7,
            mfc="white",
            mec=c["color"],
            ls="",
        )
        ax.plot(
            [r.r_gamma_v2_model] * 2,
            [r.delta_pearson_selected, r.delta_pearson_theory],
            color=c["color"],
            lw=0.6,
            ls=":",
        )
        ax.annotate(
            f"{r.cell_line} (λ = {r.selected_lambda:g})",
            (r.r_gamma_v2_model, r.delta_pearson_theory),
            xytext=OFFSETS[r.cell_line],
            textcoords="offset points",
            fontsize=8,
            ha="left" if OFFSETS[r.cell_line][0] > 0 else "center",
            va="center",
            color=c["color"],
            fontweight="bold",
        )
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.11, 0.02)
    ax.set_xlabel("interaction recovery: r(learned correction, true γ)")
    ax.set_ylabel("change in response prediction\n(Δ median Pearson vs conserved transfer)")
    ax.set_title("Recoverable is not the same as useful")
    ax.plot([], [], "ko", mec="black", label="inner-selected λ (bootstrap 95% CI)")
    ax.plot([], [], "o", mfc="white", mec="black", label="λ_theory = (3+s)/4 (≈0.88–0.90)")
    ax.legend(loc="lower right", fontsize=7.5)
    ax.text(
        0.99,
        0.97,
        "above 0 = the correction helps",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=7,
        color=style.GREY,
    )
    fig.tight_layout()
    common.save_figure(fig, "fig4_recoverability_vs_utility")


if __name__ == "__main__":
    main()
