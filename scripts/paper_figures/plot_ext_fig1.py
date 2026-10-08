# ruff: noqa: E501 (figure labels and provenance notes are kept as single strings)
"""Extended Data Fig. 1: decomposition robustness (21 predeclared variants; controlled depth experiment)."""

from __future__ import annotations

import numpy as np
import style as st
from _common import load, save

RNG = np.random.default_rng(0)  # vertical jitter only


def strip(ax, rows, v, colours, markers, xlabel, xlim, title):
    n = len(rows)
    for i, (col_name, _label) in enumerate(rows):
        y = n - 1 - i
        vals = v[col_name].to_numpy()
        jit = RNG.uniform(-0.18, 0.18, size=len(vals))
        c = colours[i]
        ax.plot(
            vals[~v.canonical],
            y + jit[~v.canonical],
            markers[i],
            color=st.tint(c, 0.6),
            mfc="white",
            ms=2.6,
            mew=0.7,
            ls="",
        )
        can = vals[v.canonical.to_numpy()][0]
        ax.plot(can, y, markers[i], color=c, ms=st.MS + 1, mec="white", mew=0.5, zorder=4)
        ax.text(
            vals.max() + (xlim[1] - xlim[0]) * 0.02,
            y,
            f"{vals.min():.1f}–{vals.max():.1f}",
            fontsize=st.FS_SMALL,
            va="center",
            color=st.INK2,
        )
    ax.set_yticks(range(n))
    ax.set_yticklabels([lab for _, lab in rows][::-1])
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xlabel(xlabel)
    ax.set_title(title, loc="left")


def main() -> None:
    v = load("ext1_variants")
    d = load("ext1_depth")
    fig = st.new_figure(st.DOUBLE, 3.9)
    ax_a = fig.add_axes([0.15, 0.58, 0.33, 0.31])
    ax_b = fig.add_axes([0.66, 0.58, 0.28, 0.31])
    ax_c = fig.add_axes([0.15, 0.1, 0.22, 0.34])
    ax_d = fig.add_axes([0.6, 0.1, 0.34, 0.34])

    strip(
        ax_a,
        [
            ("template_corrected", "Template (μ + α)"),
            ("beta_corrected", "Conserved β"),
            ("gamma_corrected", "Interaction γ"),
            ("noise", "Noise"),
        ],
        v,
        [st.TEMPLATE, st.BETA, st.GAMMA, "#9A9A9A"],
        ["s", "o", "D", "o"],
        "Share of total response energy (%), noise-corrected",
        (0, 45),
        "Energy shares, 21 variants",
    )
    strip(
        ax_b,
        [("rep_beta", "Conserved β"), ("rep_gamma", "Interaction γ")],
        v,
        [st.BETA, st.GAMMA],
        ["o", "D"],
        "Reproducible across split halves (%)",
        (0, 100),
        "Reproducibility, 21 variants",
    )

    for r in v.itertuples():
        col = st.GAMMA if r.canonical else st.tint(st.GAMMA, 0.45)
        ax_c.plot(
            [0, 1],
            [r.gamma_uncorrected, r.gamma_corrected],
            "-",
            color=col,
            lw=st.LW if r.canonical else st.LW_THIN,
            zorder=3 if r.canonical else 2,
        )
        ax_c.plot([0], [r.gamma_uncorrected], "D", color=col, mfc="white", ms=2.6)
        ax_c.plot([1], [r.gamma_corrected], "D", color=col, ms=2.6)
    ax_c.set_xticks([0, 1])
    ax_c.set_xticklabels(["uncorrected", "noise-corrected"])
    ax_c.set_xlim(-0.3, 1.3)
    ax_c.set_ylim(0, 50)
    ax_c.set_ylabel("γ share of energy (%)")
    ax_c.set_title("Noise correction halves γ in every variant", loc="left")

    ax_d.fill_between(d.cells_per_half, d.ci_lo, d.ci_hi, color=st.EMPH, alpha=0.15, lw=0)
    ax_d.plot(
        d.cells_per_half, d.median_reliability, "-o", color=st.EMPH, ms=st.MS, mec="white", mew=0.5
    )
    for r in d.itertuples():
        ax_d.text(
            r.cells_per_half,
            r.ci_hi + 0.02,
            f"{r.median_reliability:.2f}",
            fontsize=st.FS_SMALL,
            ha="center",
            va="bottom",
        )
    ax_d.set_xscale("log")
    ax_d.set_xticks(d.cells_per_half)
    ax_d.set_xticklabels([str(int(x)) for x in d.cells_per_half])
    ax_d.minorticks_off()
    ax_d.set_ylim(0, 0.55)
    ax_d.set_xlabel("Cells per split half (same 643 context × perturbation pairs)")
    ax_d.set_ylabel("Median split-half reliability")
    ax_d.set_title("Reliability rises with depth (95 % CI over pairs)", loc="left")

    from matplotlib.lines import Line2D

    h = [
        Line2D(
            [], [], marker="o", color="#777777", ls="", ms=st.MS + 1, label="canonical analysis"
        ),
        Line2D(
            [],
            [],
            marker="o",
            color="#AAAAAA",
            mfc="white",
            ls="",
            ms=2.6,
            label="other predeclared variant",
        ),
    ]
    fig.legend(handles=h, loc="upper right", bbox_to_anchor=(0.99, 0.995), ncol=2, frameon=False)
    st.panel_label(ax_a, "A", x=-0.12)
    st.panel_label(ax_b, "B", x=-0.1)
    st.panel_label(ax_c, "C", x=-0.12)
    st.panel_label(ax_d, "D", x=-0.08)
    save(
        fig,
        "ext_fig1_robustness",
        sources=["ext1_variants", "ext1_depth"],
        plot_script="scripts/paper_figures/plot_ext_fig1.py",
        upstream=["outputs/four_context_sensitivity"],
    )


if __name__ == "__main__":
    main()
