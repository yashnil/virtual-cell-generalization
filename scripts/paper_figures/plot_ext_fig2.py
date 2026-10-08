# ruff: noqa: E501 (figure labels and provenance notes are kept as single strings)
"""Extended Data Fig. 2: validity controls for the γ⊥ learning curve (N3-A; N1 comparison)."""

from __future__ import annotations

import numpy as np
import style as st
from _common import load, save


def panel_permutation(ax, p):
    ctx = st.CONTEXTS[::-1]
    for i, c in enumerate(ctx):
        for k, dy in ((20, 0.17), (100, -0.17)):
            x = p[(p.context == c) & (p.k == k)].set_index("variant")
            ok, pm = x.loc["correct anchors"], x.loc["permuted anchors"]
            y = i + dy
            st.herrorbar(
                ax,
                [ok.M3_median],
                [y],
                [ok.lo],
                [ok.hi],
                st.GAMMA,
                marker="D",
                mfc=st.GAMMA if k == 100 else "white",
                ms=2.8,
                lw=0.8,
            )
            st.herrorbar(
                ax, [pm.M3_median], [y], [pm.lo], [pm.hi], st.ZEROSHOT, marker="x", ms=3.2, lw=0.8
            )
        if i % 2 == 0:
            ax.axhspan(i - 0.5, i + 0.5, color="#F7F7F7", zorder=0, lw=0)
    st.ref_line(ax, 0, orient="v", color=st.INK2)
    ax.set_yticks(range(len(ctx)))
    ax.set_yticklabels(ctx)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.5, len(ctx) - 0.5)
    ax.set_xlim(-0.12, 0.47)
    ax.set_xlabel("γ$^{\\perp}$ recovered, E4 (median; 95 % draw interval)")
    ax.set_title("Anchor-permutation null (F3a)", loc="left")
    from matplotlib.lines import Line2D

    ax.legend(
        handles=[
            Line2D(
                [],
                [],
                marker="D",
                color=st.GAMMA,
                mfc="white",
                ls="",
                ms=3,
                label="correct anchors, k = 20",
            ),
            Line2D(
                [], [], marker="D", color=st.GAMMA, ls="", ms=3, label="correct anchors, k = 100"
            ),
            Line2D([], [], marker="x", color=st.ZEROSHOT, ls="", ms=3.5, label="permuted anchors"),
        ],
        loc="lower right",
        fontsize=st.FS_SMALL,
        bbox_to_anchor=(1.02, 0.0),
    )


def panel_shared(ax, f):
    x = f[(f.estimator == "E4")]
    ctx = st.CONTEXTS
    for j, (metric, col, mk, lab) in enumerate(
        (
            ("dM0_shared_minus_split", st.TEMPLATE, "s", "M0 (template included)"),
            ("dM3_shared_minus_split", st.GAMMA, "D", "M3 (γ$^{\\perp}$)"),
        )
    ):
        for i, c in enumerate(ctx):
            for kk, dx in ((5, -0.22), (20, 0.0), (100, 0.22)):
                v = x[(x.context == c) & (x.k == kk)][metric].iloc[0]
                ax.plot(
                    i + dx + (0.06 if j else -0.06),
                    v,
                    mk,
                    color=col,
                    mfc="white" if kk == 5 else col,
                    ms=2.6,
                    alpha=1 if kk != 100 else 0.6,
                    ls="",
                )
        ax.plot([], [], mk, color=col, ms=3, ls="", label=lab)
    st.ref_line(ax, 0, color=st.INK2)
    ax.set_xticks(range(len(ctx)))
    ax.set_xticklabels(ctx, rotation=30, ha="right")
    ax.set_ylabel("Shared − split control\n(difference in energy explained)")
    ax.set_ylim(-0.012, 0.032)
    ax.legend(loc="upper right", fontsize=st.FS_SMALL)
    ax.set_title(
        "Shared-control diagnostic (F3b), E4\nper context, left → right: k = 5 (open), 20, 100 (light)",
        loc="left",
        linespacing=1.2,
    )
    ax.text(
        0.02,
        0.06,
        "M3 identical in every cell (|Δ| < 1e-15)",
        transform=ax.transAxes,
        fontsize=st.FS_SMALL,
        color=st.GAMMA,
    )


def panel_single(ax, h):
    for c in st.CONTEXTS:
        x = h[h.context == c].sort_values("k")
        ax.plot(
            np.log10(x.k), x.E1_M0_gain_median, "-", color=st.tint(st.TEMPLATE, 0.5), lw=st.LW_THIN
        )
    med = h.groupby("k").E1_M0_gain_median.median()
    lo = h.groupby("k").E1_M0_gain_median.min()
    hi = h.groupby("k").E1_M0_gain_median.max()
    ax.fill_between(np.log10(med.index), lo, hi, color=st.TEMPLATE, alpha=0.12, lw=0)
    ax.plot(
        np.log10(med.index),
        med,
        "-s",
        color=st.TEMPLATE,
        lw=st.LW_SUMMARY,
        ms=st.MS,
        mec="white",
        mew=0.5,
    )
    st.ref_line(ax, 0, color=st.INK2)
    ax.set_xticks(np.log10(med.index))
    ax.set_xticklabels([str(k) for k in med.index])
    ax.set_ylim(-2.7, 0.6)
    ax.set_xlabel("Target anchors, k (log scale)")
    ax.set_ylabel("Full-response gain over zero-shot\n(ΔM0, target template from anchors, E1)")
    ax.set_title("Too few anchors are harmful", loc="left")
    ax.text(
        np.log10(48),
        -2.55,
        "thin: contexts (median of 200 draws)\nthick: median; band: range of contexts",
        fontsize=st.FS_SMALL,
        color=st.INK2,
        va="bottom",
        ha="right",
    )


def panel_n1n3(ax, r):
    studies = r.study.unique()
    ctx = st.CONTEXTS[::-1]
    for s_, dy, mfc in ((studies[0], 0.15, "white"), (studies[1], -0.15, st.GAMMA)):
        x = r[r.study == s_].set_index("context")
        for i, c in enumerate(ctx):
            if c in x.index:
                st.herrorbar(
                    ax,
                    [x.loc[c, "ratio_k20_kref"]],
                    [i + dy],
                    [x.loc[c, "lo"]],
                    [x.loc[c, "hi"]],
                    st.GAMMA,
                    marker="D",
                    mfc=mfc,
                    ms=2.8,
                    lw=0.8,
                )
        ax.plot([], [], "D", color=st.GAMMA, mfc=mfc, ls="", ms=3, label=s_)
    ax.axvline(0.25, color=st.INK2, lw=0.6, ls=(0, (2, 2)))
    ax.text(0.26, -0.45, "C3 bar 0.25", fontsize=st.FS_SMALL, color=st.INK2)
    st.ref_line(ax, 0, orient="v", color=st.INK2)
    ax.set_yticks(range(len(ctx)))
    ax.set_yticklabels(ctx)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(-0.35, 0.6)
    ax.set_ylim(-0.5, len(ctx) + 0.3)
    ax.set_xlabel("γ$^{\\perp}$ at k = 20 as fraction of $k_{ref}$ value")
    ax.legend(
        loc="upper left",
        fontsize=st.FS_SMALL,
        bbox_to_anchor=(-0.01, 1.03),
        ncol=2,
        frameon=True,
        facecolor="white",
        edgecolor="none",
        framealpha=1,
    )
    ax.set_title("N1 (confirmatory) and N3-A (replication) agree", loc="left")


def main() -> None:
    fig = st.new_figure(st.DOUBLE, 4.6)
    ax_a = fig.add_axes([0.08, 0.57, 0.36, 0.37])
    ax_b = fig.add_axes([0.6, 0.6, 0.38, 0.33])
    ax_c = fig.add_axes([0.08, 0.09, 0.36, 0.34])
    ax_d = fig.add_axes([0.6, 0.09, 0.38, 0.34])
    panel_permutation(ax_a, load("ext2_permutation"))
    panel_shared(ax_b, load("ext2_shared_control"))
    panel_single(ax_c, load("ext2_single_anchor"))
    panel_n1n3(ax_d, load("ext2_n1_vs_n3"))
    st.panel_label(ax_a, "A", x=-0.06)
    st.panel_label(ax_b, "B", x=-0.08)
    st.panel_label(ax_c, "C", x=-0.06)
    st.panel_label(ax_d, "D", x=-0.08)
    save(
        fig,
        "ext_fig2_n1_controls",
        sources=["ext2_permutation", "ext2_shared_control", "ext2_single_anchor", "ext2_n1_vs_n3"],
        plot_script="scripts/paper_figures/plot_ext_fig2.py",
        upstream=["outputs/n3/n3a", "outputs/n1_n4/n1"],
    )


if __name__ == "__main__":
    main()
