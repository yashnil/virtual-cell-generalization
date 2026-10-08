# ruff: noqa: E501 (figure labels and provenance notes are kept as single strings)
"""Extended Data Fig. 3: correction of the source-agreement estimate and its demotion to a baseline (N4)."""

from __future__ import annotations

import style as st
from _common import load, save

BASE = {
    "b1_magnitude": "source magnitude",
    "b2_source_reliability": "source reliability",
    "b3_source_reliable_energy": "source reliable energy",
    "b6_noise_agreement_ceiling": "noise-implied ceiling",
}


def panel_before_after(ax, h):
    # label positions repelled so that close values do not collide
    order = h.sort_values("rho_per_repeat").reset_index(drop=True)
    ypos, last = {}, -1.0
    for r in order.itertuples():
        last = max(r.rho_per_repeat, last + 0.045)
        ypos[r.context] = last
    for r in h.itertuples():
        ax.plot(
            [0, 1], [r.rho_averaged_halves, r.rho_per_repeat], "-", color=st.INK2, lw=st.LW_THIN
        )
        ax.plot(0, r.rho_averaged_halves, "o", color=st.ZEROSHOT, mfc="white", ms=st.MS)
        ax.plot(-0.12, r.rho_frozen_v1_report, "_", color=st.FAILED, ms=6, mew=1.2)
        ax.plot(1, r.rho_per_repeat, "o", color=st.EMPH, ms=st.MS)
        ax.text(1.06, ypos[r.context], r.context, fontsize=st.FS_SMALL, va="center")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["halves averaged\nover repeats (v1)", "per-repeat\nunbiased (N4)"])
    ax.set_xlim(-0.3, 1.45)
    ax.set_ylim(0, 0.9)
    ax.set_ylabel("Spearman ρ (agreement,\nheld-out transfer quality)")
    ax.text(
        -0.27,
        0.86,
        "— frozen v1 report (0.55–0.79)",
        fontsize=st.FS_SMALL,
        color=st.INK2,
        va="top",
    )
    ax.set_title("Exploratory: v1 estimator inflated", loc="left")


def panel_partial(ax, t2):
    studies = list(t2.study.unique())
    ctx = st.CONTEXTS[::-1]
    for s_, dy, mfc in ((studies[0], 0.15, "white"), (studies[1], -0.15, st.EMPH)):
        x = t2[t2.study == s_].set_index("context")
        for i, c in enumerate(ctx):
            if c in x.index:
                st.herrorbar(
                    ax,
                    [x.loc[c, "partial_spearman"]],
                    [i + dy],
                    [x.loc[c, "lo"]],
                    [x.loc[c, "hi"]],
                    st.EMPH,
                    marker="o",
                    mfc=mfc,
                    ms=2.8,
                    lw=0.8,
                )
        ax.plot([], [], "o", color=st.EMPH, mfc=mfc, ls="", ms=3, label=s_)
    st.ref_line(ax, 0, orient="v", color=st.INK2)
    ax.set_yticks(range(len(ctx)))
    ax.set_yticklabels(ctx)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(-0.2, 0.4)
    ax.set_ylim(-0.5, len(ctx) + 0.2)
    ax.set_xlabel("Partial Spearman ρ given b1–b5 (95 % CI)")
    ax.legend(loc="upper right", fontsize=st.FS_SMALL, bbox_to_anchor=(1.0, 1.04))
    ax.set_title("Incremental information is modest (T2)", loc="left")


def panel_t1(fig, gs, t1):
    ctx = st.CONTEXTS[::-1]
    axes = []
    for j, (b, lab) in enumerate(BASE.items()):
        ax = fig.add_subplot(gs[0, j])
        axes.append(ax)
        x = t1[t1.statistic == b].set_index("context")
        for i, c in enumerate(ctx):
            r = x.loc[c]
            sig_worse = r.diff_hi < 0
            col = st.EMPH if (r.diff_lo > 0 or sig_worse) else st.ZEROSHOT
            st.herrorbar(
                ax,
                [r.diff_vs_agreement],
                [i],
                [r.diff_lo],
                [r.diff_hi],
                col,
                marker="v" if sig_worse else "o",
                mfc="white" if sig_worse else None,
                ms=2.8,
                lw=0.8,
            )
        st.ref_line(ax, 0, orient="v", color=st.INK2)
        ax.set_yticks(range(len(ctx)))
        ax.set_yticklabels(ctx if j == 0 else [])
        ax.tick_params(axis="y", length=0)
        ax.set_xlim(-0.16, 0.26)
        ax.set_xticks([-0.1, 0, 0.1, 0.2])
        ax.set_ylim(-0.6, len(ctx) - 0.4)
        ax.set_title(f"vs {lab}", fontsize=st.FS_SMALL, pad=2)
    return axes


def main() -> None:
    fig = st.new_figure(st.DOUBLE, 4.3)
    ax_a = fig.add_axes([0.09, 0.58, 0.32, 0.35])
    ax_b = fig.add_axes([0.6, 0.58, 0.37, 0.35])
    panel_before_after(ax_a, load("ext3_halves_bias"))
    panel_partial(ax_b, load("ext3_partial"))
    gs = fig.add_gridspec(1, 4, left=0.09, right=0.98, top=0.37, bottom=0.12, wspace=0.12)
    axes_c = panel_t1(fig, gs, load("ext3_t1"))
    fig.text(
        0.535,
        0.03,
        "ρ(agreement) − |ρ(baseline)|, N3 six folds (paired bootstrap 95 % CI)",
        ha="center",
        fontsize=st.FS,
    )
    fig.text(
        0.09,
        0.44,
        "Agreement vs simple baselines (T1): ● agreement better · ▽ agreement worse · grey: not distinguishable",
        fontsize=st.FS_SMALL,
        color=st.INK2,
    )
    st.panel_label(ax_a, "A", x=-0.07)
    st.panel_label(ax_b, "B", x=-0.07)
    st.panel_label(axes_c[0], "C", x=-0.07, y=1.12)
    save(
        fig,
        "ext_fig3_agreement_correction",
        sources=["ext3_halves_bias", "ext3_partial", "ext3_t1"],
        plot_script="scripts/paper_figures/plot_ext_fig3.py",
        upstream=["outputs/n1_n4/n4", "outputs/n3/n4"],
    )


if __name__ == "__main__":
    main()
