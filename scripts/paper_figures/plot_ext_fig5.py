# ruff: noqa: E501 (figure labels and provenance notes are kept as single strings)
"""Extended Data Fig. 5: full N3-B source-count analysis — preregistered fixed frame and exploratory own frame.

A: preregistered fixed-reference R_γ (k = 0, 20) vs the exploratory own-frame R_γ (k = 20), per target.
B: frame-free template-removed full response, zero-shot vs 20 anchors, per target.
C: the preregistered Q1 quantity Δ(20) = R(20, 5) − R(20, 3) in the fixed frame vs the same pairing in the own frame.
"""

from __future__ import annotations

import numpy as np
import style as st
from _common import load, save
from matplotlib.lines import Line2D

M = [1, 2, 3, 4, 5]
KS = [5, 10, 20, 50, 100, 200, 743]
FS, FSS = st.FS, st.FS_SMALL


def _col(x, col, k):
    return np.array([x.loc[(k, m), col] for m in M])


def panel_frames(fig, gs, r):
    axes = []
    for i, c in enumerate(st.CONTEXTS):
        ax = fig.add_subplot(gs[0, i])
        axes.append(ax)
        x = r[r.context == c].set_index(["k", "m"])
        ref0, ref20 = _col(x, "R_gamma", 0), _col(x, "R_gamma", 20)
        lo20, hi20, own20 = (
            _col(x, "R_gamma_lo", 20),
            _col(x, "R_gamma_hi", 20),
            _col(x, "R_gamma_own", 20),
        )
        mm = np.array(M)
        on = mm >= 2
        ax.fill_between(mm[on], lo20[on], hi20[on], color=st.EMPH, alpha=0.1, lw=0)
        ax.plot(mm[on], ref0[on], "-o", color=st.ZEROSHOT, mfc="white", ms=2.8, lw=st.LW)
        ax.plot(mm[on], ref20[on], "-o", color=st.EMPH, ms=2.8, lw=st.LW)
        ax.plot(mm, own20, "-D", color=st.GAMMA, ms=2.6, lw=st.LW)
        ax.text(
            0.75,
            -0.155,
            f"m = 1, k = 0: {ref0[0]:.2f}".replace("-", "−"),
            fontsize=FSS,
            color=st.INK2,
            va="bottom",
        )
        st.ref_line(ax, 0, color=st.INK2)
        ax.set_ylim(-0.16, 0.2)
        ax.set_xlim(0.6, 5.4)
        ax.set_xticks(M)
        if i:
            ax.set_yticklabels([])
        ax.set_title(c, loc="left", fontsize=FS, fontweight="bold", pad=1.5)
    axes[0].set_ylabel("γ$^{\\perp}$ recovered")
    return axes


def panel_full_response(fig, gs, t):
    axes = []
    ylo, yhi = -0.16, 0.5
    for i, c in enumerate(st.CONTEXTS):
        ax = fig.add_subplot(gs[0, i])
        axes.append(ax)
        x = t[t.context == c]
        z = x[x.k == 0].set_index("m").loc[M]
        a = x[x.k == 20].set_index("m").loc[M]
        mm = np.array(M)
        on = (z.R_full >= ylo).to_numpy()
        ax.plot(mm[on], z.R_full[on], "-", color=st.ZEROSHOT, lw=st.LW)
        st.errorbar(
            ax,
            mm[on],
            z.R_full[on],
            z.R_full_k0_rep_min[on],
            z.R_full_k0_rep_max[on],
            st.ZEROSHOT,
            mfc="white",
            ms=2.6,
        )
        if (~on).any():
            ax.plot(
                mm[~on],
                np.full((~on).sum(), ylo + 0.01),
                "v",
                color=st.ZEROSHOT,
                ms=2.8,
                clip_on=False,
            )
            ax.text(
                1.25,
                ylo + 0.03,
                f"{z.R_full[~on].iloc[0]:.2f}".replace("-", "−"),
                fontsize=FSS,
                color=st.INK2,
                va="bottom",
            )
        ax.plot(M, a.R_full, "-", color=st.EMPH, lw=st.LW)
        st.errorbar(ax, M, a.R_full, a.R_full_lo, a.R_full_hi, st.EMPH, ms=2.6)
        st.ref_line(ax, 0)
        ax.set_ylim(ylo, yhi)
        ax.set_xlim(0.6, 5.4)
        ax.set_xticks(M)
        ax.set_yticks([0, 0.2, 0.4])
        if i:
            ax.set_yticklabels([])
        ax.set_title(c, loc="left", fontsize=FS, fontweight="bold", pad=1.5)
    return axes


def panel_delta(ax, d):
    ctx = st.CONTEXTS[::-1]
    for i, c in enumerate(ctx):
        if i % 2 == 0:
            ax.axhspan(i - 0.5, i + 0.5, color="#F7F7F7", zorder=0, lw=0)
        x = d[d.context == c].set_index("frame")
        f, o = x.loc["preregistered fixed frame"], x.loc["own frame (exploratory)"]
        st.herrorbar(
            ax,
            [f.delta20_m5_minus_m3],
            [i + 0.16],
            [f.lo],
            [f.hi],
            st.EMPH,
            marker="o",
            ms=3,
            lw=0.9,
        )
        st.herrorbar(
            ax,
            [o.delta20_m5_minus_m3],
            [i - 0.16],
            [o.lo],
            [o.hi],
            st.GAMMA,
            marker="D",
            ms=2.8,
            lw=0.9,
        )
    ax.axvline(0.03, color=st.INK2, lw=0.6, ls=(0, (2, 2)))
    ax.text(0.033, -0.45, "Q1 bar 0.03", fontsize=FSS, color=st.INK2, va="bottom")
    st.ref_line(ax, 0, orient="v", color=st.INK2)
    ax.set_yticks(range(len(ctx)))
    ax.set_yticklabels(ctx)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.5, len(ctx) + 1.6)
    ax.set_xticks([0, 0.05, 0.1, 0.15])
    ax.set_xlim(-0.05, 0.17)
    ax.set_xlabel(
        "Δ(20) = R(k = 20, m = 5) − R(k = 20, m = 3)\n(paired over 40 draws; 95 % draw interval)"
    )
    ax.set_title(
        "Preregistered source-count test (Q1) and its own-frame counterpart",
        loc="left",
        fontweight="bold",
    )
    return ax


def panel_gain(fig, gs, t):
    axes = []
    ctx = st.CONTEXTS[::-1]
    for j, k in enumerate([20, 50]):
        ax = fig.add_subplot(gs[0, j])
        axes.append(ax)
        x = t[t.k == k]
        for i, c in enumerate(ctx):
            if i % 2 == 0:
                ax.axhspan(i - 0.5, i + 0.5, color="#F7F7F7", zorder=0, lw=0)
            v = x[x.context == c].set_index("m").loc[M]
            a, b = float(v.R_gamma_own.loc[1]), float(v.R_gamma_own.loc[5])
            ax.plot(
                [a, b],
                [i, i],
                "-",
                color=st.tint(st.GAMMA, 0.45),
                lw=1.8,
                zorder=1,
                solid_capstyle="butt",
            )
            ax.plot(
                v.R_gamma_own.loc[[2, 3, 4]],
                [i] * 3,
                "|",
                color=st.GAMMA,
                ms=4.5,
                mew=0.8,
                zorder=2,
            )
            st.herrorbar(
                ax,
                [a],
                [i + 0.18],
                [v.own_lo.loc[1]],
                [v.own_hi.loc[1]],
                st.tint(st.GAMMA, 0.7),
                marker="D",
                mfc="white",
                ms=2.6,
                lw=0.5,
            )
            st.herrorbar(
                ax,
                [b],
                [i - 0.18],
                [v.own_lo.loc[5]],
                [v.own_hi.loc[5]],
                st.GAMMA,
                marker="D",
                ms=2.6,
                lw=0.5,
            )
            ax.text(
                1.03,
                i,
                f"{b - a:+.3f}".replace("-", "−"),
                fontsize=FSS,
                va="center",
                ha="left",
                color=st.INK,
                transform=ax.get_yaxis_transform(),
            )
        st.ref_line(ax, 0, orient="v")
        ax.set_xlim(-0.09, 0.28)
        ax.set_xticks([0, 0.1, 0.2])
        ax.set_ylim(-0.5, len(ctx) - 0.5)
        ax.set_yticks(range(len(ctx)))
        ax.set_yticklabels(ctx if j == 0 else [])
        ax.tick_params(axis="y", length=0)
        ax.set_title(f"k = {k}", fontsize=FS, pad=2, loc="left")
        ax.text(
            1.03,
            len(ctx) - 0.5,
            "Δ\nm 1→5",
            fontsize=FSS,
            ha="left",
            va="bottom",
            color=st.INK2,
            transform=ax.get_yaxis_transform(),
            linespacing=1.0,
        )
    return axes


def panel_heatmaps(fig, gs, t, cax):
    axes = []
    for i, c in enumerate(st.CONTEXTS):
        ax = fig.add_subplot(gs[0, i])
        axes.append(ax)
        x = (
            t[(t.context == c) & (t.k > 0)]
            .pivot(index="m", columns="k", values="R_gamma_own")
            .loc[M, KS]
        )
        im = ax.imshow(
            x.to_numpy(),
            origin="lower",
            aspect="auto",
            cmap=st.GAMMA_CMAP,
            vmin=0,
            vmax=0.42,
            interpolation="nearest",
        )
        ax.set_xticks(range(len(KS)))
        ax.set_xticklabels([str(k) for k in KS], rotation=90)
        ax.set_yticks(range(len(M)))
        ax.set_yticklabels(M if i == 0 else [])
        ax.tick_params(length=0, pad=1.5)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_title(c, fontsize=FS, fontweight="bold", pad=2)
    axes[0].set_ylabel("Source contexts, m")
    cb = fig.colorbar(im, cax=cax, extend="min")
    cb.set_label("γ$^{\\perp}$ learned\n(own frame)", fontsize=FSS)
    cb.ax.tick_params(labelsize=FSS, length=2)
    cb.outline.set_linewidth(0.5)
    return axes


def main() -> None:
    r = load("ext5_reference_frame")
    t = load("fig3_source_count")
    d = load("ext5_delta20")
    fig = st.new_figure(st.DOUBLE, 222 / 25.4)

    # Row 1: A (preregistered Q1, primary) | B (exploratory own-frame calibration gain)
    ax_a = fig.add_axes([0.1, 0.785, 0.4, 0.165])
    panel_delta(ax_a, d)
    ax_a.legend(
        handles=[
            Line2D(
                [],
                [],
                color=st.EMPH,
                marker="o",
                ls="",
                ms=3,
                label="fixed 5-source frame (preregistered Q1)",
            ),
            Line2D(
                [],
                [],
                color=st.GAMMA,
                marker="D",
                ls="",
                ms=3,
                label="subset's own frame (exploratory)",
            ),
        ],
        loc="upper right",
        ncol=1,
        frameon=True,
        facecolor="white",
        edgecolor="none",
        framealpha=1,
    )
    gs_b = fig.add_gridspec(1, 2, left=0.66, right=0.915, top=0.935, bottom=0.785, wspace=0.42)
    axes_b = panel_gain(fig, gs_b, t)
    p0, p1 = axes_b[0].get_position(), axes_b[1].get_position()
    fig.text(
        (p0.x0 + p1.x1) / 2,
        p0.y0 - 0.032,
        "γ$^{\\perp}$ from anchors, own frame (exploratory)",
        ha="center",
        va="top",
        fontsize=FS,
    )
    fig.text(
        p0.x0 - 0.06, 0.967, "Own-frame calibration gain vs m (exploratory)", fontsize=FS, va="top"
    )
    fig.legend(
        handles=[
            Line2D(
                [],
                [],
                color=st.tint(st.GAMMA, 0.7),
                marker="D",
                mfc="white",
                ls="",
                ms=3,
                label="m = 1",
            ),
            Line2D([], [], color=st.GAMMA, marker="|", ls="", ms=4.5, label="m = 2–4"),
            Line2D([], [], color=st.GAMMA, marker="D", ls="", ms=3, label="m = 5"),
        ],
        loc="upper center",
        bbox_to_anchor=((p0.x0 + p1.x1) / 2, p0.y0 - 0.052),
        ncol=3,
        columnspacing=0.8,
        handletextpad=0.2,
        handlelength=1.0,
        frameon=False,
    )

    # Row 2: C fixed vs own frame curves
    gs_c = fig.add_gridspec(1, 6, left=0.075, right=0.985, top=0.66, bottom=0.505, wspace=0.12)
    axes_c = panel_frames(fig, gs_c, r)
    fig.text(0.53, 0.475, "Number of source contexts, m", ha="center", fontsize=FS)
    h = [
        Line2D(
            [],
            [],
            color=st.ZEROSHOT,
            marker="o",
            mfc="white",
            ms=3,
            label="fixed frame, k = 0 (preregistered)",
        ),
        Line2D(
            [],
            [],
            color=st.EMPH,
            marker="o",
            ms=3,
            label="fixed frame, k = 20 (band: 95 % draw interval)",
        ),
        Line2D([], [], color=st.GAMMA, marker="D", ms=3, label="own frame, k = 20 (exploratory)"),
    ]
    fig.legend(
        handles=h,
        loc="lower center",
        bbox_to_anchor=(0.53, 0.675),
        ncol=3,
        frameon=False,
        columnspacing=1.2,
    )

    # Row 3: D frame-free full response, zero-shot vs 20 anchors
    gs_d = fig.add_gridspec(1, 6, left=0.075, right=0.985, top=0.39, bottom=0.27, wspace=0.12)
    axes_d = panel_full_response(fig, gs_d, t)
    axes_d[0].set_ylabel("Full response\nexplained (M1)")
    fig.text(0.53, 0.245, "Number of source contexts, m", ha="center", fontsize=FS)
    hd = [
        Line2D(
            [],
            [],
            color=st.ZEROSHOT,
            marker="o",
            mfc="white",
            ms=3,
            label="zero-shot (bars: range, 5 repeats)",
        ),
        Line2D(
            [], [], color=st.EMPH, marker="o", ms=3, label="k = 20 anchors (95 % draw interval)"
        ),
    ]
    fig.legend(
        handles=hd,
        loc="lower center",
        bbox_to_anchor=(0.53, 0.405),
        ncol=2,
        frameon=False,
        columnspacing=1.2,
    )

    # Row 4: E own-frame k × m heatmaps (exploratory)
    gs_e = fig.add_gridspec(1, 6, left=0.075, right=0.9, top=0.185, bottom=0.06, wspace=0.12)
    cax = fig.add_axes([0.915, 0.07, 0.01, 0.105])
    axes_e = panel_heatmaps(fig, gs_e, t, cax)
    fig.text(
        0.49,
        0.008,
        "Measured target perturbations, k   (own-frame γ⊥, exploratory; k = 0 omitted, 0 by construction; grey = below 0)",
        ha="center",
        fontsize=FS,
    )

    st.panel_label(ax_a, "A", x=-0.075, y=1.015)
    st.panel_label(axes_b[0], "B", x=-0.075, y=1.015)
    st.panel_label(axes_c[0], "C", x=-0.05, y=1.04)
    st.panel_label(axes_d[0], "D", x=-0.05, y=1.04)
    st.panel_label(axes_e[0], "E", x=-0.04, y=1.03)
    save(
        fig,
        "ext_fig5_n3_reference_frame",
        sources=["ext5_delta20", "ext5_reference_frame", "fig3_source_count"],
        plot_script="scripts/paper_figures/plot_ext_fig5.py",
        upstream=["outputs/n3/n3b", "outputs/n3/n3_decision.json"],
    )


if __name__ == "__main__":
    main()
