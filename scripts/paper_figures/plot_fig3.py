# ruff: noqa: E501 (figure labels and provenance notes are kept as single strings)
"""Figure 3: broader source panels improve zero-shot transfer across most target contexts.

N3-B source-count ablation, frame-free zero-shot quantity only (M1 at k = 0). Panel C re-displays the m = 2 and
m = 5 values of panel B per target; no new estimand. All own-frame (exploratory) analyses are in ED Fig. 5.
"""

from __future__ import annotations

import numpy as np
import style as st
from _common import load, save
from matplotlib.lines import Line2D

KS = [5, 10, 20, 50, 100, 200, 743]
M = [1, 2, 3, 4, 5]
FS, FSS = st.FS, st.FS_SMALL


def repel(vals, gap):
    """Label positions (same order as ``vals``) separated by at least ``gap``."""
    order = np.argsort(vals)
    pos = np.array(vals, dtype=float)
    last = -np.inf
    for i in order:
        pos[i] = max(vals[i], last + gap)
        last = pos[i]
    return pos


def panel_design(fig, rect):
    ax, W, H = st.schematic_axes(fig, rect)
    top = H - 2.0
    bw, bh, gap = 11.0, 3.6, 1.3
    x0 = 2.0
    y_top = top - 6.0
    ax.text(x0 + bw / 2, top, "Sources", fontsize=FS, fontweight="bold", ha="center", va="top")
    for i in range(5):
        used = i < 3
        y = y_top - i * (bh + gap) - bh
        st.box(
            ax,
            x0,
            y,
            bw,
            bh,
            fc=st.tint(st.ZEROSHOT, 0.7) if used else "white",
            ec=st.INK2 if used else "#B0B0B0",
            ls="-" if used else (0, (2, 1.5)),
        )
        ax.text(
            x0 + bw / 2,
            y + bh / 2,
            f"S{i + 1}",
            fontsize=FSS,
            ha="center",
            va="center",
            color=st.INK if used else "#9A9A9A",
        )
    yb0 = y_top - 3 * (bh + gap) + gap
    xb = x0 + bw + 1.2
    ax.plot([xb, xb + 0.8, xb + 0.8, xb], [y_top, y_top, yb0, yb0], color=st.INK, lw=st.STROKE)
    ax.text(xb + 1.5, (y_top + yb0) / 2, "m", fontsize=FS, style="italic", va="center")
    # target
    tx, tw, th = 29.5, 13.5, 15.0
    ty = y_top - th - 4.6
    ax.text(
        tx + tw / 2, top, "Held-out target", fontsize=FS, fontweight="bold", ha="center", va="top"
    )
    st.box(ax, tx, ty + th + 0.9, tw, 2.6, fc=st.tint(st.ZEROSHOT, 0.7), ec=st.INK2, radius=0.0)
    ax.text(tx + tw / 2, ty + th + 2.2, "controls", fontsize=FSS, ha="center", va="center")
    st.cell_grid(ax, tx, ty, tw, th, 7, 5, st.EMPH, rows_on=[1, 3], seed=5)
    ch = th / 7
    xk = tx + tw + 0.9
    ax.plot(
        [xk, xk + 0.8, xk + 0.8, xk],
        [ty + th - ch, ty + th - ch, ty + th - 4 * ch, ty + th - 4 * ch],
        color=st.INK,
        lw=st.STROKE,
    )
    ax.text(xk + 1.4, ty + th - 2.5 * ch, "k", fontsize=FS, style="italic", va="center")
    st.arrow(ax, (xb + 4.4, ty + th / 2), (tx - 1.2, ty + th / 2))
    ly = ty - 7.0
    ax.text(
        2.0,
        ly,
        "m = 1–5 source contexts (all 31\nsubsets; none selected)",
        fontsize=FSS,
        va="top",
        linespacing=1.15,
    )
    ax.text(
        2.0,
        ly - 6.2,
        "k = 0–743 random target anchors\n(this figure: k = 0; k > 0 in ED Fig. 5)",
        fontsize=FSS,
        va="top",
        linespacing=1.15,
    )
    return ax


def panel_zero_shot(ax, t):
    ylo, yhi = -0.16, 0.42
    z = t[t.k == 0]
    ends, names = [], []
    mm = np.array(M)
    for c in st.CONTEXTS:
        x = z[z.context == c].set_index("m").loc[M]
        on = x.R_full.to_numpy() >= ylo
        ax.plot(mm[on], x.R_full[on], "-", color=st.ZEROSHOT, lw=st.LW, zorder=2)
        st.errorbar(
            ax,
            mm[on],
            x.R_full[on],
            x.R_full_k0_rep_min[on],
            x.R_full_k0_rep_max[on],
            st.ZEROSHOT,
            mfc="white",
            ms=2.8,
        )
        ends.append(float(x.R_full.iloc[-1]))
        names.append(c)
    med = z[z.m >= 2].groupby("m").R_full.median()
    ax.plot(
        med.index,
        med.to_numpy(),
        "-",
        color=st.INK,
        lw=st.LW_SUMMARY,
        zorder=3,
        solid_capstyle="round",
    )
    pos = repel(np.array(ends), 0.03)
    for p, c in zip(pos, names, strict=True):
        ax.text(5.15, p, c, fontsize=FSS, va="center", color=st.INK2)
    m1 = z[z.m == 1].R_full
    ax.plot([1], [ylo + 0.012], "v", color=st.ZEROSHOT, ms=3.2, clip_on=False)
    ax.text(
        1.12,
        ylo + 0.01,
        f"m = 1: {m1.max():.2f} to {m1.min():.2f}\n(off-scale; unshrunk)".replace("-", "−"),
        fontsize=FSS,
        color=st.INK2,
        va="bottom",
    )
    st.ref_line(ax, 0)
    ax.set_xlim(0.75, 6.05)
    ax.set_ylim(ylo, yhi)
    ax.set_xticks(M)
    ax.set_xlabel("Number of source contexts, m")
    ax.set_ylabel("Zero-shot full response explained\n(template-removed, M1)")
    ax.set_title("Zero-shot starting point (k = 0)", loc="left")
    ax.annotate(
        "median of 6\n(descriptive)",
        xy=(2.5, float(med.loc[2] + med.loc[3]) / 2),
        xytext=(1.25, 0.3),
        fontsize=FSS,
        color=st.INK,
        arrowprops={"arrowstyle": "-", "lw": 0.5, "color": st.INK},
    )


def panel_change(ax, t):
    """Per-target zero-shot M1 at m = 2 and m = 5 (the same values as panel B), with the printed difference."""
    z = t[t.k == 0]
    ctx = st.CONTEXTS[::-1]
    for i, c in enumerate(ctx):
        if i % 2 == 0:
            ax.axhspan(i - 0.5, i + 0.5, color="#F7F7F7", zorder=0, lw=0)
        x = z[z.context == c].set_index("m")
        a, b = x.loc[2], x.loc[5]
        ax.annotate(
            "",
            xy=(b.R_full, i),
            xytext=(a.R_full, i),
            arrowprops={
                "arrowstyle": "-|>",
                "color": st.ZEROSHOT,
                "lw": 1.0,
                "mutation_scale": 7,
                "shrinkA": 2.5,
                "shrinkB": 2.5,
            },
            zorder=2,
        )
        st.herrorbar(
            ax,
            [a.R_full],
            [i],
            [a.R_full_k0_rep_min],
            [a.R_full_k0_rep_max],
            st.ZEROSHOT,
            marker="o",
            mfc="white",
            ms=3.2,
            lw=0.7,
        )
        st.herrorbar(
            ax,
            [b.R_full],
            [i],
            [b.R_full_k0_rep_min],
            [b.R_full_k0_rep_max],
            st.EMPH,
            marker="o",
            mfc="white",
            ms=3.2,
            lw=0.7,
        )
        d = b.R_full - a.R_full
        ax.text(
            1.04,
            i,
            f"{d:+.2f}".replace("-", "−"),
            transform=ax.get_yaxis_transform(),
            fontsize=FSS,
            va="center",
            ha="left",
            fontweight="bold" if d < 0 else "normal",
        )
    ax.text(
        1.04,
        len(ctx) - 0.45,
        "Δ",
        transform=ax.get_yaxis_transform(),
        fontsize=FSS,
        va="bottom",
        ha="left",
        color=st.INK2,
    )
    st.ref_line(ax, 0, orient="v")
    ax.set_yticks(range(len(ctx)))
    ax.set_yticklabels(ctx)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.5, len(ctx) - 0.5)
    ax.set_xlim(-0.15, 0.33)
    ax.set_xticks([-0.1, 0, 0.1, 0.2, 0.3])
    ax.set_xlabel("Zero-shot full response explained (M1)")
    ax.set_title("Change from m = 2 to m = 5", loc="left")
    ax.legend(
        handles=[
            Line2D(
                [], [], color=st.ZEROSHOT, marker="o", mfc="white", ls="", ms=3.2, label="m = 2"
            ),
            Line2D([], [], color=st.EMPH, marker="o", mfc="white", ls="", ms=3.2, label="m = 5"),
        ],
        loc="upper left",
        ncol=1,
        columnspacing=0.8,
        handletextpad=0.2,
        frameon=True,
        facecolor="white",
        edgecolor="none",
        framealpha=1,
    )


def main() -> None:
    t = load("fig3_source_count")
    fig = st.new_figure(st.DOUBLE, 66 / 25.4)
    ax_a = panel_design(fig, [0.02, 0.03, 0.25, 0.94])
    ax_b = fig.add_axes([0.355, 0.17, 0.26, 0.72])
    panel_zero_shot(ax_b, t)
    ax_c = fig.add_axes([0.74, 0.17, 0.19, 0.72])
    panel_change(ax_c, t)
    st.panel_label(ax_a, "A", x=-0.002, y=0.94)
    st.panel_label(ax_b, "B", x=-0.075, y=1.04)
    st.panel_label(ax_c, "C", x=-0.075, y=1.04)
    save(
        fig,
        "fig3_source_count",
        sources=["fig3_source_count"],
        plot_script="scripts/paper_figures/plot_fig3.py",
        upstream=["outputs/n3/n3b (N3-B, reports/n3_results.md §3; zero-shot M1 only)"],
    )


if __name__ == "__main__":
    main()
