# ruff: noqa: E501 (figure labels and provenance notes are kept as single strings)
"""Figure 2: context-specific perturbation effects require substantially more target data to learn (N3-A)."""

from __future__ import annotations

import numpy as np
import style as st
from _common import load, save

KREF = 743
BUDGETS = [0, 1, 2, 5, 10, 20, 50, 100, 200, KREF]


def panel_template(ax, cur):
    ylo = -1.0
    ax.axhspan(ylo - 0.2, 0, color="#F3F3F3", zorder=0, lw=0)
    ax.text(
        st.kpos(1) - 0.1,
        -0.5,
        "worse than\nzero-shot",
        fontsize=st.FS_SMALL,
        color=st.INK2,
        ha="center",
        va="center",
    )
    for c in st.CONTEXTS:
        x = cur[cur.context == c].sort_values("k")
        st.budget_line(ax, x.k, x.template_frac, st.tint(st.TEMPLATE, 0.38), lw=st.LW_THIN)
        off = x[x.template_frac < ylo]
        ax.plot(
            st.kpos(off.k),
            np.full(len(off), ylo),
            "v",
            color=st.tint(st.TEMPLATE, 0.5),
            ms=2.6,
            clip_on=False,
            zorder=2,
        )
    med = cur.groupby("k").template_frac.median()
    st.budget_line(
        ax, med.index, med.to_numpy(), st.TEMPLATE, lw=st.LW_SUMMARY, marker="s", zorder=4
    )
    st.ref_line(ax, 0)
    st.ref_line(ax, 1, ls=(0, (2, 2)))
    ax.set_ylim(ylo - 0.06, 1.12)
    ax.set_yticks([-1, -0.5, 0, 0.5, 1])
    ax.set_ylabel("Template + scale recovered\n(fraction of gain at $k_{ref}$)")
    st.budget_axis(ax, BUDGETS, KREF)
    worst = cur[cur.k > 0].template_frac.min()
    ax.text(
        st.kpos(25),
        -0.78,
        f"▼ below −1, drawn at −1\n   (k ≤ 10; minimum {worst:.0f})".replace("-", "−") + "",
        fontsize=st.FS_SMALL,
        color=st.INK2,
        ha="left",
        va="center",
    )
    kt = cur.groupby("context").k_T50.first()
    ax.text(
        st.K0_POS - 0.15,
        1.0,
        f"Half of the $k_{{ref}}$ gain\nreached at k = {kt.min()}–{kt.max()};\n≥ 86 % by k = 200\nin every context",
        fontsize=st.FS_SMALL,
        color=st.INK,
        ha="left",
        va="top",
        linespacing=1.15,
    )
    ax.set_title(
        "Template + scale: saturates early\n$\\it{normalised\\ to\\ each\\ context's\\ own}$ $k_{ref}$ $\\it{gain}$",
        loc="left",
        color=st.INK,
        linespacing=1.25,
    )


def panel_gamma(ax, cur):
    for c in st.CONTEXTS:
        x = cur[cur.context == c].sort_values("k")
        st.budget_line(ax, x.k, x.gamma_M3_median, st.tint(st.GAMMA, 0.38), lw=st.LW_THIN)
    med = cur.groupby("k").gamma_M3_median.median()
    st.budget_line(
        ax,
        med.index,
        med.to_numpy(),
        st.GAMMA,
        lw=st.LW_SUMMARY,
        marker="D",
        ms=st.MS - 0.4,
        zorder=4,
    )
    st.ref_line(ax, 0)
    ax.set_ylim(-0.02, 0.5)
    ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
    ax.set_ylabel(
        "Reliable context-specific (γ$^{\\perp}$) energy\nrecovered on held-out perturbations"
    )
    st.budget_axis(ax, BUDGETS, KREF)
    ref = cur[cur.k == KREF].gamma_M3_median
    k20 = cur[cur.k == 20].gamma_M3_median
    ax.text(
        st.kpos(KREF) - 0.02,
        0.485,
        f"$k_{{ref}}$: {ref.min():.2f}–{ref.max():.2f}, not plateaued",
        fontsize=st.FS_SMALL,
        ha="right",
        va="top",
    )
    ax.annotate(
        f"k = 20: {k20.min():.2f}–{k20.max():.2f}",
        xy=(st.kpos(20), med.loc[20]),
        xytext=(st.kpos(1.3), 0.27),
        fontsize=st.FS_SMALL,
        color=st.INK,
        arrowprops={"arrowstyle": "-", "color": st.INK2, "lw": 0.5},
    )
    ax.set_title(
        "Context-specific γ$^{\\perp}$: still rising at $k_{ref}$\n$\\it{absolute\\ scale,\\ not\\ normalised}$",
        loc="left",
        color=st.INK,
        linespacing=1.25,
    )


def legend_ab(ax_a, ax_b):
    from matplotlib.lines import Line2D

    h = [
        Line2D([], [], color="#BBBBBB", lw=st.LW_THIN),
        Line2D([], [], color="#555555", lw=st.LW_SUMMARY),
    ]
    ax_b.legend(
        h,
        ["individual held-out contexts (n = 6; see C)", "median across contexts (descriptive)"],
        loc="upper left",
        bbox_to_anchor=(0.0, 0.9),
        handlelength=2.0,
        borderaxespad=0.2,
    )


def panel_small_multiples(fig, gs, cur):
    axes = []
    for i, c in enumerate(st.CONTEXTS):
        ax = fig.add_subplot(gs[i // 3, i % 3])
        axes.append(ax)
        x = cur[cur.context == c].sort_values("k")
        xp = st.kpos(x.k)
        on = (x.k > 0).to_numpy()
        ax.fill_between(
            xp[on],
            x.gamma_M3_boot_lo[on],
            x.gamma_M3_boot_hi[on],
            color=st.GAMMA,
            alpha=st.ALPHA_CI,
            lw=0,
            zorder=1,
        )
        st.budget_line(ax, x.k, x.gamma_M3_boot_point, st.GAMMA, marker="D", ms=2.4)
        st.ref_line(ax, 0)
        ax.set_ylim(-0.08, 0.48)
        ax.set_yticks([0, 0.2, 0.4])
        st.budget_axis(ax, [0, 1, 10, 100, KREF], KREF, label=False)
        ax.set_xticklabels(["0", "1", "10", "100", "743"])
        if i % 3:
            ax.set_yticklabels([])
        if i < 3:
            ax.set_xticklabels([])
        cls = x.class_N1_rule.iloc[0]
        ax.set_title(f"{c}", loc="left", fontsize=st.FS, fontweight="bold", pad=1.5)
        ax.text(
            1.0,
            1.0,
            st.LAB[c],
            transform=ax.transAxes,
            fontsize=st.FS_SMALL,
            color=st.INK2,
            ha="right",
            va="bottom",
        )
        k20 = float(x[x.k == 20].gamma_M3_median.iloc[0])
        kr = float(x[x.k == KREF].gamma_M3_median.iloc[0])
        tag = "PASS" if cls == "PASS" else "FAIL-A"
        ax.text(
            0.03,
            0.97,
            f"k20/$k_{{ref}}$ = {k20 / kr:.2f}\n{tag} (N1 rule)",
            transform=ax.transAxes,
            fontsize=st.FS_SMALL,
            va="top",
            ha="left",
            color=st.INK if tag == "PASS" else st.INK2,
        )
    return axes


def panel_gain_fraction(fig, gs, gf):
    axes = []
    ctx = st.CONTEXTS[::-1]
    ypos = {c: i + 1.3 for i, c in enumerate(ctx)}
    for j, k in enumerate([20, 50, 100]):
        ax = fig.add_subplot(gs[0, j])
        axes.append(ax)
        x = gf[gf.k == k]
        for comp, col, mk, mfc, dy in (
            ("template_scale", st.TEMPLATE, "s", "white", 0.17),
            ("gamma_perp", st.GAMMA, "D", st.GAMMA, -0.17),
        ):
            y = x[x.component == comp].set_index("context").loc[ctx]
            yy = np.array([ypos[c] for c in ctx]) + dy
            lo = np.clip(y.lo, -0.75, None)
            st.herrorbar(ax, y.point, yy, lo, y.hi, col, marker=mk, mfc=mfc, ms=2.8, lw=0.8)
            clipped = y.lo < -0.75
            ax.plot(
                np.full(clipped.sum(), -0.75),
                yy[clipped.to_numpy()],
                marker=4,
                ls="",
                color=col,
                ms=3,
                clip_on=False,
            )
            # descriptive median across the six contexts (no interval: none exists)
            med = x[x.component == comp].point.median()
            ax.plot(med, dy * 0.9, mk, color=col, mfc=mfc, ms=3.6, mew=1.0, zorder=4)
        ax.axhline(0.72, color=st.INK2, lw=0.5)
        if k == 20:
            ax.axvline(0.25, color=st.INK2, lw=0.6, ls=(0, (2, 2)), zorder=0)
        st.ref_line(ax, 0, orient="v")
        st.ref_line(ax, 1, orient="v", ls=(0, (1, 2)))
        ax.set_xlim(-0.78, 1.08)
        ax.set_xticks([0, 0.5, 1])
        ax.set_xticklabels(["0", "0.5", "1"])
        ax.set_ylim(-0.45, len(ctx) + 0.9)
        ax.set_yticks([0] + [ypos[c] for c in ctx])
        ax.set_yticklabels((["median (descr.)"] + ctx) if j == 0 else [])
        ax.tick_params(axis="y", length=0)
        ax.set_title(f"k = {k}", fontsize=st.FS, pad=1.5)
        for i, c in enumerate(ctx):
            if i % 2 == 0:
                ax.axhspan(ypos[c] - 0.5, ypos[c] + 0.5, color="#F7F7F7", zorder=0, lw=0)
    axes[1].set_xlabel("Fraction of the $k_{ref}$ value reached, G(k)")
    return axes


def main() -> None:
    cur = load("fig2_curves")
    gf = load("fig2_gain_fraction")
    fig = st.new_figure(st.DOUBLE, 5.55)
    top = fig.add_gridspec(1, 2, left=0.085, right=0.985, top=0.94, bottom=0.6, wspace=0.32)
    ax_a = fig.add_subplot(top[0, 0])
    ax_b = fig.add_subplot(top[0, 1])
    panel_template(ax_a, cur)
    panel_gamma(ax_b, cur)
    legend_ab(ax_a, ax_b)

    gs_c = fig.add_gridspec(
        2, 3, left=0.085, right=0.565, top=0.49, bottom=0.085, wspace=0.12, hspace=0.42
    )
    axes_c = panel_small_multiples(fig, gs_c, cur)
    fig.text(
        0.325, 0.025, "Measured target perturbations, k (log scale)", ha="center", fontsize=st.FS
    )
    fig.text(
        0.035,
        0.29,
        "Reliable γ$^{\\perp}$ energy recovered (95 % CI)",
        rotation=90,
        ha="center",
        va="center",
        fontsize=st.FS,
    )

    gs_d = fig.add_gridspec(1, 3, left=0.665, right=0.985, top=0.49, bottom=0.115, wspace=0.12)
    axes_d = panel_gain_fraction(fig, gs_d, gf)
    # D legend
    axes_d[2].plot([], [], "s", color=st.TEMPLATE, mfc="white", ms=3, label="template + scale")
    axes_d[2].plot([], [], "D", color=st.GAMMA, ms=3, label="γ$^{\\perp}$")
    axes_d[2].plot(
        [], [], color=st.INK2, lw=0.6, ls=(0, (2, 2)), label="preregistered C3 bar (k = 20)"
    )
    axes_d[0].legend(
        *axes_d[2].get_legend_handles_labels(),
        loc="upper left",
        bbox_to_anchor=(-0.45, -0.17),
        ncol=3,
        columnspacing=1.0,
        handletextpad=0.3,
    )

    st.panel_label(ax_a, "A", x=-0.065)
    st.panel_label(ax_b, "B", x=-0.065)
    st.panel_label(axes_c[0], "C", x=-0.05)
    st.panel_label(axes_d[0], "D", x=-0.065)
    save(
        fig,
        "fig2_sample_complexity",
        sources=["fig2_curves", "fig2_gain_fraction"],
        plot_script="scripts/paper_figures/plot_fig2.py",
        upstream=["outputs/n3/n3a (N3-A, reports/n3_results.md §1)"],
    )


if __name__ == "__main__":
    main()
