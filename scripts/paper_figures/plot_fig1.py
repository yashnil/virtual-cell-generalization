# ruff: noqa: E501 (figure labels and provenance notes are kept as single strings)
"""Figure 1: conserved and context-specific perturbation effects differ in reproducibility (background)."""

from __future__ import annotations

import numpy as np
import style as st
from _common import load, save
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Rectangle

FS, FSS = st.FS, st.FS_SMALL


def panel_setup(fig, rect):
    """A: who knows what. All geometry in millimetres."""
    ax, W, H = st.schematic_axes(fig, rect)
    gw, gh = 22.0, 19.0  # target matrix glyph
    top = H - 3.0
    # --- column x positions (equal gaps)
    xs, xz, xc = 2.0, 42.0, 74.0
    # sources: three stacked cards
    for i in range(3):
        st.cell_grid(
            ax, xs + 2.2 * i, 12.5 - 2.2 * i + 4.4, 20.0, 17.0, 6, 6, st.ZEROSHOT, seed=11 + i
        )
    ax.text(
        xs + 12.2,
        top - 4.6,
        "Source contexts",
        fontsize=FS,
        fontweight="bold",
        ha="center",
        va="top",
    )
    ax.text(
        xs + 12.2,
        top - 8.2,
        "responses measured",
        fontsize=FSS,
        color=st.INK2,
        ha="center",
        va="top",
    )
    ax.text(
        xs + 12.2, 10.6, "perturbations × genes", fontsize=FSS, color=st.INK2, ha="center", va="top"
    )
    # header spanning the two target regimes
    ax.text(
        (xz + xc + gw) / 2,
        top,
        "Held-out target context",
        fontsize=FS,
        fontweight="bold",
        ha="center",
        va="top",
    )
    ax.plot([xz, xc + gw], [top - 3.4, top - 3.4], color=st.INK2, lw=st.STROKE)
    arrow_y = 12.5 + gh / 2
    st.arrow(ax, (xs + 27.5, arrow_y), (xz - 3.0, arrow_y))
    for x0, title, rows, note in (
        (xz, "Zero shot", [], "source responses\n+ target controls"),
        (xc, "Target calibration", [1, 3], "+ k measured target\nperturbations (anchors)"),
    ):
        ax.text(
            x0 + gw / 2, top - 4.6, title, fontsize=FS, ha="center", va="top", fontweight="bold"
        )
        st.box(ax, x0, 32.8, gw, 3.0, fc=st.tint(st.ZEROSHOT, 0.7), ec=st.INK2, radius=0.0)
        ax.text(
            x0 + gw / 2,
            34.3,
            "target controls",
            fontsize=FSS,
            color=st.INK,
            ha="center",
            va="center",
        )
        st.cell_grid(ax, x0, 12.5, gw, gh, 7, 6, st.EMPH, rows_on=rows, seed=3)
        ax.text(
            x0 + gw / 2,
            10.6,
            note,
            fontsize=FSS,
            color=st.INK2,
            ha="center",
            va="top",
            linespacing=1.15,
        )
    # anchor bracket
    ch = gh / 7
    yt, yb = 12.5 + gh - 1 * ch, 12.5 + gh - 4 * ch
    xb = xc + gw + 1.0
    ax.plot([xb, xb + 0.8, xb + 0.8, xb], [yt, yt, yb, yb], color=st.INK, lw=st.STROKE)
    ax.text(xb + 1.4, (yt + yb) / 2, "k", fontsize=FS, va="center", style="italic")
    # information-boundary key
    ky = 1.2
    items = [
        (st.tint(st.ZEROSHOT, 0.7), None, "measured, available to predictor"),
        (st.tint(st.EMPH, 0.8), None, "measured target anchor"),
        ("white", "//////", "unknown, held out for evaluation"),
    ]
    x = 2.0
    for fc, hatch, lab in items:
        ax.add_patch(Rectangle((x, ky), 2.6, 2.0, fc=fc, ec=st.INK2, lw=0.4, hatch=hatch))
        t = ax.text(x + 3.4, ky + 1.0, lab, fontsize=FSS, va="center", color=st.INK)
        fig.canvas.draw()
        bb = t.get_window_extent(renderer=fig.canvas.get_renderer())
        x = ax.transData.inverted().transform((bb.x1, 0))[0] + 3.2
    return ax


def panel_equation(fig, rect):
    ax, W, H = st.schematic_axes(fig, rect)
    r = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    y = H - 13.0
    terms = [
        (r"$\delta(c,p)$", st.INK, None),
        ("=", st.INK, None),
        (r"$\mu$", st.TEMPLATE, "global\nresponse"),
        ("+", st.INK, None),
        (r"$\alpha(c)$", st.TEMPLATE, "context-wide\ncomponent"),
        ("+", st.INK, None),
        (r"$\beta(p)$", st.BETA, "conserved\nperturbation\neffect"),
        ("+", st.INK, None),
        (r"$\gamma(c,p)$", st.GAMMA, "context-specific\nperturbation\ninteraction"),
    ]
    # measure total width first, then centre
    widths = []
    for s_, _, _ in terms:
        t = ax.text(0, y, s_, fontsize=12.5)
        bb = t.get_window_extent(renderer=r)
        widths.append(inv.transform((bb.x1, 0))[0] - inv.transform((bb.x0, 0))[0])
        t.remove()
    gap = 2.2
    x = (W - (sum(widths) + gap * (len(terms) - 1))) / 2
    centres = {}
    for (s_, col, note), w in zip(terms, widths, strict=True):
        ax.text(x, y, s_, fontsize=12.5, color=col, va="center", ha="left")
        if note:
            ax.text(
                x + w / 2,
                y - 4.2,
                note,
                fontsize=FSS,
                color=col,
                va="top",
                ha="center",
                linespacing=1.1,
            )
        centres[s_] = (x, x + w)
        x += w + gap
    ax.text(
        W / 2,
        H - 3.0,
        "Response of perturbation p in context c",
        fontsize=FS,
        ha="center",
        va="top",
        fontweight="bold",
    )
    # bracket grouping the template terms
    (a0, _), (_, a1) = centres[r"$\mu$"], centres[r"$\alpha(c)$"]
    yb = y + 4.2
    ax.plot([a0, a0, a1, a1], [yb - 0.8, yb, yb, yb - 0.8], color=st.TEMPLATE, lw=st.STROKE)
    ax.text(
        (a0 + a1) / 2,
        yb + 0.6,
        "template",
        fontsize=FSS,
        color=st.TEMPLATE,
        ha="center",
        va="bottom",
    )

    # illustration: one perturbation in four contexts (illustrative, not data)
    rng = np.random.default_rng(7)
    ctxs = ["K562", "RPE1", "HepG2", "Jurkat"]
    sw, sh, n = 11.0, 2.6, 8
    gx = 14.0
    x0 = (W - (4 * sw + 3 * 2.4)) / 2 + 4.0
    yb_ = 9.0
    beta = rng.uniform(0.2, 1.0, size=n)
    for i, c in enumerate(ctxs):
        xx = x0 + i * (sw + 2.4)
        ax.text(xx + sw / 2, yb_ + 2 * sh + 2.6, c, fontsize=FSS, ha="center", color=st.INK2)
        g = rng.uniform(0.15, 1.0, size=n)
        for j in range(n):
            ax.add_patch(
                Rectangle(
                    (xx + j * sw / n, yb_ + sh + 0.9),
                    sw / n,
                    sh,
                    fc=st.tint(st.BETA, beta[j]),
                    ec="white",
                    lw=0.3,
                )
            )
            ax.add_patch(
                Rectangle(
                    (xx + j * sw / n, yb_),
                    sw / n,
                    sh,
                    fc=st.tint(st.GAMMA, g[j]),
                    ec="white",
                    lw=0.3,
                )
            )
    ax.text(
        x0 - 1.5,
        yb_ + 1.5 * sh + 0.9,
        "β: same",
        fontsize=FSS,
        ha="right",
        va="center",
        color=st.BETA,
    )
    ax.text(
        x0 - 1.5, yb_ + sh / 2, "γ: differs", fontsize=FSS, ha="right", va="center", color=st.GAMMA
    )
    del gx
    ax.text(
        W / 2,
        0.6,
        "Decomposition framework following Molina & Zhang (2026);\nindependently re-derived here as background, not a new method.",
        fontsize=FSS,
        color=st.INK2,
        ha="center",
        va="bottom",
        style="italic",
    )
    return ax


COLOURS = {"template": st.TEMPLATE, "beta": st.BETA, "gamma": st.GAMMA, "noise": "#9A9A9A"}
MARKERS = {"template": "s", "beta": "o", "gamma": "D", "noise": "o"}


def panel_energy(ax, t):
    n = len(t)
    for i, r in enumerate(t.itertuples()):
        y = n - 1 - i
        col = COLOURS[r.component]
        ax.plot(
            [r.energy_pct_variant_min, r.energy_pct_variant_max],
            [y, y],
            color=st.tint(col, 0.45),
            lw=2.6,
            solid_capstyle="round",
            zorder=1,
        )
        ax.plot(
            r.energy_pct_corrected,
            y,
            MARKERS[r.component],
            color=col,
            ms=st.MS + 1.2,
            mec="white",
            mew=0.5,
            zorder=3,
        )
        ax.text(
            r.energy_pct_corrected,
            y + 0.3,
            f"{r.energy_pct_corrected:.1f}",
            fontsize=st.FS_SMALL,
            ha="center",
            va="bottom",
            color=st.INK,
        )
        if r.component in ("beta", "gamma"):
            ax.plot(
                r.energy_pct_uncorrected,
                y,
                MARKERS[r.component],
                color=col,
                mfc="white",
                ms=st.MS + 0.6,
                zorder=2,
            )
            ax.add_patch(
                FancyArrowPatch(
                    (r.energy_pct_uncorrected - 0.8, y),
                    (r.energy_pct_corrected + 1.0, y),
                    arrowstyle="-|>",
                    mutation_scale=5,
                    lw=0.6,
                    color=st.INK2,
                    zorder=2,
                )
            )
    g = t[t.component == "gamma"].iloc[0]
    ax.text(
        g.energy_pct_uncorrected,
        1 - 0.32,
        f"uncorrected\n{g.energy_pct_uncorrected:.1f}",
        fontsize=st.FS_SMALL,
        ha="center",
        va="top",
        color=st.INK2,
        linespacing=1.0,
    )
    ax.set_yticks(range(n))
    ax.set_yticklabels(t.label[::-1])
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 46)
    ax.set_ylim(-0.6, n - 0.3)
    ax.set_xlabel("Share of total response energy (%)")
    ax.set_title("Response energy (noise-corrected)", loc="left")


def panel_repro(ax, t):
    n = len(t)
    for i, r in enumerate(t.itertuples()):
        y = n - 1 - i
        col = COLOURS[r.component]
        if r.component == "noise":
            ax.text(
                50,
                y,
                "— (not a signal component)",
                fontsize=st.FS_SMALL,
                color=st.INK2,
                ha="center",
                va="center",
            )
            continue
        if r.component in ("beta", "gamma"):
            ax.plot(
                [r.reproducible_pct_variant_min, r.reproducible_pct_variant_max],
                [y, y],
                color=st.tint(col, 0.45),
                lw=2.6,
                solid_capstyle="round",
                zorder=1,
            )
        ax.plot(
            r.reproducible_pct,
            y,
            MARKERS[r.component],
            color=col,
            ms=st.MS + 1.2,
            mec="white",
            mew=0.5,
            zorder=3,
        )
        ax.text(
            r.reproducible_pct,
            y + 0.3,
            f"{r.reproducible_pct:.1f}",
            fontsize=st.FS_SMALL,
            ha="center",
            va="bottom",
        )
    st.ref_line(ax, 100, orient="v", ls=(0, (2, 2)))
    ax.set_yticks(range(n))
    ax.set_yticklabels([])
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 105)
    ax.set_ylim(-0.6, n - 0.3)
    ax.set_xlabel("Reproducible across split halves (%)")
    ax.set_title("Reproducibility", loc="left")


def main() -> None:
    t = load("fig1_components")
    fig = st.new_figure(st.DOUBLE, 108 / 25.4)
    ax_a = panel_setup(fig, [0.01, 0.52, 0.575, 0.44])
    ax_b = panel_equation(fig, [0.6, 0.52, 0.39, 0.44])
    ax_c = fig.add_axes([0.17, 0.115, 0.37, 0.3])
    ax_d = fig.add_axes([0.6, 0.115, 0.37, 0.3])
    panel_energy(ax_c, t)
    panel_repro(ax_d, t)
    h = [
        Line2D(
            [], [], marker="o", color="#777777", ls="", ms=st.MS + 1, label="canonical estimate"
        ),
        Line2D(
            [],
            [],
            marker="o",
            color="#777777",
            mfc="white",
            ls="",
            ms=st.MS + 0.6,
            label="before noise correction",
        ),
        Line2D([], [], color="#CCCCCC", lw=2.6, label="range over 21 preprocessing variants"),
    ]
    fig.legend(handles=h, loc="lower center", bbox_to_anchor=(0.57, -0.005), ncol=3, frameon=False)
    st.panel_label(ax_a, "A", x=0.01, y=0.975)
    st.panel_label(ax_b, "B", x=0.0, y=0.975)
    st.panel_label(ax_c, "C", x=-0.13)
    st.panel_label(ax_d, "D", x=-0.02)
    save(
        fig,
        "fig1_transfer_components",
        sources=["fig1_components"],
        plot_script="scripts/paper_figures/plot_fig1.py",
        upstream=[
            "outputs/four_context_v1/summary.json",
            "outputs/four_context_sensitivity/variant_table.csv",
        ],
    )


if __name__ == "__main__":
    main()
