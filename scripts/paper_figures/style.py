# ruff: noqa: E501 (figure labels and provenance notes are kept as single strings)
"""One visual system for every paper figure (main + Extended Data).

Semantic colour mapping (fixed; never re-used for another meaning):

=====================================  =========  ===============================================
meaning                                colour     secondary channel
=====================================  =========  ===============================================
template / scale (global adaptation)   TEMPLATE   square marker; dashed line when beside γ
conserved β                            BETA       circle marker; "β" label
context-specific γ / γ⊥                GAMMA      diamond marker; "γ" label
measurement noise                      NOISE      hatch; label
target-calibrated / emphasised source  EMPH       filled marker; bold label
zero-shot / source-only / other source ZEROSHOT   open marker; label
failed gate / non-interpretable        FAILED     italic label
positive control                       CONTROL    open square
uncertainty                            line colour at ALPHA_CI fill, or thin error bar
sequential magnitude (γ heatmap)       GAMMA_CMAP (white → vermilion → brown)
=====================================  =========  ===============================================

The categorical triple TEMPLATE / BETA / GAMMA passes the dataviz ``validate_palette.js`` checks on a white
surface for all pairs (lightness band, chroma floor, CVD ΔE ≥ 9, normal-vision ΔE ≥ 18.8, contrast ≥ 3:1).
Contexts are never given six hues: they are identified by facet titles or direct labels.
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# ------------------------------------------------------------------ colours
TEMPLATE = "#2A9D78"
BETA = "#1F6FB4"
GAMMA = "#D55E00"
NOISE = "#D0D0D0"
EMPH = "#1A1A1A"
ZEROSHOT = "#8F8F8F"
FAILED = "#BDBDBD"
CONTROL = "#3A3A3A"
INK = "#1A1A1A"  # primary text
INK2 = "#555555"  # secondary text
REF = "#B5B5B5"  # reference lines (0, gates)
ALPHA_CI = 0.22


def tint(hex_colour: str, amount: float = 0.35) -> str:
    """Blend ``hex_colour`` toward white; ``amount`` = share of the original colour kept."""
    c = np.array(mpl.colors.to_rgb(hex_colour))
    return mpl.colors.to_hex(1 - amount * (1 - c))


GAMMA_CMAP = LinearSegmentedColormap.from_list(
    "gamma_seq", ["#FFFFFF", "#FBD9C2", "#F2A270", GAMMA, "#8A3300", "#4A1B00"]
)
GAMMA_CMAP.set_under("#E3E3E3")

# Six contexts in a fixed reading order (scPertEval first, X-Atlas last).
CONTEXTS = ["K562", "RPE1", "HepG2", "Jurkat", "HCT116", "HEK293T"]
ORIGINAL4 = CONTEXTS[:4]
LAB = {
    "K562": "Replogle 2022",
    "RPE1": "Replogle 2022",
    "HepG2": "Nadig 2025",
    "Jurkat": "Nadig 2025",
    "HCT116": "X-Atlas/Orion",
    "HEK293T": "X-Atlas/Orion",
}

# ------------------------------------------------------------------ geometry (inches)
MM = 1 / 25.4
DOUBLE = 183 * MM
SINGLE = 89 * MM

FS = 7.0  # body
FS_SMALL = 6.0  # minimum used anywhere (journal floor)
FS_LABEL = 8.5  # panel letters
LW = 1.0  # data lines
LW_THIN = 0.6
LW_SUMMARY = 1.9
MS = 3.4  # marker size (points)


def apply() -> None:
    mpl.rcParams.update(
        {
            "font.family": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": FS,
            "axes.titlesize": FS,
            "axes.titleweight": "normal",
            "axes.labelsize": FS,
            "axes.labelcolor": INK,
            "axes.edgecolor": INK,
            "axes.linewidth": LW_THIN,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": False,
            "axes.titlepad": 3.0,
            "axes.labelpad": 2.0,
            "xtick.labelsize": FS_SMALL,
            "ytick.labelsize": FS_SMALL,
            "xtick.color": INK,
            "ytick.color": INK,
            "xtick.major.width": LW_THIN,
            "ytick.major.width": LW_THIN,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "xtick.minor.size": 1.5,
            "ytick.minor.size": 1.5,
            "xtick.major.pad": 1.5,
            "ytick.major.pad": 1.5,
            "legend.fontsize": FS_SMALL,
            "legend.frameon": False,
            "legend.handlelength": 1.6,
            "lines.linewidth": LW,
            "lines.markersize": MS,
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "mathtext.fallback": "stixsans",
            "text.color": INK,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.dpi": 600,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "hatch.linewidth": 0.4,
            "hatch.color": "#CFCFCF",
        }
    )


def panel_label(ax, letter: str, x: float = -0.02, y: float = 1.0, fig=None) -> None:
    """Bold upper-case panel letter at the top-left of ``ax`` (in figure coordinates)."""
    fig = fig or ax.figure
    bbox = ax.get_position()
    fig.text(
        bbox.x0 + x,
        bbox.y1 + 0.01 + (y - 1.0),
        letter,
        fontsize=FS_LABEL,
        fontweight="bold",
        va="bottom",
        ha="right",
    )


def ref_line(ax, value: float = 0.0, orient: str = "h", **kw) -> None:
    style = {"color": REF, "lw": LW_THIN, "zorder": 0, "ls": "-"}
    style.update(kw)
    (ax.axhline if orient == "h" else ax.axvline)(value, **style)


# ------------------------------------------------------------------ the target-budget axis
K0_POS = -0.55  # where k = 0 sits on the log10 axis (left of k = 1 at 0)


def kpos(k):
    """Map budgets to x positions: log10(k) for k ≥ 1, K0_POS for k = 0."""
    k = np.asarray(k, dtype=float)
    with np.errstate(divide="ignore"):
        return np.where(k > 0, np.log10(np.maximum(k, 1e-9)), K0_POS)


def budget_axis(
    ax, budgets, kref: int, *, label: bool = True, marks=(20, 50, 100), show_zero=True
) -> None:
    """Log budget axis with every tested budget ticked and k = 0 set off by a break."""
    ticks = [b for b in budgets if show_zero or b > 0]
    ax.set_xticks(kpos(ticks))
    ax.set_xticklabels([str(int(b)) for b in ticks])
    ax.minorticks_off()
    lo = K0_POS - 0.22 if show_zero else -0.15
    ax.set_xlim(lo, np.log10(kref) + 0.12)
    for m in marks:
        ax.axvline(kpos(m), color="#E6E6E6", lw=0.8, zorder=0)
    ax.axvline(kpos(kref), color="#E6E6E6", lw=0.8, zorder=0)
    if show_zero:
        # axis break between k = 0 and k = 1
        xb = (K0_POS + 0.0) / 2 - 0.02
        trans = ax.get_xaxis_transform()
        for dx in (-0.035, 0.035):
            ax.plot(
                [xb + dx - 0.03, xb + dx + 0.03],
                [-0.025, 0.025],
                transform=trans,
                color=INK,
                lw=LW_THIN,
                clip_on=False,
                zorder=5,
            )
        ax.plot(
            [xb - 0.035, xb + 0.035],
            [0, 0],
            transform=trans,
            color="white",
            lw=1.6,
            clip_on=False,
            zorder=4,
        )
    if label:
        ax.set_xlabel("Measured target perturbations, k (log scale)")


def direct_label(ax, x, y, text, color=INK, **kw) -> None:
    style = {"fontsize": FS_SMALL, "va": "center", "ha": "left", "color": color}
    style.update(kw)
    ax.text(x, y, text, **style)


def errorbar(ax, x, y, lo, hi, color, *, marker="o", mfc=None, ms=MS, lw=0.8, zorder=3, **kw):
    y = np.asarray(y, dtype=float)
    yerr = np.vstack([y - np.asarray(lo, dtype=float), np.asarray(hi, dtype=float) - y])
    return ax.errorbar(
        x,
        y,
        yerr=yerr,
        fmt=marker,
        color=color,
        mfc=mfc if mfc is not None else color,
        mec=color,
        ms=ms,
        elinewidth=lw,
        capsize=0,
        lw=0,
        zorder=zorder,
        **kw,
    )


def herrorbar(ax, x, y, lo, hi, color, *, marker="o", mfc=None, ms=MS, lw=0.9, zorder=3, **kw):
    x = np.asarray(x, dtype=float)
    xerr = np.vstack([x - np.asarray(lo, dtype=float), np.asarray(hi, dtype=float) - x])
    return ax.errorbar(
        x,
        y,
        xerr=xerr,
        fmt=marker,
        color=color,
        mfc=mfc if mfc is not None else color,
        mec=color,
        ms=ms,
        elinewidth=lw,
        capsize=0,
        lw=0,
        zorder=zorder,
        **kw,
    )


def new_figure(width: float, height: float):
    apply()
    return plt.figure(figsize=(width, height))


def budget_line(ax, k, y, color, *, lw=LW, marker=None, ms=MS, mfc=None, zorder=2, **kw):
    """Draw a learning curve without connecting k = 0 (zero-shot) across the axis break."""
    k = np.asarray(k, dtype=float)
    y = np.asarray(y, dtype=float)
    pos = kpos(k)
    on = k > 0
    ax.plot(pos[on], y[on], "-", color=color, lw=lw, zorder=zorder, **kw)
    if marker:
        ax.plot(
            pos,
            y,
            marker,
            ls="",
            color=color,
            ms=ms,
            mfc=mfc or color,
            mec="white",
            mew=0.5,
            zorder=zorder + 1,
            clip_on=True,
        )


# ------------------------------------------------------------------ schematic geometry
# All schematics share one corner radius, one stroke weight and one arrow geometry. Schematic axes are given
# equal data/physical aspect (see ``schematic_axes``) so that rounded corners are circular everywhere.
RADIUS = 0.6  # schematic units (1 unit = 1 mm at print size)
STROKE = 0.6
ARROW_MS = 8


def schematic_axes(fig, rect):
    """Axes whose data units are millimetres at print size (equal aspect), with no frame."""
    ax = fig.add_axes(rect)
    w_mm = rect[2] * fig.get_figwidth() * 25.4
    h_mm = rect[3] * fig.get_figheight() * 25.4
    ax.set_xlim(0, w_mm)
    ax.set_ylim(0, h_mm)
    ax.axis("off")
    return ax, w_mm, h_mm


def box(
    ax, x, y, w, h, *, fc="white", ec=INK2, lw=STROKE, ls="-", hatch=None, radius=RADIUS, zorder=1
):
    from matplotlib.patches import FancyBboxPatch

    p = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0,rounding_size={radius}",
        fc=fc,
        ec=ec,
        lw=lw,
        ls=ls,
        hatch=hatch,
        zorder=zorder,
    )
    ax.add_patch(p)
    return p


def arrow(ax, xy0, xy1, *, color=INK, lw=0.9, ls="-"):
    from matplotlib.patches import FancyArrowPatch

    ax.add_patch(
        FancyArrowPatch(
            xy0,
            xy1,
            arrowstyle="-|>",
            mutation_scale=ARROW_MS,
            lw=lw,
            color=color,
            ls=ls,
            shrinkA=0,
            shrinkB=0,
        )
    )


def cell_grid(ax, x, y, w, h, nrow, ncol, colour, *, rows_on=None, seed=0, unknown_hatch=True):
    """Perturbation × gene matrix glyph. Measured rows: tinted cells. Unknown rows: blank (optionally hatched)."""
    from matplotlib.patches import Rectangle

    rng = np.random.default_rng(seed)
    cw, ch = w / ncol, h / nrow
    rows_on = range(nrow) if rows_on is None else rows_on
    shade = rng.uniform(0.25, 1.0, size=(nrow, ncol))
    for i in range(nrow):
        on = i in rows_on
        for j in range(ncol):
            ax.add_patch(
                Rectangle(
                    (x + j * cw, y + h - (i + 1) * ch),
                    cw,
                    ch,
                    fc=tint(colour, shade[i, j]) if on else "white",
                    ec="white" if on else "#E4E4E4",
                    lw=0.35,
                    hatch="//////" if (unknown_hatch and not on) else None,
                    zorder=2,
                )
            )
    box(ax, x, y, w, h, fc="none", ec=INK2, radius=0.0, zorder=3)
