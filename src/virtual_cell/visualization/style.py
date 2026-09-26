"""One restrained, consistent matplotlib style for every figure in the repo.

Colours are the Okabe-Ito colour-blind-safe set. Colour is never the only
thing that separates groups: every categorical encoding here also has a hatch,
a marker or a direct label.
"""

from __future__ import annotations

import matplotlib as mpl

# Okabe-Ito
BLUE = "#0072B2"
ORANGE = "#E69F00"
GREEN = "#009E73"
VERMILION = "#D55E00"
SKY = "#56B4E9"
PURPLE = "#CC79A7"
YELLOW = "#F0E442"
BLACK = "#000000"
GREY = "#8C8C8C"
LIGHT_GREY = "#D9D9D9"

#: The four decomposition components, with the project's fixed wording.
COMPONENTS = {
    "template": {"label": "template", "color": GREY, "hatch": ""},
    "beta": {"label": "conserved effect β", "color": BLUE, "hatch": "//"},
    "gamma": {"label": "context-specific interaction γ", "color": ORANGE, "hatch": "\\\\"},
    "noise": {"label": "measurement noise", "color": "white", "hatch": ".."},
}

CONTEXTS = {
    "K562": {"color": BLUE, "marker": "o"},
    "RPE1": {"color": GREEN, "marker": "s"},
    "HepG2": {"color": VERMILION, "marker": "D"},
    "Jurkat": {"color": PURPLE, "marker": "^"},
}
CONTEXT_ORDER = ["K562", "RPE1", "HepG2", "Jurkat"]

TIERS = {
    2: {"label": "Tier 2", "color": BLUE, "hatch": "//"},
    1: {"label": "Tier 1", "color": SKY, "hatch": ".."},
    0: {"label": "Tier 0", "color": LIGHT_GREY, "hatch": ""},
}

#: Arms of the unseen-perturbation comparison.
DIRECT = {"label": "direct transfer", "color": BLUE, "marker": "o"}
PRIOR = {"label": "unseen-perturbation prior", "color": VERMILION, "marker": "s"}

GENERATORS = {
    "G0_control_resample": {"label": "G0 control resampling", "color": GREY, "hatch": ""},
    "G1_transport": {"label": "G1 control transport", "color": BLUE, "hatch": "//"},
    "G2_count_model": {"label": "G2 count model", "color": ORANGE, "hatch": "\\\\"},
}

DATASETS = {
    "arch1": {"label": "arch1", "color": BLUE, "marker": "o"},
    "kaden25rpe1": {"label": "Kaden RPE1", "color": VERMILION, "marker": "s"},
    "replogle22rpe1": {"label": "Replogle RPE1", "color": GREEN, "marker": "D"},
    "replogle22k562": {"label": "Replogle K562", "color": SKY, "marker": "^"},
    "nadig25hepg2": {"label": "Nadig HepG2", "color": PURPLE, "marker": "v"},
    "nadig25jurkat": {"label": "Nadig Jurkat", "color": ORANGE, "marker": "P"},
}

FULL_WIDTH = 7.2  # inches; readable when GitHub scales to README width
HALF_WIDTH = 3.5


def apply() -> None:
    """Install the project style (idempotent)."""
    mpl.rcParams.update(
        {
            "figure.dpi": 110,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.05,
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.titleweight": "bold",
            "axes.titlelocation": "left",
            "axes.labelsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "axes.grid": False,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "legend.frameon": False,
            "hatch.linewidth": 0.6,
            "lines.linewidth": 1.2,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )


def reference_line(ax, value: float = 0.0, *, axis: str = "y", **kw) -> None:
    """A thin dashed reference line (zero, one, analytic floor, ...)."""
    style = {"color": BLACK, "lw": 0.7, "ls": "--", "zorder": 0}
    style.update(kw)
    (ax.axhline if axis == "y" else ax.axvline)(value, **style)


def panel_label(ax, letter: str) -> None:
    ax.text(
        -0.02,
        1.04,
        letter,
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        ha="right",
        va="bottom",
    )
