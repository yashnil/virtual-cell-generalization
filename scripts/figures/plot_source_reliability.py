"""Figure 7 — Kaden source-reliability diagnostic (descriptive, not a model)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from virtual_cell.analysis.loco import spearman_brown
from virtual_cell.visualization import common, style

ME_ROWS = [
    ("arch1", "G_ARC", "all", "arch1 · all"),
    ("arch1", "G_ARC", "arc_targets", "arch1 · Arc targets"),
    ("kaden25rpe1", "G_ARC", "all", "Kaden · all"),
    ("kaden25rpe1", "G_ARC", "arc_targets", "Kaden · Arc targets"),
    ("arch1", "G_ARC", "arc_matched", "arch1 · matched Arc"),
    ("kaden25rpe1", "G_ARC", "arc_matched", "Kaden · matched Arc"),
    ("kaden25rpe1", "G_RPE", "shared_kaden_rpe1", "Kaden · shared w/ Rep. RPE1"),
    ("replogle22rpe1", "G_RPE", "shared_kaden_rpe1", "Rep. RPE1 · shared w/ Kaden"),
    ("replogle22k562", "G_4CTX", "shared_four_context", "Rep. K562 · 4-context"),
    ("replogle22rpe1", "G_4CTX", "shared_four_context", "Rep. RPE1 · 4-context"),
    ("nadig25hepg2", "G_4CTX", "shared_four_context", "Nadig HepG2 · 4-context"),
    ("nadig25jurkat", "G_4CTX", "shared_four_context", "Nadig Jurkat · 4-context"),
]

SHORT_COMPARISON = {
    "arch1 vs Kaden, natural panels (frozen m_hat inputs)": "arch1–Kaden · all",
    "arch1 vs Kaden, Arc-target subsets": "arch1–Kaden · Arc subsets",
    "arch1 vs Kaden, all shared perturbations": "arch1–Kaden · shared (39)",
    "arch1 vs Kaden, matched Arc targets": "arch1–Kaden · matched Arc (7)",
    "Kaden vs Replogle RPE1, natural panels": "Kaden–Rep. RPE1 · all",
    "Kaden vs Replogle RPE1, shared perturbations": "Kaden–Rep. RPE1 · shared (205)",
}


def main() -> None:
    style.apply()
    t = common.load_source("fig7_source_reliability")
    fig, axes = plt.subplots(
        2, 2, figsize=(style.FULL_WIDTH, 7.0), gridspec_kw={"wspace": 0.9, "hspace": 0.55}
    )

    # A — per-perturbation reliability on the identical axis
    ax = axes[0, 0]
    a = t[t.panel == "A"].set_index("dataset")
    nl = t[t.panel == "A_null"].set_index("dataset")
    order = ["arch1", "replogle22rpe1", "kaden25rpe1"]
    for i, ds in enumerate(order):
        meta = style.DATASETS[ds]
        r = a.loc[ds]
        ax.plot([r.q10, r.q90], [i, i], color=meta["color"], lw=1.2)
        ax.add_patch(
            plt.Rectangle(
                (r.q25, i - 0.18),
                r.q75 - r.q25,
                0.36,
                fc=meta["color"],
                ec="black",
                lw=0.6,
                alpha=0.85,
            )
        )
        ax.plot(r["median"], i, "|", color="black", ms=14, mew=2)
        n = nl.loc[ds]
        ax.plot(
            [n.q10, n.q90],
            [i + 0.34, i + 0.34],
            color=style.GREY,
            lw=3,
            alpha=0.7,
            solid_capstyle="butt",
        )
        ax.text(
            -0.38,
            i - 0.32,
            f"median {r['median']:.2f} · {int(r.value)} cells/pert",
            fontsize=6.3,
            va="center",
        )
    style.reference_line(ax, 0, axis="x")
    ax.set_yticks(range(len(order)), [style.DATASETS[d]["label"] for d in order])
    ax.set_xlim(-0.4, 1.0)
    ax.invert_yaxis()
    ax.set_xlabel(
        "per-perturbation split-half reliability\n"
        "(Spearman–Brown, 50 repeats; box q25–q75, line q10–q90)"
    )
    ax.set_title("Per-perturbation reliability\n(identical 8,166-gene axis)", fontsize=8.5)
    ax.plot(
        [],
        [],
        color=style.GREY,
        lw=3,
        alpha=0.7,
        label="control pseudo-perturbation null (q10–q90)",
    )
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.3), fontsize=6.5)
    style.panel_label(ax, "A")

    # B — main-effect reliability
    ax = axes[0, 1]
    b = t[t.panel == "B"]
    ys, labels = [], []
    for i, (ds, axis, panel, label) in enumerate(ME_ROWS):
        r = b[(b.dataset == ds) & (b.axis == axis) & (b.label == panel)]
        if not len(r):
            continue
        r = r.iloc[0]
        meta = style.DATASETS[ds]
        lo, hi = spearman_brown(r.q10), spearman_brown(r.q90)
        ax.errorbar(
            r["median"],
            i,
            xerr=[[r["median"] - lo], [hi - r["median"]]],
            fmt=meta["marker"],
            color=meta["color"],
            mec="black",
            ms=6,
            capsize=2,
            lw=1,
        )
        ys.append(i)
        labels.append(f"{label} (n={int(r.value)})")
    style.reference_line(ax, 0, axis="x")
    style.reference_line(ax, 1, axis="x", color=style.GREY)
    ax.set_yticks(ys, labels, fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlim(-0.2, 1.05)
    ax.set_xlabel(
        "main-effect split-half reliability (Spearman–Brown;\nbars = 5th–95th pct over repeats)"
    )
    ax.set_title("Main-effect reliability", fontsize=8.5)
    style.panel_label(ax, "B")

    # C — main-effect agreement between sources vs the noise ceiling
    ax = axes[1, 0]
    c = t[t.panel == "C"]
    c = c[~c.dataset.str.contains("four-context")].reset_index(drop=True)
    four = t[(t.panel == "C") & t.dataset.str.contains("four-context")]
    for i, r in c.iterrows():
        ax.plot([0, r.value], [i, i], color=style.LIGHT_GREY, lw=6, solid_capstyle="butt", zorder=0)
        ax.plot(
            r["median"],
            i,
            "o",
            color=style.VERMILION if "arch1" in r.dataset else style.GREEN,
            mec="black",
            ms=6,
        )
    if len(four):
        yi = len(c) + 0.5
        vals = four["median"].to_numpy()
        ax.plot([vals.min(), vals.max()], [yi, yi], color=style.BLUE, lw=2)
        ax.plot(vals, np.full(len(vals), yi), "D", color=style.BLUE, mec="black", ms=4)
        ceil = four["value"].to_numpy()
        ax.plot(
            [0, ceil.min()], [yi, yi], color=style.LIGHT_GREY, lw=6, solid_capstyle="butt", zorder=0
        )
    style.reference_line(ax, 0, axis="x")
    ylabels = [SHORT_COMPARISON.get(d, d) for d in c.dataset] + (
        ["4 research contexts (6 pairs)"] if len(four) else []
    )
    ax.set_yticks(
        list(range(len(c))) + ([len(c) + 0.5] if len(four) else []), ylabels, fontsize=6.5
    )
    ax.invert_yaxis()
    ax.set_xlim(-0.3, 1.05)
    ax.set_xlabel("Pearson between source main effects\n(grey bar = noise ceiling √(R_x R_y))")
    ax.set_title("Main-effect agreement\nbetween sources", fontsize=8.5)
    style.panel_label(ax, "C")

    # D — per-perturbation agreement on shared perturbations
    ax = axes[1, 1]
    d = t[t.panel == "D"].reset_index(drop=True)
    for i, r in d.iterrows():
        ax.plot([0, r.value], [i, i], color=style.LIGHT_GREY, lw=6, solid_capstyle="butt", zorder=0)
        color = style.VERMILION if "arch1" in r.dataset else style.GREEN
        marker = "o" if r.axis == "raw" else "s"
        lo = r["median"] - r.q10 if np.isfinite(r.q10) else 0
        hi = r.q90 - r["median"] if np.isfinite(r.q90) else 0
        ax.errorbar(
            r["median"],
            i,
            xerr=[[lo], [hi]],
            fmt=marker,
            color=color,
            mec="black",
            ms=5,
            capsize=2,
            lw=1,
        )
    style.reference_line(ax, 0, axis="x")
    lab = [
        f"{'arch1–Kaden' if 'arch1' in r.dataset else 'Kaden–Rep. RPE1'} · "
        f"{'resp.' if r.axis == 'raw' else 'β-like'} · "
        f"{'Arc' if r.label == 'Arc targets' else 'shared'} ({int(r.e25)})"
        for r in d.itertuples()
    ]
    ax.set_yticks(range(len(d)), lab, fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlim(-0.3, 1.05)
    ax.set_xlabel(
        "median per-perturbation Pearson (bootstrap 95% CI)\n(grey bar = median noise ceiling)"
    )
    ax.set_title("Per-perturbation agreement\non shared perturbations", fontsize=8.5)
    style.panel_label(ax, "D")

    common.save_figure(fig, "fig7_source_reliability")


if __name__ == "__main__":
    main()
