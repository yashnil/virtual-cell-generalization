"""Figure 6 — direct public perturbation evidence for the 300 Arc validation
targets (frozen identifier-presence tiers)."""

from __future__ import annotations

import matplotlib.pyplot as plt

from virtual_cell.visualization import common, style

PROV_HATCH = {
    "both arch1 and Kaden": "xx",
    "arch1 only": "//",
    "Kaden only": "..",
    "no public perturbation": "",
}


def main() -> None:
    style.apply()
    t = common.load_source("fig6_arc_target_support")
    total = int(t.panel_size.iloc[0])
    assert int(t.n_targets.sum()) == total == 300
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(style.FULL_WIDTH, 3.0))
    left = 0
    handles = []
    for tier in (2, 1, 0):
        n = int(t[t.tier == tier].n_targets.sum())
        meta = style.TIERS[tier]
        h = ax.barh(
            0,
            n,
            left=left,
            color=meta["color"],
            hatch=meta["hatch"],
            edgecolor="black",
            lw=0.6,
            height=0.6,
        )
        handles.append((h, f"{meta['label']}: {n} ({100 * n / total:.1f}%)"))
        if n >= 60:
            ax.text(
                left + n / 2,
                0,
                f"{meta['label']}\n{n} ({100 * n / total:.1f}%)",
                ha="center",
                va="center",
                fontsize=8,
                bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none", "alpha": 0.85},
            )
        left += n
    ax.set_xlim(0, total)
    ax.set_ylim(-0.5, 0.5)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xticks([0, 50, 100, 150, 200, 250, 300])
    ax.legend(
        *zip(*handles, strict=True),
        loc="upper center",
        bbox_to_anchor=(0.5, -0.35),
        ncol=3,
        fontsize=7.5,
    )
    ax.set_title(
        f"Arc validation panel: {total} targets by direct public perturbation evidence", fontsize=9
    )

    sup = t[t.tier > 0]
    left = 0
    handles = []
    for r in sup.itertuples():
        h = ax2.barh(
            0,
            r.n_targets,
            left=left,
            color=style.TIERS[r.tier]["color"],
            hatch=PROV_HATCH[r.provenance],
            edgecolor="black",
            lw=0.6,
            height=0.6,
        )
        handles.append((h, f"{r.provenance} ({style.TIERS[r.tier]['label']}): {r.n_targets}"))
        left += r.n_targets
    ax2.set_xlim(0, sup.n_targets.sum())
    ax2.set_ylim(-0.5, 0.5)
    ax2.set_yticks([])
    ax2.spines["left"].set_visible(False)
    ax2.set_title(
        f"The {int(sup.n_targets.sum())} supported targets, by source dataset", fontsize=9
    )
    ax2.legend(
        *zip(*handles, strict=True),
        loc="upper center",
        bbox_to_anchor=(0.5, -0.35),
        ncol=3,
        fontsize=7.5,
    )
    fig.text(
        0.01,
        0.005,
        "Tier 2 = perturbed in ≥2 public datasets; Tier 1 = exactly 1; Tier 0 = "
        "none. The frozen model gives Tier 0 zero perturbation-specific effect.",
        fontsize=6.5,
        color=style.GREY,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1), h_pad=1.5)
    common.save_figure(fig, "fig6_arc_target_support")


if __name__ == "__main__":
    main()
