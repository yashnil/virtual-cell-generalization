"""Figure 3 — gene-level vs Hallmark pathway-level zero-shot interaction
recovery, against a set-structure-preserving permutation null."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from virtual_cell.visualization import common, style


def main() -> None:
    style.apply()
    t = common.load_source("fig3_gamma_recovery")
    obs = t[t.record == "observed"]
    nulls = t[t.record == "hallmark_permuted_null"]
    tests = t[t.record == "hallmark_test"].set_index("cell_line")
    rtests = t[t.record == "reactome_test"].set_index("cell_line")
    fig, axes = plt.subplots(
        1, 2, figsize=(style.FULL_WIDTH, 3.2), gridspec_kw={"width_ratios": [1.6, 1]}
    )
    ax = axes[0]
    rng = np.random.default_rng(0)
    for i, cl in enumerate(style.CONTEXT_ORDER):
        o = obs[obs.cell_line == cl].set_index("representation")
        nv = nulls[nulls.cell_line == cl].r_gamma_normalised.to_numpy()
        ax.scatter(
            i + 0.18 + rng.uniform(-0.06, 0.06, len(nv)),
            nv,
            s=5,
            color=style.GREY,
            alpha=0.5,
            lw=0,
            label="matched-random null (100 permutations)" if i == 0 else None,
        )
        ax.plot(
            i - 0.18,
            o.loc["genes", "r_gamma_normalised"],
            marker="o",
            ms=7,
            mfc="white",
            mec="black",
            ls="",
            label="gene level" if i == 0 else None,
        )
        ax.plot(
            i,
            o.loc["hallmark", "r_gamma_normalised"],
            marker="D",
            ms=7,
            color=style.COMPONENTS["gamma"]["color"],
            mec="black",
            ls="",
            label="Hallmark pathways" if i == 0 else None,
        )
        ax.plot(
            [i - 0.18, i],
            [o.loc["genes", "r_gamma_normalised"], o.loc["hallmark", "r_gamma_normalised"]],
            color=style.GREY,
            lw=0.8,
        )
        p = tests.loc[cl, "p_value"]
        z = tests.loc[cl, "z"]
        ax.text(i, 0.74, f"z = {z:+.1f}\np = {p:.3f}", ha="center", fontsize=7)
    style.reference_line(ax, 0)
    ax.set_xticks(range(4), style.CONTEXT_ORDER)
    ax.set_xlabel("held-out context")
    ax.set_ylabel("reliability-normalised\nr(γ true, γ recovered)")
    ax.set_ylim(-0.15, 0.85)
    ax.set_title("Zero-shot γ recovery: genes vs Hallmark pathways", fontsize=9)
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.86), fontsize=7)
    style.panel_label(ax, "A")

    ax = axes[1]
    for i, cl in enumerate(style.CONTEXT_ORDER):
        o = obs[obs.cell_line == cl].set_index("representation")
        v = o.loc["reactome", "r_gamma_normalised"]
        nm, ns = rtests.loc[cl, "null_mean"], rtests.loc[cl, "null_sd"]
        ax.errorbar(i + 0.15, nm, yerr=2 * ns, fmt="none", ecolor=style.GREY, capsize=3)
        ax.plot(
            i + 0.15,
            nm,
            marker="_",
            ms=10,
            color=style.GREY,
            label="Reactome null mean ± 2 sd (20 perm.)" if i == 0 else None,
        )
        ax.plot(
            i - 0.1,
            v,
            marker="s",
            ms=6,
            color=style.SKY,
            mec="black",
            ls="",
            label="Reactome pathways" if i == 0 else None,
        )
    style.reference_line(ax, 0)
    ax.set_xticks(range(4), style.CONTEXT_ORDER)
    ax.set_ylim(-0.15, 0.85)
    ax.set_title("Reactome (independent ontology)", fontsize=9)
    ax.legend(loc="upper right", fontsize=7)
    style.panel_label(ax, "B")
    fig.tight_layout()
    common.save_figure(fig, "fig3_gamma_recovery")


if __name__ == "__main__":
    main()
