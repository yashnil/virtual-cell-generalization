"""Figures for the representation falsification battery."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
OUT = REPO_ROOT / "outputs" / "pathway_falsification_v1"
ORDER = ["K562", "RPE1", "HepG2", "Jurkat"]
CCOL = {"K562": "#1f77b4", "RPE1": "#d62728", "HepG2": "#2ca02c", "Jurkat": "#9467bd"}
NULLCOL = {"permuted": "#8C9196", "resampled": "#B8BDC2", "gaussian": "#D7DBDF"}


def main() -> None:
    obs = pd.read_csv(OUT / "observed_representations.csv")
    nulls = pd.read_csv(OUT / "null_replicates.csv")
    subsets = pd.read_csv(OUT / "source_subsets.csv")
    cross = pd.read_csv(REPO_ROOT / "outputs" / "zero_shot_v1" / "gamma_cross_context.csv")

    # --- 1. Hallmark vs null distributions -------------------------------
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.2), sharey=True)
    for ax, cl in zip(axes, ORDER, strict=True):
        for i, kind in enumerate(["permuted", "resampled", "gaussian"]):
            v = nulls[(nulls.cell_line == cl) & (nulls.representation == kind)][
                "r_gamma_normalised"
            ].dropna()
            parts = ax.violinplot(v, positions=[i], widths=0.8, showmedians=True)
            for b in parts["bodies"]:
                b.set_facecolor(NULLCOL[kind])
                b.set_alpha(0.9)
        h = obs[(obs.cell_line == cl) & (obs.representation == "hallmark")][
            "r_gamma_normalised"
        ].iloc[0]
        g = obs[(obs.cell_line == cl) & (obs.representation == "genes")]["r_gamma_normalised"].iloc[
            0
        ]
        ax.axhline(h, color=CCOL[cl], lw=2.5, label=f"Hallmark {h:.3f}")
        ax.axhline(g, color="k", lw=1.2, ls=":", label=f"genes {g:.3f}")
        ax.set_xticks(range(3))
        ax.set_xticklabels(["permuted", "resampled", "gaussian"], fontsize=8, rotation=20)
        ax.set_title(cl)
        ax.legend(fontsize=7, loc="upper left")
    axes[0].set_ylabel("reliability-normalised gamma recovery")
    fig.suptitle(
        "Hallmark vs matched nulls — biology beats aggregation in 3 of 4 contexts; "
        "HepG2 sits inside the null"
    )
    fig.tight_layout()
    fig.savefig(OUT / "fig1_hallmark_vs_null.png", dpi=140)
    plt.close(fig)

    # --- 2. representation comparison ------------------------------------
    fig, ax = plt.subplots(figsize=(11, 4.6))
    reps = ["genes", "permuted-null", "gaussian-null", "reactome", "hallmark"]
    x = np.arange(len(ORDER))
    w = 0.16
    for j, rep in enumerate(reps):
        vals = []
        for cl in ORDER:
            if rep.endswith("-null"):
                kind = rep.split("-")[0]
                vals.append(
                    nulls[(nulls.cell_line == cl) & (nulls.representation == kind)][
                        "r_gamma_normalised"
                    ].mean()
                )
            else:
                vals.append(
                    obs[(obs.cell_line == cl) & (obs.representation == rep)][
                        "r_gamma_normalised"
                    ].iloc[0]
                )
        col = {
            "genes": "#C9CED6",
            "permuted-null": "#8C9196",
            "gaussian-null": "#D7DBDF",
            "reactome": "#6F5CC4",
            "hallmark": "#D26A3A",
        }[rep]
        ax.bar(x + (j - 2) * w, vals, width=w, label=rep, color=col)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(ORDER)
    ax.set_ylabel("reliability-normalised gamma recovery")
    ax.set_title(
        "Random aggregation matches gene level exactly; only real pathways move the number"
    )
    ax.legend(fontsize=8, ncol=5)
    fig.tight_layout()
    fig.savefig(OUT / "fig2_representation_comparison.png", dpi=140)
    plt.close(fig)

    # --- 3. dependence: all 2-source subsets -----------------------------
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.0), sharey=True)
    for ax, cl in zip(axes, ORDER, strict=True):
        g = subsets[subsets.target == cl].copy()
        three = g[g.n_sources == 3]["r_gamma"].iloc[0]
        two = g[g.n_sources == 2].sort_values("r_gamma", ascending=False)
        cols = [
            "#D26A3A" if d in ("K562", "Jurkat") and cl in ("Jurkat", "K562") else "#8C9196"
            for d in two.dropped
        ]
        ax.bar(range(len(two)), two.r_gamma, color=cols)
        ax.axhline(three, color="k", ls="--", lw=1.5, label=f"all 3 sources ({three:+.3f})")
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xticks(range(len(two)))
        ax.set_xticklabels([f"drop\n{d}" for d in two.dropped], fontsize=8)
        ax.set_title(f"target {cl}")
        ax.legend(fontsize=7)
    axes[0].set_ylabel("Hallmark gamma recovery")
    fig.suptitle(
        "E. Dropping one source: Jurkat collapses without K562; K562 retains signal without Jurkat"
    )
    fig.tight_layout()
    fig.savefig(OUT / "fig3_source_dependence.png", dpi=140)
    plt.close(fig)

    # --- 4. dataset ancestry ---------------------------------------------
    origin = {"K562": "Replogle", "RPE1": "Replogle", "HepG2": "Nadig", "Jurkat": "Nadig"}
    cross["same_dataset"] = [
        origin[a] == origin[b] for a, b in zip(cross.cell_line_a, cross.cell_line_b, strict=True)
    ]
    cross["pair"] = cross.cell_line_a + "-" + cross.cell_line_b
    cross = cross.sort_values("excess_over_null", ascending=False)
    fig, ax = plt.subplots(figsize=(8, 4.4))
    cols = ["#d62728" if s else "#4C78A8" for s in cross.same_dataset]
    ax.barh(range(len(cross)), cross.excess_over_null, color=cols)
    ax.set_yticks(range(len(cross)))
    ax.set_yticklabels(cross.pair, fontsize=9)
    ax.invert_yaxis()
    ax.axvline(0, color="k", lw=1)
    ax.set_xlabel("cross-context gamma correlation, excess over the forced null")
    ax.set_title(
        "Same-dataset pairs (red) are the WEAKEST — dataset ancestry\n"
        "does not explain gamma sharing"
    )
    fig.tight_layout()
    fig.savefig(OUT / "fig4_dataset_ancestry.png", dpi=140)
    plt.close(fig)
    print("wrote 4 figures to", OUT)


if __name__ == "__main__":
    main()
