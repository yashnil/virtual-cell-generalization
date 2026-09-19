"""Figures for the four-context robustness battery."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
OUT = REPO_ROOT / "outputs" / "four_context_sensitivity"
COLORS = {"template": "#3B5BA9", "beta": "#6F5CC4", "gamma": "#D26A3A", "noise": "#C9CED6"}
DEPTHS = [15, 30, 50, 100]


def variant_label(v: dict) -> str:
    scheme = "split-ctrl" if v["control_scheme"] == "split" else "shared-ctrl"
    agg = "log(mean)" if v["aggregation"] == "log_mean" else "mean(log)"
    return f"{scheme}\n{agg}\nseed {v['seed']}"


def main() -> None:
    variants = json.loads((OUT / "variants.json").read_text())
    depth = pd.read_csv(OUT / "depth_experiment.csv")
    within = pd.read_csv(OUT / "depth_within_pair.csv")

    # --- 1. variant comparison across gene subsets -----------------------
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6), sharey=True)
    subsets = ["all", "hvg4000", "hvg2000"]
    titles = {
        "all": "all 6,640 shared genes",
        "hvg4000": "top 4,000 HVG",
        "hvg2000": "top 2,000 HVG",
    }
    labels = [variant_label(v) for v in variants]
    x = np.arange(len(variants))
    for ax, sub in zip(axes, subsets, strict=True):
        bottom = np.zeros(len(variants))
        for comp in ("template", "beta", "gamma", "noise"):
            vals = np.array(
                [v["subsets"][sub]["corrected_fractions"][comp] * 100 for v in variants]
            )
            ax.bar(x, vals, bottom=bottom, color=COLORS[comp], label=comp, width=0.72)
            for xi, (b, val) in enumerate(zip(bottom, vals, strict=True)):
                if val > 4:
                    ax.text(
                        xi,
                        b + val / 2,
                        f"{val:.1f}",
                        ha="center",
                        va="center",
                        fontsize=7,
                        color="white" if comp != "noise" else "#1F2937",
                    )
            bottom += vals
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=6.5)
        ax.set_title(titles[sub])
        ax.set_ylim(0, 100)
    axes[0].set_ylabel("% of response energy (noise-corrected)")
    axes[-1].legend(loc="upper right", fontsize=8, framealpha=0.95)
    fig.suptitle(
        "Robustness battery: component shares across every variant "
        "(independent four-context decomposition)"
    )
    fig.tight_layout()
    fig.savefig(OUT / "variant_comparison.png", dpi=140)
    plt.close(fig)

    # --- 2. beta / gamma stability summary -------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    for ax, comp in zip(axes, ("beta", "gamma"), strict=True):
        for si, sub in enumerate(subsets):
            vals = [v["subsets"][sub]["corrected_fractions"][comp] * 100 for v in variants]
            ax.scatter(
                [si] * len(vals),
                vals,
                s=45,
                alpha=0.75,
                color=COLORS[comp],
                edgecolors="k",
                linewidths=0.4,
            )
        ax.set_xticks(range(len(subsets)))
        ax.set_xticklabels([titles[s] for s in subsets], fontsize=8)
        ax.set_ylabel(f"{comp} share (%)")
        ax.set_title(f"{comp}: every variant, every feature space")
        ax.axhline(0, color="#d62728", ls="--", lw=1)
        lo = min(
            v["subsets"][s]["corrected_fractions"][comp] * 100 for v in variants for s in subsets
        )
        hi = max(
            v["subsets"][s]["corrected_fractions"][comp] * 100 for v in variants for s in subsets
        )
        ax.set_ylim(0, max(hi * 1.25, 5))
        ax.annotate(
            f"range {lo:.2f}–{hi:.2f}%", xy=(0.03, 0.92), xycoords="axes fraction", fontsize=9
        )
    fig.suptitle("Neither component approaches zero under any reasonable variant")
    fig.tight_layout()
    fig.savefig(OUT / "beta_gamma_stability.png", dpi=140)
    plt.close(fig)

    # --- 3. controlled depth experiment ----------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.4))
    med = depth.groupby("n_cells")["reliability"].median()
    rng = np.random.default_rng(0)
    los, his = [], []
    for n in DEPTHS:
        v = depth.loc[depth.n_cells == n, "reliability"].to_numpy()
        boot = [np.median(rng.choice(v, size=len(v), replace=True)) for _ in range(2000)]
        los.append(np.percentile(boot, 2.5))
        his.append(np.percentile(boot, 97.5))
    axes[0].errorbar(
        DEPTHS,
        med.loc[DEPTHS],
        yerr=[med.loc[DEPTHS] - los, his - med.loc[DEPTHS]],
        marker="o",
        lw=2,
        capsize=5,
        color="#1f77b4",
    )
    axes[0].set_xscale("log")
    axes[0].set_xticks(DEPTHS)
    axes[0].set_xticklabels(DEPTHS)
    axes[0].set_xlabel("cells per half")
    axes[0].set_ylabel("median split-half reliability")
    axes[0].set_title("Same 643 pairs, re-estimated at each depth\n(95% bootstrap CI)")

    piv = depth.pivot_table(
        index="cell_line", columns="n_cells", values="reliability", aggfunc="median"
    )
    for cl in piv.index:
        axes[1].plot(DEPTHS, piv.loc[cl, DEPTHS], marker="o", label=cl)
    axes[1].set_xscale("log")
    axes[1].set_xticks(DEPTHS)
    axes[1].set_xticklabels(DEPTHS)
    axes[1].set_xlabel("cells per half")
    axes[1].set_ylabel("median reliability")
    axes[1].set_title("Every context improves with depth")
    axes[1].legend(fontsize=8)

    sample = within.sample(n=min(250, len(within)), random_state=0)
    for _, row in sample.iterrows():
        axes[2].plot(DEPTHS, [row[str(d)] for d in DEPTHS], color="#888888", alpha=0.12, lw=0.8)
    axes[2].plot(
        DEPTHS,
        [within[str(d)].median() for d in DEPTHS],
        color="#d62728",
        lw=2.5,
        marker="o",
        label="median",
    )
    axes[2].set_xscale("log")
    axes[2].set_xticks(DEPTHS)
    axes[2].set_xticklabels(DEPTHS)
    axes[2].set_xlabel("cells per half")
    axes[2].set_ylabel("reliability")
    axes[2].set_title("Within-pair trajectories (250 sampled)\n99.5% improve from n=15 to n=100")
    axes[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "depth_experiment.png", dpi=140)
    plt.close(fig)

    # --- 4. seed stability ------------------------------------------------
    seeds = [
        v for v in variants if v["control_scheme"] == "shared" and v["aggregation"] == "mean_log"
    ]
    fig, ax = plt.subplots(figsize=(7, 4.2))
    for comp in ("template", "beta", "gamma", "noise"):
        vals = [v["subsets"]["all"]["corrected_fractions"][comp] * 100 for v in seeds]
        ax.scatter(
            [comp] * len(vals),
            vals,
            s=60,
            color=COLORS[comp],
            edgecolors="k",
            linewidths=0.4,
            zorder=3,
        )
        ax.annotate(
            f"range {max(vals) - min(vals):.3f} pp",
            xy=(comp, max(vals)),
            xytext=(0, 9),
            textcoords="offset points",
            ha="center",
            fontsize=8,
        )
    ax.set_ylabel("% of response energy")
    ax.set_title(f"Seed stability: {len(seeds)} independent seeds, all 6,640 genes")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "seed_stability.png", dpi=140)
    plt.close(fig)
    print("wrote 4 figures to", OUT)


if __name__ == "__main__":
    main()
