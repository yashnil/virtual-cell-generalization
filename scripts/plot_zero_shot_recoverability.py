"""Figures for the zero-shot recoverability diagnostic."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
OUT = REPO_ROOT / "outputs" / "zero_shot_v1"
ORDER = ["K562", "RPE1", "HepG2", "Jurkat"]
CCOL = {"K562": "#1f77b4", "RPE1": "#d62728", "HepG2": "#2ca02c", "Jurkat": "#9467bd"}
BASELINES = ["source_mean", "basal_affine", "basal_simplex", "nearest_basal"]


def main() -> None:
    met = pd.read_csv(OUT / "loco_metrics.csv")
    rel = pd.read_csv(OUT / "loco_reliability.csv")
    gam = pd.read_csv(OUT / "loco_gamma.csv")
    desc = pd.read_csv(OUT / "transferability_descriptives.csv")
    pairs = pd.read_csv(OUT / "gamma_cross_context.csv")

    # --- 1. transfer performance by held-out context ---------------------
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.4))
    sub = met[met.baseline.isin(BASELINES)]
    for ax, metric, title in zip(
        axes,
        ["pearson", "pearson_template_removed", "energy_explained"],
        ["Pearson vs held-out response", "Template-removed Pearson", "Response energy explained"],
        strict=True,
    ):
        data, labels, colors = [], [], []
        for cl in ORDER:
            for b in BASELINES:
                v = sub[(sub.cell_line == cl) & (sub.baseline == b)][metric].dropna()
                data.append(v)
                labels.append(f"{cl}\n{b.replace('_', ' ')}")
                colors.append(CCOL[cl])
        bp = ax.boxplot(data, showfliers=False, patch_artist=True, widths=0.65)
        for patch, c in zip(bp["boxes"], colors, strict=True):
            patch.set_facecolor(c)
            patch.set_alpha(0.55)
        ax.set_xticklabels(labels, fontsize=5.5, rotation=90)
        ax.axhline(0, color="k", lw=0.8, ls="--")
        ax.set_title(title, fontsize=10)
        if metric == "energy_explained":
            ax.set_ylim(-2.0, 1.0)
    axes[0].set_ylabel("per-perturbation value")
    fig.suptitle(
        "Zero-shot conserved transfer by held-out context "
        "(energy explained is mostly NEGATIVE: right direction, wrong magnitude)"
    )
    fig.tight_layout()
    fig.savefig(OUT / "fig1_transfer_by_context.png", dpi=140)
    plt.close(fig)

    # --- 2. observed vs reliability-normalised ---------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))
    sm = rel[rel.baseline == "source_mean"]
    x = np.arange(len(ORDER))
    obs = [np.nanmedian(sm[sm.cell_line == c]["r_half_mean"]) for c in ORDER]
    dis = [np.nanmedian(sm[sm.cell_line == c]["r_disattenuated"]) for c in ORDER]
    ceil = [np.nanmedian(sm[sm.cell_line == c]["ceiling"]) for c in ORDER]
    axes[0].bar(x - 0.26, obs, width=0.25, label="observed (vs one half)", color="#4C78A8")
    axes[0].bar(x, ceil, width=0.25, label="ceiling = sqrt(rho_half)", color="#C9CED6")
    axes[0].bar(x + 0.26, dis, width=0.25, label="reliability-normalised", color="#6F5CC4")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(ORDER)
    axes[0].set_ylabel("median Pearson")
    axes[0].set_title("Conserved transfer: observed vs reliability-normalised")
    axes[0].legend(fontsize=8)
    axes[0].axhline(1.0, color="k", lw=0.8, ls=":")

    for cl in ORDER:
        g = sm[sm.cell_line == cl]
        axes[1].scatter(
            g["rho_half"],
            g["r_half_mean"],
            s=5,
            alpha=0.25,
            color=CCOL[cl],
            edgecolors="none",
            label=cl,
        )
    grid = np.linspace(0.01, 1, 200)
    axes[1].plot(grid, np.sqrt(grid), "k--", lw=1.5, label=r"ceiling $\sqrt{\rho}$")
    axes[1].set_xlabel(r"target split-half reliability $\rho_{half}$")
    axes[1].set_ylabel("observed r (prediction vs half)")
    axes[1].set_title("Nothing may exceed the measurement ceiling")
    axes[1].legend(fontsize=8, markerscale=2)
    fig.tight_layout()
    fig.savefig(OUT / "fig2_observed_vs_normalised.png", dpi=140)
    plt.close(fig)

    # --- 3. transfer error vs target reliability -------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))
    d = desc.copy()
    d["stratum"] = pd.cut(
        d.target_rho_half, [-1, 0.1, 0.3, 0.6, 1.0], labels=["<0.1", "0.1-0.3", "0.3-0.6", ">0.6"]
    )
    grp = d.groupby("stratum", observed=True)
    axes[0].bar(range(4), grp["pearson"].median(), color="#4C78A8", label="raw r", width=0.38)
    axes[0].bar(
        np.arange(4) + 0.4,
        grp["r_disattenuated"].median(),
        color="#6F5CC4",
        label="reliability-normalised",
        width=0.38,
    )
    axes[0].set_xticks(np.arange(4) + 0.2)
    axes[0].set_xticklabels(["<0.1", "0.1-0.3", "0.3-0.6", ">0.6"])
    axes[0].set_xlabel(r"target reliability $\rho_{half}$")
    axes[0].set_ylabel("median Pearson")
    axes[0].set_title("Apparent failure at low reliability is mostly measurement noise")
    axes[0].legend(fontsize=8)

    for cl in ORDER:
        g = d[d.cell_line == cl]
        axes[1].scatter(
            g["target_rho_half"],
            g["mse"] if "mse" in g else g["obs_norm"],
            s=5,
            alpha=0.2,
            color=CCOL[cl],
            edgecolors="none",
            label=cl,
        )
    axes[1].set_xlabel(r"target reliability $\rho_{half}$")
    axes[1].set_ylabel(r"$\|\delta_{target}\|$")
    axes[1].set_title("Response magnitude vs reliability")
    axes[1].legend(fontsize=8, markerscale=2)
    fig.tight_layout()
    fig.savefig(OUT / "fig3_error_vs_reliability.png", dpi=140)
    plt.close(fig)

    # --- 4. gamma recoverability by context ------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
    gb = gam[gam.baseline.isin(["basal_affine", "nearest_basal"])]
    x = np.arange(len(ORDER))
    for off, b, col in ((-0.2, "basal_affine", "#6F5CC4"), (0.2, "nearest_basal", "#D26A3A")):
        med = [np.nanmedian(gb[(gb.cell_line == c) & (gb.baseline == b)]["r_gamma"]) for c in ORDER]
        axes[0].bar(x + off, med, width=0.38, label=b.replace("_", " "), color=col)
    ceil_g = [np.nanmedian(gam[gam.cell_line == c]["gamma_ceiling_full"]) for c in ORDER]
    axes[0].plot(x, ceil_g, "k_", markersize=26, label=r"ceiling $\sqrt{\rho_{full}}$")
    axes[0].axhline(0, color="k", lw=0.8)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(ORDER)
    axes[0].set_ylabel(r"median $r(\gamma_{true}, \hat{\gamma})$")
    axes[0].set_title("Zero-shot gamma recovery: real for K562/Jurkat, absent for RPE1/HepG2")
    axes[0].legend(fontsize=8)

    for cl in ORDER:
        v = gb[(gb.cell_line == cl) & (gb.baseline == "basal_affine")]["r_gamma"].dropna()
        axes[1].hist(v, bins=60, histtype="step", lw=1.6, label=cl, color=CCOL[cl])
    axes[1].axvline(0, color="k", lw=1, ls="--")
    axes[1].set_xlabel(r"$r(\gamma_{true}, \hat{\gamma})$ per perturbation")
    axes[1].set_ylabel("perturbations")
    axes[1].set_title("Distribution of gamma recovery (basal-affine)")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig4_gamma_recoverability.png", dpi=140)
    plt.close(fig)

    # --- 5. source agreement vs transfer success -------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))
    for cl in ORDER:
        g = desc[desc.cell_line == cl]
        axes[0].scatter(
            g["source_agreement"],
            g["pearson"],
            s=5,
            alpha=0.25,
            color=CCOL[cl],
            edgecolors="none",
            label=cl,
        )
        axes[1].scatter(
            g["source_agreement"],
            g["r_disattenuated"],
            s=5,
            alpha=0.25,
            color=CCOL[cl],
            edgecolors="none",
            label=cl,
        )
    rho_raw = desc["pearson"].corr(desc["source_agreement"], method="spearman")
    rho_dis = desc["r_disattenuated"].corr(desc["source_agreement"], method="spearman")
    axes[0].set_title(f"Source agreement vs transfer success (Spearman {rho_raw:+.3f})")
    axes[1].set_title(f"...after reliability normalisation (Spearman {rho_dis:+.3f})")
    for ax, yl in zip(axes, ["raw Pearson", "reliability-normalised Pearson"], strict=True):
        ax.set_xlabel("mean pairwise correlation among the 3 source responses")
        ax.set_ylabel(yl)
        ax.legend(fontsize=8, markerscale=2)
        ax.axhline(0, color="k", lw=0.6, ls=":")
    fig.suptitle("Source agreement is computable at inference time and is the strongest predictor")
    fig.tight_layout()
    fig.savefig(OUT / "fig5_source_agreement.png", dpi=140)
    plt.close(fig)

    # --- 6. basal similarity vs gamma similarity -------------------------
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(
        pairs["basal_similarity"],
        pairs["median_gamma_r"],
        s=80,
        color="#4C78A8",
        edgecolors="k",
        zorder=3,
    )
    for _, r in pairs.iterrows():
        ax.annotate(
            f"{r['cell_line_a']}-{r['cell_line_b']}",
            (r["basal_similarity"], r["median_gamma_r"]),
            textcoords="offset points",
            xytext=(6, 5),
            fontsize=8,
        )
    ax.axhline(
        pairs["null_expectation"].iloc[0],
        color="#d62728",
        ls="--",
        label=r"null $-1/(C-1)$ forced by $\sum_c \gamma = 0$",
    )
    ax.set_xlabel("basal (control) similarity")
    ax.set_ylabel(r"median cross-context $r(\gamma)$")
    rho = pairs["basal_similarity"].corr(pairs["excess_over_null"], method="spearman")
    ax.set_title(f"Basal similarity vs gamma sharing (Spearman {rho:+.3f}, n=6 pairs)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig6_basal_vs_gamma_similarity.png", dpi=140)
    plt.close(fig)

    # --- 7. example classes ----------------------------------------------
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    reliable = desc[desc.target_rho_half > 0.5]
    ax.scatter(
        desc["target_rho_half"],
        desc["r_disattenuated"],
        s=4,
        alpha=0.12,
        color="#999999",
        edgecolors="none",
        label="all pairs",
    )
    top = reliable.nlargest(8, "r_disattenuated")
    bot = reliable.nsmallest(8, "r_disattenuated")
    noisy = desc[desc.target_rho_half < 0.05].nlargest(8, "obs_norm")
    for frame, col, lab in (
        (top, "#6F5CC4", "transferable"),
        (bot, "#D26A3A", "reliably context-specific"),
        (noisy, "#C9A227", "large but unreliable"),
    ):
        ax.scatter(
            frame["target_rho_half"],
            frame["r_disattenuated"],
            s=55,
            color=col,
            edgecolors="k",
            linewidths=0.5,
            label=lab,
            zorder=3,
        )
        for _, r in frame.iterrows():
            ax.annotate(
                f"{r['perturbation']}",
                (r["target_rho_half"], r["r_disattenuated"]),
                textcoords="offset points",
                xytext=(4, 3),
                fontsize=6,
            )
    ax.axhline(0, color="k", lw=0.8, ls="--")
    ax.set_xlabel(r"target reliability $\rho_{half}$")
    ax.set_ylabel("reliability-normalised transfer r")
    ax.set_title("Three distinct perturbation classes")
    ax.legend(fontsize=8, markerscale=1.2)
    fig.tight_layout()
    fig.savefig(OUT / "fig7_example_classes.png", dpi=140)
    plt.close(fig)
    print("wrote 7 figures to", OUT)


if __name__ == "__main__":
    main()
