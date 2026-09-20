"""Figures for the pathway residual model v1."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "outputs" / "pathway_residual_v1"
ORDER = ["K562", "RPE1", "HepG2", "Jurkat"]
CCOL = {"K562": "#1f77b4", "RPE1": "#d62728", "HepG2": "#2ca02c", "Jurkat": "#9467bd"}


def main() -> None:
    folds = pd.read_csv(OUT / "hallmark_folds.csv")
    sweep = pd.read_csv(OUT / "oracle_lambda_sweep.csv")
    diag = pd.read_csv(OUT / "forced_m2_diagnostics.csv")
    conf = pd.read_csv(OUT / "confidence_strata.csv")
    nullc = pd.read_csv(OUT / "null_comparison.csv")

    # --- 1. baseline vs corrected ----------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))
    x = np.arange(len(ORDER))
    for ax, (b, c, lab) in zip(
        axes,
        [
            ("base_pearson", "corr_pearson", "median Pearson"),
            ("base_energy_explained", "corr_energy_explained", "energy explained"),
        ],
        strict=True,
    ):
        bv = [folds.loc[folds.cell_line == cl, b].iloc[0] for cl in ORDER]
        cv = [folds.loc[folds.cell_line == cl, c].iloc[0] for cl in ORDER]
        ax.bar(x - 0.2, bv, width=0.38, label="baseline B", color="#C9CED6")
        ax.bar(x + 0.2, cv, width=0.38, label="B + lambda*R_hat", color="#D26A3A")
        for i, cl in enumerate(ORDER):
            lam = folds.loc[folds.cell_line == cl, "selected_lambda"].iloc[0]
            fam = folds.loc[folds.cell_line == cl, "selected_family"].iloc[0]
            ax.annotate(
                f"{fam}\nlam={lam:g}",
                (i, max(bv[i], cv[i])),
                textcoords="offset points",
                xytext=(0, 5),
                ha="center",
                fontsize=7,
            )
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(ORDER)
        ax.set_ylabel(lab)
        ax.legend(fontsize=8)
    fig.suptitle("The learned pathway correction never improves the outer target")
    fig.tight_layout()
    fig.savefig(OUT / "fig1_baseline_vs_corrected.png", dpi=140)
    plt.close(fig)

    # --- 2. oracle lambda sweep ------------------------------------------
    fig, axes = plt.subplots(1, 4, figsize=(17, 3.9), sharey=False)
    for ax, cl in zip(axes, ORDER, strict=True):
        g = sweep[sweep.cell_line == cl]
        for fam, sub in g.groupby("family"):
            ax.plot(sub.lam, sub.pearson, marker="o", ms=4, label=fam)
        base = g[g.lam == 0].pearson.iloc[0]
        ax.axhline(base, color="k", ls=":", lw=1.2)
        ax.set_title(cl)
        ax.set_xlabel("lambda")
        ax.legend(fontsize=6)
    axes[0].set_ylabel("median Pearson (outer)")
    fig.suptitle(
        "Oracle lambda sweep (EVALUATION-ONLY): lambda = 0 is essentially "
        "optimal everywhere — max attainable gain +0.003"
    )
    fig.tight_layout()
    fig.savefig(OUT / "fig2_oracle_lambda.png", dpi=140)
    plt.close(fig)

    # --- 3. what did R_hat learn -----------------------------------------
    fig, ax = plt.subplots(figsize=(8.5, 4.4))
    x = np.arange(len(ORDER))
    g = [diag.loc[diag.cell_line == cl, "r_Rhat_gamma"].iloc[0] for cl in ORDER]
    b = [diag.loc[diag.cell_line == cl, "r_Rhat_beta"].iloc[0] for cl in ORDER]
    d = [diag.loc[diag.cell_line == cl, "r_detHat_gamma"].iloc[0] for cl in ORDER]
    ax.bar(x - 0.26, g, width=0.25, label=r"$r(\hat R,\ \gamma)$ learned", color="#D26A3A")
    ax.bar(x, b, width=0.25, label=r"$r(\hat R,\ \beta$-like$)$ learned", color="#6F5CC4")
    ax.bar(
        x + 0.26,
        d,
        width=0.25,
        label=r"$r(\hat\gamma_{det},\ \gamma)$ deterministic",
        color="#C9CED6",
    )
    ax.axhline(0, color="k", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels(ORDER)
    ax.set_ylabel("median correlation")
    ax.set_title(
        "The correction finds gamma in K562/Jurkat — but predicts beta with the "
        "WRONG SIGN,\nwhich cancels the benefit (forced M2)"
    )
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig3_residual_diagnostics.png", dpi=140)
    plt.close(fig)

    # --- 4. confidence + null --------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))
    for cl in ORDER:
        g = conf[conf.cell_line == cl].sort_values("quartile")
        axes[0].plot(g.quartile, g.base_r, marker="o", color=CCOL[cl], label=f"{cl} base")
        axes[0].plot(g.quartile, g.corr_r, marker="s", ls="--", color=CCOL[cl], alpha=0.6)
    axes[0].set_xlabel("source-agreement quartile")
    axes[0].set_ylabel("median Pearson")
    axes[0].set_title("Confidence stays monotone (solid = baseline, dashed = corrected)")
    axes[0].legend(fontsize=7)

    dn = nullc[nullc.metric == "delta_pearson"]
    xx = np.arange(len(ORDER))
    obs = [dn.loc[dn.cell_line == cl, "observed"].iloc[0] for cl in ORDER]
    nm = [dn.loc[dn.cell_line == cl, "null_mean"].iloc[0] for cl in ORDER]
    nsd = [dn.loc[dn.cell_line == cl, "null_sd"].iloc[0] for cl in ORDER]
    axes[1].bar(xx - 0.2, obs, width=0.38, label="Hallmark", color="#D26A3A")
    axes[1].bar(
        xx + 0.2, nm, width=0.38, yerr=nsd, capsize=4, label="matched-random", color="#8C9196"
    )
    axes[1].axhline(0, color="k", lw=1)
    axes[1].set_xticks(xx)
    axes[1].set_xticklabels(ORDER)
    axes[1].set_ylabel("change in median Pearson")
    axes[1].set_title("Hallmark modelling gain is indistinguishable from random\n(both are ~zero)")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig4_confidence_and_null.png", dpi=140)
    plt.close(fig)
    print("wrote 4 figures to", OUT)


if __name__ == "__main__":
    main()
