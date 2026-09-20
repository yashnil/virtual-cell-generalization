"""Figures for the transferability / confidence model v1."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "outputs" / "transferability_v1"
ORDER = ["K562", "RPE1", "HepG2", "Jurkat"]
SCOL = {
    "M0": "#D26A3A",
    "learned": "#6F5CC4",
    "random": "#C9CED6",
    "source_reliability": "#4C78A8",
    "source_magnitude": "#2ca02c",
}


def main() -> None:
    rc = pd.read_csv(OUT / "risk_coverage.csv")
    prio = pd.read_csv(OUT / "prioritisation.csv")
    cal = pd.read_csv(OUT / "calibration.csv")
    rk = pd.read_csv(OUT / "ranking.csv")

    # --- 1. risk-coverage on the normalised metric -----------------------
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.0), sharey=False)
    for ax, cl in zip(axes, ORDER, strict=True):
        for s in ["M0", "learned", "source_reliability", "source_magnitude", "random"]:
            g = rc[(rc.cell_line == cl) & (rc.score == s)].sort_values("coverage")
            ax.plot(
                g.coverage * 100,
                g.median_r_disatt,
                marker="o",
                ms=4,
                color=SCOL[s],
                label=s,
                lw=2 if s == "M0" else 1.2,
            )
        ax.invert_xaxis()
        ax.set_xlabel("coverage (%)")
        ax.set_title(cl)
    axes[0].set_ylabel("median reliability-normalised similarity")
    axes[0].legend(fontsize=7)
    fig.suptitle(
        "Selective prediction: raw source agreement (M0) is monotone in every context; "
        "the learned model does not beat it"
    )
    fig.tight_layout()
    fig.savefig(OUT / "fig1_risk_coverage.png", dpi=140)
    plt.close(fig)

    # --- 2. risk-coverage on D -------------------------------------------
    fig, axes = plt.subplots(1, 4, figsize=(17, 3.9))
    for ax, cl in zip(axes, ORDER, strict=True):
        for s in ["M0", "learned", "random"]:
            g = rc[(rc.cell_line == cl) & (rc.score == s)].sort_values("coverage")
            ax.plot(g.coverage * 100, g.median_D, marker="o", ms=4, color=SCOL[s], label=s)
        ax.axhline(1.0, color="k", ls=":", lw=1)
        ax.invert_xaxis()
        ax.set_xlabel("coverage (%)")
        ax.set_title(cl)
    axes[0].set_ylabel("median D  (lower = better)")
    axes[0].legend(fontsize=7)
    fig.suptitle(
        "D falls monotonically as coverage tightens, in all four contexts "
        "(dotted line: no better than predicting zero)"
    )
    fig.tight_layout()
    fig.savefig(OUT / "fig2_risk_coverage_D.png", dpi=140)
    plt.close(fig)

    # --- 3. M0 vs learned -------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.3))
    piv = rk.pivot_table(index="cell_line", columns="score", values="spearman_vs_negD")
    x = np.arange(len(ORDER))
    for off, s in ((-0.26, "M0"), (0.0, "M1"), (0.26, "learned")):
        axes[0].bar(
            x + off,
            [piv.loc[cl, s] for cl in ORDER],
            width=0.25,
            label=s,
            color=SCOL.get(s, "#8C9196"),
        )
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(ORDER)
    axes[0].set_ylabel(r"Spearman(score, $-D$)")
    axes[0].set_title("Ranking quality: the ladder adds nothing")
    axes[0].legend(fontsize=8)

    w = (
        rc[rc.score.isin(["M0", "learned"])]
        .pivot_table(index=["cell_line", "coverage"], columns="score", values="median_r_disatt")
        .reset_index()
    )
    w["diff"] = w["M0"] - w["learned"]
    for cl in ORDER:
        g = w[w.cell_line == cl].sort_values("coverage")
        axes[1].plot(g.coverage * 100, g["diff"], marker="o", ms=4, label=cl)
    axes[1].axhline(0, color="k", lw=1)
    axes[1].invert_xaxis()
    axes[1].set_xlabel("coverage (%)")
    axes[1].set_ylabel("M0 − learned")
    axes[1].set_title("M0 minus learned (positive = simple statistic wins)")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig3_m0_vs_learned.png", dpi=140)
    plt.close(fig)

    # --- 4. prioritisation ------------------------------------------------
    fig, axes = plt.subplots(1, 4, figsize=(17, 3.9))
    for ax, cl in zip(axes, ORDER, strict=True):
        for s in ["source_magnitude", "random", "M0"]:
            g = prio[(prio.cell_line == cl) & (prio.score == s)].sort_values("budget")
            lab = "rank by expected error" if s == "source_magnitude" else s
            ax.plot(
                g.budget * 100,
                g.captured_fraction * 100,
                marker="o",
                ms=4,
                color=SCOL[s],
                label=lab,
            )
        ax.plot([0, 50], [0, 50], "k:", lw=1)
        ax.set_xlabel("experiment budget (%)")
        ax.set_title(cl)
    axes[0].set_ylabel("% of reproducible error captured")
    axes[0].legend(fontsize=7)
    fig.suptitle(
        "Prioritisation: targeting the LEAST-CONFIDENT captures LESS than random — "
        "trustworthiness and error magnitude are different objectives"
    )
    fig.tight_layout()
    fig.savefig(OUT / "fig4_prioritisation.png", dpi=140)
    plt.close(fig)

    # --- 5. calibration ---------------------------------------------------
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    for cl in ORDER:
        g = cal[(cal.cell_line == cl) & (cal.score == "M0")].sort_values("bin")
        ax.plot(g.mean_confidence, g.mean_quality, marker="o", label=cl)
    ax.set_xlabel("mean source agreement in bin")
    ax.set_ylabel("observed reliability-normalised quality")
    ax.set_title("Calibration is monotone in every context (5 equal-count bins)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig5_calibration.png", dpi=140)
    plt.close(fig)
    print("wrote 5 figures to", OUT)


if __name__ == "__main__":
    main()
