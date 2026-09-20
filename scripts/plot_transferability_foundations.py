"""Figures for the transferability foundations study."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
OUT = REPO_ROOT / "outputs" / "transferability_foundations_v1"
ORDER = ["K562", "RPE1", "HepG2", "Jurkat"]
CCOL = {"K562": "#1f77b4", "RPE1": "#d62728", "HepG2": "#2ca02c", "Jurkat": "#9467bd"}


def main() -> None:
    diag = pd.read_csv(OUT / "template_diagnostics.csv")
    perf = pd.read_csv(OUT / "template_performance.csv")
    bsum = pd.read_csv(OUT / "source_agreement_summary.csv")
    feats = pd.read_csv(OUT / "source_agreement_features.csv")
    path = pd.read_csv(OUT / "pathway_gamma.csv")
    dcand = pd.read_csv(OUT / "candidate_d_targets.csv")
    gene_gamma = pd.read_csv(REPO_ROOT / "outputs" / "zero_shot_v1" / "loco_gamma.csv")

    # --- 1. template / scale correction ----------------------------------
    order = [
        "zero",
        "direct_basal",
        "global_scalar",
        "ridge_subspace",
        "scale_only",
        "oracle_alpha_EVAL_ONLY",
        "scale_plus_oracle_alpha_EVAL_ONLY",
    ]
    fig, ax = plt.subplots(figsize=(13, 4.8))
    x = np.arange(len(order))
    w = 0.2
    for i, cl in enumerate(ORDER):
        vals = [
            perf[(perf.cell_line == cl) & (perf.estimator == e)]["energy_explained"].iloc[0]
            for e in order
        ]
        ax.bar(x + (i - 1.5) * w, vals, width=w, label=cl, color=CCOL[cl])
    ax.axhline(0, color="k", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [e.replace("_EVAL_ONLY", "\n(oracle)").replace("_", " ") for e in order], fontsize=8
    )
    ax.set_ylabel("response energy explained")
    ax.set_ylim(-1.2, 0.6)
    ax.set_title(
        "Scale calibration — not template offset — is what fixes negative energy "
        "explained\n(direct_basal is off-scale at -6 to -17)"
    )
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig1_template_scale.png", dpi=140)
    plt.close(fig)

    # --- 2. basal vs alpha diagnostics -----------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    axes[0].bar(diag.cell_line, diag.r_alpha_basaldev, color=[CCOL[c] for c in diag.cell_line])
    axes[0].axhline(0, color="k", lw=1)
    axes[0].set_ylabel("r(alpha, basal deviation)")
    axes[0].set_title("Gene-wise alignment is ~0 and inconsistent in sign")
    axes[1].bar(
        diag.cell_line,
        diag.alpha_frac_in_source_basal_span,
        color=[CCOL[c] for c in diag.cell_line],
    )
    axes[1].set_ylabel("fraction of ||alpha||^2 in the source basal span")
    axes[1].set_ylim(0, 1)
    axes[1].set_title("alpha lies almost entirely OUTSIDE\nthe basal-deviation subspace")
    axes[2].bar(diag.cell_line, diag.global_scalar_k, color=[CCOL[c] for c in diag.cell_line])
    axes[2].axhline(0, color="k", lw=1)
    axes[2].set_ylabel("fitted global scalar k")
    axes[2].set_title("Sign of the global scalar is not consistent")
    fig.suptitle("A. Basal control expression does not encode the context response template")
    fig.tight_layout()
    fig.savefig(OUT / "fig2_basal_alpha.png", dpi=140)
    plt.close(fig)

    # --- 3. source agreement ---------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.4))
    norm = bsum[bsum.target == "success_norm"]
    raw = bsum[bsum.target == "success_raw"]
    xx = np.arange(len(ORDER))
    for off, frame, lab, col in (
        (-0.2, raw, "vs raw success", "#4C78A8"),
        (0.2, norm, "vs reliability-normalised", "#6F5CC4"),
    ):
        vals = [frame[frame.cell_line == cl]["spearman"].iloc[0] for cl in ORDER]
        los = [frame[frame.cell_line == cl]["sp_lo"].iloc[0] for cl in ORDER]
        his = [frame[frame.cell_line == cl]["sp_hi"].iloc[0] for cl in ORDER]
        axes[0].bar(
            xx + off,
            vals,
            width=0.38,
            label=lab,
            color=col,
            yerr=[np.array(vals) - np.array(los), np.array(his) - np.array(vals)],
            capsize=4,
        )
    partial = bsum[bsum.target == "partial_all_confounds"]
    pv = [partial[partial.cell_line == cl]["spearman"].iloc[0] for cl in ORDER]
    axes[0].plot(xx, pv, "k_", markersize=26, label="partial (all confounds)")
    axes[0].set_xticks(xx)
    axes[0].set_xticklabels(ORDER)
    axes[0].set_ylabel("Spearman")
    axes[0].set_title("B. Source agreement predicts transfer success in every fold")
    axes[0].legend(fontsize=8)
    axes[0].axhline(0, color="k", lw=0.8)

    for cl in ORDER:
        g = feats[feats.cell_line == cl]
        axes[1].scatter(
            g.source_agreement,
            g.success_norm,
            s=4,
            alpha=0.2,
            color=CCOL[cl],
            edgecolors="none",
            label=cl,
        )
    axes[1].set_xlabel("source agreement (sources only)")
    axes[1].set_ylabel("reliability-normalised transfer success")
    axes[1].legend(fontsize=8, markerscale=2)
    axes[1].set_title("Pooled relationship")
    fig.tight_layout()
    fig.savefig(OUT / "fig3_source_agreement.png", dpi=140)
    plt.close(fig)

    # --- 4. pathway vs gene-level gamma ----------------------------------
    gene = (
        gene_gamma[gene_gamma.baseline.isin(["nearest_basal", "basal_affine"])]
        .groupby(["cell_line", "baseline"])["r_gamma"]
        .median()
        .reset_index()
    )
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
    for ax, coll in zip(axes, ["hallmark", "reactome"], strict=True):
        pc = path[path.collection == coll]
        xx = np.arange(len(ORDER))
        gvals = [
            gene[(gene.cell_line == cl) & (gene.baseline == "basal_affine")]["r_gamma"].iloc[0]
            for cl in ORDER
        ]
        pvals = [
            pc[(pc.cell_line == cl) & (pc.baseline == "basal_affine")]["r_gamma"].iloc[0]
            for cl in ORDER
        ]
        ceil = [pc[pc.cell_line == cl]["gamma_ceiling_full"].iloc[0] for cl in ORDER]
        ax.bar(xx - 0.2, gvals, width=0.38, label="gene level", color="#C9CED6")
        ax.bar(xx + 0.2, pvals, width=0.38, label=f"{coll} pathway level", color="#D26A3A")
        ax.plot(xx + 0.2, ceil, "k_", markersize=22, label="pathway ceiling")
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xticks(xx)
        ax.set_xticklabels(ORDER)
        ax.set_ylabel(r"median $r(\gamma_{true}, \hat\gamma)$")
        ax.set_title(f"C. {coll.title()} ({int(pc.n_pathways.iloc[0])} sets)")
        ax.legend(fontsize=8)
    fig.suptitle("Pathway-resolution gamma is far more recoverable than gene-level gamma")
    fig.tight_layout()
    fig.savefig(OUT / "fig4_pathway_gamma.png", dpi=140)
    plt.close(fig)

    # --- 5. candidate D ---------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.4))
    for cl in ORDER:
        v = dcand[dcand.cell_line == cl]["D_unexplained_fraction"].dropna()
        axes[0].hist(
            np.clip(v, -0.5, 3.0), bins=70, histtype="step", lw=1.6, label=cl, color=CCOL[cl]
        )
    axes[0].axvline(1.0, color="k", ls="--", lw=1.2)
    axes[0].annotate(
        "worse than predicting zero →", xy=(1.05, axes[0].get_ylim()[1] * 0.85), fontsize=8
    )
    axes[0].set_xlabel("D_unexplained_fraction")
    axes[0].set_ylabel("perturbations")
    axes[0].set_title("Candidate D: reliable unexplained fraction")
    axes[0].legend(fontsize=8)

    for cl in ORDER:
        g = dcand[dcand.cell_line == cl]
        axes[1].scatter(
            g.reliable_energy,
            np.clip(g.D_unexplained_fraction, -0.5, 3.0),
            s=4,
            alpha=0.2,
            color=CCOL[cl],
            edgecolors="none",
            label=cl,
        )
    axes[1].axhline(1.0, color="k", ls="--", lw=1)
    axes[1].set_xscale("symlog")
    axes[1].set_xlabel("reliable target energy <h1,h2>")
    axes[1].set_ylabel("D_unexplained_fraction")
    axes[1].set_title("The per-pair ratio destabilises as reliable energy → 0")
    axes[1].legend(fontsize=8, markerscale=2)
    fig.tight_layout()
    fig.savefig(OUT / "fig5_candidate_d.png", dpi=140)
    plt.close(fig)
    print("wrote 5 figures to", OUT)


if __name__ == "__main__":
    main()
