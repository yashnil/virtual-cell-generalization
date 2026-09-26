"""Competition-v2 C1 figure set (A-E): PNG + SVG, each with its source table.

A. 300-target direct-evidence coverage, V1 vs C0 vs C1
B. public score comparison (local anchors), C0 / C0 without X-Atlas / C1a / C1b
C. source agreement vs held-out transfer quality (H1, K562, CD4)
D. perturbation-signature diversity on the Arc bundles, V1 vs C0 vs C1
E. X-Atlas ablation: with vs without, same backbone, same folds

Outputs: ``reports/competition_v2/figures/c1_*.{png,svg}`` and
``reports/competition_v2/figures/sources/c1_*.csv``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from virtual_cell.visualization import style  # noqa: E402

BASE = ROOT / "outputs" / "competition_v2" / "c1_license_clean"
FIG = ROOT / "reports" / "competition_v2" / "figures"
SRC = FIG / "sources"
COLORS = {
    "V1": style.GREY,
    "C0": style.VERMILION,
    "C0_withX": style.VERMILION,
    "C0_noX": style.ORANGE,
    "C1a": style.BLUE,
    "C1b": style.SKY,
    "C1": style.BLUE,
}


def save(fig, name: str, table: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    SRC.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "svg"):
        fig.savefig(
            FIG / f"{name}.{ext}",
            dpi=300 if ext == "png" else None,
            metadata={"Date": None} if ext == "svg" else None,
        )
    table.to_csv(SRC / f"{name}.csv", index=False)
    plt.close(fig)
    print(f"  wrote {name}")


def fig_a() -> None:
    cov = json.loads((BASE / "coverage_summary.json").read_text())
    bins = ["0", "1", "2", "3+"]
    table = pd.DataFrame(
        [
            {"model": m, "contexts": b, "targets": cov[k][b]}
            for m, k in [
                ("V1", "v1"),
                ("C0 (with X-Atlas)", "c0_expanded"),
                ("C1 (license-clean)", "c1_license_clean"),
            ]
            for b in bins
        ]
    )
    fig, ax = plt.subplots(figsize=(style.HALF_WIDTH * 1.5, 2.8))
    x = np.arange(len(bins))
    for i, (m, c) in enumerate(
        zip(table.model.unique(), [COLORS["V1"], COLORS["C0"], COLORS["C1"]], strict=True)
    ):
        v = table[table.model == m].targets.to_numpy()
        ax.bar(x + (i - 1) * 0.27, v, 0.27, color=c, label=m, edgecolor="black", lw=0.4)
        for xi, vi in zip(x, v, strict=True):
            ax.text(xi + (i - 1) * 0.27, vi + 3, str(vi), ha="center", fontsize=6)
    ax.set_xticks(x, [f"{b} contexts" for b in bins])
    ax.set_ylabel("Arc validation targets (of 300)")
    ax.set_title("A. Direct-evidence coverage of the 300 Arc targets", fontsize=9)
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    save(fig, "c1_A_target_coverage", table)


def fig_b() -> None:
    comp = pd.read_csv(BASE / "comparison_public.csv")
    dec = json.loads((BASE / "decision.json").read_text())
    rows = []
    for fold in ["H1", "K562"]:
        c1b_arm = dec["c1b_nested_selection"][fold]["outer_arm"]
        for label, arm in [
            ("C0_withX", "C0_withX"),
            ("C0_noX", "C0_noX"),
            ("C1a", "C1a"),
            ("C1b", c1b_arm),
        ]:
            r = comp[(comp.fold == fold) & (comp.arm == arm)].iloc[0]
            rows.append(
                {
                    "fold": fold,
                    "model": label,
                    "arm": arm,
                    "local_overall": r.local_overall,
                    "pds_raw": r.raw_pds_cosine,
                    "scaled_pds": r.scaled_pds_cosine,
                    "scaled_fid": r.scaled_de_wilcoxon_direction_fidelity_yield_raw,
                }
            )
    table = pd.DataFrame(rows)
    fig, axes = plt.subplots(1, 2, figsize=(style.FULL_WIDTH, 2.8))
    models = ["C0_withX", "C0_noX", "C1a", "C1b"]
    labels = ["C0 (X-Atlas; BLOCKED)", "C0 w/o X-Atlas", "C1a", "C1b (nested)"]
    for ax, metric, title in [
        (axes[0], "local_overall", "local Overall (0 = mean response)"),
        (axes[1], "pds_raw", "raw PDS (0.5 = chance)"),
    ]:
        x = np.arange(2)
        for i, m in enumerate(models):
            v = [
                table[(table.fold == f) & (table.model == m)][metric].iloc[0]
                for f in ["H1", "K562"]
            ]
            ax.bar(
                x + (i - 1.5) * 0.2,
                v,
                0.2,
                color=COLORS[m],
                edgecolor="black",
                lw=0.4,
                label=labels[i],
                hatch="//" if m == "C0_withX" else None,
            )
        ax.set_xticks(x, ["H1 held out", "K562 held out"])
        ax.set_title(title, fontsize=8)
        style.reference_line(ax, 0.5 if metric == "pds_raw" else 0.0)
    axes[1].set_ylim(0.45, None)
    axes[0].legend(fontsize=6, frameon=False, loc="best")
    fig.suptitle(
        "B. Public leave-one-atlas-out scores (V1 has no score on these folds; "
        "its official hidden Overall is -0.062)",
        fontsize=8,
    )
    fig.tight_layout()
    save(fig, "c1_B_public_scores", table)


def fig_c() -> None:
    frames = []
    for fold in ["H1", "K562", "CD4"]:
        pt = pd.read_csv(BASE / "folds" / fold / "per_target_agreement.csv")
        pt = pt[pt.agreement.notna() & pt.transfer_cosine_c1a.notna()]
        frames.append(
            pt.assign(fold=fold)[
                ["fold", "target", "agreement", "transfer_cosine_c1a", "n_sources_c1"]
            ]
        )
    table = pd.concat(frames)
    fig, axes = plt.subplots(1, 3, figsize=(style.FULL_WIDTH, 2.5), sharey=True)
    for ax, fold in zip(axes, ["H1", "K562", "CD4"], strict=True):
        d = table[table.fold == fold]
        ax.scatter(d.agreement, d.transfer_cosine_c1a, s=5, color=style.BLUE, alpha=0.6, lw=0)
        rho = stats.spearmanr(d.agreement, d.transfer_cosine_c1a)[0] if len(d) > 2 else np.nan
        ax.set_title(f"{fold} held out (n={len(d)}, ρ={rho:.2f})", fontsize=8)
        ax.set_xlabel("source agreement (mean pairwise cosine)")
        style.reference_line(ax, 0)
    axes[0].set_ylabel("held-out transfer cosine (C1a)")
    fig.suptitle("C. Source agreement vs held-out transfer quality", fontsize=9)
    fig.tight_layout()
    save(fig, "c1_C_agreement_vs_transfer", table)


def fig_d() -> None:
    path = BASE / "mechanistic" / "per_context.csv"
    if not path.exists():
        print("  D skipped: C1 bundle diversity not computed")
        return
    pc = pd.read_csv(path)
    pt = pd.read_csv(BASE / "mechanistic" / "per_target.csv")
    null = json.loads((BASE / "mechanistic" / "summary.json").read_text())[
        "reproducibility_null_p99"
    ]
    table = pc[
        [
            "bundle",
            "context",
            "deviation_effective_rank",
            "pairwise_cosine_effect_mean",
            "effect_norm_median",
        ]
    ]
    fig, axes = plt.subplots(1, 2, figsize=(style.FULL_WIDTH, 2.6))
    x = np.arange(3)
    for i, b in enumerate(["V1", "C0", "C1"]):
        v = pc[pc.bundle == b].deviation_effective_rank.to_numpy()
        axes[0].bar(
            x + (i - 1) * 0.27, v, 0.27, color=COLORS[b], edgecolor="black", lw=0.4, label=b
        )
    axes[0].set_xticks(x, ["A", "B", "C"])
    axes[0].set_ylabel("effective rank of target deviations")
    axes[0].legend(fontsize=7, frameon=False)
    bins = np.linspace(-0.2, 1, 49)
    for b in ["V1", "C0", "C1"]:
        axes[1].hist(pt[f"repro_{b}"], bins=bins, histtype="step", color=COLORS[b], label=b)
    axes[1].axvline(null, color="black", ls=":", lw=0.8)
    axes[1].set_xlabel("cross-context reproducibility of each target's signature")
    axes[1].set_ylabel("targets")
    fig.suptitle("D. Perturbation-signature diversity on the Arc A/B/C bundles", fontsize=9)
    fig.tight_layout()
    save(fig, "c1_D_signature_diversity", table)
    pt[["target", "n_usable_contexts_c1", "repro_V1", "repro_C0", "repro_C1"]].to_csv(
        SRC / "c1_D_signature_diversity_per_target.csv", index=False
    )


def fig_e() -> None:
    xa = pd.read_csv(BASE / "xatlas_ablation.csv")
    fig, axes = plt.subplots(1, 3, figsize=(style.FULL_WIDTH, 2.5))
    folds = list(xa.fold.unique())
    x = np.arange(len(folds))
    for ax, metric, title in [
        (axes[0], "mean_sources", "direct sources per target"),
        (axes[1], "pds_raw", "PDS (raw)"),
        (axes[2], "local_overall", "local Overall"),
    ]:
        for i, (arm, lab) in enumerate(
            [("C0_noX", "without X-Atlas"), ("C0_withX", "with X-Atlas")]
        ):
            v = [xa[(xa.fold == f) & (xa.arm == arm)][metric].iloc[0] for f in folds]
            ax.bar(
                x + (i - 0.5) * 0.35,
                v,
                0.35,
                color=COLORS[arm],
                edgecolor="black",
                lw=0.4,
                label=lab,
            )
        ax.set_xticks(x, [f.replace(" (effect-level)", "\n(effect)") for f in folds], fontsize=7)
        ax.set_title(title, fontsize=8)
    axes[1].set_ylim(0.45, None)
    axes[0].legend(fontsize=6, frameon=False)
    fig.suptitle(
        "E. X-Atlas ablation (same C0 backbone and weights, same public folds)", fontsize=9
    )
    fig.tight_layout()
    save(fig, "c1_E_xatlas_ablation", xa)


def main() -> None:
    style.apply()
    for fn in (fig_a, fig_b, fig_c, fig_d, fig_e):
        fn()


if __name__ == "__main__":
    main()
