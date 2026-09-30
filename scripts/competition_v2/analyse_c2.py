"""C2 decisions from the public fold scores, exactly as predeclared.

* ``--stage generator``: null qualification, G3 justification, rule J at a = 1, the
  selected generator (``generator_decision.json``).
* ``--stage amplitude``: the global amplitude grid under G0 and the selected generator,
  rule J per grid value, ``a*``, the nested estimate, and whether per-target amplitude
  may run (``amplitude_decision.json``).
* ``--stage final``: the C2 verdict over every scored candidate (``c2_decision.json``)
  plus the report tables and figures.

Rulers: primary = the frozen C1 fold anchors (G0-emitted mean response, split half).
Sensitivity = the same mean response emitted through G1c, or G1ci if that is selected
(amendment 1). MSE is clamped at 0 inside Overall, and unclamped for the member test.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from virtual_cell.competition_v2 import evaluation  # noqa: E402
from virtual_cell.visualization import style  # noqa: E402

OUT = ROOT / "outputs" / "competition_v2" / "c2_calibration"
FIG = ROOT / "reports" / "competition_v2" / "figures"
FOLDS = ["H1", "K562"]
MEMBERS = evaluation.MEMBERS
SHORT = {
    "pds_cosine": "PDS",
    "expr_mse_unbiased_capped_norm": "MSE",
    "de_wilcoxon_lfc_nmae": "NMAE",
    "de_wilcoxon_direction_fidelity_yield_raw": "FID",
    "de_wilcoxon_direction_reach_raw": "REACH",
    "de_wilcoxon_sig_jaccard": "JAC",
}
NON_PDS = ["MSE", "NMAE", "FID", "REACH", "JAC"]
QUAL_FRACTION = 0.01
C1_ARM = "G0_a1.00"


def load_raw(phases=(1, 2, 3)) -> dict[str, pd.DataFrame]:
    out = {}
    for f in FOLDS:
        frames = []
        for k in phases:
            p = OUT / "folds" / f / f"phase{k}_scores_raw.csv"
            if p.exists():
                frames.append(pd.read_csv(p, index_col=0))
        raw = pd.concat(frames)
        out[f] = raw[~raw.index.duplicated(keep="first")]
    return out


def scale(raw: pd.DataFrame, anchor: str) -> pd.DataFrame:
    """Scaled members on a fixed ruler; ``MSE`` unclamped, ``Overall`` with MSE clamped."""
    b = raw.loc[anchor, MEMBERS].astype(float)
    r = raw.loc["ANCHOR_split_half", MEMBERS].astype(float)
    s = (raw[MEMBERS].astype(float) - b) / (r - b)
    s.columns = [SHORT[m] for m in MEMBERS]
    clamped = s.copy()
    clamped["MSE"] = clamped["MSE"].clip(0, 1)
    s["Overall"] = clamped[[SHORT[m] for m in MEMBERS]].mean(axis=1)
    return s


def scaled_all(raw: dict, anchor: str) -> dict[str, pd.DataFrame]:
    return {f: scale(raw[f], anchor) for f in FOLDS}


def fold_mean(scaled: dict, arm: str) -> pd.Series:
    return pd.concat([scaled[f].loc[arm] for f in FOLDS], axis=1).mean(axis=1)


def null_table() -> pd.DataFrame:
    frames = []
    for f in FOLDS:
        t = pd.read_csv(OUT / "folds" / f / "phase1_null.csv")
        t.insert(0, "fold", f)
        frames.append(t)
    return pd.concat(frames, ignore_index=True)


def qualification(nulls: pd.DataFrame) -> dict[str, bool]:
    q = {}
    for gen, grp in nulls.groupby("generator"):
        limit = np.floor(QUAL_FRACTION * grp.tested_genes)
        q[gen] = bool(len(grp) == len(FOLDS) and (grp.null_fp_de_median <= limit).all())
    return q


def rule_j(arm, gen, primary, sens, qualified) -> dict:
    c1 = fold_mean(primary, C1_ARM)
    cand = fold_mean(primary, arm)
    better = [m for m in NON_PDS if cand[m] > c1[m]]
    pds_ret = {f: float(primary[f].loc[arm, "PDS"] / primary[f].loc[C1_ARM, "PDS"]) for f in FOLDS}
    crit = {
        "1_overall_up": bool(cand["Overall"] > c1["Overall"]),
        "2_two_members_up": len(better) >= 2,
        "3_pds_95pct_each_fold": all(v >= 0.95 for v in pds_ret.values()),
        "4_null_qualified": bool(qualified.get(gen, False)),
        "5_sensitivity_overall_up": bool(
            fold_mean(sens, arm)["Overall"] > fold_mean(sens, C1_ARM)["Overall"]
        ),
    }
    return {
        "arm": arm,
        "generator": gen,
        "pass": all(crit.values()),
        "criteria": crit,
        "members_up": better,
        "pds_retention": pds_ret,
        "mean_scaled": cand.round(4).to_dict(),
        "delta_vs_c1": (cand - c1).round(4).to_dict(),
    }


def generator_of(arm: str) -> str:
    return arm.split("_")[0]


def sens_anchor(selected: str) -> str:
    return "ANCHOR_mean_response_G1ci" if selected == "G1ci" else "ANCHOR_mean_response_G1c"


# --------------------------------------------------------------------------- stages
def stage_generator() -> dict:
    raw = load_raw((1,))
    nulls = null_table()
    q = qualification(nulls)
    g2 = nulls[nulls.generator == "G2"].set_index("fold")
    struct = {
        f: pd.read_csv(OUT / "folds" / f / "phase1_scores_raw.csv", index_col=0) for f in FOLDS
    }
    g2_var = {f: float(struct[f].loc["G2_null", "cells_var_ratio"]) for f in FOLDS}
    g3_justified = all(v < 0.8 for v in g2_var.values()) and all(
        "G3_a1.00" in raw[f].index for f in FOLDS
    )
    primary = scaled_all(raw, "ANCHOR_mean_response")
    cands = ["G1c", "G1ci", "G1b", "G2"] + (["G3"] if g3_justified else [])
    results = {}
    for gen in cands:
        sens = scaled_all(raw, sens_anchor(gen))
        results[gen] = rule_j(f"{gen}_a1.00", gen, primary, sens, q)
    qualified = [g for g in cands if q.get(g)]
    passing = [g for g in qualified if results[g]["pass"]]
    pool = passing or qualified
    if pool:
        selected = max(pool, key=lambda g: results[g]["mean_scaled"]["Overall"])
    else:
        selected = "G0"
    decision = {
        "null_qualification": q,
        "null_limit_fraction_of_tested": QUAL_FRACTION,
        "g2_null_variance_ratio": g2_var,
        "g3_justified": g3_justified,
        "rule_j_at_a1": results,
        "qualified": qualified,
        "passing_rule_j_at_a1": passing,
        "selected_generator": selected,
        "selection_basis": "passes rule J"
        if passing
        else "qualified, highest Overall (amplitude study only)",
        "c1_mean_scaled": fold_mean(primary, C1_ARM).round(4).to_dict(),
        "g2_null_by_fold": g2.null_fp_de_median.to_dict(),
    }
    (OUT / "generator_decision.json").write_text(json.dumps(decision, indent=2))
    print(
        json.dumps(
            {
                k: decision[k]
                for k in ("null_qualification", "passing_rule_j_at_a1", "selected_generator")
            },
            indent=1,
        )
    )
    return decision


def stage_amplitude() -> dict:
    gdec = json.loads((OUT / "generator_decision.json").read_text())
    best = gdec["selected_generator"]
    raw = load_raw((1, 2))
    q = gdec["null_qualification"] | {"G0": qualification(null_table()).get("G0", False)}
    primary = scaled_all(raw, "ANCHOR_mean_response")
    sens = scaled_all(raw, sens_anchor(best))
    grid = {}
    for gen in sorted({"G0", best}):
        for a in [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]:
            arm = f"{gen}_a{a:.2f}"
            if all(arm in raw[f].index for f in FOLDS):
                grid[arm] = rule_j(arm, gen, primary, sens, q)
    best_arms = {a: r for a, r in grid.items() if generator_of(a) == best}
    passing = [a for a, r in best_arms.items() if r["pass"]]
    if passing:
        a_star_arm = max(passing, key=lambda a: best_arms[a]["mean_scaled"]["Overall"])
    else:
        a_star_arm = max(best_arms, key=lambda a: best_arms[a]["mean_scaled"]["Overall"])
    a_star = float(a_star_arm.split("_a")[1])
    one = f"{best}_a1.00"
    consistent = a_star != 1.0 and all(
        primary[f].loc[a_star_arm, "Overall"] > primary[f].loc[one, "Overall"] for f in FOLDS
    )
    nested = {}
    for sel, ev in [("H1", "K562"), ("K562", "H1")]:
        pick = max(best_arms, key=lambda a: primary[sel].loc[a, "Overall"])
        nested[f"select_{sel}_score_{ev}"] = {
            "picked": pick,
            "overall": float(primary[ev].loc[pick, "Overall"]),
            "c1_overall": float(primary[ev].loc[C1_ARM, "Overall"]),
        }
    decision = {
        "generator": best,
        "grid": grid,
        "passing": passing,
        "a_star_arm": a_star_arm,
        "a_star": a_star,
        "a_star_passes_rule_j": a_star_arm in passing,
        "run_per_target": bool(consistent),
        "consistency_rule": "a* != 1 and a* beats a = 1 on primary Overall in both folds",
        "nested": nested,
    }
    (OUT / "amplitude_decision.json").write_text(json.dumps(decision, indent=2))
    print(
        json.dumps(
            {
                k: decision[k]
                for k in ("generator", "passing", "a_star_arm", "run_per_target", "nested")
            },
            indent=1,
        )
    )
    return decision


def stage_final() -> dict:
    gdec = json.loads((OUT / "generator_decision.json").read_text())
    adec = json.loads((OUT / "amplitude_decision.json").read_text())
    best = gdec["selected_generator"]
    raw = load_raw()
    q = gdec["null_qualification"] | {"G0": qualification(null_table()).get("G0", False)}
    primary = scaled_all(raw, "ANCHOR_mean_response")
    sens = scaled_all(raw, sens_anchor(best))
    arms = [
        a
        for a in raw["H1"].index
        if a in raw["K562"].index
        and not a.startswith(("ANCHOR", "IDEAL"))
        and not a.endswith("_null")
        and a != C1_ARM
    ]
    verdicts = {a: rule_j(a, generator_of(a), primary, sens, q) for a in arms}
    passing = [a for a, v in verdicts.items() if v["pass"]]
    chosen = max(passing, key=lambda a: verdicts[a]["mean_scaled"]["Overall"]) if passing else None
    decision = {
        "passing": passing,
        "c2_candidate": chosen,
        "c2_pass": chosen is not None,
        "verdicts": verdicts,
        "generator_decision": gdec["selected_generator"],
        "a_star": adec["a_star"],
    }
    (OUT / "c2_decision.json").write_text(json.dumps(decision, indent=2))
    tables(raw, primary, sens)
    figures(raw, primary, best)
    print(json.dumps({"passing": passing, "c2_candidate": chosen}, indent=1))
    return decision


# --------------------------------------------------------------------------- outputs
def tables(raw, primary, sens) -> None:
    rows = []
    for f in FOLDS:
        for arm in raw[f].index:
            row = {"fold": f, "arm": arm}
            row.update({k: v for k, v in primary[f].loc[arm].items()})
            row["Overall_sensitivity"] = sens[f].loc[arm, "Overall"]
            for col in [
                "predicted_de_median",
                "real_de_median",
                "directional_precision_pooled",
                "directional_recall_median",
                "de_overlap_median",
                "effect_norm_ratio_median",
                "library_bias",
                "genes_detected_bias",
                "pseudobulk_dispersion_ratio",
                "cells_var_ratio",
                "cells_fano_ratio",
                "cells_logcpm_var_ratio",
                "cells_zero_fraction",
                "cells_real_zero_fraction",
                "cells_cosine_distance",
                "cells_real_cosine_distance",
            ]:
                row[col] = raw[f].loc[arm].get(col, np.nan)
            rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "scores_scaled_all.csv", index=False)
    means = table.groupby("arm")[["PDS", *NON_PDS, "Overall", "Overall_sensitivity"]].mean()
    means.to_csv(OUT / "scores_scaled_fold_mean.csv")
    null_table().to_csv(OUT / "null_qualification.csv", index=False)


def figures(raw, primary, best) -> None:
    style.apply()
    (FIG / "sources").mkdir(parents=True, exist_ok=True)
    # B: mean-response ceiling vs generated cells, per member, C1 (a = 1)
    arms = ["IDEAL_a1.00", "G0_a1.00", "G1c_a1.00", "G1ci_a1.00", "G1b_a1.00", "G2_a1.00"]
    arms = [a for a in arms if all(a in primary[f].index for f in FOLDS)]
    cols = ["PDS", *NON_PDS, "Overall"]
    src = pd.DataFrame({a: fold_mean(primary, a)[cols] for a in arms}).T
    src.to_csv(FIG / "sources" / "c2_B_mean_vs_generator.csv")
    fig, ax = plt.subplots(figsize=(style.FULL_WIDTH, 3.4))
    w = 0.8 / len(arms)
    palette = [
        style.BLACK,
        style.GREY,
        style.BLUE,
        style.SKY,
        style.GREEN,
        style.ORANGE,
        style.PURPLE,
    ]
    for i, a in enumerate(arms):
        ax.bar(np.arange(len(cols)) + i * w, src.loc[a], w, label=a, color=palette[i])
    style.reference_line(ax, 0)
    ax.set_xticks(np.arange(len(cols)) + 0.4 - w / 2, cols)
    ax.set_ylabel("local scaled (primary ruler), mean H1/K562")
    ax.set_title("C1 mean prediction: idealised pseudobulk vs generated cells", fontsize=9)
    ax.legend(fontsize=6, ncol=3)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(FIG / f"c2_B_mean_vs_generator.{ext}", dpi=200)
    plt.close(fig)
    # C: amplitude tradeoff
    rows = []
    for gen in sorted({"G0", best, "IDEAL"}):
        for a in [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]:
            arm = f"{gen}_a{a:.2f}"
            if all(arm in primary[f].index for f in FOLDS):
                rows.append({"generator": gen, "a": a, **fold_mean(primary, arm)[cols].to_dict()})
    amp = pd.DataFrame(rows)
    amp.to_csv(FIG / "sources" / "c2_C_amplitude.csv", index=False)
    fig, axes = plt.subplots(1, len(cols), figsize=(style.FULL_WIDTH, 2.4), sharex=True)
    for ax, c in zip(axes, cols, strict=True):
        for gen, grp in amp.groupby("generator"):
            ax.plot(grp.a, grp[c], marker="o", ms=3, label=gen)
        ax.set_title(c, fontsize=8)
        ax.tick_params(labelsize=6)
    axes[0].legend(fontsize=6)
    fig.supxlabel("global amplitude a", fontsize=8)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(FIG / f"c2_C_amplitude.{ext}", dpi=200)
    plt.close(fig)
    # D: null false-positive DE per generator
    nulls = null_table()
    nulls.to_csv(FIG / "sources" / "c2_D_null_de.csv", index=False)
    fig, ax = plt.subplots(figsize=(style.HALF_WIDTH * 1.4, 2.8))
    piv = nulls.pivot(index="generator", columns="fold", values="null_fp_de_median")
    piv.plot.bar(ax=ax, color=[style.BLUE, style.ORANGE])
    ax.set_yscale("symlog", linthresh=10)
    ax.set_ylabel("median false-positive DE genes\n(zero effect vs reference)")
    ax.set_title("Null DE qualification", fontsize=9)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(FIG / f"c2_D_null_de.{ext}", dpi=200)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["generator", "amplitude", "final"], required=True)
    args = parser.parse_args()
    {"generator": stage_generator, "amplitude": stage_amplitude, "final": stage_final}[args.stage]()


if __name__ == "__main__":
    main()
