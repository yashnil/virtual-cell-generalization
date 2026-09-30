"""C3 decision (rule O over the capacity ladder), tables and figures R1-R5.

    uv run python scripts/competition_v2/analyse_c3.py

Reads ``outputs/competition_v2/c3_fusion/`` (mean stage + VCC folds) and writes
``c3_decision.json``, ``vcc_scaled.csv`` and ``reports/competition_v2/figures/c3_*``.
The ruler is the frozen local one (G0-emitted mean response, split half), as in C1/C2.
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
from scipy import stats as sstats  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyse_c2 as a2  # noqa: E402

from virtual_cell.visualization import style  # noqa: E402

OUT = ROOT / "outputs" / "competition_v2" / "c3_fusion"
FIG = ROOT / "reports" / "competition_v2" / "figures"
FOLDS = ["H1", "K562"]
LADDER = {2: ["R1", "R2", "Q"], 3: ["S1", "S2"], 5: ["P", "C3d_m", "C3d_s"], 6: ["K"]}
NON_PDS = ["MSE", "NMAE", "FID", "REACH"]
MATERIAL = 0.02


def vcc_scaled() -> dict[str, pd.DataFrame]:
    out = {}
    for f in FOLDS:
        raw = pd.read_csv(OUT / "folds" / f / "scores_raw.csv", index_col=0)
        out[f] = a2.scale(raw, "ANCHOR_mean_response")
    return out


def rule_o(arm: str, vcc: dict, mean: pd.DataFrame) -> dict:
    c1 = pd.concat([vcc[f].loc["C1a"] for f in FOLDS], axis=1).mean(axis=1)
    cand = pd.concat([vcc[f].loc[arm] for f in FOLDS], axis=1).mean(axis=1)
    gains = {f: float(vcc[f].loc[arm, "Overall"] - vcc[f].loc["C1a", "Overall"]) for f in FOLDS}
    total = sum(gains.values())
    m = mean.set_index(["heldout", "candidate"])
    fold_mean = lambda col, c: float(np.mean([m.loc[(h, c), col] for h in ["H1", "K562", "CD4"]]))  # noqa: E731
    pds = {f: float(vcc[f].loc[arm, "PDS"] / vcc[f].loc["C1a", "PDS"]) for f in FOLDS}
    pds["CD4_effect"] = float(
        m.loc[("CD4", arm), "effect_pds"] / m.loc[("CD4", "C1a"), "effect_pds"]
    )
    deltas = {k: float(cand[k] - c1[k]) for k in NON_PDS}
    crit = {
        "1_overall_up": bool(cand["Overall"] > c1["Overall"]),
        "2_pds_95pct_all": all(v >= 0.95 for v in pds.values()),
        "3_cosine_up": fold_mean("cosine", arm) > fold_mean("cosine", "C1a"),
        "4_sign_acc_up": fold_mean("sign_acc", arm) > fold_mean("sign_acc", "C1a"),
        "5_both_folds_up": all(v > 0 for v in gains.values()),
        "6_no_fold_over_75pct": bool(total > 0 and all(v / total <= 0.75 for v in gains.values())),
        "7_license_clean": True,
        "8_member_up_no_material_loss": any(v > 0 for v in deltas.values())
        and all(v >= -MATERIAL for v in deltas.values()),
    }
    return {
        "arm": arm,
        "pass": all(crit.values()),
        "criteria": crit,
        "overall": float(cand["Overall"]),
        "overall_c1": float(c1["Overall"]),
        "fold_gain": gains,
        "pds_retention": pds,
        "member_delta": deltas,
        "cosine_delta": fold_mean("cosine", arm) - fold_mean("cosine", "C1a"),
        "sign_acc_delta": fold_mean("sign_acc", arm) - fold_mean("sign_acc", "C1a"),
    }


def main() -> None:
    vcc = vcc_scaled()
    mean = pd.read_csv(OUT / "candidates_mean_level.csv")
    trig = json.loads((OUT / "triggers.json").read_text())
    arms = [a for a in vcc["H1"].index if a in vcc["K562"].index and not a.startswith("ANCHOR")]
    verdicts = {a: rule_o(a, vcc, mean) for a in arms if a != "C1a"}
    selected, rung_used = None, None
    for rung, cands in LADDER.items():
        passing = [c for c in cands if c in verdicts and verdicts[c]["pass"]]
        if passing:
            selected = max(passing, key=lambda c: verdicts[c]["overall"])
            rung_used = rung
            break
    rows = []
    for f in FOLDS:
        for a in arms:
            rows.append({"fold": f, "arm": a, **vcc[f].loc[a].to_dict()})
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "vcc_scaled.csv", index=False)
    decision = {
        "c3_pass": selected is not None,
        "selected": selected,
        "rung": rung_used,
        "verdicts": verdicts,
        "triggers": trig,
        "infeasible": {"C3b": "no control profile for the CD4 DE source"},
    }
    (OUT / "c3_decision.json").write_text(json.dumps(decision, indent=2, default=float))
    figures(vcc, mean, arms)
    print(json.dumps({"c3_pass": decision["c3_pass"], "selected": selected}, indent=1))
    for a, v in verdicts.items():
        fails = [k for k, ok in v["criteria"].items() if not ok]
        print(f"{a:6s} overall {v['overall']:+.4f} (C1 {v['overall_c1']:+.4f}) fails {fails}")


def figures(vcc, mean, arms) -> None:
    style.apply()
    src = FIG / "sources"
    src.mkdir(parents=True, exist_ok=True)

    def save(fig, name):
        fig.tight_layout()
        for ext in ("png", "svg"):
            fig.savefig(FIG / f"{name}.{ext}", dpi=200)
        plt.close(fig)

    # R1 donor-source transfer matrix
    tm = pd.read_csv(OUT / "source_transfer_matrix.csv")
    tm.to_csv(src / "c3_R1_transfer_matrix.csv", index=False)
    fig, axes = plt.subplots(1, 3, figsize=(style.FULL_WIDTH, 2.6))
    for ax, (col, title) in zip(
        axes,
        [
            ("cosine", "response cosine"),
            ("sign_acc", "top-200 sign accuracy"),
            ("norm_ratio_median", "norm ratio ‖p‖/‖t‖"),
        ],
        strict=True,
    ):
        piv = tm.pivot(index="source", columns="heldout", values=col).reindex(
            index=["H1", "K562", "CD4"], columns=["H1", "K562", "CD4"]
        )
        ax.imshow(piv.to_numpy(dtype=float), cmap="viridis")
        for i in range(3):
            for j in range(3):
                v = piv.iloc[i, j]
                if np.isfinite(v):
                    ax.text(j, i, f"{v:.3f}", ha="center", va="center", fontsize=7, color="white")
        ax.set_xticks(range(3), piv.columns, fontsize=7)
        ax.set_yticks(range(3), piv.index, fontsize=7)
        ax.set_xlabel("held-out atlas", fontsize=7)
        ax.set_title(title, fontsize=8)
    axes[0].set_ylabel("donor source", fontsize=7)
    save(fig, "c3_R1_transfer_matrix")

    # R2 predicted vs true magnitude
    rows = []
    fig, axes = plt.subplots(1, 3, figsize=(style.FULL_WIDTH, 2.6))
    for ax, h in zip(axes, ["H1", "K562", "CD4"], strict=True):
        d = np.load(OUT / f"diagnostic_{h}.npz", allow_pickle=True)  # our own file
        t, p = d["truth"].astype(float), d["fused_c1a"].astype(float)
        excl = np.isin(d["genes"], d["targets"])
        v = d["truth_mask"] & ~excl[None, :] & (p != 0)
        pn = np.sqrt((np.where(v, p, 0) ** 2).sum(1))
        tn = np.sqrt((np.where(v, t, 0) ** 2).sum(1))
        sig = tn * np.sqrt(1 - np.nan_to_num(d["noise_fraction"]))
        ok = (pn > 0) & (tn > 0)
        rows += [
            {"heldout": h, "pred_norm": a, "true_norm": b, "true_signal_norm": c}
            for a, b, c in zip(pn[ok], tn[ok], sig[ok], strict=True)
        ]
        ax.scatter(tn[ok], pn[ok], s=5, alpha=0.5, label="vs measured truth")
        ax.scatter(sig[ok], pn[ok], s=5, alpha=0.5, label="vs noise-corrected")
        lim = max(pn[ok].max(), tn[ok].max())
        ax.plot([0, lim], [0, lim], color=style.GREY, lw=0.8)
        ax.set_title(f"held out {h}", fontsize=8)
        ax.set_xlabel("true response norm", fontsize=7)
    axes[0].set_ylabel("C1a predicted norm", fontsize=7)
    axes[0].legend(fontsize=6)
    pd.DataFrame(rows).to_csv(src / "c3_R2_magnitude.csv", index=False)
    save(fig, "c3_R2_magnitude")

    # R3 sign accuracy by consensus
    sc = pd.read_csv(OUT / "sign_consensus.csv")
    sc.to_csv(src / "c3_R3_sign_consensus.csv", index=False)
    fig, ax = plt.subplots(figsize=(style.HALF_WIDTH * 1.4, 2.6))
    piv = sc.pivot(index="heldout", columns="class", values="sign_acc")[
        ["agree", "single", "conflict"]
    ]
    piv.plot.bar(ax=ax, color=[style.GREEN, style.GREY, style.VERMILION])
    ax.axhline(0.5, color="black", lw=0.8, ls="--")
    ax.set_ylim(0.45, 0.72)
    ax.set_ylabel("top-200 sign accuracy")
    ax.set_title("Sign accuracy by source consensus", fontsize=9)
    save(fig, "c3_R3_sign_consensus")

    # R4 C1 vs C3 VCC metrics
    cols = ["PDS", "MSE", "NMAE", "FID", "REACH", "JAC", "Overall"]
    show = [a for a in ["C1a", "R1", "Q", "S1", "S2", "P", "C3d_m", "C3d_s"] if a in arms]
    rows = [{"fold": f, "arm": a, **vcc[f].loc[a, cols].to_dict()} for f in FOLDS for a in show]
    pd.DataFrame(rows).to_csv(src / "c3_R4_vcc.csv", index=False)
    fig, axes = plt.subplots(1, 2, figsize=(style.FULL_WIDTH, 3.0), sharey=True)
    w = 0.8 / len(show)
    for ax, f in zip(axes, FOLDS, strict=True):
        for i, a in enumerate(show):
            ax.bar(np.arange(len(cols)) + i * w, vcc[f].loc[a, cols], w, label=a)
        ax.set_xticks(np.arange(len(cols)) + 0.4 - w / 2, cols, fontsize=7)
        style.reference_line(ax, 0)
        ax.set_title(f"held out {f}", fontsize=8)
    axes[0].set_ylabel("local scaled (C1 ruler)")
    axes[1].legend(fontsize=6, ncol=2)
    save(fig, "c3_R4_vcc")

    # R5 source weights selected across folds
    rows = []
    for h in ["H1", "K562", "CD4"]:
        spec = json.loads((OUT / f"specs_{h}.json").read_text())
        for cand, wts in spec["weights"].items():
            for s, v in (wts or {}).items():
                rows.append({"heldout": h, "candidate": cand, "source": s, "weight": v})
        for cand, sc_ in spec["scales"].items():
            for s, v in (sc_ or {}).items():
                rows.append({"heldout": h, "candidate": cand + "_scale", "source": s, "weight": v})
    wt = pd.DataFrame(rows)
    wt.to_csv(src / "c3_R5_weights.csv", index=False)
    fig, axes = plt.subplots(1, 3, figsize=(style.FULL_WIDTH, 2.6), sharey=True)
    for ax, h in zip(axes, ["H1", "K562", "CD4"], strict=True):
        piv = wt[wt.heldout == h].pivot(index="candidate", columns="source", values="weight")
        piv.plot.bar(ax=ax, legend=h == "H1")
        ax.set_title(f"held out {h}", fontsize=8)
        ax.tick_params(labelsize=6)
    axes[0].set_ylabel("normalised weight / scale a_s")
    save(fig, "c3_R5_weights")

    # perturbation-dependence of donor quality
    pt = pd.read_csv(OUT / "per_target_source_cosine.csv")
    out = {}
    for h, g in pt.groupby("heldout"):
        ok = g.rel_diff.notna()
        out[h] = {
            "n_targets": int(len(g)),
            "fraction_first_source_better": float((g.cos_diff > 0).mean()),
            "cos_diff_sd": float(g.cos_diff.std()),
            "spearman_cos_diff_vs_rel_diff": float(
                sstats.spearmanr(g.cos_diff[ok], g.rel_diff[ok]).statistic
            )
            if ok.sum() > 10
            else None,
            "n_with_rel": int(ok.sum()),
        }
    (OUT / "perturbation_dependence.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
