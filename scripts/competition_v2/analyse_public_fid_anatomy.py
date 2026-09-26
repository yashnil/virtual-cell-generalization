"""Aggregate the public FID anatomy (section 16) and amplitude ladder (section 17).

Reads ``outputs/competition_v2/public_fid_anatomy/<fold>/`` and writes, next to them:

* ``anatomy_by_tier.csv`` — per fold x pure-tier arm: mean precision, mean yield,
  mean FID, median n_pred / n_real, REACH, Jaccard, NMAE;
* ``anatomy_by_magnitude.csv`` — the same for the Arc-prevalence mixed panel at a = 1,
  stratified by assigned tier and by true-effect-norm tercile;
* ``ladder_raw.csv`` / ``ladder_scaled.csv`` — six members per fold x amplitude, raw and
  scaled by the **local** anchors (b = oracle mean response, r = one split-half);
* ``fid_decomposition.json`` — the sign-vs-yield verdict inputs.

Diagnostic only; nothing is selected here.

Reproduce: ``uv run python scripts/competition_v2/analyse_public_fid_anatomy.py``
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "outputs" / "competition_v2" / "public_fid_anatomy"
FOLDS = ["replogle22k562", "replogle22rpe1", "nadig25hepg2", "nadig25jurkat"]
MEMBERS = {
    "pds_cosine": True,
    "expr_mse_unbiased_capped_norm": False,
    "de_wilcoxon_direction_fidelity_yield_raw": True,
    "de_wilcoxon_direction_reach_raw": True,
    "de_wilcoxon_sig_jaccard": True,
    "de_wilcoxon_lfc_nmae": False,
}


def summarise(frame: pd.DataFrame) -> dict:
    defined = frame.dropna(subset=["fid"])
    called = defined[defined.n_pred > 0]
    return {
        "n": len(frame),
        "fid": defined.fid.mean(),
        # precision is only defined where something was called; report both the
        # conditional mean and the share of perturbations that called nothing at all
        "precision_where_called": called.precision.mean(),
        "share_calling_nothing": float((defined.n_pred == 0).mean()),
        "yield": defined["yield"].mean(),
        "n_pred_median": frame.n_pred.median(),
        "n_real_median": frame.n_real.median(),
        "pooled_precision": called.k.sum() / max(called.n_pred.sum(), 1),
        "pooled_ratio_pred_over_real": frame.n_pred.sum() / max(frame.n_real.sum(), 1),
        "reach": frame.reach.mean(),
        "jaccard": frame.jaccard.mean(),
        "nmae": frame.nmae.mean(),
        "pred_norm_median": frame.pred_norm.median(),
        "true_norm_median": frame.true_norm.median(),
    }


def main() -> None:
    per = pd.concat([pd.read_csv(BASE / f / "per_perturbation.csv").assign(fold=f) for f in FOLDS])
    scores = pd.concat([pd.read_csv(BASE / f / "scores_raw.csv").assign(fold=f) for f in FOLDS])

    tier_rows = []
    for (fold, arm), g in per[
        per.arm.isin(["T2", "T1", "T0", "ANCHOR_mean_response", "ANCHOR_split_half"])
    ].groupby(["fold", "arm"]):
        tier_rows.append({"fold": fold, "arm": arm, **summarise(g)})
    by_tier = pd.DataFrame(tier_rows)
    by_tier.to_csv(BASE / "anatomy_by_tier.csv", index=False)

    mix = per[per.arm == "MIX_a1.0"].copy()
    mix["magnitude_tercile"] = mix.groupby("fold").true_norm.transform(
        lambda s: pd.qcut(s, 3, labels=["low", "mid", "high"])
    )
    mag_rows = []
    for keys, g in mix.groupby(["fold", "tier", "magnitude_tercile"], observed=True):
        mag_rows.append(
            dict(zip(["fold", "tier", "magnitude_tercile"], keys, strict=True), **summarise(g))
        )
    for keys, g in mix.groupby(["tier", "magnitude_tercile"], observed=True):
        mag_rows.append(
            dict(zip(["tier", "magnitude_tercile"], keys, strict=True), fold="ALL", **summarise(g))
        )
    pd.DataFrame(mag_rows).to_csv(BASE / "anatomy_by_magnitude.csv", index=False)

    ladder = scores[scores.arm.str.startswith("MIX_a")].copy()
    ladder.to_csv(BASE / "ladder_raw.csv", index=False)
    scaled = []
    for fold, g in ladder.groupby("fold"):
        anchors = scores[scores.fold == fold].set_index("arm")
        b, r = anchors.loc["ANCHOR_mean_response"], anchors.loc["ANCHOR_split_half"]
        for _, row in g.iterrows():
            out = {"fold": fold, "amplitude": row.amplitude}
            for m in MEMBERS:
                span = r[m] - b[m]
                val = (row[m] - b[m]) / span if span != 0 else np.nan
                if m == "expr_mse_unbiased_capped_norm":
                    val = float(np.clip(val, 0, 1))  # official clamp
                out[m] = val
            out["overall"] = np.mean([out[m] for m in MEMBERS])
            scaled.append(out)
    scaled = pd.DataFrame(scaled)
    scaled.to_csv(BASE / "ladder_scaled.csv", index=False)
    ladder_mean = scaled.groupby("amplitude")[list(MEMBERS) + ["overall"]].mean()
    ladder_mean.to_csv(BASE / "ladder_scaled_mean.csv")

    mix_all = pd.concat(
        [per[per.arm == f"MIX_a{a}"].assign(a=a) for a in (0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0)]
    )
    ladder_anatomy = pd.DataFrame(
        [{"amplitude": a, **summarise(g)} for a, g in mix_all.groupby("a")]
    )
    ladder_anatomy.to_csv(BASE / "ladder_anatomy.csv", index=False)

    dual = BASE / "replogle22k562__dualmoment"
    generator = None
    if (dual / "scores_raw.csv").exists():
        g1 = scores[scores.fold == "replogle22k562"].set_index("arm")
        dm = pd.read_csv(dual / "scores_raw.csv").set_index("arm")
        dm_per = pd.read_csv(dual / "per_perturbation.csv")
        g1_per = per[per.fold == "replogle22k562"]
        rows = []
        for arm in dm.index.intersection(g1.index):
            rows.append(
                {
                    "arm": arm,
                    **{f"g1__{m}": g1.loc[arm, m] for m in MEMBERS},
                    **{f"dual__{m}": dm.loc[arm, m] for m in MEMBERS},
                    "g1__yield": summarise(g1_per[g1_per.arm == arm])["yield"],
                    "dual__yield": summarise(dm_per[dm_per.arm == arm])["yield"],
                    "g1__precision": summarise(g1_per[g1_per.arm == arm])["precision_where_called"],
                    "dual__precision": summarise(dm_per[dm_per.arm == arm])[
                        "precision_where_called"
                    ],
                }
            )
        generator = pd.DataFrame(rows)
        generator.to_csv(BASE / "generator_ablation_k562.csv", index=False)

    verdict = {
        "pure_tier_means_over_folds": by_tier.groupby("arm")[
            [
                "fid",
                "precision_where_called",
                "yield",
                "share_calling_nothing",
                "pooled_ratio_pred_over_real",
            ]
        ]
        .mean()
        .to_dict(orient="index"),
        "mixed_a1_over_folds": summarise(mix),
        "ladder_scaled_mean": ladder_mean.to_dict(orient="index"),
        "ladder_anatomy": ladder_anatomy.set_index("amplitude")[
            ["fid", "precision_where_called", "yield"]
        ].to_dict(orient="index"),
    }
    (BASE / "fid_decomposition.json").write_text(json.dumps(verdict, indent=2, default=float))
    pd.set_option("display.width", 200)
    print(
        by_tier.groupby("arm")[
            [
                "fid",
                "precision_where_called",
                "yield",
                "share_calling_nothing",
                "n_pred_median",
                "n_real_median",
                "reach",
                "jaccard",
                "nmae",
            ]
        ]
        .mean()
        .round(3)
    )
    print(
        ladder_anatomy[
            [
                "amplitude",
                "fid",
                "precision_where_called",
                "yield",
                "n_pred_median",
                "pred_norm_median",
            ]
        ].round(3)
    )
    print(ladder_mean.round(3))
    print(
        pd.DataFrame(mag_rows)
        .query("fold == 'ALL'")[
            ["tier", "magnitude_tercile", "n", "fid", "precision_where_called", "yield"]
        ]
        .round(3)
    )
    if generator is not None:
        print(generator.round(3).T)


if __name__ == "__main__":
    main()
