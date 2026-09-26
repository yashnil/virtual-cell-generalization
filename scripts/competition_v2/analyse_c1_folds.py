"""Apply the predeclared C1 selection and pass rule; write the comparison tables.

Reads ``outputs/competition_v2/c1_license_clean/folds/{H1,K562,CD4}`` and writes to
``outputs/competition_v2/c1_license_clean/``:

* ``decision.json`` — nested f* per outer fold, the five pass-rule criteria, the
  selected variant, and the build gate;
* ``comparison_public.csv`` — every arm x fold, raw + local-scaled members + DE anatomy;
* ``xatlas_ablation.csv``, ``promoter_ablation.csv``, ``agreement_transfer.csv``.

Everything is fixed by ``reports/competition_v2/c1_predeclaration.md`` except the
**build gate**, written here before any full-fold result existed. The Arc candidate is
built only if the selected C1 variant beats the local mean-response baseline on public
data: mean local Overall over {H1, K562} > 0, and raw PDS > 0.5 in both folds.
It is a floor, not a selection criterion.

    uv run python scripts/competition_v2/analyse_c1_folds.py
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from virtual_cell.competition_v2 import licensing  # noqa: E402
from virtual_cell.competition_v2.evaluation import MEMBERS  # noqa: E402

BASE = ROOT / "outputs" / "competition_v2" / "c1_license_clean"
FOLDS = BASE / "folds"
PREDECL = ROOT / "reports" / "competition_v2" / "c1_predeclaration.md"
PREDECL_SHA = "fd4efa6f6f69871ec9600ca894465755e726f6f15cc33f45f8d1be4eeda32957"
FULL = ["H1", "K562"]
S2_GRID = [0.0, 0.25, 0.5]
S1_GRID = [0.25, 0.5, 0.75]
MATERIAL = -0.02
DE_CHECK = {
    "NMAE": "de_wilcoxon_lfc_nmae",
    "FID": "de_wilcoxon_direction_fidelity_yield_raw",
    "REACH": "de_wilcoxon_direction_reach_raw",
    "JAC": "de_wilcoxon_sig_jaccard",
}
ANATOMY = [
    "predicted_de_median",
    "real_de_median",
    "directional_precision_pooled",
    "de_yield_median",
]


def load(fold: str):
    raw = pd.read_csv(FOLDS / fold / "scores_raw.csv", index_col="arm")
    scaled = pd.read_csv(FOLDS / fold / "scores_scaled_local.csv", index_col="arm")
    summary = json.loads((FOLDS / fold / "summary.json").read_text())
    if summary["predeclaration_sha256"] != PREDECL_SHA:
        raise SystemExit(f"{fold}: scored under a different predeclaration")
    return raw, scaled, summary


def nested(scaled: dict, prefix: str, grid: list[float]):
    out = {}
    for outer in FULL:
        inner = [f for f in FULL if f != outer]
        inner_score = {
            g: np.mean([scaled[i].loc[f"{prefix}{g}", "overall"] for i in inner]) for g in grid
        }
        best = max(grid, key=lambda g: inner_score[g])  # first max in grid order
        out[outer] = {
            "selected": best,
            "inner_overall": inner_score,
            "outer_arm": f"{prefix}{best}",
        }
    final_score = {
        g: np.mean([scaled[f].loc[f"{prefix}{g}", "overall"] for f in FULL]) for g in grid
    }
    final = max(grid, key=lambda g: final_score[g])
    return out, final, final_score


def main() -> None:
    if hashlib.sha256(PREDECL.read_bytes()).hexdigest() != PREDECL_SHA:
        raise SystemExit("predeclaration changed")
    raw, scaled, summ = {}, {}, {}
    for f in FULL:
        raw[f], scaled[f], summ[f] = load(f)
    cd4 = pd.read_csv(FOLDS / "CD4" / "scores_effect.csv", index_col="arm")
    cd4_summary = json.loads((FOLDS / "CD4" / "summary.json").read_text())

    s2, s2_final, s2_final_score = nested(scaled, "S2_f", S2_GRID)
    s1, s1_final, s1_final_score = nested(scaled, "S1_l", S1_GRID)

    def c1b(fold, member="overall", table=None):
        table = scaled if table is None else table
        return table[fold].loc[s2[fold]["outer_arm"], member]

    crit = {}
    c1b_mean = float(np.mean([c1b(f) for f in FULL]))
    c1a_mean = float(np.mean([scaled[f].loc["C1a", "overall"] for f in FULL]))
    crit["1_mean_overall"] = {"C1b": c1b_mean, "C1a": c1a_mean, "pass": c1b_mean > c1a_mean}
    pds_folds = {
        f: {"C1b": float(c1b(f, "pds_cosine", raw)), "C1a": float(raw[f].loc["C1a", "pds_cosine"])}
        for f in FULL
    }
    pds_folds["CD4"] = {
        "C1b": float(cd4.loc[f"S2_f{s2_final}", "pds_effect"]),
        "C1a": float(cd4.loc["C1a", "pds_effect"]),
    }
    better = sum(v["C1b"] > v["C1a"] for v in pds_folds.values())
    need = math.ceil(0.75 * len(pds_folds))
    crit["2_pds_folds"] = {
        "folds": pds_folds,
        "better": better,
        "required": need,
        "pass": better >= need,
    }
    de = {
        name: float(np.mean([c1b(f, m) - scaled[f].loc["C1a", m] for f in FULL]))
        for name, m in DE_CHECK.items()
    }
    crit["3_no_material_worsening"] = {
        "mean_scaled_delta": de,
        "threshold": MATERIAL,
        "pass": all(v >= MATERIAL for v in de.values()),
    }
    per_fold = {f: float(c1b(f) - scaled[f].loc["C1a", "overall"]) for f in FULL}
    crit["4_not_one_context"] = {
        "overall_delta": per_fold,
        "pass": all(v > 0 for v in per_fold.values()),
    }
    try:
        licensing.assert_sources_allowed(["K562", "H1", "CD4", "GENCODE"])
        lic = True
    except licensing.LicenseError:
        lic = False
    crit["5_no_blocked_data"] = {"c1_sources": ["K562", "H1", "CD4", "GENCODE"], "pass": lic}
    passed = all(c["pass"] for c in crit.values())
    selected = f"C1b_S2_f{s2_final}" if passed else "C1a"
    sel_arm = f"S2_f{s2_final}" if passed else "C1a"
    sel_overall = float(np.mean([scaled[f].loc[sel_arm, "overall"] for f in FULL]))
    gate = {
        "mean_local_overall": sel_overall,
        "pds_raw": {f: float(raw[f].loc[sel_arm, "pds_cosine"]) for f in FULL},
    }
    gate["pass"] = sel_overall > 0 and all(v > 0.5 for v in gate["pds_raw"].values())
    decision = {
        "predeclaration_sha256": PREDECL_SHA,
        "c1b_nested_selection": s2,
        "c1b_final_floor": s2_final,
        "c1b_final_score": s2_final_score,
        "s1_diagnostic_nested": s1,
        "s1_diagnostic_final": s1_final,
        "s1_diagnostic_final_score": s1_final_score,
        "criteria": crit,
        "c1b_passes": passed,
        "selected_variant": selected,
        "selected_arm": sel_arm,
        "build_gate": gate,
    }
    (BASE / "decision.json").write_text(json.dumps(decision, indent=2, default=float))

    # ---------------------------------------------------------------- comparison table
    rows = []
    for f in FULL:
        for arm in raw[f].index:
            r = {"fold": f, "arm": arm}
            r.update({f"raw_{m}": raw[f].loc[arm, m] for m in MEMBERS})
            r.update({f"scaled_{m}": scaled[f].loc[arm, m] for m in MEMBERS})
            r["local_overall"] = scaled[f].loc[arm, "overall"]
            r.update({k: raw[f].loc[arm, k] for k in ANATOMY if k in raw[f].columns})
            rows.append(r)
    comp = pd.DataFrame(rows)
    comp.to_csv(BASE / "comparison_public.csv", index=False)

    # ---------------------------------------------------------------- X-Atlas ablation
    xa = []
    for f in FULL:
        for arm, cov, ms in [
            ("C0_noX", "coverage_c1_any", "mean_sources_c1"),
            ("C0_withX", "coverage_withX_any", "mean_sources_withX"),
        ]:
            xa.append(
                {
                    "fold": f,
                    "arm": arm,
                    "n_targets": summ[f]["n_targets"],
                    "coverage": summ[f][cov],
                    "mean_sources": summ[f][ms],
                    "pds_raw": raw[f].loc[arm, "pds_cosine"],
                    "local_overall": scaled[f].loc[arm, "overall"],
                    **{f"scaled_{k}": scaled[f].loc[arm, m] for k, m in DE_CHECK.items()},
                    **{k: raw[f].loc[arm, k] for k in ANATOMY},
                }
            )
    for arm, cov, ms in [
        ("C0_noX", "coverage_c1_any", "mean_sources_c1"),
        ("C0_withX", "coverage_withX_any", "mean_sources_withX"),
    ]:
        xa.append(
            {
                "fold": "CD4 (effect-level)",
                "arm": arm,
                "n_targets": cd4_summary["n_targets"],
                "coverage": cd4_summary[cov],
                "mean_sources": cd4_summary[ms],
                "pds_raw": cd4.loc[arm, "pds_effect"],
            }
        )
    pd.DataFrame(xa).to_csv(BASE / "xatlas_ablation.csv", index=False)

    # ---------------------------------------------------------------- promoter ablation
    pr = []
    for f in FULL:
        on, off = scaled[f].loc["C1a"], scaled[f].loc["C1a_no_promoter"]
        pr.append(
            {
                "fold": f,
                "promoter_pairs_on_panel": summ[f]["promoter_pairs"],
                "pds_raw_on": raw[f].loc["C1a", "pds_cosine"],
                "pds_raw_off": raw[f].loc["C1a_no_promoter", "pds_cosine"],
                **{f"delta_scaled_{m}": on[m] - off[m] for m in MEMBERS},
                "delta_overall": on["overall"] - off["overall"],
            }
        )
    pd.DataFrame(pr).to_csv(BASE / "promoter_ablation.csv", index=False)

    # ---------------------------------------------------------------- agreement vs transfer
    at = []
    for f in FULL + ["CD4"]:
        pt = pd.read_csv(FOLDS / f / "per_target_agreement.csv")
        if f in FULL:
            per = pd.read_csv(FOLDS / f / "per_perturbation.csv")
            pds = per[per.arm == "C1a"].set_index("target").pds
            pt["pds_count_c1a"] = pt.target.map(pds)
        else:
            pt["pds_count_c1a"] = pt["pds_C1a"]
        ok = pt.agreement.notna() & pt.transfer_cosine_c1a.notna()
        rho, p = stats.spearmanr(pt.agreement[ok], pt.transfer_cosine_c1a[ok])
        rho2, p2 = stats.spearmanr(pt.agreement[ok], pt.pds_count_c1a[ok])
        terc = pd.qcut(pt.agreement[ok], 3, labels=["low", "mid", "high"])
        by = pt[ok].groupby(terc, observed=True).transfer_cosine_c1a.median()
        at.append(
            {
                "fold": f,
                "n_with_agreement": int(ok.sum()),
                "n_targets": len(pt),
                "spearman_agreement_vs_transfer_cosine": rho,
                "p_transfer": p,
                "spearman_agreement_vs_pds": rho2,
                "p_pds": p2,
                **{f"median_transfer_cosine_{k}": v for k, v in by.items()},
            }
        )
    pd.DataFrame(at).to_csv(BASE / "agreement_transfer.csv", index=False)
    print(
        json.dumps(
            {
                k: decision[k]
                for k in ["c1b_final_floor", "c1b_passes", "selected_variant", "build_gate"]
            },
            indent=2,
            default=float,
        )
    )
    for name, c in crit.items():
        print(name, c["pass"])


if __name__ == "__main__":
    main()
