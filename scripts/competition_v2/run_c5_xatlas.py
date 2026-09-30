"""C5: C1 + X-Atlas, gated on written permission. Prepared, not run while PENDING.

Predeclaration: ``reports/competition_v2/c5_predeclaration.md`` (SHA-256 checked).
The permission gate is checked **before any X-Atlas file is opened**: the status file must
say APPROVED and ``licensing.STATUS`` must list HCT116 and HEK293T as GREEN. Otherwise
the script prints the gate state and exits 0 without reading data.

    uv run python scripts/competition_v2/run_c5_xatlas.py [--jobs 5]

Outputs (only once permitted): ``outputs/competition_v2/c5_xatlas/``.
"""

from __future__ import annotations

import os

os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import argparse  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import re  # noqa: E402
import sys  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from virtual_cell.competition_v2 import licensing  # noqa: E402

PREDECL = ROOT / "reports" / "competition_v2" / "c5_predeclaration.md"
PREDECL_SHA = "0241d20463ef88567663248a60701f3b33b2cb4d04554e1b62f6df3091bea425"
STATUS_FILE = ROOT / "reports" / "competition_v2" / "xatlas_permission_status.md"
OUT = ROOT / "outputs" / "competition_v2" / "c5_xatlas"
XATLAS = ("HCT116", "HEK293T")


def permission_status(text: str | None = None) -> str:
    """The status recorded in ``xatlas_permission_status.md`` (PENDING/APPROVED/DENIED)."""
    text = STATUS_FILE.read_text() if text is None else text
    m = re.search(r"\|\s*\*\*status\*\*\s*\|\s*\*\*(\w+)\*\*\s*\|", text)
    if not m:
        raise ValueError("no status row in the X-Atlas permission file")
    return m.group(1).upper()


def gate(status: str | None = None, table: dict | None = None) -> tuple[bool, str]:
    """``(open, reason)``: open only if APPROVED **and** both lines are GREEN in code."""
    status = permission_status() if status is None else status
    table = licensing.STATUS if table is None else table
    if status != "APPROVED":
        return False, f"X-Atlas permission is {status}"
    bad = [n for n in XATLAS if table.get(n) != licensing.GREEN]
    if bad:
        return False, f"APPROVED in the status file but not GREEN in licensing.STATUS: {bad}"
    return True, "APPROVED"


def run(jobs: int) -> None:
    """The predeclared comparison C1a vs C5 on the H1 / K562 / CD4 folds."""
    import run_c1_public_folds as c1f
    import run_c3_fusion as c3
    import run_c4_kolf as c4

    from virtual_cell.competition_v2 import evaluation, fusion_c3, generator

    OUT.mkdir(parents=True, exist_ok=True)
    green, cd4, xatlas = c1f.load_green(), c1f.load_cd4(), c1f.load_xatlas()

    def predictors(atlas, with_x):
        srcs = {k: v for k, v in green.items() if k != atlas}
        if with_x:
            srcs.update(xatlas)
        use_cd4 = cd4 if atlas != "CD4" else None
        licensing.assert_sources_allowed(list(srcs) + (["CD4"] if use_cd4 is not None else []))
        return srcs, use_cd4

    mean_rows = []
    for atlas in ("H1", "K562", "CD4"):
        targets, genes = c3.panel(atlas)
        t_eff, t_mask = c3.truth(atlas, targets, genes, green, cd4)
        valid = t_mask & ~np.isin(genes, targets)[None, :]
        for arm, with_x in (("C1a", False), ("C5", True)):
            comp = fusion_c3.components(*predictors(atlas, with_x), targets, genes)
            pred = fusion_c3.combine(comp, "log2fc", None)
            anym = np.zeros_like(valid)
            for _, _, m in fusion_c3.source_list(comp):
                anym |= m
            per, sse, sst = fusion_c3.row_metrics(pred, t_eff, valid & anym)
            row = {"heldout": atlas, "arm": arm, **fusion_c3.summarise(per, sse, sst)}
            if atlas == "CD4":
                from virtual_cell.arc import metrics as M

                row["effect_pds"] = float(
                    M.pds_cosine(pred, t_eff, exclude=np.isin(genes, targets)).mean()
                )
            mean_rows.append(row)
    mean = pd.DataFrame(mean_rows)
    mean.to_csv(OUT / "mean_level.csv", index=False)

    import analyse_c2 as a2

    vcc = {}
    for name in ("H1", "K562"):
        fold = c1f.fold_h1(False) if name == "H1" else c1f.fold_k562(False)
        scorer = evaluation.FoldScorer(fold, jobs=jobs)
        pairs = generator.promoter_pairs(
            generator.gencode_tss(c1f.RAW / "gencode.v47.annotation.gtf.gz"),
            fold.targets,
            fold.genes,
        )
        rows = {}
        (rows["ANCHOR_split_half"], _), _ = scorer.split_half(generator.SEED + 2)
        real_cpm = np.stack([(b / b.sum(1, keepdims=True)).mean(0) for b in fold.real])
        real_bulk = np.stack([b.sum(0) / b.sum() for b in fold.real])
        n = len(fold.targets)
        moments = {
            "ANCHOR_mean_response": (
                np.tile(real_cpm.mean(0), (n, 1)),
                np.tile(real_bulk.mean(0), (n, 1)),
            )
        }
        for arm, with_x in (("C1a", False), ("C5", True)):
            comp = fusion_c3.components(*predictors(name, with_x), fold.targets, fold.genes)
            eff = fusion_c3.fused(comp, fold.controls)
            moments[arm] = generator.expected_moments(
                eff, fold.controls, fold.targets, fold.genes, pairs=pairs
            )
        for arm, (p_cpm, p_bulk) in moments.items():
            st = scorer.emit(
                arm, p_cpm / p_cpm.sum(1, keepdims=True), p_bulk / p_bulk.sum(1, keepdims=True)
            )
            rows[arm], _ = scorer.members(st)
        raw = pd.DataFrame(rows).T
        c2 = pd.read_csv(c4.C2_FOLDS / name / "phase1_scores_raw.csv", index_col=0)
        if (
            max(
                abs(float(raw.loc["C1a", m]) - float(c2.loc["G0_a1.00", m]))
                for m in evaluation.MEMBERS
            )
            > 1e-12
        ):
            raise SystemExit("C1a reproduction failed")
        raw.to_csv(OUT / f"scores_raw_{name}.csv")
        vcc[name] = a2.scale(raw, "ANCHOR_mean_response")

    ml = mean.set_index(["heldout", "arm"])
    gains = {f: float(vcc[f].loc["C5", "Overall"] - vcc[f].loc["C1a", "Overall"]) for f in vcc}
    total = sum(gains.values())
    pds = {f: float(vcc[f].loc["C5", "PDS"] / vcc[f].loc["C1a", "PDS"]) for f in vcc}
    pds["CD4"] = float(ml.loc[("CD4", "C5"), "effect_pds"] / ml.loc[("CD4", "C1a"), "effect_pds"])
    cos = {
        a: float(np.mean([ml.loc[(h, a), "cosine"] for h in ("H1", "K562", "CD4")]))
        for a in ("C1a", "C5")
    }
    improved = [f for f in vcc if gains[f] > 0] + (
        ["CD4"]
        if ml.loc[("CD4", "C5"), "effect_pds"] > ml.loc[("CD4", "C1a"), "effect_pds"]
        else []
    )
    crit = {
        "1_mean_overall_up": float(np.mean(list(gains.values()))) > 0,
        "2_pds_98pct_each_fold": all(v >= 0.98 for v in pds.values()),
        "3_cosine_not_down": cos["C5"] >= cos["C1a"],
        "4_two_folds_improve": len(improved) >= 2,
        "5_no_fold_over_75pct": bool(total > 0 and all(g / total <= 0.75 for g in gains.values())),
    }
    decision = {
        "c5_pass": all(crit.values()),
        "criteria": crit,
        "fold_gain": gains,
        "pds_retention": pds,
        "cosine": cos,
        "improved_folds": improved,
    }
    (OUT / "c5_decision.json").write_text(json.dumps(decision, indent=2))
    print(json.dumps(decision, indent=1))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=5)
    args = parser.parse_args()
    if hashlib.sha256(PREDECL.read_bytes()).hexdigest() != PREDECL_SHA:
        sys.exit("C5 predeclaration changed after it was frozen")
    ok, reason = gate()
    if not ok:
        print(f"C5 gated, not run: {reason}")
        return 0
    run(args.jobs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
