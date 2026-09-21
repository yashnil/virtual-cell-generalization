"""External generalization diagnostic: WHY did the internal result fail on arch1?

``arch1`` changed two things at once relative to the source benchmark: the
cellular context, and the study/target regime. Either could explain the
collapse. ``kaden25rpe1`` separates them: it is an independent
transcription-factor screen in **RPE1**, and Replogle RPE1 is already a source
context. So Kaden holds the cell line fixed and varies only the study and the
target regime.

    arch1   : new context  + new study/regime  -> failed
    kaden   : SAME context + new study/regime  -> ?

If Kaden's unseen-perturbation prediction works, the arch1 failure is about
cellular context. If Kaden also fails, the failure is about interpolating
biological priors across independent perturbation studies, and context is not
the operative variable.

**Nothing is tuned here.** The protocol is a transcription of the frozen arch1
script, and this script refuses to score Kaden unless it first reproduces the
frozen arch1 numbers exactly.

Usage::

    uv run python scripts/run_external_validation.py
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from virtual_cell.data import scperteval
from virtual_cell.modelling.external_benchmark import (
    ExternalResult,
    run_external_benchmark,
    target_reliability,
)

ROOT = Path(__file__).resolve().parents[1]
DESIGN = ROOT / "data" / "splits" / "four_context_v1"
CANONICAL = ROOT / "data" / "processed" / "four_context_v1"
RAW = ROOT / "data" / "raw" / "scperteval"
INVENTORY = ROOT / "data" / "provenance" / "scperteval" / "public_label_inventory.json"
FROZEN = ROOT / "outputs" / "unseen_perturbation_v1"
OUTDIR = ROOT / "outputs" / "external_validation_v1"

#: Tolerance for the arch1 reproduction check. Anything looser would let a
#: changed protocol slip through; anything tighter trips on float summation
#: order alone.
REPRODUCTION_TOLERANCE = 1e-9


def rule(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def load_sources() -> tuple[np.ndarray, np.ndarray, list[str], list[str], list[str], set[str]]:
    contexts = (DESIGN / "contexts.txt").read_text().split()
    perts = (DESIGN / "shared_perturbations.txt").read_text().split()
    genes = (DESIGN / "shared_genes.txt").read_text().split()
    delta = np.load(CANONICAL / "delta_tensor.npy").astype(np.float64)
    control = np.load(CANONICAL / "control_means.npy").astype(np.float64)
    with INVENTORY.open() as fh:
        inventory = json.load(fh)
    seen_anywhere: set[str] = set()
    for name in contexts:
        seen_anywhere |= set(inventory[name]["perturbations"])
    return delta, control, contexts, perts, genes, seen_anywhere


def score(dataset: str, **kwargs: object) -> ExternalResult:
    path = RAW / f"{dataset}_processed_complete.h5ad"
    if not path.is_file():
        raise FileNotFoundError(f"{path} is missing")
    return run_external_benchmark(dataset, path=path, root=ROOT, **kwargs)


def check_reproduces_arch1(result: ExternalResult) -> None:
    """Refuse to go further unless the protocol is provably the frozen one."""
    frozen = pd.read_csv(FROZEN / "arch1_external.csv").set_index("estimator")
    got = result.unseen.set_index("estimator")
    columns = ["beta_pearson", "beta_unexplained", "response_pearson", "response_unexplained"]
    diff = (got[columns] - frozen[columns]).abs().max().max()
    print(f"  max |difference| against frozen arch1_external.csv: {diff:.3e}")

    frozen_seen = pd.read_csv(FROZEN / "arch1_seen_perturbation_reference.csv")
    seen_diff = result.seen[frozen_seen.columns].to_numpy() - frozen_seen.to_numpy()
    seen_max = float(np.abs(seen_diff).max())
    print(f"  max |difference| against frozen seen-perturbation reference: {seen_max:.3e}")

    with (FROZEN / "arch1_summary.json").open() as fh:
        summary = json.load(fh)
    for name, got_value, want in (
        ("n_unseen", result.n_unseen, summary["n_unseen_perturbations"]),
        ("n_shared_genes", result.n_shared_genes, summary["n_shared_genes"]),
        ("main_effect_estimator", result.main_effect_estimator, summary["main_effect_estimator"]),
    ):
        if got_value != want:
            raise AssertionError(f"arch1 reproduction failed on {name}: {got_value!r} != {want!r}")

    if max(diff, seen_max) > REPRODUCTION_TOLERANCE:
        raise AssertionError(
            "The external protocol no longer reproduces the frozen arch1 result "
            f"(max difference {max(diff, seen_max):.3e}). Kaden was not scored."
        )
    print("  REPRODUCED. The protocol is the frozen one; Kaden may be scored.")


def describe(result: ExternalResult) -> None:
    print(f"\n  eligible unseen perturbations : {result.n_unseen}")
    print(f"  directly measured in sources  : {result.n_seen}")
    print(f"  shared gene axis              : {result.n_shared_genes}")
    print(f"  STRING coverage, training     : {result.prior_coverage_train:.3f}")
    print(f"  STRING coverage, unseen set   : {result.prior_coverage_unseen:.3f}")

    print("\n  feasible context main effect:")
    print(result.main_effect.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print(
        f"  diagnostic: with a perfect scalar s*={result.oracle_scalar:.4f}, "
        f"{100 * result.oracle_scalar_residual:.1f}% still unexplained"
    )

    print("\n  UNSEEN perturbations (frozen estimators):")
    print(result.unseen.to_string(index=False, float_format=lambda v: f"{v:+.4f}"))

    if not result.seen.empty:
        row = result.seen.iloc[0]
        print(f"\n  DIRECTLY MEASURED perturbations ({int(row['n_seen'])}), same context and axis:")
        print(
            f"    raw conserved transfer      r={row['raw_pearson']:+.4f}  "
            f"unexplained={row['raw_unexplained']:.4f}"
        )
        print(
            f"    scaled (s={row['scale']:.3f})           r={row['scaled_pearson']:+.4f}  "
            f"unexplained={row['scaled_unexplained']:.4f}"
        )


def confidence_report(result: ExternalResult) -> pd.DataFrame:
    frame = result.confidence
    rows = []
    for score_name in ("neighbour_agreement", "support_distance"):
        ok = np.isfinite(frame[score_name]) & np.isfinite(frame["beta_pearson"])
        rho = float("nan")
        if ok.sum() >= 20:
            rho = float(
                stats.spearmanr(frame.loc[ok, score_name], frame.loc[ok, "beta_pearson"]).statistic
            )
        rows.append(
            {
                "dataset": result.dataset,
                "score": score_name,
                "n": int(ok.sum()),
                "spearman": rho,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    started = time.time()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    delta, control, contexts, perts, genes, seen_anywhere = load_sources()
    shared = {
        "delta": delta,
        "control_means": control,
        "source_genes": genes,
        "source_perturbations": perts,
        "seen_anywhere": seen_anywhere,
    }
    print(f"source tensor {delta.shape} over {contexts}")

    rule("0. PROTOCOL CHECK — REPRODUCE THE FROZEN arch1 RESULT")
    arch1 = score("arch1", **shared)
    check_reproduces_arch1(arch1)

    rule("STAGE A. kaden25rpe1 — SAME CONTEXT (RPE1), INDEPENDENT STUDY")
    print("  Replogle RPE1 is a source context, so the cell line is NOT held out.")
    print("  What changes is the study and the target regime (transcription factors).")
    kaden = score("kaden25rpe1", **shared)
    describe(kaden)

    kaden.main_effect.to_csv(OUTDIR / "kaden_main_effect.csv", index=False)
    kaden.unseen.to_csv(OUTDIR / "kaden_unseen.csv", index=False)
    kaden.seen.to_csv(OUTDIR / "kaden_seen.csv", index=False)
    kaden.confidence.to_csv(OUTDIR / "kaden_confidence.csv", index=False)

    rule("THE CONTROL THAT DECIDES WHETHER EITHER BENCHMARK CAN ANSWER ANYTHING")
    print("  A target whose own responses do not reproduce caps every method at")
    print("  sqrt(rho), so a negative result there carries no information.\n")
    rel_rows = []
    for result, subsets in (
        (arch1, {"unseen": arch1.unseen_labels, "measured": None}),
        (kaden, {"unseen": kaden.unseen_labels, "measured": None}),
    ):
        path = RAW / f"{result.dataset}_processed_complete.h5ad"
        target_genes = set(scperteval.read_var_names(path))
        axis = [g for g in genes if g in target_genes]
        labels = set(map(str, scperteval.cell_labels(path)))
        measured = sorted(labels & set(perts))
        for label, plist in (("unseen", subsets["unseen"]), ("measured", measured)):
            if not plist:
                continue
            frame = target_reliability(path, genes=axis, perturbations=plist, max_perturbations=300)
            frame.insert(0, "subset", label)
            frame.insert(0, "dataset", result.dataset)
            rel_rows.append(frame)
            rel_rows[-1].to_csv(OUTDIR / f"reliability_{result.dataset}_{label}.csv", index=False)
            print(
                f"  {result.dataset:<14} {label:<9} n={len(frame):<5} "
                f"median cells={frame['n_cells'].median():>6.0f}  "
                f"Spearman-Brown={frame['spearman_brown'].median():+.4f}  "
                f"ceiling={np.sqrt(max(frame['spearman_brown'].median(), 0)):.4f}  "
                f"frac>0.2={float((frame['spearman_brown'] > 0.2).mean()):.3f}"
            )
    reliability = pd.concat(rel_rows, ignore_index=True)
    reliability.to_csv(OUTDIR / "target_reliability.csv", index=False)

    rule("CONFIDENCE SCORES ON BOTH EXTERNAL CONTEXTS")
    conf = pd.concat([confidence_report(arch1), confidence_report(kaden)], ignore_index=True)
    print(conf.to_string(index=False, float_format=lambda v: f"{v:+.4f}"))
    conf.to_csv(OUTDIR / "confidence_comparison.csv", index=False)

    rule("SIDE BY SIDE: WHAT CHANGED, AND WHAT BROKE")

    def knn(r: ExternalResult) -> pd.Series:
        return r.unseen.set_index("estimator").loc["U2_knn"]

    def ridge(r: ExternalResult) -> pd.Series:
        return r.unseen.set_index("estimator").loc["U3_ridge"]

    compare = pd.DataFrame(
        [
            {
                "dataset": r.dataset,
                "context_changed": r.dataset != "kaden25rpe1",
                "study_changed": True,
                "n_unseen": r.n_unseen,
                "U0_beta_unexplained": 1.0,
                "knn_beta_pearson": knn(r)["beta_pearson"],
                "knn_beta_unexplained": knn(r)["beta_unexplained"],
                "ridge_beta_pearson": ridge(r)["beta_pearson"],
                "ridge_beta_unexplained": ridge(r)["beta_unexplained"],
                "measured_transfer_pearson": (
                    float(r.seen.iloc[0]["scaled_pearson"]) if not r.seen.empty else np.nan
                ),
                "measured_transfer_unexplained": (
                    float(r.seen.iloc[0]["scaled_unexplained"]) if not r.seen.empty else np.nan
                ),
                "main_effect_oracle_residual": r.oracle_scalar_residual,
            }
            for r in (arch1, kaden)
        ]
    )
    print(compare.to_string(index=False, float_format=lambda v: f"{v:+.4f}"))
    compare.to_csv(OUTDIR / "external_comparison.csv", index=False)

    with (OUTDIR / "summary.json").open("w") as fh:
        json.dump(
            {
                "arch1_reproduced": True,
                "kaden_n_unseen": kaden.n_unseen,
                "kaden_n_seen": kaden.n_seen,
                "kaden_shared_genes": kaden.n_shared_genes,
                "kaden_knn_beta_unexplained": float(knn(kaden)["beta_unexplained"]),
                "kaden_ridge_beta_unexplained": float(ridge(kaden)["beta_unexplained"]),
                "kaden_main_effect_oracle_residual": kaden.oracle_scalar_residual,
                "runtime_minutes": (time.time() - started) / 60.0,
            },
            fh,
            indent=2,
            sort_keys=True,
        )
    print(f"\nWrote {OUTDIR}   [{(time.time() - started) / 60:.1f} min]")


if __name__ == "__main__":
    main()
