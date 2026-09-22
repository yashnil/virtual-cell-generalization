"""Arc A/B/C dry run — output-space mapping and a candidate bundle (sections J, K).

Everything the model contains was frozen on public data before this script ran:
the main-effect estimator (``M3b_basal_shrunk``), the tier shrinkages
(Tier 2 = 0.50, Tier 1 = 0.25, Tier 0 = 0 exactly) and the generator. Nothing
here is tuned, and no Arc perturbation outcome exists to tune on — the three
context files hold control cells only.

What the script does:

  J. builds ``beta_hat`` for the 86 Arc targets with direct public evidence,
     from the two datasets that actually perturb them, maps the prediction onto
     the official 18,533-gene panel, and applies the unsupported-gene rule;
  K. generates 300 x 400 cells for each of A, B and C, writes them to a single
     ``.h5ad``, and validates it locally with ``vcc prep --dry-run``.

**No submission is made and nothing is uploaded.**

Reproduce: ``uv run python scripts/run_arc_dry_run.py``
"""

from __future__ import annotations

import json
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd

from virtual_cell.arc import bundle, generate
from virtual_cell.arc.panel import Support, build_panel_map
from virtual_cell.data import arc2026, scperteval
from virtual_cell.modelling import context_main_effect as cme
from virtual_cell.modelling import mean_response as mr

ROOT = Path(__file__).resolve().parents[1]
CONTROLS = ROOT / "data" / "raw" / "arc2026" / "controls"
RAW = ROOT / "data" / "raw" / "scperteval"
SPLITS = ROOT / "data" / "splits"
OUTDIR = ROOT / "outputs" / "arc_count_space_v1"
BUNDLE_DIR = ROOT / "outputs" / "arc_dry_run_v1"

#: The public datasets that perturb any Arc target. They are the only possible
#: source of a direct ``beta_hat`` for this panel, and also the only public
#: contexts whose perturbation panel overlaps Arc's at all, so they are used
#: for ``m_hat`` too. Stated as a limitation in the report: their main effects
#: average over their own panels, not over Arc's 300 targets.
BETA_SOURCES = ("arch1", "kaden25rpe1")

#: Frozen on public leave-one-context-out folds by
#: ``scripts/run_arc_count_space_baseline.py``. Unanimous across all four folds.
MAIN_EFFECT = "M3b_basal_shrunk"
TIER_WEIGHTS = {2: 0.50, 1: 0.25, 0: 0.0}

CELLS_PER_PERT = 400
SEED = 20260921
CHUNK = 20_000


def rule(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


# --------------------------------------------------------------------------
# source responses
# --------------------------------------------------------------------------


def source_pseudobulk(name: str, genes: list[str]) -> scperteval.ContextPseudobulk:
    path = RAW / f"{name}_processed_complete.h5ad"
    labels = set(scperteval.cell_labels(path).tolist()) - {scperteval.CONTROL_LABEL}
    return scperteval.pseudobulk(
        path, genes=genes, perturbations=sorted(labels), chunk_size=CHUNK, context=name
    )


def main_effect_tensor(mains: list[np.ndarray], arc_main: np.ndarray) -> np.ndarray:
    """Shape the per-source main effects for the frozen ``m_hat`` estimators.

    The estimators in :mod:`virtual_cell.modelling.context_main_effect` read a
    ``(context, perturbation, gene)`` tensor but only ever through
    ``delta[c].mean(axis=0)``. The public sources have *different* perturbation
    panels, so no shared perturbation axis exists — and none is needed. Giving
    each context a single row equal to its own main effect makes that mean the
    identity, which is exactly the quantity the estimators want.
    """
    return np.stack([*mains, arc_main])[:, None, :]


# --------------------------------------------------------------------------
# control cells
# --------------------------------------------------------------------------


def load_arc_controls(context: str, n_genes: int) -> np.ndarray:
    """Dense float64 control matrix for one Arc context."""
    path = CONTROLS / f"context_{context}.h5ad"
    out = np.zeros((arc2026.read_shape(path)[0], n_genes), dtype=np.float64)
    for start, chunk in arc2026.stream_row_chunks(path, chunk_size=CHUNK):
        out[start : start + chunk.shape[0]] = chunk.toarray()
    return out


def log1p_cp10k_profile(counts: np.ndarray) -> np.ndarray:
    lib = counts.sum(axis=1, keepdims=True)
    return np.log1p(
        bundle.TARGET_SUM * np.divide(counts, lib, out=np.zeros_like(counts), where=lib > 0)
    ).mean(axis=0)


# --------------------------------------------------------------------------


def main() -> None:
    started = time.time()
    BUNDLE_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    rule("J1. THE OFFICIAL OUTPUT SPACE")
    arc_genes = list(arc2026.load_gene_names(CONTROLS))
    arc_targets = list(arc2026.load_pert_counts(CONTROLS)["target_gene"])
    manifest = arc2026.load_manifest(CONTROLS)
    contexts = list(manifest.contexts)
    print(f"  panel: {len(arc_genes)} genes, {len(arc_targets)} targets, contexts {contexts}")

    support = pd.read_csv(SPLITS / "arc_target_support_v1.csv")
    tier_of = dict(zip(support["arc_target"], support["support_tier"], strict=True))
    sources_of = {
        row.arc_target: [s for s in str(row.datasets_perturbed).split("|") if s]
        for row in support.itertuples()
    }
    for tier in (2, 1, 0):
        print(f"  Tier {tier}: {sum(1 for t in arc_targets if tier_of[t] == tier)} targets")

    rule("J2. SOURCE RESPONSES FOR THE ARC PANEL")
    gene_sets = []
    for name in BETA_SOURCES:
        path = RAW / f"{name}_processed_complete.h5ad"
        gene_sets.append(set(scperteval.read_var_names(path)))
    response_genes = sorted(set(arc_genes).intersection(*gene_sets))
    print(f"  response space (panel n all beta sources): {len(response_genes)} genes")

    pseudobulks = {}
    for name in BETA_SOURCES:
        t0 = time.time()
        pseudobulks[name] = source_pseudobulk(name, response_genes)
        pb = pseudobulks[name]
        print(
            f"  {name:14s} {len(pb.perturbations):>5} perturbations, "
            f"{pb.n_control_cells:>7,} control cells  ({time.time() - t0:.0f}s)"
        )

    mains = [pb.delta.mean(axis=0) for pb in pseudobulks.values()]
    betas = {}
    for name, pb in pseudobulks.items():
        centred = pb.delta - pb.delta.mean(axis=0, keepdims=True)
        betas[name] = dict(zip(pb.perturbations, centred, strict=True))

    beta_hat = np.zeros((len(arc_targets), len(response_genes)), dtype=np.float64)
    weight = np.zeros(len(arc_targets))
    for i, target in enumerate(arc_targets):
        available = [s for s in sources_of[target] if s in betas and target in betas[s]]
        if not available:
            continue
        beta_hat[i] = np.mean([betas[s][target] for s in available], axis=0)
        weight[i] = TIER_WEIGHTS[tier_of[target]]
    print(f"  beta_hat non-zero for {int((weight > 0).sum())} of {len(arc_targets)} targets")
    drift = float(np.linalg.norm((weight[:, None] * beta_hat).mean(axis=0)))
    print(f"  panel-mean drift induced by the tiered beta: ||mean_p w*beta|| = {drift:.5f}")

    rule("J2b. WHAT IS m_hat ACTUALLY AVERAGING OVER? (diagnostic, not used)")
    # m_c is the mean response over ARC's 300 targets in the Arc context. Each
    # public source can only offer the mean over ITS OWN panel, and those panels
    # are not Arc's. This measures the size of that substitution, because if the
    # two differ a lot then m_hat is estimating the wrong quantity, however well
    # the estimator family was selected.
    diag_rows = []
    for name, pb in pseudobulks.items():
        full = pb.delta.mean(axis=0)
        rows = [i for i, p in enumerate(pb.perturbations) if p in set(arc_targets)]
        on_arc = pb.delta[rows].mean(axis=0) if rows else np.zeros_like(full)
        denom = np.linalg.norm(full) * np.linalg.norm(on_arc)
        diag_rows.append(
            {
                "source": name,
                "n_perturbations_all": len(pb.perturbations),
                "n_perturbations_on_arc_panel": len(rows),
                "main_effect_norm_all": float(np.linalg.norm(full)),
                "main_effect_norm_arc_subset": float(np.linalg.norm(on_arc)),
                "cosine_all_vs_arc_subset": float(full @ on_arc / denom) if denom > 0 else np.nan,
            }
        )
    diagnostics = pd.DataFrame(diag_rows)
    print(diagnostics.to_string(index=False))

    def cosine(a: np.ndarray, b: np.ndarray) -> float:
        denom = np.linalg.norm(a) * np.linalg.norm(b)
        return float(a @ b / denom) if denom > 0 else np.nan

    cross = cosine(mains[0], mains[1])

    # Do the sources disagree about the perturbation response, or only about
    # which perturbations they measured? Averaging each source over the SAME
    # perturbations separates the two, and decides whether a matched-panel
    # m_hat is worth building.
    names = list(pseudobulks)
    shared = sorted(
        set(pseudobulks[names[0]].perturbations)
        & set(pseudobulks[names[1]].perturbations)
        & set(arc_targets)
    )
    matched = []
    for name in names:
        pb = pseudobulks[name]
        index = {p: i for i, p in enumerate(pb.perturbations)}
        matched.append(pb.delta[[index[p] for p in shared]].mean(axis=0))
    cross_matched = cosine(matched[0], matched[1]) if shared else np.nan
    cross_arc_panel = cosine(
        *(
            pseudobulks[name]
            .delta[
                [i for i, p in enumerate(pseudobulks[name].perturbations) if p in set(arc_targets)]
            ]
            .mean(axis=0)
            for name in names
        )
    )
    print(
        f"  cosine over each source's OWN panel:            {cross:.4f}\n"
        f"  cosine over each source's Arc-target subset:    {cross_arc_panel:.4f}"
        f"  ({diagnostics['n_perturbations_on_arc_panel'].tolist()} targets)\n"
        f"  cosine over the {len(shared)} targets BOTH measure:        {cross_matched:.4f}"
    )
    fitted_scale = mr.fit_main_effect_scale(
        main_effect_tensor(mains, np.zeros(len(response_genes))),
        tuple(range(len(BETA_SOURCES))),
        cme.source_pooled_mean,
    )
    print(f"  leave-one-source-out scale the estimator fits: {fitted_scale:.4f}")
    diagnostics.to_csv(BUNDLE_DIR / "main_effect_diagnostic.csv", index=False)

    rule("J3. PANEL MAPPING")
    panel_map = build_panel_map(arc_genes, response_genes)
    counts = panel_map.counts()
    print("  " + counts.to_string().replace("\n", "\n  "))
    predicted_mask = panel_map.mask(Support.PREDICTED)
    #: Where each response-space gene sits on the panel axis, so the Arc basal
    #: profile can be restated in the response space the sources share.
    response_positions = pd.Index(arc_genes).get_indexer(pd.Index(response_genes))
    if (response_positions < 0).any():
        raise SystemExit("response space contains a gene absent from the panel")
    print(
        f"  {int(predicted_mask.sum())} panel genes can carry a response; "
        f"{int((~predicted_mask).sum())} keep their target-context control distribution"
    )

    rule("K. GENERATE AND WRITE THE DRY-RUN BUNDLE")
    out_path = BUNDLE_DIR / "arc_dry_run_v1.h5ad"
    per_context_rows = []
    writer = bundle.SubmissionWriter(out_path, arc_genes)
    try:
        for context in contexts:
            t0 = time.time()
            ctrl = load_arc_controls(context, len(arc_genes))
            basal_panel = log1p_cp10k_profile(ctrl)
            basal_response = basal_panel[response_positions]
            # Source and Arc basal profiles on the shared response space, for
            # the basal-similarity weighting inside m_hat.
            control_means = np.stack(
                [pseudobulks[n].control_mean for n in BETA_SOURCES] + [basal_response]
            )
            tensor = main_effect_tensor(mains, np.zeros(len(response_genes)))
            m_hat = mr.MAIN_EFFECT_ESTIMATORS[MAIN_EFFECT](
                tensor,
                tuple(range(len(BETA_SOURCES))),
                control_means=control_means,
                target=len(BETA_SOURCES),
            )
            delta_hat = m_hat[None, :] + weight[:, None] * beta_hat
            lfc = bundle.panel_log2_fold_change(panel_map, basal_panel, delta_hat)

            for i, target in enumerate(arc_targets):
                block = generate.transport_controls(
                    ctrl, lfc[i], CELLS_PER_PERT, rng=rng, smoothing=0.5
                )
                writer.append(block, target_gene=target, context=context)
            per_context_rows.append(
                {
                    "context": context,
                    "control_cells": int(ctrl.shape[0]),
                    "m_hat_norm": float(np.linalg.norm(m_hat)),
                    "median_abs_lfc": float(np.median(np.abs(lfc))),
                    "max_abs_lfc": float(np.abs(lfc).max()),
                    "seconds": round(time.time() - t0, 1),
                }
            )
            print(
                f"  {context}: {CELLS_PER_PERT * len(arc_targets):,} cells written "
                f"in {time.time() - t0:.0f}s  (running nnz {writer.nnz:,})"
            )
            del ctrl
    finally:
        writer.close()
    print(f"  wrote {out_path} ({out_path.stat().st_size / 2**30:.1f} GiB)")
    pd.DataFrame(per_context_rows).to_csv(BUNDLE_DIR / "per_context.csv", index=False)

    rule("K2. LOCAL VALIDATION OF THE WRITTEN BUNDLE")
    report = bundle.inspect_bundle(out_path, arc_genes)
    for key, value in report.as_dict().items():
        print(f"  {key:32s} {value}")
    group_sizes = report.cells_per_group
    checks = {
        "n_cells == 360000": report.n_cells == 360_000,
        "n_genes == 18533": report.n_genes == 18_533,
        "gene order matches gene_names.csv": report.gene_order_matches,
        "raw integer counts": report.integer_valued,
        "non-negative": report.non_negative,
        "finite": report.finite,
        "exactly 400 cells per (context, perturbation)": bool((group_sizes == 400).all()),
        "300 perturbations in every context": all(
            v == 300 for v in report.n_perturbations_per_context.values()
        ),
        "contexts are exactly A, B, C": report.contexts == ("A", "B", "C"),
        "no control cells emitted": "non-targeting"
        not in set(group_sizes.index.get_level_values(1)),
        "under max_counts_per_cell (1,000,000)": report.max_counts_per_cell
        < generate.MAX_COUNTS_PER_CELL,
        "under max_nnz (4,750,000,000)": report.nnz < 4_750_000_000,
        "under max_cell_dim (400,000)": report.n_cells <= 400_000,
    }
    for label, ok in checks.items():
        print(f"  [{'PASS' if ok else 'FAIL'}] {label}")

    rule("K3. vcc prep --dry-run")
    cmd = [
        "vcc",
        "prep",
        str(out_path),
        "-g",
        str(CONTROLS / "gene_names.csv"),
        "--perts",
        str(CONTROLS / "pert_counts.csv"),
        "--dry-run",
        "--json",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    print("  $ " + " ".join(cmd))
    print(f"  exit code: {proc.returncode}")
    print(proc.stdout.strip()[:4000] or "(no stdout)")
    if proc.stderr.strip():
        print("  stderr:\n" + proc.stderr.strip()[:2000])
    (BUNDLE_DIR / "vcc_prep_dry_run.json").write_text(proc.stdout)

    summary = {
        "bundle": str(out_path),
        "bundle_bytes": out_path.stat().st_size,
        "main_effect_estimator": MAIN_EFFECT,
        "tier_weights": {str(k): v for k, v in TIER_WEIGHTS.items()},
        "beta_sources": list(BETA_SOURCES),
        "response_space_genes": len(response_genes),
        "panel_support": {k: int(v) for k, v in counts.items()},
        "targets_with_beta": int((weight > 0).sum()),
        "panel_mean_drift": drift,
        "main_effect_diagnostic": diagnostics.to_dict(orient="records"),
        "source_main_effect_cosine": cross,
        "source_main_effect_cosine_arc_subset": cross_arc_panel,
        "source_main_effect_cosine_matched": cross_matched,
        "n_matched_targets": len(shared),
        "fitted_main_effect_scale": fitted_scale,
        "report": report.as_dict(),
        "checks": checks,
        "vcc_prep_returncode": proc.returncode,
        "submitted": False,
        "elapsed_seconds": round(time.time() - started, 1),
    }
    (BUNDLE_DIR / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    np.save(BUNDLE_DIR / "beta_hat.npy", beta_hat.astype(np.float32))
    (BUNDLE_DIR / "response_genes.txt").write_text("\n".join(response_genes) + "\n")
    print(f"\nwritten to {BUNDLE_DIR}  ({summary['elapsed_seconds']}s)")
    print("NO SUBMISSION WAS MADE.")


if __name__ == "__main__":
    main()
