"""Signature diversity and cell structure of the C1 Arc bundle, beside V1 and C0.

Reuses the definitions of ``scripts/competition_v2/compare_predictions_mechanistic.py``
(imported, unchanged) on the C1 ``prediction.h5ad``. V1 and C0 rows are read from that
phase's frozen outputs (``outputs/competition_v2/mechanistic_comparison/``), never
recomputed or overwritten. Cross-context reproducibility uses the same V1 Tier-0 null.

Writes ``outputs/competition_v2/c1_license_clean/mechanistic/``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "competition_v2"))
import compare_predictions_mechanistic as mech  # noqa: E402

PRIOR = ROOT / "outputs" / "competition_v2" / "mechanistic_comparison"
OUT = ROOT / "outputs" / "competition_v2" / "c1_license_clean" / "mechanistic"
C1 = ROOT / "outputs" / "competition_v2" / "c1_license_clean" / "prediction.h5ad"
COVERAGE = ROOT / "data" / "splits" / "arc_target_support_license_clean_v1.csv"


def diversity(delta: np.ndarray) -> list[dict]:
    rows = []
    for ci, c in enumerate(mech.CONTEXTS):
        d = delta[ci]
        dev = d - d.mean(axis=0, keepdims=True)
        unit = d / np.maximum(np.linalg.norm(d, axis=1, keepdims=True), 1e-12)
        udev = dev / np.maximum(np.linalg.norm(dev, axis=1, keepdims=True), 1e-12)
        off = ~np.eye(len(d), dtype=bool)
        sv = np.linalg.svd(dev, compute_uv=False)
        p = sv**2 / (sv**2).sum()
        rows.append(
            {
                "context": c,
                "effect_norm_median": float(np.median(np.linalg.norm(d, axis=1))),
                "deviation_norm_median": float(np.median(np.linalg.norm(dev, axis=1))),
                "pairwise_cosine_effect_mean": float((unit @ unit.T)[off].mean()),
                "pairwise_cosine_deviation_mean": float((udev @ udev.T)[off].mean()),
                "deviation_effective_rank": float(np.exp(-(p * np.log(p + 1e-300)).sum())),
                "between_target_variance_share": float((dev**2).sum() / max((d**2).sum(), 1e-12)),
            }
        )
    return rows


def reproducibility(delta: np.ndarray) -> np.ndarray:
    dev = delta - delta.mean(axis=1, keepdims=True)
    unit = dev / np.maximum(np.linalg.norm(dev, axis=2, keepdims=True), 1e-12)
    return np.mean([(unit[i] * unit[j]).sum(axis=1) for i, j in [(0, 1), (0, 2), (1, 2)]], axis=0)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(mech.CONTROLS / "gene_names.csv").gene_name.astype(str).to_numpy()
    targets = pd.read_csv(mech.CONTROLS / "pert_counts.csv").target_gene.astype(str).tolist()
    keep = ~np.isin(genes, targets)
    ctrl_profile = {}
    for c in mech.CONTEXTS:
        x = sparse.csr_matrix(ad.read_h5ad(mech.CONTROLS / f"context_{c}.h5ad").X)
        ctrl_profile[c] = mech.bulk(np.asarray(x.sum(axis=0)).ravel())
    s = mech.summarise(C1, genes, targets)
    delta = np.stack(
        [
            mech.bulk(s["sums"][ci * 300 : (ci + 1) * 300]) - ctrl_profile[c]
            for ci, c in enumerate(mech.CONTEXTS)
        ]
    )[:, :, keep]
    per_context = pd.DataFrame(diversity(delta))
    per_context.insert(0, "bundle", "C1")
    prior = pd.read_csv(PRIOR / "per_context.csv")
    pd.concat([prior, per_context]).to_csv(OUT / "per_context.csv", index=False)

    structure = []
    for ci, c in enumerate(mech.CONTEXTS):
        rows = (s["group"] // 300) == ci
        structure.append(
            {
                "bundle": "C1",
                "context": c,
                "library_median": float(np.median(s["lib"][rows])),
                "detected_median": float(np.median(s["detected"][rows])),
                "heterogeneity_median": float(np.median(s["hetero"][rows])),
            }
        )
    pd.concat([pd.read_csv(PRIOR / "structure.csv"), pd.DataFrame(structure)]).to_csv(
        OUT / "structure.csv", index=False
    )

    null = json.loads((PRIOR / "summary.json").read_text())["reproducibility_null_p99"]
    prior_deltas = np.load(PRIOR / "deltas.npz")
    repro = {name: reproducibility(prior_deltas[name]) for name in prior_deltas.files}
    repro["C1"] = reproducibility(delta)
    cov = pd.read_csv(COVERAGE).set_index("arc_target").reindex(targets)
    per_target = pd.DataFrame(
        {
            "target": targets,
            "n_usable_contexts_c1": cov.n_usable_contexts_c1.to_numpy(),
            **{f"repro_{k}": v for k, v in repro.items()},
            "c1_effect_norm_A": np.linalg.norm(delta[0], axis=1),
        }
    )
    per_target.to_csv(OUT / "per_target.csv", index=False)
    summary = {
        "reproducibility_null_p99": null,
        "target_specific_count": {k: int((v > null).sum()) for k, v in repro.items()},
        "reproducibility_median": {k: float(np.median(v)) for k, v in repro.items()},
        "c1_target_specific_by_support": {
            str(n): int(((repro["C1"] > null) & (per_target.n_usable_contexts_c1 == n)).sum())
            for n in sorted(per_target.n_usable_contexts_c1.unique())
        },
        "c1_unsupported_effect_norm_median": float(
            np.median(per_target.c1_effect_norm_A[per_target.n_usable_contexts_c1 == 0])
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(pd.read_csv(OUT / "per_context.csv").to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
