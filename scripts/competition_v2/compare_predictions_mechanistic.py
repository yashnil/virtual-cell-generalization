"""Mechanistic comparison of two Arc A/B/C prediction bundles (section 12 A-D, F, G).

Reads the emitted **counts** of each bundle (no model internals), so both are measured
in the scorer's own space. Per (context, target):

* pseudobulk effect ``delta = log1p(5e4 * P / sum P) - control profile`` (PDS space),
  with all 300 panel targets removed, exactly as ``pds_cosine`` does;
* effect norm, and norm of the deviation from the context's panel-mean prediction;
* per-cell library size, detected genes, within-group heterogeneity.

"Carries perturbation-specific signal" is defined by **reproducibility across contexts**,
which is indifferent to how noisy a generator's emission is: for each target, the mean
cosine between its deviation-from-panel-mean in A, B and C. Emission noise is
independent between contexts and averages to 0; a target-specific signal carried into
every context does not. The null is V1's Tier-0 targets (one shared expected profile,
so their deviations are pure sampling noise); a target is target-specific when its
reproducibility exceeds that null's 99th percentile. (A norm threshold is not usable:
the two generators have noise floors differing by ~3x.)

Writes ``outputs/competition_v2/mechanistic_comparison/``.

Reproduce: ``uv run python scripts/competition_v2/compare_predictions_mechanistic.py``
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[2]
CONTROLS = ROOT / "data" / "raw" / "arc2026" / "controls"
V1_TIERS = ROOT / "data" / "splits" / "arc_target_support_v1.csv"
COVERAGE = ROOT / "data" / "splits" / "arc_target_support_competition_v2.csv"
OUT = ROOT / "outputs" / "competition_v2" / "mechanistic_comparison"
BUNDLES = {
    "V1": ROOT / "outputs" / "arc_dry_run_v1" / "arc_dry_run_v1.h5ad",
    "C0": ROOT / "outputs" / "competition_v2" / "atlasshift_c0" / "prediction.h5ad",
}
TS = 5.0e4
CONTEXTS = ("A", "B", "C")


def bulk(p: np.ndarray) -> np.ndarray:
    return np.log1p(TS * p / p.sum(axis=-1, keepdims=True))


def summarise(path: Path, genes: np.ndarray, targets: list[str]) -> dict:
    """One streaming pass: per-group count sums and per-cell structure."""
    a = ad.read_h5ad(path, backed="r")
    obs = a.obs[["target_gene", "context"]].astype(str)
    if not np.array_equal(a.var_names.astype(str), genes):
        raise ValueError(f"{path.name}: gene order differs from gene_names.csv")
    a.file.close()
    t_idx = pd.Index(targets).get_indexer(obs.target_gene)
    c_idx = pd.Index(CONTEXTS).get_indexer(obs.context)
    group = c_idx * len(targets) + t_idx
    n_groups = len(CONTEXTS) * len(targets)
    sums = np.zeros((n_groups, len(genes)))
    lib = np.zeros(len(obs))
    detected = np.zeros(len(obs))
    hetero = np.zeros(len(obs))
    with h5py.File(path, "r") as f:
        x = f["X"]
        indptr = x["indptr"][:]
        step = 4000  # a multiple of the 400-cell group size
        for left in range(0, len(obs), step):
            right = min(left + step, len(obs))
            lo, hi = indptr[left], indptr[right]
            block = sparse.csr_matrix(
                (x["data"][lo:hi], x["indices"][lo:hi], indptr[left : right + 1] - lo),
                shape=(right - left, len(genes)),
            )
            g = group[left:right]
            member = sparse.csr_matrix(
                (np.ones(right - left), (g, np.arange(right - left))),
                shape=(n_groups, right - left),
            )
            sums += (member @ block).toarray()
            lib[left:right] = np.asarray(block.sum(axis=1)).ravel()
            detected[left:right] = np.diff(block.indptr)
            logcp = block.multiply(1e4 / np.maximum(lib[left:right], 1)[:, None]).tocsr()
            logcp.data = np.log1p(logcp.data)
            dense = logcp.toarray()
            for gid in np.unique(g):
                rows = g == gid
                centre = dense[rows].mean(axis=0)
                hetero[left:right][rows] = np.linalg.norm(dense[rows] - centre, axis=1)
    return {"sums": sums, "lib": lib, "detected": detected, "hetero": hetero, "group": group}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundles", nargs="*", default=list(BUNDLES))
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(CONTROLS / "gene_names.csv").gene_name.astype(str).to_numpy()
    targets = pd.read_csv(CONTROLS / "pert_counts.csv").target_gene.astype(str).tolist()
    keep = ~np.isin(genes, targets)
    tiers = pd.read_csv(V1_TIERS).set_index("arc_target").support_tier.reindex(targets)
    tier0 = (tiers == 0).to_numpy()

    ctrl_profile, ctrl_cells = {}, {}
    for c in CONTEXTS:
        ctrl = ad.read_h5ad(CONTROLS / f"context_{c}.h5ad")
        x = sparse.csr_matrix(ctrl.X)
        ctrl_profile[c] = bulk(np.asarray(x.sum(axis=0)).ravel())
        lib = np.asarray(x.sum(axis=1)).ravel()
        ctrl_cells[c] = {
            "library_median": float(np.median(lib)),
            "detected_median": float(np.median(np.diff(x.indptr))),
        }

    per_target, per_context, structure, deltas = [], [], [], {}
    for name in args.bundles:
        s = summarise(BUNDLES[name], genes, targets)
        delta = np.stack(
            [
                bulk(s["sums"][ci * 300 : (ci + 1) * 300]) - ctrl_profile[c]
                for ci, c in enumerate(CONTEXTS)
            ]
        )[:, :, keep]
        deltas[name] = delta
        for ci, c in enumerate(CONTEXTS):
            d = delta[ci]
            dev = d - d.mean(axis=0, keepdims=True)
            unit = d / np.maximum(np.linalg.norm(d, axis=1, keepdims=True), 1e-12)
            udev = dev / np.maximum(np.linalg.norm(dev, axis=1, keepdims=True), 1e-12)
            cos = unit @ unit.T
            cos_dev = udev @ udev.T
            off = ~np.eye(300, dtype=bool)
            sv = np.linalg.svd(dev, compute_uv=False)
            p = sv**2 / (sv**2).sum()
            per_context.append(
                {
                    "bundle": name,
                    "context": c,
                    "effect_norm_median": float(np.median(np.linalg.norm(d, axis=1))),
                    "deviation_norm_median": float(np.median(np.linalg.norm(dev, axis=1))),
                    "pairwise_cosine_effect_mean": float(cos[off].mean()),
                    "pairwise_cosine_deviation_mean": float(cos_dev[off].mean()),
                    "deviation_effective_rank": float(np.exp(-(p * np.log(p + 1e-300)).sum())),
                    "between_target_variance_share": float(
                        (dev**2).sum() / max((d**2).sum(), 1e-12)
                    ),
                }
            )
            for ti, t in enumerate(targets):
                per_target.append(
                    {
                        "bundle": name,
                        "context": c,
                        "target": t,
                        "v1_tier": int(tiers.iloc[ti]),
                        "effect_norm": float(np.linalg.norm(d[ti])),
                        "deviation_norm": float(np.linalg.norm(dev[ti])),
                    }
                )
            rows = (s["group"] // 300) == ci
            structure.append(
                {
                    "bundle": name,
                    "context": c,
                    "library_median": float(np.median(s["lib"][rows])),
                    "detected_median": float(np.median(s["detected"][rows])),
                    "heterogeneity_median": float(np.median(s["hetero"][rows])),
                    "control_library_median": ctrl_cells[c]["library_median"],
                    "control_detected_median": ctrl_cells[c]["detected_median"],
                }
            )

    pt = pd.DataFrame(per_target)
    repro = {}
    for name, delta in deltas.items():
        dev = delta - delta.mean(axis=1, keepdims=True)
        unit = dev / np.maximum(np.linalg.norm(dev, axis=2, keepdims=True), 1e-12)
        pairs = [(0, 1), (0, 2), (1, 2)]
        repro[name] = np.mean([(unit[i] * unit[j]).sum(axis=1) for i, j in pairs], axis=0)
    np.savez_compressed(OUT / "deltas.npz", **deltas)
    null_p99 = float(np.percentile(repro["V1"][tier0], 99)) if "V1" in repro else np.nan
    pt["cross_context_reproducibility"] = [
        repro[b][targets.index(t)] for b, t in zip(pt.bundle, pt.target, strict=True)
    ]
    pt["reproducibility_null_p99"] = null_p99
    pt["target_specific"] = pt.cross_context_reproducibility > null_p99
    # cross-context adaptation: does the same target get the same direction in A/B/C?
    adaptation = []
    for name, delta in deltas.items():
        for (i, a), (j, b) in [((0, "A"), (1, "B")), ((0, "A"), (2, "C")), ((1, "B"), (2, "C"))]:
            da = delta[i] - delta[i].mean(axis=0)
            db = delta[j] - delta[j].mean(axis=0)
            num = (da * db).sum(axis=1)
            den = np.linalg.norm(da, axis=1) * np.linalg.norm(db, axis=1)
            adaptation.append(
                {
                    "bundle": name,
                    "pair": f"{a}-{b}",
                    "same_target_deviation_cosine_median": float(
                        np.median(num / np.maximum(den, 1e-12))
                    ),
                }
            )
    pt.to_csv(OUT / "per_target.csv", index=False)
    pd.DataFrame(per_context).to_csv(OUT / "per_context.csv", index=False)
    pd.DataFrame(structure).to_csv(OUT / "structure.csv", index=False)
    pd.DataFrame(adaptation).to_csv(OUT / "adaptation.csv", index=False)
    summary = {
        "reproducibility_null_p99": null_p99,
        "target_specific_count": {
            b: int(g.target_specific.sum() / len(CONTEXTS)) for b, g in pt.groupby("bundle")
        },
        "target_specific_by_v1_tier": {
            f"{b}|tier{t}": int(g.target_specific.sum() / len(CONTEXTS))
            for (b, t), g in pt.groupby(["bundle", "v1_tier"])
        },
        "reproducibility_median": {b: float(np.median(r)) for b, r in repro.items()},
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(pd.DataFrame(per_context).to_string(index=False))
    print(pd.DataFrame(structure).to_string(index=False))
    print(pd.DataFrame(adaptation).to_string(index=False))
    print(summary)


if __name__ == "__main__":
    main()
