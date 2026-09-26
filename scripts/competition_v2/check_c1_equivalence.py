"""End-to-end check: our backbone, fed C0's five prepared sources with C0's weights,
regenerates rows of the frozen C0 bundle (``prediction.h5ad``) exactly.

Deterministic subset: targets 0, 137, 299 in each of contexts A, B, C (400 cells each).
Writes ``outputs/competition_v2/c1_license_clean/equivalence_c0_bundle.json``.
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
sys.path.insert(0, str(ROOT / "src"))

from virtual_cell.competition_v2 import fusion, generator, sources  # noqa: E402

C0 = ROOT / "outputs" / "competition_v2" / "atlasshift_c0"
DATA = C0 / "data"
NAMES = [
    "K562_GWPS_CPM_full_statistics.npz",
    "HCT116_full_statistics.npz",
    "HEK293T_full_statistics.npz",
    "H1_2025_full_statistics.npz",
]
PICK = [0, 137, 299]


def main() -> None:
    targets = pd.read_csv(DATA / "pert_counts.csv").target_gene.astype(str).to_numpy()
    genes = pd.read_csv(DATA / "gene_names.csv").gene_name.astype(str).to_numpy()
    srcs = [sources.load_source(DATA / n) for n in NAMES]
    with np.load(DATA / "CD4_DE_statistics.npz") as d:
        cd4 = {k: d[k] for k in d.files}
    pairs = pd.read_csv(DATA / "official_pairs.csv")
    sub = targets[PICK]
    pred = ad.read_h5ad(C0 / "prediction.h5ad", backed="r")
    rows = []
    for ci, ctx in enumerate("ABC"):
        ctrl = ad.read_h5ad(DATA / f"context_{ctx}.h5ad")
        template, depths, mean, bulk = generator.control_template(
            ctrl.X, 400, generator.POOL_K, generator.SEED + ci
        )
        controls = {"log2fc": mean, "bulk_delta": bulk}
        effects = fusion.fused_effects(srcs, [2, 1, 1, 2], cd4, 0.5, sub, genes, controls)
        p_cpm, p_bulk = generator.expected_moments(effects, controls, sub, genes, pairs=pairs)
        for k, ti in enumerate(PICK):
            ours = generator.dual_moment_counts(
                template,
                p_cpm[k],
                p_bulk[k],
                depths=depths,
                seed=generator.seed_for(f"{ctx}:{targets[ti]}"),
            )
            start = ci * 300 * 400 + ti * 400
            theirs = sparse.csr_matrix(pred.X[start : start + 400]).toarray()
            labels = pred.obs.iloc[start : start + 400]
            rows.append(
                {
                    "context": ctx,
                    "target": str(targets[ti]),
                    "obs_labels_match": bool(
                        (labels.target_gene.astype(str) == targets[ti]).all()
                        and (labels.context.astype(str) == ctx).all()
                    ),
                    "identical_counts": bool(np.array_equal(ours, theirs)),
                    "max_abs_diff": int(np.abs(ours.astype(np.int64) - theirs).max()),
                    "total_counts": int(ours.sum()),
                }
            )
            print(rows[-1], flush=True)
    pred.file.close()
    out = {
        "subset": rows,
        "all_identical": all(r["identical_counts"] and r["obs_labels_match"] for r in rows),
    }
    path = ROOT / "outputs" / "competition_v2" / "c1_license_clean" / "equivalence_c0_bundle.json"
    path.write_text(json.dumps(out, indent=2))
    print("all identical:", out["all_identical"])


if __name__ == "__main__":
    main()
