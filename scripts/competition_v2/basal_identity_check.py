"""Is any hidden Arc context basally one of the public source cell lines? (section 7)

Correlates the mean control expression (log1p CPM, gene panel intersection, genes above
1 CPM in either profile) of each official control context A/B/C with the control
profile of every prepared AtlasShift source. The A-B-C cross-correlations are the
yardstick: a source that matched a context far better than the contexts match each
other would mean that context is (or is very close to) that cell line. That would be a
strategic fact, not leakage: the source is an independent public experiment.

Reproduce: ``uv run python scripts/competition_v2/basal_identity_check.py``
"""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[2]
CONTROLS = ROOT / "data" / "raw" / "arc2026" / "controls"
PREPARED = ROOT / "outputs" / "competition_v2" / "atlasshift_c0" / "data"
OUT = ROOT / "outputs" / "competition_v2" / "basal_identity.json"
SOURCES = {
    "K562_GWPS": "K562_GWPS_CPM_full_statistics.npz",
    "HCT116": "HCT116_full_statistics.npz",
    "HEK293T": "HEK293T_full_statistics.npz",
    "H1_2025": "H1_2025_full_statistics.npz",
}


def main() -> None:
    genes = pd.read_csv(CONTROLS / "gene_names.csv").gene_name.astype(str).to_numpy()
    profiles = {}
    for c in "ABC":
        x = sparse.csr_matrix(ad.read_h5ad(CONTROLS / f"context_{c}.h5ad").X)
        lib = np.asarray(x.sum(axis=1)).ravel()
        cpm = np.asarray(x.multiply(1e6 / lib[:, None]).mean(axis=0)).ravel()
        profiles[c] = pd.Series(cpm, index=genes)
    for name, fname in SOURCES.items():
        path = PREPARED / fname
        if not path.exists():
            continue
        with np.load(path, allow_pickle=False) as d:
            s = pd.Series(d["global_control_mean_cpm"].astype(float), index=d["genes"].astype(str))
        profiles[name] = s[~s.index.duplicated()]
    names = list(profiles)
    table = pd.DataFrame(index=names, columns=names, dtype=float)
    for a in names:
        for b in names:
            common = profiles[a].index.intersection(profiles[b].index).intersection(genes)
            pa, pb = profiles[a][common], profiles[b][common]
            # renormalise to the shared panel so gene-panel filtering cannot bias it
            pa, pb = 1e6 * pa / pa.sum(), 1e6 * pb / pb.sum()
            keep = (pa > 1) | (pb > 1)
            table.loc[a, b] = np.corrcoef(np.log1p(pa[keep]), np.log1p(pb[keep]))[0, 1]
    OUT.write_text(json.dumps(table.round(4).to_dict(), indent=2))
    print(table.round(3).to_string())


if __name__ == "__main__":
    main()
