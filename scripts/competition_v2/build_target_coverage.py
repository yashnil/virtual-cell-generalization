"""Arc 300-target direct-evidence coverage across the competition-v2 public source universe.

Section 9. Rows are the 300 official validation targets; columns are every approved
public perturbation source. A target counts as *directly measured* in a source when
that source applied a perturbation of that exact gene symbol, and as *usable* when it
also passes the source's minimum-evidence rule:

* single-cell atlases: >= 20 perturbed cells (AtlasShift ``minimum_cells``);
* CD4 DE statistics: >= 20 cells **and** the publisher quality flags (AtlasShift rule);
* scPertEval v1 datasets: present in the frozen public label inventory (the v1 rule).

Sources are counted by **cell context** for the tier comparison (K562 essential and
K562 GWPS are one context; H1 train/val/test are one context), because tiers in v1
counted independent contexts, not files.

Inputs are read-only. The frozen ``data/splits/arc_target_support_v1.csv`` is never
modified; this writes ``data/splits/arc_target_support_competition_v2.csv`` and
``outputs/competition_v2/target_coverage_summary.json``.

Reproduce: ``uv run python scripts/competition_v2/build_target_coverage.py``
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
ARC = ROOT / "data" / "raw" / "arc2026" / "controls" / "pert_counts.csv"
V1_TIERS = ROOT / "data" / "splits" / "arc_target_support_v1.csv"
INVENTORY = ROOT / "data" / "provenance" / "scperteval" / "public_label_inventory.json"
PREPARED = ROOT / "outputs" / "competition_v2" / "atlasshift_c0" / "data"
OUT_CSV = ROOT / "data" / "splits" / "arc_target_support_competition_v2.csv"
OUT_JSON = ROOT / "outputs" / "competition_v2" / "target_coverage_summary.json"
MIN_CELLS = 20

#: source column -> (cell context, technology)
V1_SOURCES = {
    "arch1": ("H1 hESC", "CRISPRi Perturb-seq (VCC 2025 train = arch1)"),
    "kaden25rpe1": ("RPE1", "CRISPRi Perturb-seq"),
    "nadig25hepg2": ("HepG2", "CRISPRi Perturb-seq"),
    "nadig25jurkat": ("Jurkat", "CRISPRi Perturb-seq"),
    "replogle22k562": ("K562", "CRISPRi Perturb-seq (essential)"),
    "replogle22rpe1": ("RPE1", "CRISPRi Perturb-seq (essential)"),
    "wessels23": ("K562/other", "CROP-seq variants"),
}
NEW_SOURCES = {
    "K562_GWPS": ("K562", "CRISPRi Perturb-seq (genome-wide)", "K562_GWPS_CPM_full_statistics.npz"),
    "XAtlas_HCT116": (
        "HCT116",
        "CRISPRi Perturb-seq (X-Atlas/Orion)",
        "HCT116_full_statistics.npz",
    ),
    "XAtlas_HEK293T": (
        "HEK293T",
        "CRISPRi Perturb-seq (X-Atlas/Orion)",
        "HEK293T_full_statistics.npz",
    ),
    "H1_2025_full": (
        "H1 hESC",
        "CRISPRi Perturb-seq (VCC 2025 train+val+test)",
        "H1_2025_full_statistics.npz",
    ),
}


def main() -> None:
    targets = pd.read_csv(ARC).target_gene.astype(str).tolist()
    frame = pd.DataFrame({"arc_target": targets})
    v1 = pd.read_csv(V1_TIERS).set_index("arc_target")
    frame["v1_tier"] = frame.arc_target.map(v1.support_tier).astype(int)
    meta = {}

    inventory = json.loads(INVENTORY.read_text())
    for name, (ctx, tech) in V1_SOURCES.items():
        perturbed = set(inventory[name]["perturbations"])
        frame[f"{name}__measured"] = frame.arc_target.isin(perturbed)
        frame[f"{name}__usable"] = frame[f"{name}__measured"]
        frame[f"{name}__n_cells"] = np.nan
        meta[name] = {"context": ctx, "technology": tech}

    for name, (ctx, tech, fname) in NEW_SOURCES.items():
        with np.load(PREPARED / fname, allow_pickle=False) as d:
            n = dict(zip(d["targets"].astype(str), d["n_cells"].astype(int), strict=True))
            splits = (
                dict(zip(d["targets"].astype(str), d["public_2025_split"].astype(str), strict=True))
                if "public_2025_split" in d
                else None
            )
        cells = frame.arc_target.map(lambda t, n=n: n.get(t, 0))
        frame[f"{name}__n_cells"] = cells
        frame[f"{name}__measured"] = cells > 0
        frame[f"{name}__usable"] = cells >= MIN_CELLS
        if splits is not None:
            frame[f"{name}__split"] = frame.arc_target.map(splits).fillna("")
        meta[name] = {"context": ctx, "technology": tech}

    with np.load(PREPARED / "CD4_DE_statistics.npz", allow_pickle=False) as d:
        cd_targets = d["targets"].astype(str)
        avail = d["available"].any(axis=0)
        usable = (d["available"] & d["quality_pass"] & (d["n_cells"] >= MIN_CELLS)).any(axis=0)
        ncells = np.where(d["available"], d["n_cells"], 0).max(axis=0)
        conditions = d["conditions"].astype(str).tolist()
    lookup = {t: i for i, t in enumerate(cd_targets)}
    idx = frame.arc_target.map(lambda t: lookup.get(t, -1)).to_numpy()
    frame["CD4_DE__measured"] = np.where(idx >= 0, avail[np.maximum(idx, 0)], False)
    frame["CD4_DE__usable"] = np.where(idx >= 0, usable[np.maximum(idx, 0)], False)
    frame["CD4_DE__n_cells"] = np.where(idx >= 0, ncells[np.maximum(idx, 0)], 0)
    meta["CD4_DE"] = {
        "context": "primary CD4 T cell",
        "technology": f"CRISPRi Perturb-seq, publisher DE statistics ({', '.join(conditions)})",
    }

    sources = list(meta)
    contexts = sorted({m["context"] for m in meta.values()})
    for ctx in contexts:
        cols = [f"{s}__usable" for s in sources if meta[s]["context"] == ctx]
        frame[f"ctx::{ctx}"] = frame[cols].any(axis=1)
    ctx_cols = [f"ctx::{c}" for c in contexts]
    frame["n_usable_contexts"] = frame[ctx_cols].sum(axis=1)
    frame["n_usable_sources"] = frame[[f"{s}__usable" for s in sources]].sum(axis=1)
    frame["usable_contexts"] = frame[ctx_cols].apply(
        lambda r: "|".join(c.removeprefix("ctx::") for c, v in r.items() if v), axis=1
    )
    atlasshift_sources = ["K562_GWPS", "XAtlas_HCT116", "XAtlas_HEK293T", "H1_2025_full", "CD4_DE"]
    frame["atlasshift_has_signal"] = frame[[f"{s}__usable" for s in atlasshift_sources]].any(axis=1)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(OUT_CSV, index=False)

    def bucket(series: pd.Series) -> dict[str, int]:
        return {
            "0": int((series == 0).sum()),
            "1": int((series == 1).sum()),
            "2": int((series == 2).sum()),
            "3+": int((series >= 3).sum()),
        }

    h1_split = frame.get("H1_2025_full__split", pd.Series("", index=frame.index))
    summary = {
        "n_targets": len(frame),
        "v1_tiers": {str(k): int(v) for k, v in frame.v1_tier.value_counts().sort_index().items()},
        "competition_v2_by_usable_context": bucket(frame.n_usable_contexts),
        "competition_v2_by_usable_source": bucket(frame.n_usable_sources),
        "targets_with_any_direct_evidence_v1": int((frame.v1_tier > 0).sum()),
        "targets_with_any_direct_evidence_v2": int((frame.n_usable_contexts > 0).sum()),
        "targets_with_atlasshift_signal": int(frame.atlasshift_has_signal.sum()),
        "per_source_usable": {s: int(frame[f"{s}__usable"].sum()) for s in sources},
        "per_source_measured": {s: int(frame[f"{s}__measured"].sum()) for s in sources},
        "per_context_usable": {c: int(frame[f"ctx::{c}"].sum()) for c in contexts},
        "h1_2025_full_by_split": {
            k: int(v) for k, v in h1_split[frame["H1_2025_full__usable"]].value_counts().items()
        },
        "h1_2025_gain_over_arch1": int(
            (frame["H1_2025_full__usable"] & ~frame["arch1__usable"]).sum()
        ),
        "v1_tier0_rescued": int(((frame.v1_tier == 0) & (frame.n_usable_contexts > 0)).sum()),
        "v1_tier0_rescued_by_source": {
            s: int(((frame.v1_tier == 0) & frame[f"{s}__usable"]).sum()) for s in sources
        },
        "sources": meta,
        "min_cells": MIN_CELLS,
    }
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
