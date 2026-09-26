"""License-clean 300-target direct-evidence coverage (C1, section 4).

For every Arc validation target and every license-GREEN direct-perturbation source:
measured (any cells / any DE row), observations, context, quality and whether C1 may use
it. Usable follows the frozen backbone: >= 20 perturbed cells (single-cell sources), or
>= 1 culture condition that passes every publisher quality flag with >= 20 cells (CD4).

Kaden 2025 RPE1 is license-GREEN but excluded from C1 on scientific grounds (see the
license register), so its columns are recorded and its ``usable`` is always False.

Writes ``data/splits/arc_target_support_license_clean_v1.csv`` and
``outputs/competition_v2/c1_license_clean/coverage_summary.json``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from virtual_cell.competition_v2 import licensing, sources  # noqa: E402

SRC = ROOT / "outputs" / "competition_v2" / "c1_license_clean" / "sources"
ARC = ROOT / "data" / "raw" / "arc2026" / "controls"
V2 = ROOT / "data" / "splits" / "arc_target_support_competition_v2.csv"
OUT = ROOT / "data" / "splits" / "arc_target_support_license_clean_v1.csv"
SUMMARY = ROOT / "outputs" / "competition_v2" / "c1_license_clean" / "coverage_summary.json"


def bins(n: pd.Series) -> dict:
    return {
        "0": int((n == 0).sum()),
        "1": int((n == 1).sum()),
        "2": int((n == 2).sum()),
        "3+": int((n >= 3).sum()),
    }


def main() -> None:
    targets = pd.read_csv(ARC / "pert_counts.csv").target_gene.astype(str).to_numpy()
    table = pd.DataFrame({"arc_target": targets})
    for key, file, context in [
        ("K562_GWPS", "K562_GWPS_CPM_full_statistics.npz", "K562"),
        ("H1_2025_full", "H1_2025_full_statistics.npz", "H1 hESC"),
    ]:
        with np.load(SRC / file) as d:
            n = dict(zip(d["targets"].astype(str), d["n_cells"].astype(int), strict=True))
            split = (
                dict(zip(d["targets"].astype(str), d["public_2025_split"].astype(str), strict=True))
                if "public_2025_split" in d.files
                else {}
            )
        cells = np.array([n.get(t, 0) for t in targets])
        table[f"{key}__measured"] = cells > 0
        table[f"{key}__n_cells"] = cells
        table[f"{key}__context"] = context
        table[f"{key}__quality"] = [
            f"n_cells={c}" + (f"; split={split[t]}" if t in split else "") if c else ""
            for t, c in zip(targets, cells, strict=True)
        ]
        table[f"{key}__usable"] = cells >= sources.MINIMUM_CELLS
    with np.load(SRC / "CD4_DE_statistics.npz") as d:
        idx = {t: i for i, t in enumerate(d["targets"].astype(str))}
        conds = d["conditions"].astype(str)
        avail, qual, ncell = d["available"], d["quality_pass"], d["n_cells"]
    usable_c = avail & qual & (ncell >= sources.MINIMUM_CELLS)
    rows = [idx.get(t, -1) for t in targets]
    table["CD4_DE__measured"] = [r >= 0 and bool(avail[:, r].any()) for r in rows]
    table["CD4_DE__n_cells"] = [int(ncell[:, r].max()) if r >= 0 else 0 for r in rows]
    table["CD4_DE__context"] = "primary CD4 T cell"
    table["CD4_DE__quality"] = [
        ";".join(
            f"{c}:{'pass' if usable_c[k, r] else ('fail' if avail[k, r] else 'absent')}"
            for k, c in enumerate(conds)
        )
        if r >= 0
        else ""
        for r in rows
    ]
    table["CD4_DE__usable"] = [r >= 0 and bool(usable_c[:, r].any()) for r in rows]
    v2 = pd.read_csv(V2).set_index("arc_target")
    table["Kaden_RPE1__measured"] = v2.loc[targets, "kaden25rpe1__measured"].to_numpy()
    table["Kaden_RPE1__n_cells"] = v2.loc[targets, "kaden25rpe1__n_cells"].fillna(0).astype(int)
    table["Kaden_RPE1__context"] = "RPE1"
    table["Kaden_RPE1__quality"] = "excluded: not in C0 backbone; weak reproducibility"
    table["Kaden_RPE1__usable"] = False
    c1_sources = ["K562_GWPS", "H1_2025_full", "CD4_DE"]
    licensing.assert_sources_allowed(["K562", "H1", "CD4"])
    table["n_usable_contexts_c1"] = table[[f"{s}__usable" for s in c1_sources]].sum(axis=1)
    table["usable_contexts_c1"] = [
        "|".join(
            ctx
            for s, ctx in zip(c1_sources, ["K562", "H1 hESC", "primary CD4 T cell"], strict=True)
            if row[f"{s}__usable"]
        )
        for _, row in table.iterrows()
    ]
    table["n_measured_green_incl_kaden"] = table[
        [f"{s}__measured" for s in c1_sources + ["Kaden_RPE1"]]
    ].sum(axis=1)
    table["n_usable_contexts_c0"] = v2.loc[targets, "n_usable_contexts"].to_numpy()
    table["v1_tier"] = v2.loc[targets, "v1_tier"].to_numpy()
    table.to_csv(OUT, index=False)

    lost = table[table.n_usable_contexts_c1 == 0]
    only_x = [t for t in lost.arc_target]
    summary = {
        "n_targets": len(table),
        "v1": {"0": 214, "1": 79, "2": 7, "3+": 0},
        "c0_expanded": bins(table.n_usable_contexts_c0),
        "c1_license_clean": bins(table.n_usable_contexts_c1),
        "c1_any_direct_evidence": int((table.n_usable_contexts_c1 > 0).sum()),
        "c1_percent_coverage": round(100 * float((table.n_usable_contexts_c1 > 0).mean()), 2),
        "c1_measured_any_green_incl_kaden": int((table.n_measured_green_incl_kaden > 0).sum()),
        "per_source_usable": {s: int(table[f"{s}__usable"].sum()) for s in c1_sources},
        "per_source_measured": {
            s: int(table[f"{s}__measured"].sum()) for s in c1_sources + ["Kaden_RPE1"]
        },
        "unsupported_after_removing_xatlas": only_x,
        "unsupported_v1_tier_counts": lost.v1_tier.value_counts().sort_index().to_dict(),
        "lost_contexts_vs_c0_mean": float(
            (table.n_usable_contexts_c0 - table.n_usable_contexts_c1).mean()
        ),
    }
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY.write_text(json.dumps(summary, indent=2, default=int))
    print(
        json.dumps(
            {k: v for k, v in summary.items() if k != "unsupported_after_removing_xatlas"},
            indent=2,
            default=int,
        )
    )
    print("unsupported:", len(only_x), only_x)


if __name__ == "__main__":
    main()
