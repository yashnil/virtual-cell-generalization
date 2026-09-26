"""Competition track, v2: expanded public perturbation evidence for the Arc 2026 challenge.

Logically separate from the frozen research modules. Three provenance classes are
kept apart and must never be merged silently:

* **our method** — the frozen v1 model (``virtual_cell.arc``, ``virtual_cell.modelling``);
* **external public baselines** — AtlasShift, vendored verbatim under
  ``third_party/atlasshift`` and only ever *run*, never imported here;
* **later hybrids** — anything combining the two, which must name both parents.

Report: ``reports/competition_v2/competitive_baseline_expansion_v1.md``.
"""

from __future__ import annotations

ATLASSHIFT_REPO = "https://github.com/kaipengm2/Virtual-Cell-Challenge-2026"
ATLASSHIFT_COMMIT = "d24ce4fdae0cd7cba3ba29546cd8737094eae98a"
XATLAS_REPO = "slaf-project/X-Atlas-Orion"
XATLAS_REVISION = "598aa5442a4ee9a8f1a383a843ab5ce142904d04"

#: Provenance label of every candidate the competition track can emit.
PROVENANCE = {
    "V1": "ours: arc_count_space_baseline_v1 (frozen, submitted 2026-09-26)",
    "C0_ATLASSHIFT_REPRODUCTION": f"external: AtlasShift @ {ATLASSHIFT_COMMIT[:12]} (MIT)",
}
