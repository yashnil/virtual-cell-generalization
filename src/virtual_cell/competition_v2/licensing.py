"""Machine-readable Challenge-use status of every competition-v2 data source.

Mirrors ``reports/competition_v2/data_license_register.md`` (evidence, quotes and URLs
live there). Any C1 input must pass :func:`assert_sources_allowed`. A status changes only
when the register changes with new written evidence. X-Atlas becomes APPROVED only on
explicit written permission from the rights holder
(``reports/competition_v2/xatlas_permission_status.md``).
"""

from __future__ import annotations

GREEN = "GREEN"
BLOCKED = "BLOCKED_PENDING_PERMISSION"
UNKNOWN = "UNKNOWN"

#: Keyed by the short source names the C1 code uses.
STATUS: dict[str, str] = {
    "H1": GREEN,  # VCC 2025 H1: CC0 1.0 (Arc Virtual Cell Atlas page)
    "K562": GREEN,  # Replogle 2022 GWPS: CC BY 4.0 (figshare+ API)
    "CD4": GREEN,  # Marson GWCD4i: "MIT License" (CZI Virtual Cells Platform listing)
    "GENCODE": GREEN,  # GENCODE v47: EMBL-EBI terms / Ensembl "no restrictions"
    "KADEN_RPE1": GREEN,  # CC BY 4.0 (Zenodo) — license-GREEN, scientifically excluded
    "KOLF2.1J_iPSC": GREEN,  # C4: Figshare+ 10.25452/figshare.plus.27261219.v1, CC BY 4.0
    "VIPERTURB_K562": GREEN,  # C4: Zenodo 10.5281/zenodo.18460279, CC BY 4.0 (not downloaded)
    "JURKAT_GSE249595": UNKNOWN,  # C4: GEO, no dataset license; paper CC BY-NC-ND 4.0
    "HCT116": BLOCKED,  # X-Atlas/Orion: CC BY-NC-SA 4.0
    "HEK293T": BLOCKED,  # X-Atlas/Orion: CC BY-NC-SA 4.0
}


class LicenseError(RuntimeError):
    pass


def assert_sources_allowed(names) -> None:
    """Raise unless every named source is GREEN (UNKNOWN and BLOCKED both fail)."""
    bad = {n: STATUS.get(n, UNKNOWN) for n in names if STATUS.get(n, UNKNOWN) != GREEN}
    if bad:
        raise LicenseError(f"non-GREEN sources in a license-clean input: {bad}")
