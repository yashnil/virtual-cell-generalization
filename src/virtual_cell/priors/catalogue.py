"""The registry of prior sources, with provenance and leakage classification.

A prior may be used in the unseen-perturbation benchmark only if it is listed
here and its :attr:`PriorSource.leakage` is acceptable for the regime. The
classification is a claim about the data, argued in the matching file under
``data/provenance/``; it is not a label applied for convenience.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

__all__ = ["Leakage", "PriorSource", "CATALOGUE", "source", "available"]


class Leakage(StrEnum):
    """Whether a source could carry perturbation-response information."""

    #: Contains no perturbation outcome of any kind.
    NONE = "none"
    #: Is itself a perturbation outcome, but of a different readout, assay and
    #: cell panel. Usable only under the restriction named in ``leakage_note``.
    PARTIAL = "partial"
    #: Derived from the perturbation responses under study. Unusable as a
    #: feature in any held-out-perturbation regime.
    DISQUALIFYING = "disqualifying"


@dataclass(frozen=True)
class PriorSource:
    """One audited source of gene-level prior information."""

    name: str
    family: str
    version: str
    licence: str
    provenance: str
    leakage: Leakage
    leakage_note: str
    paths: tuple[Path, ...]

    def missing(self, root: Path) -> list[Path]:
        return [p for p in self.paths if not (root / p).is_file()]


CATALOGUE: dict[str, PriorSource] = {
    "basal": PriorSource(
        name="basal",
        family="target-gene basal expression",
        version="derived from the four-context control cells",
        licence="n/a (derived)",
        provenance="data/provenance/scperteval/",
        leakage=Leakage.NONE,
        leakage_note=(
            "Read from control cells only. A control profile carries no "
            "perturbation response by construction; this is the same class of "
            "input the challenge permits from the Arc controls."
        ),
        paths=(),
    ),
    "hallmark": PriorSource(
        name="hallmark",
        family="pathway membership",
        version="MSigDB h.all.v2024.1.Hs",
        licence="CC BY 4.0 (MSigDB terms)",
        provenance="data/provenance/msigdb/msigdb.md",
        leakage=Leakage.NONE,
        leakage_note=(
            "Curated gene sets published 2024, fixed. No perturbation "
            "measurement from any context in this study."
        ),
        paths=(Path("data/raw/msigdb/h.all.v2024.1.Hs.symbols.gmt"),),
    ),
    "reactome": PriorSource(
        name="reactome",
        family="pathway membership",
        version="MSigDB c2.cp.reactome.v2024.1.Hs",
        licence="CC BY 4.0 (MSigDB terms)",
        provenance="data/provenance/msigdb/msigdb.md",
        leakage=Leakage.NONE,
        leakage_note="As hallmark: curated, fixed, no perturbation outcome.",
        paths=(Path("data/raw/msigdb/c2.cp.reactome.v2024.1.Hs.symbols.gmt"),),
    ),
    "string": PriorSource(
        name="string",
        family="protein association network",
        version="STRING v12.0, organism 9606",
        licence="CC BY 4.0",
        provenance="data/provenance/string/string.md",
        leakage=Leakage.NONE,
        leakage_note=(
            "Aggregate of literature, curated databases, co-occurrence and "
            "co-expression, published 2023 and fixed. Its co-expression "
            "channel is a prior over gene relatedness, not the outcome of a "
            "knockdown."
        ),
        paths=(
            Path("data/raw/string/9606.protein.links.v12.0.txt.gz"),
            Path("data/raw/string/9606.protein.info.v12.0.txt.gz"),
        ),
    ),
    "depmap": PriorSource(
        name="depmap",
        family="CRISPR gene effect",
        version="DepMap 24Q4 Public",
        licence="CC BY 4.0",
        provenance="data/provenance/depmap/depmap.md",
        leakage=Leakage.PARTIAL,
        leakage_note=(
            "Gene effect IS a perturbation outcome -- knockout viability. It "
            "is usable because the readout (fitness, not transcriptome), the "
            "assay (knockout, not CRISPRi) and the cell panel all differ. The "
            "cell panel is the part that cannot be verified, since Arc's "
            "contexts are unidentified and inferring them is prohibited. The "
            "mitigation is binding: only gene-level summaries pooled across "
            "ALL DepMap lines may be used, never a per-line column, so no "
            "context-specific quantity is ever read."
        ),
        paths=(
            Path("data/raw/depmap/CRISPRGeneEffect.csv"),
            Path("data/raw/depmap/Model.csv"),
        ),
    ),
}


def source(name: str) -> PriorSource:
    if name not in CATALOGUE:
        raise KeyError(f"Unknown prior source {name!r}; known: {sorted(CATALOGUE)}")
    return CATALOGUE[name]


def available(root: Path | str = ".") -> dict[str, bool]:
    """Which catalogued sources have all their files present under ``root``."""
    root = Path(root)
    return {name: not src.missing(root) for name, src in CATALOGUE.items()}
