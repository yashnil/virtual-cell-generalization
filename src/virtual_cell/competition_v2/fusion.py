"""Source fusion, source agreement, and agreement-controlled shrinkage.

Fusion attribution: AtlasShift ``fuse_source_centered`` / ``add_cd4_family`` (@ ``d24ce4f``,
MIT). A source that did not measure a (target, gene) contributes to neither numerator
nor denominator, so fusion is a per-cell weighted mean over the sources that measured it.

* **C0 weights** (upstream): K562 2, HCT116 1, HEK293T 1, H1 2, CD4 0.5.
* **C1a** (ours): every usable GREEN source weighs 1 — an equal per-cell mean.

Agreement and shrinkage are ours (competition track, C1b); both are predeclared in
``reports/competition_v2/c1_predeclaration.md`` before any C1b arm was scored.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np
from scipy import stats as sstats

from virtual_cell.competition_v2.atlas import cd4_effect, cd4_to_bulk_delta, source_effect
from virtual_cell.competition_v2.sources import MINIMUM_CELLS, SourceStats

#: Upstream C0 fusion weights, by source name as stored in the statistics files.
C0_WEIGHTS = {
    "K562_GWPS_CPM": 2.0,
    "HCT116": 1.0,
    "HEK293T": 1.0,
    "H1_2025_public": 2.0,
    "CD4": 0.5,
}


def fuse_centered(
    sources: list[SourceStats],
    weights: list[float],
    targets: np.ndarray,
    genes: np.ndarray,
    *,
    space: str,
    center: bool = True,
    minimum_cells: int = MINIMUM_CELLS,
    pseudocount_cpm: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Weighted per-cell mean of source-centred effects -> ``(effect, weight mass)``.

    Each source is centred over **its own** retained targets (not over the panel), then
    looked up for the panel targets.
    """
    if len(sources) != len(weights) or any(w < 0 for w in weights):
        raise ValueError("invalid source weights")
    targets = np.asarray(targets).astype(str)
    num = np.zeros((len(targets), len(genes)))
    den = np.zeros_like(num)
    for source, weight in zip(sources, weights, strict=True):
        full, mask = source_effect(
            source,
            source.targets,
            genes,
            space=space,
            minimum_cells=minimum_cells,
            pseudocount_cpm=pseudocount_cpm,
            center=center,
        )
        lookup = {t: i for i, t in enumerate(source.targets)}
        found = np.asarray([lookup.get(t, -1) for t in targets])
        rows = np.flatnonzero(found >= 0)
        num[rows] += weight * full[found[rows]]
        den[rows] += weight * mask[found[rows]]
    effect = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    return effect.astype(np.float32), den.astype(np.float32)


def add_cd4(
    effect: np.ndarray,
    mass: np.ndarray,
    cd4: np.ndarray,
    available: np.ndarray,
    weight: float,
    *,
    space: str,
    control_probability: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Fold the CD4 DE family into a fused effect as one more weighted source."""
    if weight < 0:
        raise ValueError("CD4 weight must be nonnegative")
    if space == "log2fc":
        moved = cd4
    elif space == "bulk_delta":
        moved = cd4_to_bulk_delta(cd4, control_probability)
    else:
        raise ValueError(space)
    den = mass + weight * available
    out = np.divide(
        effect * mass + weight * np.where(available, moved, 0),
        den,
        out=np.zeros_like(effect),
        where=den > 0,
    )
    return out, den


# ---------------------------------------------------------------------------
# source agreement (C1b) — predeclared statistic
# ---------------------------------------------------------------------------


def panel_source_vectors(
    sources: list[SourceStats],
    cd4_stats: dict | None,
    targets: np.ndarray,
    genes: np.ndarray,
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Per-source centred ``log2fc`` effects on the panel: ``{name: (effect, mask)}``.

    ``log2fc`` is the only space every GREEN source (including the CD4 DE table) has.
    Centering is the frozen per-source definition (over the source's own targets).
    """
    targets = np.asarray(targets).astype(str)
    out = {}
    for source in sources:
        full, mask = source_effect(source, source.targets, genes, space="log2fc")
        lookup = {t: i for i, t in enumerate(source.targets)}
        found = np.asarray([lookup.get(t, -1) for t in targets])
        eff = np.zeros((len(targets), len(genes)), dtype=np.float32)
        msk = np.zeros((len(targets), len(genes)), dtype=bool)
        ok = found >= 0
        eff[ok], msk[ok] = full[found[ok]], mask[found[ok]]
        out[source.name] = (eff, msk)
    if cd4_stats is not None:
        eff, avail = cd4_effect(cd4_stats, targets, genes)
        out["CD4"] = (eff, avail)
    return out


def source_agreement(
    vectors: dict[str, tuple[np.ndarray, np.ndarray]],
    exclude_genes: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Mean pairwise cosine between the usable sources' centred effects, per target.

    For each pair of sources that both measured target ``p``, the cosine is taken over the
    genes both measured, minus ``exclude_genes`` (the panel target genes, so the shared
    knockdown of the target itself cannot manufacture agreement). A pair where either
    vector has zero norm on that set is skipped. Targets with fewer than two usable
    sources get ``NaN``. Returns ``(agreement, n_sources)``.
    """
    names = list(vectors)
    keep = ~np.asarray(exclude_genes, dtype=bool)
    n_t = next(iter(vectors.values()))[0].shape[0]
    covered = np.stack([vectors[n][1][:, keep].any(axis=1) for n in names])
    n_sources = covered.sum(axis=0)
    total = np.zeros(n_t)
    pairs = np.zeros(n_t)
    for a, b in combinations(range(len(names)), 2):
        ea, ma = vectors[names[a]]
        eb, mb = vectors[names[b]]
        both = ma[:, keep] & mb[:, keep]
        va = np.where(both, ea[:, keep], 0).astype(np.float64)
        vb = np.where(both, eb[:, keep], 0).astype(np.float64)
        na, nb = np.linalg.norm(va, axis=1), np.linalg.norm(vb, axis=1)
        ok = covered[a] & covered[b] & (na > 0) & (nb > 0)
        cos = np.zeros(n_t)
        cos[ok] = (va[ok] * vb[ok]).sum(axis=1) / (na[ok] * nb[ok])
        total += cos
        pairs += ok
    agreement = np.divide(total, pairs, out=np.full(n_t, np.nan), where=pairs > 0)
    return agreement, n_sources


def agreement_percentile(agreement: np.ndarray) -> np.ndarray:
    """Empirical CDF position of each target's agreement among targets that have one.

    Midrank ``(rank - 0.5) / n``; targets without an agreement get 0.5 (uninformative).
    Uses only the prediction panel's own source statistics — no held-out response.
    """
    out = np.full(len(agreement), 0.5)
    ok = np.isfinite(agreement)
    if ok.sum():
        out[ok] = (sstats.rankdata(agreement[ok]) - 0.5) / ok.sum()
    return out


def shrinkage_factor(
    rule: str, n_targets: int, *, scalar: float = 1.0, floor: float = 0.0, agreement=None
) -> np.ndarray:
    """Per-target multiplier ``lambda_p`` for the predeclared family.

    * ``S0`` — ``1``;
    * ``S1`` — one global ``scalar``;
    * ``S2`` — ``floor + (1 - floor) * F(agreement_p)``, ``F`` = :func:`agreement_percentile`.
    """
    if rule == "S0":
        return np.ones(n_targets)
    if rule == "S1":
        if not 0 <= scalar <= 1:
            raise ValueError("S1 scalar must lie in [0, 1]")
        return np.full(n_targets, float(scalar))
    if rule == "S2":
        if agreement is None or not 0 <= floor <= 1:
            raise ValueError("S2 needs agreement and a floor in [0, 1]")
        return floor + (1 - floor) * agreement_percentile(np.asarray(agreement))
    raise ValueError(rule)


def fused_effects(
    sources: list[SourceStats],
    weights: list[float],
    cd4_stats: dict | None,
    cd4_weight: float,
    targets: np.ndarray,
    genes: np.ndarray,
    controls: dict[str, np.ndarray],
) -> dict[str, np.ndarray]:
    """The backbone's fused effect in both spaces, CD4 folded in as one more source.

    ``controls`` maps each space to the destination control composition the CD4 log2FC
    is carried through (``log2fc``: mean-CPM composition; ``bulk_delta``: pseudobulk).
    """
    if cd4_stats is not None:
        cd4, available = cd4_effect(cd4_stats, targets, genes)
    out = {}
    for space in ("log2fc", "bulk_delta"):
        if sources:
            effect, mass = fuse_centered(sources, weights, targets, genes, space=space)
        else:
            effect = np.zeros((len(targets), len(genes)), dtype=np.float32)
            mass = np.zeros_like(effect)
        if cd4_stats is not None:
            effect, mass = add_cd4(
                effect,
                mass,
                cd4,
                available,
                cd4_weight,
                space=space,
                control_probability=controls[space],
            )
        out[space] = effect
    return out
