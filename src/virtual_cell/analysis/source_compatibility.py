"""N5 source-compatibility quantities (reports/n5_protocol.md §2–4).

Every quantity is computed per perturbation from centred response fields. Source
and target are independent experiments, so ⟨S, T⟩ is unbiased for the latent cross
product, while reliable energies come from disjoint split halves.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

EPS = 1e-12


def centre(x: np.ndarray) -> np.ndarray:
    """Remove the mean over perturbations (axis 0): the template."""
    return x - x.mean(axis=0, keepdims=True)


def rowdot(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.einsum("pg,pg->p", a, b)


def rowwise_pearson(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    ac = a - a.mean(axis=1, keepdims=True)
    bc = b - b.mean(axis=1, keepdims=True)
    den = np.linalg.norm(ac, axis=1) * np.linalg.norm(bc, axis=1)
    return np.where(den > EPS, rowdot(ac, bc) / np.maximum(den, EPS), np.nan)


def reliable_energy(halves_a: Sequence[np.ndarray], halves_b: Sequence[np.ndarray]) -> np.ndarray:
    """Per-perturbation E_r⟨ã, b̃⟩ over repeats (each half centred over perturbations)."""
    return np.mean(
        [rowdot(centre(a), centre(b)) for a, b in zip(halves_a, halves_b, strict=True)], axis=0
    )


def split_half_reliability(halves_a, halves_b) -> np.ndarray:
    """Per-perturbation E_r Pearson(ã_p, b̃_p)."""
    return np.mean(
        [rowwise_pearson(centre(a), centre(b)) for a, b in zip(halves_a, halves_b, strict=True)],
        axis=0,
    )


def latent_cosine_terms(S_full, T_full, S_a, S_b, T_a, T_b) -> dict[str, np.ndarray]:
    """Per-perturbation numerator and the two reliable-energy denominators of C_S."""
    return {
        "num": rowdot(centre(S_full), centre(T_full)),
        "den_S": reliable_energy(S_a, S_b),
        "den_T": reliable_energy(T_a, T_b),
    }


def pooled_latent_cosine(num: np.ndarray, den_S: np.ndarray, den_T: np.ndarray) -> float:
    d = float(np.sum(den_S)) * float(np.sum(den_T))
    return float(np.sum(num) / np.sqrt(d)) if d > 0 else float("nan")


def pooled_raw_cosine(S_full: np.ndarray, T_full: np.ndarray) -> dict[str, np.ndarray]:
    Sc, Tc = centre(S_full), centre(T_full)
    return {"num": rowdot(Sc, Tc), "ss": rowdot(Sc, Sc), "tt": rowdot(Tc, Tc)}


def directional_accuracy(
    S_full, T_a: Sequence[np.ndarray], T_b: Sequence[np.ndarray], top: int = 200
):
    """Per perturbation: sign agreement of S̃ with T̃b on the top-|T̃a| genes (mean over repeats)."""
    Sc = centre(S_full)
    vals = []
    for a, b in zip(T_a, T_b, strict=True):
        ac, bc = centre(a), centre(b)
        idx = np.argpartition(-np.abs(ac), top - 1, axis=1)[:, :top]
        s = np.take_along_axis(Sc, idx, axis=1)
        t = np.take_along_axis(bc, idx, axis=1)
        vals.append(np.mean(np.sign(s) == np.sign(t), axis=1))
    return np.mean(vals, axis=0)


def matched_comparator(
    source_reliability: dict[str, np.ndarray], source_cells: dict[str, np.ndarray], reference: str
) -> dict:
    """Predeclared source-only rule: closest median reliability to the reference."""
    ref_rel = float(np.nanmedian(source_reliability[reference]))
    ref_cells = float(np.median(source_cells[reference]))
    rows = []
    for name, rel in source_reliability.items():
        if name == reference:
            continue
        rows.append(
            (
                abs(float(np.nanmedian(rel)) - ref_rel),
                abs(float(np.median(source_cells[name])) - ref_cells),
                name,
            )
        )
    rows.sort()
    return {
        "reference": reference,
        "reference_median_reliability": ref_rel,
        "matched": rows[0][2],
        "candidates": [{"source": n, "abs_rel_gap": a, "abs_cell_gap": b} for a, b, n in rows],
    }


def fe_regression(
    y: np.ndarray,
    pert: np.ndarray,
    source: np.ndarray,
    covariates: np.ndarray,
    sources: Sequence[str],
    reference: str,
) -> dict[str, float]:
    """OLS with perturbation fixed effects; returns β_S − β_reference for each source."""
    others = [s for s in sources if s != reference]
    X = np.column_stack([(source == s).astype(float) for s in others] + [covariates])
    ok = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    y, X, pert = y[ok], X[ok], pert[ok]
    # within-perturbation demeaning
    _, inv = np.unique(pert, return_inverse=True)
    counts = np.bincount(inv)

    def demean(v):
        if v.ndim == 1:
            return v - (np.bincount(inv, weights=v) / counts)[inv]
        return np.column_stack([demean(v[:, j]) for j in range(v.shape[1])])

    coef, *_ = np.linalg.lstsq(demean(X), demean(y), rcond=None)
    return {s: float(coef[i]) for i, s in enumerate(others)}


def cluster_bootstrap_fe(
    y, pert, source, covariates, sources, reference, *, n_boot: int, rng: np.random.Generator
) -> dict[str, tuple[float, float]]:
    """Perturbation-cluster bootstrap of the fixed-effect contrasts."""
    uniq = np.unique(pert)
    rows_of = {p: np.flatnonzero(pert == p) for p in uniq}
    draws = {s: [] for s in sources if s != reference}
    for _ in range(n_boot):
        pick = rng.choice(uniq, len(uniq))
        idx = np.concatenate([rows_of[p] for p in pick])
        new_pert = np.repeat(np.arange(len(pick)), [len(rows_of[p]) for p in pick])
        est = fe_regression(y[idx], new_pert, source[idx], covariates[idx], sources, reference)
        for s, v in est.items():
            draws[s].append(v)
    return {
        s: (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))) for s, v in draws.items()
    }
