"""C3: source-weighted, source-scaled fusion of direct responses (ours, competition track).

Predeclared in ``reports/competition_v2/c3_predeclaration.md``. The object changed is only
the fused mean effect. :func:`combine` reproduces :func:`fusion.fused_effects` **bit for
bit** at equal weights and unit scales (``tests/test_competition_v2_c3.py``), so C1a is
this module's zero point. Weights redistribute mass among the sources that measured a
(target, gene) cell. A cell with one source keeps that source's value, so no weight
choice shrinks a response toward zero.

Also here: sign-consensus classes, source-only reliability (split-half for cell sources,
standard-error based for the CD4 DE table) and the mean-level metrics.
"""

from __future__ import annotations

import numpy as np
from scipy import sparse

from virtual_cell.competition_v2 import fusion
from virtual_cell.competition_v2.atlas import cd4_effect, source_effect
from virtual_cell.competition_v2.sources import SourceStats

CD4 = "CD4"
TOP_K = 200
MIN_RELIABILITY_CELLS = 40
MIN_CPM = 5.0

# ---------------------------------------------------------------------------- fusion


def components(
    srcs: dict[str, SourceStats], cd4: dict | None, targets: np.ndarray, genes: np.ndarray
) -> dict:
    """Per-source panel-aligned effects: ``{space: {name: (effect, mask)}}`` plus CD4.

    Cell sources are looked up exactly as :func:`fusion.fuse_centered` does. CD4 is kept
    in log2fc and carried into ``bulk_delta`` by :func:`fusion.add_cd4` at combine time.
    """
    targets = np.asarray(targets).astype(str)
    out: dict = {"log2fc": {}, "bulk_delta": {}}
    for space in out:
        for name, s in srcs.items():
            full, mask = source_effect(s, s.targets, genes, space=space)
            lookup = {t: i for i, t in enumerate(s.targets)}
            found = np.asarray([lookup.get(t, -1) for t in targets])
            eff = np.zeros((len(targets), len(genes)), dtype=np.float32)
            msk = np.zeros((len(targets), len(genes)), dtype=bool)
            ok = found >= 0
            eff[ok], msk[ok] = full[found[ok]], mask[found[ok]]
            out[space][name] = (eff, msk)
    if cd4 is not None:
        eff, avail = cd4_effect(cd4, targets, genes)
        out[CD4] = (eff, avail)
    return out


def _w(weight, n: int):
    return weight if np.isscalar(weight) else np.asarray(weight, dtype=np.float64)[:, None]


def combine(
    comp: dict,
    space: str,
    control_probability: np.ndarray | None,
    *,
    weights: dict | None = None,
    scales: dict | None = None,
    factor: np.ndarray | None = None,
) -> np.ndarray:
    """Weighted per-cell mean of (scaled) source effects in one space.

    ``weights[name]`` is a scalar or a per-target array; ``scales[name]`` a scalar.
    ``factor`` (targets x genes) multiplies the result (C3d consensus factor).

    A (target, gene) cell that some source measured but whose weighted mass is zero
    (every measuring source has weight 0) takes the equal-weight value, so weighting only
    redistributes among the measuring sources and never removes a response.
    """
    weights = weights or {}
    scales = scales or {}
    parts = comp[space]
    n_t, n_g = next(iter(parts.values()))[0].shape if parts else comp[CD4][0].shape
    num = np.zeros((n_t, n_g))
    den = np.zeros_like(num)
    for name, (eff, msk) in parts.items():
        w = _w(weights.get(name, 1.0), n_t)
        a = scales.get(name, 1.0)
        term = eff if a == 1.0 else (eff.astype(np.float64) * a)
        num += w * term
        den += w * msk
    effect = np.divide(num, den, out=np.zeros_like(num), where=den > 0).astype(np.float32)
    mass = den.astype(np.float32)
    wmass = den
    if CD4 in comp:
        eff, avail = comp[CD4]
        a = scales.get(CD4, 1.0)
        cd4 = eff if a == 1.0 else (eff.astype(np.float64) * a).astype(np.float32)
        w = weights.get(CD4, 1.0)
        if np.isscalar(w):
            effect, wmass = fusion.add_cd4(
                effect, mass, cd4, avail, w, space=space, control_probability=control_probability
            )
        else:
            w = np.asarray(w, dtype=np.float64)[:, None]
            moved = cd4 if space == "log2fc" else fusion.cd4_to_bulk_delta(cd4, control_probability)
            d = mass + w * avail
            wmass = d
            effect = np.divide(
                effect * mass + w * np.where(avail, moved, 0),
                d,
                out=np.zeros(d.shape),
                where=d > 0,
            )
    if weights:
        lost = wmass <= 0
        if lost.any():
            equal = combine(comp, space, control_probability, scales=scales)
            effect = np.where(lost, equal, effect)
    if factor is not None:
        effect = effect * factor
    return np.asarray(effect)


def fused(comp, controls, **kw) -> dict[str, np.ndarray]:
    """Both spaces, as :func:`fusion.fused_effects` returns them."""
    return {
        space: combine(comp, space, controls.get(space) if controls else None, **kw)
        for space in ("log2fc", "bulk_delta")
    }


def source_list(comp) -> list[tuple[str, np.ndarray, np.ndarray]]:
    """``(name, log2fc effect, availability)`` for every source, CD4 included."""
    out = [(n, e, m) for n, (e, m) in comp["log2fc"].items()]
    if CD4 in comp:
        out.append((CD4, *comp[CD4]))
    return out


def consensus_classes(comp) -> np.ndarray:
    """Per (target, gene): 0 none, 1 single source, 2 agree, 3 conflict (log2fc signs)."""
    srcs = source_list(comp)
    signs = np.stack([np.where(m, np.sign(e), 0) for _, e, m in srcs])
    n = (signs != 0).sum(axis=0)
    pos = (signs > 0).sum(axis=0)
    neg = (signs < 0).sum(axis=0)
    out = np.zeros(n.shape, dtype=np.int8)
    out[n == 1] = 1
    out[(n >= 2) & ((pos == n) | (neg == n))] = 2
    out[(n >= 2) & (pos > 0) & (neg > 0)] = 3
    return out


def consensus_factor(classes: np.ndarray, conflict: float) -> np.ndarray:
    return np.where(classes == 3, conflict, 1.0).astype(np.float32)


# ---------------------------------------------------------------------------- reliability


def split_half(
    cells: sparse.csr_matrix,
    labels: np.ndarray,
    targets,
    control_label: str,
    *,
    seed: int,
    min_cells: int = MIN_RELIABILITY_CELLS,
) -> dict:
    """Spearman-Brown split-half reliability and noise energy per target, from raw cells.

    Each half's log2fc (mean per-cell CPM vs the control mean, pseudocount 1) is centred
    over targets within that half. Genes: control mean CPM ≥ 5. Returns
    ``{"genes": mask, "rel": {t: r_sb}, "noise": {t: ‖A−B‖²/4}, "energy": {t: ‖(A+B)/2‖²}}``.
    """
    cells = sparse.csr_matrix(cells, dtype=np.float64)
    lib = np.asarray(cells.sum(axis=1)).ravel()
    cpm = sparse.diags(1e6 / np.where(lib > 0, lib, 1)) @ cells
    ctrl = np.asarray(cpm[labels == control_label].mean(axis=0)).ravel()
    keep = ctrl >= MIN_CPM
    rng = np.random.default_rng(seed)
    kept, halves = [], []
    for t in targets:
        rows = np.flatnonzero(labels == t)
        if len(rows) < min_cells:
            continue
        rows = rng.permutation(rows)
        h = len(rows) // 2
        pair = [
            np.asarray(cpm[rows[:h]].mean(axis=0)).ravel()[keep],
            np.asarray(cpm[rows[h:]].mean(axis=0)).ravel()[keep],
        ]
        halves.append([np.log2((m + 1) / (ctrl[keep] + 1)) for m in pair])
        kept.append(str(t))
    a = np.stack([x[0] for x in halves])
    b = np.stack([x[1] for x in halves])
    a -= a.mean(axis=0, keepdims=True)
    b -= b.mean(axis=0, keepdims=True)
    rel, noise, energy = {}, {}, {}
    for i, t in enumerate(kept):
        r = np.corrcoef(a[i], b[i])[0, 1]
        rel[t] = float(max(0.0, 2 * r / (1 + r))) if np.isfinite(r) and r > -1 else 0.0
        noise[t] = float(((a[i] - b[i]) ** 2).sum() / 4)
        energy[t] = float((((a[i] + b[i]) / 2) ** 2).sum())
    return {"genes": keep, "rel": rel, "noise": noise, "energy": energy}


def cd4_reliability(stats: dict, minimum_cells: int = 20) -> dict:
    """Per target: mean over usable conditions of ``1 − mean(SE²)/var(lfc)``, in [0, 1]."""
    usable = stats["available"] & stats["quality_pass"] & (stats["n_cells"] >= minimum_cells)
    out = {}
    for k, t in enumerate(np.asarray(stats["targets"]).astype(str)):
        vals = []
        for c in range(stats["log2fc"].shape[0]):
            if not usable[c, k]:
                continue
            lfc = stats["log2fc"][c, k].astype(float)
            se = stats["lfcSE"][c, k].astype(float)
            ok = np.isfinite(lfc) & np.isfinite(se)
            if ok.sum() < 10 or lfc[ok].var() <= 0:
                continue
            vals.append(float(np.clip(1 - (se[ok] ** 2).mean() / lfc[ok].var(), 0, 1)))
        if vals:
            out[t] = float(np.mean(vals))
    return out


# ---------------------------------------------------------------------------- mean metrics


def row_metrics(pred: np.ndarray, truth: np.ndarray, valid: np.ndarray, top_k: int = TOP_K):
    """Per-target cosine, Pearson, norm ratio, top-k sign accuracy, and sums for energy."""
    n = len(pred)
    out = {k: np.full(n, np.nan) for k in ("cosine", "pearson", "norm_ratio", "sign_acc")}
    sse = sst = 0.0
    for i in range(n):
        v = valid[i] & np.isfinite(truth[i])
        p, t = pred[i][v].astype(float), truth[i][v].astype(float)
        if not v.any() or not np.any(p) or not np.any(t):
            continue
        np_, nt = np.linalg.norm(p), np.linalg.norm(t)
        out["cosine"][i] = p @ t / (np_ * nt)
        if p.std() > 0 and t.std() > 0:
            out["pearson"][i] = np.corrcoef(p, t)[0, 1]
        out["norm_ratio"][i] = np_ / nt
        top = np.argsort(-np.abs(t))[:top_k]
        top = top[p[top] != 0]
        if len(top):
            out["sign_acc"][i] = np.mean(np.sign(p[top]) == np.sign(t[top]))
        sse += ((t - p) ** 2).sum()
        sst += (t**2).sum()
    return out, sse, sst


def summarise(per: dict, sse: float, sst: float) -> dict:
    return {
        "cosine": float(np.nanmean(per["cosine"])),
        "pearson": float(np.nanmean(per["pearson"])),
        "norm_ratio_median": float(np.nanmedian(per["norm_ratio"])),
        "sign_acc": float(np.nanmean(per["sign_acc"])),
        "energy_explained": float(1 - sse / sst) if sst > 0 else np.nan,
        "n_targets": int(np.isfinite(per["cosine"]).sum()),
    }
