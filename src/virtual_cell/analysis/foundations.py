r"""Transferability foundations: template recovery, pathway scores, candidate D targets.

Low-capacity, leakage-audited machinery for the foundations study. Nothing here
is a learned model beyond a *single global scalar* or a heavily regularised map
on a rank-<=2 basal subspace, because four contexts (three sources per fold)
cannot support more.

Template algebra
----------------
With ``delta[c,p] = mu + alpha_c + beta_p + gamma[c,p]``, ``sum_c alpha_c = 0``
and ``sum_c gamma[c,p] = 0``, the source mean over the three sources ``S`` is

``A[p] = mu - alpha_{c*}/3 + beta_p - gamma[c*,p]/3``

Splitting it into a context-level template and a perturbation-level part:

* ``source_template = mean_p A[p] = mu - alpha_{c*}/3``
* ``A_centred[p] = A[p] - source_template = beta_p - gamma[c*,p]/3``

The *true* target template is ``T = mu + alpha_{c*} = mean_p delta[c*,p]``, so

.. math::
    T - \text{source\_template} = \tfrac{4}{3}\,\alpha_{c*}

Hence a template-corrected prediction is
``A[p] + (4/3) * alpha_hat``, and the ``alpha_hat = 0`` case is exactly the
uncorrected source mean. This factor of 4/3 is easy to miss and is asserted in
the tests.

Everything an estimator may see
-------------------------------
Control (basal) profiles of **all** contexts including the held-out one, and
perturbation responses of the **source** contexts only. Never the target's
perturbation responses.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

EPS = 1e-12


# --------------------------------------------------------------------------
# A. basal template / scale recovery
# --------------------------------------------------------------------------


def basal_deviation(control_means: np.ndarray, index: int, reference: Sequence[int]) -> np.ndarray:
    """Basal profile of ``index`` minus the mean basal profile of ``reference``.

    Uses control cells only, so it is computable for the held-out context at
    inference time.
    """
    b = np.asarray(control_means, dtype=np.float64)
    return b[index] - b[list(reference)].mean(axis=0)


def source_alphas(delta: np.ndarray, sources: Sequence[int]) -> np.ndarray:
    """Within-source context templates ``alpha_s`` from the source set alone.

    ``alpha_s = mean_p delta[s,p] - mean_{s' in S} mean_p delta[s',p]``, i.e. the
    three-context decomposition's alpha. No target response is touched.
    """
    d = np.asarray(delta, dtype=np.float64)[list(sources)]
    per_context = d.mean(axis=1)
    return per_context - per_context.mean(axis=0, keepdims=True)


def fit_global_scale(alphas: np.ndarray, deviations: np.ndarray) -> float:
    """One global scalar ``k`` minimising ``sum_s ||alpha_s - k dev_s||^2``.

    A single degree of freedom — the lowest-capacity non-trivial estimator that
    can map a basal deviation onto a response template.
    """
    a = np.asarray(alphas, dtype=np.float64).ravel()
    d = np.asarray(deviations, dtype=np.float64).ravel()
    denom = float(d @ d)
    return float(a @ d / denom) if denom > EPS else 0.0


def ridge_subspace_map(
    alphas: np.ndarray,
    deviations: np.ndarray,
    target_deviation: np.ndarray,
    *,
    n_components: int = 2,
    ridge: float = 1.0,
) -> np.ndarray:
    """Project basal deviations onto their own low-rank subspace, then ridge-map.

    With three sources the deviation matrix has rank <= 2, so this has at most
    two free parameters per output dimension before regularisation, and the ridge
    penalty shrinks it further. Fitted on sources only.
    """
    X = np.asarray(deviations, dtype=np.float64)
    Y = np.asarray(alphas, dtype=np.float64)
    centre = X.mean(axis=0)
    Xc = X - centre
    # right singular vectors span the deviation subspace
    _, _, vt = np.linalg.svd(Xc, full_matrices=False)
    k = min(n_components, vt.shape[0])
    basis = vt[:k]  # (k, G)
    coords = Xc @ basis.T  # (n_sources, k)
    gram = coords.T @ coords + ridge * np.eye(k)
    weights = np.linalg.solve(gram, coords.T @ Y)  # (k, G)
    target_coords = (np.asarray(target_deviation, dtype=np.float64) - centre) @ basis.T
    return target_coords @ weights


def fit_response_shrinkage(delta: np.ndarray, sources: Sequence[int]) -> float:
    r"""Optimal scalar shrinkage of the perturbation-specific prediction.

    Source-mean transfer predicts ``beta_p`` but cannot predict ``gamma``, so the
    *variance-optimal* point prediction is a shrunk version of it. The factor is
    estimated by an **inner leave-one-source-out** loop: each source is held out
    in turn, its centred response is regressed on the centred mean of the
    remaining sources, and the pooled least-squares scalar is returned.

    Only source contexts are ever touched, so this is available at inference.

    Caveat, stated because it biases the estimate: an inner fold has two
    remaining sources rather than three, so its gamma attenuation is ``-1/2``
    rather than ``-1/3``. The resulting shrinkage is therefore slightly
    conservative (closer to 1) than the outer-fold optimum.
    """
    d = np.asarray(delta, dtype=np.float64)
    num = 0.0
    den = 0.0
    src = list(sources)
    for held in src:
        others = [s for s in src if s != held]
        inner = d[others].mean(axis=0)
        inner_c = inner - inner.mean(axis=0, keepdims=True)
        y = d[held]
        y_c = y - y.mean(axis=0, keepdims=True)
        num += float(np.sum(y_c * inner_c))
        den += float(np.sum(inner_c * inner_c))
    return float(num / den) if den > EPS else 1.0


def apply_scale_correction(
    source_mean_pred: np.ndarray, scale: float, alpha_hat: np.ndarray | None = None
) -> np.ndarray:
    """Shrink the perturbation-specific part, optionally correcting the template.

    The source template ``mean_p A[p]`` is kept (and shifted by ``(4/3) alpha_hat``
    if given) while only the centred, perturbation-specific component is scaled.
    """
    A = np.asarray(source_mean_pred, dtype=np.float64)
    template = A.mean(axis=0, keepdims=True)
    centred = A - template
    if alpha_hat is not None:
        template = template + (4.0 / 3.0) * np.asarray(alpha_hat, dtype=np.float64)[None, :]
    return template + scale * centred


@dataclass(frozen=True)
class TemplateEstimate:
    """One estimate of the held-out context template correction."""

    name: str
    alpha_hat: np.ndarray
    detail: Mapping[str, float]


def template_estimators(
    delta: np.ndarray,
    control_means: np.ndarray,
    target: int,
    sources: Sequence[int],
    *,
    shrink: float = 1.0,
) -> list[TemplateEstimate]:
    """All permitted low-capacity estimators of ``alpha`` for the held-out context.

    ``zero`` reproduces the uncorrected source mean. ``direct`` uses the raw basal
    deviation with no scaling. ``global_scalar`` fits one coefficient on the
    sources. ``ridge_subspace`` fits a rank-<=2 regularised map on the sources.
    """
    all_contexts = sorted([target, *sources])
    dev_target = basal_deviation(control_means, target, all_contexts)
    alphas_S = source_alphas(delta, sources)
    devs_S = np.vstack([basal_deviation(control_means, s, sources) for s in sources])

    k = fit_global_scale(alphas_S, devs_S)
    out = [
        TemplateEstimate("zero", np.zeros_like(dev_target), {}),
        TemplateEstimate("direct_basal", dev_target, {}),
        TemplateEstimate("global_scalar", k * shrink * dev_target, {"k": k, "shrink": shrink}),
        TemplateEstimate(
            "ridge_subspace",
            ridge_subspace_map(alphas_S, devs_S, dev_target),
            {},
        ),
    ]
    return out


def apply_template_correction(source_mean_pred: np.ndarray, alpha_hat: np.ndarray) -> np.ndarray:
    """``A[p] + (4/3) alpha_hat`` — see the module docstring for the 4/3."""
    return np.asarray(source_mean_pred, dtype=np.float64) + (4.0 / 3.0) * np.asarray(
        alpha_hat, dtype=np.float64
    )


def oracle_alpha_evaluation_only(delta: np.ndarray, target: int) -> np.ndarray:
    """True ``alpha_{c*}`` from the four-context decomposition.

    **Evaluation-only.** Uses the held-out responses and is a ceiling, never a
    usable estimator.
    """
    d = np.asarray(delta, dtype=np.float64)
    per_context = d.mean(axis=1)
    return per_context[target] - per_context.mean(axis=0)


# --------------------------------------------------------------------------
# C. pathway aggregation
# --------------------------------------------------------------------------


def read_gmt(path: str | Path) -> dict[str, list[str]]:
    """Parse a GMT file into ``{set_name: [genes]}``, exactly as released."""
    sets: dict[str, list[str]] = {}
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 3:
            continue
        sets[fields[0]] = [g for g in fields[2:] if g]
    return sets


def pathway_membership(
    gene_sets: Mapping[str, Sequence[str]], genes: Sequence[str], *, min_genes: int = 10
) -> tuple[list[str], np.ndarray]:
    """Binary membership matrix over the frozen gene space.

    Returns ``(names, M)`` with ``M`` of shape ``(n_sets, n_genes)``. Sets with
    fewer than ``min_genes`` present are dropped — a predeclared property of our
    gene space, not a tuning knob.
    """
    index = {g: i for i, g in enumerate(genes)}
    names, rows = [], []
    for name, members in gene_sets.items():
        hits = [index[g] for g in members if g in index]
        if len(hits) < min_genes:
            continue
        row = np.zeros(len(genes))
        row[hits] = 1.0
        names.append(name)
        rows.append(row)
    if not rows:
        raise ValueError("No gene set met the minimum overlap.")
    return names, np.vstack(rows)


def pathway_scores(responses: np.ndarray, membership: np.ndarray) -> np.ndarray:
    """Mean response over each set's present genes.

    Predeclared aggregation: the unweighted mean. Deliberately the simplest
    choice; no weighting, ranking or enrichment statistic is used, and the
    definition was fixed before any pathway result was computed.
    """
    R = np.asarray(responses, dtype=np.float64)
    M = np.asarray(membership, dtype=np.float64)
    counts = M.sum(axis=1)
    if np.any(counts <= 0):
        raise ValueError("A gene set has no member genes.")
    flat = R.reshape(-1, R.shape[-1])
    scores = (flat @ M.T) / counts[None, :]
    return scores.reshape(*R.shape[:-1], M.shape[0])


# --------------------------------------------------------------------------
# D. candidate transferability targets
# --------------------------------------------------------------------------


def reliable_energy(half_a: np.ndarray, half_b: np.ndarray) -> np.ndarray:
    r"""Unbiased estimate of ``||latent||^2`` from two independent halves.

    With ``h_i = L + e_i`` and independent, zero-mean ``e_i``,
    ``E<h1, h2> = ||L||^2`` because the noise cross-term vanishes. Unlike
    ``||h||^2`` this is **not** inflated by measurement noise, which is exactly
    why a raw response norm is a biased transferability target.
    """
    a = np.asarray(half_a, dtype=np.float64)
    b = np.asarray(half_b, dtype=np.float64)
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch {a.shape} vs {b.shape}")
    return np.einsum("...g,...g->...", a, b)


def reliable_residual_energy(
    half_a: np.ndarray, half_b: np.ndarray, prediction: np.ndarray
) -> np.ndarray:
    r"""Unbiased estimate of ``||latent - prediction||^2``.

    ``E<h1 - A, h2 - A> = ||L - A||^2`` for a prediction ``A`` independent of the
    target noise — which source-only predictions are by construction.
    """
    a = np.asarray(half_a, dtype=np.float64)
    b = np.asarray(half_b, dtype=np.float64)
    p = np.asarray(prediction, dtype=np.float64)
    return np.einsum("...g,...g->...", a - p, b - p)


def pooled_unexplained_fraction(
    half_a: np.ndarray, half_b: np.ndarray, prediction: np.ndarray
) -> float:
    r"""Ratio of summed reliable residual energy to summed reliable energy.

    ``sum_p <h1-A, h2-A> / sum_p <h1, h2>`` rather than the mean of per-pair
    ratios.

    **This is the numerically safe aggregate.** The per-pair ratio in
    :func:`candidate_targets` divides by an estimate of ``||L_p||^2`` that is
    itself noisy and can approach or cross zero for unreliable perturbations, so
    its mean is unstable and can even go negative at high noise. Pooling moves
    the noisy quantity into a sum before the division, which is stable. Use the
    per-pair value for ranking individual perturbations (with a reliability
    filter) and the pooled value for any headline number.
    """
    rel_y = reliable_energy(half_a, half_b)
    rel_r = reliable_residual_energy(half_a, half_b, prediction)
    total = float(np.sum(rel_y))
    return float(np.sum(rel_r) / total) if abs(total) > EPS else np.nan


def candidate_targets(
    half_a: np.ndarray, half_b: np.ndarray, prediction: np.ndarray, *, min_energy: float = 0.0
) -> pd.DataFrame:
    """Every candidate ``D`` definition, side by side.

    * ``D_unexplained_fraction`` -- reliable residual energy over reliable target
      energy. 0 means transfer explains all reproducible signal, 1 means it
      explains none, and **>1 means the prediction is worse than predicting
      zero**. Undefined when reliable energy is non-positive.
    * ``D_explained_fraction`` -- its complement, the "success" orientation.
    * ``D_raw_residual_fraction`` -- the same ratio computed on raw (noisy)
      quantities. Included as the *naive* comparator: it is biased upward by
      measurement noise.
    * ``D_raw_gamma_norm_proxy`` -- the raw residual norm. Included to be
      **rejected**: it is depth-biased and conflates noise with biology.

    Edge case: the per-pair fraction divides by a *noisy* estimate of
    ``||L_p||^2``. For unreliable perturbations that denominator approaches zero
    and the ratio becomes unstable (and can be negative), so per-pair values need
    a reliability filter. :func:`pooled_unexplained_fraction` is the stable
    aggregate.
    """
    rel_y = reliable_energy(half_a, half_b)
    rel_r = reliable_residual_energy(half_a, half_b, prediction)
    pooled = (np.asarray(half_a, dtype=np.float64) + np.asarray(half_b, dtype=np.float64)) / 2.0
    resid = pooled - np.asarray(prediction, dtype=np.float64)
    usable = rel_y > min_energy
    frac = np.full(rel_y.shape, np.nan)
    frac[usable] = rel_r[usable] / rel_y[usable]
    return pd.DataFrame(
        {
            "reliable_energy": rel_y,
            "reliable_residual_energy": rel_r,
            "D_unexplained_fraction": frac,
            "D_explained_fraction": 1.0 - frac,
            "D_raw_residual_fraction": np.einsum("...g,...g->...", resid, resid)
            / np.maximum(np.einsum("...g,...g->...", pooled, pooled), EPS),
            "D_raw_gamma_norm_proxy": np.linalg.norm(resid, axis=-1),
            "usable": usable,
        }
    )


# --------------------------------------------------------------------------
# B. confound-controlled association
# --------------------------------------------------------------------------


def partial_spearman(x: np.ndarray, y: np.ndarray, covariates: np.ndarray) -> float:
    """Spearman correlation of ``x`` and ``y`` after linearly removing covariates.

    Ranks are taken first, covariates are regressed out of both rank vectors by
    ordinary least squares with an intercept, and the residuals are correlated.
    Rows with any non-finite value are dropped.
    """
    from scipy import stats

    X = np.asarray(x, dtype=np.float64)
    Y = np.asarray(y, dtype=np.float64)
    C = np.atleast_2d(np.asarray(covariates, dtype=np.float64))
    if C.shape[0] != X.shape[0]:
        C = C.T
    ok = np.isfinite(X) & np.isfinite(Y) & np.all(np.isfinite(C), axis=1)
    if ok.sum() < 10:
        return np.nan
    xr = stats.rankdata(X[ok])
    yr = stats.rankdata(Y[ok])
    Cr = np.column_stack([stats.rankdata(C[ok, j]) for j in range(C.shape[1])])
    design = np.column_stack([np.ones(ok.sum()), Cr])
    bx, *_ = np.linalg.lstsq(design, xr, rcond=None)
    by, *_ = np.linalg.lstsq(design, yr, rcond=None)
    rx = xr - design @ bx
    ry = yr - design @ by
    if rx.std() < EPS or ry.std() < EPS:
        return np.nan
    return float(np.corrcoef(rx, ry)[0, 1])
