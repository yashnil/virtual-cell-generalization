r"""Low-capacity pathway residual model, v1.

Three conceptually separate pieces, kept separate on purpose:

1. **scale-calibrated conserved transfer** — the baseline ``B``,
2. **pathway-level context-specific correction** — the learned residual ``R_hat``,
3. **source-agreement confidence** — carried alongside, never fed to the model
   as a target.

Lives in its own namespace so that no frozen analysis module has to change.

Residual algebra (read this before interpreting any result)
-----------------------------------------------------------
At pathway level with the four-context decomposition
``Y[c,p] = mu + alpha_c + beta_p + gamma[c,p]``, a source set ``S`` of size 3 and
a fitted scale ``s``, the baseline is

``B = mean_p A + s (A[p] - mean_p A)`` with ``A[p] = mean_{c in S} Y[c,p]``

so, using ``mean_S alpha = -alpha_c/3`` and ``mean_S gamma = -gamma[c,p]/3``,

.. math::
    R = Y - B = \tfrac{4}{3}\alpha_c + (1-s)\beta_p + (1 + \tfrac{s}{3})\gamma[c,p]

**R is therefore NOT the ANOVA gamma.** It is a mixture of three things:

* ``(4/3) alpha_c`` — a context-wide template offset, removed by the
  prediction-time centring over perturbations,
* ``(1-s) beta_p`` — **the part of the conserved effect the shrinkage
  deliberately left unpredicted**. With ``s ~ 0.44`` this is over half of beta,
  and it is the reason a model could appear to "work" merely by undoing the
  shrinkage rather than by learning any interaction,
* ``(1 + s/3) gamma[c,p]`` — the interaction, mildly amplified.

Centring ``R`` over perturbations removes the alpha term exactly (and leaves the
other two, since ``sum_p beta_p = 0`` makes the beta term zero-mean but not
zero). Because of the beta contamination this module reports, alongside the
headline prediction metrics, the correlation of ``R_hat`` with ``gamma`` and with
``beta`` separately, so that "learned interaction" and "undid the shrinkage" can
be told apart.

Leakage contract
----------------
For outer target ``t`` nothing derived from ``Y[t]`` may enter features, scaling,
model selection, the scale estimate or the shrinkage. Basal control profiles of
all contexts, including ``t``, are permitted. Enforced by tests that overwrite
``Y[t]`` with noise and require bit-identical features and predictions.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

EPS = 1e-12

#: Predeclared shrinkage grid for the interaction correction.
LAMBDA_GRID: tuple[float, ...] = (0.0, 0.25, 0.5, 0.75, 1.0)

#: Predeclared ridge penalty grid.
ALPHA_GRID: tuple[float, ...] = (1.0, 10.0, 100.0, 1000.0, 10000.0)

FEATURE_NAMES: tuple[str, ...] = (
    "source_mean",
    "source_mean_centred",
    "source_sd",
    "source_range",
    "source_sign_agreement",
    "weighted_source",
    "weighted_source_centred",
    "weighted_minus_mean",
    "source_agreement",
    "source_reliability",
    "source_cells",
    "source_magnitude",
    "target_gene_basal",
    "basal_pathway_deviation",
    "target_basal_pathway",
    "gene_in_pathway",
    "pathway_size",
)


# --------------------------------------------------------------------------
# scale calibration and baseline
# --------------------------------------------------------------------------


def fit_scale(Y: np.ndarray, contexts: Sequence[int]) -> float:
    r"""Scalar shrinkage of the perturbation-specific part, from ``contexts`` only.

    Inner leave-one-context-out within the supplied set: each member is predicted
    from the mean of the others and the pooled least-squares scalar is returned.
    Never sees any context outside ``contexts``.

    With only two contexts available each inner fold has a single predictor
    context, whose gamma attenuation differs from the three-source case. The
    estimate is then cruder; that mismatch is real and is reported rather than
    hidden.
    """
    ctx = list(contexts)
    if len(ctx) < 2:
        return 1.0
    num = den = 0.0
    for held in ctx:
        others = [c for c in ctx if c != held]
        inner = np.asarray(Y, dtype=np.float64)[others].mean(axis=0)
        inner_c = inner - inner.mean(axis=0, keepdims=True)
        y_c = np.asarray(Y, dtype=np.float64)[held]
        y_c = y_c - y_c.mean(axis=0, keepdims=True)
        num += float(np.sum(y_c * inner_c))
        den += float(np.sum(inner_c * inner_c))
    return float(num / den) if den > EPS else 1.0


def baseline(Y: np.ndarray, sources: Sequence[int], scale: float) -> np.ndarray:
    """Scale-calibrated conserved transfer ``B`` for one target, shape ``(P, K)``."""
    A = np.asarray(Y, dtype=np.float64)[list(sources)].mean(axis=0)
    template = A.mean(axis=0, keepdims=True)
    return template + scale * (A - template)


def centre_predictions(R_hat: np.ndarray) -> np.ndarray:
    """Remove the per-pathway mean over perturbations, using predictions only."""
    R = np.asarray(R_hat, dtype=np.float64)
    return R - R.mean(axis=0, keepdims=True)


# --------------------------------------------------------------------------
# inference-available features
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class GeneLevelStats:
    """Per-(context, perturbation) quantities computed once from gene space.

    All are properties of a context's own cells; using them for a *source*
    context is always legal, and none are read for the outer target.
    """

    pair_corr: np.ndarray = field(repr=False)  # (C, C, P) gene-level agreement
    reliability: np.ndarray = field(repr=False)  # (C, P)
    cells: np.ndarray = field(repr=False)  # (C, P)
    magnitude: np.ndarray = field(repr=False)  # (C, P)
    gene_basal: np.ndarray = field(repr=False)  # (C, P) target-gene basal expression


def build_features(
    Y: np.ndarray,
    target: int,
    sources: Sequence[int],
    *,
    stats: GeneLevelStats,
    weights: np.ndarray,
    basal_pathway: np.ndarray,
    in_pathway: np.ndarray,
    pathway_size: np.ndarray,
) -> np.ndarray:
    """Design matrix of shape ``(P * K, F)`` for predicting ``target`` from ``sources``.

    ``Y[target]`` is never read. ``basal_pathway`` (control profiles projected to
    pathways) and ``stats.gene_basal`` *are* read for the target, because basal
    control expression is available at inference.
    """
    Ys = np.asarray(Y, dtype=np.float64)[list(sources)]  # (S, P, K)
    n_p, n_k = Ys.shape[1], Ys.shape[2]
    src = list(sources)

    A = Ys.mean(axis=0)
    Ac = A - A.mean(axis=0, keepdims=True)
    sd = Ys.std(axis=0)
    rng_ = Ys.max(axis=0) - Ys.min(axis=0)
    sign_agree = np.abs(A) / (np.abs(Ys).mean(axis=0) + EPS)
    wA = np.tensordot(np.asarray(weights, dtype=np.float64), Ys, axes=(0, 0))
    wAc = wA - wA.mean(axis=0, keepdims=True)

    pairs = [(a, b) for i, a in enumerate(src) for b in src[i + 1 :]]
    agree = np.nanmean(np.vstack([stats.pair_corr[a, b] for a, b in pairs]), axis=0)
    rel = np.nanmean(stats.reliability[src], axis=0)
    cells = stats.cells[src].mean(axis=0)
    mag = stats.magnitude[src].mean(axis=0)
    gene_basal = stats.gene_basal[target]

    basal_dev = basal_pathway[target] - basal_pathway[src].mean(axis=0)  # (K,)
    tgt_basal = basal_pathway[target]  # (K,)

    ones_p = np.ones((n_p, 1))
    ones_k = np.ones((1, n_k))
    columns = [
        A,
        Ac,
        sd,
        rng_,
        sign_agree,
        wA,
        wAc,
        wAc - Ac,
        agree[:, None] * ones_k,
        rel[:, None] * ones_k,
        cells[:, None] * ones_k,
        mag[:, None] * ones_k,
        gene_basal[:, None] * ones_k,
        ones_p * basal_dev[None, :],
        ones_p * tgt_basal[None, :],
        np.asarray(in_pathway, dtype=np.float64),
        ones_p * np.asarray(pathway_size, dtype=np.float64)[None, :],
    ]
    X = np.stack([c.reshape(-1) for c in columns], axis=1)
    return np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)


class Standardiser:
    """Mean/SD standardisation fitted on training rows only."""

    def __init__(self) -> None:
        self.mean_: np.ndarray | None = None
        self.scale_: np.ndarray | None = None

    def fit(self, X: np.ndarray) -> Standardiser:
        Xa = np.asarray(X, dtype=np.float64)
        self.mean_ = Xa.mean(axis=0)
        sd = Xa.std(axis=0)
        self.scale_ = np.where(sd > EPS, sd, 1.0)
        return self

    def transform(self, X: np.ndarray) -> np.ndarray:
        if self.mean_ is None or self.scale_ is None:
            raise RuntimeError("Standardiser used before fit().")
        return (np.asarray(X, dtype=np.float64) - self.mean_) / self.scale_


# --------------------------------------------------------------------------
# models
# --------------------------------------------------------------------------


def fit_ridge(X: np.ndarray, y: np.ndarray, alpha: float) -> np.ndarray:
    """Closed-form ridge with an unpenalised intercept. Returns ``(F + 1,)``."""
    Xa = np.asarray(X, dtype=np.float64)
    ya = np.asarray(y, dtype=np.float64).ravel()
    n, f = Xa.shape
    Xb = np.hstack([np.ones((n, 1)), Xa])
    penalty = np.eye(f + 1) * alpha
    penalty[0, 0] = 0.0
    return np.linalg.solve(Xb.T @ Xb + penalty, Xb.T @ ya)


def predict_ridge(X: np.ndarray, coef: np.ndarray) -> np.ndarray:
    Xa = np.asarray(X, dtype=np.float64)
    return coef[0] + Xa @ coef[1:]


def fit_deterministic(feature_column: np.ndarray, y: np.ndarray) -> float:
    """One global scalar: least-squares fit of ``y`` on a single feature column."""
    x = np.asarray(feature_column, dtype=np.float64).ravel()
    ya = np.asarray(y, dtype=np.float64).ravel()
    den = float(x @ x)
    return float(x @ ya / den) if den > EPS else 0.0


# --------------------------------------------------------------------------
# metrics
# --------------------------------------------------------------------------


def rowwise_pearson(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    x = np.asarray(a, dtype=np.float64)
    y = np.asarray(b, dtype=np.float64)
    xc = x - x.mean(axis=-1, keepdims=True)
    yc = y - y.mean(axis=-1, keepdims=True)
    num = np.einsum("...g,...g->...", xc, yc)
    den = np.sqrt(np.einsum("...g,...g->...", xc, xc) * np.einsum("...g,...g->...", yc, yc))
    out = np.full(np.shape(num), np.nan)
    ok = den > EPS
    out[ok] = num[ok] / den[ok]
    return out


def prediction_metrics(observed: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    """Per-perturbation metrics across pathways, summarised by the median."""
    Y = np.asarray(observed, dtype=np.float64)
    P = np.asarray(predicted, dtype=np.float64)
    if Y.shape != P.shape:
        raise ValueError(f"shape mismatch {Y.shape} vs {P.shape}")
    err = Y - P
    sse = np.einsum("pk,pk->p", err, err)
    sst = np.einsum("pk,pk->p", Y, Y)
    cos_den = np.linalg.norm(Y, axis=1) * np.linalg.norm(P, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        energy = np.where(sst > EPS, 1.0 - sse / sst, np.nan)
        cosine = np.where(cos_den > EPS, np.einsum("pk,pk->p", Y, P) / cos_den, np.nan)
    return {
        "pearson": float(np.nanmedian(rowwise_pearson(Y, P))),
        "cosine": float(np.nanmedian(cosine)),
        "mse": float(np.nanmedian(sse / Y.shape[1])),
        "energy_explained": float(np.nanmedian(energy)),
    }


def per_perturbation_pearson(observed: np.ndarray, predicted: np.ndarray) -> np.ndarray:
    return rowwise_pearson(np.asarray(observed), np.asarray(predicted))


# --------------------------------------------------------------------------
# nested leave-one-context-out
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class FoldSpec:
    """Everything needed to build features and a baseline for one pseudo-target."""

    target: int
    sources: tuple[int, ...]


@dataclass(frozen=True)
class Selection:
    """Outcome of inner pseudo-LOCO model selection."""

    family: str
    alpha: float
    lam: float
    inner_scores: dict[str, float] = field(default_factory=dict)


def assemble_rows(
    Y: np.ndarray,
    spec: FoldSpec,
    *,
    stats: GeneLevelStats,
    weight_fn,
    basal_pathway: np.ndarray,
    in_pathway: np.ndarray,
    pathway_size: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Features, residual target, baseline and fitted scale for one pseudo-target.

    ``weight_fn(target, sources)`` supplies basal-derived source weights; it must
    not consult any perturbation response.
    """
    scale = fit_scale(Y, spec.sources)
    B = baseline(Y, spec.sources, scale)
    X = build_features(
        Y,
        spec.target,
        spec.sources,
        stats=stats,
        weights=weight_fn(spec.target, spec.sources),
        basal_pathway=basal_pathway,
        in_pathway=in_pathway,
        pathway_size=pathway_size,
    )
    R = (np.asarray(Y, dtype=np.float64)[spec.target] - B).reshape(-1)
    return X, R, B, scale


def _predict_residual(family: str, X: np.ndarray, model, shape: tuple[int, int]) -> np.ndarray:
    if family == "M0":
        return np.zeros(shape)
    if family == "M1":
        col = FEATURE_NAMES.index("weighted_minus_mean")
        return (model * X[:, col]).reshape(shape)
    return predict_ridge(X, model).reshape(shape)


def _fit_family(family: str, X: np.ndarray, R: np.ndarray, alpha: float):
    if family == "M0":
        return None
    if family == "M1":
        return fit_deterministic(X[:, FEATURE_NAMES.index("weighted_minus_mean")], R)
    return fit_ridge(X, R, alpha)


def select_by_inner_loco(
    Y: np.ndarray,
    sources: Sequence[int],
    *,
    stats: GeneLevelStats,
    weight_fn,
    basal_pathway: np.ndarray,
    in_pathway: np.ndarray,
    pathway_size: np.ndarray,
    families: Sequence[str] = ("M0", "M1", "M2"),
    alphas: Sequence[float] = ALPHA_GRID,
    lambdas: Sequence[float] = LAMBDA_GRID,
) -> Selection:
    """Choose family, ridge penalty and shrinkage using source contexts only.

    Each source context is held out in turn as a pseudo-target; the model is
    trained on rows from the remaining source contexts and scored on the held-out
    one. **The outer target is never referenced.** The predeclared inner
    criterion is the median per-perturbation Pearson of the corrected prediction.
    """
    src = list(sources)
    prepared = {}
    for c in src:
        spec = FoldSpec(c, tuple(s for s in src if s != c))
        prepared[c] = assemble_rows(
            Y,
            spec,
            stats=stats,
            weight_fn=weight_fn,
            basal_pathway=basal_pathway,
            in_pathway=in_pathway,
            pathway_size=pathway_size,
        )

    n_p, n_k = np.asarray(Y).shape[1], np.asarray(Y).shape[2]
    scores: dict[str, float] = {}
    best = None
    for family in families:
        for alpha in alphas if family == "M2" else (0.0,):
            per_lambda = {lam: [] for lam in lambdas}
            for held in src:
                train = [c for c in src if c != held]
                Xtr = np.vstack([prepared[c][0] for c in train])
                Rtr = np.concatenate([prepared[c][1] for c in train])
                scaler = Standardiser().fit(Xtr)
                model = _fit_family(family, scaler.transform(Xtr), Rtr, alpha)
                Xva, _, Bva, _ = prepared[held]
                R_hat = centre_predictions(
                    _predict_residual(family, scaler.transform(Xva), model, (n_p, n_k))
                )
                for lam in lambdas:
                    pred = Bva + lam * R_hat
                    per_lambda[lam].append(prediction_metrics(np.asarray(Y)[held], pred)["pearson"])
            for lam, vals in per_lambda.items():
                score = float(np.mean(vals))
                scores[f"{family}|alpha={alpha:g}|lam={lam:g}"] = score
                if best is None or score > best[0]:
                    best = (score, family, alpha, lam)
    assert best is not None
    return Selection(family=best[1], alpha=best[2], lam=best[3], inner_scores=scores)


def fit_and_predict_outer(
    Y: np.ndarray,
    target: int,
    sources: Sequence[int],
    selection: Selection,
    *,
    stats: GeneLevelStats,
    weight_fn,
    basal_pathway: np.ndarray,
    in_pathway: np.ndarray,
    pathway_size: np.ndarray,
) -> dict[str, object]:
    """Refit on all source pseudo-targets and predict the outer target once."""
    src = list(sources)
    Xtr, Rtr = [], []
    for c in src:
        spec = FoldSpec(c, tuple(s for s in src if s != c))
        X, R, _, _ = assemble_rows(
            Y,
            spec,
            stats=stats,
            weight_fn=weight_fn,
            basal_pathway=basal_pathway,
            in_pathway=in_pathway,
            pathway_size=pathway_size,
        )
        Xtr.append(X)
        Rtr.append(R)
    Xtr = np.vstack(Xtr)
    Rtr = np.concatenate(Rtr)
    scaler = Standardiser().fit(Xtr)
    model = _fit_family(selection.family, scaler.transform(Xtr), Rtr, selection.alpha)

    scale = fit_scale(Y, src)
    B = baseline(Y, src, scale)
    X_out = build_features(
        Y,
        target,
        src,
        stats=stats,
        weights=weight_fn(target, src),
        basal_pathway=basal_pathway,
        in_pathway=in_pathway,
        pathway_size=pathway_size,
    )
    n_p, n_k = np.asarray(Y).shape[1], np.asarray(Y).shape[2]
    R_hat = centre_predictions(
        _predict_residual(selection.family, scaler.transform(X_out), model, (n_p, n_k))
    )
    return {
        "baseline": B,
        "residual_hat": R_hat,
        "prediction": B + selection.lam * R_hat,
        "scale": scale,
        "model": model,
        "scaler": scaler,
        "n_train_rows": int(Xtr.shape[0]),
    }
