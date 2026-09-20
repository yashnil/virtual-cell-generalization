r"""Pathway residual model v2 — clean (beta-free) interaction target.

v1 is frozen and its conclusion stands: *the scale-calibrated pathway residual
correction did not improve outer zero-shot response prediction.* v2 changes
**exactly one thing** — the training target — because v1's diagnostics identified
a specific algebraic contamination as the mechanism of failure.

Everything else is reused verbatim from :mod:`virtual_cell.modelling.pathway_residual`:
the features, the standardiser, the ridge, the baseline and its scale fitting,
the nested-LOCO structure, the shrinkage grid and the leakage contract.

Why the target changes
----------------------
v1 trained on ``Y - s*A``-style residuals whose centred form is

``(1 - s) beta_p + (1 + s/3) gamma[c,p]``

which at ``s ~ 0.55`` leaves nearly half the conserved effect in the target. v1
showed the model duly learned that beta component **with the wrong sign** and
cancelled its own genuine gamma gain.

The clean target uses the **unshrunk** transfer ``A[p] = mean_{s in S} Y[s,p]``:

.. math::
    R_1 = Y[c] - A = (\alpha_c - \overline{\alpha}_{S})
                   + (\gamma[c,p] - \overline{\gamma}_{S,p})

``beta_p`` appears with coefficient 1 in both terms and **cancels exactly, for
any source-set size, before any centring**. Centring over perturbations then
removes the constant ``alpha`` term, leaving

.. math::
    \mathrm{centre}_p(R_1) = \gamma[c,p] - \overline{\gamma}_{S,p}

* with all three other contexts as sources (the outer case),
  ``= (4/3) gamma[c,p]`` — a clean multiple of the interaction;
* with two sources ``{a, b}`` (the inner pseudo-target case),
  ``= gamma[c,p] - (gamma[a,p] + gamma[b,p])/2``.

So the inner training target is beta-free but carries the *other* contexts'
interactions, and its scale differs from the outer target's. That mismatch is
inherent to four contexts and is absorbed, imperfectly, by the shrinkage.

Theory coefficient
------------------
With the v1 baseline ``B = template + s * centre_p(A)``, the algebra gives
``B = mu - alpha_c/3 + s(beta_p - gamma[c,p]/3)``, so the gamma deficit relative
to the truth is ``(1 + s/3) gamma``. A correction predicting ``(4/3) gamma``
should therefore enter with

.. math::
    \lambda_{theory} = \tfrac{3}{4}\left(1 + \tfrac{s}{3}\right) = \tfrac{3 + s}{4}

reported as a diagnostic only. It is never tuned against the outer target.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from virtual_cell.modelling.pathway_residual import (
    ALPHA_GRID,
    LAMBDA_GRID,
    FoldSpec,
    Selection,
    Standardiser,
    _fit_family,
    _predict_residual,
    baseline,
    build_features,
    centre_predictions,
    fit_scale,
    prediction_metrics,
    rowwise_pearson,
)

__all__ = [
    "ALPHA_GRID",
    "LAMBDA_GRID",
    "FoldSpec",
    "Selection",
    "assemble_rows",
    "bootstrap_delta",
    "clean_gamma_target",
    "fit_and_predict_outer",
    "lambda_theory",
    "select_by_inner_loco",
    "unshrunk_transfer",
]

EPS = 1e-12


def unshrunk_transfer(Y: np.ndarray, sources: Sequence[int]) -> np.ndarray:
    """``A[p] = mean over source contexts`` — no scale applied."""
    return np.asarray(Y, dtype=np.float64)[list(sources)].mean(axis=0)


def clean_gamma_target(Y: np.ndarray, target: int, sources: Sequence[int]) -> np.ndarray:
    """Beta-free interaction target ``centre_p(Y[target] - A)``.

    Centring uses the **observed** responses of ``target``, so this is legitimate
    for a training pseudo-target and is **evaluation-only** for the real outer
    target.
    """
    A = unshrunk_transfer(Y, sources)
    return centre_predictions(np.asarray(Y, dtype=np.float64)[target] - A)


def lambda_theory(scale: float) -> float:
    """``(3 + s) / 4`` — see the module docstring for the derivation."""
    return (3.0 + float(scale)) / 4.0


def assemble_rows(
    Y: np.ndarray,
    spec: FoldSpec,
    *,
    stats,
    weight_fn,
    basal_pathway: np.ndarray,
    in_pathway: np.ndarray,
    pathway_size: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Features (identical to v1), **clean** gamma target, baseline and scale."""
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
    target = clean_gamma_target(Y, spec.target, spec.sources).reshape(-1)
    return X, target, B, scale


def select_by_inner_loco(
    Y: np.ndarray,
    sources: Sequence[int],
    *,
    stats,
    weight_fn,
    basal_pathway: np.ndarray,
    in_pathway: np.ndarray,
    pathway_size: np.ndarray,
    families: Sequence[str] = ("M0", "M1", "M2"),
    alphas: Sequence[float] = ALPHA_GRID,
    lambdas: Sequence[float] = LAMBDA_GRID,
) -> Selection:
    """Inner pseudo-LOCO selection on the clean target. Outer target never read."""
    src = list(sources)
    prepared = {
        c: assemble_rows(
            Y,
            FoldSpec(c, tuple(s for s in src if s != c)),
            stats=stats,
            weight_fn=weight_fn,
            basal_pathway=basal_pathway,
            in_pathway=in_pathway,
            pathway_size=pathway_size,
        )
        for c in src
    }
    n_p, n_k = np.asarray(Y).shape[1], np.asarray(Y).shape[2]
    scores: dict[str, float] = {}
    best = None
    for family in families:
        for alpha in alphas if family == "M2" else (0.0,):
            per_lambda: dict[float, list[float]] = {lam: [] for lam in lambdas}
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
                    per_lambda[lam].append(
                        prediction_metrics(np.asarray(Y)[held], Bva + lam * R_hat)["pearson"]
                    )
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
    stats,
    weight_fn,
    basal_pathway: np.ndarray,
    in_pathway: np.ndarray,
    pathway_size: np.ndarray,
) -> dict[str, object]:
    """Refit on all source pseudo-targets, predict the outer target once."""
    src = list(sources)
    Xtr, Rtr = [], []
    for c in src:
        X, R, _, _ = assemble_rows(
            Y,
            FoldSpec(c, tuple(s for s in src if s != c)),
            stats=stats,
            weight_fn=weight_fn,
            basal_pathway=basal_pathway,
            in_pathway=in_pathway,
            pathway_size=pathway_size,
        )
        Xtr.append(X)
        Rtr.append(R)
    Xtr_all = np.vstack(Xtr)
    Rtr_all = np.concatenate(Rtr)
    scaler = Standardiser().fit(Xtr_all)
    model = _fit_family(selection.family, scaler.transform(Xtr_all), Rtr_all, selection.alpha)

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
        "gamma_hat": R_hat,
        "prediction": B + selection.lam * R_hat,
        "scale": scale,
        "lambda_theory": lambda_theory(scale),
        "prediction_theory": B + lambda_theory(scale) * R_hat,
        "model": model,
        "n_train_rows": int(Xtr_all.shape[0]),
    }


def bootstrap_delta(
    observed: np.ndarray,
    baseline_pred: np.ndarray,
    corrected_pred: np.ndarray,
    *,
    n_boot: int = 2000,
    seed: int = 0,
    alpha: float = 0.05,
) -> dict[str, float]:
    """Bootstrap CI over perturbations for the change in median Pearson."""
    Y = np.asarray(observed, dtype=np.float64)
    rb = rowwise_pearson(Y, np.asarray(baseline_pred, dtype=np.float64))
    rc = rowwise_pearson(Y, np.asarray(corrected_pred, dtype=np.float64))
    ok = np.isfinite(rb) & np.isfinite(rc)
    rb, rc = rb[ok], rc[ok]
    point = float(np.median(rc) - np.median(rb))
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, rb.size, size=(n_boot, rb.size))
    boot = np.median(rc[idx], axis=1) - np.median(rb[idx], axis=1)
    return {
        "delta_median_pearson": point,
        "lo": float(np.percentile(boot, 100 * alpha / 2)),
        "hi": float(np.percentile(boot, 100 * (1 - alpha / 2))),
        "frac_improved": float(np.mean(rc > rb)),
        "n": int(rb.size),
    }
