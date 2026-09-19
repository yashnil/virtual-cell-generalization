"""Four-component ANOVA decomposition of pseudobulk perturbation responses.

Our own implementation of the decomposition used by Molina & Zhang (2026),
written to match their released reference implementation exactly. See
``reports/molina_zhang_reproduction_spec.md`` for the frozen specification and
for the provenance of every convention encoded here.

The model, on the balanced tensor of pseudobulk response vectors
``D[c, p] in R^G`` over cell lines ``c`` and perturbations ``p``::

    delta_{c,p,g} = mu_g + alpha_{c,g} + beta_{p,g} + gamma_{c,p,g}

estimated by plain cell means::

    mu    = mean_{c,p} D[c, p]
    alpha = mean_p D[c, p] - mu
    beta  = mean_c D[c, p] - mu
    gamma = D - mu - alpha - beta

This is the classical two-way balanced ANOVA, so the estimates satisfy exact
reconstruction and the zero-sum (sum-to-zero) side conditions, and the four
components are mutually orthogonal under the sum-of-squares convention in
:func:`sums_of_squares`.

Two conventions matter and are easy to get wrong; both follow the reference
implementation rather than any textbook default:

* **Sum-of-squares denominator.** Totals are *uncentered* per-observation
  squared norms averaged over observations, ``mean_{c,p} ||D[c,p]||^2`` — not
  ``np.var``. The ``mu`` component therefore carries a real share of the total.
* **Template.** The reported "template" (cell-line) share is ``mu + alpha``
  combined, not ``alpha`` alone.

Noise is not part of the algebraic decomposition. It is estimated separately by
split-half resampling over cells (:func:`split_half_signal`), which measures the
*reproducible* part of each component; whatever the four reproducible components
do not explain is reported as noise.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np

#: Cell-line order used by the reference implementation's figure scripts.
REFERENCE_CELL_LINES: tuple[str, ...] = ("k562", "rpe1", "hepg2", "jurkat")

#: Minimum cells per perturbation used by the reference decomposition figure.
REFERENCE_MIN_CELLS = 10

#: Split-half resamples and seed used by the reference implementation.
REFERENCE_N_SPLITS = 50
REFERENCE_SEED = 42


class DecompositionError(ValueError):
    """Raised when an input violates a decomposition assumption."""


# --------------------------------------------------------------------------
# tensor construction
# --------------------------------------------------------------------------


def shared_perturbations(cl_deltas: Mapping[str, Mapping[str, np.ndarray]]) -> list[str]:
    """Perturbations measured in *every* cell line, sorted by name.

    Sorting is part of the specification: it makes the perturbation axis a
    deterministic function of the input, independent of dict insertion order.
    """
    if not cl_deltas:
        raise DecompositionError("No cell lines provided.")
    shared = set.intersection(*[set(cl_deltas[cl]) for cl in cl_deltas])
    if not shared:
        raise DecompositionError("Cell lines share no perturbations.")
    return sorted(shared)


def build_response_tensor(
    cl_deltas: Mapping[str, Mapping[str, np.ndarray]],
    *,
    cell_lines: Sequence[str] | None = None,
    perturbations: Sequence[str] | None = None,
    dtype: np.dtype | str = np.float64,
) -> tuple[np.ndarray, list[str], list[str]]:
    """Stack per-cell-line delta dicts into a balanced ``(C, P, G)`` tensor.

    Parameters
    ----------
    cl_deltas:
        ``cl_deltas[cell_line][perturbation] -> (G,) delta vector``.
    cell_lines:
        Explicit cell-line order. Defaults to sorted order. Pass
        :data:`REFERENCE_CELL_LINES` to match the reference figure scripts.
    perturbations:
        Explicit perturbation order. Defaults to the shared set, sorted. Passing
        an explicit list is how a *frozen* perturbation set is enforced.
    dtype:
        Accumulation dtype. Defaults to float64; the reference implementation
        builds the tensor in float64 in its figure path and float32 in its
        module path, and float64 is the safer default for the sums of squares.

    Returns
    -------
    (tensor, cell_lines, perturbations)
    """
    cls = list(cell_lines) if cell_lines is not None else sorted(cl_deltas)
    missing_cls = [c for c in cls if c not in cl_deltas]
    if missing_cls:
        raise DecompositionError(f"Cell lines absent from the input: {missing_cls}")

    perts = (
        list(perturbations)
        if perturbations is not None
        else shared_perturbations({c: cl_deltas[c] for c in cls})
    )
    if not perts:
        raise DecompositionError("Empty perturbation set.")

    for cl in cls:
        absent = [p for p in perts if p not in cl_deltas[cl]]
        if absent:
            raise DecompositionError(
                f"Unbalanced design: cell line {cl!r} is missing "
                f"{len(absent)} perturbation(s), e.g. {absent[:5]}."
            )

    widths = {int(np.asarray(cl_deltas[c][p]).shape[-1]) for c in cls for p in perts}
    if len(widths) != 1:
        raise DecompositionError(f"Inconsistent response dimensionality: {sorted(widths)}")

    tensor = np.array(
        [[np.asarray(cl_deltas[cl][p], dtype=dtype).ravel() for p in perts] for cl in cls],
        dtype=dtype,
    )
    return tensor, cls, perts


def _validate_tensor(D: np.ndarray) -> None:
    if D.ndim != 3:
        raise DecompositionError(f"Expected a 3-D (C, P, G) tensor, got shape {D.shape}.")
    if min(D.shape) == 0:
        raise DecompositionError(f"Tensor has an empty axis: shape {D.shape}.")
    if not np.all(np.isfinite(D)):
        raise DecompositionError("Tensor contains NaN or infinite values.")


# --------------------------------------------------------------------------
# the decomposition
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Decomposition:
    """Estimated components of a balanced two-way response decomposition."""

    mu: np.ndarray = field(repr=False)  # (G,)
    alpha: np.ndarray = field(repr=False)  # (C, G)
    beta: np.ndarray = field(repr=False)  # (P, G)
    gamma: np.ndarray = field(repr=False)  # (C, P, G)
    cell_lines: tuple[str, ...] = ()
    perturbations: tuple[str, ...] = ()

    @property
    def shape(self) -> tuple[int, int, int]:
        return self.gamma.shape  # type: ignore[return-value]

    def reconstruct(self) -> np.ndarray:
        """Rebuild the response tensor from the components (exact by construction)."""
        return self.mu[None, None, :] + self.alpha[:, None, :] + self.beta[None, :, :] + self.gamma


def decompose(
    D: np.ndarray,
    *,
    cell_lines: Sequence[str] = (),
    perturbations: Sequence[str] = (),
) -> Decomposition:
    """Estimate ``mu``, ``alpha``, ``beta`` and ``gamma`` from a ``(C, P, G)`` tensor."""
    _validate_tensor(D)
    mu = D.mean(axis=(0, 1))
    alpha = D.mean(axis=1) - mu
    beta = D.mean(axis=0) - mu
    gamma = D - mu[None, None, :] - alpha[:, None, :] - beta[None, :, :]
    return Decomposition(
        mu=mu,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        cell_lines=tuple(cell_lines),
        perturbations=tuple(perturbations),
    )


# --------------------------------------------------------------------------
# sums of squares
# --------------------------------------------------------------------------


def total_sum_of_squares(D: np.ndarray) -> float:
    """``mean_{c,p} ||D[c,p]||^2`` — the reference total.

    Note this is *uncentered*: it is the mean squared norm of each observed
    response vector, not a variance about the grand mean.
    """
    _validate_tensor(D)
    return float(np.mean(np.sum(D**2, axis=2)))


def sums_of_squares(dec: Decomposition) -> dict[str, float]:
    """Per-component sums of squares on the reference convention.

    Each component is averaged over the axes it does *not* span, so the four
    values sum exactly to :func:`total_sum_of_squares` of the reconstructed
    tensor (the components are mutually orthogonal in a balanced design).
    """
    return {
        "mu": float(np.sum(dec.mu**2)),
        "alpha": float(np.mean(np.sum(dec.alpha**2, axis=1))),
        "beta": float(np.mean(np.sum(dec.beta**2, axis=1))),
        "gamma": float(np.mean(np.sum(dec.gamma**2, axis=2))),
    }


def uncorrected_fractions(dec: Decomposition, total: float | None = None) -> dict[str, float]:
    """Component shares of the total sum of squares, with no noise correction.

    These are *not* the numbers reported in the paper: without the split-half
    correction every component is inflated by its own measurement noise, and
    ``gamma`` absorbs almost all of it.
    """
    ss = sums_of_squares(dec)
    denom = total if total is not None else sum(ss.values())
    if denom <= 0:
        raise DecompositionError("Total sum of squares is not positive.")
    return {k: v / denom for k, v in ss.items()}


# --------------------------------------------------------------------------
# projective template removal
# --------------------------------------------------------------------------


def project_out_template(D: np.ndarray) -> np.ndarray:
    """Remove each cell line's own template *direction* from its responses.

    ``eps[c,p] = D[c,p] - (D[c,p] . That_c) That_c`` where ``That_c`` is the unit
    vector along ``mean_p D[c,p]``. This is a projection, not a subtraction, so
    it removes the template axis rather than a fixed offset.
    """
    _validate_tensor(D)
    out = np.zeros_like(D)
    for c in range(D.shape[0]):
        template = D[c].mean(axis=0)
        norm = float(np.linalg.norm(template))
        unit = template / (norm + 1e-30)
        out[c] = D[c] - np.outer(D[c] @ unit, unit)
    return out


# --------------------------------------------------------------------------
# split-half noise correction
# --------------------------------------------------------------------------


def cross_half_signal(D1: np.ndarray, D2: np.ndarray) -> dict[str, float]:
    """Reproducible sum of squares per component, from two independent halves.

    Each half is decomposed on its own, then matching components are compared by
    dot product. Independent measurement noise has zero expected cross-product,
    so these estimate the *signal* sum of squares of each component on the same
    scale as :func:`sums_of_squares`.
    """
    _validate_tensor(D1)
    _validate_tensor(D2)
    if D1.shape != D2.shape:
        raise DecompositionError(f"Half shapes differ: {D1.shape} vs {D2.shape}.")

    a = decompose(D1)
    b = decompose(D2)
    return {
        "mu": float(np.dot(a.mu, b.mu)),
        "alpha": float(np.mean(np.sum(a.alpha * b.alpha, axis=1))),
        "beta": float(np.mean(np.sum(a.beta * b.beta, axis=1))),
        "gamma": float(np.mean(np.sum(a.gamma * b.gamma, axis=2))),
    }


def split_half_delta_tensors(
    cell_data: Mapping[str, Mapping[str, object]],
    *,
    cell_lines: Sequence[str],
    perturbations: Sequence[str],
    rng: np.random.RandomState,
) -> tuple[np.ndarray, np.ndarray]:
    """Build two delta tensors from disjoint halves of each perturbation's cells.

    ``cell_data[cl]`` must provide ``"ctrl_mean"`` (G,) and ``"pert_cells"``
    mapping perturbation -> ``(n_cells, G)``.

    Following the reference implementation, both halves are referenced against
    the **full-data** control mean (the controls are not split), and an odd cell
    is dropped so the two halves are exactly equal in size.
    """
    n_c, n_p = len(cell_lines), len(perturbations)
    G = int(np.asarray(cell_data[cell_lines[0]]["ctrl_mean"]).shape[-1])
    D1 = np.zeros((n_c, n_p, G))
    D2 = np.zeros((n_c, n_p, G))
    for ci, cl in enumerate(cell_lines):
        ctrl_mean = np.asarray(cell_data[cl]["ctrl_mean"], dtype=np.float64).ravel()
        pert_cells = cell_data[cl]["pert_cells"]  # type: ignore[index]
        for pi, p in enumerate(perturbations):
            cells = np.asarray(pert_cells[p], dtype=np.float64)  # type: ignore[index]
            perm = rng.permutation(len(cells))
            half = len(cells) // 2
            if half == 0:
                raise DecompositionError(
                    f"Cell line {cl!r} perturbation {p!r} has {len(cells)} cell(s); "
                    "at least 2 are needed to split."
                )
            D1[ci, pi] = cells[perm[:half]].mean(axis=0) - ctrl_mean
            D2[ci, pi] = cells[perm[half : 2 * half]].mean(axis=0) - ctrl_mean
    return D1, D2


def split_half_signal(
    cell_data: Mapping[str, Mapping[str, object]],
    *,
    cell_lines: Sequence[str],
    perturbations: Sequence[str],
    n_splits: int = REFERENCE_N_SPLITS,
    seed: int = REFERENCE_SEED,
    project_template: bool = False,
) -> dict[str, float]:
    """Average reproducible sums of squares over ``n_splits`` split-half resamples.

    A single ``RandomState(seed)`` is threaded through every resample, matching
    the reference implementation, so the whole sequence is reproducible.

    With ``project_template=True`` each half is template-projected independently
    before decomposition, which is how the template-removed residual bar is
    computed.
    """
    if n_splits < 1:
        raise DecompositionError("n_splits must be at least 1.")
    rng = np.random.RandomState(seed)
    accumulated: list[dict[str, float]] = []
    for _ in range(n_splits):
        D1, D2 = split_half_delta_tensors(
            cell_data, cell_lines=cell_lines, perturbations=perturbations, rng=rng
        )
        if project_template:
            D1 = project_out_template(D1)
            D2 = project_out_template(D2)
        accumulated.append(cross_half_signal(D1, D2))
    return {k: float(np.mean([s[k] for s in accumulated])) for k in accumulated[0]}


def noise_corrected_fractions(signal: Mapping[str, float], total: float) -> dict[str, float]:
    """Paper-facing component shares: template, conserved, interaction, noise.

    ``template`` is ``mu + alpha`` combined, and ``noise`` is whatever share of
    the total the four reproducible components leave unexplained. The four
    returned fractions sum to 1 by construction.
    """
    if total <= 0:
        raise DecompositionError("Total sum of squares is not positive.")
    reproducible = sum(signal[k] for k in ("mu", "alpha", "beta", "gamma"))
    return {
        "template": (signal["mu"] + signal["alpha"]) / total,
        "beta": signal["beta"] / total,
        "gamma": signal["gamma"] / total,
        "noise": (total - reproducible) / total,
    }


def beta_fraction_per_perturbation(dec: Decomposition) -> dict[str, float]:
    """Per-perturbation transferable share ``var(beta_p) / (var(beta_p) + mean_c var(gamma_cp))``.

    A descriptive "how conserved is this perturbation" statistic, following the
    reference implementation's ``beta_fraction``. Uses variance across genes.
    """
    names = dec.perturbations or tuple(str(i) for i in range(dec.beta.shape[0]))
    out = {}
    for j, name in enumerate(names):
        b = float(np.var(dec.beta[j]))
        g = float(np.mean([np.var(dec.gamma[c, j]) for c in range(dec.gamma.shape[0])]))
        out[name] = b / (b + g + 1e-30)
    return out
