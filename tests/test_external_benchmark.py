"""Checks on the external-validation protocol.

The important one is :func:`target_reliability`. Stage A of the external
diagnostic turned on it: ``kaden25rpe1`` produced a clean-looking negative
result that meant nothing, because its own responses did not reproduce. These
tests pin the statistic's behaviour at both ends -- a context with real,
repeatable effects must score high, and a context that is pure noise must score
near zero -- so that a future benchmark cannot be accepted on a number this
function got wrong.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from virtual_cell.modelling import external_benchmark as eb
from virtual_cell.priors import features as pf

GENES = [f"G{i}" for i in range(40)]


def _write(tmp_path, effects: dict[str, np.ndarray], n_cells: int, noise: float, seed: int = 0):
    """A minimal scPertEval-shaped .h5ad with known per-perturbation effects."""
    rng = np.random.default_rng(seed)
    blocks, labels = [], []
    base = rng.normal(5.0, 0.5, size=len(GENES))
    for name, effect in [*effects.items(), ("control", np.zeros(len(GENES)))]:
        block = base[None, :] + effect[None, :] + rng.normal(0, noise, size=(n_cells, len(GENES)))
        blocks.append(block)
        labels += [name] * n_cells
    adata = ad.AnnData(
        # The loaders require CSR, as the real scPertEval files are.
        X=sparse.csr_matrix(np.vstack(blocks).astype(np.float32)),
        obs=pd.DataFrame({"perturbation": labels}, index=[f"c{i}" for i in range(len(labels))]),
        var=pd.DataFrame(index=pd.Index(GENES)),
    )
    path = tmp_path / "ctx_processed_complete.h5ad"
    adata.write_h5ad(path)
    return path


def test_reliability_is_high_when_effects_are_real_and_repeatable(tmp_path) -> None:
    rng = np.random.default_rng(1)
    effects = {f"P{i}": rng.normal(0, 2.0, size=len(GENES)) for i in range(6)}
    path = _write(tmp_path, effects, n_cells=400, noise=0.5)
    frame = eb.target_reliability(path, genes=GENES, perturbations=sorted(effects))
    assert frame["spearman_brown"].median() > 0.9
    assert frame["reliability_ceiling"].median() > 0.9


def test_reliability_collapses_when_the_context_is_noise(tmp_path) -> None:
    """The kaden25rpe1 failure mode: plenty of cells, no reproducible effect."""
    effects = {f"P{i}": np.zeros(len(GENES)) for i in range(6)}
    path = _write(tmp_path, effects, n_cells=400, noise=3.0, seed=2)
    frame = eb.target_reliability(path, genes=GENES, perturbations=sorted(effects))
    assert abs(frame["spearman_brown"].median()) < 0.3


def test_reliability_reports_cells_and_effect_size(tmp_path) -> None:
    rng = np.random.default_rng(3)
    effects = {f"P{i}": rng.normal(0, 2.0, size=len(GENES)) for i in range(4)}
    path = _write(tmp_path, effects, n_cells=120, noise=0.5)
    frame = eb.target_reliability(path, genes=GENES, perturbations=sorted(effects))
    assert frame["n_cells"].tolist() == [120] * 4
    assert (frame["delta_norm"] > 0).all()
    assert list(frame["perturbation"]) == sorted(effects)


def test_reliability_subsamples_deterministically(tmp_path) -> None:
    rng = np.random.default_rng(4)
    effects = {f"P{i}": rng.normal(0, 2.0, size=len(GENES)) for i in range(8)}
    path = _write(tmp_path, effects, n_cells=60, noise=0.5)
    a = eb.target_reliability(path, genes=GENES, perturbations=sorted(effects), max_perturbations=3)
    b = eb.target_reliability(path, genes=GENES, perturbations=sorted(effects), max_perturbations=3)
    assert len(a) == 3
    assert list(a["perturbation"]) == list(b["perturbation"])
    assert np.allclose(a["spearman_brown"], b["spearman_brown"])


def test_prepare_matrix_standardises_and_flags_missingness() -> None:
    block = pf.FeatureBlock(
        "t",
        pd.Index(["A", "B", "C"]),
        np.array([[1.0, 10.0], [3.0, 30.0], [0.0, 0.0]]),
        np.array([True, True, False]),
    )
    out = eb.prepare_matrix(block)
    assert out.shape == (3, 3)
    # The covered rows are standardised against their own mean and spread.
    assert out[:2, 0].mean() == pytest.approx(0.0)
    # An uncovered row is filled with the covered mean, i.e. zero after centring,
    # and the indicator column records that it was filled.
    assert out[2, :2].tolist() == [0.0, 0.0]
    assert out[:, -1].tolist() == [1.0, 1.0, 0.0]


def test_the_frozen_constants_are_the_ones_arch1_was_scored_with() -> None:
    """These are not tunable. A change here invalidates every external result."""
    assert eb.PRIOR_FAMILY == "string"
    assert eb.KNN_K == 25
    assert eb.RIDGE_ALPHA == 10.0
    assert eb.PCA_COMPONENTS == 64
