"""N3 leakage controls (reports/n3_protocol.md §5)."""

from __future__ import annotations

import numpy as np
import pytest

from virtual_cell.analysis import calibration_budget as cb
from virtual_cell.modelling import pathway_residual as pr

P, G, TARGET = 200, 80, 4


def _world(rng):
    beta = rng.normal(size=(P, G))
    D = np.stack([beta + 0.3 * rng.normal(size=(P, G)) for _ in range(6)])
    fit = D[TARGET] + 0.1 * rng.normal(size=(P, G))
    return D, fit


@pytest.mark.parametrize("sources", [[0, 1, 2, 3, 5], [1, 5]])
def test_target_rows_and_non_anchors_never_read(sources):
    rng = np.random.default_rng(0)
    D, fit = _world(rng)
    T, K = np.arange(150, 200), np.arange(0, 30)

    def run(D_, fit_):
        sv = cb.make_source_view(D_, sources, pr.fit_scale(D_, sources))
        return cb.predict_all(sv, T, K, fit_[K])[0]

    ref = run(D, fit)
    D2, fit2 = D.copy(), fit.copy()
    D2[TARGET] = rng.normal(size=(P, G))
    mask = np.ones(P, bool)
    mask[K] = False
    fit2[mask] = rng.normal(size=(mask.sum(), G))
    new = run(D2, fit2)
    for name in ref:
        np.testing.assert_array_equal(ref[name], new[name])


def test_single_source_scale_is_one_and_estimators_run():
    rng = np.random.default_rng(1)
    D, fit = _world(rng)
    assert pr.fit_scale(D, [2]) == 1.0
    sv = cb.make_source_view(D, [2], 1.0)
    preds, _ = cb.predict_all(sv, np.arange(150, 200), np.arange(20), fit[:20])
    assert all(np.isfinite(v).all() for v in preds.values())
