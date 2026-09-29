"""C2 realisation generators: identity, unbiasedness, dispersion and the ideal DE ceiling."""

from __future__ import annotations

import numpy as np
import pytest

from virtual_cell.arc import metrics as M
from virtual_cell.competition_v2 import realisation


def _controls(rng, n=600, g=40):
    mean = rng.gamma(0.8, 3.0, size=g)
    depth = rng.integers(800, 3000, size=n)
    comp = mean / mean.sum()
    lam = depth[:, None] * comp[None, :] * rng.gamma(2.0, 0.5, size=(n, g))
    return rng.poisson(lam).astype(np.int64)


def test_effect_ratio_is_one_where_control_is_zero():
    r = realisation.effect_ratio(np.array([0.2, 0.3, 0.5]), np.array([0.4, 0.0, 0.6]))
    assert r[1] == 1.0
    np.testing.assert_allclose(r[[0, 2]], [0.5, 0.5 / 0.6])


def test_g1_identity_at_unit_ratio():
    rng = np.random.default_rng(0)
    x = _controls(rng)
    c = (x / x.sum(1, keepdims=True)).mean(0)
    out = realisation.g1_counts(x, np.ones(x.shape[1]), c, rng=rng)
    np.testing.assert_array_equal(out, x)


def test_g1_is_unbiased_for_the_mean_cpm_moment():
    rng = np.random.default_rng(1)
    x = _controls(rng, n=4000)
    c = (x / x.sum(1, keepdims=True)).mean(0)
    ratio = np.where(np.arange(x.shape[1]) % 2 == 0, 0.5, 1.8)
    target = c * ratio
    target /= target.sum()
    ratio = realisation.effect_ratio(target, c)
    out = realisation.g1_counts(x, ratio, c, rng=rng)
    got = (out / out.sum(1, keepdims=True)).mean(0)
    np.testing.assert_allclose(got, target, rtol=0.06, atol=2e-4)
    # up-regulated genes can switch on in cells that had none
    up = np.flatnonzero(ratio > 1)
    g = up[np.argmax((x[:, up] == 0).sum(0))]
    assert (x[:, g] == 0).sum() > 0
    assert (out[:, g] > 0).sum() > (x[:, g] > 0).sum()


def test_g1_rejects_bad_inputs():
    rng = np.random.default_rng(2)
    x = _controls(rng, n=10, g=5)
    with pytest.raises(ValueError):
        realisation.g1_counts(x, -np.ones(5), np.ones(5) / 5, rng=rng)
    with pytest.raises(ValueError):
        realisation.g1_counts(x, np.ones(4), np.ones(5) / 5, rng=rng)


def test_g2_preserves_depth_exactly():
    rng = np.random.default_rng(3)
    d = np.array([100, 2000, 50_000])
    out = realisation.g2_counts(d, np.array([1.0, 2.0, 7.0]), rng=rng)
    np.testing.assert_array_equal(out.sum(1), d)


def test_fit_dispersion_recovers_nb_and_poisson():
    rng = np.random.default_rng(4)
    n = 20_000
    depth = rng.integers(1000, 3000, size=n)
    comp = np.array([0.01, 0.02, 0.97])
    phi = np.array([0.5, 0.0, 0.0])
    x = realisation.g3_counts(depth, comp, phi, rng=rng)
    est = realisation.fit_dispersion(x, comp)
    assert abs(est[0] - 0.5) < 0.05
    assert est[1] < 0.01


def test_ideal_de_null_calls_nothing_and_shift_is_called():
    rng = np.random.default_rng(5)
    ref = _controls(rng, n=3000, g=60)
    cpm = ref / ref.sum(1, keepdims=True)
    tested = (1e6 * cpm).mean(0) > M.MIN_CPM
    null = cpm.mean(0)
    shifted = null.copy()
    shifted[0] *= 4
    shifted /= shifted.sum()
    table = realisation.ideal_de_table(np.stack([null, shifted]), ref, tested)
    assert not table.significant[0].any()
    assert table.significant[1][np.flatnonzero(tested).tolist().index(0)]
    assert table.lfc[1][np.flatnonzero(tested).tolist().index(0)] > 1.5


def test_cell_diagnostics_shapes():
    rng = np.random.default_rng(6)
    x = _controls(rng, n=150, g=30)
    d = realisation.cell_diagnostics(x, rng=rng)
    assert d["gene_var"].shape == (30,)
    assert 0 <= d["zero_fraction"] <= 1
    assert 0 <= d["cosine_distance_median"] <= 1
