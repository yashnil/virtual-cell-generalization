"""N4 checks (reports/n1_n4_protocol.md §3)."""

from __future__ import annotations

import numpy as np

from virtual_cell.analysis import agreement_null as an


def _world(P=200, G=300):
    rng = np.random.default_rng(0)
    m = rng.normal(size=(P, G)) * rng.uniform(0.05, 0.5, (P, 1))
    return m, np.full(P, 0.01), np.full((3, P), 0.02), np.full(P, 0.02)


def test_exchangeable_data_falls_inside_its_own_null():
    m, tau2, s2, st = _world()
    kw = dict(scale=0.5, min_signal=-1e9)
    obs = an.simulate_null(m, tau2, s2, st, rng=np.random.default_rng(1), **kw)
    null = [
        an.simulate_null(m, tau2, s2, st, rng=np.random.default_rng(i + 2), **kw) for i in range(40)
    ]
    lo, hi = np.percentile(null, [1, 99])
    assert lo <= obs <= hi


def test_quality_uses_per_repeat_energies():
    rng = np.random.default_rng(3)
    L = rng.normal(size=(50, 400))
    h1 = [L + rng.normal(size=L.shape) for _ in range(5)]
    h2 = [L + rng.normal(size=L.shape) for _ in range(5)]
    _, sig = an.quality_d(h1, h2, np.zeros_like(L), min_signal=-1e9)
    # unbiased for ||L||^2 despite unit noise per half
    np.testing.assert_allclose(sig.mean(), np.sum(L**2, axis=1).mean(), rtol=0.05)


def test_noise_variance_recovers_planted_sigma():
    rng = np.random.default_rng(4)
    L = rng.normal(size=(30, 2000))
    sigma2_full = 0.25
    # each half has twice the full-depth noise variance
    h1 = [L + np.sqrt(2 * sigma2_full) * rng.normal(size=L.shape) for _ in range(3)]
    h2 = [L + np.sqrt(2 * sigma2_full) * rng.normal(size=L.shape) for _ in range(3)]
    np.testing.assert_allclose(an.noise_variance(h1, h2).mean(), sigma2_full, rtol=0.05)
