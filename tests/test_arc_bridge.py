"""Checks on the Arc bridge: the panel mapping, the count generators, and the
local reimplementation of the six scored ``vcc2026`` metrics.

The metric tests are transcription checks. The published reference states
several values analytically -- a zero-effect prediction scores exactly 0.5 on
``pds_cosine``, a zero-fold-change prediction scores exactly 1 on
``de_wilcoxon_lfc_nmae``, a submission calling nothing scores 0 on direction
fidelity -- and each of those is asserted here. If this module has misread the
specification, these are where it shows.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from virtual_cell.arc import generate, metrics, panel

PANEL = pd.Index(["A", "B", "C", "D", "E"])


@pytest.fixture
def rng() -> np.random.Generator:
    return np.random.default_rng(20260920)


# --------------------------------------------------------------------------
# panel mapping
# --------------------------------------------------------------------------


def test_panel_map_categorises_every_gene() -> None:
    pm = panel.build_panel_map(PANEL, ["C", "A", "Z"])
    assert pm.n_panel == 5
    counts = pm.counts()
    assert counts[panel.Support.PREDICTED] == 2
    assert counts[panel.Support.UNMEASURED] == 3
    assert counts.sum() == 5
    # "Z" is not in the panel and must not leak in.
    assert list(pm.panel_genes) == list(PANEL)


def test_panel_map_separates_measured_from_responsive() -> None:
    pm = panel.build_panel_map(PANEL, ["A", "B", "C"], responsive=np.array([True, False, True]))
    counts = pm.counts()
    assert counts[panel.Support.PREDICTED] == 2
    assert counts[panel.Support.MEASURED_NO_RESPONSE] == 1
    assert counts[panel.Support.UNMEASURED] == 2


def test_project_response_places_values_and_fills_the_rest() -> None:
    pm = panel.build_panel_map(PANEL, ["C", "A"])
    out = panel.project_response(np.array([7.0, 3.0]), pm, fill=np.nan)
    assert out[0] == 3.0  # "A" is source position 1
    assert out[2] == 7.0  # "C" is source position 0
    assert np.isnan(out[[1, 3, 4]]).all()


def test_project_response_rejects_a_mismatched_width() -> None:
    pm = panel.build_panel_map(PANEL, ["A", "B"])
    with pytest.raises(ValueError, match="expected 2"):
        panel.project_response(np.zeros(3), pm)


def test_project_response_handles_a_stack_of_responses() -> None:
    pm = panel.build_panel_map(PANEL, ["A", "B"])
    out = panel.project_response(np.array([[1.0, 2.0], [3.0, 4.0]]), pm)
    assert out.shape == (2, 5)
    assert out[:, 0].tolist() == [1.0, 3.0]
    assert out[:, 2:].sum() == 0.0


def test_panel_map_rejects_duplicate_panel_genes() -> None:
    with pytest.raises(ValueError, match="unique"):
        panel.build_panel_map(pd.Index(["A", "A"]), ["A"])


# --------------------------------------------------------------------------
# generators
# --------------------------------------------------------------------------


@pytest.fixture
def controls(rng: np.random.Generator) -> np.ndarray:
    """A control block with the sparsity real single-cell counts have.

    Density matters here: the generators' failure modes are about genes a cell
    genuinely failed to catch, and a dense fixture cannot exhibit them.
    """
    shape = (300, 60)
    detected = rng.binomial(1, 0.35, size=shape)
    return rng.poisson(rng.gamma(0.9, 12.0, size=shape) * detected).astype(np.int64)


def test_generators_emit_non_negative_integer_counts(
    controls: np.ndarray, rng: np.random.Generator
) -> None:
    lfc = rng.normal(0.0, 0.3, size=controls.shape[1])
    cpm = controls.mean(axis=0) + 1.0
    for out in (
        generate.resample_controls(controls, 50, rng=rng),
        generate.transport_controls(controls, lfc, 50, rng=rng),
        generate.count_model(controls, cpm, 50, rng=rng),
    ):
        assert out.shape == (50, 60)
        assert np.issubdtype(out.dtype, np.integer)
        assert (out >= 0).all()


def test_transport_preserves_each_drawn_cell_library_size(
    controls: np.ndarray, rng: np.random.Generator
) -> None:
    lfc = rng.normal(0.0, 1.0, size=controls.shape[1])
    out = generate.transport_controls(controls, lfc, 200, rng=rng)
    # Every emitted total must be a total some control cell actually had.
    assert set(out.sum(axis=1)).issubset(set(controls.sum(axis=1).tolist()))


def test_transport_with_no_effect_matches_the_control_composition(
    controls: np.ndarray, rng: np.random.Generator
) -> None:
    out = generate.transport_controls(controls, np.zeros(controls.shape[1]), 3000, rng=rng)
    got = out.sum(axis=0) / out.sum()
    want = controls.sum(axis=0) / controls.sum()
    assert np.abs(got - want).max() < 0.01


def test_transport_moves_composition_in_the_requested_direction(
    controls: np.ndarray, rng: np.random.Generator
) -> None:
    lfc = np.zeros(controls.shape[1])
    lfc[0] = 2.0
    out = generate.transport_controls(controls, lfc, 2000, rng=rng)
    before = controls.sum(axis=0)[0] / controls.sum()
    after = out.sum(axis=0)[0] / out.sum()
    assert after > 3.0 * before


def test_count_model_tracks_the_requested_composition(
    controls: np.ndarray, rng: np.random.Generator
) -> None:
    want = np.linspace(1.0, 10.0, controls.shape[1])
    want = want / want.sum()
    out = generate.count_model(controls, want, 4000, rng=rng)
    got = out.sum(axis=0) / out.sum()
    assert np.corrcoef(got, want)[0, 1] > 0.99


def test_count_model_rejects_an_empty_composition(
    controls: np.ndarray, rng: np.random.Generator
) -> None:
    with pytest.raises(ValueError, match="sums to zero"):
        generate.count_model(controls, np.zeros(controls.shape[1]), 10, rng=rng)


def test_transport_smoothing_protects_the_detection_rate(
    rng: np.random.Generator,
) -> None:
    """Redrawing from a cell's own empirical composition is a second round of
    sampling loss; blending in the pooled composition is what prevents it.

    This builds its own control block rather than using the shared fixture,
    because the effect is a function of counts per gene. Real data carries
    roughly one count per gene per cell, which is where a cell's uncaught genes
    dominate its composition; at the shared fixture's much higher depth the
    pooled composition saturates instead and the test would measure nothing.
    """
    n_genes, n_cells = 2000, 250
    comp = rng.lognormal(0.0, 2.0, size=n_genes)
    comp /= comp.sum()
    ctrl = np.vstack(
        [rng.multinomial(int(lib), comp) for lib in rng.integers(1800, 2600, size=n_cells)]
    ).astype(np.int64)

    want = (ctrl > 0).mean()
    assert 0.2 < want < 0.4, "the fixture must have realistic sparsity"

    zero = np.zeros(n_genes)
    unsmoothed = (generate.transport_controls(ctrl, zero, 600, rng=rng, smoothing=0.0) > 0).mean()
    smoothed = (generate.transport_controls(ctrl, zero, 600, rng=rng, smoothing=0.5) > 0).mean()

    assert unsmoothed < want, "an unsmoothed redraw must lose detected genes"
    assert abs(smoothed - want) < abs(unsmoothed - want)


def test_transport_rejects_an_out_of_range_smoothing(
    controls: np.ndarray, rng: np.random.Generator
) -> None:
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        generate.transport_controls(
            controls, np.zeros(controls.shape[1]), 5, rng=rng, smoothing=1.5
        )


def test_transport_rejects_a_mismatched_effect_width(
    controls: np.ndarray, rng: np.random.Generator
) -> None:
    with pytest.raises(ValueError, match="must have 60 entries"):
        generate.transport_controls(controls, np.zeros(7), 5, rng=rng)


# --------------------------------------------------------------------------
# pds_cosine
# --------------------------------------------------------------------------


def test_pds_is_one_for_a_perfect_prediction(rng: np.random.Generator) -> None:
    real = rng.normal(size=(12, 30))
    assert metrics.pds_cosine(real, real) == pytest.approx(1.0)


def test_pds_is_exactly_one_half_for_a_zero_effect_prediction(
    rng: np.random.Generator,
) -> None:
    """The reference states this analytically: a zero effect ties with everything."""
    real = rng.normal(size=(9, 30))
    got = metrics.pds_cosine(np.zeros_like(real), real)
    assert got == pytest.approx(0.5)


def test_pds_is_exactly_one_half_for_a_single_shared_profile(
    rng: np.random.Generator,
) -> None:
    """A shared profile ranks the panel, but identically for every row, so the
    per-perturbation scores sweep the whole range and their mean -- which is
    what the metric reports -- is exactly one half."""
    real = rng.normal(size=(9, 30))
    shared = np.tile(rng.normal(size=30), (9, 1))
    per_pert = metrics.pds_cosine(shared, real)
    assert per_pert.min() < 0.5 < per_pert.max()
    assert per_pert.mean() == pytest.approx(0.5)


def test_pds_excludes_the_masked_panel_genes(rng: np.random.Generator) -> None:
    """A prediction carrying signal only in excluded coordinates learns nothing."""
    real = rng.normal(size=(8, 10))
    exclude = np.zeros(10, dtype=bool)
    exclude[:4] = True
    pred = np.zeros_like(real)
    pred[:, :4] = real[:, :4]
    assert metrics.pds_cosine(pred, real, exclude=exclude) == pytest.approx(0.5)


def test_pds_needs_at_least_two_perturbations() -> None:
    with pytest.raises(ValueError, match="at least two"):
        metrics.pds_cosine(np.ones((1, 4)), np.ones((1, 4)))


# --------------------------------------------------------------------------
# expression error
# --------------------------------------------------------------------------


def test_expr_mse_is_exactly_one_when_the_prediction_is_the_control(
    rng: np.random.Generator,
) -> None:
    """With no sampling correction, emitting the control makes N and D identical."""
    real = rng.normal(size=(6, 20))
    ctrl = rng.normal(size=20)
    got = metrics.expr_mse_unbiased_capped(
        np.tile(ctrl, (6, 1)),
        real,
        ctrl,
        pred_dispersion=np.zeros(6),
        real_dispersion=np.zeros(6),
        ctrl_dispersion=0.0,
    )
    assert got.value == pytest.approx(1.0)


def test_expr_mse_is_zero_for_a_perfect_noise_free_prediction(
    rng: np.random.Generator,
) -> None:
    real = rng.normal(size=(6, 20))
    ctrl = rng.normal(size=20)
    got = metrics.expr_mse_unbiased_capped(
        real,
        real,
        ctrl,
        pred_dispersion=np.zeros(6),
        real_dispersion=np.zeros(6),
        ctrl_dispersion=0.0,
    )
    assert got.value == pytest.approx(0.0)


def test_expr_mse_rho_never_exceeds_one(rng: np.random.Generator) -> None:
    real = rng.normal(size=(5, 15))
    ctrl = rng.normal(size=15)
    got = metrics.expr_mse_unbiased_capped(
        np.tile(ctrl, (5, 1)),
        real,
        ctrl,
        pred_dispersion=np.full(5, 1e-6),
        real_dispersion=np.full(5, 1e-6),
        ctrl_dispersion=0.0,
    )
    assert 0.0 <= got.rho <= 1.0


def test_expr_mse_excludes_each_perturbations_own_target_gene() -> None:
    real = np.zeros((2, 4))
    real[0, 0] = 100.0
    real[1, 1] = 100.0
    ctrl = np.zeros(4)
    got = metrics.expr_mse_unbiased_capped(
        real,
        real,
        ctrl,
        pred_dispersion=np.zeros(2),
        real_dispersion=np.zeros(2),
        ctrl_dispersion=0.0,
        target_gene=np.array([0, 1]),
    )
    # Every non-excluded coordinate is zero on both sides, so D collapses.
    assert got.denominator.sum() == pytest.approx(0.0)


# --------------------------------------------------------------------------
# the differential-expression family
# --------------------------------------------------------------------------


def _table(sig: np.ndarray, lfc: np.ndarray, pval: np.ndarray | None = None) -> metrics.DETable:
    p_adj = np.where(sig, 0.001, 0.9)
    return metrics.DETable(
        tested=np.ones(lfc.shape[1], dtype=bool),
        pval=p_adj if pval is None else pval,
        p_adj=p_adj,
        lfc=lfc,
    )


def test_sig_jaccard_is_one_when_the_sets_match() -> None:
    sig = np.array([[True, False, True, False]])
    t = _table(sig, np.ones((1, 4)))
    assert metrics.sig_jaccard(t, t)[0] == pytest.approx(1.0)


def test_sig_jaccard_is_zero_when_the_sets_are_disjoint() -> None:
    lfc = np.ones((1, 4))
    real = _table(np.array([[True, True, False, False]]), lfc)
    pred = _table(np.array([[False, False, True, True]]), lfc)
    assert metrics.sig_jaccard(pred, real)[0] == pytest.approx(0.0)


def test_sig_jaccard_defines_an_empty_union_as_one() -> None:
    empty = _table(np.zeros((1, 4), dtype=bool), np.ones((1, 4)))
    assert metrics.sig_jaccard(empty, empty)[0] == pytest.approx(1.0)


def test_direction_fidelity_is_zero_when_the_submission_calls_nothing() -> None:
    lfc = np.ones((1, 6))
    real = _table(np.array([[True] * 4 + [False] * 2]), lfc)
    pred = _table(np.zeros((1, 6), dtype=bool), lfc)
    assert metrics.direction_fidelity_yield(pred, real)[0] == pytest.approx(0.0)


def test_direction_fidelity_is_one_for_a_perfect_call_set() -> None:
    lfc = np.array([[1.0, -1.0, 1.0, -1.0]])
    t = _table(np.array([[True, True, False, False]]), lfc)
    assert metrics.direction_fidelity_yield(t, t)[0] == pytest.approx(1.0)


def test_direction_fidelity_penalises_calling_too_few_genes() -> None:
    lfc = np.array([[1.0, 1.0, 1.0, 1.0]])
    real = _table(np.array([[True, True, True, True]]), lfc)
    pred = _table(np.array([[True, False, False, False]]), lfc)
    # One correct call out of the four the reference found.
    assert metrics.direction_fidelity_yield(pred, real)[0] == pytest.approx(0.25)


def test_direction_fidelity_drops_a_perturbation_with_nothing_on_either_side() -> None:
    lfc = np.ones((1, 4))
    empty = _table(np.zeros((1, 4), dtype=bool), lfc)
    assert np.isnan(metrics.direction_fidelity_yield(empty, empty)[0])


def test_direction_reach_is_one_for_a_perfect_prediction() -> None:
    lfc = np.array([[1.0, -1.0, 2.0, -2.0, 1.5, -1.5, 1.0, -1.0, 2.0, -2.0, 1.0, -1.0]])
    t = _table(np.ones((1, 12), dtype=bool), lfc)
    assert metrics.direction_reach(t, t)[0] == pytest.approx(1.0)


def test_direction_reach_collapses_when_every_direction_is_wrong() -> None:
    lfc = np.array([[1.0] * 12])
    real = _table(np.ones((1, 12), dtype=bool), lfc)
    pred = _table(np.ones((1, 12), dtype=bool), -lfc)
    assert metrics.direction_reach(pred, real)[0] == pytest.approx(0.0)


def test_direction_reach_tolerates_one_miss_only_once_deep_enough() -> None:
    """At a purity floor of 0.9 the shallowest tolerant depth is ten."""
    n = 10
    lfc = np.ones((1, n))
    real = _table(np.ones((1, n), dtype=bool), lfc)
    wrong_first = lfc.copy()
    wrong_first[0, 0] = -1.0
    # Rank the miss first by giving it the smallest p-value.
    pval = np.full((1, n), 0.01)
    pval[0, 0] = 1e-6
    pred = _table(np.ones((1, n), dtype=bool), wrong_first, pval=pval)
    assert metrics.direction_reach(pred, real)[0] == pytest.approx(1.0)


def test_lfc_nmae_is_one_when_the_submission_predicts_no_change() -> None:
    """The reference states this: numerator equals denominator term by term."""
    rng = np.random.default_rng(3)
    lfc = rng.normal(size=(1, 30))
    real = _table(np.ones((1, 30), dtype=bool), lfc)
    pred = _table(np.ones((1, 30), dtype=bool), np.zeros((1, 30)))
    assert metrics.lfc_nmae(pred, real)[0] == pytest.approx(1.0)


def test_lfc_nmae_is_zero_for_a_perfect_prediction() -> None:
    rng = np.random.default_rng(4)
    lfc = rng.normal(size=(1, 30))
    t = _table(np.ones((1, 30), dtype=bool), lfc)
    assert metrics.lfc_nmae(t, t)[0] == pytest.approx(0.0)


def test_lfc_nmae_skips_a_perturbation_whose_gate_is_too_small() -> None:
    lfc = np.ones((1, 20))
    sig = np.zeros((1, 20), dtype=bool)
    sig[0, :5] = True
    real = _table(sig, lfc)
    assert np.isnan(metrics.lfc_nmae(real, real)[0])


# --------------------------------------------------------------------------
# the DE table itself, and score scaling
# --------------------------------------------------------------------------


def test_de_table_gate_reads_the_control_group_only(rng: np.random.Generator) -> None:
    ctrl = np.zeros((40, 6), dtype=np.int64)
    ctrl[:, :3] = rng.poisson(50, size=(40, 3))
    pert = ctrl.copy()
    table = metrics.de_table([pert], ctrl)
    # The three all-zero control genes sit below the 5 CPM floor.
    assert table.tested.tolist() == [True, True, True, False, False, False]
    assert table.lfc.shape == (1, 3)


def test_de_table_recovers_a_planted_up_regulation(rng: np.random.Generator) -> None:
    """The table is read on per-cell CPM, so lifting one gene necessarily pushes
    the rest down. That compositional shift is the metric's real behaviour, not
    an artefact to be normalised away, so it is asserted rather than avoided."""
    ctrl = rng.poisson(40, size=(120, 5)).astype(np.int64)
    pert = ctrl.copy()
    pert[:, 0] = rng.poisson(160, size=120)
    table = metrics.de_table([pert], ctrl)
    assert table.significant[0, 0]
    assert table.lfc[0, 0] > 1.0
    assert table.lfc[0, 0] == table.lfc[0].max()
    assert (table.lfc[0, 1:] < 0.0).all()


def test_de_table_accepts_an_externally_supplied_gate(rng: np.random.Generator) -> None:
    ctrl = rng.poisson(40, size=(30, 5)).astype(np.int64)
    gate = np.array([True, False, True, False, True])
    table = metrics.de_table([ctrl], ctrl, tested=gate)
    assert table.tested.tolist() == gate.tolist()
    assert table.lfc.shape == (1, 3)


def test_scale_score_anchors_the_baseline_at_zero_and_the_replicate_at_one() -> None:
    assert metrics.scale_score(0.5, 0.5, 0.95) == pytest.approx(0.0)
    assert metrics.scale_score(0.95, 0.5, 0.95) == pytest.approx(1.0)
    # Lower-is-better members invert naturally: r sits below b.
    assert metrics.scale_score(1.0, 1.0, 0.4) == pytest.approx(0.0)
    assert metrics.scale_score(0.7, 1.0, 0.4) == pytest.approx(0.5)


def test_scale_score_refuses_a_degenerate_anchor_pair() -> None:
    with pytest.raises(ValueError, match="undefined"):
        metrics.scale_score(0.5, 0.5, 0.5)
