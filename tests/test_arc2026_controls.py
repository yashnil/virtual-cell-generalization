"""Invariants of the official Arc Virtual Cell Challenge 2026 control bundle.

Every expectation here is derived from the downloaded files themselves or from
``manifest.json``; nothing is hard-coded from the challenge website. The bundle
is git-ignored raw data, so the whole module is skipped when it is absent.

These tests are read-only. They never write to ``data/`` and they never treat
A/B/C as perturbation-response training data.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from virtual_cell.data import arc2026
from virtual_cell.data.io import DataIntegrityError

REPO_ROOT = Path(__file__).resolve().parents[1]
CONTROLS_DIR = REPO_ROOT / arc2026.CONTROLS_SUBDIR

pytestmark = pytest.mark.skipif(
    not (CONTROLS_DIR / arc2026.MANIFEST_FILE).exists(),
    reason=f"Official Arc 2026 control bundle not present under {CONTROLS_DIR}",
)


@pytest.fixture(scope="module")
def manifest() -> arc2026.Arc2026Manifest:
    return arc2026.load_manifest(CONTROLS_DIR)


@pytest.fixture(scope="module")
def gene_names() -> pd.Index:
    return arc2026.load_gene_names(CONTROLS_DIR)


@pytest.fixture(scope="module")
def pert_counts() -> pd.DataFrame:
    return arc2026.load_pert_counts(CONTROLS_DIR)


@pytest.fixture(scope="module")
def contexts(manifest) -> tuple[str, ...]:
    return manifest.contexts


@pytest.fixture(scope="module")
def var_names(contexts) -> dict[str, pd.Index]:
    return {
        context: arc2026.read_var_names(arc2026.context_path(CONTROLS_DIR, context))
        for context in contexts
    }


@pytest.fixture(scope="module")
def obs_tables(contexts) -> dict[str, pd.DataFrame]:
    return {
        context: arc2026.read_obs(arc2026.context_path(CONTROLS_DIR, context))
        for context in contexts
    }


# --- bundle completeness --------------------------------------------------


def test_all_expected_official_files_exist(manifest):
    assert arc2026.missing_bundle_files(CONTROLS_DIR, manifest.contexts) == []


def test_manifest_declares_three_contexts(manifest):
    assert len(manifest.contexts) == len(set(manifest.contexts))
    assert manifest.n_genes > 0
    assert manifest.n_constructs > 0
    assert manifest.cells_per_pert > 0
    assert set(manifest.per_context) == set(manifest.contexts)


def test_manifest_matches_its_provenance_copy(manifest):
    """The working copy and the recorded provenance copy must not drift apart."""
    provenance = REPO_ROOT / "data" / "provenance" / "arc2026_manifest.json"
    if not provenance.exists():
        pytest.skip("No provenance copy of the manifest recorded.")
    recorded = arc2026.Arc2026Manifest.from_dict(json.loads(provenance.read_text()))
    assert dict(manifest.raw) == dict(recorded.raw)


# --- side-car tables ------------------------------------------------------


def test_gene_names_csv_has_no_duplicates(gene_names):
    assert gene_names.is_unique
    assert not gene_names.isna().any()


def test_gene_names_count_matches_manifest(gene_names, manifest):
    assert len(gene_names) == manifest.n_genes


def test_pert_counts_schema_and_uniqueness(pert_counts, manifest):
    assert manifest.pert_col in pert_counts.columns
    targets = pert_counts[manifest.pert_col]
    assert targets.is_unique
    assert not targets.isna().any()
    assert len(pert_counts) == manifest.n_constructs


def test_pert_counts_excludes_the_control_label(pert_counts, manifest):
    assert manifest.control_label not in set(pert_counts[manifest.pert_col])


# --- per-context structure ------------------------------------------------


def test_context_files_load_and_report_manifest_shapes(contexts, manifest):
    for context in contexts:
        n_cells, n_genes = arc2026.read_shape(arc2026.context_path(CONTROLS_DIR, context))
        assert n_cells == manifest.control_cells(context)
        assert n_genes == manifest.n_genes


def test_context_labels_are_exactly_the_manifest_labels(obs_tables, manifest):
    for context, obs in obs_tables.items():
        labels = obs[manifest.context_col].astype(str).unique().tolist()
        assert labels == [context], f"context {context} carries labels {labels}"


def test_contexts_contain_control_cells_only(obs_tables, manifest):
    for context, obs in obs_tables.items():
        observed = set(obs[manifest.pert_col].astype(str))
        assert observed == {manifest.control_label}, f"context {context}: {sorted(observed)}"


def test_no_duplicate_obs_names(obs_tables):
    for context, obs in obs_tables.items():
        assert obs.index.is_unique, f"context {context} has duplicate cell ids"


def test_no_duplicate_var_names(var_names):
    for context, names in var_names.items():
        assert names.is_unique, f"context {context} has duplicate gene names"


def test_cell_ids_are_disjoint_across_contexts(obs_tables):
    seen: dict[str, str] = {}
    for context, obs in obs_tables.items():
        for cell_id in obs.index:
            assert cell_id not in seen, (
                f"cell id {cell_id!r} appears in both {seen[cell_id]} and {context}"
            )
            seen[cell_id] = context


def test_x_is_sparse_csr(contexts):
    for context in contexts:
        assert arc2026.x_encoding(arc2026.context_path(CONTROLS_DIR, context)) == "csr_matrix"


# --- gene-space invariants ------------------------------------------------


def test_gene_sets_are_identical_across_contexts(var_names, contexts):
    reference = set(var_names[contexts[0]])
    for context in contexts[1:]:
        assert set(var_names[context]) == reference


def test_gene_order_is_identical_across_contexts(var_names, contexts):
    reference = var_names[contexts[0]]
    for context in contexts[1:]:
        assert arc2026.same_labels(var_names[context], reference)


def test_var_names_match_gene_names_csv_exactly_in_order(var_names, gene_names):
    for context, names in var_names.items():
        assert arc2026.same_labels(names, gene_names), f"context {context} gene order differs"


def test_same_labels_is_order_and_value_sensitive():
    base = pd.Index(["A", "B", "C"], dtype="string")
    assert arc2026.same_labels(base, pd.Index(["A", "B", "C"], dtype=object))
    assert not arc2026.same_labels(base, pd.Index(["A", "C", "B"], dtype=object))
    assert not arc2026.same_labels(base, pd.Index(["A", "b", "C"], dtype=object))
    assert not arc2026.same_labels(base, pd.Index(["A", "B ", "C"], dtype=object))
    assert not arc2026.same_labels(base, pd.Index(["A", "B"], dtype=object))


# --- raw count invariants -------------------------------------------------


@pytest.mark.parametrize("context", ["A", "B", "C"])
def test_raw_counts_are_finite_nonnegative_and_integer_valued(context, manifest):
    if context not in manifest.contexts:
        pytest.skip(f"Context {context} is not in this manifest.")
    path = arc2026.context_path(CONTROLS_DIR, context)
    seen_any = False
    for _, chunk in arc2026.stream_row_chunks(path):
        values = chunk.data
        if not values.size:
            continue
        seen_any = True
        assert np.all(np.isfinite(values)), f"context {context}: non-finite counts"
        assert values.min() >= 0, f"context {context}: negative counts"
        assert np.array_equal(values, np.round(values)), f"context {context}: non-integer counts"
    assert seen_any, f"context {context}: no stored counts read"


def test_stream_row_chunks_covers_every_cell(manifest):
    context = manifest.contexts[0]
    path = arc2026.context_path(CONTROLS_DIR, context)
    n_cells, n_genes = arc2026.read_shape(path)
    seen = 0
    for start, chunk in arc2026.stream_row_chunks(path, chunk_size=4096):
        assert start == seen
        assert chunk.shape[1] == n_genes
        seen += chunk.shape[0]
    assert seen == n_cells


def test_audit_rejects_a_wrong_expected_context(manifest):
    right, wrong = manifest.contexts[0], manifest.contexts[1]
    with pytest.raises(DataIntegrityError, match="context label mismatch"):
        arc2026.audit_context_file(
            arc2026.context_path(CONTROLS_DIR, right), expected_context=wrong
        )


# --- guide panel and manifest arithmetic ----------------------------------


def test_every_guide_has_cells_per_pert_cells(obs_tables, manifest):
    for context, obs in obs_tables.items():
        counts = obs["ntc_id"].value_counts()
        assert len(counts) == manifest.n_ntc_ids(context)
        assert set(counts) == {manifest.cells_per_pert}, f"context {context}: {sorted(set(counts))}"


def test_guide_panel_is_identical_across_contexts(obs_tables, contexts):
    reference = sorted(obs_tables[contexts[0]]["ntc_id"].astype(str).unique())
    for context in contexts[1:]:
        assert sorted(obs_tables[context]["ntc_id"].astype(str).unique()) == reference


def test_manifest_cell_arithmetic_is_self_consistent(manifest):
    for context in manifest.contexts:
        per_context = manifest.per_context[context]
        if "ground_truth_cells" not in per_context:
            pytest.skip("Manifest does not report ground_truth_cells.")
        expected = manifest.control_cells(context) + manifest.n_constructs * manifest.cells_per_pert
        assert int(per_context["ground_truth_cells"]) == expected
        assert manifest.control_cells(context) == (
            manifest.n_ntc_ids(context) * manifest.cells_per_pert
        )


# --- provenance -----------------------------------------------------------


def test_official_files_match_their_recorded_checksums(manifest):
    """Guards against the raw bundle silently changing under the audit."""
    checksums = REPO_ROOT / "data" / "provenance" / "arc2026_controls_sha256.txt"
    if not checksums.exists():
        pytest.skip("No recorded checksums for the control bundle.")
    recorded = arc2026.parse_checksum_file(checksums)
    assert recorded, "Checksum file parsed to no entries."
    for relative, digest in recorded.items():
        path = REPO_ROOT / relative
        assert path.exists(), f"Recorded file is missing: {relative}"
        assert arc2026.sha256sum(path) == digest, f"{relative} does not match its recorded digest"


def test_panel_excludes_the_highest_abundance_transcript_families(gene_names):
    """Verified from the downloaded panel; a change here means the panel changed.

    The official panel carries no cytoplasmic or mitochondrial ribosomal protein
    genes and no mitochondrial rRNA, which is why absolute expression levels are
    not comparable with unfiltered public 10x data.
    """
    composition = arc2026.panel_composition(gene_names)
    assert composition["cytoplasmic_ribosomal"] == 0
    assert composition["mitochondrial_ribosomal"] == 0
    assert composition["mitochondrial_rrna"] == 0
    assert composition["mitochondrial_protein_coding"] > 0
    assert composition["total"] == len(gene_names)


def test_all_target_genes_are_in_the_expression_panel(pert_counts, gene_names, manifest):
    """A CRISPRi target must be measurable in the panel it is scored on."""
    missing = sorted(set(pert_counts[manifest.pert_col]) - set(gene_names))
    assert missing == [], f"target genes absent from the panel: {missing}"
