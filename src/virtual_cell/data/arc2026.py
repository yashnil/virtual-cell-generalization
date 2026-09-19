"""The official Arc Virtual Cell Challenge 2026 control bundle.

This module knows the on-disk layout of the official validation bundle
(``context_{A,B,C}.h5ad``, ``gene_names.csv``, ``pert_counts.csv``,
``manifest.json``) and provides **read-only** loading and auditing helpers.

Two properties drive the design:

* The control matrices are large (18,400 cells x 18,533 genes, roughly 1e8
  stored entries each). Every statistic here is computed by streaming row
  chunks of the CSR matrix straight out of HDF5, so the full matrix is never
  densified and never fully resident in memory.
* The contexts A/B/C are *evaluation inputs*. They contain unperturbed control
  cells only; their perturbation responses are hidden. Nothing in this module
  produces perturbation-response training targets.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem
from scipy import sparse

from virtual_cell.data.io import DataIntegrityError

#: Path of the control bundle relative to the repository root.
CONTROLS_SUBDIR = Path("data") / "raw" / "arc2026" / "controls"

CONTEXT_FILE_TEMPLATE = "context_{context}.h5ad"
GENE_NAMES_FILE = "gene_names.csv"
PERT_COUNTS_FILE = "pert_counts.csv"
MANIFEST_FILE = "manifest.json"

GENE_NAMES_COLUMN = "gene_name"

#: Row chunk size used when streaming the count matrices.
DEFAULT_CHUNK_SIZE = 2_000


# --------------------------------------------------------------------------
# bundle layout
# --------------------------------------------------------------------------


def context_path(controls_dir: str | Path, context: str) -> Path:
    """Path of the ``.h5ad`` file holding the control cells of ``context``."""
    return Path(controls_dir) / CONTEXT_FILE_TEMPLATE.format(context=context)


def bundle_paths(controls_dir: str | Path, contexts: Sequence[str]) -> dict[str, Path]:
    """All expected bundle paths, keyed by a short name.

    Keys are ``manifest``, ``gene_names``, ``pert_counts`` and one entry per
    context named ``context_<label>``.
    """
    controls_dir = Path(controls_dir)
    paths = {
        "manifest": controls_dir / MANIFEST_FILE,
        "gene_names": controls_dir / GENE_NAMES_FILE,
        "pert_counts": controls_dir / PERT_COUNTS_FILE,
    }
    for context in contexts:
        paths[f"context_{context}"] = context_path(controls_dir, context)
    return paths


def missing_bundle_files(controls_dir: str | Path, contexts: Sequence[str]) -> list[Path]:
    """Expected bundle files that are not present on disk."""
    return [p for p in bundle_paths(controls_dir, contexts).values() if not p.exists()]


# --------------------------------------------------------------------------
# manifest and side-car tables
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Arc2026Manifest:
    """Typed view of the official ``manifest.json``."""

    season: str
    partition: str
    panel_id: str
    contexts: tuple[str, ...]
    pert_col: str
    context_col: str
    control_label: str
    n_genes: int
    n_constructs: int
    cells_per_pert: int
    per_context: Mapping[str, Mapping[str, int]]
    raw: Mapping[str, object] = field(repr=False, default_factory=dict)

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> Arc2026Manifest:
        missing = [
            key
            for key in (
                "season",
                "partition",
                "panel_id",
                "contexts",
                "pert_col",
                "context_col",
                "control_label",
                "n_genes",
                "n_constructs",
                "cells_per_pert",
                "per_context",
            )
            if key not in payload
        ]
        if missing:
            raise DataIntegrityError(f"manifest.json is missing keys: {missing}")
        return cls(
            season=str(payload["season"]),
            partition=str(payload["partition"]),
            panel_id=str(payload["panel_id"]),
            contexts=tuple(str(c) for c in payload["contexts"]),  # type: ignore[union-attr]
            pert_col=str(payload["pert_col"]),
            context_col=str(payload["context_col"]),
            control_label=str(payload["control_label"]),
            n_genes=int(payload["n_genes"]),  # type: ignore[arg-type]
            n_constructs=int(payload["n_constructs"]),  # type: ignore[arg-type]
            cells_per_pert=int(payload["cells_per_pert"]),  # type: ignore[arg-type]
            per_context=dict(payload["per_context"]),  # type: ignore[arg-type]
            raw=dict(payload),
        )

    def control_cells(self, context: str) -> int:
        """Expected number of control cells for ``context`` per the manifest."""
        return int(self.per_context[context]["control_cells"])

    def n_ntc_ids(self, context: str) -> int:
        """Expected number of distinct non-targeting guide ids for ``context``."""
        return int(self.per_context[context]["n_ntc_ids"])


def load_manifest(controls_dir: str | Path) -> Arc2026Manifest:
    """Read and validate ``manifest.json`` from the control bundle."""
    path = Path(controls_dir) / MANIFEST_FILE
    if not path.exists():
        raise FileNotFoundError(f"Manifest not found: {path}")
    with path.open() as handle:
        payload = json.load(handle)
    return Arc2026Manifest.from_dict(payload)


def load_gene_names(controls_dir: str | Path) -> pd.Index:
    """Read ``gene_names.csv`` as an Index, preserving file order exactly."""
    path = Path(controls_dir) / GENE_NAMES_FILE
    frame = pd.read_csv(path, dtype=str)
    if list(frame.columns) != [GENE_NAMES_COLUMN]:
        raise DataIntegrityError(
            f"{GENE_NAMES_FILE} should have exactly one column {GENE_NAMES_COLUMN!r}, "
            f"found {list(frame.columns)}."
        )
    return pd.Index(frame[GENE_NAMES_COLUMN], name=GENE_NAMES_COLUMN)


def load_pert_counts(controls_dir: str | Path) -> pd.DataFrame:
    """Read ``pert_counts.csv``, preserving file order and columns exactly."""
    path = Path(controls_dir) / PERT_COUNTS_FILE
    return pd.read_csv(path, dtype=str)


#: Gene families reported by :func:`panel_composition`, as regex patterns.
#:
#: The first few are the transcripts that normally dominate a 10x library. The
#: official 2026 panel excludes them, which is why its expression profile is much
#: flatter than raw 10x data; see ``reports/arc2026_controls_audit.md``.
PANEL_GENE_FAMILIES: Mapping[str, str] = {
    "cytoplasmic_ribosomal": r"^RP[LS]\d",
    "mitochondrial_ribosomal": r"^MRP[LS]",
    "mitochondrial_rrna": r"^MT-RNR",
    "mitochondrial_protein_coding": r"^MT-(?!RNR)",
    "histone": r"^HIST",
    "lincrna": r"^LINC",
    "clone_based": r"^A[CLP]\d{6}\.\d",
    "antisense": r"-AS\d$",
}


def panel_composition(
    gene_names: pd.Index, *, families: Mapping[str, str] = PANEL_GENE_FAMILIES
) -> pd.Series:
    """Count how many panel genes fall into each named gene family."""
    names = [str(g) for g in gene_names]
    counts = {
        family: sum(1 for name in names if re.search(pattern, name))
        for family, pattern in families.items()
    }
    counts["total"] = len(names)
    return pd.Series(counts, name="n_genes")


def sha256sum(path: str | Path, *, block_size: int = 1 << 20) -> str:
    """SHA-256 of a file, read in blocks so large matrices stay out of memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_checksum_file(path: str | Path) -> dict[str, str]:
    """Parse a ``shasum``-style file into ``{relative path: digest}``."""
    entries: dict[str, str] = {}
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        digest, _, name = line.partition(" ")
        entries[name.strip().lstrip("*")] = digest
    return entries


# --------------------------------------------------------------------------
# memory-safe streaming over the stored CSR matrices
# --------------------------------------------------------------------------


def read_obs(path: str | Path) -> pd.DataFrame:
    """Read only the ``obs`` table of an ``.h5ad`` file."""
    with h5py.File(Path(path), "r") as handle:
        return read_elem(handle["obs"])


def read_var_names(path: str | Path) -> pd.Index:
    """Read only the ``var`` index (gene names) of an ``.h5ad`` file."""
    with h5py.File(Path(path), "r") as handle:
        return pd.Index(read_elem(handle["var"]).index)


def read_shape(path: str | Path) -> tuple[int, int]:
    """Read the stored shape of ``X`` without loading any data."""
    with h5py.File(Path(path), "r") as handle:
        group = handle["X"]
        if isinstance(group, h5py.Dataset):
            return tuple(int(n) for n in group.shape)  # type: ignore[return-value]
        return tuple(int(n) for n in group.attrs["shape"])  # type: ignore[return-value]


def x_encoding(path: str | Path) -> str:
    """Encoding type of ``X`` (``csr_matrix``, ``csc_matrix`` or ``array``)."""
    with h5py.File(Path(path), "r") as handle:
        group = handle["X"]
        if isinstance(group, h5py.Dataset):
            return "array"
        return str(group.attrs.get("encoding-type", "unknown"))


def x_dtype(path: str | Path) -> np.dtype:
    """dtype of the stored values of ``X``."""
    with h5py.File(Path(path), "r") as handle:
        group = handle["X"]
        dataset = group if isinstance(group, h5py.Dataset) else group["data"]
        return np.dtype(dataset.dtype)


def stream_row_chunks(
    path: str | Path, *, chunk_size: int = DEFAULT_CHUNK_SIZE
) -> Iterator[tuple[int, sparse.csr_matrix]]:
    """Yield ``(start_row, chunk)`` CSR blocks of ``X`` straight from HDF5.

    Only ``chunk_size`` rows are resident at a time, so the full matrix is never
    materialised. Raises :class:`DataIntegrityError` if ``X`` is not CSR.
    """
    path = Path(path)
    with h5py.File(path, "r") as handle:
        group = handle["X"]
        if isinstance(group, h5py.Dataset) or group.attrs.get("encoding-type") != "csr_matrix":
            raise DataIntegrityError(
                f"{path.name}: expected a CSR-encoded X, found {x_encoding(path)!r}."
            )
        n_rows, n_cols = (int(n) for n in group.attrs["shape"])
        indptr = group["indptr"][:]
        data = group["data"]
        indices = group["indices"]
        for start in range(0, n_rows, chunk_size):
            stop = min(start + chunk_size, n_rows)
            lo, hi = int(indptr[start]), int(indptr[stop])
            chunk = sparse.csr_matrix(
                (
                    data[lo:hi],
                    indices[lo:hi],
                    indptr[start : stop + 1] - indptr[start],
                ),
                shape=(stop - start, n_cols),
            )
            yield start, chunk


# --------------------------------------------------------------------------
# per-context audit
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class ContextAudit:
    """Descriptive audit of one official control context.

    ``raw_mean`` and ``normalized_mean`` are per-gene pseudobulk profiles: the
    mean of raw counts and the mean of ``log1p(1e4 * x / library_size)``
    respectively, matching :func:`virtual_cell.preprocessing.pseudobulk.mean_expression`.
    """

    context: str
    path: Path
    n_cells: int
    n_genes: int
    x_encoding: str
    x_dtype: str
    n_stored: int
    n_nonzero: int
    n_explicit_zeros: int
    sparsity: float
    min_value: float
    max_value: float
    all_finite: bool
    all_nonnegative: bool
    all_integer: bool
    obs_columns: tuple[str, ...]
    var_columns: tuple[str, ...]
    obs_names_unique: bool
    var_names_unique: bool
    n_duplicate_obs_names: int
    n_duplicate_var_names: int
    cell_ids: pd.Index = field(repr=False)
    library_sizes: np.ndarray = field(repr=False)
    genes_detected: np.ndarray = field(repr=False)
    raw_mean: pd.Series = field(repr=False)
    normalized_mean: pd.Series = field(repr=False)
    obs_value_counts: Mapping[str, pd.Series] = field(repr=False, default_factory=dict)

    @property
    def is_sparse(self) -> bool:
        return self.x_encoding in {"csr_matrix", "csc_matrix"}

    def stats_row(self) -> dict[str, object]:
        """Flat summary suitable for a one-row-per-context table."""
        lib, det = self.library_sizes, self.genes_detected
        return {
            "context": self.context,
            "n_cells": self.n_cells,
            "n_genes": self.n_genes,
            "x_encoding": self.x_encoding,
            "x_dtype": self.x_dtype,
            "sparsity": self.sparsity,
            "n_stored": self.n_stored,
            "n_nonzero": self.n_nonzero,
            "n_explicit_zeros": self.n_explicit_zeros,
            "min_count": self.min_value,
            "max_count": self.max_value,
            "all_finite": self.all_finite,
            "all_nonnegative": self.all_nonnegative,
            "all_integer": self.all_integer,
            "lib_min": float(lib.min()),
            "lib_median": float(np.median(lib)),
            "lib_mean": float(lib.mean()),
            "lib_max": float(lib.max()),
            "genes_detected_min": float(det.min()),
            "genes_detected_median": float(np.median(det)),
            "genes_detected_mean": float(det.mean()),
            "genes_detected_max": float(det.max()),
        }


def audit_context_file(
    path: str | Path,
    *,
    expected_context: str | None = None,
    context_col: str = "context",
    target_sum: float = 1e4,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> ContextAudit:
    """Audit one control ``.h5ad`` in a single streaming pass.

    The file is opened read-only and is never modified. If ``expected_context``
    is given, the file's context label must equal it exactly.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Context file not found: {path}")

    obs = read_obs(path)
    with h5py.File(path, "r") as handle:
        var = read_elem(handle["var"])
    var_names = pd.Index(var.index)
    n_cells, n_genes = read_shape(path)

    if context_col not in obs.columns:
        raise DataIntegrityError(
            f"{path.name}: obs is missing the context column {context_col!r}; "
            f"columns are {list(obs.columns)}."
        )
    labels = obs[context_col].astype(str).unique().tolist()
    if len(labels) != 1:
        raise DataIntegrityError(
            f"{path.name}: expected exactly one context label, found {labels}."
        )
    label = labels[0]
    if expected_context is not None and label != expected_context:
        raise DataIntegrityError(
            f"{path.name}: context label mismatch, file says {label!r}, "
            f"expected {expected_context!r}."
        )

    library = np.zeros(n_cells, dtype=np.float64)
    detected = np.zeros(n_cells, dtype=np.int64)
    gene_sum = np.zeros(n_genes, dtype=np.float64)
    normalized_sum = np.zeros(n_genes, dtype=np.float64)
    n_stored = 0
    n_nonzero = 0
    min_stored = np.inf
    max_stored = -np.inf
    all_finite = True
    all_nonnegative = True
    all_integer = True

    for start, chunk in stream_row_chunks(path, chunk_size=chunk_size):
        stop = start + chunk.shape[0]
        values = chunk.data
        n_stored += values.size
        if values.size:
            finite = np.isfinite(values)
            if not finite.all():
                all_finite = False
                values = values[finite]
            if values.size:
                chunk_min = float(values.min())
                chunk_max = float(values.max())
                min_stored = min(min_stored, chunk_min)
                max_stored = max(max_stored, chunk_max)
                if chunk_min < 0:
                    all_nonnegative = False
                if all_integer and not np.array_equal(values, np.round(values)):
                    all_integer = False
            n_nonzero += int(np.count_nonzero(chunk.data))

        library[start:stop] = np.asarray(chunk.sum(axis=1)).ravel()
        counted = chunk.copy()
        counted.eliminate_zeros()
        detected[start:stop] = np.diff(counted.indptr)
        gene_sum += np.asarray(chunk.sum(axis=0)).ravel()

        lib_chunk = library[start:stop]
        scale = np.zeros_like(lib_chunk)
        nonzero_rows = lib_chunk > 0
        scale[nonzero_rows] = target_sum / lib_chunk[nonzero_rows]
        normalized = sparse.diags(scale) @ chunk.astype(np.float64)
        normalized = sparse.csr_matrix(normalized)
        normalized.data = np.log1p(normalized.data)
        normalized_sum += np.asarray(normalized.sum(axis=0)).ravel()

    n_total = n_cells * n_genes
    min_value = float(min_stored) if n_stored == n_total else min(0.0, float(min_stored))
    dup_obs = int(obs.index.duplicated().sum())
    dup_var = int(var_names.duplicated().sum())

    return ContextAudit(
        context=label,
        path=path,
        n_cells=n_cells,
        n_genes=n_genes,
        x_encoding=x_encoding(path),
        x_dtype=str(x_dtype(path)),
        n_stored=n_stored,
        n_nonzero=n_nonzero,
        n_explicit_zeros=n_stored - n_nonzero,
        sparsity=1.0 - n_nonzero / n_total,
        min_value=min_value,
        max_value=float(max_stored),
        all_finite=all_finite,
        all_nonnegative=all_nonnegative,
        all_integer=all_integer,
        obs_columns=tuple(obs.columns),
        var_columns=tuple(var.columns),
        obs_names_unique=bool(obs.index.is_unique),
        var_names_unique=bool(var_names.is_unique),
        n_duplicate_obs_names=dup_obs,
        n_duplicate_var_names=dup_var,
        cell_ids=pd.Index(obs.index),
        library_sizes=library,
        genes_detected=detected,
        raw_mean=pd.Series(gene_sum / n_cells, index=var_names, name=label),
        normalized_mean=pd.Series(normalized_sum / n_cells, index=var_names, name=label),
        obs_value_counts={
            column: obs[column].value_counts() for column in obs.columns
        },
    )


def audit_summary_table(audits: Mapping[str, ContextAudit]) -> pd.DataFrame:
    """One-row-per-context table of the core audit statistics."""
    frame = pd.DataFrame([audit.stats_row() for audit in audits.values()])
    return frame.set_index("context") if not frame.empty else frame


def basal_mean_table(
    audits: Mapping[str, ContextAudit], *, normalized: bool = True
) -> pd.DataFrame:
    """Genes x contexts pseudobulk table assembled from per-context audits.

    Requires identical gene order across contexts; that invariant is checked by
    :func:`check_cross_context_invariants`.
    """
    columns = {
        label: (audit.normalized_mean if normalized else audit.raw_mean)
        for label, audit in audits.items()
    }
    return pd.DataFrame(columns)


#: Quantiles reported for the per-cell depth distributions.
DEPTH_QUANTILES = (0.005, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)

#: Library-size thresholds used to report low-depth cell fractions.
LOW_DEPTH_THRESHOLDS = (1_000, 2_000, 5_000)


def depth_quantile_table(
    audits: Mapping[str, ContextAudit],
    *,
    quantiles: Sequence[float] = DEPTH_QUANTILES,
    low_depth_thresholds: Sequence[int] = LOW_DEPTH_THRESHOLDS,
) -> pd.DataFrame:
    """Per-context quantiles of library size and genes detected.

    Also reports the fraction of cells below each library-size threshold, which
    is how a context with an unusually shallow tail shows up.
    """
    rows = {}
    for label, audit in audits.items():
        lib = audit.library_sizes
        det = audit.genes_detected
        row: dict[str, float] = {}
        for q in quantiles:
            row[f"lib_q{q:g}"] = float(np.quantile(lib, q))
        for q in quantiles:
            row[f"genes_detected_q{q:g}"] = float(np.quantile(det, q))
        row["lib_cv"] = float(lib.std() / lib.mean()) if lib.mean() else float("nan")
        for threshold in low_depth_thresholds:
            row[f"n_cells_lib_lt_{threshold}"] = int((lib < threshold).sum())
            row[f"frac_cells_lib_lt_{threshold}"] = float((lib < threshold).mean())
        rows[label] = row
    return pd.DataFrame(rows).T.rename_axis("context")


def top_basal_difference_genes(table: pd.DataFrame, n: int = 25) -> pd.DataFrame:
    """Genes with the largest spread of basal expression across contexts.

    The statistic is purely descriptive: ``max - min`` of the per-context
    pseudobulk value, with the per-context values kept alongside it. No test,
    ranking model, or identity inference is involved.
    """
    spread = table.max(axis=1) - table.min(axis=1)
    out = table.copy()
    out["spread"] = spread
    out["argmax_context"] = table.idxmax(axis=1)
    out["argmin_context"] = table.idxmin(axis=1)
    return out.sort_values("spread", ascending=False).head(n)


# --------------------------------------------------------------------------
# cross-context invariants
# --------------------------------------------------------------------------


def same_labels(left: pd.Index, right: pd.Index) -> bool:
    """Exact element-wise, order-sensitive comparison of two label Indexes.

    ``pandas.Index.equals`` returns False when one Index has the ``string``
    extension dtype and the other plain ``object`` dtype even if every label is
    identical, which is exactly the case between an ``.h5ad`` ``var`` index and a
    CSV column. Both sides are cast to ``object`` so only the labels and their
    order are compared. No whitespace, case, or unicode normalisation is applied.
    """
    left_values = np.asarray(left, dtype=object)
    right_values = np.asarray(right, dtype=object)
    if left_values.shape != right_values.shape:
        return False
    return bool(np.array_equal(left_values, right_values))


@dataclass(frozen=True)
class InvariantCheck:
    """Outcome of one named invariant check."""

    name: str
    passed: bool
    detail: str = ""


def check_cross_context_invariants(
    audits: Mapping[str, ContextAudit],
    *,
    manifest: Arc2026Manifest | None = None,
    gene_names: pd.Index | None = None,
) -> list[InvariantCheck]:
    """Check every invariant that must hold across the official contexts.

    Returns one :class:`InvariantCheck` per invariant. Nothing is raised, so the
    caller can report all failures at once.
    """
    checks: list[InvariantCheck] = []
    labels = list(audits)
    reference = audits[labels[0]]
    ref_genes = pd.Index(reference.raw_mean.index)

    if manifest is not None:
        checks.append(
            InvariantCheck(
                "context labels match the manifest",
                labels == list(manifest.contexts),
                f"files report {labels}, manifest lists {list(manifest.contexts)}",
            )
        )

    for label, audit in audits.items():
        checks.append(
            InvariantCheck(
                f"context {label}: label stored exactly as expected",
                audit.context == label,
                f"stored label {audit.context!r}",
            )
        )
        checks.append(
            InvariantCheck(
                f"context {label}: obs_names unique",
                audit.obs_names_unique,
                f"{audit.n_duplicate_obs_names} duplicate cell ids",
            )
        )
        checks.append(
            InvariantCheck(
                f"context {label}: var_names unique",
                audit.var_names_unique,
                f"{audit.n_duplicate_var_names} duplicate gene names",
            )
        )
        checks.append(
            InvariantCheck(
                f"context {label}: counts finite, non-negative and integer-valued",
                audit.all_finite and audit.all_nonnegative and audit.all_integer,
                f"finite={audit.all_finite} nonneg={audit.all_nonnegative} "
                f"integer={audit.all_integer} min={audit.min_value} max={audit.max_value}",
            )
        )
        genes = pd.Index(audit.raw_mean.index)
        checks.append(
            InvariantCheck(
                f"context {label}: gene set identical to context {labels[0]}",
                set(genes) == set(ref_genes),
                f"{len(set(genes) ^ set(ref_genes))} genes differ",
            )
        )
        checks.append(
            InvariantCheck(
                f"context {label}: gene order identical to context {labels[0]}",
                same_labels(genes, ref_genes),
                "",
            )
        )
        if gene_names is not None:
            checks.append(
                InvariantCheck(
                    f"context {label}: var_names match gene_names.csv exactly in order",
                    same_labels(genes, pd.Index(gene_names)),
                    f"{len(genes)} vs {len(gene_names)} genes",
                )
            )
        if manifest is not None:
            expected_cells = manifest.control_cells(label)
            checks.append(
                InvariantCheck(
                    f"context {label}: cell count matches manifest",
                    audit.n_cells == expected_cells,
                    f"{audit.n_cells} vs {expected_cells}",
                )
            )
            checks.append(
                InvariantCheck(
                    f"context {label}: gene count matches manifest",
                    audit.n_genes == manifest.n_genes,
                    f"{audit.n_genes} vs {manifest.n_genes}",
                )
            )
            pert_values = set(
                audit.obs_value_counts.get(manifest.pert_col, pd.Series(dtype=int)).index
            )
            checks.append(
                InvariantCheck(
                    f"context {label}: only control cells "
                    f"({manifest.pert_col} == {manifest.control_label!r})",
                    pert_values == {manifest.control_label},
                    f"observed {sorted(pert_values)}",
                )
            )
            ntc = audit.obs_value_counts.get("ntc_id")
            if ntc is not None:
                checks.append(
                    InvariantCheck(
                        f"context {label}: non-targeting guide count matches manifest",
                        len(ntc) == manifest.n_ntc_ids(label),
                        f"{len(ntc)} vs {manifest.n_ntc_ids(label)}",
                    )
                )
                checks.append(
                    InvariantCheck(
                        f"context {label}: every guide has cells_per_pert cells",
                        bool((ntc == manifest.cells_per_pert).all()),
                        f"min={int(ntc.min())} max={int(ntc.max())} "
                        f"expected={manifest.cells_per_pert}",
                    )
                )
            expected_total = (
                manifest.control_cells(label) + manifest.n_constructs * manifest.cells_per_pert
            )
            ground_truth = int(manifest.per_context[label].get("ground_truth_cells", -1))
            checks.append(
                InvariantCheck(
                    f"context {label}: manifest cell arithmetic is self-consistent",
                    ground_truth == expected_total,
                    f"ground_truth_cells={ground_truth}, "
                    f"control_cells + n_constructs * cells_per_pert={expected_total}",
                )
            )

    guide_panels = {
        label: tuple(audit.obs_value_counts["ntc_id"].sort_index().index)
        for label, audit in audits.items()
        if "ntc_id" in audit.obs_value_counts
    }
    if len(guide_panels) == len(labels) and labels:
        reference_panel = guide_panels[labels[0]]
        checks.append(
            InvariantCheck(
                "non-targeting guide panel identical across contexts",
                all(panel == reference_panel for panel in guide_panels.values()),
                f"{len(reference_panel)} guides",
            )
        )

    for i, left in enumerate(labels):
        for right in labels[i + 1 :]:
            shared = set(audits[left].cell_ids) & set(audits[right].cell_ids)
            checks.append(
                InvariantCheck(
                    f"cell ids disjoint between {left} and {right}",
                    not shared,
                    f"{len(shared)} shared cell ids",
                )
            )
    return checks


def cell_id_overlaps(audits: Mapping[str, ContextAudit]) -> pd.DataFrame:
    """Pairwise count of shared cell identifiers between contexts."""
    labels = list(audits)
    frame = pd.DataFrame(0, index=labels, columns=labels, dtype=int)
    for i, left in enumerate(labels):
        for right in labels[i:]:
            shared = len(set(audits[left].cell_ids) & set(audits[right].cell_ids))
            frame.loc[left, right] = shared
            frame.loc[right, left] = shared
    return frame


def subsample_context(
    path: str | Path, *, n_cells: int, seed: int = 0, genes: pd.Index | None = None
) -> ad.AnnData:
    """Load a random subset of cells from a control context.

    Used for exploratory embeddings where loading all 18,400 cells x 18,533
    genes at once is wasteful. Restrict to ``genes`` to bound memory further.
    """
    path = Path(path)
    obs = read_obs(path)
    rng = np.random.default_rng(seed)
    n_cells = min(n_cells, len(obs))
    rows = np.sort(rng.choice(len(obs), size=n_cells, replace=False))

    var_names = read_var_names(path)
    keep = None
    if genes is not None:
        keep = var_names.get_indexer(pd.Index(genes))
        if (keep < 0).any():
            raise DataIntegrityError("Requested genes are not all present in the file.")

    blocks = []
    for start, chunk in stream_row_chunks(path):
        local = [r - start for r in rows if start <= r < start + chunk.shape[0]]
        if not local:
            continue
        block = chunk[local]
        if keep is not None:
            block = block[:, keep]
        blocks.append(block)

    X = sparse.vstack(blocks, format="csr") if blocks else sparse.csr_matrix((0, len(var_names)))
    adata = ad.AnnData(
        X=X,
        obs=obs.iloc[rows].copy(),
        var=pd.DataFrame(index=var_names[keep] if keep is not None else var_names),
    )
    return adata
