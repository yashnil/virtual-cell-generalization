"""The scPertEval standardized four-context perturbation datasets.

Layout, read-only loading, local validation, and memory-safe pseudobulk for the
four CRISPRi Perturb-seq contexts used by the independent four-context
decomposition (K562, RPE1, HepG2, Jurkat).

Audited specification: ``reports/scperteval_four_context_data_spec.md``.

These files are *not* the Arc challenge data and are never mixed with it. They
are also **not** a Molina & Zhang reproduction fixture: the preprocessing
differs and is documented in the spec.

Every routine here streams CSR row chunks straight out of HDF5. No full cell
matrix is ever densified, and callers are expected to process one context at a
time.
"""

from __future__ import annotations

import hashlib
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

# Reuse the audited CSR streamer rather than duplicating it.
from virtual_cell.data.arc2026 import stream_row_chunks
from virtual_cell.data.io import DataIntegrityError

#: Path of the scPertEval bundle relative to the repository root.
DATA_SUBDIR = Path("data") / "raw" / "scperteval"

FILENAME_TEMPLATE = "{dataset}_processed_complete.h5ad"
BASE_URL = "https://storage.googleapis.com/scperteval/processed"

#: The literal control label used by scPertEval.
CONTROL_LABEL = "control"

#: The obs column holding the perturbation label.
PERTURBATION_KEY = "perturbation"

#: Combination perturbations are joined with this character (none expected here).
COMBINATION_DELIMITER = "+"

DEFAULT_CHUNK_SIZE = 4_000


@dataclass(frozen=True)
class ScPertEvalDataset:
    """One audited scPertEval dataset, with the values recorded before download."""

    name: str
    cell_line: str
    expected_bytes: int
    md5_base64: str
    publication: str
    doi: str
    accession: str

    @property
    def filename(self) -> str:
        return FILENAME_TEMPLATE.format(dataset=self.name)

    @property
    def url(self) -> str:
        return f"{BASE_URL}/{self.filename}"

    def path(self, data_dir: str | Path) -> Path:
        return Path(data_dir) / self.filename


#: The four audited contexts, in the frozen context order used everywhere.
#: Values are from reports/scperteval_four_context_data_spec.md, recorded
#: before any download, and are verified against the local files on load.
DATASETS: tuple[ScPertEvalDataset, ...] = (
    ScPertEvalDataset(
        "replogle22k562",
        "K562",
        2_430_512_332,
        "zzGVsPPynk6Inru+fcMVZQ==",
        "Replogle et al., Cell 185(14):2559-2575.e28 (2022)",
        "10.1016/j.cell.2022.05.013",
        "figshare+ 10.25452/figshare.plus.20029387",
    ),
    ScPertEvalDataset(
        "replogle22rpe1",
        "RPE1",
        1_877_364_555,
        "S2VTPbbJ8p75C8i3aD1l3g==",
        "Replogle et al., Cell 185(14):2559-2575.e28 (2022)",
        "10.1016/j.cell.2022.05.013",
        "figshare+ 10.25452/figshare.plus.20029387",
    ),
    ScPertEvalDataset(
        "nadig25hepg2",
        "HepG2",
        1_236_448_196,
        "AmgRZbNFVwWbOYlZII8yLw==",
        "Nadig et al., Nat. Genet. 57(5):1228-1237 (2025)",
        "10.1038/s41588-025-02169-3",
        "GEO GSE264667",
    ),
    ScPertEvalDataset(
        "nadig25jurkat",
        "Jurkat",
        2_004_474_709,
        "2lpruoqyqbt+9hARAQcAbg==",
        "Nadig et al., Nat. Genet. 57(5):1228-1237 (2025)",
        "10.1038/s41588-025-02169-3",
        "GEO GSE264667",
    ),
)

#: Frozen context order (dataset names) for the four-context design.
CONTEXTS: tuple[str, ...] = tuple(d.name for d in DATASETS)

DATASETS_BY_NAME: Mapping[str, ScPertEvalDataset] = {d.name: d for d in DATASETS}

#: scPertEval's documented normalisation target (normalize_total(target_sum)).
TARGET_SUM = 1e4


def dataset(name: str) -> ScPertEvalDataset:
    if name not in DATASETS_BY_NAME:
        raise KeyError(f"Unknown scPertEval dataset {name!r}; known: {list(DATASETS_BY_NAME)}")
    return DATASETS_BY_NAME[name]


def missing_files(data_dir: str | Path, names: Sequence[str] = CONTEXTS) -> list[Path]:
    """Expected files that are not present on disk."""
    return [p for p in (dataset(n).path(data_dir) for n in names) if not p.exists()]


def incomplete_files(data_dir: str | Path, names: Sequence[str] = CONTEXTS) -> list[Path]:
    """Expected files that are absent **or** not yet at their audited byte size.

    A partially downloaded file exists but is unreadable as HDF5, so presence
    alone is not a usable readiness signal.
    """
    out = []
    for name in names:
        ds = dataset(name)
        path = ds.path(data_dir)
        if not path.exists() or path.stat().st_size != ds.expected_bytes:
            out.append(path)
    return out


# --------------------------------------------------------------------------
# checksums
# --------------------------------------------------------------------------


def file_digests(path: str | Path, *, block_size: int = 1 << 22) -> dict[str, object]:
    """Size, MD5 (base64, as GCS reports it) and SHA-256, in a single read."""
    import base64

    md5 = hashlib.md5()  # noqa: S324 - integrity check against GCS, not security
    sha = hashlib.sha256()
    size = 0
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            md5.update(block)
            sha.update(block)
            size += len(block)
    return {
        "size_bytes": size,
        "md5_base64": base64.b64encode(md5.digest()).decode(),
        "sha256": sha.hexdigest(),
    }


# --------------------------------------------------------------------------
# metadata reads (no expression data touched)
# --------------------------------------------------------------------------


def read_shape(path: str | Path) -> tuple[int, int]:
    with h5py.File(Path(path), "r") as f:
        g = f["X"]
        if isinstance(g, h5py.Dataset):
            return tuple(int(v) for v in g.shape)  # type: ignore[return-value]
        return tuple(int(v) for v in g.attrs["shape"])  # type: ignore[return-value]


def x_encoding(path: str | Path) -> str:
    with h5py.File(Path(path), "r") as f:
        g = f["X"]
        return "array" if isinstance(g, h5py.Dataset) else str(g.attrs.get("encoding-type", "?"))


def x_dtype(path: str | Path) -> str:
    with h5py.File(Path(path), "r") as f:
        g = f["X"]
        return str((g if isinstance(g, h5py.Dataset) else g["data"]).dtype)


def _decode(values: Iterable) -> list[str]:
    return [v.decode() if isinstance(v, bytes) else str(v) for v in values]


def read_var_names(path: str | Path) -> pd.Index:
    with h5py.File(Path(path), "r") as f:
        var = f["var"]
        node = var[var.attrs.get("_index", "_index")]
        values = node["values"] if isinstance(node, h5py.Group) else node
        return pd.Index(_decode(values[:]), name="gene_name")


def read_perturbation_labels(path: str | Path) -> tuple[list[str], np.ndarray]:
    """Return ``(categories, codes)`` for ``obs['perturbation']``.

    Only the categorical encoding is read; the per-cell codes array is small
    (one int per cell), so this never touches expression data.
    """
    with h5py.File(Path(path), "r") as f:
        node = f["obs"][PERTURBATION_KEY]
        if isinstance(node, h5py.Group) and "categories" in node:
            return _decode(node["categories"][:]), node["codes"][:]
        values = node["values"] if isinstance(node, h5py.Group) else node
        labels = np.asarray(_decode(values[:]))
        cats = sorted(set(labels.tolist()))
        lookup = {c: i for i, c in enumerate(cats)}
        return cats, np.array([lookup[v] for v in labels], dtype=np.int32)


def cell_labels(path: str | Path) -> np.ndarray:
    """Per-cell perturbation label as a string array."""
    cats, codes = read_perturbation_labels(path)
    return np.asarray(cats, dtype=object)[codes]


def structure_report(path: str | Path) -> dict[str, object]:
    """Everything a local validation needs, read without touching ``X`` values."""
    path = Path(path)
    with h5py.File(path, "r") as f:
        obs_cols = [str(c) for c in f["obs"].attrs.get("column-order", [])]
        var_cols = [str(c) for c in f["var"].attrs.get("column-order", [])]
        extras = {
            k: list(f[k].keys())
            for k in ("layers", "obsm", "varm", "obsp", "varp", "uns")
            if k in f
        }
    n_cells, n_genes = read_shape(path)
    cats, codes = read_perturbation_labels(path)
    counts = np.bincount(codes, minlength=len(cats))
    perts = [c for c in cats if c != CONTROL_LABEL]
    pert_counts = np.array([counts[cats.index(p)] for p in perts])
    return {
        "path": str(path),
        "n_cells": n_cells,
        "n_genes": n_genes,
        "x_encoding": x_encoding(path),
        "x_dtype": x_dtype(path),
        "obs_columns": obs_cols,
        "var_columns": var_cols,
        "extras": extras,
        "n_labels": len(cats),
        "n_perturbations": len(perts),
        "has_control": CONTROL_LABEL in cats,
        "n_control_cells": int(counts[cats.index(CONTROL_LABEL)]) if CONTROL_LABEL in cats else 0,
        "n_perturbed_cells": int(pert_counts.sum()),
        "min_cells_per_pert": int(pert_counts.min()) if len(pert_counts) else 0,
        "median_cells_per_pert": float(np.median(pert_counts))
        if len(pert_counts)
        else float("nan"),
        "max_cells_per_pert": int(pert_counts.max()) if len(pert_counts) else 0,
        "n_combination_perturbations": sum(COMBINATION_DELIMITER in p for p in perts),
        "var_names_unique": bool(read_var_names(path).is_unique),
    }


# --------------------------------------------------------------------------
# frozen balanced design
# --------------------------------------------------------------------------


def shared_perturbations(label_sets: Mapping[str, Iterable[str]]) -> list[str]:
    """Sorted perturbations present in every context, excluding the control label.

    Depends only on *identifier presence*. No expression value is consulted.
    """
    sets = [set(v) - {CONTROL_LABEL} for v in label_sets.values()]
    if not sets:
        raise DataIntegrityError("No contexts provided.")
    shared = set.intersection(*sets)
    if not shared:
        raise DataIntegrityError("Contexts share no perturbations.")
    return sorted(shared)


def shared_genes(gene_sets: Mapping[str, Iterable[str]]) -> list[str]:
    """Sorted genes present in every context. Identifier presence only."""
    sets = [set(v) for v in gene_sets.values()]
    if not sets:
        raise DataIntegrityError("No contexts provided.")
    shared = set.intersection(*sets)
    if not shared:
        raise DataIntegrityError("Contexts share no genes.")
    return sorted(shared)


# --------------------------------------------------------------------------
# memory-safe pseudobulk
# --------------------------------------------------------------------------


def _column_selector(var_names: pd.Index, genes: Sequence[str]) -> sparse.csr_matrix:
    """(n_vars x n_selected) 0/1 matrix selecting ``genes`` in the given order."""
    idx = var_names.get_indexer(pd.Index(genes))
    if (idx < 0).any():
        missing = [g for g, i in zip(genes, idx, strict=True) if i < 0]
        raise DataIntegrityError(f"{len(missing)} requested genes absent, e.g. {missing[:5]}")
    n_sel = len(genes)
    return sparse.csr_matrix(
        (np.ones(n_sel, dtype=np.float64), (idx, np.arange(n_sel))),
        shape=(len(var_names), n_sel),
    )


@dataclass(frozen=True)
class ContextPseudobulk:
    """Pseudobulk means for one context on a frozen gene set."""

    context: str
    genes: tuple[str, ...] = field(repr=False)
    perturbations: tuple[str, ...] = field(repr=False)
    control_mean: np.ndarray = field(repr=False)  # (G,)
    perturbation_means: np.ndarray = field(repr=False)  # (P, G)
    cell_counts: np.ndarray = field(repr=False)  # (P,)
    n_control_cells: int = 0
    n_nonzero_seen: int = 0

    @property
    def delta(self) -> np.ndarray:
        """``perturbation_mean - control_mean``, shape (P, G)."""
        return self.perturbation_means - self.control_mean[None, :]


def pseudobulk(
    path: str | Path,
    *,
    genes: Sequence[str],
    perturbations: Sequence[str],
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    context: str | None = None,
) -> ContextPseudobulk:
    """Streaming per-perturbation and control means on a frozen gene set.

    One pass over the file. The full matrix is never densified: each row chunk is
    restricted to the frozen genes by a sparse selector and accumulated with a
    sparse group-indicator product.
    """
    path = Path(path)
    var_names = read_var_names(path)
    selector = _column_selector(var_names, genes)
    labels = cell_labels(path)

    n_p, n_g = len(perturbations), len(genes)
    row_of = {p: i for i, p in enumerate(perturbations)}
    # group id per cell: 0..n_p-1 for a frozen perturbation, n_p for control, -1 to ignore
    group = np.full(len(labels), -1, dtype=np.int64)
    for i, lab in enumerate(labels):
        if lab == CONTROL_LABEL:
            group[i] = n_p
        else:
            j = row_of.get(lab)
            if j is not None:
                group[i] = j

    sums = np.zeros((n_p + 1, n_g), dtype=np.float64)
    counts = np.zeros(n_p + 1, dtype=np.int64)
    nnz_seen = 0

    for start, chunk in stream_row_chunks(path, chunk_size=chunk_size):
        g_chunk = group[start : start + chunk.shape[0]]
        keep = g_chunk >= 0
        if not keep.any():
            continue
        rows = np.flatnonzero(keep)
        sub = chunk[rows] @ selector  # (n_kept, G), still sparse
        nnz_seen += int(sub.nnz)
        ids = g_chunk[rows]
        indicator = sparse.csr_matrix(
            (np.ones(len(ids)), (ids, np.arange(len(ids)))), shape=(n_p + 1, len(ids))
        )
        sums += np.asarray((indicator @ sub).todense())
        counts += np.bincount(ids, minlength=n_p + 1)

    if counts[n_p] == 0:
        raise DataIntegrityError(f"{path.name}: no control cells found.")
    empty = [perturbations[i] for i in np.flatnonzero(counts[:n_p] == 0)]
    if empty:
        raise DataIntegrityError(f"{path.name}: {len(empty)} frozen perturbations have no cells.")

    return ContextPseudobulk(
        context=context or path.stem,
        genes=tuple(genes),
        perturbations=tuple(perturbations),
        control_mean=sums[n_p] / counts[n_p],
        perturbation_means=sums[:n_p] / counts[:n_p, None],
        cell_counts=counts[:n_p],
        n_control_cells=int(counts[n_p]),
        n_nonzero_seen=nnz_seen,
    )


def load_perturbation_cells(
    path: str | Path,
    *,
    genes: Sequence[str],
    perturbations: Sequence[str],
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> tuple[sparse.csr_matrix, np.ndarray]:
    """Cells belonging to the frozen perturbations, on the frozen gene set.

    Returns ``(X, group)`` where ``X`` is a CSR matrix of the retained cells and
    ``group[i]`` is the index into ``perturbations`` for row ``i``. Used by the
    split-half analysis, which needs cell-level access. Loads **one context** —
    callers must not hold more than one at a time.
    """
    path = Path(path)
    selector = _column_selector(read_var_names(path), genes)
    labels = cell_labels(path)
    row_of = {p: i for i, p in enumerate(perturbations)}
    group_all = np.array([row_of.get(lab, -1) for lab in labels], dtype=np.int64)

    blocks, groups = [], []
    for start, chunk in stream_row_chunks(path, chunk_size=chunk_size):
        g_chunk = group_all[start : start + chunk.shape[0]]
        rows = np.flatnonzero(g_chunk >= 0)
        if not len(rows):
            continue
        blocks.append((chunk[rows] @ selector).astype(np.float32))
        groups.append(g_chunk[rows])
    if not blocks:
        raise DataIntegrityError(f"{path.name}: no cells for the frozen perturbations.")
    return sparse.vstack(blocks, format="csr"), np.concatenate(groups)


def library_size_check(
    path: str | Path, *, n_cells: int = 2_000, seed: int = 0, chunk_size: int = DEFAULT_CHUNK_SIZE
) -> dict[str, float]:
    """Check that ``expm1(X)`` row sums look like ``normalize_total(TARGET_SUM)``.

    scPertEval documents ``X`` as ``log1p(CP10K)``. Because genes were filtered
    *after* normalisation, row sums are expected slightly *below* TARGET_SUM
    rather than exactly equal; this reports the distribution so the deviation can
    be stated rather than assumed.
    """
    path = Path(path)
    total_cells = read_shape(path)[0]
    rng = np.random.default_rng(seed)
    wanted = np.sort(rng.choice(total_cells, size=min(n_cells, total_cells), replace=False))
    sums = []
    for start, chunk in stream_row_chunks(path, chunk_size=chunk_size):
        stop = start + chunk.shape[0]
        local = wanted[(wanted >= start) & (wanted < stop)] - start
        if not len(local):
            continue
        block = chunk[local].copy()
        block.data = np.expm1(block.data)
        sums.append(np.asarray(block.sum(axis=1)).ravel())
    values = np.concatenate(sums)
    return {
        "n_cells_sampled": int(values.size),
        "min": float(values.min()),
        "p05": float(np.percentile(values, 5)),
        "median": float(np.median(values)),
        "mean": float(values.mean()),
        "p95": float(np.percentile(values, 95)),
        "max": float(values.max()),
        "target_sum": TARGET_SUM,
        "median_over_target": float(np.median(values) / TARGET_SUM),
    }
