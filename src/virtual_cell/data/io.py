"""Loading single-cell AnnData contexts with integrity checks.

A "context" in this project is one cellular background (a cell line, cell type,
or experimental system) stored as a single ``.h5ad`` file. The loader enforces the
data assumptions the rest of the pipeline relies on so that violations fail
loudly at load time rather than silently corrupting downstream statistics.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse

DEFAULT_CONTEXT_KEY = "context"


class DataIntegrityError(ValueError):
    """Raised when an AnnData object violates a project data assumption."""


def _matrix_values(X) -> np.ndarray:
    """Return the stored (non-zero for sparse) values of ``X`` as a flat array."""
    if sparse.issparse(X):
        return np.asarray(X.data).ravel()
    return np.asarray(X).ravel()


def check_counts(adata: ad.AnnData, *, require_integer: bool = True) -> None:
    """Check that ``adata.X`` looks like a raw UMI count matrix.

    Checks that the matrix exists, is numeric, contains no NaN/inf values, has no
    negative entries, and (optionally) is integer-valued.
    """
    if adata.X is None:
        raise DataIntegrityError("adata.X is missing; expected a count matrix.")

    values = _matrix_values(adata.X)
    if not np.issubdtype(values.dtype, np.number):
        raise DataIntegrityError(f"adata.X has non-numeric dtype {values.dtype}.")
    if values.size == 0:
        return
    if not np.all(np.isfinite(values)):
        raise DataIntegrityError("adata.X contains NaN or infinite values.")
    if values.min() < 0:
        raise DataIntegrityError("adata.X contains negative values; expected counts.")
    if require_integer and not np.issubdtype(values.dtype, np.integer):
        if not np.allclose(values, np.round(values)):
            raise DataIntegrityError(
                "adata.X contains non-integer values; expected raw counts. "
                "Pass require_integer=False to load normalized data."
            )


def check_context_labels(
    adata: ad.AnnData,
    *,
    context_key: str = DEFAULT_CONTEXT_KEY,
    expected_context: str | None = None,
) -> str:
    """Check the context label column and return the single context label.

    Every file is expected to describe exactly one context. The label is compared
    with ``expected_context`` (if given) using exact string equality so that
    labels are never silently normalised or remapped.
    """
    if context_key not in adata.obs.columns:
        raise DataIntegrityError(
            f"obs is missing the context column {context_key!r}; "
            f"available columns: {list(adata.obs.columns)}"
        )
    labels = adata.obs[context_key]
    if labels.isna().any():
        raise DataIntegrityError(f"obs[{context_key!r}] contains missing values.")
    unique = labels.astype(str).unique().tolist()
    if len(unique) != 1:
        raise DataIntegrityError(
            f"Expected exactly one context label in obs[{context_key!r}], found {unique}."
        )
    label = unique[0]
    if expected_context is not None and label != expected_context:
        raise DataIntegrityError(
            f"Context label mismatch: file says {label!r}, expected {expected_context!r}."
        )
    return label


def check_integrity(
    adata: ad.AnnData,
    *,
    context_key: str | None = DEFAULT_CONTEXT_KEY,
    expected_context: str | None = None,
    require_integer: bool = True,
) -> None:
    """Run all integrity checks on an AnnData object.

    Parameters
    ----------
    adata:
        Dataset to check.
    context_key:
        Name of the ``obs`` column holding the context label. Pass ``None`` to
        skip the context-label checks.
    expected_context:
        If given, the file's context label must equal this string exactly.
    require_integer:
        If True, ``adata.X`` must be integer-valued (raw counts).
    """
    if adata.n_obs == 0:
        raise DataIntegrityError("Dataset contains no cells.")
    if adata.n_vars == 0:
        raise DataIntegrityError("Dataset contains no genes.")
    if adata.var_names.has_duplicates:
        dupes = adata.var_names[adata.var_names.duplicated()].unique().tolist()
        raise DataIntegrityError(f"Duplicate gene names detected: {dupes[:10]}")
    if adata.obs_names.has_duplicates:
        raise DataIntegrityError("Duplicate cell identifiers detected in obs_names.")
    check_counts(adata, require_integer=require_integer)
    if context_key is not None:
        check_context_labels(
            adata, context_key=context_key, expected_context=expected_context
        )


def load_context(
    path: str | Path,
    *,
    context_key: str | None = DEFAULT_CONTEXT_KEY,
    expected_context: str | None = None,
    require_integer: bool = True,
) -> ad.AnnData:
    """Load one single-cell context from ``.h5ad`` and run integrity checks.

    Parameters
    ----------
    path:
        Path to an ``.h5ad`` file describing a single cellular context.
    context_key:
        ``obs`` column holding the context label (``None`` skips label checks).
    expected_context:
        Optional exact context label the file must carry.
    require_integer:
        Require ``adata.X`` to hold integer-valued raw counts.

    Returns
    -------
    anndata.AnnData
        The loaded dataset, unmodified apart from what ``anndata`` does on read.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    DataIntegrityError
        If any integrity check fails.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Dataset not found: {path}")
    if path.suffix != ".h5ad":
        raise ValueError(f"Expected an .h5ad file, got: {path.name}")

    adata = ad.read_h5ad(path)
    check_integrity(
        adata,
        context_key=context_key,
        expected_context=expected_context,
        require_integer=require_integer,
    )
    return adata


def load_contexts(
    paths: Mapping[str, str | Path] | Iterable[str | Path],
    *,
    context_key: str = DEFAULT_CONTEXT_KEY,
    require_integer: bool = True,
) -> dict[str, ad.AnnData]:
    """Load several contexts and key them by their context label.

    ``paths`` may be a mapping ``{expected_label: path}`` (labels are then
    verified against the files) or a plain iterable of paths (labels are read
    from the files). Duplicate labels across files raise an error.
    """
    if isinstance(paths, Mapping):
        items = [(label, Path(p)) for label, p in paths.items()]
    else:
        items = [(None, Path(p)) for p in paths]

    contexts: dict[str, ad.AnnData] = {}
    for expected, path in items:
        adata = load_context(
            path,
            context_key=context_key,
            expected_context=expected,
            require_integer=require_integer,
        )
        label = str(adata.obs[context_key].iloc[0])
        if label in contexts:
            raise DataIntegrityError(
                f"Context label {label!r} appears in more than one file."
            )
        contexts[label] = adata
    return contexts
