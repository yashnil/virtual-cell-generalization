r"""Assembling a submission: response -> multiplicative effect -> counts -> ``.h5ad``.

:mod:`virtual_cell.arc.panel` maps a response onto the 18,533-gene axis and
:mod:`virtual_cell.arc.generate` turns a profile into cells. Between and after
them sit three jobs this module does, and nothing else does:

1. **Translate.** The model predicts a pseudobulk ``delta`` in the space the
   public data is standardised in — the per-cell mean of ``log1p(CP10K)``. The
   generators want a *multiplicative per-gene effect* on composition. The two
   are related but not equal, and the conversion is stated here explicitly
   rather than improvised at each call site.

2. **Decide the unsupported genes.** Most of the panel carries no predicted
   response for most perturbations. The frozen rule is that such a gene keeps
   its **target-context control distribution** — which is the ``lfc = 0`` case
   of the transport generator, and is *not* the same as writing zero counts.
   :func:`panel_log2_fold_change` is where that rule lives.

3. **Write 360,000 cells without holding them.** A full submission is roughly
   two billion nonzeros; materialising it as one sparse matrix costs more RAM
   than packaging it does. :class:`SubmissionWriter` appends CSR blocks
   straight into the HDF5 file as they are generated.

The fold-change convention
--------------------------
``delta`` is a difference of mean-``log1p(CP10K)`` profiles, so
``exp(delta) - 1`` is not a fold change and ``2**delta`` is not either. What is
well defined is the ratio of the two *profiles*:

.. math::
    \text{lfc}_g = \log_2 \frac{\operatorname{expm1}(b_g + \delta_g) + f}
                              {\operatorname{expm1}(b_g) + f}

with ``b`` the control profile in the same space and ``f`` a floor below one
UMI's worth of signal, which keeps a gene the controls never detected from
producing an unbounded multiplier. Because ``delta`` is a mean of per-cell log
ratios, this is a geometric-mean fold change — the right object to push through
a multiplicative transport, and *not* the same as the ratio of arithmetic
means. Whether the generated cells actually realise the intended mean response
is therefore an empirical question, measured rather than assumed; see
``reports/arc_count_space_baseline_v1.md``.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from virtual_cell.arc.panel import PanelMap, Support

__all__ = [
    "TARGET_SUM",
    "FLOOR_FRACTION",
    "LFC_CLIP",
    "log2_fold_change",
    "predicted_profile",
    "panel_log2_fold_change",
    "SubmissionWriter",
    "BundleReport",
    "inspect_bundle",
]

#: The normalisation the public responses are defined in (``log1p(CP10K)``).
TARGET_SUM = 1.0e4

#: The floor ``f``, as a fraction of ``TARGET_SUM``. ``1e-6`` is one hundredth
#: of a single UMI in a 10,000-UMI cell: small enough not to damp a real
#: effect, large enough that a gene at exactly zero cannot divide by nothing.
FLOOR_FRACTION = 1.0e-6

#: Fold changes are clipped here before reaching a generator. A 2**10 swing is
#: far outside anything CRISPRi pseudobulk produces, so the clip is a guard
#: against a degenerate profile rather than a modelling choice.
LFC_CLIP = 10.0


def log2_fold_change(
    control_profile: np.ndarray,
    delta: np.ndarray,
    *,
    target_sum: float = TARGET_SUM,
    floor_fraction: float = FLOOR_FRACTION,
    clip: float = LFC_CLIP,
) -> np.ndarray:
    """Per-gene multiplicative effect implied by a pseudobulk ``delta``.

    ``control_profile`` is the target context's control profile in
    mean-``log1p`` space, ``delta`` the predicted response in the same space
    (``(G,)`` or ``(P, G)``). Returns an array of the same shape as ``delta``.
    """
    b = np.asarray(control_profile, dtype=np.float64)
    d = np.asarray(delta, dtype=np.float64)
    if d.shape[-1] != b.shape[-1]:
        raise ValueError(f"delta has {d.shape[-1]} genes, control profile has {b.shape[-1]}")
    floor = floor_fraction * target_sum
    basal = np.expm1(b)
    perturbed = np.expm1(b + d)
    lfc = np.log2((np.clip(perturbed, 0.0, None) + floor) / (np.clip(basal, 0.0, None) + floor))
    return np.clip(lfc, -clip, clip)


def predicted_profile(
    control_profile: np.ndarray,
    delta: np.ndarray,
    *,
    target_sum: float = TARGET_SUM,
) -> np.ndarray:
    """Predicted composition on the ``target_sum`` scale, renormalised to it.

    This is what ``G2`` draws from. Renormalisation matters: a predicted
    response that raised every gene would otherwise smuggle in a library-size
    change the model never claimed.
    """
    b = np.asarray(control_profile, dtype=np.float64)
    d = np.asarray(delta, dtype=np.float64)
    profile = np.clip(np.expm1(b + d), 0.0, None)
    total = profile.sum(axis=-1, keepdims=True)
    return np.divide(profile * target_sum, total, out=np.zeros_like(profile), where=total > 0)


def panel_log2_fold_change(
    panel_map: PanelMap,
    control_profile_panel: np.ndarray,
    delta_source: np.ndarray,
    *,
    target_sum: float = TARGET_SUM,
    floor_fraction: float = FLOOR_FRACTION,
    clip: float = LFC_CLIP,
) -> np.ndarray:
    """Lift a source-space response to the panel and state the unsupported rule.

    Panel genes outside :attr:`Support.PREDICTED` receive ``lfc = 0``. Under
    the transport generator that means the gene is emitted with the target
    context's own control distribution — its real counts, its real detection
    rate, its real dispersion. It is **not** zero-filled, and it is not
    predicted to be unchanged in any stronger sense than "we have no evidence
    that it changes".
    """
    d = np.atleast_2d(np.asarray(delta_source, dtype=np.float64))
    if d.shape[-1] != len(panel_map.source_genes):
        raise ValueError(
            f"delta has {d.shape[-1]} genes, source space has {len(panel_map.source_genes)}"
        )
    basal = np.asarray(control_profile_panel, dtype=np.float64)
    if basal.shape != (panel_map.n_panel,):
        raise ValueError(f"control profile must be ({panel_map.n_panel},)")

    take = panel_map.mask(Support.PREDICTED)
    cols = panel_map.source_index[take]
    out = np.zeros((d.shape[0], panel_map.n_panel), dtype=np.float64)
    out[:, take] = log2_fold_change(
        basal[take],
        d[:, cols],
        target_sum=target_sum,
        floor_fraction=floor_fraction,
        clip=clip,
    )
    return out


# --------------------------------------------------------------------------
# writing
# --------------------------------------------------------------------------


class SubmissionWriter:
    """Append CSR blocks of generated cells into an ``.h5ad``, one block at a time.

    Usage::

        with SubmissionWriter(path, genes) as w:
            for context, pert, block in ...:
                w.append(block, target_gene=pert, context=context)

    ``indptr`` is written as ``int64`` unconditionally: it holds *cumulative*
    nonzero counts, a full submission sits close enough to ``2**31`` that a run
    fitting in ``int32`` one day would overflow the next, and an overflowed
    offset array is silently wrong rather than loudly broken. ``indices`` holds
    *column* numbers, bounded by 18,533, so it stays ``int32`` — writing it as
    ``int64`` would add four bytes to every one of two billion nonzeros and
    grow the file by eight gigabytes for no information.
    """

    def __init__(
        self,
        path: str | Path,
        genes: Sequence[str],
        *,
        pert_col: str = "target_gene",
        context_col: str = "context",
        dtype: str = "float32",
        chunk: int = 1 << 20,
    ) -> None:
        import h5py

        self.path = Path(path)
        self.genes = pd.Index(list(genes), dtype=object)
        if self.genes.has_duplicates:
            raise ValueError("gene names must be unique")
        self.pert_col = pert_col
        self.context_col = context_col
        self._dtype = dtype
        self._n_genes = len(self.genes)
        self._n_cells = 0
        self._nnz = 0
        self._perts: list[str] = []
        self._contexts: list[str] = []

        self._file = h5py.File(self.path, "w")
        grp = self._file.create_group("X")
        grp.attrs["encoding-type"] = "csr_matrix"
        grp.attrs["encoding-version"] = "0.1.0"
        grp.create_dataset("data", shape=(0,), maxshape=(None,), dtype=dtype, chunks=(chunk,))
        grp.create_dataset("indices", shape=(0,), maxshape=(None,), dtype="int32", chunks=(chunk,))
        indptr = grp.create_dataset(
            "indptr", shape=(1,), maxshape=(None,), dtype="int64", chunks=(min(chunk, 1 << 16),)
        )
        indptr[0] = 0

    # -- context manager -------------------------------------------------
    def __enter__(self) -> SubmissionWriter:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    # -- writing ---------------------------------------------------------
    def append(self, block: np.ndarray, *, target_gene: str, context: str) -> None:
        """Append one perturbation's cells, all carrying the same labels."""
        dense = np.asarray(block)
        if dense.ndim != 2 or dense.shape[1] != self._n_genes:
            raise ValueError(f"block must be (n_cells, {self._n_genes}), got {dense.shape}")
        if np.any(dense < 0):
            raise ValueError("counts must be non-negative")
        if not np.all(np.isfinite(dense)):
            raise ValueError("counts must be finite")
        csr = sparse.csr_matrix(dense)

        grp = self._file["X"]
        data, indices, indptr = grp["data"], grp["indices"], grp["indptr"]
        n_new = int(csr.nnz)
        data.resize((self._nnz + n_new,))
        indices.resize((self._nnz + n_new,))
        if n_new:
            data[self._nnz :] = csr.data.astype(self._dtype)
            indices[self._nnz :] = csr.indices.astype(np.int32)
        rows = csr.shape[0]
        indptr.resize((self._n_cells + rows + 1,))
        indptr[self._n_cells + 1 :] = csr.indptr[1:].astype(np.int64) + self._nnz

        self._nnz += n_new
        self._n_cells += rows
        self._perts.extend([target_gene] * rows)
        self._contexts.extend([context] * rows)

    def close(self) -> None:
        """Finalise ``obs``, ``var`` and the root attributes, then close."""
        from anndata.io import write_elem

        f = self._file
        f["X"].attrs["shape"] = np.array([self._n_cells, self._n_genes], dtype="int64")
        obs = pd.DataFrame(
            {
                self.pert_col: pd.Categorical(self._perts),
                self.context_col: pd.Categorical(self._contexts),
            },
            index=pd.Index([f"cell_{i}" for i in range(self._n_cells)], dtype=object),
        )
        var = pd.DataFrame(index=self.genes)
        write_elem(f, "obs", obs)
        write_elem(f, "var", var)
        for name in ("layers", "obsm", "varm", "obsp", "varp", "uns"):
            write_elem(f, name, {})
        f.attrs["encoding-type"] = "anndata"
        f.attrs["encoding-version"] = "0.1.0"
        f.close()

    @property
    def n_cells(self) -> int:
        return self._n_cells

    @property
    def nnz(self) -> int:
        return self._nnz


# --------------------------------------------------------------------------
# reading back
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class BundleReport:
    """Everything a local pre-submission check needs, read without densifying."""

    n_cells: int
    n_genes: int
    nnz: int
    dtype: str
    gene_order_matches: bool
    integer_valued: bool
    non_negative: bool
    finite: bool
    max_counts_per_cell: int
    median_counts_per_cell: float
    min_counts_per_cell: int
    median_genes_detected: float
    density: float
    cells_per_group: pd.Series
    contexts: tuple[str, ...]
    n_perturbations_per_context: dict[str, int]

    def as_dict(self) -> dict[str, object]:
        return {
            "n_cells": self.n_cells,
            "n_genes": self.n_genes,
            "nnz": self.nnz,
            "dtype": self.dtype,
            "gene_order_matches": self.gene_order_matches,
            "integer_valued": self.integer_valued,
            "non_negative": self.non_negative,
            "finite": self.finite,
            "max_counts_per_cell": self.max_counts_per_cell,
            "median_counts_per_cell": self.median_counts_per_cell,
            "min_counts_per_cell": self.min_counts_per_cell,
            "median_genes_detected": self.median_genes_detected,
            "density": self.density,
            "contexts": list(self.contexts),
            "n_perturbations_per_context": self.n_perturbations_per_context,
        }


def inspect_bundle(
    path: str | Path,
    expected_genes: Iterable[str],
    *,
    pert_col: str = "target_gene",
    context_col: str = "context",
    chunk_rows: int = 5_000,
) -> BundleReport:
    """Stream a written bundle and report the properties ``vcc prep`` will check.

    ``chunk_rows`` is deliberately small: a full submission carries roughly
    5,700 nonzeros per cell, so 20,000 rows would be a 1.4 GB block and the
    inspection would cost more memory than the writing did.
    """
    import h5py

    from virtual_cell.data.arc2026 import stream_row_chunks

    path = Path(path)
    expected = pd.Index(list(expected_genes), dtype=object)
    with h5py.File(path, "r") as f:
        shape = tuple(int(v) for v in f["X"].attrs["shape"])
        nnz = int(f["X/data"].shape[0])
        dtype = str(f["X/data"].dtype)
        var = f["var"]
        node = var[var.attrs.get("_index", "_index")]
        values = node["values"] if isinstance(node, h5py.Group) else node
        genes = pd.Index([v.decode() if isinstance(v, bytes) else str(v) for v in values[:]])
        obs = f["obs"]

        def _col(name: str) -> np.ndarray:
            n = obs[name]
            if isinstance(n, h5py.Group) and "categories" in n:
                cats = np.array(
                    [c.decode() if isinstance(c, bytes) else str(c) for c in n["categories"][:]],
                    dtype=object,
                )
                return cats[n["codes"][:]]
            raw = n["values"] if isinstance(n, h5py.Group) else n
            return np.array(
                [v.decode() if isinstance(v, bytes) else str(v) for v in raw[:]], dtype=object
            )

        perts = _col(pert_col)
        contexts = _col(context_col)

    totals: list[np.ndarray] = []
    detected: list[np.ndarray] = []
    integer_valued = True
    non_negative = True
    finite = True
    for _start, chunk in stream_row_chunks(path, chunk_size=chunk_rows):
        data = chunk.data
        if data.size:
            finite = finite and bool(np.all(np.isfinite(data)))
            non_negative = non_negative and bool(np.all(data >= 0))
            integer_valued = integer_valued and bool(np.all(data == np.floor(data)))
        totals.append(np.asarray(chunk.sum(axis=1)).ravel())
        detected.append(np.diff(chunk.indptr))
    lib = np.concatenate(totals) if totals else np.zeros(0)
    det = np.concatenate(detected) if detected else np.zeros(0)

    frame = pd.DataFrame({"context": contexts, "pert": perts})
    per_context = {
        str(c): int(frame.loc[frame["context"] == c, "pert"].nunique())
        for c in sorted(set(contexts.tolist()))
    }
    return BundleReport(
        n_cells=shape[0],
        n_genes=shape[1],
        nnz=nnz,
        dtype=dtype,
        gene_order_matches=bool(genes.equals(expected)),
        integer_valued=integer_valued,
        non_negative=non_negative,
        finite=finite,
        max_counts_per_cell=int(lib.max()) if lib.size else 0,
        median_counts_per_cell=float(np.median(lib)) if lib.size else 0.0,
        # NOT ``lib.min(initial=0)``: ``initial`` seeds the reduction, so that
        # form reports 0 for every non-empty bundle and would hide an empty cell
        # rather than reveal one.
        min_counts_per_cell=int(lib.min()) if lib.size else 0,
        median_genes_detected=float(np.median(det)) if det.size else 0.0,
        density=float(nnz / (shape[0] * shape[1])) if shape[0] and shape[1] else 0.0,
        cells_per_group=frame.groupby(["context", "pert"], observed=True).size(),
        contexts=tuple(sorted(set(contexts.tolist()))),
        n_perturbations_per_context=per_context,
    )
