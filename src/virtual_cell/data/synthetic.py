"""Synthetic control-cell contexts for testing the pipeline without real data.

These datasets are deliberately simple (independent Poisson counts) and carry
no biology. They exist so that loaders, summaries, splits, and tests can be
developed and exercised before real Perturb-seq data or Arc credentials are
available. Nothing derived from them should be interpreted scientifically.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

DEFAULT_CONTEXTS: tuple[str, ...] = ("A", "B", "C")


def gene_names(n_genes: int) -> list[str]:
    return [f"GENE_{i:03d}" for i in range(n_genes)]


def make_synthetic_context(
    context: str,
    *,
    n_cells: int,
    n_genes: int,
    rng: np.random.Generator,
    context_effect: float = 1.0,
    base_rate: np.ndarray | None = None,
    context_key: str = "context",
) -> ad.AnnData:
    """Build one synthetic context of Poisson-distributed control-cell counts.

    Parameters
    ----------
    context:
        Context label stored in ``obs[context_key]`` (preserved exactly).
    n_cells, n_genes:
        Matrix shape.
    rng:
        NumPy random generator; controls all randomness.
    context_effect:
        Multiplicative scaling of the per-gene rates for this context.
    base_rate:
        Per-gene Poisson rates. If ``None`` a fresh vector is drawn from
        ``Uniform(0.1, 3.0)``.
    """
    if base_rate is None:
        base_rate = rng.uniform(0.1, 3.0, size=n_genes)
    base_rate = np.asarray(base_rate, dtype=float)
    if base_rate.shape != (n_genes,):
        raise ValueError(f"base_rate must have shape ({n_genes},), got {base_rate.shape}")

    counts = rng.poisson(lam=base_rate * context_effect, size=(n_cells, n_genes))

    obs = pd.DataFrame(
        {context_key: context},
        index=[f"{context}_cell_{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=pd.Index(gene_names(n_genes), name="gene"))

    adata = ad.AnnData(X=sparse.csr_matrix(counts), obs=obs, var=var)
    adata.uns["synthetic"] = {
        "context_effect": float(context_effect),
        "generator": "poisson",
        "note": "Synthetic data with no biological content; for pipeline testing only.",
    }
    return adata


def make_synthetic_contexts(
    contexts: Sequence[str] = DEFAULT_CONTEXTS,
    *,
    n_cells: int = 200,
    n_genes: int = 100,
    seed: int = 42,
    context_effect_step: float = 0.20,
    shared_base_rate: bool = True,
    context_modulation_sd: float = 0.3,
) -> dict[str, ad.AnnData]:
    """Build several synthetic contexts with mildly different baselines.

    Each context is scaled by ``1 + index * context_effect_step`` so that
    library sizes differ modestly between contexts.

    With ``shared_base_rate=True`` (default) all contexts share one per-gene
    rate vector, modulated per context by log-normal noise with standard
    deviation ``context_modulation_sd`` on the log scale. This mimics the fact
    that real cell types have strongly correlated basal profiles with a subset
    of context-specific genes. With ``shared_base_rate=False`` every context
    draws independent rates (the original behaviour), which yields essentially
    uncorrelated basal profiles.
    """
    rng = np.random.default_rng(seed)
    shared = rng.uniform(0.1, 3.0, size=n_genes) if shared_base_rate else None
    out: dict[str, ad.AnnData] = {}
    for idx, context in enumerate(contexts):
        base_rate = None
        if shared is not None:
            base_rate = shared * np.exp(rng.normal(0.0, context_modulation_sd, size=n_genes))
        out[context] = make_synthetic_context(
            context,
            n_cells=n_cells,
            n_genes=n_genes,
            rng=rng,
            context_effect=1.0 + idx * context_effect_step,
            base_rate=base_rate,
        )
    return out


def write_synthetic_contexts(
    out_dir: str | Path,
    contexts: Sequence[str] = DEFAULT_CONTEXTS,
    *,
    n_cells: int = 200,
    n_genes: int = 100,
    seed: int = 42,
) -> dict[str, Path]:
    """Generate synthetic contexts and write ``context_<label>.h5ad`` files."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    adatas = make_synthetic_contexts(contexts, n_cells=n_cells, n_genes=n_genes, seed=seed)
    paths: dict[str, Path] = {}
    for label, adata in adatas.items():
        path = out_dir / f"context_{label}.h5ad"
        adata.write_h5ad(path)
        paths[label] = path
    return paths
