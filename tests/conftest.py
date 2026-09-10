from pathlib import Path

import pytest

from virtual_cell.data.synthetic import make_synthetic_contexts, write_synthetic_contexts

REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="session")
def synthetic_contexts():
    """In-memory synthetic contexts A/B/C (small, fixed seed)."""
    return make_synthetic_contexts(("A", "B", "C"), n_cells=50, n_genes=30, seed=0)


@pytest.fixture(scope="session")
def synthetic_paths(tmp_path_factory) -> dict[str, Path]:
    """Synthetic contexts written to a temporary directory."""
    out = tmp_path_factory.mktemp("synthetic")
    return write_synthetic_contexts(out, ("A", "B", "C"), n_cells=50, n_genes=30, seed=0)
