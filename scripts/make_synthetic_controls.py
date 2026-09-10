"""Generate synthetic control-cell contexts A/B/C as ``.h5ad`` files.

Usage (from the repository root)::

    uv run python scripts/make_synthetic_controls.py
    uv run python scripts/make_synthetic_controls.py --n-cells 500 --seed 7

The data are Poisson noise with no biological content and exist only so the
pipeline can be developed and tested before real data are available.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from virtual_cell.data.synthetic import DEFAULT_CONTEXTS, write_synthetic_contexts

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = REPO_ROOT / "data" / "raw" / "synthetic"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--contexts", nargs="+", default=list(DEFAULT_CONTEXTS))
    parser.add_argument("--n-cells", type=int, default=200)
    parser.add_argument("--n-genes", type=int, default=100)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    paths = write_synthetic_contexts(
        args.out_dir,
        args.contexts,
        n_cells=args.n_cells,
        n_genes=args.n_genes,
        seed=args.seed,
    )
    for label, path in paths.items():
        shown = path.relative_to(REPO_ROOT) if path.is_relative_to(REPO_ROOT) else path
        print(f"{label}: wrote {shown}")


if __name__ == "__main__":
    main()
