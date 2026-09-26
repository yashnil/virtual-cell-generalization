"""Write every figure-source table in ``data/figure_sources/`` from frozen artifacts.

Usage: ``uv run python scripts/figures/extract_figure_sources.py [name ...]``
(no names = all). Deterministic: rerunning reproduces every CSV byte for byte;
only the ``generated`` date and ``git_head`` in the provenance sidecars change.
"""

from __future__ import annotations

import sys

from virtual_cell.visualization import common
from virtual_cell.visualization.sources import EXTRACTORS


def main(names: list[str]) -> None:
    for name in names or list(EXTRACTORS):
        ex = EXTRACTORS[name]()
        if not common.numeric_is_finite(ex.table, ex.allow_nan):
            raise SystemExit(f"{name}: non-finite values outside the allowed columns")
        path = common.write_source(
            ex.name,
            ex.table,
            sources=ex.sources,
            report=ex.report,
            figure_script=ex.figure_script,
            notes=ex.notes,
        )
        print(f"  {path.relative_to(common.ROOT)}  ({len(ex.table)} rows)")


if __name__ == "__main__":
    main(sys.argv[1:])
