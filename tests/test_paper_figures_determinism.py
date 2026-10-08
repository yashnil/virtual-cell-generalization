"""Deterministic paper-figure rebuilds: stable SVG IDs and sidecars without volatile fields.

Two consecutive ``scripts/paper_figures/run_all.sh`` runs must leave the git working tree
unchanged, so no tracked artefact may embed uuid-based SVG IDs, a wall-clock date or the live
git HEAD/dirty state.
"""

from __future__ import annotations

import io
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "paper_figures"))

import _common  # noqa: E402
import style  # noqa: E402

VOLATILE = {"generated", "git_head"}


def _svg_bytes() -> bytes:
    fig = style.new_figure(2.0, 1.5)
    ax = fig.add_subplot(111)
    ax.plot([0, 1, 2], [0, 1, 0], marker="o")
    ax.bar([0.5, 1.5], [0.3, 0.6], hatch="///")
    buf = io.BytesIO()
    fig.savefig(buf, format="svg", metadata={"Date": None})
    plt.close(fig)
    return buf.getvalue()


def test_style_sets_fixed_svg_hashsalt():
    style.apply()
    assert matplotlib.rcParams["svg.hashsalt"] == style.SVG_HASHSALT


def test_svg_output_is_byte_identical_across_renders():
    first, second = _svg_bytes(), _svg_bytes()
    assert b"clip-path=" in first and b"<pattern" in first
    assert first == second


def test_write_source_omits_volatile_fields(tmp_path, monkeypatch):
    monkeypatch.setattr(_common, "ROOT", tmp_path)
    monkeypatch.setattr(_common, "SRC_DIR", tmp_path / "src")
    monkeypatch.setattr(_common, "_MANIFEST", {})
    table = pd.DataFrame({"x": [1.0, 2.0]})
    kw = dict(sources=[], reports=[], build_script="b.py", plot_script="p.py", notes="n")
    _common.write_source("t", table, **kw)
    side = tmp_path / "src" / "t.provenance.json"
    first = side.read_bytes()
    _common.write_source("t", table, **kw)
    assert side.read_bytes() == first
    assert not VOLATILE & json.loads(first).keys()


@pytest.mark.parametrize(
    "path",
    sorted((ROOT / "data/figure_sources/paper").glob("*.provenance.json"))
    + sorted((ROOT / "reports/paper_figures").glob("*.manifest.json")),
    ids=lambda p: p.name,
)
def test_tracked_paper_sidecars_have_no_volatile_fields(path):
    assert not VOLATILE & json.loads(path.read_text()).keys()
