"""C5 permission gate: X-Atlas can never be read while permission is not APPROVED."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

from virtual_cell.competition_v2 import licensing

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def c5():
    spec = importlib.util.spec_from_file_location(
        "run_c5_xatlas", ROOT / "scripts" / "competition_v2" / "run_c5_xatlas.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_current_status_is_parsed_and_gate_is_closed(c5):
    status = c5.permission_status()
    assert status in {"PENDING", "APPROVED", "DENIED"}
    if status != "APPROVED":
        assert c5.gate() == (False, f"X-Atlas permission is {status}")


def test_gate_requires_both_status_and_license_table(c5):
    assert c5.gate("PENDING", licensing.STATUS)[0] is False
    assert c5.gate("DENIED", licensing.STATUS)[0] is False
    blocked = dict(licensing.STATUS)
    assert c5.gate("APPROVED", blocked)[0] is False  # status alone is not enough
    green = {**blocked, "HCT116": licensing.GREEN, "HEK293T": licensing.GREEN}
    assert c5.gate("APPROVED", green) == (True, "APPROVED")


def test_status_row_parser(c5):
    assert c5.permission_status("| **status** | **APPROVED** |") == "APPROVED"
    with pytest.raises(ValueError):
        c5.permission_status("no status here")


def test_main_refuses_while_pending(c5, capsys):
    if c5.permission_status() == "APPROVED":
        pytest.skip("permission recorded; the gate is legitimately open")
    import sys

    argv, sys.argv = sys.argv, ["run_c5_xatlas.py"]
    try:
        assert c5.main() == 0
    finally:
        sys.argv = argv
    assert "C5 gated, not run" in capsys.readouterr().out
    assert not (ROOT / "outputs" / "competition_v2" / "c5_xatlas" / "c5_decision.json").exists()
