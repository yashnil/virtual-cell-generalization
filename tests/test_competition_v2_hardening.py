"""Release-day hardening: the vcc CLI compatibility gate and robust official-bundle discovery.

Synthetic, tiny bundles: no data dependency.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from virtual_cell.competition_v2 import panel, vcc_compat

GENES = ["G1", "G2", "G3", "G4"]
CTX = ["D", "E", "F"]
PREP_FLAGS = (
    "-g, --genes",
    "--perts",
    "--contexts",
    "--cells-per-pert",
    "--expected-gene-dim",
    "--verify-targets / --no-verify-targets",
    "--check-cell-counts / --no-check-cell-counts",
    "--dry-run",
    "--json",
)
HELP_OK = "\n".join(f"  {o} TEXT   something" for o in PREP_FLAGS)


# ----------------------------------------------------------------------------- CLI gate


def _info(**kw) -> vcc_compat.CliInfo:
    base = dict(
        executable="/x/vcc",
        interpreter="/x/python",
        version="0.2.0",
        package_version="0.2.0",
        prep_help=HELP_OK,
        run_prep_params=list(vcc_compat.REQUIRED_RUN_PREP_KWARGS) + ["force"],
        has_read_h5ad=True,
    )
    base.update(kw)
    return vcc_compat.CliInfo(**base)


def test_compatible_cli_metadata_passes_and_is_reported():
    r = vcc_compat.require(_info())
    assert r.compatible and r.tested_version and r.version == "0.2.0"
    assert not r.problems and not r.warnings


def test_untested_but_capable_cli_passes_with_explicit_warning():
    r = vcc_compat.require(_info(version="0.3.1", package_version="0.3.1"))
    assert r.compatible and not r.tested_version
    assert r.warnings and "untested" in r.warnings[0]


@pytest.mark.parametrize(
    "kw, needle",
    [
        ({"prep_help": HELP_OK.replace("--contexts", "--context-set")}, "--contexts"),
        (
            {"prep_help": HELP_OK.replace("--expected-gene-dim", "--gene-dim")},
            "--expected-gene-dim",
        ),
        ({"run_prep_params": ["input_path", "genes_path"]}, "required_contexts"),
        ({"has_read_h5ad": False}, "read_h5ad"),
        ({"version": None}, "version"),
        ({"package_version": "0.1.9"}, "two installations"),
        ({"probe_errors": ["'vcc' is not on PATH"]}, "not on PATH"),
    ],
)
def test_incompatible_cli_fails_cleanly(kw, needle):
    with pytest.raises(vcc_compat.VccCompatibilityError) as exc:
        vcc_compat.require(_info(**kw))
    msg = str(exc.value)
    assert needle in msg and "What to do" in msg


def test_installed_cli_probe_is_readable():
    if shutil.which("vcc") is None:
        pytest.skip("vcc CLI not installed")
    report = vcc_compat.check(vcc_compat.probe())
    assert report.version is not None
    assert isinstance(report.compatible, bool)


# ----------------------------------------------------------------------------- bundle


def _controls(path: Path, label: str, n: int = 6, genes=GENES):
    a = ad.AnnData(
        sparse.csr_matrix(np.ones((n, len(genes)), dtype=np.float32)),
        obs=pd.DataFrame(
            {
                "target_gene": pd.Categorical(["non-targeting"] * n),
                "context": pd.Categorical([label] * n),
            },
            index=pd.Index([f"{label}{i}" for i in range(n)], dtype=object),
        ),
        var=pd.DataFrame(index=pd.Index(genes, dtype=object)),
    )
    a.write_h5ad(path)


def _bundle(d: Path, *, names=None, manifest=None) -> Path:
    d.mkdir(parents=True, exist_ok=True)
    names = names or {c: f"context_{c}.h5ad" for c in CTX}
    for c in CTX:
        _controls(d / names[c], c)
    pd.DataFrame({"gene_name": GENES}).to_csv(d / "gene_names.csv", index=False)
    pd.DataFrame({"target_gene": ["T2", "T1", "T3"]}).to_csv(d / "pert_counts.csv", index=False)
    m = {"contexts": CTX, "cells_per_pert": 4, "n_genes": 4, "n_constructs": 3}
    (d / "manifest.json").write_text(json.dumps(manifest if manifest is not None else m))
    return d


def test_validation_era_filenames_still_load(tmp_path):
    p = panel.load_panel(_bundle(tmp_path / "b"), checksums=False)
    assert p.contexts == tuple(CTX)
    assert set(p.discovery.values()) == {"validation-era filename"}
    assert list(p.targets) == ["T2", "T1", "T3"]


def test_manifest_named_alternate_filenames(tmp_path):
    names = {c: f"final_{c.lower()}_controls.h5ad" for c in CTX}
    m = {"contexts": [{"label": c, "file": names[c]} for c in CTX], "cells_per_pert": 4}
    p = panel.load_panel(_bundle(tmp_path / "b", names=names, manifest=m), checksums=False)
    assert {c: f.name for c, f in p.control_files.items()} == names
    assert set(p.discovery.values()) == {"named by manifest"}
    m2 = {"contexts": CTX, "cells_per_pert": 4, "files": names}
    p2 = panel.load_panel(_bundle(tmp_path / "c", names=names, manifest=m2), checksums=False)
    assert {c: f.name for c, f in p2.control_files.items()} == names


def test_alternate_filenames_identified_by_unique_obs_context(tmp_path):
    names = {"D": "ctrl_1.h5ad", "E": "ctrl_2.h5ad", "F": "ctrl_3.h5ad"}
    p = panel.load_panel(_bundle(tmp_path / "b", names=names), checksums=False)
    assert {c: f.name for c, f in p.control_files.items()} == names
    assert all(v.startswith("unique file with obs") for v in p.discovery.values())


def test_ambiguous_files_are_not_guessed(tmp_path):
    names = {"D": "ctrl_1.h5ad", "E": "ctrl_2.h5ad", "F": "ctrl_3.h5ad"}
    d = _bundle(tmp_path / "b", names=names)
    _controls(d / "ctrl_extra.h5ad", "E")  # a second file claiming context E
    with pytest.raises(panel.PanelError, match="not guessing"):
        panel.load_panel(d, checksums=False)


@pytest.mark.parametrize(
    "manifest, needle",
    [
        ({"cells_per_pert": 4}, "'contexts' is missing"),
        ({"contexts": [], "cells_per_pert": 4}, "'contexts' is missing"),
        ({"contexts": CTX, "cells_per_pert": "400"}, "positive integer"),
        ({"contexts": CTX, "cells_per_pert": 0}, "positive integer"),
        ({"contexts": ["D", "D", "F"], "cells_per_pert": 4}, "duplicate context"),
        ({"contexts": [{"file": "x.h5ad"}], "cells_per_pert": 4}, "no single label"),
        ({"contexts": [7], "cells_per_pert": 4}, "neither a label"),
        ({"contexts": CTX, "cells_per_pert": 4, "n_genes": "many"}, "'n_genes' must be an integer"),
        ({"contexts": CTX, "cells_per_pert": 4, "pert_col": ""}, "'pert_col'"),
    ],
)
def test_malformed_manifest_fails_with_a_diagnostic(tmp_path, manifest, needle):
    d = _bundle(tmp_path / "b", manifest=manifest)
    with pytest.raises(panel.PanelError) as exc:
        panel.load_panel(d, checksums=False)
    msg = str(exc.value)
    assert needle in msg and "schema not understood" in msg and "found keys" in msg


def test_manifest_file_outside_bundle_or_missing_is_rejected(tmp_path):
    d = _bundle(tmp_path / "b")
    bad = {"contexts": [{"label": "D", "file": "../escape.h5ad"}, "E", "F"], "cells_per_pert": 4}
    (d / "manifest.json").write_text(json.dumps(bad))
    with pytest.raises(panel.PanelError, match="outside the bundle"):
        panel.load_panel(d, checksums=False)
    bad["contexts"][0]["file"] = "nope.h5ad"
    (d / "manifest.json").write_text(json.dumps(bad))
    with pytest.raises(panel.PanelError, match="missing or not an .h5ad"):
        panel.load_panel(d, checksums=False)


def test_invalid_json_manifest_fails_cleanly(tmp_path):
    d = _bundle(tmp_path / "b")
    (d / "manifest.json").write_text("{not json")
    with pytest.raises(panel.PanelError, match="not valid JSON"):
        panel.load_panel(d, checksums=False)
