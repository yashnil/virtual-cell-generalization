"""Package a ``run_final_panel.py`` prediction into a .vcc. NEVER submits.

    uv run python scripts/competition_v2/package_final_panel.py \
        --output-dir outputs/final/c1 --controls-dir <official_bundle>

Every panel value (contexts, cells per perturbation, gene count) is passed explicitly.
The installed ``vcc`` 0.2.0 defaults to the validation panel (``--contexts A,B,C``,
400 cells, 18,533 genes), and those defaults would reject a D/E/F bundle.

1. Lossless compaction with the vendored ``compact.py`` (unmodified; its
   ``--cells-per-target`` is the panel value).
2. ``vcc prep --dry-run --json`` with the panel's ``gene_names.csv``,
   ``pert_counts.csv``, ``--contexts``, ``--cells-per-pert`` and
   ``--expected-gene-dim``. It must exit 0.
3. ``vcc.prep`` packaging through ``vcc_pack_panel.py``, which runs under the vcc CLI's
   interpreter and reuses the vendored pack.py disk-mapped reader unmodified.
4. A manifest with the SHA-256 of every artefact.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from virtual_cell.competition_v2 import panel  # noqa: E402

UPSTREAM = ROOT / "third_party" / "atlasshift"


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--controls-dir", required=True)
    args = ap.parse_args()
    out = Path(args.output_dir).resolve()
    controls = Path(args.controls_dir).resolve()
    prov = json.loads((out / "provenance.json").read_text())
    if not prov.get("validation", {}).get("all_pass"):
        raise SystemExit("local invariants did not pass; refusing to package")
    pred, compact = out / "prediction.h5ad", out / "prediction_compact.h5ad"
    if panel.sha256(pred) != prov["prediction"]["sha256"]:
        raise SystemExit("prediction.h5ad changed since run_final_panel.py wrote it")
    pid = (prov["panel"].get("panel_id") or "panel").replace("/", "_")
    package = out / f"c1_{pid}.vcc"
    py = UPSTREAM / ".venv" / "bin" / "python"
    pan = prov["panel"]
    contexts = ",".join(pan["contexts"])
    cells, n_genes = str(pan["cells_per_pert"]), str(pan["n_genes"])
    if not compact.exists():
        subprocess.run(
            [str(py), str(UPSTREAM / "compact.py"), str(pred), str(compact),
             "--cells-per-target", cells],
            check=True, cwd=out,
        )  # fmt: skip
    dry = subprocess.run(
        ["vcc", "prep", str(compact), "-g", str(controls / "gene_names.csv"),
         "--perts", str(controls / "pert_counts.csv"), "--contexts", contexts,
         "--cells-per-pert", cells, "--expected-gene-dim", n_genes, "--dry-run", "--json"],
        capture_output=True, text=True, check=False,
    )  # fmt: skip
    try:
        dry_report = json.loads(dry.stdout)
    except json.JSONDecodeError:
        dry_report = {"stdout": dry.stdout[-2000:], "stderr": dry.stderr[-2000:]}
    (out / "vcc_prep_dry_run.json").write_text(json.dumps(dry_report, indent=2) + "\n")
    if dry.returncode != 0:
        print(json.dumps(dry_report, indent=1)[:2000])
        raise SystemExit("vcc prep --dry-run failed; not packaging")
    print("vcc prep --dry-run: exit 0")
    if not package.exists():
        vcc_python = Path(shutil.which("vcc")).read_text().splitlines()[0].removeprefix("#!")
        subprocess.run(
            [vcc_python, str(ROOT / "scripts" / "competition_v2" / "vcc_pack_panel.py"),
             str(compact), "--genes", str(controls / "gene_names.csv"),
             "--perts", str(controls / "pert_counts.csv"), "--contexts", contexts,
             "--cells-per-pert", cells, "--expected-gene-dim", n_genes,
             "--output", str(package), "--scratch-dir", str(out / "scratch")],
            check=True, cwd=out,
        )  # fmt: skip
    manifest = {
        "submitted": False,
        "panel": prov["panel"],
        "prediction": {"path": str(pred), "sha256": panel.sha256(pred)},
        "compact": {"path": str(compact), "sha256": panel.sha256(compact)},
        "vcc_package": {
            "path": str(package),
            "sha256": panel.sha256(package),
            "bytes": package.stat().st_size,
        },
        "vcc_prep_dry_run": dry_report,
        "model_git": prov["git"],
    }
    (out / "package_manifest.json").write_text(json.dumps(manifest, indent=2, default=str) + "\n")
    print(f"packaged {package} (sha256 {manifest['vcc_package']['sha256']}); NOT submitted")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
