"""Panel-aware ``vcc.prep`` packaging adapter. Runs under the vcc CLI's own interpreter.

The vendored ``third_party/atlasshift/pack.py`` calls ``vcc.prep.run_prep`` with its
defaults: contexts A,B,C, 400 cells per perturbation and 18,533 genes. Those are
validation-panel values. This adapter reuses pack.py's disk-mapped reader **unmodified**
(imported by path) and passes the panel's own contexts, cells per perturbation and gene
count. It never submits.

    <vcc interpreter> scripts/competition_v2/vcc_pack_panel.py PRED.h5ad \
        --genes gene_names.csv --perts pert_counts.csv --contexts D,E,F \
        --cells-per-pert 400 --expected-gene-dim 18533 --output out.vcc --scratch-dir s
"""

from __future__ import annotations

import argparse
import importlib.util
import tempfile
from importlib.metadata import version
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PACK = ROOT / "third_party" / "atlasshift" / "pack.py"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("prediction", type=Path)
    ap.add_argument("--genes", type=Path, required=True)
    ap.add_argument("--perts", type=Path, required=True)
    ap.add_argument("--contexts", required=True)
    ap.add_argument("--cells-per-pert", type=int, required=True)
    ap.add_argument("--expected-gene-dim", type=int, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--scratch-dir", type=Path, required=True)
    args = ap.parse_args()
    import vcc.prep as prep

    spec = importlib.util.spec_from_file_location("atlasshift_pack", PACK)
    pack = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pack)
    if args.output.exists():
        raise FileExistsError(args.output)
    args.scratch_dir.mkdir(parents=True, exist_ok=True)
    before = pack.digest(args.prediction)
    original = prep.read_h5ad
    with tempfile.TemporaryDirectory(prefix="vcc_", dir=args.scratch_dir) as tmp:

        def reader(path):
            if Path(path).resolve() != args.prediction.resolve():
                raise ValueError("unexpected prediction input")
            return pack.mapped_adata(path, tmp)

        prep.read_h5ad = reader
        try:
            prep.run_prep(
                input_path=str(args.prediction),
                output_path=str(args.output),
                genes_path=str(args.genes),
                perts_path=str(args.perts),
                cells_per_pert=args.cells_per_pert,
                required_contexts=tuple(c for c in args.contexts.split(",") if c),
                expected_gene_dim=args.expected_gene_dim,
            )
        finally:
            prep.read_h5ad = original
    if pack.digest(args.prediction) != before:
        raise ValueError("prediction changed during packaging")
    print(f"Saved {args.output.name} with vcc-cli {version('vcc-cli')}")


if __name__ == "__main__":
    main()
