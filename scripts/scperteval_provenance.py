"""Verify the downloaded scPertEval files and write provenance records.

Checks each file's byte size and upstream MD5 against the values audited in
``reports/scperteval_four_context_data_spec.md`` (recorded before download),
computes a local SHA-256, and writes one provenance record per dataset plus a
checksum file in the same form as the Arc control bundle.

Read-only with respect to the ``.h5ad`` files.

Usage::

    uv run python scripts/scperteval_provenance.py
"""

from __future__ import annotations

import argparse
import json
from datetime import date
from pathlib import Path

from virtual_cell.data import scperteval

REPO_ROOT = Path(__file__).resolve().parents[1]

SCPERTEVAL_COMMIT = "4685f11927e887745737600170da7a655b727553"
SCPERTEVAL_REPO = "https://github.com/Virtual-Cell-Research-Community/scPertEval"

INHERITED = (
    "CRISPRi Perturb-seq (10x 3'). Replogle 2022 essential-gene panel for K562/RPE1 "
    "(not the genome-wide K562 screen); Nadig 2025 companion screens for HepG2/Jurkat. "
    "Guide-to-cell assignment, doublet/quality calls and upstream cell filtering as "
    "deposited. The >=30 cells-per-perturbation floor is inherited from upstream. "
    "Gene detection differs per dataset."
)
SCPERTEVAL_PREP = (
    "1) Label cleaning: one perturbation label per cell, all non-targeting controls "
    "collapsed to the literal label 'control', guide/plasmid suffixes stripped, "
    "combinations joined with '+', cells with a missing label dropped. "
    "2) Log-normalisation: scanpy normalize_total(target_sum=1e4) then log1p. "
    "3) Light QC: filter_cells(min_genes=200), filter_genes(min_cells=3). "
    "4) Trim to X + obs['perturbation'] + gene names; sparse float32, gzip. "
    "No HVG selection, no scaling, no PCA, no batch correction. "
    "Raw counts are NOT retained in the hosted files."
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=REPO_ROOT / scperteval.DATA_SUBDIR)
    parser.add_argument(
        "--out-dir", type=Path, default=REPO_ROOT / "data" / "provenance" / "scperteval"
    )
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    missing = scperteval.missing_files(args.data_dir)
    if missing:
        raise SystemExit(f"Missing files: {[str(p) for p in missing]}")

    retrieved = date.today().isoformat()
    records, checksum_lines, failures = [], [], []

    for ds in scperteval.DATASETS:
        path = ds.path(args.data_dir)
        print(f"hashing {ds.name} ...", flush=True)
        digests = file_or_fail(path)

        size_ok = digests["size_bytes"] == ds.expected_bytes
        md5_ok = digests["md5_base64"] == ds.md5_base64
        if not size_ok:
            failures.append(f"{ds.name}: size {digests['size_bytes']} != {ds.expected_bytes}")
        if not md5_ok:
            failures.append(f"{ds.name}: MD5 {digests['md5_base64']} != {ds.md5_base64}")

        record = {
            "dataset": ds.name,
            "cell_line": ds.cell_line,
            "filename": ds.filename,
            "source_url": ds.url,
            "source_project": SCPERTEVAL_REPO,
            "scperteval_commit": SCPERTEVAL_COMMIT,
            "retrieval_date": retrieved,
            "file_size_bytes": digests["size_bytes"],
            "expected_size_bytes": ds.expected_bytes,
            "size_matches": size_ok,
            "upstream_md5_base64": ds.md5_base64,
            "local_md5_base64": digests["md5_base64"],
            "md5_matches": md5_ok,
            "local_sha256": digests["sha256"],
            "publication": ds.publication,
            "publication_doi": ds.doi,
            "accession": ds.accession,
            "inherited_preprocessing": INHERITED,
            "scperteval_preprocessing": SCPERTEVAL_PREP,
            "local_path": str(path.relative_to(REPO_ROOT)),
        }
        records.append(record)
        (args.out_dir / f"{ds.name}.json").write_text(json.dumps(record, indent=2))
        checksum_lines.append(f"{digests['sha256']}  {path.relative_to(REPO_ROOT)}")
        print(
            f"  size={digests['size_bytes']:,} ({'OK' if size_ok else 'MISMATCH'})  "
            f"md5={'OK' if md5_ok else 'MISMATCH'}  sha256={digests['sha256'][:16]}..."
        )

    (args.out_dir / "scperteval_sha256.txt").write_text("\n".join(checksum_lines) + "\n")
    (args.out_dir / "manifest.json").write_text(
        json.dumps(
            {
                "source_project": SCPERTEVAL_REPO,
                "scperteval_commit": SCPERTEVAL_COMMIT,
                "retrieval_date": retrieved,
                "bucket": "gs://scperteval/processed/",
                "https_base": scperteval.BASE_URL,
                "contexts": list(scperteval.CONTEXTS),
                "total_bytes": sum(r["file_size_bytes"] for r in records),
                "datasets": records,
            },
            indent=2,
        )
    )

    if failures:
        raise SystemExit("CHECKSUM FAILURES:\n  " + "\n  ".join(failures))
    print(f"\nAll {len(records)} files verified. Provenance written to {args.out_dir}")


def file_or_fail(path: Path) -> dict:
    try:
        return scperteval.file_digests(path)
    except OSError as exc:  # pragma: no cover - surfaced to the user directly
        raise SystemExit(f"Could not read {path}: {exc}") from exc


if __name__ == "__main__":
    main()
