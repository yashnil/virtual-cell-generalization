"""Provenance gate for the C1 license-clean phase (read-only; any mismatch is fatal).

Checks, in order:

1. every frozen research manifest (``data/provenance/scperteval/*_freeze.txt``) and
   the V1 submission digests (``data/provenance/arc_submission_v1_sha256.txt``);
2. every raw-data checksum file (``data/provenance/**/*_sha256.txt``);
3. the frozen V1 tier file (hash pinned in ``tests/test_competition_v2.py``);
4. the vendored AtlasShift files, the C0 prepared statistics and the promoter pairs
   against ``outputs/competition_v2/atlasshift_c0/c0_manifest.json``;
5. with ``--full``: the competition-v2 URL downloads against
   ``data/provenance/competition_v2/download_log.jsonl``, the C0 prediction / compact /
   ``.vcc`` digests, and every X-Atlas/Orion file against its Hugging Face LFS etag.

    uv run python scripts/competition_v2/verify_c1_state.py [--full] [--out PATH]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PROV = ROOT / "data" / "provenance"
C0 = ROOT / "outputs" / "competition_v2" / "atlasshift_c0"
RAW_V2 = ROOT / "data" / "raw" / "competition_v2"
V1_TIERS_SHA = "2280f40b2ce874eef335dab169216d7aad476a887dd49b04845a5a5c4b3ed6f4"
DOWNLOAD_FILES = {
    "k562": "K562_gwps_raw_singlecell_01.h5ad",
    "cd4": "GWCD4i.DE_stats.h5ad",
    "gencode47": "gencode.v47.annotation.gtf.gz",
    "h1-train": "adata_Training.h5ad",
    "h1-targets": "pert_counts_Training.csv",
    "h1-genes": "gene_names.csv",
    "h1-validation": "adata_Validation.h5ad",
    "h1-test": "adata_Test.h5ad",
}


def sha256(path: Path, algo: str = "sha256") -> str:
    h = hashlib.new(algo)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024**2), b""):
            h.update(block)
    return h.hexdigest()


def git_blob_sha1(path: Path) -> str:
    data = path.read_bytes()
    return hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest()


def shasum_manifest(manifest: Path) -> dict:
    proc = subprocess.run(
        ["shasum", "-a", "256", "-c", str(manifest)],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    lines = [ln for ln in proc.stdout.splitlines() if ln.strip()]
    bad = [ln for ln in lines if not ln.endswith(": OK")] + proc.stderr.splitlines()
    return {
        "manifest": str(manifest.relative_to(ROOT)),
        "files": len(lines),
        "failed": len(bad),
        "detail": bad[:3],
    }


def check_c0_small() -> list[dict]:
    manifest = json.loads((C0 / "c0_manifest.json").read_text())
    rows = []
    for name, want in manifest["code"]["upstream_files_sha256"].items():
        path = ROOT / "third_party" / "atlasshift" / name
        if name == ".gitignore" and not path.exists():
            continue  # upstream .gitignore was not vendored
        rows.append({"item": f"third_party/atlasshift/{name}", "ok": sha256(path) == want})
    for name, want in manifest["source_data"]["prepared_statistics_sha256"].items():
        rows.append({"item": f"C0 data/{name}", "ok": sha256(C0 / "data" / name) == want})
    rows.append(
        {
            "item": "C0 data/official_pairs.csv",
            "ok": sha256(C0 / "data" / "official_pairs.csv")
            == manifest["source_data"]["official_pairs_sha256"],
        }
    )
    return rows


def check_c0_large() -> list[dict]:
    manifest = json.loads((C0 / "c0_manifest.json").read_text())
    rows = []
    for key in ["prediction", "compact", "vcc_package"]:
        entry = manifest[key]
        rows.append({"item": entry["path"], "ok": sha256(ROOT / entry["path"]) == entry["sha256"]})
    return rows


def check_downloads() -> list[dict]:
    log = [json.loads(x) for x in (PROV / "competition_v2" / "download_log.jsonl").open()]
    want = {r["source"]: r for r in log if r.get("status") == "downloaded_verified"}
    rows = []
    for source, filename in DOWNLOAD_FILES.items():
        path = RAW_V2 / filename
        ok = path.stat().st_size == want[source]["bytes"] and sha256(path) == want[source]["sha256"]
        rows.append({"item": f"data/raw/competition_v2/{filename}", "ok": ok})
    return rows


def check_xatlas() -> list[dict]:
    """Every X-Atlas file against the etag huggingface_hub recorded when it downloaded it."""
    base = RAW_V2 / "xatlas_orion"
    meta_root = base / ".cache" / "huggingface" / "download"
    rows = []
    n = n_ok = 0
    total = 0
    bad = []
    for meta in sorted(meta_root.rglob("*.metadata")):
        rel = meta.relative_to(meta_root).with_suffix("")
        path = base / rel
        lines = meta.read_text().splitlines()
        commit, etag = lines[0].strip(), lines[1].strip()
        if not path.exists():
            continue  # partial/lock remnants of superseded lance versions
        n += 1
        total += path.stat().st_size
        got = sha256(path) if len(etag) == 64 else git_blob_sha1(path)
        if got == etag and commit.startswith("598aa544"):
            n_ok += 1
        else:
            bad.append(str(rel))
    rows.append(
        {
            "item": "X-Atlas/Orion files vs HF etags (revision 598aa544)",
            "files": n,
            "bytes": total,
            "ok": n_ok == n and n > 0,
            "failed": bad[:5],
        }
    )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    started = time.time()
    freezes = [shasum_manifest(m) for m in sorted((PROV / "scperteval").glob("*_freeze.txt"))]
    freezes.append(shasum_manifest(PROV / "arc_submission_v1_sha256.txt"))
    raw = [
        shasum_manifest(m)
        for m in sorted(PROV.rglob("*_sha256.txt"))
        if m.name != "arc_submission_v1_sha256.txt"
    ]
    tiers_ok = sha256(ROOT / "data" / "splits" / "arc_target_support_v1.csv") == V1_TIERS_SHA
    c0 = check_c0_small()
    result = {
        "time": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "freeze_manifests": freezes,
        "raw_checksum_files": raw,
        "v1_tier_file_ok": tiers_ok,
        "c0_code_and_statistics": c0,
    }
    if args.full:
        result["competition_v2_downloads"] = check_downloads()
        result["c0_bundle"] = check_c0_large()
        result["xatlas"] = check_xatlas()
    failures = (
        sum(f["failed"] for f in freezes + raw)
        + (not tiers_ok)
        + sum(
            not r["ok"]
            for key in result
            if key.startswith(("c0", "competition", "xatlas"))
            for r in result[key]
        )
    )
    result["n_freeze_manifests"] = len(freezes) - 1
    result["n_freeze_digests"] = sum(f["files"] for f in freezes)
    result["n_raw_checksum_files"] = len(raw)
    result["n_raw_digests"] = sum(f["files"] for f in raw)
    result["failures"] = int(failures)
    result["elapsed_seconds"] = round(time.time() - started, 1)
    text = json.dumps(result, indent=2)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text)
    print(text)
    if failures:
        sys.exit("PROVENANCE CHECK FAILED")


if __name__ == "__main__":
    main()
