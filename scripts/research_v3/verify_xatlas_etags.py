"""Re-verify every local X-Atlas/Orion file against its cached Hugging Face etag.

Etags of 64 hex chars are LFS SHA-256 of the content; 40 hex chars are git blob
SHA-1 (``sha1(b"blob <size>\\0" + content)``). Writes a per-file manifest and a
summary; exits non-zero on any mismatch or missing file.

Reproduce: ``uv run python scripts/research_v3/verify_xatlas_etags.py``
"""

from __future__ import annotations

import hashlib
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "data" / "raw" / "competition_v2" / "xatlas_orion"
META = ROOT / ".cache" / "huggingface" / "download"
OUT = REPO / "data" / "provenance" / "research_v3"
REVISION = "598aa5442a4ee9a8f1a383a843ab5ce142904d04"


def check(meta: Path) -> dict:
    rel = meta.relative_to(META).with_suffix("")
    commit, etag, _ = meta.read_text().splitlines()[:3]
    path = ROOT / rel
    if not path.exists():
        return {"file": str(rel), "status": "missing"}
    size = path.stat().st_size
    if len(etag) == 64:
        h = hashlib.sha256()
    else:
        h = hashlib.sha1()
        h.update(f"blob {size}\0".encode())
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 24), b""):
            h.update(block)
    ok = h.hexdigest() == etag and commit == REVISION
    return {
        "file": str(rel),
        "bytes": size,
        "etag": etag,
        "kind": "sha256" if len(etag) == 64 else "git-sha1",
        "status": "ok" if ok else "mismatch",
    }


def main() -> None:
    metas = sorted(META.rglob("*.metadata"))
    with ThreadPoolExecutor(8) as ex:
        rows = list(ex.map(check, metas))
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "xatlas_etag_manifest.tsv", "w") as fh:
        fh.write("status\tkind\tetag\tbytes\tfile\n")
        for r in rows:
            fh.write(
                "\t".join(str(r.get(c, "")) for c in ("status", "kind", "etag", "bytes", "file"))
                + "\n"
            )
    local = {
        str(p.relative_to(ROOT)) for p in ROOT.rglob("*") if p.is_file() and ".cache" not in p.parts
    }
    summary = {
        "revision": REVISION,
        "n_metadata": len(rows),
        "n_ok": sum(r["status"] == "ok" for r in rows),
        "n_mismatch": sum(r["status"] == "mismatch" for r in rows),
        "n_missing": sum(r["status"] == "missing" for r in rows),
        "local_files_without_metadata": sorted(local - {r["file"] for r in rows}),
        "total_bytes": sum(r.get("bytes", 0) or 0 for r in rows),
    }
    (OUT / "xatlas_etag_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    sys.exit(
        0 if summary["n_ok"] == len(rows) and not summary["local_files_without_metadata"] else 1
    )


if __name__ == "__main__":
    main()
