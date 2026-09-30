"""Resumable parallel range download of a public repository file, then MD5 verification.

Each 256 MiB chunk is fetched with a ``Range`` GET (re-resolving the repository's
short-lived redirect each time) and written in place into a preallocated file. A chunk is
marked done in ``<file>.parts.json`` only after its full length arrives, so a rerun
resumes. The finished file's MD5 must equal the repository's ``supplied_md5``, and the
download is appended to ``data/provenance/competition_v2/download_log.jsonl``.

    uv run python scripts/competition_v2/c4_download.py <url> <dest> --bytes N --md5 HEX \
        --doi DOI --license LICENSE [--threads 8]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import UTC, datetime
from pathlib import Path

import requests

CHUNK = 256 << 20
ROOT = Path(__file__).resolve().parents[2]
LOG = ROOT / "data" / "provenance" / "competition_v2" / "download_log.jsonl"


def fetch(url: str, dest: Path, lo: int, hi: int) -> None:
    for attempt in range(8):
        try:
            with requests.get(
                url, headers={"Range": f"bytes={lo}-{hi}"}, stream=True, timeout=120
            ) as r:
                if r.status_code != 206:
                    raise OSError(f"HTTP {r.status_code}")
                pos = lo
                with dest.open("r+b") as f:
                    f.seek(lo)
                    for block in r.iter_content(1 << 20):
                        f.write(block)
                        pos += len(block)
                if pos != hi + 1:
                    raise OSError(f"short chunk {pos - lo} of {hi - lo + 1}")
                return
        except Exception as exc:  # noqa: BLE001 - retried, then re-raised
            if attempt == 7:
                raise
            time.sleep(2 + 4 * attempt)
            print(f"retry {lo}: {exc}", flush=True)


def md5(path: Path) -> str:
    h = hashlib.md5()  # noqa: S324 - repository checksum is MD5
    with path.open("rb") as f:
        for block in iter(lambda: f.read(64 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("url")
    p.add_argument("dest")
    p.add_argument("--bytes", type=int, required=True)
    p.add_argument("--md5", required=True)
    p.add_argument("--doi", required=True)
    p.add_argument("--license", required=True)
    p.add_argument("--threads", type=int, default=8)
    a = p.parse_args()
    dest = Path(a.dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    parts_path = dest.with_suffix(dest.suffix + ".parts.json")
    done = set(json.loads(parts_path.read_text())) if parts_path.exists() else set()
    if not dest.exists():
        with dest.open("wb") as f:
            f.truncate(a.bytes)
    if dest.stat().st_size != a.bytes:
        raise SystemExit("destination has the wrong size")
    chunks = [(lo, min(lo + CHUNK, a.bytes) - 1) for lo in range(0, a.bytes, CHUNK)]
    todo = [c for c in chunks if c[0] not in done]
    lock, t0, got = threading.Lock(), time.time(), [0]

    def work(c):
        fetch(a.url, dest, *c)
        with lock:
            done.add(c[0])
            got[0] += c[1] - c[0] + 1
            parts_path.write_text(json.dumps(sorted(done)))
            if len(done) % 20 == 0 or len(done) == len(chunks):
                rate = got[0] / (time.time() - t0) / 1e6
                print(f"{len(done)}/{len(chunks)} chunks, {rate:.1f} MB/s", flush=True)

    with ThreadPoolExecutor(a.threads) as pool:
        list(pool.map(work, todo))
    digest = md5(dest)
    ok = digest == a.md5
    entry = {
        "utc": datetime.now(UTC).isoformat(),
        "url": a.url,
        "path": str(dest.relative_to(ROOT)) if dest.is_relative_to(ROOT) else str(dest),
        "bytes": a.bytes,
        "md5": digest,
        "md5_expected": a.md5,
        "md5_ok": ok,
        "doi": a.doi,
        "license": a.license,
    }
    with LOG.open("a") as f:
        f.write(json.dumps(entry) + "\n")
    print(json.dumps(entry))
    if not ok:
        raise SystemExit("MD5 mismatch")
    parts_path.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
