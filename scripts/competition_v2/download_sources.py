"""Download the public AtlasShift source universe with resumable, checksum-verified transfers.

Competition track only (``reports/competition_v2/``). Sources and checksums are read
from the vendored, pinned AtlasShift ``sources.json`` (third_party/atlasshift, commit
recorded in third_party/atlasshift_commit.txt). X-Atlas/Orion has no published file
checksum in that repository; it is fetched from Hugging Face at a pinned dataset
revision, whose LFS SHA-256 values are verified by ``huggingface_hub``.

Usage:
    uv run --with huggingface_hub python scripts/competition_v2/download_sources.py \
        [--only k562 cd4 ...]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ROOT / "third_party" / "atlasshift" / "sources.json"
RAW = ROOT / "data" / "raw" / "competition_v2"
LOG = ROOT / "data" / "provenance" / "competition_v2" / "download_log.jsonl"
XATLAS_REPO = "slaf-project/X-Atlas-Orion"
XATLAS_REVISION = "598aa5442a4ee9a8f1a383a843ab5ce142904d04"
XATLAS_LINES = ["HCT116", "HEK293T"]
SAFETY_BYTES = 400 * 1024**3


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(16 * 1024**2), b""):
            h.update(block)
    return h.hexdigest()


def log(record: dict) -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    record["time"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    with LOG.open("a") as f:
        f.write(json.dumps(record) + "\n")
    print(json.dumps(record), flush=True)


def fetch_url(name: str, spec: dict) -> None:
    path = RAW / spec["filename"]
    if path.exists() and path.stat().st_size == spec["bytes"]:
        digest = sha256(path)
        if digest == spec["sha256"]:
            log({"source": name, "status": "present_verified", "sha256": digest})
            return
        raise ValueError(f"{name}: existing file checksum mismatch")
    pending = path.with_name(path.name + ".partial")
    for attempt in range(20):
        result = subprocess.run(
            [
                "curl",
                "-L",
                "--fail",
                "-s",
                "-S",
                # No curl-internal --retry: curl fixes the "-C -" offset once at launch,
                # so an internal retry restarts from that stale offset and truncates the
                # partial file. Retries happen in this loop, which re-resolves the offset.
                "--connect-timeout",
                "30",
                "-C",
                "-",
                "-o",
                str(pending),
                spec["url"],
            ],
        )
        if result.returncode == 0 and pending.stat().st_size == spec["bytes"]:
            break
        log({"source": name, "status": "retry", "attempt": attempt, "rc": result.returncode})
        time.sleep(min(30 * (attempt + 1), 300))
    else:
        raise RuntimeError(f"{name}: download failed")
    digest = sha256(pending)
    if pending.stat().st_size != spec["bytes"] or digest != spec["sha256"]:
        raise ValueError(f"{name}: checksum mismatch {digest}")
    pending.rename(path)
    log({"source": name, "status": "downloaded_verified", "bytes": spec["bytes"], "sha256": digest})


def fetch_xatlas(line: str) -> None:
    from huggingface_hub import snapshot_download

    target = RAW / "xatlas_orion"
    for attempt in range(20):
        try:
            snapshot_download(
                XATLAS_REPO,
                repo_type="dataset",
                revision=XATLAS_REVISION,
                allow_patterns=[f"data/{line}/*", "LICENSE.md", "README.md"],
                local_dir=target,
                max_workers=8,
            )
            break
        except Exception as exc:  # network errors are retried; files resume
            log(
                {
                    "source": f"xatlas_{line}",
                    "status": "retry",
                    "attempt": attempt,
                    "error": repr(exc)[:300],
                }
            )
            time.sleep(min(30 * (attempt + 1), 300))
    else:
        raise RuntimeError(f"xatlas {line}: download failed")
    size = sum(p.stat().st_size for p in (target / "data" / line).rglob("*") if p.is_file())
    log(
        {
            "source": f"xatlas_{line}",
            "status": "downloaded_hf_verified",
            "revision": XATLAS_REVISION,
            "bytes": size,
        }
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", nargs="*", default=None)
    args = parser.parse_args()
    specs = json.loads(SOURCES.read_text())
    jobs = {name: (lambda n=name: fetch_url(n, specs[n])) for name in specs}
    jobs.update({f"xatlas_{line}": (lambda x=line: fetch_xatlas(x)) for line in XATLAS_LINES})
    if args.only:
        jobs = {k: v for k, v in jobs.items() if k in args.only}
    RAW.mkdir(parents=True, exist_ok=True)
    free = shutil.disk_usage(RAW).free
    planned = sum(s["bytes"] for n, s in specs.items() if n in jobs) + 134 * 10**9
    log({"status": "start", "free_bytes": free, "planned_bytes": planned, "jobs": list(jobs)})
    if free - planned < SAFETY_BYTES:
        sys.exit("STOP: disk safety margin would fall below 400 GiB")
    with ThreadPoolExecutor(max_workers=len(jobs)) as pool:
        futures = {name: pool.submit(run) for name, run in jobs.items()}
        failed = []
        for name, future in futures.items():
            try:
                future.result()
            except Exception as exc:
                failed.append(name)
                log({"source": name, "status": "failed", "error": repr(exc)[:500]})
    log({"status": "done", "failed": failed})
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
