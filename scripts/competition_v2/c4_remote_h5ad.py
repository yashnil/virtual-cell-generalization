"""Inspect a remote .h5ad through HTTP range requests, without downloading it (C4 audit).

Reads only metadata and small arrays: obs columns (categorical categories + codes),
var names, the X encoding and a small sample of X values. The byte ranges are fetched
with ``Range`` GETs against the original repository URL, so the audit reads the
authoritative file, not a mirror.

    uv run python scripts/competition_v2/c4_remote_h5ad.py <url> --out <json> [--obs COL ...]
"""

from __future__ import annotations

import argparse
import io
import json
from pathlib import Path

import h5py
import numpy as np
import requests

BLOCK = 1 << 22  # 4 MiB


class RangeFile(io.RawIOBase):
    """A seekable read-only file over HTTP range GETs, with a block cache."""

    def __init__(self, url: str, size: int | None = None):
        self.url, self.pos, self.cache = url, 0, {}
        self.session = requests.Session()
        self.fetched = 0
        if size is None:
            r = self.session.get(url, headers={"Range": "bytes=0-0"}, timeout=60)
            size = int(r.headers["Content-Range"].split("/")[1])
        self.size = size

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, offset, whence=0):
        self.pos = {0: offset, 1: self.pos + offset, 2: self.size + offset}[whence]
        return self.pos

    def _block(self, i: int) -> bytes:
        if i not in self.cache:
            lo, hi = i * BLOCK, min((i + 1) * BLOCK, self.size) - 1
            for _ in range(5):
                r = self.session.get(self.url, headers={"Range": f"bytes={lo}-{hi}"}, timeout=120)
                if r.status_code == 206 and len(r.content) == hi - lo + 1:
                    break
            else:
                raise OSError(f"range fetch failed at block {i}")
            self.cache[i] = r.content
            self.fetched += len(r.content)
        return self.cache[i]

    def readinto(self, b):
        n = min(len(b), self.size - self.pos)
        if n <= 0:
            return 0
        out, pos = bytearray(), self.pos
        while len(out) < n:
            i = pos // BLOCK
            blk = self._block(i)
            off = pos - i * BLOCK
            take = blk[off : off + n - len(out)]
            out += take
            pos += len(take)
        b[:n] = out
        self.pos += n
        return n


def _strings(ds) -> np.ndarray:
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in ds[()]])


def read_obs_column(f: h5py.File, col: str):
    node = f["obs"][col]
    if isinstance(node, h5py.Group) and "categories" in node:
        cats = _strings(node["categories"])
        codes = node["codes"][()]
        return cats, codes
    return None, _strings(node) if node.dtype.kind in "OSU" else node[()]


def inspect(url: str, obs_cols: list[str], size: int | None = None) -> dict:
    fh = RangeFile(url, size)
    with h5py.File(fh, "r") as f:
        obs = f["obs"]
        info = {
            "url": url,
            "bytes": fh.size,
            "root_keys": list(f.keys()),
            "obs_columns": list(obs.keys()),
            "n_obs": int(f["obs"][obs.attrs.get("_index", "_index")].shape[0])
            if obs.attrs.get("_index", "_index") in obs
            else None,
        }
        var = f["var"]
        vindex = var.attrs.get("_index", "_index")
        info["var_columns"] = list(var.keys())
        genes = _strings(var[vindex])
        info["n_vars"] = len(genes)
        info["var_names_head"] = genes[:10].tolist()
        x = f["X"]
        if isinstance(x, h5py.Group):
            info["X_encoding"] = dict(x.attrs)
            data = x["data"]
            info["X_dtype"] = str(data.dtype)
            info["X_nnz"] = int(data.shape[0])
            sample = data[: min(200_000, data.shape[0])]
        else:
            info["X_encoding"] = {"dense": list(x.shape)}
            info["X_dtype"] = str(x.dtype)
            sample = x[:50].ravel()
        sample = np.asarray(sample, dtype=np.float64)
        info["X_sample"] = {
            "integer_fraction": float(np.mean(sample == np.round(sample))),
            "min": float(sample.min()),
            "max": float(sample.max()),
            "median": float(np.median(sample)),
        }
        info["layers"] = list(f["layers"].keys()) if "layers" in f else []
        info["uns_keys"] = list(f["uns"].keys()) if "uns" in f else []
        columns = {}
        for col in obs_cols:
            if col not in obs:
                continue
            cats, codes = read_obs_column(f, col)
            if cats is not None:
                counts = np.bincount(codes[codes >= 0], minlength=len(cats))
                columns[col] = {"categories": cats.tolist(), "counts": counts.tolist()}
            else:
                vals, counts = np.unique(codes, return_counts=True)
                columns[col] = {"categories": [str(v) for v in vals], "counts": counts.tolist()}
        info["obs_values"] = columns
        info["genes"] = genes.tolist()
    info["bytes_fetched"] = fh.fetched
    return info


def _jsonable(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, bytes):
        return o.decode()
    return str(o)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("url")
    parser.add_argument("--out", required=True)
    parser.add_argument("--obs", nargs="*", default=[])
    parser.add_argument("--list-only", action="store_true")
    args = parser.parse_args()
    if args.list_only:
        fh = RangeFile(args.url)
        with h5py.File(fh, "r") as f:
            print(
                json.dumps(
                    {
                        "obs": list(f["obs"].keys()),
                        "var": list(f["var"].keys()),
                        "root": list(f.keys()),
                    }
                )
            )
        return
    info = inspect(args.url, args.obs)
    Path(args.out).write_text(json.dumps(info, default=_jsonable))
    summary = {k: v for k, v in info.items() if k not in ("genes", "obs_values")}
    summary["obs_value_counts"] = {k: len(v["categories"]) for k, v in info["obs_values"].items()}
    print(json.dumps(summary, indent=1, default=str))


if __name__ == "__main__":
    main()
