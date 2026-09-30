"""Effect -> expected composition -> integer counts (the C0-style emission).

Attribution: the mapping (``desired_mean``), the CRISPRi promoter-neighbour cap, the
4-cell control template and the dual-moment integer emitter are AtlasShift's
(``model.py`` / ``predict.py`` / ``prepare.pairs`` @ ``d24ce4f``, MIT). This is our own
implementation of the same algorithm; ``tests/test_competition_v2_c1.py`` requires the
integer output to be identical to upstream for the same seed. No new generator is built
in C1.
"""

from __future__ import annotations

import gzip
import hashlib
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from virtual_cell.competition_v2.atlas import BULK_SCALE

#: Frozen C0 constants.
AMPLITUDE = {"log2fc": 0.6, "bulk_delta": 0.3}
CLIP = 3.0
PROMOTER_FRACTION = 0.15
PROMOTER_WINDOW = 5000
SEED = 20260910
POOL_K = 4


# ---------------------------------------------------------------------------
# effect -> expected composition
# ---------------------------------------------------------------------------


def expected_composition(
    control_probability: np.ndarray,
    effect: np.ndarray,
    *,
    space: str,
    amplitude: float,
    clip: float = CLIP,
) -> np.ndarray:
    """Apply a (shrunk, clipped) effect to the destination control; renormalise rows."""
    e = np.clip(amplitude * effect, -clip, clip)
    if space == "log2fc":
        out = control_probability * np.exp2(e)
    elif space == "bulk_delta":
        out = np.expm1(np.maximum(np.log1p(BULK_SCALE * control_probability) + e, 0))
    else:
        raise ValueError(space)
    total = out.sum(axis=-1, keepdims=True)
    if (total <= 0).any() or not np.isfinite(out).all():
        raise ValueError("invalid expected composition")
    return out / total


def gencode_tss(gtf_path: Path) -> pd.DataFrame:
    """Gene-level TSS table (``+`` strand: start, ``-`` strand: end) from a GENCODE GTF."""
    rows = []
    with gzip.open(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "gene":
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]*)"', f[8]))
            if "gene_name" in attrs:
                rows.append(
                    {
                        "gene": attrs["gene_name"],
                        "chromosome": f[0],
                        "strand": f[6],
                        "tss": int(f[3] if f[6] == "+" else f[4]),
                    }
                )
    return pd.DataFrame(rows)


def promoter_pairs(tss: pd.DataFrame, targets, genes, window: int = PROMOTER_WINDOW):
    """(target, neighbour, distance) for panel genes whose TSS is within ``window`` bp.

    Symbols duplicated in the annotation are dropped entirely, as upstream.
    """
    unique = tss[~tss.gene.duplicated(keep=False)].set_index("gene")
    panel = unique.loc[unique.index.intersection(pd.Index(genes))]
    rows = []
    for t in targets:
        if t not in unique.index:
            continue
        a = unique.loc[t]
        near = panel[(panel.chromosome == a.chromosome) & ((panel.tss - a.tss).abs() <= window)]
        for name, b in near.iterrows():
            if name != t:
                rows.append(
                    {"target": t, "neighbor": name, "distance": abs(int(a.tss) - int(b.tss))}
                )
    return pd.DataFrame(rows, columns=["target", "neighbor", "distance"])


def promoter_cap(
    probability: np.ndarray,
    control: np.ndarray,
    targets,
    genes,
    pairs: pd.DataFrame,
    fraction: float = PROMOTER_FRACTION,
):
    """Cap each promoter neighbour at ``control * (f + (1-f) * ramp(d))``; renormalise.

    ``ramp(d) = clip(log10(max(d, 500) / 500), 0, 1)``: 15 % of control at <= 500 bp,
    rising log-linearly to 100 % at 5 kb. Only ever lowers a value.
    """
    if not 0 < fraction <= 1:
        raise ValueError("invalid remaining fraction")
    ti = {str(t): i for i, t in enumerate(targets)}
    gi = {str(g): i for i, g in enumerate(genes)}
    out = np.asarray(probability, dtype=np.float64).copy()
    changed = []
    for row in pairs.itertuples(index=False):
        if row.target not in ti or row.neighbor not in gi or row.target == row.neighbor:
            continue
        if not 0 <= row.distance <= PROMOTER_WINDOW:
            raise ValueError("unexpected promoter distance")
        i, j = ti[row.target], gi[row.neighbor]
        ramp = np.clip(np.log(max(float(row.distance), 500) / 500) / np.log(10), 0, 1)
        ceiling = control[j] * (fraction + (1 - fraction) * ramp)
        if out[i, j] > ceiling:
            changed.append(
                {
                    "target": row.target,
                    "neighbor": row.neighbor,
                    "distance": float(row.distance),
                    "before": float(out[i, j]),
                    "after_before_normalization": float(ceiling),
                }
            )
            out[i, j] = ceiling
    out /= out.sum(axis=1, keepdims=True)
    return out, changed


# ---------------------------------------------------------------------------
# control template and integer emission
# ---------------------------------------------------------------------------


def seed_for(key: str, seed: int = SEED) -> int:
    return (int.from_bytes(hashlib.sha256(key.encode()).digest()[:4], "little") + seed) % 2**32


def control_template(raw: sparse.spmatrix, cells: int, pool_k: int = POOL_K, seed: int = SEED):
    """Smoothed template: ``cells * pool_k`` controls, sorted by depth, averaged in groups.

    Returns ``(template, depths, mean_cpm_composition, bulk_composition)``.
    """
    raw = sparse.csr_matrix(raw, dtype=np.float64)
    if (
        not np.isfinite(raw.data).all()
        or (raw.data < 0).any()
        or (raw.data != np.floor(raw.data)).any()
    ):
        raise ValueError("controls must be raw nonnegative integer counts")
    library = np.asarray(raw.sum(axis=1)).ravel()
    if (library <= 0).any() or len(library) < cells * pool_k:
        raise ValueError("need positive-depth controls and >= cells * pool_k donors")
    n_genes = raw.shape[1]
    mean = np.zeros(n_genes)
    for left in range(0, len(library), 256):
        mean += (raw[left : left + 256].toarray() / library[left : left + 256, None]).sum(axis=0)
    mean /= len(library)
    bulk = np.asarray(raw.sum(axis=0)).ravel()
    bulk /= bulk.sum()
    chosen = np.random.default_rng(seed).choice(len(library), cells * pool_k, replace=False)
    chosen = chosen[np.argsort(library[chosen], kind="stable")]
    template = raw[chosen].toarray() / library[chosen, None]
    template = template.reshape(cells, pool_k, n_genes).mean(axis=1)
    depths = np.rint(library[chosen].reshape(cells, pool_k).mean(axis=1)).astype(np.int64)
    if depths.min() < 1 or depths.max() > 1_000_000:
        raise ValueError("invalid template depths")
    return template, depths, mean, bulk


def repair_column_totals(integer: np.ndarray, expected: np.ndarray) -> np.ndarray:
    """Move single counts between cells so every gene's total is its rounded expectation,
    keeping every cell's depth fixed."""
    expected = np.asarray(expected, dtype=np.float64)
    if integer.shape != expected.shape or integer.ndim != 2:
        raise ValueError("count and expectation axes differ")
    if not np.isfinite(expected).all() or (expected < 0).any() or (integer < 0).any():
        raise ValueError("invalid expectations")
    depth = integer.sum(axis=1, dtype=np.int64)
    columns = expected.sum(axis=0)
    goal = np.floor(columns).astype(np.int64)
    spare = int(depth.sum() - goal.sum())
    if spare < 0 or spare > len(goal):
        raise ValueError("expected and integer grand totals differ")
    if spare:
        goal[np.argsort(columns - goal, kind="stable")[-spare:]] += 1
    diff = integer.sum(axis=0, dtype=np.int64) - goal
    capacity = np.zeros(len(integer), dtype=np.int64)
    for g in np.flatnonzero(diff > 0):
        excess = int(diff[g])
        up = integer[:, g] - np.floor(expected[:, g]).astype(np.int64)
        rows = np.flatnonzero(up > 0)
        if up[rows].sum() < excess:
            raise ValueError("rounding fell below its integer floor")
        rows = rows[np.argsort(-(integer[rows, g] - expected[rows, g]), kind="stable")]
        if excess <= len(rows):
            integer[rows[:excess], g] -= 1
            capacity[rows[:excess]] += 1
        else:
            for r in rows:
                take = min(excess, int(up[r]))
                integer[r, g] -= take
                capacity[r] += take
                excess -= take
                if not excess:
                    break
    short = np.flatnonzero(diff < 0)
    for g in short[np.argsort(diff[short], kind="stable")]:
        need = int(-diff[g])
        while need:
            rows = np.flatnonzero(capacity > 0)
            if not len(rows):
                raise AssertionError("column repair exhausted row capacity")
            take = min(need, len(rows))
            pick = rows[np.argsort(-(expected[rows, g] - integer[rows, g]), kind="stable")[:take]]
            integer[pick, g] += 1
            capacity[pick] -= 1
            need -= take
    if capacity.any() or (integer < 0).any():
        raise AssertionError("column repair lost counts")
    if not np.array_equal(integer.sum(axis=1, dtype=np.int64), depth):
        raise AssertionError("column repair changed a depth")
    if not np.array_equal(integer.sum(axis=0, dtype=np.int64), goal):
        raise AssertionError("column repair missed its totals")
    return integer


def _project_bulk(bulk: np.ndarray, p: np.ndarray, z: np.ndarray, *, bracket: bool):
    """Project the requested bulk composition onto what depth-tilting can reach.

    With equal depths (``bracket=False``) the caller passes ``bulk = p``; the bisection
    still runs inside ``[0, 1]``, exactly as upstream.
    """
    lower = (z.min() + 0.01 * (1 - z.min())) * p
    upper = (z.max() - 0.01 * (z.max() - 1)) * p
    lo, hi = 0.0, 1.0
    if bracket:
        for _ in range(1024):
            if np.clip(bulk * hi, lower, upper).sum() >= 1:
                break
            hi *= 2
            if not np.isfinite(hi):
                raise ValueError("bulk support cannot meet the feasible bounds")
        else:
            raise ValueError("cannot bracket the bulk projection")
    for _ in range(70):
        mid = (lo + hi) / 2
        if np.clip(bulk * mid, lower, upper).sum() < 1:
            lo = mid
        else:
            hi = mid
    return np.clip(bulk * ((lo + hi) / 2), lower, upper)


def dual_moment_counts(
    template, probability, bulk_probability, *, depths, seed, iterations=100, tolerance=2e-4
) -> np.ndarray:
    """Integer cells whose mean composition and depth-weighted composition both hit targets.

    Rescale the template columns to the mean-CPM target, tilt each gene along relative
    depth ``z`` toward the pseudobulk target, alternate to convergence, floor, place each
    cell's residual counts by systematic sampling, then repair gene totals.
    """
    template = np.asarray(template)
    p = np.asarray(probability, dtype=np.float64)
    requested = np.asarray(bulk_probability, dtype=np.float64)
    raw_depths = np.asarray(depths)
    if template.ndim != 2 or not all(template.shape):
        raise ValueError("template must be a nonempty cell-by-gene matrix")
    if p.shape != (template.shape[1],) or requested.shape != p.shape:
        raise ValueError("moment axes differ from the template")
    for v in (template, p, requested, raw_depths):
        if not np.isfinite(v).all() or (v < 0).any():
            raise ValueError("inputs must be finite and nonnegative")
    if (
        (raw_depths != np.floor(raw_depths)).any()
        or (raw_depths < 1).any()
        or (raw_depths > 1_000_000).any()
    ):
        raise ValueError("depths must be integers in [1, 1e6]")
    if (template.sum(axis=1) <= 0).any():
        raise ValueError("template rows must have positive mass")
    if not np.isclose(p.sum(), 1, rtol=0, atol=1e-8) or not np.isclose(
        requested.sum(), 1, rtol=0, atol=1e-8
    ):
        raise ValueError("moment compositions must sum to one")
    if (
        type(iterations) is not int
        or iterations < 1
        or not np.isfinite(tolerance)
        or tolerance <= 0
    ):
        raise ValueError("invalid iterations or tolerance")

    x = np.asarray(template, dtype=np.float64).copy()
    depths = np.broadcast_to(np.asarray(depths, dtype=np.int64), (len(x),))
    z = depths / depths.mean()
    desired = len(x) * p
    varied = np.ptp(depths) != 0
    bulk = _project_bulk(requested.copy() if varied else p.copy(), p, z, bracket=varied)
    if float(np.abs(bulk - requested).sum()) > 0.03:
        raise ValueError(f"requested bulk projection too large: {np.abs(bulk - requested).sum()}")
    ratio = np.divide(bulk, p, out=np.ones_like(p), where=p > 0)
    x /= np.maximum(x.sum(axis=1, keepdims=True), 1e-30)
    missing = (x.sum(axis=0) == 0) & (desired > 0)
    x[:, missing] = p[missing]
    x = 0.999 * x + 0.001 * p[None, :]
    for step in range(iterations):
        col = x.sum(axis=0)
        x *= np.divide(desired, col, out=np.zeros_like(col), where=col > 0)
        first = z @ x
        mean = np.divide(first, desired, out=np.ones_like(first), where=desired > 0)
        second = z * z @ x
        var = np.maximum(
            np.divide(second, desired, out=np.zeros_like(second), where=desired > 0) - mean * mean,
            0,
        )
        tilt = np.divide(ratio - mean, var, out=np.zeros_like(mean), where=var > 1e-12)
        tilt = np.clip(
            tilt,
            -0.95 / np.maximum(z.max() - mean, 1e-12),
            0.95 / np.maximum(mean - z.min(), 1e-12),
        )
        x *= 1 + tilt[None, :] * (z[:, None] - mean[None, :])
        x /= np.maximum(x.sum(axis=1, keepdims=True), 1e-30)
        if step % 5 == 4:
            cpm_err = float(np.abs(x.mean(axis=0) - p).sum())
            bulk_err = float(np.abs(z @ x / len(x) - bulk).sum())
            if max(cpm_err, bulk_err) < tolerance:
                break
    cpm_err = float(np.abs(x.mean(axis=0) - p).sum())
    bulk_err = float(np.abs(z @ x / len(x) - bulk).sum())
    if max(cpm_err, bulk_err) > 1e-3:
        raise ValueError(f"moment fitting failed: CPM {cpm_err}, bulk {bulk_err}")
    expected = x * depths[:, None]
    integer = np.floor(expected).astype(np.int32)
    frac = expected - integer
    rng = np.random.default_rng(seed)
    for i in range(len(x)):
        residual = int(depths[i]) - int(integer[i].sum())
        if residual:
            cum = np.cumsum(frac[i])
            cum *= residual / cum[-1]
            np.add.at(
                integer[i],
                np.searchsorted(cum, np.arange(residual) + rng.random(), side="right"),
                1,
            )
    integer = repair_column_totals(integer, expected)
    if (integer < 0).any() or not np.array_equal(integer.sum(axis=1), depths):
        raise AssertionError("emission lost row depths")
    return integer


def expected_moments(
    effects: dict[str, np.ndarray],
    controls: dict[str, np.ndarray],
    targets,
    genes,
    *,
    shrink: np.ndarray | None = None,
    pairs: pd.DataFrame | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Fused effects -> the two expected compositions the emitter must hit.

    ``shrink`` (per target, C1b) multiplies the fused effect before the frozen amplitude
    and clip; ``pairs`` enables the promoter-neighbour cap. Returns ``(p_cpm, p_bulk)``.
    """
    out = []
    for space in ("log2fc", "bulk_delta"):
        effect = effects[space]
        if shrink is not None:
            effect = (np.asarray(shrink, dtype=np.float64)[:, None] * effect).astype(effect.dtype)
        prob = expected_composition(
            controls[space], effect, space=space, amplitude=AMPLITUDE[space]
        )
        if pairs is not None:
            prob, _ = promoter_cap(prob, controls[space], targets, genes, pairs)
        out.append(prob)
    return out[0], out[1]


class CountWriter:
    """Append integer count blocks to an ``.h5ad`` with int64 CSR pointers.

    Layout follows the C0 bundle (context-major, then target, then cell), so the obs
    order matches ``pert_counts.csv`` within each context.
    """

    def __init__(self, path: Path, targets, genes, contexts, cells: int):
        import anndata as ad
        import h5py

        self.path = Path(path)
        self.rows = self.nnz = 0
        self.nobs = len(targets) * len(contexts) * cells
        self.ngenes = len(genes)
        obs = pd.DataFrame(
            {
                "target_gene": np.tile(np.repeat(targets, cells), len(contexts)),
                "context": np.repeat(contexts, len(targets) * cells),
            },
            index=[f"{c}_{t}_{i}" for c in contexts for t in targets for i in range(cells)],
        )
        ad.AnnData(
            sparse.csr_matrix((self.nobs, len(genes)), dtype=np.float32),
            obs=obs,
            var=pd.DataFrame(index=list(genes)),
        ).write_h5ad(self.path)
        self.file = h5py.File(self.path, "r+")
        group = self.file["X"]
        for name in ("data", "indices", "indptr"):
            del group[name]
        kw = dict(shape=(0,), maxshape=(None,), chunks=(1 << 20,), compression="lzf", shuffle=True)
        self.data = group.create_dataset("data", dtype="float32", **kw)
        self.indices = group.create_dataset("indices", dtype="int32", **kw)
        self.indptr = group.create_dataset("indptr", shape=(self.nobs + 1,), dtype="int64")
        self.indptr[0] = 0

    def append(self, block: sparse.csr_matrix) -> None:
        block = sparse.csr_matrix(block)
        if block.shape[1] != self.ngenes or self.rows + block.shape[0] > self.nobs:
            raise ValueError("unexpected block shape")
        block.eliminate_zeros()
        d = block.data
        if not np.isfinite(d).all() or (d < 0).any() or (d != np.floor(d)).any():
            raise ValueError("expected finite nonnegative integer counts")
        if (np.asarray(block.sum(axis=1)).ravel() > 1_000_000).any():
            raise ValueError("cell depth exceeds 1,000,000")
        end = self.nnz + block.nnz
        if end > 4_750_000_000:
            raise ValueError("sparse element limit exceeded")
        self.data.resize((end,))
        self.indices.resize((end,))
        self.data[self.nnz : end] = d.astype(np.float32)
        self.indices[self.nnz : end] = block.indices
        self.indptr[self.rows + 1 : self.rows + block.shape[0] + 1] = (
            block.indptr[1:].astype(np.int64) + self.nnz
        )
        self.rows += block.shape[0]
        self.nnz = end

    def close(self) -> None:
        complete = self.rows == self.nobs and int(self.indptr[-1]) == self.nnz
        self.file.close()
        if not complete:
            raise ValueError("incomplete prediction output")
