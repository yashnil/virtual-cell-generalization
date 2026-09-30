"""Panel-agnostic C1: the frozen C1 build, driven by an official panel and the registry.

This reproduces ``scripts/competition_v2/build_c1_candidate.py --stage emit`` exactly,
with every validation-panel constant replaced by a panel value:

* targets, genes and cells per perturbation come from ``pert_counts.csv``,
  ``gene_names.csv`` and ``manifest.json``;
* contexts come from the manifest, in manifest order, and the control template seed is
  ``SEED + context index``, exactly as the validation build used for A/B/C;
* sources come from the registry's qualified C1 sources, in ``c1_order``, with their
  registered fusion weights.

C1 behaviour is unchanged. The same functions (``fusion.fused_effects``,
``generator.expected_moments``, ``control_template``, ``dual_moment_counts`` and
``CountWriter``) are called with the same arguments. A target with no qualified source
gets the frozen C1 fallback automatically (a zero fused effect, so the context control
composition with the promoter cap). No new fallback is introduced.

The one panel-dependent part of C1 is the frozen retained-target definition: the
official targets, then the H1 2025 training targets. It is the row set of the K562 and
CD4 statistics and the set they are centred over. :func:`prepare_sources` reuses the
frozen statistics when their rows already equal it (the validation panel), and
otherwise rebuilds them from raw with the unchanged C1 preparation code.
"""

from __future__ import annotations

import json
import multiprocessing as mp
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from virtual_cell.competition_v2 import fusion, generator, sources
from virtual_cell.competition_v2.panel import Panel, Registry, list_sha, resolve, sha256

_G: dict = {}


# ----------------------------------------------------------------------------- sources


def retained_for_panel(panel: Panel, registry: Registry) -> list[str]:
    extra = pd.read_csv(resolve(registry.raw["retained_extra_targets"])).target_gene.astype(str)
    return sources.retained_targets([str(t) for t in panel.targets], extra.tolist())


def _stats_rows(entry, path: Path) -> list[str]:
    with np.load(path) as d:
        return [str(t) for t in d["targets"]]


def prepare_sources(panel: Panel, registry: Registry, out_dir: Path, log=print) -> dict:
    """Statistics paths for every qualified C1 source, reused or rebuilt for this panel.

    Returns ``{name: {"path", "sha256", "action", "rows"}}``.
    """
    out_dir = Path(out_dir)
    retained = retained_for_panel(panel, registry)
    result = {}
    for e in registry.c1_sources():
        frozen = resolve(e.info["frozen_statistics"])
        if not e.info.get("panel_dependent"):
            path, action = frozen, "reused (panel-independent)"
        elif frozen.exists() and _stats_rows(e, frozen) == retained:
            path, action = frozen, "reused (rows equal the panel's retained set)"
        else:
            path = out_dir / "sources" / f"{e.name}_panel_statistics.npz"
            path.parent.mkdir(parents=True, exist_ok=True)
            if path.exists() and _stats_rows(e, path) == retained:
                action = "reused (prepared earlier for this panel)"
            else:
                t0 = time.time()
                raw = resolve(e.info["raw_path"])
                if e.info["prepare"] == "k562":
                    stats, _cells = sources.prepare_k562(raw, retained, log=log)
                elif e.info["prepare"] == "cd4":
                    stats = sources.prepare_cd4(raw, retained)
                else:
                    raise ValueError(f"{e.name}: no preparation rule for a panel-dependent source")
                np.savez_compressed(path, **stats)
                action = f"prepared from raw for this panel ({time.time() - t0:.0f}s)"
        if e.info.get("panel_dependent") and _stats_rows(e, path) != retained:
            raise AssertionError(f"{e.name}: statistics rows differ from the retained set")
        result[e.name] = {"path": str(path), "sha256": sha256(path), "action": action}
        log(f"source {e.name}: {action}")
    return result


def load_c1_inputs(registry: Registry, prepared: dict):
    """``(cell_sources, weights, cd4_stats, cd4_weight, stats_by_name)`` in c1_order."""
    cell, weights, cd4, cd4_w, by_name = [], [], None, 0.0, {}
    for e in registry.c1_sources():
        path = Path(prepared[e.name]["path"])
        if e.kind == "cell_statistics":
            s = sources.load_source(path)
            cell.append(s)
            weights.append(float(e.info["fusion_weight"]))
            by_name[e.name] = s
        elif e.kind == "cd4_de":
            with np.load(path) as d:
                cd4 = {k: d[k] for k in d.files}
            cd4_w = float(e.info["fusion_weight"])
            by_name[e.name] = cd4
        else:
            raise ValueError(e.kind)
    return cell, weights, cd4, cd4_w, by_name


# ----------------------------------------------------------------------------- emission


def _emit_one(i: int):
    g = _G
    counts = generator.dual_moment_counts(
        g["template"],
        g["p_cpm"][i],
        g["p_bulk"][i],
        depths=g["depths"],
        seed=generator.seed_for(f"{g['ctx']}:{g['targets'][i]}"),
    )
    return sparse.csr_matrix(counts.astype(np.float32))


def context_moments(panel: Panel, ci: int, ctx: str, cell, weights, cd4, cd4_w, pairs):
    """The C1 expected moments and control template for one context (frozen recipe)."""
    import anndata as ad

    ctrl = ad.read_h5ad(panel.control_files[ctx])
    if not np.array_equal(ctrl.var_names.astype(str), panel.genes):
        raise ValueError("control gene order differs from gene_names.csv")
    template, depths, mean, bulk = generator.control_template(
        ctrl.X, panel.cells_per_pert, generator.POOL_K, generator.SEED + ci
    )
    controls = {"log2fc": mean, "bulk_delta": bulk}
    effects = fusion.fused_effects(cell, weights, cd4, cd4_w, panel.targets, panel.genes, controls)
    p_cpm, p_bulk = generator.expected_moments(
        effects, controls, panel.targets, panel.genes, shrink=None, pairs=pairs
    )
    return template, depths, p_cpm, p_bulk, bulk, effects


def emit(panel: Panel, registry: Registry, prepared: dict, out_path: Path, *, jobs: int, log=print):
    """Write the prediction ``.h5ad`` for every context of the panel. Returns provenance."""
    out_path = Path(out_path)
    if out_path.exists():
        raise FileExistsError(out_path)
    cell, weights, cd4, cd4_w, _ = load_c1_inputs(registry, prepared)
    gtf = resolve(registry.raw["promoter_annotation"])
    pairs = generator.promoter_pairs(generator.gencode_tss(gtf), panel.targets, panel.genes)
    pending = out_path.with_name(out_path.stem + ".partial.h5ad")
    pending.unlink(missing_ok=True)
    writer = generator.CountWriter(
        pending, panel.targets, panel.genes, list(panel.contexts), panel.cells_per_pert
    )
    started = time.time()
    moments = {}
    for ci, ctx in enumerate(panel.contexts):
        template, depths, p_cpm, p_bulk, bulk, _ = context_moments(
            panel, ci, ctx, cell, weights, cd4, cd4_w, pairs
        )
        moments[ctx] = (p_cpm, p_bulk, bulk)
        _G.clear()
        _G.update(
            template=template, depths=depths, p_cpm=p_cpm, p_bulk=p_bulk, ctx=ctx,
            targets=panel.targets,
        )  # fmt: skip
        with mp.get_context("fork").Pool(jobs) as pool:
            for k, block in enumerate(pool.imap(_emit_one, range(len(panel.targets)), chunksize=2)):
                writer.append(block)
                if (k + 1) % 50 == 0:
                    log(f"{ctx}: {k + 1}/{len(panel.targets)} ({time.time() - started:.0f}s)")
    writer.close()
    pending.rename(out_path)
    np.savez_compressed(
        out_path.with_name("expected_moments.npz"),
        targets=panel.targets,
        genes=panel.genes,
        **{
            f"{c}_{k}": v
            for c, (a, b, bl) in moments.items()
            for k, v in [("p_cpm", a), ("p_bulk", b), ("ctrl_bulk", bl)]
        },
    )
    pairs.to_csv(out_path.with_name("promoter_pairs.csv"), index=False)
    return {
        "seconds": round(time.time() - started, 1),
        "nnz": int(writer.nnz),
        "promoter_pairs": int(len(pairs)),
        "seed": generator.SEED,
        "template_seeds": {c: generator.SEED + i for i, c in enumerate(panel.contexts)},
        "cell_seed_rule": "seed_for(f'{context}:{target}')",
        "sources": [e.name for e in registry.c1_sources()],
        "weights": {e.name: float(e.info["fusion_weight"]) for e in registry.c1_sources()},
    }


# ----------------------------------------------------------------------------- invariants


def validate(panel: Panel, path: Path) -> dict:
    """Local invariants on the written prediction; all must hold before packaging."""
    import h5py

    from virtual_cell.arc import bundle, generate

    report = bundle.inspect_bundle(path, panel.genes, pert_col="target_gene", context_col="context")
    sizes = report.cells_per_group
    n_t, n_c = len(panel.targets), len(panel.contexts)
    with h5py.File(path, "r") as f:
        data = f["X/data"]
        zeros = sum(
            int((data[i : i + (1 << 26)] == 0).sum()) for i in range(0, data.shape[0], 1 << 26)
        )
        obs = f["obs"]

        def col(name):
            node = obs[name]
            cats = np.array(
                [c.decode() if isinstance(c, bytes) else str(c) for c in node["categories"][()]]
            )
            return cats[node["codes"][()]]

        obs_targets, obs_ctx = col("target_gene"), col("context")
    expected_targets = np.tile(np.repeat(panel.targets, panel.cells_per_pert), n_c)
    expected_ctx = np.repeat(np.asarray(panel.contexts), n_t * panel.cells_per_pert)
    checks = {
        f"n_cells == {panel.n_cells}": report.n_cells == panel.n_cells,
        f"n_genes == {len(panel.genes)}": report.n_genes == len(panel.genes),
        "gene order matches gene_names.csv": report.gene_order_matches,
        "raw integer counts": report.integer_valued,
        "non-negative": report.non_negative,
        "finite": report.finite,
        f"exactly {panel.cells_per_pert} cells per (context, perturbation)": bool(
            (sizes == panel.cells_per_pert).all()
        ),
        f"{n_t} perturbations in every context": all(
            v == n_t for v in report.n_perturbations_per_context.values()
        ),
        f"contexts are exactly {list(panel.contexts)}": tuple(report.contexts)
        == tuple(sorted(panel.contexts))
        or tuple(report.contexts) == tuple(panel.contexts),
        "obs context order is the manifest order": bool(np.array_equal(obs_ctx, expected_ctx)),
        "obs target order is the pert_counts.csv order": bool(
            np.array_equal(obs_targets, expected_targets)
        ),
        "no control cells emitted": panel.control_label not in set(sizes.index.get_level_values(1)),
        "under max_counts_per_cell (1,000,000)": report.max_counts_per_cell
        < generate.MAX_COUNTS_PER_CELL,
        "under max_nnz (4,750,000,000)": report.nnz < 4_750_000_000,
        "target set equals pert_counts.csv": set(sizes.index.get_level_values(1))
        == {str(t) for t in panel.targets},
        "no explicit stored zeros": zeros == 0,
    }
    return {
        "checks": checks,
        "all_pass": all(checks.values()),
        "nnz": int(report.nnz),
        "max_counts_per_cell": float(report.max_counts_per_cell),
        "shape": [int(report.n_cells), int(report.n_genes)],
        "target_sha256_in_prediction": list_sha(panel.targets),
    }


def write_json(path: Path, obj) -> None:
    Path(path).write_text(json.dumps(obj, indent=2, default=str) + "\n")
