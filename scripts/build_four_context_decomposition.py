"""Independent four-context response decomposition, version 1 (canonical run).

Phases 2-7 of the frozen protocol in
``reports/scperteval_four_context_data_spec.md``:

  2. local validation of the downloaded scPertEval files
  3. freeze the balanced design (shared perturbations x shared genes)
  4. canonical pseudobulk and delta construction
  5. canonical decomposition with invariant verification
  6. split-half reliability (primary scheme)
  7. numbers and diagnostic figures

This is an **independent four-context decomposition**, NOT a Molina & Zhang
reproduction; that track is blocked. Their published values are never used as a
target and never tuned toward.

No HVG selection, no scaling, no PCA, no batch correction, no Arc-gene
intersection, no outcome-dependent filtering. ``X`` is already
``log1p(CP10K)`` and is never re-normalised.

Memory: one context's cell matrix is resident at a time; the half-mean tensor is
a disk-backed memmap. The full matrices are never densified.

Usage::

    uv run python scripts/build_four_context_decomposition.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from datetime import date
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import sparse  # noqa: E402

from virtual_cell.data import scperteval  # noqa: E402
from virtual_cell.decomposition import anova  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
DESIGN_VERSION = "four_context_v1"
SCRIPT_VERSION = "1.0.0"

# Expectations from the remote audit. These are checked, never forced.
EXPECTED_N_PERTS = 1_264
EXPECTED_N_GENES = 6_640

CONTEXT_COLORS = {
    "replogle22k562": "#1f77b4",
    "replogle22rpe1": "#d62728",
    "nadig25hepg2": "#2ca02c",
    "nadig25jurkat": "#9467bd",
}


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


# --------------------------------------------------------------------------
# Phase 2
# --------------------------------------------------------------------------


def phase2_validate(data_dir: Path) -> dict[str, dict]:
    print("\n" + "=" * 70)
    print("PHASE 2 - local validation")
    print("=" * 70)
    reports = {}
    for ds in scperteval.DATASETS:
        path = ds.path(data_dir)
        rep = scperteval.structure_report(path)
        rep["cell_line"] = ds.cell_line
        lib = scperteval.library_size_check(path)
        rep["library_size_check"] = lib
        reports[ds.name] = rep
        print(f"\n{ds.name}  ({ds.cell_line})")
        print(f"  shape                 : ({rep['n_cells']:,}, {rep['n_genes']:,})")
        print(f"  X                     : {rep['x_encoding']} / {rep['x_dtype']}")
        print(f"  obs columns           : {rep['obs_columns']}")
        print(f"  var columns           : {rep['var_columns']}")
        print(f"  layers/obsm/uns etc.  : {rep['extras']}")
        print(
            f"  control cells         : {rep['n_control_cells']:,} "
            f"(control label present: {rep['has_control']})"
        )
        print(f"  perturbations         : {rep['n_perturbations']:,}")
        print(
            f"  cells/pert min/med/max: {rep['min_cells_per_pert']} / "
            f"{rep['median_cells_per_pert']:.0f} / {rep['max_cells_per_pert']}"
        )
        print(f"  combination perts     : {rep['n_combination_perturbations']}")
        print(f"  var_names unique      : {rep['var_names_unique']}")
        print(
            f"  expm1(X) row sums     : median={lib['median']:.1f} "
            f"(target {lib['target_sum']:.0f}, ratio {lib['median_over_target']:.4f}), "
            f"p05={lib['p05']:.1f} p95={lib['p95']:.1f}"
        )

    problems = []
    for name, rep in reports.items():
        if rep["x_encoding"] != "csr_matrix":
            problems.append(f"{name}: X is {rep['x_encoding']}, expected csr_matrix")
        if rep["x_dtype"] != "float32":
            problems.append(f"{name}: X dtype {rep['x_dtype']}, expected float32")
        if rep["obs_columns"] != [scperteval.PERTURBATION_KEY]:
            problems.append(f"{name}: unexpected obs columns {rep['obs_columns']}")
        if rep["var_columns"]:
            problems.append(f"{name}: unexpected var columns {rep['var_columns']}")
        for key, val in rep["extras"].items():
            if val:
                problems.append(f"{name}: unexpected {key} entries {val}")
        if not rep["has_control"]:
            problems.append(f"{name}: no control cells")
        if rep["n_combination_perturbations"]:
            problems.append(f"{name}: {rep['n_combination_perturbations']} combination perts")
        if not rep["var_names_unique"]:
            problems.append(f"{name}: duplicate var_names")
    if problems:
        raise SystemExit("PHASE 2 FAILED:\n  " + "\n  ".join(problems))
    print("\nPhase 2: all structural checks passed.")
    return reports


# --------------------------------------------------------------------------
# Phase 3
# --------------------------------------------------------------------------


def phase3_freeze(
    data_dir: Path, splits_dir: Path, provenance: dict
) -> tuple[list[str], list[str]]:
    print("\n" + "=" * 70)
    print("PHASE 3 - freeze the balanced design")
    print("=" * 70)

    label_sets, gene_sets = {}, {}
    for ds in scperteval.DATASETS:
        path = ds.path(data_dir)
        cats, _ = scperteval.read_perturbation_labels(path)
        label_sets[ds.name] = cats
        gene_sets[ds.name] = list(scperteval.read_var_names(path))
        print(f"  {ds.name:16s} labels={len(cats):>5,}  genes={len(gene_sets[ds.name]):>6,}")

    perts = scperteval.shared_perturbations(label_sets)
    genes = scperteval.shared_genes(gene_sets)
    print(f"\n  shared perturbations : {len(perts):,}  (audit expected {EXPECTED_N_PERTS:,})")
    print(f"  shared genes         : {len(genes):,}  (audit expected {EXPECTED_N_GENES:,})")

    mismatches = []
    if len(perts) != EXPECTED_N_PERTS:
        mismatches.append(f"shared perturbations {len(perts)} != expected {EXPECTED_N_PERTS}")
    if len(genes) != EXPECTED_N_GENES:
        mismatches.append(f"shared genes {len(genes)} != expected {EXPECTED_N_GENES}")
    if mismatches:
        raise SystemExit(
            "PHASE 3 STOP - locally derived design disagrees with the audit:\n  "
            + "\n  ".join(mismatches)
        )

    splits_dir.mkdir(parents=True, exist_ok=True)
    contexts_txt = "\n".join(scperteval.CONTEXTS) + "\n"
    perts_txt = "\n".join(perts) + "\n"
    genes_txt = "\n".join(genes) + "\n"
    (splits_dir / "contexts.txt").write_text(contexts_txt)
    (splits_dir / "shared_perturbations.txt").write_text(perts_txt)
    (splits_dir / "shared_genes.txt").write_text(genes_txt)

    manifest = {
        "design_version": DESIGN_VERSION,
        "generated": date.today().isoformat(),
        "generation_script": "scripts/build_four_context_decomposition.py",
        "script_version": SCRIPT_VERSION,
        "counts": {
            "n_contexts": len(scperteval.CONTEXTS),
            "n_shared_perturbations": len(perts),
            "n_shared_genes": len(genes),
        },
        "contexts": list(scperteval.CONTEXTS),
        "cell_lines": {d.name: d.cell_line for d in scperteval.DATASETS},
        "ordering_rules": {
            "contexts": "frozen declaration order in virtual_cell.data.scperteval.CONTEXTS",
            "perturbations": "sorted() ascending over the four-way intersection of obs labels",
            "genes": "sorted() ascending over the four-way intersection of var_names",
        },
        "selection_rule": (
            "Intersections depend ONLY on identifier presence (obs['perturbation'] "
            "categories and var_names). No expression value is consulted. There is "
            "no expression-dependent feature selection: no HVG, no variance filter, "
            "no outcome-dependent filtering."
        ),
        "control_label": scperteval.CONTROL_LABEL,
        "perturbation_key": scperteval.PERTURBATION_KEY,
        "source_dataset_sha256": {r["dataset"]: r["local_sha256"] for r in provenance["datasets"]},
        "scperteval_commit": provenance["scperteval_commit"],
        "file_sha256": {
            "contexts.txt": sha256_text(contexts_txt),
            "shared_perturbations.txt": sha256_text(perts_txt),
            "shared_genes.txt": sha256_text(genes_txt),
        },
    }
    import yaml

    (splits_dir / "design_manifest.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False))
    print(f"\n  frozen design written to {splits_dir}")
    for k, v in manifest["file_sha256"].items():
        print(f"    {k:28s} sha256={v[:16]}...")
    return perts, genes


# --------------------------------------------------------------------------
# Phase 4
# --------------------------------------------------------------------------


def phase4_pseudobulk(
    data_dir: Path, perts: list[str], genes: list[str], processed_dir: Path
) -> tuple[np.ndarray, dict[str, scperteval.ContextPseudobulk]]:
    print("\n" + "=" * 70)
    print("PHASE 4 - canonical pseudobulk")
    print("=" * 70)
    processed_dir.mkdir(parents=True, exist_ok=True)
    pbs, deltas = {}, []
    for ds in scperteval.DATASETS:
        t0 = time.time()
        pb = scperteval.pseudobulk(
            ds.path(data_dir), genes=genes, perturbations=perts, context=ds.name
        )
        pbs[ds.name] = pb
        deltas.append(pb.delta)
        print(
            f"  {ds.name:16s} control_cells={pb.n_control_cells:>7,} "
            f"perts={len(pb.perturbations):,} cells(min/med/max)="
            f"{pb.cell_counts.min()}/{np.median(pb.cell_counts):.0f}/{pb.cell_counts.max()} "
            f"nnz_seen={pb.n_nonzero_seen:,}  [{time.time() - t0:.1f}s]"
        )
    D = np.stack(deltas).astype(np.float64)
    print(f"\n  delta tensor: {D.shape}  (contexts x perturbations x genes)")
    np.save(processed_dir / "delta_tensor.npy", D.astype(np.float32))
    np.save(
        processed_dir / "control_means.npy",
        np.stack([pbs[c].control_mean for c in scperteval.CONTEXTS]).astype(np.float32),
    )
    np.save(
        processed_dir / "cell_counts.npy",
        np.stack([pbs[c].cell_counts for c in scperteval.CONTEXTS]),
    )
    return D, pbs


# --------------------------------------------------------------------------
# Phase 5
# --------------------------------------------------------------------------


def phase5_decompose(D: np.ndarray, perts: list[str]) -> tuple[anova.Decomposition, dict]:
    print("\n" + "=" * 70)
    print("PHASE 5 - canonical decomposition")
    print("=" * 70)
    dec = anova.decompose(D, cell_lines=scperteval.CONTEXTS, perturbations=perts)

    recon = float(np.abs(dec.reconstruct() - D).max())
    zs = {
        "alpha_sum_over_contexts": float(np.abs(dec.alpha.sum(axis=0)).max()),
        "beta_sum_over_perts": float(np.abs(dec.beta.sum(axis=0)).max()),
        "gamma_sum_over_contexts": float(np.abs(dec.gamma.sum(axis=0)).max()),
        "gamma_sum_over_perts": float(np.abs(dec.gamma.sum(axis=1)).max()),
    }
    ss = anova.sums_of_squares(dec)
    total = anova.total_sum_of_squares(D)
    partition_err = abs(sum(ss.values()) - total) / total

    n_c, n_p, _ = D.shape
    parts = {
        "mu": np.broadcast_to(dec.mu, D.shape),
        "alpha": np.broadcast_to(dec.alpha[:, None, :], D.shape),
        "beta": np.broadcast_to(dec.beta[None, :, :], D.shape),
        "gamma": dec.gamma,
    }
    names = list(parts)
    max_cross = 0.0
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            cross = abs(float(np.sum(parts[a] * parts[b])) / (n_c * n_p)) / total
            max_cross = max(max_cross, cross)

    print(f"  exact reconstruction      : max|err| = {recon:.3e}")
    for k, v in zs.items():
        print(f"  zero-sum {k:26s}: max|err| = {v:.3e}")
    print(f"  SS partition              : rel err = {partition_err:.3e}")
    print(f"  max normalised cross term : {max_cross:.3e}")

    tol = 1e-9
    failures = [k for k, v in zs.items() if v > 1e-8 * max(1.0, float(np.abs(D).max()))]
    if recon > 1e-8 * max(1.0, float(np.abs(D).max())):
        failures.append("reconstruction")
    if partition_err > tol:
        failures.append("ss_partition")
    if max_cross > tol:
        failures.append("orthogonality")
    if failures:
        raise SystemExit(f"PHASE 5 FAILED invariants: {failures}")

    fractions = {k: v / total for k, v in ss.items()}
    fractions["template(mu+alpha)"] = fractions["mu"] + fractions["alpha"]
    print("\n  Uncorrected energy shares (denominator = mean_{c,p} ||delta||^2):")
    for k in ("mu", "alpha", "template(mu+alpha)", "beta", "gamma"):
        print(f"    {k:22s} {fractions[k] * 100:6.2f}%")
    print(f"    {'total response energy':22s} {total:.4f}")

    diagnostics = {
        "reconstruction_max_abs_err": recon,
        "zero_sum": zs,
        "ss_partition_rel_err": partition_err,
        "max_normalised_cross_term": max_cross,
        "sums_of_squares": ss,
        "total_ss": total,
        "fractions": fractions,
    }
    return dec, diagnostics


# --------------------------------------------------------------------------
# Phase 6
# --------------------------------------------------------------------------


def phase6_split_half(
    data_dir: Path,
    perts: list[str],
    genes: list[str],
    total_ss: float,
    processed_dir: Path,
    *,
    n_splits: int,
    seed: int,
    keep_memmap: bool,
) -> dict:
    print("\n" + "=" * 70)
    print(f"PHASE 6 - split-half reliability ({n_splits} resamples, seed {seed})")
    print("=" * 70)
    n_c, n_p, n_g = len(scperteval.CONTEXTS), len(perts), len(genes)

    mm_path = processed_dir / "_half_means.npy"
    half = np.lib.format.open_memmap(
        mm_path, mode="w+", dtype=np.float32, shape=(n_splits, 2, n_c, n_p, n_g)
    )
    print(f"  half-mean memmap: {half.nbytes / 1e9:.1f} GB at {mm_path.name}")

    counts_by_context = np.zeros((n_c, n_p), dtype=np.int64)
    for ci, ds in enumerate(scperteval.DATASETS):
        t0 = time.time()
        X, group = scperteval.load_perturbation_cells(
            ds.path(data_dir), genes=genes, perturbations=perts
        )
        order = np.argsort(group, kind="stable")
        X = X[order]
        group = group[order]
        bounds = np.searchsorted(group, np.arange(n_p + 1))
        counts_by_context[ci] = np.diff(bounds)
        ctrl = np.load(processed_dir / "control_means.npy")[ci].astype(np.float64)
        print(
            f"  {ds.name:16s} cells={X.shape[0]:,} nnz={X.nnz:,} "
            f"({X.data.nbytes / 1e9:.2f} GB)  [{time.time() - t0:.1f}s load]",
            flush=True,
        )

        t0 = time.time()
        for s in range(n_splits):
            rng = np.random.default_rng([seed, ci, s])
            rows_a, cols_a, rows_b, cols_b = [], [], [], []
            na = np.zeros(n_p)
            nb = np.zeros(n_p)
            for p in range(n_p):
                lo, hi = bounds[p], bounds[p + 1]
                idx = np.arange(lo, hi)
                rng.shuffle(idx)
                h = len(idx) // 2
                rows_a.append(np.full(h, p))
                cols_a.append(idx[:h])
                rows_b.append(np.full(h, p))
                cols_b.append(idx[h : 2 * h])
                na[p] = nb[p] = h
            for rows, cols, n_half, slot in (
                (rows_a, cols_a, na, 0),
                (rows_b, cols_b, nb, 1),
            ):
                r = np.concatenate(rows)
                c = np.concatenate(cols)
                ind = sparse.csr_matrix(
                    (np.ones(len(r), dtype=np.float32), (r, c)), shape=(n_p, X.shape[0])
                )
                means = np.asarray((ind @ X).todense()) / n_half[:, None]
                half[s, slot, ci] = (means - ctrl).astype(np.float32)
            if (s + 1) % 10 == 0:
                print(f"    split {s + 1}/{n_splits}  [{time.time() - t0:.0f}s]", flush=True)
        del X
        half.flush()

    print("\n  computing cross-half signal per resample ...", flush=True)
    per_split, per_cp_corr = [], np.zeros((n_splits, n_c, n_p))
    for s in range(n_splits):
        D1 = np.asarray(half[s, 0], dtype=np.float64)
        D2 = np.asarray(half[s, 1], dtype=np.float64)
        per_split.append(anova.cross_half_signal(D1, D2))
        a = D1 - D1.mean(axis=2, keepdims=True)
        b = D2 - D2.mean(axis=2, keepdims=True)
        num = (a * b).sum(axis=2)
        den = np.sqrt((a**2).sum(axis=2) * (b**2).sum(axis=2))
        per_cp_corr[s] = np.where(den > 0, num / np.where(den > 0, den, 1.0), 0.0)

    signal = {k: float(np.mean([p[k] for p in per_split])) for k in per_split[0]}
    fractions = anova.noise_corrected_fractions(signal, total_ss)
    per_split_fracs = pd.DataFrame(
        [anova.noise_corrected_fractions(p, total_ss) for p in per_split]
    )

    print("\n  Noise-corrected shares (mean over resamples, +- sd):")
    for k in ("template", "beta", "gamma", "noise"):
        print(
            f"    {k:10s} {fractions[k] * 100:6.2f}%  "
            f"(sd {per_split_fracs[k].std() * 100:.3f}, "
            f"min {per_split_fracs[k].min() * 100:.2f}, max {per_split_fracs[k].max() * 100:.2f})"
        )

    reliability = per_cp_corr.mean(axis=0)
    result = {
        "n_splits": n_splits,
        "seed": seed,
        "signal": signal,
        "fractions": fractions,
        "per_split_fractions": per_split_fracs.to_dict("list"),
        "per_split_signal": [{k: float(v) for k, v in p.items()} for p in per_split],
        "reliability_by_context_pert": reliability,
        "cells_by_context_pert": counts_by_context,
    }
    np.save(processed_dir / "split_half_reliability.npy", reliability)
    if not keep_memmap:
        del half
        mm_path.unlink()
        print(f"\n  removed {mm_path.name} (regenerable)")
    return result


# --------------------------------------------------------------------------
# Phase 7 figures
# --------------------------------------------------------------------------


def phase7_figures(
    dec: anova.Decomposition,
    diag: dict,
    split: dict,
    perts: list[str],
    out_dir: Path,
) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    figs = []
    labels = [scperteval.dataset(c).cell_line for c in scperteval.CONTEXTS]

    # 1. component energy bars
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    f = diag["fractions"]
    keys = ["mu", "alpha", "beta", "gamma"]
    axes[0].bar(
        keys, [f[k] * 100 for k in keys], color=["#3B5BA9", "#5B8FD4", "#6F5CC4", "#D26A3A"]
    )
    axes[0].set_ylabel("% of response energy")
    axes[0].set_title("Uncorrected component energy")
    for i, k in enumerate(keys):
        axes[0].text(i, f[k] * 100, f"{f[k] * 100:.1f}%", ha="center", va="bottom")

    sf = split["fractions"]
    ks = ["template", "beta", "gamma", "noise"]
    axes[1].bar(ks, [sf[k] * 100 for k in ks], color=["#3B5BA9", "#6F5CC4", "#D26A3A", "#C9CED6"])
    sd = pd.DataFrame(split["per_split_fractions"])
    axes[1].errorbar(
        ks,
        [sf[k] * 100 for k in ks],
        yerr=[sd[k].std() * 100 for k in ks],
        fmt="none",
        ecolor="k",
        capsize=4,
    )
    axes[1].set_ylabel("% of response energy")
    axes[1].set_title(f"Noise-corrected ({split['n_splits']} split-half resamples)")
    for i, k in enumerate(ks):
        axes[1].text(i, sf[k] * 100, f"{sf[k] * 100:.1f}%", ha="center", va="bottom")
    fig.suptitle("Independent four-context decomposition (NOT a Molina & Zhang reproduction)")
    fig.tight_layout()
    p = out_dir / "component_energy.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    figs.append(p)

    # 2/3. beta and gamma magnitude distributions
    beta_norm = np.linalg.norm(dec.beta, axis=1)
    gamma_norm = np.linalg.norm(dec.gamma, axis=2)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    axes[0].hist(beta_norm, bins=60, color="#6F5CC4")
    axes[0].set_xlabel(r"$\|\beta_p\|_2$")
    axes[0].set_ylabel("perturbations")
    axes[0].set_title(f"Conserved effect magnitude (n={len(perts):,})")
    for ci, lab in enumerate(labels):
        axes[1].hist(
            gamma_norm[ci],
            bins=60,
            histtype="step",
            lw=1.5,
            label=lab,
            color=CONTEXT_COLORS[scperteval.CONTEXTS[ci]],
        )
    axes[1].set_xlabel(r"$\|\gamma_{c,p}\|_2$")
    axes[1].set_ylabel("perturbations")
    axes[1].set_title("Interaction magnitude by context")
    axes[1].legend()
    fig.tight_layout()
    p = out_dir / "beta_gamma_magnitude.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    figs.append(p)

    # 4/5. reliability distribution and vs cells per perturbation
    rel = split["reliability_by_context_pert"]
    cells = split["cells_by_context_pert"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    for ci, lab in enumerate(labels):
        axes[0].hist(
            rel[ci],
            bins=60,
            histtype="step",
            lw=1.5,
            label=lab,
            color=CONTEXT_COLORS[scperteval.CONTEXTS[ci]],
        )
    axes[0].set_xlabel("split-half Pearson r of $\\delta_{c,p}$")
    axes[0].set_ylabel("perturbations")
    axes[0].set_title("Per-perturbation reliability")
    axes[0].legend()
    for ci, lab in enumerate(labels):
        axes[1].scatter(
            cells[ci],
            rel[ci],
            s=5,
            alpha=0.35,
            label=lab,
            color=CONTEXT_COLORS[scperteval.CONTEXTS[ci]],
            edgecolors="none",
        )
    axes[1].set_xscale("log")
    axes[1].set_xlabel("cells per perturbation")
    axes[1].set_ylabel("split-half Pearson r")
    axes[1].set_title("Reliability vs sampling depth")
    axes[1].legend(markerscale=3)
    fig.tight_layout()
    p = out_dir / "split_half_reliability.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    figs.append(p)

    # 6. context x perturbation heatmap of ||gamma||
    order = np.argsort(gamma_norm.mean(axis=0))
    fig, ax = plt.subplots(figsize=(12, 3.2))
    im = ax.imshow(gamma_norm[:, order], aspect="auto", cmap="magma", interpolation="nearest")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_xlabel(f"perturbations, sorted by mean $\\|\\gamma\\|$ (n={len(perts):,})")
    ax.set_title(r"Interaction magnitude $\|\gamma_{c,p}\|_2$")
    fig.colorbar(im, ax=ax, label=r"$\|\gamma_{c,p}\|_2$")
    fig.tight_layout()
    p = out_dir / "gamma_heatmap.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    figs.append(p)

    # 7. conserved vs context-specific examples
    bf = pd.Series(anova.beta_fraction_per_perturbation(dec))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    axes[0].hist(bf.to_numpy(), bins=60, color="#4C78A8")
    axes[0].set_xlabel(
        r"$\beta$ fraction $=\mathrm{var}(\beta_p)/"
        r"(\mathrm{var}(\beta_p)+\overline{\mathrm{var}(\gamma_{c,p})})$"
    )
    axes[0].set_ylabel("perturbations")
    axes[0].set_title("Transferability of individual perturbations")
    top = bf.sort_values(ascending=False)
    picks = list(top.head(8).items()) + list(top.tail(8).items())
    names = [f"{k}" for k, _ in picks]
    vals = [v for _, v in picks]
    colors = ["#6F5CC4"] * 8 + ["#D26A3A"] * 8
    axes[1].barh(range(len(names)), vals, color=colors)
    axes[1].set_yticks(range(len(names)))
    axes[1].set_yticklabels(names, fontsize=7)
    axes[1].invert_yaxis()
    axes[1].set_xlabel(r"$\beta$ fraction")
    axes[1].set_title("Most conserved (purple) vs most context-specific (orange)")
    fig.tight_layout()
    p = out_dir / "beta_fraction.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    figs.append(p)
    return figs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=REPO_ROOT / scperteval.DATA_SUBDIR)
    parser.add_argument(
        "--splits-dir", type=Path, default=REPO_ROOT / "data" / "splits" / DESIGN_VERSION
    )
    parser.add_argument(
        "--processed-dir", type=Path, default=REPO_ROOT / "data" / "processed" / DESIGN_VERSION
    )
    parser.add_argument("--out-dir", type=Path, default=REPO_ROOT / "outputs" / DESIGN_VERSION)
    parser.add_argument("--n-splits", type=int, default=50)
    parser.add_argument("--seed", type=int, default=anova.REFERENCE_SEED)
    parser.add_argument("--keep-memmap", action="store_true")
    args = parser.parse_args()

    t_start = time.time()
    provenance = json.loads(
        (REPO_ROOT / "data" / "provenance" / "scperteval" / "manifest.json").read_text()
    )

    reports = phase2_validate(args.data_dir)
    perts, genes = phase3_freeze(args.data_dir, args.splits_dir, provenance)
    D, pbs = phase4_pseudobulk(args.data_dir, perts, genes, args.processed_dir)
    dec, diag = phase5_decompose(D, perts)
    split = phase6_split_half(
        args.data_dir,
        perts,
        genes,
        diag["total_ss"],
        args.processed_dir,
        n_splits=args.n_splits,
        seed=args.seed,
        keep_memmap=args.keep_memmap,
    )
    figs = phase7_figures(dec, diag, split, perts, args.out_dir)

    summary = {
        "design_version": DESIGN_VERSION,
        "script_version": SCRIPT_VERSION,
        "generated": date.today().isoformat(),
        "runtime_seconds": time.time() - t_start,
        "validation": reports,
        "design": {
            "contexts": list(scperteval.CONTEXTS),
            "n_shared_perturbations": len(perts),
            "n_shared_genes": len(genes),
        },
        "pseudobulk": {
            c: {
                "n_control_cells": pbs[c].n_control_cells,
                "cells_per_pert_min": int(pbs[c].cell_counts.min()),
                "cells_per_pert_median": float(np.median(pbs[c].cell_counts)),
                "cells_per_pert_max": int(pbs[c].cell_counts.max()),
                "total_perturbed_cells": int(pbs[c].cell_counts.sum()),
            }
            for c in scperteval.CONTEXTS
        },
        "decomposition": {k: v for k, v in diag.items() if k != "fractions"}
        | {"fractions": diag["fractions"]},
        "split_half": {
            k: v
            for k, v in split.items()
            if k not in ("reliability_by_context_pert", "cells_by_context_pert")
        },
        "figures": [str(p.relative_to(REPO_ROOT)) for p in figs],
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))

    # committed manifest for the derived pseudobulk
    derived = {
        "design_version": DESIGN_VERSION,
        "script_version": SCRIPT_VERSION,
        "generated": date.today().isoformat(),
        "generation_parameters": {
            "n_splits": args.n_splits,
            "seed": args.seed,
            "genes": "all shared genes, no HVG selection",
            "normalisation": "none applied; X is already scPertEval log1p(CP10K)",
            "delta": "perturbation_mean - control_mean on the frozen gene set",
        },
        "artifacts": {
            name: hashlib.sha256((args.processed_dir / name).read_bytes()).hexdigest()
            for name in (
                "delta_tensor.npy",
                "control_means.npy",
                "cell_counts.npy",
                "split_half_reliability.npy",
            )
            if (args.processed_dir / name).exists()
        },
        "source_dataset_sha256": {r["dataset"]: r["local_sha256"] for r in provenance["datasets"]},
    }
    import yaml

    (REPO_ROOT / "data" / "provenance" / "scperteval" / "derived_manifest.yaml").write_text(
        yaml.safe_dump(derived, sort_keys=False)
    )

    print(f"\nTotal runtime: {(time.time() - t_start) / 60:.1f} min")
    print(f"Wrote {len(figs)} figures and summary.json to {args.out_dir}")


if __name__ == "__main__":
    main()
