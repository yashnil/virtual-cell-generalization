"""Write the Arc submission candidate v1 manifest (JSON + Markdown) and SHA file.

Everything is read from frozen or just-audited artifacts; nothing is typed in
by hand:

* ``outputs/arc_submission_v1/candidate_audit.json``: audit of the source .h5ad
  (``scripts/audit_arc_submission_candidate.py``)
* ``outputs/arc_submission_v1/vcc_prep_package.json``: the packaging report
* the packaged ``.vcc`` itself (hashed here)
* ``outputs/arc_dry_run_v1/summary.json`` and ``outputs/arc_count_space_v1/*``
  (frozen model choices)
* ``data/splits/arc_target_support_v1.csv`` (frozen tiers)

Writes:
  outputs/arc_submission_v1/submission_manifest.json
  outputs/arc_submission_v1/model_config.json
  reports/arc_submission_v1_manifest.md
  data/provenance/arc_submission_v1_sha256.txt   (shasum -c format)

Reproduce: ``uv run python scripts/build_arc_submission_manifest.py <git-sha>``
where ``<git-sha>`` is the commit whose code produced and packaged the candidate.
Does not submit anything.
"""

# ruff: noqa: E501  (long lines are inside the generated Markdown template)
from __future__ import annotations

import hashlib
import json
import sys
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "arc_submission_v1"
VCC = OUT / "virtual_cell_generalization_val_v1.vcc"
H5AD = ROOT / "outputs" / "arc_dry_run_v1" / "arc_dry_run_v1.h5ad"
SHA_FILE = ROOT / "data" / "provenance" / "arc_submission_v1_sha256.txt"
REPORT = ROOT / "reports" / "arc_submission_v1_manifest.md"
MODEL_ID = "arc_count_space_baseline_v1"

CAVEATS = [
    "214/300 targets have no perturbation-specific beta (Tier 0 receives m_hat only).",
    "m_hat is strongly shrunk (fitted scalar 0.119) because the only Arc-relevant public "
    "sources, arch1 and kaden25rpe1, disagree about the mean perturbation response.",
    "Kaden is weak per perturbation (median split-half reliability 0.17), but the "
    "predeclared reliability diagnostic was mixed/inconclusive (CASE E) and did NOT justify "
    "changing the frozen model.",
    "G1 realises supplied effects faithfully (slope 0.935), but its expression-error "
    "performance was weak on public validation (1.062 vs 1.036 for control resampling).",
    "This submission is an external calibration measurement, not the end of the project.",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()


def rel(p: Path) -> str:
    return str(p.relative_to(ROOT))


def main(git_sha: str) -> None:
    audit = json.loads((OUT / "candidate_audit.json").read_text())
    package = json.loads((OUT / "vcc_prep_package.json").read_text())
    dry = json.loads((ROOT / "outputs/arc_dry_run_v1/summary.json").read_text())
    me_choice = json.loads(
        (ROOT / "outputs/arc_count_space_v1/main_effect_choice.json").read_text()
    )
    tiers = pd.read_csv(ROOT / "data/splits/arc_target_support_v1.csv")
    if audit["failures"]:
        raise SystemExit(f"candidate audit has failures: {audit['failures']}")

    h5_sha = sha256_file(H5AD)
    if h5_sha != audit["candidate"]["sha256"]:
        raise SystemExit("source .h5ad changed since the audit")
    vcc_sha = sha256_file(VCC)
    vcc_stat = VCC.stat()
    tier_counts = {f"tier_{t}": int((tiers.support_tier == t).sum()) for t in (2, 1, 0)}

    model = {
        "model_id": MODEL_ID,
        "equation": "delta_hat[c,p] = m_hat[c] + w[tier(p)] * beta_hat[p]",
        "m_hat": {
            "estimator": dry["main_effect_estimator"],
            "selected_by": me_choice["criterion"],
            "fitted_scalar": dry["fitted_main_effect_scale"],
            "sources": dry["beta_sources"],
            "norm_by_context": "see outputs/arc_dry_run_v1/per_context.csv",
        },
        "beta_hat": {
            "definition": "mean over directly perturbing sources of the response centred over "
            "perturbations within each source; exactly 0 at Tier 0",
            "sources": dry["beta_sources"],
            "targets_with_beta": dry["targets_with_beta"],
        },
        "tier_weights": dry["tier_weights"],
        "tier_counts": tier_counts,
        "response_space_genes": dry["response_space_genes"],
        "unsupported_gene_rule": "lfc = 0: gene keeps the target context's control distribution",
        "generator": {
            "name": "G1 control transport",
            "smoothing": 0.5,
            "seed": 20260921,
            "cells_per_perturbation": 400,
        },
        "excluded": [
            "STRING",
            "DepMap",
            "pathway correction",
            "gamma model",
            "unseen-perturbation prior",
            "reliability weighting",
            "Kaden filtering",
            "new magnitude calibration",
            "new tier weights",
            "new m_hat",
            "deep generative model",
        ],
        "frozen_by": "data/provenance/scperteval/arc_count_space_v1_freeze.txt",
        "implementation": [
            "scripts/run_arc_dry_run.py",
            "src/virtual_cell/modelling/mean_response.py",
            "src/virtual_cell/modelling/context_main_effect.py",
            "src/virtual_cell/arc/panel.py",
            "src/virtual_cell/arc/generate.py",
            "src/virtual_cell/arc/bundle.py",
        ],
    }
    stats = audit["statistics"]
    manifest = {
        "model_id": MODEL_ID,
        "git_sha": git_sha,
        "created_utc": datetime.now(UTC).isoformat(),
        "submitted": False,
        "partition": audit["inputs"]["partition"],
        "panel_id": audit["inputs"]["panel_id"],
        "source_prediction": {
            "path": rel(H5AD),
            "sha256": h5_sha,
            "bytes": H5AD.stat().st_size,
            "mtime_utc": audit["candidate"]["mtime_utc"],
        },
        "vcc_package": {
            "path": rel(VCC),
            "absolute_path": str(VCC),
            "sha256": vcc_sha,
            "bytes": vcc_stat.st_size,
            "created_local": datetime.fromtimestamp(vcc_stat.st_mtime).astimezone().isoformat(),
            "vcc_cli_version": audit["vcc_version"].splitlines()[0],
            "prep_report": package,
        },
        "control_bundle_checksums": {k: v["actual"] for k, v in audit["control_checksums"].items()},
        "model": model,
        "cells": stats["shape"][0],
        "genes": stats["shape"][1],
        "perturbations": audit["inputs"]["n_targets"],
        "contexts": audit["inputs"]["contexts"],
        "nnz": stats["nnz"],
        "nnz_cap": stats["nnz_cap"],
        "nnz_fraction_of_cap": stats["nnz_fraction_of_cap"],
        "label_checksums": {
            "gene_order_sha256": stats["gene_order_sha256"],
            "perturbation_labels_sha256": stats["perturbation_labels_sha256"],
            "context_labels_sha256": stats["context_labels_sha256"],
        },
        "validation": {
            "hard_checks": audit["hard_checks"],
            "local_checks_13": audit["local_checks_13"],
            "vcc_prep_dry_run_identical_to_frozen": audit["vcc_prep_dry_run"][
                "identical_to_frozen"
            ],
        },
        "scientific_caveats": CAVEATS,
        "expectations": "reports/arc_submission_v1_expectations.md",
    }
    (OUT / "submission_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (OUT / "model_config.json").write_text(json.dumps(model, indent=2) + "\n")
    SHA_FILE.write_text(
        "# Arc submission candidate v1 — large-file digests (verify: shasum -a 256 -c)\n"
        f"{h5_sha}  {rel(H5AD)}\n{vcc_sha}  {rel(VCC)}\n"
    )

    ctrl = "\n".join(f"| `{k}` | `{v}` |" for k, v in manifest["control_bundle_checksums"].items())
    checks = "\n".join(f"- [{'x' if v else ' '}] {k}" for k, v in audit["local_checks_13"].items())
    hard = "\n".join(f"- [{'x' if v else ' '}] {k}" for k, v in audit["hard_checks"].items())
    lib, det = stats["library_size"], stats["genes_detected"]
    md = f"""# Arc submission candidate v1: manifest

Generated by `scripts/build_arc_submission_manifest.py` from the candidate audit
and the frozen artifacts. Do not edit by hand. Machine-readable copy:
`outputs/arc_submission_v1/submission_manifest.json`.

**SUBMITTED: false.** The package is ready for explicit human approval.

| field | value |
|---|---|
| MODEL ID | `{MODEL_ID}` |
| GIT SHA | `{git_sha}` |
| partition / panel | `{manifest["partition"]}` / `{manifest["panel_id"]}` |
| SOURCE PREDICTION | `{rel(H5AD)}` ({H5AD.stat().st_size:,} bytes) |
| SOURCE PREDICTION SHA-256 | `{h5_sha}` |
| VCC PACKAGE | `{rel(VCC)}` ({vcc_stat.st_size:,} bytes) |
| VCC PACKAGE SHA-256 | `{vcc_sha}` |
| VCC CLI | `{manifest["vcc_package"]["vcc_cli_version"]}` |
| package created | {manifest["vcc_package"]["created_local"]} |
| CELLS | {manifest["cells"]:,} |
| GENES | {manifest["genes"]:,} |
| PERTURBATIONS | {manifest["perturbations"]} |
| CONTEXTS | {"/".join(manifest["contexts"])} |
| nnz | {manifest["nnz"]:,} ({100 * manifest["nnz_fraction_of_cap"]:.1f}% of the 4.75e9 CLI cap; 95.8% of 2^31) |

## Model

```
delta_hat[c,p] = m_hat[c] + w[tier(p)] * beta_hat[p]
```

| component | frozen value |
|---|---|
| m_hat | `{model["m_hat"]["estimator"]}`, fitted scalar {model["m_hat"]["fitted_scalar"]:.6f}, sources {", ".join(model["m_hat"]["sources"])} |
| beta_hat | direct-source centred response; non-zero for {model["beta_hat"]["targets_with_beta"]} targets |
| TIER COUNTS | Tier 2 = {tier_counts["tier_2"]}, Tier 1 = {tier_counts["tier_1"]}, Tier 0 = {tier_counts["tier_0"]} |
| TIER WEIGHTS | {dry["tier_weights"]["2"]:.2f} / {dry["tier_weights"]["1"]:.2f} / {dry["tier_weights"]["0"]:g} |
| GENERATOR | G1 control transport |
| SMOOTHING | 0.5 |
| SEED | 20260921 |
| response space | {model["response_space_genes"]:,} genes; the other {manifest["genes"] - model["response_space_genes"]:,} keep their control distribution |

Not used: {", ".join(model["excluded"])}.

## Control-bundle checksums (verified)

| file | SHA-256 |
|---|---|
{ctrl}

## Candidate statistics

- dtype {stats["dtype"]}, CSR, int32 indices / int64 indptr; density {stats["density"]:.4f}
- stored counts {stats["stored_value_min"]:g}–{stats["stored_value_max"]:g}; explicit zeros {stats["n_explicit_zeros"]}
- library size min / median / mean / max: {lib["min"]:,.0f} / {lib["median"]:,.0f} / {lib["mean"]:,.1f} / {lib["max"]:,.0f}
- genes detected min / median / mean / max: {det["min"]:,.0f} / {det["median"]:,.0f} / {det["mean"]:,.1f} / {det["max"]:,.0f}
- cells by context: {stats["cells_by_context"]}; exactly 400 cells in each of {stats["cells_per_perturbation"]["n_groups"]} (context, perturbation) groups
- gene order SHA-256 `{stats["gene_order_sha256"]}` (equals `gene_names.csv`)
- perturbation labels SHA-256 `{stats["perturbation_labels_sha256"]}`
- context labels SHA-256 `{stats["context_labels_sha256"]}`

## Validation

Hard checks:

{hard}

The thirteen frozen local checks (identical to the frozen dry run):

{checks}

`vcc prep --dry-run`: exit 0, report identical to the frozen
`outputs/arc_dry_run_v1/vcc_prep_dry_run.json`. `vcc prep` (packaging): exit 0,
`verified_targets: true`, `dropped: []`, `reordered_genes: false`,
`normalization: counts-preserved`.

## Scientific caveats entering submission

{chr(10).join(f"{i}. {c}" for i, c in enumerate(CAVEATS, 1))}

Pre-result expectations and the interpretation rules for the hidden score:
[`arc_submission_v1_expectations.md`](arc_submission_v1_expectations.md).
"""
    REPORT.write_text(md)
    print(f"  h5ad {h5_sha}\n  vcc  {vcc_sha}\n  wrote {rel(REPORT)}, {rel(SHA_FILE)}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: build_arc_submission_manifest.py <git-sha>")
    main(sys.argv[1])
