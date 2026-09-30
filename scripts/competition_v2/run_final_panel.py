"""Build the frozen C1 prediction for any official panel (final D/E/F, or validation A/B/C).

    uv run python scripts/competition_v2/run_final_panel.py \
        --controls-dir <official_bundle> \
        --source-registry configs/source_registry.yaml \
        --output-dir outputs/final/c1

Steps:

1. inspect the panel (manifest, targets, genes, contexts, cells per perturbation);
2. verify the gene axis of every control file;
3. audit source coverage (C1 sources and audit-only sources);
4. build the C1 mean responses (sources reused or re-prepared for this panel);
5. run the frozen C1 generator;
6. write ``prediction.h5ad``;
7. check the local invariants;
8. write provenance (``provenance.json``, ``panel_audit.md``, ``coverage.csv``);
9. print the packaging instructions.

**It never submits.** It refuses to overwrite an existing prediction.
"""

from __future__ import annotations

# ruff: noqa: E501  (Markdown table rows in the audit writer)
import os

os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import argparse  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from datetime import UTC, datetime  # noqa: E402
from pathlib import Path  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from virtual_cell.competition_v2 import final, panel, vcc_compat  # noqa: E402


def git_state() -> dict:
    def run(*a):
        return subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True).stdout.strip()

    return {"head": run("rev-parse", "HEAD"), "dirty": bool(run("status", "--porcelain"))}


def audit_markdown(
    p: panel.Panel, summary: dict, prepared: dict, emit: dict | None, val: dict | None, prov: dict
) -> str:
    s = p.summary()
    lines = [
        f"# Panel audit: {s['panel_id'] or 'unnamed panel'} ({s['partition'] or 'no partition'})",
        "",
        f"Generated {prov['utc']} by `scripts/competition_v2/run_final_panel.py`. "
        "Nothing was submitted.",
        "",
        "## Panel",
        "",
        "| field | value |",
        "|---|---|",
        f"| contexts | {', '.join(s['contexts'])} |",
        f"| targets | {s['n_targets']} |",
        f"| genes | {s['n_genes']} |",
        f"| cells per perturbation | {s['cells_per_pert']} |",
        f"| prediction cells | {s['n_cells']:,} |",
        f"| target list SHA-256 (pert_counts.csv order) | `{s['target_list_sha256']}` |",
        f"| gene list SHA-256 | `{s['gene_list_sha256']}` |",
        "",
        "Control files:",
        "",
        "| context | file | identified by |",
        "|---|---|---|",
        *[
            f"| {c} | `{Path(f).name}` | {s['control_file_discovery'][c]} |"
            for c, f in s["control_files"].items()
        ],
        "",
        "File checksums:",
        "",
        "| file | SHA-256 |",
        "|---|---|",
        *[f"| `{k}` | `{v}` |" for k, v in s["file_sha256"].items()],
        "",
        "## VCC CLI compatibility",
        "",
        f"- detected version: **{prov['vcc_cli']['version']}** "
        f"({'tested' if prov['vcc_cli']['tested_version'] else 'UNTESTED'}; tested: 0.2.0)",
        f"- executable: `{prov['vcc_cli']['info']['executable']}`",
        f"- required `vcc prep` options and `run_prep` keywords present: "
        f"{prov['vcc_cli']['compatible']}",
        *[f"- WARNING: {w}" for w in prov["vcc_cli"]["warnings"]],
        "",
        "## Source coverage (C1 sources)",
        "",
        f"C1 sources, in fusion order: {', '.join(summary['c1_sources'])}. "
        "A target is covered by a source with ≥ 20 cells (CD4: usable in ≥ 1 QC-passing condition).",
        "",
        "| sources per target | 0 | 1 | 2 | 3+ |",
        "|---|---|---|---|---|",
        "| targets | "
        + " | ".join(str(v) for v in summary["c1_source_count_distribution"].values())
        + " |",
        "",
        f"**Unsupported targets ({len(summary['unsupported_targets'])})**, emitted with the frozen C1 "
        "fallback (context control composition + promoter cap):",
        "",
        ", ".join(summary["unsupported_targets"]) or "none",
        "",
        "Audit-only sources (never fused):",
        "",
        "| source | license | qualification | measure | panel targets covered |",
        "|---|---|---|---|---|",
        *[
            f"| {k} | {v['license_status']} | {v['qualification']} | {v['measure']} | "
            f"{summary['audit_only_coverage'].get(k, 0)} |"
            for k, v in summary["audit_sources"].items()
        ],
        "",
        "## Source versions",
        "",
        "| source | statistics | SHA-256 | action |",
        "|---|---|---|---|",
        *[
            f"| {k} | `{v['path']}` | `{v['sha256'][:16]}…` | {v['action']} |"
            for k, v in prepared.items()
        ],
        "",
    ]
    if emit:
        lines += [
            "## Prediction",
            "",
            "| field | value |",
            "|---|---|",
            f"| file | `{prov['prediction']['path']}` |",
            f"| SHA-256 | `{prov['prediction']['sha256']}` |",
            f"| dimensions | {val['shape'][0]:,} × {val['shape'][1]:,} |",
            f"| nnz | {val['nnz']:,} |",
            "| generator | frozen C1 (4-cell template, dual-moment emitter, amplitude 0.6 / 0.3, clip 3) |",
            f"| generator seed | {emit['seed']}; template seeds {emit['template_seeds']}; cells `{emit['cell_seed_rule']}` |",
            f"| fusion | equal weights {emit['weights']} |",
            f"| promoter pairs | {emit['promoter_pairs']} |",
            f"| emission seconds | {emit['seconds']} |",
            f"| model git | `{prov['git']['head']}` (dirty: {prov['git']['dirty']}) |",
            "",
            "## Raw-count and format checks",
            "",
            *[f"- [{'PASS' if ok else 'FAIL'}] {k}" for k, ok in val["checks"].items()],
            "",
            f"**All checks pass: {val['all_pass']}**",
        ]
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--controls-dir", required=True)
    ap.add_argument("--source-registry", default=str(ROOT / "configs" / "source_registry.yaml"))
    ap.add_argument("--output-dir", required=True)
    ap.add_argument(
        "--cells-per-pert", type=int, default=None, help="only if there is no manifest.json"
    )
    ap.add_argument("--jobs", type=int, default=9)
    ap.add_argument("--audit-only", action="store_true", help="stop after step 3")
    args = ap.parse_args()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    def log(m):
        print(m, flush=True)

    try:  # step 0: fail before any work if packaging could not run later
        cli = vcc_compat.require()
    except vcc_compat.VccCompatibilityError as exc:
        log(f"[0] VCC CLI INCOMPATIBLE\n{exc}")
        return 2
    log(f"[0] vcc CLI {cli.version} ({'tested' if cli.tested_version else 'UNTESTED'}): "
        "required prep/package capabilities present")  # fmt: skip
    for w in cli.warnings:
        log(f"  WARNING: {w}")
    try:
        p = panel.load_panel(args.controls_dir, cells_per_pert=args.cells_per_pert)  # steps 1-2
    except panel.PanelError as exc:
        log(f"[1-2] OFFICIAL BUNDLE NOT UNDERSTOOD\n{exc}")
        return 3
    for c, how in p.discovery.items():
        log(f"  context {c}: {p.control_files[c].name} ({how})")
    log(f"[1-2] panel {p.contexts} | {len(p.targets)} targets | {len(p.genes)} genes | "
        f"{p.cells_per_pert} cells/pert | gene axes verified")  # fmt: skip
    reg = panel.load_registry(args.source_registry)
    prepared = final.prepare_sources(p, reg, out, log=log)  # needed for coverage and step 4
    _, _, _, _, stats = final.load_c1_inputs(reg, prepared)
    table, summary = panel.coverage_audit(p.targets, reg, stats)  # step 3
    table.to_csv(out / "coverage.csv")
    log(f"[3] coverage 0/1/2/3+ = {list(summary['c1_source_count_distribution'].values())}; "
        f"unsupported: {summary['unsupported_targets']}")  # fmt: skip
    prov = {
        "utc": datetime.now(UTC).isoformat(timespec="seconds"),
        "panel": p.summary(),
        "registry": {"path": str(reg.path), "sha256": panel.sha256(reg.path)},
        "sources": prepared,
        "coverage": summary,
        "git": git_state(),
        "vcc_cli": cli.as_dict(),
        "submitted": False,
    }
    emit_info = val = None
    if not args.audit_only:
        pred = out / "prediction.h5ad"
        emit_info = final.emit(p, reg, prepared, pred, jobs=args.jobs, log=log)  # steps 4-6
        log(f"[4-6] wrote {pred} ({emit_info['seconds']}s, nnz {emit_info['nnz']:,})")
        val = final.validate(p, pred)  # step 7
        for k, ok in val["checks"].items():
            log(f"  [{'PASS' if ok else 'FAIL'}] {k}")
        prov.update(
            emit=emit_info,
            validation=val,
            prediction={
                "path": str(pred),
                "sha256": panel.sha256(pred),
                "bytes": pred.stat().st_size,
            },
        )
    prov["elapsed_seconds"] = round(time.time() - t0, 1)
    final.write_json(out / "provenance.json", prov)  # step 8
    (out / "panel_audit.md").write_text(audit_markdown(p, summary, prepared, emit_info, val, prov))
    log(
        f"[8] wrote {out / 'provenance.json'} and {out / 'panel_audit.md'} ({prov['elapsed_seconds']}s)"
    )
    if val is not None and not val["all_pass"]:
        log("LOCAL INVARIANTS FAILED: do not package.")
        return 1
    if val is not None:  # step 9
        log(
            "\n[9] Packaging (not run here; nothing is submitted):\n"
            f"  uv run python scripts/competition_v2/package_final_panel.py "
            f"--output-dir {out} --controls-dir {args.controls_dir}\n"
            "  This compacts the prediction, runs `vcc prep --dry-run`, then builds the .vcc with vcc.prep."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
