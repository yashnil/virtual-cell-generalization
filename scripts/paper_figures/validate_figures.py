# ruff: noqa: E501 (figure labels and provenance notes are kept as single strings)
"""Validate the frozen upstream data, the figure-source tables and the figure outputs.

1. Re-hash every entry of every frozen manifest the paper-figure pipeline reads from, and confirm it matches.
2. For every source table: its CSV hash and every frozen input hash in its provenance sidecar match the current files.
3. For every figure manifest: the SVG/PNG/PDF and every source table hash match the current files.
4. Every figure has a caption file, and every source table named by a manifest exists.
5. No sidecar carries volatile fields (wall-clock date, git HEAD/dirty state) that would make a rebuild non-idempotent.

Writes ``reports/paper_figures/validation.json``; exits non-zero on any failure.
"""

from __future__ import annotations

import json
import sys

from _common import OUT_DIR, ROOT, SRC_DIR, sha256

FROZEN_MANIFESTS = [
    "outputs/n1_n4/manifest_sha256.txt",
    "outputs/n3/manifest_sha256.txt",
    "outputs/n5/manifest_sha256.txt",
    "outputs/n6/manifest_sha256.txt",
    "data/provenance/scperteval/decomposition_phase_freeze.txt",
    "data/provenance/scperteval/canonical_v1_freeze.txt",
    "data/provenance/research_v3/n5_freeze_sha256.txt",
]

VOLATILE_KEYS = {"generated", "git_head"}

FIGURES = [
    ("fig1_transfer_components", "fig1_caption.md"),
    ("fig2_sample_complexity", "fig2_caption.md"),
    ("fig3_source_count", "fig3_caption.md"),
    ("fig4_source_compatibility", "fig4_caption.md"),
    ("ext_fig1_robustness", "ext_fig1_caption.md"),
    ("ext_fig2_n1_controls", "ext_fig2_caption.md"),
    ("ext_fig3_agreement_correction", "ext_fig3_caption.md"),
    ("ext_fig4_n6_reliability", "ext_fig4_caption.md"),
    ("ext_fig5_n3_reference_frame", "ext_fig5_caption.md"),
]


def check_frozen() -> dict:
    out = {}
    for rel in FROZEN_MANIFESTS:
        man = ROOT / rel
        base = man.parent if rel.startswith("outputs/") else ROOT
        n, bad, missing = 0, [], []
        for line in man.read_text().splitlines():
            parts = line.split()
            if len(parts) != 2 or len(parts[0]) != 64:
                continue
            p = base / parts[1].removeprefix("./")
            n += 1
            if not p.exists():
                missing.append(parts[1])
            elif sha256(p) != parts[0]:
                bad.append(parts[1])
        out[rel] = {"entries": n, "mismatch": bad, "missing": missing}
    return out


def check_sources() -> dict:
    out = {}
    for side in sorted(SRC_DIR.glob("*.provenance.json")):
        rec = json.loads(side.read_text())
        csv = ROOT / rec["figure_source"]
        errs = [f"volatile field {k}" for k in sorted(VOLATILE_KEYS & rec.keys())]
        if sha256(csv) != rec["sha256"]:
            errs.append("table hash")
        for inp in rec["frozen_inputs"]:
            if sha256(ROOT / inp["path"]) != inp["sha256"]:
                errs.append(f"input {inp['path']}")
        out[csv.name] = errs
    return out


def check_figures() -> dict:
    out = {}
    for name, cap in FIGURES:
        man = json.loads((OUT_DIR / f"{name}.manifest.json").read_text())
        errs = [f"volatile field {k}" for k in sorted(VOLATILE_KEYS & man.keys())]
        for ext, rec in man["outputs"].items():
            if sha256(ROOT / rec["path"]) != rec["sha256"]:
                errs.append(f"{ext} hash")
        for src in man["figure_sources"]:
            p = ROOT / src["path"]
            if not p.exists():
                errs.append(f"missing {src['path']}")
            elif sha256(p) != src["sha256"]:
                errs.append(f"source {src['path']}")
        if not (OUT_DIR / cap).exists():
            errs.append(f"missing caption {cap}")
        out[name] = errs
    return out


def main() -> None:
    res = {"frozen": check_frozen(), "sources": check_sources(), "figures": check_figures()}
    fail = 0
    for rel, v in res["frozen"].items():
        ok = not v["mismatch"] and not v["missing"]
        fail += not ok
        print(
            f"{'OK ' if ok else 'BAD'} frozen {rel}: {v['entries']} entries, {len(v['mismatch'])} mismatch, {len(v['missing'])} missing"
        )
    for k, v in res["sources"].items():
        fail += bool(v)
        if v:
            print("BAD source", k, v)
    print(
        f"{'OK ' if not any(res['sources'].values()) else 'BAD'} {len(res['sources'])} source tables (table + frozen-input hashes, no volatile fields)"
    )
    for k, v in res["figures"].items():
        fail += bool(v)
        if v:
            print("BAD figure", k, v)
    print(
        f"{'OK ' if not any(res['figures'].values()) else 'BAD'} {len(res['figures'])} figures (SVG/PNG/PDF + source hashes + caption, no volatile fields)"
    )
    (OUT_DIR / "validation.json").write_text(json.dumps(res, indent=2) + "\n")
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()
