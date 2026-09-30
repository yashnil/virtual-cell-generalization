"""Panel-agnostic inputs for the final round: the official bundle, the source registry,
and the target x source coverage audit.

Nothing here assumes the validation panel. Contexts, targets (in file order), genes,
cells per perturbation and control files are read from the official files
(``manifest.json``, ``pert_counts.csv``, ``gene_names.csv``, ``context_<label>.h5ad``),
and every file is checksummed. The source registry (``configs/source_registry.yaml``)
is cross-checked against :mod:`licensing`, which stays the license authority: a
non-GREEN source can never be selected for C1.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from virtual_cell.competition_v2 import licensing
from virtual_cell.competition_v2.sources import MINIMUM_CELLS

ROOT = Path(__file__).resolve().parents[3]
QUALIFIED = "QUALIFIED"


class PanelError(ValueError):
    """The official bundle is incomplete or inconsistent."""


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(16 << 20), b""):
            h.update(block)
    return h.hexdigest()


def list_sha(values) -> str:
    """SHA-256 of a list of strings joined by newlines (order-sensitive)."""
    return hashlib.sha256("\n".join(map(str, values)).encode()).hexdigest()


# ----------------------------------------------------------------------------- panel


@dataclass
class Panel:
    controls_dir: Path
    contexts: tuple[str, ...]
    control_files: dict[str, Path]
    targets: np.ndarray
    genes: np.ndarray
    cells_per_pert: int
    control_label: str
    pert_col: str
    context_col: str
    manifest: dict
    checksums: dict[str, str] = field(default_factory=dict)
    discovery: dict[str, str] = field(default_factory=dict)

    @property
    def n_cells(self) -> int:
        return len(self.contexts) * len(self.targets) * self.cells_per_pert

    def summary(self) -> dict:
        return {
            "controls_dir": str(self.controls_dir),
            "partition": self.manifest.get("partition"),
            "panel_id": self.manifest.get("panel_id"),
            "contexts": list(self.contexts),
            "n_targets": int(len(self.targets)),
            "n_genes": int(len(self.genes)),
            "cells_per_pert": int(self.cells_per_pert),
            "n_cells": int(self.n_cells),
            "target_list_sha256": list_sha(self.targets),
            "gene_list_sha256": list_sha(self.genes),
            "file_sha256": dict(self.checksums),
            "control_files": {c: str(p) for c, p in self.control_files.items()},
            "control_file_discovery": dict(self.discovery),
        }


def _read_var_names(path: Path) -> np.ndarray:
    import h5py

    with h5py.File(path, "r") as f:
        var = f["var"]
        node = var[var.attrs.get("_index", "_index")]
        values = node["values"] if isinstance(node, h5py.Group) else node
        return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in values[()]])


def _read_obs_column(path: Path, col: str) -> np.ndarray:
    import h5py

    with h5py.File(path, "r") as f:
        node = f["obs"][col]
        if isinstance(node, h5py.Group):
            cats = np.array(
                [c.decode() if isinstance(c, bytes) else str(c) for c in node["categories"][()]]
            )
            return cats[node["codes"][()]]
        return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in node[()]])


MANIFEST_SCHEMA_HELP = (
    "expected manifest.json: 'contexts' = a non-empty list of unique labels, or of objects "
    "with a label ('label' | 'context' | 'name') and optionally a file ('file' | 'filename' "
    "| 'path' | 'h5ad'); 'cells_per_pert' = a positive integer; optional per-context files "
    "under 'files' {label: filename} or 'per_context' {label: {file: ...}}; optional "
    "strings 'pert_col', 'context_col', 'control_label'; optional integers 'n_genes', "
    "'n_constructs'"
)
_LABEL_KEYS = ("label", "context", "name")
_FILE_KEYS = ("file", "filename", "path", "h5ad")


def _is_int(v) -> bool:
    return isinstance(v, int) and not isinstance(v, bool)


def parse_manifest(manifest: dict) -> dict:
    """Validate the manifest fields the pipeline needs; raise a diagnostic if not understood.

    Returns ``{"contexts", "cells_per_pert", "files", "pert_col", "context_col",
    "control_label"}``, where ``files`` maps each context to a manifest-named file or None.
    Nothing is guessed: an unknown or malformed schema is an error.
    """
    problems: list[str] = []
    if not isinstance(manifest, dict):
        raise PanelError(f"manifest.json is not a JSON object. {MANIFEST_SCHEMA_HELP}")
    raw = manifest.get("contexts")
    contexts: list[str] = []
    files: dict[str, str | None] = {}
    if isinstance(raw, dict):
        raw = [{"label": k, **(v if isinstance(v, dict) else {"file": v})} for k, v in raw.items()]
    if not isinstance(raw, list) or not raw:
        problems.append("'contexts' is missing or not a non-empty list")
    else:
        for item in raw:
            if isinstance(item, str) and item:
                contexts.append(item)
                files[item] = None
            elif isinstance(item, dict):
                labels = [item[k] for k in _LABEL_KEYS if isinstance(item.get(k), str)]
                names = [item[k] for k in _FILE_KEYS if isinstance(item.get(k), str)]
                if len(set(labels)) != 1:
                    problems.append(f"context entry {item!r} has no single label")
                    continue
                if len(set(names)) > 1:
                    problems.append(f"context entry {item!r} names several files")
                    continue
                contexts.append(labels[0])
                files[labels[0]] = names[0] if names else None
            else:
                problems.append(f"context entry {item!r} is neither a label nor an object")
        if len(set(contexts)) != len(contexts):
            problems.append(f"duplicate context labels {contexts}")
    for label, name in (
        (manifest.get("files") or {}).items() if isinstance(manifest.get("files"), dict) else []
    ):
        if label in files and isinstance(name, str):
            if files[label] not in (None, name):
                problems.append(f"context {label!r} names two different files")
            files[label] = name
    per = manifest.get("per_context")
    if isinstance(per, dict):
        for label, meta in per.items():
            if label not in files and contexts:
                problems.append(f"'per_context' has unknown context {label!r}")
                continue
            if isinstance(meta, dict):
                names = [meta[k] for k in _FILE_KEYS if isinstance(meta.get(k), str)]
                if len(set(names)) > 1:
                    problems.append(f"'per_context'[{label!r}] names several files")
                elif names and label in files:
                    if files[label] not in (None, names[0]):
                        problems.append(f"context {label!r} names two different files")
                    files[label] = names[0]
    cells = manifest.get("cells_per_pert")
    if not _is_int(cells) or cells <= 0:
        problems.append(f"'cells_per_pert' must be a positive integer, got {cells!r}")
    for key in ("pert_col", "context_col", "control_label"):
        if key in manifest and not (isinstance(manifest[key], str) and manifest[key]):
            problems.append(f"{key!r} must be a non-empty string")
    for key in ("n_genes", "n_constructs"):
        if key in manifest and not _is_int(manifest[key]):
            problems.append(f"{key!r} must be an integer")
    if problems:
        raise PanelError(
            "manifest.json schema not understood:\n  - "
            + "\n  - ".join(problems)
            + f"\nfound keys: {sorted(manifest)}\n{MANIFEST_SCHEMA_HELP}"
        )
    return {
        "contexts": tuple(contexts),
        "cells_per_pert": int(cells),
        "files": files,
        "pert_col": manifest.get("pert_col", "target_gene"),
        "context_col": manifest.get("context_col", "context"),
        "control_label": manifest.get("control_label", "non-targeting"),
    }


def discover_control_files(
    d: Path, contexts, named: dict, context_col: str
) -> tuple[dict[str, Path], dict[str, str]]:
    """Resolve each context's control ``.h5ad`` without guessing.

    The order of precedence is: (1) the file the manifest names; (2) the validation-era
    ``context_<label>.h5ad``; (3) the *unique* ``.h5ad`` in the bundle whose obs
    ``context_col`` is uniformly the label. Zero or several candidates is an error.
    """
    d = Path(d).resolve()
    files, how, errors = {}, {}, []
    for c in contexts:
        name = named.get(c)
        if name:
            p = (d / name).resolve()
            if d not in p.parents:
                errors.append(f"{c}: manifest file {name!r} is outside the bundle directory")
            elif not p.exists() or p.suffix != ".h5ad":
                errors.append(f"{c}: manifest file {name!r} is missing or not an .h5ad")
            else:
                files[c], how[c] = p, "named by manifest"
        elif (d / f"context_{c}.h5ad").exists():
            files[c], how[c] = d / f"context_{c}.h5ad", "validation-era filename"
    unresolved = [c for c in contexts if c not in files]
    if unresolved and not errors:
        taken = set(files.values())
        labels = {}
        for p in sorted(d.glob("*.h5ad")):
            if p.resolve() in {f.resolve() for f in taken}:
                continue
            try:
                vals = set(_read_obs_column(p, context_col))
            except Exception as exc:  # noqa: BLE001 - reported, never guessed around
                labels[p] = f"unreadable ({type(exc).__name__})"
                continue
            labels[p] = next(iter(vals)) if len(vals) == 1 else f"mixed {sorted(vals)[:5]}"
        for c in unresolved:
            hits = [p for p, lab in labels.items() if lab == c]
            if len(hits) == 1:
                files[c], how[c] = hits[0], f"unique file with obs {context_col!r} == {c!r}"
            else:
                seen = {p.name: lab for p, lab in labels.items()}
                errors.append(
                    f"{c}: {'no' if not hits else len(hits)} candidate .h5ad files whose obs "
                    f"{context_col!r} is uniformly {c!r} (not guessing); "
                    f"unassigned files and their context labels: {seen}"
                )
    if len(set(files.values())) != len(files):
        errors.append("two contexts resolve to the same file")
    if errors:
        raise PanelError("cannot identify control files:\n  - " + "\n  - ".join(errors))
    return files, how


def load_panel(
    controls_dir: str | Path,
    *,
    cells_per_pert: int | None = None,
    checksums: bool = True,
    verify_controls: bool = True,
) -> Panel:
    """Read and verify an official Arc control bundle; the files are the authority.

    Without ``manifest.json``, contexts are the sorted ``context_<label>.h5ad`` files and
    ``cells_per_pert`` must be given. A ``cells_per_pert`` argument that disagrees with the
    manifest is an error, not an override.
    """
    d = Path(controls_dir)
    manifest_path = d / "manifest.json"
    manifest = {}
    if manifest_path.exists():
        try:
            manifest = json.loads(manifest_path.read_text())
        except json.JSONDecodeError as exc:
            raise PanelError(f"manifest.json is not valid JSON: {exc}") from exc
    if manifest_path.exists():
        spec = parse_manifest(manifest)
        contexts = spec["contexts"]
        mcells = spec["cells_per_pert"]
        if cells_per_pert is not None and int(cells_per_pert) != mcells:
            raise PanelError(f"cells_per_pert {cells_per_pert} disagrees with manifest {mcells}")
        cells_per_pert = mcells
        named = spec["files"]
        pert_col, context_col = spec["pert_col"], spec["context_col"]
        control_label = spec["control_label"]
    else:
        contexts = tuple(
            sorted(
                m.group(1)
                for p in d.glob("context_*.h5ad")
                if (m := re.fullmatch(r"context_(.+)\.h5ad", p.name))
            )
        )
        if cells_per_pert is None:
            raise PanelError("no manifest.json: cells_per_pert must be supplied")
        named = {}
        pert_col, context_col, control_label = "target_gene", "context", "non-targeting"
    if not contexts:
        raise PanelError("no contexts found")
    for f in ("gene_names.csv", "pert_counts.csv"):
        if not (d / f).exists():
            raise PanelError(f"missing bundle files: {[str(d / f)]}")
    files, discovery = discover_control_files(d, contexts, named, context_col)
    needed = [d / "gene_names.csv", d / "pert_counts.csv", *files.values()]
    missing = [str(p) for p in needed if not p.exists()]
    if missing:
        raise PanelError(f"missing bundle files: {missing}")

    genes = pd.read_csv(d / "gene_names.csv").iloc[:, 0].astype(str).to_numpy()
    pc = pd.read_csv(d / "pert_counts.csv")
    if pert_col not in pc.columns:
        raise PanelError(f"pert_counts.csv lacks {pert_col!r}")
    targets = pc[pert_col].astype(str).to_numpy()
    if len(set(targets)) != len(targets):
        raise PanelError("duplicate targets in pert_counts.csv")
    if control_label in set(targets):
        raise PanelError("the control label is listed as a target")
    if len(set(genes)) != len(genes):
        raise PanelError("duplicate genes in gene_names.csv")
    count_cols = [c for c in pc.columns if c != pert_col and pd.api.types.is_numeric_dtype(pc[c])]
    if count_cols:
        per = pc[count_cols[0]].to_numpy()
        if not (per == cells_per_pert).all():
            raise PanelError(
                f"pert_counts.csv column {count_cols[0]!r} is not uniformly {cells_per_pert}; "
                "the C1 writer needs one cells-per-perturbation value"
            )
    if manifest:
        if "n_genes" in manifest and int(manifest["n_genes"]) != len(genes):
            raise PanelError("gene_names.csv length differs from manifest n_genes")
        if "n_constructs" in manifest and int(manifest["n_constructs"]) != len(targets):
            raise PanelError("pert_counts.csv length differs from manifest n_constructs")
    if verify_controls:
        for c, p in files.items():
            if not np.array_equal(_read_var_names(p), genes):
                raise PanelError(
                    f"{p.name}: gene axis differs from gene_names.csv (order included)"
                )
            ctx = _read_obs_column(p, context_col)
            if not (ctx == c).all():
                raise PanelError(f"{p.name}: obs {context_col!r} is not uniformly {c!r}")
            pert = _read_obs_column(p, pert_col)
            if not (pert == control_label).all():
                raise PanelError(f"{p.name}: contains non-control cells")
    sums = {}
    if checksums:
        for p in needed + ([manifest_path] if manifest else []):
            sums[p.name] = sha256(p)
    return Panel(
        d, contexts, files, targets, genes, int(cells_per_pert), control_label, pert_col,
        context_col, manifest, sums, discovery,
    )  # fmt: skip


# ----------------------------------------------------------------------------- registry


@dataclass
class SourceEntry:
    name: str
    licensing_key: str
    license_status: str
    qualification: str
    role: str
    kind: str
    info: dict

    @property
    def c1_usable(self) -> bool:
        return (
            self.role == "c1_source"
            and self.qualification == QUALIFIED
            and licensing.STATUS.get(self.licensing_key, licensing.UNKNOWN) == licensing.GREEN
        )


@dataclass
class Registry:
    path: Path
    entries: list[SourceEntry]
    raw: dict

    def c1_sources(self) -> list[SourceEntry]:
        out = sorted(
            (e for e in self.entries if e.role == "c1_source"), key=lambda e: e.info["c1_order"]
        )
        for e in out:
            if not e.c1_usable:
                raise licensing.LicenseError(
                    f"{e.name} is registered as a C1 source but is not GREEN + QUALIFIED "
                    f"(licensing: {licensing.STATUS.get(e.licensing_key, licensing.UNKNOWN)}, "
                    f"qualification: {e.qualification})"
                )
        licensing.assert_sources_allowed([e.licensing_key for e in out])
        return out

    def audit_sources(self) -> list[SourceEntry]:
        return [e for e in self.entries if e.role != "c1_source"]


def load_registry(path: str | Path) -> Registry:
    """Load the registry and cross-check every license against :mod:`licensing`."""
    path = Path(path)
    raw = yaml.safe_load(path.read_text())
    entries = []
    for s in raw["sources"]:
        code = licensing.STATUS.get(s["licensing_key"], licensing.UNKNOWN)
        if s["license_status"] != code:
            raise licensing.LicenseError(
                f"registry says {s['name']} is {s['license_status']}, licensing.STATUS says {code}"
            )
        entries.append(
            SourceEntry(
                s["name"],
                s["licensing_key"],
                s["license_status"],
                s["qualification"],
                s["role"],
                s["kind"],
                s,
            )  # fmt: skip
        )
    names = [e.name for e in entries]
    if len(set(names)) != len(names):
        raise ValueError("duplicate source names in the registry")
    return Registry(path, entries, raw)


def resolve(rel: str | None) -> Path | None:
    return None if rel is None else (Path(rel) if Path(rel).is_absolute() else ROOT / rel)


# ----------------------------------------------------------------------------- coverage


def c1_usable_targets(entry: SourceEntry, stats) -> set[str]:
    """Targets a C1 source can use, by the frozen C1 rule (>= 20 cells / CD4 QC)."""
    if entry.kind == "cell_statistics":
        return {str(t) for t in stats.usable(MINIMUM_CELLS)}
    if entry.kind == "cd4_de":
        usable = stats["available"] & stats["quality_pass"] & (stats["n_cells"] >= MINIMUM_CELLS)
        return {str(t) for t in np.asarray(stats["targets"]).astype(str)[usable.any(axis=0)]}
    raise ValueError(entry.kind)


def audit_cells(entry: SourceEntry) -> dict[str, float] | None:
    """Per-target measured cells for an audit-only source (``None`` if not auditable)."""
    p = resolve(entry.info.get("coverage_path"))
    if p is None or not p.exists():
        return None
    if entry.kind == "coverage_obs_json":
        gt = json.loads(p.read_text())["obs_values"]["gene_target"]
        return dict(zip(gt["categories"], map(float, gt["counts"]), strict=True))
    if entry.kind == "coverage_manifest":
        m = pd.read_csv(p, sep=r"\s+")
        return dict(zip(m.gene.astype(str), m.cell_count.astype(float), strict=True))
    if entry.kind == "coverage_counts_csv":
        m = pd.read_csv(p)
        return dict(zip(m.gene_target.astype(str), m.n_cells.astype(float), strict=True))
    if entry.kind == "coverage_guide_library":
        names = [ln[1:].strip() for ln in gzip.open(p, "rt") if ln.startswith(">")]
        return {re.sub(r"_\d+$", "", n): np.nan for n in names}  # library presence, no cells
    return None


def coverage_audit(targets, registry: Registry, c1_stats: dict) -> tuple[pd.DataFrame, dict]:
    """Target x source availability for C1 sources and audit sources, plus summary counts.

    Only C1 sources count toward ``n_c1_sources``. Audit sources are reported with their
    license / qualification status and are never fused.
    """
    targets = [str(t) for t in targets]
    table = pd.DataFrame(index=pd.Index(targets, name="target"))
    c1 = registry.c1_sources()
    for e in c1:
        usable = c1_usable_targets(e, c1_stats[e.name])
        table[f"c1:{e.name}"] = [t in usable for t in targets]
    table["n_c1_sources"] = table[[f"c1:{e.name}" for e in c1]].sum(axis=1).astype(int)
    audit_meta = {}
    for e in registry.audit_sources():
        cells = audit_cells(e)
        if cells is None:
            continue
        col = f"audit:{e.name}"
        if e.kind == "coverage_guide_library":
            table[col] = [t in cells for t in targets]
        else:
            table[col] = [cells.get(t, 0.0) >= MINIMUM_CELLS for t in targets]
        audit_meta[e.name] = {
            "license_status": e.license_status,
            "qualification": e.qualification,
            "measure": "guide library presence"
            if e.kind == "coverage_guide_library"
            else ">= 20 cells",
        }
    counts = table.n_c1_sources.clip(upper=3).value_counts().reindex([0, 1, 2, 3], fill_value=0)
    summary = {
        "n_targets": len(targets),
        "c1_sources": [e.name for e in c1],
        "c1_source_count_distribution": {
            "0": int(counts[0]),
            "1": int(counts[1]),
            "2": int(counts[2]),
            "3+": int(counts[3]),
        },
        "unsupported_targets": [t for t in targets if table.loc[t, "n_c1_sources"] == 0],
        "audit_sources": audit_meta,
        "audit_only_coverage": {
            c.split(":", 1)[1]: int(table[c].sum()) for c in table.columns if c.startswith("audit:")
        },
    }
    return table, summary
