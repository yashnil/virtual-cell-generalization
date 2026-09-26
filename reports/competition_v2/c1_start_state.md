# C1 license-clean phase: start state

Recorded 2026-09-26, before any C1 code was written or run. Competition track.
**Nothing is submitted in this phase.**

## Git

| field | value |
|---|---|
| branch | `competition/c1-license-clean` |
| HEAD | `d375f935b6b5baf18539cb363d62c3d2dba86284` ("Complete competitive baseline expansion and AtlasShift audit", 2026-09-26 13:10 −0700) |
| working tree | clean (`git status --porcelain` empty) |

## Provenance verification at start

Script: `scripts/competition_v2/verify_c1_state.py --full`. Machine record:
`outputs/competition_v2/c1_license_clean/state_full_start.json` (2026-09-26T13:37 −0700).
**0 failures.**

| check | items | result |
|---|---|---|
| frozen research manifests (`data/provenance/scperteval/*_freeze.txt`) | 16 manifests, 285 digests | all OK |
| V1 submission digests (`arc_submission_v1_sha256.txt`: dry-run h5ad, submitted `.vcc`) | 2 | OK |
| raw-data checksum files (arc2026 controls, DepMap, Feng, MSigDB, scPertEval ×2, STRING) | 7 files, 20 digests | all OK |
| frozen V1 tier file `data/splits/arc_target_support_v1.csv` | SHA-256 `2280f40b…3ed6f4` | unchanged |
| vendored AtlasShift files vs `c0_manifest.json` | 9 files | unchanged |
| C0 prepared statistics (5 npz) + `official_pairs.csv` | 6 | unchanged |
| C0 `prediction.h5ad` / `prediction_compact.h5ad` / `.vcc` | 3 | unchanged (`155e2440…`, `d47a1518…`, `bd33b1ad…`) |
| competition-v2 URL downloads vs `download_log.jsonl` (K562 GWPS, CD4, GENCODE v47, H1 train/val/test, H1 targets/genes) | 8 files | all SHA-256 match |
| X-Atlas/Orion local copy vs Hugging Face LFS etags (revision `598aa544…`) | 8,689 files, 133.5 GB | all match |

## Competition-v2 source datasets present (`data/raw/competition_v2/`, 233 GiB)

| file | bytes | SHA-256 |
|---|---|---|
| `K562_gwps_raw_singlecell_01.h5ad` | 65,830,941,948 | `b697ef7f…a9de2c` |
| `GWCD4i.DE_stats.h5ad` | 16,786,240,107 | `c355f535…cfbb62` |
| `adata_Training.h5ad` (H1 2025) | 15,482,497,461 | `a0997710…19c02b` |
| `adata_Validation.h5ad` (H1 2025) | 6,928,967,541 | `376f0bab…77bfe1` |
| `adata_Test.h5ad` (H1 2025) | 11,950,739,168 | `ba5ce667…831f03` |
| `pert_counts_Training.csv` / `gene_names.csv` (H1 2025) | 2,824 / 116,023 | `633d202b…` / `e29ff351…` |
| `gencode.v47.annotation.gtf.gz` | 58,953,903 | `df11938c…a3db24` |
| `xatlas_orion/data/HCT116` | 50,544,250,046 | per-file LFS SHA-256 (verified) |
| `xatlas_orion/data/HEK293T` | 82,965,837,799 | per-file LFS SHA-256 (verified) |

Full 64-character hashes: `data/provenance/competition_v2/download_log.jsonl`.

## License evidence entering the phase

Taken from the prior phase (`reports/competition_v2/atlas_shift_source_audit.md` §5):

* K562 GWPS is CC BY 4.0.
* X-Atlas/Orion is CC BY-NC-SA 4.0, with the local `LICENSE.md` full text.
* H1 2025 and CD4: no license file in their buckets; unverified.

This phase re-audits every source in `data_license_register.md`.

Entering policy:

* X-Atlas HCT116 and HEK293T: BLOCKED_PENDING_PERMISSION.
* H1: GREEN_FOR_VCC_USE unless contradicted.
* CD4: GREEN if the official source shows MIT.
* K562: audit before use.

## Artifacts that must not change in this phase

* frozen research artifacts (the 16 manifests above);
* the submitted V1 model and bundle;
* the C0 reproduction outputs (`outputs/competition_v2/atlasshift_c0/`);
* the original tier file.

The end-of-phase re-verification is recorded in `license_clean_c1_v1.md` §21.
