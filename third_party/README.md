# third_party/

External code vendored **verbatim** for the competition track only. Nothing in
`src/virtual_cell/` imports from here; competition-v2 code runs it as an
external program so our scientific modules never absorb it silently.

## atlasshift/

| field | value |
|---|---|
| upstream | https://github.com/kaipengm2/Virtual-Cell-Challenge-2026 |
| commit | `d24ce4fdae0cd7cba3ba29546cd8737094eae98a` (2026-09-23, "update readme") |
| license | MIT (`atlasshift/LICENSE`, "Copyright (c) 2026 Virtual Cell Challenge 2026 solution contributors") |
| retrieved | 2026-09-25 |
| modifications | none (the upstream `.git` directory was removed; files are byte-identical) |
| method name | "A Simple Atlas Transfer Method (AtlasShift)" |
| reported score | 0.1545618019, public leaderboard #82, 2026-09-10 (upstream README; not reproduced by us) |

External **data** retain their own terms; see
`reports/competition_v2/atlas_shift_source_audit.md`. The runtime environment is
`atlasshift/.venv` (Python 3.13, upstream `requirements.txt` pins, plus `slafdb`),
which is git-ignored.
