# Current champion: C1a license-clean atlas

Updated 2026-09-30. Four successor phases (C2–C4, plus the pending X-Atlas C5) were each
judged against a rule predeclared before its results existed. None replaced C1.

## Champion

| field | value |
|---|---|
| model | **C1a**: equal-weight, source-centred fusion of direct perturbation responses, with the promoter-neighbour cap, emitted as counts by the frozen C1 (C0-style dual-moment) generator |
| origin | our reimplementation of the AtlasShift backbone (kaipengm2/Virtual-Cell-Challenge-2026 @ `d24ce4f`, MIT, attributed). It is bit-exact against upstream on shared inputs |
| package | `outputs/competition_v2/c1_license_clean/c1_license_clean_val.vcc`, SHA-256 `fcc4e2508798805d93b9f1bb3a6957ca7cc8fbaa9f8318296161e4184bc43ea7` |
| code | commit `5a2831419ba45547046f0fec84dac86dbc21ec52` (package built at `d375f935`, with the code then uncommitted and byte-unchanged since) |
| submitted | yes, as submission #2 (validation phase) |

## Official hidden-validation result

The scores are user-reported from the leaderboard; the entry id is pending (see
`c1_official_result.md`).

| Overall | rank | PDS | MSE | NMAE | FID | REACH | JAC |
|---|---|---|---|---|---|---|---|
| **0.1394** | **370 / 1207** | 0.602 | 0.057 | 0.120 | −0.019 | 0.076 | 0.001 |

For comparison, V1 (the sparse tiered model): Overall −0.062, rank 883, PDS 0.022.

## Sources (all GREEN)

| source | context | license |
|---|---|---|
| VCC 2025 H1 (train / validation / test) | H1 hESC | CC0 1.0 |
| K562 GWPS (Replogle 2022) | K562 | CC BY 4.0 |
| CD4 DE (Marson GWCD4i) | primary CD4 T cells, 3 conditions | MIT (CZI listing; caveat recorded) |
| GENCODE v47 (promoter cap only) | — | no restrictions |

Excluded:

* **X-Atlas HCT116 / HEK293T:** BLOCKED, permission PENDING.
* **Kaden RPE1:** GREEN license, excluded on scientific grounds.
* **KOLF2.1J:** GREEN, failed C4 qualification.
* **Jurkat GSE249595:** license UNKNOWN.
* **VIPerturb-seq:** GREEN, not downloaded (low incremental value).

## Coverage

**287 / 300** Arc targets have direct evidence. By number of usable contexts:
0 / 1 / 2 / 3+ = 13 / 53 / 210 / 24.

The 13 unsupported targets are emitted as the context control plus the promoter cap:
ABCD1, ANKRD52, CAPRIN2, EPHB2, KIF21B, NICN1, PARP3, PBLD, PHF19, SEMA4F, SLC44A1, SNN,
TAF4.

## Generator

The frozen C1 generator is unchanged since C0:

* the 4-cell-averaged control template;
* the dual-moment integer emitter;
* the frozen amplitude log2fc 0.6 / bulk 0.3;
* clip 3.

## Rejected successor phases

| phase | what was tried | predeclared rule | why it failed | report |
|---|---|---|---|---|
| **C2** generator / amplitude | real-donor generators (G1c / G1ci / G1b), multinomial (G2), global amplitude 0.5–2.0, per-target amplitude | rule J: Overall up, ≥ 2 other members up, PDS ≥ 95 %, generator must pass the null-DE test | The C1 generator's "defects" are score-positive: FID rewards over-calling at ~55 % precision, and the MSE noise credit rewards its averaged template. Realistic generators matched real cell heterogeneity and raised PDS 15–30 %, yet lost Overall. The only Overall gain, C1 generator at a = 1.5, failed the null test and came from K562 alone | `c2_expression_de_calibration.md` |
| **C3** source fusion / reweighting | reliability weights, nested global weights, per-source scales, per-target reliability weights, sign-consensus down-weighting | rule O: Overall up in both folds, PDS ≥ 95 %, cosine up, **sign accuracy up**, no fold > 75 % of the gain, no material member loss | Direction, not magnitude, is the bottleneck (cosine 0.03–0.085). With three atlases, global weights are unidentifiable (transfer is symmetric within pairs). The best candidate, per-target reliability weights, gave +0.003 Overall but did not raise sign accuracy. Scaling shrank the effect and collapsed PDS | `c3_mean_response_fusion.md` |
| **C4** new direct sources | KOLF2.1J iPSC (GREEN, downloaded); Jurkat (UNKNOWN); VIPerturb-seq (GREEN, not downloaded) | source rule L: Overall up in both folds, PDS preserved, cosine up, sign accuracy up, GREEN | KOLF adds coverage (≥ 1 source 287 → 298) but its signal is pluripotency-specific. It helped H1 and hurt K562 and CD4, so cosine fell and 100 % of the gain came from H1. Jurkat's license is unknown. VIPerturb-seq would add ~27 % cells to an already covered context | `c4_new_direct_evidence_audit.md` |
| **C5** X-Atlas | predeclared only (`c5_predeclaration.md`) | see there | **Not run.** Permission PENDING | `xatlas_permission_status.md` |

## What would change the champion

Only a predeclared candidate that passes its rule on the same public folds.

* The one change prepared is **C5**: C1 plus X-Atlas, with everything else fixed.
* It may run only after written permission is recorded as APPROVED.
* Nothing is tuned on the leaderboard.
