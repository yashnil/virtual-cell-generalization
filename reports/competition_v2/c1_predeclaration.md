# C1 predeclaration: folds, arms, agreement statistic, shrinkage family, pass rule

**Written 2026-09-26, before any C1a/C1b arm was scored on any fold.** The scoring script
(`scripts/competition_v2/run_c1_public_folds.py`) records this file's SHA-256 in every
output, and the final report checks that the hash has not changed. Nothing below may be
edited after scoring starts. A later change requires a new, separately dated
predeclaration.

## 1. Source universe (C1)

GREEN direct-response sources only, from `reports/competition_v2/data_license_register.md`:

* K562 GWPS (CC BY 4.0);
* VCC 2025 H1, all three splits (CC0 1.0);
* CD4 genome-wide DE statistics (MIT, as listed by the CZI Virtual Cells Platform).

X-Atlas HCT116 / HEK293T are BLOCKED_PENDING_PERMISSION. They enter **only** the
X-Atlas ablation arms (`C0_withX`), which are reported under §13 and never used for C1
selection. Kaden 2025 RPE1 is license-GREEN but not part of the C0 backbone. It is
excluded as a source because of weak reproducibility in our frozen results and an
unverified CRISPRa/CRISPRi modality.

Per-source statistics come from our own raw-file pass
(`outputs/competition_v2/c1_license_clean/sources/`), which is verified bit-identical to
the upstream prepared files.

## 2. Public folds (leave-one-atlas-out)

A held-out source is never loaded as a predictor. Its responses are read only by the
scorer, and never by agreement, centering, shrinkage or selection.

| fold | truth | panel | gene axis | predictors (C1) | scored as |
|---|---|---|---|---|---|
| **H1** | H1 2025 train cells (≤ 1,000 per target; 4,000 reference controls, 4,000 template controls; the same seeded draw as the C0 public-H1 harness) | the 150 H1-train targets | official ∩ H1 | K562, CD4 | full local cell-eval (6 members), local anchors |
| **K562** | K562 GWPS cells (≤ 1,000 per target; 8,000 random controls split 4,000 / 4,000 by `rng(SEED+1)`) | K562 retained targets with ≥ 20 cells | official ∩ K562 | H1, CD4 | full local cell-eval, local anchors |
| **CD4** | CD4 centred log2FC (per-condition source-mean centering, conditions averaged) | retained targets with ≥ 1 quality-passing condition | official ∩ CD4 | K562, H1 | effect-level PDS only (`pds_cosine` on log2FC effects) |

The folds are HCT116 and HEK293T for no C1 purpose, because they are BLOCKED.

## 3. Arms (every full-cell fold)

| arm | definition |
|---|---|
| `ANCHOR_mean_response` | ORACLE local baseline b: the held-out panel's true mean response through the same emitter |
| `ANCHOR_split_half` | ORACLE local replicate r |
| `G0_control` | emitter with zero effect |
| **`C1a`** | equal-weight (every GREEN source weight 1, CD4 weight 1) centred fusion. Frozen amplitude 0.6 / 0.3, clip 3, promoter cap on |
| `C1a_no_promoter` | C1a without the promoter cap (§14) |
| `S1_l{0.25,0.5,0.75}` | C1a × global λ (diagnostic control, see §6) |
| **`S2_f{0,0.25,0.5}`** | C1a × agreement-dependent λ_p (the C1b family) |
| `C0_noX` | upstream weights (K562 2, H1 2, CD4 0.5), X-Atlas removed, promoter on |
| `C0_withX` | upstream weights including X-Atlas HCT116 1 / HEK293T 1 (X-Atlas ablation only) |

Scaling: s = (u − b)/(r − b) per member. MSE is clamped to [0, 1] as officially, and
local Overall is the mean of the six members.

## 4. Source-agreement statistic (exact; no alternative will be substituted)

For target p, let S_p be the predictor sources usable for p in that fold (≥ 20 cells, or a
CD4 quality-passing condition). For each source, v_s is its **centred `log2fc` effect**
under the frozen per-source centering. `log2fc` is the only space the CD4 DE table
has. For each pair (s, t) in S_p:

cos_st = ⟨v_s, v_t⟩ / (‖v_s‖ ‖v_t‖), computed over the genes both sources measured,
minus **every panel target gene** (the PDS exclusion set, which includes p itself).
A pair where either vector has zero norm on that set is skipped.

**agreement a_p = mean of cos_st over the pairs.** It is undefined (NaN) when |S_p| < 2.
Implementation: `virtual_cell.competition_v2.fusion.source_agreement`.

## 5. Shrinkage family (extremely small)

β_final[p] = λ_p · β_equal_mean[p]. λ_p multiplies the fused effect in both spaces,
before the frozen amplitude and clip. The promoter cap is applied after it, unchanged.
The direction of an effect is never altered.

* **S0**: λ = 1 (= C1a).
* **S1**: λ ∈ {0.25, 0.5, 0.75}, one global scalar.
* **S2**: λ_p = f + (1 − f) · F(a_p), f ∈ {0, 0.25, 0.5}. Here F is the empirical CDF
  midrank ((rank − 0.5)/n) of a_p among the panel targets with a defined agreement, and
  F = 0.5 when a_p is undefined. This is monotone non-decreasing in agreement. The one
  parameter f is a floor, and the mean λ is (1 + f)/2, which S1 matches at 0.5 and 0.75.

No neural nets, trees or other nonlinear fitting.

## 6. Selection (nested; public folds only)

* **C1b := S2.** Only S2 uses source agreement. S1 is a diagnostic control: it separates
  "any shrinkage helps" from "agreement allocation helps". It is **not** eligible to
  replace C1a in this phase, because a global amplitude change is not the authorised
  modification.
* For outer fold k ∈ {H1, K562}, f\*_k is the f that maximises local Overall on the other
  full-cell fold (inner holdout). The outer-fold score of C1b is S2 at f\*_k. CD4 cannot
  select f, because per-target positive scaling leaves effect-level cosine unchanged.
* If C1b passes §7, the Arc candidate uses the f that maximises mean local Overall over
  both full-cell folds.

## 7. Pass rule (C1b replaces C1a only if ALL hold; strict inequalities)

1. Mean outer-fold local Overall over {H1, K562}: C1b > C1a.
2. PDS: C1b > C1a in at least ⌈0.75 × 3⌉ = **3 of the 3** PDS folds (H1, K562 raw
   `pds_cosine`; CD4 effect-level PDS). A tie is not "better".
3. No material worsening: for each of NMAE, FID, REACH and JAC, the mean over the
   full-cell folds of (C1b − C1a) local-scaled is ≥ −0.02.
4. Not driven by one context: C1b − C1a Overall > 0 in **both** H1 and K562.
5. No BLOCKED source in any C1b input (asserted in code).

Otherwise **C1a** is the C1 candidate.

**Expectation, stated in advance.** A positive per-target scale leaves an effect's cosine
unchanged. S2 is therefore exactly PDS-neutral at the effect level (CD4 fold: a tie, by
construction) and nearly PDS-neutral in counts. Criterion 2 is expected to fail unless the
count-level nonlinearity produces strict PDS gains in all three folds. Any benefit of S2
would have to show in MSE and the DE members. This expectation is recorded so that a
failure of criterion 2 cannot be reinterpreted afterwards.

## 8. Promoter correction

Retained exactly as frozen (5 kb window, 15 % floor, log ramp) **only if** GENCODE v47 is
GREEN in the license register. It is reported as the `C1a` vs `C1a_no_promoter` difference.
It is not retuned.

## 9. X-Atlas ablation (not used for any C1 decision)

`C0_withX` vs `C0_noX` on the H1 and K562 folds, plus the CD4 effect-level fold. The
ablation reports target coverage, mean sources per target, PDS, Overall and the DE
members.
