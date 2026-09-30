# C2 predeclaration, amendment 1 (before any full-fold arm was scored)

Written 2026-09-28, after a **12-target K562 smoke run** (a debugging run only: it is not
used for any selection, and it covers 3 % of that fold's targets). The frozen
predeclaration `c2_predeclaration.md` (`dcb0df82…`) is unchanged. This amendment only
**adds** one generator variant and one sensitivity anchor. It does not relax anything.

## Why

With G1's 400 donors shared by every target, the smoke null arm scored raw MSE **2.75**
(G0 null: 1.18). The mechanism is in the metric: `expr_mse_unbiased_capped` credits
sampling noise by `rho · min(C_pred, C_real)`, with `rho = min(1, b / spread)`. Here `b`
is the submission's own spread across perturbations. Shared donors make every
prediction carry the *same* sampling noise, so `b → 0`, `rho → 0`, and none of the noise
is credited. G0 escapes this because its template averages 1,600 donors, so its noise
is smaller.

## Added

* **G1ci**: G1c with **independent donors per target**. For each target, 400 cells are
  drawn without replacement from the template pool, seeded by
  `seed_for(prefix + "donors:" + target)`. Otherwise it is identical to G1c.
* It is scored in phase 1 like the other generators: a null arm and an a = 1 arm.
* A second sensitivity anchor, `ANCHOR_mean_response_G1ci`, is reported beside
  `…_G1c`.
* The sensitivity ruler of rule J (criterion 5) is the G1ci anchor if G1ci is the
  selected generator, and the G1c anchor otherwise.

Everything else, including rule J, the qualification threshold, the grids and the
G3 gate, is unchanged.
