# DepMap 24Q4 Public — CRISPR gene effect

Acquired 2026-09-20 for the unseen-perturbation prior audit.

| | |
|---|---|
| source | Figshare article 27993248, *DepMap 24Q4 Public* |
| published | 2024-12-10 |
| licence | CC BY 4.0 |
| citation | Broad Institute DepMap, 24Q4 Public release |

## Files

| file | bytes | upstream md5 | verified |
|---|---|---|---|
| `CRISPRGeneEffect.csv` | 428,678,699 | `6edf7ade09b9b34199210b559d4745d3` | yes |
| `Model.csv` | 645,696 | `675210d17675f3517b0ce39a3c274f16` | yes |

`CRISPRGeneEffect.csv` is a cell-line x gene matrix of Chronos gene-effect
scores. `Model.csv` maps DepMap model identifiers to cell-line names and
lineages.

## Leakage assessment — read this before using DepMap as a feature

**Partial risk, and it is not zero.** DepMap gene effect *is itself a
perturbation outcome*: it is the fitness consequence of a CRISPR knockout. Three
distinctions keep it usable, and all three must hold:

1. **Different readout.** DepMap measures viability, not a transcriptome. It
   carries no information about which genes move in response to a knockdown,
   which is what `beta_p` is.
2. **Different assay.** Knockout (Cas9), not the CRISPRi knockdown the
   perturbation datasets use.
3. **Different cell lines.** Public cancer lines, not Arc's A/B/C.

Point 3 is the one that could fail. Arc's contexts are unidentified, and this
project is forbidden from inferring their identities, so we cannot verify that
A, B or C is absent from DepMap. The mitigation is that we never use a *per-cell-line*
DepMap column as a context feature: only gene-level summaries pooled across all
DepMap lines enter the model, so no context-specific quantity is ever read.
**Using a DepMap column selected by similarity to an Arc context would be an
identity inference and is prohibited.**
