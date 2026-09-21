# STRING v12.0 — Homo sapiens protein association network

Acquired 2026-09-20 for the unseen-perturbation prior audit.

| | |
|---|---|
| source | https://stringdb-downloads.org/download/ |
| organism | 9606 (Homo sapiens) |
| version | v12.0 |
| licence | CC BY 4.0 |
| citation | Szklarczyk et al., *Nucleic Acids Research* (2023) |

## Files

| file | bytes | upstream last-modified |
|---|---|---|
| `9606.protein.links.v12.0.txt.gz` | 83,164,437 | 2023-05-16 |
| `9606.protein.info.v12.0.txt.gz` | 1,970,090 | 2023-05-31 |

`protein.links` carries `protein1 protein2 combined_score` over Ensembl protein
identifiers; `protein.info` supplies the preferred gene name needed to map them
onto the HGNC symbols the perturbation data uses.

## Leakage assessment

**Cannot leak Arc perturbation outcomes.** STRING is an aggregate of prior
literature, co-expression, co-occurrence and database curation, published in
2023 and fixed. It contains no measurement from the Arc 2026 panel and no
transcriptional perturbation response from any context in this study. Its
co-expression channel is derived from public expression compendia, which is a
prior over gene relatedness, not an outcome of a knockdown.
