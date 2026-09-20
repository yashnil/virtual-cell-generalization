# MSigDB gene sets

Source: Broad Institute GSEA/MSigDB data mirror
`https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/`

Release: **MSigDB 2024.1.Hs** (human, gene symbols)
Retrieved: 2026-09-20

| file | collection | sets |
|---|---|---|
| `h.all.v2024.1.Hs.symbols.gmt` | Hallmark (H) | 50 |
| `c2.cp.reactome.v2024.1.Hs.symbols.gmt` | Canonical pathways, Reactome (C2:CP:REACTOME) | 1,736 |

Identifier space: HGNC gene symbols, matching the scPertEval `var_names` space.

Gene-set definitions are used exactly as released. **They were not filtered,
edited, reweighted or selected to improve any result.** The only restriction
applied downstream is the predeclared minimum overlap with the frozen
6,640-gene response space, which is a property of our data and not of the
gene sets.

Citation: Liberzon et al., *The Molecular Signatures Database Hallmark Gene Set
Collection*, Cell Systems 1:417-425 (2015); Subramanian et al., PNAS
102:15545-15550 (2005).
