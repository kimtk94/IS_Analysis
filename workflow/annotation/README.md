# Brain eQTL colocalization adapter

`brain_eqtl_colocalization.py` joins each configured GIGASTROKE fine-mapped
locus to **full regional** eQTL statistics by GRCh build-aware chromosome,
position, REF and ALT (including allele swaps), then runs a single-causal-variant
approximate-Bayes-factor colocalization independently for every dataset.

The config supplies a locus TSV (`locus_id`, `gene`, `genome_build`, `chromosome`, `start`,
`end`, `regional_statistics`), priors, coverage thresholds, and datasets. Each
dataset points to one registry JSON and a local TSV with `gene`, `chromosome`,
`position`, `ref`, `alt`, `beta`, and `se`. Inputs must already be on the same
build declared by the registry/locus workflow; liftover belongs upstream.

Bulk tissue and cell subtype labels occupy different result columns. Use
`evidence_family`, rather than row count, when aggregating evidence. Controlled
datasets without an approved local input must set
`availability: CONTROLLED_ACCESS_REQUIRED`; the output is then
`NOT_RUN_ACCESS_REQUIRED`, never PASS. Partial coverage is retained explicitly
and should not be interpreted like complete regional coverage.
