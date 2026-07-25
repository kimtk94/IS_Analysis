# GIGASTROKE outcome analysis stage

`gigastroke_outcome_adapter.py` consumes local files created by
`scripts/colab_download_gigastroke_gwas.py`; it contains no downloader and must be
run only after that setup stage. Run it with:

```bash
python workflow/gigastroke_outcome_adapter.py --config config/gigastroke_outcomes.json
```

The configuration names one EUR discovery dataset and an explicit list of EAS
replication/subtype datasets. Every dataset declares its phenotype, ancestry,
role, source build, input, and source-column aliases. Copy the fixture config as
a template and point the inputs at the downloader's `GCST*_*.tsv.gz` files.

The configured UCSC chain and target FASTA are accepted only when their SHA-256
digests match. Each output row contains both untouched source coordinates,
alleles and variant ID and normalized GRCh38 coordinates/alleles. The adapter
writes a canonical TSV, a reasoned rejection TSV, and `dataset_manifest.json`.
Unmapped, multi-mapped, duplicate, and target-reference-mismatch records never
enter the canonical analysis input. The manifest records both builds and the
verified chain/reference paths and digests.

Production runs should use an authoritative GRCh37-to-GRCh38 chain and matching
GRCh38 FASTA from the same controlled reference release. Record their published
digests in config rather than copying the tiny synthetic review fixtures.
