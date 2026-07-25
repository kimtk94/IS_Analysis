# GIGASTROKE outcome analysis stage

`gigastroke_outcome_adapter.py` consumes local files created by
`scripts/colab_download_gigastroke_gwas.py`; it contains no downloader and must be
run only after that setup stage.

> **The command using `tests/fixtures/gigastroke/config.json` is not a real-data
> run.** Everything under `tests/fixtures/gigastroke/` is a tiny synthetic test
> input committed solely for smoke tests.

## Real-data run (Colab/Drive)

Run these as separate stages. All commands and generated data use
`/content/drive/MyDrive/IS_Analysis_V2` as the project root.

### 1. Download stage

```bash
export PROJECT_ROOT="/content/drive/MyDrive/IS_Analysis_V2"
cd "$PROJECT_ROOT"
python3 scripts/colab_download_gigastroke_gwas.py \
  --project-root "$PROJECT_ROOT" --ancestry ALL
```

This downloads the eight EUR/EAS ischemic-stroke and subtype files beneath
`data/rawdata/outcome/gigastroke_gwas_catalog/`. It does not run the adapter.

### 2. Create the real-data config

Obtain an authoritative GRCh37-to-GRCh38 UCSC chain and its matching GRCh38
FASTA during environment setup, then provide their local paths. The helper
discovers the downloaded server-controlled suffixes and writes absolute input
paths plus the locally calculated reference SHA-256 values:

```bash
python3 scripts/configure_gigastroke_outcomes.py \
  --project-root "$PROJECT_ROOT" \
  --chain "$PROJECT_ROOT/references/grch37_to_grch38.chain" \
  --reference "$PROJECT_ROOT/references/GRCh38.fa" \
  --output "$PROJECT_ROOT/config/gigastroke_outcomes.json" \
  --eur-discovery gigastroke_is_EUR \
  --eas-outcomes gigastroke_is_EAS gigastroke_las_EAS \
                 gigastroke_ces_EAS gigastroke_svs_EAS
```

The calculated digests make accidental reference replacement detectable; they
do not replace checking the files against the reference provider's published
checksums when the reference assets are first installed.

### 3. Analysis stage

```bash
python3 workflow/gigastroke_outcome_adapter.py \
  --config "$PROJECT_ROOT/config/gigastroke_outcomes.json"
```

Real canonical outputs and the dataset manifest are then written under
`data/standardized/outcome/gigastroke/`.

The configuration names one EUR discovery dataset and an explicit list of EAS
replication/subtype datasets. Every dataset declares its phenotype, ancestry,
role, source build, input, and source-column aliases. Use the example config as
the template; never use the fixture config for a scientific run.

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

## Ancestry-specific LD and fine-mapping stage

`ancestry_ld_finemap.py` is the independent Stage 04 LD/fine-mapping workflow.
It requires both EUR and EAS reference contracts, including panel name, build,
ancestry, sample provenance, matrix, and allele metadata. It will not fall back
to EUR LD for an EAS locus. Start from `config/ld_finemap.example.json` and run:

```bash
python3 workflow/ancestry_ld_finemap.py --config config/ld_finemap.json
```

The configured summary-statistics file must contain every variant in each locus,
not only significant instruments. Instrument p-value filtering is used solely
for ancestry-matched clumping. Fine-mapping uses all aligned locus variants,
creates ancestry-specific PIPs/credible sets first, and then applies the
predeclared cross-ancestry integration method. The QC table gates missingness,
allele alignment, positive-definite LD, effective sample size, and locus
coverage. The committed synthetic fixture deliberately gives EUR and EAS
different variant sets and matrices.
