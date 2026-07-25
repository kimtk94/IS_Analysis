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

### 0. Reference setup (run once)

The adapter uses the official UCSC hg19 source FASTA, hg38 target FASTA, and
matching chain. Download and verify them
with the committed setup helper:

```bash
export PROJECT_ROOT="/content/drive/MyDrive/IS_Analysis_V2"
cd "$PROJECT_ROOT"
bash scripts/colab_download_grch38_liftover_references.sh
```

The helper downloads `hg38ToHg19.over.chain.gz` from the UCSC hg38 liftOver
directory, `hg19.fa.gz` from UCSC hg19 bigZips, and `hg38.fa.gz` from UCSC hg38
bigZips, then verifies all three against the MD5 lists published in those same directories. UCSC chain
filenames are `targetToQuery`: `hg38ToHg19` is therefore the chain whose query
coordinates are hg19/GRCh37 and whose target coordinates are hg38/GRCh38. The
helper writes the decompressed files expected by the config commands:

```text
/content/drive/MyDrive/IS_Analysis_V2/references/grch37_to_grch38.chain
/content/drive/MyDrive/IS_Analysis_V2/references/GRCh37.fa
/content/drive/MyDrive/IS_Analysis_V2/references/GRCh38.fa
```

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

Use the chain and matching FASTA installed in stage 0. The config helper
discovers the downloaded server-controlled suffixes and writes absolute input
paths plus the locally calculated reference SHA-256 values:

```bash
python3 scripts/configure_gigastroke_outcomes.py \
  --project-root "$PROJECT_ROOT" \
  --chain "$PROJECT_ROOT/references/grch37_to_grch38.chain" \
  --source-reference "$PROJECT_ROOT/references/GRCh37.fa" \
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

Some GIGASTROKE files contain only `effect_allele` and `other_allele`, with no
REF/ALT or variant ID. The adapter first tries a build-matching
`chromosome:position:ref:alt` source variant ID. If none exists, it compares
both reported alleles with the verified GRCh37 FASTA and assigns REF only when
exactly one allele matches the source reference. It never assumes that effect
allele means ALT. Rows for which REF cannot be determined unambiguously are
excluded as `source_reference_allele_mismatch`.

Production runs should use an authoritative GRCh37-to-GRCh38 chain and matching
GRCh38 FASTA from the same controlled reference release. Record their published
digests in config rather than copying the tiny synthetic review fixtures.
