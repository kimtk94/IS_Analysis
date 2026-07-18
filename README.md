# IS_Analysis

This repository prepares paired UKB-PPP EUR/EAS pQTL exposures for downstream
analysis. The production workflow is batch-oriented and restartable: it obtains
10 paired genes at a time, validates every archive, filters instruments, and
writes EUR and EAS outputs separately.

## Real-data batch workflow

Run this in the Colab/Drive data-setup environment, not during PR review. Create
a tab-separated download manifest with one row per source archive:

```tsv
ancestry	gene_symbol	source_file	url	expected_size_bytes	sha256
EUR	ALPHA	ALPHA_P12345_OID1_v1_PANEL.tar	https://example.invalid/eur.tar		
EAS	ALPHA	ALPHA_P12345_OID1_v1_PANEL.tar	https://example.invalid/eas.tar		
```

`ancestry`, `gene_symbol`, `source_file`, and `url` are required. `expected_size_bytes`
and `sha256` are optional but strongly recommended when the source provides them.
Every selected gene must have at least one EUR and one EAS row; genes without a
pair are intentionally excluded.

First create and review the plan without downloading:

```bash
python3 scripts/ukb_ppp_batch_manifest_runner_fast.py \
  --download-manifest data/metadata/ukb_ppp_download_manifest.tsv
```

Then download, validate, and process all batches. The default batch size is 10.

```bash
python3 scripts/ukb_ppp_batch_manifest_runner_fast.py \
  --download-manifest data/metadata/ukb_ppp_download_manifest.tsv \
  --batch-size 10 \
  --p-threshold 5e-8 \
  --run \
  --delete-raw-after-processing \
  --stop-on-error
```

To safely resume a single failed batch, use `--only-batch batch_003`. To only
download and validate archives, add `--download-only` with `--run`.

`--delete-raw-after-processing` is the production disk-space cleanup step. It
deletes an archive only after both ancestry batch outputs exist and the R
per-source status records that archive as successfully processed (including a
valid empty result after filtering). If an R run, output, or source status is
incomplete, the archive is retained and the batch is marked
`completed_raw_retained` rather than being deleted.

### Outputs and evidence

- Raw archives: `data/rawdata/pqtl/selected_targets/{EUR,EAS}/`
- Canonical filtered results: `results/exposure_batches/EUR/exposure_batch_###.tsv`
  and `results/exposure_batches/EAS/exposure_batch_###.tsv`
- Batch state: `results/qc/batch_pipeline/batch_manifest.tsv`
- Per-batch download/verification evidence:
  `results/qc/batch_pipeline/downloads/batch_###.tsv`
- Raw deletion/retention evidence (when cleanup is requested):
  `results/qc/batch_pipeline/raw_cleanup/batch_###.tsv`
- Per-gene processing state:
  `results/exposure_batches/{EUR,EAS}/logs/batch_###_gene_status.tsv`

The R preparation stage applies cis filtering when gene coordinates are
available, then retains `pval < --p-threshold` and `F_stat > 10`. It writes a
schema-preserving empty TSV when a batch has no qualifying instruments.

## Review smoke test

The committed fixture-only validation never downloads real data:

```bash
bash scripts/codex_smoke_test.sh
```
