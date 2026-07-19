# IS_Analysis

This repository prepares paired UKB-PPP EUR/EAS pQTL exposures for downstream
analysis. The production workflow is batch-oriented and restartable: it obtains
10 paired genes at a time, validates every archive, filters instruments, and
writes EUR and EAS outputs separately.

## Real-data batch workflow

Clone the code in Colab if desired, but mount Drive and run all data, manifests,
downloads, and outputs under `/content/drive/MyDrive/IS_Analysis_V2/`. Do not
use the ephemeral Colab clone as the data workspace.

The manifest is generated from explicit Synapse parent folders; do not type gene
symbols or archive rows manually. It contains `ancestry`, the gene symbol inferred
from the archive filename, source URL/ID, the explicit `synapse_parent_id`, and
available size/checksum metadata. A `synapse_id` column makes the runner download
through `synapseclient` rather than treating the human-readable Synapse URL as a
direct archive URL.
Every selected gene must have at least one EUR and one EAS row; genes without a
pair are intentionally excluded. This TSV is also the **raw-data lifecycle
manifest**: when cleanup is enabled, the runner writes `pipeline_batch_id`,
`raw_lifecycle`, `raw_cleanup_at`, and `raw_cleanup_reason` back into this same
file. Keep it on Drive and do not delete or replace it while a run is active.

### Build the manifest from Synapse metadata

Do not manually fill archive sizes, hashes, or genes. In the data-setup
environment, explicitly provide the UKB-PPP **EUR** and **EAS** parent folders;
the builder enumerates their `.tar` files and derives the gene symbol from the
filename:

```bash
export SYNAPSE_AUTH_TOKEN='your-personal-access-token'
CODE_ROOT=/content/IS_Analysis
WORK_ROOT=/content/drive/MyDrive/IS_Analysis_V2
mkdir -p "${WORK_ROOT}/data/metadata"
python3 "${CODE_ROOT}/scripts/build_ukb_ppp_download_manifest.py" \
  --synapse-parent EUR:syn51365303 \
  --synapse-parent EAS:syn51365306 \
  --output "${WORK_ROOT}/data/metadata/ukb_ppp_download_manifest.tsv"
```

This records each file's Synapse ID, canonical Synapse URL, byte size, and MD5
without downloading the archive. For a reproducible or offline review, export
the folder metadata from Synapse and use `--synapse-metadata-file`; exported
rows must contain an ancestry and `synapse_parent_id`. The folder-query and
subsequent downloads require `synapseclient` and Synapse authentication; they
are setup/runtime operations, not review checks.

First create and review the plan without downloading:

```bash
python3 "${CODE_ROOT}/scripts/ukb_ppp_batch_manifest_runner_fast.py" \
  --base "${WORK_ROOT}/data/rawdata/pqtl/selected_targets" \
  --qc-dir "${WORK_ROOT}/results/qc/batch_pipeline" \
  --outdir "${WORK_ROOT}/results/exposure_batches" \
  --download-manifest "${WORK_ROOT}/data/metadata/ukb_ppp_download_manifest.tsv"
```

Then download, validate, and process all batches. The default batch size is 10.

```bash
python3 "${CODE_ROOT}/scripts/ukb_ppp_batch_manifest_runner_fast.py" \
  --base "${WORK_ROOT}/data/rawdata/pqtl/selected_targets" \
  --qc-dir "${WORK_ROOT}/results/qc/batch_pipeline" \
  --outdir "${WORK_ROOT}/results/exposure_batches" \
  --download-manifest "${WORK_ROOT}/data/metadata/ukb_ppp_download_manifest.tsv" \
  --batch-size 10 \
  --p-threshold 5e-8 \
  --run \
  --delete-raw-after-processing \
  --stop-on-error
```

To safely resume a single failed batch, use `--only-batch batch_003`. To only
download and validate archives, add `--download-only` with `--run`.

`--delete-raw-after-processing` is the production disk-space cleanup step. It
requires `--download-manifest`; raw-only discovery mode cannot delete files.
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
