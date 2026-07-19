# IS_Analysis

This repository prepares paired UKB-PPP EUR/EAS pQTL exposures for downstream
analysis. The production workflow is batch-oriented and restartable: it obtains
10 paired genes at a time, validates every archive, filters instruments, and
writes EUR and EAS outputs separately.

## Real-data batch workflow

Run the production workflow in Google Colab. Clone code to
`/content/IS_Analysis_V2`, mount Drive, and run all data, manifests, downloads,
and outputs under `/content/drive/MyDrive/IS_Analysis_V2/`. Do not use the
ephemeral Colab clone as the data workspace.

The manifest is generated from explicit Synapse parent folders; do not type gene
symbols or archive rows manually. Initial creation records only `ancestry`, the
gene symbol inferred from the archive filename, source URL/ID, and the explicit
`synapse_parent_id`. Before each 10-gene batch, the runner retrieves size and
checksum metadata only for that batch's Synapse files (maximum 8 concurrent
requests), persists it in the manifest, then downloads and verifies the files.
Every selected gene must have at least one EUR and one EAS row; genes without a
pair are intentionally excluded. This TSV is also the **raw-data lifecycle
manifest**: when cleanup is enabled, the runner writes `pipeline_batch_id`,
`raw_lifecycle`, `raw_cleanup_at`, and `raw_cleanup_reason` back into this same
file. Keep it on Drive and do not delete or replace it while a run is active.

### Build the manifest from Synapse metadata

Do not manually fill archive sizes, hashes, or genes. In the data-setup
environment, explicitly provide the UKB-PPP **EUR** and **EAS** parent folders;
the builder enumerates their `.tar` files and derives the gene symbol from the
filename. Execute the following **Colab cells in order**.

**1. Python cell — mount Google Drive.**

```python
from google.colab import drive

drive.mount("/content/drive")
```

**2. Bash cell — install all runtime dependencies, including `synapseclient`.**

```bash
%%bash
set -euo pipefail

CODE_ROOT="/content/IS_Analysis_V2"
cd "${CODE_ROOT}"
bash scripts/setup_codex_env.sh
```

`scripts/setup_codex_env.sh` installs Python packages from `requirements.txt`,
which includes `synapseclient>=4.9` (the non-deprecated Synapse child-listing
API), and installs the R dependencies needed for the
batch preparation stage. Run this only in the user-run Colab setup environment,
never as part of fixture-only PR review; see [AGENTS.md](AGENTS.md).

**3. Python cell — set the Synapse personal access token without printing it.**

```python
import os
from getpass import getpass

os.environ["SYNAPSE_AUTH_TOKEN"] = getpass("Synapse personal access token: ")
```

**4. Bash cell — create the EUR/EAS manifest in Drive.**

```bash
%%bash
set -euo pipefail

CODE_ROOT="/content/IS_Analysis_V2"
WORK_ROOT="/content/drive/MyDrive/IS_Analysis_V2"
SCRIPT="${CODE_ROOT}/scripts/build_ukb_ppp_download_manifest.py"

if [[ ! -f "${SCRIPT}" ]]; then
  echo "[ERROR] 스크립트를 찾을 수 없습니다: ${SCRIPT}" >&2
  echo "먼저 GitHub 저장소가 ${CODE_ROOT}에 clone되었는지 확인하세요." >&2
  exit 1
fi

: "${SYNAPSE_AUTH_TOKEN:?SYNAPSE_AUTH_TOKEN을 먼저 설정하세요.}"
mkdir -p "${WORK_ROOT}/data/metadata"

python3 "${SCRIPT}" \
  --synapse-parent "EUR:syn51365303" \
  --synapse-parent "EAS:syn51365306" \
  --output "${WORK_ROOT}/data/metadata/ukb_ppp_download_manifest.tsv"

echo "[DONE] Manifest 생성 완료"
echo "${WORK_ROOT}/data/metadata/ukb_ppp_download_manifest.tsv"
```

This records each file's Synapse ID, canonical Synapse URL, and parent ID without
downloading archives or making a per-file checksum request. The runner lazily
retrieves size/MD5 metadata for the active batch and emits progress every 100
metadata lookups. For a reproducible or offline review, export
the folder metadata from Synapse and use `--synapse-metadata-file`; exported
rows must contain an ancestry and `synapse_parent_id`. The folder-query and
subsequent downloads require `synapseclient` and Synapse authentication; they
are setup/runtime operations, not review checks.

First create and review the plan without downloading:

```bash
%%bash
set -euo pipefail

CODE_ROOT="/content/IS_Analysis_V2"
WORK_ROOT="/content/drive/MyDrive/IS_Analysis_V2"
cd "${CODE_ROOT}"

python3 "${CODE_ROOT}/scripts/ukb_ppp_batch_manifest_runner_fast.py" \
  --base "${WORK_ROOT}/data/rawdata/pqtl/selected_targets" \
  --qc-dir "${WORK_ROOT}/results/qc/batch_pipeline" \
  --outdir "${WORK_ROOT}/results/exposure_batches" \
  --download-manifest "${WORK_ROOT}/data/metadata/ukb_ppp_download_manifest.tsv"
```

Then download, validate, and process all batches. The default batch size is 10.

```bash
%%bash
set -euo pipefail

CODE_ROOT="/content/IS_Analysis_V2"
WORK_ROOT="/content/drive/MyDrive/IS_Analysis_V2"
cd "${CODE_ROOT}"

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

The runner flushes current batch and phase messages to stdout. If Colab's
`%%bash` output remains buffered in the notebook UI, use `python3 -u` in place
of `python3` in the same command.

### Reuse existing raw archives

If the required archives already exist at a **different location**, keep the
current run/staging directory in `--base` and specify that separate location
with `--existing-raw-base`. The runner checks the ancestry-specific source path
(`EUR/` or `EAS/`), requires a non-empty valid tar, and creates a symlink in
`--base` only for valid archives. Raw cleanup later removes the staging symlink,
not the original archive. A staged archive skips its Synapse metadata request
and download, then receives normal tar/checksum validation in the batch workflow
(including any checksum already stored in the manifest). For example, when raw
data are in the Colab clone but run outputs/staging are in Drive:

```bash
--base "${WORK_ROOT}/data/rawdata/pqtl/selected_targets" \
--existing-raw-base "${CODE_ROOT}/data/rawdata/pqtl/selected_targets"
```

Google Drive FUSE may reject symlink creation. In that case the runner reports
the fallback and atomically copies only the selected archive into `--base`; its
source archive remains unchanged unless optional source deletion is requested.

This check happens at the start of every selected batch. If an archive is absent
or invalid in `--existing-raw-base`, the runner follows the normal Synapse
metadata/download path for that file. This means existing raw data and newly
downloaded data can be mixed safely within the same 10-gene batch.

With `--delete-raw-after-processing`, the default cleanup removes only the
staging symlink and preserves the separate original. To delete the separate
original as well after its 10-gene batch is successfully processed, explicitly
add `--delete-existing-raw-after-processing`. That option requires both
`--existing-raw-base` and `--delete-raw-after-processing` and records the source
deletion in the batch cleanup evidence.

Keep `--download-manifest` when lifecycle tracking or raw cleanup is required;
the manifest source filename must match the existing archive filename. Without a
manifest, the runner can discover paired EUR/EAS raw archives, but cannot use
`--delete-raw-after-processing` because it has no lifecycle record.

To safely resume a single failed batch, use `--only-batch batch_003`. To only
download and validate archives, add `--download-only` with `--run`.

### Restart behavior

The runner writes the batch status after every terminal batch state. On a later
`--run`, it restores that state from `batch_manifest.tsv` and automatically skips
`completed`, `completed_raw_deleted`, and `completed_raw_retained` batches. This
prevents already processed data from being downloaded or analysed again after a
Colab interruption. Failed or incomplete batches remain eligible for the next
run. Use `--only-batch batch_003` to resume one failed batch, or add
`--rerun-completed` only when deliberately rerunning completed work.

Before retrying a batch saved as `running`, metadata/download failed, or
processing failed, the runner deletes only stale **derived** batch outputs and
per-gene status files, then recreates them. Verified raw archives are retained
and reused; they are not discarded until normal successful raw cleanup. Each
reset is recorded in
`results/qc/batch_pipeline/partial_cleanup/batch_###.tsv`.

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
- Run progress: `results/qc/batch_pipeline/batch_progress.tsv` (batch number,
  current phase, and terminal status)
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
