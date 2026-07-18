# Configuration

Runtime configuration files for the IS_Analysis pipeline belong in this directory.


## Codex/container setup

Before validating the UKB-PPP batch runner in a fresh Codex or container environment, run:

```bash
bash scripts/setup_codex_env.sh
```

The setup script installs Python dependencies from `requirements.txt`, ensures `Rscript` is available via `r-base-core` when `apt-get` exists, installs the R package `data.table` required by `scripts/01_prepare_exposure_fast.R`, and includes `openpyxl` for Excel output from `scripts/00_run_full_audit_final.py`.


## Final legacy audit configuration

The self-contained final audit entrypoint is `scripts/00_run_full_audit_final.py`. The committed `config/audit_config.json` uses `"project_root": "."` so it resolves paths from the repository root in Codex/GitHub containers. If running from a Google Drive checkout, first `cd` to that project directory or pass `--config /path/to/audit_config.json` with paths adjusted for that copy. The audit writes evidence files only and does not execute MR or mutate raw data.

## Codex/PR smoke test

After the setup phase has prepared dependencies, review-time validation should run:

```bash
bash scripts/codex_smoke_test.sh
```

This smoke test intentionally does not install packages or access external data. It checks the already-prepared Python/R dependencies, materializes tiny tar archives under `${SMOKE_ROOT:-/tmp/is_analysis_smoke_fixture}` from committed TSV fixture content, writes a runtime fixture config there, and runs `scripts/00_run_full_audit_final.py` against that runtime config.

## Google Colab layout

For Colab validation, keep the real project and full data on Drive at `/content/drive/MyDrive/IS_Analysis_V2`, but keep smoke-test scratch data on the Colab VM local filesystem:

```bash
from google.colab import drive
drive.mount('/content/drive')
%cd /content/drive/MyDrive/IS_Analysis_V2

bash scripts/setup_codex_env.sh
bash scripts/colab_smoke_test.sh
python scripts/00_run_full_audit_final.py --config config/audit_config_colab_drive.json
```

`config/audit_config_colab_drive.json` points `project_root` at `/content/drive/MyDrive/IS_Analysis_V2` for real manifests/raw data/results. `scripts/colab_smoke_test.sh` sets `SMOKE_ROOT=/content/is_analysis_smoke_fixture`, so synthetic tar files and fixture audit outputs are created outside Drive and can be discarded after the smoke test.
