# Configuration

Runtime configuration files for the IS_Analysis pipeline belong in this directory.


## Codex/container setup

Before validating the UKB-PPP batch runner in a fresh Codex or container environment, run:

```bash
bash scripts/setup_codex_env.sh
```

The setup script installs Python dependencies from `requirements.txt`, ensures `Rscript` is available via `r-base-core` when `apt-get` exists, installs the R package `data.table` required by `scripts/01_prepare_exposure_fast.R`, and includes `openpyxl` for Excel output from `scripts/00_run_full_audit_final.py`.


## Final legacy audit configuration

The self-contained final audit entrypoint is `scripts/00_run_full_audit_final.py`. By default it expects `config/audit_config.json`; pass `--config /path/to/audit_config.json` when running against a Google Drive project copy. The audit writes evidence files only and does not execute MR or mutate raw data.
