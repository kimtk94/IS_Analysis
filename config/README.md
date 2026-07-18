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

This smoke test intentionally does not install packages or access external data. It checks the already-prepared Python/R dependencies, materializes tiny tar archives under `/tmp` from committed TSV fixture content, and runs `scripts/00_run_full_audit_final.py` against the synthetic fixture config under `tests/fixtures`.
