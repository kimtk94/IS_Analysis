# Configuration

Runtime configuration files for the IS_Analysis pipeline belong in this directory.


## Codex/container setup

Before validating the UKB-PPP batch runner in a fresh Codex or container environment, run:

```bash
bash scripts/setup_codex_env.sh
```

The setup script installs Python dependencies from `requirements.txt`, ensures `Rscript` is available via `r-base-core` when `apt-get` exists, and installs the R package `data.table` required by `scripts/01_prepare_exposure_fast.R`.
