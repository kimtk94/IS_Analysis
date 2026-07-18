# Agent instructions for IS_Analysis

## Codex review verification

- Do not run `scripts/setup_codex_env.sh` during code review or PR review.
- Do not install packages during the agent/review phase.
- Run `bash scripts/codex_smoke_test.sh` for the committed fixture smoke test.
- If dependencies are missing, report `ENVIRONMENT_NOT_READY`; do not retry pip or apt installation.
- Use only the committed fixtures under `tests/fixtures` for review-time validation; runtime tar files may be materialized under `/tmp` or another local `SMOKE_ROOT`.
- Do not download or inspect full raw pQTL/GWAS data during PR review. Colab download helpers are for user-run setup only.

## Setup phase

- `bash scripts/setup_codex_env.sh` is intended only for environment setup, where network access is available.
