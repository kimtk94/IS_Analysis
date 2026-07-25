# Configuration

`audit_config.json` and `audit_config_colab_drive.json` configure the
self-contained evidence-only audit (`scripts/00_run_full_audit_final.py`).

For the real-data 10-gene download, verification, and ancestry-separated
exposure workflow, see the repository [README](../README.md). The required
input is a user-supplied UKB-PPP download manifest; credentials and URLs are
not committed to this repository.

## Fixture smoke test

Review-time validation uses only committed fixtures and does not install
packages or access external data:

```bash
bash scripts/codex_smoke_test.sh
```

In a prepared Colab setup environment, run the final audit with:

```bash
python3 scripts/00_run_full_audit_final.py --config config/audit_config_colab_drive.json
```

## GIGASTROKE outcomes

`gigastroke_outcomes.example.json` enumerates all downloaded EUR/EAS ischemic
stroke and subtype datasets. Copy it to `gigastroke_outcomes.json`, select the
single EUR discovery and desired EAS replication/subtype IDs under `selection`,
and replace the reference paths and placeholder SHA-256 values with verified
local reference artifacts. Analysis is run separately from the Colab downloader;
see `workflow/README.md`.
