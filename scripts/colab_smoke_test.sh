#!/usr/bin/env bash
# Google Colab wrapper for the install-free smoke test.
# Real project data can stay on Drive; this test materializes only tiny fixtures
# under Colab-local storage by default.
set -euo pipefail

export SMOKE_ROOT="${SMOKE_ROOT:-/content/is_analysis_smoke_fixture}"
exec bash scripts/codex_smoke_test.sh
