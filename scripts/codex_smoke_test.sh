#!/usr/bin/env bash
# Install-free Codex/PR-review smoke test for the legacy audit scripts.
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python3}"
FIXTURE_CONFIG="${FIXTURE_CONFIG:-tests/fixtures/audit_config.json}"

report_env_not_ready() {
  echo "ENVIRONMENT_NOT_READY: $*" >&2
  exit 2
}

echo "[TEST] Dependency check"
"${PYTHON_BIN}" - <<'PY' || report_env_not_ready "Missing pandas/openpyxl. Run scripts/setup_codex_env.sh during setup, not during review."
import importlib.util
import sys

missing = [name for name in ("pandas", "openpyxl") if importlib.util.find_spec(name) is None]
if missing:
    raise SystemExit("Missing Python packages: " + ", ".join(missing))

import pandas
import openpyxl

print("[OK] Python:", sys.version.split()[0])
print("[OK] pandas:", pandas.__version__)
print("[OK] openpyxl:", openpyxl.__version__)
PY

command -v Rscript >/dev/null 2>&1 || report_env_not_ready "Missing Rscript. Run scripts/setup_codex_env.sh during setup, not during review."
Rscript -e '
library(data.table)
cat("[OK] R:", R.version.string, "\n")
cat("[OK] data.table:", as.character(packageVersion("data.table")), "\n")
' || report_env_not_ready "Missing R package data.table. Run scripts/setup_codex_env.sh during setup, not during review."

echo "[TEST] Python syntax"
"${PYTHON_BIN}" -m py_compile \
  scripts/00_run_full_audit_final.py \
  scripts/ukb_ppp_batch_manifest_runner_fast.py \
  scripts/audit_common.py

echo "[TEST] Fixture config JSON"
"${PYTHON_BIN}" -m json.tool "${FIXTURE_CONFIG}" >/dev/null

echo "[TEST] Materialize tiny raw tar fixtures"
rm -rf /tmp/is_analysis_smoke_fixture /tmp/is_analysis_audit_fixture_outputs
"${PYTHON_BIN}" - <<'PY'
import tarfile
from io import BytesIO
from pathlib import Path

member = Path("tests/fixtures/raw_member.tsv").read_bytes()
for ancestry in ("EUR", "EAS"):
    path = Path("/tmp/is_analysis_smoke_fixture/raw") / ancestry / "ALPHA_P12345_OID1_v1_PANEL.tar"
    path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(path, "w") as archive:
        info = tarfile.TarInfo("ALPHA.tsv")
        info.size = len(member)
        archive.addfile(info, BytesIO(member))
PY

"${PYTHON_BIN}" - <<'PY'
import tarfile
from pathlib import Path

for path in sorted(Path("/tmp/is_analysis_smoke_fixture/raw").glob("*/*.tar")):
    with tarfile.open(path, "r") as archive:
        if archive.getnames() != ["ALPHA.tsv"]:
            raise SystemExit(f"Unexpected tar members in {path}")
        print("[OK] fixture tar:", path)
PY

echo "[TEST] Audit fixture"
"${PYTHON_BIN}" scripts/00_run_full_audit_final.py --config "${FIXTURE_CONFIG}" --fail-on-review

echo "[OK] Codex smoke test completed"
