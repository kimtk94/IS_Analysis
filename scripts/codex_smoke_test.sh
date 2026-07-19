#!/usr/bin/env bash
# Install-free Codex/PR-review smoke test for the legacy audit scripts.
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python3}"
FIXTURE_CONFIG="${FIXTURE_CONFIG:-tests/fixtures/audit_config.json}"
SMOKE_ROOT="${SMOKE_ROOT:-/tmp/is_analysis_smoke_fixture}"
SMOKE_RAW_DIR="${SMOKE_RAW_DIR:-${SMOKE_ROOT}/raw}"
SMOKE_OUTPUT_DIR="${SMOKE_OUTPUT_DIR:-${SMOKE_ROOT}/audit_outputs}"
RUNTIME_CONFIG="${RUNTIME_CONFIG:-${SMOKE_ROOT}/audit_config.runtime.json}"
export FIXTURE_CONFIG SMOKE_ROOT SMOKE_RAW_DIR SMOKE_OUTPUT_DIR RUNTIME_CONFIG

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
  scripts/build_ukb_ppp_download_manifest.py \
  scripts/synapse_metadata.py \
  scripts/ukb_ppp_batch_manifest_runner_fast.py \
  scripts/colab_download_gigastroke_gwas.py

echo "[TEST] Fixture config JSON"
"${PYTHON_BIN}" -m json.tool "${FIXTURE_CONFIG}" >/dev/null

echo "[TEST] Build manifest from committed Synapse metadata fixture"
MANIFEST_FIXTURE_OUTPUT="${SMOKE_ROOT}/ukb_ppp_download_manifest.tsv"
mkdir -p "${SMOKE_ROOT}"
"${PYTHON_BIN}" scripts/build_ukb_ppp_download_manifest.py \
  --synapse-metadata-file tests/fixtures/synapse_metadata.tsv \
  --output "${MANIFEST_FIXTURE_OUTPUT}"
"${PYTHON_BIN}" - "${MANIFEST_FIXTURE_OUTPUT}" <<'PY'
import csv
import sys

with open(sys.argv[1], newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
assert [row["synapse_id"] for row in rows] == ["syn1001", "syn1002"]
assert [row["expected_size_bytes"] for row in rows] == ["", ""]
assert [row["ancestry"] for row in rows] == ["EUR", "EAS"]
assert [row["synapse_parent_id"] for row in rows] == ["syn51365303", "syn51365306"]
assert all(row["url"].startswith("https://www.synapse.org/Synapse:syn") for row in rows)
print("[OK] Synapse metadata manifest fixture")
PY

echo "[TEST] Materialize tiny raw tar fixtures under ${SMOKE_RAW_DIR}"
rm -rf "${SMOKE_ROOT}"
mkdir -p "${SMOKE_ROOT}"
"${PYTHON_BIN}" - <<'PY'
import os
import tarfile
from io import BytesIO
from pathlib import Path

member = Path("tests/fixtures/raw_member.tsv").read_bytes()
raw_dir = Path(os.environ["SMOKE_RAW_DIR"])
for ancestry in ("EUR", "EAS"):
    path = raw_dir / ancestry / "ALPHA_P12345_OID1_v1_PANEL.tar"
    path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(path, "w") as archive:
        info = tarfile.TarInfo("ALPHA.tsv")
        info.size = len(member)
        archive.addfile(info, BytesIO(member))
PY

"${PYTHON_BIN}" - <<'PY'
import os
import tarfile
from pathlib import Path

raw_dir = Path(os.environ["SMOKE_RAW_DIR"])
for path in sorted(raw_dir.glob("*/*.tar")):
    with tarfile.open(path, "r") as archive:
        if archive.getnames() != ["ALPHA.tsv"]:
            raise SystemExit(f"Unexpected tar members in {path}")
        print("[OK] fixture tar:", path)
PY

echo "[TEST] Build runtime fixture config at ${RUNTIME_CONFIG}"
"${PYTHON_BIN}" - <<'PY'
import json
import os
from pathlib import Path

fixture_config = Path(os.environ["FIXTURE_CONFIG"])
runtime_config = Path(os.environ["RUNTIME_CONFIG"])
config = json.loads(fixture_config.read_text(encoding="utf-8"))
config["raw_dirs"] = [os.environ["SMOKE_RAW_DIR"]]
config["audit_output_dir"] = os.environ["SMOKE_OUTPUT_DIR"]
runtime_config.parent.mkdir(parents=True, exist_ok=True)
runtime_config.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
PY
"${PYTHON_BIN}" -m json.tool "${RUNTIME_CONFIG}" >/dev/null

echo "[TEST] Audit fixture"
"${PYTHON_BIN}" scripts/00_run_full_audit_final.py --config "${RUNTIME_CONFIG}" --fail-on-review

echo "[TEST] Same-gene multi-archive exposure preservation"
DUPLICATE_RAW_DIR="${SMOKE_ROOT}/duplicate_gene_raw"
DUPLICATE_OUT_DIR="${SMOKE_ROOT}/duplicate_gene_output"
DUPLICATE_GENE_FILE="${SMOKE_ROOT}/duplicate_gene.txt"
export DUPLICATE_RAW_DIR DUPLICATE_OUT_DIR DUPLICATE_GENE_FILE
"${PYTHON_BIN}" - <<'PY'
import os
import tarfile
from io import BytesIO
from pathlib import Path

member = Path("tests/fixtures/raw_member.tsv").read_bytes()
raw_dir = Path(os.environ["DUPLICATE_RAW_DIR"])
for filename in ("ALPHA_P12345_OID1_v1_PANEL.tar", "ALPHA_Q99999_OID2_v1_PANEL.tar"):
    path = raw_dir / "EUR" / filename
    path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(path, "w") as archive:
        info = tarfile.TarInfo("ALPHA.tsv")
        info.size = len(member)
        archive.addfile(info, BytesIO(member))
Path(os.environ["DUPLICATE_GENE_FILE"]).write_text("ALPHA\n", encoding="utf-8")
PY
Rscript scripts/01_prepare_exposure_fast.R \
  --gene-file "${DUPLICATE_GENE_FILE}" \
  --batch-id duplicate_gene \
  --outdir "${DUPLICATE_OUT_DIR}" \
  --rawdir "${DUPLICATE_RAW_DIR}" \
  --ancestries EUR \
  --no-cis-filter
"${PYTHON_BIN}" - "${DUPLICATE_OUT_DIR}/exposure_duplicate_gene.tsv" <<'PY'
import csv
import sys

with open(sys.argv[1], newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
assert len(rows) == 4, f"Expected 4 rows from two archives, found {len(rows)}"
assert {row["source_file"] for row in rows} == {
    "ALPHA_P12345_OID1_v1_PANEL.tar",
    "ALPHA_Q99999_OID2_v1_PANEL.tar",
}
exposures = {row["id.exposure"] for row in rows}
assert len(exposures) == 2, f"Expected distinct exposure IDs per archive, found {exposures}"
assert all(row["gene_symbol"] == "ALPHA" and row["ancestry"] == "EUR" for row in rows)
print("[OK] Both same-gene archives retained with distinct exposure IDs")
PY

echo "[OK] Codex smoke test completed"
