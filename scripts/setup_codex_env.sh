#!/usr/bin/env bash
# Prepare a Codex/container environment for the legacy IS_Analysis runners.
#
# This script covers the dependencies that blocked local validation in Codex:
#   - pandas for the Python batch/audit runners
#   - openpyxl for final audit Excel output
#   - synapseclient for the user-run Synapse metadata/download workflow
#   - Rscript plus the data.table package for scripts/01_prepare_exposure_fast.R
#
# Usage from the repository root:
#   bash scripts/setup_codex_env.sh
#
# Optional environment variables:
#   PYTHON_BIN=python3
#   REQUIREMENTS_FILE=requirements.txt
#   CRAN_REPO=https://cloud.r-project.org
#   PIP_INDEX_URL=https://pypi.org/simple
#   PIP_ROOT_USER_ACTION=ignore

set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python3}"
REQUIREMENTS_FILE="${REQUIREMENTS_FILE:-requirements.txt}"
CRAN_REPO="${CRAN_REPO:-https://cloud.r-project.org}"
PIP_INDEX_URL="${PIP_INDEX_URL:-https://pypi.org/simple}"
PIP_ROOT_USER_ACTION="${PIP_ROOT_USER_ACTION:-ignore}"

if [[ ! -f "${REQUIREMENTS_FILE}" ]]; then
  echo "[ERROR] Missing Python requirements file: ${REQUIREMENTS_FILE}" >&2
  exit 1
fi

python_deps_available() {
  "${PYTHON_BIN}" - <<'PY' >/dev/null 2>&1
import pandas
import openpyxl
import synapseclient
PY
}

r_deps_available() {
  command -v Rscript >/dev/null 2>&1 && \
    Rscript -e 'if (!requireNamespace("data.table", quietly = TRUE)) quit(status = 1)' >/dev/null 2>&1
}

if python_deps_available && r_deps_available; then
  echo "[OK] Dependencies already available; skipping installation"
else
  echo "[SETUP] Installing dependencies"
  if ! python_deps_available; then
    echo "[INFO] Installing Python requirements from ${REQUIREMENTS_FILE}"
    "${PYTHON_BIN}" -m pip install --root-user-action="${PIP_ROOT_USER_ACTION}" --index-url "${PIP_INDEX_URL}" --no-cache-dir --upgrade pip
    "${PYTHON_BIN}" -m pip install --root-user-action="${PIP_ROOT_USER_ACTION}" --index-url "${PIP_INDEX_URL}" --no-cache-dir -r "${REQUIREMENTS_FILE}"
  else
    echo "[OK] Python dependencies already available"
  fi

  if ! command -v Rscript >/dev/null 2>&1; then
    if command -v apt-get >/dev/null 2>&1; then
      echo "[INFO] Rscript not found; installing r-base-core with apt-get"
      export DEBIAN_FRONTEND=noninteractive
      apt-get update
      apt-get install -y --no-install-recommends r-base-core r-cran-data.table
    else
      echo "[ERROR] Rscript is not installed and apt-get is unavailable." >&2
      echo "[ERROR] Install R manually, then re-run this script to install R packages." >&2
      exit 1
    fi
  else
    echo "[INFO] Rscript found: $(command -v Rscript)"
  fi

  if ! Rscript -e 'if (!requireNamespace("data.table", quietly = TRUE)) quit(status = 1)' >/dev/null 2>&1; then
    echo "[INFO] Installing R package: data.table"
    Rscript -e "install.packages('data.table', repos='${CRAN_REPO}')"
  else
    echo "[OK] R package already available: data.table"
  fi
fi

echo "[INFO] Verifying required executables and packages"
"${PYTHON_BIN}" - <<'PY'
import pandas
import openpyxl
import synapseclient

print(f"[OK] pandas {pandas.__version__}")
print(f"[OK] openpyxl {openpyxl.__version__}")
print(f"[OK] synapseclient {synapseclient.__version__}")
PY
Rscript -e 'cat("[OK] R ", as.character(getRversion()), "\n", sep = ""); packageVersion("data.table")'

echo "[OK] Codex environment setup complete"
