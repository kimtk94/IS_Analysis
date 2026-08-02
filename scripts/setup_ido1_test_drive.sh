#!/usr/bin/env bash
# Prepare a Drive-backed workspace for the 15-gene IDO1 focused test.
set -euo pipefail

WORK_ROOT="${WORK_ROOT:-/content/drive/MyDrive/IS_Analysis_V3}"
SOURCE_WORK_ROOT="${SOURCE_WORK_ROOT:-/content/drive/MyDrive/IS_Analysis_V2}"

drive_root="/content/drive/MyDrive"
if [[ ! -d "${drive_root}" ]]; then
  echo "[ERROR] Google Drive is not mounted at ${drive_root}" >&2
  exit 1
fi

directories=(
  "${WORK_ROOT}/data/metadata"
  "${WORK_ROOT}/data/reference"
  "${WORK_ROOT}/data/rawdata/pqtl/selected_targets/EUR"
  "${WORK_ROOT}/data/rawdata/pqtl/selected_targets/EAS"
  "${WORK_ROOT}/results/qc/ido1_test_pipeline"
  "${WORK_ROOT}/results/test/ido1/exposure_batches"
  "${WORK_ROOT}/results/test/ido1/standardized/pqtl"
  "${WORK_ROOT}/results/test/ido1/instrument_candidates"
)

for directory in "${directories[@]}"; do
  mkdir -p "${directory}"
done

# Reuse only the small setup inputs from V2 when V3 does not have them. Raw
# archives and derived results are never copied by this setup helper.
declare -A setup_inputs=(
  ["data/metadata/ukb_ppp_download_manifest.tsv"]="download manifest"
  ["data/reference/gene_coordinates_hg38.tsv"]="GRCh38 gene coordinates"
)

for relative_path in "${!setup_inputs[@]}"; do
  destination="${WORK_ROOT}/${relative_path}"
  source="${SOURCE_WORK_ROOT}/${relative_path}"
  if [[ -s "${destination}" ]]; then
    echo "[KEEP] ${setup_inputs[${relative_path}]}: ${destination}"
  elif [[ -s "${source}" ]]; then
    cp "${source}" "${destination}"
    echo "[COPY] ${setup_inputs[${relative_path}]}: ${source} -> ${destination}"
  else
    echo "[WARN] Missing ${setup_inputs[${relative_path}]}; create it before running the test: ${destination}" >&2
  fi
done

env_file="${WORK_ROOT}/ido1_test.env"
cat >"${env_file}" <<EOF
export WORK_ROOT='${WORK_ROOT}'
export IDO1_TEST_QC_DIR='${WORK_ROOT}/results/qc/ido1_test_pipeline'
export IDO1_TEST_OUTDIR='${WORK_ROOT}/results/test/ido1/exposure_batches'
export IDO1_TEST_STANDARDIZED_DIR='${WORK_ROOT}/results/test/ido1/standardized/pqtl'
export IDO1_TEST_INSTRUMENT_DIR='${WORK_ROOT}/results/test/ido1/instrument_candidates'
EOF

echo "[OK] Drive-backed IDO1 test workspace is ready: ${WORK_ROOT}"
echo "[OK] Environment file: ${env_file}"
echo "[NEXT] source '${env_file}'"
