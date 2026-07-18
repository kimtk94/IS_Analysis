#!/usr/bin/env bash
# Download the selected deCODE/Ferkingstad 2021 aptamer pQTL summary statistics.
# Intended for Google Colab setup after mounting Drive and cd'ing to the project.
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/content/drive/MyDrive/IS_Analysis_V2}"
OUTDIR="${OUTDIR:-data/rawdata/pqtl/decode_ferkingstad_2021/aptamer_sumstats}"
DECODE_DOWNLOAD_TOKEN="${DECODE_DOWNLOAD_TOKEN:-}"

if [[ -z "${DECODE_DOWNLOAD_TOKEN}" ]]; then
  echo "[ERROR] Set DECODE_DOWNLOAD_TOKEN before running this script." >&2
  echo "[ERROR] Example: export DECODE_DOWNLOAD_TOKEN='your-token-from-decode-download-url'" >&2
  exit 2
fi

cd "${PROJECT_ROOT}"
mkdir -p "${OUTDIR}"

files=(
  "9253_52_ABO_BGAT.txt.gz"
  "9253_52_ABO_BGAT.txt.md5sum"
  "18381_16_ALDH2_ALDH_E2.txt.gz"
  "18381_16_ALDH2_ALDH_E2.txt.md5sum"
  "2190_55_F11_Coagulation_Factor_XI.txt.gz"
  "2190_55_F11_Coagulation_Factor_XI.txt.md5sum"
  "3065_65_FGF5_FGF_5.txt.gz"
  "3065_65_FGF5_FGF_5.txt.md5sum"
  "6276_16_FURIN_Furin.txt.gz"
  "6276_16_FURIN_Furin.txt.md5sum"
  "4496_60_MMP12_MMP_12.txt.gz"
  "4496_60_MMP12_MMP_12.txt.md5sum"
  "10419_1_SCARA5_SCAR5.txt.gz"
  "10419_1_SCARA5_SCAR5.txt.md5sum"
  "13552_7_SWAP70_SWP70.txt.gz"
  "13552_7_SWAP70_SWP70.txt.md5sum"
  "3336_50_TFPI_TFPI.txt.gz"
  "3336_50_TFPI_TFPI.txt.md5sum"
  "8002_27_TMPRSS5_Spinesin.txt.gz"
  "8002_27_TMPRSS5_Spinesin.txt.md5sum"
)

for file in "${files[@]}"; do
  url="https://download.decode.is/s3/download?token=${DECODE_DOWNLOAD_TOKEN}&file=${file}"
  echo "[DOWNLOAD] ${file}"
  wget -c -O "${OUTDIR}/${file}" "${url}"
done

ls -lh "${OUTDIR}"
