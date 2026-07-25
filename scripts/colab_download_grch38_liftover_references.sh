#!/usr/bin/env bash
# User-run setup stage: download the matching UCSC hg38 chain and FASTA.
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/content/drive/MyDrive/IS_Analysis_V2}"
REFERENCE_DIR="${REFERENCE_DIR:-${PROJECT_ROOT}/references}"
CHAIN_BASE="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver"
FASTA_BASE="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips"

mkdir -p "${REFERENCE_DIR}"
tmp_dir="$(mktemp -d)"
trap 'rm -rf "${tmp_dir}"' EXIT

echo "[DOWNLOAD] UCSC hg38-to-hg19 chain and published MD5 list"
curl --fail --location --retry 3 -o "${tmp_dir}/hg38ToHg19.over.chain.gz" \
  "${CHAIN_BASE}/hg38ToHg19.over.chain.gz"
curl --fail --location --retry 3 -o "${tmp_dir}/chain.md5sum.txt" \
  "${CHAIN_BASE}/md5sum.txt"

echo "[DOWNLOAD] UCSC hg38 FASTA and published MD5 list"
curl --fail --location --retry 3 -o "${tmp_dir}/hg38.fa.gz" \
  "${FASTA_BASE}/hg38.fa.gz"
curl --fail --location --retry 3 -o "${tmp_dir}/fasta.md5sum.txt" \
  "${FASTA_BASE}/md5sum.txt"

echo "[VERIFY] UCSC-published MD5 checksums"
(
  cd "${tmp_dir}"
  awk '$2 == "hg38ToHg19.over.chain.gz" {print}' chain.md5sum.txt > chain.selected.md5
  awk '$2 == "hg38.fa.gz" {print}' fasta.md5sum.txt > fasta.selected.md5
  test -s chain.selected.md5 && test -s fasta.selected.md5
  md5sum --check chain.selected.md5
  md5sum --check fasta.selected.md5
)

# UCSC chain names are targetToQuery.  Thus hg38ToHg19 has hg38 as target
# and hg19/GRCh37 as query, which is the direction required by this adapter.
gzip -dc "${tmp_dir}/hg38ToHg19.over.chain.gz" > "${REFERENCE_DIR}/grch37_to_grch38.chain"
gzip -dc "${tmp_dir}/hg38.fa.gz" > "${REFERENCE_DIR}/GRCh38.fa"

cat > "${REFERENCE_DIR}/GIGASTROKE_REFERENCE_SOURCES.txt" <<EOF
chain_url=${CHAIN_BASE}/hg38ToHg19.over.chain.gz
chain_md5_url=${CHAIN_BASE}/md5sum.txt
reference_url=${FASTA_BASE}/hg38.fa.gz
reference_md5_url=${FASTA_BASE}/md5sum.txt
downloaded_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
EOF

echo "[OK] ${REFERENCE_DIR}/grch37_to_grch38.chain"
echo "[OK] ${REFERENCE_DIR}/GRCh38.fa"
