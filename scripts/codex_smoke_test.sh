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
bash -n scripts/colab_download_grch38_liftover_references.sh
"${PYTHON_BIN}" -m py_compile \
  scripts/00_run_full_audit_final.py \
  scripts/build_ukb_ppp_download_manifest.py \
  scripts/synapse_metadata.py \
  scripts/ukb_ppp_batch_manifest_runner_fast.py \
  scripts/colab_download_gigastroke_gwas.py \
  scripts/configure_gigastroke_outcomes.py \
  workflow/gigastroke_outcome_adapter.py \
  workflow/causal_checkpoint_analysis.py

echo "[TEST] Independent causal-analysis checkpoints"
CAUSAL_OUT="/tmp/is_analysis_causal_checkpoint"
rm -rf "${CAUSAL_OUT}"
for stage in harmonization mr sensitivity colocalization; do
  "${PYTHON_BIN}" workflow/causal_checkpoint_analysis.py \
    --config tests/fixtures/causal_checkpoint/config.json --stage "${stage}"
done
"${PYTHON_BIN}" - "${CAUSAL_OUT}" <<'PY'
import csv
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
with (root / "harmonization/variants.tsv").open(newline="", encoding="utf-8") as handle:
    harmonized = {row["variant_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
assert harmonized["rs-direct"]["allele_state"] == "DIRECT"
assert harmonized["rs-swapped"]["allele_state"] == "SWAPPED"
assert harmonized["rs-complement"]["allele_state"] == "COMPLEMENT"
assert harmonized["rs-palindromic"]["reason"] == "PALINDROMIC_AMBIGUOUS"
assert harmonized["rs-mismatch"]["allele_state"] == "MISMATCH"
with (root / "mr/estimates.tsv").open(newline="", encoding="utf-8") as handle:
    estimates = {row["method"]: row for row in csv.DictReader(handle, delimiter="\t")}
assert estimates["IVW"]["status"] == "SUCCESS"
assert all(estimates[name]["status"] == "SUCCESS" for name in ("WEIGHTED_MEDIAN", "MR_EGGER", "WEIGHTED_MODE"))
with (root / "sensitivity/results.tsv").open(newline="", encoding="utf-8") as handle:
    sensitivity = list(csv.DictReader(handle, delimiter="\t"))
assert {row["analysis"] for row in sensitivity} >= {"COCHRAN_Q", "EGGER_INTERCEPT", "LEAVE_ONE_OUT", "STEIGER_DIRECTIONALITY", "REVERSE_MR"}
with (root / "colocalization/results.tsv").open(newline="", encoding="utf-8") as handle:
    coloc = list(csv.DictReader(handle, delimiter="\t"))
assert len(coloc) == 2 and all(row["common_snps"] == "5" for row in coloc)
assert json.loads((root / "colocalization/multiple_signal_status.json").read_text())["status"] == "READY_SUSIE"
assert all(json.loads((root / stage / "status.json").read_text())["status"] == "SUCCESS" for stage in ("harmonization", "mr", "sensitivity", "colocalization"))
print("[OK] allele states, MR selection, conditional sensitivity, regional coloc, and checkpoint statuses")
PY

echo "[TEST] Single-instrument Wald and reasoned insufficient-instrument states"
"${PYTHON_BIN}" - <<'PY'
import json
from pathlib import Path

source = Path("tests/fixtures/causal_checkpoint/config.json")
config = json.loads(source.read_text())
config.update(instrument_exposure="exposure_single.tsv", instrument_outcome="outcome_single.tsv",
              output_dir="/tmp/is_analysis_causal_single")
Path("tests/fixtures/causal_checkpoint/config.single.runtime.json").write_text(json.dumps(config))
PY
rm -rf /tmp/is_analysis_causal_single
"${PYTHON_BIN}" workflow/causal_checkpoint_analysis.py --config tests/fixtures/causal_checkpoint/config.single.runtime.json --stage harmonization
"${PYTHON_BIN}" workflow/causal_checkpoint_analysis.py --config tests/fixtures/causal_checkpoint/config.single.runtime.json --stage mr
"${PYTHON_BIN}" workflow/causal_checkpoint_analysis.py --config tests/fixtures/causal_checkpoint/config.single.runtime.json --stage sensitivity
rm tests/fixtures/causal_checkpoint/config.single.runtime.json
"${PYTHON_BIN}" - <<'PY'
import csv
from pathlib import Path

root = Path("/tmp/is_analysis_causal_single")
with (root / "mr/estimates.tsv").open(newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
assert rows[0]["method"] == "WALD_RATIO" and rows[0]["status"] == "SUCCESS"
assert all(row["status"] == "NOT_RUN_INSUFFICIENT_INSTRUMENTS" for row in rows[1:])
with (root / "sensitivity/results.tsv").open(newline="") as handle:
    sensitivity = list(csv.DictReader(handle, delimiter="\t"))
assert next(row for row in sensitivity if row["analysis"] == "COCHRAN_Q")["status"] == "NOT_RUN_INSUFFICIENT_INSTRUMENTS"
print("[OK] single SNP uses Wald ratio and unsupported analyses retain reasoned states")
PY

echo "[TEST] GIGASTROKE outcome adapter fixture"
GIGASTROKE_OUT="${SMOKE_ROOT}/gigastroke"
"${PYTHON_BIN}" workflow/gigastroke_outcome_adapter.py \
  --config tests/fixtures/gigastroke/config.json --output-dir "${GIGASTROKE_OUT}"
"${PYTHON_BIN}" - "${GIGASTROKE_OUT}" <<'PY'
import csv
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
with (root / "gigastroke_is_EUR.canonical.tsv").open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
assert len(rows) == 3
assert rows[0]["source_chromosome"] == "1" and rows[0]["chromosome"] == "1"
assert rows[0]["source_variant_id"] == "rs-good"
derived = next(row for row in rows if row["source_position"] == "6")
assert derived["source_variant_id"] == ""
assert (derived["source_ref"], derived["source_alt"], derived["ref"], derived["alt"]) == ("A", "T", "A", "T")
indel = next(row for row in rows if row["source_variant_id"] == "rs-indel")
assert (indel["source_position"], indel["source_ref"], indel["source_alt"]) == ("11", "AA", "A")
assert (indel["position"], indel["ref"], indel["alt"]) == ("1", "AA", "A")
with (root / "gigastroke_is_EUR.rejected.tsv").open(newline="", encoding="utf-8") as handle:
    rejects = list(csv.DictReader(handle, delimiter="\t"))
assert {row["reason"] for row in rejects} == {"unmapped", "multi_mapped", "duplicate", "reference_allele_mismatch", "source_reference_allele_mismatch"}
manifest = json.loads((root / "dataset_manifest.json").read_text(encoding="utf-8"))
assert {(item["ancestry"], item["role"], item["source_build"], item["target_build"]) for item in manifest} == {
    ("EUR", "discovery", "GRCh37", "GRCh38"), ("EAS", "replication_subtype", "GRCh37", "GRCh38")}
assert all(len(item["chain_sha256"]) == 64 and len(item["source_reference_sha256"]) == 64 and len(item["target_reference_sha256"]) == 64 for item in manifest)
print("[OK] liftover, normalization, provenance, selection, and reasoned exclusions")
PY

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

member = Path("tests/fixtures/pqtl_locus_member.tsv").read_bytes()
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
  --gene-coordinate-file tests/fixtures/gene_coordinates_grch38.tsv \
  --standardized-dir "${SMOKE_ROOT}/standardized/pqtl" \
  --instrument-dir "${SMOKE_ROOT}/instrument_candidates"
"${PYTHON_BIN}" - "${DUPLICATE_OUT_DIR}/exposure_duplicate_gene.tsv" <<'PY'
import csv
import sys

with open(sys.argv[1], newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
assert len(rows) == 6, f"Expected all three significant cis SNPs from two archives, found {len(rows)}"
assert {row["source_file"] for row in rows} == {
    "ALPHA_P12345_OID1_v1_PANEL.tar",
    "ALPHA_Q99999_OID2_v1_PANEL.tar",
}
exposures = {row["id.exposure"] for row in rows}
assert len(exposures) == 2, f"Expected distinct exposure IDs per archive, found {exposures}"
assert all(row["gene_symbol"] == "ALPHA" and row["ancestry"] == "EUR" for row in rows)
print("[OK] Both same-gene archives retained with distinct exposure IDs")
PY

echo "[TEST] Canonical full-summary retains every valid locus SNP"
"${PYTHON_BIN}" - "${SMOKE_ROOT}/standardized/pqtl/EUR/duplicate_gene/ALPHA__ALPHA_P12345_OID1_v1_PANEL.tsv" <<'PY'
import csv
import sys

with open(sys.argv[1], newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
required = {"dataset_id", "protein_id", "assay_id", "gene_symbol", "ancestry", "genome_build",
            "chromosome", "position", "rsid", "ref", "alt", "effect_allele", "other_allele",
            "beta", "se", "p_value", "eaf", "sample_size", "source_archive"}
assert required == set(rows[0]), set(rows[0])
assert len(rows) == 4, f"Full summary must retain cis and trans locus rows, found {len(rows)}"
assert {row["rsid"] for row in rows} == {"rs101", "rs102", "rs103", "rs104"}
assert all(row["dataset_id"] == "UKB-PPP" and row["genome_build"] == "GRCh38" for row in rows)
print("[OK] canonical full-summary schema and complete locus retention")
PY

echo "[TEST] Missing gene coordinate fails cis instrument selection"
MISSING_OUT="${SMOKE_ROOT}/missing_coordinate_output"
if Rscript scripts/01_prepare_exposure_fast.R \
  --gene-file "${DUPLICATE_GENE_FILE}" --batch-id missing_coordinate \
  --outdir "${MISSING_OUT}" --rawdir "${DUPLICATE_RAW_DIR}" --ancestries EUR \
  --gene-coordinate-file tests/fixtures/gene_coordinates_missing.tsv \
  --standardized-dir "${SMOKE_ROOT}/missing_standardized" \
  --instrument-dir "${SMOKE_ROOT}/missing_instruments"; then
  echo "[ERROR] Missing coordinate unexpectedly produced valid instruments" >&2
  exit 1
fi
test ! -s "${SMOKE_ROOT}/missing_instruments/EUR/exposure_missing_coordinate.tsv" || \
  test "$(wc -l < "${SMOKE_ROOT}/missing_instruments/EUR/exposure_missing_coordinate.tsv")" -eq 1
echo "[OK] Missing coordinate rejected; no valid instrument rows emitted"

echo "[OK] Codex smoke test completed"
