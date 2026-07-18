#!/usr/bin/env python3
"""Patch the current fast exposure script to replace empty data.table writes.

This legacy helper mirrors the 2026-07-17 empty-output patch and is kept for
traceability. Prefer editing ``scripts/01_prepare_exposure_fast.R`` directly for
new development.
"""
from __future__ import annotations

from pathlib import Path
import re

CANDIDATES = [Path("01_prepare_exposure_fast.R"), Path("scripts/01_prepare_exposure_fast.R")]
r_path = next((path for path in CANDIDATES if path.exists()), None)
if r_path is None:
    raise SystemExit("[ERROR] 01_prepare_exposure_fast.R not found")

text = r_path.read_text(encoding="utf-8")
backup = r_path.with_suffix(r_path.suffix + ".bak_empty_output_patch")
backup.write_text(text, encoding="utf-8")

helper = r'''
write_empty_exposure <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  empty_dt <- data.table(
    gene_symbol = character(),
    ancestry = character(),
    source_file = character(),
    chr = character(),
    pos_hg38 = integer(),
    effect_allele = character(),
    other_allele = character(),
    beta = numeric(),
    se = numeric(),
    SNP = character(),
    pval = numeric(),
    samplesize = numeric(),
    eaf = numeric(),
    exposure = character(),
    id.exposure = character(),
    beta.exposure = numeric(),
    se.exposure = numeric(),
    effect_allele.exposure = character(),
    other_allele.exposure = character(),
    eaf.exposure = numeric(),
    pval.exposure = numeric(),
    samplesize.exposure = numeric(),
    is_cis = logical(),
    F_stat = numeric()
  )

  fwrite(empty_dt, path, sep = "\t")
}
'''

if "write_empty_exposure <- function" not in text:
    if "parse_variant_id <- function" in text:
        text = text.replace("parse_variant_id <- function", helper + "\n\nparse_variant_id <- function", 1)
    else:
        text = helper + "\n\n" + text

text = re.sub(
    r'fwrite\s*\(\s*data\.table\s*\(\s*\)\s*,\s*([A-Za-z_][A-Za-z0-9_]*)\s*,\s*sep\s*=\s*["\']\\\\t["\']\s*\)',
    r'write_empty_exposure(\1)',
    text,
)
text = re.sub(
    r'fwrite\s*\(\s*data\.table\s*\(\s*\)\s*,\s*([A-Za-z_][A-Za-z0-9_]*)\s*,\s*sep\s*=\s*["\']\\t["\']\s*\)',
    r'write_empty_exposure(\1)',
    text,
)

r_path.write_text(text, encoding="utf-8")

print("[OK] patched:", r_path)
print("[OK] backup:", backup)

remain = [line for line in text.splitlines() if "fwrite(data.table()" in line]
if remain:
    print("[WARN] Remaining fwrite(data.table()) lines:")
    for line in remain:
        print(line)
else:
    print("[OK] No remaining fwrite(data.table()) lines")
