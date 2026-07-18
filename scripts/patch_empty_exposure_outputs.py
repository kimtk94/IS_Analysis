#!/usr/bin/env python3
"""Patch exposure preparation scripts to write schema-preserving empty TSVs.

This legacy recovery helper is retained for reproducibility of the 2026-07-17
repair workflow. The active implementation in ``scripts/01_prepare_exposure_fast.R``
already includes ``empty_exposure_dt()`` and ``write_empty_exposure()``.
"""
from __future__ import annotations

from pathlib import Path

CANDIDATES = [
    Path("01_prepare_exposure_fast.R"),
    Path("scripts/01_prepare_exposure_fast.R"),
    Path("scripts/01_prepare_exposure.R"),
]

script = next((path for path in CANDIDATES if path.exists()), None)
if script is None:
    raise SystemExit("[ERROR] Cannot find 01_prepare_exposure_fast.R or scripts/01_prepare_exposure.R")

text = script.read_text(encoding="utf-8")

helper = r'''
# -----------------------------
# Empty output schema helpers
# -----------------------------

empty_exposure_dt <- function() {
  data.table(
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
}

write_empty_exposure <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(empty_exposure_dt(), path, sep = "\t")
}
'''

if "empty_exposure_dt <- function()" not in text:
    marker = "# -----------------------------\n# Column helpers"
    if marker in text:
        text = text.replace(marker, helper + "\n\n" + marker, 1)
    else:
        text = helper + "\n\n" + text

text = text.replace('fwrite(data.table(), gene_out, sep = "\\t")', 'write_empty_exposure(gene_out)')
text = text.replace('fwrite(data.table(), batch_out, sep = "\\t")', 'write_empty_exposure(batch_out)')
text = text.replace('fwrite(data.table(), out, sep = "\\t")', 'write_empty_exposure(out)')

backup = script.with_suffix(script.suffix + ".bak")
backup.write_text(script.read_text(encoding="utf-8"), encoding="utf-8")
script.write_text(text, encoding="utf-8")

print("[OK] patched:", script)
print("[OK] backup:", backup)
