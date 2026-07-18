#!/usr/bin/env python3
"""Validate restored project structure and QC manifest coverage."""
from __future__ import annotations

from pathlib import Path

from audit_common import MINIMUM_REVIEW_FILES, PROJECT_ROOT, QC_DIR, QC_REQUIRED_GLOBS, REQUIRED_DIRS, file_record, write_tsv


def audit() -> int:
    rows = []
    failures = 0
    for directory in REQUIRED_DIRS:
        exists = (PROJECT_ROOT / directory).is_dir()
        failures += 0 if exists else 1
        rows.append({"check": "required_directory", "target": directory, "status": "pass" if exists else "fail", "detail": "directory present" if exists else "missing directory"})
    for file_path in MINIMUM_REVIEW_FILES:
        record = file_record(file_path)
        ok = record.exists and record.size_bytes > 0
        failures += 0 if ok else 1
        rows.append({"check": "minimum_review_file", "target": file_path, "status": "pass" if ok else "fail", "detail": f"{record.size_bytes} bytes" if record.exists else "missing file"})
    for pattern in QC_REQUIRED_GLOBS:
        matches = sorted(str(path.relative_to(PROJECT_ROOT)) for path in QC_DIR.glob(pattern))
        ok = bool(matches)
        failures += 0 if ok else 1
        rows.append({"check": "qc_tsv_family", "target": pattern, "status": "pass" if ok else "fail", "detail": ",".join(matches) if matches else "no matching TSV"})
    write_tsv(QC_DIR / "audit_summary.tsv", rows, ["check", "target", "status", "detail"])
    return failures


if __name__ == "__main__":
    raise SystemExit(1 if audit() else 0)
