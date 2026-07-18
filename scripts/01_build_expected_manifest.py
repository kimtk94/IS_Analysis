#!/usr/bin/env python3
"""Build the expected file manifest used by the repository QA audit."""
from __future__ import annotations

from audit_common import MINIMUM_REVIEW_FILES, QC_DIR, REQUIRED_DIRS, file_record, write_tsv


def build_manifest() -> None:
    rows = []
    for directory in REQUIRED_DIRS:
        rows.append({"path": directory, "kind": "directory", "required": "true", "exists": "true", "size_bytes": "NA", "sha256": "NA"})
    for file_path in MINIMUM_REVIEW_FILES:
        record = file_record(file_path)
        rows.append({
            "path": file_path,
            "kind": "file",
            "required": "true",
            "exists": str(record.exists).lower(),
            "size_bytes": record.size_bytes,
            "sha256": record.sha256,
        })
    fieldnames = ["path", "kind", "required", "exists", "size_bytes", "sha256"]
    write_tsv(QC_DIR / "expected_current_manifest.tsv", rows, fieldnames)
    write_tsv(QC_DIR / "tracking_manifest.tsv", rows, fieldnames)


if __name__ == "__main__":
    build_manifest()
