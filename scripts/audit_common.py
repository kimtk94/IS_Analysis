#!/usr/bin/env python3
"""Shared helpers for IS_Analysis audit scripts."""
from __future__ import annotations

import csv
import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

PROJECT_ROOT = Path(__file__).resolve().parents[1]
QC_DIR = PROJECT_ROOT / "results" / "qc"
REQUIRED_DIRS = ["scripts", "R", "config", "results/qc", "Backup", "scripts/Backup"]
MINIMUM_REVIEW_FILES = [
    "scripts/00_run_full_audit.py",
    "scripts/audit_common.py",
    "scripts/01_build_expected_manifest.py",
    "scripts/02_audit_raw_downloads.py",
    "scripts/01_prepare_exposure_fast.R",
    "scripts/patch_empty_exposure_outputs.py",
    "scripts/patch_empty_output_current.py",
    "scripts/patch_runner_output_exists.py",
    "scripts/Backup/01_prepare_exposure.R",
]
QC_REQUIRED_GLOBS = [
    "*current*manifest*.tsv",
    "*audit*summary*.tsv",
    "*tracking*manifest*.tsv",
]


@dataclass(frozen=True)
class FileRecord:
    path: str
    exists: bool
    size_bytes: int
    sha256: str


def project_path(path: str | Path) -> Path:
    return PROJECT_ROOT / path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: str | Path) -> FileRecord:
    rel = Path(path)
    absolute = project_path(rel)
    if not absolute.exists() or not absolute.is_file():
        return FileRecord(str(rel), False, 0, "")
    return FileRecord(str(rel), True, absolute.stat().st_size, sha256_file(absolute))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))
