#!/usr/bin/env python3
"""Self-contained integrity audit for the legacy IS_Analysis pipeline.

This file consolidates the former audit_common.py and 01-06 audit scripts.
It does not execute MR or modify raw/results data. It writes audit evidence only.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import tarfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

pd: Any = None


TRUE_TEXT = {"TRUE", "T", "1", "YES", "Y"}
DELETED_PIPELINE_STATES = {
    "deleted_after_processed",
    "deleted_after_qc",
    "released",
}
DELETED_RETENTION_STATES = {
    "deleted_after_processed",
    "deleted_after_qc",
}
ANCESTRY_ALIASES = {
    "EUR": "EUR",
    "EUROPEAN": "EUR",
    "EUROPEAN_DISCOVERY": "EUR",
    "EAS": "EAS",
    "EAST_ASIAN": "EAS",
    "EASTASIAN": "EAS",
    "AFR": "AFR",
    "AFRICAN": "AFR",
    "CSA": "CSA",
    "CENTRAL_SOUTH_ASIAN": "CSA",
    "MID": "MID",
    "MIDDLE_EAST": "MID",
    "AMR": "AMR",
    "AMERICAN": "AMR",
    "COMBINED": "COMBINED",
}
SOURCE_COLUMNS = ["source_file", "source_file_name", "file_name", "archive"]


def ensure_pandas() -> Any:
    try:
        import pandas
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "[ERROR] Missing Python dependency: pandas. "
            "Run `bash scripts/setup_codex_env.sh` or `python -m pip install -r requirements.txt`."
        ) from exc
    return pandas


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def stamp_now() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_bool(value: object) -> bool:
    return str(value).strip().upper() in TRUE_TEXT


def normalize_source_file(value: object) -> str:
    return Path(str(value).strip()).name


def normalize_ancestry(value: object) -> str:
    text = str(value).strip().upper().replace(" ", "_").replace("-", "_")
    return ANCESTRY_ALIASES.get(text, text)


def root_path(config: dict[str, Any], value: object) -> Path:
    path = Path(str(value))
    return path if path.is_absolute() else config["_root"] / path


def load_config(config_path: str | None) -> dict[str, Any]:
    scripts_dir = Path(__file__).resolve().parent
    path = Path(config_path) if config_path else scripts_dir.parent / "config" / "audit_config.json"
    if not path.exists():
        raise FileNotFoundError(f"Audit config not found: {path}")
    config = json.loads(path.read_text(encoding="utf-8"))
    config["_config_path"] = path.resolve()
    config["_root"] = Path(config["project_root"]).resolve()
    return config


def read_table(path: Path, **kwargs: Any) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path, dtype=str, **kwargs).fillna("")
    if suffix == ".csv":
        return pd.read_csv(path, dtype=str, **kwargs).fillna("")
    return pd.read_csv(path, sep="\t", dtype=str, **kwargs).fillna("")


def atomic_write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    frame.to_csv(temporary, sep="\t", index=False)
    temporary.replace(path)


def write_current_and_snapshot(
    frame: pd.DataFrame,
    output_dir: Path,
    stem: str,
    stamp: str,
) -> tuple[Path, Path]:
    current = output_dir / f"{stem}.current.tsv"
    snapshot = output_dir / f"{stem}_{stamp}.tsv"
    atomic_write_tsv(frame, current)
    atomic_write_tsv(frame, snapshot)
    return current, snapshot


def first_existing_column(frame: pd.DataFrame, candidates: Iterable[str]) -> str | None:
    mapping = {str(column).lower(): str(column) for column in frame.columns}
    for candidate in candidates:
        if candidate.lower() in mapping:
            return mapping[candidate.lower()]
    return None


def parse_source_filename(source_file: object) -> dict[str, str]:
    name = normalize_source_file(source_file)
    stem = re.sub(r"\.(tar\.gz|tgz|tar|gz|zip)$", "", name, flags=re.I)
    parts = stem.split("_")
    result = {
        "filename_gene": parts[0] if parts else "",
        "filename_uniprot": "",
        "filename_oid": "",
        "filename_version": "",
        "filename_panel": "",
    }
    version_index: int | None = None
    for index, token in enumerate(parts):
        if (
            not result["filename_uniprot"]
            and re.fullmatch(r"[A-NR-Z0-9][A-Z0-9]{5,9}(?:-\d+)?", token, re.I)
            and not token.upper().startswith("OID")
        ):
            result["filename_uniprot"] = token
        if re.fullmatch(r"OID\d+", token, re.I):
            result["filename_oid"] = token
        if re.fullmatch(r"v\d+", token, re.I):
            result["filename_version"] = token
            version_index = index
    if version_index is not None and version_index + 1 < len(parts):
        result["filename_panel"] = "_".join(parts[version_index + 1 :])
    return result


def choose_manifest(config: dict[str, Any]) -> pd.DataFrame:
    candidates = [
        root_path(config, config["existing_master_manifest"]),
        root_path(config, config["tracking_manifest"]),
    ]
    for path in candidates:
        frame = read_table(path)
        if frame.empty:
            continue
        source_column = first_existing_column(frame, SOURCE_COLUMNS)
        if source_column:
            if source_column != "source_file":
                frame = frame.rename(columns={source_column: "source_file"})
            frame["manifest_origin"] = str(path)
            return frame
    raise FileNotFoundError("No non-empty manifest with a source-file column was found.")


def build_expected_manifest(config: dict[str, Any]) -> pd.DataFrame:
    expected = choose_manifest(config).copy()
    for column in ["batch_id", "dataset_id", "ancestry"]:
        if column not in expected.columns:
            expected[column] = ""
    expected["source_file"] = expected["source_file"].map(normalize_source_file)
    expected["ancestry"] = expected["ancestry"].map(normalize_ancestry)
    parsed = pd.DataFrame(
        expected["source_file"].map(parse_source_filename).tolist(),
        index=expected.index,
    )
    for column in parsed.columns:
        expected[column] = parsed[column]
    expected["unit_key"] = (
        expected["dataset_id"].astype(str).str.strip()
        + "||"
        + expected["ancestry"].astype(str).str.strip()
        + "||"
        + expected["source_file"]
    )
    expected["batch_source_key"] = (
        expected["batch_id"].astype(str).str.strip()
        + "||"
        + expected["ancestry"].astype(str).str.strip()
        + "||"
        + expected["source_file"]
    )
    expected["duplicate_unit_key_count"] = expected.groupby("unit_key")[
        "unit_key"
    ].transform("size")
    expected["expected_manifest_status"] = "PASS"
    expected.loc[
        expected["source_file"].eq(""), "expected_manifest_status"
    ] = "FAILED_EMPTY_SOURCE"
    expected.loc[
        expected["duplicate_unit_key_count"].astype(int).gt(1),
        "expected_manifest_status",
    ] = "FAILED_DUPLICATE_UNIT"
    return expected


def declared_deleted(row: pd.Series) -> bool:
    if normalize_bool(row.get("do_not_redownload", "")):
        return True
    if str(row.get("pipeline_status", "")).strip().lower() in DELETED_PIPELINE_STATES:
        return True
    return str(row.get("raw_retention_status", "")).strip().lower() in DELETED_RETENTION_STATES


def raw_candidates(config: dict[str, Any], row: pd.Series) -> list[Path]:
    filename = normalize_source_file(row["source_file"])
    ancestry = normalize_ancestry(row.get("ancestry", ""))
    candidates: list[Path] = []
    local_path = str(row.get("local_path", "")).strip()
    if local_path:
        candidate = Path(local_path)
        candidates.append(candidate if candidate.is_absolute() else config["_root"] / candidate)
    for raw_dir in config.get("raw_dirs", []):
        base = root_path(config, raw_dir)
        if ancestry:
            candidates.append(base / ancestry / filename)
        candidates.append(base / filename)
    unique: list[Path] = []
    seen: set[str] = set()
    for candidate in candidates:
        key = str(candidate)
        if key not in seen:
            unique.append(candidate)
            seen.add(key)
    return unique


def inspect_tar(path: Path) -> dict[str, object]:
    result: dict[str, object] = {
        "tar_check_status": "NOT_CHECKED",
        "tar_member_count": "",
        "tar_first_member": "",
        "tar_error": "",
    }
    try:
        with tarfile.open(path, "r:*") as archive:
            members = archive.getmembers()
        result["tar_check_status"] = "PASS" if members else "EMPTY"
        result["tar_member_count"] = len(members)
        result["tar_first_member"] = members[0].name if members else ""
    except Exception as exc:
        result["tar_check_status"] = "FAIL"
        result["tar_error"] = f"{type(exc).__name__}: {exc}"
    return result


def audit_raw_downloads(config: dict[str, Any], expected: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    minimum_size = int(config.get("minimum_raw_size_bytes", 1))
    for _, source in expected.iterrows():
        existing = next((path for path in raw_candidates(config, source) if path.exists()), None)
        row: dict[str, object] = {
            "batch_source_key": source["batch_source_key"],
            "source_file": source["source_file"],
            "ancestry": source.get("ancestry", ""),
            "raw_path": "",
            "raw_exists": existing is not None,
            "raw_size_bytes": 0,
            "raw_size_pass": False,
            "source_sha256": "",
            "declared_deleted": declared_deleted(source),
            "tar_check_status": "NOT_AVAILABLE",
            "tar_member_count": "",
            "tar_first_member": "",
            "tar_error": "",
        }
        if existing:
            size = existing.stat().st_size
            row.update(
                raw_path=str(existing),
                raw_size_bytes=size,
                raw_size_pass=size >= minimum_size,
            )
            if config.get("calculate_sha256", False):
                row["source_sha256"] = sha256_file(existing)
            if config.get("check_tar_integrity", False):
                row.update(inspect_tar(existing))
        if not row["raw_exists"]:
            row["raw_audit_status"] = (
                "RAW_MISSING_DECLARED_DELETED"
                if row["declared_deleted"]
                else "RAW_MISSING_UNRECOVERED"
            )
        elif not row["raw_size_pass"] or row["tar_check_status"] in {"FAIL", "EMPTY"}:
            row["raw_audit_status"] = "RAW_PRESENT_INVALID"
        else:
            row["raw_audit_status"] = "RAW_PRESENT"
        rows.append(row)
    return pd.DataFrame(rows)


SITE_COLUMN_CANDIDATES = {
    "source_file": SOURCE_COLUMNS,
    "gene": ["gene", "symbol", "protein", "protein_name"],
    "uniprot": ["uniprot", "uniprot_id", "accession"],
    "oid": ["oid", "olink_id", "assay_id"],
    "version": ["version", "assay_version"],
    "panel": ["panel", "panel_name", "product"],
}


def normalized_metadata(value: object) -> str:
    return str(value).strip().lower().replace(" ", "_").replace("-", "_")


def audit_site_metadata(config: dict[str, Any], expected: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "batch_source_key",
        "source_file",
        "filename_gene",
        "filename_uniprot",
        "filename_oid",
        "filename_version",
        "filename_panel",
    ]
    base = expected[columns].copy()
    site_path = root_path(config, config["site_metadata_file"])
    site = read_table(site_path)
    if site.empty:
        base["site_metadata_available"] = False
        base["metadata_match_status"] = "NOT_AVAILABLE"
        base["metadata_mismatch_fields"] = ""
        base["site_metadata_file"] = str(site_path)
        return base

    configured = config.get("site_metadata_column_map", {})
    resolved = {
        name: configured.get(name) or first_existing_column(site, candidates)
        for name, candidates in SITE_COLUMN_CANDIDATES.items()
    }
    if not resolved["source_file"]:
        raise ValueError("Site metadata needs a source-file column mapping.")
    normalized = pd.DataFrame(
        {"source_file": site[resolved["source_file"]].map(normalize_source_file)}
    )
    for field in ["gene", "uniprot", "oid", "version", "panel"]:
        normalized[f"site_{field}"] = (
            site[resolved[field]].astype(str).str.strip() if resolved[field] else ""
        )
    merged = base.merge(
        normalized.drop_duplicates("source_file", keep="last"),
        on="source_file",
        how="left",
    )
    statuses: list[str] = []
    mismatches: list[str] = []
    for _, row in merged.iterrows():
        if not str(row.get("site_gene", "")).strip():
            statuses.append("NOT_FOUND")
            mismatches.append("source_file")
            continue
        mismatch: list[str] = []
        comparisons = [
            ("gene", str(row.filename_gene).upper(), str(row.site_gene).upper()),
            ("uniprot", str(row.filename_uniprot).upper(), str(row.site_uniprot).upper()),
            ("oid", str(row.filename_oid).upper(), str(row.site_oid).upper()),
            ("version", str(row.filename_version).lower(), str(row.site_version).lower()),
            ("panel", normalized_metadata(row.filename_panel), normalized_metadata(row.site_panel)),
        ]
        for field, left, right in comparisons:
            if left and right and left != right:
                mismatch.append(field)
        statuses.append("PASS" if not mismatch else "MISMATCH")
        mismatches.append("|".join(mismatch))
    merged["site_metadata_available"] = True
    merged["metadata_match_status"] = statuses
    merged["metadata_mismatch_fields"] = mismatches
    merged["site_metadata_file"] = str(site_path)
    return merged


def count_output_rows_and_sources(
    path: Path,
    required_groups: list[list[str]],
    scan_rows: bool,
) -> dict[str, object]:
    result: dict[str, object] = {
        "output_path": str(path),
        "output_exists": path.exists(),
        "output_size_bytes": path.stat().st_size if path.exists() else 0,
        "output_status": "NOT_RUN",
        "output_n_rows": 0,
        "output_n_cols": 0,
        "output_columns": "",
        "missing_required_groups": "",
        "source_column": "",
        "observed_sources": set(),
        "output_error": "",
    }
    if not path.exists():
        return result
    try:
        header = pd.read_csv(path, sep="\t", nrows=0)
        columns = [str(column) for column in header.columns]
        lowercase = {column.lower() for column in columns}
        missing = [
            "/".join(group)
            for group in required_groups
            if not any(str(alias).lower() in lowercase for alias in group)
        ]
        source_column = first_existing_column(header, SOURCE_COLUMNS)
        result.update(
            output_n_cols=len(columns),
            output_columns="|".join(columns),
            missing_required_groups="|".join(missing),
            source_column=source_column or "",
        )
        row_count = 0
        observed: set[str] = set()
        if scan_rows:
            usecols = [source_column] if source_column else None
            for chunk in pd.read_csv(
                path,
                sep="\t",
                dtype=str,
                usecols=usecols,
                chunksize=200_000,
            ):
                row_count += len(chunk)
                if source_column:
                    observed.update(
                        chunk[source_column].dropna().map(normalize_source_file).tolist()
                    )
        result["output_n_rows"] = row_count
        result["observed_sources"] = observed
        if missing or not source_column:
            result["output_status"] = "FAILED_SCHEMA"
        elif row_count == 0:
            result["output_status"] = "FAILED_EMPTY_UNEXPLAINED"
        else:
            result["output_status"] = "SUCCESS_NONEMPTY"
    except Exception as exc:
        result["output_status"] = "FAILED_RUNTIME"
        result["output_error"] = f"{type(exc).__name__}: {exc}"
    return result


def output_path(config: dict[str, Any], directory_key: str, batch_id: str) -> Path:
    return root_path(config, config[directory_key]) / config["batch_output_pattern"].format(
        batch_id=batch_id
    )


def audit_prepare_outputs(
    config: dict[str, Any],
    expected: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    source_rows: list[dict[str, object]] = []
    batch_rows: list[dict[str, object]] = []
    required_groups = config["required_exposure_columns_any"]
    scan_rows = bool(config.get("read_output_rows", True))
    batch_ids = sorted(set(expected["batch_id"].astype(str).str.strip()) - {""})
    for batch_id in batch_ids:
        batch_expected = expected[
            expected["batch_id"].astype(str).str.strip().eq(batch_id)
        ]
        original = count_output_rows_and_sources(
            output_path(config, "original_exposure_dir", batch_id),
            required_groups,
            scan_rows,
        )
        recovery = count_output_rows_and_sources(
            output_path(config, "recovery_exposure_dir", batch_id),
            required_groups,
            scan_rows,
        )
        original_sources = original.pop("observed_sources")
        recovery_sources = recovery.pop("observed_sources")
        batch_rows.append(
            {
                "batch_id": batch_id,
                "expected_source_count": len(batch_expected),
                "original_observed_source_count": len(original_sources),
                "recovery_observed_source_count": len(recovery_sources),
                **{f"original_{key}": value for key, value in original.items()},
                **{f"recovery_{key}": value for key, value in recovery.items()},
            }
        )
        for _, source in batch_expected.iterrows():
            filename = normalize_source_file(source["source_file"])
            in_original = filename in original_sources
            in_recovery = filename in recovery_sources
            complete = in_original or in_recovery
            source_rows.append(
                {
                    "batch_source_key": source["batch_source_key"],
                    "batch_id": batch_id,
                    "source_file": filename,
                    "complete_in_original_output": in_original,
                    "complete_in_recovery_output": in_recovery,
                    "prepare_complete": complete,
                    "prepare_completion_source": (
                        "ORIGINAL" if in_original else "RECOVERY" if in_recovery else "NONE"
                    ),
                    "prepare_qc_status": "SUCCESS_NONEMPTY" if complete else "REVIEW",
                    "original_output_status": original["output_status"],
                    "recovery_output_status": recovery["output_status"],
                }
            )
    return pd.DataFrame(source_rows), pd.DataFrame(batch_rows)


def merge_audit(
    base: pd.DataFrame,
    audit: pd.DataFrame,
    columns: list[str],
) -> pd.DataFrame:
    if audit.empty:
        return base
    selected = ["batch_source_key"] + [column for column in columns if column in audit.columns]
    return base.merge(
        audit[selected].drop_duplicates("batch_source_key"),
        on="batch_source_key",
        how="left",
    )


def build_master_manifest(
    expected: pd.DataFrame,
    raw: pd.DataFrame,
    metadata: pd.DataFrame,
    prepare: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    master = merge_audit(
        expected,
        raw,
        [
            "raw_path",
            "raw_exists",
            "raw_size_bytes",
            "raw_size_pass",
            "source_sha256",
            "declared_deleted",
            "tar_check_status",
            "raw_audit_status",
        ],
    )
    master = merge_audit(
        master,
        metadata,
        [
            "site_metadata_available",
            "metadata_mismatch_fields",
            "metadata_match_status",
        ],
    )
    master = merge_audit(
        master,
        prepare,
        [
            "complete_in_original_output",
            "complete_in_recovery_output",
            "prepare_complete",
            "prepare_completion_source",
            "prepare_qc_status",
            "original_output_status",
            "recovery_output_status",
        ],
    )
    prepare_ok = master["prepare_complete"].fillna(False).map(normalize_bool)
    raw_present = master["raw_audit_status"].eq("RAW_PRESENT")
    deleted_verified = (
        master["raw_audit_status"].eq("RAW_MISSING_DECLARED_DELETED") & prepare_ok
    )
    master["raw_final_state"] = master["raw_audit_status"]
    master.loc[deleted_verified, "raw_final_state"] = "RAW_DELETED_PROCESSED_VERIFIED"
    raw_ok = raw_present | deleted_verified
    expected_ok = master["expected_manifest_status"].eq("PASS")
    metadata_ok = master["metadata_match_status"].isin(["PASS", "NOT_AVAILABLE"])
    master["final_status"] = "REVIEW"
    master.loc[expected_ok & raw_ok & metadata_ok & prepare_ok, "final_status"] = "READY"
    master["recommended_action"] = "review"
    master.loc[~expected_ok, "recommended_action"] = "fix_manifest"
    master.loc[expected_ok & ~raw_ok, "recommended_action"] = "redownload_or_check_raw"
    master.loc[expected_ok & raw_ok & ~metadata_ok, "recommended_action"] = "review_site_metadata"
    master.loc[
        expected_ok & raw_ok & metadata_ok & ~prepare_ok,
        "recommended_action",
    ] = "reprepare_or_review_log"
    master.loc[master["final_status"].eq("READY"), "recommended_action"] = "proceed_next_step"
    reasons: list[str] = []
    for _, row in master.iterrows():
        reason: list[str] = []
        if row["expected_manifest_status"] != "PASS":
            reason.append(str(row["expected_manifest_status"]))
        if row["raw_final_state"] not in {"RAW_PRESENT", "RAW_DELETED_PROCESSED_VERIFIED"}:
            reason.append(f"RAW:{row['raw_final_state']}")
        if row["metadata_match_status"] not in {"PASS", "NOT_AVAILABLE"}:
            reason.append(f"METADATA:{row['metadata_match_status']}")
        if not normalize_bool(row.get("prepare_complete", False)):
            reason.append("PREPARE:INCOMPLETE")
        reasons.append("|".join(reason))
    master["audit_reason"] = reasons
    summary = pd.DataFrame(
        [
            ("expected_sources", len(master)),
            ("duplicate_units", int((master["duplicate_unit_key_count"].astype(int) > 1).sum())),
            ("raw_present", int(master["raw_final_state"].eq("RAW_PRESENT").sum())),
            (
                "raw_deleted_processed_verified",
                int(master["raw_final_state"].eq("RAW_DELETED_PROCESSED_VERIFIED").sum()),
            ),
            ("raw_missing_unrecovered", int((~raw_ok).sum())),
            ("metadata_review", int((~metadata_ok).sum())),
            ("prepare_complete", int(prepare_ok.sum())),
            ("prepare_incomplete", int((~prepare_ok).sum())),
            ("final_ready", int(master["final_status"].eq("READY").sum())),
            ("final_review", int(master["final_status"].ne("READY").sum())),
        ],
        columns=["metric", "value"],
    )
    return master, summary


def create_recovery_inputs(master: pd.DataFrame, output_dir: Path, stamp: str) -> Path:
    recovery_dir = output_dir / f"generated_recovery_inputs_{stamp}"
    recovery_dir.mkdir(parents=True, exist_ok=True)
    actions = [
        ("reprepare_or_review_log", "reprepare_sources"),
        ("redownload_or_check_raw", "redownload_sources"),
        ("review_site_metadata", "metadata_review_sources"),
    ]
    for action, name in actions:
        atomic_write_tsv(
            master[master["recommended_action"].eq(action)],
            recovery_dir / f"{name}.tsv",
        )
    reprepare = master[master["recommended_action"].eq("reprepare_or_review_log")]
    for batch_id, group in reprepare.groupby("batch_id"):
        values = [
            normalize_source_file(value)
            for value in group["source_file"]
            if normalize_source_file(value)
        ]
        (recovery_dir / f"{batch_id}_recovery_source_files.txt").write_text(
            "\n".join(values) + ("\n" if values else ""),
            encoding="utf-8",
        )
    return recovery_dir


def write_excel_report(output_dir: Path, stamp: str, sheets: dict[str, pd.DataFrame]) -> Path:
    path = output_dir / f"audit_report_{stamp}.xlsx"
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for sheet_name, frame in sheets.items():
            content = frame if not frame.empty else pd.DataFrame({"message": ["No rows"]})
            content.head(1_048_000).to_excel(writer, sheet_name=sheet_name[:31], index=False)
    return path


def write_run_manifest(
    config: dict[str, Any],
    output_dir: Path,
    stamp: str,
    status: str,
    summary: pd.DataFrame,
) -> Path:
    path = output_dir / f"audit_run_{stamp}.json"
    payload = {
        "run_id": stamp,
        "status": status,
        "started_or_completed_at_utc": utc_now(),
        "config_path": str(config["_config_path"]),
        "config_sha256": sha256_file(config["_config_path"]),
        "code_path": str(Path(__file__).resolve()),
        "code_sha256": sha256_file(Path(__file__).resolve()),
        "summary": dict(zip(summary["metric"], summary["value"].astype(int))),
    }
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    temporary.replace(path)
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=None)
    parser.add_argument("--create-recovery-inputs", action="store_true")
    parser.add_argument("--fail-on-review", action="store_true")
    args = parser.parse_args()

    global pd
    pd = ensure_pandas()

    config = load_config(args.config)
    stamp = stamp_now()
    output_dir = root_path(config, config["audit_output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    expected = build_expected_manifest(config)
    raw = audit_raw_downloads(config, expected)
    metadata = audit_site_metadata(config, expected)
    prepare, batches = audit_prepare_outputs(config, expected)
    master, summary = build_master_manifest(expected, raw, metadata, prepare)

    write_current_and_snapshot(expected, output_dir, "expected_manifest", stamp)
    write_current_and_snapshot(raw, output_dir, "raw_download_audit", stamp)
    write_current_and_snapshot(metadata, output_dir, "site_metadata_audit", stamp)
    write_current_and_snapshot(prepare, output_dir, "prepare_source_audit", stamp)
    write_current_and_snapshot(batches, output_dir, "batch_audit", stamp)
    write_current_and_snapshot(master, output_dir, "pipeline_master_manifest", stamp)
    write_current_and_snapshot(summary, output_dir, "audit_summary", stamp)
    atomic_write_tsv(raw[~raw["raw_audit_status"].eq("RAW_PRESENT")], output_dir / "download_review.current.tsv")
    atomic_write_tsv(metadata[~metadata["metadata_match_status"].isin(["PASS", "NOT_AVAILABLE"])], output_dir / "metadata_review.current.tsv")
    atomic_write_tsv(prepare[~prepare["prepare_qc_status"].eq("SUCCESS_NONEMPTY")], output_dir / "prepare_review.current.tsv")
    atomic_write_tsv(master[master["final_status"].ne("READY")], output_dir / "final_review.current.tsv")
    atomic_write_tsv(master[master["final_status"].eq("READY")], output_dir / "ready_for_next_step.current.tsv")

    excel = write_excel_report(
        output_dir,
        stamp,
        {
            "summary": summary,
            "master": master,
            "batch_audit": batches,
            "raw_audit": raw,
            "metadata_audit": metadata,
            "prepare_audit": prepare,
        },
    )
    recovery_dir = None
    if args.create_recovery_inputs:
        recovery_dir = create_recovery_inputs(master, output_dir, stamp)
    final_review = int(summary.loc[summary["metric"].eq("final_review"), "value"].iloc[0])
    status = "READY_FOR_NEXT_STEP" if final_review == 0 else "REVIEW_REQUIRED"
    run_manifest = write_run_manifest(config, output_dir, stamp, status, summary)

    print(summary.to_string(index=False))
    print("STATUS:", status)
    print("[OK] Excel report:", excel)
    print("[OK] Run manifest:", run_manifest)
    if recovery_dir:
        print("[OK] Recovery inputs:", recovery_dir)
    if args.fail_on_review and final_review:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
