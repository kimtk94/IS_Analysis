#!/usr/bin/env python3
"""Audit UKB-PPP rawdata files saved in Google Drive.

This script checks the actual files under ``data/rawdata/pqtl/selected_targets``
(or another ancestry-aware ``selected_targets`` directory), not only the download
manifest. It produces TSV outputs that can be used for final gene prioritization
and download cleanup.

Recommended Colab usage:
    python scripts/Rawdatacheck.py \
      --base /content/drive/MyDrive/IS_Analysis/data/rawdata/pqtl/selected_targets \
      --manifest /content/drive/MyDrive/IS_Analysis/results/qc/ukb_ppp_selected_download_manifest.tsv \
      --outdir /content/drive/MyDrive/IS_Analysis/results/qc
"""
from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit UKB-PPP selected target rawdata files.")
    parser.add_argument("--base", required=True, help="selected_targets directory containing EUR/EAS/etc folders")
    parser.add_argument("--manifest", default=None, help="Optional download manifest TSV")
    parser.add_argument("--outdir", required=True, help="Output directory for audit TSV files")
    parser.add_argument("--delete-temp", action="store_true", help="Delete temporary/broken Synapse files after writing reports")
    return parser.parse_args()


def read_manifest(path: Path | None) -> tuple[dict[str, dict[str, str]], list[dict[str, str]]]:
    if not path or not path.exists():
        return {}, []
    rows = []
    by_name = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append(row)
            by_name[row.get("name", "")] = row
    return by_name, rows


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    args = parse_args()
    base = Path(args.base)
    outdir = Path(args.outdir)
    manifest_path = Path(args.manifest) if args.manifest else None

    if not base.exists():
        raise SystemExit(f"[ERROR] base directory does not exist: {base}")

    manifest_by_name, _manifest_rows = read_manifest(manifest_path)

    rows: list[dict[str, Any]] = []
    for ancestry_dir in sorted(path for path in base.iterdir() if path.is_dir()):
        ancestry = ancestry_dir.name
        for path in sorted(ancestry_dir.iterdir(), key=lambda item: item.name):
            if not path.is_file():
                continue
            name = path.name
            size = path.stat().st_size
            is_temp_synapse = ".synapse_download_" in name
            is_tar = name.endswith(".tar")
            gene_symbol = name.split("_")[0] if "_" in name else ""
            clean_tar_name = re.sub(r"\.synapse_download_.*$", "", name)
            manifest_row = manifest_by_name.get(clean_tar_name) or manifest_by_name.get(name) or {}
            expected_size = int(manifest_row.get("dataFileSizeBytes") or 0)
            size_matches_manifest = bool(expected_size and size == expected_size)
            is_valid_raw_tar = bool(is_tar and (not is_temp_synapse) and size > 0)
            rows.append({
                "ancestry": ancestry,
                "file_name": name,
                "clean_tar_name": clean_tar_name,
                "gene_symbol": gene_symbol,
                "size_bytes": size,
                "size_gb": round(size / 1024**3, 6),
                "is_tar": is_tar,
                "is_temp_synapse": is_temp_synapse,
                "is_zero_size": size == 0,
                "is_valid_raw_tar": is_valid_raw_tar,
                "expected_size_from_manifest": expected_size,
                "size_matches_manifest": size_matches_manifest,
                "manifest_status": manifest_row.get("status", ""),
                "manifest_message": manifest_row.get("message", ""),
            })

    valid = [row for row in rows if row["is_valid_raw_tar"]]
    invalid = [row for row in rows if not row["is_valid_raw_tar"]]

    all_fields = [
        "ancestry", "file_name", "clean_tar_name", "gene_symbol", "size_bytes", "size_gb",
        "is_tar", "is_temp_synapse", "is_zero_size", "is_valid_raw_tar",
        "expected_size_from_manifest", "size_matches_manifest", "manifest_status", "manifest_message",
    ]
    write_tsv(outdir / "ukb_ppp_rawdata_file_audit.tsv", rows, all_fields)
    write_tsv(outdir / "ukb_ppp_rawdata_invalid_or_temp_files.tsv", invalid, all_fields)

    gene_rows = sorted({(row["ancestry"], row["gene_symbol"]) for row in valid})
    write_tsv(
        outdir / "ukb_ppp_rawdata_valid_gene_list.tsv",
        [{"ancestry": ancestry, "gene_symbol": gene} for ancestry, gene in gene_rows],
        ["ancestry", "gene_symbol"],
    )

    group: dict[tuple[str, str], list[str]] = defaultdict(list)
    for row in valid:
        group[(row["ancestry"], row["gene_symbol"])].append(row["file_name"])
    dup_rows = []
    for (ancestry, gene), files in sorted(group.items()):
        if len(files) > 1:
            dup_rows.append({
                "ancestry": ancestry,
                "gene_symbol": gene,
                "n_files": len(files),
                "file_names": "; ".join(files),
            })
    write_tsv(
        outdir / "ukb_ppp_rawdata_duplicate_gene_files.tsv",
        dup_rows,
        ["ancestry", "gene_symbol", "n_files", "file_names"],
    )

    summary = []
    by_ancestry = sorted({row["ancestry"] for row in rows})
    for ancestry in by_ancestry:
        sub = [row for row in rows if row["ancestry"] == ancestry]
        valid_sub = [row for row in sub if row["is_valid_raw_tar"]]
        invalid_sub = [row for row in sub if not row["is_valid_raw_tar"]]
        summary.append({
            "ancestry": ancestry,
            "all_files": len(sub),
            "valid_tar_files": len(valid_sub),
            "valid_unique_genes": len({row["gene_symbol"] for row in valid_sub}),
            "invalid_or_temp_files": len(invalid_sub),
            "zero_size_files": sum(1 for row in sub if row["is_zero_size"]),
            "temp_synapse_files": sum(1 for row in sub if row["is_temp_synapse"]),
            "valid_total_gb": round(sum(row["size_bytes"] for row in valid_sub) / 1024**3, 3),
            "manifest_size_match_files": sum(1 for row in valid_sub if row["size_matches_manifest"]),
        })
    write_tsv(outdir / "ukb_ppp_rawdata_audit_summary.tsv", summary, list(summary[0].keys()) if summary else ["ancestry"])

    print("[INFO] Rawdata audit complete")
    for item in summary:
        print(item)
    print("[INFO] Outputs:")
    for filename in [
        "ukb_ppp_rawdata_audit_summary.tsv",
        "ukb_ppp_rawdata_file_audit.tsv",
        "ukb_ppp_rawdata_valid_gene_list.tsv",
        "ukb_ppp_rawdata_invalid_or_temp_files.tsv",
        "ukb_ppp_rawdata_duplicate_gene_files.tsv",
    ]:
        print(" -", outdir / filename)

    if args.delete_temp:
        targets = [base / row["ancestry"] / row["file_name"] for row in invalid if row["is_temp_synapse"] or row["is_zero_size"]]
        print(f"[WARN] deleting {len(targets)} temporary/zero-size files")
        for path in targets:
            if path.exists():
                path.unlink()
                print("deleted", path)


if __name__ == "__main__":
    main()
