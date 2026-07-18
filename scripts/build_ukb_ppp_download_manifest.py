#!/usr/bin/env python3
"""Build a UKB-PPP download manifest from Synapse file metadata.

The input target TSV identifies each required EUR/EAS archive.  Supply either a
Synapse metadata export (suitable for review and reproducible manifest builds)
or a Synapse parent folder, which is queried with ``synapseclient`` in a
user-run data-setup environment.  Archive bytes are never downloaded by this
command.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import os
from pathlib import Path
from typing import Iterable

REQUIRED_TARGET_COLUMNS = {"ancestry", "gene_symbol", "source_file"}
OUTPUT_COLUMNS = [
    "ancestry", "gene_symbol", "source_file", "url", "expected_size_bytes",
    "sha256", "md5", "synapse_id",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV has no header: {path}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def first_value(row: dict[str, str], names: Iterable[str]) -> str:
    lowered = {key.lower(): value for key, value in row.items()}
    for name in names:
        value = lowered.get(name.lower(), "")
        if value:
            return value
    return ""


def normalize_metadata(row: dict[str, str]) -> dict[str, str]:
    """Normalize common Synapse export/API field names to manifest metadata."""
    synapse_id = first_value(row, ("synapse_id", "id", "entityId", "entity_id"))
    filename = first_value(row, ("source_file", "name", "fileName", "file_name"))
    size = first_value(row, ("expected_size_bytes", "dataFileSizeBytes", "contentSize", "content_size"))
    md5 = first_value(row, ("md5", "contentMd5", "content_md5"))
    sha256 = first_value(row, ("sha256", "contentSha256", "content_sha256"))
    return {"synapse_id": synapse_id, "source_file": Path(filename).name, "expected_size_bytes": size, "md5": md5, "sha256": sha256}


def metadata_from_parent(parent_id: str, token: str) -> list[dict[str, str]]:
    """Fetch file-handle metadata only; do not materialize Synapse file bytes."""
    if importlib.util.find_spec("synapseclient") is None:
        raise RuntimeError("synapseclient is required for --synapse-parent-id; install it during the user-run setup phase.")
    import synapseclient

    syn = synapseclient.Synapse()
    syn.login(authToken=token or None, silent=not bool(token))
    rows = []
    for child in syn.getChildren(parent_id, includeTypes=["file"]):
        entity = syn.restGET(f"/repo/v1/entity/{child['id']}")
        handle_id = entity.get("dataFileHandleId")
        if not handle_id:
            continue
        handle = syn.restGET(f"/file/v1/fileHandle/{handle_id}")
        rows.append(normalize_metadata({
            "id": child["id"], "name": child["name"], "contentSize": str(handle.get("contentSize", "")),
            "contentMd5": handle.get("contentMd5", ""), "contentSha256": handle.get("contentSha256", ""),
        }))
    return rows


def build_manifest(targets: list[dict[str, str]], metadata: list[dict[str, str]]) -> list[dict[str, str]]:
    missing = REQUIRED_TARGET_COLUMNS - set(targets[0] if targets else {})
    if missing:
        raise ValueError("Target TSV is missing columns: " + ", ".join(sorted(missing)))
    by_id = {row["synapse_id"]: row for row in metadata if row["synapse_id"]}
    by_name: dict[str, list[dict[str, str]]] = {}
    for row in metadata:
        by_name.setdefault(row["source_file"], []).append(row)
    output = []
    for target in targets:
        ancestry = target["ancestry"].upper()
        gene = target["gene_symbol"].upper()
        name = Path(target["source_file"]).name
        candidates = [by_id[target["synapse_id"]]] if target.get("synapse_id", "") in by_id else by_name.get(name, [])
        if len(candidates) != 1:
            detail = "not found" if not candidates else "ambiguous; add synapse_id to the target TSV"
            raise ValueError(f"Synapse metadata for {ancestry}/{gene}/{name}: {detail}")
        item = candidates[0]
        if not item["expected_size_bytes"].isdigit() or int(item["expected_size_bytes"]) <= 0:
            raise ValueError(f"Synapse metadata has no positive content size for {item['synapse_id']} ({name})")
        output.append({
            "ancestry": ancestry, "gene_symbol": gene, "source_file": name,
            "url": f"https://www.synapse.org/Synapse:{item['synapse_id']}",
            "expected_size_bytes": item["expected_size_bytes"], "sha256": item["sha256"],
            "md5": item["md5"], "synapse_id": item["synapse_id"],
        })
    if len({(row["ancestry"], row["gene_symbol"], row["source_file"]) for row in output}) != len(output):
        raise ValueError("Target TSV contains duplicate ancestry/gene/source_file rows")
    return output


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", required=True, help="TSV with ancestry, gene_symbol, source_file; optional synapse_id.")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--synapse-metadata-file", help="TSV exported from Synapse containing id/name and file-handle metadata.")
    source.add_argument("--synapse-parent-id", help="Synapse folder ID to query for file metadata.")
    parser.add_argument("--output", required=True, help="Output download manifest TSV.")
    parser.add_argument("--synapse-token", default=os.environ.get("SYNAPSE_AUTH_TOKEN", ""), help="Synapse personal access token; defaults to SYNAPSE_AUTH_TOKEN.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    targets = read_tsv(Path(args.targets))
    raw_metadata = read_tsv(Path(args.synapse_metadata_file)) if args.synapse_metadata_file else metadata_from_parent(args.synapse_parent_id, args.synapse_token)
    rows = build_manifest(targets, [normalize_metadata(row) for row in raw_metadata])
    write_tsv(Path(args.output), rows)
    fingerprint = hashlib.sha256(Path(args.output).read_bytes()).hexdigest()
    print(f"[OK] Wrote {len(rows)} manifest rows: {args.output}")
    print(f"[OK] Manifest SHA-256: {fingerprint}")


if __name__ == "__main__":
    main()
