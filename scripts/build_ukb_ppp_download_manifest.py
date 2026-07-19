#!/usr/bin/env python3
"""Build a UKB-PPP download manifest from Synapse file metadata.

The command enumerates archive metadata from explicitly supplied Synapse parent
folders and derives gene symbols from archive filenames. No gene target list or
archive bytes are required. Folder queries use ``synapseclient`` in a user-run
data-setup environment.
"""
from __future__ import annotations

import argparse
import asyncio
import csv
import hashlib
import importlib.util
import os
from pathlib import Path
from typing import Iterable

REQUIRED_TARGET_COLUMNS = {"ancestry", "gene_symbol", "source_file"}
OUTPUT_COLUMNS = [
    "ancestry", "gene_symbol", "source_file", "url", "expected_size_bytes",
    "sha256", "md5", "synapse_id", "synapse_parent_id",
]
ANCESTRIES = {
    "EUR": {"label": "European_discovery", "parent_id": "syn51365303", "aliases": {"EUR", "EUROPEAN", "EUROPEAN_DISCOVERY", "EUROPE"}},
    "AFR": {"label": "African", "parent_id": "syn51365304", "aliases": {"AFR", "AFRICAN", "AFRICA"}},
    "CSA": {"label": "Central_South_Asian", "parent_id": "syn51365305", "aliases": {"CSA", "CENTRAL_SOUTH_ASIAN", "CENTRAL-SOUTH-ASIAN", "SOUTH_ASIAN"}},
    "EAS": {"label": "East_Asian", "parent_id": "syn51365306", "aliases": {"EAS", "EAST_ASIAN", "EAST-ASIAN", "EASTASIAN"}},
    "MID": {"label": "Middle_East", "parent_id": "syn51365307", "aliases": {"MID", "MIDDLE_EAST", "MIDDLE-EAST", "MIDDLEEAST"}},
    "COMBINED": {"label": "Combined", "parent_id": "syn51365308", "aliases": {"COMBINED", "ALL", "META"}},
    "AMR": {"label": "American", "parent_id": "syn51500434", "aliases": {"AMR", "AMERICAN", "AMERICA"}},
}


def canonical_ancestry(value: str) -> str:
    """Return the configured ancestry code, accepting documented aliases."""
    normalized = value.upper().strip().replace(" ", "_")
    for code, details in ANCESTRIES.items():
        if normalized == code or normalized in details["aliases"]:
            return code
    raise ValueError(f"Unsupported ancestry: {value}. Supported codes: {', '.join(ANCESTRIES)}")


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
    ancestry = first_value(row, ("ancestry", "ancestry_code"))
    parent_id = first_value(row, ("synapse_parent_id", "parent_id", "parentId"))
    return {
        "synapse_id": synapse_id, "source_file": Path(filename).name,
        "expected_size_bytes": size, "md5": md5, "sha256": sha256,
        "ancestry": canonical_ancestry(ancestry) if ancestry else "", "synapse_parent_id": parent_id,
    }


def metadata_from_parent(parent_id: str, token: str, ancestry: str = "") -> list[dict[str, str]]:
    """Fetch file-handle metadata only; do not materialize Synapse file bytes."""
    if importlib.util.find_spec("synapseclient") is None:
        raise RuntimeError("synapseclient is required for --synapse-parent-id; install it during the user-run setup phase.")
    import synapseclient
    from synapseclient.api import get_children, get_entity
    from synapseclient.api.file_services import get_file_handle_for_download_async

    syn = synapseclient.Synapse()
    syn.login(authToken=token or None, silent=not bool(token))

    async def collect_metadata() -> list[dict[str, str]]:
        rows = []
        async for child in get_children(parent=parent_id, include_types=["file"], synapse_client=syn):
            entity = await get_entity(child["id"], synapse_client=syn)
            handle_id = entity.get("dataFileHandleId")
            if not handle_id:
                continue
            download_info = await get_file_handle_for_download_async(
                str(handle_id), child["id"], synapse_client=syn,
            )
            handle = download_info.get("fileHandle", {})
            rows.append(normalize_metadata({
                "id": child["id"], "name": child["name"], "contentSize": str(handle.get("contentSize", "")),
                "contentMd5": handle.get("contentMd5", ""), "contentSha256": handle.get("contentSha256", ""),
                "ancestry": ancestry, "synapse_parent_id": parent_id,
            }))
        return rows

    return asyncio.run(collect_metadata())


def parse_parent_spec(value: str) -> tuple[str, str]:
    """Parse the explicit ANCESTRY:SYNAPSE_PARENT_ID CLI value."""
    try:
        ancestry, parent_id = value.split(":", 1)
    except ValueError as error:
        raise ValueError(f"Invalid --synapse-parent {value!r}; use ANCESTRY:syn123") from error
    ancestry, parent_id = canonical_ancestry(ancestry), parent_id.strip()
    if not parent_id.startswith("syn") or not parent_id[3:].isdigit():
        raise ValueError(f"Invalid Synapse parent ID in --synapse-parent: {parent_id}")
    return ancestry, parent_id


def metadata_from_explicit_parents(parent_specs: list[str], token: str) -> list[dict[str, str]]:
    """Fetch metadata from user-specified parent IDs, retaining each ID in output."""
    parsed = [parse_parent_spec(value) for value in parent_specs]
    if len({ancestry for ancestry, _parent_id in parsed}) != len(parsed):
        raise ValueError("Specify exactly one --synapse-parent value per ancestry")
    rows = []
    for ancestry, parent_id in parsed:
        rows.extend(metadata_from_parent(parent_id, token, ancestry))
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
        ancestry = canonical_ancestry(target["ancestry"])
        gene = target["gene_symbol"].upper()
        name = Path(target["source_file"]).name
        candidates = [by_id[target["synapse_id"]]] if target.get("synapse_id", "") in by_id else by_name.get(name, [])
        candidates = [item for item in candidates if not item["ancestry"] or item["ancestry"] == ancestry]
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
            "md5": item["md5"], "synapse_id": item["synapse_id"], "synapse_parent_id": item["synapse_parent_id"],
        })
    if len({(row["ancestry"], row["gene_symbol"], row["source_file"]) for row in output}) != len(output):
        raise ValueError("Target TSV contains duplicate ancestry/gene/source_file rows")
    return output


def build_manifest_from_metadata(metadata: list[dict[str, str]]) -> list[dict[str, str]]:
    """Build all EUR/EAS manifest rows from explicitly scoped folder metadata."""
    targets = []
    for item in metadata:
        if item["ancestry"] not in {"EUR", "EAS"} or not item["source_file"].endswith(".tar"):
            continue
        gene = item["source_file"].split("_", 1)[0].upper()
        if not gene:
            raise ValueError(f"Cannot infer gene symbol from {item['source_file']}")
        targets.append({"ancestry": item["ancestry"], "gene_symbol": gene, "source_file": item["source_file"], "synapse_id": item["synapse_id"]})
    if not targets:
        raise ValueError("No EUR/EAS .tar archives found in the supplied Synapse metadata")
    return build_manifest(targets, metadata)


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
    parser.add_argument("--targets", help="Legacy TSV with ancestry, gene_symbol, source_file; optional synapse_id.")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--synapse-metadata-file", help="TSV exported from Synapse containing id/name and file-handle metadata.")
    source.add_argument("--synapse-parent-id", help="Synapse folder ID to query for file metadata.")
    source.add_argument("--synapse-parent", action="append", help="Explicit parent folder: ANCESTRY:syn123. Repeat for EUR and EAS.")
    parser.add_argument("--output", required=True, help="Output download manifest TSV.")
    parser.add_argument("--synapse-token", default=os.environ.get("SYNAPSE_AUTH_TOKEN", ""), help="Synapse personal access token; defaults to SYNAPSE_AUTH_TOKEN.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    targets = read_tsv(Path(args.targets)) if args.targets else []
    if args.synapse_parent and args.targets:
        raise ValueError("--synapse-parent builds the full manifest; do not supply --targets")
    if args.synapse_parent_id and not args.targets:
        raise ValueError("--synapse-parent-id requires --targets; use explicit --synapse-parent ANCESTRY:syn123 for a full manifest")
    if args.synapse_metadata_file:
        raw_metadata = read_tsv(Path(args.synapse_metadata_file))
    elif args.synapse_parent:
        raw_metadata = metadata_from_explicit_parents(args.synapse_parent, args.synapse_token)
    else:
        raw_metadata = metadata_from_parent(args.synapse_parent_id, args.synapse_token)
    metadata = [normalize_metadata(row) for row in raw_metadata]
    rows = build_manifest(targets, metadata) if targets else build_manifest_from_metadata(metadata)
    write_tsv(Path(args.output), rows)
    fingerprint = hashlib.sha256(Path(args.output).read_bytes()).hexdigest()
    print(f"[OK] Wrote {len(rows)} manifest rows: {args.output}")
    print(f"[OK] Manifest SHA-256: {fingerprint}")


if __name__ == "__main__":
    main()
