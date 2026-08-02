#!/usr/bin/env python3
"""Build a GRCh38 gene-coordinate table from a UKB-PPP download manifest.

This is a user-run setup helper. It queries the Ensembl REST symbol lookup API;
the batch runner itself remains offline and never invents gene coordinates.
"""
from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download-manifest", required=True, help="UKB-PPP TSV containing gene_symbol.")
    parser.add_argument("--output", required=True, help="Output GRCh38 coordinate TSV.")
    parser.add_argument("--missing-output", help="Optional newline-delimited unresolved symbol report.")
    parser.add_argument("--request-delay", type=float, default=0.1, help="Seconds between Ensembl requests.")
    parser.add_argument("--retries", type=int, default=3)
    return parser.parse_args()


def read_genes(path: Path) -> list[str]:
    if not path.is_file():
        raise SystemExit(f"[ERROR] Download manifest not found: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "gene_symbol" not in (reader.fieldnames or []):
            raise SystemExit("[ERROR] Download manifest is missing gene_symbol")
        return sorted({row["gene_symbol"].strip().upper() for row in reader if row["gene_symbol"].strip()})


def lookup_gene(symbol: str, retries: int) -> dict[str, object] | None:
    encoded = urllib.parse.quote(symbol, safe="")
    url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{encoded}?content-type=application/json"
    for attempt in range(retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=30) as response:
                return json.load(response)
        except urllib.error.HTTPError as error:
            if error.code == 400:
                return None
            if error.code != 429 and error.code < 500:
                raise
        except (urllib.error.URLError, TimeoutError):
            pass
        if attempt < retries:
            time.sleep(2**attempt)
    return None


def main() -> None:
    args = parse_args()
    if args.request_delay < 0 or args.retries < 0:
        raise SystemExit("[ERROR] --request-delay and --retries must be non-negative")
    genes = read_genes(Path(args.download_manifest))
    if not genes:
        raise SystemExit("[ERROR] Download manifest contains no gene symbols")

    rows: list[dict[str, object]] = []
    missing: list[str] = []
    for number, gene in enumerate(genes, start=1):
        record = lookup_gene(gene, args.retries)
        if not record or not all(record.get(key) is not None for key in ("seq_region_name", "start", "end")):
            missing.append(gene)
        else:
            rows.append({
                "gene_symbol": gene,
                "chr": record["seq_region_name"],
                "start": record["start"],
                "end": record["end"],
                "genome_build": "GRCh38",
            })
        if number % 100 == 0 or number == len(genes):
            print(f"[INFO] Ensembl lookup: {number}/{len(genes)}; resolved={len(rows)}; missing={len(missing)}", flush=True)
        if number < len(genes):
            time.sleep(args.request_delay)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gene_symbol", "chr", "start", "end", "genome_build"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(output)

    missing_output = Path(args.missing_output) if args.missing_output else output.with_suffix(".missing.txt")
    missing_output.write_text("".join(f"{gene}\n" for gene in missing), encoding="utf-8")
    print(f"[OK] Wrote {len(rows)} coordinates: {output}")
    if missing:
        raise SystemExit(f"[ERROR] {len(missing)} gene symbol(s) unresolved; see {missing_output}")


if __name__ == "__main__":
    main()
