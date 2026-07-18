#!/usr/bin/env python3
"""Download selected GIGASTROKE GWAS Catalog summary statistics for Colab setup.

This script is intended for data setup in Google Colab/Drive, not for PR review.
It discovers the summary-statistics gzip file in each GCST directory and downloads
it with wget -c into the project raw-data directory.
"""
from __future__ import annotations

import argparse
import os
import re
import subprocess
from pathlib import Path
from urllib.parse import urljoin
from urllib.request import urlopen

TARGETS = {
    "EUR": {
        "GCST90104540": "ischemic_stroke_EUR",
        "GCST90104542": "large_artery_stroke_EUR",
        "GCST90104541": "cardioembolic_stroke_EUR",
        "GCST90104543": "small_vessel_stroke_EUR",
    },
    "EAS": {
        "GCST90104545": "ischemic_stroke_EAS",
        "GCST90104547": "large_artery_stroke_EAS",
        "GCST90104546": "cardioembolic_stroke_EAS",
        "GCST90104548": "small_vessel_stroke_EAS",
    },
}
BASE_ROOT = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90104001-GCST90105000/"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=os.environ.get("PROJECT_ROOT", "/content/drive/MyDrive/IS_Analysis_V2"))
    parser.add_argument("--outdir", default="data/rawdata/outcome/gigastroke_gwas_catalog")
    parser.add_argument("--ancestry", choices=["EUR", "EAS", "ALL"], default="ALL")
    parser.add_argument("--dry-run", action="store_true", help="Print selected URLs and outputs without downloading.")
    return parser.parse_args()


def summary_urls(gcst: str) -> list[str]:
    url = f"{BASE_ROOT}{gcst}/"
    html = urlopen(url).read().decode("utf-8", errors="ignore")
    hrefs = re.findall(r'href="([^"]+)"', html)
    files = []
    for href in hrefs:
        if href.startswith(("?", "/")) or href == "../":
            continue
        if href.endswith((".tsv.gz", ".h.tsv.gz", ".txt.gz", ".gz")):
            files.append(urljoin(url, href))
    return sorted(set(files))


def output_extension(filename: str) -> str:
    if filename.endswith(".tsv.gz"):
        return ".tsv.gz"
    if filename.endswith(".txt.gz"):
        return ".txt.gz"
    return ".gz"


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root)
    outdir = project_root / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    ancestries = ["EUR", "EAS"] if args.ancestry == "ALL" else [args.ancestry]
    for ancestry in ancestries:
        for gcst, label in TARGETS[ancestry].items():
            url = f"{BASE_ROOT}{gcst}/"
            print(f"\n[INFO] {gcst} {label}")
            print(f"[INFO] URL: {url}")
            files = summary_urls(gcst)
            if not files:
                print("[WARN] gz summary file not found")
                continue
            file_url = files[0]
            original_name = file_url.rsplit("/", 1)[-1]
            outpath = outdir / f"{gcst}_{label}{output_extension(original_name)}"
            print(f"[DOWNLOAD] {original_name}")
            print(f"[SAVE AS] {outpath}")
            if not args.dry_run:
                subprocess.run(["wget", "-c", "-O", str(outpath), file_url], check=True)

    print("\n[INFO] Done")


if __name__ == "__main__":
    main()
