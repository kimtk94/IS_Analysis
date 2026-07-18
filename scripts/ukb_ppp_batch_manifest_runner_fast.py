#!/usr/bin/env python3
"""Run restartable UKB-PPP exposure batches from real EUR/EAS source archives.

Each batch contains paired EUR/EAS data for a fixed number of genes (10 by
default).  With ``--download-manifest``, the runner downloads only the source
archives needed for the current batch, validates them, runs the exposure filter
once per ancestry, and records all download and processing evidence.  The
canonical outputs are kept separate by ancestry; no EUR/EAS rows are mixed in a
batch result.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import os
import re
import shutil
import subprocess
import tarfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ANCESTRIES = ("EUR", "EAS")
REQUIRED_DOWNLOAD_COLUMNS = {"ancestry", "gene_symbol", "source_file", "url"}


def now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def safe_gene(gene: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", gene)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_atomic(df: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    df.to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, path)


def ensure_pandas() -> Any:
    if importlib.util.find_spec("pandas") is None:
        raise SystemExit("[ERROR] Missing Python dependency: pandas.")
    import pandas
    return pandas


def scan_valid(base: Path, pd: Any) -> Any:
    rows = []
    for ancestry in ANCESTRIES:
        folder = base / ancestry
        if not folder.exists():
            continue
        for path in sorted(folder.iterdir(), key=lambda item: item.name):
            if not path.is_file():
                continue
            rows.append({
                "ancestry": ancestry,
                "gene_symbol": path.name.split("_")[0].upper() if "_" in path.name else "",
                "source_file": path.name,
                "file_path": str(path),
                "size_bytes": path.stat().st_size,
                "valid_tar": path.suffix == ".tar" and path.stat().st_size > 0 and ".part" not in path.name,
            })
    return pd.DataFrame(rows)


def read_download_manifest(path: Path, pd: Any) -> Any:
    manifest = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    missing = REQUIRED_DOWNLOAD_COLUMNS - set(manifest.columns)
    if missing:
        raise SystemExit(f"[ERROR] Download manifest is missing columns: {', '.join(sorted(missing))}")
    manifest = manifest.copy()
    manifest["ancestry"] = manifest["ancestry"].str.upper().str.strip()
    manifest["gene_symbol"] = manifest["gene_symbol"].str.upper().str.strip()
    manifest["source_file"] = manifest["source_file"].map(lambda value: Path(value).name)
    manifest = manifest[manifest["ancestry"].isin(ANCESTRIES) & (manifest["gene_symbol"] != "")].copy()
    if manifest.empty:
        raise SystemExit("[ERROR] Download manifest has no EUR/EAS source rows")
    if manifest.duplicated(["ancestry", "gene_symbol", "source_file"]).any():
        raise SystemExit("[ERROR] Download manifest has duplicate ancestry/gene/source_file rows")
    return manifest


def paired_genes_from_manifest(manifest: Any) -> list[str]:
    groups = manifest.groupby("gene_symbol")["ancestry"].agg(lambda values: set(values))
    return sorted(gene for gene, ancestries in groups.items() if set(ANCESTRIES).issubset(ancestries))


def paired_genes_from_raw(audit: Any) -> list[str]:
    if audit.empty:
        return []
    valid = audit[audit["valid_tar"]]
    groups = valid.groupby("gene_symbol")["ancestry"].agg(lambda values: set(values))
    return sorted(gene for gene, ancestries in groups.items() if set(ANCESTRIES).issubset(ancestries))


def verify_archive(path: Path, expected_size: str, expected_sha256: str) -> tuple[bool, str, str]:
    if not path.exists() or path.stat().st_size == 0:
        return False, "missing_or_empty", ""
    if expected_size:
        try:
            if path.stat().st_size != int(expected_size):
                return False, "size_mismatch", ""
        except ValueError:
            return False, "invalid_expected_size", ""
    if not tarfile.is_tarfile(path):
        return False, "not_a_tar", ""
    try:
        with tarfile.open(path) as archive:
            archive.getmembers()
    except (tarfile.TarError, OSError) as error:
        return False, f"tar_invalid: {error}", ""
    observed_sha256 = sha256_file(path) if expected_sha256 else ""
    if expected_sha256 and observed_sha256.lower() != expected_sha256.lower():
        return False, "sha256_mismatch", observed_sha256
    return True, "verified", observed_sha256


def download_one(row: dict[str, str], destination: Path) -> tuple[bool, str]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    completed = destination.exists() and ".part" not in destination.name
    if completed:
        return True, "already_present"
    curl = shutil.which("curl")
    if not curl:
        return False, "curl_not_found"
    command = [curl, "--fail", "--location", "--retry", "3", "--continue-at", "-", "--output", str(destination), row["url"]]
    result = subprocess.run(command, check=False, text=True, capture_output=True)
    if result.returncode:
        detail = (result.stderr or result.stdout).strip().replace("\n", " ")[:500]
        return False, f"download_failed({result.returncode}): {detail}"
    return True, "downloaded"


def processing_status_allows_cleanup(status_path: Path, raw_files: list[Path]) -> tuple[bool, str]:
    """Require a successful terminal R status for every archive before deletion."""
    if not status_path.exists():
        return False, "missing_gene_status"
    with status_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    successful = {"completed", "completed_existing", "no_variants_after_filter"}
    status_by_tar = {
        str(Path(row["tar_file"]).resolve()): row.get("status", "")
        for row in rows
        if row.get("tar_file")
    }
    incomplete = [
        path.name for path in raw_files
        if status_by_tar.get(str(path.resolve())) not in successful
    ]
    if incomplete:
        return False, "incomplete_processing_status: " + ",".join(incomplete)
    return True, "all_sources_processed"


def remove_raw_files(raw_files: list[Path], status_path: Path) -> list[dict[str, str]]:
    """Delete only files backed by successful per-source processing status."""
    allowed, reason = processing_status_allows_cleanup(status_path, raw_files)
    if not allowed:
        return [{"raw_file": str(path), "action": "retained", "reason": reason} for path in raw_files]
    rows = []
    for path in raw_files:
        try:
            path.unlink()
            rows.append({"raw_file": str(path), "action": "deleted", "reason": "all_sources_processed"})
        except OSError as error:
            rows.append({"raw_file": str(path), "action": "delete_failed", "reason": str(error)})
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", default="data/rawdata/pqtl/selected_targets", help="Raw archive root containing EUR/ and EAS/.")
    parser.add_argument("--download-manifest", help="TSV: ancestry, gene_symbol, source_file, url; optional expected_size_bytes, sha256.")
    parser.add_argument("--outdir", default="results/exposure_batches")
    parser.add_argument("--qc-dir", default="results/qc/batch_pipeline")
    parser.add_argument("--batch-size", type=int, default=10)
    parser.add_argument("--only-batch", help="Comma-separated batch IDs, e.g. batch_001,batch_002.")
    parser.add_argument("--max-batches", type=int)
    parser.add_argument("--p-threshold", default="5e-8")
    parser.add_argument("--rscript", default="scripts/01_prepare_exposure_fast.R")
    parser.add_argument("--tmpdir", default="/content/ukbppp_tmp")
    parser.add_argument("--download-only", action="store_true")
    parser.add_argument(
        "--delete-raw-after-processing",
        action="store_true",
        help="Delete a batch's verified raw archives only after both ancestry outputs and every per-source status succeed.",
    )
    parser.add_argument("--run", action="store_true", help="Run downloads and exposure preparation. Without this, write the plan only.")
    parser.add_argument("--stop-on-error", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.batch_size < 1:
        raise SystemExit("[ERROR] --batch-size must be positive")
    pd = ensure_pandas()
    base, outdir, qc_dir = Path(args.base), Path(args.outdir), Path(args.qc_dir)
    manifest_path = qc_dir / "batch_manifest.tsv"
    download_manifest = read_download_manifest(Path(args.download_manifest), pd) if args.download_manifest else None
    raw_audit = scan_valid(base, pd)
    raw_audit.to_csv(qc_dir / "raw_file_audit.tsv", sep="\t", index=False)
    genes = paired_genes_from_manifest(download_manifest) if download_manifest is not None else paired_genes_from_raw(raw_audit)
    if not genes:
        raise SystemExit("[ERROR] No EUR/EAS paired genes found. Supply --download-manifest or populate --base.")

    rows: list[dict[str, object]] = []
    for offset in range(0, len(genes), args.batch_size):
        batch_genes = genes[offset : offset + args.batch_size]
        batch_id = f"batch_{offset // args.batch_size + 1:03d}"
        rows.append({"batch_id": batch_id, "n_genes": len(batch_genes), "genes": ",".join(batch_genes), "status": "pending", "raw_cleanup": "not_requested", **{f"{anc.lower()}_output": str(outdir / anc / f"exposure_{batch_id}.tsv") for anc in ANCESTRIES}})
    batch_df = pd.DataFrame(rows)
    selected_batches = batch_df
    if args.only_batch:
        wanted = {value.strip() for value in args.only_batch.split(",") if value.strip()}
        selected_batches = batch_df[batch_df["batch_id"].isin(wanted)].copy()
    if args.max_batches is not None:
        selected_batches = selected_batches.head(args.max_batches).copy()
    if selected_batches.empty:
        raise SystemExit("[ERROR] No batches selected")
    write_atomic(batch_df, manifest_path)

    for index, batch in selected_batches.iterrows():
        batch_id, batch_genes = str(batch.batch_id), str(batch.genes).split(",")
        print(f"\n===== {batch_id}: {len(batch_genes)} genes =====")
        source_rows: list[dict[str, str]] = []
        if download_manifest is not None:
            selected = download_manifest[download_manifest["gene_symbol"].isin(batch_genes)]
            source_rows = selected.to_dict("records")
        else:
            for gene in batch_genes:
                for ancestry in ANCESTRIES:
                    source_rows.append({"gene_symbol": gene, "ancestry": ancestry, "source_file": "", "url": ""})

        audit_rows = []
        raw_files_by_ancestry: dict[str, list[Path]] = {ancestry: [] for ancestry in ANCESTRIES}
        download_failed = False
        for source in source_rows:
            ancestry, gene = source["ancestry"], source["gene_symbol"]
            destination = base / ancestry / source["source_file"] if source["source_file"] else base / ancestry
            action, message = "not_requested", "existing_raw_mode"
            if args.run and download_manifest is not None:
                ok, message = download_one(source, destination)
                action = "download" if ok else "download_failed"
                if not ok:
                    download_failed = True
            if source["source_file"]:
                raw_files_by_ancestry[ancestry].append(destination)
                ok, verification, observed_sha = verify_archive(destination, source.get("expected_size_bytes", ""), source.get("sha256", ""))
            else:
                matches = list((base / ancestry).glob(f"{gene}_*.tar"))
                raw_files_by_ancestry[ancestry].extend(matches)
                ok, verification, observed_sha = bool(matches), "raw_archive_found" if matches else "missing_raw_archive", ""
            audit_rows.append({"timestamp": now(), "batch_id": batch_id, "gene_symbol": gene, "ancestry": ancestry, "source_file": source["source_file"], "destination": str(destination), "action": action, "download_message": message, "verified": ok, "verification": verification, "observed_sha256": observed_sha})
            if not ok:
                download_failed = True
        download_audit = qc_dir / "downloads" / f"{batch_id}.tsv"
        write_atomic(pd.DataFrame(audit_rows), download_audit)
        if download_failed:
            batch_df.loc[index, "status"] = "download_or_verification_failed"
            write_atomic(batch_df, manifest_path)
            if args.stop_on_error:
                raise SystemExit(f"[ERROR] {batch_id} download/verification failed; see {download_audit}")
            continue
        if not args.run or args.download_only:
            batch_df.loc[index, "status"] = "download_verified" if args.run else "planned"
            write_atomic(batch_df, manifest_path)
            continue

        gene_file = qc_dir / "gene_batches" / f"{batch_id}.txt"
        gene_file.parent.mkdir(parents=True, exist_ok=True)
        gene_file.write_text("\n".join(batch_genes) + "\n", encoding="utf-8")
        failed = False
        for ancestry in ANCESTRIES:
            output = outdir / ancestry / f"exposure_{batch_id}.tsv"
            log = qc_dir / "processing_logs" / f"{batch_id}_{ancestry}.log"
            log.parent.mkdir(parents=True, exist_ok=True)
            command = ["Rscript", args.rscript, "--gene-file", str(gene_file), "--batch-id", batch_id, "--outdir", str(outdir / ancestry), "--batch-output", str(output), "--rawdir", str(base), "--tmpdir", args.tmpdir, "--p-threshold", str(args.p_threshold), "--ancestries", ancestry]
            result = subprocess.run(command, check=False, text=True, capture_output=True)
            log.write_text(result.stdout + result.stderr, encoding="utf-8")
            if result.returncode or not output.exists():
                failed = True
                print(f"[ERROR] {batch_id} {ancestry}; see {log}")
        if failed:
            batch_df.loc[index, "status"] = "processing_failed"
            batch_df.loc[index, "raw_cleanup"] = "retained_processing_failed"
        elif args.delete_raw_after_processing:
            cleanup_rows = []
            for ancestry in ANCESTRIES:
                status_path = outdir / ancestry / "logs" / f"{batch_id}_gene_status.tsv"
                for row in remove_raw_files(raw_files_by_ancestry[ancestry], status_path):
                    cleanup_rows.append({"batch_id": batch_id, "ancestry": ancestry, **row})
            cleanup_path = qc_dir / "raw_cleanup" / f"{batch_id}.tsv"
            write_atomic(pd.DataFrame(cleanup_rows), cleanup_path)
            cleanup_failed = any(row["action"] != "deleted" for row in cleanup_rows)
            batch_df.loc[index, "status"] = "completed_raw_retained" if cleanup_failed else "completed_raw_deleted"
            batch_df.loc[index, "raw_cleanup"] = str(cleanup_path)
        else:
            batch_df.loc[index, "status"] = "completed"
            batch_df.loc[index, "raw_cleanup"] = "retained_flag_not_set"
        write_atomic(batch_df, manifest_path)
        if failed and args.stop_on_error:
            raise SystemExit(f"[ERROR] {batch_id} processing failed")

    print(f"[INFO] Batch manifest: {manifest_path}")


if __name__ == "__main__":
    main()
