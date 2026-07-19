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
import time
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
COMPLETED_BATCH_STATUSES = {"completed", "completed_raw_deleted", "completed_raw_retained"}
RESTART_CLEANUP_STATUSES = {"running", "metadata_fetch_failed", "download_or_verification_failed", "processing_failed"}


def now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def elapsed(started: float) -> str:
    seconds = int(time.monotonic() - started)
    return f"{seconds // 60:02d}:{seconds % 60:02d}"


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


def hydrate_selected_synapse_metadata(
    sources: list[dict[str, str]], manifest: Any, manifest_path: Path, base: Path,
) -> list[dict[str, str]]:
    """Look up metadata only for selected files not already valid in ``base``."""
    reusable = 0

    def has_reusable_raw(row: dict[str, str]) -> bool:
        source_file = row.get("source_file", "")
        path = base / row.get("ancestry", "") / source_file
        return bool(source_file and path.exists() and path.stat().st_size > 0 and tarfile.is_tarfile(path))

    pending = [
        row["synapse_id"] for row in sources
        if row.get("synapse_id", "") and (not row.get("expected_size_bytes", "") or not row.get("md5", ""))
        and not has_reusable_raw(row)
    ]
    reusable = sum(1 for row in sources if has_reusable_raw(row))
    if reusable:
        print(f"[INFO] Reusing {reusable} valid existing raw archive(s); skipping Synapse metadata lookup", flush=True)
    if not pending:
        return sources
    from synapse_metadata import fetch_file_metadata_many, synapse_token

    def progress(completed: int, total: int) -> None:
        print(f"[INFO] Synapse metadata: {completed}/{total} files", flush=True)

    metadata = fetch_file_metadata_many(pending, token=synapse_token(), max_concurrency=8, progress=progress)
    for source in sources:
        details = metadata.get(source.get("synapse_id", ""))
        if not details:
            continue
        source.update(details)
        mask = manifest["synapse_id"].eq(source["synapse_id"])
        for column, value in details.items():
            if column not in manifest.columns:
                manifest[column] = ""
            manifest.loc[mask, column] = value
    write_atomic(manifest, manifest_path)
    return sources


def stage_existing_raw_archives(sources: list[dict[str, str]], base: Path, existing_base: Path) -> dict[str, Path]:
    """Link valid archives from a separate existing-data root into this run's base.

    Prefer a symlink so cleanup does not touch the original archive. Google
    Drive FUSE does not support symlinks, so fall back to an atomic local copy.
    """
    staged: dict[str, Path] = {}
    for source in sources:
        source_file = source.get("source_file", "")
        ancestry = source.get("ancestry", "")
        if not source_file or not ancestry:
            continue
        destination = base / ancestry / source_file
        candidate = existing_base / ancestry / source_file
        if destination.is_symlink() and not destination.exists():
            destination.unlink()
        if destination.exists() or not candidate.exists() or candidate.stat().st_size == 0:
            continue
        if not tarfile.is_tarfile(candidate):
            print(f"[WARN] Existing raw archive is not a valid tar and will not be staged: {candidate}", flush=True)
            continue
        destination.parent.mkdir(parents=True, exist_ok=True)
        try:
            destination.symlink_to(candidate.resolve())
        except OSError as error:
            temporary = destination.with_suffix(destination.suffix + ".part")
            try:
                shutil.copy2(candidate, temporary)
                os.replace(temporary, destination)
            except OSError:
                temporary.unlink(missing_ok=True)
                raise error
        staged[str(destination)] = candidate
    return staged


def prioritize_existing_raw_batches(batches: Any, manifest: Any, existing_base: Path, pd: Any) -> Any:
    """Run batches backed by a separate raw-data root before download-only batches."""
    coverage = []
    for number, (index, batch) in enumerate(batches.iterrows(), start=1):
        genes = str(batch["genes"]).split(",")
        sources = manifest[manifest["gene_symbol"].isin(genes)]
        available = 0
        for row in sources.to_dict("records"):
            path = existing_base / row["ancestry"] / row["source_file"]
            if path.exists() and path.stat().st_size > 0 and tarfile.is_tarfile(path):
                available += 1
        coverage.append({"index": index, "existing_raw_sources": available, "source_count": len(sources)})
        if number % 100 == 0 or number == len(batches):
            print(f"[INFO] Existing raw batch scan: {number}/{len(batches)}", flush=True)
    coverage_df = pd.DataFrame(coverage).set_index("index")
    ranked = batches.join(coverage_df)
    ranked = ranked.sort_values(
        ["existing_raw_sources", "source_count", "batch_id"], ascending=[False, False, True], kind="stable",
    )
    reusable = int((ranked["existing_raw_sources"] > 0).sum())
    complete = int((ranked["existing_raw_sources"] == ranked["source_count"]).sum())
    print(f"[INFO] Existing raw priority: {complete} fully available, {reusable} partially/fully available batches first", flush=True)
    return ranked


def prioritize_genes_from_existing_raw(genes: list[str], manifest: Any, existing_base: Path) -> list[str]:
    """Create the initial stable plan with fully available raw gene pairs first."""
    available, remaining = [], []
    for number, gene in enumerate(genes, start=1):
        sources = manifest[manifest["gene_symbol"].eq(gene)]
        present = {
            row["ancestry"] for row in sources.to_dict("records")
            if (existing_base / row["ancestry"] / row["source_file"]).exists()
            and (existing_base / row["ancestry"] / row["source_file"]).stat().st_size > 0
            and tarfile.is_tarfile(existing_base / row["ancestry"] / row["source_file"])
        }
        (available if set(ANCESTRIES).issubset(present) else remaining).append(gene)
        if number % 100 == 0 or number == len(genes):
            print(f"[INFO] Existing raw gene scan: {number}/{len(genes)}", flush=True)
    print(f"[INFO] Initial batch plan: {len(available)} fully available paired genes first", flush=True)
    return available + remaining


def restore_batch_state(batch_df: Any, manifest_path: Path, pd: Any) -> Any:
    """Carry terminal state forward so a restarted run does not repeat batches."""
    if not manifest_path.exists():
        return batch_df
    previous = pd.read_csv(manifest_path, sep="\t", dtype=str).fillna("")
    if not {"batch_id", "genes", "status"}.issubset(previous.columns):
        print(f"[WARN] Ignoring incompatible existing batch manifest: {manifest_path}")
        return batch_df
    previous = previous.drop_duplicates("batch_id", keep="last").set_index("batch_id")
    for index, row in batch_df.iterrows():
        batch_id = row["batch_id"]
        if batch_id not in previous.index or previous.at[batch_id, "genes"] != row["genes"]:
            continue
        for column in ("status", "raw_cleanup"):
            if column in previous.columns:
                batch_df.loc[index, column] = previous.at[batch_id, column]
    return batch_df


def print_batch_state(batch_df: Any, label: str) -> None:
    """Print a compact checkpoint summary suitable for Colab stdout."""
    counts = batch_df["status"].value_counts().sort_index()
    summary = ", ".join(f"{status}={count}" for status, count in counts.items())
    print(f"[INFO] {label}: {len(batch_df)} batches ({summary})", flush=True)


def clean_incomplete_batch_outputs(batch_id: str, outdir: Path, qc_dir: Path, pd: Any) -> None:
    """Remove stale derived outputs before retrying an interrupted batch.

    Verified raw archives are retained as input evidence and can be reused.
    """
    candidates = []
    for ancestry in ANCESTRIES:
        candidates.extend([
            outdir / ancestry / f"exposure_{batch_id}.tsv",
            outdir / ancestry / "logs" / f"{batch_id}_gene_status.tsv",
        ])
    rows = []
    for path in candidates:
        if path.exists():
            path.unlink()
            rows.append({"path": str(path), "action": "deleted_stale_derived_output"})
    cleanup_path = qc_dir / "partial_cleanup" / f"{batch_id}.tsv"
    write_atomic(pd.DataFrame(rows, columns=["path", "action"]), cleanup_path)
    print(f"[INFO] Reset incomplete derived outputs for {batch_id}; retained verified raw archives", flush=True)


def record_raw_cleanup_in_manifest(manifest: Any, cleanup_rows: list[dict[str, str]], batch_id: str) -> None:
    """Persist the raw lifecycle state in the source manifest before continuing."""
    for column in ("pipeline_batch_id", "raw_lifecycle", "raw_cleanup_at", "raw_cleanup_reason"):
        if column not in manifest.columns:
            manifest[column] = ""
    for row in cleanup_rows:
        source_file = Path(row["raw_file"]).name
        mask = (manifest["ancestry"] == row["ancestry"]) & (manifest["source_file"] == source_file)
        lifecycle = "deleted_after_processed" if row["action"] in {"deleted", "deleted_existing_source"} else "retained_after_processing"
        manifest.loc[mask, "pipeline_batch_id"] = batch_id
        manifest.loc[mask, "raw_lifecycle"] = lifecycle
        manifest.loc[mask, "raw_cleanup_at"] = now()
        manifest.loc[mask, "raw_cleanup_reason"] = row["reason"]


def paired_genes_from_manifest(manifest: Any) -> list[str]:
    groups = manifest.groupby("gene_symbol")["ancestry"].agg(lambda values: set(values))
    return sorted(gene for gene, ancestries in groups.items() if set(ANCESTRIES).issubset(ancestries))


def paired_genes_from_raw(audit: Any) -> list[str]:
    if audit.empty:
        return []
    valid = audit[audit["valid_tar"]]
    groups = valid.groupby("gene_symbol")["ancestry"].agg(lambda values: set(values))
    return sorted(gene for gene, ancestries in groups.items() if set(ANCESTRIES).issubset(ancestries))


def verify_archive(path: Path, expected_size: str, expected_sha256: str, expected_md5: str = "") -> tuple[bool, str, str]:
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
    if expected_md5:
        digest = hashlib.md5(usedforsecurity=False)
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                digest.update(chunk)
        if digest.hexdigest().lower() != expected_md5.lower():
            return False, "md5_mismatch", observed_sha256
    return True, "verified", observed_sha256


def download_one(row: dict[str, str], destination: Path) -> tuple[bool, str]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    completed = destination.exists() and ".part" not in destination.name
    if completed:
        return True, "already_present"
    synapse_id = row.get("synapse_id", "").strip()
    if synapse_id:
        if importlib.util.find_spec("synapseclient") is None:
            return False, "synapseclient_not_installed"
        import synapseclient
        from synapseclient.operations import FileOptions, get
        try:
            syn = synapseclient.Synapse()
            syn.login(authToken=os.environ.get("SYNAPSE_AUTH_TOKEN") or None, silent=not bool(os.environ.get("SYNAPSE_AUTH_TOKEN")))
            downloaded = Path(get(
                synapse_id=synapse_id,
                file_options=FileOptions(download_location=str(destination.parent), if_collision="overwrite.local"),
                synapse_client=syn,
            ).path)
            if downloaded.resolve() != destination.resolve():
                os.replace(downloaded, destination)
            return True, "downloaded_from_synapse"
        except Exception as error:  # Synapse client errors vary by version and authentication state.
            return False, f"synapse_download_failed: {str(error).replace(chr(10), ' ')[:500]}"
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
    parser.add_argument("--existing-raw-base", help="Optional separate raw archive root to validate and symlink into --base per batch.")
    parser.add_argument("--delete-existing-raw-after-processing", action="store_true", help="Delete staged originals in --existing-raw-base only after successful processing and staging cleanup.")
    parser.add_argument("--download-manifest", help="TSV: ancestry, gene_symbol, source_file, url; optional expected_size_bytes, sha256, md5, synapse_id.")
    parser.add_argument("--outdir", default="results/exposure_batches")
    parser.add_argument("--qc-dir", default="results/qc/batch_pipeline")
    parser.add_argument("--batch-size", type=int, default=10)
    parser.add_argument("--only-batch", help="Comma-separated batch IDs, e.g. batch_001,batch_002.")
    parser.add_argument("--max-batches", type=int)
    parser.add_argument("--rerun-completed", action="store_true", help="Run batches marked completed in an existing batch manifest.")
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
    if args.delete_raw_after_processing and not args.download_manifest:
        raise SystemExit("[ERROR] --delete-raw-after-processing requires --download-manifest for raw lifecycle tracking")
    if args.delete_existing_raw_after_processing and (not args.existing_raw_base or not args.delete_raw_after_processing):
        raise SystemExit("[ERROR] --delete-existing-raw-after-processing requires --existing-raw-base and --delete-raw-after-processing")
    pd = ensure_pandas()
    base, outdir, qc_dir = Path(args.base), Path(args.outdir), Path(args.qc_dir)
    existing_raw_base = Path(args.existing_raw_base) if args.existing_raw_base else None
    if existing_raw_base is None:
        print("[INFO] Existing raw mode: OFF (all missing sources use normal download handling)", flush=True)
    else:
        print(f"[INFO] Existing raw mode: ON (source={existing_raw_base}; staging={base})", flush=True)
    qc_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = qc_dir / "batch_manifest.tsv"
    download_manifest = read_download_manifest(Path(args.download_manifest), pd) if args.download_manifest else None
    raw_audit = scan_valid(base, pd)
    raw_audit.to_csv(qc_dir / "raw_file_audit.tsv", sep="\t", index=False)
    genes = paired_genes_from_manifest(download_manifest) if download_manifest is not None else paired_genes_from_raw(raw_audit)
    if not genes:
        raise SystemExit("[ERROR] No EUR/EAS paired genes found. Supply --download-manifest or populate --base.")
    if existing_raw_base is not None and download_manifest is not None and not manifest_path.exists():
        print("[INFO] Existing raw mode: building initial 10-gene plan from archive inventory", flush=True)
        genes = prioritize_genes_from_existing_raw(genes, download_manifest, existing_raw_base)

    rows: list[dict[str, object]] = []
    for offset in range(0, len(genes), args.batch_size):
        batch_genes = genes[offset : offset + args.batch_size]
        batch_id = f"batch_{offset // args.batch_size + 1:03d}"
        rows.append({"batch_id": batch_id, "n_genes": len(batch_genes), "genes": ",".join(batch_genes), "status": "pending", "raw_cleanup": "not_requested", **{f"{anc.lower()}_output": str(outdir / anc / f"exposure_{batch_id}.tsv") for anc in ANCESTRIES}})
    batch_df = pd.DataFrame(rows)
    batch_df = restore_batch_state(batch_df, manifest_path, pd)
    print("=" * 80, flush=True)
    print("UKB-PPP EXPOSURE PIPELINE", flush=True)
    print(f"Total={len(batch_df)} | completed={int(batch_df['status'].isin(COMPLETED_BATCH_STATUSES).sum())} | batch_size={args.batch_size} | cleanup={'ON' if args.delete_raw_after_processing else 'OFF'}", flush=True)
    print("=" * 80, flush=True)
    print_batch_state(batch_df, "Current batch state")
    selected_batches = batch_df
    if args.only_batch:
        wanted = {value.strip() for value in args.only_batch.split(",") if value.strip()}
        selected_batches = batch_df[batch_df["batch_id"].isin(wanted)].copy()
    if args.run and not args.rerun_completed:
        completed = selected_batches["status"].isin(COMPLETED_BATCH_STATUSES)
        skipped = int(completed.sum())
        if skipped:
            print(f"[INFO] Skipping {skipped} completed batch(es); use --rerun-completed to override", flush=True)
        selected_batches = selected_batches[~completed].copy()
    if existing_raw_base is not None and download_manifest is not None:
        selected_batches = prioritize_existing_raw_batches(selected_batches, download_manifest, existing_raw_base, pd)
        complete_raw = int((selected_batches["existing_raw_sources"] == selected_batches["source_count"]).sum())
        partial_raw = int(((selected_batches["existing_raw_sources"] > 0) & (selected_batches["existing_raw_sources"] < selected_batches["source_count"])).sum())
        print(f"[QUEUE] raw complete={complete_raw} | raw partial={partial_raw} | download required={len(selected_batches) - complete_raw}", flush=True)
    if args.max_batches is not None:
        selected_batches = selected_batches.head(args.max_batches).copy()
    if selected_batches.empty:
        raise SystemExit("[ERROR] No batches selected")
    selected_ids = selected_batches["batch_id"].tolist()
    preview = ", ".join(selected_ids[:10])
    suffix = "" if len(selected_ids) <= 10 else f", ... (+{len(selected_ids) - 10} more)"
    print(f"[INFO] Batches selected for this run: {preview}{suffix}", flush=True)
    execution_plan_path = qc_dir / "execution_plan.tsv"
    execution_plan = selected_batches.copy()
    execution_plan.insert(0, "execution_position", range(1, len(execution_plan) + 1))
    write_atomic(execution_plan, execution_plan_path)
    print(f"[INFO] Full execution plan: {execution_plan_path}", flush=True)
    selected_batches = selected_batches.drop(columns=["existing_raw_sources", "source_count"], errors="ignore")
    write_atomic(batch_df, manifest_path)
    progress_path = qc_dir / "batch_progress.tsv"
    progress_rows: list[dict[str, str]] = []

    def record_progress(batch_id: str, number: int, phase: str, detail: str = "") -> None:
        progress_rows.append({
            "timestamp": now(), "batch_id": batch_id, "batch_number": str(number),
            "batch_total": str(len(selected_batches)), "phase": phase, "detail": detail,
        })
        write_atomic(pd.DataFrame(progress_rows), progress_path)

    for position, (index, batch) in enumerate(selected_batches.iterrows(), start=1):
        batch_started = time.monotonic()
        batch_id, batch_genes = str(batch.batch_id), str(batch.genes).split(",")
        previous_status = str(batch_df.loc[index, "status"])
        print(f"\n===== Batch {position}/{len(selected_batches)}: {batch_id} ({len(batch_genes)} genes; previous={previous_status}) =====", flush=True)
        if args.run and previous_status in RESTART_CLEANUP_STATUSES:
            clean_incomplete_batch_outputs(batch_id, outdir, qc_dir, pd)
        batch_df.loc[index, "status"] = "running"
        write_atomic(batch_df, manifest_path)
        record_progress(batch_id, position, "started", f"{len(batch_genes)} genes")
        source_rows: list[dict[str, str]] = []
        if download_manifest is not None:
            selected = download_manifest[download_manifest["gene_symbol"].isin(batch_genes)]
            source_rows = selected.to_dict("records")
        else:
            for gene in batch_genes:
                for ancestry in ANCESTRIES:
                    source_rows.append({"gene_symbol": gene, "ancestry": ancestry, "source_file": "", "url": ""})

        staged_existing: dict[str, Path] = {}
        if existing_raw_base is not None:
            stage_started = time.monotonic()
            print(f"[INFO] Batch {position}/{len(selected_batches)}: checking {len(source_rows)} source archive(s) in existing raw base", flush=True)
            staged_existing = stage_existing_raw_archives(source_rows, base, existing_raw_base)
            if staged_existing:
                copied = sum(not Path(path).is_symlink() for path in staged_existing)
                copy_detail = f"; copied={copied} because symlinks are unsupported" if copied else ""
                staged_bytes = sum(Path(path).stat().st_size for path in staged_existing)
                print(f"[1/5 STAGING ] reused={len(staged_existing)} | {staged_bytes / 1e9:.1f} GB | {elapsed(stage_started)}{copy_detail}", flush=True)
            else:
                print(f"[1/5 STAGING ] reused=0 | {elapsed(stage_started)}", flush=True)
            record_progress(batch_id, position, "existing_raw_checked", f"staged={len(staged_existing)} of {len(source_rows)}")

        audit_rows = []
        raw_files_by_ancestry: dict[str, list[Path]] = {ancestry: [] for ancestry in ANCESTRIES}
        download_failed = False
        if args.run and download_manifest is not None:
            print(f"[INFO] Batch {position}/{len(selected_batches)}: looking up selected Synapse metadata", flush=True)
            try:
                source_rows = hydrate_selected_synapse_metadata(source_rows, download_manifest, Path(args.download_manifest), base)
            except (RuntimeError, ValueError, OSError) as error:
                batch_df.loc[index, "status"] = "metadata_fetch_failed"
                write_atomic(batch_df, manifest_path)
                record_progress(batch_id, position, "metadata_fetch_failed", str(error))
                print(f"[ERROR] {batch_id} Synapse metadata lookup failed: {error}", flush=True)
                if args.stop_on_error:
                    raise SystemExit(f"[ERROR] {batch_id} metadata lookup failed")
                continue
        verify_started = time.monotonic()
        print(f"[INFO] Batch {position}/{len(selected_batches)}: verifying {len(source_rows)} source archives", flush=True)
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
                ok, verification, observed_sha = verify_archive(destination, source.get("expected_size_bytes", ""), source.get("sha256", ""), source.get("md5", ""))
            else:
                matches = list((base / ancestry).glob(f"{gene}_*.tar"))
                raw_files_by_ancestry[ancestry].extend(matches)
                ok, verification, observed_sha = bool(matches), "raw_archive_found" if matches else "missing_raw_archive", ""
            audit_rows.append({"timestamp": now(), "batch_id": batch_id, "gene_symbol": gene, "ancestry": ancestry, "source_file": source["source_file"], "destination": str(destination), "action": action, "download_message": message, "verified": ok, "verification": verification, "observed_sha256": observed_sha})
            if not ok:
                download_failed = True
        download_audit = qc_dir / "downloads" / f"{batch_id}.tsv"
        write_atomic(pd.DataFrame(audit_rows), download_audit)
        verified = sum(1 for row in audit_rows if row["verified"])
        print(f"[3/5 VERIFY  ] verified={verified}/{len(audit_rows)} | invalid={len(audit_rows) - verified} | {elapsed(verify_started)}", flush=True)
        if download_failed:
            batch_df.loc[index, "status"] = "download_or_verification_failed"
            write_atomic(batch_df, manifest_path)
            record_progress(batch_id, position, "download_or_verification_failed", str(download_audit))
            if args.stop_on_error:
                raise SystemExit(f"[ERROR] {batch_id} download/verification failed; see {download_audit}")
            continue
        if not args.run or args.download_only:
            batch_df.loc[index, "status"] = "download_verified" if args.run else "planned"
            write_atomic(batch_df, manifest_path)
            record_progress(batch_id, position, str(batch_df.loc[index, "status"]), str(download_audit))
            continue

        gene_file = qc_dir / "gene_batches" / f"{batch_id}.txt"
        gene_file.parent.mkdir(parents=True, exist_ok=True)
        gene_file.write_text("\n".join(batch_genes) + "\n", encoding="utf-8")
        failed = False
        for ancestry in ANCESTRIES:
            print(f"[INFO] Batch {position}/{len(selected_batches)}: processing {ancestry}", flush=True)
            output = outdir / ancestry / f"exposure_{batch_id}.tsv"
            log = qc_dir / "processing_logs" / f"{batch_id}_{ancestry}.log"
            log.parent.mkdir(parents=True, exist_ok=True)
            command = ["Rscript", args.rscript, "--gene-file", str(gene_file), "--batch-id", batch_id, "--outdir", str(outdir / ancestry), "--batch-output", str(output), "--rawdir", str(base), "--tmpdir", args.tmpdir, "--p-threshold", str(args.p_threshold), "--ancestries", ancestry]
            result = subprocess.run(command, check=False, text=True, capture_output=True)
            log.write_text(result.stdout + result.stderr, encoding="utf-8")
            if result.returncode or not output.exists():
                failed = True
                print(f"[ERROR] {batch_id} {ancestry}; see {log}", flush=True)
        if failed:
            batch_df.loc[index, "status"] = "processing_failed"
            batch_df.loc[index, "raw_cleanup"] = "retained_processing_failed"
        elif args.delete_raw_after_processing:
            cleanup_rows = []
            for ancestry in ANCESTRIES:
                status_path = outdir / ancestry / "logs" / f"{batch_id}_gene_status.tsv"
                for row in remove_raw_files(raw_files_by_ancestry[ancestry], status_path):
                    cleanup_rows.append({"batch_id": batch_id, "ancestry": ancestry, **row})
            if args.delete_existing_raw_after_processing:
                for row in list(cleanup_rows):
                    original = staged_existing.get(row["raw_file"])
                    if row["action"] != "deleted" or original is None:
                        continue
                    try:
                        original.unlink()
                        cleanup_rows.append({
                            "batch_id": batch_id, "ancestry": row["ancestry"], "raw_file": str(original), "action": "deleted_existing_source",
                            "reason": "staged_source_processed",
                        })
                    except OSError as error:
                        cleanup_rows.append({
                            "batch_id": batch_id, "ancestry": row["ancestry"], "raw_file": str(original), "action": "delete_existing_source_failed",
                            "reason": str(error),
                        })
            cleanup_path = qc_dir / "raw_cleanup" / f"{batch_id}.tsv"
            write_atomic(pd.DataFrame(cleanup_rows), cleanup_path)
            record_raw_cleanup_in_manifest(download_manifest, cleanup_rows, batch_id)
            write_atomic(download_manifest, Path(args.download_manifest))
            cleanup_failed = any(row["action"] not in {"deleted", "deleted_existing_source"} for row in cleanup_rows)
            batch_df.loc[index, "status"] = "completed_raw_retained" if cleanup_failed else "completed_raw_deleted"
            batch_df.loc[index, "raw_cleanup"] = str(cleanup_path)
        else:
            batch_df.loc[index, "status"] = "completed"
            batch_df.loc[index, "raw_cleanup"] = "retained_flag_not_set"
        write_atomic(batch_df, manifest_path)
        record_progress(batch_id, position, str(batch_df.loc[index, "status"]), str(batch_df.loc[index, "raw_cleanup"]))
        print(f"[DONE] {position}/{len(selected_batches)} {batch_id} | status={batch_df.loc[index, 'status']} | elapsed={elapsed(batch_started)}", flush=True)
        if failed and args.stop_on_error:
            raise SystemExit(f"[ERROR] {batch_id} processing failed")

    print_batch_state(batch_df, "Final batch state")
    print(f"[INFO] Batch manifest: {manifest_path}", flush=True)
    print(f"[INFO] Batch progress: {progress_path}", flush=True)


if __name__ == "__main__":
    main()
