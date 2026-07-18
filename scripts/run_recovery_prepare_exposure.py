#!/usr/bin/env python3
"""Run exact-source recovery batches for UKB-PPP exposure preparation.

The recovery input should contain the audit columns produced by the exposure
batch recovery audit, including ``batch_id``, ``gene``, ``source_file_key``,
``observed_output_bool``, and ``source_process_status``. The runner writes one
gene list and one exact source-file list per batch, then invokes the fast R
preparation script with ``--source-file-list`` so only missing source files are
processed.
"""
from __future__ import annotations

import argparse
import subprocess
from datetime import datetime
from pathlib import Path



def latest_file(folder: Path, pattern: str) -> Path:
    files = sorted(folder.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files matched: {folder}/{pattern}")
    return files[-1]


def read_existing_source_files(path: Path) -> set[str]:
    """Return source_file values already present in an existing recovery output."""
    if not path.exists() or path.stat().st_size == 0:
        return set()
    try:
        header = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
        if "source_file" not in header:
            return set()
        source_files: set[str] = set()
        for chunk in pd.read_csv(path, sep="\t", dtype=str, usecols=["source_file"], chunksize=200_000):
            chunk = chunk.fillna("")
            source_files |= set(chunk["source_file"].astype(str))
        source_files.discard("")
        return source_files
    except Exception as exc:  # defensive: legacy recovery should continue past malformed outputs
        print(f"[WARN] failed reading existing recovery output: {path} | {exc}")
        return set()


def write_lines(path: Path, values: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(values) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run exact recovery prepare-exposure batches.")
    parser.add_argument("--root", default="/content/drive/MyDrive/IS_Analysis")
    parser.add_argument("--recovery-sources", default=None)
    parser.add_argument("--rscript", default="01_prepare_exposure_fast.R")
    parser.add_argument("--rawdir", default="data/rawdata/pqtl/selected_targets")
    parser.add_argument("--outdir", default="results/exposure_batches_recovery_exact")
    parser.add_argument("--tmpdir", default="/content/ukbppp_tmp")
    parser.add_argument("--work-batch-dir", default="results/qc/recovery_gene_batches_exact")
    parser.add_argument("--report-dir", default="results/qc/recovery_run_reports")
    parser.add_argument("--ancestries", default="EUR")
    parser.add_argument("--p-threshold", default="5e-8")
    parser.add_argument("--max-batches", type=int, default=None)
    parser.add_argument("--only-batch", default=None)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--stop-on-error", action="store_true")
    parser.add_argument("--copy-to-local", action="store_true", default=True)
    parser.add_argument("--no-copy-to-local", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    import pandas as pd

    root = Path(args.root)

    if args.recovery_sources:
        recovery_path = root / args.recovery_sources
        if not recovery_path.exists():
            recovery_path = Path(args.recovery_sources)
    else:
        recovery_path = latest_file(
            root / "results/qc/exposure_batch_recovery_audit",
            "recovery_needed_sources_*.tsv",
        )

    if not recovery_path.exists():
        raise FileNotFoundError(recovery_path)

    print("[INFO] recovery sources:", recovery_path)
    rec = pd.read_csv(recovery_path, sep="\t", dtype=str).fillna("")

    required = ["batch_id", "gene", "source_file_key", "observed_output_bool", "source_process_status"]
    missing = [column for column in required if column not in rec.columns]
    if missing:
        raise ValueError(f"Missing columns in recovery sources: {missing}")

    rec = rec[rec["observed_output_bool"].astype(str).str.upper().isin(["FALSE", "0", "NO", "N", ""])].copy()
    rec = rec[rec["source_process_status"].isin(["not_attempted_raw_available", "prepare_failed__standardize_failed"])].copy()

    if args.only_batch:
        keep = {value.strip() for value in args.only_batch.split(",") if value.strip()}
        rec = rec[rec["batch_id"].isin(keep)].copy()

    if rec.empty:
        print("[INFO] No exact recovery targets selected.")
        return

    batch_ids = sorted(rec["batch_id"].dropna().unique())
    if args.max_batches is not None:
        batch_ids = batch_ids[: args.max_batches]
        rec = rec[rec["batch_id"].isin(batch_ids)].copy()

    work_dir = root / args.work_batch_dir
    outdir = root / args.outdir
    report_dir = root / args.report_dir
    work_dir.mkdir(parents=True, exist_ok=True)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "logs").mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)

    print("[INFO] exact recovery source files:", len(rec))
    print("[INFO] exact recovery batches:", len(batch_ids))
    print("[INFO] output dir:", outdir)

    rows = []
    for index, batch_id in enumerate(batch_ids, start=1):
        sub = rec[rec["batch_id"] == batch_id].copy()
        genes = sorted({value for value in sub["gene"].astype(str).str.upper() if value})
        source_files = sorted({value for value in sub["source_file_key"].astype(str) if value})

        expected_output = outdir / f"exposure_{batch_id}.tsv"
        existing_recovery_sources = read_existing_source_files(expected_output)
        original_source_n = len(source_files)
        source_files = [value for value in source_files if value not in existing_recovery_sources]

        if len(source_files) == 0:
            print(f"\n===== [{index}/{len(batch_ids)}] {batch_id} =====")
            print("[SKIP] all recovery source_files already exist in recovery output")
            rows.append({
                "batch_id": batch_id,
                "action": "skip_existing_recovery_output_complete",
                "n_genes": len(genes),
                "n_source_files": original_source_n,
                "n_source_files_to_run": 0,
                "gene_file": "",
                "source_file_list": "",
                "expected_output": str(expected_output.relative_to(root)),
                "returncode": 0,
                "output_exists": expected_output.exists(),
                "message": "All target source_files already recovered.",
            })
            continue

        sub = sub[sub["source_file_key"].astype(str).isin(source_files)].copy()
        genes = sorted({value for value in sub["gene"].astype(str).str.upper() if value})

        gene_file = work_dir / f"{batch_id}_recovery_genes.txt"
        source_file_list = work_dir / f"{batch_id}_recovery_source_files.txt"
        write_lines(gene_file, genes)
        write_lines(source_file_list, source_files)

        log_file = outdir / "logs" / f"{batch_id}.recovery_exact.log"
        cmd = [
            "Rscript", str(root / args.rscript),
            "--gene-file", str(gene_file),
            "--source-file-list", str(source_file_list),
            "--batch-id", batch_id,
            "--outdir", str(outdir),
            "--rawdir", str(root / args.rawdir),
            "--tmpdir", args.tmpdir,
            "--ancestries", args.ancestries,
            "--p-threshold", str(args.p_threshold),
        ]
        if not args.no_copy_to_local:
            cmd.append("--copy-to-local")

        print(f"\n===== [{index}/{len(batch_ids)}] {batch_id} =====")
        print("[INFO] recovery genes:", ",".join(genes))
        print("[INFO] exact source files:", len(source_files))
        print("[INFO] first source files:", ",".join(source_files[:5]))
        print("[CMD]", " ".join(cmd))

        if args.dry_run or not args.run:
            rows.append({
                "batch_id": batch_id,
                "action": "dry_run",
                "n_genes": len(genes),
                "n_source_files": original_source_n,
                "n_source_files_to_run": len(source_files),
                "gene_file": str(gene_file.relative_to(root)),
                "source_file_list": str(source_file_list.relative_to(root)),
                "expected_output": str(expected_output.relative_to(root)),
                "returncode": "",
                "output_exists": expected_output.exists(),
                "message": "Dry-run only",
            })
            continue

        with log_file.open("w", encoding="utf-8") as handle:
            res = subprocess.run(cmd, cwd=str(root), stdout=handle, stderr=subprocess.STDOUT, text=True, check=False)

        output_exists = expected_output.exists() and expected_output.stat().st_size > 0
        if res.returncode == 0 and output_exists:
            action = "success"
            message = "completed"
        elif res.returncode == 0 and not output_exists:
            action = "output_missing"
            message = "Rscript returned 0 but output missing"
        else:
            action = "failed"
            message = f"Rscript failed. See {log_file}"

        print("[RESULT]", batch_id, action)
        rows.append({
            "batch_id": batch_id,
            "action": action,
            "n_genes": len(genes),
            "n_source_files": original_source_n,
            "n_source_files_to_run": len(source_files),
            "gene_file": str(gene_file.relative_to(root)),
            "source_file_list": str(source_file_list.relative_to(root)),
            "expected_output": str(expected_output.relative_to(root)),
            "returncode": res.returncode,
            "output_exists": output_exists,
            "log_file": str(log_file.relative_to(root)),
            "message": message,
        })

        if action != "success" and args.stop_on_error:
            break

    report = pd.DataFrame(rows)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_path = report_dir / f"recovery_prepare_exposure_exact_report_{timestamp}.tsv"
    report.to_csv(report_path, sep="\t", index=False)

    print("\n===== RECOVERY EXACT REPORT =====")
    print(report.to_string(index=False))
    print("[OK] report:", report_path)


if __name__ == "__main__":
    main()
