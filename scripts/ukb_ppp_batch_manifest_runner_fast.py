#!/usr/bin/env python3
"""Build and optionally run UKB-PPP paired EUR/EAS exposure batches.

The runner scans ancestry-aware rawdata folders, writes the current EUR/EAS
complete-pair gene list, creates fixed-size gene batch files, maintains an
atomic batch manifest, and can invoke the R exposure preparation script batch by
batch. It is intended for legacy recovery/Colab workflows; v2 work should move
this behavior behind explicit pipeline stages and output contracts.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import re
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Any


def now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def shash(value: str, n: int = 8) -> str:
    return hashlib.md5(value.encode()).hexdigest()[:n]


def safe_gene(gene: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", gene)


def write_atomic(df: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    df.to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, path)


def scan_valid(base: Path, pd: Any) -> Any:
    rows = []
    for ancestry in ["EUR", "EAS"]:
        folder = base / ancestry
        if not folder.exists():
            continue
        for path in sorted(folder.iterdir(), key=lambda item: item.name):
            if not path.is_file():
                continue
            name = path.name
            size = path.stat().st_size
            is_temp = ".synapse_download_" in name
            is_tar = name.endswith(".tar")
            gene = name.split("_")[0].upper() if "_" in name else ""
            rows.append({
                "ancestry": ancestry,
                "gene_symbol": gene,
                "file_name": name,
                "file_path": str(path),
                "size_bytes": size,
                "is_valid_tar": bool(is_tar and not is_temp and size > 0),
                "is_temp_synapse": is_temp,
            })
    return pd.DataFrame(rows)


def batch_fname(batch_id: str, genes: list[str], max_len: int = 180) -> str:
    label = "_".join(safe_gene(gene) for gene in genes)
    filename = f"{batch_id}_{label}.txt"
    if len(filename) > max_len:
        filename = f"{batch_id}_{genes[0]}_to_{genes[-1]}_{shash(label)}.txt"
    return filename


def load_manifest(path: Path, pd: Any) -> Any:
    if path.exists() and path.stat().st_size > 0:
        return pd.read_csv(path, sep="\t")
    return pd.DataFrame()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build and optionally run UKB-PPP paired EUR/EAS exposure batches.")
    parser.add_argument("--base", default="data/rawdata/pqtl/selected_targets")
    parser.add_argument("--outdir", default="results/qc/pair_priority")
    parser.add_argument("--batch-dir", default="results/qc/pair_priority/gene_batches_4")
    parser.add_argument("--manifest", default="results/qc/pair_priority/ukb_ppp_pair_batch_manifest.tsv")
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument("--rscript", default="scripts/01_prepare_exposure_fast.R")
    parser.add_argument("--exposure-outdir", default="results/exposure_batches")
    parser.add_argument("--tmpdir", default="/content/ukbppp_tmp")
    parser.add_argument("--p-threshold", default="5e-8")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--rewrite-batch-files", action="store_true")
    parser.add_argument("--stop-on-error", action="store_true")
    parser.add_argument("--force-rerun-completed", action="store_true")
    parser.add_argument("--no-eur-first", action="store_true")
    parser.add_argument("--no-copy-to-local", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    import pandas as pd

    base = Path(args.base)
    outdir = Path(args.outdir)
    batch_dir = Path(args.batch_dir)
    manifest_path = Path(args.manifest)
    exposure_outdir = Path(args.exposure_outdir)

    outdir.mkdir(parents=True, exist_ok=True)
    batch_dir.mkdir(parents=True, exist_ok=True)
    exposure_outdir.mkdir(parents=True, exist_ok=True)
    (exposure_outdir / "logs").mkdir(parents=True, exist_ok=True)

    print(f"[{now()}] scan rawdata: {base}")
    audit = scan_valid(base, pd)
    audit.to_csv(outdir / "ukb_ppp_pair_current_file_audit.tsv", sep="\t", index=False)

    valid = audit[audit["is_valid_tar"]].copy() if not audit.empty else audit
    genes_by = valid.groupby("ancestry")["gene_symbol"].apply(lambda values: set(values.dropna().astype(str))).to_dict() if not valid.empty else {}
    eur = genes_by.get("EUR", set())
    eas = genes_by.get("EAS", set())
    paired_genes = sorted(eur & eas)

    (outdir / "ukb_ppp_complete_pair_genes.current.txt").write_text("\n".join(paired_genes) + "\n", encoding="utf-8")
    n_batches = (len(paired_genes) + args.batch_size - 1) // args.batch_size
    summary = pd.DataFrame([{
        "timestamp": now(),
        "valid_EUR_genes": len(eur),
        "valid_EAS_genes": len(eas),
        "complete_pair_genes": len(paired_genes),
        "batch_size": args.batch_size,
        "n_batches": n_batches,
    }])
    summary.to_csv(outdir / "ukb_ppp_pair_batch_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))

    old = load_manifest(manifest_path, pd)
    old_map = {str(row.batch_id): row for _, row in old.iterrows()} if not old.empty and "batch_id" in old.columns else {}

    rows = []
    for offset in range(0, len(paired_genes), args.batch_size):
        genes = paired_genes[offset : offset + args.batch_size]
        batch_id = f"batch_{offset // args.batch_size + 1:03d}"
        batch_file = batch_dir / batch_fname(batch_id, genes)
        expected_output = exposure_outdir / f"exposure_{batch_id}.tsv"
        log_file = exposure_outdir / "logs" / f"{batch_id}.runner.log"

        old_row = old_map.get(batch_id)
        status = "pending"
        started = ""
        completed = ""
        returncode = ""
        message = ""
        if old_row is not None:
            status = str(getattr(old_row, "status", "pending"))
            started = str(getattr(old_row, "started_at", ""))
            completed = str(getattr(old_row, "completed_at", ""))
            returncode = str(getattr(old_row, "returncode", ""))
            message = str(getattr(old_row, "message", ""))

        output_exists = expected_output.exists() and expected_output.stat().st_size > 0
        if output_exists and status in {"pending", "running", "failed", "output_missing", "completed", "completed_existing"}:
            status = "completed_existing"

        rows.append({
            "batch_id": batch_id,
            "n_genes": len(genes),
            "genes": ",".join(genes),
            "batch_file": str(batch_file),
            "batch_file_ok": batch_file.exists(),
            "expected_output": str(expected_output),
            "output_exists": output_exists,
            "status": status,
            "started_at": started,
            "completed_at": completed,
            "returncode": returncode,
            "log_file": str(log_file),
            "message": message,
        })

    manifest = pd.DataFrame(rows)
    for _, row in manifest.iterrows():
        batch_file = Path(row.batch_file)
        if args.rewrite_batch_files or not batch_file.exists():
            batch_file.write_text("\n".join(str(row.genes).split(",")) + "\n", encoding="utf-8")
    if not manifest.empty:
        manifest["batch_file_ok"] = manifest["batch_file"].map(lambda value: Path(value).exists())
    write_atomic(manifest, manifest_path)
    print(f"[{now()}] manifest updated: {manifest_path}")

    if not args.run:
        return

    for idx, row in manifest.iterrows():
        batch_id = str(row.batch_id)
        expected_output = Path(row.expected_output)
        output_exists = expected_output.exists() and expected_output.stat().st_size > 0
        if not args.force_rerun_completed and str(row.status) in {"completed", "completed_existing"} and output_exists:
            print(f"[SKIP] {batch_id}")
            continue

        batch_file = Path(row.batch_file)
        if not batch_file.exists():
            manifest.loc[idx, "status"] = "batch_file_error"
            manifest.loc[idx, "message"] = "Batch file missing"
            write_atomic(manifest, manifest_path)
            if args.stop_on_error:
                raise SystemExit(f"[ERROR] missing {batch_file}")
            continue

        cmd = [
            "Rscript", args.rscript,
            "--gene-file", str(batch_file),
            "--batch-id", batch_id,
            "--outdir", str(exposure_outdir),
            "--rawdir", str(base),
            "--tmpdir", args.tmpdir,
            "--p-threshold", str(args.p_threshold),
        ]
        if not args.no_eur_first:
            cmd.append("--eur-first")
        if not args.no_copy_to_local:
            cmd.append("--copy-to-local")

        log_file = Path(row.log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)
        manifest.loc[idx, "status"] = "running"
        manifest.loc[idx, "started_at"] = now()
        manifest.loc[idx, "message"] = "Started"
        write_atomic(manifest, manifest_path)

        print("[RUN]", batch_id, " ".join(cmd))
        with log_file.open("w", encoding="utf-8") as handle:
            result = subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT, text=True, check=False)

        output_exists = expected_output.exists() and expected_output.stat().st_size > 0
        manifest.loc[idx, "returncode"] = result.returncode
        manifest.loc[idx, "completed_at"] = now()
        manifest.loc[idx, "output_exists"] = output_exists

        if result.returncode != 0:
            manifest.loc[idx, "status"] = "failed"
            manifest.loc[idx, "message"] = f"Rscript failed. See {log_file}"
            write_atomic(manifest, manifest_path)
            print("[FAILED]", batch_id)
            if args.stop_on_error:
                raise SystemExit(1)
        else:
            if output_exists:
                manifest.loc[idx, "status"] = "completed"
                manifest.loc[idx, "message"] = "Completed"
                print("[DONE]", batch_id)
            else:
                manifest.loc[idx, "status"] = "output_missing"
                manifest.loc[idx, "message"] = f"Output missing: {expected_output}"
                print("[WARN]", batch_id, "output missing")
                write_atomic(manifest, manifest_path)
                if args.stop_on_error:
                    raise SystemExit(1)
            write_atomic(manifest, manifest_path)

    print(f"[{now()}] done. manifest={manifest_path}")


if __name__ == "__main__":
    main()
