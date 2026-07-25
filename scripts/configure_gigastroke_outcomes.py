#!/usr/bin/env python3
"""Materialize a runnable real-data config from downloaded GIGASTROKE files."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root", type=Path,
        default=Path(os.environ.get("PROJECT_ROOT", "/content/drive/MyDrive/IS_Analysis_V2")),
        help="Project root (default: PROJECT_ROOT or /content/drive/MyDrive/IS_Analysis_V2).",
    )
    parser.add_argument("--chain", required=True, type=Path)
    parser.add_argument("--source-reference", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--template", type=Path, default=Path("config/gigastroke_outcomes.example.json"))
    parser.add_argument("--output", type=Path, default=Path("config/gigastroke_outcomes.json"))
    parser.add_argument("--eur-discovery", default="gigastroke_is_EUR")
    parser.add_argument("--eas-outcomes", nargs="+", default=[
        "gigastroke_is_EAS", "gigastroke_las_EAS", "gigastroke_ces_EAS", "gigastroke_svs_EAS"])
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    raw = project_root / "data/rawdata/outcome/gigastroke_gwas_catalog"
    chain, source_reference, reference = args.chain.resolve(), args.source_reference.resolve(), args.reference.resolve()
    for label, path in (("raw-data directory", raw), ("chain", chain),
                        ("source reference", source_reference), ("target reference", reference)):
        if not path.exists():
            raise SystemExit(f"{label} does not exist: {path}")

    config = json.loads(args.template.read_text(encoding="utf-8"))
    config["output_dir"] = str(project_root / "data/standardized/outcome/gigastroke")
    config["liftover"] = {
        "chain": str(chain), "chain_sha256": digest(chain),
        "source_reference": str(source_reference), "source_reference_sha256": digest(source_reference),
        "target_reference": str(reference), "target_reference_sha256": digest(reference),
    }
    config["selection"] = {
        "eur_discovery": args.eur_discovery,
        "eas_replication_subtypes": args.eas_outcomes,
    }
    for dataset in config["datasets"]:
        # The downloader preserves this GCST + phenotype + ancestry stem but the
        # server controls whether the final suffix is .tsv.gz, .h.tsv.gz, etc.
        stem = Path(dataset["input"]).name.split(".tsv.gz", 1)[0]
        matches = sorted(raw.glob(stem + "*"))
        if len(matches) != 1:
            raise SystemExit(f"expected exactly one downloaded file for {dataset['id']} ({stem}*); found {len(matches)}")
        dataset["input"] = str(matches[0].resolve())

    configured = {item["id"] for item in config["datasets"]}
    selected = {args.eur_discovery, *args.eas_outcomes}
    if selected - configured:
        raise SystemExit("unknown selected dataset IDs: " + ", ".join(sorted(selected - configured)))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    print(f"[OK] wrote {args.output} with verified local paths and SHA-256 values")


if __name__ == "__main__":
    main()
