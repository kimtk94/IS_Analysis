#!/usr/bin/env python3
"""Estimate a batch size likely to contain at least one gene with instruments."""
from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("status_files", nargs="+", help="Per-archive *_gene_status.tsv files or directories containing them.")
    parser.add_argument("--target-probability", type=float, default=0.95, help="Desired probability of >=1 successful gene (default: 0.95).")
    parser.add_argument("--ancestry-rule", choices=("any", "all"), default="any", help="Count a gene when any observed ancestry succeeds, or require all observed ancestries.")
    parser.add_argument("--conservative", action="store_true", help="Use the 95%% Wilson lower bound instead of the observed success rate.")
    return parser.parse_args()


def discover(inputs: list[str]) -> list[Path]:
    paths: set[Path] = set()
    for value in inputs:
        path = Path(value)
        if path.is_dir():
            paths.update(path.rglob("*_gene_status.tsv"))
        elif path.is_file():
            paths.add(path)
        else:
            raise SystemExit(f"[ERROR] Status path not found: {path}")
    if not paths:
        raise SystemExit("[ERROR] No *_gene_status.tsv files found")
    return sorted(paths)


def positive_rows(value: str) -> bool:
    try:
        return int(float(value)) > 0
    except (TypeError, ValueError):
        return False


def gene_outcomes(paths: list[Path], ancestry_rule: str) -> dict[str, bool]:
    outcomes: dict[str, dict[str, bool]] = defaultdict(dict)
    for path in paths:
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                gene = row.get("gene_symbol", "").strip().upper()
                ancestry = row.get("ancestry", "").strip().upper()
                if not gene or not ancestry or row.get("standardization_status") != "completed":
                    continue
                success = row.get("instrument_selection_status") == "completed" and positive_rows(row.get("n_filtered_rows", ""))
                outcomes[gene][ancestry] = outcomes[gene].get(ancestry, False) or success
    return {
        gene: (any(values.values()) if ancestry_rule == "any" else all(values.values()))
        for gene, values in outcomes.items() if values
    }


def wilson_lower(successes: int, total: int, z: float = 1.959963984540054) -> float:
    proportion = successes / total
    denominator = 1 + z * z / total
    center = proportion + z * z / (2 * total)
    margin = z * math.sqrt(proportion * (1 - proportion) / total + z * z / (4 * total * total))
    return max(0.0, (center - margin) / denominator)


def required_size(rate: float, target: float) -> int:
    return math.ceil(math.log1p(-target) / math.log1p(-rate))


def main() -> None:
    args = parse_args()
    if not 0 < args.target_probability < 1:
        raise SystemExit("[ERROR] --target-probability must be between 0 and 1")
    outcomes = gene_outcomes(discover(args.status_files), args.ancestry_rule)
    total, successes = len(outcomes), sum(outcomes.values())
    if not total:
        raise SystemExit("[ERROR] No genes with completed standardization were found")
    observed = successes / total
    rate = wilson_lower(successes, total) if args.conservative else observed
    print(f"genes_evaluated\t{total}")
    print(f"genes_with_instruments\t{successes}")
    print(f"observed_gene_success_rate\t{observed:.6f}")
    print(f"rate_used\t{rate:.6f}")
    print(f"target_probability\t{args.target_probability:.6f}")
    if rate <= 0:
        print("recommended_batch_size\tNA")
        raise SystemExit("[ERROR] No positive success-rate estimate; run more pilot genes")
    print(f"recommended_batch_size\t{required_size(rate, args.target_probability)}")
    print("note\tProbabilistic estimate only; it cannot guarantee a successful gene in every batch")


if __name__ == "__main__":
    main()
