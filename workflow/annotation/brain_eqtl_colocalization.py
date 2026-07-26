#!/usr/bin/env python3
"""Harmonize GIGASTROKE regional statistics with brain eQTLs and run coloc ABF.

The program is install-free and performs no downloads.  Inputs must contain full
regional statistics (not only significant QTLs).  Unavailable controlled data
produce an explicit NOT_RUN status rather than evidence of absence.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

FIELDS = ["locus_id", "gene", "bulk_tissue", "cell_type_or_subtype", "qtl_dataset",
          "evidence_family", "common_snps", "pp4", "pp4_over_pp3_pp4", "prior_p1",
          "prior_p2", "prior_p12", "coverage_status", "coverage_fraction", "status", "reason"]


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def resolve(root: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (root / path).resolve()


def variant_key(row: dict[str, str]) -> tuple[str, int, str, str]:
    chrom = row["chromosome"].removeprefix("chr").upper()
    return chrom, int(row["position"]), row["ref"].upper(), row["alt"].upper()


def log_abf(beta: float, se: float, prior_variance: float) -> float:
    variance = se * se
    shrink = prior_variance / (variance + prior_variance)
    return .5 * (math.log1p(-shrink) + shrink * (beta / se) ** 2)


def logsum(values: list[float]) -> float:
    if not values:
        return -math.inf
    top = max(values)
    return top + math.log(sum(math.exp(value - top) for value in values))


def coloc(gwas: list[tuple[dict, dict]], priors: dict, variances: dict) -> tuple[float, float]:
    a = [log_abf(float(x[0]["beta"]), float(x[0]["se"]), variances["gwas"]) for x in gwas]
    b = [log_abf(float(x[1]["beta"]), float(x[1]["se"]), variances["eqtl"]) for x in gwas]
    l1, l2, l12 = logsum(a), logsum(b), logsum([x + y for x, y in zip(a, b)])
    # H3 is two distinct causal variants; subtract same-variant pairs in linear space.
    distinct = math.exp(l1 + l2) - math.exp(l12)
    logs = [0.0, math.log(priors["p1"]) + l1, math.log(priors["p2"]) + l2,
            math.log(priors["p1"] * priors["p2"]) + math.log(max(distinct, 1e-300)),
            math.log(priors["p12"]) + l12]
    denominator = logsum(logs)
    pp3, pp4 = math.exp(logs[3] - denominator), math.exp(logs[4] - denominator)
    return pp4, pp4 / (pp3 + pp4) if pp3 + pp4 else 0.0


def harmonize(gwas_rows: list[dict], eqtl_rows: list[dict]) -> list[tuple[dict, dict]]:
    indexed = {variant_key(row): row for row in gwas_rows}
    result = []
    for eqtl in eqtl_rows:
        key = variant_key(eqtl)
        gwas = indexed.get(key)
        if gwas is not None:
            result.append((gwas, eqtl))
            continue
        swapped = (key[0], key[1], key[3], key[2])
        gwas = indexed.get(swapped)
        if gwas is not None:
            adjusted = dict(eqtl)
            adjusted["beta"] = str(-float(eqtl["beta"]))
            result.append((gwas, adjusted))
    return result


def run(config_path: Path, output_override: Path | None = None) -> None:
    config = json.loads(config_path.read_text(encoding="utf-8")); root = config_path.parent
    out = output_override or resolve(root, config["output"]); out.parent.mkdir(parents=True, exist_ok=True)
    loci = rows(resolve(root, config["gigastroke_loci"])); priors = config["priors"]
    variances = config.get("prior_variances", {"gwas": .04, "eqtl": .0225})
    results = []
    for dataset in config["datasets"]:
        registry = json.loads(resolve(root, dataset["registry"]).read_text(encoding="utf-8"))
        resolution = registry["resolution"]
        unavailable = dataset.get("availability") in {"CONTROLLED_ACCESS_REQUIRED", "UNAVAILABLE"} or not dataset.get("input")
        locus_builds = {locus.get("genome_build", config.get("genome_build", "")) for locus in loci}
        if not unavailable and locus_builds != {registry["genome_build"]}:
            raise SystemExit(f"genome build mismatch for {registry['dataset_id']}: loci={sorted(locus_builds)}, eQTL={registry['genome_build']}")
        for locus in loci:
            base = {"locus_id": locus["locus_id"], "gene": dataset.get("gene", locus.get("gene", "")),
                    "bulk_tissue": dataset.get("tissue_or_cell_type", "") if resolution == "bulk_tissue" else "",
                    "cell_type_or_subtype": dataset.get("tissue_or_cell_type", "") if resolution != "bulk_tissue" else "",
                    "qtl_dataset": registry["dataset_id"], "evidence_family": registry["evidence_family"],
                    "prior_p1": priors["p1"], "prior_p2": priors["p2"], "prior_p12": priors["p12"]}
            if unavailable:
                status = "NOT_RUN_ACCESS_REQUIRED" if dataset.get("availability") == "CONTROLLED_ACCESS_REQUIRED" else "NOT_RUN_INPUT_UNAVAILABLE"
                results.append({**base, "common_snps": 0, "pp4": "", "pp4_over_pp3_pp4": "",
                                "coverage_status": "UNAVAILABLE", "coverage_fraction": 0, "status": status,
                                "reason": dataset.get("reason", "regional summary statistics were not supplied")})
                continue
            gwas_all = rows(resolve(root, locus["regional_statistics"])); eqtl_all = rows(resolve(root, dataset["input"]))
            chrom, start, end = locus["chromosome"].removeprefix("chr"), int(locus["start"]), int(locus["end"])
            gwas = [r for r in gwas_all if r["chromosome"].removeprefix("chr") == chrom and start <= int(r["position"]) <= end]
            eqtl = [r for r in eqtl_all if r.get("gene", base["gene"]) == base["gene"] and r["chromosome"].removeprefix("chr") == chrom and start <= int(r["position"]) <= end]
            common = harmonize(gwas, eqtl); fraction = len(common) / len(gwas) if gwas else 0
            minimum = int(config.get("minimum_common_snps", 50)); threshold = float(config.get("minimum_coverage", .8))
            coverage = "COMPLETE" if fraction >= threshold else "PARTIAL"
            if len(common) < minimum:
                results.append({**base, "common_snps": len(common), "pp4": "", "pp4_over_pp3_pp4": "",
                                "coverage_status": coverage, "coverage_fraction": f"{fraction:.6g}",
                                "status": "NOT_RUN_INSUFFICIENT_OVERLAP", "reason": f"requires at least {minimum} common SNPs"})
            else:
                pp4, ratio = coloc(common, priors, variances)
                results.append({**base, "common_snps": len(common), "pp4": f"{pp4:.8g}",
                                "pp4_over_pp3_pp4": f"{ratio:.8g}", "coverage_status": coverage,
                                "coverage_fraction": f"{fraction:.6g}", "status": "SUCCESS", "reason": ""})
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=FIELDS, lineterminator="\n"); writer.writeheader(); writer.writerows(results)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", type=Path); args = parser.parse_args(); run(args.config, args.output)
