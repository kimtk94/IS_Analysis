#!/usr/bin/env python3
"""Ancestry-specific LD clumping and single-causal fine-mapping stage.

The stage is intentionally install-free.  It consumes locus-wide summary
statistics and precomputed, ancestry-matched LD assets; it never queries or
substitutes a reference panel.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def resolve(root: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (root / path).resolve()


def positive_definite(matrix: list[list[float]], tolerance: float) -> bool:
    """Cholesky test without numpy; tolerance permits harmless rounding."""
    n = len(matrix)
    lower = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            value = matrix[i][j] - sum(lower[i][k] * lower[j][k] for k in range(j))
            if i == j:
                if value <= tolerance:
                    return False
                lower[i][j] = math.sqrt(value)
            else:
                lower[i][j] = value / lower[j][j]
    return True


def load_ld(matrix_path: Path, variants_path: Path, tolerance: float):
    variants = read_tsv(variants_path)
    metadata = {row["variant_id"]: row for row in variants}
    with matrix_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)[1:]
        rows, matrix = [], []
        for row in reader:
            rows.append(row[0]); matrix.append([float(x) for x in row[1:]])
    if rows != header or len(matrix) != len(header) or any(len(x) != len(header) for x in matrix):
        raise ValueError("LD matrix row/column variant IDs are inconsistent")
    if set(header) != set(metadata):
        raise ValueError("LD variant metadata does not exactly cover the matrix")
    symmetric = all(abs(matrix[i][j] - matrix[j][i]) <= tolerance for i in range(len(matrix)) for j in range(len(matrix)))
    diagonal = all(abs(matrix[i][i] - 1.0) <= tolerance for i in range(len(matrix)))
    return header, matrix, metadata, symmetric and diagonal and positive_definite(matrix, tolerance)


def alignment(effect: str, other: str, ref: str, alt: str):
    complement = str.maketrans("ACGT", "TGCA")
    pairs = {(alt, ref): ("DIRECT", 1), (ref, alt): ("SWAP", -1),
             (alt.translate(complement), ref.translate(complement)): ("COMPLEMENT", 1),
             (ref.translate(complement), alt.translate(complement)): ("SWAP_COMPLEMENT", -1)}
    return pairs.get((effect.upper(), other.upper()))


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def run(config_path: Path, output_override: Path | None = None) -> None:
    config = json.loads(config_path.read_text(encoding="utf-8")); root = config_path.parent
    panels = config["ld_reference_panels"]
    for ancestry in ("EUR", "EAS"):
        panel = panels.get(ancestry, {})
        required = {"name", "genome_build", "ancestry", "sample_provenance", "matrix", "variants"}
        if required - panel.keys() or panel.get("ancestry") != ancestry:
            raise SystemExit(f"invalid or incomplete {ancestry} LD reference contract")
    summaries = read_tsv(resolve(root, config["full_summary_statistics"]))
    loci = {row["locus"] for row in summaries}
    if len(loci) != 1:
        raise SystemExit("each stage unit must contain full summary statistics for exactly one locus")
    outdir = output_override.resolve() if output_override else resolve(root, config["output_dir"])
    outdir.mkdir(parents=True, exist_ok=True)
    method = config["methods"]; gates_cfg = config["gates"]
    results, instruments, leads, gates = [], [], [], []
    per_ancestry = {}
    for ancestry, panel in panels.items():
        ids, ld, meta, pd = load_ld(resolve(root, panel["matrix"]), resolve(root, panel["variants"]), gates_cfg["matrix_tolerance"])
        index = {variant: i for i, variant in enumerate(ids)}
        rows = [row for row in summaries if row["ancestry"] == ancestry]
        if not rows:
            continue
        builds_ok = all(row["genome_build"] == panel["genome_build"] for row in rows)
        aligned, failures = [], 0
        for row in rows:
            if row["variant_id"] not in meta:
                continue
            match = alignment(row["effect_allele"], row["other_allele"], meta[row["variant_id"]]["ref"], meta[row["variant_id"]]["alt"])
            if not match:
                failures += 1; continue
            item = dict(row); item["z"] = float(row["beta"]) / float(row["se"]) * match[1]; item["alignment"] = match[0]
            aligned.append(item)
        missing_rate = (len(rows) - sum(row["variant_id"] in index for row in rows)) / len(rows)
        coverage = len(aligned) / len(rows)
        min_n = min(float(row["effective_sample_size"]) for row in aligned) if aligned else 0
        passed = (builds_ok and pd and failures == 0 and missing_rate <= gates_cfg["max_variant_missing_rate"]
                  and min_n >= gates_cfg["min_effective_sample_size"] and coverage >= gates_cfg["min_locus_coverage"])
        gate = {"locus": rows[0]["locus"], "ancestry": ancestry, "ld_panel": panel["name"],
                "variant_missing_rate": f"{missing_rate:.6g}", "allele_alignment_failures": failures,
                "ld_positive_definite": pd, "effective_sample_size_min": f"{min_n:g}",
                "locus_coverage": f"{coverage:.6g}", "status": "PASS" if passed else "FAIL"}
        gates.append(gate)
        if not passed:
            continue
        # Instrument filtering is separate from fine-mapping and only this branch uses p-values.
        candidates = sorted((row for row in aligned if float(row["p_value"]) <= method["instrument_p_threshold"]), key=lambda x: float(x["p_value"]))
        kept = []
        for row in candidates:
            i = index[row["variant_id"]]
            if all(ld[i][index[prior["variant_id"]]] ** 2 < method["clump_r2"] for prior in kept):
                kept.append(row)
        for rank, row in enumerate(kept, 1):
            instruments.append({"locus": row["locus"], "signal_id": f"{row['locus']}:{ancestry}:I{rank}", "ancestry": ancestry,
                                "variant_id": row["variant_id"], "p_value": row["p_value"], "ld_panel": panel["name"]})
        if kept:
            leads.append(instruments[-len(kept)])
        # Fine-mapping uses every aligned locus row, including non-significant variants.
        prior_variance = method["prior_effect_variance"]
        weights = {}
        for row in aligned:
            variance = float(row["se"]) ** 2
            shrink = prior_variance / (prior_variance + variance)
            weights[row["variant_id"]] = math.sqrt(variance / (variance + prior_variance)) * math.exp(0.5 * row["z"] ** 2 * shrink)
        total = sum(weights.values()); pips = {key: value / total for key, value in weights.items()}
        ordered = sorted(pips, key=pips.get, reverse=True); cumulative = 0.0; credible = []
        for variant in ordered:
            credible.append(variant); cumulative += pips[variant]
            if cumulative >= method["credible_set_probability"]: break
        signal = f"{rows[0]['locus']}:{ancestry}:S1"
        for variant in ordered:
            results.append({"locus": rows[0]["locus"], "signal_id": signal, "ancestry": ancestry,
                            "credible_set_size": len(credible), "variant_id": variant,
                            "variant_pip": f"{pips[variant]:.12g}", "in_credible_set": variant in credible,
                            "integrated": False, "integration_method": "", "ld_panel": panel["name"]})
        per_ancestry[ancestry] = pips
    # Prespecified cross-ancestry product of PIPs over shared variants.
    if method["cross_ancestry_integration"] != "normalized_pip_product":
        raise SystemExit("unsupported cross-ancestry integration method")
    shared = set(per_ancestry.get("EUR", {})) & set(per_ancestry.get("EAS", {}))
    products = {v: per_ancestry["EUR"][v] * per_ancestry["EAS"][v] for v in shared}
    denominator = sum(products.values())
    for variant, value in sorted(products.items(), key=lambda x: x[1], reverse=True):
        results.append({"locus": results[0]["locus"], "signal_id": f"{results[0]['locus']}:CROSS:S1", "ancestry": "CROSS",
                        "credible_set_size": len(shared), "variant_id": variant, "variant_pip": f"{value / denominator:.12g}",
                        "in_credible_set": True, "integrated": True, "integration_method": "normalized_pip_product", "ld_panel": "EUR+EAS"})
    result_fields = ["locus", "signal_id", "ancestry", "credible_set_size", "variant_id", "variant_pip", "in_credible_set", "integrated", "integration_method", "ld_panel"]
    write_tsv(outdir / "fine_mapping_results.tsv", result_fields, results)
    write_tsv(outdir / "independent_instruments.tsv", ["locus", "signal_id", "ancestry", "variant_id", "p_value", "ld_panel"], instruments)
    write_tsv(outdir / "lead_only_instruments.tsv", ["locus", "signal_id", "ancestry", "variant_id", "p_value", "ld_panel"], leads)
    write_tsv(outdir / "ld_qc_gates.tsv", ["locus", "ancestry", "ld_panel", "variant_missing_rate", "allele_alignment_failures", "ld_positive_definite", "effective_sample_size_min", "locus_coverage", "status"], gates)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-dir", type=Path); args = parser.parse_args(); run(args.config, args.output_dir)
