#!/usr/bin/env python3
"""Convert downloaded GIGASTROKE files into build-aware canonical outcomes.

This is deliberately an analysis stage: it never performs network access.  Inputs
must already have been produced by ``scripts/colab_download_gigastroke_gwas.py``.
Only Python's standard library is required so the adapter can be fixture-tested.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import sqlite3
from collections import defaultdict
from pathlib import Path

CANONICAL = [
    "dataset_id", "ancestry", "outcome", "role", "source_build", "target_build",
    "source_chromosome", "source_position", "source_ref", "source_alt",
    "source_variant_id", "chromosome", "position", "ref", "alt",
    "effect_allele", "other_allele", "beta", "se", "p_value", "eaf", "sample_size",
]


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(encoding="utf-8", newline="")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def chromosome(value: str) -> str:
    value = value.strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    value = value.upper()
    return "MT" if value in {"M", "MTDNA"} else value


class Fasta:
    def __init__(self, path: Path):
        self.path = path
        self.index_path = Path(str(path) + ".fai")
        if not self.index_path.exists():
            self._build_index()
        self.index = {}
        with self.index_path.open(encoding="ascii") as handle:
            for line in handle:
                name, length, offset, line_bases, line_width = line.rstrip().split("\t")[:5]
                self.index[chromosome(name)] = tuple(map(int, (length, offset, line_bases, line_width)))
        self.handle = self.path.open("rb")

    def _build_index(self) -> None:
        entries = []
        with self.path.open("rb") as handle:
            name = None; length = offset = line_bases = line_width = 0
            while True:
                line_start = handle.tell(); line = handle.readline()
                if not line:
                    if name is not None: entries.append((name, length, offset, line_bases, line_width))
                    break
                if line.startswith(b">"):
                    if name is not None: entries.append((name, length, offset, line_bases, line_width))
                    name = line[1:].split()[0].decode("ascii")
                    length = 0; offset = handle.tell(); line_bases = line_width = 0
                else:
                    bases = len(line.rstrip(b"\r\n"))
                    if bases and not line_bases: line_bases, line_width = bases, len(line)
                    length += bases
        temp = Path(str(self.index_path) + ".tmp")
        temp.write_text("".join("\t".join(map(str, item)) + "\n" for item in entries), encoding="ascii")
        temp.replace(self.index_path)

    def bases(self, chrom: str, pos: int, length: int) -> str:
        entry = self.index.get(chromosome(chrom))
        if not entry or pos < 1 or pos + length - 1 > entry[0]: return ""
        _, offset, line_bases, line_width = entry
        start = pos - 1; remaining = length; chunks = []
        while remaining:
            in_line = start % line_bases
            take = min(remaining, line_bases - in_line)
            self.handle.seek(offset + (start // line_bases) * line_width + in_line)
            chunks.append(self.handle.read(take))
            start += take; remaining -= take
        return b"".join(chunks).decode("ascii").upper()


class Chain:
    """Minimal UCSC chain reader returning every mapping for a 1-based base."""
    def __init__(self, path: Path):
        self.blocks: dict[str, list[tuple[int, int, str, int, int, str]]] = defaultdict(list)
        with path.open(encoding="ascii") as handle:
            current = None
            q_cursor = t_cursor = 0
            for raw in handle:
                line = raw.strip()
                if not line:
                    current = None
                    continue
                if line.startswith("chain "):
                    f = line.split()
                    _, _, t_name, _, t_strand, t_start, _, q_name, q_size, q_strand, q_start, _, _ = f
                    if t_strand != "+":
                        raise ValueError("target-negative UCSC chains are unsupported")
                    current = (chromosome(q_name), int(q_size), q_strand, chromosome(t_name))
                    q_cursor, t_cursor = int(q_start), int(t_start)
                    continue
                if current is None:
                    continue
                f = [int(x) for x in line.split()]
                size = f[0]
                q_name, q_size, q_strand, t_name = current
                self.blocks[q_name].append((q_cursor, q_cursor + size, t_name, t_cursor, q_size, q_strand))
                if len(f) == 3:
                    t_cursor += size + f[1]
                    q_cursor += size + f[2]

    def map(self, chrom: str, position: int) -> list[tuple[str, int, str]]:
        q0 = position - 1
        found = []
        for start, end, target, target_start, q_size, strand in self.blocks.get(chromosome(chrom), []):
            oriented = q0 if strand == "+" else q_size - q0 - 1
            if start <= oriented < end:
                found.append((target, target_start + oriented - start + 1, strand))
        return sorted(set(found))


def orient_allele(allele: str, strand: str) -> str:
    if strand == "+":
        return allele
    return allele.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def normalize(fasta: Fasta, chrom: str, pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    ref, alt = ref.upper(), alt.upper()
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    # VCF-style left alignment for repeat indels.
    while len(ref) != len(alt) and pos > 1:
        previous = fasta.bases(chrom, pos - 1, 1)
        if not previous or ref[-1] != previous or alt[-1] != previous:
            break
        ref, alt, pos = previous + ref[:-1], previous + alt[:-1], pos - 1
    return pos, ref, alt


def pick(row: dict[str, str], names: list[str], required: bool = True) -> str:
    for name in names:
        if name in row and row[name] not in (None, ""):
            return row[name]
    if required:
        raise ValueError("missing column/value: " + "/".join(names))
    return ""


def source_variant(row: dict[str, str], columns: dict[str, list[str]], source_fasta: Fasta) -> tuple[str, str, str]:
    """Return source variant ID/ref/alt, deriving alleles from a coordinate ID.

    GIGASTROKE/GWAS Catalog files do not consistently expose separate REF and
    ALT columns. Their variant_id/hm_variant_id commonly encodes
    chromosome_position_ref_alt (colon and underscore separators are accepted).
    We only accept a parsed ID when its chromosome and position agree with the
    row, avoiding silent use of a harmonized ID from a different build.
    """
    variant_id = pick(row, columns.get("variant_id", []), required=False)
    ref = pick(row, columns.get("ref", []), required=False).upper()
    alt = pick(row, columns.get("alt", []), required=False).upper()
    if ref and alt:
        return variant_id, ref, alt

    chrom = chromosome(pick(row, columns["chromosome"]))
    position = int(pick(row, columns["position"]))
    candidates = columns.get("ref_alt_variant_id", []) + columns["variant_id"]
    for name in dict.fromkeys(candidates):
        encoded = row.get(name, "")
        fields = re.split(r"[:_]", encoded)
        if len(fields) < 4:
            continue
        encoded_chrom, encoded_pos, encoded_ref, encoded_alt = fields[-4:]
        if (chromosome(encoded_chrom) == chrom and encoded_pos.isdigit()
                and int(encoded_pos) == position
                and re.fullmatch(r"[ACGTN]+", encoded_ref.upper())
                and re.fullmatch(r"[ACGTN]+", encoded_alt.upper())):
            return variant_id, encoded_ref.upper(), encoded_alt.upper()
    # The published GIGASTROKE tables can contain only effect/other alleles. In
    # that case the GRCh37 FASTA, not effect direction, determines REF and ALT.
    effect = pick(row, columns["effect_allele"]).upper()
    other = pick(row, columns["other_allele"]).upper()
    if not all(re.fullmatch(r"[ACGTN]+", allele) for allele in (effect, other)) or effect == other:
        raise ValueError("source_reference_allele_mismatch")
    effect_matches = source_fasta.bases(chrom, position, len(effect)) == effect
    other_matches = source_fasta.bases(chrom, position, len(other)) == other
    if effect_matches == other_matches:
        raise ValueError("source_reference_allele_mismatch")
    return (variant_id, effect, other) if effect_matches else (variant_id, other, effect)


def run(config_path: Path, output_override: Path | None = None, dataset_id: str | None = None) -> None:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    root = config_path.parent
    resolve = lambda p: (root / p).resolve() if not Path(p).is_absolute() else Path(p)
    chain_path = resolve(config["liftover"]["chain"])
    source_reference_path = resolve(config["liftover"]["source_reference"])
    reference_path = resolve(config["liftover"]["target_reference"])
    for path, key in ((chain_path, "chain_sha256"), (source_reference_path, "source_reference_sha256"),
                      (reference_path, "target_reference_sha256")):
        expected = config["liftover"][key].lower()
        if sha256(path) != expected:
            raise SystemExit(f"checksum verification failed for {path}")
    chain, source_fasta, fasta = Chain(chain_path), Fasta(source_reference_path), Fasta(reference_path)
    outdir = output_override.resolve() if output_override else resolve(config["output_dir"])
    outdir.mkdir(parents=True, exist_ok=True)
    selected = {config["selection"]["eur_discovery"], *config["selection"]["eas_replication_subtypes"]}
    if dataset_id and dataset_id not in selected:
        raise SystemExit(f"--dataset-id is not selected in config: {dataset_id}")
    indexed = {item["id"]: item for item in config["datasets"]}
    missing = selected - indexed.keys()
    if missing:
        raise SystemExit("selected dataset IDs are not configured: " + ", ".join(sorted(missing)))
    discovery = indexed[config["selection"]["eur_discovery"]]
    if discovery["ancestry"] != "EUR" or discovery["role"] != "discovery":
        raise SystemExit("selection.eur_discovery must identify an EUR discovery dataset")
    if any(indexed[item]["ancestry"] != "EAS" or indexed[item]["role"] not in {"replication", "replication_subtype"}
           for item in config["selection"]["eas_replication_subtypes"]):
        raise SystemExit("EAS selections must identify replication or replication_subtype datasets")
    manifests = []
    for dataset in config["datasets"]:
        if dataset["id"] not in selected or (dataset_id and dataset["id"] != dataset_id):
            continue
        source = resolve(dataset["input"])
        output = outdir / f"{dataset['id']}.canonical.tsv"; reject = outdir / f"{dataset['id']}.rejected.tsv"
        database = outdir / f".{dataset['id']}.duplicates.sqlite3"
        if database.exists(): database.unlink()
        accepted_count = rejected_count = 0; counts = defaultdict(int)
        connection = sqlite3.connect(database)
        connection.execute("CREATE TABLE seen (variant TEXT PRIMARY KEY)")
        with open_text(source) as handle, output.open("w", newline="", encoding="utf-8") as accepted_handle, reject.open("w", newline="", encoding="utf-8") as rejected_handle:
            accepted_writer = csv.DictWriter(accepted_handle, delimiter="\t", fieldnames=CANONICAL, lineterminator="\n"); accepted_writer.writeheader()
            reject_fields = ["dataset_id", "source_line", "reason", "source_chromosome", "source_position", "source_ref", "source_alt", "source_variant_id"]
            rejected_writer = csv.DictWriter(rejected_handle, delimiter="\t", fieldnames=reject_fields, lineterminator="\n"); rejected_writer.writeheader()
            reader = csv.DictReader(handle, delimiter=dataset.get("delimiter", "\t"))
            for line_no, row in enumerate(reader, 2):
                try:
                    c = chromosome(pick(row, dataset["columns"]["chromosome"])); p = int(pick(row, dataset["columns"]["position"]))
                    vid, ref, alt = source_variant(row, dataset["columns"], source_fasta)
                    mappings = chain.map(c, p)
                    if len(mappings) != 1:
                        raise ValueError("unmapped" if not mappings else "multi_mapped")
                    tc, tp, strand = mappings[0]
                    tref, talt = orient_allele(ref, strand), orient_allele(alt, strand)
                    tp, tref, talt = normalize(fasta, tc, tp, tref, talt)
                    if fasta.bases(tc, tp, len(tref)) != tref:
                        raise ValueError("reference_allele_mismatch")
                    key = (tc, tp, tref, talt)
                    try:
                        connection.execute("INSERT INTO seen VALUES (?)", ("|".join(map(str, key)),))
                    except sqlite3.IntegrityError:
                        raise ValueError("duplicate")
                    col = dataset["columns"]
                    accepted_writer.writerow(dict(zip(CANONICAL, [dataset["id"], dataset["ancestry"], dataset["outcome"], dataset["role"], dataset["source_build"], "GRCh38", c, p, ref, alt, vid, tc, tp, tref, talt, pick(row, col["effect_allele"]), pick(row, col["other_allele"]), pick(row, col["beta"]), pick(row, col["se"]), pick(row, col["p_value"]), pick(row, col.get("eaf", []), False), pick(row, col.get("sample_size", []), False)])))
                    accepted_count += 1
                except (ValueError, KeyError) as error:
                    reason = str(error)
                    if reason not in {"unmapped", "multi_mapped", "duplicate", "reference_allele_mismatch", "source_reference_allele_mismatch"}:
                        reason = "invalid_input:" + reason
                    rejected_writer.writerow({"dataset_id": dataset["id"], "source_line": line_no, "reason": reason,
                                     "source_chromosome": row.get(dataset["columns"]["chromosome"][0], ""), "source_position": row.get(dataset["columns"]["position"][0], ""),
                                     "source_ref": row.get(dataset["columns"]["ref"][0], ""), "source_alt": row.get(dataset["columns"]["alt"][0], ""), "source_variant_id": row.get(dataset["columns"]["variant_id"][0], "")})
                    rejected_count += 1; counts[reason] += 1
        connection.close(); database.unlink(missing_ok=True)
        manifests.append({"dataset_id": dataset["id"], "ancestry": dataset["ancestry"], "outcome": dataset["outcome"], "role": dataset["role"], "source_build": dataset["source_build"], "target_build": "GRCh38", "input": str(source), "output": str(output), "accepted": accepted_count, "rejected": rejected_count, "rejection_counts": dict(counts), "chain": str(chain_path), "chain_sha256": sha256(chain_path), "source_reference": str(source_reference_path), "source_reference_sha256": sha256(source_reference_path), "target_reference": str(reference_path), "target_reference_sha256": sha256(reference_path)})
    manifest_path = outdir / "dataset_manifest.json"
    previous = json.loads(manifest_path.read_text()) if dataset_id and manifest_path.exists() else []
    merged = {item["dataset_id"]: item for item in previous}
    merged.update({item["dataset_id"]: item for item in manifests})
    manifest_path.write_text(json.dumps(list(merged.values()), indent=2) + "\n", encoding="utf-8")


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-dir", type=Path, help="Override config output_dir (useful for immutable fixture runs).")
    parser.add_argument("--dataset-id", help="Process only one selected dataset (for resumable batch execution).")
    args = parser.parse_args(); run(args.config, args.output_dir, args.dataset_id)
