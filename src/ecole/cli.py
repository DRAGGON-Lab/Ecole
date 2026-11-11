"""Build the established E. coli hybrid benchmark from verified local source shards.

No model or ATCG-FM dependency is needed. Downloads are separate from this builder.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

from ecole.experiment import prepare_coordinate_splits, write_fasta_lines
from ecole.fasta import iter_fasta
from ecole.records import SequenceRecord
from ecole.transforms import reverse_complement

SPLITS = ("train", "validation_within", "test_within", "validation_clade", "test_clade")
ACCESSION_RE = re.compile(r"(?:RS_|GB_)?(GC[AF]_\d+\.\d+)")
IUPAC = set("ACGTRYSWKMBDHVN-")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for part in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(part)
    return digest.hexdigest()


def accession_from_header(identifier: str, description: str) -> str | None:
    result = ACCESSION_RE.search(identifier + " " + description)
    return result.group(1) if result else None


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _eligible(rows: list[dict[str, str]], split: str, minimum: int, config: dict) -> list[dict[str, str]]:
    result = [
        row for row in rows
        if row["split"] == split
        and float(row["completeness"]) >= config["minimum_completeness"]
        and float(row["contamination"]) <= config["maximum_contamination"]
        and int(float(row["longest_contig"])) >= minimum
    ]
    return sorted(result, key=lambda row: (
        float(row["selection_rank"]) if row.get("selection_rank") else 1e12,
        -float(row["quality_score"]), row["accession"],
    ))


def _select(found: dict[str, SequenceRecord], rows: list[dict[str, str]], count: int) -> list[str]:
    chosen, groups = [], set()
    for row in rows:
        accession = row["accession"]
        group = str(row["clade_group"])
        if accession in found and group not in groups:
            chosen.append(accession)
            groups.add(group)
            if len(chosen) == count:
                break
    return chosen


def _collect(
    split: str, minimum: int, config: dict, source: Path,
    metadata: list[dict[str, str]], shards: list[dict[str, str]],
) -> tuple[list[str], dict[str, SequenceRecord], dict[str, dict[str, str]], list[dict]]:
    eligible = _eligible(metadata, split, minimum, config)
    candidates = {row["accession"] for row in eligible}
    group_by_accession = {row["accession"]: row["clade_group"] for row in eligible}
    found: dict[str, SequenceRecord] = {}
    used: list[dict] = []
    for index, file_id in enumerate(config["shard_ids"][split]):
        path = source / f"{split}_shard_{index:05d}.fa.gz"
        if not path.is_file():
            raise FileNotFoundError(f"missing verified input shard: {path}")
        expected = [row for row in shards if row["split"] == split and int(row["shard_index"]) == index]
        if len(expected) != 1:
            raise ValueError(f"missing or duplicate shard manifest entry for {split} {index}")
        digest = sha256(path)
        if digest != expected[0]["sha256"]:
            raise ValueError(f"source checksum mismatch: {path}")
        used.append({"split": split, "shard_index": index, "file_id": file_id, "sha256": digest})
        groups_seen: set[str] = set()
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            for record in iter_fasta(
                handle, source=str(path),
                record_filter=lambda identifier, description: accession_from_header(identifier, description) in candidates,
            ):
                accession = accession_from_header(record.identifier, record.description)
                if accession not in candidates or len(record.sequence) < minimum:
                    continue
                previous = found.get(accession)
                if previous is None or len(record.sequence) > len(previous.sequence):
                    found[accession] = record
                groups_seen.add(group_by_accession[accession])
                if len(groups_seen) >= config["scan_group_limit"]:
                    break
        chosen = _select(found, eligible, config["genomes_per_role"])
        if len(chosen) == config["genomes_per_role"]:
            return chosen, found, {r["accession"]: r for r in eligible}, used
    raise ValueError(f"{split}: insufficient eligible ANI groups ({len(_select(found, eligible, config['genomes_per_role']))})")


def _window(accession: str, record: SequenceRecord, length: int, role: str, seed: int) -> tuple[SequenceRecord, int]:
    available = len(record.sequence) - length
    if available < 0:
        raise ValueError(f"contig too short for {accession}")
    offset = int(hashlib.sha256(f"{seed}:{role}:{accession}".encode()).hexdigest()[:16], 16) % (available + 1)
    return SequenceRecord(
        identifier=f"{accession}__{record.identifier}__source_{offset}_{offset+length}",
        sequence=record.sequence[offset:offset+length].upper(),
        metadata={"accession": accession, "contig": record.identifier,
                  "source_start": str(offset), "source_stop": str(offset+length)},
    ), offset


def _streams(records: tuple[SequenceRecord, ...] | list[SequenceRecord],
             offset: int, role: str) -> list[SequenceRecord]:
    output = []
    for record in records:
        start = offset + int(record.metadata["start"])
        stop = offset + int(record.metadata["stop"])
        output.append(SequenceRecord(
            identifier=f"{record.metadata['accession']}__{record.metadata['contig']}__{role}__{start}_{stop}",
            sequence=record.sequence,
            metadata={**record.metadata, "role": role, "source_start": str(start), "source_stop": str(stop)},
        ))
    return output


def build(config_path: Path, source: Path, destination: Path) -> dict:
    """Construct independently verifiable, source-ordered FASTA partitions."""
    config = json.loads(config_path.read_text())
    metadata_path = source / "accession_manifest.csv"
    shards_path = source / "fasta_shard_manifest.csv"
    metadata, shards = _rows(metadata_path), _rows(shards_path)
    usable = config["within_usable_bases"]
    guard = config["guard_bases"]
    within_source = usable + 2 * guard
    heldout = config["heldout_bases"]
    chosen = {}
    found = {}
    selected_rows = {}
    files = []
    for role, minimum in (("train", within_source), ("val", heldout), ("test", heldout)):
        chosen[role], found[role], selected_rows[role], used = _collect(
            role, minimum, config, source, metadata, shards,
        )
        files.extend(used)

    role_groups = {role: {selected_rows[role][a]["clade_group"] for a in chosen[role]}
                   for role in ("train", "val", "test")}
    if role_groups["train"] & role_groups["val"] or role_groups["train"] & role_groups["test"] or role_groups["val"] & role_groups["test"]:
        raise ValueError("ANI clade leakage detected")

    splits: dict[str, list[SequenceRecord]] = defaultdict(list)
    coordinates = []
    for accession in chosen["train"]:
        record, offset = _window(accession, found["train"][accession], within_source, "within", config["seed"])
        prepared = prepare_coordinate_splits(
            [record], stream_length=config["stream_length"],
            minimum_stream_length=config["minimum_stream_length"],
            train_fraction=config["train_fraction"],
            validation_fraction=config["validation_fraction"],
            guard_bases=guard,
        )
        for internal, output in (("train", "train"), ("validation", "validation_within"), ("test", "test_within")):
            splits[output].extend(_streams(prepared.records(internal), offset, output))
        intervals = {item.split: [offset + item.start, offset + item.stop] for item in prepared.coordinates}
        coordinates.append({"accession": accession, "contig": found["train"][accession].identifier,
                            "window_start": offset, "window_stop": offset + within_source,
                            "train": intervals["train"], "validation_within": intervals["validation"],
                            "test_within": intervals["test"]})

    for role, output_name in (("val", "validation_clade"), ("test", "test_clade")):
        for accession in chosen[role]:
            original = found[role][accession]
            record, offset = _window(accession, original, heldout, output_name, config["seed"])
            for begin in range(0, len(record.sequence), config["stream_length"]):
                end = min(begin + config["stream_length"], len(record.sequence))
                if end - begin < config["minimum_stream_length"]:
                    continue
                start, stop = offset + begin, offset + end
                splits[output_name].append(SequenceRecord(
                    identifier=f"{accession}__{original.identifier}__{output_name}__{start}_{stop}",
                    sequence=record.sequence[begin:end],
                    metadata={"accession": accession, "contig": original.identifier, "role": output_name,
                              "source_start": str(start), "source_stop": str(stop)},
                ))

    canonical_hashes: dict[str, str] = {}
    for split in SPLITS:
        for record in splits[split]:
            invalid = set(record.sequence) - IUPAC
            if invalid:
                raise ValueError(f"invalid IUPAC sequence {record.identifier}: {sorted(invalid)}")
            key = hashlib.sha256(min(record.sequence, reverse_complement(record.sequence)).encode()).hexdigest()
            previous = canonical_hashes.setdefault(key, split)
            if previous != split:
                raise ValueError(f"cross-split reverse-complement duplicate: {record.identifier}")

    counts = {split: sum(len(record.sequence) for record in splits[split]) for split in SPLITS}
    size = config["genomes_per_role"]
    expected = {
        "train": int(usable * config["train_fraction"]) * size,
        "validation_within": int(usable * config["validation_fraction"]) * size,
        "test_within": (usable - int(usable * config["train_fraction"]) - int(usable * config["validation_fraction"])) * size,
        "validation_clade": heldout * size, "test_clade": heldout * size,
    }
    if counts != expected:
        raise ValueError(f"unexpected base counts: {counts} != {expected}")

    if destination.exists() and any(destination.iterdir()):
        raise FileExistsError(f"output directory must be empty: {destination}")
    destination.mkdir(parents=True, exist_ok=True)
    split_sha = {}
    fingerprint = hashlib.sha256()
    for split in SPLITS:
        path = destination / f"{split}.fa.gz"
        with path.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
                for line in write_fasta_lines(splits[split]):
                    zipped.write(line.encode())
        digest = sha256(path)
        split_sha[split] = digest
        fingerprint.update(split.encode())
        fingerprint.update(bytes.fromhex(digest))

    def role_metadata(role: str) -> list[dict]:
        return [
            {"accession": a, "clade_group": str(selected_rows[role][a]["clade_group"]),
             "quality_score": float(selected_rows[role][a]["quality_score"]),
             "completeness": float(selected_rows[role][a]["completeness"]),
             "contamination": float(selected_rows[role][a]["contamination"])}
            for a in chosen[role]
        ]

    result = {
        "schema_version": 2, "dataset_id": config["dataset_id"],
        "source_dataset": config["source_dataset"],
        "source_folder_id": config["source_folder_id"], "seed": config["seed"],
        "stream_length": config["stream_length"],
        "minimum_stream_length": config["minimum_stream_length"],
        "guard_bases": guard,
        "tokenizer_contract": "FixedAlphabetTokenizer IUPAC; BOS included, EOS excluded; one target per base; no reverse-complement augmentation",
        "split_policy": {
            "within": f"per-genome contiguous {config['train_fraction']:.0%}/{config['validation_fraction']:.0%}/{1-config['train_fraction']-config['validation_fraction']:.0%} with {guard}-base guards",
            "clade": "inherited ANI-99 clade-group holdout",
        },
        "roles": {"within": role_metadata("train"), "validation_clade": role_metadata("val"),
                  "test_clade": role_metadata("test")},
        "coordinate_contract": coordinates, "base_counts": counts,
        "stream_counts": {split: len(splits[split]) for split in SPLITS},
        "split_sha256": split_sha, "dataset_fingerprint": fingerprint.hexdigest(),
        "source_files": files, "source_manifest_sha256": {
            "accessions": sha256(metadata_path), "shards": sha256(shards_path),
        },
        "duplicate_policy": "canonical forward/reverse-complement SHA-256; cross-partition duplicates rejected",
        "producer": {"repository": "https://github.com/DRAGGON-Lab/Ecole", "configuration": config_path.name},
    }
    (destination / "manifest.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def verify(destination: Path) -> dict:
    manifest = json.loads((destination / "manifest.json").read_text())
    if manifest.get("schema_version") != 2:
        raise ValueError("unsupported manifest schema")
    digest = hashlib.sha256()
    for split in SPLITS:
        path = destination / f"{split}.fa.gz"
        file_hash = sha256(path)
        if file_hash != manifest["split_sha256"][split]:
            raise ValueError(f"checksum mismatch for {split}")
        digest.update(split.encode())
        digest.update(bytes.fromhex(file_hash))
        total = 0
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            for record in iter_fasta(handle):
                if set(record.sequence) - IUPAC:
                    raise ValueError(f"invalid DNA alphabet in {split}")
                total += len(record.sequence)
        if total != manifest["base_counts"][split]:
            raise ValueError(f"base count mismatch for {split}")
    if digest.hexdigest() != manifest["dataset_fingerprint"]:
        raise ValueError("dataset fingerprint mismatch")
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="ecole", description="Build and verify E. coli genome benchmarks")
    commands = parser.add_subparsers(dest="command", required=True)
    command = commands.add_parser("build", help="build the historical hybrid benchmark")
    command.add_argument("--config", type=Path, default=Path("configs/ecoli_hybrid_v1.json"))
    command.add_argument("--source", type=Path, required=True)
    command.add_argument("--output", type=Path, required=True)
    check = commands.add_parser("verify", help="verify generated FASTA and manifest")
    check.add_argument("--dataset", type=Path, required=True)
    arguments = parser.parse_args(argv)
    result = (build(arguments.config, arguments.source, arguments.output) if arguments.command == "build"
              else verify(arguments.dataset))
    print(json.dumps({"dataset_id": result["dataset_id"], "dataset_fingerprint": result["dataset_fingerprint"],
                      "base_counts": result["base_counts"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
