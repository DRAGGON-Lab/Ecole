import csv
import gzip
import hashlib
import json
import random
from pathlib import Path

import pytest

from ecole.cli import build, verify
from ecole.experiment import prepare_coordinate_splits
from ecole.records import SequenceRecord


def test_guarded_coordinates_and_bounded_streams():
    prepared = prepare_coordinate_splits(
        [SequenceRecord("genome", "ACGT" * 100)],
        stream_length=64, minimum_stream_length=1, guard_bases=8,
    )
    intervals = prepared.coordinates
    assert intervals[0].stop + 8 == intervals[1].start
    assert intervals[1].stop + 8 == intervals[2].start
    assert all(len(record.sequence) <= 64 for record in prepared.train)


def test_end_to_end_source_contract(tmp_path: Path):
    source = tmp_path / "sources"
    output = tmp_path / "outputs"
    source.mkdir()
    rows = []
    shard_rows = []
    for role in ("train", "val", "test"):
        accession = {"train": "GCF_000000001.1", "val": "GCF_000000002.1",
                     "test": "GCF_000000003.1"}[role]
        rows.append(dict(split=role, accession=accession, completeness="99",
                         contamination="0", longest_contig="2000", quality_score="90",
                         selection_rank="1", clade_group=role))
        path = source / f"{role}_shard_00000.fa.gz"
        with gzip.open(path, "wt") as handle:
            handle.write(f">{accession} sample\n" + "".join(random.Random({"train": 1, "val": 2, "test": 3}[role]).choices("ACGT", k=2000)) + "\n")
        shard_rows.append(dict(split=role, shard_index="0",
                               sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    for path, fields, values in (
        (source / "accession_manifest.csv", rows[0].keys(), rows),
        (source / "fasta_shard_manifest.csv", shard_rows[0].keys(), shard_rows),
    ):
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(values)
    config = dict(dataset_id="smoke", source_dataset="fixture", source_folder_id="",
                  seed=17, genomes_per_role=1, stream_length=64, minimum_stream_length=1,
                  guard_bases=8, within_usable_bases=500, heldout_bases=100,
                  train_fraction=0.8, validation_fraction=0.1,
                  minimum_completeness=95, maximum_contamination=5,
                  scan_group_limit=24, manifest_ids={},
                  shard_ids={"train": ["fake"], "val": ["fake"], "test": ["fake"]})
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config))
    result = build(config_path, source, output)
    assert result["base_counts"] == {
        "train": 400, "validation_within": 50, "test_within": 50,
        "validation_clade": 100, "test_clade": 100,
    }
    assert verify(output)["dataset_fingerprint"] == result["dataset_fingerprint"]
    with pytest.raises(FileExistsError):
        build(config_path, source, output)
    with (output / "test_clade.fa.gz").open("ab") as handle:
        handle.write(b"tampered")
    with pytest.raises(ValueError, match="checksum"):
        verify(output)
