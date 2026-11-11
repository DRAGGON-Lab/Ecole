"""Download the fixed historical Google Drive source inventory when explicitly requested."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from ecole.cli import sha256


def download(config_path: Path, destination: Path, *, force: bool = False) -> None:
    """Download source manifests/shards; verify all shards against the original manifest."""
    try:
        import gdown
    except ImportError as error:
        raise RuntimeError("Install the optional downloader: pip install -e '.[download]'") from error
    import csv

    config = json.loads(config_path.read_text())
    destination.mkdir(parents=True, exist_ok=True)

    def fetch(file_id: str, path: Path) -> None:
        if force or not path.exists():
            temporary = path.with_name(path.name + ".part")
            temporary.unlink(missing_ok=True)
            downloaded = gdown.download(id=file_id, output=str(temporary), quiet=False)
            if not downloaded or not temporary.is_file():
                raise RuntimeError(f"download failed: {file_id}")
            temporary.replace(path)

    for key, name in (("accessions", "accession_manifest.csv"),
                      ("shards", "fasta_shard_manifest.csv")):
        fetch(config["manifest_ids"][key], destination / name)
    with (destination / "fasta_shard_manifest.csv").open(newline="") as handle:
        expected = {(row["split"], int(row["shard_index"])): row for row in csv.DictReader(handle)}
    for split, identifiers in config["shard_ids"].items():
        for index, file_id in enumerate(identifiers):
            path = destination / f"{split}_shard_{index:05d}.fa.gz"
            row = expected.get((split, index))
            if row is None:
                raise ValueError(f"missing manifest entry for {split}/{index}")
            if path.exists() and sha256(path) != row["sha256"]:
                path.unlink()
            fetch(file_id, path)
            if sha256(path) != row["sha256"]:
                raise ValueError(f"checksum mismatch for {split}/{index}: {path}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="ecole-download")
    parser.add_argument("--config", type=Path, default=Path("configs/ecoli_hybrid_v1.json"))
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)
    download(args.config, args.source, force=args.force)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
