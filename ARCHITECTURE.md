# Architecture

## System boundary

`Ecole (producer) -> versioned dataset directory + manifest -> ATCG-FM (consumer)`

Ecole is a model-independent data compiler. It must not import `atcg`, PyTorch or any ATCG model runtime. ATCG-FM must not require Ecole at training time; it reads released bytes and their metadata.

## Components

1. `download.py`: explicitly downloads the pinned historical source files, checks source-shard SHA-256, and does not mutate model code.
2. `cli.py`: selects genome accessions, chooses source windows, prepares guarded coordinate splits, checks duplicate reverse complements, writes deterministic compressed FASTA, and produces/validates manifest schema 2.
3. `fasta.py`, `records.py`, `experiment.py`, `transforms.py`: extracted deterministic sequence primitives with original behavior, maintained independently.
4. `configs/`: immutable input recipes, versions, source identifiers and selection thresholds.
5. `tests/`: offline synthetic fixtures for split and artifact contracts.

## Scientific data model

- **Genome identity**: exact assembly accession including version.
- **ANI group**: metadata-provided `clade_group`; this is a grouping assertion, not an ANI calculation performed by Ecole.
- **Sequence**: source contig sequence, original header, local and global half-open coordinates.
- **Role**: train / within-genome held-out / clade-held-out.
- **Release**: named dataset, immutable content checksums, source manifest checksums, selection parameters, and per-split coordinates.

Within-genome test windows are unseen coordinates of a training accession; they do **not** imply source-genome independence. Clade test windows are from distinct metadata-provided ANI-99 groups. A repeat/homology leakage analysis across partitions remains necessary for rigorous benchmark claims.

## Compatibility contract

Output names remain exactly `train.fa.gz`, `validation_within.fa.gz`, `test_within.fa.gz`, `validation_clade.fa.gz`, and `test_clade.fa.gz`. `manifest.json` retains `schema_version=2`, `split_sha256`, `dataset_fingerprint`, `roles`, `coordinate_contract`, `base_counts`, and `stream_counts`. Additive `producer` metadata is allowed. Do not reinterpret `dataset_fingerprint`: it is the SHA-256 over each ordered split name and that split's compressed-file SHA-256 bytes.

ATCG-FM's baseline experiment uses the IUPAC fixed alphabet with BOS and no EOS, segment length 128, and ordered streams. Encoding and training policies belong to ATCG-FM rather than the producer.

## Operational policy

- Raw inputs stay in `data/raw/` (ignored by Git).
- Never overwrite a populated output directory; releases are immutable.
- Validate file hashes and sequence counts before consumption.
- Pin source and package versions; log manifest fingerprint for every training run.
- Prefer separate service accounts: Ecole publishing credentials can write; ATCG training credentials read only.
- The reference implementation is CPU-first and needs no network access during build or validation.
- Download is opt-in and uses external Google Drive access; no credentials live in the repository.
- A scheduler/container is optional; no persistent service, auth endpoint or frontend is required.

## Risks / follow-ups

- Source-shard selection stops scanning each shard after 24 distinct qualifying groups, matching the historical notebook's bounded scan heuristic; ensure this is scientifically appropriate when expanding the dataset.
- Input groups and FASTA headers must be reconciled to exact assemblies, versions, and taxonomy.
- Fixed ID download source may be unavailable; future releases should prefer accession-based NCBI/GTDB retrieval with versioned checksums.
- Legacy gzip compression was timestamp-dependent; the new producer uses zero-mtime gzip for repeatable *new* releases, so it does not promise historical byte-identical hashes.
- The script enforces exact canonical duplicate rejection, not an exhaustive approximate homology search.

See `ADR-001.md` for separation rationale.
