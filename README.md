# Ecole

**Ecole** is DRAGGON Lab's independently versioned, reproducible genome-dataset production repository. Its initial supported workflow is the **E. coli-related hybrid benchmark** formerly prepared inside [ATCG-FM](https://github.com/DRAGGON-Lab/ATCG-FM).

**Boundaries:** Ecole owns acquisition manifests, genome selection, ANI-aware split assignments, sequence extraction, QC and data-release manifests. ATCG-FM owns tokenization, model code, data loading, training, evaluation and checkpoints. Ecole has no PyTorch dependency.

## Quick start

Requires Python 3.12+; no Docker or GPU required.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[download,dev]'
python -m ecole.download --config configs/ecoli_hybrid_v1.json --source data/raw/ecoli-hybrid-v1
ecole build --config configs/ecoli_hybrid_v1.json --source data/raw/ecoli-hybrid-v1 --output data/outputs/stage-small
ecole verify --dataset data/outputs/stage-small
pytest
```

The source download uses the **fixed historical Google Drive IDs**, and the builder verifies source-shard SHA-256 hashes against the upstream shard manifest. Download credentials, access quotas and permissions depend on the source. Data is not bundled with this repository.

### Outputs

```text
data/outputs/stage-small/
  train.fa.gz                 # 16M bases
  validation_within.fa.gz    # 2M bases
  test_within.fa.gz          # 2M bases
  validation_clade.fa.gz     # 2M bases
  test_clade.fa.gz           # 2M bases
  manifest.json              # schema_version 2; accession groups, coordinates, checksums
```

ATCG-FM can consume these files using its existing `parse_fasta`, `FixedAlphabetTokenizer` and `OrderedCausalStreamDataset` APIs. See [ARCHITECTURE.md](ARCHITECTURE.md) for the contract and [docs/ATCG-FM.md](docs/ATCG-FM.md) for the migration.

### Quality and scientific scope

The default benchmark uses a deterministic seed (17), 16 ANI-99-distinct source genomes per role, minimum estimated completeness 95%, maximum contamination 5%, and 128-base guards at within-genome split boundaries. Its fixed source lists and derived hashes are captured in `configs/ecoli_hybrid_v1.json` and each output manifest. The split groups are inherited from the source manifests: Ecole does **not** recompute ANI-99 from genomes.

This is an *E. coli-related* benchmark selected from the existing bacterial corpus, **not** a verified all-*E. coli* species-wide survey. Taxonomic coverage, cross-split sequence homology, and ANI-group quality beyond the supplied metadata need independent validation before stronger biological claims.

## Repository

- `configs/`: versioned dataset recipes and source file IDs
- `src/ecole/`: FASTA processing, deterministic split logic, builder, verifier, downloader
- `tests/`: synthetic, offline contract tests
- `docs/`: producer/consumer interface and provenance policy
- `ARCHITECTURE.md`, `PRODUCT.md`, `AGENT.md`, `ADR-001.md`: maintenance and engineering decisions

No large raw data, credentials, model checkpoints or experimental artifacts are committed.

## Maintainer workflow

ChatGPT coordinates architecture, decisions and scoped work. Codex implements small PRs against the contracts in `ARCHITECTURE.md` and `AGENT.md`. New benchmark versions must never silently modify an existing published dataset fingerprint.

## Storage

Keep dataset bytes in an object store (the current lab GCS location is `gs://draggon-lab-data`), with read-only training identities and checksum verification. Make publication a separate reviewed operation; do not run destructive storage synchronization from the builder.

## Container policy

The production tool is a Python CLI. Containers are optional for portability and not needed for local CPU development. Dataset ingestion and training operate independently.

Licensed under MIT. See [LICENSE](LICENSE).
