# ATCG-FM consumption and migration

ATCG-FM owns model tokenization and training. Ecole supplies five gzip FASTA files and `manifest.json` (`schema_version=2`).

## Consumption contract

```python
import gzip
from pathlib import Path
from atcg.sequence import FixedAlphabetTokenizer, OrderedCausalStreamDataset, parse_fasta

dataset_dir = Path("/path/to/ecole-release")
with gzip.open(dataset_dir / "train.fa.gz", "rt") as handle:
    records = parse_fasta(handle)
tokenizer = FixedAlphabetTokenizer()
dataset = OrderedCausalStreamDataset(
    records, tokenizer, segment_length=128, gradient_horizon=1,
    include_bos=True, include_eos=False,
)
```

Before any model run execute `ecole verify --dataset /path/to/ecole-release`; consumers without Ecole installed must independently verify the checksums in `manifest.json`. Record the manifest's `dataset_fingerprint` in the experiment configuration and checkpoint metadata.

## Migration phases

1. Move `notebooks/00_prepare_ecoli_dataset.ipynb` to Ecole as a **legacy reference**; production logic lives in the Python package.
2. Update ATCG-FM training notebooks to take an explicit `ECOLE_DATASET_DIR` and read only the five outputs and manifest. Do not automatically rebuild data during model training.
3. Mark old ATCG-FM dataset notebook as superseded only after validating a real Ecole release; remove it later in a separate PR.
4. Keep generic `atcg.sequence` tokenization and runtime dataset classes in ATCG-FM; only producer-specific workflows belong in Ecole.

## Reproducibility note

The new deterministic gzip writer has stable timestamps. Historical notebook gzip outputs may have different file hashes even when sequence contents are identical. A parity study should compare accession lists, coordinates, FASTA records and base counts rather than require byte-identical gzip.
