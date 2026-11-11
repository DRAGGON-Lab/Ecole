# Product: Ecole

## Problem

Genomic foundation-model experiments should not implement or silently revise their own dataset selection logic. ATCG-FM currently carries an E. coli benchmark creation notebook and generic preparation helpers. That makes data releases difficult to version independently of training code.

## Users

DRAGGON Lab bioinformaticians, dataset curators, reproducibility reviewers, and genomic foundation-model developers.

## V1 goals

- Reproduce the documented E. coli-related hybrid benchmark from verified inputs.
- Preserve five split identities and the established manifest schema.
- Capture accession versions, metadata-defined ANI groups, source coordinates and SHA-256 hashes.
- Generate locally and verify offline without GPU dependencies.
- Hand off versioned dataset directories to ATCG-FM.

## Non-goals

Ecole does not train or evaluate neural networks, implement model tokenizers, manage checkpoints, or perform an exhaustive E. coli taxonomy/ANI recuration in the first release.

## Acceptance criteria

The pipeline builds all five expected split outputs using the source inventory, computes the expected quotas, rejects cross-role group overlap and cross-partition exact reverse-complement duplicates, and verifies its artifacts. Tests run on tiny synthetic data without a network connection.

## Future work

Accession-native NCBI/GTDB download, explicit release versioning, near-duplicate homology screening, independent taxonomic verification, portable Parquet / token-shard export and publishing to `gs://draggon-lab-data`. Each requires its own reviewable PR.
