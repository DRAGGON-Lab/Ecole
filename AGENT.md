# Agent coordination

## Mission

Keep Ecole the reproducible, model-neutral source of truth for microbial dataset creation, while ATCG-FM is the training consumer.

## ChatGPT / Codex

ChatGPT owns architecture, scientific assumptions, scope and reviews. Codex implements small, reviewable changes, runs tests and documents the exact command/results. Stop to consult ChatGPT before altering the split contract, source identities, fingerprint definition or publication procedure.

## Working rules

- Preserve the current CLI and manifest schema unless changing them in an explicit versioned ADR.
- Never alter or delete published artifacts in place.
- Make deterministic transforms pure and test their edge cases.
- Do not add model, GPU, ATCG runtime or cloud credentials dependencies.
- Keep source checksums, full accession versions, ANI grouping and sequence coordinates visible.
- Reject leakage and discrepancies; never silently discard missing assemblies or invalid inputs.
- Do not upload multi-gigabyte FASTA or pretrained checkpoints into Git.
- Update README, ARCHITECTURE, tests and manifest documentation with code changes.
- Run `python -m pytest` and `python -m compileall -q src`; report limitations if full-data reconstruction is unavailable.
- Keep each PR scoped to one concern; do not refactor unrelated ATCG-FM model code.
- Do not claim biological or dataset reproducibility until actual source files and generated manifests have been validated.

## Handoff

Every implementation PR states the dataset contract affected, exact inputs used, tests run, generated manifest hashes, any release incompatibility, and remaining scientific uncertainty. A change to genome selection, source version or split protocol requires an ADR and a new dataset version.
