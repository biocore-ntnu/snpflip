# Changelog

## 0.1.0 — Rust implementation

- Replace the Python/pandas implementation with a standalone Rust CLI.
- Preserve every BIM row, its order, and the literal six input fields.
- Fix list headers, joins that multiply rows, and deduplication that deletes rows.
- Handle zero, negative, and out-of-range positions without indexing the wrong base.
- Distinguish orientation, palindromicity, mismatches, unsupported variants,
  invalid alleles, and missing data with detailed reasons.
- Retain every record in schema-versioned `.annotated.tsv` with 1-based positions.
- Exclude missing/repeated IDs from action lists and reference assignments.
- Support gzip/BGZF inputs, bounded streaming and indexed reference reads,
  arbitrary contigs, explicit aliases, and an optional human chromosome preset.
- Add `.ref_allele`, category lists, stderr summaries, nonzero errors, literal
  output paths, overwrite checks, temporary outputs, and rollback on publish errors.
- Keep only compact per-row state and re-read the BIM to write outputs, so memory
  is about 45 bytes per record (0.9 GB instead of 9 GB for 20M variants).
- Check output paths before reading any input.
- Accept `-` for stdin and pipes for either input; a non-regular BIM is spooled
  to a temporary file. Piped FASTA input no longer loses its first bytes.
- Report the sequence name and position of invalid FASTA bases.
- Apply `--chromosome-set human` codes to BIM names only.
- Add Rust regression tests and cross-platform CI/build workflows.

## Earlier Python releases

- 0.0.6 (2017-12-01): attempted duplicate-output fix.
- 0.0.5 (2016-01-03): chromosome overlap error message.
- 0.0.4 (2015-12-10): large BIM fix and chromosome-name help.
- 0.0.3 (2015-10-20): reference N handling.
- 0.0.2 (2015-10-14): extra FASTA chromosome handling and version flag.
