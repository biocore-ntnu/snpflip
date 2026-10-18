# snpflip

Check PLINK BIM allele orientation against a reference genome. Writes headerless
variant ID lists, reference allele assignments, and a report containing **every
input row in its original order**.

This Rust successor to Python snpflip 0.0.6 needs no Python, pandas, htslib, or
write access to the reference directory.

## Install

Prebuilt binaries for each release are on the
[releases page](https://github.com/biocore-ntnu/snpflip/releases), with
`SHA256SUMS`. To build with Rust 1.91 or newer:

```sh
cargo install snpflip --locked          # from crates.io
cargo install --path . --locked         # from a checkout
snpflip --help
```

Or run `cargo build --release --locked` for `target/release/snpflip`
(`snpflip.exe` on Windows). The release-build workflow prepares archives for Linux
x86-64 (static musl), macOS Apple Silicon and Intel, and Windows x86-64.

## Usage

```sh
snpflip -b cohort.bim -f reference.fa -o qc/cohort

# Human numeric sex/mitochondrial chromosome codes:
snpflip -b cohort.bim.gz -f reference.fa.gz -o qc/cohort --chromosome-set human

# Included example:
cargo run -- -b examples/example.bim -f examples/example.fa -o qc/example \
  --chromosome-set human
```

The original `-b`/`--bim-file`, `-f`/`--fasta-genome`, `-o`/`--output-prefix`,
`--help`, and `-v`/`--version` options retain their meanings. The default prefix
is `snpflip_output/<BIM basename>` with compression and BIM extensions removed.
Use `--force` to replace existing outputs.

BIM and FASTA inputs can be plain, gzip, concatenated gzip, or BGZF; compression
is detected from contents. BIM requires exactly six whitespace-separated fields
per nonempty line. IDs and field spellings are preserved, including leading zeros
and literal `NA` or `null`. Positions must be signed 64-bit integers and genetic
distances finite numbers. Blank BIM lines are skipped. Malformed rows fail with
line context.

Either input may be `-` for standard input, and pipes or process substitution
such as `-b <(zcat cohort.bim.gz)` work. snpflip reads the BIM more than once to
keep memory low, so a BIM that is not a regular file is copied, decompressed,
to an unnamed temporary file (in `TMPDIR`) while it is read. A regular BIM file
that changes during a run fails. A piped FASTA is always streamed.

## Interpretation

Forward/reverse mean **compatible with the supplied reference**. They do not
prove the genome build or assay strand convention. Investigate high mismatch
rates. This tool does not perform liftover or infer strand conventions from
frequencies, LD, or probe sequences.

Only two single-base A/C/G/T alleles and an A/C/G/T reference base can receive a
forward/reverse classification. Case is ignored for comparison; equal alleles
such as A/A are allowed. Missing/ambiguous alleles (`0`, `.`, IUPAC ambiguity
codes) stay unresolved. Multi-base sequence alleles and symbolic `I`, `D`, `-`,
or `*` are unsupported non-SNPs. Other allele codes are invalid. Indels are never
classified by comparing them to a single reference base.

| Category | Meaning |
| --- | --- |
| `forward` | An allele matches reference; neither matches its complement |
| `reverse` | An allele matches the complement; neither matches reference |
| `palindromic` | A/T or C/G compatible with both orientations |
| `mismatch` | Valid SNP alleles match neither orientation |
| `non_snp` | Unsupported sequence or symbolic non-SNP allele |
| `missing` | Missing/ambiguous allele or unavailable/unresolved reference |
| `invalid` | Invalid allele code |

Palindromic SNPs are never automatically flipped. Palindromicity is a separate
flag: A/T against C is both a mismatch and palindromic; A/T against N is missing
and palindromic. Allele problems take precedence in `category`, while `reason`
retains allele and reference problems. Invalid positions can also have a missing
chromosome.

## Outputs

Every run creates these files, including empty lists:

| Suffix | Contents |
| --- | --- |
| `.annotated.tsv` | Every original row plus annotations |
| `.forward`, `.reverse` | IDs compatible with each orientation |
| `.palindromic`, `.mismatch`, `.non_snp`, `.missing`, `.invalid` | IDs in each unresolved category |
| `.ambiguous` | Union of all unresolved categories, in input order |
| `.ref_allele` | Two columns: ID and forward reference base, for resolved rows |

**Every ID list and `.ref_allele` excludes all repeated IDs and missing IDs (`.`).**
This applies even when repeated records agree: ID lists cannot uniquely target
those rows. They remain in the report with `id_status=duplicate` or `missing`;
the stderr summary counts withheld rows. Resolve those IDs before relying on
these lists to process every variant. Lists have no headers and use LF newlines.

Annotation schema version **1** begins with the six original BIM fields, with
**1-based positions**, then:

| Column | Meaning |
| --- | --- |
| `schema_version` | `1` |
| `row` | 1-based record ordinal, excluding blank lines |
| `input_line` | Physical BIM line number |
| `reference_chromosome` | Actual FASTA name, or empty if absent |
| `reference`, `reference_rev` | Uppercase reference and IUPAC complement, or empty |
| `strand` | `forward`, `reverse`, or legacy-compatible `ambiguous` |
| `category` | Detailed category above |
| `reason` | Semicolon-separated reasons; empty for resolved rows |
| `palindromic` | `true` or `false`, independent of category |
| `id_status` | `unique`, `duplicate`, or `missing` |

Reasons are `invalid_allele`, `unsupported_non_snp`,
`missing_or_ambiguous_allele`, `missing_chromosome`, `zero_position`,
`negative_position`, `out_of_range`, `reference_n`, `ambiguous_reference`,
`palindromic`, and `reference_mismatch`.

The annotation is **not a replacement BIM**: it has a header and extra columns.
Never substitute it for the BIM accompanying a BED file.

## Chromosome names

The default `--chromosome-set generic` removes one optional lowercase `chr`
prefix and treats `M`/`MT` as aliases. Arbitrary complete contig names are
supported: `2L`, `chrEBV`, `chrI`, `NC_000001.11`, etc. `2L` never becomes `2`.
Other name comparisons are case-sensitive.

`--chromosome-set human` additionally maps `23 -> X`, `24 -> Y`, `25 -> X`,
`26 -> MT`, and `XY`/`PAR1`/`PAR2` to X. PAR coordinates must be expressed on X
in the selected build; Y coordinates are not converted. These codes apply to BIM
names only, so FASTA contigs such as `PAR1` keep their own names. Numeric names such as
`23` remain ordinary contigs in generic mode. PLINK sex-chromosome codes depend
on the species' chromosome set; see [PLINK chromosome sets](https://www.cog-genomics.org/plink/1.9/input#chr_set).

`--chromosome-map aliases.tsv` takes two columns with no header:

```text
# BIM chromosome    FASTA name
1                   NC_000001.11
2                   NC_000002.12
```

Mappings apply once to BIM names, never to FASTA names, and replace the
preset's codes for the mapped names (only the `chr`/`M` normalization applies to
the target); they are not chained. Duplicate source entries fail. Colliding FASTA names
such as both `1` and `chr1`, or `M` and `MT`, fail instead of silently choosing one.

## Reference access

`--reference-mode auto` uses an existing `<reference>.fai` for plain FASTA, or
`.fai` plus `<reference>.gzi` for BGZF. Otherwise it scans the reference once with
a bounded buffer, retaining only requested bases. Ordinary gzip streams when
only a `.fai` is present. No indexes are created or modified.

`--reference-mode stream` always scans. `--reference-mode indexed` requires
indexes and a regular FASTA file, and fails otherwise. `--fasta-index path/to/index.fai` selects an
alternative existing FASTA index and cannot be combined with stream mode. The
BGZF `.gzi` must still sit beside its reference.

Indexes must match the exact reference contents. Short/invalid queries fail, but
not every stale index is detectable. Use stream mode if index provenance is
uncertain. Indexed reads combine nearby sites into windows of at most 64 KiB.
Memory scales with BIM records, not chromosome length: about 45 bytes per record,
or roughly 0.9 GB for 20 million variants.

Reference sequences must contain only IUPAC DNA codes (any case). Anything else,
including whitespace inside a sequence, fails with the sequence name and
position, since it usually means the wrong file. Stream mode checks every base;
indexed mode only reads the windows around requested sites.

## PLINK workflow

Inspect the summary and unresolved rows first. Flip resolved reverse SNPs, remove
uniquely identified unresolved variants, then assign reference alleles:

```sh
plink2 --bfile cohort --flip qc/cohort.reverse --exclude qc/cohort.ambiguous \
  --make-pgen --out cohort.stranded
plink2 --pfile cohort.stranded --ref-allele qc/cohort.ref_allele 2 1 \
  --make-pgen --out cohort.forward
```

For PLINK 1.9, use `--make-bed`/`--bfile` and
`--a2-allele qc/cohort.ref_allele 2 1` in step two. `--a1-allele` does not mean
reference assignment. See [PLINK 1.9 allele ordering](https://www.cog-genomics.org/plink/1.9/data#keep_allele_order)
and [PLINK 2 reference assignment](https://www.cog-genomics.org/plink/2.0/data#ref_allele).
Missing/repeated IDs receive no action from these lists.

Related tools: [PLINK 2 `--ref-from-fa`](https://www.cog-genomics.org/plink/2.0/data#ref_allele)
for reference assignment, and [bcftools `+fixref`](https://samtools.github.io/bcftools/howtos/plugin.fixref.html)
for VCF diagnostics/correction. Understand the data's build and strand convention
before transforming it.

## Errors and output replacement

Successful analysis exits 0, including analysis with unresolved rows. Category
and reason counts, reference access mode, and withheld-ID counts go to stderr.
Stdout is empty except for help/version. Runtime failures exit 1; CLI usage
errors exit 2. No chromosome overlap is a runtime failure.

Output paths are checked before any input is read, and again before writing.
Outputs are written during a final pass over the BIM into temporary files in the
destination directory, then flushed. Existing files require
`--force`; symlinks, directories, and input files are never replaced. Publication
errors roll back new files and restore previous outputs. Individual publications
are atomic, but the whole set is not a crash-atomic transaction. Do not read the
files while a run is writing.

## Development and migration

```sh
cargo test --locked
cargo fmt --check
cargo clippy --locked --all-targets -- -D warnings
```

Tests cover all 64 canonical reference/allele combinations, allele ordering and
complementation, regression bugs, literal IDs, ambiguous data, chromosome aliases,
gzip/BGZF, stream/index parity, output rollback, and CLI errors. CI tests Linux,
macOS, Windows, and the minimum Rust version. The release-build workflow runs
manually or on a `v*` tag; tags matching the crate version publish a GitHub
release with binaries and `SHA256SUMS`.

The Python package, executable, and tests have been replaced. Existing PyPI
installations are separate and must be replaced on PATH. `.annotated_bim` becomes
schema-versioned `.annotated.tsv` with 1-based coordinates instead of the old
internal 0-based coordinates. `.reverse`/`.ambiguous` retain their names but omit
headers and unsafe IDs. See [CHANGELOG.md](CHANGELOG.md).

License: GPL-3.0, preserving the original project's license.

![The Strand Mixture: an old Norwegian tobacco tin](IMG_0401.jpg)
