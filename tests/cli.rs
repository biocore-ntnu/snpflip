use std::{
    collections::BTreeMap,
    fs,
    io::Write,
    path::{Path, PathBuf},
    process::{Command, Output},
};

use flate2::{Compression, write::GzEncoder};
use snpflip::reference::sidecar;
use tempfile::{TempDir, tempdir};

struct Fixture {
    _root: TempDir,
    fasta: PathBuf,
    bim: PathBuf,
    prefix: PathBuf,
}

impl Fixture {
    fn new(fasta: &str, bim: &str) -> Self {
        let root = tempdir().unwrap();
        let fa = root.path().join("reference.fa");
        let bi = root.path().join("input.bim");
        let prefix = root.path().join("output");
        fs::write(&fa, fasta).unwrap();
        fs::write(&bi, bim).unwrap();
        Self {
            _root: root,
            fasta: fa,
            bim: bi,
            prefix,
        }
    }
    fn run(&self, extra: &[&str]) -> Output {
        Command::new(env!("CARGO_BIN_EXE_snpflip"))
            .arg("-f")
            .arg(&self.fasta)
            .arg("-b")
            .arg(&self.bim)
            .arg("-o")
            .arg(&self.prefix)
            .args(extra)
            .output()
            .unwrap()
    }
    fn success(&self, extra: &[&str]) -> String {
        let output = self.run(extra);
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(output.stdout.is_empty());
        String::from_utf8(output.stderr).unwrap()
    }
    fn failure(&self, extra: &[&str]) -> String {
        let output = self.run(extra);
        assert!(!output.status.success());
        assert!(output.stdout.is_empty());
        String::from_utf8(output.stderr).unwrap()
    }
    fn output(&self, suffix: &str) -> String {
        fs::read_to_string(sidecar(&self.prefix, suffix)).unwrap()
    }
    fn rows(&self) -> Vec<BTreeMap<String, String>> {
        let text = self.output(".annotated.tsv");
        let mut lines = text.lines();
        let headers: Vec<_> = lines.next().unwrap().split('\t').collect();
        lines
            .map(|line| {
                let values: Vec<_> = line.split('\t').collect();
                assert_eq!(headers.len(), values.len());
                headers
                    .iter()
                    .zip(values)
                    .map(|(h, v)| (h.to_string(), v.to_string()))
                    .collect()
            })
            .collect()
    }
    fn assert_no_outputs(&self) {
        assert!(!sidecar(&self.prefix, ".annotated.tsv").exists());
        assert!(!sidecar(&self.prefix, ".reverse").exists());
    }
    fn index(&self) {
        let index = noodles_fasta::fs::index(&self.fasta).unwrap();
        noodles_fasta::fai::fs::write(sidecar(&self.fasta, ".fai"), &index).unwrap();
    }
}

fn gzip(path: &Path, members: &[&[u8]]) {
    let mut bytes = Vec::new();
    for member in members {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(member).unwrap();
        bytes.extend(encoder.finish().unwrap());
    }
    fs::write(path, bytes).unwrap();
}

#[test]
fn supplied_example_has_safe_lists_and_explicit_categories() {
    let f = Fixture::new(
        include_str!("../examples/example.fa"),
        include_str!("../examples/example.bim"),
    );
    f.success(&["--chromosome-set", "human"]);
    assert_eq!(f.output(".reverse"), "snp3\nsnp4\n");
    assert_eq!(f.output(".forward"), "snp1\nsnp6\n");
    assert_eq!(f.output(".non_snp"), "esv5\n");
    assert_eq!(f.output(".mismatch"), "snp2\n");
    assert_eq!(f.output(".missing"), "inv1\nsnp7\n");
    assert_eq!(
        f.output(".ref_allele"),
        "snp1\tA\nsnp3\tT\nsnp4\tC\nsnp6\tA\n"
    );
    let rows = f.rows();
    assert_eq!(rows.len(), 8);
    assert_eq!(rows[1]["palindromic"], "true");
    assert_eq!(rows[1]["reason"], "palindromic;reference_mismatch");
    assert_eq!(rows[7]["reference_chromosome"], "chrX");
}

#[test]
fn preserves_duplicates_order_and_all_original_fields() {
    let bim = "2 . 0 3 A G\n1 . 0 3 A G\n1 same 0 1 A C\n1 same 0 1 A C\n1 other 0 1 A C\n";
    let f = Fixture::new(">chr1\nACG\n>chr2\nACG\n", bim);
    f.success(&[]);
    let rows = f.rows();
    assert_eq!(rows.len(), 5);
    for (i, (row, input)) in rows.iter().zip(bim.lines()).enumerate() {
        let fields: Vec<_> = input.split_whitespace().collect();
        for (key, expected) in [
            "chromosome",
            "snp_name",
            "genetic_distance",
            "position",
            "allele_1",
            "allele_2",
        ]
        .iter()
        .zip(fields)
        {
            assert_eq!(row[*key], expected);
        }
        assert_eq!(row["row"], (i + 1).to_string());
        assert_eq!(row["schema_version"], "1");
    }
    assert_eq!(f.output(".forward"), "other\n");
}

#[test]
fn ids_are_literal_strings_including_numeric_and_na_tokens() {
    for ids in [["00123", "00456"], ["NA", "null"]] {
        let f = Fixture::new(
            ">1\nAG\n",
            &format!("1 {} 0 1 A C\n1 {} 0 2 A C\n", ids[0], ids[1]),
        );
        f.success(&[]);
        assert_eq!(f.output(".forward"), format!("{}\n", ids[0]));
        assert_eq!(f.output(".reverse"), format!("{}\n", ids[1]));
    }
}

#[test]
fn repeated_and_missing_ids_never_generate_actions() {
    let f = Fixture::new(
        ">1\nAGC\n",
        "1 dup 0 1 A C\n1 dup 0 2 A C\n1 dup 0 3 A T\n1 . 0 2 A C\n",
    );
    let summary = f.success(&[]);
    for suffix in [".forward", ".reverse", ".ambiguous", ".ref_allele"] {
        assert_eq!(f.output(suffix), "");
    }
    assert!(summary.contains("withheld from ID lists and reference assignments: 4"));
    assert_eq!(f.rows()[0]["id_status"], "duplicate");
    assert_eq!(f.rows()[3]["id_status"], "missing");
}

#[test]
fn invalid_positions_are_missing_without_negative_indexing() {
    let f = Fixture::new(
        ">1\nACG\n",
        "1 zero 0 0 A G\n1 negative 0 -1 A G\n1 past 0 4 A G\n1 last 0 3 A G\n",
    );
    f.success(&[]);
    let rows = f.rows();
    for (row, reason) in rows
        .iter()
        .zip(["zero_position", "negative_position", "out_of_range"])
    {
        assert_eq!(row["category"], "missing");
        assert_eq!(row["reason"], reason);
        assert_eq!(row["reference"], "");
    }
    assert_eq!(rows[3]["reference"], "G");
    assert_eq!(f.output(".forward"), "last\n");
}

#[test]
fn iupac_reference_codes_are_unresolved_and_have_complements() {
    let bases = "NRYSWKMBDHV";
    let bim: String = bases
        .chars()
        .enumerate()
        .map(|(i, _)| format!("1 id{i} 0 {} A C\n", i + 1))
        .collect();
    let f = Fixture::new(&format!(">1\n{bases}\n"), &bim);
    f.success(&[]);
    let rows = f.rows();
    assert_eq!(rows[0]["reason"], "reference_n");
    assert_eq!(rows[1]["reference_rev"], "Y");
    assert!(rows.iter().all(|r| r["category"] == "missing"));
    assert_eq!(f.output(".reverse"), "");
}

#[test]
fn lowercase_alleles_and_reference_are_normalized_only_for_comparison() {
    let f = Fixture::new(">1\nag\n", "1 forward 0.00 0001 a c\n1 reverse 0 2 a c\n");
    f.success(&[]);
    assert_eq!(f.output(".reverse"), "reverse\n");
    assert_eq!(f.output(".forward"), "forward\n");
    let row = &f.rows()[0];
    assert_eq!(row["allele_1"], "a");
    assert_eq!(row["position"], "0001");
    assert_eq!(row["genetic_distance"], "0.00");
    assert_eq!(row["reference"], "A");
}

#[test]
fn unsupported_missing_and_invalid_alleles_are_never_flipped() {
    let f = Fixture::new(
        ">1\nA\n",
        "1 indel 0 1 AA G\n1 missing 0 1 0 A\n1 unknown 0 1 N A\n1 invalid 0 1 A Z\n1 mono 0 1 A A\n1 symbolic 0 1 I D\n",
    );
    f.success(&[]);
    assert_eq!(
        f.rows()
            .iter()
            .map(|r| r["category"].as_str())
            .collect::<Vec<_>>(),
        [
            "non_snp", "missing", "missing", "invalid", "forward", "non_snp"
        ]
    );
    assert_eq!(f.output(".forward"), "mono\n");
    assert_eq!(f.output(".reverse"), "");
    assert_eq!(
        f.output(".ambiguous"),
        "indel\nmissing\nunknown\ninvalid\nsymbolic\n"
    );
}

#[test]
fn palindromes_keep_mismatch_and_missing_reference_information() {
    let f = Fixture::new(
        ">1\nACN\n",
        "1 pal 0 1 A T\n1 mismatch 0 2 A T\n1 missing 0 3 A T\n",
    );
    f.success(&[]);
    let rows = f.rows();
    assert_eq!(rows[0]["category"], "palindromic");
    assert_eq!(rows[1]["category"], "mismatch");
    assert_eq!(rows[2]["category"], "missing");
    assert!(rows.iter().all(|r| r["palindromic"] == "true"));
    assert!(rows[1]["reason"].contains("reference_mismatch"));
    assert!(rows[2]["reason"].contains("reference_n"));
    assert_eq!(f.output(".reverse"), "");
    assert_eq!(f.output(".ref_allele"), "");
}

#[test]
fn malformed_bim_fails_with_line_context_and_no_outputs() {
    for (bim, message) in [
        ("1 x 0 1 A\n", "expected six BIM fields"),
        ("1 x 0 1 A C extra\n", "expected six BIM fields"),
        ("1 x 0 1.5 A C\n", "position must be an integer"),
        (
            "1 x 0 999999999999999999999 A C\n",
            "position must be an integer",
        ),
        ("1 x NaN 1 A C\n", "genetic distance must be finite"),
        ("", "BIM contains no records"),
    ] {
        let f = Fixture::new(">1\nA\n", bim);
        assert!(f.failure(&[]).contains(message));
        f.assert_no_outputs();
    }
}

#[test]
fn no_overlap_is_a_nonzero_error_on_stderr() {
    let f = Fixture::new(">1\nA\n", "2 id 0 1 A C\n");
    assert!(f.failure(&[]).contains("no chromosome overlap"));
    f.assert_no_outputs();
}

#[test]
fn arbitrary_contigs_are_not_truncated_or_ignored() {
    let f = Fixture::new(
        ">chr2L description\nA\n>chr2\nC\n>chrEBV\nG\n>chrI\nT\n>scaffold.1\nA\n",
        "2L fly 0 1 A C\n2 human 0 1 A C\nEBV virus 0 1 A C\nI yeast 0 1 A C\nscaffold.1 scaffold 0 1 A C\n",
    );
    f.success(&[]);
    assert_eq!(
        f.rows()
            .iter()
            .map(|r| r["reference"].as_str())
            .collect::<Vec<_>>(),
        ["A", "C", "G", "T", "A"]
    );
    assert_eq!(f.rows()[0]["reference_chromosome"], "chr2L");
}

#[test]
fn colliding_fasta_names_fail_even_when_unused() {
    for fasta in [
        ">1\nA\n>chr1\nC\n",
        ">1\nA\n>M\nC\n>MT\nG\n",
        ">1\nA\n>1\nC\n",
    ] {
        let f = Fixture::new(fasta, "1 id 0 1 A C\n");
        assert!(f.failure(&[]).contains("alias collision"));
        f.assert_no_outputs();
    }
}

#[test]
fn numeric_chromosome_aliases_require_the_human_preset() {
    let bim = "1 a 0 1 A C\n23 x 0 1 A C\n24 y 0 1 A C\n25 xy 0 1 A C\n26 mt 0 1 A C\nPAR1 p1 0 1 A C\nPAR2 p2 0 1 A C\nXY xytext 0 1 A C\n";
    let f = Fixture::new(">chr1\nA\n>chrX\nA\n>chrY\nG\n>chrM\nT\n", bim);
    f.success(&[]);
    assert_eq!(f.rows()[1]["reason"], "missing_chromosome");
    f.success(&["--chromosome-set", "human", "--force"]);
    assert_eq!(
        f.rows()
            .iter()
            .map(|r| r["reference_chromosome"].as_str())
            .collect::<Vec<_>>(),
        [
            "chr1", "chrX", "chrY", "chrX", "chrM", "chrX", "chrX", "chrX"
        ]
    );
    let autosome = Fixture::new(">chr23\nA\n", "23 cattle 0 1 A C\n");
    autosome.success(&[]);
    assert_eq!(autosome.rows()[0]["reference_chromosome"], "chr23");
}

#[test]
fn explicit_chromosome_mapping_is_validated_and_applies_to_bim_only() {
    let f = Fixture::new(">NC_000001.11\nA\n", "1 id 0 1 A C\n");
    let map = f._root.path().join("map.tsv");
    fs::write(&map, "# BIM FASTA\n1 NC_000001.11\n").unwrap();
    f.success(&["--chromosome-map", map.to_str().unwrap()]);
    assert_eq!(f.rows()[0]["reference_chromosome"], "NC_000001.11");
    fs::write(&map, "1 NC_000001.11\n1 different\n").unwrap();
    assert!(
        f.failure(&["--chromosome-map", map.to_str().unwrap(), "--force"])
            .contains("duplicate alias")
    );
}

#[test]
fn gzip_and_concatenated_members_work_without_sidecars_or_extensions() {
    let f = Fixture::new("", "");
    gzip(&f.fasta, &[b">1\nA", b"G\n>unused\nN\n"]);
    gzip(&f.bim, &[b"1 first 0 1 A C\n", b"1 second 0 2 A C\n"]);
    let summary = f.success(&[]);
    assert!(summary.contains("Reference access: stream"));
    assert_eq!(f.output(".reverse"), "second\n");
    assert!(!sidecar(&f.fasta, ".fai").exists());
    assert!(!sidecar(&f.fasta, ".gzi").exists());
}

#[test]
fn stream_and_indexed_outputs_agree_across_windows_and_line_boundaries() {
    let bases: String = (0..140000).map(|i| char::from(b"ACGT"[i % 4])).collect();
    let mut fasta = ">chr1 description\r\n".to_owned();
    for chunk in bases.as_bytes().chunks(61) {
        fasta.push_str(std::str::from_utf8(chunk).unwrap());
        fasta.push_str("\r\n");
    }
    fasta.push_str(">chr2\r\nA\r\n");
    let bim: String = [140000, 1, 65536, 65537, 2, 4096, 4100, 65536, 140001, 0, -1]
        .into_iter()
        .enumerate()
        .map(|(i, p)| format!("1 id{i} 0 {p} A C\n"))
        .collect();
    let f = Fixture::new(&fasta, &bim);
    f.success(&["--reference-mode", "stream"]);
    let expected = f.output(".annotated.tsv");
    f.index();
    let summary = f.success(&["--force"]);
    assert!(summary.contains("Reference access: indexed"));
    assert_eq!(f.output(".annotated.tsv"), expected);
}

#[test]
fn bgzf_indexed_and_streaming_access_agree_across_compressed_blocks() {
    let bases: String = (0..100000).map(|i| char::from(b"ACGT"[i % 4])).collect();
    let text = format!(">chr1\n{bases}\n");
    let f = Fixture::new(
        &text,
        "1 last 0 100000 A C\n1 first 0 1 A C\n1 middle 0 50000 A C\n",
    );
    f.index();
    let mut writer = noodles_bgzf::io::Writer::new(Vec::new());
    let mut pairs = Vec::new();
    let mut offset = 0;
    for chunk in text.as_bytes().chunks(32000) {
        if offset > 0 {
            pairs.push((writer.position(), offset));
        }
        writer.write_all(chunk).unwrap();
        writer.flush().unwrap();
        offset += chunk.len() as u64;
    }
    fs::write(&f.fasta, writer.finish().unwrap()).unwrap();
    noodles_bgzf::gzi::fs::write(
        sidecar(&f.fasta, ".gzi"),
        &noodles_bgzf::gzi::Index::from(pairs),
    )
    .unwrap();
    let summary = f.success(&[]);
    assert!(summary.contains("Reference access: indexed"));
    let expected = f.output(".annotated.tsv");
    f.success(&["--reference-mode", "stream", "--force"]);
    assert_eq!(f.output(".annotated.tsv"), expected);
}

#[test]
fn ordinary_gzip_with_only_a_fai_uses_streaming() {
    let f = Fixture::new(">1\nA\n", "1 id 0 1 A C\n");
    f.index();
    gzip(&f.fasta, &[b">1\nA\n"]);
    assert!(f.success(&[]).contains("Reference access: stream"));
    assert!(
        f.failure(&["--reference-mode", "indexed", "--force"])
            .contains("BGZF")
    );
}

#[test]
fn custom_fai_path_and_missing_or_short_indexes() {
    let f = Fixture::new(">1\nACG\n", "1 id 0 3 A C\n");
    assert!(
        f.failure(&["--reference-mode", "indexed"])
            .contains("reading FASTA index")
    );
    f.index();
    let custom = f._root.path().join("index elsewhere.fai");
    fs::rename(sidecar(&f.fasta, ".fai"), &custom).unwrap();
    assert!(
        f.success(&["--fasta-index", custom.to_str().unwrap()])
            .contains("Reference access: indexed")
    );
    fs::write(&f.fasta, ">1\nA\n").unwrap();
    assert!(
        f.failure(&["--fasta-index", custom.to_str().unwrap(), "--force"])
            .contains("short FASTA query")
    );
    assert!(
        f.failure(&[
            "--fasta-index",
            custom.to_str().unwrap(),
            "--reference-mode",
            "stream"
        ])
        .contains("cannot be combined")
    );
}

#[test]
fn corrupt_compression_and_malformed_fasta_leave_no_outputs() {
    for fasta in ["", "ACGT\n", ">1\nAC GT\n", ">1\nACZ\n"] {
        let f = Fixture::new(fasta, "1 id 0 1 A C\n");
        f.failure(&[]);
        f.assert_no_outputs();
    }
    let f = Fixture::new(">1\nACGT\n>2\nACGTA\nCZ\n", "1 id 0 1 A C\n");
    assert!(f.failure(&[]).contains("invalid FASTA base 'Z' at 2:7"));
    let f = Fixture::new("", "1 id 0 1 A C\n");
    gzip(&f.fasta, &[b">1\nA\n"]);
    let bytes = fs::read(&f.fasta).unwrap();
    fs::write(&f.fasta, &bytes[..bytes.len() - 5]).unwrap();
    f.failure(&[]);
    f.assert_no_outputs();
}

#[test]
fn human_preset_does_not_rename_fasta_contigs() {
    let f = Fixture::new(">X\nA\n>PAR1\nC\n", "25 par 0 1 A G\nX x 0 1 A G\n");
    f.success(&["--chromosome-set", "human"]);
    let rows = f.rows();
    assert_eq!(rows[0]["reference_chromosome"], "X");
    assert_eq!(rows[1]["reference_chromosome"], "X");
}

#[test]
fn outputs_are_checked_before_inputs_are_read() {
    let f = Fixture::new(">1\nA\n", "1 id 0 1 A C\n");
    f.success(&[]);
    fs::write(&f.bim, "malformed\n").unwrap();
    assert!(f.failure(&[]).contains("already exists"));
}

#[test]
fn paths_are_literal_and_existing_outputs_require_force() {
    let mut f = Fixture::new(">1\nA\n", "1 id 0 1 A C\n");
    f.prefix = f
        ._root
        .path()
        .join("space ; $(echo injected)")
        .join("result");
    f.success(&[]);
    let original = f.output(".annotated.tsv");
    assert!(f.failure(&[]).contains("already exists"));
    assert_eq!(f.output(".annotated.tsv"), original);
    fs::write(&f.bim, "1 changed 0 1 A G\n").unwrap();
    f.success(&["--force"]);
    assert_eq!(f.output(".forward"), "changed\n");
    let updated = f.output(".annotated.tsv");
    fs::write(&f.bim, "broken\n").unwrap();
    f.failure(&["--force"]);
    assert_eq!(f.output(".annotated.tsv"), updated);
}

#[test]
fn output_preflight_prevents_partial_results_and_input_overwrites() {
    let mut f = Fixture::new(">1\nA\n", "1 id 0 1 A C\n");
    fs::create_dir(sidecar(&f.prefix, ".missing")).unwrap();
    assert!(f.failure(&["--force"]).contains("non-regular output"));
    f.assert_no_outputs();
    fs::remove_dir(sidecar(&f.prefix, ".missing")).unwrap();
    let protected = sidecar(&f.prefix, ".reverse");
    fs::rename(&f.bim, &protected).unwrap();
    f.bim = protected;
    assert!(f.failure(&["--force"]).contains("would replace an input"));
    assert_eq!(fs::read_to_string(&f.bim).unwrap(), "1 id 0 1 A C\n");
}

#[test]
fn help_version_and_default_prefix() {
    let output = Command::new(env!("CARGO_BIN_EXE_snpflip"))
        .arg("--help")
        .output()
        .unwrap();
    assert!(output.status.success());
    assert!(String::from_utf8_lossy(&output.stdout).contains("--chromosome-set"));
    let output = Command::new(env!("CARGO_BIN_EXE_snpflip"))
        .arg("--version")
        .output()
        .unwrap();
    assert!(output.status.success());
    assert_eq!(String::from_utf8(output.stdout).unwrap(), "snpflip 0.1.0\n");
    for short in ["-v", "-V"] {
        assert!(
            Command::new(env!("CARGO_BIN_EXE_snpflip"))
                .arg(short)
                .output()
                .unwrap()
                .status
                .success()
        );
    }
    let f = Fixture::new(">1\nA\n", "1 id 0 1 A C\n");
    let output = Command::new(env!("CARGO_BIN_EXE_snpflip"))
        .current_dir(f._root.path())
        .arg("-f")
        .arg(&f.fasta)
        .arg("-b")
        .arg(&f.bim)
        .output()
        .unwrap();
    assert!(output.status.success());
    assert!(
        f._root
            .path()
            .join("snpflip_output/input.annotated.tsv")
            .exists()
    );
    assert_eq!(
        snpflip::default_prefix(Path::new("/data/sample.bim.gz")),
        PathBuf::from("snpflip_output/sample")
    );
}

#[test]
fn blank_lines_keep_physical_line_numbers_and_runtime_exit_codes() {
    let f = Fixture::new(">1\nA\n", "\n1 valid 0 1 A C\n\n2 unknown 0 0 A T\n");
    f.success(&[]);
    let rows = f.rows();
    assert_eq!(rows[0]["input_line"], "2");
    assert_eq!(rows[1]["input_line"], "4");
    assert_eq!(rows[1]["row"], "2");
    assert!(rows[1]["reason"].contains("zero_position"));
    assert!(rows[1]["reason"].contains("missing_chromosome"));
    assert_eq!(f.run(&[]).status.code(), Some(1));
    assert_eq!(
        Command::new(env!("CARGO_BIN_EXE_snpflip"))
            .arg("--unknown")
            .output()
            .unwrap()
            .status
            .code(),
        Some(2)
    );
}

#[test]
fn indexed_catalog_collisions_and_offset_overflow_are_errors() {
    let f = Fixture::new(">1\nA\n>chr1\nC\n", "1 id 0 1 A C\n");
    f.index();
    assert!(f.failure(&[]).contains("alias collision"));
    fs::write(
        sidecar(&f.fasta, ".fai"),
        "1\t3\t18446744073709551615\t1\t2\n",
    )
    .unwrap();
    assert!(f.failure(&[]).contains("offset overflow"));
    f.assert_no_outputs();
}

#[cfg(unix)]
#[test]
fn readonly_reference_directories_work_and_output_symlinks_are_rejected() {
    use std::os::unix::fs::{PermissionsExt, symlink};
    let mut f = Fixture::new(">1\nA\n", "1 id 0 1 A C\n");
    let directory = f._root.path().join("readonly");
    fs::create_dir(&directory).unwrap();
    let reference = directory.join("ref.fa");
    fs::rename(&f.fasta, &reference).unwrap();
    f.fasta = reference;
    fs::set_permissions(&directory, fs::Permissions::from_mode(0o555)).unwrap();
    let result = f.run(&[]);
    fs::set_permissions(&directory, fs::Permissions::from_mode(0o755)).unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(!sidecar(&f.fasta, ".fai").exists());
    let reverse = sidecar(&f.prefix, ".reverse");
    fs::remove_file(&reverse).unwrap();
    symlink(&f.bim, &reverse).unwrap();
    assert!(f.failure(&["--force"]).contains("non-regular output"));
    assert_eq!(fs::read_to_string(&f.bim).unwrap(), "1 id 0 1 A C\n");
}

#[cfg(unix)]
fn run_with_stdin(args: &[&std::ffi::OsStr], stdin: &[u8]) -> Output {
    use std::process::Stdio;
    let mut child = Command::new(env!("CARGO_BIN_EXE_snpflip"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let mut input = child.stdin.take().unwrap();
    let bytes = stdin.to_vec();
    let writer = std::thread::spawn(move || input.write_all(&bytes));
    let output = child.wait_with_output().unwrap();
    writer.join().unwrap().unwrap();
    output
}

#[cfg(unix)]
#[test]
fn piped_inputs_match_regular_files() {
    let bim = "1 dup 0 1 A C\n1 b 0 2 A C\n\n2 dup 0 1 A G\n1 c 0 3 G T\n";
    let f = Fixture::new(">1\nACG\n>2\nT\n", bim);
    f.success(&[]);
    let expected: Vec<_> = [".annotated.tsv", ".reverse", ".forward", ".ref_allele"]
        .iter()
        .map(|suffix| f.output(suffix))
        .collect();
    let mut compressed = GzEncoder::new(Vec::new(), Compression::default());
    compressed.write_all(bim.as_bytes()).unwrap();
    let compressed = compressed.finish().unwrap();
    for (name, bim_arg, fasta_arg, stdin) in [
        (
            "dash",
            "-",
            f.fasta.to_str().unwrap(),
            bim.as_bytes().to_vec(),
        ),
        ("pipe", "/dev/stdin", f.fasta.to_str().unwrap(), compressed),
        (
            "fasta",
            f.bim.to_str().unwrap(),
            "-",
            b">1\nACG\n>2\nT\n".to_vec(),
        ),
    ] {
        let prefix = f._root.path().join(name).join("out");
        let output = run_with_stdin(
            &[
                "-b".as_ref(),
                bim_arg.as_ref(),
                "-f".as_ref(),
                fasta_arg.as_ref(),
                "-o".as_ref(),
                prefix.as_os_str(),
            ],
            &stdin,
        );
        assert!(
            output.status.success(),
            "{name}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        for (suffix, expected) in [".annotated.tsv", ".reverse", ".forward", ".ref_allele"]
            .iter()
            .zip(&expected)
        {
            assert_eq!(
                &fs::read_to_string(sidecar(&prefix, suffix)).unwrap(),
                expected,
                "{name} {suffix}"
            );
        }
    }
}

#[cfg(unix)]
#[test]
fn stdin_misuse_is_reported() {
    let f = Fixture::new(">1\nA\n", "1 id 0 1 A C\n");
    let output = run_with_stdin(
        &["-b".as_ref(), "-".as_ref(), "-f".as_ref(), "-".as_ref()],
        b"",
    );
    assert!(String::from_utf8_lossy(&output.stderr).contains("only one input"));
    let output = run_with_stdin(
        &[
            "-b".as_ref(),
            f.bim.as_os_str(),
            "-f".as_ref(),
            "-".as_ref(),
            "-o".as_ref(),
            f.prefix.as_os_str(),
            "--reference-mode".as_ref(),
            "indexed".as_ref(),
        ],
        b">1\nA\n",
    );
    assert!(String::from_utf8_lossy(&output.stderr).contains("needs a regular file"));
    assert_eq!(
        snpflip::default_prefix(Path::new("-")),
        PathBuf::from("snpflip_output/stdin")
    );
}
