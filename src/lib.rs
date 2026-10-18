//! Internal library behind the `snpflip` command-line tool. It exists for the
//! binary and its tests; the API is not stable and may change in any release.

pub mod aliases;
pub mod bim;
pub mod classify;
pub mod input;
pub mod output;
pub mod reference;

use std::path::PathBuf;

use anyhow::{Result, ensure};

use aliases::{Aliases, ChromosomeSet};
use reference::ReferenceMode;

pub struct Config {
    pub bim: PathBuf,
    pub fasta: PathBuf,
    pub prefix: PathBuf,
    pub chromosome_set: ChromosomeSet,
    pub chromosome_map: Option<PathBuf>,
    pub reference_mode: ReferenceMode,
    pub fasta_index: Option<PathBuf>,
    pub force: bool,
}

pub fn run(config: &Config) -> Result<String> {
    ensure!(
        !config.prefix.as_os_str().is_empty() && config.prefix.file_name().is_some(),
        "output prefix must include a file name"
    );
    if let Some(path) = &config.fasta_index {
        ensure!(
            !matches!(config.reference_mode, ReferenceMode::Stream),
            "--fasta-index cannot be combined with --reference-mode stream"
        );
        ensure!(
            path.is_file(),
            "FASTA index {} does not exist",
            path.display()
        );
    }
    let stdin_inputs = [&config.bim, &config.fasta]
        .into_iter()
        .chain(&config.chromosome_map)
        .filter(|path| input::is_stdin(path))
        .count();
    ensure!(
        stdin_inputs <= 1,
        "only one input can be read from standard input (-)"
    );
    let mut protected = vec![
        config.bim.clone(),
        config.fasta.clone(),
        reference::sidecar(&config.fasta, ".fai"),
        reference::sidecar(&config.fasta, ".gzi"),
    ];
    protected.extend(config.chromosome_map.iter().cloned());
    protected.extend(config.fasta_index.iter().cloned());
    // Fail before reading anything large if the outputs cannot be written.
    output::preflight(&config.prefix, config.force, &protected)?;
    let aliases = Aliases::new(config.chromosome_set, config.chromosome_map.as_deref())?;
    let bim = bim::index(&config.bim, &aliases)?;
    let reference = reference::read(
        &bim,
        &config.fasta,
        config.fasta_index.as_deref(),
        config.reference_mode,
        &aliases,
    )?;
    let summary = output::write(&config.prefix, config.force, &protected, &bim, &reference)?;
    Ok(summary.to_string())
}

pub fn default_prefix(bim: &std::path::Path) -> PathBuf {
    if input::is_stdin(bim) {
        return PathBuf::from("snpflip_output/stdin");
    }
    let mut name = bim.file_name().unwrap_or(bim.as_os_str()).to_os_string();
    if bim.extension().is_some_and(|s| s == "gz" || s == "bgz") {
        name = std::path::Path::new(&name)
            .file_stem()
            .unwrap_or(&name)
            .to_os_string();
    }
    let stem = std::path::Path::new(&name).file_stem().unwrap_or(&name);
    PathBuf::from("snpflip_output").join(stem)
}
