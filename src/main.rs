use std::{path::PathBuf, process::ExitCode};

use clap::Parser;
use snpflip::{Config, aliases::ChromosomeSet, reference::ReferenceMode};

#[derive(Parser)]
#[command(
    version,
    disable_version_flag = true,
    about = "Check PLINK BIM allele orientation against a reference genome",
    after_help = "Forward/reverse indicate compatibility with the supplied reference, not proof of genome build.\nPalindromic SNPs are never automatically flipped. Repeated and missing IDs are annotated but\nexcluded from ID-based outputs. See README.md for the output schema and PLINK examples."
)]
struct Args {
    #[arg(short = 'v', short_alias = 'V', long, action = clap::ArgAction::Version)]
    version: Option<bool>,
    /// Reference FASTA (plain, gzip, or BGZF; - for stdin)
    #[arg(short = 'f', long = "fasta-genome", value_name = "FA")]
    fasta: PathBuf,
    /// PLINK BIM variant information (plain or gzip; - for stdin)
    #[arg(short = 'b', long = "bim-file", value_name = "BIM")]
    bim: PathBuf,
    /// Output prefix [default: snpflip_output/<BIM basename>]
    #[arg(short = 'o', long = "output-prefix", value_name = "PREFIX")]
    prefix: Option<PathBuf>,
    /// Human additionally maps 23/24/25/26 and XY/PAR1/PAR2 to X/Y/X/MT
    #[arg(long, value_enum, default_value = "generic")]
    chromosome_set: ChromosomeSet,
    /// Two-column BIM chromosome -> FASTA name mapping (no header)
    #[arg(long, value_name = "TSV")]
    chromosome_map: Option<PathBuf>,
    /// Auto uses existing .fai (and .gzi for BGZF); no indexes are created
    #[arg(long, value_enum, default_value = "auto")]
    reference_mode: ReferenceMode,
    /// Use an existing FASTA index at a different path
    #[arg(long, value_name = "FAI")]
    fasta_index: Option<PathBuf>,
    /// Replace existing regular output files
    #[arg(long)]
    force: bool,
}

fn main() -> ExitCode {
    let args = Args::parse();
    let prefix = args
        .prefix
        .unwrap_or_else(|| snpflip::default_prefix(&args.bim));
    let config = Config {
        bim: args.bim,
        fasta: args.fasta,
        prefix,
        chromosome_set: args.chromosome_set,
        chromosome_map: args.chromosome_map,
        reference_mode: args.reference_mode,
        fasta_index: args.fasta_index,
        force: args.force,
    };
    match snpflip::run(&config) {
        Ok(summary) => {
            eprint!("{summary}");
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("snpflip: {error:#}");
            ExitCode::FAILURE
        }
    }
}
