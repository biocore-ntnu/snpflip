use std::{collections::HashMap, path::Path};

use anyhow::{Context, Result, ensure};
use clap::ValueEnum;

use crate::input::open_text;

#[derive(Clone, Copy, Debug, Default, ValueEnum)]
pub enum ChromosomeSet {
    #[default]
    Generic,
    Human,
}

pub struct Aliases {
    preset: ChromosomeSet,
    custom: HashMap<String, String>,
}

impl Aliases {
    pub fn new(preset: ChromosomeSet, path: Option<&Path>) -> Result<Self> {
        let mut custom = HashMap::new();
        if let Some(path) = path {
            let mut reader = open_text(path)?;
            let mut text = String::new();
            let mut line = 0;
            loop {
                text.clear();
                if reader
                    .read_line(&mut text)
                    .with_context(|| format!("reading {}", path.display()))?
                    == 0
                {
                    break;
                }
                line += 1;
                if text.trim().is_empty() || text.trim_start().starts_with('#') {
                    continue;
                }
                let fields: Vec<_> = text.split_whitespace().collect();
                ensure!(
                    fields.len() == 2,
                    "{}:{line}: alias map requires BIM chromosome and FASTA name",
                    path.display()
                );
                ensure!(
                    !custom.contains_key(fields[0]),
                    "{}:{line}: duplicate alias {:?}",
                    path.display(),
                    fields[0]
                );
                custom.insert(fields[0].to_owned(), fields[1].to_owned());
            }
        }
        Ok(Self { preset, custom })
    }

    /// An explicit mapping wins; otherwise the preset's PLINK chromosome codes
    /// apply. Those codes describe BIM files only and never rename FASTA contigs.
    pub fn bim_name(&self, name: &str) -> String {
        if let Some(target) = self.custom.get(name) {
            return self.fasta_name(target);
        }
        let name = normalize(name);
        match (self.preset, name) {
            (ChromosomeSet::Human, "23" | "25" | "XY" | "PAR1" | "PAR2") => "X",
            (ChromosomeSet::Human, "24") => "Y",
            (ChromosomeSet::Human, "26") => "MT",
            _ => name,
        }
        .to_owned()
    }

    pub fn fasta_name(&self, name: &str) -> String {
        normalize(name).to_owned()
    }
}

fn normalize(name: &str) -> &str {
    let name = name
        .strip_prefix("chr")
        .filter(|s| !s.is_empty())
        .unwrap_or(name);
    if matches!(name, "M" | "MT") {
        "MT"
    } else {
        name
    }
}
