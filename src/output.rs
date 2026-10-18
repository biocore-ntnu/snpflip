use std::{
    collections::BTreeMap,
    fmt, fs,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail, ensure};
use tempfile::{NamedTempFile, TempPath};

use crate::{
    bim::{Index, Record},
    classify::{Annotation, Category, IdStatus, classify, complement},
    reference::{Lookup, ReferenceData, sidecar},
};

const ANNOTATED: usize = 0;
const AMBIGUOUS: usize = 1;
const REF_ALLELE: usize = 2;
const CATEGORIES: usize = 3;

fn destinations(prefix: &Path) -> Vec<PathBuf> {
    let mut suffixes = vec![
        ".annotated.tsv".to_owned(),
        ".ambiguous".to_owned(),
        ".ref_allele".to_owned(),
    ];
    suffixes.extend(Category::ALL.map(|c| format!(".{}", c.name())));
    suffixes.iter().map(|s| sidecar(prefix, s)).collect()
}

/// Refuse outputs that are not regular files, would replace an input, or
/// exist without `force`. Runs before any work and again before writing.
pub fn preflight(prefix: &Path, force: bool, protected: &[PathBuf]) -> Result<Vec<PathBuf>> {
    let destinations = destinations(prefix);
    let protected: Vec<_> = protected
        .iter()
        .filter_map(|p| fs::canonicalize(p).ok())
        .collect();
    for path in &destinations {
        match fs::symlink_metadata(path) {
            Ok(metadata) => {
                ensure!(
                    metadata.is_file(),
                    "refusing to replace non-regular output {}",
                    path.display()
                );
                ensure!(
                    !protected.contains(&fs::canonicalize(path)?),
                    "output {} would replace an input file",
                    path.display()
                );
                ensure!(
                    force,
                    "output {} already exists; use --force to replace outputs",
                    path.display()
                );
            }
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
            Err(error) => {
                return Err(error).with_context(|| format!("checking output {}", path.display()));
            }
        }
    }
    Ok(destinations)
}

/// Re-read the BIM and write every output in one streaming pass.
pub fn write(
    prefix: &Path,
    force: bool,
    protected: &[PathBuf],
    bim: &Index,
    reference: &ReferenceData,
) -> Result<Summary> {
    let destinations = preflight(prefix, force, protected)?;
    let parent = prefix
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    fs::create_dir_all(parent)
        .with_context(|| format!("creating output directory {}", parent.display()))?;
    let mut files = destinations
        .iter()
        .map(|_| NamedTempFile::new_in(parent).context("creating temporary output"))
        .collect::<Result<Vec<_>>>()?;
    let mut summary = Summary::new(reference);
    {
        let mut writers: Vec<_> = files
            .iter_mut()
            .map(|file| BufWriter::new(file.as_file_mut()))
            .collect();
        writeln!(
            writers[ANNOTATED],
            "chromosome\tsnp_name\tgenetic_distance\tposition\tallele_1\tallele_2\tschema_version\trow\tinput_line\treference_chromosome\treference\treference_rev\tstrand\tcategory\treason\tpalindromic\tid_status"
        )?;
        let mut row = 0;
        bim.for_each(|record| {
            bim.check(row, &record)?;
            let lookup = &reference.lookups[row];
            let name = reference.names[bim.rows[row].contig as usize].as_deref();
            let annotation = classify(
                record.alleles(),
                lookup,
                name.is_some(),
                bim.id_status(record.id()),
            );
            row += 1;
            write_row(&mut writers, row, &record, lookup, name, &annotation)?;
            summary.add(&annotation);
            Ok(())
        })?;
        ensure!(row == bim.rows.len(), "BIM changed while it was being read");
        for writer in writers {
            writer.into_inner().map_err(|e| e.into_error())?;
        }
    }
    for file in &files {
        file.as_file().sync_all()?;
    }
    publish(destinations.into_iter().zip(files).collect(), force)?;
    Ok(summary)
}

fn write_row(
    writers: &mut [BufWriter<&mut fs::File>],
    row: usize,
    record: &Record,
    lookup: &Lookup,
    name: Option<&str>,
    annotation: &Annotation,
) -> Result<()> {
    let reference = lookup
        .base
        .map(|b| char::from(b).to_string())
        .unwrap_or_default();
    let reverse = lookup
        .base
        .and_then(complement)
        .map(|b| char::from(b).to_string())
        .unwrap_or_default();
    let strand = if annotation.category.resolved() {
        annotation.category.name()
    } else {
        "ambiguous"
    };
    let [chromosome, id, distance, position, a, b] = record.fields;
    writeln!(
        writers[ANNOTATED],
        "{chromosome}\t{id}\t{distance}\t{position}\t{a}\t{b}\t1\t{row}\t{}\t{}\t{reference}\t{reverse}\t{strand}\t{}\t{}\t{}\t{}",
        record.line,
        name.unwrap_or(""),
        annotation.category.name(),
        annotation.reasons.join(";"),
        annotation.palindromic,
        annotation.id_status.name()
    )?;
    if annotation.id_status == IdStatus::Unique {
        writeln!(writers[CATEGORIES + annotation.category as usize], "{id}")?;
        if annotation.category.resolved() {
            let base = lookup.base.expect("resolved reference");
            writeln!(writers[REF_ALLELE], "{id}\t{}", char::from(base))?;
        } else {
            writeln!(writers[AMBIGUOUS], "{id}")?;
        }
    }
    Ok(())
}

/// All files are prepared before publishing. Restore previous outputs on an
/// ordinary commit error; individual renames are atomic, the whole set is not.
fn publish(staged: Vec<(PathBuf, NamedTempFile)>, force: bool) -> Result<()> {
    let mut backups: Vec<(PathBuf, TempPath)> = Vec::new();
    let mut published = Vec::new();
    let result: Result<()> = (|| {
        if force {
            for (path, _) in &staged {
                if let Ok(metadata) = fs::symlink_metadata(path) {
                    // Outputs may have changed since preflight; never move anything
                    // but a regular file (Windows can rename a directory over a file).
                    ensure!(
                        metadata.is_file(),
                        "refusing to replace non-regular output {}",
                        path.display()
                    );
                    let backup = NamedTempFile::new_in(path.parent().unwrap_or(Path::new(".")))?
                        .into_temp_path();
                    fs::rename(path, &backup)
                        .with_context(|| format!("backing up {}", path.display()))?;
                    backups.push((path.clone(), backup));
                }
            }
        }
        for (path, file) in staged {
            file.persist_noclobber(&path)
                .with_context(|| format!("publishing {}", path.display()))?;
            published.push(path);
        }
        Ok(())
    })();
    if let Err(error) = result {
        let mut failures = Vec::new();
        for path in published {
            if let Err(error) = fs::remove_file(&path) {
                failures.push(format!("removing {}: {error}", path.display()));
            }
        }
        for (path, backup) in backups {
            if let Err(error) = fs::rename(&backup, &path) {
                let backup_name = backup.display().to_string();
                let _ = backup.keep();
                failures.push(format!(
                    "restoring {}: {error}; backup retained at {backup_name}",
                    path.display()
                ));
            }
        }
        if !failures.is_empty() {
            bail!("{error:#}; rollback also failed: {}", failures.join("; "));
        }
        return Err(error);
    }
    Ok(())
}

pub struct Summary {
    access: &'static str,
    rows: usize,
    categories: [usize; Category::ALL.len()],
    reasons: BTreeMap<&'static str, usize>,
    missing: usize,
    duplicate: usize,
}

impl Summary {
    fn new(reference: &ReferenceData) -> Self {
        Self {
            access: reference.access,
            rows: 0,
            categories: [0; Category::ALL.len()],
            reasons: BTreeMap::new(),
            missing: 0,
            duplicate: 0,
        }
    }

    fn add(&mut self, annotation: &Annotation) {
        self.rows += 1;
        self.categories[annotation.category as usize] += 1;
        for &reason in &annotation.reasons {
            *self.reasons.entry(reason).or_default() += 1;
        }
        match annotation.id_status {
            IdStatus::Missing => self.missing += 1,
            IdStatus::Duplicate => self.duplicate += 1,
            IdStatus::Unique => {}
        }
    }
}

impl fmt::Display for Summary {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "Reference access: {}", self.access)?;
        writeln!(f, "Input rows: {}", self.rows)?;
        for category in Category::ALL {
            writeln!(
                f,
                "{}: {}",
                category.name(),
                self.categories[category as usize]
            )?;
        }
        for (reason, count) in &self.reasons {
            writeln!(f, "reason.{reason}: {count}")?;
        }
        writeln!(
            f,
            "Rows withheld from ID lists and reference assignments: {} (missing ID: {}, duplicate ID: {})",
            self.missing + self.duplicate,
            self.missing,
            self.duplicate
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn failed_publish_removes_new_files_and_restores_previous_outputs() {
        for force in [false, true] {
            let root = tempfile::tempdir().unwrap();
            let first = root.path().join("first");
            let blocked = root.path().join("blocked");
            if force {
                fs::write(&first, "previous").unwrap();
            }
            fs::create_dir(&blocked).unwrap();
            let mut temporary = NamedTempFile::new_in(root.path()).unwrap();
            temporary.write_all(b"replacement").unwrap();
            let another = NamedTempFile::new_in(root.path()).unwrap();
            assert!(
                publish(
                    vec![(first.clone(), temporary), (blocked.clone(), another)],
                    force
                )
                .is_err()
            );
            if force {
                assert_eq!(fs::read_to_string(&first).unwrap(), "previous");
            } else {
                assert!(!first.exists());
            }
            assert!(blocked.is_dir());
            assert_eq!(
                fs::read_dir(root.path()).unwrap().count(),
                if force { 2 } else { 1 }
            );
        }
    }
}
