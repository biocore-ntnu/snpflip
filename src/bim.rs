use std::{
    collections::{HashMap, HashSet},
    fs::File,
    hash::{DefaultHasher, Hash, Hasher},
    io::{BufRead, BufReader, BufWriter, Seek, SeekFrom, Write},
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, ensure};

use crate::{
    aliases::Aliases,
    classify::IdStatus,
    input::{is_regular, open_text},
};

/// One BIM line, borrowing its original tokens (ID spelling, allele case, and
/// coordinate spelling are preserved).
pub struct Record<'a> {
    pub fields: [&'a str; 6],
    pub position: i64,
    pub line: usize,
}

impl Record<'_> {
    pub fn chromosome(&self) -> &str {
        self.fields[0]
    }

    pub fn id(&self) -> &str {
        self.fields[1]
    }

    pub fn alleles(&self) -> [&str; 2] {
        [self.fields[4], self.fields[5]]
    }
}

/// The per-row data kept in memory; everything else is re-read from the BIM.
#[derive(Clone, Copy)]
pub struct Row {
    pub position: i64,
    pub contig: u32,
}

pub struct Index {
    pub rows: Vec<Row>,
    /// Normalized chromosome keys, indexed by `Row::contig`.
    pub contigs: Vec<String>,
    pub contig_ids: HashMap<String, u32>,
    raw: HashMap<String, u32>,
    duplicates: HashSet<String>,
    path: PathBuf,
    /// Decompressed copy of a BIM that cannot be reopened (a pipe or stdin).
    spool: Option<File>,
}

impl Index {
    /// Call `f` for every record again, from the original file or the spool.
    pub fn for_each(&self, f: impl FnMut(Record) -> Result<()>) -> Result<()> {
        let reader: Box<dyn BufRead> = match &self.spool {
            Some(file) => {
                let mut file = file.try_clone()?;
                file.seek(SeekFrom::Start(0))?;
                Box::new(BufReader::new(file))
            }
            None => open_text(&self.path)?,
        };
        scan(reader, &self.path, None, f)
    }

    fn contig(&mut self, raw: &str, aliases: &Aliases) -> u32 {
        if let Some(&id) = self.raw.get(raw) {
            return id;
        }
        let key = aliases.bim_name(raw);
        let id = match self.contig_ids.get(&key) {
            Some(&id) => id,
            None => {
                let id = self.contigs.len() as u32;
                self.contig_ids.insert(key.clone(), id);
                self.contigs.push(key);
                id
            }
        };
        self.raw.insert(raw.to_owned(), id);
        id
    }

    pub fn id_status(&self, id: &str) -> IdStatus {
        if id == "." {
            IdStatus::Missing
        } else if self.duplicates.contains(id) {
            IdStatus::Duplicate
        } else {
            IdStatus::Unique
        }
    }

    /// Fail if a later pass sees a different BIM than the one indexed.
    pub fn check(&self, row: usize, record: &Record) -> Result<()> {
        let same = self.rows.get(row).is_some_and(|r| {
            r.position == record.position && self.raw.get(record.chromosome()) == Some(&r.contig)
        });
        ensure!(
            same,
            "BIM line {} changed while it was being read",
            record.line
        );
        Ok(())
    }
}

fn hash(id: &str) -> u64 {
    let mut hasher = DefaultHasher::new();
    id.hash(&mut hasher);
    hasher.finish()
}

/// Call `f` for every non-blank BIM line, validating its fields, and copy each
/// line to `copy` if given.
fn scan(
    mut reader: Box<dyn BufRead>,
    path: &Path,
    mut copy: Option<&mut BufWriter<File>>,
    mut f: impl FnMut(Record) -> Result<()>,
) -> Result<()> {
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
        if let Some(copy) = copy.as_mut() {
            copy.write_all(text.as_bytes())
                .context("writing temporary copy of BIM; set TMPDIR to a larger disk")?;
        }
        let mut fields = [""; 6];
        let mut count = 0;
        for token in text.split_whitespace() {
            if count < 6 {
                fields[count] = token;
            }
            count += 1;
        }
        if count == 0 {
            continue;
        }
        ensure!(
            count == 6,
            "{}:{line}: expected six BIM fields, found {count}",
            path.display()
        );
        let position = fields[3].parse::<i64>().with_context(|| {
            format!(
                "{}:{line}: position must be an integer, found {:?}",
                path.display(),
                fields[3]
            )
        })?;
        let distance = fields[2].parse::<f64>().with_context(|| {
            format!(
                "{}:{line}: invalid genetic distance {:?}",
                path.display(),
                fields[2]
            )
        })?;
        ensure!(
            distance.is_finite(),
            "{}:{line}: genetic distance must be finite",
            path.display()
        );
        f(Record {
            fields,
            position,
            line,
        })?;
    }
    Ok(())
}

/// Index the BIM compactly. IDs are only hashed; a second pass compares the
/// actual strings behind repeated hashes, so a collision never marks an ID
/// as duplicate. A BIM that cannot be reopened is spooled to an unnamed
/// temporary file while it is read.
pub fn index(path: &Path, aliases: &Aliases) -> Result<Index> {
    let spool = if is_regular(path) {
        None
    } else {
        Some(tempfile::tempfile().context("creating temporary copy of BIM")?)
    };
    let mut copy = spool
        .as_ref()
        .map(|file| file.try_clone().map(BufWriter::new))
        .transpose()?;
    let mut index = Index {
        rows: Vec::new(),
        contigs: Vec::new(),
        contig_ids: HashMap::new(),
        raw: HashMap::new(),
        duplicates: HashSet::new(),
        path: path.to_owned(),
        spool: None,
    };
    let mut hashes = Vec::new();
    scan(open_text(path)?, path, copy.as_mut(), |record| {
        let contig = index.contig(record.chromosome(), aliases);
        index.rows.push(Row {
            position: record.position,
            contig,
        });
        if record.id() != "." {
            hashes.push(hash(record.id()));
        }
        Ok(())
    })?;
    if let Some(copy) = copy {
        copy.into_inner()
            .map_err(|e| e.into_error())
            .context("writing temporary copy of BIM; set TMPDIR to a larger disk")?;
    }
    index.spool = spool;
    ensure!(
        !index.rows.is_empty(),
        "{}: BIM contains no records",
        path.display()
    );
    ensure!(
        u32::try_from(index.rows.len()).is_ok(),
        "{}: more than {} BIM records",
        path.display(),
        u32::MAX
    );
    hashes.sort_unstable();
    let candidates: Vec<u64> = hashes
        .chunk_by(|a, b| a == b)
        .filter(|chunk| chunk.len() > 1)
        .map(|chunk| chunk[0])
        .collect();
    drop(hashes);
    if !candidates.is_empty() {
        let mut counts: HashMap<String, usize> = HashMap::new();
        let mut row = 0;
        index.for_each(|record| {
            index.check(row, &record)?;
            row += 1;
            if record.id() != "." && candidates.binary_search(&hash(record.id())).is_ok() {
                *counts.entry(record.id().to_owned()).or_default() += 1;
            }
            Ok(())
        })?;
        index.duplicates = counts
            .into_iter()
            .filter(|(_, n)| *n > 1)
            .map(|(id, _)| id)
            .collect();
    }
    Ok(index)
}
