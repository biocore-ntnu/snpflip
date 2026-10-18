use std::{
    collections::HashMap,
    fs::File,
    io::{BufReader, Read},
    ops::Range,
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, ensure};
use clap::ValueEnum;
use noodles_core::{Position, Region};
use noodles_fasta as fasta;

use crate::{
    aliases::Aliases,
    bim::Index,
    classify::complement,
    input::{is_gzip, is_regular, open_text},
};

#[derive(Clone, Copy, Debug, Default, ValueEnum)]
pub enum ReferenceMode {
    #[default]
    Auto,
    Stream,
    Indexed,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReferenceStatus {
    Found,
    MissingChromosome,
    ZeroPosition,
    NegativePosition,
    OutOfRange,
}

#[derive(Clone, Copy, Debug)]
pub struct Lookup {
    pub base: Option<u8>,
    pub status: ReferenceStatus,
}

pub struct ReferenceData {
    pub lookups: Vec<Lookup>,
    /// FASTA sequence name for each BIM contig the reference contains.
    pub names: Vec<Option<String>>,
    pub access: &'static str,
}

pub fn sidecar(path: &Path, extension: &str) -> PathBuf {
    let mut name = path.as_os_str().to_os_string();
    name.push(extension);
    name.into()
}

/// BIM rows grouped by contig and sorted by position, plus the lookups being filled.
struct Scan<'a> {
    bim: &'a Index,
    aliases: &'a Aliases,
    order: Vec<u32>,
    groups: Vec<Range<usize>>,
    catalog: HashMap<String, String>,
    lookups: Vec<Lookup>,
    names: Vec<Option<String>>,
    overlap: usize,
}

impl<'a> Scan<'a> {
    fn new(bim: &'a Index, aliases: &'a Aliases) -> Self {
        let lookups = bim
            .rows
            .iter()
            .map(|row| Lookup {
                base: None,
                status: match row.position {
                    0 => ReferenceStatus::ZeroPosition,
                    p if p < 0 => ReferenceStatus::NegativePosition,
                    _ => ReferenceStatus::MissingChromosome,
                },
            })
            .collect();
        let mut order: Vec<u32> = (0..bim.rows.len() as u32).collect();
        order.sort_unstable_by_key(|&i| {
            let row = bim.rows[i as usize];
            (row.contig, row.position)
        });
        let mut groups = vec![0..0; bim.contigs.len()];
        let mut start = 0;
        for chunk in
            order.chunk_by(|&x, &y| bim.rows[x as usize].contig == bim.rows[y as usize].contig)
        {
            groups[bim.rows[chunk[0] as usize].contig as usize] = start..start + chunk.len();
            start += chunk.len();
        }
        Self {
            bim,
            aliases,
            order,
            groups,
            catalog: HashMap::new(),
            lookups,
            names: vec![None; bim.contigs.len()],
            overlap: 0,
        }
    }

    /// Register a FASTA sequence and return the slots of `order` on it that
    /// have positive positions.
    fn sequence(&mut self, name: &str) -> Result<Range<usize>> {
        ensure!(!name.is_empty(), "FASTA sequence name is empty");
        let key = self.aliases.fasta_name(name);
        if let Some(previous) = self.catalog.insert(key.clone(), name.to_owned()) {
            anyhow::bail!(
                "FASTA chromosome alias collision: {previous:?} and {name:?} both map to {key:?}; use an unambiguous reference or different chromosome preset"
            );
        }
        let Some(&contig) = self.bim.contig_ids.get(&key) else {
            return Ok(0..0);
        };
        self.overlap += 1;
        self.names[contig as usize] = Some(name.to_owned());
        let range = self.groups[contig as usize].clone();
        let first = range.start
            + self.order[range.clone()]
                .partition_point(|&i| self.bim.rows[i as usize].position <= 0);
        for &i in &self.order[first..range.end] {
            self.lookups[i as usize].status = ReferenceStatus::OutOfRange;
        }
        Ok(first..range.end)
    }

    fn position(&self, slot: usize) -> u64 {
        self.bim.rows[self.order[slot] as usize].position as u64
    }

    fn found(&mut self, slot: usize, base: u8) {
        self.lookups[self.order[slot] as usize] = Lookup {
            base: Some(base.to_ascii_uppercase()),
            status: ReferenceStatus::Found,
        };
    }
}

/// Reject anything outside IUPAC DNA: such a reference is usually the wrong file.
/// `first` is the 1-based position of `bases[0]`.
fn validate_bases(bases: &[u8], name: &str, first: u64) -> Result<()> {
    if let Some(i) = bases.iter().position(|&base| complement(base).is_none()) {
        anyhow::bail!(
            "invalid FASTA base '{}' at {name}:{}; expected IUPAC DNA without embedded whitespace",
            bases[i].escape_ascii(),
            first + i as u64
        );
    }
    Ok(())
}

pub fn read(
    bim: &Index,
    path: &Path,
    index_path: Option<&Path>,
    mode: ReferenceMode,
    aliases: &Aliases,
) -> Result<ReferenceData> {
    // Pipes and stdin can only be streamed, and must be opened exactly once.
    let regular = is_regular(path);
    let default_index = sidecar(path, ".fai");
    let index_path = index_path.unwrap_or(&default_index);
    let indexed = match mode {
        ReferenceMode::Indexed => {
            ensure!(
                regular,
                "{}: indexed reference access needs a regular file; use --reference-mode stream",
                path.display()
            );
            true
        }
        ReferenceMode::Stream => false,
        ReferenceMode::Auto => {
            regular && index_path.is_file() && (!is_gzip(path)? || sidecar(path, ".gzi").is_file())
        }
    };
    let compressed = indexed && is_gzip(path)?;
    let mut scan = Scan::new(bim, aliases);
    if indexed {
        read_indexed(&mut scan, path, index_path, compressed)?;
    } else {
        read_stream(&mut scan, path)?;
    }
    ensure!(
        !scan.catalog.is_empty(),
        "{}: FASTA contains no sequences",
        path.display()
    );
    ensure!(
        scan.overlap > 0,
        "no chromosome overlap between BIM and FASTA; check chromosome names, --chromosome-set, and --chromosome-map"
    );
    Ok(ReferenceData {
        lookups: scan.lookups,
        names: scan.names,
        access: if indexed { "indexed" } else { "stream" },
    })
}

fn read_stream(scan: &mut Scan, path: &Path) -> Result<()> {
    let mut reader = fasta::io::Reader::new(open_text(path)?);
    let mut definition = fasta::record::Definition::default();
    let mut buffer = [0; 65536];
    while reader
        .read_definition(&mut definition)
        .with_context(|| format!("reading FASTA definition in {}", path.display()))?
        != 0
    {
        let name =
            std::str::from_utf8(definition.name().as_ref()).context("FASTA names must be UTF-8")?;
        let range = scan.sequence(name)?;
        let mut cursor = range.start;
        let mut offset = 0u64;
        let mut sequence = reader.sequence_reader();
        loop {
            let n = sequence
                .read(&mut buffer)
                .with_context(|| format!("reading sequence {name:?}"))?;
            if n == 0 {
                break;
            }
            let bases = &buffer[..n];
            validate_bases(bases, name, offset + 1)?;
            let end = offset + n as u64;
            while cursor < range.end {
                let position = scan.position(cursor);
                if position > end {
                    break;
                }
                scan.found(cursor, bases[(position - offset - 1) as usize]);
                cursor += 1;
            }
            offset = end;
        }
    }
    Ok(())
}

fn read_indexed(scan: &mut Scan, path: &Path, index_path: &Path, compressed: bool) -> Result<()> {
    let index = fasta::fai::fs::read(index_path).with_context(|| {
        format!(
            "reading FASTA index {}; use --reference-mode stream without an index",
            index_path.display()
        )
    })?;
    // Validate the complete catalog before querying: index lookup must never choose
    // an arbitrary record when reference names or normalized aliases collide.
    let mut sequences = Vec::new();
    for entry in index.as_ref() {
        ensure!(
            entry.line_width() >= entry.line_base_count(),
            "{}: invalid FASTA index line width",
            index_path.display()
        );
        let last = entry.length().saturating_sub(1);
        let bases = entry.line_base_count().get();
        ensure!(
            (last / bases)
                .checked_mul(entry.line_width().get())
                .and_then(|offset| offset.checked_add(last % bases))
                .and_then(|offset| offset.checked_add(entry.position()))
                .is_some(),
            "{}: FASTA index byte offset overflow",
            index_path.display()
        );
        let name = std::str::from_utf8(entry.name().as_ref())
            .context("FASTA names must be UTF-8")?
            .to_owned();
        let range = scan.sequence(&name)?;
        if !range.is_empty() {
            sequences.push((name, range, entry.length()));
        }
    }
    let inner = if compressed {
        fasta::io::BufReader::Bgzf(noodles_bgzf::io::indexed_reader::Builder::default().build_from_path(path)
            .context("opening indexed BGZF FASTA; .gzi is required, ordinary gzip needs --reference-mode stream")?)
    } else {
        fasta::io::BufReader::Uncompressed(BufReader::new(File::open(path)?))
    };
    let mut reader = fasta::io::IndexedReader::new(inner, index);
    for (name, range, length) in sequences {
        let mut cursor = range.start;
        while cursor < range.end {
            let start = scan.position(cursor);
            if start > length {
                break;
            }
            let mut next = cursor + 1;
            let mut end = start;
            // Coalesce dense positions, but never allocate an entire chromosome.
            while next < range.end {
                let position = scan.position(next);
                if position > length || position - start >= 65536 || position - end > 4096 {
                    break;
                }
                end = position;
                next += 1;
            }
            let region = Region::new(
                name.as_str(),
                Position::try_from(usize::try_from(start)?)?
                    ..=Position::try_from(usize::try_from(end)?)?,
            );
            let record = reader.query(&region).with_context(|| {
                format!("querying {region}; check that the FASTA index matches the reference")
            })?;
            let bases: &[u8] = record.sequence().as_ref();
            ensure!(
                bases.len() as u64 == end - start + 1,
                "short FASTA query for {region}; index may be stale"
            );
            validate_bases(bases, &name, start)?;
            for slot in cursor..next {
                scan.found(slot, bases[(scan.position(slot) - start) as usize]);
            }
            cursor = next;
        }
    }
    Ok(())
}
