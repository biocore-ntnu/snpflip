use std::{
    fs::{self, File},
    io::{BufRead, BufReader, Read},
    path::Path,
};

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;

pub fn is_stdin(path: &Path) -> bool {
    path.as_os_str() == "-"
}

/// Only regular files can be reopened or seeked; pipes and stdin cannot.
pub fn is_regular(path: &Path) -> bool {
    !is_stdin(path) && fs::metadata(path).is_ok_and(|m| m.is_file())
}

pub fn is_gzip(path: &Path) -> Result<bool> {
    let mut reader =
        BufReader::new(File::open(path).with_context(|| format!("opening {}", path.display()))?);
    Ok(reader.fill_buf()?.starts_with(&[0x1f, 0x8b]))
}

/// Open a path, or standard input for `-`. Compression is detected by content,
/// including concatenated gzip and BGZF members.
pub fn open_text(path: &Path) -> Result<Box<dyn BufRead>> {
    let raw: Box<dyn Read> = if is_stdin(path) {
        Box::new(std::io::stdin())
    } else {
        Box::new(File::open(path).with_context(|| format!("opening {}", path.display()))?)
    };
    let mut reader = BufReader::new(raw);
    if reader.fill_buf()?.starts_with(&[0x1f, 0x8b]) {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(reader))))
    } else {
        Ok(Box::new(reader))
    }
}
