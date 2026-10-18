use crate::reference::{Lookup, ReferenceStatus};

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum Category {
    Forward,
    Reverse,
    Palindromic,
    Mismatch,
    NonSnp,
    Missing,
    Invalid,
}

impl Category {
    pub const ALL: [Self; 7] = [
        Self::Forward,
        Self::Reverse,
        Self::Palindromic,
        Self::Mismatch,
        Self::NonSnp,
        Self::Missing,
        Self::Invalid,
    ];
    pub fn name(self) -> &'static str {
        match self {
            Self::Forward => "forward",
            Self::Reverse => "reverse",
            Self::Palindromic => "palindromic",
            Self::Mismatch => "mismatch",
            Self::NonSnp => "non_snp",
            Self::Missing => "missing",
            Self::Invalid => "invalid",
        }
    }
    pub fn resolved(self) -> bool {
        matches!(self, Self::Forward | Self::Reverse)
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum IdStatus {
    Unique,
    Duplicate,
    Missing,
}

impl IdStatus {
    pub fn name(self) -> &'static str {
        match self {
            Self::Unique => "unique",
            Self::Duplicate => "duplicate",
            Self::Missing => "missing",
        }
    }
}

#[derive(Debug)]
pub struct Annotation {
    pub category: Category,
    pub reasons: Vec<&'static str>,
    pub palindromic: bool,
    pub id_status: IdStatus,
}

pub fn complement(base: u8) -> Option<u8> {
    Some(match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'N' => b'N',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        _ => return None,
    })
}

fn canonical(allele: &[u8]) -> Option<u8> {
    match allele {
        [base] if b"ACGT".contains(&base.to_ascii_uppercase()) => Some(base.to_ascii_uppercase()),
        _ => None,
    }
}

fn is_any(allele: &[u8], options: &[&[u8]]) -> bool {
    options.iter().any(|o| allele.eq_ignore_ascii_case(o))
}

/// `on_reference` says whether the row's chromosome exists in the FASTA.
pub fn classify(
    alleles: [&str; 2],
    lookup: &Lookup,
    on_reference: bool,
    id_status: IdStatus,
) -> Annotation {
    let [a, b] = alleles.map(str::as_bytes);
    let bases = canonical(a).zip(canonical(b));
    let palindromic = bases.is_some_and(|(a, b)| a != b && complement(a) == Some(b));
    let mut reasons = Vec::new();
    let allele_problem = if bases.is_some() {
        None
    } else if [a, b].iter().any(|s| {
        !is_any(s, &[b"0", b".", b"-", b"*", b"I", b"D"])
            && !s.iter().all(|&c| complement(c).is_some())
    }) {
        Some((Category::Invalid, "invalid_allele"))
    } else if [a, b]
        .iter()
        .any(|s| s.len() > 1 || is_any(s, &[b"I", b"D", b"-", b"*"]))
    {
        Some((Category::NonSnp, "unsupported_non_snp"))
    } else {
        Some((Category::Missing, "missing_or_ambiguous_allele"))
    };
    if let Some((_, reason)) = allele_problem {
        reasons.push(reason);
    }
    let reference_problem = match lookup.status {
        ReferenceStatus::Found => match lookup.base {
            Some(base) if b"ACGT".contains(&base) => None,
            Some(b'N') => Some("reference_n"),
            _ => Some("ambiguous_reference"),
        },
        ReferenceStatus::MissingChromosome => Some("missing_chromosome"),
        ReferenceStatus::ZeroPosition => Some("zero_position"),
        ReferenceStatus::NegativePosition => Some("negative_position"),
        ReferenceStatus::OutOfRange => Some("out_of_range"),
    };
    if let Some(reason) = reference_problem {
        reasons.push(reason);
    }
    if !on_reference
        && matches!(
            lookup.status,
            ReferenceStatus::ZeroPosition | ReferenceStatus::NegativePosition
        )
    {
        reasons.push("missing_chromosome");
    }
    if palindromic {
        reasons.push("palindromic");
    }
    let category = if let Some((category, _)) = allele_problem {
        category
    } else if reference_problem.is_some() {
        Category::Missing
    } else {
        let reference = lookup.base.expect("reference checked");
        let (a, b) = bases.expect("alleles checked");
        let forward = [a, b].contains(&reference);
        let reverse = [a, b].contains(&complement(reference).expect("canonical reference"));
        match (forward, reverse) {
            (true, false) => Category::Forward,
            (false, true) => Category::Reverse,
            (true, true) => Category::Palindromic,
            (false, false) => {
                reasons.push("reference_mismatch");
                Category::Mismatch
            }
        }
    };
    Annotation {
        category,
        reasons,
        palindromic,
        id_status,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn category(reference: u8, a: u8, b: u8) -> Category {
        let lookup = Lookup {
            base: Some(reference),
            status: ReferenceStatus::Found,
        };
        let [a, b] = [a, b].map(|c| char::from(c).to_string());
        classify([&a, &b], &lookup, true, IdStatus::Unique).category
    }
    #[test]
    fn all_canonical_combinations_obey_orientation_and_allele_order() {
        for &r in b"ACGT" {
            for &a in b"ACGT" {
                for &b in b"ACGT" {
                    let result = category(r, a, b);
                    let forward = a == r || b == r;
                    let reverse = Some(a) == complement(r) || Some(b) == complement(r);
                    let expected = match (forward, reverse) {
                        (true, false) => Category::Forward,
                        (false, true) => Category::Reverse,
                        (true, true) => Category::Palindromic,
                        _ => Category::Mismatch,
                    };
                    assert_eq!(
                        result, expected,
                        "ref={} alleles={}/{}",
                        r as char, a as char, b as char
                    );
                    assert_eq!(result, category(r, b, a));
                    let complemented = category(r, complement(a).unwrap(), complement(b).unwrap());
                    let expected_complement = match expected {
                        Category::Forward => Category::Reverse,
                        Category::Reverse => Category::Forward,
                        other => other,
                    };
                    assert_eq!(complemented, expected_complement);
                }
            }
        }
    }
    #[test]
    fn complements_are_an_involution_for_all_iupac_bases() {
        for &b in b"ACGTRYSWKMBDHVN" {
            assert_eq!(complement(complement(b).unwrap()), Some(b));
        }
    }
}
