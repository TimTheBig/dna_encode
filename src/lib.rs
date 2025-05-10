use std::{fmt::Display, str::FromStr};

/// The nucleotides, top bottom order
pub struct DNACode {
    top: Vec<Nucleotide>,
    bottom: Vec<Nucleotide>
}

impl Display for DNACode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for nucleotide in &self.top {
            write!(f, "{}", nucleotide.to_char_dna())?
        }

        write!(f, "\n")?;

        for nucleotide in &self.bottom {
            write!(f, "{}", nucleotide.to_char_dna())?
        }

        Ok(())
    }
}

impl DNACode {
    // bytes could be from bincode, which will be an optional feature
    pub fn serialize_slice(bytes: &[u8]) -> Self {
        RNACode::serialize_slice(bytes).into()
    }

    /// (Top, Bottom)
    pub fn values(&self) -> (&[Nucleotide], &[Nucleotide]) {
        (&self.top, &self.bottom)
    }
}

impl From<RNACode> for DNACode {
    fn from(rna: RNACode) -> Self {
        let mut bottom = Vec::with_capacity(rna.values.len());

        for nucleotide in &rna.values {
            bottom.push(nucleotide.complementary());
        }

        Self { top: rna.values, bottom }
    }
}

#[derive(Debug, thiserror::Error)]
pub enum DNACodeParseError {
    #[error("Top and Bottom chains have different lengths")]
    TopBottomLenDiff,
    /// Missing newline separating top and bottom chains
    #[error("Missing newline separating top and bottom chains")]
    NoNewline,
    // InvalidChar at usize
    #[error(transparent)]
    InvalidChar(#[from] InvalidChar),
}

#[derive(Debug, thiserror::Error)]
pub struct InvalidChar(char);

impl Display for InvalidChar {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Invalid character '{}'", self.0)
    }
}

impl FromStr for DNACode {
    type Err = DNACodeParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (top_rna, bottom_rna) = s.split_once('\n')
            .ok_or(DNACodeParseError::NoNewline)?;
        if top_rna.len() != bottom_rna.len() {
            return Err(DNACodeParseError::TopBottomLenDiff);
        }

        Ok(Self {
            top: top_rna.parse::<RNACode>()?.values,
            bottom: bottom_rna.parse::<RNACode>()?.values,
        })
    }
}

pub struct RNACode {
    /// Half the nucleotides
    values: Vec<Nucleotide>
}

impl Display for RNACode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for nucleotide in &self.values {
            write!(f, "{:R<}", nucleotide)?//.to_char_rna())?
        }

        Ok(())
    }
}

impl RNACode {
    // bytes could be from bincode, which will be an optional feature
    pub fn serialize_slice(bytes: &[u8]) -> Self {
        // four nucleotide(two bits) per byte
        let mut output = Vec::with_capacity(bytes.len() * 4);

        for byte in bytes {
            // get all the bits from least to most significant
            // bit indices
            const MASK_0: u8 = 0b_1;
            const MASK_1: u8 = 0b_10;
            const MASK_2: u8 = 0b_100;
            const MASK_3: u8 = 0b_1000;
            const MASK_4: u8 = 0b_10000;
            const MASK_5: u8 = 0b_100000;
            const MASK_6: u8 = 0b_1000000;
            const MASK_7: u8 = 0b_10000000;

            let bit_pairs = [
                ((MASK_0 & byte) > 0, (MASK_1 & byte) > 0),
                ((MASK_2 & byte) > 0, (MASK_3 & byte) > 0),
                ((MASK_4 & byte) > 0, (MASK_5 & byte) > 0),
                ((MASK_6 & byte) > 0, (MASK_7 & byte) > 0),
            ];

            for pair in bit_pairs {
                output.push(match pair {
                    (false, false) => Nucleotide::A,
                    (false, true) => Nucleotide::C,
                    (true, false) => Nucleotide::G,
                    (true, true) => Nucleotide::T,
                })
            }
        }

        Self { values: output }
    }

    pub fn values(&self) -> &[Nucleotide] {
        &self.values
    }
}

impl FromStr for RNACode {
    type Err = InvalidChar;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut values = Vec::with_capacity(s.len());
        for char in s.chars() {
            values.push(char.try_into()?)
        }

        Ok(RNACode { values })
    }
}

/// A Nucleotide base of DNA or RNA, [wikipedia: Nucleotide base](https://en.wikipedia.org/wiki/Nucleotide_base)
#[repr(u8)]
#[derive(Debug, PartialEq, Eq)]
pub enum Nucleotide {
    /// Adenine = 00
    A = 0b_00,
    /// Cytosine = 01
    C = 0b_01,
    /// Guanine = 10
    G = 0b_10,
    /// Thymine = 11, (U) Uracil in RNA and will be displayed as such
    T = 0b_11,
}

impl Display for Nucleotide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if f.fill() == 'D' {
            match self {
                Nucleotide::A => write!(f, "A"),
                Nucleotide::C => write!(f, "C"),
                Nucleotide::G => write!(f, "G"),
                Nucleotide::T => write!(f, "T"),
            }
        } else if f.fill() == 'R' {
            match self {
                Nucleotide::A => write!(f, "A"),
                Nucleotide::C => write!(f, "C"),
                Nucleotide::G => write!(f, "G"),
                Nucleotide::T => write!(f, "U"),
            }
        } else {
            panic!("fill for Nucleotide Display must be 'D'(DNA) or 'R'(RNA) like `{{:D<}}` or `{{:R<}}`")
        }
    }
}

impl Nucleotide {
    pub const fn to_bin(&self) -> [bool; 2] {
        match self {
            Nucleotide::A => [false, false],
            Nucleotide::C => [false, true],
            Nucleotide::G => [true, false],
            Nucleotide::T => [true, true],
        }
    }

    #[inline]
    pub const fn to_char_dna(&self) -> char {
        match self {
            Nucleotide::A => 'A',
            Nucleotide::C => 'C',
            Nucleotide::G => 'G',
            Nucleotide::T => 'T',
        }
    }

    #[inline]
    pub const fn to_char_rna(&self) -> char {
        match self {
            Nucleotide::A => 'A',
            Nucleotide::C => 'C',
            Nucleotide::G => 'G',
            Nucleotide::T => 'U',
        }
    }

    fn complementary(&self) -> Nucleotide {
        match self {
            Nucleotide::A => Nucleotide::T,
            Nucleotide::C => Nucleotide::G,
            Nucleotide::G => Nucleotide::C,
            Nucleotide::T => Nucleotide::A,
        }
    }
}

impl TryFrom<char> for Nucleotide {
    type Error = InvalidChar;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c {
            'A' | 'a' => Ok(Nucleotide::A),
            'C' | 'c' => Ok(Nucleotide::C),
            'G' | 'g' => Ok(Nucleotide::G),
            'T' | 't' => Ok(Nucleotide::T),
            'U' | 'u' => Ok(Nucleotide::T),
            _ => Err(InvalidChar(c))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{RNACode, DNACode, Nucleotide};

    #[test]
    fn test_rna() {
        let rna: RNACode = "GCGCUGUCGA".parse().unwrap();
        assert_eq!(rna.to_string(), "GCGCUGUCGA");

        eprintln!("{}", RNACode::serialize_slice(include_bytes!("./lib.rs")));
    }

    #[test]
    fn test_dna() {
        let dna: DNACode = "GCGCTGTCGA\nCGCGACGGCT".parse().unwrap();
        assert_eq!(dna.to_string(), "GCGCTGTCGA\nCGCGACGGCT");

        eprintln!("{}", DNACode::serialize_slice(include_bytes!("./lib.rs")));
    }

    #[test]
    fn test_dna_twice_rna() {
        let dna = DNACode::serialize_slice(include_bytes!("./lib.rs"));
        let dna_len = dna.bottom.len() + dna.top.len();

        assert_eq!(dna_len / 2, RNACode::serialize_slice(include_bytes!("./lib.rs")).values.len());
    }

    #[test]
    fn test_nucleotide() {
        assert_eq!(format!("{:R<}", Nucleotide::T), "U");
        assert_eq!(format!("{:D<}", Nucleotide::T), "T");

        assert_eq!(Nucleotide::try_from('A').unwrap(), Nucleotide::try_from('a').unwrap());
    }

    #[test]
    #[should_panic]
    fn test_nucleotide_format() {
        // T/U is different for rna and dna so this must be specified with fmt fill
        let _ = format!("{}", Nucleotide::T);
    }
}
