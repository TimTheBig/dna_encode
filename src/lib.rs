use std::{fmt::{Debug, Display}, str::FromStr};

/// Create a new `DNACode`
/// ```
/// # use dna_encode::dna;
/// # assert_eq!(
/// dna![A, T, G],
/// // or
/// dna![
///     A, T, G;
///     T, A, C;
/// ],
/// # );
/// ```
#[macro_export]
macro_rules! dna {
    ($($nu:ident),+ $(,)?) => {
        $crate::DNACode::from((
            vec![$( $crate::dna!(stringify!($nu).as_bytes()[0] as char) ),+],
            // find complementary in const
            vec![
                $( $crate::dna!(find complementary,
                    ($crate::dna!(stringify!($nu).as_bytes()[0] as char))
                ) ),+
            ]
        ))
    };
    ($($nu1:ident),+ $(,)?; $($nu2:ident),+ $(,)? $(;)?) => {{
        let dna = $crate::DNACode::from((
            vec![$( $crate::dna!(stringify!($nu1).as_bytes()[0] as char) ),+],
            vec![$( $crate::dna!(stringify!($nu2).as_bytes()[0] as char) ),+]
        ));
        // check complementary
        #[cfg(debug_assertions)]
        { assert!(dna.is_valid()); }

        dna
    }};
    ($nu:expr) => {
        const { $crate::Nucleotide::from_caps($nu) }
    };
    (find complementary, $nu:expr) => {
        const { $nu.complementary() }
    };
}

/// Create a new `RNACode`
/// ```
/// # use dna_encode::rna;
/// rna![A, U];
/// ```
#[macro_export]
macro_rules! rna {
    ($($nu:ident),+ $(,)?) => {
        $crate::RNACode::from(vec![$( $crate::rna!(stringify!($nu).as_bytes()[0] as char) ),+]) // so user can input U or T
    };
    ($nu:expr) => {
        const { $crate::Nucleotide::from_caps($nu) }
    };
}

#[derive(PartialEq, Eq)]
pub struct DNACode {
    top: RNACode,
    bottom: RNACode,
}

impl Display for DNACode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for nucleotide in &self.top.values {
            write!(f, "{}", nucleotide.to_char_dna())?
        }

        write!(f, "\n")?;

        for nucleotide in &self.bottom.values {
            write!(f, "{}", nucleotide.to_char_dna())?
        }

        Ok(())
    }
}

impl Debug for DNACode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DNACode")
            .field("top", &self.top.values)
            .field("bottom", &self.bottom.values)
            .finish()
    }
}

impl DNACode {
    #[must_use]
    pub fn serialize_slice(bytes: &[u8]) -> Self {
        RNACode::serialize_slice(bytes).into()
    }

    /// (Top, Bottom)
    pub fn values(&self) -> (&[Nucleotide], &[Nucleotide]) {
        (&self.top.values, &self.bottom.values)
    }

    /// If top and bottom have the same len and complementary nucleotides
    pub fn is_valid(&self) -> bool {
        if self.top.values.len() != self.bottom.values.len() {
            return false;
        }

        let valid_bottom = self.top.complementary();

        valid_bottom == self.bottom
    }
}

impl From<RNACode> for DNACode {
    fn from(rna: RNACode) -> Self {
        let bottom = rna.complementary();

        Self { top: rna, bottom }
    }
}

#[derive(Debug, thiserror::Error)]
pub enum DNACodeParseError {
    /// Top and Bottom chains have different lengths
    #[error("Top and Bottom chains have different lengths, {} != {}", .top_len, .bottom_len)]
    TopBottomLenDiff {
        top_len: usize,
        bottom_len: usize,
    },
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
            return Err(DNACodeParseError::TopBottomLenDiff {
                top_len: top_rna.len(),
                bottom_len: bottom_rna.len(),
            });
        }

        Ok(Self {
            top: top_rna.parse::<RNACode>()?,
            bottom: bottom_rna.parse::<RNACode>()?,
        })
    }
}

impl From<(Vec<Nucleotide>, Vec<Nucleotide>)> for DNACode {
    fn from((top_codes, bottom_codes): (Vec<Nucleotide>, Vec<Nucleotide>)) -> Self {
        let dna = Self { top: RNACode { values: top_codes }, bottom: RNACode { values: bottom_codes } };

        #[cfg(debug_assertions)]
        { dna.is_valid(); }

        dna
    }
}

impl From<[Vec<Nucleotide>; 2]> for DNACode {
    fn from([top_codes, bottom_codes]: [Vec<Nucleotide>; 2]) -> Self {
        let dna = Self { top: RNACode { values: top_codes }, bottom: RNACode { values: bottom_codes } };

        #[cfg(debug_assertions)]
        { dna.is_valid(); }

        dna
    }
}

#[derive(Debug, PartialEq, Eq)]
pub struct RNACode {
    /// Half the nucleotides
    values: Vec<Nucleotide>
}

impl Display for RNACode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for nucleotide in &self.values {
            write!(f, "{:R<}", nucleotide)?;//.to_char_rna())?
        }

        Ok(())
    }
}

impl RNACode {
    #[must_use]
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
                });
            }
        }

        Self { values: output }
    }

    pub fn values(&self) -> &[Nucleotide] {
        &self.values
    }

    pub fn complementary(&self) -> Self {
        let mut complementary = Vec::with_capacity(self.values.len());

        for nucleotide in &self.values {
            complementary.push(nucleotide.complementary());
        }

        Self { values: complementary }
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

impl From<Vec<Nucleotide>> for RNACode {
    fn from(codes: Vec<Nucleotide>) -> Self {
        Self { values: codes }
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

    pub const fn complementary(&self) -> Nucleotide {
        match self {
            Nucleotide::A => Nucleotide::T,
            Nucleotide::C => Nucleotide::G,
            Nucleotide::G => Nucleotide::C,
            Nucleotide::T => Nucleotide::A,
        }
    }

    pub const fn from_caps(c: char) -> Nucleotide {
        match c {
            'A' => Nucleotide::A,
            'C' => Nucleotide::C,
            'G' => Nucleotide::G,
            'T' => Nucleotide::T,
            'U' => Nucleotide::T,
            _ => panic!("input must be a valid Nucleotide varent")
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
    use super::{RNACode, DNACode, Nucleotide, rna, dna};

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
        let dna_len = dna.bottom.values.len() + dna.top.values.len();

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

    #[test]
    fn test_macros() {
        assert_eq!(rna![A, U].to_string(), "AU");

        assert_eq!(dna![A, T, C].to_string(), "ATC\nTAG");

        assert_eq!(dna![A, T; T, A;].to_string(), "AT\nTA");
        assert_eq!(dna![A, U; T, A].to_string(), "AT\nTA");
    }

    #[test]
    #[should_panic]
    fn test_dna_macro_is_complement() {
        let _ = dna![A, T, G;T, A, T];
    }
}
