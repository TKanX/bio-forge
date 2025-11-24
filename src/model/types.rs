//! Atomic and residue classification primitives shared across the `bio-forge` model layer.
//!
//! The module centralizes periodic-table metadata, residue naming conventions, and helper
//! utilities that are reused by topology builders, IO translators, and downstream
//! operations such as solvation and hydrogen placement.

use nalgebra::Point3;
use std::fmt;
use std::str::FromStr;

/// Cartesian coordinate storing atomic positions in ångströms.
///
/// The alias ensures geometric calculations share the exact numeric representation across
/// IO readers, topology builders, and structure editing utilities.
pub type Point = Point3<f64>;

/// Periodic table entries supported by `bio-forge`.
///
/// Each variant encodes its atomic number in the `repr(u8)` discriminant and provides
/// conveniences for retrieving standard symbols, heavy-atom classification, and atomic
/// masses used when exporting coordinates or computing physical properties.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Deserialize)]
#[repr(u8)]
pub enum Element {
    /// Hydrogen (Z = 1).
    H = 1,
    /// Helium (Z = 2).
    He = 2,
    /// Lithium (Z = 3).
    Li = 3,
    /// Beryllium (Z = 4).
    Be = 4,
    /// Boron (Z = 5).
    B = 5,
    /// Carbon (Z = 6).
    C = 6,
    /// Nitrogen (Z = 7).
    N = 7,
    /// Oxygen (Z = 8).
    O = 8,
    /// Fluorine (Z = 9).
    F = 9,
    /// Neon (Z = 10).
    Ne = 10,
    /// Sodium (Z = 11).
    Na = 11,
    /// Magnesium (Z = 12).
    Mg = 12,
    /// Aluminum (Z = 13).
    Al = 13,
    /// Silicon (Z = 14).
    Si = 14,
    /// Phosphorus (Z = 15).
    P = 15,
    /// Sulfur (Z = 16).
    S = 16,
    /// Chlorine (Z = 17).
    Cl = 17,
    /// Argon (Z = 18).
    Ar = 18,
    /// Potassium (Z = 19).
    K = 19,
    /// Calcium (Z = 20).
    Ca = 20,
    /// Scandium (Z = 21).
    Sc = 21,
    /// Titanium (Z = 22).
    Ti = 22,
    /// Vanadium (Z = 23).
    V = 23,
    /// Chromium (Z = 24).
    Cr = 24,
    /// Manganese (Z = 25).
    Mn = 25,
    /// Iron (Z = 26).
    Fe = 26,
    /// Cobalt (Z = 27).
    Co = 27,
    /// Nickel (Z = 28).
    Ni = 28,
    /// Copper (Z = 29).
    Cu = 29,
    /// Zinc (Z = 30).
    Zn = 30,
    /// Gallium (Z = 31).
    Ga = 31,
    /// Germanium (Z = 32).
    Ge = 32,
    /// Arsenic (Z = 33).
    As = 33,
    /// Selenium (Z = 34).
    Se = 34,
    /// Bromine (Z = 35).
    Br = 35,
    /// Krypton (Z = 36).
    Kr = 36,
    /// Rubidium (Z = 37).
    Rb = 37,
    /// Strontium (Z = 38).
    Sr = 38,
    /// Yttrium (Z = 39).
    Y = 39,
    /// Zirconium (Z = 40).
    Zr = 40,
    /// Niobium (Z = 41).
    Nb = 41,
    /// Molybdenum (Z = 42).
    Mo = 42,
    /// Technetium (Z = 43).
    Tc = 43,
    /// Ruthenium (Z = 44).
    Ru = 44,
    /// Rhodium (Z = 45).
    Rh = 45,
    /// Palladium (Z = 46).
    Pd = 46,
    /// Silver (Z = 47).
    Ag = 47,
    /// Cadmium (Z = 48).
    Cd = 48,
    /// Indium (Z = 49).
    In = 49,
    /// Tin (Z = 50).
    Sn = 50,
    /// Antimony (Z = 51).
    Sb = 51,
    /// Tellurium (Z = 52).
    Te = 52,
    /// Iodine (Z = 53).
    I = 53,
    /// Xenon (Z = 54).
    Xe = 54,
    /// Cesium (Z = 55).
    Cs = 55,
    /// Barium (Z = 56).
    Ba = 56,
    /// Lanthanum (Z = 57).
    La = 57,
    /// Cerium (Z = 58).
    Ce = 58,
    /// Praseodymium (Z = 59).
    Pr = 59,
    /// Neodymium (Z = 60).
    Nd = 60,
    /// Promethium (Z = 61).
    Pm = 61,
    /// Samarium (Z = 62).
    Sm = 62,
    /// Europium (Z = 63).
    Eu = 63,
    /// Gadolinium (Z = 64).
    Gd = 64,
    /// Terbium (Z = 65).
    Tb = 65,
    /// Dysprosium (Z = 66).
    Dy = 66,
    /// Holmium (Z = 67).
    Ho = 67,
    /// Erbium (Z = 68).
    Er = 68,
    /// Thulium (Z = 69).
    Tm = 69,
    /// Ytterbium (Z = 70).
    Yb = 70,
    /// Lutetium (Z = 71).
    Lu = 71,
    /// Hafnium (Z = 72).
    Hf = 72,
    /// Tantalum (Z = 73).
    Ta = 73,
    /// Tungsten (Z = 74).
    W = 74,
    /// Rhenium (Z = 75).
    Re = 75,
    /// Osmium (Z = 76).
    Os = 76,
    /// Iridium (Z = 77).
    Ir = 77,
    /// Platinum (Z = 78).
    Pt = 78,
    /// Gold (Z = 79).
    Au = 79,
    /// Mercury (Z = 80).
    Hg = 80,
    /// Thallium (Z = 81).
    Tl = 81,
    /// Lead (Z = 82).
    Pb = 82,
    /// Bismuth (Z = 83).
    Bi = 83,
    /// Polonium (Z = 84).
    Po = 84,
    /// Astatine (Z = 85).
    At = 85,
    /// Radon (Z = 86).
    Rn = 86,
    /// Francium (Z = 87).
    Fr = 87,
    /// Radium (Z = 88).
    Ra = 88,
    /// Actinium (Z = 89).
    Ac = 89,
    /// Thorium (Z = 90).
    Th = 90,
    /// Protactinium (Z = 91).
    Pa = 91,
    /// Uranium (Z = 92).
    U = 92,
    /// Neptunium (Z = 93).
    Np = 93,
    /// Plutonium (Z = 94).
    Pu = 94,
    /// Americium (Z = 95).
    Am = 95,
    /// Curium (Z = 96).
    Cm = 96,
    /// Berkelium (Z = 97).
    Bk = 97,
    /// Californium (Z = 98).
    Cf = 98,
    /// Einsteinium (Z = 99).
    Es = 99,
    /// Fermium (Z = 100).
    Fm = 100,
    /// Mendelevium (Z = 101).
    Md = 101,
    /// Nobelium (Z = 102).
    No = 102,
    /// Lawrencium (Z = 103).
    Lr = 103,
    /// Rutherfordium (Z = 104).
    Rf = 104,
    /// Dubnium (Z = 105).
    Db = 105,
    /// Seaborgium (Z = 106).
    Sg = 106,
    /// Bohrium (Z = 107).
    Bh = 107,
    /// Hassium (Z = 108).
    Hs = 108,
    /// Meitnerium (Z = 109).
    Mt = 109,
    /// Darmstadtium (Z = 110).
    Ds = 110,
    /// Roentgenium (Z = 111).
    Rg = 111,
    /// Copernicium (Z = 112).
    Cn = 112,
    /// Nihonium (Z = 113).
    Nh = 113,
    /// Flerovium (Z = 114).
    Fl = 114,
    /// Moscovium (Z = 115).
    Mc = 115,
    /// Livermorium (Z = 116).
    Lv = 116,
    /// Tennessine (Z = 117).
    Ts = 117,
    /// Oganesson (Z = 118).
    Og = 118,
    /// Placeholder for unsupported or user-defined elements.
    Unknown = 0,
}

/// Canonical polymer residue names recognized by `bio-forge`.
///
/// Variants cover the 20 standard amino acids, the common nucleic acid bases (both RNA and
/// DNA prefixed forms), inosine, and crystallographic water. Helper methods report whether
/// a residue belongs to protein or nucleic acid polymers so that operations can branch on
/// biopolymer type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Deserialize)]
pub enum StandardResidue {
    /// Alanine residue (ALA).
    ALA,
    /// Arginine residue (ARG).
    ARG,
    /// Asparagine residue (ASN).
    ASN,
    /// Aspartate residue (ASP).
    ASP,
    /// Cysteine residue (CYS).
    CYS,
    /// Glutamine residue (GLN).
    GLN,
    /// Glutamate residue (GLU).
    GLU,
    /// Glycine residue (GLY).
    GLY,
    /// Histidine residue (HIS).
    HIS,
    /// Isoleucine residue (ILE).
    ILE,
    /// Leucine residue (LEU).
    LEU,
    /// Lysine residue (LYS).
    LYS,
    /// Methionine residue (MET).
    MET,
    /// Phenylalanine residue (PHE).
    PHE,
    /// Proline residue (PRO).
    PRO,
    /// Serine residue (SER).
    SER,
    /// Threonine residue (THR).
    THR,
    /// Tryptophan residue (TRP).
    TRP,
    /// Tyrosine residue (TYR).
    TYR,
    /// Valine residue (VAL).
    VAL,
    /// Adenosine residue (A).
    A,
    /// Cytidine residue (C).
    C,
    /// Guanosine residue (G).
    G,
    /// Uridine residue (U).
    U,
    /// Inosine residue (I).
    I,
    /// Deoxyadenosine residue (DA).
    DA,
    /// Deoxycytidine residue (DC).
    DC,
    /// Deoxyguanosine residue (DG).
    DG,
    /// Deoxythymidine residue (DT).
    DT,
    /// Deoxyinosine residue (DI).
    DI,
    /// Water molecule (HOH).
    HOH,
}

/// Bond multiplicity definitions used when describing topology edges.
///
/// Distinguishes classical single, double, triple, and aromatic delocalized interactions so
/// that exported formats and force-field derivations can preserve chemical intent.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Deserialize)]
pub enum BondOrder {
    /// Classical single bond.
    Single,
    /// Double bond with two shared electron pairs.
    Double,
    /// Triple bond with three shared electron pairs.
    Triple,
    /// Aromatic bond representing delocalized pi systems.
    Aromatic,
}

/// High-level classification for residues appearing in a structure.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ResidueCategory {
    /// Standard polymer residue (protein or nucleic acid).
    Standard,
    /// Heterogen group such as ligands or modified residues.
    Hetero,
    /// Free ion such as `Na+` or `Cl-`.
    Ion,
}

/// Topological context of a residue within a polymer chain.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ResiduePosition {
    /// No specific positional annotation is available.
    None,
    /// Residue resides in the middle of a polymer segment.
    Internal,
    /// Residue occupies the N-terminus of a polypeptide.
    NTerminal,
    /// Residue occupies the C-terminus of a polypeptide.
    CTerminal,
    /// Residue lies at the 5' terminus of a nucleic acid.
    FivePrime,
    /// Residue lies at the 3' terminus of a nucleic acid.
    ThreePrime,
}

/// Checks whether the residue belongs to the 20 canonical amino acids.
///
/// Used by cleaning and repair operations to decide when to apply protein-specific
/// heuristics (rotamer fitting, backbone rules, etc.).
///
/// # Returns
///
/// `true` if the residue is part of the protein alphabet.
impl StandardResidue {
    pub fn is_protein(self) -> bool {
        matches!(
            self,
            StandardResidue::ALA
                | StandardResidue::ARG
                | StandardResidue::ASN
                | StandardResidue::ASP
                | StandardResidue::CYS
                | StandardResidue::GLN
                | StandardResidue::GLU
                | StandardResidue::GLY
                | StandardResidue::HIS
                | StandardResidue::ILE
                | StandardResidue::LEU
                | StandardResidue::LYS
                | StandardResidue::MET
                | StandardResidue::PHE
                | StandardResidue::PRO
                | StandardResidue::SER
                | StandardResidue::THR
                | StandardResidue::TRP
                | StandardResidue::TYR
                | StandardResidue::VAL
        )
    }

    /// Checks whether the residue belongs to the standard nucleic acid alphabet.
    ///
    /// Includes RNA (`A`, `C`, `G`, `U`, `I`) and DNA (`DA`, `DC`, `DG`, `DT`, `DI`) forms so
    /// that sugar-phosphate handling can branch appropriately.
    ///
    /// # Returns
    ///
    /// `true` if the residue should be treated as a nucleotide.
    pub fn is_nucleic(self) -> bool {
        matches!(
            self,
            StandardResidue::A
                | StandardResidue::C
                | StandardResidue::G
                | StandardResidue::U
                | StandardResidue::I
                | StandardResidue::DA
                | StandardResidue::DC
                | StandardResidue::DG
                | StandardResidue::DT
                | StandardResidue::DI
        )
    }
}

impl BondOrder {
    /// Converts the bond order to a numeric multiplicity.
    ///
    /// The returned value feeds file writers and geometry algorithms that require
    /// floating-point weights while still preserving the qualitative order.
    ///
    /// # Returns
    ///
    /// The multiplicity expressed as an `f64` (1.0, 2.0, 3.0, or 1.5 for aromatic).
    pub fn value(&self) -> f64 {
        match self {
            BondOrder::Single => 1.0,
            BondOrder::Double => 2.0,
            BondOrder::Triple => 3.0,
            BondOrder::Aromatic => 1.5,
        }
    }
}

impl fmt::Display for BondOrder {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.value())
    }
}

impl FromStr for BondOrder {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "1" | "1.0" | "Single" => Ok(BondOrder::Single),
            "2" | "2.0" | "Double" => Ok(BondOrder::Double),
            "3" | "3.0" | "Triple" => Ok(BondOrder::Triple),
            "1.5" | "Aromatic" => Ok(BondOrder::Aromatic),
            _ => Err(format!("Invalid bond order: {}", s)),
        }
    }
}

impl ResidueCategory {
    /// Returns the human-readable label for the category.
    ///
    /// Primarily used by reporting layers when generating summaries or warnings in IO logs.
    ///
    /// # Returns
    ///
    /// The canonical descriptive name of the category.
    pub fn name(&self) -> &'static str {
        match self {
            ResidueCategory::Standard => "Standard Residue",
            ResidueCategory::Hetero => "Hetero Residue",
            ResidueCategory::Ion => "Ion",
        }
    }
}

impl fmt::Display for ResidueCategory {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.name())
    }
}

impl FromStr for ResidueCategory {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "Standard" => Ok(ResidueCategory::Standard),
            "Hetero" => Ok(ResidueCategory::Hetero),
            "Ion" => Ok(ResidueCategory::Ion),
            _ => Err(format!("Invalid residue category: {}", s)),
        }
    }
}

impl ResiduePosition {
    /// Returns a descriptive name for the positional annotation.
    ///
    /// This helper keeps UI messaging and diagnostic logs consistent when reporting
    /// polymer endpoints.
    ///
    /// # Returns
    ///
    /// A static string describing the position.
    pub fn name(&self) -> &'static str {
        match self {
            ResiduePosition::None => "None",
            ResiduePosition::Internal => "Internal",
            ResiduePosition::NTerminal => "N-Terminal",
            ResiduePosition::CTerminal => "C-Terminal",
            ResiduePosition::FivePrime => "5'-Terminal",
            ResiduePosition::ThreePrime => "3'-Terminal",
        }
    }
}

impl fmt::Display for ResiduePosition {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.name())
    }
}

impl FromStr for ResiduePosition {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "None" => Ok(ResiduePosition::None),
            "Internal" => Ok(ResiduePosition::Internal),
            "NTerminal" => Ok(ResiduePosition::NTerminal),
            "CTerminal" => Ok(ResiduePosition::CTerminal),
            "FivePrime" => Ok(ResiduePosition::FivePrime),
            "ThreePrime" => Ok(ResiduePosition::ThreePrime),
            _ => Err(format!("Invalid residue position: {}", s)),
        }
    }
}

impl Element {
    /// Retrieves the IUPAC element symbol.
    ///
    /// This is used across file writers to ensure the element column stays normalized
    /// regardless of how the atom was originally labeled.
    ///
    /// # Returns
    ///
    /// A static string slice containing the two-character symbol (or `"Unknown"`).
    pub fn symbol(&self) -> &'static str {
        match self {
            Element::H => "H",
            Element::He => "He",
            Element::Li => "Li",
            Element::Be => "Be",
            Element::B => "B",
            Element::C => "C",
            Element::N => "N",
            Element::O => "O",
            Element::F => "F",
            Element::Ne => "Ne",
            Element::Na => "Na",
            Element::Mg => "Mg",
            Element::Al => "Al",
            Element::Si => "Si",
            Element::P => "P",
            Element::S => "S",
            Element::Cl => "Cl",
            Element::Ar => "Ar",
            Element::K => "K",
            Element::Ca => "Ca",
            Element::Sc => "Sc",
            Element::Ti => "Ti",
            Element::V => "V",
            Element::Cr => "Cr",
            Element::Mn => "Mn",
            Element::Fe => "Fe",
            Element::Co => "Co",
            Element::Ni => "Ni",
            Element::Cu => "Cu",
            Element::Zn => "Zn",
            Element::Ga => "Ga",
            Element::Ge => "Ge",
            Element::As => "As",
            Element::Se => "Se",
            Element::Br => "Br",
            Element::Kr => "Kr",
            Element::Rb => "Rb",
            Element::Sr => "Sr",
            Element::Y => "Y",
            Element::Zr => "Zr",
            Element::Nb => "Nb",
            Element::Mo => "Mo",
            Element::Tc => "Tc",
            Element::Ru => "Ru",
            Element::Rh => "Rh",
            Element::Pd => "Pd",
            Element::Ag => "Ag",
            Element::Cd => "Cd",
            Element::In => "In",
            Element::Sn => "Sn",
            Element::Sb => "Sb",
            Element::Te => "Te",
            Element::I => "I",
            Element::Xe => "Xe",
            Element::Cs => "Cs",
            Element::Ba => "Ba",
            Element::La => "La",
            Element::Ce => "Ce",
            Element::Pr => "Pr",
            Element::Nd => "Nd",
            Element::Pm => "Pm",
            Element::Sm => "Sm",
            Element::Eu => "Eu",
            Element::Gd => "Gd",
            Element::Tb => "Tb",
            Element::Dy => "Dy",
            Element::Ho => "Ho",
            Element::Er => "Er",
            Element::Tm => "Tm",
            Element::Yb => "Yb",
            Element::Lu => "Lu",
            Element::Hf => "Hf",
            Element::Ta => "Ta",
            Element::W => "W",
            Element::Re => "Re",
            Element::Os => "Os",
            Element::Ir => "Ir",
            Element::Pt => "Pt",
            Element::Au => "Au",
            Element::Hg => "Hg",
            Element::Tl => "Tl",
            Element::Pb => "Pb",
            Element::Bi => "Bi",
            Element::Po => "Po",
            Element::At => "At",
            Element::Rn => "Rn",
            Element::Fr => "Fr",
            Element::Ra => "Ra",
            Element::Ac => "Ac",
            Element::Th => "Th",
            Element::Pa => "Pa",
            Element::U => "U",
            Element::Np => "Np",
            Element::Pu => "Pu",
            Element::Am => "Am",
            Element::Cm => "Cm",
            Element::Bk => "Bk",
            Element::Cf => "Cf",
            Element::Es => "Es",
            Element::Fm => "Fm",
            Element::Md => "Md",
            Element::No => "No",
            Element::Lr => "Lr",
            Element::Rf => "Rf",
            Element::Db => "Db",
            Element::Sg => "Sg",
            Element::Bh => "Bh",
            Element::Hs => "Hs",
            Element::Mt => "Mt",
            Element::Ds => "Ds",
            Element::Rg => "Rg",
            Element::Cn => "Cn",
            Element::Nh => "Nh",
            Element::Fl => "Fl",
            Element::Mc => "Mc",
            Element::Lv => "Lv",
            Element::Ts => "Ts",
            Element::Og => "Og",
            Element::Unknown => "Unknown",
        }
    }

    /// Indicates whether the element is treated as a heavy atom.
    ///
    /// Hydrogen is the only light atom; every other element (including `Unknown`) is flagged
    /// as heavy for the purpose of filtering coordinate exports.
    ///
    /// # Returns
    ///
    /// `true` if the element is not hydrogen.
    pub fn is_heavy_atom(&self) -> bool {
        !matches!(self, Element::H)
    }

    /// Returns the standard atomic mass in unified atomic mass units.
    ///
    /// Values follow IUPAC tables and feed mass-based calculations such as center-of-mass or
    /// density estimations in solvation operations.
    ///
    /// # Returns
    ///
    /// The atomic mass expressed as `f64`; `0.0` for `Unknown`.
    pub fn atomic_mass(&self) -> f64 {
        match self {
            Element::H => 1.00794,
            Element::He => 4.002602,
            Element::Li => 6.941,
            Element::Be => 9.012182,
            Element::B => 10.811,
            Element::C => 12.0107,
            Element::N => 14.0067,
            Element::O => 15.9994,
            Element::F => 18.9984032,
            Element::Ne => 20.1797,
            Element::Na => 22.98976928,
            Element::Mg => 24.3050,
            Element::Al => 26.9815386,
            Element::Si => 28.0855,
            Element::P => 30.973762,
            Element::S => 32.065,
            Element::Cl => 35.453,
            Element::Ar => 39.948,
            Element::K => 39.0983,
            Element::Ca => 40.078,
            Element::Sc => 44.955912,
            Element::Ti => 47.867,
            Element::V => 50.9415,
            Element::Cr => 51.9961,
            Element::Mn => 54.938045,
            Element::Fe => 55.845,
            Element::Co => 58.933195,
            Element::Ni => 58.6934,
            Element::Cu => 63.546,
            Element::Zn => 65.38,
            Element::Ga => 69.723,
            Element::Ge => 72.64,
            Element::As => 74.92160,
            Element::Se => 78.96,
            Element::Br => 79.904,
            Element::Kr => 83.798,
            Element::Rb => 85.4678,
            Element::Sr => 87.62,
            Element::Y => 88.90585,
            Element::Zr => 91.224,
            Element::Nb => 92.90638,
            Element::Mo => 95.96,
            Element::Tc => 98.0,
            Element::Ru => 101.07,
            Element::Rh => 102.90550,
            Element::Pd => 106.42,
            Element::Ag => 107.8682,
            Element::Cd => 112.411,
            Element::In => 114.818,
            Element::Sn => 118.710,
            Element::Sb => 121.760,
            Element::Te => 127.60,
            Element::I => 126.90447,
            Element::Xe => 131.293,
            Element::Cs => 132.9054519,
            Element::Ba => 137.327,
            Element::La => 138.90547,
            Element::Ce => 140.116,
            Element::Pr => 140.90765,
            Element::Nd => 144.242,
            Element::Pm => 145.0,
            Element::Sm => 150.36,
            Element::Eu => 151.964,
            Element::Gd => 157.25,
            Element::Tb => 158.92535,
            Element::Dy => 162.500,
            Element::Ho => 164.93032,
            Element::Er => 167.259,
            Element::Tm => 168.93421,
            Element::Yb => 173.054,
            Element::Lu => 174.9668,
            Element::Hf => 178.49,
            Element::Ta => 180.94788,
            Element::W => 183.84,
            Element::Re => 186.207,
            Element::Os => 190.23,
            Element::Ir => 192.217,
            Element::Pt => 195.084,
            Element::Au => 196.966569,
            Element::Hg => 200.59,
            Element::Tl => 204.3833,
            Element::Pb => 207.2,
            Element::Bi => 208.98040,
            Element::Po => 209.0,
            Element::At => 210.0,
            Element::Rn => 222.0,
            Element::Fr => 223.0,
            Element::Ra => 226.0,
            Element::Ac => 227.0,
            Element::Th => 232.03806,
            Element::Pa => 231.03588,
            Element::U => 238.02891,
            Element::Np => 237.0,
            Element::Pu => 244.0,
            Element::Am => 243.0,
            Element::Cm => 247.0,
            Element::Bk => 247.0,
            Element::Cf => 251.0,
            Element::Es => 252.0,
            Element::Fm => 257.0,
            Element::Md => 258.0,
            Element::No => 259.0,
            Element::Lr => 262.0,
            Element::Rf => 267.0,
            Element::Db => 268.0,
            Element::Sg => 271.0,
            Element::Bh => 272.0,
            Element::Hs => 270.0,
            Element::Mt => 276.0,
            Element::Ds => 281.0,
            Element::Rg => 280.0,
            Element::Cn => 285.0,
            Element::Nh => 284.0,
            Element::Fl => 289.0,
            Element::Mc => 288.0,
            Element::Lv => 293.0,
            Element::Ts => 294.0,
            Element::Og => 294.0,
            Element::Unknown => 0.0,
        }
    }
}

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.symbol())
    }
}

impl FromStr for Element {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Ok(num) = s.parse::<u8>() {
            match num {
                0 => Ok(Element::Unknown),
                1 => Ok(Element::H),
                2 => Ok(Element::He),
                3 => Ok(Element::Li),
                4 => Ok(Element::Be),
                5 => Ok(Element::B),
                6 => Ok(Element::C),
                7 => Ok(Element::N),
                8 => Ok(Element::O),
                9 => Ok(Element::F),
                10 => Ok(Element::Ne),
                11 => Ok(Element::Na),
                12 => Ok(Element::Mg),
                13 => Ok(Element::Al),
                14 => Ok(Element::Si),
                15 => Ok(Element::P),
                16 => Ok(Element::S),
                17 => Ok(Element::Cl),
                18 => Ok(Element::Ar),
                19 => Ok(Element::K),
                20 => Ok(Element::Ca),
                21 => Ok(Element::Sc),
                22 => Ok(Element::Ti),
                23 => Ok(Element::V),
                24 => Ok(Element::Cr),
                25 => Ok(Element::Mn),
                26 => Ok(Element::Fe),
                27 => Ok(Element::Co),
                28 => Ok(Element::Ni),
                29 => Ok(Element::Cu),
                30 => Ok(Element::Zn),
                31 => Ok(Element::Ga),
                32 => Ok(Element::Ge),
                33 => Ok(Element::As),
                34 => Ok(Element::Se),
                35 => Ok(Element::Br),
                36 => Ok(Element::Kr),
                37 => Ok(Element::Rb),
                38 => Ok(Element::Sr),
                39 => Ok(Element::Y),
                40 => Ok(Element::Zr),
                41 => Ok(Element::Nb),
                42 => Ok(Element::Mo),
                43 => Ok(Element::Tc),
                44 => Ok(Element::Ru),
                45 => Ok(Element::Rh),
                46 => Ok(Element::Pd),
                47 => Ok(Element::Ag),
                48 => Ok(Element::Cd),
                49 => Ok(Element::In),
                50 => Ok(Element::Sn),
                51 => Ok(Element::Sb),
                52 => Ok(Element::Te),
                53 => Ok(Element::I),
                54 => Ok(Element::Xe),
                55 => Ok(Element::Cs),
                56 => Ok(Element::Ba),
                57 => Ok(Element::La),
                58 => Ok(Element::Ce),
                59 => Ok(Element::Pr),
                60 => Ok(Element::Nd),
                61 => Ok(Element::Pm),
                62 => Ok(Element::Sm),
                63 => Ok(Element::Eu),
                64 => Ok(Element::Gd),
                65 => Ok(Element::Tb),
                66 => Ok(Element::Dy),
                67 => Ok(Element::Ho),
                68 => Ok(Element::Er),
                69 => Ok(Element::Tm),
                70 => Ok(Element::Yb),
                71 => Ok(Element::Lu),
                72 => Ok(Element::Hf),
                73 => Ok(Element::Ta),
                74 => Ok(Element::W),
                75 => Ok(Element::Re),
                76 => Ok(Element::Os),
                77 => Ok(Element::Ir),
                78 => Ok(Element::Pt),
                79 => Ok(Element::Au),
                80 => Ok(Element::Hg),
                81 => Ok(Element::Tl),
                82 => Ok(Element::Pb),
                83 => Ok(Element::Bi),
                84 => Ok(Element::Po),
                85 => Ok(Element::At),
                86 => Ok(Element::Rn),
                87 => Ok(Element::Fr),
                88 => Ok(Element::Ra),
                89 => Ok(Element::Ac),
                90 => Ok(Element::Th),
                91 => Ok(Element::Pa),
                92 => Ok(Element::U),
                93 => Ok(Element::Np),
                94 => Ok(Element::Pu),
                95 => Ok(Element::Am),
                96 => Ok(Element::Cm),
                97 => Ok(Element::Bk),
                98 => Ok(Element::Cf),
                99 => Ok(Element::Es),
                100 => Ok(Element::Fm),
                101 => Ok(Element::Md),
                102 => Ok(Element::No),
                103 => Ok(Element::Lr),
                104 => Ok(Element::Rf),
                105 => Ok(Element::Db),
                106 => Ok(Element::Sg),
                107 => Ok(Element::Bh),
                108 => Ok(Element::Hs),
                109 => Ok(Element::Mt),
                110 => Ok(Element::Ds),
                111 => Ok(Element::Rg),
                112 => Ok(Element::Cn),
                113 => Ok(Element::Nh),
                114 => Ok(Element::Fl),
                115 => Ok(Element::Mc),
                116 => Ok(Element::Lv),
                117 => Ok(Element::Ts),
                118 => Ok(Element::Og),
                _ => Ok(Element::Unknown),
            }
        } else {
            match s {
                "H" => Ok(Element::H),
                "He" => Ok(Element::He),
                "Li" => Ok(Element::Li),
                "Be" => Ok(Element::Be),
                "B" => Ok(Element::B),
                "C" => Ok(Element::C),
                "N" => Ok(Element::N),
                "O" => Ok(Element::O),
                "F" => Ok(Element::F),
                "Ne" => Ok(Element::Ne),
                "Na" => Ok(Element::Na),
                "Mg" => Ok(Element::Mg),
                "Al" => Ok(Element::Al),
                "Si" => Ok(Element::Si),
                "P" => Ok(Element::P),
                "S" => Ok(Element::S),
                "Cl" => Ok(Element::Cl),
                "Ar" => Ok(Element::Ar),
                "K" => Ok(Element::K),
                "Ca" => Ok(Element::Ca),
                "Sc" => Ok(Element::Sc),
                "Ti" => Ok(Element::Ti),
                "V" => Ok(Element::V),
                "Cr" => Ok(Element::Cr),
                "Mn" => Ok(Element::Mn),
                "Fe" => Ok(Element::Fe),
                "Co" => Ok(Element::Co),
                "Ni" => Ok(Element::Ni),
                "Cu" => Ok(Element::Cu),
                "Zn" => Ok(Element::Zn),
                "Ga" => Ok(Element::Ga),
                "Ge" => Ok(Element::Ge),
                "As" => Ok(Element::As),
                "Se" => Ok(Element::Se),
                "Br" => Ok(Element::Br),
                "Kr" => Ok(Element::Kr),
                "Rb" => Ok(Element::Rb),
                "Sr" => Ok(Element::Sr),
                "Y" => Ok(Element::Y),
                "Zr" => Ok(Element::Zr),
                "Nb" => Ok(Element::Nb),
                "Mo" => Ok(Element::Mo),
                "Tc" => Ok(Element::Tc),
                "Ru" => Ok(Element::Ru),
                "Rh" => Ok(Element::Rh),
                "Pd" => Ok(Element::Pd),
                "Ag" => Ok(Element::Ag),
                "Cd" => Ok(Element::Cd),
                "In" => Ok(Element::In),
                "Sn" => Ok(Element::Sn),
                "Sb" => Ok(Element::Sb),
                "Te" => Ok(Element::Te),
                "I" => Ok(Element::I),
                "Xe" => Ok(Element::Xe),
                "Cs" => Ok(Element::Cs),
                "Ba" => Ok(Element::Ba),
                "La" => Ok(Element::La),
                "Ce" => Ok(Element::Ce),
                "Pr" => Ok(Element::Pr),
                "Nd" => Ok(Element::Nd),
                "Pm" => Ok(Element::Pm),
                "Sm" => Ok(Element::Sm),
                "Eu" => Ok(Element::Eu),
                "Gd" => Ok(Element::Gd),
                "Tb" => Ok(Element::Tb),
                "Dy" => Ok(Element::Dy),
                "Ho" => Ok(Element::Ho),
                "Er" => Ok(Element::Er),
                "Tm" => Ok(Element::Tm),
                "Yb" => Ok(Element::Yb),
                "Lu" => Ok(Element::Lu),
                "Hf" => Ok(Element::Hf),
                "Ta" => Ok(Element::Ta),
                "W" => Ok(Element::W),
                "Re" => Ok(Element::Re),
                "Os" => Ok(Element::Os),
                "Ir" => Ok(Element::Ir),
                "Pt" => Ok(Element::Pt),
                "Au" => Ok(Element::Au),
                "Hg" => Ok(Element::Hg),
                "Tl" => Ok(Element::Tl),
                "Pb" => Ok(Element::Pb),
                "Bi" => Ok(Element::Bi),
                "Po" => Ok(Element::Po),
                "At" => Ok(Element::At),
                "Rn" => Ok(Element::Rn),
                "Fr" => Ok(Element::Fr),
                "Ra" => Ok(Element::Ra),
                "Ac" => Ok(Element::Ac),
                "Th" => Ok(Element::Th),
                "Pa" => Ok(Element::Pa),
                "U" => Ok(Element::U),
                "Np" => Ok(Element::Np),
                "Pu" => Ok(Element::Pu),
                "Am" => Ok(Element::Am),
                "Cm" => Ok(Element::Cm),
                "Bk" => Ok(Element::Bk),
                "Cf" => Ok(Element::Cf),
                "Es" => Ok(Element::Es),
                "Fm" => Ok(Element::Fm),
                "Md" => Ok(Element::Md),
                "No" => Ok(Element::No),
                "Lr" => Ok(Element::Lr),
                "Rf" => Ok(Element::Rf),
                "Db" => Ok(Element::Db),
                "Sg" => Ok(Element::Sg),
                "Bh" => Ok(Element::Bh),
                "Hs" => Ok(Element::Hs),
                "Mt" => Ok(Element::Mt),
                "Ds" => Ok(Element::Ds),
                "Rg" => Ok(Element::Rg),
                "Cn" => Ok(Element::Cn),
                "Nh" => Ok(Element::Nh),
                "Fl" => Ok(Element::Fl),
                "Mc" => Ok(Element::Mc),
                "Lv" => Ok(Element::Lv),
                "Ts" => Ok(Element::Ts),
                "Og" => Ok(Element::Og),
                "Unknown" => Ok(Element::Unknown),
                _ => Ok(Element::Unknown),
            }
        }
    }
}

impl fmt::Display for StandardResidue {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let code = match self {
            StandardResidue::ALA => "ALA",
            StandardResidue::ARG => "ARG",
            StandardResidue::ASN => "ASN",
            StandardResidue::ASP => "ASP",
            StandardResidue::CYS => "CYS",
            StandardResidue::GLN => "GLN",
            StandardResidue::GLU => "GLU",
            StandardResidue::GLY => "GLY",
            StandardResidue::HIS => "HIS",
            StandardResidue::ILE => "ILE",
            StandardResidue::LEU => "LEU",
            StandardResidue::LYS => "LYS",
            StandardResidue::MET => "MET",
            StandardResidue::PHE => "PHE",
            StandardResidue::PRO => "PRO",
            StandardResidue::SER => "SER",
            StandardResidue::THR => "THR",
            StandardResidue::TRP => "TRP",
            StandardResidue::TYR => "TYR",
            StandardResidue::VAL => "VAL",
            StandardResidue::A => "A",
            StandardResidue::C => "C",
            StandardResidue::G => "G",
            StandardResidue::U => "U",
            StandardResidue::I => "I",
            StandardResidue::DA => "DA",
            StandardResidue::DC => "DC",
            StandardResidue::DG => "DG",
            StandardResidue::DT => "DT",
            StandardResidue::DI => "DI",
            StandardResidue::HOH => "HOH",
        };
        write!(f, "{}", code)
    }
}

impl FromStr for StandardResidue {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "ALA" => Ok(StandardResidue::ALA),
            "ARG" => Ok(StandardResidue::ARG),
            "ASN" => Ok(StandardResidue::ASN),
            "ASP" => Ok(StandardResidue::ASP),
            "CYS" => Ok(StandardResidue::CYS),
            "GLN" => Ok(StandardResidue::GLN),
            "GLU" => Ok(StandardResidue::GLU),
            "GLY" => Ok(StandardResidue::GLY),
            "HIS" => Ok(StandardResidue::HIS),
            "ILE" => Ok(StandardResidue::ILE),
            "LEU" => Ok(StandardResidue::LEU),
            "LYS" => Ok(StandardResidue::LYS),
            "MET" => Ok(StandardResidue::MET),
            "PHE" => Ok(StandardResidue::PHE),
            "PRO" => Ok(StandardResidue::PRO),
            "SER" => Ok(StandardResidue::SER),
            "THR" => Ok(StandardResidue::THR),
            "TRP" => Ok(StandardResidue::TRP),
            "TYR" => Ok(StandardResidue::TYR),
            "VAL" => Ok(StandardResidue::VAL),
            "A" => Ok(StandardResidue::A),
            "C" => Ok(StandardResidue::C),
            "G" => Ok(StandardResidue::G),
            "U" => Ok(StandardResidue::U),
            "I" => Ok(StandardResidue::I),
            "DA" => Ok(StandardResidue::DA),
            "DC" => Ok(StandardResidue::DC),
            "DG" => Ok(StandardResidue::DG),
            "DT" => Ok(StandardResidue::DT),
            "DI" => Ok(StandardResidue::DI),
            "HOH" => Ok(StandardResidue::HOH),
            _ => Err(format!("Invalid standard residue: {}", s)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn element_symbol_returns_correct_value() {
        assert_eq!(Element::H.symbol(), "H");
        assert_eq!(Element::C.symbol(), "C");
        assert_eq!(Element::O.symbol(), "O");
        assert_eq!(Element::Fe.symbol(), "Fe");
        assert_eq!(Element::U.symbol(), "U");
        assert_eq!(Element::Unknown.symbol(), "Unknown");
    }

    #[test]
    fn element_is_heavy_atom_identifies_correctly() {
        assert!(!Element::H.is_heavy_atom());
        assert!(Element::C.is_heavy_atom());
        assert!(Element::O.is_heavy_atom());
        assert!(Element::Fe.is_heavy_atom());
        assert!(Element::Unknown.is_heavy_atom());
    }

    #[test]
    fn element_atomic_mass_returns_correct_value() {
        assert_eq!(Element::H.atomic_mass(), 1.00794);
        assert_eq!(Element::C.atomic_mass(), 12.0107);
        assert_eq!(Element::O.atomic_mass(), 15.9994);
        assert_eq!(Element::Fe.atomic_mass(), 55.845);
        assert_eq!(Element::Unknown.atomic_mass(), 0.0);
    }

    #[test]
    fn element_display_formats_correctly() {
        assert_eq!(format!("{}", Element::H), "H");
        assert_eq!(format!("{}", Element::C), "C");
        assert_eq!(format!("{}", Element::Unknown), "Unknown");
    }

    #[test]
    fn element_from_str_parses_valid_symbols() {
        assert_eq!(Element::from_str("H").unwrap(), Element::H);
        assert_eq!(Element::from_str("C").unwrap(), Element::C);
        assert_eq!(Element::from_str("O").unwrap(), Element::O);
        assert_eq!(Element::from_str("Fe").unwrap(), Element::Fe);
        assert_eq!(Element::from_str("Unknown").unwrap(), Element::Unknown);
    }

    #[test]
    fn element_from_str_parses_valid_numbers() {
        assert_eq!(Element::from_str("1").unwrap(), Element::H);
        assert_eq!(Element::from_str("6").unwrap(), Element::C);
        assert_eq!(Element::from_str("8").unwrap(), Element::O);
        assert_eq!(Element::from_str("26").unwrap(), Element::Fe);
        assert_eq!(Element::from_str("0").unwrap(), Element::Unknown);
        assert_eq!(Element::from_str("119").unwrap(), Element::Unknown);
    }

    #[test]
    fn element_from_str_handles_invalid_input() {
        assert_eq!(Element::from_str("Invalid").unwrap(), Element::Unknown);
        assert_eq!(Element::from_str("Zz").unwrap(), Element::Unknown);
        assert_eq!(Element::from_str("-1").unwrap(), Element::Unknown);
    }

    #[test]
    fn bond_order_value_returns_correct_f64() {
        assert_eq!(BondOrder::Single.value(), 1.0);
        assert_eq!(BondOrder::Double.value(), 2.0);
        assert_eq!(BondOrder::Triple.value(), 3.0);
        assert_eq!(BondOrder::Aromatic.value(), 1.5);
    }

    #[test]
    fn bond_order_display_formats_correctly() {
        assert_eq!(format!("{}", BondOrder::Single), "1");
        assert_eq!(format!("{}", BondOrder::Double), "2");
        assert_eq!(format!("{}", BondOrder::Triple), "3");
        assert_eq!(format!("{}", BondOrder::Aromatic), "1.5");
    }

    #[test]
    fn bond_order_from_str_parses_valid_inputs() {
        assert_eq!(BondOrder::from_str("1").unwrap(), BondOrder::Single);
        assert_eq!(BondOrder::from_str("1.0").unwrap(), BondOrder::Single);
        assert_eq!(BondOrder::from_str("Single").unwrap(), BondOrder::Single);
        assert_eq!(BondOrder::from_str("2").unwrap(), BondOrder::Double);
        assert_eq!(BondOrder::from_str("2.0").unwrap(), BondOrder::Double);
        assert_eq!(BondOrder::from_str("Double").unwrap(), BondOrder::Double);
        assert_eq!(BondOrder::from_str("3").unwrap(), BondOrder::Triple);
        assert_eq!(BondOrder::from_str("3.0").unwrap(), BondOrder::Triple);
        assert_eq!(BondOrder::from_str("Triple").unwrap(), BondOrder::Triple);
        assert_eq!(BondOrder::from_str("1.5").unwrap(), BondOrder::Aromatic);
        assert_eq!(
            BondOrder::from_str("Aromatic").unwrap(),
            BondOrder::Aromatic
        );
    }

    #[test]
    fn bond_order_from_str_handles_invalid_input() {
        assert!(BondOrder::from_str("4").is_err());
        assert!(BondOrder::from_str("Quadruple").is_err());
        assert!(BondOrder::from_str("1.2").is_err());
        assert!(BondOrder::from_str("").is_err());
    }

    #[test]
    fn residue_category_name_returns_correct_string() {
        assert_eq!(ResidueCategory::Standard.name(), "Standard Residue");
        assert_eq!(ResidueCategory::Hetero.name(), "Hetero Residue");
        assert_eq!(ResidueCategory::Ion.name(), "Ion");
    }

    #[test]
    fn residue_category_display_formats_correctly() {
        assert_eq!(format!("{}", ResidueCategory::Standard), "Standard Residue");
        assert_eq!(format!("{}", ResidueCategory::Hetero), "Hetero Residue");
        assert_eq!(format!("{}", ResidueCategory::Ion), "Ion");
    }

    #[test]
    fn residue_category_from_str_parses_valid_inputs() {
        assert_eq!(
            ResidueCategory::from_str("Standard").unwrap(),
            ResidueCategory::Standard
        );
        assert_eq!(
            ResidueCategory::from_str("Hetero").unwrap(),
            ResidueCategory::Hetero
        );
        assert_eq!(
            ResidueCategory::from_str("Ion").unwrap(),
            ResidueCategory::Ion
        );
    }

    #[test]
    fn residue_category_from_str_handles_invalid_input() {
        assert!(ResidueCategory::from_str("Unknown").is_err());
        assert!(ResidueCategory::from_str("").is_err());
        assert!(ResidueCategory::from_str("standard").is_err());
    }

    #[test]
    fn residue_position_name_returns_correct_string() {
        assert_eq!(ResiduePosition::None.name(), "None");
        assert_eq!(ResiduePosition::Internal.name(), "Internal");
        assert_eq!(ResiduePosition::NTerminal.name(), "N-Terminal");
        assert_eq!(ResiduePosition::CTerminal.name(), "C-Terminal");
        assert_eq!(ResiduePosition::FivePrime.name(), "5'-Terminal");
        assert_eq!(ResiduePosition::ThreePrime.name(), "3'-Terminal");
    }

    #[test]
    fn residue_position_display_formats_correctly() {
        assert_eq!(format!("{}", ResiduePosition::None), "None");
        assert_eq!(format!("{}", ResiduePosition::Internal), "Internal");
        assert_eq!(format!("{}", ResiduePosition::NTerminal), "N-Terminal");
        assert_eq!(format!("{}", ResiduePosition::CTerminal), "C-Terminal");
        assert_eq!(format!("{}", ResiduePosition::FivePrime), "5'-Terminal");
        assert_eq!(format!("{}", ResiduePosition::ThreePrime), "3'-Terminal");
    }

    #[test]
    fn residue_position_from_str_parses_valid_inputs() {
        assert_eq!(
            ResiduePosition::from_str("None").unwrap(),
            ResiduePosition::None
        );
        assert_eq!(
            ResiduePosition::from_str("Internal").unwrap(),
            ResiduePosition::Internal
        );
        assert_eq!(
            ResiduePosition::from_str("NTerminal").unwrap(),
            ResiduePosition::NTerminal
        );
        assert_eq!(
            ResiduePosition::from_str("CTerminal").unwrap(),
            ResiduePosition::CTerminal
        );
        assert_eq!(
            ResiduePosition::from_str("FivePrime").unwrap(),
            ResiduePosition::FivePrime
        );
        assert_eq!(
            ResiduePosition::from_str("ThreePrime").unwrap(),
            ResiduePosition::ThreePrime
        );
    }

    #[test]
    fn residue_position_from_str_handles_invalid_input() {
        assert!(ResiduePosition::from_str("Unknown").is_err());
        assert!(ResiduePosition::from_str("").is_err());
        assert!(ResiduePosition::from_str("internal").is_err());
    }

    #[test]
    fn standard_residue_display_formats_correctly() {
        assert_eq!(format!("{}", StandardResidue::ALA), "ALA");
        assert_eq!(format!("{}", StandardResidue::ARG), "ARG");
        assert_eq!(format!("{}", StandardResidue::ASN), "ASN");
        assert_eq!(format!("{}", StandardResidue::ASP), "ASP");
        assert_eq!(format!("{}", StandardResidue::CYS), "CYS");
        assert_eq!(format!("{}", StandardResidue::GLN), "GLN");
        assert_eq!(format!("{}", StandardResidue::GLU), "GLU");
        assert_eq!(format!("{}", StandardResidue::GLY), "GLY");
        assert_eq!(format!("{}", StandardResidue::HIS), "HIS");
        assert_eq!(format!("{}", StandardResidue::ILE), "ILE");
        assert_eq!(format!("{}", StandardResidue::LEU), "LEU");
        assert_eq!(format!("{}", StandardResidue::LYS), "LYS");
        assert_eq!(format!("{}", StandardResidue::MET), "MET");
        assert_eq!(format!("{}", StandardResidue::PHE), "PHE");
        assert_eq!(format!("{}", StandardResidue::PRO), "PRO");
        assert_eq!(format!("{}", StandardResidue::SER), "SER");
        assert_eq!(format!("{}", StandardResidue::THR), "THR");
        assert_eq!(format!("{}", StandardResidue::TRP), "TRP");
        assert_eq!(format!("{}", StandardResidue::TYR), "TYR");
        assert_eq!(format!("{}", StandardResidue::VAL), "VAL");
        assert_eq!(format!("{}", StandardResidue::A), "A");
        assert_eq!(format!("{}", StandardResidue::C), "C");
        assert_eq!(format!("{}", StandardResidue::G), "G");
        assert_eq!(format!("{}", StandardResidue::U), "U");
        assert_eq!(format!("{}", StandardResidue::I), "I");
        assert_eq!(format!("{}", StandardResidue::DA), "DA");
        assert_eq!(format!("{}", StandardResidue::DC), "DC");
        assert_eq!(format!("{}", StandardResidue::DG), "DG");
        assert_eq!(format!("{}", StandardResidue::DT), "DT");
        assert_eq!(format!("{}", StandardResidue::DI), "DI");
        assert_eq!(format!("{}", StandardResidue::HOH), "HOH");
    }

    #[test]
    fn standard_residue_from_str_parses_valid_inputs() {
        assert_eq!(
            StandardResidue::from_str("ALA").unwrap(),
            StandardResidue::ALA
        );
        assert_eq!(
            StandardResidue::from_str("ARG").unwrap(),
            StandardResidue::ARG
        );
        assert_eq!(
            StandardResidue::from_str("ASN").unwrap(),
            StandardResidue::ASN
        );
        assert_eq!(
            StandardResidue::from_str("ASP").unwrap(),
            StandardResidue::ASP
        );
        assert_eq!(
            StandardResidue::from_str("CYS").unwrap(),
            StandardResidue::CYS
        );
        assert_eq!(
            StandardResidue::from_str("GLN").unwrap(),
            StandardResidue::GLN
        );
        assert_eq!(
            StandardResidue::from_str("GLU").unwrap(),
            StandardResidue::GLU
        );
        assert_eq!(
            StandardResidue::from_str("GLY").unwrap(),
            StandardResidue::GLY
        );
        assert_eq!(
            StandardResidue::from_str("HIS").unwrap(),
            StandardResidue::HIS
        );
        assert_eq!(
            StandardResidue::from_str("ILE").unwrap(),
            StandardResidue::ILE
        );
        assert_eq!(
            StandardResidue::from_str("LEU").unwrap(),
            StandardResidue::LEU
        );
        assert_eq!(
            StandardResidue::from_str("LYS").unwrap(),
            StandardResidue::LYS
        );
        assert_eq!(
            StandardResidue::from_str("MET").unwrap(),
            StandardResidue::MET
        );
        assert_eq!(
            StandardResidue::from_str("PHE").unwrap(),
            StandardResidue::PHE
        );
        assert_eq!(
            StandardResidue::from_str("PRO").unwrap(),
            StandardResidue::PRO
        );
        assert_eq!(
            StandardResidue::from_str("SER").unwrap(),
            StandardResidue::SER
        );
        assert_eq!(
            StandardResidue::from_str("THR").unwrap(),
            StandardResidue::THR
        );
        assert_eq!(
            StandardResidue::from_str("TRP").unwrap(),
            StandardResidue::TRP
        );
        assert_eq!(
            StandardResidue::from_str("TYR").unwrap(),
            StandardResidue::TYR
        );
        assert_eq!(
            StandardResidue::from_str("VAL").unwrap(),
            StandardResidue::VAL
        );
        assert_eq!(StandardResidue::from_str("A").unwrap(), StandardResidue::A);
        assert_eq!(StandardResidue::from_str("C").unwrap(), StandardResidue::C);
        assert_eq!(StandardResidue::from_str("G").unwrap(), StandardResidue::G);
        assert_eq!(StandardResidue::from_str("U").unwrap(), StandardResidue::U);
        assert_eq!(StandardResidue::from_str("I").unwrap(), StandardResidue::I);
        assert_eq!(
            StandardResidue::from_str("DA").unwrap(),
            StandardResidue::DA
        );
        assert_eq!(
            StandardResidue::from_str("DC").unwrap(),
            StandardResidue::DC
        );
        assert_eq!(
            StandardResidue::from_str("DG").unwrap(),
            StandardResidue::DG
        );
        assert_eq!(
            StandardResidue::from_str("DT").unwrap(),
            StandardResidue::DT
        );
        assert_eq!(
            StandardResidue::from_str("DI").unwrap(),
            StandardResidue::DI
        );
        assert_eq!(
            StandardResidue::from_str("HOH").unwrap(),
            StandardResidue::HOH
        );
    }

    #[test]
    fn standard_residue_from_str_handles_invalid_input() {
        assert!(StandardResidue::from_str("XYZ").is_err());
        assert!(StandardResidue::from_str("Ala").is_err());
        assert!(StandardResidue::from_str("").is_err());
        assert!(StandardResidue::from_str("ALA ").is_err());
        assert!(StandardResidue::from_str("INVALID").is_err());
    }

    #[test]
    fn standard_residue_is_protein_classifies_correctly() {
        assert!(StandardResidue::ALA.is_protein());
        assert!(StandardResidue::GLY.is_protein());
        assert!(!StandardResidue::HOH.is_protein());
        assert!(!StandardResidue::A.is_protein());
    }

    #[test]
    fn standard_residue_is_nucleic_classifies_correctly() {
        assert!(StandardResidue::A.is_nucleic());
        assert!(StandardResidue::DT.is_nucleic());
        assert!(!StandardResidue::HOH.is_nucleic());
        assert!(!StandardResidue::GLY.is_nucleic());
    }
}
