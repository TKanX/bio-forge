//! PDB structure reader that normalizes chain ordering, residue templates, and metadata.
//!
//! The parser ingests legacy PDB records, applies `IoContext` aliasing to residues, filters
//! alternate locations by occupancy, and emits a fully linked [`Structure`] with terminal
//! classifications and optional unit-cell vectors.

use crate::io::context::IoContext;
use crate::io::error::Error;
use crate::model::{
    atom::Atom,
    chain::Chain,
    residue::Residue,
    structure::Structure,
    types::{Element, Point, ResidueCategory, ResiduePosition, StandardResidue},
};
use std::collections::{BTreeMap, HashMap};
use std::io::BufRead;
use std::path::Path;
use std::str::FromStr;

/// Composite key that uniquely identifies residues during the parsing pass.
///
/// Keeps chain ID, sequence number, and optional insertion code grouped together for
/// deterministic ordering while buffering atoms.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct ResKey {
    /// Chain identifier token exactly as seen in the PDB file.
    chain_id: String,
    /// Residue sequence number parsed from columns 23‑26.
    res_seq: i32,
    /// Insertion code used to disambiguate residues sharing the same `res_seq`.
    i_code: Option<char>,
}

/// Temporary residue buffer that holds atoms prior to canonicalization.
///
/// The parser collects all atoms for a residue, tracks whether it originated from an
/// `HETATM` record, and remembers the original residue name so the [`IoContext`] can map it
/// to standardized nomenclature.
struct TempResidue {
    /// Raw residue name as it appeared in the source file.
    raw_name: String,
    /// Indicates if the residue was described via an `HETATM` entry.
    is_hetatm: bool,
    /// Atom table keyed by atom name with occupancy values for altloc filtering.
    atoms: HashMap<String, (f64, Atom)>,
}

/// Parses a legacy PDB stream into a [`Structure`] using the supplied IO context.
///
/// The routine supports unit-cell (`CRYST1`) records, alternate locations (keeps the highest
/// occupancy), residue aliasing via [`IoContext`], and terminal classification for polymers.
///
/// # Arguments
///
/// * `reader` - Any buffered reader that yields PDB lines.
/// * `context` - Lookup tables and alias mappings that normalize residue names.
///
/// # Returns
///
/// A populated [`Structure`] containing chains, residues, atoms, and optional box vectors.
///
/// # Errors
///
/// Returns [`Error`] when encountering malformed numeric fields, unknown standard residues,
/// inconsistent unit-cell parameters, or IO failures from the underlying reader.
///
/// # Examples
///
/// ```
/// use bio_forge::io::{read_pdb_structure, IoContext};
/// use std::io::Cursor;
///
/// let pdb = "\
/// ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00 20.00           N\n\
/// END\n";
/// let mut cursor = Cursor::new(pdb.as_bytes());
/// let structure = read_pdb_structure(&mut cursor, &IoContext::new_default()).unwrap();
/// assert_eq!(structure.chain_count(), 1);
/// assert_eq!(structure.residue_count(), 1);
/// ```
pub fn read<R: BufRead>(reader: R, context: &IoContext) -> Result<Structure, Error> {
    let mut structure = Structure::new();

    let mut chain_order: Vec<String> = Vec::new();
    let mut chain_map: HashMap<String, BTreeMap<ResKey, TempResidue>> = HashMap::new();

    let mut line_num = 0;

    for line in reader.lines() {
        line_num += 1;
        let line = line.map_err(|e| Error::from_io(e, None))?;

        if line.starts_with("CRYST1") {
            structure.box_vectors = Some(parse_cryst1(&line, line_num)?);
            continue;
        }

        let is_atom = line.starts_with("ATOM  ");
        let is_hetatm = line.starts_with("HETATM");

        if is_atom || is_hetatm {
            parse_atom_record(&line, line_num, is_hetatm, &mut chain_order, &mut chain_map)?;
        }
    }

    for chain_id in chain_order {
        if let Some(residues) = chain_map.remove(&chain_id) {
            let mut chain = Chain::new(&chain_id);

            for (res_key, temp_res) in residues {
                let (canonical_name, std_enum) = context.classify_residue(&temp_res.raw_name);

                let category = determine_category(
                    temp_res.is_hetatm,
                    std_enum,
                    temp_res.atoms.len(),
                    &temp_res.raw_name,
                    None,
                )?;

                let mut residue = Residue::new(
                    res_key.res_seq,
                    res_key.i_code,
                    canonical_name.as_str(),
                    std_enum,
                    category,
                );

                let mut sorted_atoms: Vec<Atom> =
                    temp_res.atoms.into_values().map(|v| v.1).collect();
                sorted_atoms.sort_by(|a, b| a.name.cmp(&b.name));

                for atom in sorted_atoms {
                    residue.add_atom(atom);
                }

                chain.add_residue(residue);
            }
            structure.add_chain(chain);
        }
    }

    calculate_residue_positions(&mut structure);

    Ok(structure)
}

/// Parses an `ATOM`/`HETATM` line and updates the buffered residue map.
///
/// # Arguments
///
/// * `line` - Raw PDB record line.
/// * `line_num` - Current line number for diagnostics.
/// * `is_hetatm` - Indicates whether the record originated from `HETATM`.
/// * `chain_order` - Preserves the encounter order of chains.
/// * `chain_map` - Aggregates temporary residues keyed by [`ResKey`].
///
/// # Returns
///
/// [`Ok`] on successful parsing; [`Error`] if numeric fields are malformed or lines are too short.
fn parse_atom_record(
    line: &str,
    line_num: usize,
    is_hetatm: bool,
    chain_order: &mut Vec<String>,
    chain_map: &mut HashMap<String, BTreeMap<ResKey, TempResidue>>,
) -> Result<(), Error> {
    if line.len() < 54 {
        return Err(Error::parse("PDB", None, line_num, "Atom record too short"));
    }

    let atom_field = &line[12..16];
    let atom_name = atom_field.trim().to_string();
    let _alt_loc = line.chars().nth(16).unwrap_or(' ');
    let res_name = line[17..20].trim().to_string();
    let chain_id = line.chars().nth(21).unwrap_or(' ').to_string();
    let res_seq_str = &line[22..26];
    let i_code_char = line.chars().nth(26).unwrap_or(' ');

    let x_str = &line[30..38];
    let y_str = &line[38..46];
    let z_str = &line[46..54];

    let occ_str = if line.len() >= 60 {
        &line[54..60]
    } else {
        "1.00"
    };
    let element_str = if line.len() >= 78 {
        &line[76..78]
    } else {
        "  "
    };

    let res_seq = res_seq_str
        .trim()
        .parse::<i32>()
        .map_err(|_| Error::parse("PDB", None, line_num, "Invalid residue sequence number"))?;

    let i_code = if i_code_char == ' ' {
        None
    } else {
        Some(i_code_char)
    };

    let x = x_str
        .trim()
        .parse::<f64>()
        .map_err(|_| Error::parse("PDB", None, line_num, "Invalid X coordinate"))?;
    let y = y_str
        .trim()
        .parse::<f64>()
        .map_err(|_| Error::parse("PDB", None, line_num, "Invalid Y coordinate"))?;
    let z = z_str
        .trim()
        .parse::<f64>()
        .map_err(|_| Error::parse("PDB", None, line_num, "Invalid Z coordinate"))?;
    let pos = Point::new(x, y, z);

    let occupancy = occ_str.trim().parse::<f64>().unwrap_or(1.0);

    let element = if element_str.trim().is_empty() {
        parse_element_from_name(atom_field)
    } else {
        Element::from_str(element_str.trim()).unwrap_or(Element::Unknown)
    };

    if !chain_map.contains_key(&chain_id) {
        chain_map.insert(chain_id.clone(), BTreeMap::new());
        chain_order.push(chain_id.clone());
    }

    let residues = chain_map.get_mut(&chain_id).unwrap();
    let res_key = ResKey {
        chain_id: chain_id.clone(),
        res_seq,
        i_code,
    };

    let temp_res = residues.entry(res_key).or_insert_with(|| TempResidue {
        raw_name: res_name,
        is_hetatm,
        atoms: HashMap::new(),
    });

    match temp_res.atoms.get(&atom_name) {
        Some((old_occ, _)) => {
            if occupancy > *old_occ {
                temp_res.atoms.insert(
                    atom_name.clone(),
                    (occupancy, Atom::new(&atom_name, element, pos)),
                );
            }
        }
        None => {
            temp_res.atoms.insert(
                atom_name.clone(),
                (occupancy, Atom::new(&atom_name, element, pos)),
            );
        }
    }

    Ok(())
}

/// Converts a `CRYST1` record into orthogonal box vectors.
///
/// # Arguments
///
/// * `line` - Raw `CRYST1` text line.
/// * `line_num` - Current line number for error messages.
///
/// # Returns
///
/// A 3×3 matrix of cell vectors expressed in ångströms.
///
/// # Errors
///
/// Emits [`Error::InconsistentData`] when cell edges are zero or the line is underspecified.
fn parse_cryst1(line: &str, line_num: usize) -> Result<[[f64; 3]; 3], Error> {
    if line.len() < 54 {
        return Err(Error::parse(
            "PDB",
            None,
            line_num,
            "CRYST1 record too short",
        ));
    }

    let a = line[6..15].trim().parse::<f64>().unwrap_or(0.0);
    let b = line[15..24].trim().parse::<f64>().unwrap_or(0.0);
    let c = line[24..33].trim().parse::<f64>().unwrap_or(0.0);
    let alpha = line[33..40]
        .trim()
        .parse::<f64>()
        .unwrap_or(90.0)
        .to_radians();
    let beta = line[40..47]
        .trim()
        .parse::<f64>()
        .unwrap_or(90.0)
        .to_radians();
    let gamma = line[47..54]
        .trim()
        .parse::<f64>()
        .unwrap_or(90.0)
        .to_radians();

    if a == 0.0 || b == 0.0 || c == 0.0 {
        return Err(Error::inconsistent_data(
            "PDB",
            None,
            "Invalid unit cell dimensions",
        ));
    }

    let cos_a = alpha.cos();
    let cos_b = beta.cos();
    let cos_g = gamma.cos();
    let sin_g = gamma.sin();

    let v1_x = a;
    let v1_y = 0.0;
    let v1_z = 0.0;

    let v2_x = b * cos_g;
    let v2_y = b * sin_g;
    let v2_z = 0.0;

    let term = (cos_a - cos_b * cos_g) / sin_g;
    let v3_x = c * cos_b;
    let v3_y = c * term;
    let v3_z = c * (1.0 - cos_b * cos_b - term * term).sqrt();

    Ok([[v1_x, v1_y, v1_z], [v2_x, v2_y, v2_z], [v3_x, v3_y, v3_z]])
}

/// Infers an element symbol from an atom name when columns 77‑78 are blank.
///
/// Strips non-alphabetic characters, tries two-letter, then single-letter lookups, and
/// ultimately falls back to [`Element::Unknown`].
fn parse_element_from_name(field: &str) -> Element {
    fn parse_single(ch: char) -> Option<Element> {
        if !ch.is_ascii_alphabetic() {
            return None;
        }

        let symbol = ch.to_ascii_uppercase().to_string();
        Element::from_str(&symbol)
            .ok()
            .filter(|el| *el != Element::Unknown)
    }

    fn parse_pair(first: char, second: char) -> Option<Element> {
        if !(first.is_ascii_alphabetic() && second.is_ascii_alphabetic()) {
            return None;
        }

        let symbol = format!(
            "{}{}",
            first.to_ascii_uppercase(),
            second.to_ascii_lowercase()
        );
        Element::from_str(&symbol)
            .ok()
            .filter(|el| *el != Element::Unknown)
    }

    let letters: Vec<(usize, char)> = field
        .char_indices()
        .filter(|(_, ch)| ch.is_ascii_alphabetic())
        .collect();

    if letters.is_empty() {
        return Element::Unknown;
    }

    if letters.len() >= 2 {
        let (first_idx, first_char) = letters[0];
        let (second_idx, second_char) = letters[1];
        let contiguous = second_idx == first_idx + first_char.len_utf8();
        if let Some(el) = (first_idx == 0 && contiguous)
            .then(|| parse_pair(first_char, second_char))
            .flatten()
        {
            return el;
        }
    }

    if let Some(el) = parse_single(letters[0].1) {
        return el;
    }

    if letters.len() >= 2 {
        for window in letters.windows(2) {
            let (first_idx, first_char) = window[0];
            let (second_idx, second_char) = window[1];
            let contiguous = second_idx == first_idx + first_char.len_utf8();
            if let Some(el) = contiguous
                .then(|| parse_pair(first_char, second_char))
                .flatten()
            {
                return el;
            }
            if let Some(el) = parse_single(second_char) {
                return el;
            }
        }
    }

    Element::Unknown
}

/// Assigns a residue category using template matches or heuristics.
///
/// Standard templates override everything else. Pure `ATOM` residues without templates
/// trigger an error to protect canonical polymer expectations, while `HETATM` residues are
/// differentiated by atom count to detect ions.
///
/// # Returns
///
/// The chosen [`ResidueCategory`] or an [`Error::UnknownStandardResidue`] for unmapped
/// standard residues.
fn determine_category(
    is_hetatm: bool,
    std_enum: Option<StandardResidue>,
    atom_count: usize,
    res_name: &str,
    path: Option<&Path>,
) -> Result<ResidueCategory, Error> {
    if std_enum.is_some() {
        return Ok(ResidueCategory::Standard);
    }

    if !is_hetatm {
        return Err(Error::unknown_standard_residue(
            res_name,
            path.map(|p| p.to_path_buf()),
        ));
    }

    if atom_count == 1 {
        Ok(ResidueCategory::Ion)
    } else {
        Ok(ResidueCategory::Hetero)
    }
}

/// Updates residue positional annotations (N/C terminus, 5'/3', etc.).
///
/// The algorithm locates the first/last polymer residues within each chain (excluding water
/// and heterogens) and marks termini based on whether the polymer is protein or nucleic.
fn calculate_residue_positions(structure: &mut Structure) {
    for chain in structure.iter_chains_mut() {
        let polymer_indices: Vec<usize> = chain
            .iter_residues()
            .enumerate()
            .filter(|(_, res)| {
                res.category == ResidueCategory::Standard
                    && res.standard_name != Some(StandardResidue::HOH)
            })
            .map(|(i, _)| i)
            .collect();

        if polymer_indices.is_empty() {
            continue;
        }

        let first_idx = *polymer_indices.first().unwrap();
        let last_idx = *polymer_indices.last().unwrap();

        for (i, residue) in chain.iter_residues_mut().enumerate() {
            if residue.category != ResidueCategory::Standard
                || residue.standard_name == Some(StandardResidue::HOH)
            {
                residue.position = ResiduePosition::None;
                continue;
            }

            let is_protein = matches!(
                residue.standard_name,
                Some(std) if std.is_protein()
            );
            let is_nucleic = matches!(
                residue.standard_name,
                Some(std) if std.is_nucleic()
            );

            if i == first_idx {
                if is_protein {
                    residue.position = ResiduePosition::NTerminal;
                } else if is_nucleic {
                    residue.position = ResiduePosition::FivePrime;
                } else {
                    residue.position = ResiduePosition::Internal;
                }
            } else if i == last_idx {
                if is_protein {
                    residue.position = ResiduePosition::CTerminal;
                } else if is_nucleic {
                    residue.position = ResiduePosition::ThreePrime;
                } else {
                    residue.position = ResiduePosition::Internal;
                }
            } else if polymer_indices.contains(&i) {
                residue.position = ResiduePosition::Internal;
            } else {
                residue.position = ResiduePosition::None;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::context::IoContext;
    use crate::model::types::{Element, ResidueCategory, ResiduePosition, StandardResidue};
    use std::io::Cursor;

    fn parse_structure(pdb: &str) -> Structure {
        let mut cursor = Cursor::new(pdb.as_bytes());
        let context = IoContext::new_default();
        read(&mut cursor, &context).expect("PDB should parse")
    }

    fn parse_result(pdb: &str) -> Result<Structure, Error> {
        let mut cursor = Cursor::new(pdb.as_bytes());
        let context = IoContext::new_default();
        read(&mut cursor, &context)
    }

    #[test]
    fn read_parses_standard_polymer_and_box_vectors() {
        const PDB_DATA: &str = "\
            CRYST1   10.000   12.000   15.000  90.00  90.00  90.00 P 1           1\n\
            ATOM      1  N   ALA A   1      12.546  11.406   2.324  1.00 20.00           N\n\
            ATOM      2  CA  ALA A   1      13.123  12.345   3.210  1.00 20.00           C\n\
            ATOM      3  C   ALA A   1      14.456  11.987   4.123  1.00 20.00           C\n\
            ATOM      4  O   ALA A   1      15.123  12.456   4.987  1.00 20.00           O\n\
            ATOM      5  N   GLY A   2      14.789  10.654   4.890  1.00 20.00           N\n\
            ATOM      6  CA  GLY A   2      15.234  10.123   5.789  1.00 20.00           C\n\
            ATOM      7  C   GLY A   2      16.678  10.567   6.123  1.00 20.00           C\n\
            ATOM      8  O   GLY A   2      17.345  10.890   5.345  1.00 20.00           O\n\
            END\n";

        let structure = parse_structure(PDB_DATA);

        assert_eq!(structure.chain_count(), 1);
        assert_eq!(structure.residue_count(), 2);

        let box_vectors = structure
            .box_vectors
            .as_ref()
            .expect("CRYST1 record should set box vectors");
        assert!((box_vectors[0][0] - 10.0).abs() < 1e-6);
        assert!((box_vectors[1][1] - 12.0).abs() < 1e-6);
        assert!((box_vectors[2][2] - 15.0).abs() < 1e-6);

        let chain = structure.chain("A").expect("chain A exists");
        let residues: Vec<_> = chain.iter_residues().collect();
        assert_eq!(residues.len(), 2);

        let ala = residues[0];
        assert_eq!(ala.name, "ALA");
        assert_eq!(ala.standard_name, Some(StandardResidue::ALA));
        assert_eq!(ala.category, ResidueCategory::Standard);
        assert_eq!(ala.position, ResiduePosition::NTerminal);

        let gly = residues[1];
        assert_eq!(gly.name, "GLY");
        assert_eq!(gly.standard_name, Some(StandardResidue::GLY));
        assert_eq!(gly.category, ResidueCategory::Standard);
        assert_eq!(gly.position, ResiduePosition::CTerminal);
    }

    #[test]
    fn read_aliases_water_and_applies_occupancy_filter() {
        const PDB_DATA: &str = "\
            HETATM    1  O   WAT B   5       0.000   0.000   0.000  0.30 20.00           O\n\
            HETATM    2  O   WAT B   5       1.000   1.000   1.000  0.80 20.00           O\n\
            HETATM    3  NA  NA  B   6       5.000   5.000   5.000  1.00 20.00          NA\n";

        let structure = parse_structure(PDB_DATA);

        let chain = structure.chain("B").expect("chain B exists");
        assert_eq!(chain.residue_count(), 2);

        let water = chain.residue(5, None).expect("water residue exists");
        assert_eq!(water.name, "HOH");
        assert_eq!(water.standard_name, Some(StandardResidue::HOH));
        assert_eq!(water.category, ResidueCategory::Standard);
        let oxygen = water.atom("O").expect("oxygen atom kept");
        assert!((oxygen.pos.x - 1.0).abs() < 1e-6);
        assert!((oxygen.pos.y - 1.0).abs() < 1e-6);
        assert!((oxygen.pos.z - 1.0).abs() < 1e-6);

        let ion = chain.residue(6, None).expect("ion residue exists");
        assert_eq!(ion.category, ResidueCategory::Ion);
        assert_eq!(ion.atom_count(), 1);
        assert_eq!(ion.atom("NA").unwrap().name, "NA");
    }

    #[test]
    fn read_handles_scrambled_chain_and_residue_records() {
        const PDB_DATA: &str = "\
            ATOM      5  CA  GLY B   2      14.000  10.000   9.000  1.00 15.00           C\n\
            ATOM      6  C   GLY B   2      14.500  10.500   9.500  1.00 15.00           C\n\
            ATOM      1  N   ALA A   1      11.000  12.000   5.000  1.00 10.00           N\n\
            ATOM      2  CA  ALA A   1      11.500  12.500   5.500  1.00 10.00           C\n\
            ATOM      3  N   GLY B   1      13.000  11.000   8.000  1.00 15.00           N\n\
            ATOM      4  CA  GLY B   1      13.500  11.500   8.500  1.00 15.00           C\n\
            ATOM      7  N   ALA A   2      12.000  13.000   6.000  1.00 10.00           N\n\
            ATOM      8  CA  ALA A   2      12.500  13.500   6.500  1.00 10.00           C\n";

        let structure = parse_structure(PDB_DATA);
        let chain_ids: Vec<_> = structure.iter_chains().map(|c| c.id.clone()).collect();
        assert_eq!(
            chain_ids,
            vec!["B", "A"],
            "chain order should follow first appearance"
        );

        let chain_b = structure.chain("B").unwrap();
        let resid_b: Vec<_> = chain_b.iter_residues().map(|r| r.id).collect();
        assert_eq!(resid_b, vec![1, 2]);

        let chain_a = structure.chain("A").unwrap();
        let resid_a: Vec<_> = chain_a.iter_residues().map(|r| r.id).collect();
        assert_eq!(resid_a, vec![1, 2]);
    }

    #[test]
    fn read_assigns_terminal_positions_for_proteins_and_nucleic_acids() {
        const PDB_DATA: &str = "\
            ATOM      1  N   GLY P   1       1.000   2.000   3.000  1.00 10.00           N\n\
            ATOM      2  CA  GLY P   1       1.500   2.500   3.500  1.00 10.00           C\n\
            ATOM      3  C   GLY P   1       2.000   3.000   4.000  1.00 10.00           C\n\
            HETATM    4  O   WAT P   2       5.000   5.000   5.000  1.00 20.00           O\n\
            ATOM      5  N   SER P   3       2.500   3.500   4.500  1.00 10.00           N\n\
            ATOM      6  CA  SER P   3       3.000   4.000   5.000  1.00 10.00           C\n\
            ATOM      7  C   SER P   3       3.500   4.500   5.500  1.00 10.00           C\n\
            ATOM      8  N   LEU P   4       4.000   5.000   6.000  1.00 10.00           N\n\
            ATOM      9  CA  LEU P   4       4.500   5.500   6.500  1.00 10.00           C\n\
            ATOM     10  C   LEU P   4       5.000   6.000   7.000  1.00 10.00           C\n\
            ATOM     11  P   DA  N   1       6.000   1.000   1.000  1.00 20.00           P\n\
            ATOM     12  O5' DA  N   1       6.500   1.500   1.500  1.00 20.00           O\n\
            ATOM     13  P   DT  N   2       7.000   2.000   2.000  1.00 20.00           P\n\
            ATOM     14  O5' DT  N   2       7.500   2.500   2.500  1.00 20.00           O\n";

        let structure = parse_structure(PDB_DATA);

        let protein = structure.chain("P").expect("protein chain present");
        let residues: Vec<_> = protein.iter_residues().collect();
        assert_eq!(residues[0].position, ResiduePosition::NTerminal);
        assert_eq!(
            residues[1].position,
            ResiduePosition::None,
            "water should not be classified as polymer"
        );
        assert_eq!(residues[2].position, ResiduePosition::Internal);
        assert_eq!(residues[3].position, ResiduePosition::CTerminal);

        let nucleic = structure.chain("N").expect("nucleic chain present");
        let n_res: Vec<_> = nucleic.iter_residues().collect();
        assert_eq!(n_res[0].position, ResiduePosition::FivePrime);
        assert_eq!(n_res[1].position, ResiduePosition::ThreePrime);
    }

    #[test]
    fn read_categorizes_residues_based_on_templates_and_atom_counts() {
        const PDB_DATA: &str = "\
            ATOM      1  N   GLY C   1       9.000   9.000   9.000  1.00 12.00           N\n\
            ATOM      2  CA  GLY C   1       9.500   9.500   9.500  1.00 12.00           C\n\
            HETATM    3  O   WAT C   2       5.000   5.000   5.000  1.00 20.00           O\n\
            HETATM    4  C1  LIG C 301       1.000   1.000   1.000  1.00 30.00           C\n\
            HETATM    5  O1  LIG C 301       1.500   1.500   1.500  1.00 30.00           O\n\
            HETATM    6  NA  NA  C 401       2.000   2.000   2.000  1.00 20.00          NA\n";

        let structure = parse_structure(PDB_DATA);
        let chain = structure.chain("C").unwrap();

        let protein = chain.residue(1, None).unwrap();
        assert_eq!(protein.category, ResidueCategory::Standard);

        let water = chain.residue(2, None).unwrap();
        assert_eq!(water.category, ResidueCategory::Standard);
        assert_eq!(water.standard_name, Some(StandardResidue::HOH));

        let ligand = chain.residue(301, None).unwrap();
        assert_eq!(ligand.category, ResidueCategory::Hetero);
        assert_eq!(ligand.atom_count(), 2);

        let ion = chain.residue(401, None).unwrap();
        assert_eq!(ion.category, ResidueCategory::Ion);
        assert_eq!(ion.atom_count(), 1);
    }

    #[test]
    fn read_prefers_atoms_with_highest_occupancy_and_retains_hydrogens() {
        const PDB_DATA: &str = "\
            ATOM      1  N   GLY D   1       0.000   0.000   0.000  1.00 12.00           N\n\
            ATOM      2  CA AGLY D   1       1.000   0.000   0.000  0.40 12.00           C\n\
            ATOM      3  CA BGLY D   1       2.000   0.000   0.000  0.80 12.00           C\n\
            ATOM      4  CA CGLY D   1       3.000   0.000   0.000  0.80 12.00           C\n\
            ATOM      5  H1  GLY D   1       0.500   0.500   0.500  1.00 12.00           H\n";

        let structure = parse_structure(PDB_DATA);
        let residue = structure.chain("D").unwrap().residue(1, None).unwrap();

        let ca = residue.atom("CA").unwrap();
        assert!(
            (ca.pos.x - 2.0).abs() < 1e-6,
            "highest occupancy coordinate should be retained"
        );
        assert_eq!(
            residue.atom_count(),
            3,
            "duplicate atom names must collapse via occupancy"
        );
        assert!(
            residue.atom("H1").is_some(),
            "hydrogen atoms should be preserved for standard residues"
        );
    }

    #[test]
    fn read_supports_residues_with_insertion_codes() {
        const PDB_DATA: &str = "\
            ATOM      1  N   SER E  10A      0.000   0.000   0.000  1.00 10.00           N\n\
            ATOM      2  CA  SER E  10A      0.500   0.500   0.500  1.00 10.00           C\n\
            ATOM      3  N   SER E  10       1.000   1.000   1.000  1.00 10.00           N\n\
            ATOM      4  CA  SER E  10       1.500   1.500   1.500  1.00 10.00           C\n\
            ATOM      5  N   SER E  11       2.000   2.000   2.000  1.00 10.00           N\n\
            ATOM      6  CA  SER E  11       2.500   2.500   2.500  1.00 10.00           C\n";

        let structure = parse_structure(PDB_DATA);
        let chain = structure.chain("E").expect("chain E exists");

        let residues: Vec<_> = chain
            .iter_residues()
            .map(|r| (r.id, r.insertion_code))
            .collect();
        assert_eq!(residues, vec![(10, None), (10, Some('A')), (11, None)]);

        let base = chain.residue(10, None).expect("base residue present");
        let insertion = chain
            .residue(10, Some('A'))
            .expect("insertion residue present");

        assert_ne!(
            base.atom("CA").unwrap().pos,
            insertion.atom("CA").unwrap().pos
        );
    }

    #[test]
    fn read_errors_on_unknown_standard_atom_record() {
        const PDB_DATA: &str = "\
            ATOM      1  N   UNK C   1       0.000   0.000   0.000  1.00 20.00           N\n";

        let err = parse_result(PDB_DATA).expect_err("unknown residue should fail");

        match err {
            Error::UnknownStandardResidue { name, path } => {
                assert_eq!(name, "UNK");
                assert!(path.is_none());
            }
            other => panic!("expected UnknownStandardResidue error, got {other:?}"),
        }
    }

    #[test]
    fn parse_element_handles_common_atom_name_patterns() {
        assert_eq!(parse_element_from_name(" CA "), Element::C);
        assert_eq!(parse_element_from_name(" N  "), Element::N);
        assert_eq!(parse_element_from_name(" C1 "), Element::C);
        assert_eq!(parse_element_from_name("1HG1"), Element::H);
        assert_eq!(parse_element_from_name(" HG "), Element::H);
        assert_eq!(parse_element_from_name("FE  "), Element::Fe);
        assert_eq!(parse_element_from_name("ZN  "), Element::Zn);
        assert_eq!(parse_element_from_name("BR  "), Element::Br);
        assert_eq!(parse_element_from_name("CL  "), Element::Cl);
        assert_eq!(parse_element_from_name("HG  "), Element::Hg);
        assert_eq!(parse_element_from_name("Se  "), Element::Se);
        assert_eq!(parse_element_from_name("XE  "), Element::Xe);
    }
}
