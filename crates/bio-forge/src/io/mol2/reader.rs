//! Parses Tripos MOL2 templates into [`Template`] instances.
//!
//! The reader focuses on ligand-sized building blocks that define atoms and bonds without
//! coordinates. It validates identifiers, ensures counts match declarations, and converts
//! textual bond orders into the canonical `BondOrder` enum so that templates integrate with
//! repair and solvation pipelines.

use crate::io::error::Error;
use crate::model::template::Template;
use crate::model::types::BondOrder;
use std::collections::{BTreeMap, HashSet};
use std::io::BufRead;

/// Identifier used in diagnostics to reference the Tripos MOL2 format.
const FORMAT: &str = "MOL2";

/// Reads a MOL2 ligand template from any buffered reader.
///
/// The parser walks the `@<TRIPOS>` sections in order, collecting atom names and bond
/// records and validating that identifiers, counts, and names meet the stricter
/// requirements imposed by `bio-forge` templates. Duplicate names or references to missing
/// atoms are rejected with descriptive parse errors, making it easier to diagnose corrupt
/// ligand libraries.
///
/// # Arguments
///
/// * `reader` - A buffered text source positioned at the start of a MOL2 record.
///
/// # Returns
///
/// A fully populated [`Template`] containing atom names and bond topology.
///
/// # Errors
///
/// Returns [`Error::Parse`] when the MOL2 syntax is invalid
/// and [`Error::InconsistentData`] when declared
/// counts or references do not match the parsed content.
///
/// # Examples
///
/// ```
/// use bio_forge::io::read_mol2_template;
/// use std::io::Cursor;
///
/// const LIGAND: &str = "\
///     @<TRIPOS>MOLECULE\n\
///     ETH\n\
///     2 1 0 0 0\n\
///     SMALL\n\
///     NO_CHARGES\n\
///     @<TRIPOS>ATOM\n\
///         1 C1 0.0 0.0 0.0 C.3\n\
///         2 C2 1.5 0.0 0.0 C.3\n\
///     @<TRIPOS>BOND\n\
///         1 1 2 1\n";
///
/// let template = read_mol2_template(Cursor::new(LIGAND.as_bytes())).unwrap();
/// assert_eq!(template.name, "ETH");
/// assert!(template.has_bond("C1", "C2"));
/// ```
pub fn read<R: BufRead>(reader: R) -> Result<Template, Error> {
    let mut section = Section::None;
    let mut molecule_lines_seen = 0usize;

    let mut template_name: Option<String> = None;
    let mut expected_atoms: Option<usize> = None;
    let mut expected_bonds: Option<usize> = None;

    let mut atom_names: BTreeMap<usize, String> = BTreeMap::new();
    let mut seen_atom_names: HashSet<String> = HashSet::new();
    let mut bond_records: Vec<(usize, usize, BondOrder)> = Vec::new();

    for (idx, line_res) in reader.lines().enumerate() {
        let line = line_res.map_err(|e| Error::from_io(e, None))?;
        let line_number = idx + 1;
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if let Some(header) = trimmed.strip_prefix("@<TRIPOS>") {
            section = match header {
                "MOLECULE" => {
                    molecule_lines_seen = 0;
                    Section::Molecule
                }
                "ATOM" => Section::Atom,
                "BOND" => Section::Bond,
                _ => Section::None,
            };
            continue;
        }

        match section {
            Section::Molecule => {
                molecule_lines_seen += 1;
                match molecule_lines_seen {
                    1 => {
                        if trimmed == "****" {
                            return Err(Error::parse(
                                FORMAT,
                                None,
                                line_number,
                                "MOL2 molecule name must not be empty",
                            ));
                        }
                        template_name = Some(trimmed.to_string());
                    }
                    2 => {
                        let (atoms, bonds) = parse_expected_counts(trimmed, line_number)?;
                        expected_atoms = Some(atoms);
                        expected_bonds = Some(bonds);
                    }
                    _ => {}
                }
            }
            Section::Atom => {
                let tokens: Vec<&str> = trimmed.split_whitespace().collect();
                if tokens.len() < 2 {
                    return Err(Error::parse(
                        FORMAT,
                        None,
                        line_number,
                        "ATOM record must include at least id and name",
                    ));
                }
                let atom_id = tokens[0]
                    .parse::<usize>()
                    .map_err(|_| Error::parse(FORMAT, None, line_number, "Invalid atom id"))?;
                if atom_id == 0 {
                    return Err(Error::parse(
                        FORMAT,
                        None,
                        line_number,
                        "Atom id must be positive",
                    ));
                }
                if atom_names.contains_key(&atom_id) {
                    return Err(Error::parse(
                        FORMAT,
                        None,
                        line_number,
                        format!("Duplicate atom id {atom_id}"),
                    ));
                }

                let atom_name = tokens[1];
                if atom_name == "****" {
                    return Err(Error::parse(
                        FORMAT,
                        None,
                        line_number,
                        "Atom name must not be ****",
                    ));
                }
                if !seen_atom_names.insert(atom_name.to_string()) {
                    return Err(Error::parse(
                        FORMAT,
                        None,
                        line_number,
                        format!("Duplicate atom name '{atom_name}'"),
                    ));
                }

                atom_names.insert(atom_id, atom_name.to_string());
            }
            Section::Bond => {
                let tokens: Vec<&str> = trimmed.split_whitespace().collect();
                if tokens.len() < 4 {
                    return Err(Error::parse(
                        FORMAT,
                        None,
                        line_number,
                        "BOND record must include id, endpoints, and bond type",
                    ));
                }

                let origin = tokens[1].parse::<usize>().map_err(|_| {
                    Error::parse(FORMAT, None, line_number, "Invalid origin atom id")
                })?;
                let target = tokens[2].parse::<usize>().map_err(|_| {
                    Error::parse(FORMAT, None, line_number, "Invalid target atom id")
                })?;
                if origin == 0 || target == 0 {
                    return Err(Error::parse(
                        FORMAT,
                        None,
                        line_number,
                        "Bond endpoints must reference positive atom ids",
                    ));
                }
                let order = parse_bond_order_token(tokens[3], line_number)?;
                bond_records.push((origin, target, order));
            }
            Section::None => {}
        }
    }

    let name = template_name.ok_or_else(|| {
        Error::parse(
            FORMAT,
            None,
            0,
            "Missing @<TRIPOS>MOLECULE section with molecule name",
        )
    })?;

    if atom_names.is_empty() {
        return Err(Error::parse(
            FORMAT,
            None,
            0,
            "Missing or empty @<TRIPOS>ATOM section",
        ));
    }

    if let Some(expected) = expected_atoms.filter(|&expected| expected != atom_names.len()) {
        return Err(Error::inconsistent_data(
            FORMAT,
            None,
            format!("Declared {expected} atoms but parsed {}", atom_names.len()),
        ));
    }

    if let Some(expected) = expected_bonds.filter(|&expected| expected != bond_records.len()) {
        return Err(Error::inconsistent_data(
            FORMAT,
            None,
            format!(
                "Declared {expected} bonds but parsed {}",
                bond_records.len()
            ),
        ));
    }

    let mut bonds = Vec::with_capacity(bond_records.len());
    for (origin, target, order) in bond_records {
        let a1 = atom_names.get(&origin).ok_or_else(|| {
            Error::inconsistent_data(
                FORMAT,
                None,
                format!("Bond references unknown atom id {origin}"),
            )
        })?;
        let a2 = atom_names.get(&target).ok_or_else(|| {
            Error::inconsistent_data(
                FORMAT,
                None,
                format!("Bond references unknown atom id {target}"),
            )
        })?;
        bonds.push((a1.clone(), a2.clone(), order));
    }

    let atom_name_list: Vec<String> = atom_names.values().cloned().collect();

    Ok(Template::new(name, atom_name_list, bonds))
}

/// Parser state corresponding to the current `@<TRIPOS>` section.
enum Section {
    /// No active section; lines are ignored until the next header.
    None,
    /// `@<TRIPOS>MOLECULE` block containing metadata and expected counts.
    Molecule,
    /// `@<TRIPOS>ATOM` block listing atom identifiers and labels.
    Atom,
    /// `@<TRIPOS>BOND` block listing bond endpoints and orders.
    Bond,
}

/// Parses the count line that follows `@<TRIPOS>MOLECULE`.
///
/// The routine extracts the declared atom and bond counts, validating that both are
/// positive integers so later comparisons can detect mismatches.
///
/// # Arguments
///
/// * `line` - The raw count line string with whitespace-separated integers.
/// * `line_number` - The one-based line number for diagnostics.
///
/// # Returns
///
/// A tuple with the expected atom and bond counts.
///
/// # Errors
///
/// Returns [`Error::Parse`] when either field is missing or
/// cannot be converted into a `usize`.
fn parse_expected_counts(line: &str, line_number: usize) -> Result<(usize, usize), Error> {
    let fields: Vec<&str> = line.split_whitespace().collect();
    if fields.len() < 2 {
        return Err(Error::parse(
            FORMAT,
            None,
            line_number,
            "Molecule counts line must specify atom and bond counts",
        ));
    }
    let atoms = fields[0]
        .parse::<usize>()
        .map_err(|_| Error::parse(FORMAT, None, line_number, "Invalid atom count"))?;
    let bonds = fields[1]
        .parse::<usize>()
        .map_err(|_| Error::parse(FORMAT, None, line_number, "Invalid bond count"))?;
    Ok((atoms, bonds))
}

/// Converts a MOL2 bond order token into a [`BondOrder`].
///
/// MOL2 encodes aromatic bonds as `ar`, so the helper lowercases the token and maps both
/// textual and numeric representations to the internal enum.
///
/// # Arguments
///
/// * `token` - The bond order string found in the BOND section.
/// * `line_number` - The one-based line number for diagnostics.
///
/// # Returns
///
/// The normalized [`BondOrder`] variant for the provided token.
///
/// # Errors
///
/// Returns [`Error::Parse`] for unsupported or malformed
/// tokens.
fn parse_bond_order_token(token: &str, line_number: usize) -> Result<BondOrder, Error> {
    match token.to_ascii_lowercase().as_str() {
        "1" | "single" => Ok(BondOrder::Single),
        "2" | "double" => Ok(BondOrder::Double),
        "3" | "triple" => Ok(BondOrder::Triple),
        "ar" | "aromatic" => Ok(BondOrder::Aromatic),
        other => Err(Error::parse(
            FORMAT,
            None,
            line_number,
            format!("Unsupported bond type '{other}'"),
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::error::Error;
    use std::io::Cursor;

    const BENZENE: &str = "\
        @<TRIPOS>MOLECULE\n\
        BNZ\n\
        12 12 0 0 0\n\
        SMALL\n\
        NO_CHARGES\n\
        ****\n\
        \n\
        @<TRIPOS>ATOM\n\
            1 C1        4.5700  21.9710   1.5490  C.ar    1 BNZ       0.000\n\
            2 C2        3.9410  23.1640   1.8820  C.ar    1 BNZ       0.000\n\
            3 C3        4.6620  24.3520   1.8890  C.ar    1 BNZ       0.000\n\
            4 C4        6.0260  24.3190   1.5870  C.ar    1 BNZ       0.000\n\
            5 C5        6.6430  23.1280   1.2540  C.ar    1 BNZ       0.000\n\
            6 C6        5.9080  21.9520   1.2180  C.ar    1 BNZ       0.000\n\
            7 H1        4.0060  21.0500   1.5490  H       1 BNZ       0.000\n\
            8 H2        2.8910  23.1680   2.1360  H       1 BNZ       0.000\n\
            9 H3        4.1750  25.2870   2.1240  H       1 BNZ       0.000\n\
            10 H4        6.6030  25.2320   1.6140  H       1 BNZ       0.000\n\
            11 H5        7.6980  23.1130   1.0220  H       1 BNZ       0.000\n\
            12 H6        6.3820  21.0250   0.9320  H       1 BNZ       0.000\n\
        \n\
        @<TRIPOS>BOND\n\
            1    1    2  ar\n\
            2    2    3  ar\n\
            3    3    4  ar\n\
            4    4    5  ar\n\
            5    5    6  ar\n\
            6    6    1  ar\n\
            7    1    7   1\n\
            8    2    8   1\n\
            9    3    9   1\n\
            10    4   10   1\n\
            11    5   11   1\n\
            12    6   12   1\n";

    #[test]
    fn read_parses_benzene_template() {
        let mut cursor = Cursor::new(BENZENE.as_bytes());
        let template = read(&mut cursor).expect("benzene should parse");

        assert_eq!(template.name, "BNZ");
        assert_eq!(template.atom_count(), 12);
        assert_eq!(template.bond_count(), 12);
        assert!(template.has_atom("C1"));
        assert!(template.has_atom("H6"));
        assert!(template.has_bond("C1", "C2"));
        assert!(template.has_bond("C3", "H3"));
    }

    #[test]
    fn read_errors_on_duplicate_atom_names() {
        const DUPLICATE: &str = "\
            @<TRIPOS>MOLECULE\n\
            LIG\n\
            2 1 0 0 0\n\
            SMALL\n\
            NO_CHARGES\n\
            @<TRIPOS>ATOM\n\
                1 C1 0 0 0 C.3\n\
                2 C1 1 0 0 C.3\n\
            @<TRIPOS>BOND\n\
                1 1 2 1\n";

        let mut cursor = Cursor::new(DUPLICATE.as_bytes());
        let err = read(&mut cursor).expect_err("duplicate atom names should fail");

        match err {
            Error::Parse { details, .. } => {
                assert!(details.contains("Duplicate atom name"));
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }

    #[test]
    fn read_errors_on_unknown_bond_atom() {
        const UNKNOWN: &str = "\
            @<TRIPOS>MOLECULE\n\
            LIG\n\
            2 1 0 0 0\n\
            SMALL\n\
            NO_CHARGES\n\
            @<TRIPOS>ATOM\n\
                1 C1 0 0 0 C.3\n\
                2 C2 1 0 0 C.3\n\
            @<TRIPOS>BOND\n\
                1 1 3 1\n";

        let mut cursor = Cursor::new(UNKNOWN.as_bytes());
        let err = read(&mut cursor).expect_err("unknown atom id should fail");

        match err {
            Error::InconsistentData { details, .. } => {
                assert!(details.contains("unknown atom id"));
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }

    #[test]
    fn read_errors_on_count_mismatch() {
        const COUNT: &str = "\
            @<TRIPOS>MOLECULE\n\
            LIG\n\
            3 0 0 0 0\n\
            SMALL\n\
            NO_CHARGES\n\
            @<TRIPOS>ATOM\n\
                1 C1 0 0 0 C.3\n\
                2 C2 1 0 0 C.3\n";

        let mut cursor = Cursor::new(COUNT.as_bytes());
        let err = read(&mut cursor).expect_err("count mismatch should fail");

        match err {
            Error::InconsistentData { details, .. } => {
                assert!(details.contains("Declared 3 atoms"));
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }
}
