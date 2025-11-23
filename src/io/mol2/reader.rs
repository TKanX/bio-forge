use crate::io::error::Error;
use crate::model::template::Template;
use crate::model::types::BondOrder;
use std::collections::{BTreeMap, HashSet};
use std::io::BufRead;

const FORMAT: &str = "MOL2";

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

    if let Some(expected) = expected_atoms {
        if expected != atom_names.len() {
            return Err(Error::inconsistent_data(
                FORMAT,
                None,
                format!("Declared {expected} atoms but parsed {}", atom_names.len()),
            ));
        }
    }

    if let Some(expected) = expected_bonds {
        if expected != bond_records.len() {
            return Err(Error::inconsistent_data(
                FORMAT,
                None,
                format!(
                    "Declared {expected} bonds but parsed {}",
                    bond_records.len()
                ),
            ));
        }
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

enum Section {
    None,
    Molecule,
    Atom,
    Bond,
}

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
