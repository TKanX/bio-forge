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

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct ResKey {
    chain_id: String,
    res_seq: i32,
    i_code: Option<char>,
}

struct TempResidue {
    raw_name: String,
    is_hetatm: bool,
    atoms: HashMap<String, (f64, Atom)>,
}

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

    let atom_name = line[12..16].trim().to_string();
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
        parse_element_from_name(&atom_name)
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

fn parse_element_from_name(name: &str) -> Element {
    let name = name.trim();
    let mut symbol = String::new();
    for c in name.chars() {
        if c.is_alphabetic() {
            symbol.push(c);
        } else if !symbol.is_empty() {
            break;
        }
    }

    if let Ok(el) = Element::from_str(&symbol) {
        return el;
    }
    if !symbol.is_empty() {
        let first = &symbol[0..1];
        if let Ok(el) = Element::from_str(first) {
            return el;
        }
    }

    Element::Unknown
}

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
