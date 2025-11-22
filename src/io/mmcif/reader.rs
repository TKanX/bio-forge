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

#[derive(Default)]
struct AtomSiteIndices {
    group_pdb: Option<usize>,
    auth_atom_id: Option<usize>,
    label_atom_id: Option<usize>,
    auth_comp_id: Option<usize>,
    label_comp_id: Option<usize>,
    auth_asym_id: Option<usize>,
    label_asym_id: Option<usize>,
    auth_seq_id: Option<usize>,
    label_seq_id: Option<usize>,
    pdbx_pdb_ins_code: Option<usize>,
    cartn_x: Option<usize>,
    cartn_y: Option<usize>,
    cartn_z: Option<usize>,
    occupancy: Option<usize>,
    type_symbol: Option<usize>,
}

enum ParserState {
    Base,
    InLoopHeader,
    InAtomSiteLoop,
    InOtherLoop,
}

pub fn read<R: BufRead>(reader: R, context: &IoContext) -> Result<Structure, Error> {
    let mut structure = Structure::new();

    let mut chain_order: Vec<String> = Vec::new();
    let mut chain_map: HashMap<String, BTreeMap<ResKey, TempResidue>> = HashMap::new();

    let mut state = ParserState::Base;
    let mut atom_indices = AtomSiteIndices::default();
    let mut current_loop_headers = Vec::new();
    let mut line_num = 0;

    let mut cell_params = HashMap::new();

    for line in reader.lines() {
        line_num += 1;
        let line = line.map_err(|e| Error::from_io(e, None))?;
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let tokens = tokenize_mmcif_line(trimmed);
        if tokens.is_empty() {
            continue;
        }

        if tokens[0] == "loop_" {
            state = ParserState::InLoopHeader;
            current_loop_headers.clear();
            continue;
        }

        if tokens[0].starts_with("_cell.") {
            if tokens.len() >= 2 {
                cell_params.insert(tokens[0].clone(), tokens[1].clone());
            }
            continue;
        }

        match state {
            ParserState::Base => {}
            ParserState::InLoopHeader => {
                if tokens[0].starts_with('_') {
                    current_loop_headers.push(tokens[0].clone());
                } else {
                    if current_loop_headers
                        .iter()
                        .any(|h| h.starts_with("_atom_site."))
                    {
                        state = ParserState::InAtomSiteLoop;
                        atom_indices = map_atom_site_indices(&current_loop_headers);
                        process_atom_line(
                            &tokens,
                            &atom_indices,
                            line_num,
                            &mut chain_order,
                            &mut chain_map,
                        )?;
                    } else {
                        state = ParserState::InOtherLoop;
                    }
                }
            }
            ParserState::InAtomSiteLoop => {
                if tokens[0].starts_with('_') || tokens[0] == "loop_" {
                    state = if tokens[0] == "loop_" {
                        ParserState::InLoopHeader
                    } else {
                        ParserState::Base
                    };
                } else {
                    process_atom_line(
                        &tokens,
                        &atom_indices,
                        line_num,
                        &mut chain_order,
                        &mut chain_map,
                    )?;
                }
            }
            ParserState::InOtherLoop => {
                if tokens[0].starts_with('_') || tokens[0] == "loop_" {
                    state = if tokens[0] == "loop_" {
                        ParserState::InLoopHeader
                    } else {
                        ParserState::Base
                    };
                }
            }
        }
    }

    if let Some(box_vectors) = process_cell_parameters(&cell_params) {
        structure.box_vectors = Some(box_vectors);
    }

    build_structure(structure, chain_order, chain_map, context)
}

fn tokenize_mmcif_line(line: &str) -> Vec<String> {
    let mut tokens = Vec::new();
    let mut current = String::new();
    let mut in_quote = None;

    for c in line.chars() {
        match in_quote {
            Some(q) => {
                if c == q {
                    in_quote = None;
                } else {
                    current.push(c);
                }
            }
            None => {
                if c.is_whitespace() {
                    if !current.is_empty() {
                        tokens.push(current.clone());
                        current.clear();
                    }
                } else if c == '\'' || c == '"' {
                    if !current.is_empty() {
                        current.push(c);
                    } else {
                        in_quote = Some(c);
                    }
                } else {
                    current.push(c);
                }
            }
        }
    }
    if !current.is_empty() {
        tokens.push(current);
    }
    tokens
}

fn map_atom_site_indices(headers: &[String]) -> AtomSiteIndices {
    let mut indices = AtomSiteIndices::default();
    for (i, header) in headers.iter().enumerate() {
        match header.as_str() {
            "_atom_site.group_PDB" => indices.group_pdb = Some(i),
            "_atom_site.auth_atom_id" => indices.auth_atom_id = Some(i),
            "_atom_site.label_atom_id" => indices.label_atom_id = Some(i),
            "_atom_site.auth_comp_id" => indices.auth_comp_id = Some(i),
            "_atom_site.label_comp_id" => indices.label_comp_id = Some(i),
            "_atom_site.auth_asym_id" => indices.auth_asym_id = Some(i),
            "_atom_site.label_asym_id" => indices.label_asym_id = Some(i),
            "_atom_site.auth_seq_id" => indices.auth_seq_id = Some(i),
            "_atom_site.label_seq_id" => indices.label_seq_id = Some(i),
            "_atom_site.pdbx_PDB_ins_code" => indices.pdbx_pdb_ins_code = Some(i),
            "_atom_site.Cartn_x" => indices.cartn_x = Some(i),
            "_atom_site.Cartn_y" => indices.cartn_y = Some(i),
            "_atom_site.Cartn_z" => indices.cartn_z = Some(i),
            "_atom_site.occupancy" => indices.occupancy = Some(i),
            "_atom_site.type_symbol" => indices.type_symbol = Some(i),
            _ => {}
        }
    }
    indices
}

fn token<'a>(tokens: &'a [String], idx: usize, line_num: usize) -> Result<&'a str, Error> {
    tokens.get(idx).map(|s| s.as_str()).ok_or_else(|| {
        Error::parse(
            "mmCIF",
            None,
            line_num,
            "Atom record shorter than _atom_site definition",
        )
    })
}

fn optional_token<'a>(
    tokens: &'a [String],
    idx: Option<usize>,
    line_num: usize,
) -> Result<Option<&'a str>, Error> {
    if let Some(idx) = idx {
        token(tokens, idx, line_num).map(Some)
    } else {
        Ok(None)
    }
}

fn parse_coordinate(value: &str, axis: &str, line_num: usize) -> Result<f64, Error> {
    f64::from_str(value).map_err(|_| {
        Error::parse(
            "mmCIF",
            None,
            line_num,
            format!("Invalid {axis} coordinate"),
        )
    })
}

fn process_atom_line(
    tokens: &[String],
    indices: &AtomSiteIndices,
    line_num: usize,
    chain_order: &mut Vec<String>,
    chain_map: &mut HashMap<String, BTreeMap<ResKey, TempResidue>>,
) -> Result<(), Error> {
    let atom_name_idx = indices
        .auth_atom_id
        .or(indices.label_atom_id)
        .ok_or_else(|| {
            Error::parse(
                "mmCIF",
                None,
                line_num,
                "_atom_site loop is missing atom identifier columns",
            )
        })?;
    let res_name_idx = indices
        .auth_comp_id
        .or(indices.label_comp_id)
        .ok_or_else(|| {
            Error::parse(
                "mmCIF",
                None,
                line_num,
                "_atom_site loop is missing residue identifier columns",
            )
        })?;
    let chain_id_idx = indices
        .auth_asym_id
        .or(indices.label_asym_id)
        .ok_or_else(|| {
            Error::parse(
                "mmCIF",
                None,
                line_num,
                "_atom_site loop is missing chain identifier columns",
            )
        })?;
    let seq_id_idx = indices
        .auth_seq_id
        .or(indices.label_seq_id)
        .ok_or_else(|| {
            Error::parse(
                "mmCIF",
                None,
                line_num,
                "_atom_site loop is missing residue sequence columns",
            )
        })?;

    let x_idx = indices.cartn_x.ok_or_else(|| {
        Error::parse(
            "mmCIF",
            None,
            line_num,
            "_atom_site.Cartn_x column is required",
        )
    })?;
    let y_idx = indices.cartn_y.ok_or_else(|| {
        Error::parse(
            "mmCIF",
            None,
            line_num,
            "_atom_site.Cartn_y column is required",
        )
    })?;
    let z_idx = indices.cartn_z.ok_or_else(|| {
        Error::parse(
            "mmCIF",
            None,
            line_num,
            "_atom_site.Cartn_z column is required",
        )
    })?;

    let required_indices = [
        atom_name_idx,
        res_name_idx,
        chain_id_idx,
        seq_id_idx,
        x_idx,
        y_idx,
        z_idx,
    ];
    if let Some(max_idx) = required_indices.iter().copied().max() {
        if tokens.len() <= max_idx {
            return Err(Error::parse(
                "mmCIF",
                None,
                line_num,
                "Atom record is shorter than declared _atom_site headers",
            ));
        }
    }

    let group_pdb = optional_token(tokens, indices.group_pdb, line_num)?;
    let is_hetatm = matches!(group_pdb, Some(val) if val.eq_ignore_ascii_case("HETATM"));

    let atom_name = token(tokens, atom_name_idx, line_num)?;
    let res_name = token(tokens, res_name_idx, line_num)?;
    let chain_id_raw = token(tokens, chain_id_idx, line_num)?;
    let seq_id_str = token(tokens, seq_id_idx, line_num)?;
    let ins_code_str = optional_token(tokens, indices.pdbx_pdb_ins_code, line_num)?;

    let x_str = token(tokens, x_idx, line_num)?;
    let y_str = token(tokens, y_idx, line_num)?;
    let z_str = token(tokens, z_idx, line_num)?;

    let occ_str = optional_token(tokens, indices.occupancy, line_num)?;
    let elem_str = optional_token(tokens, indices.type_symbol, line_num)?;

    if matches!(x_str, "." | "?") || matches!(y_str, "." | "?") || matches!(z_str, "." | "?") {
        return Ok(());
    }

    let x = parse_coordinate(x_str, "X", line_num)?;
    let y = parse_coordinate(y_str, "Y", line_num)?;
    let z = parse_coordinate(z_str, "Z", line_num)?;
    let pos = Point::new(x, y, z);

    let res_seq = if matches!(seq_id_str, "." | "?") {
        1
    } else {
        seq_id_str.parse::<i32>().unwrap_or(1)
    };

    let i_code = ins_code_str.and_then(|code| {
        if matches!(code, "." | "?") {
            None
        } else {
            code.chars().next()
        }
    });

    let occupancy = occ_str
        .filter(|occ| !matches!(*occ, "." | "?"))
        .and_then(|occ| f64::from_str(occ).ok())
        .unwrap_or(1.0);

    let element = elem_str
        .filter(|elem| !matches!(*elem, "." | "?"))
        .and_then(|elem| Element::from_str(elem).ok())
        .unwrap_or(Element::Unknown);

    let chain_key = match chain_id_raw {
        "." | "?" => "?".to_string(),
        other => other.to_string(),
    };
    if !chain_map.contains_key(&chain_key) {
        chain_map.insert(chain_key.clone(), BTreeMap::new());
        chain_order.push(chain_key.clone());
    }

    let residues = chain_map.get_mut(&chain_key).unwrap();
    let res_key = ResKey {
        chain_id: chain_key.clone(),
        res_seq,
        i_code,
    };

    let temp_res = residues.entry(res_key).or_insert_with(|| TempResidue {
        raw_name: res_name.to_string(),
        is_hetatm,
        atoms: HashMap::new(),
    });

    let atom_key = atom_name.to_string();
    let candidate = Atom::new(atom_name, element, pos);

    match temp_res.atoms.get(&atom_key) {
        Some((old_occ, _)) if occupancy <= *old_occ => {}
        _ => {
            temp_res.atoms.insert(atom_key, (occupancy, candidate));
        }
    }

    Ok(())
}

fn build_structure(
    mut structure: Structure,
    chain_order: Vec<String>,
    mut chain_map: HashMap<String, BTreeMap<ResKey, TempResidue>>,
    context: &IoContext,
) -> Result<Structure, Error> {
    for chain_id in chain_order {
        if let Some(residues) = chain_map.remove(&chain_id) {
            let mut chain = Chain::new(&chain_id);

            for (res_key, temp_res) in residues {
                let (canonical_name, std_enum) = context.classify_residue(&temp_res.raw_name);

                let category = determine_category(&temp_res, std_enum)?;

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

fn determine_category(
    temp_res: &TempResidue,
    std_enum: Option<StandardResidue>,
) -> Result<ResidueCategory, Error> {
    if std_enum.is_some() {
        Ok(ResidueCategory::Standard)
    } else if !temp_res.is_hetatm {
        Err(Error::unknown_standard_residue(&temp_res.raw_name, None))
    } else if temp_res.atoms.len() == 1 {
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

            let is_protein = matches!(residue.standard_name, Some(std) if std.is_protein());
            let is_nucleic = matches!(residue.standard_name, Some(std) if std.is_nucleic());

            if i == first_idx {
                residue.position = if is_protein {
                    ResiduePosition::NTerminal
                } else if is_nucleic {
                    ResiduePosition::FivePrime
                } else {
                    ResiduePosition::Internal
                };
            } else if i == last_idx {
                residue.position = if is_protein {
                    ResiduePosition::CTerminal
                } else if is_nucleic {
                    ResiduePosition::ThreePrime
                } else {
                    ResiduePosition::Internal
                };
            } else if polymer_indices.contains(&i) {
                residue.position = ResiduePosition::Internal;
            } else {
                residue.position = ResiduePosition::None;
            }
        }
    }
}

fn process_cell_parameters(params: &HashMap<String, String>) -> Option<[[f64; 3]; 3]> {
    let get_f64 = |key: &str, def: f64| -> f64 {
        params
            .get(key)
            .and_then(|s| f64::from_str(s).ok())
            .unwrap_or(def)
    };

    let a = get_f64("_cell.length_a", 0.0);
    let b = get_f64("_cell.length_b", 0.0);
    let c = get_f64("_cell.length_c", 0.0);

    if a == 0.0 || b == 0.0 || c == 0.0 {
        return None;
    }

    let alpha = get_f64("_cell.angle_alpha", 90.0).to_radians();
    let beta = get_f64("_cell.angle_beta", 90.0).to_radians();
    let gamma = get_f64("_cell.angle_gamma", 90.0).to_radians();

    let cos_a = alpha.cos();
    let cos_b = beta.cos();
    let cos_g = gamma.cos();
    let sin_g = gamma.sin();

    let v1_x = a;
    let v2_x = b * cos_g;
    let v2_y = b * sin_g;

    let term = (cos_a - cos_b * cos_g) / sin_g;
    let v3_x = c * cos_b;
    let v3_y = c * term;
    let v3_z = c * (1.0 - cos_b * cos_b - term * term).sqrt();

    Some([[v1_x, 0.0, 0.0], [v2_x, v2_y, 0.0], [v3_x, v3_y, v3_z]])
}
