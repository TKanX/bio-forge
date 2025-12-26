//! mmCIF structure reader that reconstructs full chains, residues, and metadata.
//!
//! The parser tokenizes loop-based `_atom_site` tables, applies context-specific residue
//! aliasing, honors alternate locations using occupancy weights, and emits [`Structure`]
//! instances with categorized residues, terminal annotations, and optional unit-cell
//! vectors derived from `_cell.*` entries.

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

/// Composite key that uniquely identifies residues in a chain during parsing.
///
/// The key bundles author-provided chain IDs, sequence numbers, and optional insertion codes
/// so buffered atoms remain deterministically ordered until they are converted into
/// [`Residue`] objects.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct ResKey {
    chain_id: String,
    res_seq: i32,
    i_code: Option<char>,
}

/// Temporary aggregation of atoms prior to residue canonicalization.
///
/// Each entry tracks the raw residue name, whether it originated from `HETATM`, and the
/// best-occupancy atom positions keyed by atom name.
struct TempResidue {
    raw_name: String,
    is_hetatm: bool,
    atoms: HashMap<String, (f64, Atom)>,
}

/// Column index bookkeeping for `_atom_site` loop headers.
///
/// Fields store the zero-based column number for each required or optional header so that
/// row parsing can remain allocation-free once the header is mapped.
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

/// DFA states for the mmCIF tokenizer.
///
/// Tracks whether the parser is currently outside loops, consuming headers, or streaming
/// `_atom_site` rows versus unrelated loops so the correct handlers are invoked.
enum ParserState {
    Base,
    InLoopHeader,
    InAtomSiteLoop,
    InOtherLoop,
}

/// Parses mmCIF text into a normalized [`Structure`] using the supplied context.
///
/// The reader walks through `_cell.*` scalar entries, tokenizes `_atom_site` loops with
/// quoted-field awareness, collapses alternate locations by occupancy, and applies
/// [`IoContext`] aliasing plus template lookups to classify residues.
///
/// # Arguments
///
/// * `reader` - Any buffered reader that yields mmCIF text.
/// * `context` - Alias and template tables that normalize residue names and metadata.
///
/// # Returns
///
/// A fully populated [`Structure`] with chain ordering, residue categories, positions, and
/// optional unit-cell vectors derived from `_cell.*` parameters.
///
/// # Errors
///
/// Returns [`Error`] when encountering malformed loop headers, truncated `_atom_site`
/// records, unknown standard residues, or IO failures reported by `reader`.
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
                } else if current_loop_headers
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

/// Splits an mmCIF line into tokens while respecting quoted/semicolon blocks.
///
/// Handles both single- and double-quoted strings plus bare tokens separated by
/// whitespace so that `_atom_site` loops can be parsed without external crates.
///
/// # Arguments
///
/// * `line` - Raw line text pulled from the mmCIF stream.
///
/// # Returns
///
/// A vector of token strings preserving quoted content.
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

/// Maps `_atom_site` loop headers to column indices.
///
/// The function records the zero-based index of canonical headers so subsequent row
/// parsing can directly access coordinates, atom IDs, and optional attributes.
///
/// # Arguments
///
/// * `headers` - Ordered header strings collected immediately after `loop_`.
///
/// # Returns
///
/// A populated [`AtomSiteIndices`] structure with optional indices for each field.
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

/// Retrieves a mandatory token value by index with bounds checking.
///
/// # Arguments
///
/// * `tokens` - Slice of tokens for the current mmCIF row.
/// * `idx` - Zero-based index requested.
/// * `line_num` - Line number for contextual error reporting.
///
/// # Returns
///
/// The token as a `&str` or an [`Error::Parse`] if the index is missing.
fn token(tokens: &[String], idx: usize, line_num: usize) -> Result<&str, Error> {
    tokens.get(idx).map(|s| s.as_str()).ok_or_else(|| {
        Error::parse(
            "mmCIF",
            None,
            line_num,
            "Atom record shorter than _atom_site definition",
        )
    })
}

/// Retrieves an optional token, yielding `None` when the header is absent.
///
/// # Arguments
///
/// * `tokens` - Token slice for the processed row.
/// * `idx` - Optional column index for the field.
/// * `line_num` - Source line for diagnostics.
///
/// # Returns
///
/// `Ok(Some(&str))` when the index exists, `Ok(None)` if the header was omitted, or an
/// [`Error::Parse`] if the index exceeds the row length.
fn optional_token(
    tokens: &[String],
    idx: Option<usize>,
    line_num: usize,
) -> Result<Option<&str>, Error> {
    if let Some(idx) = idx {
        token(tokens, idx, line_num).map(Some)
    } else {
        Ok(None)
    }
}

/// Converts a coordinate token into an `f64`, enforcing proper numeric content.
///
/// # Arguments
///
/// * `value` - Text representation of the coordinate.
/// * `axis` - Axis name (X/Y/Z) for error messages.
/// * `line_num` - Source line used when building parse errors.
///
/// # Returns
///
/// Parsed coordinate as `f64`.
///
/// # Errors
///
/// [`Error::Parse`] when the token cannot be interpreted as a floating-point number.
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

/// Processes a single `_atom_site` row and appends it to the residue buffer.
///
/// Performs column lookups, handles missing values ("." and "?"), collapses chain IDs,
/// and stores highest-occupancy atoms grouped by residue and atom name.
///
/// # Arguments
///
/// * `tokens` - Tokenized row from the `_atom_site` loop.
/// * `indices` - Column indices resolved from the header.
/// * `line_num` - Source line for contextual errors.
/// * `chain_order` - Mutable list capturing encounter order of chains.
/// * `chain_map` - Aggregation of temporary residues keyed by [`ResKey`].
///
/// # Returns
///
/// [`Ok`] if the row was accepted or [`Error`] if required columns are missing/malformed.
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
    if required_indices
        .iter()
        .copied()
        .max()
        .is_some_and(|max_idx| tokens.len() <= max_idx)
    {
        return Err(Error::parse(
            "mmCIF",
            None,
            line_num,
            "Atom record is shorter than declared _atom_site headers",
        ));
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

/// Converts the buffered residue map into a [`Structure`].
///
/// Chains are created in first-seen order, residues are canonicalized through the
/// [`IoContext`], and polymer termini are annotated at the end.
///
/// # Arguments
///
/// * `structure` - Empty [`Structure`] to populate.
/// * `chain_order` - Encounter order recorded during parsing.
/// * `chain_map` - Buffered residues grouped by chain/key.
/// * `context` - IO context used for residue aliasing and templates.
///
/// # Returns
///
/// A fully populated [`Structure`] on success.
///
/// # Errors
///
/// Propagates [`Error`] from residue classification or IO context lookups.
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

/// Chooses a [`ResidueCategory`] from template matches or heuristics.
///
/// Standard residues are enforced via template lookup; unmapped standard entries produce an
/// error, while heterogens fall back to atom-count heuristics to distinguish ions.
///
/// # Arguments
///
/// * `temp_res` - Buffered residue metadata and atoms.
/// * `std_enum` - Optional standardized residue label returned by the context.
///
/// # Returns
///
/// [`ResidueCategory::Standard`] for templated residues, [`ResidueCategory::Ion`] for
/// single-atom heterogens, [`ResidueCategory::Hetero`] otherwise, or an error for unknown
/// standard residues.
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

/// Labels residues within each chain as terminal, internal, or unspecified.
///
/// Polymers are identified by standard templates (excluding water), with the first/last
/// residues marked as N/C or 5'/3' termini depending on polymer type.
///
/// # Arguments
///
/// * `structure` - Structure whose residues require positional annotations.
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

/// Converts `_cell.length_*` and `_cell.angle_*` parameters into box vectors.
///
/// # Arguments
///
/// * `params` - Key/value map collected while parsing cell metadata.
///
/// # Returns
///
/// `Some` orthogonalized 3×3 cell matrix if lengths are non-zero; otherwise `None`.
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::context::IoContext;
    use crate::model::types::{ResidueCategory, ResiduePosition, StandardResidue};
    use std::io::Cursor;

    const ATOM_SITE_HEADER: &str = "\
        loop_\n\
        _atom_site.group_PDB\n\
        _atom_site.auth_atom_id\n\
        _atom_site.auth_comp_id\n\
        _atom_site.auth_asym_id\n\
        _atom_site.auth_seq_id\n\
        _atom_site.pdbx_PDB_ins_code\n\
        _atom_site.Cartn_x\n\
        _atom_site.Cartn_y\n\
        _atom_site.Cartn_z\n\
        _atom_site.occupancy\n\
        _atom_site.type_symbol\n";

    fn parse_structure(cif: &str) -> Structure {
        let mut cursor = Cursor::new(cif.as_bytes());
        let context = IoContext::new_default();
        read(&mut cursor, &context).expect("mmCIF should parse")
    }

    fn parse_result(cif: &str) -> Result<Structure, Error> {
        let mut cursor = Cursor::new(cif.as_bytes());
        let context = IoContext::new_default();
        read(&mut cursor, &context)
    }

    #[test]
    fn read_parses_standard_polymer_and_box_vectors() {
        let cell_block = "\
            _cell.length_a 10.0\n\
            _cell.length_b 12.0\n\
            _cell.length_c 15.0\n\
            _cell.angle_alpha 90.0\n\
            _cell.angle_beta 90.0\n\
            _cell.angle_gamma 90.0\n";
        let rows = "\
            ATOM N ALA A 1 ? 12.546 11.406 2.324 1.00 N\n\
            ATOM CA ALA A 1 ? 13.123 12.345 3.210 1.00 C\n\
            ATOM N GLY A 2 ? 14.789 10.654 4.890 1.00 N\n\
            ATOM CA GLY A 2 ? 15.234 10.123 5.789 1.00 C\n";
        let cif = format!(
            "data_polymer\n{cell}\n{header}{rows}",
            cell = cell_block,
            header = ATOM_SITE_HEADER,
            rows = rows
        );

        let structure = parse_structure(&cif);
        assert_eq!(structure.chain_count(), 1);
        assert_eq!(structure.residue_count(), 2);

        let box_vectors = structure.box_vectors.expect("unit cell parsed");
        assert!((box_vectors[0][0] - 10.0).abs() < 1e-6);
        assert!((box_vectors[1][1] - 12.0).abs() < 1e-6);
        assert!((box_vectors[2][2] - 15.0).abs() < 1e-6);

        let chain = structure.chain("A").expect("chain A exists");
        let residues: Vec<_> = chain.iter_residues().collect();
        assert_eq!(residues.len(), 2);
        assert_eq!(residues[0].name, "ALA");
        assert_eq!(residues[0].position, ResiduePosition::NTerminal);
        assert_eq!(residues[1].name, "GLY");
        assert_eq!(residues[1].position, ResiduePosition::CTerminal);
    }

    #[test]
    fn read_aliases_water_and_applies_occupancy_filter() {
        let rows = "\
            HETATM O WAT B 5 ? 0.000 0.000 0.000 0.30 O\n\
            HETATM O WAT B 5 ? 1.000 1.000 1.000 0.80 O\n\
            HETATM NA NA B 6 ? 5.000 5.000 5.000 1.00 NA\n";
        let cif = format!(
            "data_water\n{header}{rows}",
            header = ATOM_SITE_HEADER,
            rows = rows
        );

        let structure = parse_structure(&cif);
        let chain = structure.chain("B").expect("chain B present");
        assert_eq!(chain.residue_count(), 2);

        let water = chain.residue(5, None).expect("water present");
        assert_eq!(water.name, "HOH");
        assert_eq!(water.category, ResidueCategory::Standard);
        let oxygen = water.atom("O").expect("oxygen retained");
        assert!((oxygen.pos.x - 1.0).abs() < 1e-6);

        let ion = chain.residue(6, None).expect("ion present");
        assert_eq!(ion.category, ResidueCategory::Ion);
        assert_eq!(ion.atom_count(), 1);
    }

    #[test]
    fn read_handles_scrambled_chain_and_residue_records() {
        let rows = "\
            ATOM CA GLY B 2 ? 14.000 10.000 9.000 1.00 C\n\
            ATOM N GLY B 1 ? 13.000 11.000 8.000 1.00 N\n\
            ATOM CA GLY B 1 ? 13.500 11.500 8.500 1.00 C\n\
            ATOM N ALA A 1 ? 11.000 12.000 5.000 1.00 N\n\
            ATOM CA ALA A 1 ? 11.500 12.500 5.500 1.00 C\n\
            ATOM N ALA A 2 ? 12.000 13.000 6.000 1.00 N\n\
            ATOM CA ALA A 2 ? 12.500 13.500 6.500 1.00 C\n\
            ATOM N GLY B 2 ? 13.500 10.500 8.500 1.00 N\n";
        let cif = format!(
            "data_order\n{header}{rows}",
            header = ATOM_SITE_HEADER,
            rows = rows
        );

        let structure = parse_structure(&cif);
        let chain_ids: Vec<_> = structure.iter_chains().map(|c| c.id.clone()).collect();
        assert_eq!(chain_ids, vec!["B", "A"]);

        let chain_b = structure.chain("B").unwrap();
        let resid_b: Vec<_> = chain_b.iter_residues().map(|r| r.id).collect();
        assert_eq!(resid_b, vec![1, 2]);

        let chain_a = structure.chain("A").unwrap();
        let resid_a: Vec<_> = chain_a.iter_residues().map(|r| r.id).collect();
        assert_eq!(resid_a, vec![1, 2]);
    }

    #[test]
    fn read_assigns_terminal_positions_for_polymers() {
        let rows = "\
            ATOM N GLY P 1 ? 1.000 2.000 3.000 1.00 N\n\
            ATOM CA GLY P 1 ? 1.500 2.500 3.500 1.00 C\n\
            HETATM O HOH P 2 ? 5.000 5.000 5.000 1.00 O\n\
            ATOM N SER P 3 ? 2.500 3.500 4.500 1.00 N\n\
            ATOM CA SER P 3 ? 3.000 4.000 5.000 1.00 C\n\
            ATOM N LEU P 4 ? 4.000 5.000 6.000 1.00 N\n\
            ATOM CA LEU P 4 ? 4.500 5.500 6.500 1.00 C\n\
            ATOM P DA N 1 ? 6.000 1.000 1.000 1.00 P\n\
            ATOM \"O5'\" DA N 1 ? 6.500 1.500 1.500 1.00 O\n\
            ATOM P DT N 2 ? 7.000 2.000 2.000 1.00 P\n\
            ATOM \"O5'\" DT N 2 ? 7.500 2.500 2.500 1.00 O\n";
        let cif = format!(
            "data_positions\n{header}{rows}",
            header = ATOM_SITE_HEADER,
            rows = rows
        );

        let structure = parse_structure(&cif);

        let protein = structure.chain("P").expect("protein chain present");
        let residues: Vec<_> = protein.iter_residues().collect();
        assert_eq!(residues[0].position, ResiduePosition::NTerminal);
        assert_eq!(residues[1].position, ResiduePosition::None);
        assert_eq!(residues[2].position, ResiduePosition::Internal);
        assert_eq!(residues[3].position, ResiduePosition::CTerminal);

        let nucleic = structure.chain("N").expect("nucleic chain present");
        let nuc_res: Vec<_> = nucleic.iter_residues().collect();
        assert_eq!(nuc_res[0].position, ResiduePosition::FivePrime);
        assert_eq!(nuc_res[1].position, ResiduePosition::ThreePrime);
    }

    #[test]
    fn read_categorizes_residues_based_on_templates_and_atom_counts() {
        let rows = "\
            ATOM N GLY C 1 ? 9.000 9.000 9.000 1.00 N\n\
            ATOM CA GLY C 1 ? 9.500 9.500 9.500 1.00 C\n\
            HETATM O HOH C 2 ? 5.000 5.000 5.000 1.00 O\n\
            HETATM C1 LIG C 301 ? 1.000 1.000 1.000 1.00 C\n\
            HETATM O1 LIG C 301 ? 1.500 1.500 1.500 1.00 O\n\
            HETATM NA NA C 401 ? 2.000 2.000 2.000 1.00 NA\n";
        let cif = format!(
            "data_categories\n{header}{rows}",
            header = ATOM_SITE_HEADER,
            rows = rows
        );

        let structure = parse_structure(&cif);
        let chain = structure.chain("C").unwrap();

        let protein = chain.residue(1, None).unwrap();
        assert_eq!(protein.category, ResidueCategory::Standard);

        let water = chain.residue(2, None).unwrap();
        assert_eq!(water.standard_name, Some(StandardResidue::HOH));
        assert_eq!(water.category, ResidueCategory::Standard);

        let ligand = chain.residue(301, None).unwrap();
        assert_eq!(ligand.category, ResidueCategory::Hetero);
        assert_eq!(ligand.atom_count(), 2);

        let ion = chain.residue(401, None).unwrap();
        assert_eq!(ion.category, ResidueCategory::Ion);
        assert_eq!(ion.atom_count(), 1);
    }

    #[test]
    fn read_prefers_atoms_with_highest_occupancy_and_retains_hydrogens() {
        let rows = "\
            ATOM N GLY D 1 ? 0.000 0.000 0.000 1.00 N\n\
            ATOM CA GLY D 1 ? 1.000 0.000 0.000 0.40 C\n\
            ATOM CA GLY D 1 ? 2.000 0.000 0.000 0.80 C\n\
            ATOM CA GLY D 1 ? 3.000 0.000 0.000 0.80 C\n\
            ATOM H1 GLY D 1 ? 0.500 0.500 0.500 1.00 H\n";
        let cif = format!(
            "data_occupancy\n{header}{rows}",
            header = ATOM_SITE_HEADER,
            rows = rows
        );

        let structure = parse_structure(&cif);
        let residue = structure.chain("D").unwrap().residue(1, None).unwrap();

        let ca = residue.atom("CA").unwrap();
        assert!((ca.pos.x - 2.0).abs() < 1e-6);
        assert_eq!(residue.atom_count(), 3);
        assert!(residue.atom("H1").is_some());
    }

    #[test]
    fn read_supports_residues_with_insertion_codes() {
        let rows = "\
            ATOM N SER E 10 ? 0.000 0.000 0.000 1.00 N\n\
            ATOM CA SER E 10 ? 0.500 0.500 0.500 1.00 C\n\
            ATOM N SER E 10 A 1.000 1.000 1.000 1.00 N\n\
            ATOM CA SER E 10 A 1.500 1.500 1.500 1.00 C\n\
            ATOM N SER E 11 ? 2.000 2.000 2.000 1.00 N\n\
            ATOM CA SER E 11 ? 2.500 2.500 2.500 1.00 C\n";
        let cif = format!(
            "data_insertions\n{header}{rows}",
            header = ATOM_SITE_HEADER,
            rows = rows
        );

        let structure = parse_structure(&cif);
        let chain = structure.chain("E").unwrap();
        let residues: Vec<_> = chain
            .iter_residues()
            .map(|r| (r.id, r.insertion_code))
            .collect();
        assert_eq!(residues, vec![(10, None), (10, Some('A')), (11, None)]);

        let base = chain.residue(10, None).unwrap();
        let insertion = chain.residue(10, Some('A')).unwrap();
        assert_ne!(
            base.atom("CA").unwrap().pos,
            insertion.atom("CA").unwrap().pos
        );
    }

    #[test]
    fn read_errors_on_unknown_standard_atom_record() {
        let rows = "\
            ATOM N UNK C 1 ? 0.000 0.000 0.000 1.00 N\n";
        let cif = format!(
            "data_error\n{header}{rows}",
            header = ATOM_SITE_HEADER,
            rows = rows
        );

        let err = parse_result(&cif).expect_err("unknown residue should fail");
        match err {
            Error::UnknownStandardResidue { name, path } => {
                assert_eq!(name, "UNK");
                assert!(path.is_none());
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }
}
