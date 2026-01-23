//! Reconstructs standard residues to match reference templates before topology building.
//!
//! The repair pipeline removes stray atoms, regenerates missing heavy atoms (including OXT on
//! C-terminal proteins and OP3 on 5'-phosphorylated nucleic acids), and aligns completions
//! using template geometry to ensure downstream hydrogenation and topology steps operate on
//! canonical coordinates.

use crate::db;
use crate::model::{
    atom::Atom,
    residue::Residue,
    structure::Structure,
    types::{Element, Point, ResidueCategory, ResiduePosition},
};
use crate::ops::error::Error;
use crate::utils::parallel::*;
use nalgebra::{Matrix3, Rotation3, Unit, Vector3};
use std::collections::HashSet;

/// Standard C-O bond length in carboxylate groups (Å).
const CARBOXYL_CO_BOND_LENGTH: f64 = 1.25;
/// O-C-OXT angle in carboxylate groups (degrees).
const CARBOXYL_OCO_ANGLE_DEG: f64 = 126.0;
/// Standard P-O bond length in phosphate groups (Å).
const PHOSPHATE_PO_BOND_LENGTH: f64 = 1.48;

/// Repairs every standard residue in a structure by invoking the internal repair logic.
///
/// Non-standard residues (heterogens, ions, solvent) are left untouched to avoid tampering
/// with ligands. Residue iteration happens first to avoid borrowing issues when mutating.
///
/// # Arguments
///
/// * `structure` - Mutable structure whose standard residues will be normalized.
///
/// # Returns
///
/// `Ok(())` when all residues are processed successfully.
///
/// # Errors
///
/// Propagates errors such as missing templates or alignment issues.
pub fn repair_structure(structure: &mut Structure) -> Result<(), Error> {
    structure
        .par_residues_mut()
        .filter(|r| r.category == ResidueCategory::Standard)
        .try_for_each(repair_residue)
}

/// Cleans and rebuilds an individual residue using its template definition.
///
/// Removes atoms absent from the template, calculates rigid alignment using shared anchors,
/// and adds back missing heavy atoms (including terminal `OXT` for proteins and `OP3` for
/// 5'-phosphorylated nucleic acids when applicable).
///
/// # Arguments
///
/// * `residue` - Residue to be normalized.
///
/// # Returns
///
/// `Ok(())` when the residue matches the template after repair.
///
/// # Errors
///
/// Returns [`Error::MissingInternalTemplate`] if the residue name lacks a template,
/// or [`Error::AlignmentFailed`] if no anchor atoms remain for alignment.
fn repair_residue(residue: &mut Residue) -> Result<(), Error> {
    let template_name = residue.name.clone();
    let template =
        db::get_template(&template_name).ok_or_else(|| Error::MissingInternalTemplate {
            res_name: template_name.to_string(),
        })?;

    let status = detect_terminal_status(residue);

    let valid_names = build_valid_names(template, &status);

    clean_invalid_atoms(residue, &valid_names);

    let (align_pairs, missing_atoms) = collect_alignment_data(residue, template, &status);

    if align_pairs.is_empty() {
        return Err(Error::alignment_failed(
            &*residue.name,
            residue.id,
            "No matching heavy atoms found for alignment",
        ));
    }

    let transform = calculate_transform(&align_pairs)?;
    synthesize_missing_template_atoms(residue, missing_atoms, &transform);

    synthesize_terminal_atoms(residue, &status);

    Ok(())
}

/// Encapsulates the terminal status of a residue for conditional processing.
struct TerminalStatus {
    /// True if this is a C-terminal protein residue.
    is_protein_c_term: bool,
    /// True if this is a 5'-terminal nucleic acid residue.
    is_nucleic_5prime: bool,
    /// True if the 5'-terminal has a phosphate group (P atom present).
    has_5prime_phosphate: bool,
}

/// Detects the terminal status of a residue based on its type and position.
///
/// # Arguments
///
/// * `residue` - Residue whose terminal status is to be determined.
///
/// # Returns
///
/// `TerminalStatus` struct indicating terminal properties.
fn detect_terminal_status(residue: &Residue) -> TerminalStatus {
    let is_protein = residue.standard_name.is_some_and(|s| s.is_protein());
    let is_nucleic = residue.standard_name.is_some_and(|s| s.is_nucleic());

    TerminalStatus {
        is_protein_c_term: is_protein && residue.position == ResiduePosition::CTerminal,
        is_nucleic_5prime: is_nucleic && residue.position == ResiduePosition::FivePrime,
        has_5prime_phosphate: is_nucleic
            && residue.position == ResiduePosition::FivePrime
            && residue.has_atom("P"),
    }
}

/// Builds the set of atom names that should exist in the final residue.
///
/// Includes all template atoms plus terminal-specific atoms, with adjustments
/// for 5'-terminal nucleic acids without phosphate groups.
///
/// # Arguments
///
/// * `template` - Template view for the residue.
/// * `status` - Terminal status of the residue.
///
/// # Returns
///
/// `HashSet<String>` containing valid atom names.
fn build_valid_names(template: db::TemplateView, status: &TerminalStatus) -> HashSet<String> {
    let mut names = HashSet::new();

    for (name, _, _) in template.heavy_atoms() {
        names.insert(name.to_string());
    }

    for (name, _, _) in template.hydrogens() {
        names.insert(name.to_string());
    }

    if status.is_protein_c_term {
        names.insert("OXT".to_string());
    }

    if status.is_nucleic_5prime {
        if status.has_5prime_phosphate {
            names.insert("OP3".to_string());
        } else {
            names.remove("P");
            names.remove("OP1");
            names.remove("OP2");
        }
    }

    names
}

/// Removes atoms from the residue that are not in the valid name set.
///
/// # Arguments
///
/// * `residue` - Residue to be cleaned.
/// * `valid_names` - Set of valid atom names to retain.
fn clean_invalid_atoms(residue: &mut Residue, valid_names: &HashSet<String>) {
    let atoms_to_remove: Vec<String> = residue
        .atoms()
        .iter()
        .filter(|a| !valid_names.contains(a.name.as_str()))
        .map(|a| a.name.to_string())
        .collect();

    for name in atoms_to_remove {
        residue.remove_atom(&name);
    }
}

/// Collects alignment pairs and identifies missing heavy atoms from the template.
///
/// # Arguments
///
/// * `residue` - Residue being repaired.
/// * `template` - Template view for the residue.
/// * `status` - Terminal status of the residue.
///
/// # Returns
///
/// Tuple containing a vector of alignment pairs and a vector of missing atoms.
fn collect_alignment_data(
    residue: &Residue,
    template: db::TemplateView,
    status: &TerminalStatus,
) -> (Vec<(Point, Point)>, Vec<(String, Element, Point)>) {
    let mut align_pairs = Vec::new();
    let mut missing_atoms = Vec::new();

    for (name, element, tmpl_pos) in template.heavy_atoms() {
        if status.is_nucleic_5prime
            && !status.has_5prime_phosphate
            && matches!(name, "P" | "OP1" | "OP2")
        {
            continue;
        }

        if let Some(atom) = residue.atom(name) {
            align_pairs.push((atom.pos, tmpl_pos));
        } else {
            missing_atoms.push((name.to_string(), element, tmpl_pos));
        }
    }

    (align_pairs, missing_atoms)
}

/// Rigid transformation consisting of rotation and translation.
struct Transform {
    rotation: Matrix3<f64>,
    translation: Vector3<f64>,
}

impl Transform {
    /// Applies the transformation to a point in template space.
    ///
    /// # Arguments
    ///
    /// * `point` - Point in template space to be transformed.
    ///
    /// # Returns
    ///
    /// Transformed `Point` in residue space.
    fn apply(&self, point: Point) -> Point {
        Point::from(self.rotation * point.coords + self.translation)
    }
}

/// Computes the best-fit rigid transform mapping template positions to residue coordinates.
///
/// Handles special cases:
/// - Single point: pure translation
/// - Two points: single-axis rotation + translation
/// - Three or more points: full Kabsch SVD alignment
///
/// # Arguments
///
/// * `pairs` - Slice of point pairs (residue position, template position).
///
/// # Returns
///
/// `Ok(Transform)` containing the computed rotation and translation.
///
/// # Errors
///
/// Returns [`Error::AlignmentFailed`] if SVD computation fails.
fn calculate_transform(pairs: &[(Point, Point)]) -> Result<Transform, Error> {
    let n = pairs.len();

    let center_res = pairs.iter().map(|p| p.0.coords).sum::<Vector3<f64>>() / n as f64;
    let center_tmpl = pairs.iter().map(|p| p.1.coords).sum::<Vector3<f64>>() / n as f64;

    if n == 1 {
        return Ok(Transform {
            rotation: Matrix3::identity(),
            translation: center_res - center_tmpl,
        });
    }

    if n == 2 {
        let v_res = pairs[1].0 - pairs[0].0;
        let v_tmpl = pairs[1].1 - pairs[0].1;

        let rotation =
            Rotation3::rotation_between(&v_tmpl, &v_res).unwrap_or_else(Rotation3::identity);

        return Ok(Transform {
            rotation: rotation.into_inner(),
            translation: center_res - rotation * center_tmpl,
        });
    }

    let mut cov = Matrix3::zeros();

    for (p_res, p_tmpl) in pairs {
        let v_res = p_res.coords - center_res;
        let v_tmpl = p_tmpl.coords - center_tmpl;

        cov += v_res * v_tmpl.transpose();
    }

    let svd = cov.svd(true, true);

    let u = svd
        .u
        .ok_or_else(|| Error::alignment_failed("N/A", 0, "SVD U matrix computation failed"))?;
    let v_t = svd
        .v_t
        .ok_or_else(|| Error::alignment_failed("N/A", 0, "SVD V_T matrix computation failed"))?;

    let mut rotation = u * v_t;
    if rotation.determinant() < 0.0 {
        let mut correction = Matrix3::identity();
        correction[(2, 2)] = -1.0;
        rotation = u * correction * v_t;
    }

    Ok(Transform {
        rotation,
        translation: center_res - rotation * center_tmpl,
    })
}

/// Synthesizes missing template atoms by applying the SVD transform.
///
/// # Arguments
///
/// * `residue` - Residue to which missing atoms will be added.
/// * `missing_atoms` - Vector of missing atom data (name, element, template position).
/// * `transform` - Rigid transform to map template positions to residue space.
fn synthesize_missing_template_atoms(
    residue: &mut Residue,
    missing_atoms: Vec<(String, Element, Point)>,
    transform: &Transform,
) {
    for (name, element, tmpl_pos) in missing_atoms {
        let new_pos = transform.apply(tmpl_pos);
        residue.add_atom(Atom::new(&name, element, new_pos));
    }
}

/// Dispatches terminal atom synthesis based on terminal status.
///
/// # Arguments
///
/// * `residue` - Residue to be modified.
/// * `status` - Terminal status of the residue.
fn synthesize_terminal_atoms(residue: &mut Residue, status: &TerminalStatus) {
    if status.is_protein_c_term {
        synthesize_oxt(residue);
    }

    if status.has_5prime_phosphate {
        synthesize_op3(residue);
    }
}

/// Synthesizes OXT atom for C-terminal protein residues using local carboxyl geometry.
///
/// The OXT position is calculated entirely in world-space using the actual positions
/// of C, O, and CA atoms, ensuring correct carboxylate geometry regardless of global
/// residue orientation.
///
/// # Arguments
///
/// * `residue` - Residue to which OXT will be added.
fn synthesize_oxt(residue: &mut Residue) {
    if residue.has_atom("OXT") {
        return;
    }

    let (c, o, ca) = match (residue.atom("C"), residue.atom("O"), residue.atom("CA")) {
        (Some(c), Some(o), Some(ca)) => (c.pos, o.pos, ca.pos),
        _ => return,
    };

    let v_co = (o - c).normalize();
    let v_cca = (ca - c).normalize();

    let normal = v_co.cross(&v_cca);
    if normal.norm() < 1e-10 {
        return;
    }
    let normal = Unit::new_normalize(normal);

    let angle = CARBOXYL_OCO_ANGLE_DEG.to_radians();
    let rotation = Rotation3::from_axis_angle(&normal, angle);
    let oxt_direction = rotation * v_co;

    let oxt_pos = c + oxt_direction * CARBOXYL_CO_BOND_LENGTH;

    residue.add_atom(Atom::new("OXT", Element::O, oxt_pos));
}

/// Synthesizes OP3 atom for 5'-terminal phosphorylated nucleic acid residues.
///
/// The OP3 position is calculated using tetrahedral geometry around the phosphorus
/// atom, placing OP3 opposite to the centroid of the other three oxygen ligands.
///
/// # Arguments
///
/// * `residue` - Residue to which OP3 will be added.
fn synthesize_op3(residue: &mut Residue) {
    if residue.has_atom("OP3") {
        return;
    }

    let (p, op1, op2, o5) = match (
        residue.atom("P"),
        residue.atom("OP1"),
        residue.atom("OP2"),
        residue.atom("O5'"),
    ) {
        (Some(p), Some(op1), Some(op2), Some(o5)) => (p.pos, op1.pos, op2.pos, o5.pos),
        _ => return,
    };

    let centroid = Point::from((op1.coords + op2.coords + o5.coords) / 3.0);

    let direction = (p - centroid).normalize();

    let op3_pos = p + direction * PHOSPHATE_PO_BOND_LENGTH;

    residue.add_atom(Atom::new("OP3", Element::O, op3_pos));
}
