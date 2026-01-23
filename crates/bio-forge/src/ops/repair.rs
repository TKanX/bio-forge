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

/// Alignment pair mapping residue position to template position.
type AlignmentPairs = Vec<(Point, Point)>;

/// Missing atom data: (name, element, template_position).
type MissingAtoms = Vec<(String, Element, Point)>;

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
) -> (AlignmentPairs, MissingAtoms) {
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

    let u = svd.u.ok_or_else(|| {
        Error::alignment_failed("", 0, "SVD decomposition failed: U matrix unavailable")
    })?;
    let v_t = svd.v_t.ok_or_else(|| {
        Error::alignment_failed("", 0, "SVD decomposition failed: V^T matrix unavailable")
    })?;

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{
        atom::Atom,
        chain::Chain,
        residue::Residue,
        types::{Element, Point, ResidueCategory, ResiduePosition, StandardResidue},
    };

    fn add_atom_from_template(
        residue: &mut Residue,
        template: db::TemplateView<'_>,
        atom_name: &str,
    ) {
        let (_, element, pos) = template
            .heavy_atoms()
            .find(|(name, _, _)| *name == atom_name)
            .unwrap_or_else(|| panic!("template atom {atom_name} missing"));
        residue.add_atom(Atom::new(atom_name, element, pos));
    }

    fn add_hydrogen_from_template(
        residue: &mut Residue,
        template: db::TemplateView<'_>,
        atom_name: &str,
    ) {
        let (_, pos, _) = template
            .hydrogens()
            .find(|(name, _, _)| *name == atom_name)
            .unwrap_or_else(|| panic!("template hydrogen {atom_name} missing"));
        residue.add_atom(Atom::new(atom_name, Element::H, pos));
    }

    fn standard_residue(name: &str, id: i32, std: StandardResidue) -> Residue {
        Residue::new(id, None, name, Some(std), ResidueCategory::Standard)
    }

    fn distance(a: Point, b: Point) -> f64 {
        (a - b).norm()
    }

    fn angle_deg(a: Point, center: Point, b: Point) -> f64 {
        let v1 = (a - center).normalize();
        let v2 = (b - center).normalize();
        v1.dot(&v2).clamp(-1.0, 1.0).acos().to_degrees()
    }

    #[test]
    fn repair_residue_rebuilds_missing_heavy_atoms_and_cleans_extras() {
        let template = db::get_template("ALA").expect("template ALA");
        let mut residue = standard_residue("ALA", 1, StandardResidue::ALA);
        residue.position = ResiduePosition::Internal;

        add_atom_from_template(&mut residue, template, "N");
        add_atom_from_template(&mut residue, template, "CA");
        add_hydrogen_from_template(&mut residue, template, "HA");
        residue.add_atom(Atom::new("FAKE", Element::C, Point::new(5.0, 5.0, 5.0)));

        repair_residue(&mut residue).expect("repair succeeds");

        for (name, _, _) in template.heavy_atoms() {
            assert!(residue.has_atom(name), "missing heavy atom {name}");
        }
        assert!(
            residue.has_atom("HA"),
            "valid hydrogen removed unexpectedly"
        );
        assert!(
            !residue.has_atom("FAKE"),
            "extraneous atom should be removed"
        );
    }

    #[test]
    fn repair_residue_adds_oxt_for_cterm_protein() {
        let template = db::get_template("ALA").expect("template ALA");
        let mut residue = standard_residue("ALA", 10, StandardResidue::ALA);
        residue.position = ResiduePosition::CTerminal;

        add_atom_from_template(&mut residue, template, "C");
        add_atom_from_template(&mut residue, template, "CA");
        add_atom_from_template(&mut residue, template, "O");

        repair_residue(&mut residue).expect("repair succeeds");

        let oxt = residue.atom("OXT").expect("OXT should be synthesized");
        assert_eq!(oxt.element, Element::O);
    }

    #[test]
    fn repair_residue_adds_op3_for_5prime_with_phosphate() {
        let template = db::get_template("DA").expect("template DA");
        let mut residue = standard_residue("DA", 1, StandardResidue::DA);
        residue.position = ResiduePosition::FivePrime;

        add_atom_from_template(&mut residue, template, "P");
        add_atom_from_template(&mut residue, template, "OP1");
        add_atom_from_template(&mut residue, template, "OP2");
        add_atom_from_template(&mut residue, template, "O5'");
        add_atom_from_template(&mut residue, template, "C5'");
        add_atom_from_template(&mut residue, template, "C4'");

        repair_residue(&mut residue).expect("repair succeeds");

        assert!(residue.has_atom("P"), "phosphorus should be retained");
        assert!(residue.has_atom("OP1"), "OP1 should be retained");
        assert!(residue.has_atom("OP2"), "OP2 should be retained");
        assert!(residue.has_atom("O5'"), "O5' should be retained");
        let op3 = residue
            .atom("OP3")
            .expect("OP3 should be synthesized for 5'-phosphate");
        assert_eq!(op3.element, Element::O);
    }

    #[test]
    fn repair_residue_excludes_phosphate_for_5prime_without_p() {
        let template = db::get_template("DA").expect("template DA");
        let mut residue = standard_residue("DA", 1, StandardResidue::DA);
        residue.position = ResiduePosition::FivePrime;

        add_atom_from_template(&mut residue, template, "O5'");
        add_atom_from_template(&mut residue, template, "C5'");
        add_atom_from_template(&mut residue, template, "C4'");
        add_atom_from_template(&mut residue, template, "C3'");

        repair_residue(&mut residue).expect("repair succeeds");

        assert!(!residue.has_atom("P"), "P should not be synthesized");
        assert!(!residue.has_atom("OP1"), "OP1 should not be synthesized");
        assert!(!residue.has_atom("OP2"), "OP2 should not be synthesized");
        assert!(!residue.has_atom("OP3"), "OP3 should not be synthesized");
        assert!(residue.has_atom("O5'"), "O5' should be retained");
    }

    #[test]
    fn repair_residue_3prime_nucleic_preserves_o3() {
        let template = db::get_template("DA").expect("template DA");
        let mut residue = standard_residue("DA", 10, StandardResidue::DA);
        residue.position = ResiduePosition::ThreePrime;

        add_atom_from_template(&mut residue, template, "C3'");
        add_atom_from_template(&mut residue, template, "O3'");
        add_atom_from_template(&mut residue, template, "C4'");
        add_atom_from_template(&mut residue, template, "C5'");

        repair_residue(&mut residue).expect("repair succeeds");

        assert!(
            residue.has_atom("O3'"),
            "O3' should be present for 3' terminal"
        );
    }

    #[test]
    fn repair_residue_errors_when_no_alignment_atoms_survive() {
        let mut residue = standard_residue("ALA", 2, StandardResidue::ALA);
        residue.add_atom(Atom::new("FAKE", Element::C, Point::origin()));

        let err = repair_residue(&mut residue).expect_err("should fail without anchors");
        match err {
            Error::AlignmentFailed { .. } => {}
            other => panic!("unexpected error: {other:?}"),
        }
    }

    #[test]
    fn repair_structure_updates_standard_residues_only() {
        let template = db::get_template("GLY").expect("template GLY");
        let mut standard = standard_residue("GLY", 5, StandardResidue::GLY);
        add_atom_from_template(&mut standard, template, "N");
        add_atom_from_template(&mut standard, template, "CA");

        let mut hetero = Residue::new(20, None, "LIG", None, ResidueCategory::Hetero);
        hetero.add_atom(Atom::new("XX", Element::C, Point::new(-1.0, 0.0, 0.0)));

        let mut chain = Chain::new("A");
        chain.add_residue(standard);
        chain.add_residue(hetero);

        let mut structure = Structure::new();
        structure.add_chain(chain);

        repair_structure(&mut structure).expect("repair succeeds");

        let chain = structure.chain("A").expect("chain A");
        let fixed = chain.residue(5, None).unwrap();
        for (name, _, _) in template.heavy_atoms() {
            assert!(fixed.has_atom(name), "missing atom {name} after repair");
        }

        let hetero_after = chain.residue(20, None).unwrap();
        assert!(
            hetero_after.has_atom("XX"),
            "hetero residue should remain untouched"
        );
    }

    #[test]
    fn oxt_geometry_has_correct_distance_and_angle() {
        let template = db::get_template("ALA").expect("template ALA");
        let mut residue = standard_residue("ALA", 1, StandardResidue::ALA);
        residue.position = ResiduePosition::CTerminal;

        for (name, element, pos) in template.heavy_atoms() {
            residue.add_atom(Atom::new(name, element, pos));
        }

        repair_residue(&mut residue).expect("repair succeeds");

        let c = residue.atom("C").expect("C").pos;
        let o = residue.atom("O").expect("O").pos;
        let oxt = residue.atom("OXT").expect("OXT").pos;

        let c_oxt_dist = distance(c, oxt);
        assert!(
            (c_oxt_dist - CARBOXYL_CO_BOND_LENGTH).abs() < 0.1,
            "C-OXT distance {c_oxt_dist:.3} should be ~{CARBOXYL_CO_BOND_LENGTH} Å"
        );

        let o_c_oxt_angle = angle_deg(o, c, oxt);
        assert!(
            (o_c_oxt_angle - CARBOXYL_OCO_ANGLE_DEG).abs() < 3.0,
            "O-C-OXT angle {o_c_oxt_angle:.1}° should be ~{CARBOXYL_OCO_ANGLE_DEG}°"
        );
    }

    #[test]
    fn op3_geometry_has_correct_distance_and_tetrahedral_angles() {
        let template = db::get_template("DA").expect("template DA");
        let mut residue = standard_residue("DA", 1, StandardResidue::DA);
        residue.position = ResiduePosition::FivePrime;

        for (name, element, pos) in template.heavy_atoms() {
            residue.add_atom(Atom::new(name, element, pos));
        }

        repair_residue(&mut residue).expect("repair succeeds");

        let p = residue.atom("P").expect("P").pos;
        let op1 = residue.atom("OP1").expect("OP1").pos;
        let op2 = residue.atom("OP2").expect("OP2").pos;
        let o5 = residue.atom("O5'").expect("O5'").pos;
        let op3 = residue.atom("OP3").expect("OP3").pos;

        let p_op3_dist = distance(p, op3);
        assert!(
            (p_op3_dist - PHOSPHATE_PO_BOND_LENGTH).abs() < 0.1,
            "P-OP3 distance {p_op3_dist:.3} should be ~{PHOSPHATE_PO_BOND_LENGTH} Å"
        );

        let tetrahedral_angle = 109.5;
        let tolerance = 10.0;

        let op1_p_op3 = angle_deg(op1, p, op3);
        let op2_p_op3 = angle_deg(op2, p, op3);
        let o5_p_op3 = angle_deg(o5, p, op3);

        assert!(
            (op1_p_op3 - tetrahedral_angle).abs() < tolerance,
            "OP1-P-OP3 angle {op1_p_op3:.1}° should be ~{tetrahedral_angle}°"
        );
        assert!(
            (op2_p_op3 - tetrahedral_angle).abs() < tolerance,
            "OP2-P-OP3 angle {op2_p_op3:.1}° should be ~{tetrahedral_angle}°"
        );
        assert!(
            (o5_p_op3 - tetrahedral_angle).abs() < tolerance,
            "O5'-P-OP3 angle {o5_p_op3:.1}° should be ~{tetrahedral_angle}°"
        );
    }
}
