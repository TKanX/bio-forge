//! Reconstructs standard residues so they match reference templates before topology building.
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
    types::{Element, ResidueCategory, ResiduePosition},
};
use crate::ops::error::Error;
use nalgebra::{Matrix3, Point3, Rotation3, Vector3};
use std::collections::HashSet;

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
    let mut targets = Vec::new();
    for chain in structure.iter_chains() {
        for residue in chain.iter_residues() {
            if residue.category == ResidueCategory::Standard {
                targets.push((chain.id.clone(), residue.id, residue.insertion_code));
            }
        }
    }

    for (chain_id, res_id, ic) in targets {
        if let Some(residue) = structure.find_residue_mut(&chain_id, res_id, ic) {
            repair_residue(residue)?;
        }
    }

    Ok(())
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
/// [`Error::AlignmentFailed`] if no anchor atoms remain, or other errors bubbled from
/// [`calculate_transform`].
fn repair_residue(residue: &mut Residue) -> Result<(), Error> {
    let template_name = residue.name.clone();
    let template_view =
        db::get_template(&template_name).ok_or_else(|| Error::MissingInternalTemplate {
            res_name: template_name.clone(),
        })?;

    let mut valid_names: HashSet<String> = HashSet::new();

    for (name, _, _) in template_view.heavy_atoms() {
        valid_names.insert(name.to_string());
    }
    for (name, _, _) in template_view.hydrogens() {
        valid_names.insert(name.to_string());
    }

    let is_protein_c_term = residue.standard_name.is_some_and(|s| s.is_protein())
        && residue.position == ResiduePosition::CTerminal;

    let is_nucleic_5prime = residue.standard_name.is_some_and(|s| s.is_nucleic())
        && residue.position == ResiduePosition::FivePrime;
    let has_5prime_phosphate = is_nucleic_5prime && residue.has_atom("P");

    if is_protein_c_term {
        valid_names.insert("OXT".to_string());
    }

    if is_nucleic_5prime {
        if has_5prime_phosphate {
            valid_names.insert("OP3".to_string());
        } else {
            valid_names.remove("P");
            valid_names.remove("OP1");
            valid_names.remove("OP2");
        }
    }

    let atoms_to_remove: Vec<String> = residue
        .atoms()
        .iter()
        .filter(|a| !valid_names.contains(&a.name))
        .map(|a| a.name.clone())
        .collect();

    for name in atoms_to_remove {
        residue.remove_atom(&name);
    }

    let mut align_pairs = Vec::new();
    let mut missing_heavy_atoms = Vec::new();

    for (name, element, tmpl_pos) in template_view.heavy_atoms() {
        if is_nucleic_5prime && !has_5prime_phosphate && matches!(name, "P" | "OP1" | "OP2") {
            continue;
        }

        if let Some(atom) = residue.atom(name) {
            align_pairs.push((atom.pos, tmpl_pos));
        } else {
            missing_heavy_atoms.push((name.to_string(), element, tmpl_pos));
        }
    }

    if is_protein_c_term {
        let tmpl_oxt_pos = calculate_template_oxt(template_view);

        if let Some(oxt) = residue.atom("OXT") {
            align_pairs.push((oxt.pos, tmpl_oxt_pos));
        } else {
            missing_heavy_atoms.push(("OXT".to_string(), Element::O, tmpl_oxt_pos));
        }
    }

    if has_5prime_phosphate {
        let tmpl_op3_pos = calculate_template_op3(template_view);

        if let Some(op3) = residue.atom("OP3") {
            align_pairs.push((op3.pos, tmpl_op3_pos));
        } else {
            missing_heavy_atoms.push(("OP3".to_string(), Element::O, tmpl_op3_pos));
        }
    }

    if align_pairs.is_empty() {
        return Err(Error::alignment_failed(
            &residue.name,
            residue.id,
            "No matching heavy atoms found for alignment",
        ));
    }

    let (rotation, translation) = calculate_transform(&align_pairs)?;

    for (name, element, tmpl_pos) in missing_heavy_atoms {
        let new_pos = rotation * tmpl_pos + translation;
        residue.add_atom(Atom::new(&name, element, new_pos));
    }

    Ok(())
}

/// Synthesizes an `OXT` position for C-terminal proteins using backbone vectors.
///
/// # Arguments
///
/// * `view` - Template view providing heavy atom coordinates.
///
/// # Returns
///
/// Estimated `OXT` coordinate; falls back to the origin if required atoms are absent.
fn calculate_template_oxt(view: db::TemplateView) -> Point3<f64> {
    let get_pos = |n| {
        view.heavy_atoms()
            .find(|(name, _, _)| *name == n)
            .map(|(_, _, pos)| pos)
    };

    let p_c = get_pos("C");
    let p_ca = get_pos("CA");
    let p_o = get_pos("O");

    if let (Some(c), Some(ca), Some(o)) = (p_c, p_ca, p_o) {
        let v_c_o = (o - c).normalize();
        let v_c_ca = (ca - c).normalize();

        let dir_oxt = -(v_c_o + v_c_ca).normalize();

        return c + dir_oxt * 1.25;
    }

    Point3::origin()
}

/// Synthesizes an `OP3` position for 5'-terminal phosphorylated nucleic acids.
///
/// The OP3 oxygen completes the tetrahedral coordination around phosphorus at 5' termini,
/// positioned opposite to the centroid of OP1, OP2, and O5' relative to P.
///
/// # Arguments
///
/// * `view` - Template view providing heavy atom coordinates.
///
/// # Returns
///
/// Estimated `OP3` coordinate; falls back to the origin if required atoms are absent.
fn calculate_template_op3(view: db::TemplateView) -> Point3<f64> {
    let get_pos = |n| {
        view.heavy_atoms()
            .find(|(name, _, _)| *name == n)
            .map(|(_, _, pos)| pos)
    };

    let p_p = get_pos("P");
    let p_op1 = get_pos("OP1");
    let p_op2 = get_pos("OP2");
    let p_o5 = get_pos("O5'");

    if let (Some(p), Some(op1), Some(op2), Some(o5)) = (p_p, p_op1, p_op2, p_o5) {
        let centroid = (op1.coords + op2.coords + o5.coords) / 3.0;
        let direction = (p.coords - centroid).normalize();

        return p + direction * 1.48;
    }

    Point3::origin()
}

/// Computes the best-fit rigid transform mapping template positions to residue coordinates.
///
/// Handles one- and two-point special cases before applying Kabsch alignment for larger sets.
///
/// # Arguments
///
/// * `pairs` - Corresponding residue/template coordinate pairs.
///
/// # Returns
///
/// Transformation `(rotation, translation)` aligning the template to the residue.
///
/// # Errors
///
/// Returns [`Error::AlignmentFailed`] when the SVD cannot produce a valid rotation.
fn calculate_transform(
    pairs: &[(Point3<f64>, Point3<f64>)],
) -> Result<(Matrix3<f64>, Vector3<f64>), Error> {
    let n = pairs.len();

    let center_res = pairs.iter().map(|p| p.0.coords).sum::<Vector3<f64>>() / n as f64;
    let center_tmpl = pairs.iter().map(|p| p.1.coords).sum::<Vector3<f64>>() / n as f64;

    if n == 1 {
        let translation = center_res - center_tmpl;
        return Ok((Matrix3::identity(), translation));
    }

    if n == 2 {
        let v_res = pairs[1].0 - pairs[0].0;
        let v_tmpl = pairs[1].1 - pairs[0].1;

        let rotation =
            Rotation3::rotation_between(&v_tmpl, &v_res).unwrap_or_else(Rotation3::identity);

        let translation = center_res - rotation * center_tmpl;

        return Ok((rotation.into_inner(), translation));
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
        .ok_or_else(|| Error::alignment_failed("N/A", 0, "SVD U failed"))?;
    let v_t = svd
        .v_t
        .ok_or_else(|| Error::alignment_failed("N/A", 0, "SVD V_T failed"))?;

    let mut rotation = u * v_t;

    if rotation.determinant() < 0.0 {
        let mut correction = Matrix3::identity();
        correction[(2, 2)] = -1.0;
        rotation = u * correction * v_t;
    }

    let translation = center_res - rotation * center_tmpl;

    Ok((rotation, translation))
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
}
