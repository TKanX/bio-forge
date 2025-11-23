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

    if is_protein_c_term {
        valid_names.insert("OXT".to_string());
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
