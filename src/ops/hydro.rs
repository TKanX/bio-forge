use crate::db;
use crate::model::{
    atom::Atom,
    residue::Residue,
    structure::Structure,
    types::{Element, Point, ResidueCategory, ResiduePosition, StandardResidue},
};
use crate::ops::error::Error;
use nalgebra::{Matrix3, Rotation3, Vector3};
use rand::Rng;
use std::collections::HashSet;

const DISULFIDE_SG_THRESHOLD: f64 = 2.2;
const N_TERM_PKA: f64 = 8.0;
const C_TERM_PKA: f64 = 3.1;

#[derive(Debug, Clone)]
pub struct HydroConfig {
    pub target_ph: Option<f64>,
    pub remove_existing_h: bool,
    pub his_strategy: HisStrategy,
}

impl Default for HydroConfig {
    fn default() -> Self {
        Self {
            target_ph: None,
            remove_existing_h: true,
            his_strategy: HisStrategy::HbNetwork,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HisStrategy {
    DirectHID,
    DirectHIE,
    Random,
    HbNetwork,
}

pub fn add_hydrogens(structure: &mut Structure, config: &HydroConfig) -> Result<(), Error> {
    let mut targets = Vec::new();
    for (c_idx, chain) in structure.iter_chains().enumerate() {
        for (r_idx, residue) in chain.iter_residues().enumerate() {
            if residue.category == ResidueCategory::Standard {
                targets.push((c_idx, r_idx));
            }
        }
    }

    mark_disulfide_bridges(structure);

    for (c_idx, r_idx) in targets {
        let new_name = determine_protonation_state(structure, c_idx, r_idx, config);

        let chain = structure.iter_chains_mut().nth(c_idx).unwrap();
        let residue = chain.iter_residues_mut().nth(r_idx).unwrap();

        if let Some(name) = new_name {
            residue.name = name;
        }

        if config.remove_existing_h {
            residue.strip_hydrogens();
        }

        construct_hydrogens_for_residue(residue, config)?;
    }

    Ok(())
}

fn determine_protonation_state(
    structure: &Structure,
    c_idx: usize,
    r_idx: usize,
    config: &HydroConfig,
) -> Option<String> {
    let chain = structure.iter_chains().nth(c_idx)?;
    let residue = chain.iter_residues().nth(r_idx)?;
    let std = residue.standard_name?;

    if std == StandardResidue::CYS && residue.name == "CYX" {
        return None;
    }

    if let Some(ph) = config.target_ph {
        return match std {
            StandardResidue::ASP => Some(if ph < 3.9 {
                "ASH".to_string()
            } else {
                "ASP".to_string()
            }),
            StandardResidue::GLU => Some(if ph < 4.2 {
                "GLH".to_string()
            } else {
                "GLU".to_string()
            }),
            StandardResidue::LYS => Some(if ph > 10.5 {
                "LYN".to_string()
            } else {
                "LYS".to_string()
            }),
            StandardResidue::ARG => Some(if ph > 12.5 {
                "ARN".to_string()
            } else {
                "ARG".to_string()
            }),
            StandardResidue::CYS => Some(if ph > 8.3 {
                "CYM".to_string()
            } else {
                "CYS".to_string()
            }),
            StandardResidue::TYR => Some(if ph > 10.0 {
                "TYM".to_string()
            } else {
                "TYR".to_string()
            }),
            StandardResidue::HIS => {
                if ph < 6.0 {
                    Some("HIP".to_string())
                } else {
                    Some(select_neutral_his(structure, residue, &config.his_strategy))
                }
            }
            _ => None,
        };
    }

    if std == StandardResidue::HIS {
        if matches!(residue.name.as_str(), "HIS" | "HID" | "HIE" | "HIP") {
            return Some(select_neutral_his(structure, residue, &config.his_strategy));
        }
    }

    None
}

fn mark_disulfide_bridges(structure: &mut Structure) {
    let mut cys_sulfurs = Vec::new();
    for (c_idx, chain) in structure.iter_chains().enumerate() {
        for (r_idx, residue) in chain.iter_residues().enumerate() {
            if residue
                .standard_name
                .map(|s| s == StandardResidue::CYS)
                .unwrap_or(false)
            {
                if let Some(sg) = residue.atom("SG") {
                    cys_sulfurs.push((c_idx, r_idx, sg.pos));
                }
            }
        }
    }

    if cys_sulfurs.len() < 2 {
        return;
    }

    let threshold_sq = DISULFIDE_SG_THRESHOLD * DISULFIDE_SG_THRESHOLD;
    let mut disulfide_residues: HashSet<(usize, usize)> = HashSet::new();
    for i in 0..cys_sulfurs.len() {
        for j in (i + 1)..cys_sulfurs.len() {
            let (ci, ri, pos_i) = &cys_sulfurs[i];
            let (cj, rj, pos_j) = &cys_sulfurs[j];
            if (*pos_i - *pos_j).norm_squared() <= threshold_sq {
                disulfide_residues.insert((*ci, *ri));
                disulfide_residues.insert((*cj, *rj));
            }
        }
    }

    if disulfide_residues.is_empty() {
        return;
    }

    for (c_idx, chain) in structure.iter_chains_mut().enumerate() {
        for (r_idx, residue) in chain.iter_residues_mut().enumerate() {
            if disulfide_residues.contains(&(c_idx, r_idx)) && residue.name != "CYX" {
                residue.name = "CYX".to_string();
            }
        }
    }
}

fn select_neutral_his(structure: &Structure, residue: &Residue, strategy: &HisStrategy) -> String {
    match strategy {
        HisStrategy::DirectHID => "HID".to_string(),
        HisStrategy::DirectHIE => "HIE".to_string(),
        HisStrategy::Random => {
            let mut rng = rand::rng();
            if rng.random_bool(0.5) {
                "HID".to_string()
            } else {
                "HIE".to_string()
            }
        }
        HisStrategy::HbNetwork => optimize_his_network(structure, residue),
    }
}

fn optimize_his_network(structure: &Structure, residue: &Residue) -> String {
    let nd1 = residue.atom("ND1");
    let ne2 = residue.atom("NE2");

    let has_acceptor_near = |atom: &Atom| -> bool {
        for chain in structure.iter_chains() {
            for other_res in chain.iter_residues() {
                if other_res.id == residue.id && chain.id == "Unknown" {
                    continue;
                }
                if other_res == residue {
                    continue;
                }

                for other_atom in other_res.iter_atoms() {
                    if matches!(other_atom.element, Element::O | Element::N) {
                        if atom.distance_squared(other_atom) < 3.5 * 3.5 {
                            return true;
                        }
                    }
                }
            }
        }
        false
    };

    let nd1_interaction = nd1.map(|a| has_acceptor_near(a)).unwrap_or(false);
    let ne2_interaction = ne2.map(|a| has_acceptor_near(a)).unwrap_or(false);

    match (nd1_interaction, ne2_interaction) {
        (true, false) => "HID".to_string(),
        (false, true) => "HIE".to_string(),
        _ => "HID".to_string(),
    }
}

fn construct_hydrogens_for_residue(
    residue: &mut Residue,
    config: &HydroConfig,
) -> Result<(), Error> {
    let template_name = residue.name.clone();

    let template_view =
        db::get_template(&template_name).ok_or_else(|| Error::MissingInternalTemplate {
            res_name: template_name.clone(),
        })?;

    let existing_atoms: HashSet<String> = residue.atoms().iter().map(|a| a.name.clone()).collect();

    for (h_name, h_tmpl_pos, anchors_iter) in template_view.hydrogens() {
        if existing_atoms.contains(h_name) {
            continue;
        }

        let anchors: Vec<&str> = anchors_iter.collect();
        if let Ok(pos) = reconstruct_geometry(residue, h_tmpl_pos, &anchors) {
            residue.add_atom(Atom::new(h_name, Element::H, pos));
        } else {
            return Err(Error::incomplete_for_hydro(
                &residue.name,
                residue.id,
                anchors.first().copied().unwrap_or("?"),
            ));
        }
    }

    match residue.position {
        ResiduePosition::NTerminal if residue.standard_name.map_or(false, |s| s.is_protein()) => {
            construct_n_term_hydrogens(residue, n_term_should_be_protonated(config))?;
        }
        ResiduePosition::CTerminal if residue.standard_name.map_or(false, |s| s.is_protein()) => {
            construct_c_term_hydrogen(residue, c_term_should_be_protonated(config))?;
        }
        ResiduePosition::ThreePrime if residue.standard_name.map_or(false, |s| s.is_nucleic()) => {
            construct_3_prime_hydrogen(residue)?;
        }
        ResiduePosition::FivePrime if residue.standard_name.map_or(false, |s| s.is_nucleic()) => {
            if !residue.has_atom("P") && residue.has_atom("O5'") {
                construct_5_prime_hydrogen(residue)?;
            }
        }
        _ => {}
    }

    Ok(())
}

fn n_term_should_be_protonated(config: &HydroConfig) -> bool {
    config.target_ph.map(|ph| ph < N_TERM_PKA).unwrap_or(true)
}

fn c_term_should_be_protonated(config: &HydroConfig) -> bool {
    config.target_ph.map(|ph| ph < C_TERM_PKA).unwrap_or(false)
}

fn reconstruct_geometry(
    residue: &Residue,
    target_tmpl_pos: Point,
    anchor_names: &[&str],
) -> Result<Point, ()> {
    let template_view = db::get_template(&residue.name).ok_or(())?;

    let mut residue_pts = Vec::new();
    let mut template_pts = Vec::new();

    for name in anchor_names {
        let r_atom = residue.atom(name).ok_or(())?;
        residue_pts.push(r_atom.pos);

        let t_pos = template_view
            .heavy_atoms()
            .find(|(n, _, _)| n == name)
            .map(|(_, _, p)| p)
            .ok_or(())?;
        template_pts.push(t_pos);
    }

    let (rot, trans) = calculate_transform(&residue_pts, &template_pts).ok_or(())?;

    Ok(rot * target_tmpl_pos + trans)
}

fn construct_n_term_hydrogens(residue: &mut Residue, protonated: bool) -> Result<(), Error> {
    residue.remove_atom("H");
    residue.remove_atom("H1");
    residue.remove_atom("H2");
    residue.remove_atom("H3");

    let n_pos = residue
        .atom("N")
        .ok_or_else(|| Error::incomplete_for_hydro(&residue.name, residue.id, "N"))?
        .pos;
    let ca_pos = residue
        .atom("CA")
        .ok_or_else(|| Error::incomplete_for_hydro(&residue.name, residue.id, "CA"))?
        .pos;

    let v_ca_n = (n_pos - ca_pos).normalize();
    let bond_len = 1.0;
    let angle = 109.5_f64.to_radians();
    let sin_a = angle.sin();
    let cos_a = angle.cos();

    let up = if v_ca_n.x.abs() < 0.9 {
        Vector3::x()
    } else {
        Vector3::y()
    };
    let v_perp = v_ca_n.cross(&up).normalize();
    let base_dir = v_ca_n.scale(cos_a) + v_perp.scale(sin_a);

    let phases = [0.0_f64, 120.0, 240.0];
    let mut candidates: Vec<Vector3<f64>> = phases
        .iter()
        .map(|deg| {
            let rot_axis = Rotation3::from_axis_angle(
                &nalgebra::Unit::new_normalize(v_ca_n),
                deg.to_radians(),
            );
            rot_axis * base_dir
        })
        .collect();

    candidates.sort_by(|a, b| {
        a.dot(&v_ca_n)
            .partial_cmp(&b.dot(&v_ca_n))
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let target_count = if protonated { 3 } else { 2 };
    let names = ["H1", "H2", "H3"];
    for (idx, dir) in candidates.into_iter().take(target_count).enumerate() {
        residue.add_atom(Atom::new(
            names[idx],
            Element::H,
            n_pos + dir.scale(bond_len),
        ));
    }

    Ok(())
}

fn construct_c_term_hydrogen(residue: &mut Residue, protonated: bool) -> Result<(), Error> {
    if !protonated {
        residue.remove_atom("HOXT");
        residue.remove_atom("HXT");
        return Ok(());
    }

    residue.remove_atom("HXT");
    if residue.has_atom("HOXT") {
        return Ok(());
    }

    let c_pos = residue
        .atom("C")
        .ok_or_else(|| Error::incomplete_for_hydro(&residue.name, residue.id, "C"))?
        .pos;
    let oxt_pos = residue
        .atom("OXT")
        .ok_or_else(|| Error::incomplete_for_hydro(&residue.name, residue.id, "OXT"))?
        .pos;

    let direction = oxt_pos - c_pos;
    if direction.norm_squared() < 1e-6 {
        return Err(Error::incomplete_for_hydro(
            &residue.name,
            residue.id,
            "OXT",
        ));
    }
    let dir = direction.normalize();
    let h_pos = oxt_pos + dir.scale(0.97);

    residue.add_atom(Atom::new("HOXT", Element::H, h_pos));
    Ok(())
}

fn construct_3_prime_hydrogen(residue: &mut Residue) -> Result<(), Error> {
    if residue.has_atom("HO3'") {
        return Ok(());
    }

    let o3 = residue
        .atom("O3'")
        .ok_or_else(|| Error::incomplete_for_hydro(&residue.name, residue.id, "O3'"))?
        .pos;
    let c3 = residue
        .atom("C3'")
        .ok_or_else(|| Error::incomplete_for_hydro(&residue.name, residue.id, "C3'"))?
        .pos;
    let c4 = residue
        .atom("C4'")
        .or_else(|| residue.atom("C2'"))
        .map(|a| a.pos)
        .unwrap_or_else(|| c3 + Vector3::x());

    let v_c3_o3 = (o3 - c3).normalize();

    let v_c4_c3 = (c3 - c4).normalize();
    let normal = v_c3_o3.cross(&v_c4_c3).normalize();

    let h_dir = (v_c3_o3 + normal).normalize();
    let h_pos = o3 + h_dir.scale(0.96);

    residue.add_atom(Atom::new("HO3'", Element::H, h_pos));
    Ok(())
}

fn construct_5_prime_hydrogen(residue: &mut Residue) -> Result<(), Error> {
    if residue.has_atom("HO5'") {
        return Ok(());
    }

    let o5 = residue
        .atom("O5'")
        .ok_or_else(|| Error::incomplete_for_hydro(&residue.name, residue.id, "O5'"))?
        .pos;
    let c5 = residue
        .atom("C5'")
        .ok_or_else(|| Error::incomplete_for_hydro(&residue.name, residue.id, "C5'"))?
        .pos;

    let v_c5_o5 = (o5 - c5).normalize();
    let aux = if v_c5_o5.z.abs() < 0.9 {
        Vector3::z()
    } else {
        Vector3::x()
    };
    let perp = v_c5_o5.cross(&aux).normalize();

    let h_dir = (v_c5_o5 + perp).normalize();
    let h_pos = o5 + h_dir.scale(0.96);

    residue.add_atom(Atom::new("HO5'", Element::H, h_pos));
    Ok(())
}

fn calculate_transform(r_pts: &[Point], t_pts: &[Point]) -> Option<(Matrix3<f64>, Vector3<f64>)> {
    let n = r_pts.len();
    if n != t_pts.len() || n == 0 {
        return None;
    }

    let r_center = r_pts.iter().map(|p| p.coords).sum::<Vector3<f64>>() / n as f64;
    let t_center = t_pts.iter().map(|p| p.coords).sum::<Vector3<f64>>() / n as f64;

    if n == 1 {
        return Some((Matrix3::identity(), r_center - t_center));
    }

    if n == 2 {
        let v_r = r_pts[1] - r_pts[0];
        let v_t = t_pts[1] - t_pts[0];
        let rot = Rotation3::rotation_between(&v_t, &v_r).unwrap_or_else(Rotation3::identity);
        let trans = r_center - rot * t_center;
        return Some((rot.into_inner(), trans));
    }

    let mut cov = Matrix3::zeros();
    for (p_r, p_t) in r_pts.iter().zip(t_pts.iter()) {
        let v_r = p_r.coords - r_center;
        let v_t = p_t.coords - t_center;
        cov += v_r * v_t.transpose();
    }

    let svd = cov.svd(true, true);
    let u = svd.u?;
    let v_t = svd.v_t?;

    let mut rot = u * v_t;
    if rot.determinant() < 0.0 {
        let mut corr = Matrix3::identity();
        corr[(2, 2)] = -1.0;
        rot = u * corr * v_t;
    }

    let trans = r_center - rot * t_center;
    Some((rot, trans))
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

    fn residue_from_template(name: &str, std: StandardResidue, id: i32) -> Residue {
        let template = db::get_template(name).unwrap_or_else(|| panic!("missing template {name}"));
        let mut residue = Residue::new(id, None, name, Some(std), ResidueCategory::Standard);
        residue.position = ResiduePosition::Internal;
        for (atom_name, element, pos) in template.heavy_atoms() {
            residue.add_atom(Atom::new(atom_name, element, pos));
        }
        residue
    }

    fn structure_with_residue(residue: Residue) -> Structure {
        let mut chain = Chain::new("A");
        chain.add_residue(residue);
        let mut structure = Structure::new();
        structure.add_chain(chain);
        structure
    }

    fn structure_with_residues(residues: Vec<Residue>) -> Structure {
        let mut chain = Chain::new("A");
        for residue in residues {
            chain.add_residue(residue);
        }
        let mut structure = Structure::new();
        structure.add_chain(chain);
        structure
    }

    fn n_terminal_residue(id: i32) -> Residue {
        let mut residue = residue_from_template("ALA", StandardResidue::ALA, id);
        residue.position = ResiduePosition::NTerminal;
        residue
    }

    fn c_terminal_residue(id: i32) -> Residue {
        let mut residue = residue_from_template("ALA", StandardResidue::ALA, id);
        residue.position = ResiduePosition::CTerminal;
        let c_pos = residue.atom("C").expect("C atom").pos;
        let o_pos = residue.atom("O").expect("O atom").pos;
        let offset = c_pos - o_pos;
        let oxt_pos = c_pos + offset;
        residue.add_atom(Atom::new("OXT", Element::O, oxt_pos));
        residue
    }

    #[test]
    fn titratable_templates_exist_in_database() {
        let expected = [
            "ASP", "ASH", "GLU", "GLH", "LYS", "LYN", "ARG", "ARN", "CYS", "CYM", "TYR", "TYM",
            "HID", "HIE", "HIP",
        ];

        for name in expected {
            assert!(
                db::get_template(name).is_some(),
                "template {name} should exist"
            );
        }
    }

    #[test]
    fn determine_protonation_state_tracks_pka_thresholds() {
        let structure =
            structure_with_residue(residue_from_template("ASP", StandardResidue::ASP, 1));
        let mut config = HydroConfig {
            target_ph: Some(2.5),
            ..HydroConfig::default()
        };
        assert_eq!(
            determine_protonation_state(&structure, 0, 0, &config),
            Some("ASH".to_string())
        );

        config.target_ph = Some(5.0);
        assert_eq!(
            determine_protonation_state(&structure, 0, 0, &config),
            Some("ASP".to_string())
        );

        let structure =
            structure_with_residue(residue_from_template("LYS", StandardResidue::LYS, 2));
        config.target_ph = Some(11.0);
        assert_eq!(
            determine_protonation_state(&structure, 0, 0, &config),
            Some("LYN".to_string())
        );

        config.target_ph = Some(7.0);
        assert_eq!(
            determine_protonation_state(&structure, 0, 0, &config),
            Some("LYS".to_string())
        );
    }

    #[test]
    fn determine_protonation_state_respects_his_strategy() {
        let mut residue = residue_from_template("HID", StandardResidue::HIS, 3);
        residue.name = "HIS".to_string();
        let structure = structure_with_residue(residue);

        let config = HydroConfig {
            target_ph: Some(7.0),
            his_strategy: HisStrategy::DirectHIE,
            ..HydroConfig::default()
        };

        assert_eq!(
            determine_protonation_state(&structure, 0, 0, &config),
            Some("HIE".to_string())
        );

        let mut acid_config = HydroConfig::default();
        acid_config.target_ph = Some(5.5);
        assert_eq!(
            determine_protonation_state(&structure, 0, 0, &acid_config),
            Some("HIP".to_string())
        );
    }

    #[test]
    fn n_terminal_defaults_to_protonated_without_ph() {
        let residue = n_terminal_residue(40);
        let mut structure = structure_with_residue(residue);

        add_hydrogens(&mut structure, &HydroConfig::default()).expect("hydrogenation succeeds");

        let residue = structure.find_residue("A", 40, None).unwrap();
        assert!(residue.has_atom("H1"));
        assert!(residue.has_atom("H2"));
        assert!(residue.has_atom("H3"));
    }

    #[test]
    fn n_terminal_deprotonates_above_pka() {
        let residue = n_terminal_residue(41);
        let mut structure = structure_with_residue(residue);
        let mut config = HydroConfig::default();
        config.target_ph = Some(9.0);

        add_hydrogens(&mut structure, &config).expect("hydrogenation succeeds");

        let residue = structure.find_residue("A", 41, None).unwrap();
        assert!(residue.has_atom("H1"));
        assert!(residue.has_atom("H2"));
        assert!(!residue.has_atom("H3"));
    }

    #[test]
    fn c_terminal_protonates_under_acidic_ph() {
        let residue = c_terminal_residue(50);
        let mut structure = structure_with_residue(residue);
        let mut config = HydroConfig::default();
        config.target_ph = Some(2.5);

        add_hydrogens(&mut structure, &config).expect("hydrogenation succeeds");

        let residue = structure.find_residue("A", 50, None).unwrap();
        assert!(residue.has_atom("HOXT"));
    }

    #[test]
    fn c_terminal_remains_deprotonated_at_physiological_ph() {
        let residue = c_terminal_residue(51);
        let mut structure = structure_with_residue(residue);

        add_hydrogens(&mut structure, &HydroConfig::default()).expect("hydrogenation succeeds");

        let residue = structure.find_residue("A", 51, None).unwrap();
        assert!(!residue.has_atom("HOXT"));
    }

    #[test]
    fn determine_protonation_state_skips_already_marked_cyx() {
        let mut residue = residue_from_template("CYS", StandardResidue::CYS, 25);
        residue.name = "CYX".to_string();
        let structure = structure_with_residue(residue);
        let mut config = HydroConfig::default();
        config.target_ph = Some(9.0);

        assert_eq!(determine_protonation_state(&structure, 0, 0, &config), None);
    }

    #[test]
    fn add_hydrogens_populates_internal_lysine_side_chain() {
        let mut residue = residue_from_template("LYS", StandardResidue::LYS, 10);
        residue.add_atom(Atom::new("FAKE", Element::H, Point::origin()));
        let mut structure = structure_with_residue(residue);

        add_hydrogens(&mut structure, &HydroConfig::default()).expect("hydrogenation succeeds");

        let residue = structure.find_residue("A", 10, None).unwrap();
        assert!(residue.has_atom("HZ1"));
        assert!(residue.has_atom("HZ2"));
        assert!(residue.has_atom("HZ3"));
        assert!(!residue.has_atom("FAKE"));
    }

    #[test]
    fn add_hydrogens_relabels_asp_under_acidic_ph() {
        let residue = residue_from_template("ASP", StandardResidue::ASP, 15);
        let mut structure = structure_with_residue(residue);
        let mut config = HydroConfig::default();
        config.target_ph = Some(2.0);

        add_hydrogens(&mut structure, &config).expect("hydrogenation succeeds");

        let residue = structure.find_residue("A", 15, None).unwrap();
        assert_eq!(residue.name, "ASH");
        assert!(residue.has_atom("HD2"));
    }

    #[test]
    fn construct_hydrogens_errors_when_anchor_missing() {
        let template = db::get_template("ALA").expect("template ALA");
        let mut residue = Residue::new(
            20,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        residue.position = ResiduePosition::Internal;

        let (name, element, pos) = template.heavy_atoms().next().unwrap();
        residue.add_atom(Atom::new(name, element, pos));

        let err = construct_hydrogens_for_residue(&mut residue, &HydroConfig::default())
            .expect_err("should fail");
        assert!(matches!(err, Error::IncompleteResidueForHydro { .. }));
    }

    #[test]
    fn close_cysteines_are_relabelled_to_cyx_and_skip_hydrogens() {
        let cys1 = residue_from_template("CYS", StandardResidue::CYS, 30);
        let mut cys2 = residue_from_template("CYS", StandardResidue::CYS, 31);

        let sg1 = cys1.atom("SG").expect("SG in cys1").pos;
        let sg2 = cys2.atom("SG").expect("SG in cys2").pos;
        let desired = sg1 + Vector3::new(0.5, 0.0, 0.0);
        let offset = desired - sg2;
        for atom in cys2.iter_atoms_mut() {
            atom.translate_by(&offset);
        }

        let mut structure = structure_with_residues(vec![cys1, cys2]);
        add_hydrogens(&mut structure, &HydroConfig::default()).expect("hydrogenation");

        let res1 = structure.find_residue("A", 30, None).unwrap();
        let res2 = structure.find_residue("A", 31, None).unwrap();
        assert_eq!(res1.name, "CYX");
        assert_eq!(res2.name, "CYX");
        assert!(!res1.has_atom("HG"));
        assert!(!res2.has_atom("HG"));
    }
}
