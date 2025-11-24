use crate::db;
use crate::model::{
    atom::Atom,
    chain::Chain,
    residue::Residue,
    structure::Structure,
    types::{Element, Point, ResidueCategory, StandardResidue},
};
use crate::ops::error::Error;
use nalgebra::{Rotation3, Vector3};
use rand::rngs::StdRng;
use rand::seq::{IndexedRandom, SliceRandom};
use rand::{Rng, SeedableRng};
use std::collections::HashMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Cation {
    Na,
    K,
    Mg,
    Ca,
    Li,
    Zn,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Anion {
    Cl,
    Br,
    I,
    F,
}

#[derive(Debug, Clone)]
pub struct SolvateConfig {
    pub margin: f64,
    pub water_spacing: f64,
    pub vdw_cutoff: f64,
    pub remove_existing: bool,
    pub cations: Vec<Cation>,
    pub anions: Vec<Anion>,
    pub target_charge: i32,
    pub rng_seed: Option<u64>,
}

impl Default for SolvateConfig {
    fn default() -> Self {
        Self {
            margin: 10.0,
            water_spacing: 3.1,
            vdw_cutoff: 2.4,
            remove_existing: true,
            cations: vec![Cation::Na],
            anions: vec![Anion::Cl],
            target_charge: 0,
            rng_seed: None,
        }
    }
}

impl Cation {
    pub fn element(&self) -> Element {
        match self {
            Cation::Na => Element::Na,
            Cation::K => Element::K,
            Cation::Mg => Element::Mg,
            Cation::Ca => Element::Ca,
            Cation::Li => Element::Li,
            Cation::Zn => Element::Zn,
        }
    }

    pub fn charge(&self) -> i32 {
        match self {
            Cation::Na | Cation::K | Cation::Li => 1,
            Cation::Mg | Cation::Ca | Cation::Zn => 2,
        }
    }

    pub fn name(&self) -> &'static str {
        match self {
            Cation::Na => "NA",
            Cation::K => "K",
            Cation::Mg => "MG",
            Cation::Ca => "CA",
            Cation::Li => "LI",
            Cation::Zn => "ZN",
        }
    }
}

impl Anion {
    pub fn element(&self) -> Element {
        match self {
            Anion::Cl => Element::Cl,
            Anion::Br => Element::Br,
            Anion::I => Element::I,
            Anion::F => Element::F,
        }
    }

    pub fn charge(&self) -> i32 {
        -1
    }

    pub fn name(&self) -> &'static str {
        match self {
            Anion::Cl => "CL",
            Anion::Br => "BR",
            Anion::I => "I",
            Anion::F => "F",
        }
    }
}

pub fn solvate_structure(structure: &mut Structure, config: &SolvateConfig) -> Result<(), Error> {
    if config.remove_existing {
        structure.retain_residues(|_chain_id, res| {
            let is_water = res.standard_name == Some(StandardResidue::HOH);
            let is_ion = res.category == ResidueCategory::Ion;
            !is_water && !is_ion
        });
        structure.prune_empty_chains();
    }

    let solvent_chain_id = next_solvent_chain_id(structure);
    let mut rng = build_rng(config);

    let (min_bound, max_bound) = calculate_bounds(structure);
    let size = max_bound - min_bound;

    let box_dim = size
        + Vector3::new(
            config.margin * 2.0,
            config.margin * 2.0,
            config.margin * 2.0,
        );

    structure.box_vectors = Some([
        [box_dim.x, 0.0, 0.0],
        [0.0, box_dim.y, 0.0],
        [0.0, 0.0, box_dim.z],
    ]);

    let target_origin = Vector3::new(config.margin, config.margin, config.margin);
    let translation = target_origin - min_bound.coords;

    translate_structure(structure, &translation);

    let grid = SpatialGrid::new(structure, 4.0);

    let mut solvent_chain = Chain::new(&solvent_chain_id);
    let mut water_positions = Vec::new();

    let water_tmpl = db::get_template("HOH").ok_or(Error::MissingInternalTemplate {
        res_name: "HOH".to_string(),
    })?;

    let tmpl_o_pos = water_tmpl
        .heavy_atoms()
        .find(|(n, _, _)| *n == "O")
        .map(|(_, _, p)| p)
        .unwrap_or(Point::origin());

    let mut z = config.water_spacing / 2.0;
    let mut res_id_counter = 1;

    while z < box_dim.z {
        let mut y = config.water_spacing / 2.0;
        while y < box_dim.y {
            let mut x = config.water_spacing / 2.0;
            while x < box_dim.x {
                let candidate_pos = Point::new(x, y, z);

                if !grid.has_clash(&candidate_pos, config.vdw_cutoff) {
                    let rotation = Rotation3::from_axis_angle(
                        &Vector3::y_axis(),
                        rng.random_range(0.0..std::f64::consts::TAU),
                    ) * Rotation3::from_axis_angle(
                        &Vector3::x_axis(),
                        rng.random_range(0.0..std::f64::consts::TAU),
                    );

                    let mut residue = Residue::new(
                        res_id_counter,
                        None,
                        "HOH",
                        Some(StandardResidue::HOH),
                        ResidueCategory::Standard,
                    );

                    let final_o_pos = candidate_pos;
                    residue.add_atom(Atom::new("O", Element::O, final_o_pos));

                    for (h_name, h_pos, _) in water_tmpl.hydrogens() {
                        let local_vec = h_pos - tmpl_o_pos;
                        let rotated_vec = rotation * local_vec;
                        residue.add_atom(Atom::new(h_name, Element::H, final_o_pos + rotated_vec));
                    }

                    solvent_chain.add_residue(residue);
                    water_positions.push(res_id_counter);
                    res_id_counter += 1;
                }

                x += config.water_spacing;
            }
            y += config.water_spacing;
        }
        z += config.water_spacing;
    }

    replace_with_ions(
        structure,
        &mut solvent_chain,
        &mut water_positions,
        config,
        &mut rng,
    )?;

    if !solvent_chain.is_empty() {
        structure.add_chain(solvent_chain);
    }

    Ok(())
}

fn calculate_bounds(structure: &Structure) -> (Point, Point) {
    let mut min = Point::new(f64::MAX, f64::MAX, f64::MAX);
    let mut max = Point::new(f64::MIN, f64::MIN, f64::MIN);
    let mut count = 0;

    for atom in structure.iter_atoms() {
        min.x = min.x.min(atom.pos.x);
        min.y = min.y.min(atom.pos.y);
        min.z = min.z.min(atom.pos.z);
        max.x = max.x.max(atom.pos.x);
        max.y = max.y.max(atom.pos.y);
        max.z = max.z.max(atom.pos.z);
        count += 1;
    }

    if count == 0 {
        return (Point::origin(), Point::origin());
    }

    (min, max)
}

fn translate_structure(structure: &mut Structure, vec: &Vector3<f64>) {
    for atom in structure.iter_atoms_mut() {
        atom.translate_by(vec);
    }
}

fn calculate_solute_charge(structure: &Structure) -> i32 {
    let mut charge = 0;
    for chain in structure.iter_chains() {
        for residue in chain.iter_residues() {
            if let Some(tmpl) = db::get_template(&residue.name) {
                charge += tmpl.charge();
            } else if residue.category == ResidueCategory::Ion {
                match residue.name.as_str() {
                    "NA" | "K" | "LI" => charge += 1,
                    "MG" | "CA" | "ZN" => charge += 2,
                    "CL" | "BR" | "I" | "F" => charge -= 1,
                    _ => {}
                }
            }
        }
    }
    charge
}

fn replace_with_ions(
    structure: &Structure,
    solvent_chain: &mut Chain,
    water_indices: &mut Vec<i32>,
    config: &SolvateConfig,
    rng: &mut impl Rng,
) -> Result<(), Error> {
    if config.cations.is_empty() && config.anions.is_empty() {
        return Ok(());
    }

    let current_charge = calculate_solute_charge(structure);
    let mut charge_diff = config.target_charge - current_charge;

    water_indices.shuffle(rng);

    let mut attempts = 0;
    let max_attempts = water_indices.len();

    while charge_diff != 0 && attempts < max_attempts {
        if let Some(res_id) = water_indices.pop() {
            let residue = solvent_chain.residue_mut(res_id, None).unwrap();
            let pos = residue.atom("O").unwrap().pos;

            if charge_diff < 0 {
                if let Some(anion) = config.anions.choose(rng) {
                    *residue = create_anion_residue(res_id, *anion, pos);
                    charge_diff -= anion.charge();
                } else {
                    break;
                }
            } else if let Some(cation) = config.cations.choose(rng) {
                *residue = create_cation_residue(res_id, *cation, pos);
                charge_diff -= cation.charge();
            } else {
                break;
            }
        }
        attempts += 1;
    }

    if charge_diff != 0 {
        if water_indices.is_empty() {
            return Err(Error::BoxTooSmall);
        }

        return Err(Error::IonizationFailed {
            details: format!(
                "Could not reach target charge {}. Remaining diff: {}. Check if proper ion types are provided.",
                config.target_charge, charge_diff
            ),
        });
    }

    Ok(())
}

fn create_cation_residue(id: i32, cation: Cation, pos: Point) -> Residue {
    let mut res = Residue::new(id, None, cation.name(), None, ResidueCategory::Ion);
    res.add_atom(Atom::new(cation.name(), cation.element(), pos));
    res
}

fn create_anion_residue(id: i32, anion: Anion, pos: Point) -> Residue {
    let mut res = Residue::new(id, None, anion.name(), None, ResidueCategory::Ion);
    res.add_atom(Atom::new(anion.name(), anion.element(), pos));
    res
}

fn build_rng(config: &SolvateConfig) -> StdRng {
    if let Some(seed) = config.rng_seed {
        StdRng::seed_from_u64(seed)
    } else {
        StdRng::from_os_rng()
    }
}

fn next_solvent_chain_id(structure: &Structure) -> String {
    const BASE_ID: &str = "W";
    if structure.chain(BASE_ID).is_none() {
        return BASE_ID.to_string();
    }

    let mut index = 1;
    loop {
        let candidate = format!("{}{}", BASE_ID, index);
        if structure.chain(&candidate).is_none() {
            return candidate;
        }
        index += 1;
    }
}

struct SpatialGrid {
    cell_size: f64,
    cells: HashMap<(isize, isize, isize), Vec<Point>>,
}

impl SpatialGrid {
    fn new(structure: &Structure, cell_size: f64) -> Self {
        let mut cells: HashMap<(isize, isize, isize), Vec<Point>> = HashMap::new();

        for atom in structure.iter_atoms() {
            if atom.element == Element::H {
                continue;
            }

            let idx = Self::get_index(atom.pos, cell_size);
            cells.entry(idx).or_default().push(atom.pos);
        }

        Self { cell_size, cells }
    }

    fn get_index(pos: Point, size: f64) -> (isize, isize, isize) {
        (
            (pos.x / size).floor() as isize,
            (pos.y / size).floor() as isize,
            (pos.z / size).floor() as isize,
        )
    }

    fn has_clash(&self, pos: &Point, cutoff: f64) -> bool {
        let center_idx = Self::get_index(*pos, self.cell_size);
        let cutoff_sq = cutoff * cutoff;

        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let idx = (center_idx.0 + dx, center_idx.1 + dy, center_idx.2 + dz);
                    if let Some(atoms) = self.cells.get(&idx) {
                        for atom_pos in atoms {
                            if nalgebra::distance_squared(pos, atom_pos) < cutoff_sq {
                                return true;
                            }
                        }
                    }
                }
            }
        }
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{
        atom::Atom,
        chain::Chain,
        residue::Residue,
        structure::Structure,
        types::{Element, Point, ResidueCategory, StandardResidue},
    };

    #[test]
    fn removes_existing_solvent_and_repositions_solute() {
        let mut structure = Structure::new();

        let mut chain_a = Chain::new("A");
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        residue.add_atom(Atom::new("CA", Element::C, Point::new(1.0, 2.0, 3.0)));
        residue.add_atom(Atom::new("CB", Element::C, Point::new(3.0, 4.0, 5.0)));
        chain_a.add_residue(residue);
        structure.add_chain(chain_a);

        let mut solvent_chain = Chain::new("W");
        let mut existing_water = Residue::new(
            999,
            None,
            "HOH",
            Some(StandardResidue::HOH),
            ResidueCategory::Standard,
        );
        existing_water.add_atom(Atom::new("O", Element::O, Point::new(20.0, 20.0, 20.0)));
        solvent_chain.add_residue(existing_water);
        structure.add_chain(solvent_chain);

        let mut ion_chain = Chain::new("I");
        let mut ion = Residue::new(1000, None, "NA", None, ResidueCategory::Ion);
        ion.add_atom(Atom::new("NA", Element::Na, Point::new(25.0, 25.0, 25.0)));
        ion_chain.add_residue(ion);
        structure.add_chain(ion_chain);

        let config = SolvateConfig {
            margin: 5.0,
            water_spacing: 6.0,
            vdw_cutoff: 1.5,
            remove_existing: true,
            cations: vec![],
            anions: vec![],
            target_charge: 0,
            rng_seed: Some(42),
        };

        solvate_structure(&mut structure, &config).expect("solvation should succeed");

        let solute_chain = structure.chain("A").expect("solute chain");
        let mut min_coords = (f64::MAX, f64::MAX, f64::MAX);
        for atom in solute_chain.iter_atoms() {
            min_coords.0 = min_coords.0.min(atom.pos.x);
            min_coords.1 = min_coords.1.min(atom.pos.y);
            min_coords.2 = min_coords.2.min(atom.pos.z);
        }

        assert!((min_coords.0 - config.margin).abs() < 1e-6);
        assert!((min_coords.1 - config.margin).abs() < 1e-6);
        assert!((min_coords.2 - config.margin).abs() < 1e-6);

        let box_vectors = structure.box_vectors.expect("box vectors");
        assert!((box_vectors[0][0] - 12.0).abs() < 1e-6);
        assert!((box_vectors[1][1] - 12.0).abs() < 1e-6);
        assert!((box_vectors[2][2] - 12.0).abs() < 1e-6);

        let has_legacy_ids = structure
            .iter_chains()
            .flat_map(|chain| chain.iter_residues())
            .any(|res| res.id == 999 || res.id == 1000);
        assert!(!has_legacy_ids);

        let solvent_residues: Vec<_> = structure
            .iter_chains()
            .filter(|chain| chain.id.starts_with('W'))
            .flat_map(|chain| chain.iter_residues())
            .filter(|res| res.standard_name == Some(StandardResidue::HOH))
            .collect();
        assert!(!solvent_residues.is_empty());
    }

    #[test]
    fn populates_expected_number_of_waters_for_uniform_grid() {
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");
        let mut residue = Residue::new(
            1,
            None,
            "GLY",
            Some(StandardResidue::GLY),
            ResidueCategory::Standard,
        );
        residue.add_atom(Atom::new("CA", Element::C, Point::origin()));
        chain.add_residue(residue);
        structure.add_chain(chain);

        let config = SolvateConfig {
            margin: 4.0,
            water_spacing: 4.0,
            vdw_cutoff: 1.0,
            remove_existing: true,
            cations: vec![],
            anions: vec![],
            target_charge: 0,
            rng_seed: Some(7),
        };

        solvate_structure(&mut structure, &config).expect("solvation should succeed");

        let water_count = structure
            .iter_chains()
            .flat_map(|chain| chain.iter_residues())
            .filter(|res| res.standard_name == Some(StandardResidue::HOH))
            .count();

        assert_eq!(water_count, 8);
    }

    #[test]
    fn replaces_waters_with_anions_to_match_target_charge() {
        let lys_charge = db::get_template("LYS").expect("LYS template").charge();
        assert!(
            lys_charge > 0,
            "Test expects positively charged LYS template"
        );

        let mut structure = Structure::new();
        let mut chain = Chain::new("A");
        let mut residue = Residue::new(
            1,
            None,
            "LYS",
            Some(StandardResidue::LYS),
            ResidueCategory::Standard,
        );
        residue.add_atom(Atom::new("NZ", Element::N, Point::origin()));
        chain.add_residue(residue);
        structure.add_chain(chain);

        let config = SolvateConfig {
            margin: 4.0,
            water_spacing: 4.0,
            vdw_cutoff: 1.0,
            remove_existing: true,
            cations: vec![],
            anions: vec![Anion::Cl],
            target_charge: 0,
            rng_seed: Some(17),
        };

        solvate_structure(&mut structure, &config).expect("solvation should succeed");

        let ion_residues: Vec<_> = structure
            .iter_chains()
            .flat_map(|chain| chain.iter_residues())
            .filter(|res| res.category == ResidueCategory::Ion)
            .collect();

        assert_eq!(ion_residues.len() as i32, lys_charge);
        assert!(ion_residues.iter().all(|res| res.name == "CL"));
    }

    #[test]
    fn returns_box_too_small_when_insufficient_waters_for_target_charge() {
        let gly_charge = db::get_template("GLY").expect("GLY template").charge();
        assert_eq!(gly_charge, 0, "GLY should be neutral for this test");

        let mut structure = Structure::new();
        let mut chain = Chain::new("A");
        let mut residue = Residue::new(
            1,
            None,
            "GLY",
            Some(StandardResidue::GLY),
            ResidueCategory::Standard,
        );
        residue.add_atom(Atom::new("CA", Element::C, Point::origin()));
        chain.add_residue(residue);
        structure.add_chain(chain);

        let config = SolvateConfig {
            margin: 2.0,
            water_spacing: 7.0,
            vdw_cutoff: 0.1,
            remove_existing: true,
            cations: vec![Cation::Na],
            anions: vec![],
            target_charge: 2,
            rng_seed: Some(5),
        };

        let result = solvate_structure(&mut structure, &config);
        assert!(matches!(result, Err(Error::BoxTooSmall)));
    }
}
