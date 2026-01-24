//! Protonation-aware hydrogen placement for protein and nucleic acid structures.
//!
//! This module inspects residue templates, predicts protonation states from pH or hydrogen
//! bonding networks, relabels residues accordingly, and rebuilds missing hydrogens while
//! respecting polymer termini, nucleic acid priming, and disulfide bridges.

use crate::db;
use crate::model::{
    atom::Atom,
    grid::Grid,
    residue::Residue,
    structure::Structure,
    types::{Element, Point, ResidueCategory, ResiduePosition, StandardResidue},
};
use crate::ops::error::Error;
use crate::utils::parallel::*;
use nalgebra::{Matrix3, Rotation3, Vector3};
use rand::Rng;
use std::collections::HashSet;

/// Maximum sulfur–sulfur distance (Å) used to detect disulfide bridges.
const DISULFIDE_SG_THRESHOLD: f64 = 2.2;
/// Henderson–Hasselbalch breakpoint for protonated N-termini.
const N_TERM_PKA: f64 = 8.0;
/// Henderson–Hasselbalch breakpoint for protonated C-termini.
const C_TERM_PKA: f64 = 3.1;
/// Henderson–Hasselbalch breakpoint for the second dissociation of terminal phosphate.
const PHOSPHATE_PKA2: f64 = 6.5;
/// Standard sp³ tetrahedral bond angle (degrees).
const SP3_ANGLE: f64 = 109.5;
/// Standard N-H bond length (Å).
const NH_BOND_LENGTH: f64 = 1.01;
/// Standard O-H bond length (Å).
const OH_BOND_LENGTH: f64 = 0.96;
/// Carboxylic acid O-H bond length (Å).
const COOH_BOND_LENGTH: f64 = 0.97;

/// Parameters controlling hydrogen addition behavior.
///
/// `HydroConfig` can target a specific solution pH, remove pre-existing hydrogens, and
/// choose how neutral histidine tautomers are assigned.
#[derive(Debug, Clone)]
pub struct HydroConfig {
    /// Optional solvent pH value used for titration decisions.
    pub target_ph: Option<f64>,
    /// Whether to strip all existing hydrogens before reconstruction.
    pub remove_existing_h: bool,
    /// Strategy for setting neutral histidine tautomer labels.
    pub his_strategy: HisStrategy,
}

impl Default for HydroConfig {
    /// Provides biologically reasonable defaults (physiological pH, removal of old hydrogens,
    /// and hydrogen-bond-aware histidine selection).
    fn default() -> Self {
        Self {
            target_ph: None,
            remove_existing_h: true,
            his_strategy: HisStrategy::HbNetwork,
        }
    }
}

/// Strategies for choosing neutral histidine labels.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HisStrategy {
    /// Force all neutral histidines to HID (delta-protonated).
    DirectHID,
    /// Force all neutral histidines to HIE (epsilon-protonated).
    DirectHIE,
    /// Randomly choose between HID and HIE.
    Random,
    /// Analyze hydrogen-bond networks to select the tautomer most likely to bond.
    HbNetwork,
}

/// Adds hydrogens to all standard residues in-place, updating protonation states when needed.
///
/// Disulfide bridges are detected prior to titration and preserved by relabeling cysteines to
/// `CYX`. Residues lacking required anchor atoms will trigger descriptive errors so upstream
/// cleanup steps can correct the issues.
///
/// # Arguments
///
/// * `structure` - Mutable structure whose residues will be protonated and hydrated.
/// * `config` - Hydrogenation configuration, including pH target and histidine strategy.
///
/// # Returns
///
/// `Ok(())` when hydrogenation succeeds.
///
/// # Errors
///
/// Returns [`Error::MissingInternalTemplate`] when no template is found or
/// [`Error::IncompleteResidueForHydro`] when required anchor atoms are missing.
pub fn add_hydrogens(structure: &mut Structure, config: &HydroConfig) -> Result<(), Error> {
    mark_disulfide_bridges(structure);

    let acceptor_grid = if config.his_strategy == HisStrategy::HbNetwork {
        Some(build_acceptor_grid(structure))
    } else {
        None
    };

    structure
        .par_chains_mut()
        .enumerate()
        .try_for_each(|(c_idx, chain)| {
            chain
                .par_residues_mut()
                .enumerate()
                .try_for_each(|(r_idx, residue)| {
                    if residue.category != ResidueCategory::Standard {
                        return Ok(());
                    }

                    let new_name = determine_protonation_state(
                        residue,
                        config,
                        acceptor_grid.as_ref(),
                        Some((c_idx, r_idx)),
                    );

                    if let Some(name) = new_name {
                        residue.name = name.into();
                    }

                    if config.remove_existing_h {
                        residue.strip_hydrogens();
                    }

                    construct_hydrogens_for_residue(residue, config)
                })
        })
}

/// Builds a spatial grid of all nitrogen and oxygen atoms in the structure.
///
/// # Arguments
///
/// * `structure` - Structure from which to extract acceptor atoms.
///
/// # Returns
///
/// A `Grid` containing positions and residue indices of all N, O, and F atoms.
fn build_acceptor_grid(structure: &Structure) -> Grid<(usize, usize)> {
    let atoms: Vec<(Point, (usize, usize))> = structure
        .iter_chains()
        .enumerate()
        .flat_map(|(c_idx, chain)| {
            chain
                .iter_residues()
                .enumerate()
                .flat_map(move |(r_idx, residue)| {
                    residue
                        .iter_atoms()
                        .filter(|a| matches!(a.element, Element::N | Element::O | Element::F))
                        .map(move |a| (a.pos, (c_idx, r_idx)))
                })
        })
        .collect();
    Grid::new(atoms, 3.5)
}

/// Predicts the protonation-induced residue rename for a given polymer residue.
///
/// Applies residue-specific pKa thresholds, handles histidine tautomer selection, and
/// respects previously tagged disulfide cysteines.
///
/// # Arguments
///
/// * `residue` - Residue under evaluation.
/// * `config` - Hydrogenation configuration containing pH and histidine options.
/// * `grid` - Optional grid of acceptor atoms for HIS network analysis.
/// * `indices` - Optional (chain_idx, res_idx) to exclude self-interactions.
///
/// # Returns
///
/// `Some(new_name)` when the residue should be relabeled; otherwise `None`.
fn determine_protonation_state(
    residue: &Residue,
    config: &HydroConfig,
    grid: Option<&Grid<(usize, usize)>>,
    indices: Option<(usize, usize)>,
) -> Option<String> {
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
                    Some(select_neutral_his(
                        residue,
                        config.his_strategy,
                        grid,
                        indices,
                    ))
                }
            }
            _ => None,
        };
    }

    if std == StandardResidue::HIS && matches!(residue.name.as_str(), "HIS" | "HID" | "HIE" | "HIP")
    {
        return Some(select_neutral_his(
            residue,
            config.his_strategy,
            grid,
            indices,
        ));
    }

    None
}

/// Identifies cysteine pairs forming disulfide bonds and renames them to `CYX`.
///
/// # Arguments
///
/// * `structure` - Mutable structure containing residues to scan and relabel.
fn mark_disulfide_bridges(structure: &mut Structure) {
    let cys_sulfurs: Vec<(Point, (usize, usize))> = structure
        .par_chains_mut()
        .enumerate()
        .flat_map(|(c_idx, chain)| {
            chain
                .par_residues_mut()
                .enumerate()
                .filter_map(move |(r_idx, residue)| {
                    if matches!(residue.standard_name, Some(StandardResidue::CYS)) {
                        residue.atom("SG").map(|sg| (sg.pos, (c_idx, r_idx)))
                    } else {
                        None
                    }
                })
        })
        .collect();

    if cys_sulfurs.len() < 2 {
        return;
    }

    let grid = Grid::new(cys_sulfurs.clone(), DISULFIDE_SG_THRESHOLD + 0.5);

    let disulfide_residues: HashSet<(usize, usize)> = cys_sulfurs
        .par_iter()
        .flat_map_iter(|(pos, (c_idx, r_idx))| {
            grid.neighbors(pos, DISULFIDE_SG_THRESHOLD)
                .exact()
                .filter_map(move |(_, &(neighbor_c, neighbor_r))| {
                    if *c_idx == neighbor_c && *r_idx == neighbor_r {
                        None
                    } else {
                        Some([(*c_idx, *r_idx), (neighbor_c, neighbor_r)])
                    }
                })
        })
        .flatten()
        .collect();

    if disulfide_residues.is_empty() {
        return;
    }

    structure
        .par_chains_mut()
        .enumerate()
        .for_each(|(c_idx, chain)| {
            chain
                .par_residues_mut()
                .enumerate()
                .for_each(|(r_idx, residue)| {
                    if disulfide_residues.contains(&(c_idx, r_idx)) && residue.name != "CYX" {
                        residue.name = "CYX".into();
                    }
                });
        });
}

/// Selects a neutral histidine tautomer using the configured strategy.
///
/// # Arguments
///
/// * `residue` - Histidine residue being relabeled.
/// * `strategy` - Selection strategy from [`HisStrategy`].
/// * `grid` - Optional grid of acceptor atoms for HIS network analysis.
/// * `indices` - Optional (chain_idx, res_idx) to exclude self-interactions.
///
/// # Returns
///
/// Either `"HID"` or `"HIE"`.
fn select_neutral_his(
    residue: &Residue,
    strategy: HisStrategy,
    grid: Option<&Grid<(usize, usize)>>,
    indices: Option<(usize, usize)>,
) -> String {
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
        HisStrategy::HbNetwork => optimize_his_network(residue, grid, indices),
    }
}

/// Determines the best histidine tautomer by inspecting nearby hydrogen-bond acceptors.
///
/// # Arguments
///
/// * `residue` - Histidine residue being evaluated.
/// * `grid` - Grid of acceptor atoms (N/O/F) in the structure.
/// * `indices` - (chain_idx, res_idx) of the current residue.
///
/// # Returns
///
/// The histidine label (`"HID"` or `"HIE"`) that maximizes compatibility with neighbors.
fn optimize_his_network(
    residue: &Residue,
    grid: Option<&Grid<(usize, usize)>>,
    indices: Option<(usize, usize)>,
) -> String {
    let (grid, self_indices) = match (grid, indices) {
        (Some(g), Some(i)) => (g, i),
        _ => return "HIE".to_string(),
    };

    let score_hid = calculate_h_bond_score(residue, "ND1", "CG", "CE1", grid, self_indices);

    let score_hie = calculate_h_bond_score(residue, "NE2", "CD2", "CE1", grid, self_indices);

    if score_hid > score_hie {
        "HID".to_string()
    } else {
        "HIE".to_string()
    }
}

/// Calculates a geometric score for a potential hydrogen bond network.
///
/// Computes the hypothetical hydrogen position for an sp2 nitrogen and sums the
/// scores of valid H-bonds formed with nearby acceptors. A valid H-bond must
/// satisfy both distance (< 2.7Å H...A) and angle (> 90° N-H...A) constraints.
///
/// # Arguments
///
/// * `residue` - The residue containing the donor nitrogen.
/// * `n_name` - Name of the nitrogen atom (donor).
/// * `c1_name` - Name of the first neighbor carbon.
/// * `c2_name` - Name of the second neighbor carbon.
/// * `grid` - Spatial index of potential acceptors.
/// * `self_idx` - Chain and residue index of the current residue to exclude self-interactions.
///
/// # Returns
///
/// A positive floating-point score, where higher indicates better H-bonding.
fn calculate_h_bond_score(
    residue: &Residue,
    n_name: &str,
    c1_name: &str,
    c2_name: &str,
    grid: &Grid<(usize, usize)>,
    self_idx: (usize, usize),
) -> f64 {
    let n = match residue.atom(n_name) {
        Some(a) => a,
        None => return 0.0,
    };
    let c1 = match residue.atom(c1_name) {
        Some(a) => a,
        None => return 0.0,
    };
    let c2 = match residue.atom(c2_name) {
        Some(a) => a,
        None => return 0.0,
    };

    let v1 = (c1.pos - n.pos).normalize();
    let v2 = (c2.pos - n.pos).normalize();
    let bisector = (v1 + v2).normalize();
    let h_dir = -bisector;
    let h_pos = n.pos + h_dir;

    let mut score = 0.0;

    for (a_pos, &idx) in grid.neighbors(&n.pos, 3.5).exact() {
        if idx == self_idx {
            continue;
        }

        let h_a_vec = a_pos - h_pos;
        let dist_sq = h_a_vec.norm_squared();

        if dist_sq > 2.7 * 2.7 {
            continue;
        }

        let h_a_dir = h_a_vec.normalize();
        let cos_theta = h_dir.dot(&h_a_dir);

        if cos_theta > 0.0 {
            score += (1.0 / dist_sq) * (cos_theta * cos_theta);
        }
    }

    score
}

/// Rebuilds hydrogens for a single residue using template geometry and terminal rules.
///
/// # Arguments
///
/// * `residue` - Residue to augment with hydrogens.
/// * `config` - Hydrogenation configuration influencing terminal protonation.
///
/// # Returns
///
/// `Ok(())` when all hydrogens were added or already present.
///
/// # Errors
///
/// Returns [`Error::MissingInternalTemplate`] or [`Error::IncompleteResidueForHydro`] when
/// required template data or anchor atoms are missing.
fn construct_hydrogens_for_residue(
    residue: &mut Residue,
    config: &HydroConfig,
) -> Result<(), Error> {
    let template_name = residue.name.clone();

    let template_view =
        db::get_template(&template_name).ok_or_else(|| Error::MissingInternalTemplate {
            res_name: template_name.to_string(),
        })?;

    let existing_atoms: HashSet<String> =
        residue.atoms().iter().map(|a| a.name.to_string()).collect();

    let rotation_override = if residue.standard_name == Some(StandardResidue::HOH) {
        let mut rng = rand::rng();
        Some(
            Rotation3::from_axis_angle(
                &Vector3::y_axis(),
                rng.random_range(0.0..std::f64::consts::TAU),
            ) * Rotation3::from_axis_angle(
                &Vector3::x_axis(),
                rng.random_range(0.0..std::f64::consts::TAU),
            ),
        )
    } else {
        None
    };

    for (h_name, h_tmpl_pos, anchors_iter) in template_view.hydrogens() {
        if existing_atoms.contains(h_name) {
            continue;
        }

        let anchors: Vec<&str> = anchors_iter.collect();
        if let Ok(pos) = reconstruct_geometry(residue, h_tmpl_pos, &anchors, rotation_override) {
            residue.add_atom(Atom::new(h_name, Element::H, pos));
        } else {
            return Err(Error::incomplete_for_hydro(
                &*residue.name,
                residue.id,
                anchors.first().copied().unwrap_or("?"),
            ));
        }
    }

    match residue.position {
        ResiduePosition::NTerminal if residue.standard_name.is_some_and(|s| s.is_protein()) => {
            construct_n_term_hydrogens(residue, n_term_should_be_protonated(config))?;
        }
        ResiduePosition::CTerminal if residue.standard_name.is_some_and(|s| s.is_protein()) => {
            construct_c_term_hydrogen(residue, c_term_should_be_protonated(config))?;
        }
        ResiduePosition::ThreePrime if residue.standard_name.is_some_and(|s| s.is_nucleic()) => {
            construct_3_prime_hydrogen(residue)?;
        }
        ResiduePosition::FivePrime if residue.standard_name.is_some_and(|s| s.is_nucleic()) => {
            if residue.has_atom("P") {
                construct_5_prime_phosphate_hydrogens(residue, config)?;
            } else if residue.has_atom("O5'") {
                construct_5_prime_hydrogen(residue)?;
            }
        }
        _ => {}
    }

    Ok(())
}

/// Evaluates whether an N-terminus should remain protonated under the configured pH.
///
/// # Arguments
///
/// * `config` - Hydrogenation configuration that may provide a target pH.
///
/// # Returns
///
/// `true` when the pH is below the configured threshold or unspecified.
fn n_term_should_be_protonated(config: &HydroConfig) -> bool {
    config.target_ph.map(|ph| ph < N_TERM_PKA).unwrap_or(true)
}

/// Evaluates whether a C-terminus should remain protonated under the configured pH.
///
/// # Arguments
///
/// * `config` - Hydrogenation configuration potentially specifying pH.
///
/// # Returns
///
/// `true` when the pH is below the acidic cutoff; otherwise `false`.
fn c_term_should_be_protonated(config: &HydroConfig) -> bool {
    config.target_ph.map(|ph| ph < C_TERM_PKA).unwrap_or(false)
}

/// Aligns template coordinates to residue anchors to predict a hydrogen position.
///
/// # Arguments
///
/// * `residue` - Residue containing measured anchor atoms.
/// * `target_tmpl_pos` - Target hydrogen position from the template.
/// * `anchor_names` - Atom names used to determine the rigid-body transform.
/// * `rotation_override` - Optional rotation to use if the system is under-constrained (e.g. water).
///
/// # Returns
///
/// `Ok(point)` containing the placed coordinate or `Err(())` if anchors are missing.
fn reconstruct_geometry(
    residue: &Residue,
    target_tmpl_pos: Point,
    anchor_names: &[&str],
    rotation_override: Option<Rotation3<f64>>,
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

    let (mut rot, mut trans) = calculate_transform(&residue_pts, &template_pts).ok_or(())?;

    if let (Some(override_rot), 1) = (rotation_override, anchor_names.len()) {
        rot = override_rot.into_inner();
        trans = residue_pts[0].coords - rot * template_pts[0].coords;
    }

    Ok(rot * target_tmpl_pos + trans)
}

/// Rebuilds the N-terminal amine hydrogens using tetrahedral geometry.
///
/// # Arguments
///
/// * `residue` - Residue whose terminal hydrogens are being reconstructed.
/// * `protonated` - Whether to place three hydrogens (`true`) or two (`false`).
///
/// # Returns
///
/// `Ok(())` when hydrogens are successfully placed.
///
/// # Errors
///
/// Returns [`Error::IncompleteResidueForHydro`] if the N or CA atoms are missing.
fn construct_n_term_hydrogens(residue: &mut Residue, protonated: bool) -> Result<(), Error> {
    residue.remove_atom("H");
    residue.remove_atom("H1");
    residue.remove_atom("H2");
    residue.remove_atom("H3");

    let n_pos = residue
        .atom("N")
        .ok_or_else(|| Error::incomplete_for_hydro(&*residue.name, residue.id, "N"))?
        .pos;
    let ca_pos = residue
        .atom("CA")
        .ok_or_else(|| Error::incomplete_for_hydro(&*residue.name, residue.id, "CA"))?
        .pos;

    let (x, y, z) = build_sp3_frame(n_pos, ca_pos, None);

    let theta = SP3_ANGLE.to_radians();
    let sin_theta = theta.sin();
    let cos_theta = theta.cos();

    let phases = [0.0_f64, 120.0, 240.0];
    let target_count = if protonated { 3 } else { 2 };
    let names = ["H1", "H2", "H3"];

    for (idx, phase) in phases.iter().take(target_count).enumerate() {
        let phi = phase.to_radians();
        let h_local = Vector3::new(sin_theta * phi.cos(), sin_theta * phi.sin(), -cos_theta);
        let h_global = x * h_local.x + y * h_local.y + z * h_local.z;
        let h_pos = n_pos + h_global * NH_BOND_LENGTH;
        residue.add_atom(Atom::new(names[idx], Element::H, h_pos));
    }

    Ok(())
}

/// Rebuilds the carboxylate proton at the C-terminus when protonated.
///
/// # Arguments
///
/// * `residue` - Residue to modify.
/// * `protonated` - Whether the terminus should include the `HOXT` hydrogen.
///
/// # Returns
///
/// `Ok(())` after either removing or adding hydrogens as needed.
///
/// # Errors
///
/// Returns [`Error::IncompleteResidueForHydro`] if the `C` or `OXT` atoms are missing.
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
        .ok_or_else(|| Error::incomplete_for_hydro(&*residue.name, residue.id, "C"))?
        .pos;
    let oxt_pos = residue
        .atom("OXT")
        .ok_or_else(|| Error::incomplete_for_hydro(&*residue.name, residue.id, "OXT"))?
        .pos;

    let c_oxt_dist = (oxt_pos - c_pos).norm();
    if c_oxt_dist < 1e-6 {
        return Err(Error::incomplete_for_hydro(
            &*residue.name,
            residue.id,
            "OXT",
        ));
    }

    let reference_pos = residue
        .atom("CA")
        .or_else(|| residue.atom("O"))
        .map(|a| a.pos);

    let h_pos = place_hydroxyl_hydrogen(
        oxt_pos,
        c_pos,
        reference_pos,
        COOH_BOND_LENGTH,
        SP3_ANGLE,
        60.0,
    );

    residue.add_atom(Atom::new("HOXT", Element::H, h_pos));
    Ok(())
}

/// Adds the 3'-terminal hydroxyl hydrogen to nucleic acid residues when absent.
///
/// # Arguments
///
/// * `residue` - Residue representing a nucleic acid terminal unit.
///
/// # Returns
///
/// `Ok(())` when the hydrogen is added or already present.
///
/// # Errors
///
/// Returns [`Error::IncompleteResidueForHydro`] if required sugar atoms are missing.
fn construct_3_prime_hydrogen(residue: &mut Residue) -> Result<(), Error> {
    if residue.has_atom("HO3'") {
        return Ok(());
    }

    let o3_pos = residue
        .atom("O3'")
        .ok_or_else(|| Error::incomplete_for_hydro(&*residue.name, residue.id, "O3'"))?
        .pos;
    let c3_pos = residue
        .atom("C3'")
        .ok_or_else(|| Error::incomplete_for_hydro(&*residue.name, residue.id, "C3'"))?
        .pos;

    let reference_pos = residue
        .atom("C4'")
        .or_else(|| residue.atom("C2'"))
        .map(|a| a.pos);

    let h_pos = place_hydroxyl_hydrogen(
        o3_pos,
        c3_pos,
        reference_pos,
        OH_BOND_LENGTH,
        SP3_ANGLE,
        180.0,
    );

    residue.add_atom(Atom::new("HO3'", Element::H, h_pos));
    Ok(())
}

/// Adds the 5'-terminal hydroxyl hydrogen for nucleic acid residues lacking phosphates.
///
/// # Arguments
///
/// * `residue` - Residue whose 5' terminus is being hydrated.
///
/// # Returns
///
/// `Ok(())` when the hydrogen is added or already exists.
///
/// # Errors
///
/// Returns [`Error::IncompleteResidueForHydro`] if sugar anchor atoms are missing.
fn construct_5_prime_hydrogen(residue: &mut Residue) -> Result<(), Error> {
    if residue.has_atom("HO5'") {
        return Ok(());
    }

    let o5_pos = residue
        .atom("O5'")
        .ok_or_else(|| Error::incomplete_for_hydro(&*residue.name, residue.id, "O5'"))?
        .pos;
    let c5_pos = residue
        .atom("C5'")
        .ok_or_else(|| Error::incomplete_for_hydro(&*residue.name, residue.id, "C5'"))?
        .pos;

    let reference_pos = residue.atom("C4'").map(|a| a.pos);

    let h_pos = place_hydroxyl_hydrogen(
        o5_pos,
        c5_pos,
        reference_pos,
        OH_BOND_LENGTH,
        SP3_ANGLE,
        180.0,
    );

    residue.add_atom(Atom::new("HO5'", Element::H, h_pos));
    Ok(())
}

/// Adds hydrogens to 5'-terminal phosphate groups based on pH.
///
/// At physiological pH (≥6.5), the terminal phosphate carries two negative charges and
/// requires no protons. Below this threshold, one proton is added to OP3 with proper
/// sp³ tetrahedral geometry: P-OP3-HOP3 ≈ 109.5°.
///
/// # Arguments
///
/// * `residue` - Nucleic acid residue with a 5'-terminal phosphate.
/// * `config` - Hydrogenation configuration containing pH settings.
///
/// # Returns
///
/// `Ok(())` when hydrogens are appropriately placed or removed.
///
/// # Errors
///
/// Returns [`Error::IncompleteResidueForHydro`] if phosphate atoms are missing.
fn construct_5_prime_phosphate_hydrogens(
    residue: &mut Residue,
    config: &HydroConfig,
) -> Result<(), Error> {
    let ph = config.target_ph.unwrap_or(7.4);

    if ph >= PHOSPHATE_PKA2 {
        residue.remove_atom("HOP3");
        residue.remove_atom("HOP2");
        return Ok(());
    }

    if residue.has_atom("HOP3") {
        return Ok(());
    }

    let op3_pos = residue
        .atom("OP3")
        .ok_or_else(|| Error::incomplete_for_hydro(&*residue.name, residue.id, "OP3"))?
        .pos;
    let p_pos = residue
        .atom("P")
        .ok_or_else(|| Error::incomplete_for_hydro(&*residue.name, residue.id, "P"))?
        .pos;

    let reference_pos = residue
        .atom("OP1")
        .or_else(|| residue.atom("OP2"))
        .map(|a| a.pos);

    let h_pos = place_hydroxyl_hydrogen(
        op3_pos,
        p_pos,
        reference_pos,
        OH_BOND_LENGTH,
        SP3_ANGLE,
        180.0,
    );

    residue.add_atom(Atom::new("HOP3", Element::H, h_pos));
    Ok(())
}

/// Builds an sp³ local coordinate frame centered at `center` with primary axis along
/// `center - attached`.
///
/// # Arguments
///
/// * `center` - Center atom position (e.g., oxygen).
/// * `attached` - Position of the atom bonded to center (e.g., carbon).
/// * `reference` - Optional third atom for defining the xy-plane orientation.
///
/// # Returns
///
/// A tuple of orthonormal vectors `(x, y, z)` where `z` points from attached toward center.
fn build_sp3_frame(
    center: Point,
    attached: Point,
    reference: Option<Point>,
) -> (Vector3<f64>, Vector3<f64>, Vector3<f64>) {
    let z = (center - attached).normalize();

    let ref_vec = reference
        .map(|r| (r - attached).normalize())
        .unwrap_or_else(|| {
            if z.x.abs() < z.y.abs() && z.x.abs() < z.z.abs() {
                Vector3::x()
            } else if z.y.abs() < z.z.abs() {
                Vector3::y()
            } else {
                Vector3::z()
            }
        });

    let x = (ref_vec - z * z.dot(&ref_vec)).normalize();

    let y = z.cross(&x);

    (x, y, z)
}

/// Places a hydroxyl hydrogen using sp³ tetrahedral geometry.
///
/// The hydrogen is placed at `bond_length` from `o_pos`, with the angle
/// `attached-O-H` equal to `bond_angle`, and rotated by `dihedral_offset`
/// around the `attached-O` axis.
///
/// # Arguments
///
/// * `o_pos` - Position of the oxygen atom.
/// * `attached_pos` - Position of the atom bonded to oxygen.
/// * `reference_pos` - Optional reference for determining dihedral orientation.
/// * `bond_length` - O-H bond length (typically 0.96 Å).
/// * `bond_angle` - The angle attached-O-H in degrees (typically 109.5°).
/// * `dihedral_offset` - Rotation around the attached-O axis in degrees.
///
/// # Returns
///
/// The calculated hydrogen position.
fn place_hydroxyl_hydrogen(
    o_pos: Point,
    attached_pos: Point,
    reference_pos: Option<Point>,
    bond_length: f64,
    bond_angle: f64,
    dihedral_offset: f64,
) -> Point {
    let (x, y, z) = build_sp3_frame(o_pos, attached_pos, reference_pos);

    let theta = bond_angle.to_radians();
    let phi = dihedral_offset.to_radians();

    let sin_theta = theta.sin();
    let cos_theta = theta.cos();

    let h_local = Vector3::new(sin_theta * phi.cos(), sin_theta * phi.sin(), -cos_theta);

    let h_global = x * h_local.x + y * h_local.y + z * h_local.z;

    o_pos + h_global * bond_length
}

/// Computes the optimal rigid transform mapping template anchor points to residue atoms.
///
/// Uses Kabsch alignment with safeguards for one- and two-point configurations.
///
/// # Arguments
///
/// * `r_pts` - Coordinates from the residue.
/// * `t_pts` - Corresponding coordinates from the template.
///
/// # Returns
///
/// `Some((rotation, translation))` when alignment succeeds; otherwise `None`.
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

    fn five_prime_residue_with_phosphate(id: i32) -> Residue {
        let template = db::get_template("DA").unwrap();
        let mut residue = Residue::new(
            id,
            None,
            "DA",
            Some(StandardResidue::DA),
            ResidueCategory::Standard,
        );
        residue.position = ResiduePosition::FivePrime;
        for (atom_name, element, pos) in template.heavy_atoms() {
            residue.add_atom(Atom::new(atom_name, element, pos));
        }
        let p_pos = residue.atom("P").unwrap().pos;
        let op1_pos = residue.atom("OP1").unwrap().pos;
        let op2_pos = residue.atom("OP2").unwrap().pos;
        let o5_pos = residue.atom("O5'").unwrap().pos;
        let centroid = (op1_pos.coords + op2_pos.coords + o5_pos.coords) / 3.0;
        let direction = (p_pos.coords - centroid).normalize();
        let op3_pos = p_pos + direction * 1.48;
        residue.add_atom(Atom::new("OP3", Element::O, op3_pos));
        residue
    }

    fn five_prime_residue_without_phosphate(id: i32) -> Residue {
        let template = db::get_template("DA").unwrap();
        let mut residue = Residue::new(
            id,
            None,
            "DA",
            Some(StandardResidue::DA),
            ResidueCategory::Standard,
        );
        residue.position = ResiduePosition::FivePrime;
        for (atom_name, element, pos) in template.heavy_atoms() {
            if !matches!(atom_name, "P" | "OP1" | "OP2") {
                residue.add_atom(Atom::new(atom_name, element, pos));
            }
        }
        residue
    }

    fn three_prime_residue(id: i32) -> Residue {
        let template = db::get_template("DA").unwrap();
        let mut residue = Residue::new(
            id,
            None,
            "DA",
            Some(StandardResidue::DA),
            ResidueCategory::Standard,
        );
        residue.position = ResiduePosition::ThreePrime;
        for (atom_name, element, pos) in template.heavy_atoms() {
            residue.add_atom(Atom::new(atom_name, element, pos));
        }
        residue
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
            determine_protonation_state(
                structure
                    .iter_chains()
                    .next()
                    .unwrap()
                    .iter_residues()
                    .next()
                    .unwrap(),
                &config,
                None,
                None
            ),
            Some("ASH".to_string())
        );

        config.target_ph = Some(5.0);
        assert_eq!(
            determine_protonation_state(
                structure
                    .iter_chains()
                    .next()
                    .unwrap()
                    .iter_residues()
                    .next()
                    .unwrap(),
                &config,
                None,
                None
            ),
            Some("ASP".to_string())
        );

        let structure =
            structure_with_residue(residue_from_template("LYS", StandardResidue::LYS, 2));
        config.target_ph = Some(11.0);
        assert_eq!(
            determine_protonation_state(
                structure
                    .iter_chains()
                    .next()
                    .unwrap()
                    .iter_residues()
                    .next()
                    .unwrap(),
                &config,
                None,
                None
            ),
            Some("LYN".to_string())
        );

        config.target_ph = Some(7.0);
        assert_eq!(
            determine_protonation_state(
                structure
                    .iter_chains()
                    .next()
                    .unwrap()
                    .iter_residues()
                    .next()
                    .unwrap(),
                &config,
                None,
                None
            ),
            Some("LYS".to_string())
        );
    }

    #[test]
    fn determine_protonation_state_respects_his_strategy() {
        let mut residue = residue_from_template("HID", StandardResidue::HIS, 3);
        residue.name = "HIS".into();
        let structure = structure_with_residue(residue);

        let config = HydroConfig {
            target_ph: Some(7.0),
            his_strategy: HisStrategy::DirectHIE,
            ..HydroConfig::default()
        };

        assert_eq!(
            determine_protonation_state(
                structure
                    .iter_chains()
                    .next()
                    .unwrap()
                    .iter_residues()
                    .next()
                    .unwrap(),
                &config,
                None,
                None
            ),
            Some("HIE".to_string())
        );

        let mut acid_config = HydroConfig::default();
        acid_config.target_ph = Some(5.5);
        assert_eq!(
            determine_protonation_state(
                structure
                    .iter_chains()
                    .next()
                    .unwrap()
                    .iter_residues()
                    .next()
                    .unwrap(),
                &acid_config,
                None,
                None
            ),
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
        residue.name = "CYX".into();
        let structure = structure_with_residue(residue);
        let mut config = HydroConfig::default();
        config.target_ph = Some(9.0);

        assert_eq!(
            determine_protonation_state(
                structure
                    .iter_chains()
                    .next()
                    .unwrap()
                    .iter_residues()
                    .next()
                    .unwrap(),
                &config,
                None,
                None
            ),
            None
        );
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

    #[test]
    fn five_prime_phosphate_deprotonated_at_physiological_ph() {
        let residue = five_prime_residue_with_phosphate(60);
        let mut structure = structure_with_residue(residue);

        add_hydrogens(&mut structure, &HydroConfig::default()).expect("hydrogenation succeeds");

        let residue = structure.find_residue("A", 60, None).unwrap();
        assert!(residue.has_atom("OP3"), "OP3 should remain");
        assert!(
            !residue.has_atom("HOP3"),
            "HOP3 should not exist at neutral pH"
        );
        assert!(
            !residue.has_atom("HOP2"),
            "HOP2 should not exist at neutral pH"
        );
    }

    #[test]
    fn five_prime_phosphate_protonated_below_pka() {
        let residue = five_prime_residue_with_phosphate(61);
        let mut structure = structure_with_residue(residue);
        let mut config = HydroConfig::default();
        config.target_ph = Some(5.5);

        add_hydrogens(&mut structure, &config).expect("hydrogenation succeeds");

        let residue = structure.find_residue("A", 61, None).unwrap();
        assert!(residue.has_atom("OP3"), "OP3 should remain");
        assert!(residue.has_atom("HOP3"), "HOP3 should be added below pKa");
    }

    #[test]
    fn five_prime_without_phosphate_gets_ho5() {
        let residue = five_prime_residue_without_phosphate(62);
        let mut structure = structure_with_residue(residue);

        add_hydrogens(&mut structure, &HydroConfig::default()).expect("hydrogenation succeeds");

        let residue = structure.find_residue("A", 62, None).unwrap();
        assert!(
            residue.has_atom("HO5'"),
            "HO5' should be added for 5'-OH terminus"
        );
        assert!(!residue.has_atom("P"), "phosphorus should not exist");
    }

    #[test]
    fn three_prime_nucleic_gets_ho3() {
        let residue = three_prime_residue(70);
        let mut structure = structure_with_residue(residue);

        add_hydrogens(&mut structure, &HydroConfig::default()).expect("hydrogenation succeeds");

        let residue = structure.find_residue("A", 70, None).unwrap();
        assert!(
            residue.has_atom("HO3'"),
            "HO3' should be added for 3' terminal"
        );
    }

    #[test]
    fn n_terminal_h_has_tetrahedral_geometry() {
        let residue = n_terminal_residue(98);
        let mut structure = structure_with_residue(residue);

        add_hydrogens(&mut structure, &HydroConfig::default()).expect("hydrogenation succeeds");

        let residue = structure.find_residue("A", 98, None).unwrap();
        let n = residue.atom("N").expect("N").pos;
        let ca = residue.atom("CA").expect("CA").pos;
        let h1 = residue.atom("H1").expect("H1").pos;
        let h2 = residue.atom("H2").expect("H2").pos;
        let h3 = residue.atom("H3").expect("H3").pos;

        let n_h1_dist = distance(n, h1);
        let n_h2_dist = distance(n, h2);
        let n_h3_dist = distance(n, h3);
        assert!(
            (n_h1_dist - NH_BOND_LENGTH).abs() < 0.1,
            "N-H1 distance {n_h1_dist:.3} should be ~{NH_BOND_LENGTH} Å"
        );
        assert!(
            (n_h2_dist - NH_BOND_LENGTH).abs() < 0.1,
            "N-H2 distance {n_h2_dist:.3} should be ~{NH_BOND_LENGTH} Å"
        );
        assert!(
            (n_h3_dist - NH_BOND_LENGTH).abs() < 0.1,
            "N-H3 distance {n_h3_dist:.3} should be ~{NH_BOND_LENGTH} Å"
        );

        let ca_n_h1_angle = angle_deg(ca, n, h1);
        let ca_n_h2_angle = angle_deg(ca, n, h2);
        let ca_n_h3_angle = angle_deg(ca, n, h3);
        assert!(
            (ca_n_h1_angle - SP3_ANGLE).abs() < 5.0,
            "CA-N-H1 angle {ca_n_h1_angle:.1}° should be ~{SP3_ANGLE}°"
        );
        assert!(
            (ca_n_h2_angle - SP3_ANGLE).abs() < 5.0,
            "CA-N-H2 angle {ca_n_h2_angle:.1}° should be ~{SP3_ANGLE}°"
        );
        assert!(
            (ca_n_h3_angle - SP3_ANGLE).abs() < 5.0,
            "CA-N-H3 angle {ca_n_h3_angle:.1}° should be ~{SP3_ANGLE}°"
        );

        let h1_n_h2_angle = angle_deg(h1, n, h2);
        let h2_n_h3_angle = angle_deg(h2, n, h3);
        let h1_n_h3_angle = angle_deg(h1, n, h3);
        assert!(
            (h1_n_h2_angle - SP3_ANGLE).abs() < 5.0,
            "H1-N-H2 angle {h1_n_h2_angle:.1}° should be ~{SP3_ANGLE}°"
        );
        assert!(
            (h2_n_h3_angle - SP3_ANGLE).abs() < 5.0,
            "H2-N-H3 angle {h2_n_h3_angle:.1}° should be ~{SP3_ANGLE}°"
        );
        assert!(
            (h1_n_h3_angle - SP3_ANGLE).abs() < 5.0,
            "H1-N-H3 angle {h1_n_h3_angle:.1}° should be ~{SP3_ANGLE}°"
        );
    }

    #[test]
    fn c_terminal_hoxt_has_tetrahedral_geometry() {
        let residue = c_terminal_residue(99);
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: Some(2.0),
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).expect("hydrogenation succeeds");

        let residue = structure.find_residue("A", 99, None).unwrap();
        let c = residue.atom("C").expect("C").pos;
        let oxt = residue.atom("OXT").expect("OXT").pos;
        let hoxt = residue.atom("HOXT").expect("HOXT").pos;

        let oxt_hoxt_dist = distance(oxt, hoxt);
        assert!(
            (oxt_hoxt_dist - COOH_BOND_LENGTH).abs() < 0.1,
            "OXT-HOXT distance {oxt_hoxt_dist:.3} should be ~{COOH_BOND_LENGTH} Å"
        );

        let c_oxt_hoxt_angle = angle_deg(c, oxt, hoxt);
        assert!(
            (c_oxt_hoxt_angle - SP3_ANGLE).abs() < 5.0,
            "C-OXT-HOXT angle {c_oxt_hoxt_angle:.1}° should be ~{SP3_ANGLE}°"
        );
    }

    #[test]
    fn five_prime_ho5_has_tetrahedral_geometry() {
        let residue = five_prime_residue_without_phosphate(80);
        let mut structure = structure_with_residue(residue);

        add_hydrogens(&mut structure, &HydroConfig::default()).expect("hydrogenation succeeds");

        let residue = structure.find_residue("A", 80, None).unwrap();
        let c5 = residue.atom("C5'").expect("C5'").pos;
        let o5 = residue.atom("O5'").expect("O5'").pos;
        let ho5 = residue.atom("HO5'").expect("HO5'").pos;

        let o5_ho5_dist = distance(o5, ho5);
        assert!(
            (o5_ho5_dist - OH_BOND_LENGTH).abs() < 0.1,
            "O5'-HO5' distance {o5_ho5_dist:.3} should be ~{OH_BOND_LENGTH} Å"
        );

        let c5_o5_ho5_angle = angle_deg(c5, o5, ho5);
        assert!(
            (c5_o5_ho5_angle - SP3_ANGLE).abs() < 5.0,
            "C5'-O5'-HO5' angle {c5_o5_ho5_angle:.1}° should be ~{SP3_ANGLE}°"
        );
    }

    #[test]
    fn three_prime_ho3_has_tetrahedral_geometry() {
        let residue = three_prime_residue(81);
        let mut structure = structure_with_residue(residue);

        add_hydrogens(&mut structure, &HydroConfig::default()).expect("hydrogenation succeeds");

        let residue = structure.find_residue("A", 81, None).unwrap();
        let c3 = residue.atom("C3'").expect("C3'").pos;
        let o3 = residue.atom("O3'").expect("O3'").pos;
        let ho3 = residue.atom("HO3'").expect("HO3'").pos;

        let o3_ho3_dist = distance(o3, ho3);
        assert!(
            (o3_ho3_dist - OH_BOND_LENGTH).abs() < 0.1,
            "O3'-HO3' distance {o3_ho3_dist:.3} should be ~{OH_BOND_LENGTH} Å"
        );

        let c3_o3_ho3_angle = angle_deg(c3, o3, ho3);
        assert!(
            (c3_o3_ho3_angle - SP3_ANGLE).abs() < 5.0,
            "C3'-O3'-HO3' angle {c3_o3_ho3_angle:.1}° should be ~{SP3_ANGLE}°"
        );
    }

    #[test]
    fn five_prime_phosphate_hop3_has_tetrahedral_geometry() {
        let residue = five_prime_residue_with_phosphate(82);
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: Some(5.5),
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).expect("hydrogenation succeeds");

        let residue = structure.find_residue("A", 82, None).unwrap();
        let p = residue.atom("P").expect("P").pos;
        let op3 = residue.atom("OP3").expect("OP3").pos;
        let hop3 = residue.atom("HOP3").expect("HOP3").pos;

        let op3_hop3_dist = distance(op3, hop3);
        assert!(
            (op3_hop3_dist - OH_BOND_LENGTH).abs() < 0.1,
            "OP3-HOP3 distance {op3_hop3_dist:.3} should be ~{OH_BOND_LENGTH} Å"
        );

        let p_op3_hop3_angle = angle_deg(p, op3, hop3);
        assert!(
            (p_op3_hop3_angle - SP3_ANGLE).abs() < 5.0,
            "P-OP3-HOP3 angle {p_op3_hop3_angle:.1}° should be ~{SP3_ANGLE}°"
        );
    }
}
