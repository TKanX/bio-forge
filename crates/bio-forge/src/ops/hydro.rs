//! Protonation-aware hydrogen placement for protein and nucleic acid structures.
//!
//! This module inspects residue templates, predicts protonation states from pH or hydrogen
//! bonding networks, relabels residues accordingly, and rebuilds missing hydrogens while
//! respecting polymer termini, nucleic acid priming, disulfide bridges, and salt bridges.

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

/// Henderson–Hasselbalch breakpoint for histidine double-protonation (HIP).
const HIS_HIP_PKA: f64 = 6.0;
/// Henderson–Hasselbalch breakpoint for aspartate protonation (ASH).
const ASP_PKA: f64 = 3.9;
/// Henderson–Hasselbalch breakpoint for glutamate protonation (GLH).
const GLU_PKA: f64 = 4.2;
/// Henderson–Hasselbalch breakpoint for lysine deprotonation (LYN).
const LYS_PKA: f64 = 10.5;
/// Henderson–Hasselbalch breakpoint for arginine deprotonation (ARN).
const ARG_PKA: f64 = 12.5;
/// Henderson–Hasselbalch breakpoint for cysteine deprotonation (CYM).
const CYS_PKA: f64 = 8.3;
/// Henderson–Hasselbalch breakpoint for tyrosine deprotonation (TYM).
const TYR_PKA: f64 = 10.0;
/// Henderson–Hasselbalch breakpoint for protonated N-termini.
const N_TERM_PKA: f64 = 8.0;
/// Henderson–Hasselbalch breakpoint for protonated C-termini.
const C_TERM_PKA: f64 = 3.1;
/// Henderson–Hasselbalch breakpoint for the second dissociation of terminal phosphate.
const PHOSPHATE_PKA2: f64 = 6.5;
/// Default assumed pH for terminal state decisions when no explicit pH is specified.
const DEFAULT_TERMINAL_PH: f64 = 7.0;
/// Maximum sulfur–sulfur distance (Å) for disulfide bridge detection.
const DISULFIDE_SG_THRESHOLD: f64 = 2.2;
/// Maximum nitrogen–oxygen distance (Å) for HIS-carboxylate salt bridge detection.
const SALT_BRIDGE_DISTANCE: f64 = 4.0;
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
/// `HydroConfig` can target a specific solution pH, remove pre-existing hydrogens,
/// choose how neutral histidine tautomers are assigned, and enable/disable salt bridge detection.
#[derive(Debug, Clone)]
pub struct HydroConfig {
    /// Optional solvent pH value used for titration decisions.
    pub target_ph: Option<f64>,
    /// Whether to strip all existing hydrogens before reconstruction.
    pub remove_existing_h: bool,
    /// Strategy for selecting neutral histidine tautomers (HID/HIE).
    pub his_strategy: HisStrategy,
    /// Whether to protonate histidine to HIP when forming salt bridges with
    /// nearby carboxylate groups (ASP⁻/GLU⁻/C-terminal COO⁻).
    pub his_salt_bridge_protonation: bool,
}

impl Default for HydroConfig {
    /// Provides biologically reasonable defaults (physiological pH, removal of old hydrogens,
    /// hydrogen-bond-aware histidine selection, and HIS salt bridge detection).
    fn default() -> Self {
        Self {
            target_ph: None,
            remove_existing_h: true,
            his_strategy: HisStrategy::HbNetwork,
            his_salt_bridge_protonation: true,
        }
    }
}

/// Strategies for selecting neutral histidine tautomers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HisStrategy {
    /// Force all neutral histidines to HID (δ-protonated).
    DirectHID,
    /// Force all neutral histidines to HIE (ε-protonated).
    DirectHIE,
    /// Randomly choose between HID and HIE with equal probability.
    Random,
    /// Analyze hydrogen-bond networks to select the tautomer most likely to form
    /// favorable interactions with nearby acceptors.
    HbNetwork,
}

/// Adds hydrogens to all standard residues in-place, updating protonation states when needed.
///
/// This function implements a multi-phase pipeline:
///
/// 1. **Disulfide detection** — Identifies CYS pairs forming S-S bonds and relabels to CYX.
/// 2. **Non-HIS protonation** — Applies pKa-based titration to ASP, GLU, LYS, ARG, CYS, TYR
///    (only when `target_ph` is specified).
/// 3. **HIS protonation** — Determines HIS state via pH thresholds, salt bridge detection,
///    and tautomer strategy.
/// 4. **Hydrogen construction** — Builds hydrogens according to template geometry and
///    terminal-specific rules.
///
/// # Arguments
///
/// * `structure` - Mutable structure whose residues will be protonated and hydrated.
/// * `config` - Hydrogenation configuration controlling pH, strategy, and options.
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

    let acceptor_grid =
        if config.his_strategy == HisStrategy::HbNetwork && config.target_ph.is_some() {
            Some(build_acceptor_grid(structure))
        } else {
            None
        };

    if config.target_ph.is_some() {
        apply_non_his_protonation(structure, config.target_ph.unwrap());
    }

    let carboxylate_grid = if config.his_salt_bridge_protonation {
        Some(build_carboxylate_grid(structure, config.target_ph))
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

                    if residue.standard_name == Some(StandardResidue::HIS) {
                        if let Some(new_name) = determine_his_protonation(
                            residue,
                            config,
                            acceptor_grid.as_ref(),
                            carboxylate_grid.as_ref(),
                            (c_idx, r_idx),
                        ) {
                            residue.name = new_name.into();
                        }
                    }

                    if config.remove_existing_h {
                        residue.strip_hydrogens();
                    }

                    construct_hydrogens_for_residue(residue, config)
                })
        })
}

/// Applies pH-based protonation to all non-HIS titratable residues.
///
/// CYX (disulfide-bonded cysteine) is never modified.
///
/// # Arguments
///
/// * `structure` - Mutable structure whose residues will be protonated.
/// * `ph` - Target pH for protonation decisions.
fn apply_non_his_protonation(structure: &mut Structure, ph: f64) {
    structure.par_residues_mut().for_each(|residue| {
        if residue.category != ResidueCategory::Standard {
            return;
        }

        if residue.name == "CYX" {
            return;
        }

        let new_name = match residue.standard_name {
            Some(StandardResidue::ASP) => Some(if ph < ASP_PKA { "ASH" } else { "ASP" }),
            Some(StandardResidue::GLU) => Some(if ph < GLU_PKA { "GLH" } else { "GLU" }),
            Some(StandardResidue::LYS) => Some(if ph > LYS_PKA { "LYN" } else { "LYS" }),
            Some(StandardResidue::ARG) => Some(if ph > ARG_PKA { "ARN" } else { "ARG" }),
            Some(StandardResidue::CYS) => Some(if ph > CYS_PKA { "CYM" } else { "CYS" }),
            Some(StandardResidue::TYR) => Some(if ph > TYR_PKA { "TYM" } else { "TYR" }),
            _ => None,
        };

        if let Some(name) = new_name {
            residue.name = name.into();
        }
    });
}

/// Builds a spatial grid of carboxylate oxygen atoms for salt bridge detection.
///
/// ASH/GLH (protonated carboxyls) are excluded because neutral COOH groups
/// cannot form ionic salt bridges.
///
/// # Arguments
///
/// * `structure` - Structure from which to extract carboxylate oxygens.
/// * `target_ph` - Optional pH used to determine C-terminal protonation.
///
/// # Returns
///
/// Spatial grid of carboxylate oxygen positions mapped to (chain_idx, residue_idx).
fn build_carboxylate_grid(structure: &Structure, target_ph: Option<f64>) -> Grid<(usize, usize)> {
    let c_term_deprotonated = c_terminus_is_deprotonated(target_ph);

    let atoms: Vec<(Point, (usize, usize))> = structure
        .iter_chains()
        .enumerate()
        .flat_map(|(c_idx, chain)| {
            chain
                .iter_residues()
                .enumerate()
                .flat_map(move |(r_idx, residue)| {
                    let mut positions = Vec::new();

                    // ASP⁻ carboxylate oxygens
                    if residue.name == "ASP" {
                        if let Some(od1) = residue.atom("OD1") {
                            positions.push((od1.pos, (c_idx, r_idx)));
                        }
                        if let Some(od2) = residue.atom("OD2") {
                            positions.push((od2.pos, (c_idx, r_idx)));
                        }
                    }

                    // GLU⁻ carboxylate oxygens
                    if residue.name == "GLU" {
                        if let Some(oe1) = residue.atom("OE1") {
                            positions.push((oe1.pos, (c_idx, r_idx)));
                        }
                        if let Some(oe2) = residue.atom("OE2") {
                            positions.push((oe2.pos, (c_idx, r_idx)));
                        }
                    }

                    // C-terminal COO⁻ (only if deprotonated)
                    if residue.position == ResiduePosition::CTerminal
                        && residue.standard_name.is_some_and(|s| s.is_protein())
                        && c_term_deprotonated
                    {
                        if let Some(o) = residue.atom("O") {
                            positions.push((o.pos, (c_idx, r_idx)));
                        }
                        if let Some(oxt) = residue.atom("OXT") {
                            positions.push((oxt.pos, (c_idx, r_idx)));
                        }
                    }

                    positions
                })
        })
        .collect();

    Grid::new(atoms, SALT_BRIDGE_DISTANCE + 0.5)
}

/// Determines the protonation state for a histidine residue.
///
/// # Decision Tree
///
/// 1. **pH < 6.0** → HIP (doubly protonated, +1 charge)
/// 2. **No pH AND no salt bridge detection** → `None` (preserve user-defined name)
/// 3. **Salt bridge detected** → HIP
/// 4. **No pH** → `None` (salt bridge didn't trigger, preserve name)
/// 5. **pH ≥ 6.0, no salt bridge** → Apply HisStrategy (HID/HIE)
///
/// # Arguments
///
/// * `residue` - Histidine residue to evaluate.
/// * `config` - Hydrogenation configuration.
/// * `acceptor_grid` - Optional spatial grid of hydrogen bond acceptors.
/// * `carboxylate_grid` - Optional spatial grid of carboxylate oxygens.
/// * `self_indices` - Tuple of (chain_idx, residue_idx) for the current residue.
///
/// # Returns
///
/// `Some(new_name)` when the residue should be renamed, `None` to preserve current name.
fn determine_his_protonation(
    residue: &Residue,
    config: &HydroConfig,
    acceptor_grid: Option<&Grid<(usize, usize)>>,
    carboxylate_grid: Option<&Grid<(usize, usize)>>,
    self_indices: (usize, usize),
) -> Option<String> {
    if let Some(ph) = config.target_ph {
        if ph < HIS_HIP_PKA {
            return Some("HIP".to_string());
        }
    }

    if config.target_ph.is_none() && !config.his_salt_bridge_protonation {
        return None;
    }

    if config.his_salt_bridge_protonation {
        if let Some(grid) = carboxylate_grid {
            if his_forms_salt_bridge(residue, grid, self_indices) {
                return Some("HIP".to_string());
            }
        }
    }

    if config.target_ph.is_none() {
        return None;
    }

    Some(select_neutral_his(
        residue,
        config.his_strategy,
        acceptor_grid,
        self_indices,
    ))
}

/// Detects if a histidine forms a salt bridge with nearby carboxylate groups.
///
/// Checks if either ND1 or NE2 nitrogen is within [`SALT_BRIDGE_DISTANCE`] of
/// any carboxylate oxygen in the grid.
///
/// # Arguments
///
/// * `residue` - Histidine residue to evaluate.
/// * `carboxylate_grid` - Spatial grid of carboxylate oxygen positions.
/// * `self_indices` - Tuple of (chain_idx, residue_idx) for the current residue.
///
/// # Returns
///
/// `true` if a salt bridge is detected, `false` otherwise.
fn his_forms_salt_bridge(
    residue: &Residue,
    carboxylate_grid: &Grid<(usize, usize)>,
    self_indices: (usize, usize),
) -> bool {
    // Check ND1
    if let Some(nd1) = residue.atom("ND1") {
        for (_, &idx) in carboxylate_grid
            .neighbors(&nd1.pos, SALT_BRIDGE_DISTANCE)
            .exact()
        {
            if idx != self_indices {
                return true;
            }
        }
    }

    // Check NE2
    if let Some(ne2) = residue.atom("NE2") {
        for (_, &idx) in carboxylate_grid
            .neighbors(&ne2.pos, SALT_BRIDGE_DISTANCE)
            .exact()
        {
            if idx != self_indices {
                return true;
            }
        }
    }

    false
}

/// Selects a neutral histidine tautomer using the configured strategy.
///
/// # Arguments
///
/// * `residue` - Histidine residue to evaluate.
/// * `strategy` - Strategy for selecting the tautomer.
/// * `acceptor_grid` - Optional spatial grid of hydrogen bond acceptors.
/// * `self_indices` - Tuple of (chain_idx, residue_idx) for the current residue.
///
/// # Returns
///
/// `"HID"` or `"HIE"` depending on the selected tautomer.
fn select_neutral_his(
    residue: &Residue,
    strategy: HisStrategy,
    acceptor_grid: Option<&Grid<(usize, usize)>>,
    self_indices: (usize, usize),
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
        HisStrategy::HbNetwork => optimize_his_network(residue, acceptor_grid, self_indices),
    }
}

/// Determines the best histidine tautomer by inspecting nearby hydrogen-bond acceptors.
///
/// Computes hypothetical H-bond scores for both ND1 (HID) and NE2 (HIE) protonation
/// and returns the tautomer with the higher score.
///
/// # Arguments
///
/// * `residue` - Histidine residue to evaluate.
/// * `acceptor_grid` - Optional spatial grid of hydrogen bond acceptors.
/// * `self_indices` - Tuple of (chain_idx, residue_idx) for the current residue.
///
/// # Returns
///
/// `"HID"` or `"HIE"` depending on which tautomer has a better H-bonding score.
fn optimize_his_network(
    residue: &Residue,
    acceptor_grid: Option<&Grid<(usize, usize)>>,
    self_indices: (usize, usize),
) -> String {
    let grid = match acceptor_grid {
        Some(g) => g,
        None => return "HIE".to_string(),
    };

    let score_hid = calculate_h_bond_score(residue, "ND1", "CG", "CE1", grid, self_indices);
    let score_hie = calculate_h_bond_score(residue, "NE2", "CD2", "CE1", grid, self_indices);

    if score_hid > score_hie {
        "HID".to_string()
    } else {
        "HIE".to_string()
    }
}

/// Calculates a geometric score for a potential hydrogen bond from an sp² nitrogen.
///
/// The score considers both distance (< 2.7 Å) and angle (> 90°) criteria for
/// valid hydrogen bonds.
///
/// # Arguments
///
/// * `residue` - Histidine residue containing the nitrogen.
/// * `n_name` - Name of the nitrogen atom (e.g., "ND1" or "NE2").
/// * `c1_name` - Name of the first carbon anchor atom.
/// * `c2_name` - Name of the second carbon anchor atom.
/// * `grid` - Spatial grid of potential acceptor atoms.
/// * `self_idx` - Tuple of (chain_idx, residue_idx) for the current residue.
///
/// # Returns
///
/// Hydrogen bond score as a floating-point value.
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

/// Identifies cysteine pairs forming disulfide bonds and renames them to CYX.
///
/// # Arguments
///
/// * `structure` - Mutable structure to analyze and modify.
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

/// Builds a spatial grid of all nitrogen and oxygen atoms for H-bond network analysis.
///
/// # Arguments
///
/// * `structure` - Structure from which to extract acceptor atoms.
///
/// # Returns
///
/// Spatial grid of acceptor atom (N, O, F) positions mapped to (chain_idx, residue_idx).
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

/// Rebuilds hydrogens for a single residue using template geometry and terminal rules.
///
/// # Arguments
///
/// * `residue` - Mutable residue to which hydrogens will be added.
/// * `config` - Hydrogenation configuration controlling pH and options.
///
/// # Returns
///
/// `Ok(())` when hydrogen construction succeeds.
///
/// # Errors
///
/// Returns [`Error::MissingInternalTemplate`] when no template is found or
/// [`Error::IncompleteResidueForHydro`] when required anchor atoms are missing.
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
            construct_n_term_hydrogens(residue, n_term_is_protonated(config.target_ph))?;
        }
        ResiduePosition::CTerminal if residue.standard_name.is_some_and(|s| s.is_protein()) => {
            construct_c_term_hydrogen(residue, c_term_is_protonated(config.target_ph))?;
        }
        ResiduePosition::ThreePrime if residue.standard_name.is_some_and(|s| s.is_nucleic()) => {
            construct_3_prime_hydrogen(residue)?;
        }
        ResiduePosition::FivePrime if residue.standard_name.is_some_and(|s| s.is_nucleic()) => {
            if residue.has_atom("P") {
                construct_5_prime_phosphate_hydrogens(residue, config.target_ph)?;
            } else if residue.has_atom("O5'") {
                construct_5_prime_hydrogen(residue)?;
            }
        }
        _ => {}
    }

    Ok(())
}
