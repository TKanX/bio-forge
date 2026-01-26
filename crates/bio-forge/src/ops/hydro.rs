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
    /// Provides biologically reasonable defaults (no protonation state changes,
    /// removal of old hydrogens, hydrogen-bond-aware histidine selection,
    /// and HIS salt bridge detection).
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

    let acceptor_grid = if config.his_strategy == HisStrategy::HbNetwork
        && config.target_ph.is_some_and(|ph| ph >= HIS_HIP_PKA)
    {
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

                    if let Some(StandardResidue::HIS) = residue.standard_name
                        && let Some(new_name) = determine_his_protonation(
                            residue,
                            config,
                            acceptor_grid.as_ref(),
                            carboxylate_grid.as_ref(),
                            (c_idx, r_idx),
                        )
                    {
                        residue.name = new_name.into();
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
        .par_chains()
        .enumerate()
        .flat_map(|(c_idx, chain)| {
            chain
                .par_residues()
                .enumerate()
                .flat_map_iter(move |(r_idx, residue)| {
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
    if let Some(ph) = config.target_ph
        && ph < HIS_HIP_PKA
    {
        return Some("HIP".to_string());
    }

    if config.target_ph.is_none() && !config.his_salt_bridge_protonation {
        return None;
    }

    if config.his_salt_bridge_protonation
        && let Some(grid) = carboxylate_grid
        && his_forms_salt_bridge(residue, grid, self_indices)
    {
        return Some("HIP".to_string());
    }

    config.target_ph?;

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
        .par_chains()
        .enumerate()
        .flat_map(|(c_idx, chain)| {
            chain
                .par_residues()
                .enumerate()
                .flat_map_iter(move |(r_idx, residue)| {
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

/// Returns the effective pH used for terminal protonation state decisions.
#[inline]
fn effective_terminal_ph(target_ph: Option<f64>) -> f64 {
    target_ph.unwrap_or(DEFAULT_TERMINAL_PH)
}

/// Determines if a C-terminus should be considered deprotonated (COO⁻).
#[inline]
fn c_terminus_is_deprotonated(target_ph: Option<f64>) -> bool {
    effective_terminal_ph(target_ph) >= C_TERM_PKA
}

/// Evaluates whether an N-terminus should be protonated (NH₃⁺).
#[inline]
fn n_term_is_protonated(target_ph: Option<f64>) -> bool {
    effective_terminal_ph(target_ph) <= N_TERM_PKA
}

/// Evaluates whether a C-terminus should be protonated (COOH).
#[inline]
fn c_term_is_protonated(target_ph: Option<f64>) -> bool {
    effective_terminal_ph(target_ph) < C_TERM_PKA
}

/// Rebuilds the N-terminal amine hydrogens using tetrahedral geometry.
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

/// Adds the 3'-terminal hydroxyl hydrogen to nucleic acid residues.
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
/// At physiological pH (≥6.5), the terminal phosphate carries two negative charges
/// and requires no protons. Below this threshold, one proton is added to OP3.
fn construct_5_prime_phosphate_hydrogens(
    residue: &mut Residue,
    target_ph: Option<f64>,
) -> Result<(), Error> {
    let ph = effective_terminal_ph(target_ph);

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

/// Aligns template coordinates to residue anchors to predict a hydrogen position.
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

/// Builds an sp³ local coordinate frame centered at `center` with primary axis along
/// `center - attached`.
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

    fn his_near_asp(his_id: i32, asp_id: i32, distance: f64) -> Structure {
        let his = residue_from_template("HID", StandardResidue::HIS, his_id);
        let mut asp = residue_from_template("ASP", StandardResidue::ASP, asp_id);

        let his_nd1 = his.atom("ND1").expect("ND1").pos;
        let asp_od1 = asp.atom("OD1").expect("OD1").pos;
        let offset = his_nd1 + Vector3::new(distance, 0.0, 0.0) - asp_od1;
        for atom in asp.iter_atoms_mut() {
            atom.translate_by(&offset);
        }

        structure_with_residues(vec![his, asp])
    }

    fn his_near_glu(his_id: i32, glu_id: i32, distance: f64) -> Structure {
        let his = residue_from_template("HID", StandardResidue::HIS, his_id);
        let mut glu = residue_from_template("GLU", StandardResidue::GLU, glu_id);

        let his_ne2 = his.atom("NE2").expect("NE2").pos;
        let glu_oe1 = glu.atom("OE1").expect("OE1").pos;
        let offset = his_ne2 + Vector3::new(distance, 0.0, 0.0) - glu_oe1;
        for atom in glu.iter_atoms_mut() {
            atom.translate_by(&offset);
        }

        structure_with_residues(vec![his, glu])
    }

    fn his_near_c_term(his_id: i32, c_term_id: i32, distance: f64) -> Structure {
        let his = residue_from_template("HID", StandardResidue::HIS, his_id);
        let mut c_term = c_terminal_residue(c_term_id);

        let his_nd1 = his.atom("ND1").expect("ND1").pos;
        let c_term_oxt = c_term.atom("OXT").expect("OXT").pos;
        let offset = his_nd1 + Vector3::new(distance, 0.0, 0.0) - c_term_oxt;
        for atom in c_term.iter_atoms_mut() {
            atom.translate_by(&offset);
        }

        structure_with_residues(vec![his, c_term])
    }

    fn his_isolated(id: i32) -> Structure {
        let his = residue_from_template("HID", StandardResidue::HIS, id);
        structure_with_residue(his)
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
    fn asp_protonates_to_ash_below_pka() {
        let residue = residue_from_template("ASP", StandardResidue::ASP, 1);
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: Some(2.5),
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(res.name, "ASH", "ASP should become ASH below pKa 3.9");
    }

    #[test]
    fn asp_remains_deprotonated_above_pka() {
        let residue = residue_from_template("ASP", StandardResidue::ASP, 1);
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: Some(7.4),
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(res.name, "ASP", "ASP should remain ASP above pKa 3.9");
    }

    #[test]
    fn asp_preserves_original_name_without_ph() {
        let residue = residue_from_template("ASP", StandardResidue::ASP, 1);
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: None,
            his_salt_bridge_protonation: false,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(res.name, "ASP", "ASP should be preserved without pH");
    }

    #[test]
    fn glu_protonates_to_glh_below_pka() {
        let residue = residue_from_template("GLU", StandardResidue::GLU, 1);
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: Some(3.0),
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(res.name, "GLH", "GLU should become GLH below pKa 4.2");
    }

    #[test]
    fn glu_remains_deprotonated_above_pka() {
        let residue = residue_from_template("GLU", StandardResidue::GLU, 1);
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: Some(7.4),
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(res.name, "GLU", "GLU should remain GLU above pKa 4.2");
    }

    #[test]
    fn glu_preserves_original_name_without_ph() {
        let residue = residue_from_template("GLU", StandardResidue::GLU, 1);
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: None,
            his_salt_bridge_protonation: false,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(res.name, "GLU", "GLU should be preserved without pH");
    }

    #[test]
    fn lys_deprotonates_to_lyn_above_pka() {
        let residue = residue_from_template("LYS", StandardResidue::LYS, 1);
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: Some(11.0),
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(res.name, "LYN", "LYS should become LYN above pKa 10.5");
    }

    #[test]
    fn lys_remains_protonated_below_pka() {
        let residue = residue_from_template("LYS", StandardResidue::LYS, 1);
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: Some(7.4),
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(res.name, "LYS", "LYS should remain LYS below pKa 10.5");
    }

    #[test]
    fn lys_preserves_original_name_without_ph() {
        let residue = residue_from_template("LYS", StandardResidue::LYS, 1);
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: None,
            his_salt_bridge_protonation: false,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(res.name, "LYS", "LYS should be preserved without pH");
    }

    #[test]
    fn cys_deprotonates_to_cym_above_pka() {
        let residue = residue_from_template("CYS", StandardResidue::CYS, 1);
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: Some(9.0),
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(res.name, "CYM", "CYS should become CYM above pKa 8.3");
    }

    #[test]
    fn cys_remains_protonated_below_pka() {
        let residue = residue_from_template("CYS", StandardResidue::CYS, 1);
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: Some(7.4),
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(res.name, "CYS", "CYS should remain CYS below pKa 8.3");
    }

    #[test]
    fn cys_preserves_original_name_without_ph() {
        let residue = residue_from_template("CYS", StandardResidue::CYS, 1);
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: None,
            his_salt_bridge_protonation: false,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(res.name, "CYS", "CYS should be preserved without pH");
    }

    #[test]
    fn cyx_is_preserved_regardless_of_ph() {
        let mut residue = residue_from_template("CYS", StandardResidue::CYS, 1);
        residue.name = "CYX".into();
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: Some(9.0),
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(res.name, "CYX", "CYX should be preserved regardless of pH");
    }

    #[test]
    fn tyr_deprotonates_to_tym_above_pka() {
        let residue = residue_from_template("TYR", StandardResidue::TYR, 1);
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: Some(11.0),
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(res.name, "TYM", "TYR should become TYM above pKa 10.0");
    }

    #[test]
    fn tyr_remains_protonated_below_pka() {
        let residue = residue_from_template("TYR", StandardResidue::TYR, 1);
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: Some(7.4),
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(res.name, "TYR", "TYR should remain TYR below pKa 10.0");
    }

    #[test]
    fn tyr_preserves_original_name_without_ph() {
        let residue = residue_from_template("TYR", StandardResidue::TYR, 1);
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: None,
            his_salt_bridge_protonation: false,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(res.name, "TYR", "TYR should be preserved without pH");
    }

    #[test]
    fn arg_deprotonates_to_arn_above_pka() {
        let residue = residue_from_template("ARG", StandardResidue::ARG, 1);
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: Some(13.0),
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(res.name, "ARN", "ARG should become ARN above pKa 12.5");
    }

    #[test]
    fn arg_remains_protonated_below_pka() {
        let residue = residue_from_template("ARG", StandardResidue::ARG, 1);
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: Some(7.4),
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(res.name, "ARG", "ARG should remain ARG below pKa 12.5");
    }

    #[test]
    fn arg_preserves_original_name_without_ph() {
        let residue = residue_from_template("ARG", StandardResidue::ARG, 1);
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: None,
            his_salt_bridge_protonation: false,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(res.name, "ARG", "ARG should be preserved without pH");
    }

    #[test]
    fn coo_grid_includes_asp_oxygens_when_deprotonated() {
        let residue = residue_from_template("ASP", StandardResidue::ASP, 1);
        let structure = structure_with_residue(residue);

        let grid = build_carboxylate_grid(&structure, Some(7.4));
        let asp = structure.find_residue("A", 1, None).unwrap();
        let od1 = asp.atom("OD1").unwrap().pos;

        let neighbors: Vec<_> = grid.neighbors(&od1, 0.1).exact().collect();
        assert!(!neighbors.is_empty(), "ASP OD1 should be in COO⁻ grid");
    }

    #[test]
    fn coo_grid_excludes_ash_oxygens_when_protonated() {
        let mut residue = residue_from_template("ASP", StandardResidue::ASP, 1);
        residue.name = "ASH".into();
        let structure = structure_with_residue(residue);

        let grid = build_carboxylate_grid(&structure, Some(7.4));
        let ash = structure.find_residue("A", 1, None).unwrap();
        let od1 = ash.atom("OD1").unwrap().pos;

        let neighbors: Vec<_> = grid.neighbors(&od1, 0.1).exact().collect();
        assert!(
            neighbors.is_empty(),
            "ASH oxygens should NOT be in COO⁻ grid"
        );
    }

    #[test]
    fn coo_grid_includes_glu_oxygens_when_deprotonated() {
        let residue = residue_from_template("GLU", StandardResidue::GLU, 1);
        let structure = structure_with_residue(residue);

        let grid = build_carboxylate_grid(&structure, Some(7.4));
        let glu = structure.find_residue("A", 1, None).unwrap();
        let oe1 = glu.atom("OE1").unwrap().pos;

        let neighbors: Vec<_> = grid.neighbors(&oe1, 0.1).exact().collect();
        assert!(!neighbors.is_empty(), "GLU OE1 should be in COO⁻ grid");
    }

    #[test]
    fn coo_grid_excludes_glh_oxygens_when_protonated() {
        let mut residue = residue_from_template("GLU", StandardResidue::GLU, 1);
        residue.name = "GLH".into();
        let structure = structure_with_residue(residue);

        let grid = build_carboxylate_grid(&structure, Some(7.4));
        let glh = structure.find_residue("A", 1, None).unwrap();
        let oe1 = glh.atom("OE1").unwrap().pos;

        let neighbors: Vec<_> = grid.neighbors(&oe1, 0.1).exact().collect();
        assert!(
            neighbors.is_empty(),
            "GLH oxygens should NOT be in COO⁻ grid"
        );
    }

    #[test]
    fn coo_grid_includes_c_term_oxygens_at_neutral_ph() {
        let residue = c_terminal_residue(1);
        let structure = structure_with_residue(residue);

        let grid = build_carboxylate_grid(&structure, Some(7.4));
        let c_term = structure.find_residue("A", 1, None).unwrap();
        let oxt = c_term.atom("OXT").unwrap().pos;

        let neighbors: Vec<_> = grid.neighbors(&oxt, 0.1).exact().collect();
        assert!(
            !neighbors.is_empty(),
            "C-term OXT should be in COO⁻ grid at pH 7.4"
        );
    }

    #[test]
    fn coo_grid_excludes_c_term_oxygens_below_pka() {
        let residue = c_terminal_residue(1);
        let structure = structure_with_residue(residue);

        let grid = build_carboxylate_grid(&structure, Some(2.0));
        let c_term = structure.find_residue("A", 1, None).unwrap();
        let oxt = c_term.atom("OXT").unwrap().pos;

        let neighbors: Vec<_> = grid.neighbors(&oxt, 0.1).exact().collect();
        assert!(
            neighbors.is_empty(),
            "C-term OXT should NOT be in COO⁻ grid below pKa 3.1"
        );
    }

    #[test]
    fn coo_grid_uses_default_ph_when_target_ph_unset() {
        let residue = c_terminal_residue(1);
        let structure = structure_with_residue(residue);

        let grid = build_carboxylate_grid(&structure, None);
        let c_term = structure.find_residue("A", 1, None).unwrap();
        let oxt = c_term.atom("OXT").unwrap().pos;

        let neighbors: Vec<_> = grid.neighbors(&oxt, 0.1).exact().collect();
        assert!(
            !neighbors.is_empty(),
            "C-term OXT should be in COO⁻ grid with default pH 7.0"
        );
    }

    #[test]
    fn his_becomes_hip_below_pka_threshold() {
        let mut structure = his_isolated(1);
        let config = HydroConfig {
            target_ph: Some(5.5),
            his_salt_bridge_protonation: false,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(res.name, "HIP", "HIS should become HIP below pKa 6.0");
    }

    #[test]
    fn his_does_not_become_hip_above_pka_threshold() {
        let mut structure = his_isolated(1);
        let config = HydroConfig {
            target_ph: Some(7.4),
            his_salt_bridge_protonation: false,
            his_strategy: HisStrategy::DirectHIE,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(
            res.name, "HIE",
            "HIS should become HIE (not HIP) above pKa 6.0"
        );
    }

    #[test]
    fn his_becomes_hip_when_nd1_near_asp_carboxylate() {
        let mut structure = his_near_asp(1, 2, 3.5);
        let config = HydroConfig {
            target_ph: Some(7.4),
            his_salt_bridge_protonation: true,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let his = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(
            his.name, "HIP",
            "HIS near ASP should become HIP via salt bridge"
        );
    }

    #[test]
    fn his_becomes_hip_when_ne2_near_glu_carboxylate() {
        let mut structure = his_near_glu(1, 2, 3.5);
        let config = HydroConfig {
            target_ph: Some(7.4),
            his_salt_bridge_protonation: true,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let his = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(
            his.name, "HIP",
            "HIS near GLU should become HIP via salt bridge"
        );
    }

    #[test]
    fn his_becomes_hip_when_near_c_term_carboxylate() {
        let mut structure = his_near_c_term(1, 2, 3.5);
        let config = HydroConfig {
            target_ph: Some(7.4),
            his_salt_bridge_protonation: true,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let his = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(
            his.name, "HIP",
            "HIS near C-term COO⁻ should become HIP via salt bridge"
        );
    }

    #[test]
    fn his_remains_neutral_when_beyond_salt_bridge_threshold() {
        let mut structure = his_near_asp(1, 2, 10.0);
        let config = HydroConfig {
            target_ph: Some(7.4),
            his_salt_bridge_protonation: true,
            his_strategy: HisStrategy::DirectHIE,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let his = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(
            his.name, "HIE",
            "HIS beyond salt bridge distance should remain neutral"
        );
    }

    #[test]
    fn his_salt_bridge_detected_without_ph() {
        let mut structure = his_near_asp(1, 2, 3.5);
        let config = HydroConfig {
            target_ph: None,
            his_salt_bridge_protonation: true,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let his = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(
            his.name, "HIP",
            "salt bridge should be detected even without pH"
        );
    }

    #[test]
    fn his_salt_bridge_skipped_when_option_disabled() {
        let mut structure = his_near_asp(1, 2, 3.5);
        let config = HydroConfig {
            target_ph: Some(7.4),
            his_salt_bridge_protonation: false,
            his_strategy: HisStrategy::DirectHIE,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let his = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(
            his.name, "HIE",
            "salt bridge should be ignored when disabled"
        );
    }

    #[test]
    fn his_uses_direct_hid_strategy() {
        let mut structure = his_isolated(1);
        let config = HydroConfig {
            target_ph: Some(7.4),
            his_salt_bridge_protonation: false,
            his_strategy: HisStrategy::DirectHID,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(res.name, "HID", "DirectHID strategy should produce HID");
    }

    #[test]
    fn his_uses_direct_hie_strategy() {
        let mut structure = his_isolated(1);
        let config = HydroConfig {
            target_ph: Some(7.4),
            his_salt_bridge_protonation: false,
            his_strategy: HisStrategy::DirectHIE,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(res.name, "HIE", "DirectHIE strategy should produce HIE");
    }

    #[test]
    fn his_uses_random_strategy_produces_valid_tautomer() {
        let mut structure = his_isolated(1);
        let config = HydroConfig {
            target_ph: Some(7.4),
            his_salt_bridge_protonation: false,
            his_strategy: HisStrategy::Random,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert!(
            res.name == "HID" || res.name == "HIE",
            "random strategy should produce either HID or HIE, got {}",
            res.name
        );
    }

    #[test]
    fn his_uses_hb_network_strategy_defaults_to_hie_without_neighbors() {
        let mut structure = his_isolated(1);
        let config = HydroConfig {
            target_ph: Some(7.4),
            his_salt_bridge_protonation: false,
            his_strategy: HisStrategy::HbNetwork,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert!(
            res.name == "HID" || res.name == "HIE",
            "HbNetwork strategy should produce HID or HIE"
        );
    }

    #[test]
    fn his_preserves_hid_name_without_ph_and_no_salt_bridge() {
        let mut residue = residue_from_template("HID", StandardResidue::HIS, 1);
        residue.name = "HID".into();
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: None,
            his_salt_bridge_protonation: false,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(
            res.name, "HID",
            "HID should be preserved without pH and no salt bridge"
        );
    }

    #[test]
    fn his_preserves_hie_name_without_ph_and_no_salt_bridge() {
        let mut residue = residue_from_template("HIE", StandardResidue::HIS, 1);
        residue.name = "HIE".into();
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: None,
            his_salt_bridge_protonation: false,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(
            res.name, "HIE",
            "HIE should be preserved without pH and no salt bridge"
        );
    }

    #[test]
    fn his_preserves_hip_name_without_ph_and_no_salt_bridge() {
        let mut residue = residue_from_template("HIP", StandardResidue::HIS, 1);
        residue.name = "HIP".into();
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: None,
            his_salt_bridge_protonation: false,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(
            res.name, "HIP",
            "HIP should be preserved without pH and no salt bridge"
        );
    }

    #[test]
    fn his_preserves_name_without_ph_when_no_salt_bridge_found() {
        let mut residue = residue_from_template("HID", StandardResidue::HIS, 1);
        residue.name = "HID".into();
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: None,
            his_salt_bridge_protonation: true,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let res = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(
            res.name, "HID",
            "HID should be preserved when salt bridge not found"
        );
    }

    #[test]
    fn his_acidic_ph_overrides_salt_bridge_check() {
        let mut structure = his_near_asp(1, 2, 3.5);
        let config = HydroConfig {
            target_ph: Some(5.0),
            his_salt_bridge_protonation: true,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let his = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(
            his.name, "HIP",
            "acidic pH should result in HIP (priority 1)"
        );
    }

    #[test]
    fn his_salt_bridge_overrides_strategy_selection() {
        let mut structure = his_near_asp(1, 2, 3.5);
        let config = HydroConfig {
            target_ph: Some(7.4),
            his_salt_bridge_protonation: true,
            his_strategy: HisStrategy::DirectHIE,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let his = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(
            his.name, "HIP",
            "salt bridge should override strategy (priority 3 > 5)"
        );
    }

    #[test]
    fn his_strategy_applies_only_after_salt_bridge_miss() {
        let mut structure = his_near_asp(1, 2, 10.0);
        let config = HydroConfig {
            target_ph: Some(7.4),
            his_salt_bridge_protonation: true,
            his_strategy: HisStrategy::DirectHID,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let his = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(
            his.name, "HID",
            "strategy should apply when salt bridge not found"
        );
    }

    #[test]
    fn add_hydrogens_populates_internal_lysine_side_chain() {
        let mut residue = residue_from_template("LYS", StandardResidue::LYS, 10);
        residue.add_atom(Atom::new("FAKE", Element::H, Point::origin()));
        let mut structure = structure_with_residue(residue);

        add_hydrogens(&mut structure, &HydroConfig::default()).unwrap();

        let residue = structure.find_residue("A", 10, None).unwrap();
        assert!(residue.has_atom("HZ1"));
        assert!(residue.has_atom("HZ2"));
        assert!(residue.has_atom("HZ3"));
        assert!(
            !residue.has_atom("FAKE"),
            "existing H should be removed by default"
        );
    }

    #[test]
    fn add_hydrogens_keeps_existing_h_when_configured() {
        let mut residue = residue_from_template("ALA", StandardResidue::ALA, 1);
        residue.add_atom(Atom::new("HX", Element::H, Point::origin()));
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            remove_existing_h: false,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let residue = structure.find_residue("A", 1, None).unwrap();
        assert!(
            residue.has_atom("HX"),
            "existing H should be kept when configured"
        );
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
    fn close_cysteines_are_relabeled_to_cyx() {
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
        add_hydrogens(&mut structure, &HydroConfig::default()).unwrap();

        let res1 = structure.find_residue("A", 30, None).unwrap();
        let res2 = structure.find_residue("A", 31, None).unwrap();
        assert_eq!(res1.name, "CYX");
        assert_eq!(res2.name, "CYX");
    }

    #[test]
    fn close_cysteines_skip_hg_hydrogen() {
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
        add_hydrogens(&mut structure, &HydroConfig::default()).unwrap();

        let res1 = structure.find_residue("A", 30, None).unwrap();
        let res2 = structure.find_residue("A", 31, None).unwrap();
        assert!(!res1.has_atom("HG"), "CYX should not have HG");
        assert!(!res2.has_atom("HG"), "CYX should not have HG");
    }

    #[test]
    fn distant_cysteines_remain_unchanged() {
        let cys1 = residue_from_template("CYS", StandardResidue::CYS, 30);
        let mut cys2 = residue_from_template("CYS", StandardResidue::CYS, 31);

        let offset = Vector3::new(20.0, 0.0, 0.0);
        for atom in cys2.iter_atoms_mut() {
            atom.translate_by(&offset);
        }

        let mut structure = structure_with_residues(vec![cys1, cys2]);
        let config = HydroConfig {
            target_ph: Some(7.4),
            ..HydroConfig::default()
        };
        add_hydrogens(&mut structure, &config).unwrap();

        let res1 = structure.find_residue("A", 30, None).unwrap();
        let res2 = structure.find_residue("A", 31, None).unwrap();
        assert_eq!(res1.name, "CYS", "distant CYS should remain CYS");
        assert_eq!(res2.name, "CYS", "distant CYS should remain CYS");
    }

    #[test]
    fn n_terminal_defaults_to_protonated_without_ph() {
        let residue = n_terminal_residue(40);
        let mut structure = structure_with_residue(residue);

        add_hydrogens(&mut structure, &HydroConfig::default()).unwrap();

        let residue = structure.find_residue("A", 40, None).unwrap();
        assert!(residue.has_atom("H1"));
        assert!(residue.has_atom("H2"));
        assert!(
            residue.has_atom("H3"),
            "N-term should have 3 H at default pH 7.0"
        );
    }

    #[test]
    fn n_terminal_remains_protonated_below_pka() {
        let residue = n_terminal_residue(40);
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: Some(7.0),
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let residue = structure.find_residue("A", 40, None).unwrap();
        assert!(residue.has_atom("H1"));
        assert!(residue.has_atom("H2"));
        assert!(
            residue.has_atom("H3"),
            "N-term should have 3 H below pKa 8.0"
        );
    }

    #[test]
    fn n_terminal_deprotonates_above_pka() {
        let residue = n_terminal_residue(41);
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: Some(9.0),
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let residue = structure.find_residue("A", 41, None).unwrap();
        assert!(residue.has_atom("H1"));
        assert!(residue.has_atom("H2"));
        assert!(
            !residue.has_atom("H3"),
            "N-term should have only 2 H above pKa 8.0"
        );
    }

    #[test]
    fn n_terminal_h_has_tetrahedral_bond_lengths() {
        let residue = n_terminal_residue(98);
        let mut structure = structure_with_residue(residue);

        add_hydrogens(&mut structure, &HydroConfig::default()).unwrap();

        let residue = structure.find_residue("A", 98, None).unwrap();
        let n = residue.atom("N").expect("N").pos;
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
    }

    #[test]
    fn n_terminal_h_has_tetrahedral_bond_angles() {
        let residue = n_terminal_residue(98);
        let mut structure = structure_with_residue(residue);

        add_hydrogens(&mut structure, &HydroConfig::default()).unwrap();

        let residue = structure.find_residue("A", 98, None).unwrap();
        let n = residue.atom("N").expect("N").pos;
        let ca = residue.atom("CA").expect("CA").pos;
        let h1 = residue.atom("H1").expect("H1").pos;
        let h2 = residue.atom("H2").expect("H2").pos;
        let h3 = residue.atom("H3").expect("H3").pos;

        let ca_n_h1_angle = angle_deg(ca, n, h1);
        let ca_n_h2_angle = angle_deg(ca, n, h2);
        let ca_n_h3_angle = angle_deg(ca, n, h3);
        let h1_n_h2_angle = angle_deg(h1, n, h2);
        let h2_n_h3_angle = angle_deg(h2, n, h3);
        let h1_n_h3_angle = angle_deg(h1, n, h3);

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
    fn c_terminal_defaults_to_deprotonated_without_ph() {
        let residue = c_terminal_residue(51);
        let mut structure = structure_with_residue(residue);

        add_hydrogens(&mut structure, &HydroConfig::default()).unwrap();

        let residue = structure.find_residue("A", 51, None).unwrap();
        assert!(
            !residue.has_atom("HOXT"),
            "C-term should be deprotonated at default pH 7.0"
        );
    }

    #[test]
    fn c_terminal_protonates_under_acidic_ph() {
        let residue = c_terminal_residue(50);
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: Some(2.5),
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let residue = structure.find_residue("A", 50, None).unwrap();
        assert!(
            residue.has_atom("HOXT"),
            "C-term should be protonated below pKa 3.1"
        );
    }

    #[test]
    fn c_terminal_remains_deprotonated_at_physiological_ph() {
        let residue = c_terminal_residue(51);
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: Some(7.4),
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let residue = structure.find_residue("A", 51, None).unwrap();
        assert!(
            !residue.has_atom("HOXT"),
            "C-term should be deprotonated at pH 7.4"
        );
    }

    #[test]
    fn c_terminal_hoxt_has_tetrahedral_bond_length() {
        let residue = c_terminal_residue(99);
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: Some(2.0),
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let residue = structure.find_residue("A", 99, None).unwrap();
        let oxt = residue.atom("OXT").expect("OXT").pos;
        let hoxt = residue.atom("HOXT").expect("HOXT").pos;

        let oxt_hoxt_dist = distance(oxt, hoxt);
        assert!(
            (oxt_hoxt_dist - COOH_BOND_LENGTH).abs() < 0.1,
            "OXT-HOXT distance {oxt_hoxt_dist:.3} should be ~{COOH_BOND_LENGTH} Å"
        );
    }

    #[test]
    fn c_terminal_hoxt_has_tetrahedral_bond_angle() {
        let residue = c_terminal_residue(99);
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: Some(2.0),
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let residue = structure.find_residue("A", 99, None).unwrap();
        let c = residue.atom("C").expect("C").pos;
        let oxt = residue.atom("OXT").expect("OXT").pos;
        let hoxt = residue.atom("HOXT").expect("HOXT").pos;

        let c_oxt_hoxt_angle = angle_deg(c, oxt, hoxt);
        assert!(
            (c_oxt_hoxt_angle - SP3_ANGLE).abs() < 5.0,
            "C-OXT-HOXT angle {c_oxt_hoxt_angle:.1}° should be ~{SP3_ANGLE}°"
        );
    }

    #[test]
    fn five_prime_phosphate_deprotonated_at_physiological_ph() {
        let residue = five_prime_residue_with_phosphate(60);
        let mut structure = structure_with_residue(residue);

        add_hydrogens(&mut structure, &HydroConfig::default()).unwrap();

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
        let config = HydroConfig {
            target_ph: Some(5.5),
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let residue = structure.find_residue("A", 61, None).unwrap();
        assert!(residue.has_atom("OP3"), "OP3 should remain");
        assert!(
            residue.has_atom("HOP3"),
            "HOP3 should be added below pKa 6.5"
        );
    }

    #[test]
    fn five_prime_without_phosphate_gets_ho5() {
        let residue = five_prime_residue_without_phosphate(62);
        let mut structure = structure_with_residue(residue);

        add_hydrogens(&mut structure, &HydroConfig::default()).unwrap();

        let residue = structure.find_residue("A", 62, None).unwrap();
        assert!(
            residue.has_atom("HO5'"),
            "HO5' should be added for 5'-OH terminus"
        );
        assert!(!residue.has_atom("P"), "phosphorus should not exist");
    }

    #[test]
    fn five_prime_ho5_has_tetrahedral_geometry() {
        let residue = five_prime_residue_without_phosphate(80);
        let mut structure = structure_with_residue(residue);

        add_hydrogens(&mut structure, &HydroConfig::default()).unwrap();

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
    fn five_prime_phosphate_hop3_has_tetrahedral_geometry() {
        let residue = five_prime_residue_with_phosphate(82);
        let mut structure = structure_with_residue(residue);
        let config = HydroConfig {
            target_ph: Some(5.5),
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

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

    #[test]
    fn three_prime_nucleic_gets_ho3() {
        let residue = three_prime_residue(70);
        let mut structure = structure_with_residue(residue);

        add_hydrogens(&mut structure, &HydroConfig::default()).unwrap();

        let residue = structure.find_residue("A", 70, None).unwrap();
        assert!(
            residue.has_atom("HO3'"),
            "HO3' should be added for 3' terminal"
        );
    }

    #[test]
    fn three_prime_ho3_has_tetrahedral_geometry() {
        let residue = three_prime_residue(81);
        let mut structure = structure_with_residue(residue);

        add_hydrogens(&mut structure, &HydroConfig::default()).unwrap();

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
    fn hydro_config_defaults_to_no_ph() {
        let config = HydroConfig::default();
        assert!(config.target_ph.is_none(), "default should have no pH");
    }

    #[test]
    fn hydro_config_defaults_to_remove_existing_h() {
        let config = HydroConfig::default();
        assert!(config.remove_existing_h, "default should remove existing H");
    }

    #[test]
    fn hydro_config_defaults_to_hb_network_strategy() {
        let config = HydroConfig::default();
        assert_eq!(
            config.his_strategy,
            HisStrategy::HbNetwork,
            "default should use HbNetwork"
        );
    }

    #[test]
    fn hydro_config_defaults_to_salt_bridge_enabled() {
        let config = HydroConfig::default();
        assert!(
            config.his_salt_bridge_protonation,
            "default should enable salt bridge detection"
        );
    }

    #[test]
    fn effective_terminal_ph_returns_target_when_set() {
        assert_eq!(effective_terminal_ph(Some(5.5)), 5.5);
        assert_eq!(effective_terminal_ph(Some(9.0)), 9.0);
    }

    #[test]
    fn effective_terminal_ph_returns_default_when_unset() {
        assert_eq!(effective_terminal_ph(None), DEFAULT_TERMINAL_PH);
    }

    #[test]
    fn c_terminus_is_deprotonated_at_and_above_pka() {
        assert!(c_terminus_is_deprotonated(Some(7.4)));
        assert!(c_terminus_is_deprotonated(Some(4.0)));
        assert!(c_terminus_is_deprotonated(Some(C_TERM_PKA)));
        assert!(c_terminus_is_deprotonated(None));
    }

    #[test]
    fn c_terminus_is_protonated_below_pka() {
        assert!(!c_terminus_is_deprotonated(Some(2.0)));
        assert!(!c_terminus_is_deprotonated(Some(3.0)));
        assert!(!c_terminus_is_deprotonated(Some(C_TERM_PKA - 0.1)));
    }

    #[test]
    fn n_terminus_is_protonated_at_and_below_pka() {
        assert!(n_term_is_protonated(Some(7.0)));
        assert!(n_term_is_protonated(Some(N_TERM_PKA)));
        assert!(n_term_is_protonated(None));
    }

    #[test]
    fn n_terminus_is_deprotonated_above_pka() {
        assert!(!n_term_is_protonated(Some(9.0)));
        assert!(!n_term_is_protonated(Some(N_TERM_PKA + 0.1)));
    }

    #[test]
    fn c_term_protonation_boundary_is_exclusive() {
        assert!(c_terminus_is_deprotonated(Some(C_TERM_PKA)));
        assert!(!c_term_is_protonated(Some(C_TERM_PKA)));
    }

    #[test]
    fn acceptor_grid_includes_nitrogen_atoms() {
        let residue = residue_from_template("ALA", StandardResidue::ALA, 1);
        let structure = structure_with_residue(residue);

        let grid = build_acceptor_grid(&structure);
        let ala = structure.find_residue("A", 1, None).unwrap();
        let n = ala.atom("N").unwrap().pos;

        let neighbors: Vec<_> = grid.neighbors(&n, 0.1).exact().collect();
        assert!(!neighbors.is_empty(), "N should be in acceptor grid");
    }

    #[test]
    fn acceptor_grid_includes_oxygen_atoms() {
        let residue = residue_from_template("ALA", StandardResidue::ALA, 1);
        let structure = structure_with_residue(residue);

        let grid = build_acceptor_grid(&structure);
        let ala = structure.find_residue("A", 1, None).unwrap();
        let o = ala.atom("O").unwrap().pos;

        let neighbors: Vec<_> = grid.neighbors(&o, 0.1).exact().collect();
        assert!(!neighbors.is_empty(), "O should be in acceptor grid");
    }

    #[test]
    fn acceptor_grid_excludes_carbon_atoms() {
        let residue = residue_from_template("ALA", StandardResidue::ALA, 1);
        let structure = structure_with_residue(residue);

        let grid = build_acceptor_grid(&structure);
        let ala = structure.find_residue("A", 1, None).unwrap();
        let ca = ala.atom("CA").unwrap().pos;

        let neighbors: Vec<_> = grid.neighbors(&ca, 0.1).exact().collect();
        assert!(neighbors.is_empty(), "CA should NOT be in acceptor grid");
    }

    #[test]
    fn full_pipeline_preserves_user_protonation_without_ph() {
        let mut his = residue_from_template("HID", StandardResidue::HIS, 1);
        his.name = "HID".into();
        let mut ash = residue_from_template("ASP", StandardResidue::ASP, 2);
        ash.name = "ASH".into();
        let mut structure = structure_with_residues(vec![his, ash]);

        let config = HydroConfig {
            target_ph: None,
            his_salt_bridge_protonation: false,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let his = structure.find_residue("A", 1, None).unwrap();
        let ash = structure.find_residue("A", 2, None).unwrap();
        assert_eq!(his.name, "HID", "user-specified HID should be preserved");
        assert_eq!(ash.name, "ASH", "user-specified ASH should be preserved");
    }

    #[test]
    fn full_pipeline_applies_all_pka_rules_with_ph() {
        let asp = residue_from_template("ASP", StandardResidue::ASP, 1);
        let glu = residue_from_template("GLU", StandardResidue::GLU, 2);
        let lys = residue_from_template("LYS", StandardResidue::LYS, 3);
        let mut structure = structure_with_residues(vec![asp, glu, lys]);

        let config = HydroConfig {
            target_ph: Some(7.4),
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let asp = structure.find_residue("A", 1, None).unwrap();
        let glu = structure.find_residue("A", 2, None).unwrap();
        let lys = structure.find_residue("A", 3, None).unwrap();
        assert_eq!(asp.name, "ASP", "ASP should remain at pH 7.4");
        assert_eq!(glu.name, "GLU", "GLU should remain at pH 7.4");
        assert_eq!(lys.name, "LYS", "LYS should remain at pH 7.4");
    }

    #[test]
    fn full_pipeline_detects_salt_bridge_and_converts_his() {
        let mut structure = his_near_asp(1, 2, 3.5);
        let config = HydroConfig {
            target_ph: Some(7.4),
            his_salt_bridge_protonation: true,
            ..HydroConfig::default()
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let his = structure.find_residue("A", 1, None).unwrap();
        let asp = structure.find_residue("A", 2, None).unwrap();
        assert_eq!(his.name, "HIP", "HIS should become HIP via salt bridge");
        assert_eq!(asp.name, "ASP", "ASP should remain deprotonated");
    }

    #[test]
    fn full_pipeline_with_all_options_disabled() {
        let mut his = residue_from_template("HID", StandardResidue::HIS, 1);
        his.name = "HID".into();
        let mut structure = structure_with_residue(his);

        let config = HydroConfig {
            target_ph: None,
            remove_existing_h: false,
            his_salt_bridge_protonation: false,
            his_strategy: HisStrategy::DirectHIE,
        };

        add_hydrogens(&mut structure, &config).unwrap();

        let his = structure.find_residue("A", 1, None).unwrap();
        assert_eq!(
            his.name, "HID",
            "HID should be preserved with all options disabled"
        );
    }
}
