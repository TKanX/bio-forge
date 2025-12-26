//! Deserializable representations of residue templates shipped with `bio-forge`.
//!
//! The schema matches the TOML documents in `templates/` and is only used during startup
//! to populate the in-memory template store. All structs deny unknown fields so that the
//! template corpus fails fast when typos or unsupported keys are introduced.

use crate::model::types::{BondOrder, Element, StandardResidue};
use serde::Deserialize;

/// Top-level template document describing metadata, heavy atoms, hydrogens, and bonds.
#[derive(Debug, Deserialize, Clone)]
#[serde(deny_unknown_fields)]
pub struct ResidueTemplateFile {
    /// Human-readable metadata shared with the rest of the crate.
    pub info: TemplateInfo,
    /// Heavy atoms included in the residue definition.
    #[serde(default)]
    pub atoms: Vec<TemplateHeavyAtom>,
    /// Hydrogen placement hints that will be materialized by the hydrogenation pass.
    #[serde(default)]
    pub hydrogens: Vec<TemplateHydrogen>,
    /// Bond connectivity between named atoms.
    #[serde(default)]
    pub bonds: Vec<TemplateBond>,
}

/// Core metadata describing the residue.
#[derive(Debug, Deserialize, Clone)]
#[serde(deny_unknown_fields)]
pub struct TemplateInfo {
    /// Internal template identifier (matches filename).
    pub name: String,
    /// Mapping to the canonical residue enum.
    pub standard_name: StandardResidue,
    /// Net integer charge for the residue.
    pub charge: i32,
}

/// Heavy atom definition with elemental identity and reference coordinates.
#[derive(Debug, Deserialize, Clone)]
#[serde(deny_unknown_fields)]
pub struct TemplateHeavyAtom {
    /// Atom name unique within the residue.
    pub name: String,
    /// Chemical element associated with the atom.
    pub element: Element,
    /// Reference coordinates (Å) used when seeding structures.
    pub pos: [f64; 3],
}

/// Hydrogen definition with anchor metadata for idealized placement.
#[derive(Debug, Deserialize, Clone)]
#[serde(deny_unknown_fields)]
pub struct TemplateHydrogen {
    /// Hydrogen atom name unique within the residue.
    pub name: String,
    /// Reference coordinates (Å) relative to the heavy-atom frame.
    pub pos: [f64; 3],
    /// List of heavy-atom anchors that determine bonding context.
    pub anchors: Vec<String>,
}

/// Bond between two named atoms and its order.
#[derive(Debug, Deserialize, Clone)]
#[serde(deny_unknown_fields)]
pub struct TemplateBond {
    /// First atom name.
    pub a1: String,
    /// Second atom name.
    pub a2: String,
    /// Bond order in canonical form.
    pub order: BondOrder,
}
