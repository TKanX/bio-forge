//! Shared error types returned by the high-level operations modules.
//!
//! Every variant maps to a specific biological invariant violation (missing templates,
//! failed alignments, protonation issues, etc.) so downstream callers can display precise
//! remediation guidance.

use thiserror::Error;

/// Error conditions surfaced by the operations layer.
#[derive(Debug, Error)]
pub enum Error {
    /// Internal template lookup failed for a standard residue.
    #[error("internal template not found for standard residue '{res_name}'")]
    MissingInternalTemplate { res_name: String },

    /// Least-squares alignment between residue coordinates and template failed.
    #[error("alignment failed for residue '{res_name}' ({res_id}): {reason}")]
    AlignmentFailed {
        res_name: String,
        res_id: i32,
        reason: String,
    },

    /// Hydrogen addition could not proceed because a required anchor atom is absent.
    #[error(
        "cannot add hydrogens to residue '{res_name}' ({res_id}): missing anchor atom '{atom_name}'"
    )]
    IncompleteResidueForHydro {
        res_name: String,
        res_id: i32,
        atom_name: String,
    },

    /// Simulation bounding box cannot accommodate requested solvent parameters.
    #[error("simulation box is too small for the requested solvation parameters")]
    BoxTooSmall,

    /// Replacement of waters with ions could not reach the requested charge balance.
    #[error("ionization failed: {details}")]
    IonizationFailed { details: String },

    /// No heterogen template was available for the residue.
    #[error("missing hetero topology template for residue '{res_name}'")]
    MissingHeteroTemplate { res_name: String },

    /// Residue is missing a heavy atom mandated by the template topology.
    #[error(
        "topology mismatch: Residue '{res_name}' ({res_id}) is missing atom '{atom_name}' required by template"
    )]
    TopologyAtomMissing {
        res_name: String,
        res_id: i32,
        atom_name: String,
    },
}

impl Error {
    /// Helper for constructing an [`Error::AlignmentFailed`] variant.
    ///
    /// # Arguments
    ///
    /// * `res_name` - Residue name to include in the message.
    /// * `res_id` - PDB/author residue identifier.
    /// * `reason` - Free-form explanation of the failure.
    pub fn alignment_failed(
        res_name: impl Into<String>,
        res_id: i32,
        reason: impl Into<String>,
    ) -> Self {
        Self::AlignmentFailed {
            res_name: res_name.into(),
            res_id,
            reason: reason.into(),
        }
    }

    /// Helper for constructing an [`Error::IncompleteResidueForHydro`] variant.
    ///
    /// # Arguments
    ///
    /// * `res_name` - Residue label.
    /// * `res_id` - Residue identifier.
    /// * `atom_name` - Anchor atom that is missing.
    pub fn incomplete_for_hydro(
        res_name: impl Into<String>,
        res_id: i32,
        atom_name: impl Into<String>,
    ) -> Self {
        Self::IncompleteResidueForHydro {
            res_name: res_name.into(),
            res_id,
            atom_name: atom_name.into(),
        }
    }

    /// Helper for constructing an [`Error::TopologyAtomMissing`] variant.
    ///
    /// # Arguments
    ///
    /// * `res_name` - Residue name as reported to the user.
    /// * `res_id` - Residue identifier.
    /// * `atom_name` - The absent atom that triggered the mismatch.
    pub fn topology_atom_missing(
        res_name: impl Into<String>,
        res_id: i32,
        atom_name: impl Into<String>,
    ) -> Self {
        Self::TopologyAtomMissing {
            res_name: res_name.into(),
            res_id,
            atom_name: atom_name.into(),
        }
    }
}
