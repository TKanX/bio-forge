use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error("internal template not found for standard residue '{res_name}'")]
    MissingInternalTemplate { res_name: String },

    #[error("alignment failed for residue '{res_name}' ({res_id}): {reason}")]
    AlignmentFailed {
        res_name: String,
        res_id: i32,
        reason: String,
    },

    #[error(
        "cannot add hydrogens to residue '{res_name}' ({res_id}): missing anchor atom '{atom_name}'"
    )]
    IncompleteResidueForHydro {
        res_name: String,
        res_id: i32,
        atom_name: String,
    },

    #[error("simulation box is too small for the requested solvation parameters")]
    BoxTooSmall,

    #[error("ionization failed: {details}")]
    IonizationFailed { details: String },

    #[error("missing hetero topology template for residue '{res_name}'")]
    MissingHeteroTemplate { res_name: String },

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
