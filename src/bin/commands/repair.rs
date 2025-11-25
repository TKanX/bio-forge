use anyhow::{Context, Result};
use clap::Args;

use bio_forge::Structure;
use bio_forge::ops::repair_structure;

use crate::commands::run_with_spinner;

/// Repairs standard residues by rebuilding missing heavy atoms and termini.
#[derive(Debug, Default, Args)]
pub struct RepairArgs {}

/// Invokes the repair pipeline on the provided structure.
pub fn run(structure: &mut Structure, _args: &RepairArgs) -> Result<()> {
    run_with_spinner("Repairing structure", || {
        repair_structure(structure).context("Failed to repair structure")
    })
}
