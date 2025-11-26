use anyhow::{Context, Result};
use clap::Args;

use bio_forge::Structure;
use bio_forge::ops::{CleanConfig, clean_structure};

use crate::commands::{build_name_set, run_with_spinner};

/// Removes solvent, ions, hydrogens, or custom residues from a structure.
#[derive(Debug, Default, Args)]
pub struct CleanArgs {
    /// Remove crystallographic water residues (HOH).
    #[arg(long)]
    pub water: bool,
    /// Remove ions from the structure.
    #[arg(long)]
    pub ions: bool,
    /// Strip all hydrogen atoms.
    #[arg(long)]
    pub hydrogens: bool,
    /// Remove hetero residues.
    #[arg(long)]
    pub hetero: bool,
    /// Residue names to keep regardless of other filters.
    #[arg(long = "keep", value_name = "RES_NAME")]
    pub keep: Vec<String>,
    /// Residue names to remove regardless of other filters.
    #[arg(long = "remove", value_name = "RES_NAME")]
    pub remove: Vec<String>,
}

/// Applies cleanup rules to the provided structure.
pub fn run(structure: &mut Structure, args: &CleanArgs) -> Result<()> {
    run_with_spinner("Cleaning structure", || {
        let config = CleanConfig {
            remove_water: args.water,
            remove_ions: args.ions,
            remove_hydrogens: args.hydrogens,
            remove_hetero: args.hetero,
            keep_residue_names: build_name_set(&args.keep),
            remove_residue_names: build_name_set(&args.remove),
        };

        clean_structure(structure, &config).context("Failed to clean structure")?;
        Ok(())
    })
}
