use anyhow::{Context, Result};
use clap::Args;

use bio_forge::ops::TopologyBuilder;
use bio_forge::{Structure, Topology};

use crate::commands::run_with_spinner;

/// Builds a bonded topology from the prepared structure.
#[derive(Debug, Args)]
pub struct TopologyArgs {
    /// Disulfide bond cutoff distance (Å).
    #[arg(long = "ss-cutoff", default_value_t = 2.2)]
    pub ss_cutoff: f64,
}

/// Generates a topology that can be written with CONECT/_struct_conn records.
pub fn run(structure: Structure, args: &TopologyArgs) -> Result<Topology> {
    run_with_spinner("Building topology", || {
        TopologyBuilder::new()
            .disulfide_cutoff(args.ss_cutoff)
            .build(structure)
            .context("Failed to build topology")
    })
}
