use anyhow::{Context, Result, bail};
use clap::Args;

use bio_forge::Structure;
use bio_forge::ops::{Anion, Cation, SolvateConfig, solvate_structure};

use crate::commands::{estimate_structure_charge, run_with_spinner};

/// Builds a solvent box populated with ions.
#[derive(Debug, Args)]
pub struct SolvateArgs {
    /// Margin (Å) to add around the solute before packing water.
    #[arg(long, default_value_t = 10.0)]
    pub margin: f64,
    /// Water grid spacing (Å).
    #[arg(long, default_value_t = 3.1)]
    pub spacing: f64,
    /// Positive ion species to consider when swapping waters.
    #[arg(long, value_name = "ELEMENT", default_value = "Na")]
    pub cation: String,
    /// Negative ion species to consider when swapping waters.
    #[arg(long, value_name = "ELEMENT", default_value = "Cl")]
    pub anion: String,
    /// Automatically neutralize the solute by targeting zero total charge.
    #[arg(long, conflicts_with = "target_charge")]
    pub neutralize: bool,
    /// Explicit total charge to reach after solvation.
    #[arg(
        long = "target-charge",
        value_name = "INT",
        conflicts_with = "neutralize"
    )]
    pub target_charge: Option<i32>,
    /// Random seed used for deterministic ion placement.
    #[arg(long, value_name = "INT")]
    pub seed: Option<u64>,
}

/// Applies the solvation routine with configured padding, ions, and RNG seed.
pub fn run(structure: &mut Structure, args: &SolvateArgs) -> Result<()> {
    run_with_spinner("Solvating structure", || {
        let mut config = SolvateConfig {
            margin: args.margin,
            water_spacing: args.spacing,
            cations: vec![parse_cation(&args.cation)?],
            anions: vec![parse_anion(&args.anion)?],
            rng_seed: args.seed,
            ..SolvateConfig::default()
        };

        let solute_charge = estimate_structure_charge(structure);
        config.target_charge = if let Some(target) = args.target_charge {
            target
        } else if args.neutralize {
            0
        } else {
            solute_charge
        };

        solvate_structure(structure, &config).context("Failed to solvate structure")
    })
}

fn parse_cation(value: &str) -> Result<Cation> {
    match value.trim().to_ascii_uppercase().as_str() {
        "NA" => Ok(Cation::Na),
        "K" => Ok(Cation::K),
        "MG" => Ok(Cation::Mg),
        "CA" => Ok(Cation::Ca),
        "LI" => Ok(Cation::Li),
        "ZN" => Ok(Cation::Zn),
        other => bail!(
            "Unsupported cation '{}'. Choose from Na, K, Mg, Ca, Li, Zn.",
            other
        ),
    }
}

fn parse_anion(value: &str) -> Result<Anion> {
    match value.trim().to_ascii_uppercase().as_str() {
        "CL" => Ok(Anion::Cl),
        "BR" => Ok(Anion::Br),
        "I" => Ok(Anion::I),
        "F" => Ok(Anion::F),
        other => bail!("Unsupported anion '{}'. Choose from Cl, Br, I, F.", other),
    }
}
