use anyhow::{Context, Result};
use clap::{Args, ValueEnum};

use bio_forge::Structure;
use bio_forge::ops::{HisStrategy, HydroConfig, add_hydrogens};

use crate::commands::run_with_spinner;

/// Adds hydrogens according to pH- and histidine-aware rules.
#[derive(Debug, Args)]
pub struct HydroArgs {
    /// Target pH used for protonation decisions.
    #[arg(long = "ph")]
    pub ph: Option<f64>,
    /// Preserve existing hydrogens instead of stripping before rebuilding.
    #[arg(long = "no-strip")]
    pub no_strip: bool,
    /// Strategy for neutral histidine tautomer assignment.
    #[arg(long = "his", value_enum, default_value = "network")]
    pub his: HistidineStrategy,
    /// Disable salt bridge detection for HIS → HIP conversion.
    /// By default, histidines near carboxylate groups (ASP⁻/GLU⁻/C-term COO⁻) become HIP.
    #[arg(long = "no-his-salt-bridge")]
    pub no_his_salt_bridge: bool,
}

/// CLI exposure of [`HisStrategy`].
#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum HistidineStrategy {
    #[value(name = "hid")]
    Hid,
    #[value(name = "hie")]
    Hie,
    #[value(name = "random")]
    Random,
    #[value(name = "network")]
    Network,
}

impl From<HistidineStrategy> for HisStrategy {
    fn from(value: HistidineStrategy) -> Self {
        match value {
            HistidineStrategy::Hid => HisStrategy::DirectHID,
            HistidineStrategy::Hie => HisStrategy::DirectHIE,
            HistidineStrategy::Random => HisStrategy::Random,
            HistidineStrategy::Network => HisStrategy::HbNetwork,
        }
    }
}

/// Hydrogenates the structure using the configured options.
pub fn run(structure: &mut Structure, args: &HydroArgs) -> Result<()> {
    run_with_spinner("Adding hydrogens", || {
        let config = HydroConfig {
            target_ph: args.ph,
            remove_existing_h: !args.no_strip,
            his_strategy: args.his.into(),
            his_salt_bridge_protonation: !args.no_his_salt_bridge,
        };

        add_hydrogens(structure, &config).context("Failed to add hydrogens")
    })
}
