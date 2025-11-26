use std::collections::HashSet;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use anyhow::{Context, Result, bail};
use clap::Args;

use bio_forge::io::read_mol2_template;
use bio_forge::ops::TopologyBuilder;
use bio_forge::{Structure, Template, Topology};

use crate::commands::run_with_spinner;

/// Builds a bonded topology from the prepared structure.
#[derive(Debug, Args)]
pub struct TopologyArgs {
    /// Disulfide bond cutoff distance (Å).
    #[arg(long = "ss-cutoff", default_value_t = 2.2)]
    pub ss_cutoff: f64,
    /// Additional hetero-residue templates (Tripos MOL2) to satisfy ligands (repeatable).
    #[arg(long = "hetero-template", value_name = "FILE")]
    pub hetero_templates: Vec<PathBuf>,
}

/// Generates a topology that can be written with CONECT/_struct_conn records.
pub fn run(structure: Structure, args: &TopologyArgs) -> Result<Topology> {
    let hetero_templates = load_mol2_templates(&args.hetero_templates)?;
    let ss_cutoff = args.ss_cutoff;

    run_with_spinner("Building topology", move || {
        let builder = hetero_templates.into_iter().fold(
            TopologyBuilder::new().disulfide_cutoff(ss_cutoff),
            |builder, template| builder.add_hetero_template(template),
        );

        builder.build(structure).context("Failed to build topology")
    })
}

fn load_mol2_templates(paths: &[PathBuf]) -> Result<Vec<Template>> {
    let mut templates = Vec::new();
    let mut seen_names = HashSet::new();

    for path in paths {
        let file = File::open(path)
            .with_context(|| format!("Failed to open hetero template (MOL2) {}", path.display()))?;
        let reader = BufReader::new(file);
        let template = read_mol2_template(reader).with_context(|| {
            format!(
                "Failed to parse hetero template from {}. Ensure the file is a valid MOL2 ligand definition.",
                path.display()
            )
        })?;

        if !seen_names.insert(template.name.clone()) {
            bail!(
                "Duplicate hetero template '{}' provided. Supply each residue definition once.",
                template.name
            );
        }

        templates.push(template);
    }

    Ok(templates)
}
