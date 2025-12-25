use std::{
    fmt,
    io::{self, Write},
};

use anyhow::{Context, Result};
use clap::Args;
use nalgebra::Vector3;
use prettytable::{Table, format, row};

use bio_forge::{Chain, ResidueCategory, StandardResidue, Structure};

use crate::commands::{estimate_structure_charge, run_with_spinner};

/// Report-only command that inspects a structure and mirrors the input stream.
#[derive(Debug, Default, Args)]
pub struct InfoArgs {}

/// Computes and prints structure statistics without mutating the structure.
pub fn run(structure: &Structure, _args: &InfoArgs) -> Result<()> {
    let (chain_reports, box_metrics, total_charge) =
        run_with_spinner("Analyzing structure", || {
            let reports = collect_chain_reports(structure);
            let box_metrics = calculate_box_metrics(structure);
            let charge = estimate_structure_charge(structure);
            Ok((reports, box_metrics, charge))
        })?;

    print_tables(&chain_reports, box_metrics.as_ref(), total_charge)?;
    Ok(())
}

fn collect_chain_reports(structure: &Structure) -> Vec<ChainReport> {
    structure
        .iter_chains()
        .map(|chain| ChainReport {
            id: chain.id.to_string(),
            residues: chain.residue_count(),
            atoms: chain.iter_atoms().count(),
            polymer: classify_chain(chain).to_string(),
        })
        .collect()
}

fn classify_chain(chain: &Chain) -> PolymerType {
    if chain.is_empty() {
        return PolymerType::Empty;
    }

    let mut protein = false;
    let mut nucleic = false;
    let mut solvent = false;
    let mut other = false;

    const LABEL_PROTEIN: &str = "Protein";
    const LABEL_NUCLEIC: &str = "Nucleic";
    const LABEL_SOLVENT: &str = "Solvent";
    const LABEL_OTHER: &str = "Other";

    for residue in chain.iter_residues() {
        match residue.category {
            ResidueCategory::Standard => {
                if let Some(std) = residue.standard_name {
                    if std.is_protein() {
                        protein = true;
                    } else if std.is_nucleic() {
                        nucleic = true;
                    } else if std == StandardResidue::HOH {
                        solvent = true;
                    } else {
                        other = true;
                    }
                } else {
                    other = true;
                }
            }
            ResidueCategory::Ion => other = true,
            ResidueCategory::Hetero => other = true,
        }
    }

    let mut components: Vec<&'static str> = Vec::new();
    if protein {
        components.push(LABEL_PROTEIN);
    }
    if nucleic {
        components.push(LABEL_NUCLEIC);
    }
    if solvent {
        components.push(LABEL_SOLVENT);
    }
    if other {
        components.push(LABEL_OTHER);
    }

    match components.as_slice() {
        [LABEL_PROTEIN] => PolymerType::Protein,
        [LABEL_NUCLEIC] => PolymerType::Nucleic,
        [LABEL_SOLVENT] => PolymerType::Solvent,
        [] => PolymerType::Mixed(vec![LABEL_OTHER]),
        _ => PolymerType::Mixed(components),
    }
}

fn print_tables(
    reports: &[ChainReport],
    box_metrics: Option<&BoxMetrics>,
    total_charge: i32,
) -> Result<()> {
    let mut stderr = io::stderr().lock();

    print_boxed_label(&mut stderr, "BioForge Structure Report")?;
    writeln!(&mut stderr)?;

    let mut chain_table = Table::new();
    print_boxed_label(&mut stderr, "Chain Breakdown")?;
    chain_table.set_format(*format::consts::FORMAT_BOX_CHARS);
    chain_table.set_titles(row!["Chain", "Residues", "Atoms", "Polymer Type"]);
    for report in reports {
        chain_table.add_row(row![
            report.id,
            report.residues,
            report.atoms,
            report.polymer
        ]);
    }
    chain_table
        .print(&mut stderr)
        .context("Failed to render chain summary")?;
    writeln!(&mut stderr)?;

    let mut summary_table = Table::new();
    print_boxed_label(&mut stderr, "Structure Summary")?;
    summary_table.set_format(*format::consts::FORMAT_BOX_CHARS);
    summary_table.set_titles(row!["Metric", "Value"]);

    if let Some(metrics) = box_metrics {
        summary_table.add_row(row![
            "Box Lengths (Å)",
            format!(
                "a = {length_a:.2}, b = {length_b:.2}, c = {length_c:.2}",
                length_a = metrics.length_a,
                length_b = metrics.length_b,
                length_c = metrics.length_c
            )
        ]);
        summary_table.add_row(row![
            "Box Angles (°)",
            format!(
                "α = {alpha:.2}, β = {beta:.2}, γ = {gamma:.2}",
                alpha = metrics.alpha,
                beta = metrics.beta,
                gamma = metrics.gamma
            )
        ]);
    } else {
        summary_table.add_row(row!["Box", "Not specified"]);
    }

    summary_table.add_row(row!["Estimated Total Charge", total_charge]);
    summary_table
        .print(&mut stderr)
        .context("Failed to render structure summary")?;

    Ok(())
}

fn print_boxed_label<W: Write>(writer: &mut W, title: &str) -> io::Result<()> {
    let inner = format!(" {title} ");
    let width = inner.chars().count();
    writeln!(writer, "╭{}╮", "─".repeat(width))?;
    writeln!(writer, "│{}│", inner)?;
    writeln!(writer, "╰{}╯", "─".repeat(width))?;
    Ok(())
}

#[derive(Debug)]
struct ChainReport {
    id: String,
    residues: usize,
    atoms: usize,
    polymer: String,
}

#[derive(Debug)]
enum PolymerType {
    Protein,
    Nucleic,
    Solvent,
    Mixed(Vec<&'static str>),
    Empty,
}

impl fmt::Display for PolymerType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PolymerType::Protein => write!(f, "Protein"),
            PolymerType::Nucleic => write!(f, "Nucleic"),
            PolymerType::Solvent => write!(f, "Solvent"),
            PolymerType::Mixed(components) => {
                if components.is_empty() {
                    write!(f, "Mixed")
                } else {
                    write!(f, "Mixed ({})", components.join(", "))
                }
            }
            PolymerType::Empty => write!(f, "Empty"),
        }
    }
}

#[derive(Debug)]
struct BoxMetrics {
    length_a: f64,
    length_b: f64,
    length_c: f64,
    alpha: f64,
    beta: f64,
    gamma: f64,
}

fn calculate_box_metrics(structure: &Structure) -> Option<BoxMetrics> {
    let vectors = structure.box_vectors?;
    let a = Vector3::from(vectors[0]);
    let b = Vector3::from(vectors[1]);
    let c = Vector3::from(vectors[2]);

    let length_a = a.norm();
    let length_b = b.norm();
    let length_c = c.norm();

    if length_a.abs() < f64::EPSILON
        || length_b.abs() < f64::EPSILON
        || length_c.abs() < f64::EPSILON
    {
        return None;
    }

    let alpha = angle_degrees(b, c);
    let beta = angle_degrees(a, c);
    let gamma = angle_degrees(a, b);

    Some(BoxMetrics {
        length_a,
        length_b,
        length_c,
        alpha,
        beta,
        gamma,
    })
}

fn angle_degrees(v1: Vector3<f64>, v2: Vector3<f64>) -> f64 {
    let denom = v1.norm() * v2.norm();
    if denom.abs() < f64::EPSILON {
        return 0.0;
    }

    let cos_angle = (v1.dot(&v2) / denom).clamp(-1.0, 1.0);
    cos_angle.acos().to_degrees()
}
