use std::collections::HashSet;
use std::fmt;
use std::fs::File;
use std::io::{self as stdio, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Duration;

use anyhow::{Context, Result, anyhow, bail};
use clap::ValueEnum;
use indicatif::{ProgressBar, ProgressStyle};
use is_terminal::IsTerminal;

use bio_forge::io::{
    IoContext, read_mmcif_structure, read_pdb_structure, write_mmcif_structure,
    write_mmcif_topology, write_pdb_structure, write_pdb_topology,
};
use bio_forge::templates;
use bio_forge::{ResidueCategory, Structure, Topology};

pub mod clean;
pub mod hydro;
pub mod info;
pub mod repair;
pub mod solvate;
pub mod topology;
pub mod transform;

/// Formats supported by the CLI when reading or writing structures.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum StructureFormat {
    /// Legacy PDB format.
    #[value(name = "pdb")]
    Pdb,
    /// mmCIF format.
    #[value(name = "mmcif")]
    Mmcif,
}

impl StructureFormat {
    /// Attempts to infer a format from a file path extension.
    pub fn from_path(path: &Path) -> Option<Self> {
        path.extension()
            .and_then(|ext| ext.to_str())
            .and_then(Self::from_extension)
    }

    fn from_extension(ext: &str) -> Option<Self> {
        match ext.to_ascii_lowercase().as_str() {
            "pdb" | "ent" => Some(Self::Pdb),
            "cif" | "mmcif" => Some(Self::Mmcif),
            _ => None,
        }
    }
}

impl fmt::Display for StructureFormat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            StructureFormat::Pdb => write!(f, "PDB"),
            StructureFormat::Mmcif => write!(f, "mmCIF"),
        }
    }
}

/// Aggregated IO parameters shared by every subcommand.
#[derive(Debug, Clone, Default)]
pub struct IoParameters {
    pub input: Option<PathBuf>,
    pub output: Option<PathBuf>,
    pub input_format: Option<StructureFormat>,
    pub output_format: Option<StructureFormat>,
}

/// Loads a structure from the configured input source.
pub fn load_input(params: &IoParameters) -> Result<(Structure, StructureFormat)> {
    let format = resolve_input_format(params)?;
    let io_context = IoContext::new_default();

    let structure = if let Some(path) = &params.input {
        let file = File::open(path)
            .with_context(|| format!("Failed to open input file {}", path.display()))?;
        let reader = BufReader::new(file);
        read_structure(reader, format, &io_context)
            .with_context(|| format!("Failed to parse {} input from {}", format, path.display()))?
    } else {
        let stdin = stdio::stdin();
        if stdin.is_terminal() {
            bail!(
                "No --input provided and stdin is a TTY. Provide -i/--input or pipe a structure into bioforge."
            );
        }
        let reader = BufReader::new(stdin.lock());
        read_structure(reader, format, &io_context)
            .with_context(|| format!("Failed to parse {} input from stdin", format))?
    };

    Ok((structure, format))
}

/// Saves a structure to the configured output destination.
pub fn save_output(structure: &Structure, params: &IoParameters) -> Result<()> {
    let format = resolve_output_format(params)?;
    write_structure(structure, params.output.as_deref(), format)
}

/// Saves a topology to the configured output destination.
pub fn save_topology(topology: &Topology, params: &IoParameters) -> Result<()> {
    let format = resolve_output_format(params)?;
    write_topology(topology, params.output.as_deref(), format)
}

fn resolve_input_format(params: &IoParameters) -> Result<StructureFormat> {
    if let Some(explicit) = params.input_format {
        Ok(explicit)
    } else if let Some(path) = &params.input {
        StructureFormat::from_path(path).ok_or_else(|| {
            anyhow!(
                "Unable to infer input format from '{}'. Please specify --format.",
                path.display()
            )
        })
    } else {
        Ok(StructureFormat::Pdb)
    }
}

fn resolve_output_format(params: &IoParameters) -> Result<StructureFormat> {
    if let Some(explicit) = params.output_format {
        return Ok(explicit);
    }

    if let Some(path) = &params.output {
        StructureFormat::from_path(path).ok_or_else(|| {
            anyhow!(
                "Unable to infer output format from '{}'. Please specify --out-format.",
                path.display()
            )
        })
    } else {
        Ok(StructureFormat::Pdb)
    }
}

fn read_structure<R: BufRead>(
    reader: R,
    format: StructureFormat,
    ctx: &IoContext,
) -> Result<Structure> {
    match format {
        StructureFormat::Pdb => read_pdb_structure(reader, ctx).map_err(anyhow::Error::new),
        StructureFormat::Mmcif => read_mmcif_structure(reader, ctx).map_err(anyhow::Error::new),
    }
}

fn write_structure(
    structure: &Structure,
    output: Option<&Path>,
    format: StructureFormat,
) -> Result<()> {
    match output {
        Some(path) => {
            let file = File::create(path)
                .with_context(|| format!("Failed to create output file {}", path.display()))?;
            let mut writer = BufWriter::new(file);
            write_structure_with_format(&mut writer, structure, format).with_context(|| {
                format!("Failed to write {} output to {}", format, path.display())
            })?;
            writer.flush().context("Failed to flush output writer")?
        }
        None => {
            let stdout = stdio::stdout();
            let handle = stdout.lock();
            let mut writer = BufWriter::new(handle);
            write_structure_with_format(&mut writer, structure, format)
                .with_context(|| format!("Failed to write {} output to stdout", format))?;
            writer.flush().context("Failed to flush stdout")?;
        }
    }
    Ok(())
}

fn write_topology(
    topology: &Topology,
    output: Option<&Path>,
    format: StructureFormat,
) -> Result<()> {
    match output {
        Some(path) => {
            let file = File::create(path)
                .with_context(|| format!("Failed to create output file {}", path.display()))?;
            let mut writer = BufWriter::new(file);
            write_topology_with_format(&mut writer, topology, format).with_context(|| {
                format!("Failed to write {} topology to {}", format, path.display())
            })?;
            writer.flush().context("Failed to flush output writer")?
        }
        None => {
            let stdout = stdio::stdout();
            let handle = stdout.lock();
            let mut writer = BufWriter::new(handle);
            write_topology_with_format(&mut writer, topology, format)
                .with_context(|| format!("Failed to write {} topology to stdout", format))?;
            writer.flush().context("Failed to flush stdout")?;
        }
    }
    Ok(())
}

fn write_structure_with_format<W: Write>(
    writer: &mut W,
    structure: &Structure,
    format: StructureFormat,
) -> Result<()> {
    match format {
        StructureFormat::Pdb => {
            write_pdb_structure(writer, structure).map_err(anyhow::Error::new)?
        }
        StructureFormat::Mmcif => {
            write_mmcif_structure(writer, structure).map_err(anyhow::Error::new)?
        }
    }
    Ok(())
}

fn write_topology_with_format<W: Write>(
    writer: &mut W,
    topology: &Topology,
    format: StructureFormat,
) -> Result<()> {
    match format {
        StructureFormat::Pdb => write_pdb_topology(writer, topology).map_err(anyhow::Error::new)?,
        StructureFormat::Mmcif => {
            write_mmcif_topology(writer, topology).map_err(anyhow::Error::new)?
        }
    }
    Ok(())
}

/// Wraps long-running operations with a spinner rendered to stderr.
pub fn run_with_spinner<T, F>(message: &str, work: F) -> Result<T>
where
    F: FnOnce() -> Result<T>,
{
    let spinner = ProgressBar::new_spinner();
    let style = ProgressStyle::with_template("{spinner:.green} {msg}")
        .unwrap_or_else(|_| ProgressStyle::default_spinner());
    spinner.set_style(style);
    spinner.enable_steady_tick(Duration::from_millis(80));
    spinner.set_message(message.to_string());

    let result = work();

    match &result {
        Ok(_) => spinner.finish_with_message(format!("{} ✓", message)),
        Err(_) => spinner.abandon_with_message(format!("{} ✗", message)),
    }

    result
}

/// Estimates the integer charge of a structure using residue templates and known ions.
pub fn estimate_structure_charge(structure: &Structure) -> i32 {
    let mut charge = 0;
    for chain in structure.iter_chains() {
        for residue in chain.iter_residues() {
            match residue.category {
                ResidueCategory::Standard => {
                    if let Some(template) = templates::get(&residue.name) {
                        charge += template.charge();
                    }
                }
                ResidueCategory::Ion => {
                    charge += match residue.name.as_str() {
                        "NA" | "K" | "LI" => 1,
                        "MG" | "CA" | "ZN" => 2,
                        "CL" | "BR" | "I" | "F" => -1,
                        _ => 0,
                    };
                }
                ResidueCategory::Hetero => {}
            }
        }
    }
    charge
}

/// Normalizes residue name lists to uppercase hash sets.
pub fn build_name_set(values: &[String]) -> HashSet<String> {
    values
        .iter()
        .map(|v| v.trim().to_ascii_uppercase())
        .collect()
}

/// Returns true when stdout is a TTY and no explicit output file was supplied.
pub fn interactive_stdout_requested(params: &IoParameters) -> bool {
    params.output.is_none() && stdio::stdout().is_terminal()
}

/// Ensures commands do not dump structured output directly into an interactive terminal.
pub fn ensure_noninteractive_stdout(command: &str, params: &IoParameters) -> Result<()> {
    if interactive_stdout_requested(params) {
        bail!(
            "Refusing to stream {command} results to an interactive terminal. Use -o/--output or pipe the command into a file."
        );
    }
    Ok(())
}
