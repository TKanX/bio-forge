use std::path::PathBuf;

use anyhow::Result;
use clap::{Parser, Subcommand};

mod commands;

use commands::{IoParameters, StructureFormat};
use commands::{clean, hydro, info, repair, solvate, topology, transform};

#[derive(Parser, Debug)]
#[command(
    name = "bioforge",
    about = "A command-line tool for the automated repair, preparation, and topology construction of biological macromolecules.",
    version,
    author,
    arg_required_else_help = true
)]
struct Cli {
    /// Input file path. When omitted, stdin is used.
    #[arg(short, long, value_name = "FILE", global = true)]
    input: Option<PathBuf>,
    /// Output file path. When omitted, stdout is used.
    #[arg(short, long, value_name = "FILE", global = true)]
    output: Option<PathBuf>,
    /// Force the input format (pdb or mmcif).
    #[arg(long = "format", value_enum, global = true)]
    input_format: Option<StructureFormat>,
    /// Force the output format (pdb or mmcif).
    #[arg(long = "out-format", value_enum, global = true)]
    output_format: Option<StructureFormat>,
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Inspect the structure without modifying the data stream.
    Info(info::InfoArgs),
    /// Remove solvent, ions, hydrogens, or selected residues.
    Clean(clean::CleanArgs),
    /// Rebuild missing atoms and termini using templates.
    Repair(repair::RepairArgs),
    /// Add hydrogens using titration-aware heuristics.
    Hydro(hydro::HydroArgs),
    /// Solvate and optionally ionize the system.
    Solvate(solvate::SolvateArgs),
    /// Apply rotations, translations, and centering transforms.
    Transform(transform::TransformArgs),
    /// Build an explicit bonding topology.
    Topology(topology::TopologyArgs),
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let io_params = IoParameters {
        input: cli.input.clone(),
        output: cli.output.clone(),
        input_format: cli.input_format,
        output_format: cli.output_format,
    };

    match cli.command {
        Command::Topology(args) => {
            commands::ensure_noninteractive_stdout("topology", &io_params)?;
            let (structure, _) = commands::load_input(&io_params)?;
            let topology = topology::run(structure, &args)?;
            commands::save_topology(&topology, &io_params)?;
        }
        Command::Info(args) => {
            let (structure, _) = commands::load_input(&io_params)?;
            info::run(&structure, &args)?;
            if !commands::interactive_stdout_requested(&io_params) {
                commands::save_output(&structure, &io_params)?;
            }
        }
        Command::Clean(args) => {
            commands::ensure_noninteractive_stdout("clean", &io_params)?;
            let (mut structure, _) = commands::load_input(&io_params)?;
            clean::run(&mut structure, &args)?;
            commands::save_output(&structure, &io_params)?;
        }
        Command::Repair(args) => {
            commands::ensure_noninteractive_stdout("repair", &io_params)?;
            let (mut structure, _) = commands::load_input(&io_params)?;
            repair::run(&mut structure, &args)?;
            commands::save_output(&structure, &io_params)?;
        }
        Command::Hydro(args) => {
            commands::ensure_noninteractive_stdout("hydro", &io_params)?;
            let (mut structure, _) = commands::load_input(&io_params)?;
            hydro::run(&mut structure, &args)?;
            commands::save_output(&structure, &io_params)?;
        }
        Command::Solvate(args) => {
            commands::ensure_noninteractive_stdout("solvate", &io_params)?;
            let (mut structure, _) = commands::load_input(&io_params)?;
            solvate::run(&mut structure, &args)?;
            commands::save_output(&structure, &io_params)?;
        }
        Command::Transform(args) => {
            commands::ensure_noninteractive_stdout("transform", &io_params)?;
            let (mut structure, _) = commands::load_input(&io_params)?;
            transform::run(&mut structure, &args)?;
            commands::save_output(&structure, &io_params)?;
        }
    }

    Ok(())
}
