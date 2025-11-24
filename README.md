# BioForge

**BioForge** is a pure-Rust toolkit for automated preparation of biological macromolecules. It reads experimental structures (PDB/mmCIF), reconciles them with high-quality residue templates, repairs missing atoms, assigns hydrogens and termini, builds topologies, and optionally solvates the system with water and ions—all without leaving the Rust type system.

## Highlights

- **Template-driven accuracy** – curated TOML templates for standard amino acids, nucleotides, and water guarantee reproducible coordinates, charges, and bonding.
- **Rich structure model** – lightweight `Atom`, `Residue`, `Chain`, and `Structure` types backed by `nalgebra` make geometric operations trivial.
- **Format interoperability** – buffered readers/writers for PDB, mmCIF, and MOL2 plus error types that surface precise parsing diagnostics.
- **Preparation pipeline** – cleaning, repairing, protonating, solvation, coordinate transforms, and topology reconstruction share a common `ops::Error` so workflows compose cleanly.
- **Rust-first ergonomics** – no FFI, no global mutable state beyond the lazily-loaded template store, and edition 2024 guarantees modern language features.

## Processing Pipeline at a Glance

1. **Load** – `io::read_pdb_structure` or `io::read_mmcif_structure` parses coordinates with help from `IoContext` alias resolution.
2. **Clean** – `ops::clean_structure` removes waters, ions, hetero residues, or arbitrary residue names via `CleanConfig`.
3. **Repair** – `ops::repair_structure` realigns residues to their templates and rebuilds missing heavy atoms (including OXT on C-termini).
4. **Hydrogenate** – `ops::add_hydrogens` infers protonation states (configurable pH and histidine strategy) and reconstructs hydrogens from template anchors.
5. **Solvate/Ionize** – `ops::solvate_structure` creates a periodic box, packs water on a configurable lattice, and swaps molecules for ions to satisfy a target charge.
6. **Topology** – `ops::TopologyBuilder` replays template bond definitions, peptide-link detection, nucleic backbone connectivity, and disulfide heuristics to emit a `Topology` object.
7. **Write** – `io::write_pdb_structure` / `io::write_mmcif_structure` serialize the processed structure; the corresponding `write_*_topology` helpers emit CONECT or `struct_conn` records.

## Quick Start

BioForge is currently distributed as a library crate. Add it to your `Cargo.toml` dependencies:

```toml
[dependencies]
bio-forge = "0.1.0"
```

### Example: Preparing a PDB Structure

```rust
use std::{fs::File, io::{BufReader, BufWriter}};

use bio_forge::{
    io::{
        read_pdb_structure,
        write_pdb_structure,
        write_pdb_topology,
        IoContext,
    },
    ops::{
        add_hydrogens, clean_structure, repair_structure, solvate_structure,
        CleanConfig, HydroConfig, SolvateConfig, TopologyBuilder,
    },
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let ctx = IoContext::new_default();
    let input = BufReader::new(File::open("input.pdb")?);
    let mut structure = read_pdb_structure(input, &ctx)?;

    clean_structure(&mut structure, &CleanConfig::water_only())?;
    repair_structure(&mut structure)?;
    add_hydrogens(&mut structure, &HydroConfig::default())?;
    solvate_structure(&mut structure, &SolvateConfig::default())?;

    let topology = TopologyBuilder::new().build(structure.clone())?;

    write_pdb_structure(BufWriter::new(File::create("prepared.pdb")?), &structure)?;
    write_pdb_topology(BufWriter::new(File::create("prepared-topology.pdb")?), &topology)?;
    Ok(())
}
```

> Prefer mmCIF? Swap in `read_mmcif_structure` / `write_mmcif_structure`. Need to process ligands? Parse them via `io::read_mol2_template` and feed the resulting `Template` into `TopologyBuilder::add_hetero_template`.

## Documentation

- [API Documentation](https://docs.rs/bio-forge) – comprehensive reference for public types and functions.
- [Architecture Overview](ARCHITECTURE.md) – detailed explanation of the internal design and algorithms used in BioForge.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
