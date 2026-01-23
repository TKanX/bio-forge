<h1><img src="https://raw.githubusercontent.com/TKanX/bio-forge/main/assets/logo.svg" alt="BioForge Logo" width="50em"> BioForge</h1>

**BioForge** is a pure-Rust toolkit for automated preparation of biological macromolecules. It reads experimental structures (PDB/mmCIF), reconciles them with high-quality residue templates, repairs missing atoms, assigns hydrogens and termini, builds topologies, and optionally solvates the system with water and ions—all without leaving the Rust type system.

## Highlights

- **Template-driven accuracy** – Curated TOML templates for standard amino acids, nucleotides, and water guarantee reproducible coordinates, charges, and bonding.
- **High performance** – Multithreaded processing (via rayon) handles million-atom systems in milliseconds; single-pass parsing, in-place mutation, and zero-copy serialization minimize overhead.
- **Rich structure model** – Lightweight `Atom`, `Residue`, `Chain`, and `Structure` types backed by `nalgebra` make geometric operations trivial.
- **Format interoperability** – Buffered readers/writers for PDB, mmCIF, and MOL2 plus error types that surface precise parsing diagnostics.
- **Preparation pipeline** – Cleaning, repairing, protonating, solvation, coordinate transforms, and topology reconstruction share a common `ops::Error` so workflows compose cleanly.
- **WebAssembly support** – Full-featured WASM bindings for modern JavaScript bundlers (Vite, webpack, Rollup); ideal for browser-based molecular viewers and web applications.
- **Rust-first ergonomics** – No FFI, no global mutable state beyond the lazily-loaded template store, and edition 2024 guarantees modern language features.

## Processing Pipeline

```
Load → Clean → Repair → Hydrogenate → Solvate → Topology → Write
```

1. **Load** – `io::read_pdb_structure` or `io::read_mmcif_structure` parses coordinates with `IoContext` alias resolution.
2. **Clean** – `ops::clean_structure` removes waters, ions, hetero residues, or arbitrary residue names via `CleanConfig`.
3. **Repair** – `ops::repair_structure` realigns residues to templates and rebuilds missing heavy atoms (OXT on C-termini, OP3 on 5'-phosphorylated nucleic acids).
4. **Hydrogenate** – `ops::add_hydrogens` infers protonation states (configurable pH and histidine strategy) and reconstructs hydrogens from template anchors.
5. **Solvate** – `ops::solvate_structure` creates a periodic box, packs water on a configurable lattice, and swaps molecules for ions to satisfy a target charge.
6. **Topology** – `ops::TopologyBuilder` emits bond connectivity with peptide-link detection, nucleic backbone connectivity, and disulfide heuristics.
7. **Write** – `io::write_pdb_structure` / `io::write_mmcif_structure` serialize the processed structure; `write_*_topology` helpers emit CONECT or `struct_conn` records.

## Quick Start

### For CLI Users

Install the latest BioForge CLI binary from the [releases page](https://github.com/TKanX/bio-forge/releases) or via `cargo`:

```bash
cargo install bio-forge
```

Once the `bioforge` binary is installed, you can repair a structure in a single step:

```bash
bioforge repair -i input.pdb -o repaired.pdb
```

Explore the complete preparation pipeline in the [user manual](MANUAL.md) and browse the [examples directory](https://github.com/TKanX/bio-forge/tree/main/examples) for runnable walkthroughs.

### For Library Developers (Rust)

BioForge is also available as a library crate. Add it to your `Cargo.toml` dependencies:

```toml
[dependencies]
bio-forge = "0.4.0"
```

#### Example: Preparing a PDB Structure

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

### For Library Developers (JavaScript/TypeScript)

Install via npm:

```bash
npm install bio-forge-wasm
```

Prepare a structure with the following code:

```typescript
import { Structure } from "bio-forge-wasm";

const pdb = await fetch("https://files.rcsb.org/view/1UBQ.pdb").then((r) =>
  r.text()
);
const structure = Structure.fromPdb(pdb);

structure.clean({ removeWater: true });
structure.repair();
structure.addHydrogens({ hisStrategy: "network" });

const topology = structure.toTopology();
console.log(`Bonds: ${topology.bondCount}`);
```

## Documentation

| Resource                                                          | Description                    |
| ----------------------------------------------------------------- | ------------------------------ |
| [CLI Manual](MANUAL.md)                                           | Command-line usage and options |
| [JS/TS API](https://github.com/TKanX/bio-forge/blob/main/API.md)  | WebAssembly bindings reference |
| [Rust API](https://docs.rs/bio-forge)                             | Library documentation          |
| [Architecture](ARCHITECTURE.md)                                   | Internal design and algorithms |
| [Examples](https://github.com/TKanX/bio-forge/tree/main/examples) | Runnable walkthroughs          |

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
