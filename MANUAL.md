# BioForge CLI Manual

BioForge ships with a composable command-line interface that mirrors the library pipeline: ingest structures, clean and repair them, add hydrogens, solvate, transform, and emit both coordinates and explicit topologies. This document explains the command surface in detail and provides examples so you can adopt the CLI confidently in scripts or interactive shells.

---

## Getting Started

### Installation

#### Option 1: Precompiled Binaries

Download the latest release for your platform from the [GitHub Releases](https://github.com/TKanX/bio-forge/releases).

#### Option 2: Cargo Install (from source)

```bash
cargo install bioforge
```

### Input & Output Basics

- **Formats** – BioForge understands PDB (`.pdb`, `.ent`) and mmCIF (`.cif`, `.mmcif`). Formats are auto-detected from file extensions, but you can override them with `--format` and `--out-format`.
- **Streaming** – Every subcommand can read from stdin (`-i` omitted) and write to stdout (`-o` omitted). Non-interactive safeguards prevent dumping structured data straight into a terminal; either redirect to a file or pipe into another command.
- **Context sharing** – Subcommands accept the same IO flags so you can combine them consistently.

### Global Flags

| Flag                        | Description                                                                    | Default                   |
| --------------------------- | ------------------------------------------------------------------------------ | ------------------------- |
| `-i, --input <FILE>`        | Structure to read. When absent, stdin is used.                                 | stdin                     |
| `-o, --output <FILE>`       | Destination for the resulting structure/topology. When absent, stdout is used. | stdout                    |
| `--format <pdb\|mmcif>`     | Force input parsing format.                                                    | Auto (or PDB when stdin)  |
| `--out-format <pdb\|mmcif>` | Force output serialization format.                                             | Auto (or PDB when stdout) |

---

## Command Reference

Each subcommand focuses on a single stage. Combine them to produce richer pipelines.

### `info` – Analyze the structure without mutation

```bash
bioforge info -i prepared.pdb
```

- Computes per-chain statistics (residue count, atom count, polymer class).
- Reports unit cell vectors/angles when available.
- Estimates total charge using residue templates and common ions.

### `clean` – Remove solvent, ions, hydrogens, or custom residues

```bash
bioforge clean -i raw.pdb -o cleaned.pdb --water --ions --remove NAG --keep LIG
```

Options:

| Flag             | Purpose                                                            |
| ---------------- | ------------------------------------------------------------------ |
| `--water`        | Drop crystallographic water (`HOH`).                               |
| `--ions`         | Remove metal and monatomic ions.                                   |
| `--hydrogens`    | Strip all hydrogen atoms.                                          |
| `--hetero`       | Drop hetero residues.                                              |
| `--keep <RES>`   | Protect specific residues from removal (may repeat).               |
| `--remove <RES>` | Forcibly remove residues regardless of other filters (may repeat). |

### `repair` – Rebuild missing atoms and termini

```bash
bioforge repair -i cleaned.pdb -o repaired.pdb
```

- Aligns each standard residue to its template and fills in missing heavy atoms, including peptide termini.
- Ideal immediately after `clean` to ensure the structure is chemically complete before protonation.

### `hydro` – Add hydrogens with titration awareness

```bash
bioforge hydro -i repaired.pdb -o protonated.pdb --ph 7.0 --his network
```

Options:

| Flag                                | Purpose                                                                   |
| ----------------------------------- | ------------------------------------------------------------------------- |
| `--ph <value>`                      | Target pH for protonation decisions.                                      |
| `--no-strip`                        | Keep existing hydrogens instead of stripping before rebuild.              |
| `--his <hid\|hie\|random\|network>` | Histidine tautomer strategy. Defaults to `network` (hydrogen-bond-aware). |

### `solvate` – Build a solvent box and add ions

```bash
bioforge solvate -i protonated.pdb -o solvated.pdb --margin 12 --spacing 3.0 --cation Na --anion Cl --neutralize --seed 42
```

Options:

| Flag                    | Purpose                                                          |
| ----------------------- | ---------------------------------------------------------------- |
| `--margin <Å>`          | Padding around the solute before packing waters (default 10 Å).  |
| `--spacing <Å>`         | Lattice spacing for initial water grid (default 3.1 Å).          |
| `--cation <element>`    | Cation species swapped into the solvent (Na, K, Mg, Ca, Li, Zn). |
| `--anion <element>`     | Anion species (Cl, Br, I, F).                                    |
| `--neutralize`          | Target zero net charge by adding/removing ions.                  |
| `--target-charge <int>` | Explicit charge goal (conflicts with `--neutralize`).            |
| `--seed <int>`          | RNG seed for deterministic ion placement.                        |

### `transform` – Apply centering, rotation, and translation

```bash
bioforge transform -i solvated.pdb -o boxed.pdb --center --rotate-z 90 --translate 0,0,5
```

Options:

| Flag                                     | Purpose                                          |
| ---------------------------------------- | ------------------------------------------------ |
| `--center`                               | Move geometric center to the origin.             |
| `--center-mass`                          | Move center of mass to the origin.               |
| `--rotate-x/--rotate-y/--rotate-z <deg>` | Rotate around axes (applied in X → Y → Z order). |
| `--translate <x,y,z>`                    | Translate by Cartesian vector (Å).               |

### `topology` – Emit bonded connectivity

```bash
bioforge topology -i boxed.pdb -o boxed-topology.pdb --out-format pdb --ss-cutoff 2.1
```

- Builds a `Topology` object using peptide, nucleic, and disulfide heuristics.
- Outputs either CONECT records (PDB) or `_struct_conn` categories (mmCIF).

Option:

| Flag                       | Purpose                                                                 |
| -------------------------- | ----------------------------------------------------------------------- |
| `--ss-cutoff <Å>`          | Maximum S–S distance used to infer disulfide bonds (default 2.2 Å).     |
| `--hetero-template <FILE>` | Include a Tripos MOL2 ligand template for hetero residues (repeatable). |

#### MOL2 template requirements

- **Molecule name in `@<TRIPOS>MOLECULE` must match the hetero residue name in the structure.** The builder uses this name to locate the correct template.
- **Atom labels in the MOL2 file must match the atom names in the structure.** Any mismatch will raise a topology atom missing error during bond graph generation.
- **Atom names within a single MOL2 file must be unique.** The parser enforces this and rejects duplicates to prevent ambiguous bonding.
- **Residue-internal atom names must also be unique within the structure for the same residue.** Duplicates in the coordinates will be rejected earlier in the pipeline.

---

## Workflow Examples

Please refer to the [examples directory](https://github.com/TKanX/bio-forge/tree/main/examples) for end-to-end usage scenarios demonstrating both single-command-per-stage and streaming pipeline approaches.

---

## Tips & Best Practices

- **Clarity over mutation** – `info` never mutates the structure, so you can insert it anywhere to inspect intermediate states.
- **Format overrides** – When piping between commands with mismatched defaults, specify `--format`/`--out-format` explicitly to avoid accidental PDB/mmCIF flips.
- **Determinism** – Provide `--seed` to `solvate` whenever you need reproducible water/ion placement.
- **Performance** – Piping avoids temporary files, but large systems may benefit from writing intermediate snapshots for debugging.

With these commands, you can automate structure preparation pipelines entirely from the terminal while reusing the same algorithms that power BioForge's Rust API.
