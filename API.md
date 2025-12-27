# BioForge JavaScript/TypeScript API

WebAssembly bindings for molecular structure manipulation in browser and Node.js environments.

## Installation

```bash
npm install bio-forge
```

## Quick Start

```typescript
import { Structure } from "bio-forge";

const pdb = await fetch("https://files.rcsb.org/view/1UBQ.pdb").then((r) =>
  r.text()
);
const structure = Structure.fromPdb(pdb);

structure.clean({ removeWater: true });
structure.repair();
structure.addHydrogens({ hisStrategy: "network" });

const topology = structure.toTopology();
console.log(`Bonds: ${topology.bondCount}`);

const output = topology.toPdb();
```

---

## Classes

### Structure

Main class for molecular structure manipulation.

#### Factory Methods

| Method           | Parameters            | Returns     | Description               |
| ---------------- | --------------------- | ----------- | ------------------------- |
| `fromPdb`        | `content: string`     | `Structure` | Parse PDB format string   |
| `fromPdbBytes`   | `content: Uint8Array` | `Structure` | Parse PDB format bytes    |
| `fromMmcif`      | `content: string`     | `Structure` | Parse mmCIF format string |
| `fromMmcifBytes` | `content: Uint8Array` | `Structure` | Parse mmCIF format bytes  |

#### Properties

| Property       | Type     | Description        |
| -------------- | -------- | ------------------ |
| `chainCount`   | `number` | Number of chains   |
| `residueCount` | `number` | Number of residues |
| `atomCount`    | `number` | Number of atoms    |

#### Instance Methods

| Method  | Parameters | Returns     | Description        |
| ------- | ---------- | ----------- | ------------------ |
| `clone` | —          | `Structure` | Create a deep copy |

#### Export Methods

| Method         | Parameters                                        | Returns      | Description               |
| -------------- | ------------------------------------------------- | ------------ | ------------------------- |
| `toPdb`        | —                                                 | `string`     | Serialize to PDB format   |
| `toPdbBytes`   | —                                                 | `Uint8Array` | Serialize to PDB bytes    |
| `toMmcif`      | —                                                 | `string`     | Serialize to mmCIF format |
| `toMmcifBytes` | —                                                 | `Uint8Array` | Serialize to mmCIF bytes  |
| `toTopology`   | `config?: TopologyConfig, templates?: Template[]` | `Topology`   | Build topology with bonds |

#### Mutation Methods

| Method         | Parameters               | Returns | Description                |
| -------------- | ------------------------ | ------- | -------------------------- |
| `clean`        | `config?: CleanConfig`   | `void`  | Remove unwanted components |
| `repair`       | —                        | `void`  | Reconstruct missing atoms  |
| `addHydrogens` | `config?: HydroConfig`   | `void`  | Add hydrogen atoms         |
| `solvate`      | `config?: SolvateConfig` | `void`  | Add water box and ions     |

#### Transform Methods

| Method           | Parameters                        | Returns | Description                              |
| ---------------- | --------------------------------- | ------- | ---------------------------------------- |
| `translate`      | `x: number, y: number, z: number` | `void`  | Translate by vector (Å)                  |
| `centerGeometry` | —                                 | `void`  | Center geometric centroid at origin      |
| `centerMass`     | —                                 | `void`  | Center of mass at origin                 |
| `rotateX`        | `radians: number`                 | `void`  | Rotate around X axis                     |
| `rotateY`        | `radians: number`                 | `void`  | Rotate around Y axis                     |
| `rotateZ`        | `radians: number`                 | `void`  | Rotate around Z axis                     |
| `rotateEuler`    | `x: number, y: number, z: number` | `void`  | Rotate using Euler angles (XYZ, radians) |

#### Query Methods

| Method | Parameters | Returns         | Description                  |
| ------ | ---------- | --------------- | ---------------------------- |
| `info` | —          | `StructureInfo` | Get comprehensive statistics |

---

### Topology

Structure with bond connectivity information.

#### Factory Methods

| Method          | Parameters                                                              | Returns    | Description          |
| --------------- | ----------------------------------------------------------------------- | ---------- | -------------------- |
| `fromStructure` | `structure: Structure, config?: TopologyConfig, templates?: Template[]` | `Topology` | Build from Structure |

#### Properties

| Property    | Type        | Description                 |
| ----------- | ----------- | --------------------------- |
| `structure` | `Structure` | Underlying structure (copy) |
| `bondCount` | `number`    | Number of bonds             |

#### Export Methods

| Method         | Parameters | Returns      | Description                     |
| -------------- | ---------- | ------------ | ------------------------------- |
| `toPdb`        | —          | `string`     | PDB with CONECT records         |
| `toPdbBytes`   | —          | `Uint8Array` | PDB with CONECT records (bytes) |
| `toMmcif`      | —          | `string`     | mmCIF with struct_conn          |
| `toMmcifBytes` | —          | `Uint8Array` | mmCIF with struct_conn (bytes)  |

---

### Template

Residue template for custom ligand topology.

#### Factory Methods

| Method          | Parameters            | Returns    | Description              |
| --------------- | --------------------- | ---------- | ------------------------ |
| `fromMol2`      | `content: string`     | `Template` | Parse MOL2 format string |
| `fromMol2Bytes` | `content: Uint8Array` | `Template` | Parse MOL2 format bytes  |

#### Properties

| Property | Type     | Description                  |
| -------- | -------- | ---------------------------- |
| `name`   | `string` | Template name (residue code) |

---

## Configuration Interfaces

### CleanConfig

| Property             | Type       | Default | Description                                            |
| -------------------- | ---------- | ------- | ------------------------------------------------------ |
| `removeWater`        | `boolean`  | `false` | Remove water molecules (HOH)                           |
| `removeIons`         | `boolean`  | `false` | Remove ion residues                                    |
| `removeHydrogens`    | `boolean`  | `false` | Remove hydrogen atoms                                  |
| `removeHetero`       | `boolean`  | `false` | Remove hetero residues (ligands)                       |
| `removeResidueNames` | `string[]` | `[]`    | Specific residue names to remove                       |
| `keepResidueNames`   | `string[]` | `[]`    | Specific residue names to keep (overrides other rules) |

### HydroConfig

| Property          | Type                                      | Default     | Description                             |
| ----------------- | ----------------------------------------- | ----------- | --------------------------------------- |
| `targetPh`        | `number \| undefined`                     | `undefined` | Target pH for protonation decisions     |
| `removeExistingH` | `boolean`                                 | `true`      | Remove existing hydrogens before adding |
| `hisStrategy`     | `'hid' \| 'hie' \| 'random' \| 'network'` | `'network'` | Histidine tautomer strategy             |

### SolvateConfig

| Property         | Type                                                 | Default     | Description                              |
| ---------------- | ---------------------------------------------------- | ----------- | ---------------------------------------- |
| `margin`         | `number`                                             | `10.0`      | Margin around solute (Å)                 |
| `waterSpacing`   | `number`                                             | `3.1`       | Water grid spacing (Å)                   |
| `vdwCutoff`      | `number`                                             | `2.4`       | Minimum solvent-solute distance (Å)      |
| `removeExisting` | `boolean`                                            | `true`      | Remove existing solvent before solvating |
| `cations`        | `Array<'Na' \| 'K' \| 'Mg' \| 'Ca' \| 'Li' \| 'Zn'>` | `['Na']`    | Cation species                           |
| `anions`         | `Array<'Cl' \| 'Br' \| 'I' \| 'F'>`                  | `['Cl']`    | Anion species                            |
| `targetCharge`   | `number`                                             | `0`         | Target net charge after solvation        |
| `rngSeed`        | `number \| undefined`                                | `undefined` | Random seed for reproducible placement   |

### TopologyConfig

| Property          | Type     | Default | Description               |
| ----------------- | -------- | ------- | ------------------------- |
| `disulfideCutoff` | `number` | `2.2`   | Disulfide bond cutoff (Å) |

---

## Data Structures

### StructureInfo

```typescript
interface StructureInfo {
  chainCount: number;
  residueCount: number;
  atomCount: number;
  boxLengths?: [number, number, number]; // [a, b, c] in Å
  boxAngles?: [number, number, number]; // [α, β, γ] in degrees
  chains: ChainInfo[];
}
```

### ChainInfo

```typescript
interface ChainInfo {
  id: string;
  residueCount: number;
  atomCount: number;
  /**
   * Polymer classifications present in this chain.
   * Multiple types can coexist (e.g., protein with bound ligand).
   * Empty array indicates an empty chain.
   */
  polymerTypes: Array<"protein" | "nucleic" | "solvent" | "hetero">;
}
```

---

## Usage Examples

### Basic Processing

```typescript
const structure = Structure.fromPdb(pdbContent);

structure.clean({ removeWater: true, removeIons: true });
structure.repair();
structure.addHydrogens({ hisStrategy: "hid", targetPh: 7.4 });

const output = structure.toPdb();
```

### Solvation

```typescript
const structure = Structure.fromPdb(pdbContent);
structure.clean({ removeWater: true });

structure.solvate({
  margin: 12.0,
  cations: ["K"],
  anions: ["Cl"],
  targetCharge: 0,
  rngSeed: 42,
});

const output = structure.toPdb();
```

### Topology with Custom Ligand

```typescript
const ligandMol2 = await fetch("ligand.mol2").then((r) => r.text());
const ligandTemplate = Template.fromMol2(ligandMol2);

const structure = Structure.fromPdb(pdbContent);
const topology = structure.toTopology({ disulfideCutoff: 2.5 }, [
  ligandTemplate,
]);

console.log(`Total bonds: ${topology.bondCount}`);
const output = topology.toPdb();
```

### Structure Analysis

```typescript
const structure = Structure.fromPdb(pdbContent);
const info = structure.info();

console.log(`Chains: ${info.chainCount}`);
console.log(`Residues: ${info.residueCount}`);
console.log(`Atoms: ${info.atomCount}`);

if (info.boxLengths) {
  console.log(`Box: ${info.boxLengths.join(" × ")} Å`);
}

for (const chain of info.chains) {
  const types =
    chain.polymerTypes.length > 0 ? chain.polymerTypes.join(", ") : "empty";
  console.log(`  Chain ${chain.id}: ${types} (${chain.residueCount} residues)`);
}
```

### Geometric Transformations

```typescript
const structure = Structure.fromPdb(pdbContent);

structure.centerGeometry();
structure.rotateZ(Math.PI / 2);
structure.translate(10.0, 0.0, 0.0);

const output = structure.toPdb();
```

### Working with Binary Data

```typescript
const response = await fetch("structure.pdb");
const buffer = await response.arrayBuffer();
const bytes = new Uint8Array(buffer);

const structure = Structure.fromPdbBytes(bytes);
const outputBytes = structure.toPdbBytes();

// Download file
const blob = new Blob([outputBytes], { type: "chemical/x-pdb" });
const url = URL.createObjectURL(blob);
```

---

## Memory Management

All classes (`Structure`, `Topology`, `Template`) expose a `free()` method to manually release WASM memory.

```typescript
const structure = Structure.fromPdb(pdbContent);
// ... use structure ...
structure.free(); // Explicitly release memory
```

**When to call `free()`:**

- Long-running applications processing many structures
- Memory-constrained environments
- After finished using large structures

**When `free()` is optional:**

- Short-lived scripts where the page/process will terminate
- Modern browsers with adequate memory headroom

> **Warning:** Accessing an object after calling `free()` will throw an error.

---

## Platform Support

The WASM package is built with `wasm-pack --target bundler` and optimized for use with modern JavaScript bundlers.

| Platform          | Support                          |
| ----------------- | -------------------------------- |
| Browser (Bundler) | ✅ Vite, webpack, Rollup, Parcel |
| Node.js (Bundler) | ✅ Via bundler with WASM loader  |

> **Note:** Direct script-tag usage and native Node.js/Deno require rebuilding with `--target web` or `--target nodejs`.
