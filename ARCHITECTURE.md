# BioForge Architecture

## 1. System Overview

BioForge is organized around five cooperating domains: data templates (`templates/` and `db`), structural models (`model`), IO codecs (`io`), molecular operations (`ops`), and public exports (`lib.rs`). Templates encode chemically faithful residue definitions, IO modules translate external files into internal structures, and ops mutate those structures into simulation-ready systems.

```mermaid
flowchart TD
    Templates --> DB
    DB --> Model
    IO --> Model
    Model --> Ops
    Ops --> Writers
    Writers --> Outputs
```

- **Templates** – static TOML descriptors for every supported residue or solvent species.
- **DB** – the template loader/store that exposes `TemplateView` references to the rest of the crate.
- **IO** – parsers and serializers (PDB, mmCIF, MOL2) that hydrate or persist `Structure` instances.
- **Model** – neutral data classes (`Atom`, `Residue`, `Chain`, `Structure`, `Topology`) shared across modules.
- **Ops** – transformation pipelines that clean, repair, protonate, solvate, and connect structures.
- **Writers** – emitters that translate the final structure/topology back to biochemical formats.

## 2. IO and Data Modules

The IO and data layers orchestrate template lookup, alias resolution, and structure instantiation before ops run.

```mermaid
flowchart LR
    Readers -->|alias map| IoCtx
    IoCtx -->|canonical name| ModelStruct
    Readers -->|atoms| ModelStruct
    ModelStruct -->|residue ids| Store
    Store -->|TemplateView| ModelStruct
    ModelStruct --> Writers
```

- **Readers** – streaming parsers (`io::pdb::reader`, `io::mmcif::reader`, `io::mol2::reader`) that emit `Structure` or `Template` values.
- **IoCtx** – `IoContext` maps aliases to canonical residue names and standard enums for consistent categorization.
- **ModelStruct** – instances of `Structure`, `Chain`, and `Residue` that accumulate atomic data and metadata.
- **Store** – singleton `db::store::DataStore` holding every pre-parsed template accessible via `TemplateView`.
- **Writers** – formatters that consume `Structure`/`Topology` pairs to write PDB or mmCIF output.

## 3. Ops Pipelines

Each ops submodule follows a focused pipeline. Flowcharts capture the order of decisions and mutations.

### 3.1 Clean Pipeline (`ops::clean`)

```mermaid
flowchart TD
    Start --> Config
    Config --> StripH
    StripH --> ResidueFilter
    ResidueFilter --> ChainPrune
    ChainPrune --> End
```

- **Start** – receives a mutable `Structure` plus `CleanConfig` options.
- **Config** – interprets flags (water-only, ion removal, name allow/deny lists).
- **StripH** – optionally removes hydrogen atoms residue-by-residue.
- **ResidueFilter** – retains or drops residues based on category, standard name, and explicit keep/remove sets.
- **ChainPrune** – prunes chains emptied by prior filtering to keep the structure minimal.
- **End** – yields a cleaned structure ready for downstream repair.

### 3.2 Repair Pipeline (`ops::repair`)

```mermaid
flowchart TD
    CollectTargets --> FetchTemplate
    FetchTemplate --> AlignPairs
    AlignPairs --> SVD
    SVD --> SynthesizeAtoms
    SynthesizeAtoms --> Cleanup
```

- **CollectTargets** – iterates standard residues to identify those eligible for repair.
- **FetchTemplate** – pulls `TemplateView` for each residue name, failing early if missing.
- **AlignPairs** – pairs existing heavy atoms with template coordinates for alignment anchors.
- **SVD** – computes rotation/translation via Kabsch/SVD, handling single/dual point fallbacks.
- **SynthesizeAtoms** – recreates missing heavy atoms (terminal OXT for peptides, OP3 for 5'-phosphorylated nucleic acids) using transformed template positions or tetrahedral geometry.
- **Cleanup** – removes atoms not present in the template to ensure canonical composition.

### 3.3 Hydro Pipeline (`ops::hydro`)

```mermaid
flowchart TD
    ScanResidues --> MarkDisulfide
    MarkDisulfide --> Protonation
    Protonation --> StripOldH
    StripOldH --> GeometryBuild
    GeometryBuild --> TerminalAdjust
    TerminalAdjust --> EndHydro
```

- **ScanResidues** – gathers standard residues and tracks their chain indices for ordered processing.
- **MarkDisulfide** – detects SG–SG pairs below the disulfide threshold and renames residues to `CYX`.
- **Protonation** – applies pH-driven heuristics plus HIS strategy selection to decide residue names.
- **StripOldH** – removes existing hydrogens if `remove_existing_h` is enabled.
- **GeometryBuild** – reconstructs hydrogens using template anchors, calling `reconstruct_geometry` and `calculate_transform` internally.
- **TerminalAdjust** – adds terminal hydrogens (N-term H1/H2/H3, C-term HOXT, nucleic 5' HO5' or pH-dependent phosphate HOP3, nucleic 3' HO3') respecting protonation states.
- **EndHydro** – structure now contains geometrically sound hydrogens.

### 3.4 Solvate Pipeline (`ops::solvate`)

```mermaid
flowchart TD
    Preclean --> Bounds
    Bounds --> Translate
    Translate --> GridPlace
    GridPlace --> RandomRotate
    RandomRotate --> ReplaceIons
    ReplaceIons --> Finalize
```

- **Preclean** – optionally removes existing solvent/ions before new placement.
- **Bounds** – computes bounding box and box vectors with configured margins.
- **Translate** – recenters solute so the solvent box origin is at margin offsets.
- **GridPlace** – populates a 3D lattice with HOH templates, skipping clashes via `SpatialGrid` checks.
- **RandomRotate** – orients each water copy with random Euler rotations for diversity.
- **ReplaceIons** – swaps selected waters with ions until the target charge is met, using RNG for distribution and reporting failure if insufficient.
- **Finalize** – appends the solvent chain and updates the structure’s periodic box.

### 3.5 Topology Builder (`ops::topology`)

```mermaid
flowchart TD
    PrepareOffsets --> IntraStandard
    IntraStandard --> HeteroTemplates
    HeteroTemplates --> TerminalBonds
    TerminalBonds --> InterResidue
    InterResidue --> DisulfideScan
    DisulfideScan --> EmitTopology
```

- **PrepareOffsets** – computes atom-index offsets for every residue to map local indexes to global topology indexes.
- **IntraStandard** – connects intra-residue bonds per template definitions, including hydrogen anchors.
- **HeteroTemplates** – injects bonds for hetero residues via user-provided `Template`s or errors if absent.
- **TerminalBonds** – adds special-case bonds for terminal atoms (H1/H2/H3, HOXT for peptides; P–OP3, OP3–HOP3, O5'–HO5', O3'–HO3' for nucleic acids).
- **InterResidue** – detects peptide and nucleic linkages by measuring atom distances against cutoffs.
- **DisulfideScan** – adds bonds between cystine sulfurs within the disulfide cutoff.
- **EmitTopology** – produces the final `Topology` pairing the structure with collected bonds.

### 3.6 Transform Utilities (`ops::transform`)

```mermaid
flowchart TD
    InputStructure --> ChooseOp
    ChooseOp --> Translate
    ChooseOp --> CenterGeom
    ChooseOp --> CenterMass
    ChooseOp --> RotateAxes
    RotateAxes --> UpdateBox
    Translate --> Output
    CenterGeom --> Output
    CenterMass --> Output
    UpdateBox --> Output
```

- **InputStructure** – mutable structure reference supplied by callers.
- **ChooseOp** – selects translation, geometric centering, mass centering, or Euler-based rotations.
- **Translate / CenterGeom / CenterMass** – adjust atomic coordinates directly, reusing vector math utilities.
- **RotateAxes** – applies axis or Euler rotations via `Rotation3`.
- **UpdateBox** – rotates the periodic box vectors to remain aligned with atom coordinates.
- **Output** – returns the transformed structure for subsequent ops or IO.

## 4. Algorithm Deep Dives

### 4.1 Template Alignment via SVD

- **Data selection** – `repair::calculate_transform` collects matched atom pairs between the residue and template; if only one or two matches exist, it falls back to translation-only or single-axis rotation.
- **Covariance build** – subtracts centroids and accumulates a 3×3 covariance matrix representing the correlation between residue and template frames.
- **SVD / Kabsch** – decomposes the covariance matrix with nalgebra’s SVD; determinant checks ensure right-handed rotations by flipping the final row when necessary.
- **Translation synthesis** – multiplies the template centroid by the rotation and subtracts from the residue centroid to obtain the translation vector.
- **Degradation path** – when the template lacks enough anchors, the algorithm returns an alignment error so the caller can report `AlignmentFailed` rather than emitting a distorted residue.

### 4.2 Hydrogen Reconstruction Geometry

- **Anchor selection** – each template hydrogen lists one or more anchor atoms; missing anchors trigger `IncompleteResidueForHydro` errors to avoid guesswork.
- **Rigid transform** – `reconstruct_geometry` retrieves the residue-specific transform (rotation + translation) derived from current heavy atoms and applies it to the template hydrogen coordinate.
- **Randomization** – none is applied for standard hydrogens, ensuring deterministic placement; terminals use evenly spaced tetrahedral vectors sorted by dot product to preserve orientation.
- **Terminal logic** – N-termini place up to three hydrogens arranged around the N–CA axis; C-termini add HOXT to OXT; nucleic 5'-terminals either add HO5' (no phosphate) or pH-dependent HOP3 (with phosphate, below pKₐ₂ ≈ 6.5); nucleic 3'-terminals always add HO3'.

### 4.3 Ion Replacement and Degradation Handling

- **Charge tracking** – `replace_with_ions` computes the total solute charge by summing template charges, accounting for keep ions, then compares against the requested total.
- **Candidate selection** – water residue IDs are shuffled and used as replacement slots; each swap decreases the charge delta by the ion’s charge.
- **Failure detection** – if the list depletes before the delta reaches zero, the function returns `IonizationFailed` (with details) or `BoxTooSmall` when no solvent exists, preventing silent mismatch.
- **Random fairness** – RNG with optional seed ensures reproducibility for tests or deterministic setups.

### 4.4 Topology Inference and Fallbacks

- **Template-driven bonds** – intra-residue connectivity is strictly template-backed; missing atoms raise `TopologyAtomMissing` unless the name is an optional terminal hydrogen.
- **Distance-based linkages** – peptide and nucleic linkages measure squared distances before adding bonds, so chains lacking adjacent residues simply skip the step.
- **Disulfide detection** – sulfurs closer than `disulfide_bond_cutoff` are connected; raising the cutoff enables permissive bonding while still preventing unrealistic links.
- **Hetero template injection** – user-supplied `Template`s cover ligands; omission surfaces `MissingHeteroTemplate` to encourage explicit topology definitions.

### 4.5 Error Propagation Strategy

- **thiserror-based enums** – `io::error::Error` and `ops::error::Error` capture precise failure contexts (source format, residue name, line numbers).
- **Eager validation** – template loading panics on duplicate names or malformed TOML via `include_str!`, guaranteeing consistent runtime data.
- **Graceful exits** – when hydrogens or alignments cannot be reconstructed, the caller receives descriptive errors rather than partial mutations, enabling upstream retry or user messaging.
