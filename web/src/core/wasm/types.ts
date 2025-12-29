/**
 * @file WASM type definitions
 *
 * TypeScript interfaces for the BioForge WASM bindings.
 */

// ============================================================================
// Structure Information Types
// ============================================================================

/** Chain polymer classification */
export type PolymerType = "protein" | "nucleic" | "solvent" | "hetero";

/** Per-chain statistics */
export interface ChainInfo {
  readonly id: string;
  readonly residueCount: number;
  readonly atomCount: number;
  readonly polymerTypes: PolymerType[];
}

/** Complete structure statistics */
export interface StructureInfo {
  readonly chainCount: number;
  readonly residueCount: number;
  readonly atomCount: number;
  readonly boxLengths?: [number, number, number];
  readonly boxAngles?: [number, number, number];
  readonly chains: ChainInfo[];
}

// ============================================================================
// Configuration Types
// ============================================================================

/** Clean operation configuration */
export interface CleanConfig {
  removeWater?: boolean;
  removeIons?: boolean;
  removeHydrogens?: boolean;
  removeHetero?: boolean;
  removeResidueNames?: string[];
  keepResidueNames?: string[];
}

/** Histidine tautomer selection strategy */
export type HisStrategy = "hid" | "hie" | "random" | "network";

/** Hydrogen addition configuration */
export interface HydroConfig {
  targetPh?: number;
  removeExistingH?: boolean;
  hisStrategy?: HisStrategy;
}

/** Cation species for solvation */
export type CationSpecies = "Na" | "K" | "Mg" | "Ca" | "Li" | "Zn";

/** Anion species for solvation */
export type AnionSpecies = "Cl" | "Br" | "I" | "F";

/** Solvation configuration */
export interface SolvateConfig {
  margin?: number;
  waterSpacing?: number;
  vdwCutoff?: number;
  removeExisting?: boolean;
  cations?: CationSpecies[];
  anions?: AnionSpecies[];
  targetCharge?: number;
  rngSeed?: number;
}

/** Topology building configuration */
export interface TopologyConfig {
  disulfideCutoff?: number;
}

// ============================================================================
// WASM Class Interfaces
// ============================================================================

/** WASM Structure class interface */
export interface WasmStructure {
  free(): void;
  clone(): WasmStructure;
  info(): StructureInfo;
  toPdb(): string;
  toPdbBytes(): Uint8Array;
  toMmcif(): string;
  toMmcifBytes(): Uint8Array;
  clean(config?: CleanConfig): void;
  repair(): void;
  addHydrogens(config?: HydroConfig): void;
  solvate(config?: SolvateConfig): void;
  toTopology(config?: TopologyConfig, templates?: WasmTemplate[]): WasmTopology;
  readonly chainCount: number;
  readonly residueCount: number;
  readonly atomCount: number;
}

/** WASM Topology class interface */
export interface WasmTopology {
  free(): void;
  readonly structure: WasmStructure;
  readonly bondCount: number;
  toPdb(): string;
  toPdbBytes(): Uint8Array;
  toMmcif(): string;
  toMmcifBytes(): Uint8Array;
}

/** WASM Template class interface */
export interface WasmTemplate {
  free(): void;
  readonly name: string;
}

// ============================================================================
// WASM Module Interface
// ============================================================================

/** The complete WASM module interface */
export interface WasmModule {
  Structure: {
    fromPdb(content: string): WasmStructure;
    fromPdbBytes(content: Uint8Array): WasmStructure;
    fromMmcif(content: string): WasmStructure;
    fromMmcifBytes(content: Uint8Array): WasmStructure;
  };
  Topology: {
    fromStructure(
      structure: WasmStructure,
      config?: TopologyConfig,
      templates?: WasmTemplate[]
    ): WasmTopology;
  };
  Template: {
    fromMol2(content: string): WasmTemplate;
    fromMol2Bytes(content: Uint8Array): WasmTemplate;
  };
}
