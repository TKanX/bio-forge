/**
 * @file Structure operations
 *
 * High-level functions for manipulating molecular structures.
 */

import type {
  WasmStructure,
  WasmTopology,
  WasmTemplate,
  StructureInfo,
  CleanConfig,
  HydroConfig,
  SolvateConfig,
  TopologyConfig,
} from "./types";

/** Supported structure formats */
export type StructureFormat = "pdb" | "mmcif";

/**
 * Get structure information.
 *
 * @param structure - WASM structure instance
 * @returns Structure statistics
 */
export function getStructureInfo(structure: WasmStructure): StructureInfo {
  return structure.info();
}

/**
 * Apply clean operation to structure.
 *
 * @param structure - WASM structure instance (mutated in place)
 * @param config - Clean configuration
 */
export function cleanStructure(
  structure: WasmStructure,
  config: CleanConfig
): void {
  structure.clean(config);
}

/**
 * Apply repair operation to structure.
 *
 * @param structure - WASM structure instance (mutated in place)
 */
export function repairStructure(structure: WasmStructure): void {
  structure.repair();
}

/**
 * Add hydrogens to structure.
 *
 * @param structure - WASM structure instance (mutated in place)
 * @param config - Hydrogen configuration
 */
export function addHydrogens(
  structure: WasmStructure,
  config: HydroConfig
): void {
  structure.addHydrogens(config);
}

/**
 * Solvate structure with water box and ions.
 *
 * @param structure - WASM structure instance (mutated in place)
 * @param config - Solvation configuration
 */
export function solvateStructure(
  structure: WasmStructure,
  config: SolvateConfig
): void {
  structure.solvate(config);
}

/**
 * Build topology from structure.
 *
 * @param structure - WASM structure instance
 * @param config - Topology configuration
 * @param templates - Optional templates for hetero residues
 * @returns Topology with bond information
 */
export function buildTopology(
  structure: WasmStructure,
  config?: TopologyConfig,
  templates?: WasmTemplate[]
): WasmTopology {
  return structure.toTopology(config, templates);
}

export type { WasmStructure, WasmTopology, WasmTemplate, StructureInfo };
