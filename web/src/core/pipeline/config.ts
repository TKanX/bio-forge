/**
 * @file Pipeline configuration types
 *
 * Type definitions and default values for all pipeline step configurations.
 */

import type {
  CleanConfig,
  HydroConfig,
  HisStrategy,
  SolvateConfig,
  TopologyConfig,
  CationSpecies,
  AnionSpecies,
} from "../wasm";

// ============================================================================
// Clean Configuration
// ============================================================================

/** Full clean configuration with all fields required */
export interface CleanSettings {
  removeWater: boolean;
  removeIons: boolean;
  removeHydrogens: boolean;
  removeHetero: boolean;
  removeResidueNames: string[];
  keepResidueNames: string[];
}

/** Default clean settings */
export const DEFAULT_CLEAN_SETTINGS: CleanSettings = {
  removeWater: true,
  removeIons: false,
  removeHydrogens: false,
  removeHetero: false,
  removeResidueNames: [],
  keepResidueNames: [],
};

/** Convert settings to WASM config */
export function toCleanConfig(settings: CleanSettings): CleanConfig {
  return { ...settings };
}

// ============================================================================
// Hydro Configuration
// ============================================================================

/** Full hydrogen configuration with all fields required */
export interface HydroSettings {
  targetPh: number | undefined;
  removeExistingH: boolean;
  hisStrategy: HisStrategy;
}

/** Default hydrogen settings */
export const DEFAULT_HYDRO_SETTINGS: HydroSettings = {
  targetPh: undefined,
  removeExistingH: true,
  hisStrategy: "network",
};

/** Convert settings to WASM config */
export function toHydroConfig(settings: HydroSettings): HydroConfig {
  return {
    targetPh: settings.targetPh,
    removeExistingH: settings.removeExistingH,
    hisStrategy: settings.hisStrategy,
  };
}

// ============================================================================
// Solvate Configuration
// ============================================================================

/** Full solvation configuration with all fields required */
export interface SolvateSettings {
  margin: number;
  waterSpacing: number;
  vdwCutoff: number;
  removeExisting: boolean;
  cations: CationSpecies[];
  anions: AnionSpecies[];
  targetCharge: number;
  rngSeed: number | undefined;
}

/** Default solvation settings */
export const DEFAULT_SOLVATE_SETTINGS: SolvateSettings = {
  margin: 10.0,
  waterSpacing: 3.1,
  vdwCutoff: 2.4,
  removeExisting: true,
  cations: ["Na"],
  anions: ["Cl"],
  targetCharge: 0,
  rngSeed: undefined,
};

/** Convert settings to WASM config */
export function toSolvateConfig(settings: SolvateSettings): SolvateConfig {
  return {
    margin: settings.margin,
    waterSpacing: settings.waterSpacing,
    vdwCutoff: settings.vdwCutoff,
    removeExisting: settings.removeExisting,
    cations: settings.cations,
    anions: settings.anions,
    targetCharge: settings.targetCharge,
    rngSeed: settings.rngSeed,
  };
}

// ============================================================================
// Topology Configuration
// ============================================================================

/** Full topology configuration with all fields required */
export interface TopologySettings {
  disulfideCutoff: number;
}

/** Default topology settings */
export const DEFAULT_TOPOLOGY_SETTINGS: TopologySettings = {
  disulfideCutoff: 2.2,
};

/** Convert settings to WASM config */
export function toTopologyConfig(settings: TopologySettings): TopologyConfig {
  return { disulfideCutoff: settings.disulfideCutoff };
}

// ============================================================================
// Pipeline Configuration
// ============================================================================

/** Complete pipeline configuration */
export interface PipelineConfig {
  clean: { enabled: boolean; settings: CleanSettings };
  repair: { enabled: boolean };
  hydro: { enabled: boolean; settings: HydroSettings };
  solvate: { enabled: boolean; settings: SolvateSettings };
  topology: { enabled: boolean; settings: TopologySettings };
}

/** Default pipeline configuration */
export const DEFAULT_PIPELINE_CONFIG: PipelineConfig = {
  clean: { enabled: true, settings: DEFAULT_CLEAN_SETTINGS },
  repair: { enabled: true },
  hydro: { enabled: false, settings: DEFAULT_HYDRO_SETTINGS },
  solvate: { enabled: false, settings: DEFAULT_SOLVATE_SETTINGS },
  topology: { enabled: false, settings: DEFAULT_TOPOLOGY_SETTINGS },
};

export type { CleanConfig, HydroConfig, SolvateConfig, TopologyConfig };
export type { HisStrategy, CationSpecies, AnionSpecies };
