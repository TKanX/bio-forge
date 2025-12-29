/**
 * @file Pipeline executor
 *
 * Core logic for executing the structure processing pipeline.
 */

import {
  cleanStructure,
  repairStructure,
  addHydrogens,
  solvateStructure,
  buildTopology,
  getStructureInfo,
  type WasmTopology,
  type StructureInfo,
} from "../wasm";
import type { FileEntry, TemplateEntry } from "../file";
import {
  type PipelineConfig,
  toCleanConfig,
  toHydroConfig,
  toSolvateConfig,
  toTopologyConfig,
} from "./config";

// ============================================================================
// Types
// ============================================================================

/** Result from successful pipeline execution */
export interface PipelineResult {
  /** Updated structure info */
  info: StructureInfo;
  /** Bond count (if topology was built) */
  bondCount?: number;
}

/** Pipeline execution error */
export class PipelineError extends Error {
  constructor(
    message: string,
    public readonly step?: string
  ) {
    super(message);
    this.name = "PipelineError";
  }
}

// ============================================================================
// Executor
// ============================================================================

/**
 * Execute the processing pipeline on a file's structure.
 *
 * Mutates the structure in-place. The file's structure object is modified directly.
 *
 * @param file - File entry with structure to process
 * @param config - Pipeline configuration
 * @param templates - Optional template entries for topology building
 * @returns Pipeline result with updated info
 * @throws PipelineError if any step fails
 */
export function executePipeline(
  file: FileEntry,
  config: PipelineConfig,
  templates?: TemplateEntry[]
): PipelineResult {
  const structure = file.structure;
  let topology: WasmTopology | null = null;

  try {
    // Step 1: Clean
    if (config.clean.enabled) {
      cleanStructure(structure, toCleanConfig(config.clean.settings));
    }

    // Step 2: Repair
    if (config.repair.enabled) {
      repairStructure(structure);
    }

    // Step 3: Hydrogenate
    if (config.hydro.enabled) {
      addHydrogens(structure, toHydroConfig(config.hydro.settings));
    }

    // Step 4: Solvate
    if (config.solvate.enabled) {
      solvateStructure(structure, toSolvateConfig(config.solvate.settings));
    }

    // Step 5: Build topology
    let bondCount: number | undefined;
    if (config.topology.enabled) {
      const wasmTemplates = templates?.map((t) => t.template);

      topology = buildTopology(
        structure,
        toTopologyConfig(config.topology.settings),
        wasmTemplates
      );
      bondCount = topology.bondCount;
    }

    const info = getStructureInfo(structure);

    return { info, bondCount };
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new PipelineError(message);
  } finally {
    // Clean up topology (structure is NOT freed - it's still used)
    // Note: Templates are NOT freed here - they're owned by the store
    topology?.free();
  }
}

/**
 * Yield control to the event loop.
 * Use between processing files to keep UI responsive.
 */
export function yieldToEventLoop(): Promise<void> {
  return new Promise((resolve) => setTimeout(resolve, 0));
}

/**
 * Execute pipeline on multiple files with progress callback.
 *
 * @param files - Array of file entries to process
 * @param config - Pipeline configuration
 * @param templates - Optional template entries for topology building
 * @param onProgress - Callback for progress updates
 * @returns Map of results/errors by file ID
 */
export async function executeBatch(
  files: FileEntry[],
  config: PipelineConfig,
  templates?: TemplateEntry[],
  onProgress?: (
    id: string,
    status: "processing" | "completed" | "error"
  ) => void
): Promise<Map<string, PipelineResult | PipelineError>> {
  const results = new Map<string, PipelineResult | PipelineError>();

  for (const file of files) {
    onProgress?.(file.id, "processing");
    await yieldToEventLoop();

    try {
      const result = executePipeline(file, config, templates);
      results.set(file.id, result);
      onProgress?.(file.id, "completed");
    } catch (error) {
      const pipelineError =
        error instanceof PipelineError
          ? error
          : new PipelineError(
              error instanceof Error ? error.message : "Unknown error"
            );
      results.set(file.id, pipelineError);
      onProgress?.(file.id, "error");
    }
  }

  return results;
}
