/**
 * @file File type definitions
 *
 * Types for file management in the application.
 */

import type {
  StructureFormat,
  StructureInfo,
  WasmStructure,
  WasmTemplate,
} from "../wasm";

// ============================================================================
// File Status
// ============================================================================

/** Processing status for a file */
export type FileStatus =
  | "pending" // Just uploaded, not yet analyzed
  | "ready" // Analyzed and ready for processing
  | "processing" // Currently being processed
  | "completed" // Successfully processed
  | "error"; // Processing failed

// ============================================================================
// File Entry
// ============================================================================

/**
 * Managed file entry.
 */
export interface FileEntry {
  /** Unique identifier */
  readonly id: string;
  /** Original filename */
  readonly name: string;
  /** Original file size in bytes */
  readonly size: number;
  /** Structure format (pdb or mmcif) */
  readonly format: StructureFormat;
  /** Original file bytes (for initial Mol* render) */
  readonly rawBytes: Uint8Array;
  /** Live WASM structure object (for all operations) */
  structure: WasmStructure;
  /** Current processing status */
  status: FileStatus;
  /** Error message if status is 'error' */
  error?: string;
  /** Cached structure information (updated after operations) */
  info?: StructureInfo;
  /** Version counter (incremented on structure mutation) */
  version: number;
}

// ============================================================================
// Template Entry
// ============================================================================

/**
 * MOL2 template for topology building.
 */
export interface TemplateEntry {
  /** Unique identifier */
  readonly id: string;
  /** Residue name from template */
  readonly name: string;
  /** Original file size in bytes */
  readonly size: number;
  /** Original file bytes */
  readonly rawBytes: Uint8Array;
  /** Live WASM template object */
  template: WasmTemplate;
}
