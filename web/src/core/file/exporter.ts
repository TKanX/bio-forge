/**
 * @file File export utilities
 *
 * Functions for downloading and exporting processed files.
 */

import JSZip from "jszip";
import type { StructureFormat } from "../wasm";
import type { FileEntry } from "./types";

// ============================================================================
// Single File Export
// ============================================================================

/**
 * Download binary content as a file.
 *
 * @param bytes - File bytes
 * @param filename - Download filename
 * @param mimeType - MIME type (defaults to application/octet-stream)
 */
export function downloadFile(
  bytes: Uint8Array,
  filename: string,
  mimeType = "application/octet-stream"
): void {
  const buffer = new Uint8Array(bytes);
  const blob = new Blob([buffer], { type: mimeType });
  const url = URL.createObjectURL(blob);

  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  link.click();

  URL.revokeObjectURL(url);
}

/**
 * Get file extension for format.
 */
export function getExtension(format: StructureFormat): string {
  return format === "pdb" ? ".pdb" : ".cif";
}

/**
 * Generate output filename from input filename.
 *
 * @param inputName - Original filename
 * @param format - Output format
 * @param suffix - Optional suffix (defaults to "_forged")
 */
export function generateOutputName(
  inputName: string,
  format: StructureFormat,
  suffix = "_forged"
): string {
  const baseName = inputName.replace(/\.[^.]+$/, "");
  return `${baseName}${suffix}${getExtension(format)}`;
}

/**
 * Get serialized bytes for a file entry.
 *
 * Uses topology serialization (with bond connectivity) when available,
 * otherwise falls back to structure-only serialization.
 *
 * @param file - File entry to serialize
 * @param format - Output format
 * @returns Serialized bytes
 */
function serializeFile(file: FileEntry, format: StructureFormat): Uint8Array {
  if (file.topology) {
    return format === "pdb"
      ? file.topology.toPdbBytes()
      : file.topology.toMmcifBytes();
  }
  return format === "pdb"
    ? file.structure.toPdbBytes()
    : file.structure.toMmcifBytes();
}

/**
 * Export a single file entry.
 *
 * @param file - File entry to export
 * @param targetFormat - Output format
 */
export function exportFileEntry(
  file: FileEntry,
  targetFormat: StructureFormat
): void {
  const bytes = serializeFile(file, targetFormat);
  const filename = generateOutputName(file.name, targetFormat);
  downloadFile(bytes, filename);
}

// ============================================================================
// Batch Export
// ============================================================================

/**
 * Export multiple files as a ZIP archive.
 *
 * @param files - File entries to export
 * @param targetFormat - Output format for all files
 * @param zipName - Name of the ZIP file (defaults to "bioforge-output.zip")
 */
export async function exportFilesAsZip(
  files: FileEntry[],
  targetFormat: StructureFormat,
  zipName = "bioforge-output.zip"
): Promise<void> {
  const zip = new JSZip();

  for (const file of files) {
    const bytes = serializeFile(file, targetFormat);
    const filename = generateOutputName(file.name, targetFormat);
    zip.file(filename, bytes);
  }

  const blob = await zip.generateAsync({ type: "blob" });
  const url = URL.createObjectURL(blob);

  const link = document.createElement("a");
  link.href = url;
  link.download = zipName;
  link.click();

  URL.revokeObjectURL(url);
}

/**
 * Export files - single file direct download, multiple files as ZIP.
 *
 * @param files - File entries to export
 * @param targetFormat - Output format
 */
export async function exportFiles(
  files: FileEntry[],
  targetFormat: StructureFormat
): Promise<void> {
  if (files.length === 0) {
    return;
  }

  if (files.length === 1) {
    exportFileEntry(files[0], targetFormat);
  } else {
    await exportFilesAsZip(files, targetFormat);
  }
}
