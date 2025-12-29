/**
 * @file File module exports
 */

export type { FileEntry, FileStatus, TemplateEntry } from "./types";

export {
  getFormatFromExtension,
  isStructureFile,
  isTemplateFile,
  readFileAsBytes,
  generateFileId,
  createFileEntry,
  createTemplateEntry,
  parseUploadedFiles,
  parseStructureBytes,
  parseTemplateBytes,
} from "./parser";

export {
  downloadFile,
  getExtension,
  generateOutputName,
  exportFileEntry,
  exportFilesAsZip,
  exportFiles,
} from "./exporter";
