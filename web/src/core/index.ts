/**
 * @file Core module exports
 */

export {
  initWasm,
  isWasmReady,
  getWasm,
  requireWasm,
  getStructureInfo,
  cleanStructure,
  repairStructure,
  addHydrogens,
  solvateStructure,
  buildTopology,
  exportTopology,
} from "./wasm";

export type {
  WasmModule,
  WasmStructure,
  WasmTopology,
  WasmTemplate,
  StructureFormat,
  StructureInfo,
  ChainInfo,
  PolymerType,
} from "./wasm";

export {
  DEFAULT_PIPELINE_CONFIG,
  DEFAULT_CLEAN_SETTINGS,
  DEFAULT_HYDRO_SETTINGS,
  DEFAULT_SOLVATE_SETTINGS,
  DEFAULT_TOPOLOGY_SETTINGS,
  executePipeline,
  executeBatch,
  yieldToEventLoop,
  PipelineError,
} from "./pipeline";

export type {
  PipelineConfig,
  CleanSettings,
  HydroSettings,
  SolvateSettings,
  TopologySettings,
  HisStrategy,
  CationSpecies,
  AnionSpecies,
  PipelineResult,
} from "./pipeline";

export {
  isStructureFile,
  isTemplateFile,
  readFileAsBytes,
  generateFileId,
  createFileEntry,
  createTemplateEntry,
  parseUploadedFiles,
  parseStructureBytes,
  parseTemplateBytes,
  downloadFile,
  getExtension,
  generateOutputName,
  exportFileEntry,
  exportFilesAsZip,
  exportFiles,
} from "./file";

export type { FileEntry, FileStatus, TemplateEntry } from "./file";
