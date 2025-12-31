/**
 * @file File store
 *
 * State management for uploaded files and templates.
 */

import { create } from "zustand";
import type {
  FileEntry,
  FileStatus,
  TemplateEntry,
  StructureInfo,
  WasmTopology,
} from "@/core";

// ============================================================================
// Types
// ============================================================================

interface FileState {
  /** List of uploaded structure files */
  files: FileEntry[];
  /** Set of selected file IDs */
  selectedIds: Set<string>;
  /** List of uploaded templates */
  templates: TemplateEntry[];
}

interface FileActions {
  // File management
  addFiles: (files: FileEntry[]) => void;
  removeFile: (id: string) => void;
  removeFiles: (ids: string[]) => void;
  clearFiles: () => void;

  // File updates
  updateFileStatus: (id: string, status: FileStatus, error?: string) => void;
  updateFileInfo: (id: string, info: StructureInfo) => void;
  updateFileResult: (
    id: string,
    info: StructureInfo,
    topology?: WasmTopology
  ) => void;
  clearFileTopology: (id: string) => void;

  // Selection
  selectFile: (id: string) => void;
  deselectFile: (id: string) => void;
  toggleFileSelection: (id: string) => void;
  selectAllFiles: () => void;
  clearSelection: () => void;
  removeSelectedFiles: () => void;

  // Bulk status operations
  removeCompletedFiles: () => void;

  // Templates
  addTemplates: (templates: TemplateEntry[]) => void;
  removeTemplate: (id: string) => void;
  clearTemplates: () => void;
}

type FileStore = FileState & FileActions;

// ============================================================================
// Helpers
// ============================================================================

/**
 * Safely free a file's WASM resources (structure and topology).
 */
function freeFileResources(file: FileEntry): void {
  file.structure.free();
  file.topology?.free();
}

// ============================================================================
// Store
// ============================================================================

export const useFileStore = create<FileStore>((set, get) => ({
  files: [],
  selectedIds: new Set(),
  templates: [],

  addFiles: (files) =>
    set((state) => ({
      files: [...state.files, ...files],
    })),

  removeFile: (id) => {
    // Query: capture resource reference before state update
    const file = get().files.find((f) => f.id === id);

    // Mutate: pure immutable state update
    set((state) => {
      const newSelectedIds = new Set(state.selectedIds);
      newSelectedIds.delete(id);
      return {
        files: state.files.filter((f) => f.id !== id),
        selectedIds: newSelectedIds,
      };
    });

    // Cleanup: free WASM resources after state update
    if (file) {
      freeFileResources(file);
    }
  },

  removeFiles: (ids) => {
    // Query: capture resource references before state update
    const idSet = new Set(ids);
    const filesToRemove = get().files.filter((f) => idSet.has(f.id));

    // Mutate: pure immutable state update
    set((state) => {
      const newSelectedIds = new Set(state.selectedIds);
      ids.forEach((id) => newSelectedIds.delete(id));
      return {
        files: state.files.filter((f) => !idSet.has(f.id)),
        selectedIds: newSelectedIds,
      };
    });

    // Cleanup: free WASM resources after state update
    filesToRemove.forEach((f) => freeFileResources(f));
  },

  clearFiles: () => {
    // Query: capture resource references before state update
    const filesToClear = get().files;

    // Mutate: pure immutable state update
    set({ files: [], selectedIds: new Set() });

    // Cleanup: free WASM resources after state update
    filesToClear.forEach((f) => freeFileResources(f));
  },

  updateFileStatus: (id, status, error) =>
    set((state) => ({
      files: state.files.map((f) =>
        f.id === id
          ? {
              ...f,
              status,
              error,
              // Increment version when completed to trigger re-render
              version: status === "completed" ? f.version + 1 : f.version,
            }
          : f
      ),
    })),

  updateFileInfo: (id, info) =>
    set((state) => ({
      files: state.files.map((f) =>
        f.id === id ? { ...f, info, status: "ready" as const } : f
      ),
    })),

  updateFileResult: (id, info, topology) => {
    // Query: capture old topology reference before state update
    const oldTopology = get().files.find((f) => f.id === id)?.topology;

    // Mutate: pure immutable state update
    set((state) => ({
      files: state.files.map((f) =>
        f.id !== id
          ? f
          : {
              ...f,
              info,
              topology,
              status: "completed" as const,
              version: f.version + 1,
            }
      ),
    }));

    // Cleanup: free old topology if it has been replaced or cleared
    if (oldTopology && oldTopology !== topology) {
      oldTopology.free();
    }
  },

  clearFileTopology: (id) => {
    // Query: capture topology reference before state update
    const oldTopology = get().files.find((f) => f.id === id)?.topology;

    // Mutate: pure immutable state update
    set((state) => ({
      files: state.files.map((f) =>
        f.id !== id ? f : { ...f, topology: undefined }
      ),
    }));

    // Cleanup: free topology after state update
    oldTopology?.free();
  },

  selectFile: (id) =>
    set((state) => {
      const newSelectedIds = new Set(state.selectedIds);
      newSelectedIds.add(id);
      return { selectedIds: newSelectedIds };
    }),

  deselectFile: (id) =>
    set((state) => {
      const newSelectedIds = new Set(state.selectedIds);
      newSelectedIds.delete(id);
      return { selectedIds: newSelectedIds };
    }),

  toggleFileSelection: (id) => {
    const { selectedIds } = get();
    if (selectedIds.has(id)) {
      get().deselectFile(id);
    } else {
      get().selectFile(id);
    }
  },

  selectAllFiles: () =>
    set((state) => ({
      selectedIds: new Set(state.files.map((f) => f.id)),
    })),

  clearSelection: () => set({ selectedIds: new Set() }),

  removeSelectedFiles: () => {
    const { selectedIds } = get();
    get().removeFiles([...selectedIds]);
  },

  removeCompletedFiles: () => {
    // Query: capture completed files before state update
    const completedFiles = get().files.filter((f) => f.status === "completed");
    const completedIds = new Set(completedFiles.map((f) => f.id));

    // Mutate: pure immutable state update
    set((state) => ({
      files: state.files.filter((f) => f.status !== "completed"),
      selectedIds: new Set(
        [...state.selectedIds].filter((id) => !completedIds.has(id))
      ),
    }));

    // Cleanup: free WASM resources after state update
    completedFiles.forEach((f) => freeFileResources(f));
  },

  addTemplates: (templates) =>
    set((state) => ({
      templates: [...state.templates, ...templates],
    })),

  removeTemplate: (id) => {
    // Query: capture template reference before state update
    const template = get().templates.find((t) => t.id === id);

    // Mutate: pure immutable state update
    set((state) => ({
      templates: state.templates.filter((t) => t.id !== id),
    }));

    // Cleanup: free WASM resource after state update
    template?.template.free();
  },

  clearTemplates: () => {
    // Query: capture template references before state update
    const templatesToClear = get().templates;

    // Mutate: pure immutable state update
    set({ templates: [] });

    // Cleanup: free WASM resources after state update
    templatesToClear.forEach((t) => t.template.free());
  },
}));

// ============================================================================
// Selectors
// ============================================================================

/** Get files by status */
export const selectFilesByStatus = (status: FileStatus) => (state: FileStore) =>
  state.files.filter((f) => f.status === status);

/** Get selected files */
export const selectSelectedFiles = (state: FileStore) =>
  state.files.filter((f) => state.selectedIds.has(f.id));

/** Get files ready for processing (ready or completed) */
export const selectProcessableFiles = (state: FileStore) =>
  state.files.filter((f) => f.status === "ready" || f.status === "completed");

/** Check if all files are selected */
export const selectAllSelected = (state: FileStore) =>
  state.files.length > 0 && state.selectedIds.size === state.files.length;

/** Get completed file count */
export const selectCompletedCount = (state: FileStore) =>
  state.files.filter((f) => f.status === "completed").length;
