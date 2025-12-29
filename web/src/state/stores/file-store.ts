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

  removeFile: (id) =>
    set((state) => {
      // Free WASM resource before removing
      const file = state.files.find((f) => f.id === id);
      file?.structure.free();

      const newSelectedIds = new Set(state.selectedIds);
      newSelectedIds.delete(id);
      return {
        files: state.files.filter((f) => f.id !== id),
        selectedIds: newSelectedIds,
      };
    }),

  removeFiles: (ids) =>
    set((state) => {
      const idSet = new Set(ids);

      // Free WASM resources before removing
      state.files
        .filter((f) => idSet.has(f.id))
        .forEach((f) => f.structure.free());

      const newSelectedIds = new Set(state.selectedIds);
      ids.forEach((id) => newSelectedIds.delete(id));
      return {
        files: state.files.filter((f) => !idSet.has(f.id)),
        selectedIds: newSelectedIds,
      };
    }),

  clearFiles: () =>
    set((state) => {
      // Free all WASM resources
      state.files.forEach((f) => f.structure.free());
      return {
        files: [],
        selectedIds: new Set(),
      };
    }),

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

  removeCompletedFiles: () =>
    set((state) => {
      const completedFiles = state.files.filter(
        (f) => f.status === "completed"
      );
      // Free WASM resources
      completedFiles.forEach((f) => f.structure.free());

      const completedIds = new Set(completedFiles.map((f) => f.id));
      const newSelectedIds = new Set(
        [...state.selectedIds].filter((id) => !completedIds.has(id))
      );
      return {
        files: state.files.filter((f) => f.status !== "completed"),
        selectedIds: newSelectedIds,
      };
    }),

  addTemplates: (templates) =>
    set((state) => ({
      templates: [...state.templates, ...templates],
    })),

  removeTemplate: (id) =>
    set((state) => {
      // Free WASM resource before removing
      const template = state.templates.find((t) => t.id === id);
      template?.template.free();

      return {
        templates: state.templates.filter((t) => t.id !== id),
      };
    }),

  clearTemplates: () =>
    set((state) => {
      // Free all WASM resources
      state.templates.forEach((t) => t.template.free());
      return { templates: [] };
    }),
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
