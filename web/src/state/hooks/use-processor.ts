/**
 * @file useProcessor hook
 *
 * Hook for managing the pipeline processing workflow.
 */

import { useCallback, useMemo } from "react";
import { useShallow } from "zustand/react/shallow";
import {
  useFileStore,
  usePipelineStore,
  useUIStore,
  showSuccess,
  showError,
  showWarning,
} from "../stores";
import {
  executePipeline,
  yieldToEventLoop,
  PipelineError,
  type FileEntry,
  type PipelineConfig,
} from "@/core";

/**
 * Hook providing pipeline processing functionality.
 */
export function useProcessor() {
  /** Batch file store selectors with shallow comparison */
  const { files, selectedIds, templates, updateFileStatus, updateFileInfo } =
    useFileStore(
      useShallow((s) => ({
        files: s.files,
        selectedIds: s.selectedIds,
        templates: s.templates,
        updateFileStatus: s.updateFileStatus,
        updateFileInfo: s.updateFileInfo,
      }))
    );

  /** Batch pipeline store selectors with shallow comparison */
  const {
    cleanEnabled,
    repairEnabled,
    hydroEnabled,
    solvateEnabled,
    topologyEnabled,
    cleanConfig,
    hydroConfig,
    solvateConfig,
    topologyConfig,
  } = usePipelineStore(
    useShallow((s) => ({
      cleanEnabled: s.cleanEnabled,
      repairEnabled: s.repairEnabled,
      hydroEnabled: s.hydroEnabled,
      solvateEnabled: s.solvateEnabled,
      topologyEnabled: s.topologyEnabled,
      cleanConfig: s.cleanConfig,
      hydroConfig: s.hydroConfig,
      solvateConfig: s.solvateConfig,
      topologyConfig: s.topologyConfig,
    }))
  );

  /** Batch UI store selectors */
  const { isProcessing, setProcessing } = useUIStore(
    useShallow((s) => ({
      isProcessing: s.isProcessing,
      setProcessing: s.setProcessing,
    }))
  );

  // ============================================================================
  // Derived State - Files by status
  // ============================================================================

  /** Files that can be processed (ready or completed) */
  const processableFiles = useMemo(
    () => files.filter((f) => f.status === "ready" || f.status === "completed"),
    [files]
  );

  /** Files that have been processed successfully */
  const completedFiles = useMemo(
    () => files.filter((f) => f.status === "completed"),
    [files]
  );

  /** Selected files that can be processed */
  const selectedProcessable = useMemo(
    () => processableFiles.filter((f) => selectedIds.has(f.id)),
    [processableFiles, selectedIds]
  );

  /** Selected files that are completed */
  const selectedCompleted = useMemo(
    () => completedFiles.filter((f) => selectedIds.has(f.id)),
    [completedFiles, selectedIds]
  );

  /** Files to process (selected if any, otherwise all processable) */
  const filesToProcess = useMemo(
    () =>
      selectedProcessable.length > 0 ? selectedProcessable : processableFiles,
    [selectedProcessable, processableFiles]
  );

  /** Files to download (selected completed if any, otherwise all completed) */
  const filesToDownload = useMemo(
    () => (selectedCompleted.length > 0 ? selectedCompleted : completedFiles),
    [selectedCompleted, completedFiles]
  );

  /** Build pipeline config from flat store properties */
  const config = useMemo<PipelineConfig>(
    () => ({
      clean: { enabled: cleanEnabled, settings: cleanConfig },
      repair: { enabled: repairEnabled },
      hydro: { enabled: hydroEnabled, settings: hydroConfig },
      solvate: { enabled: solvateEnabled, settings: solvateConfig },
      topology: { enabled: topologyEnabled, settings: topologyConfig },
    }),
    [
      cleanEnabled,
      repairEnabled,
      hydroEnabled,
      solvateEnabled,
      topologyEnabled,
      cleanConfig,
      hydroConfig,
      solvateConfig,
      topologyConfig,
    ]
  );

  /** Whether at least one pipeline step is enabled */
  const hasAnyStepEnabled =
    cleanEnabled ||
    repairEnabled ||
    hydroEnabled ||
    solvateEnabled ||
    topologyEnabled;

  /** Whether processing is currently possible */
  const canProcess =
    hasAnyStepEnabled && filesToProcess.length > 0 && !isProcessing;

  /** Whether download is currently possible */
  const canDownload = filesToDownload.length > 0;

  /**
   * Internal: Execute pipeline on given files.
   * Mutates structures in-place and updates file info.
   */
  const executeOnFiles = useCallback(
    async (targetFiles: FileEntry[]) => {
      if (targetFiles.length === 0) {
        showWarning("No files ready to process");
        return;
      }

      setProcessing(true);
      let successCount = 0;

      for (const file of targetFiles) {
        updateFileStatus(file.id, "processing");
        await yieldToEventLoop();

        try {
          const result = executePipeline(
            file,
            config,
            templates.length > 0 ? templates : undefined
          );

          updateFileInfo(file.id, result.info);
          updateFileStatus(file.id, "completed");
          successCount++;
        } catch (error) {
          const message =
            error instanceof PipelineError
              ? error.message
              : error instanceof Error
                ? error.message
                : "Processing failed";
          updateFileStatus(file.id, "error", message);
          showError(`Failed to process ${file.name}: ${message}`);
        }
      }

      setProcessing(false);

      if (successCount > 0) {
        showSuccess(
          `Processed ${successCount} file${successCount !== 1 ? "s" : ""}`
        );
      }
    },
    [templates, config, updateFileStatus, updateFileInfo, setProcessing]
  );

  /**
   * Process files (selected if any, otherwise all processable).
   */
  const processFiles = useCallback(async () => {
    await executeOnFiles(filesToProcess);
  }, [executeOnFiles, filesToProcess]);

  /**
   * Process only selected files.
   */
  const processSelected = useCallback(async () => {
    if (selectedProcessable.length === 0) {
      showWarning("No selected files ready to process");
      return;
    }
    await executeOnFiles(selectedProcessable);
  }, [executeOnFiles, selectedProcessable]);

  /**
   * Process all processable files (regardless of selection).
   */
  const processAll = useCallback(async () => {
    await executeOnFiles(processableFiles);
  }, [executeOnFiles, processableFiles]);

  /**
   * Download files (selected if any, otherwise all completed).
   */
  const downloadFiles = useCallback(
    async (format: "pdb" | "mmcif" = "pdb", asZip = false) => {
      if (filesToDownload.length === 0) {
        showWarning("No completed files to download");
        return;
      }

      const { exportFiles, exportFilesAsZip } = await import("@/core");

      if (asZip || filesToDownload.length > 1) {
        await exportFilesAsZip(filesToDownload, format);
      } else {
        await exportFiles(filesToDownload, format);
      }
    },
    [filesToDownload]
  );

  /**
   * Download only selected completed files.
   */
  const downloadSelected = useCallback(
    async (format: "pdb" | "mmcif" = "pdb", asZip = false) => {
      if (selectedCompleted.length === 0) {
        showWarning("No selected completed files to download");
        return;
      }

      const { exportFiles, exportFilesAsZip } = await import("@/core");

      if (asZip || selectedCompleted.length > 1) {
        await exportFilesAsZip(selectedCompleted, format);
      } else {
        await exportFiles(selectedCompleted, format);
      }
    },
    [selectedCompleted]
  );

  /**
   * Download all completed files (regardless of selection).
   */
  const downloadAll = useCallback(
    async (format: "pdb" | "mmcif" = "pdb", asZip = false) => {
      if (completedFiles.length === 0) {
        showWarning("No completed files to download");
        return;
      }

      const { exportFiles, exportFilesAsZip } = await import("@/core");

      if (asZip || completedFiles.length > 1) {
        await exportFilesAsZip(completedFiles, format);
      } else {
        await exportFiles(completedFiles, format);
      }
    },
    [completedFiles]
  );

  // ============================================================================
  // Return
  // ============================================================================

  return {
    // State
    isProcessing,
    hasAnyStepEnabled,
    canProcess,
    canDownload,

    // Derived file lists
    processableFiles,
    completedFiles,
    filesToProcess,
    filesToDownload,

    // Actions
    processFiles,
    processSelected,
    processAll,
    downloadFiles,
    downloadSelected,
    downloadAll,
  };
}
