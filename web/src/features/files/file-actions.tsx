/**
 * @file File actions
 *
 * Action bar for file operations (upload, process, download, clear).
 */

"use client";

import { useRef } from "react";
import { useShallow } from "zustand/react/shallow";
import { Button, IconButton, Tooltip } from "@/ui/primitives";
import {
  UploadIcon,
  PlayIcon,
  DownloadIcon,
  TrashIcon,
  RefreshIcon,
} from "@/ui/icons";
import { useFileStore, useFileUpload, showError } from "@/state";
import { useProcessor } from "@/state/hooks";
import { ACCEPTED_STRUCTURE_EXTENSIONS } from "@/lib/constants";

// ============================================================================
// Component
// ============================================================================

export function FileActions() {
  const { files, selectedIds, clearFiles, selectAllFiles, clearSelection } =
    useFileStore(
      useShallow((s) => ({
        files: s.files,
        selectedIds: s.selectedIds,
        clearFiles: s.clearFiles,
        selectAllFiles: s.selectAllFiles,
        clearSelection: s.clearSelection,
      }))
    );

  const { handleFiles } = useFileUpload();

  const {
    isProcessing,
    canProcess,
    canDownload,
    hasAnyStepEnabled,
    processFiles,
    downloadFiles,
  } = useProcessor();

  const fileInputRef = useRef<HTMLInputElement>(null);

  const hasFiles = files.length > 0;
  const hasSelected = selectedIds.size > 0;
  const allSelected = hasFiles && selectedIds.size === files.length;

  const handleInputChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const fileList = Array.from(e.target.files || []);
    if (fileList.length > 0) {
      handleFiles(fileList);
    }
    e.target.value = "";
  };

  const handleProcess = async () => {
    if (!canProcess) return;

    try {
      await processFiles();
    } catch (error) {
      showError(error instanceof Error ? error.message : "Processing failed");
    }
  };

  const handleDownload = async () => {
    if (!canDownload) return;

    try {
      await downloadFiles();
    } catch (error) {
      showError(error instanceof Error ? error.message : "Download failed");
    }
  };

  const handleToggleSelectAll = () => {
    if (allSelected) {
      clearSelection();
    } else {
      selectAllFiles();
    }
  };

  return (
    <div className="p-3 sm:p-4 border-b border-border flex items-center gap-2 sm:gap-3 flex-wrap">
      {/* Upload */}
      <Button
        variant="primary"
        size="sm"
        onClick={() => fileInputRef.current?.click()}
        className="px-2 sm:px-3"
      >
        <UploadIcon className="size-4" />
        <span className="hidden xs:inline">Upload</span>
      </Button>
      <input
        ref={fileInputRef}
        type="file"
        accept={ACCEPTED_STRUCTURE_EXTENSIONS.join(",")}
        multiple
        onChange={handleInputChange}
        className="hidden"
      />

      <div className="h-6 w-px bg-border hidden xs:block" />

      {/* Process */}
      <Tooltip
        content={
          !hasAnyStepEnabled ? "Enable at least one pipeline step" : undefined
        }
      >
        <Button
          variant="secondary"
          size="sm"
          onClick={handleProcess}
          disabled={!canProcess || isProcessing}
          className="px-2 sm:px-3"
        >
          {isProcessing ? (
            <RefreshIcon className="size-4 animate-spin" />
          ) : (
            <PlayIcon className="size-4" />
          )}
          <span className="hidden sm:inline">
            {hasSelected ? `Process (${selectedIds.size})` : "Process All"}
          </span>
        </Button>
      </Tooltip>

      {/* Download */}
      <Button
        variant="secondary"
        size="sm"
        onClick={handleDownload}
        disabled={!canDownload}
        className="px-2 sm:px-3"
      >
        <DownloadIcon className="size-4" />
        <span className="hidden sm:inline">
          {hasSelected ? `Download (${selectedIds.size})` : "Download All"}
        </span>
      </Button>

      {/* Spacer */}
      <div className="flex-1" />

      {/* Select all / Clear */}
      {hasFiles && (
        <>
          <button
            type="button"
            onClick={handleToggleSelectAll}
            className="text-xs sm:text-sm text-muted-foreground hover:text-foreground transition-colors hidden xs:block"
          >
            {allSelected ? "Deselect" : "Select all"}
          </button>

          <Tooltip content="Clear all files" position="bottom">
            <IconButton
              size="sm"
              variant="ghost"
              onClick={clearFiles}
              className="text-muted-foreground hover:text-error hover:bg-error/10"
            >
              <TrashIcon className="size-4" />
            </IconButton>
          </Tooltip>
        </>
      )}

      {/* File count */}
      <span className="text-xs sm:text-sm text-muted-foreground tabular-nums">
        {files.length}
        <span className="hidden xs:inline">
          {" "}
          file{files.length !== 1 ? "s" : ""}
        </span>
      </span>
    </div>
  );
}
