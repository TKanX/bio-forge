/**
 * @file File list
 *
 * Container component for displaying all uploaded files.
 */

"use client";

import { AnimatePresence } from "framer-motion";
import { useFileStore, useFileUpload } from "@/state";
import { EmptyState, Dropzone } from "@/ui/patterns";
import { FileIcon, UploadIcon } from "@/ui/icons";
import { ACCEPTED_STRUCTURE_EXTENSIONS } from "@/lib/constants";
import { FileCard } from "./file-card";
import { FileActions } from "./file-actions";

// ============================================================================
// Component
// ============================================================================

export function FileList() {
  const files = useFileStore((s) => s.files);
  const isEmpty = files.length === 0;

  return (
    <div className="h-full flex flex-col overflow-hidden">
      {/* Header with actions */}
      <FileActions />

      {/* File list or empty state */}
      <div className="flex-1 overflow-y-auto p-4">
        {isEmpty ? (
          <FileEmptyState />
        ) : (
          <div className="space-y-3">
            <AnimatePresence mode="popLayout">
              {files.map((file) => (
                <FileCard key={file.id} file={file} />
              ))}
            </AnimatePresence>
          </div>
        )}
      </div>
    </div>
  );
}

// ============================================================================
// Empty State
// ============================================================================

function FileEmptyState() {
  return (
    <EmptyState
      icon={<FileIcon className="size-8" />}
      title="No files uploaded"
      description="Upload structure files to start processing"
      action={
        <div className="mt-4 w-full max-w-md">
          <UploadDropzone />
        </div>
      }
      className="h-full"
    />
  );
}

// ============================================================================
// Upload Dropzone (for empty state)
// ============================================================================

function UploadDropzone() {
  const { handleFiles } = useFileUpload();

  return (
    <Dropzone
      onDrop={handleFiles}
      accept={ACCEPTED_STRUCTURE_EXTENSIONS}
      multiple
      className="min-h-30"
    >
      <div className="flex flex-col items-center justify-center gap-3 p-6 text-center">
        <div className="p-3 rounded-full bg-primary/10 text-primary">
          <UploadIcon className="size-6" />
        </div>
        <div className="space-y-1">
          <p className="text-sm font-medium text-foreground">
            Drop files here or click to browse
          </p>
          <p className="text-xs text-muted-foreground">
            Supports PDB (.pdb, .ent) and mmCIF (.cif, .mmcif)
          </p>
        </div>
      </div>
    </Dropzone>
  );
}
