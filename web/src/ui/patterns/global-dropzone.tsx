/**
 * @file Global dropzone overlay
 *
 * Full-screen drag-and-drop overlay that appears when dragging files
 * anywhere on the page. Provides a reliable drop target for file uploads.
 */

"use client";

import { useState, useEffect, useCallback, type ReactNode } from "react";
import { motion, AnimatePresence } from "framer-motion";
import { UploadIcon } from "@/ui/icons";
import { cn } from "@/lib";
import { filterFiles } from "./dropzone";

// ============================================================================
// Types
// ============================================================================

export interface GlobalDropzoneProps {
  /** Callback when files are dropped */
  onDrop: (files: File[]) => void;
  /** Accepted file extensions (e.g., [".pdb", ".mol2"]) */
  accept?: string[];
  /** Allow multiple files */
  multiple?: boolean;
  /** Disabled state */
  disabled?: boolean;
  /** Custom overlay content */
  children?: ReactNode;
  /** Additional class names for overlay */
  className?: string;
}

// ============================================================================
// Hooks
// ============================================================================

/**
 * Hook to track global drag state with proper enter/leave counting.
 */
function useGlobalDragState(disabled: boolean) {
  const [isDragging, setIsDragging] = useState(false);
  const [, setDragCounter] = useState(0);

  useEffect(() => {
    if (disabled) return;

    const handleDragEnter = (e: DragEvent) => {
      if (!e.dataTransfer?.types?.includes("Files")) return;

      e.preventDefault();
      setDragCounter((prev) => {
        const next = prev + 1;
        if (next === 1) setIsDragging(true);
        return next;
      });
    };

    const handleDragLeave = (e: DragEvent) => {
      e.preventDefault();
      setDragCounter((prev) => {
        const next = prev - 1;
        if (next <= 0) {
          setIsDragging(false);
          return 0;
        }
        return next;
      });
    };

    const handleDragOver = (e: DragEvent) => {
      e.preventDefault();
      if (e.dataTransfer) {
        e.dataTransfer.dropEffect = "copy";
      }
    };

    const handleDrop = (e: DragEvent) => {
      e.preventDefault();
      setIsDragging(false);
      setDragCounter(0);
    };

    const handleDragEnd = () => {
      setIsDragging(false);
      setDragCounter(0);
    };

    document.addEventListener("dragenter", handleDragEnter);
    document.addEventListener("dragleave", handleDragLeave);
    document.addEventListener("dragover", handleDragOver);
    document.addEventListener("drop", handleDrop);
    window.addEventListener("blur", handleDragEnd);

    return () => {
      document.removeEventListener("dragenter", handleDragEnter);
      document.removeEventListener("dragleave", handleDragLeave);
      document.removeEventListener("dragover", handleDragOver);
      document.removeEventListener("drop", handleDrop);
      window.removeEventListener("blur", handleDragEnd);
    };
  }, [disabled]);

  const reset = useCallback(() => {
    setIsDragging(false);
    setDragCounter(0);
  }, []);

  return { isDragging, reset };
}

// ============================================================================
// Component
// ============================================================================

export function GlobalDropzone({
  onDrop,
  accept,
  multiple = true,
  disabled = false,
  children,
  className,
}: GlobalDropzoneProps) {
  const { isDragging, reset } = useGlobalDragState(disabled);

  const handleDrop = useCallback(
    (e: React.DragEvent) => {
      e.preventDefault();
      e.stopPropagation();
      reset();

      if (disabled) return;

      const files = Array.from(e.dataTransfer?.files || []);
      const filtered = filterFiles(files, accept);

      if (filtered.length > 0) {
        onDrop(multiple ? filtered : [filtered[0]]);
      }
    },
    [accept, disabled, multiple, onDrop, reset]
  );

  const handleDragOver = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
  }, []);

  return (
    <AnimatePresence>
      {isDragging && !disabled && (
        <motion.div
          initial={{ opacity: 0 }}
          animate={{ opacity: 1 }}
          exit={{ opacity: 0 }}
          transition={{ duration: 0.15 }}
          className={cn(
            "fixed inset-0 z-50",
            "bg-background/80 backdrop-blur-sm",
            "flex items-center justify-center",
            className
          )}
          onDragOver={handleDragOver}
          onDrop={handleDrop}
        >
          {children || <DefaultOverlayContent accept={accept} />}
        </motion.div>
      )}
    </AnimatePresence>
  );
}

// ============================================================================
// Default Overlay Content
// ============================================================================

interface DefaultOverlayContentProps {
  accept?: string[];
}

function DefaultOverlayContent({ accept }: DefaultOverlayContentProps) {
  const formatList = accept?.length
    ? accept.map((ext) => ext.replace(".", "").toUpperCase()).join(", ")
    : null;

  return (
    <motion.div
      initial={{ scale: 0.95, opacity: 0 }}
      animate={{ scale: 1, opacity: 1 }}
      exit={{ scale: 0.95, opacity: 0 }}
      transition={{ duration: 0.15 }}
      className="pointer-events-none"
    >
      <div className="flex flex-col items-center gap-4 p-8 rounded-2xl border-2 border-dashed border-primary bg-primary/5">
        <div className="p-4 rounded-full bg-primary/20">
          <UploadIcon className="size-8 text-primary" />
        </div>
        <div className="text-center">
          <p className="text-lg font-semibold text-foreground">
            Drop files to upload
          </p>
          {formatList && (
            <p className="text-sm text-muted-foreground mt-1">
              Supported formats: {formatList}
            </p>
          )}
        </div>
      </div>
    </motion.div>
  );
}
