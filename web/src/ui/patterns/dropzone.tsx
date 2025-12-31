/**
 * @file Dropzone component
 *
 * File drop zone with drag-and-drop support.
 */

"use client";

import {
  useState,
  useRef,
  type ReactNode,
  type DragEvent,
  type ChangeEvent,
} from "react";
import { motion, AnimatePresence } from "framer-motion";
import { cn } from "@/lib";

// ============================================================================
// Types
// ============================================================================

export interface DropzoneProps {
  /** Callback when files are dropped */
  onDrop: (files: File[]) => void;
  /** Accepted file extensions (e.g., [".pdb", ".mol2"]) */
  accept?: string[];
  /** Allow multiple files */
  multiple?: boolean;
  /** Disabled state */
  disabled?: boolean;
  /** Custom content */
  children?: ReactNode;
  /** Additional class names */
  className?: string;
}

// ============================================================================
// Component
// ============================================================================

export function Dropzone({
  onDrop,
  accept,
  multiple = true,
  disabled = false,
  children,
  className,
}: DropzoneProps) {
  const [isDragging, setIsDragging] = useState(false);
  const inputRef = useRef<HTMLInputElement>(null);

  const acceptString = accept?.join(",");

  const handleDragEnter = (e: DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    if (!disabled) {
      setIsDragging(true);
    }
  };

  const handleDragLeave = (e: DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    setIsDragging(false);
  };

  const handleDragOver = (e: DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
  };

  const handleDrop = (e: DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    setIsDragging(false);

    if (disabled) return;

    const files = Array.from(e.dataTransfer.files);
    const filtered = filterFiles(files, accept);

    if (filtered.length > 0) {
      onDrop(multiple ? filtered : [filtered[0]]);
    }
  };

  const handleInputChange = (e: ChangeEvent<HTMLInputElement>) => {
    const files = Array.from(e.target.files || []);
    if (files.length > 0) {
      onDrop(files);
    }
    e.target.value = "";
  };

  const handleClick = () => {
    if (!disabled) {
      inputRef.current?.click();
    }
  };

  return (
    <div
      role="button"
      tabIndex={disabled ? -1 : 0}
      onClick={handleClick}
      onKeyDown={(e) => e.key === "Enter" && handleClick()}
      onDragEnter={handleDragEnter}
      onDragLeave={handleDragLeave}
      onDragOver={handleDragOver}
      onDrop={handleDrop}
      className={cn(
        "relative rounded-lg border-2 border-dashed",
        "transition-colors duration-150",
        "cursor-pointer",
        disabled
          ? "opacity-50 cursor-not-allowed border-border bg-muted"
          : isDragging
            ? "border-primary bg-primary/5"
            : "border-border hover:border-primary/50 hover:bg-muted/50",
        className
      )}
    >
      {/* Hidden file input */}
      <input
        ref={inputRef}
        type="file"
        accept={acceptString}
        multiple={multiple}
        disabled={disabled}
        onChange={handleInputChange}
        className="hidden"
      />

      {/* Content */}
      {children || <DefaultContent isDragging={isDragging} accept={accept} />}

      {/* Drag overlay */}
      <AnimatePresence>
        {isDragging && (
          <motion.div
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            exit={{ opacity: 0 }}
            className="absolute inset-0 rounded-lg bg-primary/10 flex items-center justify-center"
          >
            <div className="text-primary font-medium">Drop files here</div>
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
}

// ============================================================================
// Default Content
// ============================================================================

interface DefaultContentProps {
  isDragging: boolean;
  accept?: string[];
}

function DefaultContent({ isDragging, accept }: DefaultContentProps) {
  return (
    <div className="flex flex-col items-center justify-center gap-3 p-8 text-center">
      <div
        className={cn(
          "p-3 rounded-full transition-colors duration-150",
          isDragging
            ? "bg-primary/20 text-primary"
            : "bg-muted text-muted-foreground"
        )}
      >
        <svg
          className="size-6"
          viewBox="0 0 24 24"
          fill="none"
          stroke="currentColor"
          strokeWidth="2"
        >
          <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4" />
          <polyline points="17 8 12 3 7 8" />
          <line x1="12" y1="3" x2="12" y2="15" />
        </svg>
      </div>

      <div className="space-y-1">
        <p className="text-sm font-medium text-foreground">
          Drop files here or click to browse
        </p>
        {accept && (
          <p className="text-xs text-muted-foreground">
            Supported formats: {accept.join(", ")}
          </p>
        )}
      </div>
    </div>
  );
}

// ============================================================================
// Helpers
// ============================================================================

/**
 * Filter files by accepted extensions.
 *
 * @param files - Array of files to filter
 * @param accept - Array of accepted extensions (e.g., [".pdb", ".mol2"])
 * @returns Filtered array of files matching accepted extensions
 */
export function filterFiles(files: File[], accept?: string[]): File[] {
  if (!accept || accept.length === 0) return files;

  return files.filter((file) => {
    const ext = "." + file.name.split(".").pop()?.toLowerCase();
    return accept.some((a) => a.toLowerCase() === ext);
  });
}
