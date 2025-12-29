/**
 * @file useFileUpload hook
 *
 * Unified hook for handling file upload.
 */

import { useCallback } from "react";
import { useShallow } from "zustand/react/shallow";
import { useFileStore, showError, showInfo, showSuccess } from "../stores";
import { parseUploadedFiles } from "@/core";

/**
 * Hook providing file upload functionality.
 *
 * Files are parsed immediately during upload, so they arrive
 * in the store already in "ready" state with structure info.
 *
 * @returns Object with `handleFiles` function to process uploaded files.
 */
export function useFileUpload() {
  const { addFiles, addTemplates } = useFileStore(
    useShallow((s) => ({
      addFiles: s.addFiles,
      addTemplates: s.addTemplates,
    }))
  );

  /**
   * Handle file drop/upload.
   *
   * @param files - Array of File objects from input or drop event
   * @param options - Configuration options
   * @param options.showMessages - Whether to show toast messages (default: true)
   */
  const handleFiles = useCallback(
    async (files: File[], options?: { showMessages?: boolean }) => {
      const showMessages = options?.showMessages !== false;
      const { structures, templates, errors } = await parseUploadedFiles(files);

      // Show errors
      if (showMessages) {
        for (const { name, error } of errors) {
          showError(`${name}: ${error}`);
        }
      }

      // Add templates (already parsed with WASM Template objects)
      if (templates.length > 0) {
        addTemplates(templates);
        if (showMessages) {
          showInfo(
            `Added ${templates.length} template${templates.length !== 1 ? "s" : ""}`
          );
        }
      }

      // Add files (already parsed with WASM Structure objects and info)
      if (structures.length > 0) {
        addFiles(structures);
        if (showMessages) {
          showSuccess(
            `Added ${structures.length} file${structures.length !== 1 ? "s" : ""}`
          );
        }
      }

      return { structures, templates, errors };
    },
    [addFiles, addTemplates]
  );

  return { handleFiles };
}
