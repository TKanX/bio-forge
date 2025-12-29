/**
 * @file Global Providers
 *
 * Application-wide providers including WASM initialization and toast container.
 */

"use client";

import { type ReactNode } from "react";
import { ToastContainer } from "@/ui/patterns";
import { useWasmInit } from "@/state/hooks";

// ============================================================================
// Component
// ============================================================================

export function Providers({ children }: { children: ReactNode }) {
  // Initialize WASM module on app load
  useWasmInit();

  return (
    <>
      {children}
      <ToastContainer />
    </>
  );
}
