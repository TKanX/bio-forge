/**
 * @file useWasmInit hook
 *
 * Hook for initializing the WASM module.
 */

import { useEffect, useRef } from "react";
import { useUIStore, showSuccess, showError } from "../stores";
import { initWasm, isWasmReady } from "@/core";

/**
 * Hook for initializing WASM module on mount.
 */
export function useWasmInit() {
  const setWasmReady = useUIStore((s) => s.setWasmReady);
  const wasmReady = useUIStore((s) => s.wasmReady);
  const initRef = useRef(false);

  useEffect(() => {
    if (initRef.current || isWasmReady()) {
      if (isWasmReady()) {
        setWasmReady(true);
      }
      return;
    }
    initRef.current = true;

    initWasm()
      .then(() => {
        setWasmReady(true);
        showSuccess("BioForge WASM loaded", 2000);
      })
      .catch((err) => {
        initRef.current = false;
        showError(`Failed to load WASM: ${err.message}`, 0);
      });
  }, [setWasmReady]);

  return wasmReady;
}
