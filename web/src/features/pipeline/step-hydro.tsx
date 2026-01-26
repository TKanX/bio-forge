/**
 * @file Hydro step
 *
 * Pipeline step for adding hydrogens.
 */

"use client";

import { useShallow } from "zustand/react/shallow";
import { Checkbox, Input, Select } from "@/ui/primitives";
import { BeakerIcon } from "@/ui/icons";
import { usePipelineStore } from "@/state";
import { DEFAULT_HYDRO_SETTINGS, type HisStrategy } from "@/core/pipeline";
import { StepWrapper } from "./step-wrapper";

// ============================================================================
// Options
// ============================================================================

const HIS_STRATEGY_OPTIONS = [
  { value: "network", label: "Network (H-bond aware)" },
  { value: "hid", label: "HID (δ-protonated)" },
  { value: "hie", label: "HIE (ε-protonated)" },
  { value: "random", label: "Random" },
] as const;

// ============================================================================
// Component
// ============================================================================

export function StepHydro() {
  const { enabled, setEnabled, config, setConfig } = usePipelineStore(
    useShallow((s) => ({
      enabled: s.hydroEnabled,
      setEnabled: s.setHydroEnabled,
      config: s.hydroConfig,
      setConfig: s.setHydroConfig,
    }))
  );

  return (
    <StepWrapper
      icon={<BeakerIcon className="size-4" />}
      title="Add Hydrogens"
      enabled={enabled}
      onToggle={setEnabled}
      onReset={() => setConfig(DEFAULT_HYDRO_SETTINGS)}
    >
      {/* pH setting */}
      <Input
        label="Target pH"
        description="Omit to preserve original residue names (no auto-protonation)"
        type="number"
        step="0.1"
        min={0}
        max={14}
        placeholder="Leave empty to preserve names"
        value={config.targetPh ?? ""}
        onChange={(e) =>
          setConfig({
            targetPh: e.target.value ? parseFloat(e.target.value) : undefined,
          })
        }
      />

      {/* Remove existing H option */}
      <Checkbox
        label="Remove existing hydrogens"
        checked={config.removeExistingH}
        onChange={(e) => setConfig({ removeExistingH: e.target.checked })}
      />

      {/* Histidine strategy */}
      <Select
        label="Histidine strategy"
        description="Only applies when pH is specified and HIS is neutral (pH ≥ 6.0)"
        value={config.hisStrategy}
        onChange={(e) =>
          setConfig({
            hisStrategy: e.target.value as HisStrategy,
          })
        }
        options={[...HIS_STRATEGY_OPTIONS]}
      />

      {/* Salt bridge detection */}
      <Checkbox
        label="Detect salt bridges"
        description="HIS → HIP near carboxylate groups (ASP⁻/GLU⁻/C-term COO⁻). pH-independent geometric detection."
        checked={config.hisSaltBridgeProtonation}
        onChange={(e) =>
          setConfig({ hisSaltBridgeProtonation: e.target.checked })
        }
      />
    </StepWrapper>
  );
}
