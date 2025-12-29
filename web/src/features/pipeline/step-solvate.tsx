/**
 * @file Solvate step
 *
 * Pipeline step for solvation (adding water box and ions).
 */

"use client";

import { useShallow } from "zustand/react/shallow";
import { Checkbox, Input, Select } from "@/ui/primitives";
import { DropletIcon } from "@/ui/icons";
import { usePipelineStore } from "@/state";
import {
  DEFAULT_SOLVATE_SETTINGS,
  type CationSpecies,
  type AnionSpecies,
} from "@/core/pipeline";
import { StepWrapper } from "./step-wrapper";

// ============================================================================
// Options
// ============================================================================

const CATION_OPTIONS = [
  { value: "Na", label: "Na⁺" },
  { value: "K", label: "K⁺" },
  { value: "Mg", label: "Mg²⁺" },
  { value: "Ca", label: "Ca²⁺" },
  { value: "Li", label: "Li⁺" },
  { value: "Zn", label: "Zn²⁺" },
] as const;

const ANION_OPTIONS = [
  { value: "Cl", label: "Cl⁻" },
  { value: "Br", label: "Br⁻" },
  { value: "I", label: "I⁻" },
  { value: "F", label: "F⁻" },
] as const;

// ============================================================================
// Component
// ============================================================================

export function StepSolvate() {
  const { enabled, setEnabled, config, setConfig } = usePipelineStore(
    useShallow((s) => ({
      enabled: s.solvateEnabled,
      setEnabled: s.setSolvateEnabled,
      config: s.solvateConfig,
      setConfig: s.setSolvateConfig,
    }))
  );

  return (
    <StepWrapper
      icon={<DropletIcon className="size-4" />}
      title="Solvate"
      enabled={enabled}
      onToggle={setEnabled}
      onReset={() => setConfig(DEFAULT_SOLVATE_SETTINGS)}
    >
      {/* Box parameters */}
      <div className="grid grid-cols-2 gap-3">
        <Input
          label="Margin (Å)"
          type="number"
          step={0.5}
          min={0}
          value={config.margin}
          onChange={(e) =>
            setConfig({ margin: parseFloat(e.target.value) || 10 })
          }
        />
        <Input
          label="Water spacing (Å)"
          type="number"
          step={0.1}
          min={0}
          value={config.waterSpacing}
          onChange={(e) =>
            setConfig({ waterSpacing: parseFloat(e.target.value) || 3.1 })
          }
        />
      </div>

      {/* VDW cutoff */}
      <Input
        label="VDW cutoff (Å)"
        type="number"
        step={0.1}
        min={0}
        value={config.vdwCutoff}
        onChange={(e) =>
          setConfig({ vdwCutoff: parseFloat(e.target.value) || 2.4 })
        }
      />

      {/* Remove existing */}
      <Checkbox
        label="Remove existing solvent"
        checked={config.removeExisting}
        onChange={(e) => setConfig({ removeExisting: e.target.checked })}
      />

      {/* Ion selection */}
      <div className="grid grid-cols-2 gap-3">
        <Select
          label="Cation"
          value={config.cations[0] || "Na"}
          onChange={(e) =>
            setConfig({ cations: [e.target.value as CationSpecies] })
          }
          options={[...CATION_OPTIONS]}
        />
        <Select
          label="Anion"
          value={config.anions[0] || "Cl"}
          onChange={(e) =>
            setConfig({ anions: [e.target.value as AnionSpecies] })
          }
          options={[...ANION_OPTIONS]}
        />
      </div>

      {/* Target charge */}
      <Input
        label="Target charge"
        type="number"
        value={config.targetCharge}
        onChange={(e) =>
          setConfig({ targetCharge: parseInt(e.target.value) || 0 })
        }
      />

      {/* RNG seed */}
      <Input
        label="RNG seed (optional)"
        type="number"
        placeholder="Leave empty for random"
        value={config.rngSeed ?? ""}
        onChange={(e) =>
          setConfig({
            rngSeed: e.target.value ? parseInt(e.target.value) : undefined,
          })
        }
      />
    </StepWrapper>
  );
}
