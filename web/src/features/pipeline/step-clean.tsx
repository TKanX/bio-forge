/**
 * @file Clean step
 *
 * Pipeline step for cleaning structures.
 */

"use client";

import { useState, useEffect } from "react";
import { useShallow } from "zustand/react/shallow";
import { Checkbox, Input } from "@/ui/primitives";
import { SparklesIcon } from "@/ui/icons";
import { usePipelineStore } from "@/state";
import { DEFAULT_CLEAN_SETTINGS } from "@/core/pipeline";
import { StepWrapper } from "./step-wrapper";

// ============================================================================
// Helpers
// ============================================================================

function parseCommaSeparated(value: string): string[] {
  return value
    .split(",")
    .map((s) => s.trim())
    .filter(Boolean);
}

// ============================================================================
// Component
// ============================================================================

export function StepClean() {
  const { enabled, setEnabled, config, setConfig } = usePipelineStore(
    useShallow((s) => ({
      enabled: s.cleanEnabled,
      setEnabled: s.setCleanEnabled,
      config: s.cleanConfig,
      setConfig: s.setCleanConfig,
    }))
  );

  const [removeInput, setRemoveInput] = useState(
    config.removeResidueNames.join(", ")
  );
  const [keepInput, setKeepInput] = useState(
    config.keepResidueNames.join(", ")
  );

  useEffect(() => {
    setRemoveInput(config.removeResidueNames.join(", "));
  }, [config.removeResidueNames]);

  useEffect(() => {
    setKeepInput(config.keepResidueNames.join(", "));
  }, [config.keepResidueNames]);

  const handleRemoveBlur = () => {
    setConfig({ removeResidueNames: parseCommaSeparated(removeInput) });
  };

  const handleKeepBlur = () => {
    setConfig({ keepResidueNames: parseCommaSeparated(keepInput) });
  };

  return (
    <StepWrapper
      icon={<SparklesIcon className="size-4" />}
      title="Clean"
      enabled={enabled}
      onToggle={setEnabled}
      onReset={() => setConfig(DEFAULT_CLEAN_SETTINGS)}
    >
      {/* Removal options */}
      <Checkbox
        label="Remove water molecules"
        checked={config.removeWater}
        onChange={(e) => setConfig({ removeWater: e.target.checked })}
      />
      <Checkbox
        label="Remove ions"
        checked={config.removeIons}
        onChange={(e) => setConfig({ removeIons: e.target.checked })}
      />
      <Checkbox
        label="Remove hydrogens"
        checked={config.removeHydrogens}
        onChange={(e) => setConfig({ removeHydrogens: e.target.checked })}
      />
      <Checkbox
        label="Remove hetero residues"
        checked={config.removeHetero}
        onChange={(e) => setConfig({ removeHetero: e.target.checked })}
      />

      {/* Custom residue filters */}
      <div className="pt-2 border-t border-border/50 space-y-3">
        <Input
          label="Remove residue names"
          placeholder="NAG, MAN (comma separated)"
          value={removeInput}
          onChange={(e) => setRemoveInput(e.target.value)}
          onBlur={handleRemoveBlur}
        />

        <Input
          label="Keep residue names"
          placeholder="LIG, HEM (comma separated)"
          value={keepInput}
          onChange={(e) => setKeepInput(e.target.value)}
          onBlur={handleKeepBlur}
        />
      </div>
    </StepWrapper>
  );
}
