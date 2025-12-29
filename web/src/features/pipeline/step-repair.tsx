/**
 * @file Repair step
 *
 * Pipeline step for repairing missing atoms.
 */

"use client";

import { WrenchIcon } from "@/ui/icons";
import { usePipelineStore } from "@/state";
import { StepWrapper } from "./step-wrapper";

// ============================================================================
// Component
// ============================================================================

export function StepRepair() {
  const enabled = usePipelineStore((s) => s.repairEnabled);
  const setEnabled = usePipelineStore((s) => s.setRepairEnabled);

  return (
    <StepWrapper
      variant="simple"
      icon={<WrenchIcon className="size-4" />}
      title="Repair"
      enabled={enabled}
      onToggle={setEnabled}
    />
  );
}
