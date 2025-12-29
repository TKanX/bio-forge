/**
 * @file Topology step
 *
 * Pipeline step for topology generation with template management.
 */

"use client";

import { useCallback } from "react";
import { useShallow } from "zustand/react/shallow";
import { Input } from "@/ui/primitives";
import { Dropzone } from "@/ui/patterns";
import { LayersIcon, TrashIcon, FileIcon } from "@/ui/icons";
import {
  usePipelineStore,
  useFileStore,
  showSuccess,
  showError,
} from "@/state";
import { DEFAULT_TOPOLOGY_SETTINGS } from "@/core/pipeline";
import { parseUploadedFiles } from "@/core/file";
import { ACCEPTED_TEMPLATE_EXTENSIONS, formatFileSize } from "@/lib";
import { StepWrapper } from "./step-wrapper";

// ============================================================================
// Component
// ============================================================================

export function StepTopology() {
  const { enabled, setEnabled, config, setConfig } = usePipelineStore(
    useShallow((s) => ({
      enabled: s.topologyEnabled,
      setEnabled: s.setTopologyEnabled,
      config: s.topologyConfig,
      setConfig: s.setTopologyConfig,
    }))
  );

  const { templates, addTemplates, removeTemplate, clearTemplates } =
    useFileStore(
      useShallow((s) => ({
        templates: s.templates,
        addTemplates: s.addTemplates,
        removeTemplate: s.removeTemplate,
        clearTemplates: s.clearTemplates,
      }))
    );

  const handleReset = () => {
    setConfig(DEFAULT_TOPOLOGY_SETTINGS);
    clearTemplates();
  };

  const handleTemplateDrop = useCallback(
    async (files: File[]) => {
      try {
        const { templates: parsedTemplates } = await parseUploadedFiles(files);
        if (parsedTemplates.length > 0) {
          addTemplates(parsedTemplates);
          showSuccess(`Added ${parsedTemplates.length} template(s)`);
        }
      } catch (error) {
        showError(
          error instanceof Error ? error.message : "Failed to add templates"
        );
      }
    },
    [addTemplates]
  );

  return (
    <StepWrapper
      icon={<LayersIcon className="size-4" />}
      title="Topology"
      enabled={enabled}
      onToggle={setEnabled}
      onReset={handleReset}
    >
      {/* Disulfide cutoff */}
      <Input
        label="Disulfide cutoff (Å)"
        type="number"
        step={0.1}
        min={0}
        value={config.disulfideCutoff}
        onChange={(e) =>
          setConfig({
            disulfideCutoff: parseFloat(e.target.value) || 2.2,
          })
        }
      />

      {/* Template dropzone */}
      <div className="space-y-2">
        <label className="text-xs font-medium text-muted-foreground">
          Custom templates (MOL2)
        </label>
        <Dropzone
          onDrop={handleTemplateDrop}
          accept={ACCEPTED_TEMPLATE_EXTENSIONS}
          multiple
          className="py-4"
        >
          <div className="flex flex-col items-center gap-2 text-center">
            <div className="p-2 rounded-full bg-muted text-muted-foreground">
              <FileIcon className="size-4" />
            </div>
            <p className="text-xs text-muted-foreground">
              Drop MOL2 template files here
            </p>
          </div>
        </Dropzone>
      </div>

      {/* Template list */}
      {templates.length > 0 && (
        <div className="space-y-1.5">
          {templates.map((template) => (
            <TemplateItem
              key={template.id}
              name={template.name}
              size={template.size}
              onRemove={() => removeTemplate(template.id)}
            />
          ))}
        </div>
      )}
    </StepWrapper>
  );
}

// ============================================================================
// Template Item
// ============================================================================

interface TemplateItemProps {
  name: string;
  size: number;
  onRemove: () => void;
}

function TemplateItem({ name, size, onRemove }: TemplateItemProps) {
  return (
    <div className="flex items-center gap-2 px-2 py-1.5 rounded bg-muted/50 text-sm">
      <FileIcon className="size-3.5 text-muted-foreground shrink-0" />
      <span className="flex-1 truncate">{name}</span>
      <span className="text-xs text-muted-foreground shrink-0">
        {formatFileSize(size)}
      </span>
      <button
        type="button"
        onClick={onRemove}
        className="p-0.5 rounded hover:bg-error/10 hover:text-error transition-colors"
      >
        <TrashIcon className="size-3.5" />
      </button>
    </div>
  );
}
