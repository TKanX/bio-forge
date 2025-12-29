/**
 * @file App Page
 *
 * Main application page for structure processing.
 */

"use client";

import { AppLayout } from "@/layouts";
import { AppHeader } from "@/features/header";
import { PipelinePanel } from "@/features/pipeline";
import { FileList } from "@/features/files";

// ============================================================================
// Page
// ============================================================================

export default function AppPage() {
  return (
    <AppLayout header={<AppHeader />} sidebar={<PipelinePanel />}>
      <FileList />
    </AppLayout>
  );
}
