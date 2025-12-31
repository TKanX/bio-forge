/**
 * @file App layout
 *
 * Layout for the main application page with sidebar and content area.
 */

"use client";

import { type ReactNode } from "react";
import { useUIStore, useFileUpload } from "@/state";
import { useKeyboardShortcuts } from "@/state/hooks";
import { GlobalDropzone } from "@/ui/patterns";
import { ACCEPTED_STRUCTURE_EXTENSIONS, cn } from "@/lib";

// ============================================================================
// Types
// ============================================================================

interface AppLayoutProps {
  /** Header component */
  header: ReactNode;
  /** Sidebar component */
  sidebar: ReactNode;
  /** Main content */
  children: ReactNode;
}

// ============================================================================
// Component
// ============================================================================

export function AppLayout({ header, sidebar, children }: AppLayoutProps) {
  const sidebarOpen = useUIStore((s) => s.sidebarOpen);
  const { handleFiles } = useFileUpload();

  // Enable global keyboard shortcuts
  useKeyboardShortcuts();

  return (
    <div className="h-screen flex flex-col overflow-hidden bg-background">
      {/* Global drop zone overlay */}
      <GlobalDropzone
        onDrop={handleFiles}
        accept={ACCEPTED_STRUCTURE_EXTENSIONS}
      />

      {/* Header */}
      {header}

      {/* Main content area */}
      <div className="flex-1 flex overflow-hidden">
        {/* Sidebar */}
        <aside
          className={cn(
            "shrink-0 border-r border-border bg-card/50",
            "transition-all duration-300 ease-in-out",
            "w-80",
            // Mobile responsive
            "max-lg:absolute max-lg:top-14 max-lg:sm:top-16 max-lg:bottom-0 max-lg:left-0 max-lg:z-30",
            !sidebarOpen && "max-lg:-translate-x-full"
          )}
        >
          {sidebar}
        </aside>

        {/* Mobile sidebar backdrop */}
        {sidebarOpen && (
          <div
            className="lg:hidden fixed inset-0 bg-black/50 z-20"
            onClick={() => useUIStore.getState().toggleSidebar()}
          />
        )}

        {/* Content */}
        <main className="flex-1 overflow-hidden">{children}</main>
      </div>
    </div>
  );
}
