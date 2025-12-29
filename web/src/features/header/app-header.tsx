/**
 * @file App header
 *
 * Application header with logo, navigation, and action buttons.
 */

"use client";

import Link from "next/link";
import {
  Button,
  Dropdown,
  IconButton,
  type DropdownItem,
} from "@/ui/primitives";
import {
  PlayIcon,
  DownloadIcon,
  LoaderIcon,
  GithubIcon,
  BoxIcon,
  ArchiveIcon,
  FileTextIcon,
  MenuIcon,
  HomeIcon,
} from "@/ui/icons";
import { showError, useUIStore } from "@/state";
import { useProcessor } from "@/state/hooks";
import { GITHUB_URL } from "@/lib";

// ============================================================================
// Component
// ============================================================================

export function AppHeader() {
  const toggleSidebar = useUIStore((s) => s.toggleSidebar);

  const {
    isProcessing,
    canProcess,
    canDownload,
    hasAnyStepEnabled,
    filesToProcess,
    completedFiles,
    processFiles,
    downloadFiles,
  } = useProcessor();

  const handleRun = async () => {
    if (!canProcess) {
      if (!hasAnyStepEnabled) {
        showError("Enable at least one pipeline step");
      } else {
        showError("No files ready to process");
      }
      return;
    }

    try {
      await processFiles();
    } catch (error) {
      showError(error instanceof Error ? error.message : "Processing failed");
    }
  };

  const handleDownload = async (value: string) => {
    if (!canDownload) {
      showError("No completed files to download");
      return;
    }

    try {
      const [format, asZip] = value.split("-") as ["pdb" | "mmcif", string?];
      await downloadFiles(format, asZip === "zip");
    } catch (error) {
      showError(error instanceof Error ? error.message : "Download failed");
    }
  };

  const downloadItems: DropdownItem[] = [
    {
      label: "Download as PDB",
      value: "pdb",
      icon: <FileTextIcon className="size-4" />,
    },
    {
      label: "Download as mmCIF",
      value: "mmcif",
      icon: <FileTextIcon className="size-4" />,
    },
    {
      label: "Download as PDB (ZIP)",
      value: "pdb-zip",
      icon: <ArchiveIcon className="size-4" />,
    },
    {
      label: "Download as mmCIF (ZIP)",
      value: "mmcif-zip",
      icon: <ArchiveIcon className="size-4" />,
    },
  ];

  return (
    <header className="h-14 sm:h-16 border-b border-border bg-card/80 backdrop-blur-sm flex items-center justify-between px-3 sm:px-6 shrink-0">
      {/* Left - Menu button (mobile), Logo and nav */}
      <div className="flex items-center gap-2 sm:gap-6">
        {/* Mobile menu button */}
        <IconButton
          variant="ghost"
          size="sm"
          onClick={toggleSidebar}
          className="lg:hidden"
          aria-label="Toggle pipeline settings"
        >
          <MenuIcon className="size-5" />
        </IconButton>

        <Link href="/" className="flex items-center gap-2 sm:gap-3 group">
          <div className="size-8 sm:size-9 rounded-lg bg-primary/10 flex items-center justify-center group-hover:bg-primary/20 transition-colors">
            <BoxIcon className="size-4 sm:size-5 text-primary" />
          </div>
          <span className="text-lg sm:text-xl font-bold hidden xs:inline">
            <span className="text-primary">Bio</span>Forge
          </span>
        </Link>

        <nav className="hidden md:flex items-center gap-1">
          <Link href="/">
            <Button variant="ghost" size="sm">
              <HomeIcon className="size-4" />
              Home
            </Button>
          </Link>
          <a href={GITHUB_URL} target="_blank" rel="noopener noreferrer">
            <Button variant="ghost" size="sm">
              <GithubIcon className="size-4" />
              GitHub
            </Button>
          </a>
        </nav>
      </div>

      {/* Right - Actions */}
      <div className="flex items-center gap-2 sm:gap-3">
        {/* Run button */}
        <Button
          onClick={handleRun}
          disabled={!canProcess}
          size="sm"
          className="px-2 sm:px-3"
        >
          {isProcessing ? (
            <>
              <LoaderIcon className="size-4 animate-spin" />
              <span className="hidden sm:inline">Processing...</span>
            </>
          ) : (
            <>
              <PlayIcon className="size-4" />
              <span className="hidden sm:inline">Run</span>
              <span className="sm:hidden">{filesToProcess.length}</span>
              <span className="hidden sm:inline">
                ({filesToProcess.length})
              </span>
            </>
          )}
        </Button>

        {/* Download button */}
        <Dropdown
          trigger={
            <Button
              variant="secondary"
              disabled={!canDownload}
              size="sm"
              className="px-2 sm:px-3"
            >
              <DownloadIcon className="size-4" />
              <span className="hidden sm:inline">
                Download ({completedFiles.length})
              </span>
              <span className="sm:hidden">{completedFiles.length}</span>
            </Button>
          }
          items={downloadItems}
          onSelect={handleDownload}
        />
      </div>
    </header>
  );
}
