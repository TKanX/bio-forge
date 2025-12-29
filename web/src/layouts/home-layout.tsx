/**
 * @file Home layout
 *
 * Layout for the landing page.
 */

import { type ReactNode } from "react";

// ============================================================================
// Types
// ============================================================================

interface HomeLayoutProps {
  children: ReactNode;
}

// ============================================================================
// Component
// ============================================================================

export function HomeLayout({ children }: HomeLayoutProps) {
  return <div className="min-h-screen bg-background">{children}</div>;
}
