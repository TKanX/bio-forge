/**
 * @file Empty state component
 *
 * Empty state placeholder with action.
 */

import { type ReactNode } from "react";
import { cn } from "@/lib";

// ============================================================================
// Types
// ============================================================================

export interface EmptyStateProps {
  /** Icon to display */
  icon?: ReactNode;
  /** Title text */
  title: string;
  /** Description text */
  description?: string;
  /** Action button */
  action?: ReactNode;
  /** Additional class names */
  className?: string;
}

// ============================================================================
// Component
// ============================================================================

export function EmptyState({
  icon,
  title,
  description,
  action,
  className,
}: EmptyStateProps) {
  return (
    <div
      className={cn(
        "flex flex-col items-center justify-center text-center p-8",
        className
      )}
    >
      {/* Icon */}
      {icon && (
        <div className="p-4 rounded-full bg-muted text-muted-foreground mb-4">
          {icon}
        </div>
      )}

      {/* Title */}
      <h3 className="text-lg font-medium text-foreground">{title}</h3>

      {/* Description */}
      {description && (
        <p className="text-sm text-muted-foreground mt-1 max-w-sm">
          {description}
        </p>
      )}

      {/* Action */}
      {action && <div className="mt-4">{action}</div>}
    </div>
  );
}
