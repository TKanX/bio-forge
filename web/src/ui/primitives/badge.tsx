/**
 * @file Badge component
 *
 * Small status indicator primitive with optional icon support.
 */

import { forwardRef, type HTMLAttributes, type ReactNode } from "react";
import { cn } from "@/lib";

// ============================================================================
// Types
// ============================================================================

export type BadgeVariant =
  | "default"
  | "primary"
  | "success"
  | "warning"
  | "error"
  | "info"
  | "protein"
  | "nucleic"
  | "solvent"
  | "hetero";

export interface BadgeProps extends HTMLAttributes<HTMLSpanElement> {
  variant?: BadgeVariant;
  icon?: ReactNode;
}

// ============================================================================
// Styles
// ============================================================================

const variantStyles: Record<BadgeVariant, string> = {
  default: "bg-muted text-muted-foreground",
  primary: "bg-primary/20 text-primary",
  success: "bg-success/20 text-success",
  warning: "bg-warning/20 text-warning",
  error: "bg-error/20 text-error",
  info: "bg-info/20 text-info",
  protein: "bg-protein/20 text-protein",
  nucleic: "bg-nucleic/20 text-nucleic",
  solvent: "bg-solvent/20 text-solvent",
  hetero: "bg-hetero/20 text-hetero",
};

// ============================================================================
// Component
// ============================================================================

export const Badge = forwardRef<HTMLSpanElement, BadgeProps>(
  ({ variant = "default", icon, className, children, ...props }, ref) => {
    return (
      <span
        ref={ref}
        className={cn(
          "inline-flex items-center justify-center gap-1",
          "px-2 py-0.5 rounded-md",
          "text-xs font-medium",
          variantStyles[variant],
          className
        )}
        {...props}
      >
        {icon && <span className="size-3 shrink-0">{icon}</span>}
        {children}
      </span>
    );
  }
);

Badge.displayName = "Badge";
