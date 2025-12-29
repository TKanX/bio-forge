/**
 * @file Input component
 *
 * Text input primitive with label support.
 */

import { forwardRef, type InputHTMLAttributes } from "react";
import { cn } from "@/lib";

// ============================================================================
// Types
// ============================================================================

export interface InputProps extends InputHTMLAttributes<HTMLInputElement> {
  label?: string;
  error?: string;
}

// ============================================================================
// Component
// ============================================================================

export const Input = forwardRef<HTMLInputElement, InputProps>(
  ({ label, error, className, id, ...props }, ref) => {
    const inputId = id ?? label?.toLowerCase().replace(/\s+/g, "-");

    return (
      <div className="space-y-1.5">
        {label && (
          <label
            htmlFor={inputId}
            className="block text-xs font-medium text-muted-foreground"
          >
            {label}
          </label>
        )}
        <input
          ref={ref}
          id={inputId}
          className={cn(
            // Base styles
            "w-full h-9 px-3",
            "rounded-lg border border-border bg-input",
            "text-sm text-foreground placeholder:text-muted-foreground",
            "transition-colors duration-150",
            // Focus
            "focus:outline-none focus:ring-2 focus:ring-primary focus:border-transparent",
            // Disabled
            "disabled:opacity-50 disabled:cursor-not-allowed",
            // Error
            error && "border-error focus:ring-error",
            className
          )}
          {...props}
        />
        {error && <p className="text-xs text-error">{error}</p>}
      </div>
    );
  }
);

Input.displayName = "Input";
