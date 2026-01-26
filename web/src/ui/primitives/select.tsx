/**
 * @file Select component
 *
 * Dropdown select primitive with label support.
 */

import { forwardRef, type SelectHTMLAttributes } from "react";
import { cn } from "@/lib";

// ============================================================================
// Types
// ============================================================================

export interface SelectOption {
  value: string;
  label: string;
  disabled?: boolean;
}

export interface SelectProps extends Omit<
  SelectHTMLAttributes<HTMLSelectElement>,
  "children"
> {
  label?: string;
  description?: string;
  options: SelectOption[];
  placeholder?: string;
}

// ============================================================================
// Component
// ============================================================================

export const Select = forwardRef<HTMLSelectElement, SelectProps>(
  ({ label, description, options, placeholder, className, id, ...props }, ref) => {
    const selectId = id ?? label?.toLowerCase().replace(/\s+/g, "-");

    return (
      <div className="space-y-1.5">
        {label && (
          <label
            htmlFor={selectId}
            className="block text-xs font-medium text-muted-foreground"
          >
            {label}
          </label>
        )}
        <div className="relative">
          <select
            ref={ref}
            id={selectId}
            className={cn(
              // Base styles
              "w-full h-9 px-3 pr-8",
              "rounded-lg border border-border bg-input",
              "text-sm text-foreground",
              "transition-colors duration-150",
              "appearance-none cursor-pointer",
              // Focus
              "focus:outline-none focus:ring-2 focus:ring-primary focus:border-transparent",
              // Disabled
              "disabled:opacity-50 disabled:cursor-not-allowed",
              className
            )}
            {...props}
          >
            {placeholder && (
              <option value="" disabled>
                {placeholder}
              </option>
            )}
            {options.map((option) => (
              <option
                key={option.value}
                value={option.value}
                disabled={option.disabled}
              >
                {option.label}
              </option>
            ))}
          </select>
          {/* Chevron */}
          <svg
            className="absolute right-2.5 top-1/2 -translate-y-1/2 size-4 text-muted-foreground pointer-events-none"
            viewBox="0 0 24 24"
            fill="none"
            stroke="currentColor"
            strokeWidth="2"
            strokeLinecap="round"
            strokeLinejoin="round"
          >
            <path d="m6 9 6 6 6-6" />
          </svg>
        </div>
        {description && (
          <p className="text-xs text-muted-foreground/80">{description}</p>
        )}
      </div>
    );
  }
);

Select.displayName = "Select";
