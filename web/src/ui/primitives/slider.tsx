/**
 * @file Slider component
 *
 * Range slider input for numeric values.
 */

"use client";

import { forwardRef, useId, type ComponentPropsWithoutRef } from "react";
import { cn } from "@/lib";

// ============================================================================
// Types
// ============================================================================

export interface SliderProps extends Omit<
  ComponentPropsWithoutRef<"input">,
  "type" | "size"
> {
  /** Label text */
  label?: string;
  /** Current value display */
  valueDisplay?: string;
  /** Size variant */
  size?: "sm" | "md";
}

// ============================================================================
// Styles
// ============================================================================

const sizeStyles = {
  sm: "h-1",
  md: "h-2",
} as const;

// ============================================================================
// Component
// ============================================================================

export const Slider = forwardRef<HTMLInputElement, SliderProps>(
  ({ label, valueDisplay, size = "md", className, id, ...props }, ref) => {
    const generatedId = useId();
    const inputId = id || generatedId;

    return (
      <div className="space-y-2">
        {/* Label row */}
        {(label || valueDisplay) && (
          <div className="flex items-center justify-between text-sm">
            {label && (
              <label htmlFor={inputId} className="text-foreground">
                {label}
              </label>
            )}
            {valueDisplay && (
              <span className="text-muted-foreground font-mono">
                {valueDisplay}
              </span>
            )}
          </div>
        )}

        {/* Slider */}
        <input
          ref={ref}
          type="range"
          id={inputId}
          className={cn(
            "w-full appearance-none rounded-full bg-border cursor-pointer",
            // Track
            "[&::-webkit-slider-runnable-track]:rounded-full",
            "[&::-webkit-slider-runnable-track]:bg-border",
            // Thumb
            "[&::-webkit-slider-thumb]:appearance-none",
            "[&::-webkit-slider-thumb]:size-4",
            "[&::-webkit-slider-thumb]:rounded-full",
            "[&::-webkit-slider-thumb]:bg-primary",
            "[&::-webkit-slider-thumb]:border-2",
            "[&::-webkit-slider-thumb]:border-background",
            "[&::-webkit-slider-thumb]:shadow-md",
            "[&::-webkit-slider-thumb]:transition-transform",
            "[&::-webkit-slider-thumb]:duration-150",
            "[&::-webkit-slider-thumb]:hover:scale-110",
            // Firefox
            "[&::-moz-range-track]:rounded-full",
            "[&::-moz-range-track]:bg-border",
            "[&::-moz-range-thumb]:appearance-none",
            "[&::-moz-range-thumb]:size-4",
            "[&::-moz-range-thumb]:rounded-full",
            "[&::-moz-range-thumb]:bg-primary",
            "[&::-moz-range-thumb]:border-2",
            "[&::-moz-range-thumb]:border-background",
            "[&::-moz-range-thumb]:shadow-md",
            // Focus
            "focus:outline-none",
            "focus-visible:[&::-webkit-slider-thumb]:ring-2",
            "focus-visible:[&::-webkit-slider-thumb]:ring-ring",
            "focus-visible:[&::-webkit-slider-thumb]:ring-offset-2",
            // Disabled
            "disabled:opacity-50",
            "disabled:cursor-not-allowed",
            sizeStyles[size],
            className
          )}
          {...props}
        />
      </div>
    );
  }
);

Slider.displayName = "Slider";
