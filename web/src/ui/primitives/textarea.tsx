/**
 * @file Textarea component
 *
 * Multi-line text input with auto-resize.
 */

"use client";

import {
  forwardRef,
  useId,
  useRef,
  useEffect,
  type ComponentPropsWithoutRef,
} from "react";
import { cn } from "@/lib";

// ============================================================================
// Types
// ============================================================================

export interface TextareaProps extends ComponentPropsWithoutRef<"textarea"> {
  /** Label text */
  label?: string;
  /** Error message */
  error?: string;
  /** Helper text */
  hint?: string;
  /** Auto-resize height based on content */
  autoResize?: boolean;
}

// ============================================================================
// Component
// ============================================================================

export const Textarea = forwardRef<HTMLTextAreaElement, TextareaProps>(
  (
    {
      label,
      error,
      hint,
      autoResize = false,
      className,
      id,
      onChange,
      ...props
    },
    ref
  ) => {
    const generatedId = useId();
    const textareaId = id || generatedId;
    const internalRef = useRef<HTMLTextAreaElement>(null);

    // Merge refs
    const textareaRef = (ref ||
      internalRef) as React.RefObject<HTMLTextAreaElement>;

    // Auto-resize logic
    useEffect(() => {
      if (!autoResize || !textareaRef.current) return;

      const textarea = textareaRef.current;
      textarea.style.height = "auto";
      textarea.style.height = `${textarea.scrollHeight}px`;
    }, [autoResize, props.value, textareaRef]);

    const handleChange = (e: React.ChangeEvent<HTMLTextAreaElement>) => {
      if (autoResize) {
        const textarea = e.target;
        textarea.style.height = "auto";
        textarea.style.height = `${textarea.scrollHeight}px`;
      }
      onChange?.(e);
    };

    const hasError = !!error;

    return (
      <div className="space-y-1.5">
        {/* Label */}
        {label && (
          <label
            htmlFor={textareaId}
            className="block text-sm font-medium text-foreground"
          >
            {label}
          </label>
        )}

        {/* Textarea */}
        <textarea
          ref={textareaRef}
          id={textareaId}
          onChange={handleChange}
          aria-invalid={hasError}
          aria-describedby={
            hasError
              ? `${textareaId}-error`
              : hint
                ? `${textareaId}-hint`
                : undefined
          }
          className={cn(
            "w-full min-h-80px px-3 py-2",
            "rounded-md border bg-background",
            "text-sm text-foreground placeholder:text-muted-foreground",
            "transition-colors duration-150",
            "resize-y",
            autoResize && "resize-none overflow-hidden",
            // Focus
            "focus:outline-none focus:ring-2 focus:ring-ring focus:ring-offset-2",
            // Error state
            hasError
              ? "border-error focus:ring-error/20"
              : "border-input hover:border-input-hover",
            // Disabled
            "disabled:opacity-50 disabled:cursor-not-allowed",
            className
          )}
          {...props}
        />

        {/* Error message */}
        {hasError && (
          <p id={`${textareaId}-error`} className="text-xs text-error">
            {error}
          </p>
        )}

        {/* Hint text */}
        {hint && !hasError && (
          <p
            id={`${textareaId}-hint`}
            className="text-xs text-muted-foreground"
          >
            {hint}
          </p>
        )}
      </div>
    );
  }
);

Textarea.displayName = "Textarea";
