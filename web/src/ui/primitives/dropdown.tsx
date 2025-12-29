/**
 * @file Dropdown component
 *
 * Dropdown menu with trigger button.
 */

"use client";

import {
  useState,
  useRef,
  useEffect,
  type ReactNode,
  type ComponentPropsWithoutRef,
} from "react";
import { motion, AnimatePresence } from "framer-motion";
import { cn } from "@/lib";

// ============================================================================
// Types
// ============================================================================

export interface DropdownItem {
  label: string;
  value: string;
  icon?: ReactNode;
  disabled?: boolean;
  destructive?: boolean;
}

export interface DropdownProps {
  /** Trigger element */
  trigger: ReactNode;
  /** Menu items */
  items: DropdownItem[];
  /** Callback when item is selected */
  onSelect: (value: string) => void;
  /** Horizontal alignment of dropdown menu */
  align?: "start" | "end";
  /** Vertical position of dropdown menu */
  side?: "top" | "bottom";
  /** Additional class names for menu */
  className?: string;
}

// ============================================================================
// Component
// ============================================================================

export function Dropdown({
  trigger,
  items,
  onSelect,
  align = "end",
  side = "bottom",
  className,
}: DropdownProps) {
  const [open, setOpen] = useState(false);
  const containerRef = useRef<HTMLDivElement>(null);

  // Close on click outside
  useEffect(() => {
    if (!open) return;

    const handleClickOutside = (event: MouseEvent) => {
      if (
        containerRef.current &&
        !containerRef.current.contains(event.target as Node)
      ) {
        setOpen(false);
      }
    };

    const handleEscape = (event: KeyboardEvent) => {
      if (event.key === "Escape") {
        setOpen(false);
      }
    };

    document.addEventListener("mousedown", handleClickOutside);
    document.addEventListener("keydown", handleEscape);

    return () => {
      document.removeEventListener("mousedown", handleClickOutside);
      document.removeEventListener("keydown", handleEscape);
    };
  }, [open]);

  const handleSelect = (item: DropdownItem) => {
    if (item.disabled) return;
    onSelect(item.value);
    setOpen(false);
  };

  return (
    <div ref={containerRef} className="relative inline-block">
      {/* Trigger */}
      <div onClick={() => setOpen(!open)} className="cursor-pointer">
        {trigger}
      </div>

      {/* Menu */}
      <AnimatePresence>
        {open && (
          <motion.div
            initial={{ opacity: 0, y: side === "top" ? 8 : -8, scale: 0.96 }}
            animate={{ opacity: 1, y: 0, scale: 1 }}
            exit={{ opacity: 0, y: side === "top" ? 8 : -8, scale: 0.96 }}
            transition={{ duration: 0.15 }}
            className={cn(
              "absolute z-100 min-w-40",
              "rounded-lg border border-border bg-card shadow-lg",
              "py-1",
              align === "start" ? "left-0" : "right-0",
              side === "top" ? "bottom-full mb-1" : "top-full mt-1",
              className
            )}
          >
            {items.map((item) => (
              <button
                key={item.value}
                type="button"
                onClick={() => handleSelect(item)}
                disabled={item.disabled}
                className={cn(
                  "w-full flex items-center gap-2 px-3 py-2",
                  "text-sm text-left",
                  "transition-colors duration-150",
                  item.disabled
                    ? "opacity-50 cursor-not-allowed"
                    : item.destructive
                      ? "text-error hover:bg-error/10"
                      : "text-foreground hover:bg-card-hover"
                )}
              >
                {item.icon && (
                  <span className="size-4 shrink-0">{item.icon}</span>
                )}
                <span>{item.label}</span>
              </button>
            ))}
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
}

// ============================================================================
// IconButton (commonly used as dropdown trigger)
// ============================================================================

export interface IconButtonProps extends ComponentPropsWithoutRef<"button"> {
  /** Size variant */
  size?: "sm" | "md" | "lg";
  /** Visual variant */
  variant?: "ghost" | "outline";
}

const sizeStyles = {
  sm: "size-7 text-sm",
  md: "size-8 text-base",
  lg: "size-9 text-lg",
} as const;

const variantStyles = {
  ghost: "hover:bg-muted",
  outline: "border border-border hover:bg-muted",
} as const;

export function IconButton({
  size = "md",
  variant = "ghost",
  className,
  children,
  ...props
}: IconButtonProps) {
  return (
    <button
      type="button"
      className={cn(
        "inline-flex items-center justify-center rounded-md",
        "transition-colors duration-150",
        "disabled:opacity-50 disabled:pointer-events-none",
        sizeStyles[size],
        variantStyles[variant],
        className
      )}
      {...props}
    >
      {children}
    </button>
  );
}
