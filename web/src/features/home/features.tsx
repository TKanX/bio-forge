/**
 * @file Features section
 *
 * Feature highlights grid.
 */

"use client";

import { motion } from "framer-motion";
import { ZapIcon, LayersIcon, BoxIcon, FileTextIcon } from "@/ui/icons";

// ============================================================================
// Data
// ============================================================================

const features = [
  {
    icon: ZapIcon,
    title: "Blazing Fast",
    description:
      "Millisecond-level processing of million-atom systems with parallel computation via Rayon.",
  },
  {
    icon: LayersIcon,
    title: "Complete Pipeline",
    description:
      "Clean, repair, protonate, solvate, and build topologies—all in one unified workflow.",
  },
  {
    icon: BoxIcon,
    title: "WebAssembly Powered",
    description:
      "Run natively in your browser. No server uploads, no data leaves your machine.",
  },
  {
    icon: FileTextIcon,
    title: "Format Flexibility",
    description:
      "Seamless support for PDB and mmCIF formats with precise parsing diagnostics.",
  },
];

// ============================================================================
// Component
// ============================================================================

export function Features() {
  return (
    <section className="py-24 px-6">
      <div className="max-w-6xl mx-auto">
        <motion.div
          initial={{ opacity: 0, y: 20 }}
          whileInView={{ opacity: 1, y: 0 }}
          viewport={{ once: true }}
          transition={{ duration: 0.6 }}
          className="text-center mb-16"
        >
          <h2 className="text-3xl md:text-4xl font-bold mb-4">Why BioForge?</h2>
          <p className="text-muted-foreground text-lg max-w-xl mx-auto">
            Built from the ground up in Rust for maximum performance and
            reliability.
          </p>
        </motion.div>

        <div className="grid md:grid-cols-2 gap-6">
          {features.map((feature, i) => (
            <motion.div
              key={feature.title}
              initial={{ opacity: 0, y: 20 }}
              whileInView={{ opacity: 1, y: 0 }}
              viewport={{ once: true }}
              transition={{ duration: 0.6, delay: i * 0.1 }}
              className="p-6 rounded-2xl bg-card border border-border hover:border-primary/30 transition-colors"
            >
              <div className="size-12 rounded-xl bg-primary/10 flex items-center justify-center mb-4">
                <feature.icon className="size-6 text-primary" />
              </div>
              <h3 className="text-xl font-semibold mb-2">{feature.title}</h3>
              <p className="text-muted-foreground">{feature.description}</p>
            </motion.div>
          ))}
        </div>
      </div>
    </section>
  );
}
