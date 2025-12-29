/**
 * @file Hero section
 *
 * Landing page hero with title and CTA.
 */

"use client";

import { useState } from "react";
import { useRouter } from "next/navigation";
import { motion } from "framer-motion";
import { Button } from "@/ui/primitives";
import { ChevronRightIcon, DownloadIcon, LoaderIcon } from "@/ui/icons";
import { RELEASES_URL } from "@/lib";

// ============================================================================
// Component
// ============================================================================

export function Hero() {
  const router = useRouter();
  const [isLaunching, setIsLaunching] = useState(false);

  const handleLaunch = () => {
    setIsLaunching(true);
    router.push("/app");
  };

  return (
    <section className="relative min-h-screen flex flex-col items-center justify-center px-6 overflow-hidden">
      {/* Background effects */}
      <div className="absolute inset-0 -z-10">
        <div className="absolute top-1/4 left-1/4 w-96 h-96 bg-primary/20 rounded-full blur-[120px]" />
        <div className="absolute bottom-1/4 right-1/4 w-80 h-80 bg-accent/20 rounded-full blur-[100px]" />
      </div>

      {/* Logo and title */}
      <motion.div
        initial={{ opacity: 0, y: 30 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.8, ease: "easeOut" }}
        className="text-center"
      >
        <motion.div
          initial={{ scale: 0.8, opacity: 0 }}
          animate={{ scale: 1, opacity: 1 }}
          transition={{ duration: 0.6, delay: 0.2 }}
          className="mb-8"
        >
          <h1 className="text-7xl md:text-9xl font-bold tracking-tight">
            <span className="bg-linear-to-r from-primary via-primary to-accent bg-clip-text text-transparent">
              Bio
            </span>
            <span className="text-foreground">Forge</span>
          </h1>
        </motion.div>

        <motion.p
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.6, delay: 0.4 }}
          className="text-xl md:text-2xl text-muted-foreground max-w-2xl mx-auto mb-12"
        >
          High-performance molecular structure preparation toolkit.
          <br />
          <span className="text-foreground/80">Clean. Repair. Forge.</span>
        </motion.p>

        {/* CTA buttons */}
        <motion.div
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.6, delay: 0.6 }}
          className="flex flex-col sm:flex-row gap-4 justify-center"
        >
          <Button
            size="lg"
            className="gap-2 text-base px-8"
            onClick={handleLaunch}
            disabled={isLaunching}
          >
            {isLaunching ? (
              <>
                <LoaderIcon className="size-5 animate-spin" />
                Loading...
              </>
            ) : (
              <>
                Launch Web App
                <ChevronRightIcon className="size-5" />
              </>
            )}
          </Button>
          <a href={RELEASES_URL} target="_blank" rel="noopener noreferrer">
            <Button
              variant="secondary"
              size="lg"
              className="gap-2 text-base px-8"
            >
              <DownloadIcon className="size-5" />
              Download CLI
            </Button>
          </a>
        </motion.div>
      </motion.div>

      {/* Scroll indicator */}
      <motion.div
        initial={{ opacity: 0 }}
        animate={{ opacity: 1 }}
        transition={{ delay: 1.2 }}
        className="absolute bottom-8"
      >
        <motion.div
          animate={{ y: [0, 8, 0] }}
          transition={{ duration: 2, repeat: Infinity }}
          className="w-6 h-10 rounded-full border-2 border-muted-foreground/30 flex items-start justify-center p-2"
        >
          <motion.div className="w-1.5 h-1.5 rounded-full bg-muted-foreground" />
        </motion.div>
      </motion.div>
    </section>
  );
}
