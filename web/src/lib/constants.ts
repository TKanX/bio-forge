/**
 * @file Constants
 *
 * Application-wide constants and configuration values.
 */

/** Application name */
export const APP_NAME = "BioForge";

/** GitHub repository URL */
export const GITHUB_URL = "https://github.com/TKanX/bio-forge";

/** npm package URL */
export const NPM_URL = "https://www.npmjs.com/package/bio-forge";

/** Documentation URLs */
export const DOCS = {
  manual: `${GITHUB_URL}/blob/main/MANUAL.md`,
  api: `${GITHUB_URL}/blob/main/API.md`,
  rustDocs: "https://docs.rs/bio-forge",
  examples: `${GITHUB_URL}/tree/main/examples`,
} as const;

/** Maximum chains to show in collapsed view */
export const MAX_VISIBLE_CHAINS = 10;

/** Animation durations (ms) */
export const ANIMATION = {
  fast: 150,
  base: 200,
  slow: 300,
} as const;

/** File drop zone accepted extensions */
export const ACCEPTED_STRUCTURE_EXTENSIONS = [".pdb", ".ent", ".cif", ".mmcif"];
export const ACCEPTED_TEMPLATE_EXTENSIONS = [".mol2"];
export const ACCEPTED_EXTENSIONS = [
  ...ACCEPTED_STRUCTURE_EXTENSIONS,
  ...ACCEPTED_TEMPLATE_EXTENSIONS,
];

// Alias for backward compatibility
export const SUPPORTED_EXTENSIONS = ACCEPTED_EXTENSIONS;

/** MIME types for file downloads */
export const MIME_TYPES = {
  pdb: "chemical/x-pdb",
  mmcif: "chemical/x-mmcif",
  text: "text/plain",
  zip: "application/zip",
} as const;
