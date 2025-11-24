//! Canonical error type for all structure and template IO operations.
//!
//! This module wraps parser, serializer, and filesystem failures into a single
//! `Error` enum that higher-level operations can bubble up or convert into
//! user-facing diagnostics with uniform wording.

use std::fmt;
use std::path::PathBuf;
use thiserror::Error;

/// Errors that can occur while reading or writing biomolecular data.
///
/// The enum captures I/O failures, structured parser issues, unknown residues, and
/// integrity mismatches so callers can inspect the variant and react accordingly.
#[derive(Debug, Error)]
pub enum Error {
    /// Wrapper around operating-system level I/O failures.
    ///
    /// Includes both filesystem and stream sources, optionally carrying the file path for
    /// richer error messages.
    #[error(
        "I/O error for {path_desc}: {source}",
        path_desc = PathDisplay(path)
    )]
    Io {
        /// Path to the file involved in the failed operation, if any.
        path: Option<PathBuf>,
        /// Underlying error emitted by the standard library.
        #[source]
        source: std::io::Error,
    },

    /// Indicates that an input line could not be parsed into the expected record.
    ///
    /// Exposes the textual format, source path, failing line number, and an explanatory
    /// detail string to assist with debugging malformed files.
    #[error(
        "failed to parse {format} {path_desc}: {details} (line {line_number})",
        path_desc = PathDisplay(path)
    )]
    Parse {
        /// Name of the textual format (e.g., `"PDB"`, `"mmCIF"`).
        format: &'static str,
        /// Path to the offending file, if known.
        path: Option<PathBuf>,
        /// One-based line number where parsing failed.
        line_number: usize,
        /// Human-readable description of what went wrong.
        details: String,
    },

    /// Raised when an atom record references a residue not present in template libraries.
    #[error(
        "unknown standard residue name '{name}' in {path_desc} could not be resolved",
        path_desc = PathDisplay(path)
    )]
    UnknownStandardResidue { name: String, path: Option<PathBuf> },

    /// Reports logical inconsistencies such as mismatched atom counts or invalid records.
    #[error(
        "inconsistent data in {format} {path_desc}: {details}",
        path_desc = PathDisplay(path)
    )]
    InconsistentData {
        /// Name of the textual format being processed.
        format: &'static str,
        /// Related file path when available.
        path: Option<PathBuf>,
        /// Summary of the detected inconsistency.
        details: String,
    },
}

impl Error {
    /// Constructs an [`Error::Io`] variant from a standard I/O error.
    ///
    /// # Arguments
    ///
    /// * `source` - The original `std::io::Error` emitted by the OS or runtime.
    /// * `path` - Optional file path associated with the operation.
    ///
    /// # Returns
    ///
    /// A ready-to-use `Error` that preserves the source error for chaining.
    pub fn from_io(source: std::io::Error, path: Option<PathBuf>) -> Self {
        Self::Io { path, source }
    }

    /// Builds a [`Error::Parse`] variant with consistent messaging.
    ///
    /// # Arguments
    ///
    /// * `format` - Name of the textual format being parsed.
    /// * `path` - Optional path pointing to the input file.
    /// * `line_number` - Line where the failure occurred (1-indexed).
    /// * `details` - Additional context about the parsing problem.
    ///
    /// # Returns
    ///
    /// An `Error::Parse` carrying the supplied metadata.
    pub fn parse(
        format: &'static str,
        path: Option<PathBuf>,
        line_number: usize,
        details: impl Into<String>,
    ) -> Self {
        Self::Parse {
            format,
            path,
            line_number,
            details: details.into(),
        }
    }

    /// Generates an [`Error::UnknownStandardResidue`] for unresolved template names.
    ///
    /// # Arguments
    ///
    /// * `name` - Residue identifier that failed to look up.
    /// * `path` - Optional path to the source file.
    ///
    /// # Returns
    ///
    /// An error variant signalling the missing residue definition.
    pub fn unknown_standard_residue(name: impl Into<String>, path: Option<PathBuf>) -> Self {
        Self::UnknownStandardResidue {
            name: name.into(),
            path,
        }
    }

    /// Creates an [`Error::InconsistentData`] describing logical mismatches.
    ///
    /// # Arguments
    ///
    /// * `format` - Name of the textual format being processed.
    /// * `path` - Optional file path, if applicable.
    /// * `details` - Explanation of the inconsistency.
    ///
    /// # Returns
    ///
    /// An error variant pointing to structural or semantic discrepancies.
    pub fn inconsistent_data(
        format: &'static str,
        path: Option<PathBuf>,
        details: impl Into<String>,
    ) -> Self {
        Self::InconsistentData {
            format,
            path,
            details: details.into(),
        }
    }
}

/// Lightweight formatter for optional paths used in error messages.
///
/// When a path is present it prints `file '<path>'`; otherwise it emits `stream source` so
/// error messages remain grammatically consistent.
struct PathDisplay<'a>(&'a Option<PathBuf>);

impl<'a> fmt::Display for PathDisplay<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.0 {
            Some(p) => write!(f, "file '{}'", p.display()),
            None => write!(f, "stream source"),
        }
    }
}
