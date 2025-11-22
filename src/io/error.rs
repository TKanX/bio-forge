use std::path::PathBuf;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error("I/O error for file '{path}': {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("failed to parse {format} file '{path}': {details} (line {line_number})")]
    Parse {
        format: &'static str,
        path: PathBuf,
        line_number: usize,
        details: String,
    },

    #[error("unknown standard residue name '{name}' in file '{path}' could not be resolved")]
    UnknownStandardResidue { name: String, path: PathBuf },

    #[error("inconsistent data in {format} file '{path}': {details}")]
    InconsistentData {
        format: &'static str,
        path: PathBuf,
        details: String,
    },
}

impl Error {
    pub fn from_io(source: std::io::Error, path: impl Into<PathBuf>) -> Self {
        Self::Io {
            path: path.into(),
            source,
        }
    }

    pub fn parse(
        format: &'static str,
        path: PathBuf,
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

    pub fn unknown_standard_residue(name: impl Into<String>, path: PathBuf) -> Self {
        Self::UnknownStandardResidue {
            name: name.into(),
            path,
        }
    }

    pub fn inconsistent_data(
        format: &'static str,
        path: PathBuf,
        details: impl Into<String>,
    ) -> Self {
        Self::InconsistentData {
            format,
            path,
            details: details.into(),
        }
    }
}
