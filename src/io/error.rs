use std::fmt;
use std::path::PathBuf;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error(
        "I/O error for {path_desc}: {source}",
        path_desc = PathDisplay(path)
    )]
    Io {
        path: Option<PathBuf>,
        #[source]
        source: std::io::Error,
    },

    #[error(
        "failed to parse {format} {path_desc}: {details} (line {line_number})",
        path_desc = PathDisplay(path)
    )]
    Parse {
        format: &'static str,
        path: Option<PathBuf>,
        line_number: usize,
        details: String,
    },

    #[error(
        "unknown standard residue name '{name}' in {path_desc} could not be resolved",
        path_desc = PathDisplay(path)
    )]
    UnknownStandardResidue { name: String, path: Option<PathBuf> },

    #[error(
        "inconsistent data in {format} {path_desc}: {details}",
        path_desc = PathDisplay(path)
    )]
    InconsistentData {
        format: &'static str,
        path: Option<PathBuf>,
        details: String,
    },
}

impl Error {
    pub fn from_io(source: std::io::Error, path: Option<PathBuf>) -> Self {
        Self::Io { path, source }
    }

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

    pub fn unknown_standard_residue(name: impl Into<String>, path: Option<PathBuf>) -> Self {
        Self::UnknownStandardResidue {
            name: name.into(),
            path,
        }
    }

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

struct PathDisplay<'a>(&'a Option<PathBuf>);

impl<'a> fmt::Display for PathDisplay<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.0 {
            Some(p) => write!(f, "file '{}'", p.display()),
            None => write!(f, "stream source"),
        }
    }
}
