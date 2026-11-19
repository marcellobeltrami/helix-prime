//! Error types for the `helix-core` bioinformatics and 3D structure engine.

use std::fmt;

/// Unified error type returned by fallible operations across `helix-core`.
#[derive(Debug)]
pub enum HelixError {
    /// Failed to parse a biological sequence, alignment, annotation, or 3D structure file.
    Parse {
        format: &'static str,
        line: Option<usize>,
        message: String,
    },
    /// Invalid alignment state or parameters (e.g., empty input or out-of-bounds range).
    Alignment(String),
    /// Invalid 3D structure data (e.g., no valid ATOM/C-alpha records).
    Structure(String),
    /// Network or remote API error when querying UniProt or AlphaFold DB.
    Database {
        service: &'static str,
        message: String,
    },
    /// Underlying filesystem I/O failure.
    Io(std::io::Error),
}

impl HelixError {
    /// Constructs a parse error with format name, optional 1-based line number, and detail message.
    pub fn parse(format: &'static str, line: Option<usize>, message: impl Into<String>) -> Self {
        Self::Parse {
            format,
            line,
            message: message.into(),
        }
    }

    /// Constructs a remote database error for UniProt or AlphaFold DB requests.
    pub fn database(service: &'static str, message: impl Into<String>) -> Self {
        Self::Database {
            service,
            message: message.into(),
        }
    }
}

impl fmt::Display for HelixError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parse {
                format,
                line: Some(line_no),
                message,
            } => write!(f, "{format} parse error at line {line_no}: {message}"),
            Self::Parse {
                format,
                line: None,
                message,
            } => write!(f, "{format} parse error: {message}"),
            Self::Alignment(msg) => write!(f, "Alignment error: {msg}"),
            Self::Structure(msg) => write!(f, "3D structure error: {msg}"),
            Self::Database { service, message } => write!(f, "{service} API error: {message}"),
            Self::Io(err) => write!(f, "I/O error: {err}"),
        }
    }
}

impl std::error::Error for HelixError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(err) => Some(err),
            _ => None,
        }
    }
}

impl From<std::io::Error> for HelixError {
    fn from(err: std::io::Error) -> Self {
        Self::Io(err)
    }
}

/// Convenience result alias for `helix-core` operations.
pub type Result<T> = std::result::Result<T, HelixError>;
