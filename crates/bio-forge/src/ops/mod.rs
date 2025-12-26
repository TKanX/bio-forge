//! High-level operations that clean, repair, solvate, and otherwise prepare structures.
//!
//! This module groups together the public entry points for structure processing:
//! cleaning, hydrogen addition, repairs, solvation, coordinate transforms, and
//! topology reconstruction. Each submodule exposes a cohesive API and shares a
//! common error type so downstream consumers can compose workflows easily.

mod clean;
mod error;
mod hydro;
mod repair;
mod solvate;
mod topology;
mod transform;

pub use clean::{CleanConfig, clean_structure};

pub use repair::repair_structure;

pub use hydro::{HisStrategy, HydroConfig, add_hydrogens};

pub use solvate::{Anion, Cation, SolvateConfig, solvate_structure};

pub use transform::Transform;

pub use topology::TopologyBuilder;

pub use error::Error;
