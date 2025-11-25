//! # BioForge
//!
//! **BioForge** is a pure-Rust molecular preparation engine that ingests experimental macromolecular files, reconciles them with curated residue templates, and emits topology-aware structures ready for simulation or analysis. The crate favors deterministic workflows, strong typing, and clean error surfaces so pipelines remain auditable from parsing to solvation.
//!
//! ## Features
//!
//! - **High-fidelity templates** – Embedded TOML templates for proteins, nucleic acids, and solvent ensure reproducible coordinates, elements, and bonding across runs.
//! - **Ergonomic structure model** – Lightweight `Atom`, `Residue`, `Chain`, and `Structure` types backed by `nalgebra` power geometric manipulation and metadata-aware queries.
//! - **Versatile I/O** – Buffered readers and writers for PDB, mmCIF, and MOL2 share a common context that normalizes residue aliases and captures precise diagnostics.
//! - **Compositional operations** – Cleaning, repairing, hydrogenation, solvation, topology building, and transforms live under `ops`, producing interoperable results with a unified error type.
//! - **Topology reconstruction** – `Topology` and `Bond` exports allow downstream force-field assignment or connectivity validation without re-parsing raw coordinates.

mod db;
mod model;

pub mod io;
pub mod ops;
pub mod templates;

pub use model::atom::Atom;
pub use model::chain::Chain;
pub use model::residue::Residue;
pub use model::structure::Structure;
pub use model::template::Template;
pub use model::topology::{Bond, Topology};
pub use model::types::{
    BondOrder, Element, Point, ResidueCategory, ResiduePosition, StandardResidue,
};
