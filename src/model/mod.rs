//! Core data structures modeling biological macromolecules.
//!
//! This module defines the foundational types for representing atoms, residues, chains,
//! structures, and topologies. These types form the backbone of `bio-forge` and are
//! consumed and mutated by I/O parsers, operations pipelines, and export routines.

pub mod atom;
pub mod chain;
pub mod grid;
pub mod residue;
pub mod structure;
pub mod template;
pub mod topology;
pub mod types;
