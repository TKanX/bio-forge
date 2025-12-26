//! IO front-end exposing structure parsing and export utilities for common biomolecular formats.
//!
//! The module re-exports format-specific readers and writers so applications can import PDB,
//! mmCIF, or MOL2 data into `bio-forge` structures, enrich them via the operations pipeline,
//! and export updated coordinates or topologies without touching lower-level submodules.

mod context;
mod error;
mod mmcif;
mod mol2;
mod pdb;

pub use pdb::reader::read as read_pdb_structure;
pub use pdb::writer::{
    write_structure as write_pdb_structure, write_topology as write_pdb_topology,
};

pub use mmcif::reader::read as read_mmcif_structure;
pub use mmcif::writer::{
    write_structure as write_mmcif_structure, write_topology as write_mmcif_topology,
};

pub use mol2::reader::read as read_mol2_template;

pub use context::IoContext;

pub use error::Error;
