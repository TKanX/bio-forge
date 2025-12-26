//! Rule-driven structure cleanup routines that strip solvent, ions, hydrogens, and
//! user-specified residues prior to downstream modeling steps.
//!
//! The cleaners operate directly on mutable [`Structure`]
//! instances and honor fine-grained controls such as keep/remove lists so workflows can
//! standardize inputs before repair, hydrogenation, or topology building.

use crate::model::structure::Structure;
use crate::model::types::{ResidueCategory, StandardResidue};
use crate::ops::error::Error;
use std::collections::HashSet;

/// Configuration switches describing which components should be removed during cleaning.
///
/// Every flag maps to a different biological filter (water, ions, heterogens, etc.), and the
/// string sets allow explicit whitelist/blacklist overrides to preserve ligands of interest.
#[derive(Debug, Clone, Default)]
pub struct CleanConfig {
    /// Strip crystallographic waters (HOH) when `true`.
    pub remove_water: bool,
    /// Strip ionic residues (category `Ion`) when `true`.
    pub remove_ions: bool,
    /// Remove hydrogen atoms before structural refinement.
    pub remove_hydrogens: bool,
    /// Remove heterogen residues (category `Hetero`).
    pub remove_hetero: bool,
    /// Case-sensitive residue names to always remove, regardless of category.
    pub remove_residue_names: HashSet<String>,
    /// Case-sensitive residue names to always keep, overriding other rules.
    pub keep_residue_names: HashSet<String>,
}

impl CleanConfig {
    /// Convenience constructor that only removes waters.
    ///
    /// Useful for workflows that only need to strip solvent boxes before further processing.
    ///
    /// # Returns
    ///
    /// A [`CleanConfig`] with `remove_water` enabled and all other fields at their defaults.
    pub fn water_only() -> Self {
        Self {
            remove_water: true,
            ..Default::default()
        }
    }

    /// Convenience constructor that removes waters and ions simultaneously.
    ///
    /// Often used for preparing apo structures prior to topology generation.
    ///
    /// # Returns
    ///
    /// A [`CleanConfig`] enabling `remove_water` and `remove_ions`.
    pub fn water_and_ions() -> Self {
        Self {
            remove_water: true,
            remove_ions: true,
            ..Default::default()
        }
    }
}

/// Applies the cleaning rules to a mutable structure in-place.
///
/// Hydrogens can be stripped prior to solvation or protonation, and residues matching the
/// configured filters are removed with chain bookkeeping handled automatically.
///
/// # Arguments
///
/// * `structure` - Mutable structure that will be filtered.
/// * `config` - Cleaning switches describing which components to remove or keep.
///
/// # Returns
///
/// `Ok(())` when the structure is processed successfully.
///
/// # Errors
///
/// Currently never returns [`Error`] variants but reserves the signature for future
/// validation failures to stay compatible with other ops APIs.
pub fn clean_structure(structure: &mut Structure, config: &CleanConfig) -> Result<(), Error> {
    structure.par_retain_residues_mut(|_chain_id, residue| {
        if config.keep_residue_names.contains(residue.name.as_str()) {
            if config.remove_hydrogens {
                residue.strip_hydrogens();
            }
            return true;
        }

        if config.remove_residue_names.contains(residue.name.as_str()) {
            return false;
        }

        if config.remove_water && residue.standard_name == Some(StandardResidue::HOH) {
            return false;
        }

        if config.remove_ions && residue.category == ResidueCategory::Ion {
            return false;
        }

        if config.remove_hetero && residue.category == ResidueCategory::Hetero {
            return false;
        }

        if config.remove_hydrogens {
            residue.strip_hydrogens();
        }

        true
    });

    structure.prune_empty_chains();

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{
        atom::Atom,
        chain::Chain,
        residue::Residue,
        structure::Structure,
        types::{Element, Point},
    };

    fn make_structure(
        residues: Vec<(String, ResidueCategory, Option<StandardResidue>)>,
    ) -> Structure {
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");

        for (idx, (name, category, standard)) in residues.into_iter().enumerate() {
            let mut residue = Residue::new(idx as i32 + 1, None, &name, standard, category);
            residue.add_atom(Atom::new("C", Element::C, Point::origin()));
            chain.add_residue(residue);
        }

        structure.add_chain(chain);
        structure
    }

    #[test]
    fn removes_hydrogens_when_flag_enabled() {
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");
        let mut residue = Residue::new(
            1,
            None,
            "GLY",
            Some(StandardResidue::GLY),
            ResidueCategory::Standard,
        );
        residue.add_atom(Atom::new("CA", Element::C, Point::origin()));
        residue.add_atom(Atom::new("HA", Element::H, Point::new(1.0, 0.0, 0.0)));
        chain.add_residue(residue);
        structure.add_chain(chain);

        let config = CleanConfig {
            remove_hydrogens: true,
            ..Default::default()
        };

        clean_structure(&mut structure, &config).unwrap();

        let residue = structure.chain("A").unwrap().residue(1, None).unwrap();
        assert_eq!(residue.atom_count(), 1);
        assert!(residue.atom("HA").is_none());
    }

    #[test]
    fn removes_water_and_ions_when_configured() {
        let mut structure = make_structure(vec![
            (
                "HOH".to_string(),
                ResidueCategory::Standard,
                Some(StandardResidue::HOH),
            ),
            ("NA".to_string(), ResidueCategory::Ion, None),
            (
                "GLY".to_string(),
                ResidueCategory::Standard,
                Some(StandardResidue::GLY),
            ),
        ]);

        let config = CleanConfig {
            remove_water: true,
            remove_ions: true,
            ..Default::default()
        };

        clean_structure(&mut structure, &config).unwrap();

        let chain = structure.chain("A").unwrap();
        assert_eq!(chain.residue_count(), 1);
        assert_eq!(chain.residue(3, None).unwrap().name, "GLY");
    }

    #[test]
    fn removes_named_residues_and_hetero_categories() {
        let mut structure = make_structure(vec![
            ("LIG".to_string(), ResidueCategory::Hetero, None),
            ("SO4".to_string(), ResidueCategory::Hetero, None),
            (
                "ALA".to_string(),
                ResidueCategory::Standard,
                Some(StandardResidue::ALA),
            ),
        ]);

        let config = CleanConfig {
            remove_hetero: true,
            remove_residue_names: HashSet::from(["SO4".to_string()]),
            ..Default::default()
        };

        clean_structure(&mut structure, &config).unwrap();

        let chain = structure.chain("A").unwrap();
        let names: Vec<_> = chain.iter_residues().map(|res| res.name.as_str()).collect();
        assert_eq!(names, vec!["ALA"]);
    }

    #[test]
    fn keep_list_overrides_removal_rules() {
        let mut structure = make_structure(vec![(
            "HOH".to_string(),
            ResidueCategory::Standard,
            Some(StandardResidue::HOH),
        )]);

        let config = CleanConfig {
            remove_water: true,
            keep_residue_names: HashSet::from(["HOH".to_string()]),
            ..Default::default()
        };

        clean_structure(&mut structure, &config).unwrap();

        assert_eq!(structure.chain("A").unwrap().residue_count(), 1);
    }

    #[test]
    fn prunes_empty_chains_after_residue_removal() {
        let mut structure = make_structure(vec![(
            "HOH".to_string(),
            ResidueCategory::Standard,
            Some(StandardResidue::HOH),
        )]);
        let mut chain_b = Chain::new("B");
        let mut residue_b = Residue::new(
            10,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        residue_b.add_atom(Atom::new("CA", Element::C, Point::new(2.0, 0.0, 0.0)));
        chain_b.add_residue(residue_b);
        structure.add_chain(chain_b);

        let config = CleanConfig {
            remove_water: true,
            ..Default::default()
        };

        clean_structure(&mut structure, &config).unwrap();

        assert!(structure.chain("A").is_none());
        assert!(structure.chain("B").is_some());
    }
}
