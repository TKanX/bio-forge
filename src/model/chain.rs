use super::residue::Residue;
use std::fmt;

#[derive(Debug, Clone, PartialEq)]
pub struct Chain {
    pub id: String,
    residues: Vec<Residue>,
}

impl Chain {
    pub fn new(id: &str) -> Self {
        Self {
            id: id.to_string(),
            residues: Vec::new(),
        }
    }

    pub fn add_residue(&mut self, residue: Residue) {
        debug_assert!(
            self.residue(residue.id).is_none(),
            "Attempted to add a duplicate residue ID '{}' to chain '{}'",
            residue.id,
            self.id
        );
        self.residues.push(residue);
    }

    pub fn residue(&self, id: i32) -> Option<&Residue> {
        self.residues.iter().find(|r| r.id == id)
    }

    pub fn residue_mut(&mut self, id: i32) -> Option<&mut Residue> {
        self.residues.iter_mut().find(|r| r.id == id)
    }

    pub fn residues(&self) -> &[Residue] {
        &self.residues
    }

    pub fn residue_count(&self) -> usize {
        self.residues.len()
    }

    pub fn is_empty(&self) -> bool {
        self.residues.is_empty()
    }

    pub fn iter_residues(&self) -> std::slice::Iter<'_, Residue> {
        self.residues.iter()
    }

    pub fn iter_residues_mut(&mut self) -> std::slice::IterMut<'_, Residue> {
        self.residues.iter_mut()
    }

    pub fn iter_atoms(&self) -> impl Iterator<Item = &super::atom::Atom> {
        self.residues.iter().flat_map(|r| r.iter_atoms())
    }

    pub fn iter_atoms_mut(&mut self) -> impl Iterator<Item = &mut super::atom::Atom> {
        self.residues.iter_mut().flat_map(|r| r.iter_atoms_mut())
    }
}

impl fmt::Display for Chain {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Chain {{ id: \"{}\", residues: {} }}",
            self.id,
            self.residue_count()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::atom::Atom;
    use crate::model::types::{Element, Point, ResidueCategory, StandardResidue};

    #[test]
    fn chain_new_creates_correct_chain() {
        let chain = Chain::new("A");

        assert_eq!(chain.id, "A");
        assert!(chain.residues.is_empty());
    }

    #[test]
    fn chain_add_residue_adds_residue_correctly() {
        let mut chain = Chain::new("A");
        let residue = Residue::new(
            1,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );

        chain.add_residue(residue);

        assert_eq!(chain.residue_count(), 1);
        assert!(chain.residue(1).is_some());
        assert_eq!(chain.residue(1).unwrap().name, "ALA");
    }

    #[test]
    fn chain_residue_returns_correct_residue() {
        let mut chain = Chain::new("A");
        let residue = Residue::new(
            1,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        chain.add_residue(residue);

        let retrieved = chain.residue(1);

        assert!(retrieved.is_some());
        assert_eq!(retrieved.unwrap().id, 1);
        assert_eq!(retrieved.unwrap().name, "ALA");
    }

    #[test]
    fn chain_residue_returns_none_for_nonexistent_residue() {
        let chain = Chain::new("A");

        let retrieved = chain.residue(999);

        assert!(retrieved.is_none());
    }

    #[test]
    fn chain_residue_mut_returns_correct_mutable_residue() {
        let mut chain = Chain::new("A");
        let residue = Residue::new(
            1,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        chain.add_residue(residue);

        let retrieved = chain.residue_mut(1);

        assert!(retrieved.is_some());
        assert_eq!(retrieved.unwrap().id, 1);
    }

    #[test]
    fn chain_residue_mut_returns_none_for_nonexistent_residue() {
        let mut chain = Chain::new("A");

        let retrieved = chain.residue_mut(999);

        assert!(retrieved.is_none());
    }

    #[test]
    fn chain_residues_returns_correct_slice() {
        let mut chain = Chain::new("A");
        let residue1 = Residue::new(
            1,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        let residue2 = Residue::new(
            2,
            "GLY",
            Some(StandardResidue::GLY),
            ResidueCategory::Standard,
        );
        chain.add_residue(residue1);
        chain.add_residue(residue2);

        let residues = chain.residues();

        assert_eq!(residues.len(), 2);
        assert_eq!(residues[0].id, 1);
        assert_eq!(residues[1].id, 2);
    }

    #[test]
    fn chain_residue_count_returns_correct_count() {
        let mut chain = Chain::new("A");

        assert_eq!(chain.residue_count(), 0);

        let residue = Residue::new(
            1,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        chain.add_residue(residue);

        assert_eq!(chain.residue_count(), 1);
    }

    #[test]
    fn chain_is_empty_returns_true_for_empty_chain() {
        let chain = Chain::new("A");

        assert!(chain.is_empty());
    }

    #[test]
    fn chain_is_empty_returns_false_for_non_empty_chain() {
        let mut chain = Chain::new("A");
        let residue = Residue::new(
            1,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        chain.add_residue(residue);

        assert!(!chain.is_empty());
    }

    #[test]
    fn chain_iter_residues_iterates_correctly() {
        let mut chain = Chain::new("A");
        let residue1 = Residue::new(
            1,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        let residue2 = Residue::new(
            2,
            "GLY",
            Some(StandardResidue::GLY),
            ResidueCategory::Standard,
        );
        chain.add_residue(residue1);
        chain.add_residue(residue2);

        let mut ids = Vec::new();
        for residue in chain.iter_residues() {
            ids.push(residue.id);
        }

        assert_eq!(ids, vec![1, 2]);
    }

    #[test]
    fn chain_iter_residues_mut_iterates_correctly() {
        let mut chain = Chain::new("A");
        let residue = Residue::new(
            1,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        chain.add_residue(residue);

        for residue in chain.iter_residues_mut() {
            residue.position = crate::model::types::ResiduePosition::Internal;
        }

        assert_eq!(
            chain.residue(1).unwrap().position,
            crate::model::types::ResiduePosition::Internal
        );
    }

    #[test]
    fn chain_iter_atoms_iterates_over_all_atoms() {
        let mut chain = Chain::new("A");
        let mut residue1 = Residue::new(
            1,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        let mut residue2 = Residue::new(
            2,
            "GLY",
            Some(StandardResidue::GLY),
            ResidueCategory::Standard,
        );

        let atom1 = Atom::new("CA1", Element::C, Point::new(0.0, 0.0, 0.0));
        let atom2 = Atom::new("CB1", Element::C, Point::new(1.0, 0.0, 0.0));
        let atom3 = Atom::new("CA2", Element::C, Point::new(2.0, 0.0, 0.0));

        residue1.add_atom(atom1);
        residue1.add_atom(atom2);
        residue2.add_atom(atom3);

        chain.add_residue(residue1);
        chain.add_residue(residue2);

        let mut atom_names = Vec::new();
        for atom in chain.iter_atoms() {
            atom_names.push(atom.name.clone());
        }

        assert_eq!(atom_names, vec!["CA1", "CB1", "CA2"]);
    }

    #[test]
    fn chain_iter_atoms_mut_iterates_over_all_atoms_mutably() {
        let mut chain = Chain::new("A");
        let mut residue = Residue::new(
            1,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        let atom = Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0));
        residue.add_atom(atom);
        chain.add_residue(residue);

        for atom in chain.iter_atoms_mut() {
            atom.translate_by(&nalgebra::Vector3::new(1.0, 0.0, 0.0));
        }

        let translated_atom = chain.residue(1).unwrap().atom("CA").unwrap();
        assert!((translated_atom.pos.x - 1.0).abs() < 1e-10);
    }

    #[test]
    fn chain_iter_atoms_returns_empty_iterator_for_empty_chain() {
        let chain = Chain::new("A");

        let count = chain.iter_atoms().count();

        assert_eq!(count, 0);
    }

    #[test]
    fn chain_iter_atoms_mut_returns_empty_iterator_for_empty_chain() {
        let mut chain = Chain::new("A");

        let count = chain.iter_atoms_mut().count();

        assert_eq!(count, 0);
    }

    #[test]
    fn chain_display_formats_correctly() {
        let mut chain = Chain::new("A");
        let residue = Residue::new(
            1,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        chain.add_residue(residue);

        let display = format!("{}", chain);
        let expected = "Chain { id: \"A\", residues: 1 }";

        assert_eq!(display, expected);
    }

    #[test]
    fn chain_display_formats_empty_chain_correctly() {
        let chain = Chain::new("B");

        let display = format!("{}", chain);
        let expected = "Chain { id: \"B\", residues: 0 }";

        assert_eq!(display, expected);
    }

    #[test]
    fn chain_clone_creates_identical_copy() {
        let mut chain = Chain::new("A");
        let residue = Residue::new(
            1,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        chain.add_residue(residue);

        let cloned = chain.clone();

        assert_eq!(chain, cloned);
        assert_eq!(chain.id, cloned.id);
        assert_eq!(chain.residues, cloned.residues);
    }

    #[test]
    fn chain_partial_eq_compares_correctly() {
        let mut chain1 = Chain::new("A");
        let mut chain2 = Chain::new("A");
        let residue = Residue::new(
            1,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        chain1.add_residue(residue.clone());
        chain2.add_residue(residue);

        let chain3 = Chain::new("B");

        assert_eq!(chain1, chain2);
        assert_ne!(chain1, chain3);
    }

    #[test]
    fn chain_with_multiple_residues_and_atoms() {
        let mut chain = Chain::new("A");

        for i in 1..=3 {
            let mut residue =
                Residue::new(i, &format!("RES{}", i), None, ResidueCategory::Standard);
            let atom = Atom::new(
                &format!("ATOM{}", i),
                Element::C,
                Point::new(i as f64, 0.0, 0.0),
            );
            residue.add_atom(atom);
            chain.add_residue(residue);
        }

        assert_eq!(chain.residue_count(), 3);
        assert_eq!(chain.iter_atoms().count(), 3);

        let residue = chain.residue(2).unwrap();
        assert_eq!(residue.name, "RES2");
        assert_eq!(residue.atom("ATOM2").unwrap().name, "ATOM2");
    }
}
