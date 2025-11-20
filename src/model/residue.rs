use super::atom::Atom;
use super::types::{ResidueCategory, ResiduePosition, StandardResidue};
use std::fmt;

#[derive(Debug, Clone, PartialEq)]
pub struct Residue {
    pub id: i32,
    pub name: String,
    pub standard_name: Option<StandardResidue>,
    pub category: ResidueCategory,
    pub position: ResiduePosition,
    atoms: Vec<Atom>,
}

impl Residue {
    pub fn new(
        id: i32,
        name: &str,
        standard_name: Option<StandardResidue>,
        category: ResidueCategory,
    ) -> Self {
        Self {
            id,
            name: name.to_string(),
            standard_name,
            category,
            position: ResiduePosition::None,
            atoms: Vec::new(),
        }
    }

    pub fn is_standard(&self) -> bool {
        self.standard_name.is_some()
    }

    pub fn add_atom(&mut self, atom: Atom) {
        debug_assert!(
            self.atom(&atom.name).is_none(),
            "Attempted to add a duplicate atom name '{}' to residue '{}'",
            atom.name,
            self.name
        );
        self.atoms.push(atom);
    }

    pub fn remove_atom(&mut self, name: &str) -> Option<Atom> {
        if let Some(index) = self.atoms.iter().position(|a| a.name == name) {
            Some(self.atoms.remove(index))
        } else {
            None
        }
    }

    pub fn atom(&self, name: &str) -> Option<&Atom> {
        self.atoms.iter().find(|a| a.name == name)
    }

    pub fn atom_mut(&mut self, name: &str) -> Option<&mut Atom> {
        self.atoms.iter_mut().find(|a| a.name == name)
    }

    pub fn has_atom(&self, name: &str) -> bool {
        self.atom(name).is_some()
    }

    pub fn atoms(&self) -> &[Atom] {
        &self.atoms
    }

    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    pub fn iter_atoms(&self) -> std::slice::Iter<'_, Atom> {
        self.atoms.iter()
    }

    pub fn iter_atoms_mut(&mut self) -> std::slice::IterMut<'_, Atom> {
        self.atoms.iter_mut()
    }

    pub fn strip_hydrogens(&mut self) {
        self.atoms
            .retain(|a| a.element != crate::model::types::Element::H);
    }
}

impl fmt::Display for Residue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(std_name) = self.standard_name {
            write!(
                f,
                "Residue {{ id: {}, name: \"{}\" ({}), category: {}, atoms: {} }}",
                self.id,
                self.name,
                std_name,
                self.category,
                self.atom_count()
            )
        } else {
            write!(
                f,
                "Residue {{ id: {}, name: \"{}\", category: {}, atoms: {} }}",
                self.id,
                self.name,
                self.category,
                self.atom_count()
            )
        }
    }
}
