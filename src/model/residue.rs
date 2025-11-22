use super::atom::Atom;
use super::types::{ResidueCategory, ResiduePosition, StandardResidue};
use std::fmt;

#[derive(Debug, Clone, PartialEq)]
pub struct Residue {
    pub id: i32,
    pub insertion_code: Option<char>,
    pub name: String,
    pub standard_name: Option<StandardResidue>,
    pub category: ResidueCategory,
    pub position: ResiduePosition,
    atoms: Vec<Atom>,
}

impl Residue {
    pub fn new(
        id: i32,
        insertion_code: Option<char>,
        name: &str,
        standard_name: Option<StandardResidue>,
        category: ResidueCategory,
    ) -> Self {
        Self {
            id,
            insertion_code,
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
        let insertion_code_str = self
            .insertion_code
            .map(|c| format!(" (ic: {})", c))
            .unwrap_or_default();

        if let Some(std_name) = self.standard_name {
            write!(
                f,
                "Residue {{ id: {}{}, name: \"{}\" ({}), category: {}, atoms: {} }}",
                self.id,
                insertion_code_str,
                self.name,
                std_name,
                self.category,
                self.atom_count()
            )
        } else {
            write!(
                f,
                "Residue {{ id: {}{}, name: \"{}\", category: {}, atoms: {} }}",
                self.id,
                insertion_code_str,
                self.name,
                self.category,
                self.atom_count()
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::types::{Element, Point};

    #[test]
    fn residue_new_creates_correct_residue() {
        let residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );

        assert_eq!(residue.id, 1);
        assert_eq!(residue.insertion_code, None);
        assert_eq!(residue.name, "ALA");
        assert_eq!(residue.standard_name, Some(StandardResidue::ALA));
        assert_eq!(residue.category, ResidueCategory::Standard);
        assert_eq!(residue.position, ResiduePosition::None);
        assert!(residue.atoms.is_empty());
    }

    #[test]
    fn residue_new_with_insertion_code() {
        let residue = Residue::new(
            1,
            Some('A'),
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );

        assert_eq!(residue.id, 1);
        assert_eq!(residue.insertion_code, Some('A'));
    }

    #[test]
    fn residue_new_with_none_standard_name() {
        let residue = Residue::new(2, None, "UNK", None, ResidueCategory::Hetero);

        assert_eq!(residue.id, 2);
        assert_eq!(residue.name, "UNK");
        assert_eq!(residue.standard_name, None);
        assert_eq!(residue.category, ResidueCategory::Hetero);
    }

    #[test]
    fn residue_is_standard_returns_true_for_standard_residue() {
        let residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        assert!(residue.is_standard());
    }

    #[test]
    fn residue_is_standard_returns_false_for_non_standard_residue() {
        let residue = Residue::new(2, None, "UNK", None, ResidueCategory::Hetero);
        assert!(!residue.is_standard());
    }

    #[test]
    fn residue_add_atom_adds_atom_correctly() {
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        let atom = Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0));

        residue.add_atom(atom);

        assert_eq!(residue.atom_count(), 1);
        assert!(residue.has_atom("CA"));
        assert_eq!(residue.atom("CA").unwrap().name, "CA");
    }

    #[test]
    fn residue_remove_atom_removes_existing_atom() {
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        let atom = Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0));
        residue.add_atom(atom);

        let removed = residue.remove_atom("CA");

        assert!(removed.is_some());
        assert_eq!(removed.unwrap().name, "CA");
        assert_eq!(residue.atom_count(), 0);
        assert!(!residue.has_atom("CA"));
    }

    #[test]
    fn residue_remove_atom_returns_none_for_nonexistent_atom() {
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );

        let removed = residue.remove_atom("NONEXISTENT");

        assert!(removed.is_none());
    }

    #[test]
    fn residue_atom_returns_correct_atom() {
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        let atom = Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0));
        residue.add_atom(atom);

        let retrieved = residue.atom("CA");

        assert!(retrieved.is_some());
        assert_eq!(retrieved.unwrap().name, "CA");
    }

    #[test]
    fn residue_atom_returns_none_for_nonexistent_atom() {
        let residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );

        let retrieved = residue.atom("NONEXISTENT");

        assert!(retrieved.is_none());
    }

    #[test]
    fn residue_atom_mut_returns_correct_mutable_atom() {
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        let atom = Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0));
        residue.add_atom(atom);

        let retrieved = residue.atom_mut("CA");

        assert!(retrieved.is_some());
        assert_eq!(retrieved.unwrap().name, "CA");
    }

    #[test]
    fn residue_atom_mut_returns_none_for_nonexistent_atom() {
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );

        let retrieved = residue.atom_mut("NONEXISTENT");

        assert!(retrieved.is_none());
    }

    #[test]
    fn residue_has_atom_returns_true_for_existing_atom() {
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        let atom = Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0));
        residue.add_atom(atom);

        assert!(residue.has_atom("CA"));
    }

    #[test]
    fn residue_has_atom_returns_false_for_nonexistent_atom() {
        let residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );

        assert!(!residue.has_atom("NONEXISTENT"));
    }

    #[test]
    fn residue_atoms_returns_correct_slice() {
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        let atom1 = Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0));
        let atom2 = Atom::new("CB", Element::C, Point::new(1.0, 0.0, 0.0));
        residue.add_atom(atom1);
        residue.add_atom(atom2);

        let atoms = residue.atoms();

        assert_eq!(atoms.len(), 2);
        assert_eq!(atoms[0].name, "CA");
        assert_eq!(atoms[1].name, "CB");
    }

    #[test]
    fn residue_atom_count_returns_correct_count() {
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );

        assert_eq!(residue.atom_count(), 0);

        let atom = Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0));
        residue.add_atom(atom);

        assert_eq!(residue.atom_count(), 1);
    }

    #[test]
    fn residue_is_empty_returns_true_for_empty_residue() {
        let residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );

        assert!(residue.is_empty());
    }

    #[test]
    fn residue_is_empty_returns_false_for_non_empty_residue() {
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        let atom = Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0));
        residue.add_atom(atom);

        assert!(!residue.is_empty());
    }

    #[test]
    fn residue_iter_atoms_iterates_correctly() {
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        let atom1 = Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0));
        let atom2 = Atom::new("CB", Element::C, Point::new(1.0, 0.0, 0.0));
        residue.add_atom(atom1);
        residue.add_atom(atom2);

        let mut names = Vec::new();
        for atom in residue.iter_atoms() {
            names.push(atom.name.clone());
        }

        assert_eq!(names, vec!["CA", "CB"]);
    }

    #[test]
    fn residue_iter_atoms_mut_iterates_correctly() {
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        let atom1 = Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0));
        residue.add_atom(atom1);

        for atom in residue.iter_atoms_mut() {
            atom.translate_by(&nalgebra::Vector3::new(1.0, 0.0, 0.0));
        }

        assert!((residue.atom("CA").unwrap().pos.x - 1.0).abs() < 1e-10);
    }

    #[test]
    fn residue_strip_hydrogens_removes_hydrogen_atoms() {
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        let carbon = Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0));
        let hydrogen1 = Atom::new("HA", Element::H, Point::new(1.0, 0.0, 0.0));
        let hydrogen2 = Atom::new("HB", Element::H, Point::new(2.0, 0.0, 0.0));
        residue.add_atom(carbon);
        residue.add_atom(hydrogen1);
        residue.add_atom(hydrogen2);

        residue.strip_hydrogens();

        assert_eq!(residue.atom_count(), 1);
        assert!(residue.has_atom("CA"));
        assert!(!residue.has_atom("HA"));
        assert!(!residue.has_atom("HB"));
    }

    #[test]
    fn residue_strip_hydrogens_preserves_non_hydrogen_atoms() {
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        let carbon = Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0));
        let nitrogen = Atom::new("N", Element::N, Point::new(1.0, 0.0, 0.0));
        let oxygen = Atom::new("O", Element::O, Point::new(2.0, 0.0, 0.0));
        residue.add_atom(carbon);
        residue.add_atom(nitrogen);
        residue.add_atom(oxygen);

        residue.strip_hydrogens();

        assert_eq!(residue.atom_count(), 3);
        assert!(residue.has_atom("CA"));
        assert!(residue.has_atom("N"));
        assert!(residue.has_atom("O"));
    }

    #[test]
    fn residue_display_formats_standard_residue_correctly() {
        let residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );

        let display = format!("{}", residue);
        let expected =
            "Residue { id: 1, name: \"ALA\" (ALA), category: Standard Residue, atoms: 0 }";

        assert_eq!(display, expected);
    }

    #[test]
    fn residue_display_formats_residue_with_insertion_code_correctly() {
        let residue = Residue::new(
            1,
            Some('A'),
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );

        let display = format!("{}", residue);
        let expected =
            "Residue { id: 1 (ic: A), name: \"ALA\" (ALA), category: Standard Residue, atoms: 0 }";

        assert_eq!(display, expected);
    }

    #[test]
    fn residue_display_formats_non_standard_residue_correctly() {
        let residue = Residue::new(2, None, "UNK", None, ResidueCategory::Hetero);

        let display = format!("{}", residue);
        let expected = "Residue { id: 2, name: \"UNK\", category: Hetero Residue, atoms: 0 }";

        assert_eq!(display, expected);
    }

    #[test]
    fn residue_display_includes_atom_count() {
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        let atom1 = Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0));
        let atom2 = Atom::new("CB", Element::C, Point::new(1.0, 0.0, 0.0));
        residue.add_atom(atom1);
        residue.add_atom(atom2);

        let display = format!("{}", residue);
        let expected =
            "Residue { id: 1, name: \"ALA\" (ALA), category: Standard Residue, atoms: 2 }";

        assert_eq!(display, expected);
    }

    #[test]
    fn residue_clone_creates_identical_copy() {
        let mut residue = Residue::new(
            1,
            Some('A'),
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        let atom = Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0));
        residue.add_atom(atom);
        residue.position = ResiduePosition::Internal;

        let cloned = residue.clone();

        assert_eq!(residue, cloned);
        assert_eq!(residue.id, cloned.id);
        assert_eq!(residue.insertion_code, cloned.insertion_code);
        assert_eq!(residue.name, cloned.name);
        assert_eq!(residue.standard_name, cloned.standard_name);
        assert_eq!(residue.category, cloned.category);
        assert_eq!(residue.position, cloned.position);
        assert_eq!(residue.atoms, cloned.atoms);
    }

    #[test]
    fn residue_partial_eq_compares_correctly() {
        let mut residue1 = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        let mut residue2 = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        let atom = Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0));
        residue1.add_atom(atom.clone());
        residue2.add_atom(atom);

        let residue3 = Residue::new(
            2,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );

        assert_eq!(residue1, residue2);
        assert_ne!(residue1, residue3);
    }
}
