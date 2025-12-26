//! Fundamental atom representation comprising name, chemical element, and Cartesian position.
//!
//! This module defines the smallest structural unit used throughout `bio-forge`. Atoms are
//! instantiated by IO readers, manipulated by topology operations, and rendered back into
//! biomolecular formats. Distance helpers and translation utilities keep vector math inside
//! the type, ensuring consistent use of the chosen coordinate system.

use super::types::{Element, Point};
use smol_str::SmolStr;
use std::fmt;

/// Labeled atom with immutable element identity and mutable position.
///
/// The struct is shared across residue, chain, and structure builders. Keeping the element
/// metadata close to the coordinate allows downstream algorithms (e.g., heavy-atom filters
/// or hydrogen placement) to reason locally without traversing additional tables.
#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    /// Atom name as it appears in crystallographic or modeling files (e.g., `CA`).
    pub name: SmolStr,
    /// Chemical element derived from the periodic table definitions.
    pub element: Element,
    /// Cartesian coordinates measured in ångströms.
    pub pos: Point,
}

impl Atom {
    /// Creates a new atom from a name, element, and position.
    ///
    /// Caller controls ownership of the label string while the element enforces chemical
    /// consistency. The position is copied as-is; no normalization is performed.
    ///
    /// # Arguments
    ///
    /// * `name` - Atom label such as `"CA"` or `"OXT"`.
    /// * `element` - `Element` variant describing the chemical identity.
    /// * `pos` - `Point` describing the Cartesian coordinates in ångströms.
    ///
    /// # Returns
    ///
    /// A fully initialized `Atom` instance.
    pub fn new(name: &str, element: Element, pos: Point) -> Self {
        Self {
            name: SmolStr::new(name),
            element,
            pos,
        }
    }

    /// Computes the squared Euclidean distance to another atom.
    ///
    /// Prefer this when comparing relative distances or feeding cutoffs, as it avoids the
    /// costly square-root step while remaining in ångström squared units.
    ///
    /// # Arguments
    ///
    /// * `other` - Reference atom to measure against.
    ///
    /// # Returns
    ///
    /// The squared distance as `f64`.
    pub fn distance_squared(&self, other: &Atom) -> f64 {
        nalgebra::distance_squared(&self.pos, &other.pos)
    }

    /// Computes the Euclidean distance to another atom.
    ///
    /// This is the fully realized length in ångströms and is suitable for reporting or
    /// geometry calculations that require actual bond lengths.
    ///
    /// # Arguments
    ///
    /// * `other` - Reference atom to measure against.
    ///
    /// # Returns
    ///
    /// The distance in ångströms as `f64`.
    pub fn distance(&self, other: &Atom) -> f64 {
        nalgebra::distance(&self.pos, &other.pos)
    }

    /// Translates the atom by an arbitrary vector.
    ///
    /// The operation mutates the underlying position and is commonly used during rigid-body
    /// transforms or when applying simulation displacements.
    ///
    /// # Arguments
    ///
    /// * `vector` - Translation expressed as a `nalgebra::Vector3<f64>` in ångströms.
    pub fn translate_by(&mut self, vector: &nalgebra::Vector3<f64>) {
        self.pos += vector;
    }
}

impl fmt::Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Atom {{ name: \"{}\", element: {}, pos: [{:.3}, {:.3}, {:.3}] }}",
            self.name, self.element, self.pos.x, self.pos.y, self.pos.z
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn atom_new_creates_correct_atom() {
        let pos = Point::new(1.0, 2.0, 3.0);
        let atom = Atom::new("C1", Element::C, pos);

        assert_eq!(atom.name, "C1");
        assert_eq!(atom.element, Element::C);
        assert_eq!(atom.pos, pos);
    }

    #[test]
    fn atom_distance_squared_calculates_correctly() {
        let atom1 = Atom::new("A", Element::H, Point::new(0.0, 0.0, 0.0));
        let atom2 = Atom::new("B", Element::H, Point::new(3.0, 4.0, 0.0));

        let dist_sq = atom1.distance_squared(&atom2);
        assert!((dist_sq - 25.0).abs() < 1e-10);
    }

    #[test]
    fn atom_distance_calculates_correctly() {
        let atom1 = Atom::new("A", Element::H, Point::new(0.0, 0.0, 0.0));
        let atom2 = Atom::new("B", Element::H, Point::new(3.0, 4.0, 0.0));

        let dist = atom1.distance(&atom2);
        assert!((dist - 5.0).abs() < 1e-10);
    }

    #[test]
    fn atom_distance_squared_zero_for_same_position() {
        let pos = Point::new(1.5, -2.3, 4.7);
        let atom1 = Atom::new("A", Element::O, pos);
        let atom2 = Atom::new("B", Element::O, pos);

        let dist_sq = atom1.distance_squared(&atom2);
        assert!((dist_sq - 0.0).abs() < 1e-10);
    }

    #[test]
    fn atom_distance_zero_for_same_position() {
        let pos = Point::new(1.5, -2.3, 4.7);
        let atom1 = Atom::new("A", Element::O, pos);
        let atom2 = Atom::new("B", Element::O, pos);

        let dist = atom1.distance(&atom2);
        assert!((dist - 0.0).abs() < 1e-10);
    }

    #[test]
    fn atom_translate_by_updates_position_correctly() {
        let mut atom = Atom::new("Test", Element::N, Point::new(1.0, 2.0, 3.0));
        let vector = nalgebra::Vector3::new(0.5, -1.0, 2.5);

        atom.translate_by(&vector);

        assert!((atom.pos.x - 1.5).abs() < 1e-10);
        assert!((atom.pos.y - 1.0).abs() < 1e-10);
        assert!((atom.pos.z - 5.5).abs() < 1e-10);
    }

    #[test]
    fn atom_translate_by_with_zero_vector_no_change() {
        let mut atom = Atom::new("Test", Element::N, Point::new(1.0, 2.0, 3.0));
        let original_pos = atom.pos;
        let vector = nalgebra::Vector3::new(0.0, 0.0, 0.0);

        atom.translate_by(&vector);

        assert_eq!(atom.pos, original_pos);
    }

    #[test]
    fn atom_display_formats_correctly() {
        let atom = Atom::new("CA", Element::C, Point::new(1.234, -5.678, 9.012));

        let display = format!("{}", atom);
        let expected = "Atom { name: \"CA\", element: C, pos: [1.234, -5.678, 9.012] }";

        assert_eq!(display, expected);
    }

    #[test]
    fn atom_display_with_unknown_element() {
        let atom = Atom::new("UNK", Element::Unknown, Point::new(0.0, 0.0, 0.0));

        let display = format!("{}", atom);
        let expected = "Atom { name: \"UNK\", element: Unknown, pos: [0.000, 0.000, 0.000] }";

        assert_eq!(display, expected);
    }

    #[test]
    fn atom_clone_creates_identical_copy() {
        let atom = Atom::new("CloneTest", Element::Fe, Point::new(7.89, -1.23, 4.56));
        let cloned = atom.clone();

        assert_eq!(atom, cloned);
        assert_eq!(atom.name, cloned.name);
        assert_eq!(atom.element, cloned.element);
        assert_eq!(atom.pos, cloned.pos);
    }

    #[test]
    fn atom_partial_eq_compares_correctly() {
        let atom1 = Atom::new("Test", Element::O, Point::new(1.0, 2.0, 3.0));
        let atom2 = Atom::new("Test", Element::O, Point::new(1.0, 2.0, 3.0));
        let atom3 = Atom::new("Different", Element::O, Point::new(1.0, 2.0, 3.0));
        let atom4 = Atom::new("Test", Element::N, Point::new(1.0, 2.0, 3.0));
        let atom5 = Atom::new("Test", Element::O, Point::new(1.1, 2.0, 3.0));

        assert_eq!(atom1, atom2);
        assert_ne!(atom1, atom3);
        assert_ne!(atom1, atom4);
        assert_ne!(atom1, atom5);
    }

    #[test]
    fn atom_distance_with_negative_coordinates() {
        let atom1 = Atom::new("A", Element::H, Point::new(-1.0, -2.0, -3.0));
        let atom2 = Atom::new("B", Element::H, Point::new(1.0, 2.0, 3.0));

        let dist_sq = atom1.distance_squared(&atom2);
        assert!((dist_sq - 56.0).abs() < 1e-10);

        let dist = atom1.distance(&atom2);
        assert!((dist - (56.0_f64).sqrt()).abs() < 1e-10);
    }

    #[test]
    fn atom_translate_by_multiple_times_accumulates() {
        let mut atom = Atom::new("Test", Element::C, Point::new(0.0, 0.0, 0.0));

        atom.translate_by(&nalgebra::Vector3::new(1.0, 0.0, 0.0));
        atom.translate_by(&nalgebra::Vector3::new(0.0, 2.0, 0.0));
        atom.translate_by(&nalgebra::Vector3::new(0.0, 0.0, 3.0));

        assert!((atom.pos.x - 1.0).abs() < 1e-10);
        assert!((atom.pos.y - 2.0).abs() < 1e-10);
        assert!((atom.pos.z - 3.0).abs() < 1e-10);
    }
}
