//! Graph-based description of bonded connectivity for `Structure` instances.
//!
//! The topology module stores canonicalized atom-to-atom bonds, exposes iterators for
//! residue-level analysis, and provides helper utilities used by operations such as repair,
//! hydrogen completion, and solvation to reason about neighboring atoms.

use super::structure::Structure;
use super::types::BondOrder;
use std::fmt;

/// Undirected bond connecting two atoms within a structure.
///
/// Bonds store canonical atom indices (ascending order) so equality, hashing, and sorting
/// remain stable regardless of the order in which the connection was created.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Bond {
    /// Index of the first atom (always the lesser index after canonicalization).
    pub a1_idx: usize,
    /// Index of the second atom (greater-or-equal to `a1_idx`).
    pub a2_idx: usize,
    /// Chemical multiplicity assigned to the bond.
    pub order: BondOrder,
}

impl Bond {
    /// Creates a new bond while canonicalizing the endpoint ordering.
    ///
    /// The smaller atom index is stored in `a1_idx` to keep hashing and equality symmetric,
    /// ensuring duplicated bonds collapse in sets and maps.
    ///
    /// # Arguments
    ///
    /// * `idx1` - Index of one bonded atom within the owning `Structure`.
    /// * `idx2` - Index of the partner atom.
    /// * `order` - Chemical bond order describing multiplicity or aromaticity.
    ///
    /// # Returns
    ///
    /// A `Bond` whose indices are sorted so `a1_idx <= a2_idx`.
    pub fn new(idx1: usize, idx2: usize, order: BondOrder) -> Self {
        if idx1 <= idx2 {
            Self {
                a1_idx: idx1,
                a2_idx: idx2,
                order,
            }
        } else {
            Self {
                a1_idx: idx2,
                a2_idx: idx1,
                order,
            }
        }
    }
}

/// Bond graph overlay for a [`Structure`].
///
/// A `Topology` pairs structural coordinates with explicit bonds, enabling neighbor queries,
/// validation routines, and format writers that require connectivity information.
#[derive(Debug, Clone)]
pub struct Topology {
    structure: Structure,
    bonds: Vec<Bond>,
}

impl Topology {
    /// Builds a topology from a structure and its associated bonds.
    ///
    /// Indices are assumed to reference atoms within the provided structure. A debug assert
    /// validates this assumption in development builds to catch mismatched templates.
    ///
    /// # Arguments
    ///
    /// * `structure` - Fully instantiated structure containing all atoms.
    /// * `bonds` - Canonical bond list describing connectivity.
    ///
    /// # Returns
    ///
    /// A `Topology` ready for neighbor queries and downstream processing.
    pub fn new(structure: Structure, bonds: Vec<Bond>) -> Self {
        debug_assert!(
            bonds.iter().all(|b| b.a2_idx < structure.atom_count()),
            "Bond index out of bounds"
        );
        Self { structure, bonds }
    }

    /// Exposes the underlying structure.
    ///
    /// # Returns
    ///
    /// Immutable reference to the wrapped [`Structure`].
    pub fn structure(&self) -> &Structure {
        &self.structure
    }

    /// Returns all bonds present in the topology.
    ///
    /// # Returns
    ///
    /// Slice containing every [`Bond`], preserving insertion order.
    pub fn bonds(&self) -> &[Bond] {
        &self.bonds
    }

    /// Counts the number of stored bonds.
    ///
    /// # Returns
    ///
    /// Total number of bonds in the topology.
    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    /// Counts the atoms tracked by the underlying structure.
    ///
    /// # Returns
    ///
    /// Number of atoms derived from the wrapped structure.
    pub fn atom_count(&self) -> usize {
        self.structure.atom_count()
    }

    /// Iterates over all bonds incident to a specific atom.
    ///
    /// # Arguments
    ///
    /// * `atom_idx` - Index of the atom whose incident bonds should be returned.
    ///
    /// # Returns
    ///
    /// Iterator yielding references to [`Bond`] instances connected to `atom_idx`.
    pub fn bonds_of(&self, atom_idx: usize) -> impl Iterator<Item = &Bond> {
        self.bonds
            .iter()
            .filter(move |b| b.a1_idx == atom_idx || b.a2_idx == atom_idx)
    }

    /// Enumerates the neighboring atom indices for the provided atom.
    ///
    /// # Arguments
    ///
    /// * `atom_idx` - Index of the atom whose neighbors will be traversed.
    ///
    /// # Returns
    ///
    /// Iterator producing the indices of atoms bonded to `atom_idx`.
    pub fn neighbors_of(&self, atom_idx: usize) -> impl Iterator<Item = usize> + '_ {
        self.bonds_of(atom_idx).map(move |b| {
            if b.a1_idx == atom_idx {
                b.a2_idx
            } else {
                b.a1_idx
            }
        })
    }
}

impl fmt::Display for Topology {
    /// Formats the topology by reporting the atom and bond counts.
    ///
    /// This user-friendly summary is leveraged in logs and debugging output.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Topology {{ atoms: {}, bonds: {} }}",
            self.atom_count(),
            self.bond_count()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::atom::Atom;
    use crate::model::chain::Chain;
    use crate::model::residue::Residue;
    use crate::model::types::{Element, Point, ResidueCategory, StandardResidue};

    #[test]
    fn bond_new_creates_bond_with_canonical_ordering() {
        let bond = Bond::new(5, 2, BondOrder::Single);

        assert_eq!(bond.a1_idx, 2);
        assert_eq!(bond.a2_idx, 5);
        assert_eq!(bond.order, BondOrder::Single);
    }

    #[test]
    fn bond_new_preserves_order_when_already_canonical() {
        let bond = Bond::new(2, 5, BondOrder::Double);

        assert_eq!(bond.a1_idx, 2);
        assert_eq!(bond.a2_idx, 5);
        assert_eq!(bond.order, BondOrder::Double);
    }

    #[test]
    fn bond_new_handles_same_indices() {
        let bond = Bond::new(3, 3, BondOrder::Triple);

        assert_eq!(bond.a1_idx, 3);
        assert_eq!(bond.a2_idx, 3);
        assert_eq!(bond.order, BondOrder::Triple);
    }

    #[test]
    fn topology_new_creates_topology_correctly() {
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        residue.add_atom(Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0)));
        residue.add_atom(Atom::new("CB", Element::C, Point::new(1.0, 0.0, 0.0)));
        chain.add_residue(residue);
        structure.add_chain(chain);

        let bonds = vec![Bond::new(0, 1, BondOrder::Single)];
        let topology = Topology::new(structure, bonds);

        assert_eq!(topology.atom_count(), 2);
        assert_eq!(topology.bond_count(), 1);
    }

    #[test]
    fn topology_structure_returns_correct_reference() {
        let structure = Structure::new();
        let topology = Topology::new(structure, Vec::new());

        let retrieved = topology.structure();

        assert_eq!(retrieved.chain_count(), 0);
    }

    #[test]
    fn topology_bonds_returns_correct_slice() {
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        residue.add_atom(Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0)));
        residue.add_atom(Atom::new("CB", Element::C, Point::new(1.0, 0.0, 0.0)));
        residue.add_atom(Atom::new("N", Element::N, Point::new(2.0, 0.0, 0.0)));
        chain.add_residue(residue);
        structure.add_chain(chain);

        let bonds = vec![
            Bond::new(0, 1, BondOrder::Single),
            Bond::new(1, 2, BondOrder::Double),
        ];
        let topology = Topology::new(structure, bonds);

        let retrieved = topology.bonds();

        assert_eq!(retrieved.len(), 2);
        assert_eq!(retrieved[0].order, BondOrder::Single);
        assert_eq!(retrieved[1].order, BondOrder::Double);
    }

    #[test]
    fn topology_bond_count_returns_correct_count() {
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        residue.add_atom(Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0)));
        residue.add_atom(Atom::new("CB", Element::C, Point::new(1.0, 0.0, 0.0)));
        chain.add_residue(residue);
        structure.add_chain(chain);

        let bonds = vec![Bond::new(0, 1, BondOrder::Single)];
        let topology = Topology::new(structure, bonds);

        assert_eq!(topology.bond_count(), 1);
    }

    #[test]
    fn topology_bond_count_returns_zero_for_empty_topology() {
        let structure = Structure::new();
        let topology = Topology::new(structure, Vec::new());

        assert_eq!(topology.bond_count(), 0);
    }

    #[test]
    fn topology_atom_count_returns_correct_count() {
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        residue.add_atom(Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0)));
        residue.add_atom(Atom::new("CB", Element::C, Point::new(1.0, 0.0, 0.0)));
        chain.add_residue(residue);
        structure.add_chain(chain);

        let topology = Topology::new(structure, Vec::new());

        assert_eq!(topology.atom_count(), 2);
    }

    #[test]
    fn topology_atom_count_returns_zero_for_empty_structure() {
        let structure = Structure::new();
        let topology = Topology::new(structure, Vec::new());

        assert_eq!(topology.atom_count(), 0);
    }

    #[test]
    fn topology_bonds_of_returns_correct_bonds() {
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        residue.add_atom(Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0)));
        residue.add_atom(Atom::new("CB", Element::C, Point::new(1.0, 0.0, 0.0)));
        residue.add_atom(Atom::new("N", Element::N, Point::new(2.0, 0.0, 0.0)));
        chain.add_residue(residue);
        structure.add_chain(chain);

        let bonds = vec![
            Bond::new(0, 1, BondOrder::Single),
            Bond::new(1, 2, BondOrder::Double),
        ];
        let topology = Topology::new(structure, bonds);

        let bonds_of_1: Vec<_> = topology.bonds_of(1).collect();

        assert_eq!(bonds_of_1.len(), 2);
        assert!(bonds_of_1.iter().any(|b| b.a1_idx == 0 && b.a2_idx == 1));
        assert!(bonds_of_1.iter().any(|b| b.a1_idx == 1 && b.a2_idx == 2));
    }

    #[test]
    fn topology_bonds_of_returns_empty_for_atom_with_no_bonds() {
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        residue.add_atom(Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0)));
        residue.add_atom(Atom::new("CB", Element::C, Point::new(1.0, 0.0, 0.0)));
        chain.add_residue(residue);
        structure.add_chain(chain);

        let bonds = vec![Bond::new(0, 1, BondOrder::Single)];
        let topology = Topology::new(structure, bonds);

        let bonds_of_5: Vec<_> = topology.bonds_of(5).collect();

        assert!(bonds_of_5.is_empty());
    }

    #[test]
    fn topology_neighbors_of_returns_correct_neighbors() {
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        residue.add_atom(Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0)));
        residue.add_atom(Atom::new("CB", Element::C, Point::new(1.0, 0.0, 0.0)));
        residue.add_atom(Atom::new("N", Element::N, Point::new(2.0, 0.0, 0.0)));
        residue.add_atom(Atom::new("O", Element::O, Point::new(3.0, 0.0, 0.0)));
        chain.add_residue(residue);
        structure.add_chain(chain);

        let bonds = vec![
            Bond::new(0, 1, BondOrder::Single),
            Bond::new(1, 2, BondOrder::Double),
            Bond::new(1, 3, BondOrder::Triple),
        ];
        let topology = Topology::new(structure, bonds);

        let neighbors: Vec<_> = topology.neighbors_of(1).collect();

        assert_eq!(neighbors.len(), 3);
        assert!(neighbors.contains(&0));
        assert!(neighbors.contains(&2));
        assert!(neighbors.contains(&3));
    }

    #[test]
    fn topology_neighbors_of_returns_empty_for_isolated_atom() {
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        residue.add_atom(Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0)));
        residue.add_atom(Atom::new("CB", Element::C, Point::new(1.0, 0.0, 0.0)));
        chain.add_residue(residue);
        structure.add_chain(chain);

        let bonds = vec![Bond::new(0, 1, BondOrder::Single)];
        let topology = Topology::new(structure, bonds);

        let neighbors: Vec<_> = topology.neighbors_of(5).collect();

        assert!(neighbors.is_empty());
    }

    #[test]
    fn topology_display_formats_correctly() {
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        residue.add_atom(Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0)));
        residue.add_atom(Atom::new("CB", Element::C, Point::new(1.0, 0.0, 0.0)));
        chain.add_residue(residue);
        structure.add_chain(chain);

        let bonds = vec![Bond::new(0, 1, BondOrder::Single)];
        let topology = Topology::new(structure, bonds);

        let display = format!("{}", topology);
        let expected = "Topology { atoms: 2, bonds: 1 }";

        assert_eq!(display, expected);
    }

    #[test]
    fn topology_display_formats_empty_topology_correctly() {
        let structure = Structure::new();
        let topology = Topology::new(structure, Vec::new());

        let display = format!("{}", topology);
        let expected = "Topology { atoms: 0, bonds: 0 }";

        assert_eq!(display, expected);
    }

    #[test]
    fn topology_clone_creates_identical_copy() {
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        residue.add_atom(Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0)));
        residue.add_atom(Atom::new("CB", Element::C, Point::new(1.0, 0.0, 0.0)));
        chain.add_residue(residue);
        structure.add_chain(chain);

        let bonds = vec![Bond::new(0, 1, BondOrder::Single)];
        let topology = Topology::new(structure, bonds);

        let cloned = topology.clone();

        assert_eq!(topology.atom_count(), cloned.atom_count());
        assert_eq!(topology.bond_count(), cloned.bond_count());
        assert_eq!(topology.bonds(), cloned.bonds());
    }

    #[test]
    fn bond_partial_eq_compares_correctly() {
        let bond1 = Bond::new(0, 1, BondOrder::Single);
        let bond2 = Bond::new(1, 0, BondOrder::Single);
        let bond3 = Bond::new(0, 2, BondOrder::Single);
        let bond4 = Bond::new(0, 1, BondOrder::Double);

        assert_eq!(bond1, bond2);
        assert_ne!(bond1, bond3);
        assert_ne!(bond1, bond4);
    }

    #[test]
    fn bond_hash_considers_canonical_ordering() {
        use std::collections::HashSet;

        let bond1 = Bond::new(0, 1, BondOrder::Single);
        let bond2 = Bond::new(1, 0, BondOrder::Single);
        let mut set = HashSet::new();

        set.insert(bond1);
        set.insert(bond2);

        assert_eq!(set.len(), 1);
    }

    #[test]
    fn topology_with_complex_structure() {
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");
        let mut residue = Residue::new(
            1,
            None,
            "HOH",
            Some(StandardResidue::HOH),
            ResidueCategory::Standard,
        );
        residue.add_atom(Atom::new("O", Element::O, Point::new(0.0, 0.0, 0.0)));
        residue.add_atom(Atom::new("H1", Element::H, Point::new(0.96, 0.0, 0.0)));
        residue.add_atom(Atom::new("H2", Element::H, Point::new(-0.24, 0.93, 0.0)));
        chain.add_residue(residue);
        structure.add_chain(chain);

        let bonds = vec![
            Bond::new(0, 1, BondOrder::Single),
            Bond::new(0, 2, BondOrder::Single),
        ];
        let topology = Topology::new(structure, bonds);

        assert_eq!(topology.atom_count(), 3);
        assert_eq!(topology.bond_count(), 2);

        let o_neighbors: Vec<_> = topology.neighbors_of(0).collect();
        assert_eq!(o_neighbors.len(), 2);
        assert!(o_neighbors.contains(&1));
        assert!(o_neighbors.contains(&2));

        let h1_neighbors: Vec<_> = topology.neighbors_of(1).collect();
        assert_eq!(h1_neighbors, vec![0]);

        let h2_neighbors: Vec<_> = topology.neighbors_of(2).collect();
        assert_eq!(h2_neighbors, vec![0]);
    }

    #[test]
    fn topology_bonds_of_with_aromatic_bond() {
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        residue.add_atom(Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0)));
        residue.add_atom(Atom::new("CB", Element::C, Point::new(1.0, 0.0, 0.0)));
        chain.add_residue(residue);
        structure.add_chain(chain);

        let bonds = vec![Bond::new(0, 1, BondOrder::Aromatic)];
        let topology = Topology::new(structure, bonds);

        let bonds_of_0: Vec<_> = topology.bonds_of(0).collect();

        assert_eq!(bonds_of_0.len(), 1);
        assert_eq!(bonds_of_0[0].order, BondOrder::Aromatic);
    }
}
