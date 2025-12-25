//! Representation of multi-chain biomolecular assemblies with geometric helpers.
//!
//! The `Structure` type aggregates polymer chains, provides fast lookup utilities, and
//! offers derived properties such as geometric centers or mass-weighted centroids. It is
//! the central container consumed by IO readers, cleaning operations, and solvation tools.

use super::chain::Chain;
use super::grid::Grid;
use super::residue::Residue;
use super::types::Point;
use crate::utils::parallel::*;
use std::fmt;

/// High-level biomolecular assembly composed of zero or more chains.
///
/// A `Structure` wraps individual chains, tracks optional periodic box vectors, and offers
/// convenience iterators for traversing chains, residues, and atoms alongside contextual
/// metadata. Builders and operations mutate the structure to clean, solvate, or analyze
/// biological systems.
#[derive(Debug, Clone, Default)]
pub struct Structure {
    /// Internal collection of polymer chains preserving insertion order.
    chains: Vec<Chain>,
    /// Optional periodic box represented as crystallographic basis vectors.
    pub box_vectors: Option<[[f64; 3]; 3]>,
}

impl Structure {
    /// Creates an empty structure with no chains or box vectors.
    ///
    /// # Returns
    ///
    /// A new `Structure` identical to `Structure::default()`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Appends a chain to the structure, asserting unique chain IDs in debug builds.
    ///
    /// The chain is inserted at the end of the current collection and becomes visible to
    /// iterator methods immediately.
    ///
    /// # Arguments
    ///
    /// * `chain` - Chain instance whose `id` must be unique within the structure.
    pub fn add_chain(&mut self, chain: Chain) {
        debug_assert!(
            self.chain(&chain.id).is_none(),
            "Attempted to add a duplicate chain ID '{}'",
            chain.id
        );
        self.chains.push(chain);
    }

    /// Removes and returns a chain by identifier if it exists.
    ///
    /// # Arguments
    ///
    /// * `id` - Chain identifier to search for.
    ///
    /// # Returns
    ///
    /// `Some(Chain)` when a chain with the provided ID is present, otherwise `None`.
    pub fn remove_chain(&mut self, id: &str) -> Option<Chain> {
        if let Some(index) = self.chains.iter().position(|c| c.id == id) {
            Some(self.chains.remove(index))
        } else {
            None
        }
    }

    /// Drops every chain from the structure, leaving box vectors untouched.
    pub fn clear(&mut self) {
        self.chains.clear();
    }

    /// Retrieves an immutable chain by identifier.
    ///
    /// # Arguments
    ///
    /// * `id` - Chain identifier to search for.
    ///
    /// # Returns
    ///
    /// `Some(&Chain)` if found, otherwise `None`.
    pub fn chain(&self, id: &str) -> Option<&Chain> {
        self.chains.iter().find(|c| c.id == id)
    }

    /// Retrieves a mutable chain by identifier.
    ///
    /// # Arguments
    ///
    /// * `id` - Chain identifier to search for.
    ///
    /// # Returns
    ///
    /// `Some(&mut Chain)` if found, otherwise `None`.
    pub fn chain_mut(&mut self, id: &str) -> Option<&mut Chain> {
        self.chains.iter_mut().find(|c| c.id == id)
    }

    /// Finds a residue using chain ID, residue number, and optional insertion code.
    ///
    /// # Arguments
    ///
    /// * `chain_id` - Identifier of the chain to search.
    /// * `residue_id` - Numeric residue index (typically PDB `resSeq`).
    /// * `insertion_code` - Optional insertion code differentiating duplicate IDs.
    ///
    /// # Returns
    ///
    /// `Some(&Residue)` when the residue is located, otherwise `None`.
    pub fn find_residue(
        &self,
        chain_id: &str,
        residue_id: i32,
        insertion_code: Option<char>,
    ) -> Option<&Residue> {
        self.chain(chain_id)
            .and_then(|c| c.residue(residue_id, insertion_code))
    }

    /// Finds a mutable residue reference using chain and residue identifiers.
    ///
    /// # Arguments
    ///
    /// * `chain_id` - Identifier of the chain to search.
    /// * `residue_id` - Numeric residue index.
    /// * `insertion_code` - Optional insertion code to disambiguate residues.
    ///
    /// # Returns
    ///
    /// `Some(&mut Residue)` when located, otherwise `None`.
    pub fn find_residue_mut(
        &mut self,
        chain_id: &str,
        residue_id: i32,
        insertion_code: Option<char>,
    ) -> Option<&mut Residue> {
        self.chain_mut(chain_id)
            .and_then(|c| c.residue_mut(residue_id, insertion_code))
    }

    /// Sorts chains lexicographically by their identifier.
    pub fn sort_chains_by_id(&mut self) {
        self.chains.sort_by(|a, b| a.id.cmp(&b.id));
    }

    /// Returns the number of chains currently stored.
    ///
    /// # Returns
    ///
    /// Chain count as `usize`.
    pub fn chain_count(&self) -> usize {
        self.chains.len()
    }

    /// Counts all residues across every chain.
    ///
    /// # Returns
    ///
    /// Total residue count as `usize`.
    pub fn residue_count(&self) -> usize {
        self.chains.iter().map(|c| c.residue_count()).sum()
    }

    /// Counts all atoms across every chain.
    ///
    /// # Returns
    ///
    /// Total atom count as `usize`.
    pub fn atom_count(&self) -> usize {
        self.chains.iter().map(|c| c.iter_atoms().count()).sum()
    }

    /// Indicates whether the structure contains zero chains.
    ///
    /// # Returns
    ///
    /// `true` if no chains are present.
    pub fn is_empty(&self) -> bool {
        self.chains.is_empty()
    }

    /// Provides an iterator over immutable chains.
    ///
    /// # Returns
    ///
    /// `std::slice::Iter<'_, Chain>` spanning all chains in insertion order.
    pub fn iter_chains(&self) -> std::slice::Iter<'_, Chain> {
        self.chains.iter()
    }

    /// Provides an iterator over mutable chains.
    ///
    /// # Returns
    ///
    /// `std::slice::IterMut<'_, Chain>` for in-place modification of chains.
    pub fn iter_chains_mut(&mut self) -> std::slice::IterMut<'_, Chain> {
        self.chains.iter_mut()
    }

    /// Provides a parallel iterator over immutable chains.
    ///
    /// # Returns
    ///
    /// A parallel iterator yielding `&Chain`.
    #[cfg(feature = "parallel")]
    pub fn par_chains(&self) -> impl IndexedParallelIterator<Item = &Chain> {
        self.chains.par_iter()
    }

    /// Provides a parallel iterator over immutable chains (internal fallback).
    #[cfg(not(feature = "parallel"))]
    pub(crate) fn par_chains(&self) -> impl IndexedParallelIterator<Item = &Chain> {
        self.chains.par_iter()
    }

    /// Provides a parallel iterator over mutable chains.
    ///
    /// # Returns
    ///
    /// A parallel iterator yielding `&mut Chain`.
    #[cfg(feature = "parallel")]
    pub fn par_chains_mut(&mut self) -> impl IndexedParallelIterator<Item = &mut Chain> {
        self.chains.par_iter_mut()
    }

    /// Provides a parallel iterator over mutable chains (internal fallback).
    #[cfg(not(feature = "parallel"))]
    pub(crate) fn par_chains_mut(&mut self) -> impl IndexedParallelIterator<Item = &mut Chain> {
        self.chains.par_iter_mut()
    }

    /// Provides a parallel iterator over immutable residues across all chains.
    ///
    /// # Returns
    ///
    /// A parallel iterator yielding `&Residue`.
    #[cfg(feature = "parallel")]
    pub fn par_residues(&self) -> impl ParallelIterator<Item = &Residue> {
        self.chains.par_iter().flat_map(|c| c.par_residues())
    }

    /// Provides a parallel iterator over immutable residues across all chains (internal fallback).
    #[cfg(not(feature = "parallel"))]
    pub(crate) fn par_residues(&self) -> impl ParallelIterator<Item = &Residue> {
        self.chains.par_iter().flat_map(|c| c.par_residues())
    }

    /// Provides a parallel iterator over mutable residues across all chains.
    ///
    /// # Returns
    ///
    /// A parallel iterator yielding `&mut Residue`.
    #[cfg(feature = "parallel")]
    pub fn par_residues_mut(&mut self) -> impl ParallelIterator<Item = &mut Residue> {
        self.chains
            .par_iter_mut()
            .flat_map(|c| c.par_residues_mut())
    }

    /// Provides a parallel iterator over mutable residues across all chains (internal fallback).
    #[cfg(not(feature = "parallel"))]
    pub(crate) fn par_residues_mut(&mut self) -> impl ParallelIterator<Item = &mut Residue> {
        self.chains
            .par_iter_mut()
            .flat_map(|c| c.par_residues_mut())
    }

    /// Provides a parallel iterator over immutable atoms across all chains.
    ///
    /// # Returns
    ///
    /// A parallel iterator yielding `&Atom`.
    #[cfg(feature = "parallel")]
    pub fn par_atoms(&self) -> impl ParallelIterator<Item = &super::atom::Atom> {
        self.chains
            .par_iter()
            .flat_map(|c| c.par_residues().flat_map(|r| r.par_atoms()))
    }

    /// Provides a parallel iterator over immutable atoms across all chains (internal fallback).
    #[cfg(not(feature = "parallel"))]
    pub(crate) fn par_atoms(&self) -> impl ParallelIterator<Item = &super::atom::Atom> {
        self.chains
            .par_iter()
            .flat_map(|c| c.par_residues().flat_map(|r| r.par_atoms()))
    }

    /// Provides a parallel iterator over mutable atoms across all chains.
    ///
    /// # Returns
    ///
    /// A parallel iterator yielding `&mut Atom`.
    #[cfg(feature = "parallel")]
    pub fn par_atoms_mut(&mut self) -> impl ParallelIterator<Item = &mut super::atom::Atom> {
        self.chains
            .par_iter_mut()
            .flat_map(|c| c.par_residues_mut().flat_map(|r| r.par_atoms_mut()))
    }

    /// Provides a parallel iterator over mutable atoms across all chains (internal fallback).
    #[cfg(not(feature = "parallel"))]
    pub(crate) fn par_atoms_mut(&mut self) -> impl ParallelIterator<Item = &mut super::atom::Atom> {
        self.chains
            .par_iter_mut()
            .flat_map(|c| c.par_residues_mut().flat_map(|r| r.par_atoms_mut()))
    }

    /// Iterates over immutable atoms across all chains.
    ///
    /// # Returns
    ///
    /// An iterator yielding `&Atom` in chain/residue order.
    pub fn iter_atoms(&self) -> impl Iterator<Item = &super::atom::Atom> {
        self.chains.iter().flat_map(|c| c.iter_atoms())
    }

    /// Iterates over mutable atoms across all chains.
    ///
    /// # Returns
    ///
    /// An iterator yielding `&mut Atom` in chain/residue order.
    pub fn iter_atoms_mut(&mut self) -> impl Iterator<Item = &mut super::atom::Atom> {
        self.chains.iter_mut().flat_map(|c| c.iter_atoms_mut())
    }

    /// Retains residues that satisfy a predicate, removing all others.
    ///
    /// The predicate receives the chain ID and a residue reference, enabling
    /// context-sensitive filtering.
    ///
    /// # Arguments
    ///
    /// * `f` - Closure returning `true` to keep the residue.
    pub fn retain_residues<F>(&mut self, mut f: F)
    where
        F: FnMut(&str, &Residue) -> bool,
    {
        for chain in &mut self.chains {
            let chain_id = chain.id.clone();
            chain.retain_residues(|residue| f(&chain_id, residue));
        }
    }

    /// Retains residues that satisfy a predicate, removing all others (Mutable version).
    ///
    /// The predicate receives the chain ID and a mutable residue reference.
    ///
    /// # Arguments
    ///
    /// * `f` - Closure returning `true` to keep the residue.
    pub fn retain_residues_mut<F>(&mut self, mut f: F)
    where
        F: FnMut(&str, &mut Residue) -> bool,
    {
        for chain in &mut self.chains {
            let chain_id = chain.id.clone();
            chain.retain_residues_mut(|residue| f(&chain_id, residue));
        }
    }

    /// Retains residues that satisfy a predicate, removing all others (Parallel version).
    ///
    /// This method processes chains in parallel when the `parallel` feature is enabled.
    /// The predicate must be thread-safe (`Sync` + `Send`) and immutable (`Fn`).
    ///
    /// # Arguments
    ///
    /// * `f` - Thread-safe closure returning `true` to keep the residue.
    #[cfg(feature = "parallel")]
    pub fn par_retain_residues<F>(&mut self, f: F)
    where
        F: Fn(&str, &Residue) -> bool + Sync + Send,
    {
        self.chains.par_iter_mut().for_each(|chain| {
            let chain_id = chain.id.clone();
            chain.retain_residues(|residue| f(&chain_id, residue));
        });
    }

    /// Retains residues that satisfy a predicate, removing all others (Sequential fallback).
    #[cfg(not(feature = "parallel"))]
    pub fn par_retain_residues<F>(&mut self, f: F)
    where
        F: Fn(&str, &Residue) -> bool + Sync + Send,
    {
        self.retain_residues(f);
    }

    /// Retains residues that satisfy a predicate, removing all others (Parallel Mutable version).
    ///
    /// # Arguments
    ///
    /// * `f` - Thread-safe closure returning `true` to keep the residue.
    #[cfg(feature = "parallel")]
    pub fn par_retain_residues_mut<F>(&mut self, f: F)
    where
        F: Fn(&str, &mut Residue) -> bool + Sync + Send,
    {
        self.chains.par_iter_mut().for_each(|chain| {
            let chain_id = chain.id.clone();
            chain.retain_residues_mut(|residue| f(&chain_id, residue));
        });
    }

    /// Retains residues that satisfy a predicate, removing all others (Sequential Mutable fallback).
    #[cfg(not(feature = "parallel"))]
    pub fn par_retain_residues_mut<F>(&mut self, f: F)
    where
        F: Fn(&str, &mut Residue) -> bool + Sync + Send,
    {
        self.retain_residues_mut(f);
    }

    /// Removes any chain that became empty after residue pruning.
    pub fn prune_empty_chains(&mut self) {
        self.chains.retain(|chain| !chain.is_empty());
    }

    /// Iterates over atoms while including chain and residue context.
    ///
    /// # Returns
    ///
    /// An iterator yielding triples `(&Chain, &Residue, &Atom)` for every atom.
    pub fn iter_atoms_with_context(
        &self,
    ) -> impl Iterator<Item = (&Chain, &Residue, &super::atom::Atom)> {
        self.chains.iter().flat_map(|chain| {
            chain.iter_residues().flat_map(move |residue| {
                residue.iter_atoms().map(move |atom| (chain, residue, atom))
            })
        })
    }

    /// Computes the geometric center of all atom coordinates.
    ///
    /// Falls back to the origin when the structure contains no atoms.
    ///
    /// # Returns
    ///
    /// A `Point` located at the unweighted centroid.
    pub fn geometric_center(&self) -> Point {
        let mut sum = nalgebra::Vector3::zeros();
        let mut count = 0;

        for atom in self.iter_atoms() {
            sum += atom.pos.coords;
            count += 1;
        }

        if count > 0 {
            Point::from(sum / (count as f64))
        } else {
            Point::origin()
        }
    }

    /// Computes the mass-weighted center of all atoms.
    ///
    /// Uses element atomic masses and returns the origin when the total mass is below
    /// numerical tolerance.
    ///
    /// # Returns
    ///
    /// A `Point` representing the center of mass.
    pub fn center_of_mass(&self) -> Point {
        let mut total_mass = 0.0;
        let mut weighted_sum = nalgebra::Vector3::zeros();

        for atom in self.iter_atoms() {
            let mass = atom.element.atomic_mass();
            weighted_sum += atom.pos.coords * mass;
            total_mass += mass;
        }

        if total_mass > 1e-9 {
            Point::from(weighted_sum / total_mass)
        } else {
            Point::origin()
        }
    }

    /// Constructs a spatial grid indexing all atoms in the structure.
    ///
    /// The grid stores `(chain_idx, residue_idx, atom_idx)` tuples, allowing efficient
    /// retrieval of atoms within a spatial neighborhood.
    ///
    /// # Arguments
    ///
    /// * `cell_size` - The side length of each spatial bin.
    ///
    /// # Returns
    ///
    /// A [`Grid`] containing all atoms in the structure.
    pub fn spatial_grid(&self, cell_size: f64) -> Grid<(usize, usize, usize)> {
        let items: Vec<_> = self
            .par_chains()
            .enumerate()
            .flat_map(|(c_idx, chain)| {
                chain
                    .par_residues()
                    .enumerate()
                    .flat_map_iter(move |(r_idx, residue)| {
                        residue
                            .iter_atoms()
                            .enumerate()
                            .map(move |(a_idx, atom)| (atom.pos, (c_idx, r_idx, a_idx)))
                    })
            })
            .collect();
        Grid::new(items, cell_size)
    }
}

impl fmt::Display for Structure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Structure {{ chains: {}, residues: {}, atoms: {} }}",
            self.chain_count(),
            self.residue_count(),
            self.atom_count()
        )
    }
}

impl FromIterator<Chain> for Structure {
    fn from_iter<T: IntoIterator<Item = Chain>>(iter: T) -> Self {
        Self {
            chains: iter.into_iter().collect(),
            box_vectors: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::atom::Atom;
    use crate::model::types::{Element, ResidueCategory, StandardResidue};

    fn make_residue(id: i32, name: &str) -> Residue {
        Residue::new(
            id,
            None,
            name,
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        )
    }

    #[test]
    fn structure_new_creates_empty_structure() {
        let structure = Structure::new();

        assert!(structure.is_empty());
        assert_eq!(structure.chain_count(), 0);
        assert_eq!(structure.residue_count(), 0);
        assert_eq!(structure.atom_count(), 0);
        assert!(structure.box_vectors.is_none());
    }

    #[test]
    fn structure_default_creates_empty_structure() {
        let structure = Structure::default();

        assert!(structure.is_empty());
        assert!(structure.box_vectors.is_none());
    }

    #[test]
    fn structure_add_chain_adds_chain_correctly() {
        let mut structure = Structure::new();
        let chain = Chain::new("A");

        structure.add_chain(chain);

        assert_eq!(structure.chain_count(), 1);
        assert!(structure.chain("A").is_some());
        assert_eq!(structure.chain("A").unwrap().id, "A");
    }

    #[test]
    fn structure_remove_chain_removes_existing_chain() {
        let mut structure = Structure::new();
        let chain = Chain::new("A");
        structure.add_chain(chain);

        let removed = structure.remove_chain("A");

        assert!(removed.is_some());
        assert_eq!(removed.unwrap().id, "A");
        assert_eq!(structure.chain_count(), 0);
        assert!(structure.chain("A").is_none());
    }

    #[test]
    fn structure_remove_chain_returns_none_for_nonexistent_chain() {
        let mut structure = Structure::new();

        let removed = structure.remove_chain("NONEXISTENT");

        assert!(removed.is_none());
    }

    #[test]
    fn structure_clear_removes_all_chains() {
        let mut structure = Structure::new();
        structure.add_chain(Chain::new("A"));
        structure.add_chain(Chain::new("B"));

        structure.clear();

        assert!(structure.is_empty());
        assert_eq!(structure.chain_count(), 0);
    }

    #[test]
    fn structure_chain_returns_correct_chain() {
        let mut structure = Structure::new();
        let chain = Chain::new("A");
        structure.add_chain(chain);

        let retrieved = structure.chain("A");

        assert!(retrieved.is_some());
        assert_eq!(retrieved.unwrap().id, "A");
    }

    #[test]
    fn structure_chain_returns_none_for_nonexistent_chain() {
        let structure = Structure::new();

        let retrieved = structure.chain("NONEXISTENT");

        assert!(retrieved.is_none());
    }

    #[test]
    fn structure_chain_mut_returns_correct_mutable_chain() {
        let mut structure = Structure::new();
        let chain = Chain::new("A");
        structure.add_chain(chain);

        let retrieved = structure.chain_mut("A");

        assert!(retrieved.is_some());
        assert_eq!(retrieved.unwrap().id, "A");
    }

    #[test]
    fn structure_chain_mut_returns_none_for_nonexistent_chain() {
        let mut structure = Structure::new();

        let retrieved = structure.chain_mut("NONEXISTENT");

        assert!(retrieved.is_none());
    }

    #[test]
    fn structure_find_residue_finds_correct_residue() {
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");
        let residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        chain.add_residue(residue);
        structure.add_chain(chain);

        let found = structure.find_residue("A", 1, None);

        assert!(found.is_some());
        assert_eq!(found.unwrap().id, 1);
        assert_eq!(found.unwrap().name, "ALA");
    }

    #[test]
    fn structure_find_residue_returns_none_for_nonexistent_chain() {
        let structure = Structure::new();

        let found = structure.find_residue("NONEXISTENT", 1, None);

        assert!(found.is_none());
    }

    #[test]
    fn structure_find_residue_returns_none_for_nonexistent_residue() {
        let mut structure = Structure::new();
        let chain = Chain::new("A");
        structure.add_chain(chain);

        let found = structure.find_residue("A", 999, None);

        assert!(found.is_none());
    }

    #[test]
    fn structure_find_residue_mut_finds_correct_mutable_residue() {
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");
        let residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        chain.add_residue(residue);
        structure.add_chain(chain);

        let found = structure.find_residue_mut("A", 1, None);

        assert!(found.is_some());
        assert_eq!(found.unwrap().id, 1);
    }

    #[test]
    fn structure_sort_chains_by_id_sorts_correctly() {
        let mut structure = Structure::new();
        structure.add_chain(Chain::new("C"));
        structure.add_chain(Chain::new("A"));
        structure.add_chain(Chain::new("B"));

        structure.sort_chains_by_id();

        let ids: Vec<&str> = structure.iter_chains().map(|c| c.id.as_str()).collect();
        assert_eq!(ids, vec!["A", "B", "C"]);
    }

    #[test]
    fn structure_chain_count_returns_correct_count() {
        let mut structure = Structure::new();

        assert_eq!(structure.chain_count(), 0);

        structure.add_chain(Chain::new("A"));
        assert_eq!(structure.chain_count(), 1);

        structure.add_chain(Chain::new("B"));
        assert_eq!(structure.chain_count(), 2);
    }

    #[test]
    fn structure_residue_count_returns_correct_count() {
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");
        chain.add_residue(Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        ));
        chain.add_residue(Residue::new(
            2,
            None,
            "GLY",
            Some(StandardResidue::GLY),
            ResidueCategory::Standard,
        ));
        structure.add_chain(chain);

        assert_eq!(structure.residue_count(), 2);
    }

    #[test]
    fn structure_atom_count_returns_correct_count() {
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

        assert_eq!(structure.atom_count(), 2);
    }

    #[test]
    fn structure_is_empty_returns_true_for_empty_structure() {
        let structure = Structure::new();

        assert!(structure.is_empty());
    }

    #[test]
    fn structure_is_empty_returns_false_for_non_empty_structure() {
        let mut structure = Structure::new();
        structure.add_chain(Chain::new("A"));

        assert!(!structure.is_empty());
    }

    #[test]
    fn structure_iter_chains_iterates_correctly() {
        let mut structure = Structure::new();
        structure.add_chain(Chain::new("A"));
        structure.add_chain(Chain::new("B"));

        let mut ids = Vec::new();
        for chain in structure.iter_chains() {
            ids.push(chain.id.clone());
        }

        assert_eq!(ids, vec!["A", "B"]);
    }

    #[test]
    fn structure_iter_chains_mut_iterates_correctly() {
        let mut structure = Structure::new();
        structure.add_chain(Chain::new("A"));

        for chain in structure.iter_chains_mut() {
            chain.id = "MODIFIED".into();
        }

        assert_eq!(structure.chain("MODIFIED").unwrap().id, "MODIFIED");
    }

    #[test]
    fn structure_par_chains_iterates_correctly() {
        let mut structure = Structure::new();
        structure.add_chain(Chain::new("A"));
        structure.add_chain(Chain::new("B"));

        let ids: Vec<String> = structure.par_chains().map(|c| c.id.to_string()).collect();

        assert_eq!(ids, vec!["A", "B"]);
    }

    #[test]
    fn structure_par_chains_mut_iterates_correctly() {
        let mut structure = Structure::new();
        structure.add_chain(Chain::new("A"));
        structure.add_chain(Chain::new("B"));

        structure.par_chains_mut().for_each(|c| {
            c.id = format!("{}_MOD", c.id).into();
        });

        assert_eq!(structure.chain("A_MOD").unwrap().id, "A_MOD");
        assert_eq!(structure.chain("B_MOD").unwrap().id, "B_MOD");
    }

    #[test]
    fn structure_par_residues_iterates_correctly() {
        let mut structure = Structure::new();
        let mut chain_a = Chain::new("A");
        chain_a.add_residue(make_residue(1, "ALA"));
        let mut chain_b = Chain::new("B");
        chain_b.add_residue(make_residue(2, "GLY"));
        structure.add_chain(chain_a);
        structure.add_chain(chain_b);

        let count = structure.par_residues().count();
        assert_eq!(count, 2);

        let names: Vec<String> = structure
            .par_residues()
            .map(|r| r.name.to_string())
            .collect();
        assert!(names.contains(&"ALA".to_string()));
        assert!(names.contains(&"GLY".to_string()));
    }

    #[test]
    fn structure_par_residues_mut_iterates_correctly() {
        let mut structure = Structure::new();
        let mut chain_a = Chain::new("A");
        chain_a.add_residue(make_residue(1, "ALA"));
        let mut chain_b = Chain::new("B");
        chain_b.add_residue(make_residue(2, "GLY"));
        structure.add_chain(chain_a);
        structure.add_chain(chain_b);

        structure.par_residues_mut().for_each(|r| {
            r.name = format!("{}_MOD", r.name).into();
        });

        let chain_a = structure.chain("A").unwrap();
        assert_eq!(chain_a.residue(1, None).unwrap().name, "ALA_MOD");

        let chain_b = structure.chain("B").unwrap();
        assert_eq!(chain_b.residue(2, None).unwrap().name, "GLY_MOD");
    }

    #[test]
    fn structure_par_atoms_iterates_correctly() {
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");
        let mut r1 = make_residue(1, "ALA");
        r1.add_atom(Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0)));
        let mut r2 = make_residue(2, "GLY");
        r2.add_atom(Atom::new("N", Element::N, Point::new(1.0, 0.0, 0.0)));

        chain.add_residue(r1);
        chain.add_residue(r2);
        structure.add_chain(chain);

        let count = structure.par_atoms().count();
        assert_eq!(count, 2);

        let names: Vec<String> = structure.par_atoms().map(|a| a.name.to_string()).collect();
        assert!(names.contains(&"CA".to_string()));
        assert!(names.contains(&"N".to_string()));
    }

    #[test]
    fn structure_par_atoms_mut_iterates_correctly() {
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");
        let mut r1 = make_residue(1, "ALA");
        r1.add_atom(Atom::new("CA", Element::C, Point::new(0.0, 0.0, 0.0)));
        chain.add_residue(r1);
        structure.add_chain(chain);

        structure.par_atoms_mut().for_each(|a| {
            a.pos.x += 10.0;
        });

        let atom = structure
            .chain("A")
            .unwrap()
            .residue(1, None)
            .unwrap()
            .atom("CA")
            .unwrap();
        assert!((atom.pos.x - 10.0).abs() < 1e-6);
    }

    #[test]
    fn structure_iter_atoms_iterates_over_all_atoms() {
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

        let mut atom_names = Vec::new();
        for atom in structure.iter_atoms() {
            atom_names.push(atom.name.clone());
        }

        assert_eq!(atom_names, vec!["CA", "CB"]);
    }

    #[test]
    fn structure_iter_atoms_mut_iterates_over_all_atoms_mutably() {
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
        chain.add_residue(residue);
        structure.add_chain(chain);

        for atom in structure.iter_atoms_mut() {
            atom.translate_by(&nalgebra::Vector3::new(1.0, 0.0, 0.0));
        }

        let translated_atom = structure
            .find_residue("A", 1, None)
            .unwrap()
            .atom("CA")
            .unwrap();
        assert!((translated_atom.pos.x - 1.0).abs() < 1e-10);
    }

    #[test]
    fn structure_retain_residues_filters_using_chain_context() {
        let mut structure = Structure::new();
        let mut chain_a = Chain::new("A");
        chain_a.add_residue(make_residue(1, "ALA"));
        chain_a.add_residue(make_residue(2, "GLY"));
        let mut chain_b = Chain::new("B");
        chain_b.add_residue(make_residue(3, "SER"));
        structure.add_chain(chain_a);
        structure.add_chain(chain_b);

        structure.retain_residues(|chain_id, residue| chain_id == "A" && residue.id == 1);

        let chain_a = structure.chain("A").unwrap();
        assert_eq!(chain_a.residue_count(), 1);
        assert!(chain_a.residue(1, None).is_some());
        assert!(structure.chain("B").unwrap().is_empty());
    }

    #[test]
    fn structure_retain_residues_mut_filters_and_modifies() {
        let mut structure = Structure::new();
        let mut chain_a = Chain::new("A");
        chain_a.add_residue(make_residue(1, "ALA"));
        chain_a.add_residue(make_residue(2, "GLY"));
        let mut chain_b = Chain::new("B");
        chain_b.add_residue(make_residue(3, "SER"));
        structure.add_chain(chain_a);
        structure.add_chain(chain_b);

        structure.retain_residues_mut(|chain_id, residue| {
            if chain_id == "A" && residue.id == 1 {
                residue.name = format!("{}_MOD", residue.name).into();
                true
            } else {
                false
            }
        });

        let chain_a = structure.chain("A").unwrap();
        assert_eq!(chain_a.residue_count(), 1);
        assert_eq!(chain_a.residue(1, None).unwrap().name, "ALA_MOD");
        assert!(structure.chain("B").unwrap().is_empty());
    }

    #[test]
    fn structure_par_retain_residues_filters_correctly() {
        let mut structure = Structure::new();
        let mut chain_a = Chain::new("A");
        chain_a.add_residue(make_residue(1, "ALA"));
        chain_a.add_residue(make_residue(2, "GLY"));
        let mut chain_b = Chain::new("B");
        chain_b.add_residue(make_residue(3, "SER"));
        structure.add_chain(chain_a);
        structure.add_chain(chain_b);

        structure.par_retain_residues(|chain_id, residue| chain_id == "A" && residue.id == 1);

        let chain_a = structure.chain("A").unwrap();
        assert_eq!(chain_a.residue_count(), 1);
        assert!(chain_a.residue(1, None).is_some());
        assert!(structure.chain("B").unwrap().is_empty());
    }

    #[test]
    fn structure_par_retain_residues_mut_filters_and_modifies() {
        let mut structure = Structure::new();
        let mut chain_a = Chain::new("A");
        chain_a.add_residue(make_residue(1, "ALA"));
        chain_a.add_residue(make_residue(2, "GLY"));
        let mut chain_b = Chain::new("B");
        chain_b.add_residue(make_residue(3, "SER"));
        structure.add_chain(chain_a);
        structure.add_chain(chain_b);

        structure.par_retain_residues_mut(|chain_id, residue| {
            if chain_id == "A" && residue.id == 1 {
                residue.name = format!("{}_MOD", residue.name).into();
                true
            } else {
                false
            }
        });

        let chain_a = structure.chain("A").unwrap();
        assert_eq!(chain_a.residue_count(), 1);
        assert_eq!(chain_a.residue(1, None).unwrap().name, "ALA_MOD");
        assert!(structure.chain("B").unwrap().is_empty());
    }

    #[test]
    fn structure_prune_empty_chains_removes_them() {
        let mut structure = Structure::new();
        let mut chain_a = Chain::new("A");
        chain_a.add_residue(make_residue(1, "ALA"));
        let chain_b = Chain::new("B");
        structure.add_chain(chain_a);
        structure.add_chain(chain_b);

        structure.prune_empty_chains();

        assert!(structure.chain("A").is_some());
        assert!(structure.chain("B").is_none());
        assert_eq!(structure.chain_count(), 1);
    }

    #[test]
    fn structure_iter_atoms_with_context_provides_correct_context() {
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
        chain.add_residue(residue);
        structure.add_chain(chain);

        let mut contexts = Vec::new();
        for (chain, residue, atom) in structure.iter_atoms_with_context() {
            contexts.push((chain.id.clone(), residue.id, atom.name.clone()));
        }

        assert_eq!(contexts, vec![("A".into(), 1, "CA".into())]);
    }

    #[test]
    fn structure_geometric_center_calculates_correctly() {
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
        residue.add_atom(Atom::new("CB", Element::C, Point::new(2.0, 0.0, 0.0)));
        chain.add_residue(residue);
        structure.add_chain(chain);

        let center = structure.geometric_center();

        assert!((center.x - 1.0).abs() < 1e-10);
        assert!((center.y - 0.0).abs() < 1e-10);
        assert!((center.z - 0.0).abs() < 1e-10);
    }

    #[test]
    fn structure_geometric_center_returns_origin_for_empty_structure() {
        let structure = Structure::new();

        let center = structure.geometric_center();

        assert_eq!(center, Point::origin());
    }

    #[test]
    fn structure_center_of_mass_calculates_correctly() {
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        residue.add_atom(Atom::new("H", Element::H, Point::new(0.0, 0.0, 0.0)));
        residue.add_atom(Atom::new("C", Element::C, Point::new(2.0, 0.0, 0.0)));
        chain.add_residue(residue);
        structure.add_chain(chain);

        let com = structure.center_of_mass();

        let expected_x = (0.0 * 1.00794 + 2.0 * 12.0107) / (1.00794 + 12.0107);
        assert!((com.x - expected_x).abs() < 1e-3);
        assert!((com.y - 0.0).abs() < 1e-10);
        assert!((com.z - 0.0).abs() < 1e-10);
    }

    #[test]
    fn structure_center_of_mass_returns_origin_for_empty_structure() {
        let structure = Structure::new();

        let com = structure.center_of_mass();

        assert_eq!(com, Point::origin());
    }

    #[test]
    fn structure_display_formats_correctly() {
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
        chain.add_residue(residue);
        structure.add_chain(chain);

        let display = format!("{}", structure);
        let expected = "Structure { chains: 1, residues: 1, atoms: 1 }";

        assert_eq!(display, expected);
    }

    #[test]
    fn structure_display_formats_empty_structure_correctly() {
        let structure = Structure::new();

        let display = format!("{}", structure);
        let expected = "Structure { chains: 0, residues: 0, atoms: 0 }";

        assert_eq!(display, expected);
    }

    #[test]
    fn structure_from_iterator_creates_structure_correctly() {
        let chains = vec![Chain::new("A"), Chain::new("B")];
        let structure: Structure = chains.into_iter().collect();

        assert_eq!(structure.chain_count(), 2);
        assert!(structure.chain("A").is_some());
        assert!(structure.chain("B").is_some());
        assert!(structure.box_vectors.is_none());
    }

    #[test]
    fn structure_clone_creates_identical_copy() {
        let mut structure = Structure::new();
        structure.add_chain(Chain::new("A"));
        structure.box_vectors = Some([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);

        let cloned = structure.clone();

        assert_eq!(structure.chain_count(), cloned.chain_count());
        assert_eq!(structure.box_vectors, cloned.box_vectors);
    }

    #[test]
    fn spatial_grid_bins_and_neighbor_queries_work() {
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
        residue.add_atom(Atom::new("CB", Element::C, Point::new(1.5, 0.0, 0.0)));
        chain.add_residue(residue);
        structure.add_chain(chain);

        let grid = structure.spatial_grid(1.0);

        let big_center = structure.geometric_center();
        let all_count = grid.neighbors(&big_center, 1e6).count();
        assert_eq!(all_count, structure.atom_count());

        let center = Point::new(0.0, 0.0, 0.0);
        let neighbors: Vec<_> = grid.neighbors(&center, 0.1).collect();
        assert_eq!(neighbors.len(), 1);

        let &(c_idx, r_idx, a_idx) = neighbors[0];
        let chain_ref = structure.iter_chains().nth(c_idx).unwrap();
        let residue_ref = chain_ref.iter_residues().nth(r_idx).unwrap();
        let atom_ref = residue_ref.iter_atoms().nth(a_idx).unwrap();
        assert_eq!(atom_ref.name, "CA");

        let coarse = grid.neighbors(&center, 2.0).count();
        assert_eq!(coarse, 2);

        let exact: Vec<_> = grid.neighbors(&center, 1.0).exact().collect();
        assert_eq!(exact.len(), 1);
    }

    #[test]
    fn spatial_grid_empty_structure_is_empty() {
        let structure = Structure::new();

        let grid = structure.spatial_grid(1.0);
        assert_eq!(grid.neighbors(&Point::origin(), 1.0).count(), 0);
        assert_eq!(grid.neighbors(&Point::origin(), 1.0).exact().count(), 0);
    }

    #[test]
    fn spatial_grid_dense_packing_and_indices_are_consistent() {
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");
        let mut residue = Residue::new(
            1,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );

        for i in 0..50 {
            residue.add_atom(Atom::new(
                &format!("X{}", i),
                Element::C,
                Point::new(0.1, 0.1, 0.1),
            ));
        }
        chain.add_residue(residue);
        structure.add_chain(chain);

        let grid = structure.spatial_grid(1.0);
        let center = Point::new(0.1, 0.1, 0.1);
        assert_eq!(grid.neighbors(&center, 1e6).count(), structure.atom_count());

        let center = Point::new(0.1, 0.1, 0.1);
        let count = grid.neighbors(&center, 0.5).count();
        assert_eq!(count, 50);

        for &(c, r, a) in grid.neighbors(&center, 0.5) {
            let chain_ref = structure.iter_chains().nth(c).unwrap();
            let residue_ref = chain_ref.iter_residues().nth(r).unwrap();
            let atom_ref = residue_ref.iter_atoms().nth(a).unwrap();
            assert!(atom_ref.name.starts_with('X'));
        }
    }
}
