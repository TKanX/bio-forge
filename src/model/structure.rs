use super::chain::Chain;
use super::residue::Residue;
use super::types::Point;
use std::fmt;

#[derive(Debug, Clone, Default)]
pub struct Structure {
    chains: Vec<Chain>,
    pub box_vectors: Option<[[f64; 3]; 3]>,
}

impl Structure {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_chain(&mut self, chain: Chain) {
        debug_assert!(
            self.chain(&chain.id).is_none(),
            "Attempted to add a duplicate chain ID '{}'",
            chain.id
        );
        self.chains.push(chain);
    }

    pub fn remove_chain(&mut self, id: &str) -> Option<Chain> {
        if let Some(index) = self.chains.iter().position(|c| c.id == id) {
            Some(self.chains.remove(index))
        } else {
            None
        }
    }

    pub fn clear(&mut self) {
        self.chains.clear();
    }

    pub fn chain(&self, id: &str) -> Option<&Chain> {
        self.chains.iter().find(|c| c.id == id)
    }

    pub fn chain_mut(&mut self, id: &str) -> Option<&mut Chain> {
        self.chains.iter_mut().find(|c| c.id == id)
    }

    pub fn find_residue(
        &self,
        chain_id: &str,
        residue_id: i32,
        insertion_code: Option<char>,
    ) -> Option<&Residue> {
        self.chain(chain_id)
            .and_then(|c| c.residue(residue_id, insertion_code))
    }

    pub fn find_residue_mut(
        &mut self,
        chain_id: &str,
        residue_id: i32,
        insertion_code: Option<char>,
    ) -> Option<&mut Residue> {
        self.chain_mut(chain_id)
            .and_then(|c| c.residue_mut(residue_id, insertion_code))
    }

    pub fn sort_chains_by_id(&mut self) {
        self.chains.sort_by(|a, b| a.id.cmp(&b.id));
    }

    pub fn chain_count(&self) -> usize {
        self.chains.len()
    }

    pub fn residue_count(&self) -> usize {
        self.chains.iter().map(|c| c.residue_count()).sum()
    }

    pub fn atom_count(&self) -> usize {
        self.chains.iter().map(|c| c.iter_atoms().count()).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.chains.is_empty()
    }

    pub fn iter_chains(&self) -> std::slice::Iter<'_, Chain> {
        self.chains.iter()
    }

    pub fn iter_chains_mut(&mut self) -> std::slice::IterMut<'_, Chain> {
        self.chains.iter_mut()
    }

    pub fn iter_atoms(&self) -> impl Iterator<Item = &super::atom::Atom> {
        self.chains.iter().flat_map(|c| c.iter_atoms())
    }

    pub fn iter_atoms_mut(&mut self) -> impl Iterator<Item = &mut super::atom::Atom> {
        self.chains.iter_mut().flat_map(|c| c.iter_atoms_mut())
    }

    pub fn iter_atoms_with_context(
        &self,
    ) -> impl Iterator<Item = (&Chain, &Residue, &super::atom::Atom)> {
        self.chains.iter().flat_map(|chain| {
            chain.iter_residues().flat_map(move |residue| {
                residue.iter_atoms().map(move |atom| (chain, residue, atom))
            })
        })
    }

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
            chain.id = "MODIFIED".to_string();
        }

        assert_eq!(structure.chain("MODIFIED").unwrap().id, "MODIFIED");
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

        assert_eq!(contexts, vec![("A".to_string(), 1, "CA".to_string())]);
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
}
