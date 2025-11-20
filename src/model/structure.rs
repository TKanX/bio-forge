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

    pub fn find_residue(&self, chain_id: &str, residue_id: i32) -> Option<&Residue> {
        self.chain(chain_id).and_then(|c| c.residue(residue_id))
    }

    pub fn find_residue_mut(&mut self, chain_id: &str, residue_id: i32) -> Option<&mut Residue> {
        self.chain_mut(chain_id)
            .and_then(|c| c.residue_mut(residue_id))
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
