//! Bond topology builder used to reconstruct intra- and inter-residue connectivity.
//!
//! The routines in this module consume a fully prepared [`Structure`] and emit a
//! [`Topology`] that encodes all covalent bonds. Internal templates are used for
//! standard residues, while callers can provide additional hetero templates.
//! Beyond template-driven intra-residue bonds, the builder also infers peptide,
//! nucleic-backbone, terminal, and disulfide bonds using geometric thresholds.

use crate::db;
use crate::model::{
    structure::Structure,
    template::Template,
    topology::{Bond, Topology},
    types::{BondOrder, ResidueCategory, ResiduePosition},
};
use crate::ops::error::Error;
use std::collections::HashMap;

/// Builder responsible for creating [`Topology`] objects from a [`Structure`].
///
/// The builder can augment the internal template database with additional
/// hetero templates and tweak geometric cutoffs for disulfide, peptide, and
/// nucleic bonds before running [`TopologyBuilder::build`].
pub struct TopologyBuilder {
    hetero_templates: HashMap<String, Template>,
    disulfide_bond_cutoff: f64,
    peptide_bond_cutoff: f64,
    nucleic_bond_cutoff: f64,
}

impl Default for TopologyBuilder {
    fn default() -> Self {
        Self {
            hetero_templates: HashMap::new(),
            disulfide_bond_cutoff: 2.2,
            peptide_bond_cutoff: 1.5,
            nucleic_bond_cutoff: 1.8,
        }
    }
}

impl TopologyBuilder {
    /// Creates a builder that uses default cutoffs and no extra templates.
    ///
    /// Default cutoffs were chosen to capture typical covalent geometry in
    /// Ångström. They can be overridden through the builder-style setters.
    pub fn new() -> Self {
        Self::default()
    }

    /// Registers a new template that will be used for hetero residues.
    ///
    /// # Arguments
    ///
    /// * `template` - Template describing atoms and bonds for a residue that
    ///   is classified as `ResidueCategory::Hetero`.
    pub fn add_hetero_template(mut self, template: Template) -> Self {
        self.hetero_templates
            .insert(template.name.clone(), template);
        self
    }

    /// Configures the maximum SG···SG distance allowed for disulfide bonds.
    ///
    /// # Arguments
    ///
    /// * `cutoff` - Maximum distance in Ångström for connecting cysteine SG
    ///   atoms.
    pub fn disulfide_cutoff(mut self, cutoff: f64) -> Self {
        self.disulfide_bond_cutoff = cutoff;
        self
    }

    /// Builds a [`Topology`] for the provided structure.
    ///
    /// All intra-residue bonds are taken from templates and terminal rules,
    /// while inter-residue bonds rely on distance checks using the configured
    /// cutoffs.
    ///
    /// # Errors
    ///
    /// Returns [`Error`] when a required template or atom is missing.
    pub fn build(self, structure: Structure) -> Result<Topology, Error> {
        let mut bonds = Vec::new();

        self.build_intra_residue(&structure, &mut bonds)?;

        self.build_inter_residue(&structure, &mut bonds)?;

        Ok(Topology::new(structure, bonds))
    }

    /// Populates bonds that lie within each residue using templates and
    /// terminal-specific heuristics.
    fn build_intra_residue(
        &self,
        structure: &Structure,
        bonds: &mut Vec<Bond>,
    ) -> Result<(), Error> {
        let mut global_atom_offset = 0;

        for chain in structure.iter_chains() {
            for residue in chain.iter_residues() {
                let atom_count = residue.atom_count();

                if residue.category == ResidueCategory::Ion {
                    global_atom_offset += atom_count;
                    continue;
                }

                if residue.category == ResidueCategory::Standard {
                    let tmpl_name = &residue.name;
                    let tmpl_view = db::get_template(tmpl_name).ok_or_else(|| {
                        Error::MissingInternalTemplate {
                            res_name: tmpl_name.clone(),
                        }
                    })?;

                    for (a1_name, a2_name, order) in tmpl_view.bonds() {
                        self.try_add_bond(
                            residue,
                            global_atom_offset,
                            a1_name,
                            a2_name,
                            order,
                            bonds,
                        )?;
                    }

                    for (h_name, _, anchors) in tmpl_view.hydrogens() {
                        match anchors.into_iter().next() {
                            Some(anchor) if residue.has_atom(h_name) => {
                                self.try_add_bond(
                                    residue,
                                    global_atom_offset,
                                    h_name,
                                    anchor,
                                    BondOrder::Single,
                                    bonds,
                                )?;
                            }
                            _ => {}
                        }
                    }

                    self.handle_terminal_intra_bonds(residue, global_atom_offset, bonds)?;
                } else if residue.category == ResidueCategory::Hetero {
                    let tmpl = self.hetero_templates.get(&residue.name).ok_or_else(|| {
                        Error::MissingHeteroTemplate {
                            res_name: residue.name.clone(),
                        }
                    })?;

                    for (a1_name, a2_name, order) in tmpl.bonds() {
                        self.try_add_bond(
                            residue,
                            global_atom_offset,
                            a1_name,
                            a2_name,
                            *order,
                            bonds,
                        )?;
                    }
                }

                global_atom_offset += atom_count;
            }
        }
        Ok(())
    }

    /// Attempts to add a bond and reports informative errors when atoms are
    /// missing.
    fn try_add_bond(
        &self,
        residue: &crate::model::residue::Residue,
        offset: usize,
        name1: &str,
        name2: &str,
        order: BondOrder,
        bonds: &mut Vec<Bond>,
    ) -> Result<(), Error> {
        let idx1 = residue.iter_atoms().position(|a| a.name == name1);
        let idx2 = residue.iter_atoms().position(|a| a.name == name2);

        match (idx1, idx2) {
            (Some(i1), Some(i2)) => {
                bonds.push(Bond::new(offset + i1, offset + i2, order));
                Ok(())
            }
            (None, _) if self.is_optional_terminal_atom(residue, name1) => Ok(()),
            (_, None) if self.is_optional_terminal_atom(residue, name2) => Ok(()),
            (None, _) => Err(Error::topology_atom_missing(
                &residue.name,
                residue.id,
                name1,
            )),
            (_, None) => Err(Error::topology_atom_missing(
                &residue.name,
                residue.id,
                name2,
            )),
        }
    }

    /// Returns whether the provided atom is optional for the residue position.
    fn is_optional_terminal_atom(
        &self,
        residue: &crate::model::residue::Residue,
        atom_name: &str,
    ) -> bool {
        let is_protein = residue.standard_name.is_some_and(|std| std.is_protein());
        let is_nucleic = residue.standard_name.is_some_and(|std| std.is_nucleic());

        match residue.position {
            ResiduePosition::NTerminal if is_protein => atom_name == "H",
            ResiduePosition::CTerminal if is_protein => matches!(atom_name, "HXT" | "HOXT"),
            ResiduePosition::FivePrime if is_nucleic => atom_name == "HO5'",
            ResiduePosition::ThreePrime if is_nucleic => atom_name == "HO3'",
            _ => false,
        }
    }

    /// Adds missing bonds for termini (protein and nucleic acid) that are not
    /// explicitly listed in templates.
    fn handle_terminal_intra_bonds(
        &self,
        residue: &crate::model::residue::Residue,
        offset: usize,
        bonds: &mut Vec<Bond>,
    ) -> Result<(), Error> {
        if residue.position == ResiduePosition::NTerminal
            && residue.standard_name.is_some_and(|s| s.is_protein())
        {
            for h_name in ["H1", "H2", "H3"] {
                if let (Some(h_idx), Some(n_idx)) = (
                    residue.iter_atoms().position(|a| a.name == h_name),
                    residue.iter_atoms().position(|a| a.name == "N"),
                ) {
                    bonds.push(Bond::new(offset + h_idx, offset + n_idx, BondOrder::Single));
                }
            }
        }

        if residue.position == ResiduePosition::CTerminal
            && residue.standard_name.is_some_and(|s| s.is_protein())
        {
            let c_idx = residue.iter_atoms().position(|a| a.name == "C");
            let oxt_idx = residue.iter_atoms().position(|a| a.name == "OXT");

            if let (Some(c_idx), Some(oxt_idx)) = (c_idx, oxt_idx) {
                bonds.push(Bond::new(
                    offset + c_idx,
                    offset + oxt_idx,
                    BondOrder::Single,
                ));

                for h_name in ["HXT", "HOXT"] {
                    if let Some(h_idx) = residue.iter_atoms().position(|a| a.name == h_name) {
                        bonds.push(Bond::new(
                            offset + oxt_idx,
                            offset + h_idx,
                            BondOrder::Single,
                        ));
                    }
                }
            }
        }

        if residue.position == ResiduePosition::FivePrime
            && residue.standard_name.is_some_and(|s| s.is_nucleic())
        {
            let ho5_idx = residue.iter_atoms().position(|a| a.name == "HO5'");
            let o5_idx = residue.iter_atoms().position(|a| a.name == "O5'");

            if let (Some(h_idx), Some(o_idx)) = (ho5_idx, o5_idx) {
                bonds.push(Bond::new(offset + h_idx, offset + o_idx, BondOrder::Single));
            }
        }

        if residue.position == ResiduePosition::ThreePrime
            && residue.standard_name.is_some_and(|s| s.is_nucleic())
        {
            let ho3_idx = residue.iter_atoms().position(|a| a.name == "HO3'");
            let o3_idx = residue.iter_atoms().position(|a| a.name == "O3'");

            if let (Some(h_idx), Some(o_idx)) = (ho3_idx, o3_idx) {
                bonds.push(Bond::new(offset + h_idx, offset + o_idx, BondOrder::Single));
            }
        }

        Ok(())
    }

    /// Connects adjacent residues through peptide, nucleic-backbone, and
    /// disulfide bonds based on geometry.
    fn build_inter_residue(
        &self,
        structure: &Structure,
        bonds: &mut Vec<Bond>,
    ) -> Result<(), Error> {
        let mut residue_offsets: Vec<Vec<usize>> = Vec::new();
        let mut current_offset = 0;

        for chain in structure.iter_chains() {
            let mut chain_offsets = Vec::new();
            for residue in chain.iter_residues() {
                chain_offsets.push(current_offset);
                current_offset += residue.atom_count();
            }
            residue_offsets.push(chain_offsets);
        }

        for (c_idx, chain) in structure.iter_chains().enumerate() {
            let residues: Vec<_> = chain.iter_residues().collect();
            if residues.len() < 2 {
                continue;
            }

            for i in 0..residues.len() - 1 {
                let curr = residues[i];
                let next = residues[i + 1];

                if curr.category != ResidueCategory::Standard
                    || next.category != ResidueCategory::Standard
                {
                    continue;
                }

                let curr_offset = residue_offsets[c_idx][i];
                let next_offset = residue_offsets[c_idx][i + 1];

                if let (Some(std1), Some(std2)) = (curr.standard_name, next.standard_name) {
                    if std1.is_protein() && std2.is_protein() {
                        self.connect_atoms_if_close(
                            AtomLocator::new(curr, curr_offset, "C"),
                            AtomLocator::new(next, next_offset, "N"),
                            self.peptide_bond_cutoff,
                            BondOrder::Single,
                            bonds,
                        );
                    } else if std1.is_nucleic() && std2.is_nucleic() {
                        self.connect_atoms_if_close(
                            AtomLocator::new(curr, curr_offset, "O3'"),
                            AtomLocator::new(next, next_offset, "P"),
                            self.nucleic_bond_cutoff,
                            BondOrder::Single,
                            bonds,
                        );
                    }
                }
            }
        }

        let mut sulfur_atoms = Vec::new();

        for (c_idx, chain) in structure.iter_chains().enumerate() {
            for (r_idx, residue) in chain.iter_residues().enumerate() {
                match residue.name.as_str() {
                    "CYX" | "CYM" => {
                        if let Some(sg) = residue.atom("SG") {
                            let offset = residue_offsets[c_idx][r_idx]
                                + residue.iter_atoms().position(|a| a.name == "SG").unwrap();
                            sulfur_atoms.push((offset, sg.pos));
                        }
                    }
                    _ => {}
                }
            }
        }

        let cutoff_sq = self.disulfide_bond_cutoff * self.disulfide_bond_cutoff;
        for i in 0..sulfur_atoms.len() {
            for j in (i + 1)..sulfur_atoms.len() {
                let (idx1, pos1) = sulfur_atoms[i];
                let (idx2, pos2) = sulfur_atoms[j];

                if nalgebra::distance_squared(&pos1, &pos2) <= cutoff_sq {
                    bonds.push(Bond::new(idx1, idx2, BondOrder::Single));
                }
            }
        }

        Ok(())
    }

    /// Adds a bond when the specified atoms are within the provided cutoff.
    fn connect_atoms_if_close(
        &self,
        first: AtomLocator<'_>,
        second: AtomLocator<'_>,
        cutoff: f64,
        order: BondOrder,
        bonds: &mut Vec<Bond>,
    ) {
        if let (Some(idx1), Some(idx2)) = (
            first
                .residue
                .iter_atoms()
                .position(|a| a.name == first.atom_name),
            second
                .residue
                .iter_atoms()
                .position(|a| a.name == second.atom_name),
        ) {
            let p1 = first.residue.atoms()[idx1].pos;
            let p2 = second.residue.atoms()[idx2].pos;

            if nalgebra::distance_squared(&p1, &p2) <= cutoff * cutoff {
                bonds.push(Bond::new(first.offset + idx1, second.offset + idx2, order));
            }
        }
    }
}

/// Utility that couples a residue reference with its global atom offset.
struct AtomLocator<'a> {
    residue: &'a crate::model::residue::Residue,
    offset: usize,
    atom_name: &'a str,
}

impl<'a> AtomLocator<'a> {
    /// Creates a locator describing a specific atom inside a residue.
    fn new(residue: &'a crate::model::residue::Residue, offset: usize, atom_name: &'a str) -> Self {
        Self {
            residue,
            offset,
            atom_name,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{
        atom::Atom,
        chain::Chain,
        residue::Residue,
        structure::Structure,
        template::Template,
        types::{Element, Point, ResidueCategory, ResiduePosition},
    };
    use nalgebra::Vector3;

    fn structure_from_residues(residues: Vec<Residue>) -> Structure {
        let mut chain = Chain::new("A");
        for residue in residues {
            chain.add_residue(residue);
        }

        let mut structure = Structure::new();
        structure.add_chain(chain);
        structure
    }

    fn standard_residue(name: &str, id: i32, position: ResiduePosition) -> Residue {
        let template = db::get_template(name)
            .unwrap_or_else(|| panic!("missing internal template for test residue '{name}'"));

        let mut residue = Residue::new(
            id,
            None,
            name,
            Some(template.standard_name()),
            ResidueCategory::Standard,
        );
        residue.position = position;

        for (atom_name, element, point) in template.heavy_atoms() {
            residue.add_atom(Atom::new(atom_name, element, point));
        }

        for (atom_name, point, _) in template.hydrogens() {
            residue.add_atom(Atom::new(atom_name, Element::H, point));
        }

        residue
    }

    fn translate_residue(residue: &mut Residue, offset: Vector3<f64>) {
        for atom in residue.iter_atoms_mut() {
            atom.translate_by(&offset);
        }
    }

    fn global_atom_index(
        topology: &Topology,
        chain_id: &str,
        residue_id: i32,
        atom_name: &str,
    ) -> usize {
        let mut idx = 0;
        for chain in topology.structure().iter_chains() {
            for residue in chain.iter_residues() {
                for atom in residue.iter_atoms() {
                    if chain.id == chain_id && residue.id == residue_id && atom.name == atom_name {
                        return idx;
                    }
                    idx += 1;
                }
            }
        }

        panic!(
            "atom '{}' not found in chain '{}' residue '{}'",
            atom_name, chain_id, residue_id
        );
    }

    fn has_bond(topology: &Topology, idx1: usize, idx2: usize, order: BondOrder) -> bool {
        let (a1, a2) = if idx1 <= idx2 {
            (idx1, idx2)
        } else {
            (idx2, idx1)
        };
        topology
            .bonds()
            .iter()
            .any(|bond| bond.a1_idx == a1 && bond.a2_idx == a2 && bond.order == order)
    }

    #[test]
    fn build_creates_peptide_bond_between_adjacent_proteins() {
        let residue1 = standard_residue("GLY", 1, ResiduePosition::NTerminal);
        let mut residue2 = standard_residue("ALA", 2, ResiduePosition::Internal);

        let c_pos = residue1.atom("C").unwrap().pos;
        let n_pos = residue2.atom("N").unwrap().pos;
        let target_n = c_pos + Vector3::new(1.33, 0.0, 0.0);
        translate_residue(&mut residue2, target_n - n_pos);

        let structure = structure_from_residues(vec![residue1, residue2]);
        let topology = TopologyBuilder::new()
            .build(structure)
            .expect("build topology");

        let c_idx = global_atom_index(&topology, "A", 1, "C");
        let n_idx = global_atom_index(&topology, "A", 2, "N");

        assert!(has_bond(&topology, c_idx, n_idx, BondOrder::Single));
    }

    #[test]
    fn build_creates_nucleic_backbone_bond() {
        let residue1 = standard_residue("DA", 1, ResiduePosition::FivePrime);
        let mut residue2 = standard_residue("DT", 2, ResiduePosition::ThreePrime);

        let o3_pos = residue1.atom("O3'").unwrap().pos;
        let p_pos = residue2.atom("P").unwrap().pos;
        let target_p = o3_pos + Vector3::new(0.0, 0.0, 1.6);
        translate_residue(&mut residue2, target_p - p_pos);

        let structure = structure_from_residues(vec![residue1, residue2]);
        let topology = TopologyBuilder::new()
            .build(structure)
            .expect("build topology");

        let o3_idx = global_atom_index(&topology, "A", 1, "O3'");
        let p_idx = global_atom_index(&topology, "A", 2, "P");

        assert!(has_bond(&topology, o3_idx, p_idx, BondOrder::Single));
    }

    #[test]
    fn build_adds_terminal_protein_bonds() {
        let mut residue = standard_residue("ALA", 1, ResiduePosition::CTerminal);
        let c_pos = residue.atom("C").unwrap().pos;
        let oxt_pos = c_pos + Vector3::new(1.24, 0.0, 0.0);
        let hxt_pos = oxt_pos + Vector3::new(0.96, 0.0, 0.0);

        residue.add_atom(Atom::new("OXT", Element::O, oxt_pos));
        residue.add_atom(Atom::new("HXT", Element::H, hxt_pos));

        let structure = structure_from_residues(vec![residue]);
        let topology = TopologyBuilder::new()
            .build(structure)
            .expect("build topology");

        let c_idx = global_atom_index(&topology, "A", 1, "C");
        let oxt_idx = global_atom_index(&topology, "A", 1, "OXT");
        let hxt_idx = global_atom_index(&topology, "A", 1, "HXT");

        assert!(has_bond(&topology, c_idx, oxt_idx, BondOrder::Single));
        assert!(has_bond(&topology, oxt_idx, hxt_idx, BondOrder::Single));
    }

    #[test]
    fn build_errors_when_required_standard_atom_missing() {
        let mut residue = standard_residue("GLY", 1, ResiduePosition::Internal);
        assert!(residue.remove_atom("CA").is_some(), "expected CA atom");

        let structure = structure_from_residues(vec![residue]);
        let err = TopologyBuilder::new().build(structure).unwrap_err();

        match err {
            Error::TopologyAtomMissing { atom_name, .. } => assert_eq!(atom_name, "CA"),
            other => panic!("unexpected error: {other:?}"),
        }
    }

    #[test]
    fn build_errors_when_hetero_template_missing() {
        let mut residue = Residue::new(1, None, "LIG", None, ResidueCategory::Hetero);
        residue.add_atom(Atom::new("C1", Element::C, Point::origin()));

        let structure = structure_from_residues(vec![residue]);
        let err = TopologyBuilder::new().build(structure).unwrap_err();

        match err {
            Error::MissingHeteroTemplate { res_name } => assert_eq!(res_name, "LIG"),
            other => panic!("unexpected error: {other:?}"),
        }
    }

    #[test]
    fn build_uses_hetero_template_for_hetero_residue() {
        let template = Template::new(
            "LIG",
            vec!["C1".into(), "O1".into(), "H1".into()],
            vec![
                ("C1".into(), "O1".into(), BondOrder::Double),
                ("C1".into(), "H1".into(), BondOrder::Single),
            ],
        );

        let mut residue = Residue::new(1, None, "LIG", None, ResidueCategory::Hetero);
        residue.add_atom(Atom::new("C1", Element::C, Point::new(0.0, 0.0, 0.0)));
        residue.add_atom(Atom::new("O1", Element::O, Point::new(1.2, 0.0, 0.0)));
        residue.add_atom(Atom::new("H1", Element::H, Point::new(-1.0, 0.0, 0.0)));

        let structure = structure_from_residues(vec![residue]);
        let builder = TopologyBuilder::new().add_hetero_template(template);
        let topology = builder.build(structure).expect("build topology");

        let c_idx = global_atom_index(&topology, "A", 1, "C1");
        let o_idx = global_atom_index(&topology, "A", 1, "O1");
        let h_idx = global_atom_index(&topology, "A", 1, "H1");

        assert!(has_bond(&topology, c_idx, o_idx, BondOrder::Double));
        assert!(has_bond(&topology, c_idx, h_idx, BondOrder::Single));
    }

    #[test]
    fn build_creates_disulfide_bond_when_sg_within_cutoff() {
        let residue1 = standard_residue("CYX", 1, ResiduePosition::Internal);
        let mut residue2 = standard_residue("CYX", 2, ResiduePosition::Internal);

        let sg1_pos = residue1.atom("SG").unwrap().pos;
        let sg2_pos = residue2.atom("SG").unwrap().pos;
        let target = sg1_pos + Vector3::new(2.0, 0.0, 0.0);
        translate_residue(&mut residue2, target - sg2_pos);

        let structure = structure_from_residues(vec![residue1, residue2]);
        let topology = TopologyBuilder::new()
            .build(structure)
            .expect("build topology");

        let sg1_idx = global_atom_index(&topology, "A", 1, "SG");
        let sg2_idx = global_atom_index(&topology, "A", 2, "SG");

        assert!(has_bond(&topology, sg1_idx, sg2_idx, BondOrder::Single));
    }
}
