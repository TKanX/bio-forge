//! Bond topology builder used to reconstruct intra- and inter-residue connectivity.
//!
//! The routines in this module consume a fully prepared [`Structure`] and emit a
//! [`Topology`] that encodes all covalent bonds. Internal templates are used for
//! standard residues, while callers can provide additional hetero templates.
//! Beyond template-driven intra-residue bonds, the builder also infers peptide,
//! nucleic-backbone, terminal, and disulfide bonds using geometric thresholds.

use crate::db;
use crate::model::{
    grid::Grid,
    structure::Structure,
    template::Template,
    topology::{Bond, Topology},
    types::{BondOrder, ResidueCategory, ResiduePosition},
};
use crate::ops::error::Error;
use crate::utils::parallel::*;
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
    /// # Arguments
    ///
    /// * `structure` - Structure for which to build the bond topology.
    ///
    /// # Returns
    ///
    /// A `Result` containing the built [`Topology`] or an [`Error`].
    ///
    /// # Errors
    ///
    /// Returns [`Error`] when a required template or atom is missing.
    pub fn build(self, structure: Structure) -> Result<Topology, Error> {
        let mut chain_offsets = Vec::with_capacity(structure.chain_count());
        let mut current_offset = 0;
        for chain in structure.iter_chains() {
            chain_offsets.push(current_offset);
            current_offset += chain.atom_count();
        }

        assert_eq!(chain_offsets.len(), structure.chain_count());

        let hetero_templates = &self.hetero_templates;
        let peptide_cutoff = self.peptide_bond_cutoff;
        let nucleic_cutoff = self.nucleic_bond_cutoff;
        let disulfide_cutoff = self.disulfide_bond_cutoff;

        let (mut bonds, sulfurs) = structure
            .par_chains()
            .zip(chain_offsets)
            .map(|(chain, chain_start_offset)| {
                let mut local_bonds = Vec::new();
                let mut local_sulfurs = Vec::new();
                let mut residue_offset = chain_start_offset;

                let residues: Vec<_> = chain.iter_residues().collect();

                for (i, residue) in residues.iter().enumerate() {
                    let atom_count = residue.atom_count();

                    Self::build_intra_residue_for_residue(
                        residue,
                        residue_offset,
                        hetero_templates,
                        &mut local_bonds,
                    )?;

                    if i < residues.len() - 1 {
                        let next_residue = residues[i + 1];
                        let next_offset = residue_offset + atom_count;

                        Self::build_backbone_bond(
                            residue,
                            residue_offset,
                            next_residue,
                            next_offset,
                            peptide_cutoff,
                            nucleic_cutoff,
                            &mut local_bonds,
                        );
                    }

                    match residue.name.as_str() {
                        "CYX" | "CYM" => {
                            if let Some(sg_idx) = residue.iter_atoms().position(|a| a.name == "SG")
                            {
                                let sg_pos = residue.atoms()[sg_idx].pos;
                                local_sulfurs.push((sg_pos, residue_offset + sg_idx));
                            }
                        }
                        _ => {}
                    }

                    residue_offset += atom_count;
                }

                Ok((local_bonds, local_sulfurs))
            })
            .try_reduce(
                || (Vec::new(), Vec::new()),
                |mut a, b| {
                    a.0.extend(b.0);
                    a.1.extend(b.1);
                    Ok(a)
                },
            )?;

        if !sulfurs.is_empty() {
            let grid = Grid::new(
                sulfurs.iter().map(|(p, i)| (*p, *i)),
                disulfide_cutoff + 0.5,
            );

            let disulfide_bonds: Vec<Bond> = sulfurs
                .par_iter()
                .flat_map(|(pos, idx1)| {
                    grid.neighbors(pos, disulfide_cutoff)
                        .exact()
                        .filter_map(|idx2| {
                            if *idx1 < *idx2 {
                                Some(Bond::new(*idx1, *idx2, BondOrder::Single))
                            } else {
                                None
                            }
                        })
                        .collect::<Vec<_>>()
                })
                .collect();

            bonds.extend(disulfide_bonds);
        }

        bonds.par_sort_unstable();
        bonds.dedup();

        Ok(Topology::new(structure, bonds))
    }

    /// Helper to generate intra-residue bonds for a single residue.
    fn build_intra_residue_for_residue(
        residue: &crate::model::residue::Residue,
        offset: usize,
        hetero_templates: &HashMap<String, Template>,
        bonds: &mut Vec<Bond>,
    ) -> Result<(), Error> {
        if residue.category == ResidueCategory::Ion {
            return Ok(());
        }

        if residue.category == ResidueCategory::Standard {
            let tmpl_name = &residue.name;
            let tmpl_view =
                db::get_template(tmpl_name).ok_or_else(|| Error::MissingInternalTemplate {
                    res_name: tmpl_name.to_string(),
                })?;

            for (a1_name, a2_name, order) in tmpl_view.bonds() {
                Self::try_add_bond(residue, offset, a1_name, a2_name, order, bonds)?;
            }

            Self::handle_terminal_intra_bonds(residue, offset, bonds)?;
        } else if residue.category == ResidueCategory::Hetero {
            let tmpl = hetero_templates.get(residue.name.as_str()).ok_or_else(|| {
                Error::MissingHeteroTemplate {
                    res_name: residue.name.to_string(),
                }
            })?;

            for (a1_name, a2_name, order) in tmpl.bonds() {
                Self::try_add_bond(residue, offset, a1_name, a2_name, *order, bonds)?;
            }
        }

        Ok(())
    }

    /// Helper to generate backbone bonds between two residues.
    fn build_backbone_bond(
        curr: &crate::model::residue::Residue,
        curr_offset: usize,
        next: &crate::model::residue::Residue,
        next_offset: usize,
        peptide_cutoff: f64,
        nucleic_cutoff: f64,
        bonds: &mut Vec<Bond>,
    ) {
        if curr.category != ResidueCategory::Standard || next.category != ResidueCategory::Standard
        {
            return;
        }

        if let (Some(std1), Some(std2)) = (curr.standard_name, next.standard_name) {
            if std1.is_protein() && std2.is_protein() {
                Self::connect_atoms_if_close(
                    AtomLocator::new(curr, curr_offset, "C"),
                    AtomLocator::new(next, next_offset, "N"),
                    peptide_cutoff,
                    BondOrder::Single,
                    bonds,
                );
            } else if std1.is_nucleic() && std2.is_nucleic() {
                Self::connect_atoms_if_close(
                    AtomLocator::new(curr, curr_offset, "O3'"),
                    AtomLocator::new(next, next_offset, "P"),
                    nucleic_cutoff,
                    BondOrder::Single,
                    bonds,
                );
            }
        }
    }

    /// Attempts to add a bond and reports informative errors when atoms are
    /// missing.
    fn try_add_bond(
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
            (None, _) if Self::is_optional_terminal_atom(residue, name1) => Ok(()),
            (_, None) if Self::is_optional_terminal_atom(residue, name2) => Ok(()),
            (None, _) => Err(Error::topology_atom_missing(
                &*residue.name,
                residue.id,
                name1,
            )),
            (_, None) => Err(Error::topology_atom_missing(
                &*residue.name,
                residue.id,
                name2,
            )),
        }
    }

    /// Returns whether the provided atom is optional for the residue position.
    fn is_optional_terminal_atom(
        residue: &crate::model::residue::Residue,
        atom_name: &str,
    ) -> bool {
        let is_protein = residue.standard_name.is_some_and(|std| std.is_protein());
        let is_nucleic = residue.standard_name.is_some_and(|std| std.is_nucleic());

        match residue.position {
            ResiduePosition::NTerminal if is_protein => atom_name == "H",
            ResiduePosition::CTerminal if is_protein => matches!(atom_name, "HXT" | "HOXT"),
            ResiduePosition::FivePrime if is_nucleic => {
                if residue.has_atom("P") {
                    matches!(atom_name, "OP3" | "HOP3" | "HOP2")
                } else {
                    matches!(
                        atom_name,
                        "P" | "OP1" | "OP2" | "OP3" | "HOP3" | "HOP2" | "HO5'"
                    )
                }
            }
            ResiduePosition::ThreePrime if is_nucleic => atom_name == "HO3'",
            _ => false,
        }
    }

    /// Adds missing bonds for termini (protein and nucleic acid) that are not
    /// explicitly listed in templates.
    fn handle_terminal_intra_bonds(
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
            if let (Some(p_idx), Some(op3_idx)) = (
                residue.iter_atoms().position(|a| a.name == "P"),
                residue.iter_atoms().position(|a| a.name == "OP3"),
            ) {
                bonds.push(Bond::new(
                    offset + p_idx,
                    offset + op3_idx,
                    BondOrder::Single,
                ));

                if let Some(hop3_idx) = residue.iter_atoms().position(|a| a.name == "HOP3") {
                    bonds.push(Bond::new(
                        offset + op3_idx,
                        offset + hop3_idx,
                        BondOrder::Single,
                    ));
                }
            }

            if let (Some(ho5_idx), Some(o5_idx)) = (
                residue.iter_atoms().position(|a| a.name == "HO5'"),
                residue.iter_atoms().position(|a| a.name == "O5'"),
            ) {
                bonds.push(Bond::new(
                    offset + ho5_idx,
                    offset + o5_idx,
                    BondOrder::Single,
                ));
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

    /// Adds a bond when the specified atoms are within the provided cutoff.
    fn connect_atoms_if_close(
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
    use std::collections::HashSet;

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

    fn five_prime_residue_with_phosphate_and_op3(id: i32) -> Residue {
        let mut residue = standard_residue("DA", id, ResiduePosition::FivePrime);
        let p_pos = residue.atom("P").unwrap().pos;
        let op1_pos = residue.atom("OP1").unwrap().pos;
        let op2_pos = residue.atom("OP2").unwrap().pos;
        let o5_pos = residue.atom("O5'").unwrap().pos;
        let centroid = (op1_pos.coords + op2_pos.coords + o5_pos.coords) / 3.0;
        let direction = (p_pos.coords - centroid).normalize();
        let op3_pos = p_pos + direction * 1.48;
        residue.add_atom(Atom::new("OP3", Element::O, op3_pos));
        residue
    }

    fn five_prime_residue_without_phosphate(id: i32) -> Residue {
        let template = db::get_template("DA").unwrap();
        let mut residue = Residue::new(
            id,
            None,
            "DA",
            Some(template.standard_name()),
            ResidueCategory::Standard,
        );
        residue.position = ResiduePosition::FivePrime;

        for (atom_name, element, point) in template.heavy_atoms() {
            if !matches!(atom_name, "P" | "OP1" | "OP2") {
                residue.add_atom(Atom::new(atom_name, element, point));
            }
        }

        for (atom_name, point, _) in template.hydrogens() {
            residue.add_atom(Atom::new(atom_name, Element::H, point));
        }

        let o5_pos = residue.atom("O5'").unwrap().pos;
        let c5_pos = residue.atom("C5'").unwrap().pos;
        let h_dir = (o5_pos - c5_pos).normalize();
        let ho5_pos = o5_pos + h_dir * 0.96;
        residue.add_atom(Atom::new("HO5'", Element::H, ho5_pos));

        residue
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

    #[test]
    fn build_avoids_duplicate_bonds_for_standard_residue() {
        let residue = standard_residue("ALA", 1, ResiduePosition::Internal);
        let structure = structure_from_residues(vec![residue]);
        let topology = TopologyBuilder::new()
            .build(structure)
            .expect("build topology");

        let mut seen = HashSet::new();
        for bond in topology.bonds() {
            let key = (bond.a1_idx, bond.a2_idx, bond.order);
            assert!(seen.insert(key), "duplicate bond detected: {key:?}");
        }
    }

    #[test]
    fn build_sets_expected_water_degree() {
        let residue = standard_residue("HOH", 1, ResiduePosition::Internal);
        let structure = structure_from_residues(vec![residue]);
        let topology = TopologyBuilder::new()
            .build(structure)
            .expect("build topology");

        let o_idx = global_atom_index(&topology, "A", 1, "O");
        let neighbors: Vec<_> = topology.neighbors_of(o_idx).collect();

        assert_eq!(neighbors.len(), 2, "water oxygen should have two neighbors");

        let mut uniq = HashSet::new();
        for idx in neighbors {
            assert!(uniq.insert(idx), "neighbor duplicated for water oxygen");
        }
    }

    #[test]
    fn build_creates_p_op3_bond_for_5prime_phosphate() {
        let residue = five_prime_residue_with_phosphate_and_op3(1);
        let structure = structure_from_residues(vec![residue]);
        let topology = TopologyBuilder::new()
            .build(structure)
            .expect("build topology");

        let p_idx = global_atom_index(&topology, "A", 1, "P");
        let op3_idx = global_atom_index(&topology, "A", 1, "OP3");

        assert!(has_bond(&topology, p_idx, op3_idx, BondOrder::Single));
    }

    #[test]
    fn build_creates_op3_hop3_bond_when_protonated() {
        let mut residue = five_prime_residue_with_phosphate_and_op3(1);
        let op3_pos = residue.atom("OP3").unwrap().pos;
        let p_pos = residue.atom("P").unwrap().pos;
        let h_dir = (op3_pos - p_pos).normalize();
        let hop3_pos = op3_pos + h_dir * 0.96;
        residue.add_atom(Atom::new("HOP3", Element::H, hop3_pos));

        let structure = structure_from_residues(vec![residue]);
        let topology = TopologyBuilder::new()
            .build(structure)
            .expect("build topology");

        let op3_idx = global_atom_index(&topology, "A", 1, "OP3");
        let hop3_idx = global_atom_index(&topology, "A", 1, "HOP3");

        assert!(has_bond(&topology, op3_idx, hop3_idx, BondOrder::Single));
    }

    #[test]
    fn build_creates_o5_ho5_bond_for_5prime_hydroxyl() {
        let residue = five_prime_residue_without_phosphate(1);
        let structure = structure_from_residues(vec![residue]);
        let topology = TopologyBuilder::new()
            .build(structure)
            .expect("build topology");

        let o5_idx = global_atom_index(&topology, "A", 1, "O5'");
        let ho5_idx = global_atom_index(&topology, "A", 1, "HO5'");

        assert!(has_bond(&topology, o5_idx, ho5_idx, BondOrder::Single));
    }

    #[test]
    fn optional_terminal_atoms_includes_nucleic_5prime_special_atoms() {
        let residue_with_p = five_prime_residue_with_phosphate_and_op3(1);

        assert!(TopologyBuilder::is_optional_terminal_atom(
            &residue_with_p,
            "OP3"
        ));
        assert!(TopologyBuilder::is_optional_terminal_atom(
            &residue_with_p,
            "HOP3"
        ));
        assert!(TopologyBuilder::is_optional_terminal_atom(
            &residue_with_p,
            "HOP2"
        ));
        assert!(!TopologyBuilder::is_optional_terminal_atom(
            &residue_with_p,
            "P"
        ));
        assert!(!TopologyBuilder::is_optional_terminal_atom(
            &residue_with_p,
            "OP1"
        ));

        let residue_no_p = five_prime_residue_without_phosphate(2);
        assert!(TopologyBuilder::is_optional_terminal_atom(
            &residue_no_p,
            "P"
        ));
        assert!(TopologyBuilder::is_optional_terminal_atom(
            &residue_no_p,
            "OP1"
        ));
        assert!(TopologyBuilder::is_optional_terminal_atom(
            &residue_no_p,
            "OP2"
        ));
        assert!(TopologyBuilder::is_optional_terminal_atom(
            &residue_no_p,
            "HO5'"
        ));
    }
}
