//! mmCIF writer utilities that encode structures and topologies into loop-based records.
//!
//! The serializer emits `data_` headers, optional `_cell.*` metadata, `_atom_site` rows with
//! consistent numbering, and `_struct_conn` loops reconstructed from [`Topology`] bonds so
//! downstream crystallography pipelines can round-trip `bio-forge` structures.

use crate::io::error::Error;
use crate::model::{
    atom::Atom, residue::Residue, structure::Structure, topology::Topology, types::BondOrder,
};
use std::collections::HashMap;
use std::io::Write;

/// Serializes a [`Structure`] into mmCIF format with optional cell metadata.
///
/// The writer emits a `data_` header, `_cell.*` entries when box vectors exist, and a
/// complete `_atom_site` loop that captures atom identities, entity IDs, and coordinates.
///
/// # Arguments
///
/// * `writer` - Destination that implements [`Write`].
/// * `structure` - Source structure providing atoms and unit-cell information.
///
/// # Returns
///
/// [`Ok`] on success or [`Error`] if any IO operation fails.
pub fn write_structure<W: Write>(writer: W, structure: &Structure) -> Result<(), Error> {
    let mut ctx = WriterContext::new(writer);

    ctx.write_header()?;

    ctx.write_cell(structure.box_vectors)?;

    ctx.write_entity_poly_seq(structure)?;

    ctx.write_atoms(structure)?;

    Ok(())
}

/// Serializes a [`Topology`] into mmCIF, including `_struct_conn` bond loops.
///
/// Coordinates are emitted via [`write_structure`] semantics; additionally, bonds are
/// converted into `_struct_conn` records with distance and order annotations.
///
/// # Arguments
///
/// * `writer` - Output sink implementing [`Write`].
/// * `topology` - Topology containing a structure and bond list to serialize.
///
/// # Returns
///
/// [`Ok`] when writing succeeds or [`Error`] if IO fails or bonds reference missing atoms.
pub fn write_topology<W: Write>(writer: W, topology: &Topology) -> Result<(), Error> {
    let mut ctx = WriterContext::new(writer);
    let structure = topology.structure();

    ctx.write_header()?;

    ctx.write_cell(structure.box_vectors)?;

    ctx.write_entity_poly_seq(structure)?;

    ctx.write_atoms(structure)?;

    ctx.write_connections(topology)?;

    Ok(())
}

/// Stateful helper that tracks atom numbering and writes mmCIF sections.
struct WriterContext<W> {
    writer: W,
    current_atom_id: usize,
    atom_index_to_id: HashMap<usize, usize>,
    residue_label_map: HashMap<(String, i32, Option<char>), String>,
}

impl<W: Write> WriterContext<W> {
    /// Creates a writer context with cleared indices and serial counters.
    ///
    /// # Arguments
    ///
    /// * `writer` - Sink receiving the emitted mmCIF text.
    fn new(writer: W) -> Self {
        Self {
            writer,
            current_atom_id: 1,
            atom_index_to_id: HashMap::new(),
            residue_label_map: HashMap::new(),
        }
    }

    /// Writes the `data_` block header and separating comment line.
    ///
    /// # Returns
    ///
    /// [`Ok`] when the header is flushed; [`Error`] if writing fails.
    fn write_header(&mut self) -> Result<(), Error> {
        writeln!(self.writer, "data_bio_forge_export")
            .and_then(|_| writeln!(self.writer, "#"))
            .map_err(|e| Error::from_io(e, None))
    }

    /// Emits `_cell.*` parameters derived from optional box vectors.
    ///
    /// # Arguments
    ///
    /// * `box_vectors` - Optional orthogonal matrix describing the unit cell.
    fn write_cell(&mut self, box_vectors: Option<[[f64; 3]; 3]>) -> Result<(), Error> {
        if let Some(vectors) = box_vectors {
            let v1 = nalgebra::Vector3::from(vectors[0]);
            let v2 = nalgebra::Vector3::from(vectors[1]);
            let v3 = nalgebra::Vector3::from(vectors[2]);

            let a = v1.norm();
            let b = v2.norm();
            let c = v3.norm();

            let alpha = v2.angle(&v3).to_degrees();
            let beta = v1.angle(&v3).to_degrees();
            let gamma = v1.angle(&v2).to_degrees();

            writeln!(self.writer, "_cell.entry_id           bio_forge_export")
                .map_err(|e| Error::from_io(e, None))?;
            writeln!(self.writer, "_cell.length_a           {:.3}", a)
                .map_err(|e| Error::from_io(e, None))?;
            writeln!(self.writer, "_cell.length_b           {:.3}", b)
                .map_err(|e| Error::from_io(e, None))?;
            writeln!(self.writer, "_cell.length_c           {:.3}", c)
                .map_err(|e| Error::from_io(e, None))?;
            writeln!(self.writer, "_cell.angle_alpha        {:.2}", alpha)
                .map_err(|e| Error::from_io(e, None))?;
            writeln!(self.writer, "_cell.angle_beta         {:.2}", beta)
                .map_err(|e| Error::from_io(e, None))?;
            writeln!(self.writer, "_cell.angle_gamma        {:.2}", gamma)
                .map_err(|e| Error::from_io(e, None))?;
            writeln!(self.writer, "_cell.Z_PDB              1")
                .map_err(|e| Error::from_io(e, None))?;
            writeln!(self.writer, "#").map_err(|e| Error::from_io(e, None))?;
        }
        Ok(())
    }

    /// Writes the `_entity_poly_seq` loop for polymer chains.
    ///
    /// # Arguments
    ///
    /// * `structure` - Structure whose polymer residues will be serialized.
    fn write_entity_poly_seq(&mut self, structure: &Structure) -> Result<(), Error> {
        let mut buffer = Vec::new();
        let mut next_entity_id = 1usize;
        let mut entity_ids: HashMap<smol_str::SmolStr, usize> = HashMap::new();

        for chain in structure.iter_chains() {
            let chain_id = chain.id.clone();
            let entity_id = *entity_ids.entry(chain_id).or_insert_with(|| {
                let val = next_entity_id;
                next_entity_id += 1;
                val
            });

            let polymer_residues: Vec<_> = chain
                .iter_residues()
                .filter(|r| {
                    r.standard_name
                        .is_some_and(|s| s.is_protein() || s.is_nucleic())
                })
                .collect();

            if !polymer_residues.is_empty() {
                for (i, residue) in polymer_residues.iter().enumerate() {
                    writeln!(buffer, "{} {} {} n", entity_id, i + 1, residue.name)
                        .map_err(|e| Error::from_io(e, None))?;
                }
            }
        }

        if !buffer.is_empty() {
            writeln!(self.writer, "loop_").map_err(|e| Error::from_io(e, None))?;
            writeln!(self.writer, "_entity_poly_seq.entity_id")
                .map_err(|e| Error::from_io(e, None))?;
            writeln!(self.writer, "_entity_poly_seq.num").map_err(|e| Error::from_io(e, None))?;
            writeln!(self.writer, "_entity_poly_seq.mon_id")
                .map_err(|e| Error::from_io(e, None))?;
            writeln!(self.writer, "_entity_poly_seq.hetero")
                .map_err(|e| Error::from_io(e, None))?;
            self.writer
                .write_all(&buffer)
                .map_err(|e| Error::from_io(e, None))?;
            writeln!(self.writer, "#").map_err(|e| Error::from_io(e, None))?;
        }

        Ok(())
    }

    /// Writes the `_atom_site` loop and assigns mmCIF atom IDs.
    ///
    /// # Arguments
    ///
    /// * `structure` - Structure whose atoms will be serialized.
    fn write_atoms(&mut self, structure: &Structure) -> Result<(), Error> {
        writeln!(self.writer, "loop_").map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_atom_site.group_PDB").map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_atom_site.id").map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_atom_site.type_symbol").map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_atom_site.label_atom_id").map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_atom_site.label_alt_id").map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_atom_site.label_comp_id").map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_atom_site.label_asym_id").map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_atom_site.label_entity_id").map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_atom_site.label_seq_id").map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_atom_site.pdbx_PDB_ins_code")
            .map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_atom_site.Cartn_x").map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_atom_site.Cartn_y").map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_atom_site.Cartn_z").map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_atom_site.occupancy").map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_atom_site.B_iso_or_equiv").map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_atom_site.auth_seq_id").map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_atom_site.auth_comp_id").map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_atom_site.auth_asym_id").map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_atom_site.auth_atom_id").map_err(|e| Error::from_io(e, None))?;

        self.atom_index_to_id.clear();
        self.residue_label_map.clear();
        let mut entity_ids: HashMap<smol_str::SmolStr, usize> = HashMap::new();
        let mut next_entity_id = 1usize;
        let mut global_atom_index = 0usize;

        for chain in structure.iter_chains() {
            let chain_id = chain.id.clone();
            let entity_id = *entity_ids.entry(chain_id.clone()).or_insert_with(|| {
                let val = next_entity_id;
                next_entity_id += 1;
                val
            });

            let mut polymer_seq_id = 0;

            for residue in chain.iter_residues() {
                let is_polymer = residue
                    .standard_name
                    .is_some_and(|s| s.is_protein() || s.is_nucleic());
                let label_seq_id = if is_polymer {
                    polymer_seq_id += 1;
                    polymer_seq_id.to_string()
                } else {
                    ".".to_string()
                };

                self.residue_label_map.insert(
                    (chain_id.to_string(), residue.id, residue.insertion_code),
                    label_seq_id.clone(),
                );

                for atom in residue.iter_atoms() {
                    let group_pdb = if is_polymer { "ATOM" } else { "HETATM" };

                    let atom_id = self.current_atom_id;

                    self.write_atom_record(
                        group_pdb,
                        atom,
                        residue,
                        &chain_id,
                        entity_id,
                        &label_seq_id,
                    )?;

                    self.atom_index_to_id.insert(global_atom_index, atom_id);
                    global_atom_index += 1;
                    self.current_atom_id += 1;
                }
            }
        }
        writeln!(self.writer, "#").map_err(|e| Error::from_io(e, None))?;
        Ok(())
    }

    /// Formats a single `_atom_site` row with quoted identifiers.
    ///
    /// # Arguments
    ///
    /// * `group_pdb` - Either `"ATOM"` or `"HETATM"`.
    /// * `atom` - Atom providing coordinates and element symbol.
    /// * `residue` - Residue metadata for labels and auth fields.
    /// * `chain_id` - Parent chain identifier string.
    /// * `entity_id` - Numeric entity identifier assigned to the chain.
    /// * `label_seq_id` - Sequential residue index for polymers, or `.` for others.
    fn write_atom_record(
        &mut self,
        group_pdb: &str,
        atom: &Atom,
        residue: &Residue,
        chain_id: &str,
        entity_id: usize,
        label_seq_id: &str,
    ) -> Result<(), Error> {
        let atom_id = self.current_atom_id;
        let type_symbol = atom.element.symbol();
        let label_atom_id = quote_string(&atom.name);
        let label_comp_id = quote_string(&residue.name);
        let label_asym_id = quote_string(chain_id);
        let ins_code = residue
            .insertion_code
            .map(|c| c.to_string())
            .unwrap_or_else(|| "?".to_string());

        let auth_seq_id = residue.id.to_string();
        let auth_comp_id = label_comp_id.clone();
        let auth_asym_id = label_asym_id.clone();
        let auth_atom_id = label_atom_id.clone();

        writeln!(
            self.writer,
            "{group_pdb} {atom_id} {type_symbol} {label_atom_id} . {label_comp_id} {label_asym_id} {entity_id} {label_seq_id} {ins_code} {x:.3} {y:.3} {z:.3} 1.00 0.00 {auth_seq_id} {auth_comp_id} {auth_asym_id} {auth_atom_id}",
            group_pdb = group_pdb,
            atom_id = atom_id,
            type_symbol = type_symbol,
            label_atom_id = label_atom_id,
            label_comp_id = label_comp_id,
            label_asym_id = label_asym_id,
            entity_id = entity_id,
            label_seq_id = label_seq_id,
            ins_code = ins_code,
            x = atom.pos.x,
            y = atom.pos.y,
            z = atom.pos.z,
            auth_seq_id = auth_seq_id,
            auth_comp_id = auth_comp_id,
            auth_asym_id = auth_asym_id,
            auth_atom_id = auth_atom_id
        )
        .map_err(|e| Error::from_io(e, None))
    }

    /// Serializes topology bonds into `_struct_conn` records with distances.
    ///
    /// # Arguments
    ///
    /// * `topology` - Topology whose bonds will be emitted.
    ///
    /// # Returns
    ///
    /// [`Ok`] if all bonds were written or [`Error::InconsistentData`] when atom indices are
    /// missing because coordinates were not emitted beforehand.
    fn write_connections(&mut self, topology: &Topology) -> Result<(), Error> {
        if topology.bond_count() == 0 {
            return Ok(());
        }

        writeln!(self.writer, "loop_").map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_struct_conn.id").map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_struct_conn.conn_type_id").map_err(|e| Error::from_io(e, None))?;

        writeln!(self.writer, "_struct_conn.ptnr1_label_atom_id")
            .map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_struct_conn.ptnr1_label_alt_id")
            .map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_struct_conn.ptnr1_label_comp_id")
            .map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_struct_conn.ptnr1_label_asym_id")
            .map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_struct_conn.ptnr1_label_seq_id")
            .map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_struct_conn.ptnr1_PDB_ins_code")
            .map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_struct_conn.ptnr1_symmetry")
            .map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_struct_conn.ptnr1_auth_asym_id")
            .map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_struct_conn.ptnr1_auth_comp_id")
            .map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_struct_conn.ptnr1_auth_seq_id")
            .map_err(|e| Error::from_io(e, None))?;

        writeln!(self.writer, "_struct_conn.ptnr2_label_atom_id")
            .map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_struct_conn.ptnr2_label_alt_id")
            .map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_struct_conn.ptnr2_label_comp_id")
            .map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_struct_conn.ptnr2_label_asym_id")
            .map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_struct_conn.ptnr2_label_seq_id")
            .map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_struct_conn.ptnr2_PDB_ins_code")
            .map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_struct_conn.ptnr2_symmetry")
            .map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_struct_conn.ptnr2_auth_asym_id")
            .map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_struct_conn.ptnr2_auth_comp_id")
            .map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_struct_conn.ptnr2_auth_seq_id")
            .map_err(|e| Error::from_io(e, None))?;

        writeln!(self.writer, "_struct_conn.pdbx_dist_value")
            .map_err(|e| Error::from_io(e, None))?;
        writeln!(self.writer, "_struct_conn.pdbx_value_order")
            .map_err(|e| Error::from_io(e, None))?;

        let atom_lookup: Vec<(&crate::model::chain::Chain, &Residue, &Atom)> =
            topology.structure().iter_atoms_with_context().collect();

        for (conn_idx, bond) in topology.bonds().iter().enumerate() {
            let _ = self.atom_index_to_id.get(&bond.a1_idx).ok_or_else(|| {
                Error::inconsistent_data(
                    "mmCIF",
                    None,
                    format!(
                        "bond references atom index {} that was not written",
                        bond.a1_idx
                    ),
                )
            })?;
            let _ = self.atom_index_to_id.get(&bond.a2_idx).ok_or_else(|| {
                Error::inconsistent_data(
                    "mmCIF",
                    None,
                    format!(
                        "bond references atom index {} that was not written",
                        bond.a2_idx
                    ),
                )
            })?;

            let (chain1, res1, atom1) = atom_lookup[bond.a1_idx];
            let (chain2, res2, atom2) = atom_lookup[bond.a2_idx];

            let label_seq_1 = self
                .residue_label_map
                .get(&(chain1.id.to_string(), res1.id, res1.insertion_code))
                .map(|s| s.as_str())
                .unwrap_or("?");
            let label_seq_2 = self
                .residue_label_map
                .get(&(chain2.id.to_string(), res2.id, res2.insertion_code))
                .map(|s| s.as_str())
                .unwrap_or("?");

            let conn_id = format!("conn_{:04}", conn_idx + 1);
            let conn_type_id = "covale";
            let symmetry = "1_555";
            let order_str = match bond.order {
                BondOrder::Single => "SING",
                BondOrder::Double => "DOUB",
                BondOrder::Triple => "TRIP",
                BondOrder::Aromatic => "AROM",
            };
            let dist = atom1.distance(atom2);

            let ins1 = res1
                .insertion_code
                .map(|c| c.to_string())
                .unwrap_or_else(|| "?".to_string());
            let ins2 = res2
                .insertion_code
                .map(|c| c.to_string())
                .unwrap_or_else(|| "?".to_string());

            writeln!(
                self.writer,
                "{conn_id} {conn_type_id} {pt1_atom} . {pt1_res} {pt1_asym} {pt1_seq} {pt1_ins} {symmetry} {pt1_auth_asym} {pt1_auth_res} {pt1_auth_seq} {pt2_atom} . {pt2_res} {pt2_asym} {pt2_seq} {pt2_ins} {symmetry} {pt2_auth_asym} {pt2_auth_res} {pt2_auth_seq} {dist:.3} {order_str}",
                conn_id = conn_id,
                conn_type_id = conn_type_id,
                pt1_atom = quote_string(&atom1.name),
                pt1_res = quote_string(&res1.name),
                pt1_asym = quote_string(&chain1.id),
                pt1_seq = label_seq_1,
                pt1_ins = ins1,
                symmetry = symmetry,
                pt1_auth_asym = quote_string(&chain1.id),
                pt1_auth_res = quote_string(&res1.name),
                pt1_auth_seq = res1.id,
                pt2_atom = quote_string(&atom2.name),
                pt2_res = quote_string(&res2.name),
                pt2_asym = quote_string(&chain2.id),
                pt2_seq = label_seq_2,
                pt2_ins = ins2,
                pt2_auth_asym = quote_string(&chain2.id),
                pt2_auth_res = quote_string(&res2.name),
                pt2_auth_seq = res2.id,
                dist = dist,
                order_str = order_str
            )
            .map_err(|e| Error::from_io(e, None))?;
        }
        writeln!(self.writer, "#").map_err(|e| Error::from_io(e, None))?;

        Ok(())
    }
}

/// Wraps strings containing whitespace or quotes with CIF-safe quoting.
///
/// Empty strings become `?`, single quotes trigger double-quote wrapping, and all other
/// cases fall back to single quotes.
///
/// # Arguments
///
/// * `s` - Raw string to sanitize.
///
/// # Returns
///
/// Quoted string acceptable for mmCIF fields.
fn quote_string(s: &str) -> String {
    if s.is_empty() {
        return "?".to_string();
    }
    if !s.contains(char::is_whitespace) && !s.contains('\'') && !s.contains('"') {
        return s.to_string();
    }
    if s.contains('\'') && !s.contains('"') {
        return format!("\"{}\"", s);
    }
    format!("'{}'", s)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::atom::Atom;
    use crate::model::chain::Chain;
    use crate::model::residue::Residue;
    use crate::model::topology::{Bond, Topology};
    use crate::model::types::{BondOrder, Element, Point, ResidueCategory, StandardResidue};

    fn create_atom(name: &str, element: Element) -> Atom {
        Atom::new(name, element, Point::new(0.0, 0.0, 0.0))
    }

    fn create_residue(
        id: i32,
        name: &str,
        std: Option<StandardResidue>,
        cat: ResidueCategory,
    ) -> Residue {
        Residue::new(id, None, name, std, cat)
    }

    fn build_test_structure() -> Structure {
        let mut structure = Structure::new();
        structure.box_vectors = Some([[10.0, 0.0, 0.0], [0.0, 11.0, 0.0], [0.0, 0.0, 12.0]]);

        let mut chain = Chain::new("A");

        let mut gly = create_residue(
            1,
            "GLY",
            Some(StandardResidue::GLY),
            ResidueCategory::Standard,
        );
        gly.add_atom(Atom::new("N", Element::N, Point::new(0.0, 0.0, 0.0)));
        gly.add_atom(Atom::new("CA", Element::C, Point::new(1.0, 0.0, 0.0)));

        let mut lig = create_residue(2, "LIG", None, ResidueCategory::Hetero);
        lig.add_atom(Atom::new("C1", Element::C, Point::new(4.0, 5.0, 6.0)));

        chain.add_residue(gly);
        chain.add_residue(lig);
        structure.add_chain(chain);

        structure
    }

    #[test]
    fn write_header_and_cell_emits_correct_metadata() {
        let structure = build_test_structure();
        let mut buffer = Vec::new();

        write_structure(&mut buffer, &structure).expect("structure write failed");
        let output = String::from_utf8(buffer).expect("invalid UTF-8");

        assert!(output.contains("data_bio_forge_export"));
        assert!(output.contains("_cell.length_a           10.000"));
        assert!(output.contains("_cell.length_b           11.000"));
        assert!(output.contains("_cell.length_c           12.000"));
        assert!(output.contains("_cell.angle_alpha        90.00"));
    }

    #[test]
    fn write_entity_poly_seq_generates_loop_for_polymers_only() {
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");

        let ala = create_residue(
            1,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        let val = create_residue(
            2,
            "VAL",
            Some(StandardResidue::VAL),
            ResidueCategory::Standard,
        );
        let hoh = create_residue(
            3,
            "HOH",
            Some(StandardResidue::HOH),
            ResidueCategory::Hetero,
        );

        chain.add_residue(ala);
        chain.add_residue(val);
        chain.add_residue(hoh);
        structure.add_chain(chain);

        let mut buffer = Vec::new();
        write_structure(&mut buffer, &structure).expect("structure write failed");
        let output = String::from_utf8(buffer).expect("invalid UTF-8");

        assert!(output.contains("_entity_poly_seq.entity_id"));
        assert!(output.contains("_entity_poly_seq.mon_id"));

        assert!(output.contains("1 1 ALA n"));
        assert!(output.contains("1 2 VAL n"));
        assert!(!output.contains("1 3 HOH"));
    }

    #[test]
    fn label_seq_id_is_sequential_and_skips_non_polymers() {
        let mut structure = Structure::new();
        let mut chain = Chain::new("A");

        let mut r1 = create_residue(
            10,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        r1.add_atom(create_atom("CA", Element::C));

        let mut r2 = create_residue(20, "LIG", None, ResidueCategory::Hetero);
        r2.add_atom(create_atom("X", Element::C));

        let mut r3 = create_residue(
            30,
            "VAL",
            Some(StandardResidue::VAL),
            ResidueCategory::Standard,
        );
        r3.add_atom(create_atom("CA", Element::C));

        chain.add_residue(r1);
        chain.add_residue(r2);
        chain.add_residue(r3);
        structure.add_chain(chain);

        let mut buffer = Vec::new();
        write_structure(&mut buffer, &structure).expect("structure write failed");
        let output = String::from_utf8(buffer).expect("invalid UTF-8");

        let lines: Vec<&str> = output.lines().collect();
        let atom_lines: Vec<&str> = lines
            .iter()
            .filter(|l| l.starts_with("ATOM") || l.starts_with("HETATM"))
            .cloned()
            .collect();

        assert_eq!(atom_lines.len(), 3);

        let ala_parts: Vec<&str> = atom_lines[0].split_whitespace().collect();
        assert_eq!(ala_parts[5], "ALA");
        assert_eq!(ala_parts[8], "1");
        assert_eq!(ala_parts[15], "10");

        let lig_parts: Vec<&str> = atom_lines[1].split_whitespace().collect();
        assert_eq!(lig_parts[5], "LIG");
        assert_eq!(lig_parts[8], ".");
        assert_eq!(lig_parts[15], "20");

        let val_parts: Vec<&str> = atom_lines[2].split_whitespace().collect();
        assert_eq!(val_parts[5], "VAL");
        assert_eq!(val_parts[8], "2");
        assert_eq!(val_parts[15], "30");
    }

    #[test]
    fn multiple_chains_are_assigned_distinct_entity_ids() {
        let mut structure = Structure::new();

        let mut chain_a = Chain::new("A");
        let mut res_a = create_residue(
            1,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        res_a.add_atom(create_atom("CA", Element::C));
        chain_a.add_residue(res_a);
        structure.add_chain(chain_a);

        let mut chain_b = Chain::new("B");
        let mut res_b = create_residue(
            1,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        res_b.add_atom(create_atom("CA", Element::C));
        chain_b.add_residue(res_b);
        structure.add_chain(chain_b);

        let mut buffer = Vec::new();
        write_structure(&mut buffer, &structure).expect("structure write failed");
        let output = String::from_utf8(buffer).expect("invalid UTF-8");

        let atom_lines: Vec<&str> = output.lines().filter(|l| l.starts_with("ATOM")).collect();

        let parts_a: Vec<&str> = atom_lines[0].split_whitespace().collect();
        assert_eq!(parts_a[6], "A");
        assert_eq!(parts_a[6], "A");
        assert_eq!(parts_a[7], "1");

        let parts_b: Vec<&str> = atom_lines[1].split_whitespace().collect();
        assert_eq!(parts_b[6], "B");
        assert_eq!(parts_b[7], "2");
    }

    #[test]
    fn write_topology_emits_struct_conn_records() {
        let structure = build_test_structure();
        let topology = Topology::new(structure.clone(), vec![Bond::new(0, 1, BondOrder::Single)]);

        let mut buffer = Vec::new();
        write_topology(&mut buffer, &topology).expect("topology write failed");
        let output = String::from_utf8(buffer).expect("invalid UTF-8");

        assert!(output.contains("_struct_conn.id"));

        let conn_line = output
            .lines()
            .find(|l| l.starts_with("conn_0001"))
            .expect("connection line found");
        let tokens: Vec<&str> = conn_line.split_whitespace().collect();

        assert_eq!(tokens[0], "conn_0001");
        assert_eq!(tokens[2], "N");
        assert_eq!(tokens[4], "GLY");
        assert_eq!(tokens[6], "1");
        assert_eq!(tokens[11], "1");
        assert_eq!(tokens[12], "CA");
        assert_eq!(tokens[22], "1.000");
        assert_eq!(tokens[23], "SING");
    }

    #[test]
    fn write_connections_returns_error_when_atom_missing() {
        let structure = build_test_structure();
        let topology = Topology::new(structure.clone(), vec![Bond::new(0, 1, BondOrder::Single)]);

        let mut ctx = WriterContext::new(Vec::new());

        let err = ctx.write_connections(&topology).expect_err("should fail");
        match err {
            Error::InconsistentData { details, .. } => {
                assert!(details.contains("bond references atom index"));
            }
            _ => panic!("wrong error type"),
        }
    }
}
