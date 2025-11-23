use crate::io::error::Error;
use crate::model::{
    atom::Atom, residue::Residue, structure::Structure, topology::Topology, types::ResidueCategory,
};
use std::collections::HashMap;
use std::io::Write;

pub fn write_structure<W: Write>(writer: W, structure: &Structure) -> Result<(), Error> {
    let mut ctx = WriterContext::new(writer);

    ctx.write_cryst1(structure.box_vectors)?;

    ctx.write_atoms(structure)?;

    ctx.write_end()?;

    Ok(())
}

pub fn write_topology<W: Write>(writer: W, topology: &Topology) -> Result<(), Error> {
    let mut ctx = WriterContext::new(writer);
    let structure = topology.structure();

    ctx.write_cryst1(structure.box_vectors)?;

    ctx.write_atoms(structure)?;

    ctx.write_connects(topology)?;

    ctx.write_end()?;

    Ok(())
}

struct WriterContext<W> {
    writer: W,
    current_serial: usize,
    atom_index_to_serial: HashMap<usize, usize>,
}

impl<W: Write> WriterContext<W> {
    fn new(writer: W) -> Self {
        Self {
            writer,
            current_serial: 1,
            atom_index_to_serial: HashMap::new(),
        }
    }

    fn write_cryst1(&mut self, box_vectors: Option<[[f64; 3]; 3]>) -> Result<(), Error> {
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

            writeln!(
                self.writer,
                "CRYST1{:9.3}{:9.3}{:9.3}{:7.2}{:7.2}{:7.2} P 1           1",
                a, b, c, alpha, beta, gamma
            )
            .map_err(|e| Error::from_io(e, None))?;
        }
        Ok(())
    }

    fn write_atoms(&mut self, structure: &Structure) -> Result<(), Error> {
        let mut global_idx = 0;

        for chain in structure.iter_chains() {
            for residue in chain.iter_residues() {
                for atom in residue.iter_atoms() {
                    let record_type = match residue.standard_name {
                        Some(std) if std.is_protein() || std.is_nucleic() => "ATOM  ",
                        _ => "HETATM",
                    };

                    let serial = self.current_serial;

                    self.atom_index_to_serial.insert(global_idx, serial);

                    self.write_atom_record(record_type, serial, atom, residue, &chain.id)?;

                    self.current_serial += 1;
                    global_idx += 1;
                }
            }

            if let Some(last_standard) = chain
                .iter_residues()
                .rev()
                .find(|res| res.category == ResidueCategory::Standard)
            {
                let serial = self.current_serial;
                self.write_ter_record(serial, last_standard, &chain.id)?;
                self.current_serial += 1;
            }
        }
        Ok(())
    }

    fn write_atom_record(
        &mut self,
        record_type: &str,
        serial: usize,
        atom: &Atom,
        residue: &Residue,
        chain_id: &str,
    ) -> Result<(), Error> {
        let atom_name = if atom.name.len() >= 4 {
            format!("{:<4}", &atom.name[0..4])
        } else {
            format!(" {:<3}", atom.name)
        };

        let res_name = if residue.name.len() > 3 {
            &residue.name[0..3]
        } else {
            &residue.name
        };

        let element_str = format!("{:>2}", atom.element.symbol().to_uppercase());

        writeln!(
            self.writer,
            "{:6}{:5} {:4}{:1}{:3} {:1}{:4}{:1}   {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}          {:2}",
            record_type,
            serial % 100000,
            atom_name,
            ' ',
            res_name,
            chain_id.chars().next().unwrap_or(' '),
            residue.id % 10000,
            residue.insertion_code.unwrap_or(' '),
            atom.pos.x,
            atom.pos.y,
            atom.pos.z,
            1.00,
            0.00,
            element_str
        )
        .map_err(|e| Error::from_io(e, None))
    }

    fn write_ter_record(
        &mut self,
        serial: usize,
        residue: &Residue,
        chain_id: &str,
    ) -> Result<(), Error> {
        let res_name = if residue.name.len() > 3 {
            &residue.name[0..3]
        } else {
            &residue.name
        };

        writeln!(
            self.writer,
            "TER   {:5}      {:3} {:1}{:4}{:1}",
            serial % 100000,
            res_name,
            chain_id.chars().next().unwrap_or(' '),
            residue.id % 10000,
            residue.insertion_code.unwrap_or(' ')
        )
        .map_err(|e| Error::from_io(e, None))
    }

    fn write_connects(&mut self, topology: &Topology) -> Result<(), Error> {
        let mut adjacency: HashMap<usize, Vec<usize>> = HashMap::new();

        for bond in topology.bonds() {
            let s1 = *self.atom_index_to_serial.get(&bond.a1_idx).ok_or_else(|| {
                Error::inconsistent_data(
                    "PDB",
                    None,
                    format!(
                        "bond references atom index {} that was not written",
                        bond.a1_idx
                    ),
                )
            })?;
            let s2 = *self.atom_index_to_serial.get(&bond.a2_idx).ok_or_else(|| {
                Error::inconsistent_data(
                    "PDB",
                    None,
                    format!(
                        "bond references atom index {} that was not written",
                        bond.a2_idx
                    ),
                )
            })?;

            adjacency.entry(s1).or_default().push(s2);
            adjacency.entry(s2).or_default().push(s1);
        }

        let mut serials: Vec<_> = adjacency.keys().copied().collect();
        serials.sort();

        for src_serial in serials {
            let targets = adjacency.get(&src_serial).unwrap();
            let mut targets = targets.clone();
            targets.sort();
            targets.dedup();

            for chunk in targets.chunks(4) {
                write!(self.writer, "CONECT{:5}", src_serial)
                    .map_err(|e| Error::from_io(e, None))?;
                for target in chunk {
                    write!(self.writer, "{:5}", target).map_err(|e| Error::from_io(e, None))?;
                }
                writeln!(self.writer).map_err(|e| Error::from_io(e, None))?;
            }
        }

        Ok(())
    }

    fn write_end(&mut self) -> Result<(), Error> {
        writeln!(self.writer, "END   ").map_err(|e| Error::from_io(e, None))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::atom::Atom;
    use crate::model::chain::Chain;
    use crate::model::residue::Residue;
    use crate::model::topology::{Bond, Topology};
    use crate::model::types::{BondOrder, Element, Point, ResidueCategory, StandardResidue};

    fn assert_cryst1_line(line: &str, params: (f64, f64, f64, f64, f64, f64)) {
        assert!(line.starts_with("CRYST1"));
        let (a, b, c, alpha, beta, gamma) = params;
        assert!((parse_float(&line[6..15]) - a).abs() < 1e-3);
        assert!((parse_float(&line[15..24]) - b).abs() < 1e-3);
        assert!((parse_float(&line[24..33]) - c).abs() < 1e-3);
        assert!((parse_float(&line[33..40]) - alpha).abs() < 1e-2);
        assert!((parse_float(&line[40..47]) - beta).abs() < 1e-2);
        assert!((parse_float(&line[47..54]) - gamma).abs() < 1e-2);
    }

    fn assert_atom_line(
        line: &str,
        record: &str,
        serial: usize,
        atom_name: &str,
        res_name: &str,
        chain_id: char,
        res_seq: i32,
        insertion_code: char,
        coords: (f64, f64, f64),
        element: &str,
    ) {
        assert!(line.len() >= 78, "line too short: {line}");
        assert_eq!(&line[0..6], record);
        assert_eq!(line[6..11].trim(), serial.to_string());
        assert_eq!(line[12..16].trim(), atom_name);
        assert_eq!(line[17..20].trim(), res_name);
        assert_eq!(line.chars().nth(21).unwrap(), chain_id);
        assert_eq!(line[22..26].trim(), res_seq.to_string());
        assert_eq!(line.chars().nth(26).unwrap(), insertion_code);
        assert!((parse_float(&line[30..38]) - coords.0).abs() < 1e-3);
        assert!((parse_float(&line[38..46]) - coords.1).abs() < 1e-3);
        assert!((parse_float(&line[46..54]) - coords.2).abs() < 1e-3);
        assert_eq!(line[76..78].trim(), element);
    }

    fn assert_ter_line(
        line: &str,
        serial: usize,
        res_name: &str,
        chain_id: char,
        res_seq: i32,
        insertion_code: char,
    ) {
        assert!(line.starts_with("TER   "));
        assert_eq!(line[6..11].trim(), serial.to_string());
        assert_eq!(line[17..20].trim(), res_name);
        assert_eq!(line.chars().nth(21).unwrap(), chain_id);
        assert_eq!(line[22..26].trim(), res_seq.to_string());
        assert_eq!(line.chars().nth(26).unwrap(), insertion_code);
    }

    fn assert_conect_line(line: &str, source: usize, targets: &[usize]) {
        assert!(line.starts_with("CONECT"));
        let tokens: Vec<_> = line.split_whitespace().collect();
        assert_eq!(tokens[0], "CONECT");
        assert_eq!(tokens[1].parse::<usize>().unwrap(), source);
        let parsed_targets: Vec<_> = tokens[2..]
            .iter()
            .map(|tok| tok.parse::<usize>().unwrap())
            .collect();
        assert_eq!(parsed_targets, targets);
    }

    fn parse_float(slice: &str) -> f64 {
        slice.trim().parse::<f64>().expect("valid float")
    }

    #[test]
    fn write_structure_emits_cryst1_atoms_ter_and_end() {
        let mut structure = Structure::new();
        structure.box_vectors = Some([[10.0, 0.0, 0.0], [0.0, 11.0, 0.0], [0.0, 0.0, 12.0]]);

        let mut chain = Chain::new("A");

        let mut gly = Residue::new(
            1,
            None,
            "GLY",
            Some(StandardResidue::GLY),
            ResidueCategory::Standard,
        );
        gly.add_atom(Atom::new("N", Element::N, Point::new(1.0, 2.0, 3.0)));
        gly.add_atom(Atom::new("CA", Element::C, Point::new(1.5, 2.5, 3.5)));

        let mut lig = Residue::new(2, None, "LIG", None, ResidueCategory::Hetero);
        lig.add_atom(Atom::new("C1", Element::C, Point::new(4.0, 5.0, 6.0)));

        chain.add_residue(gly);
        chain.add_residue(lig);
        structure.add_chain(chain);

        let mut buffer = Vec::new();
        write_structure(&mut buffer, &structure).expect("writer should succeed");

        let output = String::from_utf8(buffer).expect("valid UTF-8");
        let lines: Vec<&str> = output.lines().collect();
        assert_eq!(lines.len(), 6, "unexpected number of lines: {lines:?}");

        assert_cryst1_line(lines[0], (10.0, 11.0, 12.0, 90.0, 90.0, 90.0));
        assert_atom_line(
            lines[1],
            "ATOM  ",
            1,
            "N",
            "GLY",
            'A',
            1,
            ' ',
            (1.0, 2.0, 3.0),
            "N",
        );
        assert_atom_line(
            lines[2],
            "ATOM  ",
            2,
            "CA",
            "GLY",
            'A',
            1,
            ' ',
            (1.5, 2.5, 3.5),
            "C",
        );
        assert_atom_line(
            lines[3],
            "HETATM",
            3,
            "C1",
            "LIG",
            'A',
            2,
            ' ',
            (4.0, 5.0, 6.0),
            "C",
        );
        assert_ter_line(lines[4], 4, "GLY", 'A', 1, ' ');
        assert_eq!(lines[5], "END   ");
    }

    #[test]
    fn write_structure_without_box_starts_with_atom_records() {
        let mut structure = Structure::new();
        let mut chain = Chain::new("B");

        let mut ser = Residue::new(
            7,
            Some('A'),
            "SER",
            Some(StandardResidue::SER),
            ResidueCategory::Standard,
        );
        ser.add_atom(Atom::new("OG", Element::O, Point::new(-1.0, 0.5, 2.0)));
        chain.add_residue(ser);
        structure.add_chain(chain);

        let mut buffer = Vec::new();
        write_structure(&mut buffer, &structure).expect("writer should succeed");

        let output = String::from_utf8(buffer).expect("valid UTF-8");
        let mut lines = output.lines();

        let first_line = lines.next().expect("at least one line");
        assert!(first_line.starts_with("ATOM"));
        assert!(lines.any(|line| line == "END   "));
    }

    #[test]
    fn non_polymer_standard_residue_uses_hetatm_record() {
        let mut structure = Structure::new();
        let mut chain = Chain::new("W");

        let mut water = Residue::new(
            42,
            None,
            "HOH",
            Some(StandardResidue::HOH),
            ResidueCategory::Standard,
        );
        water.add_atom(Atom::new("O", Element::O, Point::new(0.0, 0.0, 0.0)));
        chain.add_residue(water);
        structure.add_chain(chain);

        let mut buffer = Vec::new();
        write_structure(&mut buffer, &structure).expect("writer should succeed");

        let output = String::from_utf8(buffer).expect("valid UTF-8");
        let first_line = output.lines().next().expect("at least one line");
        assert!(
            first_line.starts_with("HETATM"),
            "line should be HETATM: {first_line}"
        );
    }

    #[test]
    fn write_topology_emits_conect_records() {
        let mut structure = Structure::new();
        let mut chain = Chain::new("C");
        let mut ala = Residue::new(
            10,
            None,
            "ALA",
            Some(StandardResidue::ALA),
            ResidueCategory::Standard,
        );
        ala.add_atom(Atom::new("N", Element::N, Point::new(0.0, 0.0, 0.0)));
        ala.add_atom(Atom::new("CA", Element::C, Point::new(1.0, 0.0, 0.0)));
        chain.add_residue(ala);
        structure.add_chain(chain);

        let topology = Topology::new(structure.clone(), vec![Bond::new(0, 1, BondOrder::Single)]);

        let mut buffer = Vec::new();
        write_topology(&mut buffer, &topology).expect("topology writer succeeds");

        let output = String::from_utf8(buffer).expect("valid UTF-8");
        let conect_lines: Vec<&str> = output
            .lines()
            .filter(|line| line.starts_with("CONECT"))
            .collect();

        assert_eq!(conect_lines.len(), 2);
        assert_conect_line(conect_lines[0], 1, &[2]);
        assert_conect_line(conect_lines[1], 2, &[1]);
    }

    #[test]
    fn write_connects_returns_error_when_serial_missing() {
        let mut structure = Structure::new();
        let mut chain = Chain::new("D");
        let mut gly = Residue::new(
            5,
            None,
            "GLY",
            Some(StandardResidue::GLY),
            ResidueCategory::Standard,
        );
        gly.add_atom(Atom::new("N", Element::N, Point::new(0.0, 0.0, 0.0)));
        gly.add_atom(Atom::new("CA", Element::C, Point::new(1.0, 1.0, 0.0)));
        chain.add_residue(gly);
        structure.add_chain(chain);

        let topology = Topology::new(structure.clone(), vec![Bond::new(0, 1, BondOrder::Single)]);

        let mut ctx = WriterContext::new(Vec::new());
        ctx.write_atoms(&structure).expect("atoms should write");
        ctx.atom_index_to_serial.clear();

        let err = ctx
            .write_connects(&topology)
            .expect_err("missing serial map should error");

        match err {
            Error::InconsistentData { details, .. } => {
                assert!(details.contains("bond references atom index"));
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }
}
