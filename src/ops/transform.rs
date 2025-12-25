//! Geometric transformations for molecular structures.
//!
//! This module provides utilities for translating, centering, and rotating structures.

use crate::model::structure::Structure;
use crate::model::types::Point;
use crate::utils::parallel::*;
use nalgebra::{Rotation3, Vector3};

/// Collection of geometric transformation operations for structures.
///
/// The `Transform` type groups static methods that mutate structure coordinates in place.
/// These operations include translations, centering (by geometry or mass), and rotations
/// about principal axes or via Euler angles.
pub struct Transform;

impl Transform {
    /// Translates all atoms by the specified displacement vector.
    ///
    /// # Arguments
    ///
    /// * `structure` - Mutable structure whose atoms will be displaced.
    /// * `x` - Translation along the x-axis in ångströms.
    /// * `y` - Translation along the y-axis in ångströms.
    /// * `z` - Translation along the z-axis in ångströms.
    pub fn translate(structure: &mut Structure, x: f64, y: f64, z: f64) {
        let translation = Vector3::new(x, y, z);
        structure.par_residues_mut().for_each(|residue| {
            for atom in residue.iter_atoms_mut() {
                atom.translate_by(&translation);
            }
        });
    }

    /// Centers the structure's geometric centroid at the target point.
    ///
    /// When `target` is `None`, the structure is centered at the origin.
    ///
    /// # Arguments
    ///
    /// * `structure` - Mutable structure to be centered.
    /// * `target` - Optional target point; defaults to the origin.
    pub fn center_geometry(structure: &mut Structure, target: Option<Point>) {
        let current_center = structure.geometric_center();
        let target_point = target.unwrap_or(Point::origin());
        let translation = target_point - current_center;

        structure.par_residues_mut().for_each(|residue| {
            for atom in residue.iter_atoms_mut() {
                atom.translate_by(&translation);
            }
        });
    }

    /// Centers the structure's center of mass at the target point.
    ///
    /// Mass weighting uses atomic masses from element definitions. When `target` is
    /// `None`, the structure is centered at the origin.
    ///
    /// # Arguments
    ///
    /// * `structure` - Mutable structure to be centered.
    /// * `target` - Optional target point; defaults to the origin.
    pub fn center_mass(structure: &mut Structure, target: Option<Point>) {
        let current_com = structure.center_of_mass();
        let target_point = target.unwrap_or(Point::origin());
        let translation = target_point - current_com;

        structure.par_residues_mut().for_each(|residue| {
            for atom in residue.iter_atoms_mut() {
                atom.translate_by(&translation);
            }
        });
    }

    /// Rotates the structure about the x-axis by the specified angle.
    ///
    /// # Arguments
    ///
    /// * `structure` - Mutable structure to be rotated.
    /// * `radians` - Rotation angle in radians.
    pub fn rotate_x(structure: &mut Structure, radians: f64) {
        let rotation = Rotation3::from_axis_angle(&Vector3::x_axis(), radians);
        Self::apply_rotation(structure, rotation);
    }

    /// Rotates the structure about the y-axis by the specified angle.
    ///
    /// # Arguments
    ///
    /// * `structure` - Mutable structure to be rotated.
    /// * `radians` - Rotation angle in radians.
    pub fn rotate_y(structure: &mut Structure, radians: f64) {
        let rotation = Rotation3::from_axis_angle(&Vector3::y_axis(), radians);
        Self::apply_rotation(structure, rotation);
    }

    /// Rotates the structure about the z-axis by the specified angle.
    ///
    /// # Arguments
    ///
    /// * `structure` - Mutable structure to be rotated.
    /// * `radians` - Rotation angle in radians.
    pub fn rotate_z(structure: &mut Structure, radians: f64) {
        let rotation = Rotation3::from_axis_angle(&Vector3::z_axis(), radians);
        Self::apply_rotation(structure, rotation);
    }

    /// Rotates the structure using Euler angles (XYZ convention).
    ///
    /// # Arguments
    ///
    /// * `structure` - Mutable structure to be rotated.
    /// * `x_rad` - Rotation about x-axis in radians.
    /// * `y_rad` - Rotation about y-axis in radians.
    /// * `z_rad` - Rotation about z-axis in radians.
    pub fn rotate_euler(structure: &mut Structure, x_rad: f64, y_rad: f64, z_rad: f64) {
        let rotation = Rotation3::from_euler_angles(x_rad, y_rad, z_rad);
        Self::apply_rotation(structure, rotation);
    }

    /// Applies a rotation matrix to all atoms and box vectors.
    fn apply_rotation(structure: &mut Structure, rotation: Rotation3<f64>) {
        structure.par_residues_mut().for_each(|residue| {
            for atom in residue.iter_atoms_mut() {
                atom.pos = rotation * atom.pos;
            }
        });

        if let Some(box_vecs) = structure.box_vectors {
            let v1 = Vector3::from(box_vecs[0]);
            let v2 = Vector3::from(box_vecs[1]);
            let v3 = Vector3::from(box_vecs[2]);

            let v1_rot = rotation * v1;
            let v2_rot = rotation * v2;
            let v3_rot = rotation * v3;

            structure.box_vectors = Some([v1_rot.into(), v2_rot.into(), v3_rot.into()]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Transform;
    use crate::model::{
        atom::Atom,
        chain::Chain,
        residue::Residue,
        structure::Structure,
        types::{Element, Point, ResidueCategory, StandardResidue},
    };

    fn structure_with_points(points: &[Point]) -> Structure {
        let mut chain = Chain::new("A");
        let mut residue = Residue::new(
            1,
            None,
            "GLY",
            Some(StandardResidue::GLY),
            ResidueCategory::Standard,
        );

        for (idx, point) in points.iter().enumerate() {
            let name = format!("C{}", idx);
            residue.add_atom(Atom::new(&name, Element::C, *point));
        }

        chain.add_residue(residue);
        let mut structure = Structure::new();
        structure.add_chain(chain);
        structure
    }

    fn assert_point_close(actual: &Point, expected: &Point) {
        assert!((actual.x - expected.x).abs() < 1e-6);
        assert!((actual.y - expected.y).abs() < 1e-6);
        assert!((actual.z - expected.z).abs() < 1e-6);
    }

    #[test]
    fn translate_moves_all_atoms_by_vector() {
        let mut structure =
            structure_with_points(&[Point::new(0.0, 0.0, 0.0), Point::new(1.0, 2.0, 3.0)]);

        Transform::translate(&mut structure, 5.0, -2.0, 1.5);

        let mut atoms = structure.iter_atoms();
        assert_point_close(&atoms.next().unwrap().pos, &Point::new(5.0, -2.0, 1.5));
        assert_point_close(&atoms.next().unwrap().pos, &Point::new(6.0, 0.0, 4.5));
    }

    #[test]
    fn center_geometry_moves_geometric_center_to_target() {
        let mut structure =
            structure_with_points(&[Point::new(2.0, 0.0, 0.0), Point::new(4.0, 0.0, 0.0)]);

        Transform::center_geometry(&mut structure, Some(Point::new(10.0, 0.0, 0.0)));

        let center = structure.geometric_center();
        assert_point_close(&center, &Point::new(10.0, 0.0, 0.0));
    }

    #[test]
    fn center_mass_moves_center_of_mass_to_origin_by_default() {
        let mut structure =
            structure_with_points(&[Point::new(2.0, 0.0, 0.0), Point::new(4.0, 0.0, 0.0)]);

        {
            let chain = structure.chain_mut("A").unwrap();
            let residue = chain.residue_mut(1, None).unwrap();
            residue
                .iter_atoms_mut()
                .enumerate()
                .for_each(|(idx, atom)| {
                    atom.element = if idx == 0 { Element::H } else { Element::O };
                });
        }

        Transform::center_mass(&mut structure, None);

        let com = structure.center_of_mass();
        assert_point_close(&com, &Point::origin());
    }

    #[test]
    fn rotate_z_rotates_atoms_about_origin() {
        let mut structure =
            structure_with_points(&[Point::new(1.0, 0.0, 0.0), Point::new(0.0, 2.0, 0.0)]);

        Transform::rotate_z(&mut structure, std::f64::consts::FRAC_PI_2);

        let mut atoms = structure.iter_atoms();
        assert_point_close(&atoms.next().unwrap().pos, &Point::new(0.0, 1.0, 0.0));
        assert_point_close(&atoms.next().unwrap().pos, &Point::new(-2.0, 0.0, 0.0));
    }

    #[test]
    fn rotate_euler_updates_box_vectors() {
        let mut structure = structure_with_points(&[Point::new(1.0, 0.0, 0.0)]);
        structure.box_vectors = Some([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]);

        Transform::rotate_euler(&mut structure, 0.0, 0.0, std::f64::consts::FRAC_PI_2);

        let box_vectors = structure.box_vectors.unwrap();
        assert_point_close(&Point::from(box_vectors[0]), &Point::new(0.0, 1.0, 0.0));
        assert_point_close(&Point::from(box_vectors[1]), &Point::new(-2.0, 0.0, 0.0));
        assert_point_close(&Point::from(box_vectors[2]), &Point::new(0.0, 0.0, 3.0));
    }
}
