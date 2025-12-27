//! Internal database API exposing read-only views over residue templates.
//!
//! Callers obtain [`TemplateView`] handles keyed by template name, enabling topology and IO
//! layers to iterate atoms, hydrogens, and bonds without cloning the underlying schema.

mod loader;
mod schema;
mod store;

use crate::model::types::{BondOrder, Element, Point, StandardResidue};

/// Retrieves a template by its canonical name.
///
/// # Arguments
///
/// * `name` - Template identifier such as `"ALA"` or `"HOH"`.
///
/// # Returns
///
/// `Some(TemplateView)` when the template exists, otherwise `None`.
pub fn get_template(name: &str) -> Option<TemplateView<'_>> {
    store::get_store()
        .templates_by_name
        .get(name)
        .map(TemplateView::new)
}

/// Lightweight wrapper granting read-only access to a stored template.
#[derive(Debug, Clone, Copy)]
pub struct TemplateView<'a> {
    inner: &'a store::InternalTemplate,
}

impl<'a> TemplateView<'a> {
    /// Creates a new view from the internal store entry.
    ///
    /// # Arguments
    ///
    /// * `inner` - Reference to the cached template.
    pub fn new(inner: &'a store::InternalTemplate) -> Self {
        Self { inner }
    }

    /// Returns the template's display name.
    ///
    /// # Returns
    ///
    /// The literal name string stored in the schema.
    pub fn name(&self) -> &'a str {
        &self.inner.schema.info.name
    }

    /// Reports the standard residue enum associated with the template.
    ///
    /// # Returns
    ///
    /// The [`StandardResidue`] discriminant stored in metadata.
    pub fn standard_name(&self) -> StandardResidue {
        self.inner.schema.info.standard_name
    }

    /// Returns the net integer charge of the residue.
    ///
    /// # Returns
    ///
    /// Signed charge in electrons.
    pub fn charge(&self) -> i32 {
        self.inner.schema.info.charge
    }

    /// Iterates heavy atoms with their elements and reference coordinates.
    ///
    /// # Returns
    ///
    /// An iterator yielding `(name, element, Point)` tuples preserving declaration order.
    pub fn heavy_atoms(&self) -> impl Iterator<Item = (&'a str, Element, Point)> {
        self.inner
            .schema
            .atoms
            .iter()
            .map(|a| (a.name.as_str(), a.element, Point::from(a.pos)))
    }

    /// Iterates hydrogens, their coordinates, and anchor names.
    ///
    /// # Returns
    ///
    /// An iterator producing `(name, Point, anchor_iter)` tuples where `anchor_iter`
    /// traverses the heavy-atom anchor names as `&str`.
    pub fn hydrogens(
        &self,
    ) -> impl Iterator<Item = (&'a str, Point, impl Iterator<Item = &'a str>)> {
        self.inner.schema.hydrogens.iter().map(|h| {
            (
                h.name.as_str(),
                Point::from(h.pos),
                h.anchors.iter().map(|s| s.as_str()),
            )
        })
    }

    /// Iterates bonds as name pairs plus their bond order.
    ///
    /// # Returns
    ///
    /// An iterator over `(atom1, atom2, BondOrder)` tuples mirroring the schema definition.
    pub fn bonds(&self) -> impl Iterator<Item = (&'a str, &'a str, BondOrder)> {
        self.inner
            .schema
            .bonds
            .iter()
            .map(|b| (b.a1.as_str(), b.a2.as_str(), b.order))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::types::{BondOrder, Element};

    fn create_mock_template(
        name: &str,
        standard_name: StandardResidue,
        charge: i32,
        atoms: Vec<schema::TemplateHeavyAtom>,
        hydrogens: Vec<schema::TemplateHydrogen>,
        bonds: Vec<schema::TemplateBond>,
    ) -> store::InternalTemplate {
        let schema = schema::ResidueTemplateFile {
            info: schema::TemplateInfo {
                name: name.to_string(),
                standard_name,
                charge,
            },
            atoms,
            hydrogens,
            bonds,
        };
        store::InternalTemplate { schema }
    }

    fn create_simple_mock_template() -> store::InternalTemplate {
        let atoms = vec![
            schema::TemplateHeavyAtom {
                name: "CA".to_string(),
                element: Element::C,
                pos: [0.0, 0.0, 0.0],
            },
            schema::TemplateHeavyAtom {
                name: "CB".to_string(),
                element: Element::C,
                pos: [1.0, 0.0, 0.0],
            },
        ];

        let hydrogens = vec![schema::TemplateHydrogen {
            name: "HA".to_string(),
            pos: [0.5, 1.0, 0.0],
            anchors: vec!["CA".to_string()],
        }];

        let bonds = vec![schema::TemplateBond {
            a1: "CA".to_string(),
            a2: "CB".to_string(),
            order: BondOrder::Single,
        }];

        create_mock_template("TEST", StandardResidue::ALA, 0, atoms, hydrogens, bonds)
    }

    #[test]
    fn get_template_returns_none_for_nonexistent_template() {
        let result = get_template("NONEXISTENT");
        assert!(result.is_none());
    }

    #[test]
    fn get_template_returns_some_for_existing_template() {
        let result = get_template("ALA");
        assert!(result.is_some());
    }

    #[test]
    fn template_view_name_returns_correct_name() {
        let mock_template = create_simple_mock_template();
        let view = TemplateView::new(&mock_template);

        assert_eq!(view.name(), "TEST");
    }

    #[test]
    fn template_view_standard_name_returns_correct_standard_name() {
        let mock_template = create_simple_mock_template();
        let view = TemplateView::new(&mock_template);

        assert_eq!(view.standard_name(), StandardResidue::ALA);
    }

    #[test]
    fn template_view_charge_returns_correct_charge() {
        let mock_template = create_simple_mock_template();
        let view = TemplateView::new(&mock_template);

        assert_eq!(view.charge(), 0);
    }

    #[test]
    fn template_view_charge_returns_correct_charge_for_charged_template() {
        let atoms = vec![schema::TemplateHeavyAtom {
            name: "N".to_string(),
            element: Element::N,
            pos: [0.0, 0.0, 0.0],
        }];
        let mock_template =
            create_mock_template("LYS", StandardResidue::LYS, 1, atoms, vec![], vec![]);
        let view = TemplateView::new(&mock_template);

        assert_eq!(view.charge(), 1);
    }

    #[test]
    fn template_view_heavy_atoms_returns_correct_iterator() {
        let mock_template = create_simple_mock_template();
        let view = TemplateView::new(&mock_template);

        let heavy_atoms: Vec<_> = view.heavy_atoms().collect();

        assert_eq!(heavy_atoms.len(), 2);

        assert_eq!(heavy_atoms[0].0, "CA");
        assert_eq!(heavy_atoms[0].1, Element::C);
        assert_eq!(heavy_atoms[0].2, Point::new(0.0, 0.0, 0.0));

        assert_eq!(heavy_atoms[1].0, "CB");
        assert_eq!(heavy_atoms[1].1, Element::C);
        assert_eq!(heavy_atoms[1].2, Point::new(1.0, 0.0, 0.0));
    }

    #[test]
    fn template_view_heavy_atoms_returns_empty_iterator_for_template_without_atoms() {
        let mock_template =
            create_mock_template("EMPTY", StandardResidue::GLY, 0, vec![], vec![], vec![]);
        let view = TemplateView::new(&mock_template);

        let heavy_atoms: Vec<_> = view.heavy_atoms().collect();

        assert_eq!(heavy_atoms.len(), 0);
    }

    #[test]
    fn template_view_hydrogens_returns_correct_iterator() {
        let mock_template = create_simple_mock_template();
        let view = TemplateView::new(&mock_template);

        let hydrogens: Vec<_> = view.hydrogens().collect();

        assert_eq!(hydrogens.len(), 1);

        let (name, pos, _) = &hydrogens[0];
        assert_eq!(*name, "HA");
        assert_eq!(*pos, Point::new(0.5, 1.0, 0.0));
    }

    #[test]
    fn template_view_hydrogens_returns_empty_iterator_for_template_without_hydrogens() {
        let atoms = vec![schema::TemplateHeavyAtom {
            name: "CA".to_string(),
            element: Element::C,
            pos: [0.0, 0.0, 0.0],
        }];
        let mock_template =
            create_mock_template("NO_H", StandardResidue::GLY, 0, atoms, vec![], vec![]);
        let view = TemplateView::new(&mock_template);

        let hydrogens: Vec<_> = view.hydrogens().collect();

        assert_eq!(hydrogens.len(), 0);
    }

    #[test]
    fn template_view_hydrogens_handles_multiple_anchors() {
        let atoms = vec![
            schema::TemplateHeavyAtom {
                name: "CA".to_string(),
                element: Element::C,
                pos: [0.0, 0.0, 0.0],
            },
            schema::TemplateHeavyAtom {
                name: "CB".to_string(),
                element: Element::C,
                pos: [1.0, 0.0, 0.0],
            },
        ];

        let hydrogens = vec![schema::TemplateHydrogen {
            name: "HA".to_string(),
            pos: [0.5, 1.0, 0.0],
            anchors: vec!["CA".to_string(), "CB".to_string()],
        }];

        let mock_template = create_mock_template(
            "MULTI_ANCHOR",
            StandardResidue::ALA,
            0,
            atoms,
            hydrogens,
            vec![],
        );
        let view = TemplateView::new(&mock_template);

        let hydrogens_vec: Vec<_> = view.hydrogens().collect();
        assert_eq!(hydrogens_vec.len(), 1);

        let (name, pos, _) = &hydrogens_vec[0];
        assert_eq!(*name, "HA");
        assert_eq!(*pos, Point::new(0.5, 1.0, 0.0));
    }

    #[test]
    fn template_view_bonds_returns_correct_iterator() {
        let mock_template = create_simple_mock_template();
        let view = TemplateView::new(&mock_template);

        let bonds: Vec<_> = view.bonds().collect();

        assert_eq!(bonds.len(), 1);

        let (a1, a2, order) = bonds[0];
        assert_eq!(a1, "CA");
        assert_eq!(a2, "CB");
        assert_eq!(order, BondOrder::Single);
    }

    #[test]
    fn template_view_bonds_returns_empty_iterator_for_template_without_bonds() {
        let atoms = vec![schema::TemplateHeavyAtom {
            name: "CA".to_string(),
            element: Element::C,
            pos: [0.0, 0.0, 0.0],
        }];
        let mock_template =
            create_mock_template("NO_BONDS", StandardResidue::GLY, 0, atoms, vec![], vec![]);
        let view = TemplateView::new(&mock_template);

        let bonds: Vec<_> = view.bonds().collect();

        assert_eq!(bonds.len(), 0);
    }

    #[test]
    fn template_view_bonds_handles_different_bond_orders() {
        let atoms = vec![
            schema::TemplateHeavyAtom {
                name: "C1".to_string(),
                element: Element::C,
                pos: [0.0, 0.0, 0.0],
            },
            schema::TemplateHeavyAtom {
                name: "C2".to_string(),
                element: Element::C,
                pos: [1.0, 0.0, 0.0],
            },
            schema::TemplateHeavyAtom {
                name: "N1".to_string(),
                element: Element::N,
                pos: [2.0, 0.0, 0.0],
            },
        ];

        let bonds = vec![
            schema::TemplateBond {
                a1: "C1".to_string(),
                a2: "C2".to_string(),
                order: BondOrder::Single,
            },
            schema::TemplateBond {
                a1: "C2".to_string(),
                a2: "N1".to_string(),
                order: BondOrder::Double,
            },
            schema::TemplateBond {
                a1: "C1".to_string(),
                a2: "N1".to_string(),
                order: BondOrder::Triple,
            },
        ];

        let mock_template =
            create_mock_template("BONDS", StandardResidue::ALA, 0, atoms, vec![], bonds);
        let view = TemplateView::new(&mock_template);

        let bonds_vec: Vec<_> = view.bonds().collect();

        assert_eq!(bonds_vec.len(), 3);
        assert_eq!(bonds_vec[0].2, BondOrder::Single);
        assert_eq!(bonds_vec[1].2, BondOrder::Double);
        assert_eq!(bonds_vec[2].2, BondOrder::Triple);
    }

    #[test]
    fn template_view_bonds_handles_aromatic_bonds() {
        let atoms = vec![
            schema::TemplateHeavyAtom {
                name: "C1".to_string(),
                element: Element::C,
                pos: [0.0, 0.0, 0.0],
            },
            schema::TemplateHeavyAtom {
                name: "C2".to_string(),
                element: Element::C,
                pos: [1.0, 0.0, 0.0],
            },
        ];

        let bonds = vec![schema::TemplateBond {
            a1: "C1".to_string(),
            a2: "C2".to_string(),
            order: BondOrder::Aromatic,
        }];

        let mock_template =
            create_mock_template("AROMATIC", StandardResidue::PHE, 0, atoms, vec![], bonds);
        let view = TemplateView::new(&mock_template);

        let bonds_vec: Vec<_> = view.bonds().collect();

        assert_eq!(bonds_vec.len(), 1);
        assert_eq!(bonds_vec[0].2, BondOrder::Aromatic);
    }

    #[test]
    fn template_view_clone_creates_identical_copy() {
        let mock_template = create_simple_mock_template();
        let view1 = TemplateView::new(&mock_template);
        let view2 = view1;

        assert_eq!(view1.name(), view2.name());
        assert_eq!(view1.standard_name(), view2.standard_name());
        assert_eq!(view1.charge(), view2.charge());

        let atoms1: Vec<_> = view1.heavy_atoms().collect();
        let atoms2: Vec<_> = view2.heavy_atoms().collect();
        assert_eq!(atoms1, atoms2);

        let bonds1: Vec<_> = view1.bonds().collect();
        let bonds2: Vec<_> = view2.bonds().collect();
        assert_eq!(bonds1, bonds2);
    }

    #[test]
    fn template_view_debug_formatting() {
        let mock_template = create_simple_mock_template();
        let view = TemplateView::new(&mock_template);

        let debug_str = format!("{:?}", view);
        assert!(debug_str.contains("TemplateView"));
    }

    #[test]
    fn get_template_with_real_database_templates() {
        let templates_to_test = vec![
            "ALA", "GLY", "VAL", "LEU", "ILE", "MET", "PHE", "PRO", "SER", "THR", "TYR", "CYS",
            "ASN", "GLN", "ASP", "GLU", "LYS", "ARG", "HID", "HIE", "HIP", "HOH",
        ];

        for template_name in templates_to_test {
            let result = get_template(template_name);
            assert!(
                result.is_some(),
                "Template '{}' should exist",
                template_name
            );

            let template = result.unwrap();
            assert_eq!(template.name(), template_name);
        }
    }

    #[test]
    fn template_view_with_real_template_has_expected_properties() {
        if let Some(ala_template) = get_template("ALA") {
            assert_eq!(ala_template.name(), "ALA");
            assert_eq!(ala_template.standard_name(), StandardResidue::ALA);
            assert_eq!(ala_template.charge(), 0);

            let heavy_atoms: Vec<_> = ala_template.heavy_atoms().collect();
            assert!(!heavy_atoms.is_empty());

            let bonds: Vec<_> = ala_template.bonds().collect();
            assert!(!bonds.is_empty());
        }
    }

    #[test]
    fn template_view_with_real_template_water_properties() {
        if let Some(hoh_template) = get_template("HOH") {
            assert_eq!(hoh_template.name(), "HOH");
            assert_eq!(hoh_template.standard_name(), StandardResidue::HOH);
            assert_eq!(hoh_template.charge(), 0);

            let heavy_atoms: Vec<_> = hoh_template.heavy_atoms().collect();
            assert_eq!(heavy_atoms.len(), 1);
            assert_eq!(heavy_atoms[0].1, Element::O);

            let hydrogens: Vec<_> = hoh_template.hydrogens().collect();
            assert!(!hydrogens.is_empty());
        }
    }
}
