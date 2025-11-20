use super::types::BondOrder;
use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Template {
    pub name: String,
    atom_names: Vec<String>,
    bonds: Vec<(String, String, BondOrder)>,
}

impl Template {
    pub fn new<S: Into<String>>(
        name: S,
        atom_names: Vec<String>,
        bonds: Vec<(String, String, BondOrder)>,
    ) -> Self {
        debug_assert!(
            bonds
                .iter()
                .all(|(a1, a2, _)| { atom_names.contains(a1) && atom_names.contains(a2) }),
            "Bond in template '{}' refers to an atom name that does not exist in the atom list.",
            name.into()
        );

        Self {
            name: name.into(),
            atom_names,
            bonds,
        }
    }

    pub fn has_bond(&self, name1: &str, name2: &str) -> bool {
        self.bonds
            .iter()
            .any(|(a1, a2, _)| (a1 == name1 && a2 == name2) || (a1 == name2 && a2 == name1))
    }

    pub fn has_atom(&self, name: &str) -> bool {
        self.atom_names.contains(&name.to_string())
    }

    pub fn atom_names(&self) -> &[String] {
        &self.atom_names
    }

    pub fn bonds(&self) -> &[(String, String, BondOrder)] {
        &self.bonds
    }

    pub fn atom_count(&self) -> usize {
        self.atom_names.len()
    }

    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }
}

impl fmt::Display for Template {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Template {{ name: \"{}\", atoms: {}, bonds: {} }}",
            self.name,
            self.atom_count(),
            self.bond_count()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn template_new_creates_correct_template() {
        let atom_names = vec!["C1".to_string(), "C2".to_string(), "O1".to_string()];
        let bonds = vec![
            ("C1".to_string(), "C2".to_string(), BondOrder::Single),
            ("C2".to_string(), "O1".to_string(), BondOrder::Double),
        ];

        let template = Template::new("LIG", atom_names.clone(), bonds.clone());

        assert_eq!(template.name, "LIG");
        assert_eq!(template.atom_names, atom_names);
        assert_eq!(template.bonds, bonds);
    }

    #[test]
    fn template_new_with_empty_atom_names() {
        let atom_names = Vec::new();
        let bonds = Vec::new();

        let template = Template::new("EMPTY", atom_names, bonds);

        assert_eq!(template.name, "EMPTY");
        assert_eq!(template.atom_count(), 0);
        assert_eq!(template.bond_count(), 0);
    }

    #[test]
    fn template_new_with_empty_bonds() {
        let atom_names = vec!["C1".to_string(), "C2".to_string()];
        let bonds = Vec::new();

        let template = Template::new("NO_BONDS", atom_names.clone(), bonds);

        assert_eq!(template.name, "NO_BONDS");
        assert_eq!(template.atom_names, atom_names);
        assert_eq!(template.bonds, Vec::new());
    }

    #[test]
    fn template_new_with_string_name() {
        let atom_names = vec!["N1".to_string()];
        let bonds = Vec::new();

        let template = Template::new(String::from("ATP"), atom_names.clone(), bonds);

        assert_eq!(template.name, "ATP");
        assert_eq!(template.atom_names, atom_names);
    }

    #[test]
    fn template_has_bond_returns_true_for_existing_bond() {
        let atom_names = vec!["C1".to_string(), "C2".to_string(), "O1".to_string()];
        let bonds = vec![
            ("C1".to_string(), "C2".to_string(), BondOrder::Single),
            ("C2".to_string(), "O1".to_string(), BondOrder::Double),
        ];
        let template = Template::new("LIG", atom_names, bonds);

        assert!(template.has_bond("C1", "C2"));
        assert!(template.has_bond("C2", "O1"));
    }

    #[test]
    fn template_has_bond_returns_true_for_reverse_order() {
        let atom_names = vec!["C1".to_string(), "C2".to_string()];
        let bonds = vec![("C1".to_string(), "C2".to_string(), BondOrder::Single)];
        let template = Template::new("LIG", atom_names, bonds);

        assert!(template.has_bond("C2", "C1"));
    }

    #[test]
    fn template_has_bond_returns_false_for_nonexistent_bond() {
        let atom_names = vec!["C1".to_string(), "C2".to_string(), "O1".to_string()];
        let bonds = vec![("C1".to_string(), "C2".to_string(), BondOrder::Single)];
        let template = Template::new("LIG", atom_names, bonds);

        assert!(!template.has_bond("C1", "O1"));
        assert!(!template.has_bond("O1", "C2"));
    }

    #[test]
    fn template_has_bond_returns_false_for_empty_template() {
        let template = Template::new("EMPTY", Vec::new(), Vec::new());

        assert!(!template.has_bond("C1", "C2"));
    }

    #[test]
    fn template_has_atom_returns_true_for_existing_atom() {
        let atom_names = vec!["C1".to_string(), "C2".to_string(), "O1".to_string()];
        let bonds = Vec::new();
        let template = Template::new("LIG", atom_names, bonds);

        assert!(template.has_atom("C1"));
        assert!(template.has_atom("C2"));
        assert!(template.has_atom("O1"));
    }

    #[test]
    fn template_has_atom_returns_false_for_nonexistent_atom() {
        let atom_names = vec!["C1".to_string(), "C2".to_string()];
        let bonds = Vec::new();
        let template = Template::new("LIG", atom_names, bonds);

        assert!(!template.has_atom("O1"));
        assert!(!template.has_atom("N1"));
    }

    #[test]
    fn template_has_atom_returns_false_for_empty_template() {
        let template = Template::new("EMPTY", Vec::new(), Vec::new());

        assert!(!template.has_atom("C1"));
    }

    #[test]
    fn template_atom_names_returns_correct_slice() {
        let atom_names = vec!["C1".to_string(), "C2".to_string(), "O1".to_string()];
        let bonds = Vec::new();
        let template = Template::new("LIG", atom_names.clone(), bonds);

        let names = template.atom_names();

        assert_eq!(names, atom_names.as_slice());
        assert_eq!(names.len(), 3);
        assert_eq!(names[0], "C1");
        assert_eq!(names[1], "C2");
        assert_eq!(names[2], "O1");
    }

    #[test]
    fn template_atom_names_returns_empty_slice_for_empty_template() {
        let template = Template::new("EMPTY", Vec::new(), Vec::new());

        let names = template.atom_names();

        assert_eq!(names, &[] as &[String]);
        assert_eq!(names.len(), 0);
    }

    #[test]
    fn template_bonds_returns_correct_slice() {
        let atom_names = vec!["C1".to_string(), "C2".to_string(), "O1".to_string()];
        let bonds = vec![
            ("C1".to_string(), "C2".to_string(), BondOrder::Single),
            ("C2".to_string(), "O1".to_string(), BondOrder::Double),
        ];
        let template = Template::new("LIG", atom_names, bonds.clone());

        let template_bonds = template.bonds();

        assert_eq!(template_bonds, bonds.as_slice());
        assert_eq!(template_bonds.len(), 2);
        assert_eq!(
            template_bonds[0],
            ("C1".to_string(), "C2".to_string(), BondOrder::Single)
        );
        assert_eq!(
            template_bonds[1],
            ("C2".to_string(), "O1".to_string(), BondOrder::Double)
        );
    }

    #[test]
    fn template_bonds_returns_empty_slice_for_empty_template() {
        let template = Template::new("EMPTY", Vec::new(), Vec::new());

        let bonds = template.bonds();

        assert_eq!(bonds, &[]);
        assert_eq!(bonds.len(), 0);
    }

    #[test]
    fn template_atom_count_returns_correct_count() {
        let atom_names = vec![
            "C1".to_string(),
            "C2".to_string(),
            "O1".to_string(),
            "N1".to_string(),
        ];
        let bonds = Vec::new();
        let template = Template::new("LIG", atom_names, bonds);

        assert_eq!(template.atom_count(), 4);
    }

    #[test]
    fn template_atom_count_returns_zero_for_empty_template() {
        let template = Template::new("EMPTY", Vec::new(), Vec::new());

        assert_eq!(template.atom_count(), 0);
    }

    #[test]
    fn template_bond_count_returns_correct_count() {
        let atom_names = vec!["C1".to_string(), "C2".to_string(), "O1".to_string()];
        let bonds = vec![
            ("C1".to_string(), "C2".to_string(), BondOrder::Single),
            ("C2".to_string(), "O1".to_string(), BondOrder::Double),
            ("C1".to_string(), "O1".to_string(), BondOrder::Single),
        ];
        let template = Template::new("LIG", atom_names, bonds);

        assert_eq!(template.bond_count(), 3);
    }

    #[test]
    fn template_bond_count_returns_zero_for_empty_template() {
        let template = Template::new("EMPTY", Vec::new(), Vec::new());

        assert_eq!(template.bond_count(), 0);
    }

    #[test]
    fn template_display_formats_correctly() {
        let atom_names = vec!["C1".to_string(), "C2".to_string(), "O1".to_string()];
        let bonds = vec![
            ("C1".to_string(), "C2".to_string(), BondOrder::Single),
            ("C2".to_string(), "O1".to_string(), BondOrder::Double),
        ];
        let template = Template::new("LIG", atom_names, bonds);

        let display = format!("{}", template);
        let expected = "Template { name: \"LIG\", atoms: 3, bonds: 2 }";

        assert_eq!(display, expected);
    }

    #[test]
    fn template_display_formats_empty_template_correctly() {
        let template = Template::new("EMPTY", Vec::new(), Vec::new());

        let display = format!("{}", template);
        let expected = "Template { name: \"EMPTY\", atoms: 0, bonds: 0 }";

        assert_eq!(display, expected);
    }

    #[test]
    fn template_display_formats_template_with_special_characters() {
        let atom_names = vec!["C1".to_string()];
        let bonds = Vec::new();
        let template = Template::new("LIG-123", atom_names, bonds);

        let display = format!("{}", template);
        let expected = "Template { name: \"LIG-123\", atoms: 1, bonds: 0 }";

        assert_eq!(display, expected);
    }

    #[test]
    fn template_clone_creates_identical_copy() {
        let atom_names = vec!["C1".to_string(), "C2".to_string()];
        let bonds = vec![("C1".to_string(), "C2".to_string(), BondOrder::Single)];
        let template = Template::new("LIG", atom_names.clone(), bonds.clone());

        let cloned = template.clone();

        assert_eq!(template, cloned);
        assert_eq!(template.name, cloned.name);
        assert_eq!(template.atom_names, cloned.atom_names);
        assert_eq!(template.bonds, cloned.bonds);
    }

    #[test]
    fn template_partial_eq_compares_correctly() {
        let atom_names = vec!["C1".to_string(), "C2".to_string()];
        let bonds = vec![("C1".to_string(), "C2".to_string(), BondOrder::Single)];

        let template1 = Template::new("LIG", atom_names.clone(), bonds.clone());
        let template2 = Template::new("LIG", atom_names.clone(), bonds.clone());
        let template3 = Template::new("DIF", atom_names.clone(), bonds.clone());

        assert_eq!(template1, template2);
        assert_ne!(template1, template3);
    }

    #[test]
    fn template_with_complex_molecule() {
        let atom_names = vec![
            "C1".to_string(),
            "C2".to_string(),
            "O1".to_string(),
            "H11".to_string(),
            "H12".to_string(),
            "H13".to_string(),
            "H21".to_string(),
            "H22".to_string(),
            "H1".to_string(),
        ];
        let bonds = vec![
            ("C1".to_string(), "C2".to_string(), BondOrder::Single),
            ("C2".to_string(), "O1".to_string(), BondOrder::Single),
            ("C1".to_string(), "H11".to_string(), BondOrder::Single),
            ("C1".to_string(), "H12".to_string(), BondOrder::Single),
            ("C1".to_string(), "H13".to_string(), BondOrder::Single),
            ("C2".to_string(), "H21".to_string(), BondOrder::Single),
            ("C2".to_string(), "H22".to_string(), BondOrder::Single),
            ("O1".to_string(), "H1".to_string(), BondOrder::Single),
        ];

        let template = Template::new("ETOH", atom_names, bonds);

        assert_eq!(template.name, "ETOH");
        assert_eq!(template.atom_count(), 9);
        assert_eq!(template.bond_count(), 8);

        assert!(template.has_bond("C1", "C2"));
        assert!(template.has_bond("C2", "O1"));
        assert!(template.has_bond("O1", "H1"));

        assert!(template.has_atom("C1"));
        assert!(template.has_atom("H1"));
        assert!(!template.has_atom("C3"));
    }

    #[test]
    fn template_with_aromatic_bonds() {
        let atom_names = vec!["C1".to_string(), "C2".to_string(), "C3".to_string()];
        let bonds = vec![
            ("C1".to_string(), "C2".to_string(), BondOrder::Aromatic),
            ("C2".to_string(), "C3".to_string(), BondOrder::Aromatic),
            ("C3".to_string(), "C1".to_string(), BondOrder::Aromatic),
        ];

        let template = Template::new("BENZENE_RING", atom_names, bonds);

        assert_eq!(template.bond_count(), 3);
        assert!(template.has_bond("C1", "C2"));
        assert!(template.has_bond("C2", "C3"));
        assert!(template.has_bond("C3", "C1"));
    }

    #[test]
    fn template_with_triple_bond() {
        let atom_names = vec!["N1".to_string(), "N2".to_string()];
        let bonds = vec![("N1".to_string(), "N2".to_string(), BondOrder::Triple)];

        let template = Template::new("N2", atom_names, bonds);

        assert_eq!(template.bond_count(), 1);
        assert!(template.has_bond("N1", "N2"));
        assert!(!template.has_bond("N1", "N3"));
    }
}
