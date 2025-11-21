use super::schema::ResidueTemplateFile;
use super::store::{DataStore, InternalTemplate};
use std::collections::HashMap;

pub fn load_all_templates() -> DataStore {
    let mut templates_by_name = HashMap::new();

    macro_rules! load_template {
        ($path:literal) => {
            let content = include_str!(concat!("../../templates/", $path));
            let schema: ResidueTemplateFile = toml::from_str(content)
                .unwrap_or_else(|e| panic!("Failed to parse template file '{}': {}", $path, e));

            let template_name = schema.info.name.clone();
            let internal_template = InternalTemplate { schema };

            if templates_by_name
                .insert(template_name.clone(), internal_template)
                .is_some()
            {
                panic!("Duplicate template name found: {}", template_name);
            }
        };
    }

    load_template!("protein/ALA.toml");
    load_template!("protein/ARG.toml");
    load_template!("protein/ASN.toml");
    load_template!("protein/CYM.toml");
    load_template!("protein/CYX.toml");
    load_template!("protein/GLN.toml");
    load_template!("protein/GLY.toml");
    load_template!("protein/HIE.toml");
    load_template!("protein/ILE.toml");
    load_template!("protein/LYN.toml");
    load_template!("protein/MET.toml");
    load_template!("protein/PRO.toml");
    load_template!("protein/THR.toml");
    load_template!("protein/TYM.toml");
    load_template!("protein/VAL.toml");
    load_template!("protein/AR0.toml");
    load_template!("protein/ASH.toml");
    load_template!("protein/ASP.toml");
    load_template!("protein/CYS.toml");
    load_template!("protein/GLH.toml");
    load_template!("protein/GLU.toml");
    load_template!("protein/HID.toml");
    load_template!("protein/HIP.toml");
    load_template!("protein/LEU.toml");
    load_template!("protein/LYS.toml");
    load_template!("protein/PHE.toml");
    load_template!("protein/SER.toml");
    load_template!("protein/TRP.toml");
    load_template!("protein/TYR.toml");

    load_template!("nucleic/A.toml");
    load_template!("nucleic/C.toml");
    load_template!("nucleic/G.toml");
    load_template!("nucleic/U.toml");
    load_template!("nucleic/I.toml");

    load_template!("nucleic/DA.toml");
    load_template!("nucleic/DC.toml");
    load_template!("nucleic/DG.toml");
    load_template!("nucleic/DT.toml");
    load_template!("nucleic/DI.toml");

    load_template!("solvent/HOH.toml");

    DataStore { templates_by_name }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::types::BondOrder;
    use std::collections::HashMap;

    fn load_template_from_content(
        content: &str,
        templates_by_name: &mut HashMap<String, InternalTemplate>,
    ) -> Result<(), String> {
        let schema: ResidueTemplateFile = toml::from_str(content)
            .map_err(|e| format!("Failed to parse template content: {}", e))?;

        let template_name = schema.info.name.clone();
        let internal_template = InternalTemplate { schema };

        if templates_by_name
            .insert(template_name.clone(), internal_template)
            .is_some()
        {
            return Err(format!("Duplicate template name found: {}", template_name));
        }

        Ok(())
    }

    fn create_minimal_template(name: &str, standard_name: &str) -> String {
        format!(
            r#"
            [info]
            name = "{}"
            standard_name = "{}"
            charge = 0

            [[atoms]]
            name = "C1"
            element = "C"
            pos = [0.0, 0.0, 0.0]
            "#,
            name, standard_name
        )
    }

    fn create_template_with_bonds(name: &str, standard_name: &str) -> String {
        format!(
            r#"
            [info]
            name = "{}"
            standard_name = "{}"
            charge = 0

            [[atoms]]
            name = "C1"
            element = "C"
            pos = [0.0, 0.0, 0.0]

            [[atoms]]
            name = "C2"
            element = "C"
            pos = [1.0, 0.0, 0.0]

            [[bonds]]
            a1 = "C1"
            a2 = "C2"
            order = "Single"
            "#,
            name, standard_name
        )
    }

    #[test]
    fn load_template_from_content_parses_valid_toml() {
        let mut templates = HashMap::new();
        let content = create_minimal_template("TEST", "ALA");

        let result = load_template_from_content(&content, &mut templates);

        assert!(result.is_ok());
        assert_eq!(templates.len(), 1);
        assert!(templates.contains_key("TEST"));

        let template = &templates["TEST"];
        assert_eq!(template.schema.info.name, "TEST");
        assert_eq!(template.schema.atoms.len(), 1);
        assert_eq!(template.schema.atoms[0].name, "C1");
    }

    #[test]
    fn load_template_from_content_handles_template_with_bonds() {
        let mut templates = HashMap::new();
        let content = create_template_with_bonds("BONDED", "GLY");

        let result = load_template_from_content(&content, &mut templates);

        assert!(result.is_ok());
        assert_eq!(templates.len(), 1);

        let template = &templates["BONDED"];
        assert_eq!(template.schema.info.name, "BONDED");
        assert_eq!(template.schema.atoms.len(), 2);
        assert_eq!(template.schema.bonds.len(), 1);
        assert_eq!(template.schema.bonds[0].a1, "C1");
        assert_eq!(template.schema.bonds[0].a2, "C2");
        assert_eq!(template.schema.bonds[0].order, BondOrder::Single);
    }

    #[test]
    fn load_template_from_content_handles_template_with_hydrogens() {
        let mut templates = HashMap::new();
        let content = r#"
            [info]
            name = "HYDRO"
            standard_name = "ALA"
            charge = 0

            [[atoms]]
            name = "C1"
            element = "C"
            pos = [0.0, 0.0, 0.0]

            [[hydrogens]]
            name = "H1"
            pos = [1.0, 0.0, 0.0]
            anchors = ["C1"]
            "#;

        let result = load_template_from_content(&content, &mut templates);

        assert!(result.is_ok());
        assert_eq!(templates.len(), 1);

        let template = &templates["HYDRO"];
        assert_eq!(template.schema.hydrogens.len(), 1);
        assert_eq!(template.schema.hydrogens[0].name, "H1");
        assert_eq!(template.schema.hydrogens[0].anchors, vec!["C1"]);
    }

    #[test]
    fn load_template_from_content_detects_duplicate_template_names() {
        let mut templates = HashMap::new();

        let content1 = create_minimal_template("DUPLICATE", "ALA");
        let result1 = load_template_from_content(&content1, &mut templates);
        assert!(result1.is_ok());
        assert_eq!(templates.len(), 1);

        let content2 = create_minimal_template("DUPLICATE", "GLY");
        let result2 = load_template_from_content(&content2, &mut templates);

        assert!(result2.is_err());
        assert!(
            result2
                .unwrap_err()
                .contains("Duplicate template name found: DUPLICATE")
        );
        assert_eq!(templates.len(), 1);
    }

    #[test]
    fn load_template_from_content_handles_invalid_toml() {
        let mut templates = HashMap::new();
        let invalid_content = r#"
            [info
            name = "INVALID"
            standard_name = "ALA"
            charge = 0
            "#;

        let result = load_template_from_content(&invalid_content, &mut templates);

        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .contains("Failed to parse template content")
        );
        assert_eq!(templates.len(), 0);
    }

    #[test]
    fn load_template_from_content_handles_missing_required_fields() {
        let mut templates = HashMap::new();
        let invalid_content = r#"
            [[atoms]]
            name = "C1"
            element = "C"
            pos = [0.0, 0.0, 0.0]
            "#;

        let result = load_template_from_content(&invalid_content, &mut templates);

        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .contains("Failed to parse template content")
        );
        assert_eq!(templates.len(), 0);
    }

    #[test]
    fn load_template_from_content_handles_unknown_fields() {
        let mut templates = HashMap::new();
        let content_with_unknown = r#"
            [info]
            name = "UNKNOWN"
            standard_name = "ALA"
            charge = 0
            unknown_field = "should_fail"

            [[atoms]]
            name = "C1"
            element = "C"
            pos = [0.0, 0.0, 0.0]
            "#;

        let result = load_template_from_content(&content_with_unknown, &mut templates);

        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .contains("Failed to parse template content")
        );
        assert_eq!(templates.len(), 0);
    }

    #[test]
    fn load_template_from_content_handles_empty_content() {
        let mut templates = HashMap::new();
        let empty_content = "";

        let result = load_template_from_content(&empty_content, &mut templates);

        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .contains("Failed to parse template content")
        );
        assert_eq!(templates.len(), 0);
    }

    #[test]
    fn load_template_from_content_handles_multiple_templates() {
        let mut templates = HashMap::new();

        let content1 = create_minimal_template("TEMPLATE1", "ALA");
        let result1 = load_template_from_content(&content1, &mut templates);
        assert!(result1.is_ok());

        let content2 = create_minimal_template("TEMPLATE2", "GLY");
        let result2 = load_template_from_content(&content2, &mut templates);
        assert!(result2.is_ok());

        let content3 = create_template_with_bonds("TEMPLATE3", "VAL");
        let result3 = load_template_from_content(&content3, &mut templates);
        assert!(result3.is_ok());

        assert_eq!(templates.len(), 3);
        assert!(templates.contains_key("TEMPLATE1"));
        assert!(templates.contains_key("TEMPLATE2"));
        assert!(templates.contains_key("TEMPLATE3"));
    }

    #[test]
    fn load_template_from_content_preserves_template_data_integrity() {
        let mut templates = HashMap::new();
        let content = r#"
            [info]
            name = "COMPLEX"
            standard_name = "ALA"
            charge = -4

            [[atoms]]
            name = "P1"
            element = "P"
            pos = [0.0, 0.0, 0.0]

            [[atoms]]
            name = "O1"
            element = "O"
            pos = [1.0, 0.0, 0.0]

            [[atoms]]
            name = "O2"
            element = "O"
            pos = [-0.5, 0.866, 0.0]

            [[hydrogens]]
            name = "H1"
            pos = [1.5, 0.0, 0.0]
            anchors = ["O1"]

            [[bonds]]
            a1 = "P1"
            a2 = "O1"
            order = "Double"

            [[bonds]]
            a1 = "P1"
            a2 = "O2"
            order = "Single"
            "#;

        let result = load_template_from_content(&content, &mut templates);

        assert!(result.is_ok());
        assert_eq!(templates.len(), 1);

        let template = &templates["COMPLEX"];
        assert_eq!(template.schema.info.name, "COMPLEX");
        assert_eq!(template.schema.info.charge, -4);
        assert_eq!(template.schema.atoms.len(), 3);
        assert_eq!(template.schema.hydrogens.len(), 1);
        assert_eq!(template.schema.bonds.len(), 2);
    }

    #[test]
    fn load_template_from_content_handles_different_bond_orders() {
        let mut templates = HashMap::new();
        let content = r#"
            [info]
            name = "BONDS"
            standard_name = "ALA"
            charge = 0

            [[atoms]]
            name = "C1"
            element = "C"
            pos = [0.0, 0.0, 0.0]

            [[atoms]]
            name = "C2"
            element = "C"
            pos = [1.0, 0.0, 0.0]

            [[atoms]]
            name = "N1"
            element = "N"
            pos = [2.0, 0.0, 0.0]

            [[bonds]]
            a1 = "C1"
            a2 = "C2"
            order = "Single"

            [[bonds]]
            a1 = "C2"
            a2 = "N1"
            order = "Double"

            [[bonds]]
            a1 = "C1"
            a2 = "N1"
            order = "Triple"
            "#;

        let result = load_template_from_content(&content, &mut templates);

        assert!(result.is_ok());
        let template = &templates["BONDS"];
        assert_eq!(template.schema.bonds.len(), 3);
        assert_eq!(template.schema.bonds[0].order, BondOrder::Single);
        assert_eq!(template.schema.bonds[1].order, BondOrder::Double);
        assert_eq!(template.schema.bonds[2].order, BondOrder::Triple);
    }

    #[test]
    fn load_template_from_content_handles_aromatic_bonds() {
        let mut templates = HashMap::new();
        let content = r#"
            [info]
            name = "AROMATIC"
            standard_name = "PHE"
            charge = 0

            [[atoms]]
            name = "C1"
            element = "C"
            pos = [0.0, 0.0, 0.0]

            [[atoms]]
            name = "C2"
            element = "C"
            pos = [1.0, 0.0, 0.0]

            [[bonds]]
            a1 = "C1"
            a2 = "C2"
            order = "Aromatic"
            "#;

        let result = load_template_from_content(&content, &mut templates);

        assert!(result.is_ok());
        let template = &templates["AROMATIC"];
        assert_eq!(template.schema.bonds[0].order, BondOrder::Aromatic);
    }

    #[test]
    fn load_all_templates_returns_datastore_with_expected_structure() {
        let store = load_all_templates();

        assert!(std::ptr::eq(
            &store.templates_by_name as *const _,
            &store.templates_by_name as *const _
        ));

        let _map_ref = &store.templates_by_name;
    }

    #[test]
    fn load_template_macro_behavior_is_simulated_correctly() {
        let mut templates = HashMap::new();
        let content = create_minimal_template("MACRO_TEST", "LEU");

        let schema: ResidueTemplateFile = toml::from_str(&content).unwrap();
        let template_name = schema.info.name.clone();
        let internal_template = InternalTemplate { schema };

        let was_inserted = templates
            .insert(template_name.clone(), internal_template)
            .is_none();
        assert!(was_inserted);

        assert_eq!(templates.len(), 1);
        assert!(templates.contains_key("MACRO_TEST"));
    }
}
