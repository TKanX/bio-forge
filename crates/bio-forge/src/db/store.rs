//! In-memory cache for residue templates loaded from TOML resources.
//!
//! The store is initialized once on demand and exposes simple lookup helpers for the rest
//! of the crate. Tests interact with the same structures to verify loading logic.

use super::loader;
use super::schema::ResidueTemplateFile;
use std::collections::HashMap;
use std::sync::OnceLock;

/// Wrapper around a parsed template, preserving the original schema for inspection.
#[derive(Debug, Clone)]
pub struct InternalTemplate {
    /// Deserialized template schema retained for reuse across queries.
    pub schema: ResidueTemplateFile,
}

/// Holds every template indexed by its unique name.
pub struct DataStore {
    /// Map from template name (e.g., `"ALA"`) to the parsed template definition.
    pub templates_by_name: HashMap<String, InternalTemplate>,
}

static STORE: OnceLock<DataStore> = OnceLock::new();

/// Returns the lazily initialized template store singleton.
///
/// The loader walks every TOML file under `templates/` exactly once and caches the result
/// for subsequent topology-building operations.
///
/// # Returns
///
/// A reference to the initialized [`DataStore`].
pub fn get_store() -> &'static DataStore {
    STORE.get_or_init(loader::load_all_templates)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::types::StandardResidue;
    use std::collections::HashMap;

    fn create_mock_template(
        name: &str,
        standard_name: StandardResidue,
        charge: i32,
    ) -> ResidueTemplateFile {
        ResidueTemplateFile {
            info: crate::db::schema::TemplateInfo {
                name: name.to_string(),
                standard_name,
                charge,
            },
            atoms: vec![],
            hydrogens: vec![],
            bonds: vec![],
        }
    }

    #[test]
    fn internal_template_new_creates_correct_template() {
        let schema = create_mock_template("ALA", StandardResidue::ALA, 0);
        let template = InternalTemplate {
            schema: schema.clone(),
        };

        assert_eq!(template.schema.info.name, "ALA");
        assert_eq!(template.schema.info.standard_name, StandardResidue::ALA);
        assert_eq!(template.schema.info.charge, 0);
    }

    #[test]
    fn internal_template_clone_creates_identical_copy() {
        let schema = create_mock_template("GLY", StandardResidue::GLY, 0);
        let template = InternalTemplate { schema };
        let cloned = template.clone();

        assert_eq!(template.schema.info.name, cloned.schema.info.name);
        assert_eq!(
            template.schema.info.standard_name,
            cloned.schema.info.standard_name
        );
        assert_eq!(template.schema.info.charge, cloned.schema.info.charge);
    }

    #[test]
    fn internal_template_debug_formats_correctly() {
        let schema = create_mock_template("VAL", StandardResidue::VAL, 0);
        let template = InternalTemplate { schema };

        let debug_str = format!("{:?}", template);
        assert!(debug_str.contains("InternalTemplate"));
        assert!(debug_str.contains("VAL"));
    }

    #[test]
    fn data_store_new_creates_empty_store() {
        let store = DataStore {
            templates_by_name: HashMap::new(),
        };

        assert!(store.templates_by_name.is_empty());
        assert_eq!(store.templates_by_name.len(), 0);
    }

    #[test]
    fn data_store_with_templates_stores_correctly() {
        let mut templates = HashMap::new();
        let ala_template = InternalTemplate {
            schema: create_mock_template("ALA", StandardResidue::ALA, 0),
        };
        let gly_template = InternalTemplate {
            schema: create_mock_template("GLY", StandardResidue::GLY, 0),
        };

        templates.insert("ALA".to_string(), ala_template);
        templates.insert("GLY".to_string(), gly_template);

        let store = DataStore {
            templates_by_name: templates,
        };

        assert_eq!(store.templates_by_name.len(), 2);
        assert!(store.templates_by_name.contains_key("ALA"));
        assert!(store.templates_by_name.contains_key("GLY"));
        assert!(!store.templates_by_name.contains_key("VAL"));
    }

    #[test]
    fn data_store_get_template_returns_correct_template() {
        let mut templates = HashMap::new();
        let ala_template = InternalTemplate {
            schema: create_mock_template("ALA", StandardResidue::ALA, 0),
        };
        templates.insert("ALA".to_string(), ala_template);

        let store = DataStore {
            templates_by_name: templates,
        };

        let retrieved = store.templates_by_name.get("ALA");
        assert!(retrieved.is_some());
        assert_eq!(retrieved.unwrap().schema.info.name, "ALA");
        assert_eq!(
            retrieved.unwrap().schema.info.standard_name,
            StandardResidue::ALA
        );
    }

    #[test]
    fn data_store_get_template_returns_none_for_nonexistent() {
        let store = DataStore {
            templates_by_name: HashMap::new(),
        };

        let retrieved = store.templates_by_name.get("NONEXISTENT");
        assert!(retrieved.is_none());
    }

    #[test]
    fn get_store_returns_singleton_instance() {
        let store1 = get_store();
        let store2 = get_store();

        assert!(std::ptr::eq(store1, store2));
    }
}
