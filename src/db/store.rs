use super::loader;
use super::schema::ResidueTemplateFile;
use std::collections::HashMap;
use std::sync::OnceLock;

#[derive(Debug, Clone)]
pub struct InternalTemplate {
    pub schema: ResidueTemplateFile,
}

pub struct DataStore {
    pub templates_by_name: HashMap<String, InternalTemplate>,
}

static STORE: OnceLock<DataStore> = OnceLock::new();

pub fn get_store() -> &'static DataStore {
    STORE.get_or_init(loader::load_all_templates)
}
