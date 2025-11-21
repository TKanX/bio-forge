use crate::model::types::{BondOrder, Element, StandardResidue};
use serde::Deserialize;

#[derive(Debug, Deserialize, Clone)]
#[serde(deny_unknown_fields)]
pub struct ResidueTemplateFile {
    pub info: TemplateInfo,
    #[serde(default)]
    pub atoms: Vec<TemplateHeavyAtom>,
    #[serde(default)]
    pub hydrogens: Vec<TemplateHydrogen>,
    #[serde(default)]
    pub bonds: Vec<TemplateBond>,
}

#[derive(Debug, Deserialize, Clone)]
#[serde(deny_unknown_fields)]
pub struct TemplateInfo {
    pub name: String,
    pub standard_name: StandardResidue,
    pub charge: i32,
}

#[derive(Debug, Deserialize, Clone)]
#[serde(deny_unknown_fields)]
pub struct TemplateHeavyAtom {
    pub name: String,
    pub element: Element,
    pub pos: [f64; 3],
}

#[derive(Debug, Deserialize, Clone)]
#[serde(deny_unknown_fields)]
pub struct TemplateHydrogen {
    pub name: String,
    pub pos: [f64; 3],
    pub anchors: [String; 3],
}

#[derive(Debug, Deserialize, Clone)]
#[serde(deny_unknown_fields)]
pub struct TemplateBond {
    pub a1: String,
    pub a2: String,
    pub order: BondOrder,
}
