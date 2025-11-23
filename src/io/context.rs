use crate::model::types::StandardResidue;
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct IoContext {
    alias_map: HashMap<String, String>,
    standard_map: HashMap<String, StandardResidue>,
}

impl IoContext {
    pub fn new_default() -> Self {
        let mut alias_map = HashMap::new();
        let mut standard_map = HashMap::new();

        macro_rules! register_standard {
            ($canonical:expr, $enum_val:expr) => {
                alias_map.insert($canonical.to_string(), $canonical.to_string());
                standard_map.insert($canonical.to_string(), $enum_val);
            };
        }

        macro_rules! register_alias {
            ($alias:expr, $canonical:expr) => {
                alias_map.insert($alias.to_string(), $canonical.to_string());
            };
        }

        register_standard!("ALA", StandardResidue::ALA);
        register_standard!("ARG", StandardResidue::ARG);
        register_standard!("ARN", StandardResidue::ARG);
        register_standard!("ASN", StandardResidue::ASN);
        register_standard!("ASP", StandardResidue::ASP);
        register_standard!("ASH", StandardResidue::ASP);
        register_standard!("CYS", StandardResidue::CYS);
        register_standard!("CYM", StandardResidue::CYS);
        register_standard!("CYX", StandardResidue::CYS);
        register_standard!("GLN", StandardResidue::GLN);
        register_standard!("GLU", StandardResidue::GLU);
        register_standard!("GLH", StandardResidue::GLU);
        register_standard!("GLY", StandardResidue::GLY);
        register_standard!("HID", StandardResidue::HIS);
        register_standard!("HIE", StandardResidue::HIS);
        register_standard!("HIP", StandardResidue::HIS);
        register_standard!("ILE", StandardResidue::ILE);
        register_standard!("LEU", StandardResidue::LEU);
        register_standard!("LYS", StandardResidue::LYS);
        register_standard!("LYN", StandardResidue::LYS);
        register_standard!("MET", StandardResidue::MET);
        register_standard!("PHE", StandardResidue::PHE);
        register_standard!("PRO", StandardResidue::PRO);
        register_standard!("SER", StandardResidue::SER);
        register_standard!("THR", StandardResidue::THR);
        register_standard!("TRP", StandardResidue::TRP);
        register_standard!("TYR", StandardResidue::TYR);
        register_standard!("TYM", StandardResidue::TYR);
        register_standard!("VAL", StandardResidue::VAL);

        register_standard!("DA", StandardResidue::DA);
        register_standard!("DC", StandardResidue::DC);
        register_standard!("DG", StandardResidue::DG);
        register_standard!("DT", StandardResidue::DT);
        register_standard!("DI", StandardResidue::DI);

        register_standard!("A", StandardResidue::A);
        register_standard!("C", StandardResidue::C);
        register_standard!("G", StandardResidue::G);
        register_standard!("U", StandardResidue::U);
        register_standard!("I", StandardResidue::I);

        register_standard!("HOH", StandardResidue::HOH);

        register_alias!("AIB", "ALA");
        register_alias!("ALM", "ALA");
        register_alias!("AYA", "ALA");
        register_alias!("BNN", "ALA");
        register_alias!("CHG", "ALA");
        register_alias!("CSD", "ALA");
        register_alias!("DAL", "ALA");
        register_alias!("DHA", "ALA");
        register_alias!("DNP", "ALA");
        register_alias!("FLA", "ALA");
        register_alias!("HAC", "ALA");
        register_alias!("MAA", "ALA");
        register_alias!("PRR", "ALA");
        register_alias!("TIH", "ALA");
        register_alias!("TPQ", "ALA");

        register_alias!("ACL", "ARG");
        register_alias!("AGM", "ARG");
        register_alias!("ARM", "ARG");
        register_alias!("DAR", "ARG");
        register_alias!("HAR", "ARG");
        register_alias!("HMR", "ARG");

        register_alias!("AR0", "ARN");

        register_alias!("MEN", "ASN");

        register_alias!("2AS", "ASP");
        register_alias!("ASA", "ASP");
        register_alias!("ASB", "ASP");
        register_alias!("ASK", "ASP");
        register_alias!("ASL", "ASP");
        register_alias!("ASQ", "ASP");
        register_alias!("BHD", "ASP");
        register_alias!("DAS", "ASP");
        register_alias!("DSP", "ASP");
        register_alias!("IAS", "ASP");

        register_alias!("BCS", "CYS");
        register_alias!("BUC", "CYS");
        register_alias!("C5C", "CYS");
        register_alias!("C6C", "CYS");
        register_alias!("CAS", "CYS");
        register_alias!("CCS", "CYS");
        register_alias!("CEA", "CYS");
        register_alias!("CME", "CYS");
        register_alias!("CSO", "CYS");
        register_alias!("CSP", "CYS");
        register_alias!("CSS", "CYS");
        register_alias!("CSW", "CYS");
        register_alias!("CSX", "CYS");
        register_alias!("CY1", "CYS");
        register_alias!("CY3", "CYS");
        register_alias!("CYG", "CYS");
        register_alias!("CYQ", "CYS");
        register_alias!("DCY", "CYS");
        register_alias!("EFC", "CYS");
        register_alias!("OCS", "CYS");
        register_alias!("PEC", "CYS");
        register_alias!("PR3", "CYS");
        register_alias!("PYX", "CYS");
        register_alias!("SCH", "CYS");
        register_alias!("SCS", "CYS");
        register_alias!("SCY", "CYS");
        register_alias!("SHC", "CYS");
        register_alias!("SMC", "CYS");
        register_alias!("SOC", "CYS");

        register_alias!("5HP", "GLU");
        register_alias!("CGU", "GLU");
        register_alias!("DGL", "GLU");
        register_alias!("GGL", "GLU");
        register_alias!("GMA", "GLU");
        register_alias!("PCA", "GLU");

        register_alias!("GLP", "GLH");

        register_alias!("DGN", "GLN");

        register_alias!("GL3", "GLY");
        register_alias!("GLZ", "GLY");
        register_alias!("GSC", "GLY");
        register_alias!("MPQ", "GLY");
        register_alias!("MSA", "GLY");
        register_alias!("NMC", "GLY");
        register_alias!("SAR", "GLY");

        register_alias!("HIS", "HID");
        register_alias!("3AH", "HID");
        register_alias!("DHI", "HID");
        register_alias!("HIC", "HID");
        register_alias!("MHS", "HID");
        register_alias!("NEM", "HID");
        register_alias!("NEP", "HID");

        register_alias!("DIL", "ILE");
        register_alias!("IIL", "ILE");

        register_alias!("BUG", "LEU");
        register_alias!("CLE", "LEU");
        register_alias!("DLE", "LEU");
        register_alias!("MK8", "LEU");
        register_alias!("MLE", "LEU");
        register_alias!("NLE", "LEU");
        register_alias!("NLN", "LEU");
        register_alias!("NLP", "LEU");

        register_alias!("5OW", "LYS");
        register_alias!("ALY", "LYS");
        register_alias!("DLY", "LYS");
        register_alias!("KCX", "LYS");
        register_alias!("LLP", "LYS");
        register_alias!("LLY", "LYS");
        register_alias!("LYM", "LYS");
        register_alias!("LYZ", "LYS");
        register_alias!("SHR", "LYS");
        register_alias!("TRG", "LYS");

        register_alias!("CXM", "MET");
        register_alias!("FME", "MET");
        register_alias!("MSE", "MET");
        register_alias!("OMT", "MET");

        register_alias!("DAH", "PHE");
        register_alias!("DPN", "PHE");
        register_alias!("HPQ", "PHE");
        register_alias!("PHI", "PHE");
        register_alias!("PHL", "PHE");

        register_alias!("DPR", "PRO");
        register_alias!("HYP", "PRO");

        register_alias!("DSN", "SER");
        register_alias!("MIS", "SER");
        register_alias!("OAS", "SER");
        register_alias!("SAC", "SER");
        register_alias!("SEL", "SER");
        register_alias!("SEP", "SER");
        register_alias!("SET", "SER");
        register_alias!("SVA", "SER");

        register_alias!("ALO", "THR");
        register_alias!("BMT", "THR");
        register_alias!("DTH", "THR");
        register_alias!("TPO", "THR");

        register_alias!("DTR", "TRP");
        register_alias!("HTR", "TRP");
        register_alias!("LTR", "TRP");
        register_alias!("TPL", "TRP");
        register_alias!("TRO", "TRP");

        register_alias!("DTY", "TYR");
        register_alias!("IYR", "TYR");
        register_alias!("PAQ", "TYR");
        register_alias!("PTR", "TYR");
        register_alias!("STY", "TYR");
        register_alias!("TYB", "TYR");
        register_alias!("TYI", "TYR");
        register_alias!("TYQ", "TYR");
        register_alias!("TYS", "TYR");
        register_alias!("TYY", "TYR");

        register_alias!("APP", "ASH");

        register_alias!("DIV", "VAL");
        register_alias!("DVA", "VAL");
        register_alias!("MVA", "VAL");

        register_alias!("WAT", "HOH");
        register_alias!("SOL", "HOH");
        register_alias!("TIP", "HOH");
        register_alias!("TIP3", "HOH");
        register_alias!("TP3", "HOH");
        register_alias!("SPC", "HOH");

        let self_register_list = [
            "ALA", "ARN", "ARG", "ASH", "ASN", "ASP", "CYM", "CYS", "CYX", "GLH", "GLN", "GLU",
            "GLY", "HID", "HIE", "HIP", "ILE", "LEU", "LYN", "LYS", "MET", "PHE", "PRO", "SER",
            "THR", "TRP", "TYM", "TYR", "VAL", "A", "C", "DA", "DC", "DG", "DI", "DT", "G", "I",
            "U", "HOH",
        ];

        for name in self_register_list {
            alias_map.insert(name.to_string(), name.to_string());
        }

        Self {
            alias_map,
            standard_map,
        }
    }

    pub fn resolve_name<'a>(&'a self, name: &'a str) -> &'a str {
        self.alias_map.get(name).map(|s| s.as_str()).unwrap_or(name)
    }

    pub fn map_to_standard(&self, name: &str) -> Option<StandardResidue> {
        self.standard_map.get(name).copied()
    }

    pub fn add_alias(&mut self, alias: impl Into<String>, canonical: impl Into<String>) {
        self.alias_map.insert(alias.into(), canonical.into());
    }

    pub fn classify_residue(&self, raw_name: &str) -> (String, Option<StandardResidue>) {
        let canonical = self.resolve_name(raw_name);
        let standard = self.map_to_standard(canonical);
        (canonical.to_string(), standard)
    }
}

impl Default for IoContext {
    fn default() -> Self {
        Self::new_default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::types::StandardResidue;

    #[test]
    fn io_context_new_default_creates_context_with_mappings() {
        let context = IoContext::new_default();

        assert!(context.alias_map.contains_key("ALA"));
        assert!(context.standard_map.contains_key("ALA"));
        assert!(context.alias_map.contains_key("HOH"));
        assert!(context.standard_map.contains_key("HOH"));

        assert!(context.alias_map.contains_key("AIB"));
        assert!(context.alias_map.contains_key("WAT"));
    }

    #[test]
    fn io_context_default_creates_same_as_new_default() {
        let context1 = IoContext::new_default();
        let context2 = IoContext::default();

        assert_eq!(context1.alias_map.len(), context2.alias_map.len());
        assert_eq!(context1.standard_map.len(), context2.standard_map.len());

        assert_eq!(context1.alias_map.get("ALA"), context2.alias_map.get("ALA"));
        assert_eq!(
            context1.standard_map.get("ALA"),
            context2.standard_map.get("ALA")
        );
    }

    #[test]
    fn io_context_clone_creates_identical_copy() {
        let context = IoContext::new_default();
        let cloned = context.clone();

        assert_eq!(context.alias_map, cloned.alias_map);
        assert_eq!(context.standard_map, cloned.standard_map);
    }

    #[test]
    fn io_context_debug_formats_correctly() {
        let context = IoContext::new_default();
        let debug_str = format!("{:?}", context);

        assert!(debug_str.contains("IoContext"));
        assert!(debug_str.contains("alias_map"));
        assert!(debug_str.contains("standard_map"));
    }

    #[test]
    fn resolve_name_returns_canonical_for_standard_residue() {
        let context = IoContext::new_default();

        assert_eq!(context.resolve_name("ALA"), "ALA");
        assert_eq!(context.resolve_name("GLY"), "GLY");
        assert_eq!(context.resolve_name("HOH"), "HOH");
    }

    #[test]
    fn resolve_name_returns_canonical_for_alias() {
        let context = IoContext::new_default();

        assert_eq!(context.resolve_name("AIB"), "ALA");
        assert_eq!(context.resolve_name("WAT"), "HOH");
        assert_eq!(context.resolve_name("SOL"), "HOH");
        assert_eq!(context.resolve_name("DAL"), "ALA");
    }

    #[test]
    fn resolve_name_returns_original_for_unknown_name() {
        let context = IoContext::new_default();

        assert_eq!(context.resolve_name("UNKNOWN"), "UNKNOWN");
        assert_eq!(context.resolve_name("XYZ123"), "XYZ123");
    }

    #[test]
    fn map_to_standard_returns_correct_enum_for_standard_residues() {
        let context = IoContext::new_default();

        assert_eq!(context.map_to_standard("ALA"), Some(StandardResidue::ALA));
        assert_eq!(context.map_to_standard("GLY"), Some(StandardResidue::GLY));
        assert_eq!(context.map_to_standard("ARG"), Some(StandardResidue::ARG));
        assert_eq!(context.map_to_standard("HOH"), Some(StandardResidue::HOH));
        assert_eq!(context.map_to_standard("DA"), Some(StandardResidue::DA));
        assert_eq!(context.map_to_standard("A"), Some(StandardResidue::A));
    }

    #[test]
    fn map_to_standard_returns_none_for_aliases() {
        let context = IoContext::new_default();

        assert_eq!(context.map_to_standard("AIB"), None);
        assert_eq!(context.map_to_standard("WAT"), None);
        assert_eq!(context.map_to_standard("ARN"), Some(StandardResidue::ARG));
    }

    #[test]
    fn map_to_standard_returns_none_for_unknown_names() {
        let context = IoContext::new_default();

        assert_eq!(context.map_to_standard("UNKNOWN"), None);
        assert_eq!(context.map_to_standard("XYZ"), None);
    }

    #[test]
    fn add_alias_adds_new_alias_mapping() {
        let mut context = IoContext::new_default();

        context.add_alias("TEST_ALIAS", "ALA");

        assert_eq!(context.resolve_name("TEST_ALIAS"), "ALA");
        assert_eq!(context.map_to_standard("TEST_ALIAS"), None);
    }

    #[test]
    fn add_alias_overwrites_existing_alias() {
        let mut context = IoContext::new_default();

        assert_eq!(context.resolve_name("AIB"), "ALA");

        context.add_alias("AIB", "GLY");

        assert_eq!(context.resolve_name("AIB"), "GLY");
    }

    #[test]
    fn add_alias_with_string_types() {
        let mut context = IoContext::new_default();

        context.add_alias("STR_ALIAS", "GLY");
        assert_eq!(context.resolve_name("STR_ALIAS"), "GLY");

        context.add_alias(String::from("OWNED_ALIAS"), String::from("ALA"));
        assert_eq!(context.resolve_name("OWNED_ALIAS"), "ALA");
    }

    #[test]
    fn context_handles_case_sensitivity() {
        let context = IoContext::new_default();

        assert_eq!(context.resolve_name("ala"), "ala");
        assert_eq!(context.resolve_name("ALA"), "ALA");
        assert_eq!(context.map_to_standard("ala"), None);
        assert_eq!(context.map_to_standard("ALA"), Some(StandardResidue::ALA));
    }

    #[test]
    fn classify_residue_returns_canonical_and_standard() {
        let context = IoContext::new_default();

        let (name, standard) = context.classify_residue("WAT");
        assert_eq!(name, "HOH");
        assert_eq!(standard, Some(StandardResidue::HOH));

        let (name, standard) = context.classify_residue("ALA");
        assert_eq!(name, "ALA");
        assert_eq!(standard, Some(StandardResidue::ALA));
    }

    #[test]
    fn classify_residue_handles_unknowns() {
        let context = IoContext::new_default();

        let (name, standard) = context.classify_residue("LIG");
        assert_eq!(name, "LIG");
        assert!(standard.is_none());
    }
}
