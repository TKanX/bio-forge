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
