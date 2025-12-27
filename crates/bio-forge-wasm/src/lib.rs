//! WebAssembly bindings for the BioForge molecular preparation library.
//!
//! This crate provides JavaScript/TypeScript bindings for BioForge's core functionality,
//! enabling molecular structure manipulation directly in the browser or other WASM-capable environments.

use bio_forge::io::{
    IoContext, read_mmcif_structure, read_mol2_template, read_pdb_structure, write_mmcif_structure,
    write_mmcif_topology, write_pdb_structure, write_pdb_topology,
};
use bio_forge::ops::{
    Anion as CoreAnion, Cation as CoreCation, CleanConfig as CoreCleanConfig,
    HisStrategy as CoreHisStrategy, HydroConfig as CoreHydroConfig,
    SolvateConfig as CoreSolvateConfig, TopologyBuilder, Transform as CoreTransform,
    add_hydrogens as core_add_hydrogens, clean_structure as core_clean_structure,
    repair_structure as core_repair_structure, solvate_structure as core_solvate_structure,
};
use bio_forge::{
    Chain as CoreChain, ResidueCategory, StandardResidue, Structure as CoreStructure,
    Template as CoreTemplate, Topology as CoreTopology,
};
use serde::{Deserialize, Serialize};
use std::io::{BufWriter, Cursor};
use tsify::Tsify;
use wasm_bindgen::prelude::*;

// ============================================================================
// Initialization
// ============================================================================

/// Initializes panic hook for better error messages in browser console.
///
/// This function is automatically called when the WASM module is loaded.
/// It sets up `console_error_panic_hook` to provide readable panic messages
/// in the browser's developer console instead of cryptic WASM error codes.
#[wasm_bindgen(start)]
pub fn init() {
    console_error_panic_hook::set_once();
}

// ============================================================================
// Error Handling
// ============================================================================

/// Converts any displayable error into a JavaScript error.
fn to_js_error<E: std::fmt::Display>(e: E) -> JsError {
    JsError::new(&e.to_string())
}

// ============================================================================
// Configuration: CleanConfig
// ============================================================================

/// Configuration for structure cleaning operations.
#[derive(Debug, Clone, Default, Serialize, Deserialize, Tsify)]
#[tsify(into_wasm_abi, from_wasm_abi)]
#[serde(rename_all = "camelCase")]
pub struct CleanConfig {
    /// Remove water molecules (HOH). Default: `false`
    #[serde(default)]
    pub remove_water: bool,
    /// Remove ion residues. Default: `false`
    #[serde(default)]
    pub remove_ions: bool,
    /// Remove hydrogen atoms. Default: `false`
    #[serde(default)]
    pub remove_hydrogens: bool,
    /// Remove hetero residues (ligands). Default: `false`
    #[serde(default)]
    pub remove_hetero: bool,
    /// Specific residue names to remove. Default: `[]`
    #[serde(default)]
    pub remove_residue_names: Vec<String>,
    /// Specific residue names to keep (overrides other rules). Default: `[]`
    #[serde(default)]
    pub keep_residue_names: Vec<String>,
}

impl From<CleanConfig> for CoreCleanConfig {
    fn from(cfg: CleanConfig) -> Self {
        CoreCleanConfig {
            remove_water: cfg.remove_water,
            remove_ions: cfg.remove_ions,
            remove_hydrogens: cfg.remove_hydrogens,
            remove_hetero: cfg.remove_hetero,
            remove_residue_names: cfg.remove_residue_names.into_iter().collect(),
            keep_residue_names: cfg.keep_residue_names.into_iter().collect(),
        }
    }
}

// ============================================================================
// Configuration: HydroConfig
// ============================================================================

/// Configuration for hydrogen addition.
#[derive(Debug, Clone, Serialize, Deserialize, Tsify)]
#[tsify(into_wasm_abi, from_wasm_abi)]
#[serde(rename_all = "camelCase")]
pub struct HydroConfig {
    /// Target pH for protonation state decisions. Default: `undefined` (neutral)
    #[serde(default)]
    pub target_ph: Option<f64>,
    /// Remove existing hydrogens before adding new ones. Default: `true`
    #[serde(default = "default_true")]
    pub remove_existing_h: bool,
    /// Histidine tautomer strategy: `"hid"`, `"hie"`, `"random"`, `"network"`. Default: `"network"`
    #[tsify(type = "\"hid\" | \"hie\" | \"random\" | \"network\"")]
    #[serde(default = "default_his_strategy")]
    pub his_strategy: String,
}

fn default_true() -> bool {
    true
}

fn default_his_strategy() -> String {
    "network".to_string()
}

impl Default for HydroConfig {
    fn default() -> Self {
        Self {
            target_ph: None,
            remove_existing_h: true,
            his_strategy: default_his_strategy(),
        }
    }
}

impl From<HydroConfig> for CoreHydroConfig {
    fn from(cfg: HydroConfig) -> Self {
        let his_strategy = match cfg.his_strategy.to_lowercase().as_str() {
            "hid" => CoreHisStrategy::DirectHID,
            "hie" => CoreHisStrategy::DirectHIE,
            "random" => CoreHisStrategy::Random,
            _ => CoreHisStrategy::HbNetwork,
        };
        CoreHydroConfig {
            target_ph: cfg.target_ph,
            remove_existing_h: cfg.remove_existing_h,
            his_strategy,
        }
    }
}

// ============================================================================
// Configuration: SolvateConfig
// ============================================================================

/// Configuration for solvation.
#[derive(Debug, Clone, Serialize, Deserialize, Tsify)]
#[tsify(into_wasm_abi, from_wasm_abi)]
#[serde(rename_all = "camelCase")]
pub struct SolvateConfig {
    /// Margin around solute (Å). Default: `10.0`
    #[serde(default = "default_margin")]
    pub margin: f64,
    /// Water grid spacing (Å). Default: `3.1`
    #[serde(default = "default_water_spacing")]
    pub water_spacing: f64,
    /// Minimum solvent-solute distance (Å). Default: `2.4`
    #[serde(default = "default_vdw_cutoff")]
    pub vdw_cutoff: f64,
    /// Remove existing solvent before solvating. Default: `true`
    #[serde(default = "default_true")]
    pub remove_existing: bool,
    /// Cation species: `"Na"`, `"K"`, `"Mg"`, `"Ca"`, `"Li"`, `"Zn"`. Default: `["Na"]`
    #[tsify(type = "Array<\"Na\" | \"K\" | \"Mg\" | \"Ca\" | \"Li\" | \"Zn\">")]
    #[serde(default = "default_cations")]
    pub cations: Vec<String>,
    /// Anion species: `"Cl"`, `"Br"`, `"I"`, `"F"`. Default: `["Cl"]`
    #[tsify(type = "Array<\"Cl\" | \"Br\" | \"I\" | \"F\">")]
    #[serde(default = "default_anions")]
    pub anions: Vec<String>,
    /// Target net charge after solvation. Default: `0`
    #[serde(default)]
    pub target_charge: i32,
    /// Random seed for reproducible placement.
    #[serde(default)]
    pub rng_seed: Option<u64>,
}

fn default_margin() -> f64 {
    10.0
}
fn default_water_spacing() -> f64 {
    3.1
}
fn default_vdw_cutoff() -> f64 {
    2.4
}
fn default_cations() -> Vec<String> {
    vec!["Na".to_string()]
}
fn default_anions() -> Vec<String> {
    vec!["Cl".to_string()]
}

impl Default for SolvateConfig {
    fn default() -> Self {
        Self {
            margin: default_margin(),
            water_spacing: default_water_spacing(),
            vdw_cutoff: default_vdw_cutoff(),
            remove_existing: true,
            cations: default_cations(),
            anions: default_anions(),
            target_charge: 0,
            rng_seed: None,
        }
    }
}

impl From<SolvateConfig> for CoreSolvateConfig {
    fn from(cfg: SolvateConfig) -> Self {
        let cations = cfg
            .cations
            .iter()
            .filter_map(|s| match s.to_uppercase().as_str() {
                "NA" => Some(CoreCation::Na),
                "K" => Some(CoreCation::K),
                "MG" => Some(CoreCation::Mg),
                "CA" => Some(CoreCation::Ca),
                "LI" => Some(CoreCation::Li),
                "ZN" => Some(CoreCation::Zn),
                _ => None,
            })
            .collect();
        let anions = cfg
            .anions
            .iter()
            .filter_map(|s| match s.to_uppercase().as_str() {
                "CL" => Some(CoreAnion::Cl),
                "BR" => Some(CoreAnion::Br),
                "I" => Some(CoreAnion::I),
                "F" => Some(CoreAnion::F),
                _ => None,
            })
            .collect();
        CoreSolvateConfig {
            margin: cfg.margin,
            water_spacing: cfg.water_spacing,
            vdw_cutoff: cfg.vdw_cutoff,
            remove_existing: cfg.remove_existing,
            cations,
            anions,
            target_charge: cfg.target_charge,
            rng_seed: cfg.rng_seed,
        }
    }
}

// ============================================================================
// Configuration: TopologyConfig
// ============================================================================

/// Configuration for topology building.
#[derive(Debug, Clone, Serialize, Deserialize, Tsify)]
#[tsify(into_wasm_abi, from_wasm_abi)]
#[serde(rename_all = "camelCase")]
pub struct TopologyConfig {
    /// Disulfide bond cutoff (Å). Default: `2.2`
    #[serde(default = "default_disulfide_cutoff")]
    pub disulfide_cutoff: f64,
}

fn default_disulfide_cutoff() -> f64 {
    2.2
}

impl Default for TopologyConfig {
    fn default() -> Self {
        Self {
            disulfide_cutoff: default_disulfide_cutoff(),
        }
    }
}

// ============================================================================
// Data Structures: Info
// ============================================================================

/// Comprehensive structure information.
#[derive(Debug, Serialize, Deserialize, Tsify)]
#[tsify(into_wasm_abi)]
#[serde(rename_all = "camelCase")]
pub struct StructureInfo {
    /// Number of chains.
    pub chain_count: usize,
    /// Total number of residues.
    pub residue_count: usize,
    /// Total number of atoms.
    pub atom_count: usize,
    /// Box dimensions `[a, b, c]` in Å, if present.
    pub box_lengths: Option<[f64; 3]>,
    /// Box angles `[α, β, γ]` in degrees, if present.
    pub box_angles: Option<[f64; 3]>,
    /// Per-chain information.
    pub chains: Vec<ChainInfo>,
}

/// Per-chain information.
#[derive(Debug, Serialize, Deserialize, Tsify)]
#[tsify(into_wasm_abi)]
#[serde(rename_all = "camelCase")]
pub struct ChainInfo {
    /// Chain identifier.
    pub id: String,
    /// Number of residues in this chain.
    pub residue_count: usize,
    /// Number of atoms in this chain.
    pub atom_count: usize,
    /// Polymer classifications present in this chain.
    /// Possible values: `"protein"`, `"nucleic"`, `"solvent"`, `"hetero"`.
    /// An empty array indicates an empty chain.
    #[tsify(type = "Array<\"protein\" | \"nucleic\" | \"solvent\" | \"hetero\">")]
    pub polymer_types: Vec<String>,
}

// ============================================================================
// Class: Template
// ============================================================================

/// A residue or ligand template for topology building.
#[wasm_bindgen]
pub struct Template {
    inner: CoreTemplate,
}

#[wasm_bindgen]
impl Template {
    /// Parses a MOL2 string into a Template.
    #[wasm_bindgen(js_name = fromMol2)]
    pub fn from_mol2(content: &str) -> Result<Template, JsError> {
        let cursor = Cursor::new(content);
        let inner = read_mol2_template(cursor).map_err(to_js_error)?;
        Ok(Template { inner })
    }

    /// Parses a MOL2 byte array into a Template.
    #[wasm_bindgen(js_name = fromMol2Bytes)]
    pub fn from_mol2_bytes(content: &[u8]) -> Result<Template, JsError> {
        let cursor = Cursor::new(content);
        let inner = read_mol2_template(cursor).map_err(to_js_error)?;
        Ok(Template { inner })
    }

    /// Returns the template name (residue code).
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        self.inner.name.clone()
    }
}

// ============================================================================
// Class: Topology
// ============================================================================

/// A structure with bond connectivity information.
#[wasm_bindgen]
pub struct Topology {
    inner: CoreTopology,
}

#[wasm_bindgen]
impl Topology {
    /// Creates a Topology from a Structure.
    #[wasm_bindgen(js_name = fromStructure)]
    pub fn from_structure(
        structure: &Structure,
        config: Option<TopologyConfig>,
        templates: Option<Vec<Template>>,
    ) -> Result<Topology, JsError> {
        let cfg = config.unwrap_or_default();
        let mut builder = TopologyBuilder::new().disulfide_cutoff(cfg.disulfide_cutoff);

        if let Some(tpls) = templates {
            for t in tpls {
                builder = builder.add_hetero_template(t.inner);
            }
        }

        let inner = builder
            .build(structure.inner.clone())
            .map_err(to_js_error)?;
        Ok(Topology { inner })
    }

    /// Returns a copy of the underlying structure.
    #[wasm_bindgen(getter)]
    pub fn structure(&self) -> Structure {
        Structure {
            inner: self.inner.structure().clone(),
        }
    }

    /// Returns the number of bonds.
    #[wasm_bindgen(js_name = bondCount, getter)]
    pub fn bond_count(&self) -> usize {
        self.inner.bonds().len()
    }

    /// Serializes to PDB format with CONECT records.
    #[wasm_bindgen(js_name = toPdb)]
    pub fn to_pdb(&self) -> Result<String, JsError> {
        let mut buf = Vec::new();
        write_pdb_topology(BufWriter::new(&mut buf), &self.inner).map_err(to_js_error)?;
        String::from_utf8(buf).map_err(to_js_error)
    }

    /// Serializes to PDB format as bytes.
    #[wasm_bindgen(js_name = toPdbBytes)]
    pub fn to_pdb_bytes(&self) -> Result<Vec<u8>, JsError> {
        let mut buf = Vec::new();
        write_pdb_topology(BufWriter::new(&mut buf), &self.inner).map_err(to_js_error)?;
        Ok(buf)
    }

    /// Serializes to mmCIF format with struct_conn records.
    #[wasm_bindgen(js_name = toMmcif)]
    pub fn to_mmcif(&self) -> Result<String, JsError> {
        let mut buf = Vec::new();
        write_mmcif_topology(BufWriter::new(&mut buf), &self.inner).map_err(to_js_error)?;
        String::from_utf8(buf).map_err(to_js_error)
    }

    /// Serializes to mmCIF format as bytes.
    #[wasm_bindgen(js_name = toMmcifBytes)]
    pub fn to_mmcif_bytes(&self) -> Result<Vec<u8>, JsError> {
        let mut buf = Vec::new();
        write_mmcif_topology(BufWriter::new(&mut buf), &self.inner).map_err(to_js_error)?;
        Ok(buf)
    }
}

// ============================================================================
// Class: Structure
// ============================================================================

/// A molecular structure for manipulation and export.
#[wasm_bindgen]
pub struct Structure {
    inner: CoreStructure,
}

#[wasm_bindgen]
impl Structure {
    // -------------------------------------------------------------------------
    // Factory Methods
    // -------------------------------------------------------------------------

    /// Parses a PDB string into a Structure.
    #[wasm_bindgen(js_name = fromPdb)]
    pub fn from_pdb(content: &str) -> Result<Structure, JsError> {
        let ctx = IoContext::new_default();
        let cursor = Cursor::new(content);
        let inner = read_pdb_structure(cursor, &ctx).map_err(to_js_error)?;
        Ok(Structure { inner })
    }

    /// Parses a PDB byte array into a Structure.
    #[wasm_bindgen(js_name = fromPdbBytes)]
    pub fn from_pdb_bytes(content: &[u8]) -> Result<Structure, JsError> {
        let ctx = IoContext::new_default();
        let cursor = Cursor::new(content);
        let inner = read_pdb_structure(cursor, &ctx).map_err(to_js_error)?;
        Ok(Structure { inner })
    }

    /// Parses an mmCIF string into a Structure.
    #[wasm_bindgen(js_name = fromMmcif)]
    pub fn from_mmcif(content: &str) -> Result<Structure, JsError> {
        let ctx = IoContext::new_default();
        let cursor = Cursor::new(content);
        let inner = read_mmcif_structure(cursor, &ctx).map_err(to_js_error)?;
        Ok(Structure { inner })
    }

    /// Parses an mmCIF byte array into a Structure.
    #[wasm_bindgen(js_name = fromMmcifBytes)]
    pub fn from_mmcif_bytes(content: &[u8]) -> Result<Structure, JsError> {
        let ctx = IoContext::new_default();
        let cursor = Cursor::new(content);
        let inner = read_mmcif_structure(cursor, &ctx).map_err(to_js_error)?;
        Ok(Structure { inner })
    }

    // -------------------------------------------------------------------------
    // Instance Methods
    // -------------------------------------------------------------------------

    /// Creates a deep copy of the structure.
    #[wasm_bindgen(js_name = clone)]
    pub fn clone_structure(&self) -> Structure {
        Structure {
            inner: self.inner.clone(),
        }
    }

    // -------------------------------------------------------------------------
    // Export Methods
    // -------------------------------------------------------------------------

    /// Serializes to PDB format.
    #[wasm_bindgen(js_name = toPdb)]
    pub fn to_pdb(&self) -> Result<String, JsError> {
        let mut buf = Vec::new();
        write_pdb_structure(BufWriter::new(&mut buf), &self.inner).map_err(to_js_error)?;
        String::from_utf8(buf).map_err(to_js_error)
    }

    /// Serializes to PDB format as bytes.
    #[wasm_bindgen(js_name = toPdbBytes)]
    pub fn to_pdb_bytes(&self) -> Result<Vec<u8>, JsError> {
        let mut buf = Vec::new();
        write_pdb_structure(BufWriter::new(&mut buf), &self.inner).map_err(to_js_error)?;
        Ok(buf)
    }

    /// Serializes to mmCIF format.
    #[wasm_bindgen(js_name = toMmcif)]
    pub fn to_mmcif(&self) -> Result<String, JsError> {
        let mut buf = Vec::new();
        write_mmcif_structure(BufWriter::new(&mut buf), &self.inner).map_err(to_js_error)?;
        String::from_utf8(buf).map_err(to_js_error)
    }

    /// Serializes to mmCIF format as bytes.
    #[wasm_bindgen(js_name = toMmcifBytes)]
    pub fn to_mmcif_bytes(&self) -> Result<Vec<u8>, JsError> {
        let mut buf = Vec::new();
        write_mmcif_structure(BufWriter::new(&mut buf), &self.inner).map_err(to_js_error)?;
        Ok(buf)
    }

    /// Builds a Topology from this structure.
    #[wasm_bindgen(js_name = toTopology)]
    pub fn to_topology(
        &self,
        config: Option<TopologyConfig>,
        templates: Option<Vec<Template>>,
    ) -> Result<Topology, JsError> {
        Topology::from_structure(self, config, templates)
    }

    // -------------------------------------------------------------------------
    // Mutation Methods
    // -------------------------------------------------------------------------

    /// Removes unwanted components (water, ions, hydrogens, hetero residues).
    #[wasm_bindgen]
    pub fn clean(&mut self, config: Option<CleanConfig>) -> Result<(), JsError> {
        let cfg = config.unwrap_or_default();
        core_clean_structure(&mut self.inner, &cfg.into()).map_err(to_js_error)
    }

    /// Reconstructs missing atoms from templates.
    #[wasm_bindgen]
    pub fn repair(&mut self) -> Result<(), JsError> {
        core_repair_structure(&mut self.inner).map_err(to_js_error)
    }

    /// Adds hydrogen atoms.
    #[wasm_bindgen(js_name = addHydrogens)]
    pub fn add_hydrogens(&mut self, config: Option<HydroConfig>) -> Result<(), JsError> {
        let cfg = config.unwrap_or_default();
        core_add_hydrogens(&mut self.inner, &cfg.into()).map_err(to_js_error)
    }

    /// Solvates the structure with water and ions.
    #[wasm_bindgen]
    pub fn solvate(&mut self, config: Option<SolvateConfig>) -> Result<(), JsError> {
        let cfg = config.unwrap_or_default();
        core_solvate_structure(&mut self.inner, &cfg.into()).map_err(to_js_error)
    }

    // -------------------------------------------------------------------------
    // Transform Methods
    // -------------------------------------------------------------------------

    /// Translates all atoms by the given vector (x, y, z in Å).
    #[wasm_bindgen]
    pub fn translate(&mut self, x: f64, y: f64, z: f64) {
        CoreTransform::translate(&mut self.inner, x, y, z);
    }

    /// Centers the geometric centroid at the origin.
    #[wasm_bindgen(js_name = centerGeometry)]
    pub fn center_geometry(&mut self) {
        CoreTransform::center_geometry(&mut self.inner, None);
    }

    /// Centers the center of mass at the origin.
    #[wasm_bindgen(js_name = centerMass)]
    pub fn center_mass(&mut self) {
        CoreTransform::center_mass(&mut self.inner, None);
    }

    /// Rotates around the X axis (angle in radians).
    #[wasm_bindgen(js_name = rotateX)]
    pub fn rotate_x(&mut self, radians: f64) {
        CoreTransform::rotate_x(&mut self.inner, radians);
    }

    /// Rotates around the Y axis (angle in radians).
    #[wasm_bindgen(js_name = rotateY)]
    pub fn rotate_y(&mut self, radians: f64) {
        CoreTransform::rotate_y(&mut self.inner, radians);
    }

    /// Rotates around the Z axis (angle in radians).
    #[wasm_bindgen(js_name = rotateZ)]
    pub fn rotate_z(&mut self, radians: f64) {
        CoreTransform::rotate_z(&mut self.inner, radians);
    }

    /// Rotates using Euler angles (XYZ convention, all in radians).
    #[wasm_bindgen(js_name = rotateEuler)]
    pub fn rotate_euler(&mut self, x: f64, y: f64, z: f64) {
        CoreTransform::rotate_euler(&mut self.inner, x, y, z);
    }

    // -------------------------------------------------------------------------
    // Query Methods
    // -------------------------------------------------------------------------

    /// Returns comprehensive structure statistics.
    #[wasm_bindgen]
    pub fn info(&self) -> Result<StructureInfo, JsError> {
        let chains: Vec<ChainInfo> = self
            .inner
            .iter_chains()
            .map(|c| ChainInfo {
                id: c.id.to_string(),
                residue_count: c.residue_count(),
                atom_count: c.atom_count(),
                polymer_types: classify_chain(c),
            })
            .collect();

        let (box_lengths, box_angles) = self
            .inner
            .box_vectors
            .map(|m| {
                let a = [m[0][0], m[0][1], m[0][2]];
                let b = [m[1][0], m[1][1], m[1][2]];
                let c = [m[2][0], m[2][1], m[2][2]];
                let la = vec_norm(a);
                let lb = vec_norm(b);
                let lc = vec_norm(c);
                let alpha = angle_deg(b, c);
                let beta = angle_deg(a, c);
                let gamma = angle_deg(a, b);
                (Some([la, lb, lc]), Some([alpha, beta, gamma]))
            })
            .unwrap_or((None, None));

        Ok(StructureInfo {
            chain_count: self.inner.chain_count(),
            residue_count: self.inner.residue_count(),
            atom_count: self.inner.atom_count(),
            box_lengths,
            box_angles,
            chains,
        })
    }

    /// Returns the number of chains.
    #[wasm_bindgen(js_name = chainCount, getter)]
    pub fn chain_count(&self) -> usize {
        self.inner.chain_count()
    }

    /// Returns the number of residues.
    #[wasm_bindgen(js_name = residueCount, getter)]
    pub fn residue_count(&self) -> usize {
        self.inner.residue_count()
    }

    /// Returns the number of atoms.
    #[wasm_bindgen(js_name = atomCount, getter)]
    pub fn atom_count(&self) -> usize {
        self.inner.atom_count()
    }
}

// ============================================================================
// Helper Functions
// ============================================================================

/// Classifies a chain into one or more polymer types.
fn classify_chain(chain: &CoreChain) -> Vec<String> {
    if chain.is_empty() {
        return vec![];
    }

    let mut protein = false;
    let mut nucleic = false;
    let mut solvent = false;
    let mut hetero = false;

    for res in chain.iter_residues() {
        match res.category {
            ResidueCategory::Standard => {
                if let Some(std) = res.standard_name {
                    if std.is_protein() {
                        protein = true;
                    } else if std.is_nucleic() {
                        nucleic = true;
                    } else if std == StandardResidue::HOH {
                        solvent = true;
                    } else {
                        hetero = true;
                    }
                } else {
                    hetero = true;
                }
            }
            ResidueCategory::Ion => solvent = true,
            ResidueCategory::Hetero => hetero = true,
        }
    }

    let mut types = Vec::new();
    if protein {
        types.push("protein".to_string());
    }
    if nucleic {
        types.push("nucleic".to_string());
    }
    if solvent {
        types.push("solvent".to_string());
    }
    if hetero {
        types.push("hetero".to_string());
    }
    types
}

fn vec_norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn vec_dot(v1: [f64; 3], v2: [f64; 3]) -> f64 {
    v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]
}

fn angle_deg(v1: [f64; 3], v2: [f64; 3]) -> f64 {
    let denom = vec_norm(v1) * vec_norm(v2);
    if denom < f64::EPSILON {
        return 0.0;
    }
    let cos_angle = (vec_dot(v1, v2) / denom).clamp(-1.0, 1.0);
    cos_angle.acos().to_degrees()
}
