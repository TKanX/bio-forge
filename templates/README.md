# Template Catalog

This catalog lists all residue templates available in BioForge, with their usage and pH ranges for protonation states.

## Protein Templates

These templates cover all standard amino acids plus their titration-dependent labels. Alternate names (ASH, GLH, LYN, etc.) describe the protonation state that the hydrogenation pass will swap in at the given pH.

| Protonated Residue Name | Standard Residue Name | Charge | pH Range   | Description                                                                              |
| ----------------------- | --------------------- | ------ | ---------- | ---------------------------------------------------------------------------------------- |
| **"ALA"**               | `ALA`                 | 0      | (−∞, +∞)   | Alanine with a compact methyl side chain used as a neutral hydrophobic scaffold.         |
| **"ARG"**               | `ARG`                 | +1     | (−∞, 12.5] | Arginine in its guanidinium cation form; dominant across physiological pH.               |
| **"ARN"**               | `ARG`                 | 0      | (12.5, +∞) | Neutral arginine reached after guanidinium deprotonation at very high pH.                |
| **"ASN"**               | `ASN`                 | 0      | (−∞, +∞)   | Asparagine with a polar carboxamide side chain.                                          |
| **"ASP"**               | `ASP`                 | −1     | [3.9, +∞)  | Aspartate with a deprotonated carboxylate group.                                         |
| **"ASH"**               | `ASP`                 | 0      | (−∞, 3.9)  | Aspartic acid with a protonated carboxyl, favored in strongly acidic media.              |
| **"CYS"**               | `CYS`                 | 0      | (−∞, 8.3]  | Cysteine with a neutral thiol, used below the thiol pKₐ.                                 |
| **"CYM"**               | `CYS`                 | −1     | (8.3, +∞)  | Cysteine thiolate after sulfur deprotonation.                                            |
| **"CYX"**               | `CYS`                 | 0      | (−∞, +∞)   | Disulfide-bonded cystine label assigned whenever SG–SG distance indicates a bridge.      |
| **"GLN"**               | `GLN`                 | 0      | (−∞, +∞)   | Glutamine with a polar carboxamide side chain.                                           |
| **"GLU"**               | `GLU`                 | −1     | [4.2, +∞)  | Glutamate with a deprotonated γ-carboxylate.                                             |
| **"GLH"**               | `GLU`                 | 0      | (−∞, 4.2)  | Glutamic acid with a protonated γ-carboxyl.                                              |
| **"GLY"**               | `GLY`                 | 0      | (−∞, +∞)   | Glycine, the smallest amino acid lacking a side chain carbon.                            |
| **"HID"**               | `HIS`                 | 0      | [6.0, +∞)  | Neutral histidine tautomer protonated on ND1 (δ) as selected by the configured strategy. |
| **"HIE"**               | `HIS`                 | 0      | [6.0, +∞)  | Neutral histidine tautomer protonated on NE2 (ε) as selected by the configured strategy. |
| **"HIP"**               | `HIS`                 | +1     | (−∞, 6.0)  | Histidine imidazolium carrying two protons in acidic conditions.                         |
| **"ILE"**               | `ILE`                 | 0      | (−∞, +∞)   | Isoleucine with a β-branched hydrophobic side chain.                                     |
| **"LEU"**               | `LEU`                 | 0      | (−∞, +∞)   | Leucine with a γ-branched hydrophobic side chain.                                        |
| **"LYS"**               | `LYS`                 | +1     | (−∞, 10.5] | Lysine with a protonated ε-ammonium group.                                               |
| **"LYN"**               | `LYS`                 | 0      | (10.5, +∞) | Neutral lysine formed after ε-amino deprotonation at very high pH.                       |
| **"MET"**               | `MET`                 | 0      | (−∞, +∞)   | Methionine featuring a thioether side chain.                                             |
| **"PHE"**               | `PHE`                 | 0      | (−∞, +∞)   | Phenylalanine with an aromatic phenyl ring.                                              |
| **"PRO"**               | `PRO`                 | 0      | (−∞, +∞)   | Proline with a pyrrolidine ring that constrains backbone φ.                              |
| **"SER"**               | `SER`                 | 0      | (−∞, +∞)   | Serine with a hydroxymethyl side chain.                                                  |
| **"THR"**               | `THR`                 | 0      | (−∞, +∞)   | Threonine with a β-branched hydroxyethyl group.                                          |
| **"TRP"**               | `TRP`                 | 0      | (−∞, +∞)   | Tryptophan bearing an indole aromatic system.                                            |
| **"TYR"**               | `TYR`                 | 0      | (−∞, 10.0] | Tyrosine with an intact phenolic OH.                                                     |
| **"TYM"**               | `TYR`                 | −1     | (10.0, +∞) | Tyrosinate, the phenolate form of tyrosine favored under strongly basic conditions.      |
| **"VAL"**               | `VAL`                 | 0      | (−∞, +∞)   | Valine with a branched hydrophobic isopropyl side chain.                                 |

_Neutral histidines (`HID` / `HIE`) share the same pH interval; BioForge chooses between them using `HydroConfig.his_strategy` (default: hydrogen-bond network analysis)._

## Nucleic Acid Templates

RNA entries use ribose sugars (A, C, G, U, I) while DNA entries carry the `D` prefix. All nucleotides carry a −1 charge from the phosphate and remain valid at any pH because the hydrogenation pass does not titrate nucleic acids.

| Protonated Residue Name | Standard Residue Name | Charge | pH Range | Description                                                     |
| ----------------------- | --------------------- | ------ | -------- | --------------------------------------------------------------- |
| **"A"**                 | `A`                   | −1     | (−∞, +∞) | RNA adenosine residue with adenine base and 2′-hydroxyl ribose. |
| **"C"**                 | `C`                   | −1     | (−∞, +∞) | RNA cytidine residue carrying cytosine and ribose.              |
| **"G"**                 | `G`                   | −1     | (−∞, +∞) | RNA guanosine residue with guanine and ribose.                  |
| **"U"**                 | `U`                   | −1     | (−∞, +∞) | RNA uridine residue with uracil and ribose.                     |
| **"I"**                 | `I`                   | −1     | (−∞, +∞) | RNA inosine residue (hypoxanthine base) commonly used in tRNA.  |
| **"DA"**                | `DA`                  | −1     | (−∞, +∞) | DNA deoxyadenosine residue lacking the ribose 2′-hydroxyl.      |
| **"DC"**                | `DC`                  | −1     | (−∞, +∞) | DNA deoxycytidine residue.                                      |
| **"DG"**                | `DG`                  | −1     | (−∞, +∞) | DNA deoxyguanosine residue.                                     |
| **"DT"**                | `DT`                  | −1     | (−∞, +∞) | DNA thymidine residue.                                          |
| **"DI"**                | `DI`                  | −1     | (−∞, +∞) | DNA deoxyinosine residue for inosine-containing duplexes.       |

## Solvent Templates

| Protonated Residue Name | Standard Residue Name | Charge | pH Range | Description                                       |
| ----------------------- | --------------------- | ------ | -------- | ------------------------------------------------- |
| **"HOH"**               | `HOH`                 | 0      | (−∞, +∞) | Water molecule used for explicit solvent records. |
