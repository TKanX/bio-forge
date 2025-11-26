# End-to-End Preparation Example (single command per stage)

![1CRN Structure](./images/structure.png)

## Inputs

- **Structure**: `1CRN.pdb` — Crambin crystal structure.

## Commands

```bash
bioforge clean    -i 1CRN.pdb -o 1CRN.clean.pdb --water --ions
bioforge repair   -i 1CRN.clean.pdb -o 1CRN.repair.pdb
bioforge hydro    -i 1CRN.repair.pdb -o 1CRN.h.pdb --ph 7.4
bioforge solvate  -i 1CRN.h.pdb -o 1CRN.solvated.pdb --neutralize
bioforge topology -i 1CRN.solvated.pdb -o 1CRN.topo.cif --out-format mmcif
```

## Results

- Sequentially cleans, repairs, protonates, solvates, and builds the topology for `1CRN.pdb`.
- Produces the final topology file `1CRN.topo.cif`.
