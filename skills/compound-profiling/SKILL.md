---
name: compound-profiling
description: Molecular property calculation, ADMET prediction, drug-likeness scoring. RDKit descriptors and medchem filters.
metadata:
    skill-author: Albert Ying
---

# Compound profiling

## When to use

- Computing molecular descriptors for a compound set
- Applying drug-likeness and medchem filters
- ADMET property prediction
- Comparing compound profiles across series

## Property calculation

```python
import datamol as dm
from rdkit.Chem import Descriptors, rdMolDescriptors

mol = dm.to_mol("CCO")
profile = {
    "MW": Descriptors.MolWt(mol),
    "LogP": Descriptors.MolLogP(mol),
    "TPSA": Descriptors.TPSA(mol),
    "HBA": Descriptors.NumHAcceptors(mol),
    "HBD": Descriptors.NumHDonors(mol),
    "RotBonds": Descriptors.NumRotatableBonds(mol),
    "Rings": rdMolDescriptors.CalcNumRings(mol),
    "AromaticRings": rdMolDescriptors.CalcNumAromaticRings(mol),
    "FractionCSP3": rdMolDescriptors.CalcFractionCSP3(mol),
    "QED": Descriptors.qed(mol),
}
```

## Drug-likeness rules

| Rule | Criteria |
|------|----------|
| Lipinski | MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10 |
| Veber | RotBonds ≤ 10, TPSA ≤ 140 |
| Lead-like | MW 250-350, LogP ≤ 3.5, RotBonds ≤ 7 |
| Fragment-like | MW ≤ 300, LogP ≤ 3, HBD ≤ 3, HBA ≤ 3 |

## QED score

Quantitative Estimate of Drug-likeness (0-1). Compounds > 0.5 are typically drug-like. Use `Descriptors.qed(mol)` from RDKit.
