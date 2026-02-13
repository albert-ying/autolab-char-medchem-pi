---
name: virtual-screening
description: Computational virtual screening with RDKit and datamol. Molecular filtering, fingerprint similarity, substructure search, docking preparation.
metadata:
    skill-author: Albert Ying
---

# Virtual screening

## When to use

- Filtering compound libraries by drug-likeness
- Similarity or substructure searches against a query
- Preparing molecules for docking
- Building focused screening libraries

## Library filtering

```python
import datamol as dm
from rdkit.Chem import Descriptors, FilterCatalog

# Load and standardize
mols = dm.read_sdf("library.sdf")
mols = [dm.standardize_mol(m) for m in mols if m is not None]

# Lipinski filter
def passes_lipinski(mol):
    return (Descriptors.MolWt(mol) <= 500
            and Descriptors.MolLogP(mol) <= 5
            and Descriptors.NumHDonors(mol) <= 5
            and Descriptors.NumHAcceptors(mol) <= 10)

hits = [m for m in mols if passes_lipinski(m)]

# PAINS filter
params = FilterCatalog.FilterCatalogParams()
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
catalog = FilterCatalog.FilterCatalog(params)
clean = [m for m in hits if not catalog.HasMatch(m)]
```

## Similarity search

```python
from rdkit.Chem import AllChem, DataStructs

query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=2048)
results = []
for mol in library:
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    sim = DataStructs.TanimotoSimilarity(query_fp, fp)
    if sim >= 0.7:
        results.append((mol, sim))
results.sort(key=lambda x: -x[1])
```
