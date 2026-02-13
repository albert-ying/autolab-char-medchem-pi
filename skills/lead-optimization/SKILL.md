---
name: lead-optimization
description: Structure-activity relationship analysis, matched molecular pairs, property optimization for drug candidates.
metadata:
    skill-author: Albert Ying
---

# Lead optimization

## When to use

- Analyzing SAR trends across a compound series
- Identifying matched molecular pairs for property optimization
- Balancing potency vs ADMET properties
- Prioritizing compounds for synthesis

## SAR analysis

```python
import datamol as dm
import pandas as pd
from rdkit.Chem import AllChem, Descriptors

# Compute properties across a series
df["MW"] = df["smiles"].apply(lambda s: Descriptors.MolWt(dm.to_mol(s)))
df["LogP"] = df["smiles"].apply(lambda s: Descriptors.MolLogP(dm.to_mol(s)))
df["TPSA"] = df["smiles"].apply(lambda s: Descriptors.TPSA(dm.to_mol(s)))
df["HBA"] = df["smiles"].apply(lambda s: Descriptors.NumHAcceptors(dm.to_mol(s)))
df["HBD"] = df["smiles"].apply(lambda s: Descriptors.NumHDonors(dm.to_mol(s)))

# Activity cliffs: high similarity, large activity difference
from rdkit.Chem import DataStructs
fps = [AllChem.GetMorganFingerprintAsBitVect(dm.to_mol(s), 2) for s in df["smiles"]]
for i in range(len(fps)):
    for j in range(i+1, len(fps)):
        sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
        if sim > 0.8:
            delta = abs(df.iloc[i]["pIC50"] - df.iloc[j]["pIC50"])
            if delta > 1.0:  # activity cliff
                print(f"Cliff: {i} vs {j}, sim={sim:.2f}, delta={delta:.1f}")
```

## Key principles

- Always check Lipinski/Veber before advancing a series.
- Log activity (pIC50) is preferred over raw IC50 for SAR analysis.
- Track MW, LogP, TPSA, and rotatable bonds across optimization rounds.
- Flag compounds with MW > 450 or LogP > 4 early — hard to optimize later.
