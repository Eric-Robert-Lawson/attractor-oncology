# MU Ask Pro — Query 1 Response
**Date:** 2026-04-03  
**Query:** Candidate 1 solvation structure — LiFSI 1.2 mol/kg in 2-MeTHF:DME 0.60:0.40  
**Tool used:** Ask Pro (1 of 10 queries used)  
**Platform:** SES AI Molecular Universe  

---

```
For a lithium electrolyte consisting of LiFSI salt at 
1.2 mol/kg in a mixed solvent of 2-methyltetrahydrofuran 
(2-MeTHF) and 1,2-dimethoxyethane (DME) at a mole 
fraction ratio of 0.60:0.40 (2-MeTHF:DME):

1. What are the predicted SSIP, CIP, and AGG coordination 
   fractions around Li+ at 25°C?

2. What is the dominant coordination species and its 
   percentage of the total population?

3. How many coordination species have greater than 5% 
   population (the number of significant configs)?

4. What is the predicted ionic conductivity at 25°C 
   and at -20°C?

Please give numerical fractions not qualitative descriptions.
```

---

## MU Response — Full Text

MU stated it cannot give defensible exact predicted numbers for the specific 
formulation (1.2 mol/kg LiFSI in 2-MeTHF:DME = 0.60:0.40) because it does 
not have a live solvation/transport simulation endpoint for that exact 
composition. MU explicitly declined to invent values.

The closest quantitative data retrievable are analog systems in LiFSI ether 
electrolytes, used as anchors only, not as exact outputs for the mixture.

---

## Anchor Data Returned by MU

### Ion Speciation (from Raman proxies, 25°C)

Raman band assignments for LiFSI ether systems:
- 727 cm⁻¹ → SSIP
- 735 cm⁻¹ → CIP  
- 747 cm⁻¹ → AGG

| Species | 1.0 mol/kg LiFSI in DME | 1.0 mol/kg LiFSI in DEE |
|---------|------------------------|------------------------|
| SSIP | 61.1% | 34.5% |
| CIP | 35.4% | 55.3% |
| AGG | 3.5% | 10.2% |

DME endpoint: SSIP-dominated (61.1% SSIP)  
DEE endpoint: CIP-dominated (55.3% CIP)

### Ionic Conductivity

| Condition | 1.0 mol/kg LiFSI in DME | 1.0 mol/kg LiFSI in DEE |
|-----------|------------------------|------------------------|
| 25°C | 16.69 mS/cm | 7.86 mS/cm |
| −20°C | 8.44 mS/cm | 3.14 mS/cm |

---

## MU's Stated Limitations

- Exact SSIP/CIP/AGG fractions for the 2-MeTHF/DME mixture: **not available**
- Dominant species percentage for the mixture: **not available**
- Number of coordination species with >5% population: **not available**
- Conductivities for the exact mixture: **not available**

MU stated the mixture should lie between the DME-rich SSIP-dominant limit 
and the weaker-solvating ether CIP-richer limit, but declined to assign 
specific numerical values without labelling them as unsupported estimates.

MU offered to provide a clearly labelled interpolation estimate using the 
DME and DEE endpoints as bounds.

---

## MU's Position on the Mixture

> "Your electrolyte should lie between the DME-rich SSIP-dominant limit 
> and the weaker-solvating ether CIP-richer limit, but any single numerical 
> values I gave for your exact mixture would be an unsupported estimate, 
> not a validated prediction."

---

## References Cited by MU

1. Pham, T.D., Lee, K.-K. *Simultaneous Stabilization of the Solid/Cathode 
   Electrolyte Interface in Lithium Metal Batteries by a New Weakly Solvating 
   Electrolyte*. Small, 2021.  
   DOI: [10.1002/smll.202100133](https://doi.org/10.1002/smll.202100133)

2. Li, Z., Chen, X., Li, W., et al. *High-Concentrated Binary-Salt Ether 
   Electrolytes for High-Voltage Lithium Metal Batteries with Ni-Rich Cathode*. 
   ACS Applied Materials & Interfaces, 2024.  
   DOI: [10.1021/acsami.4c06491](https://doi.org/10.1021/acsami.4c06491)

---

## What This Response Means for the SCE Framework

### What MU confirmed

- LiFSI is recognised as a valid salt in the platform
- SSIP/CIP/AGG speciation framework is the correct output format
- The mixture will sit between the DME endpoint (SSIP-dominant) and a 
  weaker ether endpoint (CIP-richer)
- This directional prediction is consistent with the SCE framework's 
  prediction that the 2-MeTHF:DME mixture at 0.60:0.40 should have 
  intermediate solvation character

### What MU did not confirm

- Exact config fractions for the specific mixture
- SCE value for the specific mixture
- Whether the mixture hits SCE* = 1.466

### What the anchor data implies about SCE

Using the DME pure endpoint fractions as a lower bound estimate:

```
SSIP: 0.611    ln(0.611) = -0.492    contribution = 0.301
CIP:  0.354    ln(0.354) = -1.039    contribution = 0.368
AGG:  0.035    ln(0.035) = -3.352    contribution = 0.117

SCE(DME pure, 1.0M) ≈ 0.301 + 0.368 + 0.117 = 0.786
```

Note: This three-species SSIP/CIP/AGG model underestimates SCE because 
it collapses all AGG into one category and all SSIP into one category, 
ignoring the within-category config diversity. The true SCE from a full 
config population distribution (which distinguishes, e.g., 
(2DME,0FSI), (1DME,1FSI), (0DME,2FSI) etc.) will be higher than 0.786.

The SCE values in the framework dataset (e.g., DME at 1.0M = 1.2396) 
were computed from the full config population, not from the three-species 
SSIP/CIP/AGG summary. This is why MU's three-species anchor data cannot 
be directly compared to the framework's SCE values without a correction 
for within-category diversity.

### Critical note for interpretation

The SSIP/CIP/AGG fractions from Raman speciation are a **compressed 
representation** of the full config population. SCE computed from 
SSIP/CIP/AGG alone will always be lower than SCE computed from the 
full (n_solvent, n_anion) distribution. The relationship is:

```
SCE_full = SCE_SSIP/CIP/AGG + within-category entropy contribution
```

The within-category contribution is typically 0.3–0.6 SCE units for 
ether-based electrolytes based on the framework dataset. This correction 
must be applied before comparing MU Raman-proxy outputs to framework 
SCE values.

---

## Next Action

Use Pro Query 2 to request the interpolation estimate MU offered, with 
explicit weighting for 2-MeTHF donor number (~17) relative to DME (~20) 
and DEE (~15). This will give a labelled estimate of SSIP/CIP/AGG 
fractions for the exact 2-MeTHF:DME 0.60:0.40 mixture, from which 
a corrected SCE estimate can be computed.

---

## Query Log

| Query | Tool | Status | Queries remaining |
|-------|------|--------|-------------------|
| Candidate 1 solvation structure | Ask Pro | Complete | Pro: 9/10 |
| Candidate 1 interpolation estimate | Ask Pro | Pending | Pro: 9/10 |
| Candidate 2 solvation structure | Ask Pro | Not started | — |
| Arctic solvation structure | Ask Pro | Not started | — |
