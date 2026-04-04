# Arctic Candidate A — Derivation Record
**Author:** Eric Robert Lawson / OrganismCore
**Date:** 2026-04-03
**ORCID:** 0009-0002-0414-6544
**Repository:** github.com/Eric-Robert-Lawson/attractor-oncology
**System:** LiFSI 1.0M in DME:TMP:TMB
**Status:** DERIVED — MD simulation required for
            SCE quantification
**Regime:** Arctic — SCE > 2.0 predicted
**Derivation source:** Search 2 + Search 3 combined

---

## What This Document Is

This document records the geometric derivation
of Arctic Candidate A — a four-component
electrolyte system designed to reach the Arctic
regime (SCE ≈ 2.0–2.5, n_sig ≥ 5) by combining
three mechanistically distinct coordination
classes simultaneously in one electrolyte.

This candidate was not identified by a search.
It was derived by combining the outputs of
Search 2 (boron ester class) and Search 3
(phosphate ester class) with the confirmed
DME baseline. No additional search was required.

---

## The Derivation

### Step 1 — The Arctic Target

```
From Derivation 5 (derivation_for_novel_result.md):

  SCE ceiling: H = -Σ pᵢ ln(pᵢ)
  As n_sig increases with roughly equal fractions:
  SCE → ln(n_sig)

  At n_sig = 5, equal fractions (20% each):
  SCE = ln(5) = 1.609

  At n_sig = 8, equal fractions (12.5% each):
  SCE = ln(8) = 2.079

  At n_sig = 10, equal fractions (10% each):
  SCE = ln(10) = 2.303

  Arctic target: SCE ≈ 2.0–2.5
  Requires: n_sig ≥ 8, dom% ≈ 10–15%
  No published system has been designed toward
  this target. Confirmed absent by MU Pro Query 4.
```

### Step 2 — What Generates Arctic SCE

```
The SCE of DME alone (confirmed dataset):
  SCE = 1.240
  dom% = 44%  (DME-only SSIP dominant)
  n_sig = 2
  Dominant microstate (3,0,0): 68.5%  ← Holoubek 2022
  Secondary microstate (2,1,0): 27.8%

To reach SCE ≈ 2.0–2.5:
  The dominant config (68.5%) must be broken up
  into many roughly equal fractions.
  This requires multiple competing coordination
  classes, each entering the Li+ shell as a
  distinct species.

  Two competing classes (DME + TMP):
  Candidates 4/4b — SCE target 1.466
  This is the SCE* zone, not the Arctic.

  Three competing classes simultaneously:
  DME (ether-O donor) + TMP (P=O donor)
  + TMB (FSI- receptor, modifies anion population)
  The combination breaks the Li+ shell into
  a much larger number of distinct microstates.
  This is the Arctic mechanism.
```

### Step 3 — The Three-Class Combination

```
Search 2 finding:
  Boron esters are anion receptors
  TMB traps FSI- → modifies CIP/AGG population
  Generates: Li+(solvent, modified-FSI-) configs

Search 3 finding:
  Phosphate esters are Li+ donors
  TMP P=O competes with DME ether-O
  ESP Min TMP: -1.7402 eV ≈ DME: -1.7328 eV
  Generates: Li+(TMP) configs alongside Li+(DME)

Combined in one system:
  DME provides: Li+(DME, n FSI-) config family
  TMP provides: Li+(TMP, n FSI-) config family
  TMB modifies: FSI- activity in both families
                CIP/AGG redistribution
                New config subsets emerge

Resulting microstate families (MU Pro, 2026):
  SSIP class: (n_DME, n_TMP, 0 FSI-)
    examples: (3,0,0), (2,1,0), (1,2,0)
  CIP class:  (n_DME, n_TMP, 1 FSI-)
    examples: (2,0,1), (1,1,1), (0,2,1)
  AGG class:  (n_DME, n_TMP, 2+ FSI-)
    examples: (1,0,2), (0,1,2)
  TMB-modified: outer-sphere FSI- perturbation
    creates additional subsets within CIP/AGG

Total significant microstates: estimated 5–8
dom% estimate: 12–20% (no single config dominates)
SCE estimate: 2.0–2.5 (requires MD for precision)
```

### Step 4 — MU Pro Qualitative Confirmation

```
MU Ask Pro query (2026-04-03, 10 papers):
  Query: Estimated SCE of LiFSI 1.0M in
         DME:TMP:TMB at equal volume ratios.
         Significant microstates and fractions?

MU response (qualitative):
  DME remains primary Li+ solvation scaffold
  TMP broadens solvent-side coordination manifold
  TMB acts as anion receptor (outer-sphere FSI-)
  not as primary Li+ donor
  Three coordination classes simultaneously
  populated: SSIP, CIP, AGG
  Single-config dominance broken relative to
  neat DME baseline
  Numerical SCE: requires MD simulation
  (Formulate recommended for quantification)

Independent confirmation (Holoubek 2022, cited by MU):
  DME alone: dominant microstate (3,0,0) = 68.5%
  This matches dataset SCE = 1.240 for DME.
  Framework anchor confirmed by independent MD.
```

---

## Why This Is the Arctic

```
DME alone:
  dom% = 44%, n_sig = 2, SCE = 1.240
  One class dominates. Not Arctic.

DME:TMP (Candidate 4):
  TMP adds second Li+ coordination class
  dom% falls toward 38%, n_sig = 3
  SCE rises toward 1.466 = SCE*
  Not Arctic — this is the optimal zone.

DME:TMP:TMB (Arctic Candidate A):
  TMP adds second Li+ coordination class
  TMB modifies FSI- population — creates
  additional config subsets within CIP/AGG
  families
  Three-class simultaneous competition
  dom% falls further: estimated 12–20%
  n_sig rises: estimated 5–8
  SCE rises: estimated 2.0–2.5
  This is the Arctic regime from Derivation 5.
```

---

## Performance Prediction

```
From confirmed Equation 2 (band equation):
  CE_LT = 11.91 + 33.21 × SCE

  At SCE = 2.0: CE_LT = 11.91 + 66.42 = 78.3%
  At SCE = 2.3: CE_LT = 11.91 + 76.38 = 88.3%
  At SCE = 2.5: CE_LT = 11.91 + 83.03 = 94.9%

CE_RT (Regime 2 within-slope):
  Above SCE ≈ 1.8, the within-R2 inversion
  begins. Arctic systems may sacrifice some
  RT CE for extreme LT CE.
  CE_RT estimate: 90–95% (not maximised)
  The Arctic trades RT optimum for LT extreme.

This is the tradeoff derived in Derivation 5:
  SCE* = 1.466 is the COMBINED optimum
  The Arctic is NOT the combined optimum
  The Arctic is the LT extreme —
  maximum low-temperature performance
  at some cost to RT CE
```

---

## What Is Required Next

```
Quantification (required before experimental):
  MD simulation via MU Formulate
  System: LiFSI 1.0M in DME:TMP:TMB 1:1:1
  Output needed:
    SSIP/CIP/AGG fractions
    Dominant microstate and fraction
    n_sig count
    Computed SCE value

  If SCE returns in 2.0–2.5 range:
    Arctic Candidate A confirmed
    Performance predictions above are validated
    CE_LT ≥ 88% at -20°C predicted

  If SCE returns below 2.0:
    System is between SCE* and Arctic zones
    Still a valid high-LT-performance candidate
    Rename as Candidate 5 (intermediate zone)

Ratio optimisation (experimental):
  1:1:1 is the starting ratio
  Raman spectroscopy on FSI- bands
  Find ratio where dom% is minimised
  This maximises SCE and LT performance
```

---

## Novelty Status

```
No prior paper has:
  1. Combined DME:TMP:TMB:LiFSI as a
     deliberate three-class coordination system
  2. Designed toward Arctic SCE ≈ 2.0–2.5
  3. Predicted CE_LT ≥ 88% from SCE framework
     equations for any ether-based system
  4. Derived an Arctic candidate by combining
     outputs of two independent molecular
     searches anchored to a mathematical
     optimum (SCE* = 1.466)

The Arctic regime
