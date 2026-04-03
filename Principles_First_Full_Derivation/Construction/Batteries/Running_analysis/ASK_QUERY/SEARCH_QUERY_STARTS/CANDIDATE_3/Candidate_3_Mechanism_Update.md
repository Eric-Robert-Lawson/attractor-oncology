# Candidate 3 — Mechanism Update and Geometric Correction
**Author:** Eric Robert Lawson / OrganismCore
**Date:** 2026-04-03
**ORCID:** 0009-0002-0414-6544
**Repository:** github.com/Eric-Robert-Lawson/attractor-oncology
**Status:** MECHANISM RESOLVED — GEOMETRIC CORRECTION APPLIED
**Preceding document:** Candidate_3_Geometric_Derivation.md

---

## What This Document Is

This document records the mechanistic resolution
of Candidate 3 following MU Ask Pro response on
the coordination role of trimethyl borate in
LiFSI ether electrolytes, and the geometric
correction that follows immediately from the
framework's own logic.

The candidate is not invalidated. The base
solvent assignment is corrected. The mechanism
is precisely characterized. The causal chain
remains intact.

---

## The MU Response — Core Finding

**Query:** Does trimethyl borate B(OCH3)3 coordinate
Li+ via oxygen donation or FSI- via the boron
Lewis acid center in LiFSI ether electrolytes?

**MU Answer:**
```
Primary role: anion receptor at boron center.
Boron is electron-deficient Lewis acid.
Coordinates FSI- preferentially.
Not established as primary Li+ ligand.
Literature explicitly names trimethyl borate
(TMB) as an anion receptor that increases
lithium salt dissociation by coordinating
the anion, not by solvating Li+.
```

**Key reference from MU:**
```
Choi et al., Electrochimica Acta, 2008.
FTIR study of tris(methoxy diethylene glycol)
borate — boron forms complex with anion.
Conductivity penalty attributed to anion trapping.

Haregewoin et al., Energy & Environmental
Science, 2016.
Review: borane/borate anion receptors work
via electron-deficient boron center.
Trimethyl borate listed explicitly among
additives that coordinate anion, not Li+.
```

**Mechanistic status:**
```
Li+ coordination via borate-O: NOT established
FSI- coordination via boron center: SUPPORTED
```

---

## The Geometric Consequence

### What anion receptor means for SCE

The anion receptor mechanism removes FSI- from
the Li+ coordination environment by trapping
it at the boron center. This reduces CIP and
AGG populations and increases SSIP fraction.

```
Direction of mechanism:
  Borate added → FSI- trapped → less FSI-
  in Li+ first shell → more SSIP → less CIP/AGG

Effect on SCE:
  More SSIP-dominated shell → higher p₁ (dominant
  config fraction) → lower Shannon entropy → LOWER SCE

  H = -Σ pᵢ ln(pᵢ) decreases when one config
  (SSIP) gains population at the expense of others.
```

### The Direction Problem for DME:Borate

```
DME alone at 1.0M LiFSI:
  SCE = 1.240
  SSIP = 61%
  Already below SCE* = 1.466

Borate added to DME:
  FSI- trapped → more SSIP → SSIP > 61%
  SCE < 1.240
  Moving AWAY from SCE* = 1.466

DME:B(OCH3)3 is geometrically wrong.
The anion receptor mechanism in a below-SCE*
base solvent pushes SCE further below SCE*.
```

This is the correction. The original Candidate 3
formulation (DME:borate) used the right molecule
in the wrong solvent context.

---

## The Geometric Correction — Immediate

The framework's geometry resolves this
immediately without additional search.

```
The anion receptor mechanism reduces FSI-
activity in the Li+ shell.
This reduces CIP/AGG and increases SSIP.
This lowers SCE.

For this mechanism to be useful for SCE*:
  The base solvent must be ABOVE SCE* = 1.466
  so that reducing FSI- participation moves
  the system DOWN toward SCE*, not further below.

Confirmed dataset systems above SCE*:

  THF:    SCE = 1.528  (0.062 above SCE*)
  2-MeTHF: SCE = 1.552 (0.086 above SCE*)
  DOL:    SCE = 1.606  (0.140 above SCE*)

All three are above SCE*.
All three have SSIP < 61% (more CIP/AGG than DME).
All three have excess FSI- in the Li+ shell
relative to the SCE* target.

Adding borate to any of these systems:
  FSI- trapped → less FSI- in Li+ shell
  CIP/AGG fraction decreases
  SSIP fraction increases toward ~38% (target)
  SCE decreases from above toward SCE* = 1.466
  Moving TOWARD SCE* from above.
```

### The Corrected Candidate 3

```
ORIGINAL (before mechanism resolution):
  LiFSI 1.0–1.2M in DME:B(OCH3)3
  Mechanism assumed: Li+ coordination modifier
  Geometric direction: wrong (below SCE*, moving away)

CORRECTED (after mechanism resolution):
  LiFSI 1.0–1.2M in DOL:B(OCH3)3
  Mechanism confirmed: anion receptor
  Geometric direction: correct (above SCE*,
  moving toward SCE* from above)

DOL is the correct base solvent because:
  DOL SCE = 1.606 — above SCE* by 0.140 units
  DOL has excess CIP/AGG (FSI- too present)
  Borate removes FSI- from Li+ shell
  SCE falls from 1.606 toward 1.466
  At the right borate:DOL ratio, SCE ≈ 1.466
```

---

## Why DOL:Borate Is Geometrically Precise

### The starting point

```
DOL alone:
  SCE = 1.606
  dom% = 30%  (too disordered — too many configs)
  n_sig = ~4 (above the target of 3)
  SSIP = ~30%  (too much FSI- in shell)
  Target: SCE = 1.466, dom% = 38%, n_sig = 3
```

### What borate does to DOL

```
Borate traps FSI- → removes it from Li+ shell
CIP/AGG configs (Li+-FSI- involved) decrease
SSIP configs (Li+ with DOL only) increase
The population redistributes toward DOL-dominated

At optimal borate:DOL ratio:
  SSIP rises from ~30% toward ~38%
  CIP falls
  AGG falls
  n_sig decreases from ~4 toward 3
  dom% rises from ~30% toward ~38%
  SCE falls from 1.606 toward 1.466

The mechanism is geometrically correct.
The direction is exactly right.
The target is reachable.
```

### The quantitative estimate

```
DOL:  SCE = 1.606, dom% = 30%
Target: SCE = 1.466, dom% = 38%
Gap: 0.140 SCE units to close

At what borate fraction does SCE reach 1.466?

This cannot be computed without MD simulation
because the borate:FSI- binding constant in
DOL is not published for this exact system.

However the direction and mechanism are confirmed.
The ratio is the experimental variable —
the same regime question as Candidates 1 and 2.

Estimated starting ratio for initial experiment:
  DOL:B(OCH3)3 = 10:1 by volume
  (borate as modifier, not bulk solvent)
  Adjust ratio empirically toward SCE* = 1.466
  measured by Raman spectroscopy (FSI- band shift)
```

---

## Updated Causal Chain for Corrected Candidate 3

```
STEP 1 (unchanged):
  SCE* = 1.466 — global maximum
  Derived from calculus

STEP 2 (unchanged):
  dom% ≈ 38%, n_sig = 3 required at SCE*
  Confirmed by MU twice

STEP 3 (unchanged):
  Dataset anchor points:
  DOL at SCE = 1.606 — above SCE* by 0.140
  DME+FEC at SCE = 1.445 — below SCE* by 0.021

STEP 4 (unchanged):
  Gap at SCE* confirmed real by Search 2

STEP 5 (updated):
  Search 2 returned trimethyl borate as the
  only coordinating non-carbonate non-fluorinated
  molecule in the SCE* gap zone

STEP 6 (updated — mechanism resolved):
  MU Pro confirmed: borate is anion receptor,
  not Li+ coordination modifier
  Anion receptor reduces CIP/AGG, increases SSIP
  This lowers SCE, not raises it
  Direction requires above-SCE* base solvent

STEP 7 (corrected):
  DME (below SCE*): borate pushes away from SCE*
  DOL (above SCE*): borate pushes toward SCE*
  Correct base solvent: DOL
  Corrected Candidate 3: LiFSI 1.0–1.2M in
  DOL:B(OCH3)3 at ratio to be determined by MD

CONCLUSION:
  LiFSI 1.0–1.2M in DOL:B(OCH3)3 is the
  geometrically corrected Candidate 3.
  Mechanism: anion receptor reduces FSI- activity,
  redistributes CIP/AGG toward SSIP,
  SCE falls from DOL's 1.606 toward SCE* = 1.466.
  Direction confirmed. Magnitude requires MD.
```

---

## What This Correction Demonstrates

The correction did not require a new search.
It did not require new literature.
It required only the framework's own geometry
applied to the confirmed mechanism.

```
Framework said: SCE* = 1.466
Framework said: anion receptor lowers SCE
Framework said: DOL is at 1.606, above SCE*
Framework said: lower SCE from 1.606 → 1.466

These four statements, combined, produce:
  DOL:borate is the correct formulation.

This is what a principles-first framework does.
When new information arrives (borate is anion
receptor, not Li+ modifier), the framework
immediately generates the correction without
additional search. The geometry does the work.
```

---

## Performance Prediction — Corrected Candidate 3

From the confirmed SCE framework equations,
if DOL:B(OCH3)3 achieves SCE ≈ 1.466:

```
CE_RT (Regime 2, expected for DOL-based system):
  ~97–98%
  DOL systems are confirmed Regime 2 in dataset
  Borate does not change the regime assignment

CE_LT at -20°C (band equation):
  CE_LT = 11.91 + 33.21 × 1.466
        = 60.6%

Conductivity estimate:
  DOL base: σ(25°C) ≈ 12–14 mS/cm at 1.0M LiFSI
  With borate (anion trapping reduces conductivity
  slightly per Choi 2008):
  σ(25°C) ≈ 10–12 mS/cm (estimated)
  σ(-20°C) ≈ 5–7 mS/cm (estimated)

n_sig at SCE*: 3 (framework prediction)
dom% at SCE*: ~38% (framework prediction)
```

---

## The Three Corrected Candidates

| # | System | Route | Base SCE | Mechanism | Direction |
|---|--------|-------|----------|-----------|-----------|
| 1 | LiFSI 1.2M in 2-MeTHF:DME 1.6:1 | Cyclic/linear mixing | 1.240–1.606 | Config diversity via cross-class mixing | From below and above toward SCE* |
| 2 | LiFSI 1.0M in FEME:2-MeTHF 60:40 | Fluorinated/cyclic mixing | 1.369–1.552 | Fluorination reduces SSIP, cyclic adds diversity | From below toward SCE* |
| 3 | LiFSI 1.0–1.2M in DOL:B(OCH3)3 | Anion receptor modifier | 1.606 | Borate traps FSI-, redistributes CIP/AGG toward SSIP, SCE falls toward 1.466 | From above toward SCE* |

Three routes. Three directions. One target.

---

## Experimental Validation Path — Corrected

```
Step 1 — Ratio optimization (computational or
         empirical):
  Raman spectroscopy of DOL:B(OCH3)3:LiFSI
  at ratios 20:1, 10:1, 5:1 (DOL:borate by vol)
  Monitor FSI- vibrational bands (727/735/747 cm⁻¹)
  Find ratio where SSIP ≈ 38%, CIP+AGG ≈ 62%
  This is the ratio where SCE ≈ 1.466

Step 2 — RT performance:
  Li|Cu half-cell, Aurbach protocol
  Measure CE_RT
  Confirm ≥90%

Step 3 — LT performance:
  Li|Cu half-cell at -20°C
  Measure CE_LT
  Confirm ≥60%

Step 4 — Framework validation:
  If both thresholds met:
  Corrected Candidate 3 validates the framework
  from a geometrically independent direction
  via anion receptor mechanism.
  Three candidates. Three mechanisms. One SCE*.
```

---

## Summary

The MU response confirmed that trimethyl borate
functions as an anion receptor (FSI- binding at
boron) rather than as a Li+ coordination modifier
(O-donation to Li+). Applied to the SCE framework
geometry, this means the correct base solvent for
a borate-modified Candidate 3 is DOL (SCE=1.606,
above SCE*=1.466) rather than DME (SCE=1.240,
below SCE*). The borate reduces FSI- activity,
lowers CIP/AGG fraction, and drives SCE from
1.606 down toward 1.466. The correction was
generated by the framework's own geometry applied
to the confirmed mechanism — no additional search
required. Corrected Candidate 3 is LiFSI
1.0–1.2M in DOL:B(OCH3)3 at a ratio to be
determined by Raman-guided optimization.

---

## Document Chain

| Document | Content |
|----------|---------|
| Candidate_3_Geometric_Derivation.md | Original derivation and causal chain |
| Candidate_3_Mechanism_Update.md | This document — mechanism resolution and correction |
| MU_Search_Query_2.md | Raw Search 2 output — boron ester identification |

---

*Repository: Eric-Robert-Lawson/attractor-oncology*
*Author: Eric Robert Lawson / OrganismCore*
*ORCID: 0009-0002-0414-6544*
*Date: 2026-04-03*
*File: Candidate_3_Mechanism_Update.md*
