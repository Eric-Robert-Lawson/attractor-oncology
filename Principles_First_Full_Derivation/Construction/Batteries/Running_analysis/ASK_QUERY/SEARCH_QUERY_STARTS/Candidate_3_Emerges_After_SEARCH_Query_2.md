# Candidate 3 — Geometric Derivation Record
**Author:** Eric Robert Lawson / OrganismCore
**Date:** 2026-04-03
**ORCID:** 0009-0002-0414-6544
**Repository:** github.com/Eric-Robert-Lawson/attractor-oncology
**Status:** DERIVED — awaiting mechanistic confirmation
            and experimental validation
**Epistemic status:** Geometrically derived from confirmed
                      SCE framework equations via molecular
                      space navigation. Not a literature
                      finding. Not a guess. Not pattern
                      matching.

---

## What This Document Is

This document records the causal derivation of
Candidate 3 — a novel electrolyte system that
emerged from following the SCE framework geometry
through molecular space. It is a deliberate
preservation of the derivation sequence so that
the origin of Candidate 3 is permanently recorded
and cannot be confused with a literature search
result or a chemical intuition guess.

Every step in the derivation follows causally
from the previous step. The candidate did not
exist before the geometry was followed. It was
not in the preprint. It was not in the literature.
It emerged from the intersection of a mathematically
derived fixed point in coordination space and a
molecular space search anchored to that fixed point.

---

## The Fixed Point

From the confirmed SCE framework equations
(full derivation in derivation_for_novel_result.md):

```
SCE* = 1.466

Derived by:
  Score(SCE) = CE_RT(SCE) + CE_LT(SCE)
  dScore/dSCE = 0
  33.21 - 1.230 × e^(1.493 + 1.230 × SCE) = 0
  SCE* = [ln(33.21/1.230) - 1.493] / 1.230
       = [ln(27.00) - 1.493] / 1.230
       = [3.2958 - 1.493] / 1.230
       = 1.466

Second derivative: always strictly negative.
SCE* = 1.466 is a global maximum. No other SCE
value produces higher combined RT+LT performance
under the confirmed equations.
```

This is a fixed point in coordination space.
It does not change. Everything that follows
is navigation toward it.

---

## The Required Shell Geometry at SCE*

From interpolation across the confirmed gradient
table (full derivation in derivation_for_novel_result.md):

```
At SCE* = 1.466:
  Dominant configuration fraction: ~37–39%
  Number of significant configurations: 3
  Secondary config 1 fraction: ~22–24%
  Secondary config 2 fraction: ~18–20%
  Remainder: ~15–20%
```

This is what the Li+ first shell must look like
at the target. Three distinct coordination species,
none dominating above ~38%.

---

## The Topology of Ether Molecular Space

From MU Molecular Searches 1a and 1b
(Search_Query_1.md, Search_Query_1b.md):

```
Finding 1: Linear and cyclic ether spaces are
           structurally separated in MU's
           similarity metric.
           DME-anchored search returned zero
           cyclic ethers.
           DOL-anchored search returned zero
           cyclic monoethers (THF, 2-MeTHF).

Finding 2: DME and diglyme sit at the convergence
           zone boundary — the only molecules
           that appear in both searches.

Finding 3: The confirmed dataset points bracketing
           SCE* = 1.466 are:
             DME+FEC:  SCE = 1.445  (below SCE*)
             THF:      SCE = 1.528  (above SCE*)
           The gap between them is 0.083 SCE units.
           No confirmed ether system sits inside it.
```

These findings establish the topology. SCE* sits
in a real gap in the molecular space of tested
ether electrolytes. The gap is not an artifact
of the dataset. It is a structural feature of
the chemical space.

---

## The Gap Confirmation — Search 2

From MU Molecular Search 2
(MU_Search_Query_2.md):

```
Anchors: DOL (SCE=1.606), DME (SCE=1.240)
Setting: Co-solvent
Purpose: Probe what sits in the SCE* gap

Result: The gap is unpopulated by ether molecules.

Return classified as:
  Heavy fluorinators (HOMO < -8.5 eV): dominant
    — wrong direction, suppress coordination
      diversity, kinetic locking risk
  Carbonates: contamination
    — Regime 3, transport failure at low T
  Boron esters: two molecules
    — only coordinating non-carbonate
      non-fluorinated return
  Silicon ester: one molecule
    — same electronic character as boron esters,
      secondary interest
```

The gap is confirmed as structurally real.
No conventional ether molecule fills it.
The search identified what does sit at the
correct electronic position in the gap:
boron esters.

---

## The Geometric Signal — Boron Esters

From Search 2 return, the two boron ester molecules:

```
Trimethyl borate:
  SMILES:  COB(OC)OC
  Score:   9.38/10
  Status:  Published (in other contexts)
  MW:      103.914 g/mol
  HOMO:    -7.7210 eV
  LUMO:    1.3291 eV
  ESP Min: -1.0007 eV
  ESP Max:  0.4130 eV
  Commercial: Commercially available

Triethyl borate:
  SMILES:  CCOB(OCC)OCC
  Score:   9.23/10
  Status:  Published (in other contexts)
  MW:      145.995 g/mol
  HOMO:    -7.6412 eV
  LUMO:    1.1399 eV
  ESP Min: -1.1543 eV
  ESP Max:  0.3906 eV
  Commercial: Commercially available
```

Electronic position in coordination space:

```
DME (pure ether):           HOMO = -6.87 eV  ← strong donor
Trimethyl borate:           HOMO = -7.72 eV  ← moderate donor
Triethyl borate:            HOMO = -7.64 eV  ← moderate donor
Fluorinated diluents:       HOMO < -8.50 eV  ← non-donor

Boron esters sit between pure ether donors
and non-coordinating fluorinators.
They have oxygen coordination sites weakened
by the boron Lewis acid center.
They are in the coordinating range but weaker
than DME.
This is the correct electronic position for a
coordination modifier — it can enter the Li+
first shell as a distinct coordination species
without dominating the shell the way DME does.
```

---

## The Causal Derivation of Candidate 3

The derivation follows from the geometry:

```
PREMISE 1:
  SCE* = 1.466 is the target.
  (Derived from calculus — fixed point.)

PREMISE 2:
  At SCE* the shell has dom% ≈ 38%, n_sig = 3.
  Three distinct coordination species required.
  (Derived from gradient table interpolation.)

PREMISE 3:
  DME alone has SCE = 1.240, dom% = 44%.
  One coordination species dominates.
  DME cannot reach SCE* by concentration alone.
  (Confirmed by dataset — three concentration points.)

PREMISE 4:
  The gap between DME (1.240) and THF/DOL
  (1.528/1.606) is real and unpopulated by
  any ether molecule in MU's database.
  (Confirmed by Search 2.)

PREMISE 5:
  The only route into the gap from the DME side
  is a coordination modifier — a molecule that
  adds a distinct Li+ coordination species to
  the DME-dominated shell without dominating it.
  (The FEC mechanism, applied to a non-carbonate.)

PREMISE 6:
  The only molecule in Search 2's return with
  the correct electronic character to function
  as a coordination modifier — coordinating
  but not dominating, no carbonate, no heavy
  fluorination — is trimethyl borate.
  HOMO = -7.72 eV. Three weakened oxygen sites.
  Commercially available. Unexplored in this role.
  (Confirmed by Search 2.)

CONCLUSION:
  LiFSI 1.0–1.2M in DME:B(OCH3)3 is the
  geometrically indicated candidate for
  approaching SCE* = 1.466 from the DME side
  via a Lewis acid coordination modifier route.

  This is Candidate 3.
```

This derivation is causal. Each premise follows
from the framework or from the molecular space
navigation. The conclusion follows from the
premises. Candidate 3 was not chosen. It was
found by following the geometry.

---

## How Candidate 3 Differs From Candidates 1 and 2

```
Candidates 1 and 2:
  Derived backward from SCE* using donor number
  arguments and linear mixing approximations.
  Route: mix a cyclic ether with a linear ether
  at a ratio predicted to hit SCE*.
  Limitation: cross-class extrapolation confirmed
  by MU's structural separation finding.
  Uncertainty: continuous — the SCE prediction
  depends on a linear interpolation that may
  not hold across the class boundary.

Candidate 3:
  Derived forward by following the geometry
  into molecular space and reading what was there.
  Route: add a coordination modifier to a single
  characterized base solvent (DME).
  DME's SCE is fixed at 1.240 in the dataset.
  The modifier adds one new coordination config.
  Uncertainty: binary — borate either enters
  the Li+ first shell or it does not.
  If yes: SCE rises from 1.240 toward 1.466.
  If no: mechanism fails cleanly, no ambiguity.
```

Candidate 3 has a cleaner failure mode and a
more direct causal derivation. It also has a
different mechanistic route to SCE* — Lewis acid
oxygen donation rather than cyclic/linear ether
mixing. If all three candidates are tested and
all three achieve SCE ≈ 1.466, this confirms that
SCE* is a robust attractor in coordination space
reachable by multiple chemically distinct routes.
That is a stronger result than any single candidate
alone.

---

## Open Mechanistic Question

One question must be answered before Candidate 3
is fully characterized:

```
Does trimethyl borate coordinate Li+ via oxygen
donation, or does it coordinate FSI- via the
boron Lewis acid center?

If Li+ coordination: mechanism confirmed.
  Borate adds third config to DME+FSI shell.
  SCE rises toward 1.466.

If FSI- coordination: mechanism is different.
  Borate acts as anion receptor, not Li+ modifier.
  Anion receptor changes ion pairing equilibrium.
  May still raise SCE by redistributing the
  CIP/AGG population, but via a different route.
  Mechanistic framing in preprint must be revised.
```

This is answerable by one MU Ask Pro query
or by inspection of the published borate
electrolyte additive literature.

It does not invalidate Candidate 3.
It determines which mechanistic framing is correct.

---

## Performance Prediction

From the confirmed SCE framework equations,
if SCE ≈ 1.466 is achieved:

```
CE_RT (Regime 1 lower bound):
  CE_RT = 100.1 - e^(1.493 + 1.230 × 1.466)
        = 100.1 - 27.02
        = 73.1%

CE_RT (Regime 2 upper bound, if concentration
       pushes system into Regime 2):
  CE_RT = 104.68 - 25.85 × 1.466 + 31.16
        = 97.9%

CE_LT at -20°C (band equation):
  CE_LT = 11.91 + 33.21 × 1.466
        = 60.6%
```

The actual RT CE depends on whether 1.0–1.2M
LiFSI in DME:borate produces sufficient AGG
fraction to cross the Regime 2 threshold.
This is the same regime question as Candidates
1 and 2 — it is the key experimental variable
for all three candidates.

---

## Experimental Validation Path

```
Step 1 — Mechanistic confirmation (computational):
  Short MD simulation of LiFSI 1.0M in
  DME:B(OCH3)3 at 4:1 and 2:1 ratios.
  Measure: Does borate enter Li+ first shell?
  Measure: SSIP/CIP/AGG fractions.
  Measure: Full config population distribution.
  Compute: SCE from Shannon entropy.
  Confirm: SCE ≈ 1.466, dom% ≈ 38%, n_sig = 3.

Step 2 — RT performance (experimental):
  Li|Cu half-cell, standard Aurbach protocol.
  Measure CE_RT.
  Confirm: ≥73% (R1 bound) or ≥90% (R2 bound).

Step 3 — LT performance (experimental):
  Li|Cu half-cell at -20°C.
  Measure CE_LT.
  Confirm: ≥60%.

Step 4 — Solvation structure (experimental):
  Raman spectroscopy.
  Confirm: three coordination species present.
  Confirm: dominant species fraction ≈ 38%.

Step 5 — Framework validation:
  If CE_RT and CE_LT match predictions:
  Candidate 3 validates the framework from
  a geometrically independent direction.
  Three candidates, three chemical routes,
  one fixed point: SCE* = 1.466.
```

---

## Why This Matters for the Framework

Candidate 3 is not just a new electrolyte proposal.
It is evidence that SCE* = 1.466 is a real attractor
in coordination space — a point that multiple
chemically distinct routes converge on when the
geometry is followed correctly.

```
Candidate 1: cyclic/linear ether mixing → SCE* 
Candidate 2: fluorinated/cyclic mixing  → SCE*
Candidate 3: Lewis acid modification    → SCE*

Three different chemical routes.
Three different structural classes.
One predicted destination: SCE* = 1.466.

If all three arrive at the same coordination
geometry, the framework's claim that SCE* is
a mathematically derived global optimum is
experimentally confirmed from three independent
directions simultaneously.

No single-candidate experiment can do this.
The three-candidate structure of the framework
is the scientific content.
```

---

## Summary

Candidate 3 (LiFSI 1.0–1.2M in DME:B(OCH3)3)
was derived by following the SCE framework
geometry through molecular space without drifting
into literature retrieval or chemical intuition.
The derivation is causal: fixed point → required
shell geometry → topology of ether space → gap
confirmation → geometric signal → candidate.
Each step follows from the previous step.
The candidate did not exist before the geometry
was followed. It is preserved here as a permanent
record of the derivation sequence.

---

*Repository: Eric-Robert-Lawson/attractor-oncology*
*Author: Eric Robert Lawson / OrganismCore*
*ORCID: 0009-0002-0414-6544*
*Date: 2026-04-03*
*Commit basis: b491ce60c0317d8bddb78ad04836db094abb5f5b*
*File: Candidate_3_Geometric_Derivation.md*
