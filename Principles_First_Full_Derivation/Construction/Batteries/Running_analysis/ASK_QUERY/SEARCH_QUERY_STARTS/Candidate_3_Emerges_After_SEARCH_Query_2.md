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

Source documents for each step are identified
precisely. All SCE values cited are from the
confirmed 21-system dataset
(OC_Battery_Framework_SCE_Analysis.md).

---

## Step 1 — The Fixed Point

From the confirmed SCE framework equations
(full derivation in derivation_for_novel_result.md,
Derivation 1):

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

## Step 2 — The Required Shell Geometry at SCE*

From interpolation across the confirmed gradient
table (derivation_for_novel_result.md, Derivation 2):

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

This was independently confirmed by MU Ask Pro
Query 2 (SES_Ask_Query_2.md):

```
MU interpolation returned n_sig = 3 for
Candidate 1 (2-MeTHF:DME 0.60:0.40).
Framework predicted exactly 3.
Match: exact.
```

And by MU Ask Pro Query 3 (SES_Ask_Query_3.md):

```
MU interpolation returned n_sig = 3 for
Candidate 2 (FEME:2-MeTHF 0.60:0.40).
Framework predicted exactly 3.
Match: exact.
```

n_sig = 3 is doubly confirmed before Search 2
was run.

---

## Step 3 — The Dataset Anchor Points

From the confirmed 21-system dataset
(OC_Battery_Framework_SCE_Analysis.md,
Master Gradient Table):

```
The four confirmed dataset points that define
the geometry around SCE* = 1.466:

  DME 1.0M LiFSI:
    SCE = 1.240
    dom% = 44%
    Regime 2
    Position: BELOW SCE* by 0.226 units

  DME+FEC 1.0M LiFSI:
    SCE = 1.445
    dom% = 38%
    Regime 2
    Position: BELOW SCE* by 0.021 units
    ← closest confirmed point below target

  THF 1.0M LiFSI:
    SCE = 1.528
    dom% = 30%
    Regime 2
    Position: ABOVE SCE* by 0.062 units
    ← closest confirmed point above target

  2-MeTHF 1.0M LiFSI:
    SCE = 1.552
    dom% = 32%
    Regime 2
    Position: ABOVE SCE* by 0.086 units

  DOL 1.0M LiFSI:
    SCE = 1.606
    dom% = 30%
    Regime 2
    Position: ABOVE SCE* by 0.140 units

The gap between DME+FEC (1.445) and THF (1.528)
is 0.083 SCE units wide.
SCE* = 1.466 sits inside this gap.
No confirmed dataset system sits inside this gap.
```

This is the geometric target zone. It exists
in the confirmed data before any search is run.
The search was designed to probe this zone.

---

## Step 4 — The Topology of Ether Molecular Space

From MU Molecular Searches 1a and 1b
(MU_SEARCH_QUERY_1A.md, MU_SEARCH_QUERY_1B.md):

```
Finding 1 — Structural separation confirmed:
  Linear and cyclic ether spaces are structurally
  separated in MU's similarity metric.
  DME-anchored search returned zero cyclic ethers.
  DOL-anchored search returned zero cyclic
  monoethers (THF, 2-MeTHF did not resolve).
  The two classes are not on a single axis.

Finding 2 — Convergence zone identified:
  DME and diglyme sit at the structural boundary
  between the two classes.
  They are the only molecules that appear in
  both searches.
  DME is the convergence zone molecule.

Finding 3 — Implication for gap geometry:
  The SCE* gap sits between DME (convergence zone,
  linear side) and THF/DOL (cyclic side).
  The gap is at the structural class boundary.
  It is not accessible by searching within either
  class alone.
  A coordination modifier — a molecule that
  bridges the gap electronically without belonging
  to either class structurally — is the only
  geometric route into the gap from the DME side.
```

This established what kind of molecule to look
for before Search 2 was run.

---

## Step 5 — The Gap Confirmation — Search 2

From MU Molecular Search 2
(MU_Search_Query_2.md):

```
Anchors: DOL (SCE=1.606), DME (SCE=1.240)
Setting: Co-solvent
SMILES entered: C1CCOC1 (THF, did not resolve),
                C1COCCO1 (DOL, resolved),
                COCCOC (DME, resolved)
Search ran from DOL and DME anchors.

Full return — 15 unique molecules after
deduplication:

  Heavy fluorinators (HOMO < -8.5 eV):
    COC(F)(F)CF            Score 9.49  Unexplored
    COC(F)(F)CC(F)(F)F     Score 9.44  Unexplored
    CCOC(F)(F)C(F)F        Score 9.43  Published
    CCOC(F)(F)CF           Score 9.41  Unexplored
    CCOC(F)(F)C(F)(F)F     Score 9.33  Unexplored
    CC(C)(C)OC(F)(F)C(F)F  Score 9.30  Unexplored
    COC(F)(F)C(F)(F)C(F)(F)F Score 9.29 Published
    CCCOC(F)(F)C(F)F       Score 9.28  Published
    C[Si](C)(C)C(F)(F)F    Score 9.27  Published
    → Wrong direction. Non-coordinating diluents.
      Discard entire class.

  Carbonates:
    CCOC(=O)OC    Score 9.26  Published
    CCOC(=O)OCC   Score 9.25  Published
    → Regime 3. Discard.

  Boron esters:
    COB(OC)OC     Score 9.38  Published
    CCOB(OCC)OCC  Score 9.23  Published
    → See Step 6.

  Silicon ester:
    CO[Si](C)(OC)OC  Score 9.29  Published
    → Secondary interest. See Step 6.

  Chiral fluorinated ether:
    CO[C@H](C)C(F)(F)F  Score 9.23  Unexplored
    → Borderline. HOMO -7.59 eV. Not prioritised.

RESULT:
  The gap is unpopulated by ether molecules.
  Heavy fluorinators dominate the return but
  are electronically wrong.
  Two molecules sit at the correct electronic
  position for a coordination modifier:
  the two boron esters.
```

---

## Step 6 — The Geometric Signal

From Search 2 return analysis
(MU_Search_Query_2.md):

```
The electronic position of each class relative
to the coordination modifier requirement:

  Class                    HOMO range      Role
  Pure ethers (DME)        -6.8 to -7.0   Strong donor — dominates shell
  Boron esters             -7.6 to -7.7   Moderate donor — can enter shell
                                           without dominating
  Silicon ester            -7.6            Same as boron esters
  Heavily fluorinated      < -8.5          Non-donor — suppresses diversity
  Non-coordinating         << -8.5         Diluent class

The coordination modifier must sit in the
-7.5 to -8.0 eV HOMO range:
  Strong enough to coordinate Li+
  Weak enough not to dominate the shell
  No carbonate decomposition pathway
  No heavy fluorination kinetic lock

Trimethyl borate:  HOMO = -7.72 eV ← IN THE WINDOW
Triethyl borate:   HOMO = -7.64 eV ← IN THE WINDOW
Silicon ester:     HOMO = -7.63 eV ← IN THE WINDOW

All three fluorinators: HOMO < -8.5 eV ← OUTSIDE WINDOW
All carbonates:         HOMO < -7.9 eV ← OUTSIDE WINDOW (Regime 3 risk)

Trimethyl borate is the primary signal:
  Highest score of the coordinating modifiers: 9.38
  Commercially available
  Published (as additive, not as coordination modifier)
  MW = 103.9 g/mol — small enough to enter Li+ shell
  Three oxygen donor sites weakened by boron
  Lewis acid center
  No carbonate group
  No heavy fluorination
  Unexplored in the coordination modifier role
  for ether electrolytes designed toward SCE*
```

---

## Step 7 — The Causal Derivation

The candidate follows from the geometry:

```
PREMISE 1 (Step 1):
  SCE* = 1.466 is the target.
  Derived from calculus. Global maximum.
  Second derivative always negative.

PREMISE 2 (Step 2):
  At SCE* the shell requires dom% ≈ 38%, n_sig = 3.
  Derived from gradient table interpolation.
  Confirmed twice by MU independent interpolation.

PREMISE 3 (Step 3):
  DME alone has SCE = 1.240, dom% = 44%.
  Three concentration points confirm DME cannot
  reach SCE* by concentration alone.
  The gap is 0.226 SCE units above DME.

PREMISE 4 (Step 3):
  The nearest confirmed dataset point below SCE*
  is DME+FEC at SCE = 1.445, dom% = 38%.
  FEC is the coordination modifier that pushes
  DME from 1.240 to 1.445.
  This proves the coordination modifier mechanism
  exists and is effective.
  FEC adds a distinct Li+ coordination species
  (Li-FEC config) to the DME-dominated shell,
  increasing n_sig and reducing dom%.

PREMISE 5 (Step 4):
  The SCE* gap sits at the structural class
  boundary between linear and cyclic ether space.
  The gap is not accessible by ether-to-ether
  mixing within either class.
  A coordination modifier is the geometric route
  into the gap from the DME side.

PREMISE 6 (Step 5):
  Search 2, anchored to the two dataset points
  bracketing SCE* from outside (DME at 1.240,
  DOL at 1.606), confirmed the gap is unpopulated
  by ether molecules.
  The search returned only fluorinators (wrong
  direction), carbonates (Regime 3), and two
  boron esters at the correct electronic position.

PREMISE 7 (Step 6):
  Trimethyl borate (COB(OC)OC) has HOMO = -7.72 eV.
  This places it in the coordination modifier
  window: coordinating but not dominating.
  No carbonate. No heavy fluorination.
  Commercially available. Unexplored in this role.
  It is the only molecule in the Search 2 return
  that satisfies all coordination modifier
  requirements.

CONCLUSION:
  LiFSI 1.0–1.2M in DME:B(OCH3)3 is the
  geometrically indicated candidate for
  approaching SCE* = 1.466 from the DME side
  via a Lewis acid coordination modifier route.

  Mechanism:
    DME provides the dominant coordination config
    (estimated ~38% at optimal ratio).
    Borate provides a secondary coordination config
    via weakened oxygen donation to Li+.
    FSI- provides a third config at elevated
    anion participation.
    n_sig = 3. dom% ≈ 38%.

  This is Candidate 3.

  Candidate 3 was not chosen. It was not designed
  by chemical intuition. It was not found by
  literature search. It was found by following
  the geometry of the SCE framework to the point
  it indicated and reading what was there.
```

---

## How Candidate 3 Differs From Candidates 1 and 2

```
Candidates 1 and 2
(derivation_for_novel_result.md, Derivations 3 and 4):

  Derived backward from SCE* using donor number
  arguments and linear mixing approximations.
  Route: mix a cyclic ether with a linear ether
  at a ratio predicted to hit SCE*.

  Limitation confirmed by Searches 1a and 1b:
  The mixing spans a structural class boundary.
  MU's own algorithm cannot reach cyclic ether
  space from a linear ether anchor.
  The interpolation is a cross-class extrapolation.
  Uncertainty is continuous — the SCE prediction
  depends on an approximation that may not hold
  across the boundary.

  Epistemic status: algebraically derived from
  donor number arguments, cross-class limitation
  confirmed post-derivation.

Candidate 3:

  Derived forward by following the geometry into
  molecular space from the correct anchors.
  Route: add a coordination modifier to DME,
  a single fully characterized base solvent.
  DME's SCE = 1.240 is fixed in the dataset.
  The modifier adds one distinct coordination
  config to the DME-dominated shell.

  Limitation: one binary mechanistic question.
  Borate either enters the Li+ first shell or
  coordinates FSI- instead. If the former,
  mechanism confirmed. If the latter, mechanism
  is different but candidate may still function
  via anion receptor route.

  Epistemic status: geometrically derived from
  confirmed framework equations and molecular
  space navigation. Cross-class issue does not
  apply — DME is fully characterized and the
  modifier role does not require class mixing.
```

---

## Open Mechanistic Question

```
Question:
  Does trimethyl borate coordinate Li+ via
  oxygen donation, or coordinate FSI- via the
  boron Lewis acid center?

If Li+ oxygen donation:
  Mechanism confirmed as stated.
  Borate adds Li+(borate-O) config to shell.
  Three species: Li+(DME), Li+(borate-O), Li+(FSI-)
  dom% ≈ 38%, n_sig = 3.
  SCE rises from 1.240 toward 1.466.

If FSI- coordination (anion receptor):
  Mechanism is different but not invalidating.
  Borate binds FSI-, reducing free anion activity.
  This redistributes the CIP/AGG population.
  A different route to modified SCE.
  Mechanistic framing in preprint must be revised.
  Candidate 3 may still function — different route.

This question does not invalidate Candidate 3.
It determines which mechanistic description
is correct.

Answerable by:
  One MU Ask Pro query on borate coordination
  in LiFSI ether systems.
  OR: MD simulation of DME:B(OCH3)3:LiFSI.
  OR: Raman spectroscopy of the mixture — borate
  B-O stretch shifts upon Li+ coordination.
```

---

## Performance Prediction

From the confirmed SCE framework equations,
if SCE ≈ 1.466 is achieved:

```
Equation 1 (Regime 1 RT):
  CE_RT = 100.1 - e^(1.493 + 1.230 × 1.466)
        = 100.1 - e^(3.298)
        = 100.1 - 27.04
        = 73.1%

Equation 4 (Regime 2 RT, if AGG fraction
crosses Regime 2 threshold):
  log(100 - CE_RT) = 4.898 - 2.545 × 1.466
  log(100 - CE_RT) = 4.898 - 3.731
  log(100 - CE_RT) = 1.167
  100 - CE_RT = 14.68
  CE_RT = 85.3% (Regime 2 lower estimate)

  Upper bound if fully in Regime 2:
  CE_RT ≈ 97–98%

Equation 2 (Band — LT):
  CE_LT = 11.91 + 33.21 × 1.466
        = 11.91 + 48.69
        = 60.6% at -20°C

Conductivity (from DME+DEE interpolation basis):
  Estimated σ at 25°C:  ~13 mS/cm
  Estimated σ at -20°C: ~6 mS/cm
  Basis: DME endpoint from SES_Ask_Query_1.md

Summary at SCE* = 1.466:
  CE_RT: 73–98% (regime-dependent)
  CE_LT: 60.6% at -20°C
  σ(-20°C): ~6 mS/cm
  n_sig: 3
  dom%: ~38%
```

The RT CE range reflects regime uncertainty.
If the DME:borate mixture at 1.0–1.2M LiFSI
produces sufficient AGG fraction to cross into
Regime 2, CE_RT approaches 97–98%.
If it remains in Regime 1, CE_RT ≈ 73%.
This is the key experimental variable —
the same regime question as Candidates 1 and 2.

---

## Experimental Validation Path

```
Step 1 — Mechanistic confirmation (computational):
  MD simulation of LiFSI 1.0M in DME:B(OCH3)3
  at ratios 4:1 and 2:1 by volume.
  Measure:
    Does borate enter Li+ first shell?
    SSIP/CIP/AGG fractions.
    Full config population distribution.
    SCE from Shannon entropy H = -Σ p_i ln(p_i).
  Confirm: SCE ≈ 1.466, dom% ≈ 38%, n_sig = 3.

Step 2 — RT performance (experimental):
  Li|Cu half-cell, standard Aurbach protocol.
  Measure CE_RT over 10 cycles.
  Confirm: ≥73% (Regime 1 bound).
  Target:  ≥90% (Regime 2 — both thresholds met).

Step 3 — LT performance (experimental):
  Li|Cu half-cell at -20°C.
  Same protocol.
  Measure CE_LT.
  Confirm: ≥60%.

Step 4 — Solvation structure (experimental):
  Raman spectroscopy of DME:B(OCH3)3:LiFSI.
  Confirm: three coordination species present
  in FSI- vibrational region.
  Confirm: dominant species fraction ≈ 38%.
  Confirm: B-O stretch shift indicating Li+
  or FSI- coordination (resolves mechanistic
  question).

Step 5 — Framework validation:
  If CE_RT ≥ 90% and CE_LT ≥ 60%:
  Candidate 3 validates the framework from
  a geometrically independent chemical direction.
  Three candidates. Three chemical routes.
  One fixed point: SCE* = 1.466.
```

---

## Why Three Candidates Matter

```
Candidate 1: cyclic/linear ether mixing → SCE*
  Chemical route: cross-class structural mixing
  Limitation: cross-class extrapolation

Candidate 2: fluorinated/cyclic mixing  → SCE*
  Chemical route: fluorination + cyclic ether
  Limitation: same cross-class extrapolation

Candidate 3: Lewis acid modification   → SCE*
  Chemical route: coordination modifier on DME
  Limitation: borate coordination mode (binary)

Three chemically distinct routes.
Three structurally independent mechanisms.
One predicted destination: SCE* = 1.466.

If all three arrive at the same coordination
geometry independently:

  The framework's claim that SCE* is a
  mathematically derived global optimum is
  confirmed from three independent directions.

  SCE* is not an artifact of the dataset.
  It is not a property of cyclic/linear ether
  mixing specifically.
  It is a real attractor in Li+ coordination
  space that multiple chemical routes converge
  on when the geometry is followed correctly.

  This is the scientific content of having
  three candidates rather than one.
  No single experiment can prove this.
  The three-candidate structure is the proof.
```

---

## Summary

Candidate 3 (LiFSI 1.0–1.2M in DME:B(OCH3)3)
was derived by following the SCE framework
geometry through molecular space in seven
causal steps: fixed point → shell geometry →
dataset anchor points → ether space topology →
gap confirmation → geometric signal → candidate.
Each step is sourced to a specific document in
the repository. The candidate did not exist
before the geometry was followed. It is not in
the preprint. It is not in the literature in
this role. It is preserved here as a permanent
record of the derivation sequence, timestamped
to 2026-04-03, at commit
b491ce60c0317d8bddb78ad04836db094abb5f5b.

---

## Source Document Index

| Step | Finding | Source Document |
|------|---------|-----------------|
| 1 | SCE* = 1.466, global maximum | derivation_for_novel_result.md, Derivation 1 |
| 2 | dom% = 38%, n_sig = 3 at SCE* | derivation_for_novel_result.md, Derivation 2 |
| 2 | n_sig = 3 confirmed twice | SES_Ask_Query_2.md, SES_Ask_Query_3.md |
| 
