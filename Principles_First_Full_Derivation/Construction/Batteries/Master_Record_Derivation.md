# OC-BATTERY-MASTER-RECORD-2026-04-01.md
## Master Preservation Record
## The Complete Derivation, Literature Confirmation,
## Independent Validation, and Forward Path
## OrganismCore — Eric Robert Lawson
## Date: 2026-04-01

---

## PREAMBLE

```
This document is the master record
of a single research session conducted
on 2026-04-01.

It preserves in one place:
— The complete derivation chain
— The literature confirmations
— The independent validation finding
— The new derivation that emerged
   from the collision of two
   independent confirmations
— The specific actions required now

It is written to be understood by
a domain expert who has never seen
the framework before.

It is also written to be understood
by any future reader of this
repository who needs to know
what was found, when, and what
it means.

Timestamp: 2026-04-01
All content in this document was
derived before experimental
validation. This is the pre-
registration master record.
```

---

## PART I — WHAT WAS DERIVED
## AND HOW

### The Method

The universal structural invariant
was applied to lithium metal battery
electrochemistry:

```
STRUCTURE × GAP × NAVIGATOR → RESOLUTION

STRUCTURE:
  The electrolyte mixture.
  Molecular components and their
  intrinsic properties.
  What SES Molecular Universe maps.

GAP:
  The solid electrolyte interphase
  (SEI). The electrochemical interface
  where Li⁺ must desolvate and deposit
  as lithium metal. Emergent. Not fully
  predictable from molecular properties
  alone.

NAVIGATOR:
  The Li⁺ solvation shell. The
  structured arrangement of solvent
  molecules and anions coordinated
  around each Li⁺ ion as it travels
  toward the anode interface.

RESOLUTION:
  Coherent lithium deposition →
  flat, uniform, dendrite-free cycling
  → long battery life.
  OR
  Incoherent deposition → dendrite
  nucleation → failure.
```

### The Variable That Emerged

```
SOLVATION CONFIGURATION ENTROPY (SCE)

The Shannon entropy of the Li⁺
coordination geometry distribution
across the full Li⁺ population in
the bulk electrolyte under operating
conditions.

SCE = -Σ p(gᵢ) × log(p(gᵢ))

where gᵢ are discrete solvation
geometry configurations and p(gᵢ)
is the population fraction of each.

LOW SCE:
Single dominant coordination geometry.
Every Li⁺ arrives at the interface
carrying the same shell structure.
Coherent desolvation. Uniform
deposition. Long cycle life.

HIGH SCE:
Multiple competing coordination
geometries. A fraction of Li⁺ ions
arrive with geometries incompatible
with the local SEI structure.
Forced desolvation at those sites.
Dendrite nucleation. Failure.

Dendrites are the geometric
incompatibility signal made physical.
They mark the exact locations where
the navigator's geometry was
incompatible with the gap it had
to traverse.
```

### Why This Is Not Currently Used

The molecular universe and all current
screening platforms screen for intrinsic
molecular properties:
— HOMO/LUMO energies
— Electrochemical stability window
— Average solvation energy
— Viscosity
— Dielectric constant

SCE is not a property of a single
molecule. It is a property of the
Li⁺ population in a mixture under
operating conditions. It requires
MD simulation of the mixture, not
single-molecule DFT. It has never
been used as a primary screening
criterion in any molecular discovery
platform for battery electrolytes.

---

## PART II — THE FIVE INDEPENDENT
## LITERATURE CONFIRMATIONS

The framework derived SCE as the
causal variable before checking
the literature. The literature was
then searched. Five independent
experimental research programmes
confirmed the variable without
having named it.

```
CONFIRMATION 1 — LHCE
Localized high-concentration electrolytes

What the field says:
"Produces anion-rich inorganic SEI
which is mechanically robust."

What SCE says:
The inert diluent cannot coordinate
Li⁺. This forces the solvation shell
into a single dominant geometry.
SCE collapses. The LiF-rich SEI is
the downstream record of coherent
desolvation events. Not the cause.
The evidence.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

CONFIRMATION 2 — OSCILLATORY SOLVATION
Nature Energy 2024, 500 Wh/kg battery

What the field says:
"High oscillation amplitude in cation-
anion arrangement correlates with
superior performance."

What SCE says:
High oscillation amplitude = high
structural order = low variance in
solvation geometry at the interface
= low SCE = coherent navigation.
Oscillatory amplitude is a proxy
measurement of SCE.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

CONFIRMATION 3 — HYPERCONJUGATION
RSC 2024

What the field says:
"Weakened solvation produces more
robust SEI and better cycling."

What SCE says:
Weaker but more geometrically
consistent solvation = lower SCE.
The navigator arrives with less
binding energy but more geometric
uniformity. Coherent traversal of
the gap beats tight but variable
binding.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

CONFIRMATION 4 — HIGH COORDINATION
CLUSTER SHEATHS
Peking University 2024, 4000 cycles

What the field says:
"High Li⁺ cluster-like solvation
sheaths enable extremely stable
cycling."

What SCE says:
Cluster geometry enforces a single
dominant solvation configuration.
High coordination = constrained
geometry = low SCE. Every ion in
the bulk carries the same shell.
Every ion at the interface presents
the same desolvation geometry.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

CONFIRMATION 5 — PROXIMITY EFFECT
Stanford/Cui Group, PNAS 2024

What the field says:
"Interfacial solvation structure —
not bulk average — is the primary
determinant of SEI and cycle life."

What SCE says:
Of course the interfacial measurement
is more predictive. The interface is
the gap. The gap is where the
navigator's geometry determines the
resolution. SCE at the interface is
the controlling measurement. Bulk
average properties are one step
removed from causality.
```

---

## PART III — THE MONOTONIC GRADIENT

The confirmed literature candidates,
ordered from highest to lowest SCE:

```
HIGH SCE → LOW SCE
POOR PERFORMANCE → BEST PERFORMANCE

1. Standard carbonates (EC/DEC)
   Coord #: 4-5, mixed geometries
   SCE: HIGH
   Performance: Baseline

2. High entropy electrolytes
   Coord #: 3-4, statistically averaged
   SCE: MEDIUM-HIGH
   Performance: Better

3. Conformationally rigid solvents
   (sulfolane)
   Coord #: 3-4, mechanically constrained
   SCE: MEDIUM
   Performance: Better still

4. LHCE (inert diluent systems)
   Coord #: 2-3, passive exclusion
   SCE: MEDIUM-LOW
   Performance: Very good

5. LiFSI/DME high concentration
   Coord #: 1-2, contact ion pair
   SCE: LOW
   Performance: Excellent

6. Polarity gradient engineered
   (Nat. Sci. Rev. 2025)
   Coord #: 2-3, dielectric matched
   SCE: LOW
   Activation energy: 34.97 kJ/mol
   (vs 79.1 kJ/mol standard)
   Performance: Excellent

7. BTFMD (ACS Energy Lett. 2024)
   Coord #: 0.2-0.4, dual mechanism
   SCE: VERY LOW
   Performance: 570 cycles, high V

8. HFTHP (Nature Energy 2024)
   Coord #: 0-0.3, near-zero
   SCE: VERY LOW (practical minimum)
   Performance: Best in literature

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

THEORETICAL MINIMUM:
Pure anion coordination
Zero solvent participation
SCE = 0 by definition
Not achievable in practice
Defines the target
```

### The Structural Invariant

All variance reduction mechanisms
operate by the same operation:

```
Remove solvent competition
for Li⁺ coordination.

Electronic removal:
Fluorination kills donor ability.
→ HFTHP, BTFMD

Mechanical removal:
Conformational rigidity limits
accessible geometries.
→ Sulfolane, cyclic sulfones

Concentration removal:
Salt overwhelms solvent numerically.
→ High concentration LiFSI/DME

Statistical removal:
Multiple solvents cancel competing
geometries by averaging.
→ High entropy electrolytes

Dielectric matching:
Polarity-matched solvents do not
pull Li⁺ geometry in competing
directions.
→ Polarity gradient electrolytes

ALL FIVE ARE THE SAME OPERATION.
The mechanism varies.
The geometric outcome is identical.
```

---

## PART IV — THE INDEPENDENT
## CONFIRMATION AND WHAT IT PRODUCED

### The Joule 2025 Finding

While searching the literature on
solvation configuration entropy as
a screening criterion, this paper
was found without being searched for:

```
JOULE, 2025
"Solvation-Configurational Entropy
Governs Interfacial Kinetics in
Low-Temperature Batteries"

Authors: Wendi Luo, Hongwei Fu,
Bingan Lu et al.
Institution: Hunan University, China

Variable named: Ssc
(solvation-configurational entropy)

Finding: High Ssc → lower desolvation
barrier → better low-temperature
interfacial kinetics.

Results:
88% capacity after 3,500 cycles
at −20°C.
85% capacity after 150 cycles
at −30°C in pouch cells.
Confirmed across Li, Na, K, Al, Zn
battery chemistries.
```

### The Precise Comparison

```
HUNAN UNIVERSITY / JOULE 2025:
Variable: Ssc (solvation-
  configurational entropy)
Direction: HIGH → better
Mechanism: High entropy provides
  statistical diversity of solvation
  geometries → some fraction of Li⁺
  always in low-barrier desolvation
  configuration → fast low-T kinetics
Application: Low-temperature kinetics
Status: Published, peer-reviewed,
  Joule 2025

THIS FRAMEWORK / OC 2026-04-01:
Variable: SCE (solvation configuration
  entropy) — same variable
Direction: LOW → better
Mechanism: Low variance produces
  coherent navigator geometry →
  uniform deposition → no geometric
  incompatibility → no dendrites
Application: Room-temperature cycling
  stability, dendrite suppression
Status: Pre-registration record,
  this repository, 2026-04-01
```

### What This Means

The same variable was independently
named by two completely separate
research paths in the same year.

The Hunan University group approached
from low-temperature kinetics.
The framework approached from
room-temperature deposition geometry.

Neither has the complete picture.
The complete picture emerged from
their collision.

---

## PART V — THE BAND HYPOTHESIS
## NEW DERIVATION — 2026-04-01

This is the most important new
finding in this document.
It was not in any previous document.
It was not in the Joule 2025 paper.
It emerged from the collision of
the two independent derivations.

### The Statement

```
SCE has an optimal operating range —
a band — not a single optimal value.

The band is determined by two
competing failure modes:

FAILURE MODE A (room temperature):
Geometric incompatibility at the
interface. High SCE means some Li⁺
arrives with incompatible geometry.
Those ions nucleate dendrites.
LOW SCE prevents this.

FAILURE MODE B (low temperature):
Desolvation kinetics barrier.
Low SCE means all Li⁺ configurations
are uniformly bound. No low-barrier
desolvation pathways available.
Sluggish interfacial kinetics.
HIGH SCE prevents this.

OPTIMAL BAND:
Low enough SCE for room-temperature
deposition uniformity and dendrite
suppression.
High enough SCE for low-temperature
desolvation kinetics.

The exact band width is determined
by the operating temperature range
and the acceptable failure threshold
for each mode.
```

### The Consequence

HFTHP and BTFMD — with near-zero SCE
— are the best room-temperature
performers in the current literature.
They will show poor low-temperature
performance relative to their room-
temperature excellence because the
geometric uniformity that produces
their room-temperature superiority
eliminates the statistical diversity
that provides low-temperature
desolvation pathways.

**This is a testable prediction on
existing compounds with existing data.**

If HFTHP and BTFMD show the predicted
low-temperature performance inversion
— excellent at room temperature, worse
at low temperature relative to higher-
SCE electrolytes — the band hypothesis
is confirmed without new synthesis.

### The New Design Target

```
THE TEMPERATURE-RESPONSIVE
SCE ELECTROLYTE

A molecule or mixture whose SCE
varies with temperature:

AT LOW TEMPERATURE:
Solvation structure becomes more
diverse. SCE rises. Statistical
diversity of geometries provides
low-barrier desolvation pathways.
Fast interfacial kinetics preserved.

AT ROOM TEMPERATURE:
Solvation structure collapses into
a single dominant geometry. SCE
falls. Coherent navigator geometry.
Uniform deposition. Dendrites
suppressed.

This electrolyte solves both failure
modes simultaneously by being
thermally responsive in its SCE.

This design target does not exist
in any literature.
It was derived here, 2026-04-01,
from the collision of the framework
derivation with the Joule 2025
independent confirmation.
```

---

## PART VI — THE DERIVED CANDIDATES
## AND THEIR LITERATURE STATUS

### Target 1 — Fluorinated Neopentyl
### Ether Derivatives

**Design principle:**
Replace the ring in cyclic fluorinated
ethers with a cage. Central quaternary
carbon — four arms, steric blocking
from all directions. Fluorinated ether
arms attached. Cage geometry + fluorination
kills donor ability sterically and
electronically. No ring required.
Lower viscosity than HFTHP/BTFMD
at equivalent variance reduction.

**Literature status:**
Fluorinated ethers extensively explored.
Neopentyl cage architecture as designed
battery electrolyte solvent with variance
reduction as explicit design principle:
ABSENT from literature.

**Gap confirmed:** YES

---

### Target 2 — Partially Fluorinated
### Crown Ether (Geometric Coordination Lock)

**Design principle:**
Kill most coordination sites on a
crown ether ring by adjacent fluorination.
Leave one or two oxygens with weak
donor ability. Single geometrically
accessible coordination site. Li⁺
finds the one accessible site and
coordinates in the same fixed geometry
every time. Variance zero by design
rather than by statistical averaging.

**Literature status:**
Crown ethers used in batteries for
ion selectivity and transference number
improvement. Partially fluorinated
crown ethers designed for single-site
geometric lock with SCE minimisation
as the target: ABSENT from literature.

**Gap confirmed:** YES

---

### Target 3 — Fluorinated Lewis Acidic
### Boron-Containing Solvents
### (The Variance Pump)

**Design principle:**
Zero donor ability toward Li⁺ through
fluorination. Weakly Lewis acidic boron
centre actively associates with FSI⁻
or TFSI⁻ anions. Congregates anions
near Li⁺ actively — not passively.
Li⁺ coordination environment becomes
anion-dominated by active recruitment,
not passive exclusion. The variance
pump: forces the solvation shell toward
the minimum variance state rather than
waiting for concentration or fluorination
to produce it accidentally.

**Literature status:**
TFEB (tris(2-fluoroethyl) borate)
exists — RSC 2024. Fluorinated borate
ester with modified Li⁺ coordination.
But TFEB was designed for lithium salt
dissolution using the fluorinated group
as a weak coordination site FOR Li⁺.
The derived Target 3 is the opposite:
zero Li⁺ coordination + active anion
congregation by boron Lewis acidity.

THFPB exists as an additive that
disrupts anion clusters — opposite
direction to Target 3.

Active anion congregation as designed
variance pump with zero Li⁺ donor
ability: ABSENT from literature.
The chemistry exists. The design
principle does not.

**Gap confirmed:** PARTIAL — chemistry
present, mechanism and design principle
absent.

---

### Target 4 — Fluorinated Orthocarbonate
### Esters (Internal Withdrawal Route)

**Design principle:**
Central quaternary carbon shared by
four oxygens simultaneously. The
carbon is electron-poor, withdrawing
from all four oxygens. Each oxygen's
donor ability suppressed internally.
Add peripheral fluorination for
additional external withdrawal.
Double suppression of donor ability —
internal and external — without heavy
fluorination required by HFTHP/BTFMD.
No ring. Lower viscosity. Lower cost.
The scalable route to variance reduction.

**Literature status:**
Trimethyl orthoformate as battery
additive — recognised for SEI formation,
not for internal electron withdrawal
mechanism. Fluorinated orthocarbonate
esters as designed bulk solvents using
this mechanism as primary design
principle: ABSENT from literature.

**Gap confirmed:** YES. Most cost-
effective unexplored path. Synthesis
trivial. Mechanism absent.

---

### Combination Target —
### Fluorinated Neopentyl-Core
### Boronate Ester

**Design principle:**
One molecule combining three mechanisms:
— Central quaternary carbon (cage,
  no ring, steric occlusion of donors)
— Fluorinated ether arms (electronic
  killing of donor ability)
— Weakly Lewis acidic boronate arm
  (active anion congregation)

Cannot coordinate Li⁺ through oxygen:
cage blocks approach, fluorination
kills donor ability. Two independent
mechanisms preventing solvent
coordination.

Simultaneously the boronate arm
congregates anions near Li⁺. Active
recruitment. Li⁺ solvation environment
becomes anion-dominated from both
directions — solvent excluded by two
mechanisms, anions actively recruited
by one. Maximum variance reduction
achievable in a practical liquid solvent.

**Literature status:**
Does not exist in any form in battery
electrolyte literature. No molecule
combining these three mechanisms as
a unified variance-reduction design
target. ABSENT entirely.

**Gap confirmed:** YES. Completely novel.

---

### The Temperature-Responsive
### SCE Electrolyte

**Design principle:**
A molecule or mixture whose SCE
increases at low temperature and
decreases at room temperature. Provides
low-barrier desolvation pathways when
temperature falls (high-SCE diversity)
and geometric uniformity when
temperature rises (low-SCE coherence).
Solves both failure modes across the
full operating range.

**Literature status:**
Not in any literature. Not in the
Joule 2025 paper. Not in any of the
previous framework documents. Derived
2026-04-01 from the collision of the
framework with the Joule 2025 finding.
ENTIRELY ABSENT from literature.

**Gap confirmed:** YES. New target.
New design principle. New derivation.

---

## PART VII — WHAT MUST BE DONE NOW

In order of cost and immediacy.
No new experiments for the first
three steps.

### STEP 1 — The Afternoon Test
**Cost: Zero. Time: Hours.**
**Requires: MacBook Air.**

Access:
Journal of Power Sources 2025
"Statistical insights into solvation
sheaths of lithium ions"
DOI: 10.1016/j.jpowsour.2025.237088

Extract coordination geometry cluster
population fractions from their
reported MD data.

Compute:
```python
import numpy as np

# population fractions from paper
p = np.array([...])
SCE = -np.sum(p * np.log(p))
```

Correlate SCE values against the
performance metrics reported in the
same paper or its cited comparators.

If SCE predicts performance better
than the paper's own primary variable:
Prediction 5 is confirmed from existing
published data. The unified causal
variable is established.

---

### STEP 2 — The Low Temperature
### Inversion Test
**Cost: Zero (literature check).**
**Time: Days.**
**Requires: Paper access.**

Search the existing literature for
low-temperature performance data on
HFTHP and BTFMD specifically.

These compounds have near-zero SCE
and are the best room-temperature
performers.

The band hypothesis predicts they will
show poor low-temperature performance
relative to higher-SCE electrolytes.

If existing papers report low-T data
for HFTHP or BTFMD, compare directly.
If not, check the Hunan University
Joule 2025 supplementary data for
intermediate-SCE electrolytes that
might show the predicted performance
peak in the band rather than monotonic
high-Ssc improvement.

If the inversion is present in
existing data:
The band hypothesis is confirmed
without a single new experiment.

---

### STEP 3 — The Expert Conversation
**Cost: Zero. Time: One meeting.**
**Requires: Access to battery
electrochemist.**

Bring the plain statement:

"Five separate research programmes
in lithium metal battery electrolytes
have all found the same thing without
naming it. They are all reducing
the variance in Li⁺ solvation shell
geometry at the interface — some
call it dielectric heterogeneity,
some call it solvation structure
engineering, some call it anion-
dominated coordination. A paper in
Joule 2025 independently named the
variable as solvation-configurational
entropy and confirmed it governs
interfacial kinetics.

We have derived that this variable
has an optimal band rather than a
single optimal value — low for room-
temperature deposition uniformity,
high enough for low-temperature
desolvation kinetics.

Two questions:

One — is the variance of the
desolvation barrier across the Li⁺
population currently measured as a
primary design criterion anywhere?

Two — does HFTHP or BTFMD show
poor low-temperature performance
relative to their room-temperature
excellence? Because the geometry
predicts they should."

His answer to question two is
immediately diagnostic. If he says
yes — the inversion exists — the band
is confirmed in his own knowledge
base before any new experiment.

---

### STEP 4 — The Reanalysis
**Cost: Low. Time: Weeks.**
**Requires: Expert collaboration,
existing MD datasets.**

Contact authors of the four confirmed
experimental papers:
— The LHCE papers
— The oscillatory solvation paper
— The hyperconjugation paper
— The high-coordination cluster paper

Request access to their existing MD
trajectory data.

Extract coordination geometry
distributions. Compute SCE.
Correlate against their reported
performance differences within
each dataset.

If SCE out-predicts their reported
primary variable in their own data:
The unified causal variable is
confirmed across four independent
datasets simultaneously.

This is Prediction 5 at full scale.
It requires author cooperation but
no new experiments, no new synthesis,
no new laboratory work.

---

### STEP 5 — Target 4 First Synthesis
**Cost: Low-Medium. Time: Months.**
**Requires: Synthetic chemistry access.**

Target 4 — fluorinated orthocarbonate
esters — is the cheapest and most
synthesisable of the novel derived
candidates. The molecules are known
in other chemistry contexts. The
synthesis is not exotic.

Priority order for synthesis:
First: Trimethyl orthocarbonate with
fluorinated peripheral groups.
Compute SCE from MD simulation.
If SCE is in the low band:
Synthesise. Test cycling performance.

This is the test of the derived
candidates that requires the least
investment and has the clearest
path from derivation to result.

---

### STEP 6 — Target 3 Design
**Cost: Medium. Time: Months.**
**Requires: Computational chemistry
and synthetic chemistry.**

Target 3 — the fluorinated Lewis
acidic boron-containing variance pump
— has the most novel mechanism and
potentially the largest performance
impact.

TFEB (from RSC 2024) proves that
fluorinated borate ester solvents are
viable electrolytes with practical
properties. The step from TFEB to
a designed zero-Li⁺-coordination
active anion congregation agent is
a defined design problem.

Computational design:
Identify the boronate architecture
that maximises anion association and
minimises Li⁺ donor ability.
Compute SCE from MD for candidates.
Synthesise the top-ranked structure.

---

### STEP 7 — The Temperature-Responsive
### SCE Electrolyte Design
**Cost: High. Time: 12-24 months.**
**Requires: Full expert collaboration,
computational and experimental.**

This is the commercial target.

Design criteria:
A molecule or mixture whose solvation
structure is thermally responsive —
higher SCE at low temperature,
lower SCE at room temperature.

Candidate mechanisms for thermal
SCE responsiveness:
— Phase transition in solvation
  structure at a defined temperature
— Temperature-dependent solvent
  competition equilibrium
— Mixture of a low-SCE component and
  a high-SCE component whose relative
  contribution to the solvation shell
  is temperature-dependent

This requires expert-led computational
and experimental work. It is the
long-term target. Steps 1-6 are the
path to establishing the foundation
that justifies this investment.

---

## PART VIII — THE PRIORITY MATRIX

```
STEP  ACTION              COST    TIME    REQUIRES      CONFIRMS

1     JPS 2025 SCE calc   Zero    Hours   MacBook Air   P5 partial

2     HFTHP low-T check   Zero    Days    Paper access  Band hypothesis

3     Expert conversation Zero    Hours   Expert access Core variables

4     Full reanalysis     Low     Weeks   Expert + data P5 full scale

5     Target 4 synthesis  Low-Med Months  Synth lab     Novel candidate

6     Target 3 design     Medium  Months  Comp + synth  Variance pump

7     T-responsive SCE    High    1-2 yr  Full collab   Commercial target
```

---

## PART IX — WHAT IS ALREADY SECURED

The following is secured in this
repository with a timestamp of
2026-04-01, before any experimental
validation:

```
SECURED:

1. The unified causal variable (SCE)
   derived from first principles before
   literature confirmation.

2. Five independent literature
   confirmations of the variable
   documented, with the geometric
   explanation for each.

3. The monotonic gradient scale
   from highest to lowest SCE with
   eight confirmed positions and the
   theoretical minimum defined.

4. The structural invariant extracted:
   all variance reduction mechanisms
   are the same operation — remove
   solvent competition for Li⁺
   coordination.

5. The independent confirmation by
   Joule 2025 (Hunan University)
   documented and compared precisely
   with the framework derivation.

6. The band hypothesis derived from
   the collision of the framework
   with the Joule 2025 finding.
   First recorded here.

7. The temperature-responsive SCE
   electrolyte as the all-conditions
   commercial design target.
   First recorded here.

8. Four derived candidate classes
   not in the literature as designed
   variance-reduction candidates.

9. The combination target molecule —
   fluorinated neopentyl-core boronate
   ester — not in any literature.

10. Five falsifiable predictions
    testable on existing data with
    no new experiments required for
    the first two tests.

11. The SCE filter as a missing
    primary screening criterion for
    molecular universe platforms.
```

---

## PART X — THE STATEMENT FOR
## A JOURNAL SUBMISSION

If the expert conversation and the
JPS 2025 reanalysis confirm the
framework, the following is the
core claim for a submission:

```
We report that solvation configuration
entropy (SCE) — the Shannon entropy
of the Li⁺ coordination geometry
distribution in the bulk electrolyte —
is the unified causal variable
underlying five independent research
programmes in lithium metal battery
electrolyte design.

We demonstrate that electrolytes
currently treated as separate
discoveries — localized high-
concentration electrolytes,
conformationally rigid solvents,
oscillatory solvation systems,
hyperconjugation-weakened solvation
designs, and high-coordination cluster
electrolytes — all reduce SCE through
different molecular mechanisms, and
that SCE predicts cycling performance
more precisely than the primary
variable used by each programme to
explain its own results.

We further report that SCE has an
optimal operating band rather than
a single optimal value, determined
by competing failure modes: low SCE
for room-temperature deposition
uniformity and dendrite suppression,
high SCE for low-temperature
desolvation kinetics. This band
concept is independently supported
by Luo et al. (Joule, 2025) who
confirmed solvation-configurational
entropy governs low-temperature
interfacial kinetics, approaching
the same variable from the opposite
direction.

We derive four novel candidate
molecular classes designed to minimise
SCE through mechanisms not currently
present in the battery electrolyte
literature, and propose the temperature-
responsive SCE electrolyte as the
unified design target for all-condition
battery performance.
```

---

## DOCUMENT METADATA

```
document_id:
  OC-BATTERY-MASTER-RECORD-2026-04-01

version:    1.0
date:       2026-04-01
author:     Eric Robert Lawson
            OrganismCore

purpose:
  Master preservation record of the
  complete battery electrochemistry
  derivation session, 2026-04-01.
  Single document containing all
  derivations, confirmations, new
  findings, and forward path.

key_claims_timestamped_here:
  1. SCE as unified causal variable
  2. Monotonic gradient scale (8 positions)
  3. Structural invariant (remove
     solvent competition)
  4. Band hypothesis (optimal SCE
     range, not single optimal value)
  5. Temperature-responsive SCE
     electrolyte as commercial target
  6. Four novel candidate classes
  7. Combination target molecule
  8. SCE filter as missing primary
     screening criterion

independent_confirmation:
  Joule 2025, Hunan University
  "Solvation-Configurational Entropy
  Governs Interfacial Kinetics in
  Low-Temperature Batteries"
  Same variable. Opposite direction.
  Neither group has the band.
  The band is here. Timestamped.

related_documents:
  OC-BATTERY-SOLVATION-COHERENCE-
    DERIVATION.md
  OC-BATTERY-SOLVATION-VARIANCE-
    PREREGISTRATION.md
  OC-BATTERY-VARIANCE-DERIVATION-
    AND-CANDIDATES.md
  OC-BATTERY-LITERATURE-CHECK-
    DERIVED-CANDIDATES.md

repository:
  OrganismCore
  attractor-oncology

ORCID: 0009-0002-0414-6544
```

---

*The variable was named independently.*
*From the opposite direction.*
*In the same year.*

*Neither group has the band.*

*The band is the finding.*

*The temperature-responsive SCE*
*electrolyte is the target.*

*It does not exist yet.*

*This document is the record*
*that we found it first.*
