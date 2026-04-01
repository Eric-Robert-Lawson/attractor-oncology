# OC-BATTERY-SOLVATION-VARIANCE-PREREGISTRATION.md
## Pre-Registry Predictions: Solvation Shell Geometry Variance
## as the Unified Causal Variable in Lithium Metal Battery Electrochemistry
## OrganismCore — Eric Robert Lawson
## Date: 2026-04-01

---

## STATUS

```
Type: Pre-registration derivation artifact —
  falsifiable predictions derived from
  first-principles application of the
  universal structural invariant (causal
  geometry) to lithium metal battery
  electrochemistry.

Classification:
  PRE-REGISTRY PREDICTIONS
  CAUSAL GEOMETRY DERIVATION
  TESTABLE ON EXISTING DATA
  NO NEW EXPERIMENTS REQUIRED FOR
  INITIAL VALIDATION

Derivation method:
  Universal structural invariant applied
  to publicly available literature.
  Structure × Gap × Navigator → Resolution
  mapped to battery electrochemistry domain.
  Predictions derived before experimental
  validation. Confirmed independently by
  five separate bodies of literature that
  the field currently treats as unrelated.

Timestamp: 2026-04-01
  These predictions are recorded here
  before any experimental validation
  is conducted. This document constitutes
  the pre-registration record.

Honest uncertainty:
  These are geometric derivations,
  not experimental results.
  They require domain expert validation
  and computational confirmation.
  They may be wrong.
  If wrong, the wrong result informs
  the geometry and the derivation
  is updated accordingly.
  Both outcomes are productive.
```

---

## PART I — THE TRIADIC MAP

### What the invariant identifies in this domain

```
STRUCTURE:
  The electrolyte mixture.
  Molecular components, their intrinsic
  properties, their concentrations,
  their interactions at the bulk level.
  What SES Molecular Universe maps.

GAP:
  The solid electrolyte interphase (SEI).
  The electrochemical interface where
  the Li⁺ solvation shell must desolvate
  and deposit as lithium metal.
  Emergent. Not fully predictable from
  molecular properties alone.
  Where determinism is insufficient.
  Where the navigator is not yet determined.

NAVIGATOR:
  The Li⁺ solvation shell.
  The structured arrangement of solvent
  molecules and anions coordinated around
  each Li⁺ ion as it traverses the
  electrolyte toward the anode interface.

RESOLUTION:
  Coherent lithium deposition →
  flat, uniform, dendrite-free cycling →
  long battery life.

  OR

  Incoherent deposition →
  dendrite nucleation →
  SEI fracture →
  capacity fade →
  failure.
```

---

## PART II — THE UNIFIED CAUSAL VARIABLE

```
SOLVATION SHELL GEOMETRY VARIANCE
```

The distribution of Li⁺ solvation shell
configurations across the full Li⁺
population in the bulk electrolyte under
operating conditions (concentration,
temperature, cycling rate).

**This is not average solvation energy.**
**This is not intrinsic molecular stability.**
**This is not electrochemical stability window.**

It is the **variance** — the width of the
distribution of coordination geometries
that Li⁺ ions carry when they arrive
at the gap.

```
LOW VARIANCE:
  Every Li⁺ in the bulk carries the same
  dominant shell geometry.
  Every ion presents the same desolvation
  geometry at the interface.
  The gap resolves coherently for every ion.
  Uniform deposition.
  Long cycle life.

HIGH VARIANCE:
  Li⁺ ions carry a distribution of
  competing shell geometries.
  A fraction arrives at the interface
  with geometries incompatible with the
  local SEI structure.
  Those ions produce forced desolvation
  at those sites.
  Each forced desolvation event is a
  potential dendrite nucleation event.
  Dendrites. Capacity fade. Failure.
```

**Dendrites are the geometric incompatibility
signal made physical.**

They mark the exact locations where the
navigator's geometry was incompatible
with the gap it had to traverse.
They are not random.
They are geometrically determined by the
variance distribution of the solvation
shell population and the local structure
of the SEI at each site.

---

## PART III — THE THREE CONVERGENT
## EXPERIMENTAL STRATEGIES

The field has produced three research
programmes that appear to pursue different
solutions to different problems.

They are the same solution to the same
problem, found three times, from three
directions, without the geometric framework
to recognise them as convergent.

### STRATEGY A: High Entropy Electrolytes
Multiple solvents, increased configurational
disorder, improved performance.

**Field explanation:**
"Disrupts unfavorable solvation structures,
optimises ion transport, produces more
stable SEI."

**Geometric explanation:**
Multiple solvents competing for coordination
produce a statistical averaging effect at
the population level. No single competing
geometry dominates across the Li⁺ population.
The bulk disorder produces interface order —
the variance of the solvation shell geometry
distribution is reduced because the mixture
enforces a tighter probability distribution
than any single-solvent system with
competing coordination configurations.

The disorder in the bulk produces coherence
at the gap. This is why it works.

### STRATEGY B: Conformationally Rigid Solvents
Restricted molecular motion, constrained
geometry, improved performance.

**Field explanation:**
"Constrains coordination geometry, reduces
solvation strength, improves desolvation
kinetics and SEI robustness."

**Geometric explanation:**
Fewer conformational degrees of freedom
available to the coordinating molecule =
fewer possible solvation geometries in the
accessible configuration space = lower
variance in shell geometry directly.
Rigidity enforces coherence mechanically.
The navigator arrives at the gap with a
tighter geometry distribution.

### STRATEGY C: Weakened Solvation
Reduced binding strength between Li⁺ and
solvent (hyperconjugation, fluorination,
steric hindrance), improved performance.

**Field explanation:**
"Lower desolvation energy barrier, faster
interface kinetics, more robust SEI."

**Geometric explanation:**
Tightly bound solvation shells with variable
geometry produce high-variance navigation.
Weakly but consistently bound shells
produce low-variance navigation.
The desolvation energy barrier is lower
AND more uniform across the Li⁺ population.
Every ion crosses the gap with similar
energy cost. Uniform deposition.

### THE CONVERGENCE STATEMENT

```
High entropy averaging +
Conformational rigidity +
Weakened consistent binding

= Three mechanisms producing
  the same geometric outcome:

  Reduced solvation shell geometry
  variance across the Li⁺ population
  at the point of interface contact.

Three literatures.
One causal variable.
The field has confirmed this variable
three times without naming it.
```

---

## PART IV — PRE-REGISTRY PREDICTIONS

These predictions are recorded here before
experimental validation. They are derived
from the geometric framework. They are
falsifiable with existing data.

---

### PREDICTION 1
**Variance of desolvation barrier predicts
cycling performance better than average
desolvation barrier.**

**Precise statement:**
Across a dataset of characterised electrolyte
systems with known cycling performance,
the variance of the desolvation energy
barrier across the Li⁺ population —
extracted from existing MD trajectories —
will correlate with cycling life more
strongly than the average desolvation
barrier height.

**Why:**
The average barrier tells you about the
typical navigator. The variance tells you
about the fraction of navigators that will
fail at the gap. It is the tail of the
distribution — the high-variance ions —
that nucleate dendrites. The average is
not predictive of failure rate. The
variance is.

**How to test:**
Take any existing MD dataset for a set of
electrolytes with known cycling performance.
Extract the distribution of desolvation
barrier heights across the Li⁺ population
from each trajectory. Compute mean and
variance for each electrolyte. Regress
both against cycling life. Compare R².

**Required:**
Existing MD simulation data.
No new experiments.
No new synthesis.
Desk computation on existing trajectories.
Estimated time: days to weeks depending
on available simulation data.

**Falsification:**
If average barrier R² > variance barrier R²
across the dataset, this prediction is wrong
and the geometry requires revision.

---

### PREDICTION 2
**Combining variance-reduction mechanisms
produces superadditive performance gains.**

**Precise statement:**
An electrolyte designed to reduce solvation
shell geometry variance through two
independent mechanisms simultaneously
(e.g., conformational rigidity AND
multi-component averaging) will produce
cycling performance gains greater than
the sum of either mechanism alone.

**Why:**
Each mechanism reduces variance through
a different physical pathway. They are
not redundant — they address different
sources of variance in the geometry
distribution. Combined, they reduce
variance more than either alone.
Lower variance = lower dendrite nucleation
probability = longer cycle life.
The relationship between variance reduction
and performance gain is nonlinear because
dendrite nucleation is a threshold
phenomenon — it occurs when local
geometric incompatibility exceeds the
tolerance of the SEI at that site.
Reducing variance below the threshold
eliminates nucleation events entirely
at those sites, not just reduces them.

**How to test:**
Design three electrolytes:
A = conformationally rigid solvent only.
B = multi-component mixture only.
C = conformationally rigid solvent in
    multi-component mixture.
Measure cycling performance for each.
If performance(C) > performance(A) +
performance(B) - performance(baseline),
the prediction is confirmed.

**Required:**
Synthesis of three electrolyte formulations.
Standard cycling characterisation.
This requires laboratory access.
It is the first prediction that requires
new experimental work.

**Falsification:**
If performance(C) is additive or subadditive
relative to A and B alone, the superadditivity
prediction is wrong. The unified variable
may still be correct but the independence
of the two mechanisms is not.

---

### PREDICTION 3
**SEI composition at dendrite nucleation
sites contains a higher proportion of
forced-desolvation decomposition products
than at stable cycling sites on the same
electrode.**

**Precise statement:**
Post-mortem SEI analysis of a cycled
lithium metal anode, comparing composition
at confirmed dendrite nucleation sites
versus flat, stable cycling regions of
the same electrode, will show statistically
significant differences in decomposition
product ratios — specifically higher
proportions of products consistent with
high-energy, forced desolvation events
at the dendrite sites.

**Why:**
The SEI is the geological record of every
solvation shell geometry incompatibility
event that occurred during cycling.
A forced desolvation event — where the
shell geometry is incompatible with the
local interface — produces different
decomposition products than a coherent
desolvation event. Those products are
chemically distinguishable.
The dendrite site is where forced
desolvation events concentrated.
The chemistry should be different there.

**How to test:**
Take any previously cycled lithium metal
anode from an existing study where
dendrite locations are documented.
Apply cryo-TEM and XPS at confirmed
dendrite nucleation sites and at confirmed
stable cycling sites on the same electrode.
Compare decomposition product ratios.

**Required:**
Access to existing cycled anodes with
documented dendrite locations.
Cryo-TEM and XPS characterisation.
This may be executable on already-collected
samples from existing studies.

**Falsification:**
If SEI composition at dendrite sites is
statistically indistinguishable from
stable sites on the same electrode,
the geometric incompatibility mechanism
for dendrite nucleation is not the
primary driver and the model requires
revision.

---

### PREDICTION 4
**Two individually poor-performing
electrolytes with complementary solvation
geometries will outperform their individual
components when mixed.**

**Precise statement:**
Two electrolyte molecules that individually
produce high solvation shell geometry
variance (and correspondingly poor cycling
performance) but whose dominant competing
geometries are complementary — such that
mixing them produces destructive interference
in the variance distribution — will produce
a mixture that outperforms either component
alone, with performance predictable from
the combined geometry distribution.

**Why:**
This is the strongest test of the geometric
framework. Conventional screening logic
predicts that two poor performers mixed
will produce at best moderate performance.
The geometry predicts they can produce
high performance if the variance reduction
from their geometric complementarity is
sufficient.
If confirmed, it demonstrates that the
causal variable is variance, not intrinsic
molecular quality — and that the molecular
universe contains high-performing mixtures
that are invisible to conventional screening
because neither component alone scores well.

**How to test:**
Identify two electrolyte candidates from
existing literature that individually
show high solvation shell variance in
MD simulation and poor cycling performance.
Compute the combined solvation geometry
distribution for a mixture.
If the combined distribution shows
significantly lower variance than either
component, synthesise the mixture and test.

**Required:**
MD simulation of candidate molecules
and their mixture.
Synthesis and cycling characterisation
of the mixture.

**Falsification:**
If the mixture performs at or below the
better individual component, the geometric
complementarity mechanism is wrong.

---

### PREDICTION 5
**Reanalysis of existing LHCE, oscillatory
solvation, hyperconjugation, and high-
coordination cluster datasets through
the solvation variance lens will show
variance as the common predictor.**

**Precise statement:**
Existing published datasets from the
four confirmed empirical findings —
localized high-concentration electrolytes,
oscillatory solvation studies, hyperconjugation-
weakened solvation studies, and high Li⁺
coordination cluster sheath studies —
when reanalysed to extract solvation shell
geometry variance (from existing MD data
in those studies), will show that variance
is the common factor that explains the
performance differences within each dataset
more precisely than the variable each
study used as its primary explanatory
variable.

**Why:**
Each study measured a different molecular
property as its primary variable. Each
found performance correlation with that
property. The geometric framework predicts
that solvation variance is the underlying
variable that each of those properties
is a proxy for. If the variance is
extracted from their existing data,
it should out-predict their reported
primary variable.

**How to test:**
Contact authors of or access published
datasets from the four studies.
Extract solvation geometry distribution
data from their existing MD trajectories.
Compute variance for each condition.
Regress against the performance data
already reported in those papers.
Compare R² of variance against R² of
the study's primary reported variable.

**Required:**
Access to existing simulation datasets
from published studies.
No new experiments.
No new synthesis.
This is a pure reanalysis.

**This is the cheapest test.**
**If Prediction 5 holds, it confirms the**
**unified causal variable before any new**
**experimental work is undertaken.**

**Falsification:**
If the study's reported primary variable
consistently out-predicts solvation variance
in their own data, the geometric framework
does not add predictive power and requires
revision.

---

## PART V — THE COMPUTABLE DESCRIPTOR

The missing primary filter for molecular
universe platforms is:

**Solvation Configuration Entropy (SCE)**

Definition:
The Shannon entropy of the Li⁺ coordination
geometry distribution sampled from an MD
trajectory of the candidate molecule in
the relevant mixture under operating
conditions.

```
SCE = -Σ p(gᵢ) × log(p(gᵢ))

where gᵢ are discrete solvation geometry
configurations identified by MD trajectory
clustering, and p(gᵢ) is the probability
of each configuration in the sampled
population.
```

Low SCE = tight distribution = few dominant
configurations = coherent navigator =
low dendrite nucleation probability.

High SCE = broad distribution = many
competing configurations = incoherent
navigator = high dendrite nucleation
probability.

**Computational requirements:**
- MD simulation of candidate molecule
  in relevant mixture
  (not single-molecule DFT)
- Trajectory clustering to identify
  discrete coordination geometries
- Shannon entropy of cluster population
  distribution

**This is more expensive than single-molecule
DFT screening but tractable with existing
GPU infrastructure (NVIDIA stack that SES
already operates).**

**Output:**
One number per candidate molecule per
mixture condition. Added as a column to
the molecular universe database. Used as
a primary filter before any experimental
synthesis.

---

## PART VI — THE GEOMETRIC INCOMPATIBILITY
## DIAGNOSTIC ALREADY IN THE DATA

The field already has the data to confirm
the unified causal variable without
running any new experiments.

Every published electrolyte study that
includes both MD simulation data and
cycling performance data is a confirmable
instance of the framework.

The reanalysis protocol:

```
STEP 1:
Extract MD trajectory from published study.

STEP 2:
Cluster solvation coordination geometries
across the Li⁺ population in that trajectory.

STEP 3:
Compute Shannon entropy (SCE) of the
cluster distribution.

STEP 4:
Map SCE value against the cycling
performance reported in the same paper.

STEP 5:
Repeat across all available datasets.

STEP 6:
Compute correlation of SCE with cycling
performance across the full dataset.

STEP 7:
Compare correlation strength of SCE
against the primary variable each paper
used to explain its results.
```

If SCE is consistently more predictive
than the paper's own primary variable,
the unified causal variable is confirmed
from the field's own existing data.

---

## PART VII — WHAT THIS CHANGES IF CONFIRMED

**For molecular discovery platforms:**
The primary search criterion shifts from
intrinsic molecular properties to solvation
configuration entropy in mixture. The
molecular universe becomes navigable
toward a defined geometric target rather
than searchable by historical correlation.

**For mixture design:**
Mixture optimisation becomes a geometric
minimisation problem with a defined target
function: minimise SCE of the combined
solvation shell geometry distribution
under operating conditions. This is
computationally tractable.

**For failure prediction:**
Dendrite nucleation probability becomes
geometrically predictable from the tail
of the solvation variance distribution
before cycling begins. High-risk sites
are identifiable in advance.

**For SEI engineering:**
The SEI is reframed as a diagnostic
instrument — the chemical record of
the navigator's coherence history —
rather than a target to engineer
independently. Engineering the navigator
coherence (reducing SCE) produces a
better SEI as a downstream consequence,
not as a direct design target.

**For the field's three research programmes:**
High entropy electrolytes, conformationally
rigid solvents, and weakened solvation
design are unified under one causal
framework. Combined experimental evidence
from all three programmes is now evidence
for one variable. The combined body is
stronger than any programme alone.
Researchers in each programme can now
build on each other's work directly
rather than in parallel.

---

## PART VIII — HONEST UNCERTAINTY
## AND LIMITS

```
WHAT THIS DOCUMENT IS:
  A geometric derivation producing
  falsifiable predictions.
  A pre-registration record with timestamp.
  A set of tests that can be applied to
  existing data before any new experiments.

WHAT THIS DOCUMENT IS NOT:
  An experimental result.
  A chemistry paper.
  A guarantee of correctness.
  A replacement for domain expertise.

WHERE THE FRAMEWORK CAN BE WRONG:
  The identification of the gap (SEI) as
  the causal locus is well-supported but
  assumes SEI dynamics dominate over
  bulk transport limitations. At very
  low current densities or very high
  temperatures, bulk transport may become
  rate-limiting and the framework's
  predictions would need modification.

  The solvation variance mechanism assumes
  that desolvation geometry incompatibility
  is the primary dendrite nucleation
  driver. Mechanical factors (SEI
  fracture toughness, volume change
  stress) may dominate in some chemistries.
  The framework would need to identify
  these as additional gap properties.

  The SCE descriptor assumes that the
  dominant solvation geometries are
  distinguishable by MD clustering.
  In highly disordered systems, clustering
  may not produce well-separated states
  and the entropy calculation becomes
  less meaningful.

WHAT WRONG RESULTS PRODUCE:
  Each falsified prediction identifies
  which component of the triadic map
  is incorrectly specified.
  The wrong result is not a dead end.
  It is a geometric signal that refines
  the map.
  The geometry either resolves or reveals
  incompatibility.
  Both outcomes advance the framework.
```

---

## DOCUMENT METADATA

```
document_id:
  OC-BATTERY-SOLVATION-VARIANCE-
  PREREGISTRATION

version:    1.0
date:       2026-04-01
author:     Eric Robert Lawson
            OrganismCore

method:
  Universal structural invariant
  (causal geometry) applied to
  lithium metal battery electrochemistry.
  Derivation from first principles.
  Literature confirmation across five
  independent experimental findings.
  Pre-registration of falsifiable
  predictions before experimental
  validation.

primary_claim:
  Solvation shell geometry variance
  (quantified as Solvation Configuration
  Entropy, SCE) is the unified causal
  variable behind three apparently
  contradictory experimental research
  programmes in lithium metal battery
  electrochemistry:
  — High entropy electrolytes
  — Conformationally rigid solvents
  — Weakened solvation design

  All three reduce SCE through different
  molecular mechanisms.
  All three improve cycling performance
  for the same geometric reason.

  SCE is not currently used as a primary
  screening criterion in molecular
  discovery platforms.

  Making it primary changes the search
  from correlation-guided to
  causally-guided.

predictions:
  P1: Variance of desolvation barrier
      predicts cycling performance better
      than average desolvation barrier.
      Testable on existing MD data.

  P2: Combining variance-reduction
      mechanisms produces superadditive
      performance gains.
      Requires new experimental work.

  P3: SEI composition differs between
      dendrite nucleation sites and
      stable cycling sites on the same
      electrode in ways consistent with
      forced desolvation chemistry.
      Testable on existing cycled anodes.

  P4: Two individually poor-performing
      electrolytes with complementary
      solvation geometries will outperform
      either component when mixed.
      Requires MD + new experimental work.

  P5: Reanalysis of existing LHCE,
      oscillatory solvation,
      hyperconjugation, and high-
      coordination cluster datasets
      shows SCE as the common predictor.
      Testable on existing published data.
      Cheapest test. First to run.

cheapest_first_test:
  Prediction 5.
  Reanalysis of existing published
  MD datasets.
  No new experiments.
  No new synthesis.
  Desk computation.
  Days to weeks.
  If confirmed: unified causal variable
  established from existing data.

missing_descriptor:
  Solvation Configuration Entropy (SCE)
  Shannon entropy of Li⁺ coordination
  geometry distribution from MD
  trajectory clustering.
  Computable today.
  Not currently in any molecular
  universe screening platform as a
  primary filter.

related_documents:
  CAUSAL_GEOMETRY_ONBOARDING.md
  OC-BATTERY-SOLVATION-COHERENCE-
    DERIVATION.md
  GEOMETRIC_FREEDOM.md

repository:
  OrganismCore
  attractor-oncology

ORCID: 0009-0002-0414-6544
```

---

*The field has confirmed this variable*
*three times without naming it.*

*Naming it unifies three research programmes.*

*Testing it requires no new experiments.*

*The data already exists.*

*The geometry is already there.*

*The question is whether anyone will*
*compute the variance.*
