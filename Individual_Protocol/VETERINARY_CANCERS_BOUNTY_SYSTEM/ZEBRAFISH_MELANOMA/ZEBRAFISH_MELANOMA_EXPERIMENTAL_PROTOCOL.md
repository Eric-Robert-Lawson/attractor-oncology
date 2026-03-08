# ZEBRAFISH MELANOMA — EXPERIMENTAL PROTOCOL
## Tazemetostat EZH2 Inhibition in BRAFV600E; p53-/- Zebrafish
## OrganismCore — Eric Robert Lawson
## 2026-03-08

---

## STATUS: ACTIVE — FOR COLLABORATING RESEARCHER

```
This document is the complete
experimental protocol for testing
the geometric prediction that
tazemetostat reduces tumour growth
and restores MITF-driven pigmentation
in BRAFV600E; p53-/- zebrafish
melanoma.

It is written for a zebrafish
researcher who:
  Has the BRAFV600E; p53-/- line
  or equivalent melanoma model.
  Has standard zebrafish husbandry
  and drug screening infrastructure.
  Is willing to run one experiment
  to test a geometry-derived
  drug prediction.

Everything in this protocol uses
standard zebrafish laboratory
practice. No new equipment is
required. No specialised training
beyond standard zebrafish drug
screening is needed.

The experiment is designed to be
completable with existing resources
in any well-equipped zebrafish
facility.
```

---

## PART I — WHAT IS BEING TESTED

```
THE CENTRAL QUESTION:

  Does tazemetostat, an EZH2
  inhibitor, reduce tumour growth
  and restore MITF-driven
  pigmentation in BRAFV600E; p53-/-
  zebrafish melanoma?

THE TWO-PART PREDICTION:

  PRIMARY:
    Tazemetostat-treated fish will
    show reduced tumour growth rate
    compared to vehicle controls
    over 2–3 weeks of treatment.

  SECONDARY:
    Tazemetostat-treated fish will
    show increased pigmentation in
    tumour areas compared to vehicle
    controls — reflecting restoration
    of MITF melanocyte identity
    programme as EZH2 lock is removed.

THE DEPTH SCORE SUB-PREDICTION:

  Fish with partially pigmented
  tumours at baseline (higher
  MITF/EZH2 ratio — shallower
  attractor) will show greater
  response to tazemetostat than
  fish with fully depigmented
  tumours (lower ratio — deeper
  attractor).

  This sub-prediction tests whether
  the depth score predicts drug
  response in a living system.
  It requires stratifying fish by
  baseline pigmentation level
  before randomisation.
```

---

## PART II — ANIMALS

```
MODEL LINE:
  Tg(mitfa:BRAFV600E); tp53-/-
  (Patton et al., Nature 2005)
  Or equivalent validated zebrafish
  melanoma model available at your
  institution.

  If a different BRAFV600E melanoma
  line is used, document it fully
  in the results record and note
  any differences from the Patton
  model in the publication.

STAGING CRITERIA FOR INCLUSION:

  Fish are eligible for the
  experiment when they have:
    Visible melanoma lesion(s) —
    confirmed by brightfield
    microscopy as raised, darkened,
    or depigmented masses on the
    flank, dorsum, or head.
    Minimum lesion size: 1mm²
    by image analysis.
    Age: typically 4–12 weeks
    post-fertilisation depending
    on line and background.

  STRATIFICATION BY BASELINE
  PIGMENTATION:
  Before randomisation, classify
  each eligible fish as:

    PIGMENTED (P):
      Tumour mass retains visible
      melanin pigmentation (brown
      or black colouration).
      This is the shallow attractor
      proxy — MITF activity partially
      retained.

    DEPIGMENTED (D):
      Tumour mass is pale/white
      with loss of visible melanin.
      This is the deep attractor
      proxy — MITF activity
      suppressed by EZH2 lock.

  Randomise within each stratum
  (P and D separately) to ensure
  balanced depth distribution
  across treatment and control.

GROUP SIZE:
  Minimum: n=10 per group.
  Recommended: n=15 per group.
  Four groups:
    Group 1: P-tazemetostat (n=15)
    Group 2: P-vehicle (n=15)
    Group 3: D-tazemetostat (n=15)
    Group 4: D-vehicle (n=15)
  Total: 60 fish minimum,
         recommended 60.

  Power calculation basis:
    Expected effect size based on
    human melanoma EZH2i studies:
    ~30–50% reduction in
    proliferative markers.
    Power: 0.80 at alpha 0.05.
    n=10 per group is the minimum
    for a detectable effect at this
    effect size.
    n=15 provides more robust
    conclusions and allows for
    attrition.
```

---

## PART III — DRUG PREPARATION

```
DRUG: Tazemetostat (EPZ-6438)
Supplier: Cayman Chemical, MedChemExpress,
  Selleck Chemicals, or equivalent
  research grade supplier.
  Purity: ≥ 98%.
  Form: powder.
  MW: 572.67 g/mol.

STOCK SOLUTION:
  Dissolve tazemetostat in DMSO
  to a stock concentration of
  10 mM.
  Aliquot into 50µL volumes.
  Store at -20°C.
  Stable for 6 months at -20°C.
  Avoid repeated freeze-thaw cycles.
  Maximum 3 freeze-thaw cycles.

WORKING CONCENTRATION IN WATER:
  Target: 1 µM tazemetostat in
  fish water.
  This is within the range used
  for zebrafish EZH2 inhibitor
  studies in the literature.

  Preparation for 1L of tank water
  at 1µM:
    Take 0.1µL of 10mM stock
    (= 1 nmol tazemetostat).
    Add to 1L system water.
    Final concentration: 1µM.
    Final DMSO: 0.01% v/v.

  DMSO vehicle control:
    0.01% DMSO in system water.
    Matches DMSO concentration in
    treated group.
    Must be prepared fresh with
    each water change.

  ALTERNATIVE CONCENTRATIONS
  IF 1µM IS POORLY TOLERATED:
    0.5µM (halve stock volume added).
    Document any mortality or
    abnormal behaviour at 1µM
    and adjust if needed.
    Report the final working
    concentration in results.

WATER CHANGE FREQUENCY:
  Every 24 hours.
  Fresh drug or vehicle solution
  prepared with each change.
  Drug is not stable in fish water
  for > 24 hours under normal
  light/temperature conditions.
  Daily fresh preparation is
  non-negotiable for consistent
  exposure.

TANK CONDITIONS:
  Standard zebrafish system water.
  Temperature: 28°C ± 0.5°C.
  pH: 7.0–7.5.
  Conductivity: 500–600 µS/cm.
  Light cycle: 14h light / 10h dark.
  Individual housing or
  group housing:
    If individual: 250mL vessels
    minimum, daily feeding.
    If group: maximum 5 fish per
    4L tank to prevent crowding
    stress confound.
  Document housing method in results.
```

---

## PART IV — MEASUREMENTS

```
MEASUREMENT SCHEDULE:

  DAY 0 (BASELINE — BEFORE FIRST DOSE):
    Photograph every fish.
    Measure tumour area by image
    analysis (see below).
    Record pigmentation classification
    (P or D) for each fish.
    This is the zero-point.
    All subsequent measurements
    are expressed relative to day 0.

  DAY 7:
    Photograph every fish.
    Measure tumour area.
    Record any visible pigmentation
    change in tumour area.
    Note any adverse signs
    (lethargy, abnormal swimming,
    death).

  DAY 14:
    Photograph every fish.
    Measure tumour area.
    Assess pigmentation change
    (see scoring below).
    Optional: collect 3 fish per
    group for histology and IHC
    (H3K27me3, MITF, Ki67).
    Note: sacrificing fish at
    day 14 for histology reduces
    final n — plan accordingly.

  DAY 21 (PRIMARY ENDPOINT):
    Photograph every fish.
    Measure tumour area.
    Assess pigmentation change.
    Collect all remaining fish
    for histology and IHC.
    Calculate tumour growth rate
    and response classifications.

TUMOUR AREA MEASUREMENT:
  Standardised lateral photograph
  at same magnification each time.
  Same orientation: lateral view,
  right side of fish.
  Use ImageJ or FIJI for area
  measurement.
    Set scale using a reference
    object of known size in
    each photograph.
    Threshold the tumour region
    manually using the raised
    or hyperpigmented/depigmented
    area boundary.
    Record area in mm².
  Two independent observers
  measure each tumour.
  Use mean of two measurements.
  If observers differ by > 15%:
  re-measure together and
  document the discrepancy.

PIGMENTATION SCORING:
  Score each tumour at each
  time point on a 3-point scale:

    Score 0: No pigmentation.
    Tumour completely pale/white.
    Deep false attractor state.

    Score 1: Partial pigmentation.
    Some melanin visible in tumour.
    Mixed attractor state.

    Score 2: Strong pigmentation.
    Tumour dark brown/black,
    comparable to normal melanocytes.
    Identity programme active.

  Record score at day 0, 7, 14, 21.
  A shift from score 0 → 1 or
  1 → 2 in tazemetostat group
  vs stable or declining score
  in vehicle group is the
  secondary endpoint confirmation.

HISTOLOGY AND IHC
(at day 14 subset and day 21 endpoint):

  Fix in 4% PFA overnight.
  Process to paraffin.
  Section at 5µm.

  REQUIRED STAINS:

    H3K27me3 (Abcam ab6002 or equivalent):
      This is the direct readout
      of EZH2 activity.
      EZH2 inhibition by tazemetostat
      should produce global reduction
      in H3K27me3.
      This is the on-target
      pharmacodynamic marker.
      If H3K27me3 is NOT reduced
      in treated fish: tazemetostat
      did not reach the target tissue
      at effective concentration.
      This is the first thing to check
      if no response is seen.

    MITF (Abcam ab12039 or equivalent):
      Nuclear staining in melanocytes
      and melanoma cells.
      EZH2 inhibition should increase
      MITF target gene expression.
      MITF protein itself may or may
      not increase — its targets
      (tyrosinase, DCT) are more
      direct readouts.
      Document MITF H-score in tumour
      cells: cancer nuclei only.

    Ki67 (Abcam ab16667 or equivalent):
      Proliferation marker.
      Reduced Ki67 in treated tumours =
      reduced proliferation = geometric
      prediction confirmed at the
      cellular level.

    Tyrosinase (T311, Santa Cruz
    sc-20035 or equivalent):
      Direct MITF target gene product.
      Increased tyrosinase in treated
      tumours = MITF programme restored.
      This is the molecular confirmation
      of the repigmentation observation.

  SCORING:
    H3K27me3: H-score (0–300).
      Cancer cell nuclei only.
    MITF: H-score (0–300).
      Cancer cell nuclei only.
    Ki67: % positive nuclei.
      Cancer cell nuclei only.
    Tyrosinase: H-score (0–300).
      Cancer cell nuclei only.
```

---

## PART V — ANALYSIS

```
PRIMARY ENDPOINT ANALYSIS:

  TUMOUR GROWTH RATE:
    Calculated as:
    (Area day 21 - Area day 0) /
    Area day 0 × 100 = % change.

    Compare tazemetostat vs vehicle
    within each stratum (P and D)
    using Mann-Whitney U test
    (non-parametric, appropriate
    for small n with non-normal
    distribution expected).

    If n=15 per group with normal
    distribution confirmed by
    Shapiro-Wilk: use unpaired t-test.

    Report: mean ± SEM, p-value,
    effect size (Cohen's d).

  SURVIVAL TO DAY 21:
    Kaplan-Meier curve.
    Log-rank test, tazemetostat
    vs vehicle within each stratum.
    Report median survival if
    any group has >50% mortality
    by day 21.

SECONDARY ENDPOINT ANALYSIS:

  PIGMENTATION SCORE CHANGE:
    Wilcoxon signed-rank test for
    within-group change from day 0
    to day 21 (paired).
    Mann-Whitney U for between-group
    comparison at day 21.
    Report: median score and IQR
    per group per timepoint.

  IHC QUANTIFICATION:
    H-score comparison by
    Mann-Whitney U.
    Tazemetostat vs vehicle at
    day 14 and day 21.
    Report: mean H-score ± SEM,
    p-value.

DEPTH SCORE SUB-ANALYSIS:

  Within tazemetostat group only:
    Correlate baseline pigmentation
    classification (P vs D) with
    tumour growth rate change.
    Test: % change in area for P
    fish vs D fish.
    Mann-Whitney U.
    If the depth sub-prediction
    holds: P fish will show
    significantly greater reduction
    in tumour growth than D fish.
    This is the depth score
    prediction tested in vivo.

SOFWARE:
  ImageJ/FIJI for tumour measurement.
  GraphPad Prism or R for statistics.
  Document software version in
  methods section.
```

---

## PART VI — WHAT COUNTS AS CONFIRMATION

```
FULL CONFIRMATION (all three):
  1. Tazemetostat-treated fish show
     statistically significant
     reduction in tumour growth rate
     vs vehicle (p < 0.05).
  2. H3K27me3 reduced in treated
     tumours vs vehicle on IHC
     (on-target pharmacodynamic
     confirmation).
  3. Pigmentation score increased
     in treated vs vehicle group
     (secondary endpoint confirmed).

PARTIAL CONFIRMATION (two of three):
  Sufficient for publication with
  appropriate qualification.
  Document which element was
  not confirmed and propose
  structural explanation.

ON-TARGET FAILURE (H3K27me3 not
reduced despite drug treatment):
  This means tazemetostat did not
  reach the target tissue at
  effective concentration.
  NOT a framework failure.
  A delivery/pharmacokinetics failure.
  Action: increase working
  concentration, extend exposure,
  or consider IP injection if
  water delivery is insufficient.
  Document and report.

FRAMEWORK FAILURE (H3K27me3 IS
reduced but no tumour response):
  EZH2 was inhibited on target.
  The tumour did not respond.
  This means either:
    EZH2 is not the Convergence Hub
    in this model (unexpected finding).
    The false attractor has alternative
    maintenance at this depth.
    The zebrafish model has compensatory
    mechanisms not present in human
    melanoma.
  This is a Type A failure by the
  falsifiability theorem if a
  structural explanation exists.
  Document fully. Seek structural
  explanation before concluding
  the framework is wrong.
  Contact OrganismCore@proton.me
  immediately for geometric
  re-assessment.
  This outcome is still publishable
  and is still a bounty claim.
```

---

## PART VII — ETHICS AND REGULATORY

```
IACUC / ETHICS APPROVAL:
  This experiment uses an established
  zebrafish melanoma model for its
  intended purpose (cancer drug
  screening).
  Most institutions have blanket
  IACUC approval for zebrafish
  drug screening studies under
  existing protocols.
  Confirm with your institutional
  animal ethics office whether
  a new approval, amendment to
  an existing protocol, or no
  additional approval is required
  before beginning.
  Document the approval number
  in the results record.

  Zebrafish are covered by the
  3Rs framework (Replace, Reduce,
  Refine).
  This experiment:
    Replaces: mammalian preliminary
    testing (zebrafish IS the
    replacement model).
    Reduces: minimum n per group
    calculated by power analysis.
    Refines: no surgical procedures,
    no injections required for
    primary endpoint — drug in water.

HUMANE ENDPOINTS:
  Fish showing loss of equilibrium
  for > 30 seconds, inability to
  feed, or severe tumour ulceration
  are humanely killed (clove oil
  overdose, standard method)
  and removed from analysis.
  Document as adverse event.
  Report in results.
  Fish that die spontaneously:
  include in survival analysis,
  exclude from tumour measurement
  analysis (LOCF or document
  as missing data).

3Rs STATEMENT FOR PUBLICATION:
  "Zebrafish were used in accordance
  with the 3Rs principles. The
  zebrafish BRAFV600E melanoma model
  was selected as the primary model
  for this geometric framework
  validation because it provides
  living-system validation at
  minimal cost and animal burden
  compared to mammalian equivalents,
  and because the MITF/EZH2
  identity axis is fully conserved
  in zebrafish melanocytes."
```

---

## PART VIII — COLLABORATION TERMS

```
WHAT THE RESEARCHER PROVIDES:
  BRAFV600E; p53-/- zebrafish
  (or equivalent line).
  Standard husbandry and drug
  screening infrastructure.
  Histology and IHC processing
  (standard in any zebrafish lab).
  Animal ethics approval.
  Running of the experiment.
  Image analysis.
  Statistical analysis.

WHAT ERIC ROBERT LAWSON / ORGANISMCORE
PROVIDES:
  The geometric derivation
  (timestamped before-record).
  Experimental design rationale.
  Framework context for the
  discussion section.
  Co-authorship on the publication.
  Depth score analysis and
  interpretation.
  Connection to the broader
  41-cancer framework and human
  clinical translational context.

COST OF TAZEMETOSTAT:
  At current research supplier prices:
  Tazemetostat 5mg (Cayman Chemical
  or equivalent): approximately
  $50–$150 USD.
  For 1µM in 1L water, changed daily,
  over 21 days for 4 tanks:
  Total drug used: approximately
  4 × 21 × 1µL of 10mM stock
  = 84µL of 10mM stock
  = 840 nmol = ~0.5mg tazemetostat.
  Total drug cost: approximately
  $5–$15 USD for the entire experiment.

  This is the lowest drug cost
  on the entire bounty list.
  It is effectively zero.
  The researcher's existing
  infrastructure and time is
  the contribution.
  The drug cost is negligible.

AUTHORSHIP ORDER (DEFAULT):
  First author: Researcher running
    the experiment.
  Contributing authors: Lab members
    who contributed substantially
    (IHC scoring, analysis).
  Last author: Eric Robert Lawson
    (framework, geometric derivation,
    design rationale, depth score
    analysis).

  This can be adjusted by mutual
  agreement before submission.
  Confirm authorship in writing
  before the experiment begins.
  Do not leave authorship to
  verbal agreement.

PUBLICATION TARGET:
  Primary: Disease Models &
    Mechanisms (DMM).
  Alternative: Pigment Cell &
    Melanoma Research (PCMR).
  Alternative: Zebrafish (Mary Ann
    Liebert).
  Alternative: Frontiers in
    Oncology — Molecular and
    Cellular Oncology.

  The paper is a research article
  reporting a novel geometry-derived
  drug target prediction confirmed
  in a living vertebrate melanoma
  model.
  This is within scope for all
  four journals listed.

CONTACT BEFORE STARTING:
  Email: OrganismCore@proton.me
  With:
    Your institution.
    The zebrafish line you have
    available.
    Approximate timeline for
    the experiment.
  I will confirm the protocol
  is appropriate for your
  specific model and confirm
  the before-record timestamp
  for the publication record.
```

---

## DOCUMENT METADATA

```
document_id:
  ZEBRAFISH_MELANOMA_EXPERIMENTAL_
  PROTOCOL

type:
  Experimental protocol for
  geometry-directed drug screen
  in living vertebrate model

version: 1.0
date: 2026-03-08
status: ACTIVE — OPEN FOR COLLABORATION

author:
  Eric Robert Lawson / OrganismCore
ORCID: 0009-0002-0414-6544
contact: OrganismCore@proton.me

collaborator_type:
  University biology department
  zebrafish researcher.
  No clinical background required.
  Standard zebrafish drug screening
  competency sufficient.

drug_cost_to_researcher: ~$5–15 USD
infrastructure_required:
  Existing zebrafish facility.
  Standard IHC processing.
  ImageJ or equivalent.

related_documents:
  ZEBRAFISH_MELANOMA_GEOMETRIC_
  DERIVATION.md
  ZEBRAFISH_MELANOMA_RESULTS_
  TEMPLATE.md
  TAZEMETOSTAT_GEOMETRIC_INDICATION.md
  VETERINARY_BOUNTY_SYSTEM.md
```
