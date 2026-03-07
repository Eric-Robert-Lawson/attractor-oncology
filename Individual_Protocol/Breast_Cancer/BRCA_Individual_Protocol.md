# BRCA INDIVIDUAL PATIENT GEOMETRIC ANALYSIS PROTOCOL
## A Principles-First Derivation Protocol for Individual Breast Cancer Cases
## OrganismCore — Eric Robert Lawson
## 2026-03-07

---

## PREAMBLE — THE EPISTEMOLOGICAL COMMITMENT

```
This protocol exists because of one principle:

  The geometry must reveal itself.
  You do not go in knowing the answer.
  You go in knowing the method.

The population analysis found FOXA1 and EZH2
as the primary axis of breast cancer geometry.
That was a result.
It was not assumed going in.
The top_mover_scan found the top movers.
The saddle point scan found the saddle point.
The identity TF direction test found the
direction. The lock emerged from the data.

For an individual patient, the same commitment
applies. The patient's data may confirm the
FOXA1/EZH2 axis. It may reveal a different
dominant signal in this specific case — a
deeper lock, a secondary axis, an unexpected
switch gene suppression.

Whatever the data shows is what the analysis
reports.

This is the difference between:
  "Your FOXA1 is X and your EZH2 is Y"
  (imposing a known axis onto data)
and:
  "The geometry of your data reveals this
  structure, and here is what that structure
  means" (letting the geometry speak).

The first approach is a lookup table.
The second approach is this protocol.

The population map is the coordinate system.
The patient's data is placed into it.
The geometry is what the placement reveals.
Not what was assumed before the placement.
```

---

## PART I — WHAT THIS PROTOCOL IS NOT

```
This protocol is NOT:

  A lookup table.
  You do not match the patient's ER status
  to a known subtype and report the subtype's
  known biology. That is what the clinical
  report already does.

  A FOXA1/EZH2 ratio calculator.
  FOXA1/EZH2 is what the population analysis
  found. It is a result. For this patient,
  the primary axis is whatever the geometry
  reveals. It will often be FOXA1/EZH2 because
  that axis was confirmed across 7,500 patients.
  But the protocol does not assume this in
  advance. It discovers it from the data.

  A subtype classifier.
  The clinical team has already classified
  the subtype. That classification is one
  of the inputs. It is not the output.
  The output is the geometry underneath the
  subtype classification — what is driving it,
  how deep it is, what the lock is, and what
  the drug target geometry looks like.

  A prognosis tool.
  This protocol does not predict survival.
  It describes geometric position and
  geometric structure. What that structure
  implies for treatment vulnerability is
  stated as a geometric observation, not a
  clinical prediction.
```

---

## PART II — THE DATA LANDSCAPE

```
What data a breast cancer patient actually
has determines what geometric resolution
is achievable. The protocol adapts to the
data. The data does not determine whether
the protocol runs — it determines the
precision of the output.

Three data tiers exist in the real world.

─────────────────────────────────────────
TIER 1 — UNIVERSAL
What every diagnosed patient has.

  From the pathology report:
    ER status
    PR status
    HER2 status
    Ki-67 percentage
    Tumour grade (Nottingham 1-3)
    Histological type (IDC / ILC / other)
    Tumour size and nodal status
    (if surgical pathology available)

  What this tier enables:
    Geometric REGION placement.
    Not a point. A region.
    The region is bounded by what these
    proxies can constrain.
    The constraint is real and informative.
    The patient is currently receiving
    less information from this data
    than the geometry can provide.

  What this tier cannot provide:
    Direct gene expression values.
    Switch gene panel profiling.
    Exact depth score.
    Lock type confirmation.
    Claudin-low vs TNBC disambiguation.

─────────────────────────────────────────
TIER 2 — AVAILABLE TO SOME
Genomic risk score results.

  OncotypeDX (provides):
    Recurrence score (0-100)
    Individual gene group scores
    (ER group, PR group, HER2 group,
     proliferation group, invasion group)
    Some platforms include individual
    gene expression estimates.

  Prosigna / PAM50 (provides):
    Subtype assignment
    Risk of recurrence (ROR) score
    Centroid distance scores to each
    PAM50 subtype (LumA, LumB, HER2e,
    Basal — sometimes in detailed report)

  What this tier enables:
    Refined geometric placement.
    The centroid distances (if available)
    directly describe the patient's
    position relative to known population
    centroids in gene expression space.
    Proliferation group score = partial
    depth signal.
    PR group score = partial axis signal.

  What to request:
    "I would like the detailed gene group
    scores from my OncotypeDX report,
    not just the recurrence score."
    "I would like the centroid distance
    scores from my Prosigna report."
    Both are in the assay output.
    They are often not in the patient summary.
    They can be requested from the
    ordering physician.

─────────────────────────────────────────
TIER 3 — AVAILABLE TO A MINORITY
Full transcriptome expression data.

  Sources:
    Tempus xR (23,000+ gene expression)
    Caris MI Profile with WTS
    Foundation Medicine with RNA panel
    Research-grade RNA-seq (if enrolled
    in a research study)

  What this tier enables:
    Full geometric analysis.
    Every step of the protocol runs.
    FOXA1 and EZH2 directly measured.
    All switch genes measurable.
    Depth score computable.
    Lock type identifiable.
    Position on the population map precise.

  How to get it:
    The data exists if any of these tests
    were run. It is the patient's data.
    It is accessible under HIPAA.

    Request language (exact):
      "I am requesting the complete gene
      expression table from my [test name]
      analysis. I am requesting the
      normalised expression values for all
      genes as a CSV or text file under
      my HIPAA right of access to my
      medical records (45 CFR § 164.524).
      Please provide within 30 days."

    If the laboratory or provider refuses:
      Send the request in writing.
      Cite the regulation explicitly.
      Escalate to the laboratory director.
      The right is unambiguous.
      The data belongs to the patient.
```

---

## PART III — THE GEOMETRIC ANALYSIS
## PROCEDURE BY DATA TIER

```
The procedure below is written in the
order the analysis actually runs.
The steps apply at all tiers.
What differs is the precision of the
output at each step, not whether the
step runs.

Every step that can be executed with the
available data is executed.
Every step that cannot is documented as
a limitation, not skipped silently.
```

### Step 1 — Establish the Reference Frame

```
Before any patient data is examined:

  STATE which population reference
  geometry applies.

  For breast cancer:
    Population analysis: BRCA-S8b
    Dataset: GSE176078
    Normal reference populations:
      Mature Luminal
      Luminal Progenitors
      Myoepithelial
    Cancer populations mapped:
      LumA, LumB, HER2e, Basal/TNBC,
      Claudin-low (geometric classifier),
      ILC (TCGA histology annotation)
    Primary axis confirmed:
      FOXA1/EZH2 (result of saddle point
      scan — not assumed in advance)
    Depth score formula:
      Derived from BRCA-S8b
      (switch gene composite vs FA panel
       composite, normalised to mature
       luminal reference)

  This reference frame is the map.
  The patient's data will be placed on it.
  The reference frame does not change
  based on the patient's data.
  The patient's position on the map is
  what is discovered.

  WRITE THIS DOWN before looking at
  any patient data.
  The reference frame is fixed before
  the analysis begins.
  This is the before-document principle
  applied to the individual patient case.
```

### Step 2 — Examine the Data Without Conclusions

```
Receive the patient data.
Look at what is present.
Do not interpret yet.

For Tier 1 (IHC values):
  List every value present:
    ER: [value and units]
    PR: [value and units]
    HER2: [IHC score and/or ISH result]
    Ki-67: [percentage]
    Grade: [1/2/3]
    Histotype: [IDC/ILC/other]
  Note anything absent or ambiguous.
  Note units explicitly.
  (ER 80% is not the same as ER H-score 80)

For Tier 2 (Genomic score):
  List every component score present.
  Note which components are missing.
  Note the assay platform.

For Tier 3 (Expression table):
  Confirm format: TPM / RPKM / raw counts /
  normalised log2.
  Confirm gene identifier type:
    HGNC symbol / Ensembl ID / Entrez ID
  Confirm coverage: are the reference genes
  present in the data?
  Check the following genes are present
  and have non-zero values:
    FOXA1, GATA3, ESR1, PGR, SPDEF,
    EZH2, MKI67, CDH1, CLDN3, CLDN4
  If any are absent: note it as a
  coverage limitation.

Do not proceed to Step 3 until Step 2
is fully documented.
What exists in the data is the foundation.
What does not exist is equally important.
```

### Step 3 — Place on the Primary Axis

```
The primary axis of breast cancer geometry
is the FOXA1/EZH2 axis.
This was the result of the population
saddle point scan. It is the confirmed axis.

For each data tier, the axis placement
procedure is:

─────────────────────────────────────────
TIER 3 (direct measurement):

  Compute the ratio:
    FOXA1_value / EZH2_value
  (using whatever normalisation is present
   in the expression table)

  Place on the population reference scale:
    LumA median: 9.38
    LumB median: 8.10
    HER2e median: 3.34
    TNBC median: 0.52
    Claudin-low median: 0.10

  Note:
    The population medians are in RNA-space
    from GSE176078 / TCGA-BRCA.
    The patient's expression table must be
    in a comparable space for direct
    comparison. If the normalisation method
    differs, the ratio is directionally
    informative but not directly comparable
    to population medians.
    State this explicitly if it applies.

  State the patient's ratio and where it
  falls relative to the population scale.
  Do not round to the nearest subtype label.
  State the position precisely.
  "Patient ratio: 6.2. This places the
  patient between the LumB median (8.10)
  and the HER2e median (3.34), closer to
  the LumB zone."

─────────────────────────────────────────
TIER 2 (partial axis signal):

  OncotypeDX ER group score and PR group
  score provide partial axis information.

    ER group score = partial FOXA1-axis proxy
    (ESR1, BCL2, SCUBE2 are in the ER group)

    PR group score = partial EZH2-competition
    proxy (PGR is the EZH2-sensitive target)

    Proliferation group score = partial
    depth signal (Ki-67, STK15, Survivin,
    CCNB1, MYBL2 are proliferation markers)

  A high ER group score with a low PR group
  score is the signature of EZH2 beginning
  to compete at the FOXA1 axis — this is
  the LumB geometry signal.

  Prosigna centroid distances (if available)
  give direct geometric position in the
  PAM50 subtype space.

  State what the available scores imply
  about axis position.
  State the confidence level honestly.
  "Based on OncotypeDX group scores, the
  axis position is consistent with the
  [LumA / LumB / boundary] region. Direct
  ratio computation is not possible from
  this data tier."

─────────────────────────────────────────
TIER 1 (proxy inference):

  ER and PR are transcriptionally downstream
  of FOXA1. They are proxies, not direct
  measurements.

  The inference chain is:
    ER expression requires FOXA1 binding
    to the ESR1 promoter.
    PR expression requires FOXA1 + GATA3
    binding to the PGR promoter.
    PGR is more sensitive to EZH2 competition
    than ESR1.
    Therefore:
      ER high + PR high → FOXA1 dominant,
        EZH2 competition minimal
        → high ratio region (LumA zone)
      ER high + PR low + Ki-67 high →
        FOXA1 partially active, EZH2
        competition at PGR visible
        → intermediate ratio region (LumB zone)
      ER low or negative → FOXA1 suppressed
        → low ratio region (HER2e, TNBC,
          or CL zone)

  This inference is directional. It is not
  a ratio computation. It places the patient
  in a region, not at a point.

  State the region. State the basis for it.
  State what would sharpen it.

  "Based on ER [value] and PR [value] with
  Ki-67 [value]%, the axis position is
  consistent with the [region] zone. This
  is a proxy inference from downstream
  markers, not a direct measurement of
  FOXA1 or EZH2 expression. The LumA/LumB
  ambiguity [is / is not] resolvable from
  Tier 1 data alone."
```

### Step 4 — Compute or Estimate Attractor Depth

```
Attractor depth measures how far the cells
have moved from the normal reference state.
Deeper attractor = stronger lock = harder
to unlock = different drug priority.

This is NOT a prognosis measure.
It is a geometric position measure.
The clinical implications (drug priority,
lock type) flow from the geometry.
Prognosis is not the output.

─────────────────────────────────────────
TIER 3 (depth score computable):

  The depth score from the population
  analysis uses:

  SWITCH GENE PANEL (suppressed in cancer):
    FOXA1, GATA3, ESR1, CDH1,
    KRT18, KRT8, SPDEF

  FA MARKER PANEL (elevated in cancer):
    EZH2, MKI67, TOP2A, PCNA

  Computation:
    switch_score = mean of z-scores of
    switch genes relative to mature luminal
    reference (negative = suppressed)

    fa_score = mean of z-scores of FA
    markers relative to mature luminal
    reference (positive = elevated)

    depth_score = fa_score - switch_score
    (higher = deeper attractor)

  For an individual patient without a
  paired normal, use the population
  mature luminal mean and SD from BRCA-S8b
  as the reference.

  A patient-level depth score places the
  patient on the population depth
  distribution established in BRCA-S8b.

  Report the depth score and its percentile
  relative to the population distribution.
  Report which panel members contributed
  most to the score (which switch genes
  are most suppressed, which FA markers
  are most elevated).

─────────────────────────────────────────
TIER 2 (partial depth signal):

  OncotypeDX proliferation group score
  is a partial FA signal.
  Higher proliferation group = greater FA
  marker elevation = deeper attractor signal.

  OncotypeDX invasion group score (CTSL2,
  MMP11) is a partial mesenchymal depth
  signal.

  Prosigna ROR score is a partial depth
  proxy — higher ROR = deeper attractor.

  Combine available signals into a
  qualitative depth estimate:
    Shallow / Moderate / Deep / Very Deep

  State the basis and the confidence.

─────────────────────────────────────────
TIER 1 (depth indicators only):

  Grade 1 → Shallow attractor indicators
  Grade 2 → Moderate depth indicators
  Grade 3 → Deep attractor indicators

  Ki-67 < 10% → Shallow FA signal
  Ki-67 10-20% → Moderate FA signal
  Ki-67 > 20% → Deep FA signal
  Ki-67 > 30% → Very deep FA signal

  PR absent with ER present → EZH2
  competition visible at PGR promoter.
  This is a depth signal within the
  luminal zone specifically. It does not
  indicate overall depth — it indicates
  the EZH2 lock is active at a specific
  target.

  Combine into qualitative depth estimate.
  State basis. State limitations.
```

### Step 5 — Identify the Attractor Type

```
The four attractor types from
Attractor_Geometry_Axioms.md are:

  TYPE I — BLOCKED APPROACH
  Cells have the correct lineage identity
  but cannot complete differentiation.
  Identity TF is present but suppressed.
  The block is above the correct valley.

  TYPE II — WRONG VALLEY
  Cells are in an incorrect attractor.
  An oncogenic signal drives and maintains
  the wrong identity programme.
  Identity TF of the correct lineage is
  absent or very low.

  TYPE III — OVERSHOT IDENTITY
  Cells are in the correct valley but the
  valley floor has been removed.
  The terminal differentiation programme
  is present but uncontrolled (wrong copy
  number, wrong stoichiometry, wrong
  downstream wiring).

  TYPE IV — ROOT LOCK
  Cells are arrested at the pre-commitment
  node. They have not adopted any mature
  lineage identity. They are stem-like
  in a specific structural sense — not
  by surface marker but by transcriptional
  programme.

─────────────────────────────────────────
HOW TO CLASSIFY FROM DATA:

  The Identity TF Direction Test:
    What is the expression level of the
    lineage-defining transcription factor
    (FOXA1 for breast epithelial lineage)?

    HIGH → The lineage identity is present.
           The cancer has not escaped the
           luminal lineage.
           Type I (blocked) or Type III
           (overshot) depending on
           whether the valley floor is
           intact.

    LOW → The lineage identity is absent
          or suppressed.
          Type II (wrong valley — if an
          alternative oncogenic programme
          is dominant) or Type IV
          (root lock — if no mature lineage
          programme is visible at all).

  For breast cancer specifically:
    FOXA1 HIGH + EZH2 HIGH:
      Identity present but under active
      competitive suppression.
      TYPE I — Blocked Approach.

    FOXA1 LOW + HER2 (ERBB2) HIGH:
      Lineage identity overridden by
      HER2 amplification programme.
      TYPE II — Wrong Valley.

    FOXA1 LOW + EZH2 HIGH + claudin LOW
    + immune markers HIGH:
      Pre-commitment arrest.
      No mature lineage programme.
      TYPE IV — Root Lock (Claudin-low).

    FOXA1 PRESENT + CDH1 ABSENT:
      ILC geometry. Composite type.
      Luminal identity (FOXA1) present
      but structural programme (CDH1,
      E-cadherin) disrupted.

─────────────────────────────────────────
TIER 1 CLASSIFICATION:

  Use ER, HER2, and histotype as proxies
  for the Identity TF Direction Test.

  ER HIGH → FOXA1 programme active
  HER2 POSITIVE → HER2 programme dominant
  ER LOW + HER2 LOW → FOXA1 suppressed
  ILC histotype → CDH1 loss confirmed
                  (structural disruption)

  State the attractor type.
  State the basis.
  State the confidence.
  A Tier 1 attractor type classification
  is directionally reliable but cannot
  resolve:
    TYPE I depth (shallow vs deep)
    TYPE IV vs deep TYPE I in TNBC
    These require additional data.
```

### Step 6 — Identify the Epigenetic Lock

```
The lock is what maintains the attractor.
Identifying the lock is not identifying
the driver mutation. It is identifying
the epigenetic maintenance mechanism
that keeps the cells in the false attractor.
The lock is the drug target.

─────────────────────────────────────────
LOCK TYPES IN BRCA (from BRCA-S8b):

  LumA:
    EZH2 lock MINIMAL.
    CDK4/6 cell cycle lock prominent
    (CDKN1A suppressed).
    Primary drug target: CDK4/6 inhibitors.

  LumB:
    EZH2 lock ACTIVE at PGR and GATA3
    promoters. DNMT3A/HDAC2 co-elevation
    (r=+0.267 vs r=+0.071 in LumA).
    Primary drug target: HDACi + ET.
    Secondary: EZH2i where EZH2 elevation
    is the dominant signal.

  HER2-enriched:
    HER2 amplification drives false identity.
    EZH2 secondary lock in the deep fraction
    (CDH3-high, AR-low, EZH2 +118%).
    Primary drug target: anti-HER2.
    Secondary in deep fraction: EZH2i.

  TNBC:
    EZH2 dominant lock.
    FOXA1 near-completely suppressed.
    EZH2 +189% vs Mature Luminal.
    Primary drug target: tazemetostat.
    Sequence: EZH2i → fulvestrant
    (unlock FOXA1 programme → restore ET
    sensitivity → apply ET).

  Claudin-low:
    TYPE IV — Root Lock.
    Not an epigenetic lock on a specific
    lineage gene. An arrest at the
    pre-commitment node.
    Immune infiltration prominent.
    Drug target: immunotherapy (anti-TIGIT,
    anti-PD-1) to target the immune
    context of the uncommitted state.
    FOXP3/CD8A ratio is the specific
    immune predictor (HR=2.212 from
    BRCA-S8h CS-LIT-20).

  ILC:
    Structural lock. CDH1 (E-cadherin)
    lost. Luminal identity (FOXA1)
    preserved. The lock is structural
    (CDH1 loss changes downstream
    β-catenin signalling) rather than
    epigenetic in the primary axis sense.
    Drug implication: fulvestrant preferred
    over aromatase inhibitors specifically
    in ILC geometry (CS-LIT-18).

─────────────────────────────────────────
HOW TO IDENTIFY THE LOCK FROM PATIENT DATA:

  TIER 3:
    Measure EZH2 directly.
    Measure HDAC1, HDAC2, DNMT3A directly.
    Measure KDM6A (EZH2 eraser — if low,
    EZH2 is unopposed).
    The dominant elevated epigenetic
    regulator is the lock.
    Confirm against the attractor type
    from Step 5.

  TIER 2:
    OncotypeDX does not directly measure
    epigenetic regulators.
    Lock inference from subtype + depth:
    LumB pattern → EZH2/HDAC lock inferred.
    Requires Tier 3 for confirmation.

  TIER 1:
    Lock inference only.
    PR absent with ER present → EZH2
    competition at PGR inferred.
    Histotype ILC → CDH1 structural lock.
    TNBC → EZH2 dominant lock inferred.
    All Tier 1 lock identifications are
    inferences. State them as such.
    "The lock type is inferred from
    proxy data. Direct EZH2 measurement
    would confirm or revise this inference."
```

### Step 7 — Locate on the Drug Map

```
The drug map is a consequence of the
geometry, not an independent classification.
Steps 3 through 6 produce the geometric
description. The drug map follows from it.

This step does not produce treatment
recommendations. It produces geometric
observations about drug target relevance
given the identified position, depth,
type, and lock.

─────────────────────────────────────────
FOR EACH DRUG TARGET RELEVANT TO BRCA:

  CDK4/6 inhibitors (palbociclib,
  ribociclib, abemaciclib):
    Relevant when: LumA geometry,
    CDKN1A suppressed (depth signal),
    cell cycle lock dominant.
    Geometric basis: CDKN1A loss enables
    unregulated CDK4/6 activity.
    Observation: "CDK4/6i geometry is
    [present / absent / uncertain]
    based on [evidence from patient data]."

  HDAC inhibitors (entinostat):
    Relevant when: LumB geometry,
    DNMT3A/HDAC2 co-elevation present,
    PR-absent pattern.
    Geometric basis: HDAC lock at GATA3/PGR
    promoters.
    Observation: "HDACi + ET geometry is
    [present / absent / uncertain]."
    Note: TFF1/ESR1 decoupling is the
    specific patient selector for this
    target (CS-LIT-23). If available
    in Tier 3 data, compute and report.

  EZH2 inhibitors (tazemetostat):
    Relevant when: TNBC geometry with
    EZH2 dominant lock, or HER2-deep
    fraction (CDH3-high, AR-low, EZH2 high).
    Geometric basis: EZH2 maintaining
    FOXA1 suppression.
    Critical sequence: EZH2i must precede
    fulvestrant to first unlock the FOXA1
    programme, then exploit restored ET
    sensitivity (CS-LIT-16).
    Observation: "Tazemetostat → fulvestrant
    sequence geometry is [present / absent /
    uncertain]. EZH2 elevation [confirmed
    directly / inferred from TNBC pattern /
    not assessable from available data]."

  Immunotherapy (anti-TIGIT, anti-PD-1):
    Relevant when: Claudin-low geometry
    (TYPE IV).
    Patient selector: FOXP3/CD8A ratio
    (HR=2.212).
    SKYLINE trial (NCT06175390): tiragolumab
    + atezolizumab.
    Observation: "Immunotherapy geometry
    [present — claudin-low pattern confirmed
    / uncertain — claudin-low not resolvable
    from available data / absent]."

  Fulvestrant vs aromatase inhibitors:
    ILC-specific distinction.
    Relevant when: ILC histotype confirmed
    and ER-positive.
    Geometric basis: CDH1 loss creates
    specific ER signalling context where
    fulvestrant outperforms AI.
    Observation: "ILC-specific fulvestrant
    preference geometry [applies / does not
    apply / not assessable]."

─────────────────────────────────────────
FOR EACH DRUG TARGET:

  State: relevant / not relevant / uncertain
  State the geometric basis
  State what data would resolve uncertainty
  Do not use clinical recommendation language
  Use geometric observation language throughout
```

### Step 8 — Produce the Report

```
The report has a fixed structure.
Every section is always present.
Absent sections are not omitted —
they are present with a stated reason
for why they cannot be completed from
the available data.

─────────────────────────────────────────
REPORT STRUCTURE:

SECTION 1 — DATA RECEIVED
  List what was provided.
  List format, completeness, any quality notes.
  State data tier.

SECTION 2 — REFERENCE FRAME
  State which population analysis applies.
  State the reference populations.
  State the primary axis and where it came from.
  (This section is identical for every BRCA
  patient. The map does not change patient
  to patient.)

SECTION 3 — PRIMARY AXIS POSITION
  State the result of Step 3.
  Position on the FOXA1/EZH2 axis.
  Where this falls relative to population
  medians.
  What region of the attractor landscape
  this corresponds to.
  Confidence level (Tier 1: region only /
  Tier 2: refined region / Tier 3: precise).

SECTION 4 — ATTRACTOR DEPTH
  State the result of Step 4.
  Depth score (Tier 3) or depth estimate
  (Tier 1/2).
  What population percentile (if computable).
  What the depth implies geometrically.

SECTION 5 — ATTRACTOR TYPE
  State the attractor type from Step 5.
  State the basis (which data supports it).
  State confidence.
  If ILC: note composite classification.
  If TNBC/CL ambiguity: state explicitly.

SECTION 6 — LOCK IDENTIFICATION
  State the lock type from Step 6.
  Direct measurement or inference.
  What the lock implies for drug target
  geometry.

SECTION 7 — DRUG TARGET GEOMETRY
  State each drug target observation
  from Step 7.
  Relevant / Not relevant / Uncertain.
  Geometric basis for each.
  What data would resolve uncertain cases.

SECTION 8 — THE KEY QUESTION FOR THE
CLINICAL TEAM
  One question only.
  Derived from the most clinically
  significant geometric finding.
  Written in language the patient can
  bring to their oncologist.
  Uses the patient's own clinical data
  as the basis.
  Does not include geometric framework
  language unless the patient wants it.

  Format:
  "Based on [patient's own values], the
  geometric analysis suggests [finding].
  The question for your clinical team is:
  [specific clinical question]."

SECTION 9 — HONEST UNCERTAINTY
  What is established from available data.
  What cannot be established.
  What additional data would add.
  What this analysis cannot determine
  (clinical factors, drug interactions,
  comorbidities, prior treatment history).

SECTION 10 — WHAT THE GEOMETRY SHOWS
  Two to three sentences of plain language.
  Not biology jargon.
  Not framework language.
  What the geometry found, stated for
  a person who is not a scientist.
  This is the last section.
  It is the one the patient reads first.

─────────────────────────────────────────
TONE RULES FOR THE REPORT:

  State positions, not conclusions.
  "The geometry is consistent with X"
  not "You have X."

  State observations, not predictions.
  "The drug target geometry is present"
  not "This drug will work for you."

  State questions, not recommendations.
  "The question for your team is whether..."
  not "You should ask about..."

  State uncertainty where it exists.
  Never imply precision you do not have.
  Never imply certainty about clinical
  outcomes from geometric observations.

  State the founding principle if the
  patient asks why this is offered:
    "I do not want to take people's money
    and promise them bullshit."
```

---

## PART IV — THE UNIVERSAL PRINCIPLE
## THAT GOVERNS ALL CANCER TYPES

```
Everything in Parts I through III is
specific to BRCA.

The method is not specific to BRCA.

The method is:
  Establish the reference frame
  (from the population analysis for
  that cancer type).

  Place the patient's data in the
  reference frame.

  Let the geometry reveal the position,
  depth, type, and lock.

  Report what the geometry shows.

For any of the 19 other cancer types
with completed population analyses in
this repository, the same steps apply.
Step 1 names a different population
analysis. Steps 3-6 use the genes and
reference values from that analysis.
Step 7 uses the drug predictions from
that analysis.

The method does not change.
The reference frame changes.
The genes change.
The geometry reveals itself in every case.

FOXA1 and EZH2 are the answer the geometry
gave for breast cancer.
For AML, the geometry gave different answers.
For GBM, different again.
The method was the same.
The method is what this protocol teaches.
The answers are what the population analyses
found.
An individual patient analysis applies
the answers to one person at a time.
It does not re-derive the answers.
It does not assume them either.
It places one person on the map that
the answers define, and reads where
they land.
```

---

## PART V — WHEN TO STOP

```
Stop and document without producing a
geometric report if:

  The patient's cancer type is not BRCA
  and the patient presents as a breast
  cancer patient but the data suggests
  a different primary. Refer to the
  correct cancer type protocol.

  The available data is insufficient to
  place the patient in any region of the
  landscape with confidence. This means:
  no ER status, no HER2 status, no grade.
  Without these, the Tier 1 placement
  cannot be done. Ask for the complete
  pathology report.

  The patient is seeking a second opinion
  on a clinical decision that has already
  been made and implemented. The geometric
  analysis of a post-treatment sample
  describes the post-treatment geometry,
  not the pre-treatment geometry. State
  this distinction clearly.

  The data is ambiguous in a way that
  cannot be resolved and where the
  ambiguity is clinically significant.
  Do not produce a report that could
  be misread as confident when the
  underlying data is genuinely ambiguous.
  Produce a limited report that states
  what is clear, what is ambiguous,
  and what would resolve it.
```

---

## DOCUMENT METADATA

```
document_id:    BRCA_INDIVIDUAL_GEOMETRIC_PROTOCOL
folder:         Individual_Protocol/Breast_Cancer/
type:           Cancer-type individual patient
                protocol
cancer_type:    BRCA (Breast Cancer)
version:        1.0
date:           2026-03-07
status:         ACTIVE

parent_document:
  Individual_Patient_Geometric_Analysis_Protocol.md

population_analysis:
  BRCA-S8b (BRCA_Cross_Subtype_Script1.py)
  BRCA-S8h (BRCA_Cross_Subtype_Literature_Check.md)

axiom_document:
  Cancer_Research/Attractor_Geometry_Axioms.md

governing_principle:
  "I do not want to take people's money
   and promise them bullshit."

note:
  FOXA1 and EZH2 appear in this document
  because the BRCA population analysis found
  them as the result of the saddle point scan.
  They are the answer the geometry gave.
  They are not imposed. They are named here
  because this document applies the result
  to individual patients.
  The method that found them — the saddle
  point scan, the top mover scan, the identity
  TF direction test — is what is universal.
  Not the specific genes.
```
