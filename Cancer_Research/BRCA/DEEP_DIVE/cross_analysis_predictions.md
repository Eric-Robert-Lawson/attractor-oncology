# BREAST CANCER — CROSS-SUBTYPE ANALYSIS
## Before-Document: Predictions Locked Before Any Script Runs
## OrganismCore — Document BRCA-S8a
## Date: 2026-03-05

---

## DOCUMENT METADATA

```
document_id:        BRCA-S8a
series:             BRCA Deep Dive — Cross-Subtype Analysis
folder:             Cancer_Research/BRCA/DEEP_DIVE/Cross_Subtype/
type:               BEFORE-DOCUMENT
                    All predictions locked before
                    any cross-subtype script runs.
                    No cross-subtype data examined
                    prior to this document.
date:               2026-03-05
author:             Eric Robert Lawson
                    OrganismCore
protocol:           Workflow_Protocol.md v2.0
ethics:             Individual_Protocol/ETHICS.md
prerequisite_docs:
    BRCA_Subtypes.md          (orientation)
    BRCA-S2c  Luminal A       literature check
    BRCA-S5c-LC Luminal B     literature check
    BRCA-S3e  HER2-enriched   literature check
    BRCA-S4e/f TNBC           literature checks ×2
    BRCA-S6e  ILC             literature check
    BRCA-S7e/i Claudin-low    literature checks ×2
    Attractor_Geometry_Axioms.md
    (Document 90 — TYPE 1 through TYPE 4 taxonomy)
status:             PREDICTIONS LOCKED
                    This document is complete when
                    the last prediction is stated.
                    It is not updated after script runs.
                    Results go in BRCA-S8b.
governing_principle:
    "I do not want to take people's money
     and promise them bullshit."
     — Eric Robert Lawson, March 4, 2026
```

---

## PREAMBLE — WHAT THIS ANALYSIS IS

```
Six breast cancer subtypes have now been
individually analysed:

    Luminal A    (BRCA-S1a through S2c)
    Luminal B    (BRCA-S5a through S5c-LC)
    HER2-enriched (BRCA-S3a through S3e)
    TNBC/Basal   (BRCA-S4a through S4f)
    ILC          (BRCA-S6a through S6e)
    Claudin-low  (BRCA-S7a through S7l)

Each was treated as a separate disease.
Each has its own cell of origin.
Each has its own attractor type.
Each has its own drug map.

That is correct.
Breast cancer is not one disease.

AND —

The six analyses now share a common reference
coordinate system. Each was measured against the
same normal cell populations (Mature Luminal,
Luminal Progenitors from GSE176078). Each used the
same framework. Each produced a depth score, an
attractor classification, a switch gene profile,
and a drug map.

The cross-subtype analysis does one thing:
places all six subtypes on a single geometric map
and asks what that map reveals that the individual
analyses could not.

This document contains the predictions for
what that map will show.

Nothing from the cross-subtype data has been
examined. Every prediction below is derived from:
    1. The confirmed individual subtype findings
    2. The attractor geometry taxonomy
    3. The developmental biology of the
       mammary ductal hierarchy
    4. The confirmed cross-subtype geometric
       ordering already observed in Script 1
       results (HER2 log: FOXA1/GATA3/EZH2
       ordering across LumA/HER2/TNBC)

The predictions are specific.
They are stated before data.
They are falsifiable.
Wrong predictions will be documented
alongside confirmed predictions.
```

---

## PART I — THE FRAMEWORK INPUT
### What Each Individual Analysis Established

---

```
This is the complete input to the
cross-subtype predictions.
These findings are confirmed.
They are not predictions.
They are the foundation.

─────────────────────────────────────────
LUMINAL A (BRCA-S1/S2):
─────────────────────────────────────────

  Attractor type:   TYPE 1-L (Locked Luminal)
  Cell of origin:   Mature luminal epithelial cell
  Primary lock:     CDK4/CDK6 — cell cycle arrest
                    axis dismantled
  Switch genes:     FOXA1 (+37% vs mature luminal —
                    ELEVATED not suppressed in LumA SC)
                    GATA3 (+34% — elevated)
                    ESR1  (flat — present but
                           uncoupled from outputs)
  Key novel finding: LumA is NOT a cell that has lost
                    luminal identity. It has AMPLIFIED
                    it — and lost the exit mechanism.
                    FOXA1/GATA3 are present; the
                    CDK4/6 brake is missing.
                    TGFB→SMAD→CDKN1A circuit broken.
                    NCOA1/NCOA2 coactivator gap between
                    ESR1 and its downstream targets.
  Depth axis:       CDKN1A (p21) depletion is the
                    primary depth variable.
                    TFF1/TFF3 as output markers.
  Drug target:      CDK4/6i + ET (confirmed convergent)
                    NCOA axis (partially novel)
  EZH2:             +13% vs mature luminal — minimal

─────────────────────────────────────────
LUMINAL B (BRCA-S5):
─────────────────────────────────────────

  Attractor type:   TYPE 1-L (Deeper Locked Luminal)
  Cell of origin:   Luminal progenitor
                    (less committed than LumA)
  Primary lock:     HDAC1/2 + DNMT3A co-complex
                    Chromatin-level ER output blockade
  Switch genes:     ESR1 +64.4% vs mature luminal
                    (HIGHER than normal)
                    TFF1 −82.9% vs LumA despite intact
                    ESR1 → chromatin-level silencing
                    GATA3 expressed but lower than LumA
  Key novel finding: DNMT3A-HDAC2 coupling r=+0.267
                    in LumB vs r=+0.071 in LumA
                    (3.8× stronger — confirmed p=5.68e-56)
                    LumB is not a less-differentiated LumA.
                    It retains and AMPLIFIES luminal
                    identity while LOSING ER output.
                    The receptor is present; the output
                    is blocked at the chromatin level.
  Depth axis:       CDKN1A depletion (same as LumA but
                    deeper — 21-fold range across bulk)
                    HDAC1/2 elevation as secondary
  Drug target:      CDK4/6i + ET + HDACi (entinostat)
                    Novel specificity: LumB-specific
                    HDAC1/2 co-elevation not in prior
                    trial selection
  EZH2:             Intermediate elevation (higher than
                    LumA, lower than TNBC)

────���────────────────────────────────────
HER2-ENRICHED (BRCA-S3):
─────────────────────────────────────────

  Attractor type:   TYPE 1-A (Amplicon Arrest)
                    Novel type — slope arrest by
                    copy number event rather than
                    epigenetic lock or identity loss
  Cell of origin:   Luminal progenitor
                    (ER-negative — not terminally
                    differentiated)
  Primary lock:     ERBB2 amplicon (17q12)
                    Not a TF suppression.
                    Not a chromatin lock.
                    A constitutive proliferative signal
                    from copy number amplification.
  Switch genes:     FOXA1 retained (+12.7% vs normal —
                    elevated but less than LumA)
                    GATA3 −44.7% vs normal (suppressed)
                    ESR1 −38.0% vs LumA in bulk
  Key novel finding: FOXA1 is retained in HER2-enriched
                    (luminal scaffold intact).
                    EZH2 +176% (intermediate —
                    between LumA +13% and TNBC +224%).
                    PCA distance from Mature Luminal:
                    0.81 (HER2) vs 0.92 (LumA) vs 3.50 (TNBC)
                    HER2 is GEOMETRICALLY CLOSER to
                    mature luminal than LumA SC.
                    (PCA paradox — confirmed in log)
                    CDH3 as dominant novel gained gene.
                    AR-low subpopulation (24.1%) is the
                    pre-resistant deep end.
  Depth axis:       ERBB3 / CDH1 / AR composite
                    (inverse — lower = deeper = worse OS)
  Drug target:      Trastuzumab primary.
                    EZH2i for deep subpopulation.
                    CDH3 ADC (BC3195) for pre-resistant.
  EZH2:             +176% — intermediate

─────────────────────────────────────────
TNBC / BASAL-LIKE (BRCA-S4):
─────────────────────────────────────────

  Attractor type:   TYPE 2 (Wrong Valley)
                    Full luminal identity loss.
                    Basal/mesenchymal programme engaged.
                    Cell is in a different valley,
                    not merely blocked on the approach
                    to the luminal valley.
  Cell of origin:   Luminal progenitor OR
                    myoepithelial progenitor
                    (BRCA1-mutant → composite type 1→2)
  Primary lock:     EZH2 / PRC2 complex
                    H3K27me3 silencing of FOXA1/GATA3
  Switch genes:     FOXA1 −80.7% (p=8.34e-162)
                    GATA3 −53.4% (p=2.30e-104)
                    ESR1  −96.7% (p=0.00)
                    SOX10 +1323% (neural crest programme —
                    cross-lineage marker of false attractor)
  Key novel finding: AR as continuous depth biomarker
                    (r=−0.547). EED:EZH2 ratio as
                    tazemetostat patient selector.
                    pCR/DRFS paradox explained by
                    single depth score variable.
                    Ferroptosis as novel depth-targeted
                    therapy direction (ZEB1→GPX4 axis).
                    Composite type: BRCA1 loss (Type 1)
                    precipitates Type 2 fall into
                    basal attractor.
  Depth axis:       AR (continuous, r=−0.547)
                    ZEB1/ZEB2 (depth correlates)
                    SOX10 (false attractor membership)
  Drug target:      Tazemetostat → Fulvestrant sequence
                    PARPi + EZH2i combination
                    Ferroptosis (GPX4 inhibition) novel
  EZH2:             +224% to +270% — highest among
                    classical subtypes

─────────────────────────────────────────
ILC (BRCA-S6):
─────────────────────────────────────────

  Attractor type:   TYPE 1-C (Cohesion Loss)
                    Unique geometry — not a depth
                    arrest along the luminal axis.
                    A structural loss (CDH1/E-cadherin)
                    while luminal identity is retained
                    and in fact HYPERACTIVATED.
  Cell of origin:   Lobular epithelial cell
                    (different normal starting point
                    from ductal subtypes)
  Primary lock:     CDH1 loss (mutation ~65%,
                    methylation ~35%)
                    EZH2 as depth modifier within ILC
  Switch genes:     FOXA1 ABOVE normal (hyperactivated)
                    GATA3 ABOVE normal (hyperactivated)
                    ESR1 present and active
                    CDH1 ABSENT (the defining event)
  Key novel finding: Luminal TFs HYPERACTIVATED above
                    normal baseline in ILC — opposite
                    of all other subtypes.
                    ILC and TNBC are GEOMETRIC OPPOSITES:
                    ILC = maximum luminal identity +
                    CDH1 loss
                    TNBC = minimum luminal identity +
                    SOX10 cross-lineage programme
                    MKI67 as dominant within-ILC
                    survival axis (HR=3.2 in METABRIC)
                    EZH2-high ILC: HR=2.656 (novel OS signal)
                    Composite escape ILC = pleomorphic ILC
  Depth axis:       CDH1 loss depth (orthogonal to ESR1)
                    MKI67 proliferative escape
                    EZH2 composite escape
  Drug target:      ET + CDK4/6i (primary)
                    EZH2i dual mechanism (novel ILC-specific)
                    Anti-HER2 has NO geometric basis for ILC
  EZH2:             Elevated — ILC-specific mechanism

─────────────────────────────────────────
CLAUDIN-LOW (BRCA-S7):
─────────────────────────────────────────

  Attractor type:   TYPE 4 (Root Lock)
                    First confirmed TYPE 4 in repository.
                    Below commitment threshold.
                    Both luminal AND basal identity
                    simultaneously absent.
  Cell of origin:   Mammary stem cell (memory-low)
                    OR de-differentiated luminal/basal
                    (memory-high — via EMT)
  Primary lock:     Pre-commitment programme
                    (EZH2 elevated, BRD4 predicted)
                    No single dominant epigenetic lock —
                    the lock IS the pre-commitment state
  Switch genes:     FOXA1 ~−100% (bilateral absence)
                    GATA3 ~−100% (bilateral absence)
                    ESR1  ~−100% (bilateral absence)
                    TP63/KRT5 also absent (not merely
                    luminal loss — basal identity also
                    absent)
                    VIM/FN1/SNAI1/ZEB1/CD44 elevated
  Key novel finding: TYPE 4 Root Lock taxonomy itself —
                    novel structural classification.
                    Memory-low / memory-high two-axis
                    structure (independent of depth axis).
                    CT antigen de-repression tracks
                    memory-low status (p=8.86e-09 TCGA,
                    p=0.0004 METABRIC — confirmed ×2).
                    OS null in untreated patients =
                    immune cancellation mechanism.
                    Memory-low = greatest predicted
                    anti-TIGIT benefit.
                    FOXP3/CD8A ratio clinical-grade IHC
                    biomarker already validated.
                    Anti-PD-1 monotherapy contraindicated.
  Depth axis:       Depth score (HR≈2-3 in ESR1-low CL
                    confirmed ×3 cohorts)
                    Memory axis independent of depth
  Drug target:      Anti-TIGIT → anti-PD-1 sequence
                    (memory-low, FOXP3-high only)
                    CT antigen targeting (GAGE) step 3
                    CD44/TGFβ + BETi (depth-high)
  EZH2:             Elevated — consistent with TNBC
                    direction but not dominant node
```

---

## PART II — THE CROSS-SUBTYPE PREDICTIONS
### Locked Before Any Script Runs

---

### AXIS 1 PREDICTIONS — THE DEPTH SPECTRUM

---

```
PREDICTION CS-1 (LOCKED 2026-03-05):

THE BREAST CANCER DEPTH AXIS IS
A CONTINUOUS SPECTRUM, NOT DISCRETE CLASSES.

When all six subtypes are placed on a single
geometric axis measuring luminal identity
retention (FOXA1 + GATA3 + ESR1 composite),
they will form a continuous ordered spectrum:

MOST DIFFERENTIATED → LEAST DIFFERENTIATED:

  ILC > LumA > LumB > HER2 > TNBC > Claudin-low

WITH THE FOLLOWING SPECIFIC ORDERING RULES:

  R1: ILC will have the HIGHEST luminal TF
      composite of any cancer subtype —
      ABOVE the Mature Luminal normal.
      (Confirmed individually for FOXA1, GATA3
      in ILC script. Cross-subtype positioning
      not yet confirmed.)

  R2: LumA will have ELEVATED FOXA1/GATA3
      (relative to Mature Luminal) — confirmed.

  R3: LumB will have lower FOXA1/GATA3 than LumA
      but still positive relative to a deeper
      reference (Luminal Progenitor).

  R4: HER2-enriched will have FOXA1 retained
      (confirmed: +12.7%) but GATA3 suppressed
      (confirmed: −44.7%) — a SPLIT TF profile
      unique to this subtype.

  R5: TNBC will show strong suppression of all
      three luminal TFs simultaneously (confirmed).

  R6: Claudin-low will show the deepest
      suppression of all three simultaneously
      (confirmed — ~100% for all three TFs).

FALSIFICATION:
  If LumB luminal TF profile is NOT between
  LumA and HER2 when expressed on the same
  scale and reference — prediction fails.

  If ILC does not exceed normal Mature Luminal
  on FOXA1/GATA3 composite — prediction fails.

  If HER2 FOXA1 retention is not distinguishable
  from TNBC FOXA1 suppression on the unified axis
  — prediction fails.

─────────────────────────────────────────

PREDICTION CS-2 (LOCKED 2026-03-05):

EZH2 ELEVATION TRACKS THE DEPTH SPECTRUM
MONOTONICALLY ACROSS ALL SUBTYPES.

EZH2 % elevation above Mature Luminal will
follow the same ordering as the luminal TF
suppression axis:

  Predicted EZH2 ordering (ascending):
    LumA:          +13%     (confirmed)
    ILC:           intermediate-to-low
                   (higher than LumA, lower than LumB)
    LumB:          intermediate (higher than LumA)
    HER2:          +176%    (confirmed)
    TNBC:          +224%    (confirmed)
    Claudin-low:   highest  (deepest attractor)

  SPECIFIC QUANTITATIVE PREDICTION:
    EZH2 in claudin-low will exceed TNBC.
    Both claudin-low and TNBC will exceed HER2.
    HER2 will exceed LumB.
    LumB will exceed LumA.
    LumA will have the lowest EZH2 elevation of
    all ER-negative or partially ER-negative subtypes.

  ILC is the exception to test:
    ILC has hyperactivated luminal TFs but still
    shows EZH2 elevation (confirmed in individual
    analysis). This tests whether EZH2 tracks
    luminal TF suppression (in which case ILC
    should have low EZH2) or whether EZH2 tracks
    something else in ILC (CDH1-depth independently
    of TF axis). PREDICTION: ILC EZH2 elevation
    will be LOWER than LumB on the unified scale,
    because the ILC false attractor is defined by
    CDH1 loss not by TF suppression — EZH2 in ILC
    is a secondary depth modifier, not the primary
    lock as in TNBC.

FALSIFICATION:
  If EZH2 in LumA exceeds EZH2 in HER2 — fails.
  If EZH2 in claudin-low is not the maximum —
  partially fails (directional ordering remains
  but magnitude prediction fails).
  If ILC EZH2 exceeds TNBC — fails (would suggest
  ILC has a deeper EZH2-driven lock than the wrong-
  valley subtypes, which contradicts the geometry).

─────────────────────────────────────────

PREDICTION CS-3 (LOCKED 2026-03-05):

THE PCA DISTANCE FROM MATURE LUMINAL IS
MONOTONICALLY ORDERED WITH CLINICAL PROGNOSIS.

Using the same PCA coordinate system established
in HER2 Script 1 (Mature Luminal as origin):

  Predicted PCA distance ordering:
    ILC:          < LumA  (hyperactivated luminal TFs
                           — closest to or above normal)
    LumA:         ~0.92   (confirmed)
    LumB:         0.92 < LumB < HER2  (between LumA
                           and HER2 on this axis)
    HER2:         0.81    (confirmed — closer than LumA
                           due to PCA paradox)
    TNBC:         3.50    (confirmed — maximally distant)
    Claudin-low:  > TNBC  (below commitment — should
                           be furthest of all)

  NOTE ON HER2 PCA PARADOX:
    The HER2 script confirmed that HER2 is
    GEOMETRICALLY CLOSER to Mature Luminal than
    LumA SC in PCA space (0.81 vs 0.92). This
    is because LumA SC has elevated FOXA1/GATA3
    (moving away from the Mature Luminal reference
    in the TF-elevation direction) while HER2
    retains FOXA1 partially but does not elevate
    it. The cross-subtype analysis must account for
    this. ILC is expected to be even closer to
    or above Mature Luminal on PC1 (due to
    hyperactivated TFs).

  CLINICAL CORRELATION PREDICTION:
    PCA distance from Mature Luminal will
    POSITIVELY correlate with clinical aggressiveness
    and NEGATIVELY correlate with prognosis
    across the six subtypes.
    5-year OS ordering (worst to best):
    Claudin-low ≈ TNBC > HER2 (untreated) > LumB > LumA > ILC
    (ILC has late recurrence but excellent early OS)

FALSIFICATION:
  If LumA is more distant from Mature Luminal
  than TNBC in PCA — fails completely.
  If Claudin-low PCA distance does not exceed
  TNBC — partially fails.
```

---

### AXIS 2 PREDICTIONS — THE ATTRACTOR TYPE MAP

---

```
PREDICTION CS-4 (LOCKED 2026-03-05):

THE SIX SUBTYPES OCCUPY FOUR DISTINCT
ATTRACTOR TYPES WITH ONE NOVEL INTERNAL DIVISION.

Classification predicted:

  TYPE 1 (Blocked Approach variants):
    TYPE 1-L Shallow (LumA):
      Differentiation TFs ELEVATED above normal.
      Exit mechanism (CDK4/6 brake) MISSING.
      Cell is above the valley floor, cannot exit.
      The most differentiated false attractor.

    TYPE 1-L Deep (LumB):
      Differentiation TFs present but output BLOCKED
      at chromatin level (HDAC1/2 + DNMT3A).
      Cell is in the valley, cannot read its own
      identity programme out into function.
      Deeper than LumA — the TF is present but
      the downstream is silenced.

    TYPE 1-A (HER2-enriched):
      Luminal scaffold (FOXA1) partially retained.
      Amplicon (ERBB2) drives constitutive
      proliferative signal that overrides
      the differentiation checkpoint.
      Copy number arrest rather than TF or
      chromatin arrest.
      Unique mechanism within TYPE 1 family.

    TYPE 1-C (ILC):
      Luminal TFs HYPERACTIVATED above normal.
      Cell is over-committed to luminal identity.
      The structural programme (CDH1-mediated
      cohesion) has failed — not the TF programme.
      The geometry is the inverse of TYPE 2/4.
      Unique mechanism: cohesion loss, not
      identity loss.

  TYPE 2 (Wrong Valley):
    TNBC/Basal:
      Luminal identity lost.
      Basal/mesenchymal programme active.
      EZH2 PRC2 lock maintaining wrong-valley state.
      Classic TYPE 2 geometry.

  TYPE 4 (Root Lock):
    Claudin-low:
      Below commitment threshold.
      Both luminal AND basal identity absent.
      Pre-commitment programme active.
      Memory-low/high internal axis.

  TYPE 3 (Correct Valley, Floor Removed):
    NOT OBSERVED in breast cancer.
    PREDICTION: No breast cancer subtype will
    show TYPE 3 geometry in this analysis.
    TYPE 3 requires ELEVATED differentiation TFs
    with a DISMANTLED arrest axis — e.g., LumA
    if CDK4/6 were the ONLY change and TFs were
    fully normal. In LumA, TFs are elevated
    above normal — this is closer to TYPE 1-L
    than TYPE 3. The distinction holds.

FALSIFICATION:
  If claudin-low shows any residual basal
  identity programme (TP63 or KRT5 elevation)
  — TYPE 4 classification would need revision.
  If TNBC shows ANY luminal TF elevation rather
  than suppression — TYPE 2 classification fails.
  If ILC shows luminal TF suppression rather than
  hyperactivation — TYPE 1-C fails entirely.

─────────────────────────────────────────

PREDICTION CS-5 (LOCKED 2026-03-05):

THE LOCK MECHANISM CHANGES SYSTEMATICALLY
ACROSS THE DEPTH SPECTRUM.

Predicted lock mechanism by depth:

  SHALLOW END (LumA):
    Lock type:    Kinase (CDK4/CDK6)
    Mechanism:    Cell cycle arrest axis.
                  Not epigenetic — enzymatic.
    Drug class:   CDK4/6 inhibitors (approved)
    EZH2 role:    Minimal — not the primary lock

  SHALLOW-INTERMEDIATE (LumB):
    Lock type:    Chromatin (HDAC1/2 + DNMT3A)
    Mechanism:    Epigenetic silencing of ER outputs
                  despite intact receptor.
    Drug class:   HDACi + CDK4/6i combination
    EZH2 role:    Rising but not dominant

  AMPLICON (HER2):
    Lock type:    Copy number (ERBB2 amplicon)
    Mechanism:    Constitutive kinase signalling
                  overriding differentiation checkpoint.
    Drug class:   Anti-HER2 (trastuzumab — approved)
    EZH2 role:    Intermediate — rising with depth
                  within subtype (deep HER2 fraction)

  DEEP WRONG VALLEY (TNBC):
    Lock type:    Epigenetic (EZH2/PRC2 dominant)
    Mechanism:    H3K27me3 silencing of FOXA1/GATA3/ESR1
                  loci. Full luminal identity erasure.
    Drug class:   EZH2 inhibitor (tazemetostat)
                  → conversion → fulvestrant
    EZH2 role:    DOMINANT CONVERGENCE NODE

  ROOT LOCK (Claudin-low):
    Lock type:    Pre-commitment programme
    Mechanism:    Cells never committed to somatic
                  identity — or retreated below
                  commitment threshold. No single
                  dominant lock — the pre-commitment
                  STATE is the lock.
    Drug class:   Treg depletion (anti-TIGIT) to
                  unmask immune recognition of
                  exposed germline antigens.
                  Commitment forcing — novel.
    EZH2 role:    Present but not dominant node
                  in the way it is in TNBC

  COHESION LOSS (ILC):
    Lock type:    Structural (CDH1 loss)
    Mechanism:    E-cadherin protein absent.
                  Tight junction architecture lost.
                  Cell identity programme intact —
                  structural programme destroyed.
    Drug class:   ET + CDK4/6i (luminal target valid)
                  EZH2i for deep escape subtype
    EZH2 role:    Depth modifier, not primary lock

FALSIFICATION:
  If CDK4/6 is the dominant target in TNBC
  — the lock mechanism progression fails.
  If EZH2 is not a meaningful drug target for
  LumB (i.e., if HDAC1/2 finding is spurious)
  — the chromatin lock prediction fails.
  If anti-HER2 has no geometric basis (i.e., if
  ERBB2 amplicon is not linked to proliferative
  arrest in the individual analysis) — fails.
```

---

### AXIS 3 PREDICTIONS — THE CROSS-SUBTYPE DRUG MAP

---

```
PREDICTION CS-6 (LOCKED 2026-03-05):

EZH2 INHIBITION IS A CROSS-SUBTYPE TARGET
BUT WITH DEPTH-STRATIFIED PRIORITY.

EZH2 is elevated in EVERY breast cancer subtype
in this analysis. But its role differs by depth:

  In LumA:    EZH2 elevation minimal (+13%).
              EZH2 inhibition alone would be
              insufficient as monotherapy.
              Not a primary target.

  In LumB:    EZH2 elevation intermediate.
              EZH2i was WITHDRAWN as a standalone
              prediction in the LumB literature check
              (CONVERGENT with emerging literature).
              May have role in combination with HDACi.
              Not first-line.

  In HER2:    EZH2 +176%. Rises within subtype.
              EZH2i + trastuzumab predicted for the
              deep HER2 fraction (AR-low, CDH1-low,
              ERBB3-low deep-end cells).
              Specific patient selection criterion:
              AR-low AND ERBB3-low AND CDH1-low.
              Not for all HER2-enriched.

  In TNBC:    EZH2 +224%. Dominant convergence node.
              EZH2i → fulvestrant sequence is the
              primary novel drug prediction for TNBC.
              EED:EZH2 ratio as patient selector
              (partially novel — literature emerging).
              PARPi + EZH2i combination for
              BRCA1-mutant composite type.
              EZH2i has highest priority here.

  In ILC:     EZH2 elevated as composite escape marker.
              EZH2-high ILC: HR=2.656.
              EZH2i dual mechanism in ILC is novel:
              (1) reduces proliferative escape
              (2) may partially restore CDH1
                  (H3K27me3 at CDH1 promoter in
                  methylation-silenced ILC cases)
              Patient selection: EZH2-high AND
              MKI67-high ILC.

  In CL:      EZH2 present but not dominant.
              BRD4 predicted as the more important
              epigenetic target (BETi predicted
              — Script 5 Priority 5, not yet tested).

CROSS-SUBTYPE RULE DERIVED:
  EZH2 inhibition priority follows the depth axis:
  TNBC > HER2 (deep fraction) > ILC (deep fraction)
  > LumB (combination only) > LumA (not primary)
  > Claudin-low (BRD4 more important)

FALSIFICATION:
  If EZH2 elevation in LumA equals or exceeds
  TNBC — the depth-stratified priority rule fails.
  If EZH2i shows no predicted utility in HER2
  deep fraction — partially fails.

─────────────────────────────────────────

PREDICTION CS-7 (LOCKED 2026-03-05):

THE LUMINAL IDENTITY AXIS (FOXA1/GATA3/ESR1)
IS THE MASTER VARIABLE OF THE BREAST CANCER
TREATMENT LANDSCAPE.

Specifically:

  ESR1 status determines endocrine therapy
  eligibility. Confirmed standard of care.
  (Convergent — not novel.)

  FOXA1 status does MORE than determine ER
  binding. FOXA1 level distinguishes:
    ILC/LumA (FOXA1 high or hyper) →
    CDK4/6i sensitive, ET primary target
    HER2 (FOXA1 partially retained) →
    trastuzumab primary, EZH2i for deep fraction
    TNBC (FOXA1 suppressed) →
    EZH2i to restore FOXA1 → then ET
    Claudin-low (FOXA1 absent) →
    commitment forcing required before any
    luminal-directed therapy can engage

  FOXA1 is a THERAPEUTIC READINESS MARKER:
    FOXA1 present → endocrine therapy can engage
    FOXA1 partially present → can be restored
    by EZH2i (TNBC prediction: tazemetostat
    restores FOXA1/GATA3)
    FOXA1 absent (claudin-low) → requires
    commitment forcing before restoration

  NOVEL CROSS-SUBTYPE PREDICTION:
    FOXA1 IHC score will be MORE PREDICTIVE
    of EZH2 inhibitor response across ALL
    TNBC/claudin-low patients than EZH2
    expression level alone.
    REASON: EZH2 elevation identifies the LOCK
    is present. FOXA1 suppression confirms
    the LOCK IS ON THE RIGHT TARGET
    (FOXA1 locus is H3K27me3-methylated by EZH2).
    A TNBC patient with high EZH2 but FOXA1
    partially retained (shallow TNBC) will
    respond to EZH2i with faster FOXA1 restoration
    than a deep TNBC with maximal FOXA1 suppression.
    Composite selector: EZH2-high + FOXA1-low
    (not EZH2 alone).

FALSIFICATION:
  If FOXA1 does not predict EZH2i response
  better than EZH2 alone in existing TNBC
  datasets — fails.
  If ILC FOXA1-hyperactivation does not
  distinguish ILC from other subtypes on the
  same scale — fails.

────────────────────��────────────────────

PREDICTION CS-8 (LOCKED 2026-03-05):

THE CROSS-SUBTYPE DRUG MAP HAS A SINGLE
UNDERLYING LOGIC — DEPTH-GUIDED TARGETING.

The complete drug map for breast cancer,
derived from attractor geometry, follows
one principle:

  The target is the mechanism maintaining
  the specific depth of false attractor the
  cell is in.

  Match the drug to the lock at that depth.
  Not to the surface phenotype.
  Not to the receptor status alone.
  To the LOCK at that DEPTH.

Cross-subtype drug map (predicted):

  DEPTH 0–1 (Normal):
    No intervention required.

  DEPTH 1 (LumA):
    Lock: CDK4/6 kinase
    Drug: CDK4/6i + ET
    Logic: Remove the cell cycle block that
           prevents terminal exit from the
           attractor. The TFs are present.
           The luminal programme is active.
           Release the brake.

  DEPTH 2 (LumB):
    Lock: HDAC1/2 + DNMT3A (chromatin)
    Drug: HDACi (entinostat) + CDK4/6i + ET
    Logic: HDAC1/2-DNMT3A co-complex is
           silencing ER target genes at
           the chromatin level despite
           intact receptor. Dissolve the
           chromatin lock THEN release the
           cycle brake THEN block the receptor.
           Sequence: HDACi → CDK4/6i + ET.

  DEPTH 3 (HER2):
    Lock: ERBB2 amplicon (copy number)
    Drug: Anti-HER2 (trastuzumab primary)
          + EZH2i for deep fraction
    Logic: The amplicon is the constitutive
           signal that prevents exit. Block
           the amplified kinase. For the deep
           fraction where EZH2 is also rising:
           dissolve the secondary epigenetic
           lock before the anti-HER2 can
           fully engage FOXA1-retained programme.

  DEPTH 4 (TNBC):
    Lock: EZH2/PRC2 (dominant epigenetic)
    Drug: Tazemetostat → Fulvestrant sequence
          + PARPi combination for BRCA1 subtype
          + Ferroptosis (GPX4i) for deep-end ZEB1-high
    Logic: EZH2 has silenced FOXA1/GATA3/ESR1.
           Dissolve the H3K27me3 marks (tazemetostat).
           Allow FOXA1 to re-bind its loci.
           GATA3 and ESR1 re-activate.
           NOW target the re-expressed ESR1
           (fulvestrant). The sequence matters —
           ESR1 must be present to be targeted.
           Converting TNBC to a luminal-like state
           before applying endocrine therapy.

  DEPTH 5 (Claudin-low):
    Lock: Pre-commitment programme
          (below the level where any luminal
          target is present to engage)
    Drug: Memory-low subtype:
            Anti-TIGIT (Treg depletion)
            → anti-PD-1 (checkpoint release)
            → CT antigen targeting (GAGE)
          Memory-high/both:
            Commitment forcing (BETi, FOXA1
            activator — to restore any identity
            before depth-targeting can apply)
          Depth-high CL:
            CD44/TGFβ targeting
    Logic: Luminal programme is below the
           commitment threshold. Cannot target
           ESR1 because ESR1 doesn't exist
           in these cells. Cannot target FOXA1
           because FOXA1 doesn't exist.
           Must either: (A) exploit the
           pre-commitment state's immune
           vulnerability (CT antigens exposed),
           or (B) force commitment to ANY
           identity before applying subtype-
           specific therapy.

  ILC (orthogonal axis):
    Lock: CDH1 structural loss
    Drug: ET + CDK4/6i (luminal target valid)
          EZH2i for deep escape (MKI67-high,
          EZH2-high ILC)
    Logic: ILC is not on the same depth axis
           as the ductal subtypes. It has the
           MOST luminal identity of any breast
           cancer. The treatment target is
           the tumour's own hyperactivated
           luminal programme (ET) plus the
           cell cycle escape mechanism (CDK4/6i).
           EZH2i for the subset where epigenetic
           escape is producing a secondary
           depth within ILC.

FALSIFICATION:
  If any TNBC subpopulation shows treatment
  response to ET without prior EZH2i — the
  sequence prediction is wrong (ESR1 may be
  residually expressed in shallow TNBC).
  If commitment forcing in claudin-low fails
  to produce any luminal-like gene expression
  — the deepest drug prediction is wrong.
  If HDACi shows no differential benefit in
  LumB vs LumA in clinical data — the chromatin
  lock prediction for LumB is wrong.
```

---

### AXIS 4 PREDICTIONS — NOVEL CROSS-SUBTYPE FINDINGS

---

```
PREDICTION CS-9 (LOCKED 2026-03-05):

THE SINGLE MOST IMPORTANT CLINICAL FINDING
FROM THE CROSS-SUBTYPE ANALYSIS WILL BE:

FOXA1 AS A UNIVERSAL BREAST CANCER ATTRACTOR
DEPTH BIOMARKER — MORE INFORMATIVE THAN
PAM50 SUBTYPE ALONE.

Reasoning:
  PAM50 assigns a discrete label (LumA, LumB, etc).
  FOXA1 protein level is continuous.
  FOXA1 tracks the depth axis continuously
  across all ductal subtypes.

  The predicted clinical utility:
    FOXA1 IHC level within ANY subtype predicts
    which depth-guided drug is most appropriate.

    Within ER+ breast cancer:
      FOXA1-high, CDKN1A-high → LumA geometry
      → CDK4/6i + ET, low intensity
      FOXA1-high, TFF1-low → LumB geometry
      → HDACi + CDK4/6i + ET, HDACi priority
      FOXA1-intermediate + ERBB2-amplified
      → HER2 geometry → trastuzumab primary

    Within TNBC:
      FOXA1-absent + EZH2-high
      → deep TNBC → tazemetostat first
      FOXA1-partial + EZH2-intermediate
      → shallow TNBC → may be convertible
        more rapidly by EZH2i

    Within claudin-low:
      FOXA1-absent + GATA3-absent
      → memory-low TYPE 4 → anti-TIGIT candidate
      FOXA1-partial + GATA3-partial
      → memory-high TYPE 4 → different immune logic

  A single FOXA1 IHC result, combined with
  EZH2 IHC and CDH1 IHC, produces a
  three-marker geometric classification of
  breast cancer depth that is more treatment-
  relevant than PAM50 alone.

  THREE-MARKER PANEL PREDICTION:
    FOXA1 + EZH2 + CDH1 IHC determines:
      FOXA1 high, EZH2 low, CDH1 present
      → LumA geometry → CDK4/6i + ET
      FOXA1 high, EZH2 intermediate, CDH1 absent
      → ILC geometry → ET + CDK4/6i ± EZH2i
      FOXA1 intermediate, EZH2 high, CDH1 present
      → HER2 or deep LumB geometry → EZH2i priority
      FOXA1 absent, EZH2 very high, CDH1 present
      → TNBC geometry → tazemetostat sequence
      FOXA1 absent, GATA3 absent, CDH1 present
      → claudin-low TYPE 4 → anti-TIGIT candidate
        (if FOXP3/CD8A high)

FALSIFICATION:
  If FOXA1 IHC does not correlate with PAM50
  depth ordering in TCGA-BRCA bulk data —
  the three-marker panel prediction fails.
  If the ordering is present but non-monotonic —
  partially fails.

─────────────────────────────────────────

PREDICTION CS-10 (LOCKED 2026-03-05):

THE PROGNOSIS GRADIENT ACROSS BREAST CANCER
FOLLOWS THE ATTRACTOR DEPTH SPECTRUM,
NOT THE CLINICAL RECEPTOR STATUS.

The current clinical staging system uses:
  ER/PR status (endocrine therapy eligibility)
  HER2 status (anti-HER2 eligibility)
  Ki-67 (proliferation proxy)
  Grade (differentiation proxy)

These four variables capture PARTS of the
attractor depth axis but not its full structure.

Specific predictions:

  CS-10a: Within ER+ breast cancer, CDKN1A
          (p21) level will predict OS BETTER
          than Ki-67 alone.
          Reason: CDKN1A is the primary depth
          variable in LumA/LumB (confirmed in
          individual analyses). Ki-67 captures
          proliferation but not the underlying
          attractor depth that causes the
          CDKN1A depletion.

  CS-10b: Within TNBC, the composite TNBC
          depth score (AR + ZEB1 + SOX10
          composite) will predict pCR rate
          better than individual markers.
          Reason: Individual markers capture
          one dimension of depth. The composite
          captures the full false attractor
          membership.

  CS-10c: Within HER2-enriched, the three-marker
          deep-end composite (AR-low + ERBB3-low
          + CDH1-low) will predict trastuzumab
          resistance better than ERBB2 amplicon
          level alone.
          Reason: The amplicon level determines
          subtype membership. The depth within
          the HER2 false attractor determines
          treatment response.

  CS-10d: Within ILC, MKI67 level will predict
          late recurrence (>10 years) better
          than grade alone.
          Reason: ILC late recurrence is driven
          by the small MKI67-high proliferative
          escape subpopulation. Grade is a
          population-level descriptor. MKI67
          tracks the escape cells.

FALSIFICATION:
  These are testable using existing TCGA-BRCA
  clinical data with survival endpoints.
  If CDKN1A is not more predictive than Ki-67
  in ER+ TCGA-BRCA — CS-10a fails.
  If composite TNBC depth does not outpredict
  any individual marker — CS-10b fails.

─────────────────────────────────────────

PREDICTION CS-11 (LOCKED 2026-03-05):

THE METASTATIC TRANSITION IS A DEPTH SHIFT
WITHIN THE EXISTING ATTRACTOR LANDSCAPE —
NOT A TRANSITION TO A NEW ATTRACTOR TYPE.

Specifically:
  Metastatic LumA cells will show LOWER FOXA1
  and HIGHER EZH2 than primary LumA cells
  from the same tumour. They will geometrically
  resemble LumB or early HER2-like geometry.

  Metastatic LumB cells will show LOWER FOXA1,
  LOWER TFF1, HIGHER SOX10, and will
  geometrically approach the TNBC attractor.

  Metastatic TNBC cells will show HIGHER ZEB1/
  SNAI1, LOWER AR, and will geometrically
  approach claudin-low geometry.

  The metastatic cell is a DEEPER VERSION of
  its primary tumour's attractor type —
  not a new type.

  PRACTICAL CONSEQUENCE:
    Drug resistance and metastasis are
    GEOMETRIC PHENOMENA — both involve
    depth increase in the existing false
    attractor.
    The drug target for metastatic disease
    should be the DEPTH-APPROPRIATE target
    for where the metastatic cell NOW sits —
    which may be one level deeper than the
    primary tumour.

    A patient whose primary was LumA and whose
    metastasis is now geometrically LumB-like
    should have HDACi added to the CDK4/6i + ET
    backbone — not just dose escalation of
    the CDK4/6i.

    This is the cross-subtype mechanism for
    the clinical phenomenon of acquired
    endocrine resistance: the cell has
    shifted one attractor depth level.

FALSIFICATION:
  If metastatic LumA cells do NOT show
  higher EZH2 or lower FOXA1 than primary
  LumA cells in matched datasets — fails.
  If metastatic cells are not geometrically
  between their primary type and the next
  deeper type — fails.
  (Test: Doc 107 — BRCA chemo resistance
  dataset GSE161533 planned in TODO.)

─────────────────────────────────────────

PREDICTION CS-12 (LOCKED 2026-03-05):

ILC IS THE STRUCTURAL INVERSE OF
THE DUCTAL DEPTH SPECTRUM —
AND THIS INVERSION HAS A CLINICAL CONSEQUENCE
THAT IS NOT CURRENTLY RECOGNISED IN GUIDELINES.

ILC is the only breast cancer subtype where:
  Luminal TFs are hyperactivated above normal.
  CDH1/E-cadherin is lost.
  The patient's tumour is MORE committed to
  luminal identity than a normal cell,
  but has lost the ability to form cohesive
  architecture.

The clinical consequence NOT in current guidelines:

  Anti-HER2 therapy has NO GEOMETRIC BASIS
  for ILC as a class (confirmed ILC lit check).
  The ILC biology is driven by hyperactivated
  luminal TFs — the opposite of HER2-enriched
  biology. When ILC is HER2+ by IHC, it is
  because of an incidental co-occurring
  amplification, not because the ILC false
  attractor shares geometry with the HER2-enriched
  false attractor. Anti-HER2 targets the
  amplicon-driven constitutive signal.
  ILC does not have an amplicon-driven
  constitutive signal in its core geometry.
  Anti-HER2 in ILC will not dissolve the
  ILC false attractor.

  The ILC-specific treatment principle
  (not in guidelines):
    ET is the PRIMARY target because the
    hyperactivated FOXA1/GATA3/ESR1 programme
    is the most targetable feature of ILC.
    Blocking ESR1 directly (fulvestrant over
    aromatase inhibitors) may be more effective
    in ILC because the luminal programme is
    more active, not less.
    CDK4/6i addresses the CCND1-driven
    proliferative component.
    EZH2i specifically for MKI67-high
    EZH2-high ILC where epigenetic escape
    has produced a secondary depth within
    the ILC false attractor.

FALSIFICATION:
  If ILC clinical data shows equivalent
  anti-HER2 benefit to HER2-amplified ductal
  tumours — CS-12 needs revision.
  If fulvestrant is not superior to AIs
  in ILC outcomes (specifically in ESR1-high,
  FOXA1-hyperactive ILC) — partial failure.
  (NOTE: The ILC fulvestrant prediction
  has some clinical support from the FALCON
  trial subgroup analysis — partially
  convergent already.)

─────────────────────────────────────────

PREDICTION CS-13 (LOCKED 2026-03-05):

THE TNBC PARADOX (HIGH pCR RATE +
HIGH EARLY RECURRENCE RATE) IS EXPLAINED
BY THE DEPTH AXIS.

Confirmed individually in TNBC analysis:
The pCR/DRFS paradox is explained by the
depth score as a single variable.

Cross-subtype prediction:
  The TNBC paradox is the predictable
  consequence of the TYPE 2 attractor geometry
  in the context of conventional chemotherapy.

  Shallow TNBC (low depth score, AR-intermediate,
  ZEB1-intermediate, FOXA1-partially-retained):
    These cells are sensitive to chemotherapy.
    Anthracycline/taxane kills proliferating cells.
    Low depth → cells are closer to a committed
    state → more sensitive to cytotoxic agents
    that exploit cell cycle machinery.
    pCR rate is HIGH in this group.
    Post-pCR, few residual cells → low recurrence.

  Deep TNBC (high depth score, AR-absent,
  ZEB1-high, FOXA1-absent, SOX10-high):
    These cells are in a deep false attractor.
    They cycle differently and more slowly.
    Chemotherapy kills the shallow cells first.
    The deep cells SURVIVE because they are
    less proliferative and more EMT-like.
    pCR rate is LOW in this group.
    But the surviving deep cells recur EARLY —
    they are aggressive and chemo-resistant.

  The clinical paradox (TNBC has high pCR rates
  AND high early recurrence rates
  simultaneously):
    Both are true — for DIFFERENT patient
    subpopulations within TNBC.
    The shallow TNBC achieves pCR and rarely
    recurs. The deep TNBC fails to achieve pCR
    and recurs early. When averaged across all
    TNBC, the population appears paradoxically
    both chemo-sensitive and aggressive.

  CROSS-SUBTYPE IMPLICATION:
    This same paradox does NOT appear in LumA
    (shallow attractor — chemotherapy is rarely
    used because hormone therapy suffices).
    It does NOT appear in HER2-enriched
    (the amplicon is directly targeted by
    trastuzumab — depth matters less when the
    primary lock is a copy number event).
    It DOES appear in claudin-low (where pCR
    rate is low and recurrence is early —
    consistent with deep attractor geometry
    where no chemotherapy can cross the
    commitment threshold).

FALSIFICATION:
  If depth score does not separate TNBC into
  high-pCR (low depth) and low-pCR (high depth)
  subgroups — fails.
  (Already partially confirmed in individual
  TNBC analysis — GSE25066 DRFS log-rank
  p<0.0001. Cross-subtype test extends this
  to the full breast cancer landscape.)
```

---

### AXIS 5 PREDICTIONS — THE INDIVIDUAL PATIENT PROTOCOL

---

```
PREDICTION CS-14 (LOCKED 2026-03-05):

THE CROSS-SUBTYPE ANALYSIS WILL PRODUCE
THE SINGLE MOST CLINICALLY USEFUL DOCUMENT
IN THE BRCA CHAPTER — A DEPTH-GUIDED
TREATMENT ALGORITHM APPLICABLE TO EVERY
BREAST CANCER PATIENT.

The output of this analysis — a unified
geometric classification with a depth-guided
drug map across all six subtypes — will be
the primary reference document for every
individual patient analysis in the BRCA
individual patient protocol.

When a breast cancer patient sends their data:
  Step 1: Determine subtype by geometry
          (FOXA1 + EZH2 + CDH1 + ERBB2 IHC)
  Step 2: Determine depth within subtype
          (subtype-specific depth score)
  Step 3: Apply depth-guided drug map
          from this document
  Step 4: Generate clinical questions
          specific to their position on the map

This is the individual patient protocol
for breast cancer.
Every prediction in this analysis either
confirms, refines, or falsifies the drug map
that individual patient reports will use.

WHAT MAKES THIS PREDICTION FALSIFIABLE:
  The utility claim is falsifiable only by
  outcome data from individual patients.
  The structural claim — that a unified
  depth-guided map exists and is derivable
  from the cross-subtype analysis — is
  confirmed or falsified by whether the
  predictions CS-1 through CS-13 hold.

  If ≥ 8 of 13 structural predictions confirm:
    The map is valid and the individual patient
    protocol is grounded.

  If < 8 confirm:
    The depth spectrum is not continuous
    or not as ordered as predicted.
    The map requires revision before use
    in individual patient reports.
    The revision is documented and the
    revised map is used.

─────────────────────────────────────────

PREDICTION CS-15 (LOCKED 2026-03-05):

THE THREE-MARKER IHC PANEL
(FOXA1 + EZH2 + CDH1) CLASSIFIES EVERY
BREAST CANCER PATIENT INTO A DEPTH TIER
THAT IS MORE TREATMENT-RELEVANT THAN
PAM50 SUBTYPE ALONE — AND DOES SO USING
ASSAYS THAT ALREADY EXIST IN EVERY
PATHOLOGY LAB.

This is the most clinically translatable
prediction in the entire BRCA chapter.

  FOXA1 IHC:  Already available — used in
              routine clinical breast pathology
              to distinguish ER+ subtypes
              and characterise ILC.
  EZH2 IHC:   Already available — used in
              multiple cancer types as a
              prognostic marker and research
              tool. No new assay development.
  CDH1 IHC:   Already available — standard
              clinical IHC in breast cancer
              to identify ILC (CDH1 loss
              is the diagnostic IHC criterion
              for ILC).

  The three tests together cost less than
  an additional $200 per patient in a
  standard pathology workflow.
  They are available in every oncology
  centre that performs breast cancer
  pathology.
  They do not require RNA-seq.
  They do not require genomic profiling.
  They do not require any new technology.

  THE PATIENT SELECTION TABLE PREDICTED:

  FOXA1  EZH2    CDH1  → Geometry
  ─────────────────────────────────────────
  High   Low     Pres  → LumA: CDK4/6i+ET
  High   High    Abs   → ILC:  ET+CDK4/6i
                               ±EZH2i (if MKI67 high)
  Int    High    Pres  → LumB/HER2: check ERBB2
                         If ERBB2 amplified → anti-HER2
                         If ERBB2 normal → HDACi+CDK4/6i+ET
  Abs    V.High  Pres  → TNBC: tazemetostat sequence
  Abs    Abs     Pres  → Claudin-low: FOXP3/CD8A IHC
                         next → anti-TIGIT if FOXP3 high

  (Abbreviations: Pres=present, Abs=absent,
   Int=intermediate, V.High=very high)

  Adding FOXP3 and CD8A (already validated as
  clinical IHC markers in TNBC) completes the
  claudin-low patient selection algorithm for
  anti-TIGIT therapy.

  Five IHC markers.
  Five existing clinical-grade antibodies.
  A geometric classification of breast cancer
  depth applicable to any pathology lab today.

FALSIFICATION:
  If the three-marker panel does not
  distinguish the six subtype geometries
  when applied to TCGA-BRCA with known PAM50
  subtype calls — the panel utility fails.
  If FOXA1 is not separating LumA from TNBC
  in bulk TCGA data — fails (very unlikely
  given confirmed individual analyses, but
  testable).
```

---

## PART III — WHAT THE SCRIPT MUST PRODUCE
### Technical Requirements for BRCA-S8b

---

```
The cross-subtype script (BRCA-S8b) must
produce the following outputs to test
the 15 predictions above:

REQUIRED OUTPUTS:

1. UNIFIED DEPTH AXIS PLOT
   All six subtypes placed on a single axis
   using the luminal TF composite score
   (FOXA1 + GATA3 + ESR1 mean z-score,
   normalised to Mature Luminal = 0).
   Reference populations:
     Mature Luminal (n=1265, GSE176078)
     Luminal Progenitor (n=1992, GSE176078)
   Predicted order: ILC > LumA > LumB >
                    HER2 > TNBC > Claudin-low
   Tests: CS-1

2. EZH2 CROSS-SUBTYPE ELEVATION TABLE
   EZH2 mean expression per subtype
   normalised to Mature Luminal.
   All six subtypes + normal reference.
   Predicted order: CL > TNBC > HER2 >
                    LumB > ILC > LumA
   Tests: CS-2

3. PCA GEOMETRY — SIX POPULATIONS
   Same PCA coordinate system as HER2 Script 1
   (already partially done — all six subtype
   centroids computed in same PCA space).
   Centroid distances from Mature Luminal.
   Predicted order: ILC < LumA < HER2 <
                    LumB < TNBC < CL
   Tests: CS-3

4. ATTRACTOR TYPE CLASSIFICATION TABLE
   Formal classification of each subtype
   using the Attractor Geometry Axioms.
   Diagnostic test applied: Identity TF
   Direction Test (suppressed vs elevated).
   Tests: CS-4

5. LOCK MECHANISM TABLE
   Primary lock mechanism per subtype.
   Drug class matched to lock.
   Tests: CS-5

6. EZH2 DEPTH-STRATIFIED DRUG PRIORITY TABLE
   EZH2i priority score per subtype.
   Patient selection criterion per subtype.
   Tests: CS-6

7. FOXA1 THERAPEUTIC READINESS ANALYSIS
   FOXA1 level as continuous predictor
   across all six subtypes.
   Correlation with EZH2.
   Correlation with clinical depth proxy (MKI67).
   Tests: CS-7, CS-9

8. DEPTH-GUIDED DRUG MAP — UNIFIED TABLE
   One table. All six subtypes.
   Depth level, primary lock,
   first-line drug, second-line drug,
   patient selection tool, validation status.
   Tests: CS-8

9. FOXA1 + EZH2 + CDH1 THREE-MARKER PANEL TEST
   Simulate the three-marker IHC classification
   using TCGA-BRCA expression data as IHC proxy.
   Show separation of PAM50 subtypes using
   the three-marker rule.
   Tests: CS-15

10. METASTASIS DEPTH SHIFT PREDICTION SETUP
    Summary of what Doc 107 will test.
    Not run in this script — flags for Script 5.
    Tests: CS-11 setup only

DATASET:
  Primary: GSE176078 (scRNA-seq, already used)
           All six subtype cell populations
           identified using existing metadata.
  Secondary: TCGA-BRCA (bulk RNA-seq)
           PAM50 subtype calls available.
           Used for CS-9, CS-10, CS-15.

REFERENCE POPULATIONS (from GSE176078):
  Normal reference:   Mature Luminal (n=1265)
  Progenitor:         Luminal Progenitors (n=1992)
  Cancer populations:
    Cancer LumA SC    (n=7742)
    Cancer Basal SC   (n=4312)
    Cancer Her2 SC    (confirmed in metadata)
    ILC cells         (identified in ILC script)
    Claudin-low cells (identified in S7 series)
    LumB cells        (identified in LumB script)
```

---

## PART IV — PREDICTION SCORECARD TEMPLATE
### To Be Filled in BRCA-S8b After Script Runs

---

```
CROSS-SUBTYPE PREDICTION SCORECARD
To be completed in BRCA-S8b (results document)
Date locked: 2026-03-05

CS-1:  Depth spectrum continuous and ordered    □ ✓  □ ✗  □ ?
CS-2:  EZH2 tracks depth monotonically         □ ✓  □ ✗  □ ?
CS-3:  PCA distance tracks prognosis           □ ✓  □ ✗  □ ?
CS-4:  Six subtypes → four attractor types     □ ✓  □ ✗  □ ?
CS-5:  Lock mechanism changes with depth       □ ✓  □ ✗  □ ?
CS-6:  EZH2i priority depth-stratified         □ ✓  □ ✗  □ ?
CS-7:  FOXA1 = master therapeutic variable     □ ✓  □ ✗  □ ?
CS-8:  Unified depth-guided drug map valid     □ ✓  □ ✗  □ ?
CS-9:  FOXA1 better predictor than EZH2 alone  □ ✓  □ ✗  □ ?
CS-10: Depth better than standard markers      □ ✓  □ ✗  □ ?
CS-11: Metastasis = depth shift (setup only)   □ N/A for this script
CS-12: ILC is structural inverse — clinical Δ  □ ✓  □ ✗  □ ?
CS-13: TNBC paradox explained by depth         □ ✓  □ ✗  □ ?
CS-14: Map clinically usable (≥8/13 confirmed) □ ✓  □ ✗  □ ?
CS-15: Three-marker IHC panel separates types  □ ✓  □ ✗  □ ?

OVERALL SCORECARD:    __ / 14 tested
INDIVIDUAL PROTOCOL:  Ready / Needs revision
```

---

## DOCUMENT STATUS

```
document:     BRCA_Cross_Subtype_Before.md
folder:       Cancer_Research/BRCA/DEEP_DIVE/Cross_Subtype/
document_id:  BRCA-S8a
version:      1.0
date:         2026-03-05
author:       Eric Robert Lawson
              OrganismCore
status:       COMPLETE — PREDICTIONS LOCKED
              Do not modify after this date.
              Results document: BRCA-S8b
              Script:           BRCA_Cross_Subtype_Script.py

predictions_count:  15
predictions_locked: 2026-03-05
data_examined:      NONE (before-document)

governing_documents:
  Individual_Protocol/ETHICS.md
  Attractor_Geometry_Axioms.md (Document 90)
  Workflow_Protocol.md v2.0

founding_principle:
  "I do not want to take people's money
   and promise them bullshit."
   — Eric Robert Lawson, March 4, 2026
```

---

*"The map is drawn one subtype at a time.*
*The cross-subtype analysis is the moment*
*the six maps become one.*
*That map is what every breast cancer patient*
*in the individual protocol will be located on.*
*Get it right."*

— Eric Robert Lawson, OrganismCore, 2026-03-05
