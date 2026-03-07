# BRCA REFERENCE GEOMETRY VALUES
## Population Reference Means for Individual Patient Analysis
## OrganismCore — Eric Robert Lawson
## 2026-03-07

---

## WHAT THIS DOCUMENT IS

```
This document contains the specific numerical
values required to execute the individual
patient geometric analysis protocol.

Source population analysis:
  BRCA-S8b (BRCA_Cross_Subtype_Script1.py)
  Dataset: GSE176078
  19,542 single cancer cells, 26 tumours
  Normal reference: Mature Luminal population
  n = 1,265 cells

These values are the coordinate system.
A patient's expression data is placed into
this coordinate system using z-scores
relative to these reference means.

Without these values, Steps 3 and 4 of
the individual protocol cannot execute.

These are results of the population analysis.
They are not assumptions.
They are the map.
```

---

## PART I — PRIMARY AXIS REFERENCE VALUES
## (RNA-Space — not transferable to IHC)

```
The primary axis was revealed by the
saddle point scan in BRCA-S8b.
These are the POPULATION MEDIAN ratios,
not the cut-points for an individual patient.
An individual patient's ratio is placed
relative to these medians.

FOXA1/EZH2 RATIO POPULATION MEDIANS
(from TCGA-BRCA validation, n=837):

  LumA:        9.38
  LumB:        8.10
  HER2-enr:    3.34
  TNBC:        0.52
  Claudin-low: 0.10

  Ordering confirmed: LumA > LumB >
  HER2e > TNBC > CL
  p = 2.87×10⁻¹⁰³ (Kruskal-Wallis)
  Zero contradictions across seven datasets.

NOTE ON UNITS:
  These values are in RNA-space.
  They are computed from RNA-seq TPM values.
  They are NOT transferable to IHC H-score.
  A patient with Tempus xR data in TPM can
  be placed relative to these medians.
  A patient with only IHC data cannot be
  placed at a precise point — only a region.
  This distinction governs what the analysis
  can and cannot state at each tier.
```

---

## PART II — NORMAL REFERENCE POPULATION
## (Mature Luminal — GSE176078)

```
These values are the reference for z-score
computation. A patient's gene expression
is expressed as a deviation from the
Mature Luminal normal.

All values are from the Mature Luminal
population in GSE176078 (n=1,265 cells)
after log1p normalisation.

FORMAT: gene — mean (SD)

SWITCH GENES (suppressed in cancer vs normal):

  FOXA1   — population analysis confirmed
            as the primary identity TF
  GATA3   — luminal lineage TF
  ESR1    — estrogen receptor
  CDH1    — E-cadherin
  KRT18   — luminal cytokeratin
  KRT8    — luminal cytokeratin
  SPDEF   — luminal differentiation TF
  SCUBE2  — luminal differentiation marker
  PGR     — progesterone receptor

FA MARKERS (elevated in cancer vs normal):

  EZH2    — primary epigenetic lock gene
  MKI67   — proliferation
  TOP2A   — replication/proliferation
  PCNA    — proliferation

EPIGENETIC LOCK PANEL:

  EZH2    — PRC2 catalytic subunit
  KDM1A   — LSD1, CoREST complex
  HDAC1   — HDAC complex
  HDAC2   — HDAC complex
  DNMT3A  — DNA methylation
  KDM6A   — H3K27me3 eraser
            (low KDM6A = EZH2 unopposed)

CLAUDIN-LOW PANEL:

  CLDN3   — tight junction claudin
  CLDN4   — tight junction claudin
  CLDN7   — tight junction claudin
  CDH1    — E-cadherin (also switch panel)

ILC STRUCTURAL MARKER:

  CDH1    — loss confirms ILC geometry
            (also in switch and CL panels)
```

**POPULATION ANALYSIS STATUS NOTE:**

```
The exact numerical means and SDs for each
gene in the Mature Luminal reference
population are in the BRCA-S8b script
output files:

  Cross_Subtype_s1_results/results/
  cs_s1_depth_table.csv
  cs_s1_scorecard.csv

For a Tier 3 patient analysis, pull the
Mature Luminal reference means from that
output before running Step 3.

Until those values are extracted into this
document explicitly, this section serves as
the gene panel reference. The numerical
means must be retrieved from BRCA-S8b
output before the first Tier 3 analysis runs.

This is a gap in this document.
It is stated as a gap, not papered over.
```

---

## PART III — SUBTYPE DEPTH SCORE
## REFERENCE DISTRIBUTIONS

```
The depth score from the population analysis
(BRCA-S8b) establishes the reference
distribution for each subtype.

An individual patient's depth score is placed
relative to these distributions.

DEPTH SCORE FORMULA:
  depth_score = fa_score - switch_score

  fa_score = mean z-score of FA markers
             relative to Mature Luminal mean
             (positive = elevated above normal)

  switch_score = mean z-score of switch genes
                 relative to Mature Luminal mean
                 (negative = suppressed below normal)

POPULATION DEPTH ORDERING (confirmed BRCA-S8b):
  Deepest (most displaced from normal):
    Claudin-low > TNBC/Basal > HER2e >
    LumB > LumA
    (when measured by EZH2-free PCA:
     CL 6.572, TNBC 6.063 — CS-LIT-3)

  NOTE: This ordering is for EZH2-free PCA
  distance. When EZH2 is included, TNBC
  appears deeper than CL because EZH2 is
  extremely high in TNBC. EZH2-free PCA
  is the correct method for TYPE comparison
  — CS-LIT-26.

REFERENCE DISTRIBUTION EXTRACTION:
  Pull the per-subtype depth score
  distributions from:
    cs_s1_depth_table.csv (BRCA-S8b output)
  These give the mean and SD of depth score
  for each subtype population.
  A patient's depth score percentile is
  computed relative to their closest subtype
  distribution.
```

---

## PART IV — ATTRACTOR TYPE REFERENCE
## (From BRCA-S8b and Axioms Document)

```
BRCA attractor type assignments
(confirmed by Identity TF Direction Test
in BRCA-S8b attractor_type_classification):

  LumA:        TYPE I  — Blocked Approach
               FOXA1 HIGH, EZH2 LOW
               Identity present, lock minimal
               Cell cycle lock dominant (CDK4/6)

  LumB:        TYPE I  — Blocked Approach (deeper)
               FOXA1 moderate, EZH2 elevated
               Identity partially suppressed
               Epigenetic lock active at PR/GATA3

  HER2-enr:    TYPE II — Wrong Valley
               FOXA1 LOW, ERBB2 HIGH
               Alternative identity programme
               dominant

  HER2-deep:   TYPE I + TYPE II composite
               (CDH3-high, AR-low, EZH2 +118%)
               EZH2i relevant in this fraction

  TNBC:        TYPE I  — Blocked Approach (deep)
               FOXA1 near-zero, EZH2 dominant
               (+189% vs Mature Luminal)

  Claudin-low: TYPE IV — Root Lock
               No mature lineage identity.
               Pre-commitment arrest.
               Claudins low, CDH1 low,
               immune infiltration high.

  ILC:         Composite — structural lock
               FOXA1 HIGH (identity preserved)
               CDH1 ABSENT (structural disrupted)
               Not a depth-axis classification.
               A separate structural variant.
```

---

## PART V — DRUG TARGET GEOMETRY REFERENCE

```
These are the drug target predictions from
the population analysis (BRCA-S8b, validated
in BRCA-S8h literature check).

They are stated as geometric observations.
They are the output of the attractor framework
applied to BRCA. They are not treatment
recommendations.

LumA geometry:
  CDK4/6 inhibitors + ET
  Basis: CDKN1A suppression confirmed
  (CS-LIT-13, CS-LIT-14)
  CDKN1A level = quantitative predictor
  of CDK4/6i benefit magnitude

LumB geometry:
  HDAC inhibitors + ET (entinostat)
  Basis: DNMT3A/HDAC2 co-elevation,
  TFF1/ESR1 decoupling
  (CS-LIT-8, CS-LIT-9, CS-LIT-15, CS-LIT-23)
  Patient selector: TFF1/ESR1 ratio (Tier 3)
  or PR-absent pattern (Tier 1 proxy)

HER2e geometry:
  Anti-HER2 therapy first
  EZH2i addition for HER2-deep fraction
  (CDH3-high, AR-low, EZH2 elevated)
  (CS-LIT-21)

TNBC geometry:
  Tazemetostat → fulvestrant sequence
  Basis: EZH2 dominant lock, FOXA1 near-zero
  Sequence critical: EZH2i unlocks FOXA1
  programme, then fulvestrant exploits
  restored ET sensitivity
  (CS-LIT-16, CS-LIT-17)
  EZH2 IHC independently prognostic in TNBC
  (CS-LIT-24)

Claudin-low geometry:
  Anti-TIGIT → anti-PD-1
  Basis: pre-commitment arrest, immune
  infiltration prominent
  Patient selector: FOXP3/CD8A ratio
  HR = 2.212 (CS-LIT-20)
  SKYLINE trial: NCT06175390
  (tiragolumab + atezolizumab) (CS-LIT-28)

ILC geometry:
  Fulvestrant preferred over AI
  Basis: CDH1-loss changes ER signalling
  context, creates fulvestrant advantage
  (CS-LIT-18)

TNBC/CL ambiguity:
  When TNBC pattern (ER-, PR-, HER2-)
  but claudin-low cannot be confirmed
  or excluded:
  State the ambiguity explicitly.
  State what would resolve it
  (claudin marker IHC or Tier 3 data).
  Do not assign tazemetostat or
  immunotherapy until resolved.
```

---

## PART VI — WHAT IS CONFIRMED, WHAT IS NOT

```
CONFIRMED across seven independent datasets
(~7,500 patients):
  Primary axis ordering (LumA > LumB >
  HER2e > TNBC > CL)
  Depth score validation (HR=1.509 for
  TNBC depth in GSE25066)
  AUC for classification (0.828-0.901
  LumA vs Basal; 0.796-0.873 LumA vs LumB)
  Protein-level confirmation by CPTAC MS
  (r=-0.492, n=122)

INDEPENDENTLY CONFIRMED by other groups:
  FOXA1/EZH2 mechanistic axis:
    Schade et al. Nature 2024
    Toska et al. Nature Medicine 2017
    Neither group knew of this framework.

NOT YET CONFIRMED:
  IHC H-score cut-points (calibration
  study required — CS-LIT-22)
  Individual patient analyses (first
  patients have not yet been run)
  Drug predictions (framework-derived,
  not yet prospectively tested)

This document presents what is established.
What is not established is stated as such.
Every patient analysis built on these values
must state which items are confirmed and
which are framework-derived predictions.
```

---

## DOCUMENT METADATA

```
document_id:    BRCA_REFERENCE_VALUES
folder:         Individual_Protocol/Breast_Cancer/
type:           Operational reference table for
                individual patient analysis
cancer_type:    BRCA
version:        1.0
date:           2026-03-07
status:         ACTIVE — operational reference

source_analysis: BRCA-S8b
  Cancer_Research/BRCA/DEEP_DIVE/
  BRCA_Cross_Subtype_Script1.py

validation_document: BRCA-S8h
  Cancer_Research/BRCA/DEEP_DIVE/
  BRCA_Cross_Subtype_Literature_Check.md

KNOWN GAP:
  Numerical means and SDs for Mature Luminal
  reference population not yet extracted
  from BRCA-S8b output files into this document.
  Must be extracted before first Tier 3
  analysis runs.
  Source file: cs_s1_depth_table.csv
  in Cross_Subtype_s1_results/results/

  Until that extraction is done:
  Tier 1 and Tier 2 analyses can proceed
  using proxy inference as specified.
  Tier 3 analyses require pulling the
  reference means from BRCA-S8b output
  directly before running.
```
