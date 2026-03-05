# FOXA1/EZH2 RATIO — SCRIPT 1 BEFORE-DOCUMENT
## Predictions Locked Before Any Validation Script Runs
## OrganismCore — Document RATIO-S1a
## Date: 2026-03-05

---

## DOCUMENT METADATA

```
document_id:        RATIO-S1a
series:             FOXA1/EZH2 Ratio Validation
folder:             Cancer_Research/BRCA/DEEP_DIVE/FOXA1_EZH2_RATIO/
type:               BEFORE-DOCUMENT
                    All predictions stated and locked
                    before any validation script runs.
                    This document cannot be modified
                    after Script 1 runs.
date:               2026-03-05
author:             Eric Robert Lawson
                    OrganismCore
protocol_version:   Workflow_Protocol.md v2.0
status:             LOCKED

precursor_documents:
  BRCA-S8c-PLAIN  (Script 1 Plain Account — ratio derived)
  BRCA-S8e-PLAIN  (Script 2 Plain Update)
  BRCA-S8g        (Script 3 Results and Reasoning)
  BRCA-S8h        (Cross-Subtype Literature Check)
  foxa1_ezh2_ratio_reasoning_artifact.md
                  (Full reasoning artifact for the ratio)

derivation_source:
  GSE176078 (Wu et al. 2021, Nature Genetics)
  Single-cell RNA-seq, 100,064 cells, 26 patients
  Cross-subtype Script 1 — BRCA-S8c
  The ratio was derived from attractor geometry
  applied to 19,542 cancer cells across five subtypes.
  Predictions were locked before the script ran.
  The ordering was confirmed exactly.

ratio_values_to_be_validated:
  LumA:  9.38
  LumB:  8.10
  HER2:  3.34
  TNBC:  0.52
  CL:    0.10
  Order: LumA > LumB > HER2 > TNBC > CL
```

---

## PREAMBLE

```
The FOXA1/EZH2 ratio was derived from single-cell
RNA-seq geometry in the cross-subtype analysis.
It correctly ordered all five measurable breast
cancer subtypes in the predicted sequence across
19,542 individual cancer cells.

That is one dataset. One technology. One derivation.

Before this ratio can be proposed as a diagnostic
instrument, it must be tested in independent
datasets using independent measurement technologies.

This document locks the predictions before any
validation script runs.

Three datasets are available for validation without
any new data collection:

  DATASET 1 — GSE176078 (already in pipeline):
    Same scRNA-seq data the ratio was derived from.
    New analysis: per-patient ratio, formally scored.
    Purpose: formal concordance, not re-derivation.
    What this adds: per-patient classification vs
    PAM50. Cohen's kappa. Sensitivity/specificity.

  DATASET 2 — TCGA-BRCA RPPA (new download):
    ~800 breast cancer patients.
    Actual protein quantification (RPPA).
    Not RNA. Protein.
    PAM50 labels available and matched.
    FOXA1 protein: confirmed in TCPA panel.
    EZH2 protein: confirmed in TCPA panel.
    Purpose: first protein-level test of the ratio.
    This is the critical tier.

  DATASET 3 — METABRIC (already in pipeline):
    1,980 patients, mRNA (microarray z-scores).
    FOXA1 and EZH2 both confirmed present
    in prior scripts (Claudin-low Script 4 log).
    PAM50 equivalent: CLAUDIN_SUBTYPE column.
    Purpose: mRNA-level large-cohort confirmation.
    n=1,980 — the largest available public cohort.

The three datasets span:
  — Two measurement technologies
    (RNA: scRNA-seq + bulk microarray;
     Protein: RPPA)
  — Three independent patient cohorts
  — n=19,542 cells + n=~800 + n=1,980

Together they constitute the strongest public-data
validation package available before any prospective
IHC study.
```

---

## PART I — WHAT IS BEING VALIDATED

### I.1 — The primary claim

```
PRIMARY CLAIM:
  The FOXA1/EZH2 ratio, when computed from gene or
  protein expression data for a breast cancer sample,
  orders the five major breast cancer subtypes
  (LumA, LumB, HER2, TNBC/Basal, Claudin-low) in
  the sequence:

    LumA > LumB > HER2 > TNBC > CL

  with statistically significant separation between
  groups, across independent datasets and measurement
  technologies.

This is the claim being tested.
Nothing more. Nothing less.
```

### I.2 — The secondary claims

```
SECONDARY CLAIM 1:
  The ratio, applied as a continuous variable with
  the derived cut-points (>8 / ~3 / ~0.5 / <0.2),
  correctly classifies individual patient samples
  into their PAM50-assigned subtype at a concordance
  level (Cohen's kappa) exceeding chance.

SECONDARY CLAIM 2:
  In TCGA-BRCA with survival data, patients stratified
  by FOXA1/EZH2 ratio show survival differences
  consistent with subtype-stratified survival.
  Ratio-guided stratification should be at least
  as informative as PAM50 for survival separation.

SECONDARY CLAIM 3:
  The ratio ordering is more robust than either
  FOXA1 alone or EZH2 alone as a subtype classifier.
  The ratio provides additional information beyond
  either individual protein.
```

### I.3 — What would constitute failure

```
FAILURE CRITERION 1 (primary):
  If the ratio in TCGA RPPA protein data does NOT
  order subtypes in the sequence LumA > LumB > HER2
  > TNBC — that is, if two or more adjacent subtypes
  are inverted — the primary claim fails in protein
  data.

  Note: CL vs TNBC ordering is EXPECTED to be
  sensitive to measurement technology due to the
  EZH2 confound established in BRCA-S8e-PLAIN.
  A CL/TNBC inversion in protein data would be
  noted and explained, not scored as a primary
  failure.

FAILURE CRITERION 2 (secondary):
  If Cohen's kappa for ratio-based classification
  vs PAM50 is < 0.20 (slight agreement or less),
  the secondary concordance claim fails.

FAILURE CRITERION 3 (survival):
  If ratio quartiles show no survival difference
  (logrank p > 0.20) in TCGA-BRCA with OS data,
  the survival prediction fails.

BIOLOGICAL FAILURE (unfalsifiable without this test):
  If EZH2 protein in TCGA RPPA is NOT highest in
  TNBC/Basal, this would directly contradict the
  graded elevation claim confirmed in BRCA-S8h
  (CS-LIT-4). This would be a significant finding
  requiring explanation.
```

---

## PART II — LOCKED PREDICTIONS

### II.1 — Dataset 1: GSE176078 per-patient formal concordance

```
PREDICTION R1-A:
  When the FOXA1/EZH2 ratio is computed per patient
  (median FOXA1 mRNA / median EZH2 mRNA across all
  cancer cells for each patient), the per-patient
  ratio orders the five subtypes in the sequence
  LumA > LumB > HER2 > TNBC > CL.

  Kruskal-Wallis p < 0.001 across groups.
  All pairwise comparisons LumA>LumB, LumB>HER2,
  HER2>TNBC: p < 0.05.

PREDICTION R1-B:
  Using the derived cut-points applied to per-patient
  ratios, Cohen's kappa for ratio-based vs PAM50
  classification ≥ 0.40 (moderate agreement or
  better) in the GSE176078 cohort.

  Rationale: The ratio was derived from this data,
  so concordance should be high. Kappa < 0.40 would
  indicate the ratio ordering does not translate
  from per-cell to per-patient level cleanly.

PREDICTION R1-C:
  The ratio outperforms either FOXA1 alone or EZH2
  alone as a subtype classifier (higher kappa for
  ratio than for either single gene).
  The ratio is not reducible to either component.
```

### II.2 — Dataset 2: TCGA RPPA protein-level validation

```
PREDICTION R2-A (PRIMARY):
  The FOXA1/EZH2 protein ratio in TCGA RPPA data
  orders the PAM50 subtypes in the sequence:
  LumA > LumB > HER2 > Basal

  With Kruskal-Wallis p < 0.001 across the four
  major subtypes (LumA, LumB, HER2, Basal).

  With pairwise Mann-Whitney:
    LumA vs LumB: p < 0.05
    LumB vs HER2: p < 0.05
    HER2 vs Basal: p < 0.05

  This is the most important prediction in this
  document. If protein-level data confirms the
  RNA-derived ordering, the ratio has crossed from
  RNA observation to protein-level validation.

PREDICTION R2-B:
  EZH2 protein in TCGA RPPA is highest in
  Basal/TNBC, followed by HER2, then LumB, then
  LumA — in that order.
  This replicates CS-LIT-4 (CONFIRMED) at the
  protein level in an independent cohort.

PREDICTION R2-C:
  FOXA1 protein in TCGA RPPA is highest in LumA,
  followed by LumB, then HER2, then Basal — in that
  order.
  This is the expected protein-level gradient
  consistent with FOXA1's role as a luminal
  identity marker.

PREDICTION R2-D (survival):
  TCGA-BRCA patients in the upper quartile of
  FOXA1/EZH2 protein ratio have better overall
  survival than patients in the lower quartile.
  Logrank p < 0.05.
  This is a test of whether the ratio encodes
  prognostic information independent of subtype.

PREDICTION R2-E (concordance):
  Cohen's kappa for protein ratio-based
  classification vs PAM50 ≥ 0.30 (fair agreement
  or better).

  Rationale: RPPA protein does not map identically
  to mRNA. Some compression is expected. Kappa ≥ 0.30
  in protein data would constitute meaningful
  validation. Kappa < 0.20 would indicate the
  protein ratio does not reliably classify subtypes.

PREDICTION R2-F (ratio vs single protein):
  The FOXA1/EZH2 ratio classifies PAM50 subtypes
  more accurately (higher kappa) than either FOXA1
  protein alone or EZH2 protein alone.
  The ratio is not reducible to either component
  at the protein level.
```

### II.3 — Dataset 3: METABRIC mRNA large-cohort confirmation

```
PREDICTION R3-A:
  The FOXA1/EZH2 mRNA ratio in METABRIC (n=1,980)
  orders the CLAUDIN_SUBTYPE groups in the sequence:
  LumA > LumB > HER2 > Basal/TNBC

  Kruskal-Wallis p < 0.001 across groups.
  All pairwise comparisons between adjacent subtypes
  significant at p < 0.01 given the large n.

PREDICTION R3-B:
  The LumA/LumB separation is statistically
  significant (Mann-Whitney p < 0.001) in METABRIC.
  This is the therapeutically critical cut — LumA
  vs LumB is the clinically most important subtype
  distinction for the ratio's treatment stratification
  logic (kinase lock vs chromatin lock).

PREDICTION R3-C:
  Median FOXA1/EZH2 ratio in METABRIC is higher in
  OS-survivors at 10 years than in non-survivors.
  Patients in the upper quartile of ratio have
  better 10-year OS than patients in the lower
  quartile.
  Logrank p < 0.05.
  (METABRIC has genuine 10-year follow-up data.
  This tests whether the ratio encodes long-window
  survival information in a large cohort.)

PREDICTION R3-D:
  The mRNA ratio in METABRIC provides the same
  directional ordering as the protein ratio in
  TCGA RPPA and the scRNA-seq ratio in GSE176078.
  Three datasets, three technologies, same ordering.
  If this prediction is confirmed, the ordering
  principle is robust to measurement technology.
```

---

## PART III — THE EZH2 CONFOUND PREDICTION

```
The EZH2-free PCA analysis (BRCA-S8e-PLAIN,
confirmed in BRCA-S8h CS-LIT-3) established:

  When EZH2 is included in the measurement,
  TNBC appears deeper than CL.
  When EZH2 is removed, CL is correctly deeper.

  Reason: In TNBC, EZH2 is the depth mechanism.
  EZH2 elevation and identity loss are the same
  biology measured from two angles.
  In CL, EZH2 is not the depth mechanism.
  CL's EZH2 is only moderately elevated (+67%).

For the FOXA1/EZH2 ratio specifically:

PREDICTION R4-A:
  The CL/TNBC ordering in the ratio will be
  dataset-dependent and technology-dependent.

  In scRNA-seq (GSE176078): CL ratio < TNBC ratio.
  (Already confirmed: 0.10 < 0.52)
  This is the expected direction — lower ratio means
  less luminal identity and/or more EZH2 suppression.
  CL has both (very low FOXA1, moderate EZH2).
  TNBC has high EZH2 specifically.
  The ratio correctly places CL below TNBC here.

  In RPPA protein data (TCGA):
  The CL/TNBC comparison in protein ratio will
  depend on whether the RPPA captures FOXA1's
  near-total absence in CL (-97.8%) more strongly
  than EZH2's moderate elevation in CL (+67%).
  If FOXA1 protein absence dominates: CL ratio
  will be very low, correctly below TNBC.
  If EZH2 RPPA signal compresses CL more than
  TNBC: the ordering may approach equivalence.

  PREDICTION: CL ratio ≤ TNBC ratio in all three
  datasets. The FOXA1 near-total absence (-97.8%)
  dominates the ratio regardless of EZH2 level.

PREDICTION R4-B:
  The LumA/LumB separation is the most robust
  and reproducible pairwise comparison across all
  three datasets and both technologies.
  This is the therapeutically most important
  comparison and should show the least technology-
  dependent variation because:
    — Both subtypes have expressed FOXA1
    — Both subtypes have moderate EZH2
    — The ratio difference is driven by chromatin
      lock effects on relative expression
    — This is not dominated by near-total absence
      of either protein (as in TNBC/CL)
```

---

## PART IV — THE SCRIPT DESIGN

### IV.1 — What Script 1 will do

```
Script 1 runs three components in sequence.
Each component is independent.
Results from one component do not influence
the analysis of another.

COMPONENT A — GSE176078 PER-PATIENT CONCORDANCE:

  Input:
    GSE176078 scRNA-seq (already downloaded —
    same data used in cross-subtype Script 1)
    Patient PAM50 labels from sample metadata

  Analysis:
    1. For each patient: compute median FOXA1 mRNA
       and median EZH2 mRNA across all cancer cells
    2. Compute per-patient ratio = FOXA1/EZH2
    3. Test ordering across PAM50 groups
       (Kruskal-Wallis + pairwise Mann-Whitney)
    4. Apply cut-points to classify each patient:
       >8.0  → LumA predicted
       5.0–8.0 → LumB predicted
       1.0–5.0 → HER2 predicted
       0.2–1.0 → TNBC predicted
       <0.2  → CL predicted
    5. Compute Cohen's kappa vs PAM50 label
    6. Compute sensitivity/specificity per subtype
    7. Compare ratio kappa vs FOXA1-alone kappa
       and EZH2-alone kappa

  Note on cut-points:
    The cut-points above are the first proposed
    clinical cut-points derived from the ratio values.
    They are centred on the midpoints between
    adjacent subtype values:
      LumA=9.38 / LumB=8.10 → cut at 8.0
      LumB=8.10 / HER2=3.34 → cut at 5.0
      HER2=3.34 / TNBC=0.52 → cut at 1.0
      TNBC=0.52 / CL=0.10   → cut at 0.2
    These cut-points are LOCKED before the script
    runs. They cannot be adjusted after seeing
    the concordance results.
    If the concordance is poor, the cut-points
    are recorded as suboptimal and a revised set
    is computed post-hoc (clearly labelled as
    data-derived, not pre-specified).

─────────────────────────────────────────────────

COMPONENT B — TCGA RPPA PROTEIN VALIDATION:

  Input:
    TCGA-BRCA RPPA protein data
    Source: TCPA portal (tcpaportal.org)
    File: BRCA_RPPA_data.csv or equivalent
    FOXA1 antibody: confirmed in TCPA BRCA panel
    EZH2 antibody: confirmed in TCPA BRCA panel

    TCGA-BRCA clinical data with PAM50 labels
    Source: cBioPortal API
    Study: brca_tcga_pan_can_atlas_2018
    PAM50 column: SUBTYPE or equivalent

    TCGA-BRCA survival data (OS_MONTHS, OS_STATUS)
    Source: same cBioPortal clinical download

  Data acquisition steps:
    1. Download RPPA data from TCPA portal
       URL: https://www.tcpaportal.org/tcpa/download.html
       Select: BRCA cohort, Level 4 normalized data
    2. Download PAM50 + survival via cBioPortal API
       (same API used in all prior scripts)
    3. Match samples by TCGA barcode (first 15 chars)
    4. Extract FOXA1 and EZH2 rows from RPPA matrix
    5. Compute ratio per sample

  Analysis:
    1. Distribution of FOXA1, EZH2, and ratio
       by PAM50 subtype (violin plots + boxplots)
    2. Kruskal-Wallis test across subtypes
    3. Pairwise Mann-Whitney with Bonferroni correction
    4. Cohen's kappa for ratio-based classification
    5. KM curves by ratio quartile (OS endpoint)
    6. Logrank p for Q1 vs Q4 survival comparison
    7. Comparison: ratio kappa vs FOXA1-alone kappa
       vs EZH2-alone kappa
    8. Spearman r between RPPA ratio and RNA ratio
       in the overlapping TCGA samples (if any
       TCGA samples also appear in scRNA-seq data
       — unlikely but checked)

─────────────────────────────────────────────────

COMPONENT C — METABRIC mRNA LARGE COHORT:

  Input:
    METABRIC expression data (already downloaded)
    brca_metabric_mrna_median_all_sample_Zscores
    n=1,980 samples with FOXA1 and EZH2 confirmed
    present (verified in Claudin-low Script 4 log)

    METABRIC clinical data (already downloaded)
    CLAUDIN_SUBTYPE column for subtype labels
    OS_MONTHS and OS_STATUS for survival

  Analysis:
    1. Distribution of FOXA1, EZH2, and ratio
       by CLAUDIN_SUBTYPE (violin + box)
    2. Kruskal-Wallis + pairwise Mann-Whitney
    3. Cohen's kappa for ratio-based classification
       vs CLAUDIN_SUBTYPE
    4. KM curves by ratio quartile (OS endpoint)
    5. Logrank p Q1 vs Q4
    6. 10-year OS comparison: upper vs lower
       ratio quartile
    7. Comparison: ratio kappa vs each single gene

─────────────────────────────────────────────────

OUTPUTS:

  For each component:
    — Ordered ratio values per subtype (table)
    — Violin/boxplot of ratio by subtype
    — Pairwise p-value matrix
    — Cohen's kappa table (ratio, FOXA1-alone,
      EZH2-alone)
    — KM curve by quartile
    — Summary verdict for each prediction

  Combined summary:
    — Cross-dataset ordering comparison table
    — Technology-dependence assessment
    — Concordance comparison table
    — Final verdict: how many of the locked
      predictions were confirmed
```

### IV.2 — Technical notes for Script 1

```
Z-SCORE HANDLING:
  METABRIC data is already z-scored (Median
  All Sample Z-Scores). The ratio of two z-scores
  is a valid relative measure but is NOT directly
  comparable to the ratio of raw or log2 expression
  values from scRNA-seq. The ordering test (is
  LumA > LumB > HER2 > TNBC?) is valid regardless
  of the absolute values. The specific numeric
  values (9.38, 8.10, etc.) are scRNA-seq derived
  and will NOT match METABRIC z-score ratios.
  The ordering is what is being validated. The
  absolute values are the scRNA-seq observation.

RPPA NORMALIZATION:
  TCPA Level 4 data is median-centered and
  normalized. The ratio of two RPPA measurements
  is valid as a relative quantification. The
  numeric values will differ from scRNA-seq.
  Again: ordering is the test.

CUT-POINT APPLICATION:
  The pre-specified cut-points (>8.0, 5.0-8.0,
  1.0-5.0, 0.2-1.0, <0.2) were derived from
  scRNA-seq values. They will not be optimal for
  RPPA or METABRIC z-score units.
  For the concordance analysis:
    — Apply the pre-specified cut-points first
      (locked, reported as-is).
    — Compute optimal cut-points post-hoc using
      ROC analysis for each dataset (clearly
      labelled as data-derived).
    — Report both. The pre-specified cut-points
      show robustness. The optimised cut-points
      show what is achievable in each technology.

FOXA1 IN TCPA RPPA:
  The TCPA BRCA RPPA panel was confirmed to
  include both FOXA1 and EZH2 in the search
  conducted 2026-03-05. If either protein is
  absent from the downloaded panel, this must
  be reported immediately and Component B must
  be scored as NOT TESTABLE — not as FAILED.
  A missing protein is a data availability issue,
  not a biological falsification.

SAMPLE SIZE NOTE:
  TCGA RPPA: approximately 800 samples with both
  RPPA and PAM50 data after matching.
  This provides adequate power for subtype
  ordering tests and moderate power for
  concordance and survival analysis.
  Subtype-specific n: LumA ~300, LumB ~100,
  HER2 ~50, Basal ~100, Normal ~50.
  HER2 is the smallest subtype. Power for
  pairwise LumB vs HER2 comparison will be
  the most limited.
```

---

## PART V — WHAT THE THREE TIERS PROVE IF CONFIRMED

```
IF COMPONENT A CONFIRMS (GSE176078 per-patient):
  The ratio ordering translates from the per-cell
  level (19,542 cells) to the per-patient level
  (26 patients classified correctly).
  The cut-points are applicable at patient level.
  Formal concordance vs PAM50 is established
  in the derivation dataset.

IF COMPONENT B CONFIRMS (TCGA RPPA protein):
  The ratio ordering holds at the PROTEIN level.
  This is the critical tier.
  RNA-derived ratios that fail at the protein level
  are not translatable to IHC.
  RNA-derived ratios that hold at the protein level
  provide strong grounds for IHC translation.
  A confirmed protein-level ordering in ~800 patients
  is the most important single result this script
  can produce.

IF COMPONENT C CONFIRMS (METABRIC mRNA):
  The ratio ordering holds in the largest available
  public cohort (n=1,980) using a completely
  independent measurement platform (Illumina
  microarray, different from both scRNA-seq and RPPA).
  Technology-independence is established at the
  mRNA level.

IF ALL THREE CONFIRM:
  The ratio ordering principle is validated across:
    — Three independent patient cohorts
    — Two measurement technologies (RNA and protein)
    — n = 19,542 cells + ~800 patients + 1,980 patients
    — scRNA-seq + RPPA + bulk microarray

  This is the public-data validation package that
  justifies a prospective IHC study.
  This is the evidence base that makes the
  FOXA1/EZH2 ratio a publishable finding.
  This is the foundation for a clinical team to
  take the prospective concordance study seriously.

IF COMPONENT B FAILS (protein level):
  The ratio is an RNA-level observation that does
  not translate to protein.
  The IHC diagnostic tool as proposed is not
  supported.
  The finding is downgraded from a potential
  diagnostic instrument to a biological observation
  about transcript levels.
  This would be recorded honestly and the reasons
  investigated (protein/RNA discordance in
  specific subtypes, RPPA technical limitations,
  FOXA1 antibody specificity in the TCPA panel).
  A protein-level failure does NOT eliminate
  the finding — it redefines its scope.
```

---

## PART VI — DOCUMENT SERIES STATUS

```
FOXA1/EZH2 RATIO VALIDATION SERIES:

  RATIO-S1a:  before_script1_ratio_validation.md
              COMPLETE — LOCKED [THIS DOCUMENT]

  RATIO-S1b:  Script 1 (three-component validation)
              TO BE WRITTEN AND RUN

  RATIO-S1c:  script1_results_and_reasoning.md
              AFTER Script 1 runs

  RATIO-S2a:  before_script2.md
              Survival and cut-point optimisation
              (if Component B confirms)

  RATIO-S2b:  Script 2
              Extended survival analysis,
              multivariate Cox, subtype-specific
              cut-point refinement

  RATIO-S2c:  script2_results_and_reasoning.md

  RATIO-S3a:  literature_check.md
              After both scripts complete

  RATIO-FINAL: ratio_clinical_proposal.md
              Trial design and clinical proposal
              based on validated results
```

---

## STATUS BLOCK

```
document:           RATIO-S1a
type:               Before-Document — Script 1
status:             COMPLETE AND LOCKED
date:               2026-03-05
author:             Eric Robert Lawson / OrganismCore

predictions_locked:
  R1-A, R1-B, R1-C  — GSE176078 per-patient
  R2-A, R2-B, R2-C, R2-D, R2-E, R2-F — TCGA RPPA
  R3-A, R3-B, R3-C, R3-D — METABRIC mRNA
  R4-A, R4-B — EZH2 confound predictions

primary_prediction:
  R2-A — FOXA1/EZH2 protein ratio in TCGA RPPA
  orders LumA > LumB > HER2 > Basal at
  Kruskal-Wallis p < 0.001.
  This is the most important single test.

critical_failure_definition:
  Protein ratio (TCGA RPPA) does NOT order
  subtypes in correct sequence — two or more
  adjacent subtypes inverted.

data_status:
  GSE176078:  already downloaded, in pipeline
  METABRIC:   already downloaded, in pipeline
  TCGA RPPA:  requires download from TCPA portal
              URL: tcpaportal.org/tcpa/download.html
              File: BRCA Level 4 RPPA data
  TCGA PAM50: requires cBioPortal API call
              Study: brca_tcga_pan_can_atlas_2018

repository:         https://github.com/Eric-Robert-Lawson/
                    attractor-oncology
orcid:              https://orcid.org/0009-0002-0414-6544
contact:            OrganismCore@proton.me
```

---

*"The ratio was derived from geometry.*
*The predictions are locked.*
*The data will say what it says.*
*That is the only way to do this."*

— Eric Robert Lawson, March 5, 2026
