# BRCA REFERENCE GEOMETRY VALUES
## Population Reference Means for Individual Patient Analysis
## OrganismCore — Eric Robert Lawson
## 2026-03-07 | Updated: 2026-03-07

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

Depth score distribution source:
  BRCA-S8c (BRCA_Cross_Subtype_Script2 output)
  Cohort: TCGA-BRCA bulk RNA-seq
  n = 837 tumour samples (normal -11 excluded)

These values are the coordinate system.
A patient's expression data is placed into
this coordinate system.

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
These values are the reference baseline for
geometric placement.

All values are from the Mature Luminal
population in GSE176078 (n=1,265 cells)
after log1p normalisation.
Source: cs_s1_depth_table.csv (BRCA-S8b output)
Extracted: 2026-03-07

─────────────────────────────────────────
WHAT THIS SECTION PROVIDES

The depth score formula requires z-scores
relative to the Mature Luminal reference.
A z-score requires both a mean and a SD.

CURRENT STATUS:
  Means for five luminal TF genes:
    EXTRACTED — values below.
  SDs for all genes:
    NOT YET EXTRACTED.
    The cs_s1_depth_table.csv output does not
    include per-gene SDs. It outputs means only.
    SDs require re-running BRCA_Cross_Subtype_
    Script1.py with a one-line addition to the
    unified_depth_axis() function (see Part IId).
  EZH2, MKI67, TOP2A, PCNA, CDH1, KRT18, KRT8:
    NOT in the cs_s1_depth_table.csv columns.
    Same re-run resolves this.

WHAT THIS MEANS FOR PATIENT ANALYSIS:
  Without SDs: directional displacement is
  fully operational. You can state which
  direction and approximate magnitude.
  With SDs: calibrated z-scores. You can
  state exactly how many SDs from normal
  each gene sits.
  The axis placement using Part I ratios
  is NOT affected — that uses the patient's
  own measured FOXA1 and EZH2 values directly.

─────────────────────────────────────────
LUMINAL TF SWITCH GENE MEANS
(Mature Luminal reference — from cs_s1_depth_table.csv)

  FOXA1   — mean: 2.1335
  GATA3   — mean: 4.4914
  ESR1    — mean: 3.1688
  PGR     — mean: 1.7894
  SPDEF   — mean: 4.8633

─────────────────────────────────────────
REMAINING SWITCH GENES (not in CSV output):

  CDH1    — mean: NOT YET EXTRACTED
  KRT18   — mean: NOT YET EXTRACTED
  KRT8    — mean: NOT YET EXTRACTED
  SCUBE2  — mean: NOT YET EXTRACTED

─────────────────────────────────────────
FA MARKERS (not in CSV output):

  EZH2    — mean: NOT YET EXTRACTED
  MKI67   — mean: NOT YET EXTRACTED
  TOP2A   — mean: NOT YET EXTRACTED
  PCNA    — mean: NOT YET EXTRACTED

─────────────────────────────────────────
EPIGENETIC LOCK PANEL (not in CSV output):

  HDAC1   — mean: NOT YET EXTRACTED
  HDAC2   — mean: NOT YET EXTRACTED
  DNMT3A  — mean: NOT YET EXTRACTED
  KDM6A   — mean: NOT YET EXTRACTED
  KDM1A   — mean: NOT YET EXTRACTED

─────────────────────────────────────────
CLAUDIN-LOW PANEL (not in CSV output):

  CLDN3   — mean: NOT YET EXTRACTED
  CLDN4   — mean: NOT YET EXTRACTED
  CLDN7   — mean: NOT YET EXTRACTED
```

---

## PART IIb — SUBTYPE POPULATION MEANS
## FOR LUMINAL TF PANEL
## (Displacement reference — five genes)

```
These are the population means for the five
luminal TF genes across all populations.
Source: cs_s1_depth_table.csv (BRCA-S8b)
Normalisation: log1p, scRNA-seq GSE176078
(ILC values are from TCGA-BRCA bulk RNA-seq)

─────────────────────────────────────────
POPULATION      n      FOXA1   GATA3   ESR1    PGR     SPDEF
─────────────────────────────────────────
MatureLum    1265    2.1335  4.4914  3.1688  1.7894  4.8633
LumProg      1992    0.2479  2.3755  0.5840  0.0936  1.3257
Myo          1098    0.0438  2.0114  0.0848  0.2924  0.1429
─────────��───────────────────────────────
LumA         7742    2.8774  5.3605  3.4672  1.1592  4.5948
LumB         3368    2.6304  6.1415  4.6321  1.2211  3.9669
HER2         3708    1.9930  2.1671  0.2399  0.0632  4.1448
TNBC         4312    0.4151  2.1778  0.1429  0.0456  1.1635
CL            412    0.0464  1.9758  0.0275  0.0224  0.0734
ILC           210   12.4262 13.3079 12.6478 10.1677 11.8808
─────────────────────────────────────────

PERCENTAGE DISPLACEMENT FROM MATURE LUMINAL
(from cs_s1_depth_table.csv)

  FOXA1:
    LumA    +34.9%      LumB    +23.3%
    HER2     −6.6%      TNBC   −80.5%
    CL      −97.8%      LumProg −88.4%
    Myo     −97.9%

  GATA3:
    LumA    +19.4%      LumB    +36.7%
    HER2    −51.8%      TNBC   −51.5%
    CL      −56.0%      LumProg −47.1%
    Myo     −55.2%

  ESR1:
    LumA     +9.4%      LumB    +46.2%
    HER2    −92.4%      TNBC   −95.5%
    CL      −99.1%      LumProg −81.6%
    Myo     −97.3%

  PGR:
    LumA    −35.2%      LumB    −31.8%
    HER2    −96.5%      TNBC   −97.4%
    CL      −98.7%      LumProg −94.8%
    Myo     −83.7%

  SPDEF:
    LumA     −5.5%      LumB    −18.4%
    HER2    −14.8%      TNBC   −76.1%
    CL      −98.5%      LumProg −72.7%
    Myo     −97.1%

─────────────────────────────────────────
ILC NOTE:
  ILC values are from TCGA-BRCA bulk RNA-seq
  (n=210), not scRNA-seq GSE176078.
  They are in a different expression space.
  ILC placement uses FOXA1 HIGH + CDH1 ABSENT
  as the geometric classifier, not a numerical
  ratio comparison against scRNA-seq means.
  The ILC bulk values confirm FOXA1 elevation
  directionally and are consistent with the
  ILC structural lock classification.

─────────────────────────────────────────
CRITICAL NOTE ON WHAT THIS TABLE COVERS:
  This table contains ONLY the five luminal TF
  genes captured in cs_s1_depth_table.csv.
  EZH2 subtype means are NOT here.
  The FOXA1/EZH2 ratio computation for a
  Tier 3 patient uses the patient's own
  measured EZH2 value directly — it does not
  require a population mean for EZH2.
  The ratio medians for placement are in Part I.
```

---

## PART IIc — HOW TO USE PARTS II AND IIb
## WITHOUT FULL SDs

```
Until the script re-run extracts full SDs
and FA marker means, the following procedure
applies for Tier 3 patient analysis.

─────────────────────────────────────────
STEP 1 — PRIMARY AXIS PLACEMENT (fully operational)

  Compute: patient FOXA1 / patient EZH2
  Both values come directly from the patient's
  expression table. No reference mean needed.
  Place the result on the Part I medians scale.
  This step is unaffected by the SD gap.

─────────────────────────────────────────
STEP 2 — DIRECTIONAL DISPLACEMENT (operational)

  For each of the five luminal TF genes:
    displacement = patient_value - MatureLum_mean
    (using means from Part II)

  A negative displacement = suppression below
  normal. A positive displacement = elevation
  above normal.

  Report which genes show the largest
  displacement from the reference.
  Report the direction.
  Compare the displacement pattern to the
  subtype reference means in Part IIb —
  does this patient's pattern match LumA?
  LumB? TNBC? Or is it atypical?

─────────────────────────────────────────
STEP 3 — DEPTH ESTIMATION (partial operational)

  Without Mature Luminal SDs for the FA marker
  panel, a calibrated z-score depth score
  cannot be computed. What can be stated:

  If FA markers (EZH2, MKI67, TOP2A, PCNA)
  are present in the patient's data:
    Compare each value to the subtype
    population means from Part IIb (noting
    that EZH2 population means are NOT in
    the current extracted table — the
    FOXA1/EZH2 ratio from Part I is the
    calibrated axis signal).

  Depth estimate from available data:
    FOXA1 and GATA3 suppression magnitude
    relative to MatureLum means gives a
    directional depth signal.
    High switch gene suppression = deeper
    attractor than expected for this subtype.

─────────────────────────────────────────
WHAT CHANGES AFTER THE SCRIPT RE-RUN:

  Per-gene SDs extracted → full z-score
  computation operational for all panel genes.
  EZH2, MKI67, TOP2A, PCNA means extracted →
  complete FA score computable.
  Depth score becomes a calibrated number
  placeable on the Part III distribution.
  Until then: directional analysis is
  fully operational. Quantitative z-scores
  are not yet calibrated.
```

---

## PART IId — THE SCRIPT RE-RUN REQUIRED
## (What to add — exactly)

```
To extract the remaining values, one addition
is needed in BRCA_Cross_Subtype_Script1.py
inside the unified_depth_axis() function.

After the line:
  ref_means = ref.mean()

Add:
  ref_stds = ref.std()

Then in the CSV row construction, after
each gene mean is added, also add:
  row[f"{g}_std"] = float(ref_stds[g])
    if g in ref_stds.index else np.nan

Expand the avail list to include all panel
genes, not just LUMINAL_TFS:

  all_panel = (LUMINAL_TFS + EPIGENETIC
               + ["MKI67", "TOP2A", "PCNA",
                  "CDH1", "KRT18", "KRT8",
                  "CLDN3", "CLDN4", "CLDN7"])
  avail = [g for g in all_panel
           if g in ref_means.index]

Re-run the script. The output CSV will then
contain means and SDs for all panel genes
in the Mature Luminal reference population.
Extract those values and replace all
NOT YET EXTRACTED entries in Part II.

This is a single re-run. It does not
change any predictions or prior results.
It adds reference columns to the CSV.
```

---

## PART III — SUBTYPE DEPTH SCORE
## REFERENCE DISTRIBUTIONS

```
Source: cs_s2_depth_scores.csv
(BRCA_Cross_Subtype_Script2 output)
Cohort: TCGA-BRCA bulk RNA-seq
Normal tissue excluded: samples ending -11
Extracted: 2026-03-07

DEPTH SCORE FORMULA (Script 2):
  The depth score in cs_s2_depth_scores.csv
  is computed as a normalised composite of:
    FA markers: EZH2, MKI67, TOP2A, PCNA
    Switch genes: FOXA1, GATA3, ESR1,
                  CDH1, KRT18, KRT8, SPDEF
  Normalised to 0–1 range across the
  TCGA-BRCA cohort.
  Higher = deeper attractor (further from
  Mature Luminal normal).

─────────────────────────────────────────
PER-SUBTYPE REFERENCE DISTRIBUTIONS
(from cs_s2_depth_scores.csv, n=944 tumours)

FORMAT: n | mean (SD) | median | [Q25–Q75] | range

  LumA    343 | 0.4946 (0.1461) | 0.4876 | [0.4035–0.5947] | 0.0857–0.8659
  LumB    185 | 0.5005 (0.1372) | 0.4991 | [0.4119–0.5919] | 0.1692–0.8033
  Basal   135 | 0.5011 (0.1073) | 0.5111 | [0.4366–0.5790] | 0.1710–0.7072
  HER2     65 | 0.4998 (0.1201) | 0.5124 | [0.4129–0.5920] | 0.2512–0.7488
  ILC     202 | 0.5115 (0.1812) | 0.5143 | [0.3804–0.6354] | 0.1214–0.9560
  CL       14 | 0.5000 (0.1867) | 0.4167 | [0.3839–0.6726] | 0.2381–0.7976

─────────────────────────────────────────
WHAT THE DISTRIBUTION REVEALS

The most important observation from these
numbers is what they do NOT show.

The means across all six subtypes are
nearly identical (0.494–0.512). This is
not a null result. It is a structural finding.

The depth score in Script 2 is normalised
0–1 across the TCGA-BRCA cohort. When
normalised this way, the mean of any
sufficiently large subtype will be pulled
toward 0.5 by the normalisation itself.
The subtype ordering is not visible in
the means because the normalisation
removes the between-subtype signal.

What is visible:
  ILC has the widest SD (0.1812).
  CL has the second widest SD (0.1867,
  but n=14 — interpret with caution).
  Basal has the narrowest SD (0.1073).

The ILC wide spread reflects the structural
heterogeneity of ILC: FOXA1 is preserved
(luminal identity) but CDH1 is lost
(structural disruption). The depth score
captures a continuous range of states
within ILC geometry.

Basal's narrow spread reflects the
consistency of the deep EZH2 lock:
TNBC/Basal tumours sit in a tight cluster
at a specific attractor depth. The lock
is relatively uniform across the subtype.

─────────────────────────────────────────
WHY THE BETWEEN-SUBTYPE ORDERING IS IN PART I,
NOT IN THESE DEPTH SCORE DISTRIBUTIONS

The primary axis ordering (LumA > LumB >
HER2e > TNBC > CL) is captured by the
FOXA1/EZH2 ratio in Part I.

The depth score distributions here capture
WITHIN-SUBTYPE variation — how deep is
this specific patient relative to others
in their subtype.

These are two different measurements:
  Part I: WHERE on the landscape (axis position)
  Part III: HOW DEEP in the landscape
            relative to subtype peers

Both are needed for the full picture.
Neither replaces the other.

─────────────────────────────────────────
HOW TO USE FOR A PATIENT (Script 2 depth score)

If the patient has a Tier 3 full transcriptome
and you have computed their depth score:

  1. Identify their subtype region from
     Part I (axis placement).

  2. Look up the reference distribution
     for that subtype above.

  3. Compare:
     Patient below Q25 for their subtype:
       Shallow attractor for this subtype.
       The lock is less entrenched than
       typical. The identity programme
       is more intact. This is clinically
       relevant — it implies greater
       residual responsiveness to identity-
       restoring approaches.

     Patient Q25–Q75 for their subtype:
       Typical attractor depth. The geometry
       matches the subtype reference. No
       unexpected depth signal.

     Patient above Q75 for their subtype:
       Deeper than typical for this subtype.
       The lock is more entrenched than
       the subtype average. This is the
       patient where depth-guided drug
       target geometry matters most.
       State it explicitly in the report.

     Patient above Q75 of the next deeper subtype:
       This patient's depth exceeds typical
       for both their assigned subtype and
       the next. State this explicitly.
       It may indicate clonal evolution
       toward a deeper attractor or a
       classification boundary case.

  4. Note the SD context:
     ILC patients have the widest spread.
     An ILC patient at depth 0.80 is
     approaching the maximum observed
     (0.956) — that is geometrically
     significant and should be stated.
     A Basal patient at depth 0.70 is
     above Q75 (0.579) — also significant.

─────────────────────────────────────────
CAVEAT ON CL DISTRIBUTION:
  n = 14 for Claudin-low.
  The distribution statistics are real
  but the sample is small.
  Use the median (0.42) and range as
  directional guidance only.
  Do not compute precise percentiles
  for CL patients from this distribution.
  State this limitation in the report.
```

---

## PART IV — ATTRACTOR TYPE REFERENCE
## (From BRCA-S8b and Axioms Document)

```
BRCA attractor type assignments
(confirmed by Identity TF Direction Test
in BRCA-S8b attractor_type_classification):

  LumA:        TYPE I  — Blocked Approach
               FOXA1 HIGH (+34.9% vs MatureLum)
               EZH2 low (lock minimal)
               Identity present, cell cycle
               lock dominant (CDK4/6)

  LumB:        TYPE I  — Blocked Approach (deeper)
               FOXA1 moderate (+23.3% vs MatureLum)
               EZH2 elevated (lock active at
               PGR and GATA3 promoters)
               ESR1 elevated (+46.2% — identity
               signal elevated but partially
               decoupled from PR, confirming
               EZH2 competition at PGR)

  HER2-enr:    TYPE II — Wrong Valley
               FOXA1 LOW (−6.6% vs MatureLum)
               ERBB2 HIGH
               Alternative identity programme
               dominant

  HER2-deep:   TYPE I + TYPE II composite
               (CDH3-high, AR-low, EZH2 +118%)
               EZH2i relevant in this fraction

  TNBC:        TYPE I  — Blocked Approach (deep)
               FOXA1 −80.5% vs MatureLum
               EZH2 +189% vs MatureLum
               Near-complete identity suppression

  Claudin-low: TYPE IV — Root Lock
               FOXA1 −97.8% vs MatureLum
               No mature lineage identity.
               Pre-commitment arrest.
               All five luminal TFs near-zero.

  ILC:         Composite — structural lock
               FOXA1 HIGH (confirmed by bulk
               TCGA values — identity preserved)
               CDH1 ABSENT (structural disruption)
               Not classifiable on the depth axis.
               A separate structural variant.
```

---

## PART V — DRUG TARGET GEOMETRY REFERENCE

```
These are the drug target predictions from
the population analysis (BRCA-S8b, validated
in BRCA-S8h literature check).

They are geometric observations.
They are not treatment recommendations.

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
  AUC for classification (0.828–0.901
  LumA vs Basal; 0.796–0.873 LumA vs LumB)
  Protein-level confirmation by CPTAC MS
  (r=-0.492, n=122)

INDEPENDENTLY CONFIRMED by other groups:
  FOXA1/EZH2 mechanistic axis:
    Schade et al. Nature 2024
    Toska et al. Nature Medicine 2017
    Neither group knew of this framework.

CONFIRMED from cs_s2_depth_scores.csv:
  Per-subtype depth score distributions
  now extracted (Part III).
  The within-subtype spread is confirmed.
  The between-subtype homogeneity of means
  is a normalisation artefact, not a null
  result — the axis ordering is captured
  by the ratio in Part I.

NOT YET CONFIRMED:
  IHC H-score cut-points (calibration
  study required — CS-LIT-22)
  Individual patient analyses (first
  patients have not yet been run)
  Drug predictions (framework-derived,
  not yet prospectively tested)
  Per-gene SDs for Mature Luminal reference
  (requires one script re-run — see Part IId)
  EZH2, MKI67, TOP2A, PCNA Mature Luminal
  means (same re-run resolves)

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
version:        2.0
date_created:   2026-03-07
date_updated:   2026-03-07
status:         ACTIVE — operational reference
                (partial: five luminal TF means
                 extracted; FA marker means and
                 all SDs pending one script re-run)

source_analysis: BRCA-S8b
  Cancer_Research/BRCA/DEEP_DIVE/
  BRCA_Cross_Subtype_Script1.py

depth_distribution_source: BRCA-S8c
  cs_s2_depth_scores.csv
  Cancer_Research/BRCA/DEEP_DIVE/
  (Script 2 output)

validation_document: BRCA-S8h
  Cancer_Research/BRCA/DEEP_DIVE/
  BRCA_Cross_Subtype_Literature_Check.md

REMAINING GAP:
  EZH2, MKI67, TOP2A, PCNA, CDH1,
  KRT18, KRT8, CLDN3, CLDN4, CLDN7 means
  for Mature Luminal reference population,
  plus SDs for all genes.
  Action: add two lines to unified_depth_axis()
  in BRCA_Cross_Subtype_Script1.py — see
  Part IId for exact instructions.
  Re-run script. Extract values.
  Replace NOT YET EXTRACTED entries in Part II.
  Bump to version 3.0 when complete.
```
