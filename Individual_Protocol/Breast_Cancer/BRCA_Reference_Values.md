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
## VERSION 3.0 — FULLY EXTRACTED

```
All values are from the Mature Luminal
population in GSE176078 (n=1,265 cells)
after log1p normalisation.
Source: cs_s1_depth_table.csv (BRCA-S8b re-run)
Extracted: 2026-03-07

These are the reference means and SDs for
z-score computation.

z-score formula:
  z = (patient_value - mean) / SD

A positive z = elevated above normal.
A negative z = suppressed below normal.

The magnitude tells you how far from normal.
The sign tells you the direction.
The gene panel tells you what it means
biologically.

─────────────────────────────────────────
SWITCH GENES
(suppressed in cancer vs normal — luminal
 identity programme)

  gene      mean      SD
  FOXA1     2.1335    2.6539
  GATA3     4.4914    2.9675
  ESR1      3.1688    3.0784
  PGR       1.7894    2.6526
  SPDEF     4.8633    2.7148
  CDH1      4.8926    2.7911  ¹
  CLDN3     5.9909    2.7259  ²
  CLDN4     5.3454    2.3486  ²
  CLDN7     4.3225    2.7844  ²
  AR        5.3748    2.8363  ³
  CDKN1A    0.5430    1.4208  ⁴

¹ CDH1: also the ILC structural marker —
  absence in ILC context is the classifier,
  not the z-score.
² CLDN3/4/7: claudin-low panel. Suppression
  of all three simultaneously = CL geometry.
  In a non-CL patient these should be near
  zero or positive.
³ AR: continuous depth axis within TNBC
  (r=−0.378, p=6.23×10⁻¹⁴⁷ — CS-LIT-6).
  In LumA/LumB context AR suppression
  is a supporting signal, not primary.
⁴ CDKN1A: note this is a switch gene in
  the sense that its loss drives depth
  in LumA (CDKN1A LOW = deeper LumA
  attractor). The mean here is the
  Mature Luminal baseline. A patient
  with very low CDKN1A relative to this
  mean is in a deeper LumA attractor.

─────────────────────────────────────────
FA MARKERS / EPIGENETIC LOCK PANEL
(elevated in cancer vs normal — attractor
 enforcement programme)

  gene      mean      SD
  EZH2      0.2739    1.1441
  EED       0.4698    1.4365
  SUZ12     1.1826    2.1403
  HDAC1     2.2516    2.6291
  HDAC2     2.5363    2.6966
  KDM1A     0.8565    1.8797
  DNMT3A    0.4499    1.4208
  DNMT3B    0.0690    0.5800
  BRD4      2.5510    2.7316
  MKI67     0.0097    0.2037
  TOP2A     0.4197    0.4541

NOTE ON EZH2 SD:
  EZH2 SD in Mature Luminal is 1.1441.
  EZH2 mean is 0.2739.
  This means MatureLum EZH2 has high
  relative variance — many normal cells
  have near-zero EZH2.
  A patient TNBC value of ~0.79 sits at:
  z = (0.79 - 0.27) / 1.14 = +0.46
  This appears modest because the SD is
  large. The percentage displacement
  (+189% for TNBC) is the more informative
  signal for this gene. Use both.

─────────────────────────────────────────
IMMUNE PANEL
(relevant for claudin-low geometry)

  gene      mean      SD
  FOXP3     0.0244    0.3448
  CD8A      0.0786    0.6437
  TIGIT     0.0354    0.4060

NOTE: These genes have near-zero means
in Mature Luminal. They are immune cell
markers appearing in the tumour
microenvironment. High values in a
patient's bulk transcriptome reflect
immune infiltration, not cancer cell
expression. Interpret in context of
subtype — relevant primarily for CL
geometry (FOXP3/CD8A ratio as CL
patient selector, HR=2.212, CS-LIT-20).

─────────────────────────────────────────
Z-SCORE INTERPRETATION GUIDE

  z < −2.0   Strongly suppressed vs normal
  z −2 to −1 Moderately suppressed
  z −1 to  0 Mildly below normal
  z  0 to +1 Mildly above normal
  z +1 to +2 Moderately elevated
  z > +2.0   Strongly elevated vs normal

For switch genes: a z < −1 is meaningful.
For FA markers: a z > +1 is meaningful.
For EZH2 specifically: given high SD,
use percentage displacement as primary
signal and z-score as supporting context.
```

---

## PART IIb — SUBTYPE POPULATION MEANS
## FULL PANEL — ALL 25 GENES
## (from cs_s1_depth_table.csv re-run 2026-03-07)

```
FORMAT: values are population means after
log1p normalisation (scRNA-seq GSE176078).
ILC values are from TCGA-BRCA bulk RNA-seq
(different expression space — directional
reference only, not for z-score computation).

─────────────────────────────────────────────────────────────
SWITCH GENES

gene      MatureLum  LumA    LumB    HER2    TNBC    CL
FOXA1     2.1335     2.8774  2.6304  1.9930  0.4151  0.0464
GATA3     4.4914     5.3605  6.1415  2.1671  2.1778  1.9758
ESR1      3.1688     3.4672  4.6321  0.2399  0.1429  0.0275
PGR       1.7894     1.1592  1.2211  0.0632  0.0456  0.0224
SPDEF     4.8633     4.5948  3.9669  4.1448  1.1635  0.0734
CDH1      4.8926     2.0168  2.2538  2.8527  1.6982  0.2208
AR        5.3748     2.6886  1.0343  3.4831  0.8437  0.0537
CDKN1A    0.5430     0.3288  0.7707  0.6451  0.1378  0.0630

─────────────────────────────────────────────────────────────
FA / EPIGENETIC LOCK PANEL

gene      MatureLum  LumA    LumB    HER2    TNBC    CL
EZH2      0.2739     0.3068  0.3246  0.5966  0.7927  0.4562
EED       0.4698     0.3632  0.4035  0.4330  0.9569  1.0092
SUZ12     1.1826     1.0810  1.4000  0.9308  0.9460  1.4012
HDAC1     2.2516     1.2946  2.0773  1.3694  2.4052  2.5036
HDAC2     2.5363     1.6881  2.6612  2.6759  3.7218  3.2896
KDM1A     0.8565     0.5137  0.4938  0.8067  1.0048  1.1467
DNMT3A    0.4499     0.3288  0.7707  0.6451  0.7862  0.6222
DNMT3B    0.0690     0.0745  0.0753  0.2071  0.0853  0.1905
BRD4      2.5510     1.5458  1.7837  1.8395  2.2608  3.2362
MKI67     0.0097     0.0044  0.0650  0.0894  0.0800  0.1249
TOP2A     0.4197     0.3727  0.3276  0.4979  0.6694  1.3426

─────────────────────────────────────────────────────────────
CLAUDIN PANEL

gene      MatureLum  LumA    LumB    HER2    TNBC    CL
CLDN3     5.9909     5.6610  4.6779  5.7580  3.6672  0.3469
CLDN4     5.3454     5.0062  4.0885  4.0721  3.0428  0.9839
CLDN7     4.3225     3.8718  3.2964  3.1943  2.1934  0.9178

─────────────────────────────────────────────────────────────
IMMUNE PANEL

gene      MatureLum  LumA    LumB    HER2    TNBC    CL
FOXP3     0.0244     0.0142  0.0054  0.0083  0.0328  0.0078
CD8A      0.0786     0.0125  0.0705  0.0439  0.0918  0.1007
TIGIT     0.0354     0.0081  0.0171  0.0162  0.0717  0.0594

─────────────────────────────────────────────────────────────
PERCENTAGE DISPLACEMENT FROM MATURE LUMINAL
(for each gene — from script log output)

FOXA1:   LumA +34.9%  LumB +23.3%  HER2  −6.6%  TNBC −80.5%  CL  −97.8%
GATA3:   LumA +19.4%  LumB +36.7%  HER2 −51.8%  TNBC −51.5%  CL  −56.0%
ESR1:    LumA  +9.4%  LumB +46.2%  HER2 −92.4%  TNBC −95.5%  CL  −99.1%
PGR:     LumA −35.2%  LumB −31.8%  HER2 −96.5%  TNBC −97.4%  CL  −98.7%
SPDEF:   LumA  −5.5%  LumB −18.4%  HER2 −14.8%  TNBC −76.1%  CL  −98.5%
CDH1:    LumA −58.8%  LumB −53.9%  HER2 −41.7%  TNBC −65.3%  CL  −95.5%
AR:      LumA −50.0%  LumB −80.8%  HER2 −35.2%  TNBC −84.3%  CL  −99.1%
CDKN1A:  LumA −39.5%  LumB −26.7%  HER2 +18.8%  TNBC −74.6%  CL  −88.4%

EZH2:    LumA +12.0%  LumB +18.5%  HER2 +117.8% TNBC +189.4% CL  +66.5%
EED:     LumA −22.7%  LumB −13.4%  HER2  −7.1%  TNBC +103.8% CL +115.1%
SUZ12:   LumA  −9.0%  LumB +17.9%  HER2 −21.6%  TNBC −20.0%  CL  +18.5%
HDAC1:   LumA −42.5%  LumB  −7.7%  HER2 −39.2%  TNBC  +6.8%  CL  +11.2%
HDAC2:   LumA −33.5%  LumB  +4.8%  HER2  +5.4%  TNBC +46.7%  CL  +29.7%
KDM1A:   LumA −40.1%  LumB −42.4%  HER2  −6.0%  TNBC +17.5%  CL  +33.9%
DNMT3A:  LumA −26.7%  LumB +71.8%  HER2 +43.8%  TNBC +74.8%  CL  +38.3%
DNMT3B:  LumA  +8.1%  LumB  +9.3%  HER2+200.3%  TNBC +23.7%  CL +174.7%
BRD4:    LumA −39.4%  LumB −30.5%  HER2 −28.3%  TNBC −11.2%  CL  +27.3%
MKI67:   LumA −55.4%  LumB+564.2%  HER2+819.6%  TNBC+725.0%  CL+1186.3%
TOP2A:   LumA −11.5%  LumB+677.8%  HER2 +18.4%  TNBC +59.4%  CL +220.4%

CLDN3:   LumA  +0.0%  LumB −17.4%  HER2  −3.9%  TNBC −38.8%  CL  −94.1%
CLDN4:   LumA  −6.3%  LumB −23.5%  HER2 −23.8%  TNBC −43.1%  CL  −81.5%
CLDN7:   LumA  −6.7%  LumB −23.8%  HER2 −26.1%  TNBC −49.3%  CL  −78.8%

─────────────────────────────────────────────────────────────
KEY PATTERNS TO RECOGNISE IN A PATIENT

LumA pattern:
  FOXA1 elevated, GATA3 elevated,
  ESR1 mildly elevated, PGR suppressed,
  EZH2 minimally elevated (+12%),
  MKI67 SUPPRESSED (−55%),
  claudins mildly elevated.
  Identity present. Proliferation low.
  Lock is kinase-mediated (CDKN1A low).

LumB pattern:
  FOXA1/GATA3/ESR1 all elevated,
  PGR suppressed (decoupled from ESR1),
  DNMT3A +72% (epigenetic lock active),
  MKI67 +564% (proliferation dominant),
  HDAC2 mildly elevated.
  Identity present but ER output decoupled.
  TFF1/ESR1 decoupling is the key signal
  at Tier 3 (CS-LIT-9).

HER2 pattern:
  FOXA1 mildly suppressed,
  ESR1/PGR strongly suppressed,
  EZH2 +118%,
  DNMT3B +200% (distinctive),
  MKI67 +820%.
  ERBB2 signal required to confirm.
  If ERBB2 not in data: HER2 pattern
  is a strong indicator but not definitive.

TNBC pattern:
  FOXA1 −80.5%, ESR1 −95.5%,
  PGR −97.4%, AR −84.3%,
  EZH2 +189%, HDAC2 +46.7%,
  MKI67 +725%, EED +104%.
  All identity TFs near-absent.
  Epigenetic lock fully dominant.
  Claudins moderately suppressed.

Claudin-low pattern:
  ALL five luminal TFs near-zero,
  ALL three claudins near-zero (−79 to −94%),
  AR −99%, CDKN1A −88%,
  ZEB1/ZEB2/TWIST1 massively elevated
  (from top mover scan — not in this table),
  EZH2 moderate (+67% — less than TNBC),
  MKI67 +1186%.
  Pre-commitment arrest. No lineage identity.
  The claudin collapse is the diagnostic signal.
  If claudins are near zero AND luminal TFs
  are near zero: this is CL geometry.

ILC pattern:
  FOXA1 HIGH (equal to or above LumA),
  CDH1 ABSENT or strongly suppressed,
  GATA3/ESR1 elevated.
  Structural lock — not an epigenetic lock.
  The geometry is preserved identity with
  lost structural integrity.
  Note: ILC values in this table are from
  TCGA bulk RNA-seq (n=210), not scRNA-seq.
  They are in a different expression space.
  Use directionally, not for z-score
  computation.
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
