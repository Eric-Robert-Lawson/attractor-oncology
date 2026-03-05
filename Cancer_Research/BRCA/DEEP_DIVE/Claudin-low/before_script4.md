# CLAUDIN-LOW — SCRIPT 4 BEFORE-DOCUMENT
## External Dataset Validation — METABRIC + GSE96058
## Predictions Locked Before Script 4 Runs
## OrganismCore — Document BRCA-S7j
## Date: 2026-03-05

---

## DOCUMENT METADATA

```
document_id:        BRCA-S7j
series:             BRCA Deep Dive — Claudin-Low
folder:             Cancer_Research/BRCA/DEEP_DIVE/Claudin_Low/
type:               BEFORE-DOCUMENT (Script 4)
date:               2026-03-05
author:             Eric Robert Lawson
                    OrganismCore
protocol_version:   Workflow_Protocol.md v2.0
attractor_type:     TYPE 4 — ROOT LOCK
                    (ATTRACTOR_GEOMETRY_AXIOMS.md v2.0)
based_on:           BRCA-S7d (script2_results_and_reasoning.md)
                    BRCA-S7h (script3_results_and_reasoning.md)
                    BRCA-S7i (literature_check_2.md)
status:             LOCKED — predictions cannot change
                    after this document is committed
```

---

## PART I — WHY SCRIPT 4 EXISTS

Three scripts and two literature checks were run on TCGA-BRCA.
The central statistical limitation of all three scripts is
identical: TCGA-BRCA claudin-low has 33 events in 268 samples.
Every survival result is underpowered. The confirmed findings
(depth score HR=3.112, composite Treg:effector ratio HR=2.212,
S3-P8, S3-P9, S3-P10) are directionally sound but statistically
fragile in TCGA. External validation is the only path to
clinical credibility.

Two datasets provide that validation:

**METABRIC (Molecular Taxonomy of Breast Cancer International
Consortium):**
- n ≈ 1,980 patients (full cohort)
- Median follow-up > 7 years (versus TCGA ~3–4 years)
- Microarray expression (Illumina HT-12 v3)
- Published survival with >100 expected events in any
  claudin-low-equivalent population
- Available via cBioPortal and EMBL-EBI ArrayExpress
  (E-MTAB-365 / E-MTAB-4439)
- No PAM50 calls in original dataset; claudin-low
  classification must be derived by geometry

**GSE96058 (SCAN-B — Sweden Cancerome Analysis Network —
Breast):**
- n ≈ 3,273 primary breast cancers
- RNA-seq (Illumina HiSeq 2000)
- Median follow-up > 5 years
- Published claudin-low PAM50-like calls available in
  the supplementary data (Brueffer et al. 2018 JCO Precis
  Oncol) — allows direct comparison of geometry-first
  classification vs canonical published classification
- Available via GEO

Script 4 runs both datasets. METABRIC is the primary
validation dataset (more events, longer follow-up).
GSE96058 is the secondary validation dataset (canonical
claudin-low classification available for direct comparison).

---

## PART II — WHAT SCRIPT 4 TESTS

Script 4 has four analysis blocks, each carrying forward
predictions from prior documents:

```
BLOCK 1 — CORE DEPTH SCORE VALIDATION
  Primary question: Does the depth score predict OS in
  claudin-low in METABRIC and GSE96058?
  Predictions: EXT-P1 (from BRCA-S7d)
               EXT-P4 (from BRCA-S7d)

BLOCK 2 — INDIVIDUAL GENE SURVIVAL REPLICATION
  Primary question: Does CLDN3 predict OS in METABRIC?
  Does MKI67 predict OS?
  Predictions: EXT-P2 (from BRCA-S7d)

BLOCK 3 — SCRIPT 3 FEATURE REPLICATION IN METABRIC
  Primary question: Do the memory-group differences
  (CT antigen, depth, Treg ratio) replicate in METABRIC?
  Predictions: EXT-P5 (from BRCA-S7h)
               EXT-P6 (from BRCA-S7h)
               EXT-P7 (from BRCA-S7h)

BLOCK 4 — GSE96058 CANONICAL CLASSIFICATION COMPARISON
  Primary question: Does geometry-first depth score
  perform as well or better in canonical claudin-low
  vs geometry-derived claudin-low?
  Predictions: EXT-P4 (variant)
               EXT-P8 (new — stated below)
```

---

## PART III — COMPLETE PREDICTION SET
### All predictions locked 2026-03-05 before Script 4 runs

---

### BLOCK 1 — DEPTH SCORE SURVIVAL (METABRIC)

---

#### EXT-P1 — DEPTH SCORE PREDICTS OS IN METABRIC CLAUDIN-LOW
**Origin:** BRCA-S7d (locked before Script 2 ran)
**Direction:** HR > 2.0, p < 0.05 for depth score tertile
in METABRIC ESR1-low claudin-low

**Basis:**
TCGA Stratum B produced HR=3.112, p=0.064 with 20 total
events across the tertile extremes. METABRIC has at minimum
3–4× more claudin-low events given its larger cohort size
and longer follow-up. With 60–80 events expected in the
claudin-low ESR1-low equivalent, a true HR of 3.0 would
produce p < 0.001. The prediction is deliberately
conservative (HR > 2.0, p < 0.05) to allow for:
  — Microarray vs RNA-seq platform noise
  — Different clinical era (METABRIC: 1977–2005)
  — Slightly different patient selection

**Population definition for METABRIC:**
Primary (Stratum B equivalent): geometry score ≥ 7 AND
ESR1 expression below METABRIC cohort median.
Secondary (Stratum C equivalent): geometry score ≥ 7
AND IHC ER-negative (where available).
Both will be tested and compared.

**Confidence: HIGH**

**Falsification criterion:**
Not confirmed if p > 0.10 or HR < 1.5 in both populations.

---

#### EXT-P1b — DEPTH SCORE PREDICTS OS IN GSE96058 CLAUDIN-LOW
**Direction:** HR > 2.0, p < 0.05 in ESR1-low claudin-low
in GSE96058

**Basis:** Same reasoning as EXT-P1. GSE96058 has n≈3,273
with published claudin-low classification available as
ground truth comparison.

**Confidence: HIGH**

---

### BLOCK 2 — INDIVIDUAL GENE SURVIVAL REPLICATION (METABRIC)

---

#### EXT-P2 — CLDN3 PREDICTS OS IN METABRIC CLAUDIN-LOW
**Origin:** BRCA-S7d
**Direction:** HR < 0.80 for CLDN3-high vs CLDN3-low
(CLDN3-low = worse OS), p < 0.10

**Basis:**
TCGA confirmed CLDN3 as the strongest depth correlate
(r=-0.641). CLDN3-low directionally predicted worse OS in
TCGA but was underpowered. In METABRIC with more events the
direction should confirm.

**Confidence: MODERATE**

---

#### EXT-P2b — MKI67 PREDICTS OS IN METABRIC CLAUDIN-LOW
**Direction:** HR > 1.3 for MKI67-high vs MKI67-low, p < 0.10

**Basis:**
MKI67 +41% vs normal in Script 1. MKI67-high = higher
proliferation = worse OS. This is the most straightforward
survival predictor and should confirm easily in METABRIC.

**Confidence: MODERATE-HIGH**

---

#### EXT-P3 — IMMUNE SCORE DOES NOT PREDICT OS IN METABRIC
### (NEGATIVE PREDICTION — PRE-CHECKPOINT ERA COHORT)
**Direction:** Total immune score (FOXP3+PDCD1+TIGIT+LAG3)
does NOT predict OS, p > 0.15

**Basis:**
TCGA S2-P3: immune score was null. METABRIC is pre-checkpoint
era (recruited 1977–2005). Immune infiltrate in pre-checkpoint
patients provides no survival benefit because it is Treg-
suppressed (Morel 2017). The null should replicate.

This is a negative prediction — confirming a null is as
important as confirming a positive result. If immune score
DOES predict OS in METABRIC, the framework's immune
cancellation interpretation (BRCA-S7h Part I, 1.5) would
need revision.

**Confidence: MODERATE**

---

### BLOCK 3 — SCRIPT 3 FEATURE REPLICATION (METABRIC)

---

#### EXT-P5 — COMPOSITE TREG:EFFECTOR RATIO PREDICTS OS
### IN METABRIC CLAUDIN-LOW
**Origin:** BRCA-S7h (locked after Script 3 ran)
**Direction:** HR > 1.5, p < 0.05 for composite ratio
(FOXP3+TIGIT)/(CD8A+GZMB+PRF1) tertile in METABRIC
ESR1-low claudin-low

**Basis:**
TCGA Stratum B HR=2.212, p=0.171 — underpowered. METABRIC
provides the power to confirm. The IHC literature (BMC Cancer
2021, METABRIC-validated) already confirmed the CD8/FOXP3
ratio as prognostic in TNBC in METABRIC using IHC data.
This prediction tests whether the RNA-seq version of the
same ratio in the claudin-low subset replicates.

**Note on METABRIC platform:** METABRIC uses microarray
(Illumina HT-12 v3). FOXP3, CD8A, GZMB, TIGIT, PRF1 must
all be present in the array. Availability will be confirmed
in Step 1 of the script. If any gene is missing, the
available subset ratio will be used with a logged note.

**Confidence: MODERATE-HIGH**
(IHC version already validated in METABRIC — RNA proxy
expected to replicate the direction)

---

#### EXT-P6 — MEMORY-LOW HAS HIGHER CT ANTIGEN IN METABRIC
**Origin:** BRCA-S7h
**Direction:** CT antigen composite (GAGE family, CT45
family, STRA8, DPPA2) higher in memory-low claudin-low
vs memory-high, p < 0.001

**Basis:**
TCGA p = 8.86e-09. METABRIC microarray likely has GAGE
family probes. This is the most robust finding from Script 3
and should replicate strongly. The p-value threshold is
tightened (p < 0.001 rather than p < 0.05) to reflect the
expected effect size in a larger dataset.

**Note:** METABRIC microarray probe availability for GAGE
family genes must be confirmed. If GAGE family probes are
absent, the test will use available CT antigen genes only
(CT45, STRA8, DPPA2) and the threshold will be adjusted
accordingly.

**Confidence: MODERATE-HIGH**
Conditional on CT antigen gene availability in METABRIC array.

---

#### EXT-P7 — MEMORY-LOW DOES NOT HAVE WORSE OS THAN
### MEMORY-HIGH IN METABRIC (NEGATIVE PREDICTION)
**Origin:** BRCA-S7h (null result S3-P7 replication)
**Direction:** HR ≈ 1.0, p > 0.30 for memory-low vs
memory-high OS in METABRIC claudin-low

**Basis:**
TCGA S3-P7: HR=1.043, p=0.878 — genuine biological null
interpreted as immune cancellation. If the immune
cancellation interpretation is correct, this null should
replicate in any pre-checkpoint-era cohort (METABRIC
is pre-checkpoint). If the null does NOT replicate
(memory-low shows significantly worse OS in METABRIC),
the immune cancellation interpretation would require
revision — possibly the null in TCGA was a power artefact.

This is the highest-information prediction in Script 4:
  — If replicated (null again): immune cancellation
    interpretation is correct. Pre-checkpoint = no
    Treg removal = CT antigen recognition stays suppressed
    = OS equivalence regardless of depth.
  — If not replicated (memory-low worse OS): the depth
    axis predicts OS directly independent of the immune
    mechanism, and the S3-P7 null in TCGA was a power
    artefact.

Both outcomes are scientifically informative. The
framework predicts the null replicates.

**Confidence: MODERATE**

---

### BLOCK 4 — GSE96058 CANONICAL CLASSIFICATION COMPARISON

---

#### EXT-P4 — DEPTH HR IN CANONICAL CLAUDIN-LOW ≥ DEPTH HR
### IN GEOMETRY-DERIVED CLAUDIN-LOW (GSE96058)
**Origin:** BRCA-S7d
**Direction:** HR(canonical CL) ≥ HR(geometry CL)

**Basis:**
Published claudin-low calls (Brueffer et al. 2018) in
GSE96058 represent a purer claudin-low population than
the geometry-first 10-gene classifier. Purer populations
produce larger depth HRs — established by the
Stratum A → B → C gradient in TCGA (HR: 1.26 → 3.11 → 3.25).
If the canonical classifier identifies a purer set,
its depth HR should be at least as large.

**Confidence: MODERATE**

---

#### EXT-P8 — GEOMETRY-FIRST DEPTH SCORE IDENTIFIES CLAUDIN-LOW
### SAMPLES THAT CANONICAL CLASSIFIER MISSES, WITH EQUIVALENT
### OR BETTER PROGNOSIS SEPARATION (NEW PREDICTION)
**Direction:** Geometry-first claudin-low set overlaps
≥ 70% with canonical claudin-low, and samples in the
geometry set NOT in the canonical set have intermediate
depth scores (shallower than core canonical claudin-low)

**Basis:**
The TCGA LumA contamination analysis (BRCA-S7d) showed
that the geometry-first classifier captures a gradient
including some LumA-like samples at the shallow end.
The canonical classifier (Brueffer et al.) is expected
to identify the deep core. Geometry-first should include
the canonical set plus additional shallow samples. The
non-canonical geometry-first samples should have lower
depth scores (shallower — closer to the LumA edge).

**Confidence: MODERATE**

---

#### EXT-P9 — DEPTH SCORE PREDICTS OS IN CANONICAL CLAUDIN-LOW
### IN GSE96058 (p < 0.01)
**Direction:** HR > 2.0, p < 0.01 in ESR1-low canonical
claudin-low in GSE96058

**Basis:**
GSE96058 n ≈ 3,273 with ≈ 100–200 expected claudin-low
samples and long follow-up. With canonical classification
providing purity and large n providing power, p < 0.01
is expected for the depth score. If TCGA HR=3.112 holds,
this should confirm easily.

**Confidence: HIGH** (conditional on canonical CL n ≥ 50
with sufficient events)

---

## PART IV — COMPLETE PREDICTION REFERENCE TABLE

| ID | Prediction | Dataset | Direction | Confidence |
|----|-----------|---------|-----------|------------|
| EXT-P1 | Depth score HR > 2.0, p < 0.05 in METABRIC CL | METABRIC | Depth-high = worse | HIGH |
| EXT-P1b | Depth score HR > 2.0, p < 0.05 in GSE96058 CL | GSE96058 | Depth-high = worse | HIGH |
| EXT-P2 | CLDN3-low predicts worse OS in METABRIC | METABRIC | HR < 0.80, p < 0.10 | MODERATE |
| EXT-P2b | MKI67-high predicts worse OS in METABRIC | METABRIC | HR > 1.3, p < 0.10 | MOD-HIGH |
| EXT-P3 | Immune score NULL in METABRIC (negative pred.) | METABRIC | p > 0.15 | MODERATE |
| EXT-P4 | Canonical CL depth HR ≥ geometry CL depth HR | GSE96058 | Canonical ≥ geometry | MODERATE |
| EXT-P5 | Composite Treg:effector ratio HR > 1.5, p < 0.05 | METABRIC | High ratio = worse | MOD-HIGH |
| EXT-P6 | Memory-low CT antigen > memory-high, p < 0.001 | METABRIC | Memory-low higher | MOD-HIGH |
| EXT-P7 | Memory-low OS = memory-high OS (null replication) | METABRIC | HR ≈ 1.0, p > 0.30 | MODERATE |
| EXT-P8 | Geometry CL overlaps ≥ 70% with canonical CL | GSE96058 | Overlap ≥ 70% | MODERATE |
| EXT-P9 | Depth score HR > 2.0, p < 0.01 in canonical CL | GSE96058 | Depth-high = worse | HIGH |

---

## PART V — DATA SOURCES AND ACQUISITION

```
METABRIC:
  Primary source: cBioPortal
    https://www.cbioportal.org/study/summary?id=brca_metabric
    Expression: z-scored microarray
    Clinical: overall survival, disease-specific survival,
              ER/PR/HER2 IHC status, PAM50 IntClust

  Secondary source: EMBL-EBI ArrayExpress
    E-MTAB-365 (discovery cohort, n=997)
    E-MTAB-4439 (validation cohort, n=989)

  Key files required:
    data_mrna_illumina_microarray_zscores.txt
    data_clinical_patient.txt
    (survival columns: OS_MONTHS, OS_STATUS)

GSE96058:
  Primary source: GEO
    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058
    Platform: GPL11154 (Illumina HiSeq 2000)
    Expression: FPKM / TPM (will use log2-transformed)
    Clinical: supplement from Brueffer et al. 2018
              JCO Precision Oncology
    Claudin-low calls: from Brueffer 2018 Table S1

  Key files required:
    GSE96058_gene_expression_3273_samples_and_136_replicates.csv.gz
    GSE96058_survival_and_PAM50.tsv (if available as supplement)
    Brueffer 2018 Table S1 (claudin-low calls)

Both datasets are publicly available without access
restrictions. Script 4 attempts download from primary
sources and falls back to secondary sources if primary
URLs are unavailable.
```

---

## PART VI — PLATFORM CONSIDERATIONS

METABRIC uses Illumina HT-12 v3 microarray. Gene-level
probes exist for all claudin/mesenchymal depth score genes
and most immune genes. The following must be verified at
runtime and logged:

```
MUST BE PRESENT (script will abort if missing):
  CLDN3, CLDN4, CLDN7, CDH1, ESR1 (negative depth axis)
  VIM, FN1, SNAI1, ZEB1, CD44 (positive depth axis)

SHOULD BE PRESENT (script continues with available subset):
  FOXP3, CD8A, GZMB, PRF1, TIGIT, PDCD1 (immune panel)
  FOXA1, SPDEF, GATA3 (lineage memory panel)
  GAGE1, GAGE2D, GAGE4, CT45A3, CT45A4, STRA8, DPPA2
  (CT antigen panel — may be absent in microarray)

NOTE ON CT ANTIGEN GENES IN METABRIC:
  GAGE family genes are expressed at very low levels.
  Microarray probes for GAGE1, GAGE2D etc. may have
  low signal or may not be present in the HT-12 v3
  probe set. If CT antigen genes are absent or have
  flat variance across the METABRIC cohort, EXT-P6
  will be tested using available CT antigen genes only.
  If no CT antigen genes are available, EXT-P6 will
  be recorded as DATA MISSING, not NOT CONFIRMED.
```

---

## PART VII — POPULATION CONSTRUCTION IN EXTERNAL DATASETS

Because METABRIC does not have RNA-seq (no clean 01/11
sample type suffix), population construction differs from
TCGA:

```
METABRIC claudin-low classification:
  Step 1: Compute cohort-wide medians for all 10 geometry
  signature genes across all METABRIC tumour samples.

  Step 2: Score each sample:
    +1 for each POS marker above cohort median
    +1 for each NEG marker below cohort median
    Max score = 10

  Step 3: Classify as claudin-low if score ≥ 7
  (same threshold as TCGA)

  Step 4: Stratum B equivalent:
    ESR1 expression below METABRIC claudin-low set median
    AND geometry score ≥ 7

  Step 5: IHC Stratum (where ER IHC data available):
    ER-negative by IHC AND geometry score ≥ 7

GSE96058 claudin-low classification:
  Step 1: Apply geometry classifier (same as above)
  Step 2: Extract canonical claudin-low calls from
  Brueffer 2018 Table S1
  Step 3: Compute overlap: geometry-CL ∩ canonical-CL,
  geometry-CL only, canonical-CL only
  Step 4: Run depth score survival in all three sets
  and compare HRs
```

---

## STATUS BLOCK

```
document:           BRCA-S7j (before_script4.md)
folder:             Cancer_Research/BRCA/DEEP_DIVE/Claudin_Low/
status:             LOCKED
date:               2026-03-05
author:             Eric Robert Lawson
                    OrganismCore
script_4_status:    NOT YET RUN
all_predictions:    stated before any Script 4 output is seen
primary_test:       EXT-P1 (METABRIC depth HR)
                    EXT-P5 (Treg:effector ratio OS in METABRIC)
highest_info_test:  EXT-P7 (memory-low OS null replication —
                    immune cancellation or power artefact)
next_document:      BRCA-S7k
                    Script 4 results and reasoning artifact
datasets:           METABRIC (primary)
                    GSE96058 (secondary, canonical CL calls)
```
