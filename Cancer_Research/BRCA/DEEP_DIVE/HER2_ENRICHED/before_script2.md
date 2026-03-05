# HER2-ENRICHED — SCRIPT 2 BEFORE-DOCUMENT
## Predictions Locked Before Script 2 Runs
## OrganismCore — Document BRCA-S3c
## Date: 2026-03-05

---

## DOCUMENT METADATA

```
document_id:        BRCA-S3c
series:             BRCA Deep Dive — HER2-Enriched
folder:             Cancer_Research/BRCA/DEEP_DIVE/HER2_Enriched/
type:               BEFORE-DOCUMENT
                    All predictions stated before
                    any data is loaded.
                    This document cannot be modified
                    after Script 2 runs.
date:               2026-03-05
author:             Eric Robert Lawson
                    OrganismCore
protocol_version:   Workflow_Protocol.md v2.0
precursor:          BRCA-S3b (Script 1 reasoning artifact)
dataset_planned:    TCGA-BRCA (PAM50 HER2-enriched subset)
                      Bulk RNA-seq
                      Clinical outcomes (OS, RFS)
                      PAM50 subtype calls
                    GSE (trastuzumab response dataset TBD
                      by diagnostic — see Part III)
status:             LOCKED — predictions cannot change
                    after this document is committed
```

---

## PART I — WHAT SCRIPT 1 ESTABLISHED

```
From BRCA-S3b (Script 1 reasoning artifact):

CONFIRMED:
  Type 1 geometry (Slope Arrest) — unambiguous
  FOXA1 retained (+12.7%) — luminal programme still running
  ESR1 -92.5%, PGR -95.2% — ER/PR circuit epigenetically severed
  SOX10 -96.8%, KRT5 -64.6% — no Type 2 / neural crest programme
  EZH2 +176% — intermediate between LumA (+13%) and TNBC (+224%)
  MKI67 +1151% — most proliferative subtype in dataset
  STARD3 +51.9%, MIEN1 +37.0% — amplicon readouts in scRNA-seq
  PCA distance HER2 → Luminal: 0.81 (closer than LumA at 0.92)

DEPTH AXIS INVERSION (novel):
  Deeper HER2 cells LOSE the HER2 programme:
    ERBB3  r=-0.264  (strongest negative correlation with depth)
    CDH1   r=-0.247
    AR     r=-0.246
    AKT1   r=-0.234
    ERBB2  r=-0.172
  No strong positive correlations with depth.
  The HER2 attractor is unstable at its deep end.
  Deep = pre-resistant = phenotypically undefined.

FRAMEWORK LESSON CONFIRMED:
  Copy-number amplification is attenuated in log-normalised
  scRNA-seq. STARD3 > ERBB2 as amplicon readout.
  All future predictions on copy-number events adjusted.
```

---

## PART II — WHAT SCRIPT 2 MUST DO

```
Script 2 has FOUR independent tasks:

TASK 1 — TCGA-BRCA BULK VALIDATION
  Load TCGA-BRCA expression (HiSeqV2, log2 RSEM+1)
  Load TCGA-BRCA clinical matrix (PAM50, survival)
  Isolate PAM50 = HER2-enriched subset
  Confirm amplicon signature: ERBB2, STARD3, GRB7, MIEN1
  Confirm ERBB3 suppression relative to LumA
  Confirm EZH2 elevation relative to LumA
  Confirm FOXA1 retention relative to TNBC (Basal)
  Compute depth score:
    Depth = inverse of mean(ERBB3, CDH1, AR) normalised
    (low ERBB3 + low CDH1 + low AR = deep = pre-resistant)
  Test depth score vs OS (log-rank, Cox PH)
  Test depth score vs RFS where available

TASK 2 — TRASTUZUMAB RESISTANCE DATASET
  Identify and download a dataset with:
    HER2+ tumours
    Pre-treatment expression profiling
    Trastuzumab (or trastuzumab-containing) response annotation
    pCR or clinical response as endpoint
  Candidates (diagnostic will confirm which is accessible):
    GSE37946 — trastuzumab response, HER2+ neoadjuvant
    GSE50948 — NeoSphere trial, HER2+ neoadjuvant
    GSE55348 — lapatinib + trastuzumab, HER2+
    GSE66399 — PAMELA trial, HER2+ dual blockade
  Test: ERBB3 expression as pCR predictor in HER2+
  Test: Depth score as resistance predictor

TASK 3 — EZH2 / ESR1 AXIS IN BULK
  In TCGA-BRCA HER2-enriched:
    Test r(EZH2, ESR1) — do EZH2-high HER2 tumours
    have deeper ESR1 silencing?
  In available trastuzumab dataset:
    Test whether EZH2-high HER2+ tumours are
    more resistant (lower pCR)

TASK 4 — FOXA1 ENDOCRINE SENSITIVITY
  In TCGA-BRCA HER2-enriched:
    Test r(FOXA1, ESR1) in HER2+ tumours
    Confirm FOXA1 retention across the cohort
    Identify FOXA1-high vs FOXA1-low subgroups
    Test FOXA1-high vs FOXA1-low for survival difference
```

---

## PART III — TECHNICAL PREREQUISITES
### Must be resolved before predictions are stated

```
The diagnostic script must confirm:

1. TCGA-BRCA expression file:
   Candidate URL:
     https://tcga.xenahubs.net/download/
     TCGA.BRCA.sampleMap/HiSeqV2.gz
   Fallback:
     https://pancanatlas.xenahubs.net/download/
     TCGA.BRCA.sampleMap/HiSeqV2.gz
   Expected: ~1200 samples, ~20000 genes, log2 RSEM+1

2. TCGA-BRCA clinical matrix:
   Candidate URL:
     https://tcga.xenahubs.net/download/
     TCGA.BRCA.sampleMap/BRCA_clinicalMatrix
   Contains: PAM50, ER/PR/HER2 IHC, survival

3. TCGA-BRCA survival table:
   Candidate URL:
     https://pancanatlas.xenahubs.net/download/
     Survival_SupplementalTable_S1_20171025_xena_sp
   Contains: OS, DSS, DFI, PFI with events and times

4. Trastuzumab response dataset:
   Priority order (diagnostic to confirm accessibility):
     1. GSE37946
     2. GSE50948
     3. GSE55348
     4. GSE66399
   Minimum required:
     HER2+ annotated samples
     Pre-treatment expression
     pCR or response annotation

These are confirmed or alternatives substituted
BEFORE Script 2 runs. Script 2 does not attempt
URLs that have not been verified by the diagnostic.
```

---

## PART IV — PREDICTIONS
### All stated 2026-03-05 before Script 2 runs

---

### S2-P1 — ERBB3 AS TRASTUZUMAB RESPONSE BIOMARKER

```
PREDICTION:
  In a HER2+ neoadjuvant trastuzumab dataset,
  ERBB3 expression will be HIGHER in pCR patients
  than in non-pCR (residual disease) patients.

DIRECTION: ERBB3-high → pCR (trastuzumab sensitive)
           ERBB3-low  → RD (trastuzumab resistant)

GEOMETRY BASIS:
  ERBB3 r=-0.264 with depth in scRNA-seq.
  Deeper HER2 cells lose ERBB3 first.
  ERBB2/ERBB3 heterodimer is the oncogenic pair.
  Loss of ERBB3 dismantles the signalling complex
  that trastuzumab relies on to trigger ADCC
  and block downstream PI3K activation.

MAGNITUDE PREDICTION:
  ERBB3 mean in pCR group > ERBB3 mean in RD group
  Difference statistically significant (p < 0.05)
  ERBB3 AUC as pCR predictor > 0.60

FALSIFICATION CRITERION:
  ERBB3 does not differ between pCR and RD, OR
  ERBB3-low correlates with pCR (opposite direction)
```

---

### S2-P2 — DEPTH SCORE PREDICTS OS IN TCGA HER2-ENRICHED

```
PREDICTION:
  A depth score computed as:
    Depth = 1 - norm(mean[ERBB3, CDH1, AR])
    (low expression of all three = high depth = resistant)
  will stratify TCGA-BRCA HER2-enriched patients
  into distinct survival groups.

DIRECTION: depth-high (low ERBB3/CDH1/AR) → worse OS
           depth-low  (high ERBB3/CDH1/AR) → better OS

GEOMETRY BASIS:
  Depth axis inversion confirmed in Script 1.
  Deeper cells shed ERBB3, CDH1, AR simultaneously.
  This is the pre-resistant subpopulation.
  In a bulk tumour, higher depth score = higher fraction
  of pre-resistant cells = worse clinical outcome.

MAGNITUDE PREDICTION:
  Log-rank p < 0.05 separating depth-high vs depth-low
  Hazard ratio > 1.5 for depth-high group
  Median OS separation > 12 months

FALSIFICATION CRITERION:
  No survival separation by depth score, OR
  depth-high associates with better OS
```

---

### S2-P3 — EZH2 HIGH IN HER2 CORRELATES WITH DEEPER ESR1 SILENCING

```
PREDICTION:
  In TCGA-BRCA HER2-enriched bulk tumours:
  r(EZH2, ESR1) will be NEGATIVE.
  EZH2-high HER2 tumours will have lower ESR1.

GEOMETRY BASIS:
  EZH2 is the epigenetic gate that silences ESR1/PGR.
  This was confirmed in TNBC (EZH2 r=+0.274 with depth,
  depth = lower ESR1).
  In HER2, EZH2 +176% drives ESR1 -92.5%.
  Within HER2, the EZH2/ESR1 inverse relationship
  should be visible as a negative correlation.

SECONDARY PREDICTION:
  EZH2-high HER2 tumours will also show lower FOXA1
  (FOXA1 is the upstream pioneer TF; if EZH2 is high
  enough, it will suppress FOXA1 as well, not just ESR1)
  r(EZH2, FOXA1) negative in HER2-enriched subset

MAGNITUDE PREDICTION:
  r(EZH2, ESR1) < -0.15 in HER2-enriched
  p < 0.05

FALSIFICATION CRITERION:
  r(EZH2, ESR1) positive or near-zero in HER2-enriched
```

---

### S2-P4 — STARD3 IS A MORE STABLE AMPLICON READOUT THAN ERBB2 IN BULK

```
PREDICTION:
  In TCGA-BRCA HER2-enriched vs LumA comparison:
  STARD3 fold-change will be comparable to or larger
  than ERBB2 fold-change.
  Both will show clear amplicon signal in bulk RNA-seq
  (unlike scRNA-seq where ERBB2 was attenuated).

SECONDARY PREDICTION:
  In bulk RNA-seq, ERBB2 will show the expected large
  fold-change (>5x, >+400%) that was absent in scRNA-seq.
  This confirms the scRNA-seq lesson: the attenuation
  was a normalisation artefact, not a biological finding.

GEOMETRY BASIS:
  Script 1 established that scRNA-seq attenuates
  copy-number-driven fold-change.
  Bulk RNA-seq does not apply per-cell normalisation.
  The full amplicon signal should be visible in bulk.

MAGNITUDE PREDICTION:
  ERBB2: >+400% in HER2-enriched vs LumA in TCGA
  STARD3: >+100% in same comparison

FALSIFICATION CRITERION:
  ERBB2 fold-change in bulk < +100%
  (would indicate dataset or subtyping problem)
```

---

### S2-P5 — FOXA1 RETENTION SEPARATES HER2 FROM TNBC IN BULK

```
PREDICTION:
  In TCGA-BRCA:
  FOXA1 mean expression:
    HER2-enriched > Basal-like (TNBC)
  The retention confirmed in scRNA-seq holds in bulk.

SECONDARY PREDICTION:
  Within HER2-enriched, FOXA1-high tumours will show
  higher ESR1 (partial ER+) and better response
  to endocrine therapy where annotated.
  r(FOXA1, ESR1) POSITIVE within HER2-enriched.
  (Unlike TNBC where FOXA1 is absent and ESR1 is near-zero
  for the entire subtype — within-subtype correlation
  is not meaningful when both genes are near-zero.)

GEOMETRY BASIS:
  FOXA1 retention is the mechanistic basis for why
  EZH2 inhibition can restore ESR1 in HER2+ but not TNBC.
  FOXA1 is the pioneer TF that rebinds ESR1 enhancers
  once H3K27me3 is cleared by EZH2 inhibition.
  If FOXA1 is present (HER2), restoration is possible.
  If FOXA1 is absent (TNBC), restoration requires
  both EZH2 inhibition AND FOXA1 re-expression —
  a harder combinatorial requirement.

MAGNITUDE PREDICTION:
  FOXA1 in HER2-enriched: mean > 8.0 (log2 RSEM+1)
  FOXA1 in Basal-like: mean < 6.0
  r(FOXA1, ESR1) within HER2-enriched > +0.20

FALSIFICATION CRITERION:
  FOXA1 not significantly different between HER2 and Basal,
  OR r(FOXA1, ESR1) negative in HER2-enriched
```

---

### S2-P6 — DEEP HER2 CELLS MAP GRADE 3 IN BULK

```
PREDICTION:
  In TCGA-BRCA HER2-enriched:
  Depth score (ERBB3-low, CDH1-low, AR-low) will
  correlate positively with histological grade.
  Grade 3 tumours will have higher depth scores
  than Grade 1/2 tumours.

SECONDARY PREDICTION:
  MKI67 / Ki-67 will correlate positively with depth score.
  r(depth, MKI67) > 0 in HER2-enriched.
  (Deeper cells in scRNA-seq had marginally higher KRT5
  but mainly shed differentiation markers — in bulk,
  the dedifferentiated deep fraction should associate
  with proliferation markers at the tumour level.)

GEOMETRY BASIS:
  The depth axis identifies the most dedifferentiated
  fraction within HER2-enriched.
  Dedifferentiation and proliferation are coupled
  in the wrong-slope geometry — cells that shed
  ERBB3/CDH1/AR are also more likely to be cycling.
  Grade 3 is histologically defined by high mitotic rate
  and nuclear pleomorphism — both consistent with
  the depth-high, phenotypically undefined population.

MAGNITUDE PREDICTION:
  Mean depth score: Grade 3 > Grade 2 > Grade 1
  p < 0.05 (ANOVA or Kruskal-Wallis across grades)

FALSIFICATION CRITERION:
  No grade association with depth score
```

---

### S2-P7 — AR-LOW WITHIN HER2 IS A DISTINCT SUBPOPULATION

```
PREDICTION:
  Within TCGA-BRCA HER2-enriched, AR expression
  will be bimodally or continuously distributed.
  AR-low HER2+ tumours will have:
    Lower ERBB3
    Lower CDH1
    Higher EZH2
    Worse OS than AR-high HER2+ tumours

GEOMETRY BASIS:
  AR r=-0.246 with depth in scRNA-seq.
  AR-low cells are the deep, pre-resistant fraction.
  AR tracks with the luminal differentiation programme
  in HER2 (AR is a luminal marker in this context,
  opposite to its role in TNBC/LAR where it defines
  a luminal-within-basal subtype).
  AR-low HER2+ is the furthest from luminal identity
  within the HER2 attractor.

SECONDARY PREDICTION:
  AR-low HER2+ overlaps with the CDH1-low, ERBB3-low
  population identified in S2-P1 and S2-P2.
  These are the same cells viewed from different angles
  of the depth axis.
  A composite three-marker score (ERBB3 + CDH1 + AR)
  is more predictive than any single marker alone.

MAGNITUDE PREDICTION:
  Log-rank p < 0.05 separating AR-high vs AR-low
  in HER2-enriched OS analysis
  Hazard ratio > 1.3 for AR-low group

FALSIFICATION CRITERION:
  No OS difference between AR-high and AR-low
  in HER2-enriched, OR AR-low associates with better OS
```

---

## PART V — WHAT SCRIPT 2 DOES NOT PREDICT

```
Script 2 does not predict:

1. Absolute pCR rates in any dataset.
   We predict direction and biomarker utility,
   not absolute response frequencies.

2. The specific trastuzumab dataset to be used.
   The diagnostic determines which dataset is
   accessible. The predictions hold regardless
   of which qualifying dataset is used.

3. Mechanism of ERBB3 downregulation.
   We observe the correlation. We predict the
   clinical consequence. The molecular mechanism
   (transcriptional, post-translational, or
   epigenetic) is not predicted here.

4. Whether EZH2 inhibition plus trastuzumab
   will show clinical benefit.
   We predict the biomarker correlation.
   The clinical efficacy question requires
   a clinical trial — outside this framework.

5. ERBB2 amplification status in individual samples.
   We are working with expression data only.
   FISH confirmation is not modelled here.
```

---

## PART VI — COMPLETE PREDICTION REFERENCE

```
S2-P1:  ERBB3-high predicts trastuzumab pCR
        Direction: ERBB3-high → pCR, ERBB3-low → RD
        Test:      pre-treatment HER2+ dataset, pCR endpoint

S2-P2:  Depth score (ERBB3/CDH1/AR inverse) predicts OS
        Direction: depth-high → worse OS
        Test:      TCGA-BRCA HER2-enriched, log-rank + Cox

S2-P3:  r(EZH2, ESR1) negative in HER2-enriched bulk
        Direction: EZH2-high → lower ESR1
        Test:      TCGA-BRCA HER2-enriched Pearson r

S2-P4:  ERBB2 fold-change visible in bulk (>+400%)
        STARD3 co-elevated (>+100%)
        Test:      TCGA-BRCA HER2 vs LumA comparison

S2-P5:  FOXA1 retained in HER2 vs Basal in bulk
        r(FOXA1, ESR1) positive within HER2-enriched
        Test:      TCGA-BRCA subtype comparison

S2-P6:  Depth score correlates with Grade 3 in bulk
        Direction: higher depth → higher grade
        Test:      TCGA-BRCA HER2-enriched, grade annotation

S2-P7:  AR-low HER2+ has worse OS than AR-high HER2+
        Direction: AR-low → worse OS
        Test:      TCGA-BRCA HER2-enriched, log-rank
```

---

## STATUS BLOCK

```
document_id:    BRCA-S3c
status:         LOCKED
locked_date:    2026-03-05
predictions:    7 (S2-P1 through S2-P7)
next_action:    Run diagnostic script to confirm
                data URLs, then run Script 2
next_document:  BRCA-S3d (Script 2 reasoning artifact)
```

---

*OrganismCore — Eric Robert Lawson — 2026-03-05*
