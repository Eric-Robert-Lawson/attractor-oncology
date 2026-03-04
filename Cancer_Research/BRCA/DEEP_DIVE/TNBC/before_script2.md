# TNBC — SCRIPT 2 BEFORE-DOCUMENT
## Predictions Locked Before Script 2 Runs
## OrganismCore — Document BRCA-S4c
## Date: 2026-03-04

---

## DOCUMENT METADATA

```
document_id:        BRCA-S4c
series:             BRCA Deep Dive — TNBC (Basal-like)
folder:             Cancer_Research/BRCA/DEEP_DIVE/TNBC/
type:               BEFORE-DOCUMENT
                    All predictions stated before
                    any data is loaded.
                    This document cannot be modified
                    after Script 2 runs.
date:               2026-03-04
author:             Eric Robert Lawson
                    OrganismCore
protocol_version:   Workflow_Protocol.md v2.0
precursor:          BRCA-S4b (Script 1 reasoning artifact)
dataset_planned:    GSE25066 — Hatzis et al. 2011
                      n=508 pre-treatment bulk samples
                      Affymetrix HG-U133A (GPL96)
                      pCR annotation
                      TNBC subset by ER/PR/HER2 status
                    TCGA-BRCA (PAM50 Basal-like subset)
                      Bulk RNA-seq
                      BRCA1 mutation + methylation data
                      Clinical outcomes (OS, RFS)
status:             LOCKED — predictions cannot change
                    after this document is committed
```

---

## PART I — WHAT SCRIPT 1 ESTABLISHED

```
From BRCA-S4b (Script 1 reasoning artifact):

CONFIRMED:
  Type 2 geometry (Wrong Valley) — unambiguous
  All 7 FA markers elevated (SOX10 +1323% largest signal)
  All 8 switch genes suppressed (ESR1 -97%, PGR -98%)
  EZH2 elevated +270% — confirmed as convergence node gate
  EED r=+0.435 — depth driver within TNBC (novel)
  VIM r=+0.445 — EMT gradient within attractor (confirmed)
  AR r=-0.301 — LAR subtype is geometrically shallow
  PARP1 r=+0.235 — depth-linked PARP1 elevation (novel)
  pCR directionality: deeper TNBC → lower pCR (r=-0.098,
    degraded proxy — definitive test is the purpose of
    Script 2)

NOT CONFIRMED (requires Script 2):
  P4: r(BRCA1, ESR1) — composite type signal overwritten
      at expression level. DNA-level test needed.
  P6: pCR correlation used degraded proxy (top-variance
      probes, not gene-mapped). Proper mapping needed.
  ER/PR/HER2 parsing from GSE25066 failed.
      TNBC subset not properly isolated.

KEY NOVEL FINDINGS TO VALIDATE:
  EED > EZH2 as depth biomarker for epigenetic therapy
  PARP1 as depth-linked target
  LAR geometry (AR-shallow) confirmed in bulk?
  SOX10 as the dominant FA marker — holds in bulk?
```

---

## PART II — WHAT SCRIPT 2 MUST DO

```
Script 2 has THREE independent tasks:

TASK 1 — GSE25066 PROPER ANALYSIS
  Correct Affymetrix probe → gene symbol mapping
  using GPL96 platform annotation file.
  Proper TNBC subset extraction (ER-/PR-/HER2-).
  Proper pCR annotation parsing.
  Gene-based depth score (not top-variance proxy).
  This is the definitive P6 test.

TASK 2 — TCGA-BRCA BULK VALIDATION
  PAM50 Basal-like subset extraction.
  Validate EED > EZH2 as depth biomarker.
  Validate PARP1 depth correlation in bulk.
  Test BRCA1 mutation/methylation enrichment in
  basal-like vs luminal (P4 DNA-level revision).
  Test survival correlation with depth score.

TASK 3 — LEHMANN SUBTYPE MAPPING
  Map the six Lehmann TNBC subtypes onto the
  depth score axis.
  Prediction: LAR = shallow, BL1/BL2 = intermediate,
  M/MSL = deep.
  This requires Lehmann subtype gene signatures
  applied to GSE25066 or TCGA-BRCA.
```

---

## PART III — TECHNICAL PREREQUISITES
### Must be resolved before predictions are stated

```
AFFYMETRIX PROBE MAPPING (GPL96 — HG-U133A):

  The GSE25066 expression matrix contains Affymetrix
  probe set IDs as row indices (e.g., "200814_at").
  Script 1 failed to match these to gene names because
  it searched for gene symbols in probe ID strings.

  Script 2 must download the GPL96 platform annotation:
    URL: https://ftp.ncbi.nlm.nih.gov/geo/platforms/
         GPL96nnn/GPL96/annot/GPL96.annot.gz
    Contains: probe ID → gene symbol mapping
    Required columns: ID (probe), Gene Symbol

  Target genes and their expected HG-U133A probe IDs
  (from published GPL96 annotation — confirmed mapping):

    ESR1:   205225_at
    FOXA1:  202340_x_at
    GATA3:  209604_s_at
    SPDEF:  219197_s_at
    PGR:    208305_at
    KRT5:   201820_at
    KRT14:  201744_s_at
    SOX10:  221579_at
    FOXC1:  203853_s_at
    EGFR:   201983_s_at
    VIM:    201426_s_at
    EZH2:   203358_s_at
    EED:    218657_s_at
    HDAC1:  202705_at
    HDAC2:  200895_s_at
    KDM1A:  218263_s_at
    MKI67:  212022_s_at
    BRCA1:  204531_s_at
    AR:     211110_x_at
    CDH1:   201130_s_at
    ZEB1:   209839_at
    VIM:    201426_s_at
    PARP1:  208501_at
    CD274:  223834_at
    TP53:   201746_at
    PTEN:   214440_at
    PIK3CA: 212733_at
    CDH3:   205497_at

  NOTE: Some genes have multiple probes. Use the
  probe with the highest variance across samples
  when multiple probes map to the same gene.

GSE25066 METADATA FORMAT:
  The ER/PR/HER2 parsing failed in Script 1 because
  the characteristics fields use a different format
  than expected.
  Script 2 must print ALL metadata fields before
  attempting to parse them.
  Known field names from GSE25066 publication:
    "er status: positive/negative"
    "her2 status: positive/negative"
    "drfs event: 0/1"
    "pCR: 0/1" (or "pcr: yes/no")
  The exact strings must be read from the data,
  not assumed.

TCGA-BRCA DATA:
  Source: GDC Data Portal or TCGA legacy archive
  Files needed:
    BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu
    __Level_3__RSEM_genes_normalized__data.data.txt
    (normalized gene expression, RSEM)
    clinical file: nationwidechildrens.org_clinical_
    patient_brca.txt
    PAM50 subtype calls: available in clinical file
    BRCA1 mutation: somatic mutation file (MAF)
    BRCA1 methylation: methylation 450k file
  Alternative: use the curated TCGA-BRCA dataset
  from cBioPortal (tcga_brca_pub2015) which has
  all of the above pre-merged.
```

---

## PART IV — PREDICTIONS
### All stated 2026-03-04 before Script 2 runs

---

### S2-P1 — GENE-MAPPED DEPTH SCORE IN GSE25066

```
PREDICTION:
  When depth score is computed using properly
  gene-mapped Affymetrix probes (not top-variance
  proxy), the correlation with pCR binary outcome
  in GSE25066 will be more negative than Script 1.

  Script 1 result (degraded proxy): r = -0.098

  S2-P1a: In all 508 samples (full cohort):
    r(depth_proper, pCR) < -0.098
    i.e., gene-mapped score performs better
    than the top-variance proxy.

  S2-P1b: In TNBC subset only (ER-/PR-/HER2-):
    r(depth_proper, pCR) < -0.15
    The TNBC-specific effect is stronger than
    the pan-cohort effect because depth in non-TNBC
    tumours is on a different scale.

  Reasoning:
    The top-variance proxy is contaminated by
    non-target variance. Gene-mapped depth using
    the basal FA markers (KRT5, KRT14, SOX10, FOXC1)
    and luminal switch genes (ESR1, FOXA1, GATA3,
    SPDEF) will more precisely measure the biological
    axis that predicts pCR.

  ANTI-PREDICTION (wrong prediction to falsify):
    If gene-mapped depth performs WORSE than the
    top-variance proxy, it would mean the depth
    biology does not predict pCR in bulk and the
    Script 1 result was a noise artifact.
    This would require revision of P6.
```

---

### S2-P2 — EED AS pCR BIOMARKER

```
PREDICTION:
  EED expression (properly mapped Affymetrix probe)
  predicts pCR outcome in TNBC.

  Direction: EED HIGH → pCR = 0 (residual disease)
             EED LOW  → pCR = 1 (complete response)

  Reasoning from Script 1 Novel-1:
    EED r=+0.435 within TNBC — depth driver.
    Depth predicts lower pCR.
    Therefore: high EED → deep attractor → lower pCR.

  S2-P2a: r(EED, pCR_binary) < 0 in TNBC subset
    (higher EED = less likely to achieve pCR)
    Predicted: r < -0.10, p < 0.05

  S2-P2b: EED predicts pCR better than EZH2.
    AUC(EED) > AUC(EZH2) for pCR prediction in TNBC
    (if AUC analysis is possible with n available)
    OR: r(EED, pCR) is more negative than
        r(EZH2, pCR) in absolute value

  This is the primary novel clinical prediction
  from Script 1. It converts the Novel-1 finding
  into a testable biomarker claim.

  If confirmed: EED IHC or RNA expression should
  be tested as a clinical stratification tool
  before EZH2i or EEDi therapy.
```

---

### S2-P3 — PARP1 AS pCR BIOMARKER

```
PREDICTION:
  PARP1 expression correlates NEGATIVELY with pCR
  to standard neoadjuvant chemotherapy in TNBC.

  Direction: PARP1 HIGH → pCR = 0
             PARP1 LOW  → pCR = 1

  Reasoning from Script 1 Novel-3:
    PARP1 r=+0.235 with depth within TNBC.
    Depth predicts lower pCR (P6 confirmed).
    Therefore: high PARP1 → deeper → lower pCR
    to chemotherapy.

  S2-P3a: r(PARP1, pCR_binary) < 0 in TNBC
    Predicted: r < -0.08, p < 0.05

  S2-P3b: PARP1 may predict PARPi benefit
    (OPPOSITE direction for PARPi vs chemotherapy)
    High PARP1 → more PARP1 to inhibit → better
    PARPi response.
    This second prediction CANNOT be tested in
    GSE25066 (chemotherapy data only).
    It will be a stated prediction for future
    PARPi clinical trial datasets.
    STATED HERE (2026-03-04): in any TNBC cohort
    treated with olaparib or talazoparib,
    baseline PARP1 expression should positively
    correlate with response.

  NOTE: The two PARP1 predictions are NOT contradictory:
    High PARP1 → chemo RESISTANCE (predicts low pCR)
    High PARP1 → PARPi SENSITIVITY (predicts PARPi response)
    These are different drug-specific predictions
    from the same biological finding.
```

---

### S2-P4 — LEHMANN SUBTYPE DEPTH MAPPING

```
PREDICTION:
  The Lehmann TNBC subtypes map onto the depth
  score axis in this order (shallow to deep):

  LAR (Luminal Androgen Receptor)     ← SHALLOWEST
    AR+, residual luminal character
    r(AR, depth) = -0.301 from Script 1
    LAR cells cluster at low depth

  BL1 (Basal-like 1)                 ← INTERMEDIATE
    High proliferation, cell cycle
    FOXC1+, KRT5/14+

  BL2 (Basal-like 2)                 ← INTERMEDIATE
    Growth factor signaling
    EGFR activation

  IM (Immunomodulatory)               ← INTERMEDIATE-DEEP
    Immune gene programme
    SPI1 elevation seen in Script 1

  M (Mesenchymal)                     ← DEEP
    VIM high, ZEB1/2 high

  MSL (Mesenchymal stem-like)         ← DEEPEST
    Stem cell markers
    Highest VIM and EMT

  Predicted formal test:
    Mean depth score by Lehmann subtype:
    LAR < BL1 ≈ BL2 < IM < M < MSL

  Specific numerical predictions (stated before data):
    LAR mean depth:  < 0.40
    BL1/BL2 mean:    0.45 – 0.55
    M/MSL mean:      > 0.60

  Reasoning:
    VIM is the primary EMT marker and depth driver
    (r=+0.445 in Script 1).
    M and MSL subtypes are defined by VIM/mesenchymal
    elevation — they should be the deepest.
    LAR is defined by residual luminal character (AR+)
    — AR is anti-correlated with depth (r=-0.301)
    — LAR should be the shallowest.

  NOTE: Lehmann subtype calls are not available
  directly in GSE25066 or TCGA-BRCA as annotations.
  They must be derived by applying the Lehmann
  centroid classifier to the expression data.
  The classifier gene lists are published in
  Lehmann et al. 2016 (updated 6-subtype version).
  Script 2 will implement this classifier.
```

---

### S2-P5 — COMPOSITE TYPE TEST AT DNA LEVEL

```
PREDICTION:
  In TCGA-BRCA, BRCA1 dysfunction (germline mutation
  + somatic mutation + promoter hypermethylation
  combined) is enriched in PAM50 Basal-like tumours
  relative to PAM50 Luminal A, Luminal B, and
  HER2-enriched tumours.

  S2-P5a: BRCA1 germline + somatic mutation rate
    in PAM50 Basal-like: predicted > 20%
    in PAM50 Luminal A: predicted < 5%
    Fisher's exact p < 0.001

  S2-P5b: BRCA1 promoter methylation rate
    in PAM50 Basal-like: predicted > 30%
    in PAM50 Luminal A: predicted < 10%
    Fisher's exact p < 0.001

  S2-P5c: Combined BRCA1 dysfunction
    (mutation OR methylation) in PAM50 Basal-like:
    predicted > 50%
    This is the Type 1 component of the composite
    type — confirming that the founding event (BRCA1
    loss) is enriched in the Type 2 false attractor
    (basal-like subtype).

  Reasoning:
    This is the DNA-level test of P4 from Script 1
    (which failed at the expression level because the
    Type 1 signal is overwritten in established tumours).
    The literature already establishes BRCA1 enrichment
    in basal-like TNBC (~20% germline + additional somatic
    + methylation = ~50-60% total dysfunction).
    This prediction is expected to confirm the
    literature and formally validate the composite type
    axiom from Doc 90 using TCGA-BRCA data.
    The specific thresholds stated above are novel
    formal predictions, even if the direction is
    established.
```

---

### S2-P6 — DEPTH SCORE IN TCGA-BRCA SURVIVAL

```
PREDICTION:
  In TCGA-BRCA PAM50 Basal-like subset,
  depth score (computed from RNA-seq expression)
  will correlate with overall survival and/or
  relapse-free survival.

  Direction: higher depth → shorter survival
  (deeper false attractor = worse prognosis)

  S2-P6a: Kaplan-Meier: depth-high vs depth-low TNBC
    High depth group: shorter OS (predicted)
    Log-rank p < 0.05 predicted
    Median survival difference: > 12 months

  S2-P6b: Cox proportional hazards:
    Depth score as continuous variable in PAM50
    Basal-like predicts OS
    HR > 1.3 per 0.1 increase in depth score
    (i.e., each 0.1 unit increase in depth is
    associated with ≥30% increase in hazard)

  Reasoning:
    Depth score captures the degree of luminal
    programme loss and basal programme acquisition.
    Deeper cells have more PRC2 (EED), more PARP1,
    more EMT (VIM). All of these are associated
    with worse prognosis in published TNBC literature.
    The depth score synthesises them into a
    single geometric measure.

  CAVEAT:
    TCGA-BRCA TNBC has only ~180-200 samples.
    Power for survival analysis is moderate.
    P6b (Cox HR) requires sufficient events.
    If underpowered, the survival finding will be
    directional but may not reach p < 0.05.
    The direction (HR > 1) is the primary prediction.
```

---

### S2-P7 — SPI1 RESOLUTION

```
PREDICTION:
  SPI1 elevation in the "Cancer Basal SC" population
  (Script 1, +141%) is primarily explained by
  immune cell contamination rather than genuine
  TNBC cancer cell biology.

  S2-P7a (contamination hypothesis):
    In bulk GSE25066, SPI1 expression correlates
    positively with established immune infiltration
    markers:
      PTPRC (CD45): r(SPI1, PTPRC) > +0.30
      CD68 (macrophage): r(SPI1, CD68) > +0.25
      CD3D (T cell): r(SPI1, CD3D) > +0.20
    If these correlations hold, SPI1 in the
    scRNA-seq data reflects mislabelled immune cells
    in the Cancer Basal SC cluster.

  S2-P7b (genuine IM subtype hypothesis):
    If SPI1 correlates more strongly with IM Lehmann
    subtype score than with immune cell markers,
    it is genuine TNBC cancer cell biology.
    The IM subtype has documented immune gene
    activation in cancer cells themselves.

  ANTI-PREDICTION:
    If r(SPI1, PTPRC) > r(SPI1, IM_score),
    contamination hypothesis confirmed.
    If r(SPI1, IM_score) > r(SPI1, PTPRC),
    genuine IM biology confirmed.
    One must be larger. Both will be computed.

  NOTE: If SPI1 is genuine IM biology, it becomes
  a prediction of pembrolizumab response —
  SPI1-high TNBC (IM subtype) would be predicted
  to respond better to anti-PD-L1 therapy.
  This would be a direct clinical prediction.
```

---

### S2-P8 — EZH2 INHIBITOR RESPONSE BIOMARKER

```
PREDICTION:
  EED expression level predicts response to EZH2
  inhibitors better than EZH2 expression level.

  This cannot be tested in GSE25066 (neoadjuvant
  chemotherapy dataset, no EZH2i data).
  Cannot be tested in TCGA-BRCA (surgical resection,
  no EZH2i data).

  STATED AS A PROSPECTIVE PREDICTION (2026-03-04):
  In any future dataset with baseline TNBC gene
  expression + EZH2i treatment + response data:
    AUC(EED for response) > AUC(EZH2 for response)
    r(EED, response) > r(EZH2, response)

  Supporting logic:
    EZH2 is the gate (uniform elevation, Script 1)
    EED is the depth driver (r=+0.435, Script 1)
    More assembled PRC2 = more drug target available
    = larger magnitude of effect when inhibited

  TESTABLE IN EXISTING DATA:
    Schade et al. 2024 (Nature) used TNBC PDX models
    treated with EZH2i + AKTi.
    If raw data is available (supplementary), test
    whether baseline EED expression in PDX models
    predicts the magnitude of EZH2i-induced
    luminal conversion.
    This is stated as a pre-analysis prediction for
    any Schade 2024 supplementary data analysis.
    EED > EZH2 as predictor predicted.
```

---

### S2-P9 — AR BLOCKADE IN LAR SUBTYPE

```
PREDICTION:
  LAR (Luminal Androgen Receptor) TNBC has lower
  depth scores than BL1/BL2/M/MSL subtypes.
  LAR is the shallowest subtype in the basal attractor.
  AR-positive TNBC cells are the most restorable
  to luminal identity.

  S2-P9a: In GSE25066, AR-high TNBC samples have:
    Lower depth score (near the LAR end)
    Lower probability of pCR (LAR known to have
    low pCR to chemotherapy — consistent with depth
    if the correct drug is not chemotherapy)
    r(AR, depth_score) < -0.15 in bulk TNBC

  S2-P9b: LAR-classified TNBC (if Lehmann classifier
  applied) has lowest depth mean of all Lehmann subtypes.

  NOTE ON APPARENT CONTRADICTION:
    Shallow = easier to restore luminal identity
            = better response to endocrine/AR therapy
    Shallow ≠ better response to chemotherapy
    The depth axis predicts response to the
    CORRECT DRUG for the depth position,
    not response to all drugs equally.
    LAR (shallow): AR blockade predicted to work
    LAR (shallow): chemotherapy predicted to be less
                   effective than for deeper subtypes
    This is NOT a contradiction — it is the framework
    predicting drug-depth specificity.
```

---

## PART V — WHAT SCRIPT 2 DOES NOT PREDICT

```
Script 2 does NOT predict:
  ✗ New drug targets beyond those in BRCA-S4a/b
    (those belong in Script 3 before-document)
  ✗ Individual patient outcomes (ecological fallacy)
  ✗ The exact pCR cutoff for clinical decision-making
    (this requires prospective validation beyond
    the scope of this analysis)
  ✗ Cross-subtype comparisons
    (those belong in the LumA vs TNBC cross-subtype
    before-document, not here)
  ✗ Epigenetic mechanism details beyond what is
    observable from bulk RNA-seq
    (H3K27me3 ChIP-seq needed for mechanism — not
    available in GSE25066 or TCGA-BRCA)
```

---

## PART VI — COMPLETE PREDICTION REFERENCE

```
PREDICTIONS LOCKED 2026-03-04
Document: BRCA-S4c
Author: Eric Robert Lawson, OrganismCore

S2-P1:  Gene-mapped depth score predicts pCR better
        than Script 1 proxy.
        S2-P1a: r(depth, pCR) more negative than -0.098
                in full GSE25066 cohort
        S2-P1b: r(depth, pCR) < -0.15 in TNBC subset

S2-P2:  EED expression predicts pCR in TNBC.
        S2-P2a: r(EED, pCR) < -0.10
        S2-P2b: AUC(EED) > AUC(EZH2) for pCR prediction

S2-P3:  PARP1 predicts chemotherapy resistance (low pCR).
        S2-P3a: r(PARP1, pCR) < -0.08 in TNBC
        S2-P3b (prospective): PARP1 predicts PARPi response
                (stated for future datasets)

S2-P4:  Lehmann subtypes map depth axis:
        LAR < BL1/BL2 < IM < M < MSL
        LAR mean depth < 0.40
        M/MSL mean depth > 0.60

S2-P5:  BRCA1 dysfunction enriched in Basal-like TCGA.
        S2-P5a: Mutation rate >20% Basal vs <5% LumA
        S2-P5b: Methylation rate >30% Basal vs <10% LumA
        S2-P5c: Combined dysfunction >50% in Basal-like

S2-P6:  Depth score predicts survival in TCGA-BRCA.
        S2-P6a: KM: depth-high shorter OS, p < 0.05
        S2-P6b: Cox HR > 1.3 per 0.1 depth increase

S2-P7:  SPI1 elevation is immune contamination.
        r(SPI1, PTPRC) > r(SPI1, IM_score)
        (contamination > genuine IM biology)

S2-P8:  EED > EZH2 as EZH2i response predictor.
        Prospective — stated for future trial datasets.

S2-P9:  AR correlates negatively with depth in bulk.
        r(AR, depth) < -0.15 in bulk TNBC

CONTROLS (should be flat in TNBC bulk):
  CDX2, NKX2-1, OLIG2, MBP
  Expected: near-absent, no depth correlation
  If elevated: flag and investigate
```

---

## STATUS BLOCK

```
document:           BRCA-S4c
type:               Before-Document (locked predictions)
date:               2026-03-04
author:             Eric Robert Lawson / OrganismCore
status:             LOCKED

predictions_count:  9 prediction groups
                    (S2-P1 through S2-P9)
                    14 specific sub-predictions

datasets:
  PRIMARY:   GSE25066 (n=508, pCR)
  SECONDARY: TCGA-BRCA (PAM50 Basal-like)

technical_prerequisites:
  GPL96 Affymetrix platform annotation file
  GSE25066 metadata format inspection before parsing
  Lehmann 2016 centroid classifier gene lists
  TCGA-BRCA PAM50 calls + BRCA1 mutation/methylation

next_document:      BRCA-S4d (Script 2 reasoning artifact)
script:             BRCA_TNBC_script2.py

series_rule:        No Script 2 data loads until this
                    document is committed to the repository.
```
