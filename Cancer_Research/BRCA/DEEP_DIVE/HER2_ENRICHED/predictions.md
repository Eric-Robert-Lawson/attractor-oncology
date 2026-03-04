# HER2-ENRICHED BREAST CANCER — BEFORE DOCUMENT
## Predictions Locked Before Script 1 Runs
## OrganismCore — Document BRCA-S3a
## Date: 2026-03-04

---

## DOCUMENT METADATA

```
document_id:        BRCA-S3a
series:             BRCA Deep Dive — HER2-Enriched
folder:             Cancer_Research/BRCA/DEEP_DIVE/HER2_ENRICHED/
type:               BEFORE-DOCUMENT
                    All predictions stated before
                    any data is loaded.
                    This document cannot be modified
                    after Script 1 runs.
date:               2026-03-04
author:             Eric Robert Lawson
                    OrganismCore
protocol_version:   Workflow_Protocol.md v2.0
dataset_planned:    GSE176078 — Wu et al. 2021
                    (same scRNA-seq dataset as LumA and TNBC —
                    Cancer Her2 SC cells extracted separately)
                    n=3,708 Cancer Her2 SC cells available
                    Reference: Mature Luminal (n=1,265)
                               Luminal Progenitors (n=1,992)
status:             LOCKED — predictions cannot change
                    after this document is committed
precursor_documents:
  ATTRACTOR_GEOMETRY_AXIOMS.md (Doc 90)
  BRCA_Subtypes.md (BRCA_Subtype_Orientation)
  BRCA-S2b (LumA Script 2 reasoning artifact)
  BRCA-S2c (LumA literature check)
  BRCA-S4b (TNBC Script 1 reasoning artifact)
```

---

## PART I — ATTRACTOR TYPE ASSIGNMENT
## Before any biology is stated

```
Per ATTRACTOR_GEOMETRY_AXIOMS.md (Doc 90),
the first step before any prediction is to
assign the attractor type using the diagnostic
algorithm.

QUESTION 1: Cell of origin?
  The normal cell is the LUMINAL PROGENITOR —
  a cell partway through luminal differentiation.
  The correct terminal destination is the mature
  luminal epithelial cell (ESR1+, FOXA1+, GATA3+).

QUESTION 2: Are identity TFs of the correct
  terminal state expressed?
  PARTIAL — HER2-enriched retains luminal
  structural markers (CK7, CK18) but is
  ER-NEGATIVE by definition (PAM50 criterion).
  ESR1 is absent or very low.
  FOXA1 may be partially retained.
  GATA3 is reduced.
  The luminal identity programme is INCOMPLETE —
  not fully expressed, not fully absent.

QUESTION 3: Is the cell structurally within
  the correct valley but the floor is removed?
  PARTIAL — the cell retains some luminal
  markers but the ER axis is off.
  The floor removal hypothesis does not fully fit
  because the cell is not near the terminal luminal
  state — it is arrested mid-differentiation.

QUESTION 4: Is the cell expressing identity TFs
  of a DIFFERENT normal cell type?
  NO — HER2-enriched does not adopt a basal,
  neural crest, or mesenchymal identity.
  KRT5/KRT14 are not elevated.
  SOX10, ZEB1, ZEB2, VIM are not the dominant signal.
  The cell has not fallen into a wrong valley.

ATTRACTOR TYPE DIAGNOSTIC CONCLUSION:

  This subtype does not fit cleanly into
  Type 1, Type 2, or Type 3 as originally defined.

  The HER2-enriched cell is a LUMINAL PROGENITOR
  that has been intercepted mid-differentiation
  by a constitutive signaling amplicon (ERBB2
  on 17q12). The PI3K/AKT/mTOR axis is
  constitutively activated. Proliferation is
  locked on. The cell cannot descend to the
  correct luminal terminal state because the
  ERBB2 amplicon is generating a force vector
  that holds it on the Waddington slope —
  unable to complete the descent.

PRIMARY TYPE ASSIGNMENT: TYPE 1 VARIANT —
  "SLOPE ARREST"

  The cell is not blocked at the top of the
  valley (classic Type 1 — blocked approach).
  The cell is not in a wrong valley (Type 2).
  The cell is not inside the correct valley
  with the floor removed (Type 3).

  The cell is ARRESTED ON THE SLOPE OF DESCENT
  toward the correct luminal terminal state.
  The Waddington landscape is intact.
  The ERBB2 amplicon is the external force
  preventing descent completion.

  PREDICTED GEOMETRIC SIGNATURE:
    — Partial luminal identity retained
      (CK7/CK18 present, FOXA1 partial)
    — ER axis absent (ESR1 very low)
    — HER2/ERBB2 massively elevated
    — GRB7 co-elevated (17q12 amplicon)
    — Proliferative scaffold elevated
      (MKI67, TOP2A, CCNB1)
    — Basal identity markers NOT elevated
      (KRT5, SOX10 should be near-normal)
    — Luminal switch genes PARTIALLY retained
      not fully suppressed as in TNBC

  THIS IS THE KEY STRUCTURAL DISTINCTION
  FROM TNBC:
    TNBC: switch genes near-completely ERASED
          (ESR1 -97%, PGR -98%)
    HER2: switch genes PARTIALLY RETAINED
          (some luminal programme still present)

  If this holds in the data, HER2-enriched
  confirms a new attractor geometry type not
  fully captured by the original three axioms.
  This would require an update to
  ATTRACTOR_GEOMETRY_AXIOMS.md after results
  are locked.
```

---

## PART II — CROSS-SUBTYPE CONTEXT
## What prior analyses established

```
From BRCA Deep Dive series:

LumA (BRCA-S2b):
  Type 3 geometry — correct valley, floor removed
  FOXA1/GATA3/ESR1 retained but GRADED
  Depth axis: CDKN1A / SMAD3 / TGF-β throttle
  EZH2 NOT the primary convergence node in LumA
  GATA3 as depth biomarker (novel)

TNBC (BRCA-S4b):
  Type 2 geometry — wrong valley
  All switch genes near-completely erased
  SOX10 +1323% — neural crest programme active
  EZH2 +270% — convergence node gate confirmed
  EED > EZH2 as depth driver (novel)
  pCR directionality confirmed

Cross-subtype prediction from BRCA-S4b:
  LumA and TNBC are geometrically OPPOSITE ends
  of the breast cancer landscape.
  HER2-enriched should occupy the MIDDLE —
  partial luminal retention, partial identity loss,
  and a completely different depth axis
  driven by copy number rather than
  transcription factor loss.
```

---

## PART III — GENE PANEL AND RATIONALE
## What to measure and why

```
CONFIRMED FALSE ATTRACTOR MARKERS (from prior subtypes):
  ESR1   — luminal switch gene (near-zero expected in HER2)
  FOXA1  — luminal switch gene (partial retention predicted)
  GATA3  — luminal switch gene (partial retention predicted)
  PGR    — luminal switch gene (near-zero expected in HER2)

HER2-SPECIFIC MARKERS:
  ERBB2  — the false attractor driver
           Expected: MASSIVELY elevated vs mature luminal
  GRB7   — co-amplified on 17q12
           Expected: elevated (co-amplicon signal)
  EGFR   — HER family member, cross-talk
           Expected: elevated (less so than in TNBC)

PROLIFERATIVE SCAFFOLD:
  MKI67  — proliferation index
           Expected: elevated (Grade 3 biology)
  TOP2A  — co-amplified in some HER2 tumors (17q21)
           Expected: elevated

EPIGENETIC AXIS:
  EZH2   — convergence node in TNBC
           Expected: elevated but LESS than in TNBC
           Rationale: EZH2 is a downstream target of
           PI3K/AKT signaling. If ERBB2 activates
           PI3K/AKT, EZH2 should be elevated.
           But the primary driver here is ERBB2 copy
           number — not epigenetic gate-keeping as
           in TNBC.
  EED    — PRC2 component (novel biomarker from TNBC)
           Expected: elevated but less than EZH2
           relative shift (opposite of TNBC)

PI3K AXIS:
  PIK3CA — mutated ~24% of HER2-enriched
           Expression proxy: not directly gene-level
           but downstream targets will show activation

BASAL/WRONG-VALLEY MARKERS (negative controls):
  KRT5   — basal keratin
           Expected: LOW (not in wrong valley)
  SOX10  — neural crest (TNBC dominant signal)
           Expected: LOW — KEY TEST
           If SOX10 is NOT elevated, this confirms
           HER2-enriched is NOT a Type 2 geometry.
  VIM    — mesenchymal/EMT
           Expected: LOW to moderate (not dominant)
  ZEB1   — EMT TF (dominant in TNBC)
           Expected: LOW

CROSS-CANCER CONTROLS:
  SPI1   — AML switch gene (should be low in breast)
  MBP    — GBM switch gene (should be low in breast)
  CDX2   — CRC switch gene (should be low in breast)
```

---

## PART IV — PREDICTIONS
## All stated 2026-03-04 before Script 1 runs

```
PREDICTION 1 — THE PARTIAL LUMINAL IDENTITY

  H1: FOXA1 is SUPPRESSED relative to Mature Luminal
      but NOT near-zero.
      Expected range: -30% to -60% suppression.
      TNBC was -81%. HER2 should be LESS suppressed.
      If FOXA1 in HER2 is > -70% suppression,
      prediction is confirmed.
      If FOXA1 in HER2 approaches TNBC levels (-81%),
      prediction is falsified.

  H2: GATA3 is SUPPRESSED relative to Mature Luminal
      but NOT near-zero.
      Expected: -20% to -50% suppression.
      TNBC was -53%. HER2 should be near or above that
      boundary — partial overlap with TNBC geometry
      is possible at the GATA3 level.

  H3: ESR1 is near-zero (ER-negative definition).
      Expected: >-85% suppression.
      This is NOT a distinguishing finding — it is
      expected by clinical definition.
      But the DEGREE of suppression relative to
      TNBC (-97%) is the question.
      Prediction: ESR1 suppression is LESS complete
      than in TNBC.

  H4: PGR is near-zero (PR-negative definition).
      Expected: >-85% suppression.
      Same logic as ESR1.
```

```
PREDICTION 2 — THE ERBB2 AMPLITUDE

  H5: ERBB2 is the LARGEST elevated gene in the panel.
      Expected: >+500% elevation vs Mature Luminal.
      This is the defining event — the copy number
      amplification should dominate the expression
      signal above all other elevated markers.

  H6: GRB7 is co-elevated with ERBB2.
      Expected: >+200% elevation.
      GRB7 sits on the 17q12 amplicon with ERBB2.
      Co-elevation confirms the copy number
      mechanism rather than transcriptional
      regulation alone.

  FALSIFICATION: If ERBB2 is NOT the largest
  elevated gene in the panel, and instead
  SOX10 or another basal marker dominates,
  the Type 1 Variant geometry is falsified
  and HER2-enriched may be closer to TNBC
  geometry than predicted.
```

```
PREDICTION 3 — THE BASAL MARKER TEST

  H7: SOX10 is NOT elevated in Cancer Her2 SC
      relative to Mature Luminal.
      Expected: SOX10 change < +50%.
      In TNBC: SOX10 was +1323%.
      If HER2-enriched shows SOX10 elevation
      above +200%, the cell has partially adopted
      a neural crest identity — Type 2 component
      is present and composite type evaluation
      is required.

  H8: KRT5 is NOT elevated.
      Expected: KRT5 change < +100%.
      In TNBC: KRT5 was +508%.
      KRT5 elevation would indicate basal identity
      programme activation.

  H9: ZEB1/VIM are NOT dominant.
      EMT markers should not be the top movers
      in HER2-enriched.
      In TNBC: ZEB1 +1024%, ZEB2 +1036%, VIM +370%.
      If VIM > +200% in HER2-enriched, there is
      an EMT component present.
```

```
PREDICTION 4 — EZH2 AXIS IN HER2

  H10: EZH2 IS elevated relative to Mature Luminal.
       Expected: +100% to +200% elevation.
       Rationale: ERBB2 → PI3K/AKT → EZH2 is a
       known signaling axis. HER2 amplification
       drives EZH2 indirectly via AKT.
       But the elevation should be LESS than
       in TNBC (+270%) because the primary
       driver is signaling, not epigenetic
       gate-keeping.

  H11: EZH2 elevation in HER2 is NOT correlated
       with depth in the same way as TNBC.
       In TNBC: EED r=+0.435 within TNBC cells.
       In HER2: EZH2 should show weaker
       within-subtype depth correlation because
       the depth axis is ERBB2-driven (copy number
       heterogeneity), not epigenetically driven.

  H12: EED/EZH2 ratio in HER2 is LOWER than in TNBC.
       In TNBC: EED was a stronger depth driver
       than EZH2.
       In HER2: this relationship should not hold —
       ERBB2 dominates the depth axis.
```

```
PREDICTION 5 — THE PROLIFERATIVE AXIS

  H13: MKI67 is elevated relative to Mature Luminal.
       Expected: >+100% elevation.
       HER2-enriched is almost always Grade 3.
       High proliferation is a defining feature.

  H14: TOP2A is elevated.
       Expected: >+150% elevation.
       TOP2A is frequently co-amplified at 17q21
       (near the ERBB2 amplicon) in HER2+ tumors.
       Co-elevation of TOP2A with ERBB2 would
       confirm the amplicon's reach beyond 17q12.
```

```
PREDICTION 6 — THE GEOMETRIC COMPARISON
               (CROSS-SUBTYPE ORDERING)

  H15: The suppression of FOXA1 in HER2 falls
       BETWEEN LumA and TNBC:
         LumA:  FOXA1 graded, partially retained
         HER2:  FOXA1 moderately suppressed
         TNBC:  FOXA1 -81% (near-absent)

  H16: The suppression of GATA3 in HER2 is
       SIMILAR TO OR SLIGHTLY GREATER THAN TNBC:
         LumA:  GATA3 retained with depth gradient
         HER2:  GATA3 ~-40% to -60%
         TNBC:  GATA3 -53%

  H17: The elevation of EZH2 in HER2 is
       LESS THAN TNBC:
         LumA:  EZH2 not a primary convergence node
         HER2:  EZH2 +100% to +200%
         TNBC:  EZH2 +270%

  H18: The elevation of ERBB2 in HER2 exceeds
       ANY single marker elevation seen in
       either LumA or TNBC.
       TNBC's largest signal: SOX10 +1323%.
       Prediction: ERBB2 in HER2 exceeds +500%.
       If ERBB2 > SOX10 (in absolute % terms
       within their respective subtypes),
       this confirms that copy number amplification
       produces a qualitatively different class
       of attractor force than transcription
       factor re-expression.
```

```
PREDICTION 7 — LUMINAL PROGENITOR AS REFERENCE

  H19: Cancer Her2 SC is CLOSER to Luminal
       Progenitors than Cancer Basal SC (TNBC) is.
       PCA geometry prediction:
         TNBC should be maximally separated from
         Luminal Progenitors on PC1.
         HER2 should occupy intermediate PCA space —
         further from Mature Luminal than LumA,
         but NOT as far as TNBC.
       This is the geometric test of the
       "slope arrest" hypothesis.

  H20: The PCA distance from Cancer Her2 SC
       to Mature Luminal is LESS than the
       PCA distance from Cancer Basal SC
       to Mature Luminal.
       Quantitative prediction:
         If TNBC centroid-to-Luminal distance = 1.0,
         HER2 centroid-to-Luminal distance < 0.7.
```

---

## PART V — WHAT SCRIPT 1 DOES NOT PREDICT

```
Script 1 uses the SAME scRNA-seq dataset
(GSE176078 — Wu et al. 2021) used for LumA and TNBC.
The data is already downloaded and verified.

What Script 1 CANNOT test:
  - Bulk RNA validation (requires TCGA-BRCA Script 2)
  - Clinical outcome correlations (Script 2)
  - pCR rates in HER2-targeted therapy (not in scRNA-seq)
  - ERBB2 copy number (expression only, not genomic)
  - GRB7 — only if present in the gene matrix
    (single-cell dropout is high for low-expression genes)
  - TP53 mutation status (expression proxy only)
  - Treatment response (this dataset has no treatment data)

What Script 1 CAN test:
  - Cell state geometry (Cancer Her2 SC vs Mature Luminal)
  - Attractor type confirmation or falsification
  - ERBB2 expression amplitude (copy number → expression)
  - Partial vs complete luminal identity loss
  - SOX10/KRT5 test (basal identity absent or present)
  - EZH2 elevation and within-subtype depth correlation
  - PCA geometry (slope arrest vs wrong valley vs floor removed)
  - Cross-subtype comparison (LumA vs HER2 vs TNBC)
    using the same dataset and reference populations
```

---

## PART VI — FALSIFICATION CRITERIA

```
The TYPE 1 VARIANT "slope arrest" geometry is FALSIFIED if:

  FALSIFICATION 1:
    SOX10 elevation in Cancer Her2 SC > +500%.
    Interpretation: HER2-enriched has adopted
    neural crest identity — Type 2 geometry present.
    This would mean HER2-enriched is a composite
    Type 1 → Type 2, analogous to TNBC.

  FALSIFICATION 2:
    FOXA1 suppression in Cancer Her2 SC > -75%.
    Interpretation: Luminal identity is near-erased,
    not partially retained. The slope arrest
    geometry prediction fails.

  FALSIFICATION 3:
    ERBB2 is NOT elevated or is elevated < +100%.
    Interpretation: The copy number event is not
    detectable at the expression level in this
    scRNA-seq dataset (possible due to dropout).
    In this case: Script 1 result is inconclusive,
    not falsified — require TCGA bulk confirmation.

  FALSIFICATION 4:
    EZH2 elevation in HER2 EXCEEDS TNBC (+270%).
    Interpretation: EZH2 is a more dominant
    convergence node in HER2 than in TNBC —
    the epigenetic gate-keeping hypothesis
    is the primary driver, not ERBB2.
    This would require major revision of the
    geometric model.

  PARTIAL FALSIFICATION:
    If predictions H15–H18 (cross-subtype ordering)
    are violated — specifically if HER2-enriched
    falls geometrically CLOSER to TNBC than to LumA —
    the slope arrest framing is incorrect and
    HER2-enriched shares more with TNBC than predicted.
```

---

## PART VII — WHAT SCRIPT 1 NEEDS TO RUN

```
DATA REQUIRED (all already available locally):
  scRNA-seq: GSE176078 Wu_etal_2021_BRCA_scRNASeq/
    count_matrix_sparse.mtx     — confirmed OK
    count_matrix_genes.tsv      — confirmed OK
    count_matrix_barcodes.tsv   — confirmed OK
    metadata.csv                — confirmed OK

CELL TYPES TO EXTRACT:
  TUMOR:      Cancer Her2 SC        (n=3,708)
  REFERENCE 1: Mature Luminal       (n=1,265)
  REFERENCE 2: Luminal Progenitors  (n=1,992)
  COMPARISON:  Cancer Basal SC      (n=4,312)   [TNBC]
               Cancer LumA SC       (n=7,742)   [LumA]

  The TNBC and LumA cells are included in Script 1
  for the cross-subtype comparison in PART VI.
  This is the first time all three subtypes are
  analyzed in a single script against the same
  reference population in the same dataset.

GENE PANEL:
  ERBB2, GRB7
  FOXA1, GATA3, ESR1, PGR
  EZH2, EED
  KRT5, SOX10, VIM, ZEB1
  MKI67, TOP2A
  EGFR, AR
  SPI1, MBP, CDX2  (cross-cancer controls)
```

---

## STATUS BLOCK

```
status:         LOCKED
                This document cannot be modified
                after Script 1 runs.
next_document:  BRCA-S3b
                (HER2-Enriched Script 1 reasoning artifact)
dataset:        GSE176078 — no download required
script:         BRCA_HER2_script1.py
date_locked:    2026-03-04
```
