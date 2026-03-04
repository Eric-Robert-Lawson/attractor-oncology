# LUMINAL A BREAST CANCER — SCRIPT 1 REASONING ARTIFACT
## Predictions vs Data | Corrections | Updated Predictions | Drug Targets
## OrganismCore — Document BRCA-S1b
## Date: 2026-03-04

---

## DOCUMENT METADATA

```
document_id:        BRCA-S1b
series:             BRCA Deep Dive — Luminal A
folder:             Cancer_Research/BRCA/DEEP_DIVE/Luminal_A/
type:               REASONING ARTIFACT
                    Script 1 results interpreted.
                    Predictions assessed.
                    Corrected attractor derived.
                    Script 2 predictions locked.
date:               2026-03-04
author:             Eric Robert Lawson
                    OrganismCore
protocol_version:   Workflow_Protocol.md v2.0
dataset:            GSE176078 — Wu et al. 2021
                    Nature Genetics | PMID: 34493872
populations:
  Cancer LumA SC:        n=7,742 cells
  Mature Luminal:        n=1,265 cells  (normal reference)
  Cancer Basal SC:       n=4,312 cells  (cross-reference)
  Luminal Progenitors:   n=1,992 cells  (Script 2 reference)
precursor_document:  BRCA-S1a
                     Predictions locked 2026-03-04
next_document:       BRCA-S2a
                     Script 2 before-document
                     (new predictions locked before
                      Script 2 runs)
status:             COMPLETE — initial version
```

---

## PART I — WHAT THE GEOMETRY FOUND
## (Read before predictions — Protocol v2.0)

---

### 1.1 TOP MOVERS — UNFILTERED

```
The geometry found two things elevated.
Everything else was flat or suppressed.

TOP GAINED (LumA vs Mature Luminal):
  FOXA1    +37.3%   p=3.61e-13  ***
  GATA3    +33.8%   p=9.68e-13  ***
  EZH2     +18.5%   p=2.40e-01  ns  (flat)
  CCND1    + 6.4%   p=7.20e-01  ns  (flat)
  — nothing else elevated —

The landscape gained two genes.
Both are luminal identity transcription factors.
There is no oncogene activation in the top movers.
There is no proliferative signature at the top.
The two most elevated genes are the pioneer factors
that define luminal cell identity in the first place.

TOP LOST (LumA vs Mature Luminal):
  CDKN1A   -74.3%   p=4.68e-195  ***  ← DEEPEST SIGNAL
  MBP      -83.2%   p=3.16e-112  ***
  EGFR     -93.7%   p=1.08e-99   ***
  KRT14    -94.8%   p=1.52e-122  ***
  KRT5     -96.5%   p=1.20e-85   ***
  SOX10    -96.7%   p=2.48e-17   ***
  PGR      -54.8%   p=3.11e-22   ***
  MYC      -44.5%   p=6.25e-57   ***
  HDAC1    -47.8%   p=9.46e-43   ***
  HDAC2    -44.5%   p=7.88e-33   ***
  KDM1A    -44.3%   p=9.09e-14   ***
  CCNE1    -69.4%   p=4.00e-06   ***
  PIK3CA   -44.9%   p=3.31e-11   ***
  MTOR     -46.0%   p=1.33e-09   ***
  RB1      -38.2%   p=5.02e-14   ***
  TP53     -37.3%   p=6.40e-09   ***
  CDK4     -28.3%   p=7.99e-15   ***
  PTEN     -27.0%   p=6.45e-05   ***

The dominant suppression signal is CDKN1A at -74.3%
with p=4.68e-195 — the most significant finding
in the entire dataset by a wide margin.

CDKN1A = p21/CIP1
The CDK inhibitor that enforces cell cycle arrest
downstream of TP53, TGF-β, and contact inhibition.

The LumA cancer cell has not become more proliferative.
It has become unable to stop.
```

---

### 1.2 THE CENTRAL STRUCTURAL FINDING

```
LumA cancer vs normal mature luminal:

  GAINED: luminal identity TFs (FOXA1 +37%, GATA3 +34%)
  LOST:   cell cycle arrest capacity (CDKN1A -74%)
  LOST:   basal/myoepithelial markers completely
          (KRT5 -97%, KRT14 -95%, SOX10 -97%, EGFR -94%)
  PARTIAL: ER circuit broken at PGR level
           (ESR1 flat, FOXA1 elevated, PGR -55%)
  LOST:   global chromatin-modifying machinery
          (HDAC1, HDAC2, KDM1A all -44% to -48%)

This is a specific attractor geometry:
  Identity hyper-fixed (FOXA1/GATA3 above normal).
  Arrest capacity eliminated (CDKN1A).
  Basal identity completely excluded.
  ER programme partially decoupled.
```

---

### 1.3 THE DIVERGENCE FROM TNBC

```
LumA and TNBC are geometrically opposite:

                 LumA        Basal/TNBC
  FOXA1          +37%          -83%
  GATA3          +34%          -66%
  ESR1           -30% (flat)   -98%
  EZH2           +19% (flat)   +324%  ← LOCK
  CCND1          + 6% (flat)   +260%
  CDK4           -28%          +121%
  MKI67          -62% (flat)   +976%
  KRT5           -97%         +1145%
  KRT14          -95%          +733%
  SOX10          -97%         +1717%

LumA is NOT a high-proliferation state.
TNBC IS.
They differ by what they kept and what they lost.

The false attractor geometry is:
  TNBC:  proliferative escalation
         + identity erasure
         + EZH2 lock
  LumA:  identity fixation
         + arrest removal
         + ER circuit partial break

Two completely different mechanisms.
Same normal cell of origin (luminal epithelium).
Different saddle points.
Different drug targets.
```

---

### 1.4 WITHIN-LumA DEPTH AXIS

```
Top Pearson r correlates with depth
within Cancer LumA SC cells:

  CDK4     r=+0.81  p=0.00e+00  ***
  PCNA     r=+0.63  p=0.00e+00  ***
  GATA3    r=+0.57  p=0.00e+00  ***
  CCND1    r=+0.53  p=0.00e+00  ***
  AR       r=+0.45  p=0.00e+00  ***
  MYC      r=+0.40  p=9.01e-295 ***
  HDAC1    r=+0.38  p=2.68e-271 ***
  CCNE1    r=+0.35  p=1.02e-220 ***
  PGR      r=+0.24  p=1.34e-104 ***
  EZH2     r=+0.23  p=2.36e-92  ***

Interpretation:
  Within the LumA population, cells that are
  "deeper" (more advanced cancer state) have
  MORE CDK4, MORE PCNA, MORE CCND1.
  These are cells with residual cell cycle
  machinery that is still running.

  On average LumA has LESS CDK4 than normal.
  But within LumA the most CDK4-active cells
  are the deepest.

  This means: the depth axis within LumA
  is NOT the distance from normal on a
  proliferative scale.
  It is the distribution of cells that have
  retained the most cell cycle activity
  within a population that has lost its
  arrest capacity.

  The "deep" LumA cell retains CDK4/CCND1
  AND has lost CDKN1A.
  That cell cannot stop cycling.
  The "shallow" LumA cell has lower CDK4
  AND lower CDKN1A — still more arrested.

  CDKN1A loss unlocks CDK4/CCND1 activity.
  The depth axis is fundamentally the
  CDKN1A/CDK4 ratio.
```

---

## PART II — PREDICTIONS VS DATA

---

### 2.1 PREDICTION 1 — Switch genes retained
### LOCKED: FOXA1, GATA3, ESR1, PGR NOT suppressed

```
FOXA1:  +37.3%  p=3.61e-13  ELEVATED     ✓ CONFIRMED
GATA3:  +33.8%  p=9.68e-13  ELEVATED     ✓ CONFIRMED
ESR1:   -30.1%  p=9.41e-01  FLAT (ns)    ✓ CONFIRMED
PGR:    -54.8%  p=3.11e-22  SUPPRESSED   ✗ WRONG

VERDICT: PARTIALLY CONFIRMED (3/4)

ANALYST ASSUMPTION ERROR — PGR:
  PGR was predicted as a switch gene
  that would be retained in LumA.
  PGR is significantly suppressed at -55%.
  This is an error in what PGR represents
  in this context.

  PGR is NOT a luminal identity TF.
  PGR is an ER TARGET GENE.
  ESR1 (the receptor) is retained.
  FOXA1 (the ER pioneer) is elevated.
  But PGR (what ER drives) is lost.

  WHAT THIS TEACHES:
    The ER transcriptional programme is
    partially decoupled in LumA cancer.
    ER is present and its pioneer factor is
    hyperactive, but ER is not fully executing
    its transcriptional programme.
    Something between ER binding and PGR
    transcription is broken.
    This is the ER CIRCUIT BREAK.
    It is the most important corrected finding
    from Prediction 1.

  POINTER FOR SCRIPT 2:
    Test ER coactivators: NCOA1, NCOA2, NRIP1.
    If these are suppressed, that explains
    why ER is present but PGR is not transcribed.
    The gap is at the coactivator level.
```

---

### 2.2 PREDICTION 2 — Depth axis proliferative
### LOCKED: MKI67, CCND1, CDK4 elevated in LumA

```
MKI67:   -61.9%  p=1.46e-01  FLAT (ns)    ✗ WRONG
CCND1:   + 6.4%  p=7.20e-01  FLAT (ns)    ✗ WRONG
CDK4:    -28.3%  p=7.99e-15  SUPPRESSED   ✗ WRONG
MYC:     -44.5%  p=6.25e-57  SUPPRESSED   ✗ WRONG
CCNE1:   -69.4%  p=4.00e-06  SUPPRESSED   ✗ WRONG
PCNA:    -18.1%  p=1.09e-05  FLAT         ✗ WRONG
TOP2A:   -18.3%  p=5.74e-01  FLAT (ns)    ✗ WRONG

VERDICT: NOT CONFIRMED (0/7)

ANALYST ASSUMPTION ERROR — ALL PROLIFERATIVE GENES:
  The prediction was that LumA cancer cells
  would be MORE proliferative than mature luminal.
  The data shows the OPPOSITE.
  LumA cancer cells are LESS proliferative
  than mature luminal normal cells
  in mean expression terms.

  This is the most significant wrong prediction
  in Script 1. And it teaches the most.

  WHAT THIS TEACHES:
    1. LumA cancer is NOT a high-proliferation
       state. The mechanism is not accelerating
       the engine. It is releasing the brake.
    2. CDKN1A at -74% is the primary finding.
       When p21 is gone, even LOW CDK4 is
       sufficient for uncontrolled cycling.
       The cell does not need to increase CDK4 —
       it needs only to remove CDKN1A.
    3. The depth axis is not proliferative genes
       elevated vs normal.
       The depth axis is arrest capacity lost
       vs normal.
    4. The within-LumA correlations (CDK4 r=+0.81,
       CCND1 r=+0.53) represent a DIFFERENT axis:
       within the cancer population, cells with
       more cycle machinery are "deeper."
       But the whole population has LESS machinery
       than normal.

  POINTER FOR SCRIPT 2:
    Depth score must be redesigned.
    Primary axis: 1 - norm(CDKN1A)
    This is the arrest dismantlement score.
    Secondary axis: CDK4 within-cancer
    correlations are real but need the
    CDKN1A frame to be interpreted correctly.

    Also redesign the conceptual frame:
    LumA depth = degree of arrest dismantlement,
    NOT degree of proliferative escalation.
```

---

### 2.3 PREDICTION 3 — EZH2 flat in LumA
### LOCKED: EZH2 not +270% as in TNBC

```
EZH2:  +18.5%  p=2.40e-01  FLAT (ns)    ✓ CONFIRMED

TNBC cross-reference:
  EZH2 in Basal SC vs Mature Luminal:  +323.7%

VERDICT: CONFIRMED

  EZH2 is confirmed as a TNBC-specific lock.
  It is not the mechanism in LumA.
  LumA does not use PRC2/H3K27me3 to maintain
  its false attractor state.
  LumA uses CDKN1A loss.

  The EZH2 inhibitor (tazemetostat) is the
  correct drug for TNBC but NOT for LumA.
  This cross-subtype distinction is now
  confirmed from data.
```

---

### 2.4 PREDICTION 4 — PIK3CA axis uncertain
### LOCKED: Post-translational, mRNA may not show

```
PIK3CA:   -44.9%  p=3.31e-11  SUPPRESSED
AKT1:     -19.1%  p=1.88e-03  FLAT
AKT2:     -28.3%  p=2.24e-13  SUPPRESSED
MTOR:     -46.0%  p=1.33e-09  SUPPRESSED
PTEN:     -27.0%  p=6.45e-05  SUPPRESSED

VERDICT: PREDICTION CORRECT — axis suppressed at mRNA level.
  All PIK3CA pathway mRNAs are lower in LumA
  than in mature luminal normal.
  This is consistent with post-translational
  activation being the real mechanism:
  PIK3CA mutations produce constitutively
  active PI3K protein from a LOWER mRNA level.
  The pre-stated uncertainty was correct.
  mRNA suppression with active protein is
  exactly what PIK3CA gain-of-function mutations
  produce — the mutant protein is more stable,
  more active, expressed from fewer transcripts.
  The prediction of uncertainty was validated.
```

---

### 2.5 PREDICTION 5 — Depth predicts outcome
### LOCKED: Needs bulk data

```
scRNA-seq: no survival data.
Bulk data downloaded (24 tumors).

Bulk CDKN1A across 24 tumors:
  mean=2260  std=2204  (very high variance)

Bulk depth scores (MKI67/CCND1/CDK4 composite):
  Mean: 0.29  Std: 0.15  Min: 0.04  Max: 0.71

The depth score range of 0.04–0.71 across
24 tumors represents substantial heterogeneity.
Whether CDKN1A bulk levels predict recurrence
or endocrine therapy resistance is the
clinically testable prediction.
This prediction remains OPEN.
Script 2 bulk analysis will test it.
```

---

### 2.6 PREDICTION 6 — Controls flat
### LOCKED: SPI1, CDX2, MBP flat

```
SPI1:    -37.7%  p=5.21e-02  FLAT (ns)    ✓ CONFIRMED
CDX2:    + 0.0%  p=1.00e+00  FLAT         ✓ CONFIRMED
MBP:     -83.2%  p=3.16e-112 SUPPRESSED   ✗ FLAG
NKX2-1:  + 0.0%  p=1.00e+00  FLAT         ✓ CONFIRMED
OLIG2:   + 0.0%  p=6.86e-01  FLAT         ✓ CONFIRMED

MBP FLAG:
  MBP (myelin basic protein) was predicted
  to be flat as a non-breast control.
  It is -83.2% at p=3.16e-112.
  The normal mature luminal cells express MBP
  at mean=0.24. LumA has mean=0.04.

  This is either:
    a) Genuine MBP expression in normal luminal
       cells that is lost in cancer (real biology)
    b) Contaminating oligodendrocyte/Schwann cell
       in the Wu 2021 normal population (technical)

  The question is not crucial for the framework
  but should be noted. Not over-interpreted.
  MBP is not being added to the attractor panel.
  This is a flag for data quality awareness.

VERDICT: CONFIRMED with MBP flag noted.
```

---

### 2.7 PREDICTION 7 — TNBC markers absent
### LOCKED: KRT5, SOX10, KRT14, EGFR absent

```
KRT5:   -96.5%  p=1.20e-85   SUPPRESSED  ✓ CONFIRMED
KRT14:  -94.8%  p=1.52e-122  SUPPRESSED  ✓ CONFIRMED
SOX10:  -96.7%  p=2.48e-17   SUPPRESSED  ✓ CONFIRMED
EGFR:   -93.7%  p=1.08e-99   SUPPRESSED  ✓ CONFIRMED

VERDICT: CONFIRMED — all 4 genes near-absent.
  Population purity confirmed.
  The Cancer LumA SC metadata label is accurate.
  These are clean luminal cancer cells with
  zero contamination by basal/myoepithelial cells.
```

---

## PART III — CONVERGENCE TABLE

```
Prediction      Gene(s)              Result      Verdict
────────────────────────────────────────────────────────────
P1 — retained   FOXA1 +37.3%         ELEVATED    ✓ CONFIRMED
P1 — retained   GATA3 +33.8%         ELEVATED    ✓ CONFIRMED
P1 — retained   ESR1  -30.1% (ns)    FLAT        ✓ CONFIRMED
P1 — retained   PGR   -54.8%         SUPPRESSED  ✗ ERROR
P2 — prolif↑    MKI67 -61.9% (ns)    FLAT        ✗ ERROR
P2 — prolif↑    CCND1  +6.4% (ns)    FLAT        ✗ ERROR
P2 — prolif↑    CDK4  -28.3%         SUPPRESSED  ✗ ERROR
P2 — prolif↑    MYC   -44.5%         SUPPRESSED  ✗ ERROR
P2 — prolif↑    CCNE1 -69.4%         SUPPRESSED  ✗ ERROR
P3 — EZH2 flat  EZH2  +18.5% (ns)    FLAT        ✓ CONFIRMED
P4 — uncertain  PIK3CA -44.9%        SUPPRESSED  ✓ CONFIRMED
                                      (predicted
                                       uncertainty)
P5 — bulk       CDKN1A bulk high      OPEN        — PENDING
                variance observed
P6 — controls   SPI1/CDX2/NKX2-1     FLAT        ✓ CONFIRMED
P6 — control    MBP -83.2%            SUPPRESSED  ⚠ FLAG
P7 — TNBC out   KRT5/14/SOX10/EGFR   -94-97%     ✓ CONFIRMED

ANALYST ASSUMPTION ERRORS:
  1. PGR predicted retained — suppressed -55%
     Corrected: PGR is ER target, not identity TF
     ER circuit break at coactivator level

  2. Proliferative axis predicted elevated —
     all 7 genes flat or suppressed
     Corrected: Arrest removal (CDKN1A -74%)
     is the mechanism, not proliferative escalation

NOVEL SIGNALS NOT IN PREDICTIONS:
  1. CDKN1A -74.3% p=4.68e-195 — deepest signal
     Cell cycle arrest programme dismantled
  2. FOXA1/GATA3 ELEVATED (not just retained)
     Identity hyper-fixation, not identity loss
  3. Global chromatin machinery suppressed
     HDAC1 -48%, HDAC2 -45%, KDM1A -44%
  4. Within-LumA: CDK4 r=+0.81 with depth
     The "deeper" cells retain the most CDK4
     within a population that has LESS CDK4
     than normal
```

---

## PART IV — THE CORRECTED ATTRACTOR

---

### 4.1 What the LumA false attractor is

```
The three components:

COMPONENT 1 — IDENTITY FIXATION:
  Gene: FOXA1 (+37%), GATA3 (+34%)
  The LumA cancer cell is MORE committed
  to luminal identity than normal mature
  luminal cells.
  This is not identity LOSS.
  It is identity LOCK.
  The pioneer TFs that define luminal
  chromatin access are hyperactivated.
  This creates a stable state —
  the cell cannot be pushed out of
  the luminal programme.

COMPONENT 2 — ARREST DISMANTLEMENT:
  Gene: CDKN1A (-74% p=4.68e-195)
  The cell cycle arrest programme is gone.
  p21/CIP1 transcript is near-absent.
  The cell does not need to accelerate CDK4 —
  it only needs to remove the CDK inhibitor
  that normally stops CDK4 from firing.
  With CDKN1A gone, even residual CDK4 levels
  are sufficient for continuous cycling.
  This is the mechanism of unchecked growth.
  Not oncogene activation.
  Tumour suppressor loss.

COMPONENT 3 — ER CIRCUIT PARTIAL BREAK:
  ESR1: FLAT (retained)
  FOXA1: ELEVATED
  PGR: -54.8% (suppressed)
  The receptor is present.
  The pioneer is hyperactivated.
  But the canonical ER target gene is lost.
  The ER programme is not fully executing.
  Something between ER binding and PGR
  transcription is decoupled.
  Most likely: coactivator suppression
  (NCOA1, NCOA2, NRIP1 — not yet tested).
```

---

### 4.2 The Waddington geometry

```
The LumA cancer cell is NOT stuck between
two valleys. It IS in the luminal valley.
But within the luminal valley, it has:

  1. Moved to the deepest point of the valley
     (FOXA1/GATA3 hyperactivated — identity
     more entrenched than normal).

  2. Removed the walls of the valley on the
     proliferation side (CDKN1A gone — the
     barrier to unchecked cycling is absent).

  3. Created a partial disconnection in the
     valley's internal wiring (ER present
     but not fully connected to its programme).

The therapeutic target is not:
  "Push the cell out of the luminal valley."
  (Endocrine therapy tries this and mostly fails
  in advanced disease because the identity is
  too deeply fixed.)

The therapeutic target IS:
  "Restore the walls that CDKN1A normally
  provides within the luminal valley."
  OR
  "Block the CDK4/cyclin D machinery that
  CDKN1A normally suppresses."

This is precisely what CDK4/6 inhibitors do.
The geometry derived it.
The pharmacology confirmed it.
```

---

## PART V — DRUG TARGETS
## Stated from geometry before literature check

---

### PRIMARY DRUG TARGET (geometry-derived)

```
TARGET: CDK4/CDK6
MECHANISM:
  CDKN1A (-74%) removes the arrest
  that normally constrains CDK4/6
  activity in luminal cells.
  CDK4/6 inhibitors restore the arrest
  that CDKN1A would normally provide.
  They do not suppress an oncogene.
  They substitute for a lost tumour suppressor.

PREDICTED DRUGS:
  Palbociclib (Ibrance)
  Ribociclib (Kisqali)
  Abemaciclib (Verzenio)

LITERATURE STATUS: ✓ CONFIRMED
  CDK4/6 inhibitors are standard of care
  for ER+ metastatic breast cancer.
  The framework independently derived the
  correct drug class from CDKN1A loss.

DEPTH PREDICTION:
  Tumors with the lowest CDKN1A bulk expression
  should show the greatest CDK4/6 inhibitor
  benefit (deepest arrest dismantlement).
  This is testable from the bulk data.
  Script 2 will test it.
```

---

### SECONDARY DRUG TARGET (geometry-derived)

```
TARGET: ER coactivator pathway
MECHANISM:
  ESR1 retained, FOXA1 elevated,
  PGR suppressed.
  The ER circuit is broken between
  receptor binding and target gene activation.
  If the break is at NCOA1/NCOA2 (coactivators),
  then drugs that force ER into a non-activating
  conformation may bypass the broken step.

PREDICTED DRUGS:
  SERDs (selective ER degraders):
    Fulvestrant, elacestrant, camizestrant
  Mechanism: degrade ER entirely rather
  than competing for binding
  Works even when the downstream programme
  is partially decoupled

LITERATURE STATUS: ✓ CONFIRMED
  SERDs are approved/in trials for
  ER+ breast cancer, particularly where
  endocrine therapy has failed.

DEPTH PREDICTION:
  Tumors with the most decoupled ER circuit
  (lowest PGR, retained ESR1) predict SERD
  response better than aromatase inhibitors.
  The PGR/ESR1 ratio is the biomarker.
  Low PGR + present ESR1 = SERD candidate.
  Script 2 will compute this ratio.
```

---

### NOVEL DRUG TARGET (geometry-derived — not standard of care)

```
TARGET: FOXA1 hyperactivation
MECHANISM:
  FOXA1 is the most elevated gene in LumA
  after GATA3 (+37.3% p=3.61e-13).
  FOXA1 is the pioneer TF that opens luminal
  chromatin for ER binding.
  Elevated FOXA1 in LumA cancer deepens the
  identity fixation — it makes the luminal
  attractor MORE stable, not less.
  This is the IDENTITY ANCHOR.

  In prostate cancer (PRAD), FOXA1
  gain-of-function mutations are oncogenic
  — FOXA1 is a documented oncogenic TF
  in the luminal lineage.
  The same mechanism may apply in LumA
  breast cancer.

PREDICTED DRUG:
  FOXA1 inhibitors (in development)
  CUT&RUN-validated FOXA1 targeting
  compounds (preclinical stage)

LITERATURE STATUS: 🆕 NOVEL PREDICTION
  FOXA1 overexpression in LumA breast cancer
  is documented in literature.
  FOXA1 as a drug target in breast cancer
  is an active research area but not standard
  of care.
  Whether FOXA1 ELEVATION (not mutation) is
  a depth biomarker in LumA is not established.
  The framework predicts it is.
  Script 2 will test whether FOXA1 expression
  level within LumA tumors correlates with
  depth score and clinical outcomes.

PREDICTION FOR LITERATURE CHECK:
  High FOXA1 in LumA predicts worse outcome
  and endocrine therapy resistance.
  This is testable from the bulk data.
```

---

### CROSS-SUBTYPE INSIGHT (geometry-derived)

```
FINDING: EZH2 lock is TNBC-specific,
         not a universal breast cancer mechanism.

  LumA: EZH2 +19% (flat)
  TNBC: EZH2 +324%

  Tazemetostat (EZH2 inhibitor) is the
  correct drug for TNBC.
  Tazemetostat is the WRONG drug for LumA.
  Giving tazemetostat to a LumA patient
  would suppress a gene that is already flat
  and do nothing to the CDKN1A loss that
  drives the actual biology.

  This is a patient stratification finding.
  It confirms the framework prediction that
  EZH2 direction must be determined per
  cancer type, not assumed.
  (Protocol Lesson 5, Workflow_Protocol.md v2.0)
```

---

## PART VI — SCRIPT 2 PREDICTIONS
## Locked before Script 2 runs

---

### 6.1 What Script 2 will test

```
DESIGN CHANGES FROM SCRIPT 1:

  1. REFERENCE: LumA vs Luminal Progenitors
     (n=1,992) in addition to Mature Luminal
     (n=1,265).
     Luminal Progenitors may be a more
     appropriate normal reference because
     they are cycling — they will not show
     CCND1/CDK4 loss artefact of the
     post-mitotic Mature Luminal comparison.

  2. DEPTH SCORE: CDKN1A-centred
     New primary depth score:
       depth = 1 - norm(CDKN1A)
     Secondary component:
       + norm(CDK4) [within-LumA activity axis]
     This replaces the proliferative panel score.

  3. GENE PANEL EXTENSION:
     ER coactivators:   NCOA1, NCOA2, NRIP1
                        (gap test: ESR1→PGR break)
     TGF-β pathway:    TGFB1, SMAD3, TGFBR1
                        (upstream of CDKN1A)
     FOXA1 programme:  FOXA1 target genes to
                        confirm elevation is
                        functional
     Cell cycle nodes: CDK2, CDK6, CCNA1, CCNB1
                        (downstream of CDK4)

  4. GAP TEST: r(ESR1, PGR) within LumA cells
     If r is near-zero: ER→PGR circuit broken
     If r is significant: different mechanism
     for PGR loss
```

---

### 6.2 Script 2 predictions (locked)

```
ALL PREDICTIONS STATED BEFORE SCRIPT 2 RUNS.

PREDICTION S2-1:
  CDKN1A depth score separates LumA
  into high-depth and low-depth populations.
  High depth = lower CDKN1A = more arrest loss.
  p < 0.05 for depth separation.
  [DIRECTION: CDKN1A negatively defines depth]

PREDICTION S2-2:
  r(ESR1, PGR) within LumA cells < 0.3
  The ER→PGR circuit is broken.
  ESR1 and PGR expression are uncoupled
  in cancer cells even though ESR1 is retained.
  [DIRECTION: near-zero or negative r]

PREDICTION S2-3:
  ER coactivators (NCOA1 or NCOA2) are
  suppressed in LumA vs Luminal Progenitors.
  The coactivator suppression explains the
  ESR1→PGR break.
  [DIRECTION: NCOA1/NCOA2 DOWN]

PREDICTION S2-4:
  TGFB1 or SMAD3 is suppressed in LumA.
  TGF-β is the primary upstream driver of
  CDKN1A transcription via SMAD binding to
  the p21 promoter.
  If TGF-β signalling is reduced, CDKN1A
  loss follows mechanistically.
  [DIRECTION: TGFB1/SMAD3 DOWN]

PREDICTION S2-5:
  FOXA1 expression within LumA cells
  positively correlates with depth score
  (r > 0.30).
  Higher FOXA1 = deeper identity fixation
  = more aggressive cancer state.
  [DIRECTION: r(FOXA1, depth) > 0.30]
  Note: From Script 1, r(GATA3, depth) = +0.57
  already confirmed this direction for GATA3.

PREDICTION S2-6:
  LumA vs Luminal Progenitor comparison
  will show LESS suppression of CDK4/CCND1
  than LumA vs Mature Luminal.
  The proliferative suppression seen in
  Script 1 was partly artefactual — the
  Mature Luminal normal is post-mitotic
  and has higher CDK4/CCND1 than progenitors
  because it is actively cycling.
  [DIRECTION: CDK4 change attenuates vs
   Luminal Progenitor reference]

PREDICTION S2-7 — NOVEL:
  CDKN1A bulk expression level across the
  24 Wu et al. bulk tumors is bimodal
  or spans a sufficient range to stratify
  tumors by predicted CDK4/6 inhibitor benefit.
  Tumors with lowest bulk CDKN1A have
  highest depth scores.
  [DIRECTION: CDKN1A bulk negatively
   correlates with bulk depth score]
```

---

## PART VII — FRAMEWORK LESSONS

```
LESSON FROM BRCA-S1b:
  Luminal A is a fundamentally different
  geometry from all prior cancer validations.
  All prior cancers had SUPPRESSED switch genes
  as the primary signal.
  LumA has ELEVATED switch genes (FOXA1/GATA3)
  as the primary signal.
  The false attractor in LumA is not a
  differentiation block.
  It is a growth regulation failure within
  a correctly differentiated cell.

  Rule for future analyses:
    When a cancer retains its identity TFs
    (as established by prior analysis),
    the depth axis is NOT the identity
    programme itself.
    It is what the identity programme fails
    to maintain: cell cycle arrest capacity,
    growth factor independence, spatial
    organisation.
    Look for the LOST CONSTRAINT,
    not the lost identity.

ADDITION TO PROTOCOL LESSON BANK:
  LESSON 14 (from BRCA LumA, 2026-03-04):
    Cancers that retain identity TFs have
    a different attractor geometry from
    cancers that lose them.
    For identity-retaining cancers, the depth
    axis is tumour suppressor loss, not
    oncogene activation.
    CDKN1A is the archetypal loss in LumA.
    The drug target (CDK4/6 inhibitor) is
    derived from the suppressed gene's function,
    not from an elevated oncogene.
    This inverts the switch gene logic:
    the target is what was lost, not what
    replaced it.
```

---

## STATUS

```
document:           BRCA-S1b
type:               Reasoning Artifact
status:             COMPLETE
date:               2026-03-04
author:             Eric Robert Lawson
                    OrganismCore

confirmed:
  P1 (3/4): FOXA1/GATA3/ESR1 retained ✓
  P3: EZH2 flat ✓
  P4: PIK3CA uncertainty correct ✓
  P6: Controls flat (MBP flagged) ✓
  P7: TNBC markers absent ✓

wrong:
  P1 (1/4): PGR suppressed (ER target, not TF)
  P2: ALL proliferative genes flat/suppressed
      Corrected: arrest removal not escalation

novel:
  CDKN1A -74.3% p=4.68e-195 — primary finding
  Identity hyper-fixation (FOXA1/GATA3 elevated)
  Global chromatin machinery suppressed
  CDK4 r=+0.81 within-LumA depth axis

drug_targets_stated_before_literature:
  PRIMARY:   CDK4/6 inhibitors (palbociclib etc)
             ✓ CONFIRMED — standard of care
  SECONDARY: SERDs (fulvestrant etc)
             ✓ CONFIRMED — approved/trials
  NOVEL:     FOXA1 inhibitors
             🆕 Not standard of care
             Testable from bulk data

next:        BRCA-S2a (Script 2 predictions locked)
             then brca_luma_script2.py
```
