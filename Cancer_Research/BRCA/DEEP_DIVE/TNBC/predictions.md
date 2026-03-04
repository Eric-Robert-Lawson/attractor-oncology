# TNBC — SCRIPT 1 BEFORE-DOCUMENT
## Predictions Locked Before Data Loads
## OrganismCore — Document BRCA-S4a
## Date: 2026-03-04

---

## DOCUMENT METADATA

```
document_id:        BRCA-S4a
series:             BRCA Deep Dive — TNBC (Basal-like)
folder:             Cancer_Research/BRCA/DEEP_DIVE/TNBC/
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
                    (same scRNA-seq dataset as LumA —
                    basal-like cells extracted separately)
                    + GSE25066 — Hatzis et al. 2011
                    (508 pre-treatment bulk samples,
                    pCR annotation — clinical endpoint)
status:             LOCKED — predictions cannot change
                    after this document is committed
precursor_documents:
  ATTRACTOR_GEOMETRY_AXIOMS.md (Doc 90)
  BRCA_Subtype_Orientation.md
  BRCA-S2b (LumA Script 2 artifact)
  BRCA-S2c (LumA literature check)
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
  specifically the BRCA1-regulated luminal progenitor
  in the breast ductal-lobular hierarchy.
  The correct terminal destination is the mature
  luminal epithelial cell.

QUESTION 2: Are identity TFs of the correct
  terminal state expressed?
  NO — TNBC is ER-negative, PR-negative.
  ESR1 and PGR are absent.
  FOXA1 expression is markedly reduced.
  GATA3 is low or absent.
  The luminal identity programme is NOT expressed.

QUESTION 4: Is the cell expressing identity TFs
  of a DIFFERENT normal cell type?
  YES — KRT5, KRT14, SOX10, FOXC1, EGFR.
  These are markers of:
    — Myoepithelial / basal cells (KRT5/14)
    — Neural crest (SOX10)
    — Basal progenitors (FOXC1)
  None of these are the correct terminal identity
  for a luminal progenitor.

PRIMARY TYPE ASSIGNMENT: TYPE 2 — WRONG VALLEY

BUT — COMPOSITE TYPE EVALUATION REQUIRED:

  The founding event in a large fraction of TNBC
  is BRCA1 loss in a LUMINAL PROGENITOR.
  BRCA1 is required for luminal progenitor identity.
  Without BRCA1:
    Stage 1: Luminal progenitor cannot complete
             differentiation (TYPE 1 — BLOCKED APPROACH)
    Stage 2: Progenitor falls into basal/stem-like
             false attractor (TYPE 2 — WRONG VALLEY)

  PREDICTION: TNBC is a COMPOSITE TYPE 1 → TYPE 2.

  This is the first test of the composite type
  axiom stated in Doc 90 Prediction C.

  If confirmed: the correct therapeutic strategy
  addresses both stages:
    Type 1: PARP inhibitors (exploit BRCA1 defect)
    Type 2: EZH2 inhibitors (dissolve basal attractor)
  Their combination is in active trials.
  The composite type predicts why.
```

---

## PART II — BIOLOGICAL GROUNDING
## The lineage and its normal architecture

```
THE NORMAL JOURNEY THIS CANCER INTERRUPTS:

  Mammary stem cell
    ↓
  Luminal progenitor   ← BRCA1 required here
    ↓                    for identity maintenance
  Committed luminal    ← Where TNBC blocks (Type 1)
    ↓
  Mature luminal       ← Correct terminal state
  epithelial cell        (ESR1+, FOXA1+, GATA3+)

THE FALSE ATTRACTOR STATE (Type 2 destination):

  When BRCA1 fails in the luminal progenitor,
  the cell cannot take the luminal path.
  It activates the basal/myoepithelial programme
  — the nearest accessible stable state in the
  Waddington landscape.

  The basal false attractor is maintained by:
    EZH2 — epigenetic silencer of luminal programme
            (PRC2 complex, H3K27me3 at luminal TF loci)
    SOX10 — neural crest / basal identity TF
    FOXC1 — basal progenitor identity TF
    KRT5/14 — structural markers of basal identity
    EGFR — growth factor receptor of basal state

  These are the FA markers — the identity genes of
  the wrong valley the cell has fallen into.

THE SADDLE POINT:

  The saddle point between the normal luminal
  valley and the basal false attractor is maintained
  by BRCA1-dependent chromatin organisation.

  In normal luminal progenitors:
    BRCA1 maintains chromatin accessibility at
    luminal TF loci (ESR1, FOXA1, GATA3).
    When BRCA1 is lost, these loci become inaccessible.
    EZH2 deposits H3K27me3 at these loci.
    The luminal programme is epigenetically silenced.
    The cell is locked in the basal false attractor.

  The saddle point is therefore:
    BRCA1 loss → EZH2 gain → luminal TF silencing
    → basal programme activation → stable false
    attractor with EZH2 as the convergence node.
```

---

## PART III — PREDICTIONS
## All stated 2026-03-04 before data loads

---

### P1 — SWITCH GENES (suppressed in TNBC vs normal)

```
Switch genes are the luminal identity TFs that
were silenced when the luminal progenitor fell
into the basal false attractor.
These are suppressed in TNBC relative to:
  — Normal luminal progenitors
  — Mature luminal cells (the correct terminal state)

PREDICTED SUPPRESSED (switch genes):

  ESR1      Estrogen receptor alpha
            Master luminal TF — drives luminal identity
            Expected: near-absent in TNBC
            CONFIRMED by clinical definition (ER-)
            Direction: DOWN  p predicted: <0.001

  FOXA1     Pioneer TF — opens chromatin for ER binding
            Required for luminal programme initiation
            Expected: suppressed — BRCA1 loss prevents
            chromatin accessibility at FOXA1 targets
            Direction: DOWN  p predicted: <0.001

  GATA3     Luminal differentiation TF
            Co-expressed with ER/FOXA1 in normal luminal
            Expected: suppressed — downstream of FOXA1
            in the luminal programme
            Direction: DOWN  p predicted: <0.001

  SPDEF     Luminal/secretory TF
            Downstream target of FOXA1
            Expected: suppressed
            Direction: DOWN  p predicted: <0.01

  PGR       Progesterone receptor
            ER target gene — confirmed absent (PR-)
            Direction: DOWN  p predicted: <0.001

NOTE: These suppressions are clinically confirmed
by the triple-negative definition.
The question is not IF they are suppressed but
HOW MUCH and whether the depth of suppression
correlates with clinical outcome (pCR).
```

---

### P2 — FALSE ATTRACTOR MARKERS (elevated in TNBC)

```
FA markers are the identity genes of the basal
false attractor — the wrong valley the cell
has fallen into.

PREDICTED ELEVATED (FA markers):

  KRT5      Basal cytokeratin 5
            Marker of myoepithelial / basal cells
            Expected: highly elevated in TNBC
            Direction: UP  p predicted: <0.001

  KRT14     Basal cytokeratin 14
            Co-marker with KRT5 of basal identity
            Expected: highly elevated
            Direction: UP  p predicted: <0.001

  SOX10     Neural crest / basal identity TF
            Elevated specifically in TNBC/basal-like
            Not expressed in luminal subtypes
            Expected: strongly elevated
            Direction: UP  p predicted: <0.001

  FOXC1     Basal progenitor TF
            Drives aggressive basal programme
            Expected: elevated vs luminal normal
            Direction: UP  p predicted: <0.01

  EGFR      Epidermal growth factor receptor
            Overexpressed in basal-like TNBC
            Not amplified — overexpressed
            Expected: elevated vs normal luminal
            Direction: UP  p predicted: <0.01

  VIM       Vimentin — mesenchymal marker
            Rises with depth in TNBC
            Expected: elevated, correlates with depth
            Direction: UP (graded with depth)
```

---

### P3 — CONVERGENCE NODE

```
PREDICTED CONVERGENCE NODE: EZH2

Reasoning:
  EZH2 is the catalytic subunit of PRC2.
  PRC2 deposits H3K27me3 at luminal TF gene loci.
  EZH2 elevation in TNBC simultaneously:
    — Silences ESR1/FOXA1/GATA3 (luminal TFs)
    — Maintains KRT5/KRT14/SOX10 (basal identity)
    — Suppresses differentiation across the genome
    — Maintains chromatin in the basal false attractor

  This is the single node whose inhibition is
  predicted to dissolve the basal false attractor.

  Direction: EZH2 ELEVATED in TNBC vs normal luminal
  Predicted r with depth score: strongly positive

  CONFIRMED EXTERNALLY:
    Schade et al. Nature 635, 755-763 (2024)
    Harvard / Dana-Farber / Ludwig Center
    AKT and EZH2 inhibitors dissolve TNBC
    false attractor and produce luminal conversion.
    Independent derivation confirms EZH2 as
    convergence node.
    This was documented in Doc 71 and Doc 83.
    The prediction is pre-confirmed by this
    independent work.

  PREDICTION: EZH2 is the top positive depth correlate
  within TNBC. r(EZH2, depth) > 0.30 predicted.
```

---

### P4 — THE COMPOSITE TYPE PREDICTION

```
PREDICTION: TNBC shows evidence of BOTH Type 1 AND
Type 2 geometry embedded in the expression data.

Specifically:

  TYPE 1 SIGNAL:
    BRCA1 expression itself should be reduced
    in the cancer cells vs normal luminal progenitors
    (somatic BRCA1 loss or promoter methylation
    in non-germline cases — ~50-60% of TNBC have
    BRCA1 dysfunction by one mechanism or another)

    The luminal progenitor programme was NOT
    expressed at full capacity even before the
    basal programme was activated.
    This should be visible as a DOUBLE SUPPRESSION:
      Luminal TFs suppressed (Type 1 signal)
      AND basal TFs elevated (Type 2 signal)
    where the luminal suppression EXCEEDS what
    would be seen from EZH2 alone.

  TYPE 2 SIGNAL:
    EZH2 elevated — convergence node active
    Basal identity markers elevated
    Epigenetic lock on luminal loci (H3K27me3 proxy:
    look for HDAC1/2 and EZH2 co-elevation)

  THE CRITICAL TEST:
    r(BRCA1, ESR1) within TNBC cells.
    If > 0 and significant: BRCA1 expression
    correlates with residual luminal programme —
    confirms Type 1 component is present.
    If near zero: BRCA1 and luminal identity are
    uncoupled in these cells — Type 1 component
    already completed and only Type 2 remains.

    Prediction: r(BRCA1, ESR1) in TNBC > 0.15
    (partial Type 1 signal still detectable)
```

---

### P5 — DEPTH SCORE PREDICTION

```
PREDICTED DEPTH SCORE CONSTRUCTION:

  Component 1: Basal identity elevation
    norm(mean of KRT5, KRT14, SOX10, FOXC1)
    Higher = more committed to basal false attractor

  Component 2: Luminal identity suppression
    1 - norm(mean of ESR1, FOXA1, GATA3)
    Higher = further from correct luminal identity

  Depth = (Component 1 + Component 2) / 2

PREDICTED DEPTH CORRELATES (top within TNBC):
  EZH2:   r > +0.30  (convergence node rises with depth)
  VIM:    r > +0.25  (mesenchymal character with depth)
  MKI67:  r > +0.20  (proliferation with depth)
  BRCA1:  r < -0.15  (residual BRCA1 falls with depth —
                      deeper = more BRCA1 dysfunction)
  ESR1:   r < -0.30  (luminal TF falls with depth)
  FOXA1:  r < -0.25  (pioneer TF falls with depth)

ANTI-PREDICTED (should NOT correlate with depth):
  CDX2, SPI1, NKX2-1 — non-breast lineage controls
  These should be flat (near-absent in all TNBC)

NOTE ON TYPE 3 CONTRAST:
  In LumA: luminal TFs ELEVATED, correlated with depth
  In TNBC: luminal TFs SUPPRESSED, anti-correlated
  This is the Type 2 vs Type 3 inversion.
  The direction of the identity TF depth correlation
  is the cleanest between-type diagnostic.
```

---

### P6 — pCR PREDICTION (bulk GSE25066)

```
PREDICTED RELATIONSHIP BETWEEN DEPTH AND pCR:

  Hypothesis:
    Deeper TNBC tumors (higher basal false attractor
    commitment, lower residual luminal identity)
    will have LOWER pCR rates to standard
    neoadjuvant chemotherapy.

  Reasoning:
    Deeper basal attractor = more EZH2-mediated
    epigenetic lock = more stable false attractor
    = more resistant to perturbation by cytotoxics.
    Shallower attractor = cells closer to a
    transition state = more vulnerable to
    chemotherapy-induced apoptosis.
    This is the Waddington depth → therapy
    vulnerability prediction.

  ALTERNATIVE HYPOTHESIS (to be tested):
    Deeper attractor = more proliferative (higher MKI67)
    = more sensitive to chemotherapy (anti-metabolites
    and taxanes kill cycling cells preferentially).
    Some evidence in the literature that high TILs
    (immune infiltration) and high Ki-67 predict
    pCR. Both are consistent with the basal false
    attractor being the context.

  The two hypotheses generate OPPOSITE predictions.
  The data will discriminate.

  FORMAL PREDICTION: depth score is negatively
  correlated with pCR rate in TNBC (deeper = less
  likely to achieve pCR).
  r(depth, pCR_binary) < 0 predicted.
  p < 0.05 predicted in GSE25066 TNBC subset.
```

---

### P7 — EPIGENETIC PREDICTION

```
Based on Lesson 5 (Workflow_Protocol v2.0):
  EZH2 direction must be determined from data.
  In BRCA (bulk analysis): EZH2 elevated (gain of
  function lock) — confirmed.
  In LumA: EZH2 near-neutral (+19% vs progenitor, ns)
  In TNBC: EZH2 is the convergence node.

PREDICTION: EZH2 ELEVATED in TNBC vs normal luminal

  EZH2 direction: UP
  Magnitude predicted: >50% vs normal luminal
  p predicted: <0.001
  r with depth score: > +0.30

  Supporting reasoning:
    The Schade 2024 paper (Nature) demonstrates
    EZH2 inhibition produces luminal conversion
    in TNBC. This confirms EZH2 is active and
    maintaining the epigenetic lock.
    The framework independently derived EZH2 as
    the convergence node from the geometry.
    Both converge on EZH2 ELEVATED as the lock.

ADDITIONAL EPIGENETIC:
  HDAC1/2:   Predicted ELEVATED (co-epigenetic
             lock with EZH2 — maintains basal
             chromatin state)
  KDM1A:     Predicted ELEVATED (LSD1 demethylates
             H3K4me2 at luminal TF loci —
             additional silencing layer)
  TET2:      Predicted SUPPRESSED or neutral
             (DNA demethylase — its loss would
             increase methylation at luminal TF
             promoters, supporting silencing)
```

---

### P8 — DRUG TARGET PREDICTIONS
## Before data, before literature (beyond what is known)

```
All stated 2026-03-04 before Script 1 runs.

DRUG TARGET 1 — EZH2 INHIBITORS
  Tazemetostat (EZH2 inhibitor)
  Mechanism: dissolve epigenetic lock on luminal
  TF loci → allow partial luminal reconversion
  → push cells toward the luminal/senescent state
  Schade 2024 (Nature) independently confirmed.
  This is pre-confirmed by external literature.
  Status: ✓ CONVERGENT (already known)

DRUG TARGET 2 — PARP INHIBITORS
  Olaparib, talazoparib
  Mechanism: exploit BRCA1 deficiency (Type 1
  component of the composite type geometry)
  Synthetic lethality: BRCA1-deficient cells
  cannot repair double-strand breaks via
  homologous recombination.
  PARP inhibitors trap PARP at single-strand
  breaks, converting them to DSBs, which kill
  BRCA1-deficient cells selectively.
  Status: ✓ CONVERGENT (approved in BRCA1-mutated
  TNBC — confirms Type 1 component logic)

DRUG TARGET 3 — EZH2 + PARP COMBINATION
  Predicted from composite type axiom (Doc 90):
  Type 1 + Type 2 composite = both type-specific
  drugs should be synergistic.
  EZH2i dissolves the Type 2 false attractor.
  PARPi exploits the Type 1 BRCA1 defect.
  Together: the cell cannot maintain the false
  attractor (EZH2i) AND cannot repair the
  resulting DNA damage (PARPi).
  PREDICTION: EZH2i + PARPi combination is
  synergistic in TNBC beyond what either drug
  achieves alone.
  Status: 🆕 NOVEL derivation from composite
  type geometry — to be checked in literature

DRUG TARGET 4 — DEPTH-STRATIFIED IMMUNOTHERAPY
  Pembrolizumab (anti-PD-L1) is now standard
  with neoadjuvant chemotherapy in TNBC
  (KEYNOTE-522).
  PREDICTION: Depth score stratifies pembrolizumab
  benefit. Shallower tumors (closer to luminal,
  more partial differentiation) may have higher
  TIL infiltration and greater immunotherapy
  response. Deeper tumors (more committed to
  basal false attractor) may be more immune-excluded
  and less pembrolizumab-sensitive.
  Status: 🆕 NOVEL — depth as immunotherapy
  stratification variable not established.

DRUG TARGET 5 — AKT INHIBITORS
  Schade 2024 found AKT + EZH2 inhibition
  synergistic in TNBC (PI3K/AKT pathway
  co-activates EZH2 in basal cells).
  PREDICTION: AKT1 elevated in deeper TNBC cells
  (r(AKT1, depth) > 0.15 predicted).
  Status: ⚠ PARTIALLY PRE-KNOWN (Schade 2024)
  but AKT as depth marker within TNBC is not
  yet established.
```

---

### P9 — INTERNAL TNBC HETEROGENEITY

```
TNBC contains at least 6 molecular subtypes
(Lehmann et al. 2011, 2016).
The GSE176078 scRNA-seq data does not have
Lehmann subtype calls.
The depth score will be applied across ALL
basal-like cells in the dataset.

PREDICTION: Depth score will reveal heterogeneity
within the TNBC population that partially maps
to known Lehmann subtypes.

Specifically:
  Deeper cells (high basal FA markers, low luminal
  switch genes, high EZH2): BL1 or BL2 subtypes
  Intermediate cells (partial luminal character):
  LAR (luminal androgen receptor) subtype

  AR (androgen receptor):
    The LAR subtype is AR-positive.
    PREDICTION: AR expression in TNBC cells is
    negatively correlated with depth.
    Cells that retained partial luminal character
    (shallower in the basal attractor) may retain
    AR as a residual luminal feature.
    r(AR, depth) < -0.15 predicted.
    (This is the opposite of LumA where
    r(AR, depth) = +0.285)

  VIM (vimentin):
    Mesenchymal/claudin-low feature.
    PREDICTION: r(VIM, depth) > +0.20.
    Deeper cells acquire partial EMT character.
```

---

## PART IV — WHAT THIS ANALYSIS CANNOT TELL US

```
The scRNA-seq data (GSE176078) contains:
  Basal-like cells: n=4,312
  Normal references: mature luminal, luminal progenitor

This data was not collected specifically for
TNBC analysis. The "basal SC" annotation is the
closest available population.

LIMITATIONS:
  1. No Lehmann subtype annotation — cannot test
     BL1/BL2/LAR directly
  2. No pCR annotation — the pCR prediction
     (P6) requires GSE25066
  3. No immune cell data from the cancer compartment
     — cannot test TIL depth correlation directly
  4. The normal references are from the same dataset
     — matched donors, not separate cohorts

GSE25066 IS REQUIRED for P6 (pCR prediction).
It is independent of GSE176078.
The script will load both separately.
They will be analyzed in separate steps.
```

---

## PART V — COMPLETE PREDICTION REFERENCE

```
PREDICTIONS LOCKED 2026-03-04
Document: BRCA-S4a
Author: Eric Robert Lawson, OrganismCore

ATTRACTOR TYPE: COMPOSITE TYPE 1 → TYPE 2

P1 — SWITCH GENES (all DOWN in TNBC):
  ESR1, FOXA1, GATA3, SPDEF, PGR
  p < 0.001 for all predicted

P2 — FA MARKERS (all UP in TNBC):
  KRT5, KRT14, SOX10, FOXC1, EGFR, VIM

P3 — CONVERGENCE NODE:
  EZH2 ELEVATED, r(EZH2, depth) > +0.30
  (pre-confirmed by Schade 2024 / Nature)

P4 — COMPOSITE TYPE TEST:
  r(BRCA1, ESR1) within TNBC > +0.15
  (partial Type 1 signal detectable)

P5 — DEPTH SCORE:
  norm(KRT5+KRT14+SOX10+FOXC1) +
  (1 - norm(ESR1+FOXA1+GATA3)) / 2
  Top correlates: EZH2, VIM, MKI67 (positive)
                  ESR1, FOXA1, BRCA1 (negative)

P6 — pCR PREDICTION (GSE25066):
  r(depth, pCR_binary) < 0
  (deeper = lower pCR probability)

P7 — EPIGENETIC:
  EZH2 UP >50% vs normal luminal
  HDAC1/2 elevated
  KDM1A elevated

P8 — DRUG TARGETS:
  ✓ EZH2 inhibitors (tazemetostat)
  ✓ PARP inhibitors (olaparib)
  🆕 EZH2 + PARP combination (from composite type)
  🆕 Depth-stratified pembrolizumab
  ⚠ AKT inhibitors (partially pre-known)

P9 — HETEROGENEITY:
  r(AR, depth) < -0.15 (LAR cells shallower)
  r(VIM, depth) > +0.20 (EMT cells deeper)

CONTROLS (should be flat in TNBC):
  CDX2, SPI1, NKX2-1, NKX3-1, OLIG2
  If any are elevated: analyst assumption error
  Record and process per Wrong Prediction Protocol
```

---

## STATUS BLOCK

```
document:           BRCA-S4a
type:               Before-Document (locked predictions)
date:               2026-03-04
author:             Eric Robert Lawson / OrganismCore
status:             LOCKED

predictions_count:  9 prediction groups
attractor_type:     Composite Type 1 → Type 2
datasets:           GSE176078 (scRNA-seq)
                    GSE25066 (bulk, pCR annotation)
next_document:      BRCA-S4b (Script 1 reasoning artifact)
script:             BRCA_TNBC_script1.py
```
