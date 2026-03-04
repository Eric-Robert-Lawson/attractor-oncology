# TRIPLE-NEGATIVE BREAST CANCER — BEFORE DOCUMENT
## Predictions Locked Before Script 1 Runs
## OrganismCore — Document BRCA-S4a
## Date: 2026-03-04

---

## DOCUMENT METADATA

```
document_id:        BRCA-S4a
series:             BRCA Deep Dive — TNBC (Basal-like)
folder:             Cancer_Research/BRCA/DEEP_DIVE/TNBC/
type:               BEFORE DOCUMENT
                    Predictions only.
                    No data loaded.
                    No results present.
                    This document is locked
                    before Script 1 runs.
date_locked:        2026-03-04
author:             Eric Robert Lawson
                    OrganismCore
protocol_version:   Workflow_Protocol.md v2.0
dataset:            GSE176078 — Wu et al. 2021
                    Nature Genetics | PMID: 34493872
                    100,064 cells — 26 primary breast tumors
                    SAME DATASET AS LUMINAL A DEEP DIVE
                    NO NEW DOWNLOAD REQUIRED
population_of_interest:
                    Cancer Basal SC   (TNBC cancer cells)
normal_reference:
                    Luminal Progenitors  (primary — cell of
                                         origin per current
                                         consensus)
                    Mature Luminal       (secondary — confirm
                                         TNBC identity is absent)
cross_reference:
                    Cancer LumA SC       (confirm LumA/TNBC
                                         geometry are opposite)
additional_dataset:
                    GSE25066 — Hatzis et al. 2011
                    Lancet | PMID: 21641858
                    508 pre-treatment ER-/HER2- bulk RNA-seq
                    pCR annotated — the clinical endpoint test
                    To be incorporated in Script 2 if
                    depth score is established in Script 1
precursor_documents:
                    BRCA_Subtype_Orientation.md
                    ATTRACTOR_GEOMETRY_AXIOMS.md (Doc 90)
                    BRCA-S1a/S1b/S2a/S2b/S2c (LumA series)
next_document:      BRCA-S4b
                    Script 1 results + reasoning artifact
status:             LOCKED — predictions immutable
                    after this date
```

---

## CRITICAL CONTEXT — WHAT IS ALREADY KNOWN

```
Two sources of prior knowledge constrain these predictions.
Neither is allowed to determine them.
Both are stated explicitly so that any influence is visible.

SOURCE 1 — THE PRIOR BULK BRCA ANALYSIS (Documents 75, 82, 83):

  From Document 82 extended (BRCA drug target exploration):
    EZH2 elevated in TNBC/basal SC: +269.7%
    EZH2 is the identified convergence node
    Tazemetostat predicted as drug target
    Confirmed: Schade et al. Nature 2024

  From Document 75 secondary comparison:
    Cancer Basal SC vs Mature Luminal:
      KRT5/14, SOX10, FOXC1 elevated in Basal SC
      FOXA1, GATA3, ESR1 suppressed in Basal SC
    This established the wrong-identity signature.

  What this means for these predictions:
    The prior analysis was bulk-level at the
    population mean. It confirmed the geometry exists.
    The deep dive asks finer questions:
      What is the depth axis WITHIN TNBC?
      What predicts pCR response?
      What is the composite type structure?
    These questions were not asked before.

SOURCE 2 — THE ATTRACTOR GEOMETRY AXIOMS (Document 90):

  Prediction A stated in Document 90:
    TNBC may be a Type 1 → Type 2 composite.
    Type 1: BRCA1 loss blocks luminal progenitor
            from completing differentiation
    Type 2: cell falls into basal false attractor
    Therapeutic implication: combination of
    Type 1 drug (PARP inhibitor) + Type 2 drug
    (EZH2 inhibitor) should show synergy.

  This prediction was stated before this document.
  It is the primary structural hypothesis
  being tested by Script 1.

THE RULE:
  The prior analysis confirmed EZH2 as the
  convergence node. That finding is used as the
  reference point, not as the answer.
  The deep dive must go further:
    How deep is the basin?
    What determines depth within TNBC?
    Where exactly is the saddle point?
    Is the Type 1 element detectable
    in the expression data?
```

---

## SECTION I — BIOLOGICAL GROUNDING

---

### I.1 The normal cell of origin

```
CURRENT CONSENSUS (2024):
  Most TNBC / basal-like cancers arise from a
  BRCA1-regulated LUMINAL PROGENITOR.

  Not the mature luminal cell (LumA's origin).
  Not the myoepithelial cell (despite basal markers).
  The LUMINAL PROGENITOR — the cell partway
  through the luminal differentiation journey
  that BRCA1 normally supervises.

  Why this matters for the geometry:
    A luminal progenitor is already on the path
    toward luminal identity.
    It has partial luminal TF expression.
    Without BRCA1, it cannot complete the journey.
    This is the TYPE 1 component.
    The block at the luminal progenitor → mature
    luminal saddle point.

    But instead of stalling in the progenitor state
    (which would produce a Type 1 cancer like AML),
    the BRCA1-null luminal progenitor activates the
    basal transcriptional programme.
    It falls into the basal false attractor.
    This is the TYPE 2 component.

THE COMPOSITE SEQUENCE:
  Luminal progenitor
    → BRCA1 loss (Type 1 block installed)
    → Cannot complete luminal differentiation
    → Basal programme activated (Type 2 entry)
    → Basal false attractor stabilized by EZH2

  Each stage is detectable in the expression data.
  Script 1 tests whether both stages are visible.
```

---

### I.2 The normal differentiation pathway

```
MAMMARY STEM CELL
  ↓
LUMINAL PROGENITOR  ← TNBC cell of origin
  |
  | ← BRCA1 required for this transition
  ↓
MATURE LUMINAL CELL  ← Normal endpoint (LumA territory)

LUMINAL PROGENITOR
  |
  | ← BRCA1 loss blocks this path
  | ← Basal programme activates instead
  ↓
BASAL FALSE ATTRACTOR  ← TNBC destination

NORMAL MYOEPITHELIAL CELL  ← TNBC resembles this
                              but does NOT originate
                              from it (per consensus)

SADDLE POINT 1 (Type 1 element):
  Luminal progenitor → Mature luminal transition
  Gate: BRCA1 / luminal commitment TFs
  When BRCA1 is lost: gate closes, cell cannot pass
  Switch genes at this saddle:
    BRCA1 itself
    Luminal commitment TFs that BRCA1 enables:
    ESR1, FOXA1, GATA3 — the LumA identity TFs
    that are ACTIVATED as the cell commits to luminal
    These should be SUPPRESSED in TNBC relative to
    Luminal Progenitor (not just vs mature luminal)

SADDLE POINT 2 (Type 2 element):
  The energy barrier between the luminal progenitor
  state and the basal false attractor.
  Once crossed, EZH2 stabilizes the new state.
  The convergence node: EZH2.
```

---

### I.3 What TNBC should look like in the data

```
BASED ON TYPE 1 → TYPE 2 COMPOSITE PREDICTION:

  VERSUS LUMINAL PROGENITOR (Type 1 signal):
    ESR1    SUPPRESSED  (luminal commitment failed)
    FOXA1   SUPPRESSED  (pioneer factor absent)
    GATA3   SUPPRESSED  (luminal TF absent)
    KRT8    SUPPRESSED  (luminal structural gene)
    These being suppressed vs progenitor = the cell
    could not complete luminal commitment = TYPE 1.

  VERSUS MATURE LUMINAL (Type 2 signal — already known):
    KRT5/14  ELEVATED   (basal programme active)
    SOX10    ELEVATED   (basal/neural crest identity)
    FOXC1    ELEVATED   (basal-like TF)
    EGFR     ELEVATED   (basal growth signalling)
    VIM      ELEVATED   (mesenchymal marker)
    EZH2     ELEVATED   (convergence node)

  THE COMPOSITE SIGNAL:
    ESR1/FOXA1/GATA3 suppressed vs PROGENITOR
    = Type 1 element (blocked before luminal commitment)
    KRT5/SOX10/EZH2 elevated vs BOTH references
    = Type 2 element (stabilized in basal state)

  If BOTH are detected, the composite type is confirmed.
  If ONLY Type 2 is detected:
    Either the Type 1 element is too subtle
    or the consensus cell-of-origin is wrong
    (myoepithelial, not luminal progenitor).

THE DECISIVE TEST:
  r(ESR1, EZH2) within TNBC cancer cells.
  If negative: more EZH2 → less ESR1
    EZH2 is actively suppressing luminal programme
    Supports Type 2 mechanism (EZH2 as convergence node
    maintaining the false attractor by silencing
    the luminal identity TFs)
  If near zero: EZH2 and ESR1 are independent
    EZH2 may not be directly suppressing ESR1
    Alternative mechanism for luminal suppression
```

---

## SECTION II — PREDICTIONS
## Locked 2026-03-04 before any data loads

---

### II.1 Switch gene predictions (Type 1 element)

```
PREDICTED SUPPRESSED in TNBC vs Luminal Progenitor:

  P1-SW-1: ESR1   SUPPRESSED  vs Luminal Progenitor
    Role: master luminal TF, BRCA1-regulated
    Reasoning: BRCA1 loss prevents luminal commitment
               The first step of commitment is ESR1 activation
               ESR1 should be LOWER in TNBC even than
               the progenitor it came from
    Direction: DOWN vs progenitor
    Predicted magnitude: >50% suppression

  P1-SW-2: FOXA1  SUPPRESSED  vs Luminal Progenitor
    Role: pioneer factor, opens chromatin for ESR1 binding
    Reasoning: Without BRCA1, FOXA1 cannot establish
               the luminal chromatin state
    Direction: DOWN vs progenitor
    Predicted magnitude: >50% suppression

  P1-SW-3: GATA3  SUPPRESSED  vs Luminal Progenitor
    Role: luminal differentiation TF
    Reasoning: same as FOXA1 — luminal commitment TF
               that BRCA1 enables
    Direction: DOWN vs progenitor
    Predicted magnitude: >30% suppression
    NOTE: GATA3 has lower baseline in progenitors than
          FOXA1 — smaller absolute change expected

  P1-SW-4: KRT8   SUPPRESSED  vs Luminal Progenitor
    Role: luminal structural cytokeratin
    Reasoning: Luminal structural identity lost in TNBC
    Direction: DOWN vs progenitor

  EXPECTED: All four show a gradient:
    Mature Luminal > Luminal Progenitor > TNBC
    This gradient confirms the Type 1 element:
    the cell is FURTHER from luminal identity
    than even the progenitor it came from.
```

---

### II.2 False attractor predictions (Type 2 element)

```
PREDICTED ELEVATED in TNBC vs BOTH references:

  P2-FA-1: KRT5   ELEVATED  vs both references
    Role: basal cytokeratin 5 — basal programme marker
    Reasoning: confirmed in prior bulk analysis
    Direction: UP vs both
    Predicted magnitude: >500% vs mature luminal
                         >200% vs luminal progenitor

  P2-FA-2: KRT14  ELEVATED  vs both references
    Role: basal cytokeratin 14
    Direction: UP vs both

  P2-FA-3: SOX10  ELEVATED  vs both references
    Role: neural crest / basal identity TF
    Reasoning: confirmed in prior bulk analysis
    Direction: UP vs both

  P2-FA-4: EZH2   ELEVATED  vs both references
    Role: PRC2 catalytic subunit — convergence node
    Reasoning: confirmed +269.7% in prior analysis
               confirmed Schade et al. Nature 2024
    Direction: UP vs both
    This is the most confident prediction in the panel.

  P2-FA-5: FOXC1  ELEVATED  vs both references
    Role: basal-like identity TF
    Reasoning: FOXC1 is a documented TNBC marker
    Direction: UP vs both

  P2-FA-6: VIM    ELEVATED  vs both references
    Role: vimentin — mesenchymal/EMT marker
    Reasoning: basal false attractor has partial EMT
    Direction: UP vs both
    NOTE: VIM elevation with depth would suggest
          EMT is a depth axis, not just a binary marker
```

---

### II.3 Epigenetic prediction

```
EZH2: ELEVATED (predicted above)
  Direction: UP — gain of function lock
  This is the Type 2 convergence node.
  EZH2 is elevated to stabilize the basal false
  attractor by silencing the luminal programme.
  This is the GAIN OF FUNCTION epigenetic lock
  (same class as BRCA, different cancer).

  CONTRAST WITH LumA:
    LumA: EZH2 FLAT (ns)
    TNBC: EZH2 ELEVATED (+270%)
  This contrast is among the strongest cross-subtype
  geometric distinctions in the breast cancer series.
  LumA is Type 3 (arrest dismantlement — no epigenetic lock).
  TNBC is Type 2 (false attractor — EZH2 as the lock).

DNMT3A / TET2:
  Direction: UNCERTAIN
  Not predicting confidently.
  Will read from the data.

HDAC1/2:
  Direction: SUPPRESSED in TNBC predicted
  Reasoning: LumA showed HDAC1/2 suppressed.
  TNBC may show same pattern for different reason
  (basal programme requires chromatin opening,
  not closing via HDACs).
  Weak prediction. Will read from data.
```

---

### II.4 Depth axis prediction

```
The depth axis in TNBC is expected to differ
fundamentally from LumA.

LumA depth axis: CDKN1A (arrest dismantlement)
  Cells with less p21 = deeper in the arrest-removed
  valley = more arrest-gone

TNBC depth axis prediction:
  The depth axis should reflect DEGREE OF BASAL
  PROGRAMME ACTIVATION, not arrest removal.
  EZH2 level should positively correlate with depth.
  KRT5/14/SOX10 should correlate with depth.
  ESR1/FOXA1 should NEGATIVELY correlate with depth.

  P3-DEPTH-1:
    r(EZH2, depth) > +0.40 within TNBC
    Higher EZH2 = deeper basal false attractor
    = more luminal programme suppressed

  P3-DEPTH-2:
    r(ESR1, depth) < -0.25 within TNBC
    More ESR1 = shallower (less committed to basal)
    = cells closer to the luminal progenitor origin

  P3-DEPTH-3:
    EZH2 and ESR1 are ANTI-CORRELATED within TNBC:
    r(EZH2, ESR1) < -0.20
    This is the decisive test of whether EZH2 is
    actively suppressing the luminal programme
    or is independently elevated.

  P3-DEPTH-4:
    r(KRT5, depth) > +0.30 within TNBC
    KRT5 rises with depth — deeper basal commitment
    correlates with more basal structural gene expression

  P3-DEPTH-5:
    pCR prediction (if GSE25066 depth is computable):
    Low depth (more ESR1, less EZH2, less KRT5) =
    more residual luminal-progenitor character =
    more genomic instability sensitivity =
    higher pCR to chemotherapy
    High depth (more EZH2, less ESR1, more KRT5) =
    fully committed basal attractor =
    lower pCR = more resistant to chemotherapy
    This is the most clinically significant prediction.
```

---

### II.5 Composite type test predictions

```
THE DECISIVE COMPOSITE TYPE TESTS:

  P4-COMP-1: TYPE 1 GRADIENT VISIBLE
    ESR1/FOXA1/GATA3 should be lower in TNBC
    than in Luminal Progenitors.
    FULL GRADIENT:
      Mature Luminal >> Luminal Progenitor > TNBC
    If TNBC < Luminal Progenitor for all three:
      Type 1 element CONFIRMED in expression data.
    If TNBC ≈ Luminal Progenitor for all three:
      Type 1 element not visible at expression level.
      BRCA1 loss may be functional not expressional.

  P4-COMP-2: TYPE 2 ELEMENT INDEPENDENT OF TYPE 1
    EZH2/KRT5/SOX10 elevated vs Luminal Progenitor
    independently of ESR1/FOXA1 suppression.
    Both signals present simultaneously.

  P4-COMP-3: r(ESR1, EZH2) < -0.15 within TNBC
    The anti-correlation confirms EZH2 suppresses
    ESR1 within TNBC — the convergence node is
    actively maintaining the false attractor by
    silencing the luminal programme in real time.
    This is the single most important correlation
    in the entire TNBC analysis.

  P4-COMP-4: MKI67 DEPTH CORRELATION
    r(MKI67, depth) positive within TNBC
    Deeper basal attractor = more proliferative
    This contrasts with LumA where MKI67 was flat
    (LumA is not a proliferation problem at the
    population mean — TNBC is)
```

---

### II.6 Drug target predictions

```
All stated before data. Stated 2026-03-04.

DRUG 1 — EZH2 INHIBITOR (tazemetostat)
  Mechanism: EZH2 is the convergence node.
  Inhibiting EZH2 dissolves the false attractor
  by removing the epigenetic silencing of the
  luminal programme.
  Status: CONFIRMED (Schade et al. Nature 2024)
          (Framework derived independently)
  Depth prediction: Higher EZH2 = greater benefit
                    from EZH2 inhibition
  Novel extension: EZH2 level at diagnosis should
                   predict magnitude of tazemetostat
                   response — analogous to CDKN1A
                   in LumA predicting CDK4/6i benefit

DRUG 2 — PARP INHIBITOR (olaparib, talazoparib)
  Mechanism: BRCA1 dysfunction is the Type 1 block.
  PARP inhibitors exploit the BRCA1-null state.
  The cell cannot repair DNA double-strand breaks.
  PARP inhibition is synthetic lethal with BRCA1
  deficiency.
  Status: CONFIRMED (approved for BRCA1-mutated TNBC)
  This is the Type 1 drug — it exploits the block,
  not the false attractor.

DRUG 3 — EZH2i + PARPi COMBINATION
  Mechanism: Addresses both stages of composite type.
  Type 1 stage: PARPi exploits the BRCA1 block
  Type 2 stage: EZH2i dissolves the false attractor
  Framework prediction: synergy between the two
  drugs because they address different stages.
  Status: 🆕 NOVEL COMBINATION PREDICTION
  The composite type taxonomy predicts this synergy.
  Single-drug logic from either type alone would
  not predict it.
  Literature check (S4c) will assess whether
  this combination is already in trials.

DRUG 4 — ANTI-PD-L1 / PEMBROLIZUMAB
  Mechanism: TNBC has high TIL infiltration.
  The deepest basal attractors may have the highest
  immune infiltration (EMT and basal programme
  create an inflammatory microenvironment).
  Depth score should correlate with TIL markers
  and predicted pembrolizumab response.
  Status: Pembrolizumab approved for PD-L1+ TNBC
          (KEYNOTE-522)
  Depth prediction: Depth score within TNBC may
                    stratify pembrolizumab benefit
                    better than PD-L1 IHC alone.
                    🆕 NOVEL STRATIFICATION PREDICTION

DRUG 5 — SACITUZUMAB GOVITECAN (TROP2 ADC)
  Mechanism: TROP2 is a surface antigen on TNBC.
  Depth prediction: Does TROP2 expression vary with
                    attractor depth within TNBC?
                    High depth (deep basal) = more TROP2?
                    This is an open question to read from
                    the depth correlation data.
  Status: Approved second-line metastatic TNBC
  Novel question: TROP2 as depth-correlated target
```

---

### II.7 Control gene predictions

```
PREDICTED FLAT in TNBC vs Luminal references:

  SPI1 (PU.1)    — myeloid TF, not breast
  CDX2           — intestinal TF, not breast
  NKX2-1         — lung TF, not breast
  OLIG2          — oligodendrocyte TF, not breast

  LumA identity TFs should be ABSENT in TNBC:
    FOXA1  ABSENT (confirmed prior analysis)
    GATA3  ABSENT (confirmed prior analysis)
    ESR1   ABSENT (confirmed prior analysis)
    PGR    ABSENT (confirmed — TNBC by definition)

  These absent LumA markers in TNBC confirm the
  Type 2 geometry — wrong valley, correct identity
  programme absent.
```

---

## SECTION III — DATASET STRUCTURE CHECK

```
Before Script 1 runs, confirm:

REQUIRED FROM GSE176078 (Wu et al. 2021):
  Population: "Cancer Basal SC"
    — This is the TNBC population
    — Already labeled in metadata.csv
    — No re-classification needed

  Normal reference 1: "Luminal Progenitors"
    — Primary reference for Type 1 test
    — Already labeled in metadata.csv

  Normal reference 2: "Mature Luminal"
    — Secondary reference for cross-comparison
    — Already labeled in metadata.csv

  Cross-reference: "Cancer LumA SC"
    — For LumA vs TNBC geometry confirmation
    — Already labeled in metadata.csv

FILE SOURCES (same as LumA scripts — reuse cache):
  count_matrix_sparse.mtx
  count_matrix_genes.tsv
  count_matrix_barcodes.tsv
  metadata.csv

  IF ALREADY CACHED FROM LUMINAL A SCRIPTS:
    No re-download required.
    Script 1 checks for cache and skips download.

LIBRARY SIZE CHECK (Protocol v2.0 requirement):
  Report library sizes per cell in each population.
  Flag any cell with library size > 5x median
  or < median/5 as outlier.
  TNBC cells have higher proliferation — expect
  slightly higher library sizes than LumA.
  Do not flag proliferation-related library size
  differences as technical artefacts.
  Flag genuine extreme outliers only.

GENE NAME FORMAT:
  Human (all caps): BRCA1, EZH2, KRT5
  Confirmed from LumA scripts.
  No mouse gene name issues expected.
```

---

## SECTION IV — SCRIPT 1 OUTPUT STRUCTURE

```
Per Protocol v2.0 — discovery frame enforced.
Output must be in this order:

SECTION 1: TOP MOVERS (unfiltered)
  Top 20 genes GAINED in TNBC vs Luminal Progenitor
  Top 20 genes LOST in TNBC vs Luminal Progenitor
  Read cold. No prediction panel imposed.

SECTION 2: TOP MOVERS vs MATURE LUMINAL
  Top 20 GAINED
  Top 20 LOST
  For comparison with Luminal Progenitor results.

SECTION 3: TYPE 1 GRADIENT TEST
  ESR1, FOXA1, GATA3, KRT8 across three populations:
  Mature Luminal > Luminal Progenitor > TNBC?
  Is the gradient confirmed?

SECTION 4: PCA GEOMETRY
  PC1 variance explained
  PC1 top loadings
  Does PC1 separate TNBC from normal references?
  What do the loadings encode — basal or luminal axis?

SECTION 5: DEPTH SCORE
  Blind depth score from top movers
  EZH2-centred component
  ESR1/FOXA1 inverse component
  Depth correlations within TNBC

SECTION 6: PREDICTION PANEL CHECK
  Where do all predictions above land in the landscape?
  Confirmed / not confirmed / inverted / unexpected

SECTION 7: COMPOSITE TYPE ASSESSMENT
  r(ESR1, EZH2) within TNBC — the decisive test
  Type 1 gradient present/absent
  Type 2 markers present/absent
  Composite type confirmed/rejected/partial

SECTION 8: TNBC vs LUMA CROSS-COMPARISON
  Top movers in TNBC that are opposite in LumA
  Geometry contrast: Type 2 vs Type 3

A script that leads with the prediction panel
has failed to implement the discovery frame.
```

---

## SECTION V — THE COMPLETE PREDICTION SUMMARY

```
ALL PREDICTIONS LOCKED 2026-03-04
BEFORE ANY TNBC DATA IS LOADED

TYPE 1 ELEMENT (BRCA1 block):
  P1-SW-1: ESR1   DOWN vs Luminal Progenitor  >50%
  P1-SW-2: FOXA1  DOWN vs Luminal Progenitor  >50%
  P1-SW-3: GATA3  DOWN vs Luminal Progenitor  >30%
  P1-SW-4: KRT8   DOWN vs Luminal Progenitor

TYPE 2 ELEMENT (Basal false attractor):
  P2-FA-1: KRT5   UP   vs both references     >500% vs mature
  P2-FA-2: KRT14  UP   vs both references
  P2-FA-3: SOX10  UP   vs both references
  P2-FA-4: EZH2   UP   vs both references     confirmed prior
  P2-FA-5: FOXC1  UP   vs both references
  P2-FA-6: VIM    UP   vs both references

EPIGENETIC:
  P-EPI-1: EZH2 elevated, gain-of-function lock

DEPTH AXIS:
  P3-DEPTH-1: r(EZH2,  depth) > +0.40 within TNBC
  P3-DEPTH-2: r(ESR1,  depth) < -0.25 within TNBC
  P3-DEPTH-3: r(EZH2,  ESR1)  < -0.20 within TNBC
  P3-DEPTH-4: r(KRT5,  depth) > +0.30 within TNBC
  P3-DEPTH-5: pCR prediction (Script 2 with GSE25066)

COMPOSITE TYPE:
  P4-COMP-1: ESR1/FOXA1/GATA3 lower in TNBC
             than in Luminal Progenitors
  P4-COMP-2: EZH2/KRT5/SOX10 elevated vs progenitor
             independently
  P4-COMP-3: r(ESR1, EZH2) < -0.15 within TNBC
  P4-COMP-4: r(MKI67, depth) > 0 within TNBC

DRUG TARGETS:
  DRUG 1: EZH2 inhibitor (tazemetostat)    ✓ prior confirmation
  DRUG 2: PARP inhibitor (olaparib)        ✓ approved
  DRUG 3: EZH2i + PARPi combination       🆕 novel from composite type
  DRUG 4: Pembrolizumab depth stratified  🆕 novel stratification
  DRUG 5: Sacituzumab govitecan (TROP2)   ? depth correlation open

CONTROLS (all predicted FLAT):
  SPI1, CDX2, NKX2-1, OLIG2

LUMA CONTRAST (all predicted ABSENT/SUPPRESSED):
  FOXA1, GATA3, ESR1, PGR, CDKN1A (LumA markers)
```

---

## STATUS

```
document:           BRCA-S4a
date_locked:        2026-03-04
predictions:        26 total
  Type 1 element:   4 predictions
  Type 2 element:   6 predictions
  Epigenetic:       1 prediction
  Depth axis:       5 predictions (1 deferred to S2)
  Composite type:   4 predictions
  Drug targets:     5 predictions (2 novel, 3 confirmed)
  Controls:         4 predictions

author:             Eric Robert Lawson
                    OrganismCore
status:             LOCKED
```
