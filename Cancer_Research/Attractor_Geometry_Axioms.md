# ATTRACTOR GEOMETRY AXIOMS
## The Invariant Structural Types of Cancer as False Attractors
## Empirically Derived from 14 Cancer Validations
## OrganismCore — Document 90
## Date: 2026-03-04

---

## PREAMBLE — HOW THIS DOCUMENT CAME TO EXIST

```
This document was not written before the cancer series.
It could not have been.

The three geometric types described here were not
assumed as a theoretical framework and then tested.
They emerged from the data across 14 independent
cancer analyses and were first recognized as a
unified taxonomy during the Luminal A deep dive
(BRCA-S2b, 2026-03-04).

That distinction matters.

A taxonomy imposed before data is a hypothesis.
A taxonomy that emerges from data across 14
independent cancer types is an empirical finding.

This document records the finding.
It does not claim the finding is complete.
It explicitly invites revision as more cancers
are analyzed.

The axioms stated here are the strongest kind:
  Derived, not assumed.
  Convergent across multiple independent datasets.
  Each type confirmed by drug target derivation
  that matched published pharmacology.
  Each type falsifiable by future analysis.

Read this document before writing any new
cancer before-document.
The three types determine what to look for
before the data is loaded.
```

---

## DOCUMENT METADATA

```
document_id:        ATTRACTOR_GEOMETRY_AXIOMS
document_number:    90
type:               Framework axiom document
                    Cross-cancer structural invariants
date_derived:       2026-03-04
                    (taxonomy first named explicitly
                    during BRCA LumA S2b analysis)
date_written:       2026-03-04
author:             Eric Robert Lawson
                    OrganismCore
protocol_version:   Workflow_Protocol.md v2.0
derived_from:       14 cancer validations:
                    AML, CML, CRC, GBM, BRCA (bulk),
                    LUAD, B-ALL, T-ALL, CLL, MM,
                    MDS, PAAD, PRAD, STAD,
                    + BRCA LumA deep dive
status:             LIVING DOCUMENT
                    Update as new types are identified
                    or existing types are refined
                    Version 1.0 — initial derivation
supersedes:         No prior document
                    This taxonomy was implicit in the
                    series but never formally stated
```

---

## THE CENTRAL CLAIM

```
Every cancer analyzed in this repository is a false
attractor in the Waddington epigenetic landscape.

This is not in dispute across the series.
What this document adds is the observation that
false attractors are not all geometrically identical.

They differ in one critical structural dimension:

  THE RELATIONSHIP BETWEEN THE CANCER CELL'S
  POSITION IN THE LANDSCAPE AND THE CORRECT
  DIFFERENTIATION VALLEY.

That relationship takes three distinct forms.
These forms are not arbitrary categories.
They are geometrically distinct configurations
of the Waddington landscape.
Each configuration has a distinct:
  — Diagnostic signature in expression data
  — Mechanism of attractor stabilization
  — Location of the saddle point
  — Drug logic for dissolution

The three types are stated below as axioms —
meaning: propositions that, if true, predict
the structure of every cancer analysis that
applies this framework.

They are not assumed.
They are the distilled residue of 14 validations.
```

---

## AXIOM I — THE BLOCKED APPROACH
### (Type 1 — Stuck Above the Valley)

---

### I.1 Geometric statement

```
The cancer cell originated from a progenitor that
was in transit toward a normal differentiation
attractor and could not complete the journey.

The cell is not in the correct valley.
The cell is not in a different valley.
The cell is ABOVE the valley — stalled on the slope
before reaching the stable differentiated state.

The valley exists.
The path exists.
The cell cannot take it.

                    Normal
                   attractor
      Cancer          ↓
      cell →    ████████████
      (stuck)   ██ normal  ██
                ████████████
                  (valley)

The false attractor is the stalled-progenitor state.
It is stable because the transition to the normal
attractor is blocked — not because the cancer cell
has acquired a new stable identity.
```

### I.2 Expression signature

```
WHAT IS SUPPRESSED:
  The differentiation TFs required to cross the
  saddle point into the terminal identity valley.
  These are the SWITCH GENES.
  They are the genes whose activation would dissolve
  the attractor by completing the journey.

  Pattern: the switch genes are the master TFs of
  the NEXT STAGE in the differentiation sequence —
  the stage the cancer cell cannot reach.

WHAT IS ELEVATED:
  The identity genes of the PROGENITOR STATE where
  the cell is stuck.
  These are the FALSE ATTRACTOR MARKERS.
  They are the surface markers, TFs, and
  structural genes of the stalled state.

  Pattern: the FA markers are the canonical
  markers of the progenitor population from which
  the cancer derives — expressed at high levels
  because the cell is indefinitely maintained in
  that state.

KEY DIAGNOSTIC:
  The switch genes are SUPPRESSED vs both the
  normal progenitor reference AND the normal
  terminal cell reference.
  The FA markers are ELEVATED vs the terminal
  cell reference.
  Identity TFs of the WRONG (progenitor) state
  are retained.
  Identity TFs of the CORRECT (terminal) state
  are absent.
```

### I.3 Confirmed examples

```
AML (Acute Myeloid Leukemia):
  Normal journey: HSC → GMP → granulocyte
  Block:          blast cannot cross HSC→GMP saddle
  Switch genes:   SPI1 (PU.1), KLF4, IRF8
                  (GMP-stage TFs — suppressed)
  FA markers:     CD34, HOXA9
                  (progenitor identity — retained)
  Drug:           CRISPRa activation of SPI1/KLF4
                  (confirmed: independent CRISPR work)

CML (Chronic Myeloid Leukemia):
  Normal journey: HSC → myeloid lineage commitment
  Block:          BCR-ABL locks HSC in progenitor
  Switch gene:    GATA1/2 (myeloid commitment)
  Drug:           Imatinib (BCR-ABL inhibitor)
                  (confirmed: standard of care)

B-ALL (B-cell Acute Lymphoblastic Leukemia):
  Normal journey: HSC → pro-B → pre-B → mature B
  Block:          cannot cross pro-B saddle
  Switch gene:    PAX5 (B-cell commitment — suppressed)
  FA markers:     CD34, CD10 (progenitor retained)
  Drug:           Confirmed by literature

MDS (Myelodysplastic Syndrome):
  Normal journey: GMP → promyelocyte → myelocyte
  Block:          one stage LATER than AML
                  (different saddle point in same
                  lineage — proves each clinical
                  entity has a distinct geometry)
  Switch gene:    ELANE (neutrophil elastase —
                  effector-level block, not TF)
  Drug:           Confirmed by literature
```

### I.4 Drug logic

```
For Type 1 cancers, the therapeutic objective is:

  COMPLETE THE INTERRUPTED JOURNEY.

The cell is close to the correct valley.
It needs to cross a saddle point it cannot
cross on its own.

The minimal control set is the set of TFs or
effectors that, when activated, push the cell
over the saddle point and into the normal
terminal attractor.

Once in the terminal attractor, normal
differentiation programs execute.
The cell may arrest, senesce, or terminally
differentiate — all of which end its malignant
activity.

Drug classes:
  — TF activation (CRISPRa, small molecule activators)
  — Kinase inhibition (removes the block maintaining
    the progenitor state — BCR-ABL in CML)
  — Differentiation inducers (ATRA in APL — the
    best existing example of Type 1 therapy)
  — Corepressor inhibition (removes epigenetic
    silencing of switch genes)

The diagnostic prediction:
  For Type 1 cancers, the most important biomarker
  is the expression level of the switch gene.
  Lower switch gene = deeper in the false attractor
  = larger minimal control set required.
```

---

## AXIOM II — THE WRONG VALLEY
### (Type 2 — Stuck in an Incorrect Attractor)

---

### II.1 Geometric statement

```
The cancer cell has crossed into and stabilized in
a different attractor basin — one that does not
correspond to any normal terminal differentiation
state of the cell's lineage of origin.

The cell is not stalled on the path to the correct
valley. It has arrived at a different valley entirely.

The normal differentiation attractor exists.
The cell is not near it.
The cell is in a false attractor that has its own
geometry, its own identity programme, and its own
stabilizing mechanisms.

      Normal           FALSE
     attractor       attractor
        ↓               ↓
  ████████████    ████████████
  ██ normal  ██   ██ cancer ██
  ████████████    ████████████
       ↑
  Cell should
  be here but
  is here ──────────────────→

The false attractor has basin depth.
The deeper the basin, the larger the intervention
required to push the cell out.
```

### II.2 Expression signature

```
WHAT IS SUPPRESSED:
  The identity TFs of the CORRECT normal terminal state.
  The lineage-specific differentiation markers.
  These genes are not just absent — they are
  actively maintained in a suppressed state by the
  convergence node (the epigenetic or transcriptional
  lock that stabilizes the false attractor).

WHAT IS ELEVATED:
  The identity genes of the FALSE ATTRACTOR STATE.
  These may correspond to:
    — A different normal cell type (trans-identity)
    — A developmentally primitive state
      (de-differentiation)
    — A hybrid state with no normal counterpart
    — A stress-response or wound-healing programme

KEY DIAGNOSTIC:
  The switch genes are SUPPRESSED specifically
  relative to the cell's correct normal counterpart.
  The FA markers are ELEVATED and correspond to
  a different normal cell type's identity programme
  or a primitive stem-like state.
  There is a CONVERGENCE NODE — a single gene or
  pathway that simultaneously maintains switch gene
  suppression and FA marker elevation.
  The convergence node is the drug target.
```

### II.3 Confirmed examples

```
GBM (Glioblastoma Multiforme):
  Normal cell:    mature neuron / astrocyte
  False attractor: oligoprogenitor-like state
  Switch genes:   neural differentiation TFs
                  (suppressed)
  FA markers:     OLIG2 (oligoprogenitor identity —
                  the wrong valley)
  Convergence:    OLIG2 (drug target)
  Drug:           CT-179 (OLIG2 inhibitor)
                  Phase 1 trial Oct 2025 ✓

TNBC (Triple-Negative Breast Cancer):
  Normal cell:    luminal progenitor (BRCA1-regulated)
  False attractor: basal/stem-like state
  Switch genes:   luminal progenitor TFs (suppressed)
  FA markers:     KRT5/14, SOX10, FOXC1
                  (basal identity — wrong valley)
  Convergence:    EZH2 (confirmed, Schade et al.
                  Nature 2024 ✓ independent)
  Drug:           Tazemetostat (EZH2 inhibitor)

CLL (Chronic Lymphocytic Leukemia):
  Normal cell:    mature B lymphocyte
  False attractor: apoptosis-resistant B cell state
  Convergence:    BCL2 (maintains false attractor
                  against programmed death)
  Drug:           Venetoclax (BCL2 inhibitor)
                  FDA approved ✓

MM (Multiple Myeloma):
  Special case — within-terminal-state false attractor
  Normal cell:    plasma cell (already terminal)
  False attractor: immortal plasma cell
                  (arrested within the terminal state
                  against senescence/death)
  Switch gene:    IRF8 (commitment marker — suppressed)
  Drug:           Proteasome inhibitors / IMiDs
                  (confirmed standard of care)
```

### II.4 Drug logic

```
For Type 2 cancers, the therapeutic objective is:

  DISSOLVE THE FALSE ATTRACTOR.

The cell is not near the correct valley.
Pushing it toward the correct valley requires
crossing a large energy barrier.
The more tractable intervention is to destabilize
the false attractor itself — remove the convergence
node that is maintaining it — and allow the cell
to find ANY more stable state (which may be
quiescence, senescence, apoptosis, or
differentiation depending on what is available).

Drug classes:
  — Convergence node inhibitors (EZH2i, BCL2i,
    OLIG2i — direct targeting of the node that
    holds the false attractor stable)
  — Epigenetic reprogramming (demethylation,
    HAT activators — dissolve the epigenetic lock)
  — Synthetic lethality exploitation (the false
    attractor creates dependencies that the normal
    cell does not have — PARP inhibitors in BRCA1-
    mutated TNBC exploit this)

Basin depth prediction (from Document 71):
  Deeper false attractor = larger minimal control set.
  Velocity coherence ratio predicts MCS size before
  running the full analysis.
  Deep basin (score > 1.5) → 5+ gene MCS
  Shallow basin (score < 1.2) → 1-2 gene MCS

The diagnostic prediction:
  For Type 2 cancers, the convergence node expression
  level is the depth marker.
  EZH2 level in TNBC.
  BCL2 level in CLL.
  OLIG2 level in GBM.
  Higher convergence node = deeper attractor =
  more drug required to dissolve it.
```

---

## AXIOM III — THE OVERSHOT IDENTITY
### (Type 3 — Correct Valley, Floor Removed)

---

### III.1 Geometric statement

```
The cancer cell is in the CORRECT differentiation
valley. Its identity programme is intact — or
more than intact: it has been driven past the
normal terminal state deeper into the identity
programme than normal cells go.

But the ARREST mechanisms that normally form the
floor of the valley — preventing indefinite descent
and anchoring the cell at the correct depth —
have been dismantled.

The cell is in the right valley.
It has gone too far into it.
And it cannot stop.

      Normal          
     attractor       
        ↓               
  ████████████    
  ██ normal  ██   ← Normal cell arrested here
  ████████████    
        ↓               
  ░░░░░░░░░░░░    ← Walls removed (TGF-β/CDKN1A gone)
  ░░ cancer ░░    ← Cancer cell keeps descending
  ░░░░░░░░░░░░    
        ↓
   (no floor)

The false attractor here is unusual: it is
located WITHIN the correct differentiation valley,
deeper than the normal terminal state.
The cell has the correct identity but has lost the
mechanism that says "stop here."
```

### III.2 Expression signature

```
WHAT IS SUPPRESSED:
  NOT the identity TFs of the correct state.
  These are RETAINED or ELEVATED.
  The suppression is in the ARREST AXIS:
    — Tumour suppressor pathway components
    — Cell cycle inhibitors (CDKIs)
    — TGF-β receptor or downstream transducers
    — Senescence pathway genes
  The arrest axis is dismantled at multiple levels
  simultaneously — not just one gene suppressed but
  the full pathway from receptor to effector.

WHAT IS ELEVATED:
  The identity genes of the CORRECT terminal state —
  often elevated BEYOND normal terminal cell levels.
  Compared to progenitors: massive elevation.
  Compared to mature normal: retained or slightly
  reduced but above progenitor by orders of magnitude.
  There is NO alternative identity programme active.
  The false attractor markers ARE the correct
  identity markers — just in excess.

  Compared to progenitor:        FOXA1 +1038%
                                 GATA3  +286%
                                 ESR1   +606%
  Compared to mature normal:     FOXA1  +37%
                                 GATA3  +34%
                                 ESR1   flat

KEY DIAGNOSTIC — THE INVERSION TEST:
  In Type 1 and Type 2 cancers:
    Identity TFs of the CORRECT state are SUPPRESSED
  In Type 3 cancers:
    Identity TFs of the CORRECT state are ELEVATED
    (or retained) relative to BOTH references

  This inversion is the diagnostic signature.
  If you see a cancer where the differentiation
  TFs are elevated rather than suppressed,
  you are in Type 3 geometry.
  Do not interpret this as "the framework is not
  working." Interpret it as a different type of
  false attractor.

ADDITIONAL KEY DIAGNOSTIC — THE LIGAND PARADOX:
  In Type 3, the ligand for the arrest pathway
  is often ELEVATED (the cell is trying to arrest
  itself) while the receptor or transducer is absent.
  The cell produces the arrest signal.
  It cannot receive it.
  This is receptor-level resistance to self-arrest.
  TGFB1 +54% with TGFBR2 -97% in LumA is the
  confirmed example.
  The arrest machinery on the receiving end is gone.
  The signal being sent is not.
```

### III.3 Confirmed examples

```
Luminal A Breast Cancer (BRCA LumA deep dive):
  Normal cell:    mature luminal epithelial cell
  False attractor: hyper-committed luminal state
                  with arrest axis removed
  Identity TFs:   FOXA1/GATA3/ESR1 all ELEVATED
                  (not suppressed — Type 3 signature)
  Arrest axis:    TGFBR2 -97% / SMAD3 -44% /
                  CDKN1A -69% (full axis dismantled)
  Ligand paradox: TGFB1 +54% (ligand elevated,
                  receptor absent)
  Circuit gap:    r(SMAD3, CDKN1A) = +0.041
                  (pathway uncoupled, not just depleted)
  Drug:           CDK4/6 inhibitors (compensate for
                  absent p21) ✓ standard of care
                  TGFBR2 restoration (causal target) 🆕

CANDIDATE — Invasive Lobular Carcinoma (ILC):
  Not yet analyzed in deep dive series.
  Preliminary structural reasoning:
  CDH1 (E-cadherin) loss as the founding event
  in a cell that retains full luminal identity
  (ER+, FOXA1+, GATA3+).
  The cell did not acquire a wrong identity.
  It lost a structural component (adhesion) while
  maintaining its full transcriptional identity.
  This may be Type 3 with the architectural
  arrest rather than the cell cycle arrest
  as the dismantled component.
  To be confirmed in BRCA-S6 series.

CANDIDATE — Low-grade Prostate Cancer (early PRAD):
  Prior PRAD analysis identified NKX3-1 suppression
  as the switch gene.
  But the luminal prostate identity TFs
  (FOXA1, AR) are retained.
  The early-stage geometry may be Type 3 with the
  PTEN-driven arrest axis as the dismantled component.
  AR/FOXA1 retention in early PRAD is consistent
  with identity overshoot rather than identity loss.
  Requires re-examination with the Type 3 lens.
```

### III.4 Drug logic

```
For Type 3 cancers, the therapeutic objective is:

  RESTORE THE WALLS — RE-ESTABLISH THE ARREST AXIS.

The cell is in the right valley.
It does not need to be pushed anywhere.
It needs to be stopped.

Two strategies:
  
  COMPENSATORY (current clinical approach):
    Replace the missing arrest mechanism with a drug
    that performs the same function.
    CDK4/6 inhibitors substitute for absent p21 /
    CDKN1A. The drug does what CDKN1A cannot.
    Limitation: the drug must be administered
    continuously. Remove the drug and cycling resumes.
    The underlying arrest axis remains absent.
    This explains CDK4/6 inhibitor resistance
    development: the cell finds alternative ways to
    cycle that bypass CDK4/6, which it can do because
    it was already adapted to low p21.

  CAUSAL (what the framework predicts is needed):
    Restore the missing component of the arrest axis
    so the cell can re-arrest itself.
    TGFBR2 restoration in LumA: the TGF-β ligand
    is already elevated. Restoring the receptor
    allows the cell's own TGF-β to re-engage
    SMAD3 → CDKN1A → self-arrest.
    This is a durable intervention because the
    cell regains autonomous arrest capacity.
    No continuous drug administration required.
    Limitation: no TGFBR2 restoration pharmacology
    currently exists at clinical stage.

DEPTH PREDICTION for Type 3:
  Depth = degree of arrest axis dismantlement.
  CDKN1A level (bulk CV = 0.976 in LumA).
  Lower CDKN1A = more dismantlement = deeper.
  Deeper = greater CDK4/6 inhibitor benefit
          (compensatory strategy: more reliance
          on the drug because less endogenous arrest)
  Deeper = greater TGFBR2 restoration benefit
          (causal strategy: the arrest axis has
          been more completely removed and must
          be more completely restored)

Drug classes:
  — CDK4/6 inhibitors (compensatory — confirmed)
  — Receptor restoration for arrest pathway
    (causal — novel, not yet at clinical stage)
  — SMAD activators (restore transducer function)
  — Coactivator stabilization (for ER circuit
    in LumA — addresses coactivator depletion)
  — Endocrine therapy + SERD
    (addresses ER programme throttling — confirmed)
```

---

## PART II — THE CROSS-TYPE COMPARISON TABLE

```
Property              Type 1           Type 2           Type 3
                    BLOCKED          WRONG VALLEY     CORRECT VALLEY
                    APPROACH         (False Identity)  FLOOR REMOVED
──────────────────────────────────────────────────────────────────────
Cell position       Above correct    Inside wrong     Inside correct
in landscape        valley           valley           valley, too deep

Identity TFs of     ABSENT           ABSENT           PRESENT or
correct lineage                                       OVER-EXPRESSED

Identity TFs of     PRESENT          ABSENT           ABSENT
progenitor state    (retained)       (different       (no alternative
                                     identity         identity — just
                                     expressed)       over-committed)

FA markers          = Progenitor     = Wrong-identity = Correct identity
are...              identity markers identity markers genes, elevated

What is the         The saddle       The convergence  The arrest axis
convergence node?   point TF that    node maintaining (CDKN1A, TGFBR2,
                    the cell cannot  the false        SMAD3, or
                    activate         identity (EZH2,  equivalent)
                                     BCL2, OLIG2...)

Ligand paradox?     No               No               YES — arrest
                                                      ligand elevated,
                                                      receptor absent

Therapeutic goal    Complete the     Dissolve the     Restore the
                    journey          false attractor  walls / arrest

Drug logic          Activate switch  Inhibit          Compensate for
                    gene / remove    convergence node absent arrest
                    the block        (epigenetic or   (CDK4/6i) OR
                                     anti-apoptotic)  restore arrest
                                                      axis causally

Resistance          Block reforms    New convergence  Alternative
mechanism           (alternative     node replaces    cycling pathway
                    lock forms)      inhibited node   bypasses CDK4/6

Depth metric        Switch gene      Convergence node Arrest axis
                    expression level expression /     expression level
                    (lower = deeper) basin depth      (lower = deeper)

Canonical examples  AML, MDS,        TNBC (EZH2),    LumA (CDKN1A/
                    B-ALL, CML       GBM (OLIG2),    TGFBR2 axis)
                                     CLL (BCL2),
                                     MM

Drug confirmation   CRISPRa (AML)    Venetoclax (CLL) Palbociclib (LumA)
status              Imatinib (CML)   Tazemetostat     All confirmed ✓
                    ATRA (APL)       (TNBC) ✓
                    All confirmed ✓  CT-179 (GBM) ✓
```

---

## PART III — THE COMPOSITE TYPE

```
Some cancers may involve more than one type.
This was first identified in the TNBC structural
reasoning during the LumA post-analysis.

TNBC may be a Type 1 → Type 2 CASCADE:

  The founding event in BRCA1-mutated TNBC:
    BRCA1 loss in a LUMINAL PROGENITOR.
    BRCA1 is required for luminal progenitor
    identity maintenance.
    Without BRCA1, the luminal progenitor cannot
    complete differentiation (Type 1 block).
    But instead of stalling indefinitely in the
    progenitor state, it falls into the
    basal/stem-like false attractor (Type 2).

  The geometry:
    Stage 1: Type 1 — block above luminal valley
    Stage 2: Type 2 — falls into basal false attractor

  The therapeutic implication:
    Type 1 logic alone would predict:
      restore BRCA1 / push toward luminal
    Type 2 logic alone would predict:
      dissolve the basal false attractor via EZH2
    The correct answer may be both:
      EZH2 inhibition (Type 2: dissolve the
      basal attractor) + PARP inhibition (exploits
      the BRCA1 defect that caused the Type 1 block)
    This is EXACTLY the combination being tested
    in clinical trials.
    The composite type predicts combination therapy.
    Single-type logic would not.

ILC MAY BE A TYPE 3 VARIANT:
  CDH1 loss = structural arrest component removed
  (adhesion/architecture rather than cell cycle)
  Luminal identity retained (ER+, FOXA1+, GATA3+)
  Type 3 with architectural dismantlement
  rather than cell cycle dismantlement.
  To be confirmed in BRCA-S6 analysis.

RULE FOR COMPOSITE TYPES:
  If a cancer does not fit cleanly into one type,
  ask: does it involve a SEQUENCE of type transitions?
  Map the sequence before writing the before-document.
  The therapeutic strategy should address each stage
  of the sequence, not just the final attractor state.
```

---

## PART IV — DIAGNOSTIC ALGORITHM

```
When beginning a new cancer analysis, before loading
any data, determine the type using this algorithm:

QUESTION 1:
  What is the normal cell of origin?
  What is the terminal differentiation state
  it should reach?

QUESTION 2:
  Is the cancer cell expressing the identity TFs
  of the correct terminal state?

  YES → Go to QUESTION 3
  NO  → Go to QUESTION 4

QUESTION 3 (correct identity TFs present):
  Are the identity TFs elevated ABOVE normal
  terminal cell levels? (vs progenitor: massively
  elevated? vs mature: retained or slightly elevated?)

  YES → TYPE 3 — Correct valley, floor removed
        Look for: arrest axis genes (CDKIs, TGF-β
        receptors, senescence pathway)
        Look for: ligand paradox (arrest ligand
        elevated + receptor absent)
        Drug logic: restore the walls

  NO  → INDETERMINATE — The cell is in the
        correct valley at the correct depth.
        May not be a cancer.
        May be a Type 3 at very shallow depth.
        Re-examine the arrest axis genes.

QUESTION 4 (correct identity TFs absent):
  Is the cell expressing identity TFs of a
  DIFFERENT normal cell type?

  YES → Likely TYPE 2 — Wrong valley
        Look for: convergence node maintaining
        the wrong identity
        Look for: basin depth markers
        Drug logic: dissolve the false attractor

  NO  → Is the cell expressing the identity TFs
        of an EARLIER progenitor state?

        YES → TYPE 1 — Blocked approach
              Look for: the saddle point TF(s)
              the cell cannot activate
              Drug logic: complete the journey

        NO  → COMPOSITE or UNKNOWN type
              Map the lineage history
              Look for multi-stage sequence
              Apply composite type analysis

SHORTCUT — THE IDENTITY TF DIRECTION TEST:
  In the data, before any interpretation:
  Are the differentiation TFs of the correct
  lineage SUPPRESSED or ELEVATED?

  SUPPRESSED → Type 1 or Type 2
  ELEVATED   → Type 3
               (highest probability rule)
```

---

## PART V — PREDICTIONS GENERATED BY THE AXIOMS

```
These are predictions that the three-type taxonomy
makes about analyses not yet performed.
Stated here before those analyses run.
Dated 2026-03-04.

PREDICTION A — ILC is Type 3 with architectural axis:
  Invasive Lobular Carcinoma will show:
    Luminal identity TFs retained (ER+, FOXA1+, GATA3+)
    CDH1 absent (founding event — architectural arrest
    component removed, not cell-cycle arrest)
    Adjacent non-CDH1 structural genes may also
    show arrest-like suppression
    Drug logic: CDH1 restoration strategies or
    downstream stabilization of p120-catenin
    signalling
  This is a Type 3 variant with adhesion/architecture
  as the dismantled arrest component rather than
  cell cycle.

PREDICTION B — Early PRAD may be Type 3:
  If re-examined with the Type 3 lens:
    AR and FOXA1 are retained in early PRAD
    PTEN loss removes the arrest floor (PI3K/AKT)
    NKX3-1 suppression may be secondary to
    AR over-activation rather than primary switch
    gene suppression
  This would reclassify early PRAD geometry.
  Implications: AR inhibition (compensatory) vs
  PTEN/PI3K pathway restoration (causal) as
  the correct drug logic distinction.

PREDICTION C — Composite type predicts combination therapy:
  Any cancer where the composite type analysis
  reveals a Type 1 → Type 2 sequence should
  show synergy between:
    — The Type 1 drug (remove the block)
    — The Type 2 drug (dissolve the false attractor)
  Synergy that is not explained by either drug alone.
  TNBC: PARP inhibitor (Type 1) + EZH2 inhibitor
  (Type 2) — this combination is in active trials.
  The composite type predicts the combination.
  The single-type analysis would not.

PREDICTION D — No Type 3 cancer will show classic
  oncogene activation as the primary driver:
  Type 3 geometry requires the identity programme
  to be intact or over-expressed. Classic oncogene
  activation (MYC amplification, RAS mutation,
  HER2 amplification) disrupts identity programmes
  and would push toward Type 2 geometry.
  Therefore: if a cancer is Type 3, its primary
  driver should be TUMOR SUPPRESSOR LOSS rather
  than ONCOGENE ACTIVATION.
  LumA: TGFBR2 loss, SMAD3/4 loss, CDKN1A loss.
  All tumor suppressors. No oncogene activation.
  This is a testable cross-type prediction.

PREDICTION E — Drug resistance mechanisms differ by type:
  Type 1 resistance: alternative lock reforms
                     (new epigenetic silencing of
                     switch gene)
  Type 2 resistance: new convergence node replaces
                     inhibited node (EZH2i resistance
                     via alternative PRC2 subunit)
  Type 3 resistance: alternative cycling pathway
                     bypasses restored/compensated
                     arrest (CDK4/6i resistance via
                     CCNE1 amplification, which drives
                     CDK2 instead of CDK4/6)
  These resistance mechanisms are TYPE-SPECIFIC.
  Predicting them before resistance develops allows
  rational combination therapy design.
```

---

## PART VI — FRAMEWORK LESSONS UPDATED

```
LESSON 17 (from cross-cancer taxonomy, 2026-03-04):
  Not all false attractors are geometrically identical.
  The three types require different analytical frames
  before data is loaded.
  The Identity TF Direction Test (suppressed vs elevated)
  is the fastest pre-data discriminator.
  Apply it before writing any before-document.
  If differentiation TFs are elevated:
    You are likely in Type 3.
    Switch the frame before predicting.

LESSON 18 (from cross-cancer taxonomy, 2026-03-04):
  The ligand paradox (arrest ligand elevated +
  receptor absent) is the defining diagnostic
  of Type 3 geometry beyond simple TF direction.
  If the arrest ligand is elevated but the pathway
  is suppressed — the cell is producing the signal
  to stop but cannot hear it.
  This finding changes the drug logic from
  "compensate for absent arrest" to
  "restore the missing receptor."
  Look for this pattern in every Type 3 candidate.

LESSON 19 (from cross-cancer taxonomy, 2026-03-04):
  Composite type cancers predict combination therapy.
  If a cancer shows evidence of Type 1 → Type 2
  sequence, the correct treatment is the combination
  of both type-specific drugs.
  Single-drug strategies that address only one stage
  will have limited efficacy in composite type cancers.
  Map the type sequence before designing the
  therapeutic prediction.
```

---

## STATUS AND VERSION HISTORY

```
document:           ATTRACTOR_GEOMETRY_AXIOMS.md
document_number:    90
version:            1.0
date:               2026-03-04
author:             Eric Robert Lawson
                    OrganismCore
status:             LIVING DOCUMENT

version_history:
  1.0  2026-03-04  Initial derivation
                   Three types named and specified
                   Derived from 14 cancer validations
                   Confirmed examples listed per type
                   Cross-type comparison table built
                   Diagnostic algorithm stated
                   Novel predictions generated (A-E)
                   Framework lessons 17-19 added

to_be_updated_when:
  - New cancer analysis produces a Type 3 example
    beyond LumA (update confirmed examples)
  - ILC analysis completes (confirm/revise Type 3
    architectural variant prediction)
  - Composite type analysis produces new example
  - A fourth structural type is identified
  - Any prediction A-E is confirmed or falsified

instructions_for_update:
  Do not revise confirmed examples retroactively.
  Add new examples with date.
  If a prediction is confirmed: mark ✓ + date + source.
  If a prediction is falsified: mark ✗ + date +
    record what the data showed instead +
    derive the corrected principle.
  Wrong predictions are as important as confirmed ones.
  The document is a living record, not a polished claim.

next_analysis_to_test_these_axioms:
  TNBC (BRCA-S4a)
    — Test composite type hypothesis (Type 1 → Type 2)
    — Test whether BRCA1 loss is the Type 1 block
      that precedes the Type 2 basal false attractor
    — Test whether PARP + EZH2 combination is
      predicted by composite type logic
```
