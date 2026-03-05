# LUMINAL B BREAST CANCER — BEFORE DOCUMENT
## Predictions Locked Before Any Data Loads
## OrganismCore — Document BRCA-S5a
## Date: 2026-03-05

---

## DOCUMENT METADATA

```
document_id:        BRCA-S5a
series:             BRCA Deep Dive — Luminal B
folder:             Cancer_Research/BRCA/DEEP_DIVE/LUMINAL_B/
type:               BEFORE-DOCUMENT
                    All predictions stated before
                    any data is loaded.
                    This document cannot be modified
                    after Script 1 runs.
date:               2026-03-05
author:             Eric Robert Lawson
                    OrganismCore
protocol_version:   Workflow_Protocol.md v2.0
dataset_planned:    GSE176078 — Wu et al. 2021
                    (same scRNA-seq dataset as LumA and TNBC —
                    Luminal B cancer cells extracted by
                    celltype_subset annotation)
                    + TCGA-BRCA PAM50 LumB subset
                    (bulk RNA-seq, n~200-250)
                    for depth score / clinical correlation
status:             LOCKED — predictions cannot change
                    after this document is committed
precursor_documents:
  ATTRACTOR_GEOMETRY_AXIOMS.md (Doc 90)
  BRCA_Subtype_Orientation.md
  BRCA-S1a (LumA predictions)
  BRCA-S1b (LumA Script 1 artifact)
  BRCA-S2a (LumA Script 2 predictions)
  BRCA-S2b (LumA Script 2 artifact)
  BRCA-S2c (LumA Literature Check)
  BRCA-S4a (TNBC predictions)
```

---

## WHAT IS ALREADY KNOWN FROM PRIOR ANALYSES

```
This is NOT a blank slate.
Four subtypes have been completed before this one.
The following are established facts, not predictions:

FROM LUMINAL A (BRCA-S1b, BRCA-S2b):

  LUMINAL A GEOMETRY:
    FOXA1:   +37.3% vs Mature Luminal (p=3.61e-13)
             +1038% vs Luminal Progenitor
    GATA3:   +33.8% vs Mature Luminal (p=9.68e-13)
    ESR1:    -30.1% vs Mature Luminal (ns — retained)
    CDKN1A:  -74.3% vs Mature Luminal (p=4.68e-195)
    TGFBR2:  -97%   vs Luminal Progenitor
    EZH2:    +18.5% vs Mature Luminal (ns — flat)
    PGR:     -54.8% vs Mature Luminal (p=3.11e-22)

  LUMINAL A ATTRACTOR TYPE: TYPE 3
    Identity programme RETAINED and AMPLIFIED.
    Cell is in the correct Waddington valley.
    The problem is GROWTH REGULATION, not identity.
    The depth axis in LumA is the ARREST axis (CDKN1A).
    Deeper LumA = less p21 = more CDK4 activity.

  LUMINAL A DEPTH SCORE:
    1 - norm(CDKN1A) + norm(CDK4) / 2
    CDK4     r=+0.808  within LumA
    GATA3    r=+0.384  within LumA
    TGFBR2   -97% vs progenitor — upstream of CDKN1A

FROM TNBC (BRCA-S4a — predictions only, not yet run):
  TNBC geometry: Composite Type 1 → Type 2
  Luminal identity TFs SUPPRESSED
  Basal false attractor ACTIVATED
  EZH2 elevated as convergence node

THE STRUCTURAL QUESTION THIS ANALYSIS ANSWERS:
  Luminal A = terminally differentiated luminal cell
              that RETAINED identity but lost arrest
              (Type 3 — in the right valley,
               growth regulation failed)

  Luminal B = luminal progenitor that PARTIALLY
              differentiated but never fully committed
              to the mature luminal programme
              AND acquired proliferative drive

  These are different positions on the same landscape.

  The Luminal B analysis answers:
    1. How much luminal identity is retained vs LumA?
    2. Is there a measurable differentiation gap —
       FOXA1/GATA3/ESR1 intermediate between
       LumA and TNBC?
    3. Is the depth axis the same (arrest/CDKN1A)?
       Or different (partial identity loss)?
    4. Does EZH2 rise above the LumA flat baseline?
       LumA: EZH2 +18.5% ns
       TNBC: EZH2 large elevation (predicted)
       LumB: EZH2 intermediate? This is the key test.
    5. Does the proliferative drive (Ki-67 high)
       correspond to a different attractor geometry
       or simply a different depth on the same axis?

THIS IS THE CROSS-SUBTYPE STRUCTURAL TEST:
  Luminal B sits between Luminal A and TNBC on
  the breast Waddington landscape.
  Every prediction below derives from that position.
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
  a less committed cell than the mature luminal
  cell that is the LumA cell of origin.
  The correct terminal destination is the mature
  luminal epithelial cell (same destination as LumA,
  different starting point).

QUESTION 2: Are identity TFs of the correct
  terminal state expressed?
  PARTIAL — Luminal B is ER+ (ESR1 expressed).
  FOXA1 and GATA3 are expressed — but at lower
  levels than in LumA.
  The luminal identity programme is PRESENT but
  INCOMPLETELY EXECUTED.

QUESTION 3: Is the cell progressing toward the
  correct terminal state?
  NO — the cell is proliferating, not differentiating.
  Ki-67 high (>20%). Grade 2-3.
  The cell has arrested partway through luminal
  differentiation and is cycling instead of
  completing the luminal programme.

ATTRACTOR TYPE ASSIGNMENT:

  PRIMARY TYPE: TYPE 3 — SAME VALLEY, DIFFERENT DEPTH
  The cell is in the luminal valley (ER+, FOXA1+,
  GATA3+). It has NOT fallen into a different valley
  (not TNBC, not basal-like). This is the same
  Waddington valley as LumA.

  CRITICAL DISTINCTION FROM LumA (STATED AS PREDICTION):
  The LumA geometry showed identity AMPLIFIED above
  mature luminal (FOXA1 +37%, GATA3 +34%) and
  arrest dismantled (CDKN1A -74%).
  The cell had OVERSHOOT on identity AND arrested
  growth control simultaneously.

  LumB PREDICTED GEOMETRY: PARTIAL EXECUTION
  Identity TFs will be LOWER than LumA:
    FOXA1:  LumB < LumA
    GATA3:  LumB < LumA
    ESR1:   LumB < LumA (ER lower — clinically confirmed
            by lower PR expression in LumB vs LumA)
  Arrest axis will be compromised but differently:
    The depth axis in LumB may be DUAL:
    Component 1: Identity gap (FOXA1/GATA3 lower)
    Component 2: Proliferative drive (Ki-67/CCND1 higher)
    Unlike LumA where identity was amplified and
    arrest was dismantled independently.

  COMPOSITE EVALUATION:
  Is LumB a composite type?
    PARTIAL TYPE 1 COMPONENT:
      The cell never fully differentiated from the
      luminal progenitor. Progenitor identity partially
      retained (partial Type 1 — differentiation
      incomplete).
    TYPE 3 DOMINANT:
      The cell is in the correct valley but at a
      different depth. Not a wrong-valley problem.

  FINAL TYPE ASSIGNMENT: TYPE 3 (DOMINANT)
  with PARTIAL TYPE 1 COMPONENT
  (incomplete differentiation from progenitor state)

  This is a DIFFERENT COMPOSITE than TNBC.
  TNBC: Type 1 → Type 2 (progenitor → wrong valley)
  LumB: Type 1 partial → Type 3 (progenitor →
        correct valley, incomplete execution)
```

---

## PART II — BIOLOGICAL GROUNDING
## The lineage and its normal architecture

```
THE NORMAL JOURNEY THIS CANCER INTERRUPTS:

  Mammary stem cell
    ↓
  Luminal progenitor   ← WHERE LumB ORIGINATES
    ↓                    (LumA originates one step lower)
  Committed luminal    ← Where LumB partially stalls
    ↓                    (incomplete commitment)
  Mature luminal       ← Correct terminal state
  epithelial cell        (ESR1+, FOXA1+, GATA3+ HIGH)

THE DIFFERENCE BETWEEN LumA AND LumB CELLS OF ORIGIN:

  LumA cell of origin: MATURE LUMINAL CELL
    Fully differentiated. High ESR1, very high FOXA1,
    very high GATA3. Slow cycling (Ki-67 low).
    The arrest programme (CDKN1A) is the primary
    barrier that was dismantled.

  LumB cell of origin: LUMINAL PROGENITOR
    Partially differentiated. ESR1 expressed but
    lower. FOXA1 and GATA3 present but at progenitor
    levels, not mature luminal levels.
    BOTH the differentiation completion AND the
    arrest programme are compromised.

THE FALSE ATTRACTOR STATE IN LumB:
  LumB is NOT a false attractor in the Type 2 sense.
  The cell is in the luminal valley.
  The problem is the depth within that valley —
  the cell is higher up the slope (less committed)
  than LumA, closer to the progenitor starting point.

  The "false attractor" in LumB is a proliferating
  luminal progenitor state:
    ER+, FOXA1+, GATA3+ — but at progenitor levels
    Ki-67 HIGH (cycling)
    CCND1 HIGH (cyclin D1 — drives G1/S transition)
    p21 (CDKN1A) LOW (arrest dismantled as in LumA)
    CDK4/6 ACTIVE (downstream of CCND1, unchecked
    without p21)

  EZH2 STATUS (KEY PREDICTION):
    LumA: EZH2 flat (+18.5%, ns)
    TNBC: EZH2 highly elevated
    LumB: EZH2 INTERMEDIATE — predicted elevated above
    LumA but below TNBC.
    The partial differentiation block requires some
    epigenetic silencing of the mature luminal programme
    (EZH2 must be doing something) but the cell has
    not undergone the complete epigenetic lock seen
    in TNBC.

THE SADDLE POINT:
  The saddle point in LumB is HIGHER in the landscape
  than in LumA — the cell is closer to the luminal
  progenitor branch point.
  The transition from normal progenitor → LumB false
  state is a shallower barrier than in TNBC.
  This is why LumB is:
    — More chemotherapy-responsive than LumA
      (less committed = more vulnerable to cytotoxics)
    — Less endocrine-responsive as monotherapy than LumA
      (lower ER circuit fidelity)
    — More genomically unstable than LumA
      (less committed identity = less chromatin stability)
```

---

## PART III — PREDICTIONS
## All stated 2026-03-05 before data loads

---

### P1 — SWITCH GENE PANEL vs LumA AND vs TNBC

```
The key structural test for LumB:
  Do identity TFs fall between LumA and TNBC?
  Or do they match LumA?

This is the prediction that most directly tests
the "intermediate position" hypothesis.

COMPARISON 1: LumB vs Mature Luminal
  (same comparison as LumA Script 1)

  FOXA1:  PREDICTED RETAINED but LOWER than LumA
          LumA was +37.3% above mature luminal.
          LumB predicted: +10–20% above mature luminal
          OR near-equal to mature luminal
          NOT suppressed (still in luminal valley)
          Direction: UP or neutral  p predicted: ns to <0.05

  GATA3:  PREDICTED RETAINED but LOWER than LumA
          LumA was +33.8% above mature luminal.
          LumB predicted: near-equal to or slightly
          above mature luminal (not the overshoot
          seen in LumA)
          Direction: neutral to slightly UP

  ESR1:   PREDICTED RETAINED but REDUCED vs LumA
          Lower PR expression in LumB is clinically
          established (PR often lower in LumB than LumA).
          PR is an ER target gene (PGR).
          If ESR1 circuit fidelity is lower in LumB,
          PGR should be lower than in LumA.
          ESR1 itself may be near LumA level
          (ER still present by definition) but
          ESR1 circuit OUTPUT (PGR) predicted lower.
          Direction: ESR1 neutral, PGR DOWN vs LumA

  SPDEF:  PREDICTED LOWER than LumA
          Secretory/luminal TF downstream of FOXA1.
          LumB is less committed to the secretory
          mature luminal programme.
          Direction: DOWN vs LumA

  COMPARISON: LumB vs LumA (direct structural test)
    FOXA1:  LumB < LumA  (less identity overshoot)
    GATA3:  LumB < LumA  (less identity overshoot)
    ESR1:   LumB ≈ LumA  (both ER+)
    PGR:    LumB < LumA  (less ER circuit output)
    CDKN1A: LumB ≈ LumA  (both arrest-dismantled)
            OR LumB < LumA (even more depleted given
            higher Ki-67)

  These directional predictions are the primary test
  of the intermediate-position model.
  If FOXA1 and GATA3 are equal to LumA:
    LumB and LumA are at the same depth —
    the difference is only in growth control.
  If FOXA1 and GATA3 are lower in LumB than LumA:
    LumB is at a genuinely higher position on
    the slope — less differentiated, not just
    less growth-controlled.
```

---

### P2 — PROLIFERATIVE AXIS

```
Clinically, LumB is defined by HIGH Ki-67 (>20%).
This is the most robust clinical marker separating
LumB from LumA.

PREDICTED ELEVATED in LumB vs LumA:

  MKI67:  Ki-67 proliferation marker
          Expected: significantly higher in LumB than LumA
          The defining clinical distinction.
          Direction: UP vs LumA  p predicted: <0.001

  CCND1:  Cyclin D1 — drives G1/S transition
          In LumA: CCND1 was flat (+6.4%, ns)
          In LumB: CCND1 predicted elevated
          (more active cell cycle at the G1/S checkpoint)
          Direction: UP vs LumA  p predicted: <0.01

  CDK4:   Cyclin-dependent kinase 4
          In LumA: CDK4 r=+0.808 within-LumA depth axis
          In LumB: CDK4 predicted elevated vs LumA
          (more active CDK4/6 signaling = higher Ki-67)
          Direction: UP vs LumA  p predicted: <0.05

  MYC:    Master proliferation TF
          In LumA: MYC -44.5% (suppressed)
          In LumB: predicted near-neutral or elevated
          MYC is a known driver in LumB specifically
          (TCGA data shows higher MYC activity in LumB)
          Direction: UP vs LumA  p predicted: <0.05

  TOP2A:  Topoisomerase II alpha — DNA replication
          Highly cell-cycle regulated.
          Elevated in high Ki-67 tumors.
          Direction: UP vs LumA  p predicted: <0.01

  CDKN1A: p21 — arrest gene
          In LumA: CDKN1A -74.3% (the primary finding)
          In LumB: predicted EQUALLY or MORE depleted
          Given higher Ki-67 in LumB, p21 suppression
          must be at least as deep as in LumA.
          Direction: DOWN or equal vs LumA
```

---

### P3 — EZH2 POSITION (THE CRITICAL BETWEEN-SUBTYPE TEST)

```
EZH2 levels across the luminal subtypes is the
most important structural prediction in this document.

ESTABLISHED:
  LumA:   EZH2 +18.5% vs Mature Luminal (ns — flat)
          EZH2 is NOT the LumA lock.
  TNBC:   EZH2 predicted highly elevated (convergence node)
          (BRCA-S4a — not yet tested)

PREDICTION FOR LumB:

  EZH2 ELEVATED vs LumA  (not flat as in LumA)
  EZH2 LOWER than TNBC   (not the dominant lock as in TNBC)

  Reasoning:
    LumB cells are luminal progenitors that have
    not fully differentiated. Some degree of EZH2-mediated
    chromatin silencing of the mature luminal programme
    is required to maintain the progenitor state.
    In LumA (terminally differentiated), EZH2 is not
    needed to suppress differentiation — the cell
    already completed differentiation before the
    malignant transformation.
    In TNBC, EZH2 is the full epigenetic lock on
    the basal false attractor.
    In LumB: INTERMEDIATE.

  FORMAL PREDICTION:
    EZH2 in LumB > EZH2 in LumA (above flat baseline)
    EZH2 in LumB < EZH2 in TNBC (not the dominant lock)
    Magnitude vs Mature Luminal: +30–60% predicted
    p predicted: <0.05

  IF CONFIRMED:
    This is the most direct evidence that LumB
    sits between LumA and TNBC on the epigenetic
    axis of the breast Waddington landscape.
    EZH2 elevation would predict:
      CDK4/6 inhibitor response (as in LumA)
      BUT ALSO partial EZH2 inhibitor sensitivity
      (unlike LumA where EZH2i would have no target)

  IF FLAT (equal to LumA):
    LumB and LumA share the same epigenetic state.
    The difference between them is purely in
    growth regulation, not differentiation depth.
    This would simplify the therapeutic picture
    (CDK4/6 inhibitors serve both).

  IF HIGHLY ELEVATED (equal to TNBC):
    LumB has an epigenetic component beyond LumA.
    EZH2 inhibition would be predicted to benefit
    LumB specifically — a subtype-specific target
    that is currently not approved for LumB.
    This would be a novel finding.

  PREDICTION RANKING:
    Most likely: EZH2 intermediate (between LumA and TNBC)
    Second: EZH2 flat (equal to LumA)
    Third: EZH2 highly elevated (equal to TNBC)
```

---

### P4 — DEPTH SCORE CONSTRUCTION

```
The depth score for LumB must address the
DUAL NATURE of the predicted geometry:

  Component 1 — PROLIFERATIVE AXIS (same as LumA)
    1 - norm(CDKN1A) + norm(CDK4) / 2
    This is the LumA depth score.
    It should work in LumB too — but may capture
    less variance because LumB has higher baseline
    proliferation (less dynamic range within LumB
    than within LumA for the CDKN1A axis).

  Component 2 — IDENTITY AXIS (new for LumB)
    1 - norm(mean of FOXA1, GATA3, ESR1)
    Higher = less committed to mature luminal identity
    This component was FLAT in LumA because identity
    was amplified (FOXA1/GATA3 elevated, so this
    component = 0 for most LumA cells).
    In LumB: identity TFs are lower → this component
    has non-zero variance → it may contribute to depth.

  PROPOSED LumB DEPTH SCORE:
    Depth = w1 * (1 - norm(CDKN1A)) +
            w2 * (1 - norm(mean(FOXA1, GATA3))) +
            w3 * norm(MKI67)
    All three components normalized.
    Weights: w1=0.4, w2=0.4, w3=0.2
    (CDKN1A and identity axes equal weight;
    Ki-67 as secondary component)

  IF THE IDENTITY AXIS HAS ZERO VARIANCE IN LumB:
    LumB depth = LumA depth formula
    (CDKN1A-dominated, same biology)
    This would mean: LumB and LumA are in the
    same attractor at different depths.

  IF THE IDENTITY AXIS HAS SIGNIFICANT VARIANCE:
    LumB depth is a two-component score.
    This would mean: LumB contains a range of cells
    from nearly-LumA (full luminal identity, high
    FOXA1/GATA3) to partial progenitor (lower
    FOXA1/GATA3, higher MKI67).
    This is consistent with the known intratumoral
    heterogeneity of LumB tumors.

PREDICTED TOP DEPTH CORRELATES within LumB:
  CDK4:     r > +0.40  (proliferation axis — as in LumA)
  MKI67:    r > +0.30  (Ki-67 rises with depth)
  CCND1:    r > +0.25  (cyclin D1 with depth)
  EZH2:     r > +0.20  (epigenetic component)
  FOXA1:    r uncertain (UP in LumA depth; may be
            UP or DOWN in LumB depending on which
            component dominates — FLAG FOR DATA)
  CDKN1A:   r < -0.40  (p21 falls with depth — as LumA)
  ESR1:     r < -0.15  (ER circuit weaker at depth)
```

---

### P5 — ER CIRCUIT INTEGRITY

```
The ER circuit in LumA showed a BREAK:
  ESR1 retained (+, near normal)
  PGR (ER target) SUPPRESSED (-54.8%)
  NCOA1/2 (coactivators) SUPPRESSED
  TGFBR2 (upstream) -97%
  The circuit was present but the coactivator
  machinery and upstream TGF-β signalling were gone.

PREDICTION FOR LumB ER CIRCUIT:

  The ER circuit break will be MORE PRONOUNCED
  in LumB than in LumA.

  Reasoning:
    LumB is less ER-dependent than LumA clinically
    (endocrine therapy less effective as monotherapy).
    Higher Ki-67 despite ER expression indicates
    the ER-dependent growth arrest programme is
    less functional.
    PR is clinically lower in LumB than LumA
    (this is a known clinical fact — used as one
    criterion for LumA vs LumB classification in
    St. Gallen guidelines).

  SPECIFIC PREDICTIONS:
    PGR:    LOWER in LumB than LumA
            In LumA: -54.8% vs Mature Luminal
            In LumB: predicted -60 to -75% vs
            Mature Luminal (more suppressed)
            Direction: DOWN  p predicted: <0.001

    NCOA1:  Predicted SUPPRESSED or lower than LumA
            Coactivator gap more pronounced.

    TGFBR2: Predicted LOW (similar to LumA -97%)
            TGF-β blindness is a shared feature
            of the luminal cancer landscape.
            Direction: DOWN  p predicted: <0.001

    r(ESR1, PGR) within LumB:
            Predicted LOWER than in LumA
            (circuit gap larger in LumB)
            The ER→PR link is weaker.
```

---

### P6 — WHAT DISTINGUISHES LumB FROM LumA AT THE DATA LEVEL

```
This section states the DIFFERENTIAL PREDICTIONS:
What should be DIFFERENT between LumB and LumA
in the gene expression data.

These are the key discriminating tests.
If these differentials are not found, the two
subtypes are indistinguishable by this framework.

PREDICTED HIGHER IN LumB vs LumA:
  MKI67    Ki-67 (defining clinical marker)
  CCND1    Cyclin D1 (cell cycle driver)
  CDK4     CDK4 kinase (depth axis gene in LumA)
  MYC      MYC (oncogene, more active in LumB)
  TOP2A    DNA replication enzyme
  EZH2     Epigenetic lock (more active than in LumA)
  ERBB2    HER2 (Luminal B HER2+ variant has ERBB2
           amplification — even in HER2- LumB, ERBB2
           may be slightly elevated vs LumA)

PREDICTED LOWER IN LumB vs LumA:
  FOXA1    Pioneer TF (less identity overshoot)
  GATA3    Luminal identity TF (less overshoot)
  PGR      PR target of ER (less ER circuit fidelity)
  SPDEF    Secretory TF (less mature luminal commitment)
  AGR2     FOXA1 target (less secretory programme)
  CDKN1A   p21 (predicted equal or lower — both depleted)

PREDICTED EQUAL:
  ESR1     ER alpha (both ER+)
  KRT8/18  Luminal cytokeratins (both luminal)
  KRT5/14  Basal markers (absent in both)
  CDX2     Non-breast control (absent in both)

THE CRITICAL TEST (single most important):
  FOXA1: LumB < LumA?
  If YES: LumB is genuinely less differentiated.
          The intermediate position model holds.
  If NO (FOXA1 equal): LumB and LumA are in the
          same differentiation state. The difference
          is only in growth regulation (Ki-67).
          The intermediate position model fails.
  Record outcome explicitly in BRCA-S5b.
```

---

### P7 — EPIGENETIC PREDICTIONS

```
Based on the LumA epigenetic findings and the
intermediate-position model:

EZH2:
  LumA:   +18.5% vs Mature Luminal (ns, flat)
  LumB:   Predicted +30–60% vs Mature Luminal
          p predicted: <0.05
  Reasoning: partial differentiation block requires
  some EZH2 activity at mature luminal TF loci.

HDAC1/2:
  LumA:   SUPPRESSED (part of the LumA geometry —
          global epigenetic simplification)
  LumB:   Predicted FLAT OR SLIGHTLY ELEVATED vs LumA
  Reasoning: LumB is less fully committed to any
  identity, so the epigenetic landscape is less
  simplified — HDAC may be more active to maintain
  the progenitor chromatin state.

KDM1A (LSD1):
  LumA:   SUPPRESSED
  LumB:   Predicted near-LumA (both luminal valley
          cells — the LumA result may generalize)

DNMT3A/3B:
  Predicted ELEVATED in LumB vs LumA
  LumB has more extensive promoter methylation
  at mature luminal TF targets (to block completion
  of differentiation). This is a prediction of
  elevated de novo methylation activity.
  Direction: UP vs LumA  p predicted: <0.05

NOTE:
  The epigenetic predictions for LumB are less
  certain than for TNBC because the LumA result
  was surprising (EZH2 flat, HDAC suppressed).
  The field should report what it finds without
  forcing the intermediate model.
  Flag all epigenetic results as HIGH INTEREST
  regardless of direction.
```

---

### P8 — DRUG TARGET PREDICTIONS

```
All stated 2026-03-05 before Script 1 runs.

DRUG TARGET 1 — CDK4/6 INHIBITORS
  Palbociclib, ribociclib, abemaciclib
  Mechanism: Block CDK4/6 → restore G1 arrest
  compensating for CDKN1A depletion
  In LumA: CDK4/6 inhibitors are approved and
  standard of care in metastatic disease.
  In LumB: SAME MECHANISM — also approved
  and effective.
  PREDICTION: CDK4/6 inhibitors are effective in
  LumB for the same reason as LumA: CDKN1A depletion
  leaves CDK4/6 unchecked, and pharmacological
  inhibition restores the arrest signal.
  Status: ✓ CONVERGENT (already approved)

DRUG TARGET 2 — CDK4/6 INHIBITOR BENEFIT MAGNITUDE
  NOVEL PREDICTION:
  Within LumB, depth score should predict CDK4/6
  inhibitor benefit magnitude.
  Deeper LumB (lower CDKN1A, higher CDK4) =
  more CDK4/6 dependent = larger benefit from
  CDK4/6 inhibition.
  This is the same NOVEL-1 prediction made for LumA
  (BRCA-S2c NOVEL-1).
  If it holds in LumB as well: it is a pan-luminal
  prediction (not LumA-specific).
  Status: 🆕 NOVEL extension from LumA geometry

DRUG TARGET 3 — EZH2 INHIBITORS
  Tazemetostat
  In LumA: EZH2 flat → EZH2i predicted ineffective
  in LumA specifically.
  In TNBC: EZH2 elevated → EZH2i predicted effective.
  In LumB: IF EZH2 elevated (P3 prediction):
    EZH2 inhibition may have partial activity in LumB
    that it lacks in LumA.
    This would mean LumB (not LumA) is the luminal
    subtype that benefits from EZH2i.
  Status: 🆕 NOVEL — depends on P3 confirmation
  High interest if EZH2 intermediate elevation confirmed.

DRUG TARGET 4 — ENDOCRINE THERAPY + CDK4/6 COMBINATION
  The standard of care for metastatic LumB is
  endocrine therapy + CDK4/6 inhibitor.
  PREDICTION: The depth score explains WHY the
  combination is required where monotherapy fails.
  Monotherapy (endocrine alone): fails because
    the ER circuit gap (TGFBR2/coactivator depletion)
    means ER cannot re-engage the arrest programme
    even when activated by ligand manipulation.
  CDK4/6 inhibitor required: to bypass the broken
    TGF-β → p21 → CDK4/6 circuit directly.
  This is the MECHANISTIC EXPLANATION for the
  standard of care, derived from the geometry.
  Status: ✓ EXPLANATORY (derives mechanism of
  approved combination from attractor geometry)

DRUG TARGET 5 — HER2-TARGETED THERAPY (LumB HER2+ variant)
  ~15–20% of LumB is HER2+.
  In these tumors: both endocrine and anti-HER2
  therapy are standard.
  PREDICTION: In the HER2+ LumB subgroup, ERBB2
  expression will be elevated in the expression
  data. The depth score for this subgroup will be
  higher than HER2- LumB (ERBB2 amplification
  pushes cells deeper in the attractor by adding
  a second proliferative driver beyond CCND1/CDK4).
  Status: ⚠ PARTIAL — the scRNA-seq data may not
  have sufficient HER2+ LumB cells to test this.
  Flag as contingent on sample composition.
```

---

### P9 — CROSS-SUBTYPE STRUCTURAL PREDICTIONS

```
This is the cross-subtype prediction that this
analysis is specifically designed to test.
Stated here as prediction (to be tested in
subsequent cross-subtype analysis):

THE ER AXIS GRADIENT:
  ESR1 expression should form a gradient:
    LumA  >  LumB  >  HER2-enriched  >  TNBC
  Specifically:
    LumA: ESR1 near-normal (retained luminal identity)
    LumB: ESR1 lower than LumA
    HER2-enriched: ESR1 near-absent (ER-)
    TNBC: ESR1 absent

  THE FOXA1/GATA3 GRADIENT:
    LumA:  FOXA1/GATA3 highest (overshoot above normal)
    LumB:  FOXA1/GATA3 intermediate
    HER2:  FOXA1/GATA3 low
    TNBC:  FOXA1/GATA3 near-absent

  IF CONFIRMED: the breast cancer landscape has a
  single continuous axis — the luminal identity axis —
  along which all four ER-relevant subtypes can be
  ordered. LumB is the second point on this axis.
  This is the most important structural claim
  in the entire BRCA series.

  THE EZH2 GRADIENT (inverse of identity):
    LumA:  EZH2 flat
    LumB:  EZH2 intermediate
    TNBC:  EZH2 high
  If the EZH2 gradient mirrors the FOXA1/GATA3
  gradient in reverse, the epigenetic lock and
  the identity retention are inversely coupled
  across the entire luminal landscape.
  This is the breast Waddington energy surface
  with EZH2 as the barrier height and
  FOXA1/GATA3 as the valley position.

CONTROLS:
  The following should be absent in LumB
  (as they were in LumA):
    KRT5, KRT14  (basal markers — absent in LumB)
    SOX10        (neural crest / basal — absent)
    CDX2         (intestinal lineage — absent)
    NKX2-1       (lung lineage — absent)
  If any are elevated: analyst assumption error.
  Record and process per Wrong Prediction Protocol.
```

---

## PART IV — WHAT THIS ANALYSIS CANNOT TELL US

```
The scRNA-seq data (GSE176078) may or may not
contain a labeled "Cancer LumB SC" population.

CRITICAL DATA NOTE:
  In GSE176078, the cell annotations include:
    "Cancer LumA SC"   — confirmed present, n=7,742
    "Cancer Basal SC"  — confirmed present, n=4,312
  It is UNCERTAIN whether there is a
  "Cancer LumB SC" label or equivalent.

  If no LumB-specific label exists:
    Option A: Use TCGA-BRCA PAM50 LumB bulk data
              as the primary dataset.
    Option B: Extract LumB cells from GSE176078
              using ERBB2 expression and Ki-67
              analog (MKI67) within ER+ cells
              to enrich for a LumB-like population.
    Option C: Use a separate scRNA-seq dataset
              with LumB annotation.
    The script must assess this at load time
    and document the approach taken.

  The TCGA-BRCA bulk analysis will provide
  the depth score / clinical correlation
  regardless of scRNA-seq availability.

LIMITATIONS:
  1. If LumB label absent from GSE176078:
     single-cell comparison to LumA will require
     proxy definition — lower methodological rigor
  2. No pCR annotation available in bulk data
     for LumB (pCR not standard in luminal disease —
     neoadjuvant chemotherapy rarely used in LumA/B)
  3. The HER2+ LumB subgroup may be too small
     for separate analysis
  4. LumB is split by PAM50 — the clinical Ki-67
     cutoff misclassifies ~20% of cases
     (this affects any dataset using clinical
     rather than PAM50 calls)
```

---

## PART V — COMPLETE PREDICTION REFERENCE

```
PREDICTIONS LOCKED 2026-03-05
Document: BRCA-S5a
Author: Eric Robert Lawson, OrganismCore

ATTRACTOR TYPE: TYPE 3 (dominant) with PARTIAL TYPE 1
  (intermediate position — luminal progenitor origin,
  correct valley, incomplete differentiation execution)

P1 — SWITCH GENES vs LumA:
  FOXA1: LumB < LumA  (less identity overshoot)
  GATA3: LumB < LumA  (less identity overshoot)
  ESR1:  LumB ≈ LumA  (both ER+)
  PGR:   LumB < LumA  (less ER circuit output)
  SPDEF: LumB < LumA  (less secretory commitment)
  Critical test: FOXA1 direction vs LumA

P2 — PROLIFERATIVE AXIS (all UP vs LumA):
  MKI67, CCND1, CDK4, MYC, TOP2A
  CDKN1A: DOWN or equal vs LumA

P3 — EZH2 POSITION (CRITICAL TEST):
  EZH2: INTERMEDIATE between LumA (flat) and TNBC (high)
  Predicted: +30–60% vs Mature Luminal, p<0.05
  r(EZH2, depth) > +0.20 predicted

P4 — DEPTH SCORE:
  Dual component:
    1 - norm(CDKN1A) [proliferative axis]
    + 1 - norm(mean(FOXA1, GATA3)) [identity axis]
  Top correlates: CDK4, MKI67, CCND1 (positive)
                  CDKN1A, FOXA1, GATA3, PGR (negative)
  Flag FOXA1 direction as most important single result.

P5 — ER CIRCUIT:
  PGR more suppressed in LumB than LumA
  r(ESR1, PGR) within LumB < r in LumA
  TGFBR2 similarly depleted as in LumA

P6 — DIFFERENTIAL LumB vs LumA:
  Higher: MKI67, CCND1, CDK4, MYC, EZH2, ERBB2
  Lower: FOXA1, GATA3, PGR, SPDEF, AGR2
  Equal: ESR1, KRT8/18

P7 — EPIGENETIC:
  EZH2: UP vs LumA (+30–60%)
  DNMT3A/3B: UP vs LumA (de novo methylation)
  HDAC1/2: FLAT or slight UP vs LumA
  KDM1A: near-LumA

P8 — DRUG TARGETS:
  ✓ CDK4/6 inhibitors (same mechanism as LumA)
  🆕 Depth score predicts CDK4/6 benefit within LumB
  🆕 EZH2 inhibitors (if P3 confirmed — LumB-specific)
  ✓ Endocrine + CDK4/6 combination (mechanism explained)
  ⚠ HER2-targeted (LumB HER2+ variant only)

P9 — CROSS-SUBTYPE:
  ESR1 gradient: LumA > LumB > HER2 > TNBC
  FOXA1/GATA3 gradient: LumA > LumB > HER2 > TNBC
  EZH2 gradient (inverse): LumA < LumB < TNBC
  If confirmed: single continuous luminal identity axis

CONTROLS (should be absent in LumB):
  KRT5, KRT14, SOX10, CDX2, NKX2-1, NKX3-1
  If any elevated: analyst assumption error
  Record per Wrong Prediction Protocol
```

---

## STATUS BLOCK

```
document:           BRCA-S5a
type:               Before-Document (locked predictions)
date:               2026-03-05
author:             Eric Robert Lawson / OrganismCore
status:             LOCKED

predictions_count:  9 prediction groups
attractor_type:     Type 3 (dominant) + partial Type 1
datasets:           GSE176078 (scRNA-seq — LumB label TBC)
                    TCGA-BRCA PAM50 LumB bulk (primary)
next_document:      BRCA-S5b (Script 1 reasoning artifact)
script:             BRCA_LumB_script1.py

structural_position:
  LumB sits between LumA and TNBC on the breast
  Waddington landscape. Every prediction in this
  document derives from that single claim.
  The data will confirm or refute it.
  No exceptions — record all discrepancies.

cross_subtype_flag:
  After BRCA-S5b is complete, the LumA vs LumB
  comparison (FOXA1/GATA3/EZH2 gradient test)
  is the first cross-subtype analysis in the
  BRCA series. It must be done before Claudin-low
  begins.
```
