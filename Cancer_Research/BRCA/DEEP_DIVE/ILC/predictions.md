# INVASIVE LOBULAR CARCINOMA — BEFORE DOCUMENT
## Predictions Locked Before Script 1 Runs
## OrganismCore — Document BRCA-S6a
## Date: 2026-03-05

---

## DOCUMENT METADATA

```
document_id:        BRCA-S6a
series:             BRCA Deep Dive — Invasive Lobular Carcinoma (ILC)
folder:             Cancer_Research/BRCA/DEEP_DIVE/ILC/
type:               BEFORE-DOCUMENT
                    All predictions stated before
                    any data is loaded.
                    This document cannot be modified
                    after Script 1 runs.
date:               2026-03-05
author:             Eric Robert Lawson
                    OrganismCore
protocol_version:   Workflow_Protocol.md v2.0
dataset_planned:    PRIMARY: TCGA-BRCA (histologically
                    confirmed ILC samples, n~180)
                    Normal reference: TCGA-BRCA adjacent
                    normal (n~113) — same reference used
                    across all BRCA subtypes
                    SECONDARY: GSE109169 (ILC-specific)
                    TERTIARY: METABRIC (if needed —
                    ILC histology annotation available)
status:             LOCKED — predictions cannot change
                    after this document is committed
precursor_documents:
  ATTRACTOR_GEOMETRY_AXIOMS.md (Doc 90)
  BRCA_Subtypes.md (BRCA_Subtype_Orientation)
  BRCA-S3a (HER2-enriched before-document)
  BRCA-S3b (HER2-enriched Script 1 reasoning artifact)
  BRCA-S3d (HER2-enriched Script 2 reasoning artifact)
  Luminal A complete series (BRCA-S1a through S1c)
  TNBC complete series (BRCA-S4a through S4c)
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
  The normal cell is the LOBULAR EPITHELIAL CELL —
  the secretory cell within the breast lobule (acini).
  This is a DIFFERENT cell of origin from all other
  BRCA subtypes analyzed so far.
  All prior subtypes (LumA, TNBC, HER2-enriched)
  arose from cells within the DUCTAL lineage.
  ILC arises from the LOBULAR lineage.

  The lobular epithelial cell's defining property:
    — It expresses E-cadherin (CDH1) as part of its
      normal structural identity.
    — CDH1 holds lobular cells together in the acinar
      architecture. Without CDH1, lobular cells lose
      cohesion and dissociate.
    — CDH1 is therefore not merely a marker of the
      normal cell — it is the STRUCTURAL REQUIREMENT
      for the normal cell to exist as a functional unit.

  The correct terminal destination:
    A mature lobular epithelial cell that:
      — Expresses CDH1 (E-cadherin) — architectural
      — Expresses ESR1 (ER) — hormonal identity
      — Expresses FOXA1 — ER pioneer TF
      — Expresses GATA3 — luminal differentiation TF
      — Expresses PGR (PR) — ER target
      — Expresses CK8/18 — luminal cytokeratins
    This is the lobular equivalent of the mature
    luminal epithelial cell.

QUESTION 2: Are identity TFs of the correct
  terminal state expressed?
  PARTIALLY — ILC is almost always ER+.
  ESR1, FOXA1, and GATA3 are RETAINED in ILC
  (unlike HER2-enriched and TNBC where they fall).
  PGR and CK8/18 are also retained.
  The TRANSCRIPTION FACTOR identity programme
  of the normal lobular cell is INTACT.

  HOWEVER:
  CDH1 ��� the STRUCTURAL component of the normal
  lobular cell identity — is LOST.
  This is the defining event of ILC.
  In ~65% of ILC cases: biallelic CDH1 mutation.
  In the remaining ~35%: promoter methylation or
  other CDH1 silencing mechanism.
  The result is the same: E-cadherin protein
  is absent from the cell surface.

QUESTION 3: Is the cell structurally within
  the correct valley but the floor is removed?
  THIS IS THE KEY DISTINCTION.

  ILC uniquely matches a modified Type 3 geometry:
    — The TF programme of the correct terminal state
      IS present (ESR1+, FOXA1+, GATA3+)
    — The cell IS in the correct valley transcriptionally
    — BUT: the structural component that makes the
      correct valley a STABLE ATTRACTOR (CDH1)
      has been removed.
  
  When CDH1 is lost, the lobular epithelial cell
  can no longer anchor itself in the normal acinar
  architecture. The Waddington valley is present
  but the WALLS have been dissolved — not the floor.
  The cell does not escape into a different identity
  (it stays ER+, FOXA1+, GATA3+) but it LOSES
  COHESION and acquires invasive behaviour through
  p120-catenin redistribution.

  This is structurally distinct from Type 3 as
  previously defined (Type 3 = correct valley,
  floor removed = differentiation programme
  partially maintained but proliferative scaffold
  activated beneath it).

  ILC represents a NEW VARIANT:
  TYPE 3 VARIANT — "ADHESION LOCK DISSOLUTION"
  The cell is in the correct transcriptional valley
  but the structural lock (CDH1) that maintains
  architectural identity has been removed,
  enabling dissociation and single-file invasion
  WITHOUT transcription factor reprogramming.

QUESTION 4: Is the cell expressing identity TFs
  of a DIFFERENT normal cell type?
  NO — ILC does not adopt basal, mesenchymal,
  HER2-amplicon, or stem cell identity.
  KRT5/KRT14, SOX10, ERBB2, CD44(high)/CD24(low)
  are NOT characteristic of ILC.
  The cell identity remains LUMINAL throughout.
  The disease mechanism is ARCHITECTURAL, not
  transcriptional.

ATTRACTOR TYPE DIAGNOSTIC CONCLUSION:

  ILC IS STRUCTURALLY UNIQUE IN THIS FRAMEWORK.

  Primary assignment: TYPE 3 VARIANT —
  "ADHESION LOCK DISSOLUTION"

  The false attractor in ILC is not defined by
  transcription factor replacement (Type 2),
  blocked approach to a terminal state (Type 1),
  or proliferative override within the correct
  valley (classic Type 3).

  ILC is defined by the REMOVAL OF THE STRUCTURAL
  LOCK that holds the correct terminal state intact.
  The cell retains its transcriptional identity
  but loses its physical identity.

  This is the only subtype in this repository where:
    — The correct identity TF programme is retained
    — The disease is driven by STRUCTURAL rather
      than transcriptional reprogramming
    — The false attractor is a physical/architectural
      state, not a molecular identity switch

  PREDICTED GEOMETRIC SIGNATURE (pre-data):
    — ESR1: HIGH (retained — not a marker of depth)
    — FOXA1: HIGH (retained)
    — GATA3: HIGH (retained)
    — PGR: HIGH (retained)
    — CK8/18: retained
    — CDH1: LOW or ABSENT (the depth marker)
    — p120-catenin (CTNND1): redistributed
      (cytoplasmic rather than membrane)
    — VIM: LOW (not a mesenchymal transition)
    — KRT5/KRT14: LOW (not a basal transition)
    — ERBB2: LOW-MODERATE (not amplified in ILC)
    — TP53: very LOW mutation rate (~4%)
      The most genomically stable subtype.
```

---

## PART II — CROSS-SUBTYPE CONTEXT
## What prior analyses established

```
Three BRCA deep dives are now complete:
  BRCA-S1 (Luminal A):
    — Type 3 (floor removed below a ER+ identity)
    — Depth axis: ESR1/FOXA1/GATA3 fall as EZH2 rises
    — The luminal identity TFs are the depth markers
    — EZH2 is the epigenetic driver of attractor depth
    — CDK4/6 inhibitor sensitivity confirmed

  BRCA-S4 (TNBC / Basal-like):
    — Type 2 (wrong valley — basal false attractor)
    — Depth axis: ESR1/FOXA1 erasure + VIM/ZEB rise
    — BRCA1 loss is the gate to the wrong valley
    — Deep cells are more mesenchymal, more resistant
    — pCR prediction: deep = resistant (confirmed)

  BRCA-S3 (HER2-enriched):
    — Type 1 variant "Slope Arrest"
    — Depth axis inverted: CDH1/CDH3 high at top,
      decline with depth; ESR1 very low throughout
    — The HER2 amplicon arrests descent on the slope
    — EZH2 high but correlation with ESR1 not confirmed
    — FOXA1 retention distinguishes HER2 from TNBC
      (the most important single confirmed finding)
    — Deep end cells: AR-low, more basal-adjacent,
      potential TNBC crossover sub-population

  WHAT THIS MEANS FOR ILC:
    — ILC is the FOURTH GEOMETRY in the breast landscape
    — It shares the ER+ identity with LumA but has
      a completely different false attractor mechanism
    — The depth axis in ILC CANNOT be ESR1 loss —
      ESR1 is retained throughout ILC by definition
    — The depth axis must be CDH1 and its downstream
      structural consequences
    — ILC is the direct counterpoint to Type 2 (TNBC):
        TNBC = TF identity lost, structure retained
        ILC  = TF identity retained, structure lost
    — This is the breast landscape's structural
      inversion: the two extreme geometries are
      exact opposites on the identity/structure axis
```

---

## PART III — GENE PANEL AND RATIONALE
## What to measure and why

```
TIER 1 — THE DEFINING AXIS (CDH1 structural axis)
These genes define ILC vs. normal lobular identity.
All predictions derive from this axis.

  CDH1 (E-cadherin)
    RATIONALE: The founding event of ILC.
    Loss defines ILC. Must be near-absent in
    cancer vs. normal.
    EXPECTED: Massively downregulated in ILC.
    This is the anchor of the depth axis.

  CTNND1 (p120-catenin)
    RATIONALE: Downstream effector of CDH1 loss.
    When CDH1 is lost, p120-catenin redistributes
    from the membrane to the cytoplasm.
    At the RNA level, total CTNND1 expression may
    not fall dramatically, but its functional
    state changes completely.
    EXPECTED: expression may be maintained but
    functionally context-shifted. Monitor for
    any bulk expression change.

  CDH3 (P-cadherin)
    RATIONALE: From HER2 deep-dive — CDH3 was
    the major depth signal in HER2-enriched.
    CDH3 is expressed in myoepithelial cells
    and in some deep cancer cells.
    In ILC, CDH3 upregulation may signal a
    partial shift toward a less-differentiated
    or myoepithelial-adjacent state in deep cells.
    EXPECTED: LOW in ILC overall (luminal identity
    retained) but may RISE with attractor depth.

TIER 2 — LUMINAL IDENTITY MARKERS
(predicted to be RETAINED in ILC — contrast with
LumA, TNBC, and HER2-enriched where they fall)

  ESR1 (estrogen receptor alpha)
    EXPECTED: HIGH throughout ILC — this is what
    makes ILC luminal. ESR1 should NOT be the
    depth axis. Any fall in ESR1 with depth would
    indicate a Type 3 component and would be a
    novel finding.

  FOXA1 (pioneer TF for ER binding)
    EXPECTED: HIGH — retained in ILC.
    FOXA1 mutations occur in ~14% of ILC and
    alter ER binding patterns (may not affect
    bulk expression dramatically).
    Monitor for any expression change.

  GATA3 (luminal differentiation TF)
    EXPECTED: HIGH — luminal identity retained.
    Should not fall as a function of attractor depth.

  PGR (progesterone receptor)
    EXPECTED: HIGH — ILC is ER+PR+ in most cases.
    PGR is an ER target gene.
    Any fall would indicate partial ER axis
    dysfunction (would be a novel finding).

  CK8 (KRT8), CK18 (KRT18)
    EXPECTED: HIGH — luminal cytokeratins
    retained in ILC. Structural identity markers.

TIER 3 — WHAT SHOULD BE ABSENT
(markers that characterize other subtypes —
ILC should NOT express these prominently)

  KRT5, KRT14 (basal cytokeratins)
    EXPECTED: LOW — ILC is not basal-like.
    Any elevation would indicate a basal-adjacent
    false attractor component.

  VIM (vimentin — mesenchymal marker)
    EXPECTED: LOW-MODERATE.
    ILC is not a mesenchymal transition cancer.
    Single-file invasion occurs WITHOUT the
    canonical EMT programme.
    This is a critical distinction from claudin-low.

  ZEB1, ZEB2, SNAI1 (EMT transcription factors)
    EXPECTED: LOW — ILC does not use canonical
    EMT. If these are elevated, it would challenge
    the structural-not-transcriptional model.

  ERBB2 (HER2)
    EXPECTED: LOW-MODERATE (not amplified in ILC).
    Amplification rate is ~2% in ILC.
    Should be clearly lower than HER2-enriched.

  SOX10 (basal/neural crest TF)
    EXPECTED: LOW — ILC is not TNBC.

TIER 4 — EPIGENETIC AXIS
(What does EZH2 do in a subtype where ESR1 is retained?)

  EZH2 (PRC2 catalytic subunit)
    RATIONALE: EZH2 was the depth driver in
    Luminal A (ESR1 suppressor) and was high in
    HER2-enriched. In ILC, ESR1 is retained,
    so EZH2 cannot be targeting ESR1 silencing
    in the same way.
    HYPOTHESIS: EZH2 in ILC may be targeting
    CDH1 itself — EZH2-mediated H3K27me3 at the
    CDH1 promoter is a known silencing mechanism.
    EXPECTED: EZH2 elevated in ILC vs. normal,
    and correlation of EZH2 with CDH1 suppression
    rather than ESR1 suppression.
    This would represent a completely different
    EZH2 target axis from LumA.

  EED, SUZ12 (PRC2 complex members)
    EXPECTED: Co-elevated with EZH2 if PRC2
    complex is activated.

  DNMT3A (de novo DNA methyltransferase)
    RATIONALE: CDH1 promoter methylation accounts
    for ~35% of ILC cases (the non-CDH1-mutant cases).
    EXPECTED: Elevated in ILC vs. normal — DNMT3A
    may be the mechanism of CDH1 silencing in
    the methylation-driven ILC subset.

  KDM1A (LSD1 — demethylase)
    EXPECTED: Monitor — may play a role in
    maintaining the CDH1-silenced state.

TIER 5 — PROLIFERATION AND CELL CYCLE
(ILC is characteristically low-proliferation)

  MKI67 (Ki-67)
    EXPECTED: LOW — ILC is almost always Grade 1-2.
    This distinguishes it from all other subtypes
    analyzed so far (LumA Grade 1-2 but HER2 and
    TNBC are Grade 3).
    LOW Ki-67 in ILC cancer cells vs. high Ki-67
    in other subtype cancer cells would confirm
    the subtype-specific proliferative signature.

  TOP2A, PCNA, CCNB1
    EXPECTED: LOW — consistent with low-grade
    biology.

  CDK2, CDK4
    EXPECTED: Monitor — CDK4/6 inhibitors are
    used in metastatic ILC (same as LumA).

TIER 6 — PI3K SIGNALLING
(PIK3CA is the most common driver in ILC at ~48%)

  PIK3CA (not directly measurable at RNA level
    but downstream targets are)
  AKT1, MTOR, PTEN
    EXPECTED: AKT1/MTOR axis elevated (PIK3CA
    mutant tumors constitutively activate this
    pathway). PTEN may be reduced.
    This is the main proliferative driver in ILC
    given the low TP53 rate and low Ki-67.

TIER 7 — ILC-SPECIFIC MARKERS

  TBX3
    RATIONALE: TBX3 mutations are enriched in
    ILC vs. IDC. TBX3 is a transcription factor
    involved in mammary gland development.
    EXPECTED: Expression change vs. normal?
    Monitor.

  RUNX1
    RATIONALE: RUNX1 mutations enriched in ILC.
    RUNX1 normally promotes luminal differentiation.
    Mutations may alter ER programme without
    eliminating ER expression.
    EXPECTED: Expression change vs. normal? Monitor.

  CDX2, SPI1, MBP (lineage controls)
    EXPECTED: ALL LOW — these are negative
    controls. Absent in breast tissue.
```

---

## PART IV — PREDICTIONS
## All stated 2026-03-05 before Script 1 runs

```
PREDICTION FRAMEWORK NOTE:
  ILC is structurally different from all prior
  subtypes because the depth axis is NOT defined
  by transcription factor loss.
  Predictions are organized around the CDH1
  structural axis rather than the ESR1/FOXA1
  identity axis used in LumA, TNBC, and HER2.
```

### P1 — CDH1 AS THE DOMINANT CANCER-VS-NORMAL SIGNAL

```
PREDICTION:
  CDH1 will be the most significantly
  downregulated gene in ILC cancer cells
  vs. normal lobular/luminal reference.

RATIONALE:
  CDH1 loss is the defining molecular event
  of ILC (~65% mutation, ~35% methylation).
  It is the equivalent of ERBB2 amplification
  in HER2-enriched — the single founding event
  that defines the entire subtype.
  Unlike the other founding events (ESR1 loss
  in LumA, BRCA1 loss in TNBC, ERBB2 amplification
  in HER2), CDH1 loss leaves the transcription
  factor identity intact while destroying the
  structural programme.

EXPECTED RESULT:
  CDH1: -60% to -90% in ILC vs. normal
  Effect size: large (comparable to ESR1 in TNBC)
  p-value: <0.001

FALSIFICATION:
  If CDH1 is NOT the top cancer-vs-normal signal,
  the structural-loss model is challenged.
  If CDH1 expression is only modestly reduced,
  the bulk RNA-seq may be capturing a mixed
  population (CDH1-mutant and CDH1-methylated
  cells together with residual normal cells).
  This is a known limitation of bulk data for ILC.
```

### P2 — LUMINAL IDENTITY TFs ARE RETAINED (NOT THE DEPTH AXIS)

```
PREDICTION:
  ESR1, FOXA1, GATA3, and PGR will be
  HIGH in ILC cancer vs. reference —
  comparable to Luminal A, NOT reduced as
  in TNBC or HER2-enriched.
  These genes will NOT correlate with
  attractor depth in ILC.

RATIONALE:
  ILC retains luminal identity by definition.
  This is the structural inversion argument:
  in TNBC, the TF programme is destroyed;
  in ILC, the TF programme is intact.
  This prediction establishes the baseline:
  ILC is not a luminal identity loss disease.

EXPECTED RESULT:
  ESR1: comparable to LumA cancer cells
  FOXA1: comparable to LumA cancer cells
  GATA3: comparable to LumA cancer cells
  PGR: comparable to LumA cancer cells
  None of these should show a depth gradient
  in ILC (unlike LumA where EZH2→ESR1 fall
  defines the depth axis).

FALSIFICATION:
  If ESR1 or FOXA1 show a clear depth gradient
  in ILC (as they do in LumA), then ILC has
  a DUAL mechanism: both CDH1 structural loss
  AND a partial Type 3 ESR1-silencing component.
  This would be a novel and important finding
  — ILC with deeper attractor depth may be
  acquiring a LumA-type epigenetic silencing
  component on top of CDH1 loss.
  This is specifically worth watching.
```

### P3 — EZH2 TARGETS CDH1 NOT ESR1 IN ILC

```
PREDICTION:
  EZH2 will be elevated in ILC vs. normal,
  BUT its expression will correlate with
  CDH1 suppression rather than ESR1 suppression.
  This is the opposite of Luminal A where
  EZH2 drives ESR1 silencing.

RATIONALE:
  EZH2-mediated H3K27me3 at the CDH1 promoter
  is a documented epigenetic silencing mechanism.
  In the ~35% of ILC cases driven by CDH1 promoter
  methylation (non-CDH1-mutant), EZH2 and DNMT3A
  may be the primary CDH1 silencing machinery.
  If EZH2 targets CDH1 in ILC while targeting
  ESR1 in LumA, the same epigenetic enzyme is
  driving two completely different attractor
  geometries in two ER+ subtypes.
  This would be a cross-subtype framework insight
  of significant importance.

EXPECTED RESULT:
  EZH2: elevated in ILC vs. normal (2-4 fold)
  EZH2 correlation with CDH1: negative
    (high EZH2 → low CDH1)
  EZH2 correlation with ESR1: minimal or absent
    (high EZH2 does NOT predict ESR1 loss in ILC)

FALSIFICATION:
  If EZH2 does NOT correlate with CDH1 suppression,
  the epigenetic silencing mechanism for CDH1 in
  ILC is not PRC2-mediated in bulk RNA (may be
  a DNA methylation dominant mechanism — DNMT3A
  would then be the primary signal to interrogate).
  If EZH2 DOES correlate with ESR1 in ILC
  (same as LumA), then EZH2 in ILC is behaving
  like LumA epigenetics despite the different
  attractor geometry.
```

### P4 — DNMT3A ELEVATED AS CDH1 METHYLATION DRIVER

```
PREDICTION:
  DNMT3A will be elevated in ILC vs. normal
  breast tissue, consistent with CDH1 promoter
  methylation as an alternative CDH1-silencing
  mechanism to CDH1 mutation.

RATIONALE:
  ~35% of ILC does not carry CDH1 mutation —
  these cases achieve CDH1 loss through promoter
  hypermethylation. DNMT3A is the de novo
  methyltransferase most commonly associated
  with cancer-specific promoter methylation.
  If DNMT3A is systematically elevated across
  ILC (not just in the methylation subset),
  it suggests the methylation machinery is
  broadly active in ILC, possibly targeting
  multiple loci beyond CDH1.

EXPECTED RESULT:
  DNMT3A: elevated in ILC vs. normal (2-3 fold)
  If expressible with current data:
  DNMT3A elevation may correlate with lower
  CDH1 in the CDH1-mutation-negative sub-cluster.

FALSIFICATION:
  If DNMT3A is NOT elevated in ILC bulk RNA,
  CDH1 promoter methylation in the non-mutant
  subset is driven by other DNMT family members
  (DNMT1 maintenance methylation; DNMT3B)
  rather than DNMT3A de novo methylation.
```

### P5 — ILC IS LOW PROLIFERATION: Ki-67 PROXY GENES LOW

```
PREDICTION:
  ILC cancer cells will show LOW expression
  of proliferation markers (MKI67, TOP2A,
  PCNA, CCNB1) relative to the TNBC and
  HER2-enriched cancer cells analyzed in
  the same dataset.
  ILC will be more similar to Luminal A
  in proliferative signature than to TNBC
  or HER2-enriched.

RATIONALE:
  ILC is almost always Grade 1-2.
  The clinical Ki-67 staining is low (<14%
  in most ILC cases).
  The attractor geometry here is not
  proliferative-drive-dominant (unlike HER2-enriched
  where the amplicon constitutively drives
  proliferation, or TNBC where MKI67 is
  near-universally high).
  The ILC false attractor is maintained by
  architectural disruption, not by
  proliferative override.

EXPECTED RESULT:
  MKI67: LOW in ILC cancer (comparable to LumA,
    clearly lower than TNBC or HER2 cancer cells
    if cross-subtype comparison is possible)
  TOP2A, PCNA, CCNB1: LOW
  This is a distinguishing feature that should
  be cleanly visible in the data.

FALSIFICATION:
  If MKI67 is HIGH in ILC cancer cells, the
  attractor geometry involves a proliferative
  component that the framework did not predict.
  High-grade ILC exists as a rare variant —
  if the dataset contains Grade 3 ILC cases
  they would confound this prediction.
  Grade annotation in TCGA-BRCA should be
  checked before drawing conclusions.
```

### P6 — PIK3CA PATHWAY ELEVATED AS THE PROLIFERATIVE DRIVER

```
PREDICTION:
  The PI3K/AKT/mTOR pathway will be the
  dominant proliferative signalling axis
  in ILC, elevated vs. normal, and will
  be the primary drug target signal.
  This is predicted because:
    — PIK3CA mutation is the most common
      driver mutation in ILC (~48%)
    — TP53 is very rare in ILC (~4%)
    — No HER2 amplification
    — The cell has no other obvious
      constitutive growth factor pathway
      beyond PIK3CA

EXPECTED RESULT:
  AKT1: elevated in ILC vs. normal
  MTOR: elevated or pathway-activated proxies
  PTEN: reduced (PTEN loss enables PIK3CA gain)
  EGFR: NOT elevated (distinguishes from TNBC)
  ERBB2: NOT elevated (distinguishes from HER2)

FALSIFICATION:
  If the PI3K axis is not clearly elevated in
  ILC, another proliferative mechanism is dominant
  — this would require re-examination of what
  drives ILC proliferation in the absence of
  the usual suspects (HER2, TP53, EGFR).
```

### P7 — ABSENCE OF BASAL OR MESENCHYMAL PROGRAMME

```
PREDICTION:
  ILC will NOT show elevation of basal cytokeratins
  (KRT5, KRT14) or mesenchymal markers (VIM, ZEB1,
  ZEB2, SNAI1, FOXC1, SOX10).
  ILC single-file invasion is NOT canonical EMT.
  The framework predicts that ILC achieves
  dissociation through CDH1 loss ALONE, without
  the transcriptional mesenchymal transition
  that characterises claudin-low or TNBC.

RATIONALE:
  This is the structural inversion principle:
  ILC breaks cohesion structurally (CDH1 loss),
  not transcriptionally (EMT programme activation).
  If VIM and ZEB1 are LOW while CDH1 is LOW,
  that is the geometric signature of structural-
  lock dissolution without valley switching.
  It is structurally unlike anything in the
  rest of this repository.

EXPECTED RESULT:
  VIM: LOW (unlike TNBC where VIM is high)
  ZEB1, ZEB2: LOW
  SNAI1: LOW
  KRT5, KRT14: LOW
  SOX10: LOW
  All clearly lower than TNBC cancer cells

FALSIFICATION:
  If VIM or ZEB1 are elevated in ILC
  (approaching TNBC levels), ILC has an EMT
  component and the structural-lock model
  needs revision — ILC may be a partial
  mesenchymal transition subtype.
  This would link ILC more closely to claudin-low
  than currently predicted.
```

### P8 — DEPTH SCORE: CDH1 AS THE AXIS, NOT ESR1

```
PREDICTION:
  The depth score in ILC will be defined by
  CDH1 suppression (and downstream structural
  consequences) rather than ESR1 suppression.
  Cells with lower CDH1 are DEEPER in the
  ILC false attractor.
  ESR1 will NOT define attractor depth in ILC.

RATIONALE:
  In LumA: depth = ESR1/FOXA1/GATA3 fall
  In TNBC: depth = ESR1/FOXA1 erasure + VIM/ZEB rise
  In HER2: depth = CDH1/CDH3 fall (inverted)
  In ILC: depth = CDH1 fall (structural axis)
  
  The cross-subtype pattern emerging is that
  CDH1 is the STRUCTURAL DEPTH MARKER:
    — In HER2-enriched: CDH1 fall correlates
      with deeper attractor (deep-end sub-population)
    — In ILC: CDH1 loss IS the attractor
  This is converging evidence that CDH1 is a
  cross-subtype structural attractor marker.

PREDICTED DEPTH AXIS SIGNATURE:
  Low CDH1 → deep ILC attractor
  Low CDH1 + low proliferation = less likely
    to respond to standard chemotherapy
  Low CDH1 + elevated EZH2/DNMT3A = epigenetically
    locked state (may respond to EZH2i or DNMT
    inhibition — novel prediction for ILC)

FALSIFICATION:
  If CDH1 does not show continuous variation
  within the ILC cancer population (i.e., CDH1
  is universally absent with no gradient), then
  CDH1 alone cannot serve as a depth axis.
  In that scenario, a downstream marker must
  define depth — p120-catenin redistribution,
  PIK3CA pathway activity, or CDH3 elevation.
```

### P9 — DRUG TARGET PREDICTION FROM GEOMETRY ALONE
### Stated before literature check

```
PRIMARY DRUG PREDICTION:
  ENDOCRINE THERAPY (standard — ESR1 retained)
  ER is retained → endocrine therapy is the
  foundational treatment.
  Tamoxifen and aromatase inhibitors work.
  CDK4/6 inhibitors (palbociclib, ribociclib,
  abemaciclib) should work by the same
  mechanism as in LumA (ER+ proliferative
  override through CDK4/6 axis).
  This is confirmed standard of care —
  not a novel prediction.

SECONDARY DRUG PREDICTION (novel from geometry):
  EZH2 INHIBITORS (tazemetostat class)
  If EZH2 is elevated and correlates with
  CDH1 suppression in ILC, then EZH2 inhibition
  may relieve CDH1 silencing and force
  re-expression of E-cadherin — partially
  restoring the structural lock.
  This is a NOVEL GEOMETRIC PREDICTION:
  EZH2i in ILC should be tested not for
  ESR1 re-expression (as would be predicted
  in LumA) but for CDH1 re-expression.
  The therapeutic goal is structural reversion
  rather than TF identity reversion.
  This prediction is falsifiable by ILC
  cell line EZH2i experiments measuring
  CDH1 and E-cadherin protein.

TERTIARY DRUG PREDICTION (novel from geometry):
  PI3K/AKT/mTOR PATHWAY INHIBITORS
  PIK3CA mutation at 48% is the dominant
  driver mutation. CDK4/6 + PI3K/AKT/mTOR
  dual blockade is the logical extension of
  current standard of care.
  Alpelisib (PI3Kα inhibitor, SOLAR-1 trial)
  is approved for PIK3CA-mutant ER+ MBC —
  already applicable to ILC, possibly with
  enhanced benefit given the higher PIK3CA
  rate in ILC vs. IDC-LumA.
  Everolimus (mTOR inhibitor) is used in
  endocrine-resistant ER+ MBC and ILC is
  a natural target population.

QUATERNARY DRUG PREDICTION (geometric):
  ANTI-CDH1 DOES NOT WORK (structural reason)
  CDH1 is already absent in ILC.
  Any therapy requiring CDH1 surface expression
  as a target will not function in ILC.
  This is a negative prediction from geometry —
  CDH1-based ADCs or immunotherapy would
  have no target in ILC cancer cells.
  However, CDH1 RE-EXPRESSION after EZH2i
  could make such therapies applicable to
  treated ILC — a combination strategy.

PREDICTION NOT MADE:
  Anti-HER2 therapy:
  ERBB2 is not amplified in ILC.
  Trastuzumab has no geometric basis here.
  This is confirmed literature.
```

---

## PART V — WHAT SCRIPT 1 DOES NOT PREDICT

```
Script 1 (single-cell RNA-seq, GSE176078)
will not directly address:

  — CDH1 protein / E-cadherin IHC
    (bulk RNA is the wrong measurement for
    structural protein distribution; CDH1
    mRNA may be present even when protein
    is absent due to post-translational
    regulation or mutation-only loss)

  — p120-catenin SUBCELLULAR LOCALIZATION
    (cytoplasmic vs. membrane redistribution
    cannot be determined from RNA-seq)

  — CDH1 PROMOTER METHYLATION STATUS
    (requires bisulfite sequencing or EPIC
    methylation array — not in GSE176078)

  — CDH1 MUTATION STATUS per cell
    (single-cell RNA-seq does not reliably
    detect mutations with normal coverage)

  — FOXA1 MUTATION FUNCTIONAL CONSEQUENCES
    (14% of ILC has FOXA1 mutations that
    alter ER binding patterns — not detectible
    from expression alone)

  — RESPONSE TO ENDOCRINE THERAPY
    (no treatment outcome data in scRNA dataset)

  — LATE RECURRENCE BIOLOGY
    (the ILC paradox — equal early prognosis,
    worse late recurrence — requires longitudinal
    data not available in scRNA)

WHAT SCRIPT 1 CAN ADDRESS:
  — CDH1 mRNA expression in ILC cancer cells
    vs. normal reference
  — Luminal TF retention (ESR1, FOXA1, GATA3)
  — Absence of basal/mesenchymal markers
  — EZH2 elevation
  — Proliferative marker levels
  — Depth score construction from CDH1 axis
  — Cross-comparison to LumA, TNBC, HER2 if
    using the same GSE176078 dataset
    (ILC cells may not be in this dataset —
    see NOTE below)

CRITICAL NOTE ON DATASET:
  GSE176078 (Wu et al. 2021) contains:
    Cancer Her2 SC cells — used for BRCA-S3
    Cancer LumA SC cells — used for BRCA-S1
    Cancer Basal SC cells — used for BRCA-S4
  IT IS UNCLEAR whether Wu et al. 2021
  includes histologically confirmed ILC cells
  or whether the "Cancer LumA SC" population
  includes ILC samples (which are often
  Luminal A by PAM50).

  IF GSE176078 does not contain ILC:
    Script 1 must use TCGA-BRCA bulk RNA-seq
    with histological ILC annotation.
    This changes the analysis from single-cell
    to bulk — lower resolution but larger n.
    The depth score would be a pseudo-bulk
    composite rather than a per-cell measure.

  Script 1 must resolve this data question
  FIRST before any gene expression analysis.
  If GSE176078 has ILC → use it.
  If not → use TCGA-BRCA histological ILC.
```

---

## PART VI — FALSIFICATION CRITERIA

```
The framework is falsified for ILC if ALL
of the following are true simultaneously:

  F1: CDH1 is NOT significantly reduced in
      ILC cancer cells vs. normal
      (would mean the founding event is not
      detectable at the RNA level — possible
      if mutation-only with no mRNA effect,
      but bulk RNA-seq in CDH1-methylated ILC
      should show reduction)

  F2: ESR1 IS significantly reduced in ILC
      cancer vs. normal
      (would mean ILC is a Type 3 subtype
      like LumA, not the structural-lock
      dissolution geometry predicted)

  F3: VIM and ZEB1 ARE elevated in ILC cancer
      (would mean ILC uses canonical EMT
      rather than structural CDH1 loss —
      the ILC geometry would be TYPE 2,
      not a novel variant)

  F4: EZH2 is NOT elevated in ILC
      (would mean the epigenetic lock mechanism
      is not PRC2-dependent in ILC — would
      require revision of the drug target logic)

  Any single falsified prediction generates
  a lesson. All four together falsify the
  structural-lock dissolution model and would
  require framework revision.

  PARTIAL FALSIFICATION:
  If CDH1 shows no gradient within the ILC
  cancer population (universally absent with
  no continuous variation), Prediction P8
  (depth axis from CDH1 gradient) is falsified.
  The depth axis would then need to be defined
  by a downstream marker (CDH3, p120-catenin
  pathway proxy, or PIK3CA target genes).
```

---

## PART VII — WHAT SCRIPT 1 NEEDS TO RUN

```
DATASET:
  Primary attempt: GSE176078 (Wu et al. 2021)
    Check for ILC cells in metadata
    If present: extract and use directly
    If absent: fall through to TCGA

  Fallback: TCGA-BRCA
    Expression matrix: TCGA_BRCA_HiSeqV2
    Clinical file: BRCA clinical matrix
    Histology filter: histological_type ==
      "Infiltrating Lobular Carcinoma"
    Normal reference: TCGA-BRCA adjacent normal
      (same reference used in other BRCA subtypes)

GENE PANEL:
  All genes listed in PART III above.
  At minimum:
    CDH1, CTNND1, CDH3
    ESR1, FOXA1, GATA3, PGR
    KRT8, KRT18, KRT5, KRT14
    VIM, ZEB1, ZEB2, SNAI1
    EZH2, EED, SUZ12, DNMT3A
    MKI67, TOP2A, PCNA, CCNB1
    AKT1, MTOR, PTEN, PIK3CA
    ERBB2, EGFR, SOX10
    TBX3, RUNX1
    CDX2, SPI1, MBP (controls)

ANALYSIS STEPS:
  1. Dataset acquisition and metadata check
     (confirm ILC cells or TCGA ILC samples)
  2. Expression matrix construction
  3. Cancer vs. normal differential expression
     (ILC cancer vs. adjacent normal)
  4. Ranked gene list — identify top
     cancer-vs-normal signals
  5. Panel gene expression levels
     (confirm CDH1 loss, ESR1 retention,
     absence of basal/mesenchymal markers)
  6. Depth score construction:
     CDH1-based (low CDH1 = deep)
  7. Depth score correlation:
     — EZH2 correlation with CDH1
     — EZH2 correlation with ESR1
       (should be absent in ILC)
     — DNMT3A correlation with CDH1
  8. Cross-subtype comparison (if using
     GSE176078 with ILC cells present):
     ILC vs. LumA vs. TNBC vs. HER2
     on CDH1, ESR1, VIM, EZH2 axes
  9. PCA geometry:
     Does ILC occupy a distinct PCA space
     from other subtypes? Does it sit
     adjacent to LumA (shared TF identity)
     or is it structurally separated?
  10. Figure generation and log output
```

---

## STATUS BLOCK

```
document:           BRCA-S6a (ILC predictions.md)
folder:             Cancer_Research/BRCA/DEEP_DIVE/ILC/
status:             LOCKED — written before Script 1
date:               2026-03-05
author:             Eric Robert Lawson
                    OrganismCore

attractor_type_assigned:
  TYPE 3 VARIANT — "ADHESION LOCK DISSOLUTION"
  Structural-lock loss without TF reprogramming.
  Novel geometry not previously encountered
  in this repository.

predictions_locked:   9
  P1: CDH1 dominant cancer-vs-normal signal
  P2: ESR1/FOXA1/GATA3 retained (not depth axis)
  P3: EZH2 targets CDH1, not ESR1 in ILC
  P4: DNMT3A elevated as CDH1 methylation driver
  P5: Low proliferation (Ki-67 proxy genes low)
  P6: PIK3CA/AKT/mTOR elevated (dominant driver)
  P7: No basal or mesenchymal programme
  P8: Depth axis is CDH1, not ESR1
  P9: Drug targets — endocrine, EZH2i, PI3Ki

drug_targets_predicted:
  Primary (confirmed SOC):   Endocrine therapy,
                             CDK4/6 inhibitors
  Secondary (novel):         EZH2 inhibitors
                             (for CDH1 re-expression,
                             not ESR1 re-expression)
  Tertiary (novel/rational): PI3K/AKT/mTOR inhibitors
                             (alpelisib, everolimus)
  Negative prediction:       Anti-HER2 (no amplicon)
                             Anti-CDH1 targeting
                             (CDH1 already absent)

key_structural_distinction:
  ILC is the INVERSE of TNBC:
    TNBC: TF identity erased, structure retained
    ILC:  TF identity retained, structure erased
  This is the breast landscape's geometric
  inversion pair.

cross_subtype_insight:
  CDH1 loss is a cross-subtype attractor
  depth marker:
    HER2-enriched deep cells: CDH1 falling
    ILC: CDH1 absent (extreme CDH1-loss state)
  CDH1 may be the structural depth axis
  across the entire breast landscape.

next_document:      BRCA-S6b — Script 1 results
                    and reasoning artifact
                    (written after Script 1 runs)
```
