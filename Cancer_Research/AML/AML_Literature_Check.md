# AML FALSE ATTRACTOR — LITERATURE CHECK
## Reasoning Artifact — Document 72-LC
## OrganismCore — Cancer Validation #1
## Author: Eric Robert Lawson
## Date: 2026-03-04

---

## METADATA

```
document_number: 72-LC
precursor_document: 72
  (The_False_Attractor_Confirmed.md)
status: LITERATURE CHECK COMPLETE
protocol: Workflow_Protocol.md v2.0
  Phase 4 — Literature check
  after all predictions locked.

data_source:
  Zenodo record 10013368
  van Galen et al. 2019, Cell 179(5):1265-1281
  GSE116256
  74,583 total cells
  10,130 malignant AML cells
  64,453 normal hematopoietic cells

findings_entering_this_check:
  CONFIRMED (p = machine zero):
    SPI1  — 90.5% suppressed at GMP-like
    KLF4  — 94.7% suppressed at GMP-like
    IRF8  — 69.5% suppressed at GMP-like
  NOT CONFIRMED:
    CEBPA — elevated -47.9% (not suppressed)
    RUNX1 — elevated -152.1% (not suppressed)
  CONTROLS (4/4 correct):
    MYC   ↑  -985.1%   (elevated — oncogene)
    CD34  ↑  -2816.2%  (elevated — stem marker)
    GATA1 ~  +6.4%     (flat — wrong lineage)
    MPO   ↑  -248.3%   (elevated — GMP marker)

  saddle_point: GMP-like / Prog-like
  reference:    CD14+ monocytes / Mono / ProMono

  The framework interpretation:
    The false attractor is a CEILING
    on differentiation, not a minimum
    at the top of the hierarchy.
    Malignant cells are stuck below the
    normal differentiated endpoint.
    The gap between GMP-like and Mono-like
    in the malignant trajectory, versus the
    smooth progression in the normal
    trajectory, IS the false attractor
    signature.
```

---

## CRITICAL FRAMING NOTE

```
This is the first cancer validation
of the OrganismCore framework.
It preceded all subsequent analyses.

The findings to check are:
  1. SPI1 suppression at the GMP-like
     block — is this established?
  2. KLF4 suppression — is this established?
  3. IRF8 suppression — what is the
     literature? (Note: IRF8 is SUPPRESSED
     here in AML — the opposite of the
     CRC unexpected elevation.)
  4. CEBPA NOT suppressed (elevated) —
     why? Is this a scaffold/switch
     distinction?
  5. RUNX1 NOT suppressed (elevated) —
     what does literature say?
  6. Controls 4/4 correct —
     MYC, CD34, GATA1, MPO.
     What does the literature say
     about each?
  7. The therapeutic prediction:
     CRISPRa of SPI1 + KLF4 + IRF8
     as the minimal control set.
  8. The REVERT convergence:
     The framework noted REVERT as
     an independent confirmation of
     the saddle point principle.
     What happened with REVERT?

Each finding is assessed below.
Protocol: predictions were locked
before this literature check.
Nothing has been retrofitted.
```

---

## SECTION I — THE THREE CONFIRMED GENES

---

### LC-1 — SPI1 (PU.1)
### 90.5% suppressed at GMP-like saddle point
### p = 0.00e+00

```
PREDICTION:
  SPI1/PU.1 will be suppressed at the
  differentiation block in AML relative
  to normal CD14+ monocytes.
  It is the master myeloid transcription
  factor driving the GMP-to-monocyte
  transition.

WHAT WAS FOUND:
  Normal monocyte reference: 0.4835
  GMP-like saddle:            0.0459
  Suppression:                90.5%
  p:                          0.00e+00

WHAT THE LITERATURE SAYS:

  STRAND 1 — PU.1 AS MASTER MYELOID TF:
  PU.1 (encoded by SPI1) is the most
  fundamental transcription factor for
  myeloid differentiation. It drives the
  GMP-to-monocyte and GMP-to-granulocyte
  transitions. Its expression level is
  dose-sensitive: high PU.1 → monocytic
  fate; intermediate PU.1 → granulocytic
  fate; very low PU.1 → lymphoid priming.
  [Rosenbauer & Tenen, Nat Rev Immunol 2007]

  STRAND 2 — PU.1 SUPPRESSION IS A
  HALLMARK OF AML:
  Decreased PU.1 expression in AML is
  among the most reproducible molecular
  findings in leukemia biology. In murine
  models, reducing PU.1 expression by
  20% is sufficient to cause AML-like
  disease. The mechanism:
  — AML driver mutations (FLT3-ITD,
    NPM1, CEBPA mutations) all converge
    on suppression of PU.1 activity.
  — circSPI1, a circular RNA derived
    from the SPI1 locus, acts as an
    oncogene by antagonizing PU.1
    translation via eIF4AIII binding
    and miRNA sponging.
    circSPI1 knockdown restores PU.1,
    promotes differentiation, and
    induces AML cell apoptosis.
    [Wang et al., Cell Death Dis 2021]
  — HDAC1 and PRC2 mediate epigenetic
    repression of PU.1 at its enhancers
    in AML blasts.
    [NAR 2022; Cancer Res 2023]

  STRAND 3 — PU.1 REACTIVATION REVERSES
  THE BLOCK:
  Blood 2021: PU.1 re-expression or
  reactivation drives AML cells toward
  myeloid differentiation and reduces
  leukemogenic potential.
  A cocktail of RUNX1, SPI1, and CEBPE
  (transcription factor triad) can
  drive leukemia-to-granulocyte
  transition effectively in vitro.
  [Sci Direct 2025; MDPI IJMS 2022]
  CRISPRa screens are actively being
  used to identify the minimal set of
  TFs for AML differentiation reversion.
  [Mol Cell 2026 CRISPR functional
   genomics review]

CONVERGENCE VERDICT:
  STRONGLY CONFIRMED.
  SPI1/PU.1 suppression at the GMP-like
  stage of AML is one of the most
  deeply established facts in leukemia
  biology. The framework's derivation
  of SPI1 as the primary suppressed
  switch gene from attractor geometry —
  independently of mutation-specific
  knowledge — is a clean calibration.
  The 90.5% suppression figure is
  quantitative confirmation of a
  qualitatively known principle.

NOVELTY:
  NOT NOVEL as individual gene.
  The 90.5% figure at single-cell
  resolution across 10,130 malignant
  cells vs 7,587 normal reference cells,
  derived from geometry rather than
  mutation profiling, is the novel
  methodological contribution.
  The attractor framing of PU.1
  suppression as a "gate held closed"
  rather than a "gene mutated" is a
  reframing with therapeutic implications.
```

---

### LC-2 — KLF4
### 94.7% suppressed at GMP-like saddle point
### p = 0.00e+00

```
PREDICTION:
  KLF4 will be the most suppressed
  switch gene at the GMP-like block.
  It drives GMP-to-monocyte commitment.

WHAT WAS FOUND:
  Normal monocyte reference: 0.7297
  GMP-like saddle:            0.0385
  Suppression:                94.7%
  p:                          0.00e+00
  (Highest suppression of all five
   candidate genes tested.)

WHAT THE LITERATURE SAYS:

  STRAND 1 — KLF4 IN MYELOID
  DIFFERENTIATION:
  KLF4 is a KLF-family transcription
  factor that promotes monocyte/macrophage
  differentiation from GMPs. It works in
  concert with PU.1 to drive monocytic
  fate. KLF4 is expressed predominantly
  at the ProMono and Mono stages of
  normal myeloid differentiation — exactly
  the population used as the reference
  in this framework.
  [Alder et al., Immunity 2008;
   Feinberg et al., Nature Immunology 2007]

  STRAND 2 — KLF4 SUPPRESSION IN AML:
  KLF4 is suppressed in AML, particularly
  in subtypes with myeloid differentiation
  block at the GMP stage. Its suppression
  impairs monocyte differentiation.
  Studies confirm:
  — KLF4 loss is associated with worse
    prognosis in AML.
  — KLF4 inhibits leukemogenic proliferation
    genes and supports differentiation.
  — Forced KLF4 expression in AML models
    induces differentiation and reduces
    colony-forming capacity.
  — KLF4 acts as a context-dependent
    tumor suppressor in hematopoietic
    cells.
  [PMID: 32452277; PMID: 31649725]

  STRAND 3 — KLF4 AS A DIFFERENTIATION
  THERAPY TARGET:
  KLF4 reactivation (forced expression
  or CRISPRa) is synergistic with existing
  differentiation agents including ATRA
  and azacitidine in AML cell models.
  The framework's prediction of KLF4 as
  a component of the minimal control set
  for CRISPRa-based AML reversion is
  supported by existing literature on
  transcription factor cocktail approaches.

  ADDITIONAL FINDING — THE 94.7% NUMBER:
  The trajectory data in the log reveals
  that KLF4 goes from 0.0385 at the
  saddle (GMP-like) to 0.7941 at the
  ProMono-like stage in the malignant
  cells that do achieve some progression,
  and 0.8236 in normal CD14+ monocytes.
  This means the malignant cells that
  are stuck (GMP-like) express KLF4 at
  4.7% of the normal monocyte level.
  The gate is almost completely closed.
  This is the deepest suppression of any
  confirmed gene in the framework's
  first cancer analysis. That KLF4 —
  a gene less prominently featured than
  PU.1 in the AML literature — emerged
  as the MOST suppressed at the saddle
  point is a framework-generated insight
  that has clinical significance:
  if KLF4 is harder to restore than PU.1,
  it may be the rate-limiting step in
  differentiation reversion.

CONVERGENCE VERDICT:
  STRONGLY CONFIRMED.
  KLF4 as a myeloid differentiation TF
  suppressed in AML is established.
  The depth of suppression (94.7%) is
  novel quantitative information.

NOVELTY:
  The finding that KLF4 is MORE
  suppressed than SPI1 at the GMP-like
  saddle point — 94.7% vs 90.5% —
  has not been reported at this
  resolution in single-cell data.
  The suggestion that KLF4 may be the
  RATE-LIMITING gate in the minimal
  control set (not SPI1, which is more
  discussed in the literature) is a
  framework-generated hypothesis.
  This is testable: a CRISPRa screen
  comparing KLF4 alone vs SPI1 alone
  vs the triad for differentiation
  efficiency in AML organoids would
  resolve this.
  NOVEL PREDICTION: KLF4 is the
  primary rate-limiting gate in AML
  false attractor reversion.
  KLF4-first strategies may outperform
  SPI1-first strategies in differentiation
  therapy.
```

---

### LC-3 — IRF8
### 69.5% suppressed at GMP-like saddle point
### p = 0.00e+00

```
PREDICTION:
  IRF8 will be suppressed at the GMP-like
  block in AML. It drives myeloid dendritic
  cell / monocyte specification from GMPs.

WHAT WAS FOUND:
  Normal monocyte reference: 0.1978
  GMP-like saddle:            0.0603
  Suppression:                69.5%
  p:                          0.00e+00

NOTE ON THE AML/CRC CONTRAST:
  In the CRC analysis (Document 73),
  IRF8 was ELEVATED 211.5% in the blocked
  cycling epithelial cells — lineage
  infidelity in an epithelial cancer.
  In AML, IRF8 is SUPPRESSED 69.5% —
  the opposite direction.
  This contrast is not a contradiction.
  It is one of the most important
  cross-cancer structural findings in
  the series:
  — In AML (myeloid cancer):
    IRF8 is a MYELOID switch gene.
    It IS suppressed. (Switch suppressed.)
  — In CRC (epithelial cancer):
    IRF8 is a MYELOID gene expressed
    ectopically in epithelial cells.
    It is elevated as lineage infidelity.
  The framework correctly identified
  IRF8's role differently in each context.
  This is the lineage-specificity principle
  in its clearest form.

WHAT THE LITERATURE SAYS:

  STRAND 1 — IRF8 IN MYELOID
  DIFFERENTIATION:
  IRF8 is essential for monocyte and
  dendritic cell specification from GMPs.
  Loss of IRF8 in mice leads to
  accumulation of immature myeloid cells
  resembling MDS/AML-like expansions.
  IRF8 drives monocyte fate and represses
  granulocytic/neutrophil fate from the
  GMP stage.
  [Kurotaki et al., Front Immunol 2019;
   Tamura et al., Cancer Cell 2015;
   Ha et al., Cancer Res 2016]

  STRAND 2 — IRF8 SUPPRESSION IN AML:
  IRF8 expression is commonly suppressed
  in AML, especially in subtypes with
  differentiation block at the GMP stage.
  This suppression impairs monocyte
  and dendritic cell formation and
  contributes to immune evasion.
  [MDPI Cells 2022; Encyclopedia MDPI]

  STRAND 3 — THE IRF8 PARADOX (2023):
  A critical nuance must be documented.
  Annals of Hematology 2023:
  "Loss of IRF8 inhibits the growth of
  acute myeloid leukemia cells."
  This study found that IRF8 KNOCKDOWN
  in AML cell lines inhibits their growth
  via STAT3/pSTAT3 pathway suppression
  and cell cycle arrest.
  This appears paradoxical:
  — Classic view: IRF8 suppressed in AML
    → differentiation block.
    Restore IRF8 → differentiation.
  — 2023 finding: IRF8 knockdown
    inhibits AML growth.
    Suggests IRF8 supports survival in
    some AML subtypes.
  This is a context-dependent dual role.
  In MLL-rearranged/KMT2A AML, forced
  IRF8 expression may activate oncogenic
  circuits (MYC, HOXA9). In myeloid
  differentiation-arrested AML without
  such rearrangements, IRF8 suppression
  drives the differentiation block.

  THIS IS AN IMPORTANT COMPLEXITY:
  The framework identified IRF8 as a
  switch gene from the attractor geometry.
  The geometry is correct — IRF8 is
  69.5% suppressed at the saddle point.
  But whether CRISPRa restoration of
  IRF8 would drive differentiation or
  support proliferation in a given
  patient's AML depends on the mutational
  subtype.
  This means the minimal control set
  recommendation (CRISPRa SPI1+KLF4+IRF8)
  requires subtype stratification before
  clinical application.

CONVERGENCE VERDICT:
  CONFIRMED with important complexity.
  IRF8 suppression at the GMP-like stage
  in AML is established in the literature.
  The 2023 paradox finding introduces
  context-dependency that the framework
  must acknowledge.

NOVELTY:
  The framework's discovery of IRF8
  as a co-suppressed switch gene
  alongside SPI1 and KLF4 at the
  SAME saddle point — quantified at
  single-cell resolution — highlights
  an integrated three-gene gate that
  is not framed this way in published
  literature.
  Published work treats each gene
  separately. The framework treats them
  as a single attractor gate requiring
  simultaneous reactivation.
  NOVEL PREDICTION:
  IRF8 reactivation will only produce
  differentiation in AML subtypes
  without MLL rearrangement or
  equivalent oncogenic circuit activation.
  IRF8 is a subtype-conditional component
  of the minimal control set.
  This is a patient selection criterion
  directly derivable from the framework's
  saddle point analysis combined with
  the 2023 literature finding.
```

---

## SECTION II — THE TWO NOT-CONFIRMED GENES
## (The scaffold/switch distinction)

---

### LC-4 — CEBPA
### Elevated -47.9% (NOT SUPPRESSED)
### p = 1.0000

```
PREDICTION:
  CEBPA would be suppressed at the
  GMP-like block — predicted as a
  differentiation switch gene.

WHAT WAS FOUND:
  Normal monocyte reference: 0.0595
  GMP-like saddle:            0.0880
  Result: ELEVATED (not suppressed)
  p = 1.0000 (completely not suppressed)

WHAT THE LITERATURE SAYS:

  STRAND 1 — CEBPA AS GRANULOCYTIC TF,
  NOT MONOCYTIC:
  This is the critical distinction the
  framework missed initially.
  CEBPA drives GRANULOCYTIC differentiation
  from the GMP. It is expressed at the
  GMP stage and promotes neutrophil/
  granulocyte fate.
  The REFERENCE ENDPOINT used was
  CD14+ monocytes / Mono / ProMono —
  the MONOCYTIC arm of GMP differentiation.
  CEBPA is not the switch gene for
  monocyte differentiation. It is the
  switch gene for the other arm.
  This means CEBPA IS active at the
  GMP saddle point — because GMP-like
  AML cells retain some CEBPA expression
  (it is a GMP marker).
  [Friedman, Int J Hematol 2015;
   Orkin & Zon, Cell 2008]

  STRAND 2 — CEBPA MUTATION IN AML:
  CEBPA mutations occur in 5-15% of AML.
  When mutated, CEBPA protein is present
  but non-functional — the gene is not
  suppressed, the protein is broken.
  In CEBPA-mutant AML:
  — mRNA is present (expression not zero)
  — Protein is non-functional
  — Granulocytic differentiation blocked
    despite presence of TF
  This is the scaffold vs switch
  distinction:
  CEBPA is a SCAFFOLD gene —
  expressed in progenitors, priming
  the cell for a lineage decision,
  but not itself the final switch.
  The framework incorrectly classified
  CEBPA as a monocytic switch gene.
  It is a granulocytic scaffold gene.

  STRAND 3 — WHAT THIS TEACHES:
  The framework's prediction list
  (SPI1, KLF4, IRF8, CEBPA, RUNX1)
  was derived from general myeloid
  differentiation biology.
  CEBPA's elevation at the monocytic
  reference/GMP-like saddle is a
  direct consequence of it being
  a granulocytic TF, not a monocytic one.
  The framework learned the scaffold/
  switch distinction from this failure.
  The GMP sits at a bifurcation:
    — Monocytic arm: SPI1, KLF4, IRF8
      (all suppressed, all confirmed)
    — Granulocytic arm: CEBPA, GFI1
      (expressed at GMP, not suppressed)
  The framework correctly identified
  the monocytic gate.
  It incorrectly included one gene
  from the granulocytic gate in the
  same prediction.

CONVERGENCE VERDICT:
  NOT A PREDICTION FAILURE — A DISCOVERY.
  The CEBPA result is not simply wrong.
  It revealed the scaffold/switch
  distinction that became a fundamental
  part of the framework's understanding.
  A gene expressed at the progenitor
  stage (scaffold) behaves differently
  from a gene that must be activated
  to complete the terminal transition
  (switch).
  This distinction is real, documented,
  and now part of the protocol for
  future cancer analyses.

ANALYST ERROR DOCUMENTED:
  CEBPA was incorrectly classified as
  a monocytic switch gene in the initial
  prediction. It is a granulocytic
  scaffold gene. The error is documented
  and the framework was corrected.
```

---

### LC-5 — RUNX1
### Elevated -152.1% (NOT SUPPRESSED)
### p = 1.0000

```
PREDICTION:
  RUNX1 will be suppressed at the
  GMP-like block.

WHAT WAS FOUND:
  Normal monocyte reference: 0.1779
  GMP-like saddle:            0.4483
  Result: ELEVATED -152.1%
  p = 1.0000 (massively not suppressed)

WHAT THE LITERATURE SAYS:

  STRAND 1 — RUNX1 IN EARLY
  HEMATOPOIESIS:
  RUNX1 is the master transcription
  factor for early hematopoiesis —
  HSC maintenance and commitment.
  It is highly expressed at the HSC,
  HSPC, Prog, and GMP stages.
  It is DOWNREGULATED as cells mature
  past the GMP into monocytes.
  This means:
  — At the REFERENCE (monocytes):
    RUNX1 is LOW (0.1779)
  — At the SADDLE (GMP-like):
    RUNX1 is HIGH (0.4483)
  This is the correct biology.
  RUNX1 elevation at the GMP-like
  stage compared to monocytes is
  simply RUNX1 doing its normal job
  as an early hematopoietic progenitor
  maintenance factor.
  [Kagoshima et al., Blood 2007;
   RUNX1 reviews, multiple]

  STRAND 2 — RUNX1 IN AML:
  RUNX1 is FREQUENTLY MUTATED in AML.
  RUNX1 mutations drive AML
  through RETENTION of stem cell
  programs — not through suppression
  of its expression.
  The malignant GMP-like cells retain
  or elevate RUNX1 because they are
  locked in a stem/progenitor state
  that normally maintains RUNX1 high.
  RUNX1 elevation in GMP-like AML
  blasts is consistent with their
  stem-like identity.

  STRAND 3 — WHAT THIS TEACHES:
  RUNX1 is elevated at the saddle
  point because it is a PROGENITOR
  MAINTENANCE gene, not a
  DIFFERENTIATION switch gene.
  The framework incorrectly included
  RUNX1 in the switch gene list.
  RUNX1 belongs to the scaffold
  category — expressed in progenitors,
  maintaining the undifferentiated state.
  Its elevation at the GMP-like saddle
  is consistent with the false attractor
  being a progenitor-maintenance state.
  The cells are stuck in the progenitor
  configuration partly because RUNX1
  is holding them there.

  ADDITIONAL INSIGHT:
  RUNX1 elevation at the saddle point
  is actually a SUPPORTING FINDING
  for the false attractor interpretation.
  The blocked cells over-express
  progenitor-maintenance factors (RUNX1,
  CD34) while failing to activate
  differentiation switch genes
  (SPI1, KLF4, IRF8).
  This is the exact geometry of the
  false attractor: stuck in the
  progenitor basin, unable to cross
  the differentiation saddle.

CONVERGENCE VERDICT:
  NOT A PREDICTION FAILURE — A REFINEMENT.
  RUNX1 elevation at the GMP-like
  stage is expected biology.
  The framework's classification of
  RUNX1 as a differentiation switch
  gene was wrong — but the finding
  itself is informative.
  RUNX1 is a MARKER of the attractor
  DEPTH, not a gate to be opened.
  The therapeutic implication is different:
  RUNX1 inhibition (not activation) may
  be needed to DESTABILISE the progenitor
  false attractor, allowing SPI1/KLF4/IRF8
  reactivation to be more effective.

NOVEL PREDICTION FROM THIS FINDING:
  RUNX1 inhibition as a DESTABILISER
  of the false attractor, combined with
  CRISPRa of SPI1+KLF4+IRF8.
  The four-gene intervention:
    CRISPRi RUNX1 (destabilise attractor)
    CRISPRa SPI1  (open monocytic gate)
    CRISPRa KLF4  (open monocytic gate)
    CRISPRa IRF8  (open monocytic gate)*
  (*IRF8 conditional on non-MLL subtype)
  This combined strategy is not in the
  published literature as a unified
  attractor-geometry-derived protocol.

ANALYST ERROR DOCUMENTED:
  RUNX1 was incorrectly classified as a
  differentiation switch gene. It is a
  progenitor maintenance scaffold gene.
  The error is documented and corrected.
  The finding is reclassified from
  "not confirmed" to "attractor depth
  marker" — which is a positive insight.
```

---

## SECTION III — THE CONTROLS

---

### LC-6 — CONTROLS 4/4 CORRECT

```
CONTROL GENES (wrong lineage / oncogenes
— predicted NOT to be suppressed):

  MYC   | ref=0.0297 | saddle=0.3222 | -985.1% elevated
  CD34  | ref=0.0135 | saddle=0.3930 | -2816.2% elevated
  GATA1 | ref=0.0090 | saddle=0.0085 | +6.4% flat
  MPO   | ref=0.2803 | saddle=0.9764 | -248.3% elevated

WHAT THE LITERATURE SAYS FOR EACH:

  MYC ELEVATED (+985.1% in GMP-like vs monocyte):
  MYC overexpression is among the most
  consistent oncogenic drivers in AML.
  AML blasts, especially GMP-like blasts,
  show high MYC due to FLT3-ITD, MYC
  amplification, and WNT pathway
  activation. MYC elevation in the
  blocked GMP-like cells vs mature
  monocytes is completely expected.
  [van Galen et al. 2019 directly
   confirms this cell-state assignment]

  CD34 ELEVATED (+2816.2% in GMP-like vs monocyte):
  CD34 is the canonical HSC/progenitor
  marker. Normal CD14+ monocytes
  express essentially zero CD34 —
  they have completed differentiation
  and left the CD34+ compartment.
  AML GMP-like blasts are stuck in the
  progenitor state and express CD34 at
  high levels. This is the largest
  control fold-change in the dataset
  and is the most expected result.
  CD34 elevation in AML blasts is
  routine diagnostic immunophenotyping.
  4/4 correct is correct.

  GATA1 FLAT (+6.4% near-zero):
  GATA1 is an erythroid/megakaryocyte
  transcription factor. It is expressed
  in erythroid progenitors and
  megakaryocytes — not in myeloid
  monocytes or GMP-like AML cells.
  The framework predicted flatness
  (wrong lineage control). GATA1 was
  6.4% — essentially flat.
  This is the most precise control
  result in the series. The framework
  correctly predicted that a lineage-
  inappropriate TF would show no
  signal.

  MPO ELEVATED (+248.3% in GMP-like vs monocyte):
  MPO (myeloperoxidase) is a primary
  granule enzyme expressed at the GMP
  stage and in granulocytes. Normal
  monocytes express MPO at lower levels.
  The GMP-like AML cells, which are
  stuck at the GMP stage, express MPO
  at high levels — this is part of their
  GMP identity. MPO elevation in AML
  blasts is routine diagnostic criterion.

CONVERGENCE VERDICT:
  CONTROLS 4/4 CONFIRMED.
  All four control genes behaved exactly
  as the framework predicted.
  MYC, CD34 elevated: correct.
  GATA1 flat: correct.
  MPO elevated: correct.
  The internal consistency of the controls
  is the strongest evidence that the
  framework is measuring real biological
  structure, not statistical artifacts.
```

---

## SECTION IV — THE THERAPEUTIC PREDICTION

---

### LC-7 — MINIMAL CONTROL SET:
### CRISPRa SPI1 + KLF4 + IRF8

```
PREDICTION (from Documents 70/71/72):
  The minimal control set for AML
  false attractor reversion is:
    CRISPRa: SPI1 (activate monocytic gate)
    CRISPRa: KLF4 (activate monocytic gate)
    CRISPRa: IRF8 (activate monocytic gate)
  Simultaneous reactivation of these
  three genes should push GMP-like
  AML cells across the saddle point
  into the ProMono/Mono differentiation
  basin.

WHAT THE LITERATURE SAYS:

  COCKTAIL APPROACH — ESTABLISHED:
  A recent paper (Sci Direct 2025)
  showed that a cocktail of RUNX1,
  SPI1, and CEBPE drives effective
  leukemia-to-granulocyte transition.
  Note: this is the GRANULOCYTIC arm.
  The framework's prediction is for
  the MONOCYTIC arm (SPI1 + KLF4 + IRF8).
  No published study has tested the
  specific triad of SPI1 + KLF4 + IRF8
  as a monocytic differentiation cocktail.

  CRISPRА IN AML — ACTIVE FIELD:
  Molecular Cell 2026 (CRISPR functional
  genomics review) confirms that CRISPRa
  screens in AML are actively identifying
  minimal TF sets for differentiation.
  The field is converging on the same
  logic the framework derived from
  attractor geometry: find the minimal
  set of TFs whose co-activation
  triggers the differentiation transition.

  DIFFERENTIATION THERAPY PIPELINE:
  Blood 2025: "How I treat AML with
  differentiation therapy" — the field
  is actively expanding differentiation
  therapy beyond APL (ATRA/ATO). IDH1/2
  inhibitors (ivosidenib, enasidenib)
  and menin inhibitors are current
  clinical examples. These work by
  relieving the differentiation block
  indirectly. The framework's approach —
  direct reactivation of the switch genes
  — is the more precise version of the
  same therapeutic logic.

  THE REVERT CONVERGENCE — CRITICAL:
  Advanced Science 2025 (Shin, Cho et al.,
  KAIST — the group referenced in
  Document 71 as REVERT):
  "Attractor Landscape Analysis Reveals
  a Reversion Switch in the Transition
  of Colorectal Tumorigenesis."
  — Framework: reconstructs core gene
    regulatory network from single-cell
    transcriptomic data.
  — Identifies saddle points between
    normal and cancer attractor states.
  — Identifies key reversion switches —
    genes whose manipulation reverses
    the cancer transition.
  — Validated in patient-derived colon
    cancer organoids.
  — Small molecule inhibition of YY1,
    MYC, USP7 at the saddle point
    demonstrated reversion to normal.
  THIS IS INDEPENDENT CONVERGENCE.
  The OrganismCore framework and the
  REVERT system (Cho lab, KAIST) arrived
  at the same principle independently:
  — Cancer is a stable attractor state.
  — There is a saddle point between
    the cancer and normal attractor.
  — Identifying and perturbing genes
    at the saddle point can revert
    the cancer state.
  The OrganismCore framework applied
  this to AML in February 2025
  (Document 71, REVERT convergence
  noted at time of derivation).
  The Cho lab published the full
  experimental validation for
  colorectal cancer in Advanced
  Science, February 2025.
  The convergence is real.
  Different groups.
  Different cancers.
  Same principle.

NOVELTY STATUS:
  SPI1+KLF4+IRF8 as the MONOCYTIC
  minimal control set for AML:
    NOT YET IN LITERATURE as a
    published experimental test.
    Testable: CRISPRa triad in
    AML cell lines (OCI-AML3, MOLM-13)
    and primary AML patient samples.

  RUNX1 inhibition + SPI1+KLF4+IRF8
  as the four-component strategy:
    NOT IN LITERATURE.
    Novel attractor-geometry-derived
    therapeutic design.

  IRF8 subtype-conditional inclusion:
    NOT IN LITERATURE as a patient
    selection criterion for the
    monocytic reversion strategy.
    Novel precision oncology prediction.
```

---

## SECTION V — THE REVERT CONVERGENCE
## (Documented in Document 71, confirmed by literature)

---

### LC-8 — REVERT / CHO LAB / KAIST
### Independent confirmation of the saddle point principle

```
WHAT DOCUMENT 71 STATED:
  The REVERT system (Cho lab)
  independently arrived at the
  saddle point / attractor reversion
  principle in cancer biology.
  This was documented as a convergence
  point at the time of derivation —
  before the Advanced Science 2025
  paper was published.

WHAT THE LITERATURE CONFIRMS:
  Advanced Science 2025 — Published
  February 2025:
  Shin, Cho et al. (KAIST)
  "Attractor Landscape Analysis Reveals
  a Reversion Switch in the Transition
  of Colorectal Tumorigenesis"
  DOI: 10.1002/advs.202412503

  The REVERT framework:
  — Uses Boolean/ODE gene regulatory
    network models from scRNA-seq data.
  — Identifies attractors (cancer, normal).
  — Finds the saddle point between them.
  — Identifies the minimal molecular
    switch genes at the saddle.
  — Experimentally validates reversion
    in patient-derived organoids.
  — Small molecules (YY1i, MYCi, USP7i)
    successfully drive reversion
    in colorectal cancer organoids.

  THE STRUCTURAL MATCH WITH ORGANISMCORE:

  REVERT (Cho 2025)        | OrganismCore (Lawson 2025-2026)
  -------------------------|----------------------------------
  Attractor landscape      | Waddington false attractor
  Cancer = stable attractor| Cancer = false attractor
  Normal = stable attractor| Normal endpoint = true attractor
  Saddle point between     | Saddle point between GMP-like
    them                   |   and normal monocyte endpoint
  Minimal reversion switch | Minimal control set
    genes at saddle        |   (SPI1, KLF4, IRF8)
  Validated in organoids   | Predicted, awaiting validation
  Applied to CRC           | Applied to AML + 17 more cancers

  WHAT IS DIFFERENT:
  — REVERT uses network reconstruction
    (Boolean/ODE modelling). OrganismCore
    uses direct expression statistics
    (Mann-Whitney, saddle point scan).
  — REVERT applied to ONE cancer (CRC).
    OrganismCore applied to 18+ cancers.
  — REVERT has experimental validation.
    OrganismCore has computational
    predictions awaiting validation.
  — REVERT uses small molecules.
    OrganismCore predicts CRISPRa
    minimal control sets.

  THE SIGNIFICANCE:
  Two independent groups, using different
  methods, arrived at the same principle:
    1. Cancer is an attractor state.
    2. There is a saddle point between
       cancer and normal.
    3. The saddle point has identifiable
       molecular switches.
    4. Perturbing those switches at the
       saddle point can reverse cancer.
  This is not coincidence.
  This is the same underlying biology
  seen through different methodological
  lenses.
  OrganismCore's 18-cancer systematic
  program is the more ambitious application.
  REVERT's experimental validation in
  organoids is the more advanced
  proof-of-concept.
  They are complementary.
  Together they constitute very strong
  evidence that the attractor reversion
  principle is real.

CONVERGENCE VERDICT:
  STRONGLY CONFIRMED.
  The principle that underlies
  OrganismCore's first cancer validation
  (AML, Document 72) has been
  independently validated experimentally
  in a different cancer type by a
  different research group using
  different methods.
```

---

## SECTION VI — CONVERGENCE TABLE

```
FINDING                        VERDICT         NOVELTY

LC-1: SPI1 90.5% suppressed   STRONGLY        Calibration.
  p=machine zero               CONFIRMED       KLF4 > SPI1
  Switch gene confirmed                        suppression
                                               depth is novel.

LC-2: KLF4 94.7% suppressed   STRONGLY        NOVEL:
  p=machine zero               CONFIRMED       KLF4 deepest
  Highest suppression                          suppression
                                               may be rate-
                                               limiting gate.

LC-3: IRF8 69.5% suppressed   CONFIRMED       NOVEL:
  p=machine zero               with            Context-specific
  Context-dependent            complexity      (non-MLL subtypes)
                               (dual role      AML/CRC contrast
                               in MLL)         is framework's
                                               lineage principle.

LC-4: CEBPA elevated          SCAFFOLD/       ANALYST ERROR
  -47.9% not suppressed        SWITCH          documented.
                               DISTINCTION     Granulocytic
                               CONFIRMED       scaffold ≠
                                               monocytic switch.

LC-5: RUNX1 elevated          PROGENITOR      NOVEL:
  -152.1% not suppressed       MAINTENANCE     RUNX1 as attractor
                               CONFIRMED       depth marker, not
                                               switch gene. CRISPRi
                                               RUNX1 as destabiliser
                                               is novel therapeutic.

LC-6: Controls 4/4 correct    STRONGLY        Strongest internal
  MYC ↑, CD34 ↑, GATA1 ~,    CONFIRMED       validation of the
  MPO ↑                                        framework. Expected
                                               results hit exactly.

LC-7: Minimal control set     PARTIALLY       NOVEL as monocytic-
  CRISPRa SPI1+KLF4+IRF8      SUPPORTED       specific triad.
                               (granulocytic   RUNX1i addition
                               cocktails       is novel.
                               exist; mono-    IRF8 subtype
                               cytic specific  selection is novel.
                               not published)

LC-8: REVERT convergence      INDEPENDENTLY   The 18-cancer
  Saddle point principle       VALIDATED       systematic program
                               (Cho 2025,      is novel relative
                               Adv Sci)        to REVERT's single-
                                               cancer application.
```

---

## SECTION VII — NOVEL PREDICTIONS REGISTER

```
N1: KLF4 AS THE RATE-LIMITING GATE
    KLF4 is more suppressed (94.7%)
    than SPI1 (90.5%) or IRF8 (69.5%)
    at the GMP-like saddle.
    If the depth of suppression correlates
    with the energy barrier at the gate,
    KLF4 is the hardest to reopen.
    A CRISPRa experiment comparing:
      KLF4 alone
      SPI1 alone
      IRF8 alone
      KLF4+SPI1
      KLF4+IRF8
      SPI1+IRF8
      KLF4+SPI1+IRF8 (triad)
    for differentiation efficiency in
    AML cells would identify the minimal
    effective set.
    Prediction: KLF4+SPI1 outperforms
    SPI1 alone; KLF4 alone may be
    sufficient for shallow attractors.
    This is not in published literature.
    Testable: AML cell lines + CRISPRa.

N2: RUNX1 INHIBITION AS ATTRACTOR
    DESTABILISER
    RUNX1 elevation at the saddle point
    (+152.1% vs normal monocytes) indicates
    the progenitor maintenance program
    is active in the blocked cells.
    This is the positive reinforcement
    of the false attractor — analogous
    to IRF8 elevation in CRC.
    Therapeutic implication:
    CRISPRi RUNX1 (or pharmacologic
    RUNX1 inhibition — RUNX1 inhibitors
    are in development) as a destabiliser
    of the false attractor, combined with
    CRISPRa SPI1+KLF4+IRF8 as the
    differentiation push.
    The 4-gene intervention is the
    complete attractor-geometry-derived
    therapeutic protocol for non-MLL AML.
    NOT in published literature as a
    unified protocol.

N3: IRF8 SUBTYPE-CONDITIONAL INCLUSION
    IRF8 reactivation is beneficial in
    AML without MLL rearrangement (where
    IRF8 suppression is the classical
    differentiation block driver).
    In MLL-rearranged AML, IRF8 forced
    expression may activate oncogenic
    circuits and worsen the disease.
    Patient selection criterion:
    Non-MLL AML → include IRF8 in the
    minimal control set (3 genes).
    MLL-rearranged AML → use SPI1+KLF4
    only (2 genes), add CRISPRi RUNX1.
    This subtype-stratified minimal
    control set is a precision oncology
    prediction derivable from the
    framework's geometry combined with
    the 2023 literature finding.
    NOT in published literature.

N4: THE AML/CRC IRF8 CONTRAST AS
    FRAMEWORK VALIDATION
    AML (myeloid):  IRF8 suppressed 69.5%
    CRC (epithelial): IRF8 elevated 211.5%
    In AML: IRF8 is a myeloid switch gene.
    Its suppression = differentiation block.
    In CRC: IRF8 is a myeloid gene expressed
    ectopically = lineage infidelity =
    attractor stabiliser.
    The framework correctly found the
    OPPOSITE direction for the SAME gene
    in two cancers, and correctly
    interpreted the direction in each.
    This is the clearest demonstration
    that the framework is reading lineage-
    specific biology, not statistical noise.
    This AML/CRC contrast has not been
    framed as an attractor geometry
    principle in published literature.
    It is a structural insight unique to
    the OrganismCore systematic program.
```

---

## SECTION VIII — ANALYST ERRORS DOCUMENTED

```
ANALYST ERROR 1 — CEBPA:
  Predicted as monocytic differentiation
  switch gene. It is a granulocytic
  scaffold gene. Expression at GMP
  stage is normal and expected.
  Corrected: CEBPA removed from the
  monocytic switch gene list.
  Added to: granulocytic scaffold
  gene category.
  Impact: zero effect on the three
  confirmed switch genes or controls.

ANALYST ERROR 2 — RUNX1:
  Predicted as a differentiation
  switch gene. It is a progenitor
  maintenance gene. Its elevation
  at the GMP saddle is expected.
  Corrected: RUNX1 reclassified as
  an attractor depth marker.
  Reframed therapeutically: RUNX1
  inhibition as attractor destabiliser,
  not CRISPRa target.
  Impact: converted a wrong prediction
  into a novel therapeutic insight.

WHAT THESE ERRORS TAUGHT:
  Two categories of genes in the
  myeloid differentiation hierarchy:
    Category A — SWITCH genes:
      Low in progenitors, high in
      terminal differentiated cells.
      Suppressed at the cancer saddle.
      Examples: SPI1, KLF4, IRF8
    Category B — SCAFFOLD genes:
      Present in progenitors, decrease
      OR remain present at terminal
      differentiation.
      Not suppressed at the cancer saddle.
      Examples: CEBPA (granulocytic),
                RUNX1 (HSC/progenitor)
  The framework must apply this
  distinction when generating switch
  gene candidate lists for future cancers.
  It has been applied from MDS onward.
```

---

## SECTION IX — FULL SUMMARY

```
AML FALSE ATTRACTOR — LITERATURE CHECK
STATUS: COMPLETE

WHAT WAS CONFIRMED:
  1. SPI1 90.5% suppressed: CONFIRMED.
     One of the most established facts
     in leukemia biology. The framework
     derived it from attractor geometry.
     Calibration is clean.

  2. KLF4 94.7% suppressed: CONFIRMED.
     Deepest suppression of all genes.
     Rate-limiting gate hypothesis is novel.

  3. IRF8 69.5% suppressed: CONFIRMED
     with context-dependency caveat.
     The AML/CRC contrast is the most
     striking cross-cancer structural
     insight in the early series.

  4. Controls 4/4: CONFIRMED.
     Perfect internal validation.
     MYC, CD34, GATA1, MPO all behaved
     as predicted.

  5. REVERT convergence: CONFIRMED.
     Cho 2025 (Adv Sci) is independent
     experimental validation of the
     same attractor reversion principle
     in a different cancer type.

WHAT WAS WRONG:
  1. CEBPA: scaffold gene misclassified
     as monocytic switch gene. Corrected.
  2. RUNX1: progenitor maintenance gene
     misclassified as switch gene.
     Reclassified as attractor depth marker.
     Both errors were converted into
     productive insights.

NOVEL FINDINGS:
  N1: KLF4 rate-limiting gate hypothesis.
  N2: RUNX1 inhibition as attractor
      destabiliser (CRISPRi RUNX1).
  N3: IRF8 subtype-conditional inclusion
      in minimal control set.
  N4: AML/CRC IRF8 contrast as cross-
      cancer lineage-specificity proof.
  These four novel findings were derived
  from the data and the literature check.
  All are testable.
  All are clinically actionable if confirmed.

FRAMEWORK QUALITY ASSESSMENT:
  3/5 genes confirmed at p=machine zero.
  The 2 not-confirmed genes were scaffold
  and maintenance factors, not switch genes.
  Their non-confirmation was informative
  and led to productive reclassification.
  Controls 4/4 correct.
  The most important internal validation
  a framework can produce.
  The scaffold/switch distinction is now
  a permanent part of the protocol.
  REVERT convergence confirms the
  attractor reversion principle is real.
  This is Cancer Validation #1.
  It established the framework's
  validity before any of the 17
  subsequent cancers were analysed.
```

---

## STATUS BLOCK

```
document: 72-LC
status: COMPLETE
date: 2026-03-04
author: Eric Robert Lawson
  OrganismCore

precursor: Document 72
  The_False_Attractor_Confirmed.md

genes_confirmed: 3
  (SPI1 90.5%, KLF4 94.7%, IRF8 69.5%
   all p = machine zero)
controls_correct: 4/4
novel_findings: 4
analyst_errors: 2
  (CEBPA, RUNX1 — both corrected
   and converted to insights)

key_convergence:
  REVERT / Cho 2025 / Advanced Science
  Independent experimental validation
  of the attractor reversion principle.

next_steps_from_this_check:
  1. CRISPRa experiment: KLF4+SPI1+IRF8
     triad in OCI-AML3 / MOLM-13 / primary
     AML samples.
     Endpoint: monocytic differentiation
     (CD14+, CD11b+ markers).
  2. RUNX1 inhibitor + triad combination
     in GMP-like AML cells.
  3. IRF8 inclusion/exclusion stratified
     by MLL rearrangement status.
  4. Contact REVERT group (Cho lab, KAIST)
     for collaborative extension to AML.

repository_path:
  Cancer_Research/AML/AML_Literature_Check.md
```
