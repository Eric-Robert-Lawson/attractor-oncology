# MYELODYSPLASTIC NEOPLASMS (MDS) — SUBTYPE ORIENTATION DOCUMENT
## Before Any Subtype Analysis Begins
## OrganismCore | 2026-03-03 | Author: Eric Robert Lawson

---

## PURPOSE OF THIS DOCUMENT

```
This document exists before any script runs.
It contains no predictions.
It contains no epigenetic hypotheses.
It contains no depth score derivations.

What it contains:

  A complete map of the myelodysplastic neoplasm
  (MDS) molecular subtype landscape — the WHO 2022
  classification, the three molecularly-defined
  entities (SF3B1, del(5q), biallelic TP53), the
  broader genomic classes (spliceosome, epigenetic
  modifier, transcription factor, cohesion-complex
  mutations), the CHIP-to-MDS-to-AML progression
  continuum, the normal myeloid differentiation
  hierarchy that is being disrupted, and the full
  treatment landscape from supportive care through
  hypomethylating agents to transplant.

MDS HOLDS A UNIQUE STRUCTURAL POSITION
IN THE ORGANISCORE REPOSITORY.

  Unlike every other cancer in the repository,
  MDS is not a disease of ONE cell type acquiring
  a false attractor.
  MDS is a disease of the ENTIRE MYELOID
  DIFFERENTIATION HIERARCHY being disrupted
  simultaneously at the stem/progenitor level —
  affecting multiple downstream lineages at once
  (erythroid, granulocytic, megakaryocytic).

  THE MDS FALSE ATTRACTOR IS IN THE STEM CELL.

  A hematopoietic stem cell (HSC) or early
  multipotent progenitor acquires a founding
  mutation (typically TET2, DNMT3A, or ASXL1 —
  the CHIP-stage mutations).
  This cell is not yet a cancer cell — it is a
  PRENEOPLASTIC STEM CELL clonally expanding
  in the bone marrow (CHIP — clonal hematopoiesis
  of indeterminate potential).
  Additional hits convert CHIP to MDS:
    A spliceosome mutation (SF3B1, SRSF2, U2AF1)
    OR a transcription factor mutation (RUNX1,
       ETV6, GATA2)
    OR a chromosomal lesion (del(5q), monosomy 7,
       del(7q), del(20q))
    OR TP53 mutation / deletion
  The resulting MDS clone produces DYSPLASTIC,
  DYSFUNCTIONAL blood cells — cells that look
  wrong and don't work — while failing to produce
  ENOUGH cells (cytopenias).

  THE DEPTH AXIS IN MDS IS DIFFERENT FROM
  ALL OTHER CANCERS IN THE REPOSITORY.

  In PRAD, PAAD, STAD, BRCA: depth measures
  distance from a TERMINAL DIFFERENTIATED CELL
  identity (luminal prostate, acinar pancreas,
  gastric chief/parietal, luminal breast).
  The cancer cell has traveled from its normal
  terminal attractor toward a false attractor.

  In MDS: depth measures distance from a NORMAL
  STEM CELL IDENTITY — not a terminal
  differentiated cell identity.
  The MDS stem cell is STUCK at an early
  progenitor state, producing dysplastic output,
  unable to differentiate properly.
  But the dysplastic cells it produces are
  themselves ALSO in a false attractor:
  a promyelocyte that cannot complete
  granulocytic maturation;
  an erythroblast that cannot complete
  erythroid maturation;
  a megakaryocyte that cannot complete
  platelet production.

  MDS IS A FALSE ATTRACTOR OF FALSE ATTRACTORS:
  The stem cell is stuck (first false attractor).
  Its progeny are stuck (derived false attractors).
  Dysplasia = the phenotype of cells caught
  between two attractors — not enough signal to
  complete differentiation, not enough signal to
  remain at stem identity.
  Caught in the Waddington valley between
  attractors.

  THE EXISTING MDS ANALYSIS IN THE REPOSITORY:
  The completed bulk MDS analysis identified
  the switch gene and drug target in the bulk
  MDS signal.
  THIS SUBTYPE SERIES exists to ask:
  WHICH MOLECULAR SUBTYPE drives the bulk signal?
  Do different subtypes have different switch genes?
  Does the depth score vary by molecular subtype?
  Which subtype most benefits from which treatment?

  THE STRUCTURAL IMPLICATION OF THE MDS FINDING:
  Unlike AML (where the differentiation block is
  at a specific stage — myelocyte/promyelocyte —
  and is driven by a specific fusion gene like
  PML-RARA or AML1-ETO), MDS has NO SINGLE BLOCK.
  It has a GENERAL FAILURE OF DIFFERENTIATION
  arising from stem cell dysfunction.
  This means the depth score in MDS measures
  something different from depth scores in
  solid tumours:
  NOT "how far has a cell traveled from its
  terminal identity?"
  BUT "how completely has the stem cell
  lost the ability to execute differentiation
  programmes?"
  The depth axis in MDS is a DIFFERENTIATION
  EXECUTION CAPACITY axis — how much of the normal
  myeloid differentiation circuitry is still
  functional in the MDS stem cell clone?
```

---

## DOCUMENT METADATA

```
document_id:        MDS_Subtype_Orientation
series:             MDS (Myelodysplastic Neoplasms —
                    Subtypes)
folder:             Cancer_Research/MDS/Subtypes/
date:               2026-03-03
author:             Eric Robert Lawson
status:             ORIENTATION ONLY — no predictions,
                    no analysis, no epigenetic hypotheses
next_document:      MDS_SF3B1_before.md
                    (Document MDS-S1a — SF3B1-mutant
                    before-doc, the most benign and
                    best-characterised molecular subtype)
protocol_version:   Workflow_Protocol.md v2.0
existing_analysis:  See OrganismCore cancer series
                    for the completed MDS bulk analysis.
                    Switch gene and drug target identified.
                    This subtype series resolves which
                    subtype drives the bulk signal and
                    whether different subtypes have
                    different switch genes.
```

---

## SECTION I — THE NORMAL MYELOID HIERARCHY

```
The Waddington baseline for MDS is the normal
bone marrow myeloid differentiation hierarchy —
the cascade from HSC to mature blood cells that
operates continuously throughout life.

SCALE OF NORMAL MYELOPOIESIS:
  The adult human produces:
    ~200 billion red blood cells per day
    ~100 billion neutrophils per day
    ~100 billion platelets per day
  This is the most productive tissue in the body.
  Any disruption at the stem cell level propagates
  across this enormous output — even a 10%
  reduction in differentiation efficiency
  produces clinically significant cytopenias.

═══════════════════════════════════════════════════════
THE MYELOID DIFFERENTIATION HIERARCHY — FULL MAP
═══════════════════════════════════════════════════════

LEVEL 1: HEMATOPOIETIC STEM CELL (HSC)
  Location:   Bone marrow niche (endosteal and
              perivascular zones).
              Low oxygen — hypoxic niche.
  Function:   Multipotent self-renewal.
              Produces all blood cell types.
              The founding cell type of MDS:
              the MDS-initiating cell is at
              this level or the MPP level.
  Markers:    CD34+ CD38- CD90+ CD45RA-
              LT-HSC (long-term): CD34+ CD38- CD49f+
  TF network: GATA2 (HSC survival + identity)
              RUNX1 (AML1) — required for definitive
                             haematopoiesis
              GFI1 — transcriptional repressor,
                     required for HSC quiescence
              ERG, ETV6 — ETS factors for HSC
              TAL1/SCL — bHLH factor for HSC
  The normal HSC is the deepest normal attractor
  in the myeloid hierarchy — it is the SOURCE
  attractor from which all myeloid identities
  are derived by progressive restriction.
  MDS DISRUPTS THIS SOURCE ATTRACTOR.
  The MDS HSC acquires a false programme:
  it can still self-renew (the false attractor
  is STABLE) but cannot efficiently execute
  differentiation (the downstream circuits are
  disrupted by the driver mutations).

LEVEL 2: MULTIPOTENT PROGENITOR (MPP)
  Subclasses:
    MPP1 (self-renewing, myeloid-biased)
    MPP2 (myeloid biased)
    MPP3 (myeloid/lymphoid balanced)
    MPP4 (lymphoid primed)
  MDS predominantly affects MPP1-MPP3
  (myeloid-biased progenitors).
  Lymphoid output is relatively SPARED in MDS —
  B cells and T cells are usually near-normal
  in count (unlike in AML where lymphoid
  output can also be suppressed).

LEVEL 3: COMMON MYELOID PROGENITOR (CMP)
  The commitment point: CMP can only give rise
  to myeloid cells (not lymphoid).
  Divides into:
    GMP (granulocyte-monocyte progenitor)
    MEP (megakaryocyte-erythroid progenitor)
  Key TF switch at CMP level:
    HIGH PU.1 (SPI1) → GMP fate
    HIGH GATA1 → MEP fate
    THE GATA1/PU.1 ANTAGONISM:
    GATA1 and PU.1 are mutual transcriptional
    repressors — high GATA1 represses PU.1 target
    genes; high PU.1 represses GATA1 targets.
    The balance between GATA1 and PU.1 at the
    CMP level determines whether a cell becomes
    a granulocyte/monocyte (GMP) or an erythroid
    cell/megakaryocyte (MEP).
    In MDS: SRSF2 and U2AF1 mutations tend to
    SKEW TOWARD GMP (granulocytic bias with
    dysplasia), while SF3B1 mutations SKEW
    TOWARD ERYTHROID DYSPLASIA.
    This lineage-skewing is a structural feature
    of MDS subtypes that the depth score must
    capture.

LEVEL 4a: GRANULOCYTE-MONOCYTE PROGENITOR (GMP)
  Gives rise to:
    Neutrophils (most abundant granulocyte)
    Eosinophils
    Basophils / Mast cells
    Monocytes / Macrophages / Dendritic cells
  Maturation stages of granulopoiesis:
    GMP →
    Promyelocyte (granule formation begins) →
    Myelocyte (secondary granule formation) →
    Metamyelocyte →
    Band cell →
    Mature neutrophil
  Normal transit: ~7–14 days GMP to neutrophil.
  KEY TF NETWORK FOR GRANULOPOIESIS:
    CEBPα (CEBPA) — THE MASTER GRANULOCYTIC TF
      Drives GMP-to-promyelocyte transition.
      Activates: G-CSF receptor (CSF3R),
      granule genes (ELANE, MPO, PRTN3),
      and granulocytic identity.
      CEBPA LOSS → AML (acute promyelocytic-like
      differentiation block at promyelocyte).
      CEBPA HYPOMORPHY → MDS (partial
      granulopoiesis — dysplastic neutrophils
      with nuclear hyposegmentation, toxic
      granulation, Pelger-Huët anomaly).
    PU.1 (SPI1) — upstream driver of myeloid fate
      Activates CEBPA.
      In MDS: PU.1 expression is maintained but
      its downstream targets may be mis-spliced
      (if SRSF2/U2AF1 mutation) or
      epigenetically silenced (if ASXL1 mutation).
    GFI1 — repressor of monocyte programme
      in granulopoiesis.
      GFI1 mutation → monocyte-biased output.
    ELANE (Neutrophil Elastase) — a GFI1 target
      Expressed in promyelocytes.
      Dysregulation contributes to neutrophil
      maturation block in some MDS subtypes.

LEVEL 4b: MEGAKARYOCYTE-ERYTHROID PROGENITOR (MEP)
  Gives rise to:
    Erythroid lineage:
      MEP → BFU-E (burst-forming unit erythroid)
           → CFU-E (colony-forming unit erythroid)
           → Proerythroblast
           → Basophilic erythroblast
           → Polychromatophilic erythroblast
           → Orthochromatophilic erythroblast
           → Reticulocyte
           → Mature red blood cell (erythrocyte)
      Normal transit: ~7 days BFU-E to reticulocyte.
    Megakaryocyte lineage:
      MEP → Megakaryoblast → Promegakaryocyte
           → Megakaryocyte (polyploid, endomitosis)
           → Platelet release (thrombopoiesis)
  KEY TF NETWORK FOR ERYTHROPOIESIS:
    GATA1 — THE MASTER ERYTHROID TF
      Drives MEP-to-erythroid commitment and
      all subsequent erythroid maturation.
      Activates: HBB (beta-globin), HBA1/2
      (alpha-globin), KLF1, ALAS2 (haem synthesis
      rate-limiting enzyme), SLC4A1 (band 3),
      GYPA (glycophorin A).
      GATA1 LOSS → DYSERYTHROPOIESIS (the cardinal
      feature of MDS erythroid dysplasia).
      In MDS: GATA1 is frequently epigenetically
      silenced (EZH2 high in some MDS subtypes)
      OR its target genes are mis-spliced (SF3B1
      mutation specifically disrupts ABCB7 and
      TMEM14C — both required for haem synthesis).
    KLF1 — erythroid-specific TF activated by GATA1
      Activates: HBB, ALAS2, BCL2L1 (Bcl-xL —
      erythroid survival signal).
      KLF1 LOSS → erythroid maturation failure.
    EPOR (erythropoietin receptor) — EPO signalling
      drives erythroid differentiation via JAK2/STAT5.
      STAT5 activates BCL2L1 (anti-apoptotic).
      In MDS: EPO response is impaired — the
      dyserythropoietic progenitors do not
      respond normally to EPO even when EPO
      levels are high.
      This is the pathophysiology behind ESA
      resistance in MDS (exogenous EPO doesn't
      work because the differentiation circuit
      is broken downstream of the receptor).
    FOG1 (ZFPM1) — GATA1 cofactor required for
      erythroid and megakaryocytic maturation.

  THE SF3B1 → RING SIDEROBLAST MECHANISM:
    SF3B1 K700E (most common mutation) alters
    U2 snRNP splicing to use cryptic 3' splice
    sites → mis-splicing of hundreds of genes.
    The two critical mis-spliced genes in erythroid
    precursors:
    ABCB7 (ATP-binding cassette transporter B7):
      Exports iron-sulfur (Fe-S) clusters from
      mitochondria. Mis-splicing → reduced
      functional ABCB7 → Fe-S clusters
      accumulate in mitochondria.
    TMEM14C (transmembrane protein 14C):
      Required for mitochondrial haem synthesis.
      Mis-splicing → reduced TMEM14C → haem
      synthesis impaired.
    Combined effect:
      Iron enters mitochondria (via SLC25A37/
      mitoferrin-1) but CANNOT be incorporated
      into haem (TMEM14C loss) AND cannot be
      exported as Fe-S clusters (ABCB7 loss).
      Iron ACCUMULATES in the mitochondria of
      erythroid precursors → perinuclear iron
      deposits visible as Prussian blue-staining
      granules around the nucleus (≥15% = ring
      sideroblast diagnostic criterion).
    This is a PARTIAL DIFFERENTIATION ARREST:
      The erythroid cell begins differentiation
      (GATA1 activates haem synthesis pathway)
      but CANNOT COMPLETE it (mitochondrial iron
      trapping prevents haem production).
      The cell is caught BETWEEN normal erythroid
      attractor positions — it has started the
      programme but cannot finish it.
      THIS IS THE WADDINGTON INTERPRETATION
      OF RING SIDEROBLASTS:
      They are cells trapped in the valley
      between the MEP attractor and the mature
      erythrocyte attractor, held there by
      mitochondrial iron trapping.
      The false attractor in SF3B1-MDS is an
      INTERMEDIATE STATE — partially erythroid,
      partially dysfunctional — not a completely
      foreign identity (like basal-like PDAC or
      NEPC) but a STALLED programme.

THE WADDINGTON STRUCTURE OF NORMAL MYELOID
DIFFERENTIATION (SUMMARY):

  HSC (multipotent attractor)
    ↓ MYELOID COMMITMENT (RUNX1, GATA2)
  CMP (restricted myeloid attractor)
    ↓ GATA1/PU.1 balance decides fate
  GMP → GRANULOCYTE (CEBPA drives) ← normal attractor
  MEP → ERYTHROCYTE (GATA1 drives) ← normal attractor
      → MEGAKARYOCYTE (GATA1 + FLI1 drive) ← attractor

  MDS FALSE ATTRACTORS:
  CHIP-stage: HSC in a pre-MDS false attractor
    (self-renewing, clonally expanding, but
    differentiation subtly impaired)
  LOW-RISK MDS: Stable false attractor at the
    MEP/GMP level — cells initiate differentiation
    but stall at intermediate stages.
    SF3B1 → erythroid stalling (ring sideroblasts)
    del(5q) → erythroid stalling (RPS14 loss)
    Both: REVERTIBLE with the right intervention.
  HIGH-RISK MDS: False attractor at GMP level
    with escalating genetic instability —
    deeper stalling, more dysplasia, approach
    to transformation.
    SRSF2 + ASXL1 → granulocytic dysplasia +
    monocytosis (approaching CMML).
  MDS-AML transition: The GMP-level stall
    becomes a complete differentiation block
    (>20% blasts = AML by WHO definition).
    RUNX1 mutation + TP53 = rapid transition.
```

---

## SECTION II — THE WHO 2022 CLASSIFICATION

```
The 2022 WHO classification renamed
"myelodysplastic syndromes" to
"myelodysplastic neoplasms" — retaining the
MDS abbreviation. The key structural change:
increasing emphasis on MOLECULARLY DEFINED
ENTITIES over pure morphologic categories.

THREE MOLECULARLY DEFINED ENTITIES:

1. MDS WITH LOW BLASTS AND SF3B1 MUTATION
   (MDS-SF3B1)
   See Section III — the most important
   molecular entity for circuit analysis.

2. MDS WITH LOW BLASTS AND ISOLATED del(5q)
   (MDS-del(5q))
   See Section IV.

3. MDS WITH BIALLELIC TP53 INACTIVATION
   (MDS-biTP53)
   See Section VII — the worst prognosis entity.

MORPHOLOGICALLY DEFINED ENTITIES (non-molecular):
  MDS with low blasts (MDS-LB)
  MDS, hypoplastic (MDS-h)
  MDS with increased blasts 1 (MDS-IB1: 5–9%)
  MDS with increased blasts 2 (MDS-IB2: 10–19%)
  MDS with fibrosis (MDS-f)

THE BLAST THRESHOLD HIERARCHY:
  The blast percentage in bone marrow (or blood)
  defines the MDS/AML boundary:
    <5% blasts: low blast MDS (SF3B1, del(5q),
                or non-molecular MDS-LB)
    5–9% blasts: MDS-IB1 (high risk)
    10–19%: MDS-IB2 (approaching AML)
    ≥20%: AML (by WHO 2022 definition)
            NOTE: AML previously required ≥30%
            (FAB classification). The 20% threshold
            is the current standard.
  The blast percentage is the CLINICAL DEPTH
  MEASURE in MDS — analogous to Gleason score
  in PRAD (increasing blast % = increasing depth).
  The framework depth score should CORRELATE
  with blast percentage: this is the primary
  clinical validation criterion for MDS depth.

THE IPSS-M (MOLECULAR IPSS):
  Published 2022. The most accurate prognostic
  scoring system in MDS.
  Integrates:
    Blast percentage
    Haemoglobin level
    Platelet count
    Absolute neutrophil count
    Cytogenetics (complex karyotype, monosomy 7,
    del(7q), del(17p), etc.)
    31 gene mutations (including all of the
    genes discussed in this document)
  Produces 6 risk categories:
    Very low | Low | Moderate-low | Moderate-high
    | High | Very high
  The IPSS-M is effectively a CLINICAL DEPTH
  SCORE THAT ALREADY EXISTS FOR MDS.
  The framework's depth score must be compared
  to IPSS-M in the analysis.
  If the framework depth score predicts IPSS-M
  risk category from expression data alone —
  that is a convergent validation.
  If the framework depth score provides
  additional predictive power BEYOND IPSS-M —
  that is the going-further finding.
```

---

## SECTION III — SF3B1-MUTANT MDS: THE STALLED ERYTHROID PROGRAMME

```
FREQUENCY:     ~20–30% of MDS overall.
               The most common single-gene MDS
               molecular entity.
               Essentially 100% of MDS with
               ring sideroblasts ≥15%.
               WHO 2022: defined by SF3B1 mutation
               + low blasts + ≥5% ring sideroblasts
               (in absence of complex karyotype,
               del(5q), or biallelic TP53).

PROGNOSIS:     BEST in MDS.
               Median OS: ~5–10 years
               Low risk of AML transformation (~10%
               at 10 years).
               The most indolent MDS subtype.

CELL OF ORIGIN:
  A haematopoietic stem cell or early MEP-biased
  progenitor that acquired SF3B1 K700E (or
  other hotspot mutations: K666N, K666T, R625C).
  The mutation drives a SELECTIVE ADVANTAGE at
  the HSC/progenitor level (CHIP) AND simultaneously
  disrupts erythroid differentiation (MDS phenotype).
  The SF3B1 clone expands in the bone marrow
  and progressively dominates erythroid output —
  but because the granulocyte and megakaryocyte
  lineages are RELATIVELY SPARED (the critical
  mis-spliced genes are most expressed during
  erythroid differentiation), cytopenias are
  mainly anaemia.

THE WADDINGTON INTERPRETATION:
  SF3B1-MDS is the most tractable MDS subtype
  for the framework because the false attractor
  is WELL DEFINED:
  The erythroid progenitor is STUCK at the
  basophilic-to-polychromatophilic erythroblast
  transition because mitochondrial haem synthesis
  is blocked (TMEM14C mis-spliced) and iron
  cannot be exported (ABCB7 mis-spliced).
  The stalling is a CIRCUIT EXECUTION FAILURE —
  not an identity-loss failure. The erythroid
  cell has GATA1 active. It HAS the programme.
  It cannot EXECUTE the haem synthesis step
  because of the splicing defect in the
  mitochondrial machinery.
  This is a FALSE ATTRACTOR at an INTERMEDIATE
  DIFFERENTIATION POSITION — the cell is
  not at the stem cell and not at the mature
  erythrocyte. It is stuck in between.
  THE SWITCH GENE QUESTION IN SF3B1-MDS:
  What is the single gene that, if restored,
  would allow the erythroid programme to
  complete execution?
  Structural candidates:
    ABCB7 — restoring splicing → haem synthesis
    TMEM14C — restoring splicing → haem synthesis
    ALAS2 — the rate-limiting haem enzyme
             (activated by GATA1, but limited by
             iron availability due to ABCB7/TMEM14C loss)
    FOG1 (ZFPM1) — GATA1 cofactor for erythroid
                    maturation
  THE SWITCH GENE IN SF3B1-MDS IS LIKELY AT THE
  MITOCHONDRIAL IRON METABOLISM STEP —
  the bottleneck in the circuit.

LUSPATERCEPT — THE CLINICAL CONFIRMATION:
  Luspatercept (activin receptor ligand trap,
  specifically an ACVR2B/TGFBR3 antagonist):
    Activin ligands (GDF11, GDF8, activin A/B)
    signal through SMAD2/3 and INHIBIT late
    erythroid maturation — they block the
    transition from late erythroblast to
    reticulocyte.
    In SF3B1-MDS: activin signalling is
    ABNORMALLY HIGH because the dysplastic
    erythroid progenitors secrete excess
    activin/GDF signals (a consequence of
    the stalled differentiation programme —
    stalled cells send SOS signals).
    Luspatercept blocks this excess signalling
    → allows late erythroid maturation to
    complete → reduces transfusion dependence.
    LANDMARK TRIALS:
    MEDALIST (lower-risk MDS with ring
    sideroblasts): luspatercept vs. placebo.
    38% luspatercept vs. 13% placebo achieved
    transfusion independence ≥8 weeks. (NEJM 2020)
    COMMANDS (first-line, ESA-naive, lower-risk
    MDS): luspatercept vs. epoetin alfa.
    Superiority established — luspatercept
    is now preferred over ESAs for lower-risk
    MDS with ring sideroblasts (2023–2024).
    FDA APPROVED for lower-risk MDS with
    ring sideroblasts after ESA failure (2020)
    AND now as first-line (2023).
  FRAMEWORK INTERPRETATION:
    Luspatercept is the MOST DIRECT CLINICAL
    CONFIRMATION OF THE STALLED-PROGRAMME
    ATTRACTOR CONCEPT IN MDS.
    The stalled erythroid cells are sending
    activin signals because they cannot complete
    differentiation.
    Blocking those signals ALLOWS the incomplete
    programme to partially execute.
    Luspatercept is not restoring the splicing
    defect — it is REMOVING A SECONDARY BRAKING
    SIGNAL that the stalled cells have activated.
    The primary defect (splicing) remains.
    The secondary brake (activin excess) is removed.
    The analogy: in PRAD Luminal A, EZH2 is the
    PRIMARY LOCK on NKX3-1 — tazemetostat removes
    it, the circuit executes.
    In SF3B1-MDS, TMEM14C/ABCB7 mis-splicing is
    the primary defect AND activin excess is the
    SECONDARY BRAKE — luspatercept removes the
    secondary brake, allowing partial compensation.
    A true circuit restoration in SF3B1-MDS would
    require splicing correction of ABCB7/TMEM14C —
    which spliceosome modulator drugs are attempting.

SPLICEOSOME MODULATORS (NOVEL PIPELINE):
  H3B-8800 (splicing modulator, Epizyme/Syros):
    Preferentially lethal to spliceosome-mutant
    cells (SF3B1, SRSF2, U2AF1) over wildtype.
    Mechanism: modulates U2/U12 snRNP — exploits
    the vulnerability of already-compromised
    spliceosome machinery.
    Phase I data in MDS/AML: safety established;
    modest single-agent activity.
    The framework prediction: H3B-8800 kills the
    SF3B1 clone but does not differentiate it.
    Combining H3B-8800 with luspatercept may:
    H3B-8800 reduces SF3B1 clone size →
    luspatercept allows residual erythroid cells
    to differentiate → additive anaemia improvement
    and possibly MRD reduction.
    This combination is not yet in clinical trials.
    It is the structural framework prediction
    for SF3B1-MDS beyond current standard of care.
    STATUS: GOING FURTHER — not in literature.

TREATMENT:
  Lower-risk SF3B1-MDS:
    Observation if asymptomatic anaemia.
    ESA (erythropoietin stimulating agents) if
    low serum EPO (EPO <200 U/L):
    ~30–40% response rate.
    Luspatercept: First-line for SF3B1-mutant
    lower-risk MDS with ring sideroblasts.
    Superior to ESA (COMMANDS trial, 2023).
    Imetelstat (telomerase inhibitor):
    FDA approved 2024 for lower-risk, transfusion-
    dependent, ESA-refractory MDS (non-del5q).
    Some SF3B1-mutant patients respond.
    Higher-risk SF3B1-MDS (rare — if blast count
    rises to ≥5%):
    Hypomethylating agent (azacitidine/decitabine).
    Allo-HSCT if eligible.
```

---

## SECTION IV — del(5q) MDS: THE RIBOSOMAL STRESS ATTRACTOR

```
FREQUENCY:     ~5–10% of MDS.
               More common in older women
               (F:M ratio ~2:1 for the classic
               low-risk del(5q) phenotype).
               Classic presentation: macrocytic
               anaemia + normal/elevated platelets
               + hypolobulated megakaryocytes.

PROGNOSIS:     GOOD (isolated del(5q)):
               Median OS: ~6–9 years.
               Low AML transformation risk if isolated.
               WORSENS if additional mutations present:
               del(5q) + TP53 mutation = very high risk
               del(5q) + ASXL1 = high risk
               del(5q) + RUNX1 = high risk
               These "complex del(5q)" cases behave
               like higher-risk MDS despite low blast%.

CELL OF ORIGIN:
  An HSC or early erythroid progenitor that has
  lost a segment of chromosome 5q (typically
  5q31–q33 for the commonly deleted region, CDR).

THE RPS14 MECHANISM — WADDINGTON INTERPRETATION:
  RPS14 is located at 5q33 and is within the CDR.
  RPS14 haploinsufficiency → ribosome stress:
    RPS14 is required for 40S ribosomal subunit
    maturation.
    With only one functional RPS14 allele, the
    erythroid progenitor cannot produce enough
    40S ribosomes to sustain the enormous protein
    synthesis demands of erythroid maturation
    (erythroid cells must produce ~280 million
    haemoglobin molecules per cell during
    maturation — the single highest protein
    synthesis demand of any cell in the body).
    Ribosomal stress → free RPL5/RPL11
    (large ribosomal subunit proteins) accumulate
    → RPL5/RPL11 bind and INHIBIT MDM2 →
    MDM2 cannot ubiquitinate p53 → p53 stabilises
    → p53 triggers apoptosis of erythroid
    precursors.
  The del(5q) false attractor:
    Erythroid precursors CANNOT COMPLETE
    DIFFERENTIATION because they die via
    p53-mediated apoptosis before they can
    complete the ribosome-intensive maturation
    programme.
    This is an apoptotic false attractor —
    the cell is pushed toward death before it
    can reach the mature erythrocyte attractor.
    Different from SF3B1-MDS (where cells stall
    but survive — ring sideroblasts are alive).
    In del(5q)-MDS: erythroid cells die prematurely
    → anaemia.
  THE SWITCH GENE QUESTION IN del(5q)-MDS:
    RPS14 restoration would reverse the defect.
    BUT: RPS14 is on the deleted chromosome —
    it cannot be directly restored without gene
    therapy.
    The BYPASS strategy: lenalidomide.

LENALIDOMIDE — THE CLINICAL CONFIRMATION:
  Lenalidomide (thalidomide analogue, CRBN binder):
    In del(5q) MDS: lenalidomide is selectively
    cytotoxic to del(5q) haematopoietic progenitors.
    Mechanism:
      CRBN (cereblon) binding → redirects CRL4A
      ubiquitin ligase to degrade CDC25C (cell
      cycle phosphatase), CSNK1A1 (casein kinase
      1A1 — located on 5q!).
      CSNK1A1 is haploinsufficient in del(5q)
      cells (only one copy).
      Lenalidomide further depletes CSNK1A1 →
      del(5q) cells are selectively killed (one
      copy of CSNK1A1 is not enough to survive
      lenalidomide).
      Normal cells (two CSNK1A1 copies) can
      tolerate CSNK1A1 reduction → selective
      kill of del(5q) clone.
    ADDITIONALLY: lenalidomide modulates p53/MDM2
      → reduces p53-driven apoptosis in surviving
      erythroid progenitors → erythroid maturation
      can proceed → transfusion independence.
    CLINICAL RESULTS (MDS-003 and MDS-004 trials):
      ~67% transfusion independence in del(5q) MDS.
      ~45% complete cytogenetic remission.
      The only standard of care for del(5q) MDS.
      FDA approved 2005 for del(5q) MDS.
    FRAMEWORK INTERPRETATION:
      Lenalidomide is an ATTRACTOR DISSOLUTION TOOL
      for the del(5q) clone:
      It selectively kills the false-attractor cells
      (del(5q) clone) rather than restoring their
      differentiation circuit.
      This is the attractor dissolution strategy —
      the same principle as venetoclax in CLL
      (dissolving BCL2-dependent false attractor)
      but operating via a synthetic lethality
      mechanism (del(5q) haploinsufficiency for
      CSNK1A1 creates the vulnerability).
      The distinction from SF3B1-MDS:
      SF3B1-MDS stalled programme = circuit
      restoration target (luspatercept removes
      secondary brake; spliceosome modulators
      could remove primary defect).
      del(5q)-MDS = synthetic lethality target
      (lenalidomide kills the clone via CSNK1A1
      depletion; the remaining normal clone
      regenerates normal erythropoiesis).

TP53 MUTATION IN del(5q) MDS:
  ~15–20% of del(5q) MDS patients on lenalidomide
  develop TP53 mutations over time.
  The mechanism: lenalidomide increases p53
  activity → selects for TP53-mutant cells that
  have escaped p53-mediated apoptosis → the
  TP53-mutant cells resist lenalidomide AND
  have a growth advantage → AML transformation.
  Clinical implication:
    Patients on lenalidomide for del(5q) MDS
    must have serial TP53 mutational monitoring.
    TP53 mutation emergence = stop lenalidomide,
    consider HMA or transplant.
  This is an example of THERAPY-INDUCED DEPTH
  INCREASE — lenalidomide selects for deeper
  (TP53-mutant) clones, analogous to AR inhibitor
  therapy selecting for NEPC in PRAD.
  The depth score in del(5q) patients on
  lenalidomide should track over time — rising
  depth score signals TP53-clone emergence before
  clinical progression.
  This is a prospective clinical application of
  the framework in del(5q) MDS.

TREATMENT:
  Standard first-line: Lenalidomide 10 mg/day
  (21/28 days) until response or progression.
  ESAs as bridge therapy if EPO <500 U/L.
  Allo-HSCT: for complex del(5q) (with TP53 or
  other adverse mutations), higher-risk disease,
  or lenalidomide failure.
  TP53 mutation monitoring: essential — serial
  NGS every 3–6 months on lenalidomide.
```

---

## SECTION V — SPLICEOSOME-MUTANT MDS (SRSF2, U2AF1)

```
FREQUENCY:
  SRSF2:  ~10–15% of MDS
  U2AF1:  ~8–12% of MDS

PROGNOSIS:
  SRSF2: ADVERSE — especially when co-mutated
    with ASXL1 (SRSF2/ASXL1 combination is one
    of the most adverse mutation pairs in MDS,
    ~15% 2-year OS in some series).
    Strong association with CMML overlap
    (chronic myelomonocytic leukaemia).
  U2AF1: INTERMEDIATE-TO-ADVERSE.
    U2AF1 S34F/S34Y: associated with del(20q),
    monosomy 7 — worse prognosis.
    U2AF1 Q157: associated with better prognosis.

MECHANISM:
  SRSF2 (Serine-Arginine Splicing Factor 2):
    Core component of the U1 snRNP complex.
    SRSF2 P95 mutations (P95H, P95L, P95R) are
    the defining hotspot.
    Unlike SF3B1 (which uses cryptic 3' splice
    sites), SRSF2 mutations alter the RNA
    recognition motifs → changes which exons
    are included in mRNA.
    CRITICAL TARGET: EZH2 is mis-spliced by
    SRSF2 mutation → reduced EZH2 in SRSF2-mutant
    cells → H3K27me3 reduction → activation of
    PRC2-silenced genes.
    This means SRSF2-mutant MDS has LOWER EZH2
    ACTIVITY than normal, in contrast to
    SF3B1-mutant MDS (where EZH2 can be elevated).
    SRSF2 mutation also mis-splices:
      RUNX1 (alternative isoforms affecting
      myeloid differentiation TF activity)
      BCOR (epigenetic repressor)
      KANSL1 (chromatin regulator)
    GRANULOCYTIC BIAS:
      SRSF2 mutation shifts differentiation
      toward monocytes (via altered PU.1/
      CEBPA ratio) → monocytosis is a clinical
      hallmark of SRSF2-mutant MDS.
      When monocyte count >1×10⁹/L = CMML criteria
      → SRSF2-mutant MDS often co-diagnosed as
      MDS/MPN overlap or CMML.

  U2AF1 (U2 Small Nuclear RNA Auxiliary Factor 1):
    Component of the U2AF complex that recognises
    3' splice sites (AG dinucleotide).
    S34F/S34Y: alters 3' splice site recognition
    → different exon inclusion patterns from
    SRSF2 and SF3B1 mutations.
    Critical target: STRAP (serine-threonine
    kinase receptor-associated protein) → altered
    mRNA stability and ATM pathway.
    U2AF1 S34 mutations associate with monosomy 7
    — a marker of genomic instability.
    U2AF1 Q157 mutations have a different
    splicing signature and better prognosis.

WADDINGTON INTERPRETATION:
  SRSF2/U2AF1 MDS is the most complex MDS subtype
  for the framework because the differentiation
  failure is MULTI-LINEAGE and the mis-splicing
  affects TFs directly.
  The false attractor is:
    A GMP/early myeloid cell that cannot properly
    execute either granulocytic or monocytic
    differentiation — it produces dysplastic cells
    of mixed granulocyte/monocyte identity.
    SRSF2 mutation particularly drives toward
    a monocyte-biased false attractor
    (altered PU.1 downstream targets favour
    monocyte gene programme over granulocyte).
  SWITCH GENE IN SRSF2-MDS:
    CEBPα is a structural candidate:
    CEBPA drives granulocytic identity.
    If SRSF2 mutation reduces CEBPA activity
    (through mis-splicing or CEBPA-AS lncRNA
    upregulation), the GMP cannot complete
    granulocytic maturation.
    CEBPα restoration = circuit restoration
    for the granulocytic programme.
    This is a structural hypothesis to test
    in MDS-S3a.

TREATMENT:
  Lower-risk SRSF2-mutant MDS with anaemia:
    ESA (if EPO <500 U/L)
    Luspatercept (if ring sideroblasts present —
    uncommon in SRSF2-mutant but possible)
    Imetelstat (FDA approved 2024)
  Higher-risk SRSF2/ASXL1-mutant MDS:
    Hypomethylating agents (azacitidine/decitabine)
    as the current standard — BUT:
    SRSF2/ASXL1-mutant MDS has lower HMA response
    rate (~25–30%) than other high-risk MDS subtypes.
    Allo-HSCT: indicated for eligible patients.
    HDACi combinations: vorinostat + azacitidine —
    ASXL1-mutant cells are dependent on the
    histone deacetylase pathway (ASXL1 normally
    activates HDAC activity via BAP1 complex);
    ASXL1 mutation → BAP1 loss → altered histone
    deacetylation → HDACi combinations are
    mechanistically rational.
    Phase II data (AZA + vorinostat): modest
    improvement in ASXL1-mutant mCRPC — data in
    MDS emerging.
```

---

## SECTION VI — EPIGENETIC MODIFIER MDS (TET2, DNMT3A, ASXL1, IDH1/2)

```
FREQUENCY:
  TET2:   ~20–30% of MDS (most common single mutation)
  DNMT3A: ~5–10% of MDS (more common in CHIP)
  ASXL1:  ~10–15% of MDS
  IDH1:   ~3–5% of MDS
  IDH2:   ~3–5% of MDS

These mutations are rarely disease-defining on
their own (only IDH1/2 approaches this in some
AML contexts). They are predominantly FOUNDER
MUTATIONS — the first hit in the CHIP-to-MDS
progression — or CO-MUTATIONS that amplify
risk when combined with spliceosome or TF mutations.

THE CHIP-TO-MDS PROGRESSION AXIS:
  CHIP (CHIP mutations alone: TET2, DNMT3A, ASXL1):
    Clone present in peripheral blood.
    No cytopenias. No dysplasia.
    ~0.5–1% per year progression to MDS.
    No treatment — monitoring only.
  CCUS (Clonal Cytopenia of Undetermined Significance):
    CHIP mutation + cytopenia(s).
    No morphological dysplasia.
    ~15–20% per year progression to MDS.
    Close monitoring. Consider HMA trial.
  MDS (CHIP mutation + second hit):
    Morphological dysplasia + cytopenia(s).
    MDS diagnosis established.
  AML:
    ≥20% blasts. Transformation.
  THE CHIP-TO-MDS TRANSITION:
  The second hit that converts CHIP to MDS is
  typically a spliceosome mutation (SF3B1,
  SRSF2, U2AF1), a transcription factor mutation
  (RUNX1, GATA2, ETV6), or a chromosomal lesion.
  The epigenetic modifier mutations (TET2, DNMT3A,
  ASXL1) create the PERMISSIVE EPIGENETIC STATE
  that enables the second hit to establish MDS.

TET2:
  TET2 is a DNA hydroxymethylase — converts
  5-methylcytosine (5mC) to 5-hydroxymethylcytosine
  (5hmC), the first step in active DNA demethylation.
  TET2 LOSS → failure of DNA demethylation at
  loci that should be demethylated during
  myeloid differentiation → HYPERMETHYLATION of
  differentiation-associated promoters.
  GATA1 locus is a TET2 target: TET2 loss →
  GATA1 promoter hypermethylation → reduced
  GATA1 → erythroid differentiation impaired.
  EZH2 AND TET2 — THE CONNECTION:
    TET2 and EZH2 (PRC2) are in opposition:
    TET2 demethylates DNA → opens chromatin for
    differentiation TFs.
    EZH2 adds H3K27me3 → closes chromatin,
    blocking differentiation TF access.
    When TET2 is lost: EZH2 activity is relatively
    UNOPPOSED → more H3K27me3 at differentiation
    loci → deeper differentiation block.
    This is the epigenetic mechanism by which
    TET2 loss creates the permissive state for
    MDS establishment.
    Hypomethylating agents (azacitidine, decitabine)
    work in part by REVERSING THE TET2-LOSS
    PHENOTYPE: they reduce DNA methylation,
    partially compensating for TET2 loss.
    This is the mechanistic rationale for HMA
    therapy in TET2-mutant MDS.

DNMT3A:
  DNMT3A is a de novo DNA methyltransferase —
  adds methyl groups to previously unmethylated CpG
  sites during development and lineage commitment.
  DNMT3A R882 (hotspot) — loss of function.
  DNMT3A LOSS → failure to properly methylate
  HSC self-renewal genes during commitment →
  HSC remains in an inappropriately stem-like
  state → CHIP establishment.
  DNMT3A is the most common CHIP mutation in
  humans (>40% of CHIP >70 years old).
  In MDS: DNMT3A mutation is the FOUNDER
  mutation upon which a spliceosome or TF
  mutation is layered.
  DNMT3A-mutant MDS with SF3B1: the most common
  two-hit combination in lower-risk MDS.
  DNMT3A-mutant MDS with RUNX1: higher-risk,
  rapid AML transformation.

ASXL1:
  ASXL1 (additional sex combs-like 1) is a
  Polycomb group protein that normally ACTIVATES
  BAP1 histone deubiquitinase and modulates
  PRC2 (EZH2) activity.
  ASXL1 LOSS → disrupts the PRC2 regulatory
  complex → LOSS OF H3K27me3 AT SOME LOCI
  (similar to EZH2 loss-of-function) AND
  GAIN OF H3K27me3 AT OTHER LOCI (paradox).
  The paradox: ASXL1 is required for BOTH
  activating and repressing specific loci.
  The net effect: DYSREGULATED PRC2 TARGETING
  → inappropriate silencing of myeloid
  differentiation genes AND inappropriate
  activation of self-renewal genes.
  EZH2 RELATIONSHIP:
    In ASXL1-mutant MDS, EZH2 is DYSREGULATED —
    not simply elevated or reduced, but
    mistargeted. EZH2 inhibitors in ASXL1-mutant
    MDS have shown activity (preclinical) but
    clinical benefit is modest.
    The reason: inhibiting EZH2 in ASXL1-mutant
    cells removes H3K27me3 globally, including at
    loci that benefit from H3K27me3 for
    differentiation fidelity. The dysregulated
    targeting means that EZH2 inhibition has
    mixed consequences in ASXL1-mutant cells.
  ASXL1 AS THE HIGHEST-RISK CO-MUTATION:
    ASXL1 + SRSF2 = one of the worst MDS
    mutation combinations (CMML-like, rapid AML).
    ASXL1 + RUNX1 = very high risk.
    ASXL1 + STAG2 (cohesin) = high risk.

IDH1/IDH2:
  IDH1 (isocitrate dehydrogenase 1, cytoplasmic)
  IDH2 (isocitrate dehydrogenase 2, mitochondrial)
  Gain-of-function mutations produce the
  oncometabolite 2-hydroxyglutarate (2-HG).
  2-HG inhibits TET2 (the alpha-ketoglutarate-
  dependent dioxygenase) → phenocopies TET2 loss.
  2-HG ALSO inhibits KDM family histone
  demethylases → increases H3K9me2/H3K27me3.
  IDH MUTATIONS IN MDS ARE CLINICALLY ACTIONABLE:
    Ivosidenib (IDH1 inhibitor) — Phase I/II data
    in IDH1-mutant MDS: ~40% ORR.
    Enasidenib (IDH2 inhibitor) — Phase I/II data
    in IDH2-mutant MDS: ~35% ORR.
    Both FDA approved for IDH-mutant AML.
    MDS approval is in progress.
  CROSS-REPOSITORY CONNECTION:
    IDH mutation in MDS creates a 2-HG-driven
    CIMP phenotype — identical in mechanism to
    IDH1 mutation in GBM and IDH1 in PRAD (IDH
    subtype), and related to EBV-CIMP in GS STAD
    (where EBV infection drives global DNA
    methylation via a different mechanism).
    All of these are ONCOMETABOLITE/VIRAL-DRIVEN
    EPIGENETIC REPROGRAMMING states.
    Same mechanism — TET2/KDM inhibition by
    either 2-HG (IDH mutation) or EBV viral
    factors — leading to hypermethylation and
    differentiation block.
    The depth score in IDH-mutant MDS should
    show the highest H3K27me3-associated gene
    suppression pattern (because BOTH 2-HG
    inhibiting TET2 AND inhibiting KDMs
    increases H3K27me3 at differentiation loci).

TREATMENT:
  TET2-mutant lower-risk MDS:
    HMA therapy (azacitidine/decitabine) shows
    better response in TET2-mutant vs. wildtype.
    Imetelstat (FDA 2024 for lower-risk MDS):
    some TET2-mutant patients respond.
    Vitamin C (ascorbic acid) — TET2 restorer:
    Vitamin C activates TET2 enzymatic activity
    (ascorbic acid is a cofactor for TET2 as an
    Fe(II)/alpha-ketoglutarate-dependent enzyme).
    In TET2-mutant cells with residual TET2
    function (haploinsufficient but not null):
    high-dose vitamin C may partially restore
    TET2 activity.
    Phase I/II data in TET2-mutant MDS/AML:
    modest responses.
    This is a circuit restoration strategy:
    vitamin C restores partial TET2 function →
    DNA demethylation proceeds → differentiation
    genes open → programme can execute.
    FRAMEWORK PREDICTION: In TET2-haploinsufficient
    MDS, vitamin C depth score response test:
    add vitamin C → measure depth score change
    in CD34+ cells → response prediction.
    Not yet tested in this form.
  IDH1-mutant MDS:
    Ivosidenib (IDH1 inhibitor) — active.
    Azacitidine + ivosidenib: emerging combination.
  IDH2-mutant MDS:
    Enasidenib (IDH2 inhibitor) — active.
  ASXL1-mutant higher-risk MDS:
    Azacitidine ± venetoclax (VERONA trial negative
    in unselected higher-risk MDS — the ASXL1-
    selected response data is pending analysis).
    Allo-HSCT strongly recommended.
```

---

## SECTION VII — BIALLELIC TP53-MDS: THE GENOME GUARDIAN LOST

```
FREQUENCY:     ~4–5% of MDS.
               ~20–30% of therapy-related MDS
               (t-MDS — arising after cytotoxic
               chemotherapy or radiation for a
               prior malignancy).
               The most common molecular subtype
               in therapy-related MDS.

PROGNOSIS:     THE WORST IN MDS.
               Median OS: ~6–9 months without
               transplant.
               AML transformation: rapid (~6 months).
               WHO 2022: defined by BIALLELIC
               (two-hit) TP53 inactivation — either
               two separate mutations, or one mutation
               + deletion of the other allele (LOH).
               MONOALLELIC TP53 mutation alone is
               NOT a defining criterion for MDS-biTP53
               and has intermediate prognosis.
               Only biallelic = the worst entity.

WADDINGTON INTERPRETATION:
  TP53 biallelic inactivation in MDS is the
  closest parallel to TP53 mutation/deletion
  in NEPC (PRAD) and to TP53 loss in basal-like
  PDAC — it is the GENOME GUARDIAN LOSS that
  enables lineage plasticity and identity loss.
  In MDS: TP53 loss allows the stem/progenitor
  cell to accumulate chromosomal instability
  (complex karyotype: ≥3 structural abnormalities
  including del(5q), del(7q), del(17p), +8).
  The MDS-biTP53 false attractor is characterised
  by:
    Maximum chromosomal instability
    Maximum apoptosis resistance
    Maximum proliferative advantage for the
    TP53-null clone
    Minimum differentiation capacity
    This is the deepest false attractor in MDS —
    the MDS equivalent of DNPC in PRAD or
    double-negative CRC.
  THE CIRCUIT IN MDS-biTP53:
    The circuit is not merely broken — it has
    been rendered non-functional by chromosomal
    instability.
    Multiple differentiation TFs are deleted
    or disrupted by the complex karyotype.
    RUNX1 is on chromosome 21q — deleted in some
    complex karyotype cases.
    GATA2 is on chromosome 3q21 — affected by
    inv(3) or del(3q).
    ETV6 is on chromosome 12p13 — commonly
    deleted in MDS.
    Circuit restoration is not feasible in
    MDS-biTP53.
    ATTRACTOR DISSOLUTION is the strategy —
    selectively kill the TP53-null clone before
    it becomes AML.

EPRENETAPOPT (APR-246):
  APR-246 is a small molecule that binds mutant
  p53 protein and restores its wild-type
  conformation — re-activating p53 function
  in cells with GOF TP53 mutations.
  In MDS-biTP53: the mechanism is more
  complex because biallelic loss includes
  DELETION of one allele (no protein to
  restore). APR-246 works best in MUTANT/MUTANT
  biallelic TP53 (both alleles mutated, both
  proteins present and restorable) rather than
  MUTANT/DELETED (one allele deleted — no
  protein from the deleted allele).
  CLINICAL DATA:
    Phase II (azacitidine + APR-246):
    CR rate 41%, ORR 69%, median OS ~12 months.
    Phase III (NCT03745716):
    Did NOT meet primary endpoint vs. azacitidine
    alone — CR rate and OS improvement not
    statistically significant.
  FRAMEWORK INTERPRETATION:
    The Phase III failure parallels the CELLO-1
    (tazemetostat in mCRPC) failure:
    APR-246 in UNSELECTED biTP53 MDS — negative.
    APR-246 in MUTANT/MUTANT (restorable p53)
    biTP53 MDS — potentially positive.
    Patient selection by TP53 mutation TYPE
    (both alleles mutated vs. one mutant/one
    deleted) is the depth score equivalent here:
    the mutation type = the patient selection
    criterion.
    MRD monitoring by duplex sequencing post-
    transplant + APR-246 maintenance shows
    promising data (ASH 2024) — this is the
    correct patient population and timing.
  TREATMENT:
    Azacitidine or decitabine as bridge to
    transplant (HMA shows ~30–40% CR rate in
    biTP53 MDS).
    APR-246 + azacitidine: investigational, not
    standard.
    Allo-HSCT: only potentially curative option.
    Outcomes remain poor (<25% OS at 2 years
    post-transplant).
    New: decitabine + cedazuridine (oral HMA —
    ASTX727) approved for MDS generally. Being
    tested in biTP53 subtype.
    MAGROLIMAB (anti-CD47) + azacitidine:
    ENHANCE Phase III trial — NEGATIVE and HALTED
    (higher death risk in magrolimab arm vs.
    placebo, OS 15.9 vs. 18.6 months).
    SABATOLIMAB (anti-TIM-3) + azacitidine:
    STIMULUS-MDS2 Phase III — NEGATIVE.
    VENETOCLAX + azacitidine: VERONA Phase III
    in higher-risk MDS — NEGATIVE.
    Summary of the 2024 landscape in higher-risk
    MDS:
    FOUR major Phase III combination trials all
    failed in 2023–2024:
    ENHANCE (magrolimab) — negative
    STIMULUS-MDS2 (sabatolimab) — negative
    VERONA (venetoclax) — negative
    CELLO-1 (tazemetostat in CRPC — parallel) —
    negative
    The pattern: all four are UNSELECTED patient
    population trials. All four have a rational
    mechanism but failed because the patient
    selection did not use the molecular
    stratification that the mechanism requires.
    This is the STRONGEST POSSIBLE ARGUMENT
    FOR THE DEPTH SCORE FRAMEWORK:
    Four consecutive Phase III failures in MDS
    and PRAD for mechanistically rational drugs —
    all attributable to the same root cause:
    treating the wrong patients with the right drug.
    The depth score is the patient selection tool
    that would have stratified these trials.
```

---

## SECTION VIII — TRANSCRIPTION FACTOR MUTANT MDS (RUNX1, GATA2, ETV6)

```
RUNX1 MUTATION (~10% of MDS):
  RUNX1 (AML1) is the most commonly mutated
  transcription factor in myeloid malignancies.
  It is required for:
    HSC to CMP transition
    GMP specification
    Both granulocytic and megakaryocytic lineages
  RUNX1 MUTATION IN MDS:
    Loss-of-function mutations → impaired myeloid
    differentiation at multiple stages simultaneously.
    High risk of AML transformation (~30% within
    2 years of diagnosis).
    Associated with platelet dysfunction
    (inherited RUNX1 mutations cause familial
    platelet disorder with predisposition to AML).
  SWITCH GENE CONTEXT:
    RUNX1 IS a switch gene candidate in MDS —
    the existing bulk MDS analysis (OrganismCore)
    may have identified RUNX1 or a RUNX1 target
    as the switch gene.
    RUNX1 restoration = the circuit restoration
    strategy in RUNX1-mutant MDS.
    BUT: RUNX1 mutations in MDS are
    loss-of-function — unlike EZH2 in PRAD (a
    gain-of-function lock that can be inhibited),
    RUNX1 is lost and cannot be restored by
    inhibiting an opposing pathway.
    The therapeutic approach is therefore:
    Find what RUNX1 normally represses (the
    false attractor programme) and inhibit THAT.
    OR: Find what RUNX1 normally activates and
    supplement THAT.
    This is attractor dissolution from the
    opposite direction — not blocking the lock,
    but activating the alternative pathway.

GATA2 MUTATION/DELETION (~5% of MDS,
higher in young patients):
  GATA2 is required for HSC survival and
  early progenitor maintenance.
  GATA2 haploinsufficiency syndrome:
    An inherited disorder characterised by:
      MDS/AML predisposition
      Lymphoedema
      Mycobacterial infections
      Deafness
    GATA2-haploinsufficient patients accumulate
    CHIP mutations (ASXL1, RUNX1) that convert
    their haploinsufficient HSCs into MDS.
    HSC transplant is curative before MDS.
  In MDS: GATA2 loss → reduced HSC output →
  progressive cytopenias starting in young
  adulthood.
  The switch gene here IS GATA2 — but it is
  haploinsufficient (one allele lost), not
  epigenetically silenced.
  Gene therapy approaches for GATA2 deficiency
  are in development.

ETV6 MUTATION (~5% of MDS):
  ETV6 (TEL) is an ETS family transcriptional
  repressor.
  Located at 12p13 — frequently deleted in
  complex karyotype MDS.
  ETV6 regulates:
    HSC quiescence
    MEP differentiation
    Platelet production
  ETV6 loss → increased HSC cycling (loss of
  quiescence) → increased accumulation of
  driver mutations → accelerated CHIP-to-MDS.
  ETV6 mutations are also the founding events
  in childhood B-ALL (ETV6-RUNX1 fusion) —
  demonstrating the dual role of ETV6 in
  myeloid and lymphoid contexts.
```

---

## SECTION IX — CHROMOSOMAL LESIONS IN MDS

```
CYTOGENETICS IS THE BACKBONE OF MDS PROGNOSIS.
The conventional karyotype (bone marrow
metaphase cytogenetics) was the primary
stratification tool before molecular profiling.
Even in the IPSS-M era, cytogenetics remains
essential because chromosomal lesions carry
prognostic information not fully captured by
single-gene mutations.

KEY LESIONS:

del(5q) — GOOD PROGNOSIS (isolated):
  Covered in Section IV. Lenalidomide target.
  GOOD when isolated; BAD when complex.

del(20q) — GOOD TO INTERMEDIATE:
  Isolated del(20q): generally favourable.
  ASXL1 is NOT on chromosome 20q — del(20q)
  does not directly affect ASXL1.
  The mechanism of del(20q) MDS is not as
  well characterised as del(5q).
  Associated with thrombocytopenia.
  JAK2 V617F occasionally co-occurs.

+8 (trisomy 8) — INTERMEDIATE:
  Common in both MDS and AML.
  MYC is on chromosome 8q24 — trisomy 8 →
  MYC gain → proliferative advantage.
  MYC depth correlation: trisomy 8 MDS should
  have elevated MYC expression → depth-positive.
  Associated with immune dysregulation:
  oligoclonal T cells reactive against trisomy 8
  cells (unique immune surveillance mechanism
  explaining why some trisomy 8 MDS responds to
  immunosuppressive therapy like cyclosporin).

del(7q) / monosomy 7 — POOR PROGNOSIS:
  Monosomy 7 is among the highest-risk
  cytogenetic lesions in MDS.
  Mechanism: loss of chromosome 7q material →
  loss of multiple myeloid regulatory genes:
    EZH2 is on chromosome 7q36 —
    MONOSOMY 7 = EZH2 HAPLOINSUFFICIENCY.
    In MDS with monosomy 7: EZH2 is REDUCED
    (one copy lost), in contrast to MDS with
    SF3B1 or ASXL1 mutations where EZH2 may
    be elevated.
    This is the CRITICAL EZH2 DIRECTION
    DETERMINATION (Framework Pattern 4):
    MONOSOMY 7 → EZH2 LOW → EZH2 INHIBITORS
    WOULD WORSEN DISEASE (the EZH2 lesson
    from the Framework: EZH2 direction must be
    determined from data for each cancer).
    EZH2 haploinsufficiency in monosomy 7 MDS
    is the exact inverse of EZH2 overexpression
    in BRCA/PRAD/PAAD.
    MYELOID TUMOUR SUPPRESSOR ON 7q:
    CUX1 (7q22) — transcriptional repressor
    required for HSC quiescence.
    CUX1 haploinsufficiency → increased HSC
    cycling → clone expansion.
    NF1 haploinsufficiency (NF1 on 17q —
    co-deleted in complex karyotype with monosomy 7
    cases): RAS pathway activation.
  U2AF1 S34 mutations strongly associate with
  del(7q)/monosomy 7 — providing a combined
  cytogenetic + molecular signal.
  TREATMENT:
    Allo-HSCT strongly recommended for monosomy 7
    MDS — the only potentially curative option.
    HMA as bridge.

Complex karyotype (≥3 abnormalities):
  By definition: high/very-high risk.
  ~70% of complex karyotype MDS has TP53 mutation.
  Complex karyotype itself = a surrogate for
  TP53 inactivation (which enables the
  chromosomal instability that creates complexity).
  Treatment: as per MDS-biTP53 (Section VII).
  Allo-HSCT is the only curative option.
```

---

## SECTION X — THE DEPTH AXIS PROBLEM IN MDS

```
THE UNIQUE METHODOLOGICAL CHALLENGE IN MDS:
MDS is not a solid tumour. The analysis of MDS
using the framework requires understanding several
structural differences from solid tumour analysis.

DIFFERENCE 1: THE NORMAL REFERENCE IS
THE SAME TISSUE TYPE (BONE MARROW).
  In solid tumour analysis:
  Normal = adjacent tissue or GTEx sample.
  Tumour = cancer cells from the same patient.
  The depth axis = distance from normal tissue.
  In MDS:
  Normal = CD34+ HSCs from healthy donors.
  MDS = CD34+ cells from MDS patients.
  The depth axis = distance from normal HSC gene
  expression to MDS progenitor gene expression.
  The challenge: BOTH are bone marrow CD34+ cells.
  The gene expression differences are more subtle
  than in solid tumours (which have radical
  identity changes).
  The depth score in MDS must be derived from
  CD34+-enriched samples (not unsorted bone marrow)
  to get the progenitor signal without dilution
  by mature blood cells.

DIFFERENCE 2: MULTI-LINEAGE DYSPLASIA DILUTES
THE SIGNAL.
  In PRAD bulk analysis: ~60% luminal cancer cells.
  In MDS bone marrow: the sample contains
  CD34+ progenitors (~1–5% of total cells) PLUS
  maturing granulocytes, monocytes, erythroid cells,
  and megakaryocytes at various differentiation stages.
  If unsorted bulk bone marrow is used:
    The dominant signal is from MATURE cells —
    neutrophils, red cells, etc.
    The MDS progenitor signal is a minority.
  SOLUTION: Use CD34+ enriched datasets
  (GSE114922, GSE19429 — both are CD34+
  purified samples from MDS patients and healthy
  controls).
  This is mandatory for MDS, more critical than
  tumour purity filtering in PAAD.
  ALL MDS ANALYSIS IN THIS SERIES MUST USE
  CD34+ ENRICHED SAMPLES.

DIFFERENCE 3: THE EZH2 DIRECTION IS SUBTYPE-SPECIFIC.
  The existing MDS bulk analysis found a
  convergence node.
  But EZH2 direction in MDS varies by subtype:
    SF3B1-mutant MDS: EZH2 may be ELEVATED
    (gain-of-function lock on erythroid genes)
    SRSF2-mutant MDS: EZH2 may be REDUCED
    (SRSF2 mutation mis-splices EZH2)
    Monosomy 7 MDS: EZH2 is HAPLOINSUFFICIENT
    (gene deletion — one copy lost)
    ASXL1-mutant MDS: EZH2 is DYSREGULATED
    (mistargeted — not simply high or low)
  FRAMEWORK PATTERN 4 APPLIED TO MDS:
  "EZH2 is not a universal oncogene. The direction
  of EZH2 change must be determined from data for
  each cancer. It cannot be assumed."
  In MDS, it must be determined for each SUBTYPE.
  EZH2 inhibition (tazemetostat) is appropriate
  ONLY in subtypes where EZH2 is elevated as a
  gain-of-function lock.
  EZH2 inhibition in monosomy 7 MDS (where EZH2
  is already haploinsufficient) would worsen
  disease — confirmed by the GBM lesson (where
  EZH2 direction determined treatment choice).
  THE SUBTYPE SERIES RESOLVES EZH2 DIRECTION
  FOR EACH MOLECULAR CLASS.

DIFFERENCE 4: THE DEPTH AXIS IS DIFFERENTIATION
EXECUTION CAPACITY, NOT IDENTITY LOSS.
  In solid tumours: deeper = more identity lost
  (more PTF1A, NKX3-1, or other master TF
  expression lost).
  In MDS: deeper = LESS CAPACITY TO EXECUTE
  DIFFERENTIATION PROGRAMMES.
  The depth score markers should be:
  DEPTH-NEGATIVE (falling with depth):
    Genes that mark SUCCESSFUL differentiation:
    GATA1 targets (when erythroid programme failing)
    CEBPA targets (when granulocytic programme failing)
    KLF1, HBB, GYPA (erythroid)
    ELANE, MPO, G-CSF-R (granulocytic)
    Mature blood cell markers
  DEPTH-POSITIVE (rising with depth):
    Genes that mark PROGENITOR STALLING:
    CD34, CD38, HOXA genes (stem/progenitor markers)
    MYC, MCM genes (proliferation)
    BIRC5 (survivin — apoptosis resistance)
    CD47 (anti-phagocytic signal — "don't eat me")
    CD117 (c-KIT) — stem cell survival receptor
  THE CLINICAL TRANSLATION:
  The MDS depth score measures the DEGREE OF
  DIFFERENTIATION EXECUTION FAILURE in the
  CD34+ progenitor compartment.
  High depth score = less differentiation capacity
  = higher blast proportion = higher IPSS-M risk
  = worse prognosis.
  This makes the depth score a continuous-variable
  version of the IPSS-M blast percentage component.
```

---

## SECTION XI — DATA AVAILABILITY SUMMARY

```
Dataset         Accession    Samples              Notes      Power
──────────────────────────────────────────────────────────────────────
GSE19429        GEO          183 MDS patients     CD34+      HIGH
                             17 healthy controls  enriched
                             Subtypes:            Microarray
                             del(5q), +8,
                             del(7q)/−7
                             normal karyotype

GSE114922       GEO          82 MDS + controls    CD34+      HIGH
                             Splicing mutations   RNA-seq
                             (SF3B1, SRSF2,       Mutation-
                             U2AF1 labelled)      annotated

GSE63569        GEO          SF3B1 mutant vs.     CD34+      MOD
                             wildtype + controls  RNA-seq

GSE140101       GEO          MDS vs. healthy      MSC        MOD
                             mesenchymal stem     (stromal)
                             cells                RNA-seq
                             NOTE: MSC not
                             haematopoietic —
                             stromal axis only

GSE253355       GEO (2024)   Healthy human BM     scRNA-seq  HIGH
                             normal reference     (normal
                             all cell types       only)

BEAT AML        Oregon H&U   AML+MDS patients     Multi-omic MOD-HIGH
                             ex vivo drug
                             sensitivity

BeatAML         GEO/dbGaP    >400 AML/MDS         RNA-seq    HIGH
                             CD34+ bulk

IPSS-M          Web tool     Not a GEO dataset    Clinical   N/A
Calculator      (Bernard     — uses mutation
                2022)        and cytogenetic
                             data for risk score

CRITICAL NOTES:

NOTE 1 — CD34+ ENRICHMENT IS MANDATORY:
  ALL MDS DEPTH SCORE ANALYSIS MUST USE
  CD34+ ENRICHED SAMPLES.
  Unsorted bone marrow will produce a signal
  dominated by mature cells (granulocytes,
  red cells) — the progenitor MDS signal will
  be a minority and the depth score will be
  confounded by mature-cell gene expression.
  GSE19429 and GSE114922 are both CD34+-enriched
  — these are the primary datasets.
  The universal_discovery_start_script.py must
  be configured to use these datasets, NOT
  unsorted bone marrow datasets.

NOTE 2 — MUTATION ANNOTATION IS REQUIRED:
  Unlike solid tumour datasets (TCGA-PAAD,
  TCGA-PRAD) where tumour purity is the
  primary confound, in MDS the primary
  confound is MOLECULAR SUBTYPE MIXING in
  the dataset.
  GSE19429 has cytogenetic subtype labels.
  GSE114922 has molecular mutation labels
  (SF3B1, SRSF2, U2AF1).
  The analysis should be run:
  a) On the FULL DATASET (all MDS vs. normal)
     — to replicate the existing bulk analysis.
  b) On SUBTYPE-SPECIFIC SUBSETS (SF3B1-only,
     SRSF2-only, etc.) — to resolve which
     subtypes drive the bulk signal.

NOTE 3 — THE EZH2 DIRECTION TEST:
  The first analytical question in each
  before-document must be:
  "In this specific molecular subtype,
  is EZH2 ELEVATED, REDUCED, or NORMAL
  relative to healthy CD34+ controls?"
  This must be stated as a PRE-DATA PREDICTION
  and then confirmed or corrected by the analysis.
  The prediction from biology:
    SF3B1-mutant: EZH2 likely ELEVATED
    SRSF2-mutant: EZH2 likely REDUCED (mis-spliced)
    TET2-mutant: EZH2 likely ELEVATED RELATIVELY
                 (TET2 loss allows EZH2 to be
                 more dominant)
    Monosomy 7: EZH2 likely REDUCED (haploinsufficient)
    IDH-mutant: EZH2 likely ELEVATED (2-HG drives
                H3K27me3 via KDM inhibition)
  These are the pre-data predictions to be stated
  in the before-documents and confirmed by data.

NOTE 4 — IPSS-M VALIDATION:
  The framework depth score in MDS must be
  compared to IPSS-M risk category for every
  patient in the dataset.
  Expected finding: depth score correlates with
  IPSS-M risk category.
  Going-further finding: depth score provides
  additional prognostic information BEYOND
  IPSS-M in at least one subtype.
  This is the primary clinical validation and
  contribution of the MDS depth score series.
```

---

## SECTION XII — PLANNED ANALYSIS ORDER

```
ORDER:

  MDS-S1   SF3B1-mutant      GSE114922 +          HIGH
           (benign,          GSE63569
           erythroid         SF3B1 subset
           stalling,         + healthy CD34+
           luspatercept      controls
           confirmed)        REASON: The most
                             tractable subtype.
                             Well-characterised
                             erythroid stalling
                             mechanism (ABCB7/
                             TMEM14C mis-splicing).
                             Luspatercept is the
                             confirmed drug —
                             depth score should
                             predict luspatercept
                             response.
                             EZH2 direction test
                             for SF3B1 subtype.
                             Circuit status:
                             is the GATA1 circuit
                             intact (erythroid
                             programme stalled
                             but restoreable)?
                             Switch gene candidate:
                             TMEM14C, ABCB7,
                             or ALAS2?
                             Clinical output:
                             depth score as
                             luspatercept response
                             predictor — which
                             SF3B1-mutant patients
                             are deepest (most
                             stalled) and most
                             likely to benefit.

  MDS-S2   del(5q)           GSE19429             MOD-HIGH
           (ribosomal        del(5q) subset
           stress,           + healthy CD34+
           lenalidomide      controls
           confirmed)        REASON: The second
                             most tractable subtype.
                             Lenalidomide mechanism
                             well-understood
                             (CSNK1A1 synthetic
                             lethality).
                             Depth score should
                             capture RPS14 loss
                             signature (ribosomal
                             stress markers).
                             P53 signature
                             in del(5q) — depth-
                             positive?
                             TP53 mutation emergence
                             depth signal — can
                             the depth score detect
                             the TP53 clone before
                             clinical progression?
                             Clinical output:
                             depth score as TP53-
                             emergence early warning
                             tool in del(5q) patients
                             on lenalidomide.

  MDS-S3   SRSF2-mutant      GSE114922            MOD
           (granulocytic     SRSF2 subset
           bias,             + healthy CD34+
           CMML risk)        REASON: The most
                             complex spliceosome
                             subtype. Granulocytic
                             dysplasia mechanism.
                             EZH2 direction test:
                             is EZH2 REDUCED
                             (SRSF2 mis-splices EZH2)?
                             CEBPα circuit test:
                             is granulocytic
                             differentiation circuit
                             intact or broken?
                             CMML transition signal:
                             does depth score
                             separate SRSF2-MDS
                             from SRSF2-CMML?
                             Clinical output:
                             depth score as CMML
                             progression risk tool
                             in SRSF2-mutant MDS.

  MDS-S4   TET2/ASXL1/IDH    GSE114922            MOD
           epigenetic        TET2/ASXL1 subset
           modifiers         + IDH1/IDH2 subset
                             + healthy CD34+
                             REASON: Epigenetic
                             modifier MDS. EZH2
                             relationship to TET2
                             and ASXL1. IDH-mutant
                             CIMP signature.
                             Vitamin C/TET2 depth
                             prediction: does high
                             dose vitamin C reverse
                             depth signature?
                             IDH inhibitor target:
                             does ivosidenib depth
                             score predict response?
                             Cross-cancer: IDH-mutant
                             MDS vs. IDH-mutant
                             GBM — same CIMP depth
                             signature?
                             Clinical output:
                             EZH2 direction resolved
                             by epigenetic modifier
                             type; ivosidenib
                             selection tool for
                             IDH-mutant MDS.

  MDS-S5   biTP53/complex    GSE19429             MOD
           karyotype         complex karyotype
           (deepest          subset + biTP53
           false attractor)  cases + healthy
                             CD34+ controls
                             REASON: The most
                             advanced MDS state.
                             Circuit broken.
                             APR-246 Phase III
                             failure reinterpreted:
                             which TP53 mutation
                             TYPE (mut/mut vs.
                             mut/del) predicts
                             APR-246 response?
                             EZH2 in complex
                             karyotype (with EZH2
                             deletion on monosomy 7)?
                             MYC as depth-positive
                             driver in +8 cases?
                             Clinical output:
                             APR-246 patient
                             selection depth score;
                             TP53-mutation-type
                             selection criterion.

  MDS-X    Cross-subtype     After S1–S5.
           + IPSS-M          QUESTIONS:
           correlation         1. Does depth score
                                  correlate with
                                  IPSS-M risk category
                                  across all subtypes?
                               2. Which subtype's
                                  depth score is most
                                  predictive of blast
                                  count trajectory?
                               3. EZH2 direction
                                  summary across all
                                  five subtypes:
                                  is the pattern
                                  consistent with
                                  the pre-data
                                  predictions?
                               4. The 2024 Phase III
                                  failure analysis:
                                  ENHANCE, STIMULUS-
                                  MDS2, VERONA, APR-246
                                  Phase III — in each,
                                  which SUBTYPE-DEPTH
                                  combination would
                                  have predicted
                                  response?
                                  This is the unified
                                  lesson for MDS
                                  clinical development:
                                  unselected higher-
                                  risk MDS is four
                                  distinct diseases.
                                  Treating them as one
                                  = four consecutive
                                  trial failures.
```

---

## SECTION XIII — WHAT THIS DOCUMENT DOES NOT CONTAIN

```
This document contains:
  ✓ Normal myeloid differentiation hierarchy
    (HSC → MPP → CMP → GMP/MEP, with full TF
    network: GATA2/RUNX1/PU.1/CEBPA/GATA1/KLF1)
  ✓ The GATA1/PU.1 antagonism at the CMP level
    and its MDS subtype implications
  ✓ WHO 2022 classification (three molecular
    entities + morphological categories + blast
    threshold hierarchy + IPSS-M)
  ✓ SF3B1 subtype: ring sideroblast mechanism
    (ABCB7/TMEM14C mis-splicing, mitochondrial
    iron trapping, stalled erythroid programme),
    luspatercept as secondary brake removal
    (confirmed), spliceosome modulators as
    primary defect targeting (novel)
  ✓ del(5q) subtype: RPS14/ribosomal stress
    mechanism, p53/MDM2 axis, lenalidomide
    CSNK1A1 synthetic lethality (confirmed),
    TP53-emergence under lenalidomide (depth
    score monitoring prediction)
  ✓ SRSF2/U2AF1 subtypes: granulocytic bias,
    EZH2 mis-splicing by SRSF2, CMML overlap
  ✓ Epigenetic modifier mutations (TET2, DNMT3A,
    ASXL1, IDH1/2) and CHIP-to-MDS progression
  ✓ TET2/EZH2 opposition and HMA mechanism
  ✓ IDH-mutant CIMP mechanism and cross-cancer
    connection (GBM, STAD, PRAD IDH subtype)
  ✓ Chromosomal lesions with EZH2 direction
    (monosomy 7 = EZH2 haploinsufficient —
    EZH2 inhibitors would worsen disease)
  ✓ biTP53-MDS as deepest false attractor —
    APR-246 Phase III failure reinterpreted
    (unselected patient population problem)
  ✓ The 2024 Phase III failure cluster:
    ENHANCE, STIMULUS-MDS2, VERONA, APR-246
    Phase III — four consecutive failures from
    the same root cause (unselected population)
  ✓ The depth axis in MDS as differentiation
    execution capacity (vs. identity loss in
    solid tumours)
  ✓ The CD34+ enrichment requirement (mandatory
    — more critical than purity filtering in PAAD)
  ✓ EZH2 direction predictions by subtype
    (to be confirmed or corrected by data)
  ✓ Data availability (GSE19429, GSE114922,
    GSE63569, GSE253355)
  ✓ IPSS-M validation as the primary clinical
    validation criterion for MDS depth score

This document does NOT contain:
  ✗ Depth score predictions
  ✗ Specific switch gene predictions
    (beyond the structural framework analysis
    from the existing bulk MDS analysis)
  ✗ Drug target predictions beyond the
    confirmed drugs (luspatercept, lenalidomide,
    ivosidenib, enasidenib) and the structural
    novel prediction (H3B-8800 + luspatercept
    for SF3B1-MDS)
  ✗ Epigenetic mechanism hypotheses

All of the above belong in the BEFORE documents.
MDS-S1a (SF3B1 before-document) is next.
Written before any script runs.
Before any data loads.
```

---

## STATUS BLOCK

```
document:           MDS_Subtype_Orientation.md
folder:             Cancer_Research/MDS/Subtypes/
status:             COMPLETE — ORIENTATION ONLY
date:               2026-03-03
author:             Eric Robert Lawson
                    OrganismCore

subtypes_covered:
  SF3B1-mutant (~25%):    Erythroid stalling,     [1 of 5]
                          ABCB7/TMEM14C mis-splice,
                          luspatercept confirmed,
                          EZH2 likely elevated
  del(5q) (~8%):          Ribosomal stress,        [2 of 5]
                          RPS14/p53/MDM2 axis,
                          lenalidomide confirmed,
                          TP53-emergence depth signal
  SRSF2/U2AF1 (~15%):     Granulocytic dysplasia,  [3 of 5]
                          EZH2 likely reduced,
                          CMML risk
  TET2/DNMT3A/ASXL1/IDH   CHIP-to-MDS progression, [4 of 5]
  (~40% combined):        EZH2/TET2 opposition,
                          IDH CIMP cross-cancer
  biTP53/complex karyotype:Deepest false attractor, [5 of 5]
  (~5–10% of MDS,          circuit absent,
   ~25% of t-MDS)         APR-246 Phase III reinterp.

key_structural_difference_from_solid_tumours:
  MDS is a FALSE ATTRACTOR OF FALSE ATTRACTORS.
  The HSC/progenitor is in a false attractor
  (first level). Its dysplastic progeny are
  also in false attractors (derived level).
  Depth = differentiation execution capacity
  (not identity loss as in solid tumours).
  Normal reference = healthy donor CD34+ cells
  (not adjacent tissue or GTEx).
  CD34+ enrichment of samples is MANDATORY.

key_pattern_4_application:
  EZH2 DIRECTION IS SUBTYPE-SPECIFIC IN MDS:
  SF3B1-mutant → EZH2 likely ELEVATED (predicted)
  SRSF2-mutant → EZH2 likely REDUCED (predicted)
  Monosomy 7 → EZH2 HAPLOINSUFFICIENT (genomic)
  TET2-mutant → EZH2 relatively ELEVATED (predicted)
  IDH-mutant → EZH2 ELEVATED via 2-HG (predicted)
  biTP53 → EZH2 VARIABLE (complex karyotype dependent)
  EZH2 inhibition is ONLY appropriate where EZH2
  is elevated as a gain-of-function lock.
  In monosomy 7 MDS: EZH2 inhibition = worsening.
  This is Framework Pattern 4 most completely
  illustrated across a single cancer type.

key_2024_clinical_landscape_finding:
  FOUR CONSECUTIVE PHASE III FAILURES IN 2023–2024
  IN HIGHER-RISK MDS:
  ENHANCE (magrolimab + AZA) — negative, halted
  STIMULUS-MDS2 (sabatolimab + AZA) — negative
  VERONA (venetoclax + AZA) — negative
  APR-246 Phase III (eprenetapopt + AZA) — negative
  PLUS: CELLO-1 (tazemetostat in mCRPC) — negative
  All five negative results share the same root cause:
  molecularly rational drugs given to UNSELECTED
  patient populations.
  The depth score framework is the proposed solution
  to all five failures — molecular stratification
  of patient populations before trial enrolment.
  This is the most important clinical lesson that
  emerges from the MDS + PRAD subtype series
  considered together.

analyses_started:   0 (new subtype series)
existing_analysis:  Cancer_Research/MDS/ — complete
                    (OrganismCore cancer series,
                    bulk signal)
next_document:      MDS-S1a
                    SF3B1-mutant Before-Document
```
