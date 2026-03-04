# ACUTE LYMPHOBLASTIC LEUKEMIA — SUBTYPE ORIENTATION DOCUMENT
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

  A complete map of the acute lymphoblastic leukemia
  subtype landscape — what the two primary lineages
  are (B-ALL and T-ALL), why they are biologically
  distinct diseases arising from different cells in
  different anatomical compartments, what the major
  molecular subtypes within each lineage are, what the
  normal lymphoid development hierarchy looks like as
  a Waddington landscape, what public data exists to
  analyze each entity, and what the existing analysis
  in this repository captured.

ALL has a lineage problem that is unique among the
haematological malignancies in this repository:

  B-ALL and T-ALL are not subtypes of the same disease.
  They are two diseases that share a name because both
  produce immature lymphoid blasts on the blood film.
  Beyond that, the similarity ends.

  B-ALL:   Arises in the bone marrow.
           From a B-cell progenitor.
           In a child aged 2–10 years at peak.
           5-year survival in children: >90%.
           Driven by PAX5, EBF1, IKZF1 disruption.

  T-ALL:   Arises in the thymus.
           From a T-cell progenitor.
           In adolescents and young adults.
           5-year survival: ~70-80% in children,
                            ~40-50% in adults.
           Driven by NOTCH1, TAL1, TLX1/3 disruption.

  The existing analysis in Cancer_Research/ALL/ ran on
  a mixed or lineage-unspecified ALL dataset.
  This document establishes what that analysis captured
  and what the subtype analyses must add.

Within B-ALL alone, there are at least 15 molecularly
distinct subtypes — each with a different initiating
genetic event, different cell of arrest in the B-cell
hierarchy, different prognosis, and different drug
sensitivity. They all look alike on the blood film.
They are not alike in any meaningful biological sense.

This is the most molecularly fragmented disease in
the repository. The orientation document must hold
that complexity without collapsing it.
```

---

## DOCUMENT METADATA

```
document_id:        ALL_Subtype_Orientation
series:             ALL (Acute Lymphoblastic Leukemia — Subtypes)
folder:             Cancer_Research/ALL/Subtypes/
date:               2026-03-03
author:             Eric Robert Lawson
status:             ORIENTATION ONLY — no predictions,
                    no analysis, no epigenetic hypotheses
next_document:      ALL_BALL_before.md
                    (Document ALL-S1a — B-ALL before-doc)
protocol_version:   Workflow_Protocol.md v2.0
```

---

## SECTION I — THE NORMAL LYMPHOID DEVELOPMENT HIERARCHY

```
The Waddington baseline for ALL is the normal lymphoid
development hierarchy — two entirely separate programmes
running in two different anatomical compartments.

═══════════════════════════════════════════════════════
THE B-CELL HIERARCHY (bone marrow)
═══════════════════════════════════════════════════════

ALL B-CELL DEVELOPMENT OCCURS IN THE BONE MARROW:

  Haematopoietic Stem Cell (HSC)
  ↓
  Common Lymphoid Progenitor (CLP)
    CD34+, CD10+, CD127+ (IL7R)
    Committed to lymphoid lineage but not yet B or T
  ↓
  Pro-B cell  [EARLIEST B-LINEAGE CELL]
    CD19+, CD34+, CD10+, TdT+
    PAX5 first activates here — marks B-lineage commitment
    EBF1 active — opens chromatin at B-cell loci
    IKZF1 (Ikaros) active — lymphoid commitment
    Immunoglobulin heavy chain (IGH) begins D-J
    rearrangement
    NOT YET: surface immunoglobulin, light chain
  ↓
  Pre-B cell  [PRE-B RECEPTOR CHECKPOINT]
    CD19+, CD34−, CD10+, cytoplasmic μ heavy chain+
    Pre-BCR assembled (μ heavy chain + surrogate
    light chain VpreB/λ5)
    Pre-BCR signals survival and proliferation
    PAX5 high — locks B-identity, represses T/myeloid
    IGH V-DJ rearrangement completes
    Light chain (IGL/IGK) rearrangement begins
  ↓
  Immature B cell
    CD19+, surface IgM+
    Complete BCR assembled
    Self-reactive BCR → central tolerance (deletion/
    receptor editing)
    Not yet fully functional
  ↓
  Mature naive B cell
    CD19+, IgM+, IgD+
    Exits bone marrow to peripheral blood
    Awaits antigen encounter

KEY TRANSCRIPTION FACTORS THAT DEFINE B-CELL IDENTITY:

  IKZF1 (Ikaros):   Lymphoid master TF
                    Required for CLP → pro-B transition
                    Acts as tumour suppressor in B-ALL:
                    IKZF1 deletion/mutation = POOR PROGNOSIS
                    in Ph+ ALL and Ph-like ALL
                    When IKZF1 is lost, the B progenitor
                    loses its lymphoid identity driver —
                    the cell becomes more stem-like and
                    more resistant to treatment

  EBF1:             B-lineage specification TF
                    Opens chromatin at pro-B stage
                    Cooperates with PAX5 to lock B fate
                    Mutations in EBF1 are rare but seen
                    in some B-ALL subtypes

  PAX5:             The "guardian of B-cell identity"
                    Activated at pro-B stage — maintained
                    through all B-lineage stages
                    PAX5 does two things simultaneously:
                      1. Activates B-cell identity genes
                         (CD19, BLNK, CD79a)
                      2. Represses alternative fates
                         (T-cell genes, myeloid genes,
                         NK cell genes)
                    PAX5 MUTATION/DELETION: present in
                    ~30% of B-ALL — allows the cell to
                    de-repress alternative lineage genes
                    This is the B-ALL identity crisis:
                    loss of PAX5 → lineage infidelity
                    PAX5 amplification: paradoxical
                    amplification seen in some subtypes

  E2A (TCF3):       Early B-cell TF — cooperates with
                    EBF1 at pro-B stage
                    TCF3-PBX1 fusion: a B-ALL subtype
                    TCF3-HLF fusion: rare but very poor
                    prognosis B-ALL

  RAG1/RAG2:        V(D)J recombinase — expressed only
                    in pro-B and pre-B cells when IGH/
                    IGL rearrangements are occurring
                    Off-target RAG activity is a source
                    of B-ALL translocations (ETV6-RUNX1,
                    BCR-ABL1, TCF3-PBX1)

THE WADDINGTON STRUCTURE OF NORMAL B-CELL DEVELOPMENT:

  The B-cell hierarchy is a Waddington landscape where:
  - The HSC is the most "undifferentiated" position
  - The mature B cell is the terminal differentiated
    normal attractor
  - Each developmental stage is locked by specific TFs

  B-ALL FALSE ATTRACTOR STRUCTURE:
  Each B-ALL subtype arrests the cell at a specific
  developmental stage — producing a blast cell that
  has the immunophenotype of that stage:

    KMT2A-rearranged:  Arrest near pro-B / CLP
                       (the most primitive arrest)
    ETV6-RUNX1:        Arrest at late pro-B / early pre-B
    BCR-ABL1:          Arrest at pre-B
    TCF3-PBX1:         Arrest at pre-B
    Hyperdiploid:      Arrest at pre-B
    DUX4:              Arrest at pre-B / immature B
    MEF2D:             Arrest at pro-B to pre-B
    Ph-like:           Variable, often pre-B arrest

  The depth score in B-ALL is therefore not primarily
  "how far has this cell moved from normal" —
  it is "how far has this arrested progenitor drifted
  from its developmental arrest position deeper into
  a self-renewing proliferative state."
  The depth axis captures:
    Loss of the arrest-stage identity (e.g., pre-BCR
    signalling genes in pre-B ALL)
    Gain of proliferative/survival genes
    Loss of B-lineage TF activity (PAX5 target genes)
    Gain of stem-like signatures (HOXA genes, FLT3)

══════════════════════════════════════════════════════���
THE T-CELL HIERARCHY (thymus)
═══════════════════════════════════════════════════════

ALL T-CELL DEVELOPMENT OCCURS IN THE THYMUS:

  Haematopoietic Stem Cell (HSC)
  ↓
  T-lineage progenitor migrates from bone marrow
  to THYMUS via blood
  ↓
  Early Thymic Progenitor (ETP) — arrives in thymus
    CD7+, CD5−/low, CD1a−, cCD3+
    Still multipotent — can make T, NK, myeloid cells
    This is the cell of origin for ETP-ALL
  ↓
  Double Negative (DN) thymocytes
    DN2: CD44+, CD25+, CD117+, KIT+
         T-lineage commitment beginning
    DN3: CD44−, CD25+
         TCR β-chain rearrangement — the β-selection
         checkpoint (equivalent to pre-BCR checkpoint)
    DN4: CD44−, CD25−
         Passed β-selection. Proliferating.
  ↓
  Double Positive (DP) thymocyte  [THE CORTICAL STAGE]
    CD4+, CD8+ simultaneously
    TCR α-chain rearranges
    Positive selection: TCR must recognise self-MHC
    Negative selection: TCR must not react strongly
    to self-antigens
    This is the cell of origin for TLX1/TLX3-ALL
  ↓
  Single Positive (SP) thymocyte
    CD4+ or CD8+ (lineage commitment)
    Migrates to peripheral blood
    Mature naive T cell

KEY TRANSCRIPTION FACTORS THAT DEFINE T-CELL IDENTITY:

  NOTCH1:           THE master TF of T-ALL
                    ~60% of all T-ALL has activating
                    NOTCH1 mutation regardless of subtype
                    NOTCH1 is the secondary driver that
                    amplifies the proliferative signal
                    from whatever the primary subtype
                    driver is (TLX1, TAL1, HOXA, etc.)
                    NOTCH1 inhibitors (γ-secretase
                    inhibitors) have been studied —
                    limited clinical success so far due
                    to gut toxicity

  GATA3:            T-lineage identity TF — cooperates
                    with NOTCH1 at DN2/DN3 stage

  BCL11B:           T-lineage commitment TF
                    Loss of BCL11B = ETP phenotype
                    (most primitive T-ALL)

  TCF7 (TCF1):      T-cell identity and quiescence

  TAL1 (SCL):       Normally expressed in HSC and
                    erythroid cells — repressed in
                    mature T cells
                    TAL1 reactivation in T-ALL = the
                    cell de-represses a primitive
                    haematopoietic identity

  CDKN2A/B:         Deleted in >70% of T-ALL
                    (higher rate than any other ALL
                    subtype) — ubiquitous cell cycle
                    brake removal

THE WADDINGTON STRUCTURE OF NORMAL T-CELL DEVELOPMENT:

  The T-cell hierarchy is a two-organ Waddington
  landscape (bone marrow → thymus). The critical
  difference from B-cell development:

  T-CELL DEVELOPMENT IS EXTERNAL TO THE BONE MARROW.
  It requires migration from marrow to thymus,
  physical interactions with thymic stroma, and
  two checkpoints (β-selection and positive/negative
  selection) that have no marrow equivalent.

  T-ALL FALSE ATTRACTOR STRUCTURE:
  Like B-ALL, each T-ALL subtype arrests at a specific
  stage — but the developmental positions are defined
  by thymic development, not marrow development:

    ETP-ALL:         Arrest at earliest thymic
                     progenitor (before DN2)
                     Still multipotent — can express
                     myeloid markers
    HOXA-ALL:        Arrest at DN2-DN3
    TLX1/TLX3-ALL:   Arrest at DP cortical stage
    TAL1-ALL:        Arrest at late DP/SP stage

  THE KEY IMPLICATION FOR DEPTH SCORING:
  The "normal" reference for T-ALL is not blood or
  bone marrow — it is the THYMUS.
  The depth score axis runs from the arrested thymic
  stage toward a self-renewing primitive state that
  has lost thymic identity.
  This requires a thymic normal reference dataset —
  which is less available than bone marrow reference
  datasets. This is addressed in the data section.
```

---

## SECTION II — THE B-ALL MOLECULAR SUBTYPES

```
B-ALL is the most molecularly heterogeneous malignancy
in this repository. At least 15 discrete molecular
subtypes are defined by WHO 2022 / ICC 2022.

They are grouped here by MECHANISM and PROGNOSIS
rather than alphabetically, to make the Waddington
structure legible.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
GROUP 1 — FAVOURABLE PROGNOSIS SUBTYPES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

SUBTYPE 1A — HYPERDIPLOIDY (HIGH HYPERDIPLOID)
  Frequency:        25% of childhood B-ALL
                    (the most common single subtype
                    in children)
                    5% of adult B-ALL
  Chromosomes:      >50 chromosomes (51–67)
                    Specific trisomies: +4, +10, +17
                    are the most favorable combination
  Prognosis:        EXCELLENT — 5-year EFS >85%
                    in children
  Initiating event: Unknown — hyperdiploidy is the
                    founding event itself.
                    Thought to arise from a single
                    abnormal mitosis producing a
                    near-tetraploid cell that then
                    loses chromosomes to hyperdiploid
                    range.
  Why it works:     Extra chromosome copies = gene
                    dosage effects for survival genes
                    AND for methotrexate sensitivity
                    (extra chromosome 21 = extra DHFR
                    copies → paradoxically more
                    methotrexate polyglutamate
                    accumulation)
  Arrest stage:     Pre-B cell
  Normal identity:  Pre-BCR signalling retained

SUBTYPE 1B — ETV6-RUNX1 (TEL-AML1, t(12;21))
  Frequency:        25% of childhood B-ALL
                    (together with hyperdiploid,
                    these two account for ~50% of
                    all pediatric B-ALL)
                    <1% of adult B-ALL
  Translocation:    t(12;21)(p13;q22)
                    Fuses ETV6 (TEL) with RUNX1 (AML1)
                    ETV6 is a transcriptional repressor
                    RUNX1 is the master haematopoietic
                    TF of the CLP stage
                    The fusion creates a dominant
                    negative repressor of RUNX1 targets
  Prognosis:        EXCELLENT — 5-year OS >95%
                    Long-term cure rates approaching
                    those of normal children
  Key feature:      ETV6-RUNX1 B-ALL can RELAPSE very
                    late — years after apparent cure.
                    The founding clone may persist in
                    a quiescent state (possible neonatal
                    origin) and re-emerge.
  Arrest stage:     Late pro-B / early pre-B
  Structural note:  ETV6-RUNX1 is found in neonatal
                    blood spots in children who later
                    develop ALL — the translocation
                    arises in utero. A second hit is
                    required for frank ALL.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
GROUP 2 — INTERMEDIATE PROGNOSIS SUBTYPES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

SUBTYPE 2A — TCF3-PBX1 (E2A-PBX1, t(1;19))
  Frequency:        5–6% of childhood B-ALL
  Translocation:    t(1;19)(q23;p13)
                    Fuses TCF3 (E2A) with PBX1
                    PBX1 is normally repressed in
                    B-cell progenitors — the fusion
                    reactivates a HOXA-like programme
  Prognosis:        INTERMEDIATE — initially poor,
                    outcome has improved dramatically
                    with intensified therapy (now ~80%
                    EFS with contemporary regimens)
  Arrest stage:     Pre-B cell

SUBTYPE 2B — DUX4-REARRANGED
  Frequency:        5–7% of B-ALL in children
                    Higher in AYA (7–10%)
  Genetics:         IGH-DUX4 or ERG deletion (del4q28)
                    DUX4 is normally expressed only in
                    2-cell embryos (totipotent stage)
                    Its reactivation in a B-cell
                    progenitor creates a profound
                    identity conflict
  Prognosis:        FAVOURABLE-INTERMEDIATE
                    Surprisingly good response to
                    standard therapy despite unusual
                    biology
  Special feature:  DUX4-rearranged B-ALL frequently
                    has an unusual immunophenotype with
                    CD2 aberrant expression (T-cell
                    marker on a B-cell tumour)

SUBTYPE 2C — MEF2D-REARRANGED
  Frequency:        3–5% of B-ALL
                    Higher in adults than children
  Genetics:         MEF2D fusion with various partners
                    (BCL9, HNRNPUL1, SS18, DAZAP1)
                    MEF2D is a differentiation TF —
                    its rearrangement disrupts normal
                    B-cell maturation signals
  Prognosis:        INTERMEDIATE-POOR
                    Higher relapse rate than standard
                    B-ALL; not yet clearly targetable

SUBTYPE 2D — ZNF384-REARRANGED
  Frequency:        3–5% of B-ALL
  Genetics:         ZNF384 fused with TCF3, EP300,
                    CREBBP, TAF15, ARID1B
                    ZNF384 is a chromatin-associated
                    zinc finger TF
  Prognosis:        VARIABLE — depends on fusion
                    partner
  Note:             ZNF384-rearranged ALL has a
                    mixed lineage appearance — often
                    co-expresses myeloid markers
                    (CD13, CD33) alongside B-ALL
                    markers. This lineage ambiguity
                    is a Waddington identity feature:
                    the false attractor has not fully
                    committed to B identity.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
GROUP 3 — HIGH-RISK SUBTYPES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

SUBTYPE 3A — BCR-ABL1 (PHILADELPHIA CHROMOSOME, t(9;22))
  Frequency:        3–5% of childhood B-ALL
                    25–30% of adult B-ALL
                    Increases dramatically with age —
                    the most common subtype in adults
  Translocation:    t(9;22)(q34;q11)
                    Fuses BCR with ABL1
                    Creates a constitutively active
                    ABL1 tyrosine kinase
  Prognosis:        HISTORICALLY POOR — now transformed
                    by tyrosine kinase inhibitors (TKIs)
                    Imatinib + chemotherapy: ~70% EFS
                    in children (vs ~25% historically)
                    Dasatinib or ponatinib + chemo:
                    standard; stem cell transplant
                    increasingly avoidable in good
                    responders
  Cell of origin:   Likely HSC or CLP — the BCR-ABL1
                    fusion can drive either CML (myeloid
                    dominant) or ALL (lymphoid dominant)
                    depending on the cell in which it
                    arises. In ALL, the cell arrested
                    at the pre-B stage.
  IKZF1 deletion:   ~70% of BCR-ABL1 B-ALL
                    IKZF1 deletion is the key cooperating
                    event that drives myeloid to lymphoid
                    dominance and confers
                    treatment resistance.

SUBTYPE 3B — Ph-LIKE (BCR-ABL1-LIKE)
  Frequency:        10–15% of standard-risk children
                    25–30% of adolescents and adults
                    (the prevalence RISES with age —
                    opposite of hyperdiploid and
                    ETV6-RUNX1 which fall with age)
  Definition:       Has the gene expression profile of
                    BCR-ABL1 B-ALL but LACKS the
                    BCR-ABL1 fusion itself.
                    More than 70 different kinase-
                    activating alterations have been
                    identified:
  MECHANISM GROUPS:
    JAK/STAT group (~50% of Ph-like):
      CRLF2 rearrangement (to P2RY8 or IGH)
        + JAK2 mutation (~60% of CRLF2-r cases)
      EPOR rearrangement + JAK2
      JAK2 fusion (ETV6-JAK2, BCR-JAK2)
      Targeted by: ruxolitinib (JAK inhibitor)
                   baricitinib
    ABL-class group (~25% of Ph-like):
      ABL1 fusions (non-BCR-ABL1)
      ABL2 fusions
      PDGFRB fusions
      CSF1R fusions
      Targeted by: dasatinib, imatinib, ponatinib
    Other kinase group (~25%):
      FLT3 mutations
      NTRK fusions
      RET fusions
      FGFR fusions
  Prognosis:        POOR without targeted therapy
                    5-year OS as low as ~23% in historical
                    cohorts treated with chemo alone
                    TKI addition improving MRD clearance
                    Ruxolitinib + chemo: COG AALL1521
                    trial showing improved response
  IKZF1 deletion:   ~70% of Ph-like ALL
                    Same cooperating event as BCR-ABL1.
                    The IKZF1 deletion is the common
                    thread that makes both BCR-ABL1 and
                    Ph-like resistant to standard therapy.

SUBTYPE 3C — KMT2A-REARRANGED (MLL-REARRANGED)
  Frequency:        5% of childhood B-ALL overall
                    ~80% of INFANT B-ALL (<1 year old)
                    The dominant leukemia of infancy
  Genetics:         KMT2A (MLL) gene at 11q23 fused
                    with >80 different partner genes
                    Most common in B-ALL:
                      KMT2A-AFF1 (t(4;11)): infant ALL
                      KMT2A-MLLT1 (t(11;19))
                      KMT2A-MLLT3 (t(9;11))
                    KMT2A is a histone methyltransferase
                    (H3K4 methylation — active chromatin)
                    KMT2A fusion = HOXA cluster genes
                    constitutively activated
                    HOXA activation = the cell reverts
                    to a very primitive haematopoietic
                    identity
  Prognosis:        VERY POOR in infants (<1 year)
                    5-year OS ~25–35% in infant ALL
                    Better in older children (~60%)
  Cell of origin:   The most primitive arrest in B-ALL —
                    pro-B or even CLP. In infants, the
                    KMT2A rearrangement likely arises
                    in a fetal HSC or CLP.
                    This is in utero origin.
  Note:             KMT2A-AFF1 infant ALL also expresses
                    myeloid markers — the most lineage-
                    ambiguous B-ALL subtype.
                    Menin inhibitors (revumenib, ziftomenib)
                    target the KMT2A/menin interaction —
                    now in clinical trials and showing
                    early efficacy in KMT2A-r ALL.

SUBTYPE 3D — HYPODIPLOIDY (LOW HYPODIPLOIDY)
  Frequency:        ~2% of B-ALL
  Chromosomes:      <44 chromosomes
                    Near-haploid (24–31 chromosomes):
                    very rare, dismal outcome
                    Low hypodiploid (32–39 chromosomes):
                    slightly more common, still very poor
  TP53 mutation:    ~75% of low hypodiploid ALL
                    AND ~50% of these TP53 mutations
                    are GERMLINE (inherited TP53 —
                    Li-Fraumeni syndrome).
                    Hypodiploidy is the only B-ALL
                    subtype with a significant germline
                    cancer predisposition component.
  Prognosis:        WORST of all B-ALL subtypes
                    5-year OS ~25% for near-haploid
                    ~30–40% for low hypodiploid
  NOTE:             Low hypodiploid B-ALL frequently
                    doubles its chromosomes to produce
                    a near-tetraploid clone that mimics
                    hyperdiploid morphologically —
                    but the outcome is very different.
                    Distinguishing doubled-hypodiploid
                    from true hyperdiploid requires
                    SNP array (not just chromosome count).

SUBTYPE 3E — iAMP21 (INTRACHROMOSOMAL AMPLIFICATION
              OF CHROMOSOME 21)
  Frequency:        2% of B-ALL
                    Enriched in older children (9–14 yr)
  Genetics:         ≥5 copies of the RUNX1 locus on a
                    single structurally abnormal chr21
                    Mechanism: chromothripsis of chr21
  Prognosis:        POOR with standard therapy
                    GOOD with intensified therapy
                    (iAMP21 is treatable if identified)
  Note:             iAMP21 is one of the reasons
                    molecular diagnostics matter in
                    ALL — the STANDARD therapy
                    is INSUFFICIENT for this subtype,
                    but the cells look identical to
                    low-risk B-ALL on the blood film.
```

---

## SECTION III — THE T-ALL MOLECULAR SUBTYPES

```
T-ALL is defined by the stage of thymic arrest and
the TF that is aberrantly activated at that stage.
Unlike B-ALL (which has discrete translocations),
T-ALL is more often driven by overexpression of
developmental TFs through enhancer hijacking.

NOTCH1 MUTATION: THE UNIVERSAL BACKGROUND

  NOTCH1 activating mutation is present in ~60% of
  ALL T-ALL subtypes regardless of the primary driver.
  It is a secondary/cooperating event that amplifies
  proliferation — not the initiating event.
  It is the T-ALL equivalent of IKZF1 deletion in
  B-ALL: a near-universal cooperating event that
  worsens prognosis and drives resistance.

CDKN2A/B DELETION: THE NEAR-UNIVERSAL CELL CYCLE EVENT

  >70% of T-ALL has biallelic CDKN2A/B deletion —
  the highest rate of any ALL subtype.
  This removes both p16 (CDK4/6 inhibition) and
  p14ARF (p53 activation) simultaneously.
  It is essentially universal in T-ALL.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
THE FOUR MAJOR T-ALL SUBTYPES
(by developmental arrest position)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

SUBTYPE A — ETP-ALL (EARLY T-CELL PRECURSOR ALL)
  Frequency:        10–15% of T-ALL
  Immunophenotype:  CD7+, CD2+/−, CD1a−, CD8−
                    CD5 dim/negative
                    Aberrant myeloid markers: CD13+,
                    CD33+, CD34+, CD117+ (KIT)
                    This immunophenotype is diagnostic:
                    a T-cell marker (CD7) on a cell that
                    also expresses myeloid markers
  Cell of origin:   ETP — the earliest thymic progenitor
                    Has not yet committed to T-cell fate
                    Still retains multi-lineage potential
  Genetics:         FLT3 mutation: ~30% (highest of any
                    T-ALL subtype) — FLT3 is myeloid
                    DNMT3A mutation: ~10%
                    RUNX1 mutation: ~10%
                    IDH1/2 mutation: ~5%
                    NOTCH1 mutation: ~30%
                    (lower than other T-ALL subtypes —
                    reflects primitive origin)
                    KMT2A rearrangements: present
  Prognosis:        HISTORICALLY POOR
                    With intensified/risk-adapted therapy:
                    improving (~50–70% EFS in children)
                    FLT3 inhibitors (sorafenib, gilteritinib)
                    being studied — ETP-ALL with FLT3
                    mutation may benefit
  Waddington note:  ETP-ALL has a mixed
                    lymphoid/myeloid identity — the most
                    primitive false attractor in ALL.
                    The false attractor is not "B-cell
                    arrested" or "T-cell arrested" but
                    "lymphoid/myeloid progenitor arrested."
                    This is structurally analogous to
                    KMT2A infant B-ALL — both are arrests
                    at the most primitive haematopoietic
                    progenitor stage.

SUBTYPE B — HOXA-DEREGULATED T-ALL
  Frequency:        ~20% of T-ALL
  Genetics:         KMT2A rearrangements in T-ALL
                    MLLT10 (AF10) fusions
                    PICALM-MLLT10 (t(10;11))
                    SET-NUP214
                    All result in HOXA cluster
                    upregulation (HOXA5, HOXA9, HOXA10)
  Cell of origin:   DN2–DN3 thymocyte
                    (after ETP but before DP stage)
  Prognosis:        INTERMEDIATE — PICALM-MLLT10 is
                    particularly aggressive
  Connection to B-ALL:
                    KMT2A rearrangements drive HOXA
                    upregulation in BOTH KMT2A B-ALL
                    and HOXA T-ALL — the same molecular
                    engine operates in both lineages
                    through the same epigenetic mechanism
                    (H3K4 methylation of HOXA loci).
                    This is a cross-subtype structural
                    connection noted here for the
                    before-document to address.

SUBTYPE C — TLX1 AND TLX3 T-ALL
  Frequency:        TLX1 (HOX11): ~30% of T-ALL
                    TLX3 (HOX11L2): ~20% of T-ALL
                    Together, the most common T-ALL
                    subtypes in adults
  Genetics:         TLX1 overexpression:
                      t(7;10)(q35;q24) — most common
                      TLX1 is a homeodomain TF normally
                      expressed only in spleen/thymus
                      during embryonic development —
                      silenced in mature lymphocytes
                      Its reactivation arrests the cell
                      at the cortical (DP) stage
                    TLX3 overexpression:
                      t(5;14)(q35;q32) — most common
                      Similar TF family to TLX1 but
                      different chromosomal mechanism
  Cell of origin:   Cortical thymocyte (DP stage)
                    The most mature arrest in T-ALL
  Prognosis:        TLX1: FAVOURABLE — best prognosis
                    of all T-ALL subtypes
                    TLX3: INTERMEDIATE-POOR — worse
                    than TLX1 despite TF family
                    similarity
  Normal identity:  DP thymocyte markers: CD1a, CD4,
                    CD8, CD3 — retained in false attractor
                    This is the T-ALL subtype where the
                    normal identity is most preserved —
                    analogous to ETV6-RUNX1 in B-ALL.

SUBTYPE D — TAL1-REARRANGED T-ALL
  Frequency:        25–30% of T-ALL
  Genetics:         STIL-TAL1 microdeletion (~20% T-ALL)
                    TAL1 enhancer hijacking
                    TAL1 is normally expressed in HSCs
                    and erythroid cells — it is
                    COMPLETELY REPRESSED in mature
                    T cells. Its reactivation in a
                    T-cell progenitor produces profound
                    identity dysregulation.
                    TAL1 forms a complex with RUNX1,
                    GATA3, and LMO1/2 — this complex
                    drives T-ALL proliferation
  Cell of origin:   Late cortical or SP thymocyte
                    (the most mature T-ALL arrest)
  Prognosis:        INTERMEDIATE — PTEN loss is common
                    in TAL1 T-ALL and confers resistance
  Additional:       PTEN deletion: ~25% of TAL1 T-ALL
                    (PI3K/AKT activation — relevant
                    drug target)
```

---

## SECTION IV — THE CLINICAL FRAMEWORK

```
THE FUNDAMENTAL PROGNOSTIC SPLIT IN ALL:

  PAEDIATRIC vs ADULT DISEASE:

  Pediatric ALL (age <18):
    5-year OS: >90% overall
               >95% for standard risk
               ~80% for high risk
    Cure is the realistic expectation for most
    children with ALL in high-income countries.
    This is the greatest success story in
    oncology — ALL was uniformly fatal in 1960.
    Now >90% of children are cured.

  Adult ALL (age >18):
    5-year OS: ~40–50% for AYA (18–40)
                ~20–30% for adults >40
                ~15–25% for adults >60
    The prognosis deteriorates sharply with age.
    Reason: not just more aggressive disease —
    adults have:
      Higher frequency of Ph+ and Ph-like ALL
      Higher frequency of IKZF1 deletion
      Less ability to tolerate intensive chemo
      Less frequent hyperdiploidy and ETV6-RUNX1
      (the favourable childhood subtypes)

  This age-biology interaction is structurally
  important for the framework:
  The depth score in B-ALL should correlate with
  the adult-enriched subtypes (Ph+, Ph-like,
  hypodiploid) as the deepest positions in the
  false attractor landscape.

TREATMENT STRUCTURE (brief — for orientation):

  ALL treatment has three phases regardless of
  subtype:

  INDUCTION (weeks 1–5):
    Goal: morphological remission (<5% blasts)
    Drugs: vincristine, dexamethasone/prednisone,
           L-asparaginase, ± anthracycline
    MRD-negativity at end of induction is the
    strongest independent predictor of outcome.
    MRD+ after induction = treatment escalation.

  CONSOLIDATION (months 2–7):
    Goal: eradicate MRD, prevent CNS relapse
    Drugs: high-dose methotrexate, 6-mercaptopurine,
           cytarabine ± intrathecal chemotherapy
    CNS prophylaxis: IT chemotherapy (replacing
    cranial radiation in modern protocols)

  MAINTENANCE (months 8–30):
    Goal: prevent late relapse
    Drugs: oral 6-MP daily + weekly methotrexate
    ± monthly vincristine/steroid pulses
    Duration: 2–3 years total from diagnosis

  TARGETED ADDITIONS BY SUBTYPE:
    Ph+ ALL:        TKI throughout (dasatinib
                    or ponatinib)
    Ph-like ALL:    TKI if ABL-class; ruxolitinib
                    if JAK-class (trial protocols)
    KMT2A infant:   Menin inhibitors (trials)
    ETP-ALL/FLT3:   FLT3 inhibitors (trials)
    Relapsed/MRD+:  Blinatumomab (CD19xCD3 BiTE)
                    Inotuzumab ozogamicin (anti-CD22)
                    CAR-T cell therapy (tisagenlecleucel
                    for paediatric/young adult r/r B-ALL)

THE MRD AXIS:

  Minimal residual disease (MRD) detection by flow
  cytometry or PCR/NGS is the most powerful clinical
  biomarker in ALL — more powerful than any single
  genetic feature for predicting relapse.

  MRD negative at end of induction: low relapse risk
  MRD positive at end of induction: high relapse risk,
  escalate therapy

  THE DEPTH SCORE CONNECTION:
  MRD persistence reflects the cell's resistance to
  chemotherapy-induced apoptosis.
  The false attractor depth score should predict MRD
  persistence: deeper false attractor position =
  more chemotherapy-resistant = more likely to be
  MRD-positive at end of induction.
  This is the highest-value clinical output the ALL
  depth score could produce:
  A pre-treatment depth score that predicts MRD
  outcome and guides early therapy intensification.
  This hypothesis is stated here as a structural
  observation — it will be formalised in the
  before-documents.
```

---

## SECTION V — THE EXISTING ALL ANALYSIS — CONTEXT

```
The existing analysis in Cancer_Research/ALL/ ran on
an ALL dataset — the specific lineage composition
(B-ALL only, T-ALL only, or mixed) determines what
the depth axis found represents.

IF THE EXISTING ANALYSIS WAS PREDOMINANTLY B-ALL:
  The depth axis captured the B-ALL false attractor —
  likely the WNT/BCR-ABL1/IKZF1 resistance axis.
  The most likely depth-positive genes:
    MKI67, TOP2A (proliferation — universal)
    HOXA9, HOXA10 (KMT2A programme — deepest B-ALL)
    FLT3 (primitive progenitor marker)
    EZH2 (predicted epigenetic lock)
  The most likely depth-negative genes:
    PAX5 targets (CD19, BLNK, CD79a)
    Pre-BCR signalling genes
    CDX2 equivalent = PAX5 itself

IF THE EXISTING ANALYSIS WAS MIXED B-ALL + T-ALL:
  The depth axis captured the common variance between
  both lineages — likely the proliferation programme
  and loss of lymphoid identity markers.
  The two lineage-specific signals (PAX5 for B-ALL,
  NOTCH1 targets for T-ALL) would have partially
  cancelled, leaving only the universal:
    High MKI67, TOP2A
    High EZH2
    Low differentiation markers

THE SUBTYPE ANALYSES WILL:
  1. Separate B-ALL and T-ALL completely
  2. Analyze major B-ALL subtypes (Ph+/Ph-like,
     KMT2A, hyperdiploid, ETV6-RUNX1) individually
     or in appropriately grouped cohorts
  3. Analyze T-ALL subtypes (ETP, TLX1/3, TAL1/HOXA)
  4. Compare across subtypes for universal ALL markers
```

---

## SECTION VI — A NOTE ON ANALYTICAL DESIGN FOR ALL

```
ALL presents a unique analytical challenge for the
depth score framework that no other cancer in this
repository shares:

THE CHALLENGE:

  In solid tumours (RCC, BRCA, GBM, CRC):
    The normal cell is the same tissue cell
    adjacent to the tumour.
    The tumour is derived from that cell.
    The depth axis = tumour vs. normal tissue.

  In ALL (and AML, CML, CLL, MM):
    There is no "tumour adjacent normal."
    The normal reference is a specific developmental
    stage cell (pro-B, pre-B, thymocyte) that is
    not routinely biopsied alongside the leukemia.
    The depth axis must be defined as:
    "how far has the blast cell moved from the
    normal cell at its developmental arrest stage"
    OR
    "how far has the blast cell moved toward a
    primitive, self-renewing, treatment-resistant
    haematopoietic progenitor state"

  THESE ARE NOT THE SAME AXIS.
  The first (distance from arrested stage) is the
  within-subtype depth axis.
  The second (distance toward primitive progenitor)
  is the cross-subtype resistance axis.

  RESOLUTION FOR THIS FRAMEWORK:
  The framework will use BOTH axes:

    AXIS 1 (within-subtype):
      Normal reference = sorted normal B-cell or
      T-cell developmental stage cells (available
      in GEO from flow-sorted bone marrow studies)
      Depth score = distance from the arrested stage
      This axis tells you how far a given B-ALL
      blast has de-differentiated beyond its
      normal arrest position.

    AXIS 2 (cross-subtype resistance axis):
      Reference = the most primitive HSC-like state
      Depth score = how stem-like is this blast?
      This axis predicts MRD and relapse — the
      more stem-like the blast, the more resistant
      it is to chemotherapy.

  The before-document for each ALL subtype will
  specify which axis is the primary one and lock
  that decision before any data loads.

THE NORMAL REFERENCE DATASETS:

  Flow-sorted B-cell stages:
    GSE14186:     Sorted human B-cell progenitors
                  (pro-B, pre-B, immature B)
                  n=small but stage-pure
                  CRITICAL for axis 1 in B-ALL
    GSE42519:     Human bone marrow B-cell stages
                  sorted by flow cytometry
    GSE17054:     HSC through mature B-cell RNA-seq
                  Multiple sorted stages

  Flow-sorted T-cell thymic stages:
    GSE28703:     Sorted human thymic subsets
                  (ETP, DN, DP, SP stages)
                  CRITICAL for axis 1 in T-ALL
    GSE74738:     Human thymus single-cell

  Normal bone marrow (mixed):
    GSE120795:    n=11 normal bone marrow
    GSE24759:     Complete haematopoietic hierarchy
                  sorted populations — the gold
                  standard normal reference for
                  haematological malignancies
```

---

## SECTION VII — DATA AVAILABILITY SUMMARY

```
Entity          Primary Dataset     n (samples)    Subtype labels  Power
──────────────────────────────────────────────────────────────────────────
B-ALL (all)     TARGET ALL          ~1,978 total   YES             HIGH
                (phs000218)         (B + T mixed)
                GEO GSE49032        n=221 B-ALL     YES (translo)  MOD
                GSE47051            n=73 B-ALL      YES            LOW
                St Jude GSE60926    n=145 B-ALL     YES            MOD
Ph+/Ph-like     GSE49032 subset     ~50-80          YES            LOW-MOD
                COG AALL0434        via TARGET      YES            MOD
KMT2A infant    TARGET infant ALL   ~50-80          YES            LOW-MOD
Hyperdiploid    TARGET ALL subset   ~400            YES            HIGH
ETV6-RUNX1      TARGET ALL subset   ~400            YES            HIGH
T-ALL           TARGET ALL          ~250 T-ALL      YES            MOD
                St Jude GSE60927    n=72 T-ALL      YES            LOW
ETP-ALL         GSE28703 + clinical ~30-50          YES            LOW
                (limited data)
Normal B-cell   GSE24759            sorted stages   —              REFERENCE
stages          GSE14186            sorted stages   —              REFERENCE
Normal thymic   GSE28703            sorted stages   —              REFERENCE
stages

KEY NOTE ON TARGET DATABASE:
  The NCI TARGET programme (Therapeutically Applicable
  Research to Generate Effective Treatments) is the
  primary public resource for paediatric ALL.
  It contains the largest collection of uniformly
  processed paediatric B-ALL and T-ALL RNA-seq with
  molecular subtype annotation.
  Access: https://ocg.cancer.gov/programs/target
  GDC portal: phs000218 (ALL)
  UCSC Xena: TARGET ALL phase 2 cohort

KEY NOTE ON SAMPLE SIZES:
  The small n problem is most acute for:
    Hypodiploidy: ~30-40 samples in TARGET
    iAMP21: ~20-30 samples
    MEF2D: ~30-50 samples
    ETP-ALL: ~30-50 samples
  These may require cross-cohort pooling or will
  produce LOW-POWER analyses. The before-documents
  must flag this and adjust predictions accordingly.
  The N_MIN threshold from the framework should be
  applied — subtypes with n<15 are UNTESTABLE.
```

---

## SECTION VIII — PLANNED ANALYSIS ORDER

```
ALL has too many subtypes to analyze individually
with the current data availability. The framework
will group them by analytical feasibility and
clinical priority:

ORDER:

  ALL-S1   B-ALL combined    TARGET ALL B-ALL    HIGH POWER
           (all subtypes)    All subtypes pooled
                             REASON: Establish the
                             universal B-ALL depth axis
                             first. What is the common
                             false attractor programme
                             across ALL B-ALL? This
                             is the equivalent of the
                             combined TCGA analysis in
                             solid tumours — the first
                             pass that establishes the
                             overall landscape before
                             decomposing by subtype.

  ALL-S2   Ph+/Ph-like       TARGET Ph+/Ph-like   MODERATE
           B-ALL             combined
                             REASON: The most clinically
                             relevant high-risk adult
                             B-ALL group. TKI and JAK
                             inhibitor targeting.
                             IKZF1 deletion as the
                             depth amplifier.
                             The depth score here
                             predicts MRD persistence
                             and TKI resistance.

  ALL-S3   KMT2A infant      TARGET infant ALL    LOW-MOD
           B-ALL             REASON: The most primitive
                             B-ALL arrest. HOXA
                             programme. Menin inhibitor
                             connection. Lowest survival.
                             Most analogous to ETP-ALL
                             in T-cell lineage.

  ALL-S4   T-ALL combined    TARGET ALL T-ALL     MODERATE
                             All subtypes pooled
                             REASON: Establish the
                             universal T-ALL depth axis.
                             NOTCH1 and CDKN2A/B as
                             near-universal events.
                             TLX1 vs ETP as the
                             differentiation axis.

  ALL-S5   ETP-ALL           TARGET + GEO pooled  LOW
           (if n sufficient) REASON: Most primitive
                             T-ALL. FLT3-high. Most
                             connected to AML landscape.
                             If n<15 → UNTESTABLE.

  ALL-X    Cross-lineage     After S1–S4 complete
           B-ALL vs T-ALL   Questions:
                              1. Universal ALL depth
                                 markers across both
                                 lineages?
                              2. EZH2 in both?
                              3. HOXA genes: depth-
                                 positive in KMT2A
                                 AND ETP-ALL?
                                 (same primitive
                                 progenitor programme?)
                              4. MRD prediction:
                                 does depth score at
                                 diagnosis predict
                                 MRD positivity at
                                 end of induction?
```

---

## SECTION IX — WHAT THIS DOCUMENT DOES NOT CONTAIN

```
This document contains:
  ✓ B-ALL and T-ALL as separate diseases
  ✓ Normal B-cell and T-cell developmental hierarchy
  ✓ All major B-ALL subtypes (15+ grouped by mechanism)
  ✓ All major T-ALL subtypes (4 grouped by arrest stage)
  ✓ Clinical facts (survival, standard of care, MRD)
  ✓ The analytical design challenge (two depth axes)
  ✓ The IKZF1/NOTCH1 cooperating event structure
  ✓ Data availability (TARGET, GEO, St Jude)
  ✓ The MRD connection as a structural hypothesis
  ✓ The HOXA cross-lineage structural connection
    (KMT2A B-ALL ↔ HOXA/ETP T-ALL — same mechanism)
  ✓ Context for the existing ALL analysis

This document does NOT contain:
  ✗ Depth score predictions
  ✗ Specific false attractor gene predictions
  ✗ Drug target predictions
  ✗ Epigenetic lock mechanism hypotheses
  ✗ Cross-subtype structural predictions beyond
    what is already established in the literature

All of the above belong in the BEFORE documents.
ALL-S1a (B-ALL combined before-document) is next.
Written before any script runs. Before any data loads.
```

---

## STATUS BLOCK

```
document:           ALL_Subtype_Orientation.md
folder:             Cancer_Research/ALL/Subtypes/
status:             COMPLETE — ORIENTATION ONLY
date:               2026-03-03
author:             Eric Robert Lawson
                    OrganismCore

entities_covered:
  B-ALL subtypes:   Hyperdiploidy, ETV6-RUNX1,       [9 of ~15]
                    BCR-ABL1, Ph-like, KMT2A,
                    DUX4, MEF2D, ZNF384,
                    TCF3-PBX1, Hypodiploidy, iAMP21
  T-ALL subtypes:   ETP, HOXA, TLX1/3, TAL1          [4 of 4]

analyses_started:   0

existing_analysis:  Cancer_Research/ALL/ — lineage
                    composition of the analyzed dataset
                    determines interpretation.
                    If B-ALL dominant: captured the
                    WNT/BCR-ABL1 resistance axis.
                    If mixed B+T: captured universal
                    proliferation/lymphoid identity axis.
                    Decomposed by analyses in this folder.

next_document:      ALL-S1a
                    B-ALL Combined Before-Document
                    (all B-ALL subtypes pooled from
                    TARGET; predictions locked before
                    any data loads)

critical_note_1:    B-ALL and T-ALL are not subtypes
                    of the same disease. They arise
                    from different cells, in different
                    organs, with different TF programmes,
                    different drugs, different outcomes.
                    They must be analyzed as separate
                    Waddington landscapes.

critical_note_2:    The depth score framework for ALL
                    requires explicit axis definition
                    BEFORE data loads:
                    Axis 1 = distance from arrested
                    developmental stage (within-subtype)
                    Axis 2 = distance toward primitive
                    HSC-like state (resistance axis)
                    The before-document must specify
                    which axis is primary.

critical_note_3:    HOXA cluster activation is the
                    common thread connecting the deepest
                    B-ALL (KMT2A infant) and the most
                    primitive T-ALL (ETP/HOXA T-ALL).
                    The same KMT2A/H3K4me3/HOXA
                    epigenetic mechanism operates in
                    both lineages. This cross-lineage
                    connection is the most structurally
                    significant observation in this
                    orientation document.

critical_note_4:    ALL is the only cancer in this
                    repository where CURE IS ACHIEVED
                    in >90% of affected children.
                    The depth score framework here is
                    not primarily about finding a cure —
                    it is about identifying the ~10%
                    of children who will fail standard
                    therapy BEFORE they fail, so that
                    those children receive the right
                    treatment from day one.
```
