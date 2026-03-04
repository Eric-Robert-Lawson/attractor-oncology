# BLADDER CANCER (UROTHELIAL CARCINOMA) — SUBTYPE ORIENTATION DOCUMENT
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

  A complete map of the bladder cancer molecular
  subtype landscape — what the four consensus
  molecular subtypes are, where each arises in
  the normal urothelial hierarchy, what the
  initiating and cooperating genetic events are,
  what the clinical characteristics and treatment
  are, and what public data exists to analyze
  each entity.

Bladder cancer has a structural feature that makes
it unique among the solid tumours in this repository:

  THE INVASION BINARY.

  All other solid tumours in this repository exist
  on a continuous depth axis: shallow to deep,
  differentiated to primitive, low-grade to
  high-grade. The clinical severity increases
  continuously with depth.

  In bladder cancer, there is a CATEGORICAL JUMP
  between two completely different disease states
  that share a tissue of origin:

  NMIBC — Non-Muscle Invasive Bladder Cancer:
    Confined to the urothelium and lamina propria.
    Does NOT invade the detrusor muscle.
    5-year survival: ~90%.
    Treated with transurethral resection + BCG.
    Managed endoscopically — the bladder is kept.
    The main threat: RECURRENCE and PROGRESSION.
    ~70% recur after initial treatment.
    ~15–30% progress to MIBC.

  MIBC — Muscle-Invasive Bladder Cancer:
    Has invaded the detrusor muscle (pT2+).
    5-year survival: ~50% with radical cystectomy.
    Treated with neoadjuvant chemotherapy +
    radical cystectomy or bladder preservation
    with chemoradiation.
    Once lymph node positive or metastatic:
    5-year survival falls to ~5–20%.

  The NMIBC → MIBC transition is not merely
  disease progression. It is a qualitative shift
  in the cancer's Waddington position — from a
  cell that retains partial urothelial identity
  and is epithelial-confined to a cell that has
  lost enough identity to breach the basement
  membrane and invade the muscle wall.

  This INVASION BINARY overlays the molecular
  subtype classification in a complex way:
    Luminal Papillary: predominantly NMIBC
                       (rarely invades muscle)
    Luminal Infiltrated: predominantly MIBC
    Basal/Squamous:    predominantly MIBC
                       (most aggressive)
    Neuronal:          almost always MIBC
                       (most aggressive of all)

  The framework must address both dimensions:
  the subtype (molecular identity) AND the
  invasion status (Waddington position relative
  to the muscle wall).
```

---

## DOCUMENT METADATA

```
document_id:        BLCA_Subtype_Orientation
series:             BLCA (Bladder Cancer — Subtypes)
folder:             Cancer_Research/BLCA/Subtypes/
date:               2026-03-03
author:             Eric Robert Lawson
status:             ORIENTATION ONLY — no predictions,
                    no analysis, no epigenetic hypotheses
next_document:      BLCA_LumPap_before.md
                    (Document BLCA-S1a — Luminal Papillary
                    before-doc)
protocol_version:   Workflow_Protocol.md v2.0
```

---

## SECTION I — THE NORMAL UROTHELIUM

```
The Waddington baseline for bladder cancer is the
normal urothelial epithelium — a stratified
specialised epithelium that lines the bladder,
ureters, and renal pelvis.

The urothelium is called TRANSITIONAL EPITHELIUM
because it was historically thought to be
"transitional" between squamous and columnar
epithelium. This is a misnomer — the urothelium
is a specialised tissue with unique biology that
does not simply transition between two other states.

═══════════════════════════════════════════════════════
THE UROTHELIAL HIERARCHY — THREE LAYERS
═══════════════════════════════════════════════════════

BASAL LAYER (deepest — resting on basement membrane)
  Cells:      Small, round, tightly packed
              Resting on lamina propria
  Markers:    KRT5 (Keratin 5) — HIGH
              KRT14 (Keratin 14) — expressed
              CD44 — stem-like marker
              CD49f — integrin, basal adhesion
              TP63 (p63) — THE master TF of
                           basal identity
  Function:   Progenitor compartment
              Responsible for urothelial
              regeneration after injury
              Relatively slow cycling under
              normal conditions — proliferate
              rapidly in response to damage
  Waddington: The most "primitive" normal
              urothelial position — highest
              self-renewal, least differentiated
              This layer is the cell of origin
              for BASAL/SQUAMOUS bladder cancer

INTERMEDIATE LAYER
  Cells:      Polygonal, multiple layers
              Between basal and umbrella
  Markers:    KRT5 — decreasing
              KRT8, KRT18 — rising
              GATA3 — present and rising
              CD44 — diminishing
  Function:   Transit-amplifying population
              Committed to luminal fate but
              not yet fully differentiated
              Active proliferative reserve
  Waddington: The mid-trajectory position —
              committed to luminal identity
              but not yet terminally
              differentiated

UMBRELLA LAYER (most superficial — faces the urine)
  Cells:      Very large, hexagonal, sometimes
              binucleated
              The surface layer — in direct
              contact with urine
  Markers:    KRT20 (Keratin 20) — HIGH
              Uroplakins: UPK1A, UPK1B,
                          UPK2, UPK3A, UPK3B
              (Uroplakins are unique to the
              urothelium — they form the
              asymmetric unit membrane that
              makes the bladder impermeable
              to urine)
              FOXA1 — high
              GATA3 — high
              PPARG — high
  Function:   The barrier.
              Prevents urine (hypertonic,
              acidic, containing waste
              metabolites) from diffusing
              back into the bladder wall.
              The uroplakin plaque is the
              molecular seal.
              Non-proliferating — terminally
              differentiated.
              Shed periodically, replaced
              by differentiating intermediate
              cells.
  Waddington: The terminal normal attractor —
              the most differentiated position
              in the urothelial hierarchy.
              This layer is the cell of origin
              for LUMINAL PAPILLARY bladder
              cancer (the most differentiated
              false attractor).

KEY TRANSCRIPTION FACTORS OF UROTHELIAL IDENTITY:

  TP63 (p63):
    The master TF of basal urothelial identity.
    TP63 does two things:
      1. Maintains basal progenitor programme
         (KRT5, KRT14, CD44, p63 itself)
      2. Represses luminal differentiation genes
         (uroplakins, KRT20, FOXA1)
    TP63 IS THE PAX5 OF THE UROTHELIUM —
    it is the guardian of basal identity.
    In basal/squamous bladder cancer: TP63
    is MAINTAINED/ACTIVATED → the cancer cell
    retains a basal progenitor identity and
    cannot access the luminal programme.
    In luminal bladder cancer: TP63 is
    DOWNREGULATED → the luminal programme
    is de-repressed.

  FOXA1 (HNF3α):
    The pioneer TF for luminal urothelial
    differentiation.
    Opens chromatin at uroplakin and GATA3
    loci in intermediate and umbrella cells.
    FOXA1 activates: UPK genes, KRT20,
    GATA3 targets, PPARG.
    FOXA1 is the entry point into the
    luminal false attractor landscape.
    In luminal papillary cancer: FOXA1 high.
    In basal/squamous: FOXA1 low/absent.

  GATA3:
    Luminal identity TF.
    GATA3 is the urothelial equivalent of
    PAX5 in B cells — it marks committed
    luminal identity throughout the
    intermediate and umbrella layers.
    GATA3 maintains urothelial identity
    (prevents squamous or glandular
    metaplasia).
    GATA3 loss in bladder cancer:
    → the cell can access squamous or
    glandular fates (as seen in CIS and
    basal/squamous MIBC).
    IMPORTANT: GATA3 is also the primary
    immunohistochemical marker used by
    pathologists to confirm urothelial
    origin in metastatic disease.
    If GATA3 is negative in a bladder
    cancer metastasis — the tumour has lost
    urothelial identity entirely.

  PPARG (PPARγ):
    Nuclear receptor — the master regulator
    of luminal urothelial differentiation.
    PPARG activation drives:
      Uroplakin gene expression
      FOXA1 and GATA3 target activation
      Luminal papillary tumour maintenance
    PPARG is the co-driver of the luminal
    papillary false attractor alongside FGFR3.
    PPARG amplification/activation is the
    founding event of many luminal papillary
    low-grade NMIBC — independent of FGFR3
    mutation.

THE WADDINGTON STRUCTURE OF NORMAL UROTHELIUM:

  The urothelium is a Waddington landscape
  with:
    Basal cell = the most stem-like normal
                 position (highest TP63,
                 lowest uroplakin)
    Umbrella cell = the terminal normal
                    attractor (lowest TP63,
                    highest uroplakin/KRT20)

  Normal maturation runs BASAL → UMBRELLA:
    Increasing FOXA1, GATA3, PPARG
    Increasing uroplakins, KRT20
    Decreasing TP63, KRT5, KRT14, CD44

  BLADDER CANCER FALSE ATTRACTOR STRUCTURE:
    Luminal Papillary: arrested near umbrella/
                       intermediate position
                       Retains FOXA1, GATA3,
                       PPARG, KRT20
                       FGFR3 constitutively
                       active → proliferates
                       from a near-differentiated
                       position
    Luminal Infiltrated: arrested at intermediate
                         position with stromal
                         infiltration
    Basal/Squamous:    arrested at basal/
                       progenitor position
                       Retains KRT5, KRT14, TP63
                       Cannot access luminal
                       programme
    Neuronal:          arrested in a non-
                       urothelial neuronal state —
                       most radical loss of
                       urothelial identity

  THE DEPTH AXIS FOR BLADDER CANCER:
    Normal umbrella cell (KRT20+, UPK+, TP63−)
    → Luminal Papillary (FGFR3+, PPARG+)
    → Luminal Infiltrated (stromal signal)
    → Basal/Squamous (KRT5+, TP63+, KRT20−)
    → Neuronal (synaptophysin+, chromogranin+)

    The depth axis runs from most differentiated
    (luminal papillary, near umbrella cell) to
    most primitive / most identity-lost (neuronal,
    entirely outside the urothelial landscape).
    This ordering is INVERTED relative to most
    other cancers in the repository — because
    luminal papillary is the FALSE ATTRACTOR
    CLOSEST TO THE NORMAL TERMINAL DIFFERENTIATED
    STATE. Most cancers arrest cells in primitive
    progenitor states. Luminal papillary bladder
    cancer arrests cells in a near-differentiated
    state.
```

---

## SECTION II — LUMINAL PAPILLARY: THE DOMINANT NMIBC SUBTYPE

```
PREVALENCE:    ~35% of all MIBC in TCGA
               BUT represents the vast majority
               of ALL bladder cancers when NMIBC
               is included (NMIBC is ~75% of all
               new bladder cancer diagnoses)
               The most common bladder cancer
               overall — low grade NMIBC is
               predominantly luminal papillary.

CLINICAL PROFILE:
  Invasion status:  Almost exclusively NMIBC
                    (Stage Ta, T1)
                    Rarely progresses to MIBC
                    WITHOUT first acquiring
                    additional mutations
  Grade:            Low-grade or high-grade
                    Low-grade papillary (LG-Ta):
                    ~70% of NMIBC at first
                    presentation
  Histology:        Papillary architecture
                    Thin fibrovascular cores
                    Urothelium with minimal
                    atypia in low-grade cases
  Prognosis:        EXCELLENT for low-grade
                    5-year OS: >90%
                    5-year cancer-specific
                    survival: >95%
                    PROBLEM: HIGH RECURRENCE
                    ~50–70% recur after TURBT
                    alone
                    ~5–10% of LG-Ta progress
                    to MIBC
  Standard therapy:
    TURBT (transurethral resection) — always
    + Intravesical BCG (high-risk NMIBC):
      Bacillus Calmette-Guérin stimulates
      local immune response in the bladder
      Reduces recurrence by ~30–40%
      Reduces progression to MIBC
    + Surveillance cystoscopy every 3 months
      (year 1), every 6 months (year 2),
      then annually — for life.
      The cystoscopy burden is the clinical
      reality of bladder cancer — patients
      undergo repeated procedures indefinitely.
    Erdafitinib (FGFR3 inhibitor):
      Approved for FGFR-altered advanced/
      metastatic urothelial carcinoma
      (not yet standard for NMIBC, but
      the target is most prevalent here)
  Response to immunotherapy:
    POOR as single agent for NMIBC
    BCG failure → pembrolizumab approved
    (but modest efficacy)

CELL OF ORIGIN:
  Intermediate or umbrella urothelial cell
  (or a committed progenitor in the upper
  intermediate layer) — the most differentiated
  cell of origin in the repository.

DEFINING MOLECULAR EVENTS:

  FGFR3 mutation (~60% of low-grade NMIBC):
    The most common oncogenic mutation in
    bladder cancer overall.
    FGFR3 is a receptor tyrosine kinase
    for fibroblast growth factors (FGF7,
    FGF10 in urothelium).
    Activating mutations: S249C (most common),
    Y375C, R248C — create ligand-independent
    dimerization and constitutive signalling.
    FGFR3 mutation activates:
      RAS/MAPK → proliferation
      PI3K/AKT → survival
    FGFR3 in bladder cancer is the OPPOSITE
    of FGFR3 in skeletal development — in
    bone, FGFR3 activation is ANTI-proliferative
    (achondroplasia mutations cause premature
    growth plate closure). In urothelium,
    FGFR3 activation is PRO-proliferative.
    This tissue-context specificity is
    important for understanding why FGFR3
    mutation produces a LOW-GRADE, WELL-
    DIFFERENTIATED tumour in urothelium —
    the constitutive signal maintains the
    cell in a near-differentiated proliferative
    state, not a primitive one.

  FGFR3 fusion/translocation (~10% of MIBC):
    FGFR3-TACC3 fusion — high FGFR3 kinase
    activity. Present in a subset of luminal
    MIBC and targetable by erdafitinib.

  PPARG amplification/activation (~50% of
  luminal papillary NMIBC):
    PPARG locus amplification at 3p25.
    OR: RXRA mutation (S427F) which activates
    PPARG signalling.
    PPARG activation drives the luminal
    differentiation programme independently
    of FGFR3.
    PPARG is the second major driver of
    luminal papillary identity.

  STAG2 mutation (~26% of Ta NMIBC):
    Cohesin complex component.
    STAG2 mutations are enriched in NMIBC.
    Associated with increased chromosomal
    instability — a paradox for a
    well-differentiated tumour.

  ARID1A, KDM6A mutation (~20–25%):
    Chromatin remodelling genes.
    KDM6A is a histone demethylase
    (H3K27me3 demethylase).
    KDM6A LOSS → H3K27me3 elevated at
    differentiation gene loci → epigenetic
    silencing of luminal identity genes.
    KDM6A is one of the most commonly
    mutated genes in bladder cancer overall.

  TP53, RB1:    RARELY mutated in NMIBC
                LOW-GRADE luminal papillary
                Their mutation marks the
                NMIBC → MIBC transition.

CORE EXPRESSION PROGRAMME:
  DEPTH-POSITIVE (rises as tumour deepens
  within luminal papillary category):
    FGFR3 and FGFR3 pathway targets
    CCND1 (cyclin D1) — cell cycle
    MYC targets — proliferation
    E2F targets
  DEPTH-NEGATIVE (falls as lumour deepens):
    UPK1A, UPK1B, UPK2, UPK3A (uroplakins)
    KRT20
    FOXA1
    GATA3
    PPARG targets
    CDH1 (E-cadherin — epithelial cohesion)
  THE LUMINAL PAPILLARY DEPTH AXIS:
    Shallow: near-normal urothelium, high UPK,
             high GATA3, FGFR3-mutant
             → essentially a papilloma with
             a constitutively active RTK
    Deep:    Low UPK, GATA3 dimming, TP53
             acquiring mutations
             → the NMIBC → MIBC boundary
```

---

## SECTION III — BASAL/SQUAMOUS: THE DOMINANT MIBC SUBTYPE

```
PREVALENCE:    ~35% of MIBC (TCGA)
               The most common muscle-invasive
               subtype. Almost all MIBC
               patients who are basal/squamous
               present with MIBC — NMIBC does
               not usually precede it.

CLINICAL PROFILE:
  Invasion status:  Almost always MIBC (pT2–T4)
                    Frequently presents at
                    advanced stage (T3/T4)
                    with lymph node involvement
  Histology:        Often shows squamous
                    differentiation or squamous
                    morphology. CK5/6 positive
                    by IHC. Occasionally shows
                    sarcomatoid features at
                    the invasive front.
  Prognosis:        POOR
                    5-year OS: ~30–40% after
                    radical cystectomy
                    Metastatic: ~10–15%
  Standard therapy:
    Neoadjuvant cisplatin/gemcitabine → radical
    cystectomy (MVAC or GC regimen)
    BUT: basal/squamous MIBC may be
    LESS responsive to cisplatin-based NAC
    than luminal MIBC (evidence is
    conflicting — the before-document will
    address this as a prediction to test)
    Enfortumab vedotin + pembrolizumab:
    NOW FIRST-LINE STANDARD for advanced/
    metastatic urothelial carcinoma
    (EV-302 trial) — supersedes platinum
    chemotherapy in unselected patients
  Response to immunotherapy:
    HIGHER response rate for basal/squamous
    than luminal — the immune infiltrate
    is higher, PD-L1 is more often elevated,
    and the basal/squamous false attractor
    appears more "visible" to T cells.
    This is the structural paradox:
    the most aggressive subtype is also
    the most immunotherapy-responsive.
    (Structural parallel to CMS1 in CRC —
    the aggressive subtype responds to
    checkpoint blockade.)

CELL OF ORIGIN:
  Basal urothelial cell (KRT5+, KRT14+, TP63+)
  The progenitor compartment of the urothelium.
  Basal cells are the repair cells — they
  respond to injury and inflammation by
  proliferating. The persistent activation of
  basal progenitor programmes (TP63, KRT5,
  KRT14) without accessing the differentiation
  programme is the false attractor.

DEFINING MOLECULAR EVENTS:

  TP63 activation / maintained expression:
    TP63 is the master TF of basal urothelial
    identity. Its persistent expression in
    basal/squamous MIBC maintains the basal
    programme and REPRESSES the luminal
    differentiation programme.
    TP63 drives:
      KRT5, KRT14 expression
      CD44 expression (stemness)
      Repression of FOXA1, GATA3 targets
      Repression of uroplakins
    TP63 is both a marker of basal identity
    AND a driver of the false attractor engine
    in basal/squamous MIBC.
    This makes TP63 conceptually analogous to
    PAX5 in B cells — but INVERTED: in B-ALL,
    PAX5 LOSS drives the false attractor.
    In basal bladder cancer, TP63 RETENTION
    drives the false attractor. The normal
    uroplakin/luminal programme is repressed
    by TP63 maintaining basal identity.

  TP53 mutation (~50% of basal/squamous):
    NOT the same as TP63.
    TP53 mutation loses genome guardian
    function → genomic instability →
    accumulation of additional mutations.
    TP53 mutation is absent in luminal
    papillary NMIBC and present primarily
    in MIBC. Its presence marks the
    acquisition of invasive capacity.

  RB1 deletion/mutation (~20%):
    Retinoblastoma tumour suppressor —
    cell cycle checkpoint.
    RB1 loss → E2F transcription factors
    constitutively active → unrestrained
    G1/S progression.
    RB1 deletion in bladder cancer is
    enriched in basal/squamous and predicts
    poor response to cisplatin NAC.

  CDKN2A deletion (~40%):
    p16 and p14ARF — same cell cycle
    brake as in T-ALL.
    Loss is near-universal in MIBC
    regardless of subtype but enriched
    in basal/squamous.

  EGFR amplification (~10%):
    The basal/squamous subtype is enriched
    for EGFR copy number gains.
    EGFR drives the basal proliferative
    programme (analogous to FGFR3 in
    luminal but through a different receptor).

IMMUNE ARCHITECTURE:
  Basal/squamous MIBC has a COMPLEX immune
  landscape — different from simple immune
  exclusion:

  CD8+ T cells:   Present but variable
                  Some Ba/Sq tumours are
                  immune-infiltrated ("hot")
                  Others are immune-excluded
  PD-L1:          Often HIGH on tumour cells
                  and stroma
  Tumour-associated macrophages (TAMs):
                  Elevated — often M2 polarised
                  (immunosuppressive) in
                  invasive front
  Cancer-associated fibroblasts (CAFs):
                  Prominent — create stromal
                  desmoplasia at invasive front
  NK cells:       Reduced

  THE IMMUNE PARADOX:
  Basal/squamous MIBC is simultaneously:
    The highest immune infiltrate MIBC subtype
    The most aggressive and lethal subtype
  Resolution: The immune cells are present
  but EXHAUSTED or functionally impaired.
  PD-L1 on tumour + checkpoint exhaustion
  prevent killing despite T cell presence.
  Enfortumab vedotin + pembrolizumab
  overcomes this: EV kills the cancer cell
  directly via ADC while pembrolizumab
  releases the exhausted T cells.
  The combination targets both axes
  simultaneously.

CORE EXPRESSION PROGRAMME:
  DEPTH-POSITIVE (rises with depth):
    KRT5, KRT14, KRT6A/B (squamous keratins)
    TP63
    CD44, CD49f (stemness markers)
    EGFR
    ZEB1, ZEB2 (EMT drivers at invasive front)
    VIM (vimentin — mesenchymal marker)
    MMP2, MMP9 (invasion)
    CDH2 (N-cadherin — EMT)
  DEPTH-NEGATIVE (falls with depth):
    UPK1A, UPK1B, UPK2 (uroplakins — LOST)
    KRT20 (LOST)
    FOXA1 (LOST)
    GATA3 (progressively lost)
    CDH1 (E-cadherin — lost at invasive front)
    PPARG (LOST)
```

---

## SECTION IV — LUMINAL INFILTRATED: THE INTERMEDIATE SUBTYPE

```
PREVALENCE:    ~19% of MIBC (TCGA)
               The most immunologically active MIBC
               subtype.

CLINICAL PROFILE:
  Invasion status:  MIBC
  Histology:        Luminal markers present but
                    tumour stroma is heavily
                    infiltrated by immune cells
                    (T cells, B cells,
                    macrophages, fibroblasts)
  Prognosis:        INTERMEDIATE
                    5-year OS: ~45%
                    Better than basal/squamous
                    at matched stage but worse
                    than luminal papillary
  Standard therapy: Cisplatin/gemcitabine NAC
                    The immune infiltrate in
                    luminal infiltrated MIBC
                    predicts:
                    a) Better response to
                       cisplatin NAC (immune
                       cells amplify chemo
                       response)
                    b) Potential response to
                       immunotherapy (high
                       immune infiltrate)
  Response to immunotherapy:
    COMPLEX — luminal infiltrated has the
    most complex immune architecture.
    "Immune active" subtype — LIKELY to
    respond to checkpoint blockade BUT
    also contains regulatory elements
    (Tregs, M2 macrophages, TGF-β) that
    can suppress response.
    Pembrolizumab trial (KEYNOTE-052):
    luminal infiltrated had variable
    response — not cleanly predictable
    by subtype alone.

CELL OF ORIGIN:
  Intermediate urothelial cell — similar to
  luminal papillary cell of origin but with
  stromal/immune co-activation.
  The "infiltrated" designation refers to
  the tumour microenvironment, not just
  the cancer cell itself.
  CMS4 parallel: like CMS4 in CRC, the
  luminal infiltrated subtype has a heavy
  stromal/immune component in the bulk
  RNA-seq signal. The tumour cell identity
  is luminal; the microenvironment signal
  is immune/stromal. Bulk RNA-seq mixes both.

DEFINING MOLECULAR EVENTS:
  FGFR3 mutation: LOW (~10%)
                  (contrast to luminal papillary
                  where it is ~60%)
  TP53 mutation:  Intermediate (~35%)
  CDKN2A:         Frequent deletion
  TGF-β pathway:  ACTIVE
                  Drives immune infiltration
                  AND produces an
                  immunosuppressive stroma
                  TGF-β is both the
                  attractor of immune cells
                  and the mechanism of
                  their suppression
  EMT markers:    Intermediate — partial EMT
                  at the invasive front
  IMMUNE MARKERS:
    CD8+ T cells:  HIGH (by definition)
    FOXP3+ Tregs:  Elevated
    B cells:       Present (lymphoid aggregates
                   visible histologically —
                   "tertiary lymphoid structures")
    PD-L1:         Variable
    TGF-β:         HIGH (same as CMS4 in CRC)

STRUCTURAL NOTE FOR THE FRAMEWORK:
  Luminal infiltrated MIBC is the most complex
  subtype for the depth score framework.
  Like CMS4 in CRC:
    Bulk RNA-seq captures both tumour cell
    transcriptome AND immune/stromal
    transcriptome.
    The depth score in luminal infiltrated
    may partly reflect immune infiltrate
    intensity rather than tumour cell
    depth.
  This must be addressed explicitly in the
  before-document: the primary axis is
  tumour cell depth (loss of urothelial
  identity), and the immune signal is a
  secondary axis to be measured separately.
```

---

## SECTION V — NEURONAL: THE RAREST AND MOST AGGRESSIVE

```
PREVALENCE:    ~5% of MIBC (TCGA)
               May be underrepresented in TCGA
               because it is often treated as
               a separate entity (neuroendocrine
               bladder cancer) in clinical trials

CLINICAL PROFILE:
  Invasion status:  MIBC — always
  Histology:        Small cell / neuroendocrine
                    morphology in the most
                    extreme cases.
                    Express neuroendocrine
                    markers: synaptophysin,
                    chromogranin A, CD56,
                    NCAM1.
                    May have mixed
                    urothelial/neuroendocrine
                    histology.
  Prognosis:        WORST of all subtypes
                    5-year OS: <20%
                    Often presents with
                    distant metastases
  Standard therapy: Platinum-based regimens
                    (cisplatin/etoposide for
                    small cell component)
                    Immunotherapy: limited data
  Response to immunotherapy:
    POOR — neuroendocrine tumours generally
    have low mutation burden and reduced
    antigen presentation (MHC-I frequently
    lost)

CELL OF ORIGIN:
  The neuronal subtype represents the most
  extreme loss of urothelial identity in the
  landscape.
  Two hypotheses:
    a) Transdifferentiation from a urothelial
       cell — the cancer cell has undergone
       such extreme lineage plasticity that it
       has accessed a neuroendocrine programme
       entirely outside the urothelial
       Waddington landscape.
    b) Origin from a rare neuroendocrine
       progenitor cell in the urothelium
       (enteroendocrine equivalent — rare
       neuroendocrine cells are present in
       most mucosal epithelia).

DEFINING MOLECULAR EVENTS:
  TP53 mutation: ~90% — universal
  RB1 deletion:  ~90% — universal
  CDKN2A:        ~80%
  NEUROD1:       Neuroendocrine TF — elevated
  ASCL1:         Neuroendocrine TF — elevated
  RET:           Occasionally
  NOTCH1:        Often inactivated
                 (NOTCH1 maintains epithelial
                 identity and suppresses
                 neuroendocrine fate —
                 its loss releases the
                 neuroendocrine programme)
  FOXA1, GATA3:  LOST
  TP63:          LOST

THE WADDINGTON INTERPRETATION:
  Neuronal MIBC is the deepest false attractor
  in the bladder cancer landscape — not in the
  sense of being "more stemlike" but in the
  sense of being the MOST IDENTITY-LOST.
  The cell has abandoned the urothelial
  Waddington landscape entirely and entered
  a neuroendocrine false attractor.
  This is analogous to the mesenchymal shift
  in GBM recurrence — the tumour has changed
  its landscape entirely under selective
  pressure (often in response to treatment).
  TREATMENT-EMERGENT NEURONAL/NE TRANSDIFFERENTIATION:
  Some bladder cancers that are initially
  luminal or basal acquire neuronal markers
  after immunotherapy or targeted therapy —
  this transdifferentiation is analogous to
  small-cell transformation in EGFR-mutant
  lung cancer after osimertinib.
  The neuronal subtype at diagnosis is a
  primary false attractor.
  The neuronal subtype after treatment is a
  treatment-selected escape.
  Both look the same on RNA-seq.
  This distinction must be flagged in the
  before-document.
```

---

## SECTION VI — THE NMIBC AXIS: CIS AND HIGH-GRADE T1

```
Non-Muscle Invasive Bladder Cancer has an internal
depth axis that the four-subtype consensus
classification does not fully capture.

NMIBC is not one disease — it has three clinical
and molecular states:

LOW-GRADE Ta (LG-Ta):
  The most benign NMIBC.
  Papillary architecture, minimal atypia.
  Recurs frequently but almost never
  progresses to MIBC.
  Predominantly luminal papillary subtype.
  FGFR3 mutation rate: ~80%.
  Management: TURBT + surveillance.
  Treatment burden is HIGH (repeated
  cystoscopies) but cancer mortality is LOW.

HIGH-GRADE Ta / T1 (HG-Ta, HG-T1):
  More aggressive papillary or flat tumours.
  High-grade T1: has invaded the lamina propria
  but not the muscle.
  The most dangerous NMIBC — ~30–50% progress
  to MIBC.
  Molecular characteristics:
    FGFR3 mutation: ~30–40%
    TP53 mutation: beginning to accumulate
    CDKN2A deletion: ~30%
    Some acquire basal markers
  Management: BCG intravesical immunotherapy
  BCG failure: either early radical cystectomy
  or pembrolizumab/nadofaragene.

CARCINOMA IN SITU (CIS):
  Not papillary — flat, high-grade.
  The most dangerous NMIBC in terms of
  progression risk.
  CIS is molecularly closer to MIBC than to
  LG-Ta:
    TP53 often mutated
    GATA3 dimming
    KRT5 sometimes upregulating
    The cell is transitioning from luminal
    toward basal/squamous identity.
  CIS is the Waddington transition state:
  the cell is losing luminal identity but
  has not yet established basal identity.
  It is at the saddle point between the
  luminal and basal false attractors.
  Management: BCG — the most urgent BCG
  indication. CIS without BCG → ~50%
  progress to MIBC within 5 years.

THE CLINICAL URGENCY HIERARCHY:
  LG-Ta:    Recurs, rarely progresses
            → Surveillance
  HG-Ta:    Recurs, may progress
            → BCG + close surveillance
  T1-HG:    Invades lamina propria
            → BCG + consider early cystectomy
  CIS:      Flat, high-grade, unstable
            → BCG + cystectomy if BCG fails
  MIBC:     Muscle invasion confirmed
            → NAC + cystectomy or
              chemoradiation

THE DEPTH SCORE IN NMIBC:
  The depth score framework applied to
  NMIBC-specific datasets should separate:
    LG-Ta (shallow — near-normal luminal
           urothelial identity, FGFR3+)
    HG-Ta/T1 (mid depth — losing some
              identity, TP53 accumulating)
    CIS (deep NMIBC — at the transition
         to invasive potential)
  A depth score that predicts progression
  from NMIBC to MIBC would be the highest-
  clinical-value output for the bladder
  cancer framework — identifying the 15-30%
  of NMIBC patients who will progress BEFORE
  they progress.
```

---

## SECTION VII — ENFORTUMAB VEDOTIN: THE CROSS-SUBTYPE DRUG

```
The most important clinical advance in bladder
cancer in the last decade is the approval of
enfortumab vedotin + pembrolizumab (EV+P) as
first-line therapy for locally advanced or
metastatic urothelial carcinoma.

EV-302 / KEYNOTE-A39 trial (2023-2024):
  EV+P vs. platinum chemotherapy
  Median OS: 31.5 months (EV+P) vs. 16.1 months
  Hazard ratio: 0.47 — 53% reduction in death
  This is the largest survival benefit ever
  demonstrated in first-line metastatic
  bladder cancer.

ENFORTUMAB VEDOTIN — MECHANISM:
  Antibody-drug conjugate (ADC):
  Anti-NECTIN4 antibody attached to MMAE
  (monomethyl auristatin E, a microtubule
  poison) via a cleavable linker.

  NECTIN4 is the key:
    NECTIN4 is highly expressed on urothelial
    carcinoma cells across ALL subtypes —
    luminal papillary, basal/squamous,
    luminal infiltrated, neuronal.
    It is not subtype-specific.
    This is a UNIVERSAL BLADDER CANCER
    SURFACE TARGET.
    Why is NECTIN4 high on bladder cancer?
    NECTIN4 is an adhesion molecule involved
    in the normal urothelial barrier.
    Urothelial cancers maintain NECTIN4
    expression even as they lose other
    urothelial identity markers.
    NECTIN4 is the urothelial identity
    marker that SURVIVES the transition to
    cancer — it is the last vestige of
    urothelial identity in every subtype.

  MECHANISM OF KILL:
    EV binds NECTIN4 on the cancer cell
    surface → internalised → MMAE released
    intracellularly → microtubule disruption
    → mitotic arrest → apoptosis.
    BYSTANDER EFFECT: MMAE can diffuse to
    adjacent cells — kills NECTIN4-negative
    cells in the tumour microenvironment.

  THE DEPTH SCORE PREDICTION:
    NECTIN4 expression should be a depth-
    NEGATIVE gene: it should fall with
    depth as urothelial identity is lost.
    IF NECTIN4 falls with depth → the most
    deeply positioned tumours (neuronal,
    deep basal/squamous) would be LESS
    responsive to EV because they have
    lost the target.
    This is a structural hypothesis —
    the before-documents will test it.

PEMBROLIZUMAB MECHANISM:
  Anti-PD-1 checkpoint inhibitor.
  Releases exhausted CD8+ T cells.
  More effective in:
    Basal/squamous (higher immune infiltrate)
    Luminal infiltrated (highest immune
    infiltrate)
  Less effective in:
    Luminal papillary (low immune infiltrate)
    Neuronal (MHC-I loss)
```

---

## SECTION VIII — THE EXISTING BLCA ANALYSIS — CONTEXT

```
The existing analysis in Cancer_Research/BLCA/
ran on the combined TCGA-BLCA dataset.

TCGA-BLCA contains ~412 MIBC samples with:
  Luminal Papillary:   ~35% = ~144 samples
  Luminal Infiltrated: ~19% = ~78 samples
  Basal/Squamous:      ~35% = ~144 samples
  Neuronal:            ~5%  = ~20 samples
  Luminal/Stroma-rich: ~6%  = ~25 samples

NOTE: TCGA-BLCA is primarily MIBC.
The low-grade NMIBC (LG-Ta) is massively
underrepresented — most patients with LG-Ta
do not undergo cystectomy and are not in
TCGA. The TCGA dataset therefore captures
primarily the muscle-invasive disease and
is BIASED toward higher-grade disease.

WHAT THE EXISTING COMBINED ANALYSIS CAPTURED:
  The depth axis from combined TCGA-BLCA is
  the BASAL/SQUAMOUS vs. LUMINAL PAPILLARY
  axis — the two largest subgroups pulling
  in opposite directions:
    Basal/squamous: KRT5+, TP63+, KRT14+,
                    deep in the Waddington
                    landscape (basal progenitor)
    Luminal papillary: KRT20+, UPK+, FOXA1+,
                    shallow (near-differentiated)
  The combined depth axis runs from the
  umbrella cell pole (shallow) to the basal
  progenitor pole (deep).
  The depth-positive genes should be:
    KRT5, KRT14, KRT6A (basal keratins)
    TP63
    CD44
    ZEB1 (EMT marker at invasive front)
    VIM
  The depth-negative genes should be:
    UPK1A, UPK1B, UPK2, UPK3A (uroplakins)
    KRT20
    FOXA1
    GATA3
    PPARG targets

THE INVERSION NOTE:
  The depth score in BLCA is INVERTED relative
  to most other cancers in the repository.
  In most cancers: depth-positive = primitive
  (HSC-like, stem-like)
  In BLCA: depth-positive = BASAL PROGENITOR
  (most aggressive) but depth-negative =
  UMBRELLA CELL (most differentiated)
  The MOST DIFFERENTIATED false attractor
  (luminal papillary) is the SHALLOWEST
  and the LEAST DIFFERENTIATED is the DEEPEST.
  This is anatomically correct:
  The bladder cancer depth axis runs from
  the apical surface inward — the umbrella
  cell is at the surface (shallow) and the
  basal cell is at the basement membrane
  (deep). The cancer cell that has retreated
  to the basement membrane position is the
  deepest cancer.
  This anatomical depth axis corresponds to
  the clinical depth axis: MIBC is literally
  deeper in the bladder wall than NMIBC.
```

---

## SECTION IX — DATA AVAILABILITY SUMMARY

```
Dataset         Accession        n (tumour)  Normal    Subtype   Power
                                             (n)       labels
────────────────────────────────────────────────────────────────────────
TCGA-BLCA       GDC/phs000178    ~412 MIBC   19        YES       HIGH
                cBioPortal                   (adj)     (Kamoun
                                                       consensus)
IMvigor210      GSE93157         348          0        YES       MOD-HIGH
(atezolizumab                   (met UC)               (mapped)
trial)
GSE13507        GSE13507         165          23       YES       MODERATE
                                (NMIBC+MIBC)           (partial)
GSE32894        GSE32894         308          12       YES       HIGH
                                                       (Sjödahl)
GSE48075        GSE48075         ~70          ~10      YES       LOW-MOD
                                NMIBC focus
GSE124305       GSE124305        ~100         ~10      YES       MODERATE
(NAC response                   (NAC)                  + NAC
cohort)                                                response
Normal urothelium references:
  GTEx:          Bladder tissue n=21 normal
                 (mixed, not sorted by layer)
  TCGA-BLCA adj: n=19 (adjacent normal,
                 proximity-effect risk)
  GSE13507 normal: n=23 (normal urothelium —
                 best available for BLCA)

CRITICAL NOTES:
  NMIBC IS UNDERREPRESENTED IN ALL PUBLIC
  DATASETS. TCGA-BLCA is essentially all MIBC.
  GSE13507 is the best source for NMIBC
  with normal urothelium controls but is
  a microarray study, not RNA-seq.
  For NMIBC-specific depth score analysis:
  GSE48075 + GSE13507 NMIBC subset are
  the primary sources.
  For MIBC subtype analysis:
  TCGA-BLCA + IMvigor210 provide the best
  combined depth and treatment response data.
  IMvigor210 is uniquely valuable because:
    It contains treatment response data
    (atezolizumab — anti-PD-L1)
    It has consensus subtype labels
    It enables prediction testing:
    "does depth score predict immunotherapy
    response within each subtype?"
    This is the highest-value clinical
    output from the BLCA framework.

CONSENSUS SUBTYPE LABELLING:
  The Kamoun 2020 (Nature Reviews Urology)
  consensus classifier is the standard.
  It is available as an R package
  (consensusMIBC) that takes RNA-seq data
  as input and assigns subtype labels.
  This must be run on TCGA-BLCA and
  IMvigor210 before any subtype analysis —
  it is a prerequisite data step, not a
  prediction.
```

---

## SECTION X — PLANNED ANALYSIS ORDER

```
ORDER:

  BLCA-S1   Luminal Papillary   TCGA-BLCA       HIGH POWER
            (MIBC component)    LumPap subset
                                + GSE32894
            REASON: The canonical FGFR3/PPARG-
            driven bladder cancer. Represents
            the vast majority of all bladder
            cancers by number (including all
            NMIBC). The existing combined
            analysis signal is dominated by
            its contrast with basal/squamous.
            FGFR3 depth correlation and
            erdafitinib response prediction
            is the clinical output.

  BLCA-S2   Basal/Squamous      TCGA-BLCA       HIGH POWER
            MIBC                Ba/Sq subset
                                + IMvigor210
                                Ba/Sq subset
            REASON: The most aggressive and
            most prevalent MIBC subtype.
            The basal progenitor false attractor.
            Cisplatin NAC response prediction.
            Enfortumab vedotin + pembrolizumab
            response prediction from depth score
            and NECTIN4 expression.
            The NECTIN4 depth-negative hypothesis
            will be tested here first.

  BLCA-S3   Luminal Infiltrated TCGA-BLCA       MODERATE
                                LumInf subset
                                + IMvigor210
                                LumInf subset
            REASON: The immunotherapy-responsive
            MIBC. TGF-β immune exclusion axis.
            Atezolizumab/pembrolizumab response
            data from IMvigor210 enables direct
            testing of: "does depth score predict
            immunotherapy response in luminal
            infiltrated MIBC?"
            The CMS4/CRC structural parallel
            (TGF-β, immune exclusion, stromal
            contamination in bulk RNA-seq) must
            be addressed.

  BLCA-S4   NMIBC               GSE13507        LOW-MOD
            (LG-Ta and CIS)     NMIBC subset
                                + GSE48075
            REASON: The clinical question with
            the highest unmet need: which NMIBC
            patient will progress to MIBC?
            A depth score that separates LG-Ta
            (no progression risk) from CIS/HG-T1
            (high progression risk) is the most
            clinically actionable output in
            bladder cancer.
            Low n — analysis will be LOW-POWER
            and predictions adjusted accordingly.

  BLCA-X    Cross-subtype       After S1–S4
                                Questions:
                                  1. NECTIN4:
                                     depth-negative
                                     across all
                                     subtypes?
                                  2. GATA3/KRT20
                                     as universal
                                     depth-negative
                                     markers?
                                  3. Does depth
                                     score predict
                                     EV+pembrolizumab
                                     response in
                                     IMvigor210?
                                  4. CIS depth
                                     score vs. MIBC
                                     — does CIS sit
                                     at the transition
                                     point?
                                  5. KDM6A mutation
                                     and depth —
                                     does KDM6A loss
                                     shift tumours
                                     toward deeper
                                     positions?
```

---

## SECTION XI — WHAT THIS DOCUMENT DOES NOT CONTAIN

```
This document contains:
  ✓ The invasion binary (NMIBC vs. MIBC) as the
    unique structural feature of bladder cancer
  ✓ Normal urothelial hierarchy (basal → umbrella)
  ✓ All four consensus molecular subtypes with
    cells of origin and molecular drivers
  ✓ The NMIBC internal depth axis (LG-Ta → CIS)
  ✓ Enfortumab vedotin mechanism and NECTIN4
    as a cross-subtype target
  ✓ Clinical facts (survival, standard of care,
    BCG, cisplatin, EV+pembrolizumab)
  ✓ The depth score inversion note (luminal
    papillary = SHALLOWEST in the bladder cancer
    landscape despite being the most common
    cancer subtype)
  ✓ IMvigor210 as the treatment-response dataset
  ✓ Data availability (TCGA-BLCA, IMvigor210,
    GEO datasets, normal urothelium reference)
  ✓ Context for the existing combined BLCA analysis

This document does NOT contain:
  ✗ Depth score predictions
  ✗ Specific false attractor gene predictions
  ✗ Drug target predictions
  ✗ Epigenetic lock mechanism hypotheses
  ✗ Cross-subtype structural predictions beyond
    what is established in the literature

All of the above belong in the BEFORE documents.
BLCA-S1a (Luminal Papillary before-document)
is next. Written before any script runs.
Before any data loads.
```

---

## STATUS BLOCK

```
document:           BLCA_Subtype_Orientation.md
folder:             Cancer_Research/BLCA/Subtypes/
status:             COMPLETE — ORIENTATION ONLY
date:               2026-03-03
author:             Eric Robert Lawson
                    OrganismCore

subtypes_covered:
  Luminal Papillary:   FGFR3/PPARG/umbrella cell  [1 of 4]
  Basal/Squamous:      TP63/KRT5/basal progenitor  [2 of 4]
  Luminal Infiltrated: TGF-β/immune/stromal        [3 of 4]
  Neuronal:            NE transdiff/TP53/RB1       [4 of 4]

analyses_started:   0

existing_analysis:  Cancer_Research/BLCA/ — run on
                    combined TCGA-BLCA (MIBC-dominant).
                    Depth axis is the
                    luminal papillary ↔ basal/squamous
                    contrast: umbrella cell pole
                    (shallow) vs. basal progenitor
                    pole (deep). Depth-positive genes
                    should be KRT5, KRT14, TP63.
                    Depth-negative genes should be
                    UPK1A, KRT20, FOXA1, GATA3.
                    NMIBC is underrepresented — the
                    low-grade papillary pole of the
                    landscape is not fully captured.

next_document:      BLCA-S1a
                    Luminal Papillary Before-Document
                    (predictions locked before
                    TCGA-BLCA LumPap subset loads)

critical_note_1:    The depth score in BLCA is
                    INVERTED relative to most other
                    repository cancers.
                    Most cancers: deep = primitive,
                    poor prognosis.
                    BLCA: deep = basal progenitor
                    (most aggressive), shallow =
                    near-differentiated umbrella cell
                    (most favourable).
                    The before-documents must flag
                    this inversion so the drug target
                    logic is applied correctly.

critical_note_2:    NMIBC → MIBC progression prediction
                    is the highest-priority clinical
                    output for the BLCA framework.
                    A depth score that identifies
                    which NMIBC patient will progress
                    to MIBC BEFORE they progress
                    would change clinical management
                    for the ~15–30% of patients who
                    currently receive inadequate
                    initial treatment.
                    This requires NMIBC-specific data
                    (GSE13507, GSE48075) — not just
                    TCGA-BLCA which is MIBC-dominant.

critical_note_3:    NECTIN4 as a depth-negative marker
                    is a structural hypothesis to be
                    tested in BLCA-S1a and BLCA-S2a.
                    If NECTIN4 falls with depth:
                    enfortumab vedotin will be LESS
                    effective in the deepest tumours.
                    If NECTIN4 is depth-stable:
                    EV efficacy is depth-independent.
                    IMvigor210 (atezolizumab trial,
                    not EV trial) cannot directly
                    test EV response — but expression
                    of NECTIN4 in IMvigor210 can test
                    whether depth correlates with
                    NECTIN4 level.

critical_note_4:    The TP63 RETENTION (basal/squamous)
                    vs. TP63 SUPPRESSION (luminal
                    papillary) axis is the primary
                    molecular identity switch in
                    bladder cancer. It is the
                    urothelial equivalent of:
                    PAX5 in B cells (loss drives
                    false attractor)
                    CEBPA in AML (loss drives
                    granulocytic arrest)
                    CDX2 in CRC (loss marks depth)
                    TP63 RETENTION is the mechanism
                    that maintains basal identity and
                    blocks luminal differentiation in
                    the deepest false attractor.
```
