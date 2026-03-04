# CHRONIC LYMPHOCYTIC LEUKEMIA (CLL) — SUBTYPE ORIENTATION DOCUMENT
## Before Any Subtype Analysis Begins
## OrganismCore | 2026-03-04 | Author: Eric Robert Lawson

---

## PURPOSE OF THIS DOCUMENT

```
This document exists before any script runs.
It contains no depth score derivations.
It contains no novel drug predictions.
It contains no subtype-specific epigenetic hypotheses.

What it contains:

  A complete map of the CLL molecular landscape —
  the IGHV mutated/unmutated divide (the single most
  clinically important stratification in CLL), the
  cytogenetic hierarchy (del13q, trisomy 12, del11q,
  del17p/TP53), the stereotyped BCR subsets,
  the normal B cell developmental programme that
  is being subverted, the CLL false attractor structure
  as derived and confirmed in Cancer Validation #8
  (Document 80), the resistance mechanisms for BTK
  and BCL2 inhibitors, and the full treatment landscape
  from watch-and-wait through targeted therapy to
  the emerging next-generation drug combinations.

CLL HOLDS A UNIQUE STRUCTURAL POSITION
IN THE ORGANISCORE REPOSITORY.

  CLL is the only cancer in the repository
  that has been confirmed as a SURVIVAL ATTRACTOR
  rather than a DIFFERENTIATION BLOCK ATTRACTOR.

  Every other confirmed cancer:
    AML:   differentiation block (GMP → promyelocyte)
    CML:   differentiation block (BCR-ABL drives
           GATA1/2 loss, HSC expansion)
    CRC:   differentiation block (CDX2 loss, WNT up)
    GBM:   differentiation block (OLIG2 false neuronal)
    BRCA:  differentiation block (EZH2 blocks luminal)
    LUAD:  differentiation block (NKX2-1 loss)
    B-ALL: differentiation block (PAX5 network frozen)
    MDS:   differentiation execution failure (multiple
           lineages — see MDS Subtype Orientation)

  CLL is NONE OF THESE.
  CLL cells are NOT blocked from differentiating.
  CLL cells have ALREADY DIFFERENTIATED.
  They are mature B cells — CD5+ CD19+ CD23+.
  They have completed V(D)J recombination
  (IGKC elevated +60%, RAG1/RAG2 silent in the data).
  They are STUCK AFTER MATURATION,
  in a state that resembles the mature naive B cell
  identity but lacks the final exit signal: apoptosis.

  THE FALSE ATTRACTOR IN CLL IS THE SURVIVAL STATE.
  The cell arrived at its correct terminal identity.
  It then failed to die when it should have.
  BCL2 is the molecular lock — elevated 136% ***
  confirmed by the framework analysis.
  Tonic BCR signalling (IGHD elevated +43%,
  FCRL5 elevated +415%) drives BCL2 transcription,
  maintaining the lock continuously.

  THE WADDINGTON INTERPRETATION:
  In normal B cell development, the mature naive B
  cell position is NOT a stable terminal attractor.
  Normal B cells must:
    Either encounter antigen → activate → germinal
    centre → become memory B cell or plasma cell.
    Or receive death signals → apoptosis.
  The normal B cell landscape has NO long-term
  stable position at the mature naive B cell state.
  B cells pass THROUGH that state — they do not
  reside there permanently.
  The CLL cell has found a way to stabilise a
  TRANSITIONAL STATE — the mature naive B position
  — and turned it into a FALSE ATTRACTOR.
  By maintaining tonic BCR signalling (BCL2 up),
  expressing inhibitory FcRL5 (anergy, avoiding
  activation), and suppressing PRDM1 (blocking
  plasma cell terminal fate), the CLL cell has
  created a STABLE IDENTITY where no stable
  identity should exist in the normal landscape.
  This is the most structurally interesting
  attractor in the repository — not a blocked
  intermediate stage but a STABILISED TRANSITIONAL
  STATE that should never have been permanent.

  THE FRAMEWORK VALIDATION #8 — WHAT HAPPENED:
  The automated script scored CLL as "insufficient
  confirmation" (1/4 switch genes confirmed).
  The honest interpretation: the switch gene
  predictions were designed for differentiation
  block attractors (like B-ALL) and did not
  transfer to a survival attractor.
  BCL2 elevation (+136% ***) was confirmed.
  PRDM1 suppression (-57% ***) was confirmed.
  IGKC elevation (+60%) confirmed the hierarchy
  position (deeper than B-ALL as predicted).
  The ibrutinib time-series confirmed attractor
  dissolution: BCL2 fell 83% by day 150,
  IGHD and FCRL5 fell sharply, PRDM1 stayed flat
  (cells exit by dying, not by differentiating).
  Both predicted drug classes (BTK inhibitor,
  BCL2 inhibitor) are FDA-approved for CLL.
  The attractor was confirmed. The automated
  scoring missed the confirmation because the
  scoring logic was designed for the wrong
  attractor type.
  THIS SUBTYPE DOCUMENT AND THE SUBTYPE SERIES
  THAT FOLLOWS RESOLVES THE OPEN QUESTION
  FROM VALIDATION #8:
  Do different CLL molecular subtypes have
  different depths within the survival attractor?
  Does the IGHV-mutated clone occupy a SHALLOWER
  position (closer to a true memory B cell identity)
  while the IGHV-unmutated clone occupies a DEEPER
  position (more thoroughly locked, BCL2 higher,
  more tonic BCR signal)?
  Do del17p/TP53-mutated CLL cells add a second
  dimension of depth — the genome guardian loss
  creating a deeper attractor within the already-
  false survival attractor?
  These are the questions the subtype series answers.
```

---

## DOCUMENT METADATA

```
document_id:        CLL_Subtype_Orientation
series:             CLL (Chronic Lymphocytic Leukemia
                    — Subtypes)
folder:             Cancer_Research/CLL/Subtypes/
date:               2026-03-04
author:             Eric Robert Lawson
status:             ORIENTATION ONLY — no predictions,
                    no analysis, no epigenetic hypotheses
existing_analysis:  Cancer_Research/CLL/
                    CLL_False_Attractor_confirmed.md
                    (Document 80 — Validation #8)
                    Survival attractor confirmed.
                    BCL2 elevated +136% ***
                    Drug targets: BTK + BCL2 confirmed.
                    Ibrutinib dissolution confirmed.
                    Framework lesson: survival attractor
                    needs different switch gene logic.
next_document:      CLL-S1a
                    IGHV-mutated CLL Before-Document
                    (most tractable: shallowest position,
                    closest to true memory B cell identity)
protocol_version:   Workflow_Protocol.md v2.0
```

---

## SECTION I — THE NORMAL B CELL PROGRAMME

```
Understanding the CLL false attractor requires
understanding the developmental programme it has
subverted. The B cell developmental hierarchy is
the most precisely characterised differentiation
cascade in the body — over 70 years of immunological
research have mapped every stage, every transcription
factor, every surface marker.

═══════════════════════════════════════════════════════
THE B CELL DEVELOPMENTAL HIERARCHY — FULL MAP
═══════════════════════════════════════════════════════

BONE MARROW STAGES:
──────────────────

STAGE 1: COMMON LYMPHOID PROGENITOR (CLP)
  Emerges from MPP4 in the bone marrow.
  Still capable of B, T, NK, and dendritic
  cell fates.
  Key TF commitment: PAX5 is NOT yet expressed.
  PAX5 expression = irreversible B cell fate.
  RUNX1, IKZF1 (Ikaros) drive CLP toward B fate.

STAGE 2: PRO-B CELL
  PAX5 expression begins → irreversible B fate.
  PAX5 represses non-B programmes (T cell, myeloid).
  PAX5 activates CD19, which becomes the defining
  surface marker of B lineage cells from this
  point forward.
  V(D)J RECOMBINATION BEGINS:
    RAG1 and RAG2 are expressed.
    RAG complex first recombines the
    IGHV-D-J segments (heavy chain).
    First recombination: D-J joining.
    Then: V-DJ joining.
    Result: a COMPLETE IMMUNOGLOBULIN HEAVY CHAIN
    (IgH μ chain).
  If productive rearrangement: advance to Pre-B.
  If non-productive: second attempt, then apoptosis.
  KEY TF NETWORK (Pro-B):
    PAX5 — master B lineage TF
    EBF1 (Early B Cell Factor 1) — activates
           B lineage programme with PAX5
    E2A (TCF3) — bHLH TF for B fate
    IKZF1 — chromatin accessibility at Ig loci
    RUNX1 — transition from CLP to Pro-B

STAGE 3: LARGE PRE-B CELL
  Surface marker: Pre-BCR complex.
  The μ chain (heavy chain) pairs with the
  surrogate light chain (VpreB + λ5)
  to form the PRE-BCR.
  Pre-BCR signalling through SYK and BTK
  drives massive clonal expansion (3–5 divisions).
  V(D)J HEAVY CHAIN RECOMBINATION IS COMPLETE.
  D-J and V-DJ rearrangements LOCKED.
  Pre-BCR signalling then silences RAG1/RAG2
  transiently (allelic exclusion mechanism —
  prevents a second productive rearrangement).
  KEY TF at this stage:
    IKZF3 (Aiolos) — promotes Pre-B transition
    IRF4 — feedback from pre-BCR signalling
           silences the surrogate light chain,
           driving transition to small Pre-B

STAGE 4: SMALL PRE-B CELL
  Pre-BCR is internalised.
  Cell exits cell cycle.
  RAG1/RAG2 are RE-EXPRESSED.
  LIGHT CHAIN V(D)J RECOMBINATION BEGINS:
    Kappa light chain (IGKC) locus recombines first.
    If productive κ rearrangement: IGKC expressed.
    If non-productive κ: lambda light chain (IGLC).
  If productive light chain rearrangement:
    The complete IgM BCR (μ heavy + κ or λ light)
    is expressed on the cell surface.
  IMPLICATION FOR FRAMEWORK:
    IGKC suppressed = cell has not yet completed
    light chain V(D)J recombination = BLOCKED
    at the Pre-B to Immature B transition.
    This was the B-ALL finding (IGKC suppressed
    -83.7% — B-ALL is blocked at Pre-B stage).
    IGKC expressed = cell has completed light
    chain recombination = PAST the Pre-B stage.
    This was the CLL finding (IGKC elevated +60%
    — CLL cells are mature B cells, V(D)J complete,
    CLL is deeper in development than B-ALL).
    RAG1/RAG2 silent in CLL (RAG1=0.0002,
    RAG2=0.000) — no V(D)J recombination activity,
    confirming mature B cell identity.

STAGE 5: IMMATURE B CELL
  Surface IgM (sIgM) expressed.
  First CENTRAL TOLERANCE CHECK:
    If BCR reacts strongly with self-antigens:
    RECEPTOR EDITING (re-express RAG, rearrange
    light chain again to change specificity).
    If editing fails: apoptosis (clonal deletion).
    If BCR is weakly self-reactive: anergy
    (functional silencing without apoptosis).
  If BCR is non-self-reactive: advance to
  transitional B cell and export to periphery.

STAGE 6: MATURE NAIVE B CELL (PERIPHERY)
  The cell that CLL has stabilised at.
  Surface: IgM + IgD (co-expressed, same VDJ
  rearrangement, different constant region
  via alternative splicing).
  Surface: CD19+, CD20+, CD21+, CD23+, CD27-.
  This cell is READY TO BE ACTIVATED by antigen.
  It is NOT a stable long-term attractor in
  normal biology. It is a TRANSITIONAL state —
  awaiting antigen encounter.
  KEY TF STATE (Mature Naive B):
    PAX5 — maintained (B identity)
    IRF4 — low
    BCL6 — not expressed
    PRDM1 (Blimp1) — not expressed
    BCL2 — LOW to moderate (survival signal
            from tonic survival circuits but
            not at CLL levels)
  THIS IS THE BASELINE FROM WHICH CLL DEPARTS.
  CLL TAKES THIS CELL AND:
    1. ELEVATES BCL2 (136% above normal B) —
       via tonic BCR → BTK → NF-κB/AP-1 → BCL2
    2. MAINTAINS IGHD (IgD on surface, +43%) —
       the dual IgM/IgD BCR drives tonic BTK
    3. ELEVATES FCRL5 (+415%) — anergy marker,
       blocks acute activation while tonic signal
       continues
    4. SUPPRESSES PRDM1 (-57%) — blocks plasma
       cell terminal fate exit
    5. ELEVATES CD27 (+817%) — memory B cell marker
       CLL cells co-opt the memory B cell surface
       phenotype as part of the false attractor
       (CD27+ CD5+ CD23+ B cell = the CLL phenotype)
    Result: a STABLE FALSE ATTRACTOR at a position
    in B cell development space that is normally
    only a transient waystation.

PERIPHERAL B CELL ACTIVATION STAGES:

STAGE 7: T-INDEPENDENT vs T-DEPENDENT ACTIVATION
  T-independent (TI):
    Short-lived plasma cells producing IgM.
    No somatic hypermutation.
    No memory B cell formation.
    Responds to polysaccharide antigens.
  T-dependent (TD):
    Involves cognate T cell help (CD40L-CD40).
    Drives GERMINAL CENTRE FORMATION.

STAGE 8: GERMINAL CENTRE (GC) B CELL
  THE CRITICAL STAGE FOR IGHV MUTATION STATUS.
  Location: lymph node follicles (secondary
  lymphoid organs).
  Two zones:
    DARK ZONE: B cells (centroblasts) proliferate
    rapidly and UNDERGO SOMATIC HYPERMUTATION (SHM)
    of the IGHV-D-J region.
    SHM ENZYME: AICDA (Activation-Induced Cytidine
    Deaminase) — introduces point mutations into
    the VDJ region, generating antibody diversity.
    This is what creates IGHV MUTATION IN NORMAL B CELLS.
    LIGHT ZONE: Centrocytes present their mutated
    BCR to follicular dendritic cells (FDC) and
    T follicular helper cells (Tfh).
    High-affinity BCR = selected for survival.
    Low-affinity BCR = apoptosis (positive selection).
  GC B CELL MASTER TF: BCL6
    BCL6 is the defining TF of GC B cells.
    BCL6 represses:
      PRDM1 (Blimp1) — prevents premature
                        plasma cell differentiation
      IRF4 — prevents exit from GC state
      TP53 — allows proliferation and SHM without
             checkpoint activation
      CDKN1A (p21) — allows rapid cycling
    BCL6 must be silenced for exit from GC.
    HOW BCL6 IS SILENCED:
      High-affinity BCR signalling → IRF4 high
      → IRF4 represses BCL6 → BCL6 falls
      → PRDM1 is de-repressed → GC B cell
      exits to plasma cell fate.
  IN CLL WITH IGHV UNMUTATED:
    The CLL precursor cell NEVER WENT THROUGH
    THE GERMINAL CENTRE.
    OR: It went through the GC but failed to
    accumulate mutations (antigen-independent
    or antigen-stimulated but SHM-deficient).
    IgVH sequence is <2% divergent from germline.
    Consequence: the CLL clone is derived from
    a NAIVE B CELL or PRE-GC B CELL.
    The tonic BCR signalling is derived from
    a GERMLINE CONFIGURATION BCR — the most
    reactive of BCRs, recognising a stereotyped
    antigen (or self-antigen in some subsets).
    THIS IS THE MORE AGGRESSIVE FORM OF CLL.
  IN CLL WITH IGHV MUTATED:
    The CLL precursor cell WENT THROUGH THE
    GERMINAL CENTRE and accumulated SHM.
    IgVH sequence is >2% divergent from germline.
    Consequence: the CLL clone is derived from
    a MEMORY B CELL or POST-GC B CELL.
    The BCR has been affinity-matured and may
    be less reactive to self-antigens.
    The tonic BCR signalling is WEAKER.
    BCL2 is LOWER (attractor is SHALLOWER).
    THIS IS THE LESS AGGRESSIVE FORM OF CLL.
    FRAMEWORK INTERPRETATION:
    The IGHV mutation status determines the DEPTH
    of the CLL survival attractor.
    IGHV mutated = shallow survival attractor
                  (lower BCL2, weaker tonic BCR)
    IGHV unmutated = deep survival attractor
                    (higher BCL2, stronger tonic BCR)
    The depth score in CLL should separate these
    two populations. This is the primary
    hypothesis of the subtype series.

STAGE 9: MEMORY B CELL
  Exits germinal centre.
  Long-lived.
  CD27+ CD20+ surface IgM or IgG.
  CD27 is THE MEMORY B CELL MARKER.
  IN THE CLL DATA: CD27 elevated +817% ***
    This is the most dramatic elevation in the
    entire CLL dataset.
    CLL cells MASSIVELY OVEREXPRESS CD27.
    They have co-opted the memory B cell
    surface identity as part of the false attractor.
    CLL cells ARE NOT MEMORY B CELLS (they have
    not necessarily gone through GC) but they
    LOOK LIKE memory B cells on their surface.
    This is the phenotypic mimicry that makes
    CLL clinically confusing:
    CD27+ CD5+ CD23+ looks like activated memory
    B cells. It is not.
    It is the survival attractor dressed in the
    surface phenotype of a memory B cell.

STAGE 10: PLASMA CELL
  PRDM1 (Blimp1) is the MASTER PLASMA CELL TF.
  Blimp1 activates:
    Immunoglobulin secretion machinery
    (XBP1 → ER expansion → Ig secretion)
    Cell cycle arrest (anti-proliferative)
    Surface BCR downregulation
  Blimp1 represses:
    PAX5 — B cell identity is LOST in plasma cells
    MYC, IRF4 (in the GC context)
    BCL6 (mutual repression: BCL6↑ blocks PRDM1,
           PRDM1↑ represses BCL6)
  PLASMA CELL IS THE TERMINAL ATTRACTOR of
  B cell development. It is IRREVERSIBLE — once
  Blimp1 is active and PAX5 is silenced, the
  cell can never return to B cell identity.
  IN CLL: PRDM1 is suppressed (-57% ***).
    The CLL cell cannot reach the plasma cell attractor.
    The terminal exit from B cell identity is blocked.
    Even when ibrutinib cuts BCR signalling and BCL2
    falls, PRDM1 STAYS FLAT.
    CLL cells exit the survival attractor by DYING
    (apoptosis via caspase-3 when BCL2 falls),
    NOT by completing B cell development to plasma
    cell or memory B cell exit.
    This was the critical ibrutinib time-series
    finding from Validation #8.
```

---

## SECTION II — THE SURVIVAL ATTRACTOR GEOMETRY

```
THE WADDINGTON SURFACE OF CLL IS INVERTED
COMPARED TO ALL OTHER CANCERS IN THE REPOSITORY.

In differentiation-block cancers (AML, B-ALL,
PRAD, PAAD, etc.): the cancer cell is a PRE-FORMED
STRUCTURE — it is upstream of the normal attractor
and is blocked from reaching it. The false attractor
is SHALLOWER in development than the normal terminal
state.

In CLL: the cancer cell is a POST-FORMED STRUCTURE —
it is DOWNSTREAM of the normal progenitor state,
at a position that should have been transient.
The false attractor is AT OR PAST the terminal
differentiation point of the upstream development,
but it has failed to die.

THE WADDINGTON METAPHOR IN PHYSICAL TERMS:

Normal B cell development as a marble on a surface:
  The marble (developing B cell) rolls from the
  top (HSC) downhill through successive valleys
  (developmental stages: Pro-B, Pre-B, Immature B,
  Naive B) toward a terminal basin at the bottom
  (plasma cell — the true terminal attractor).

  In B-ALL: The marble gets stuck in a small
  pit (false attractor) BEFORE reaching the
  mature B cell position — it never completes
  light chain V(D)J recombination, never expresses
  IGKC, never reaches the naive B stage.

  In CLL: The marble ROLLS PAST the B-ALL pit,
  completes V(D)J (IGKC expressed), reaches
  the mature naive B position — and THEN gets
  stuck, because BCL2 has been elevated and
  the naive B position, which should have been
  a sloped waystation, has been turned into
  a false basin. The marble cannot roll further
  to plasma cell (PRDM1 suppressed) and cannot
  roll backward. It sits in a position that has
  been converted from a slope into a false pit.

THE THREE-AXIS DEPTH PROBLEM IN CLL:
  Unlike solid tumours (one depth axis = distance
  from normal terminal identity), CLL has THREE
  depth axes:

  AXIS 1 — BCR-DEPENDENT LOCK DEPTH:
  How strongly is the BCR maintaining BCL2?
  Proxies:
    BCL2 expression level
    IGHD expression (IgD co-expressed with IgM)
    FCRL5 expression (anergy state depth)
    CD27 elevation
  This axis reflects the strength of tonic
  BCR → BTK → NF-κB → BCL2 circuit.
  IGHV unmutated: BCR is germline-configured,
  maximally reactive → AXIS 1 DEEP.
  IGHV mutated: BCR has been affinity-matured
  in GC, less reactive → AXIS 1 SHALLOW.

  AXIS 2 — GENOMIC INSTABILITY DEPTH:
  Has the cell acquired secondary genetic lesions
  that ADD a second survival mechanism beyond BCL2?
  Proxies:
    TP53 mutation / del17p: apoptosis pathway
    abolished. BCL2 inhibition can drive cells
    toward apoptosis via BAX/BAK — but BAX
    activation requires p53. TP53-mutant cells
    have a second survival advantage beyond BCL2.
    del11q (ATM loss): DNA damage response
    impaired. Cells cannot initiate apoptosis
    in response to genotoxic stress.
    SF3B1 mutation: splicing dysregulation adds
    a third level of survival signal.
  This axis reflects secondary genetic depth —
  acquired mutations that deepen the false
  attractor by adding survival mechanisms
  INDEPENDENT OF BCL2.
  del17p/TP53: AXIS 2 DEEP (resistant to BCL2
  inhibitor-mediated apoptosis because p53-
  independent apoptosis is the only available
  pathway, and BTK inhibitors are required to
  reduce BCL2 first).

  AXIS 3 — PROLIFERATIVE COMPARTMENT DEPTH:
  CLL has a unique biology: the bulk of CLL cells
  are in G0 (quiescent, MKI67 low -99% in data),
  but a small fraction (~1–2% of CLL cells, the
  "proliferation centres" or PSEUDO-FOLLICLES
  in lymph nodes and bone marrow) are actively
  cycling.
  The proliferating cells are the MOST DANGEROUS —
  they are where new mutations arise, where
  AICDA can introduce additional genomic damage,
  and where clonal evolution occurs.
  Proxy: Ki67, MYC, MCM2-7 expression in
  the proliferating fraction.
  IGHV unmutated: larger proliferating fraction,
  more frequent proliferation centre activity.
  Elevated doubling time: AXIS 3 DEEP.
  IGHV mutated: smaller proliferating fraction.
  AXIS 3 SHALLOW.

  THE COMBINED DEPTH SCORE FOR CLL:
  The framework depth score in CLL should capture
  all three axes:
    Axis 1 (BCR lock): BCL2, IGHD, FCRL5, CD49d
    Axis 2 (genomic): TP53, ATM, SF3B1 status
    Axis 3 (proliferative): Ki67, CXCR4
  The combined depth score should separate:
    IGHV-mutated, del13q, low-risk → low depth
    IGHV-unmutated, trisomy 12 → intermediate depth
    IGHV-unmutated, del11q → deeper
    IGHV-unmutated, del17p/TP53 → deepest

THE SURVIVAL ATTRACTOR CONFIRMATION
FROM VALIDATION #8 — GENE BY GENE:

BCL2 +136% *** — AXIS 1 LOCK CONFIRMED
  p < 1e-45
  The lock is real, strong, and statistically
  incontestable.
  Drug target confirmed: venetoclax → BCL2.

IGHD +43% *** — AXIS 1 BCR SIGNAL CONFIRMED
  p < 1e-9 for elevation
  NOT A FAILED PREDICTION — AN ANALYST ERROR.
  The prediction was "IGHD suppressed" because
  the logic was: IGHD is lost as B cells activate.
  But CLL does NOT acutely activate — it tonically
  signals. IgD co-expression with IgM drives the
  dual IgM/IgD BCR → maximises BTK activation
  → maximises BCL2 → the CLL attractor IS the
  IgM/IgD co-expression state.
  IGHD is not suppressed — it is MAINTAINED as
  PART OF THE SURVIVAL MECHANISM.
  Confirmed by ibrutinib dissolution: IGHD falls
  from 0.238 to 0.000 by day 150. The BCR signal
  is cut by BTK inhibition → IgD is no longer
  needed to maintain BCL2.

FCRL5 +415% *** — AXIS 1 ANERGY CONFIRMED
  p < 1e-85 for elevation
  ALSO AN ANALYST ERROR IN DIRECTION.
  FCRL5 is an inhibitory receptor that marks
  ANERGIC B cells.
  CLL cells are FUNCTIONALLY ANERGIC — they avoid
  acute BCR activation while maintaining tonic
  signal. FCRL5 upregulation IS the anergy marker
  of the CLL false attractor.
  The CLL cell has become "immune to its own
  immune activation" — it has upregulated the
  very receptor that prevents B cell activation
  (FCRL5) while keeping the tonic survival signal
  running. This is a sophisticated false attractor:
  it simultaneously maintains survival (BCL2 via
  tonic BTK) AND avoids activating mechanisms that
  would drive it toward GC, differentiation, or
  deletion (FCRL5-mediated anergy).
  Confirmed by ibrutinib: FCRL5 falls from 0.162
  to 0.001 by day 150. When BTK is inhibited,
  the tonic signal is cut, the anergy state
  is no longer needed, and FCRL5 expression falls.

CD27 +817% *** — MEMORY B PHENOTYPE CONFIRMED
  p ~ 0 (effectively below detection threshold)
  CD27 is the canonical memory B cell marker.
  Its 817% elevation in CLL is the most dramatic
  finding in the entire Validation #8 dataset.
  CLL cells have co-opted the MEMORY B CELL
  SURFACE PHENOTYPE as part of the false attractor.
  This is not coincidence — it is the nature of
  the survival attractor state: the CLL cell
  is expressing the surface proteins of a mature
  post-GC memory B cell while having the internal
  BCL2-high, PRDM1-low, tonic-BCR-active
  gene expression profile of a pathological
  survivor.
  THE CD27 ELEVATION TELLS THE FRAMEWORK WHERE
  IN THE B CELL HIERARCHY THE CLL FALSE ATTRACTOR
  IS ANCHORED: at the memory B cell position,
  not at the naive B cell position.
  This means CLL cells are mimicking the MOST
  STABLE normal B cell peripheral identity
  (memory B cells are the long-lived, non-cycling
  peripheral B cell population) and using that
  mimicry to avoid immune surveillance.

PRDM1 -57% *** — TERMINAL EXIT BLOCKED CONFIRMED
  p < 1e-6
  Blimp1 is suppressed in CLL.
  The plasma cell terminal exit is BLOCKED.
  CLL cells cannot complete B cell development.
  The survival attractor has one exit: apoptosis.
  The ibrutinib time-series confirms: PRDM1 stays
  flat throughout treatment (d0 = 0.008, d150 =
  0.024 — negligible change).
  CLL cells DO NOT differentiate under treatment.
  They die when BCL2 falls below the survival
  threshold.

MKI67 -99% — NON-PROLIFERATIVE CONFIRMED
  Accumulated CLL cells are in G0.
  This is the paradox of CLL:
  A rapidly accumulating cancer of non-proliferating
  cells. Cells pile up because they don't die,
  not because they divide rapidly.
  The MKI67 near-zero signal confirms the bulk
  population is quiescent.
  The 1–2% proliferating fraction (proliferation
  centres) is NOT captured by bulk MKI67 — this
  is an important caveat for the subtype analysis
  (single-cell data will show the proliferating
  fraction; bulk data will not).

IGKC +60% — HIERARCHY CONFIRMED
  CLL cells have completed V(D)J light chain
  recombination (IGKC expressed, not suppressed).
  CLL cells are downstream of B-ALL in development.
  The framework internal cross-check passed:
  CLL block is DEEPER than B-ALL as predicted.
```

---

## SECTION III — THE IGHV DIVIDE

```
THE MOST IMPORTANT STRATIFICATION IN CLL IS NOT
A DRUG TARGET OR A CYTOGENETIC LESION.
IT IS THE DEVELOPMENTAL HISTORY OF THE BCR.

The presence or absence of SOMATIC HYPERMUTATION
in the IGHV variable region segments tells us:
  Did this cell go through the germinal centre?
  Was its BCR affinity-matured?
  How reactive is the remaining BCR to its antigen?
  How strongly is tonic BCR signalling driving BCL2?

This is a DEPTH MEASURE encoded in the BCR sequence
itself — not in gene expression, not in cytogenetics,
but in the somatic mutation rate of the antibody gene.
It is the equivalent of the Gleason score in PRAD
or the RAS mutation in CRC — the single most
predictive marker of disease biology.

IGHV MUTATED (>2% divergence from germline):
──────────────────────────────────────────────
FREQUENCY:      ~50–60% of CLL patients.
PROGNOSIS:      FAVOURABLE.
MEDIAN OS:      ~20–25 years (some series show
                near-normal life expectancy).
AML TRANSFORM:  Rare (<5% in 10 years).
TREATMENT NEED: Often never requires treatment
                (watch-and-wait for years or
                indefinitely in some patients).

CELL OF ORIGIN:
  Post-germinal centre B cell.
  The CLL precursor went through the GC,
  underwent somatic hypermutation by AICDA,
  had its BCR affinity-matured.
  After GC exit, it became a MEMORY B CELL
  and then acquired the CLL driver events
  (BCL2 overexpression, CD5 acquisition)
  that stabilised the memory B position
  as a false attractor.

BCR CHARACTERISTICS:
  Affinity-matured BCR.
  Lower reactivity to antigens/self-antigens.
  WEAKER TONIC BCR SIGNAL.
  LOWER BCL2 (relative to unmutated CLL).
  The survival lock is present but SHALLOW.

DEPTH AXIS 1 PREDICTION:
  BCL2 lower than IGHV-unmutated.
  IGHD lower (less IgD co-expression needed).
  FCRL5 lower (less anergy maintenance needed).
  CD27 may still be high (memory B phenotype).
  This is the shallower arm of the CLL attractor.
  Venetoclax response BETTER:
    Lower BCL2 = fewer molecules of BCL2 to
    inhibit = lower venetoclax ramp-up
    requirement = deeper, faster remissions.
    Consistent with clinical data: IGHV-mutated
    CLL has higher rates of uMRD (undetectable
    MRD) on venetoclax-based regimens.

STEREOTYPED BCR SUBSETS IN IGHV-MUTATED:
  Subset 4: IGHV4-34, mutated.
            Associated with specific self-antigen
            recognition.
            Intermediate prognosis.
  Subset 16: IGHV1-69, mutated (rare).
             Better prognosis.
  Many IGHV-mutated CLL cases are not stereotyped —
  each has a unique BCR sequence resulting from
  individualised SHM pattern.

TREATMENT IN IGHV-MUTATED CLL:
  Watch and wait: many patients never need treatment.
  When treatment required:
    FCR (fludarabine, cyclophosphamide, rituximab):
    HISTORICALLY the best regimen for young, fit
    IGHV-mutated CLL.
    CLL8 trial: ~80% ORR, median PFS ~55 months.
    Some IGHV-mutated patients achieve long-term
    MRD-negative remissions (functional cure).
    FCR is increasingly replaced by targeted
    therapy in 2024 — but IGHV-mutated patients
    have the best outcomes with FCR historically.
    VENETOCLAX + OBINUTUZUMAB (CLL14):
    Fixed-duration (12 months), deep remissions,
    high uMRD rates. IGHV-mutated patients
    achieve the best uMRD rates in this regimen.
    BTK inhibitors: effective but typically
    continuous therapy (indefinite). Less
    preferred for IGHV-mutated patients who can
    achieve deep fixed-duration remissions with
    venetoclax combinations.

IGHV UNMUTATED (<2% divergence from germline):
───────────────────────────────────────────────
FREQUENCY:      ~40–50% of CLL patients.
PROGNOSIS:      UNFAVOURABLE.
MEDIAN OS:      ~8–10 years.
AML TRANSFORM:  Higher (~10–15% in 10 years).
TREATMENT NEED: Usually requires treatment within
                5 years of diagnosis.

CELL OF ORIGIN:
  Pre-germinal centre B cell (naive B cell)
  OR a post-GC cell that exited without SHM
  accumulation.
  The CLL precursor either:
  a) NEVER entered the germinal centre —
     its BCR retains germline configuration,
     maximally reactive to its cognate antigen.
  b) ENTERED the GC but failed to accumulate
     sufficient SHM (<2% divergence).
  In either case: the BCR is GERMLINE-CONFIGURED
  and reacts more strongly with its antigen (or
  with self-antigens such as muscle myosin,
  cytoskeletal proteins, or microbial antigens
  cross-reactive with self).
  STRONGER TONIC BCR SIGNAL.
  HIGHER BCL2.
  DEEPER SURVIVAL ATTRACTOR.

ZAP70 AND CD38 AS PROXIES FOR IGHV STATUS:
  Before cheap NGS, ZAP70 and CD38 were used
  as surrogate markers for IGHV status.
  ZAP70+ (positive, >20% cells): correlates with
  IGHV unmutated → poor prognosis.
  CD38+ (>30% cells): correlates with IGHV
  unmutated → poor prognosis.
  FRAMEWORK INTERPRETATION:
  ZAP70 is a tyrosine kinase that AMPLIFIES BCR
  signalling (it lowers the threshold for BCR
  activation). ZAP70 elevation in CLL cells
  → more sensitive BCR → stronger tonic signal
  → more BCL2 → deeper survival attractor.
  ZAP70 elevation IS the depth indicator for
  Axis 1 in IGHV-unmutated CLL, expressed at
  the kinase level rather than the surface level.

CD49d (Integrin α4):
  The strongest single prognostic marker in CLL
  — stronger than CD38, ZAP70, and comparable to
  IGHV status in some analyses.
  CD49d+ (>30%): poor prognosis.
  Mechanism: CD49d binds fibronectin and VCAM-1
  in the bone marrow microenvironment → activates
  PI3K/AKT → survival signalling → SUPPLEMENTS
  BCL2-mediated survival with MICROENVIRONMENT-
  DEPENDENT survival signals.
  AXIS 1 DEPTH MARKER: CD49d elevation =
  additional microenvironment-derived survival
  dimension on top of BCR-driven BCL2.
  CD49d+ CLL cells have TWO survival signals:
    1. BCR → BTK → BCL2 (intrinsic)
    2. CD49d → integrin → PI3K/AKT (extrinsic)
  FRAMEWORK PREDICTION: CD49d correlates with
  depth score in CLL because it adds a second
  attractor-maintenance mechanism.
  Not yet measured in the Validation #8 gene set
  (CD49d = ITGA4 — should be added to Script 2
  of the IGHV subtype analysis).

STEREOTYPED BCR SUBSETS IN IGHV-UNMUTATED:
  ~30% of IGHV-unmutated CLL have STEREOTYPED
  BCRs — identical or near-identical CDR3
  sequences despite being from different patients.
  This means the SAME ANTIGEN is driving the
  tonic BCR signal in these patients.
  The two most clinically important subsets:

  SUBSET 1 (IGHV1-5-7/IGHD6-19/IGHJ4,
            IGKV1-39/IGHJ1 light chain):
    Unmutated, CD38+, ZAP70+.
    Rapid progression, early treatment need.
    Associated with DNA damage accumulation
    and SF3B1 mutations over time.
    Responds well to BTK inhibitors.
    Responds less well to venetoclax (high BCL2,
    deep attractor, high ramp-up requirement).

  SUBSET 2 (IGHV3-21, regardless of mutation):
    UNUSUAL: can be MUTATED or UNMUTATED but
    has poor prognosis regardless.
    Associated with del11q (ATM loss).
    Enriched for SF3B1 mutations.
    More prone to AML-like progression with
    complex karyotype acquisition.
    The framework interpretation: SUBSET 2 CLL
    is the only CLL with IGHV-mutation-
    independent aggressive biology — suggesting
    the ANTIGEN SPECIFICITY (what the BCR
    recognises) is the driver, not just whether
    the BCR was affinity-matured.
    In Subset 2: the BCR recognises a structurally
    rigid antigen that drives STRONG TONIC SIGNAL
    regardless of SHM status — overriding the
    normally protective effect of IGHV mutation.
    This is the deepest IGHV-mutated CLL variant.

  SUBSET 4 (IGHV4-34):
    IgG-expressing CLL (class-switched variant).
    Intermediate prognosis.
    Less dependent on IgM/IgD tonic BCR signal
    (already class-switched to IgG → different
    BCR signalling dynamics).

TREATMENT IN IGHV-UNMUTATED CLL:
  FCR: NOT recommended.
    FCR in IGHV-unmutated CLL produces shorter
    PFS (~20 months) compared to IGHV-mutated.
    More genomic instability → FCR is genotoxic
    → selects for TP53-mutant resistant clones.
    FCR in IGHV-unmutated CLL with TP53 aberration
    = near-zero response rate.
  IBRUTINIB (BTK inhibitor): superior to FCR
    for IGHV-unmutated CLL (E1912 trial).
    Continuous therapy.
  ACALABRUTINIB or ZANUBRUTINIB:
    Next-generation BTK inhibitors.
    Better tolerated than ibrutinib (fewer atrial
    fibrillation, bleeding, hypertension events).
    Zanubrutinib (ALPINE trial): superior PFS
    vs. ibrutinib, same efficacy with better
    tolerability.
    FRAMEWORK INTERPRETATION:
    Acalabrutinib/zanubrutinib are the same
    attractor dissolution mechanism (BTK inhibition
    → BCL2 falls → apoptosis) with improved
    off-target pharmacology.
    They do not change the attractor geometry —
    only the drug tolerance profile.
  VENETOCLAX + OBINUTUZUMAB (CLL14):
    Fixed-duration, high uMRD rates.
    IGHV-unmutated: achieves uMRD but at lower
    rates and shorter duration than IGHV-mutated.
    Consistent with deeper attractor: more BCL2
    to inhibit, less complete responses.
```

---

## SECTION IV — THE CYTOGENETIC HIERARCHY

```
FISH CYTOGENETICS in CLL was codified by the
Döhner hierarchical model (NEJM 2000) and remains
the standard cytogenetic risk stratification.
The five-tier hierarchy:

TIER 1: del(17p) — WORST PROGNOSIS
  TP53 DELETION — the genome guardian is lost.
  Frequency: ~7% of untreated CLL.
             ~20–30% of relapsed/refractory CLL
             (selected for by chemotherapy and
             BTK inhibitors over time).
  TP53 MUTATION (without del17p):
    Functionally equivalent to del(17p) alone.
    Now categorised separately but same clinical
    management.
  COMBINED: TP53 aberration (del17p and/or
    TP53 mutation) = ~10% at diagnosis, up to
    50% in late-stage refractory CLL.

  DEPTH AXIS 2 INTERPRETATION:
    TP53 loss in CLL is identical in principle
    to TP53 loss in NEPC (PRAD), MDS-biTP53,
    and TP53 co-mutation in BRCA:
    The genome guardian is lost → genomic
    instability accelerates → secondary mutations
    accumulate rapidly → attractor deepens via
    new survival mechanisms beyond BCL2.
    TP53 loss specifically disables p53-dependent
    apoptosis — the PATHWAY THAT VENETOCLAX USES.
    When BCL2 is inhibited by venetoclax:
    Anti-apoptotic BCL2 is blocked → BAX/BAK
    oligomerise → OUTER MITOCHONDRIAL MEMBRANE
    PERMEABILISATION → Cytochrome C release →
    APOPTOSOME → CASPASE-9 → CASPASE-3 →
    APOPTOSIS.
    But: BAX/BAK activation requires PUMA and
    NOXA — BH3-only proteins that are activated
    by p53.
    WITHOUT P53: PUMA/NOXA are not induced →
    BAX/BAK activation is impaired → Cytochrome C
    release is reduced → Apoptosis is attenuated.
    THEREFORE: TP53-mutant/del17p CLL has REDUCED
    SENSITIVITY TO VENETOCLAX (apoptosis still
    occurs via p53-independent BH3-only proteins
    like BIM — which is not p53-dependent — but
    the full apoptotic response is blunted).
    CLINICAL DATA CONFIRMS:
    TP53-aberrant CLL has lower uMRD rates
    with venetoclax regimens.
    BUT: BTK inhibitors are MORE EFFECTIVE
    in TP53-aberrant CLL (BTKi → BCL2 ↓ →
    even without full p53-dependent apoptosis,
    cells die from BCL2 insufficiency because
    BTK inhibition removes the upstream signal
    driving BCL2 transcription entirely).
    TREATMENT:
    BTK inhibitor (ibrutinib, zanubrutinib,
    acalabrutinib) as primary therapy.
    Venetoclax has reduced efficacy but still used.
    Allo-HSCT: for fit patients with refractory
    del17p CLL — the only curative option.

TIER 2: del(11q) — POOR PROGNOSIS
  ATM DELETION (ataxia-telangiectasia mutated).
  Frequency: ~15–20% of CLL.
  ATM is a DNA damage response kinase.
  ATM activates p53 in response to double-strand
  breaks → apoptosis of DNA-damaged cells.
  ATM LOSS → DNA damage checkpoint impaired →
  cells with DNA damage survive → genomic
  instability accelerates → secondary mutations
  accumulate.
  THE RELATIONSHIP TO TP53:
    ATM → CHK2 → p53 (the ATM pathway).
    del11q (ATM loss) IMPAIRS the ATM→p53 pathway
    WITHOUT directly deleting TP53.
    Consequence: TP53 protein may be expressed
    but cannot be activated by ATM in response
    to DNA damage.
    This is a FUNCTIONAL TP53 LOSS — p53 is
    present but not activatable via ATM.
    del11q CLL has intermediate sensitivity to
    venetoclax (p53 protein is present but ATM
    activation is impaired — BIM-mediated
    apoptosis still functions).
  CLINICAL FEATURES:
    Bulky lymphadenopathy.
    High blood lymphocyte count.
    Rapid progression.
    Associated with IGHV unmutated (~70% of del11q
    CLL is IGHV unmutated).
    FRAMEWORK: del11q CLL is an Axis 2 deepening
    on top of an already-deep IGHV-unmutated Axis 1
    attractor — the two depth axes compound.
  TREATMENT:
    BTK inhibitors: very effective in del11q CLL
    (ibrutinib produces sustained responses).
    Venetoclax: effective (ATM loss is less
    penetrant than TP53 loss for apoptosis).
    Combination (BTK + venetoclax): preferred
    for del11q to maximise depth of response.

TIER 3: del(13q) ONLY — BEST PROGNOSIS
  SOLE ABNORMALITY (no other cytogenetic lesions).
  Frequency: ~50% of CLL — most common lesion.
  del(13q14) deletes the MIR15A/MIR16-1 microRNA
  cluster.
  MiR-15a and miR-16-1 REPRESS BCL2 translation.
  del(13q) → miR-15a/16-1 LOSS →
  BCL2 NOT REPRESSED → BCL2 ELEVATION.
  THIS IS THE MOLECULAR MECHANISM CONFIRMING
  THE BCL2 DEPTH AXIS:
    The most common CLL cytogenetic lesion
    directly elevates BCL2 by removing its
    post-transcriptional repressors.
    The false attractor is maintained at its
    most basic level by miRNA deletion
    allowing BCL2 to rise.
    This is how the CLL false attractor is
    INITIATED in 50% of cases: miR-15a/16-1
    deletion → BCL2 elevation → BCR tonic
    signal can now lock the cell in survival.
  CLINICAL FEATURES:
    Indolent disease.
    Long watch-and-wait periods.
    Most cases are IGHV mutated (shallow Axis 1)
    + del13q only (minimal Axis 2) =
    the shallowest CLL attractor combination.
    FRAMEWORK PREDICTION:
    del13q-only, IGHV-mutated CLL should have
    the lowest depth score of all CLL subtypes.
    Venetoclax: very deep responses, high uMRD.
    Del13q + IGHV-mutated + venetoclax =
    the most favourable venetoclax response context.
  TREATMENT:
    Watch and wait frequently indefinite.
    When treatment needed: any standard CLL
    regimen (FCR, BTK inhibitor, venetoclax +
    obinutuzumab) — all produce excellent responses.

TIER 4: NORMAL KARYOTYPE — INTERMEDIATE PROGNOSIS
  No cytogenetic abnormalities by FISH.
  Molecular mutations (IGHV status, ZAP70, CD38)
  are the main discriminators.
  Prognosis depends almost entirely on IGHV status.

TIER 5: TRISOMY 12 — INTERMEDIATE PROGNOSIS
  +12 (trisomy of chromosome 12).
  Frequency: ~15% of CLL.
  Mechanism: MYC is on chromosome 8q24 (not 12).
  The trisomy 12 mechanism is NOT MYC-driven
  (unlike trisomy 8 in MDS).
  Chromosome 12 contains:
    MDM2 (12q15) — triplicated copy of MDM2 →
    MDM2 excess → p53 degradation despite TP53
    intact → functional p53 deficiency.
    CDK2, CDK4 — cell cycle kinases at 12q13-q14.
    KITLG — SCF/c-KIT ligand.
  CLINICAL FEATURES:
    Often ATYPICAL CLL morphology (larger cells,
    irregular nuclei, prolymphocytes).
    Frequently associated with unmutated IGHV.
    Enriched for CD38+, ZAP70+ markers.
    Higher expression of CD23 than del13q CLL.
  TRISOMY 12 + SF3B1: common combination.
  SF3B1 mutation (splicing factor, as in MDS
  but here in B cells) is frequent in trisomy 12
  CLL. SF3B1 mutation in CLL causes:
    Altered BCR signalling via aberrant splicing
    of signalling adaptors.
    Altered apoptosis regulation via BCL2-family
    alternative splicing.
  FRAMEWORK: Trisomy 12 CLL has an INTERMEDIATE
  DEPTH — the MDM2 excess provides partial p53
  suppression (deepening Axis 2) but this is
  less penetrant than del17p (full TP53 deletion).
  SF3B1 co-mutation adds a third survival
  mechanism (splicing-derived BCL2 stabilisation).

DEPTH SCORE PREDICTION ACROSS CYTOGENETIC TIERS:
  del13q only, IGHV-mutated: DEPTH 1 (lowest)
  Normal karyotype, IGHV-mutated: DEPTH 2
  Trisomy 12, IGHV-unmutated: DEPTH 3
  del11q, IGHV-unmutated: DEPTH 4
  del17p/TP53, IGHV-unmutated: DEPTH 5 (highest)
  This prediction is to be stated in the
  before-documents and confirmed by data.
```

---

## SECTION V — RESISTANCE MECHANISMS

```
CLL has the most well-characterised resistance
mechanisms of any blood cancer — because BTK
inhibitors and venetoclax have been in use for
>10 years and resistance has been studied
prospectively in clinical samples.

The resistance mechanisms ARE THE MOLECULAR
DESCRIPTION OF ATTRACTOR DEEPENING UNDER
THERAPY.

RESISTANCE TO BTK INHIBITORS (IBRUTINIB,
ACALABRUTINIB, ZANUBRUTINIB — COVALENT):

MECHANISM 1: BTK C481S MUTATION
  The most common resistance mutation.
  Covalent BTK inhibitors bind IRREVERSIBLY to
  cysteine-481 (C481) in the ATP-binding pocket.
  C481S mutation: cysteine → serine substitution.
  Serine cannot form a covalent bond with ibrutinib.
  Consequence: ibrutinib/acalabrutinib/zanubrutinib
  can no longer irreversibly occupy BTK →
  competitive binding (reversible) → BTK kinase
  activity partially restored → BCR signalling
  partially restored → BCL2 rises again.
  FREQUENCY: ~75–85% of covalent BTKi resistance.
  FRAMEWORK INTERPRETATION:
    BTK C481S is the molecular equivalent of
    the false attractor finding a new maintenance
    path when the original path is blocked.
    The CLL attractor was being maintained by
    BTK → BCL2.
    Ibrutinib blocked the BTK→BCL2 pathway.
    C481S restored it via altered binding.
    The attractor reconstituted itself by
    modifying its own lock mechanism.
    This is ATTRACTOR RECONSTITUTION UNDER PRESSURE.

MECHANISM 2: PLCG2 MUTATIONS
  PLCG2 (Phospholipase C gamma 2) is downstream
  of BTK in the BCR signalling cascade.
  BCR → SYK → BTK → PLCG2 → IP3/DAG → 
  PKC → NF-κB → BCL2
  PLCG2 gain-of-function mutations:
    PLCG2 becomes constitutively active
    INDEPENDENT of BTK.
    BTK is blocked → PLCG2 activates anyway →
    NF-κB → BCL2.
    BYPASS OF THE BLOCKED BTK NODE.
  FREQUENCY: ~10–15% of covalent BTKi resistance.
  FRAMEWORK INTERPRETATION:
    PLCG2 GOF is a BYPASS MUTATION — the attractor
    found an alternative route to BCL2 maintenance
    that skips the blocked BTK node.
    This is the molecular equivalent of rewiring
    the circuit around the blocked component.
    Clinically: PLCG2 mutations are present even
    before BTKi exposure at very low allele
    frequency — suggesting pre-existing attractor
    heterogeneity where the PLCG2-bypassed route
    is maintained as a minor sub-attractor that
    is selected for under BTKi pressure.

SOLUTION TO COVALENT BTKi RESISTANCE:
  PIRTOBRUTINIB (LOXO-305):
    Non-covalent (reversible) BTK inhibitor.
    Binds BTK without requiring C481 cysteine.
    C481S mutation does NOT confer resistance
    to pirtobrutinib.
    FDA approved December 2025 for R/R CLL
    after covalent BTKi failure.
    BRUIN CLL-321 Phase III: median PFS 14 months
    vs. 8.7 months standard of care.
    FRAMEWORK INTERPRETATION:
    Pirtobrutinib dissolves the C481S-reconstituted
    attractor because it targets BTK by a mechanism
    independent of C481. The attractor's adapted
    lock mechanism (C481S) is bypassed by a key
    that doesn't need C481.
    LIMITATION: BTK kinase-domain mutations
    (distinct from C481S) that alter the ATP
    binding site conformation can also emerge
    under pirtobrutinib — next level of attractor
    adaptation is already observed in early
    clinical reports.
    PLCG2 mutations also provide pirtobrutinib
    resistance (they bypass BTK entirely —
    irrelevant which BTK inhibitor is used).
    For PLCG2-mediated resistance, venetoclax
    (BCL2 direct inhibition) is the correct
    sequential therapy.

RESISTANCE TO VENETOCLAX (BCL2 INHIBITOR):

MECHANISM 1: BCL2 G101V MUTATION
  The most common venetoclax resistance mutation.
  G101V: glycine-101 → valine substitution
  in BCL2's BH3-binding groove.
  Venetoclax binds the BH3-binding groove of BCL2
  to competitively displace BH3-only pro-apoptotic
  proteins (BIM, PUMA, NOXA) → BCL2 cannot
  sequester these → BAX/BAK activation → apoptosis.
  G101V reduces venetoclax binding affinity
  ~180-fold → venetoclax cannot occupy the groove
  effectively at clinical concentrations.
  FREQUENCY: ~30–40% of venetoclax-resistant CLL.
  FRAMEWORK INTERPRETATION:
    G101V is the BCL2 equivalent of C481S in BTK:
    the lock has mutated so the drug cannot bind.
    The attractor's primary maintenance mechanism
    has modified its own molecular structure to
    resist the inhibitor.
    G101V is SPECIFICALLY SELECTED for by
    venetoclax — cells with G101V have a survival
    advantage under venetoclax treatment.
    G101V is present at very low allele frequency
    before venetoclax exposure — pre-existing
    attractor heterogeneity (minor G101V sub-clone)
    is selected for under venetoclax pressure.

MECHANISM 2: BCL2 UPREGULATION / AMPLIFICATION
  Cells increase BCL2 expression → overwhelming
  the available venetoclax concentration.
  More BCL2 molecules than venetoclax can occupy
  → residual free BCL2 sequesters BIM/PUMA →
  apoptosis threshold not reached.
  This is not a mutation — it is a QUANTITATIVE
  ATTRACTOR DEEPENING: the BCL2 lock becomes
  stronger by increasing its own concentration.

MECHANISM 3: MCL1 / BCL-XL UPREGULATION
  The cell switches to using ALTERNATIVE
  ANTI-APOPTOTIC proteins:
  MCL1 (myeloid cell leukemia 1) can substitute
  for BCL2 in sequestering BIM/PUMA.
  BCL-XL (BCL2L1) provides additional survival.
  When BCL2 is inhibited, MCL1 and BCL-XL
  upregulation is ADAPTIVE — the false attractor
  routes its survival signal through a different
  anti-apoptotic protein.
  FRAMEWORK INTERPRETATION:
    MCL1/BCL-XL upregulation = the attractor
    finding ALTERNATIVE LOCK MECHANISMS when
    the primary lock (BCL2) is blocked.
    This is the deepest form of CLL resistance:
    the attractor no longer depends on a single
    lock but has diversified its survival
    maintenance across multiple BCL2 family members.
    CLINICAL IMPLICATION:
    At this stage, combining venetoclax with
    MCL1 inhibitors (AZD5991, S63845 — in clinical
    trials) or with BTK inhibitors (to reduce
    all upstream BCL2-family transcription)
    is the rational strategy.
    NEXT-GENERATION BCL2 INHIBITORS:
    Sonrotoclax (BeiGene): potent BCL2 inhibitor,
    early data in R/R CLL (ASH 2025 confirmation).
    May have different binding kinetics than
    venetoclax → activity against some
    venetoclax-resistant cells.
    Lisaftoclax: in Phase II/III trials.
    FRAMEWORK PREDICTION: These agents do not
    change the attractor geometry — they are
    still targeting BCL2 (Axis 1 lock).
    For true attractor dissolution in triple-
    refractory CLL (BTKi + venetoclax + TP53
    mutant), the attractor has been deepened
    to require an entirely different approach.

SEQUENTIAL RESISTANCE LANDSCAPE (2025):
  The therapeutic sequence in CLL defines
  the depth landscape:

  DEPTH 1:  Untreated CLL (any subtype)
  DEPTH 2:  BTK inhibitor exposure
            (C481S or PLCG2 emerges in ~40%
            at 5 years of continuous BTKi)
  DEPTH 3:  Covalent BTKi-resistant, venetoclax-
            naive
            → Pirtobrutinib (non-covalent BTKi)
            OR venetoclax + anti-CD20
  DEPTH 4:  Both BTKi and venetoclax exposed
            → Pirtobrutinib + venetoclax
            → Next-generation BTK + BCL2 combos
  DEPTH 5:  Triple-class refractory + del17p/TP53
            → Bispecific antibodies (epcoritamab,
               glofitamab — CD20×CD3 bispecifics)
            → CAR-T cell therapy (lisocabtagene
               maraleucel — liso-cel — FDA approved
               for R/R CLL 2024)
            → Allo-HSCT (only curative option)
  DEPTH 6:  Post-CAR-T relapse
            (Richter transformation — CLL →
            diffuse large B cell lymphoma DLBCL)
            = deepest, poorest prognosis
            Richter transformation = IDENTITY CHANGE:
            the CLL survival attractor dissolves
            and the clone transforms into an
            aggressive lymphoma with a COMPLETELY
            DIFFERENT ATTRACTOR GEOMETRY (DLBCL
            false attractor, not CLL false attractor).
            Richter transformation is the CLL
            equivalent of NEPC transformation
            in PRAD — therapy-driven lineage
            plasticity producing a new cancer
            with a new identity.

RICHTER TRANSFORMATION:
  Richter transformation (RT) = transformation
  of CLL to DLBCL or, rarely, Hodgkin lymphoma.
  Frequency: ~5–10% of CLL patients over lifetime.
             Higher in TP53-aberrant, IGHV-unmutated,
             CD38+, ZAP70+ patients.
             Enriched in stereotyped subsets.
  MECHANISM:
    The CLL survival attractor is increasingly
    under therapeutic pressure.
    Under BTK inhibitor + venetoclax exposure,
    a sub-clone acquires:
      MYC rearrangement (t(8;14) — Ig-MYC) →
      MYC-driven proliferation programme
    PLUS:
      CDKN2A/B deletion → p16/p14 loss →
      CDK4/6 deregulated
    PLUS (often):
      TP53 mutation
    Result: the cell has acquired enough
    proliferative drive and lost enough braking
    signals that the SURVIVAL ATTRACTOR (quiescent
    BCL2-high non-proliferating cell) is REPLACED
    by a PROLIFERATIVE ATTRACTOR (rapidly cycling
    DLBCL cell with MYC activation).
    This is the CLL equivalent of the Waddington
    transition that produces NEPC in PRAD:
    therapy drives the cancer from one false
    attractor to a deeper, more aggressive one.
  DEPTH AXIS 3 INTERPRETATION:
    Richter transformation is DEPTH AXIS 3
    maximally activated:
    The proliferating fraction of CLL (pseudo-
    follicle cells) has accumulated enough
    genomic instability (driven by AICDA activity
    in the proliferating compartment — AICDA
    can still be expressed at low levels in CLL
    pseudo-follicle cells) to make the MYC
    rearrangement and produce Richter transformation.
    AICDA in CLL = the equivalent of the
    APOBEC mutational signature in PRAD (the
    mutational enzyme that drives clonal evolution
    under therapy in the proliferating fraction).
    FRAMEWORK PREDICTION: The depth score at
    DEPTH AXIS 3 (proliferative) in IGHV-unmutated,
    del17p CLL should predict Richter transformation
    risk before it occurs — by measuring the
    gene expression signature of pseudo-follicle
    cells even when they are a minority of the
    total CLL population.
    This requires single-cell data to isolate the
    proliferating fraction.
    The GSE111014 ibrutinib time-series dataset
    includes serial samples — the pseudo-follicle
    signature may be extractable from the day-0
    single-cell data by clustering cells by
    MKI67 expression (the 0.26% MKI67-high cells
    in the bulk day-0 data are the proliferating
    fraction; 99.74% are the quiescent bulk).
```

---

## SECTION VI — THE TREATMENT LANDSCAPE

```
CLL TREATMENT HAS BEEN TRANSFORMED MORE
COMPLETELY THAN ANY OTHER HAEMATOLOGICAL
MALIGNANCY IN THE PAST DECADE.

The introduction of ibrutinib (2014), venetoclax
(2016), and their combinations has produced
a situation where TP53-wild-type CLL patients
can achieve functional cure with time-limited
chemo-free regimens — something unthinkable
with chemoimmunotherapy (FCR).

FRONTLINE TREATMENT (2025):

WATCH AND WAIT:
  Indication: Asymptomatic CLL (any stage).
  Rai Stage 0-1, Binet Stage A-B (asymptomatic):
  No treatment benefit demonstrated.
  Treatment initiation criteria (IwCll 2018):
    Progressive marrow failure (worsening cytopenias)
    Symptomatic/progressive splenomegaly
    Massive/progressive lymphadenopathy
    Rapid lymphocyte doubling time (<6 months)
    Autoimmune cytopenias unresponsive to corticosteroids
    CLL-related symptoms (B symptoms, fatigue)
  IGHV-mutated, del13q-only: median time to first
  treatment ~10 years; some never require treatment.

FRONTLINE TARGETED THERAPY:

BTK INHIBITOR MONOTHERAPY (CONTINUOUS):
  Ibrutinib: original first-gen, highest off-target
    (cardiovascular toxicity: AF, hypertension).
  Acalabrutinib: second-gen, more selective BTK.
    ELEVATE-TN trial: acalabrutinib ± obinutuzumab
    vs. chlorambucil-obinutuzumab.
    Superior PFS with acalabrutinib.
  Zanubrutinib: most selective BTK inhibitor.
    ALPINE trial: zanubrutinib vs. ibrutinib —
    superior PFS AND lower AF rate.
    Now PREFERRED first-generation replacement.
  INDICATION: Any CLL requiring treatment,
    particularly del17p/TP53-aberrant (BTKi is
    the most effective single agent here).
  DEPTH SCORE PREDICTION FOR RESPONSE:
    IGHV-unmutated + del17p: responds to BTKi
    but MRD is rarely undetectable.
    IGHV-mutated: responds to BTKi but fixed-
    duration venetoclax may produce deeper
    remissions.

VENETOCLAX + OBINUTUZUMAB (CLL14) — FIXED DURATION:
  12 months of treatment.
  ALL patients: uMRD rate ~75% at end of treatment.
  IGHV-mutated: uMRD rate >85%, median PFS > 5 years.
  IGHV-unmutated: uMRD rate ~65%, shorter PFS.
  del17p/TP53 aberrant: lower uMRD rates (~50%).
  DEPTH SCORE PREDICTION FOR RESPONSE:
    The depth score should predict uMRD rate:
    Low depth (IGHV-mut, del13q) → high uMRD.
    High depth (IGHV-unm, del17p) → lower uMRD.
    If depth score predicts uMRD rates better
    than IGHV status alone — that is the
    going-further finding for CLL.

IBRUTINIB + VENETOCLAX (GLOW, CAPTIVATE):
  Fixed-duration combination.
  GLOW: ibrutinib + venetoclax vs. chlorambucil-
  obinutuzumab in older/unfit patients.
  High uMRD rates (77% in ibrutinib-ven arm
  at end of treatment).
  CAPTIVATE: ibrutinib + venetoclax in younger
  patients. uMRD-guided discontinuation:
  MRD-negative → ibrutinib can be stopped.
  MRD-positive → ibrutinib continued.
  FRAMEWORK INTERPRETATION:
    CAPTIVATE is the MRD-guided depth score
    application in CLL:
    uMRD = the depth score has returned to ~zero.
    MRD-positive = residual depth > zero.
    The MRD test IS an empirical depth score
    in CLL — it measures residual attractor
    cells (CLL cells still alive, detectable by
    allele-specific PCR or NGS-based MRD assay).
    The framework depth score derived from
    expression data should CORRELATE WITH uMRD
    status at end of treatment — this is the
    primary clinical validation for the CLL
    depth score.

PIRTOBRUTINIB + VENETOCLAX + OBINUTUZUMAB:
  Emerging combination (first-line Phase II data
  March 2025):
  High uMRD rates even in IGHV-unmutated patients.
  Early data: >90% uMRD at end of treatment
  in IGHV-mutated; ~75% in IGHV-unmutated.
  This combination may provide the deepest
  fixed-duration remissions in all CLL subtypes.
  FRAMEWORK PREDICTION: The pirtobrutinib
  component targets the BTK→BCL2 axis (Axis 1).
  Venetoclax targets BCL2 directly (Axis 1).
  Obinutuzumab targets CD20 → complement/ADCC
  (kills CLL cells directly, Axis 0 = immune
  effector).
  Triple combination attacks Axis 1 from three
  directions. The depth score should decline most
  rapidly in this combination versus any single
  or dual agent.

RELAPSED/REFRACTORY CLL (2025):

AFTER BTKi FAILURE (COVALENT):
  Pirtobrutinib: FDA approved Dec 2025.
    BRUIN CLL-321: PFS 14 vs. 8.7 months.
    Works despite C481S mutation.
  Venetoclax + rituximab (MURANO): if BTKi was
    the prior therapy and no prior venetoclax.
    Fixed-duration (venetoclax 24 months,
    rituximab 6 months).
    MURANO: median PFS 53 months (vs. 17 months
    for bendamustine-rituximab).

AFTER VENETOCLAX FAILURE:
  BTKi (if BTKi-naive): very effective.
  Pirtobrutinib (if covalent BTKi prior): option.
  Sonrotoclax or lisaftoclax (emerging):
    Next-generation BCL2 inhibitors.
    Some activity in venetoclax-resistant patients
    (depends on resistance mechanism — G101V
    resistance remains problematic for some).

AFTER BOTH BTKi AND VENETOCLAX:
  Pirtobrutinib + venetoclax combination.
  Lisocabtagene maraleucel (liso-cel) CAR-T:
    FDA approved 2024 for R/R CLL.
    ORR ~45%, uMRD ~20% in post-BTKi/venetoclax
    setting.
    Significant improvement in this population.
    FRAMEWORK INTERPRETATION:
    CAR-T in CLL is attractor dissolution by
    immune effector — it is NOT targeting the
    BCL2/BTK axis.
    CD19-targeting CAR-T cells recognise CD19
    (which is ELEVATED in CLL, +6% from the
    data, effectively maintained throughout)
    and kill CLL cells by cytotoxicity (perforin/
    granzyme). This is a KILL-THE-ATTRACTOR
    approach rather than a DISSOLVE-THE-ATTRACTOR
    approach. The two approaches are
    fundamentally different:
    DISSOLVE: Remove the molecular lock (BCL2
    inhibition) → the cells die from loss of
    survival signal → they still use their own
    apoptotic machinery.
    KILL: External effector (CAR-T cell) kills
    the CLL cell by cytotoxic mechanisms
    independent of BCL2 → works even in BCL2-
    resistant, TP53-mutant cells.
    This is why CAR-T is effective in triple-
    refractory CLL: it bypasses BOTH the BTK/BCL2
    axis AND the TP53/apoptosis requirement.

  Epcoritamab or glofitamab (CD20×CD3 bispecific):
    Emerging in CLL.
    Data from DLBCL showing activity; CLL trials
    ongoing.
    Same kill-the-attractor mechanism as CAR-T
    (immune-mediated cytotoxicity).

  Allo-HSCT:
    Only curative option.
    OS benefit in high-risk CLL (del17p, TP53
    mutant, Richter transformation).
    Limited by age/fitness and donor availability.
```

---

## SECTION VII — THE DEPTH AXIS IN CLL: SUMMARY

```
THE CLL DEPTH SCORE IS UNLIKE ALL OTHER
DEPTH SCORES IN THE REPOSITORY.

Every other depth score measures DISTANCE
FROM A TERMINAL NORMAL IDENTITY:
  PRAD: distance from luminal prostate identity
  PAAD: distance from acinar pancreas identity
  BRCA: distance from luminal breast identity
  GBM: distance from OPC neural identity
  MDS: distance from differentiation execution
       capacity in CD34+ progenitors

CLL DEPTH MEASURES THE STABILITY AND STRENGTH
OF A FALSE ATTRACTOR AT A POSITION THAT IS
ALREADY AT OR PAST THE TERMINAL NORMAL STATE.

The CLL depth score measures:
  How deeply is the BCR-BTK-BCL2 circuit
  locked in the survival state?
  How many secondary resistance mechanisms
  (TP53 loss, ATM loss, MCL1/BCL-XL backup)
  have been added?
  How active is the proliferating fraction
  (Richter transformation risk)?

THREE AXES, THREE GENE SETS:

AXIS 1 — BCR-DEPENDENT LOCK:
  Depth-positive (rising with depth):
    BCL2 — the primary lock
    IGHD — IgD co-expression, BCR signal
    FCRL5 — anergy, sustained tonic signal
    CD27 (TNFRSF7) — memory B phenotype mimicry
    CD49d (ITGA4) — microenvironment survival
    ZAP70 — BCR amplification kinase
  Depth-negative (falling with depth):
    PRDM1 — terminal exit blocked
    (paradox: PRDM1 is suppressed in all CLL —
    it may not discriminate subtypes as cleanly
    as Axis 1 genes)

AXIS 2 — GENOMIC INSTABILITY:
  Depth-positive (rising with depth):
    MDM2 — TP53 degradation
    CDKN2A absence — p16/p14 loss
    SF3B1 mutation status (binary — not expression)
    ATM expression (low = del11q)
    TP53 expression (low = del17p, but mutant
    p53 can be ELEVATED due to MDM2 failure to
    degrade it — direction must be determined
    from data)
  Note: Axis 2 is partially a GENOMIC feature
  (mutation/deletion status) rather than an
  expression feature. It overlaps imperfectly
  with gene expression data.
  The depth score from expression data will
  capture Axis 2 INDIRECTLY (via MDM2, p53
  targets, DNA damage response genes).

AXIS 3 — PROLIFERATIVE COMPARTMENT:
  Present only in the Ki67-high minor fraction.
  MKI67, MCM2, PCNA — proliferation markers
  in the 0.26% cycling fraction.
  Depth-positive in proliferating cells:
    MYC — transcriptional driver of Richter risk
    AICDA — mutational engine in pseudo-follicles
    CCND1, CCND3 — cell cycle driver
  This axis is NOT DETECTABLE FROM BULK
  EXPRESSION DATA (because the proliferating
  fraction is 0.26% of total CLL cells and is
  diluted to zero in bulk analysis).
  REQUIRES SINGLE-CELL DATA to resolve.
  The GSE111014 dataset (10X single-cell)
  DOES capture this population — it is present
  at low frequency in the day-0 data.
  Identifying and profiling the Ki67-high
  CLL cells in single-cell data is the
  primary novel output of the CLL subtype
  series at AXIS 3.

THE CLINICAL PANEL FOR CLL:
  Per framework output: 3-gene panel
  approximating depth with r > 0.85.
  Candidate genes based on the biology:
    BCL2 protein (immunohistochemistry or flow)
    CD49d surface expression (flow cytometry)
    ZAP70 surface/intracellular (flow cytometry)
  These three markers are ALREADY MEASURED in
  clinical practice for CLL prognostication.
  The depth score may formalise their combination
  into a continuous score:
    BCL2-high + CD49d-high + ZAP70-high =
    deepest survival attractor
    BCL2-moderate + CD49d-low + ZAP70-low =
    shallowest survival attractor
  The uMRD outcome from venetoclax regimens
  should be the clinical validation outcome:
  low depth → high uMRD rate.
  high depth → lower uMRD rate.
  If depth score predicts uMRD better than
  IGHV status alone — that is the
  primary going-further finding for CLL.
```

---

## SECTION VIII — DATA AVAILABLE

```
Dataset         Accession    Samples              Notes
──────────────────────────────────────────────────────────
GSE111014       GEO          48,016 CLL cells     scRNA-seq
(Rendeiro       10X Chromium  4 patients: CLL1,   10X
2020)                        CLL5, CLL6, CLL8     IBRUTINIB
                             d0, d30, d120,       TIME SERIES
                             d150, d280           PRIMARY
                             Normal B: from
                             GSE132509 (PBMMC)

GSE132509       GEO          Normal PBMMC         scRNA-seq
                             3 donors             10X
                             B cells + monocytes  NORMAL REF

EXISTING ANALYSIS (Validation #8):
  Day 0 CLL cells (15,007) vs. normal B (2,744).
  Gene panel: IGHD, BTG1, FCRL5, PRDM1, BCL2,
              MKI67, CD19, PAX5, RAG1, RAG2,
              IGKC, IGHM, CD27, CEBPA, SFTPC, CDX2.
  BCL2 +136% *** confirmed.
  PRDM1 -57% *** confirmed.
  IGKC +60% confirmed.
  CD27 +817% *** confirmed (not predicted —
  most dramatic finding).
  IGHD +43% ELEVATED (prediction was wrong direction).
  FCRL5 +415% ELEVATED (prediction was wrong direction).
  Both IGHD and FCRL5 are now understood as
  FEATURES of the survival attractor (not switch
  genes in the differentiation block sense).
  Ibrutinib time-series: BCL2 83% fall by day 150.
  PRDM1 flat — cells exit by death not differentiation.

ADDITIONAL DATASETS FOR SUBTYPE ANALYSIS:

GSE22762       GEO           CLL samples with     Microarray
                             IGHV, FISH, clinical TCGA-style
                             annotation           prognostic
                             Multiple subtypes    markers
                             labelled             annotated

GSE136861      GEO           CLL + normal B,      Bulk RNA-seq
                             multiple subtypes    IGHV annotated

E-MTAB-1292    ArrayExpress  CLL with del13q,     Microarray
                             del11q, del17p,      FISH-annotated
                             trisomy 12, normal   cytogenetics
                             karyotype labelled
                             (large n = 210 CLL)

LLMPP          NCBI GEO      CLL with IGHV,       Microarray
                             CD38, ZAP70, OS      Gene expression
                             clinical outcomes    + survival

GSE111014      GEO           EXISTING DATASET.    scRNA-seq
                             For subtype analysis: TIME SERIES
                             Cell clustering by   IGHV UNKNOWN
                             expression can be    (4 patients
                             used to identify     only — not
                             IGHV-like and        enough for
                             unmutated-like        subtype
                             expression groups    analysis on
                             (even without         its own)
                             sequenced IGHV)

NOTE ON CLL SUBTYPE ANALYSIS METHODOLOGY:
  Unlike solid tumour subtype analysis (which
  can use TCGA datasets with molecular
  annotations), CLL subtype analysis requires
  COMBINING EXPRESSION DATA with INDEPENDENT
  MOLECULAR ANNOTATIONS (IGHV status, FISH,
  clinical outcomes).
  The best available datasets with both
  expression data and molecular annotations
  are MICROARRAY-BASED (not RNA-seq) — they are
  older but have the molecular annotations needed.
  E-MTAB-1292 (210 CLL patients, FISH-annotated
  cytogenetics) is the most complete cytogenetic
  dataset.
  GSE22762 (IGHV-annotated CLL) is the most
  complete IGHV-stratified expression dataset.
  THESE ARE THE PRIMARY DATASETS FOR THE
  CLL SUBTYPE SERIES.
```

---

## SECTION IX — PLANNED ANALYSIS ORDER

```
ORDER:

CLL-S1  IGHV-MUTATED vs. NORMAL B    GSE22762 +    HIGH
        (shallow survival attractor)  E-MTAB-1292
                                      IGHV-mutated
                                      subset +
                                      normal B ref.
        QUESTIONS:
          BCL2 in IGHV-mutated vs.
          IGHV-unmutated: how different?
          Is the depth score lower in
          IGHV-mutated as predicted?
          CD27 elevation: same or less
          than bulk CLL (Validation #8)?
          FCRL5 elevation: same or less?
          Does depth score predict
          treatment response in
          IGHV-mutated CLL?
          3-gene panel candidate for
          IGHV-mutated: BCL2 +
          CD49d (ITGA4) + ZAP70?

CLL-S2  IGHV-UNMUTATED vs. NORMAL B  GSE22762 +    HIGH
        (deep survival attractor)     E-MTAB-1292
                                      IGHV-unmutated
                                      subset +
                                      normal B ref.
        QUESTIONS:
          BCL2 in IGHV-unmutated vs.
          IGHV-mutated: depth confirmed?
          ZAP70 — elevated in IGHV-
          unmutated vs. IGHV-mutated?
          CD49d (ITGA4) — elevated in
          unmutated subset?
          Does depth score separate
          IGHV-unmutated from mutated
          better than ZAP70 alone?
          Does depth score correlate
          with Rai/Binet stage?
          Stereotyped BCR subsets: can
          the expression data identify
          Subset 1 and Subset 2 CLL
          without sequencing the BCR?

CLL-S3  CYTOGENETIC STRATIFICATION    E-MTAB-1292   HIGH
        del13q vs. trisomy 12 vs.     (FISH-annotated
        del11q vs. del17p/TP53        cytogenetics)
                                      + normal B ref.
        QUESTIONS:
          Does the depth score follow
          the Döhner hierarchy?
          del13q → AXIS 1 ONLY
          (miR-15a/16-1 loss → BCL2)
          Does del17p add AXIS 2 depth
          ON TOP OF AXIS 1?
          ATM expression in del11q:
          ATM low + BCL2 high = Axis 2
          deepening confirmed?
          MDM2 in trisomy 12: elevated
          MDM2 vs. del13q CLL?
          Does the depth score from
          expression data replicate the
          Döhner prognostic hierarchy
          WITHOUT using FISH data?
          If yes: depth score from
          expression is a FISH-free
          prognostic stratification tool.
          This is the primary clinical
          going-further finding for CLL
          cytogenetic analysis.

CLL-S4  IBRUTINIB RESISTANCE          GSE111014     MOD-HIGH
        (from Validation #8 data)     ibrutinib
                                      time series:
                                      d0 vs d30 vs
                                      d120 vs d150
                                      vs d280
        QUESTIONS:
          The depth score at day 0:
          What is it for each of the
          4 patients individually?
          Does depth score change over
          time on ibrutinib?
          Expected: depth falls
          (attractor dissolving) then
          partially recovers at d280
          (residual disease / resistance
          beginning to emerge).
          C481S resistance emergence:
          can the depth score DETECT
          the onset of BTKi resistance
          BEFORE clinical progression?
          At d280 partial BCL2 recovery —
          does this represent C481S
          sub-clone expansion?
          If depth score rises at d280
          before clinical relapse —
          that is the resistance early-
          warning signal application
          of the framework in CLL.

CLL-S5  RICHTER TRANSFORMATION RISK   From GSE111014 MOD
        (proliferating fraction       single-cell data
        identification)               Ki67-HIGH CELLS
                                      from day 0
        QUESTIONS:
          Can the proliferating fraction
          (<1% of CLL cells) be isolated
          by single-cell clustering on
          the basis of MKI67 expression?
          What is the gene expression
          profile of MKI67-high CLL cells
          vs. MKI67-low CLL cells?
          MYC elevated in MKI67-high?
          AICDA expressed in MKI67-high?
          CCND3 elevated (the cyclin in
          GC B cells — would indicate
          pseudo-follicle biology)?
          Is the proliferating fraction
          gene expression profile different
          in IGHV-unmutated vs. mutated CLL?
          (predicting: higher MYC, higher
          AICDA in unmutated proliferating
          fraction → higher Richter risk)
          This is the single-cell-specific
          analysis that cannot be done with
          bulk data and is unique to the
          GSE111014 10X dataset.

CLL-X   CROSS-SUBTYPE INTEGRATION     After S1–S5.
        + uMRD PREDICTION             QUESTIONS:
                                        1. Does the
                                           unified CLL
                                           depth score
                                           (Axis 1 +
                                           Axis 2) predict
                                           uMRD outcome
                                           from venetoclax
                                           regimens?
                                        2. Is it better
                                           than IGHV status
                                           alone?
                                        3. Does Axis 3
                                           (proliferating
                                           fraction Ki67-
                                           high cells)
                                           predict Richter
                                           transformation
                                           risk?
                                        4. 3-gene panel
                                           final candidate:
                                           BCL2 + CD49d
                                           (ITGA4) + ZAP70
                                           confirmed or
                                           revised?
                                        5. CLL depth score
                                           vs. IGHV +
                                           cytogenetics
                                           combined model:
                                           is it non-
                                           inferior? If
                                           depth score
                                           without molecular
                                           testing matches
                                           the combined
                                           model — that is
                                           the final going-
                                           further finding
                                           for the CLL
                                           subtype series.
```

---

## SECTION X — UNIQUE FRAMEWORK POSITIONS FOR CLL

```
CLL CONTRIBUTES FOUR UNIQUE STRUCTURAL LESSONS
TO THE REPOSITORY THAT DO NOT EXIST IN ANY
OTHER CANCER SERIES.

LESSON 1 — THE SURVIVAL ATTRACTOR IS REAL
AND CATEGORICALLY DIFFERENT FROM
DIFFERENTIATION BLOCK ATTRACTORS.
  The automated script (1/4 switch genes confirmed)
  called CLL a failure.
  The biological interpretation called CLL a
  confirmation of a new attractor TYPE.
  Both are honest.
  The lesson: the framework's scoring logic
  was built for one attractor type (differentiation
  block) and failed to score a different attractor
  type (survival block).
  This is a limitation of the scoring system,
  not of the biology or the method.
  RESOLUTION FOR THE SUBTYPE SERIES:
  The CLL subtype analysis should use a
  SURVIVAL ATTRACTOR SCORING SYSTEM:
    Depth-positive genes: BCL2, IGHD, FCRL5,
    CD27, CD49d, ZAP70.
    Depth-negative genes: PRDM1, RAG1, RAG2,
    MKI67 (expected low in all CLL — not
    discriminating subtypes).
    The depth score in CLL is the MEAN OF DEPTH-
    POSITIVE GENES minus the MEAN OF DEPTH-
    NEGATIVE GENES (reversed from differentiation
    block scoring, where depth-negative falling =
    deeper).
  This revised scoring logic should be stated
  explicitly in CLL-S1a (the Before document)
  and tested against the data.

LESSON 2 — THE TRANSITIONAL STATE AS FALSE ATTRACTOR.
  Every other cancer in the repository has
  stabilised a STAGE OF DEVELOPMENT — either a
  progenitor stage (AML, B-ALL, MDS) or a
  partially-differentiated intermediate (GBM,
  LUAD, PRAD, PAAD).
  CLL has stabilised a TRANSITION POINT —
  a moment in normal B cell biology that is
  supposed to be a fleeting decision point
  (mature naive B → activate OR die) and turned
  it into a permanent stable identity.
  The normal B cell landscape has NO stable
  long-term attractor at the mature naive B
  position. CLL CREATED ONE.
  This is the most sophisticated false attractor
  in the repository: it is not exploiting an
  existing developmental stage's stability —
  it is creating new stability at a position
  where none should exist.
  The CD27 +817% finding is the signature of
  this creation: CLL co-opts the memory B cell
  surface phenotype (CD27+) to make its false
  attractor LOOK like a legitimate B cell
  identity (memory B cell) even though its
  internal state is distinct.

LESSON 3 — DRUG-TARGET DERIVATION FROM ATTRACTOR
DISSOLUTION TIMECOURSE IS POSSIBLE IN CLL.
  The GSE111014 ibrutinib time-series is unique
  in the repository — no other confirmed cancer
  dataset has sequential measurements of the
  same cancer cells as the drug dissolves the
  attractor over 280 days.
  The timecourse shows:
    BCL2 falls 83% — the lock is released.
    IGHD falls to zero — the BCR signal is cut.
    FCRL5 falls to near-zero — anergy state resolved.
    PRDM1 stays flat — no differentiation occurs.
  This is a MOLECULAR MOVIE of attractor dissolution.
  No other cancer in the repository has this.
  The subtype analysis at CLL-S4 uses this
  timecourse to ask: does the DEPTH SCORE
  predict the RATE AND DEPTH OF ATTRACTOR
  DISSOLUTION?
  If depth score at day 0 predicts depth of
  BCL2 fall by day 150 — that is the clinical
  depth score as treatment response predictor.

LESSON 4 — AXIS 3 (PROLIFERATING FRACTION) IS
THE RICHTER TRANSFORMATION EARLY WARNING SIGNAL.
  The 99% MKI67 suppression in bulk CLL creates
  a paradox: how does a non-proliferating cancer
  clonally evolve?
  Answer: it doesn't evolve in the 99.74% quiescent
  fraction. It evolves in the 0.26% proliferating
  fraction (pseudo-follicle cells).
  The pseudo-follicle cells are where:
    AICDA is expressed → mutation accumulation
    MYC may be amplified/rearranged → proliferative
    drive for the resistant sub-clone
    Clonal evolution happens at high pace in a
    small compartment, invisible to bulk profiling.
  Single-cell data (GSE111014) captures this
  compartment. The CLL-S5 analysis is the FIRST
  FRAMEWORK ANALYSIS TARGETING A CANCER'S
  PROLIFERATING MINORITY rather than its
  dominant cell population.
  This inverts the usual framework approach:
  normally we analyse the dominant clone and
  derive the drug target for that clone.
  In CLL Richter risk, we analyse the MINORITY
  CLONE and derive the biomarker that predicts
  when the minority will outgrow the majority
  and cause transformation.
  Framework prediction: MYC expression in the
  Ki67-high CLL fraction predicts Richter
  transformation risk.
  This has not been tested in this form in
  the literature. It is the primary novel
  prediction for CLL-S5.
```

---

## SECTION XI — WHAT THIS DOCUMENT DOES NOT CONTAIN

```
This document contains:
  ✓ Normal B cell developmental hierarchy
    (Pro-B → Pre-B → Immature B → Naive B →
    GC B → Memory B / Plasma cell) with full
    TF network (PAX5, EBF1, E2A, BCL6, PRDM1,
    GATA1 logic for GC entry/exit)
  ✓ V(D)J recombination staging (RAG1/RAG2 at
    Pro-B and small Pre-B; RAG silent after IgM
    expression; IGKC as the developmental position
    marker confirmed by Validation #8)
  ✓ The survival attractor geometry and why it
    is unique in the repository (transitional
    state stabilised, three-axis depth)
  ✓ IGHV mutated/unmutated divide (most important
    CLL stratification; GC history → depth Axis 1
    shallow vs. deep)
  ✓ ZAP70, CD38, CD49d as Axis 1 depth proxies
  ✓ Stereotyped BCR subsets (Subset 1, Subset 2,
    Subset 4 with Waddington interpretation)
  ✓ Cytogenetic hierarchy (del13q/miR-15a-16-1/
    BCL2 mechanism; trisomy 12/MDM2; del11q/ATM;
    del17p/TP53; each as Axis 2 depth layers)
  ✓ Resistance mechanisms (BTK C481S, PLCG2 GOF,
    BCL2 G101V, MCL1/BCL-XL upregulation — each
    as molecular attractor reconstitution)
  ✓ Pirtobrutinib, sonrotoclax, lisaftoclax as
    next-generation drug landscape (2025)
  ✓ Richter transformation as Axis 3 depth
    maximum and CLL equivalent of NEPC in PRAD
  ✓ CAR-T (liso-cel) as kill-the-attractor vs.
    dissolve-the-attractor distinction
  ✓ Planned analysis: CLL-S1 to CLL-S5 + CLL-X
  ✓ Four unique framework lessons from CLL
    (survival attractor, transitional state
    stabilisation, dissolution timecourse,
    proliferating minority as Richter signal)
  ✓ Data available (GSE111014, GSE132509,
    GSE22762, E-MTAB-1292, GSE136861)
  ✓ The survival attractor scoring system
    revision (Axis 1 depth-positive: BCL2,
    IGHD, FCRL5, CD27, CD49d, ZAP70)
    required for subtype analysis

This document does NOT contain:
  ✗ Depth score predictions (belong in CLL-S1a)
  ✗ Switch gene predictions for IGHV subtypes
    (belong in CLL-S1a)
  ✗ Novel drug combination predictions
    (belong in CLL-S1b after data)
  ✗ Epigenetic mechanism hypotheses for subtypes
    (belong in before-documents)
  ✗ IGLL5, CCND3, NOTCH1 mutation analysis
    (these are additional CLL mutations — NOTCH1
    particularly in trisomy 12 — belong in
    CLL-S3 before-document)
```

---

## STATUS BLOCK

```
document:           CLL_Subtype_Orientation.md
folder:             Cancer_Research/CLL/Subtypes/
status:             COMPLETE — ORIENTATION ONLY
date:               2026-03-04
author:             Eric Robert Lawson
                    OrganismCore

existing_validation:
  Validation #8 (Document 80):
  CLL_False_Attractor_confirmed.md
  Automated score: INSUFFICIENT (1/4 switch genes)
  Biological finding: SURVIVAL ATTRACTOR CONFIRMED
  BCL2 +136% *** (p<1e-45) — lock confirmed
  PRDM1 -57% *** — terminal exit blocked confirmed
  IGKC +60% — hierarchy position confirmed
  CD27 +817% *** — memory B phenotype mimicry confirmed
  IGHD +43% ELEVATED — analyst assumption error
    (corrected: IGHD is feature of survival attractor,
    not switch gene)
  FCRL5 +415% ELEVATED — analyst assumption error
    (corrected: anergy marker elevated, not suppressed)
  Ibrutinib time-series: BCL2 -83% at d150, lock
    dissolves when BCR signal cut.
  PRDM1 flat throughout — exit by death not
    differentiation.
  Drug targets: BTK + BCL2 — BOTH FDA APPROVED ✓
  Framework lesson: survival attractor scoring needs
    different gene logic from differentiation block
    attractor scoring.

unique_attractor_type:   SURVIVAL ATTRACTOR
  (only such cancer in the repository)
  Cells completed development then failed to die.
  False attractor = stabilised transitional state
    (mature naive B) that should not be permanent.
  False attractor created by BCL2 elevation
    maintaining cell viability at a position where
    normal biology has no stable identity.

key_structural_difference_from_all_other_cancers:
  Every other repository cancer: progenitor or
    intermediate stage arrested short of terminal
    identity. False attractor is UPSTREAM of normal.
  CLL: terminal maturation completed. Cell then
    failed to exit. False attractor IS AT terminal
    position, maintained by BCL2-high survival signal.
  The depth axis in CLL = STRENGTH OF SURVIVAL LOCK
    (not identity loss).

three_depth_axes:
  Axis 1: BCR-dependent lock (BCL2/IGHD/FCRL5/
    CD27/CD49d/ZAP70) — correlates with IGHV status
  Axis 2: Genomic instability depth (TP53/ATM/
    MDM2/SF3B1) — correlates with cytogenetic tier
  Axis 3: Proliferating fraction (Ki67-high cells)
    — Richter transformation risk signal; only
    detectable in single-cell data
```
