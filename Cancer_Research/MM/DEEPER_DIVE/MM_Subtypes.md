# MULTIPLE MYELOMA — SUBTYPE ORIENTATION DOCUMENT
## Before Any Subtype Analysis Begins
## OrganismCore | 2026-03-04 | Author: Eric Robert Lawson

---

## PURPOSE OF THIS DOCUMENT

```
This document exists before any script runs.
It contains no depth score predictions.
It contains no specific switch gene claims.
It contains no epigenetic hypotheses.

What it contains:

  A complete map of the multiple myeloma (MM)
  molecular subtype landscape — the IgH
  translocation classes, the hyperdiploid class,
  the high-risk cytogenetic lesions (del(17p),
  amp(1q)), the MAF and MMSET subtypes, the
  normal plasma cell differentiation circuit
  that the framework must measure, the MGUS-
  to-SMM-to-MM-to-RRMM progression continuum,
  the unique features of the MM false attractor
  that differentiate it from every other cancer
  in the repository, and the full treatment
  landscape from frontline quadruplet therapy
  through BCMA-targeting bispecifics and CAR-T.

MM HOLDS A UNIQUE STRUCTURAL POSITION IN
THE REPOSITORY.

  Unlike all other cancers in the repository:

  THE MYELOMA CELL IS A TERMINALLY DIFFERENTIATED
  CELL THAT HAS REFUSED TO DIE.

  In every other solid tumour:
    A progenitor or early cell type acquires
    mutations and FAILS TO DIFFERENTIATE.
    It is stuck at an immature state.
    It expresses progenitor markers and loses
    lineage identity.
    The depth axis measures how far it has
    traveled FROM its normal terminal identity.

  In multiple myeloma:
    The myeloma cell is a PLASMA CELL — the most
    terminally differentiated B cell in the body.
    It is ALREADY differentiated.
    It RETAINS its plasma cell identity:
    IRF4 high, BLIMP1 high, XBP1s active,
    BCMA expressed, CD138 expressed, high
    immunoglobulin (M-protein) production.
    What is wrong is NOT that it has lost its
    differentiation programme.
    What is wrong is that it has UNCOUPLED
    PROLIFERATION and SURVIVAL from the normal
    plasma cell terminal fate — which should be
    growth arrest and programmed death.

  THE MYELOMA FALSE ATTRACTOR:
    Normal plasma cells live 3–5 days, then die.
    They are post-mitotic and apoptosis-prone
    (BCL2-family balanced toward apoptosis;
    the plasma cell's job is to make antibody
    then die).
    Myeloma cells retain all plasma cell markers
    AND gain proliferative immortality —
    they have found an additional attractor
    that OVERLAPS with the plasma cell identity
    but ADDS self-renewal and apoptosis resistance.
    This is an AUGMENTED ATTRACTOR, not a
    replacement attractor.
    The myeloma cell is not a different cell type.
    It is a plasma cell that has acquired an
    additional survival programme on top of
    its existing identity.

  THE DEPTH AXIS IN MM IS THEREFORE DIFFERENT
  FROM ALL OTHER CANCERS:

    DEPTH IN PRAD, PAAD, STAD, BRCA, GBM:
    Distance FROM normal terminal identity.
    Depth = identity loss.
    Deeper = less like the normal cell type.

    DEPTH IN MDS:
    Differentiation execution capacity.
    Depth = inability to complete programmes.
    Deeper = less programme execution.

    DEPTH IN MM:
    Distance FROM normal plasma cell biology
    TOWARD an augmented survival/proliferative
    state that COEXISTS with plasma cell markers.
    Depth = added oncogenic programme ON TOP
    of retained plasma cell identity.
    Deeper = more MYC, more BCL2/MCL1,
    more proliferative signalling, more
    chromosomal instability — while RETAINING
    IRF4, BLIMP1, XBP1, BCMA, CD138, M-protein.

  THE CLINICAL CONSEQUENCE:
    All myeloma cells express plasma cell markers
    (BCMA, CD38, SLAMF7/CS1) regardless of depth.
    This is why CD38 antibodies (daratumumab),
    BCMA-targeting agents (teclistamab, CAR-T),
    and SLAMF7 antibodies (elotuzumab) work
    across ALL myeloma subtypes — they are
    targeting the retained plasma cell identity,
    not the depth-specific oncogenic programme.
    The depth-specific targets are the oncogenic
    add-ons: BCL2 (t(11;14)), NSD2/FGFR3
    (t(4;14)), MAF/MAFB (t(14;16)/(14;20)),
    TP53 loss (del(17p)), MYC activation
    (late event), CRBN vulnerability
    (lenalidomide/IMiD sensitivity).
```

---

## DOCUMENT METADATA

```
document_id:        MM_Subtype_Orientation
series:             MM (Multiple Myeloma — Subtypes)
folder:             Cancer_Research/MM/Subtypes/
date:               2026-03-04
author:             Eric Robert Lawson
status:             ORIENTATION ONLY — no predictions,
                    no analysis, no epigenetic hypotheses
next_document:      MM-S1a
                    t(11;14)/BCL2 subtype before-document
                    (the most tractable subtype —
                    venetoclax confirms BCL2 as convergence
                    node; the framework's CLL lesson applied
                    to plasma cells)
protocol_version:   Workflow_Protocol.md v2.0
existing_analysis:  Cancer_Research/MM/ — complete
                    bulk MM analysis from OrganismCore
                    series. Switch gene and drug target
                    identified. This subtype series resolves
                    which molecular subtype drives the bulk
                    signal and whether different subtypes
                    have different switch genes.
```

---

## SECTION I — NORMAL PLASMA CELL BIOLOGY

```
Understanding the normal plasma cell is the
prerequisite for understanding the myeloma
false attractor. The myeloma cell KEEPS almost
everything in this section. What it adds on top
is the oncogenic programme.

═══════════════════════════════════════════════════════
THE PLASMA CELL DIFFERENTIATION CIRCUIT
═══════════════════════════════════════════════════════

STAGE 1: NAIVE B CELL
  Location: Blood, lymph nodes, spleen, bone marrow.
  Function: Antigen surveillance.
            Expresses B cell receptor (BCR).
            Patrols for cognate antigen.
  Markers: CD19+, CD20+, PAX5+, BCL6−.
  Transcription factor identity:
    PAX5 — the B cell master TF.
    Activates all B cell genes.
    Represses non-B cell genes.
    PAX5 MUST BE SILENCED for plasma cell fate.
    PAX5 is still expressed in naive B cells.
    It is the barrier to plasma cell commitment.

STAGE 2: GERMINAL CENTRE B CELL
  Location: Germinal centre (GC) of secondary
            lymphoid organs.
  Function: Somatic hypermutation (SHM) of
            immunoglobulin variable regions for
            affinity maturation.
            Class-switch recombination (CSR) —
            IgM → IgG, IgA, or IgE switching.
            This is the stage at which IgH
            translocations occur in myeloma.
            The AID (activation-induced cytidine
            deaminase) enzyme that drives SHM
            and CSR creates DNA double-strand
            breaks — these breaks are the
            substrate for aberrant translocation
            joining the IgH locus to oncogenes
            (CCND1, FGFR3/NSD2, MAF, MAFB).
  Markers: CD19+, CD20+, BCL6+, AID high.
  TF identity:
    BCL6 — GC master TF.
    Represses MYC (prevents premature exit).
    Represses BLIMP1 (prevents premature
    plasma cell commitment).
    BCL6 must be turned OFF for plasma cell fate.
    In MM: some subtypes have BCL6 dysregulation.
    MAF (t(14;16)) and MAFB (t(14;20)) activate
    BCL6 targets and are expressed in cells
    that have partially maintained GC identity
    features — contributing to the MAF-subtype
    gene signature.

STAGE 3: THE COMMITMENT SWITCH
  (GC B cell → Plasmablast)
  The most important molecular event in normal
  plasma cell development.
  WHAT HAPPENS:
    BCL6 is downregulated.
    BLIMP1 (PRDM1) is upregulated by IRF4.
    BLIMP1 REPRESSES:
      PAX5 — B cell identity is erased.
      BCL6 — GC programme is erased.
      MYC — proliferation is halted.
      CIITA — MHC class II is downregulated.
    BLIMP1 ACTIVATES:
      IRF4 (positive feedback).
      XBP1 (via a secondary activation loop).
      J-chain (IgM polymerisation signal).
    IRF4 ACTIVATES:
      BLIMP1 (positive feedback loop with BLIMP1).
      MYC (in the MM context — not in normal PC).
      IRF4/MYC co-activation is the pathological
      feed-forward loop in myeloma.
    XBP1 SPLICING:
      ER stress during Ig production activates
      IRE1α, which splices XBP1 mRNA to produce
      XBP1s (the active transcription factor).
      XBP1s drives the UNFOLDED PROTEIN RESPONSE:
        ER expansion (more rough ER for Ig folding)
        Chaperone upregulation (BIP/GRP78, PDIA3)
        ERAD pathway (Ig folding quality control)
        Secretory pathway expansion
      XBP1s is ESSENTIAL for efficient Ig secretion.
      Without XBP1s: Ig misfolding and ER stress
      → apoptosis.
      IN MYELOMA: XBP1s activity is very high
      (M-protein production demands enormous ER
      capacity). This is the mechanistic basis for
      proteasome inhibitor sensitivity — bortezomib,
      carfilzomib, and ixazomib block proteasomal
      clearance of misfolded proteins → ER stress
      exceeds a threshold → apoptosis.
      ALL myeloma subtypes are proteasome inhibitor
      sensitive because ALL retain XBP1s-driven
      high-rate Ig production.
      This is the retained plasma cell identity
      being exploited therapeutically.

STAGE 4: PLASMABLAST
  Location: Blood (circulating), lymph nodes.
  Function: Rapid Ig production (IgM first,
            then class-switched Ig).
            Short-lived proliferating cell.
  Markers: CD138−, CD19+, CD38+, CD27+.
           Begins to express BCMA (TNFRSF17).
           BCMA is a survival receptor whose
           ligands are BAFF and APRIL (bone
           marrow stromal-derived cytokines).
           BCMA is the primary target of
           teclistamab, belantamab, ciltacabtagene
           autoleucel (cilta-cel/Carvykti), and
           idecabtagene vicleucel (ide-cel/Abecma).
           ALL MYELOMA CELLS RETAIN BCMA EXPRESSION
           because they are derived from plasmablasts
           and long-lived plasma cells — both of which
           require BCMA for niche survival signalling.

STAGE 5: LONG-LIVED PLASMA CELL (LLPC)
  Location: Bone marrow survival niches.
            Specialised microenvironmental zones
            where VCAM1+ stromal cells, osteoblasts,
            and APRIL/BAFF-producing cells provide
            survival signals.
  Function: Long-term antibody production.
            Maintain immunological memory
            without requiring repeated antigen
            stimulation.
            Live for YEARS (some estimated to
            last decades in humans).
  Markers: CD138+, CD38+, CD19−, CD20−.
           BCMA+, CD27+, CXCR4+.
           MHC class II LOW (CIITA repressed by BLIMP1).
           CD20 NEGATIVE (BLIMP1 represses CD20).
  THE SURVIVAL BALANCE IN NORMAL LLPC:
    BCL2 is expressed in LLPC (survival signal
    mediated by BAFF/APRIL via BCMA).
    MCL1 is also expressed.
    BIM (pro-apoptotic) is held in check by BCL2.
    NOXA and PUMA are expressed at low levels.
    The normal LLPC is SURVIVAL-BALANCED —
    it is not maximally anti-apoptotic (it will
    die in 3–5 days if removed from the niche).
    THE MYELOMA CELL SHIFTS THIS BALANCE:
    BCL2 is FURTHER elevated in t(11;14) myeloma
    (CCND1 co-upregulates BCL2 via the IgH enhancer
    effect at the t(11;14) locus).
    MCL1 is FURTHER elevated in non-t(11;14)
    myeloma (IL-6 → STAT3 → MCL1 transcription;
    amp(1q) → MCL1 copy number gain).
    The myeloma cell is MORE ANTI-APOPTOTIC than
    a normal LLPC — the balance has been shifted
    toward survival by the oncogenic add-on.

THE FULL PLASMA CELL CIRCUIT (SUMMARY):
  Naive B cell
    → [antigen + GC entry + BCL6 high] →
  GC B cell
    → [AID + SHM + CSR + IgH translocations here] →
  Post-GC B cell
    → [BCL6 off + BLIMP1 on + PAX5 off + XBP1s on] →
  Plasmablast / Long-lived plasma cell
    [IRF4 high, BLIMP1 high, XBP1s high,
    BCL2/MCL1 balanced, MYC OFF,
    BCMA high, CD138 high, CD38 high, CD20 off]

MYELOMA = The plasma cell state PLUS:
  MYC ON (the key BLIMP1 repression broken)
  BCL2 or MCL1 elevated (survival extended)
  Proliferative signalling (cyclin D1/2/3 elevated)
  Chromosomal instability (in high-risk subtypes)
  Bone marrow niche dependence (CXCR4, VLA-4 adhesion)
  IL-6 autocrine/paracrine signalling
  RANKL-mediated osteoclast activation (bone disease)
```

---

## SECTION II — THE MGUS-TO-MM-TO-RRMM PROGRESSION CONTINUUM

```
MULTIPLE MYELOMA HAS THE MOST WELL-CHARACTERISED
PRE-MALIGNANT SERIES IN ALL HAEMATOLOGICAL
MALIGNANCIES — COMPARABLE TO BARRETT'S OESOPHAGUS
IN SOLID TUMOURS.

MONOCLONAL GAMMOPATHY OF UNDETERMINED SIGNIFICANCE
(MGUS):
  Definition: M-protein in blood + <10% clonal
              plasma cells in bone marrow + NO
              end-organ damage (no CRAB:
              hyperCalcaemia, Renal failure,
              Anaemia, Bone lesions).
  Frequency:  ~3% of adults >50 years.
              ~5% of adults >70 years.
              >10% of adults >80 years (some series).
              Most common undetected pre-malignancy
              in the developed world.
  Progression: ~1% per year to MM or related
               malignancy (Waldenström's, AL amyloid).
  Biology:    The MGUS clone already has the
              PRIMARY IgH TRANSLOCATION (if a
              translocation subtype) OR trisomies
              (if a hyperdiploid subtype).
              The IgH translocations are FOUNDING
              EVENTS that occur in the GC B cell
              and establish the MGUS clone.
              MGUS has the first false attractor:
              a plasma cell clone that is stable,
              anti-apoptotic (BCL2/MCL1 balanced
              to survival) and slowly expanding —
              but NOT YET PROLIFERATING (MYC is
              still off, or near-off).
  Waddington: The MGUS clone sits in a very
              SHALLOW false attractor — close to
              the normal LLPC attractor, with
              only the IgH translocation effect
              distinguishing it from normal.
              The depth score in MGUS should be
              NEAR ZERO — barely distinguishable
              from normal plasma cells.

SMOLDERING MULTIPLE MYELOMA (SMM):
  Definition: M-protein ≥30 g/L OR ≥500 mg/24h
              urine M-protein OR 10–60% bone
              marrow plasma cells, WITH NO
              end-organ damage AND no biomarkers
              of malignancy.
  Progression: ~10%/year for first 5 years,
               then declining.
               High-risk SMM (based on serum
               free light chain ratio, M-protein
               level, and BM plasma cell %) can
               approach 50% progression at 2 years.
  Biology:    SMM represents the acquisition of
              SECONDARY EVENTS on top of the
              MGUS clone:
                MYC activation (enhancer hijacking
                or structural translocation — a
                LATE EVENT).
                Additional chromosomal copy number
                variations (gain of 1q, del 13,
                del 17p early events can begin here).
                Increasing bone marrow dependence
                (IL-6 signalling, CXCR4 upregulation).
              MYC is the CRITICAL SECONDARY EVENT
              that transitions MGUS to SMM/MM:
                BLIMP1 normally represses MYC in
                mature plasma cells.
                In MGUS/MM: a MYC translocation or
                enhancer hijacking (IGH, IGL, IGK
                enhancers are brought near MYC or
                vice versa) overcomes BLIMP1 repression.
                ~40–45% of MM have MYC translocation.
                ~60% of MM have native MYC enhancer
                gain-of-function (no structural variant —
                IRF4 and other TFs drive a plasma
                cell-specific MYC enhancer in an
                abnormally active state).
                MYC activation = the TRANSITION
                FROM STATIC CLONE (MGUS) TO
                PROLIFERATIVE CLONE (SMM/MM).
  Waddington: The SMM clone has moved DEEPER
              in the false attractor space:
              MYC is active → proliferative programme
              added to the plasma cell identity.
              The cell is now BOTH a plasma cell
              AND a cycling cell — two normally
              mutually exclusive states.
              BLIMP1 and MYC are simultaneously
              expressed (in normal plasma cells,
              BLIMP1 represses MYC; in myeloma,
              MYC escapes BLIMP1 repression).

SYMPTOMATIC MM (CRAB criteria):
  Definition: Clonal plasma cell disease + end-organ
              damage (CRAB) OR ≥60% BM plasma cells
              OR serum free light chain ratio >100
              OR >1 focal lesion on MRI.
  Biology:    On top of MGUS + SMM events:
                Osteolytic bone disease:
                Myeloma cells produce DKK1, RANKL,
                IL-3 → osteoclast activation and
                osteoblast suppression.
                Osteoclasts resorb bone → release
                bone matrix growth factors (IGF-1,
                TGF-β, IL-6 from matrix) → feed
                myeloma growth (the vicious cycle).
                RANKL from myeloma cells → osteoclast
                differentiation → bone resorption.
                DKK1 from myeloma cells → inhibits
                WNT signalling in osteoblasts →
                osteoblast suppression → no bone repair.
                Denosumab (anti-RANKL) and
                bisphosphonates (zolendronic acid)
                interrupt this cycle.
                Anaemia: marrow replacement + increased
                hepcidin (acute phase response) +
                RANKL-driven bone marrow fibrosis.
                Renal failure: light chain cast
                nephropathy (Bence Jones protein) or
                direct tubular toxicity.
                Hypercalcaemia: osteoclast-mediated
                calcium release.

RELAPSED/REFRACTORY MM (RRMM):
  Biology:    The critical distinction from first
              relapse to later-line disease:
              FIRST RELAPSE: often a single
              therapy-driven clonal evolution
              event (e.g., lenalidomide resistance
              via CRBN mutation or downregulation;
              bortezomib resistance via proteasome
              subunit mutation).
              LATER RELAPSE (≥3 prior lines): the
              dominant biology is CLONAL HETEROGENEITY
              and MICROENVIRONMENTAL ADAPTATION.
              Multiple subclones co-exist, some
              with del(17p), some with amp(1q),
              some with MYC translocation variants.
              The depth score at relapse reflects
              the DOMINANT SUBCLONE's molecular
              state — which may have shifted from
              diagnosis.
              THE LESSON: The depth score must be
              re-derived at each relapse, not assumed
              to be constant from diagnosis.
              Molecular re-biopsy (FISH, NGS, gene
              expression) at first and subsequent
              relapses is the clinical implementation
              of this principle.
```

---

## SECTION III — t(11;14) MM: THE BCL2 ATTRACTOR

```
FREQUENCY:   ~16–24% of MM. The most common IgH
             translocation subtype.
             More common in men (slight).
             Often associated with:
               Lymphoplasmacytic morphology
               (cells look like lymphocytes more
               than typical plasma cells)
               CD20 expression (rare in MM generally;
               common in t(11;14) MM)
               Low CD56 expression
               Light chain restriction: often λ light
               chain (vs. κ in other subtypes)

PROGNOSIS:   STANDARD RISK.
             Median OS with modern therapy: 5–7+ years.
             HOWEVER: a subset (~30%) of t(11;14) MM
             has co-occurring high-risk features
             (del(17p), amp(1q)) that worsen prognosis.
             Pure t(11;14) without secondary hits:
             among the most indolent MM subtypes.

MOLECULAR MECHANISM:
  t(11;14)(q13;q32):
  The IgH enhancer (on chromosome 14q32) is
  juxtaposed to CCND1 (cyclin D1, on chromosome
  11q13) → CCND1 overexpression.
  CCND1 normally drives G1→S cell cycle transition.
  In normal mature B cells and plasma cells, cyclin
  D1 is NOT expressed — it is a proliferative
  signal that the terminally differentiated plasma
  cell has shut off.
  t(11;14) → CCND1 permanently ON via IgH enhancer
  → constitutive G1→S progression → proliferative
  advantage in the bone marrow niche.
  The BLIMP1 circuit still partially represses
  MYC — which is why t(11;14) has relatively
  restrained proliferation (low proliferation
  index in the TC1 class of GEP profiling).

BCL2 OVEREXPRESSION IN t(11;14):
  The IgH enhancer region at t(11;14) is close
  to the BCL2 locus on chromosome 18 in 3D
  nuclear space.
  When the IgH enhancer from chromosome 14 is
  translocated to chromosome 11, it can also
  drive BCL2 upregulation via 3D chromatin
  contacts (enhancer-hijacking effect).
  ADDITIONALLY: CCND1 itself drives BCL2
  expression by activating BCL2 regulatory
  elements.
  RESULT: t(11;14) myeloma cells are the
  MOST BCL2-DEPENDENT of all MM subtypes.
  The BCL2 dependency creates the venetoclax
  vulnerability.

VENETOCLAX IN t(11;14) MM — THE CLINICAL CONFIRMATION:
  BELLINI TRIAL (venetoclax + bortezomib +
  dexamethasone vs. placebo + Vd):
    t(11;14)-positive subgroup:
    Median PFS: 36.8 months (venetoclax) vs.
    9.3 months (placebo) — HR 0.17, p=0.00041.
    BCL2-high subgroup:
    Median PFS: 30.1 months vs. 9.9 months —
    HR 0.36, p=0.00014.
    UNSELECTED MM: NO OS BENEFIT. Signal toward
    increased mortality (infections) in the
    non-t(11;14) population.
  CANOVA TRIAL (venetoclax + dexamethasone vs.
  pomalidomide + dexamethasone in t(11;14) RRMM):
    Primary endpoint not met (PFS 9.9 vs.
    5.8 months, p=0.237).
    But: post-hoc analyses confirm strong BCL2-
    high subgroup benefit.
    The CANOVA negative result is directly
    analogous to the Phase III failures in MDS:
    t(11;14) is not uniform — some t(11;14) patients
    have additional high-risk features (del(17p),
    amp(1q)) that REDUCE BCL2 dependency and
    venetoclax sensitivity.
    Unselected t(11;14) trial = mixed response.
    Depth score by BCL2 expression level within
    t(11;14) = the correct patient selection tool.
  FRAMEWORK INTERPRETATION:
    Venetoclax in t(11;14) MM is the DIRECT
    PARALLEL TO VENETOCLAX IN CLL (from the
    OrganismCore repository CLL analysis).
    In CLL: BCL2 is elevated (not by IgH
    translocation but by del(13q) and mir-15/16
    loss) → venetoclax dissolves the BCL2-
    dependent false attractor → deep remission.
    In t(11;14) MM: BCL2 is elevated by IgH
    enhancer effect on CCND1/BCL2 → venetoclax
    dissolves the BCL2 component of the false
    attractor → deep remission in BCL2-high
    patients.
    THE SWITCH GENE QUESTION:
    In CLL: BCL2 IS the convergence node.
    In t(11;14) MM: BCL2 may be the convergence
    node, OR CCND1 may be the primary convergence
    node (with BCL2 as a downstream effect).
    The framework analysis of t(11;14) MM must
    distinguish:
      Is CCND1 the top depth-correlating gene
      (the cell cycle driver)?
      Is BCL2 a downstream consequence?
      Or does BCL2 correlate independently with
      depth beyond CCND1?
    If BCL2 r > CCND1 r with depth:
    BCL2 is the convergence node →
    venetoclax is the depth-stratified target.
    If CCND1 r > BCL2 r:
    CCND1 is the convergence node →
    CDK4/6 inhibitors (palbociclib, ribociclib)
    are the depth-stratified targets, not venetoclax.
    The depth correlation between CCND1 and BCL2
    in t(11;14) bulk data is the key analytical
    question.

LENALIDOMIDE + t(11;14):
  An important note: t(11;14) MM has LOWER
  lenalidomide sensitivity than other subtypes.
  Mechanism: CRBN is a required substrate receptor
  for lenalidomide's IKZF1/3 degradation mechanism.
  t(11;14) cells have lower CRBN expression
  levels compared to other MM subtypes.
  The biological explanation: CCND1 overexpression
  alters the IRF4/CRBN regulatory axis.
  CLINICAL IMPLICATION: IMiD-based regimens
  (lenalidomide-based quadruplets) are LESS
  effective in t(11;14) MM.
  Venetoclax-based regimens are MORE effective.
  The depth score should show lower CRBN expression
  at depth in t(11;14) vs. other subtypes — if
  confirmed, this is a Going Further finding that
  explains the lenalidomide resistance mechanistically.

TREATMENT:
  Frontline: Anti-CD38 quadruplet (Dara-VRd
  or Isa-VRd) as standard per 2025 guidelines.
  t(11;14) specific: Venetoclax + bortezomib +
  dexamethasone (VenVd) — off-label but used
  in BCL2-high t(11;14) at relapse.
  Novel: Venetoclax + daratumumab + dexamethasone
  (VenDd) — emerging combination exploiting
  both CD38-targeting AND BCL2 inhibition.
  Active trials: phase III venetoclax combinations
  in t(11;14)-selected patients (ongoing as of 2026).
```

---

## SECTION IV — t(4;14) MM: THE NSD2/H3K36ME2 ATTRACTOR

```
FREQUENCY:   ~11–15% of newly diagnosed MM.
             The most studied high-risk IgH
             translocation subtype.
             Slightly more common in women.

PROGNOSIS:   HIGH RISK (intermediate-to-poor).
             Median OS without bortezomib: ~2 years.
             Median OS with bortezomib: ~4–5 years
             (bortezomib PARTIALLY overcomes the
             poor prognosis).
             IMPORTANT 2024 FINDING: The t(4;14)
             risk is HETEROGENEOUS based on NSD2
             breakpoint location within the NSD2 gene:
               "No disruption" or "early disruption"
               (5' UTR of NSD2): NSD2 protein is
               fully expressed. Median OS ~75 months.
               "Late disruption" (coding region
               of NSD2): NSD2 function is partially
               lost. Median OS ~29 months.
             This is DEPTH STRATIFICATION WITHIN
             A SUBTYPE — exactly what the framework
             depth score should capture.

MOLECULAR MECHANISM:
  t(4;14)(p16;q32):
  TWO genes are dysregulated simultaneously:
  1. FGFR3 (fibroblast growth factor receptor 3,
     on chromosome 4p16) — brought under IgH
     enhancer control → FGFR3 overexpression
     in ~70% of t(4;14) cases.
     BUT: FGFR3 is lost in ~25–30% of t(4;14)
     MM at progression (FGFR3 loss is not lethal
     to the t(4;14) clone — NSD2 is the critical
     oncogene, not FGFR3).
  2. NSD2/MMSET/WHSC1 (nuclear receptor binding
     SET domain protein 2, on chromosome 4p16) —
     IgH enhancer → NSD2 overexpression in
     VIRTUALLY ALL t(4;14) cases.
     NSD2 is the methyltransferase responsible
     for H3K36me2 (di-methylation of histone H3
     at lysine 36) — an active chromatin mark
     that promotes gene transcription.
     NSD2 OVEREXPRESSION:
       H3K36me2 globally increased across the
       t(4;14) myeloma cell genome.
       Increased H3K36me2 → MORE PERMISSIVE
       chromatin → more transcription of
       oncogenic programmes, proliferative genes.
       COMPETITION WITH H3K27me3:
         H3K36me2 and H3K27me3 are mutually
         exclusive on the same nucleosome.
         NSD2 overexpression → more H3K36me2 →
         LESS H3K27me3 at the same loci →
         LESS EZH2-MEDIATED REPRESSION.
         In t(4;14) myeloma: EZH2 activity is
         LESS EFFECTIVE than in normal plasma
         cells or other MM subtypes, because
         NSD2-driven H3K36me2 COMPETES with
         EZH2-driven H3K27me3 for the same
         chromatin.
         This means: EZH2 inhibitors (tazemetostat)
         would NOT be the primary target in
         t(4;14) MM — NSD2 is the better target.
         EZH2 inhibition in t(4;14) MM would
         further reduce H3K27me3 at loci already
         lacking it, with unpredictable effects.
         NSD2 inhibition (RK-552 and related
         compounds in Phase I/II 2024–2025) is
         the mechanistically correct approach.
  NSD2 METABOLIC MECHANISM (2024 finding):
     NSD2 consumes S-adenosylmethionine (SAM)
     as the methyl donor for H3K36me2.
     Excess NSD2 activity depletes SAM →
     creatine synthesis is impaired (creatine
     synthesis requires SAM) →
     Adenylate kinase 2 (AK2) is upregulated
     to compensate (AK2 regenerates ADP from
     AMP to maintain energy balance under
     creatine-synthesis impairment) →
     t(4;14) cells are uniquely DEPENDENT ON AK2.
     AK2 dependence → sensitivity to DNA
     replication stress → bortezomib sensitivity
     (bortezomib increases unfolded protein stress
     which further taxes the AK2-dependent
     energy balance in t(4;14) cells).
     THIS IS THE MECHANISTIC BASIS FOR BORTEZOMIB
     SENSITIVITY IN t(4;14) MM — previously not
     fully understood.
     FRAMEWORK INTERPRETATION: The NSD2 → SAM
     depletion → AK2 dependence chain is a
     metabolic false attractor maintained by
     epigenetic rewiring.
     NSD2 inhibition → SAM restored → AK2
     dependence removed → bortezomib sensitivity
     may be partially bypassed (resistance risk).
     Conversely: NSD2 inhibition + bortezomib
     combination: NSD2 inhibition reduces H3K36me2,
     partially restores H3K27me3 silencing of
     oncogenic programmes, AND SIMULTANEOUSLY
     reduces the SAM-depletion-AK2 dependence
     that currently makes bortezomib work.
     Net effect: antagonism? Or synergy?
     This is a structural novel prediction for
     the t(4;14) depth analysis.
  BORTEZOMIB RESISTANCE IN t(4;14):
     HRP2 (HDGFRP2) is a Polycomb complex member
     that normally COUNTERS NSD2 — it activates
     a MINA/ROX1 demethylase pathway that reduces
     H3K36me2 and raises H3K27me3 at specific loci.
     HRP2 LOSS (LOH, deletion, or silencing) →
     more permissive H3K36me2 landscape →
     reduced ER stress sensing (fewer UPR genes
     exposed to repression) →
     bortezomib resistance.
     EZH2 inhibition (tazemetostat) + bortezomib
     in HRP2-loss t(4;14) MM can re-sensitise
     cells to proteasome inhibition (because
     tazemetostat, by reducing H3K27me3,
     counter-intuitively re-opens ER stress
     genes that bortezomib depends on).
     This is a special case where tazemetostat
     is relevant in t(4;14) MM — but ONLY in
     the HRP2-loss context.
     FRAMEWORK PATTERN 4 APPLIED TO t(4;14):
     The EZH2 direction in t(4;14) MM is CONTEXT-
     SPECIFIC within the subtype:
     General t(4;14): EZH2 activity REDUCED
     (NSD2 competes with H3K27me3).
     EZH2 inhibitor NOT appropriate.
     t(4;14) with HRP2 loss: EZH2 may be
     PARADOXICALLY relevant for re-sensitisation.
     EZH2 direction must be measured in the
     HRP2-loss subgroup separately.

TREATMENT:
  Frontline: Bortezomib-containing quadruplet
  (Dara-VRd) — bortezomib particularly important
  due to t(4;14) bortezomib sensitivity.
  Relapse: Carfilzomib-based regimens (KRd,
  DKd) — carfilzomib shows activity.
  Novel: NSD2 inhibitors in Phase I/II (2024–2025).
  Autologous SCT: strongly recommended if eligible.
  Tandem transplant: considered for some high-risk
  t(4;14) (especially with additional high-risk
  features like del(17p) or amp(1q)).
```

---

## SECTION V — MAF SUBTYPES: t(14;16) AND t(14;20)

```
FREQUENCY:   t(14;16): ~5–7% of MM.
             t(14;20): ~1–2% of MM.
             Together: ~7–9% of MM.
             Rarer but among the highest-risk
             IgH translocation subtypes.

PROGNOSIS:   HIGH RISK.
             t(14;16): Median OS ~3–4 years
             with modern therapy.
             t(14;20): Similar to t(14;16).
             Both are classified as standard
             high-risk cytogenetics.

MOLECULAR MECHANISM:
  t(14;16)(q32;q23):
  IgH enhancer → MAF (c-MAF transcription factor)
  overexpression.
  MAF is a large Maf family bZIP transcription factor.
  NORMAL MAF FUNCTION:
    Regulates lens differentiation (normally
    expressed in lens and some immune cells).
    In B cells/plasma cells: MAF is NOT normally
    expressed at high levels.
  MAF OVEREXPRESSION IN t(14;16):
    MAF activates:
      CCND2 (cyclin D2) — the cyclin D isoform
      used in t(14;16) (vs. CCND1 in t(11;14),
      CCND2 in hyperdiploid, CCND3 in others).
      ITGB7 — integrin beta-7, promotes bone
      marrow homing and adhesion.
      DKK1 — Dickkopf 1, WNT inhibitor that
      drives osteoblast suppression (worsening
      bone disease).
      AKT/PI3K pathway genes — proliferative
      survival signalling.
      IGF1R — insulin-like growth factor receptor,
      promotes survival.
    MAF co-activates BCL6 targets — the t(14;16)
    myeloma cell has PARTIAL GC IDENTITY preserved
    (the MAF-driven programme includes some BCL6-
    driven genes), explaining the partial GC-like
    phenotype of MAF-subtype cells.
  t(14;20)(q32;q12):
  IgH enhancer → MAFB (B-Maf transcription factor)
  overexpression.
  MAFB normally expressed in monocytes/macrophages.
  In plasma cells: MAFB is normally absent.
  MAFB OVEREXPRESSION:
    Similar to MAF but with:
    Monocytic gene activation (MAFB is a monocyte
    TF — t(14;20) myeloma cells have partial
    monocyte/macrophage gene expression).
    CCND2 activation (same as MAF).
    AKT pathway activation.
  EZH2 IN MAF SUBTYPES:
    MAF and MAFB activate WNT target genes including
    some Polycomb-regulated targets.
    EZH2 in t(14;16)/t(14;20) MM is context-dependent
    but may be RELATIVELY ELEVATED compared to
    normal plasma cells — MAF activates proliferative
    programmes that EZH2 normally silences, and
    the cell must maintain EZH2 to prevent further
    identity loss (the MAF programme would drive
    the cell AWAY from plasma cell identity toward
    a more dedifferentiated state).
    EZH2 in MAF subtypes: ELEVATED (predicted).
    EZH2 inhibitors MIGHT be relevant in MAF-subtype
    MM — more so than in t(4;14), less so than in
    BRCA (where EZH2 is the dominant convergence node).
    This is a structural prediction to test in MM-S2.

MAF-SPECIFIC DRUG TARGETS:
  No MAF-specific inhibitor exists clinically.
  MAF is a transcription factor — historically
  considered undruggable.
  2024 approaches:
    Indirect targeting via cell cycle (CCND2 → CDK4/6
    inhibitors; palbociclib active in preclinical MAF MM).
    AKT/PI3K inhibitors (given AKT pathway activation).
    Bortezomib sensitivity: MAF subtypes are moderately
    bortezomib-sensitive (not as markedly as t(4;14)).
  Novel: First-in-class MAF inhibitors (targeting the
  bZIP domain) are in early-stage discovery.
  Daratumumab/BCMA-targeting: Effective regardless of
  MAF status (because these target the retained plasma
  cell identity — BCMA, CD38).

TREATMENT:
  Frontline: Anti-CD38 quadruplet (Dara-VRd or
  Isa-VRd) — same as other MM subtypes.
  Bortezomib inclusion important.
  Autologous SCT: strongly recommended.
  Maintenance: extended lenalidomide ±
  bortezomib; continuous therapy may offset
  the high-risk biology.
  Relapse: BCMA-targeting bispecifics (teclistamab)
  and CAR-T (cilta-cel, ide-cel) —
  BCMA is retained in MAF-subtype MM.
  Clinical trial enrollment recommended
  given the poor standard-of-care outcomes.
```

---

## SECTION VI — HYPERDIPLOIDY: THE MOST COMMON MM SUBTYPE

```
FREQUENCY:   ~50–55% of newly diagnosed MM.
             The MOST COMMON MM subtype.
             More common in older patients.
             Slightly more common in men.

PROGNOSIS:   STANDARD-TO-FAVORABLE RISK.
             Median OS with modern therapy: 6–10+ years.
             Among the best-prognosis MM subtypes.
             EXCEPTION: Hyperdiploid MM with co-occurring
             del(17p), amp(1q), or MYC translocation =
             "double-hit" → much worse prognosis.

MOLECULAR MECHANISM:
  DEFINITION: Gain of multiple odd-numbered
  chromosomes — typically 3, 5, 7, 9, 11, 15,
  19, and 21.
  Usually results in 48–75 chromosomes
  (normal = 46).
  This is POLYSOMY (multiple copies of entire
  chromosomes) rather than focal amplification.
  MECHANISM:
    Hyperdiploidy arises from abnormal mitosis
    during the GC B cell stage — possibly from
    impaired spindle assembly checkpoint.
    The chromosomal gains are not random:
    The preferential gain of ODD-NUMBERED
    chromosomes (and not even-numbered ones)
    suggests a specific mitotic mechanism.
    The gained chromosomes bring:
    Extra copies of CCND2 (chromosome 12) —
    Cyclin D2 overexpression drives G1→S.
    Extra copies of CDK4 (chromosome 12) —
    Further cell cycle promotion.
    Extra copies of IRF4 (chromosome 6) —
    Increased IRF4 expression in some hyperdiploid
    subtypes (but IRF4 is on chr 6, even-numbered;
    chr 6 gain is less common in canonical
    hyperdiploidy — this depends on the specific
    trisomy pattern).
  NO IgH TRANSLOCATION:
    Hyperdiploid MM is characterised by the
    ABSENCE of IgH translocations.
    The cell cycle driver is CCND2 (from chr 12
    polysomy) rather than CCND1 (from t(11;14)).
    This is a key diagnostic and biological
    distinction.
  LESS BCL2-DEPENDENT:
    Without the t(11;14)/CCND1/BCL2 mechanism,
    hyperdiploid MM cells are typically LESS
    BCL2-dependent and MORE MCL1-dependent.
    Venetoclax is generally NOT effective in
    non-t(11;14) MM (BELLINI trial — increased
    mortality in unselected non-t(11;14) MM).
    MCL1 inhibitors (AZD5991, AMG176) are the
    BCL2-family targets in hyperdiploid MM —
    but MCL1 inhibitors have cardiotoxicity
    concerns and are still in clinical trials.
  EZH2 IN HYPERDIPLOID MM:
    Hyperdiploid MM does not have a specific
    epigenetic driver analogous to NSD2 in t(4;14).
    EZH2 in hyperdiploid MM: context-dependent.
    If CCND2 is the primary driver and IRF4
    is intact: EZH2 may be in an intermediate
    state — not strongly elevated, not reduced.
    Additional hits (del(17p), amp(1q)) will
    modulate EZH2 direction.
    The depth score in hyperdiploid MM should
    reflect primarily: CCND2 level, proliferation
    index, BM plasma cell %, and the presence
    of secondary high-risk features.

TREATMENT:
  Frontline: Anti-CD38 quadruplet —
  Dara-VRd or Isa-VRd per 2025 guidelines.
  Lenalidomide-based maintenance (standard —
  hyperdiploid MM has good CRBN/lenalidomide
  sensitivity because it lacks the CCND1/CRBN
  interference of t(11;14)).
  Autologous SCT: recommended for eligible patients.
  Relapse: BCMA-targeting agents (teclistamab,
  cilta-cel) effective — BCMA retained.
  Venetoclax: NOT effective in BCL2-low
  hyperdiploid MM. Do NOT use.
  MCL1 inhibitors: under investigation.
```

---

## SECTION VII — HIGH-RISK CYTOGENETIC LESIONS

```
THREE LESIONS THAT CONFER HIGH RISK REGARDLESS
OF WHICH SUBTYPE THEY CO-OCCUR WITH:

═══════════════════════════════════════════════════════
DELETION 17p (del(17p13)) — TP53 LOSS
═══════════════════════════════════════════════════════

FREQUENCY:   ~10% of newly diagnosed MM.
             ~25–30% of RRMM (selected under
             treatment pressure).
             ~70–80% of del(17p) MM also has
             TP53 mutation on the retained allele
             → biallelic TP53 inactivation
             (analogous to MDS-biTP53).

PROGNOSIS:   VERY POOR.
             Median OS ~2–3 years even with
             modern quadruplet therapy.
             The most adverse single cytogenetic
             feature in MM.
             "Double-hit" MM: del(17p) + one other
             high-risk feature (amp(1q), t(4;14),
             t(14;16)) → median OS <2 years.
             "Triple-hit" MM: del(17p) + two others →
             median OS <1.5 years in some series.

MECHANISM:
  TP53 LOSS IN MM (same principles as MDS-biTP53):
    p53 is the genome guardian.
    TP53 haploinsufficiency → reduced response to
    DNA damage → accumulation of additional
    mutations.
    Biallelic TP53 loss → no p53 function →
    maximum chromosomal instability.
    In MM with del(17p): the myeloma clone can
    accumulate MYC translocation, amp(1q),
    immunoglobulin-driven translocation PLUS
    TP53 loss → DEEPEST ATTRACTOR IN MM.
    Del(17p) MM is MM at its greatest depth —
    the augmented survival/proliferation programme
    has become maximally dysregulated, most
    resistant, and most unstable.
  THE MULTI-HIT TP53 CONCEPT IN MM (2025):
    Like MDS-biTP53, multi-hit TP53 MM (both
    alleles inactivated) is worse than
    monoallelic del(17p) alone.
    A 2025 paper confirms: multi-hit TP53 in MM
    confers the poorest survival in the novel
    agent era.
    The same patient selection principle applies:
    monoallelic del(17p) vs. biallelic TP53
    inactivation are different disease states,
    not one.

TREATMENT:
  Frontline: Intensified quadruplet — Dara-VRd
  with planned tandem autologous SCT considered.
  Bortezomib STRONGLY included (PI therapy
  partially overcomes TP53-null proliferative
  advantage by increasing unresolvable
  proteotoxic stress regardless of TP53 status).
  Dual maintenance: bortezomib + lenalidomide.
  CAR-T / bispecifics: being moved to earlier
  lines for del(17p) MM given poor outcomes
  with standard therapy.
  Clinical trial enrollment: ESSENTIAL.
  No standard therapy produces acceptable
  long-term outcomes in del(17p)/biTP53 MM.

═══════════════════════════════════════════════════════
GAIN/AMPLIFICATION OF 1q21 (amp(1q))
═══════════════════════════════════════════════════════

FREQUENCY:   ~40% of newly diagnosed MM.
             ~60–70% of RRMM.
             The most common chromosomal
             abnormality in MM overall.
             Increases markedly with relapse —
             therapy-driven clonal selection
             for 1q-amplified subclones.

PROGNOSIS:   POOR (especially ≥4 copies of 1q21).
             3 copies (gain): intermediate-risk.
             ≥4 copies (amplification): high-risk,
             approaching del(17p) in some series.

MECHANISM:
  1q21 contains CKS1B (CDC28 protein kinase
  regulatory subunit 1B) — a cell cycle regulator
  that promotes S-phase entry and degradation
  of CDK inhibitors (p27, p21).
  CKS1B overexpression (3–4+ copies) →
  accelerated cell cycle entry → proliferative
  advantage.
  1q21 ALSO contains MCL1 (myeloid cell
  leukaemia 1) — the main anti-apoptotic BCL2
  family member in non-t(11;14) MM.
  MCL1 amplification → increased anti-apoptotic
  buffer → resistance to proteasome inhibitor-
  induced apoptosis (bortezomib/carfilzomib
  induce apoptosis partly via NOXA which
  displaces MCL1; if MCL1 is overexpressed,
  NOXA cannot efficiently displace it).
  ALSO: IL-6 receptor signalling → STAT3 →
  MCL1 transcription amplifies the 1q21-copy
  number-driven MCL1 elevation.
  EZH2 AND amp(1q):
    EZH2 is not directly on 1q21 but MCL1 is.
    MCL1 overexpression in amp(1q) MM means the
    apoptosis axis is resistant to BH3 mimetics.
    The depth score in amp(1q) MM should show
    high CKS1B and MCL1 as the dominant
    depth-positive genes.
    EZH2 in amp(1q) MM: context-dependent —
    amplification of 1q does not directly
    regulate EZH2. If t(4;14) co-occurs with
    amp(1q): NSD2 drives H3K36me2 and reduces
    EZH2 activity (as in Section IV). If t(11;14)
    co-occurs with amp(1q): BCL2 + MCL1 are BOTH
    elevated → dual anti-apoptotic protection.
    This t(11;14) + amp(1q) combination requires
    BCL2 + MCL1 co-inhibition (venetoclax + MCL1
    inhibitor) — neither alone is sufficient.
    NOVEL PREDICTION: t(11;14) + amp(1q) MM
    (double BCL2 + MCL1 elevation) requires
    venetoclax + MCL1 inhibitor combination.
    Neither venetoclax alone (BELLINI) nor MCL1
    inhibitor alone will produce durable response.
    The BELLINI partial response data may
    partly reflect amp(1q) contamination in the
    t(11;14) cohort — patients with t(11;14) +
    amp(1q) had MCL1 compensation reducing
    venetoclax response.
    DEPTH SCORE APPLICATION: Can the depth score
    in t(11;14) MM identify the amp(1q) co-occurring
    subset by MCL1/CKS1B elevation WITHIN the
    BCL2-high t(11;14) population? This would
    define the venetoclax-resistance signature
    within t(11;14) MM and provide the molecular
    selection criterion for venetoclax + MCL1
    inhibitor combination trials.

TREATMENT:
  Frontline: Daratumumab-based quadruplet
  (Dara-VRd).
  Carfilzomib preferred over bortezomib in some
  amp(1q) series (carfilzomib generates more
  irreversible UPR stress that MCL1 cannot
  fully buffer).
  Novel: MCL1 inhibitors in Phase I/II — AZD5991,
  AMG176. Cardiotoxicity monitoring required.
  Daratumumab + MCL1 inhibitor combinations
  in amp(1q)-selected trials: emerging.
  Allo-SCT: considered for eligible patients with
  amp(1q) + del(17p) (double-hit).

═══════════════════════════════════════════════════════
MYC ACTIVATION (LATE EVENT, ~70% OF MM AT PROGRESSION)
═══════════════════════════════════════════════════════

FREQUENCY:   ~40–45% by structural MYC translocation
             at diagnosis.
             ~60–70% by native MYC enhancer
             gain-of-function (no structural variant).
             MYC is ACTIVE in the vast majority of
             established symptomatic MM.
             The mechanism of MYC activation
             determines the target.

MECHANISM:
  STRUCTURAL MYC TRANSLOCATION:
    MYC brought near IgH, IgL, or IgK enhancers.
    Strong transcriptional drive → very high MYC.
    Often a SECONDARY event (on top of the primary
    IgH translocation or hyperdiploidy).
    In MGUS/SMM: MYC translocation is the
    TRANSITION EVENT that converts the slowly
    expanding pre-malignant clone into
    proliferative MM.
  NATIVE ENHANCER GAIN-OF-FUNCTION (2024):
    IRF4, SPIB, and other plasma cell TFs
    drive a plasma cell-specific MYC enhancer
    to abnormally high activity.
    This occurs WITHOUT any structural change.
    Detectable only by chromatin accessibility
    assays (ATAC-seq), histone ChIP-seq, or
    high-resolution 3D genomics.
    NOT detectable by FISH (standard clinical test).
    This explains why ~30% of MM have MYC-high
    gene expression without a detectable MYC
    translocation by FISH — the enhancer is
    epigenetically activated, not structurally rearranged.
  MYC AND BLIMP1 ANTAGONISM:
    Normal BLIMP1 should repress MYC.
    In MM: the MYC activation exceeds BLIMP1's
    repressive capacity → both BLIMP1 and MYC
    are simultaneously expressed.
    The depth axis in MM directly measures the
    degree of MYC escape from BLIMP1 repression:
    BLIMP1 expression (should be constant or high
    in plasma cells) vs. MYC expression (should
    be low; its elevation = depth signal).
    DEPTH-POSITIVE IN MM: MYC, CDKs, proliferation
    markers, MCL1, BCL2 (in t(11;14)).
    DEPTH-NEGATIVE IN MM: BLIMP1 target genes
    related to normal PC function (J-chain,
    some secretory pathway genes that become
    relatively downregulated in high-proliferation
    state), or apoptosis sensors (BIM).

IRF4-MYC CO-ADDICTION (THE CORE ONCOGENIC LOOP):
  IRF4 → activates MYC (in MM, IRF4 has
    gained the ability to activate MYC; in
    normal PCs, this connection is weaker or
    absent because BLIMP1 counters it).
  MYC → activates IRF4 (MYC drives IRF4 transcription).
  IRF4 and MYC mutually reinforce each other
  in a feed-forward loop.
  This is the CORE SURVIVAL LOOP of myeloma cells.
  LENALIDOMIDE KILLS THE LOOP:
  Lenalidomide → cereblon → degrades IKZF1/3
  (Ikaros/Aiolos) → IKZF1/3 normally activate
  IRF4 transcription → IRF4 drops → MYC drops
  → myeloma cell death.
  This is the LENALIDOMIDE MECHANISM IN MM:
  not direct BCL2 inhibition (as in venetoclax),
  not direct epigenetic reprogramming (as in
  tazemetostat), but COLLAPSE OF THE IRF4-MYC
  CO-ADDICTION LOOP via IKZF1/3 degradation.
  WHY t(11;14) MM IS LENALIDOMIDE-RESISTANT:
  t(11;14) MM has lower CRBN expression.
  Without sufficient CRBN: lenalidomide cannot
  efficiently degrade IKZF1/3 → IRF4/MYC loop
  remains intact → lenalidomide resistance.
  This connects to the depth score: in t(11;14) MM,
  CRBN expression should be DEPTH-NEGATIVE
  (lower at depth) — lower CRBN = deeper disease
  = more lenalidomide resistant.
  If confirmed, this is a DEPTH-STRATIFIED
  LENALIDOMIDE RESISTANCE PREDICTOR within t(11;14).
```

---

## SECTION VIII — EZH2 DIRECTION IN MM SUBTYPES

```
FRAMEWORK PATTERN 4: THE EZH2 DIRECTION RULE
APPLIED TO MM.

This is more complex in MM than in any other
cancer in the repository because MM has five
distinct molecular classes with different
EZH2 mechanisms.

═══════════════════════════════════════════════════════
SUMMARY TABLE: EZH2 DIRECTION BY MM SUBTYPE
═══════════════════════════════════════════════════════

SUBTYPE           EZH2 DIRECTION     MECHANISM        EZH2i
                  (PREDICTED)                         APPROPRIATE?
──────────────────────────────────────────────────────────────────
t(11;14)          ELEVATED           BCL2/CCND1       POSSIBLY
                  (relative to       drives enhanced  (if EZH2 is
                  normal LLPC)       H3K27me3 at      maintaining
                                     B cell genes to  the plasma
                                     maintain plasma  cell state
                                     cell fate while  while BCL2
                                     BCL2 increases   is the
                                     survival         survival node)
                                                      TARGET: BCL2
                                                      FIRST, not EZH2

t(4;14)           REDUCED            NSD2 drives      NO (in
                  (activity-wise,    H3K36me2         general)
                  not expression-    competing with   YES only in
                  wise — EZH2        H3K27me3;        HRP2-loss
                  protein may be     EZH2-placed      subgroup
                  present but        marks displaced  (re-sensitises
                  functionally       by H3K36me2      to bortezomib)
                  reduced)           NSD2 inhibitor   TARGET: NSD2
                                     is the target    FIRST

t(14;16)/t(14;20) RELATIVELY         MAF/MAFB drive   POSSIBLY
                  ELEVATED           WNT/proliferative (MAF requires
                                     genes; EZH2      EZH2 to
                                     likely needed to prevent full
                                     maintain plasma  GC-reversion)
                                     cell identity    TEST IN DATA
                                     against MAF-
                                     driven
                                     dedifferentiation

Hyperdiploidy     INTERMEDIATE       No specific      UNCERTAIN
                  (near-normal       EZH2 mechanism   Depends on
                  range)             unless secondary secondary
                                     hits present     features
                                                      TEST IN DATA

del(17p)/         VARIABLE —         TP53 loss allows GENERALLY
double-hit/       CONTEXT-           MYC-driven       NOT
triple-hit        DEPENDENT          H3K27me3         appropriate
                                     reprogramming    (circuit broken,
                                     (may or may not  too unstable
                                     drive EZH2 up)   for single-
                                                      target therapy)

═══════════════════════════════════════════════════════
THE MM-SPECIFIC EZH2 LESSON:
═══════════════════════════════════════════════════════

In BRCA (OrganismCore BRCA analysis):
EZH2 elevated = gain-of-function lock on luminal
identity → tazemetostat is the target.

In MM: EZH2 is elevated in some subtypes
(t(11;14), MAF), reduced in others (t(4;14)),
and intermediate in others (hyperdiploid).
The wrong EZH2 assumption in MM could produce:
  EZH2 inhibitor in t(4;14): further reduces
  H3K27me3 that is already competed away by
  NSD2 → unpredictable gene derepression →
  potential worsening.
  EZH2 inhibitor in t(11;14): may disrupt the
  EZH2-maintained plasma cell identity that is
  NEEDED to maintain BCMA/CD38/IRF4 expression
  → loss of drug targets (if BCMA falls after
  EZH2 inhibition, BCMA-targeting agents lose
  their target) → catastrophic therapeutic
  failure.
EZH2 DIRECTION MUST BE MEASURED, NOT ASSUMED.
The depth score analysis for each MM subtype
will include the EZH2 correlation as one of the
first outputs — confirming or correcting these
structural predictions.
```

---

## SECTION IX — THE IRF4 NODE: MASTER CONVERGENCE CANDIDATE

```
THE CANDIDATE FOR MASTER CONVERGENCE NODE IN MM:
IRF4.

In every cancer in the repository, the
convergence node is the gene whose expression
most consistently maintains the false attractor
and whose inhibition/modulation would dissolve
it.

WHY IRF4 IS THE MASTER CANDIDATE IN MM:

1. IRF4 IS REQUIRED FOR ALL MM SUBTYPES.
   Every MM subtype — t(11;14), t(4;14), t(14;16),
   hyperdiploid, del(17p) — retains IRF4 high
   expression.
   IRF4 knockdown is lethal to MM cells regardless
   of cytogenetic subtype.
   This is "non-oncogene addiction" — MM cells
   are addicted to IRF4 even though IRF4 is not
   itself a classically mutated oncogene.

2. IRF4 DRIVES THE MYC FEED-FORWARD LOOP.
   IRF4 → MYC → IRF4 (mutual activation).
   This loop is the CORE ONCOGENIC CIRCUIT
   superimposed on the retained plasma cell identity.
   Disrupting IRF4 collapses the loop.

3. IRF4 ACTIVATES BLIMP1.
   BLIMP1 represses PAX5, BCL6, MYC.
   The normal circuit is: IRF4 → BLIMP1 →
   MYC repressed.
   In MM: IRF4 → BLIMP1 is maintained BUT
   simultaneously IRF4 → MYC is also maintained
   (IRF4 has gained the MYC activation function
   that normal plasma cell IRF4 lacks).
   The myeloma cell is in a paradox state where
   BLIMP1 tries to repress MYC but IRF4 is
   simultaneously activating it.
   The depth of the false attractor is determined
   by the balance between BLIMP1-mediated MYC
   repression and IRF4-mediated MYC activation.
   At low depth: BLIMP1 wins, MYC low, more
   plasma cell-like, slower proliferation.
   At high depth: IRF4/MYC loop wins, MYC high,
   high proliferation, bortezomib sensitive
   (high Ig production = high UPR stress = high
   dependency on proteasome clearance).

4. LENALIDOMIDE WORKS BY INDIRECTLY TARGETING IRF4.
   The entire lenalidomide mechanism (CRBN →
   IKZF1/3 degradation → IRF4 reduction → MYC
   reduction) is the clinical confirmation that
   IRF4 is the convergence node.
   Lenalidomide does not directly inhibit IRF4.
   It depletes the TFs (IKZF1/3) that drive
   IRF4 transcription.
   The success of lenalidomide in MM (it is
   effective across most subtypes) is the
   CROSS-SUBTYPE CONFIRMATION of the IRF4 node.

5. IRF4 DEPTH CORRELATION IS THE PREDICTED
   ANALYSIS OUTPUT.
   In the depth score analysis:
   IRF4 should be among the top depth-correlating
   genes in the BULK MM analysis.
   The existing OrganismCore MM analysis may
   have already identified IRF4 or an IRF4 target
   as a top depth correlate.
   The subtype series will determine:
   Does IRF4 correlation with depth vary by
   molecular subtype?
   If IRF4 is equally depth-correlated in all
   subtypes → IRF4 is the universal MM convergence node.
   If IRF4 depth correlation is lower in t(11;14)
   (which shows BCL2/CCND1 as primary driver)
   and higher in t(4;14) or MAF subtypes →
   IRF4 may be a secondary node in some subtypes.
   This is the key analytical output of MM-S1
   through MM-S4.

DIRECT IRF4 TARGETING:
  IRF4 is a transcription factor — historically
  undruggable.
  2024 approaches:
    IRF4 degraders (PROTAC-based) — in development.
    CELMoDs (next-gen cereblon modulators) —
    more potent IKZF1/3 degradation → more
    complete IRF4 depletion:
      CC-220 (iberdomide) — clinical trials
      CC-92480 (mezigdomide) — clinical trials
      Both show activity after lenalidomide
      resistance (by more efficiently degrading
      IKZF1/3 even when CRBN is partially
      downregulated).
    BCL6 inhibition — indirect IRF4 approach:
    In MAF subtypes with BCL6 activity, BCL6
    inhibition + IRF4 disruption may collapse
    the GC-like component of the MAF programme.
```

---

## SECTION X — THE BONE MARROW NICHE AS DEPTH CO-DRIVER

```
MM IS UNIQUELY MICROENVIRONMENT-DEPENDENT AMONG
ALL CANCERS IN THE REPOSITORY.

Myeloma cells CANNOT survive outside the bone
marrow niche in vivo — ex vivo cultures die
rapidly without stromal support.
This extreme niche-dependence means the depth
score in MM must account for niche-driven
survival signals that are not intrinsic to the
myeloma cell itself.

KEY NICHE SURVIVAL SIGNALS:

IL-6 AXIS:
  Bone marrow stromal cells (BMSCs) produce IL-6.
  IL-6 → JAK1/2/TYK2 → STAT3 phosphorylation →
    STAT3 activates:
      MCL1 (anti-apoptotic) — the PRIMARY
      IL-6-driven survival gene.
      MYC (proliferative).
      CCND1/D2 (cell cycle).
      VEGF (angiogenesis — MM bone lesions are
      well-vascularised).
  MYELOMA CELL ADHESION → MORE IL-6:
    VLA-4 (integrin α4β1) on myeloma cells binds
    VCAM-1 on BMSCs → contact-mediated IL-6
    upregulation.
    This adhesion-dependent IL-6 loop is the
    "cell adhesion-mediated drug resistance"
    (CAM-DR) phenomenon — myeloma cells adherent
    to stroma are MORE resistant to ALL standard
    therapies, because IL-6 → STAT3 → MCL1
    provides an additional anti-apoptotic buffer.
  CAM-DR AND THE DEPTH SCORE:
    Gene expression data from CD138+ cells
    (the standard MM research cell type) includes
    the INTRINSIC myeloma cell programme.
    It does NOT include the stromal IL-6 signal
    (which is paracrine, not intrinsic).
    The depth score may UNDERESTIMATE depth in
    patients with high-adhesion, high-IL-6 tumours.
    CD138+ gene expression should show STAT3
    target gene elevation as an INDIRECT read-out
    of the IL-6/CAM-DR state.
    STAT3 target gene elevation in CD138+ cells
    = proxy for niche-driven IL-6 depth.
    If MCL1 expression in CD138+ cells correlates
    with depth: MCL1 is measuring BOTH intrinsic
    (amp(1q)) AND extrinsic (IL-6-driven) survival.
    The framework cannot fully disentangle these
    from bulk CD138+ data — this is a known
    limitation to record in the analysis documents.

BCMA/BAFF/APRIL AXIS:
  BCMA (B cell maturation antigen, TNFRSF17) on
  myeloma cells binds BAFF (TNFSF13B) and APRIL
  (TNFSF13) from BMSCs → PI3K/AKT/NF-κB →
  MCL1, BCL2 upregulation + plasma cell survival.
  BCMA signalling is the PRIMARY PLASMA CELL
  IDENTITY SURVIVAL SIGNAL — it is the signal
  that allows long-lived plasma cells to survive
  in the bone marrow niche for years.
  Myeloma cells hijack this signal for immortality.
  BCMA as a drug target USES the niche dependency
  AGAINST the myeloma cell:
  Teclistamab/bispecifics: redirects T cells
  toward BCMA-expressing myeloma cells.
  CAR-T cells: engineered to kill all BCMA+
  cells (myeloma) while sparing non-BCMA normal
  cells (most tissues don't express BCMA).
  Belantamab mafodotin: anti-BCMA ADC.
  ALL BCMA-TARGETING AGENTS WORK ACROSS ALL
  MM SUBTYPES because BCMA is the retained
  plasma cell niche-survival signal — it is
  depth-independent.
  The depth score does NOT predict BCMA-agent
  response. BCMA EXPRESSION LEVEL predicts it.
  A separate BCMA expression analysis should
  accompany the depth score analysis for MM.

RANKL/OPG AXIS (BONE DISEASE):
  Myeloma cells produce RANKL (receptor activator
  of NF-κB ligand) and suppress OPG (osteoprotegerin,
  the RANKL decoy receptor).
  RANKL → osteoclast differentiation → bone
  resorption → calcium release + growth factor
  release → myeloma cell growth (vicious cycle).
  DKK1 from myeloma cells → inhibits WNT in
  osteoblasts → bone repair suppressed.
  CLINICAL INTERVENTION:
  Denosumab (anti-RANKL monoclonal) — FDA approved
  for MM bone disease. Reduces SREs (skeletal
  related events) and may have direct anti-myeloma
  effects (interrupting the vicious cycle).
  Zolendronic acid (bisphosphonate) — kills
  osteoclasts → reduces bone disease + potential
  direct anti-myeloma effect.
  The depth score's bone disease connection:
  DKK1 expression and RANKL expression in CD138+
  cells should be depth-positive (higher at depth)
  → deeper myeloma = more bone disease = worse
  prognosis. This is testable.
```

---

## SECTION XI — DATA AVAILABILITY SUMMARY

```
Dataset         Accession     Samples               Notes      Power
───────────────────────────────────────────────────���─────────────────
MMRF CoMMpass  dbGaP         ~1,000 newly           RNA-seq    VERY
(IA20)         phs000748     diagnosed MM patients  + WES      HIGH
                             + longitudinal         + FISH
                             at relapse             + outcome
                             CD138+ cells sorted    data
                             NOTE: access requires
                             dbGaP application but
                             academic access widely
                             granted.

GSE167968      GEO           CD138+ MM cells        RNA-seq    HIGH
                             Multiple patients +    purified
                             response data

GSE269875      GEO (2024)    MM + normal plasma     Spatial    MOD
                             cells BM spatial       + bulk
                             transcriptomics
                             5 normal PC +
                             28 MM samples

GSE223060      GEO (2023)    41 MM patients, 53 BM  scRNA-seq  HIGH
                             samples + bulk RNA-seq + bulk
                             Normal BM reference    RNA-seq

HOVON65/GMMG-  European      Large randomised       Arrays +   HIGH
HO65 dataset   consortium    trial data incl.       molecular
                             cytogenetic subtype    subtype
                             labels                 data

UAMS GEP data  Multiple      TT2/TT3 trials         Arrays     HIGH
(Barlogie)     GEO           Large n, subtype       GEP70/GEP80
accessions     accessions    annotated, survival    prognostic
                             data                   signatures

CRITICAL NOTES FOR MM ANALYSIS:

NOTE 1 — CD138+ ENRICHMENT IS MANDATORY:
  As with MDS (CD34+ required), ALL MM analysis
  must use CD138+-sorted plasma cells.
  Unsorted bone marrow contains:
    ~80–90% normal cells (lymphocytes, stromal,
    RBCs, granulocytes)
    ~1–60% myeloma cells (depending on stage)
  The myeloma signal in unsorted BM is diluted
  by normal cells.
  CD138+ FACS-sorted cells are the standard
  in MM research — all major datasets use this.
  The universal_discovery_start_script.py must
  be configured to use CD138+ datasets only.

NOTE 2 — MOLECULAR SUBTYPE ANNOTATION REQUIRED:
  The MMRF CoMMpass dataset has FISH data for
  all major translocations (t(11;14), t(4;14),
  t(14;16), del(17p), amp(1q)) — this is the
  ideal dataset for subtype-stratified analysis.
  GEO datasets (GSE167968, GSE223060) have
  variable subtype annotation — check metadata.
  Analysis plan:
  a) Full dataset (all MM vs. normal plasma cells)
     to replicate OrganismCore bulk MM analysis.
  b) Subtype-specific (t(11;14) only, t(4;14) only,
     etc.) to resolve convergence nodes per subtype.

NOTE 3 — THE NORMAL PLASMA CELL REFERENCE:
  Normal plasma cells are RARE in normal bone
  marrow (<0.1% of total BM cells).
  GEO datasets that include normal PCs:
    GSE269875: 5 normal PC samples (spatial)
    GSE223060: normal BM reference included
    MMRF CoMMpass: does NOT include normal PC —
    comparative analysis requires EXTERNAL normal
    PC reference datasets.
  The framework normal reference in MM must be
  normal BONE MARROW PLASMA CELLS, not:
    Not B cells (wrong differentiation stage)
    Not peripheral blood PCs (niche-different)
    Not plasmablasts (incomplete differentiation)
  Using the wrong normal reference will distort
  the depth axis. This is more important in MM
  than in any other cancer in the repository.

NOTE 4 — THE EZH2 DIRECTION TEST BY SUBTYPE:
  EZH2 expression and H3K27me3-associated gene
  signatures must be assessed per molecular subtype.
  Pre-data predictions (from Section VIII):
  t(11;14): EZH2 elevated (predicted)
  t(4;14): EZH2 activity reduced (NSD2 competition)
  t(14;16)/t(14;20): EZH2 relatively elevated (predicted)
  Hyperdiploid: EZH2 near-normal (predicted)
  del(17p): EZH2 variable (predicted)
  These must be confirmed or corrected by data.

NOTE 5 — MYC AS THE DEPTH AXIS PRIMARY DRIVER:
  In the bulk MM depth score, MYC should be
  the dominant depth-positive correlate.
  The depth score is essentially measuring:
  How far has the IRF4-MYC co-activation loop
  escaped from BLIMP1-mediated MYC repression?
  LOW DEPTH: MYC low, BLIMP1 target genes high,
  more plasma cell-like, near normal LLPC.
  HIGH DEPTH: MYC high, proliferation markers
  high, BCL2 or MCL1 high, closer to the
  most aggressive RRMM state.
  The BLIMP1/MYC ratio = the MM depth axis.
  Measurable by standard IHC: Ki67 (proliferation),
  MYC protein, BCL2 protein, CD138 (plasma cell
  marker — should be HIGH and constant, not depth-
  driven). These may constitute the 3-gene MM
  clinical panel.
```

---

## SECTION XII — PLANNED ANALYSIS ORDER

```
ORDER:

  MM-S1    t(11;14)/BCL2     MMRF CoMMpass +      VERY
           (venetoclax       GSE167968             HIGH
           confirmed,        t(11;14) subset
           BCL2 node)        + normal PC ref
                             REASON: The most
                             tractable subtype.
                             BCL2 as convergence
                             node confirmed by
                             BELLINI (venetoclax
                             benefit in t(11;14)).
                             CCND1 vs. BCL2 depth
                             competition: which
                             is the convergence
                             node?
                             EZH2 direction in
                             t(11;14): elevated?
                             CRBN expression depth
                             correlation: does
                             CRBN fall with depth
                             in t(11;14)?
                             (Lenalidomide resistance
                             predictor)
                             t(11;14) + amp(1q)
                             subset: MCL1 elevation
                             = venetoclax resistance
                             signature?
                             Clinical output:
                             BCL2 vs. MCL1 depth
                             stratification within
                             t(11;14) → venetoclax
                             alone vs. venetoclax +
                             MCL1 inhibitor selection.

  MM-S2    Hyperdiploid      MMRF CoMMpass +      VERY
           (most common,     GSE167968             HIGH
           standard risk,    hyperdiploid subset
           CCND2 driver)     + normal PC ref
                             REASON: The largest
                             subgroup (~50% of MM).
                             Bulk MM signal may
                             be dominated by
                             hyperdiploid cells.
                             CCND2 as the primary
                             cell cycle driver
                             (vs. CCND1 in t(11;14)).
                             MCL1 as the primary
                             anti-apoptotic driver.
                             IRF4/MYC loop: is MYC
                             activation equally
                             present in hyperdiploid
                             vs. translocation MM?
                             CKS1B elevation (from
                             trisomy 1q when present)?
                             Clinical output:
                             CDK4/6 inhibitor
                             (palbociclib) vs. MCL1
                             inhibitor selection by
                             depth in hyperdiploid MM.

  MM-S3    t(4;14)/NSD2      MMRF CoMMpass +      HIGH
           (bortezomib-      GSE167968
           sensitive,        t(4;14) subset
           NSD2 node)        + normal PC ref
                             REASON: The clearest
                             epigenetic convergence
                             node in MM (NSD2 vs.
                             EZH2).
                             NSD2 depth correlation:
                             is NSD2 expression
                             linearly depth-correlated?
                             EZH2 direction in t(4;14):
                             is H3K36me2-associated
                             gene activation
                             measurable in bulk
                             RNA-seq (proxy: genes
                             that carry H3K36me2 marks
                             are globally elevated)?
                             HRP2 loss subset:
                             does HRP2 expression
                             anti-correlate with
                             depth in t(4;14)?
                             (HRP2 low + NSD2 high
                             = deepest t(4;14))
                             NSD2 inhibitor target
                             confirmation:
                             does NSD2 depth
                             correlation exceed
                             EZH2 depth correlation?
                             Bortezomib sensitivity
                             mechanism: does the
                             AK2 upregulation
                             appear in t(4;14)
                             depth-positive samples?
                             Clinical output:
                             NSD2 inhibitor patient
                             selection; HRP2 loss
                             as EZH2i + PI
                             re-sensitisation signal.

  MM-S4    MAF subtypes      MMRF CoMMpass +      MOD-HIGH
           t(14;16)/t(14;20) GSE167968
           (high risk,       MAF/MAFB subset
           MAF node)         + normal PC ref
                             REASON: The highest-
                             risk IgH translocation
                             class.
                             MAF/MAFB depth
                             correlation.
                             GC-like gene programme:
                             do MAF-subtype MM
                             cells show BCL6 target
                             gene elevation (partial
                             GC identity retained)?
                             EZH2 in MAF subtypes:
                             elevated (maintaining
                             plasma cell identity
                             against MAF-driven
                             GC reactivation)?
                             CCND2 (shared with
                             hyperdiploid): MAF
                             activates CCND2 just
                             as trisomy does —
                             does CCND2 depth
                             correlation in MAF
                             converge with
                             hyperdiploid CCND2?
                             AKT/PI3K activation
                             depth: is PTEN low
                             or AKT high at depth?
                             Clinical output:
                             CDK4/6 inhibitor
                             (CCND2) vs. AKT
                             inhibitor depth
                             stratification in
                             MAF MM.

  MM-S5    del(17p)/         MMRF CoMMpass +      MOD-HIGH
           double-hit/       del(17p) subset +
           triple-hit        double-hit subset
           (deepest MM       + normal PC ref
           attractor)        REASON: Most adverse
                             subgroup. TP53 circuit
                             absent. Analogous to
                             MDS-biTP53.
                             Is the depth score
                             measurable or saturated
                             in del(17p) MM?
                             MYC activation in
                             del(17p): universally
                             present?
                             BCL2 vs. MCL1: which
                             is dominant in del(17p)?
                             (Predicts venetoclax
                             vs. MCL1i selection)
                             APR-246 applicability:
                             TP53 R273/R175 hotspot
                             mutations present?
                             (APR-246 works on
                             specific mutant p53
                             proteins — not on
                             deletion/null TP53)
                             EZH2 direction: what
                             happens to H3K27me3
                             when TP53 is lost in MM?
                             (TP53 loss + MYC high
                             can reprogram H3K27me3
                             distribution)
                             Clinical output:
                             depth score vs. RRMM
                             therapy sequence;
                             bispecific/CAR-T
                             indication at first
                             relapse in del(17p) MM.

  MM-X     Cross-subtype     After S1–S5.
           + clinical        QUESTIONS:
           integration         1. IRF4 depth
                                  correlation across
                                  all subtypes:
                                  is IRF4 depth-
                                  positive universally?
                                  (Confirming IRF4 as
                                  master convergence node)
                               2. Does the MM depth
                                  score correlate with
                                  the R-ISS or IMWG
                                  clinical risk scores?
                                  (Clinical validation)
                               3. The BLIMP1/MYC ratio:
                                  can it serve as the
                                  3-gene MM clinical panel
                                  approximation together
                                  with Ki67?
                               4. The CANOVA failure
                                  reinterpreted:
                                  in t(11;14) MM, which
                                  depth-BCL2/MCL1
                                  subgroup would have
                                  responded to venetoclax
                                  in CANOVA?
                                  (Retroactive patient
                                  selection analysis)
                               5. MajesTEC-3 (teclistamab
                                  + dara): depth-
                                  independent because
                                  both target retained
                                  plasma cell markers
                                  (BCMA + CD38)?
                                  Confirmation that
                                  immune-targeting of
                                  the retained identity
                                  circumvents depth
                                  stratification.
```

---

## SECTION XIII — WHAT THIS DOCUMENT DOES NOT CONTAIN

```
This document contains:
  ✓ Normal plasma cell differentiation circuit:
    Naive B cell → GC B cell → Plasmablast →
    LLPC, with full TF network (PAX5, BCL6,
    BLIMP1, IRF4, XBP1s, KLF4 missing from
    B cell → MYC repressed in normal LLPC)
  ✓ The AID mechanism creating IgH translocations
    in GC B cells — the founding molecular event
    of translocation-subtype MM
  ✓ The MGUS → SMM → MM → RRMM continuum with
    Waddington depth at each stage
  ✓ MYC activation as the MGUS→SMM/MM transition
    event (structural translocation in ~45% vs.
    native enhancer gain-of-function in ~55%)
  ✓ The augmented attractor concept: MM is not
    identity loss but IDENTITY RETENTION + added
    oncogenic programme (depth = degree of
    oncogenic addition, not identity subtraction)
  ✓ t(11;14): CCND1 + BCL2, venetoclax confirmed
    (BELLINI), lenalidomide reduced sensitivity
    (CRBN low), CCND1 vs. BCL2 convergence node
    question, t(11;14) + amp(1q) = BCL2 + MCL1
    co-elevation requiring dual inhibition
  ✓ t(4;14): NSD2/H3K36me2 mechanism, AK2
    metabolic dependence, HRP2-loss bortezomib
    resistance, NSD2 vs. EZH2 competition —
    NSD2 inhibitor as primary target (not EZH2)
  ✓ MAF subtypes: MAF/MAFB, CCND2, AKT pathway,
    partial GC identity, EZH2 relatively elevated
  ✓ Hyperdiploidy: CCND2 driver, MCL1-dependent
    (not BCL2), venetoclax NOT appropriate,
    lenalidomide-sensitive (CRBN normal)
  ✓ del(17p)/double-hit: TP53 loss, deepest MM
    attractor, analogy to MDS-biTP53, multi-hit
    TP53 worst prognosis confirmed (2025)
  ✓ amp(1q): CKS1B + MCL1, therapy-selected
    clonal evolution, MCL1 inhibitor target
  ✓ EZH2 direction table by subtype:
    t(11;14) elevated, t(4;14) reduced, MAF elevated,
    hyperdiploid intermediate, del(17p) variable
  ✓ IRF4 as master convergence node candidate:
    universal MM, non-oncogene addiction, IRF4-MYC
    feed-forward loop, lenalidomide mechanism via
    IKZF1/3 → IRF4 → MYC collapse
  ✓ CELMoDs (iberdomide, mezigdomide) as next-gen
    lenalidomide equivalents with deeper IKZF1/3
    degradation
  ✓ Bone marrow niche: IL-6/STAT3/MCL1, VLA-4/
    VCAM-1/CAM-DR, BCMA/BAFF/APRIL niche survival,
    RANKL/DKK1 vicious cycle (bone disease)
  ✓ 2024–2025 treatment landscape: Dara-VRd and
    Isa-VRd as quadruplet frontline standard (FDA
    approved 2025), teclistamab + daratumumab
    (MajesTEC-3: paradigm-shifting — 83% 3-year
    PFS in 1-3 prior line RRMM), cilta-cel in
    earlier lines (CARTITUDE-4)
  ✓ Data available: MMRF CoMMpass (dbGaP,
    ~1,000 patients), GSE167968, GSE223060,
    GSE269875 (normal PC reference), UAMS GEP
    data (Barlogie trials)
  ✓ CD138+ enrichment mandatory (same principle
    as CD34+ in MDS)
  ✓ Normal PC reference: must be BM-resident
    plasma cells, not B cells, PBMCs, or
    plasmablasts

This document does NOT contain:
  ✗ Depth score predictions (numerical)
  ✗ Switch gene predictions beyond the structural
    biology analysis presented here
  ✗ Specific drug sequence recommendations
    beyond confirmed clinical data
  ✗ Epigenetic mechanism hypotheses for
    individual patient scenarios
  ✗ Claims about the specific depth score numbers
    the analysis will produce

All of the above belong in the BEFORE documents.
MM-S1a (t(11;14)/BCL2 before-document) is next.
Written before any script runs.
Before any data loads.
```

---

## STATUS BLOCK

```
document:           MM_Subtype_Orientation.md
folder:             Cancer_Research/MM/Subtypes/
status:             COMPLETE — ORIENTATION ONLY
date:               2026-03-04
author:             Eric Robert Lawson
                    OrganismCore

subtypes_covered:
  t(11;14)/CCND1/BCL2 (~20%):   BCL2 convergence   [1 of 6]
                                 node (venetoclax
                                 confirmed), CRBN
                                 low → lenalidomide
                                 less effective

  t(4;14)/NSD2/FGFR3 (~13%):    NSD2/H3K36me2      [2 of 6]
                                 node, AK2 metabolic
                                 dependency,
                                 bortezomib sensitive,
                                 EZH2 activity reduced

  t(14;16)/MAF; t(14;20)/MAFB   Partial GC          [3 of 6]
  (~7% combined):                identity, CCND2,
                                 AKT pathway, EZH2
                                 relatively elevated

  Hyperdiploid (~50%):           CCND2, MCL1-        [4 of 6]
                                 dependent (not BCL2),
                                 lenalidomide-sensitive,
                                 standard risk

  del(17p)/TP53 (~10%):          Deepest MM          [5 of 6]
                                 attractor, analogy
                                 to MDS-biTP53,
                                 multi-hit TP53 worst

  amp(1q) (~40%, overlapping):   CKS1B + MCL1,       [6 of 6]
                                 therapy-selected,
                                 MCL1 inhibitor target

key_structural_difference_from_all_other_cancers:
  MM IS AN AUGMENTED ATTRACTOR, NOT A REPLACEMENT.
  The myeloma cell RETAINS full plasma cell identity
  (IRF4, BLIMP1, XBP1s, BCMA, CD138, M-protein).
  It ADDS a survival/proliferation programme on top.
  Depth = degree of oncogenic addition, not
  identity subtraction.
  This inverts the depth axis relative to all
  solid tumours (PRAD, PAAD, STAD, BRCA, GBM).
  The normal reference must be BONE MARROW PLASMA
  CELLS, not distant tissue or GTEx.

key_convergence_node_candidate:
  IRF4 — non-oncogene addiction, universal across
  all MM subtypes, drives the IRF4-MYC feed-forward
  loop, confirmed as lethal dependency by lenalidomide
  mechanism (IKZF1/3 → IRF4 → MYC → cell death).
  SUBTYPE-SPECIFIC CONVERGENCE NODES:
  t(11;14): BCL2 (or CCND1) — venetoclax confirmed
  t(4;14): NSD2 — NSD2 inhibitors in trials
  MAF: MAF (undruggable directly; CDK4/6 and AKT
       as indirect nodes)
  Hyperdiploid: MCL1 — MCL1 inhibitors in trials
  del(17p): none tractable (circuit broken)

key_pattern_4_application:
  EZH2 direction is SUBTYPE-SPECIFIC in MM:
  t(11;14) → EZH2 relatively elevated (predicted)
  t(4;14) → EZH2 activity reduced (NSD2 competes)
  MAF → EZH2 relatively elevated (predicted)
  Hyperdiploid → EZH2 near-normal (predicted)
  del(17p) → EZH2 variable (predicted)
  DO NOT assume EZH2 inhibitor is appropriate
  in t(4;14) — EZH2 is already outcompeted by NSD2.
  EZH2 inhibitor in t(4;14) without HRP2-loss
  selection = potentially harmful.

key_2025_clinical_landscape_finding:
  MajesTEC-3 (teclistamab + daratumumab, relapsed
  MM): paradigm-shifting. 83.4% 3-year PFS vs.
  29.7% standard of care. BCMA + CD38 co-targeting.
  Both targets are RETAINED PLASMA CELL IDENTITY
  markers — depth-independent.
  This is the clinical confirmation that immune
  targeting of the RETAINED IDENTITY circumvents
  the depth stratification problem:
  You don't need to know the subtype or depth to
  select teclistamab + dara — BCMA and CD38 are
  expressed in ALL MM subtypes at ALL depths.
  Depth stratification becomes important when
  the TARGET IS THE ONCOGENIC ADD-ON:
  venetoclax in t(11;14) (targets BCL2 add-on),
  bortezomib-based PIs in t(4;14) (targets NSD2-
  driven AK2/UPR dependency add-on), lenalidomide
  (targets IRF4-MYC loop add-on in non-t(11;14)).
  The depth score framework's contribution to MM
  is patient selection for the ONCOGENIC ADD-ON
  targeting therapies — not for the immunotherapy
  (BCMA/CD38) class, which works across all depths.

analyses_started:   0 (new subtype series)
existing_analysis:  Cancer_Research/MM/ — complete
                    (OrganismCore cancer series,
                    bulk signal)
next_document:      MM-S1a
                    t(11;14)/BCL2 before-document
```
