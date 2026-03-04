# PANCREATIC DUCTAL ADENOCARCINOMA (PAAD/PDAC) — SUBTYPE ORIENTATION DOCUMENT
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

  A complete map of the pancreatic ductal
  adenocarcinoma molecular subtype landscape —
  what the two robust epithelial subtypes are
  (and what the two contested subtypes were),
  where each arises in the normal exocrine
  pancreatic hierarchy, what the initiating
  and cooperating genetic events are, what the
  clinical characteristics and treatment are,
  and what public data exists to analyze each.

PAAD holds a specific position in the OrganismCore
repository because it was the cancer that first
demonstrated intact circuit restoration:

  THE CIRCUIT INTEGRITY FINDING.

  The existing PAAD bulk analysis (OrganismCore
  cancer series) found that the master
  differentiation TF circuit IS INTACT in PAAD —
  unlike STAD where the circuit is broken.

  Pattern 5 from the OrganismCore Cancer Framework:
    In PAAD (and PRAD):
      The master differentiation TF circuit is intact.
      Restoring the switch gene would execute the
      differentiation programme.
      Circuit restoration is therapeutic.
    In STAD (and some other cancers):
      The circuit is broken.
      Restoring the switch gene would NOT execute
      the programme.
      Attractor dissolution is the correct strategy.

  The intact circuit finding in PAAD means:
    The switch gene (PTF1A — the master acinar
    TF identified in the existing analysis)
    sits at the top of a downstream programme
    that is still connected.
    If PTF1A expression is restored in a PAAD
    cell — the cell would re-enter the acinar
    differentiation programme and exit the
    cancer state.
    This is circuit restoration as therapy,
    not attractor dissolution.

  THIS STRUCTURAL FINDING — the intact circuit —
  is the most therapeutically important result
  in the PAAD series and is the primary context
  for all subtype analysis that follows.

  The subtype orientation document situates this
  finding: WHICH subtype has an intact circuit?
  Is circuit integrity subtype-specific?
  Does the basal-like subtype have a MORE broken
  circuit than classical? These questions motivate
  the subtype series.
```

---

## DOCUMENT METADATA

```
document_id:        PAAD_Subtype_Orientation
series:             PAAD (Pancreatic Ductal Adenocarcinoma
                    — Subtypes)
folder:             Cancer_Research/PAAD/Subtypes/
date:               2026-03-03
author:             Eric Robert Lawson
status:             ORIENTATION ONLY — no predictions,
                    no analysis, no epigenetic hypotheses
next_document:      PAAD_Classical_before.md
                    (Document PAAD-S1a — Classical
                    subtype before-doc)
protocol_version:   Workflow_Protocol.md v2.0
existing_analysis:  See OrganismCore cancer series
                    for the completed PAAD bulk analysis
                    (PTF1A switch gene identified,
                    intact circuit finding, EZH2
                    elevated — gain-of-function lock,
                    circuit restoration as therapeutic
                    strategy).
```

---

## SECTION I — THE NORMAL EXOCRINE PANCREAS

```
The Waddington baseline for pancreatic cancer
is the normal exocrine pancreatic epithelium —
the acinar-ductal unit that constitutes ~95%
of pancreatic mass and is the cell of origin
for PDAC.

THE PANCREAS IS TWO ORGANS IN ONE.
The pancreas contains two completely distinct
functional compartments:

  ENDOCRINE (~2% of mass):
    The islets of Langerhans.
    Beta cells (insulin), alpha cells (glucagon),
    delta cells (somatostatin), PP cells.
    These cells arise from a separate lineage
    (Neurogenin3+ endocrine progenitors) and are
    the cell of origin for pancreatic
    neuroendocrine tumours (pNET) — NOT PDAC.
    The endocrine compartment is IRRELEVANT to
    PDAC biology directly but IS the source of
    the ADEX artefact (see Section V).

  EXOCRINE (~95–98% of mass):
    ACINAR CELLS + DUCTAL CELLS.
    The cell of origin for PDAC.
    This is the organ of interest.

═══════════════════════════════════════════════════════
THE EXOCRINE PANCREATIC HIERARCHY
═══════════════════════════════════════════════════════

ACINAR CELLS (the terminal differentiated cells)
  Location:   The acini — grape-like clusters at the
              terminal units of the ductal tree.
              Each acinus is surrounded by a basement
              membrane and drains into a centroacinar
              ductule.
  Function:   Enzyme secretion.
              Acinar cells are among the most
              secretory cells in the body — they
              synthesise and package:
                Trypsinogen (PRSS1, PRSS2)
                Chymotrypsinogen (CTRB1, CTRB2)
                Elastase (CELA1, CELA2A, CELA3A)
                Carboxypeptidase A (CPA1, CPA2)
                Lipase (PNLIP, CELA3B)
                Amylase (AMY2A, AMY2B)
                DNase (DNASE1L3)
              These are synthesised as inactive
              zymogens and activated by enterokinase
              in the duodenum.
              The acinar cell's biology is entirely
              devoted to this enormous secretory load:
                Rough ER: massive (to make protein)
                Golgi: massive (to package)
                Zymogen granules: prominent
                Mitochondria: dense
              The acinar cell is the chief cell
              of the pancreas — exactly analogous to
              the gastric chief cell (MIST1 is shared
              between them: MIST1 is the TF for the
              secretory programme in BOTH gastric
              chief cells AND pancreatic acinar cells).
  Markers:    PTF1A — THE MASTER TF OF ACINAR IDENTITY
              NR5A2 (LRH-1) — acinar identity cofactor
              RBPJL — PTF1A complex partner
              MIST1 (BHLHA15) — secretory programme TF
              CPA1 — enzyme marker
              PRSS1/PRSS2 — trypsinogen
              AMY2A — amylase
              CELA family — elastase
  Self-renewal: VERY LOW under normal conditions.
                Acinar cells are long-lived, post-
                mitotic-like cells under homeostasis.
                Regenerate in response to injury
                (pancreatitis) — the stem-like
                population in the adult pancreas
                remains controversial but may include
                both centroacinar cells and acinar cells
                themselves via dedifferentiation.

  PTF1A — THE GUARDIAN TF OF ACINAR IDENTITY:
    PTF1A (Pancreas Transcription Factor 1a) is
    the master transcription factor for:
      a) Pancreatic specification (embryo): PTF1A
         specifies the entire pancreatic organ from
         the foregut endoderm. PTF1A knockout mice
         have no pancreas.
      b) Acinar cell maintenance (adult): PTF1A
         maintains the acinar gene expression
         programme in the adult gland.
         It forms a trimeric complex:
         PTF1A + RBPJL + E-protein
         This complex activates:
           All digestive enzyme genes
           The secretory programme
           The acinar identity programme
         And REPRESSES:
           Ductal identity genes
           Progenitor genes
           MYC targets (proliferation)
    PTF1A LOSS → immediate acinar identity
    destabilisation → ADM (acinar-to-ductal
    metaplasia) — the first step in PDAC
    formation.
    PTF1A is the CDKN2A analogue of identity —
    it is both a tumour suppressor (when lost,
    enables PDAC) AND a differentiation guardian
    (when active, enforces acinar fate).
    PTF1A IS THE SWITCH GENE IDENTIFIED IN THE
    EXISTING PAAD ANALYSIS.
    Its restoration = circuit restoration therapy
    (the intact circuit finding).

  NR5A2 (LRH-1) — THE IDENTITY COFACTOR:
    Nuclear receptor subfamily 5, group A.
    NR5A2 cooperates with PTF1A to maintain
    acinar identity.
    NR5A2 haploinsufficiency in humans is a
    GWAS risk factor for PDAC — supporting its
    tumour suppressor role.
    NR5A2 activates: PTF1A targets, CPA1,
    PRSS1, AMY2A.
    NR5A2 loss: acinar identity destabilises,
    KRAS-driven ADM is accelerated.
    NR5A2 is the second node in the acinar
    identity circuit alongside PTF1A.

CENTROACINAR CELLS (the isthmus equivalent)
  Location:   The cells that line the lumen where
              the acinus joins the ductule.
              A small population at the junction
              of acini and ducts.
  Function:   May represent a stem/progenitor
              compartment for the exocrine pancreas.
              The controversy: whether PDAC arises
              from acinar cells (via ADM), ductal
              cells (directly), or centroacinar cells
              (as a progenitor source).
              Current evidence favours ACINAR CELLS
              as the primary cell of origin for PDAC
              (via ADM and PanIN formation with KRAS
              activation), though ductal-origin models
              are not excluded.
  Markers:    SOX9 (ductal/progenitor TF)
              HES1 (Notch target)
              CA9 (carbonic anhydrase)
              Low PTF1A, low enzyme expression

DUCTAL CELLS (the conduit system)
  Location:   Line the ductal tree from main pancreatic
              duct to ductules.
              The main duct, secondary branches,
              tertiary branches, ductules.
  Function:   Transport of enzyme-containing fluid
              from acini to the duodenum.
              Secrete bicarbonate (CFTR — the
              cystic fibrosis protein lives here).
              HCO3- secretion alkalinises the
              ductal fluid to prevent premature
              enzyme activation.
  Markers:    KRT7, KRT8, KRT18, KRT19 —
              simple ductal cytokeratins
              SOX9 — THE DUCTAL TF
                     (repressed by PTF1A in acinar cells)
              HNF1B — ductal TF
              FOXA1, FOXA2 — ductal TF
              CFTR — bicarbonate secretion
              MUC1, MUC6 — ductal mucins
              CA2 — carbonic anhydrase
  Self-renewal: MODERATE — ductal cells turn over
                more readily than acinar cells.
                Main duct epithelium has an active
                progenitor compartment.

THE WADDINGTON STRUCTURE OF NORMAL EXOCRINE
PANCREAS:

  Embryonic progenitor (PTF1A+, SOX9+, NKX6.1+)
    ↓ ACINAR COMMITMENT (PTF1A high, SOX9 falls)
    ↓ Acinar differentiation (PTF1A/RBPJL/NR5A2
      complex assembled)
  ACINAR CELL (PTF1A+, CPA1+, PRSS1+, AMY2A+)
    = THE TERMINAL NORMAL ATTRACTOR

  Embryonic progenitor
    ↓ DUCTAL COMMITMENT (SOX9 maintained, PTF1A off)
  DUCTAL CELL (SOX9+, KRT7+, KRT19+, CFTR+)
    = THE SECONDARY NORMAL ATTRACTOR

  CANCER FORMATION ROUTE (the false attractor):
  Normal acinar cell (PTF1A+, KRAS wild-type)
    ↓ KRAS G12D/G12V mutation (the founding event)
    ↓ Inflammation / injury / pancreatitis
    ↓ PTF1A loss → ACINAR-TO-DUCTAL METAPLASIA (ADM)
      (the acinar cell loses its identity and
      transdifferentiates toward a ductal-progenitor
      state: PTF1A falls, SOX9 rises, CPA1 falls,
      KRT19 rises)
  ADM CELL (SOX9+, KRT19+, low PTF1A, KRAS G12D)
    = THE PRENEOPLASTIC INTERMEDIATE (Waddington
      transition state on the way to PanIN)
    ↓ Additional mutations: CDKN2A loss (PanIN-2)
    ↓ TP53 mutation, SMAD4 loss (PanIN-3 → invasive)
  CLASSICAL PDAC FALSE ATTRACTOR:
    (PTF1A gone, GATA6 present, HNF1A/HNF4A
    present, KRT19 high, ductal-like)
  OR → further identity loss →
  BASAL-LIKE/SQUAMOUS PDAC FALSE ATTRACTOR:
    (GATA6 gone, KRT5/KRT14/TP63 expressed,
    most primitive false attractor)

  THE DEPTH AXIS IN PAAD:
  Normal acinar cell (deepest identity, PTF1A+)
  → Classical PDAC (partial ductal identity,
    GATA6+, PTF1A lost)
  → Basal-like PDAC (total identity loss, GATA6−,
    TP63+, KRT5+)

  NOTE: The depth axis in PAAD is MORE COMPLEX
  than in most cancers because the normal tissue
  itself has two lineages (acinar and ductal).
  The classical subtype resembles DUCTAL cells
  (not acinar cells) — it has adopted a ductal
  false attractor identity. The acinar identity
  is what was LOST.
  The depth axis therefore measures:
    Distance from the normal acinar identity
    (how much PTF1A, CPA1, AMY2A, PRSS1 is left)
  AND separately:
    Distance from any epithelial identity at all
    (how much GATA6, HNF1A, KRT8/18/19 is left
    in the classical subtype, none of which are
    left in the basal-like subtype).
```

---

## SECTION II — THE TWO ROBUST SUBTYPES

```
Multiple classification systems have been proposed
for PDAC. Three groups deserve acknowledgment:

  Collisson et al. (2011): Classical, Quasi-
  Mesenchymal (QM), Exocrine-like.
  Three subtypes.

  Moffitt et al. (2015): Tumor Classical, Tumor
  Basal-like; Stroma Normal, Stroma Activated.
  Two tumor subtypes + two stromal subtypes.
  The cleanest system — separates tumour from
  stroma explicitly.

  Bailey et al. (2016): Pancreatic Progenitor,
  Squamous, Immunogenic, ADEX.
  Four subtypes — but two are artefacts.

  TCGA (2017): Confirmed primarily two robust
  epithelial subtypes.

THE 2024 CONSENSUS:
  Two robust tumor-intrinsic subtypes:
    CLASSICAL (= Pancreatic Progenitor,
                Moffitt Classical)
    BASAL-LIKE (= Squamous, Moffitt Basal-like,
                  Collisson Quasi-Mesenchymal)

  Two contested subtypes with artefact evidence:
    IMMUNOGENIC: Largely an immune cell
                 infiltration signal, not a
                 tumor-intrinsic programme.
    ADEX:        Largely an acinar/endocrine cell
                 contamination signal — NOT a
                 tumor cell programme.

  These are documented in Section V.

SINGLE MARKER THAT SEPARATES THEM:
  GATA6 expression by IHC or RNA.
  GATA6 HIGH → Classical (better prognosis)
  GATA6 LOW/ABSENT → Basal-like (worse prognosis)
  GATA6 is the single most clinically useful
  molecular marker in PDAC — its IHC status
  can be determined on the diagnostic biopsy.
  It is the GATA3 of PDAC — a single TF that
  marks subtype identity and predicts treatment
  response.
```

---

## SECTION III — CLASSICAL SUBTYPE: THE DIFFERENTIATED FALSE ATTRACTOR

```
FREQUENCY:     ~60% of PDAC
               (estimates range 40–65% depending
               on purity filtering and classification
               system — the majority subtype)

PROGNOSIS:     BETTER (relative to basal-like)
               Resectable: median OS 18–24 months
               Metastatic: median OS ~11–14 months
               (vs. ~6–9 months for basal-like)
               The prognosis is still extremely poor
               by any absolute standard — "better"
               in PDAC is still devastating.

CELL OF ORIGIN:
  Acinar cell that has undergone ADM — the
  acinar cell LOST its identity (PTF1A gone)
  and adopted a ductal-progenitor false attractor
  (GATA6+, HNF1A+, HNF4A+, SOX9 retained in
  some cells).
  The classical PDAC cell resembles a ductal
  cell more than an acinar cell — it has
  traveled from the terminal acinar attractor
  to a ductal-identity false attractor.
  This is the WADDINGTON TRANSITION from acinar
  to ductal-like identity, STABILISED by KRAS
  + loss of CDKN2A + TP53 mutation.
  The cells form glands (tubular/ductal
  architecture on histology) — consistent with
  ductal identity being maintained even in the
  false attractor.
  GATA6 is NOT a normal acinar cell marker —
  it is a DUCTAL IDENTITY MARKER that the
  classical PDAC cell has adopted as part of
  its ductal-false-attractor programme.
  GATA6 in PDAC = the residual ductal identity
  maintained in the false attractor.
  Classical PDAC is a DUCTAL IDENTITY FALSE
  ATTRACTOR — the most differentiated false
  attractor in PDAC, just as luminal papillary
  is in BLCA.

DEFINING MOLECULAR EVENTS:
  KRAS mutation:      ~95% — universal, founding
                      event, G12D (most common in
                      PDAC), G12V (second), G12R
                      KRAS G12D: the dominant allele
                      in classical PDAC
                      KRAS is the originating oncogene
                      — present in PanIN-1, the
                      earliest identifiable lesion,
                      before any other driver mutation
  CDKN2A loss:        ~95% — universal (homozygous
                      deletion ~40%, promoter
                      methylation ~35%, mutation ~20%)
                      p16 loss → CDK4/6 constitutively
                      active → Rb phosphorylation →
                      unrestrained G1/S progression
                      This is the cell cycle brake
                      that allows the KRAS-driven ADM
                      cell to proliferate without
                      restraint
  TP53 mutation:      ~75%
                      Late-stage driver (PanIN-3
                      → invasive transition)
                      Gain-of-function p53 alleles
                      common in PDAC (not just null)
  SMAD4 deletion/     ~50%
  mutation:           TGF-β pathway tumour suppressor
                      SMAD4 loss:
                        Allows TGF-β signalling to
                        switch from tumour suppressive
                        to tumour promoting (TGF-β
                        promotes EMT and invasion in
                        SMAD4-null cells)
                        Correlates with haematogenous
                        metastasis (liver mets in
                        SMAD4-null; peritoneal mets
                        in SMAD4-intact)
                      SMAD4 status as a clinical
                      stratification tool:
                        SMAD4-null: liver dominant
                        disease → adjuvant chemo
                        most critical
                        SMAD4-intact: locally
                        aggressive / peritoneal →
                        radiation consideration
  GATA6 expression:   HIGH — defines classical subtype
                      GATA6 is not mutated — it is
                      RETAINED/EXPRESSED in classical
                      PDAC as the residual ductal
                      identity marker.
                      GATA6 in classical PDAC drives:
                        Pancreatic progenitor genes
                        (HNF1A, HNF4A, FOXA2)
                        Epithelial adhesion
                        Differentiation markers
                      GATA6 is a DEPTH-NEGATIVE gene
                      in PAAD: falls as depth increases
                      (classical → basal-like
                      transition = GATA6 loss)
  EZH2:               ELEVATED (the existing PAAD
                      analysis finding — gain-of-
                      function epigenetic lock)
                      EZH2 drives H3K27me3 at the
                      PTF1A locus and other acinar
                      identity gene loci — epigenetically
                      silencing the return to acinar
                      identity.
                      EZH2 is the EPIGENETIC LOCK
                      that maintains the ductal false
                      attractor by preventing PTF1A
                      re-expression.
                      This is analogous to EZH2 in
                      BRCA (where EZH2 drives the
                      cancer false attractor by
                      silencing luminal differentiation
                      genes).
                      EZH2 is ELEVATED in classical
                      PDAC and is the primary
                      convergence node for circuit
                      restoration therapy.
                      TAZEMETOSTAT (EZH2 inhibitor)
                      — the predicted drug target from
                      the existing PAAD analysis —
                      is the circuit restoration agent:
                      Block EZH2 → H3K27me3 cleared
                      at PTF1A locus → PTF1A expression
                      can return → acinar differentiation
                      programme re-executed → exit from
                      false attractor.

CIRCUIT INTEGRITY — WHY PAAD IS DIFFERENT:
  In the existing analysis (Framework Pattern 5):
    PAAD has an intact circuit.
    The PTF1A → downstream acinar programme
    connections are NOT severed in classical PDAC.
    The upstream blocker (EZH2) is the ONLY
    reason the programme is not executing.
    Remove the blocker → programme executes.
  This makes PAAD structurally analogous to:
    Prostate cancer (PRAD) — also intact circuit,
    also EZH2-elevated, also amenable to circuit
    restoration.
  The contrast with STAD:
    STAD GS subtype — CDH1 lost, the downstream
    connections are severed. Restoring a single
    switch gene cannot re-execute the programme.
  The clinical implication:
    Tazemetostat in classical PDAC is not just
    a cytostatic agent — it is potentially a
    DIFFERENTIATION THERAPY that forces the
    cancer cell back to acinar identity.
    This is a fundamentally different mechanism
    from gemcitabine (killing by DNA damage)
    or FOLFIRINOX (multi-target cytotoxic).
    It addresses the FALSE ATTRACTOR MECHANISM,
    not just downstream proliferation.

TREATMENT:
  Resectable disease:
    Whipple procedure (pancreaticoduodenectomy)
    for head of pancreas.
    Distal pancreatectomy for body/tail.
    Adjuvant chemotherapy post-resection:
      Modified FOLFIRINOX — standard of care
      for fit patients (PRODIGE 24 trial)
      gemcitabine/capecitabine — alternative
      for less-fit patients (ESPAC-4 trial)
    Neoadjuvant approach (borderline resectable):
      FOLFIRINOX or gemcitabine/nab-paclitaxel
      × 4–6 months → surgery → adjuvant
  Metastatic:
    FOLFIRINOX (5-FU/leucovorin + oxaliplatin +
    irinotecan): ~31% ORR, median OS 11.1 months
    OR:
    Gemcitabine + nab-paclitaxel: ~23% ORR,
    median OS 8.5 months
    NALIRIFOX (liposomal irinotecan + 5-FU +
    oxaliplatin): NAPOLI-3 trial — superior to
    gem/nab-paclitaxel (OS 11.1 vs 9.2 months)
    Approved as first-line option (2023–2024).
    GATA6-HIGH (classical) patients respond
    BETTER to intensive chemotherapy (FOLFIRINOX,
    NALIRIFOX) than GATA6-LOW (basal-like).
    GATA6 IHC is increasingly used to guide
    chemotherapy selection in clinical practice.
    BRCA1/2 or PALB2 germline mutation (~7–9%
    of PDAC): olaparib (PARP inhibitor) as
    maintenance after platinum-based first-line
    (POLO trial — standard of care).
  Framework-derived novel prediction:
    Tazemetostat + FOLFIRINOX: the circuit
    restoration + cytotoxic combination
    targeting classical PDAC.
    Rationale:
      Tazemetostat removes EZH2 lock → PTF1A
      can return → acinar differentiation begins
      FOLFIRINOX kills cells in the transition
      (cells exiting false attractor are
      transiently vulnerable)
      Combination exploits both the attractor
      geometry AND the cytotoxic vulnerability
    Status: Not yet in clinical trials as
    a combination for classical PDAC.
    The EZH2 inhibitor in PDAC concept is
    in early trials (PDAC is not an approved
    EZH2i indication as of 2026).
    This is the framework's novel prediction
    for classical PDAC.
```

---

## SECTION IV — BASAL-LIKE / SQUAMOUS SUBTYPE: THE IDENTITY-LOSS ATTRACTOR

```
FREQUENCY:     ~30–40% of PDAC
               (range depends on classification
               system and tumour purity filtering)
               The most aggressive PDAC subtype.

PROGNOSIS:     WORST
               Resectable: median OS 10–14 months
               Metastatic: median OS 6–9 months
               The subtype-specific prognosis gap
               in PDAC is larger than in most
               solid tumours — basal-like patients
               have roughly HALF the survival of
               classical patients.

CELL OF ORIGIN:
  Two hypotheses not fully resolved:
  a) Acinar cell that underwent MORE EXTREME
     ADM — lost not just PTF1A (as in classical)
     but also ALL residual ductal-epithelial
     identity (GATA6, HNF1A, HNF4A) — and instead
     accessed a squamous/basal progenitor programme
     (TP63, KRT5, KRT14) that is not normally
     expressed in the pancreas at all.
  b) A distinct cell-of-origin in the ductal
     compartment or a rare basal progenitor
     that directly accessed the squamous programme.
  Current evidence favours (a): the basal-like
  programme in PDAC is an ABERRANT activation
  of an embryonic squamous/basal programme in
  response to extreme oncogenic stress —
  KRAS amplification, TP53 mutation, SMAD4 loss
  all occurring, destabilising epithelial identity
  so completely that the cell accesses a non-
  pancreatic identity programme.
  THIS IS THE MOST DRAMATIC FALSE ATTRACTOR
  IN THE PANCREATIC LANDSCAPE:
  The cell originated in the pancreas but has
  acquired gene expression from a squamous
  epithelial programme it was never supposed
  to express.
  Cross-repository connection:
  Basal-like PDAC is structurally analogous to:
    Basal/squamous BLCA (TP63+, KRT5+)
    Claudin-low breast cancer (EMT, basal)
    CMS4 CRC (worst prognosis, mesenchymal)
    GS gastric cancer (identity loss)
  Different tissues. Same false attractor type.
  The most identity-lost false attractor in
  each tissue has poor prognosis, immune
  exclusion, and chemotherapy resistance.
  This is Pattern 2 (Drug Target Derivation)
  applied across the repository.

DEFINING MOLECULAR EVENTS:
  KRAS mutation:       ~95% — same as classical
                       BUT: more frequent KRAS
                       amplification (KRAS copy
                       number gain) and higher
                       KRAS allele frequency
  TP53 mutation:       ~85% — higher than classical
                       (almost universal in basal-like)
  SMAD4 loss:          ~60% — higher than classical
  CDKN2A loss:         ~95% — same as classical
  GATA6:               LOST — defines basal-like
                       GATA6 loss is not mutation —
                       it is epigenetic silencing
                       (promoter methylation or
                       PRC2-mediated H3K27me3)
                       GATA6 loss → the ductal
                       identity programme collapses
                       → the cell has no anchor to
                       any differentiated pancreatic
                       programme
  TP63 (ΔNp63):        EXPRESSED — defines basal-like
                       The same TP63 that appears in
                       basal/squamous BLCA and basal
                       breast cancer.
                       In the pancreas, TP63 is NOT
                       normally expressed.
                       Its aberrant activation in
                       basal-like PDAC is the
                       acquisition of a foreign
                       programme — squamous basal
                       identity from another tissue.
  KRT5, KRT14:         EXPRESSED — squamous keratins
                       Not normally expressed in
                       exocrine pancreas.
                       These are the markers of the
                       squamous/basal false attractor.
  MYC amplification:   ~50% of basal-like PDAC
                       Higher than classical.
                       MYC drives rapid proliferation
                       and is the downstream
                       proliferative engine of the
                       basal-like false attractor.
  ARID1A mutation:     Higher in basal-like than
                       classical — chromatin
                       remodelling gene loss
                       contributing to the epigenetic
                       instability that allows
                       foreign programme (squamous)
                       access.

CIRCUIT INTEGRITY — THE CRITICAL QUESTION:
  The existing PAAD analysis found intact circuit
  in the BULK PAAD signal (which is ~60%
  classical-dominated).
  THE SUBTYPE QUESTION:
  Is the circuit intact in BASAL-LIKE PDAC?
  The structural hypothesis (to be tested in
  PAAD-S2a):
    Basal-like PDAC has a MORE BROKEN circuit
    than classical PDAC.
    In classical PDAC: PTF1A is silenced by
    EZH2 but the downstream acinar programme
    (CPA1, PRSS1, AMY2A, etc.) is still
    connected to PTF1A. If EZH2 is blocked,
    PTF1A can re-engage the programme.
    In basal-like PDAC: PTF1A is gone, GATA6
    is gone — two levels of identity are lost.
    Additionally, TP63 is aberrantly active,
    actively driving a competing squamous
    programme. Restoring PTF1A in this context
    would not simply re-execute the acinar
    programme — the basal-like false attractor
    is ACTIVELY maintained by TP63, not just
    passively maintained by EZH2 silencing.
    The therapeutic implication:
      CLASSICAL PDAC: EZH2 inhibition → PTF1A
      restoration → circuit re-execution →
      differentiation therapy (CIRCUIT
      RESTORATION).
      BASAL-LIKE PDAC: EZH2 inhibition alone
      may be INSUFFICIENT. TP63 knockdown +
      EZH2 inhibition may both be required.
      OR: ATTRACTOR DISSOLUTION (like STAD GS)
      — dissolve the false attractor rather
      than restore the differentiation circuit.
      TP63 ΔNp63 inhibition is a potential
      basal-like PDAC target.
  This structural distinction is the most
  important novel prediction in the PAAD
  subtype series.

TREATMENT:
  Chemotherapy response: POOR relative to
  classical PDAC.
  FOLFIRINOX: less benefit in basal-like
  (GATA6-low patients do worse on all regimens)
  PASS-01 trial data: GATA6-high (classical)
  responded better to gemcitabine/nab-paclitaxel
  (OS 13.9 months) than GATA6-low (basal-like).
  PANGEA trial (active 2025):
    EGFR inhibitors (erlotinib) specifically
    for basal-like PDAC (PurIST basal).
    Rationale: EGFR signalling through KRAS
    is more relevant in basal-like context where
    alternative survival signals dominate.
    Gemcitabine + erlotinib + nab-paclitaxel
    vs standard — specifically in basal-like
    subtype patients (NCT06483555).
    This is the first prospectively subtype-
    stratified trial in PDAC — the clinical
    implementation of the framework's logic.
  KRAS inhibitors:
    KRAS G12D is most common in PDAC overall
    and enriched in basal-like.
    MRTX1133 (KRAS G12D-specific covalent
    inhibitor): preclinical data shows activity
    in PDAC — Phase I/II entering 2024–2025.
    RMC-9805 (KRAS G12D ON/OFF state inhibitor):
    Phase I active.
    Resistance mechanism paradox:
    Basal-like cells are INITIALLY more sensitive
    to KRAS inhibitors (because the basal-like
    programme depends heavily on KRAS signalling)
    but quickly ESCAPE by adopting a classical-
    like state (the classical state is relatively
    KRAS-independent). This is the most striking
    finding in the recent PDAC KRAS inhibitor
    literature (GSE271300 dataset): resistance
    to KRAS inhibition involves SUBTYPE SWITCHING
    — basal-like cells escape by reverting to
    classical identity.
    This is attractor landscape dynamics in real
    time: KRAS inhibition destabilises the
    basal-like false attractor → cells slide
    to the adjacent classical false attractor →
    survival despite KRAS inhibition.
    The therapeutic implication:
    KRAS inhibitor + EZH2 inhibitor (prevents
    classical-attractor stabilisation as an
    escape route) — a rational combination
    targeting both the primary false attractor
    and the escape route simultaneously.
  AURKA inhibitors:
    AURKA is a depth-positive gene in PDAC bulk
    analysis (cross-reference from STAD where
    AURKA/ZEB2 co-expression was found).
    Basal-like PDAC may show higher AURKA
    expression than classical.
    Alisertib (Aurora A inhibitor) — the drug
    predicted in the bulk PAAD/STAD analysis —
    may be MORE relevant for basal-like PDAC
    than classical.
    This is a structural hypothesis to test
    in PAAD-S2a.
```

---

## SECTION V — THE CONTESTED SUBTYPES: ADEX AND IMMUNOGENIC

```
ADEX (Aberrantly Differentiated Endocrine-Exocrine):
  Original description (Bailey 2016):
    A subtype characterised by high expression
    of exocrine (acinar enzyme) genes and
    endocrine (islet) genes.
    "Aberrantly differentiated" — suggesting
    the cancer cells had re-activated acinar
    and endocrine programmes.
  2024 consensus re-evaluation:
    ADEX is NOT a tumour-intrinsic subtype.
    The ADEX signal arises from CONTAMINATION
    of bulk RNA-seq samples by normal/non-
    malignant exocrine acinar cells and
    endocrine islet cells that are admixed
    with the tumour.
    PDAC tumours contain ~60–80% non-tumour
    cells by volume (desmoplastic stroma,
    acinar remnants, islets, immune cells).
    Bulk RNA-seq from an ADEX-classified tumour
    is measuring:
      ~30% tumour (classical or basal-like)
      ~30% acinar cell contamination
      ~20% endocrine contamination
      ~20% stroma/immune
    The "acinar and endocrine genes" are from
    the non-tumour cells — NOT from the cancer.
    When tumour purity is controlled for (using
    ESTIMATE, TIMER, or single-cell deconvolution),
    the ADEX signal disappears — these samples
    resolve to classical PDAC.
  FRAMEWORK IMPLICATION:
    ADEX was the single most seductive false
    positive in PDAC subtyping — it appeared to
    show cancer cells with PTF1A activity
    (the switch gene!) but was actually just
    measuring neighbouring normal acinar cells
    whose PTF1A programme was intact.
    This is a critical methodological point for
    the depth score analysis:
    THE DEPTH-NEGATIVE GENES IN PAAD (PTF1A,
    CPA1, PRSS1, AMY2A) MUST BE INTERPRETED WITH
    TUMOUR PURITY IN MIND.
    A tumour with high PTF1A expression in bulk
    RNA-seq is NOT necessarily a cancer cell with
    high PTF1A — it may be measuring acinar cell
    contamination.
    TCGA-PAAD tumour purity estimates must be
    examined alongside depth score derivation.
    Low-purity samples should be flagged or
    excluded from depth score analysis.

IMMUNOGENIC:
  Original description (Bailey 2016):
    A subtype with high expression of immune
    pathway genes and immune cell markers.
  2024 consensus re-evaluation:
    IMMUNOGENIC is also NOT a tumour-intrinsic
    subtype.
    The immunogenic signal arises from IMMUNE
    CELL INFILTRATION in the tumour stroma.
    When immune cell deconvolution is applied,
    the tumour cell component of "immunogenic"
    samples resolves to either classical or
    basal-like PDAC.
    The immune infiltration is a MICROENVIRONMENT
    feature, not a tumour-intrinsic programme.
  FRAMEWORK IMPLICATION:
    Like luminal infiltrated MIBC or CMS1 CRC
    (where the immune signal in bulk RNA-seq
    partially reflects infiltrating T cells rather
    than tumour cell biology), the "immunogenic"
    PDAC signal is the IMMUNE MICROENVIRONMENT
    signal that co-exists with either classical
    or basal-like tumour cells.
    In the depth score analysis, immune marker
    genes (CD8A, FOXP3, CD4, CD274/PD-L1)
    must be treated as SEPARATE from the
    tumour-intrinsic depth axis.
    A two-axis model is required:
      AXIS 1 (tumour intrinsic):
        Depth = distance from acinar/ductal
        normal identity
        Markers: PTF1A, GATA6, CPA1, AMY2A
        (depth-negative) vs. TP63, KRT5,
        MYC (depth-positive)
      AXIS 2 (microenvironment):
        Immune infiltration score
        Markers: CD8A, FOXP3, CD274, CXCL9
        Prognostically relevant independently
        of tumour depth.

THE STROMAL SUBTYPE:
  Moffitt et al. (2015) separated the STROMAL
  signal from the tumour signal and found two
  stromal subtypes:
    STROMA NORMAL: Lower CAF activation,
    less desmoplasia — better prognosis.
    STROMA ACTIVATED: High TGF-β signalling,
    high CAF activation markers (ACTA2, FAP,
    COL1A1), dense desmoplasia — worse prognosis.
  The Moffitt stromal subtypes add prognostic
  information INDEPENDENT of the tumour subtype.
  Classical + stroma normal = best outcome.
  Basal-like + stroma activated = worst outcome.
  This is the most complete prognostic model
  for PDAC available in the 2024 literature.
  For the PAAD depth score analysis, the stromal
  signal represents a THIRD AXIS to be managed:
    Not the tumour identity axis
    Not the immune infiltration axis
    The CAF/desmoplasia axis (TGF-β, ACTA2, FAP)
  These three axes must be separated in bulk
  RNA-seq analysis to avoid conflation.
```

---

## SECTION VI — THE FOUR DRIVER MUTATIONS AND THE PanIN LADDER

```
PDAC has the most stereotyped mutational
progression of any solid tumour. Four genes
are mutated in virtually every case. Their
order of acquisition defines the PanIN ladder.

KRAS (>95%):
  THE INITIATING EVENT.
  Present in PanIN-1 — the earliest identifiable
  preneoplastic lesion.
  Present in normal-appearing acinar cells
  adjacent to tumours (field effect).
  KRAS G12D (~38%), G12V (~27%), G12R (~16%)
  are the most common alleles.
  G12D predominates in PDAC — different from
  G12C (lung cancer) where sotorasib/adagrasib
  are used.
  KRAS in the acinar cell drives:
    ADM (acinar-to-ductal metaplasia)
    — the first Waddington transition
    PanIN-1 formation (flat and papillary
    mucin-producing lesions in the duct)
    Suppression of PTF1A (indirectly —
    through RAF/MEK/ERK activation that
    opposes acinar TF binding to chromatin)
  KRAS alone is NOT sufficient for PDAC:
    Acinar-specific KRAS G12D in mice:
    → ADM → PanIN-1 → PanIN-2
    → but PDAC only with additional hit
    (e.g., pancreatitis, CDKN2A loss)

CDKN2A (~95%):
  THE SECOND HIT — accelerates PanIN-2.
  p16 (encoded by CDKN2A) inhibits CDK4/6.
  CDKN2A loss → CDK4/6 constitutively active
  → Rb phosphorylation → E2F release → G1/S
  progression unchecked.
  Timing: CDKN2A loss occurs at PanIN-1 → PanIN-2
  transition.
  Mechanism of loss:
    Homozygous deletion (~40%)
    Promoter methylation (~35%)
    Point mutation/small deletion (~20%)
  CDK4/6 inhibitors (palbociclib, ribociclib,
  abemaciclib) have limited single-agent
  activity in PDAC because:
    The KRAS → ERK axis bypasses Rb in PDAC
    CDK4/6 inhibition is not the rate-limiting
    step once KRAS is constitutively active
    BUT: in the CDKN2A-deleted context, CDK4/6
    inhibitors may restore some cell cycle
    control when combined with KRAS inhibitors.
    This combination is under investigation.

TP53 (~75%):
  THE LATE DRIVER — PanIN-2 → PanIN-3 transition.
  TP53 mutations in PDAC are frequently
  gain-of-function (GOF) alleles (R175H, R248W,
  R248Q) that not only lose tumour suppressor
  function but ACTIVELY promote oncogenic
  programmes:
    GOF TP53 activates:
      VEGF-A (angiogenesis)
      PI3K targets
      PDGFR (stromal activation)
      MDM2 (a paradoxical activator via
      gain-of-function interactions)
    GOF TP53 represses:
      Differentiation TFs
      Tumour suppressor networks
  GOF p53 isoforms are a therapeutic target:
    APR-246 (eprenetapopt) restores wild-type
    conformation to some GOF TP53 alleles.
    Phase III trials in TP53 GOF contexts
    (haematological malignancies) — Phase I/II
    in PDAC.
  Timing: TP53 mutation is also the event
  that enables invasiveness — the cell acquires
  the ability to breach the basement membrane
  at the PanIN-3 → invasive carcinoma step.

SMAD4 (~50%):
  THE FINAL DRIVER — enriched at PanIN-3 →
  invasive transition but also acquired earlier.
  SMAD4 encodes the central signal transducer
  of TGF-β signalling.
  In normal cells: TGF-β → SMAD4 → growth arrest
  and apoptosis (tumour suppressive).
  In SMAD4-null PDAC: TGF-β → SMAD4-independent
  signalling (via SMAD2/3 without SMAD4) →
  pro-invasive, pro-EMT, pro-metastatic outcomes.
  The switch from tumour-suppressive to tumour-
  promoting TGF-β signalling requires SMAD4 loss.
  SMAD4 as a stratification biomarker:
    SMAD4-intact: peritoneal metastasis dominant
    SMAD4-null: haematogenous (liver) metastasis
    This was confirmed in autopsy study of PDAC
    patients (Iacobuzio-Donahue et al.).
    Clinical implication: SMAD4 IHC on the
    diagnostic biopsy predicts METASTATIC PATTERN
    and therefore which patients are most likely
    to benefit from aggressive local control
    (radiation/surgery) vs. systemic therapy.

THE CANON: KRAS → CDKN2A → TP53 → SMAD4

  This four-gene sequence is the most
  well-characterised carcinogenic sequence
  of any solid tumour — better characterised
  than the APC → KRAS → TP53 → SMAD4 sequence
  in CRC (which is similar but uses APC where
  PDAC uses CDKN2A at the second position).

  CRC comparison:
    CRC: APC → KRAS → TP53 → SMAD4
    PDAC: [NORMAL] → KRAS → CDKN2A → TP53 → SMAD4
  The structural parallel is exact: both sequences
  end in TP53 and SMAD4 as late drivers.
  The difference: CRC begins with APC (WNT pathway)
  while PDAC begins with KRAS (RTK/MAPK pathway).
  Both end in the same loss of genome guardian (TP53)
  and the same TGF-β pathway failure (SMAD4).
  This cross-cancer convergence is a framework
  prediction: the final drivers of invasive
  capacity converge across cancer types even
  when the initiating drivers differ.
```

---

## SECTION VII — THE DESMOPLASTIC STROMA: THE UNIQUE PAAD CONTEXT

```
PDAC has the most extreme desmoplastic stroma
of any solid tumour. This stroma is:

  PHYSICAL COMPOSITION:
    ~60–90% non-tumour cells by volume
    CAFs (cancer-associated fibroblasts): ~40%
    Immune cells: ~10–15%
    Endothelial cells / vasculature: ~5%
    Normal pancreatic acinar/ductal remnants:
    variable (up to ~20% — THE ADEX ARTEFACT
    SOURCE)
    Extracellular matrix: collagen, fibronectin,
    hyaluronic acid (dense, stiff)
    The tumour cell content: ~10–30% of tissue

  BIOLOGICAL CONSEQUENCES:
    Hypovascular: the stroma compresses blood
    vessels → PDAC tumours have very low blood
    flow relative to mass → poor drug delivery
    → gemcitabine concentration in tumour cells
    is LOWER than in plasma by 10-fold.
    Hypoxic: low vascular density → HIF-1α
    activation → VEGF → paradoxically poor
    angiogenesis (VEGF inhibitors have not
    worked in PDAC for this reason).
    Immune-excluded: the stroma physically
    prevents T cell penetration of the tumour
    core. TGF-β from CAFs induces T cell
    exhaustion and excludes them from the
    invasive front. PDAC is the paradigm
    "immune desert" cancer.
    This is why immunotherapy does not work
    in PDAC: the T cells are not even getting
    to the tumour cells.

  CAF HETEROGENEITY (2024 insight):
    CAFs are not one cell type — two major
    subtypes have been identified:
    MYOFIBROBLASTIC CAFs (myCAFs):
      ACTA2 (alpha-SMA) HIGH
      Located adjacent to tumour cells
      Produce ECM (collagen, fibronectin)
      PRO-TUMOUR: form the physical desmoplastic
      barrier
    INFLAMMATORY CAFs (iCAFs):
      ACTA2 LOW
      Located away from tumour cells in stroma
      Produce cytokines: IL-6, IL-11, CXCL1,
      CXCL2, LIF
      Pro-inflammatory — drive systemic
      cachexia and local immunosuppression
      ALSO PRO-TUMOUR but via different mechanism
    Both CAF subtypes are maintained by TGF-β
    signalling from the tumour cells.
    TGF-β inhibition (galunisertib, vactosertib)
    targets both CAF subtypes.
    Phase I/II data: modest single-agent activity;
    combination with chemotherapy or immunotherapy
    under active investigation.

  THE STROMA PARADOX:
    The desmoplastic stroma appears to be BOTH
    tumour-promoting (by excluding immune cells
    and compressing vasculature) AND potentially
    tumour-restraining (by mechanically limiting
    invasion at early stages).
    Early anti-stroma strategies (hedgehog
    pathway inhibition with IPI-926/vismodegib
    to dissolve stroma) FAILED in Phase II trials
    and potentially ACCELERATED invasion/metastasis
    by removing the physical restraint.
    The lesson: the stroma cannot be simply
    dissolved. It must be REPROGRAMMED.
    Current strategies target specific pro-tumour
    CAF signals (TGF-β, IL-6, FAP) while
    preserving restraining stroma functions.
    This is directly analogous to the framework's
    distinction between ATTRACTOR DISSOLUTION
    (inappropriate in PDAC stroma — removes
    restraint) and CIRCUIT RESTORATION (reprogramming
    CAFs toward a less activated state).

  STROMA AND DEPTH SCORE:
    The depth score derivation from bulk PAAD
    RNA-seq MUST account for the stroma signal.
    The stroma contributes:
      TGF-β pathway genes (TGFB1, TGFB2,
      TGFBR1, SMAD2, SMAD3)
      ECM genes (COL1A1, COL1A2, FN1, VIM)
      FAP, ACTA2 (CAF markers — depth-confounded)
    These will CORRELATE with depth (more
    aggressive = more activated stroma) but are
    NOT tumour-intrinsic depth markers.
    The before-documents must separate:
      TUMOUR DEPTH AXIS (PTF1A, GATA6, CPA1
      falling; TP63, KRT5, MYC rising)
      from:
      STROMA ACTIVATION AXIS (ACTA2, FAP, COL1A1,
      TGFB1 — depth-correlated but not causal)
    This separation is the most technically
    challenging aspect of PAAD depth score
    derivation and the primary reason why
    subtype-pure analysis (using maximum tumour
    purity filtered samples) is required.
```

---

## SECTION VIII — EXISTING PAAD ANALYSIS IN CONTEXT

```
THE COMPLETED ANALYSIS (OrganismCore cancer series)
ran on a bulk PAAD dataset
(TCGA-PAAD or equivalent).

WHAT THE ANALYSIS FOUND:

  FINDING 1: PTF1A identified as switch gene.
    CORRECT CONTEXT:
    PTF1A is the master acinar TF — its loss
    is the founding event in ADM and the Waddington
    transition from normal acinar to cancer cell.
    Its depth-negative correlation (falls with
    depth) is expected: deeper tumours have lost
    more acinar identity (more PTF1A gone).
    Its identification as the switch gene confirms:
    the depth axis is measuring loss of acinar
    identity, which is the correct Waddington axis.
    SUBTYPE CONTEXT:
    PTF1A is relevant primarily to CLASSICAL
    PDAC (the dominant subtype at ~60%).
    In basal-like PDAC, PTF1A is also gone —
    but restoring it may not be sufficient to
    overcome the active TP63-driven squamous
    programme that is maintaining the basal-like
    false attractor.
    The circuit integrity finding (intact in bulk)
    reflects the CLASSICAL-DOMINATED signal.
    Whether the circuit is intact in PURE BASAL-
    LIKE samples requires PAAD-S2a.

  FINDING 2: EZH2 elevated (gain-of-function lock).
    CORRECT CONTEXT:
    EZH2 in PDAC silences PTF1A and other acinar
    identity genes via H3K27me3.
    This is the epigenetic lock that prevents
    circuit restoration from occurring
    spontaneously.
    EZH2 is the convergence node.
    Tazemetostat blocks EZH2 → H3K27me3 cleared
    → PTF1A locus accessible → PTF1A expression
    possible → acinar differentiation could
    re-execute.
    This is CONFIRMED as elevated in CLASSICAL
    PDAC bulk signal.
    SUBTYPE QUESTION:
    Is EZH2 ALSO elevated in basal-like PDAC?
    If EZH2 is elevated in both subtypes:
    tazemetostat addresses BOTH.
    If EZH2 is not the driver in basal-like
    (where TP63 is the active driver):
    tazemetostat alone is insufficient for
    basal-like.
    This question is resolved in PAAD-S2a.

  FINDING 3: Circuit restoration therapeutic
    strategy — NOT attractor dissolution.
    CORRECT CONTEXT:
    Unlike GS gastric cancer (broken circuit —
    attractor dissolution strategy) or BLCA
    basal/squamous (TP63 retention — dissolution
    strategy), classical PDAC has an intact
    downstream circuit that can be re-executed.
    The metaphor: the programme is paused,
    not deleted. Unblocking EZH2 presses play.
    This is the most therapeutically optimistic
    finding in the PAAD series — if correct,
    classical PDAC cells can be CONVERTED to
    acinar differentiation, not just killed.

  FINDING 4: Tazemetostat predicted as drug target.
    CONFIRMED: EZH2 inhibitor prediction in PDAC.
    Literature check status (from OrganismCore
    series): EZH2 is known to be overexpressed
    in PDAC and preclinical data supports
    tazemetostat activity. Clinical trials in
    PDAC for EZH2 inhibitors are in early stage.
    The framework-derived prediction that
    tazemetostat = circuit restoration agent
    (not just cytostatic) is the NOVEL component
    — not just that EZH2 is elevated, but that
    blocking it restores the downstream
    differentiation programme.
    THIS IS THE GOING FURTHER FINDING FOR PAAD.
```

---

## SECTION IX — DATA AVAILABILITY SUMMARY

```
Dataset         Accession         n(tumour)  Normal  Subtype  Power
                                             (n)     labels
────────────────────────────────────────────────────────────────────
TCGA-PAAD       GDC/phs000178     ~178       4 adj   YES      MOD
                                             normal  (TCGA +
                                                      Bailey
                                                      mapped)
ICGC-AU         PACA-AU           ~269       some    YES      HIGH
(Bailey         (controlled        adj              (ORIGINAL
original)       access)                             BAILEY
                                                     LABELS)
GSE62452        GEO               ~130       ~65     PARTIAL  HIGH
                (Moffitt 2015     (micro-           (Moffitt
                dataset)          array)            Classical/
                                                    Basal)
GSE16515        GEO               ~36        ~16     PARTIAL  LOW-MOD
GSE15471        GEO               ~36        ~16     PARTIAL  LOW-MOD
GTEx            —                 0          ~150    N/A      HIGH
                                             healthy         (NORMAL
                                             pancreas        ONLY)
GSE84133        GEO (scRNA)       normal     YES     Cell     MOD
                                  pancreas           type
                                  only               labels

CRITICAL NOTES:

NOTE 1 — TUMOUR PURITY:
  TCGA-PAAD has a mean tumour purity of ~35%
  (the lowest purity of any TCGA solid tumour
  cohort). Many samples have purity <20%.
  Low purity samples contain predominantly
  stroma/acinar/immune cells, NOT tumour cells.
  ALL DEPTH SCORE ANALYSIS IN PAAD MUST FILTER
  FOR TUMOUR PURITY ≥40%.
  TCGA provides ESTIMATE purity scores.
  Samples with purity <40% should be flagged
  or excluded from the depth score derivation.
  Failure to filter produces ADEX-artefact
  results — the depth-negative genes will be
  PTF1A/CPA1/PRSS1 from contaminating acinar
  cells, not from low-depth tumour cells.
  This is the primary methodological risk in
  PAAD depth score analysis.

NOTE 2 — NORMAL TISSUE:
  GTEx pancreas tissue (healthy donors) is the
  best normal reference for the PAAD depth axis.
  TCGA adjacent normal (n=4) is too few for
  reliable normal reference.
  If GTEx pancreas is not accessible in the
  script, TCGA-PAAD n=4 adjacent normal is
  the fallback. Flag this in the before-document.

NOTE 3 — SUBTYPE LABELS:
  Bailey subtype labels for TCGA-PAAD are
  available in the Bailey 2016 supplementary
  tables (Nature 2016).
  Moffitt subtype labels (Classical/Basal-like)
  are cleaner for two-subtype analysis and
  available in the Moffitt 2015 supplementary.
  The Moffitt two-subtype system is preferred
  for the subtype-stratified depth score because:
    It is more robust to tumour purity
    It separates tumour from stroma
    It is increasingly used in clinical trials
    (PANGEA trial uses PurIST — the Moffitt
    classifier adapted for IHC/low-input data)

NOTE 4 — THE PURITY-DEPTH CONFOUND:
  In PAAD bulk RNA-seq:
    Low tumour purity → high expression of
    acinar cell genes (PTF1A, CPA1, AMY2A)
    FROM CONTAMINATING ACINAR CELLS.
    High tumour purity → low expression of
    acinar cell genes (tumour cells have lost
    acinar identity).
  This INVERTS the expected depth relationship
  if purity is not controlled:
    In a correctly purity-filtered dataset:
    Depth-negative = acinar genes (PTF1A, CPA1)
    because TUMOUR cells with more identity
    loss have less PTF1A.
    In an unfiltered dataset:
    "Depth-negative" = acinar genes because
    HIGH-PURITY (high depth) TUMOUR SAMPLES
    have fewer contaminating acinar cells,
    producing an apparent depth correlation
    that is ENTIRELY ARTIFACTUAL.
  THE BEFORE-DOCUMENTS MUST STATE EXPLICITLY:
  Purity filtering is performed BEFORE the
  depth score is derived.
  The depth score is computed on tumour cells,
  not on contaminating acinar cell admixtures.
```

---

## SECTION X — PLANNED ANALYSIS ORDER

```
ORDER:

  PAAD-S1   Classical subtype   TCGA-PAAD         MOD-HIGH
            (dominant,          Classical-only
            intact circuit)     (purity-filtered,
                                Moffitt Classical
                                labels)
                                + GTEx normal
                                pancreas
                                REASON: The dominant
                                subtype (~60%).
                                Contains the intact
                                circuit finding from
                                the bulk analysis.
                                EZH2 elevation and
                                PTF1A depth axis
                                confirmed in pure
                                classical samples.
                                Tazemetostat
                                prediction validated
                                in classical-specific
                                context.
                                Clinical output:
                                depth score as
                                circuit restoration
                                therapy candidate
                                predictor — which
                                classical PDAC
                                patients are most
                                likely to benefit
                                from tazemetostat.

  PAAD-S2   Basal-like subtype  TCGA-PAAD         MOD
            (aggressive,        Basal-only
            broken circuit?)    (purity-filtered,
                                Moffitt Basal-like
                                labels)
                                + GTEx normal
                                REASON: The second
                                subtype. Circuit
                                integrity question:
                                is EZH2 the
                                convergence node in
                                basal-like as in
                                classical, or is
                                TP63 the primary
                                driver?
                                AURKA depth
                                correlation in
                                basal-like — is it
                                higher than in
                                classical?
                                Alisertib prediction
                                for basal-like.
                                KRAS inhibitor
                                resistance via
                                subtype switching:
                                test whether
                                classical markers
                                RISE when basal-like
                                depth score falls
                                (the escape mechanism).
                                Clinical output:
                                basal-like drug map
                                (TP63 inhibitor +
                                AURKA inhibitor vs.
                                classical EZH2
                                inhibitor circuit
                                restoration path).

  PAAD-S3   Cross-subtype       After S1–S2.
            + KRAS inhibitor    QUESTIONS:
            escape analysis       1. Does depth score
                                     separate Classical
                                     from Basal-like
                                     within purity-
                                     filtered samples?
                                  2. EZH2 direction:
                                     elevated in both
                                     subtypes or only
                                     classical?
                                  3. PTF1A recovery:
                                     does tazemetostat
                                     re-induce PTF1A
                                     expression in
                                     classical but not
                                     basal-like?
                                     (testable in cell
                                     line data — GSE
                                     search needed)
                                  4. SMAD4 status and
                                     depth: do SMAD4-
                                     null tumours cluster
                                     at deeper positions?
                                  5. KRAS inhibitor
                                     escape: in the
                                     GSE271300 dataset
                                     (resistance to
                                     KRAS inhibition),
                                     do escape cells
                                     show GATA6 recovery
                                     (classical identity
                                     return)?

  PAAD-X    Stroma analysis     After S1–S3.
                                STROMA AXIS:
                                  Separate the
                                  stroma activation
                                  axis (ACTA2, FAP,
                                  TGFB1) from the
                                  tumour depth axis.
                                  Test the dual-axis
                                  model:
                                  Axis 1 (tumour
                                  depth) × Axis 2
                                  (stroma activation)
                                  — does the product
                                  improve OS prediction
                                  over tumour depth
                                  alone?
                                  This is the Moffitt
                                  dual model (tumour
                                  subtype × stroma
                                  subtype) tested in
                                  the depth score
                                  framework.
```

---

## SECTION XI — WHAT THIS DOCUMENT DOES NOT CONTAIN

```
This document contains:
  ✓ Normal exocrine pancreatic hierarchy
    (acinar → centroacinar → ductal, with
    full TF network: PTF1A/NR5A2/RBPJL/MIST1)
  ✓ The ADM → PanIN → PDAC Waddington cascade
    (KRAS → CDKN2A → TP53 → SMAD4)
  ✓ The two robust subtypes (Classical and
    Basal-like) with cells of origin, molecular
    events, and treatment
  ✓ The two contested subtypes (ADEX and
    Immunogenic) with the contamination
    explanation — their invalidity as tumour-
    intrinsic subtypes
  ✓ The circuit integrity finding from the
    existing analysis placed in subtype context
    (intact circuit = classical signal; basal-
    like circuit integrity requires PAAD-S2a)
  ✓ The PTF1A/EZH2 depth axis and tazemetostat
    prediction in context
  ✓ GATA6 as the single separating marker
  ✓ The tumour purity warning (PDAC has the
    lowest tumour purity in TCGA — purity
    filtering is mandatory before depth score)
  ✓ Desmoplastic stroma biology and CAF subtypes
    (myCAFs vs. iCAFs) and the stroma paradox
  ✓ KRAS inhibitor resistance via subtype
    switching (basal-like → classical escape)
  ✓ The cross-repository structural parallels
    (basal-like PDAC ↔ basal/squamous BLCA ↔
    GS STAD ↔ CMS4 CRC — identity-loss attractor
    type across GI cancers)
  ✓ Data availability with the purity-filtering
    note

This document does NOT contain:
  ✗ Depth score predictions
  ✗ Specific false attractor gene predictions
  ✗ Drug target predictions beyond what the
    existing analysis already established
  ✗ Epigenetic lock mechanism hypotheses
  ✗ Novel finding claims

All of the above belong in the BEFORE documents.
PAAD-S1a (Classical before-document) is next.
Written before any script runs.
Before any data loads.
```

---

## STATUS BLOCK

```
document:           PAAD_Subtype_Orientation.md
folder:             Cancer_Research/PAAD/Subtypes/
status:             COMPLETE — ORIENTATION ONLY
date:               2026-03-03
author:             Eric Robert Lawson
                    OrganismCore

subtypes_covered:
  Classical (~60%):    GATA6+, PTF1A lost,   [1 of 2]
                       ductal false attractor,
                       intact circuit, EZH2 elevated
  Basal-like (~35%):   GATA6−, TP63+, KRT5+, [2 of 2]
                       identity-loss attractor,
                       circuit status TBD
  (ADEX, Immunogenic:  ARTEFACTS — not tumour-
  intrinsic subtypes — documented and excluded)

analyses_started:   0 (new subtype series)
existing_analysis:  Cancer_Research/PAAD/ — complete
                    (OrganismCore cancer series)
                    PTF1A switch gene identified
                    EZH2 elevated (gain-of-function lock)
                    Circuit intact (restoration strategy)
                    Tazemetostat predicted as drug target
                    All confirmed findings in context
                    of Classical-dominated bulk signal

next_document:      PAAD-S1a
                    Classical Subtype Before-Document
                    (predictions locked before
                    purity-filtered TCGA-PAAD Classical
                    subset loads)

critical_note_1:    TUMOUR PURITY IS THE PRIMARY
                    METHODOLOGICAL RISK IN PAAD.
                    TCGA-PAAD has the lowest mean
                    tumour purity (~35%) of any TCGA
                    solid tumour cohort. Low-purity
                    samples produce the ADEX artefact:
                    acinar cell contamination genes
                    (PTF1A, CPA1, AMY2A) appear as
                    depth-positive signals when they
                    are contamination signals.
                    ALL analysis in the PAAD subtype
                    series must begin with ESTIMATE-
                    based purity filtering (purity ≥40%)
                    before any depth score derivation.
                    This is not optional. This is the
                    difference between measuring tumour
                    biology and measuring normal tissue
                    contamination.

critical_note_2:    The circuit restoration finding
                    (Pattern 5: intact circuit in PAAD)
                    is the most therapeutically important
                    result in the PAAD series.
                    Its subtype specificity matters:
                    Classical PDAC = intact circuit
                    (EZH2 blocks PTF1A restoration but
                    the downstream programme is connected)
                    Basal-like PDAC = circuit status
                    unknown — TP63 is actively driving
                    a competing programme, not just
                    blocking a dormant one.
                    The distinction determines whether
                    PAAD-S2 concludes:
                    a) Both subtypes benefit from
                       tazemetostat (EZH2 elevated in
                       both, circuit intact in both)
                    b) Only classical benefits from
                       tazemetostat (basal-like needs
                       TP63 inhibitor as the primary
                       intervention)
                    This is the most important binary
                    outcome of the entire PAAD subtype
                    series.

critical_note_3:    The KRAS inhibitor resistance via
                    subtype switching is the most novel
                    structural insight available in the
                    PAAD literature as of 2026.
                    When basal-like PDAC cells are
                    exposed to KRAS G12D inhibitors,
                    they escape by ADOPTING CLASSICAL
                    IDENTITY (GATA6 rises, KRT5 falls).
                    This is an attractor landscape
                    transition — the basal-like false
                    attractor is destabilised and the
                    cell slides to the adjacent classical
                    false attractor.
                    The framework prediction from this:
                    KRAS inhibitor + EZH2 inhibitor
                    simultaneously — prevents the
                    classical attractor from serving as
                    an escape route (EZH2i destabilises
                    classical attractor by allowing
                    PTF1A to return → cell exits false
                    attractor entirely → differentiation
                    or death, not escape to classical).
                    This combination has not been
                    tested clinically.
                    It is the framework's primary novel
                    prediction for the PAAD subtype
                    series.

critical_note_4:    The cross-repository structural
                    identity of the basal-like pattern
                    across GI cancers:
                    GS STAD (CDH1 loss, signet ring)
                    Basal/Squamous BLCA (TP63+, KRT5+)
                    Basal-like PDAC (TP63+, GATA6−)
                    CMS4 CRC (mesenchymal, TGF-β)
                    All four are the MOST IDENTITY-LOST
                    false attractor in their respective
                    GI tissue of origin.
                    All four have the worst prognosis
                    in their cancer type.
                    All four are immune-excluded.
                    All four are chemotherapy-resistant.
                    The framework sees them as the SAME
                    WADDINGTON ATTRACTOR TYPE in
                    different tissue coordinates.
                    The drug target derivation for the
                    most identity-lost attractor type
                    should therefore show convergent
                    features across all four:
                    Identity markers that can be targeted
                    (CLDN18.2 in STAD, NECTIN4 in BLCA,
                    EGFR in PDAC basal-like, ZEB1 in
                    CMS4 CRC) represent the RESIDUAL
                    SURFACE IDENTITY that survives even
                    in the most identity-lost false
                    attractor — and each is the ADC or
                    antibody target for that tissue's
                    worst-prognosis subtype.
                    This cross-cancer convergence is
                    Framework Pattern 2 demonstrated
                    at the subtype resolution level.
```
