# STOMACH ADENOCARCINOMA (STAD) — SUBTYPE ORIENTATION DOCUMENT
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

  A complete map of the gastric cancer molecular
  subtype landscape — what the four TCGA molecular
  subtypes are, where each arises in the normal
  gastric glandular hierarchy, what the initiating
  and cooperating genetic events are, what the
  clinical characteristics and treatment are,
  and what public data exists to analyze each.

STAD holds a position in this repository that
requires explicit contextualisation:

  THE EXISTING ANALYSIS IS ALREADY DONE.

  The STAD bulk analysis (GSE32571 or equivalent)
  has already been run as part of the original
  OrganismCore cancer series. The findings from
  that analysis — ZEB2/AURKA co-expression,
  WNT5A depth correlation, ZEB2 r=+0.56 with
  depth, CDH2/VIM suppression (the ANALYST
  ASSUMPTION ERROR), the EZH2 suppression
  finding, the broken differentiation circuit —
  are documented in the OrganismCore series
  (Documents 70–84 lineage).

  This subtype orientation document therefore
  serves a different function than in other
  cancer series:

    1. It places the existing STAD findings in
       the context of the four TCGA subtypes —
       specifically, it identifies WHICH subtype
       dominated the existing bulk analysis.

    2. It maps the existing findings onto each
       subtype's biology to assess which are
       subtype-specific and which are cross-
       subtype signals.

    3. It frames the new subtype-specific analyses
       that must be run to resolve the subtype
       heterogeneity that the bulk analysis
       could not distinguish.

  THE EXISTING ANALYSIS CAPTURED:
    The combined TCGA-STAD signal is dominated
    by CIN (~50% of cases) — the most common
    subtype. ZEB2/AURKA co-expression, the
    loss of CDH1/MUC5AC (gastric surface
    markers), and the mesenchymal/EMT signal
    are consistent with CIN biology.
    However, the CDH2/VIM SUPPRESSION (the
    analyst assumption error — EMT NOT being
    the dominant bulk signal) is consistent with
    the GS subtype being a MINORITY signal:
    GS is the diffuse/mesenchymal subtype where
    CDH1 is lost (consistent) but also where
    CDH2/VIM are NOT systematically elevated
    (because the GS cell is not a classic
    mesenchymal EMT cell — it is a signet ring
    cell or poorly cohesive carcinoma cell
    with a completely different biology).
    The existing analysis correctly identified
    that the BULK STAD signal is NOT an EMT
    signal — because GS (the EMT-like subtype)
    is a minority and its transcriptome is
    DILUTED by the dominant CIN signal.

    THIS IS THE CORE SUBTYPE PROBLEM IN STAD:
    Treating all gastric cancer as one disease
    produces a depth axis that is a MIXTURE
    of four completely different biological
    programmes. The subtype-specific analyses
    will unmix this.
```

---

## DOCUMENT METADATA

```
document_id:        STAD_Subtype_Orientation
series:             STAD (Stomach Adenocarcinoma — Subtypes)
folder:             Cancer_Research/STAD/Subtypes/
date:               2026-03-03
author:             Eric Robert Lawson
status:             ORIENTATION ONLY — no predictions,
                    no analysis, no epigenetic hypotheses
next_document:      STAD_CIN_before.md
                    (Document STAD-S1a — CIN before-doc)
protocol_version:   Workflow_Protocol.md v2.0
existing_analysis:  See OrganismCore Documents 70–84
                    lineage for the completed STAD bulk
                    analysis (ZEB2/AURKA, WNT5A, broken
                    circuit finding, EZH2 suppression).
```

---

## SECTION I — THE NORMAL GASTRIC GLAND

```
The Waddington baseline for gastric cancer is the
normal gastric glandular epithelium — a complex,
regionalized tissue with a distinct cell hierarchy
unlike any other epithelium in the body.

THE STOMACH IS NOT UNIFORM.
The stomach has three anatomically and
molecularly distinct regions:
  CARDIA:    Proximate to the oesophagus.
             Mucous cells, minimal parietal cells.
             CIN-type gastric cancer and GEJ
             (gastroesophageal junction) cancer
             arises here predominantly.
  CORPUS     The main body — the acid-secreting
  (FUNDUS):  zone. Contains the full glandular
             hierarchy described below.
             Parietal cells abundant.
             This is where the normal Waddington
             baseline is best defined.
  ANTRUM     The pyloric region — mucous-secreting,
  (PYLORUS): minimal parietal cells.
             G cells (gastrin-secreting) abundant.
             Most MSI and GS subtypes arise in
             the antrum.
             Most H. pylori-driven carcinogenesis
             begins in the antrum.

═══════════════════════════════════════════════════════
THE CORPUS GASTRIC GLAND — FROM SURFACE TO BASE
═══════════════════════════════════════════════════════

THE FOVEOLAR ZONE (Surface — facing the gastric lumen)
  Cells:      SURFACE MUCOUS (FOVEOLAR) CELLS
              Columnar, mucus-filled
              Coat the gastric surface with a
              protective mucous layer
  Markers:    MUC5AC — HIGH (the defining mucin
              of the stomach surface; its loss
              is a marker of metaplastic
              transformation)
              TFF1 (trefoil factor 1) — HIGH
              (protective, also a tumour suppressor
              in stomach; TFF1 KO mice develop
              gastric cancer)
              MKI67 (Ki67): LOW — post-mitotic
  Function:   Barrier against HCl and pepsin.
              Short-lived cells — renewed every
              ~3–5 days from the isthmus below.

THE ISTHMUS ZONE (The proliferative zone)
  Cells:      STEM CELLS + TRANSIT-AMPLIFYING
              PROGENITORS
              The isthmus is the gastric equivalent
              of the colonic crypt base — where
              all cell turnover originates.
  Markers:    MKI67 (Ki67): HIGH — most mitotic
              CD44: present (stem marker)
              AXIN2: WNT target — stem activity
              SOX9: progenitor TF
              LGR5: marks a stem population
                    (context-dependent)
  Function:   Multipotent stem cell niche.
              Produces:
              UPWARD → foveolar cells
              DOWNWARD → neck cells → chief cells
              LATERAL → parietal cells

THE NECK ZONE
  Cells:      MUCOUS NECK CELLS
              Wedge-shaped, MUC6+, TFF2+
              The progenitor of the chief cell
              (through transdifferentiation in
              the base zone)
              Also the cell that undergoes
              SPEM in response to injury
  Markers:    MUC6 — HIGH (distinct from
              the surface MUC5AC)
              TFF2 (spasmolytic polypeptide) — HIGH
              MIST1 — beginning to rise
              (MIST1 levels rise as the mucous
              neck cell matures into a chief cell
              — MIST1 drives chief cell identity)

THE BASE ZONE
  Cells:      CHIEF CELLS (predominant)
              + PARIETAL CELLS
              + ENTEROENDOCRINE CELLS

  CHIEF CELLS:
    The terminally differentiated secretory
    cell of the gastric gland base.
    Produce: PEPSINOGEN (PGC) — the gastric
             protease zymogen activated to pepsin
             at low pH
    Markers: MIST1 (BHLHA15) — THE MASTER TF
                               OF CHIEF CELL
                               IDENTITY
             PGC (pepsinogen C)
             GIF (gastric intrinsic factor)
             LIPF (lipase)
             MUC6 — some
    Self-renewal: ZERO — post-mitotic.
    The most terminally differentiated cell
    in the gastric gland. Once a chief cell,
    always a chief cell (under normal conditions).
    The chief cell is the terminal Waddington
    attractor of the gastric corpus.
    MIST1 drives chief cell identity by
    activating secretory and exocytosis genes
    and simultaneously REPRESSING progenitor
    identity genes.
    MIST1 IS THE CDX2 OF THE GASTRIC CORPUS —
    the guardian TF of terminal differentiated
    gastric identity. Its loss is a marker of
    dedifferentiation in CIN and GS gastric
    cancer.

  PARIETAL CELLS:
    The HCl-secreting cells.
    One of the most energy-demanding cells
    in the body — H+/K+ ATPase pumps protons
    against enormous concentration gradients.
    Markers: ATP4A, ATP4B — the H+/K+ ATPase
             subunits (absolute parietal cell
             markers — not expressed anywhere
             else in the body)
             CLDN18 — HIGH
             (Claudin-18 is the tight junction
             protein of parietal cells and
             gastric neck cells.
             CLDN18.2 — the stomach-specific
             isoform — is both a normal gastric
             identity marker AND the target of
             zolbetuximab in GS/CIN gastric cancer.
             CLDN18.2 RETENTION in cancer cells
             is the basis for the drug target.
             CLDN18.2 LOSS (in deep tumours)
             predicts non-response to zolbetuximab —
             this is the predicted depth relationship
             to be tested.)
    Self-renewal: LOW — long-lived but not post-
                  mitotic. Some capacity for
                  renewal in response to injury.

  ENTEROENDOCRINE CELLS:
    Scattered throughout the gland base.
    G cells (antrum): secrete gastrin (GAST)
    ECL cells (corpus): secrete histamine
    D cells: secrete somatostatin (SST)
    Enterochromaffin cells: serotonin
    These cells are the neuroendocrine compartment
    of the stomach — they are relevant to the
    occasional gastric neuroendocrine tumour
    (gastric carcinoid) and to the neuronal
    features that emerge at extreme depth.

CRITICAL NORMAL MARKERS FOR THE DEPTH AXIS:

  GASTRIC IDENTITY MARKERS (normal, depth-negative
  in cancer — LOST as depth increases):
    MUC5AC:  Surface mucous (foveolar) identity
    TFF1:    Tumour suppressor + foveolar marker
    ATP4A:   Parietal cell — absolute gastric marker
    CLDN18:  Gastric tight junction marker
    MIST1:   Chief cell identity master TF
    PGC:     Chief cell secretory product
    GIF:     Chief cell secretory product

  PROGENITOR/CANCER MARKERS (depth-positive
  in cancer — RISE as depth increases):
    CDH2 (N-cadherin): Rises in some subtypes
                       but NOT in bulk STAD
                       (the ANALYST ASSUMPTION ERROR)
    VIM: Similarly context-dependent
    ZEB2: Rises with depth in CIN bulk signal
    MKI67: Rises (proliferation marker)
    AURKA: Rises with depth in CIN

THE H. PYLORI → PRECURSOR → CANCER CASCADE:

  The Correa cascade maps the Waddington
  landscape of gastric carcinogenesis:

  Normal gastric mucosa (MUC5AC+, ATP4A+, MIST1+)
    ↓  H. pylori → chronic inflammation → NF-κB
  Chronic superficial gastritis
    ↓  progressive gland loss
  Atrophic gastritis (parietal cells lost —
     ATP4A falls, acid production drops)
    ↓  injury response → chief cell transdifferentiation
  SPEM (Spasmolytic Polypeptide-Expressing
     Metaplasia: TFF2 high, chief cell markers fall,
     MUC6 persists)
    ↓  CDX2 expression begins (intestinal TF)
  Intestinal Metaplasia (CDX2+, MUC2+, goblet cells
     present — the stomach has partially
     reprogrammed to intestinal identity)
    ↓  additional mutations (TP53, KRAS, etc.)
  Dysplasia → Carcinoma

  THIS CASCADE IS THE WADDINGTON LANDSCAPE OF
  GASTRIC CANCER FORMATION:
  The false attractor in gastric cancer IS the
  product of a Waddington transition — from gastric
  identity (MIST1+, MUC5AC+, ATP4A+) through an
  intermediate metaplastic state (SPEM, IM) to
  cancer (CDX2+, MUC2+ in intestinal types, or
  CDH1-lost in diffuse type).
  The depth axis in gastric cancer ultimately
  measures how far the cell has traveled from
  the normal gastric glandular identity —
  how much MUC5AC, MIST1, and ATP4A it has
  lost and how much it has gained either
  intestinal (CDX2, MUC2) or mesenchymal/
  progenitor (ZEB2, AURKA, CDH2) character.
```

---

## SECTION II — THE FOUR TCGA SUBTYPES: STRUCTURAL MAP

```
The TCGA 2014 comprehensive molecular
characterisation of gastric adenocarcinoma
identified four molecular subtypes with
distinct biology, prognosis, and therapeutic
vulnerabilities. They are not equal in frequency
and they are NOT all on the same depth axis —
they represent different Waddington transitions
from the normal gastric gland.

SUBTYPE FREQUENCY (TCGA):
  CIN (Chromosomal Instability):  ~50%
  GS (Genomically Stable):        ~20%
  MSI (Microsatellite Instable):  ~20%
  EBV (Epstein-Barr Virus+):      ~9%

DEPTH ORDERING (THE FRAMEWORK'S STRUCTURAL CLAIM):
  The four subtypes are NOT simply ordered on
  a single depth axis. They represent different
  false attractors in the gastric Waddington
  landscape.

  EBV and MSI: The IMMUNE-INFILTRATED attractors.
    Both are hypermethylated (EBV extremely so).
    Both are immunotherapy-responsive.
    The false attractor engine is EPIGENETIC
    SILENCING (CIMP in EBV, MLH1 silencing in MSI).
    These are analogous to CRC CMS1 —
    immune hot, high mutation burden, responsive
    to checkpoint blockade, BETTER prognosis
    within their stage despite being cancer.

  CIN: The CHROMOSOMAL INSTABILITY attractor.
    Aneuploidy, TP53 mutation, RTK amplification.
    The intestinal-type Lauren pathway.
    Analogous to CRC CMS2 — chromosomally
    unstable, intermediate prognosis, RTK-driven.
    This is the dominant signal in bulk STAD analysis.

  GS: The IDENTITY-LOSS attractor.
    CDH1 or RHOA mutation / CLDN18-ARHGAP fusion.
    Diffuse type — signet ring cells or poorly
    cohesive carcinoma.
    The cell that has most catastrophically
    lost normal gastric identity.
    The deepest false attractor in the gastric
    Waddington landscape.
    Worst prognosis of all four subtypes.
    Analogous to CRC CMS4 — mesenchymal,
    stromal, immune-excluded, resistant.
    BUT — see critical note: GS is NOT a
    classic EMT cancer. It is a LOSS OF
    COHESION cancer, not a gain of
    mesenchymal identity cancer.
    This is the source of the ANALYST
    ASSUMPTION ERROR in the original STAD
    analysis (CDH2/VIM suppressed, not elevated,
    in bulk STAD signal dominated by non-GS
    subtypes mixing with GS).

LAUREN CLASSIFICATION ALIGNMENT:
  Intestinal type: CIN and MSI subtypes
  Diffuse type:    GS subtype
  Mixed type:      Overlap / transitional cases
  EBV:            Mostly intestinal but distinct

ANATOMICAL PREFERENCE:
  Cardia / GEJ:    CIN (intestinal type, HER2+)
  Corpus:          EBV (proximal gastric body)
  Antrum / Pylorus: MSI (older patients, antrum)
                    GS (antrum and body, younger)
```

---

## SECTION III — CIN SUBTYPE: THE DOMINANT SIGNAL

```
FREQUENCY:     ~50% of TCGA-STAD
               The most common gastric cancer
               molecular subtype.
               Predominantly located at GEJ
               (gastroesophageal junction)
               and gastric cardia.
               Mean age: ~63 years.
               Male predominance.

PROGNOSIS:     INTERMEDIATE
               5-year OS: ~35–45% for resectable
               Molecular subtype with best
               response to standard oxaliplatin/
               fluoropyrimidine chemotherapy.

CELL OF ORIGIN:
  Gastric surface mucous (foveolar) cell or
  intermediate progenitor — intestinal type
  in Lauren classification.
  The CIN tumour often arises via the
  IM → dysplasia route (Correa cascade),
  acquiring chromosomal abnormalities at the
  dysplasia → carcinoma transition.
  Some CIN cases (GEJ location) may arise
  without the H. pylori/IM route —
  particularly those with reflux-related
  carcinogenesis.

DEFINING MOLECULAR EVENTS:
  TP53 mutation:      ~73% — the highest
                      TP53 mutation rate of
                      any TCGA STAD subtype
                      (identifies CIN)
  Chromosomal instability: EXTENSIVE
    Arm-level amplifications and deletions
    Multiple chromosomal breakpoints
    This is the molecular feature that
    defines the subtype name.
  RTK amplifications (FOCAL):
    ERBB2 (HER2):     ~17% of CIN
                      (the highest HER2 rate
                      in any gastric subtype)
                      The clinically actionable
                      target — trastuzumab
                      benefit is CIN-enriched
    EGFR:             ~12% of CIN
    FGFR2:            ~7% of CIN
                      (Bemarituzumab target —
                      Phase III data expected)
    VEGFA:            ~8% of CIN
                      (Ramucirumab target)
    MET:              ~4% of CIN
    KRAS:             ~11% (mutations)
    PIK3CA:           ~12%
  CCND1 amplification: ~6%
    (CDK4/6 inhibitor potential — but this
    is the finding that the framework
    document FAILED TO IMPLEMENT correctly:
    CDK4/6 inhibitors work in CCND1-amplified
    CIN gastric cancer, NOT in the overall
    population. Treating "gastric cancer" as
    one entity and concluding CDK4/6 inhibitors
    don't work — the paradigm failure that
    the OrganismCore framework was explicitly
    built to address, as stated in the
    OrganismCore_Cancer_Framework.md document.)
  E-cadherin (CDH1): ~6% mutation
    (Lower than GS; when CDH1 is mutated in
    CIN it tends to be focal, not the defining
    event)

EXPRESSION PROGRAMME (the bulk STAD signal):
  DEPTH-POSITIVE (rises with depth):
    ZEB2 — the confirmed depth correlate
            from the existing STAD analysis
    AURKA — confirmed co-expressed with ZEB2
            (r=0.9871 in the existing analysis)
    MKI67 — proliferation
    TOP2A, PCNA — replication markers
    MYC targets
    E2F targets
  DEPTH-NEGATIVE (falls with depth):
    MUC5AC — the primary gastric surface marker
    TFF1 — tumour suppressor / foveolar marker
    ATP4A — parietal cell marker
    CLDN18 — gastric tight junction marker
    MIST1 / PGC — chief cell markers
    CDH1 — E-cadherin (epithelial cohesion)

TREATMENT:
  Fit resectable:
    Perioperative FLOT (5-FU/leucovorin +
    oxaliplatin + docetaxel) — standard in
    Europe (FLOT4 trial)
    Or: FOLFOX / CAPOX
  Advanced/metastatic:
    First line: FOLFOX or XELOX +
    Nivolumab (CheckMate-649 — CPS ≥5)
    OR:
    Trastuzumab + FOLFOX/CAPOX (HER2+, ~17%
    of CIN — ToGA trial, landmark approval)
    + Pembrolizumab in HER2+/PD-L1+ cases
    (KEYNOTE-811 trial — improved PFS/OS)
    Second line:
    Ramucirumab + paclitaxel (REGARD trial)
  HER2+ CIN:
    The most actionable patient group in STAD.
    Trastuzumab → trastuzumab deruxtecan (T-DXd)
    at progression (DESTINY-Gastric01/02)
    T-DXd is an ADC (trastuzumab + topoisomerase
    I inhibitor payload) — 51% ORR in HER2+
    GC after prior trastuzumab. One of the
    highest response rates in any gastrointestinal
    cancer.

THE ZEB2/AURKA FINDING IN CONTEXT:
  The existing OrganismCore STAD analysis found:
    ZEB2 r=+0.56 with depth
    AURKA r correlated with depth
    ZEB2/AURKA r=0.9871 (97.4% shared variance)
  These are CIN-subtype-enriched signals:
    ZEB2 in CIN is a marker of:
      Loss of gastric epithelial cohesion
      (not full EMT — explains why CDH2/VIM
      were NOT elevated as predicted)
      Acquisition of a partially
      mesenchymal transition state
      Contributing to invasion capacity
      WITHOUT full mesenchymal conversion
      (This is PARTIAL or INCOMPLETE EMT —
      the most dangerous form for metastasis)
    AURKA in CIN is a marker of:
      Chromosomal instability maintenance
      (Aurora A drives centrosome amplification
      and aneuploid division)
      Proliferation
      Mitotic checkpoint override
  The alisertib (Aurora A inhibitor) prediction
  from the existing analysis is CORRECT for the
  CIN subtype — AURKA drives the chromosomal
  instability programme and is a depth-correlated
  target in the dominant STAD subtype.
  This prediction needs subtype stratification
  to confirm which patients benefit: CIN > MSI >
  EBV > GS (in approximate order of AURKA
  relevance).
```

---

## SECTION IV — GS SUBTYPE: THE IDENTITY-LOSS ATTRACTOR

```
FREQUENCY:     ~20% of TCGA-STAD
               ENRICHED in late-stage disease
               (~36% of late-stage) because
               GS spreads peritoneally and
               is often diagnosed at advanced
               stage.
               Younger mean age (~55 years)
               than other subtypes.
               Strong female enrichment
               (diffuse-type gastric cancer is
               more common in women relative
               to men compared to intestinal type).
               Often presents as:
               Signet ring cell carcinoma (SRC)
               Poorly cohesive carcinoma (PCC)

PROGNOSIS:     WORST of the four subtypes.
               5-year OS: ~20–25% even at
               resectable stages.
               High peritoneal metastasis rate.
               Chemotherapy resistance.
               Almost no actionable targets
               in the current armamentarium.

CELL OF ORIGIN:
  The diffuse type does NOT reliably follow the
  Correa cascade via intestinal metaplasia.
  Two routes:
  a) CDH1 germline mutation (Hereditary Diffuse
     Gastric Cancer — HDGC):
     Caused by CDH1 germline loss-of-function.
     The second hit (somatic CDH1 mutation or
     LOH) transforms normal gastric mucous neck
     cells or stem cells into signet ring cells
     WITHOUT passing through IM or dysplasia.
     The GS false attractor arises DIRECTLY
     from a normal gastric progenitor.
  b) RHOA mutation (~25% of GS):
     Activating RHOA (Y42C, L57V, E40K mutations)
     drives cytoskeletal reorganisation and
     cell roundness/non-cohesion — producing
     the signet ring phenotype from a normal
     gastric epithelial cell.
  c) CLDN18-ARHGAP fusions (~15% of GS):
     CLDN18 (the parietal cell marker) is
     fused to ARHGAP26 or ARHGAP6 — both
     RhoGAP proteins that inactivate RHOA.
     The fusion produces the same downstream
     effect as activating RHOA mutation:
     Rho pathway disruption → cytoskeletal
     collapse → signet ring morphology.
     Paradoxically: CLDN18 is now PART OF
     AN ONCOGENIC FUSION — and the residual
     CLDN18.2 expression on the cancer cell
     makes it the target of zolbetuximab.

THE SIGNET RING CELL BIOLOGY:
  A signet ring cell is a gastric cancer cell
  that has:
    Lost CDH1 (E-cadherin) — no epithelial
    adhesion — cannot form glands
    Accumulated INTRACELLULAR MUCIN that
    pushes the nucleus to one side
    (the signet ring appearance — the mucin
    vacuole is the "ring" and the displaced
    nucleus is the "signet")
    No basement membrane attachment —
    individually infiltrating cells
    Lost collective migration — individual
    cell infiltration through stroma
    PARADOXICALLY: often MUC5AC+ or MUC6+
    (retains some gastric surface/neck mucin
    expression despite catastrophic structural
    loss of gastric identity)
  This is NOT a mesenchymal cell.
  CDH2 (N-cadherin) is NOT elevated in GS SRC.
  VIM is NOT elevated in GS SRC.
  The cell has not undergone EMT in the classic
  sense — it has undergone LOSS OF EPITHELIAL
  COHESION without GAIN OF MESENCHYMAL IDENTITY.
  THIS IS THE PRECISE EXPLANATION FOR THE
  ANALYST ASSUMPTION ERROR IN THE ORIGINAL
  STAD ANALYSIS.
  The prediction was: CDH2/VIM elevated
  (standard EMT logic).
  The data showed: CDH2/VIM SUPPRESSED.
  The correction: STAD bulk signal is NOT EMT.
  The GS subtype is EMT-like in CDH1 loss but
  NOT in CDH2/VIM gain. The CIN subtype (dominant)
  shows partial ZEB2-driven transition but not
  full mesenchymal conversion.
  The framework correctly identified this as an
  ANALYST ASSUMPTION ERROR and documented it.
  This subtype orientation document provides
  the molecular explanation for why the
  assumption was wrong.

DEFINING MOLECULAR EVENTS:
  CDH1 mutation/deletion:  ~37% of GS
                           (highest CDH1 mutation
                           rate of any TCGA subtype)
                           Loss of E-cadherin =
                           loss of epithelial cohesion
  RHOA mutation:           ~25% of GS
                           (Y42C most common)
                           Cytoskeletal disruption
                           → signet ring phenotype
  CLDN18-ARHGAP fusion:    ~15% of GS
                           (and ~5% of other subtypes)
                           CLDN18 fused to ARHGAP26/6
                           → RhoGAP dysregulation
  CDH1, RHOA, and CLDN18-  LARGELY MUTUALLY
  ARHGAP:                  EXCLUSIVE —
                           alternative paths to the
                           same cytoskeletal/cohesion
                           false attractor endpoint.
  TP53 mutation:           ~12% (LOW — contrast with
                           CIN at 73%)
                           GS maintains the genome
                           relatively stably —
                           hence "genomically stable"
                           despite catastrophic
                           phenotypic change.
  RTK amplifications:      RARE — no HER2, minimal
                           EGFR, FGFR2 amplification.
                           This is why HER2 therapy
                           does not work in GS.
                           This is why the treatment
                           landscape for GS is so poor.

TREATMENT (the hardest problem in STAD):
  Standard chemotherapy: EOX / FLOT / FOLFOX
  Response rate: POOR (worse than CIN or MSI)
  Surgical resection with clear margins:
    Technically difficult because GS spreads
    as individual cells through the stroma.
    R0 rate lower than CIN.
    Peritoneal seeding is common at diagnosis
    even when not visible radiologically.
  Zolbetuximab (IMAB362):
    Anti-CLDN18.2 monoclonal antibody.
    Phase III data (SPOTLIGHT trial, GLOW trial):
      SPOTLIGHT: Zolbetuximab + mFOLFOX6
        PFS: 10.6 vs. 8.7 months (HR 0.75)
        OS:  18.2 vs. 15.5 months (HR 0.75)
      GLOW: Zolbetuximab + CAPOX
        Similar outcomes.
    ELIGIBILITY: CLDN18.2 ≥2+ in ≥75% of tumour
    cells by IHC.
    This test enriches for GS subtypes but
    is not restricted to them.
    The DEPTH SCORE PREDICTION:
    CLDN18.2 expression should fall with
    depth as gastric identity is lost.
    If CLDN18 is depth-NEGATIVE:
    Shallow (less-advanced) GS tumours
    are the best zolbetuximab candidates.
    Deep (most-advanced) GS tumours have
    lost CLDN18 and are NOT zolbetuximab
    candidates — they need a different drug.
    This is stated here as a structural
    hypothesis — to be tested in STAD-S2a.
  Anti-FGFR2:
    Bemarituzumab (anti-FGFR2b) —
    enriched in GS (FGFR2 amplification rare
    but FGFR2b overexpression from isoform
    switching IS a GS feature).
    Phase III (FORTITUDE-101) data expected 2025.
  Immunotherapy: POOR for GS.
    GS tumours have low TMB, low PD-L1 CPS,
    and immune-excluded microenvironment.
    The immune cells cannot penetrate the
    desmoplastic/stromal barrier.

HEREDITARY DIFFUSE GASTRIC CANCER (HDGC):
  ~1–3% of all gastric cancer is HDGC.
  CDH1 germline pathogenic variant.
  Lifetime risk of diffuse gastric cancer: ~70%
  Also increased risk of lobular breast cancer
  in CDH1 germline carriers.
  Prophylactic total gastrectomy is recommended
  for CDH1 germline carriers who have confirmed
  lobular adenocarcinoma foci on biopsy.
  This is the BRCA equivalent for gastric cancer
  — a germline variant that mandates prophylactic
  surgery rather than surveillance.
  The STAD framework should note CDH1 germline
  in the clinical context.
```

---

## SECTION V — MSI SUBTYPE: THE HYPERMUTATED ATTRACTOR

```
FREQUENCY:     ~20% of TCGA-STAD
               BUT: ~20% at early stage, dropping
               to ~5% in late-stage disease.
               (MSI gastric cancer has better
               prognosis — patients are diagnosed
               earlier or survive longer, so the
               late-stage denominator is enriched
               for worse subtypes.)
               Older patients (mean age ~68 years)
               Female enrichment.
               Antrum location predominant.

PROGNOSIS:     FAVOURABLE within stage.
               5-year OS: ~50–60% for resectable.
               Excellent immunotherapy response.

CELL OF ORIGIN:
  Gastric antral epithelial cells (foveolar or
  intermediate progenitor) following the
  Correa cascade with MLH1 promoter silencing
  as the key metaplasia → cancer transition event.
  MSI gastric cancer follows the Correa cascade
  more reliably than GS — it is an intestinal-
  type Lauren cancer (columnar gland-forming)
  arising through intestinal metaplasia.

DEFINING MOLECULAR EVENTS:
  MLH1 promoter methylation: ~87% of MSI GC
    The primary mechanism.
    MLH1 encodes a mismatch repair (MMR) protein.
    Methylation silences MLH1 → MMR deficient →
    microsatellite instability (insertion-deletion
    mutations at short tandem repeat sequences
    accumulate throughout the genome).
    Cross-repository connection:
    MLH1 methylation in MSI GC is the SAME
    epigenetic silencing mechanism as in:
    CRC MSS-like tumours with MLH1 promoter
    methylation → CIMP CRC MSI
    Endometrial MSI cancer
    These are the same mechanism: CpG island
    methylation of a mismatch repair gene
    driving a hypermutator phenotype.
    The metabolic/epigenetic connection is
    confirmed: MLH1 CIMP in GC, CRC, and
    endometrial cancer are all products of
    the same aberrant DNA methylation machinery.
  High mutation burden (TMB-H):
    MSI GC has the highest TMB of STAD subtypes.
    ~1000–3000 somatic mutations per tumour
    (vs. ~100–400 in CIN).
    High TMB → high neoantigen load →
    high T cell infiltration (CD8+) →
    immune recognition → immunotherapy
    sensitivity.
  RTK mutations (not amplifications):
    Various RTK mutations occur at high rate
    (ERBB3, FGFR2, MET, EGFR) — but as
    MUTATIONS not amplifications.
    These are hypermutation by-products —
    they occur in these genes because they
    contain microsatellite-like sequences.
    Clinically: these are NOT driving the
    cancer — the high mutation burden is.
    Immunotherapy is the correct treatment,
    not RTK targeting.
  ARID1A mutation: ~44% of MSI GC
    Chromatin remodelling gene.
    Part of the SWI/SNF complex.
    ARID1A mutations are very common in MSI
    cancers (GC, CRC, endometrial, ovarian
    clear cell) — again, hypermutation by-product.

IMMUNE ARCHITECTURE:
  MSI GC is the most immune-infiltrated STAD
  subtype:
    CD8+ T cells:   HIGH — tumour-infiltrating
    CD4+ T cells:   HIGH
    B cells:        Present (lymphoid aggregates)
    PD-L1 CPS:      Often HIGH (>10)
    Interferon-γ signature: HIGH
  This is the gastric equivalent of CRC CMS1
  (the MSI-H, immune-hot, right-sided colorectal
  cancer that responds to pembrolizumab).
  SAME MECHANISM: MLH1 silencing → MSI →
  TMB-H → neoantigen → CD8+ T cell
  infiltration → PD-L1 upregulation →
  immune exhaustion → immunotherapy
  sensitivity after checkpoint release.

TREATMENT:
  Resectable:
    Perioperative FLOT / FOLFOX standard
    BUT: MSI patients may NOT benefit from
    5-FU-based adjuvant chemotherapy
    (same finding as in CRC MSI — 5-FU
    resistance in MSI).
    This is an important clinical subtlety.
    The molecular mechanism: 5-FU targets
    thymidylate synthase and creates DNA
    lesions that MMR normally repairs.
    MMR-deficient cells cannot repair these
    lesions but also trigger cell death differently.
    The clinical net result in colon MSI
    is 5-FU resistance; gastric data is mixed
    but trending in the same direction.
  Advanced/metastatic:
    Nivolumab + FOLFOX/CAPOX (CheckMate-649):
    CPS ≥5 — approved, MSI GC enriched
    for high CPS → largest benefit here.
    Pembrolizumab:
    MSI-H agnostic approval (KEYNOTE-158) —
    covers MSI gastric cancer regardless of
    CPS.
    MSI-H gastric cancer response rate
    to pembrolizumab: ~57% (vs. ~12% in MSS).
    This is among the highest response rates
    seen in any gastrointestinal MSI cancer.
    Long-term survivors with checkpoint
    blockade exist in this subtype.
```

---

## SECTION VI — EBV SUBTYPE: THE EPIGENETICALLY LOCKED ATTRACTOR

```
FREQUENCY:     ~9% of TCGA-STAD
               Small but molecularly unique.
               Predominantly male (higher
               male:female ratio than other
               subtypes).
               Predominantly proximal stomach
               (cardia, gastric body).
               EBV infects the cancer cells —
               the virus is clonally maintained
               in the tumour.

PROGNOSIS:     BEST of the four subtypes.
               5-year OS: ~55–65% for resectable.
               Excellent immunotherapy response.
               The best-prognosis STAD subtype —
               partly because it is detected
               earlier (proximal location)
               and partly because of immune
               infiltration.

THE VIRUS AND THE CANCER:
  Epstein-Barr virus (EBV) infects B cells
  normally — EBV is the cause of infectious
  mononucleosis and is associated with Burkitt
  lymphoma, Hodgkin lymphoma, and nasopharyngeal
  carcinoma.
  In gastric cancer, EBV infects GASTRIC
  EPITHELIAL CELLS and establishes latent
  infection.
  In EBV+ gastric cancer:
    EBV latency type II — expresses LMP1, LMP2A,
    EBER1/2 (small non-coding RNAs)
    EBV-encoded miRNAs downregulate immune
    surveillance genes
    BUT: simultaneously the viral antigens drive
    a high CD8+ T cell infiltrate — the immune
    system recognises the virus even if the cancer
    exploits it.

DEFINING MOLECULAR EVENTS:
  EBV infection:          Clonal (all cells positive)
                          Confirmed by EBER ISH in path.
  PIK3CA mutation:        ~80% — highest of any STAD
                          subtype
                          Many are helical domain
                          mutations (E545K, E542K)
                          Not in kinase domain as in
                          other cancers
                          PI3K → AKT pathway constitutive
                          → cell survival
  DNA hypermethylation:   EXTREME — the most methylated
  (CpG island methylator  STAD subtype
  phenotype — CIMP):      More methylation than MSI GC
                          Silences multiple tumour
                          suppressor genes:
                          CDKN2A (p16)
                          MLH1 (in some cases)
                          PTEN
                          Multiple others
  PD-L1 / PD-L2 amplif:  ~15% of EBV subtype
  (9p24.1 amplification): JAK2 co-amplified
                          High PD-L1 surface expression
                          → PD-1 axis constitutively
                          active → deep T cell
                          exhaustion → exceptional
                          response to PD-1 blockade
  JAK2 amplification:     ~15% — JAK/STAT signalling
                          constitutively active →
                          STAT3 → immune evasion
                          AND JAK inhibitor potential
  ARID1A mutation:        ~55% — highest of any STAD
                          subtype
  PIK3R1 mutation:        ~15%

EPIGENETIC LOCK — CROSS-REPOSITORY CONNECTION:
  The EBV CIMP in gastric cancer is the most
  extreme DNA methylation phenotype in the
  STAD repository.
  It is mechanistically distinct from:
    MSI CIMP: MLH1 methylation → MMR → TMB
    (the methylation drives hypermutation)
    EBV CIMP: EBV-induced methylation →
    tumour suppressor gene silencing →
    BUT: not hypermutation (TP53 mutation
    is VERY LOW in EBV GC — ~6%)
  The EBV CIMP is more analogous to:
    IDH-mutant glioma CIMP (2-HG → TET2
    inhibition → global methylation) in
    that it produces a deeply methylated
    genome WITHOUT necessarily driving
    high mutation burden.
  Two different CIMP mechanisms:
    EBV-driven: virus directly recruits
    DNMT3A/3B to CpG islands via viral
    protein interactions
    IDH-driven: oncometabolite inhibits
    TET enzymes
  Same endpoint: hypermethylated genome
  with silenced differentiation/tumour
  suppressor genes.
  Different clinical consequence:
    IDH CIMP → block in differentiation
    (no immune response)
    EBV CIMP → high immune infiltration
    (because EBV antigens drive T cells)
  This structural distinction will be
  stated in STAD-S4a (EBV before-doc).

TREATMENT:
  Surgery (resectable):
    Same as other STAD — FLOT perioperative.
    EBV status does not currently alter
    surgical approach.
  Advanced/metastatic:
    Pembrolizumab: EXCELLENT response
    The highest PD-L1 CPS scores of any STAD
    subtype (PD-L1 amplification → CPS often
    >30) predicts best single-agent checkpoint
    response in gastric cancer.
    Nivolumab + chemotherapy: also effective.
    Potential: JAK2 inhibitors (ruxolitinib,
    baricitinib) for JAK2-amplified EBV GC.
    Potential: PI3K/AKT inhibitors (PIK3CA
    mutation rate 80% — highest in all STAD)
    → ipatasertib, capivasertib.
    Phase II trials are ongoing for
    AKT inhibitors in PIK3CA-mutant gastric
    cancer, enriched in EBV subtype.
```

---

## SECTION VII — ZOLBETUXIMAB AND CLDN18.2

```
Zolbetuximab (anti-CLDN18.2 antibody) was
approved in 2024 as first-line treatment for
CLDN18.2+ HER2-negative gastric/GEJ cancer.

This is relevant to ALL four subtypes but
has differential implications by subtype:

CLDN18.2 BIOLOGY:
  CLDN18 (claudin-18) is a tight junction
  protein normally expressed in the STOMACH
  and LUNG.
  The stomach-specific isoform is CLDN18.2
  (differs from CLDN18.1 by 8 amino acids
  in the first extracellular domain).
  CLDN18.2 expression in normal tissue:
    Stomach: HIGH — parietal cells and
             neck cells of the gastric gland
    Other organs: ABSENT (except lung CLDN18.1)
  This makes CLDN18.2 a nearly gastric-specific
  surface antigen — and therefore a relatively
  safe target for antibody therapy.

CLDN18.2 IN GASTRIC CANCER:
  Retained in: ~38% of gastric cancers overall
               ~50% of GS subtype
               ~40% of CIN subtype
               ~35% of MSI subtype
               ~30% of EBV subtype
  Clinical eligibility criterion:
    CLDN18.2 IHC ≥2+ in ≥75% of tumour cells

THE ZOLBETUXIMAB MECHANISM:
  Anti-CLDN18.2 IgG1 antibody.
  Binds CLDN18.2 on gastric cancer cells.
  Kills via:
    ADCC (antibody-dependent cellular cytotoxicity)
    CDC (complement-dependent cytotoxicity)
  NOT an ADC — no cytotoxic payload.
  Pure antibody mechanism.
  This means efficacy depends entirely on
  CLDN18.2 expression level.
  If CLDN18.2 falls (depth-negative):
    Deeper tumours express less CLDN18.2
    → less cell-surface target
    → less ADCC/CDC
    → less efficacy

SPOTLIGHT trial (GS and CIN enriched):
  Zolbetuximab + mFOLFOX6 vs FOLFOX6
  PFS: 10.6 vs 8.7 months (HR 0.75, p=0.0066)
  OS: 18.2 vs 15.5 months (HR 0.75, p=0.0053)
  Approved by FDA April 2024 for CLDN18.2+/
  HER2-negative advanced gastric cancer.

GLOW trial:
  Zolbetuximab + CAPOX vs CAPOX
  Similar magnitude of benefit.
  Confirmed the approval.

THE DEPTH SCORE HYPOTHESIS (stated, not predicted):
  If CLDN18.2 is a depth-NEGATIVE gene (falls
  as depth increases — as gastric identity is
  progressively lost), then:
    Shallow depth → high CLDN18.2 → best
    zolbetuximab candidates
    Deep depth → low CLDN18.2 → poor
    zolbetuximab candidates — need alternative
    This would be a clinically actionable
    depth-stratified drug map for zolbetuximab.
  This is NOT stated as a prediction here
  (that belongs in the before-documents).
  It is stated as a structural hypothesis
  to motivate the design of STAD-S1a and
  STAD-S2a.
```

---

## SECTION VIII — THE EXISTING STAD ANALYSIS IN CONTEXT

```
THE COMPLETED ANALYSIS (OrganismCore Documents
70–84 lineage) ran on a combined STAD dataset
(GSE32571 or equivalent bulk RNA-seq/microarray).

WHAT SUBTYPES WERE IN THAT DATASET:
  TCGA-STAD approximate composition:
    CIN:  ~50% (dominant)
    GS:   ~20%
    MSI:  ~20%
    EBV:  ~9%

WHAT THE ANALYSIS FOUND AND WHAT SUBTYPE IT REFLECTS:

  FINDING 1: ZEB2 r=+0.56 with depth.
    SUBTYPE: Predominantly CIN.
    ZEB2 in CIN reflects partial EMT-like
    transcriptional shift at the invasive
    front — consistent with incomplete
    mesenchymal transition in intestinal-
    type gastric cancer.
    WNT5A r=+0.56 (also found): WNT5A is
    a non-canonical WNT ligand associated
    with invasion and aggressiveness —
    confirmed in the 2006 Cancer Research
    paper. This is a CIN-enriched signal
    (intestinal-type gastric cancer shows
    WNT5A elevation at the invasive front).

  FINDING 2: ZEB2/AURKA r=0.9871.
    SUBTYPE: Predominantly CIN.
    AURKA drives centrosome amplification
    and chromosomal instability — the
    defining feature of the CIN subtype.
    This co-expression circuit is a CIN
    signature. The alisertib (Aurora A
    inhibitor) prediction is SPECIFICALLY
    VALID FOR THE CIN SUBTYPE.

  FINDING 3: CDH2/VIM SUPPRESSED
    (ANALYST ASSUMPTION ERROR).
    CORRECT INTERPRETATION (provided here):
    The bulk STAD signal is dominated by CIN
    (not GS). CIN tumours do NOT show full
    EMT (CDH2/VIM elevation). They show PARTIAL
    ZEB2-driven transition. GS tumours (which
    have lost CDH1 but not gained CDH2/VIM)
    are a minority in the bulk signal.
    The assumption error was correctly identified
    and recorded. Its molecular explanation is
    now established.

  FINDING 4: EZH2 SUPPRESSED (not elevated).
    SUBTYPE: The CIN-dominant bulk signal.
    EZH2 suppression in the CIN subtype is
    consistent with the EZH2 context-dependence
    pattern (Pattern 4 in the OrganismCore Cancer
    Framework):
    In BRCA, PAAD, PRAD: EZH2 elevated.
    In STAD bulk: EZH2 suppressed.
    The EBV subtype likely has EZH2 ELEVATED
    (CIMP requires polycomb-mediated silencing).
    The GS subtype may have EZH2 VARIABLE.
    The CIN and MSI bulk signal: EZH2 suppressed.
    THIS REQUIRES SUBTYPE STRATIFICATION to resolve.
    The single most important EZH2 question
    for STAD subtype analysis:
    Is EZH2 elevated specifically in EBV subtype
    (where CIMP is most extreme) but suppressed
    or neutral in CIN, GS, and MSI?
    If yes — EZH2 inhibitors are relevant ONLY
    in EBV subtype gastric cancer. Not CIN.
    Not GS. Not MSI.
    This structural hypothesis will be tested
    in STAD-S4a (EBV before-doc).

  FINDING 5: BROKEN DIFFERENTIATION CIRCUIT.
    The existing analysis found that the
    MASTER DIFFERENTIATION TF CIRCUIT is not
    restoreable in bulk STAD (unlike PAAD and
    PRAD where restoring the switch gene could
    execute the differentiation programme —
    Circuit Integrity Pattern 5).
    The correct interpretation:
    For CIN and GS — the differentiation
    circuit is indeed broken. The cancer cell
    has moved so far from normal gastric
    identity that restoring a single switch
    gene cannot execute the programme.
    ATTRACTOR DISSOLUTION (rather than circuit
    restoration) is the correct strategy for
    CIN and GS — consistent with AURKA/alisertib
    (dissolving the chromosomal instability
    programme) and zolbetuximab (targeting
    a residual surface identity marker).
    For MSI and EBV — the circuit may be MORE
    intact because these tumours are less
    chromosomally disrupted. RESTORING MIST1
    or TFF1 expression via epigenetic targets
    may be more feasible in MSI and EBV than
    in CIN or GS.
    This is stated as a structural hypothesis —
    to be tested in STAD-S3a (MSI) and
    STAD-S4a (EBV).
```

---

## SECTION IX — DATA AVAILABILITY SUMMARY

```
Dataset           Accession         n (tumour)  Normal  Subtype  Power
                                                (n)     labels
─────────────────────────────────────────────────────────────────────────
TCGA-STAD         GDC/phs000178     ~443         36 adj  YES      HIGH
                  cBioPortal                     normal  (TCGA
                                                         4 sub-
                                                         types)
GSE32571          GEO (existing)    ~300+        SOME   PARTIAL  HIGH
                  [THE DATASET      (bulk        normal  (Lauren/
                  ALREADY RUN]      array)       tissue  FAB but
                                                        not TCGA
                                                        4 subtype)
GSE26253          GEO               ~432         0       PARTIAL  MOD-HIGH
                  (Asian cohort)    (survival)
GSE62254          GEO               ~300         0       YES      HIGH
                  (ACRG cohort)     RNA-seq              (ACRG 4
                                                         subtypes
                                                         — see note)
GSE134520         GEO (scRNA-seq)   ~46K cells   YES     YES      MOD
                  (normal +         from        (biopsies (cell
                  IM + early GC)    13 donors   incl.    type res.)
                                                normal)
GSE184198         GEO (scRNA-seq)   STAD cells  YES adj  YES      MOD
                                    + immune     normal
                                    stromal

IMPORTANT NOTE — THE ACRG CLASSIFICATION:
  The Asian Cancer Research Group (ACRG) 2015
  published an alternative 4-subtype classification:
    MSS/EMT: Corresponds to GS subtype
    MSI:     Corresponds to MSI subtype
    MSS/TP53+: TP53 wildtype CIN-like
    MSS/TP53−: TP53 mutant CIN-like
  The ACRG classification is used in some datasets
  (especially Asian cohorts — ACRG, TCGA has
  non-Asian overrepresentation).
  GSE62254 uses ACRG labels.
  When using GSE62254, MSS/EMT = GS.
  Cross-mapping is required.

NORMAL GASTRIC REFERENCE:
  GSE134520 (scRNA-seq): The best normal gastric
  reference for the cellular hierarchy.
  Contains cells from non-atrophic gastritis
  biopsies (the closest available normal).
  True healthy gastric mucosa without any H.
  pylori or inflammation is rare in public
  datasets (most Asian populations have
  H. pylori exposure by adulthood).

PRIMARY DATASET FOR SUBTYPE ANALYSIS:
  TCGA-STAD + GSE62254 (ACRG) + GSE26253
  provide the bulk RNA-seq power for
  subtype-stratified analysis.
  TCGA-STAD for the 4-subtype framework.
  GSE62254 ACRG for treatment response
  validation (survival curves by subtype).
  GSE134520 for normal gastric cell
  identity markers (the Waddington baseline).
```

---

## SECTION X — PLANNED ANALYSIS ORDER

```
ORDER:

  STAD-S1   CIN subtype       TCGA-STAD          HIGH POWER
            (dominant         CIN subset +
            signal)           GSE62254
                              REASON: The dominant
                              STAD subtype (~50%).
                              The source of the
                              existing ZEB2/AURKA
                              finding. The HER2+
                              actionable population.
                              Confirms alisertib
                              prediction in CIN-
                              specific context.
                              CLDN18.2 depth
                              relationship first
                              tested here.
                              Clinical output:
                              HER2+ vs HER2− drug
                              maps within CIN;
                              depth score as T-DXd
                              candidate predictor.

  STAD-S2   GS subtype        TCGA-STAD          MOD-HIGH
            (hardest          GS subset +
            problem)          GSE62254
                              REASON: Worst
                              prognosis, least
                              actionable, but
                              the molecular
                              explanation for
                              the analyst
                              assumption error.
                              CDH1/RHOA/CLDN18
                              biology.
                              Zolbetuximab depth
                              prediction (CLDN18.2
                              depth-negative
                              hypothesis tested).
                              Novel hypothesis:
                              RHOA/FAK pathway
                              depth-positive in GS
                              → FAK inhibitors
                              (defactinib) predicted
                              for deep GS.
                              Clinical output:
                              zolbetuximab patient
                              selection from depth
                              score.

  STAD-S3   MSI subtype       TCGA-STAD          MOD-HIGH
                              MSI subset +
                              GSE62254 MSI +
                              CheckMate-649 data
                              (if accessible)
                              REASON: The immune
                              attractor. CRC CMS1
                              parallel (MLH1 CIMP
                              as the shared
                              mechanism confirmed).
                              Immunotherapy
                              response prediction.
                              5-FU resistance
                              depth prediction.
                              Clinical output:
                              depth score as
                              pembrolizumab
                              response predictor
                              in MSI GC.

  STAD-S4   EBV subtype       TCGA-STAD          LOW-MOD
                              EBV subset (~40
                              samples — low n)
                              REASON: The most
                              epigenetically
                              distinct subtype.
                              PIK3CA mutation ~80%.
                              EZH2 direction in
                              EBV: elevated or
                              suppressed?
                              PI3K/AKT inhibitor
                              prediction.
                              JAK2 amplification
                              → JAK inhibitor
                              prediction.
                              PD-L1 amplification
                              → pembrolizumab best
                              response in STAD.
                              Low n — findings
                              are HYPOTHESIS-
                              GENERATING only.

  STAD-X    Cross-subtype     After S1–S4.
            + EZH2            QUESTIONS:
            resolution          1. Does depth score
                                   order subtypes
                                   by prognosis?
                                   MSI/EBV (shallowest,
                                   best) → CIN → GS
                                   (deepest, worst)?
                                2. EZH2 direction by
                                   subtype:
                                   EBV: elevated?
                                   CIN: suppressed?
                                   GS: suppressed?
                                   MSI: suppressed?
                                   This resolves the
                                   bulk suppression
                                   finding.
                                3. CLDN18.2/CLDN18
                                   across all four
                                   subtypes — does
                                   it fall uniformly
                                   with depth or only
                                   in specific
                                   subtypes?
                                4. AURKA/ZEB2 co-
                                   expression: is it
                                   CIN-specific or
                                   present in all
                                   subtypes?
                                5. MUC5AC/TFF1 as
                                   universal depth-
                                   negative: confirmed
                                   across all four?
```

---

## SECTION XI — WHAT THIS DOCUMENT DOES NOT CONTAIN

```
This document contains:
  ✓ Normal gastric gland hierarchy (foveolar →
    isthmus → neck → chief/parietal)
  ✓ The H. pylori → SPEM → IM → cancer
    Correa cascade as the Waddington map
  ✓ All four TCGA subtypes (CIN, GS, MSI, EBV)
    with frequencies, cells of origin, molecular
    events, and treatment
  ✓ The existing STAD analysis findings placed
    in subtype context (ZEB2/AURKA = CIN signal;
    CDH2/VIM suppression = explained by GS
    biology; EZH2 suppression = bulk CIN signal)
  ✓ The analyst assumption error molecular
    explanation (GS is NOT classic EMT)
  ✓ Zolbetuximab mechanism and CLDN18.2 biology
  ✓ The CRC cross-repository connections
    (MSI GC = CMS1 parallel; GS = CMS4 parallel;
    CIN = CMS2 parallel)
  ✓ The IDH-glioma/EBV CIMP structural comparison
  ✓ Data availability (TCGA, ACRG/GSE62254,
    GSE32571 existing, GSE134520 scRNA-seq)
  ✓ The ACRG classification cross-mapping note

This document does NOT contain:
  ✗ Depth score predictions
  ✗ Specific false attractor gene predictions
  ✗ Drug target predictions
  ✗ Epigenetic lock mechanism hypotheses
  ✗ Cross-subtype structural predictions beyond
    what is established in the literature

All of the above belong in the BEFORE documents.
STAD-S1a (CIN before-document) is next.
Written before any script runs.
Before any data loads.
```

---

## STATUS BLOCK

```
document:           STAD_Subtype_Orientation.md
folder:             Cancer_Research/STAD/Subtypes/
status:             COMPLETE — ORIENTATION ONLY
date:               2026-03-03
author:             Eric Robert Lawson
                    OrganismCore

subtypes_covered:
  CIN (~50%): TP53, RTK amp, chromosomal  [1 of 4]
              instability, HER2, ZEB2/AURKA
  GS (~20%):  CDH1/RHOA/CLDN18-ARHGAP,   [2 of 4]
              signet ring, broken cohesion
  MSI (~20%): MLH1 CIMP, TMB-H, CD8+     [3 of 4]
              infiltrate, pembrolizumab
  EBV (~9%):  Viral CIMP, PIK3CA, PD-L1  [4 of 4]
              amp, JAK2 amp, best prognosis

analyses_started:   0 (new subtype series)
existing_analysis:  Cancer_Research/STAD/ — complete
                    (OrganismCore Documents 70–84 lineage)
                    ZEB2/AURKA finding = CIN signal
                    CDH2/VIM suppression = analyst
                    assumption error EXPLAINED by GS
                    biology (GS ≠ classic EMT)
                    EZH2 suppression = bulk CIN signal
                    Broken differentiation circuit =
                    correct for CIN and GS; may be
                    intact for MSI and EBV

next_document:      STAD-S1a
                    CIN Subtype Before-Document
                    (predictions locked before
                    TCGA-STAD CIN subset loads)

critical_note_1:    The analyst assumption error
                    (CDH2/VIM suppressed, not elevated)
                    is now fully explained:
                    GS subtype (the only STAD subtype
                    with CDH1 loss) does NOT gain
                    CDH2/VIM — it is a LOSS OF
                    COHESION cancer, not a gain of
                    mesenchymal identity cancer.
                    The signet ring cell is CDH1-low,
                    CDH2-LOW, VIM-LOW, MUC5AC-variable.
                    This distinguishes GS from classic
                    EMT and from CMS4 CRC (which DOES
                    gain CDH2/VIM/ZEB2 in the mesenchymal
                    EMT programme).
                    The framework correctly identified
                    and documented the error. The
                    explanation was not available in
                    the original bulk analysis — it
                    required subtype stratification.

critical_note_2:    The CRC structural parallels are
                    the strongest cross-cancer structural
                    connections in the repository:
                    STAD MSI  ↔  CRC CMS1 (MLH1 CIMP
                                → TMB-H → immune hot
                                → checkpoint response)
                    STAD CIN  ↔  CRC CMS2 (chromosomal
                                instability → RTK amp
                                → intermediate prognosis
                                → chemotherapy response)
                    STAD GS   ↔  CRC CMS4 (stromal/
                                mesenchymal → worst
                                prognosis → immune
                                excluded → resistant)
                    These are not coincidences — they
                    reflect the same Waddington
                    attractor landscape operating in
                    two different GI epithelial tissues.
                    Different tissues. Same geometry.

critical_note_3:    The EZH2 direction in STAD requires
                    subtype resolution. The bulk
                    suppression finding may reflect:
                    a) CIN-dominant signal where EZH2
                       is genuinely suppressed (unlike
                       BRCA/PAAD/PRAD)
                    b) EBV subtype where EZH2 is
                       ELEVATED (CIMP requires PRC2)
                       being masked by the dominant
                       CIN signal
                    Resolution in STAD-S4a (EBV) will
                    determine whether EZH2 inhibitors
                    (tazemetostat) are relevant for a
                    minority EBV subtype patient
                    selection tool.

critical_note_4:    CLDN18 / CLDN18.2 is the most
                    clinically actionable depth-negative
                    marker hypothesis in STAD.
                    CLDN18.2 is a NORMAL PARIETAL CELL
                    MARKER — it is expressed by
                    terminally differentiated gastric
                    cells. Its retention in cancer
                    marks SHALLOW depth (near-normal
                    gastric identity). Its loss marks
                    DEEP depth (identity-lost cancer).
                    If the depth score can predict
                    CLDN18.2 expression level, the
                    depth score becomes a ZOLBETUXIMAB
                    PATIENT SELECTION TOOL — the most
                    directly clinically translatable
                    output in the STAD subtype series.
                    This is the primary clinical output
                    target for STAD-S1a and STAD-S2a.
```
