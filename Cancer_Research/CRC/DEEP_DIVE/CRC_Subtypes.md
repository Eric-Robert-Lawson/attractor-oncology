# COLORECTAL CANCER — SUBTYPE ORIENTATION DOCUMENT
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

  A complete map of the colorectal cancer subtype
  landscape — what the four Consensus Molecular Subtypes
  (CMS1–4) are, where they come from anatomically and
  cellularly, why they represent genuinely distinct
  Waddington landscapes rather than superficial
  expression variation, what public data exists to
  analyze each one, and what the existing analysis
  in this repository captured.

CRC has the most rigorously validated molecular subtype
classification system of any solid tumour in this
repository. The CMS classification (Guinney et al.,
Nature Medicine 2015) was derived from six independent
datasets across 3,000+ patients and has been
independently replicated in dozens of cohorts since.

It is not a soft classification.
It is the closest thing to an internationally accepted
standard for CRC molecular stratification that exists
in solid tumour oncology.

The four CMS subtypes are not statistical clusters.
They are biologically distinct disease states:
  CMS1: immune-driven, MSI-high, right-sided, BRAF
  CMS2: WNT/MYC-driven, canonical, left-sided, CIN
  CMS3: metabolic, KRAS-driven, mixed MSI
  CMS4: mesenchymal, TGF-β-driven, worst prognosis

The existing analysis in Cancer_Research/CRC/ ran on
the combined TCGA-COAD/READ dataset without CMS
stratification. It captured the dominant variance
across the full dataset — predominantly CMS2 (the
largest subtype, ~37%) with contributions from all.

The subtype analyses in this folder will decompose that.
```

---

## DOCUMENT METADATA

```
document_id:        CRC_Subtype_Orientation
series:             CRC (Colorectal Cancer — Subtypes)
folder:             Cancer_Research/CRC/Subtypes/
date:               2026-03-03
author:             Eric Robert Lawson
status:             ORIENTATION ONLY — no predictions,
                    no analysis, no epigenetic hypotheses
next_document:      CRC_CMS2_before.md
                    (Document CRC-S1a — CMS2 before-doc)
protocol_version:   Workflow_Protocol.md v2.0
```

---

## SECTION I — THE NORMAL COLORECTAL EPITHELIAL HIERARCHY

```
Before defining the four false attractors, the normal
colorectal epithelial hierarchy must be understood.

THE CRYPT — THE BASIC UNIT OF INTESTINAL IDENTITY:

  The colon and rectum are lined by a single layer of
  columnar epithelium organised into finger-like
  invaginations called crypts. The crypt is a
  self-renewing unit with a strict spatial hierarchy:

  CRYPT BASE (bottom):
    Intestinal stem cells (ISCs):
      LGR5+, OLFM4+, ASCL2+
      The continuously cycling cell population that
      drives all epithelial renewal.
      Every intestinal epithelial cell traces its
      lineage to an LGR5+ ISC.
      WNT signalling is highest here — required for
      ISC maintenance and proliferation.
      LGR5 itself is a WNT target gene.

    Paneth cells (small intestine only):
      LYZ+, DEFA5+
      Support ISC niche via EGF, WNT3, and Notch
      ligands.
      In the colon: "deep secretory cells" play an
      analogous niche role.

  TRANSIT AMPLIFYING ZONE (mid-crypt):
    Daughter cells of ISCs that amplify cell number
    before committing to a differentiated fate.
    Still proliferating (Ki67+).
    WNT signalling declining.
    CDX2 expression rising.
    Beginning to acquire absorptive or secretory fate.

  DIFFERENTIATED SURFACE (upper crypt and surface):

    ABSORPTIVE lineage — Colonocytes:
      CDX2+, KRT20+, CEACAM7+, SLC26A3+
      Absorb water and electrolytes.
      Non-proliferating. Shed every 3–5 days.
      The most abundant cell type (~80% of colon).

    SECRETORY lineages — three cell types:
      Goblet cells:  MUC2+, TFF3+, AGR2+
                     Secrete mucus layer — protection
                     of epithelial surface
      Enteroendocrine cells: CHGA+, NEUROG3+
                     Hormone-secreting cells
      Tuft cells:    DCLK1+, POU2F3+
                     Rare chemosensory cells

  THE WADDINGTON LANDSCAPE OF THE NORMAL CRYPT:
    The crypt IS a Waddington landscape in miniature.
    LGR5+ ISC at the base = the most stem-like state
    (capable of all fates, high WNT, high proliferation)
    Colonocyte at the surface = the terminal
    differentiated state (no WNT, no proliferation,
    absorptive function)
    The normal cell moves UPWARD through the crypt:
    ISC → transit amplifying → colonocyte → shed

    CRC REVERSAL OF DIRECTION:
    In CRC, the cancer cell moves DOWNWARD — back
    toward the ISC state or BEYOND the ISC into a
    fully de-anchored proliferative state.
    APC/WNT pathway activation is the molecular engine
    that drives this reversal.
    In normal tissue: WNT is HIGH at the base, LOW
    at the surface.
    In CRC: WNT is constitutively HIGH everywhere —
    the cancer cell is perpetually stuck at the base
    of the crypt or below it.

  CELL OF ORIGIN FOR ALL CRC SUBTYPES:
    The LGR5+ intestinal stem cell is the cell of
    origin for the majority of CRC cases.
    Evidence:
      - Mouse models: APC loss in LGR5+ cells rapidly
        produces adenomas; APC loss in differentiated
        colonocytes produces fewer, smaller tumours
      - Human CRC genetic analysis: the founder
        mutation (APC) occurs in cells with stem-like
        epigenetic marks
      - Single-cell sequencing: CRC cells map to
        a range of crypt positions but origin traces
        to ISC-like progenitors
    EXCEPTION: CMS4 mesenchymal tumours may have a
    more complex origin involving epithelial-to-
    mesenchymal transition (EMT) from progenitor cells,
    and the stromal compartment contributes heavily
    to their transcriptional signature.

  THE NORMAL IDENTITY PROGRAMME THAT IS SUPPRESSED
  IN ALL CRC SUBTYPES:
    CDX2:       Master intestinal TF — its expression
                is reduced as CRC deepens. CDX2 loss
                is a key marker of dedifferentiation.
    KRT20:      Colonocyte differentiation marker
    SLC26A3:    Chloride transporter — colonocyte
                function marker
    CEACAM7:    Surface differentiation marker
    MUC2:       Goblet cell identity (secretory lineage)
    HOXB13:     Colonic identity TF
    These are the genes expected to fall with depth
    in ALL four CMS subtypes. Their loss defines the
    depth axis from the normal crypt position.
```

---

## SECTION II — THE CONSENSUS MOLECULAR SUBTYPES: OVERVIEW

```
The CMS classification (Guinney et al., Nat Med 2015)
was the culmination of six independent research groups
each having produced their own CRC subtype systems.
The CMS framework emerged from a meta-analysis
integrating all six systems across 4,151 patients.

What the CMS system captures:
  Not just gene expression variation.
  Distinct biological programmes:
    CMS1: An immune activation programme
    CMS2: A WNT/MYC proliferative programme
    CMS3: A metabolic reprogramming programme
    CMS4: A stromal invasion programme

  These programmes are not points on a continuum.
  They are qualitatively different false attractor
  states — different molecular engines driving
  the cancer cell's persistence in a non-normal
  proliferative state.

CMS PREVALENCE (from the 2015 paper, n=3,816):
  CMS1: 14%  (MSI-immune)
  CMS2: 37%  (canonical)
  CMS3: 13%  (metabolic)
  CMS4: 23%  (mesenchymal)
  Mixed/unclassified: 13%

  CMS2 is the largest — it is what most people think
  of when they think of "colorectal cancer."
  CMS4 has the worst prognosis.
  CMS1 has a paradox: good initial prognosis but
  poor survival after metastasis.
  CMS3 is the most metabolically distinct.

GEOGRAPHIC DISTRIBUTION NOTE:
  Left-sided colon and rectum: CMS2 dominant
  Right-sided colon:           CMS1 and CMS3 enriched
  This anatomical split is not incidental —
  the right and left colon have different embryological
  origins, different microbiomes, different cell
  populations, and different molecular backgrounds.
  Right-sided CRC arises in the midgut-derived colon.
  Left-sided CRC arises in the hindgut-derived colon.
  The CMS subtype distribution reflects this embryological
  split. This is a Waddington consequence — the starting
  cell position in the landscape differs between right
  and left colon.
```

---

## SECTION III — CMS2: THE CANONICAL SUBTYPE

```
PREVALENCE: 37% — THE DOMINANT FORM

CLINICAL FACTS:
  Location:         Predominantly left-sided colon
                    and rectum
  Histology:        Well to moderately differentiated
                    adenocarcinoma. Glandular structures
                    retained. CDX2 usually maintained.
  MSI status:       Microsatellite STABLE (MSS) — the
                    mismatch repair system is intact
  Prognosis:        Intermediate — better than CMS4,
                    worse than CMS1 at early stages
                    Relapse-free survival: moderate
  Standard therapy: 5-FU/capecitabine + oxaliplatin
                    (FOLFOX / CAPOX) — standard adjuvant
                    Anti-EGFR therapy (cetuximab,
                    panitumumab) in RAS-wildtype left-sided
                    CRC — most effective in CMS2.
                    Bevacizumab (anti-VEGF).
  Response to immunotherapy: LOW — MSS tumours do not
                    respond to PD-1 inhibitors.
                    PD-1/PD-L1 blockade is NOT effective
                    in CMS2 (or CMS3/CMS4) as a single
                    agent based on current trial data.

CELL OF ORIGIN:
  LGR5+ intestinal stem cell (ISC).
  CMS2 most closely resembles the transit amplifying
  zone of the crypt — high proliferation, WNT active,
  MYC high, CDX2 still expressed.
  The false attractor is a proliferative ISC-like state
  that has blocked terminal differentiation but retained
  enough colonocyte identity to maintain glandular
  architecture.

DEFINING MOLECULAR EVENTS:
  APC mutation:     ~85% of CMS2
                    The founding event in CMS2.
                    APC is the destruction complex
                    component that normally degrades
                    β-catenin. APC mutation → β-catenin
                    accumulates → WNT target genes
                    constitutively active.
                    WNT target genes include:
                      MYC (the key downstream effector)
                      CCND1 (cyclin D1 — cell cycle)
                      AXIN2 (feedback regulator)
                      LGR5 (ISC marker — retained)
                      CD44 (stem-like marker)
  KRAS mutation:    ~30% (lower than CMS3)
                    Secondary event — accelerates
                    proliferation but not the primary
                    driver.
  TP53 mutation:    ~60%
                    Required for adenoma-to-carcinoma
                    transition.
  Chromosomal instability (CIN):
                    The dominant genomic mechanism.
                    Whole-arm copy number alterations:
                      Chr18q loss (SMAD4, DCC)
                      Chr17p loss (TP53)
                      Chr8q gain (MYC)
                      Chr20q gain
                    These copy number events amplify the
                    WNT/MYC signal further.
  SMAD4 loss:       ~15–20% (chr18q)
                    TGF-β tumour suppressor — loss
                    late event, promotes invasion.
  No MSI:           Mismatch repair is intact.
  No BRAF V600E:    Present in CMS1, absent in CMS2.

CORE EXPRESSION PROGRAMME OF THE FALSE ATTRACTOR:
  UP in CMS2 deep stratum:
    MYC and MYC targets (the dominant programme)
    CCND1 (cyclin D1)
    CDK4, CDK6
    E2F targets (proliferation genes)
    EGFR pathway components
    LGR5 (ISC marker — retained)
  DOWN in CMS2 deep stratum:
    CDX2 (intestinal master TF — progressive loss)
    KRT20 (colonocyte marker)
    CEACAM7
    SLC26A3
    Goblet cell markers (MUC2, TFF3)

WADDINGTON STRUCTURE:
  CMS2 is the canonical false attractor that most
  closely resembles the textbook picture of CRC:
  a cell stuck at the ISC/transit amplifying position
  due to constitutive WNT/MYC activation.
  The depth axis runs from:
    Normal colonocyte (CDX2+, KRT20+, non-proliferating)
    → CMS2 shallow (LGR5+, MYC+, glandular)
    → CMS2 deep (MYC high, CDX2 lost, CIN)
  This is a classical gain-of-proliferation, loss-of-
  differentiation false attractor.

AVAILABLE PUBLIC DATA:
  TCGA-COAD/READ:   CMS labels available via published
                    classifier. CMS2: ~37% of ~600 total
                    samples = ~220 CMS2 samples.
  GEO datasets:
    GSE14333:       n=290, CMS labels in supplementary
    GSE39582:       n=585 (largest published CRC cohort
                    with CMS labels — French PETACC-3
                    cohort), survival annotated. CRITICAL.
    GSE17536:       n=177, survival annotated
    GSE33113:       n=90, treatment response data
  Normal colon reference:
    GTEx colon:     Transverse + sigmoid colon
                    Multiple tissue samples.
    GSE4183:        Normal colon mucosa n=8 + adenoma
    TCGA adjacent:  TCGA-COAD normal adjacent tissue

ANALYSIS STATUS:    NOT YET STARTED
PLANNED DOCUMENT:   CRC-S1a (CMS2 before-document)
PLANNED DATASET:    GSE39582 (n=585, primary large cohort)
                    + TCGA-COAD/READ CMS2 subset
PRIORITY:           FIRST — CMS2 is the largest subtype
                    and the primary driver of the existing
                    combined analysis. The WNT/MYC axis
                    is the canonical CRC programme.
                    Anti-EGFR response is best in CMS2
                    (RAS-wildtype, left-sided) — the
                    clinical utility of depth stratification
                    in CMS2 is highest.
```

---

## SECTION IV — CMS1: THE MSI-IMMUNE SUBTYPE

```
PREVALENCE: 14%

CLINICAL FACTS:
  Location:         Predominantly right-sided colon
  Histology:        Poorly differentiated, mucinous,
                    signet-ring features common.
                    High immune infiltration — tumour-
                    infiltrating lymphocytes (TILs) are
                    a defining histological feature.
  MSI status:       MSI-HIGH (microsatellite instable)
                    Deficient mismatch repair (dMMR)
                    Most commonly: MLH1 promoter
                    hypermethylation (sporadic MSI-H)
                    Or: Lynch syndrome (germline MLH1,
                    MSH2, MSH6, PMS2 mutation — ~3% of
                    all CRC)
  Prognosis:        PARADOXICAL:
                    Stage II/III: GOOD prognosis
                    (better than CMS2/3/4 at same stage)
                    Stage IV (metastatic): POOR prognosis
                    (worse than MSS after metastasis —
                    MSI-H metastatic CRC is more aggressive
                    once it has disseminated)
  Standard therapy: Stage II: No adjuvant 5-FU benefit
                    (paradox — 5-FU actually may harm
                    MSI-H stage II patients)
                    Stage III: FOLFOX still used
                    Stage IV: PEMBROLIZUMAB (PD-1 inhibitor)
                    is NOW FIRST-LINE standard for MSI-H/
                    dMMR metastatic CRC (KEYNOTE-177 trial,
                    2020 — landmark approval).
                    This is the most clinically actionable
                    distinction in CRC:
                    MSI-H CRC → immunotherapy works
                    MSS CRC → immunotherapy does not work
  Response to immunotherapy: HIGH — CMS1 is the only
                    CRC subtype that responds robustly
                    to PD-1 inhibitors.

CELL OF ORIGIN:
  LGR5+ ISC — same as CMS2.
  BUT the carcinogenesis pathway is completely different.

  CMS1 arises via the SERRATED PATHWAY:
    Normal crypt → hyperplastic polyp →
    Sessile serrated lesion (SSL) →
    CMS1 adenocarcinoma
    (vs. CMS2 which arises via the adenoma pathway:
    normal → tubular adenoma → CRC)

  The serrated pathway features:
    BRAF V600E mutation as the initiating event
    (not APC)
    MLH1 promoter hypermethylation (CIMP-H)
    → MSI-high via epigenetic silencing of mismatch
    repair, not genetic mutation
    This means CMS1 starts with an EPIGENETIC event
    (methylation of MLH1), not a genetic mutation.
    The CIMP-H programme broadly methylates the genome —
    silencing multiple tumour suppressors simultaneously.

DEFINING MOLECULAR EVENTS:
  BRAF V600E:       ~50% of CMS1
                    THE defining mutation of the serrated
                    pathway. Constitutively activates
                    MAPK/ERK signalling.
                    CRITICAL: BRAF V600E in CRC is NOT
                    the same as BRAF V600E in melanoma.
                    Vemurafenib (BRAF inhibitor) alone
                    fails in CRC because EGFR feedback
                    reactivates MAPK.
                    ENCORAFENIB + CETUXIMAB (anti-EGFR)
                    is the approved combination for BRAF
                    V600E MSS CRC (BEACON trial) — but
                    most BRAF V600E CRC is MSI-H (CMS1).
  MLH1 hypermethylation:
                    ~75% of sporadic MSI-H CRC
                    The epigenetic silencing of the
                    mismatch repair gene that drives
                    microsatellite instability.
  MSI-H:            By definition. Thousands of insertion/
                    deletion mutations at microsatellite
                    loci. Very high tumour mutation burden
                    (TMB-high). This is why immunotherapy
                    works — the tumour produces many
                    neoantigens that T cells can recognise.
  CIMP-H:           CpG island methylator phenotype —
                    high. Widespread epigenetic silencing.
                    This is the structural connection to
                    IDH-mutant glioma (CIMP also present
                    there via 2-HG mechanism) and to the
                    TCA/αKG axis across the repository.
  APC mutation:     LOW (~15% — much less than CMS2)
                    CMS1 does NOT primarily use the
                    WNT/APC pathway. BRAF drives it instead.
  RAS mutation:     LOW (~5% — KRAS mutations are rare
                    in BRAF V600E CRC because the two
                    mutations are largely mutually exclusive)

IMMUNE ARCHITECTURE:
  CMS1 has the most complex immune architecture of
  any CRC subtype — and one of the most studied immune
  microenvironments of any cancer in this repository.

  CD8+ T cells:     HIGH — cytotoxic T cell infiltration
                    This is why immunotherapy works.
  CD4+ T cells:     HIGH — helper T cell infiltration
  FOXP3+ Tregs:     Present but outnumbered
  NK cells:         Elevated
  B cells:          Elevated
  Macrophages:      M1-polarised (pro-inflammatory)
                    in contrast to CMS4 where M2
                    (immunosuppressive) macrophages
                    dominate.
  PD-L1:            High expression on tumour cells
                    and immune cells.
                    This is the checkpoint that
                    pembrolizumab releases.

  THE STRUCTURAL PARADOX:
    CMS1 has the highest immune infiltration and yet
    the cancer still grows — it is not eliminated by
    the immune system despite high TIL density.
    The explanation: PD-L1 on the tumour cell binds
    PD-1 on the T cell and suppresses the T cell response.
    This is the checkpoint. Pembrolizumab blocks this
    interaction — releasing the T cell to kill.
    The depth score in CMS1 should therefore correlate
    with PD-L1 expression and with T cell exhaustion
    markers (LAG3, TIM3, TIGIT).

WADDINGTON STRUCTURE:
  CMS1 is the most structurally unusual false attractor
  in the CRC landscape.
  The cell of origin (LGR5+ ISC) is the same as CMS2 —
  but the initiating event (BRAF + CIMP, not APC) and
  the false attractor state (immune-active, BRAF-driven,
  highly mutated) are completely different.
  The depth axis in CMS1 runs:
    Normal serrated epithelium → SSL → MSI-H CMS1 →
    deep CMS1 (immune exclusion, T cell exhaustion)
  The deep CMS1 tumour has paradoxically escaped
  immune elimination by upregulating checkpoint
  molecules — not by reducing immune cell density.
  The immune infiltrate is HIGH but EXHAUSTED.
  This is the CMS1 depth axis: not loss of normal
  identity (which is what drives CMS2/3/4 depth) but
  progressive T cell exhaustion and checkpoint
  upregulation despite maintained tumour immunogenicity.

AVAILABLE PUBLIC DATA:
  TCGA-COAD/READ:   CMS1 ~14% of ~600 = ~84 samples.
                    Small n within TCGA alone.
  GSE39582:         Best powered — CMS labels available,
                    n=585 total, ~82 CMS1.
  GSE14333:         n=290, ~40 CMS1
  MSI-specific datasets:
    GSE20916:       MSI-H CRC specific
    GSE13294:       Includes MSI/MSS annotation

ANALYSIS STATUS:    NOT YET STARTED
PLANNED DOCUMENT:   CRC-S2a (CMS1 before-document)
PLANNED DATASET:    GSE39582 CMS1 subset
                    + TCGA-COAD CMS1 subset
                    (combined n ~160–170 for CMS1)
PRIORITY:           SECOND — CMS1 is the only CRC subtype
                    where immunotherapy works. The depth
                    score in CMS1 maps the exhaustion
                    axis. Clinically the most actionable
                    finding: depth score predicts
                    pembrolizumab response in dMMR CRC.
```

---

## SECTION V — CMS4: THE MESENCHYMAL SUBTYPE

```
PREVALENCE: 23% — THE WORST PROGNOSIS SUBTYPE

CLINICAL FACTS:
  Location:         Both right and left colon, rectum
  Histology:        Poorly differentiated to undifferentiated.
                    Stromal desmoplasia prominent.
                    Mucinous features common.
                    Invasive front with EMT morphology.
  MSI status:       Microsatellite STABLE (MSS)
  Prognosis:        WORST of all four subtypes.
                    Highest rate of metastasis.
                    Shortest relapse-free and overall
                    survival.
                    Stage IV CMS4: median OS ~12–14 months
                    (vs. ~24 months in CMS2 at Stage IV)
                    Most patients present at advanced stage.
  Standard therapy: FOLFOX / FOLFIRI + bevacizumab
                    (anti-VEGF — bevacizumab has
                    proportionally MORE benefit in CMS4
                    than in other subtypes due to high
                    VEGF/angiogenesis gene expression)
                    Anti-EGFR therapy: LESS effective
                    in CMS4 (despite RAS status) due to
                    stromal resistance mechanisms.
  Response to immunotherapy: LOW — MSS, and the immune
                    microenvironment is suppressed by
                    TGF-β-driven immune exclusion.

CELL OF ORIGIN:
  The most complex origin of any CRC subtype.
  Two models exist:

  Model 1 (epithelial model):
    LGR5+ ISC undergoes epithelial-to-mesenchymal
    transition (EMT) — driven by TGF-β signalling.
    The cell loses CDX2 and epithelial identity,
    acquires VIM, FN1, CDH2 (N-cadherin), ZEB1.
    The false attractor is a post-EMT mesenchymal state.
    This model predicts: the deepest CMS4 cells will
    have LOST all ISC/epithelial markers and gained
    mesenchymal identity markers.

  Model 2 (stromal contamination model):
    The CMS4 transcriptional signature is partly driven
    by the tumour microenvironment — specifically by
    cancer-associated fibroblasts (CAFs) and macrophages
    that infiltrate the tumour stroma.
    Bulk RNA-seq of CMS4 tumours captures CAF and
    macrophage transcriptomes in addition to the
    tumour cell transcriptome.
    Single-cell RNA-seq has shown that CMS4 tumours
    contain both:
      a. EMT-high tumour cells (true mesenchymal cancer
         cells)
      b. High CAF infiltration (stromal contamination)
    BOTH contribute to the CMS4 bulk signal.
    The framework must address this: the depth score
    in CMS4 will partly capture stromal content, not
    just tumour cell depth.

  Resolution:
    Both models are partly correct.
    CMS4 is the mesenchymal CRC state where:
    - The tumour cells have undergone at least partial EMT
    - The microenvironment is heavily stromal/fibroblastic
    - TGF-β is the dominant signalling axis (both in
      tumour cells and CAFs)
    The depth score captures the combined programme.
    This is noted explicitly so the before-document
    can address it.

DEFINING MOLECULAR EVENTS:
  TGF-β pathway activation:
                    The defining molecular programme.
                    TGFBR1/2, SMAD2/3 signalling active.
                    TGF-β drives:
                      EMT (CDH1 loss, VIM gain)
                      CAF recruitment and activation
                      Immune exclusion (T cell suppression)
                      Angiogenesis
  SMAD4 loss:       ~30% (higher than other subtypes)
                    Paradox: SMAD4 loss should block
                    TGF-β signalling — but in the
                    SMAD4-intact CMS4 tumours, TGF-β
                    drives invasion directly.
                    SMAD4 loss removes the growth
                    inhibitory arm of TGF-β while
                    leaving the pro-invasive arm active.
  KRAS mutation:    ~40%
  APC mutation:     ~40% (lower than CMS2)
  TP53 mutation:    ~60%
  WNT pathway:      Active but not dominant
                    (less APC-driven than CMS2)
  EMT master TFs:   ZEB1, ZEB2, TWIST1, SNAI1
                    These reprogram the ISC toward
                    mesenchymal identity.
  VEGF:             HIGH — angiogenesis programme active.
                    This is why bevacizumab shows
                    proportionally greater benefit.
  MMP2, MMP9, MMP14: HIGH — matrix metalloproteinases
                    drive invasive front remodelling.
  CDH1 loss:        Progressive with depth — E-cadherin
                    loss releases epithelial cohesion.
  VIM, FN1, CDH2:   HIGH — mesenchymal identity markers.

IMMUNE ARCHITECTURE:
  CMS4 has the most immunosuppressed microenvironment
  of all CRC subtypes.
  CD8+ T cells:     LOW (excluded from tumour core)
  FOXP3+ Tregs:     ELEVATED
  M2 macrophages:   HIGH (TGF-β drives M2 polarisation)
  TGF-β levels:     HIGHEST of all subtypes
  TGF-β drives immune exclusion by:
    - Suppressing CD8+ T cell infiltration
    - Converting M1 to M2 macrophages
    - Activating Tregs
    - Reducing NK cell activity
  This immune exclusion explains why immunotherapy
  fails in CMS4 even though the tumour expresses
  some checkpoint molecules.
  The immune cells cannot reach the tumour.
  TGF-β is the wall.
  TGF-β inhibitors combined with immunotherapy are
  in clinical trials specifically for CMS4 CRC — the
  rationale is to remove the TGF-β wall and allow
  immune cell entry.

WADDINGTON STRUCTURE:
  CMS4 is the deepest false attractor in the CRC
  landscape by clinical outcome — worst prognosis,
  highest metastatic rate, most immune-excluded.
  The depth axis runs:
    Normal colonocyte (CDH1+, CDX2+, non-motile)
    → CMS4 shallow (TGF-β elevated, partial EMT)
    → CMS4 deep (full EMT, mesenchymal, immune-excluded,
      metastatic programme active)
  The structural analogy within this repository:
    CMS4 is the CRC equivalent of the mesenchymal
    shift at GBM recurrence — a treatment-driven
    (and sometimes intrinsic) shift from an epithelial
    to a mesenchymal false attractor state.
    In GBM: proneural → mesenchymal (at recurrence)
    In CRC: CMS2 → CMS4-like (at metastasis or
    progression) — some evidence that CMS2 tumours
    acquire CMS4-like characteristics at metastatic sites.

AVAILABLE PUBLIC DATA:
  TCGA-COAD/READ:   CMS4 ~23% of ~600 = ~138 samples
  GSE39582:         ~135 CMS4 — best powered CMS4 dataset
  GSE14333:         ~67 CMS4
  GSE17536:         ~41 CMS4
  CAF-specific:
    GSE35602:       CAF vs tumour cell sorted fractions
                    Critical for deconvolution of
                    CMS4 stromal signal

ANALYSIS STATUS:    NOT YET STARTED
PLANNED DOCUMENT:   CRC-S3a (CMS4 before-document)
PLANNED DATASET:    GSE39582 CMS4 subset
                    + TCGA-COAD CMS4 subset
PRIORITY:           THIRD — CMS4 has the worst prognosis
                    and the greatest unmet need, but the
                    stromal contamination problem makes
                    it analytically complex. Tackle after
                    the cleaner CMS2 and CMS1 analyses
                    establish the framework for CRC.
```

---

## SECTION VI — CMS3: THE METABOLIC SUBTYPE

```
PREVALENCE: 13%

CLINICAL FACTS:
  Location:         Mixed — right and left colon
  Histology:        Moderately differentiated.
                    Intermediate between CMS2 and CMS4.
                    Sometimes mixed features.
  MSI status:       MIXED — both MSS and low MSI occur.
                    ~30% of CMS3 has low MSI.
                    NOT the high MSI of CMS1.
  Prognosis:        INTERMEDIATE — better than CMS4,
                    comparable to or slightly worse
                    than CMS2.
  Standard therapy: FOLFOX / CAPOX
                    Anti-EGFR therapy: less clear benefit
                    (KRAS mutation common → EGFR inhibitors
                    not applicable in KRAS-mutant CMS3)
  Response to immunotherapy: LOW (mostly MSS)

CELL OF ORIGIN:
  LGR5+ ISC — same starting point.
  CMS3 is the subtype where the metabolic reprogramming
  of the false attractor is the defining feature rather
  than WNT (CMS2), immune evasion (CMS1), or EMT (CMS4).

DEFINING MOLECULAR EVENTS:
  KRAS mutation:    ~68% of CMS3 — the HIGHEST KRAS
                    frequency of any CRC subtype.
                    KRAS is the master driver of CMS3.
                    KRAS G12D and G12V are most common.
                    KRAS activates:
                      RAS/RAF/MEK/ERK (proliferation)
                      PI3K/AKT/mTOR (survival)
                      RAL-GDS pathway (invasion)
                      Metabolic reprogramming:
                        Glycolysis upregulation
                        Glutamine addiction
                        Fatty acid synthesis
                        Pentose phosphate pathway
  KRAS-driven metabolic programme:
                    The entire metabolic signature of
                    CMS3 flows from KRAS:
                    KRAS → MYC upregulation → glutamine
                    metabolism increase → nucleotide
                    biosynthesis → increased proliferation
                    KRAS → LKB1 pathway suppression →
                    metabolic flexibility
                    KRAS → NRF2 activation → ROS defence
  APC mutation:     ~50% — intermediate between CMS1
                    and CMS2
  PIK3CA mutation:  ~25% — highest of any CRC subtype
                    (in addition to KRAS — both drivers
                    of PI3K/AKT active)
  Mixed epigenetics: Low or intermediate CIMP
                    (CMS3 sits between CMS1-CIMP-H and
                    CMS2-non-CIMP)
  ERBB2 (HER2):     ~5% amplification in CRC overall
                    but enriched in CMS3

KEY METABOLIC GENES (the CMS3 signature):
  FASN:             Fatty acid synthase — high in CMS3
  ACACA:            Acetyl-CoA carboxylase — lipid synthesis
  LDHA:             Lactate dehydrogenase — glycolysis
  GLS:              Glutaminase — glutamine → glutamate
  SLC1A5 (ASCT2):   Glutamine transporter — elevated
  GLUT1 (SLC2A1):   Glucose transporter — elevated
  These metabolic genes define the CMS3 depth axis:
  the deeper the CMS3 tumour, the more extreme the
  metabolic reprogramming.

WADDINGTON STRUCTURE:
  CMS3 is the metabolic false attractor:
  a cell whose Waddington position is defined primarily
  by its rewired metabolism rather than by loss of
  differentiation (CMS2), immune activation (CMS1),
  or mesenchymal transition (CMS4).
  The depth axis in CMS3:
    Normal colonocyte (oxidative phosphorylation,
    normal glucose/glutamine metabolism)
    → CMS3 shallow (KRAS mutant, partial metabolic shift)
    → CMS3 deep (full Warburg effect, glutamine-addicted,
      lipid synthesis maximal, PI3K + KRAS co-active)
  The cross-cancer structural connection:
    The metabolic reprogramming in CMS3 is structurally
    analogous to the TCA disruption in RCC and the IDH
    mutation in glioma — different entry points into the
    same metabolic vulnerability pattern.
    In RCC: TCA enzyme suppression → αKG reduction
    In IDH-mutant glioma: IDH neomorphism → 2-HG
    In CMS3 CRC: KRAS → metabolic rewiring → glutamine
    and lipid addiction
    All three converge on altered intermediate metabolism
    as a false attractor maintenance mechanism.
    This cross-cancer connection will be tested formally
    in the CMS3 before-document.

AVAILABLE PUBLIC DATA:
  TCGA-COAD/READ:   CMS3 ~13% of ~600 = ~78 samples
                    Small within TCGA alone.
  GSE39582:         ~76 CMS3 — best available for CMS3
  GSE14333:         ~38 CMS3
  KRAS-specific:
    GSE26906:       KRAS-mutant vs wildtype CRC
    GSE62080:       Metabolic gene expression CRC

ANALYSIS STATUS:    NOT YET STARTED
PLANNED DOCUMENT:   CRC-S4a (CMS3 before-document)
PLANNED DATASET:    GSE39582 CMS3 subset
                    + TCGA-COAD CMS3 subset
PRIORITY:           FOURTH — the smallest subtype by n.
                    The metabolic axis is the connecting
                    thread to RCC (TCA) and IDH-mutant
                    glioma. Analytically tractable but
                    lower n than CMS1/2/4.
```

---

## SECTION VII — THE MSI/MSS SPLIT AND ITS CLINICAL WEIGHT

```
The single most clinically important binary in CRC
molecular classification is:

  MSI-HIGH (dMMR) vs. MICROSATELLITE STABLE (MSS/pMMR)

This is not the same as the CMS classification —
but they overlap strongly:

  CMS1 = almost entirely MSI-H (~85% of CMS1)
  CMS2 = almost entirely MSS
  CMS3 = mixed (~30% MSI-low, mostly MSS)
  CMS4 = almost entirely MSS

CLINICAL WEIGHT OF THE MSI/MSS BINARY:

  MSI-H CRC (dMMR):
    ~15% of all CRC
    ~3% caused by Lynch syndrome (germline MMR mutation)
    ~12% sporadic (MLH1 hypermethylation, BRAF V600E)
    Standard test: immunohistochemistry for MMR proteins
                   (MLH1, MSH2, MSH6, PMS2) OR PCR for
                   microsatellite markers OR NGS-TMB
    Treatment implication:
      Stage II:   Do NOT give 5-FU (may worsen outcome)
      Stage IV:   PEMBROLIZUMAB first-line (KEYNOTE-177)
                  Complete response rate ~10-15%
                  2-year PFS ~48% (vs ~19% with chemo)
                  UNPRECEDENTED for metastatic CRC

  MSS CRC:
    ~85% of all CRC
    Does NOT respond to single-agent PD-1 blockade
    Current trials: combinations of TGF-β inhibitors
    + PD-1 inhibitors specifically targeting CMS4 MSS CRC
    (LSZ102, bintrafusp alfa, others)

THE TESTING IMPLICATION FOR THIS FRAMEWORK:
  The depth score in CMS1 (MSI-H) should be testable
  as a predictor of pembrolizumab response:
  Hypothesis: shallow CMS1 (low depth) = good
  pembrolizumab response (tumour still immunogenic,
  PD-1 blockade releases functioning T cells).
  Deep CMS1 = exhausted immune landscape, primary
  resistance to pembrolizumab.
  This is the highest-value clinical prediction the
  CMS1 analysis could produce.
  It is stated here as a structural hypothesis —
  the before-document will formalise predictions.
```

---

## SECTION VIII — THE SIDEDNESS AXIS — RIGHT vs. LEFT COLON

```
The right colon (cecum, ascending, transverse) and
left colon (descending, sigmoid) are embryologically
distinct:

  RIGHT colon:  Derived from midgut
                Longer, more capacious, absorptive
                Different microbiome composition
                More alkaline luminal environment
                Enriched for CMS1 and CMS3
                APC mutation rate: lower
                BRAF V600E rate: higher
                MSI-H rate: higher
                Female sex predominance at older ages

  LEFT colon and rectum:  Derived from hindgut
                Shorter, tubular, motility function
                Enriched for CMS2
                APC mutation rate: higher
                KRAS rate: lower (in CMS2)
                MSI rate: lower
                Male sex predominance

CLINICAL CONSEQUENCE:
  Anti-EGFR therapy (cetuximab, panitumumab) works
  in LEFT-sided RAS-wildtype metastatic CRC
  but NOT in right-sided CRC — regardless of RAS status.
  This was established in retrospective analyses of
  PRIME, CRYSTAL, and FIRE-3 trials.

  The mechanism is incompletely understood but reflects:
  - Different EGFR pathway dependency by side
  - Different co-mutation landscape (BRAF, PI3K)
  - Different ligand expression (amphiregulin,
    epiregulin — higher on left)
  - CMS subtype distribution (CMS2 left-predominant;
    CMS2 is most anti-EGFR responsive)

THE SIDEDNESS AXIS IN THE FRAMEWORK:
  The depth score analysis should test sidedness as
  a correlate — are depth-positive genes enriched
  for right- or left-sided tumours?
  This is not a prediction — it is a structural
  analytical step to be included in each CMS
  before-document.
```

---

## SECTION IX — THE EXISTING CRC ANALYSIS — CONTEXT

```
The existing analysis in Cancer_Research/CRC/ ran on
the combined TCGA-COAD/READ dataset without CMS
stratification.

WHAT THE EXISTING ANALYSIS CAPTURED:
  TCGA-COAD/READ contains approximately:
    CMS1: ~84 samples (14%)
    CMS2: ~222 samples (37%)
    CMS3: ~78 samples (13%)
    CMS4: ~138 samples (23%)
    Unclassified: ~78 samples (13%)

  The dominant variance in the combined dataset is
  driven by CMS2 (the largest group) — so the depth
  axis from the existing analysis is predominantly
  the WNT/MYC axis of CMS2.

  HOWEVER: the combined dataset also contains the
  immune signal of CMS1 and the mesenchymal signal
  of CMS4. These create two competing signals that
  partially cancel each other:
    CMS1: immune genes UP with depth
    CMS4: immune genes paradoxically DOWN (excluded)
  In the combined analysis, these opposing immune
  signals may have flattened the immune architecture.

  The genes that appeared robustly depth-positive
  in the existing analysis are likely:
    Universal CRC proliferation markers: MKI67, TOP2A
    WNT target genes: MYC, CCND1, LGR5
    Epigenetic lock markers: EZH2 (predicted)
    ECM/invasion markers: COL1A1, MMP genes

  The genes that appeared depth-negative:
    Colonocyte identity markers: CDX2, KRT20, SLC26A3
    Normal crypt markers: OLFM4, LGR5 (paradoxically
    lost at deepest depth despite being ISC marker —
    because deep CMS4 loses ALL epithelial markers)

THE EXISTING ANALYSIS IS NOT DISCARDED.
It correctly identified the central CRC depth axis.
The CMS subtype analyses will refine it — showing
which parts of the combined signal belong to which
subtype.
```

---

## SECTION X — DATA AVAILABILITY SUMMARY

```
Subtype  Primary Dataset    n (tumour)   CMS label   Normal ref   Power
─────────────────────────────────────────────────────────────────────────
CMS2     GSE39582 CMS2      ~216         YES         GTEx colon   HIGH
         TCGA-COAD CMS2     ~222         YES         TCGA adj     HIGH
CMS1     GSE39582 CMS1      ~82          YES         GTEx colon   MOD
         TCGA-COAD CMS1     ~84          YES         TCGA adj     MOD
CMS4     GSE39582 CMS4      ~135         YES         GTEx colon   MOD
         TCGA-COAD CMS4     ~138         YES         TCGA adj     MOD
CMS3     GSE39582 CMS3      ~76          YES         GTEx colon   LOW-MOD
         TCGA-COAD CMS3     ~78          YES         TCGA adj     LOW-MOD
Combined TCGA-COAD/READ     ~600         YES (retro) TCGA adj     HIGH
         (existing analysis)

CRITICAL NOTE ON CMS LABELS:
  CMS labels are NOT in the TCGA clinical files natively.
  They must be applied by running the CMScaller R package
  (Eide et al., 2017) or using the published classifier
  on the expression matrix.
  CMScaller is freely available and runs on bulk RNA-seq.
  This must be done as the FIRST STEP of each CMS
  subtype analysis before any script runs.
  The before-document for each subtype will specify
  that CMS assignment is a prerequisite data step.

NORMAL COLON REFERENCE:
  GTEx colon (transverse + sigmoid): gold standard.
  TCGA adjacent normal tissue: available but proximity
  to tumour may introduce field cancerisation effects —
  use as validation, not primary reference.
```

---

## SECTION XI — PLANNED ANALYSIS ORDER

```
ORDER:

  CRC-S1   CMS2     GSE39582 + TCGA CMS2        HIGH POWER
                    REASON: Largest subtype.
                    Canonical CRC. WNT/MYC axis.
                    Anti-EGFR response prediction
                    is the clinical output.
                    The existing combined analysis
                    most reflects CMS2 — this will
                    confirm and refine it.

  CRC-S2   CMS1     GSE39582 + TCGA CMS1        MODERATE
                    REASON: Immunotherapy-responsive.
                    The depth score in CMS1 maps the
                    T cell exhaustion axis.
                    Highest clinical utility:
                    pembrolizumab response prediction
                    in dMMR CRC.
                    The BRAF/CIMP epigenetic entry
                    point connects to IDH glioma.

  CRC-S3   CMS4     GSE39582 + TCGA CMS4        MODERATE
                    REASON: Worst prognosis.
                    TGF-β / EMT / immune exclusion axis.
                    Bevacizumab response prediction.
                    Stromal contamination must be
                    addressed in the before-document.

  CRC-S4   CMS3     GSE39582 + TCGA CMS3        LOW-MOD
                    REASON: KRAS/metabolic axis.
                    Connects to RCC (TCA) and IDH
                    glioma (metabolic reprogramming).
                    Smallest n — requires care.

  CRC-X    Cross-   After CMS1–4 complete
           CMS      Questions:
                      1. Universal CRC depth genes
                         across all four subtypes?
                      2. EZH2 in all four? Direction?
                      3. CDX2 loss as universal depth
                         marker (expected — but does
                         the rate differ by CMS)?
                      4. Basket trial candidates?
                         What works across CMS1–4?
                      5. Does the combined depth score
                         (existing analysis) decompose
                         cleanly into four CMS axes?
```

---

## SECTION XII — WHAT THIS DOCUMENT DOES NOT CONTAIN

```
This document contains:
  ✓ CMS1–4 subtype definitions and clinical facts
  ✓ Cells of origin (LGR5+ ISC for all; CMS4 complex)
  ✓ Published molecular drivers (literature only)
  ✓ Normal crypt hierarchy (the Waddington baseline)
  ✓ Clinical characteristics (survival, standard of care)
  ✓ The MSI/MSS binary and its treatment implications
  ✓ The sidedness axis (right vs. left colon)
  ✓ Structural connections to RCC and IDH glioma
    (CMS1-CIMP ↔ IDH glioma; CMS3-metabolism ↔ RCC TCA)
    stated as observations, not predictions
  ✓ Context for the existing combined analysis
  ✓ Data availability (GSE39582, TCGA, GTEx)
  ✓ CMS label assignment prerequisite note

This document does NOT contain:
  ✗ Depth score predictions
  ✗ Specific false attractor gene predictions
  ✗ Drug target predictions
  ✗ Epigenetic lock mechanism hypotheses
  ✗ Cross-CMS structural predictions beyond literature

All of the above belong in the BEFORE documents.
CRC-S1a (CMS2 before-document) is next.
Written before any script runs. Before any data loads.
```

---

## STATUS BLOCK

```
document:           CRC_Subtype_Orientation.md
folder:             Cancer_Research/CRC/Subtypes/
status:             COMPLETE — ORIENTATION ONLY
date:               2026-03-03
author:             Eric Robert Lawson
                    OrganismCore

subtypes_covered:
  CMS2:             Canonical / WNT-MYC               [1 of 4]
  CMS1:             MSI-Immune / BRAF-CIMP             [2 of 4]
  CMS4:             Mesenchymal / TGF-β               [3 of 4]
  CMS3:             Metabolic / KRAS                   [4 of 4]

analyses_started:   0

existing_analysis:  Cancer_Research/CRC/ — conducted on
                    combined TCGA-COAD/READ, no CMS split.
                    Depth axis predominantly reflects CMS2
                    (largest subtype).
                    Decomposed by analyses in this folder.

cms_label_note:     CMS labels must be computed from
                    expression data using CMScaller
                    before any subtype analysis begins.
                    This is a prerequisite data step,
                    not a prediction.

next_document:      CRC-S1a
                    CMS2 Before-Document
                    (predictions locked before
                    GSE39582 / TCGA CMS2 subset loads)

critical_note_1:    CMS2 drives the existing combined
                    analysis signal. CMS1 is the only
                    immunotherapy-responsive subtype.
                    CMS4 has the worst prognosis.
                    CMS3 connects to the metabolic axis
                    of the broader repository.
                    These are four separate Waddington
                    landscapes sharing one organ and one
                    cell of origin (LGR5+ ISC) but with
                    four distinct false attractors.

critical_note_2:    CMS4 analysis must account for
                    stromal contamination — the bulk
                    RNA-seq signal contains CAF and
                    macrophage transcriptomes in addition
                    to tumour cell signal. This must be
                    addressed in CRC-S3a.

critical_note_3:    The MSI-H / MSS binary is the most
                    clinically actionable split in CRC.
                    The depth score in CMS1 (MSI-H)
                    maps the T cell exhaustion axis and
                    is the highest-priority clinical
                    output of the CRC subtype series.
```
