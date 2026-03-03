# BREAST CANCER — SUBTYPE ORIENTATION DOCUMENT
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

  A complete map of the breast cancer subtype landscape —
  what the subtypes are, where they come from anatomically,
  what two parallel classification systems (histological
  and molecular) describe them, how those two systems
  relate to each other, what the clinical and molecular
  literature has established about each subtype, and
  what public data exists to analyze them.

Breast cancer has a classification problem more complex
than any other cancer type in this repository.
It has TWO classification systems that do not map
cleanly onto each other:
  1. Histological classification (what the cells look like)
  2. Molecular classification (what the cells express)

Both systems are used clinically. Neither alone is sufficient.
This document reconciles them before any analysis begins.

This document is the prerequisite to every analysis in the
BRCA/Subtypes/ folder.
```

---

## DOCUMENT METADATA

```
document_id:        BRCA_Subtype_Orientation
series:             BRCA (Breast Cancer — Subtypes)
folder:             Cancer_Research/BRCA/Subtypes/
date:               2026-03-03
author:             Eric Robert Lawson
status:             ORIENTATION ONLY — no predictions,
                    no analysis, no epigenetic hypotheses
next_document:      BRCA_LumA_before.md
                    (Document BRCA-S1a — Luminal A before-doc)
protocol_version:   Workflow_Protocol.md v2.0
```

---

## SECTION I — THE FUNDAMENTAL PROBLEM WITH "BREAST CANCER"

```
"Breast cancer" is not one disease.

This is now universally acknowledged in oncology — and yet
clinical trials, treatment protocols, and survival statistics
are still frequently reported as if it is.

The breast is composed of multiple cell types organized
into a ductal-lobular tree:
  - Luminal epithelial cells lining the ducts and lobules
  - Basal/myoepithelial cells surrounding the luminal cells
  - Stromal cells, fibroblasts, adipocytes (not cancer cells)
  - Stem/progenitor cell populations at branch points

Breast cancers arise from different cells within this tree,
travel through different molecular transitions, and arrive
at different false attractors.

The clinical consequence:
  - Luminal A: cured by hormone blockade in many patients
  - TNBC (basal-like): no hormone receptor to block,
    no HER2 to target, cytotoxic chemotherapy is the
    only approved option for most patients
  - HER2-enriched: transformed by trastuzumab after 1998
    (from 18-month median OS to near-normal in HER2+ patients)
  - Invasive lobular: grows silently, resists mammography,
    often presents at large size with diffuse spread

One protocol does not serve all of these.
The framework in this repository treats each subtype as a
separate disease with a separate normal cell of origin,
a separate Waddington landscape, and a separate false
attractor. The cross-subtype comparison comes after
each individual analysis is complete.
```

---

## SECTION II — THE TWO CLASSIFICATION SYSTEMS

```
BREAST CANCER HAS TWO PARALLEL CLASSIFICATION SYSTEMS.
Both are used clinically. They are not the same thing.
Understanding how they relate is essential before analysis.

SYSTEM 1 — HISTOLOGICAL CLASSIFICATION
  Based on: What the tumour cells look like under
            the microscope. What structures they form.
  Primary divide:
    Invasive Ductal Carcinoma (IDC):   ~75% of breast cancer
                                        Also called "invasive
                                        carcinoma of no special
                                        type" (NST) in WHO 2022
    Invasive Lobular Carcinoma (ILC):  ~10–15% of breast cancer
    Rare special types:                ~10–15% combined
      Tubular, cribriform, mucinous,
      medullary, metaplastic, etc.

SYSTEM 2 — MOLECULAR CLASSIFICATION (PAM50)
  Based on: Gene expression profile of 50 key genes
            (the PAM50 assay / Prosigna)
  Primary divide (5 subtypes):
    Luminal A
    Luminal B
    HER2-enriched
    Basal-like
    Normal-like (rare, likely technical artifact in some cases)

  Additional molecular class (not formally in PAM50):
    Claudin-low      (stem cell-like, triple-negative)

HOW THEY RELATE:
  The two systems are correlated but not equivalent.

  IDC can be ANY molecular subtype:
    Luminal A IDC    (most common)
    Luminal B IDC
    HER2-enriched IDC
    Basal-like IDC   (most common presentation of TNBC)

  ILC is ALMOST ALWAYS luminal:
    ~95% of ILC is Luminal A or Luminal B
    ILC that is TNBC or HER2-enriched is extremely rare
    This is because the molecular event defining ILC
    (CDH1/E-cadherin loss) is incompatible with the
    basal-like programme

  The special types map as follows:
    Mucinous:          Almost always Luminal A
    Tubular:           Almost always Luminal A
    Medullary:         Often basal-like / TNBC
    Metaplastic:       Almost always TNBC (basal-like
                       or claudin-low)

FOR THIS FRAMEWORK:
  The analysis will use MOLECULAR subtypes as the primary
  classification, because:
    1. Molecular subtype determines drug target
       (ER/PR → endocrine therapy, HER2 → anti-HER2,
       TNBC → chemotherapy / PARP inhibitors)
    2. Molecular subtype determines depth axis
       (the normal cell programme being replaced differs
       by molecular subtype)
    3. Waddington landscape geometry is a molecular
       phenomenon, not a histological one

  ILC will be treated as a SEPARATE ANALYSIS because:
    Although ILC is almost always Luminal A/B molecularly,
    it has a different cell of origin (lobular epithelium),
    a different architectural programme (E-cadherin loss,
    single-file invasion), and a different clinical behavior
    (diffuse spread, late recurrence, endocrine therapy
    resistance pattern). The landscape geometry of ILC
    is not the same as IDC-Luminal A despite the similar
    PAM50 assignment.
```

---

## SECTION III — THE BREAST DUCTAL-LOBULAR ARCHITECTURE

```
Before defining cells of origin, the normal breast
architecture must be understood — because the Waddington
landscape geometry is defined relative to the normal cell.

THE NORMAL BREAST DUCTAL TREE:

  Terminal Duct Lobular Unit (TDLU):
    The functional unit of the breast.
    Where milk is produced and secreted.
    Composed of:
      Luminal cells:    Line the duct lumen
                        Produce milk proteins
                        Express ER, PR, CK8/18
                        Respond to estrogen/progesterone
      Myoepithelial cells: Surround the luminal layer
                        Contractile (squeeze milk out)
                        Express CK5/6, CK14, p63, SMA
                        Do NOT express ER
                        The "basal" layer of the duct

  Stem/progenitor cells:
    Present at ductal branch points and within the TDLU
    Give rise to both luminal and myoepithelial cells
    Express CD44 (high), CD24 (low) — stem cell phenotype
    Express CK5/6, CK14 (basal markers)
    BRCA1-mutated cancers arise partly from this population

THE HIERARCHY:
  Mammary stem cell
    → Luminal progenitor → Mature luminal cell
    → Myoepithelial progenitor → Mature myoepithelial cell

  This hierarchy is the Waddington landscape of the
  normal breast. The normal attractors are:
    Mature luminal cell:      Waddington valley 1
    Mature myoepithelial cell: Waddington valley 2

  Cancer is what happens when a cell in this hierarchy
  escapes the valley it is in and enters a false attractor
  that supports malignant proliferation.

  Luminal A and B: Mature luminal cell escapes partially
  HER2-enriched:   Luminal progenitor with HER2 amplification
  Basal-like TNBC: Luminal progenitor or myoepithelial
                   escape with basal programme activated
  Claudin-low:     Mammary stem cell-like false attractor
  ILC:             Lobular epithelial cell with E-cadherin
                   loss — a different normal starting point
```

---

## SECTION IV — THE MOLECULAR SUBTYPES IN DETAIL

---

### SUBTYPE 1 — LUMINAL A

```
CLINICAL FACTS:
  Prevalence:       ~40–50% of all breast cancer
                    Largest single subtype.
  ER/PR/HER2 status: ER+, PR+ (usually high), HER2-
  Ki-67:            Low (<14% by conventional threshold)
                    Low proliferation index
  Grade:            Usually Grade 1 or 2
  5-year survival:  ~90%+ (all stages combined)
                    Best prognosis of all subtypes
  Standard of care: Endocrine therapy (tamoxifen or
                    aromatase inhibitors) — primary treatment
                    Chemotherapy: often NOT needed
                    (Oncotype DX / genomic scores guide
                    chemotherapy decisions)
                    CDK4/6 inhibitors (palbociclib,
                    ribociclib, abemaciclib) in combination
                    with endocrine therapy for metastatic disease
  Recurrence:       CAN recur late — 10–20 years after
                    initial diagnosis. Unique to luminal subtypes.
                    Late recurrence is not seen in TNBC or
                    HER2-enriched at the same rate.

CELL OF ORIGIN:
  Mature luminal epithelial cell
  The most differentiated, most committed luminal cell.
  Normal identity programme:
    ESR1 (ER):      Master TF — estrogen receptor alpha
                    Drives luminal cell identity and
                    hormone-responsive gene programmes
    FOXA1:          Pioneer TF that opens chromatin for
                    ER binding. Co-expressed with ER.
                    FOXA1 is the gatekeeper of luminal identity.
    GATA3:          TF required for luminal differentiation
                    Co-expressed with ER/FOXA1
    CK8/18:         Luminal cytokeratins
    ELF5:           Luminal alveolar TF (milk production)
    PGR (PR):       Progesterone receptor — target of ER
    CCND1:          Cell cycle, ER-responsive
  The luminal A cell is highly differentiated, slow cycling,
  estrogen-dependent for survival.

INITIATING MOLECULAR EVENTS:
  PIK3CA mutation:  ~45% of Luminal A (most common)
  MAP3K1 mutation:  ~14%
  CDH1 mutation:    In the ILC subset (~10% of luminal A)
  GATA3 mutation:   ~10%
  TP53 mutation:    ~12% (low — unlike basal-like where it
                    is nearly universal)
  Few copy number alterations (unlike HGSC or TNBC)
  Genomic stability: Relatively stable genome

KEY DISTINGUISHING FEATURES:
  - ER+ PR+ HER2- Low Ki-67
  - FOXA1/GATA3/ESR1 co-expression (the luminal identity TFs)
  - PIK3CA most common driver
  - Slow growing, highly differentiated
  - Hormone-dependent: estrogen deprivation is therapeutic
  - Late recurrence risk (10–20 year horizon)
  - CDK4/6 inhibitor-sensitive in metastatic disease

AVAILABLE PUBLIC DATA:
  TCGA-BRCA:        n~500–600 Luminal A samples
                    PAM50 calls available in clinical file
                    RNA-seq + methylation + CNV
  GSE96058:         n=3,273 total; large luminal A component
                    Clinical metadata for subtype calls
  GSE25066:         Pre-treatment samples, chemotherapy
                    response annotated

ANALYSIS STATUS:    NOT YET STARTED
PLANNED DOCUMENT:   BRCA-S1a (Luminal A before-document)
PLANNED DATASET:    TCGA-BRCA (PAM50 Luminal A subset)
                    Normal reference: GTEx breast tissue
                    or TCGA adjacent normal
```

---

### SUBTYPE 2 — LUMINAL B

```
CLINICAL FACTS:
  Prevalence:       ~20% of all breast cancer
  ER/PR/HER2 status: ER+, PR+ (often lower than Luminal A),
                    HER2- OR HER2+
                    Two variants:
                      Luminal B HER2-: ER+, HER2-, high Ki-67
                      Luminal B HER2+: ER+, HER2+, high Ki-67
  Ki-67:            High (>20% by most thresholds)
  Grade:            Usually Grade 2 or 3
  5-year survival:  ~80–85%
                    Worse than Luminal A due to higher
                    proliferation and more aggressive behavior
  Standard of care: Endocrine therapy (required — ER+)
                    PLUS chemotherapy (more often needed
                    than in Luminal A — higher Ki-67 drives
                    chemotherapy recommendation)
                    HER2-targeted therapy if HER2+
                    CDK4/6 inhibitors (metastatic)

CELL OF ORIGIN:
  Luminal epithelial cell — but less differentiated than
  Luminal A. Likely a luminal progenitor rather than
  a terminally differentiated luminal cell.
  The key distinction from Luminal A:
    Luminal A = terminally differentiated luminal cell
                that has partly escaped ER control
    Luminal B = luminal progenitor that has not fully
                committed to terminal differentiation AND
                has acquired proliferative drive
  Normal identity markers: same as Luminal A but
  expressed at lower levels (less committed to the
  mature luminal identity)

INITIATING MOLECULAR EVENTS:
  TP53 mutation:    Higher than Luminal A (~32%)
  PIK3CA mutation:  ~29%
  MAP3K1:           Less common than in Luminal A
  ERBB2 amplification: in Luminal B HER2+ variant
  RB1 alterations:  Some Luminal B tumors
  Higher copy number alterations than Luminal A
  (more genomically unstable)

KEY DISTINGUISHING FEATURES:
  - ER+ but HIGHER proliferation than Luminal A
  - FOXA1/GATA3 expressed but at lower levels
  - Less hormone-dependent than Luminal A —
    endocrine therapy less effective as monotherapy
  - Chemotherapy more often required
  - Luminal B HER2+ subset: overlap with HER2-enriched
    clinically (gets both endocrine and HER2-targeted tx)
  - Worse prognosis than Luminal A stage-for-stage

AVAILABLE PUBLIC DATA:
  TCGA-BRCA:        n~200–250 Luminal B (PAM50 calls)
  GSE96058:         Includes Luminal B with subtype calls
  NOTE: Luminal B and Luminal A are separated by PAM50
        gene expression profiling, not just Ki-67 alone.
        Ki-67 alone misclassifies ~20% of cases.

ANALYSIS STATUS:    NOT YET STARTED
PLANNED DOCUMENT:   BRCA-S2a (Luminal B before-document)
PLANNED DATASET:    TCGA-BRCA (PAM50 Luminal B subset)
STRUCTURAL NOTE:    Luminal A and Luminal B share a cell
                    of origin but differ in differentiation
                    state. This is structurally analogous
                    to ccRCC/PRCC sharing a proximal tubule
                    origin but diverging at different
                    saddle points. The cross-subtype
                    comparison of Luminal A vs. Luminal B
                    will be the first cross-BRCA analysis.
```

---

### SUBTYPE 3 — HER2-ENRICHED

```
CLINICAL FACTS:
  Prevalence:       ~10–15% of all breast cancer
  ER/PR/HER2 status: HER2+, ER-, PR-
                    IMPORTANT: Not all HER2+ tumors are
                    "HER2-enriched" by PAM50.
                    ~50% of HER2+ by IHC/FISH are actually
                    Luminal B by PAM50 (they co-express ER).
                    "HER2-enriched" by PAM50 means
                    HER2-high / ER-negative.
  Grade:            Almost always Grade 3
  5-year survival:  Pre-trastuzumab era: ~50–60%
                    Post-trastuzumab era: ~85–90%
                    Trastuzumab is one of the most
                    transformative targeted therapies in
                    oncology history.
  Standard of care: Anti-HER2 therapy MANDATORY:
                    Trastuzumab + pertuzumab (dual HER2
                    blockade, standard neoadjuvant)
                    T-DM1 (ado-trastuzumab emtansine)
                    in residual disease after neoadjuvant
                    T-DXd (trastuzumab deruxtecan) —
                    second-line metastatic
                    Tucatinib (small molecule HER2 kinase
                    inhibitor) — metastatic brain mets
                    Chemotherapy: anthracycline + taxane
                    backbone
  pCR rate:         ~40–60% with dual HER2 blockade +
                    chemotherapy (highest pCR of all
                    subtypes except TNBC in some trials)

CELL OF ORIGIN:
  Luminal progenitor cell with HER2 amplification.
  Not a fully differentiated luminal cell (ER-negative).
  Not a basal/myoepithelial cell.
  A cell that is partway through luminal differentiation
  and has been diverted by HER2 amplification into a
  different proliferative identity.
  Normal identity markers:
    GRB7:           Co-amplified with ERBB2 on 17q12
    ERBB2 (HER2):   Amplified and overexpressed
                    The amplicon at 17q12 drives the
                    false attractor — it is a copy number
                    event, not a point mutation
    CK7, CK18:      Luminal markers retained
    TP53:           Frequently mutated (>72%)

INITIATING MOLECULAR EVENTS:
  ERBB2 amplification: 17q12 amplicon (~100% of HER2-enriched)
                    This is the founding event — not a
                    point mutation but a massive gene
                    amplification driving constitutive
                    HER2/HER3/PI3K signaling
  TP53 mutation:    ~72% (second highest after TNBC)
  PIK3CA mutation:  ~24%
  GRB7 co-amplification with ERBB2 (co-amplicon)
  High genomic instability (copy number driven, like HGSC)

KEY DISTINGUISHING FEATURES:
  - ERBB2 amplification = defining event
  - TP53 co-mutation in majority (genomically unstable)
  - ER-negative (PAM50 HER2-enriched definition)
  - Responds dramatically to HER2-targeted therapy
  - Without HER2-targeted therapy: aggressive disease
  - With HER2-targeted therapy: excellent prognosis
  - The most dramatic example in breast oncology of
    a targetable false attractor switch gene (ERBB2)

AVAILABLE PUBLIC DATA:
  TCGA-BRCA:        n~80–100 HER2-enriched (PAM50)
                    (smaller sample size than luminal)
  GSE55348:         HER2-enriched enriched dataset
  GSE96058:         Includes HER2-enriched with subtype calls
  NOTE: HER2-enriched is relatively rare in TCGA-BRCA.
        Power will be moderate.

ANALYSIS STATUS:    NOT YET STARTED
PLANNED DOCUMENT:   BRCA-S3a (HER2-enriched before-document)
PLANNED DATASET:    TCGA-BRCA (PAM50 HER2-enriched subset)
```

---

### SUBTYPE 4 — BASAL-LIKE / TRIPLE-NEGATIVE (TNBC)

```
CLINICAL FACTS:
  Prevalence:       ~15–20% of all breast cancer
  ER/PR/HER2 status: ER-, PR-, HER2- (triple negative)
                    NOTE: ~80% of TNBC are basal-like by
                    PAM50. ~20% of TNBC are claudin-low
                    or HER2-enriched by PAM50.
                    And ~20% of basal-like by PAM50 are
                    NOT triple-negative by IHC.
                    The terms are overlapping but not
                    synonymous. This analysis will use
                    PAM50 basal-like as the primary class.
  Grade:            Almost always Grade 3
  Patient demographic: Higher prevalence in:
                    - Pre-menopausal women
                    - African-American women (~40% of
                      breast cancer in AA women is TNBC)
                    - BRCA1 germline mutation carriers
                      (~75% of BRCA1-associated breast
                      cancers are TNBC/basal-like)
  5-year survival:  ~77% (all stages, though dropping
                    significantly at Stage III/IV)
  Recurrence:       Recurrence is EARLY (within 3–5 years)
                    — not late like luminal subtypes.
                    If no recurrence within 5 years:
                    TNBC patients have similar long-term
                    survival to luminal patients.
                    This is the "TNBC paradox."
  Standard of care: Chemotherapy: anthracycline + taxane
                    Immune checkpoint:
                      Pembrolizumab (PD-L1+ or high TIL)
                      added to neoadjuvant chemotherapy
                      (KEYNOTE-522) — now standard
                    PARP inhibitors: olaparib/talazoparib
                      in BRCA1/2-mutated TNBC
                    Sacituzumab govitecan (TROP2 ADC) —
                      second-line metastatic
                    Antibody-drug conjugates increasingly
                    important (TROP2, LIV-1)
  pCR rate:         ~35–50% with chemotherapy ± immunotherapy
                    If pCR achieved: excellent outcomes
                    If residual disease: poor prognosis
                    This is the "neoadjuvant window" that
                    stratifies TNBC prognosis better than
                    any pre-treatment biomarker

CELL OF ORIGIN:
  THE MOST DEBATED CELL OF ORIGIN IN BREAST CANCER.
  Two competing hypotheses with evidence for both:

  HYPOTHESIS 1 — LUMINAL PROGENITOR ORIGIN:
    A luminal progenitor that fails to complete
    differentiation. The basal-like programme reflects
    a progenitor state, not a myoepithelial state.
    Evidence: BRCA1-null mammary cells are luminal
    progenitors. BRCA1 loss in a luminal progenitor
    drives basal-like gene expression.
    The cell looks "basal" but comes from the luminal
    progenitor pool.

  HYPOTHESIS 2 — MYOEPITHELIAL / BASAL CELL ORIGIN:
    The cancer arises from the myoepithelial/basal cell
    layer and retains basal markers (CK5/6, CK14, p63)
    because it came from that layer.
    Evidence: CK5/6+, CK14+, p63+ expression in TNBC
    mirrors the myoepithelial programme.

  CURRENT CONSENSUS (2024):
    Most TNBC / basal-like cancers arise from a
    BRCA1-regulated LUMINAL PROGENITOR that cannot
    complete differentiation and instead activates
    a basal-like transcriptional programme.
    The cell starts luminal but ends up in a
    false attractor that resembles a myoepithelial
    / basal state.
    For the Waddington framework: the normal starting
    cell is the luminal progenitor. The false attractor
    is a basal-like / stem-like state.

  Normal identity markers (luminal progenitor):
    BRCA1:          DNA damage repair, luminal progenitor
                    maintenance
    KRT5/6:         Basal cytokeratins (upregulated in
                    false attractor, absent in normal
                    luminal progenitor)
    SOX10:          Neural crest / basal cell TF —
                    upregulated in TNBC false attractor
    FOXC1:          Basal-like identity TF — upregulated
    EGFR:           Receptor tyrosine kinase — activated
                    in basal-like TNBC false attractor
    MKI67:          High proliferation in false attractor
    VIM:            Vimentin — mesenchymal marker, rises
                    with depth in TNBC

INITIATING MOLECULAR EVENTS:
  TP53 mutation:    ~80% of TNBC (highest of all subtypes)
  BRCA1 mutation:   ~20% germline + somatic in TNBC
  PIK3CA:           ~7% (lower than luminal)
  RB1:              ~20%
  PTEN loss:        ~35%
  CDKN2A:           ~30%
  High copy number instability (like HGSC)
  MYC amplification: ~40%

KEY DISTINGUISHING FEATURES:
  - Triple-negative (ER-, PR-, HER2-) by IHC
  - Basal cytokeratin expression (CK5/6, CK14)
  - EGFR overexpression (not mutated — overexpressed)
  - TP53 mutation near-universal
  - BRCA1 dysfunction (germline or somatic or
    promoter methylation) in large fraction
  - HIGHEST immune infiltration of all subtypes
    (TILs — tumor-infiltrating lymphocytes —
    predict pembrolizumab response)
  - Early recurrence, early death — or early cure
  - No hormone receptor to block
  - Chemotherapy is effective but resistance develops
  - Sacituzumab govitecan (TROP2 ADC) now standard
  - PARP inhibitors in BRCA-mutated patients
  - Greatest active clinical trial pipeline of all
    breast cancer subtypes as of 2026

INTERNAL HETEROGENEITY OF TNBC:
  TNBC is itself heterogeneous. Six TNBC subtypes
  have been described (Lehmann et al., 2011, 2016):
    BL1:  Basal-like 1 — cell cycle, DNA damage response
    BL2:  Basal-like 2 — growth factor signaling, glycolysis
    M:    Mesenchymal — EMT, differentiation
    MSL:  Mesenchymal stem-like — stem cell features
    IM:   Immunomodulatory — immune gene expression
    LAR:  Luminal androgen receptor — AR+, endocrine-like
  Of these, LAR (luminal androgen receptor) TNBC is
  the most clinically distinct — it expresses AR and
  may respond to androgen receptor blockade.
  This internal TNBC heterogeneity will need to be
  addressed in the TNBC analysis before-document.

AVAILABLE PUBLIC DATA:
  TCGA-BRCA:        n~140–200 basal-like (PAM50)
  GSE96058:         Includes TNBC with subtype calls
  GSE58812:         TNBC-enriched, survival annotated
  GSE25066:         Pre-treatment samples with pCR data
                    Critical for depth-score/pCR analysis
  GSE21653:         TNBC enriched
  NOTE: TNBC datasets with pCR annotations are
        the most clinically valuable — depth score
        prediction of pCR response is a primary goal.

ANALYSIS STATUS:    NOT YET STARTED
PLANNED DOCUMENT:   BRCA-S4a (TNBC/Basal-like before-doc)
PLANNED DATASET:    TCGA-BRCA (PAM50 Basal-like subset)
                    + GSE25066 (pCR annotation)
PRIORITY NOTE:      TNBC has the greatest unmet clinical
                    need of all molecular breast subtypes.
                    It is the second analysis priority
                    after Luminal A (highest sample size).
```

---

### SUBTYPE 5 — CLAUDIN-LOW

```
CLINICAL FACTS:
  Prevalence:       ~5–10% of all breast cancer
                    (not a formal PAM50 class — identified
                    by separate gene expression classifier)
  ER/PR/HER2 status: Usually triple-negative (ER-, PR-, HER2-)
                    BUT overlaps with all three negative
                    classifiers. Not all claudin-low is TNBC.
  Grade:            Grade 3
  5-year survival:  Similar to TNBC — poor at advanced stage
  Standard of care: Treated as TNBC (no specific approved
                    therapy for claudin-low as distinct subtype)
                    This is a gap — claudin-low responds
                    differently to chemotherapy than
                    conventional basal-like TNBC
                    (lower pCR rates to standard regimens)

CELL OF ORIGIN:
  MAMMARY STEM CELL — the most primitive cell in the
  breast hierarchy.
  Claudin-low tumors express stem cell markers:
    CD44 high / CD24 low  — breast cancer stem cell phenotype
    ALDH1                 — cancer stem cell marker
    Vimentin, fibronectin — mesenchymal markers
    SOX2, OCT4            — pluripotency TFs
  The false attractor in claudin-low is a stem/progenitor
  state with partial EMT (epithelial-mesenchymal transition).
  Low expression of cell-cell junction proteins
  (claudins 3, 4, 7 — hence "claudin-low").
  This is the breast cancer analogue of the undifferentiated
  state seen in other cancer types where the false attractor
  resembles a stem cell.

KEY DISTINGUISHING FEATURES:
  - Low expression of claudins 3, 4, 7 and E-cadherin
  - High immune infiltration (like TNBC/IM subtype)
  - Stem cell gene expression programme
  - EMT markers elevated
  - Lower pCR than conventional basal-like TNBC
  - No targeted therapy currently approved

AVAILABLE PUBLIC DATA:
  Not a formal PAM50 class — identified by separate
  gene expression classifier.
  Claudin-low calls available from published analyses
  of TCGA-BRCA and GSE96058.
  Small effective sample size for claudin-low-specific
  analysis — LOW POWER.

ANALYSIS STATUS:    NOT YET STARTED
PLANNED DOCUMENT:   BRCA-S5a (Claudin-low before-document)
                    AFTER Luminal A, TNBC, HER2-enriched
                    are complete.
PRIORITY:           Lower — small n, no specific clinical
                    target currently approved.
```

---

### SUBTYPE 6 — INVASIVE LOBULAR CARCINOMA (ILC)

```
CLINICAL FACTS:
  Prevalence:       ~10–15% of all breast cancer
                    Second most common histological subtype
  ER/PR/HER2 status: ~95% ER+, mostly HER2-
                    Almost always PAM50 Luminal A or B
  Grade:            Usually Grade 1 or 2
  5-year survival:  ~85–90% overall (similar to IDC-Luminal A)
                    BUT: late recurrence is HIGHER in ILC
                    than in IDC-Luminal A matched for ER+
                    This is the ILC paradox:
                    equal early prognosis, worse late prognosis
  Diagnostic challenge:
                    ILC grows as single-file strands of cells
                    (Indian file pattern) with no cohesion.
                    Mammography misses ILC more often than IDC.
                    ILC is often larger at diagnosis than
                    imaging suggests.
                    ILC spreads to unusual metastatic sites:
                    GI tract, peritoneum, meninges, ovaries
                    (unlike IDC which metastasizes to lung,
                    liver, bone, brain more typically)
  Standard of care: Same as Luminal A/B:
                    Endocrine therapy (tamoxifen or AI)
                    CDK4/6 inhibitors in metastatic disease
                    BUT: ILC responds LESS WELL to
                    neoadjuvant chemotherapy than IDC
                    (pCR rate ~5–10% vs. ~20% for IDC)
                    Chemotherapy is less effective in ILC
                    This is clinically significant.

CELL OF ORIGIN:
  LOBULAR EPITHELIAL CELL — the secretory cells within
  the breast lobule (acini), not the duct.
  This is a DIFFERENT CELL from the ductal luminal cell
  that gives rise to IDC.
  Key distinction of the lobular epithelial cell:
    - Expresses E-cadherin (CDH1) normally — this is what
      holds lobular cells together
    - Loss of E-cadherin is the DEFINING EVENT of ILC
    - When CDH1 is lost (mutation or promoter methylation),
      lobular cells lose cohesion and invade as single files

  Normal identity markers:
    CDH1 (E-cadherin):  LOST in ILC — the founding event
    ESR1 (ER):          High expression (retained)
    FOXA1:              High (luminal identity maintained)
    GATA3:              High (luminal identity maintained)
    PR:                 High
    CK8/18:             Retained
    p120-catenin:       Redistributes to cytoplasm when
                        E-cadherin is lost — diagnostic marker

INITIATING MOLECULAR EVENTS:
  CDH1 mutation/loss: ~65% of ILC (biallelic inactivation)
                    This is THE defining event of ILC.
                    No other cancer type has CDH1 loss as
                    the primary initiating event.
                    CDH1 loss creates the architectural
                    change (single-file invasion) AND
                    activates downstream survival signalling
                    through p120-catenin redistribution.
  PIK3CA:           ~48% of ILC (higher than IDC Luminal A)
  FOXA1 mutation:   ~14% (FOXA1 is the ER pioneer TF —
                    FOXA1 mutation alters ER binding
                    pattern and may contribute to
                    endocrine therapy resistance)
  ERBB2:            ~2% amplification (rare in ILC)
  TP53:             ~4% (extremely low — among the lowest
                    of any breast cancer subtype)
  TBX3, RUNX1:      Mutations enriched in ILC vs. IDC

KEY DISTINGUISHING FEATURES:
  - CDH1 loss = defining event (E-cadherin loss)
  - Single-file invasion architecture
  - p120-catenin cytoplasmic redistribution (diagnostic IHC)
  - Almost always ER+, rarely HER2+
  - Missed by mammography more often than IDC
  - Unusual metastatic sites (GI tract, peritoneum)
  - Lower pCR to chemotherapy
  - Higher late recurrence than IDC-Luminal A
  - FOXA1 mutations — may alter endocrine therapy response
  - The Waddington false attractor here is unique:
    The cell retains luminal identity (ER+, FOXA1+, GATA3+)
    but loses architectural coherence (E-cadherin gone).
    It is a false attractor characterized not by TF
    replacement but by cell adhesion programme dismantling.
    This is structurally distinct from all other subtypes.

AVAILABLE PUBLIC DATA:
  TCGA-BRCA:        n~180 ILC samples (histology field)
                    A dedicated TCGA analysis of ILC was
                    published (TCGA 2015 Nat Commun) —
                    ILC clinical and molecular data available
  GSE109169:        ILC-specific dataset
  GSE146558:        ILC-enriched
  METABRIC:         Includes ILC with histology annotation

ANALYSIS STATUS:    NOT YET STARTED
PLANNED DOCUMENT:   BRCA-S6a (ILC before-document)
PRIORITY NOTE:      ILC is structurally the most interesting
                    breast cancer subtype for this framework
                    because the false attractor is defined
                    by adhesion programme loss, not TF
                    replacement. CDH1 loss as the Waddington
                    saddle point event is a tractable and
                    novel framing that does not appear in
                    the existing attractor literature.
```

---

## SECTION V — THE METAPLASTIC SUBTYPE (RARE)

```
BRIEF NOTE — NOT A PRIMARY ANALYSIS TARGET

Metaplastic breast carcinoma (MBC):
  Prevalence:       <1% of all breast cancer
  ER/PR/HER2:       Almost always triple-negative
  Histology:        Conversion to non-glandular elements
                    (spindle cell, squamous, chondroid,
                    osseous, or mixed components)
                    The cancer cell has undergone EMT
                    and differentiated toward mesenchymal
                    lineages — the false attractor is
                    a completely different tissue type
  Prognosis:        Worse than conventional TNBC
  Cell of origin:   Disputed — likely luminal progenitor
                    undergoing extreme EMT
  Public data:      Very limited (n<50 in any single dataset)
  Analysis status:  Not planned at this time due to
                    insufficient data and extreme rarity.
                    Would be addressed as a special case
                    after all other subtypes are complete.
```

---

## SECTION VI — THE TWO CLASSIFICATION SYSTEMS RECONCILED

```
For this repository, the analysis order and primary
classification is as follows:

MOLECULAR SUBTYPE    HISTOLOGICAL OVERLAP        ANALYSIS SERIES
────────────────────────────────────────────────────────────────
Luminal A            IDC (mostly), ILC (~80% of ILC)  BRCA-S1
Luminal B            IDC (mostly), ILC (~15% of ILC)  BRCA-S2
HER2-enriched        IDC (almost exclusively)          BRCA-S3
Basal-like / TNBC    IDC (almost exclusively)          BRCA-S4
Claudin-low          IDC, rare ILC                     BRCA-S5
ILC (separate)       Lobular histology — any PAM50     BRCA-S6

ILC IS ANALYZED SEPARATELY despite being mostly
Luminal A/B because:
  1. Different cell of origin (lobular vs. ductal)
  2. Different Waddington landscape (CDH1-defined)
  3. Different clinical behavior (late recurrence,
     unusual metastases, low chemo response)
  4. Different normal identity programme (E-cadherin
     is an essential structural component, not a TF)

The Luminal A analysis will use IDC-Luminal A samples.
The ILC analysis will use histologically confirmed ILC.
Where samples overlap in TCGA-BRCA, they will be
assigned to the appropriate series by histology first,
then by PAM50.
```

---

## SECTION VII — CROSS-SUBTYPE STRUCTURAL PICTURE

```
Before any analysis runs, here is what is ALREADY KNOWN
structurally about how these subtypes relate.
This is literature, not prediction.

THE DIFFERENTIATION HIERARCHY MAPPING:

  Mammary stem cell
  └─ Luminal progenitor
     ├─ Mature luminal cell → Luminal A (least displaced)
     ├─ Intermediate luminal → Luminal B (partly displaced)
     ├─ Luminal progenitor (HER2-amplified) → HER2-enriched
     └─ Luminal progenitor (BRCA1-null, basal programme)
        → Basal-like / TNBC
        → Claudin-low (most displaced — nearest stem)
  └─ Lobular progenitor
     └─ Lobular epithelial cell → ILC (CDH1 loss)

  This is a SINGLE WADDINGTON LANDSCAPE for the breast
  with MULTIPLE FALSE ATTRACTORS accessible from
  different positions in the normal cell hierarchy.

  The key cross-subtype structural prediction
  (stated here as literature, not as a novel claim):
    The depth score for Luminal A should correlate with
    ESR1/FOXA1/GATA3 loss (the luminal identity programme).
    The depth score for TNBC should correlate with
    BRCA1/luminal progenitor identity loss and basal
    programme gain.
    The depth score for HER2-enriched should correlate
    with HER2 amplicon expression rising with depth.
    The depth score for ILC should correlate with
    CDH1 loss and p120-catenin redistribution.
    These are different axes in the same landscape.

THE ER-AXIS AS THE MASTER LUMINAL AXIS:
  ESR1 (ER) is the master luminal identity TF.
  It falls progressively across the luminal subtypes:
    Luminal A: highest ER expression
    Luminal B: intermediate ER expression
    HER2-enriched: ER negative (or very low)
    Basal-like: ER negative
  This is the main axis of the breast Waddington landscape —
  the ER axis. Displacement from ER expression = displacement
  from luminal identity.
  The cross-subtype analysis will test whether ER axis
  position predicts depth score across all luminal subtypes
  simultaneously.

THE BRCA1 CONNECTION:
  BRCA1 germline mutations strongly predispose to
  basal-like TNBC specifically (not luminal subtypes).
  This is because BRCA1 is required for LUMINAL PROGENITOR
  identity — without BRCA1, the luminal progenitor cannot
  complete differentiation and instead activates the
  basal-like programme.
  BRCA1 is therefore a gate between the normal luminal
  progenitor valley and the basal-like false attractor.
  When BRCA1 is lost, the gate opens.
  This is the most mechanistically precise Waddington
  crossing in all of breast cancer.
```

---

## SECTION VIII — DATA AVAILABILITY SUMMARY

```
Subtype        Primary Dataset     n (tumour)   Normal Ref        Power
────────────────────────────────────────────────────────────��───────────
Luminal A      TCGA-BRCA PAM50     ~500–600     TCGA adj normal   HIGH
Luminal B      TCGA-BRCA PAM50     ~200–250     TCGA adj normal   HIGH
HER2-enriched  TCGA-BRCA PAM50     ~80–100      TCGA adj normal   MOD
Basal-like     TCGA-BRCA PAM50     ~140–200     TCGA adj normal   HIGH
               + GSE25066 (pCR)    ~508         pre-treatment
Claudin-low    TCGA-BRCA subset    ~60–80       TCGA adj normal   LOW
ILC            TCGA-BRCA histology ~180         TCGA adj normal   MOD
               + GSE109169

NOTE: TCGA-BRCA contains n~1,100 total samples.
      PAM50 subtype calls are available in the clinical file.
      Adjacent normal samples (n~113) available for normal
      reference — the most important single dataset for
      normal identity programme definition.
```

---

## SECTION IX — THE PLANNED ANALYSIS ORDER

```
ORDER:

  BRCA-S1   Luminal A      TCGA-BRCA PAM50 LumA    HIGH POWER
                           REASON: Largest subtype. Highest
                           power. Establishes the luminal
                           identity baseline for all comparisons.
                           ESR1/FOXA1/GATA3 axis will anchor
                           all subsequent luminal analyses.

  BRCA-S2   TNBC           TCGA-BRCA PAM50 Basal    HIGH POWER
            (Basal-like)   + GSE25066 pCR cohort
                           REASON: Greatest clinical unmet need.
                           pCR data enables the most clinically
                           actionable depth score analysis
                           (does depth score predict pCR?).
                           Second priority after Luminal A.

  BRCA-S3   HER2-enriched  TCGA-BRCA PAM50 HER2     MODERATE
                           REASON: Most targeted-therapy-
                           transformed subtype. HER2 amplicon
                           as attractor switch gene is the
                           cleanest targetable false attractor
                           in all of breast cancer.

  BRCA-S4   ILC            TCGA-BRCA histology      MODERATE
                           + GSE109169
                           REASON: Structurally most novel
                           for this framework (CDH1-defined
                           false attractor). The most
                           scientifically interesting analysis.

  BRCA-S5   Luminal B      TCGA-BRCA PAM50 LumB     HIGH POWER
                           REASON: After Luminal A is complete,
                           Luminal B comparison answers:
                           what is the structural difference
                           between a terminally differentiated
                           luminal cell that escaped (LumA)
                           vs. a luminal progenitor that
                           never fully differentiated (LumB)?

  BRCA-S6   Claudin-low    TCGA-BRCA classifier     LOW POWER
                           REASON: Last because lowest n
                           and no specific clinical target.
                           Most speculative analysis.

CROSS-SUBTYPE ANALYSES:
  After BRCA-S1 and BRCA-S2:
    → Luminal vs. Basal: ER-axis as the master breast
      landscape axis

  After BRCA-S1 through BRCA-S5:
    → Full breast landscape cross-subtype analysis
      (analogous to RCC cross-type series)
    → Does depth score across all luminal subtypes
      predict the same drugs? Are there basket trial
      candidates across Luminal A + B + HER2-enriched?
```

---

## SECTION X — WHAT THIS DOCUMENT DOES NOT CONTAIN

```
This document contains:
  ✓ Subtype definitions (histological and molecular)
  ✓ Cells of origin (per subtype)
  ✓ Published molecular drivers (literature, not prediction)
  ✓ Clinical characteristics (survival, standard of care)
  ✓ Data availability and GEO/TCGA accessions
  ✓ The relationship between the two classification systems
  ✓ The planned analysis order with justification

This document does NOT contain:
  ✗ Depth score hypotheses
  ✗ False attractor gene predictions
  ✗ Drug target predictions
  �� Circuit topology predictions
  ✗ Epigenetic lock mechanisms
  ✗ Cross-subtype structural predictions beyond
    what is established literature

All of the above belong in the BEFORE documents
for each individual subtype analysis.
The before-document for Luminal A (BRCA-S1a) is the
next document in this series. It will be written
before any script runs and before any data is loaded.
```

---

## STATUS BLOCK

```
document:           BRCA_Subtype_Orientation.md
folder:             Cancer_Research/BRCA/Subtypes/
status:             COMPLETE — ORIENTATION ONLY
date:               2026-03-03
author:             Eric Robert Lawson
                    OrganismCore

subtypes_covered:
  molecular:        Luminal A, Luminal B, HER2-enriched,
                    Basal-like/TNBC, Claudin-low     [5 of 5]
  histological:     ILC (separate series)            [1 of 1]
  rare:             Metaplastic (noted, not planned)  [1 noted]

analyses_started:   0

next_document:      BRCA-S1a
                    Luminal A Before-Document
                    (predictions locked before
                    TCGA-BRCA PAM50 Luminal A loads)

series_rule:        No subtype analysis begins until
                    this document is read and acknowledged.
                    The order is Luminal A first.
                    No exceptions.

relationship_to_existing_BRCA_analysis:
                    The existing analysis in Cancer_Research/BRCA/
                    was conducted on TCGA-BRCA without
                    explicit PAM50 subtype stratification.
                    The analyses in this Subtypes/ folder
                    are subtype-stratified analyses using
                    the same dataset but with PAM50 calls
                    as the primary stratification variable.
                    They are not replacements for the
                    existing analysis — they are extensions
                    that go deeper into the subtype structure
                    the existing analysis identified.
```
