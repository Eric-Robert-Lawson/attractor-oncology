# OVARIAN CANCER — SUBTYPE ORIENTATION DOCUMENT
## Before Any Analysis Begins
## OrganismCore | 2026-03-03 | Author: Eric Robert Lawson

---

## PURPOSE OF THIS DOCUMENT

```
This document exists before any script runs.
It contains no predictions.
It contains no epigenetic hypotheses.
It contains no depth score derivations.

What it contains:

  A complete map of the ovarian cancer subtype landscape —
  what the subtypes are, where they come from anatomically,
  what the clinical and molecular literature has established
  about each one, and what public data exists to analyze them.

This document is the prerequisite to every analysis in this
folder. Any analyst — human or machine — should read this
document before touching any data in OV/.

It is the orientation. Not the analysis.
```

---

## DOCUMENT METADATA

```
document_id:        OV_Subtype_Orientation
series:             OV (Ovarian Cancer)
folder:             Cancer_Research/OV/
date:               2026-03-03
author:             Eric Robert Lawson
status:             ORIENTATION ONLY — no predictions, no analysis
next_document:      OV_HGSC_before.md (Document OV-1a)
protocol_version:   Workflow_Protocol.md v2.0
```

---

## SECTION I — THE FUNDAMENTAL PROBLEM WITH "OVARIAN CANCER"

```
"Ovarian cancer" is not one disease.

It is a clinical label applied to at least seven biologically
distinct malignancies that happen to present in or near the
ovary. They arise from different cells. They travel through
different molecular transitions. They respond to different
drugs. They kill patients through different mechanisms.

Treating them as one disease is the central historical error
in ovarian oncology.

The five-year survival rate for ovarian cancer overall is
approximately 49%. For high-grade serous carcinoma diagnosed
at Stage III/IV (the most common presentation), it is
approximately 29%. These numbers have improved modestly
over 30 years despite platinum-based chemotherapy remaining
the backbone of treatment.

The modest improvement is partly explained by the fact that
"ovarian cancer" clinical trials have predominantly enrolled
high-grade serous carcinoma patients, produced protocols
optimized for that subtype, and then applied those protocols
to all subtypes — including the ones with completely different
biology that do not respond to platinum at all.

The framework in this repository treats each subtype as a
separate disease with a separate normal cell of origin,
a separate Waddington landscape, and a separate false
attractor. The cross-subtype comparison comes after
each individual analysis is complete — not before.
```

---

## SECTION II — THE TWO MAJOR GROUPS

```
WHO Classification (5th Edition, 2020) divides ovarian
malignancies into two major groups:

GROUP 1 — EPITHELIAL OVARIAN CANCERS (EOC)
  ~90% of all ovarian malignancies.
  Arise from epithelial cells lining the reproductive tract.
  Five distinct subtypes with distinct cells of origin,
  distinct molecular programmes, and distinct clinical behavior.

GROUP 2 — NON-EPITHELIAL OVARIAN TUMORS
  ~10% of ovarian malignancies.
  Two major families: germ cell tumors and sex cord-stromal tumors.
  Completely different cell lineages from Group 1.
  Different patient demographics (often younger patients).
  Different treatment paradigms.

This repository will analyze Group 1 subtypes first,
as they represent the largest patient populations and
have the most available public RNA-seq data.

Group 2 will be addressed as separate analyses after
Group 1 is complete. They are not variants of EOC.
They are different diseases entirely.
```

---

## SECTION III — THE FIVE EPITHELIAL SUBTYPES

---

### SUBTYPE 1 — HIGH-GRADE SEROUS CARCINOMA (HGSC)

```
CLINICAL FACTS:
  Prevalence:       ~70% of all epithelial ovarian cancer
                    ~63% of all ovarian malignancies
  Stage at diagnosis: ~75% diagnosed at Stage III or IV
                      (peritoneal spread)
  5-year survival:  ~29% (Stage III/IV)
                    ~92% (Stage I — rare at diagnosis)
  Standard of care: Platinum + taxane chemotherapy
                    PARP inhibitors (BRCA-mutated patients)
                    Bevacizumab (anti-VEGF)
  Recurrence:       ~80% of patients recur within 2 years
                    Platinum-resistant recurrence is the
                    dominant cause of death

ANATOMICAL ORIGIN:
  Long debated. The prevailing model (Kurman & Shih, 2016)
  now places the cell of origin in the FALLOPIAN TUBE —
  specifically in the secretory cells of the fimbriated end
  (distal fallopian tube / fimbriae), not in the ovarian
  surface epithelium as was historically assumed.

  Evidence:
  - Serous tubal intraepithelial carcinoma (STIC) lesions
    found in the fallopian tubes of BRCA carriers and in
    the fallopian tubes adjacent to ovarian HGSC tumors
  - STIC and adjacent HGSC share identical TP53 mutations
    — establishing STIC as the precursor lesion
  - The secretory cell of the fallopian tube expresses
    the same PAX8, WT1, and ER profile as HGSC tumors
  - Prophylactic salpingo-oophorectomy in BRCA carriers
    finds STIC in fallopian tube in ~5–10% of specimens

NORMAL CELL IDENTITY:
  Fallopian tube secretory cell
  Key identity markers:
    PAX8       — transcription factor, secretory cell identity
    WT1        — Wilms tumor protein, serous identity marker
    OVGP1      — oviduct-specific glycoprotein 1
    CAPS       — capping protein (ciliated cell marker)
    FOXJ1      — ciliated cell identity (FOXJ1-negative = secretory)
    ESR1       — estrogen receptor, expressed in secretory cells
  The fallopian tube secretory cell is a quiescent,
  hormone-responsive, PAX8+ cell that does not normally
  proliferate autonomously.

INITIATING MOLECULAR EVENTS:
  TP53 mutation:    >95% of HGSC cases
                    This is not a late event — it is the
                    initiating event. TP53 loss permits the
                    secretory cell to survive DNA damage
                    without apoptosis and to replicate
                    with genomic instability.
  BRCA1/BRCA2:      ~25% of HGSC patients carry germline or
                    somatic BRCA1/2 mutations
                    BRCA loss impairs homologous recombination
                    — explains platinum/PARP inhibitor
                    sensitivity in this subgroup
  Genomic instability: HGSC has the highest somatic copy
                    number alteration burden of any cancer
                    type in TCGA. Few recurrent point
                    mutations beyond TP53. The driver of
                    HGSC is chromosomal, not point-mutational.

KEY DISTINGUISHING FEATURES:
  - TP53 mutated in essentially every case
  - Wild-type KRAS, BRAF (unlike LGSC)
  - Wild-type ARID1A (unlike clear cell and endometrioid)
  - Platinum-sensitive initially (~80% response rate)
  - Platinum-resistant at recurrence in majority of patients
  - BRCA-mutated patients: PARP inhibitor responsive
  - BRCA wild-type patients: fewer targeted options

AVAILABLE PUBLIC DATA:
  TCGA-OV:          n=489 tumour samples (RNA-seq + microarray)
                    TCGA's ovarian cancer dataset is overwhelmingly
                    HGSC — this was by design. TCGA-OV is the
                    reference dataset for HGSC.
                    GDC portal: phs000178 / TCGA-OV
                    UCSC Xena: available in harmonized format
  GEO datasets:
    GSE26712:       n=195 tumours, Affymetrix, subtyped
                    Includes HGSC, endometrioid, clear cell,
                    mucinous. Largest multi-subtype microarray.
    GSE18520:       n=53 tumours + n=10 normal fallopian tube
                    Affymetrix, includes normal reference
    GSE9891:        n=285 tumours, microarray, HGSC enriched
    GSE69428:       RNA-seq, serous subtype
    GSE32062:       n=260 tumours, survival annotated

ANALYSIS STATUS:    NOT YET STARTED
PLANNED DOCUMENT:   OV-1a (HGSC, Document 1 in OV series)
PLANNED DATASET:    TCGA-OV (primary) + GSE18520 (normal reference)
```

---

### SUBTYPE 2 — LOW-GRADE SEROUS CARCINOMA (LGSC)

```
CLINICAL FACTS:
  Prevalence:       ~5% of epithelial ovarian cancer
  Stage at diagnosis: Often diagnosed at advanced stage but
                      slower growing than HGSC
  5-year survival:  ~70–80% (better than HGSC despite
                    advanced stage — reflects slower biology)
  Standard of care: Surgery (critical — often the most
                    effective intervention)
                    Platinum + taxane (BUT: LGSC is largely
                    PLATINUM-RESISTANT — response rate ~3%)
                    MEK inhibitors (binimetinib, trametinib)
                    in development / approved in some settings
                    Hormonal therapy (letrozole, tamoxifen)
                    — LGSC is often hormone-sensitive
  Key insight:      LGSC is treated as "ovarian cancer" and
                    receives platinum-based chemotherapy
                    that does not work in this subtype.
                    This is one of the most clear-cut
                    examples of subtype misclassification
                    harming patients.

ANATOMICAL ORIGIN:
  Two competing models:
  1. Independent development from ovarian surface epithelium
     or inclusion cysts via adenoma → borderline tumor →
     invasive carcinoma sequence (the "Type I" pathway)
  2. Serous tubal intraepithelial carcinoma precursor
     (less well established than in HGSC)
  The majority evidence supports the adenoma → borderline
  → LGSC progression sequence. LGSC is a SLOW TRANSFORMER
  — it takes years to decades to develop from a borderline
  serous tumor.

NORMAL CELL IDENTITY:
  Ovarian surface epithelium or inclusion cyst epithelium
  Key identity markers:
    OSE markers: calretinin, CA125 (low), vimentin
    LGSC retains more differentiated serous features
    than HGSC — it is a less committed false attractor

INITIATING MOLECULAR EVENTS:
  KRAS mutation:    ~30% of LGSC
  BRAF mutation:    ~5% of LGSC
  NRAS mutation:    ~10% of LGSC
  (KRAS/BRAF/NRAS mutations are mutually exclusive)
  The MAPK/ERK pathway is the dominant driver in LGSC.
  This is fundamentally different from HGSC.
  Wild-type TP53 in >95% of cases.
  Low genomic instability (unlike HGSC).
  ERBB2 (HER2) amplification in a subset.

KEY DISTINGUISHING FEATURES:
  - Wild-type TP53 (HGSC has mutated TP53 — this is a
    definitive diagnostic distinction)
  - MAPK pathway mutations (KRAS/BRAF/NRAS)
  - Low genomic instability
  - Platinum-RESISTANT
  - MEK inhibitor-SENSITIVE (targetable)
  - Hormone receptor-positive (ER/PR) — important for
    hormonal therapy options
  - Slow growing — indolent course but incurable at
    advanced stage

AVAILABLE PUBLIC DATA:
  LGSC is rare. Public datasets are small.
  GSE73614:         LGSC samples included (subtyped)
  GSE6008:          Contains LGSC samples
  MD Anderson LGSC series: available via GEO (multiple
                    small cohorts, typically n<50 for LGSC)
  Note: Sample sizes for LGSC in any single dataset
        are typically n=20–60. Low-power analysis will
        apply. Same approach as cdRCC.

ANALYSIS STATUS:    NOT YET STARTED
PLANNED DOCUMENT:   OV-2a (LGSC, Document 2 in OV series)
PLANNED DATASET:    GSE73614 primary + GSE6008 supplementary
                    LOW POWER — results directional only
```

---

### SUBTYPE 3 — ENDOMETRIOID CARCINOMA (EnOC)

```
CLINICAL FACTS:
  Prevalence:       ~10% of epithelial ovarian cancer
  Stage at diagnosis: More often diagnosed at Stage I/II
                      than HGSC — better prognosis in part
                      because it is caught earlier
  5-year survival:  ~70–80% (Stage I/II); ~40–50% (Stage III)
  Standard of care: Platinum + taxane (responds better
                    than HGSC at advanced stage in some series)
                    Hormonal therapy (ER+/PR+ cases)
                    MSI-H cases: pembrolizumab (checkpoint)
  Association:      ~30% of EnOC co-occurs with
                    endometrial carcinoma (synchronous tumors)
                    ~15–20% arise in a background of
                    ovarian endometriosis

ANATOMICAL ORIGIN:
  Endometriosis — ectopic endometrial glands and stroma
  implanted on the ovarian surface or cortex.
  The cell of origin is the ENDOMETRIAL GLAND CELL that
  has implanted outside the uterus.
  This is the same cell type that gives rise to endometrial
  carcinoma — which is why EnOC and endometrial carcinoma
  share molecular features and sometimes occur synchronously.

  The transformation sequence:
  Endometriosis → atypical endometriosis → EnOC
  This sequence is well documented histologically.

NORMAL CELL IDENTITY:
  Endometrial gland cell (secretory phase epithelium)
  Key identity markers:
    FOXA2      — endometrial gland identity TF
    PAX2       — endometrial epithelial marker
    Vimentin   — endometrial stroma marker
    PR         — progesterone receptor (high in normal
                 endometrial glands)
    ER         — estrogen receptor

INITIATING MOLECULAR EVENTS:
  ARID1A mutation:  ~30–40% of EnOC
                    ARID1A is a SWI/SNF chromatin remodeling
                    subunit. Loss of ARID1A impairs nucleosome
                    repositioning — this is an epigenetic
                    event at the chromatin remodeling level,
                    NOT a DNA methylation event.
  PTEN mutation:    ~20% of EnOC
                    PI3K/Akt pathway activation
  PIK3CA mutation:  ~20% of EnOC
  CTNNB1 mutation:  ~16% of EnOC (β-catenin, Wnt pathway)
  MSI-H:            ~20% of EnOC (mismatch repair deficiency)
  Wild-type TP53 in low-grade; TP53 mutations in high-grade
  variant (high-grade EnOC overlaps with HGSC molecularly)

KEY DISTINGUISHING FEATURES:
  - ARID1A loss — the defining molecular feature
  - PTEN / PIK3CA co-mutations — PI3K pathway driven
  - Wnt/β-catenin active in CTNNB1-mutated cases
  - MSI-H in subset — pembrolizumab-eligible
  - ER/PR positive in most cases
  - Often co-occurs with synchronous endometrial carcinoma
  - BETTER prognosis than HGSC stage-for-stage
  - May respond to hormonal therapy

AVAILABLE PUBLIC DATA:
  GSE26712:         EnOC samples (subtyped, Affymetrix)
  GSE6008:          EnOC samples included
  TCGA endometrial: (UCEC) — not OV, but contains the same
                    cell of origin and shares molecular features
                    Useful as normal reference for EnOC

ANALYSIS STATUS:    NOT YET STARTED
PLANNED DOCUMENT:   OV-3a (EnOC, Document 3 in OV series)
```

---

### SUBTYPE 4 — CLEAR CELL CARCINOMA (OCCC)

```
CLINICAL FACTS:
  Prevalence:       ~10% of epithelial ovarian cancer globally
                    Higher prevalence in East Asian populations
                    (~25% of ovarian cancer in Japan)
  Stage at diagnosis: More often Stage I than HGSC
  5-year survival:  ~85% (Stage I) — EXCELLENT at Stage I
                    ~25% (Stage III/IV) — POOR at advanced stage
                    Why? OCCC is platinum-RESISTANT.
                    Stage I patients are cured by surgery.
                    Advanced stage patients have almost no
                    effective systemic treatment.
  Standard of care: Surgery (critical)
                    Platinum + taxane (low response rate ~15–20%)
                    No standard second-line option
                    mTOR inhibitors (temsirolimus) — modest activity
                    PD-1/PD-L1 inhibitors — being studied
                    (TMB-high / MSI-H subset eligible)
  This is one of the most treatment-refractory subtypes
  in all of oncology. An effective therapy for advanced
  OCCC does not currently exist.

ANATOMICAL ORIGIN:
  Same as EnOC: ENDOMETRIOSIS.
  ~50% of OCCC arise in a background of endometriosis.
  The same endometrial gland cell that gives rise to EnOC
  appears to give rise to OCCC via a different molecular
  transition — a different path from the same starting valley
  in the Waddington landscape.
  This is structurally analogous to ccRCC and PRCC
  sharing a proximal tubule origin but diverging at
  different saddle points.

NORMAL CELL IDENTITY:
  Endometrial gland cell (same as EnOC)
  BUT: OCCC diverges from EnOC early in its Waddington
  transition and acquires a distinct identity programme:
    HNF-1β overexpression — the defining marker of OCCC
                            (transcription factor driving
                            glycogen accumulation — the "clear
                            cell" appearance comes from
                            glycogen-laden cytoplasm)
    NAPSA              — marker of OCCC
    AMACR             — expressed in OCCC
    The OCCC false attractor is characterized by
    HNF-1β dominance — a completely different TF programme
    from the HGSC (WT1/PAX8) or EnOC (FOXA2/PAX2) attractor

INITIATING MOLECULAR EVENTS:
  ARID1A mutation:  ~50–60% of OCCC (HIGHEST of all subtypes)
                    ARID1A is the most commonly mutated gene
                    in OCCC. Same gene as EnOC — but higher
                    frequency and different co-mutations
                    leading to a different outcome.
  PIK3CA mutation:  ~33% of OCCC
  PTEN mutation:    ~8% (lower than EnOC)
  Wild-type TP53 in most cases
  HNF-1β overexpression: nearly universal
                    (mostly epigenetic/transcriptional,
                    not mutation-driven)
  Low genomic instability (unlike HGSC)
  MSI-H: rare (~2%)

KEY DISTINGUISHING FEATURES:
  - HNF-1β overexpression (the defining clinical marker)
  - ARID1A loss (highest frequency of all subtypes)
  - PIK3CA mutations (PI3K/mTOR pathway)
  - Platinum-RESISTANT (like LGSC, not like HGSC)
  - Clear cell morphology (glycogen-laden cytoplasm)
  - Higher prevalence in East Asian populations
  - Excellent Stage I prognosis; dismal Stage III/IV prognosis
  - No effective treatment at advanced stage currently
  - This is arguably the subtype most in need of a novel
    treatment framework — greatest unmet need

AVAILABLE PUBLIC DATA:
  GSE73614:         OCCC samples (subtyped)
  GSE6008:          OCCC samples
  GSE29450:         OCCC-enriched dataset
  Japanese cohorts: Multiple GEO datasets with OCCC
                    enrichment due to higher East Asian
                    prevalence (search GSE datasets from
                    Keio, Osaka, Tokyo institutions)
  Note: OCCC sample sizes in individual datasets
        are typically n=30–80. Moderate power.

ANALYSIS STATUS:    NOT YET STARTED
PLANNED DOCUMENT:   OV-4a (OCCC, Document 4 in OV series)
PRIORITY NOTE:      OCCC has the greatest unmet clinical need
                    of all epithelial ovarian subtypes.
                    It is the second priority after HGSC.
```

---

### SUBTYPE 5 — MUCINOUS CARCINOMA (MOC)

```
CLINICAL FACTS:
  Prevalence:       ~3% of epithelial ovarian cancer
                    (true primary MOC is rarer than reported —
                    many "ovarian mucinous carcinomas" are
                    actually metastases from the GI tract,
                    pancreas, or appendix masquerading as
                    primary ovarian tumors)
  Stage at diagnosis: ~75% diagnosed at Stage I
  5-year survival:  ~90%+ (Stage I — very good)
                    ~10–20% (Stage III/IV — very poor)
  Standard of care: Surgery (Stage I: surgery alone may cure)
                    Platinum + taxane (low response rate)
                    5-FU / capecitabine (GI-type regimens
                    — used in advanced stage, modest activity)
                    HER2-targeted therapy (subset with HER2+)
  Important caveat: DIAGNOSIS IS DIFFICULT.
                    A mucinous tumor in the ovary is more
                    likely to be a GI/appendix/pancreas
                    metastasis than a primary ovarian cancer.
                    Any analysis of MOC must use datasets
                    that have rigorously excluded metastatic
                    disease.

ANATOMICAL ORIGIN:
  DEBATED. This is the most uncertain of all subtypes.
  Three candidate origins:
  1. Walthard rests — transitional cell nests at the
     tubo-ovarian junction. These express mucinous markers
     and are the current leading hypothesis.
  2. Mucinous cystadenoma → borderline mucinous tumor →
     MOC progression sequence (the "Type I" pathway,
     analogous to LGSC)
  3. Teratomatous origin — from a mature cystic teratoma
     containing gastrointestinal epithelium
  The true cell of origin is unknown. This is a limitation
  for Waddington landscape analysis — if the normal cell
  is uncertain, the normal identity programme is uncertain.

NORMAL CELL IDENTITY:
  Uncertain — most likely intestinal-type mucinous epithelium
  Key features of the normal cell (if Walthard rest origin):
    CDX2       — intestinal transcription factor
    MUC2, MUC5B — mucin proteins
    CK20       — intestinal epithelial marker
  The MOC false attractor identity looks like intestinal
  mucosa — which may reflect normal cell of origin or
  may reflect the false attractor programming.

INITIATING MOLECULAR EVENTS:
  KRAS mutation:    ~50% of MOC
                    Same KRAS mutation as pancreatic and
                    colorectal cancer — supports the
                    intestinal-type origin hypothesis
  HER2 amplification: ~18% of MOC (targetable)
  BRAF mutation:    rare in MOC (unlike LGSC)
  TP53 mutation:    present in high-grade MOC
  Wild-type ARID1A, BRCA1/2

KEY DISTINGUISHING FEATURES:
  - KRAS-driven (like LGSC but different subtype entirely)
  - HER2 amplification in subset (targetable with
    trastuzumab — active clinical trials)
  - Platinum-RESISTANT
  - Often mimics GI cancers molecularly
  - Stage I: surgery may be curative
  - Stage III/IV: no effective treatment
  - Uncertain cell of origin — limits depth score derivation
  - Smallest true prevalence of all epithelial subtypes

AVAILABLE PUBLIC DATA:
  Very limited for true primary MOC.
  GSE26712:         Includes mucinous samples but small n
  GSE6008:          Includes mucinous samples
  Note: MOC datasets are small and contaminated by
        metastatic cases in older series. Any MOC analysis
        requires careful dataset curation.
  ANALYSIS PRIORITY: LOWEST among the five epithelial
                     subtypes, due to uncertain cell of
                     origin and small dataset sizes.

ANALYSIS STATUS:    NOT YET STARTED
PLANNED DOCUMENT:   OV-5a (MOC, Document 5 in OV series)
                    AFTER all other subtypes are complete.
```

---

## SECTION IV — THE TWO NON-EPITHELIAL GROUPS

---

### GROUP 2A — GERM CELL TUMORS

```
CLINICAL FACTS:
  Prevalence:       ~3–5% of all ovarian malignancies
  Patient demographic: Predominantly young women (median
                    age ~19 years). Most common ovarian
                    malignancy in women under 20.
  Key types:
    Dysgerminoma:   Most common germ cell tumor (~50% of
                    malignant GCT). Highly chemosensitive.
                    Analogous to testicular seminoma.
    Yolk sac tumor: Aggressive. Produces AFP.
    Choriocarcinoma: Extremely rare. Produces hCG.
    Immature teratoma: Graded I–III by neuroectodermal content
    Mixed GCT:      Contains multiple germ cell elements
  5-year survival:  ~90%+ (dysgerminoma, BEP chemotherapy)
                    ~75%+ (non-dysgerminoma, BEP)
  Standard of care: BEP chemotherapy (bleomycin, etoposide,
                    cisplatin) — highly effective
                    Fertility-sparing surgery where possible

CELL OF ORIGIN:
  Primordial germ cell — the embryonic precursor that
  migrates from the yolk sac to the developing gonad.
  These are the most primitive cells in the body.
  Their normal programme is pluripotency suppression
  and meiotic preparation. When that programme fails,
  they revert toward pluripotent or extraembryonic
  identity — which is what the different GCT subtypes
  represent (embryonal, yolk sac, trophoblastic, etc.)

KEY MOLECULAR FEATURES:
  KIT mutations:    ~25% of dysgerminoma
  i(12p):           Isochromosome 12p (also seen in
                    testicular GCT)
  OCT4 / NANOG:     Expressed in dysgerminoma (pluripotency
                    markers — the germ cell has reverted
                    toward pluripotent identity)
  DICER1:           Some GCT subtypes

WADDINGTON RELEVANCE:
  Germ cell tumors are the most direct expression of
  attractor biology in cancer:
    Normal germ cell → suppresses pluripotency → meiosis
    GCT → pluripotency reactivated → proliferative identity
  The "false attractor" in a dysgerminoma is an embryonic
  stem cell-like state. The depth score would measure
  re-expression of pluripotency markers relative to the
  meiotic preparation programme.
  This is a tractable analysis — but a separate series.

ANALYSIS STATUS:    NOT YET STARTED
PLANNED SERIES:     OV-GCT (separate from EOC series)
DATA AVAILABILITY:  GEO datasets available but small
                    (GCT is rare — most cohorts n<30)
```

---

### GROUP 2B — SEX CORD-STROMAL TUMORS

```
CLINICAL FACTS:
  Prevalence:       ~7% of all ovarian malignancies
  Key types:
    Adult granulosa cell tumor (AGCT): Most common.
                    Can produce estrogen. Endometrial
                    hyperplasia in 25–50% of patients.
    Juvenile granulosa cell tumor (JGCT): Rare, young patients
    Sertoli-Leydig cell tumor: Rare. Can produce androgens.
    Fibroma/thecoma: Benign or low-grade, rarely malignant
  5-year survival:  ~90% (Stage I AGCT)
                    ~50% (Stage III — AGCT is an indolent
                    recurrent disease; recurrences may
                    occur 10–20 years after initial treatment)
  Standard of care: Surgery
                    BEP or platinum-based chemotherapy
                    for advanced/recurrent disease

CELL OF ORIGIN:
  Gonadal stromal cells and sex cord precursors derived
  from embryonic gonadal mesenchyme — specifically the
  cells that give rise to granulosa cells (follicle
  supportive cells) and theca cells (androgen-producing
  stromal cells).
  These are support cells, not germ cells and not
  surface epithelium. Completely different lineage
  from all Group 1 subtypes.

KEY MOLECULAR FEATURES:
  FOXL2 C134W:      ~97% of adult granulosa cell tumors
                    This single point mutation is essentially
                    pathognomonic for AGCT. FOXL2 is the
                    master transcription factor for granulosa
                    cell identity. C134W is a gain-of-function
                    mutation that keeps granulosa cells in a
                    proliferative state.
  DICER1:           ~60% of Sertoli-Leydig cell tumors
                    (DICER1 is an RNA processing enzyme —
                    its loss disrupts miRNA biogenesis)
  SMAD3:            Mutations in some AGCT

WADDINGTON RELEVANCE:
  AGCT is a compelling attractor story:
    Normal granulosa cell: FOXL2 wild-type, follicle
                    supportive, non-proliferative, responds
                    to FSH signaling to support oocyte development
    AGCT false attractor: FOXL2 C134W locks the cell into
                    a proliferative granulosa identity that
                    does not respond to normal cycle signals
  The depth score would measure displacement from the
  normal FSH-responsive granulosa programme toward the
  constitutively proliferative FOXL2 C134W-driven state.
  This is a tractable analysis — but requires the right
  dataset. AGCT is rare and datasets are small.

ANALYSIS STATUS:    NOT YET STARTED
PLANNED SERIES:     OV-SCST (separate from EOC series)
DATA AVAILABILITY:  Very limited. Most AGCT cohorts n<30.
```

---

## SECTION V — THE CROSS-SUBTYPE STRUCTURAL PICTURE

```
Before any analysis runs, here is what is ALREADY KNOWN
structurally about how these subtypes relate to each other.
This is literature, not prediction.

THE TYPE I / TYPE II DIVISION (Kurman & Shih model):

  TYPE I CANCERS (slow, stepwise, indolent):
    LGSC, EnOC, OCCC, MOC
    Characteristics:
    - Arise from precursor lesions (borderline tumors,
      endometriosis, cystadenomas)
    - Low genomic instability
    - Specific driver mutations (KRAS, BRAF, ARID1A, PIK3CA)
    - Wild-type TP53 in majority
    - Diagnosed at lower stage more often
    - Better prognosis stage-for-stage

  TYPE II CANCERS (rapid, genomically unstable, aggressive):
    HGSC (and high-grade EnOC — disputed)
    Characteristics:
    - Arise rapidly, often without identifiable precursor
    - Extreme genomic instability
    - TP53 mutation nearly universal
    - BRCA1/2 mutations in significant subset
    - Diagnosed at advanced stage in majority
    - Worse prognosis

  THIS DIVISION IS THE FIRST CROSS-SUBTYPE STRUCTURAL FACT.
  It maps onto the Waddington framework:
    Type I: slow Waddington crossing — the cell traverses
            a gradual slope from normal to false attractor
            over years to decades. Precursor lesions are
            the cells sitting on the slope.
    Type II: rapid Waddington crossing — the cell falls
             off a cliff, not a slope. STIC → HGSC in months.
             No visible slope (no precursor that persists).

THE ENDOMETRIOSIS CONNECTION (EnOC and OCCC):

  EnOC and OCCC share a cell of origin (endometrial gland
  cell) but diverge at different saddle points.
  This is structurally identical to ccRCC and PRCC sharing
  a proximal tubule origin but diverging at different
  Waddington saddle points.
  The cross-subtype comparison of EnOC and OCCC is
  therefore the direct analogue of the ccRCC/PRCC
  cross-subtype comparison within RCC.
  They start in the same valley. They end in different
  false attractors. What made them diverge is the question.

THE SEROUS SPLIT (HGSC and LGSC):

  HGSC and LGSC are both "serous" — they produce serous
  fluid and express serous markers (CA125, WT1 in some).
  But they are not variants of the same disease.
  They arise from different cells (fallopian tube vs.
  ovarian surface/inclusion cyst), have different driver
  mutations (TP53 vs. KRAS/BRAF), and have opposite
  responses to platinum chemotherapy.
  Calling them both "serous" is a histological
  classification that does not reflect molecular reality.
  In the Waddington framework, they have different
  normal cells and different false attractors.
  They will be analyzed as entirely separate diseases.

THE MUCINOUS PROBLEM:

  MOC is the most difficult subtype for this framework.
  The cell of origin is uncertain. If we do not know the
  normal cell, we cannot define the normal identity
  programme, which means we cannot define the axis from
  normal to false attractor, which means the depth score
  cannot be derived.
  MOC will be the last subtype analyzed, after the cell
  of origin question has been addressed in the before-document.
```

---

## SECTION VI — DATA AVAILABILITY SUMMARY

```
Subtype     Primary Dataset      n (tumour)  Normal Ref        Power
─────────────────────────────────────────────────────────────────────
HGSC        TCGA-OV              ~489        GSE18520 FT       HIGH
LGSC        GSE73614 + GSE6008   ~30–60      OSE cell lines    LOW
EnOC        GSE26712 + GSE6008   ~40–60      TCGA-UCEC proxy   LOW-MOD
OCCC        GSE73614 + GSE29450  ~40–80      Endometriosis     MOD
MOC         GSE26712             ~20–40      Unknown           LOW
GCT         GEO series           <30         Primordial GC     VERY LOW
SCST        GEO series           <30         Granulosa cell    VERY LOW

NOTE: Only HGSC has HIGH power from a single dataset.
      All other subtypes will carry the LOW POWER caveat.
      Results from low-power subtypes are directional only.
      This is not a reason to skip them — it is a reason
      to state the caveat explicitly in every document.
```

---

## SECTION VII — THE PLANNED ANALYSIS ORDER

```
The following order is determined by:
  1. Data quality and sample size (highest power first)
  2. Clinical unmet need (highest unmet need prioritized)
  3. Structural relationship (related subtypes grouped)

ORDER:

  OV-1   HGSC    High-grade serous    TCGA-OV n=489    HIGH POWER
                 REASON: Largest dataset. Most common.
                 Established the baseline for all comparisons.

  OV-2   OCCC    Clear cell           GSE datasets     MODERATE
                 REASON: Greatest clinical unmet need.
                 Platinum-resistant. No effective treatment.
                 ARID1A makes it the most epigenetically
                 tractable of the rarer subtypes.

  OV-3   EnOC    Endometrioid         GSE datasets     LOW-MOD
                 REASON: Shares cell of origin with OCCC.
                 Cross-subtype comparison with OCCC will
                 answer: why do cells from the same starting
                 valley reach different false attractors?

  OV-4   LGSC    Low-grade serous     GSE datasets     LOW
                 REASON: MEK-driven, platinum-resistant,
                 distinct biology from HGSC. The contrast
                 with HGSC across the serous classification
                 is instructive.

  OV-5   MOC     Mucinous             GSE datasets     LOW
                 REASON: Uncertain origin. Last because
                 the before-document must address the
                 cell of origin problem before analysis.

  OV-GCT  Germ cell tumors           GEO datasets     VERY LOW
                 REASON: Different lineage entirely.
                 Separate analysis series.

  OV-SCST Sex cord-stromal            GEO datasets     VERY LOW
                 REASON: Different lineage entirely.
                 Separate analysis series.

CROSS-SUBTYPE ANALYSIS:
  After OV-1 through OV-3 are complete:
  → OV Cross-Type Script 1: What do HGSC, OCCC, EnOC
    share structurally? (analogous to RCC cross-type series)

  After OV-4 and OV-5:
  → OV Cross-Type Script 2: Serous division structural
    comparison (HGSC vs. LGSC)
```

---

## SECTION VIII — WHAT THIS DOCUMENT DOES NOT CONTAIN

```
This document contains:
  ✓ Subtype definitions
  ✓ Cells of origin
  ✓ Published molecular drivers (literature, not prediction)
  ✓ Clinical characteristics (survival, standard of care)
  ✓ Data availability
  ✓ Planned analysis order

This document does NOT contain:
  ✗ Depth score hypotheses
  ✗ False attractor gene predictions
  ✗ Drug target predictions
  ✗ Circuit topology hypotheses
  ✗ Epigenetic lock mechanisms
  ✗ Cross-subtype structural predictions

All of the above belong in the BEFORE documents
for each individual subtype analysis.
The before-document for HGSC (OV-1a) is the next document.
It will be written before any script runs and before
any data is loaded.
```

---

## STATUS BLOCK

```
document:           OV_Subtype_Orientation.md
status:             COMPLETE — ORIENTATION ONLY
date:               2026-03-03
author:             Eric Robert Lawson
                    OrganismCore

subtypes_covered:
  epithelial:       HGSC, LGSC, EnOC, OCCC, MOC   [5 of 5]
  non_epithelial:   GCT, SCST                       [2 of 2]

analyses_started:   0
next_document:      OV-1a — HGSC Before-Document
                    (predictions locked before TCGA-OV loads)

series_rule:        No subtype analysis begins until
                    this document is read and acknowledged.
                    The order is HGSC first.
                    No exceptions.
```
