# GLIOBLASTOMA AND DIFFUSE GLIOMA — SUBTYPE ORIENTATION DOCUMENT
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

  A complete map of the glioblastoma and diffuse glioma
  subtype landscape — what the WHO 2021 classification
  established, why it changed what "GBM" means, what
  the TCGA molecular subtypes within true GBM are,
  what the cells of origin are for each entity, what
  public data exists, and what the existing analysis
  in this repository captured.

GBM has a classification problem that is unique in this
repository: the name of the disease CHANGED in 2021.

What was called "GBM" before 2021 included both:
  1. True glioblastoma (IDH-wildtype) — the rapidly lethal
     primary brain cancer most people think of when they
     hear "GBM"
  2. Secondary glioblastoma (IDH-mutant) — a different
     disease that progressed slowly from a lower-grade
     glioma and only looked like GBM histologically

WHO 2021 reclassified (2) as "Astrocytoma, IDH-mutant,
grade 4" and removed it from the GBM category.

The word "GBM" now means ONLY IDH-wildtype glioblastoma.

The TCGA-GBM dataset predates this reclassification.
It contains both IDH-wildtype and IDH-mutant tumors.

The existing GBM analysis in Cancer_Research/GBM/ must
be interpreted in this context — and the subtype analyses
in this folder will decompose it accordingly.
```

---

## DOCUMENT METADATA

```
document_id:        GBM_Subtype_Orientation
series:             GBM (Glioblastoma and Diffuse Glioma — Subtypes)
folder:             Cancer_Research/GBM/Subtypes/
date:               2026-03-03
author:             Eric Robert Lawson
status:             ORIENTATION ONLY — no predictions,
                    no analysis, no epigenetic hypotheses
next_document:      GBM_IDHwt_before.md
                    (Document GBM-S1a — IDH-wildtype GBM
                    before-document)
protocol_version:   Workflow_Protocol.md v2.0
```

---

## SECTION I — THE CLASSIFICATION PROBLEM: WHAT CHANGED IN 2021

```
The WHO 2021 CNS Tumor Classification (5th edition)
made a single change that restructured the entire field
of brain tumor oncology:

  BEFORE 2021:
    "Glioblastoma" = any grade IV diffuse glioma
    This included IDH-mutant tumors that had progressed
    from lower-grade astrocytomas (historically called
    "secondary GBM") as well as the de novo IDH-wildtype
    primary tumors (historically called "primary GBM").
    Prognosis and biology were treated as similar
    because the histology looked similar.

  AFTER 2021 (WHO 5th Edition):
    "Glioblastoma, IDH-wildtype" = a single biological
    entity. The de novo high-grade glioma arising without
    a lower-grade precursor. Almost always in adults
    over 55. Uniformly lethal. Median OS ~15–21 months.

    "Astrocytoma, IDH-mutant, grade 4" = the former
    "secondary GBM." Now a separate entity. Different
    cell of origin, different molecular drivers, better
    prognosis (median OS ~3–4 years), different treatment
    implications.

    "Oligodendroglioma, IDH-mutant and 1p/19q-codeleted"
    = previously a subtype of GBM in some classifications.
    Now a completely separate entity.

THE PRACTICAL CONSEQUENCE FOR THIS REPOSITORY:

  TCGA-GBM was collected before 2021.
  It contains approximately:
    ~80% IDH-wildtype (true GBM by 2021 definition)
    ~20% IDH-mutant (now classified as astrocytoma
         grade 4 or other IDH-mutant glioma)

  The existing analysis in Cancer_Research/GBM/ ran on
  the full TCGA-GBM dataset, likely including both.
  This does not invalidate the analysis — but the
  depth axis found reflects a mixture of two biologically
  distinct Waddington landscapes.

  The subtype analyses in this folder will:
    1. Analyze IDH-wildtype GBM as the primary entity
       (what "GBM" means post-2021)
    2. Analyze IDH-mutant astrocytoma as a separate
       entity (completely different Waddington landscape)
    3. Address oligodendroglioma as a third entity
       (yet another separate landscape, with 1p/19q
       codeletion as the defining molecular event)

  The three entities will then be compared in a
  cross-subtype analysis after individual analyses
  are complete.
```

---

## SECTION II — THE NORMAL BRAIN CELL HIERARCHY

```
Before defining false attractors, the normal brain cell
hierarchy must be understood — the Waddington landscape
of the normal brain.

THE NORMAL GLIAL LINEAGE:

  Neural stem cell (NSC)
  ├─ Astrocyte progenitor → Mature astrocyte
  │    Function: glutamate uptake, BBB maintenance,
  │    synapse support, water homeostasis
  │    Identity: GFAP+, SOX9+, ALDH1L1+, AQP4+
  │    Non-proliferating in adult brain
  │
  ├─ Oligodendrocyte precursor cell (OPC)
  │    Identity: PDGFRA+, SOX10+, OLIG2+, NG2+
  │    Function: progenitor pool — produces
  │    oligodendrocytes that myelinate axons
  │    Slowly cycling in adult brain
  │    → Oligodendrocyte: MBP+, PLP1+, MOG+
  │
  └─ Ependymal cell (lines ventricles)
       Not directly relevant to GBM

  KEY IDENTITY MARKERS BY CELL TYPE:
    Astrocyte:      GFAP, AQP4, ALDH1L1, SOX9, S100B
    OPC:            PDGFRA, OLIG2, SOX10, NG2 (CSPG4)
    Oligodendrocyte: MBP, PLP1, MOG, OLIG1
    NSC/stem:       NES (nestin), SOX2, CD133 (PROM1)
    Neurons:        MAP2, TUBB3, RBFOX3 (NeuN)

  WHY THIS MATTERS FOR THE FRAMEWORK:
    The three main glioma entities arise from different
    positions in this hierarchy:

    GBM (IDH-wildtype): Arises from NSC or astrocyte
      progenitor. The false attractor is a dedifferentiated
      proliferative neural progenitor-like state.

    Astrocytoma (IDH-mutant): Arises from astrocytic
      lineage neural progenitor. The false attractor is
      an astrocyte-like proliferative state, slower moving
      than GBM.

    Oligodendroglioma (IDH-mutant + 1p/19q): Arises from
      oligodendrocyte precursor cell (OPC) or a progenitor
      with oligodendroglial potential. The false attractor
      retains OPC markers (PDGFRA, OLIG2) and is the
      most differentiated-looking of the three.

  THE KEY STRUCTURAL FACT:
    All three entities arise from the same broad
    neural stem/progenitor hierarchy — but they
    arrest at different positions and in different
    lineage branches. This is directly analogous to
    the RCC situation: ccRCC, PRCC, and chRCC all
    arise from the renal tubular hierarchy but diverge
    at different saddle points.

    The critical difference from RCC:
    In RCC, the four subtypes all arise from the
    kidney. In glioma, the three entities arise from
    genuinely distinct cell populations within the
    brain (NSC vs. astrocytic progenitor vs. OPC).
    The cross-glioma analysis will test whether there
    are universal markers despite these different origins.
```

---

## SECTION III — ENTITY 1: GLIOBLASTOMA, IDH-WILDTYPE (TRUE GBM)

```
CLINICAL FACTS:
  Prevalence:       ~90% of all glioblastomas
                    ~3.2 cases per 100,000 per year (USA)
                    ~14,000 new cases/year in the USA
  Age at diagnosis: Median ~64 years. Rare under 40.
                    If a young patient has apparent GBM
                    histology, IDH testing is essential —
                    it may be astrocytoma grade 4, not GBM.
  Sex:              Male > Female (~1.6:1)
  5-year survival:  <7% — one of the lowest of any cancer
  Median OS:        15–21 months with full treatment
                    (Stupp protocol + TTFields + maintenance
                    temozolomide — best current standard)
                    ~3 months untreated or best supportive
                    care only
  Standard of care:
    Newly diagnosed:
      1. Maximal safe surgical resection
      2. Concurrent radiation + temozolomide (Stupp, 2005)
         60 Gy over 6 weeks + TMZ daily
      3. Adjuvant temozolomide x 6 cycles
      4. Tumor Treating Fields (TTFields / Optune):
         Wearable device using alternating electric fields
         to disrupt mitosis. Added to maintenance TMZ.
         EF-14 trial (2017): improved median OS from
         16.0 → 20.9 months. Now standard in eligible
         patients.
      5. MGMT promoter methylation status guides TMZ
         benefit: methylated MGMT → TMZ most effective
    Recurrence:
      No universal standard.
      Options: bevacizumab (anti-VEGF, FDA approved),
               re-irradiation, clinical trials,
               re-resection.
      No approved therapy has demonstrated OS benefit
      at recurrence beyond modest PFS improvement
      with bevacizumab.
  Key unmet need:    Recurrent GBM has NO effective
                    treatment. Median OS after recurrence
                    is 6–9 months. This is one of the
                    most refractory cancers in oncology.
                    Immunotherapy has failed in every
                    major trial to date (GBM is an
                    immunosuppressive microenvironment).

CELL OF ORIGIN:
  Neural stem cell (NSC) or astrocytic progenitor cell —
  the neural progenitor population of the subventricular
  zone (SVZ) or hippocampal dentate gyrus.

  Evidence:
  - GBM cells express NSC markers: NES (nestin), SOX2,
    CD133 (PROM1), and OLIG2 in stem-like subpopulations
  - The SVZ is spatially associated with GBM origin in
    many cases — GBMs adjacent to the SVZ behave more
    aggressively and may seed ventricular spread
  - Mouse models: targeted mutations in NSCs recapitulate
    GBM more faithfully than mutations in mature astrocytes
  - Single-cell RNA-seq of GBM reveals tumor cells that
    map to NSC, OPC, astrocyte, and mesenchymal states —
    all within the SAME tumor (intratumoral heterogeneity)
    suggesting origin from a NSC that can adopt multiple
    glial identities

  Normal NSC identity programme:
    NES (nestin):   Intermediate filament — NSC marker
    SOX2:           Core pluripotency TF — NSC identity
    PROM1 (CD133):  NSC surface marker
    GFAP:           Low in NSC, high in astrocytes
    OLIG2:          Expressed in NSC and OPCs
    VIM:            Vimentin — NSC and progenitor marker
    The normal NSC is:
    - Slowly cycling (quiescent in adult brain)
    - Dependent on niche signals (VEGF, EGF, FGF)
    - Capable of producing all glial lineages
    In GBM: the NSC cannot exit the proliferative state.
    The differentiation signal is blocked.
    The result is a perpetually cycling progenitor that
    accumulates genomic instability and invades the brain.

INITIATING MOLECULAR EVENTS:
  These are well-established and form the basis of the
  WHO 2021 molecular diagnostic criteria:

  TERT promoter mutation: ~80% of IDH-wildtype GBM
                    C228T or C250T mutations.
                    Reactivates telomerase — allows
                    the cancer cell to divide indefinitely
                    without telomere shortening.
                    The earliest detectable molecular event
                    in GBM evolution.

  Chromosome 7 gain + Chromosome 10 loss (+7/−10):
                    ~80% of IDH-wildtype GBM
                    Chromosome 7 contains EGFR and MET —
                    gaining an extra copy amplifies RTK
                    signalling. Chromosome 10 contains
                    PTEN — losing a copy disables the
                    PI3K brake.
                    This single copy number event does
                    more to activate growth signalling
                    than any single point mutation.

  EGFR amplification: ~45% of IDH-wildtype GBM
                    Further amplification beyond the
                    chromosome 7 gain.
                    EGFRvIII mutation (in-frame deletion
                    of exons 2-7): ~25% — constitutively
                    active EGFR regardless of ligand.
                    This is the most studied drug target
                    in GBM, and also the most frustrating —
                    every EGFR-targeted clinical trial in
                    GBM has failed to date.

  CDKN2A/B deletion: ~60% of IDH-wildtype GBM
                    Homozygous deletion of p16INK4a
                    and p14ARF — removes both the
                    cell cycle brake (CDK4/6 inhibition)
                    and the p53 activating arm.
                    Combined with PTEN loss and EGFR
                    amplification: the three key brakes
                    on proliferation are simultaneously
                    released.

  PTEN mutation/loss: ~40% of IDH-wildtype GBM
                    PI3K/AKT/mTOR constitutive activation.

  TP53 mutation:    ~30% of IDH-wildtype GBM
                    (lower than in astrocytoma, where
                    TP53 mutation is ~80%)
                    When present, associated with the
                    proneural TCGA subtype.

  NF1 loss:         ~15% of IDH-wildtype GBM
                    RAS pathway activator.
                    Strongly associated with the
                    mesenchymal TCGA subtype.

  IDH status:       WILDTYPE — definitional.
                    IDH1/2 wildtype means the normal
                    αKG/TCA cycle is intact.
                    The TCA→αKG→EZH2 circuit studied in
                    RCC is NOT the primary mechanism here.
                    (This is a critical difference from
                    RCC and will be assessed in the
                    before-document.)

  MGMT promoter methylation:
                    ~45% of IDH-wildtype GBM
                    MGMT is a DNA repair gene that
                    reverses the damage TMZ causes.
                    Methylation = MGMT silenced =
                    TMZ works better.
                    This is currently the ONLY validated
                    predictive biomarker in GBM.

THE THREE TCGA MOLECULAR SUBTYPES WITHIN IDH-WILDTYPE GBM:

  TCGA (Verhaak et al., 2010) identified three subtypes
  using gene expression profiling. These subtypes are
  within IDH-wildtype GBM only (after IDH-mutant cases
  are removed). They are not separate diseases — they are
  expression states within the same false attractor.

  SUBTYPE A — CLASSICAL:
    Defining alteration: EGFR amplification
    Key features:
      EGFR amplified and/or mutated (EGFRvIII)
      Chromosome 7 gain / 10 loss
      NES (nestin) high
      NOTCH and Sonic Hedgehog pathway activation
      Low TP53 mutation rate
      High CDKN2A deletion rate
    Cell of origin mapping: Closest to NSC
    Expression programme: Neural stem cell-like
    % of IDH-wildtype GBM: ~36%
    Clinical note: EGFR-targeted therapies have failed
      despite this subtype being defined by EGFR.
      Resistance mechanisms include EGFRvIII
      heterogeneity, RTK switching, and BBB
      penetration issues.

  SUBTYPE B — PRONEURAL:
    Defining alteration: PDGFRA amplification / TP53
    Key features:
      PDGFRA amplified or mutated
      TP53 mutation (when present in IDH-wildtype)
      OLIG2 high
      SOX2 high
      IDH1 R132H mutation historically enriched here —
      BUT this subtype in the pre-2021 TCGA included
      IDH-mutant cases. After IDH-mutant removal,
      true IDH-wildtype proneural GBM is smaller.
    Cell of origin mapping: Closest to OPC
    Expression programme: Oligodendrocyte precursor-like
    % of IDH-wildtype GBM: ~30% (reduced after IDH
      purification)
    Clinical note: Historically appeared to have better
      prognosis — but this was entirely explained by
      the inclusion of IDH-mutant cases (which are
      less aggressive). True IDH-wildtype proneural
      GBM does not have better prognosis.

  SUBTYPE C — MESENCHYMAL:
    Defining alteration: NF1 loss
    Key features:
      NF1 loss (mutation or deletion)
      CHI3L1 (YKL-40) high
      CD44 high
      LGALS3 high
      MET expression elevated
      High macrophage / microglia infiltration
      Highest immune gene expression of the three subtypes
      High expression of mesenchymal markers (VIM, FN1)
    Cell of origin mapping: Astrocytic progenitor or
      NSC undergoing mesenchymal transition
    Expression programme: Reactive astrocyte / mesenchymal
    % of IDH-wildtype GBM: ~34%
    Clinical note: Most immune-infiltrated subtype.
      Despite this, immune checkpoint inhibitors have
      not improved outcomes in mesenchymal GBM in
      trials to date. The immune infiltrate is
      predominantly immunosuppressive macrophages
      (M2 polarized), not cytotoxic T cells.

  IMPORTANT CAVEAT ON THE THREE SUBTYPES:
    These three subtypes are not fixed cellular states —
    they are expression snapshots. GBM shows significant
    intratumoral heterogeneity at the single-cell level:
    cells within a single tumor can map to proneural,
    classical, and mesenchymal expression states
    simultaneously. Subtype assignment from bulk RNA-seq
    reflects the DOMINANT expression state.
    Recurrent GBM frequently shifts from proneural or
    classical to mesenchymal — a treatment-driven
    transition analogous to the neuroendocrine
    transformation in prostate cancer.
    The mesenchymal shift at recurrence is a Waddington
    transition within the false attractor landscape:
    one false attractor position → another false attractor
    position, driven by treatment pressure.

AVAILABLE PUBLIC DATA:
  TCGA-GBM:         n~590 total
                    ~480 IDH-wildtype by modern calls
                    RNA-seq, methylation, CNV
                    GDC portal (phs000178 / TCGA-GBM)
                    NOTE: TCGA-GBM has very few normal
                    brain tissue samples (~5 adjacent
                    normal). Normal reference must come
                    from external sources.
  CGGA (Chinese Glioma Genome Atlas):
    CGGA_693:       n=693 total gliomas
                    ~160–170 IDH-wildtype GBM
                    RNA-seq, survival annotated
                    Available at cgga.org.cn
    CGGA_325:       n=325 total, RNA-seq
                    ~80–90 IDH-wildtype GBM
  GEO datasets:
    GSE4290:        n=180 (gliomas of various grades +
                    n=23 non-tumour brain)
                    Affymetrix. CRITICAL normal reference.
                    One of the few datasets with true
                    non-tumour brain tissue.
    GSE16011:       n=276 gliomas (grades II–IV)
                    Includes WHO grade IV (GBM)
                    Non-tumour brain included
    GSE83300:       Paired primary/recurrent GBM
                    Critical for recurrence analysis
    GSE108474:      IDH-stratified glioma dataset
    GSE162631:      Single-cell RNA-seq GBM
  Normal brain reference:
    GTEx brain:     Multiple brain region RNA-seq
                    Available via GTEx portal
                    Most comprehensive normal brain
                    transcriptome available.
                    Critical — will be the primary
                    normal reference for GBM analysis.

ANALYSIS STATUS:    NOT YET STARTED
PLANNED DOCUMENT:   GBM-S1a (IDH-wildtype GBM before-doc)
PLANNED DATASET:    TCGA-GBM (IDH-wildtype subset)
                    + GSE4290 (normal brain reference)
                    + GTEx brain (normal reference)
PRIORITY:           FIRST — IDH-wildtype GBM is what
                    the existing analysis predominantly
                    captured. The highest clinical
                    unmet need. The most tractable
                    Waddington analysis within the CNS.
```

---

## SECTION IV — ENTITY 2: ASTROCYTOMA, IDH-MUTANT (GRADES 2–4)

```
CLINICAL FACTS:
  Prevalence:       ~30% of all diffuse gliomas
                    Predominantly in young adults (20–45)
  Grades:           WHO grades 2, 3, and 4
                    Grade 4 astrocytoma, IDH-mutant =
                    the former "secondary GBM"
  5-year survival:  Grade 2: ~70–80%
                    Grade 3: ~50–60%
                    Grade 4: ~30–40%
                    (ALL better than IDH-wildtype GBM)
                    The survival difference vs. GBM is
                    almost entirely explained by IDH mutation.
  Standard of care: Surgery (maximal safe resection)
                    Radiation + PCV chemotherapy
                    (procarbazine, CCNU/lomustine,
                    vincristine) for grade 3/4
                    Temozolomide (increasingly used for
                    grade 3/4 — CATNON trial data)
                    Ivosidenib (IDH1 inhibitor):
                    FDA approved 2023 for IDH1-mutant
                    low-grade glioma — first targeted
                    therapy for this entity
                    Grade 4: TMZ + radiation (Stupp-like)
                    but outcomes still much better than
                    IDH-wildtype GBM

CELL OF ORIGIN:
  Astrocytic lineage neural progenitor — distinct from
  the NSC of IDH-wildtype GBM.
  Evidence:
    - ATRX mutations (see below) affect the chromatin
      landscape of astrocytic progenitors specifically
    - TP53 mutations present in astrocytic cells
    - IDH mutation occurs in a partially committed
      astrocytic progenitor, not in the most primitive NSC
    - Mouse models: IDH mutation in GFAP+ astrocytic
      progenitors produces astrocytoma, not GBM

  Normal astrocytic progenitor identity:
    GFAP:           Astrocytic intermediate filament
                    (present but low in progenitors)
    VIM:            Vimentin — reactive astrocyte marker
    ALDH1L1:        Pan-astrocyte marker
    SOX9:           Astrocytic lineage TF
    NFIA:           Nuclear factor I A — astrocyte TF
    NES:            Nestin — expressed in progenitor state

DEFINING MOLECULAR EVENT — IDH MUTATION:
  IDH1 R132H:       ~90% of IDH-mutant astrocytomas
  IDH2 R172K:       ~5% (less common)
  The IDH mutation is the founding event.
  It converts αKG (alpha-ketoglutarate) into
  2-hydroxyglutarate (2-HG), a "oncometabolite."
  2-HG inhibits αKG-dependent dioxygenases including:
    TET2:           DNA demethylase — inhibited by 2-HG
                    → CpG island methylator phenotype (CIMP)
                    → widespread DNA hypermethylation
                    → silencing of differentiation genes
    KDM histone demethylases: inhibited by 2-HG
                    → histone hypermethylation
                    → epigenetic lock
  The IDH mutation therefore creates an EPIGENETIC
  LOCK through the TCA cycle — exactly the mechanism
  identified in RCC (TCA→αKG→epigenetic writers).
  In RCC: TCA gene suppression → reduced αKG → EZH2
  In IDH-mutant glioma: IDH mutation → 2-HG → TET2
  inhibition → DNA hypermethylation lock.
  The direction is different (loss of αKG in RCC via
  TCA suppression; gain of 2-HG in IDH-mutant glioma
  via neomorphic enzyme activity) but the mechanism
  is structurally analogous.
  This cross-cancer structural observation will be
  addressed in the before-document.

ADDITIONAL MOLECULAR EVENTS:
  TP53 mutation:    ~80% of IDH-mutant astrocytoma
                    (much higher than IDH-wildtype GBM)
  ATRX loss:        ~80% of IDH-mutant astrocytoma
                    ATRX is a chromatin remodeling factor
                    required for histone H3.3 deposition
                    at telomeres. Loss of ATRX activates
                    ALT (alternative lengthening of
                    telomeres) — the cancer cell no longer
                    needs TERT to maintain telomere length.
                    This is the opposite of IDH-wildtype GBM
                    where TERT promoter mutation reactivates
                    TERT. Two different solutions to the
                    same problem (telomere maintenance).
  CDKN2A deletion:  Grade 4 astrocytoma (IDH-mutant, G4)
                    requires CDKN2A/B homozygous deletion
                    by WHO 2021 for the grade 4 designation
                    in the absence of necrosis/MVP.
  No 1p/19q codeletion (distinguishes from oligodendroglioma)
  No TERT promoter mutation (distinguishes from GBM and
  oligodendroglioma)

KEY DISTINGUISHING FEATURES vs. GBM (IDH-wildtype):
  - IDH1/2 MUTATED (wildtype in GBM)
  - ATRX loss (intact in GBM)
  - TP53 mutation ~80% (only ~30% in GBM)
  - No EGFR amplification
  - No +7/−10 copy number signature
  - CIMP (CpG island methylator phenotype) — widespread
    DNA hypermethylation (absent in GBM)
  - Better prognosis (years, not months)
  - Younger patients (20s–40s, not 60s–70s)
  - IDH1 inhibitor (ivosidenib) active — first targeted
    therapy for this entity
  - Slow progression with late recurrence

AVAILABLE PUBLIC DATA:
  TCGA-LGG:         n=516 lower grade gliomas
                    Contains grades 2–3, mostly IDH-mutant
                    Including astrocytoma and oligodendroglioma
                    IDH status annotated in clinical file
  TCGA-GBM:         ~20% IDH-mutant (the old "secondary GBM")
                    These need to be SEPARATED from
                    IDH-wildtype before GBM analysis
  CGGA:             Contains IDH-mutant astrocytoma samples
                    IDH status annotated
  GEO:
    GSE16011:       Includes grade II/III astrocytoma
    GSE108474:      IDH-stratified, excellent for
                    astrocytoma vs. GBM comparison

ANALYSIS STATUS:    NOT YET STARTED
PLANNED DOCUMENT:   GBM-S2a (IDH-mutant astrocytoma
                    before-document)
PLANNED DATASET:    TCGA-LGG (IDH-mutant astrocytoma subset)
                    + CGGA (astrocytoma grades 2–4)
PRIORITY:           SECOND — the IDH mutation creates a
                    direct TCA→epigenetic lock mechanism
                    that connects this entity structurally
                    to the RCC work. The ivosidenib
                    approval (2023) makes this clinically
                    timely. The αKG axis is the unifying
                    thread.
```

---

## SECTION V — ENTITY 3: OLIGODENDROGLIOMA, IDH-MUTANT AND 1P/19Q-CODELETED

```
CLINICAL FACTS:
  Prevalence:       ~15% of all diffuse gliomas
  Grades:           WHO grades 2 and 3
                    No grade 4 oligodendroglioma (by
                    WHO 2021 — the old grade 4 oligo
                    is now reclassified)
  Age at diagnosis: 35–55 years (younger than GBM,
                    older than astrocytoma)
  5-year survival:  Grade 2: ~90%+
                    Grade 3: ~60–70%
                    The most favorable prognosis of all
                    diffuse gliomas — median OS can
                    exceed 15 years in grade 2
  Standard of care: Surgery
                    Radiation + PCV chemotherapy
                    (RTOG 9402 / EORTC 26951 trials:
                    PCV after radiation improved OS by
                    ~7 years in grade 2/3 oligo)
                    IDH inhibitors being studied
                    (vorasidenib — IDH1/2 inhibitor —
                    INDIGO trial 2023 showed benefit in
                    grade 2 IDH-mutant glioma including
                    oligodendroglioma)

CELL OF ORIGIN:
  Oligodendrocyte precursor cell (OPC) — the PDGFRA+,
  OLIG2+, SOX10+ progenitor of the oligodendroglial
  lineage.
  Evidence:
    - Oligodendroglioma cells retain OPC markers
      (PDGFRA, OLIG2, SOX10) — the cell of origin
      programme is partially preserved
    - 1p/19q codeletion is a structural variant that
      appears to arise in cells of OPC origin
    - Mouse models: OPC-targeted IDH1 mutation + 1p/19q
      codeletion events recapitulate oligodendroglioma

  Normal OPC identity:
    PDGFRA:         OPC marker and growth factor receptor
    OLIG2:          Master TF for oligodendroglial lineage
    OLIG1:          Co-TF with OLIG2
    SOX10:          Oligodendroglial lineage TF
    NKX2-2:         OPC identity TF
    NG2 (CSPG4):    OPC surface marker
    MBP:            Low in OPC, high in mature oligo
    The normal OPC is:
    - A slowly cycling progenitor in adult brain
    - Activated by demyelination injury to produce
      new oligodendrocytes
    - PDGFRA-dependent for survival and proliferation

DEFINING MOLECULAR EVENTS:
  IDH1/2 mutation:  Required (same 2-HG oncometabolite
                    mechanism as astrocytoma)
  1p/19q codeletion: Required (defines the entity)
                    Concurrent loss of chromosome arms
                    1p and 19q via an unbalanced
                    translocation. The genes on these
                    arms that drive transformation are
                    not fully defined — multiple candidate
                    tumor suppressors lost simultaneously.
                    CIC (1p) and FUBP1 (19q) are the
                    most recurrently mutated at the
                    remaining allele.
  TERT promoter mutation: ~90% of oligodendrogliomas
                    (high — same as IDH-wildtype GBM,
                    but different context: in oligo it
                    co-occurs with IDH mutation and 1p/19q)
  No ATRX loss     (distinguishes from astrocytoma —
                    ATRX and 1p/19q codeletion are
                    mutually exclusive)
  No TP53 mutation  (usually — distinguishes from
                    astrocytoma)
  PIK3CA:           ~10%
  NOTCH1:           ~10%

KEY DISTINGUISHING FEATURES vs. ASTROCYTOMA:
  - 1p/19q CODELETED (absent in astrocytoma)
  - TERT mutation HIGH (rare in astrocytoma)
  - ATRX INTACT (lost in astrocytoma)
  - TP53 WILDTYPE (mutated in astrocytoma)
  - OPC expression programme (PDGFRA, OLIG2, SOX10)
  - Best prognosis of all diffuse gliomas
  - PCV chemotherapy particularly effective

AVAILABLE PUBLIC DATA:
  TCGA-LGG:         Contains oligodendroglioma
                    1p/19q and IDH status annotated
                    Can be separated from astrocytoma
                    using clinical file
  CGGA:             Contains oligodendroglioma subset
  GEO:
    GSE16011:       Contains oligodendroglioma grade 2/3
    GSE108474:      IDH-stratified — includes oligo

ANALYSIS STATUS:    NOT YET STARTED
PLANNED DOCUMENT:   GBM-S3a (Oligodendroglioma before-doc)
PLANNED DATASET:    TCGA-LGG (oligo subset)
                    + CGGA (oligo subset)
PRIORITY:           THIRD — after GBM and astrocytoma.
                    Smallest of the three entities.
                    Most differentiated false attractor.
                    The OPC identity programme retained
                    in the false attractor makes this
                    the most "normal-looking" glioma —
                    and therefore the most interesting
                    test of the depth score framework.
```

---

## SECTION VI — THE MGMT METHYLATION AXIS — A SPECIAL STRUCTURAL NOTE

```
MGMT promoter methylation is the only validated predictive
biomarker in GBM. It deserves specific attention before
any analysis begins.

WHAT MGMT IS:
  O6-methylguanine-DNA methyltransferase.
  A DNA repair enzyme that removes the methyl group
  that TMZ adds to guanine. When MGMT is active, the
  cell repairs TMZ-induced damage — TMZ fails.
  When MGMT promoter is methylated (silenced), the
  cell cannot repair TMZ damage — TMZ works.

PREVALENCE:
  ~45% of IDH-wildtype GBM have MGMT promoter methylation
  ~75% of IDH-mutant astrocytoma/glioma have MGMT methylation
  (IDH mutation drives CpG island methylation broadly —
  CIMP — which includes MGMT methylation)

CLINICAL SIGNIFICANCE:
  Methylated MGMT GBM:    Median OS ~22–24 months
  Unmethylated MGMT GBM:  Median OS ~13–15 months
  This is the largest prognostic split in GBM.

WADDINGTON RELEVANCE:
  MGMT methylation in GBM is a depth-relevant marker:
  It measures the epigenetic silencing at one specific
  gene. The deeper into the false attractor, the more
  broadly the epigenetic programme is altered.
  In IDH-mutant glioma: MGMT methylation is part of
  the CIMP programme driven by 2-HG.
  In IDH-wildtype GBM: MGMT methylation is independent
  of IDH mutation — it is a separate epigenetic event.
  The depth score analysis should correlate with MGMT
  methylation status as a validation step — higher
  depth should predict methylated MGMT in IDH-wildtype
  GBM if the epigenetic lock scales with depth.
  This is stated here as a structural hypothesis to
  be tested, NOT as a prediction for the before-document.
```

---

## SECTION VII — THE EXISTING GBM ANALYSIS — CONTEXT

```
The existing analysis in Cancer_Research/GBM/ ran on
TCGA-GBM without IDH stratification.

WHAT THE EXISTING ANALYSIS CAPTURED:
  The TCGA-GBM dataset is ~80% IDH-wildtype and ~20%
  IDH-mutant. The depth axis found in the existing
  analysis is therefore predominantly driven by the
  IDH-wildtype landscape — which is the correct target.

  However:
  - The ~20% IDH-mutant cases in the dataset will have
    introduced variance unrelated to the IDH-wildtype
    false attractor. The IDH-mutant tumors have a
    different depth axis (CIMP-driven, different genes).
  - The depth axis from the combined dataset is a
    weighted average of both landscapes.
  - Genes that appear depth-positive in both landscapes
    (EZH2, MKI67, etc.) will be correctly identified.
  - Genes that are depth-positive in IDH-wildtype GBM
    but depth-negative in IDH-mutant glioma (or vice
    versa) will be diluted or hidden.

WHAT THE SUBTYPE ANALYSES WILL ADD:
  GBM-S1 (IDH-wildtype): Clean depth axis for true GBM.
    Will confirm or refine what the existing analysis found.
    Will resolve which genes are IDH-wildtype specific.

  GBM-S2 (IDH-mutant astrocytoma): A completely different
    depth axis centered on the 2-HG/CIMP mechanism.
    The TCA→αKG→EZH2 circuit from RCC will be tested
    here in a different context (IDH mutation as the
    TCA disruptor rather than gene suppression).

  GBM-S3 (Oligodendroglioma): The most differentiated
    glioma. The depth score here measures how far an
    OPC has drifted from its normal identity.

THE EXISTING ANALYSIS IS NOT DISCARDED.
It represents the combined landscape — valid as a
first-pass view of what separates brain tumours from
normal brain. The subtype analyses go deeper.
```

---

## SECTION VIII — THE RECURRENCE TRANSITION — A STRUCTURAL NOTE

```
GBM recurrence is itself a Waddington event.

At diagnosis, GBM bulk RNA-seq shows a dominant
expression subtype (classical, proneural, or mesenchymal).

At recurrence (post-TMZ and radiation):
  ~30% of proneural GBMs shift to mesenchymal subtype.
  Classical GBMs also trend toward mesenchymal on recurrence.
  Mesenchymal GBMs at recurrence are DIFFERENT from
  mesenchymal GBMs at diagnosis — the recurrence-emergent
  mesenchymal state appears driven by:
    - NF1 loss (acquired at recurrence)
    - YAP/TAZ activation
    - Treatment-induced macrophage recruitment
      remodeling the microenvironment
    - TP53 mutation accumulation

This recurrence transition is structurally analogous to
the neuroendocrine transformation in prostate cancer:
  - Treatment pressure drives the cell deeper into the
    Waddington landscape
  - The cell escapes the drug by acquiring a new
    false attractor identity
  - The new identity (mesenchymal) is less drug-sensitive

The paired primary/recurrent datasets (GSE83300) will
be used in the GBM-S1 analysis to test whether the
depth score at diagnosis predicts mesenchymal shift
at recurrence. This is stated here as a structural
opportunity, not as a prediction.
```

---

## SECTION IX — DATA AVAILABILITY SUMMARY

```
Entity         Primary Dataset      n (tumour)  Normal Ref       Power
───────────────────────────────────────────────────────────────────────
IDH-wt GBM    TCGA-GBM IDH-wt     ~480        GTEx + GSE4290   HIGH
              CGGA_693 IDH-wt     ~165        GTEx             MOD
IDH-mut Astro TCGA-LGG astro      ~200        GTEx + GSE4290   HIGH
              CGGA astro subset   ~150        GTEx             MOD
Oligo         TCGA-LGG oligo      ~170        GTEx + GSE4290   MOD
              CGGA oligo subset   ~100        GTEx             MOD
Normal brain  GTEx brain          n=2000+     —                REFERENCE
              GSE4290 normal      n=23        —                REFERENCE

CRITICAL NOTE ON NORMAL REFERENCE:
  GBM is unique in this repository: the "normal" tissue
  is living brain — not surgically removed normal adjacent
  tissue (which barely exists for brain tumours).
  The depth score axis will be defined relative to
  GTEx brain (multiple regions) as the normal reference.
  GSE4290 provides a smaller but directly GEO-accessible
  normal brain reference that can be used for validation.
  This decision — which normal reference to use — must be
  locked in the before-document before any data loads.
```

---

## SECTION X — PLANNED ANALYSIS ORDER

```
ORDER:

  GBM-S1   IDH-wildtype GBM    TCGA-GBM IDH-wt     HIGH POWER
                                + GTEx normal
                                REASON: The primary entity.
                                What the existing analysis
                                predominantly captured.
                                The highest clinical unmet
                                need in CNS oncology.
                                The three TCGA subtypes
                                (classical, proneural,
                                mesenchymal) will be tested
                                as depth positions within
                                the single IDH-wt landscape.

  GBM-S2   IDH-mutant          TCGA-LGG astro      HIGH POWER
           Astrocytoma         + CGGA
                                REASON: The 2-HG/CIMP
                                mechanism is the direct
                                structural parallel to
                                the TCA→αKG axis in RCC.
                                This analysis connects
                                the GBM work to the RCC
                                work at the mechanistic
                                level. Ivosidenib approval
                                makes this clinically timely.

  GBM-S3   Oligodendroglioma   TCGA-LGG oligo      MODERATE
                                + CGGA
                                REASON: Most differentiated
                                glioma. Tests whether the
                                depth score framework works
                                for a cancer where the
                                false attractor RETAINS
                                significant normal identity.

  GBM-X    Cross-glioma        After S1–S3 complete
           comparison          Questions:
                                  1. Universal glioma
                                     depth markers?
                                  2. EZH2 in all three?
                                     Same or different
                                     direction?
                                  3. TCA/αKG axis: present
                                     in IDH-wildtype GBM
                                     through a different
                                     mechanism?
                                  4. Basket trial candidates
                                     across all three glioma
                                     entities?
                                  5. Does the mesenchymal
                                     shift at recurrence
                                     map to a deeper
                                     depth position?
```

---

## SECTION XI — WHAT THIS DOCUMENT DOES NOT CONTAIN

```
This document contains:
  ✓ The WHO 2021 reclassification and its consequences
  ✓ The three glioma entities (GBM, astrocytoma, oligo)
  ✓ Cells of origin for each entity
  ✓ Published molecular drivers (literature, not prediction)
  ✓ Clinical characteristics (survival, standard of care)
  ✓ The three TCGA subtypes within IDH-wildtype GBM
  ✓ The MGMT methylation axis
  ✓ The recurrence transition as a Waddington event
  ✓ The IDH mutation / 2-HG / αKG structural connection
    to the RCC work (noted as observation, not prediction)
  ✓ Context for the existing combined GBM analysis
  ✓ Data availability (TCGA, CGGA, GEO, GTEx)

This document does NOT contain:
  ✗ Depth score predictions
  ✗ Specific false attractor gene predictions
  ✗ Drug target predictions
  ✗ Epigenetic lock mechanism hypotheses
  ✗ Cross-glioma structural predictions beyond literature

All of the above belong in the BEFORE documents.
GBM-S1a (IDH-wildtype GBM before-document) is next.
Written before any script runs. Before any data loads.
```

---

## STATUS BLOCK

```
document:           GBM_Subtype_Orientation.md
folder:             Cancer_Research/GBM/Subtypes/
status:             COMPLETE — ORIENTATION ONLY
date:               2026-03-03
author:             Eric Robert Lawson
                    OrganismCore

entities_covered:
  primary:          IDH-wildtype GBM              [1 of 3]
  secondary:        IDH-mutant Astrocytoma        [2 of 3]
  tertiary:         Oligodendroglioma             [3 of 3]

analyses_started:   0

existing_analysis:  Cancer_Research/GBM/ — conducted on
                    full TCGA-GBM (~80% IDH-wt, ~20%
                    IDH-mutant, no IDH stratification).
                    Valid as a combined first-pass analysis.
                    The IDH-wildtype signal dominates.
                    Decomposed by the analyses in this folder.

next_document:      GBM-S1a
                    IDH-Wildtype GBM Before-Document
                    (predictions locked before
                    TCGA-GBM IDH-wt subset loads)

critical_note_1:    "GBM" post-2021 means IDH-WILDTYPE ONLY.
                    IDH-mutant grade 4 glioma is a separate
                    disease (Astrocytoma, IDH-mutant, G4).
                    These are not subtypes of GBM.
                    They are three separate Waddington
                    landscapes in the same organ.

critical_note_2:    The IDH mutation creates a direct
                    structural connection to the RCC work:
                    IDH mutation → 2-HG → αKG-dependent
                    enzyme inhibition → epigenetic lock.
                    This is the TCA→αKG→EZH2 circuit
                    identified in RCC, operating through
                    a different entry point.
                    This connection will be tested formally
                    in GBM-S2a (astrocytoma before-document).
```
