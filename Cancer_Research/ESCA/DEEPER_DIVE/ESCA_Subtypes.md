# ESOPHAGEAL CANCER — SUBTYPE ORIENTATION DOCUMENT
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

  A complete map of the esophageal cancer subtype landscape —
  what the two primary subtypes are, where they come from
  anatomically, why they are biologically distinct diseases
  that share only an organ location, what the clinical and
  molecular literature has established about each one, and
  what public data exists to analyze them.

Esophageal cancer is arguably the clearest example in this
entire repository of a single clinical label applied to
two completely unrelated diseases.

ESCC (squamous cell carcinoma) and EAC (adenocarcinoma)
arise from different cells in different parts of the
esophagus, are driven by different molecular events,
have different geographic distributions, different risk
factors, different chemotherapy response profiles, and
different prognoses.

The only thing they share is the organ.

The existing ESCA analysis in Cancer_Research/ESCA/ was
conducted on TCGA-ESCA, which contains both subtypes.
The first task of this Subtypes/ folder is to establish
which subtype the existing analysis predominantly captured,
and then to conduct dedicated subtype analyses for each.

This document is the prerequisite to every analysis in
Cancer_Research/ESCA/Subtypes/.
```

---

## DOCUMENT METADATA

```
document_id:        ESCA_Subtype_Orientation
series:             ESCA (Esophageal Cancer — Subtypes)
folder:             Cancer_Research/ESCA/Subtypes/
date:               2026-03-03
author:             Eric Robert Lawson
status:             ORIENTATION ONLY — no predictions,
                    no analysis, no epigenetic hypotheses
next_document:      ESCA_ESCC_before.md
                    (Document ESCA-S1a — ESCC before-doc)
protocol_version:   Workflow_Protocol.md v2.0
```

---

## SECTION I — THE FUNDAMENTAL PROBLEM WITH "ESOPHAGEAL CANCER"

```
"Esophageal cancer" is the single most misleading cancer
label in the repository in terms of biological coherence.

The esophagus is a muscular tube approximately 25 cm long
running from the pharynx to the stomach. It is lined by
two completely different epithelial cell types depending
on location:

  UPPER and MIDDLE esophagus:
    Lined by stratified SQUAMOUS epithelium —
    the same cell type that lines the skin, mouth,
    throat, and cervix.
    Function: mechanical protection against food passage.
    No secretory function.

  LOWER esophagus / gastroesophageal junction:
    Lined by SQUAMOUS epithelium in normal conditions —
    but in the setting of chronic acid reflux
    (gastroesophageal reflux disease, GERD), this
    squamous epithelium undergoes METAPLASIA:
    it is replaced by columnar/glandular epithelium
    that resembles intestinal mucosa.
    This metaplastic state is called BARRETT'S ESOPHAGUS.
    Barrett's esophagus is the precursor to
    esophageal adenocarcinoma.

TWO SUBTYPES, TWO DISEASES:

  ESCC — Esophageal Squamous Cell Carcinoma
    Arises from squamous epithelial cells
    Upper/middle esophagus
    Risk factors: tobacco, alcohol, hot beverages,
                  nutritional deficiency, HPV (in some
                  geographic regions)
    Geographic distribution: China, Iran, parts of Africa
                  ("esophageal cancer belt")
    Molecularly similar to: head and neck SCC,
                  lung squamous cell carcinoma,
                  cervical squamous cell carcinoma

  EAC — Esophageal Adenocarcinoma
    Arises from metaplastic columnar cells (Barrett's)
    Lower esophagus / gastroesophageal junction
    Risk factors: GERD, obesity, Western diet, male sex
    Geographic distribution: Western countries,
                  rising incidence in UK, USA, Australia
    Molecularly similar to: gastric adenocarcinoma
                  (chromosomally unstable subtype, CIN)

These are not subtypes of the same disease.
They are two diseases that happen to present in the
same anatomical structure.
The TCGA ESCA dataset contains both —
and the existing analysis must be interpreted with
that in mind.
```

---

## SECTION II — THE ESOPHAGEAL ANATOMY AND NORMAL CELL HIERARCHY

```
Understanding the normal esophageal anatomy is essential
for defining the Waddington landscapes in each subtype.

THE SQUAMOUS ESOPHAGUS (upper and middle):

  Architecture:
    Basal layer:      Stem/progenitor cells
                      Express p63, CK5/6, CK14
                      Divide to replenish the epithelium
    Suprabasal layers: Post-mitotic, differentiating cells
                      Flatten and keratinize toward surface
    Surface layer:    Squamous, non-keratinizing cells
                      (unlike skin — the esophagus does
                      not form a true keratin layer)

  Normal identity programme of the squamous cell:
    TP63 (p63):       Master TF of squamous epithelium
                      Required for squamous identity and
                      stratification maintenance
    SOX2:             Pluripotency and squamous lineage TF
                      Co-expressed with p63 in squamous cells
                      Amplified in ESCC (paradoxical — the
                      cancer amplifies a "normal" identity
                      gene to drive proliferation)
    KRT5/6:           Basal squamous cytokeratins
    KRT14:            Squamous basal layer marker
    NOTCH1:           Squamous differentiation regulator
                      (loss of NOTCH1 permits escape from
                      differentiation — key in ESCC)
    NFE2L2 (NRF2):   Oxidative stress response
                      Mutated or amplified in ~20% of ESCC

  The ESCC false attractor is:
    A cell that retains squamous identity markers (p63, SOX2)
    but has disabled the differentiation programme (NOTCH1
    loss, TP53 loss) and activated proliferative drive.
    It is a squamous cell that cannot terminally differentiate
    and instead cycles indefinitely.

THE COLUMNAR / GLANDULAR ESOPHAGUS (Barrett's, lower):

  Normal squamous epithelium of the lower esophagus is
  replaced by columnar epithelium in Barrett's metaplasia.

  The Barrett's cell is an intermediate state —
  it is no longer squamous and not yet cancerous.
  Key identity markers:
    CDX2:             Intestinal transcription factor
                      Master regulator of intestinal
                      epithelial identity. Drives the
                      columnar metaplasia phenotype.
    MUC2:             Goblet cell mucin
    CK8/18:           Glandular/columnar cytokeratins
    TFF3:             Trefoil factor 3 — intestinal metaplasia
    VILLIN:           Intestinal brush border protein

  The Waddington transition in EAC:
    Normal squamous → Barrett's (CDX2+ columnar metaplasia)
    → Barrett's with dysplasia (TP53 mutation occurs here)
    → EAC (full false attractor)

  The Barrett's cell is the saddle point cell —
  it has already crossed from squamous identity to
  columnar identity, and EAC is the further descent
  from Barrett's into the false attractor.

  The EAC false attractor is:
    A glandular cell with intestinal-like identity (CDX2+)
    that has disabled apoptosis (TP53 mutation), activated
    RTK signaling (ERBB2, EGFR, VEGFA amplification),
    and disabled cell cycle control (CDKN2A loss).
    It is a quasi-intestinal cell that grows in the
    esophagus where intestinal cells do not belong.

TCGA COMPOSITION OF ESCA:
  TCGA-ESCA contains:
    n=89  ESCC samples (~52% of the dataset)
    n=80  EAC samples (~47% of the dataset)
    n=4   adenosquamous or mixed/other (~1%)
  Total: ~173 samples

  This is an unusually balanced split for a TCGA dataset —
  most TCGA datasets are dominated by one histotype.
  ESCA was deliberately designed to capture both.

  IMPLICATION FOR EXISTING ANALYSIS:
    The existing Cancer_Research/ESCA/ analysis ran on
    the combined TCGA-ESCA dataset.
    The depth score and gene correlations it produced
    reflect a MIXTURE of two fundamentally different
    Waddington landscapes.
    This does not invalidate the existing analysis —
    but it means the depth axis found represents the
    variance that separates BOTH subtypes from normal,
    not the within-subtype depth structure.
    The subtype analyses below will decompose this.
```

---

## SECTION III — SUBTYPE 1: ESOPHAGEAL SQUAMOUS CELL CARCINOMA (ESCC)

```
CLINICAL FACTS:
  Prevalence:       ~90% of esophageal cancer globally
                    (ESCC dominates worldwide)
                    ~30% in Western countries
                    (EAC has overtaken ESCC in the West)
  Geographic peak:  The "esophageal cancer belt":
                    Northern China, Central Asia (Iran,
                    Turkmenistan, Kazakhstan), parts of
                    sub-Saharan Africa, Northern France
                    (Normandy/Calvados region)
                    In China: ESCC is among the top 5
                    cancer killers.
  Age at diagnosis: Typically 55–70 years
  Sex:              Male >> Female (3:1)
  5-year survival:  ~15–25% overall
                    ~40–50% (Stage I/II — rare at diagnosis)
                    ~5–15% (Stage III/IV — most common
                    presentation)
  Standard of care: Stage I:   Surgery alone (esophagectomy)
                    Stage II/III: Neoadjuvant concurrent
                    chemoradiotherapy (cisplatin + 5-FU or
                    paclitaxel + carboplatin) followed by
                    surgery (CROSS trial regimen)
                    Stage IV:  Platinum + 5-FU + PD-1
                    inhibitor (nivolumab or pembrolizumab)
                    Immunotherapy:
                      CPS ≥10: pembrolizumab/nivolumab
                      now first-line standard in advanced ESCC
                      CheckMate-648, KEYNOTE-590 changed
                      the standard of care in 2021-2022
  pCR rate:         ~25–35% with neoadjuvant CRT
                    (lower than EAC)

CELL OF ORIGIN:
  Squamous epithelial progenitor cell (basal layer) of
  the upper or middle esophagus.

  The key normal identity:
    The basal squamous cell is a p63+ stem-like cell that
    continuously divides to replenish the squamous epithelium
    while its daughters differentiate and flatten toward
    the surface. Its identity is defined by:
      TP63 (p63):   Nuclear master TF — squamous identity
      SOX2:         Co-TF with p63 — squamous programme
      KRT5/6:       Structural identity
      NOTCH1:       Differentiation signal — when NOTCH1
                    activates, the basal cell exits the
                    cell cycle and differentiates upward.
                    Loss of NOTCH1 = failure to differentiate
                    = persistence in a proliferative basal state.

  In ESCC: the basal squamous cell fails to differentiate.
  The normal identity programme is partially retained
  (p63 and SOX2 remain expressed) but the differentiation
  exit is blocked (NOTCH1 lost) and genomic instability
  accumulates.

INITIATING MOLECULAR EVENTS:
  TP53 mutation:    ~85% of ESCC
                    The founding event in most ESCC.
                    Loss of TP53 allows the basal cell to
                    survive DNA damage and continue cycling.
  NOTCH1 mutation:  ~20% of ESCC
                    Loss of squamous differentiation signal.
  ZNF750 mutation:  ~15% of ESCC
                    Squamous differentiation TF — loss
                    cooperates with NOTCH1 loss.
  NFE2L2 (NRF2):   ~20% mutated or amplified
                    Oxidative stress survival pathway.
  PIK3CA:           ~15% of ESCC
  CDKN2A (p16):    ~75% of ESCC (mostly by deletion/
                    methylation — not mutation)
                    Loss of p16 = unrestricted CDK4/6
                    activity = uncontrolled cell cycling.
  SOX2 amplification: ~40% of ESCC
                    The squamous identity gene is amplified
                    as a pro-proliferative driver.
                    Same amplicon at 3q26-27 seen in lung
                    SCC and head and neck SCC.
  FGFR1 amplification: ~15% of ESCC
  EGFR amplification: ~15% of ESCC
  High copy number instability (but point-mutation driven
  more than copy-number driven compared to EAC)

TCGA MOLECULAR SUBTYPES WITHIN ESCC:
  TCGA (2017) identified three ESCC subtypes:
    ESCC1: NFE2L2/KEAP1 pathway altered — oxidative stress
    ESCC2: NOTCH1 mutations — differentiation failure
    ESCC3: Less defined — fewer copy number alterations
  These are sub-classifications within ESCC.
  For this framework: ESCC will be analyzed as a single
  entity first (n=89 in TCGA-ESCA). The three TCGA
  sub-subtypes are noted but will not be analyzed
  separately given small sample sizes per sub-group.

MOLECULAR SIMILARITY TO OTHER CANCERS:
  ESCC shares molecular programme with:
    Head and neck squamous cell carcinoma (HNSCC)
    Lung squamous cell carcinoma (LUSC)
    Cervical squamous cell carcinoma
  These share:
    p63/SOX2 co-amplification (3q26-27)
    NOTCH1 loss
    CDKN2A loss
    TP53 mutation
    PIK3CA mutation
  This is not coincidental — it reflects a shared
  cell of origin (squamous progenitor cell) and a
  shared false attractor (undifferentiated squamous
  proliferative state).
  The Waddington landscape for ESCC is structurally
  similar to HNSCC and LUSC. This is a cross-cancer
  prediction that will be tested.

KEY DISTINGUISHING FEATURES vs. EAC:
  - Squamous histology (not glandular)
  - p63 positive, CDX2 negative
  - Upper/middle esophagus location
  - TP53 + NOTCH1 as primary drivers
  - SOX2/p63 amplicon at 3q
  - Responds to cisplatin + 5-FU (moderate response)
  - Responds to PD-1 inhibitors (high immune infiltration
    in some cases — immunotherapy now standard)
  - No ERBB2 amplification (unlike EAC)
  - No Barrett's precursor (unlike EAC)

AVAILABLE PUBLIC DATA:
  TCGA-ESCA:        n=89 ESCC samples (histology: SCC)
                    RNA-seq available — usable for
                    primary ESCC analysis
  GEO datasets:
    GSE53625:       n=119 ESCC + n=179 adjacent normal
                    Largest ESCC dataset with normal
                    reference. CRITICAL for depth axis.
                    Chinese cohort. Affymetrix platform.
    GSE20347:       n=17 ESCC + n=17 adjacent normal
                    Smaller but has paired normal.
    GSE161533:      ESCC with survival annotation
    GSE75241:       Chinese ESCC cohort
  NOTE: Chinese GEO datasets are critical for ESCC
        because ESCC prevalence is highest in China
        and the largest datasets come from Chinese
        institutions. The biology is the same regardless
        of geography — the driver mutations and false
        attractor are consistent.

ANALYSIS STATUS:    NOT YET STARTED
PLANNED DOCUMENT:   ESCA-S1a (ESCC before-document)
PLANNED DATASET:    GSE53625 primary (n=119 + n=179 normal)
                    + TCGA-ESCA ESCC subset (n=89)
PRIORITY:           FIRST — ESCC is the globally dominant
                    form of esophageal cancer by prevalence.
                    GSE53625 provides the highest-powered
                    normal reference of any ESCC dataset.
```

---

## SECTION IV — SUBTYPE 2: ESOPHAGEAL ADENOCARCINOMA (EAC)

```
CLINICAL FACTS:
  Prevalence:       ~10% of esophageal cancer globally
                    ~65–70% of esophageal cancer in
                    Western countries (USA, UK, Australia)
                    Rising incidence: EAC has increased
                    ~600% in incidence in the USA since
                    the 1970s — one of the steepest
                    rising cancer incidences in oncology.
                    The rise parallels the obesity and
                    GERD epidemic in Western populations.
  Geographic peak:  UK, USA, Australia, Northern Europe
  Age at diagnosis: Typically 60–75 years
  Sex:              Male >> Female (7:1) — one of the
                    most male-predominant cancers in
                    oncology. Reason incompletely understood.
                    Adiposity, bile reflux, and sex hormone
                    effects on Barrett's metaplasia are
                    candidate explanations.
  5-year survival:  ~20% overall
                    ~45% (Stage I/II)
                    ~5% (Stage IV)
  Standard of care: Stage I:   Endoscopic resection (if
                    high-grade dysplasia or T1)
                    Stage II/III: Neoadjuvant CRT
                    (CROSS regimen: paclitaxel +
                    carboplatin + RT) then surgery
                    OR perioperative chemotherapy
                    (FLOT: docetaxel + oxaliplatin +
                    5-FU + leucovorin) — FLOT4 trial
                    Stage IV: Platinum + 5-FU ±
                    nivolumab ± trastuzumab (if HER2+)
                    HER2-positive EAC (~15–20%):
                      Trastuzumab + chemotherapy
                      (ToGA trial 2010 — pivotal)
                    PD-1 inhibitors:
                      Nivolumab (CheckMate-649)
                      now standard first-line for
                      CPS ≥5 gastric/GEJ cancers
                      (EAC/GEJ included in this trial)
  pCR rate:         ~25–30% with neoadjuvant CROSS CRT

CELL OF ORIGIN:
  Barrett's epithelial cell — the metaplastic columnar
  cell that has replaced normal squamous epithelium
  in the lower esophagus under chronic acid exposure.

  The transformation sequence (the Waddington pathway):
  Step 1: NORMAL squamous cell at gastroesophageal junction
          Chronic GERD → acid and bile reflux
          → oxidative stress and DNA damage
  Step 2: METAPLASIA — squamous cell is replaced by
          columnar/intestinal-type epithelium (Barrett's)
          CDX2 is activated — intestinal master TF
          This is NOT yet cancer — it is a cell that has
          changed identity (metaplasia)
  Step 3: DYSPLASIA — Barrett's cell acquires TP53 mutation
          Low-grade dysplasia → high-grade dysplasia
          TP53 mutation is the key gating event —
          it occurs at the dysplasia stage, not at
          the Barrett's stage
  Step 4: ADENOCARCINOMA — full false attractor entry
          Copy number amplifications, ERBB2/VEGFA/EGFR
          amplification, CDKN2A loss, SMAD4 loss

  The normal starting point for the depth score is:
  The Barrett's cell (columnar, CDX2+) — NOT the squamous
  cell. The Waddington axis runs from Barrett's (CDX2+,
  non-proliferating, non-invasive) through dysplasia
  to EAC (proliferating, invasive, genomically unstable).

  Normal identity markers of the Barrett's/columnar cell:
    CDX2:           Intestinal master TF — maintained
    MUC2:           Goblet cell mucin — maintained early
    CK8/18:         Columnar cytokeratins
    TFF3:           Trefoil factor — maintained
    HNF4A:          Hepatocyte nuclear factor — expressed
                    in columnar metaplasia
    VILLIN:         Intestinal brush border
    ESR1:           Absent (unlike gastric cardia)

INITIATING MOLECULAR EVENTS:
  TP53 mutation:    ~72% of EAC
                    Occurs at the dysplasia → EAC transition.
                    The key gatekeeping mutation.
  CDKN2A (p16):    ~75% (deletion or methylation)
                    Loss of cell cycle brake
  ERBB2 (HER2):    ~15–20% amplification
                    The most clinically targetable event
                    in EAC. ToGA/TOGA trials established
                    trastuzumab benefit.
  VEGFA:            ~10% amplification
                    Bevacizumab studied but inconsistent
  EGFR:             ~10% amplification
  SMAD4:            ~17% mutation/deletion
                    TGF-β pathway — loss promotes invasion
  ARID1A:           ~9%
  PIK3CA:           ~10%
  KRAS:             ~5%
  High copy number instability (CIN-dominant, like gastric
  CIN and HGSC — large amplicons driving oncogenes)

MOLECULAR SIMILARITY TO OTHER CANCERS:
  EAC shares molecular programme with:
    Gastric adenocarcinoma (CIN subtype)
    Colorectal carcinoma (CMS2, CIN type)
  The TCGA explicitly noted that EAC clusters with
  gastric CIN adenocarcinoma, not with ESCC.
  This reflects a shared origin in glandular/intestinal
  epithelium vs. squamous epithelium.
  The Waddington landscape for EAC is a glandular
  landscape, not a squamous landscape.

KEY DISTINGUISHING FEATURES vs. ESCC:
  - Glandular/adenocarcinoma histology
  - CDX2 positive, p63 negative
  - Lower esophagus / GEJ location
  - Barrett's precursor (identifiable, surveil-able)
  - TP53 + CDKN2A + ERBB2 amplification as drivers
  - CIN-dominant (copy number driven, not point-mutation)
  - HER2+ in ~15–20% — targetable with trastuzumab
  - Immunotherapy active (nivolumab standard in CPS≥5)
  - Male predominance (7:1) far exceeds ESCC (3:1)
  - Rising incidence in Western populations
  - No SOX2/p63 amplification (squamous markers absent)

AVAILABLE PUBLIC DATA:
  TCGA-ESCA:        n=80 EAC samples (histology: AC)
                    Primary dataset — all RNA-seq available
  GEO datasets:
    GSE13898:       EAC + Barrett's + normal squamous
                    CRITICAL: has all three states
                    (normal squamous, Barrett's, EAC)
                    allowing the full Waddington axis
                    to be reconstructed
    GSE26886:       EAC + Barrett's + normal
                    Similar to GSE13898 — another
                    three-state dataset
    GSE37203:       EAC survival annotated
    GSE19417:       Barrett's progression to EAC
  NOTE: The three-state datasets (normal squamous /
        Barrett's / EAC) are the most analytically
        powerful for EAC because they capture the
        full Waddington transition including the
        intermediate Barrett's state. This is rare
        in oncology — most cancers do not have a
        well-defined, capturable intermediate state.
        EAC is an exception.

ANALYSIS STATUS:    NOT YET STARTED
PLANNED DOCUMENT:   ESCA-S2a (EAC before-document)
PLANNED DATASET:    GSE13898 primary (3-state dataset)
                    + TCGA-ESCA EAC subset (n=80)
PRIORITY:           SECOND — after ESCC. The three-state
                    datasets make EAC particularly
                    tractable for the Waddington framework.
                    The Barrett's saddle point is
                    directly capturable from public data.
```

---

## SECTION V — THE MIXED / RARE SUBTYPES

```
ADENOSQUAMOUS CARCINOMA:
  Prevalence:       <3% of esophageal cancers
  Definition:       Contains both squamous AND glandular
                    elements in the same tumor
  Cell of origin:   Unknown / pluripotent progenitor
                    or collision tumor
  Clinical behavior: Generally poor prognosis
  Public data:      Too rare for dedicated analysis
  Status:           Not planned — insufficient data

SMALL CELL CARCINOMA OF THE ESOPHAGUS:
  Prevalence:       <1%
  Cell of origin:   Neuroendocrine cell (rare in esophagus)
  Molecular:        Similar to small cell lung cancer
                    (TP53, RB1, ASCL1)
  Status:           Not planned — insufficient data

UNDIFFERENTIATED CARCINOMA:
  Prevalence:       <1%
  Status:           Not planned — insufficient data

NOTE: These rare types will not receive dedicated analyses
in this repository due to insufficient public data.
They are noted here for completeness only.
```

---

## SECTION VI — THE TWO-DISEASE STRUCTURAL PICTURE

```
The most important structural fact about esophageal cancer
for this framework is this:

  ESCC and EAC are not subtypes of the same Waddington
  landscape. They are two separate Waddington landscapes
  that happen to be located in the same organ.

  ESCC landscape:
    Normal cell:    Squamous progenitor (p63+, SOX2+)
    Saddle point:   NOTCH1 loss + TP53 mutation in
                    basal squamous cell
    False attractor: Undifferentiated squamous proliferative
                    state (p63+, SOX2 amplified, cycling)
    Depth axis:     Loss of squamous differentiation
                    programme (ZNF750, NOTCH targets) rising
                    with depth; proliferation markers
                    (MKI67, PCNA) rising with depth

  EAC landscape:
    Normal cell:    Barrett's columnar cell (CDX2+)
    Saddle point:   TP53 mutation in Barrett's dysplasia
    False attractor: Glandular proliferative state with
                    CIN and RTK amplification
                    (ERBB2+, CDX2 retained or lost,
                    intestinal identity partially preserved)
    Depth axis:     Loss of Barrett's identity programme
                    (CDX2, TFF3, MUC2) with increasing depth;
                    RTK/proliferation markers rising

  THE KEY STRUCTURAL DIFFERENCE:
    ESCC: the false attractor RETAINS some normal identity
    markers (p63, SOX2) — it is a blocked differentiation,
    not a complete identity replacement.

    EAC: the false attractor LOSES some of the Barrett's
    identity markers as it deepens — it becomes less
    intestinal and more undifferentiated with depth.

    This is the same pattern seen in RCC:
    ccRCC retains some proximal tubule markers while PRCC
    loses them more completely with depth.

  CROSS-SUBTYPE COMPARISON:
    After ESCC and EAC individual analyses are complete,
    the cross-subtype comparison will ask:
      1. Is there a gene that is depth-positive in BOTH
         subtypes? (Universal esophageal depth marker)
      2. Is EZH2 the epigenetic lock in both?
         (Prediction: YES in EAC, UNCERTAIN in ESCC)
      3. Does the depth score predict platinum response
         in both subtypes?
      4. What is the structural analogue of the
         "chRCC inversion" in the esophageal landscape?
         (Prediction: EAC has an inverted depth axis
         relative to ESCC because normal Barrett's is
         already a "metaplastic" state — the depth axis
         direction will need careful definition)
```

---

## SECTION VII — THE BARRETT'S ESOPHAGUS OPPORTUNITY

```
Barrett's esophagus is the single most important
intermediate state in this entire repository.

In most cancers, the Waddington saddle point is a
theoretical construct — we can infer it from the gene
expression data at the extremes (normal vs. deep tumour)
but we cannot directly observe cells sitting at the
saddle point.

In EAC, we can.

Barrett's esophagus is:
  - Clinically detectable (endoscopy + biopsy)
  - Biologically defined (CDX2+, MUC2+, columnar)
  - Surveil-able (biopsies taken every 3-5 years)
  - Available in public datasets (GSE13898, GSE26886)
  - The literal saddle point cell between normal
    squamous and EAC false attractor

This means that for EAC:
  - The depth score has THREE reference points,
    not two: normal squamous / Barrett's / EAC
  - The transition from Barrett's to EAC is the
    clinically relevant window (intervention window)
  - Genes that are elevated in Barrett's but rise
    further in EAC are depth markers in the
    Barrett's → EAC transition
  - Genes that are elevated in Barrett's but fall
    in EAC may represent protective identity genes
    that are lost as the cancer deepens

The Barrett's saddle point analysis is the most
clinically actionable output of the EAC analysis:
  If depth score can be derived from Barrett's
  biopsies, it would predict progression risk.
  This is a clinical need that is currently unmet.
  Current surveillance relies on histological grading
  of dysplasia — a crude, subjective, and imperfect
  predictor of EAC progression.
  A gene expression depth score in Barrett's biopsies
  would be a quantitative, continuous risk predictor
  derived from first principles.

This is stated here as a structural opportunity,
not as a prediction for the before-document.
The before-document (ESCA-S2a) will formalize
the predictions before any data loads.
```

---

## SECTION VIII — THE EXISTING ESCA ANALYSIS — CONTEXT

```
The existing analysis in Cancer_Research/ESCA/ ran on
the full TCGA-ESCA dataset (ESCC + EAC combined, n~173).

What that analysis most likely captured:
  The depth axis that separates BOTH subtypes from their
  respective normal cells simultaneously. Because ESCC and
  EAC have different normal cells (squamous vs. Barrett's),
  the combined depth score is not a clean axis for either
  subtype individually.

  The genes most likely to have appeared as depth-positive
  in the combined analysis are those shared between both
  false attractors — likely:
    EZH2 (epigenetic lock — likely in both)
    MKI67, PCNA, TOP2A (proliferation — universal)
    MYC (transcriptional amplifier — likely in both)

  The genes most likely to have appeared as depth-negative
  are those that are normal identity markers for BOTH:
    Squamous markers (KRT4, KRT13, IVL — normal squamous)
    Columnar markers (CDX2, TFF3, MUC2 — normal Barrett's)
    Both sets would appear negative if the combined
    depth score penalizes any normal identity loss

  The ESCC-specific and EAC-specific depth structure
  will only be visible when the two subtypes are
  analyzed separately.

This is not a criticism of the existing analysis —
it is a statement of what the existing analysis
is and what it is not.
The subtype analyses will decompose it.
```

---

## SECTION IX — DATA AVAILABILITY SUMMARY

```
Subtype   Primary Dataset        n (tumour)  Normal Ref         Power
──────────────────────────────────────────────────────────────────────
ESCC      GSE53625               n=119       n=179 adj normal   HIGH
          + TCGA-ESCA SCC        n=89        TCGA adj normal    MOD
EAC       GSE13898               EAC+BE+Sq   3-state dataset    MOD
          + GSE26886             EAC+BE+Sq   3-state dataset    MOD
          + TCGA-ESCA AC         n=80        TCGA adj normal    MOD
Mixed     TCGA-ESCA (combined)   n=173       TCGA adj normal    HIGH
          (existing analysis)

LEGEND: BE = Barrett's Esophagus, Sq = normal squamous

KEY NOTE ON NORMAL REFERENCE:
  For ESCC: normal adjacent squamous esophagus is the
            reference. GSE53625 provides n=179 normal
            adjacent samples — the largest normal
            reference for any cancer in this repository.
            This is exceptional.

  For EAC:  Two options:
            1. Normal squamous (what the tissue USED to be)
            2. Barrett's (what the immediate precursor is)
            The depth axis will be defined relative to
            BARRETT'S as the reference state — because
            Barrett's is the cell that gives rise to EAC,
            not the squamous cell.
            This is a critical analytical decision that
            must be stated in the before-document.
```

---

## SECTION X — PLANNED ANALYSIS ORDER

```
ORDER:

  ESCA-S1   ESCC     GSE53625 primary         HIGH POWER
                     REASON: ESCC is the globally dominant
                     form. GSE53625 provides the best normal
                     reference of any ESCA dataset (n=179).
                     Squamous cell of origin is the cleaner
                     starting point for depth score derivation
                     because the normal cell is well-defined
                     and the false attractor is a
                     blocked-differentiation state, not a
                     metaplastic state.

  ESCA-S2   EAC      GSE13898 / GSE26886      MODERATE
                     REASON: The three-state datasets
                     (squamous / Barrett's / EAC) make EAC
                     the most analytically rich cancer for
                     the Waddington framework. The Barrett's
                     saddle point is directly observable.
                     The depth score in Barrett's biopsies
                     is the clinical output.

  ESCA-X    Cross-   After ESCC and EAC complete
            subtype  Questions:
                       1. Universal esophageal depth genes?
                       2. EZH2 in both? Same direction?
                       3. Shared drug targets (basket trial)?
                       4. What does the combined depth
                          score (existing analysis) actually
                          represent structurally?
```

---

## SECTION XI — WHAT THIS DOCUMENT DOES NOT CONTAIN

```
This document contains:
  ✓ Subtype definitions (ESCC and EAC)
  ✓ Cells of origin (squamous progenitor vs. Barrett's)
  ✓ Published molecular drivers (literature, not prediction)
  ✓ Clinical characteristics (survival, standard of care)
  ✓ Normal cell identity programmes
  ✓ The Waddington transition sequence for each subtype
  ✓ Data availability (GEO accessions, sample sizes)
  ✓ Context for the existing combined analysis
  ✓ The Barrett's saddle point structural opportunity

This document does NOT contain:
  ✗ Depth score predictions
  ✗ False attractor gene predictions (specific)
  ✗ Drug target predictions
  ✗ Epigenetic lock mechanism hypotheses
  ✗ Cross-subtype structural predictions beyond literature

All of the above belong in the BEFORE documents.
ESCA-S1a (ESCC before-document) is next.
It will be written before any script runs and before
any data is loaded.
```

---

## STATUS BLOCK

```
document:           ESCA_Subtype_Orientation.md
folder:             Cancer_Research/ESCA/Subtypes/
status:             COMPLETE — ORIENTATION ONLY
date:               2026-03-03
author:             Eric Robert Lawson
                    OrganismCore

subtypes_covered:
  primary:          ESCC, EAC                        [2 of 2]
  rare/mixed:       Adenosquamous, Small cell noted  [noted]

analyses_started:   0

existing_analysis:  Cancer_Research/ESCA/ — conducted on
                    combined TCGA-ESCA (ESCC + EAC mixed).
                    Valid as a combined landscape analysis.
                    Not a subtype-specific analysis.
                    Decomposed by the analyses in this folder.

next_document:      ESCA-S1a
                    ESCC Before-Document
                    (predictions locked before
                    GSE53625 / TCGA-ESCA ESCC loads)

critical_note:      ESCC and EAC are not subtypes of
                    the same disease.
                    They are two different diseases
                    in the same organ.
                    They will be analyzed as separate
                    Waddington landscapes.
                    Their cross-subtype comparison comes
                    after both individual analyses are
                    complete.
```
