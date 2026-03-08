# FELINE MAMMARY CARCINOMA — IHC DEPTH SCORE PROTOCOL
## FOXA1 / EZH2 / H3K27me3 Panel
## Two-Pathway Protocol: Luminal FMC and Triple-Negative FMC
## OrganismCore — Eric Robert Lawson
## 2026-03-08

---

## STATUS: ACTIVE — REQUIRED BEFORE TAZEMETOSTAT

```
THIS PROTOCOL BRANCHES AT STEP 1.

The IHC panel determines which
geometric pathway applies to this
patient before any depth score
calculation.

STEP 1 — SUBTYPE DETERMINATION:
  Run ER, PR, HER2 (standard
  diagnostic panel — likely already
  done at diagnosis).
  Add FOXA1 to this panel.
  This determines which geometry applies.

PATHWAY A — LUMINAL FMC
(ER+ and/or FOXA1+ on diagnostic panel):
  Use FOXA1/EZH2 depth score.
  R = FOXA1 H-score / EZH2 H-score.
  Treatment: tazemetostat +
  endocrine therapy sequence.

PATHWAY B — TRIPLE-NEGATIVE FMC
(ER-, PR-, HER2-, FOXA1-):
  Use EZH2 H-score as depth measure.
  EZH2 + H3K27me3 panel.
  Treatment: tazemetostat +/-
  doxorubicin.

CRITICAL: FOXA1 ANTIBODY FELINE
VALIDATION IS MANDATORY.

  No commercial FOXA1 antibody is
  formally validated for feline
  FFPE tissue. Before using FOXA1
  IHC results to guide clinical
  decisions in cats:

  MANDATORY IN-HOUSE VALIDATION:
  Run the chosen FOXA1 antibody clone
  on feline mammary tissue from an
  archived case with known ER status.
  An ER-positive feline mammary tumour
  should show nuclear FOXA1 positivity
  in luminal cells.
  An ER-negative triple-negative FMC
  should show absent or very low FOXA1.
  Normal feline mammary gland: luminal
  cells should show FOXA1 positivity.
  Myoepithelial/basal cells should
  be FOXA1-negative.

  If staining pattern is consistent
  with this expectation: proceed.
  If staining is absent across all
  tissues including normal: antibody
  is not working in feline tissue.
  Try alternative clone.
  Do not report a depth score from
  unvalidated staining.

CRITICAL: FELINE DRUG METABOLISM.

  See dosing section. Cats are not
  dogs. CYP3A metabolises tazemetostat
  but feline CYP3A activity is not
  quantitatively equivalent to canine.
  No published feline PK data.
  Conservative starting doses.
  Mandatory feline oncology specialist
  consultation before prescribing.
  More frequent monitoring than dogs.

SURGERY PROCEEDS FIRST:
  Do not delay mastectomy for IHC
  depth score results. The depth
  score guides adjuvant and
  maintenance therapy decisions,
  not the surgical decision.
  Run the panel on the surgical
  specimen in parallel with the
  routine histopathology.
  Results guide the systemic therapy
  plan at the post-surgical consultation.
```

---

## STEP 1 — SUBTYPE DETERMINATION PANEL

```
REQUIRED STAINS FOR SUBTYPE
DETERMINATION:
  ER (estrogen receptor)
  PR (progesterone receptor)
  HER2
  FOXA1

  ER and PR are standard feline FMC
  diagnostic stains in most reference
  labs. HER2 is increasingly included.
  FOXA1 is the addition for this protocol.

  Feline ER antibodies:
    Clone 6F11 (Leica/Novocastra):
    widely used and validated for
    feline mammary tissue in the
    veterinary literature.
    Clone SP1 (rabbit monoclonal):
    also used in feline FMC studies.
    Use whichever is validated at
    the performing lab.

  Feline FOXA1 — ANTIBODY CANDIDATES
  (in-house validation required):

    Option 1: D4E2 (CST #53528)
    Rabbit monoclonal. The same clone
    used in the BRCA IHC calibration
    protocol. Validated for human
    and canine tissue. Feline
    cross-reactivity expected (high
    FOXA1 sequence conservation) but
    not commercially confirmed.
    VALIDATION REQUIRED.
    Dilution (human/canine): 1:400–1:800.
    Antigen retrieval: HIER, citrate
    pH 6.0, 20 minutes.
    Try this dilution range for feline.

    Option 2: Clone 3A8
    (Thermo Fisher MA1-091)
    Mouse monoclonal. Validated for
    human, mouse, dog. No feline data.
    Try as second option if D4E2
    fails feline validation.
    Dilution: 1:200–1:500.

  FOXA1 in-house validation controls:
    Positive: feline ER+ mammary
    carcinoma (archived FFPE block).
    Negative: feline TNBC mammary
    carcinoma (archived FFPE block).
    Species control: normal feline
    mammary gland tissue (luminal
    cells should be FOXA1+,
    myoepithelial cells FOXA1-).

SUBTYPE CLASSIFICATION RESULT:

  ☐ LUMINAL FMC:
    ER+ (any staining > 1%) OR
    FOXA1+ (any nuclear staining
    ≥ 5% of neoplastic cells
    on validated staining).
    → Proceed to PATHWAY A.

  ☐ TRIPLE-NEGATIVE FMC:
    ER < 1% AND PR < 1% AND
    HER2 score 0–1+ AND
    FOXA1 < 5% of neoplastic cells.
    → Proceed to PATHWAY B.

  ☐ HER2-POSITIVE FMC:
    HER2 score 3+ or 2+ with FISH
    amplification.
    → HER2-targeted therapy primary.
    → Also run FOXA1/EZH2 depth score
    for maintenance strategy.
    Contact OrganismCore@proton.me.

  ☐ LUMINAL-AR SUBTYPE
  (within TNBC-FMC):
    ER < 1%, PR < 1%, HER2 0–1+
    BUT FOXA1 ≥ 5% (despite no ER).
    This is the Dagher 2019 "luminal-AR"
    subtype. Best prognosis within
    TNBC-FMC. Use FOXA1/EZH2 depth
    score (Pathway A) even without
    ER expression — the FOXA1 programme
    is partially active.
    Treatment: tazemetostat +
    androgen receptor-directed therapy
    if AR+ is also confirmed.
    Contact OrganismCore@proton.me.
```

---

## PATHWAY A — LUMINAL FMC
## FOXA1 / EZH2 DEPTH SCORE

```
STAINS REQUIRED FOR DEPTH SCORE
(PATHWAY A):
  FOXA1 — Identity Anchor
  EZH2 — Convergence Hub protein
  H3K27me3 — Hub activity readout
              and safety gate

PROTEIN 1 — FOXA1
(Forkhead Box A1 / HNF3-alpha)
  Nuclear staining.
  Score neoplastic luminal cells only.
  Exclude myoepithelial/basal cells
  (these are FOXA1-negative in normal
  tissue — use as internal negative
  control for FOXA1 specificity).

  ANTIBODY: D4E2 (CST #53528)
  OR clone 3A8 (Thermo Fisher MA1-091)
  — WHICHEVER PASSED IN-HOUSE FELINE
  VALIDATION.
  Dilution: as validated.
  Antigen retrieval: HIER, citrate
  pH 6.0, 20 minutes.

  HIGH H-SCORE (> 150):
    FOXA1 programme substantially active.
    Shallower luminal lock.
    R likely > 1.0.
    Endocrine therapy may be partially
    effective already. Tazemetostat
    as maintenance, not urgent first-line.

  INTERMEDIATE (80–150):
    FOXA1 partially suppressed.
    Luminal lock intermediate.
    R depends on EZH2 level.

  LOW (< 80):
    FOXA1 substantially suppressed.
    Deep luminal lock.
    R likely < 0.5.
    Tazemetostat urgent.
    Endocrine therapy after lock
    dissolution.

  DIAGNOSTIC NOTE — MELANIN:
    N/A for mammary tissue.
    Lipofuscin granules can occasionally
    appear in aged feline mammary tissue
    — these are granular cytoplasmic
    deposits, not nuclear, and do not
    interfere with nuclear H-scoring.

PROTEIN 2 — EZH2
(Enhancer of Zeste Homolog 2)
  Nuclear staining.
  Score neoplastic cells only.

  ANTIBODY: D2H1 (CST #5246).
  Rabbit monoclonal. Validated in
  multiple species including feline
  mammary tissue (cross-reactivity
  confirmed — H3K27me3 is a
  conserved PTM; EZH2 protein itself
  is highly conserved across vertebrates).
  Dilution: 1:400–1:800.
  Antigen retrieval: HIER, citrate
  pH 6.0, 20 minutes.

  HIGH H-SCORE (> 150):
    EZH2 overexpressed. Hub candidate.
    Confirm with H3K27me3.

  INTERMEDIATE (80–150):
    EZH2 present. Assess H3K27me3.

  LOW (< 80):
    EZH2 not dominant. Consider whether
    EZH2 is actually the hub in this
    case. Contact OrganismCore@proton.me.

PROTEIN 3 — H3K27me3
  THE SAFETY GATE AND HUB ACTIVITY
  READOUT. Same function as in all
  other protocols in this series.

  ANTIBODY: C36B11 (CST #9733).
  Rabbit monoclonal. Fully conserved
  PTM. Cross-species validated.

  INTERNAL POSITIVE CONTROL:
    Stromal fibroblasts, lymphocytes,
    endothelial cells in the same
    section must show H3K27me3
    positivity. If all cells including
    stroma are negative: TECHNICAL
    FAILURE. Repeat.

  H3K27me3 HIGH (> 100, > 75% positive,
  stromal control positive):
    Hub ACTIVE. Proceed to depth score.

  H3K27me3 LOW (< 60, < 25% positive,
  stromal control positive):
    Hub PARALYSED. Do not use
    tazemetostat. Rare in FMC but
    check. Contact OrganismCore@proton.me.

  H3K27me3 ALL CELLS NEGATIVE
  (including stroma):
    Technical failure. Repeat stain.

DEPTH SCORE CALCULATION (PATHWAY A):
  SAFETY GATE FIRST.
  H3K27me3 > 100 with stromal
  control positive? → PROCEED.

  R = FOXA1 H-score / EZH2 H-score.

  R > 2.0 — SHALLOW (Luminal A-like):
    MITF dominant (FOXA1 dominant).
    Standard surgery + OVH (if intact)
    + endocrine therapy (toremifene
    if available).
    Tazemetostat: reserve for relapse
    or if response to endocrine therapy
    is incomplete.

  R 1.0–2.0 — INTERMEDIATE (Luminal B-like):
    Moderate EZH2 lock.
    Surgery + OVH + endocrine therapy.
    Tazemetostat: add as adjuvant
    post-surgical with endocrine therapy.
    Monitor ER target gene activity
    (Ki-67 reduction is first measurable
    signal of EZH2 lock dissolution).

  R 0.5–1.0 — DEEP (Luminal B, deep lock):
    EZH2 approaching dominance.
    Surgery + OVH.
    Tazemetostat FIRST (2–4 weeks
    after surgery).
    Endocrine therapy AFTER tazemetostat
    has been running (lock must be
    partially dissolved before endocrine
    therapy has target).
    Monitor at month 1 for any reduction
    in tumour Ki-67 on repeat biopsy
    if available (not mandatory — clinical
    signs are primary endpoint).

  R < 0.5 — VERY DEEP:
    EZH2 dominant. FOXA1 heavily
    suppressed despite ER positivity.
    Endocrine therapy alone likely
    ineffective.
    Surgery + OVH.
    Tazemetostat immediately post-surgical.
    Endocrine therapy deferred until
    tazemetostat response established
    (month 2–3).
    Doxorubicin chemotherapy may be
    added for high-grade tumours.

SCORING (PATHWAY A):
  H-score: standard method.
  Score neoplastic luminal epithelial
  cells only.
  Exclude basal/myoepithelial cells
  (FOXA1-negative normally — provides
  internal antibody specificity check).
  Minimum 200 neoplastic cells per stain.
  Two observers, blinded, ICC > 0.7.
```

---

## PATHWAY B — TRIPLE-NEGATIVE FMC
## EZH2 / H3K27me3 DEPTH SCORE

```
STAINS REQUIRED FOR DEPTH SCORE
(PATHWAY B):
  EZH2 — Convergence Hub protein
  H3K27me3 — Hub activity readout
              and safety gate

  FOXA1 staining is not informative
  in confirmed TNBC-FMC (expected
  absent or < 5%). Run EZH2 and
  H3K27me3 only for the depth score.

PROTEIN 1 — EZH2
  Same antibody and protocol as
  Pathway A. (D2H1, CST #5246.)

  In TNBC-FMC, EZH2 is typically
  higher than in luminal FMC.
  EZH2 > 200 H-score is not uncommon
  in aggressive TNBC-FMC.

PROTEIN 2 — H3K27me3
  Same antibody and protocol as
  Pathway A. (C36B11, CST #9733.)
  Stromal internal control mandatory.

DEPTH CLASSIFICATION (PATHWAY B):

  SAFETY GATE:
  H3K27me3 > 100 with stromal
  control positive → hub active →
  PROCEED.
  H3K27me3 < 60 with stromal
  control positive → hub paralysed
  → DO NOT USE TAZEMETOSTAT.
  Contact OrganismCore@proton.me.

  DEPTH SCORE (hub confirmed):

  EZH2 > 200 AND H3K27me3 > 150:
    VERY DEEP.
    Maximum EZH2 dependence.
    Tazemetostat highest priority.
    Combine with doxorubicin if
    cat is a good surgical candidate
    and performance status allows.

  EZH2 150–200 AND H3K27me3 100–150:
    DEEP.
    Tazemetostat recommended.
    +/- doxorubicin based on
    clinical status and grade.

  EZH2 100–150 AND H3K27me3 60–100:
    INTERMEDIATE.
    Tazemetostat may benefit.
    Monitor carefully.
    Doxorubicin if high grade.

  EZH2 < 100 OR H3K27me3 < 60:
    LOW EZH2 ACTIVITY.
    EZH2 may not be the primary hub.
    Do not use tazemetostat without
    further investigation.
    Contact OrganismCore@proton.me.
    Consider standard doxorubicin
    chemotherapy as primary.
```

---

## TAZEMETOSTAT DOSING GUIDANCE
## (FELINE — OFF-LABEL, HIGHLY CONSERVATIVE)

```
CRITICAL PREAMBLE — READ BEFORE
PRESCRIBING:

  Cats are NOT small dogs.
  Feline drug metabolism is
  fundamentally different from
  canine and human:

    UGT glucuronidation: DEFICIENT
    in cats. Tazemetostat is CYP3A
    metabolised (not primarily
    glucuronidated) — this is
    favourable. But UGT deficiency
    affects other drugs the cat may
    be receiving concurrently.

    CYP3A in cats: present but
    activity not quantitatively
    comparable to human or dog.
    No published data on feline
    CYP3A activity for tazemetostat
    specifically.

    Half-life prediction: unknown.
    May be longer than in humans/dogs
    due to different CYP3A activity.
    This means accumulation risk at
    human/canine equivalent doses.

    Starting philosophy:
    START LOW. ESCALATE SLOWLY.
    Monitor intensively.
    The first feline cases are
    establishing the safety floor
    for the species.
    Document every adverse event.

  MANDATORY: Consult a veterinary
  oncologist with feline
  pharmacology experience before
  prescribing. This is not optional.
  Off-label use requires full
  informed owner consent.

HUMAN REFERENCE DOSE:
  800 mg PO BID (~70 kg adult).
  ~23 mg/kg/day or ~960 mg/m²/day.

FELINE STARTING DOSE RECOMMENDATION:
  200–300 mg/m² PO BID.
  (Approximately 20–30% of the human
  BSA dose as the initial feline
  starting point — much more
  conservative than the canine
  50–60% starting range, reflecting
  unknown feline CYP3A kinetics.)

FELINE BSA FORMULA:
  BSA (m²) = (body weight kg)^0.67
  × 0.100
  (Modified Meeh formula for cats.)

WEIGHT-BASED EXAMPLES:
  3 kg cat (BSA ≈ 0.21 m²):
    42–63 mg PO BID.
    Compounded liquid or capsule
    required — no commercial
    formulation for this dose.
    Compounding pharmacy essential.

  4 kg cat (BSA ≈ 0.26 m²):
    52–78 mg PO BID.
    Compounded formulation.

  5 kg cat (BSA ≈ 0.31 m²):
    62–93 mg PO BID.
    Compounded formulation.

  6 kg cat (BSA ≈ 0.36 m²):
    72–108 mg PO BID.
    Compounded formulation.

  All feline doses require compounding.
  200mg tablets cannot be used directly
  for any typical cat weight at the
  conservative starting dose range.
  Compounding as oral liquid or
  gelatin capsule is required.

ADMINISTRATION:
  With food. pH-dependent drug.
  Feline appetite unpredictability
  is a risk — administer with a
  small amount of palatable food
  (tuna broth, chicken broth) to
  ensure ingestion with gastric acid
  present.
  NEVER administer tazemetostat with
  omeprazole or famotidine — these
  reduce gastric acidity and may
  impair absorption.

DOSE ESCALATION:
  If well-tolerated at starting dose
  after 4 weeks of monitoring:
  may increase to 300–400 mg/m²
  PO BID.
  Maximum recommended escalation:
  400 mg/m² PO BID (< 50% of
  human BSA dose) until feline
  PK data exists.

DRUG INTERACTIONS IN CATS
(more complex than dogs):

  CYP3A inhibitors (AVOID):
    Ketoconazole — common in cats
    for dermatophytosis. Potent
    CYP3A inhibitor. If cat is on
    ketoconazole for any reason:
    do NOT add tazemetostat until
    ketoconazole is discontinued
    and cleared.
    Itraconazole — common in cats
    for fungal infections. CYP3A
    inhibitor.
    Fluconazole — moderate inhibitor.
    Erythromycin — used in cats as
    prokinetic. CYP3A inhibitor.
    Cyclosporine — used in cats for
    immune-mediated conditions.
    CYP3A substrate AND inhibitor.
    MUST be discontinued before
    tazemetostat.

  CYP3A inducers (REDUCE EFFICACY):
    Phenobarbital — common in feline
    epilepsy. Induces CYP3A.
    Rifampicin — occasionally used.
    If cat is on phenobarbital:
    tazemetostat may be rapidly
    metabolised and have insufficient
    exposure for efficacy. Discuss
    with oncologist.

  GLUCURONIDATION-AFFECTED DRUGS:
    (These are NOT interactions
    with tazemetostat itself but are
    dangerous in cats and must be
    monitored if concurrent.)
    NSAIDs: piroxicam is commonly
    used for anti-tumour effect in
    cats. Low-dose only. Cats have
    reduced ability to glucuronidate
    NSAIDs. Monitor renal function.
    Meloxicam: licensed for cats
    at very low doses. Use at
    approved feline doses only.
    Acetaminophen (paracetamol):
    ABSOLUTELY CONTRAINDICATED
    in cats. Toxic at any dose.
    Ensure no owner administers
    acetaminophen.

MONITORING SCHEDULE (MORE
INTENSIVE THAN CANINE):

  Baseline: CBC + differential,
  chemistry (ALT, ALP, AST, GGT,
  bilirubin), BUN, creatinine,
  electrolytes, total protein,
  albumin. Body weight.

  Week 1: CBC + ALT + weight.
  Week 2: Full CBC + chemistry.
  Week 4: Full CBC + chemistry.
  Month 2: Full CBC + chemistry.
  Then monthly.

  ALSO MONITOR:
  Body weight at every visit.
  Cats lose weight rapidly with
  illness. Weight loss > 10%
  from baseline: reduce dose
  and investigate.
  Appetite diary kept by owner.
  QoL scored by owner weekly.

ABORT CRITERIA (stop tazemetostat
immediately in cats):
  ALT > 3× ULN (lower threshold
  than dogs — cats have less hepatic
  reserve for novel drugs).
  Any bilirubin elevation above ULN.
  Neutrophils < 1,000/μL.
  Platelets < 75,000/μL.
  Weight loss > 15% from baseline.
  Persistent vomiting (> 2×/day
  for > 3 consecutive days).
  Complete anorexia (not eating
  for > 2 days).
  Any neurological signs.
  Marked lethargy persisting > 48h.
  Owner distress regarding quality
  of remaining life.
```

---

## IHC REPORTING TEMPLATE

```
PATIENT:
  Species: Cat
  Breed: _______________
  Age: _____ years
  Sex: ☐ FN  ☐ FI  ☐ MN  ☐ MI
  Weight: _____ kg
  BSA: _____ m²
    (= weight^0.67 × 0.100)
  Case ID: _____________
  Intact at diagnosis: ☐ Y  ☐ N
  OVH performed: ☐ Y  ☐ N
    Date: ________________

DIAGNOSIS:
  Confirmed FMC:
  ☐ Histopathology
    Gland(s) biopsied: __________
    Date: ____________________
    Grade: ☐ I  ☐ II  ☐ III
  ☐ Multiple glands involved
    Worst grade gland: __________

DIAGNOSTIC PANEL (SUBTYPE):
  ER: _____ % positive / H-score ___
  PR: _____ % positive
  HER2 score: ______
  FOXA1: ☐ Validated for feline tissue
          ☐ In-house validation
            performed — result:
            ☐ Nuclear staining in
              luminal cells of
              positive control:
              ☐ YES → proceed
              ☐ NO → antibody not
                working. Switch clone.
  FOXA1 H-score: _____ / 300
  (Only report if validated.)

SUBTYPE CLASSIFICATION:
  ☐ LUMINAL (ER+ or FOXA1+ ≥ 5%)
    → PATHWAY A
  ☐ TRIPLE-NEGATIVE (ER<1%, PR<1%,
    HER2 0-1+, FOXA1 < 5%)
    → PATHWAY B
  ☐ HER2-POSITIVE → Contact
    OrganismCore@proton.me
  ☐ LUMINAL-AR (TNBC but FOXA1+):
    → Pathway A variant. Contact
    OrganismCore@proton.me.

PATHWAY A RESULTS (if luminal):
  EZH2 H-score: _____ / 300
  H3K27me3 H-score: _____ / 300
  H3K27me3 % positive: _____%
  Stromal H3K27me3 control:
  ☐ Positive (normal) → proceed
  ☐ Negative (technical failure) →
    repeat stain
  Safety gate: ☐ PASS  ☐ FAIL
  R = FOXA1 / EZH2 = _______
  Depth:
  ☐ Shallow (R > 2.0)
  ☐ Intermediate (R 1.0–2.0)
  ☐ Deep (R 0.5–1.0)
  ☐ Very Deep (R < 0.5)

PATHWAY B RESULTS (if TNBC):
  EZH2 H-score: _____ / 300
  H3K27me3 H-score: _____ / 300
  H3K27me3 % positive: _____%
  Stromal control: ☐ POS  ☐ NEG
  Safety gate: ☐ PASS  ☐ FAIL
  EZH2 depth:
  ☐ Very Deep (EZH2 > 200, H3K27me3
    > 150)
  ☐ Deep (EZH2 150–200, H3K27me3
    100–150)
  ☐ Intermediate (EZH2 100–150,
    H3K27me3 60–100)
  ☐ Low EZH2 activity — contact
    OrganismCore@proton.me

Pathologist: ___________________
Second scorer: _________________
ICC (FOXA1 if run): ___
ICC (EZH2): ___  (H3K27me3): ___
Date: ___________________________

TREATMENT PLAN:
  Surgery: ☐ Completed  ☐ Planned
    ☐ Not feasible
  OVH: ☐ Completed  ☐ Planned
    ☐ Already spayed
  Tazemetostat:
    ☐ Planned — starting dose:
      _____ mg PO BID
      _____ mg/m²/day
      Start date: ___________
    ☐ CONTRAINDICATED —
      H3K27me3 safety gate failed
    ☐ On hold — FOXA1 validation
      not yet completed
  Endocrine therapy (luminal only):
    ☐ Toremifene — dose: _______
    ☐ After tazemetostat initiation
    ☐ Not applicable (TNBC)
  Doxorubicin:
    ☐ Planned (TNBC deep/very deep
      + good performance status)
    ☐ Not planned
  Veterinary oncologist consulted:
    ☐ Yes — name: _______________
    ☐ No — OBTAIN BEFORE PRESCRIBING

Supervising veterinarian: ________
Date: ___________________________
```

---

## DOCUMENT METADATA

```
document_id:
  FELINE_FMC_IHC_DEPTH_PROTOCOL

type:
  Patient selection IHC protocol —
  two-pathway: luminal and TNBC

version: 1.0
date: 2026-03-08
status: ACTIVE

author:
  Eric Robert Lawson / OrganismCore
ORCID: 0009-0002-0414-6544
contact: OrganismCore@proton.me

critical_warnings:
  1. FOXA1 feline validation:
     No commercial clone validated
     for feline tissue. In-house
     validation MANDATORY before
     clinical use.

  2. Feline drug metabolism:
     Cats are NOT dogs. Conservative
     starting doses (20–30% of human
     BSA equivalent, NOT 50–60% as
     for dogs). Veterinary oncologist
     mandatory. All feline doses
     require compounding.

  3. Tazemetostat abort threshold
     for ALT in cats is 3× ULN
     (NOT 5× as in dogs). Lower
     threshold reflects reduced
     feline hepatic metabolic reserve
     for novel compounds.

  4. Ketoconazole, itraconazole,
     cyclosporine, erythromycin:
     all common in feline medicine.
     All CYP3A inhibitors. All must
     be discontinued before
     tazemetostat.

related_documents:
  FELINE_FMC_GEOMETRIC_DERIVATION.md
  FELINE_FMC_CASE_REPORT_TEMPLATE.md
  CANINE_DLBCL_IHC_DEPTH_PROTOCOL.md
    (reference for H-score methodology)
```
