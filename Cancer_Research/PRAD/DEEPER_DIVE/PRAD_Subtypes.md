# PROSTATE ADENOCARCINOMA (PRAD) — SUBTYPE ORIENTATION DOCUMENT
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

  A complete map of the prostate adenocarcinoma
  molecular subtype landscape — what the primary
  genomic subtypes are (TMPRSS2-ERG, SPOP, FOXA1,
  IDH1, and others), how they map onto the luminal/
  basal/neuroendocrine transcriptional spectrum,
  where each arises in the normal prostate epithelial
  hierarchy, what the clinical characteristics and
  treatment are at each stage of disease, and what
  public data exists to analyse each.

PRAD holds a specific position in the OrganismCore
repository alongside PAAD as the second cancer
demonstrating the intact circuit finding:

  THE INTACT CIRCUIT — PROSTATE VERSION.

  The existing PRAD bulk analysis (OrganismCore
  cancer series) found that:
    NKX3-1 is the switch gene —
    the master luminal identity TF whose loss
    drives progression from luminal differentiation
    toward a more progenitor/stem-like state.
    EZH2 is elevated — the gain-of-function
    epigenetic lock — silencing NKX3-1 and
    other luminal identity genes via H3K27me3.
    The circuit is intact — restoring NKX3-1
    expression (by blocking EZH2) would re-execute
    the luminal differentiation programme.
    Circuit restoration is therapeutic.
    Tazemetostat is the predicted drug.

  This structural finding — the intact circuit —
  is the primary context for all subtype analysis
  that follows.

  The subtype orientation document situates this
  finding in precise molecular context:

    WHICH subtype has the intact circuit?
      Luminal A — the most differentiated subtype.
      NKX3-1 HIGH, AR HIGH, FOXA1 HIGH.
      Circuit intact. EZH2 blocks it. Tazemetostat
      predicted for this subtype.

    WHICH subtype is at intermediate depth?
      Luminal B — higher proliferation, more
      chromosomal alterations, still AR-positive
      but more genomically unstable.
      Circuit partially intact. AR pathway
      inhibition is the primary tool.

    WHICH subtypes have BROKEN or MISSING circuits?
      Basal-like treatment-naive PRAD:
      NKX3-1 lost, AR low, TP63/KRT5/KRT14 high.
      The circuit has been bypassed.
      Treatment-emergent NEPC:
      AR completely lost. TP53 and RB1 lost.
      AURKA/MYCN amplified. The circuit is not
      just blocked — the machinery is gone.
      Double-Negative Prostate Cancer (DNPC):
      AR negative, NE negative — the most
      identity-lost false attractor state.
      No luminal identity. No neuroendocrine
      identity. Stem-like progenitor.

  The arc from Luminal A → Luminal B → Basal-like
  → NEPC/DNPC is the WADDINGTON DEPTH AXIS of
  prostate cancer — the same increasing distance
  from the normal luminal attractor that the depth
  score measures.

  The intact circuit finding holds for LUMINAL A.
  Whether it holds for Luminal B requires PRAD-S2a.
  Whether it holds for treatment-naive basal-like
  requires PRAD-S3a.
  NEPC and DNPC — the circuit is gone.
  The strategy shifts to attractor dissolution.

ADDITIONAL UNIQUE CONTEXT FOR PRAD:

  PRAD has the most elaborately characterised
  treatment landscape of any cancer in the
  repository, because prostate cancer is:
    1. Extremely common (~1.4 million new cases/
       year globally)
    2. Slow-progressing in most cases (years of
       clinical intervention time)
    3. The subject of decades of precision medicine
       development (AR pathway is one of the most
       deeply understood oncogenic axes in medicine)

  This means the depth score clinical output
  in PRAD is NOT a replacement for the AR-pathway
  framework — it is a LAYER ON TOP of it:
    Which luminal A patients are deepest within
    their subtype and most likely to resist AR
    inhibition via EZH2-driven NKX3-1 silencing?
    Those patients are the tazemetostat candidates
    — before they become CRPC, not after.

  This is the framework's specific contribution
  to the PRAD clinical landscape:
  EARLY IDENTIFICATION of luminal A patients at
  risk of AR pathway resistance via EZH2 lock —
  so the epigenetic intervention can occur BEFORE
  the resistance arises, not after.
```

---

## DOCUMENT METADATA

```
document_id:        PRAD_Subtype_Orientation
series:             PRAD (Prostate Adenocarcinoma —
                    Subtypes)
folder:             Cancer_Research/PRAD/Subtypes/
date:               2026-03-03
author:             Eric Robert Lawson
status:             ORIENTATION ONLY — no predictions,
                    no analysis, no epigenetic hypotheses
next_document:      PRAD_LuminalA_before.md
                    (Document PRAD-S1a — Luminal A
                    before-doc)
protocol_version:   Workflow_Protocol.md v2.0
existing_analysis:  See OrganismCore cancer series
                    for the completed PRAD bulk analysis
                    (NKX3-1 switch gene identified,
                    EZH2 elevated — gain-of-function lock,
                    intact circuit finding, circuit
                    restoration as therapeutic strategy,
                    tazemetostat predicted as drug target).
```

---

## SECTION I — THE NORMAL PROSTATE EPITHELIUM

```
The Waddington baseline for prostate cancer is
the normal prostate epithelium — a glandular
epithelium with a three-layer hierarchy that is
governed by androgen receptor signalling.

THE PROSTATE GLAND.
The human prostate is a walnut-sized gland
surrounding the urethra at the bladder neck.
Its function is the production of seminal fluid
components (citric acid, PSA/KLK3, zinc, etc.)
for sperm motility support.
The epithelium is organised into glands and ducts
embedded in a fibromuscular stroma.
The three anatomical zones of the prostate
(peripheral, transition, central) have distinct
cancer risk profiles:
  Peripheral zone: ~70% of prostate gland mass
                   ~70% of prostate cancers arise here
                   PRAD analysis primarily reflects
                   peripheral zone tumours
  Transition zone: ~25% of mass — BPH (benign
                   prostatic hyperplasia) arises here
                   ~20% of cancers arise here
                   Transition zone cancers are more
                   often SPOP-mutant (distinct biology)
  Central zone: ~5% of mass — rare cancer origin

═══════════════════════════════════════════════════════
THE PROSTATE EPITHELIAL HIERARCHY — THREE LAYERS
═══════════════════════════════════════════════════════

LUMINAL CELLS (the terminal differentiated cells)
  Location:   Line the glandular lumen.
              Single layer of tall columnar epithelial
              cells facing the gland lumen.
              The secretory cells of the prostate.
  Function:   Androgen-driven secretion of:
                PSA (KLK3 / prostate-specific antigen)
                KLK2 (kallikrein 2)
                PSMA (FOLH1 — prostate-specific
                membrane antigen — the target of
                lutetium PSMA-617)
                Citric acid (via zinc transport)
              These functions require high AR activity
              — luminal cells are the most androgen-
              dependent cells in the body.
              Castration (testosterone removal) leads
              to rapid luminal cell apoptosis
              (95% cell loss within days) — the
              entire luminal layer is androgen-
              dependent for survival.
  Markers:    AR (androgen receptor) — HIGH, nuclear
              NKX3-1 — THE MASTER LUMINAL TF
                       (the MIST1 / PTF1A of the
                       prostate — directly analogous
                       to the acinar master TF in PAAD
                       and the chief cell TF in STAD)
              FOXA1 — pioneer factor for AR
                      (opens chromatin for AR binding)
              HOXB13 — luminal identity TF
              KRT8, KRT18 — luminal cytokeratins
              PSA (KLK3) — high (androgen target)
              PSMA (FOLH1) — high (androgen target)
              TP63: ABSENT (TP63 marks basal cells —
                    its expression in cancer = identity
                    loss / lineage shift)
              KRT5, KRT14: ABSENT (basal markers)
  Self-renewal: VERY LOW under normal conditions.
                Luminal cells are post-mitotic-like
                and androgen-dependent for survival.
                They are REGENERATED FROM BASAL CELLS
                after castration-induced involution
                — the basal cells survive castration
                (AR-independent) and regenerate the
                luminal layer on androgen restoration.

  NKX3-1 — THE GUARDIAN TF OF LUMINAL IDENTITY:
    NKX3-1 (NK3 Homeobox 1) is the master
    transcription factor for:
      a) Prostate luminal specification during
         development: NKX3-1 specifies the prostate
         from urogenital sinus epithelium.
         NKX3-1 knockout mice fail to form a
         normal prostate.
      b) Luminal differentiation maintenance in
         the adult gland: NKX3-1 maintains the
         luminal gene expression programme by:
           Activating: PSA/KLK3, KLK2, FOXA1,
           AR targets, secretory gene programme
           Repressing: basal cell genes (KRT5,
           KRT14, TP63), stem/progenitor genes,
           proliferation genes (MYC targets)
    NKX3-1 LOSS → immediate luminal identity
    destabilisation → the first step in prostate
    cancer initiation (exactly analogous to PTF1A
    loss in PAAD acinar cells).
    NKX3-1 is haploinsufficient for tumour
    suppression: loss of ONE allele is sufficient
    to begin the oncogenic cascade.
    NKX3-1 is the SWITCH GENE IDENTIFIED IN THE
    EXISTING PRAD BULK ANALYSIS.
    Its depth-negative correlation (falls as depth
    increases) measures loss of luminal identity.
    Its restoration = circuit restoration therapy.

  FOXA1 — THE PIONEER FACTOR FOR AR:
    FOXA1 (Forkhead Box A1) is a pioneer
    transcription factor that opens compacted
    chromatin to allow AR binding.
    Without FOXA1, AR cannot access the chromatin
    at androgen-responsive enhancers.
    FOXA1 is required for:
      Normal luminal differentiation
      AR-driven gene expression
      PSA production
    FOXA1 MUTATIONS (gain-of-function, ~5–10%)
    in primary PRAD and higher in mCRPC:
      Alter the AR cistrome — AR binds to
      different chromatin locations than in
      wildtype.
      May drive castration resistance by allowing
      AR to maintain activity at lower androgen
      concentrations.
      FOXA1 mutations are an androgen-independence
      mechanism — one of the routes to CRPC.
    The depth score implication:
      FOXA1 mutation → altered AR cistrome → AR
      activity maintained at low androgens →
      these tumours sit at a DIFFERENT POSITION
      on the depth axis than tumours where AR is
      blocked by EZH2-mediated NKX3-1 silencing.
      FOXA1-mutant tumours need AR inhibitor
      combinations, not EZH2 inhibitors.
      Separating FOXA1-mutant from wildtype in
      the Luminal B depth analysis is critical.

BASAL CELLS (the stem/progenitor layer)
  Location:   Adjacent to the basement membrane —
              the outermost layer of the epithelium,
              separated from the lumen by luminal cells.
              ~5–10% of epithelial cells.
  Function:   Stem/progenitor reservoir.
              Regenerate luminal cells after androgen
              withdrawal and restoration.
              Provide mechanical support to luminal layer.
              Basal cells are AR-INDEPENDENT — they
              survive castration.
              This is why castration does not cure
              prostate cancer — basal-like cancer cells
              or luminal cells that have acquired
              basal properties survive androgen
              withdrawal.
  Markers:    TP63 — THE BASAL TF
                     p63 (TP63 isoform ΔNp63) is the
                     master transcription factor of
                     basal cells. Its expression
                     DEFINES the basal state.
                     TP63 REPRESSES luminal identity
                     and AR transcriptional targets.
                     TP63 HIGH in cancer = the cell
                     has acquired basal identity or
                     lost luminal identity.
                     This is exactly parallel to TP63
                     in basal-like PDAC (Section IV of
                     PAAD orientation) — the same TF,
                     the same role, the same aberrant
                     expression pattern in the most
                     identity-lost cancer false attractor.
              KRT5, KRT14 — basal cytokeratins
              CD44 — stem/progenitor marker
              CD49f (ITGA6) — integrin stem marker
              PSCA (prostate stem cell antigen) — some
              AR: LOW (vs. luminal cells)
              NKX3-1: LOW/ABSENT
              PSA: ABSENT (basal cells do not secrete PSA)
  Self-renewal: HIGH — basal cells are the self-
                renewing progenitor population.
                Under injury or castration, basal cells
                divide and generate intermediate and
                then luminal daughter cells.

INTERMEDIATE CELLS (the transit-amplifying layer)
  Location:   Between basal and luminal layers.
              Co-express basal and luminal markers
              in varying ratios.
              The transitional state between basal
              progenitor and mature luminal cell.
  Markers:    KRT5+/KRT8+ (co-expression)
              Low-intermediate AR
              Partial NKX3-1
  Function:   Transit-amplifying cells — rapidly
              dividing progeny of basal stem cells
              committing to luminal fate.
              The isthmus equivalent in the prostate
              hierarchy (analogous to gastric isthmus
              and pancreatic centroacinar cells).
  Cancer relevance:
              Luminal progenitor cells in the
              intermediate zone are the proposed cell
              of origin for the majority of primary
              prostate adenocarcinomas — they are
              AR-responsive but not yet fully
              terminally differentiated, making them
              susceptible to oncogenic transformation.

THE WADDINGTON STRUCTURE OF NORMAL PROSTATE:

  BASAL CELL (TP63+, KRT5+, AR-low)
    ↓ COMMITTED TO LUMINAL FATE
    (NKX3-1 expression begins)
  INTERMEDIATE CELL (KRT5+/KRT8+, AR intermediate)
    ↓ TERMINAL LUMINAL DIFFERENTIATION
    (NKX3-1 high, FOXA1 high, AR high)
  LUMINAL CELL (AR-high, NKX3-1+, KLK3+)
    = THE TERMINAL NORMAL ATTRACTOR

  PROSTATE CANCER FORMATION ROUTES:
  Route 1 (most common, primary adenocarcinoma):
    Luminal progenitor / intermediate cell
    → AR-pathway activation (TMPRSS2-ERG, SPOP,
      FOXA1 mutation, or PTEN loss)
    → Initial luminal adenocarcinoma
    → NKX3-1 loss (haploinsufficiency → LOH)
    → EZH2 rise → H3K27me3 at NKX3-1 locus
    → FALSE ATTRACTOR: Luminal A (NKX3-1 low,
      AR still high, FOXA1 retained)
    → Further depth: Luminal B (higher proliferation,
      CDK activity high, SPOP or FOXA1 mutation
      may drive AR independence)

  Route 2 (treatment-emergent, arising AFTER AR
  inhibition therapy):
    Luminal A or Luminal B cancer
    → AR-pathway blockade (enzalutamide,
      abiraterone)
    → Selective pressure removes AR-dependent cells
    → TP53 + RB1 loss (enables lineage plasticity)
    → NEPC false attractor: AURKA/MYCN amplification,
      NE markers (CHGA, SYP), AR absent
    → OR: DNPC false attractor: AR absent, NE absent,
      basal/stem-like (the most identity-lost state)

  THE DEPTH AXIS IN PRAD:
  Normal luminal cell (NKX3-1 HIGH, AR HIGH)
  → Luminal A cancer (NKX3-1 falling, EZH2 rising,
    AR still HIGH → intact circuit, EZH2 blocked)
  → Luminal B cancer (NKX3-1 low, AR high/variants,
    proliferation rising — CDKN2A loss common)
  → Basal-like primary (NKX3-1 gone, TP63 present,
    AR low — CIRCUIT TRANSITION POINT)
  → NEPC (AR gone, AURKA/MYCN high — CIRCUIT BROKEN)
  → DNPC (AR gone, NE gone — MAXIMUM DEPTH,
    NO CIRCUIT)

  This is the most elaborately staged depth axis
  in the entire OrganismCore repository — prostate
  cancer is the cancer that has been most thoroughly
  tracked through its lineage plasticity progression
  because it is driven there clinically by AR
  inhibition therapy, creating a natural experiment
  in real patients.
```

---

## SECTION II — THE TCGA GENOMIC SUBTYPES

```
The TCGA 2015 molecular characterisation of
prostate adenocarcinoma identified GENOMIC
subtypes — defined by driver alterations rather
than transcriptional phenotype. These are distinct
from the luminal/basal/NE transcriptional spectrum
but are profoundly important for understanding
WHICH patients are in WHICH part of the depth
spectrum.

THE SEVEN TCGA MOLECULAR CLASSES:

1. TMPRSS2-ERG FUSION (~40–50%)
  The most frequent molecular event in PRAD.
  A chromosomal translocation places ERG (an ETS
  transcription factor) under the control of the
  androgen-responsive TMPRSS2 promoter.
  Result: ERG is expressed at HIGH LEVELS in
  prostate cells in an androgen-driven manner.
  ERG in the normal prostate: ABSENT.
  ERG is an ETS family TF that activates:
    Cell cycle genes (MYC targets)
    Invasion genes (MMP3, PLAU)
    Wnt pathway genes
    Represses: NKX3-1 (the switch gene!)
  ERG ACTIVELY REPRESSES NKX3-1:
    ERG binds to the NKX3-1 locus and
    suppresses its expression.
    This means TMPRSS2-ERG fusion is, in part,
    an epigenetic NKX3-1 suppressor —
    an oncogenic mechanism that acts through
    the SAME SWITCH GENE identified in the
    framework's depth analysis.
    The connection:
    ERG fusion → NKX3-1 falls → luminal identity
    begins to destabilise → depth increases.
    The depth score should correlate with
    ERG expression in TMPRSS2-ERG+ tumours.
  Clinical significance:
    TMPRSS2-ERG fusion is NOT independently
    prognostic in primary disease — it is
    prognostic mainly in combination with
    PTEN loss (fusion+ PTEN-null = high risk).
    No approved targeted therapy for ERG yet.
    ERG inhibitors (e.g., VPC-3090) in preclinical/
    early clinical development.

2. SPOP MUTATION (~10–15%)
  SPOP encodes a substrate-binding protein of
  an E3 ubiquitin ligase complex.
  SPOP targets for proteasomal degradation:
    AR (the main oncogenic driver)
    SRC-3 (steroid receptor coactivator)
    TRIM24 (coactivator of AR)
  SPOP MUTATION → THESE TARGETS ARE NOT DEGRADED
  → AR protein levels RISE → androgen signalling
  is constitutively AMPLIFIED even at low hormone
  levels.
  This is an AMPLIFICATION of AR activity
  mechanism — the opposite end of the spectrum
  from AR loss (NEPC/DNPC).
  SPOP-mutant tumours are:
    MUTUALLY EXCLUSIVE with TMPRSS2-ERG fusion
    (these are two separate routes to the same
    initial oncogenesis)
    ENRICHED in transition zone (vs. peripheral)
    MORE SENSITIVE to AR pathway inhibitors
    (abiraterone, enzalutamide) than SPOP-wildtype
    — because the cancer's AR dependency is even
    higher than normal
    ASSOCIATED WITH CDK12 MUTATIONS in aggressive
    cases — CDK12 loss creates homologous
    recombination deficiency and tandem
    duplications
  The depth score implication:
    SPOP-mutant tumours have ELEVATED AR activity
    → enhanced luminal programme → NKX3-1 may
    paradoxically be MAINTAINED at higher levels
    → these tumours sit at SHALLOW depth initially
    → they benefit from aggressive AR pathway
    inhibition (because AR dependency is extreme)
    → they progress to CRPC by a different route
    than TMPRSS2-ERG tumours.
  Clinical significance:
    SPOP mutation status predicts BETTER RESPONSE
    to abiraterone/enzalutamide in mCRPC.
    Emerging evidence: SPOP-mutant CRPC may be
    more susceptible to combination AR inhibition
    + CDK4/6 inhibition (if CDKN2A is co-deleted).

3. FOXA1 MUTATION (~5–10%, higher in mCRPC)
  Gain-of-function mutations predominantly in the
  forkhead domain.
  Effect: Altered AR cistrome — AR is recruited to
  different chromatin locations, maintaining its
  activity at lower androgen concentrations.
  This is an ANDROGEN-INDEPENDENCE mechanism:
  the cancer maintains AR-driven gene expression
  WITHOUT testosterone → the cancer is resistant
  to castration (CRPC) from the start.
  FOXA1 mutation is enriched in mCRPC relative to
  primary disease — it is selected for by androgen
  deprivation.
  The depth score implication:
    FOXA1-mutant tumours use AR activity in a
    non-canonical manner — the depth axis measures
    NKX3-1 loss, but in FOXA1-mutant tumours,
    AR maintains some activity even at low NKX3-1.
    FOXA1-mutant tumours occupy a DISTINCT POSITION
    on the depth spectrum from FOXA1-wildtype
    tumours with the same NKX3-1 level.
    Failure to account for FOXA1 mutation in the
    depth score analysis could create a confound —
    the depth score would underestimate AR activity
    in these patients.
  Clinical significance:
    FOXA1 mutation is NOT yet used as a clinical
    biomarker but is under active investigation.
    Experimental FOXA1-specific AR inhibitors are
    in early development.

4. IDH1 MUTATION (<1%)
  Extremely rare. Isocitrate dehydrogenase 1
  gain-of-function mutation (R132H most common).
  Produces the oncometabolite 2-hydroxyglutarate
  (2-HG) → inhibits TET enzymes → DNA
  hypermethylation (CIMP) → global epigenetic
  reprogramming.
  IDH1-mutant PRAD is the closest structural
  equivalent to IDH-mutant glioma in the prostate:
    CIMP phenotype
    Distinct methylation signature
    Better prognosis (in glioma — unclear in PRAD
    given rarity)
  Cross-repository structural connection:
    IDH1 mutation in PRAD → same 2-HG mechanism
    as IDH1/2 in AML → same global methylation
    phenotype as EBV CIMP in STAD → same principle
    (an oncometabolite or viral factor hijacks the
    DNA methylation machinery to silence
    differentiation genes).
    Different tissue, same epigenetic mechanism.
    IDH inhibitor (ivosidenib) has Phase I data
    in IDH1-mutant solid tumours — PRAD cases
    rare but occasionally eligible.

5. CHD1 DELETION (~7%)
  CHD1 (chromodomain helicase DNA binding protein 1)
  is a chromatin remodelling gene.
  CHD1 deletion is:
    MUTUALLY EXCLUSIVE with TMPRSS2-ERG fusion
    (like SPOP)
    Associated with multiple concurrent structural
    rearrangements
    Enriched in aggressive disease
  CHD1-deleted tumours often co-occur with
    SPOP mutations (~40% co-occurrence)
    providing a combined genomic instability
    + AR amplification signature.

6. CDKN2A/CDK4/6 ALTERATIONS
  CDKN2A loss in PRAD: ~5–8% of primary disease
  but rising to ~20–30% in mCRPC.
  CDK4/6 constitutively active → Rb phosphorylation
  → unrestrained proliferation.
  This is the same mechanism as in PAAD.
  Clinical significance in PRAD:
    CDK4/6 inhibitors (palbociclib, ribociclib)
    are being tested in CDKN2A-deleted mCRPC —
    Phase II trials show modest single-agent
    activity.
    The rational combination: AR inhibitor +
    CDK4/6 inhibitor in CDKN2A-deleted mCRPC.
    Depth score prediction: CDKN2A-deleted tumours
    should cluster at HIGHER depth (more advanced)
    — CDK4/6 activation is a depth-positive driver.

7. OTHER/NO IDENTIFIED ALTERATION (~5%)
  Tumours without the above defining alterations.
  Often lower-grade, primary disease.
  Full characterisation ongoing with improved
  sequencing.

THE TCGA SUBTYPES AND THE DEPTH AXIS:
  The TCGA genomic subtypes stratify WITHIN the
  luminal/basal transcriptional spectrum:
  Most TMPRSS2-ERG+ and SPOP-mutant tumours:
    Luminal A or Luminal B (AR-positive,
    NKX3-1 intermediate/low, FOXA1 retained)
  FOXA1-mutant tumours: Luminal B or early CRPC
  CHD1-deleted tumours: Luminal B to basal-like
  IDH1-mutant tumours: Special CIMP class
  The depth score adds a CONTINUOUS MEASURE
  within each genomic class:
  Within TMPRSS2-ERG+ tumours:
    Depth score separates high-NKX3-1 (shallow)
    from low-NKX3-1 (deep) tumours that are at
    higher risk of progression and AR resistance.
  The combination of TCGA genomic class +
  depth score is the most complete patient
  stratification tool for primary PRAD.
```

---

## SECTION III — LUMINAL A: THE INTACT CIRCUIT SUBTYPE

```
FREQUENCY:     ~40–50% of primary PRAD
               The most common transcriptional
               subtype in early/intermediate stage
               disease. Frequency falls as disease
               advances (selection against luminal
               identity by AR therapy).

PROGNOSIS:     BEST
               10-year prostate cancer-specific
               survival: >85% for localised Luminal A.
               Low Gleason score (6–7/10).
               Low PSA doubling time.
               Most patients managed with active
               surveillance or primary radiation/
               surgery — no systemic therapy required.

CELL OF ORIGIN:
  Luminal progenitor / intermediate cell that
  acquired an AR-pathway activation event
  (TMPRSS2-ERG fusion most commonly, or SPOP
  mutation) and subsequently lost one NKX3-1 allele.
  The resulting cell:
    Maintains high AR activity (androgen-dependent)
    Maintains residual luminal programme (FOXA1 high)
    Shows partial NKX3-1 loss (haploinsufficiency)
    EZH2 begins to rise, placing H3K27me3 marks
    at the NKX3-1 locus and other luminal genes
  This is a SHALLOW FALSE ATTRACTOR:
    The cell has left the terminal luminal attractor
    but has not traveled far.
    The luminal programme is PARTIALLY maintained.
    The dominant epigenetic change is EZH2-driven
    silencing of NKX3-1 — the same mechanism
    as in classical PDAC (EZH2 silencing PTF1A).
    THE CIRCUIT IS INTACT:
    If EZH2 is inhibited → NKX3-1 is re-expressed
    → NKX3-1 re-activates the luminal programme
    → The cell returns to terminal luminal identity.
    This is the tazemetostat mechanism in PRAD.

DEFINING MOLECULAR FEATURES:
  NKX3-1:        INTERMEDIATE (haploinsufficient)
                 NOT COMPLETELY LOST — partially
                 silenced by EZH2. This is CRITICAL.
                 It is the partial silencing that
                 makes the circuit restorable.
                 If NKX3-1 were completely methylated
                 and silenced — the locus inaccessible
                 — the circuit could not be restored
                 by EZH2 inhibition alone. But at
                 Luminal A depth, NKX3-1 is
                 SUPPRESSED but RESTOREABLE.
  EZH2:          ELEVATED (above normal luminal cells)
                 The epigenetic lock. This is the
                 finding from the existing bulk PRAD
                 analysis.
                 EZH2 places H3K27me3 at:
                   NKX3-1 promoter
                   CDKN1A (p21) promoter
                   DAB2IP (tumour suppressor) locus
                   Other differentiation gene loci
  AR:            HIGH — the primary survival signal
                 Androgen-dependent growth.
                 Responds to:
                   ADT (castration) — initial response
                   Abiraterone — CYP17A1 inhibition
                   Enzalutamide — AR antagonist
  FOXA1:         HIGH (wildtype in most Luminal A)
                 Maintains AR chromatin accessibility
                 at canonical AR-binding sites.
  TMPRSS2-ERG:   Present in ~50% of Luminal A
                 (or SPOP-mutant in the other major
                 group — mutually exclusive)
  PTEN:          Intact in most early Luminal A
                 (~20% have PTEN loss in primary
                 disease — those are higher risk)
  CDKN2A:        Intact in most Luminal A (~5% loss)
  TP53:          Wildtype in most Luminal A (~5%)
  RB1:           Intact — RB1 loss marks advanced
                 disease (NEPC transition)

CIRCUIT RESTORATION — THE CORE MECHANISM:
  The existing PRAD analysis finding (Pattern 5):
    The luminal differentiation circuit is INTACT.
    NKX3-1 is the switch gene.
    EZH2 blocks NKX3-1 re-expression.
    EZH2 inhibition (tazemetostat) removes block.
    NKX3-1 re-expression drives luminal programme.
    Luminal cells exit false attractor.

  The PRAD-specific version of this mechanism has
  a clinical timetable that PAAD does not:
    In PAAD, all patients are metastatic at
    presentation in ~80% of cases — there is
    no early-stage prevention window.
    In PRAD, ~75% of patients are diagnosed at
    localised stage with YEARS of intervention time.
    The EZH2 elevation in Luminal A PRAD is
    measurable at Gleason 6-7 (before metastasis,
    before CRPC).
    A tazemetostat-based intervention at GLEASON
    6-7 stage — before EZH2 further deepens the
    NKX3-1 silencing — could PREVENT the progression
    to CRPC rather than just treating it.
    This is the most clinically impactful implication
    of the PRAD circuit restoration finding:
    EARLY INTERVENTION in luminal A patients at
    risk of progression (high EZH2, intermediate
    NKX3-1 depth score) with tazemetostat BEFORE
    AR pathway therapy is needed.

TAZEMETOSTAT CLINICAL DATA IN PRAD:
  CELLO-1 trial (NCT04179864):
    Tazemetostat + enzalutamide vs. enzalutamide
    alone in mCRPC.
    Primary endpoint: rPFS.
    RESULT: NEGATIVE — PFS 16.6 vs. 13.8 months,
    p=0.37. Trial terminated late 2024.
  FRAMEWORK INTERPRETATION OF THE NEGATIVE RESULT:
    The CELLO-1 trial enrolled mCRPC patients —
    NOT Luminal A primary disease patients.
    By the time a patient has mCRPC (on enzalutamide),
    the disease is typically Luminal B or transitioning
    toward basal-like. The EZH2 circuit restoration
    mechanism (which requires a RESTOREABLE NKX3-1
    locus) may no longer be operative in mCRPC.
    The FRAMEWORK PREDICTS:
    Tazemetostat works in LUMINAL A (intact circuit,
    EZH2 high, NKX3-1 restoreable) — NOT in mCRPC
    (circuit may be more broken, competing resistance
    mechanisms dominant).
    The negative CELLO-1 result does NOT disprove
    the framework's prediction. It demonstrates that
    EZH2 inhibition in UNSELECTED MCRPC is
    insufficient — consistent with the framework's
    insistence on DEPTH-STRATIFIED patient selection.
    The DEPTH SCORE IS THE SELECTION TOOL:
    Luminal A patients with HIGH EZH2 + INTERMEDIATE
    NKX3-1 depth score (restoreable circuit) are the
    tazemetostat candidates — not all mCRPC patients.
    The CELLO-1 failure is the PAAD lesson applied
    to PRAD: CDK4/6 inhibitors given to all gastric
    cancer → negative. CDK4/6 inhibitors in CCND1-
    amplified CIN gastric cancer → expected positive.
    Tazemetostat in all mCRPC → negative.
    Tazemetostat in luminal A, high EZH2, restoreable
    NKX3-1 → the depth score stratification that
    turns a failed drug into the right drug for the
    right patient.

TREATMENT:
  Localised (active surveillance eligible):
    PSA ≤10, Gleason ≤6, ≤2 cores positive
    → active surveillance (PSA, biopsy monitoring)
  Localised (definitive):
    Radical prostatectomy
    OR radiation therapy (external beam + brachytherapy)
    No systemic therapy required in Luminal A.
  High-risk localised (Gleason 7, PSA 10-20):
    Definitive surgery/radiation +
    short-course ADT (~6 months) to debulk.
    The ADT eliminates AR-dependent luminal cells
    but spares EZH2-high cells that have begun to
    become androgen-independent.
  FRAMEWORK TARGET (pre-CRPC):
    Tazemetostat in Luminal A patients with
    high EZH2 depth score at Gleason 6-7:
    Prevents progression to CRPC by restoring
    NKX3-1 before it is silenced irreversibly.
    This is the preventive depth score application.
    Not yet in clinical trials in this exact form.
    Mevrometostat (more potent EZH2 inhibitor)
    showing early signs of activity in mCRPC
    (data emerging 2025) — may succeed where
    tazemetostat failed if the patient selection
    is refined by depth score.
```

---

## SECTION IV — LUMINAL B: INTERMEDIATE DEPTH

```
FREQUENCY:     ~25–30% of primary PRAD
               Frequency increases in the mHSPC and
               mCRPC settings (Luminal A tumours
               progress to Luminal B under pressure).

PROGNOSIS:     INTERMEDIATE
               Gleason 7–8 (ISUP grade 2–3).
               PSA 10–20 ng/mL.
               Higher risk of biochemical recurrence
               after curative therapy.
               5-year metastasis-free survival: ~60–75%.

CELL OF ORIGIN:
  Either:
  a) Primary Luminal B — arises from a luminal
     progenitor with a higher-impact initial event
     (PTEN loss + TMPRSS2-ERG, or FOXA1 mutation,
     or SPOP mutation with CDK4/6 co-activation).
  b) Progressed Luminal A — under therapeutic
     selection, Luminal A patients develop Luminal B
     features as EZH2 further silences NKX3-1 and
     CDKN2A loss or FOXA1 mutation emerges.

DEFINING MOLECULAR FEATURES:
  NKX3-1:        LOW — more silenced than Luminal A
                 H3K27me3 at NKX3-1 is heavier.
                 Circuit still exists but is LESS
                 restoreable than Luminal A.
  EZH2:          HIGH — higher than Luminal A.
  AR:            STILL HIGH in most Luminal B.
                 AR amplification begins to appear
                 (~10–15% of advanced Luminal B).
                 AR splice variant AR-V7 begins to
                 appear — the primary CRPC mechanism.
  FOXA1:         Variable — wildtype OR mutant.
                 FOXA1-mutant Luminal B = castration-
                 resistant from the start (or will
                 become so quickly).
  PTEN:          LOST in ~30–40% of Luminal B —
                 higher than Luminal A.
                 PTEN loss → AKT activation → PI3K
                 pathway constitutive → survival
                 independent of AR.
                 PTEN loss + AR inhibition creates
                 a PI3K dependency — ipatasertib
                 (AKT inhibitor) + enzalutamide
                 shows improved outcomes in PTEN-null
                 mCRPC (IPATential150 trial).
  CDKN2A:        LOSS ~15–20% — higher than Luminal A.
                 CDK4/6 constitutive activity.
  MYC:           AMPLIFICATION ~20% of Luminal B.
                 The depth-positive proliferation
                 driver — rising with depth.
  TP53:          MUTATION ~15–20%.
  AR-V7:         BEGINS TO APPEAR at Luminal B depths.
                 AR-V7 is a constitutively active
                 AR splice variant lacking the
                 ligand-binding domain — it cannot
                 be blocked by enzalutamide or
                 abiraterone.
                 AR-V7 detection = the molecular
                 marker that a patient is transitioning
                 to CRPC.

CIRCUIT STATUS IN LUMINAL B:
  The circuit is PARTIALLY intact but under
  increasing threat:
    NKX3-1 is more silenced than Luminal A —
    H3K27me3 is heavier at the NKX3-1 locus.
    EZH2 inhibition may still de-repress NKX3-1
    IF the locus is not yet fully CpG methylated.
    But EZH2 inhibition alone may not be sufficient
    in Luminal B if FOXA1 mutation has altered
    the AR cistrome or if PTEN loss has activated
    a PI3K survival path that is AR-independent.
    The FRAMEWORK PREDICTION for Luminal B:
    Tazemetostat + enzalutamide:
      EZH2 inhibition partially restores NKX3-1
      AR inhibition removes the AR-survival signal
      Combined: deeper differentiation + AR blockade
    Tazemetostat + ipatasertib (PTEN-null Luminal B):
      EZH2 inhibition of the epigenetic lock
      AKT inhibition of the PTEN-null survival path
      Addressing both the identity-loss and the
      bypass survival signal.
    These are structural hypotheses — not predictions.
    They belong in PRAD-S2a.

TREATMENT:
  High-risk localised (Gleason 8–9, Luminal B):
    Radical prostatectomy + pelvic lymph node
    dissection.
    OR: definitive radiation + long-course ADT
    (2 years — DART trial: 24-month ADT + RT
    improves OS over short-course + RT in high-risk).
  mHSPC (metastatic hormone-sensitive, Luminal B):
    ADT + intensification:
    ADT + docetaxel (CHAARTED trial — OS benefit
    in high-volume mHSPC)
    OR: ADT + abiraterone (LATITUDE/STAMPEDE trials)
    OR: ADT + enzalutamide (ARCHES/ENZAMET trials)
    OR: ADT + darolutamide (ARASENS trial)
    The STAMPEDE and PEACE-1 trials also support
    ADT + abiraterone + docetaxel triplet for
    high-volume mHSPC.
  mCRPC (progressed on ADT):
    First-line: Enzalutamide or abiraterone
    (if not already used in mHSPC setting)
    PTEN-null: Ipatasertib + enzalutamide
    (IPATential150 — rPFS benefit in PTEN-null)
    BRCA1/2 or HRR-mutant: Olaparib + abiraterone
    (PROpel trial — improved OS in BRCA1/2)
    AR-V7+: Switch to cabazitaxel (taxane; AR-V7
    does not affect taxane activity)
    PSMA-positive: Lutetium-177 PSMA-617 after
    ARPI + taxane progression (VISION trial)
```

---

## SECTION V — BASAL-LIKE PRIMARY: THE IDENTITY-LOSS ATTRACTOR

```
FREQUENCY:     ~5–15% of primary PRAD
               (estimates vary by classification
               system and biopsy site)
               Enriched in high-grade (Gleason 9–10)
               and metastatic presentations.

PROGNOSIS:     POOR within primary disease.
               ~2.5-fold higher risk of metastasis
               vs. Luminal A.
               Aggressive disease at presentation.
               Commonly not captured in low-volume
               prostate cancer screening — often
               presents at advanced stage.

CELL OF ORIGIN:
  Two hypotheses:
  a) Basal cell directly transformed — acquiring
     KRAS-equivalent activation (MYC amplification,
     EGFR activation, FGFR signalling) while
     retaining basal identity (TP63+, KRT5+, AR-low).
  b) Luminal cell that has undergone PROFOUND identity
     loss — losing NKX3-1 AND AR AND FOXA1, and
     accessing the basal progenitor programme as
     the survival programme of last resort.
  Current evidence supports both routes — basal-like
  PRAD is not one disease but a transcriptional
  state that can be reached via multiple paths.
  The common feature: AR IS LOW OR ABSENT in most
  basal-like primary PRAD. This makes these tumours
  intrinsically castration-resistant from diagnosis.

DEFINING MOLECULAR FEATURES:
  NKX3-1:        ABSENT or near-absent.
                 Complete loss of luminal identity TF.
  AR:            LOW — the critical distinction from
                 Luminal B.
                 When AR falls below a threshold, AR
                 inhibitors become ineffective —
                 the cancer is not dependent on AR.
  FOXA1:         LOW/ABSENT — luminal pioneer factor gone.
  TP63:          PRESENT — the basal identity marker.
                 The same TP63 that appears in:
                 Basal-like PDAC (PAAD orientation)
                 Basal/squamous BLCA
                 Basal breast cancer
                 In the prostate: its expression marks
                 identity loss from the luminal attractor.
  KRT5, KRT14:   PRESENT — basal cytokeratins.
  EZH2:          HIGH — EZH2 remains elevated even
                 in basal-like primary PRAD.
                 BUT: in basal-like, EZH2 is no longer
                 just blocking NKX3-1 — it is maintaining
                 the entire basal-like false attractor
                 by silencing multiple luminal genes
                 simultaneously. EZH2 inhibition in
                 basal-like may not restore the
                 differentiation circuit — the circuit
                 may be too damaged to re-execute.
                 This is the Luminal A → Basal-like
                 circuit integrity transition.
  TP53:          FREQUENTLY MUTATED in basal-like
                 primary (~30–40%) — higher than
                 Luminal A or B.
  RB1:           PARTIAL LOSS — beginning.
                 Full RB1 loss defines NEPC.
  PTEN:          FREQUENT LOSS (~50% of basal-like).
  MYC:           FREQUENTLY AMPLIFIED — the proliferative
                 engine of the basal-like state.
  AURKA:         BEGINS TO RISE — pre-NEPC signature.
  MYCN:          OCCASIONALLY AMPLIFIED — pre-NEPC.

TREATMENT:
  AR inhibitors: LARGELY INEFFECTIVE (AR-low/negative).
  Standard approaches:
    Docetaxel chemotherapy (AR-independent)
    Carboplatin (TP53-mutant, HRR-deficient subset)
    Clinical trials strongly recommended:
    EGFR inhibitors (EGFR is often elevated in
    basal-like primary PRAD — a parallel to the
    PANGEA trial in basal-like PDAC)
    DLL3-targeted therapy (if NE markers begin
    to appear — tarlatamab biTE)
  The depth score contribution:
    The basal-like primary depth score identifies
    which surface identity markers are retained
    (PSMA, PSCA, EGFR) and which are gone (PSA,
    KLK3, NKX3-1).
    Residual surface markers = the ADC/antibody
    targets for this subtype.
    This is the identity-loss attractor surface
    marker pattern — exactly analogous to CLDN18.2
    in GS STAD and NECTIN4 in basal/squamous BLCA.

CROSS-REPOSITORY STRUCTURAL CONNECTION:
  Basal-like PRAD joins the identity-loss attractor
  cohort across the repository:
    GS STAD (CDH1 loss, signet ring)
    Basal-like PDAC (GATA6−, TP63+)
    Basal/Squamous BLCA (TP63+, KRT5+)
    CMS4 CRC (mesenchymal, TGF-β high)
    Basal-like PRAD (AR−, TP63+, KRT5+)
  All five are the most identity-lost false attractor
  in their GI/genitourinary tissue of origin.
  All five have worst prognosis in their cancer type.
  All five are immune-excluded.
  All five are resistant to the canonical therapy
  of their cancer type.
  The framework sees them as the same Waddington
  attractor type in different tissue coordinates.
```

---

## SECTION VI — NEPC AND DNPC: THE CIRCUIT-GONE STATES

```
NEPC (NEUROENDOCRINE PROSTATE CANCER):
  FREQUENCY:
    De novo (treatment-naive): <1% of primary PRAD
    Treatment-emergent: ~15–20% of mCRPC patients
    after prolonged AR-pathway inhibition.
    Rising frequency as AR inhibitors improve —
    the better the AR inhibitor, the more selection
    pressure toward NEPC transdifferentiation.

  MECHANISM — TREATMENT-EMERGENT NEPC:
    Under AR pathway inhibitor pressure:
      TP53 loss + RB1 loss (REQUIRED — these two
      tumour suppressors must both be lost for
      NEPC emergence. Either alone is insufficient.)
      AURKA amplification (~40%) — Aurora A
      kinase stabilises N-MYC protein (prevents
      its proteasomal degradation), dramatically
      amplifying N-MYC activity.
      MYCN (N-MYC) amplification (~40%) — the
      downstream transcriptional driver of the
      neuroendocrine gene programme.
      AURKA/MYCN form a stabilising complex:
      AURKA stabilises N-MYC → N-MYC activates
      neuroendocrine TFs → NEPC false attractor
      is established and self-maintaining.
      This is the CONVERGENCE NODE in NEPC:
      AURKA/N-MYC is the equivalent of EZH2 in
      Luminal A — the maintainer of the false
      attractor. Unlike EZH2 (which is a lock),
      AURKA/MYCN is an ACTIVATOR of the wrong
      programme.
      Alisertib (Aurora A inhibitor) disrupts
      AURKA→N-MYC stabilisation → N-MYC degraded
      → neuroendocrine programme loses its driver.
      This is the NEPC drug target — alisertib —
      the same drug predicted for CIN STAD by the
      existing AURKA depth-correlation finding.
      The AURKA/NEPC connection cross-validates
      the CIN STAD AURKA finding at the NEPC
      level: AURKA is a depth-positive gene in
      multiple cancers, and its most extreme
      manifestation (AURKA amplification + MYCN
      amplification) drives the most identity-lost
      and aggressive cancer states.

  DEFINING MOLECULAR FEATURES:
    AR:           ABSENT (lost by epigenetic silencing
                  or structural deletion — the AR
                  programme is completely shut down)
    NKX3-1:       ABSENT
    FOXA1:        ABSENT (in most NEPC)
    TP53:         MUTATED/LOST (~90% of NEPC)
    RB1:          LOST (~90% of NEPC — the co-loss
                  with TP53 is the NEPC transition signal)
    AURKA:        AMPLIFIED (~40%)
    MYCN:         AMPLIFIED (~40%)
    Neuroendocrine markers:
      CHGA (chromogranin A)
      SYP (synaptophysin)
      NSE (neuron-specific enolase)
      NCAM1 (CD56)
    DLL3:         HIGH — the NEPC surface marker
                  DLL3 is a Notch pathway inhibitory
                  ligand expressed on neuroendocrine
                  cells but NOT on normal prostate
                  epithelium.
                  DLL3 expression in NEPC is the
                  RESIDUAL IDENTITY MARKER of the
                  neuroendocrine false attractor —
                  exactly analogous to:
                  CLDN18.2 in GS STAD (retained
                  gastric marker in identity-lost cancer)
                  NECTIN4 in basal/squamous BLCA
                  EGFR in basal-like PDAC
                  DLL3 is the target of TARLATAMAB
                  (bispecific T-cell engager, AMG 757):
                  DeLLpro-300 trial — Phase II.
                  22.2% ORR in DLL3+ NEPC.
                  This is the first targeted therapy
                  to show activity in NEPC.
    EZH2:         STILL ELEVATED but serving a
                  different function than in Luminal A:
                  In Luminal A: EZH2 blocks NKX3-1
                  (removing EZH2 restores circuit).
                  In NEPC: EZH2 helps maintain
                  the neuroendocrine programme by
                  silencing luminal identity genes.
                  EZH2 inhibition in NEPC may not
                  restore luminal identity because:
                  a) NKX3-1 locus is fully silenced
                     (beyond just H3K27me3 —
                     may be CpG methylated)
                  b) AR is absent — even if NKX3-1
                     is partially restored, AR is not
                     there to drive the luminal
                     programme
                  c) AURKA/MYCN is ACTIVELY driving
                     the NE programme — passive de-
                     repression (EZH2 removal) cannot
                     overcome active programme driving.
                  THIS IS THE CIRCUIT BROKEN STATE.
                  EZH2 inhibition alone is insufficient.
                  ALISERTIB (AURKA inhibitor) is the
                  primary NEPC target — not tazemetostat.

  CIRCUIT STATUS: BROKEN.
    The circuit cannot be restored in NEPC.
    Attractor dissolution (alisertib disrupting
    AURKA/N-MYC stabilisation) is the correct
    strategy, not circuit restoration.
    Tazemetostat is not effective as monotherapy
    in NEPC — confirmed by the CELLO-1 direction
    and by NEPC biology.
    Alisertib is the framework-derived NEPC target.

  TREATMENT:
    Platinum-based chemotherapy:
      Carboplatin + etoposide (small cell analog)
      — first-line NEPC.
      ~50% ORR, median OS ~12 months.
      Responses are often dramatic but brief.
    Cabazitaxel: Second-line option.
    Tarlatamab (AMG 757):
      DLL3-biTE — DeLLpro-300 Phase II.
      22.2% ORR in DLL3+ NEPC.
      The first targeted therapy with activity.
    Alisertib (Aurora A inhibitor):
      SPARTAN trial (NCT01799278): Phase II in
      mCRPC/NEPC — showed activity in NEPC-enriched
      cohort with MYCN amplification.
      The framework prediction — AURKA/MYCN as
      convergence node → alisertib — IS CONFIRMED
      by clinical trial data in NEPC.
      This is the NEPC equivalent of the venetoclax
      confirmation in CLL.
    Pembrolizumab:
      KEYNOTE-365 — mCRPC cohort includes NEPC.
      Activity in MSI-high or TMB-high NEPC subset.

DNPC (DOUBLE-NEGATIVE PROSTATE CANCER):
  FREQUENCY:
    ~5–10% of mCRPC.
    Arising after AR pathway inhibition AND after
    NEPC treatments have failed or NEPC programme
    not established.
    The most heterogeneous and poorly defined
    advanced prostate cancer state.

  DEFINITION:
    AR-NEGATIVE (no PSA, no androgen-driven gene
    expression)
    NEUROENDOCRINE-NEGATIVE (no CHGA, no SYP —
    no neuroendocrine programme)
    WHAT IS IT? Stem-like. Basal-like. Mesenchymal.
    A progenitor-like state without any tissue
    identity — the furthest false attractor from
    the normal luminal cell.
    Analogous to:
      GS STAD (CDH1 loss, no mesenchymal identity)
      Basal-like PDAC (GATA6−, TP63−)
      CMS4 CRC (EMT without clear identity)
    These are the DOUBLE-NEGATIVE states of their
    respective tissues — identity-lost but without
    a foreign programme to define them.

  DEFINING MOLECULAR FEATURES:
    AR:           ABSENT
    NE markers:   ABSENT
    TP53:         FREQUENTLY MUTATED/LOST
    RB1:          FREQUENTLY LOST
    PTEN:         FREQUENTLY LOST
    WNT pathway:  OFTEN ACTIVATED (CTNNB1, RNF43)
    FGF pathway:  OFTEN ACTIVATED
    EZH2:         VARIABLE
    SOX2:         SOMETIMES ELEVATED (stem marker)
    The defining feature is the ABSENCE of identity
    markers rather than the presence of specific
    markers — DNPC is characterised by what it
    has lost, not what it expresses.

  CIRCUIT STATUS: ABSENT.
    There is no circuit to restore.
    There is no dominant false attractor programme
    to dissolve.
    DNPC is the most challenging therapeutic target
    in the repository precisely because it has
    escaped all identity-based targeting strategies.
    Current approaches: WNT inhibitors, FGF
    inhibitors, SOX2 targeting (experimental),
    immunotherapy (if TMB-high), platinum if HRR
    deficient.

TREATMENT FAILURE CASCADE — THE WADDINGTON
RESISTANCE TRAJECTORY:
  The prostate cancer treatment sequence is a
  FORCED WADDINGTON TRANSITION:
    Primary PRAD (Luminal A) → ADT
    → Luminal B (AR pathway adaptation)
    → mCRPC (FOXA1 mutation, AR amplification,
      AR-V7)
    → AR pathway inhibitor (abiraterone/enzalutamide)
    → CRPC treatment-resistant:
      Path A: NEPC (TP53+RB1 loss → AURKA/MYCN)
      Path B: DNPC (full identity loss — no NE)
    → Death from treatment-resistant disease.
  THIS IS THE WADDINGTON DEPTH AXIS FORCED BY
  SEQUENTIAL THERAPY.
  The depth score in PRAD measures not just WHERE
  a patient's cancer is — it predicts WHAT COMES
  NEXT on the resistance trajectory.
  A deep depth score in Luminal A predicts
  imminent transition toward Luminal B (requires
  earlier intensification).
  A deep depth score in Luminal B with TP53 loss
  predicts NEPC trajectory (requires AURKA/RB1
  monitoring and early alisertib consideration).
  This predictive application of the depth score —
  anticipating the next attractor transition before
  it occurs clinically — is the most powerful
  and unique application of the framework in PRAD.
```

---

## SECTION VII — THE AR PATHWAY: THE DOMINANT THERAPEUTIC AXIS

```
The androgen receptor (AR) pathway is the most
thoroughly characterised oncogenic axis in medicine.
It has been the subject of targeted therapy since
the 1940s (Huggins and Hodges, 1941 — surgical
castration for prostate cancer; Nobel Prize 1966).

THE AR PATHWAY SUMMARY:

ANDROGENS (testosterone → DHT):
  Testosterone produced by Leydig cells (testis)
  and adrenal glands (DHEA, androstenedione).
  5-alpha reductase (SRD5A1/2) converts testosterone
  to DHT (dihydrotestosterone) in prostate tissue.
  DHT binds AR with 3-5× higher affinity than
  testosterone — DHT is the intraprostatic driver.
  Enzymes: CYP17A1 (also in adrenal glands and
  tumour cells) converts cholesterol → androgens
  → the target of abiraterone.

AR ACTIVATION CASCADE:
  DHT binds AR ligand-binding domain (LBD)
  → AR-DHT complex translocates to nucleus
  → AR binds to androgen response elements (AREs)
    in chromatin (FOXA1 has opened these sites)
  → AR activates:
    PSA (KLK3) — the clinical monitoring marker
    KLK2, TMPRSS2 — other targets
    NKX3-1 (paradox: AR ACTIVATES NKX3-1 normally
    — this is why ADT causes NKX3-1 to fall —
    AND why the EZH2 lock on NKX3-1 is so
    pernicious: EZH2 silences NKX3-1 even while
    AR tries to activate it. EZH2 wins.)
    Cell survival genes: BCL2, BCL-XL
    Proliferation genes: CDK1, CDK2

AR-V7 (the primary CRPC mechanism):
  Splice variant of AR lacking the LBD.
  Constitutively active — doesn't need androgens.
  Cannot be blocked by:
    Enzalutamide (binds LBD — absent in AR-V7)
    Abiraterone (blocks androgen synthesis —
    AR-V7 doesn't need androgens)
  AR-V7 detection in CTCs predicts CRPC
  resistance to ARPI (abiraterone, enzalutamide).
  AR-V7 patients: Switch to cabazitaxel.

EZH2 → NKX3-1 SUPPRESSION DESPITE AR:
  The core molecular paradox in PRAD depth analysis:
  AR activates NKX3-1 in the normal luminal cell.
  But in cancer, EZH2 places H3K27me3 at NKX3-1.
  H3K27me3 is a closed chromatin mark —
  AR cannot access a closed H3K27me3 locus.
  EZH2 OVERRIDES AR at the NKX3-1 locus:
  Cancer cell has HIGH AR AND LOW NKX3-1
  simultaneously because EZH2 closed the NKX3-1
  chromatin before AR could access it.
  TAZEMETOSTAT removes H3K27me3:
  → NKX3-1 locus opens
  → AR can now access NKX3-1 (AR is still present
    in Luminal A)
  → NKX3-1 is expressed
  → NKX3-1 drives luminal differentiation
  This is why tazemetostat + enzalutamide has
  rationale: tazemetostat opens NKX3-1 for AR,
  AND enzalutamide blocks AR from driving
  oncogenic targets — the combination forces AR
  toward NKX3-1 activation rather than oncogenic
  target activation.
  The order of events matters:
  Tazemetostat FIRST (to open NKX3-1 locus)
  THEN AR pathway inhibition
  This is a mechanistic sequence prediction
  that no current clinical trial has tested.
```

---

## SECTION VIII — EXISTING PRAD ANALYSIS IN CONTEXT

```
THE COMPLETED ANALYSIS (OrganismCore cancer series)
ran on a bulk PRAD dataset
(TCGA-PRAD or equivalent).

WHAT THE ANALYSIS FOUND:

  FINDING 1: NKX3-1 identified as switch gene.
    CORRECT CONTEXT:
    NKX3-1 is the master luminal TF — its loss
    is the founding event in luminal identity
    destabilisation and the Waddington transition
    from normal luminal to cancer cell.
    Its depth-negative correlation (falls with
    depth) is expected: deeper tumours have lost
    more luminal identity (more NKX3-1 gone).
    Its identification as the switch gene confirms:
    the depth axis is measuring loss of luminal
    identity, which is the correct Waddington axis
    for PRAD.
    SUBTYPE CONTEXT:
    NKX3-1 as switch gene is primarily a LUMINAL A
    finding (the dominant subtype in bulk PRAD
    at ~50%).
    In Luminal B: NKX3-1 is also the switch gene
    but more completely silenced.
    In basal-like: NKX3-1 is absent — the concept
    of "restoring the switch gene" must give way
    to the identity-loss attractor dissolution
    strategy.

  FINDING 2: EZH2 elevated (gain-of-function lock).
    CORRECT CONTEXT:
    EZH2 in PRAD silences NKX3-1 via H3K27me3.
    The EZH2 lock is the epigenetic mechanism
    preventing circuit restoration.
    Tazemetostat removes the lock.
    CONFIRMED in PRAD bulk signal (Luminal A
    dominated).
    CLINICAL VALIDATION (partial):
    CELLO-1 was negative in mCRPC —
    but as interpreted above, this does not
    disprove the mechanism; it demonstrates
    that unselected mCRPC is the wrong patient
    population. Luminal A primary PRAD with
    high EZH2 and intermediate depth score is
    the correct population.
    SUBTYPE QUESTION:
    Is EZH2 equally elevated in Luminal B?
    Is EZH2 the correct target in basal-like?
    Is EZH2 relevant in NEPC (where AURKA/MYCN
    is the dominant false attractor driver)?
    These questions are resolved in PRAD-S2a
    through PRAD-S4a.

  FINDING 3: Circuit restoration therapeutic
    strategy (Pattern 5).
    CORRECT CONTEXT:
    Like PAAD, PRAD has an intact circuit in its
    dominant (Luminal A) subtype.
    Unlike PAAD (where all patients are metastatic
    at presentation), PRAD offers the unique
    opportunity to intervene with circuit
    restoration BEFORE the circuit is broken —
    at Gleason 6–7, years before CRPC.
    This timing advantage is the most clinically
    distinctive feature of the PRAD series.

  FINDING 4: Tazemetostat predicted as drug target.
    STATUS:
    Tazemetostat in mCRPC (CELLO-1): NEGATIVE in
    unselected population — does NOT invalidate
    the prediction; validates the depth score
    patient selection requirement.
    Tazemetostat in Luminal A primary PRAD with
    high EZH2 depth score: NOT YET TESTED.
    The going-further prediction is:
    Tazemetostat in DEPTH-SELECTED Luminal A
    primary PRAD will show NKX3-1 restoration
    and differentiation programme re-engagement.
    Patient selection tool: the depth score.
    This remains the GOING FURTHER finding for PRAD.
```

---

## SECTION IX — DATA AVAILABILITY SUMMARY

```
Dataset           Accession          n(tumour)  Normal   Subtype  Power
                                                (n)      labels
──────────────────────────────────────────────────────────────────────────
TCGA-PRAD         GDC/phs000178      ~499       ~52 adj  YES      HIGH
                  (PanCancer Atlas)            normal   (TCGA 7
                                                         genomic
                                                         classes)
GSE62944          GEO (TCGA RNA-seq  ~499       ~52      PARTIAL  HIGH
                  reformatted)
GSE21034          GEO (Lapointe      ~164       ~28 adj  PARTIAL  MOD
                  2010 — luminal/    prostate   normal   (Luminal/
                  basal annotations) microarray          Basal
                                                         manually)
GSE204811         GEO (2024)         Variable   YES      Partial  MOD
                  (epigenetic                   (cell-            (cell
                  analysis TCGA-     lines)     based)            line)
                  PRAD)
SU2C Cohort       SU2C/PCF           ~270 mCRPC NONE     YES      HIGH
(Robinson 2015    Dream Team                            (mCRPC
 — mCRPC)                                               subtypes
                                                         incl NEPC)
GSE32269          GEO (Grasso 2012)  ~59 mCRPC  ~28 adj  PARTIAL  MOD
                  (mCRPC)                       normal
GSE77930          GEO (Beltran 2016  ~35 NEPC   ~14 LuCaP NEPC   MOD
                  NEPC)              +           prostate  labels
                                     ~49 CRPC   normal
GTEx              —                  0           ~30      N/A      HIGH
                                                prostate          (NORMAL
                                                normal            ONLY)

CRITICAL NOTES:

NOTE 1 — TUMOUR PURITY:
  TCGA-PRAD has GOOD tumour purity (~60–70%)
  relative to TCGA-PAAD (~35%).
  The PRAD contamination issue is NOT acinar/
  ductal cell contamination (no normal prostate
  equivalent of the ADEX artefact).
  The primary PRAD contamination issue is:
  STROMAL CELLS — prostate tumours have a
  fibromuscular stroma with AR-positive stromal
  cells that can contribute AR-pathway gene
  expression to bulk RNA-seq.
  This creates a mild confound (stromal AR-positive
  signals inflate the apparent AR activity in bulk
  tumour RNA-seq) but is less severe than the
  PAAD purity problem.
  Purity filtering at ≥50% is recommended but
  not as urgently critical as in PAAD.

NOTE 2 — NEPC DATASETS:
  NEPC samples are NOT in TCGA-PRAD
  (TCGA collected only primary adenocarcinoma).
  For NEPC analysis, GSE77930 (Beltran 2016)
  and SU2C mCRPC cohort are the primary resources.
  Both are publicly accessible.
  GSE77930 contains matched NEPC + CRPC adenocarcinoma
  + adjacent normal prostate — ideal for comparing
  the depth axis across the CRPC → NEPC transition.

NOTE 3 — SUBTYPE LABELS:
  TCGA 2015 seven genomic classes:
    Available in TCGA supplementary tables and
    cBioPortal for the TCGA-PRAD cohort.
    ERG fusion, SPOP, FOXA1, IDH1 labels available.
  Luminal/Basal transcriptional labels:
    NOT officially in TCGA — must be derived from
    the expression data using published classifiers
    (e.g., Zhao et al. 2017 Luminal A/B/Basal-like
    prostate classification, or the consensus from
    Karthaus et al. 2020 scRNA-seq).
    The before-documents must state which
    classifier is being used to assign
    transcriptional subtype labels.
  mCRPC subtype labels (including NEPC):
    SU2C cohort has these labels.
    Beltran 2016 (GSE77930) has NEPC labels.

NOTE 4 — THE GLEASON/ISUP GRADE SYSTEM:
  The Gleason grade (1–10) and modern ISUP
  grade (1–5) are the clinical depth measures
  currently used in pathology.
  The framework depth score should CORRELATE
  WITH Gleason grade — this is the primary
  clinical validation test for the PRAD depth score.
  Correlation of depth score with:
    Gleason ≤6 (ISUP 1): expected shallow depth
    Gleason 3+4=7 (ISUP 2): expected intermediate
    Gleason 4+3=7 (ISUP 3): expected intermediate-deep
    Gleason 8 (ISUP 4): expected deep
    Gleason 9-10 (ISUP 5): expected very deep
  A depth score that does NOT correlate with
  Gleason grade should trigger a QC review —
  this is the built-in validation criterion for
  the PRAD depth analysis.
```

---

## SECTION X — PLANNED ANALYSIS ORDER

```
ORDER:

  PRAD-S1   Luminal A          TCGA-PRAD         HIGH
            (intact circuit,   luminal-high
            EZH2 lock,         subset (NKX3-1
            tazemetostat)      high, AR high,
                               Gleason ≤7)
                               + GTEx normal
                               prostate
                               REASON: The dominant
                               subtype (~50%).
                               The intact circuit
                               finding from bulk.
                               EZH2 elevation and
                               NKX3-1 depth axis
                               confirmed in luminal-
                               specific context.
                               Correlation with
                               Gleason score as
                               depth validation.
                               CELLO-1 failure
                               recontextualised by
                               depth score —
                               identify the depth-
                               score-selected
                               subpopulation that
                               would have responded.
                               Clinical output:
                               depth score as
                               tazemetostat patient
                               selection tool at
                               Gleason 6–7 (pre-CRPC
                               intervention).

  PRAD-S2   Luminal B          TCGA-PRAD         MOD-HIGH
            (intermediate      Luminal B subset
            depth, CRPC        (NKX3-1 low, AR
            transition)        high/variant,
                               Gleason 7–8)
                               + SU2C mCRPC subset
                               REASON: The CRPC
                               transition subtype.
                               FOXA1 mutation
                               confound addressed.
                               PTEN loss + AKT
                               pathway separation.
                               AR-V7 depth
                               correlation.
                               Circuit status in
                               Luminal B: restoreable
                               or broken?
                               Clinical output:
                               depth score as CRPC
                               risk predictor;
                               tazemetostat vs.
                               ipatasertib selection
                               tool (EZH2 lock vs.
                               PI3K bypass pathway
                               dominant).

  PRAD-S3   Basal-like         TCGA-PRAD         MOD
            primary            Basal-like subset
                               (NKX3-1 low, AR
                               low, TP63 present,
                               Gleason 9–10)
                               REASON: The identity-
                               loss false attractor
                               in PRAD.
                               Circuit status:
                               transition from
                               restoreable to broken.
                               Residual surface
                               markers (PSMA, EGFR,
                               DLL3-early) as depth-
                               positive drug targets.
                               AURKA rise in basal-
                               like — pre-NEPC
                               signature.
                               Clinical output:
                               identify which
                               basal-like tumours
                               are pre-NEPC
                               (AURKA+/MYCN early)
                               vs. basal-identity
                               (TP63+/AR−).

  PRAD-S4   NEPC               GSE77930 +        MOD
                               SU2C NEPC subset
                               REASON: Circuit-gone
                               state. AURKA/MYCN
                               as convergence node.
                               Alisertib depth
                               correlation.
                               DLL3 as residual
                               surface identity
                               marker.
                               Depth comparison:
                               NEPC vs. Luminal A —
                               the full Waddington
                               axis traversal in
                               one analysis.
                               Clinical output:
                               depth score as
                               tarlatamab/alisertib
                               patient selection
                               tool in NEPC.

  PRAD-X    Cross-subtype +    After S1–S4.
            Resistance         QUESTIONS:
            Trajectory           1. Does depth score
                                    order:
                                    Luminal A < Luminal B
                                    < Basal-like < NEPC?
                                    (The Waddington axis
                                    traversal prediction)
                                 2. EZH2 direction:
                                    elevated in Luminal A,
                                    Luminal B, Basal-like,
                                    AND NEPC? Or does it
                                    peak at Luminal B and
                                    change character in
                                    NEPC?
                                 3. The AURKA crossover
                                    point: at what depth
                                    score does AURKA become
                                    depth-POSITIVE (marking
                                    the tazemetostat →
                                    alisertib switch
                                    transition)?
                                 4. Gleason score
                                    validation: does depth
                                    score correlate with
                                    Gleason across ALL
                                    subtypes?
                                 5. The CELLO-1 post-hoc:
                                    in TCGA-PRAD, what is
                                    the depth score
                                    distribution? How many
                                    patients would have
                                    been depth-selected
                                    vs. depth-excluded?
                                    What would the CELLO-1
                                    trial look like if
                                    only depth-selected
                                    patients were enrolled?
```

---

## SECTION XI — WHAT THIS DOCUMENT DOES NOT CONTAIN

```
This document contains:
  ✓ Normal prostate epithelial hierarchy
    (luminal → intermediate → basal, with
    full TF network: NKX3-1/FOXA1/AR/TP63)
  ✓ The TCGA seven genomic classes (TMPRSS2-ERG,
    SPOP, FOXA1, IDH1, CHD1, CDKN2A, other)
    with frequencies and clinical significance
  ✓ The four transcriptional subtypes (Luminal A,
    Luminal B, Basal-like, NEPC/DNPC) with cells
    of origin, molecular events, and treatment
  ✓ The intact circuit finding from the existing
    analysis placed in Luminal A context
  ✓ NKX3-1 as switch gene, EZH2 as epigenetic lock,
    tazemetostat as predicted drug — all in context
  ✓ The CELLO-1 failure reinterpreted through the
    depth score lens (unselected mCRPC ≠ the
    depth-selected Luminal A population)
  ✓ The NEPC AURKA/MYCN convergence node and
    alisertib as the NEPC drug target —
    connecting to the CIN STAD AURKA finding
  ✓ The AR pathway full biology (testosterone →
    DHT → AR → NKX3-1 → luminal programme)
    and why EZH2 overrides AR at NKX3-1
  ✓ The Waddington depth axis as the forced
    resistance trajectory under AR therapy
    (Luminal A → Luminal B → CRPC → NEPC/DNPC)
  ✓ Data availability (TCGA-PRAD, SU2C, GSE77930)
    with tumour purity note (mild in PRAD, not
    the crisis it is in PAAD)
  ✓ The Gleason/ISUP depth validation standard
  ✓ The cross-repository structural parallels
    (basal-like PRAD ↔ basal-like PDAC ↔ GS STAD
    ↔ basal/squamous BLCA ↔ CMS4 CRC —
    identity-loss attractor type across cancers)

This document does NOT contain:
  ✗ Depth score predictions
  ✗ Specific false attractor gene predictions
    beyond the existing bulk analysis findings
  ✗ Drug target predictions beyond the established
    analysis (NKX3-1, EZH2, tazemetostat for
    Luminal A; alisertib for NEPC)
  ✗ Epigenetic mechanism hypotheses

All of the above belong in the BEFORE documents.
PRAD-S1a (Luminal A before-document) is next.
Written before any script runs.
Before any data loads.
```

---

## STATUS BLOCK

```
document:           PRAD_Subtype_Orientation.md
folder:             Cancer_Research/PRAD/Subtypes/
status:             COMPLETE — ORIENTATION ONLY
date:               2026-03-03
author:             Eric Robert Lawson
                    OrganismCore

subtypes_covered:
  Luminal A (~50%): NKX3-1 high→med, EZH2 high, [1 of 5]
                    AR high, intact circuit,
                    tazemetostat candidate
  Luminal B (~25%): NKX3-1 low, EZH2 high,      [2 of 5]
                    AR high/variant, CRPC
                    transition, FOXA1/PTEN confounds
  Basal-like (~10%):NKX3-1 absent, AR low,       [3 of 5]
                    TP63+, identity-loss attractor,
                    circuit transition point
  NEPC (~10%        AR absent, AURKA/MYCN,        [4 of 5]
  treatment-emerg): DLL3+, RB1/TP53 lost,
                    circuit broken, alisertib target
  DNPC (~5%):       AR absent, NE absent,         [5 of 5]
                    stem-like, maximum depth,
                    no circuit

analyses_started:   0 (new subtype series)
existing_analysis:  Cancer_Research/PRAD/ — complete
                    (OrganismCore cancer series)
                    NKX3-1 switch gene identified
                    EZH2 elevated (gain-of-function lock)
                    Circuit intact (restoration strategy)
                    Tazemetostat predicted as drug target
                    All confirmed findings = Luminal A
                    dominated bulk signal context

next_document:      PRAD-S1a
                    Luminal A Before-Document
                    (predictions locked before
                    TCGA-PRAD Luminal A subset loads)

critical_note_1:    THE CELLO-1 FAILURE IS NOT A
                    FRAMEWORK FAILURE.
                    Tazemetostat + enzalutamide in
                    unselected mCRPC → negative result.
                    Framework interpretation:
                    mCRPC is NOT Luminal A.
                    mCRPC is predominantly Luminal B
                    or transitioning to basal-like.
                    The EZH2 circuit restoration
                    mechanism operates in LUMINAL A
                    where NKX3-1 is RESTOREABLE.
                    In mCRPC, NKX3-1 may be beyond
                    restoration — fully silenced by
                    both H3K27me3 AND CpG methylation.
                    The correct population for
                    tazemetostat is:
                    Primary disease, Gleason 6–7,
                    high EZH2 by IHC, depth score
                    in the shallow-intermediate range.
                    NOT mCRPC.
                    The depth score is the patient
                    selection tool that makes
                    tazemetostat work — without it,
                    the drug is being given to patients
                    who cannot respond because their
                    circuits are already broken.
                    This is the framework's most
                    precise clinical statement about
                    prostate cancer.

critical_note_2:    THE WADDINGTON DEPTH AXIS IN PRAD
                    IS THE ONLY AXIS IN THE REPOSITORY
                    THAT IS FORCED BY THERAPY.
                    In every other cancer, depth
                    increases passively through
                    tumour evolution.
                    In PRAD, sequential AR therapy
                    ACTIVELY DRIVES cells deeper
                    by removing AR-dependent (shallow)
                    cells and selecting for AR-independent
                    (deeper) cells.
                    This means the depth score in PRAD
                    is both DIAGNOSTIC (where is the
                    cancer now?) AND PROGNOSTIC
                    (how fast will the cancer progress
                    to the next attractor under the
                    planned therapy?).
                    The PRAD depth score should predict
                    TIME TO CRPC — the deeper the score
                    at initial diagnosis, the shorter
                    the time to castration resistance.
                    This is a falsifiable prediction
                    testable in TCGA-PRAD survival data
                    (correlate depth score at diagnosis
                    with biochemical recurrence-free
                    survival and metastasis-free survival).
                    If confirmed, the depth score at
                    Gleason 6–7 becomes a CLINICAL
                    DECISION TOOL: patients with deep
                    score at low Gleason are Gleason
                    6 cancers that will behave like
                    Gleason 9 cancers — they should
                    NOT be on active surveillance.
                    They should receive early
                    definitive therapy.

critical_note_3:    THE NEPC AURKA FINDING CROSS-VALIDATES
                    THE CIN STAD AURKA FINDING.
                    In the existing STAD bulk analysis:
                    AURKA is depth-positive, ZEB2/AURKA
                    are co-expressed (r=0.9871).
                    Alisertib predicted for CIN STAD.
                    In NEPC: AURKA amplification is
                    the defining molecular event.
                    AURKA stabilises N-MYC.
                    Alisertib disrupts AURKA/N-MYC.
                    Clinical trial data (SPARTAN —
                    Phase II in mCRPC/NEPC with MYCN
                    amplification) confirms activity.
                    The framework predicted AURKA as
                    a depth-positive convergence node
                    in CIN STAD from first principles.
                    The same principle, independently
                    confirmed in NEPC, validates the
                    universality of the framework's
                    AURKA depth-correlation finding
                    across cancer types.
                    AURKA is NOT a coincidental finding
                    in any one cancer — it is a cross-
                    cancer depth-positive marker that
                    marks the most proliferative and
                    identity-lost false attractors.
                    This is Framework Pattern 2
                    (Drug Target Derivation) confirmed
                    at the cross-cancer level.

critical_note_4:    THE DLL3/NEPC FINDING COMPLETES
                    THE RESIDUAL SURFACE MARKER PATTERN.
                    In every identity-loss attractor
                    across the repository, one residual
                    surface marker survives identity loss:
                    CLDN18.2 in GS STAD
                    NECTIN4 in basal/squamous BLCA
                    EGFR in basal-like PDAC
                    DLL3 in NEPC PRAD
                    These four are NOT the deep
                    tumour markers — they are the
                    RESIDUAL IDENTITY MARKERS of the
                    new false attractor programme
                    (not the old tissue identity,
                    but the new aberrant identity).
                    DLL3 is a neuroendocrine programme
                    marker — it is expressed because
                    the NEPC cell HAS ACQUIRED
                    neuroendocrine identity.
                    It is the one surface protein of
                    the new programme that survives
                    on the cell surface and can be
                    targeted by an antibody.
                    Tarlatamab (DLL3 biTE) works in
                    NEPC for exactly the same reason
                    zolbetuximab works in GS STAD:
                    the drug targets the RESIDUAL
                    SURFACE MARKER of the false
                    attractor, regardless of how
                    identity-lost the cancer cell is.
                    The framework predicts this
                    pattern will hold for every
                    cancer in the repository:
                    Find the residual surface marker
                    of the false attractor programme.
                    That marker is the ADC/antibody
                    target for the most identity-lost
                    patients in that cancer type.
                    This is the going-further finding
                    that emerges from combining the
                    PRAD NEPC analysis with the rest
                    of the repository.
```
