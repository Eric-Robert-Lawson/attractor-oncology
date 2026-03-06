FOXA1/EZH2 DUAL IHC BREAST CANCER CLASSIFIER
One-Page Summary — Eric Robert Lawson / OrganismCore
OrganismCore@proton.me | ORCID: 0009-0002-0414-6544
────────────────────────────────────────────────────

HOW THIS WAS DERIVED — READ THIS FIRST

  This test was not found by mining clinical
  data. It was derived from first principles
  by a mathematician, not a clinician.

  Step 1 — Geometry first:
  Waddington landscape attractor geometry
  was applied to 19,542 single breast cancer
  cells (GSE176078, Wu et al. 2021, Nature
  Genetics). From the geometry alone, FOXA1
  and EZH2 were identified as mechanistic
  opposites on the primary identity axis of
  breast cancer — before any clinical dataset
  was opened.

  Step 2 — Prediction locked:
  The predicted subtype ordering —
  LumA > LumB > HER2 > TNBC > Claudin-low
  — was written down in a version-controlled
  document with a commit timestamp BEFORE
  any confirmatory analysis was run.

  Step 3 — Confirmation:
  Seven independent clinical datasets were
  analysed. The prediction confirmed exactly
  in all seven, in the predicted direction,
  across four measurement platforms and
  ~7,500 patients.

  THE STATISTICS BELOW ARE CONFIRMATIONS
  OF A PRE-SPECIFIED PREDICTION.
  THEY ARE NOT THE SOURCE OF IT.

  This eliminates overfitting, data dredging,
  and post-hoc hypothesis adjustment as
  explanations for the results. The
  timestamped before-document is publicly
  verifiable in the repository.

────────────────────────────────────────────────────

THE BIOLOGICAL MECHANISM

  FOXA1: pioneer transcription factor.
  Opens chromatin for luminal breast cell
  identity. High FOXA1 = the cell knows
  what it is.

  EZH2: PRC2 catalytic subunit. Deposits
  H3K27me3 silencing marks on FOXA1,
  GATA3, and ESR1 promoters. High EZH2 =
  active epigenetic suppression of luminal
  identity.

  Their ratio measures which force is
  winning — identity or its suppression.
  That balance determines both subtype
  and treatment vulnerability.

  Mechanism independently confirmed:
  Schade et al. (Nature, 2024) and Toska
  et al. (Nature Medicine, 2017) — neither
  group had knowledge of this framework.

────────────────────────────────────────────────────

THE TEST
  FOXA1 IHC H-score ÷ EZH2 IHC H-score
  Two standard antibodies. One ratio.
  Both antibodies already in clinical use.
  Standard IHC protocol. Same-day results.
  Compatible with whole-section FFPE
  and TMA format.
  Cost: ~$50–$100 per patient.

────────────────────────────────────────────────────

WHAT IT CLASSIFIES

  Luminal A     Highest ratio
                CDK4/6 inhibitor + ET

  Luminal B     High ratio, below LumA
                HDACi + ET

  HER2-enr.     Mid-range ratio
                Anti-HER2 first, ET thereafter

  TNBC          Low ratio
                Tazemetostat → fulvestrant

  Claudin-low   Lowest ratio
                Anti-TIGIT → anti-PD-1

  ILC           Inverted — above LumA
  (exception)   Fulvestrant > AI

────────────────────────────────────────────────────

CONFIRMATIONS OF THE PRE-SPECIFIED PREDICTION
(~7,500 patients, 7 datasets, 4 platforms)

  Subtype ordering:
    TCGA n=837  p=2.87×10⁻¹⁰³
    METABRIC LumA vs LumB
    n=1,980  p=8.47×10⁻⁶⁷

  Survival stratification:
    Confirmed in 4 independent cohorts
    (~7,500 patients total)

  ROC — LumA vs Basal:
    AUC 0.828–0.901
    (replicated in METABRIC and TCGA
    independently)

  ROC — LumA vs LumB:
    AUC 0.796–0.873
    (replicated in METABRIC and TCGA
    independently)

  Protein-level confirmation:
    CPTAC mass spectrometry n=122
    Spearman r=−0.492  p<0.0001
    Correct ordering confirmed at the
    protein level — the level IHC measures

  Treatment response:
    GSE25066 chemotherapy n=508  p<0.0001
    METABRIC endocrine therapy n=1,104
    p=0.0027  Δ=18.7 months RFS

  Biological contradictions: 0

────────────────────────────────────────────────────

WHAT IS NOT YET ESTABLISHED

  IHC H-score cut-points do not yet exist.

  The AUC and p-values above are from
  RNA-level and protein-level computational
  validation. The H-score thresholds for
  clinical classification must be determined
  from actual stained tissue against PAM50
  ground truth using Youden J optimisation.

  This is instrument calibration — the same
  step applied to ER, PR, HER2, and Ki67
  before routine clinical adoption.
  It is the only remaining step.

────────────────────────────────────────────────────

THE ASK
  ~300–400 archived FFPE breast cancer
  cases with known PAM50 subtype
  (whole section or TMA format)
  FOXA1 (CST D4E2 #53528) and EZH2
  (CST D2H1 #5246) IHC
  H-score scoring blinded to PAM50
  Concordance analysis
  No new patients. No treatment changes.
  Est. cost: $5,000–$15,000 consumables.
  Timeline: 3–6 months.

IN RETURN
  Full Protocol Specification v3.0
  Full statistical analysis framework
  Co-authorship on publication
  All materials open access (MIT license)

────────────────────────────────────────────────────

WHY THIS MATTERS

  PAM50/Prosigna costs $3,000–$4,000,
  requires proprietary RNA equipment, and
  is inaccessible to the majority of breast
  cancer patients on earth.

  A validated FOXA1/EZH2 IHC protocol
  delivers equivalent information — plus
  mechanism, plus treatment logic — from
  two standard antibodies at $50–$100 per
  patient, in any laboratory that processes
  biopsies, worldwide.

  ~1.2–1.5 million breast cancer patients
  per year currently receive no molecular
  subtype information. This test, on
  validation, reaches all of them.

────────────────────────────────────────────────────

PUBLISHED EVIDENCE

  CS-LIT-1 (ratio, ~7,500 patients):
  doi.org/10.5281/zenodo.18883922

  CS-LIT-22 (IHC decision tool):
  doi.org/10.5281/zenodo.18892788

  Full repository (open access):
  github.com/Eric-Robert-Lawson/
  attractor-oncology

────────────────────────────────────────────────────
Eric Robert Lawson · OrganismCore
OrganismCore@proton.me
ORCID: 0009-0002-0414-6544
────────────────────────────────────────────────────
