# FOXA1/EZH2 Dual IHC Breast Cancer Classifier
**One-Page Summary** — Eric Robert Lawson / OrganismCore
`OrganismCore@proton.me` | ORCID: `0009-0002-0414-6544`

---

## How This Was Derived — Read This First

This test was **not found by mining clinical data**.
It was derived from first principles by a mathematician, not a clinician.

**Step 1 — Geometry first:**
Waddington landscape attractor geometry was applied to 19,542 single breast
cancer cells (GSE176078, Wu et al. 2021, *Nature Genetics*). From the geometry
alone, FOXA1 and EZH2 were identified as mechanistic opposites on the primary
identity axis of breast cancer — **before any clinical dataset was opened.**

**Step 2 — Prediction locked:**
The predicted subtype ordering — LumA > LumB > HER2 > TNBC > Claudin-low —
was written down in a version-controlled document with a commit timestamp
**before any confirmatory analysis was run.**

**Step 3 — Confirmation:**
Seven independent clinical datasets were analysed. The prediction confirmed
exactly in all seven, in the predicted direction, across four measurement
platforms and ~7,500 patients.

> **The statistics below are confirmations of a pre-specified prediction.**
> **They are not the source of it.**
>
> This eliminates overfitting, data dredging, and post-hoc hypothesis
> adjustment as explanations for the results. The timestamped before-document
> is publicly verifiable in the repository.

---

## The Biological Mechanism

| Protein | Role | High means |
|---------|------|------------|
| **FOXA1** | Pioneer transcription factor — opens chromatin for luminal breast cell identity | Cell identity programme is active. Luminal. |
| **EZH2** | PRC2 catalytic subunit — deposits H3K27me3 silencing marks on FOXA1, GATA3, and ESR1 promoters | Active epigenetic suppression of luminal identity. |

Their **ratio measures which force is winning** — identity or its suppression.
That balance determines both subtype and treatment vulnerability.

*Mechanism independently confirmed by Schade et al. (Nature, 2024) and Toska
et al. (Nature Medicine, 2017) — neither group had knowledge of this framework.*

---

## The Test

```
FOXA1 IHC H-score ÷ EZH2 IHC H-score
```

- Two standard antibodies. One ratio.
- Both antibodies already in routine clinical use worldwide.
- Standard IHC protocol. Same-day results.
- Compatible with whole-section FFPE and TMA format.
- **Cost: ~$50–$100 per patient.**

---

## What It Classifies

| Subtype | Ratio character | Treatment logic |
|---------|----------------|-----------------|
| Luminal A | Highest ratio | CDK4/6 inhibitor + endocrine therapy |
| Luminal B | High ratio, below LumA | HDACi + endocrine therapy |
| HER2-enriched | Mid-range ratio | Anti-HER2 first, ET thereafter |
| TNBC | Low ratio | Tazemetostat → fulvestrant |
| Claudin-low | Lowest ratio | Anti-TIGIT → anti-PD-1 sequence |
| ILC *(exception)* | Inverted — above LumA | Fulvestrant > aromatase inhibitor |

---

## Confirmations of the Pre-Specified Prediction

*~7,500 patients · 7 independent datasets · 4 measurement platforms*

| Analysis | Result |
|----------|--------|
| Subtype ordering | TCGA-BRCA n=837, p=2.87×10⁻¹⁰³ |
| LumA vs LumB separation | METABRIC n=1,980, p=8.47×10⁻⁶⁷ |
| Survival stratification | Confirmed in 4 independent cohorts (~7,500 patients) |
| ROC — LumA vs Basal | AUC 0.828–0.901 *(replicated in METABRIC and TCGA independently)* |
| ROC — LumA vs LumB | AUC 0.796–0.873 *(replicated in METABRIC and TCGA independently)* |
| Protein-level confirmation | CPTAC mass spectrometry n=122, Spearman r=−0.492, p<0.0001 |
| Treatment response — chemo | GSE25066 n=508, p<0.0001 |
| Treatment response — ET | METABRIC n=1,104, p=0.0027, Δ=18.7 months RFS |
| **Biological contradictions** | **0** |

---

## What Is Not Yet Established

> **IHC H-score cut-points do not yet exist.**

The AUC and p-values above are from RNA-level and protein-level computational
validation. The H-score thresholds for clinical classification must be
determined from actual stained tissue against PAM50 ground truth using
Youden J optimisation.

This is **instrument calibration** — the same step applied to ER, PR, HER2,
and Ki67 before routine clinical adoption. It is the only remaining step.

---

## The Ask

| Item | Detail |
|------|--------|
| Tissue | ~300–400 archived FFPE cases with known PAM50 subtype (whole section **or** TMA) |
| Stains | FOXA1 (CST D4E2 #53528) and EZH2 (CST D2H1 #5246) |
| Scoring | H-score blinded to PAM50 subtype |
| Analysis | Concordance against PAM50; cut-points derived from scratch |
| Patients | No new recruitment. No treatment changes. |
| Cost | ~$5,000–$15,000 consumables |
| Timeline | 3–6 months from tissue access to submission |

---

## In Return

- Complete **Protocol Specification v3.0** — antibody clones, antigen retrieval
  conditions (Ventana / Dako / Leica), H-score method, TMA instructions,
  claudin-low stromal contamination warning
- Full **statistical analysis framework** — ROC, Youden J cut-point
  optimisation, training/test split, bootstrapped CIs, inter-observer kappa
- Complete public computational validation dataset and scripts
  *(reproducible from public data)*
- **Co-authorship** on the resulting publication
- All materials open access (MIT license)

---

## Why This Matters

PAM50/Prosigna costs **$3,000–$4,000**, requires proprietary RNA equipment
($150k–$250k capital cost), and is inaccessible to the majority of breast
cancer patients on earth.

A validated FOXA1/EZH2 IHC protocol delivers equivalent classification
information — plus mechanism, plus treatment logic — from two standard
antibodies at **$50–$100 per patient**, in any laboratory that processes
biopsies, on the same day as the biopsy workup.

**~1.2–1.5 million breast cancer patients per year currently receive no
molecular subtype information. This test, on validation, reaches all of them.**

---

## Published Evidence

| Document | DOI / Link |
|----------|-----------|
| CS-LIT-1 — Ratio validation (~7,500 patients) | [doi.org/10.5281/zenodo.18883922](https://doi.org/10.5281/zenodo.18883922) |
| CS-LIT-22 — IHC decision tool | [doi.org/10.5281/zenodo.18892788](https://doi.org/10.5281/zenodo.18892788) |
| Full repository (open access) | [github.com/Eric-Robert-Lawson/attractor-oncology](https://github.com/Eric-Robert-Lawson/attractor-oncology) |

---

*Eric Robert Lawson · OrganismCore · OrganismCore@proton.me · ORCID: 0009-0002-0414-6544*
