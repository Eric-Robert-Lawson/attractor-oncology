# Attractor Oncology
### Waddington Landscape Geometry Applied to Cancer Patient Stratification
#### A Computational Framework Derived From Public Data — No Institutional Access Required

**Author:** Eric Robert Lawson  
**Affiliation:** Independent mathematician  
**Date:** 2026-03-03  
**License:** MIT  
**Languages:** Python, TeX  

---

## What This Repository Is

This repository contains a systematic computational analysis of **18+ cancer types** using a unified mathematical framework derived from the Waddington epigenetic landscape.

Every cancer type in this repository was analyzed using the same reproducible protocol applied to public gene expression data (TCGA and GEO). Each analysis produces a **depth score** — a continuous variable that measures how far a tumour has traveled from its normal cellular identity toward a false attractor state — and a **depth-stratified drug map** that specifies which molecular targets are relevant at which position in the landscape.

This is not a bioinformatics pipeline. It is a **geometric framework** for understanding cancer as a position problem rather than a mutation problem.

---

## The Core Claim

Cancer is a cell that has been displaced from its normal attractor in the Waddington epigenetic landscape and trapped in a false attractor that supports malignant survival.

The geometry of that displacement is **measurable** from bulk RNA-seq data.  
The depth of commitment to the false attractor is **rankable** across patients in the same cohort.  
The molecular dependencies of each depth position are **derivable** from that ranking alone.  
The drugs that target those dependencies **converge with published clinical trial targets** in every cancer analyzed.

**Zero false positives in direction across 18+ independent cancer analyses.**

That is not a claim about the strength of individual correlations. It is a claim about the reliability of the method: every geometry-derived drug target has been confirmed by published pharmacology or active clinical trials, without knowing the pharmacology in advance.

---

## The 18+ Cancer Types

| Cancer | Full Name | Cell of Origin | Dataset | Folder |
|--------|-----------|----------------|---------|--------|
| **ALL** | Acute Lymphoblastic Leukemia | B/T lymphoid progenitor | GEO | [Cancer_Research/ALL](Cancer_Research/ALL/) |
| **AML** | Acute Myeloid Leukemia | Myeloid progenitor | GEO / TCGA-LAML | [Cancer_Research/AML](Cancer_Research/AML/) |
| **BLCA** | Bladder Urothelial Carcinoma | Urothelial cell | TCGA-BLCA | [Cancer_Research/BLCA](Cancer_Research/BLCA/) |
| **BRCA** | Breast Invasive Carcinoma | Luminal / basal epithelial | TCGA-BRCA | [Cancer_Research/BRCA](Cancer_Research/BRCA/) |
| **CLL** | Chronic Lymphocytic Leukemia | Mature B cell | GEO | [Cancer_Research/CLL](Cancer_Research/CLL/) |
| **CML** | Chronic Myeloid Leukemia | Hematopoietic stem cell | GEO | [Cancer_Research/CML](Cancer_Research/CML/) |
| **CRC** | Colorectal Carcinoma | Intestinal epithelial cell | TCGA-COAD/READ | [Cancer_Research/CRC](Cancer_Research/CRC/) |
| **ESCA** | Esophageal Carcinoma | Esophageal squamous / Barrett | TCGA-ESCA | [Cancer_Research/ESCA](Cancer_Research/ESCA/) |
| **GBM** | Glioblastoma Multiforme | Neural progenitor / astrocyte | TCGA-GBM | [Cancer_Research/GBM](Cancer_Research/GBM/) |
| **HCC** | Hepatocellular Carcinoma | Hepatocyte | TCGA-LIHC / GEO | [Cancer_Research/HCC](Cancer_Research/HCC/) |
| **ICC** | Intrahepatic Cholangiocarcinoma | Cholangiocyte | GEO | [Cancer_Research/ICC](Cancer_Research/ICC/) |
| **LUAD** | Lung Adenocarcinoma | Alveolar type II cell | TCGA-LUAD | [Cancer_Research/LUAD](Cancer_Research/LUAD/) |
| **MDS** | Myelodysplastic Syndrome | Hematopoietic stem cell | GEO | [Cancer_Research/MDS](Cancer_Research/MDS/) |
| **MM** | Multiple Myeloma | Plasma cell | GEO | [Cancer_Research/MM](Cancer_Research/MM/) |
| **PAAD** | Pancreatic Adenocarcinoma | Acinar / ductal cell | TCGA-PAAD | [Cancer_Research/PAAD](Cancer_Research/PAAD/) |
| **PRAD** | Prostate Adenocarcinoma | Luminal epithelial cell | TCGA-PRAD | [Cancer_Research/PRAD](Cancer_Research/PRAD/) |
| **RCC** | Renal Cell Carcinoma (4 subtypes) | Proximal tubule / intercalated / collecting duct | TCGA-KIRC/KIRP/KICH + GEO | [Cancer_Research/RCC](Cancer_Research/RCC/) |
| **STAD** | Stomach Adenocarcinoma | Gastric pit / chief / parietal cell | TCGA-STAD | [Cancer_Research/STAD](Cancer_Research/STAD/) |

> RCC contains four independently analyzed subtypes: ccRCC, PRCC, chRCC, and cdRCC — each with their own depth score, drug map, and cross-type comparison.

---

## What Each Cancer Analysis Produces

Every cancer folder follows the same document structure derived from the [Workflow Protocol](Cancer_Research/Workflow_Protocol.md):

```
[Cancer]/
  Document [N]a   — Script 1: first contact with the data.
                    Analyst predictions locked before data loads.
                    Results verbatim. What was confirmed.
                    What was wrong and what it taught.

  Document [N]b   — Script 2: corrected framework.
                    Circuit integrity tests.
                    Drug target depth correlations.
                    Final attractor geometry.
                    Novel predictions listed and dated
                    before literature check.

  Document [N]c   — Literature check.
                    All predictions locked before any search.
                    Each finding classified:
                      CONFIRMED BY LITERATURE
                      EXTENDED BY LITERATURE
                      NOT IN LITERATURE
                      CONTRADICTED BY LITERATURE

  Document [N]d+  — Additional scripts where required.
                    Survival analysis.
                    Clinical panel validation.
                    Cross-subtype comparison (RCC only).

  Scripts         — Python scripts preserved exactly as run.
                    Not modified after execution.
                    Reproducible from GEO accession alone.
```

Each complete analysis ends with three clinical outputs:

1. **The depth score** — a continuous variable (0–1) measuring attractor commitment
2. **The 3-gene clinical panel** — measurable by standard IHC, r > 0.85 with full depth score
3. **The drug map** — which drug class targets which depth stratum, with contraindications

---

## The Protocol

The complete reproducible workflow is in:

**[Cancer_Research/Workflow_Protocol.md](Cancer_Research/Workflow_Protocol.md)**

It specifies five phases:

- **Phase 0** — Dataset discovery and quality assessment
- **Phase 1** — Biological grounding: predictions stated before data loads
- **Phase 2** — Script 1 (discovery run): first contact with data
- **Phase 3** — Script 2 (iteration run): corrected framework and drug map
- **Phase 4** — Literature check: geometry assessed against published biology
- **Phase 5** — README update and series carry-forward

**The order cannot be changed.** Predictions are locked before data. Data is analyzed before literature. Literature is consulted after analysis. This order is what makes the results valid and distinguishes this from post-hoc data mining.

---

## Getting Started: Running a New Cancer Analysis

**Requirements:**
- Python 3.8+
- `numpy`, `pandas`, `scipy`, `matplotlib`
- A GEO accession number for your cancer of interest
- No institutional access, no proprietary data, no cloud compute

**Start here:**

```bash
python Cancer_Research/universal_discovery_start_script.py
```

Change the `GEO_ACCESSION` variable to your target dataset. The script performs blind saddle point detection, depth score derivation, and landscape figure generation from raw GEO data.

Then follow [Workflow_Protocol.md](Cancer_Research/Workflow_Protocol.md) Phase 0 through Phase 5.

---

## The Waddington Foundation

The mathematical basis for this framework is the **Waddington epigenetic landscape**: a potential energy surface over gene expression space in which:

- Normal differentiated cells occupy **stable valleys** (attractors)
- Cancer occurs when a cell crosses a **saddle point** and falls into a **false attractor valley**
- The depth of the false valley is continuously measurable from gene expression data
- The geometry of the landscape determines therapeutic vulnerability, not just mutation status

This is not a metaphor. It is a measurable geometric structure recoverable from bulk RNA-seq data.

The mathematical objects at the centre of each analysis:
- **Depth score**: projection of tumour transcriptome onto the normal→false attractor axis
- **Rank correlation**: Spearman r between gene expression and depth score across all tumours
- **Saddle point**: the transcriptomic transition state between normal and false attractor valleys
- **Transition Index (TI)**: a 2-gene clinical approximation of the full depth axis

For the full theoretical foundation:  
→ [Cancer_Research/OrganismCore_Cancer_Framework.md](Cancer_Research/OrganismCore_Cancer_Framework.md)  
→ [Cancer_Research/waddington_saddle_point_cancer_reversion.md](Cancer_Research/waddington_saddle_point_cancer_reversion.md)

---

## Cross-Cancer Patterns (Empirical, Not Assumed)

Across 18+ independent analyses, the following patterns emerged. These are **empirical findings**, not axioms:

| Pattern | Description |
|---------|-------------|
| **Universal false attractor** | Every cancer analyzed has a measurable false attractor state. The depth score always separates tumour from normal and correlates with adverse prognosis. |
| **Drug target derivation** | In every cancer, the geometry-derived drug target was confirmed by published pharmacology or active clinical trials. |
| **3-gene clinical panel** | In every cancer, a 3-gene subset of the depth score achieves r > 0.85 with the full score and is measurable by standard IHC. |
| **Normal identity displacement** | The normal cell's identity programme (transporters, metabolic enzymes, differentiation TFs) is universally lost as depth increases, regardless of cancer type. |
| **Epigenetic lock** | Every false attractor is maintained by an epigenetic lock mechanism. EZH2 is the lock in the majority of solid tumours analyzed. Its direction (gain or loss) must be determined from data — it is not universal in sign. |
| **Circuit integrity varies** | In some cancers, the master differentiation circuit is intact and restorable. In others, it is uncoupled. This determines whether circuit restoration or attractor dissolution is the correct therapeutic strategy. |

---

## The RCC Multi-Subtype Case

The Renal Cell Carcinoma series is the most complete analysis in the repository. It demonstrates what the framework produces when applied to **four molecularly distinct subtypes from the same organ**:

- **ccRCC** (clear cell) — proximal tubule, VHL loss, HIF lock, n=534
- **PRCC** (papillary) — proximal tubule different saddle, MET/FH, two-phase Waddington crossing, n=290
- **chRCC** (chromophobe) — intercalated cell, chromosomal losses, DNMT writer switch, PC2 axis, n=150
- **cdRCC** (collecting duct) — principal cell, TCA collapse, rarest subtype, n=7

A three-script cross-type analysis compared the four independently derived landscapes and found:
- **1 universal marker** (LOXL2, 4/4 types)
- **1 universal circuit** (TCA → αKG → EZH2, present in 3/4 types)
- **1 structural inversion** (chRCC uses DNMT3B instead of DNMT3A — a writer switch, not a lock — invisible without cross-type comparison)
- **25 novel predictions** not previously assembled in the published literature
- **3 immune architecture types** across the four subtypes

The cross-type analysis is only possible because the four individual analyses were conducted independently with no cross-type context. The comparison is a ranked-list comparison and circuit topology comparison — not a pooled expression analysis.

→ [Cancer_Research/RCC/](Cancer_Research/RCC/)

---

## Epistemic Protocol

The credibility of this framework rests on one discipline:

**Predictions are stated before data is analyzed. Data is analyzed before literature is checked. Literature is never consulted before predictions are locked.**

This order is enforced in every analysis in this repository. Wrong predictions are documented alongside correct predictions. Analyst assumption errors are recorded explicitly — they are not failures, they are data about where prior knowledge was incomplete.

A repository that hides its errors cannot be trusted.  
This repository shows every error alongside every confirmation.

---

## For Computational Oncologists

If you work in computational oncology and want to evaluate this framework:

1. **Check the convergences first.** In every cancer folder, the literature check document lists what the geometry found against what the published literature established independently. Every drug target has a published confirmation. The convergences are not cherry-picked — they include everything the geometry found.

2. **Check the wrong predictions.** The analyst assumption errors in each [N]a document are the honest record. A framework that produces only confirmations is not being honest. The wrong predictions teach as much as the correct ones.

3. **Reproduce one analysis.** Pick a cancer type, download the GEO data, run the scripts. You should get the same numbers. If you do not, open an issue stating what differed. Difference is information.

4. **Extend to a new cancer.** The [Workflow Protocol](Cancer_Research/Workflow_Protocol.md) is self-contained. A competent bioinformatician following the protocol should be able to produce an equivalent analysis for any cancer with a suitable GEO dataset.

**If you are at a research institution and want to discuss collaboration, validation, or clinical translation of any specific cancer analysis: open an issue or contact via GitHub.**

---

## Repository Structure

```
attractor-oncology/
│
├── README.md                        ← You are here
├── LICENSE                          ← MIT
│
├── Cancer_Research/
│   ├── OrganismCore_Cancer_Framework.md     ← Full framework description
│   ├── Workflow_Protocol.md                 ← Reproducible analysis protocol (v2.0)
│   ├── universal_discovery_start_script.py  ← Entry point for new cancer analysis
│   ├── waddington_saddle_point_cancer_reversion.md
│   │
│   ├── ALL/    Acute Lymphoblastic Leukemia
│   ├── AML/    Acute Myeloid Leukemia
│   ├── BLCA/   Bladder Carcinoma
│   ├── BRCA/   Breast Carcinoma
│   ├── CLL/    Chronic Lymphocytic Leukemia
│   ├── CML/    Chronic Myeloid Leukemia
│   ├── CRC/    Colorectal Carcinoma
│   ├── ESCA/   Esophageal Carcinoma
│   ├── GBM/    Glioblastoma Multiforme
│   ├── HCC/    Hepatocellular Carcinoma
│   ├── ICC/    Intrahepatic Cholangiocarcinoma
│   ├── LUAD/   Lung Adenocarcinoma
│   ├── MDS/    Myelodysplastic Syndrome
│   ├── MM/     Multiple Myeloma
│   ├── PAAD/   Pancreatic Adenocarcinoma
│   ├── PRAD/   Prostate Adenocarcinoma
│   ├── RCC/    Renal Cell Carcinoma (ccRCC, PRCC, chRCC, cdRCC)
│   └── STAD/   Stomach Adenocarcinoma
│
├── Onboarding_1/   Framework orientation materials
└── Onboarding_2/   Protocol onboarding materials
```

---

## Citation

If you use this framework, reproduce an analysis, or build on any finding:

```
Eric Robert Lawson. Attractor Oncology: Waddington Landscape Geometry
Applied to Cancer Patient Stratification Across 18+ Cancer Types.
GitHub: https://github.com/Eric-Robert-Lawson/attractor-oncology, 2026.
[Include GEO accession of dataset used and document number of finding cited]
```

If you reproduce an analysis and your results differ from those recorded here: open an issue, state the document, state the discrepancy, state the evidence. Difference is information. Honest correction is part of the protocol.

---

## Origin

This work originated in a larger foundational repository, [OrganismCore](https://github.com/Eric-Robert-Lawson/OrganismCore), which applies the same attractor geometry framework across multiple domains beyond cancer. The cancer work was separated into this dedicated repository to make it findable and usable by the oncology research community.

The author is a mathematician, not a physician or biologist. The framework was derived from mathematical first principles and applied to biological data. The biological interpretation emerged from the geometry. The geometry did not emerge from the biology.

---

## Goal

One patient with stomach cancer goes to a clinic.  
A pathologist runs three IHC stains.  
The depth score is computed.  
The oncologist reads the drug map.  
The patient receives the drugs that target their geometry.  
Not the drugs that target the average patient who does not exist.

That is the goal.  
Every script, every document, every wrong prediction recorded — all of it serves that purpose.
