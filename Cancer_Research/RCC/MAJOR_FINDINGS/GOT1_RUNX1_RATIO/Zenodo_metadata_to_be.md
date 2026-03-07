# Zenodo Upload Metadata
## GOT1/RUNX1 Transition Index — ccRCC Depth Staging Tool

---

## TITLE

The GOT1/RUNX1 Transition Index: A Geometry-Derived 2-Gene Attractor Depth Staging Tool for Clear Cell Renal Cell Carcinoma with Overall Survival Validation in TCGA-KIRC (n=532)

---

## AUTHORS

**Eric Robert Lawson**
- Affiliation: OrganismCore (Independent Research)
- ORCID: 0009-0002-0414-6544
- Email: OrganismCore@proton.me

---

## DESCRIPTION (paste into Zenodo description field)

This document records the derivation, mechanistic basis, pre-specified predictions, computational survival validation, literature novelty confirmation, and IHC translation pathway for the GOT1/RUNX1 Transition Index (TI) as a 2-gene attractor depth staging tool for clear cell renal cell carcinoma (ccRCC).

**The Transition Index:**
TI = norm(GOT1) − norm(RUNX1)

GOT1 (aspartate aminotransferase / malate-aspartate shuttle anchor) falls continuously with attractor depth (r = −0.527) in TCGA-KIRC. RUNX1 (false attractor chromatin lock hub) rises with depth (r = +0.559). Their normalised difference captures the position of any tumour on the proximal tubule identity → false attractor axis.

**How it was derived:**
Not from survival data mining. Derived from Waddington landscape attractor geometry applied to TCGA-KIRC bulk expression data (n=534 tumour samples). The TI definition and the survival predictions were locked in a timestamped document before Script 4 (the survival analysis) was executed on 2026-03-07.

**Survival validation (pre-specified, TCGA-KIRC, n=532, events=175):**
- Cox regression (continuous −TI): HR = 6.94 [95% CI: 3.62–13.29], p = 5.09×10⁻⁹, C-index = 0.627
- KM median split: TI-High (GOT1-dominant) median OS = not reached; TI-Low (RUNX1-dominant) median OS = 1,964 days. Log-rank p = 7.43×10⁻⁶
- Q4 vs Q1 depth quartile log-rank: p = 0.0001
- Q3+Q4 vs Q1+Q2 log-rank: p = 9.60×10⁻⁵
- 20/21 individual gene OS directions confirmed as predicted
- Contradictions: 0

**Literature novelty:**
RUNX1 as a standalone prognostic marker in ccRCC is confirmed in prior literature (Cancer Research 2020, PeerJ 2019, Immunity/Inflammation and Disease 2024/2025). GOT1 as a metabolic collapse marker in ccRCC is consistent with published cancer metabolism literature (Frontiers in Oncology 2024). The GOT1/RUNX1 ratio as a geometry-derived attractor depth index with pre-specified OS validation has no prior publication. Searched directly: no paper reports a combined GOT1/RUNX1 prognostic ratio in ccRCC.

**Additional novel finding (PRCC Type 2):**
The EZH2 paradox in PRCC Type 2 — EZH2 univariate HR = 1.85 (worse OS), EZH2 multivariate adjusting for MKI67 HR = 0.19 (better OS within MKI67-high) — is confirmed in TCGA-KIRP (n=142 Type 2, events=31). This is the same attractor depth paradox mechanism confirmed in breast cancer TNBC (OrganismCore CS-LIT-10, GSE25066). Cross-cancer structural confirmation of the attractor depth framework.

**Treatment selection implication:**
The TI assigns patients to quartile-matched drug logic. Q4 (TI-Low, RUNX1-dominant) patients should not receive anti-PD-L1 monotherapy (CD274 flat/falling in Q4) or belzutifan alone (HIF-2α not the dominant lock in deep tumours). The RUNX1-high belzutifan resistance prediction (RCC-N1, locked 2026-03-02) is directly testable from LITESPARK-005/013 trial data.

**IHC translation:**
The calibration study required before clinical deployment follows the same design as the FOXA1/EZH2 breast cancer protocol (OrganismCore CS-LIT-22, DOI: 10.5281/zenodo.18892788). n=200–400 archived FFPE ccRCC cases, RUNX1 IHC + GOT1 IHC, H-score blinded to RNA-TI and outcome, Youden J cut-point derivation against RNA-TI quartile ground truth. No IHC H-score cut-points currently exist.

**Honest limits:**
Single dataset (TCGA-KIRC). RNA-level only. All drug predictions are preclinical hypotheses. GOT1 IHC requires additional validation in FFPE renal tissue before the calibration study can begin.

**Repository:** https://github.com/Eric-Robert-Lawson/attractor-oncology

---

## KEYWORDS (enter each separately in Zenodo)

- clear cell renal cell carcinoma
- ccRCC
- GOT1
- RUNX1
- Transition Index
- attractor depth
- Waddington landscape
- biomarker
- overall survival
- TCGA-KIRC
- prognosis
- 2-gene index
- epigenetic
- chromatin lock
- proximal tubule identity
- belzutifan resistance
- EZH2 paradox
- PRCC
- geometry-derived
- OrganismCore

---

## LICENSE

Creative Commons Attribution 4.0 International (CC BY 4.0)

---

## RESOURCE TYPE

Publication → Preprint

---

## RELATED IDENTIFIERS (add in Zenodo "Related works" field)

| Relation | Identifier | Notes |
|---|---|---|
| Is supplement to | https://github.com/Eric-Robert-Lawson/attractor-oncology | OrganismCore repository |
| Is related to | https://doi.org/10.5281/zenodo.18892788 | CS-LIT-22: FOXA1/EZH2 IHC tool (same methodology, breast cancer) |
| Is related to | https://doi.org/10.5281/zenodo.18883922 | CS-LIT-1: FOXA1/EZH2 main BRCA validation paper |
| References | https://peerj.com/articles/7854/ | RUNX1 in ccRCC (Qiu et al., independent confirmation) |
| References | https://aacrjournals.org/cancerres/article/80/11/2325/640655/RUNX1-Is-a-Driver-of-Renal-Cell-Carcinoma | RUNX1 as ccRCC driver (Cancer Research 2020) |
| References | https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2024.1519046/full | GOT1 in cancer metabolism (Frontiers 2024) |

---

## VERSION NOTES

v1.0 (2026-03-07): Initial release. Survival validation complete.
IHC calibration study not yet initiated.

---

## SUGGESTED ZENODO COMMUNITY SUBMISSIONS

- Zenodo > Cancer Research
- Zenodo > Bioinformatics
- Zenodo > Kidney Cancer
- Zenodo > Biomarkers
