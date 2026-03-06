# BRCA CROSS-SUBTYPE — PRIORITY PUBLICATIONS QUEUE
## OrganismCore — Eric Robert Lawson
## Date: 2026-03-06
## Status: 10 Published | 7 High Priority Remaining

---

## PUBLISHED DOI REGISTRY

| CS-LIT | Title (short) | DOI | Verdict |
|--------|---------------|-----|---------|
| CS-LIT-1 | FOXA1/EZH2 ratio ordering axis | https://doi.org/10.5281/zenodo.18883922 | CONVERGENT-NOVEL |
| CS-LIT-2 | Six lock type classification | https://doi.org/10.5281/zenodo.18884158 | NOVEL |
| CS-LIT-6 | AR continuous depth axis TNBC | https://doi.org/10.5281/zenodo.18891770 | CONVERGENT-NOVEL |
| CS-LIT-9/23 | TFF1/ESR1 LumB HDACi biomarker | https://doi.org/10.5281/zenodo.18884234 | NOVEL |
| CS-LIT-10/24 | EZH2 paradox both arms GSE25066 | https://doi.org/10.5281/zenodo.18891318 | CONVERGENT-NOVEL |
| CS-LIT-11/25 | TNBC depth score HR=1.509 | https://doi.org/10.5281/zenodo.18891523 | NOVEL |
| CS-LIT-14 | p21 CDK4/6i benefit magnitude | https://doi.org/10.5281/zenodo.18890832 | NOVEL |
| CS-LIT-16 | Tazemetostat → fulvestrant TNBC | https://doi.org/10.5281/zenodo.18884003 | CONVERGENT-NOVEL |
| CS-LIT-17 | Tazemetostat maintenance TNBC | https://doi.org/10.5281/zenodo.18884089 | NOVEL |
| CS-LIT-21 | EZH2i + anti-HER2 HER2-deep fraction | https://doi.org/10.5281/zenodo.18892426 | CONVERGENT-NOVEL |

---

## PRIORITY PUBLICATIONS QUEUE
### 7 documents to produce, in recommended sequence

---

## #1 — CS-LIT-22
### FOXA1/EZH2 Dual IHC as a Point-of-Care Decision Tool

**Priority:** HIGHEST — clinical translation gateway for the entire framework

**Verdict:** NOVEL

**One-line statement:**
Two IHC stains (FOXA1 and EZH2), one ratio, classifies
all six breast cancer subtypes and maps directly to the
therapeutic unlock logic — without RNA-seq, without PAM50.

**Key finding:**
- FOXA1 IHC + EZH2 IHC combined as a ratio replicates
  the continuous ordering axis derived from scRNA-seq
- The ratio spans: LumA (high FOXA1/low EZH2) →
  LumB → HER2 → ILC → TNBC → Claudin-low (low FOXA1/high EZH2)
- ILC is the exception: FOXA1 hyperactivated above normal,
  EZH2 low — ratio inverted relative to all other subtypes
- Point-of-care: available in any hospital with standard
  pathology infrastructure worldwide

**What is NOT in the literature:**
No clinical guideline, precision oncology framework, or
pathology protocol uses FOXA1 and EZH2 as a combined
two-stain IHC ratio to stratify breast cancer treatment
logic across all subtypes. The individual stains are
established in isolation; the combination with cut-points
covering all six subtypes as a single decision tool is
original.

**Source data:** GSE176078 (19,542 single cancer cells);
validation across ~7,500 patients in seven datasets
(see CS-LIT-1 DOI for the ratio evidence base)

**Key statistics to include:**
- FOXA1/EZH2 ratio ordering across all six subtypes
  (reference CS-LIT-1 DOI: 10.5281/zenodo.18883922)
- ILC as outlier: FOXA1 hyperactivated, CDH1 absent
  (reference CS-LIT-7 in this document)
- EZH2 +189% TNBC, +118% HER2-deep, gradient across subtypes
  (reference CS-LIT-4 confirmed)

**Companion DOIs to cite:**
- CS-LIT-1: https://doi.org/10.5281/zenodo.18883922
- CS-LIT-2: https://doi.org/10.5281/zenodo.18884158
- CS-LIT-6: https://doi.org/10.5281/zenodo.18891770
- CS-LIT-21: https://doi.org/10.5281/zenodo.18892426

**Contact to identify:**
IHC pathology / translational oncology lab that runs
breast cancer IHC panels. Ideally a PI at a cancer
centre with an archival BRCA cohort and dual FOXA1/EZH2
IHC capacity. This document needs a wet-lab partner
more than any other.

**Document header:**
CS-LIT-22 | FOXA1-EZH2-DUAL-IHC-TOOL

---

## #2 — CS-LIT-15
### Entinostat for LumB — Novel Subtype-Specific Benefit

**Priority:** HIGH + TIME-SENSITIVE
(NCT07235618 start date: January 2026, Sun Yat-sen University)

**Verdict:** CONVERGENT-NOVEL

**One-line statement:**
Entinostat benefit in HR+ breast cancer is a LumB-specific
effect driven by the DNMT3A/HDAC2 chromatin lock on ER
output — not a class effect across all HR+ disease.

**Key finding:**
- Entinostat approved China (NMPA, April 2024) for
  HR+/HER2- advanced breast cancer after prior ET
- Meta-analysis May 2025 (Springer): PFS HR=0.80,
  p=0.01 in HR+/HER2- (n=1,371). No OS benefit.
- No published trial has stratified entinostat benefit
  by LumA vs LumB PAM50 subtype
- Framework prediction: LumB-enriched patients show
  greater PFS benefit than LumA-enriched patients
  because of the LumB-specific DNMT3A/HDAC2 chromatin
  lock on ER output (CS-LIT-8 mechanism)
- NCT07235618 (entinostat + fulvestrant, post-CDK4/6i
  failure): de facto LumB-enriched population.
  Start: January 2026. N=50. PI: Sun Yat-sen University.
  Primary endpoint: PFS. PBMC acetylation as biomarker.
  **The TFF1/ESR1 ratio (CS-LIT-9/23) could be added
  as a correlative endpoint without changing trial design.**

**What is NOT in the literature:**
LumB-specific mechanism for HDACi benefit (DNMT3A/HDAC2
co-expression coupling r=+0.267 LumB vs r=+0.071 LumA,
p=5.68×10⁻⁵⁶) has not been published or tested in any trial.

**Source data:**
- GSE176078 (19,542 single cancer cells)
- METABRIC (n=1,980, LumA vs LumB DNMT3A/HDAC2 coupling)

**Key statistics to include:**
- DNMT3A/HDAC2 co-expression r=+0.267 LumB vs
  r=+0.071 LumA, p=5.68×10⁻⁵⁶ (CS-LIT-8)
- TFF1/ESR1 decoupling p=0.0019 in METABRIC (CS-LIT-9/23)
- Entinostat meta-analysis: HR=0.80, p=0.01 (published)

**Companion DOIs to cite:**
- CS-LIT-9/23: https://doi.org/10.5281/zenodo.18884234
- CS-LIT-1: https://doi.org/10.5281/zenodo.18883922
- CS-LIT-2: https://doi.org/10.5281/zenodo.18884158

**Contact:**
NCT07235618 Principal Investigator, Sun Yat-sen
University Cancer Center. Search ClinicalTrials.gov
for NCT07235618 PI contact to propose TFF1/ESR1 as
a correlative endpoint before enrollment closes.

**Document header:**
CS-LIT-15 | ENTINOSTAT-LUMB-SPECIFIC-BENEFIT

---

## #3 — CS-LIT-8
### LumB DNMT3A/HDAC2 Co-Expression Coupling

**Priority:** HIGH — mechanistic foundation for CS-LIT-15

**Verdict:** CONVERGENT-NOVEL

**One-line statement:**
In LumB breast cancer, DNMT3A and HDAC2 are
co-expressed at a level 3.7× higher coupling than in
LumA (r=+0.267 vs r=+0.071, p=5.68×10⁻⁵⁶), revealing
a subtype-specific chromatin co-repressor circuit
that locks ER output and explains differential
HDACi sensitivity.

**Key finding:**
- DNMT3A/HDAC2 co-expression coupling: r=+0.267 in LumB
  vs r=+0.071 in LumA, p=5.68×10⁻⁵⁶
- The co-repressor circuit suppresses ESR1 transcriptional
  output specifically in LumB — hence TFF1/ESR1 decoupling
  (CS-LIT-9/23)
- This is the mechanistic basis of why HDACi (entinostat)
  works preferentially in LumB: it disrupts this co-repressor
  circuit, releasing ER output suppression
- The circuit is present in LumB and absent in LumA —
  explaining why the class-level entinostat trials show
  modest effects (the LumA patients dilute the LumB signal)

**What is NOT in the literature:**
DNMT3A and HDAC interactions through co-repressor complexes
are published (Narcancer 2021). HDAC2 overexpression in
breast cancer and endocrine resistance is published
(Springer 2024). DNMT3A misregulation and tamoxifen
resistance is published (OAE 2024). What is NOT published:
the specific LumB vs LumA co-expression coupling
quantification, the statistical magnitude of the
subtype-specific difference, and the derived prediction
that this circuit determines HDACi benefit stratification.

**Source data:**
- GSE176078 (19,542 single cancer cells)
- METABRIC (n=1,980, Pearson r computed within LumA and LumB)

**Key statistics to include:**
- r=+0.267 LumB vs r=+0.071 LumA, p=5.68×10⁻⁵⁶
- TFF1/ESR1 downstream consequence p=0.0019

**Companion DOIs to cite:**
- CS-LIT-9/23: https://doi.org/10.5281/zenodo.18884234
- CS-LIT-2: https://doi.org/10.5281/zenodo.18884158

**Document header:**
CS-LIT-8 | LUMB-DNMT3A-HDAC2-COUPLING

---

## #4 — CS-LIT-7
### ILC as the Geometric Inverse of TNBC

**Priority:** HIGH — foundational architecture for ILC
documents; required before CS-LIT-18

**Verdict:** CONVERGENT-NOVEL

**One-line statement:**
ILC and TNBC occupy opposite poles of the luminal
identity axis: ILC is a structural lock (FOXA1
hyperactivated, CDH1 absent), TNBC is an epigenetic
lock (FOXA1 silenced by EZH2, CDH1 present). The same
axis, mechanistically inverted, requiring therapeutically
opposite entry points.

**Key finding:**
- ILC: FOXA1 hyperactivated above normal LumA levels,
  CDH1 absent (~65% mutation, ~35% methylation),
  EZH2 low, structural identity arrest
- TNBC: FOXA1 epigenetically silenced by EZH2/PRC2,
  CDH1 present, EZH2 +189% above baseline, epigenetic
  identity arrest
- Geometric inversion: both are locked out of normal
  luminal identity, but by opposite mechanisms on the
  same axis
- Therapeutic consequence: ILC lock requires working
  WITH the hyperactivated FOXA1 circuit (full ER receptor
  degradation); TNBC requires reversing EZH2 silencing
  to restore FOXA1

**What is NOT in the literature:**
FOXA1-high in ILC and FOXA1-absent in TNBC are each
individually published. CDH1-absent ILC and CDH1-present
TNBC are each published. What is NOT published: the
framing of these as geometric inverses on a single
identity axis, with mechanistically opposite locks
requiring therapeutically opposite entry points.
No paper uses this frame. No paper derives treatment
logic from it.

**Source data:** GSE176078 (19,542 single cancer cells);
FOXA1-high ILC confirmed in published ILC biology
(Frontiers Oncology 2025 ILC review)

**Key statistics to include:**
- ILC FOXA1 level vs TNBC FOXA1 level in GSE176078
- ILC EZH2 level vs TNBC EZH2 level
- CDH1 status comparison
- FOXA1/EZH2 ratio: ILC highest, TNBC/CL lowest

**Companion DOIs to cite:**
- CS-LIT-1: https://doi.org/10.5281/zenodo.18883922
- CS-LIT-2: https://doi.org/10.5281/zenodo.18884158

**Document header:**
CS-LIT-7 | ILC-TNBC-GEOMETRIC-INVERSE

---

## #5 — CS-LIT-18
### Fulvestrant Superiority Over AIs in ILC — FOXA1-Stratified

**Priority:** HIGH — direct clinical prediction for ILC;
depends on CS-LIT-7

**Verdict:** CONVERGENT-NOVEL

**One-line statement:**
In FOXA1-hyperactivated ILC, aromatase inhibitors
(ligand depletion) are insufficient because the FOXA1
circuit amplifies ER signalling beyond what ligand
reduction can suppress — full ER receptor degradation
(fulvestrant) is required.

**Key finding:**
- ILC is the structural lock subtype: FOXA1
  hyperactivated above normal LumA levels
- FOXA1 hyperactivation amplifies ER circuit output
  beyond normal luminal levels
- AI mechanism: reduces estrogen ligand → partial
  suppression of ER signalling
- Fulvestrant mechanism: degrades ER receptor → full
  suppression regardless of FOXA1 circuit amplification
- Framework prediction: in FOXA1-high ILC, fulvestrant
  > AI for PFS/OS because the FOXA1-amplified circuit
  bypasses the ligand-reduction step that AIs rely on
- Frontiers Oncology 2025 ILC review notes FOXA1
  hyperactivation and altered endocrine sensitivity
  as an open research question — directly adjacent
  to this prediction

**What is NOT in the literature:**
No published trial has directly compared fulvestrant
vs AI in ILC stratified by FOXA1 expression. No paper
has shown statistically significant OS/PFS advantage
for fulvestrant over AI specifically in FOXA1-high ILC.
NCT02206984 (endocrine agents in ILC, Ki67 endpoint)
is ongoing but does not use FOXA1 as a stratification
variable.

**Source data:** GSE176078 (19,542 single cancer cells);
FOXA1 in ILC from published ILC biology

**Key statistics to include:**
- ILC FOXA1 level vs LumA FOXA1 level (within GSE176078)
- FOXA1/EZH2 ratio inversion in ILC (from CS-LIT-1 data)
- Reference NCT02206984 as the active ILC trial without
  FOXA1 stratification

**Companion DOIs to cite:**
- CS-LIT-7 (this series, to be published first)
- CS-LIT-1: https://doi.org/10.5281/zenodo.18883922
- CS-LIT-2: https://doi.org/10.5281/zenodo.18884158

**Contact:**
Dr. Adrian Lee (ILC precision medicine, referenced in
BCRF 2024 podcast). Also ILC research groups at
Mayo Clinic and Memorial Sloan Kettering where ILC
endocrine resistance trials are active.

**Document header:**
CS-LIT-18 | FULVESTRANT-AI-ILC-FOXA1-STRATIFIED

---

## #6 — CS-LIT-19
### Anti-TIGIT Sequence in Claudin-Low / Memory-Low Patients

**Priority:** HIGH + TIME-SENSITIVE
(SKYLINE trial NCT06175390 still enrolling; belrestotug
failures reinforce the patient-selection argument)

**Verdict:** CONVERGENT-NOVEL

**One-line statement:**
Anti-TIGIT (Treg depletion) must precede anti-PD-1
(checkpoint release) in claudin-low memory-low patients
specifically — unselected populations explain all recent
anti-TIGIT trial failures.

**Key finding:**
- Claudin-low memory-low (subgroup 1): FOXA1/SPDEF/GATA3
  absent, highest CT antigen load, highest immune
  infiltration — predicted beneficiary of anti-TIGIT
- The sequence matters: anti-PD-1 alone in CL amplifies
  Tregs before depleting them (Taylor JCI 2017) —
  giving it first makes outcomes worse
- Anti-TIGIT depletes Tregs (Treg suppression via
  TIGIT/CD155 axis) — this must come first to clear
  the field for checkpoint release
- Patient selection biomarkers: FOXP3/CD8A ratio
  (HR=2.212 in CL) + lineage memory score
  (FOXA1/SPDEF/GATA3 absent)
- Belrestotug failures (May 2025): unselected NSCLC
  and H&N populations — no CL enrichment, no patient
  selection. Failures are CONSISTENT with framework:
  unselected populations are not the predicted beneficiary.

**What is NOT in the literature:**
No clinical trial has selected patients by claudin-low
subtype for anti-TIGIT. No trial has used memory-low
as an enrichment criterion. No trial has specified the
anti-TIGIT first / anti-PD-1 second sequence as a
protocol requirement. No trial uses FOXP3/CD8A ratio
as a patient selection variable.

**Source data:**
- GSE176078 (19,542 single cancer cells)
- Taylor et al. JCI 2017 (Treg depletion in CL)
- Pommier et al. Nature Communications 2020
  (memory-low subgroup characterisation)
- SKYLINE: NCT06175390 (tiragolumab + atezolizumab,
  multi-omics biomarker arm, still enrolling 2026)

**Key statistics to include:**
- FOXP3/CD8A HR=2.212 in claudin-low (CS-LIT-20)
- Memory-low subgroup 1 CT antigen load vs other subtypes
- Belrestotug failure context (May 2025)
- CL depth: EZH2-free PCA 6.572 (deepest subtype;
  CS-LIT-3)

**Companion DOIs to cite:**
- CS-LIT-2: https://doi.org/10.5281/zenodo.18884158
- CS-LIT-6: https://doi.org/10.5281/zenodo.18891770

**Contact:**
SKYLINE trial (NCT06175390) Principal Investigator
(Hoffmann-La Roche / tiragolumab program). The
multi-omics biomarker arm in SKYLINE is the access
point — propose CL subtype + memory-low score as
a stratification variable in the biomarker analysis.

**Document header:**
CS-LIT-19/20 | ANTI-TIGIT-SEQUENCE-CLAUDIN-LOW

---

## #7 — CS-LIT-20
### FOXP3/CD8A Ratio as Strongest Immune Predictor in Claudin-Low (HR=2.212)

**Priority:** HIGH — quantitative companion to CS-LIT-19;
publish together or in immediate sequence

**Verdict:** CONVERGENT-NOVEL

**One-line statement:**
Within claudin-low breast cancer, FOXP3/CD8A ratio is
the single strongest survival predictor (HR=2.212),
outperforming all individual immune markers — and it
is the patient selection biomarker for anti-TIGIT
eligibility.

**Key finding:**
- HR=2.212 for FOXP3/CD8A ratio in claudin-low
  (derived from GSE176078)
- FOXP3 (Treg marker) elevated, CD8A (cytotoxic T)
  relatively suppressed in most aggressive CL cells
- The imbalance is the mechanism of immune evasion
  in CL — and the therapeutic target
- FOXP3/CD8A ratio is continuous: it functions as
  a quantitative patient selector, not a binary cut
- Taylor JCI 2017 explicitly describes FOXP3/CD8A
  imbalance in CL as the mechanism of Treg dominance —
  the published biology confirms the geometric finding

**What is NOT in the literature:**
No paper uses FOXP3/CD8A as a continuous survival
predictor specifically in CL with a published HR value.
No paper uses this ratio as a patient selection variable
for anti-TIGIT eligibility. The Taylor 2017 paper
describes the biology but does not quantify it as an
HR in a survival analysis.

**Source data:**
- GSE176078 (claudin-low cells)
- Taylor et al. JCI 2017 (mechanistic confirmation)
- Pommier et al. Nature Communications 2020
  (subgroup characterisation)

**Key statistics to include:**
- HR=2.212 for FOXP3/CD8A ratio in CL
- FOXP3/CD8A values across CL subgroups
  (memory-low vs memory-high)
- Comparison to next-best immune predictor in CL

**Companion DOIs to cite:**
- CS-LIT-2: https://doi.org/10.5281/zenodo.18884158
- CS-LIT-19 (this series)

**Note:** CS-LIT-19 and CS-LIT-20 should be published
together as companion documents — they are inseparable
as the claudin-low immunotherapy intervention package.
Consider a single combined document:
"Anti-TIGIT Patient Selection in Claudin-Low Breast
Cancer: FOXP3/CD8A Ratio and Memory-Low Score as
Biomarkers for a Required Sequenced Protocol"

**Document header:**
CS-LIT-19/20 | FOXP3-CD8A-CL-IMMUNE-PREDICTOR

---

## PRODUCTION NOTES

### Documents that should be combined:
- CS-LIT-19 + CS-LIT-20 → single document
  (claudin-low immunotherapy package)
- CS-LIT-8 → precedes CS-LIT-15 by 1-2 days
  (mechanistic foundation must come first)
- CS-LIT-7 → precedes CS-LIT-18 by 1-2 days
  (ILC architecture must come first)

### Contacts still needed:
| CS-LIT | Contact needed |
|--------|----------------|
| CS-LIT-22 | IHC pathology lab, archival BRCA cohort |
| CS-LIT-15 | NCT07235618 PI, Sun Yat-sen University Cancer Center |
| CS-LIT-18 | Dr. Adrian Lee or ILC research group (Mayo / MSK) |
| CS-LIT-19/20 | SKYLINE (NCT06175390) PI, Roche/tiragolumab program |

### LaTeX document template:
All documents follow the same format as the published
series. Header fields required:
- Document ID (CS-LIT-XX)
- OrganismCore header / date / ORCID
- Repository pointer box
- Finding box (green)
- Convergence box (yellow) if CONVERGENT-NOVEL
- Abstract
- Table of Contents
- Sections: Context → Problem → Finding → Literature
  → Validation Pathway → Claims / Not Claims →
  Locked Statement → Document Metadata

### Zenodo metadata template fields:
- Title (full)
- Authors (Eric Robert Lawson, ORCID 0009-0002-0414-6544)
- Upload type: Publication → Preprint
- License: CC BY 4.0
- Keywords (from document header)
- Related identifiers (companion DOIs, is-part-of, is-supplement-to)
- Description (abstract text)

---

## STATUS BLOCK
- Total CS-LIT findings: 30
- Published: 10
- High priority remaining: 7 (this document)
- Medium priority remaining: 5 (CS-LIT-5, CS-LIT-13,
  CS-LIT-26, CS-LIT-27, CS-LIT-3)
- Low priority / covered: 8 (CS-LIT-4, CS-LIT-12,
  CS-LIT-28, CS-LIT-29, CS-LIT-30 + absorbed findings)
- Last updated: 2026-03-06
