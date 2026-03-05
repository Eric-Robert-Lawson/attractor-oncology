# CLAUDIN-LOW — SCRIPT 2 BEFORE-DOCUMENT
## Predictions Locked Before Script 2 Runs
## OrganismCore — Document BRCA-S7c
## Date: 2026-03-05

---

## DOCUMENT METADATA

```
document_id:        BRCA-S7c
series:             BRCA Deep Dive — Claudin-Low
folder:             Cancer_Research/BRCA/DEEP_DIVE/Claudin_Low/
type:               BEFORE-DOCUMENT (Script 2)
                    All predictions stated before
                    Script 2 runs. This document
                    cannot be modified after Script 2
                    begins.
date:               2026-03-05
author:             Eric Robert Lawson
                    OrganismCore
protocol_version:   Workflow_Protocol.md v2.0
attractor_type:     TYPE 4 — ROOT LOCK
                    (ATTRACTOR_GEOMETRY_AXIOMS.md v2.0)
based_on:           BRCA-S7b (script1_results_and_reasoning.md)
status:             LOCKED — predictions cannot change
                    after this document is committed
```

---

## PART I — WHAT SCRIPT 1 ESTABLISHED
### Confirmed facts entering Script 2. Not predictions.

**Geometry confirmed:**

1. Claudin-low identified from 10-gene geometry signature
   (score ≥ 7/10, ERBB2 exclusion): n=268 in TCGA-BRCA.
   Mixed purity — includes LumA-PAM50 (n=110), Basal (n=40),
   Normal-like (n=39), LumB (n=18), unassigned (n=58).
   Pure claudin-low estimated at ~40–60 samples.
   The depth axis structure is coherent even in the mixed set.

2. **Claudin loss is the defining structural fingerprint.**
   CLDN3, CLDN4, CLDN7 are all lower in the geometry-identified
   claudin-low set than in every other BRCA subtype.
   Depth correlations: r(CLDN3, depth)=-0.637, r(CLDN4,
   depth)=-0.641, r(CLDN7, depth)=-0.619. These are the
   strongest depth correlates in the analysis.

3. **The depth axis is coherent and bidirectional.**
   As depth increases: claudins fall (r≈-0.62 to -0.64),
   CDH1 falls (r=-0.459), luminal TFs fall (ESR1 r=-0.393,
   FOXA1 r=-0.430, GATA3 r=-0.507), and mesenchymal markers
   rise (VIM r=+0.609, SNAI1 r=+0.390, FN1 r=+0.351).

4. **FN1 (+26.4%, p=2.34e-49) and CDH2 (+48.2%, p=3.31e-29)**
   are the strongest confirmed mesenchymal/partial-EMT signals.
   The E→N cadherin switch is present and significant.

5. **Immune programme is the most robustly confirmed finding.**
   FOXP3 +66.7% (p=1.64e-40), PDCD1 +52.5% (p=8.59e-20),
   TIGIT +51.6% (p=6.86e-28), LAG3 +30.6% (p=9.90e-17),
   IFNG +67.8% (p=3.12e-06). Large, exhausted TIL population
   confirmed.

6. **Cancer-testis antigen de-repression in the unfiltered scan.**
   GAGE1, GAGE2D, GAGE4, GAGE12D, GAGE12J, CT45A3, CT45A4,
   STRA8, DPPA2 elevated vs. normal. Not predicted — emerged
   from data. Interpreted per Axioms v2.0 (Lesson 21) as a
   TYPE 4 epigenetic disinhibition signature.

7. **Proliferation is elevated** — MKI67 +41.0% (p=9.21e-43),
   TOP2A +44.3% (p=8.58e-44), CCNB1 +26.7% (p=7.13e-45).
   Cross-subtype rank: CL (10.51) between LumA (10.17) and
   LumB (11.37). Below TNBC (12.03).

8. **PARP1 elevated** +8.6% (p=1.53e-43). ERBB2 not amplified
   (+4.0% vs. normal, not at HER2-amplified levels). ESR1 near
   flat vs. normal in the full mixed set (depth-dependent in
   the deep subset).

9. **CLDN3 within claudin-low correlates positively with luminal
   identity markers:** r(CLDN3, FOXA1)=+0.326 (p=4.80e-08),
   r(CLDN3, GATA3)=+0.367 (p=6.04e-10). Shallow-end samples
   retain both claudins and partial luminal identity
   simultaneously.

10. **ALDH1A1 is reduced** vs. normal (-19.5%, p=1.29e-38) and
    does not correlate with the depth score (r=-0.025, ns).
    The bulk RNA-seq stem marker axis for claudin-low is
    VIM/FN1/SNAI1-driven, not ALDH1A1-driven.

---

## PART II — THE CLASSIFICATION CONTEXT FOR SCRIPT 2

The classifier impurity (268 samples, mixed PAM50 distribution)
is a structural finding, not an error. Script 2 addresses it
directly through stratified analysis.

**Three analysis strata for Script 2:**

```
STRATUM A — Full geometry set (n=268)
  All samples scoring ≥7 on the 10-gene signature
  after ERBB2 exclusion.
  The broadest stem/mesenchymal attractor space.
  Survival analysis here tests whether the
  geometry-derived depth score predicts outcomes
  regardless of PAM50 label.

STRATUM B — ESR1-low subset (n=TBD)
  Geometry set filtered to ESR1 below cohort median.
  Removes the LumA-PAM50 samples that contaminate
  the full set.
  Estimated n=80–120.
  Survival analysis here tests the purer
  claudin-low/basal-like attractor.

STRATUM C — Basal+Normal-PAM50 within geometry set
  The 40 Basal + 39 Normal-PAM50 claudin-low
  geometry samples (n=79).
  Closest to canonical published claudin-low.
  Small n — survival analysis may be underpowered
  but should be attempted.
```

All three strata should be tested. The depth score
prediction (S2-P1) should hold across all three.
If it holds in Stratum A but not Stratum C, the depth
axis is driven by the LumA impurity. If it holds in
Stratum C, the depth axis is genuinely claudin-low.

---

## PART III — WHAT SCRIPT 2 MUST DO

```
REQUIRED ANALYSES:

  1. Survival analysis: depth score tertiles vs OS
     (all three strata — primary test)

  2. Survival analysis: CLDN3 tertiles vs OS
     within claudin-low (claudin loss as survival axis)

  3. Survival analysis: immune score tertiles vs OS
     (FOXP3 + PDCD1 + TIGIT + LAG3 composite)

  4. Survival analysis: MKI67 tertiles vs OS
     (proliferation as survival axis — comparison
     to TNBC where MKI67 is a strong predictor)

  5. Cross-subtype survival: claudin-low vs LumA
     vs TNBC vs LumB vs HER2 vs Normal
     (all five subtypes in the same Kaplan-Meier plot)

  6. Stratified analysis: depth score in Stratum B
     (ESR1-low subset) vs OS — purity test

  7. PARP1-high vs PARP1-low survival within
     claudin-low (PARP inhibitor geometry test)

  8. CD274 tertile vs OS within claudin-low
     (PD-L1 as biomarker)

REQUIRED DATASETS:
  TCGA-BRCA expression: already downloaded
  TCGA-BRCA clinical: already downloaded
  TCGA pancan survival supplement: required
    (Survival_SupplementalTable_S1_20171025_xena_sp
    — confirmed available from prior scripts)

NO EXTERNAL DATASETS REQUIRED.
All analyses use the expression + clinical + survival
files already present from Script 1.
```

---

## PART IV — PREDICTIONS
### All stated 2026-03-05 before Script 2 runs

---

### S2-P1 — DEPTH SCORE PREDICTS OS WITHIN CLAUDIN-LOW
**Direction:** Depth-high claudin-low has WORSE overall survival
than depth-low claudin-low

**Geometric basis:**
The depth score constructed in Script 1 (positive markers VIM,
CD44, SNAI1, ZEB1, FN1 minus negative markers CLDN3, CLDN4,
CLDN7, CDH1, ESR1) measures how completely the cell has
abandoned all differentiation commitment. Deeper samples have
lost both claudin identity and luminal identity simultaneously.
This corresponds to the most undifferentiated, most invasive,
most immune-exhausted subset of the claudin-low population.

The structural prediction is: deeper root lock = worse clinical
behaviour. This is the Type 4 analogue of the depth score
predictions confirmed in TNBC (script2_results_and_reasoning.md
BRCA-S4d: depth predicts DRFS, log-rank p<0.0001) and tested
in ILC (S2-P2).

**Predicted result:**
Log-rank p < 0.05 for depth-high vs depth-low tertiles within
claudin-low (Stratum A). HR > 1.0 for depth-high (worse OS).
Effect should hold in Stratum B (ESR1-low subset) if the depth
axis is genuine claudin-low biology and not impurity-driven.

**Confidence: HIGH**

**Falsification criterion:**
S2-P1 is NOT confirmed if log-rank p > 0.10 in Stratum A,
OR if the depth score predicts in the wrong direction
(depth-high doing better). A null result in Stratum A but
a significant result in Stratum C would be partial confirmation
(impurity-dependent depth axis).

---

### S2-P2 — CLDN3 EXPRESSION PREDICTS OS — LOW CLDN3 = WORSE
**Direction:** CLDN3-low claudin-low has WORSE overall survival

**Geometric basis:**
CLDN3 is the strongest individual depth correlate
(r=-0.637). Within claudin-low, CLDN3 co-varies with
luminal identity retention (r(CLDN3, FOXA1)=+0.326,
r(CLDN3, GATA3)=+0.367). CLDN3-high samples are the
shallowest — retaining partial claudin and partial luminal
identity. CLDN3-low samples are the deepest — fully committed
to the root lock false attractor.

CLDN3 is a single-gene proxy for the full depth axis.
It should predict survival in the same direction as the
composite depth score, with slightly lower effect size.

This is a clinically relevant prediction: CLDN3 is a
measurable protein (IHC-testable). If CLDN3 expression
predicts OS within the claudin-low subset, it has immediate
clinical utility as a biomarker requiring no gene expression
profiling.

**Predicted result:**
Log-rank p < 0.10 for CLDN3-high vs CLDN3-low within
claudin-low. CLDN3-low = worse OS. HR > 1.0 for CLDN3-low.

**Confidence: MODERATE-HIGH**

**Falsification criterion:**
S2-P2 is NOT confirmed if log-rank p > 0.15 or HR is
in the wrong direction.

---

### S2-P3 — IMMUNE SCORE PREDICTS OS — HIGH IMMUNE = BETTER SURVIVAL
**Direction:** High immune infiltration (FOXP3/PDCD1/TIGIT/LAG3
composite) within claudin-low predicts BETTER overall survival

**Geometric basis:**
This is the most structurally nuanced prediction in this
document. The immune programme in claudin-low is large and
exhausted — FOXP3+, PDCD1+, TIGIT+, LAG3+ simultaneously.
This pattern has two possible survival interpretations:

  INTERPRETATION A (immune = better):
    High TIL burden in untreated tumours predicts better
    survival. The immune system is attempting to clear the
    tumour. Even in the absence of checkpoint therapy, high
    TIL = slower progression. This is established in TNBC
    literature (Adams et al., Denkert et al.): TIL count
    is a positive prognostic marker in TNBC.

  INTERPRETATION B (immune = neutral/worse due to exhaustion):
    The TILs are exhausted (PDCD1+, TIGIT+, LAG3+ together).
    Exhausted TILs do not kill effectively. High FOXP3
    indicates Tregs actively suppressing the response.
    The immune infiltrate may be present but non-functional.
    In this case, immune score would not predict better
    survival — or might predict worse (Treg-dominated
    immune microenvironment = immunosuppressive).

**The prediction is Interpretation A (better OS with high
immune)**, based on the TNBC precedent and the magnitude
of infiltration. The Treg/exhaustion pattern is acknowledged
as a confound.

**Predicted result:**
Log-rank p < 0.10 for immune-score-high vs immune-score-low
within claudin-low. High immune score = better OS or
longer progression-free interval.

**Confidence: MODERATE** — genuine mechanistic ambiguity
acknowledged above.

**Falsification criterion:**
S2-P3 is NOT confirmed if p > 0.15 in either direction,
OR if immune-high shows significantly WORSE survival
(which would indicate Treg-dominated immunosuppression
dominates over TIL benefit — a meaningful finding in its
own right).

---

### S2-P4 — MKI67 PREDICTS OS WITHIN CLAUDIN-LOW — HIGH MKI67 = WORSE
**Direction:** MKI67-high claudin-low has WORSE overall survival

**Geometric basis:**
MKI67 in claudin-low is +41.0% vs. normal — substantial
elevation. Cross-subtype rank: CL (10.51) between LumA
and LumB, below TNBC (12.03). Proliferation elevation in
undifferentiated/stem-like tumours is expected to be a
negative prognostic marker. Within claudin-low, the
samples with highest MKI67 are likely those with the
most proliferative stem-like behaviour.

In TNBC (BRCA-S4d), MKI67 correlated positively with
depth (r=+0.216) and depth predicted DRFS (p<0.0001).
The claudin-low prediction follows the same logic —
the most proliferative claudin-low tumours are the deepest
in the root lock false attractor and should have the
worst outcomes.

**Predicted result:**
Log-rank p < 0.10 for MKI67-high vs MKI67-low within
claudin-low. MKI67-high = worse OS. Effect size expected
to be smaller than TNBC MKI67 effect (because the
claudin-low set has lower average MKI67 than TNBC and
wider heterogeneity due to classifier impurity).

**Confidence: MODERATE**

**Falsification criterion:**
S2-P4 is NOT confirmed if p > 0.15 or HR in the wrong
direction. A null result here is expected to be more
likely than in TNBC due to the compressed MKI67
distribution in the mixed-purity claudin-low set.

---

### S2-P5 — CLAUDIN-LOW HAS WORSE OS THAN LumA, COMPARABLE TO TNBC
**Direction:** Claudin-low OS worse than LumA; not significantly
different from TNBC

**Geometric basis:**
The cross-subtype PCA showed claudin-low in an intermediate
position between LumA and TNBC (with impurity pulling it
toward LumA). The biology of the deep claudin-low subset
is more similar to TNBC: no ER, no HER2, high proliferation,
high immune infiltration. LumA has the best OS of all BRCA
subtypes; TNBC has among the worst.

The mixed-purity claudin-low set (containing LumA-PAM50
samples) may have intermediate OS between LumA and TNBC
as a consequence of the mixture. This is expected and
interpretable — it does not falsify the biology.

**Predicted result:**
Cross-subtype Kaplan-Meier: LumA > claudin-low ≥ TNBC
in OS (log-rank p < 0.05 for LumA vs claudin-low).
Claudin-low vs TNBC: not significantly different
(p > 0.10) — they occupy similar survival territory.

**Confidence: MODERATE** — mixture composition confounds
the claudin-low survival estimate.

**Falsification criterion:**
S2-P5 is NOT confirmed if claudin-low OS is equivalent
to LumA (p > 0.15). If claudin-low does as well as LumA
in TCGA, that reflects the LumA contamination dominating
the survival signal — interpretable but not the biology
of the true claudin-low population.

---

### S2-P6 — PARP1 ELEVATION DOES NOT PREDICT SURVIVAL WITHIN CLAUDIN-LOW
**Direction:** PARP1 level is NOT a significant survival
predictor within claudin-low (negative prediction)

**Geometric basis:**
PARP1 is elevated +8.6% in claudin-low vs normal
(p=1.53e-43). This is a cohort-level finding reflecting
elevated replication stress / proliferation. However,
PARP1 elevation in bulk RNA-seq is a consequence of
high proliferation across the sample set — it tracks
MKI67 (both are replication-associated). Within the
claudin-low set, PARP1 and MKI67 should co-vary.
PARP1 as an independent survival predictor (beyond
MKI67) is therefore not expected.

The PARP inhibitor drug target geometry is present
at the cohort level — but PARP1 mRNA level within
claudin-low is not expected to be a survival
biomarker for the reasons above.

**This is a negative prediction.** It distinguishes the
framework's drug target geometry conclusion (PARP1 is
elevated = PARP inhibitor geometry is present) from
the survival prediction (PARP1 level within the subtype
does not stratify outcomes — because it is confounded
by proliferation).

**Predicted result:**
Log-rank p > 0.10 for PARP1-high vs PARP1-low within
claudin-low. No significant survival difference
attributable to PARP1 independent of MKI67.

**Confidence: MODERATE (negative prediction)**

**Falsification criterion:**
S2-P6 is wrong (PARP1 DOES predict) if log-rank p < 0.05
in multivariate analysis controlling for MKI67. If PARP1
independently predicts OS beyond MKI67, that would suggest
BRCA1-pathway dysfunction is stratifying outcomes within
the claudin-low set — a finding requiring investigation.

---

### S2-P7 — DEPTH SCORE IN STRATUM B (ESR1-LOW SUBSET) PREDICTS OS
**Direction:** Depth score predicts OS in the ESR1-low
claudin-low subset with equal or greater effect size than
in the full set (Stratum A)

**Geometric basis:**
If the depth axis found in the full 268-sample set is a
genuine claudin-low biological signal (not an artefact of
LumA contamination), it should be stronger — not weaker —
in the purified ESR1-low subset. Removing the LumA-PAM50
samples from the analysis removes noise from the depth
score and should increase the signal-to-noise ratio of
the survival prediction.

If depth predicts OS in Stratum A but NOT in Stratum B,
the depth signal is LumA contamination driving the
result. That is an important methodological finding.

If depth predicts OS in Stratum B with comparable or
greater HR as in Stratum A, the depth axis is real
claudin-low biology.

**Predicted result:**
Log-rank p < 0.10 in Stratum B for depth-high vs depth-low.
HR in Stratum B ≥ HR in Stratum A. The depth signal
survives purification.

**Confidence: MODERATE** — Stratum B n will be smaller
(estimated 80–120) and power will be reduced. A null
result due to underpowering is not falsification.

**Falsification criterion:**
S2-P7 is NOT confirmed if HR in Stratum B is substantially
lower than in Stratum A (e.g., < 0.5× the Stratum A HR),
suggesting the depth signal was driven by the LumA
contamination rather than claudin-low biology.

---

### S2-P8 — CD274 (PD-L1) DOES NOT STRATIFY SURVIVAL
### WITHOUT CHECKPOINT THERAPY CONTEXT
**Direction:** CD274 expression level does NOT predict
OS within claudin-low in the TCGA cohort (negative prediction)

**Geometric basis:**
CD274 (PD-L1) is elevated in claudin-low (+6.3% vs normal,
p=0.038) but modestly. The correlation with depth score is
near zero (r=+0.002). Within the claudin-low set, CD274
varies independently of the depth structure — it is not
a depth marker; it is a TIL-response marker.

More critically: CD274 predicts survival only in the
CONTEXT of checkpoint therapy. In the TCGA cohort
(predominantly pre-checkpoint era, 2008–2013), CD274
expression level should not predict OS because the
mechanism by which PD-L1 affects outcomes (checkpoint
exhaustion release) was not available as treatment.

This negative prediction separates the drug target geometry
finding (PD-L1 elevated → checkpoint therapy rationale)
from the survival prediction (PD-L1 mRNA does not predict
OS in untreated cohort).

**Predicted result:**
Log-rank p > 0.10 for CD274-high vs CD274-low within
claudin-low in the TCGA cohort. No significant OS
stratification by CD274 alone.

**Confidence: HIGH (negative prediction)**

**Falsification criterion:**
S2-P8 is wrong if CD274 significantly predicts OS
(p < 0.05) in TCGA. If it does, this likely reflects
an indirect signal (CD274 tracking overall immune
infiltration, which has its own prognostic value
as tested in S2-P3) rather than a direct PD-L1
mechanism effect.

---

## PART V — COMPLETE PREDICTION REFERENCE TABLE

| ID | Prediction | Direction | Confidence |
|----|-----------|-----------|------------|
| S2-P1 | Depth score predicts OS within claudin-low | High depth = worse OS | HIGH |
| S2-P2 | CLDN3 predicts OS within claudin-low | Low CLDN3 = worse OS | MODERATE-HIGH |
| S2-P3 | Immune score predicts OS | High immune = better OS | MODERATE |
| S2-P4 | MKI67 predicts OS within claudin-low | High MKI67 = worse OS | MODERATE |
| S2-P5 | Claudin-low OS worse than LumA, ≈ TNBC | CL < LumA, CL ≈ TNBC | MODERATE |
| S2-P6 | PARP1 does NOT predict OS independently | Null result | MODERATE (negative) |
| S2-P7 | Depth score holds in Stratum B (ESR1-low) | HR(B) ≥ HR(A) | MODERATE |
| S2-P8 | CD274 does NOT predict OS in TCGA cohort | Null result | HIGH (negative) |

---

## PART VI — WHAT SCRIPT 2 DOES NOT PREDICT

1. **pCR or neoadjuvant response in claudin-low.**
   There is no claudin-low-specific neoadjuvant dataset
   confirmed accessible. TCGA does not capture pCR.
   The TNBC Script 2 used GSE25066 for this purpose.
   Claudin-low neoadjuvant response analysis requires a
   separate dataset not yet identified. This analysis is
   deferred to a potential Script 3.

2. **Specific median OS values.**
   The claudin-low n=268 with TCGA follow-up makes point
   estimates unreliable. Predictions are directional and
   statistical. No median OS numbers are predicted.

3. **Checkpoint therapy response directly.**
   TCGA is a pre-checkpoint era cohort. No checkpoint
   response prediction is made for TCGA. The immune
   geometry (S2-P3) is tested as a prognostic marker
   in the untreated context, not as a predictive
   biomarker for checkpoint therapy.

4. **BRCA1 mutation frequency.**
   PARP1 mRNA is available; BRCA1 mutation status
   requires the TCGA mutation file. Script 2 will
   use BRCA1 mRNA as a proxy (established precedent
   from TNBC Script 2). Actual mutation-based PARP
   inhibitor prediction is not made here.

5. **Claudin-low-specific subtypes (e.g., within-claudin-low
   molecular subtypes).** The classifier impurity makes
   within-claudin-low molecular subtype analysis unreliable
   in this dataset. Subtype characterisation is deferred
   to a dedicated analysis with a larger or purer
   claudin-low dataset (e.g., GSE96058 where published
   claudin-low calls are available).

---

## PART VII — THE METHODOLOGICAL NOTE FOR SCRIPT 2

```
Script 2 must acknowledge and directly address the
classifier impurity from Script 1.

The recommended approach:

  SECTION A — Full set analysis (n=268)
    Run all survival analyses.
    Report all results.
    Flag that this is a mixed-purity population.

  SECTION B — Stratum B analysis (ESR1-low subset)
    Re-run depth score and CLDN3 survival analyses.
    Compare HR and p-values to Section A.
    Draw the purity-test conclusion:
      If Stratum B HR ≥ Stratum A HR → depth axis
      is genuine claudin-low biology.
      If Stratum B HR < Stratum A HR → depth axis
      is partly driven by LumA contamination.

  SECTION C — Cross-subtype survival (all subtypes)
    Include claudin-low (full set) alongside
    LumA, LumB, HER2, TNBC in a single KM plot.
    This is the most publishable figure in Script 2 —
    it positions claudin-low in the breast cancer
    survival landscape for the first time in this series.

The impurity is not a reason to delay Script 2.
It is a reason to design Script 2 carefully so
the conclusions are robust to the impurity.
The stratified analysis is the correct approach.
```

---

## STATUS BLOCK

```
document:           BRCA-S7c (before_script2.md)
folder:             Cancer_Research/BRCA/DEEP_DIVE/Claudin_Low/
status:             LOCKED — predictions cannot change
                    after this document is committed
date:               2026-03-05
author:             Eric Robert Lawson
                    OrganismCore
script_2_status:    NOT YET RUN
all_predictions:    stated before any Script 2 output
                    is seen
next_document:      BRCA-S7d
                    Script 2 results and reasoning artifact
                    (written after Script 2 runs)
primary_test:       S2-P1 — depth score vs OS (Stratum A)
purity_test:        S2-P7 — depth score vs OS (Stratum B)
key_negative_preds: S2-P6 (PARP1), S2-P8 (CD274)
novel_prediction:   S2-P2 (CLDN3 as single-gene survival
                    biomarker — not in standard clinical use)
```
