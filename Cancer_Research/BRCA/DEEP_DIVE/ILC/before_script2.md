# ILC — SCRIPT 2 BEFORE-DOCUMENT
## Predictions Locked Before Script 2 Runs
## OrganismCore — Document BRCA-S6c
## Date: 2026-03-05

---

## DOCUMENT METADATA

- **Cancer type:** Invasive Lobular Carcinoma (ILC) — BRCA subtype
- **Script:** Before Script 2 (predictions locked)
- **Attractor type:** TYPE 3 VARIANT — ADHESION LOCK DISSOLUTION
- **Based on:** BRCA-S6b (script1_results_and_reasoning.md)
- **Protocol version:** v2.0 (predictions locked before data contact)
- **Status:** LOCKED — Script 2 has not run at time of writing

---

## PART I — WHAT SCRIPT 1 ESTABLISHED

These are confirmed facts entering Script 2. They are not predictions —
they are the geometric substrate from which predictions are derived.

**Confirmed geometry:**
1. Luminal TF programme is hyperactivated ABOVE normal:
   ESR1 +9.4%, FOXA1 +15.4%, GATA3 +11.9%, SPDEF +19.0%
2. CDH1 mRNA partially reduced (-14.5%); protein loss is complete
   via mutation in ~65% and methylation in ~35% of ILC
3. No basal or mesenchymal programme: all EMT/basal markers reduced
4. EZH2 elevated +8.9%; weak negative correlation with CDH1 (r=-0.147)
5. Low proliferation vs TNBC and HER2 (MKI67: ILC=9.77 vs TNBC=11.99)
6. PTEN reduced -4.3% (p=2.00e-19); PI3K context real
7. CDH1 depth axis is orthogonal to ESR1 identity axis
8. Structural inversion with TNBC confirmed in cross-subtype table
9. PCA: ILC clusters nearest LumA (distance 1.296 vs 8.314 to TNBC)
10. CCND1 elevated +3.4% — CDK4/6 context present

---

## PART II — WHAT SCRIPT 2 MUST DO

Script 2 uses the same TCGA-BRCA dataset from Script 1.
It focuses on outcomes, survival, and treatment-context validation.

**Required analyses:**
1. OS survival analysis: ESR1-high vs ESR1-low within ILC
2. OS survival analysis: CDH1-depth (high vs low) within ILC
3. OS survival analysis: ILC vs LumA (are they clinically equivalent?)
4. SPDEF as a novel stratification marker within ILC
5. EZH2 tertile survival within ILC
6. PIK3CA mRNA as a proxy for pathway context (mutation data if available)
7. PTEN-low survival within ILC
8. Cross-subtype survival: all 5 groups compared

Script 2 does NOT require external datasets. All analyses use TCGA
clinical + survival data cross-referenced with the expression matrix
already downloaded.

---

## PART III — PREDICTIONS
### All stated 2026-03-05 before Script 2 runs

---

### S2-P1 — ESR1 LEVEL PREDICTS SURVIVAL IN ILC
**Direction:** ESR1-high ILC has BETTER overall survival than ESR1-low ILC

**Geometric basis:**
ESR1 is the dominant active transcription programme in ILC. It is elevated
above normal, confirming the attractor is deeply committed to the ER
state. ER-positive cancers treated with endocrine therapy have better
outcomes when ER expression is higher. The geometric prediction is that
within ILC, ESR1 level is a continuous survival predictor — not just
ER+/ER- binary.

**Predicted result:**
Log-rank test for ESR1-high (top tertile) vs ESR1-low (bottom tertile)
within ILC: p < 0.05, hazard ratio favoring ESR1-high.

**Falsification criterion:**
P1 is NOT confirmed if log-rank p > 0.10 or HR is in the wrong direction
(ESR1-low doing better).

---

### S2-P2 — CDH1 DEPTH SCORE PREDICTS WORSE SURVIVAL IN ILC
**Direction:** Deep ILC (low CDH1) has WORSE overall survival

**Geometric basis:**
CDH1 loss is the structural dissolution event. Deeper CDH1 loss = more
complete dissolution of the adhesion constraint = more invasive potential.
Deep ILC is more likely to spread in single-file peritoneal pattern, which
is associated with worse outcomes than ductal spread patterns.

However, this prediction has LOWER confidence than S2-P1. Reasons:
- CDH1 mRNA is an incomplete proxy for protein loss (as established in S6b)
- TCGA cohort OS follow-up may be insufficient for ILC (a slow-growing cancer)
- Endocrine therapy confounds CDH1 prognosis (all ESR1+ patients are treated)

**Predicted result:**
Log-rank p < 0.10 (lower confidence threshold than P1), CDH1-low worse.

**Falsification criterion:**
P2 is NOT confirmed if log-rank p > 0.20 or HR is in the wrong direction.
A null result here is interpretable — CDH1 mRNA may be too weak a proxy.

---

### S2-P3 — ILC HAS WORSE SURVIVAL THAN LumA DESPITE SIMILAR GEOMETRY
**Direction:** ILC overall survival is WORSE than LumA in the same cohort

**Geometric basis:**
PCA shows ILC and LumA are nearest neighbors (distance 1.296). They share
the same luminal TF programme. However, ILC has CDH1 structural dissolution
that LumA does not. This enables a different metastatic pattern in ILC:
peritoneal spread, diffuse bone marrow infiltration, retroperitoneal
involvement — patterns that are harder to detect and treat than the
pulmonary/hepatic pattern of IDC.

The geometric prediction is that the structural dissolution imposes a
survival cost that is not captured by the shared identity profile.
Two cancers that look similar on identity markers can have different
outcomes if their structural geometries differ.

**Predicted result:**
Median OS or survival probability at 5 years: ILC < LumA, log-rank
p < 0.05 OR a clear trend (p < 0.15) in the direction of ILC worse.

**Falsification criterion:**
P3 is NOT confirmed if ILC and LumA are statistically indistinguishable
in survival (p > 0.20, no trend). A null result means the shared luminal
identity dominates outcomes over the structural difference.

---

### S2-P4 — SPDEF IS A NOVEL PROGNOSTIC MARKER WITHIN ILC
**Direction:** SPDEF-high ILC has BETTER survival (deeper luminal
commitment = more endocrine therapy responsive)

**Geometric basis:**
SPDEF showed the largest TF elevation in Script 1 (+19.0%, p=1.58e-28).
SPDEF is downstream of FOXA1 and GATA3 and drives terminal luminal
differentiation. Higher SPDEF = deeper terminal luminal state = more
fully committed to the ER attractor = more endocrine therapy responsive.

This is a novel prediction. SPDEF is not a standard clinical ILC biomarker.
If confirmed, it would represent a geometry-derived discovery.

**Predicted result:**
SPDEF-high (top tertile) vs SPDEF-low (bottom tertile) within ILC:
log-rank p < 0.05, HR favoring SPDEF-high.

**Falsification criterion:**
P4 is NOT confirmed if log-rank p > 0.10 in either direction.

---

### S2-P5 — PTEN-LOW WITHIN ILC PREDICTS WORSE SURVIVAL
**Direction:** PTEN-low ILC has WORSE overall survival

**Geometric basis:**
PTEN mRNA is reduced -4.3% (p=2.00e-19) in bulk ILC vs normal.
Within ILC, the subset with lowest PTEN has the most active PI3K pathway.
PI3K pathway activation in ER+ breast cancer is associated with endocrine
therapy resistance and worse outcomes. This is established biology. The
geometric prediction confirms and quantifies the PTEN axis within ILC.

**Predicted result:**
Log-rank p < 0.05, PTEN-low worse survival within ILC.

**Falsification criterion:**
P5 is NOT confirmed if log-rank p > 0.10 or HR is in the wrong direction.

---

### S2-P6 — EZH2-HIGH WITHIN ILC DOES NOT PREDICT WORSE SURVIVAL
**Direction:** EZH2 level is NOT a significant survival predictor within ILC

**Geometric basis:**
EZH2 elevation in ILC is modest (+8.9%) and its correlation with CDH1
is weak (r=-0.147). Unlike TNBC where EZH2 is strongly elevated and
associated with aggressive biology, in ILC EZH2 elevation is a partial
signal from the methylation-driven subset. The majority of ILC (mutation-only)
has near-normal EZH2. Therefore EZH2 should not stratify ILC survival.

**This is a NEGATIVE prediction.** Framework negative predictions are as
important as positive ones.

**Predicted result:**
Log-rank p > 0.10 for EZH2-high vs EZH2-low within ILC. No significant
survival difference.

**Falsification criterion:**
P6 is NOT confirmed (i.e., prediction is wrong) if log-rank p < 0.05.
If EZH2 DOES predict ILC survival, that is a novel finding requiring
explanation — likely identifying the methylation-driven subset as a
distinct clinical entity.

---

### S2-P7 — CCND1-HIGH WITHIN ILC PREDICTS BETTER SURVIVAL
**Direction:** CCND1-high ILC has BETTER survival
(paradoxical — higher cyclin D1 = more CDK4/6 inhibitor sensitive)

**Geometric basis:**
CCND1 is elevated +3.4% in ILC vs normal (p=0.0001). CCND1 amplification
and overexpression in ER+ breast cancer predicts CDK4/6 inhibitor
sensitivity. High CCND1 in the context of high ESR1 indicates a cell that
is cycling through CDK4/6-dependent proliferation pathways — which are
targetable. The geometric prediction is that CCND1-high ILC represents
the endocrine-sensitive, CDK4/6-inhibitable subgroup with better outcomes
in the treatment era.

Note: this prediction is treatment-era dependent. In a pre-CDK4/6i cohort,
CCND1-high might not predict better survival. TCGA is predominantly a
pre-treatment-era cohort. Interpret with caution.

**Predicted result:**
Trend toward CCND1-high better survival, p < 0.15 (lower confidence
given treatment-era confound). May be null in TCGA era.

**Falsification criterion:**
P7 is NOT confirmed if p > 0.20 in either direction.

---

### S2-P8 — ILC SURVIVAL IS NOT PREDICTED BY MKI67
**Direction:** MKI67 level does NOT significantly stratify ILC survival

**Geometric basis:**
MKI67 in ILC is low across the board (9.77 vs TNBC 11.99). The range
within ILC is narrow. Low-grade ILC is Grade 1-2 by definition. Survival
difference between high-MKI67 and low-MKI67 ILC may be real but should
be smaller in magnitude than the MKI67 survival split in TNBC or HER2,
because ILC has a compressed proliferation distribution.

**This is a conditional negative prediction** — MKI67 may predict
within ILC but with weaker effect size than in other subtypes.

**Predicted result:**
Log-rank p > 0.05 within ILC for MKI67-high vs MKI67-low. If significant,
effect size (HR) should be smaller than in TNBC.

**Falsification criterion:**
P8 is wrong if log-rank p < 0.02 with large HR in ILC.

---

## PART IV ��� COMPLETE PREDICTION REFERENCE

| ID | Prediction | Direction | Confidence |
|----|-----------|-----------|-----------|
| S2-P1 | ESR1 predicts OS within ILC | High ESR1 = better | HIGH |
| S2-P2 | CDH1 depth predicts OS within ILC | Low CDH1 = worse | MODERATE |
| S2-P3 | ILC worse survival than LumA | ILC < LumA | MODERATE |
| S2-P4 | SPDEF novel prognostic marker | High SPDEF = better | MODERATE (novel) |
| S2-P5 | PTEN-low predicts worse OS within ILC | Low PTEN = worse | HIGH |
| S2-P6 | EZH2 does NOT predict OS within ILC | Null result | MODERATE (negative) |
| S2-P7 | CCND1-high trends toward better OS | High CCND1 = better | LOW (treatment confound) |
| S2-P8 | MKI67 does NOT strongly predict OS in ILC | Null or weak | MODERATE (negative) |

---

## PART V — WHAT SCRIPT 2 DOES NOT PREDICT

1. **Specific survival curves or median OS values.** The TCGA ILC cohort
   (n=204) with incomplete follow-up makes point estimates unreliable.
   The predictions are directional and statistical, not numerical.

2. **PIK3CA mutation status directly.** TCGA mutation data requires a
   separate file. Script 2 will use PIK3CA mRNA and PTEN mRNA as proxies.
   If the mutation file is accessible, PIK3CA mutation vs wild-type
   survival will be tested as a bonus analysis.

3. **Neoadjuvant therapy response.** ILC is typically not treated with
   neoadjuvant chemotherapy. This analysis is not applicable.

4. **Metastatic site data.** TCGA does not routinely capture peritoneal
   vs pulmonary metastatic patterns. The framework prediction about ILC
   metastatic pattern cannot be tested in TCGA.

---

## STATUS BLOCK

- Document BRCA-S6c: COMPLETE (predictions locked)
- Script 2: NOT YET RUN
- All predictions above are stated before any Script 2 output is seen
- Timestamp: 2026-03-05
