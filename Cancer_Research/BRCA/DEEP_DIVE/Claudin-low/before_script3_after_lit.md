# CLAUDIN-LOW — SCRIPT 3 BEFORE-DOCUMENT
## Predictions Locked Before Script 3 Runs
## OrganismCore — Document BRCA-S7f
## Date: 2026-03-05

---

## DOCUMENT METADATA

```
document_id:        BRCA-S7f
series:             BRCA Deep Dive — Claudin-Low
folder:             Cancer_Research/BRCA/DEEP_DIVE/Claudin_Low/
type:               BEFORE-DOCUMENT (Script 3)
date:               2026-03-05
author:             Eric Robert Lawson
                    OrganismCore
protocol_version:   Workflow_Protocol.md v2.0
attractor_type:     TYPE 4 — ROOT LOCK
                    (ATTRACTOR_GEOMETRY_AXIOMS.md v2.0)
based_on:           BRCA-S7b (script1_results_and_reasoning.md)
                    BRCA-S7d (script2_results_and_reasoning.md)
                    BRCA-S7e (literature_check.md)
status:             LOCKED — predictions cannot change
                    after this document is committed
```

---

## PART I — WHY SCRIPT 3 EXISTS

Script 2 completed the survival analysis. The literature check
(BRCA-S7e) produced three findings that generate testable
predictions not addressed in Scripts 1 or 2:

**Finding 1 — Morel JCI 2017:**
Anti-PD-1 monotherapy is counterproductive in claudin-low because
Tregs in the claudin-low TME express high PD-1, and anti-PD-1
paradoxically makes them MORE suppressive. The key variable is
not total immune infiltration (tested in S2-P3, null result) but
the RATIO of Treg markers to effector T cell markers. The
FOXP3/CD8A ratio measures how Treg-dominated vs. effector-T-cell-
dominated the microenvironment is. This ratio was never computed
in Scripts 1 or 2.

**Finding 2 — Pommier Nature Comm 2020:**
Claudin-low contains at least three molecular subgroups with
different cells of origin: mammary stem cell-like (subgroup 1 —
true TYPE 4), luminal-derived via EMT (subgroup 2), and EMT-like
(subgroup 3). Subgroups 2 and 3 retain residual lineage memory
of the committed identity they came from. This residual memory is
measurable in expression space: samples with detectable FOXA1,
SPDEF, or GATA3 above the claudin-low set median retain partial
luminal identity (subgroup 2/3 proxy). Samples with no detectable
luminal TF expression (all three below claudin-low median) have
no lineage memory (subgroup 1 proxy — true root lock). These
populations have never been tested for different survival outcomes
within the claudin-low set.

**Finding 3 — Axioms Prediction G (novel, untested):**
CT antigen de-repression (GAGE family, CT45 family, STRA8, DPPA2)
was the dominant signal in the Script 1 unfiltered scan. The Axioms
document (Prediction G) proposes that CT antigen de-repression
should correlate with depth in TYPE 4 cancers because it reflects
the degree to which the cell has never committed to a somatic
identity programme. This correlation was never tested. If confirmed,
CT antigen level becomes both a prognostic marker and a direct
immunotherapy target (GAGE CAR-T, CT45-directed vaccines).

Script 3 addresses all three findings with locked predictions.

---

## PART II — WHAT SCRIPT 3 USES

```
DATA:
  TCGA-BRCA expression matrix (already downloaded)
  TCGA-BRCA pancan survival (already downloaded)
  Claudin-low depth scores (recomputed from Script 1
  logic — identical parameters)

POPULATIONS:
  Stratum A: full geometry set (n=268)
  Stratum B: ESR1-low subset (n=177) — primary analysis
  Stratum C: Basal+Normal-PAM50 (n=79) — confirmatory

NO NEW DATA DOWNLOADS REQUIRED.
All analyses use files already present from Scripts 1 and 2.

OUTPUT:
  Claudin_Low_s3_results/
    cl_s3_log.txt
    cl_s3_figure.png
    cl_s3_scorecard.csv
    cl_s3_immune_ratios.csv
    cl_s3_lineage_memory.csv
    cl_s3_ct_antigen.csv
```

---

## PART III — ANALYSIS STRUCTURE

Script 3 runs three independent analyses. Each is derived from
a specific literature check finding.

```
ANALYSIS A — TREG/EFFECTOR IMMUNE RATIO
  Source: Morel JCI 2017
  Question: Does the ratio of Treg markers to effector T cell
  markers predict OS in claudin-low, and does it track depth?

  Genes:
    Treg/suppressor: FOXP3, PDCD1 (PD-1 on Tregs), TIGIT
    Effector:        CD8A, GZMB, PRF1
  Ratios computed:
    FOXP3 / CD8A
    FOXP3 / GZMB
    TIGIT  / CD8A
    (FOXP3 + TIGIT) / (CD8A + GZMB + PRF1)
    [composite Treg:effector ratio]
  Tests:
    Correlation of each ratio with depth score
    Tertile survival (Stratum B — ESR1-low subset)
    Cox: ratio + depth score jointly (does ratio add
    independent prognostic information?)

ANALYSIS B — CT ANTIGEN DEPTH CORRELATION AND SURVIVAL
  Source: Axioms Prediction G (BRCA-S7e novel finding)
  Question: Does CT antigen expression correlate with depth
  within claudin-low, and does it predict OS?

  Genes (confirmed elevated in Script 1 unfiltered scan):
    GAGE1, GAGE2D, GAGE4, GAGE12D, GAGE12J
    CT45A3, CT45A4
    STRA8, DPPA2
  CT antigen composite score:
    Mean z-score of available CT antigen genes within set
  Tests:
    r(each CT antigen gene, depth score)
    r(CT composite score, depth score)
    Tertile survival for CT composite score (Stratum B)
    Correlation of CT score with FOXP3/CD8A ratio
    (does CT antigen load correlate with Treg dominance?)

ANALYSIS C — LINEAGE MEMORY SUBGROUP ANALYSIS
  Source: Pommier Nature Comm 2020
  Question: Do claudin-low samples with residual luminal
  lineage memory (Pommier subgroup 2/3 proxy) have better
  survival than samples with no lineage memory (subgroup 1
  proxy — true root lock)?

  Lineage memory score:
    Computed as mean(FOXA1, SPDEF, GATA3) within claudin-low.
    Samples above claudin-low set median = MEMORY-HIGH
    (residual luminal identity — subgroup 2/3 proxy)
    Samples below claudin-low set median = MEMORY-LOW
    (no residual luminal identity — subgroup 1 proxy)
  Tests:
    Survival: memory-high vs memory-low (Stratum B)
    CT antigen score in memory-high vs memory-low
    (does no-memory group have more CT antigen de-repression?)
    Depth score in memory-high vs memory-low
    (does memory-low = deeper depth score?)
    Treg/effector ratio in memory-high vs memory-low
    (does no-memory group have more Treg-dominated TME?)
```

---

## PART IV — PREDICTIONS
### All stated 2026-03-05 before Script 3 runs

---

### ANALYSIS A — TREG/EFFECTOR RATIO

---

#### S3-P1 — FOXP3/CD8A RATIO CORRELATES POSITIVELY WITH DEPTH SCORE
**Direction:** r(FOXP3/CD8A, depth) > 0.20 within claudin-low
(Stratum A)

**Geometric basis:**
Depth score measures how completely the cell has abandoned
differentiation commitment (claudin loss + mesenchymal gain +
luminal TF loss). The immune microenvironment in TYPE 4 cancers
is predicted by Axiom IV.2 to be exhaustion-dominant, and this
exhaustion should be more complete in deeper tumours. Deeper
tumours have more completely abandoned epithelial identity, more
CT antigen de-repression (attracting immune cells), AND a more
developed Treg-suppression mechanism. FOXP3/CD8A ratio should
therefore rise with depth: deeper claudin-low = more Treg-dominated
= higher ratio.

**Predicted result:**
r(FOXP3/CD8A ratio, depth score) > 0.20 (Stratum A, n=268).
Pearson and Spearman both positive.

**Confidence: MODERATE**

**Falsification criterion:**
S3-P1 not confirmed if r < 0.10 or negative.

---

#### S3-P2 — FOXP3/CD8A RATIO PREDICTS WORSE OS IN STRATUM B
**Direction:** High FOXP3/CD8A ratio = WORSE overall survival
(HR > 1.0)

**Geometric basis:**
The Morel JCI 2017 mechanism: the Treg-dominated microenvironment
suppresses effector T cell activity AND paradoxically amplifies
under anti-PD-1 monotherapy. Independently of therapy, higher
Treg:effector ratio = less spontaneous immune clearance = worse
untreated outcomes. This extends the Script 2 null result for total
immune score (S2-P3): the TOTAL immune score was null because it
mixed Tregs (bad) and effector T cells (good). The RATIO separates
these opposing effects.

This is the primary mechanistic prediction from the literature check
that was not testable in Scripts 1 or 2.

**Predicted result:**
Log-rank p < 0.10 for FOXP3/CD8A tertile in Stratum B (ESR1-low,
n≈176). HR > 1.0 for high-ratio group (worse OS).

**Confidence: MODERATE-HIGH**

**Falsification criterion:**
S3-P2 not confirmed if p > 0.15 or HR < 1.0 (wrong direction).
If FOXP3/CD8A-high does BETTER — that would be unexpected and would
require explanation (possibly: high CD8A in the denominator drives
the ratio, and CD8A is favourable, masking the FOXP3 effect).

---

#### S3-P3 — COMPOSITE TREG:EFFECTOR RATIO PREDICTS OS MORE STRONGLY
### THAN INDIVIDUAL IMMUNE GENES
**Direction:** The composite ratio (FOXP3+TIGIT)/(CD8A+GZMB+PRF1)
predicts OS in Stratum B with greater HR and lower p-value than
any single immune gene alone

**Geometric basis:**
A composite ratio captures the balance of the immune microenvironment
more accurately than any single gene. FOXP3 alone reflects Treg
quantity. CD8A alone reflects effector quantity. The ratio is the
functional variable — the net immunosuppressive index.

**Predicted result:**
Composite ratio HR in Stratum B > individual gene HRs for FOXP3,
CD8A, GZMB, TIGIT separately. At minimum, composite ratio p-value
is lower than the best individual gene p-value.

**Confidence: MODERATE**

**Falsification criterion:**
S3-P3 not confirmed if the composite ratio HR is not the largest
of all immune gene and ratio tests.

---

### ANALYSIS B — CT ANTIGEN

---

#### S3-P4 — CT ANTIGEN GENES CORRELATE POSITIVELY WITH DEPTH SCORE
**Direction:** r(GAGE1, depth) > 0.20; r(CT antigen composite,
depth) > 0.25 within Stratum A

**Geometric basis:**
Axioms Prediction G: CT antigen de-repression is a TYPE 4 structural
marker. The mechanism is that a cell which has never committed to
a somatic identity programme retains epigenetic access to germline
loci that committed cells silence. In the claudin-low depth axis,
deeper cells are more completely pre-commitment (more bilateral
identity absence, more mesenchymal, more stem-like). Deeper cells
should therefore have less somatic epigenetic silencing of CT loci,
producing higher CT antigen expression. The correlation should be
positive: depth rises → CT antigen rises.

This is the first direct test of Axioms Prediction G in any cancer.

**Predicted result:**
r(CT composite score, depth score) > 0.25 (Pearson, Stratum A).
At least 3 of the 9 individual CT antigen genes (GAGE1, GAGE2D,
GAGE4, GAGE12D, GAGE12J, CT45A3, CT45A4, STRA8, DPPA2) show
r > 0.20 with the depth score.

**Confidence: MODERATE**
Novel prediction — no prior data to anchor the expected magnitude.

**Falsification criterion:**
S3-P4 not confirmed if CT composite r < 0.15, or fewer than 2
individual CT genes reach r > 0.15. A null correlation would
suggest CT antigen de-repression in claudin-low is driven by a
mechanism other than the depth axis (e.g., epigenetic accident,
replication stress, unrelated to the stem lock geometry).

---

#### S3-P5 — CT ANTIGEN COMPOSITE SCORE PREDICTS WORSE OS IN STRATUM B
**Direction:** CT antigen-high claudin-low has WORSE overall survival

**Geometric basis:**
If CT antigen score tracks depth (S3-P4), and depth predicts worse
OS (confirmed in S2-P7, HR=3.112), then CT antigen score should
also predict worse OS. The CT antigen score is a depth proxy with
the additional property of being a direct immunotherapy target.
Higher CT antigen = deeper root lock = worse prognosis = more
antigen for GAGE CAR-T or vaccine to target. Both the prognostic
and therapeutic implications are aligned.

**Predicted result:**
Log-rank p < 0.10 for CT composite tertile in Stratum B. HR > 1.0
for CT-high. Effect size expected to be similar to or smaller than
depth score HR (3.112) because CT antigen score is one component
of the depth biology, not the full signal.

**Confidence: MODERATE**
Conditional on S3-P4 confirming.

**Falsification criterion:**
S3-P5 not confirmed if p > 0.15 or HR < 1.0.

---

#### S3-P6 — CT ANTIGEN SCORE CORRELATES POSITIVELY WITH
### FOXP3/CD8A RATIO
**Direction:** r(CT composite, FOXP3/CD8A) > 0.15

**Geometric basis:**
Both CT antigen de-repression and Treg dominance are predicted to
be features of deep TYPE 4 false attractors. If both track the
depth axis (S3-P1 and S3-P4), they should also correlate with each
other. Deeper tumours have both more CT antigen and more Treg
infiltration. This is a consistency test: if S3-P1 and S3-P4 both
confirm, S3-P6 should follow automatically. If S3-P6 does not
confirm when S3-P1 and S3-P4 do, the depth axis is splitting into
two independent components — immunological and epigenetic — which
would itself be a significant finding.

**Predicted result:**
r(CT composite, FOXP3/CD8A ratio) > 0.15 (Pearson, Stratum A).

**Confidence: MODERATE**

**Falsification criterion:**
Not confirmed if r < 0.10. A near-zero correlation despite both
correlating with depth would indicate independent mechanisms
converging on the same depth axis from different directions.

---

### ANALYSIS C — LINEAGE MEMORY SUBGROUPS

---

#### S3-P7 — MEMORY-LOW GROUP HAS WORSE OS THAN MEMORY-HIGH GROUP
**Direction:** Claudin-low samples with no residual luminal lineage
memory (FOXA1/SPDEF/GATA3 all below claudin-low set median) have
WORSE overall survival than memory-high samples

**Geometric basis:**
Pommier 2020 subgroup 1 (direct stem cell origin — no committed
precursor) is the deepest TYPE 4 geometry. Subgroups 2 and 3
(luminal-derived or EMT-derived) retain residual lineage memory
from the committed state they came from. This residual memory is
a remnant of the committed identity's normal arrest machinery
(partial luminal identity may retain partial CDKN1A, partial TGF-β
signalling). The memory-low group (full bilateral absence, no
luminal remnant) has therefore lost even the partial arrest capacity
that subgroups 2/3 retain. This predicts worse clinical behaviour
for memory-low.

**Predicted result:**
Log-rank p < 0.10 for memory-high vs memory-low in Stratum B.
HR > 1.0 for memory-low (worse OS). Effect size expected to be
comparable to or larger than depth score HR (3.112) because
lineage memory score targets the most biologically meaningful
depth gradient within the claudin-low set.

**Confidence: MODERATE-HIGH**

**Falsification criterion:**
S3-P7 not confirmed if p > 0.15 or HR < 1.0.

---

#### S3-P8 — MEMORY-LOW GROUP HAS HIGHER CT ANTIGEN SCORE
**Direction:** Memory-low claudin-low samples have higher CT antigen
composite score than memory-high samples

**Geometric basis:**
Memory-low = more complete pre-commitment state = less somatic
epigenetic silencing = more CT antigen de-repression. This is the
mechanistic link between Pommier subgroup 1 (direct stem origin)
and the CT antigen signature. If true, CT antigen de-repression is
not just a depth marker — it is specifically a marker of the direct-
stem-origin subtype of claudin-low, distinguishing it from the EMT-
derived subtypes.

**Predicted result:**
Two-sample t-test or Mann-Whitney U: mean CT antigen composite
score is higher in memory-low vs memory-high (p < 0.05).

**Confidence: MODERATE**
Conditional on S3-P4 confirming that CT antigen tracks depth at all.

**Falsification criterion:**
Not confirmed if p > 0.10 in either direction. If memory-HIGH has
higher CT antigen — that would be unexpected and would require
reinterpretation of what the memory score is measuring.

---

#### S3-P9 — MEMORY-LOW GROUP HAS DEEPER DEPTH SCORES
**Direction:** Memory-low claudin-low samples have higher depth
scores than memory-high samples

**Geometric basis:**
The lineage memory score (mean FOXA1/SPDEF/GATA3) should be
negatively correlated with the depth score: deeper samples have
less luminal TF expression (FOXA1 r(depth)=-0.430, GATA3
r(depth)=-0.507, established in Script 1). Memory-low = less luminal
TF = deeper depth score. This is a direct consistency check between
the two scoring approaches.

**Predicted result:**
Two-sample test: depth score mean/median is higher in memory-low
vs memory-high (p < 0.01, expected to be highly significant given
the direct relationship between FOXA1/GATA3 and the depth score
components).

**Confidence: HIGH**
This is a near-mathematical necessity given the Script 1 depth
correlations. If it does not confirm, there is an error in the
score computation.

**Falsification criterion:**
Not confirmed if memory-low depth scores are not significantly
higher (p > 0.05). This would indicate the lineage memory score
is measuring something orthogonal to the depth axis — a meaningful
finding if true.

---

#### S3-P10 — MEMORY-LOW GROUP HAS HIGHER FOXP3/CD8A RATIO
**Direction:** Memory-low claudin-low samples have higher
FOXP3/CD8A ratio than memory-high samples

**Geometric basis:**
If memory-low = deeper TYPE 4 geometry, and deeper TYPE 4 =
more Treg-dominated TME (S3-P1), then memory-low should have the
most Treg-dominated immune microenvironment of all claudin-low
subgroups. This links the Pommier subgroup biology (cell of origin)
to the Morel 2017 immune mechanism (Treg dominance): the cells that
arose directly from stem cells (no committed precursor, no lineage
memory) are also the cells that generate the most Treg-dominated
immune environment — and are therefore the most resistant to
anti-PD-1 monotherapy and the most likely to benefit from
anti-TIGIT-first sequencing.

**Predicted result:**
Memory-low FOXP3/CD8A ratio > memory-high FOXP3/CD8A ratio
(p < 0.10, two-sample test).

**Confidence: MODERATE**

**Falsification criterion:**
Not confirmed if p > 0.15 or direction reversed.

---

## PART V — COMPLETE PREDICTION REFERENCE TABLE

| ID | Prediction | Analysis | Direction | Confidence |
|----|-----------|----------|-----------|------------|
| S3-P1 | FOXP3/CD8A ratio correlates with depth | A | r > 0.20 | MODERATE |
| S3-P2 | FOXP3/CD8A ratio predicts OS in Stratum B | A | High ratio = worse OS | MODERATE-HIGH |
| S3-P3 | Composite ratio strongest immune predictor | A | HR(composite) > HR(single genes) | MODERATE |
| S3-P4 | CT antigen composite correlates with depth | B | r > 0.25 | MODERATE |
| S3-P5 | CT antigen score predicts OS in Stratum B | B | High CT = worse OS | MODERATE |
| S3-P6 | CT antigen score correlates with Treg ratio | B | r > 0.15 | MODERATE |
| S3-P7 | Memory-low has worse OS than memory-high | C | Memory-low = worse | MODERATE-HIGH |
| S3-P8 | Memory-low has higher CT antigen score | C | Memory-low > memory-high | MODERATE |
| S3-P9 | Memory-low has deeper depth scores | C | Memory-low deeper | HIGH |
| S3-P10 | Memory-low has higher FOXP3/CD8A ratio | C | Memory-low more Treg-dominant | MODERATE |

---

## PART VI — WHAT SCRIPT 3 DOES NOT PREDICT

1. **METABRIC or external dataset validation.** Script 3 uses
   TCGA-BRCA only. The external dataset validation (EXT-P1 through
   EXT-P4 from BRCA-S7d) is deferred.

2. **Mutation data.** BRCA1 mutation frequency in memory-low vs
   memory-high is not tested. Mutation data requires the TCGA MAF
   file which is not required for these analyses.

3. **Single-cell resolution of Pommier subgroups.** The lineage
   memory score is a bulk RNA-seq proxy for Pommier subgroups 1
   vs 2/3. True subgroup assignment requires single-cell RNA-seq
   or the mouse model data from Pommier 2020. This analysis is a
   bulk RNA-seq approximation.

4. **Anti-TIGIT clinical response prediction.** TCGA has no
   treatment or response data. The mechanistic prediction (anti-TIGIT
   before anti-PD-1 in high-ratio patients) cannot be tested in TCGA.

---

## STATUS BLOCK

```
document:           BRCA-S7f (before_script3.md)
folder:             Cancer_Research/BRCA/DEEP_DIVE/Claudin_Low/
status:             LOCKED
date:               2026-03-05
author:             Eric Robert Lawson
                    OrganismCore
script_3_status:    NOT YET RUN
all_predictions:    stated before any Script 3 output is seen
next_document:      BRCA-S7g
                    Script 3 results and reasoning artifact
primary_tests:      S3-P2 (FOXP3/CD8A ratio vs OS, Stratum B)
                    S3-P7 (lineage memory vs OS, Stratum B)
axioms_test:        S3-P4 (CT antigen correlates with depth —
                    first direct test of Axioms Prediction G)
lit_check_test:     S3-P2 (Morel 2017 mechanism test)
                    S3-P7 (Pommier 2020 subgroup test)
```
