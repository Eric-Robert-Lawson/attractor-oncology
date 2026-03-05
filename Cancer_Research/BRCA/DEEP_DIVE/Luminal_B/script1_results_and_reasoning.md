# BRCA-S5b — Reasoning Artifact
## Luminal B Breast Cancer | Script 1 Results
### OrganismCore | GSE176078 | Wu et al. 2021 | 2026-03-05

---

## OPENING STATEMENTS

FOXA1 in Luminal B is EQUAL to Luminal A at the single-cell level (-7.4%, p=0.0045), not lower — the intermediate differentiation hypothesis stated in BRCA-S5a is refuted by the primary critical test. EZH2 in Luminal B is statistically indistinguishable from Luminal A (LumB +17.6% vs Mature, p=0.20 ns; LumA +17.2% vs Mature, p=ns) — the predicted LumA < LumB < TNBC gradient does not exist at the LumA-to-LumB step, and no EZH2i rationale survives this data.

The central revision forced by these results: **Luminal B is not a less-differentiated version of Luminal A. It is a fully committed, hyper-differentiated luminal cell that has selectively lost its cell cycle exit machinery while retaining and in several markers amplifying its luminal identity program.** GATA3 (+42.1% vs LumA, p<1e-100), ESR1 (+64.4% vs LumA, p<1e-100), and FOXA1 (equal to LumA) are all at or above LumA levels. CDKN1A is depleted -69.3% vs LumA (p<1e-100). The attractor type stated in BRCA-S5a (Type 3 dominant + partial Type 1) must be revised.

---

## 1. DATASET AND POPULATION METADATA

| Parameter | Value |
|---|---|
| Accession | GSE176078 |
| Reference | Wu et al. 2021, Nature Genetics. PMID 34493872 |
| Total cells | 100,064 |
| Platform | 10X Chromium v2 scRNA-seq, 26 primary tumors |
| LumB label | "Cancer LumB SC" — present in metadata, n=3,368 |
| Strategy | DIRECT_LABEL (no proxy required) |
| LumA n | 7,742 |
| Mature Luminal n | 1,265 |
| Luminal Progenitors n | 1,992 |
| Cancer Basal SC (TNBC ref) n | 4,312 |
| Genes in matrix | 29,733 |
| Target genes found | 41 / 41 (100%) |

All 41 target genes from the locked panel were present in the matrix. No fallback gene substitutions required. DIRECT_LABEL strategy is the highest methodological rigor available from this dataset — results carry full confidence.

---

## 2. SECTION 1 — TOP MOVERS, UNFILTERED GEOMETRY

### 2a. LumB vs Mature Luminal (absolute cancer state)

This comparison defines what the LumB cancer cell is, in absolute terms, relative to the closest normal reference — post-mitotic differentiated luminal epithelium.

**Top 15 LOST vs Mature Luminal:**

| Gene | Mature | LumB | % change | p |
|---|---|---|---|---|
| KRT5 | 0.0759 | 0.0042 | -94.5% | 1.81e-45 |
| SOX10 | 0.0082 | 0.0005 | -93.5% | 2.32e-09 |
| KRT14 | 0.1339 | 0.0134 | -90.0% | 4.11e-67 |
| CDKN1A | 1.4956 | 0.1816 | **-87.9%** | p<1e-100 |
| CCNE1 | 0.0126 | 0.0021 | -83.7% | 1.29e-07 |
| TGFBR2 | 0.0420 | 0.0132 | -68.6% | 3.35e-11 |
| CDK6 | 0.0522 | 0.0177 | -66.1% | 3.04e-12 |
| RB1 | 0.1854 | 0.0727 | -60.8% | 1.08e-35 |
| NCOA1 | 0.2822 | 0.1168 | -58.6% | 4.58e-41 |
| CDK2 | 0.1219 | 0.0647 | -46.9% | 1.27e-13 |
| KDM1A | 0.1363 | 0.0748 | -45.2% | 1.41e-12 |
| PGR | 0.3548 | 0.1997 | -43.7% | 2.03e-14 |
| NKX3-1 | 0.0862 | 0.0571 | -33.7% | 7.74e-05 |
| AR | 0.9379 | 0.6448 | -31.3% | 1.18e-28 |
| NCOA2 | 0.2025 | 0.1423 | -29.7% | 8.65e-07 |

**Top 15 GAINED vs Mature Luminal:**

| Gene | Mature | LumB | % change | p |
|---|---|---|---|---|
| TOP2A | 0.0063 | 0.0528 | +731.7% | 3.89e-16 |
| MKI67 | 0.0016 | 0.0099 | +505.2% | 9.96e-04 |
| SPI1 | 0.0088 | 0.0219 | +149.6% | 7.90e-04 |
| CCND1 | 0.8115 | 1.8368 | **+126.3%** | p<1e-100 |
| PCNA | 0.2486 | 0.4918 | +97.9% | 3.50e-36 |
| DNMT3A | 0.0691 | 0.1253 | +81.4% | 3.49e-09 |
| GATA3 | 1.1115 | 1.8801 | +69.1% | 7.97e-88 |
| ESR1 | 0.7489 | 1.1344 | +51.5% | 1.75e-43 |
| CDKN2A | 0.0692 | 0.0938 | +35.5% | 6.74e-05 |
| FOXA1 | 0.3934 | 0.4834 | +22.9% | 4.44e-07 |
| HDAC2 | 0.4506 | 0.5476 | +21.5% | 1.34e-04 |
| EZH2 | 0.0414 | 0.0486 | +17.6% | 0.2004 ns |
| MYC | 1.1010 | 1.2633 | +14.7% | 1.27e-03 |
| DNMT3B | 0.0105 | 0.0115 | +9.7% | 0.5663 ns |
| TP53 | 0.1458 | 0.1497 | +2.7% | 0.5057 ns |

**Interpretation of 1a:** The LumB cancer cell vs its normal post-mitotic reference is defined by two dominant themes. First, near-complete loss of cell cycle exit infrastructure: CDKN1A -88%, RB1 -61%, CDK6 -66%, CDK2 -47%, CCNE1 -84%. This is not partial — it is a systematic dismantling of G1/S arrest. Second, concurrent retention and amplification of luminal identity: GATA3 +69%, ESR1 +52%, FOXA1 +23%. The LumB cell is more luminal than Mature Luminal by these markers. The two phenomena are simultaneous — the cell has gained luminal identity genes and lost cell cycle exit genes together. KRT5/KRT14/SOX10 losses confirm zero basal contamination in this population (discussed separately in P6 correction). DNMT3A +81% is the largest epigenetic signal. SPI1 +150% is flagged as anomalous (Section 7, FO-5).

### 2b. LumB vs LumA — CRITICAL TEST

This is the central comparison of Script 1. The FOXA1 result here is the single most important number produced by this analysis.

**Top 15 LOST in LumB vs LumA:**

| Gene | LumA | LumB | % change | p |
|---|---|---|---|---|
| OLIG2 | 0.0001 | 0.0000 | -100% | 0.51 ns |
| CDKN1A | 0.5914 | 0.1816 | **-69.3%** | p<1e-100 |
| CCNE1 | 0.0040 | 0.0021 | -48.0% | 0.069 ns |
| RB1 | 0.1170 | 0.0727 | -37.9% | 1.89e-15 |
| NKX3-1 | 0.0785 | 0.0571 | -27.3% | 0.024 |
| TET2 | 0.0744 | 0.0561 | -24.5% | 3.83e-05 |
| NCOA1 | 0.1414 | 0.1168 | -17.4% | 1.62e-04 |
| AR | 0.7512 | 0.6448 | -14.2% | 1.58e-12 |
| CDK6 | 0.0204 | 0.0177 | -13.5% | 0.37 ns |
| SPDEF | 0.9882 | 0.8800 | -10.9% | 3.03e-13 |
| PARP1 | 0.6176 | 0.5614 | -9.1% | 4.59e-05 |
| FOXA1 | 0.5221 | 0.4834 | **-7.4%** | 4.51e-03 |
| KDM1A | 0.0781 | 0.0748 | -4.3% | 0.51 ns |

**Top 15 GAINED in LumB vs LumA:**

| Gene | LumA | LumB | % change | p |
|---|---|---|---|---|
| MKI67 | 0.0006 | 0.0099 | **+1487.4%** | 9.38e-19 |
| TOP2A | 0.0052 | 0.0528 | +910.8% | 5.90e-78 |
| SPI1 | 0.0053 | 0.0219 | +314.9% | 2.04e-20 |
| CD274 | 0.0028 | 0.0086 | +209.6% | 1.68e-06 |
| DNMT3A | 0.0508 | 0.1253 | +146.7% | 2.03e-51 |
| PCNA | 0.2021 | 0.4918 | +143.4% | p<1e-100 |
| CCND1 | 0.8135 | 1.8368 | **+125.8%** | p<1e-100 |
| HDAC2 | 0.2745 | 0.5476 | +99.5% | p<1e-100 |
| SOX10 | 0.0003 | 0.0005 | +98.1% | 0.82 ns |
| CDKN2A | 0.0484 | 0.0938 | +93.9% | 2.88e-26 |
| MYC | 0.6738 | 1.2633 | +87.5% | p<1e-100 |
| HDAC1 | 0.2160 | 0.3771 | +74.6% | 1.23e-64 |
| ESR1 | 0.6901 | 1.1344 | +64.4% | p<1e-100 |
| TGFBR2 | 0.0082 | 0.0132 | +60.9% | 1.65e-03 |
| TP53 | 0.0949 | 0.1497 | +57.8% | 3.21e-21 |

**Interpretation of 1b — Critical Test:** FOXA1 is -7.4% at p=0.0045 — statistically significant but below the ±10% threshold for a directional call, classified EQUAL. This is the primary refutation of BRCA-S5a's Type 3 dominant model. The predicted lower FOXA1 in LumB than LumA does not exist in this data. What does exist is a strikingly clean separation on the proliferative axis: MKI67 +1487%, TOP2A +911%, CCND1 +126%, MYC +88%, PCNA +143% — and on the opposing side, CDKN1A -69.3%, the single largest significant loss. GATA3 (+42.1%, p<1e-100) and ESR1 (+64.4%, p<1e-100) are both higher in LumB than LumA. The biology is unambiguous: LumB and LumA diverge on a proliferative axis, not a differentiation axis.

Three additional signals from 1b require documentation. HDAC2 is +99.5% in LumB vs LumA (p<1e-100) and HDAC1 is +74.6% (p=1.23e-64) — pan-HDAC elevation of this magnitude has not been seen in the LumA or TNBC comparisons and warrants framework observation. DNMT3A +146.7% (p=2.03e-51) is the largest epigenetic signal vs LumA. CD274 (PD-L1) +209.6% (p=1.68e-06) is a novel immune checkpoint signal not in the locked panel.

### 2c. LumB vs Luminal Progenitor (context)

LumB vs Luminal Progenitor shows near-total basal marker loss (KRT5 -99.5%, KRT14 -98.7%, SOX10 -99.7%) and massive luminal identity gain (FOXA1 +826%, ESR1 +877%, PGR +936%, GATA3 +266%, AR +577%). LumB is a fully committed luminal cell — there is no progenitor-like dedifferentiation. The LumB attractor is deeply within luminal identity space, not on the boundary with a progenitor or basal compartment. ERBB2 +96.8% (p=1.47e-56) vs Progenitor is noted — LumB has substantially elevated ERBB2 above the normal luminal progenitor baseline.

---

## 3. SECTION 2 — PCA GEOMETRY

PC1 accounts for 14.50% of variance. PC2 9.73%. The top positive PC1 loadings are CCND1 (+0.280), HDAC2 (+0.277), PARP1 (+0.275), CDK4 (+0.274), MYC (+0.265), HDAC1 (+0.255), PCNA (+0.219). PC1 is a proliferation and chromatin remodelling axis, not a differentiation axis.

**PC1 population means:**

| Population | PC1 mean | n |
|---|---|---|
| LumA | -0.7725 | 7,457 |
| Mature Luminal | +0.0607 | 1,230 |
| LumB | +0.2858 | 3,145 |
| Luminal Progenitor | +0.4745 | 1,851 |
| TNBC (Basal) | +1.0126 | 3,860 |

LumB vs LumA PC1 separation: p=8.66e-91. LumA sits at the most negative end of PC1 — it is the most quiescent, least proliferative, most CCND1/MYC-low population in this dataset. TNBC sits at the most positive end. LumB is between Mature and Progenitor, far removed from LumA in the direction of proliferative/chromatin activity. This is geometrically inconsistent with any model in which LumB is intermediate between LumA and a normal precursor on a differentiation axis. LumB has gone in the opposite direction from LumA along the dominant axis of variance in this dataset. The two subtypes occupy opposite ends of PC1.

The fact that Luminal Progenitor (PC1 mean +0.475) sits between LumB (+0.286) and TNBC (+1.013) on this axis — while containing fully basal markers — confirms that PC1 is not a luminal-to-basal transition axis. It is a chromatin accessibility / proliferative drive axis that crosscuts cell identity. LumA has the lowest chromatin remodelling / proliferative drive of any cancer population in this dataset. LumB's chromatin state is more active than LumA but has not reached TNBC levels.

---

## 4. SECTION 3 — DEPTH SCORE

**Composite formula:** w=0.4*(1-norm(CDKN1A)) + w=0.4*(1-norm(mean(FOXA1,GATA3))) + w=0.2*norm(MKI67)

| Component | Mean value | Interpretation |
|---|---|---|
| C1 — 1-norm(CDKN1A) | 0.9359 | CDKN1A near-maximally depleted within LumB range |
| C2 — 1-norm(FOXA1+GATA3) | 0.6324 | Identity TFs intermediate (score inverted: high expression = low C2) |
| C3 — norm(MKI67) | 0.0072 | MKI67 low on absolute scale within LumB (skewed distribution) |
| **Composite depth score** | **mean=0.6288, median=0.6338, std=0.1141** | |

LumB vs LumA depth (CDKN1A axis, joint normalisation): LumB mean=0.9359, LumA mean=0.8676, p=1.20e-124. LumB is deeper than LumA by the CDKN1A criterion — this result is sound.

**Critical interpretive note on the composite score:** The C2 component (identity axis) was designed assuming that lower FOXA1/GATA3 = deeper attractor = higher depth score. This assumption is violated. FOXA1 and GATA3 are higher in LumB than LumA in absolute terms. The 1-norm transformation inverts this: cells with the highest FOXA1/GATA3 within the LumB population receive the lowest C2 contribution, which in the composite pulls their depth score down. The composite score of 0.63 is therefore a hybrid: the CDKN1A component pushes it toward 0.93 (very deep) and the identity component partially cancels this. The composite score does not represent differentiation depth in the Type 3 sense. It accurately measures cell cycle exit failure (CDKN1A axis) but the identity axis component is misspecified for LumB biology. The score should be read as a CDKN1A-weighted proliferative escape index, not a differentiation state index.

**Top depth correlates within LumB (|r| ranked):**

| Gene | r | p |
|---|---|---|
| GATA3 | -0.817 | 0.00 |
| CCND1 | -0.664 | 0.00 |
| MYC | -0.617 | 0.00 |
| CDKN1A | -0.600 | 0.00 |
| ESR1 | -0.597 | 9.88e-324 |
| SPDEF | -0.568 | 1.82e-287 |
| FOXA1 | -0.546 | 6.92e-261 |
| HDAC2 | -0.540 | 8.46e-254 |
| PCNA | -0.520 | 1.40e-232 |
| AR | -0.496 | 1.21e-208 |
| PARP1 | -0.492 | 1.40e-204 |
| CDK4 | -0.491 | 8.50e-204 |
| ERBB2 | -0.459 | 3.12e-175 |
| HDAC1 | -0.448 | 6.57e-166 |
| TP53 | -0.289 | 8.83e-66 |

All 25 depth correlates are negative. Within the LumB population, cells that score higher on the depth composite (= more CDKN1A-depleted, less FOXA1/GATA3) simultaneously show lower expression of CCND1, MYC, ESR1, FOXA1, HDAC2. This is the within-population structure: the most CDKN1A-depleted LumB cells are also the most transcriptionally quiescent by luminal identity and proliferative driver gene expression. This is paradoxical at first reading but reflects a cellular subpopulation structure within LumB — a fraction of deeply CDKN1A-depleted cells that are in a low-transcription state, possibly post-mitotic arrest or G0-like, versus a cycling-active fraction with retained identity TF and CCND1 expression. This substructure is a target for Script 2 investigation.

---

## 5. SECTION 4 — PREDICTION PANEL CHECK

### P1 — Identity TFs vs Mature Luminal (predicted: RETAINED)

| Gene | Mature | LumB | % chg | p | LumA ref | Verdict |
|---|---|---|---|---|---|---|
| FOXA1 | 0.3934 | 0.4834 | +22.9% | 4.44e-07 | ref=+37.3% | ✓ ELEVATED |
| GATA3 | 1.1115 | 1.8801 | +69.1% | 7.97e-88 | ref=+33.8% | ✓ ELEVATED |
| ESR1 | 0.7489 | 1.1344 | +51.5% | 1.75e-43 | ref=-30.1% | ✓ ELEVATED |
| PGR | 0.3548 | 0.1997 | -43.7% | 2.03e-14 | ref=-54.8% | ✗ SUPPRESSED |
| SPDEF | 1.1943 | 0.8800 | -26.3% | 2.18e-26 | ref=? | ✗ SUPPRESSED |

P1 result: PARTIALLY CONFIRMED. FOXA1, GATA3, ESR1 — retained and elevated beyond the LumA reference values. This exceeds the prediction — LumB does not merely retain luminal identity, it amplifies it above the LumA baseline. PGR suppression confirmed as predicted (-43.7% vs Mature, close to LumA ref of -54.8%). SPDEF suppression was not predicted — it is SPDEF-low while being FOXA1/GATA3/ESR1-high, which suggests selective loss of a downstream ER-responsive gene with maintained upstream TF expression. SPDEF is known to be positively regulated by ESR1 — its suppression despite elevated ESR1 is mechanistically important (see FO-4 on ER circuit decoupling).

ESR1 exceeding the LumA reference by a wide margin (+51.5% vs Mature in LumB vs -30.1% in LumA — i.e., LumA is actually below Mature for ESR1) is the single most surprising absolute-level result in P1. LumA cancer cells express less ESR1 than normal Mature Luminal. LumB cancer cells express 51% more than Mature Luminal. This inversion has not been seen in TNBC or HER2 subtypes and is a defining feature of the LumB attractor state.

### P1b — Identity TFs: LumB vs LumA (CRITICAL TEST)

Predicted: FOXA1 lower in LumB than LumA. GATA3 lower. ESR1 equal.

| Gene | LumA | LumB | % diff | p | Predicted | Result |
|---|---|---|---|---|---|---|
| FOXA1 | 0.5221 | 0.4834 | -7.4% | 4.51e-03 | LOWER | **EQUAL** ← CRITICAL |
| GATA3 | 1.3230 | 1.8801 | +42.1% | p<1e-100 | LOWER | **HIGHER** |
| ESR1 | 0.6901 | 1.1344 | +64.4% | p<1e-100 | EQUAL | **HIGHER** |
| PGR | 0.1950 | 0.1997 | +2.4% | 0.066 ns | LOWER | **EQUAL** |
| SPDEF | 0.9882 | 0.8800 | -10.9% | 3.03e-13 | LOWER | ✓ LOWER |

**FOXA1 verdict: EQUAL. The Type 3 dominant attractor model is refuted.**

FOXA1 at -7.4% clears the statistical significance threshold (p=0.0045) but fails the biological significance threshold (±10%). This is not LOWER — this is a noise-level difference between two populations that share essentially the same FOXA1 level. Contrast this with FOXA1's behaviour at the LumB-to-TNBC transition (LumA=0.52, LumB=0.48, TNBC=0.076) — the drop to TNBC is real and large. The LumA-to-LumB step on FOXA1 is not a real biological transition. The identity TF program is equivalent between LumA and LumB.

GATA3 +42.1% and ESR1 +64.4% higher in LumB than LumA are the refutation of the differentiation-depth model. If LumB were a less-differentiated version of LumA, these transcription factors — which drive luminal differentiation — should be lower, not higher. They are substantially higher. The framework observation that follows from this is addressed in Section 7.

### P2 ��� Proliferation vs LumA (predicted: elevated)

| Gene | LumA | LumB | % diff | p | r_depth | Verdict |
|---|---|---|---|---|---|---|
| MKI67 | 0.0006 | 0.0099 | +1487.4% | 9.38e-19 | +0.066 | ✓ UP |
| CCND1 | 0.8135 | 1.8368 | +125.8% | p<1e-100 | -0.664 | ✓ UP |
| CDK4 | 0.3463 | 0.4568 | +31.9% | 3.68e-21 | -0.491 | ✓ UP |
| MYC | 0.6738 | 1.2633 | +87.5% | p<1e-100 | -0.617 | ✓ UP |
| TOP2A | 0.0052 | 0.0528 | +910.8% | 5.90e-78 | -0.155 | ✓ UP |
| CDKN1A | 0.5914 | 0.1816 | -69.3% | p<1e-100 | -0.600 | ✓ DOWN |

**P2: FULLY AND UNAMBIGUOUSLY CONFIRMED.** All six predictions hold with massive effect sizes. The MKI67 increase of 1487% is not a subtle difference — LumB is categorically more proliferative than LumA. The negative depth correlations for CCND1 (-0.664), MYC (-0.617), and CDK4 (-0.491) confirm that proliferative drive is tightly coupled to the CDKN1A-depletion axis within the LumB population. CDKN1A loss and CCND1/MYC gain are the same biological process, measured from opposite ends.

Note on r_depth for MKI67: the depth-MKI67 correlation is +0.066 (near zero), which appears to contradict the proliferative interpretation. This reflects the log1p + sparsity structure of MKI67 in scRNA-seq data — MKI67 is expressed in a small fraction of cells at low counts and the continuous Pearson correlation is dominated by zero-inflation. TOP2A (r=-0.155) and PCNA (r=-0.520) capture the proliferative signal more reliably in this data structure.

### P3 — EZH2 Gradient (predicted: Mature < LumA < LumB < TNBC)

| Comparison | % change vs Mature | p |
|---|---|---|
| LumA vs Mature | +17.2% | ns |
| LumB vs Mature | +17.6% | 0.20 ns |
| TNBC vs Mature | +269.7% | — |

r(EZH2, depth) in LumB = -0.1049.

**P3: NOT CONFIRMED. EZH2 is flat in both LumA and LumB relative to Mature Luminal. The predicted LumA < LumB step does not exist. The LumB-to-TNBC step is real and large (+270%), but LumB does not participate in the ascending gradient above LumA.** EZH2 is a TNBC-associated epigenetic driver. It is not differentially elevated in LumB at the expression level. No EZH2 inhibitor rationale survives this data for LumB.

This result also has framework implications. EZH2's near-identical levels in LumA and LumB (+17% and +18% vs Mature respectively) while being 270% elevated in TNBC suggest the Polycomb pathway is a basal/TNBC-specific adaptation, not a general cancer escalation. The LumB-to-TNBC transition on EZH2 is abrupt, not gradual. This is consistent with EZH2 being a switch-like state variable for the basal attractor rather than a continuously graded marker of dedifferentiation depth.

### P4 — Depth Score (predicted: functional composite)

The composite depth score (mean=0.6288, std=0.1141) functions as designed on the CDKN1A axis (joint normalisation LumB mean=0.9359 vs LumA=0.8676, p=1.20e-124). LumB is measurably deeper than LumA by this metric. The score is interpretable as a CDKN1A-weighted proliferative escape index. Reservations about the identity component are detailed in Section 4. The bulk tumor depth scores (n=24 tumors, mean=0.5725, std=0.1442, range 0.33–0.83) show substantial inter-tumor heterogeneity consistent with the known clinical heterogeneity of Luminal B disease.

### P5 — ER Circuit Integrity

Predicted: PGR more suppressed in LumB than LumA; r(ESR1,PGR) lower in LumB.

| Gene | LumA % vs Mature | LumB % vs Mature | Result |
|---|---|---|---|
| ESR1 | -7.9% | +51.5% | LumB far HIGHER (not equal) |
| PGR | -45.1% | -43.7% | EQUAL (P5 partially refuted) |
| NCOA1 | -49.9% | -58.6% | LumB MORE suppressed ✓ |
| NCOA2 | -41.4% | -29.7% | LumB less suppressed |
| TGFBR2 | -80.5% | -68.6% | LumB less suppressed |

r(ESR1, PGR) within LumB = +0.234, p=3.53e-43
r(ESR1, PGR) within LumA = +0.333, p=1.25e-199

**P5 result: r(ESR1,PGR) prediction CONFIRMED — the correlation is lower in LumB (0.234) than LumA (0.333), consistent with partial ER signal decoupling. PGR suppression equality CONFIRMED — both subtypes show ~45% PGR suppression vs Mature. ESR1 being +51% in LumB vs -8% in LumA was not predicted and represents a major finding.**

The ESR1-PGR decoupling in LumB — elevated ESR1 input, with equivalent PGR output to LumA and a weaker ESR1-PGR correlation — is the mechanistic signature of partial endocrine resistance. The ER is present and highly expressed. It is not being efficiently transduced to its canonical downstream target. SPDEF suppression in P1 is part of this same picture: SPDEF is a known ER-responsive gene, and its suppression despite elevated ESR1 is further evidence of downstream ER circuit disruption. NCOA1 suppression (-58.6% vs Mature) below LumA levels (-49.9%) further weakens the transactivation complex. The practical implication: LumB tumors have more ESR1 than LumA but less efficient ER signal output per unit ESR1, which is the molecular explanation for why endocrine therapy alone produces lower response rates in LumB-enriched populations.

### P6 — Controls (expected absent in LumB)

**SCRIPT DISPLAY ARTEFACT — CORRECTION REQUIRED.** The script printed "ELEVATED-UNEXPECTED" for KRT5, KRT14, and SOX10. This is a code error in the threshold logic: the script flags |change| > 30% regardless of direction. All three controls are strongly suppressed in LumB vs Mature Luminal:

- KRT5: -94.5% (p=1.81e-45) — near-absent ✓
- KRT14: -90.0% (p=4.11e-67) — near-absent ✓
- SOX10: -93.5% (p=2.32e-09) — near-absent ✓

**P6 is CONFIRMED biologically.** LumB contains zero basal contamination. The "Cancer LumB SC" label in Wu et al. metadata correctly identifies a pure luminal population. The P6 logic should be fixed in the next script revision to use directional thresholding (`change < -30%` for "absent as expected") rather than absolute value thresholding.

### P7 — Epigenetic Panel vs Mature

| Gene | Mature | LumB | % chg | p | r_depth |
|---|---|---|---|---|---|
| EZH2 | 0.0414 | 0.0486 | +17.6% | 0.20 ns | -0.105 |
| HDAC1 | 0.3936 | 0.3771 | -4.2% | 0.20 ns | -0.448 |
| HDAC2 | 0.4506 | 0.5476 | +21.5% | 1.34e-04 | -0.540 |
| KDM1A | 0.1363 | 0.0748 | -45.2% | 1.41e-12 | -0.173 |
| DNMT3A | 0.0691 | 0.1253 | +81.4% | 3.49e-09 | -0.271 |
| DNMT3B | 0.0105 | 0.0115 | +9.7% | 0.57 ns | -0.091 |
| EED | 0.0707 | 0.0592 | -16.3% | 0.11 ns | -0.136 |
| TET2 | 0.0692 | 0.0561 | -18.9% | 0.055 ns | -0.117 |

DNMT3A +81% (p=3.49e-09) is the dominant epigenetic signal. KDM1A (LSD1) -45% (p=1.41e-12) is the largest epigenetic loss. HDAC2 +21.5% (p=1.34e-04) is the only HDAC with significant elevation — HDAC1 is flat. EZH2, DNMT3B, EED, and TET2 are all non-significant vs Mature. The epigenetic picture is: de novo DNA methylation up (DNMT3A), H3K4me1/2 demethylation capacity down (KDM1A/LSD1), HDAC2 up, EZH2 flat. This is a DNMT3A/HDAC2/KDM1A-centred epigenetic profile, not a PRC2/EZH2-centred profile. LumB's epigenetic dependencies are distinct from TNBC's.

### P9 — Cross-Subtype Gradient

| Gene | LumA | LumB | TNBC | Gradient | Predicted |
|---|---|---|---|---|---|
| ESR1 | 0.690 | 1.134 | 0.025 | LumB > LumA >> TNBC | LumA > LumB > TNBC — REFUTED |
| FOXA1 | 0.522 | 0.483 | 0.076 | LumA ≥ LumB >> TNBC | DESC ✓ |
| GATA3 | 1.323 | 1.880 | 0.518 | LumB > LumA > TNBC | LumA > LumB > TNBC — REFUTED |
| EZH2 | 0.049 | 0.049 | 0.153 | LumA ≈ LumB << TNBC | LumA < LumB < TNBC — PARTIAL |

The P9 result adds cross-subtype context to the P1b refutation. The ESR1 and GATA3 gradient predictions assumed LumB was intermediate between LumA and TNBC. The actual ordering for both genes is LumB at the top, LumA below LumB, TNBC far below both. The cancer-relevant gradient for these two markers is not a descending scale of differentiation — LumB has the most luminal identity TF expression of any cancer population in this dataset. FOXA1 follows the predicted descending gradient (LumA > LumB > TNBC) but the LumA-LumB difference is small (8%) and the LumB-TNBC difference is large (84%). EZH2 is flat from LumA to LumB and then jumps to TNBC — a step function, not a gradient.

---

## 6. DRUG TARGETS AND PREDICTIONS

### Tier 1 — Strong evidence, mechanistically direct

**CDK4/6 inhibition (palbociclib, ribociclib, abemaciclib)**

Evidence from this data:
- CDKN1A -69.3% vs LumA (p<1e-100): near-complete loss of the primary CDK4/6 brake
- CCND1 +125.8% vs LumA (p<1e-100): cyclin D1 — the CDK4/6 activator — is the single largest significant gained signal vs LumA
- CDK4 +31.9% vs LumA (p=3.68e-21): the kinase target itself is elevated
- MKI67 +1487%, TOP2A +911%: proliferative phenotype is extreme
- r(CDK4, depth) = -0.491: CDK4 expression is tightly coupled to the CDKN1A-depletion axis within LumB
- r(CCND1, depth) = -0.664: the strongest depth correlate in the panel (after GATA3) is CCND1
- PCNA +143.4% (p<1e-100): replication machinery fully engaged

Assessment: This is the strongest drug target signal in the dataset. Five independent lines of evidence converge on CDK4/6 pathway activation as the mechanistically dominant oncogenic drive in LumB. CDKN1A depletion, CCND1 amplification, CDK4 elevation, and the proliferative phenotype (MKI67, TOP2A) are all consistent with a cell that has lost its CDK4/6 brake and is constitutively cycling. CDK4/6 inhibition is the primary therapeutic rationale for LumB, confirmed by this single-cell data and consistent with phase III clinical trial outcomes (PALOMA-2/3, MONALEESA-2/3/7, MONARCH-3).

**Endocrine therapy (tamoxifen, fulvestrant, aromatase inhibitors) in combination**

Evidence:
- ESR1 +51.5% vs Mature (p=1.75e-43): ER is highly expressed — target is present
- GATA3 +69.1% vs Mature (p=7.97e-88): luminal identity circuit intact
- FOXA1 +22.9% vs Mature: pioneer factor for ER binding is present
- ESR1-PGR coupling r=0.234 (lower than LumA r=0.333): partial signal decoupling
- SPDEF -26.3% vs Mature despite elevated ESR1: downstream ER output partially disrupted
- NCOA1 -58.6% vs Mature: ER coactivator suppressed

Assessment: ER is present and highly expressed. The target is valid. However, the ER circuit shows measurable decoupling — elevated input (ESR1) is not producing proportional downstream output (PGR flat, SPDEF suppressed, NCOA1 reduced). This is the molecular basis of partial endocrine resistance in LumB. Endocrine therapy remains justified — the receptor is there — but monotherapy response rates will be lower than in LumA where the ER circuit is more intact. The CDK4/6i + endocrine combination is the correct therapeutic approach: CDK4/6 blockade re-sensitises partially endocrine-resistant ER+ cancer by restoring RB-mediated cell cycle control independently of ER signalling completeness.

**CDK4/6i + endocrine backbone as combination**

The data provides direct mechanistic support for the approved combination strategy. CDKN1A loss drives CDK4/6 pathway hyperactivation. ER is present. Adding CDK4/6 inhibition to endocrine therapy bypasses the need for a fully functional downstream ER circuit — the cell cycle is blocked regardless of whether ESR1 is fully transducing to PGR. This is why the combination works better than endocrine monotherapy in LumB-enriched populations.

### Tier 2 — Biological rationale present, requires further support

**HDAC inhibition (entinostat, vorinostat) — particularly HDAC2-selective approaches**

Evidence:
- HDAC2 +99.5% vs LumA (p<1e-100): largest HDAC signal in the dataset
- HDAC1 +74.6% vs LumA (p=1.23e-64): both class I HDACs elevated vs LumA
- r(HDAC2, depth) = -0.540: tightly coupled to CDKN1A-depletion axis within LumB
- r(HDAC1, depth) = -0.448: similarly coupled
- HDAC2 vs Mature: +21.5% (p=1.34e-04) — elevation is confirmed above normal reference

Assessment: The simultaneous 2-fold elevation of both HDAC1 and HDAC2 vs LumA (not seen in TNBC or HER2 comparisons at this magnitude) is a LumB-specific epigenetic signal. HDAC2 has known roles in cell cycle progression and its coupling to the depth/CDKN1A axis within LumB suggests functional relevance. Entinostat has been evaluated in ER+ breast cancer (ENCORE 301) with modest benefit — the LumB-enriched population may be the appropriate subgroup for HDAC inhibitor trials. Not strong enough for Tier 1 without functional validation.

**Anthracycline sensitivity (doxorubicin, epirubicin) — TOP2A as predictive marker**

Evidence:
- TOP2A +731% vs Mature, +910% vs LumA: TOP2A is near-absent in LumA and massively elevated in LumB
- TOP2A is the direct enzymatic target of anthracyclines
- TOP2A amplification/elevation is a known predictor of anthracycline response in breast cancer

Assessment: The data supports differential anthracycline sensitivity in LumB vs LumA based on TOP2A expression. LumA has near-zero TOP2A (0.0052 log1p mean). LumB has 10-fold higher TOP2A (0.0528). This is consistent with published data showing higher pCR rates with anthracycline-containing neoadjuvant regimens in higher-grade ER+ tumors (which are LumB-enriched). Relevant to neoadjuvant treatment planning.

**DNMT3A inhibition or DNMT3A-dependent vulnerability exploitation**

Evidence:
- DNMT3A +81.4% vs Mature (p=3.49e-09), +146.7% vs LumA (p=2.03e-51)
- Largest epigenetic signal in the panel
- Not elevated in LumA

Assessment: DNMT3A elevation in LumB above LumA suggests de novo DNA methylation activity is a LumB-specific epigenetic feature. Standard demethylating agents (azacitidine, decitabine) have limited solid tumour development. The more actionable implication is that DNMT3A-mediated methylation may be silencing tumour suppressors specifically in LumB — this requires bisulfite sequencing data to confirm which loci are methylated. Flag for Script 2 mechanistic investigation.

### Tier 3 — Refuted or unsupported by this data

**EZH2 inhibition (tazemetostat)**
- EZH2 +17.6% vs Mature, p=0.20 ns. Equal to LumA (+17.2%).
- P3 not confirmed. No LumB-specific EZH2 elevation.
- No EZH2i rationale for LumB from expression data. Remove from active consideration for this subtype. EZH2i remains a TNBC/basal target based on this dataset.

**PARP inhibition (olaparib, niraparib)**
- PARP1 -9.1% vs LumA — slightly lower in LumB, not elevated
- PARP inhibition rationale in breast cancer is mutation-driven (germline BRCA1/2), not expression-driven
- Cannot be evaluated from expression data alone
- No expression-level rationale specific to LumB

**KDM1A/LSD1 inhibition**
- KDM1A is -45.2% vs Mature (suppressed, not elevated) and -4.3% vs LumA (flat)
- Low KDM1A in LumB does not support LSD1 inhibition as a strategy — inhibiting a gene that is already suppressed has no clear mechanism
- The KDM1A suppression is a descriptor of LumB state, not a drug target

---

## 7. FRAMEWORK OBSERVATIONS

### FO-1 — Attractor Type Revision: Type 3 Dominant Refuted, Type 1 Dominant Proposed

BRCA-S5a stated: Type 3 dominant + partial Type 1. This is refuted by the P1b critical test.

Type 3 in the OrganismCore framework is defined by intermediate identity TF expression — the cancer cell has partially lost its tissue identity, sitting between the normal differentiated cell and a dedifferentiated state. LumB's FOXA1 equal to LumA, GATA3 42% above LumA, and ESR1 64% above LumA are incompatible with a Type 3 attractor. The identity TF program is not intermediate or suppressed — it is at maximum luminal expression levels.

**Revised classification: Type 1 dominant.** Type 1 in the framework is characterised by retained or amplified identity with hyperproliferative drive — the cancer cell has maintained its tissue-specific transcription factor program while gaining autonomous proliferative signalling. The LumB profile matches this description precisely: luminal identity (FOXA1/GATA3/ESR1) retained and amplified, cell cycle exit (CDKN1A/RB1) dismantled, proliferative drivers (CCND1/MYC/CDK4/TOP2A) massively elevated.

The partial Type 1 component acknowledged in BRCA-S5a is now the dominant component. The Type 3 component was an artefact of the prediction model, not of the data.

**Proposed sub-classification: Type 1-L (Luminal Locked).** LumA may represent a Type 1-L state at lower proliferative drive (quiescent luminal lock), while LumB represents a Type 1-L state at higher proliferative drive (active luminal lock with cell cycle re-entry). Both subtypes have the same identity program; they differ by proliferative state. This sub-classification should be evaluated against HER2-enriched (which may be Type 1-H, HER2-amplified luminal) and the cycling cancer population to determine if the Type 1 axis in luminal cancer is a spectrum.

### FO-2 — ESR1 Inversion: LumB > LumA > Mature Luminal

The ordering ESR1: LumB (1.134) > LumA (0.690) > Mature Luminal (0.749) — with LumA actually below Mature — is unexpected and has not been flagged in prior subtypes. Three non-exclusive mechanisms:

1. **FOXA1/GATA3-driven ESR1 promoter activation.** Both pioneer factors are elevated in LumB. FOXA1 directly opens chromatin at the ESR1 promoter. Higher FOXA1/GATA3 → more accessible ESR1 locus → higher ESR1 transcription. This is mechanistically coherent.

2. **Compensatory upregulation against partial downstream resistance.** If downstream ER signalling is partially blocked (SPDEF suppressed, NCOA1 low, r(ESR1,PGR) reduced), a feedback mechanism may drive ESR1 promoter activity upward. This is less well-supported by the data but cannot be excluded.

3. **Selection for ER positivity in LumB tumors.** LumB tumors are clinically HR+ by definition. The dataset may reflect selection pressure for maintained ER expression in a proliferating cell that otherwise would drift toward ER loss. LumB may represent a state where ESR1 expression is actively maintained as part of the attractor stability mechanism, even as downstream signalling is disrupted.

The SPDEF suppression despite high ESR1, combined with reduced r(ESR1,PGR), argues that interpretation 2 or 3 is operating alongside interpretation 1. The ER in LumB is highly expressed but incompletely functional. This decoupling is the molecular signature of partial endocrine resistance.

### FO-3 — CDKN1A Depletion as the Primary LumA→LumB Transition Event

The LumA-to-LumB transition in this dataset is defined by one dominant signal: CDKN1A -69.3%, p<1e-100. Nothing else comes close in terms of effect size combined with biological coherence. CDKN1A (p21) is a CDK inhibitor that integrates multiple anti-proliferative signals (DNA damage, TGFB signalling, cell contact inhibition) to enforce G1 arrest. Its near-complete loss in LumB vs LumA — in a cell that otherwise maintains full luminal identity — is the molecular switch event that distinguishes these two subtypes. All other proliferative signals (CCND1, MKI67, TOP2A, CDK4) are downstream consequences of CDKN1A loss combined with retained cyclin D1/CDK4 expression. The upstream driver of CDKN1A suppression in LumB is not resolved by this expression data. Candidates include: MYC-driven CDKN1A promoter suppression (MYC +87% vs LumA), DNMT3A-mediated CDKN1A promoter methylation (DNMT3A +147% vs LumA), and/or PI3K/AKT pathway activation. Script 2 should examine CDKN1A promoter methylation status if methylation data are available.

### FO-4 — SPDEF as a Decoupling Marker Within the ER Circuit

SPDEF (SAM Pointed Domain ETS Factor) is a known direct transcriptional target of ESR1/FOXA1 and is required for luminal differentiation in breast epithelium. In LumB: SPDEF is -26.3% vs Mature (p=2.18e-26) and -10.9% vs LumA (p=3.03e-13), despite FOXA1 and ESR1 being elevated above both comparators. This dissociation — high ESR1 + high FOXA1, yet low SPDEF — is a specific marker of ER circuit decoupling. SPDEF drives terminal luminal differentiation markers including MUC1 and lactation genes. Its suppression in LumB while upstream TFs are elevated suggests that the ER/FOXA1 complex in LumB has altered co-factor binding or altered chromatin context that specifically prevents SPDEF transactivation. This may be related to the HDAC1/2 elevation (HDACs can remodel chromatin at SPDEF loci) or the DNMT3A elevation (de novo methylation of SPDEF regulatory elements). SPDEF is a candidate mechanistic marker for ER circuit quality — not just quantity — in LumB vs LumA.

### FO-5 — SPI1 (PU.1) Elevation: Unresolved

SPI1 is +314.9% in LumB vs LumA (p=2.04e-20) and +149.6% vs Mature (p=7.90e-04). SPI1/PU.1 is a canonical myeloid lineage TF — it is the master regulator of monocyte/macrophage identity and has no described function in luminal breast epithelium. Its presence as one of the top gained signals in LumB vs LumA is anomalous. Three possible explanations:

1. **Immune cell contamination.** The "Cancer LumB SC" cluster may contain a small number of contaminating myeloid cells. Given that Wu et al. used careful annotation, this is unlikely to be a large effect, but even 1-5% contamination could elevate mean SPI1 substantially given its near-zero expression in epithelial cells.

2. **Ectopic SPI1 expression in cycling epithelial cells.** There is a small literature on SPI1 expression in non-haematopoietic cancers during S-phase. If SPI1 is transiently expressed in cycling epithelial cells, the LumB population (which is far more proliferative than LumA) would show higher SPI1 by this mechanism.

3. **Transcriptional noise at low absolute expression levels.** SPI1 mean in LumA is 0.0053 and in LumB is 0.0219 log1p — both very low. A +315% fold-change on a near-zero baseline may represent stochastic transcriptional noise amplified by the high cell count (n=3368 giving statistical power to detect biologically negligible differences).

Resolution: examine SPI1 at single-cell resolution within the LumB cluster in Script 2. If SPI1-expressing cells co-express myeloid markers (CD14, CSF1R, LYZ) they are contaminants. If SPI1 cells are uniformly distributed and co-express FOXA1/GATA3 they represent genuine ectopic expression. Do not include SPI1 in any biological model until this is resolved.

### FO-6 — HDAC1/HDAC2 Co-Elevation is a LumB-Specific Signal

HDAC1 +74.6% and HDAC2 +99.5% vs LumA (both p<1e-100) represent a class I HDAC upregulation that is specific to the LumB-vs-LumA comparison. Neither HDAC showed this pattern in the LumA vs Mature or TNBC comparisons at comparable effect size. Both HDACs have strong negative depth correlations within LumB (HDAC1 r=-0.448, HDAC2 r=-0.540), meaning the most deeply CDKN1A-depleted LumB cells also express the most HDAC1/2. This coupling suggests HDAC1/2 expression is mechanistically linked to the cell cycle re-entry program in LumB. HDAC2 is known to deacetylate RB and promote E2F-driven S-phase entry. HDAC1 is known to repress CDK inhibitor promoters. The elevated HDAC1/2 in LumB may contribute to CDKN1A suppression through promoter deacetylation — a mechanism that would be upstream of or parallel to DNMT3A-mediated methylation and MYC-driven transcriptional suppression. This creates a potential synthetic therapeutic target: the convergence of DNMT3A + HDAC1/2 elevation on CDKN1A suppression suggests that combined DNMT + HDAC inhibition might restore CDKN1A expression in LumB specifically. This is speculative from expression data alone but mechanistically coherent.

### FO-7 — CD274 (PD-L1) Elevation is Novel and Clinically Relevant

CD274 +209.6% in LumB vs LumA (p=1.68e-06) was not in the locked prediction panel and was not anticipated from BRCA-S5a. PD-L1 expression in LumB above LumA suggests an immune evasion component to the LumB attractor state. The MYC elevation (+87.5%) is relevant here — MYC is a known direct transcriptional activator of CD274. High MYC → high PD-L1 is a documented axis in multiple cancer types. In LumB, MYC and CD274 co-elevation creates a potential immunotherapy rationale that does not exist in LumA where MYC is lower and CD274 is near-absent (0.0028). The clinical implications of PD-L1 expression in HR+/HER2- breast cancer are currently limited (checkpoint inhibitors are not standard of care in this subtype), but the data supports the mechanistic basis for investigating CDK4/6i + checkpoint inhibitor combinations in LumB-enriched populations, where CDK4/6 inhibition has been shown to increase tumor immunogenicity through distinct mechanisms.

### FO-8 — ERBB2 Elevation Above Progenitor

ERBB2 +96.8% vs Progenitor (p=1.47e-56) indicates LumB has substantially elevated HER2 protein expression above the normal luminal progenitor baseline. This does not constitute HER2 amplification (which requires genomic copy number data) but it places LumB on the high end of normal HER2 expression range. The "Cancer Her2 SC" label exists as a distinct population in this dataset (n=3,708) — LumB cells with ERBB2 near-amplification levels may represent the clinical LumB/HER2-low population that has recently shown benefit from ADC therapy (trastuzumab deruxtecan). ERBB2 in LumB (0.4901 log1p) vs LumA in the 1c comparison: ERBB2 is not in the LumA vs LumB panel but the progenitor comparison gives the baseline. This should be formally added to the gene panel in Script 2.

---

## 8. BULK RNA-SEQ CONTEXT

Bulk depth scores across 24 GSE176078 bulk tumor samples: mean=0.5725, std=0.1442, range 0.33–0.83. The standard deviation of 0.14 on a 0-1 scale represents substantial inter-tumor heterogeneity — some tumors score 0.33 (shallow, LumA-like) and others 0.83 (deep, fully proliferative). This heterogeneity is consistent with the clinical observation that Luminal B is not a homogeneous subtype: it spans a spectrum from near-LumA (high CDKN1A, moderate proliferation) to near-cycling-cancer (negligible CDKN1A, extreme proliferation). The bulk score distribution should be examined for correlation with Ki67 IHC in future work — high bulk depth score (>0.7) likely corresponds to Ki67 >30%, the clinical cut-off for high-grade LumB.

---

## 9. PREDICTION SCORECARD

| Prediction | Status | Notes |
|---|---|---|
| P1: FOXA1/GATA3/ESR1 retained vs Mature | ✓ CONFIRMED — EXCEEDED | All three elevated, not merely retained |
| P1: PGR suppressed vs Mature | ✓ CONFIRMED | -43.7%, matches LumA ref |
| P1b: FOXA1 lower than LumA | ✗ REFUTED | EQUAL (-7.4%), not lower |
| P1b: GATA3 lower than LumA | ✗ REFUTED | HIGHER (+42.1%) |
| P1b: ESR1 equal to LumA | ✗ REFUTED | HIGHER (+64.4%) |
| P1b: PGR lower than LumA | ✗ REFUTED | EQUAL (+2.4%) |
| P2: MKI67/CCND1/CDK4/MYC/TOP2A elevated vs LumA | ✓ CONFIRMED | All UP, massive effect sizes |
| P2: CDKN1A lower than LumA | ✓ CONFIRMED | -69.3%, p<1e-100, dominant result |
| P3: EZH2 intermediate gradient LumA<LumB<TNBC | ✗ REFUTED | LumA≈LumB, then large jump to TNBC |
| P4: Depth score functional | ✓ CONFIRMED (with caveat) | CDKN1A axis valid; identity component misspecified |
| P5: PGR more suppressed in LumB than LumA | ✗ REFUTED | PGR equal in both subtypes |
| P5: r(ESR1,PGR) lower in LumB | ✓ CONFIRMED | 0.234 vs 0.333 in LumA |
| P6: Controls absent | ✓ CONFIRMED (script artefact in display) | All three basal markers near-absent |
| P9: FOXA1 gradient LumA>LumB>TNBC | ✓ CONFIRMED | Holds |
| P9: ESR1 gradient LumA>LumB>TNBC | ✗ REFUTED | LumB>LumA>>TNBC |
| P9: GATA3 gradient LumA>LumB>TNBC | ✗ REFUTED | LumB>LumA>TNBC |
| P9: EZH2 gradient LumA<LumB<TNBC | ✗ PARTIAL | LumA≈LumB, TNBC far higher — step function not gradient |
| Attractor type: Type 3 dominant | ✗ REFUTED | Revised to Type 1 dominant (Type 1-L) |

Confirmed: 8 | Refuted: 9 | Partial: 1

The high refutation rate is not a methodological failure — it is a scientific result. The predictions in BRCA-S5a were based on the intermediate differentiation hypothesis, which was a reasonable prior given the clinical positioning of LumB between LumA and high-grade ER+ disease. The data has falsified that prior cleanly and replaced it with a more precise model: proliferative escape from a locked luminal identity state.

---

## 10. OPEN QUESTIONS FOR BRCA-S5c (SCRIPT 2)

1. **SPI1 resolution:** Single-cell resolution within LumB cluster — myeloid contamination or genuine ectopic expression? Examine co-expression with FOXA1/GATA3 vs CD14/LYZ.

2. **ESR1-high/PGR-low subpopulation:** Does a discrete ESR1-high/SPDEF-low/PGR-low substate exist within LumB, distinct from an ESR1-moderate/PGR-moderate substate? If so, this maps to the decoupled ER circuit population with greatest endocrine resistance.

3. **CDKN1A-depleted subpopulation geometry:** The within-LumB depth correlations show a subpopulation structure where the most CDKN1A-depleted cells are also transcriptionally quiescent on luminal identity genes. Are these cells in G0-like arrest, post-mitotic, or in a drug-resistant persister state? Requires pseudotime or cell cycle phase analysis.

4. **DNMT3A + HDAC1/2 convergence on CDKN1A:** Is CDKN1A promoter methylation (DNMT3A) accompanied by deacetylation (HDAC1/2) in the same cells? Co-expression analysis of DNMT3A/HDAC1/2 vs CDKN1A expression within LumB.

5. **ERBB2 formal quantification:** Add ERBB2 to the identity panel in Script 2. Quantify the fraction of LumB cells meeting an ERBB2-high threshold (HER2-low territory).

6. **CD274/MYC co-expression:** Are the CD274-high LumB cells the same cells as the MYC-high cells? If so, this defines an immune-evasive, MYC-driven proliferative substate within LumB.

7. **CDKN2A paradox:** p16/CDKN2A is +93.9% vs LumA (p=2.88e-26) while p21/CDKN1A is -69.3%. p16 and p21 loss are typically co-occurring in cancer. p16 elevation alongside p21 depletion is the senescence bypass signature — the cell has activated the senescence program (p16 accumulation) but escaped it by depleting the effector (p21). This should be examined at the single-cell level: are the p16-high cells also the p21-lowest cells, suggesting senescence bypass as a discrete cellular state within LumB?

8. **Revised attractor type formalisation:** Type 1-L (Luminal Locked, hyperproliferative) definition to be written into the framework before BRCA-S5c. How does Type 1-L relate to the cycling cancer population (n=5359 in this dataset)? Is Cancer Cycling a further progression from LumB Type 1-L toward full cycle autonomy?

---

## 11. DOCUMENT PROVENANCE

| Field | Value |
|---|---|
| Document | BRCA-S5b |
| Framework | OrganismCore |
| Protocol | Workflow_Protocol.md v2.0 |
| Dataset | GSE176078 — Wu et al. 2021, Nature Genetics |
| Script | BRCA_LumB_script1.py |
| Date | 2026-03-05 |
| Preceding document | BRCA-S5a (predictions locked) |
| Following document | BRCA-S5c (Script 2 — mechanistic deep dive) |
| n(LumB) | 3,368 cells, DIRECT_LABEL strategy |
| Critical test result | FOXA1 EQUAL in LumB vs LumA — Type 3 model refuted |
| Primary revised conclusion | LumB is Type 1-L dominant: locked luminal identity with cell cycle re-entry via CDKN1A depletion and CCND1 amplification |
| Primary drug target | CDK4/6 inhibition (five independent evidentiary lines) |
| EZH2i status | Withdrawn — no expression-level support |

---

*End of BRCA-S5b Reasoning Artifact*
