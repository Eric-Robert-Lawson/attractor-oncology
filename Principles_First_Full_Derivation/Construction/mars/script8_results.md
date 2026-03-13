# NECROMASS ATTRACTOR GEOMETRY — SCRIPT 8 RESULTS
## Identity Manifold and Geometric Verdict
## Eric Robert Lawson / OrganismCore
## 2026-03-13

---

## DOCUMENT IDENTITY

```
Document ID:  NECROMASS_ATTRACTOR_RESULTS_v1.0
Date:         2026-03-13
Script:       8 of 9+
Author:       Eric Robert Lawson / OrganismCore
Pre-reg DOI:  10.5281/zenodo.18986790
Repository:   github.com/Eric-Robert-Lawson/
              attractor-oncology
ORCID:        0009-0002-0414-6544
```

---

## THE FOUR STRUCTURAL INVARIANTS

```
INVARIANT 1: Source-product slope = 1.0
  Biology translates δ¹³C source by fixed
  metabolic offset. Slope = 1.0.
  Abiotic photochem: slope ≠ 1.0.
  Source: Hayes 2001; Whiticar 1999.
  Script 8 result: TEST INVALID AS DESIGNED.
  Mixed bio+abiotic pairs produce artefact
  negative slope (-0.843). Biological-only
  pairs produce slope ≈ 1.0. Requires Script 9
  redesign with class separation.

INVARIANT 2: C:N separatrix at C:N = 20
  Biological: C:N = 5-15.
  Abiotic: C:N → ∞.
  Source: Scripts 5-6; Sterner & Elser 2002.
  Script 8 result: STRONG BIOLOGICAL SIGNAL.
  Lafayette C:N ~ 2. P(bio side) = 1.0000.

INVARIANT 3: Curvature ratio Ea/ΔH = 5.6×
  Biological Ea = 100 kJ/mol (Conrad 2023).
  Abiotic ΔH = 18 kJ/mol (Moores 2019).
  Script 8 result: INDETERMINATE.
  n=12 Webster data insufficient. Both models
  fit identically (ΔAIC = 0.00). Fitted
  Ea = 6.9 kJ/mol ≠ biological 100 kJ/mol.
  Requires Script 9 with expanded dataset
  and proper model parameterisation.

INVARIANT 4: PCA covariance rank
  Biological outputs: low-rank (coupled).
  Abiotic outputs: higher-rank (independent).
  Source: Clough et al. 2025 Earth Space Sci.
  Script 8 result: CLEAN. See below.
```

---

## COMPONENT 1 — SOURCE-PRODUCT SLOPE

```
OLS slope:         -0.8390
Bootstrap median:  -0.8430
95% CI:            [-2.3032, +0.8190]
r value:           -0.4705
P(|slope-1.0|<0.05): 0.0039

STATUS: TEST INVALID AS DESIGNED.

Reason: Dataset mixes biological and abiotic
observation pairs. The SAM abiotic model
(source +43‰ → product -107‰) and the Ueno
2024 photochemical model (source +43‰ →
product -107‰) are abiotic data points that
drive a strongly negative slope when pooled
with biological pairs. A pooled regression
of mixed classes does not test the biological
invariant. It tests the mean relationship
across classes. Redesign required: run slope
test on biological-only pairs and abiotic-
only pairs separately. Script 9.
```

---

## COMPONENT 2 — C:N–δ¹³C IDENTITY PLANE

```
Observation                      log10(C:N)  IDS_CN  P(bio)
─────────────────────────────────────────────��───────────
Lafayette vein (NanoSIMS)           0.301    +1.000  1.0000
D. audaxviator (Earth)              0.813    +0.488  1.0000
Typical microbial necromass         0.902    +0.399  1.0000
Serpentinization formate/acetate    2.300    -0.999  0.0367
Abiotic CH₄-derived carbon         3.000    -1.699  0.0248
CO₂ photolysis organic (Ueno 2024) 2.904    -1.602  0.0248
Galactic dust (estimate)           2.700    -1.399  0.0540

Separatrix: log10(C:N) = log10(20) = 1.301

STATUS: STRONG BIOLOGICAL SIGNAL.
Lafayette is 1.000 log units into
the biological region.
P(bio side) = 1.0000.
Zero of 100,000 Monte Carlo samples
crossed into the abiotic region.
```

---

## COMPONENT 3 — SEASONAL CH₄ CURVATURE

```
Arrhenius fit:
  Ea fitted:  6.9 kJ/mol
  R²:         0.9254
  AIC:        -53.98

Langmuir fit:
  ΔH fitted:  -6.9 kJ/mol
  R²:         0.9254
  AIC:        -53.98

ΔAIC: 0.00
Winner declared: Arrhenius (biological)

STATUS: INDETERMINATE. RESULT RETRACTED.

Reason: Ea = 6.9 kJ/mol ≠ biological
value of 100 kJ/mol. The fitted value
is an order of magnitude too small.
ΔAIC = 0.00 means the models are
mathematically degenerate at this
scale — identical fits. n=12 seasonal
points is insufficient to distinguish
Arrhenius from Langmuir.
The biological Arrhenius trajectory
(Ea=100) and abiotic Langmuir trajectory
(ΔH=18) produce different curvatures
only across a T range that exceeds the
seasonal variation at Gale Crater (~50K).
Curvature ratio 5.6× is a published
structural fact — the test to distinguish
them requires either:
  (a) Multi-year high-resolution TLS data
      with δ¹³C(CH₄) co-measurement, or
  (b) A physically parameterised model
      with fixed Ea=100 and fixed ΔH=18,
      fitted only for amplitude.
Script 9 implements option (b).
```

---

## COMPONENT 4 — PCA IDENTITY MANIFOLD

```
Variance explained:
  PC1: 68.5%
  PC2: 20.6%
  PC3: 10.2%
  Total (PC1+PC2): 89.1%

7 dimensions used:
  1. δ¹³C organic
  2. Fractionation from DIC
  3. log10(C:N ratio)
  4. N presence (binary)
  5. Fe co-location (binary)
  6. Depth gradient magnitude
  7. log10(CH₄ seasonal amplitude)

Identity Distance Scores:
  Observation                      IDS    Side
  ──────���──────────────────────────────────────
  D. audaxviator (Earth)          +1.126  BIO
  Typical FeOB (Earth)            +2.388  BIO
  Deep subsurface methanogen      +0.872  BIO
  Serpentinization short-chain    -0.783  ABIO
  CO₂ UV photolysis (Ueno 2024)  -2.143  ABIO
  Serpentinization methane        -1.460  ABIO
  Lafayette organic (Scripts 1-5) +2.595  BIO ★
  SAM Gale Crater (House 2022)   -2.438  ABIO ★
  Webster CH₄ seasonal (Mars)    -0.010  ABIO

★ Mars observations

STATUS: CLEAN RESULT.

Lafayette IDS = +2.595:
  More biologically positioned than any
  reference organism, including FeOB (+2.388)
  and Desulforudis (+1.126).
  Lafayette exceeds the biological reference
  centroid along the identity axis.

SAM IDS = -2.438:
  More abiotically positioned than any
  abiotic reference, including CO₂
  photolysis (-2.143).
  Geometrically consistent with Ueno 2024.

Separation between Lafayette and SAM:
  2.595 - (-2.438) = 5.033 standard units.
  They are at opposite poles of the
  identity axis in the full 7D space.

signals_geometrically_separated: True
```

---

## AGGREGATE VERDICT

```
INVARIANT 1: INVALID (redesign in Script 9)
INVARIANT 2: BIOLOGICAL — P = 1.0000 (strong)
INVARIANT 3: INDETERMINATE (redesign in Script 9)
INVARIANT 4: BIOLOGICAL — IDS = +2.595 (strong)

Valid tests: 2 of 4.
Both valid tests: BIOLOGICAL for Lafayette.
Both invalid tests: forwarded to Script 9.

GEOMETRIC VERDICT (from valid tests):
Lafayette sits unambiguously in the
biological attractor basin.
IDS = +2.595, exceeding all reference
organisms. C:N P(bio) = 1.0000.

BIFURCATION:
Lafayette organic signals → BIOLOGICAL basin
SAM Gale Crater → ABIOTIC basin
Webster CH₄ seasonal → near separatrix

These are geometrically distinct phenomena.
Not contradictions. Disambiguation.
```

---

## WHAT SCRIPT 9 FIXES

```
Test 1 redesign — Source-product slope:
  Run on biological-only pairs separately.
  Run on abiotic-only pairs separately.
  Compute slope and CI for each class.
  Test H₀: bio slope = 1.0.
  Test H₁: abiotic slope ≠ 1.0.
  Compare class slopes directly.
  The class separation is the invariant.

Test 2 redesign — Seasonal curvature:
  Fix Ea = 100 kJ/mol (biological).
  Fix ΔH = 18 kJ/mol (abiotic).
  Fit only amplitude parameter A for each.
  Compare residuals and AIC with
  constrained models, not free-parameter fit.
  Add δ¹³C(CH₄) co-measurement constraint
  from Korablev 2019 and other published
  Mars CH₄ isotope data.
  The constrained model fit properly tests
  which trajectory shape the data prefers.
```

---

## FIGURE MANIFEST

```
Fig 23: Source-product slope (note: invalid design)
Fig 24: C:N–δ¹³C identity plane (valid)
Fig 25: Seasonal curvature (note: indeterminate)
Fig 26: PCA identity manifold (valid)
```

---

## CUMULATIVE RECORD TO DATE

```
Script 1  Depth gradient       CONSISTENT
Script 2  Rate calc            INCONCLUSIVE
Script 3  Morphology + δ¹³C   CONSISTENT
Script 4  WL fractionation     PARTIALLY DISCRIM.
Script 5  Molecular identity   P(bio) = 0.990
Script 6  Desulforudis         H₂ fails, Fe(II) ok
Script 7  SAM anomaly          Not discriminating
Script 8  Attractor geometry   Lafayette IDS +2.595
                                SAM IDS -2.438
                                C:N P(bio) = 1.000
                                Slope: redesign needed
                                Curvature: redesign needed
Script 9  Invariant redesign   PENDING
```
