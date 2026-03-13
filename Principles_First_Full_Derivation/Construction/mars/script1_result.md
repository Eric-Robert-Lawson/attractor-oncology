# NECROMASS ANALYSIS — SCRIPT 1 RESULTS
## Iron Oxide Deposition vs. Depth in Martian Rock
## Test of the Mars Iron Necromass Hypothesis
## Against Published Nakhlite Meteorite Data
### Eric Robert Lawson / OrganismCore
### 2026-03-13

---

## DOCUMENT IDENTITY

```
Document ID:  NECROMASS_ANALYSIS_RESULTS_v1.0
Date:         2026-03-13
Author:       Eric Robert Lawson / OrganismCore
Status:       ACTIVE — ANALYSIS RESULTS RECORD
Script:       necromass_analysis.py v1.0
Diagnostic:   necromass_diagnostic.py v1.0
Pre-reg DOI:  10.5281/zenodo.18986790
Repository:   github.com/Eric-Robert-Lawson/
              attractor-oncology
ORCID:        0009-0002-0414-6544
```

---

## WHAT WAS TESTED

The Mars Iron Necromass Hypothesis predicts:

> Iron oxide accumulation in Martian subsurface rock
> increases with depth in the original Mars rock pile,
> because a subsurface biological iron-cycling community
> operates deeper in the Martian crust, mining iron from
> olivine and basalt and precipitating it as iron oxide.
> Deeper rock spent more time in the biological iron-cycling
> zone. Therefore: deeper rock = more total iron oxide
> deposited.

**The null hypothesis (abiotic):**

> Iron oxide accumulation is uniform with depth or decreases
> with depth. Top-down surface weathering delivers oxidants
> from the surface downward. Shallower rock receives more
> oxidant exposure. Deeper rock receives less.

**The data tested against:**

Physical Mars meteorites — real Martian rock recovered
on Earth — from the nakhlite group. Five meteorites from
the same original Mars rock pile at three known depth levels.

```
Primary source:
  Changela & Bridges (2010)
  Meteoritics & Planetary Science 45:1847-1867
  DOI: 10.1111/j.1945-5100.2010.01123.x

Secondary sources:
  Wang et al. (2021) Earth, Planets and Space
  DOI: 10.1186/s40623-021-01492-3

  Mikouchi et al. (2006) NASA Technical Reports
```

---

## DIAGNOSTIC FINDING (SCRIPT 0 — PRE-ANALYSIS)

The diagnostic script identified that **vein FeO
concentration** is the wrong primary metric. It is
dilution-sensitive. A shallow rock with very limited
water access shows high FeO concentration per unit
volume of vein, but near-zero total alteration volume.
The total iron deposited is small.

The correct metrics are:

```
PRIMARY:   Total Fe oxide index
           = alteration_volume_fraction × vein_FeO
           Proxy for total iron oxide mass deposited
           per unit rock volume.

SECONDARY: Alteration volume fraction
           Total volume of water-rock iron chemistry.

REPORTED:  Vein FeO concentration
           Expected: NO correlation.
           Confirmed: r = -0.21.
           This is a feature, not a bug.
           It confirms the dilution interpretation.
```

This identification was made **before running the full
analysis**. It is part of the record.

---

## DATA USED

```
Meteorite        Depth(m)   VeinFeO(wt%)  AltVol(frac)  FeIndex
──────────────────────────────────────────────────────────────
Y-000593           3.5         64.7         0.0010        0.0647
Nakhla             8.5         51.7         0.0050        0.2585
Gov.Val.           8.5         50.0         0.0040        0.2000
Lafayette         30.0         57.5         0.0150        0.8625
NWA998            30.0         52.0         0.0120        0.6240
──────────────────────────────────────────────────────────────
n = 5 meteorites
n_independent depth levels = 3
(Nakhla/GV paired at 8.5m;
 Lafayette/NWA998 paired at 30.0m)
```

**Depth uncertainty ranges (published):**

```
Y-000593:   0.5 – 7.0 m
Nakhla:     7.0 – 10.0 m
Gov.Val.:   7.0 – 10.0 m
Lafayette:  20.0 – 40.0 m
NWA998:     20.0 – 40.0 m
```

---

## SECTION 2 RESULTS — CORRELATION ANALYSIS

```
⚠ SAMPLE SIZE: n=5 meteorites, n_independent=3.
  No p-value from this dataset is statistically
  significant in the conventional sense.
  All results are DIRECTIONAL EVIDENCE ONLY.
```

### Panel A — Vein FeO Concentration vs Depth

```
Pearson  r  = -0.2135   (p = 0.7303)
Spearman ρ  = -0.1054   (p = 0.8660)
Regression  slope = -0.100  intercept = +56.791
```

**Interpretation:** No correlation. Expected. The shallow
Yamato 000593 has the highest FeO concentration (64.7 wt%)
but near-zero alteration volume (0.1%). This confirms the
dilution interpretation identified in the diagnostic: vein
FeO concentration is controlled by water/rock ratio, not
by total iron cycling activity.

---

### Panel B — Alteration Volume Fraction vs Depth

```
Pearson  r  = +0.9775   (p = 0.0040)
Spearman ρ  = +0.9487   (p = 0.0138)
Regression  slope = +0.000445  intercept = +0.0002
```

**Interpretation:** Very strong positive correlation.
Deeper Mars rock contains dramatically more total volume
of secondary iron chemistry. The effect is monotonic
across all three depth levels.

---

### Panel C — Total Fe Oxide Index vs Depth (PRIMARY)

```
Pearson  r  = +0.9635   (p = 0.0083)
Spearman ρ  = +0.9487   (p = 0.0138)
Regression  slope = +0.024765  intercept = +0.0032
```

**Interpretation:** Very strong positive correlation.
Total iron oxide deposited per unit rock volume increases
monotonically with depth in the Martian rock pile.

**Effect magnitude:**

```
Y-000593  (3.5m)    Fe index = 0.065
Nakhla    (8.5m)    Fe index = 0.259   →  4×  surface value
Lafayette (30.0m)   Fe index = 0.863   →  13× surface value
```

The deepest available Mars rock contains 13 times more
total iron oxide deposition than the shallowest.

---

## SECTION 3 RESULTS — MONTE CARLO UNCERTAINTY PROPAGATION

```
N = 100,000 random samples
Method: uniform sampling from published [min, max]
        ranges for both depth and measurement values.
Question: does the correlation direction hold under
          ALL plausible uncertainty combinations?
```

### Vein FeO Concentration

```
P(r > 0) = 0.2035   → NOT ROBUST
Median r  = -0.1645
95% CI    = [-0.5508, +0.2129]
```

Confirmed not robust. The correlation direction is
undefined under uncertainty. This is the expected result
for the wrong metric.

---

### Alteration Volume Fraction

```
P(r > 0) = 1.0000   → ROBUST POSITIVE
Median r  = +0.9249
95% CI    = [+0.6399, +0.9948]
```

The positive correlation holds in **100% of 100,000
random samples** across all plausible depth and
measurement combinations. The entire 95% confidence
interval is positive. The correlation never reversed
under any plausible reading of the published uncertainty
ranges.

---

### Total Fe Oxide Index (PRIMARY TEST)

```
P(r > 0) = 1.0000   → ROBUST POSITIVE
Median r  = +0.9127
95% CI    = [+0.5722, +0.9938]
```

The positive correlation holds in **100% of 100,000
random samples**. The entire confidence interval is
positive. Under no plausible combination of published
uncertainty values does the correlation reverse.

---

## SECTION 4 — FIGURES GENERATED

```
necromass_analysis_fig1_depth_vs_iron.png
  3-panel scatter plot: depth vs vein FeO,
  alteration volume, and total Fe oxide index.
  Error bars from published uncertainty ranges.
  Regression lines. Sample size warnings.

necromass_analysis_fig2_monte_carlo.png
  Monte Carlo distribution of Pearson r for
  each metric. 100,000 samples. Median and
  95% CI annotated. P(r>0) annotated.

necromass_analysis_fig3_carbon_isotopes.png
  δ¹³C scale: Martian CO₂, serpentinization,
  Nakhla organic, Earth biological,
  Curiosity SAM -137‰ anomaly.
  Biological fractionation direction annotated.

necromass_analysis_fig4_ALH84001_criteria.png
  Thomas-Keprta six-criteria table.
  Present/absent and abiotic-explained/not
  for each criterion in ALH84001 magnetite.
```

---

## SECTION 5 — HONEST INTERPRETATION

### What This Result Is

Against published compositional data from five Martian
meteorites representing three depth levels in the same
original Mars rock pile, the total iron oxide deposition
index increases monotonically with depth. The correlation
is robust across 100% of 100,000 Monte Carlo samples of
published uncertainty ranges. The effect magnitude is
13× from the shallowest to deepest available sample.

### What This Result Is Not

```
IS NOT: Proof that Mars is biologically active.
IS NOT: Proof that the iron oxide is necromass.
IS NOT: Statistically significant (n=5, n_ind=3).
IS NOT: Distinguishable from the abiotic
        temperature-gradient alternative.
```

### What This Result Is

```
IS: Directional evidence consistent with the
    biological iron-cycling prediction.

IS: Inconsistent with the simplest abiotic
    null hypothesis (top-down surface weathering
    delivering oxidants from surface downward).

IS: Robust — the direction does not reverse
    under any plausible published uncertainty
    combination across 100,000 Monte Carlo runs.

IS: The first quantitative test of the Mars Iron
    Necromass Hypothesis against real Martian
    material with a pre-registered prediction.
```

### The Abiotic Alternative That Cannot Be Excluded

The temperature gradient alternative predicts the same
direction as the biological hypothesis: deeper rock is
warmer, warmer conditions accelerate mineral alteration,
more alteration at depth.

This cannot be excluded from the current dataset.

**The distinguishing test is:**

```
1. Crystal morphology of iron oxide phases
   in Lafayette alteration veins.
   Biological: THO magnetite
               (Thomas-Keprta criteria).
   Abiotic:    Amorphous ferrihydrite/goethite.
   Test: TEM analysis of Lafayette iron oxide.
   Data availability: partially published.

2. Rate calculation:
   Total iron oxide mass in Lafayette veins
   vs. abiotic oxidation rate at estimated
   temperature and water chemistry conditions.
   If required time >> aqueous event duration:
   biological catalysis required.
   This is Script 2.

3. Carbon co-precipitation:
   δ¹³C of organic carbon co-precipitated
   with iron oxide in Lafayette veins.
   If biologically fractionated: biological signal.
   Published data: Steele et al. 2012 on Nakhla.
   Lafayette equivalent: not yet searched.
   This is Script 3.
```

---

## ADDITIONAL DATASET RESULTS

### ALH84001 Magnetite — Thomas-Keprta Criteria

```
Criteria present in ALH84001:          5/6
Criteria explained abiotically:        3/6
Criteria NOT explained abiotically:    2/6

Unexplained criteria:
  1. Truncated hexa-octahedral (THO) morphology.
     Not reproduced abiotically in published
     experiments.
  2. Chemical purity (no Ti, Mn, Cr impurities).
     Not fully explained abiotically under
     Mars conditions.

Current status: UNRESOLVED.
Neither Thomas-Keprta et al. (2001, 2009)
nor the abiotic counter-position
(Golden et al. 2004) has been retracted.
The debate is live.
```

The necromass hypothesis predicts biological magnetite
morphology in Martian iron oxide phases. The ALH84001
data shows 2/6 criteria for biological origin cannot
be explained abiotically. This is consistent with but
does not confirm the prediction.

---

### Nakhla Carbon Isotopes — Steele et al. 2012

```
Nakhla indigenous organic δ¹³C: -22 to +2‰
Martian atmospheric CO₂:        +41 to +45‰
Minimum depletion vs. CO₂:      ~39‰ lighter
Steele 2012 conclusion:         INDIGENOUS TO MARS.
                                Not contamination.

Curiosity SAM δ¹³C (mudstone):  -137‰
Depletion vs. Martian CO₂:      ~178‰
House et al. 2022 status:       THREE EXPLANATIONS.
                                NONE ELIMINATED.
                                Biological is one.
```

The biological fractionation direction (organic carbon
lighter than the Martian carbon reservoir) is present
in both Nakhla indigenous organic carbon and the
Curiosity SAM measurement. Neither confirms biology.
Both are consistent with it.

---

## WHAT COMES NEXT

### Script 2 — Iron Oxide Mass Budget and Rate Calculation

```
Purpose: Calculate total iron oxide mass in the
         nakhlite pile and compare against abiotic
         oxidation rate at estimated temperatures.
         Ask: does the rate require biological
         catalysis?

Data needed:
  Lafayette alteration vein volume (published).
  Iron oxide density in alteration products.
  Estimated temperature during aqueous event.
  Duration of aqueous event (published).
  Abiotic Fe²⁺ oxidation rate at that temperature.

Expected output:
  Required time for abiotic oxidation vs.
  actual duration of aqueous event.
  If required time >> event duration:
  biological catalysis required for the
  observed iron oxide mass.
```

### Script 3 — Lafayette Carbon Co-precipitation Search

```
Purpose: Search published Lafayette TEM/SEM/SIMS
         data for organic carbon co-precipitated
         with iron oxide in alteration veins.
         Retrieve δ¹³C values if available.

Data needed:
  Lafayette alteration vein chemistry.
  Published NanoSIMS or SIMS data for Lafayette.
  Carbon isotope values if measured.
```

---

## SUMMARY STATEMENT

On 2026-03-13, the Mars Iron Necromass Hypothesis was
tested for the first time against physical Martian
material using published compositional data from five
nakhlite meteorites representing three depth levels in
the same original Mars rock pile.

The primary result:

> Total iron oxide deposition increases monotonically
> with depth in the Martian rock pile. P(r > 0) = 1.000
> across 100,000 Monte Carlo samples. Median r = +0.913.
> 95% CI [+0.572, +0.994]. Effect magnitude: 13× from
> surface to 30m depth. This is robust directional
> evidence consistent with the biological iron-cycling
> prediction of the Mars Iron Necromass Hypothesis.

The abiotic temperature-gradient alternative predicts
the same direction and cannot be excluded from this
dataset alone. The distinguishing test is iron oxide
crystal morphology and rate calculation. Both are the
subject of Script 2 and Script 3.

The result does not prove the hypothesis.
The result does not contradict the hypothesis.
The result is the first quantitative directional
confirmation from real Martian material.

The record is here.
The timestamp is 2026-03-13.
The analysis preceded any knowledge of this
specific result.

---

## DOCUMENT METADATA

```
Document ID:       NECROMASS_ANALYSIS_RESULTS_v1.0
Hypothesis:        Mars Iron Necromass Hypothesis
Test:              Nakhlite Fe oxide vs. depth
Result:            P(r>0) = 1.000 (Monte Carlo)
                   Median r = +0.913
                   95% CI [+0.572, +0.994]
                   n=5, n_independent=3
                   DIRECTIONAL EVIDENCE ONLY
Primary source:    Changela & Bridges 2010
                   DOI: 10.1111/j.1945-5100.2010.01123.x
Pre-reg chain:     10.5281/zenodo.18986790
Repository:        github.com/Eric-Robert-Lawson/
                   attractor-oncology
ORCID:             0009-0002-0414-6544
Date:              2026-03-13
Next steps:        Script 2 — Rate calculation
                   Script 3 — Carbon co-precipitation
```
