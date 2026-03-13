# NECROMASS RATE CALCULATION — SCRIPT 2 RESULTS
## Iron Oxide Mass Budget and Rate Analysis
## Lafayette Nakhlite vs. Published Aqueous Event Parameters
### Eric Robert Lawson / OrganismCore
### 2026-03-13

---

## DOCUMENT IDENTITY

```
Document ID:  NECROMASS_RATE_RESULTS_v1.0
Date:         2026-03-13
Author:       Eric Robert Lawson / OrganismCore
Status:       ACTIVE — SCRIPT 2 RESULTS RECORD
Script:       necromass_rate_calculation.py v1.0
Pre-reg DOI:  10.5281/zenodo.18986790
Repository:   github.com/Eric-Robert-Lawson/
              attractor-oncology
ORCID:        0009-0002-0414-6544
```

---

## QUESTION ASKED

Given the total mass of iron oxide in Lafayette
nakhlite alteration veins, and given the published
duration and temperature of the single aqueous event
that produced it, could abiotic chemistry have
deposited that iron oxide in the time available?
Or does the mass require biological catalysis?

---

## PARAMETERS USED

```
Lafayette alteration:
  Alt. volume fraction : 0.015 [0.010–0.025]
  Product density      : 2.3 g/cm³
  FeO wt fraction      : 0.35 [0.32–0.575]
  Event temperature    : 85°C [50–150°C]
  Event duration       : 1e4 yr [10–1e6 yr]
  Fluid:rock ratio     : 0.05 [0.01–0.10]

Sources: Changela & Bridges 2010; Treiman 1993;
         Bridges et al. 2001; Wang et al. 2021

Fe oxide mass computed:
  12.07 kg per m³ rock (mid values)
  244.2 mol Fe per m³ rock (MC median)
  150 L water per m³ rock (mid fluid:rock ratio)
```

---

## RESULTS

### Point Estimate (Mid Values)

```
Pathway               Rate (mol/L/yr)   Time req.   Adequate?
─────────────────────────────────────────────────────────────
Fenton (abiotic)         5.762e+05     1.94e-06 yr     YES
Surface-cat (abiotic)    1.199e-01     9.35e+00 yr     YES
Biological               9.468e-03     1.18e+02 yr     YES

Event duration (mid): 1e+04 yr
```

### Monte Carlo (N = 100,000 samples)

```
Fenton (abiotic):
  P(adequate) = 1.0000
  Median time = 4.72e-05 yr (~25 minutes)
  95% CI      = [2.74e-06, 1.98e-03] yr

Surface-cat (abiotic):
  P(adequate) = 0.9965
  Median time = 8.58e+02 yr
  95% CI      = [4.00e+01, 2.17e+04] yr

Biological:
  P(adequate) = 1.0000
  Median time = 3.83 yr
  95% CI      = [0.64, 93.9] yr

Fenton / Biological time ratio:
  Median      = 1.16e-05×
  95% CI      = [2.73e-07×, 5.21e-04×]
  Interpretation: Fenton is ~86,000× FASTER
  than biological at these conditions.
  P(Fenton slower than bio) = 0.0000
```

---

## VERDICT

```
ABIOTIC FENTON PATHWAY SUFFICIENT IN ALL
PARAMETER COMBINATIONS.

The Fenton (H2O2 radiolysis) pathway can
deposit the entire observed Lafayette iron
oxide mass in approximately 25 minutes at
mid-parameter values. Biological catalysis
is not required to explain the quantity or
the rate of iron oxide deposition.

RATE CALCULATION IS INCONCLUSIVE AS A
DISCRIMINATING TEST.
Both pathways are fast enough.
The rate argument does not distinguish
biological from abiotic for this system.
```

---

## WHAT THIS MEANS

### What This Result Is Not

```
NOT: Proof that Mars has no biology.
NOT: Refutation of the necromass hypothesis.
NOT: Invalidation of Script 1 results.
NOT: Evidence against biological iron cycling.
```

### What This Result Is

```
IS: The rate argument is not a useful
    discriminating test for this system.
    The Fenton pathway is too fast.
    It can explain any reasonable iron oxide
    quantity in any reasonable event duration.

IS: An honest negative result that eliminates
    one potential line of evidence.
    Science works by eliminating tests that
    cannot distinguish hypotheses.
    This is one of those tests.

IS: A clarification of what the actual
    discriminating tests are (see below).
```

### Script 1 Stands Independently

The rate test being inconclusive does not affect
the Script 1 result:

```
Script 1 result (stands):
  Total iron oxide deposition increases
  monotonically with depth in the Martian
  rock pile. P(r>0) = 1.000. Median r = +0.913.
  95% CI [+0.572, +0.994].
  Effect magnitude: 13× from surface to 30m.
```

The rate test and the depth gradient test are
independent lines of evidence. The rate test
returning inconclusive does not affect the
depth gradient result.

---

## WHAT THE RATE TEST REVEALED

Script 2 produced an unexpected finding about
the Fenton pathway:

```
At 85°C with 0.1 mM H2O2 and 1 mM Fe2+:
  Fenton rate = 5.76e+05 mol/L/yr
  Required time = ~25 minutes

This is extraordinarily fast.

Implications:
  1. The Fenton pathway is so efficient at
     Mars-relevant temperatures that it cannot
     be used as a discriminator.
  2. Any aqueous event on Mars with U/Th/K
     in the surrounding rock and iron in
     solution will rapidly oxidise available
     iron regardless of biology.
  3. The iron oxide quantity alone is not a
     biosignature in an environment with
     active radiolysis.
  4. The DISTRIBUTION of iron oxide
     (depth gradient, Script 1) remains
     a more useful test than the quantity.
```

---

## THE ACTUAL DISCRIMINATING TESTS

Script 2 has clarified the hierarchy of tests:

```
ELIMINATED AS DISCRIMINATOR:
  ✗ Rate/quantity argument
    Fenton too fast. Inconclusive.

PARTIALLY DISCRIMINATING (Script 1):
  ✓ Iron oxide depth gradient
    P(r>0) = 1.000. Strong directional.
    Caveat: temperature gradient
    predicts same direction.

KEY DISCRIMINATING TESTS (Script 3):
  ✓ Crystal morphology of iron oxide
    in Lafayette alteration veins.
    THO magnetite = biological signal.
    Amorphous ferrihydrite = abiotic.
    This is the most direct available test.

  ✓ Carbon co-precipitation with iron oxide.
    δ¹³C of organic carbon in Lafayette veins.
    Biologically-fractionated carbon
    co-located with iron oxide = biological.
    Independent of rate argument entirely.
```

---

## FIGURES GENERATED

```
necromass_rate_fig5_time_required.png
  Time required distributions for each pathway.
  Event duration range overlaid.
  Shows all three pathways adequate.

necromass_rate_fig6_abiotic_vs_bio_ratio.png
  Fenton/Biological time ratio distribution.
  Shows Fenton is ~86,000× faster than bio.

necromass_rate_fig7_fe_mass_and_time.png
  Fe moles in Lafayette and time comparison
  bar chart vs event duration range.
```

---

## DOCUMENT METADATA

```
Document ID:      NECROMASS_RATE_RESULTS_v1.0
Test:             Iron oxide rate calculation
                  vs published event parameters
Result:           INCONCLUSIVE
                  Fenton P(adequate) = 1.000
                  Bio P(adequate) = 1.000
                  Fenton/Bio ratio = 1.16e-05×
                  (Fenton ~86,000× faster)
Interpretation:   Rate argument cannot
                  distinguish bio from abiotic.
                  Not a useful discriminator.
Script 1 status:  STANDS INDEPENDENTLY.
                  P(r>0) = 1.000 unaffected.
Next:             Script 3 — crystal morphology
                  and carbon co-precipitation
                  in Lafayette alteration veins.
Pre-reg chain:    10.5281/zenodo.18986790
Date:             2026-03-13
```
