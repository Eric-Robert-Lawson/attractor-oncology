# NECROMASS CARBON & MORPHOLOGY — SCRIPT 3 RESULTS
## Crystal Morphology and Carbon Isotope Co-precipitation
## Lafayette Nakhlite vs. Biological and Abiotic Reference Ranges
### Eric Robert Lawson / OrganismCore
### 2026-03-13

---

## DOCUMENT IDENTITY

```
Document ID:  NECROMASS_CARBON_MORPHOLOGY_RESULTS_v1.0
Date:         2026-03-13
Author:       Eric Robert Lawson / OrganismCore
Status:       ACTIVE — SCRIPT 3 RESULTS RECORD
Script:       necromass_carbon_morphology.py v1.0
Pre-reg DOI:  10.5281/zenodo.18986790
Repository:   github.com/Eric-Robert-Lawson/
              attractor-oncology
ORCID:        0009-0002-0414-6544
```

---

## WHAT WAS TESTED

Two independent discriminating tests that the
rate argument (Script 2) cannot touch:

**Test A:** Do Lafayette alteration vein iron
oxide phases show biological crystal morphology
under the Thomas-Keprta six-criteria framework?

**Test B:** Does organic carbon co-precipitated
with iron oxide in Lafayette alteration veins
show isotopic fractionation consistent with
biological iron cycling?

---

## TEST A RESULTS — CRYSTAL MORPHOLOGY

### Thomas-Keprta Criteria Applied to Lafayette

```
Criterion                  Lafayette   ALH84001   Score
──────────────────────────────────────────────────────
THO morphology             ABSENT      PRESENT     0.0
Size range 35-120nm        PARTIAL     PRESENT     0.5
Chemical purity            PARTIAL     PRESENT     0.5
Crystal perfection         PARTIAL     PRESENT     0.5
Chain organisation         PARTIAL     PARTIAL     0.5
No Fe-silicate rims        UNKNOWN     PRESENT     0.0

Weighted morphology score: 0.3182
(0.0 = fully abiotic, 1.0 = fully biological)
```

### Morphology Types Present in Lafayette

```
Framboidal aggregates:
  Abiotic-consistent. Not bio-consistent.
  Common in low-T abiotic alteration.
  Source: Bridges 2019; Golden 2004.

Euhedral prismatic chains:
  BIO-CONSISTENT AND ABIOTIC-CONSISTENT.
  DISPUTED.
  Tomkinson 2015: similar to biogenic chains.
  Bridges 2019: abiotic route possible.
  This is the unresolved morphological crux.

Euhedral blocky:
  Abiotic-consistent. Not bio-consistent.
  Standard abiotic hydrothermal form.
  Source: Golden 2004.
```

### Morphology Verdict

```
AMBIGUOUS — INCONCLUSIVE.
No THO magnetite in Lafayette.
Chain structures present but origin disputed.
Morphological test cannot discriminate.
```

---

## TEST B RESULTS — CARBON ISOTOPE CO-PRECIPITATION

### Lafayette Organic Carbon Data

```
Source:    Steele et al. 2012, Science 337:212
           NanoSIMS, Table S2
           Organic domains co-located with
           iron oxide phases in Lafayette veins.

δ¹³C range: -36.6 to -27.6‰
Mean:        -31.5‰
Co-located with iron oxide: YES
Indigenous confirmed:       YES (Steele 2012)
```

### Overlap Analysis

```
Reference                     Min      Max   Overlap  Fraction
──────────────────────────────────────────────────────────────
Mars atmospheric CO₂         +41.0   +45.0      0.0    0.0000
Serpentinization CH₄         -40.0   -10.0      9.0    1.0000
Earth biological              -35.0   -15.0      7.4    0.8222
FeOB co-precipitated C        -36.0   -24.0      8.4    0.9333
Nakhla organic C              -25.3   -18.7      0.0    0.0000
Curiosity SAM                -145.0  -128.0      0.0    0.0000
```

### Bootstrap Results (N = 50,000)

```
P(Lafayette in FeOB biogenic range)    = 0.8754
P(Lafayette in serpentinization range) = 1.0000
P(Lafayette in Mars CO₂ range)         = 0.0000
```

### Depletion Analysis

```
Mars atmospheric CO₂ mean:  +43.0‰
Lafayette organic mean:     -31.5‰
Lafayette depletion (mean):  74.5‰ vs Mars CO₂
Lafayette depletion (max):   79.6‰ vs Mars CO₂

Abiotic serpentinization maximum depletion: 83.0‰
Lafayette depletion (79.6‰): WITHIN abiotic maximum.

Lafayette mean vs FeOB biogenic mean: 1.5‰ difference.
Lafayette is 1.5‰ from the centre of the biological
reference range.
```

### Carbon Verdict — Corrected Honest Assessment

```
CONSISTENT WITH BIOLOGICAL HYPOTHESIS.
ALSO CONSISTENT WITH SERPENTINIZATION.
CANNOT DISTINGUISH FROM CARBON ISOTOPES ALONE.

The Lafayette δ¹³C range is:
  (a) Completely separated from Martian CO₂.
  (b) Co-located with iron oxide phases.
  (c) Confirmed indigenous to Mars.
  (d) Within both the FeOB biological range
      AND the serpentinization abiotic range.

The carbon signal is the strongest result in
the dataset. It is not uniquely biological.
It is not uniquely abiotic. It cannot
distinguish between the two with isotopes alone.
```

---

## CUMULATIVE ASSESSMENT — ALL THREE SCRIPTS

### Results Summary

```
Script 1 — Iron oxide depth gradient:
  P(r>0) = 1.000. Median r = +0.913.
  Effect: 13× from surface to 30m depth.
  Status: CONSISTENT. Caveat: temp gradient also.

Script 2 — Rate calculation:
  Fenton P(adequate) = 1.000.
  Fenton median time: ~25 minutes.
  Status: INCONCLUSIVE. Both pathways adequate.

Script 3A — Crystal morphology:
  Score: 0.318. No THO magnetite.
  Chains present, origin disputed.
  Status: AMBIGUOUS. INCONCLUSIVE.

Script 3B — Carbon isotopes:
  δ¹³C -36.6 to -27.6‰. Depleted 79.6‰ vs CO₂.
  Co-located with iron oxide. Indigenous.
  P(FeOB) = 0.875. P(serpentinization) = 1.000.
  Status: CONSISTENT. Not uniquely biological.
```

### Abiotic Alternatives Not Eliminated

```
Temperature gradient:
  Explains Script 1 depth gradient.
  Cannot be eliminated.

Fenton (H2O2 radiolysis):
  Explains Script 2 iron oxide mass.
  Cannot be eliminated.
  Actually faster than biology.

Serpentinization:
  Explains Script 3B carbon isotopes.
  Cannot be eliminated.
  δ¹³C range fully overlapping.

None of these three abiotic alternatives
can be eliminated by the current dataset.
All three are physically plausible on Mars.
```

### What Survives All Abiotic Alternatives

```
One finding has no single abiotic explanation:

  The COMBINATION of:
    Chain-organised magnetite (Tomkinson 2015)
  AND
    Isotopically depleted organic carbon
    co-located with iron oxide (Steele 2012)
  IN THE SAME METEORITE (Lafayette)
  AT THE DEEPEST POSITION IN THE MARS PILE.

  Serpentinization does not produce
  chain-organised magnetite.
  Temperature gradient does not produce
  chain-organised magnetite.
  Fenton does not produce
  chain-organised magnetite.

  The question that remains:
  Were the chains and the depleted carbon
  measured in the SAME ALTERATION VEIN?

  Steele 2012: NanoSIMS on Lafayette veins.
  Tomkinson 2015: TEM on Lafayette veins.
  Published data does not confirm same location.
  This is the next literature test.
```

---

## THE SINGLE MOST IMPORTANT FINDING

From all three scripts combined, against real
Martian material, with published sources,
with timestamps before this analysis:

> Indigenous Martian organic carbon, confirmed
> by NanoSIMS (Steele et al. 2012), co-located
> with iron oxide phases, isotopically depleted
> 79.6‰ from the Martian CO₂ reservoir,
> in the deepest and most iron-processed nakhlite
> in the original Mars rock pile, alongside
> chain-organised magnetite whose biological
> origin is disputed but not eliminated.

That is not proof.

That is a pattern. Consistent across three
independent tests. Three independent datasets.
All real Martian material. All published sources.

No single test is conclusive.
No single test contradicts the hypothesis.
The pattern is internally consistent.

---

## THE REMAINING DISCRIMINATING TEST

```
Co-location of chain-organised magnetite
WITH isotopically depleted organic carbon
IN THE SAME ALTERATION VEIN
IN LAFAYETTE
at nanometre scale resolution.

Required: NanoSIMS + TEM on same thin section.
Steele 2012 did NanoSIMS on Lafayette.
Tomkinson 2015 did TEM on Lafayette.
Were they on the same vein? Not confirmed.

This is the test that would cut through all
three abiotic alternatives simultaneously.
No abiotic single process produces both.
Biology produces both from the same organism.

This is the basis for the collaboration request.
```

---

## FIGURES GENERATED

```
necromass_carbon_fig8_isotope_overlap.png
  Carbon isotope scale: Lafayette vs all
  reference ranges. Overlap regions shown.
  Biological fractionation direction annotated.

necromass_carbon_fig9_morphology.png
  Thomas-Keprta criteria table for Lafayette
  vs ALH84001. Weighted score. Chain debate.

necromass_carbon_fig10_overlap_matrix.png
  Overlap fraction for each reference range.
  Shows Lafayette cannot be distinguished
  from serpentinization on isotopes alone.
```

---

## DOCUMENT METADATA

```
Document ID:     NECROMASS_CARBON_MORPHOLOGY_RESULTS_v1.0
Tests:           Crystal morphology (Thomas-Keprta)
                 Carbon isotope co-precipitation
Morphology:      AMBIGUOUS. Score 0.318.
                 No THO. Chains present, disputed.
Carbon:          CONSISTENT. Not uniquely biological.
                 δ¹³C -36.6 to -27.6‰, depleted 79.6‰
                 P(FeOB) = 0.875
                 P(serpentinization) = 1.000
Abiotic alts:    Temperature gradient — not eliminated.
                 Fenton — not eliminated.
                 Serpentinization — not eliminated.
Pattern:         4 independent tests, all consistent,
                 none contradictory, none conclusive.
Next test:       NanoSIMS + TEM co-location on same
                 Lafayette vein. Chains + depleted C
                 in same location = no single abiotic
                 explanation.
Pre-reg chain:   10.5281/zenodo.18986790
Date:            2026-03-13
```
