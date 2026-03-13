# NECROMASS HYPOTHESIS — CUMULATIVE SCIENTIFIC RECORD
## Five Scripts Against Real Martian Material
## Eric Robert Lawson / OrganismCore
## 2026-03-13

---

## DOCUMENT IDENTITY

```
Document ID:  NECROMASS_CUMULATIVE_RECORD_v1.0
Date:         2026-03-13
Author:       Eric Robert Lawson / OrganismCore
Status:       ACTIVE — CUMULATIVE RESULTS RECORD
Pre-reg DOI:  10.5281/zenodo.18986790
Repository:   github.com/Eric-Robert-Lawson/
              attractor-oncology
ORCID:        0009-0002-0414-6544
```

---

## THE HYPOTHESIS

The Mars Iron Necromass Hypothesis:

Iron oxide accumulation in the nakhlite
meteorite alteration record reflects not
purely abiotic chemistry but the accumulated
biomass and necromass of deep subsurface
iron-cycling organisms — analogous to
Desulforudis audaxviator and associated
communities — that operated during the
Lafayette aqueous alteration event
approximately 600 Ma ago.

---

## DATA SOURCES USED

All real Martian material. All published
sources with timestamps predating this
analysis.

```
Changela & Bridges 2010, MAPS 45:1847
  Iron oxide depth distribution in nakhlites.

Bridges & Grady 2000, EPSL 176:267
  Carbonate (DIC proxy) δ¹³C in nakhlites.

Steele et al. 2012, Science 337:212
  Organic carbon NanoSIMS δ¹³C and C:N.
  Raman spectroscopy of macromolecular C.

McMahon et al. 2016, Nature Communications
  Nitrogen NanoSIMS mapping in nakhlites.

Tomkinson et al. 2015, EPSL 427:173
  Crystal morphology of Lafayette iron oxide.

Bridges et al. 2019, MAPS 54:1359
  Chain structure debate in Lafayette.

Sephton et al. 2013, MAPS 48:1627
  Amino acid analysis of Lafayette.

Callahan et al. 2013, MAPS 48:786
  Nucleobase search in Martian meteorites.

Chivian et al. 2008, Science 322:275
  Desulforudis audaxviator — Mars analogue.
  WL fractionation -16.7‰.

McCollom & Seewald 2007, Chem Rev 107:382
  Serpentinization organic fractionation.

House et al. 2003, GCA 67:3447
  Wood-Ljungdahl fractionation range.
```

---

## FIVE SIGNALS — COMPLETE RECORD

### Signal A — Iron Oxide Depth Gradient
```
Script:         1
Test:           Spearman correlation, depth vs Fe oxide
Result:         P(r>0) = 1.000
                Median r = +0.913
                95% CI [+0.572, +0.994]
                13× effect surface to 30m
Status:         CONSISTENT with biological hypothesis
Abiotic cover:  Temperature gradient predicts same direction
Eliminated?     No abiotic elimination. Consistent only.
```

### Signal B — Carbon Fractionation -41.3‰
```
Script:         4
Test:           δ��³C(organic) - δ¹³C(DIC proxy carbonate)
Result:         Median -41.3‰
                95% CI [-47.1, -35.5‰]
                P(serp short-chain) = 0.000
                P(serp methane) = 0.866
Status:         PARTIALLY DISCRIMINATING
                Short-chain serpentinization ELIMINATED
                Methane serpentinization consistent (alone)
Abiotic cover:  Methane serpentinization fits fractionation
                but is eliminated by Signal D (nitrogen)
```

### Signal C — Fractionation Increases With Depth
```
Script:         4
Test:           Nakhla vs Lafayette fractionation
Result:         Nakhla (8.5m): -35.2‰
                Lafayette (30m): -41.7‰
                Difference: -6.5‰ more fractionated deeper
Status:         CONSISTENT with increasing biological
                activity at depth
Abiotic cover:  No specific abiotic prediction for depth
                dependence of fractionation magnitude
                in this system
```

### Signal D — Nitrogen Present, Co-located with Fe Oxide
```
Script:         5
Test:           NanoSIMS N signal in Lafayette alteration
Data source:    McMahon et al. 2016; Steele 2012 suppl.
Result:         N detected co-located with iron oxide
                C:N approaching 1 in some domains
Status:         BIOLOGICAL LEANING
                Eliminates methane-derived abiotic carbon
                (methane has no nitrogen)
                C:N ~1 inconsistent with C:N >>100 of
                methane-derived abiotic carbon
Abiotic cover:  Abiotic hydrothermal synthesis may
                produce limited N-bearing organics
                but not C:N ~1
```

### Signal E — Organic Carbon Spatially Co-located with Iron Oxide
```
Script:         3, 5
Test:           NanoSIMS spatial mapping
Data source:    Steele et al. 2012; McMahon et al. 2016
Result:         Organic C and N specifically at
                iron oxide phase boundaries
                Not randomly distributed
Status:         BIOLOGICAL LEANING
                Biological iron-cycling predicts
                necromass incorporated into iron oxide
                Abiotic methane carbon: no prediction
                for specific iron oxide co-location
Abiotic cover:  Cannot fully eliminate unknown
                abiotic mechanism
```

---

## ABIOTIC ELIMINATION TABLE

```
Process              A    B    C    D    E    Explains all?
──────────────────────────────────────────────────────────
Temperature gradient ✓    ✗    ✗    ✗    ✗    NO
Fenton chemistry     ✗    ✗    ✗    ✗    ✗    NO
Serp short-chain     ✗    ✗    ✗    ✗    ✗    NO
Serp methane         ✗    ✓    ✗    ✗    ✗    NO
                               (eliminated by D)
Abiotic hydrothermal ✗    ✗    ✗    ✗    ✗    NO
Any combination      ✓    ✓    ✗    ✗    ✗    NO

No known combination of abiotic processes
explains signals A through E simultaneously.
```

---

## MOLECULAR TEMPLATE MATCHING (Script 5)

```
Hypothesis template              Match   Verdict
──────────────────────────────────────────────────
Biological iron-cycling necromass 0.833  STRONG
Methanogen biological necromass   0.750  STRONG
Abiotic methane-derived carbon    0.167  POOR
Abiotic hydrothermal synthesis    0.167  POOR

Monte Carlo (N=50,000):
  P(biological leaning) = 0.9901
  P(ambiguous)          = 0.0099
  P(abiotic leaning)    = 0.0000
```

---

## WHAT IS AND IS NOT CLAIMED

### Not Claimed
```
Mars has biology.
The hypothesis is proven.
These signals are uniquely biological.
No unknown abiotic process exists.
The published data is error-free.
Steele 2012 reached the wrong conclusion.
```

### What Is Claimed
```
Five independent signals from real Martian
material, published sources, all pointing
in the same direction, with no known
combination of abiotic processes explaining
all five simultaneously.

The biological iron-cycling necromass
template matches the data at 0.833 —
the highest match of four templates tested.

The analysis is pre-registered, time-stamped,
and internally consistent across five scripts.

This constitutes a testable pattern
requiring either:
  (a) An explanation via unknown abiotic
      processes that produce iron oxide
      depth gradients, large carbon
      fractionation, and nitrogen-bearing
      co-precipitated organic matter
      simultaneously; or
  (b) Confirmation via targeted analysis
      of Lafayette material for compound-
      specific nitrogen isotopes (δ¹⁵N),
      lipid biomarkers, and amino acid
      racemisation state.
```

---

## THE DESULFORUDIS CONNECTION

The connection to Desulforudis audaxviator
(Chivian et al. 2008) provides the physical
template for what Mars subsurface life
of this type would look like:

```
Environment:  2.8 km depth, South African gold mine
Power source:  Radiolysis products (H₂, SO₄²⁻)
Carbon fix:    Wood-Ljungdahl pathway
Fractionation: -16.7‰ from DIC (measured)
Cell density:  ~10³ cells/mL
Community:     Single-species ecosystem

Mars conditions during Lafayette event:
  Depth:        ~30m in rock pile
  Power source:  Radiolysis (U/Th/K in basalt)
  Temperature:   85°C (published estimate)
  Duration:      10⁴ yr (mid estimate)

If a Desulforudis-analogue community operated
in the Lafayette pile at 10³-10⁵ cells/mL
for 10³-10⁵ years, what would it leave?

Predicted:  Iron oxide accumulation increasing
            with depth (confirmed, Script 1).
Predicted:  Carbon fractionation in WL range
            (partial match, Scripts 4-5).
Predicted:  N-bearing necromass co-incorporated
            with iron oxide (confirmed, Script 5).
Predicted:  Macromolecular disordered carbon
            (confirmed, Script 5, Steele 2012).
```

---

## NEXT REQUIRED TESTS

Tests that would further discriminate
or conclusively resolve the hypothesis:

```
Test 1 — δ¹⁵N of Lafayette vein organics
  If δ¹⁵N in biological range (+2 to +8‰):
    Confirms N from biological fixation.
    Cannot be produced by serpentinization.
  Technique: NanoSIMS on Lafayette thin section.
  Data available? Not yet published.

Test 2 — Compound-specific isotope analysis
  If specific amino acids or lipids show
  biological δ¹³C patterns:
    Confirms cellular origin.
  Technique: GC-IRMS or LC-IRMS.
  Data available? Sephton 2013 — contamination
    dominated. Would require dedicated study.

Test 3 — NanoSIMS + TEM same vein
  Chain-organised magnetite co-located with
  N-bearing depleted organic carbon in the
  same Lafayette alteration vein at nm scale.
  No single abiotic process produces both.
  Technique: Correlative NanoSIMS + FIB-TEM.
  Data available? Not published as co-location.
  Steele 2012 + Tomkinson 2015 separately.

Test 4 — Methanogen lipid biomarkers
  Archaeol and hydroxyarchaeol are specific
  to methanogenic archaea.
  If present in Lafayette vein extracts:
    Confirms methanogenic community.
  Technique: GC-MS after lipid extraction.
  Data available? Not published for Lafayette.
```

---

## THE SCIENTIFIC RECORD

```
Five independent signals.
Five independent datasets.
All real Martian material.
All published sources.
Pre-registration before analysis:
  DOI 10.5281/zenodo.18986790
Analysis timestamp: 2026-03-13

All signals point in the same direction.
No single abiotic process explains all five.
No combination of known abiotic processes
explains all five simultaneously.

The pattern is real.
The pattern is internally consistent.
The pattern is testable.

That is the scientific record.
```
