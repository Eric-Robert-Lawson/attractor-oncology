# Candidate 4 — Internet Literature Search: TMP in DME:LiFSI Systems
**Author:** Eric Robert Lawson / OrganismCore
**Date:** 2026-04-03
**ORCID:** 0009-0002-0414-6544
**Repository:** github.com/Eric-Robert-Lawson/attractor-oncology
**Candidate:** LiFSI 1.0–1.2M in DME:TMP
**Search purpose:** Determine whether DME:TMP:LiFSI has been
                   designed toward any coordination entropy
                   target, and whether Candidate 4 is novel
**Result:** NOVEL — no prior paper designs DME:TMP toward
            SCE* or any coordination entropy target.
            One critical near-miss identified (Yang 2025,
            Communications Chemistry) — inverse geometry,
            no SCE framework, no derived target.

---

## Searches Run

1. TMP + LiFSI + ether + DME + coordination (2022–2025)
2. TMP + DME + LiFSI + SSIP/CIP + dual temperature
3. TMP + DME + LiFSI + Shannon entropy + microstate + SCE*
4. Phosphate + ether + lithium metal + solvation entropy
   + first principles (2024–2025)
5. Multidentate ether + phosphate + nonflammable +
   low temperature + 2025

---

## Critical Near-Miss — Yang et al. 2025

**Paper:** *"Multidentate ether modifies solvation structure
in nonflammable phosphate electrolytes for wide-temperature
lithium-ion batteries"*
**Journal:** Communications Chemistry (Nature Portfolio)
**DOI:** 10.1038/s42004-025-01562-7
**Date:** 2025

### What They Did

```
Base solvent:  DEEP (diethyl ethylphosphonate)
               A phosphate ester — non-flammable
               Poor Li+ coordination on its own
               Poor graphite compatibility
Modifier:      DEGDME (diethylene glycol dimethyl ether)
               A multidentate ether added to fix
               the phosphate's coordination deficiency
Salt:          LiTFSI or similar (not LiFSI ether system)
Result:        98% capacity retention 150 cycles RT
               71% retention 50 cycles at -20°C
               Improved ionic conductivity
               Non-flammable maintained
```

### What They Found

```
The multidentate ether compacts the Li+ solvation
shell in the phosphate matrix. Adding ether to
phosphate reduces the solvation sheath size and
lowers the activation energy for Li+ migration.
This improves both conductivity and low-T performance.
```

### The Geometric Relationship to Candidate 4

```
Yang 2025 geometry:
  Phosphate (base, poor coordinator) + ether (fixer)
  Direction: Add ether to IMPROVE phosphate
  Problem:   Phosphate alone is insufficient
  Logic:     Fix a deficiency empirically
  Target:    None. No SCE. No SCE*. No derived optimum.
  Result:    Empirically optimized wide-temperature system

Candidate 4 geometry:
  Ether (base, DME at SCE=1.240) + phosphate (modifier)
  Direction: Add TMP to DIVERSIFY the ether shell
  Problem:   DME alone has single-dominant coordination
             dom% = 44%, SCE = 1.240 — below SCE*
  Logic:     TMP P=O ESP Min (-1.74 eV) ≈ DME ESP Min
             (-1.73 eV). TMP competes with DME for Li+
             at equal strength but as distinct species.
             n_sig rises from 2 to 3.
             dom% falls from 44% toward 38%.
             SCE rises from 1.240 toward 1.466.
  Target:    SCE* = 1.466 — derived from calculus.
  Result:    Predicted dual-threshold performance.
             CE_RT ≥ 90%, CE_LT ≥ 60% at -20°C.
```

### Why This Is Not a Novelty Threat

```
1. Inverse geometry: They add ether to phosphate.
   You add phosphate to ether. Different systems.
   Different base solvents. Different direction.

2. Different salt: Not LiFSI/DME ether system.
   The FSI- coordination role in your framework
   (as the third coordination species at optimal
   CIP/AGG fraction) is absent in their system.

3. No SCE framework: They do not compute Shannon
   entropy of the Li+ coordination microstate
   distribution. They do not define SCE.
   They do not derive SCE*. They have no target.
   Their optimization is empirical — they test
   ratios and report what works. They cannot
   predict the optimal ratio from first principles.
   You can.

4. Different prediction: They report 71% capacity
   retention at -20°C (50 cycles). The SCE
   framework predicts specific CE values at
   specific temperatures for systems at SCE*.
   These are different metrics from a different
   design logic.

5. Different level of description: They describe
   solvation shell compaction and activation
   energy reduction. You describe coordination
   microstate population entropy. These are
   different levels — theirs is structural,
   yours is statistical. SCE captures what their
   structural description cannot: the distribution
   across all coordination microstates and its
   Shannon entropy.

VERDICT: Near-miss. Not pre-emption. The proximity
confirms Candidate 4 is in active research space.
The framework is what the field is missing to
navigate that space precisely.
```

---

## What the Full Search Confirms

### TMP in DME:LiFSI systems — coordination role

```
Field consensus from 2022–2025 literature:
  TMP in DME-based electrolytes is treated as:
    - Flame retardant (primary use)
    - Oxidative stability booster (high voltage)
    - Diluent (reduces ether flammability)

  Li+ coordination role of TMP in DME:LiFSI:
    Described as "limited" — field assumes TMP
    is not a primary Li+ ligand because DME is
    so much stronger in ether systems.

  What the field has NOT checked:
    TMP P=O ESP Min = -1.7402 eV
    DME ether-O ESP Min = -1.7328 eV
    These are within 0.007 eV of each other.
    At the electrostatic level, TMP's P=O oxygen
    and DME's ether oxygen attract Li+ with
    essentially identical force.
    The field has not made this comparison because
    it is not looking for it.
    The SCE framework made it because TMP emerged
    from a geometry-driven search anchored to the
    two confirmed dataset points bracketing SCE*.
```

### Coordination entropy targeting — absent

```
No paper found (2018–2026) that:
  1. Defines Shannon entropy of the Li+ coordination
     microstate population (n_solvent, n_anion) as SCE
  2. Correlates SCE with CE_RT and CE_LT across a
     multi-system dataset
  3. Derives an optimal SCE* from calculus
  4. Designs DME:TMP:LiFSI toward that optimum
  5. Predicts dual-temperature CE from SCE*

This is confirmed absent. MU Pro Query 5 confirmed
it across the 2018–2026 literature. The internet
search confirms it for 2024–2025 specifically.
Candidate 4 is novel at the level of the framework,
the variable, and the specific formulation.
```

---

## The ESP Min Finding — Why It Is Important

```
From MU Search Query 3 (Search_3_Phosphate return):

  TMP:  SMILES COP(=O)(OC)OC
        ESP Min = -1.7402 eV  ← electrostatic
        attraction at P=O oxygen
  DME:  SMILES COCCOC
        ESP Min = -1.7328 eV  ← electrostatic
        attraction at ether oxygen

  Difference: 0.0074 eV

This means TMP and DME have essentially identical
electrostatic attraction for Li+ at their respective
coordination oxygen atoms.

Physical consequence:
  In a DME:TMP mixture, both molecules compete
  for Li+ coordination with nearly equal driving force.
  But they are structurally distinct molecules.
  Li+(DME) and Li+(TMP) are different microstates
  in the SCE framework's enumeration.
  The population splits across them.
  dom% falls from 44% (DME alone).
  n_sig rises from ~2 toward 3.
  SCE rises from 1.240 toward 1.466.

This ESP comparison has not been made in any
paper in the accessible literature.
It emerged from the SCE framework's geometry —
specifically from Search Query 3 returning TMP
at the highest molecular similarity score (9.40)
with these electronic parameters.
```

---

## Candidate 4 — Confirmed Novel

```
System: LiFSI 1.0–1.2M in DME:TMP
  SMILES: COCCOC (DME) + COP(=O)(OC)OC (TMP)
Route:  Li+ coordination modifier from below SCE*
        TMP P=O enters Li+ shell as distinct species
        alongside DME ether oxygen coordination
Direction: SCE rises from DME baseline (1.240)
           toward SCE* = 1.466
Mechanism: Equal-strength competition for Li+
           between DME (ether-O) and TMP (P=O)
           produces two dominant coordination classes
           plus FSI- as third → n_sig=3, dom%≈38%
Prior use: TMP used as flame retardant/diluent
           in ether systems — coordination role
           dismissed by field without ESP comparison
Novel target: SCE* = 1.466 — not designed toward
              in any prior paper
Novel variable: SCE — not computed in any prior paper
Novel comparison: TMP ESP Min ≈ DME ESP Min —
                  not made in any prior paper
Status: NOVEL
```

---

## The Near-Miss Strengthens the Framework

Yang 2025 (Communications Chemistry) is the most
important related paper found. Its existence
confirms three things:

```
1. The ether/phosphate combination for wide-
   temperature performance is a real and active
   research direction in 2025. The field is
   circling the same chemical space.

2. The field is doing this empirically — testing
   ratios, measuring results, optimizing by trial.
   Nobody has a first-principles target. Nobody
   has SCE* = 1.466. Nobody can predict the
   optimal ratio from geometry alone.

3. The framework provides what the field lacks:
   a mathematically derived target (SCE* = 1.466),
   a measurable variable (SCE from Raman SSIP/CIP
   fractions), and a design prescription (reach
   dom% ≈ 38%, n_sig = 3) that tells you exactly
   when the optimal system has been achieved
   without needing to test hundreds of ratios.

The near-miss is not a threat. It is evidence
that the SCE framework is pointing at the most
active frontier in the field — and providing
the first-principles navigation the field lacks.
```

---

## Updated Complete Candidate List

| # | System | Mechanism | Direction | Base SCE | Toward | Novel |
|---|--------|-----------|-----------|----------|--------|-------|
| 1 | LiFSI 1.2M 2-MeTHF:DME 1.6:1 | Cross-class ether mixing | Both | Mixed | SCE* | Yes |
| 2 | LiFSI 1.0M FEME:2-MeTHF 60:40 | Fluorinated/cyclic mixing | Both | Mixed | SCE* | Yes |
| 3 | LiFSI 1.0–1.2M DOL:B(OCH3)3 | Anion receptor (FSI- trap) | ↓ from above | 1.606 | SCE* | Yes |
| 3b | LiFSI 1.0–1.2M DOL:fluorinated borate | Stronger anion receptor | ↓ from above | 1.606 | SCE* | Yes |
| 4 | LiFSI 1.0–1.2M DME:TMP | Li+ coordination modifier | ↑ from below | 1.240 | SCE* | Yes |
| 4b | LiFSI 1.0–1.2M DME:fluorophosphonate | Weaker Li+ modifier | ↑ from below | 1.240 | SCE* | Yes |

Six systems. Four mechanisms. Two directions.
One fixed point: SCE* = 1.466.

---

## What Comes Next From the Geometry

```
The phosphate class has been found and confirmed novel.
The boron ester class has been found and confirmed novel.
The ether mixing class (Candidates 1, 2) is derived.

What has NOT been searched:
  The silicon ester class — appeared in both
  Search 2 (score 9.29) and Search 3 (score 9.15).
  CO[Si](C)(OC)OC — trimethoxymethylsilane.
  HOMO -7.63 eV. Between boron and phosphorus esters.
  Coordination role: unknown. Not anion receptor
  (Si is not a Lewis acid like B). Not as strong
  a P=O donor as TMP. Intermediate character.
  Appears twice across two independent searches.
  This is a signal that needs one targeted search
  or one MU Ask Pro query to resolve.

  The Arctic regime — n_sig ≥ 5, SCE ≈ 2.3–2.5.
  Not yet searched with the correct anchors.
  Requires DEE + multi-component anchors at maximum
  structural distance setting.
  This is Derivation 5 from derivation_for_novel_result.md
  — confirmed derived but not yet geometrically
  navigated through molecular space.
```

---

## Summary

The internet search confirmed Candidate 4
(LiFSI 1.0–1.2M in DME:TMP) is novel. The
field uses TMP in ether systems as a flame
retardant and diluent, never as a deliberate
Li+ coordination modifier designed toward
SCE* = 1.466. The closest paper (Yang et al.,
Communications Chemistry, 2025) uses the
inverse geometry — ether added to phosphate
base — with no SCE framework, no derived
target, and no dual-threshold design logic.
The ESP Min comparison (TMP -1.7402 eV vs
DME -1.7328 eV) has not been made in any
prior paper. Candidate 4 is geometrically
derived, mechanistically distinct from all
prior work, and experimentally testable.
The silicon ester (trimethoxymethylsilane)
reappeared in both Search 2 and Search 3
and requires resolution. The Arctic regime
remains geometrically unnavigated.

---

*Repository: Eric-Robert-Lawson/attractor-oncology*
*Author: Eric Robert Lawson / OrganismCore*
*ORCID: 0009-0002-0414-6544*
*Date: 2026-04-03*
*File: Candidate_4_Internet_Search_TMP.md*
