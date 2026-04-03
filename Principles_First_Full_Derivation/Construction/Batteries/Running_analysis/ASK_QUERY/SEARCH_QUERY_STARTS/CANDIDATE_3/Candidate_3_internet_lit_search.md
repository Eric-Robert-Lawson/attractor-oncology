# Candidate 3 — Internet Literature Search Results
**Author:** Eric Robert Lawson / OrganismCore
**Date:** 2026-04-03
**ORCID:** 0009-0002-0414-6544
**Repository:** github.com/Eric-Robert-Lawson/attractor-oncology
**Purpose:** Internet search to determine what is
             published on trimethyl borate in DOL/LiFSI
             ether electrolytes and whether Candidate 3
             is novel or pre-empted
**Result:** NOVEL — prior use of TMB is for different
            reasons. One new constraint identified.
            One new candidate variant generated (3b).

---

## Search Queries Run

1. Trimethyl borate DOL dioxolane LiFSI lithium metal
   battery electrolyte additive Coulombic efficiency
2. TMB LiFSI published paper 2020–2025 CE results
3. Trimethyl borate DOL LiFSI ether electrolyte
   solvation structure SSIP CIP coordination
4. Trimethyl borate DOL low temperature -20°C CE
   2022 2023 2024
5. Trimethyl borate lithium metal battery anion
   receptor solvation entropy coordination diversity
   electrolyte design 2023 2024 2025

---

## What Was Found — Paper by Paper

### Paper 1 — Ding et al. 2023 (ACS Applied Materials & Interfaces)
*"Roles of Trimethyl Borate in Constructing an
Interphase on Li Anode: Angel or Demon?"*

```
TMB studied as: SEI-forming additive on Li anode
Primary finding: TMB contributes to interphase
                 formation. At low concentration
                 (≤0.5 wt%) beneficial. At high
                 concentration harmful (excess LiF,
                 increased impedance).
Salt system:    LiPF6-based carbonate, not LiFSI ether
SCE framing:    Not used. Not computed. Not mentioned.
Coordination
entropy target: None. No SCE* concept.
Low-T data:     Not the primary focus.
```

**Relationship to Candidate 3:**
TMB has been studied as SEI former. Not as anion
receptor in DOL:LiFSI designed toward SCE* = 1.466.
The formulation, the salt, the purpose, and the
design principle are all different. Not a novelty
threat. Confirms TMB is compatible with Li metal
interface chemistry.

### Paper 2 — Yang et al. 2021 (ACS Applied Materials & Interfaces)
*"Synergy Effect of Trimethyl Borate on Protecting
High-Voltage Cathode Materials in Dual-Additive
Electrolytes"*

```
TMB studied as: Cathode protective additive
                in carbonate electrolytes with
                high-voltage Ni-rich cathodes
Primary finding: TMB + FEC dual-additive improves
                 cathode interphase stability
Salt system:    LiPF6-based, high-voltage carbonate
SCE framing:    Not used.
Low-T data:     Not reported.
```

**Relationship to Candidate 3:**
Cathode-focused study in carbonate electrolytes.
No overlap with DOL:LiFSI ether system.
Not a novelty threat.

### Paper 3 — RSC 2024 (Chemical Science)
*"In situ polymerization of 1,3-dioxolane and
formation of fluorine/boron rich SEI"*

```
System: Borate (THB — tris(hexafluoroisopropyl)
        borate) initiates DOL ring-opening
        polymerization to form PDOL gel electrolyte
Finding: Borate Lewis acid catalyzes DOL
         polymerization. Forms boron/fluorine-rich
         SEI. Li+ transference number 0.76.
         >300 cycles stable.
Critical: INTENTIONAL polymerization of DOL using
          borate Lewis acid as initiator.
```

**Relationship to Candidate 3:**
**THIS IS THE CRITICAL CONSTRAINT.**

Borate Lewis acids can catalyze DOL ring-opening
polymerization. This paper uses this intentionally.
For Candidate 3 (DOL:B(OCH3)3:LiFSI), TMB is a
weaker Lewis acid than THB (trifluoroisopropyl
borate). However the polymerization risk is real
and must be assessed before any CE measurement.

```
EXPERIMENTAL CONSTRAINT ADDED TO CANDIDATE 3:

Step 0 (before any CE measurement):
  Verify DOL:TMB:LiFSI solution stability.
  Check for viscosity increase or gelation at
  target ratios (20:1, 10:1, 5:1 DOL:TMB vol).
  If gelation occurs at any ratio:
    Either reduce TMB concentration further
    Or switch to Candidate 3b (see below)
  TMB is weaker Lewis acid than THB so may
  not trigger polymerization at low concentration.
  Must be verified, not assumed.
```

### Paper 4 — RSC Journal of Materials Chemistry A 2024
*"Fluorination promotes lithium salt dissolution
in borate esters for lithium batteries"*

```
System: Fluorinated borate esters vs unfluorinated
        (including trimethyl borate) as electrolyte
        components
Finding: Fluorination of the borate ester increases
         Lewis acid strength at boron.
         Stronger Lewis acid → stronger FSI- binding
         → greater salt dissociation → higher conductivity
         → better performance at same concentration
Key:     Fluorinated borate analogs outperform TMB
         for the anion receptor function
```

**Relationship to Candidate 3:**
This paper proves the anion receptor mechanism for
borate esters (consistent with MU Pro response) and
shows that fluorinated analogs are more effective.

**This generates Candidate 3b:**

```
Candidate 3b:
  LiFSI 1.0–1.2M in DOL:fluorinated borate ester
  Stronger anion receptor than TMB
  Reaches same SCE reduction at lower concentration
  Lower gelation risk (stronger Lewis acid at lower vol%)
  Commercially available examples:
    Trifluoroethyl borate: B(OCH2CF3)3
    Or tris(2-fluoroethyl) borate variants
  Same SCE* target: 1.466
  Same mechanism: FSI- trapped, CIP/AGG reduced,
  SSIP fraction rises toward 38%, SCE falls
  from DOL 1.606 toward 1.466
  Novel: Yes — not designed toward SCE* from
  any prior paper
```

### Paper 5 — Science China Chemistry 2025
*"Designing boron-based anion acceptors as
electrolyte additives for lithium batteries"*

```
System: General review of boron-based anion
        acceptors in lithium battery electrolytes
Finding: Confirms electron-deficient boron center
         is the active site for anion binding.
         Stronger Lewis acidity → better performance.
         Fluorinated borates > alkyl borates.
         TMB mentioned as baseline.
No SCE framing. No dual-temperature optimization.
No coordination entropy target.
```

**Relationship to Candidate 3:**
Confirms the mechanism. Confirms the fluorinated
analog hierarchy. Does not pre-empt Candidate 3.

### Paper 6 — Nature Energy 2025
*"Unified affinity paradigm for the rational design
of high-efficiency lithium metal electrolytes"*

```
Descriptor used: Normalized cation/anion-solvent
                 affinity — a single-molecule
                 binding energy descriptor
Finding: >99.8% RT CE achieved by optimizing
         affinity. Molecular-level descriptor.
SCE:     Not computed. Not mentioned. Not derived.
         No Shannon entropy of coordination
         microstate distribution.
         No dual-temperature optimization.
         No SCE* derivation.
         No band identification.
```

**Relationship to Candidate 3 and SCE framework:**
The affinity paradigm is the closest adjacent paper
to the SCE framework in the 2025 literature. It
optimizes at the molecular level. SCE optimizes at
the population-statistical level. They are one level
apart in description. Affinity tells you how strongly
a molecule binds. SCE tells you what the resulting
population distribution looks like statistically and
what its Shannon entropy is. You cannot derive SCE*
from affinity alone because affinity does not capture
the population-level diversity.

Not a novelty threat. A near-miss that strengthens
the SCE framework's distinctiveness.

---

## What the Search Confirmed

### Confirmed Novelty

```
Prior TMB use:
  SEI former (Ding 2023)
  Cathode protector (Yang 2021)
  Polymer initiator (RSC 2024 — unintended risk)

TMB in DOL:LiFSI designed toward SCE* = 1.466
via anion receptor mechanism to reduce coordination
entropy from above the target:

  NOT IN ANY PAPER FOUND.
  NOT DESIGNED TOWARD.
  NOT FRAMED THIS WAY.

Candidate 3 is novel.
```

### New Experimental Constraint

```
DOL polymerization risk from TMB Lewis acid.
Must verify stability at Step 0 before any
CE measurement.
Mitigation: low TMB concentration (<5 vol%),
or switch to Candidate 3b.
```

### New Candidate Generated — 3b

```
Candidate 3b: DOL:fluorinated borate ester:LiFSI
  Stronger anion receptor
  Lower concentration needed
  Lower polymerization risk (less vol%)
  Same SCE* = 1.466 target
  Same mechanism, stronger lever
  Published as better TMB analog (RSC 2024)
  Not used for SCE* targeting — novel
```

---

## Updated Candidate 3 Experimental Path

```
Step 0 (NEW — from internet search):
  Stability verification
  Mix DOL:B(OCH3)3:LiFSI at target ratios
  Check for gelation or viscosity change
  If stable: proceed to Step 1
  If gelation: switch to Candidate 3b

Step 1 — Ratio optimization:
  Raman spectroscopy at 20:1, 10:1, 5:1 ratios
  Monitor FSI- band (727/735/747 cm⁻¹)
  Find ratio where SSIP ≈ 38%

Step 2 — RT CE:
  Li|Cu Aurbach protocol
  Target ≥90%

Step 3 — LT CE at -20°C:
  Li|Cu protocol at -20°C
  Target ≥60%

Step 4 — Framework confirmation:
  CE_RT ≥90% + CE_LT ≥60% = band achieved
  SCE ≈ 1.466 confirmed by Raman
  Candidate 3 validates framework from above
```

---

## The Three Candidates — Final State

| # | System | Route | Base SCE | Direction | Status |
|---|--------|-------|----------|-----------|--------|
| 1 | LiFSI 1.2M in 2-MeTHF:DME 1.6:1 | Cyclic/linear mixing | Mixed | From below toward SCE* | Derived. MD needed. |
| 2 | LiFSI 1.0M in FEME:2-MeTHF 60:40 | Fluorinated/cyclic | Mixed | From below toward SCE* | Derived. MD needed. |
| 3 | LiFSI 1.0–1.2M in DOL:B(OCH3)3 | Anion receptor | 1.606 | From above toward SCE* | Derived. Step 0 needed. |
| 3b | LiFSI 1.0–1.2M in DOL:fluorinated borate | Stronger anion receptor | 1.606 | From above toward SCE* | Generated by search. Novel. |

Four systems. Three chemical routes. One SCE*.

---

## One-Paragraph Summary

The internet search confirmed that trimethyl borate
has been studied as an SEI former and cathode
protector in lithium metal batteries, but never
as a deliberate anion receptor designed to reduce
coordination entropy from above SCE* = 1.466 in
a DOL:LiFSI ether system. Candidate 3 is novel.
The search identified one new experimental
constraint — DOL polymerization risk from borate
Lewis acid catalysis — that must be verified before
CE measurement. The search also identified fluorinated
borate esters as stronger anion receptor analogs
(RSC 2024), generating Candidate 3b: DOL plus a
fluorinated borate ester at lower concentration,
same SCE* target, stronger lever, lower gelation
risk. The Nature Energy 2025 affinity paradigm
paper reappeared and remains a near-miss — it
optimizes at the molecular binding level, not at
the population-statistical entropy level. SCE
cannot be derived from affinity alone. The SCE
framework remains novel. The candidate count is
now four: Candidates 1, 2, 3, and 3b.

---

*Repository: Eric-Robert-Lawson/attractor-oncology*
*Author: Eric Robert Lawson /*
