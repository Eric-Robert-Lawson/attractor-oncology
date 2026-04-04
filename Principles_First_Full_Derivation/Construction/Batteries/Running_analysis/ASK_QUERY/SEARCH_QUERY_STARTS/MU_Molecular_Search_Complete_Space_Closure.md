# MU Molecular Search — Complete Space Closure Record
**Author:** Eric Robert Lawson / OrganismCore
**Date:** 2026-04-04
**ORCID:** 0009-0002-0414-6544
**Repository:** github.com/Eric-Robert-Lawson/attractor-oncology
**Status:** COMPLETE — MU molecular space closed from all directions
**Searches completed:** 9 (1a, 1b, 2, 3, 4, 5, 6, 7, 8, 9)

---

## What This Document Is

This document records the definitive closure of the MU
Molecular Universe search session for the SCE framework.
It records what was searched, what was found, what failed
to resolve, and what the complete geometric picture means.

It is not a search record. Each search has its own document.
This is the closure statement — the record that the
search space has been exhausted from all accessible
anchor directions.

---

## The Complete Search Record

| Search | Anchors | Setting | What Resolved | Primary Return | Status |
|--------|---------|---------|---------------|----------------|--------|
| 1a | COCCOC (DME) | Solvent | DME only (3 of 4 failed) | Linear glymes, fluorinated linear ethers | Complete |
| 1b | C1COCCO1 (DOL) | Co-solvent | DOL only (3 of 4 failed) | Cyclic DOL analogs, fluorinated cyclic ethers | Complete |
| 2 | COCCOC + C1COCCO1 (FEC+DME) | Solvent | Both | Boron esters → Candidate 3 | Complete |
| 3 | CCOCC + C1COCCO1 + COCCOC | Solvent | DOL + DME (DEE failed) | Phosphate esters → Candidate 4 | Complete |
| 4 | COP(=O)(OC)OC × 3 phosphate anchors | Solvent | All three | Phosphate neighborhood saturated | Complete |
| 5 | COC1OCCO1 + COP(=O)(OC)OC + COB1OB(OC)OB(OC)O1 + CCOCC | Solvent | 3 of 4 (DEE failed ×2) | Phosphate esters — identical to Search 3 | Complete |
| 6 | COC1OCCO1 + COP(=O)(OC)OC + COB1OB(OC)OB(OC)O1 + CCOCC | Solvent | 3 of 4 (DEE failed ×3) | Phosphate esters — identical to Search 5 | Complete |
| 7 | CCOCC alone | Solvent | Failed | No return | Complete — DEE inaccessible |
| 8 | Mixed: all prior anchors combined | Solvent | 3 of 4 | Phosphate esters — identical to all prior | Complete |
| 9 | CS(C)=O + CS(C)(=O)=O + O=S1(=O)CCCC1 + CC(=O)N(C)C | Solvent | All four | Phosphate esters — identical return | Complete |

---

## The Anchors That Failed To Resolve

```
CCOCC         DEE           — failed ×3 (Searches 3, 6, 7)
CCOC          DPE analog    — failed ×1
C1CCOC1       THF           — failed in Search 1a
C1CCOCC1      THP           — failed in Search 1b
CC1CCCO1      2-MeTHF       — failed in Search 1a
```

These molecules sit in regions of chemical space that
MU cannot navigate from. Their structural neighborhoods
are inaccessible — not empty, but outside MU's indexed
representation. What is in those neighborhoods cannot
be determined from MU searches.

---

## The Universal Attractor — TMP

```
COP(=O)(OC)OC   Trimethyl phosphate (TMP)
  Score: 9.40/10 — highest score in every search
  Returned in: Searches 2, 3, 4, 5, 6, 8, 9
  Anchor classes it was returned from:
    Ether anchors (Search 2, 3)
    Phosphate anchors (Search 4)
    Mixed anchors (Search 5, 6, 8)
    Sulfone + amide anchors (Search 9)

TMP is the central attractor of MU's
battery-relevant molecular space.
It is the highest-scoring molecule
returned from every anchor direction
that has been tried.
This is a geometric fact about MU's
representation, not a chemistry claim.
```

---

## The Three Regions Of MU's Battery-Relevant Space

### Region 1 — Ether Class
```
Primary molecules:
  COCCOC        DME         SCE = 1.240  (dataset)
  C1COCCO1      DOL         SCE = 1.606  (dataset)
  CCOCCO        Diglyme G2  SCE unknown
  CCCOCCC       Triglyme G3 SCE unknown
  CC1CCCO1      2-MeTHF     SCE = 1.552  (dataset)

Structural finding:
  Linear and cyclic ether space are SEPARATE
  primary regions in MU's similarity metric.
  They meet only at the convergence zone:
  DME and diglyme appear in both spaces.
  No linear ether anchor reaches cyclic ether
  space and vice versa (confirmed Search 1a, 1b).

Candidates generated:
  Candidate 1: LiFSI 1.2M in 2-MeTHF:DME 1.6:1
  Candidate 2: LiFSI 1.0M in FEME:2-MeTHF 60:40
```

### Region 2 — Phosphate Ester Class
```
Primary molecules:
  COP(=O)(OC)OC     TMP              ESP -1.7402 eV
  COP(=O)(F)OC      DMFP             ESP -1.7882 eV
  CCO[P@@](=O)(F)OC Unexplored       ESP -1.7227 eV
  CCOP(=O)(OC)OC    Ethyl DMP        ESP -1.8544 eV
  CCO[P@](=O)(F)OC  Unexplored       ESP -1.7474 eV
  COP(=O)(OC)OC(C)C Isopropyl DMP   ESP -2.0386 eV
  COP(C)(=O)OC      DMMP             ESP -1.9417 eV

Key finding:
  TMP ESP Min (-1.7402 eV) ≈ DME ESP Min (-1.7328 eV)
  Difference: 0.007 eV — essentially identical.
  TMP competes with DME for Li+ coordination
  at equal electrostatic force as a distinct
  structural species. Population splits.
  SCE rises from DME baseline toward SCE*.

Candidates generated:
  Candidate 4:  LiFSI 1.0–1.2M in DME:TMP
  Candidate 4b: LiFSI 1.0–1.2M in DME:fluorophosphonate
```

### Region 3 — Boron Ester Class
```
Primary molecules:
  COB(OC)OC         Trimethyl borate (TMB)
  COB1OB(OC)OB(OC)O1 Cyclic boroxine

Key finding:
  Boron center is electron-deficient Lewis acid.
  Coordinates FSI- (anion receptor), not Li+.
  Removes FSI- from Li+ shell.
  CIP/AGG fall. SSIP rises.
  SCE falls from above toward SCE*.
  Correct base solvent: DOL (above SCE*),
  not DME (below SCE*) — correction generated
  by invariant geometry, not chemistry.

Candidates generated:
  Candidate 3:  LiFSI 1.0–1.2M in DOL:TMB
  Candidate 3b: LiFSI 1.0–1.2M in DOL:fluorinated borate
```

---

## Region 4 — Discards (All Searches)

```
Class              SMILES examples     Reason discarded
Carbonates         CCOC(=O)OC          Regime 3 risk
                   COC(=O)OC
Nitriles           CC#N                Li metal incompatible
Silicon esters     CO[Si](C)(OC)OC     ESP Min too weak
                                       (-1.2956 eV vs DME -1.7328)
                                       Sterically hindered
Sulfones           CS(C)=O             Melting point too high
                   CS(C)(=O)=O         (DMSO 18.5°C, DMSO2 109°C)
                   O=S1(=O)CCCC1       (Sulfolane 27.5°C)
Amides             CC(=O)N(C)C         Oxidative stability unclear
                                       for Li metal high voltage
Heavy fluorinators FC(F)(F)C(F)(F)OC  Non-coordinating diluents
                                       not Li+ donors
```

---

## The Sulfolane ESP Finding — Noted But Discarded

```
O=S1(=O)CCCC1  Sulfolane
  ESP Min: -1.7166 eV
  DME ESP Min: -1.7328 eV
  Difference: 0.016 eV

Sulfolane's S=O oxygen attracts Li+ at
essentially the same electrostatic strength
as DME. The same mechanism that makes TMP
a Candidate 4 precursor (ESP match to DME)
would apply to sulfolane in principle.

Why it is discarded:
  Melting point 27.5°C — solid at room temperature.
  Cannot function as bulk Li metal battery solvent
  at RT without co-solvent that defeats the purpose.
  Not a candidate. Mechanism noted for completeness.

If a liquid sulfone analog with similar ESP Min
and mp below -20°C exists, it would be a
candidate by the same logic as TMP.
MU did not return such a molecule.
It may exist outside MU's indexed space.
```

---

## What Cannot Be Reached From MU

```
The following structural neighborhoods are
inaccessible from MU's anchor resolution:

1. Weak linear ether space (DEE, DPE, CPME)
   DEE failed ×3. The molecules in this
   neighborhood cannot be found via MU search.
   What is there: unknown.
   Relevance: DEE may be near SCE* (corrected
   SCE estimate 1.33–1.53). Its structural
   analogs are inaccessible.

2. Medium cyclic ether space (THF, THP, 2-MeTHF)
   All failed to resolve as anchors.
   Their structural neighborhoods cannot be
   navigated. The cyclic ether space that
   was mapped (DOL and analogs) is only the
   five-membered ring region. Six-membered
   ring cyclic ether space (THP) is closed.

3. Everything outside MU's database
   MU indexes molecules with sufficient
   property data for its scoring algorithm.
   Novel or obscure molecules not in its
   training set are invisible.
   This includes: custom synthesized molecules,
   recently reported electrolyte components,
   and anything in the DEE/THF/2-MeTHF
   structural neighborhood that MU cannot
   anchor from.
```

---

## The Arctic Regime — Not Closed By Search

```
Arctic Candidate A: LiFSI 1.0M in DME:TMP:TMB 1:1:1
Target: SCE ≈ 2.0–2.5

This candidate was NOT found by search.
It was DERIVED by combining the outputs of
Search 2 (boron ester class) and Search 3
(phosphate ester class) with the confirmed
DME baseline.

The Arctic regime (SCE > 2.0) was confirmed
absent from the published literature by
MU Ask Pro Query 4. No published system
has been designed toward SCE ≈ 2.0–2.5.

The Arctic candidate does not require
a new molecular class. It requires three
molecules from already-found classes
combined in one system:
  DME (ether class, Region 1)
  TMP (phosphate ester class, Region 2)
  TMB (boron ester class, Region 3)

All three regions were needed to derive it.
The search session produced all three regions.
The Arctic candidate is the combination
derivation from the complete triangle.

Status: Derived. MD simulation required
to quantify SCE. Not yet experimentally tested.
```

---

## Geometric Closure — The Precise Statement

```
MU's battery-relevant molecular space has been
probed exhaustively from:

  Ether anchors (linear):    DME
  Ether anchors (cyclic):    DOL
  Phosphate ester anchors:   TMP × 3 variants
  Boron ester anchors:       boroxine
  Mixed class anchors:       all combinations
  External class anchors:    sulfone + amide

Every resolving anchor from every direction
returns the same three regions:
  Ether → phosphate ester
  Ether → boron ester
  Phosphate → phosphate neighborhood
  Boron → phosphate neighborhood
  Mixed → phosphate ester
  Sulfone/amide → phosphate ester

The triangle has no external connections to
new battery-relevant coordination classes
within MU's indexed space.

THIS IS REAL GEOMETRIC CLOSURE OF MU'S SPACE.

It does not mean all chemical space is closed.
It means: within the molecular database MU
can search, anchored from any resolving class,
no new coordination mechanism exists beyond
what has been found.

Chemical space beyond MU's indexed database:
  Open. Requires different tool.
  Cannot be closed by MU searches.
  Noted as a limitation, not a failure.
```

---

## Complete Candidate Summary

| # | System | Direction | Mechanism | Base SCE | Target | MD needed |
|---|--------|-----------|-----------|----------|--------|-----------|
| 1 | LiFSI 1.2M in 2-MeTHF:DME 1.6:1 | From below | Cross-class ether mixing | 1.240 | SCE* 1.466 | Yes — exact SCE |
| 2 | LiFSI 1.0M in FEME:2-MeTHF 60:40 | From below | Fluorinated/cyclic mixing | 1.368 | SCE* 1.466 | Yes — exact SCE |
| 3 | LiFSI 1.0–1.2M in DOL:TMB | From above | Anion receptor (boron) | 1.606 | SCE* 1.466 | No — ratio by Raman |
| 3b | LiFSI 1.0–1.2M in DOL:fluorinated borate | From above | Anion receptor (stronger) | 1.606 | SCE* 1.466 | No — ratio by Raman |
| 4 | LiFSI 1.0–1.2M in DME:TMP | From below | Equal-strength Li+ donor | 1.240 | SCE* 1.466 | No — ratio by Raman |
| 4b | LiFSI 1.0–1.2M in DME:fluorophosphonate | From below | Weaker Li+ donor | 1.240 | SCE* 1.466 | No — ratio by Raman |
| A | LiFSI 1.0M in DME:TMP:TMB 1:1:1 | Past SCE* | Three-class Arctic | 1.240 | SCE 2.0–2.5 | Yes — SCE unknown |

---

## Document Index — Complete Search Record

| # | Document | Content |
|---|----------|---------|
| 1 | MU_SEARCH_QUERY_1A.md | Search 1a — linear ether space |
| 2 | MU_SEARCH_QUERY_1B.md | Search 1b — cyclic ether space |
| 3 | search_session_summary_query1.md | Session summary after 1a + 1b |
| 4 | SEARCH_QUERY_2.md | Search 2 — boron esters, Candidate 3 |
| 5 | MU_Search_Query_3_Phosphate_Corrected.md | Search 3 — phosphate esters, Candidate 4 |
| 6 | MU_Search_Query_4_Phosphate_Anchor.md | Search 4 — phosphate space saturated |
| 7 | Silicon_Ester_Resolution.md | Silicon ester dead end confirmed |
| 8 | MU_Search_Complete_Space_Closure.md | This document — searches 5–9, closure |

---

*Repository: Eric-Robert-Lawson/attractor-oncology*
*Author: Eric Robert Lawson / OrganismCore*
*ORCID: 0009-0002-0414-6544*
*Date: 2026-04-04*
*File: MU_Search_Complete_Space_Closure.md*
