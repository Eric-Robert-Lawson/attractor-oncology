# MU Search — Query 4: Phosphate Anchor Convergence
**Author:** Eric Robert Lawson / OrganismCore
**Date:** 2026-04-03
**ORCID:** 0009-0002-0414-6544
**Repository:** github.com/Eric-Robert-Lawson/attractor-oncology
**Anchors:** COP(=O)(OC)OC, COP(=O)(OC)OC(C)C, CCOP(=O)(OC)OC
**Setting:** Solvent
**Status:** COMPLETE — PHOSPHATE SPACE SATURATED
**Significance:** Convergence confirmed. Phosphate ester
                 neighborhood bounded. No new chemical
                 classes returned. First non-ether
                 anchored search in this study.

---

## Search Configuration

First search anchored entirely from the phosphate
ester class — not from ethers. Anchors are the
three phosphate esters identified in Search 3:

| SMILES | Identity | Score (Search 3) | ESP Min (eV) |
|--------|----------|------------------|--------------|
| COP(=O)(OC)OC | TMP | 9.40 | -1.7402 |
| COP(=O)(OC)OC(C)C | Isopropyl DMP | 9.27 | -2.0386 |
| CCOP(=O)(OC)OC | Ethyl DMP | 9.31 | -1.8544 |

---

## Resolved Anchors From MU

All three anchors resolved with their Search 3
property suitability scores:

| SMILES | Identity | Property Suitability |
|--------|----------|----------------------|
| COP(=O)(OC)OC | TMP | 9.40/10 |
| COP(=O)(OC)OC(C)C | Isopropyl DMP | 9.27/10 |
| CCOP(=O)(OC)OC | Ethyl DMP | 9.31/10 |

All three resolved. Compare to Search 3 where
DEE failed to resolve. The phosphate ester
class is fully represented in MU's database.

---

## Raw Results — Deduplicated

| SMILES | Score | Status | HOMO (eV) | ESP Min (eV) |
|--------|-------|--------|-----------|--------------|
| COP(=O)(F)OC | 9.36 | Published | -8.4699 | -1.7882 |
| CCOP(=O)(F)OC | 9.32 | Unexplored | -8.3440 | -1.7227 |
| CCOP(=O)(F)OC | 9.29 | Unexplored | -8.1288 | -1.7474 |
| CC#N | 9.27 | Published | -9.1441 | -1.6116 |
| CCOC(=O)OC | 9.19 | Published | -7.9355 | -1.3605 |
| COC(OC)OC | 9.19 | Published | -7.4318 | -1.4881 |
| COC(C)=O | 9.17 | Published | -7.5720 | -1.4350 |
| COC(=O)CF | 9.17 | Published | -7.9406 | -1.5753 |
| COC(=O)OC | 9.16 | Published | -7.9768 | -1.3606 |
| COP(C)(=O)OC | 9.16 | Published | -7.7632 | -1.9417 |

---

## The Key Finding — Phosphate Space Saturated

```
The return is identical to Search 3.
Same molecules. Same scores. Same ordering.

Three phosphate ester anchors returned the
same chemical classes as two ether anchors
pointing toward the phosphate region.

This is the same convergence that occurred
in Search 3 when DEE failed to resolve
and the return matched a two-anchor search.

The phosphate ester molecular neighborhood
in MU's space is bounded by TMP and its
immediate alkyl analogs. Adding isopropyl
and ethyl phosphate analogs as anchors
did not pull the search into new territory.

The phosphate ester anchor search series
is complete after one search.
```

---

## Structural Classification

### Class A — Fluorophosphonates (PRIMARY RETURN)

```
COP(=O)(F)OC     Score 9.36  Published
  HOMO -8.47 eV   ESP Min -1.7882 eV

CCOP(=O)(F)OC    Score 9.32  Unexplored
  HOMO -8.34 eV   ESP Min -1.7227 eV

CCOP(=O)(F)OC    Score 9.29  Unexplored
  HOMO -8.13 eV   ESP Min -1.7474 eV
```

Same class as Search 3 Class B. Fluorophosphonates
are the stable neighbor class of the phosphate
ester anchors in MU's space. The two unexplored
chiral isomers return again at the same scores.

### Class B — Dimethyl Methylphosphonate (PERSISTENT)

```
COP(C)(=O)OC     Score 9.16  Published
  HOMO -7.76 eV   ESP Min -1.9417 eV
```

Appeared in Search 3 and reappears here.
ESP Min -1.9417 eV — between TMP (-1.7402)
and isopropyl DMP (-2.0386). Structural
neighbor within the phosphonate family.
Not a new class. Persistent secondary signal
same as silicon ester in ether space.

### Class C — Acetals (MODERATE — BASE SOLVENTS)

```
COC(OC)OC    Score 9.19  Published  ESP Min -1.4881
COC(C)=O     Score 9.17  Published  ESP Min -1.4350
```

Same acetal class as Search 3. Strong donors.
Not modifiers — would sit above SCE* as base
solvents. Not discarded but not candidates.

### Class D — Carbonates and Esters (DISCARD)

```
CCOC(=O)OC   EMC    9.19  Regime 3 risk
COC(=O)CF    MFA    9.17  Regime 3 risk
COC(=O)OC    DMC    9.16  Regime 3 risk
```

Same discard class across all searches.

### Class E — Nitrile (DISCARD)

```
CC#N   9.27  Li metal incompatible
```

Same discard as all prior searches.

---

## What Search 4 Established

```
1. Phosphate ester space converges after
   one non-ether-anchored search.
   The boundary of the phosphate ester
   neighborhood is TMP and its alkyl analogs.
   No new chemical classes returned.

2. The fluorophosphonate class is the stable
   neighbor of the phosphate ester class —
   same relationship as silicon ester was to
   the boron ester and phosphate ester classes.

3. Dimethyl methylphosphonate (COP(C)(=O)OC)
   is a persistent secondary signal with
   ESP Min -1.9417 eV — stronger than TMP
   but same mechanism, same class.
   Not a new candidate. A stronger-lever
   variant of the TMP family.

4. The two unexplored fluorophosphonate
   chiral isomers returned again at the
   same scores. They are confirmed as the
   primary unexplored signal in phosphate
   ester space. They are Candidate 4b
   variants — weaker P=O donors than TMP,
   same from-below mechanism.

5. No new chemical class emerged from
   phosphate-anchored search. The molecular
   space bounded by ethers and phosphate
   esters has been fully mapped by
   Searches 1–4.
```

---

## Complete Search Session Map

```
Search 1a: DME anchor, ether space
           → Linear ethers. Converged.

Search 1b: DOL anchor, cyclic ether space
           → Cyclic ethers. Converged.

Search 2:  FEC+DME anchors, SCE* zone
           → Boron esters. Candidate 3.

Search 3:  DEE+DOL+DME anchors, solvent
           DEE failed to resolve.
           → Phosphate esters. Candidate 4.

Search 4:  TMP+isopropyl DMP+ethyl DMP anchors
           All three resolved.
           → Same phosphate ester neighborhood.
           Convergence confirmed.
           Phosphate space saturated.
```

The ether-phosphate-boron region of MU's
molecular space is fully mapped. All searches
have converged. No new chemical classes
remain to be found from these anchor sets.

---

## What This Means for Next Steps

```
The search series within the ether-phosphate-
boron triangle is complete.

What remains geometrically:
  Arctic Candidate A (DME:TMP:TMB:LiFSI)
  — derived without search, requires MD

  The question of whether there are chemical
  classes OUTSIDE the ether-phosphate-boron
  triangle that could reach SCE ≈ 2.0–2.5
  — this requires anchors from a completely
  different region of chemical space
  — sulfone class, carbonate class (Regime 3),
  nitrile class (Li metal incompatible) have
  all been returned and discarded
  — no viable fourth class has appeared

The geometry has stopped returning new
chemical classes from the current anchor
region. The search series is complete.
```

---

## Document Chain

| Document | Content |
|----------|---------|
| MU_SEARCH_QUERY_1A.md | Search 1a — linear ethers |
| MU_SEARCH_QUERY_1B.md | Search 1b — cyclic ethers |
| SEARCH_QUERY_2.md | Search 2 — boron esters |
| MU_Search_Query_3_Phosphate_Corrected.md | Search 3 — phosphate esters |
| MU_Search_Query_4_Phosphate_Anchor.md | This document |

---

*Repository: Eric-Robert-Lawson/attractor-oncology*
*Author: Eric Robert Lawson / OrganismCore*
*ORCID: 0009-0002-0414-6544*
*Date: 2026-04-03*
*File: MU_Search_Query_4_Phosphate_Anchor.md*
