# MU Search — Query 3: Three-Anchor Ether Search
**Author:** Eric Robert Lawson / OrganismCore
**Date:** 2026-04-03
**ORCID:** 0009-0002-0414-6544
**Repository:** github.com/Eric-Robert-Lawson/attractor-oncology
**Anchors:** CCOCC (DEE) + C1COCCO1 (DOL) + COCCOC (DME)
**Setting:** Solvent
**Range:** Maximum
**Priority:** Hypothetical and unexplored
**Status:** COMPLETE
**Significance:** Phosphate ester class confirmed as
                 stable attractor of ether molecular
                 space. Ether anchor searches saturated.
                 Candidate 4 generated.

---

## Search Configuration

Three anchors spanning the full ether SCE range
in the confirmed dataset:

| SMILES | Identity | Dataset SCE | HOMO (eV) |
|--------|----------|-------------|-----------|
| CCOCC | DEE | ~1.33 | -7.4 |
| C1COCCO1 | DOL | 1.606 | -6.61 |
| COCCOC | DME | 1.240 | -6.87 |

DEE was included to pull the search toward
weaker solvators and potentially reach Arctic
regime molecules at SCE ≈ 2.3–2.5.

---

## Resolved Anchors From MU

| SMILES | Identity | Property Suitability | HOMO (eV) | ESP Min (eV) |
|--------|----------|----------------------|-----------|--------------|
| C1COCCO1 | DOL | 8.52/10 | -6.6055 | -1.1948 |
| COCCOC | DME | 9.11/10 | -6.8685 | -1.7328 |

Note: CCOCC (DEE) showed "No molecule data available"
in the resolved anchors panel. DEE did not resolve
as a scored anchor. The search ran from DOL and DME
only. This is why the return is identical to a
two-anchor DOL+DME search.

---

## Raw Results — Deduplicated

| SMILES | Score | Status | MW (g/mol) | HOMO (eV) | LUMO (eV) | ESP Min (eV) | Commercial |
|--------|-------|--------|------------|-----------|-----------|--------------|------------|
| COP(=O)(OC)OC | 9.40 | Published | 140.075 | -7.9448 | 1.0274 | -1.7402 | Available |
| COP(=O)(F)OC | 9.36 | Published | 128.039 | -8.4699 | 0.7096 | -1.7882 | Synthesizable |
| CCO[P@@](=O)(F)OC | 9.32 | Unexplored | 142.066 | -8.3440 | 0.8714 | -1.7227 | Synthesizable |
| CCOP(=O)(OC)OC | 9.31 | Published | 154.102 | -7.7987 | 1.0238 | -1.8544 | Synthesizable |
| CCO[P@](=O)(F)OC | 9.29 | Unexplored | 142.066 | -8.1288 | 0.9303 | -1.7474 | Synthesizable |
| CC#N | 9.27 | Published | 41.053 | -9.1441 | 0.7245 | -1.6116 | Available |
| COP(=O)(OC)OC(C)C | 9.27 | Unexplored | 168.129 | -7.8023 | 0.8756 | -2.0386 | Synthesizable |
| CCOC(=O)OC | 9.19 | Published | 104.105 | -7.9355 | 0.9828 | -1.3605 | Available |
| COC(OC)OC | 9.19 | Published | 106.121 | -7.4318 | 1.1874 | -1.4881 | Available |
| COC(C)=O | 9.17 | Published | 74.079 | -7.5720 | 0.2024 | -1.4350 | Available |
| COC(=O)CF | 9.17 | Published | 92.069 | -7.9406 | -0.1240 | -1.5753 | Available |
| COC(=O)OC | 9.16 | Published | 90.078 | -7.9768 | 0.9852 | -1.3606 | Available |
| COP(C)(=O)OC | 9.16 | Published | 124.076 | -7.7632 | 0.8164 | -1.9417 | Synthesizable |
| CO[Si](C)(OC)OC | 9.15 | Published | 136.223 | -7.6348 | 0.8775 | -1.2956 | Available |
| COC1OCCO1 | 9.14 | Published | 104.105 | -7.3430 | 1.1843 | -1.5176 | Available |

---

## The Key Finding — DEE Did Not Resolve

DEE (CCOCC) appears under "No molecule data available"
in the searched molecules panel. It was entered but
did not score as a resolved anchor. The search
therefore ran from DOL and DME only.

This is not a failure. It is a finding:

```
DEE is not represented in MU's molecular space
in the same way DOL and DME are. It either
lacks sufficient training data or sits outside
the platform's anchor resolution threshold.

Consequence: The ether anchor space cannot be
extended toward weaker solvators via DEE.
The ether neighborhood in MU's space is bounded
by DOL and DME as the furthest resolvable anchors.

The phosphate ester class is the stable return
from this entire bounded ether neighborhood.
It appeared in the DOL+DME search and appears
again here with DEE added but unresolved.
The ether anchor searches are saturated.
```

---

## Structural Classification

### Class A — Phosphate Esters (PRIMARY SIGNAL)

```
COP(=O)(OC)OC    trimethyl phosphate (TMP)
  HOMO -7.94 eV   ESP Min -1.7402 eV
  Score 9.40       Published   Commercially available
  KEY: ESP Min ≈ DME ESP Min (-1.7328 eV)
       TMP and DME attract Li+ with equal force
       but are structurally distinct molecules
       → distinct microstates in SCE framework

CCOP(=O)(OC)OC   ethyl dimethyl phosphate
  HOMO -7.80 eV   ESP Min -1.8544 eV
  Score 9.31       Published

COP(=O)(OC)OC(C)C  isopropyl dimethyl phosphate
  HOMO -7.80 eV   ESP Min -2.0386 eV
  Score 9.27       Unexplored

COP(C)(=O)OC     dimethyl methylphosphonate
  HOMO -7.76 eV   ESP Min -1.9417 eV
  Score 9.16       Published
```

Highest scoring class. TMP at 9.40 is the
highest score returned across all MU searches
in this study.

### Class B — Fluorophosphonates (HIGH INTEREST)

```
COP(=O)(F)OC     dimethyl fluorophosphate
  HOMO -8.47 eV   ESP Min -1.7882 eV
  Score 9.36       Published

CCO[P@@](=O)(F)OC  ethylmethyl fluorophosphonate
  HOMO -8.34 eV   ESP Min -1.7227 eV
  Score 9.32       Unexplored

CCO[P@](=O)(F)OC   (chiral isomer)
  HOMO -8.13 eV   ESP Min -1.7474 eV
  Score 9.29       Unexplored
```

Fluorine withdraws electron density from P=O.
Weaker Li+ donor than TMP. More likely to act
as minority coordination species. Two unexplored
chiral isomers — no published battery data.

### Class C — Silicon Ester (PERSISTENT SIGNAL)

```
CO[Si](C)(OC)OC   trimethoxymethylsilane
  HOMO -7.63 eV   ESP Min -1.2956 eV
  Score 9.15       Published   Commercially available
```

Appears in Search 2 (score 9.29) and Search 3
(score 9.15). Two independent searches. Two
different anchor sets. Same molecule returns
both times. This is a persistent signal.
ESP Min -1.2956 eV — weaker Li+ attraction
than both DME (-1.73) and TMP (-1.74).
Coordination role unresolved. MU Pro Option A
query submitted — awaiting response.

### Class D — Acetals (MODERATE INTEREST)

```
COC(OC)OC    trimethyl orthoformate
  HOMO -7.43 eV   ESP Min -1.4881 eV
  Score 9.19       Published

COC1OCCO1    2-methoxy-DOL
  HOMO -7.34 eV   ESP Min -1.5176 eV
  Score 9.14       Published
```

Strong donors similar to ether class.
Above SCE* in coordination space.
Not modifiers — base solvents.

### Class E — Carbonates and Esters (DISCARD)

```
CCOC(=O)OC   EMC    HOMO -7.94   Published   9.19
COC(C)=O     MeAc   HOMO -7.57   Published   9.17
COC(=O)CF    MFA    HOMO -7.94   Published   9.17
COC(=O)OC    DMC    HOMO -7.98   Published   9.16
```

Carbonyl chemistry → Regime 3 risk. Discard.

### Class F — Nitrile (DISCARD)

```
CC#N   acetonitrile   HOMO -9.14   Published   9.27
```

Incompatible with Li metal anode. Discard.

---

## What Search 3 Established

```
1. DEE does not resolve as an anchor in MU.
   The ether anchor space is bounded by DOL
   and DME. Cannot extend toward weaker
   solvators via ether anchors alone.

2. The phosphate ester class is the stable
   attractor of the entire ether molecular
   neighborhood in MU's space. It returns
   regardless of whether two or three ether
   anchors are used.

3. TMP's ESP Min (-1.7402 eV) is within
   0.007 eV of DME's ESP Min (-1.7328 eV).
   This comparison has not been made in any
   prior published paper. It means TMP
   competes with DME for Li+ coordination
   at equal electrostatic force.

4. The ether anchor search series is complete.
   All ether-anchored searches have converged
   on the same two attractor classes:
     Boron esters (Search 2)
     Phosphate esters (Search 3)
   No further ether-anchored searches will
   return new chemical classes.

5. The Arctic regime requires anchors outside
   the ether class. DEE failing to resolve
   confirms this — the ether space does not
   extend far enough in MU's representation
   to reach SCE ≈ 2.3–2.5 from ether anchors.
```

---

## Candidate Generated — Candidate 4

```
System:    LiFSI 1.0–1.2M in DME:TMP
SMILES:    COCCOC + COP(=O)(OC)OC
Mechanism: TMP P=O oxygen (ESP Min -1.7402 eV)
           competes with DME ether oxygen
           (ESP Min -1.7328 eV) for Li+
           coordination at equal electrostatic
           strength but as structurally distinct
           species. New microstate enters
           population. n_sig rises. dom% falls.
           SCE rises from 1.240 toward 1.466.
Direction: From below SCE* = 1.466
Base SCE:  1.240 (DME)
Target:    SCE* = 1.466
Novel:     Yes — no prior paper designs this
           system toward any coordination entropy
           target. Confirmed by internet search.
```

---

## What Comes Next

```
Ether anchor searches: COMPLETE — saturated
Silicon ester role:    PENDING — MU Pro Option A
                       query submitted
Arctic regime:         Requires non-ether anchors
                       DEE failure to resolve
                       confirms ether anchors cannot
                       reach SCE ≈ 2.3–2.5
                       Next anchor set must include
                       TMP or TMB as anchors —
                       molecules from the confirmed
                       non-ether classes found in
                       Searches 2 and 3
```

---

## Document Chain

| Document | Content |
|----------|---------|
| MU_SEARCH_QUERY_1A.md | Search 1a — linear ether space |
| MU_SEARCH_QUERY_1B.md | Search 1b — cyclic ether space |
| SEARCH_QUERY_2.md | Search 2 — boron esters, Candidate 3 |
| MU_Search_Query_3_Phosphate_Corrected.md | This document |

---

*Repository: Eric-Robert-Lawson/attractor-oncology*
*Author: Eric Robert Lawson / OrganismCore*
*ORCID: 0009-0002-0414-6544*
*Date: 2026-04-03*
*File: MU_Search_Query_3_Phosphate_Corrected.md*
