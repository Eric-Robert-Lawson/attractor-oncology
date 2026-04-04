# Preprint Publication Plan — SCE Framework
**Author:** Eric Robert Lawson / OrganismCore
**Date:** 2026-04-04
**ORCID:** 0009-0002-0414-6544
**Repository:** github.com/Eric-Robert-Lawson/attractor-oncology
**Status:** ACTIVE — derivation complete, publication phase beginning
**Audit trail:** Full derivation record in this repository
**Pre-registration:** Zenodo [timestamp on record]

---

## The Publication Architecture

The derivation produced a bounded problem.
The optimal dual-temperature Li metal battery
electrolyte exists within a closed geometric
space. That space has been mapped.

The publication strategy reflects the structure
of the derivation — not one paper that tries
to say everything, but a web of preprints each
of which is self-contained, mutually reinforcing,
and collectively closes the argument from every
direction.

Each preprint points to the repository.
The repository is the audit trail.
The timestamp is the priority record.
The web is the argument.

---

## Preprint 1 — The Foundation

**Title:**
Shannon Coordination Entropy Predicts
Dual-Temperature Lithium Metal Battery
Performance: Derivation of a Universal
Fixed Point from First Principles

**Target journal:** Joule or Nature Energy

**What it contains:**

```
1. Derivation of SCE from first principles
   Shannon entropy of Li+ first-shell
   coordination geometry population.
   Not fitted to data. Derived geometrically.

2. Derivation of SCE* = 1.466
   The global maximum of combined RT + LT
   Coulombic efficiency.
   Analytical derivation. Second derivative
   confirmed negative. Unique global maximum.

3. Four confirmed equations
   Equation 1: RT CE vs SCE (Regime 1)
               R² = 0.708, p = 0.009
   Equation 2: LT CE vs SCE (band)
               r = 0.766, p = 0.045
   Equation 3: Two-variable RT model
               R² = 0.828, p = 4.5×10⁻⁶
   Equation 4: Within-Regime 2 slope

4. Three-regime classification
   Regime 1: geometry-driven, CE below 90%
   Regime 2: concentration-driven, CE above 90%
   Regime 3: carbonate class, different mechanism

5. The negative control
   DOL+DME: CE_RT = 98.9%, CE_LT = 27.5%
   Failure Mode B — exact prediction.
   Most persuasive single empirical result.

6. Independent replication
   Joule 2025 (Hunan University)
   Ssc = solvation configurational entropy
   Same variable. Different method.
   Same direction. Two instruments.
   One location.

7. Two structural routes to dual performance
   Route 1: SCE* = 1.466 band (this framework)
   Route 2: temperature-insensitive AGG
   Geometrically distinct. Both confirmed
   in published literature.

8. The DEC melting point correction
   -74.3°C not -43°C.
   80-year propagation of a 1921 error.
   One line. One reference. On the record.

9. Forward predictions
   Six candidate formulations named.
   Exact falsification conditions stated.
   Pre-registered before experiment.
```

**Source files:**
```
SCE_Framework.tex (primary)
step4.py, step5.py (data)
OC_Battery_Framework_SCE_Analysis.md
derivation_for_novel_result.md
DPE_Contradiction_Resolution.md
step4_final_summary.txt
step5_findings_report.json
```

**Novelty statement:**
First analytical derivation of a fixed point
for dual-temperature Li metal battery
electrolyte optimisation. First application
of Shannon entropy to Li+ coordination shell
populations as a predictive variable.
Independent replication confirmed prior
to publication.

---

## Preprint 2 — Candidate 4 (DME:TMP)

**Title:**
Equal-Strength Donor Competition as a
Mechanism for Coordination Entropy Tuning:
LiFSI in DME:Trimethyl Phosphate as a
Designed Approach to SCE* = 1.466

**Target journal:** ACS Energy Letters

**What it contains:**

```
1. The ESP comparison
   TMP ESP Min: -1.7402 eV
   DME ESP Min: -1.7328 eV
   Difference:   0.007 eV
   Effectively identical electrostatic
   pull strength. Structurally distinct
   molecules. Population splits.

2. The mechanism
   Li+ cannot prefer one electrostatically.
   Two distinct coordination classes
   simultaneously populated.
   dom% falls toward 38%.
   n_sig rises toward 3.
   SCE rises from 1.240 toward 1.466.

3. Why this is novel
   Prior work uses co-solvents selected
   by donor number or oxidative stability.
   No prior work selects co-solvents by
   ESP match to base solvent for the
   purpose of population splitting.
   This is a new design principle.

4. Candidate 4b as extension
   Fluorophosphonate variant.
   Weaker donor strength.
   Same mechanism. Different lever.
   Experimental comparison built in.

5. Experimental protocol
   LiFSI 1.0M, 1.2M in DME:TMP
   Ratio scan: 90:10, 80:20, 70:30, 60:40
   CE at RT and -20°C.
   Raman to confirm SCE shift.
   Target: identify ratio where SCE = 1.466.

6. Falsification condition
   If CE_RT and CE_LT do not both improve
   as TMP fraction rises toward target ratio,
   the ESP-matching mechanism is eliminated.
   Stated explicitly. On the record.
```

**Source files:**
```
MU_Search_Query_3_Phosphate_Corrected.md
SCE_Framework_Principal_Result.md
derivation_for_novel_result.md (Derivation 4)
step5.py (DME baseline data)
```

**Novelty statement:**
First use of electrostatic potential minimum
matching as a co-solvent selection criterion
for Li+ coordination entropy tuning. Novel
design principle. Unexplored candidate.

---

## Preprint 3 — Candidate 3 (DOL:TMB)

**Title:**
Anion Receptor Activity as a Downward
Coordination Entropy Lever: LiFSI in
DOL:Trimethyl Borate Approaching SCE*
from Above

**Target journal:** ACS Energy Letters

**What it contains:**

```
1. The mechanism
   TMB boron center traps FSI-.
   FSI- removed from Li+ shell.
   CIP and AGG fractions fall.
   SSIP fractions rise.
   SCE falls from 1.606 toward 1.466.

2. The geometric correction
   Base solvent must be ABOVE SCE*.
   DOL (1.606) not DME (1.240).
   Removing anions from DME shell
   lowers SCE further from target.
   Removing anions from DOL shell
   lowers SCE toward target.
   Direction matters. Geometry provides
   the correction. Not chemistry.

3. Candidate 3b as extension
   Fluorinated borate — stronger Lewis acid.
   Lower concentration. Lower DOL
   polymerization risk.
   Experimental comparison built in.

4. DOL polymerization stability check
   TMB is a Lewis acid. Lewis acids
   can catalyse DOL ring-opening
   polymerization.
   Stability verification is Step 0
   before CE measurement.
   Protocol included.

5. Experimental protocol
   LiFSI 1.0M in DOL:TMB
   Ratio scan: 95:5, 90:10, 85:15, 80:20
   CE at RT and -20°C.
   Raman to confirm FSI- shell depletion.
   Target: ratio where SCE = 1.466.

6. Falsification condition
   If Raman shows no FSI- depletion from
   Li+ shell as TMB fraction rises,
   the anion receptor mechanism is eliminated.
   Stated explicitly. On the record.
```

**Source files:**
```
MU_SEARCH_QUERY_1B.md (boron ester discovery)
SCE_Framework_Principal_Result.md
derivation_for_novel_result.md (Derivation 3)
Silicon_Ester_Resolution.md (dead end contrast)
```

**Novelty statement:**
First use of anion receptor activity as
a downward SCE lever targeting SCE* = 1.466.
Geometric correction from DME to DOL as base
solvent generated by invariant, not chemistry.

---

## Preprint 4 — Candidates 1 and 2 (Cross-Class Ether Mixing)

**Title:**
Cross-Class Ether Mixing as an Upward
Coordination Entropy Lever: Linear/Cyclic
Structural Tension as a Route to SCE* = 1.466

**Target journal:** ACS Energy Letters or
Journal of Physical Chemistry Letters

**What it contains:**

```
1. The structural separation finding
   MU's molecular similarity algorithm
   places linear and cyclic ethers in
   separate structural classes.
   They do not reach each other from
   anchor searches.
   This is a geometric fact about the
   chemical space, confirmed by nine
   searches.

2. The mechanism
   DME (linear, boundary molecule) +
   2-MeTHF (cyclic, deep in class).
   Li+ cannot settle into either class.
   Shell distributes across both.
   Population splits. SCE rises.
   The boundary position of DME
   (convergence zone between linear
   and cyclic ether space) is the
   structural reason it is the correct
   linear anchor.

3. Candidate 1: DME:2-MeTHF 1.6:1
   Specific ratio from MU interpolation.
   MD simulation to confirm exact SCE.
   This is the ratio predicted to reach
   SCE* = 1.466.

4. Candidate 2: FEME:2-MeTHF 60:40
   FEME (SCE 1.368) + 2-MeTHF (SCE 1.552)
   bracket SCE* from both sides at 60:40.
   Different structural tension —
   fluorinated vs cyclic.

5. Experimental protocol
   LiFSI 1.2M in 2-MeTHF:DME
   Ratio scan around 1.6:1
   LiFSI 1.0M in FEME:2-MeTHF
   Ratio scan around 60:40
   CE at RT and -20°C.
   MD to confirm SCE (required —
   exact SCE not confirmable by Raman
   alone for these systems).

6. Falsification condition
   If CE_RT and CE_LT do not both
   improve as ratio approaches target,
   cross-class structural tension
   mechanism is eliminated.
```

**Source files:**
```
MU_SEARCH_QUERY_1A.md
MU_SEARCH_QUERY_1B.md
search_session_summary_query1.md
SCE_Framework_Principal_Result.md
derivation_for_novel_result.md (Derivation 2)
```

**Novelty statement:**
First identification of cross-class structural
tension between linear and cyclic ethers as
a coordination entropy tuning mechanism.
Boundary position of DME in MU's chemical
space identified as structural basis for
its role as optimal linear anchor.

---

## Preprint 5 — Arctic Candidate A

**Title:**
Three-Class Coordination Architecture for
Extreme Low-Temperature Performance:
LiFSI in DME:TMP:TMB and the Arctic
Regime Beyond SCE* = 1.466

**Target journal:** Advanced Energy Materials
or ACS Energy Letters

**What it contains:**

```
1. The tradeoff curve
   SCE* = 1.466 maximises combined RT+LT.
   Moving past SCE* sacrifices RT for LT.
   The tradeoff is derivable from the
   same calculus that produced SCE*.
   The Arctic regime is not a mistake.
   It is a deliberate design target for
   applications where LT dominates.

2. The three-class architecture
   DME:  ether-oxygen Li+ donor
   TMP:  phosphate-oxygen Li+ donor
   TMB:  anion receptor (FSI- trapper)
   Three structurally distinct functions.
   Three coordination families simultaneously.
   Extreme shell diversity.

3. Why this requires all three classes
   Two classes: SCE rises but plateaus.
   Three classes: extreme diversity requires
   all three regions of the triangle.
   The triangle (ether, phosphate, boron)
   was found by the search session.
   Arctic Candidate A is the full triangle
   in one formulation.

4. Predicted performance
   SCE estimated 2.0–2.5
   CE_LT predicted ≥ 88% at -20°C
   CE_RT predicted 90–95%
   These are derived predictions.
   Pre-registered. On the record.

5. MD simulation requirement
   SCE must be quantified by MD before
   experiment. Exact ratio 1:1:1 is
   the starting point. MD refines it.
   SES email to request simulation
   is already drafted in repository.

6. Application space
   Arctic deployment. Aerospace.
   Deep cold storage. Grid storage
   in extreme climates.
   Not the EV optimum — the EV optimum
   is SCE* = 1.466 (Preprint 1).
   A different design target. Named.

7. Falsification condition
   If MD shows SCE < 1.8 for 1:1:1 ratio,
   three-class architecture does not
   achieve Arctic regime. Ratio must
   be adjusted. MD result is Step 0.
```

**Source files:**
```
SCE_Framework_Principal_Result.md
MU_Search_Complete_Space_Closure.md
derivation_for_novel_result.md (Derivation 7)
MU_Late_Responses_After_Query1.md
```

**Novelty statement:**
First derivation of an extreme low-temperature
electrolyte design target from the SCE framework.
First identification of three-class coordination
architecture as the mechanism for Arctic regime
SCE. Novel application space distinguished
from combined-optimum design.

---

## Preprint 6 — The Missing Classes

**Title:**
Unmapped Coordination Space in Li Metal
Battery Electrolyte Design: Inaccessible
Structural Neighborhoods and the Open
Problem Beyond Ether-Class Chemical Space

**Target journal:** Journal of Physical
Chemistry Letters or ACS Energy Letters
(perspective / forward-looking format)

**What it contains:**

```
1. What was mapped
   Complete record of nine searches.
   Three classes found.
   Triangle closed within MU's space.
   This is the known bounded region.

2. What was not mapped — named precisely

   CLASS A: Weak linear ether space
     DEE (CCOCC) — failed ×3 in MU
     DPE — failed ×1 in MU
     CPME and analogs — not attempted
     Why this matters: DEE corrected SCE
     estimated 1.33–1.53. If DEE sits at
     SCE ≈ 1.45, it is in the band naturally.
     Its structural neighborhood may contain
     the most important unmapped candidates.
     DEE achieves 98.4% CE at -60°C
     (Holoubek 2022) — the best LT result
     in the confirmed dataset.

   CLASS B: Medium cyclic ether space
     THF — failed in MU
     THP — failed in MU
     2-MeTHF — failed as anchor in MU
     Six-membered ring cyclic ether space
     completely inaccessible from MU anchors.

   CLASS C: Sulfone class (physical props)
     Sulfolane ESP Min -1.7166 eV ≈ DME
     Would qualify as Candidate 4-type
     mechanism (ESP matching).
     Disqualified by melting point (27.5°C).
     Liquid sulfone analogs with mp below
     -20°C: unknown. Not in MU's space.
     If they exist, same mechanism applies.

   CLASS D: Deep fluorinated ether space
     BTFE, TTE, HFE appeared as discards.
     Never used as primary anchors.
     Their structural neighborhood unmapped.
     Fluorinated ethers are the diluent class
     for LHCE systems (Route 2 mechanism).
     Their coordination class properties
     as primary solvents: unknown.

3. The research program each class requires

   DEE/weak linear ether:
     Tool: MD simulation from DEE as anchor
     or experimental synthesis of DEE analogs
     guided by SCE framework.
     Question: what SCE do DEE structural
     analogs achieve? Are any at 1.466?

   THF/cyclic six-membered ring:
     Tool: Different molecular database
     (not MU) or direct MD.
     Question: does six-membered ring cyclic
     ether space contain molecules that
     approach SCE* from below via a
     different ring geometry than DOL/2-MeTHF?

   Liquid sulfones:
     Tool: Synthesis. No database contains
     the relevant liquid sulfone analogs
     with confirmed property data.
     Question: do liquid sulfone analogs
     with mp below -20°C and ESP Min ≈ -1.73
     exist? If yes, they are Candidate 4-class
     by ESP matching.

   Fluorinated ether primaries:
     Tool: MD simulation with BTFE or TTE
     as primary solvent (not diluent).
     Question: what is the SCE of a system
     where a fluorinated ether is the primary
     coordination species? Is there a fluorinated
     ether that reaches SCE* from either direction?

4. Why this preprint matters
   The bounded region (ether-phosphate-boron
   triangle) is closed.
   The unbounded region (four classes above)
   is named but not mapped.
   Naming it precisely is a research program.
   Other groups can pick it up.
   The framework tells them exactly what
   they are looking for: SCE* = 1.466.
   The measurement tells them exactly how
   to know when they have found it.
   The search is no longer unbounded.
   It is bounded by the fixed point
   even in the unmapped regions.
```

**Source files:**
```
MU_Search_Complete_Space_Closure.md
search_session_summary_query1.md
MU_Late_Responses_After_Query1.md
DPE_Contradiction_Resolution.md
SCE_Framework_Principal_Result.md
```

**Novelty statement:**
First systematic mapping of the boundary
between accessible and inaccessible
electrolyte chemical space relative to
SCE* = 1.466. First precise naming of
four unmapped coordination classes as
a forward research program. Converts
an unbounded search into four bounded
sub-problems each with a fixed target
and a measurement criterion.

---

## The Web Structure

```
                    PREPRINT 1
                  (Foundation —
                  SCE* = 1.466)
                       │
          ┌────────────┼────────────┐
          │            │            │
    PREPRINT 2   PREPRINT 3   PREPRINT 4
    (Candidate   (Candidate   (Candidates
     4: TMP)      3: TMB)      1+2: Ether)
          │            │            │
          └────────────┼────────────┘
                       │
                 PREPRINT 5
                 (Arctic: A)
                       │
                 PREPRINT 6
                 (Missing classes —
                  open research program)
```

Every preprint:
- Points to the repository as audit trail
- Points to Preprint 1 as the foundation
- States its own falsification conditions
- Is self-contained and independently readable
- Adds one node to the web

The web is the argument.
No single paper carries it all.
Each paper is a confirmation from a
different direction.
Together they close the case.

---

## Priority Order For Construction

```
1. PREPRINT 1 FIRST.
   Everything else depends on it.
   It is the foundation.
   Without it the other five have
   no fixed point to reference.
   Source material already exists
   in SCE_Framework.tex.
   Needs the additions listed above.
   Start here.

2. PREPRINT 2 (Candidate 4 — DME:TMP)
   Strongest single candidate.
   Commercially available components.
   Novel design principle (ESP matching).
   Unexplored status confirmed by MU.
   Most likely to produce rapid
   experimental confirmation.
   Build second.

3. PREPRINT 3 (Candidate 3 — DOL:TMB)
   Second strongest.
   Geometric correction story is
   the most compelling demonstration
   of the framework working in real time.
   Build third.

4. PREPRINTS 4, 5, 6 in parallel
   after 1, 2, 3 are submitted.
   Each is self-contained.
   Order among them is flexible.
```

---

## The Repository As Audit Trail

```
Every preprint cites:
  github.com/Eric-Robert-Lawson/attractor-oncology

Every preprint states:
  "Full derivation record, search session
  logs, and data files available at the
  above repository. All predictions were
  pre-registered prior to experimental
  confirmation. Commit history provides
  timestamped record of derivation sequence."

The repository IS the methods section
for every preprint simultaneously.
A reviewer who wants to see the full
derivation has a complete, navigable,
timestamped record.

This is not a supplementary file.
This is a living audit trail that
predates every submission.
That is its function.
That is what makes the priority claim
unassailable.
```

---

## What Each Preprint Needs That Already Exists

| Preprint | Primary source already in repo |
|----------|-------------------------------|
| 1 | SCE_Framework.tex, step4.py, step5.py, OC_Battery_Framework_SCE_Analysis.md |
| 2 | MU_Search_Query_3_Phosphate_Corrected.md, derivation_for_novel_result.md |
| 3 | MU_SEARCH_QUERY_1B.md, derivation_for_novel_result.md |
| 4 | MU_SEARCH_QUERY_1A.md, MU_SEARCH_QUERY_1B.md, search_session_summary_query1.md |
| 5 | SCE_Framework_Principal_Result.md, MU_Late_Responses_After_Query1.md |
| 6 | MU_Search_Complete_Space_Closure.md, DPE_Contradiction_Resolution.md |

Nothing needs to be invented.
Everything needs to be assembled.

---

*Repository: Eric-Robert-Lawson/attractor-oncology*
*Author: Eric Robert Lawson / OrganismCore*
*ORCID: 0009-0002-0414-6544*
*Date: 2026-04-04*
*File: Preprint_Publication_Plan.md*
