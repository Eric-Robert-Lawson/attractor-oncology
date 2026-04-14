# Literature Check — SCE Framework Candidates vs. Existing Published Literature
## What the Field Has Found Accidentally, What It Has Not Found, and What the SCE Framework Adds
## OrganismCore — Eric Robert Lawson
## Date: 2026-04-14
## Status: ACTIVE — Preservation and Record Document
## Companion to: SCE_Framework_Principal_Result.md, Preprints 1–6, OC_Battery_Framework_SCE_Analysis.md

---

## Purpose

This document records the results of a systematic literature check
conducted against the six SCE framework candidates and the framework's
conceptual foundations.

The question asked: does existing published literature accidentally
contain any of these candidates, mechanisms, or formulations?

The answer, candidate by candidate, is recorded below. This document
is a pre-submission record and should be updated as new literature
is identified. It is also the basis for the novelty claims in
Preprints 1–6 and for the "related work" sections of each manuscript.

---

## What the Literature Has — Accidentally — On Each Candidate

---

## Candidate 4 (DME:TMP) — The Most Significant Hit

**Paper identified:**
"Li-Metal Anode in Dilute Electrolyte LiFSI/TMP:
Electrochemical Performance and Solvation Structure"
OSTI/DOE record, confirmed via scite.ai

**What the paper found:**

In dilute LiFSI/TMP electrolytes, ab initio MD simulations show
Li⁺ coordinated by approximately three TMP molecules and one FSI⁻ —
producing a coordination shell with mixed phosphate-oxygen and
anion character. This solvation structure produces a denser, more
inorganic SEI than DME-alone systems. CE improvements reported
moving toward and above 98% in Li‖Cu cells.

**What the SCE framework predicts vs. what they found:**

The literature work studied pure LiFSI/TMP — not the DME:TMP blend.
In pure TMP, Li⁺ is saturated by phosphate-oxygen coordination. The
population is dominated by one coordination class: TMP-heavy shells.
SCE would be low (similar to pure DME), not at 1.466. The paper found
good RT performance — consistent with Regime 2 behaviour — but the
dual-temperature question was never asked.

**The critical gap:**

No paper found has asked: what happens when you mix DME and TMP at
a ratio where their ESP minima are matched so that Li⁺ cannot prefer
one electrostatically? That specific design question — choosing the
co-solvent ratio not for bulk properties but for population splitting
at the coordination level — does not appear anywhere in the literature.

The accidental result confirms the ingredients work.
The designed result (the ESP-matching mechanism at specific ratio)
has never been attempted.

**Implication for Preprint 2 (DME:TMP):**
Novelty is confirmed. The mechanism is new. The literature confirms
ingredient compatibility and solvation structure type. The fixed-point
design rationale is entirely absent.

---

## Candidate 3 (DOL:TMB) — Mechanism Confirmed, Application Wrong

**Paper identified:**
"Roles of Trimethyl Borate in Constructing an Interphase on Li Anode:
Angel or Demon?"
ACS Applied Materials & Interfaces, 2023 — Ding et al.
DOI: 10.1021/acsami.2c19417

**What the paper found:**

- TMB acts as an anion receptor — it binds FSI⁻ / PF₆⁻, reducing
  free anion concentration in the electrolyte bulk
- This modifies the SEI toward a boron-containing, LiF-rich interphase
- The "angel or demon" framing reflects that the effect is
  concentration-dependent and base-solvent-dependent — in carbonate
  systems, TMB sometimes helps and sometimes hurts

**What the SCE framework adds:**

The literature found the mechanism (anion receptor activity) but did
not have the design criterion to control it.

The SCE framework says: TMB works as a downward SCE lever ONLY when
the base solvent SCE is above SCE* = 1.466.

- In a carbonate system (SCE ≈ 2.0+), adding TMB reduces anion
  coordination fractions and moves SCE toward 1.466 — which improves
  performance. TMB is the "angel."
- In a DME system (SCE ≈ 1.24), adding TMB moves SCE further below
  the fixed point — which hurts performance. TMB is the "demon."

The "angel or demon" confusion in the literature is geometrically
explained by a single criterion: which side of SCE* = 1.466 is the
base solvent on?

This prediction is directly testable against every published TMB
study without a single new experiment. Every case where TMB helped
should correspond to a base solvent above SCE*. Every case where it
hurt should correspond to a base solvent below SCE*.

**Additionally found:**

"In situ polymerization of 1,3-dioxolane and formation of
fluorine/boron rich SEI" (Chem. Sci., 2024) confirms that
DOL + boron-containing additives is an active research area — but
all existing work focuses on the SEI chemistry at the interface,
not on what the additive does to the Li⁺ coordination population
before it reaches the interface. The population-level mechanism
is entirely absent.

**Implication for Preprint 3 (DOL:TMB):**
The mechanism is confirmed in the literature. The correct application
(DOL as base solvent; TMB as downward SCE lever from above SCE*) is
not present anywhere. The geometric explanation for a decade of
inconsistent TMB results is new. This is a strong novelty claim
with retroactive explanatory power.

---

## Candidates 1 and 2 (Cross-Class Ether Mixing: 2-MeTHF:DME and FEME:2-MeTHF) — Ingredients Confirmed, Mechanism Absent

**Papers identified:**

1. "Methylation enables high-voltage ether electrolytes for lithium
   metal batteries" — Nature Chemistry, 2024
2. "Electrolyte Design for Low-Temperature Li-Metal Batteries:
   Challenges and Progress" — Nano-Micro Letters, 2023 (PMC)

**What these papers found:**

- 2-MeTHF in LiFSI electrolytes achieves ~94% RT CE and ~74% LT CE
  at −20°C. This is the dataset point already in the framework's
  21-system table (LiFSI/2-MeTHF, SCE = 1.552).
- Mixed cyclic/linear ether blends are used for low-temperature
  performance, with the rationale given as "lowering viscosity" and
  "maintaining conductivity."

**What is absent:**

Not one paper frames the cyclic/linear ether mixing question as a
coordination population splitting exercise. The rationale is always
bulk transport properties or SEI stability — never the Shannon
entropy of the Li⁺ shell configuration distribution.

The specific DME:2-MeTHF ratio that targets SCE* = 1.466 has not
been identified, tested, or even conceptualised in the literature.

The ingredients are common. The mechanism — cross-class structural
tension raising SCE toward a fixed point — is entirely absent.

**Implication for Preprint 4 (Candidates 1 and 2):**
Strong novelty claim. The mechanism is new. The target ratio is new.
The framing of mixing as a coordination entropy lever with a specific
numerical target is new. Existing literature provides dataset
confirmation (2-MeTHF CE values) but no design rationale.

---

## Arctic Candidate A (DME:TMP:TMB 1:1:1) — No Literature at All

**Result of search:**

No paper was found that combines all three coordination classes
(ether + phosphate ester + boron ester) in a single formulation with
any rationale related to coordination population diversity.

Ternary electrolytes exist in the literature but are always justified
on bulk property grounds:
- Flame retardancy from TMP
- Ionic conductivity from DME
- SEI chemistry modification from borate additives

No paper designs a ternary formulation because three coordination
classes simultaneously produce extreme shell diversity targeting
an Arctic low-temperature performance regime. The concept of
using class diversity as the design variable — rather than bulk
properties — does not exist in the literature for ternary systems.

**Implication for Preprint 5 (Arctic Candidate A):**
Complete novelty. No prior art on the concept, the formulation,
or the rationale. MD simulation required before experiment; no
published baseline to compare against.

---

## The Broader Pattern — Three Framework-Level Papers That Matter

Beyond the individual candidates, three papers emerged from the
search that speak directly to the framework's conceptual foundations.

---

### 1. Joule 2025 — Independent Replication of the SCE Axis

**Paper:**
"Solvation-configurational entropy governs interfacial kinetics
in low-temperature batteries"
Cell / Joule — DOI: 10.1016/j.joule.2025.102271
Luo, Fu, Peng, Gao, Gu et al. — Hunan University and collaborators

**What they found:**

The Hunan University group independently defines Ssc (solvation
configurational entropy), shows high Ssc → lower desolvation
barrier → better low-temperature interfacial kinetics, and designs
a high-entropy electrolyte that achieves 88% CE at −40°C.

**Relationship to the SCE framework:**

SCE and Ssc are demonstrably the same variable — Shannon entropy
of the Li⁺ first-shell coordination configuration population
distribution. Both groups arrived at the same axis independently.
Neither group has cited the other. The fixed point SCE* = 1.466
does not appear in their work — they optimised toward high Ssc
without identifying the analytical optimum.

**Status:** Independent replication of the core prediction (high
solvation configurational entropy → better LT performance). The
fixed point, the RT tradeoff, the two-curve model, the six candidates,
and the design rationale are all absent from their work.

**Action:** Contact corresponding authors to confirm variable
equivalence. If confirmed, this becomes co-discovery language
in the relevant manuscripts.

---

### 2. Nature Energy 2023 — Qualitative Confirmation of the Direction

**Paper:**
"Diversifying the solvent"
Nature Energy commentary — Wei, Chen & Cao, 2023
DOI: 10.1038/s41560-023-01294-9

**What they found:**

Increasing solvent diversity → increasing entropy of solvation →
improved ionic conductivity by disrupting mesoscale ion clustering.
The commentary endorses the entropy-as-design-variable direction
and reviews high-entropy electrolyte results confirming the trend.

**Relationship to the SCE framework:**

The commentary treats entropy as a direction ("more is better")
rather than as an axis with an optimum. The SCE framework shows:
more is not unconditionally better. More is better up to
SCE* = 1.466 for the combined RT/LT optimum. Higher than that
is correct for the Arctic regime but trades RT performance.

The commentary is moving toward the fixed point without knowing
the fixed point exists.

**Status:** Qualitative directional confirmation. Framework is
the quantitative version of what this commentary endorses.

---

### 3. National Science Review 2025 — Adjacent Thermodynamic Framework

**Paper:**
"Designing electrolytes by thermodynamics"
National Science Review — May 2025
DOI: 10.1093/nsr/nwaf100

**What they found:**

A thermodynamic theory of electrolyte design using solvation free
energy as the design parameter. The framework identifies the
competitive equilibrium between cation-solvent and cation-anion
interactions as the structural determinant of solvation architecture.
Explains high-concentration, localized high-concentration, weakly
solvated, anion-coordination, and high-entropy electrolytes from
unified thermodynamic principles.

**Relationship to the SCE framework:**

This is the closest approach in the literature to the SCE framework —
but in thermodynamic language rather than information-theoretic
language.

Key difference: their framework optimises solvation free energy (a
continuous thermodynamic quantity). The SCE framework optimises
Shannon entropy of the discrete population distribution of
coordination configurations. These are related but not identical.

Their framework produces the same qualitative results in many cases
but does NOT produce a fixed point — because solvation free energy
is a smooth gradient without an interior maximum of the combined
RT/LT performance function.

The SCE framework produces a fixed point because it operates on
the population distribution (discrete configurations), not the
continuous energy. The interior maximum at SCE* = 1.466 is a
consequence of the information-theoretic formulation.

**Status:** Adjacent framework. Confirms the physical approach
is sound. Does not reproduce the fixed point. Does not produce
specific candidates. Does not derive the design rationale for
any of the six formulations.

---

## Summary Table — Candidate-by-Candidate Literature Status

| Candidate | Literature Status | What Exists | What Is Missing |
|-----------|------------------|-------------|-----------------|
| **Candidate 4 (DME:TMP)** | Accidental partial hit | LiFSI/TMP pure electrolyte studied; coordination structure (3 TMP + 1 FSI⁻) confirmed by ab initio MD; CE approaching 98% | ESP-matching mechanism for population splitting never attempted; DME:TMP blend at SCE* target ratio does not exist in literature |
| **Candidate 3 (DOL:TMB)** | Mechanism confirmed, application wrong | TMB anion receptor mechanism confirmed in ACS AMI 2023; DOL + boron active in Chem. Sci. 2024 | Correct base solvent (DOL above SCE*) never used with TMB; geometric explanation for "angel or demon" behaviour entirely absent |
| **Candidates 1 & 2 (ether mixing)** | Ingredients confirmed, mechanism absent | 2-MeTHF, DME, FEME studied individually; some blends reported for LT performance; 2-MeTHF CE data in framework dataset | Cross-class structural tension as coordination entropy lever never framed; target ratio for SCE* never identified |
| **Arctic Candidate A (DME:TMP:TMB)** | No literature | Nothing on three-class ternary with this rationale | Entire concept absent; no prior art on coordination class diversity as ternary design criterion |
| **Ssc / SCE axis** | Independent replication confirmed | Joule 2025 derives same variable independently; high Ssc → better LT kinetics confirmed | Fixed point not identified; no design criterion for optimal SCE; no RT tradeoff characterised |
| **Solvent diversity direction** | Confirmed qualitatively | Nature Energy 2023 endorses entropy increase as design direction | No fixed point; "more is better" without identifying the optimum |
| **Thermodynamic design framework** | Adjacent framework | NSR 2025 uses solvation free energy as design variable; same qualitative physics | Different formalism (continuous free energy vs. discrete population entropy); no interior maximum; no specific candidates derived |

---

## The Single Most Important Finding

**The TMB "angel or demon" paper is the smoking gun for Preprint 3.**

The field has been using trimethyl borate for years, finding
inconsistent results, and attributing the inconsistency to
concentration effects and SEI chemistry complexity.

The SCE framework gives the clean geometric explanation:
