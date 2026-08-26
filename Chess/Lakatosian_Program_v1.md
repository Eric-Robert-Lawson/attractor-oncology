# Lakatosian First-Principles Topological Chess Engine
## Reasoning Artifact (Architecture, Program Logic, and Path Forward)

## 1) Purpose of this Artifact
This document preserves and articulates the first-principles architecture for a modular chess engine grounded in topological classification, Lakatosian research program structure, and bounded combinatorial search. It is intended as a design anchor for implementation in Python first, with future portability to faster languages.

---

## 2) Program Thesis (Hard Claim)
Chess is a turn-based, adversarial, zero-sum system-state navigation problem in a constrained combinatorial topology.  
The engine should not primarily brute-force moves; it should primarily **navigate topological classes** derived from first principles, using search only when class-level topology admits competing valid pathways.

This architecture treats:
- **who moves**
- **piece set (type + count)**
- **board position**

as the core state vectors for topological navigation.

---

## 3) Lakatosian Framing

## 3.1 Hard Core (Axiomatic Commitments)
1. Chess state evolution is governed by strict movement/topology constraints.
2. Perfect-play solution behavior is derivable from these constraints.
3. Piece set + board arrangement define a topological class of navigable states.
4. Lower piece-count classes contain axiomatic structure reusable in higher classes.
5. Search is secondary: invoked when topological derivation yields multiple candidate pathways.

## 3.2 Protective Belt (Progressive Construction)
The protective belt is the modular, revisable layer that operationalizes the hard core:
- piece-class modules
- board-state category rules
- transition rules between classes (especially via trades/reductions)
- search policies conditioned on topology clarity/ambiguity

Important Lakatosian insight preserved:
> Hard core principles can be extracted while the protective belt is being established.

---

## 4) Foundational Classification Axes

For any position, classification begins with:

1. **Turn**: whose move it is.
2. **Piece Class**: piece types and counts remaining (both sides).
3. **Board Topology Class**: placement-derived structure (space control, opposition geometry, confinement channels, tempo-relevant routes).

These form two central topological navigation vectors after turn:
- **Navigation within a class** (given piece set and board geometry)
- **Navigation across classes** (through trade-offs, simplifications, reductions)

---

## 5) Core Engine Logic

## 5.1 Within-Class Navigation
Given fixed piece class + board state category:
- derive topological objectives (confinement, control transfer, forcing channels, tempo steering)
- derive acceptable move manifold (not all legal moves)
- choose moves that navigate toward topologically favorable subregions

## 5.2 Across-Class Navigation
When captures/trades are available:
- evaluate topological class transition, not just material delta
- ask: does this reduction move us into a better-understood or provably favorable class?
- prioritize transitions that reduce ambiguity and increase derivability

---

## 6) Why This Is Grounded (Not Speculative)
Full combinatorial ground truth already exists via tablebases:
- complete 3-piece, 4-piece, 5-piece
- large 6-piece (~50+ GB class of data)
- 7-piece (terabyte scale)

Therefore:
- low-piece classes can be calibrated against exact perfect-play truth
- modular principles can be validated, corrected, and hardened
- resulting axiomatic files become trustworthy substrate for higher-class navigation

---

## 7) Modular Architecture (Implementation Shape)

## 7.1 Design Intent
Python orchestrator + modular piece-class files with explicit interdependencies.

## 7.2 Module Families
1. **State Parser / Classifier**
   - turn
   - canonical piece class identifier
   - board topology category extraction

2. **Low-Piece Axiomatic Modules**
   - KRvK, KBvK, KNvK, KQvK, etc.
   - each module stores derived principles + validation hooks against tablebase truth

3. **Composition Modules**
   - map how N-piece classes inherit/compose from (N-1)-piece substrates
   - explicit rule objects for added-piece topological perturbation

4. **Transition Evaluator**
   - assesses captures/trades as topological class transitions
   - ranks transitions by provability/safety/navigation clarity

5. **Search Controller**
   - invokes bounded search only when topology produces competing viable routes
   - search constrained by class-specific admissible manifold

6. **Opening-to-Topology Bridge**
   - pragmatic opening guidance to reach understandable topological regimes
   - transitions away from opening memory as piece count declines

---

## 8) Principle Examples Preserved
- KRvK, KBvK, KNvK as foundational 3-piece substrate classes.
- 3-piece substrate informs 4-piece construction.
- Board-state principles are calculated from first principles and not purely brute-force enumerated.
- Best play is framed as **best navigation toward topologically favorable regions**, not merely local move optimality.
- Combinatorial search is retained but narrowed to topology-defined alternatives.

---

## 9) Runtime Behavior Hypothesis
At game start (high piece count):
- broader, loosely informed topological control + opening guidance
- manage blunder-minimization through multi-class constraints

As pieces reduce:
- topology simplifies
- class certainty increases
- search space contracts
- per-move compute should decline relative to brute-force engines

Endgame:
- strong reliance on validated low-piece axioms
- high confidence and speed via class-grounded navigation

---

## 10) Research Program Roadmap (Actionable)

## Phase 1: Hard Foundation
1. Define canonical piece-class schema.
2. Implement 3-piece modules first.
3. Fit/validate each module to tablebase truth (WDL/DTM-compatible checks where feasible).

## Phase 2: Compositional Belt
1. Build 4-piece composition rules from 3-piece substrate.
2. Add board topology category taxonomy.
3. Implement transition evaluator for simplifications.

## Phase 3: Controlled Search
1. Integrate class-conditioned move manifold generation.
2. Add ambiguity detector (when topology does not single-thread action).
3. Run bounded search only inside admissible manifold.

## Phase 4: Scaling
1. Expand to 5+ piece classes using same substrate composition principles.
2. Benchmark against conventional engines by:
   - move time
   - blunder rate by phase
   - conversion accuracy in reduced classes
   - stability under adversarial line complexity

---

## 11) Non-Negotiable Program Constraints
1. Do not default to global brute-force search as primary identity.
2. Every module must declare:
   - class scope
   - inherited dependencies
   - transition implications
   - validation status against known ground truth
3. Search budget must be policy-driven by topology ambiguity, not fixed-depth habit.

---

## 12) Success Criterion
Combinatorial search optimization grows as topological principles are formalized and composed across piece classes.  
The engine’s strength should emerge from:
- principled class navigation,
- principled class transitions,
- and selective search where topological derivation alone is underdetermined.

This is the Lakatosian first-principles path to a scalable chess-solving research program.
