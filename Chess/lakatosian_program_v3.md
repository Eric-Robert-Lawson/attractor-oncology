# Lakatosian First-Principles Topological Chess Engine
## Artifact #3 — Implementation Method, Distinctions, and KRvK-to-Scale Construction Plan

## 1) Purpose of this Artifact
This artifact records **how implementation will be done** while preserving the critical distinctions:

1. Start from **ground truth** for low piece classes.
2. Encode class files as **topological principle engines**, not flat combinatorial lookup maps.
3. Build higher classes through **interdependencies on lesser classes**.
4. Enable **direct navigation through loaded class files** at runtime with minimal extra processing.
5. Ensure performance improves as pieces are reduced due to deterministic class-routing and lower-class certainty.

---

## 2) Critical Distinctions (Non-Negotiable)

## 2.1 What this is **not**
- Not a brute-force search-first engine.
- Not merely storing huge move tables.
- Not evaluating positions primarily by generic scalar eval.
- Not dependent on move history for class identity (state is Markovian in board state + turn).

## 2.2 What this **is**
- A first-principles topological navigation engine.
- A class-structured Lakatosian program:
  - **hard core:** topological axioms and class logic,
  - **protective belt:** modular, revisable category/transition/search policies.
- A compositional system where higher classes inherit and constrain from lower validated classes.
- A runtime that chooses moves by navigating topology classes; search appears only under structured ambiguity.

---

## 3) Grand Program Logic (as intended)

Given a position:
1. Determine `side_to_move`.
2. Determine `piece_class` (type + count on each side).
3. Determine `board_topology_category` for that class.
4. Produce class-conditioned objectives.
5. Navigate within class (or across classes via trade transitions).
6. Search only when two or more topologically valid pathways compete.

### Key preserved insight
> Moves before the current board state do not matter for class identity and topology derivation (except rule-state necessities like castling rights/en passant where relevant).  
Therefore class modules can be preloaded and directly invoked, improving speed and consistency.

---

## 4) Ground Truth First: Low-Piece Class Construction

## 4.1 Why start here
Low-piece classes (3/4/5 and partly 6/7 via available resources) provide exact truth anchors.  
This lets us calibrate principles against perfect-play reality before scaling complexity.

## 4.2 Construction sequence
1. Implement class module skeleton (e.g., `KRvK`).
2. Define topological category taxonomy for that class.
3. Define features and rule derivations for categories.
4. Query ground truth adapter (tablebase) for validation.
5. Iterate until category-driven policy aligns strongly with truth.
6. Mark module `tablebase-validated` once acceptance thresholds are met.

## 4.3 What gets validated
- Outcome class correctness (W/D/L).
- Distance-class alignment where available (e.g., DTM/DTZ style measures).
- Category-to-objective consistency (does category produce expected forcing behavior?).
- Transition correctness (captures/trades leading to expected class outcomes).

---

## 5) Class Files Must Encode Principles, Not Just Mappings

Each class file stores:
1. **Category definitions** (topological state classes within piece class).
2. **Feature extraction rules** (how categories are detected from board geometry).
3. **Objective logic** (what structural target to pursue).
4. **Admissible manifold logic** (which moves are topologically coherent).
5. **Transition logic** (when moving to another class is beneficial/harmful).
6. **Validation metadata** (coverage + truth alignment status).

This prevents “table lookup dependency syndrome” and preserves first-principles portability.

---

## 6) Interdependency Model (Lesser -> Higher)

## 6.1 Composition principle
A higher class module must explicitly consume outputs of lesser class modules:
- not by replaying full searches,
- but by importing validated principles and perturbing by added piece topology.

## 6.2 Example
`KBNvK` composes:
- substrate from `KBvK` principles,
- substrate from `KNvK` principles,
- added interaction terms for bishop+knight coupling,
- board category routing for corner forcing motifs and king-routing constraints.

## 6.3 Interdependency contract
Every higher class declares:
- `DEPENDENCIES` (lesser modules),
- `INHERITED_OBJECTIVES`,
- `ADDED_TOPOLOGY_TERMS`,
- `CONFLICT_RESOLUTION_RULES` (if substrates disagree).

---

## 7) Runtime Navigation Without Additional Processing Burden

## 7.1 Preload strategy
At engine startup:
- load registry of class modules,
- load dependency graph,
- load/attach ground-truth adapters,
- precompile category predicates and manifold filters.

## 7.2 Fast-path runtime
At move time:
1. classify state
2. direct dispatch to module/composer
3. category assignment
4. reduced manifold extraction
5. transition check
6. optional bounded search
7. return move + trace

No historical replay, no global recomputation.  
This is the direct methodological processing event you specified.

## 7.3 Why performance should improve with reductions
As pieces are removed:
- class dimensionality declines,
- category ambiguity declines,
- manifold size shrinks,
- dependence on deep search decreases,
- module certainty increases due to stronger low-piece grounding.

---

## 8) Principle Tree and Reduction Semantics

## 8.1 Principle tree
For each class, maintain tree:
- root: hard-core axioms
- intermediate: class-specific topological invariants
- leaves: category-conditioned actionable policies

## 8.2 Reduction semantics
Captures/trades are modeled as transitions in principle tree space:
- evaluate if target class has stronger derivability and safer topology,
- prefer reductions that move toward validated, lower-ambiguity classes,
- reject reductions that enter poorly understood or adversarially volatile classes unless required.

---

## 9) KRvK Starter Blueprint (Concrete Build Pattern)

## 9.1 Module objective
Implement fully validated `KRvK` as first production-quality reference module.

## 9.2 KRvK category set (initial)
1. `CUT_OFF_AVAILABLE`
2. `CUT_OFF_ESTABLISHED`
3. `OPPOSITION_FAVORABLE`
4. `OPPOSITION_NEUTRAL`
5. `KING_BOX_SHRINKING`
6. `STALEMATE_RISK`
7. `TEMPO_CRITICAL`
8. `FORCING_SEQUENCE_AVAILABLE`

## 9.3 KRvK feature set (initial)
- defender king legal mobility count
- rook line control effectiveness
- attacker king support distance
- confinement area metric (defender reachable region size)
- check forcing count
- stalemate risk indicators

## 9.4 KRvK objectives by category (example pattern)
- if `CUT_OFF_AVAILABLE`: prioritize establishing rook barrier with king support
- if `CUT_OFF_ESTABLISHED`: shrink defender region while preserving rook safety
- if `STALEMATE_RISK`: switch to waiting/topology-preserving move class
- if `FORCING_SEQUENCE_AVAILABLE`: follow forced mode unless contradiction detected

## 9.5 KRvK validation gates
- >= target WDL agreement with tablebase sample/full supported set
- bounded distance disagreement thresholds
- no catastrophic misclassification in stalemate-critical categories
- deterministic rationale emission for all chosen moves

---

## 10) Scaling Recipe (After KRvK)

1. Build `KBvK`, `KNvK` with same protocol.
2. Mark all as validated substrate modules.
3. Build first compositional module (`KBNvK`) using dependencies.
4. Expand to 4-piece families with composition-first methodology.
5. Introduce 5-piece with selective category expansion and ambiguity-gated search.
6. Keep strict “topology first, search second” invariant.

---

## 11) Search Placement (Exactly as intended)

Search is a resolver for competing valid topological pathways, not a universal planner.

### Trigger conditions
- category certainty = `AMBIGUOUS`,
- multiple top-ranked moves with conflicting transition consequences,
- tactical volatility exceeds deterministic class resolution.

### Envelope requirements
- candidate moves limited to admissible manifold,
- depth/node/time bounded by category policy,
- no unrestricted tree expansion by default.

---

## 12) Data Persistence and Research Traceability

Each move decision stores:
1. class id
2. category label + certainty
3. objective stack
4. manifold size vs legal move count
5. transition analysis
6. search envelope (if used)
7. final rationale chain

This is essential for Lakatosian iteration (diagnosing where protective belt needs refinement).

---

## 13) Acceptance Criteria for “Architecture Correctness”

The architecture is considered faithful to your intent only if:

1. Ground truth is used first to build low-piece class principles.
2. Class files are principle/categorization engines, not raw mapping dumps.
3. Higher modules declare and consume explicit lesser-class dependencies.
4. Runtime dispatch is direct and class-based (no history-dependent recomputation).
5. Move quality is produced by topological navigation, with search as constrained fallback.
6. Measured compute per move trends downward as piece count drops.

---

## 14) Immediate Implementation Checklist

1. Build module registry + dependency graph loader.
2. Implement canonical class classifier.
3. Implement `KRvK` end-to-end (features, categories, objectives, manifold, transitions).
4. Add ground-truth adapter integration and validation harness.
5. Add decision trace serialization.
6. Add ambiguity-gated bounded search policy.
7. Run KRvK benchmark suite; freeze spec once stable.
8. Repeat for `KBvK`, `KNvK`, then compose `KBNvK`.

---

## 15) Closing Preservation Statement
This artifact preserves the intended grand picture:

- The engine is a **Lakatosian, first-principles, topological navigation system**.
- It is built bottom-up from **truth-validated low-piece classes**.
- Higher classes are formed through **interdependent composition** from lesser classes.
- Class files encode **topological categories and navigational principles**, not mere combinatorial tables.
- Runtime speed and quality improve through **direct class routing** and **ambiguity-bounded search**.
- The end goal is a scalable method for chess solving that mitigates combinatorial explosion through principled topology.
