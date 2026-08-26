# Lakatosian First-Principles Topological Chess Engine
## Artifact #2 — Corrected Understanding (Preserved Program Logic)

## 1) Correction Statement
Previous framing that emphasized discrete feature extraction and fixed category assignment is not aligned with the intended Lakatosian first-principles program.

The correct architecture is:
- **not** heuristic bucketing,
- **not** static category scoring,
- **not** combinatorial state storage with lookup-driven policy,
- but a **generative topological class system** where each class object encodes:
  1. intrinsic topology of its piece configuration space,
  2. morphisms to/from neighboring piece-count classes,
  3. derivation operators for board-structure navigation,
  4. ambiguity conditions under which bounded combinatorial search is invoked.

---

## 2) Core Distinction (Essential)
A class file is **not** “evaluate position with metrics.”  
A class file is a **topological data object + derivation engine** that supports two navigation modes:

1. **Intra-class navigation**  
   (move within same piece-count topology over board structure)

2. **Inter-class navigation**  
   (move across piece-count topologies via reductions/trades/captures)

These are the two primary principle vectors.

---

## 3) What a Piece Class Object Actually Is
A piece class object (e.g., `KRvK`, `KRvKp`, `KRvKpp`) must encode:

1. **State manifold definition**  
   The legal board-structure topology induced by that piece multiset.

2. **Topological invariants**  
   Structural truths that remain meaningful under legal transformations in class.

3. **Generators (legal transformation operators)**  
   Move-induced topology transformations available from any state in manifold.

4. **Internal navigation principles**  
   Directional rules for moving toward favorable regions in the class manifold.

5. **Boundary maps**  
   What transitions move this class to adjacent classes (capture/reduction/addition).

6. **Interdependency links to lesser classes**  
   How lower classes constrain or inform interpretation of this class manifold.

7. **Ambiguity surface definition**  
   Regions where principle derivation does not uniquely determine navigation, triggering constrained search.

---

## 4) Ground Truth Role (Precise)
Ground truth (3/4/5 and beyond where available) is not used as a giant path database for runtime retrieval.  
It is used to **derive and validate class-level principles**:

- discover invariant structures,
- validate transition morphisms,
- verify boundary behavior between classes,
- test whether derived principles preserve perfect-play outcome relations.

Then runtime uses these derived principles generatively on current board structure.

---

## 5) Interdependency Principle (Lesser -> Higher)
Higher classes are not independent models.

For class `C_{n+1}`:
- it must explicitly depend on one or more `C_n` classes,
- inherit stable principles from those lesser classes,
- add perturbation logic for added material/topological generators,
- and re-evaluate manifold structure under the expanded operator set.

Example:
- `KRvK` informs `KRvKp`
- `KRvKp` informs `KRvKpp`
- each additional pawn alters manifold geometry and transition boundaries
- but interpretation remains anchored in inherited lower-class structure

This is compositional topology, not flat evaluation.

---

## 6) Runtime Logic (Correct Form)

Given board state:
1. Identify piece class manifold object.
2. Instantiate current board point in that manifold.
3. Apply intra-class derivation operators to determine navigable principle directions.
4. Evaluate available inter-class morphisms (reductions/transitions).
5. If derivation yields a unique direction: execute without broad search.
6. If derivation yields competing valid directions (ambiguity surface): run constrained combinatorial search over those derived branches only.

Search is a local resolver on a derived branch set, not global planner.

---

## 7) Why This Is Lakatosian
## Hard Core
- chess as topological-combinatorial adversarial system
- piece-class manifolds and morphisms as primary explanatory structure
- derivability from first principles over board/piece topology

## Protective Belt
- concrete module implementations
- refinement of invariants, morphisms, ambiguity surfaces
- error corrections against truth without replacing hard-core commitments

Hard core is preserved while belt improves.

---

## 8) What Must Be Avoided Going Forward
1. Treating classes as heuristic scorecards.
2. Treating categories as fixed handcrafted labels detached from manifold derivation.
3. Treating ground truth as runtime brute-force replacement.
4. Treating higher classes as isolated from lower-class dependencies.
5. Running unconstrained search before derivation.

---

## 9) Required Shape of Future Class Files
Each class file should be structured as:

1. `ClassIdentity`
   - piece multiset definition
   - adjacency in class graph

2. `ManifoldDefinition`
   - legal configuration topology for class

3. `InvariantSet`
   - derived, validated structural invariants

4. `GeneratorSet`
   - legal transformation operators

5. `IntraClassNavigationOperators`
   - principle-driven directional derivations

6. `InterClassMorphisms`
   - transitions to neighboring classes with conditions

7. `DependencyContract`
   - inherited principles from lesser classes
   - perturbation terms introduced by added material

8. `AmbiguitySurface`
   - formal trigger conditions for constrained search

9. `ValidationContract`
   - how ground truth is used to test derivations and morphisms

---

## 10) Preserved Statement of Intent
The program is a generative first-principles architecture where:
- board state is interpreted through class topology,
- navigation happens within and across piece-class manifolds,
- lower-class truth-derived principles scaffold higher classes,
- and combinatorial search is used only when topological derivation yields legitimate ambiguity.

This is the intended research program and should be the baseline for all future implementation artifacts.
