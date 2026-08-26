# Lakatosian First-Principles Topological Chess Engine
## Artifact #2 — Full Module Interface Specification (Python-First, Language-Portable)

## 1) Purpose
This artifact defines the implementation-facing interface contract for a modular topological chess engine that:
- derives action from first principles (turn, piece class, board topology),
- composes higher piece classes from lower validated classes,
- and invokes combinatorial search only when topological derivation is underdetermined.

It is designed for:
1. **immediate Python implementation**, and
2. **future migration** to faster runtimes/languages while preserving semantics.

---

## 2) Architectural Commitments (Preserved)
1. **Primary mechanism:** topological navigation through classed state spaces.
2. **Secondary mechanism:** bounded search only on competing topological pathways.
3. **Knowledge growth:** low-piece ground truth (tablebase-validated) composes into higher-piece classes.
4. **State identity:** `(turn, piece_class, board_topology_category)` governs reasoning.
5. **Trading/capture decisions:** treated as transitions across piece classes in topology space.

---

## 3) Global Data Contracts

## 3.1 Core Types

```python name=types_core.py
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Tuple, Set, Any

Color = str  # "white" | "black"
Square = int  # 0..63 (portable representation)
PieceType = str  # "K","Q","R","B","N","P"
MoveUCI = str  # e.g., "e2e4", "a7a8q"

class PhaseTag(str, Enum):
    OPENING = "opening"
    MIDDLEGAME = "middlegame"
    ENDGAME = "endgame"

class OutcomeClass(str, Enum):
    WIN = "win"
    DRAW = "draw"
    LOSS = "loss"
    UNKNOWN = "unknown"

class CertaintyClass(str, Enum):
    FORCED = "forced"             # single topological pathway
    CONSTRAINED = "constrained"   # few pathways, topology dominant
    AMBIGUOUS = "ambiguous"       # needs bounded search
    OPAQUE = "opaque"             # fallback behavior needed

@dataclass(frozen=True)
class PieceCountVector:
    # canonical material signature; include both sides
    # example keys: "wK","wQ","wR","wB","wN","wP","bK","bQ"...
    counts: Dict[str, int]

@dataclass(frozen=True)
class PieceClassId:
    # canonical class token, e.g. "KRvK", "KBNvK", "KQRPvKQRP"
    token: str
    piece_count_vector: PieceCountVector

@dataclass
class PositionState:
    fen: str
    side_to_move: Color
    piece_class_id: PieceClassId
    phase_tag: PhaseTag
    # optional cached board object from chess library
    board_obj: Any = None

@dataclass
class TopologyFeature:
    name: str
    value: float
    meta: Dict[str, Any] = field(default_factory=dict)

@dataclass
class TopologyCategory:
    # class-specific category label derived from features
    label: str
    certainty: CertaintyClass
    features: List[TopologyFeature]
    rationale: List[str]  # human-readable derivation trace
```

## 3.2 Policy / Recommendation Types

```python name=types_policy.py
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional
from types_core import MoveUCI, OutcomeClass, CertaintyClass

@dataclass
class MoveCandidate:
    move: MoveUCI
    topo_score: float
    risk_score: float
    transition_score: float
    rationale: List[str]

@dataclass
class TransitionOption:
    move: MoveUCI
    from_class: str
    to_class: str
    predicted_outcome_shift: OutcomeClass
    confidence: float
    rationale: List[str]

@dataclass
class SearchEnvelope:
    enabled: bool
    reason: str
    depth_cap: int
    node_cap: int
    time_ms_cap: int
    candidate_moves: List[MoveUCI]
    # constraints that prune search manifold
    manifold_constraints: Dict[str, Any] = field(default_factory=dict)

@dataclass
class DecisionPacket:
    chosen_move: MoveUCI
    certainty: CertaintyClass
    principal_mode: str  # topological mode label
    alternatives: List[MoveCandidate]
    used_search: bool
    search_envelope: Optional[SearchEnvelope]
    explanation: List[str]
```

---

## 4) Module Interfaces (Primary Contracts)

## 4.1 `PieceClassModule` Interface
Each piece-class module (e.g., KRvK, KBNvK, KQRPvKQRP) must implement this contract.

```python name=interfaces_piece_class.py
from abc import ABC, abstractmethod
from typing import List, Dict, Any
from types_core import PositionState, TopologyCategory, OutcomeClass
from types_policy import MoveCandidate, TransitionOption

class PieceClassModule(ABC):
    """
    Contract for topological reasoning within a fixed piece class.
    """

    MODULE_ID: str                  # e.g., "KRvK"
    VERSION: str                    # semantic version
    DEPENDENCIES: List[str]         # lower class module IDs
    VALIDATION_STATUS: str          # "unvalidated" | "partial" | "tablebase-validated"
    COVERAGE_NOTES: str             # scope boundaries
    
    @abstractmethod
    def supports(self, state: PositionState) -> bool:
        """Return True if this module applies to state.piece_class_id."""
        raise NotImplementedError

    @abstractmethod
    def extract_topology_features(self, state: PositionState) -> List:
        """
        Compute class-relevant topological features (control, confinement, opposition,
        mobility restrictions, tempo motifs, etc.).
        """
        raise NotImplementedError

    @abstractmethod
    def categorize_topology(self, state: PositionState, features: List) -> TopologyCategory:
        """
        Map features -> category label + certainty + rationale.
        """
        raise NotImplementedError

    @abstractmethod
    def derive_topological_objectives(self, state: PositionState, category: TopologyCategory) -> List[str]:
        """
        Return ordered objectives (e.g., 'cut king', 'improve opposition', 'force corner').
        """
        raise NotImplementedError

    @abstractmethod
    def admissible_move_manifold(self, state: PositionState, category: TopologyCategory) -> List[str]:
        """
        Return a reduced move set consistent with class topology.
        """
        raise NotImplementedError

    @abstractmethod
    def rank_moves_topologically(self, state: PositionState, moves: List[str], category: TopologyCategory) -> List[MoveCandidate]:
        """
        Score candidate moves using topological criteria (not generic eval first).
        """
        raise NotImplementedError

    @abstractmethod
    def evaluate_transitions(self, state: PositionState, ranked_moves: List[MoveCandidate]) -> List[TransitionOption]:
        """
        Evaluate captures/trades as class transitions in topology space.
        """
        raise NotImplementedError

    @abstractmethod
    def expected_outcome_class(self, state: PositionState, category: TopologyCategory) -> OutcomeClass:
        """
        Predict W/D/L class for optimal navigation under this module's scope.
        """
        raise NotImplementedError
```

---

## 4.2 `CompositionModule` Interface
Build higher class principles from lower class modules.

```python name=interfaces_composition.py
from abc import ABC, abstractmethod
from typing import List, Dict, Any
from types_core import PositionState, TopologyCategory

class CompositionModule(ABC):
    MODULE_ID: str
    TARGET_CLASS: str              # e.g., "KBNvK"
    SUBSTRATE_CLASSES: List[str]   # e.g., ["KBvK", "KNvK"]

    @abstractmethod
    def can_compose(self, state: PositionState) -> bool:
        raise NotImplementedError

    @abstractmethod
    def compose_principles(self, state: PositionState, substrate_reports: Dict[str, Any]) -> Dict[str, Any]:
        """
        Merge inherited lower-class principles and perturb by added-piece topology.
        """
        raise NotImplementedError

    @abstractmethod
    def generate_composed_category(self, state: PositionState, composition_report: Dict[str, Any]) -> TopologyCategory:
        raise NotImplementedError
```

---

## 4.3 `TransitionPolicy` Interface
Decide when moving across classes (captures/trades/simplifications) is desirable.

```python name=interfaces_transition_policy.py
from abc import ABC, abstractmethod
from typing import List
from types_core import PositionState
from types_policy import MoveCandidate, TransitionOption

class TransitionPolicy(ABC):
    POLICY_ID: str
    VERSION: str

    @abstractmethod
    def select_transition_bias(self, state: PositionState) -> str:
        """
        e.g., 'simplify', 'preserve tension', 'avoid liquidation', 'force reduction'
        """
        raise NotImplementedError

    @abstractmethod
    def score_transition_options(self, state: PositionState, options: List[TransitionOption]) -> List[TransitionOption]:
        raise NotImplementedError

    @abstractmethod
    def apply_transition_constraints(self, state: PositionState, ranked_moves: List[MoveCandidate], transition_options: List[TransitionOption]) -> List[MoveCandidate]:
        raise NotImplementedError
```

---

## 4.4 `SearchPolicy` Interface
Enable bounded search only if topology is ambiguous.

```python name=interfaces_search_policy.py
from abc import ABC, abstractmethod
from typing import List
from types_core import PositionState, TopologyCategory
from types_policy import SearchEnvelope, MoveCandidate

class SearchPolicy(ABC):
    POLICY_ID: str
    VERSION: str

    @abstractmethod
    def should_search(self, state: PositionState, category: TopologyCategory, ranked_moves: List[MoveCandidate]) -> bool:
        raise NotImplementedError

    @abstractmethod
    def build_envelope(self, state: PositionState, category: TopologyCategory, ranked_moves: List[MoveCandidate]) -> SearchEnvelope:
        """
        Build depth/node/time caps + candidate manifold constraints.
        """
        raise NotImplementedError
```

---

## 4.5 `GroundTruthAdapter` Interface (Tablebase calibration)
Connect low-piece modules to exact truth data for validation.

```python name=interfaces_ground_truth.py
from abc import ABC, abstractmethod
from types_core import PositionState, OutcomeClass

class GroundTruthAdapter(ABC):
    ADAPTER_ID: str
    VERSION: str

    @abstractmethod
    def supports(self, state: PositionState) -> bool:
        raise NotImplementedError

    @abstractmethod
    def query_outcome(self, state: PositionState) -> OutcomeClass:
        """W/D/L truth for supported classes."""
        raise NotImplementedError

    @abstractmethod
    def query_distance_metric(self, state: PositionState) -> int:
        """
        e.g., DTM/DTZ/other supported distance metric.
        """
        raise NotImplementedError
```

---

## 5) Central Orchestrator Contract

```python name=engine_orchestrator.py
class TopologicalEngine:
    """
    Main coordinator:
    1) classify state
    2) choose class module (or composed module)
    3) derive topology category
    4) produce admissible manifold
    5) apply transition policy
    6) optionally bounded search
    7) emit decision packet + reasoning trace
    """

    def choose_move(self, fen: str) -> "DecisionPacket":
        ...
```

## Required orchestration sequence
1. Parse `fen`, detect `side_to_move`.
2. Compute canonical `PieceClassId`.
3. Select module:
   - direct class module if available,
   - else composition module built from substrate classes.
4. Extract topology features.
5. Categorize topology (`FORCED|CONSTRAINED|AMBIGUOUS|OPAQUE`).
6. Derive objectives and admissible move manifold.
7. Rank moves topologically.
8. Evaluate transition options (trade/simplify/preserve tension).
9. Apply transition policy.
10. Invoke search policy only if needed.
11. Return `DecisionPacket` with full rationale chain.

---

## 6) Canonical Piece-Class Taxonomy Rules

## 6.1 Naming convention
- Format: `"WhitePiecesvBlackPieces"` with kings implicit but always present.
- Example tokens:
  - `KRvK`
  - `KBNvK`
  - `KQRPvKQRP`
- Canonical ordering within side:
  - `K Q R B N P` by descending piece hierarchy and multiplicity

## 6.2 Normalization
- Side-symmetric states should map to equivalent normalized IDs where possible.
- Store side-to-move separately; do not encode turn into class token.

## 6.3 Class tiering
- Tier by total piece count (including kings), then by non-king material complexity.
- Use tiers for scheduling module maturity and testing.

---

## 7) Board Topology Category Schema (Initial)

Each module defines class-specific categories; global categories provide common vocabulary:

1. `CONFINEMENT_GAIN`
2. `CONFINEMENT_STABLE`
3. `CONFINEMENT_LOSS_RISK`
4. `OPPOSITION_FAVORABLE`
5. `OPPOSITION_CONTESTED`
6. `FORCING_CHANNEL_AVAILABLE`
7. `FORCING_CHANNEL_BLOCKED`
8. `TRANSITION_OPPORTUNITY`
9. `TRANSITION_HAZARD`
10. `TEMPO_CRITICAL`

Each category record includes:
- feature vector,
- certainty class,
- explicit rationale list.

---

## 8) Feature Library (Portable primitives)

Minimum first-principles features (module-selective use):
1. legal mobility counts (both sides)
2. king-zone accessibility map
3. opposition relation metrics
4. confinement perimeter measure
5. controlled escape-square count
6. forcing-move availability count (checks, skewers, tempo gains)
7. trade/reduction availability and consequence class
8. repetition/stalemate risk flags
9. promotion race and blockage metrics (pawn classes)
10. move-to-dominance delta (control-shift forecast)

---

## 9) Decision Semantics

## 9.1 Choosing move != raw eval maximization
Primary target:
- navigate to topologically favorable region
- improve certainty class (AMBIGUOUS -> CONSTRAINED -> FORCED)
- reduce adversary’s counter-topology

## 9.2 Search trigger semantics
Search is invoked when:
- multiple near-equivalent topological pathways exist,
- transition options conflict,
- tactical volatility exceeds topology-only resolution.

Search is not invoked when:
- a forced topological mode is identified with high certainty.

---

## 10) Validation and Truth Alignment

## 10.1 Module validation ladder
1. `unvalidated`
2. `partial` (unit/regression tests + sampled truth checks)
3. `tablebase-validated` (full supported-space calibration)

## 10.2 Required validation outputs
Per module:
- class coverage declaration
- confusion matrix vs truth (`predicted W/D/L` vs ground truth)
- distance metric error summaries (where available)
- failure archetype inventory (where topology heuristic diverges)

---

## 11) Suggested Python Package Layout

```text name=package_layout.txt
topochess/
  core/
    types_core.py
    types_policy.py
    classifier.py
    topology_features.py
  modules/
    piece_classes/
      krvk.py
      kbvk.py
      knvk.py
      kbnvk.py
      ...
    composition/
      compose_4piece.py
      compose_5piece.py
      ...
  policies/
    transition_policy_default.py
    search_policy_default.py
  adapters/
    ground_truth_syzygy.py
    ground_truth_local.py
  engine/
    orchestrator.py
    tracing.py
  tests/
    test_classifier.py
    test_krvk_truth.py
    test_transition_policy.py
    ...
  docs/
    artifacts/
      LAKATOSIAN_TOPOLOGICAL_CHESS_ENGINE_REASONING_ARTIFACT.md
      LAKATOSIAN_TOPOLOGICAL_CHESS_ENGINE_ARTIFACT_2_MODULE_SPEC.md
```

---

## 12) Tracing, Explainability, and Preservation

Every `DecisionPacket` must include:
1. module selected and version
2. topology category and certainty
3. admissible manifold size vs legal move count
4. transition analysis summary
5. whether search was used and exact envelope
6. human-readable rationale chain (ordered)

This preserves first-principles accountability and supports research iteration.

---

## 13) Performance Strategy (Python-first)
1. Python orchestration, compiled kernels later (Rust/C++/Cython/Numba) for:
   - feature extraction
   - move manifold filtering
   - bounded search inner loops
2. deterministic caching at:
   - class normalization
   - topology feature vectors
   - category derivation outputs
3. strict manifold pruning before search to avoid combinatorial blow-up.

---

## 14) Minimal Viable Build Sequence

1. Implement classifier + piece class normalization.
2. Implement `KRvK`, `KBvK`, `KNvK` modules with truth adapter integration.
3. Add transition policy with simple simplify/preserve modes.
4. Add search policy with strict ambiguity trigger.
5. Add first composition module (`KBNvK` from `KBvK` + `KNvK` substrate principles).
6. Benchmark:
   - move latency
   - manifold reduction ratio
   - truth-alignment in low-piece classes.

---

## 15) Non-Negotiable Invariants
1. Engine must always attempt topological derivation first.
2. Search without topology envelope is disallowed by default.
3. Every module declares dependencies and validation status.
4. Cross-class transitions are first-class citizens (not side effects).
5. Low-piece truth is the calibration anchor for upward composition.

---

## 16) Closing Preservation Statement
This specification operationalizes the first-principles Lakatosian architecture:
- Hard core: topological classed reasoning over chess states.
- Protective belt: modular class files, composition rules, transition/search policies.
- Path forward: truth-validated low-piece modules driving scalable higher-piece navigation with constrained combinatorial search.

The objective is not brute-force best-move extraction over full space, but principled navigation through interdependent topological classes with search only where theory-preserving ambiguity remains.
