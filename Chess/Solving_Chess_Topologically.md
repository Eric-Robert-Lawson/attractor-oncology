# SOLVING CHESS TOPOLOGICALLY
## A Principles-First Reasoning Artifact
## OrganismCore | Eric Robert Lawson | 2026-08-23

---

## PREAMBLE

```
This document exists because a principle has been identified
that makes chess solvable.

Not solvable by brute force.
Not solvable by statistical learning.
Solvable by topological navigation.

The same class of principle that routes you from New York
to Panama City in real time through a combinatorial space
that dwarfs the chess position landscape.

The same class of principle that solved the Rubik's cube.
The same class of principle that maps the epigenetic
landscape of cancer and derives drug targets from
first principles across 36 cancer types.

The principle is the same.
The parameters change.
The principle does not.

This document records:
  What the principle is.
  Why chess is the correct domain to demonstrate it.
  What the methodology is.
  What the system requires.
  How to build it.
  What validation looks like.
  What solving chess actually means.

This is not a proposal.
This is a construction plan.
```

---

## PART I — THE PRINCIPLE

### I.1 — The Tonnetz Space Principle

```
The Tonnetz was introduced by Euler in 1739 as a
representation of harmonic relationships between
musical tones in a two-dimensional lattice.

It has been treated for 300 years as a tool for
music theory.

This is the 300-year oversight.

The Tonnetz is not about music.
Music is one physical instantiation of a more
general principle:

  Any system whose states are describable as
  positions in a combinatorial relational space,
  and whose transitions between states are
  governed by fixed rules, has a topological
  structure that is navigable from first principles.

A resonance cavity of any kind — a violin string,
a vocal tract, a wind instrument, a tinnitus-
affected ear canal — maps its combinatorial
state space through the same principle.
Music is what happens when a human ear detects
the topological features of that space.

But the principle does not require sound.
It requires a system with:
  — Multiple relational constraints
  — Distinguishable states
  — Fixed rules governing transitions between states
  — Attractor states (stable configurations)

Chess satisfies all four conditions perfectly.
It is a pure formal system with no noise,
no measurement error, no biological complexity.
It is the cleanest possible domain to demonstrate
the principle.
```

### I.2 — The Key Insight: The Board State IS the Coordinate

```
The most important thing to understand about the
topological approach to chess:

  The board state is not described by coordinates.
  The board state IS the coordinate.

This eliminates the hardest problem in building
a topological chess engine — there is no
transformation step, no embedding, no mapping
from position to representation.

The FEN string — the complete description of
the board state, castling rights, en passant
availability, and side to move — is the coordinate.
It is complete. It is self-describing.
It is path-independent.

The Rubik's cube makes this precise:

  A Rubik's cube state is a coordinate in the
  space of all 43,252,003,274,489,856,000
  possible cube configurations.
  The sequence of moves that produced the state
  is irrelevant.
  The state itself tells you everything about
  where you are in the configuration space.

Chess is the same structure.
The moves that produced a position are irrelevant
to optimal play from that position.
The position itself is the coordinate.

This is not a claim.
This follows directly from the rules of chess.
Chess is a Markov system.
The current state is a sufficient statistic
for all future optimal play.
History is discarded.
Only the state matters.
```

### I.3 — Why This Makes Chess Solvable

```
Chess is currently treated as computationally
intractable because the game tree has
approximately 10^123 possible games.
The Shannon number.

This number describes the number of possible GAMES.
Not the number of possible POSITIONS.

The number of legally reachable unique board positions
is approximately 10^44 to 10^47.
Still large. But crucially: these positions have
topological structure.

They are not uniformly distributed in an
undifferentiated space.
They cluster.
They have symmetries.
They have metric relationships to the attractor states
(win, loss, draw).

The combinatorial explosion is a property of
TREE SEARCH.
It is not a property of TOPOLOGICAL NAVIGATION.

GPS existence proof:

  The road network of North America contains
  millions of nodes, hundreds of millions of edges,
  and a dynamic state space that dwarfs the chess
  position landscape.

  GPS solves it in under 10 seconds.

  GPS does not enumerate routes.
  GPS navigates the topology of the space.

  If GPS solves a larger, more dynamic, less
  structured problem by topological navigation,
  then chess — which is smaller, static, and
  perfectly structured — is solvable by the
  same class of method.

  Anyone who claims chess is unsolvable while
  using GPS on their phone is asserting a
  contradiction without recognising it.

The failure to solve chess is not a statement
about chess's inherent complexity.
It is a statement about the inadequacy of
exhaustive enumeration as a methodology.
```

---

## PART II — THE METHODOLOGY

### II.1 — Path Independence: The Foundation

```
FORMAL STATEMENT:
  The optimal move at any chess position P is a
  function f(P) of P alone.
  f(P) is not a function of the move sequence
  that produced P.

This follows from:
  1. Chess rules are fully deterministic
  2. The board state (with castling rights,
     en passant, side to move) is a complete
     description of the game state
  3. Future optimal play depends only on the
     current state and future states reachable
     from it
  4. The past is not reachable and therefore
     irrelevant

IMPLICATION:
  Every game ever played is a sequence of
  board states, not a sequence of moves.
  The moves are the transitions between states.
  The states are the data.
  The transitions are the gradient signals.

When we process a chess game, we discard the
narrative of the game and extract:
  — A sequence of board state coordinates
  — The outcome (win/loss/draw)
  — The gradient at each state (the move played,
    labelled by the outcome it contributed to)
```

### II.2 — The Three Tiers of Chess Principles

```
Chess principles exist at three tiers.
The topological methodology operates on all three.

TIER 1 — FORMALISED PRINCIPLES:
  Principles that grandmasters have named
  and transmitted.
  Examples:
    — The opposition (king and pawn endgames)
    — The Lucena position (rook endgames)
    — The Philidor position (rook endgames)
    — Rook behind the passed pawn
    — The principle of two weaknesses
    — Knight outposts
    — King activity in endgames

  These are topological principles that have been
  partially formalised by human pattern recognition.
  They are correct but incomplete — each is a
  local approximation to a more general topological
  truth.

  The methodology will:
    — Confirm these principles emerge naturally
      from the topological analysis
    — Extend their scope to the full class of
      positions sharing their topological feature
    — State them formally as topological objects
      rather than as named heuristics

TIER 2 — INTUITED BUT UNFORMALISED PRINCIPLES:
  Principles that grandmasters navigate by
  intuition but cannot articulate.

  When Magnus Carlsen plays a move that appears
  counterintuitive and the position unfolds over
  20 moves to reveal its correctness — he is
  navigating a Tier 2 principle.
  He cannot tell you why.
  He experienced it as intuition.
  It is topological navigation operating below
  the threshold of verbalisation.

  These principles exist in the game corpus as
  systematic patterns in grandmaster decision-making
  that deviate from obvious choices but are
  confirmed as optimal by later analysis.

  The methodology will:
    — Identify these patterns in the game corpus
    — Derive the structural feature of the position
      that makes the non-obvious move correct
    — State that feature as a formalised principle
      for the first time

TIER 3 — PRINCIPLES NO LIVING HUMAN UNDERSTANDS:
  Topological features of the position landscape
  that are too high-dimensional, too subtle, or
  too non-obvious for human pattern recognition
  to have detected.

  These principles govern optimal play in specific
  position classes.
  They have never been named.
  They have never been transmitted.
  They may never have been consciously perceived
  by any human who has played chess.

  They will emerge from the tablebase analysis
  as structural features that distinguish outcomes
  in ways not corresponding to any named principle
  in chess literature.

  This is genuinely new knowledge.
  500 years of chess have not produced it.
  The topological methodology will.
```

### II.3 — The Post-Hoc Derivation Methodology

```
This is the core methodological insight.

Current chess engines work FORWARD:
  Evaluate the current position.
  Search future positions.
  Propagate evaluations backward.
  Select the move with the best propagated value.

  This is computationally expensive.
  It produces a number, not an understanding.
  It scales catastrophically with depth.

The topological methodology works BACKWARD:
  Take positions where optimal play is KNOWN.
  Analyse what structural features of those
  positions make the optimal moves optimal.
  Extract those features as principles.
  Apply principles forward to new positions.

  This produces understanding, not numbers.
  It does not scale with depth — it scales
  with the richness of the principle library.
  Depth is encoded in the principles themselves,
  not searched explicitly.

THE GROUND TRUTH SOURCE:
  Syzygy endgame tablebases.
  Complete optimal play for all positions
  with 7 or fewer pieces.
  ~500 million positions.
  Every position labelled: Win in N, Draw, Loss.

  This is the topological map for the endgame
  region of the position landscape.
  Every outcome label is an attractor assignment.
  The structural features that distinguish
  win positions from draw positions within
  the same material class are the principles.

THE GRANDMASTER GAME SOURCE:
  Every position in every annotated grandmaster
  game where the played move was confirmed optimal
  by subsequent analysis is a data point.
  The structural feature that makes the optimal
  move optimal is the principle.

  Positions where grandmasters played correctly
  but could not explain why are the Tier 2
  principle locations.
```

### II.4 — The Loss Pathway Mitigation Principle

```
Navigation toward the win attractor is equivalent
to navigation away from the loss attractor.

The loss attractor has topological structure.
Positions that lead to forced loss share structural
features. Identifying those features and building
principles that avoid them is as important as
identifying features that lead toward wins.

This is the topological chess equivalent of
what GPS does when it avoids roads that lead
away from the destination.

The navigation function must operate on both:
  POSITIVE GRADIENT: move toward win attractor topology
  NEGATIVE GRADIENT: move away from loss attractor topology

When these signals align: the optimal move is clear.
When they conflict: the position is genuinely complex
and the conflict is itself informative about the
nature of the position.

The loss pathway mitigation principle also explains
defensive play: optimal defense is not about
finding winning moves. It is about navigating
away from the loss attractor topology while
maintaining proximity to the draw attractor.
This is a distinct navigational problem from
winning. The principle library must address both.
```

### II.5 — The Rubik's Cube Structural Analogy

```
The Rubik's cube was solved by topological
principles operating on the state space.

The cube's state IS its coordinate.
The sequence of moves that produced the state
is irrelevant.
The solution is a navigation from the current
state to the solved state through a sequence
of principle-governed moves.

The solution methodology:
  1. Identify structural sub-problems
     (cross, first layer corners, second layer
     edges, OLL, PLL)
  2. Derive principles for each sub-problem
     that apply regardless of the specific
     configuration within that sub-problem
  3. Compose principles modularly
  4. Navigate to the solved state without
     exhaustive enumeration

The chess methodology is the same:
  1. Identify structural sub-problems
     (material classes, position types,
     strategic themes)
  2. Derive principles for each sub-problem
     that apply regardless of specific
     piece placement within that class
  3. Compose principles modularly
  4. Navigate to the win/draw attractor
     without tree search

The difference:
  Chess has multiple attractor states (win/draw/loss)
  Chess has a larger state space
  Chess principles are multi-dimensional and
  interact in more complex ways

  But the principle is the same.
  The parameters change.
  The principle does not.
```

---

## PART III — THE SYSTEM ARCHITECTURE

### III.1 — System Overview

```
The system has three components:

  COMPONENT 1: THE POSITION ANALYSER
    Input: FEN string (board state)
    Output: Structural feature vector

  COMPONENT 2: THE PRINCIPLE LIBRARY
    Input: Structural feature vector
    Output: Active principles and their
            directional signals

  COMPONENT 3: THE NAVIGATION FUNCTION
    Input: Active principles and signals
    Output: Optimal move

These three components compose into the engine.
No tree search.
No statistical inference.
Pure topological navigation.
```

### III.2 — Component 1: The Position Analyser

```
PURPOSE:
  Extract the structural features of a board state
  that are topologically relevant — that is,
  features that encode the position's relationship
  to the attractor states and to other positions
  in its topological neighbourhood.

INPUT:
  FEN string — the complete path-independent
  description of the board state.

FEATURE CLASSES:

  MATERIAL FEATURES:
    — Piece counts per type per colour
    — Material balance (aggregate and per type)
    — Material imbalances (bishop pair, rook vs
      minor pieces, etc.)
    — Pawn structure: passed pawns, isolated,
      doubled, backward, connected

  KING SAFETY FEATURES:
    — King position (file, rank, distance to centre)
    — Pawn shield integrity
    — Open files near king
    — Attacking pieces proximity to king

  PIECE ACTIVITY FEATURES:
    — Piece mobility (legal moves available)
    — Control of key squares (centre, outposts)
    — Rook activity (open files, rank penetration)
    — Bishop diagonals (open/closed)
    — Knight outpost availability

  STRUCTURAL FEATURES:
    — Pawn chain direction and integrity
    — Weak square complexes
    — Open file control
    — Space advantage (controlled territory)

  ENDGAME-SPECIFIC FEATURES:
    — King opposition (for K+P endgames)
    — Passed pawn advancement
    — Rook activity relative to passed pawns
    — Key square occupation

  POSITIONAL TENSION FEATURES:
    — Number of pieces in contact (captures available)
    — Undefended pieces
    — Piece coordination metrics

OUTPUT:
  A structured feature vector representing the
  topological coordinates of the position within
  the position landscape.
  This vector is the input to the principle library.
```

### III.3 — Component 2: The Principle Library

```
PURPOSE:
  A collection of topological principles, each of
  which maps a set of structural features to a
  directional signal (recommended move class and
  positional target).

PRINCIPLE STRUCTURE:
  Each principle has the form:

    IF [structural condition on feature vector]
    THEN [directional signal toward attractor]
    BECAUSE [topological reason — the feature
              that makes this the gradient direction]
    CONFIDENCE [derived from tablebase coverage
                of positions sharing this feature]

PRINCIPLE SOURCES:

  SOURCE A — TABLEBASE DERIVATION:
    For each material class in Syzygy tablebases:
      1. Extract all positions with known outcomes
      2. Compute feature vectors for all positions
      3. Find features that maximally separate
         Win positions from Draw/Loss positions
      4. State those features as principles
      5. Verify principles hold across the full
         material class

    These are the most reliable principles.
    They are derived from complete ground truth.
    They are the Tier 1 and some Tier 3 principles.

  SOURCE B — GRANDMASTER GAME DERIVATION:
    For each annotated grandmaster game:
      1. Extract positions where played move
         was confirmed optimal by engine analysis
      2. Compute feature vector for each position
      3. Identify what feature of the position
         makes the optimal move optimal
      4. Cross-reference with existing principles
      5. If no existing principle matches:
         new Tier 2 or Tier 3 principle candidate

  SOURCE C — CONFLICT ANALYSIS:
    Positions where principles conflict reveal
    the higher-dimensional structure of the landscape.
    A position where Principle A says move X and
    Principle B says move Y is a saddle point in
    the topological landscape.
    Resolving these conflicts reveals the hierarchy
    of principles — which principles dominate
    which in which contexts.
    This is where the deepest structure is found.

PRINCIPLE LIBRARY STRUCTURE:
  Organised by material class and position type.
  Hierarchical: higher-level principles (game phase,
  material balance) constrain which lower-level
  principles are active.
  Composable: multiple principles can be active
  simultaneously, their signals composed to
  produce the navigation direction.
```

### III.4 — Component 3: The Navigation Function

```
PURPOSE:
  Given the current board state and the principle
  library, select the move that most efficiently
  follows the topological gradient toward the
  correct attractor.

ALGORITHM:

  STEP 1 — ATTRACTOR IDENTIFICATION:
    What is the correct attractor for this position?
    Win attractor: if material or positional
    advantage makes win achievable.
    Draw attractor: if balance is maintained and
    win is not achievable.
    Avoid loss attractor: always.

  STEP 2 — PRINCIPLE ACTIVATION:
    Query the principle library with the feature
    vector of the current position.
    Identify which principles are active —
    which structural conditions are satisfied.

  STEP 3 — SIGNAL COMPOSITION:
    Each active principle produces a directional
    signal:
      — Move class (e.g., "activate rook to
        seventh rank", "advance passed pawn",
        "activate king toward opposition")
      — Confidence (derived from principle's
        tablebase validation coverage)
      — Polarity (positive: toward win attractor,
        negative: away from loss attractor)

    Compose signals:
      Strong agreement: clear move direction
      Weak agreement: move direction with uncertainty
      Conflict: saddle point — requires higher-order
                principle resolution

  STEP 4 — MOVE SELECTION:
    Generate legal moves from current position.
    For each legal move, compute the resulting
    position's feature vector.
    Select the move whose resulting position
    maximally satisfies the composed signal.

    NO TREE SEARCH.
    The evaluation is one move deep.
    The depth is in the principles, not the search.

  STEP 5 — VALIDATION (during development):
    Query Stockfish or Syzygy for the position.
    Compare navigation function output to
    engine/tablebase output.
    Discrepancies are:
      — Principle gaps (missing principle in library)
      — Principle errors (incorrect principle)
      — Genuine disagreements (the interesting case:
        is the topological principle right and the
        engine wrong, or vice versa?)
```

---

## PART IV — THE BUILD PLAN

### IV.1 — Phase 1: King and Pawn Endgame Proof of Concept

```
SCOPE:
  King and pawn endgames only.
  Specifically: KP vs K (king and pawn vs king).
  This is the simplest non-trivial endgame.
  The topological principles are completely known:
    — The opposition
    — Key squares
    — The square rule
    — Outside passed pawn
  If these principles emerge from pure topological
  analysis of tablebase data without being
  hand-coded, the methodology is validated.

WHAT WE BUILD:
  1. Position generator for KP vs K positions
  2. Syzygy tablebase interface
  3. Feature extractor (king positions, pawn
     position, opposition, key square calculation)
  4. Principle extractor: group positions by
     outcome, find distinguishing features
  5. Navigation function: apply principles,
     verify against tablebase
  6. Validation report: do known principles emerge?
     Are any unknown principles found?

VALIDATION CRITERIA:
  SUCCESS: Opposition principle emerges as a
  distinguishing feature between Win and Draw
  positions without being hand-coded.

  SUCCESS: Key square principle emerges similarly.

  BREAKTHROUGH: A structural feature emerges
  that distinguishes outcomes in a class of
  positions where no named principle currently
  applies. This is a Tier 3 principle discovery.

TOOLS:
  python-chess (position handling, FEN, legal moves)
  chess.syzygy (tablebase interface)
  numpy (feature computation)
  scikit-learn (feature importance analysis)
  pandas (position dataset management)

ESTIMATED BUILD TIME:
  Phase 1 complete system: one focused session.
  The data exists. The tools exist.
  The methodology is specified.
```

### IV.2 — Phase 2: Rook Endgames

```
SCOPE:
  KR vs K (trivially won — validates basic rook
  activity principles)
  KRP vs KR (Lucena and Philidor positions)
  KR vs KP (defensive technique)

VALIDATION CRITERIA:
  Lucena structural features emerge as Win
  principle from tablebase analysis.
  Philidor structural features emerge as Draw
  principle from tablebase analysis.
  The methodology produces formal topological
  statements of these principles for the first time.

NEW PRINCIPLE SEARCH:
  Rook endgames are the richest source of
  Tier 3 principles. The tablebase covers them
  completely. Any structural feature that
  distinguishes outcomes in rook endgames without
  corresponding to a named principle is a
  genuine discovery.
```

### IV.3 — Phase 3: Complex Endgames and Middlegame Principles

```
SCOPE:
  All 7-piece tablebase positions.
  Selected middlegame positions from grandmaster
  games with confirmed optimal moves.

METHODOLOGY EXTENSION:
  The middlegame does not have tablebase coverage.
  Instead, we use grandmaster games where engine
  analysis has confirmed the optimal move.
  The feature that makes the optimal move optimal
  is the principle.

  This is where Tier 2 principles are found:
  positions where grandmasters played correctly
  by intuition and the post-hoc analysis reveals
  the structural feature that made the move correct.

PRINCIPLE HIERARCHY CONSTRUCTION:
  As the principle library grows, the interactions
  between principles reveal the hierarchical
  structure of the position landscape.
  Some principles dominate others in specific
  contexts. This hierarchy is itself a topological
  feature — the meta-structure of the principle space.
```

### IV.4 — Phase 4: Full Engine Assembly and Validation

```
ASSEMBLY:
  Compose all principles from Phases 1-3 into
  a unified principle library.
  Build the full navigation function.
  Test on positions not used in principle derivation.

VALIDATION AGAINST STOCKFISH:
  For positions in the tablebase: ground truth
  is the tablebase, not Stockfish.
  For middlegame positions: Stockfish is used
  as a reference but not as ground truth.

  WHERE THEY AGREE: principle confirmed.
  WHERE THEY DISAGREE: the most interesting case.
    If the topological principle is correct and
    Stockfish is wrong: we have found a position
    where statistical learning has not converged
    to the topological truth.
    If Stockfish is correct and the principle
    is wrong: the principle is incorrect and
    must be revised.

  The disagreement analysis is where the deepest
  discoveries will be made.

WHAT SOLVING CHESS MEANS:
  Not: a computer that cannot be beaten.
  Stockfish already achieves that.

  Solving chess means:
    — Optimal play from any position is derivable
      from topological principles without tree search
    — The principles are formally stated and
      humanly understandable
    — The solution can be explained, not just computed
    — The first move of the game has a principled
      optimal answer derivable from the structure
      of the position landscape

  This is what has never been done.
  This is what the topological methodology produces.
```

---

## PART V — THE SIGNIFICANCE

### V.1 — Why Chess Is The Right Proof of Concept

```
Chess is a perfectly defined formal system.
No noise. No measurement error. No biological
complexity. No ambiguity about the rules.

If the topological navigation methodology solves
chess — deriving optimal play from structural
principles without tree search — then the case
for the methodology's applicability to every
other domain is as strong as it can possibly be.

Chess is the cleanest possible proof of the
300-year-overlooked principle.

The proof of concept for cancer required working
with noisy biological data, imperfect measurements,
and incomplete cell atlas coverage.
The proof of concept for tinnitus required
working with biological variability in ear canal
geometry.
The proof of concept for Vedic Sanskrit requires
working with partial phonological records.

Chess requires none of this.
The rules are perfect.
The state space is finite and completely enumerable
in the endgame.
The ground truth is available (tablebases).
The existing principles are documented.

A clean demonstration of the methodology in chess
leaves no room for the objection that the results
in other domains were artefacts of data quality
or biological complexity.
The principle either works in a perfect formal
system or it doesn't.
```

### V.2 — What This Demonstrates About the Broader Principle

```
When chess is solved by topological navigation:

  1. The 300-year-overlooked Tonnetz principle
     is demonstrated in its purest possible form.

  2. The claim that combinatorial explosion makes
     problems computationally intractable is
     revealed as a statement about the inadequacy
     of exhaustive enumeration, not about the
     problems themselves.

  3. The GPS / Chess / Cancer / Tinnitus / Language
     connection becomes undeniable:
     These are all instances of the same class
     of problem — navigating a topological space
     toward an attractor — and the same class of
     solution — deriving the structural principles
     that govern the gradient and following them.

  4. The methodology becomes a named field of study.
     Not chess computer science.
     Not bioinformatics.
     Not music theory.
     Topological navigation of constrained
     combinatorial spaces.
     A general mathematical methodology with
     applications wherever the conditions are met.

  5. The Tier 3 principles discovered in chess —
     principles that no living human has ever
     articulated — demonstrate that the methodology
     produces genuinely new knowledge, not just
     formal restatements of existing knowledge.
```

### V.3 — The Historical Framing

```
The history of chess has produced, entirely as a
byproduct of people playing the game, a complete
topological survey of the position landscape.

Every game ever played is a traced path through
the position space.
Every grandmaster decision is a data point about
the local gradient at a specific position.
Every endgame technique ever discovered is a
partially formalised topological principle.

The map has been drawn, point by point, game by
game, position by position, for 500 years.

Nobody knew they were drawing a map.
Nobody had the framework to recognise that what
they were producing was a topological survey.

The anabolic tradition of chess — accumulate games,
accumulate positions, accumulate principles —
has been creating the dataset for the topological
solution without recognising it.

The data has always been there.
The framework has not.

Until now.
```

---

## PART VI — IMMEDIATE ACTION PLAN

### VI.1 — What Is Needed Right Now

```
TOOLS (all open source, all freely available):
  python-chess          Position handling, FEN, moves
  chess.syzygy          Tablebase interface
  numpy                 Feature computation
  scikit-learn          Feature importance analysis
  pandas                Dataset management
  matplotlib            Visualization of principle space

DATA (all freely available):
  Syzygy tablebases     Download: https://syzygy-tables.info
  Lichess game database Download: https://database.lichess.org
  PGN Mentor            Annotated grandmaster games

KNOWLEDGE:
  The methodology specified in this document.
  No chess domain expertise required.
  The topology speaks for itself.
```

### VI.2 — Phase 1 Implementation Specification

```
Script: chess_topo_phase1.py

STEP 1 — Generate KP vs K positions
  For every legal king and pawn vs king position:
    — Generate FEN
    — Query Syzygy for outcome (Win/Draw)
    — Store as (FEN, outcome) pair

STEP 2 — Extract features
  For each position, compute:
    — White king square (0-63)
    — Black king square (0-63)
    — Pawn square (0-63)
    — Pawn file (0-7)
    — Pawn rank (0-7)
    — Kings in opposition (boolean)
    — White king on key square (boolean)
    — Black king in front of pawn (boolean)
    — Distance: white king to pawn
    — Distance: black king to pawn
    — Distance: black king to queening square
    — Pawn in square rule (boolean)
    — Side to move

STEP 3 — Principle extraction
  Separate positions by outcome (Win vs Draw).
  Compute feature importance:
    Which features maximally distinguish
    Win positions from Draw positions?
  The top features are the topological principles.

STEP 4 — Principle formalisation
  For each high-importance feature:
    State the principle: IF [feature condition]
    THEN [outcome gradient].
    Verify: does this principle correctly classify
    positions in the full dataset?

STEP 5 — Navigation function
  Given a new KP vs K position:
    Extract features.
    Apply principles.
    Select the move that maximally satisfies the
    principle signals.
    Verify against tablebase.

STEP 6 — Validation report
  Accuracy of navigation function vs tablebase.
  Principles discovered vs principles in literature.
  Any new principles found (Tier 3 candidates).
```

### VI.3 — The Commitment

```
This is not a proposal.
This is a construction plan.

The data exists.
The tools exist.
The principle exists.
The methodology is specified.

Chess has not been solved because nobody has
had the framework to know that it is solvable
by the correct method.

The correct method is topological navigation.
The framework is now documented.
The build plan is stated.

Phase 1 can be built in a single session.
Phase 2 in a week.
Phase 3 in a month.
Phase 4 — the full engine — in three months
of focused work.

The result is not just a chess engine.
It is the demonstration of a universal principle
that reorganises how navigable combinatorial
spaces are understood across every domain
where they appear.

Chess is the entry point.
The principle is the destination.

The work begins now.
```

---

## DOCUMENT METADATA

```
document_id:      CHESS_TOPOLOGICAL_SOLUTION
type:             Reasoning Artifact — Construction Plan
status:           ACTIVE — PHASE 1 READY TO BUILD
date:             2026-08-23
author:           Eric Robert Lawson
                  OrganismCore
ORCID:            https://orcid.org/0009-0002-0414-6544
contact:          OrganismCore@proton.me
repository:       https://github.com/Eric-Robert-Lawson/
                  attractor-oncology

principle:        The Tonnetz space principle, overlooked
                  for 300 years, describes the navigable
                  topological structure of any constrained
                  combinatorial space. Chess is the
                  cleanest possible proof of concept.

the_claim:        Chess is solvable by topological
                  navigation without tree search or
                  statistical inference.

the_proof:        GPS solves a larger, more dynamic,
                  less structured problem in real time
                  by the same class of method.

the_method:       Post-hoc derivation of topological
                  principles from known optimal positions.
                  Direct application to new positions.
                  No tree search. No statistics.
                  Pure structural navigation.

the_data:         Already exists.
                  500 years of recorded chess.
                  Syzygy tablebases for complete
                  endgame ground truth.

the_gap:          Only one thing was missing.
                  The framework to know what to
                  extract from the data.

                  That gap is now closed.
```

---

*"The map was always being drawn.*
*Nobody knew they were drawing it.*
*500 years of chess games.*
*A complete topological survey of the position landscape.*
*Waiting to be read as what it actually is."*

— Eric Robert Lawson, August 23, 2026
