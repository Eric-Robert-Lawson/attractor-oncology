# An Explainable Perfect-Play Research Programme for Chess Endgames

## 1. Programme Overview

This document describes a research programme, in the Lakatosian sense: a
core set of animating commitments (the *hard core*), a body of current,
revisable implementation choices that operationalize those commitments (the
*protective belt*), a direction for future work (the *positive heuristic*),
and a boundary that should not be abandoned under pressure (the *negative
heuristic*). Progress in a programme like this is judged by whether it is
*progressive* — whether revising the protective belt in response to problems
keeps producing new, independently checkable facts — rather than by whether
any single claim within it has already been fully realized.

The current, working result is a complete, exhaustively verified KQvK
(King+Queen vs. King) tablebase whose recorded "best move" at every position
is the outcome of an explicit, two-tier, recursively self-consistent
optimality criterion — not just moves-to-mate, but moves-to-mate *and a
stated, computed reason* for the choice among every move that ties on that
first measure. That result, and the debugging history that produced it, is
documented in full below, alongside the broader programme it sits inside.

## 2. The Hard Core

The commitments that define this programme, and that are not themselves
being tested by any single experiment below — they are the reason the
programme exists:

- **Perfect play should be explainable, not merely computed.** Wherever
  multiple moves achieve identical primary optimality, an arbitrary
  selection among them is an incomplete answer. A principled, stated,
  computable secondary criterion is a more complete one.
- **Optimality criteria should be self-consistent across an entire line of
  play, not evaluated one ply at a time.** A move is only truly justified by
  its criterion if that justification survives being checked against every
  subsequent forced choice, recursively, to the end of the game.
- **The generation of a perfect-play table is itself a first-principles
  problem, not merely a brute-force one.** The same geometric and
  positional understanding that defines what "perfect" means at a given
  position should, in principle, also inform *how much of the position
  space needs to be searched at all* to construct the table — generation
  and criterion are not separate concerns.
- **Scope reduction through mutual optimal play is real and worth pursuing
  directly.** If both sides are assumed to play according to the same
  understood principles, the portion of the state space that is ever
  actually reachable under that mutual assumption is a genuine, potentially
  much smaller target than the full space of all legal positions — and
  that gap is worth trying to characterize and exploit, not just accept.

## 3. The Protective Belt: Current Implementation

The current, specific, revisable choices that operationalize the hard core
for KQvK:

- **The specific secondary criterion:** minimize (queen side) / maximize
  (king side) the opponent's own legal-move count, summed recursively over
  every future point in the line where that side is actually to move.
  This is the programme's current best implementation of "explainable
  tie-breaking" — not a claim that no better secondary criterion could ever
  be found.
- **Search architecture:** forward, recursive, exhaustive search with no
  domain-heuristic move exclusion (every legal candidate genuinely
  considered), accelerated by a memoization/database layer so that
  transpositions into already-solved territory resolve in O(1) rather than
  by re-derivation.
- **Generation ordering:** positions processed in ascending order of
  estimated difficulty, so easier positions are solved first and accelerate
  harder ones through transposition — an initial, coarse implementation of
  the hard-core commitment to first-principles-informed generation, using
  external DTZ estimates as the ordering signal for now.
- **Symmetry:** 8-fold board symmetry (3 rotations, 4 reflections) applied
  to every solved position, with correctly re-derived move notation under
  each transform.
- **Output format:** a plain-text table recording position, turn, best
  move, mate distance, ply count, cumulative opponent-mobility count, and
  per-ply mobility trajectory — independently parseable and checkable by
  anyone, including against external ground truth like Syzygy DTZ.

## 4. Formal Definition of the Current Auxiliary Metric

For a position `p` with side to move `s`, let `cands(p)` be every legal
move. For each candidate `c`, let `v(c)` be one plus the recursively defined
value of the resulting position (0 at checkmate). Let `dir = min` if `s` is
the side with the queen, `max` otherwise.

```
best_v(p)     = dir over c in cands(p) of v(c)
tied(p)       = { c in cands(p) : v(c) == best_v(p) }
own(c)        = legal-move-count of the lone king, evaluated at c,
                regardless of whose turn c represents
bn_cum(c)     = own(c) + bn_cum(chosen continuation from c)   [0 at checkmate]
best_move(p)  = argmin_{c in tied(p)} bn_cum(c)   if s = queen side
              = argmax_{c in tied(p)} bn_cum(c)   if s = king side
```

This is well-defined by the same backward-induction argument that grounds
ordinary minimax on a finite game tree: KQvK terminates in bounded depth
under any legal play, so the recursion has a floor. Ties surviving both
criteria are treated as genuinely, provably symmetric under the stated
metrics, and resolved arbitrarily among themselves.

## 5. Engineering History: Progressive Problemshifts

Each of these was an anomaly the protective belt had to accommodate. None
of them refutes the hard core; each one, once resolved, produced a new
result checkable against independent ground truth — which is what makes
this a progressive rather than a degenerating sequence of revisions.

- **Anomaly: heuristic contamination at recursion leaves.** An early
  implementation substituted a cheap geometric distance estimate whenever
  real search ran out of depth budget, rather than reporting "unproven."
  *Revision:* depth-exhaustion now returns an explicit unproven state,
  identical in kind to "no legal moves," forcing the caller to retry at
  greater depth. *Checkable result:* previously "solved" positions that had
  relied on the estimate could now be re-verified against real search.
- **Anomaly: domain-heuristic move exclusion.** Candidate generation
  excluded certain legal moves (e.g., apparent king retreats) based on
  general chess intuition. At least one real position required exactly
  such a move to reach true optimal play, which the exclusion made
  structurally impossible to find. *Revision:* all exclusion restricted to
  genuine illegality. *Checkable result:* the previously-unreachable
  optimal line became findable and was confirmed via independent fresh
  search.
- **Anomaly: asymmetric tie-break.** The maximizing (king) side had an
  explicit tie-break rule under the secondary criterion; the minimizing
  (queen) side did not. *Revision:* added the missing, direction-correct
  branch, completing the symmetry the hard core requires (principles
  inversion by side to move).
- **Anomaly: depth-vs-cache-shortcut unfairness.** Iterative-deepening
  search accepts the first depth at which any value is found. A database
  cache hit resolves regardless of remaining depth budget; a candidate
  needing genuine fresh recursion is still bound by it. This let a
  cache-accelerated but objectively worse move win a comparison against a
  genuinely better move that had not yet been given enough depth to prove
  itself. *Revision:* a cached value is only trusted for comparison if it
  would also have been provable by real recursion within the current depth
  budget; otherwise treated as unresolved, exactly like real search running
  out of budget. *Checkable result:* reproduced the exact failure on two
  specific positions, confirmed their true values via independent fresh
  search, and confirmed the fix recovered the correct values through the
  real pipeline end to end.
- **Anomaly: unmarked search failures recorded as solved.** Reaching a
  move-count ceiling without finding real checkmate was indistinguishable
  from a genuine solve, and got written to the database as proven.
  *Revision:* an explicit proven/not-proven signal threaded through the
  full call chain, gating what gets recorded.

## 6. Positive Heuristic: Where the Programme Points Next

- **Extend to KRvK, then KPvK**, using the same architecture and criterion.
  KPvK chains naturally at the promotion boundary: a promoting pawn's
  resulting position is looked up directly in the already-solved KQvK
  table, carrying the dual-metric criterion through the material
  transition, not just the raw win/loss value.
- **Build a table of principles, not just a table of positions.** Mine the
  completed KQvK table for recurring geometric motifs — king-wedged-between
  configurations, boundary-proximity classes, race conditions between
  competing containment patterns — and test, against the table's own
  exhaustive ground truth, which motifs reliably predict the optimal move
  class without search. This directly operationalizes the hard-core
  commitment that generation and criterion are not separate concerns.
- **Generative scope reduction.** Investigate whether recognizing a
  region-level pattern — not just an exact transposition — can let table
  construction skip enumerating positions that provably cannot be reached,
  or provably cannot be optimal, under mutual principled play. The closest
  existing formal analog is Relevance-Zone reduction, a current technique
  for exhaustively solving Hex, Go, and Killall-Go via sub-pattern reuse
  rather than full-state enumeration; mapping its formalism onto chess
  board geometry is the most concrete lead for this branch of the
  programme.
- **Principles inversion as a standing generative mechanism.** The
  by-side-to-move inversion of the tie-break criterion (minimize for one
  side, maximize for the other) is already implemented and empirically
  load-bearing; treat it as a template to apply wherever a new principle is
  added, so every future criterion is defined for both sides from the
  start rather than needing a later symmetry fix.
- **Public release as a falsifiable artifact.** Publish the solver, the
  ordered position list, and the resulting table together, so the
  programme's central claim — that a stated, computed criterion can
  legibly and correctly distinguish among tied-optimal moves at scale — is
  independently checkable by anyone, not just asserted.

## 6a. How the Table of Principles, Scope Reduction, and Turn-Inversion Combine

These three items are listed separately above, but the actual claim being
made about them is stronger than any one in isolation: a compact table of
geometric principles, applied with the correct by-side inversion to model
what an optimally-playing opponent will do in response, is the mechanism by
which the reachable-under-mutual-optimal-play space could be navigated
directly — rather than discovered only after the fact by having already
solved the full space exhaustively. If that combination works, it is a
genuinely different mode of construction than "search everything, then
compress" — it is "know enough, in advance, to only need to search the part
that mutual optimal play can actually reach."

**This combination has one specific, unresolved technical requirement that
must be satisfied before it is sound, not just plausible, and it should be
named precisely rather than assumed away:** knowing which positions are
reachable under mutual optimal play ordinarily presupposes already knowing
what optimal play is — the very thing exhaustive search exists to determine.
A principle-based scoping mechanism that quietly assumes this knowledge in
advance is circular, and would not actually avoid the enumeration it claims
to avoid.

The way out of that circularity is the same one alpha-beta pruning and
retrograde tablebase construction already use, and it is not merely
theoretical for this programme — a working instance of it was built and
verified this session. The depth-gated cache fix (Section 5) never assumed
global knowledge of true optimal values; it only ever trusted a value when
it had an *incrementally provable bound* — "this is at least as good as
what real search could also establish within the budget already spent" —
and let correctness refine itself as the budget grew, rather than assuming
the answer up front. A sound version of principle-based scope reduction
would need to follow the identical discipline: a proposed geometric
principle can only be used to exclude a position from search once it comes
with a similarly incremental, checkable proof that the exclusion is safe —
not because the principle "usually holds," but because it is bounded the
same rigorous way the depth-gate bounds a cache hit. This is the concrete
bar the generative-scoping branch of the programme (Section 6, Section 11)
needs to clear to move from a compelling idea to a working result.

## 7. Negative Heuristic: What the Programme Does Not Abandon Under Pressure

- The commitment to *exhaustive, verified* correctness is not traded away
  for speed or elegance. Every acceleration technique (caching, symmetry,
  future pattern-based scoping) must be provably equivalent to full search,
  not merely usually correct — this is precisely the standard the
  depth-vs-cache bug violated, and precisely why it counted as an anomaly
  requiring resolution rather than an acceptable trade-off.
- The secondary criterion, once chosen, is applied *recursively and
  completely* — not just at the root of a search, and not just for one
  side. A criterion that only holds at the first move, or only for one
  side to move, does not satisfy the hard core.
- A found anomaly is resolved by revising the protective belt (the
  implementation), never by loosening the hard core's standard of what
  counts as "explainable" or "proven."

## 8. Novel, Checkable Predictions This Programme Makes

Claims below are stated so they can be checked by someone other than the
programme's author — this is what distinguishes a progressive programme
from an unfalsifiable one:

- For any KQvK position with multiple DTZ-optimal moves, the programme
  predicts a specific one is distinguished by strictly lower (or higher,
  for the king side) cumulative opponent mobility than every other tied
  move — checkable directly against the published table.
- Classical, search-free endgame technique (e.g., the box method) will
  diverge from this table's true optimum in a characterizable, non-random
  subset of positions, identifiable by geometric criteria stated in
  advance of running the comparison — not fitted to the divergence after
  the fact.
- If a geometric motif is proposed as predictive of optimal play, it will
  hold across all symmetric transforms of a matching position and can be
  falsified by finding a single matching position where it fails against
  the table — a concrete, immediate test for any future "table of
  principles" candidate.

## 9. Adjacent and Sibling Programmes

Situating this work relative to existing efforts — not as prior claims to
be distinguished from defensively, but as the closest existing context for
evaluating what's genuinely being added:

- **Syzygy / Nalimov / Lomonosov tablebases** pursue a different hard core:
  completeness and correctness of win/draw/loss and a single distance
  metric, with no commitment to distinguishing among ties. This programme
  depends on and does not dispute their results; it adds a criterion on
  top of an already-solved value rather than re-solving it.
- **Ciancarini & Favini's Kriegspiel tablebase work** (*Theoretical
  Computer Science*, 2010) shares a hard-core element with this programme
  — shortest distance to mate, then a secondary count-based criterion among
  ties — in a different game (imperfect-information chess). The closest
  existing intellectual relative found so far.
- **"R-mobility"** (TCEC) shares the hard-core intuition that minimizing an
  opponent's legal-move count is a meaningful chess quantity, applied
  currently as a tournament tiebreak rather than a per-position tablebase
  criterion.
- **Classical endgame theory** (box method, opposition/key-squares,
  wrong-corner maneuvers) already operationalizes "search-free geometric
  first-principles play" for several specific endgame classes, with proven
  move-count *bounds* rather than proofs of exact optimality at every
  position — the gap this programme's motif-mining direction (Section 6)
  is positioned to measure precisely.
- **Relevance-Zone reduction** is the closest sibling programme for the
  generative-scoping branch of the hard core: same underlying goal (prune
  exhaustive solving below full state enumeration via reusable patterns),
  different game class, not yet — as far as located — applied to chess.
- **The information-theoretic entropy bound** on tablebase compression
  gives this programme's scope-reduction ambitions a concrete target to
  approach rather than an open-ended one: progress is measured as movement
  toward that floor, not as its elimination.

## 10. Current Boundary Conditions

Stated plainly so the programme's current scope is legible: this is, today,
a complete and correct result for one endgame class (KQvK), with a
documented and partially-precedented plan to extend further. It depends on
already-established tablebase results for its own ground-truth validation,
shares its motivating hard core with at least one documented prior effort
in an adjacent game variant, and has not yet demonstrated the generative
scope-reduction branch of the hard core in working form — that remains the
programme's most important open frontier, not a completed result.

## 11. Open Problems / Next Empirical Tests

- Verify the depth-gated cache fix (Section 5) with the same empirical
  rigor on the maximizing (king) side that it received on the minimizing
  (queen) side — flagged as theoretically distinct and not yet separately
  tested.
- Validate the full corrected table against Syzygy DTZ at scale, from a
  completely wiped database, not just on spot-checked positions.
- Run the reordering experiment (solving in the order the landscape itself
  reveals as easiest, rather than raw DTZ order) and measure whether it
  converges after one pass or needs iterating.
- Read the Relevance-Zone paper in full and attempt a concrete mapping of
  its formalism onto chess board geometry, as the leading candidate for
  turning the generative-scoping hard-core commitment into a working
  result.
