# Boundary Board States: The Choice-Force Boundary of Perfect Play

**This is recorded as the primary objective of this project, not a secondary
analysis — stated explicitly and preserved as a fact about the project's own
priorities, not something this document should soften.** A separate, honestly
reasoned assessment of this idea's actual scope and limits is kept in §13,
clearly labeled as a distinct, independently-argued position rather than folded
into or diluting the objective stated here.

A record of a concept developed across this project's own conversation history,
preserved here precisely because it went through real revision — an initial
framing that didn't hold up, a correction, a genuine gap in an earlier draft of
this very document (the candidate-list/incremental-verification workflow in §3,
first flattened into "ruled out" when only one of its component claims actually
was), a second gap (the scoping mechanism in §10 was missing entirely from the
prior draft), a third addition (retrograde seeding from already-decisive
positions as a second, complementary candidate-generation method, now in §2), a
fourth refinement (separating generation from validation within that same
method, also in §2), a fifth addition (the walk's termination guarantee, both
its unconditional in-material form and its broader, conditional form tied to
chess's own unproven start-position status, now also in §2), a sixth addition
(open, distributed resolution as the practical path to scale, now its own
section, §12), a seventh addition (a fifth generation method, construction
from race-condition primitives — Method 5 in §2, recorded in full in
`boundary_board_state_topology.md` §8 — together with caveats on Uri's rule's
status and on wrong-colored versus opposite-colored bishops, in §2, §8 and §9),
an eighth addition (locks generalized from frozen pawn walls to duty-based
locks, composition as complete mutual interdependence, the scope of
verification as an interference closure, and a proposed local/global
distinction — recorded in `boundary_board_state_topology.md` §5a and §8f-§8j,
with the definitional discrepancies noted in §14 here), a ninth addition
(population-first primitive extraction as the methodology's own justification,
§8m; composition identified as the specific mechanism that extends a
checkers-style forward-search floor past Syzygy's own storage ceiling, §7a;
Uri's rule checked against its own published formula and found to decompose
into an already-catalogued race primitive composed with a color-lock fact,
§8n; and the Možina et al. relation refined from a rejected claim of failed
attempt into a structural-blindness claim that does hold, given their material
and goal — all in `boundary_board_state_topology.md`), and the corrected,
complete version below. Kept separate from
`syzygy_improvement_proven_results.md` and the other Claim A/B/C documents because
this is a genuinely different kind of object: those are about distance and escape
count; this is about win/draw/loss pattern-matching at the single precise point
where a losing side's last real choice exists.

> **Where this stands, and how to resume (read this first).**
>
> - **Two files.** This one is the construction pipeline and the record of the
>   idea. `boundary_board_state_topology.md` is the abstract layer: why a
>   confirmed boundary state's status cannot be invalidated from upstream,
>   start-reachability, the scaffolding strategy (§7 there), and the
>   composition-from-primitives idea (§8 there).
> - **The concept in one line.** A boundary state is a drawn position in which
>   the side to move has exactly one non-losing move (§1, condition 2).
> - **Four external anchors, each with its own honest status.** Checkers' weak
>   solution — a tight structural match to the scaffolding strategy (§9 here,
>   §7 there). Možina, Sadikov and Bratko's KBNK work — a looser family
>   resemblance, for Method 3 (§9). The zugzwang statistics of Althöfer and
>   Walter and of Haworth — related to, but a different definition from,
>   boundary states (§9). Uri's rule — a real geometric proxy for a
>   race-condition fact, but a heuristic with documented gaps, so it informs
>   candidates and never supplies certainty (Method 4; topology §8d).
> - **Standing rules.** Heuristics may order and propose candidates, never
>   supply certainty. The record of an idea and the assessment of it are kept
>   separate. What upstream discoveries can invalidate is a catalog's claim to
>   be complete, never a confirmed boundary state.
> - **The live thread.** Composition from race-condition primitives (Method 5):
>   choose one concrete placement of frozen, color-locked material, write out
>   the tempo timings for every route that breaks it, and decide whether the
>   target is a proven wide fortress or an exactly-one boundary state (topology
>   §8d, §9). The thread has since been generalized (topology §8f-§8j): a lock
>   is a *duty*, of which a frozen pawn wall is only the degenerate case; a
>   boundary state is a set of mutually dependent locks whose joint
>   compatibility is a constraint system checked as a whole; and verification
>   is scoped by an interference closure rather than run over the composed
>   position. The immediate next step is to write one primitive (the rook-pawn
>   edge race) as an assume-guarantee contract (topology §8i-§8j).
> - **Open definitional decisions (read §14 before building).** Whether the
>   definition of a boundary state is purely local or upstream-aware; what
>   "local versus global" is meant to distinguish (topology §5a, including
>   the later refinement in §5a-iii). A conflict between two predictions about
>   KBNK/KQK-type materials was flagged and then resolved (§14 item 2, kept in
>   the record with its correction). The first two are not settled, and
>   construction targets depend on them.
> - **Status of the work.** The only code in this document is the
>   `is_boundary_state()` sketch in §6, checked for syntax but not yet run
>   against real Syzygy or project data. A guarded revision is proposed in the
>   §6 addendum (syntax-checked only). Everything else is design.

---

## 1. The core concept, in the terms it was first stated

There is a boundary between "choices" made to *reach* a boundary board state, and
the forced sequence that follows it if the drawing option isn't taken. Everything
before this boundary is choice — a path being steered toward the boundary.
Everything after, if the draw isn't taken, is forced toward mate. The boundary
itself has a specific geometry: **discontinuous**, not continuous — continuity
isn't a well-defined vector for this kind of transition. There is a discrete,
localized upstream/downstream structure: certain board positions flow *toward*
these draw positions, and if a forced win/loss exists downstream, there necessarily
exists a point upstream where a choice could have been made to force a draw
instead, wherever in the game that choice actually sits.

**The three defining conditions, stated exactly as first given:**

1. White and Black can both still win based on piece combination alone, where
   reduction from sacrificing will force a draw.
2. Either Black or White is choosing a draw, and the other side has no options
   available that don't immediately put them past a boundary state into a forced
   sequence — or that also become a draw.
3. Every other upstream position for both Black and White's choices will also
   route into a draw or a loss — otherwise, if a forced win were possible, that
   choice would have been made instead.

**On which draws count.** Draws that emerge as stalemates arising from
sub-optimal play are excluded — a choice could have been made upstream not to
allow the stalemate, so these are not genuine boundary states. But a position
validated (via Syzygy or full-landscape data) as a draw through **forced
repetition**, not stalemate, is a mathematically proven boundary board position.
The final board state reached after repetition sets in, for both sides, is itself
the boundary position — the furthest downstream point, since everything before it
was the approach toward it.

**What a boundary board state is not**: it is not a property of *how* a position
was reached, and it does not require the position to have arisen from optimal
play upstream in some real game. It is a pure, self-contained property of the
position itself — the same way this project has always insisted a full landscape
must be correct for every legal position, not just ones a perfect game would
produce.

---

## 2. Candidate generation — five methods, complementary, all feeding the same verification step

This is the pragmatic, buildable core of the proposal. Five distinct ways to
generate boundary-state candidates are recorded here — a structural filter over
a material's whole space (Method 1); a targeted, seeded search starting from
positions already known to be decisive (Method 2); structural pattern analysis
over the population of states already confirmed (Method 3); and stratified
ablation over defending material, deliberately varied one piece at a time
(Method 4); and construction from race-condition primitives (Method 5). All
five produce the same thing: a candidate needing the same
exactly-one-non-losing-move check (§1, condition 2) before it counts as a
confirmed boundary state, not a replacement for that check — though for
Method 5 that check can take the form of a complete case-list argument.

### Method 1: the material-and-check-potential filter

**The insight from parsing condition 2 directly**: for condition 2 to even be
*possible* at a given position, the position must meet certain structural
requirements first. Board positions that cannot possibly satisfy condition 2 can
be excluded entirely, based on piece combination alone, and on specific
arrangements within a piece combination that don't satisfy the constraint. This
is explicitly **broader than plain insufficient material**: the necessary
condition is *sufficient material to theoretically win* **and** *real check
potential* — a position with sufficient material but literally no way to give or
block a check is pruned from the candidate list on this basis, a genuine, distinct
expansion beyond the already-recognized insufficient-material case.

**What this produces**: not the exact boundary-state set directly, but a
**candidate list guaranteed to contain the full, true set of boundary board
states, plus some non-boundary "garbage."** The garbage gets eliminated by
exhaustive verification, working through a finite list rather than the full
landscape — explicitly framed as bounded by the same limitations Syzygy itself
has (finite piece count, finite practical scale), not a way around that ceiling
but a way to work more efficiently within it.

**Incremental, partial construction**: because the candidate list is finite, it
supports iterative refinement — work through it, prove or disprove entries, and
the unresolved pool shrinks with each pass. A partial catalog, worked down from
the full candidate set, is itself useful before the whole thing is complete.

**Distributed, modular computation**: proposed as a further extension — a public,
open-source contribution network, approached modularly, where individual
candidate positions are computed to full exhaustion one at a time, potentially
distributed across many independent contributors rather than one machine.

**The batch-size-1 analogy, as originally offered in support of this**: in
`run_full_sweep.py`, a batch size of 1 often resolves a large fraction of
positions quickly, with larger batch sizes becoming safe only once the database's
own accumulated knowledge starts doing most of the work. Offered as intuition for
why incremental, candidate-by-candidate proving might behave similarly well. (§4c
records where this specific analogy was directly challenged.)

### Method 2: retrograde seeding from already-decisive positions

**The starting point**: any position with a finite, proven distance-to-mate — a
confirmed decisive win or loss. By definition, this means the losing side has
**zero** non-losing options there; if even one existed, the position wouldn't be
decisive at all, it would be drawn.

**The walk**: from that seed, proceed backward via retrograde analysis — the same
minimax-inversion this engine's own `classify()` already performs, not a forward
trace through one played-out game, and not subject to blunders, since every
position visited is reached by strict backward induction over the whole
decisive region, not by anyone's actual move choice. Track, at each ancestor,
how many non-losing options the defender has. Every position strictly inside the
decisive region has zero, by the same definitional fact as the seed. The **first
ancestor where that count becomes nonzero** is the exact transition point where
the drawn region begins.

**Why this is a candidate, not a guaranteed hit — stated precisely, not as a
hedge**: the transition could land on exactly one non-losing option (a true
boundary state, condition 2 satisfied exactly) or on several (the edge of a
wider drawn region — a real, useful position, but not itself a boundary state
under the strict definition). Distinguishing these is the same one-line
verification any candidate needs, regardless of which method produced it: count
the non-losing options directly.

**Why this is genuinely complementary to Method 1, not a replacement for it**:
Method 1 narrows the *whole* state space of a material via a structural filter,
then still requires exhaustive work across everything that survives it. Method 2
is seeded and targeted — it starts from one specific position already known to
be decisive and walks directly to one specific transition point, with no need to
first enumerate or filter a broader space at all. It's directly buildable right
now, seeded from Syzygy's own already-solved decisive positions (up to 7 pieces)
or from any material this project has already solved — reusing the existing
retrograde machinery as-is for the walk itself, with no new solving technique
required, only a new way of directing an already-verified capability.

**Refining the walk: generation and validation are two separate operations,
and only one of them needs value data.** *Generating* a retrograde ancestor —
finding every position from which some legal move reaches a given position — is
purely mechanical, rules-based, requiring no value information at all; this is
why the walk doesn't need perfect play to produce candidates in the first place,
only real chess rules run in reverse. *Validating* whether a specific traced
ancestor-to-child edge represents a step real optimal play could actually have
taken is a separate question entirely, and it does require value data —
specifically, whether the traced move is among the ancestor's own best
(tied-optimal) choices, not strictly worse than some alternative that mover had
available. Determining this correctly needs the ancestor's *entire* set of legal
children's values, not just the one edge being walked — the same fixed-point
requirement `classify()`'s own backward induction already depends on for every
position it resolves, applied here to one ancestor at a time rather than a whole
material.

Where this can be checked with certainty: a direct Syzygy WDL probe, if the
ancestor's own material — including, for retrograde steps that undo a capture,
the larger material with the extra piece restored — is itself Syzygy-covered; or
this project's own already-solved data, for any material already fully
classified. A traced move that would turn what was otherwise a preserved
advantage into a handed-away one for the ancestor's own mover — an ancestor
confirmed to have been a forced win, where the specific traced move alone
inverts that into a forced loss for its own side — is not a path real optimal
play would take, and the branch is discarded, not carried further backward. A
branch confirmed to preserve the ancestor's own advantage, or to already sit
inside the pre-boundary drawn region itself, is validated and worth continuing
to walk.

Where the ancestor's material isn't yet resolvable by either source — too many
pieces for current Syzygy coverage, not yet solved by this project — the branch
is neither confirmed nor discarded. It becomes a deferred candidate, held in the
unresolved pool until that larger material becomes checkable, rather than
assumed valid or thrown away for lack of an immediate answer.

**One further clarification worth stating explicitly, since it's easy to
conflate**: condition 2's "exactly one non-losing move" is a claim about a
position's own immediate choice, not about how many distinct ways the resulting
draw can subsequently unfold. A boundary state's single surviving move can
itself lead into a position with several different viable repetition cycles or
fortress continuations — that multiplicity doesn't affect whether the position
upstream of it qualifies as a boundary state, since condition 2 is evaluated at
that one position's own move count, not at the downstream branching of however
the draw plays out afterward.

**A structural guarantee about where the walk terminates — two distinct claims,
one unconditional and immediately useful, one conditional and currently out of
computational reach. Kept separate deliberately, since conflating them would
overstate what either one actually establishes.**

**The unconditional version, true within any single, already-fully-classified
material, no speculative premise required**: because a material's own landscape
is finite and completely known once solved, Method 2's validated walk from any
decisive position, followed backward through purely decisive ancestors, cannot
continue forever and cannot loop — it must terminate, in finitely many steps, at
one of exactly two places: a genuine transition point within that material (a
real boundary-region candidate), or the material's own outer edge — a position
whose only validated predecessors require adding a piece back, crossing into a
larger material entirely. There is no third possibility. This doesn't reduce the
cost of the walk itself — each step still needs the same full-children
validation established above — but it does guarantee the walk always produces
something useful: a real candidate, or an explicit, well-defined handoff to a
larger material, never an unproductive infinite regress within the material
already in hand.

**The broader, conditional version, tied to a premise this project doesn't get
to assume proven**: chess's own starting position has never been proven
decisive, and is widely believed — not proven — to be undecided. If that belief
holds, any validated retrograde chain connecting a truly game-reachable decisive
position back to the actual 32-piece starting array must cross a transition
point somewhere along the way, for the same reason the endpoints differ in
category and a finite chain of individual moves connects them — a real,
structurally sound argument, *conditional on an assumption nobody has proven*.
Even granting the assumption, nothing about it guarantees the transition sits
within any material currently computable by Syzygy or this project — it could
require crossing multiple, successively larger, currently-unsolvable materials
first. The wall is guaranteed to exist somewhere along the full chain; its
location within reach of anything buildable today is not guaranteed at all.

**What this offers toward the open cost-mitigation question (§13), stated
plainly rather than left to be inferred**: neither version reduces the per-step
cost of the validated walk — the full-children requirement is untouched by
either. What the unconditional version provides is confidence that the walk, run
within any material already in hand, always resolves to something well-defined
rather than running unproductively — real architectural understanding, not a
computational shortcut. Whether that understanding leads to an actual cost
reduction remains open.

### Method 3: geometric/topological population analysis over confirmed boundary states

**Explicitly a candidate-*prioritization* method, not a certainty method** — the
same standing distinction §1 of `boundary_board_state_topology.md` holds
throughout applies here without exception. Nothing in this method ever confirms
a boundary state; it only proposes which unverified candidates to try first,
using structure discovered in the population already confirmed with certainty.

**The idea, stated precisely.** Once the catalog holds enough confirmed
boundary states to work with, treat that population as data: for each
confirmed state, record structural features — king distance and coverage,
whose move, piece mobility, material combination, distance from any material's
own extremes — and look for what correlates with boundary-state-hood across
the population. A useful extension beyond single-position features: treat each
confirmed state as the root of its own local cone — its own upstream ancestry,
out to wherever it connects via the scaffolding strategy
(`boundary_board_state_topology.md` §7) — and compare the *shape* of many such
cones to each other, not just static features of the states in isolation.

**Why adding pieces to a confirmed boundary state cannot be assumed to
preserve its status — a specific mechanism, not a general caution.** A single
counterexample is enough to rule this out as a certainty method: take a
confirmed boundary state realized through stalemate (§8), and add a piece for
the opponent that blocks an existing line of attack, opening a previously
covered escape square. The stalemate — and with it the entire drawing
mechanism — can vanish, even though the added piece had no legal move of its
own and threatened nothing. Blocking, discovered effects, and escape-square
control are ordinary chess mechanics, not rare edge cases, so no category of
"inconsequential" addition can be assumed safe. Every candidate generated this
way still requires the full verification any other candidate does (§1
condition 2, checked directly) — this method changes which candidates get
tried first, never what counts as confirmed.

**Real, if loose, precedent for mining tablebase-derived data this way** — see
§9's discussion of Možina, Sadikov, and Bratko's KBNK research program, and
the honest distinction drawn there between their specific technique (ABML rule
induction, over the whole tablebase, toward human-readable strategic goals)
and this narrower one (structural feature correlation, over the boundary-state
population specifically, toward candidate prioritization). Family resemblance
in the general approach — tablebases contain more structure than raw values,
worth mining deliberately — not a claim that the specific method matches.

### Method 4: stratified ablation over defender material, testing a geometric/topological distinction directly

**A precondition this method depends on, stated first because it corrects
something in this document's own earlier build plan (§11).** Whether a
boundary state can exist in a material at all depends on that material having
*genuine drawn territory* somewhere in it. A material that is a forced win
*everywhere* against a bare king — KQvK, KRvK, KBNvK among them — has no
drawn positions to sit at the edge of, aside from stalemate reachable only
through a prior blunder, already excluded by §1's own definition. So these
materials are not merely unlikely candidates for Method 2's seeding — they are
**structurally empty of boundary states**, by the same logic that makes
insufficient material an unconditional draw in the other direction. §11's
build plan listed KQvK and KBNvK as Method 2 targets; that needs correcting,
not just noting as a caveat. **This has a satisfying, independent check
already sitting in §9**: Haworth's exhaustive enumeration found zero mutual
zugzwangs in KBNK specifically. That's not a coincidence next to this
document's own claim of zero boundary states there — both are the same
underlying fact (a total forced win has no drawn-adjacent structure of any
kind) observed from two independent directions, one already published, one
derived here.

**The distinction this method is actually built to test — named precisely,
not left as an intuition.** A position's *geometric* properties — piece
placement, king distance, coverage, material count — are computable directly
from the board, with no search required. Whether that position is a boundary
state is a *topological* property in the sense already established throughout
this document: it depends on the full recursive value of every legal child,
not on anything locally readable from the position alone (§4b, and the same
argument that has ruled out every proposed geometric shortcut in this
document's own record). The open, empirical question Method 4 targets: do
geometrically different positions that are topologically the same kind of
boundary state share any discoverable structural signature — not a shortcut
around full verification, but a pattern worth testing for.

**A real, independent, already-documented instance of exactly this
distinction — checked directly, not asserted.** Wrong-colored-bishop-with-
rook-pawn material is not unconditional insufficiency the way KB vs K or KN
vs K is — it's won or drawn depending on a genuine race between the two
kings for a specific corner, decided by tempo and distance, not by material
alone. Confirmed directly: a defending king cut off from the safe corner
loses; one that reaches it draws, and real analysis of such positions
identifies single correct moves among alternatives that surrender the win —
the same narrow-margin shape this whole document is built around. Chess
programming literature states plainly that this exact ending "need[s] to be
encoded in the evaluation function, as search would not evaluate [it]
properly" — engines get it wrong without special-cased knowledge — and there
is already a real, existing example of the kind of geometric proxy Method 4
is looking for: "Uri's rule," a Chebyshev-distance formula comparing the
defending king, the pawn, and the attacking king (with a tempo adjustment for
whose move it is) that approximates this exact topological fact from cheap,
purely geometric inputs. This is direct, human-verified, centuries-old
precedent that the distinction Method 4 targets is real, and a genuine
validation target once tooling exists: does a formal boundary-state scan
correctly reproduce the already-known wrong-bishop-rook-pawn boundary.

**A necessary caveat on that precedent's own status.** Uri's rule is a
heuristic with documented gaps, not an exact boundary: the source that gives it
lists a separate "too far" case where the stronger side wins, and notes the
heuristic depends on recognizing the draw even after the stronger side
sacrifices its bishop. It is a good informer of candidates and a good
validation *target* only when checked against the exact tablebase boundary for
that small material, not against the formula itself (see
`boundary_board_state_topology.md` §8d).

**The experimental design, stated precisely — a controlled, factorial
ablation, not an unstructured search.** Hold the attacking material fixed and
vary the defending side's material one component at a time — KQvK, then
KQvNK, KQvBK, KQvBNK — each addition tested for whether it introduces genuine
drawn territory, and if so, what boundary states emerge there. This isolates
the effect of each added piece rather than conflating several at once, the
same logic behind ablation studies generally. Each material in the sequence
still requires full construction and verification (§2's Methods 1 and 2,
unchanged) — this method decides *which materials are worth that full
construction effort*, using structure discovered from what's already
confirmed, exactly the same prioritization-not-certainty role Method 3 plays,
applied one level up: not which candidate position to try next, but which
*material* to build out next.

### Method 5: construction from race-condition primitives

**Generates candidates by construction, with a constructive argument as its
certificate.** Compose small structures whose values are already established —
the catalog in `boundary_board_state_topology.md` §8b — using frozen,
zero-legal-move pieces treated as landscape because disturbing them provably
loses, so that the composed position's value follows from the primitives' own
race-condition arguments. The full record of the idea, the six-condition proof
obligation, and the separate assessment live in that document's §8; they are
kept there because the reasoning is structural rather than procedural.

**How it differs from Methods 1-4.** The others propose candidates that a
search or lookup then confirms. This one's confirmation can be a complete
case-list argument, with a small obstacle-board solver as insurance rather than
as a gate. That holds only if the case list is genuinely complete, and the
standing rule is unchanged: nothing here is certain because it looks
convincing, only because the argument closes every case.

**Generalization recorded after Method 5 was first written** (full record in
`boundary_board_state_topology.md` §8f-§8j). Frozen, zero-legal-move pieces
are one special case of a wider class: a lock is a *duty* — a piece committed
to a post, movable in principle, but any move that abandons the duty loses
under a race condition. Composition then means several such duties held in
mutual dependence, certified by a per-resource accounting and a
non-interference argument rather than by search over the composed position.
The certificate remains a proof; what changes is where the exhaustive work
sits (exact derivation of each primitive from small material, and validation
of the composition rules on small instances). The standing rule is unchanged.

---

## 3. The actual caching-architecture proposal, restored here at full precision

An earlier draft of this document flattened this into a single rejected claim.
It deserved better, because it's actually two separable claims with different
soundness — restored here as the five-point workflow originally described:

1. **Use §1's three conditions as necessary, not sufficient, filters** to build a
   finite candidate pool for a given material — sufficient material on both
   sides, plus whatever structural setup makes condition 2 possible at all. This
   pool is guaranteed to contain every true boundary state, plus false positives
   that fail full verification.
2. **Fully, exhaustively verify each candidate individually** — no shortcut
   inside this step. This is where the real cost lives.
3. **Cache each candidate's proven verdict — boundary or not — permanently, once
   resolved.** This is the crux, and it is a *different* claim from "avoid
   caching the full landscape at all." The distinction, precisely: the original
   proposal was read, in an earlier draft, as "don't build the full landscape's
   memoized state, recompute subtrees on demand whenever revisited" — that
   version genuinely fails (§4b). This version is the opposite shape: do the
   full, normal, single memoized computation exactly as the existing pipeline
   already does, and simply choose *not to persist* the non-boundary majority of
   that computation to the long-term output file — write out only the small,
   finished, permanently-true verdict for whichever positions turn out to be
   boundary states.
4. **The candidate pool shrinks monotonically as verification proceeds** — real,
   visible, incremental progress, even though the total work to get there is
   unchanged from full classification.
5. **Each candidate's verification is largely independent of the others**,
   making this naturally parallelizable, including as a genuinely distributed
   effort — many contributors, each taking one candidate to full exhaustion, no
   coordination needed beyond merging proven verdicts into a shared cache.

**Why point 3 is consistent with how this project's own engine already works —
checked directly against the real code, not asserted**: `SolvedPositionDatabase`
in `compositional_trajectory_solver_modular_shapes.cpp` already separates exactly
these two concerns. The full computation state lives in one in-memory structure
(`map<pair<string,char>, SolvedPosition> solved`) — this has to hold everything,
the same way full classification always has to. But what gets *written to disk*
is governed by a separate `pending_export` vector, tracking only keys added or
mutated since the last write — the code's own comment states plainly this is
"what makes frequent, cheap on-disk checkpointing possible," and
`append_new_to_file()` costs "O(size of the delta), not O(total database size)."
The existing architecture already treats "what's needed to compute" and "what's
worth persisting long-term" as two different questions with two different costs.
Point 3 above is a direct extension of that same pattern.

**What this grounding does and does not establish — stated precisely, since it
matters for §11.** The existing engine confirms the *general architectural
approach* is sound: enumerating downstream from a given starting position,
memoizing fully in memory, persisting selectively. It does **not** mean the
current `compositional_trajectory_solver_modular_shapes.cpp` or
`run_full_sweep.py`, unmodified, are the tool this process actually runs on. Both
were built and tuned for a different job — classifying an entire material's full
landscape end to end. The boundary-state workflow needs its own build: forward
search from specific starting positions, early termination against a boundary
catalog (§10), and persistence of only proven verdicts rather than the full
`solved` map. That build reuses real pieces of the existing architecture — the
memoization pattern, the delta-checkpoint approach, the general retrograde-search
machinery — repurposed and extended, not the current scripts run as-is against a
different output target.

---

## 4. Assessment — kept clearly separate from the record in §§1-3

### 4a. The candidate-list filter (§2's "sufficient material + check potential")

Tested against a concrete counterexample: **king-and-pawn opposition and
zugzwang.** One side's move loses, the other's doesn't, in positions with zero
check activity anywhere nearby — opposition is a central, textbook concept in
endgame theory precisely because the deciding factor is king geometry and tempo,
not tactics. A "check potential" filter would discard exactly this class of
positions. The sound version — "no check is *ever* possible under any
continuation from here" — is a true global fact, but establishing it costs the
same as the computation the filter was meant to avoid.

### 4b. What genuinely doesn't work: recomputing unproven subtrees on demand

This is the one part of the original proposal that does fail, and it's narrower
than "the whole caching idea" — specifically, treating full-landscape memoization
as avoidable *during* the search itself, re-deriving a subtree each time it's
revisited rather than keeping it resident. This project already ran the closest
real experiment to this exact tradeoff, in the opposite direction: caching *more*
(`reachable_shapes` by canonical form, a shared cross-pass move-tuple cache, both
in `general_analyzer.py`) to avoid recomputation was tested and came back **~6.5%
slower, twice**, because `explore()` still walks the full subtree regardless of
cache hits — traversal, not value storage, is the dominant cost. Discarding the
in-memory computation early and re-deriving it later hits the identical wall from
the other side: more re-traversal makes the already-dominant cost worse. This is
exactly what §3's point 3 does *not* propose — it keeps the full in-memory
computation intact for the single pass, and only economizes on what gets written
out afterward, which is a different, and sound, question.

### 4c. The batch-size-1 analogy

A small batch covering a large fraction of positions quickly isn't caused by
*what* gets cached — `discover()`'s frontier expands based on what's reachable
from a seed, independent of batch size. The observation actually shows full
exhaustive discovery already happens fast, with full caching, from very few
seeds — evidence that the underlying graph traverses quickly once you're doing
it, not evidence that persisting less of the result would help.

### 4d. What survives, unambiguously real

The insufficient-material ceiling exclusion (condition 1) is real, narrow, and
already implemented in the engine's own per-position terminal checks (confirmed
directly against the engine's own comments on unconditional-draw terminals) —
though not yet as a pre-generation material-category skip; whether that's worth
adding is a small, separate, real question. And once a material is already fully
solved — by this project or by Syzygy — finding its boundary states is genuinely
cheap (§7): a direct scan, no new tablebase construction. That part of the
proposal's spirit is correct.

---

## 5. Relationship to this project's existing architecture

A boundary board state is a specific, checkable pattern over values `classify()`
already computes in full — "exactly one child is non-losing" — the same shape as
`analyze_longest_line.py`'s `NumTiedAlternativesThisPly == 1` check for forced
plies, applied to win/draw/loss instead of distance/escape. No new solving engine
is required for *identifying* boundary states in already-solved data — a new
*scan*. §10 covers the separate, further build this enables for materials not yet
solved.

---

## 6. Where this gets easy for already-solved material

For any material Syzygy already covers (up to 7 pieces), boundary-state
identification needs no tablebase construction — legal move generation plus WDL
probes, both standard, already-solved problems:

```python
import chess
import chess.syzygy

def is_boundary_state(board, tablebase):
    """
    Returns True if the side to move has exactly one non-losing legal move
    and every other legal move loses -- checked directly against
    already-computed Syzygy WDL values, no new computation.

    See section 7 for the collapse-cursed/blessed choice this makes, and
    section 8 for what "non-losing" is deliberately left agnostic about.
    """
    non_losing = []
    for move in board.legal_moves:
        board.push(move)
        wdl = tablebase.probe_wdl(board)
        board.pop()
        # WDL is from the perspective of the side to move AFTER the push,
        # i.e. the opponent -- so a losing value for them (-2 or -1) is a
        # non-losing outcome for the side that just moved.
        if wdl <= -1:      # opponent loses or blessed-loses -> we're fine
            non_losing.append(move)
        elif wdl == 0:
            non_losing.append(move)
        # wdl >= 1 (opponent wins or cursed-wins) -> this move loses for us
    return len(non_losing) == 1
```

**Addendum (proposed revision, not yet run against real tables; syntax-checked
only).** The sketch above counts moves that are non-losing for the mover
without first checking that the position itself is drawn. As written it would
also return True for a *won* position in which the mover has exactly one
winning move and no drawing move — an "only move to win" position, which is a
different object from the drawn positions §1 defines. A guard on the
position's own WDL closes this, and `boundary_report()` returns the
legal-move count alongside, which is the natural refinement §8 already asks
for to separate a genuine choice-point from "only one legal move existed."

```python
def is_boundary_state(board, tablebase):
    """
    True if the position is DRAWN and the side to move has exactly one
    non-losing legal move (every other legal move loses).
    """
    # Guard: the position itself must be a draw. After the cursed/blessed
    # collapse (section 7), any nonzero WDL is decisive.
    if tablebase.probe_wdl(board) != 0:
        return False
    non_losing = []
    for move in board.legal_moves:
        board.push(move)
        wdl = tablebase.probe_wdl(board)   # from the opponent's perspective
        board.pop()
        if wdl <= 0:   # opponent loses/blessed-loses, or draw -> not a loss for us
            non_losing.append(move)
    return len(non_losing) == 1


def boundary_report(board, tablebase):
    """Boundary verdict plus legal-move count (section 8 refinement)."""
    legal = board.legal_moves.count()
    return is_boundary_state(board, tablebase), legal
```

In a drawn position no move can win (that would make the position itself a
win), so the `wdl <= 0` branch there counts exactly the drawing moves. The
`legal == 1` case is what §8 flags as "no genuine choice existed." A further
guard — that both sides are materially live (condition 1) — is not yet in the
sketch; §14 item 2 gives the proposed definition.

The same logic applies to any material this project has already fully solved — a
scan over the completed database, structurally identical to how
`analyze_longest_line.py` already finds forced plies. The only place real, heavy
computation is still required is a material neither Syzygy nor this project has
solved yet, and that's not a new cost this idea introduces — it's the same cost
already being paid for KRPvK right now, for entirely separate reasons. §10 is
where that changes.

---

## 7. The 50-move rule — the full reasoning, then what was verified

**As originally worked through**: the 50-move rule implies a board position might,
geometrically, appear to avoid a forced draw — but the specific sequence of moves
required to *reach* that position might itself force a draw first, under the
50-move rule, before the position is ever actually reached. If so, such a
position is impossible to reach without violating the rule, and a
board-geometry-only analysis wouldn't be the true exhaustive list. An alternative
possibility was also raised: Syzygy might not track the 50-move rule at all, in
which case the whole consideration is irrelevant. Stated resolution at the time:
it doesn't matter which is true, since this project's own engine tracks no move
history at all, and the rule "can be applied post-hoc regardless."

**What direct verification against Syzygy's own documentation found — more
precise than either possibility above**: `probe_wdl` returns five values, not
three — `-2` (unconditional loss), `-1` (**blessed loss**, mate forceable but the
50-move rule saves the mover), `0` (draw), `1` (**cursed win**, mate forceable but
the 50-move rule would turn it into a draw), `2` (unconditional win). Syzygy
*does* structurally track the 50-move rule. But only relative to a hypothetically
fresh clock at the probed position — the documentation states unconditional
win/loss values hold "assuming 50-move counter is zero"; `probe_wdl` never reads
a real game's actual current halfmove clock.

**The reconciling fact**: this project's own engine tracks no move history or
halfmove clock at all — pure distance-to-mate, exactly matching Syzygy's
"assume a fresh clock" semantics. So the correct practice when bootstrapping from
Syzygy is to **collapse cursed-win into unconditional win, and blessed-loss into
unconditional loss** — not a workaround, but the direct, correct alignment
between two systems that already agree on the underlying question. §6's code
already implements this collapse via its `<= -1` / implicit `>= 1` thresholds.
Neither original possibility was quite right; the truth was more specific than
either.

---

## 8. Draw-mechanism taxonomy

Exactly four position-level drawing mechanisms exist under standard chess rules:
**stalemate**, **insufficient material**, **repetition** (threefold, claimable,
or fivefold, automatic), and the **50/75-move rule**. Repetition and the 50-move
rule are treated as interchangeable for this project's purposes, by deliberate
choice (§1) — both mean "the position is drawn," and which specific rule a real
game would invoke to claim it doesn't change whether the position qualifies.
Insufficient material is already the real, narrow ceiling exclusion (§4d).
Stalemate is the one mechanism needing real care: stalemate as the *realization*
of a boundary state's drawing branch (the non-losing move leaves the *opponent*
stalemated on their reply) is a real, legitimate defensive resource with genuine
precedent (well-known "stalemate tricks" in king-and-pawn endings), not something
to exclude by default. What should be flagged, not silently conflated with
genuine boundary states: a "boundary state" whose only distinguishing feature is
that no other legal move existed at all, rather than a genuine choice among real
alternatives where all-but-one specifically lose. `is_boundary_state()` above
doesn't yet distinguish these — checking `board.legal_moves` count independent of
the non-losing count is the natural refinement.

**Insufficient material itself splits into two genuinely different cases, and
conflating them is a real mistake worth naming precisely, since it was
actually made and corrected in this project's own record (§2, Method 4).**
*Unconditional* insufficiency — KB vs K, KN vs K — is drawn everywhere,
position-independent, with no possible boundary-adjacent structure of any
kind; nothing about board configuration ever changes the result. *Position-
dependent* fortress material looks similar at first glance but isn't the same
category at all: wrong-colored-bishop-with-rook-pawn material is won in many
positions and drawn in others, decided by a genuine race between the two
kings for a specific corner — exactly the narrow-margin, tempo-dependent
structure this whole document is built around, not a case falling outside
it. The test that actually separates the two: does *every* position in the
material share the same result regardless of configuration (unconditional),
or does the result depend on the specific arrangement (position-dependent,
and therefore a real source of boundary states, not an excluded category).

**Wrong-colored bishop and opposite-colored bishops are different concepts, and
only the first is the material this document's construction work is built on.**
Wrong-colored means one bishop that is the wrong color for a rook pawn's
promotion square. Opposite-colored means each side has a bishop, on opposite
colors. The second is the far better-known drawish structure in opening and
middlegame theory, but its draws are robust blockades, and the endgame
literature states that zugzwang is less common there than with same-colored
bishops. It should therefore be expected to produce wide fortresses, not
exactly-one positions; the transition into the blockade is the more plausible
place for a tempo-critical moment (`boundary_board_state_topology.md` §8d).

---

## 9. Relationship to existing chess-computing literature

**This sits inside a real, decades-old research tradition, and it already covers
this project's own material.** Althöfer and Walter, *"Weak Zugzwang: Statistics
on some Chess Endgames"*, ICCA Journal 17(2), 1994 — an exhaustive,
tablebase-based statistical study across four endgames, one of which is **KBNK**,
the exact material this project has real, verified data for. Guy Haworth's
*"Mutual Zugzwangs in Chess"* (Computer Olympiad Workshop 6, 2001) goes further:
a comprehensive, exhaustive enumeration of mutual zugzwang across essentially
every 2-to-5-man endgame, explicitly listing which materials have zero such
positions — **KBNK is on that list, by name, confirmed exhaustively: zero mutual
zugzwangs.** Both predate Syzygy (2013) by 19 and 12 years respectively — the
real research window for this territory is closer to 45-50 years, not the 13
implied by dating it from Syzygy alone.

**The precise distinction**: mutual zugzwang and weak zugzwang both require
toggling which side is to move on a fixed arrangement and comparing the result —
a two-position comparison. Traditional zugzwang (the older concept both papers
build on) asks whether the value *category* changes under that same toggle —
closest in spirit to this document's own concept, but still a comparison across
a toggled pair. A **boundary board state (§1)** requires no toggle at all: fix
the side to move as given, and count how many of *that side's own* legal moves
stay non-losing. Related to, but not identical to, anything found published under
the zugzwang name specifically.

**Honest verdict**: not operating in a vacuum, and not simply "zugzwang, already
done" either — a genuinely adjacent question to a well-established one, built on
the same real infrastructure that already proved workable at this exact
material's scale. The overlap between which positions get flagged by either lens
is presumably large but unmeasured — a real, open, directly checkable question.

**The concrete action this adds to §11's plan**: Haworth's KBNK-zero-mZZ result is
a precise, already-published, external number to check this project's own
boundary-state scan against, the moment it runs on real KBNvK data. Since a
boundary board state doesn't require the mutual property mZZ does, boundary-state
count should come back *larger* than zero even where mutual-zugzwang count is
exactly zero — a real, falsifiable prediction, not just a plausible-sounding one.

**A second precedent, for Method 3 specifically, with the distinction stated
precisely rather than left loose.** Možina, Sadikov, and Bratko (University of
Ljubljana, roughly 2005-2011) built a sustained research program mining the
KBNK tablebase itself for human-understandable strategic knowledge — using
argument-based machine learning combined with specialized minimax search to
induce a goal hierarchy (restraining the defending king, driving it toward the
correct corner), validated against real chess teachers and grandmasters. This
is real, genuine precedent for "a tablebase contains more usable structure
than its raw values alone" as a general approach. It is not a precedent for
Method 3's specific technique: their target is the *entire* tablebase, their
method is rule induction toward human-readable goals; Method 3's target is the
*boundary-state population specifically* — a category that doesn't exist in
their framework — and its method is structural feature correlation toward
candidate prioritization, not strategic explanation. Family resemblance in the
general instinct, not a mapping of the specific method — worth being exact
about that distinction rather than letting the precedent's real strength imply
more than it does.

**A third precedent, for the scaffolding strategy in §10 and §11 step 10 —
recorded in full in `boundary_board_state_topology.md` §7 and pointed to here so
this section stands complete on its own.** Checkers was weakly solved by
combining an exhaustively built retrograde endgame database with forward
proof-number search from the starting position, terminating wherever the search
reconnects to the database — described in that literature as raising the water
level of where lines terminate. This is a tight structural match to the
scaffolding strategy, unlike the looser family resemblances above. The same
literature treats applying the technique to chess at full scale as blocked by
the scale gap, which is why the narrower, topology-mapping version recorded here
does not attempt whole-game solving.

---

## 10. The scoping mechanism — where the real impact actually comes from

Everything in §§1-9 establishes *what* a boundary board state is and how to find
one in already-solved material. This section is the separate, further claim: what
a **validated catalog of them** lets you do that full classification alone
doesn't hand you for free, and why that's the actual path to impact beyond
cataloging materials that were already fully solved anyway.

**The logical justification, made explicit rather than left implicit.** A
position whose own perfect-play line routes into a known, validated boundary
state is *itself* provably a draw — not by proximity, not by heuristic, but
directly from condition 3 (§1): if a forcing win existed anywhere upstream along
that line, perfect play would have taken it instead, and the boundary state
simply wouldn't be the outcome reached at all. A validated boundary state's very
existence already certifies "nothing better was available anywhere feeding into
this" — the boundary wouldn't be a boundary otherwise. This isn't a new
assumption; it's condition 3 applied to a specific, useful consequence that
wasn't spelled out in earlier drafts.

**The mechanism this actually enables**: forward search from one specific
starting position of interest, terminating the instant it reaches an already-
known boundary state, rather than needing the rest of that branch explored or
the whole material's landscape pre-built. This is the concrete form of "the
architecture for enumerating from a single position downstream is sound" — the
existing engine's general approach (retrograde discovery seeded from starting
positions) is the right *kind* of tool, repurposed (§3) so that a validated
boundary catalog acts as an early-termination oracle during that search, the same
role a Syzygy WDL probe already plays when a real chess engine's own search hits
tablebase-covered material mid-game.

**What this does and does not achieve — stated precisely, because the
distinction is the whole point.** It does **not** make building the complete,
every-single-position full landscape of a material cheaper — full backward
induction already reuses every child's computed value automatically the moment
it's known; a separate catalog adds no further savings on top of that for the
"classify everything" goal specifically. What it *does* achieve is different and
genuinely new: determining the outcome of a **specific position, or a scoped set
of them**, without needing that material's entire landscape classified first.
That reframes the actual target — not "classify every KRPvK position," but
"determine this position's outcome, and anything else that shares its boundary."

**The connection to `--preload-from` — the same reduction logic, one level more
abstract.** `general_workflow_reference.md` already documents and verifies a
material sealing instantly against an already-solved *simpler material's full
database* the moment it reduces into one via capture or promotion. What this
section describes is the same principle applied to a **boundary-state catalog**
instead of a full sub-material table: a large material's forward search, hitting
a position whose relevant sub-structure is already a validated boundary state
(of that same material, or one it reduces into), seals immediately the same way.

**The reframed central open question, exactly as intended — not "which ceiling
do we hit first," but "how far does this actually reach."** Whether targeted
forward search with boundary-catalog early termination is tractable for
materials whose *full* exhaustive classification is not, is the real, currently
untested question this whole document has been building toward. That's the
literal path past Syzygy's own 7-piece ceiling and past whatever this project's
own engine can fully classify end to end: not a bigger version of full
classification, but a fundamentally more targeted question, asked and answered
one position (or one scoped candidate set) at a time, using a catalog built the
way §§2-3 describe. Finishing the boundary board state problem for a material
already within reach is the concrete, buildable step that stands between here
and being able to ask that question on real data instead of in the abstract.

---

## 11. Concrete next build

1. **Bootstrap from Syzygy** for any material already covered, using
   `is_boundary_state()` (§6, with the §7 collapse applied) — a pure scan, no new
   tablebase construction, real boundary states out immediately.
2. **Seed Method 2's retrograde walk (§2) from Syzygy's own decisive positions**,
   directly alongside step 1 — walk backward from known wins/losses, validating
   each traced step against the ancestor's full set of children (via Syzygy or
   this project's own data, discarding branches the ancestor's own optimal play
   wouldn't take) until reaching the first ancestor with any non-losing option,
   then verify that transition point the same way step 1's candidates are
   verified. No new solving technique required, only a new way of directing the
   existing retrograde machinery.
3. **Extend both methods to this project's own solved materials with genuine
   drawn territory** — KRPvK once complete, and the Method 4 ablation
   sequence (KQvNK, KQvBK, KQvBNK) once built — the same way
   `analyze_longest_line.py` already scans for forced plies. **Not** KQvK or
   KBNvK directly: both are structurally empty of boundary states (§2, Method
   4) — running this step against either is expected to return zero, and that
   expected-empty result is itself the first real test of that claim, not a
   wasted run.
4. **Confirm the expected-zero result for KBNvK against Haworth's
   already-published zero-mutual-zugzwang result** (§9) — not a comparison of
   two independent numbers hoping to match, but two different arguments (this
   document's structural one, Haworth's exhaustive enumeration) predicting
   the same zero from two independent directions. A non-zero result here
   would mean something in this document's own reasoning is wrong, not just
   surprising.
5. **Build the repurposed forward-search tool** (§10) — a new, dedicated
   pipeline, not the current `compositional_trajectory_solver_modular_shapes.cpp`
   or `run_full_sweep.py` run unmodified against a different output target.
   Reuses their memoization pattern and delta-checkpoint approach (§3), extended
   with boundary-catalog early termination during forward search from a
   specified starting position.
6. **Cache proven boundary-state verdicts to a small, dedicated file**, following
   §3's point 3 — persisting only the proven-verdict subset, not the whole
   landscape, using the same compute-fully/persist-selectively pattern the
   engine's own `pending_export` mechanism already demonstrates is sound.
7. **Test the scoping mechanism (§10) on a material this project can already
   fully classify**, comparing forward-search-with-early-termination against
   full classification directly — the first real, measured data point on whether
   the approach actually reaches further than full classification does, before
   trusting it on a material too large to check the answer against.
8. **For a material neither Syzygy nor this project can fully classify**: §10's
   mechanism, not full classification (§4b), is the candidate path forward — with
   step 7's own result as the evidence for whether it's real.
9. **Once steps 1-8 are working at a small, verified scale, the candidate-list
   resolution problem itself is the natural target for the distributed model in
   §12** — a coordinating registry of candidates, independent contributors each
   resolving one, permanent caching of verified verdicts (§3), following the
   GIMPS precedent rather than attempting a single, centrally-coordinated sweep.
10. **Extend confirmed boundary states outward via the scaffolding strategy in
    `boundary_board_state_topology.md` §7** — tracing further upstream than
    strictly needed for each state's own confirmation, recording which
    confirmed states turn out to connect to each other and which specific
    larger materials mark the current frontier, using that frontier as the
    prioritized target for steps 8-9 rather than an arbitrary choice of what
    to solve next.
11. **Build the primitives catalog** (`boundary_board_state_topology.md` §8b):
    for each entry, record the deciding inequality and the locked-squares
    signature, taking each small material's boundary states from its exact
    tablebase, not from a heuristic formula.
12. **Attempt one composition (Method 5) on paper**: choose a single concrete
    placement — the 4-against-4 diagonal chain, or something simpler — write out
    the tempo timings for every route that breaks the frozen structure, decide
    wide fortress versus exactly-one, then run the small obstacle-board solver
    as insurance (topology §8d).
13. **Build a population of confirmed local boundary states before attempting
    decomposition** (topology §8m, correcting this step's own earlier draft,
    which proposed starting from rook-pawn KPvK alone). Scan every material
    Syzygy or this project already covers for *local* boundary states — drawn,
    exactly-one, failing the materially-live test on one side (topology
    §5a-iii) — across as many distinct materials and entry families (topology
    §8b) as the scan reaches, not only the rook-pawn corner draw. Depends on
    the definitional decisions in §14.
14. **Discover primitives from that population and generate candidates by
    subtraction** (topology §8l, §8m). Decompose each population member by
    minimization-by-ablation into fragments with explicit interfaces; catalog
    them; look for what recurs across the population; remove one primitive
    from an entry and solve for a structure that discharges the residual
    assumptions. This requires the obstacle-board solver (walls as blocked
    squares) as required infrastructure, not insurance — Syzygy cannot
    evaluate positions with walls. If step 13's scan turns up too small a
    population for recurrence to be meaningful, fall back to the single-example
    path on rook-pawn KPvK alone (topology §8j, steps 5 onward).

---

## 12. Open, distributed resolution — the practical path forward, complementary to Syzygy

Everything in §§1-11 specifies *what* a boundary state is, *how* to generate and
validate candidates, and *why* the resulting catalog has real leverage (§10).
This section is about a separate, further question: given a candidate list that
may be large and a per-candidate verification cost that stays real regardless of
which generation method produced it (§4b), what actually gets it resolved at
scale — and why that's a genuinely different shape of problem than the one
Syzygy itself is bounded by.

**Where Syzygy's own ceiling actually sits, stated precisely.** Syzygy's limit
isn't a limit on *knowing how* to solve a material — it's a limit on *storing*
the result of solving it exhaustively. A complete table requires a value for
every legal position of that material, and the storage this demands grows
enormously with piece count; 7-piece coverage already required dedicated,
large-scale computing resources, and 8-piece coverage runs into the same wall,
harder. This is a constraint on full, simultaneous, every-position storage — not
on whether any *individual* position's value can be determined.

**Why the boundary-state catalog is a different shape of problem, not subject to
the same wall.** It doesn't require every position's value stored at once — it
requires a (hopefully much smaller) candidate list, each entry resolved by an
independent verification (§2), with the result cached permanently once found
(§3) and never needing re-derivation. That's a fundamentally different
computational shape: not "store everything about one material, all at once," but
"resolve one candidate at a time, permanently, and accumulate." Work on it
incrementally, over however long it takes — months, years, decades — without
needing the storage capacity full landscape construction would demand.

**The real, established precedent for this shape of problem, not a hypothetical
analogy.** The Great Internet Mersenne Prime Search (GIMPS) has run exactly this
model for decades: a coordinating registry supplies candidates, individual
volunteers each verify one independently (the Lucas-Lehmer test, entirely
self-contained), and a confirmed result is permanent and never re-checked. This
is worth distinguishing from a different, also-real precedent: the Lomonosov
7-piece tables were built on distributed *supercomputer* resources, but as one
coordinated, simultaneous sweep toward a single exhaustive table — closer in
shape to full-landscape construction than to candidate-by-candidate resolution.
GIMPS is the closer match to what's being proposed here, not Lomonosov.

**One honest place the analogy doesn't transfer cleanly, worth stating rather
than glossing over.** A Lucas-Lehmer test needs nothing but the candidate number
and a CPU — fully self-contained. Boundary-state candidate validation (§2,
Method 2's refinement) is not quite that clean: confirming a traced move was an
ancestor's own best choice may require checking against Syzygy's own tables for
the relevant piece count, or against this project's own growing boundary
catalog — real, if largely static, reference data a contributor needs local
access to or a way to query. A real distribution model for this needs to account
for getting that reference data to contributors, not only the compute task
itself. That's a genuine engineering question, not a reason the approach
doesn't work — GIMPS-style projects solve analogous distribution problems for
their own precomputed data routinely — but it's not free, and shouldn't be
assumed away.

**The honest scale of the claim, consistent with the calibration already
established elsewhere in this document.** Real progress on this candidate list —
for a material beyond what Syzygy or this project can currently fully classify —
would be genuine progress on a problem that specifically struggles with
exhaustive, simultaneous storage, not a small thing. But it is progress on
*that* problem, specifically — extending how far a boundary-state catalog and
the scoping mechanism (§10) reach — not, by itself, insight into chess's own
starting position's status, for the same scale-gap reasons already established
there. Keeping those two claims separate is what keeps this section honest
rather than inflated.

**Left explicitly open, not claimed as solved**: the validated-retrograde-walk
approach specified in §2 is the method documented here, not asserted as the
final or most efficient one possible. Whether other techniques — leveraging the
topology of an already-partially-built catalog, or patterns emerging once enough
candidates have been resolved — could further reduce the per-candidate cost
established in §4b is a genuine, acknowledged, unresolved question, one that
sustained distributed effort over time is arguably better positioned to surface
than any single, fixed method decided in advance.

---

## 13. What remains open

**The separately-reasoned scope assessment referenced in the preamble**: for any
material that is *already fully classified* — Syzygy-covered, or already solved
by this project's own engine — a boundary-state catalog *by itself* is a derived
filter over data already completely known. §10 is this document's own answer to
that limitation — the catalog's real leverage is on materials not yet fully
classified, via targeted forward search, not on materials already finished. This
does not by itself constitute progress toward solving chess from its starting
position — a categorically larger problem (roughly 10^43 to 10^47 legal
positions) that nothing built here touches. A real, argued position, not a
dismissal — it does not change whether building the catalog is worthwhile, which
§1's preamble records as the project's own stated priority regardless.

- **The central open question (§10, stated once, precisely, not repeated as a
  vague "scale" question)**: whether targeted forward search with boundary-
  catalog early termination is tractable for materials whose full exhaustive
  classification is not. Everything else in this list is secondary to this one.
- **Whether Method 2's transition points (§2) are typically exactly-one
  (true boundary states) or typically several (edge-of-region, not boundary
  states under the strict definition)** — a real, empirical, currently untested
  question about the actual shape of the drawn/decisive boundary, not assumed
  either way going in.
- **Whether the conditional closure argument's "wall" (§2) ever falls within
  currently-computable material, or always requires crossing materials beyond
  present reach** — the unconditional, within-material termination guarantee
  holds regardless; this is specifically about whether the broader, start-
  position-tied version ever becomes practically checkable, not just logically
  real.
- **How reference data (Syzygy tables, this project's own growing boundary
  catalog) actually gets distributed to independent contributors under §12's
  model** — flagged there as a real engineering question distinct from the
  computation itself, not yet designed.
- **Whether a method beyond the validated-retrograde-walk of §2 could further
  reduce per-candidate verification cost** — explicitly left open in §12 rather
  than assumed settled by the one method this document currently specifies.
- Whether `board.legal_moves` count alone soundly separates "genuine choice-point
  boundary state" from "only one legal move existed, and it happens to draw"
  (§8) — proposed, not yet tested against real data.
- Whether a boundary-state catalog, once built for a real material, reveals
  further structure of its own (families, clustering, anything analogous to the
  shape graph) — genuinely unexplored.
- Whether the insufficient-material ceiling check (§4d) is worth adding as a
  pre-generation material-category skip, rather than only the per-position
  terminal check it already is.
- **Whether Method 5 can produce any placement satisfying all of its conditions
  at once, and whether that placement is an exactly-one boundary state or only
  a proven wide fortress** — the discovery itself; the full list of conditions
  and the related open questions are in `boundary_board_state_topology.md` §9.
- **The items in §14 below**, which are definitional and need a decision before
  construction targets can be chosen.

---

## 14. Discrepancies and definitional questions found during onboarding

**Status of this section: assessment, not record.** Written when the
documents were re-read from a fresh start (September 28, 2026) and set against
the ideas developed in that session. Nothing here changes §§1-13; these are
points where the record is internally inconsistent or underspecified, kept
visible so that they are decided deliberately and not by drift. Every claim
marked "untested" was reasoned, not run.

**1. The definition of a boundary state exists in two versions that diverge.**
The operational definition is purely local: the position is drawn and the side
to move has exactly one non-losing move (§1 condition 2; the §6 sketch). The
narrative definition adds upstream-flavored parts: condition 3, and the
exclusion of stalemates that arise from sub-optimal play. §1 also says the
concept is "a pure, self-contained property of the position itself," which
sits uneasily with any exclusion justified by what could have been chosen
upstream. The two versions agree on most positions and diverge on
blunder-derived ones (item 2). Three tiers make the choice explicit:

- **Tier 1 — locally exactly-one:** drawn, exactly one non-losing move.
  Checkable from the position alone.
- **Tier 2 — exactly-one and reachable within its material through non-losing
  play:** Tier 1, plus a validated chain from some drawn position of the same
  material in which neither side ever needed to have already erred.
- **Tier 3 — start-reachable:** the definition in the topology document's §5.
  Currently uncheckable at scale.

Each tier contains the next. The methods in §2 and the build in §11 should
name the tier they target. Until decided, "boundary state" in this project
should be read as Tier 1. **Correction added afterwards (see item 2):** every
tier also includes condition 1, that both sides are materially live. Tier 1 as
first written omitted it; read Tier 1 as "drawn, exactly one non-losing move,
both sides materially live."

**2. Two predictions about total-win materials conflict — RESOLVED: the flag
was mistaken; the original reasoning is kept below, and the resolution follows
it.** §2 Method 4 and §11
steps 3-4 state that KQvK, KRvK and KBNvK are "structurally empty of boundary
states," and treat Haworth's zero mutual zugzwangs in KBNK as independent
confirmation of an expected zero. §9 states the opposite prediction: that the
boundary-state count should come back *larger* than zero even where the mutual
zugzwang count is zero. They cannot both hold. Reasoning suggests the §9
prediction is the one that survives at Tier 1 (untested): take black Ke5,
white Qe4 undefended (white king far away, say a1), black to move. Kxe4
reaches KvK, a draw; Kd6 and Kf6 are legal and lose to KQvK; the other king
moves are illegal. Exactly one non-losing move among several legal ones, and
the position is drawn — a Tier 1 boundary state in a material Method 4 calls
empty. Analogous positions with a hanging rook, bishop or knight should exist
in KRvK and KBNvK. This falls squarely under condition 1's "reduction from
sacrificing will force a draw." It would be excluded only at Tier 2 or above
(the position is reachable only after the winning side hangs a piece), by the
same reasoning that excludes sub-optimal stalemates. This is checkable in
minutes with a scan of KQvK, and doing so settles the question of which
prediction to keep. Note that Haworth's zero is a statement about mutual
zugzwang, which requires the toggle, and so does not by itself bear on
Tier 1 counts either way.

**Resolution (recorded after discussion with the project's author).** The flag
above was mistaken, and the mistake is kept in the record. Condition 1
requires that *both* sides can still win on piece combination alone. In a
hanging-piece position in KQvK, KRvK or KBNvK one side is a bare king, which
has no winning material, so the position fails condition 1 and is not a
boundary state. The reading that let the example through took the operational
test (drawn, exactly one non-losing move) as the whole definition, and read
condition 1's "reduction from sacrificing will force a draw" as covering the
case, without applying condition 1's first clause. §2 Method 4's claim that
these materials are empty of boundary states therefore stands, as does the §11
step 4 comparison with Haworth. What follows from it:

- The §9 prediction (a nonzero boundary count in KBNK where the mutual-zugzwang
  count is zero) is out of scope wherever a side is not materially live, and
  should be read as restricted to materials where both sides are. It is not
  edited in place; this note is the annotation.
- Condition 1 belongs in the operational check. Proposed definition: a side is
  *materially live* in a given material split if some legal position of that
  split is a win for it under optimal play (a WDL win, after the §7 collapse).
  A bare king, and a king with a single minor piece against a bare king, are
  never live. It is computable once per material split. It is not yet in the
  §6 sketch.
- The example positions are still drawn, exactly-one positions in a material
  where one side cannot win. Under the vocabulary of topology §5a-iii they are
  candidates for *local* boundary states, alongside the rook-pawn KPvK draw.
  Whether they count, given that they involve a reduction to insufficient
  material and no locking mechanism, is a question recorded there.

**3. The §6 sketch lacked a drawn-position guard.** See the §6 addendum. This
one is small and mechanical.

**4. "Local" and "global" now carry four meanings across the two documents.**
Topology §3 (a status is local and unconditional), §5 (global as
start-reachable), §6 (scoped versus global coverage), and the closure-
contribution distinction proposed in topology §5a. The last is the least
settled. It needs either a firm definition or a different pair of labels.

**5. "Insufficient material" is being used loosely for the rook-pawn KPK
draw.** That result depends on the king race and is not unconditional
insufficiency in the §8 sense. The property being pointed at appears to be that
the primitive cannot host additional material without breaking or being
broken by it; "saturated primitive" or "closed primitive" would say that
without overloading §8's taxonomy. Proposed, not adopted.

**6. Higher-order constraints are not captured by pairwise compatibility.**
See topology §8g. This is a design constraint on any future composition
tooling and is recorded here so it is not lost when the build in §11 step 12
is scoped.
