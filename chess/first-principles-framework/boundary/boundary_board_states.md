# Boundary Board States: The Choice-Force Boundary of Perfect Play

**This is recorded as the primary objective of this project, not a secondary
analysis — stated explicitly and preserved as a fact about the project's own
priorities, not something this document should soften.** A separate, honestly
reasoned assessment of this idea's actual scope and limits is kept in §13,
clearly labeled as a distinct, independently-argued position rather than folded
into or diluting the objective stated here.

**Provenance of the programme (the designer's statement).** The programme was constructed from the designer's own informal understanding of
chess, which the designer holds was complete in structure before the programme was started. The formal names used in these documents (mating
net, decoy, opposition and so on) are labels applied afterwards to map that structure onto formal chess understanding. The programme was not
built from formal chess theory, and its agreement with formal theory is read as convergence: formal theory and the programme describe the same
structure in different ways. The exactness of that agreement is derivable from the primitive calculations and from how the primitives compose,
by joint compatibility, into a calculable equilibrium; testing it is the purpose of the rediscovery suite (`boundary_board_state_topology.md` §8r).
The designer has no formal training in game theory, combinatorics, chess, systems design, search optimization or computer science, and built the
architecture as an autodidact from this informal understanding. The designer states this as a fact about the origin of the programme and not as
a supplement to any proof, and regards the framing as important to the programme.

This file is the construction pipeline and the record of the idea: what a boundary state is (§1, §14), how candidates are generated and validated (§2, §3), and how the programme is framed (§15). The other documents of the programme are:

- `boundary_board_state_topology.md`: the abstract layer, covering the status of a confirmed boundary state, start-reachability, the scaffolding strategy and composition from race-condition primitives.
- `boundary_board_state_primitive_catalog.md`: the three primitives (the mating and capture net, chained dependencies, the coverage race) and the composition rules.
- `boundary_board_states_opening_survey.md`: the anchored survey from the initial position, the closure and its certificate.
- `boundary_board_states_start_selection.md`: the choice of where to begin, with the designer's statements on choice regions, relativistic comparison, triangulation and the two senses of reachable.
- `boundary_board_states_syzygy_mining.md`: the Syzygy mining process.
- `boundary_board_states_case_study_otb.md`: a case study from over-the-board play.
- `chess_concepts_survey.md`: a survey of chess concepts read against the programme.

Kept separate from `syzygy_improvement_proven_results.md` and the other Claim A/B/C documents because
this is a genuinely different kind of object: those are about distance and escape
count; this is about win/draw/loss pattern-matching at the single precise point
where a losing side's last real choice exists.

> **Where this stands, and how to resume (read this first).**
>
> - **The concept in one line (§14 item 8).** A boundary state is a drawn position that is an
>   equilibrium: each side has a deviation that loses against perfect play, and the
>   draw is sustained by a cycle of moves of both sides (§14 items 7 and 8). The
>   number of non-losing moves at the node is a recorded width, with exactly one the
>   strict case (§1, condition 2).
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
> - **The live thread, generalized (topology §8f-§8j).** A lock
>   is a *duty*, of which a frozen pawn wall is only the degenerate case; a
>   boundary state is a set of mutually dependent locks whose joint
>   compatibility is a constraint system checked as a whole; and verification
>   is scoped by an interference closure rather than run over the composed
>   position. The immediate next step is to write one primitive (the rook-pawn
>   edge race) as an assume-guarantee contract (topology §8i-§8j).
> - **Definition (settled; §14 items 7 and 8).** Under perfect play the position is a draw,
>   because deviation from the repeated draw line is suboptimal and loses decisively on a
>   primitive race condition. Entailed by the definition: every boundary state has a
>   composition of primitive race conditions that enumerates its equilibrium, each losing
>   deviation of each side being refuted by a primitive race condition; the composition is
>   recorded on the proven entry of the boundary state catalog; a primitive missing from the
>   catalog is reducible from the joint compatibility of the others (topology §8x); and draws
>   by stalemate or insufficient material are not boundary states. What remains is
>   construction: the primitives, the composition rules, and the cost model.
>   A draw by forced stalemate or forced capture into a dead position is a boundary state when it arises with no
>   choice to deviate (item 10), so the exclusion applies to draws that could have been avoided.
> - **The primitive catalog now exists** (`boundary_board_state_primitive_catalog.md`):
>   a specification with an entry template, an index of three primitives (the earlier index of eight is reduced; topology §8ab) of which
>   P-001 (the capture net, whose royal case is the mating net) is written, and the generalization rule by which it
>   covers the fool's mate and any other position. The definition of a global
>   boundary state is restated in §14 item 7 (a cyclic equilibrium, both sides
>   live), and the loop it describes is the hard core and not a primitive (topology §8ab).
> - **Status of the work.** The only code in this document is the
>   `is_boundary_state()` sketch in §6, checked for syntax but not yet run
>   against real Syzygy or project data. No run results are kept in this
>   record. The catalog (`boundary_board_state_primitive_catalog.md`) is a specification, and no
>   proven boundary-state entry exists yet.
> - **Chained dependencies (topology §8z, §8ab, catalog §7).** Chains, which are locks read as the
>   dependency of one element on another, and the movement rule of every unit are specified at the level of primitives. P-003
>   is specified (SPEC, version 3). It derives orbits, route parity, edge dependence and the permanence of coverage from a
>   unit's movement relation without naming a piece, and its race terms cover the activation of a piece or pawn. Nothing has been run.
>   Two predictions are stated in advance and not run (catalog §7.13).
> - **The coverage race (topology §8aa, catalog §8).** P-008 is derived from P-001 and P-003 by the reduction rule: the effect of a
>   perturbation on the held terms of every engaged primitive, evaluated exactly at each position along a line and for both sides
>   from the definitions of those primitives, with the support that held terms share (coverage used once, overload) and the race to
>   open one before it is restored. The coverage is complete over the perturbations whose change meets the support closure that
>   the primitives determine. The dissolution race is in P-008 and not in P-003. SPEC, not run; three predictions are stated in advance
>   (catalog §8.13).
> - **The catalog has three primitives (topology §8ab, catalog §1 and §10.5).** The mating and capture net (P-001), chained dependencies
>   (P-003, which includes the activation of a piece or pawn, the movement relation and the frozen structure) and the coverage race
>   (P-008). Race to a square is a term, edge restriction is an input, the repetition loop is the definition of a boundary state, and
>   activation is the result of a chain resolving; none of these is a primitive. **Guard: do not add a primitive, and do not index one,
>   unless it passes the admission rule of catalog §1** (a conjunction of terms with mechanisms that decides a race under
>   perturbation; irreducible by the reduction rule; not a description of a result or of an emergent state; piece- and
>   board-agnostic; modular, with no overlap). What fails the rule is classified as a definition, input, term, mechanism, derived
>   form or instance, and recorded in the table of catalog §1. Zugzwang without capture is an open residual, not a fourth primitive.
> - **The anchored opening survey (topology §8ac, `boundary_board_states_opening_survey.md`).** The method that uses the catalog along a
>   chosen opening line from the initial position to find a local boundary-state candidate. Local means a verified equilibrium with no
>   claim about the path; global means membership of a family reachable from the initial position under the optimal backtrace constraint (item 9
>   of §14). It has its own cache (forced region, boundary found, none can exist) keyed on the position. Specified, not run.
>   Also recorded there and in §14 items 10 and 11: global boundary states are many, and avoidability is not a ground for non-global;
>   a forced terminal draw with no choice to deviate is a boundary state; a pruned structure that is closed (every move inside the
>   region or a certified loss) is itself the certificate for a local state; and the closure of the global boundary states is the
>   hard core's one failure condition, with the premise that the initial position is a draw.

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

*Note (§14 item 10): a draw realized by a forced stalemate or a forced capture into a dead position is the same as a boundary state when it arises with no choice to deviate, that is, when the equilibrium holds. The exclusion in this text is a filter against positions where the draw could have been avoided, and against cataloging such positions as discoveries. It is not an exclusion of the equilibrium case. Read "cycle" in the definition as the paradigm realization.*

*The cyclic requirement is restated as a definition, and its fit with §8's treatment of
forced stalemate is examined, in §14 item 7. That "exactly one non-losing move" is a
measure of width and not the definition is recorded in §14 item 8.*

**What a boundary board state is not**: it is not a property of *how* a position
was reached, and it does not require the position to have arisen from optimal
play upstream in some real game. It is a pure, self-contained property of the
position itself — the same way this project has always insisted a full landscape
must be correct for every legal position, not just ones a perfect game would
produce.

---

## 2. Candidate generation — five methods, complementary, all feeding the same verification step

*Methods 6 and 7 were added later and follow Method 5.*

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

**The batch-size-1 analogy, offered in support of this**: in
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
distinction, from standard endgame knowledge.** Wrong-colored-bishop-with-
rook-pawn material is not unconditional insufficiency the way KB vs K or KN
vs K is — it's won or drawn depending on a genuine race between the two
kings for a specific corner, decided by tempo and distance, not by material
alone. Standard endgame knowledge: a defending king cut off from the safe corner
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

**Standing requirement, added after a direct challenge to this method's own
scope** (full record in `boundary_board_state_topology.md` §8q). Internal
validity — pieces placed so the composition's race conditions hold — is not
the same as start-reachability: whether that same arrangement is actually
reachable via legal play from the game's own starting position. An internally
valid but unreachable composition is still a legitimate catalog entry (it
shows the structure *can* exist as a boundary state), but it is useless for
this catalog's stated forward-search-scoping purpose (§7a, §10, §8p). Any
composition work aimed at that purpose must be checked against Tier 3 (§14
item 1, "start-reachable"), not left at Tier 1's internal-validity bar.

**Generalization of Method 5** (full record in
`boundary_board_state_topology.md` §8f-§8j). Frozen, zero-legal-move pieces
are one special case of a wider class: a lock is a *duty* — a piece committed
to a post, movable in principle, but any move that abandons the duty loses
under a race condition. Composition then means several such duties held in
mutual dependence, certified by a per-resource accounting and a
non-interference argument rather than by search over the composed position.
The certificate remains a proof; what changes is where the exhaustive work
sits (exact derivation of each primitive from small material, and validation
of the composition rules on small instances). The standing rule is unchanged.

### Method 6: perturbation search — varying a parameter around a known position and watching the value cross a boundary

**Added after a worked case study** (full record, including the case study
itself, in `boundary_board_state_topology.md` §8o). Take a position of known
or engine-estimated value, of any size — not restricted to tablebase-reachable
material — and vary a single parameter: whose move it is, a tempo removed or
added, a defensive resource present or absent. Check the resulting position's
value. Where the value crosses from decisive to drawn, the position sits near
a boundary and is a candidate worth closer examination; the same scan can also
land on a decisive-to-decisive crossing, which is a different, equally real
structure and not itself a boundary-state candidate (see §8o for a worked
instance of exactly this outcome).

**Generates candidates by local search around a reference position, not by
retrograde walk or construction.** This is Method 2's retrograde-seeding
instinct — follow value changes to find the seam — applied as a single-step,
parametric probe rather than a full backward walk, which is what makes it
usable at piece counts retrograde generation and Syzygy cannot reach. Its
confirmation is whatever confirms the reference position's own value: engine
search for anything above tablebase range, or the tablebase itself within it.
That confirmation is the standing rule's heuristic tier (§1), not a
certificate — Method 6 orders and proposes candidates, including the two
different kinds of crossing it can land on, and certifies neither.

**What it costs to run, and what it buys.** Running it is cheap: one parameter
change, one engine evaluation, repeated across a batch of candidate
perturbations. What it buys, when it lands on a drawn crossing, is a candidate
worth taking to full verification (Method 1-5's own certification, or the
tablebase where the material is small enough); when it lands on a
decisive-to-decisive crossing instead, as in the worked case in §8o, it buys a
different kind of information — see that section for what, and for the limits
on reading more into it than the data supports.

### Method 7: primitive-state perturbation — the move-level change in the state of the catalog's primitives

**Added after the activation discussion** (record and assessment in
`boundary_board_state_topology.md` §8y; the primitives are in
`boundary_board_state_primitive_catalog.md`). For a position P, compute its primitive state S(P): for
each side, the status of every catalog primitive that engages (dormant, engaged with held terms and
tempo, or completed). A move m perturbs it: ΔS(P, m) = S(P·m) − S(P). The search runs over primitive
state and not over piece placement. Moves are ordered by their effect on the asymmetry between the
two sides' states, and positions in which each side's invariants hold with every deviation losing are
the candidates for boundary-state composition.

**Gap detection.** Where the proven verdict of a move (tablebase or exhaustive search) changes in a way no primitive's ΔS explains, the unexplained change is a residual interface, and by the reduction rule (topology §8x) it specifies the primitive that is missing or not jointly compatible. An engine verdict can point to where to look; the gap is established by the proof.

**Status.** The inference is deductive under the definition and is not an approximation. Take a position
whose verdict is proven, by tablebase or exhaustive search. If every catalog primitive is equalized at it
and it is nevertheless not an equilibrium, then by the definition some primitive race condition decides
it, and it is one of three: a primitive not yet in the catalog, whose joint compatibility with the known
ones carries the imbalance; a known primitive that is not jointly compatible with another in this
position and is decisive; or a catalog entry that was evaluated wrongly, which is an error in the
entry. All three are developments of the protective belt, the primitives, and none touches the hard core,
the definition of a boundary state. The residual interface fixes what the missing or incompatible
primitive must discharge. The distance from a state to equilibrium is the order in which this reduction is
carried out, and it never changes what is true. A boundary state itself is proven by its composition as a
complete case list, with exhaustive insurance on the first entries. It has not been run, because it needs
the cost model of P-001 and the written entries of the catalog (P-003 and P-008 are specified and not run).

**Perturbation classes (catalog §7.7, topology §8z).** The enumeration over the finite geometry that this method relies on is made specific by the movement relation of each unit. A move perturbs the structure in one of five ways: a unit move (its node in the navigation graph changes, and under a grading the parity of its remaining route flips); a blocker move or capture (edges of line-travelling relations whose dependency sets contain the square appear or disappear); a change of coverage or a pin (edge exposure changes); a relation change (promotion); and the capture of a node of the dependency graph (a cascade along its supporting and enabling edges). ΔS is the sum of their effects. Not run.

**Two-sided evaluation (catalog §8.8, topology §8aa).** S(P) and ΔS are computed for both sides with the colors swapped, and the candidate moves of each side perturb the other's evaluation: the disturbances that can reach a structure are generated from the other side's candidate moves inside the closure of the structure's key set, and a disturbance is available only if the disturbing side's own structure still holds afterwards. This keeps the search over the primitive state scoped and adversarial. Not run.

---

## 3. The caching-architecture proposal

The proposal is two separable claims with different soundness. It is the following
five-point workflow:

1. **Use §1's three conditions as necessary, not sufficient, filters** to build a
   finite candidate pool for a given material — sufficient material on both
   sides, plus whatever structural setup makes condition 2 possible at all. This
   pool is guaranteed to contain every true boundary state, plus false positives
   that fail full verification.
2. **Fully, exhaustively verify each candidate individually** — no shortcut
   inside this step. This is where the real cost lives.
3. **Cache each candidate's proven verdict — boundary or not — permanently, once
   resolved.** This is the crux, and it is a *different* claim from "avoid
   caching the full landscape at all." The distinction, precisely: the proposal
   is not "don't build the full landscape's
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
sketch; §14 item 2 gives the proposed definition. The sketch tests only the strict
case, width one; §14 item 8 records that exactly-one is a width measure and not the
definition.

The same logic applies to any material this project has already fully solved — a
scan over the completed database, structurally identical to how
`analyze_longest_line.py` already finds forced plies. The only place real, heavy
computation is still required is a material neither Syzygy nor this project has
solved yet, and that's not a new cost this idea introduces — it's the same cost
already being paid for KRPvK right now, for entirely separate reasons. §10 is
where that changes.

---

## 7. The 50-move rule — not followed

**The designer's statement.** The program does not follow the 50-move rule. The designer has removed it from the program's chess, which
makes the program simpler. Repetition draws are not affected: the hard core is a drawn equilibrium sustained by a repetition cycle.

**This file's reading, offered for correction.** The 75-move rule is the automatic form of the same rule and is read as removed with it.
The earlier concern, that the rule could make a position unreachable because the moves that reach it force a draw first, no longer
arises.

**What it means for Syzygy values.** `probe_wdl` returns five values: `-2` (unconditional loss), `-1` (**blessed loss**, mate forceable
but the 50-move rule saves the mover), `0` (draw), `1` (**cursed win**, mate forceable but the 50-move rule would turn it into a draw)
and `2` (unconditional win). The unconditional values hold assuming a fresh 50-move counter. Since the program does not follow the rule,
a cursed win is a win and a blessed loss is a loss, so the values are collapsed: `-2` and `-1` are losses, `1` and `2` are wins, and
`0` is the draw. This is the program's own value and not an approximation. §6's code implements the collapse through its `<= -1` and
`>= 1` thresholds, and the engine of this project tracks no move history or halfmove clock.

---

## 8. Draw-mechanism taxonomy

Exactly four position-level drawing mechanisms exist under standard chess rules:
**stalemate**, **insufficient material**, **repetition** (threefold, claimable,
or fivefold, automatic), and the **50/75-move rule**. The program follows the
first three and does not follow the 50/75-move rule (§7).
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
assumption; it's condition 3 applied to a specific, useful consequence.

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

**A generalization of this scoping mechanism, recorded in
`boundary_board_state_topology.md` §8p.** The termination-oracle role
described above uses a single catalog entry as a leaf node: forward search
seals the instant it literally reaches one. §8p records a broader proposal —
using the whole catalog as a network of reference points against which other
positions are classified by their relativistic relationship to several known
entries at once, not only by directly reaching one. This already exists in
exact, worked-out form for a narrow class of king-and-pawn structures, as
classical corresponding-squares theory, which §8p checks and cites directly.
It also carries one sharp guard rail, which §8p records precisely: failing to
find a path from an unknown position to any cataloged entry establishes
nothing on its own, since that failure could reflect the catalog's own
incompleteness rather than a genuine absence of a drawing structure.

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
    decomposition** (topology §8m; not rook-pawn KPvK alone). Scan every material
    Syzygy or this project already covers for *local* boundary states — drawn,
    a losing deviation for the side to move (width recorded; §14 item 8),
    failing the materially-live test on one side (topology
    §5a-iii) — across as many distinct materials and entry families (topology
    §8b) as the scan reaches, not only the rook-pawn corner draw. Uses the definition in §14
    items 7 and 8. A population of boundary states needs a material in which both sides can lose; king and pawn against king has a side that cannot lose, and its exactly-one positions include captures of a hanging pawn into dead draws, which are not boundary states.
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
15. **Add a forward-driven selection layer on top of steps 13-14** (topology
    §8r). Steps 13-14 are retrograde-first because exact anchors are cheapest
    there. Which materials and primitives to catalog first should be chosen by
    start-position reachability (Tier 3, §14 item 1), with forward-trace and
    cooperative-mate (helpmate-style) enumeration used to *propose* primitives
    and the exact methods used to *certify* them. This changes the order of
    work, not what counts as a certificate. Depends on a definition of forward
    perturbation that does not yet exist beyond Method 6's single-parameter
    probe.
16. **Run the rediscovery validation suite** (topology §8r): give each known
    result (rule of the square, opposition and key squares, corresponding
    squares, Uri's rule, Philidor, Lucena, Vančura, triangulation) a primitive
    decomposition and a pass criterion, and check the decomposition reproduces
    the result exactly. Cheap and external, so it comes before any claim of new
    structure.
17. **Run the mating-net report across a corpus of mating positions and, as a
    stronger test, try to recognize the net before the mate** (topology §8u).
    The four terms and the role of each piece are computed with
    no reference to piece types. Not run; recognition from a pre-mate position is
    not attempted.
18. **Define the primitive state S and the perturbation ΔS over the written
    primitives, and apply Method 7** (topology §8y). The per-ply tempo ledger of
    topology §8w is its first instance. It needs the cost model of P-001 and the
    written entries of the catalog (P-003 and P-008 are specified and not run).
19. **Run the two predictions of topology §8z** (catalog §7.13): that the load-bearing sets found by ablation (§8s) are
    the cut sets of the dependency graph read from the same position, and that orbit exclusion reproduces the known
    color-lock verdicts across the whole class. Nothing in this step has been run.
20. **Run the three predictions of the coverage race** (topology §8aa, catalog §8.13): that the key elements (the singleton
    opening sets of the held terms) coincide with the load-bearing elements found by ablation (§8s); that the evaluation of the
    primitive under perturbation, composed over a line, reproduces the verdict of exhaustive search; and that every
    verdict-changing perturbation in the ablation and deviation data meets the support closure. Count the perturbations that must
    be branched against those legal in the footprint, since that decides how much the hybrid saves. Nothing in this step has
    been run.


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

**The separately-reasoned scope assessment stated at the top of this file**: for any
material that is *already fully classified* — Syzygy-covered, or already solved
by this project's own engine — a boundary-state catalog *by itself* is a derived
filter over data already completely known. §10 is this document's own answer to
that limitation — the catalog's real leverage is on materials not yet fully
classified, via targeted forward search, not on materials already finished. This
does not by itself constitute progress toward solving chess from its starting
position — a categorically larger problem (roughly 10^43 to 10^47 legal
positions) that nothing built here touches. A real, argued position, not a
dismissal — it does not change whether building the catalog is worthwhile, which
the statement at the top of this file records as the project's own stated priority regardless.

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
name the tier they target. The tier structure is kept as a description of how much of a position's history is assumed; the definition of a boundary state itself is settled (item 8). Every
tier also includes condition 1, that both sides are materially live (item 2).

Read "exactly one non-losing move" in the tier descriptions as "a recorded width w, with w = 1 the strict case" (item 8: exactly-one is a measure of width, not the definition).

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
settled. It needs either a firm definition or a different pair of labels. The designer's statement of local and global by the path of arrival is recorded in item 9 as the firm definition.

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

**7. The cyclic requirement, and where the record already says it.**

*As stated.* A boundary state is defined by equilibrium and a forced draw under
optimal play. The draw must end cyclically: a valid boundary state needs a
composition of primitives whose narrow constraints present both sides with moves
that sustain a loop of repetition, each move returning to the same equilibrium
structure. No boundary state ends in stalemate unless that is forced under a
sequence, and no evidence of that has been seen. This follows from the definition
and was believed to be recorded already.

*Assessment — separate from the record above.*

- **It is recorded, in §1 ("On which draws count").** A position validated as a
  draw through forced repetition, not stalemate, is a boundary position, and
  stalemates arising from sub-optimal play are excluded. Two things in the record
  had not been tied to it.
- **The operational test does not encode a cycle, so the cyclic property does not
  follow from the test alone.** The operational test was: drawn, with exactly one
  non-losing move (item 8 records that this is a width measure). In king and pawn against king,
  positions that pass it include those whose saving move captures the pawn into a
  dead draw. These are not boundary states:
  capturing a pawn left en prise is a blunder by the pawn's side and ends in a dead
  draw, not an equilibrium; and the stronger side there cannot lose, so there is no
  equilibrium between two sides (§1 condition 1 again).
- **§8 and §1 appear to disagree about stalemate, and the stated definition
  reconciles them.** §8 treats stalemate as a legitimate *realization* of the
  drawing branch, not something to exclude by default, while §1 reads as excluding
  it. The reconciliation: stalemate from sub-optimal play is excluded (§1);
  stalemate forced under a sequence is allowed in principle (§8); a cycle is the
  paradigm. No non-degenerate forced-stalemate case is
  claimed here; the point is the reconciliation of §1 and §8.
- **Resulting definition, recorded as the project's definition of a global
  boundary state (restated in item 8):** a drawn position that is an equilibrium
  (each side has a losing deviation against perfect play), whose drawing branch is
  realized by a cycle of repetition sustained by moves of both sides, with both
  sides live. Positions whose draw is realized by a
  forced stalemate or by a forced capture into a dead position are terminal draws of
  the local kind (Tier 1, topology §5a-iii), not global.
  *Note (item 10): for the equilibrium case, a draw realized by a forced stalemate or a forced capture
  into a dead position is a boundary state when no side has a choice that avoids it and each side has a losing deviation. It is a terminal
  draw of the local kind only where the draw could have been avoided.*
- **A consequence worth stating.** Both sides must be able to lose, since otherwise
  the stronger side's non-losing moves are all of its moves and there is no narrow
  constraint. Both being able to lose implies both being live, so this definition
  contains §1 condition 1.
- **The open choice recorded earlier no longer arises.** Item 8 removes exactly-one
  as a requirement at any node, so a perpetual check with two interchangeable checks
  is not excluded.
- **This is the definition of a boundary state and not a catalog entry.** The repetition loop is the emergent effect of the
  primitives in equilibrium, so it belongs to the hard core. It is not a catalog entry (topology §8ab, catalog §1).

**8. "Exactly one non-losing move" is a measure of width, not the definition.**

*As stated.* Exactly one is not the criterion. A boundary state has at least one, and
potentially several, tied race conditions, and uniqueness would imply a single
primitive in the composition. A boundary state is a position of equilibrium under
optimal play: a forced draw because any deviation by the side that deviates, against a
perfect opponent, loses. Draws by insufficient material and by stalemate are not
boundary states, because they say nothing about an equilibrium of optimal play; they
show a position forcing a draw from a place where it could have been avoided.
Capturing a pawn the opponent leaves en prise is a blunder by the pawn's side, unless
the race was already lost.

*Note (item 10): the exclusion of insufficient-material and stalemate draws in this text applies where the draw could have been
avoided, and is not an exclusion of a forced terminal draw that arises with no choice to deviate.*

*Assessment — separate from the record above.*

- **Where "exactly one" came from.** It is in the document first uploaded: the
  one-line concept, the strict reading of §1 condition 2, and the `is_boundary_state()`
  sketch with `len(non_losing) == 1`. Condition 2 as originally worded does not require
  it ("no options available that don't immediately put them past a boundary state into
  a forced sequence — or that also become a draw"), and the same document keeps a
  separate, looser category for transitions with several non-losing options ("the edge
  of a wider drawn region"). The strict reading appears in the tiers, the scan and the catalog; the
  wording of condition 2 does not itself require it.
- **Move count and primitive count are different quantities.** The number of
  non-losing moves at a node describes the move tree. The number of tied primitives
  describes the structure of the equilibrium. A composition of many primitives can
  leave one equilibrium move or several, and neither follows from the number of
  primitives. What joint compatibility requires is that the equilibrium moves keep
  every invariant, not that there be exactly one of them. The earlier phrasing in this
  record that the exactly-one move is the only move that leaves every lock in place
  (topology §8g) tied the two together and is not a definition.
- **The rook-pawn draw is the paradigm.** The project treats the rook-pawn draw as
  the paradigm local primitive. A test that selects against the paradigm example is
  not the definition.
- **Definition.** A boundary state is a drawn position that is an
  equilibrium: each side has at least one deviation that loses against perfect play, so
  neither can improve and both can worsen; the draw is realized by a cycle of repetition
  sustained by the moves of both sides; and the equilibrium is accounted for by one or
  more tied race-condition primitives whose invariants the cycle's moves maintain
  (catalog §5, recorded on the proven entry as in topology §8x). Both sides being able to lose implies both are
  live. Positions in which
  only one side has a losing deviation are the local kind (topology §5a-iii).
- **Width.** The number of non-losing moves at a node is recorded as a descriptive
  parameter w ≥ 1. The case w = 1 is the strict case the original text called a
  boundary state; w > 1 is a wider equilibrium. Neither is excluded. Width is cheap to
  compute, so it stays useful for ordering candidates and scoping search, as a measure
  and not as the test.
- **Draws that are not boundary states.** Stalemate and insufficient-material draws,
  including the capture of an unprotected pawn, are terminal draws. They are excluded
  for the reasons in §1 and §8, and item 7 records the reconciliation. The exactly-one test alone does not exclude them, which is part
  of why that test is not the definition.
  *Note (item 10): the exclusion applies to terminal draws that could have been avoided. A forced terminal draw with no choice to
  deviate is a boundary state when the equilibrium holds.*
- **Consequence for the starting position.** The argument that the
  starting position is not a boundary state because it has many non-losing first moves,
  rests on the exactly-one test, which is a width measure and not the definition. The starting position is not ruled out
  by width. Whether it is an equilibrium of the kind defined here depends on things not
  established: that chess is a draw, that each side has a losing deviation, and that the
  draw is sustained by a cycle.
- **Consequences for earlier items.** The tiers of item 1 take exactly-one as Tier 1's
  test and should be read with w as a parameter. The question in item 7 of whether
  exactly-one holds at every node of the loop no longer arises. King and pawn against king, under this definition, holds only local
  candidates, and not a population of boundary states.

**9. "Local" and "global" by the path of arrival (the designer's statement; topology §8ac).** A boundary state's own status is local and
unconditional (topology §3) and is the position's. Its *classification* depends on how it is reached. A **local** boundary state is a
verified equilibrium (the primitives show it and exhaustive verification confirms it), with no requirement on the path and none that
its family reach the initial position. A **global** boundary state is a boundary state that is a member of a reachable family that is reachable from the initial
position under the optimal backtrace constraint (survey §5a).
A global state is sound against all deviation replies, and there are many of them: a deviation may lead to a different global state,
which still closes the phase between choice and force, so avoidability is not a ground for calling a state non-global. A state not shown to be in such a family is local. This is offered as the firm definition that item 4 asked for. The record's start-reachability
(topology §5, Tier 3 of item 1) is read as the reachability part of *global*, which is an open reading recorded in
`boundary_board_states_opening_survey.md` §2. Earlier uses of "local candidate" in this file, for positions that satisfy only the
local conditions of a material, are not changed. Nothing in this item changes the definition of a boundary state (items 7 and 8).

**10. Forced terminal draws are boundary states when there is no choice to deviate (the designer's ruling; topology §8ad).** A forced draw
by stalemate, or by a forced capture into a dead position, is the same as a boundary state by the definition when it arises with no
choice to deviate: the position is drawn, each side has a losing deviation, and no side has a non-losing move that avoids the terminal
draw. The classification of such draws as not boundary states rested on their being sub-optimal relative to other branches, so that
the draw could have been avoided. Where exhaustive coverage of all branches shows that some branches force stalemate or a dead
position, that is not evidence against the classification. The text of §1, item 7 and item 8 that excludes these draws is a filter
with a purpose: the boundary-state catalog is not to be filled with stalemate positions that arise from sub-optimal play, and such
positions are not to be treated as discoveries. The definition of a global boundary state in item 7 is read with the cycle of
repetition as its paradigm realization, and with the forced terminal draw as a second realization of the same equilibrium. The case
in king and pawn against king (item 7, positions in which the saving move captures the pawn into a dead draw) is excluded
for a different reason, which stands: the stronger side cannot lose there, so there is no equilibrium between two sides. This closes
the open realization question recorded in `boundary_board_states_opening_survey.md` §10.

**11. The closure claim (the designer's statement; topology §8ae).** The closure of the global boundary board states must be ensured,
otherwise the programme fails. A viable deviation in the choice phase reaches a boundary state, and a route does not pass from choice
to force without hitting one. If it is not the same one, it is a different global boundary state, by definition. If this were not true,
an opening would force a win from the initial position, and no such line is known to exist; the closure is predicated on its
nonexistence, and enumerating it would be the proof. A route out of the closure, with a certificate, would be a forced opening line.
A closure that holds makes the forward and backward trace complete and a weak solution viable. A missing primitive or a failed belt
prediction is not a failure of the programme. The full statement, with its outcomes and the premise that the initial position is a
draw, is in `boundary_board_states_opening_survey.md` §2a. The segmented progression is recorded there in §7a with only what is stated.

---

## 15. The project as a Lakatosian research programme — an honest accounting

*Scope: this section accounts for the evidence behind the primitives, methods and constructions built so far. It says nothing about the definition of a boundary state, which is settled (§14 items 7 and 8).*

**As stated.** The boundary board state catalog functions as the protective
belt of a Lakatosian research programme, surrounding whatever counts as this
project's hard core. Discovering and cataloging the primitives is the actual
point, since that catalog is what feeds the forward-search component the
whole programme needs to make up for what Syzygy cannot supply. That specific
results keep aligning with established chess theory is expected and correct,
not a weakness. Separately, the framework is being used operationally, in
real play, as a formalized understanding that improves capability.

**Assessment — separate from the record above, and held to Lakatos's own
stricter test rather than the softer one of "did the framework survive."**

The structural mapping holds. A hard core — the definitional and topological
apparatus: a boundary state as a narrow, drawn equilibrium, composable from
joint-compatible primitives under a race-condition calculus — has stood
through every challenge raised against it, revised at its edges but not
abandoned. A protective belt — the specific methods, the primitive catalog
entries, the composition rules, the definitional boundaries around "local"
and "global" — has been repeatedly adjusted in response to real anomalies:
the tier structure and condition 1's liveness requirement, the local/global
distinction rewritten twice (§5a through §5a-iii), the Možina relation
corrected and then refined again, the correction that the opening position
fails the exactly-one test (withdrawn, §14 item 8: that correction rested on the strict
reading), the correction that the puzzle's mutual-mate-flip
is not an equilibrium. Every one of these was a principled response to a
specific, checkable objection, not an evasion of one — which is exactly what
keeps a programme non-degenerating in Lakatos's terms.

Non-degenerating is not the same as progressive, and Lakatos's own criterion
for the stronger label is specific: a belt modification counts as progressive
only when it predicts a genuinely novel fact — not used to construct the
modification — that is then independently corroborated. Measured against
that bar rather than the softer one, here is an honest accounting of this
project's own history so far:

- **Genuinely progressive, by the strict test.** One clean instance: the
  proposed decomposition of Uri's rule (topology §8n) — a race primitive
  composed with a color-lock fact — predicted, before checking, that the
  bishop would not appear in the rule's governing formula. That prediction
  was checked against the independently published rule and held. That is
  real novel-fact corroboration, not restated theory.
- **A real, falsifiable prediction, stated and still untested.** Topology
  §7a's claim that composition specifically, not merely narrower scoping, is
  what lets a floor extend past Syzygy's own storage ceiling. Nothing has yet
  been built that tests it — no composed position beyond ordinary search's
  reach has actually been verified. This is exactly the kind of prediction a
  progressive programme needs to keep generating, and exactly the kind that
  needs to be *run*, not just stated, to count toward the label.
- **Valuable, but not novel-fact generation in the strict sense.** Finding
  that corresponding-squares theory already does what topology §8p proposed,
  and that Bleicher's Freezer tool already does part of what §7a proposed, is
  real and useful — it shows the programme's intuitions are pointed at
  something independently-developed theory also found. But confirming an
  idea already exists elsewhere corroborates plausibility; it is not the
  programme's own prediction being vindicated. Lakatos would read this as a
  sign the hard core is aimed at something real, not as progressive content
  in its own right.
- **The larger share of this project's work, honestly, has been protective-
  belt defense, not novel-fact generation.** Most of what's listed above as
  "non-degenerating" is exactly that: real, principled, and necessary, and
  also not yet evidence of a progressive programme by Lakatos's stricter
  standard. A programme that only ever successfully defends its core against
  anomalies, however honestly, without generating new confirmed predictions,
  is non-degenerating but not yet progressive — and that is the fair, current
  label for most of this project's history, not a criticism of the
  individual corrections, each of which was the right thing to do when it
  came up.
- **What would move the needle, concretely — the same list this document
  already converges on.** Run the population scan (topology §8m) and see
  whether primitives actually recur. Build one composed primitive and check
  it against the tablebase (§11 step 14). Test whether a composition can be
  certified at a piece count Syzygy cannot reach — the direct test of §7a's
  prediction. Each of these, run and checked rather than merely argued, is a
  novel fact in Lakatos's sense if it comes back positive. A negative result
  would be equally valuable, reported with the same honesty every correction
  in this document already has been — and would mean revising the hard core
  itself, not just the belt, which is the harder and more useful thing to
  learn early rather than late.

**On "it should align with known chess theory."** Agreed plainly, and worth
saying why: when a belt modification, built from first principles and
checked afterward, lands on something independently-developed theory already
found, that is a good sign for the hard core's underlying intuitions, not a
reason to think the modification was unnecessary. Uri's rule and
corresponding squares are both instances of this. Both are reasons for more
confidence in the programme's direction, even though neither is itself the
strict novel-fact test Lakatos's stronger label requires.

**On the operational claim.** The OTB game and the chess.com puzzle,
independently verified and recorded in
`boundary_board_states_case_study_otb.md`, are real evidence that the
intuitions motivating this programme's hard core are tracking something
genuine in actual play — confirmed multiple times, including a mate-in-one
coincidence that is as clean a piece of independent corroboration as this
project has produced anywhere. That is a different claim from the hard
core's strongest formal commitment — that boundary states compose by
provable, exhaustively verifiable rules into a usable catalog — and the two
should keep being held separately, the way that document's own §6 already
does. A programme's heuristic producing real practical value well before its
hardest formal claims are settled is not a weakness in Lakatos's account. It
is close to the normal case.

**The substrate and the chain (topology §8z).** The movement-relation substrate and the chain primitive are
developments of the protective belt, and they leave the definition of a boundary state, the hard core, untouched. By the
stricter standard above they are not yet progressive. Nothing in them has been run. They make two predictions that can
fail and have not been run (catalog §7.13): that the load-bearing sets found by ablation are the cut sets of the dependency
graph, and that orbit exclusion reproduces the known color-lock verdicts across the whole class. Positive results would be
novel facts in Lakatos's sense, and negative ones would be reported in the same way as the corrections already on record.

**The coverage race (topology §8aa).** P-008 was derived from P-001 and P-003 by the reduction rule, and it is a
development of the protective belt that leaves the hard core untouched. It is not yet progressive by the stricter standard above.
The derivation of one primitive from two others is made under the definitions and is not a novel fact, and nothing in it has been
run. Three predictions can fail (catalog §8.13): that the key elements coincide with the load-bearing elements found by ablation;
that the evaluation of the primitive under perturbation, composed over a line, reproduces the verdict of exhaustive search, a
disagreement being a residual interface that develops the catalog; and that every verdict-changing perturbation in the ablation
and deviation data meets the support closure. The primitive is evaluated exactly at each position, so the saving is the
perturbations that need no branching, and how large that is depends on the size of the closure, which has not been measured. The
designer keeps open that further primitives may reduce the residual coverage.
