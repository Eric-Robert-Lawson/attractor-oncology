# Boundary Board State Topology: Necessity, Redundancy, and Scoped Coverage

A companion to `boundary_board_states.md`, not a replacement for it. That
document is the practical construction pipeline — what a boundary state is, how
to generate and validate candidates, how to cache and scope from them. This
document is the more abstract layer sitting on top of it: what the *structure*
of an already-growing boundary-state catalog can and can't tell you, once more
than one boundary state is known and their relationships to each other and to
their own ancestors start to matter. Written the same way as everything else in
this project — a real claim tested, a real correction made, kept in the record
rather than smoothed over.

---

## 1. On heuristics — the standing principle this document holds to throughout

Stated explicitly, not left implicit: heuristics are legitimate for *ordering
and scoping* which candidates to work on first. They are never legitimate as a
basis for *certainty* about whether something is a boundary state, a necessary
one, or anything else this document discusses. Every claim below is labeled as
one or the other — provable with the same rigor as the rest of this project's
work, or explicitly a prioritization aid — and the two are never allowed to
blur into each other.

---

## 2. On the Lakatosian framing — a real correction to the previous document's hedge

`boundary_board_states.md` described the periodic-table framing as an analogy.
That undersold it, and it's worth saying so plainly rather than leaving the
softer version on record. Lakatos's structure — a hard core that stays fixed,
surrounded by a protective belt of specific claims that absorb real anomalies
and get revised without the core itself being touched — matches this
project's own process structurally, not just evocatively. The hard core has
been the minimax/perfect-play formalism itself, untouched throughout every
revision. The protective belt is every specific technical claim that's actually
been tested and revised in response to a real counterexample: the
check-potential filter, refined after the zugzwang case; the retrograde
chess-dichotomy, corrected after its circularity was found; the closure
argument, split into its conditional and unconditional forms. Each revision
responded to a real anomaly without abandoning the core. That is the actual
shape of what has happened in this conversation, not a metaphor for it.

**What this framing does not, by itself, extend to**: Lakatos describes the
structure of theoretical revision under empirical pressure. It says nothing
about whether a growing boundary-state catalog makes future *computational
verification* cheaper. That is a separate, further claim, addressed on its own
terms in §6 — real, but bounded, and not something the Lakatosian framing
settles just by being apt for the project's own methodology.

---

## 3. The correction: a boundary state's own status is local and unconditional

This is the technical core of this document, and it's worth stating precisely
before building anything on top of it.

**The claim tested**: could discovering an alternative path from some ancestor
of a confirmed boundary state B — a path that doesn't route through B — mean B
was never really a boundary state, or was only a "local" one rather than a
"global" one?

**Why this doesn't hold, worked through directly**: condition 2 (§1 of
`boundary_board_states.md`) is evaluated entirely from B's own legal moves and
the fixed, well-defined value of each resulting child. Under the minimax
framework this whole project is built on, a position's value is determined
recursively by its own children — never by its ancestors, and never by
whatever else exists elsewhere in the graph. B's own move list doesn't change
based on what some ancestor A does with *its* other options. The value of each
of B's children doesn't change either. So nothing discoverable about A, or
about any alternative route A might have, can make B's own condition-2 status
stop being true. A confirmed boundary state stays a boundary state,
permanently, exactly as §3 of `boundary_board_states.md` already established
for caching purposes — this is the same fact, viewed from the angle of "can it
be un-proven," not just "does it need re-proving."

**What this means concretely**: there is no mechanism, in this framework or any
correct extension of it, for excluding or invalidating a confirmed boundary
state based on what its ancestors turn out to have available elsewhere. If a
supposedly-validated retrograde chain leading to B turns out to be wrong, the
error was in the *validation of that specific chain* (§2's generation-versus-
validation split in `boundary_board_states.md`, already built to catch exactly
this) — not in B itself.

**A necessary refinement — the collective claim is different from the one just
addressed, and it's correct.** The correction above answers "can discovering an
alternative ancestor path invalidate *one specific, already-confirmed* B." It
does not answer a different, genuinely valid claim: whether the *complete set*
of all true transition points for a material, taken together, forms a closure
that every decisive position's validated ancestry must cross into somewhere —
not necessarily through any one pre-chosen member, but into the set as a whole.

That collective version isn't new work — it's the unconditional in-material
guarantee already established (`boundary_board_states.md` §2) applied to every
decisive position in the material simultaneously rather than one chain at a
time. Since each individual chain must terminate at some transition point
before its decisive ancestry runs out (finiteness, no loop possible), and the
complete set of transition points is, by definition, everywhere that
termination can happen, every decisive position's ancestry crosses into that
complete set somewhere. Not a separate argument — the same one, restated for
the whole material at once.

**The real, sound, practical consequence: a completeness check, not an
invalidation mechanism.** If a validated chain from a decisive position is
found that never touches anything in the current catalog before its decisive
ancestry runs out, that's proof — not a guess — that the catalog is missing a
transition point somewhere along that specific chain. Genuine, sound
gap-detection: a way to actively drive candidate discovery, distinct from
anything already in the construction pipeline's own generation methods.

**Where the earlier correction in this section still holds, precisely.**
Finding such a chain falsifies the hypothesis "the catalog is already
complete." It does not, and structurally cannot, falsify any specific,
already-confirmed member of that catalog — B's own condition-2 status remains
exactly as unconditional as established above. "Exclude" is the right word for
the completeness hypothesis. It is the wrong word for any individual confirmed
boundary state. These are two different targets, not a contradiction between
this refinement and the correction it's refining.

**One further precision, tied to an already-open question (§9, and
`boundary_board_states.md` §13).** The unconditional guarantee concerns
transition points generally — the first ancestor where a decisive position's
non-losing-option count moves from zero to nonzero — not specifically strict,
condition-2 boundary states. A routing-around chain proves some transition
point is missing from the catalog. Whether that missing point turns out to be
a true boundary state or the weaker, several-options edge case remains the
same open, empirical question already on record — this argument locates where
to look, not what will be found there.

---

## 4. What's actually being pointed at: necessity versus redundancy, a real and different property

The instinct behind "local versus global" is real — it just isn't a property
of B. It's a property of the *relationship* between B and a specific ancestor.

**The precise question**: for a given ancestor A, is B the *only* validated
path A has to a confirmed draw, or does A have another, independent validated
path — through B, or through some other confirmed boundary state entirely?

- If B is A's only validated route to safety, B is **load-bearing** for A —
  removing B (hypothetically) would leave A with no confirmed way to draw.
- If A has two or more independent validated routes to confirmed draws, B is
  **redundant** for A specifically — A doesn't depend on B in particular, even
  though B remains, on its own terms, a perfectly real boundary state.

**The graph-theoretic anchor, named precisely rather than loosely borrowed**:
the closest established concept is *dominators* from compiler theory and
directed-graph analysis, not the classic articulation point (which is defined
for undirected graphs; a retrograde game tree is directed — edges only run one
way). A node D dominates a node N if every path from the graph's root to N
passes through D. "B is load-bearing for A" is exactly "B dominates A, within
the subgraph of A's own validated safe paths." The spirit matches articulation
points closely enough to be a useful intuition; the mechanics match dominator
analysis more precisely, and that's the concept worth building any actual tool
against.

**This is provable, not heuristic — worth being explicit about, given §1**:
determining whether A has a second independent validated path is exactly as
rigorous as validating the first one (§2 of `boundary_board_states.md`) — full
verification, not approximation. "B is load-bearing for A" and "B is redundant
for A" are both certain facts once the relevant paths are checked, not guesses.
What's heuristic is only the *decision about which ancestors to check this for
first* — a scoping choice, not a claim about the answer.

**What this is actually useful for — prioritization, explicitly not
exclusion**: a load-bearing boundary state has more riding on its own
correctness than a redundant one, since more of its own ancestors' certainty
depends on it alone. That makes "is this candidate load-bearing for a lot of
already-confirmed positions" a legitimate, real signal for which unresolved
candidates are worth verifying first — a genuine use of heuristic ordering,
exactly the kind §1 already sanctions, and nothing more than that. It does not
let any confirmed boundary state be discarded, and no claim in this document
should be read as saying otherwise.

---

## 5. Start-reachability — anchoring "global" to a fixed point, not a relative one

*Clarification (§8ac, §8aq, `boundary_board_states.md` §14 item 9): "global" below is start-reachability. The designer's definition of a global boundary state is membership of a family reachable from the initial position under the optimal backtrace constraint (`boundary_board_states_opening_survey.md` §5a), sound against all deviation replies, with many global states. The reading that start-reachability is the reachability in that definition is recorded in the survey §2 and §5a.*

§4's necessity/redundancy is anchored to *some* ancestor, chosen arbitrarily —
useful for prioritization, but relative. What's being described here is
sharper: anchor "global" to one fixed, unique reference point — the game's
actual initial position — rather than to whichever ancestor happens to be
under consideration. This is a real improvement in precision over §4's framing
for the specific question of which boundary states matter for understanding
chess itself, not just for ordering verification work.

**The precise definition**: a boundary state is *start-reachable* if there
exists a validated chain connecting it back to the true, 32-piece starting
position, where at every step along that chain, both sides always had at least
one non-losing option available and the chain never required either side to
have already made a losing move to arrive there. A boundary state that's
provably valid (condition 2 holds) but only reachable via a chain that, at some
point, required passing through a position where a side had *already* played
into a position other, better non-losing choices existed to avoid — is real,
provably a boundary state, but not part of this specific, game-start-anchored
set. §4's terms still apply to it: it may still be load-bearing for some other
ancestor, just not for the start position itself.

**Why this is the concept actually worth calling "global," and why determining
it is a different problem from defining it.** Definitionally, this is precise
and well-posed — no ambiguity in what's being asked. Determining it for any
specific, already-confirmed boundary state runs into the same wall already
established two turns ago for the broader closure argument: the chain has to
be traced all the way back through the true starting position, crossing
successively larger, currently unsolvable materials along the way. The concept
is sound and sharply defined; checking it for a real, specific boundary state
is not currently within reach, for the same scale-gap reason as before, not a
new one.

**A fact that is not empirical at all, and should be stated as the near-
tautology it is, not hedged as merely plausible.** If the start position's
true value is draw, whoever is to move at any position along a line of
sustained mutual optimal play has, by definition, a value-preserving move
available — that is what "the position's value is draw" means. If both sides
always take that move, the value never changes, at any step, because nothing
has been played that could change it. So two players who always play their
own truly optimal move, starting from a drawn position, cannot reach a
decisive position — not empirically unlikely, structurally impossible, the
same way summing only even numbers can't produce an odd one. Any chain that
does connect a drawn start to a decisive position must therefore contain, at
some point, a move that was *not* the position's own value-preserving choice —
a genuine deviation, in the precise sense of departing from what that
position's own true value required.

**What this fact does not, by itself, establish — and why the correction below
still stands as its own separate claim, not weakened by the fact above.** The
deviation point guaranteed to exist is not automatically a *strict* boundary
state. It is a position where a non-losing option existed and a different,
losing one was taken instead — but that position may have had several
non-losing options, not exactly one, making it a transition point in the
broader sense (already established two turns ago) without satisfying
condition 2 specifically. So: "no start-reachable boundary state exists, under
the strict condition-2 definition" does not imply chess is forced from the
start. It would mean something narrower and still fully consistent with chess
being a draw: that every such guaranteed deviation point, wherever it occurs,
happens to have several non-losing options rather than exactly one — the
drawing mechanism staying wide, with real slack, rather than ever narrowing to
a single safe choice. That's a different, weaker finding than "no transition,"
and it doesn't license the forced-win inference. The absence of a strict
boundary state doesn't mean the absence of a transition at all — the two facts
in this section are complementary, not in tension: a transition is guaranteed;
whether it's ever the narrow, strict kind is what remains genuinely open.

**Where this leaves the actual range of outcomes — and why treating it as a
spectrum, not a dichotomy, is the more honest framing.** Not "either a forced
win is discovered, or the closure is confirmed" — a genuine, wider range:
strict start-reachable boundary states might be found, incrementally mapping
real structure; none might ever be found while the drawing mechanism turns out
to stay wide the whole way, itself a real and interesting structural fact
about the game; or the scale gap might simply make the question unanswerable
for the foreseeable future regardless of which is true. All three are learning
something real about chess's own structure that isn't currently known. None of
them require the dramatic, low-probability alternative to be the only other
option.

---

## 5a. A further proposed distinction — local versus global by closure contribution

Recorded as stated, then assessed separately, in the same way as §8. Written
after §§3-5 and in some tension with their use of the word "global"; the
tension is part of what the assessment addresses.

### 5a-i. The idea, as stated

Local boundary states are those that do not contribute to the closure of the
choice/force boundary; global boundary states are those that do. The trivial
motivating case is king and rook pawn against king with the pawn on the a- or
h-file: a localized draw from the corner race. It is characterized as a state
whose scenario would be difficult to extend — additional black pieces to
create structural immobility, or additional white pieces to complicate — so it
does not grow into a larger structure. Its importance is as a primitive: the
grander idea is structural locking mechanisms that depend on race conditions
to unlock, where unlocking flips a forced draw into a loss. A catalog of
boundary states lets that compositional nature be understood, whether a given
state is local or global.

### 5a-ii. Assessment — separate from the record above

*Correction (`boundary_board_states.md` §14 item 8): this section uses "exactly one non-losing move" as if it were the definition of a boundary state. It is a measure of width (w = 1 is the strict case), not the definition. Where a conclusion below depends on w = 1 specifically, it holds only for that strict case.*

**As currently worded, closure membership cannot be the criterion.** §3 and
Method 2 define a transition point as the first ancestor on a validated
retrograde walk where the loser's non-losing count moves from zero to nonzero.
A strict boundary state with at least two legal moves has at least one losing
alternative; that alternative leads to a decisive position; the walk from that
position arrives at the boundary state, which has a nonzero count. So every
strict boundary state that is a genuine choice-point (more than one legal
move) is already a transition point, hence already in the closure. Under that
reading every such state would be global, and the distinction would separate
nothing. This is a reading of the existing definitions and may miss a
refinement intended by the idea; if "closure of the choice/force boundary"
means something narrower than the §3 set, that narrower object needs to be
defined before the distinction can be used.

**The corner-race example argues from a different property.** What it shows is
*extensibility*: whether the lock can host more material without being broken
or breaking. That is separate from *upward connectivity*: whether non-losing
lines from richer material land on the state. The rook-pawn KPK draw is a
poor host but a very good landing point (larger drawing lines trade down into
it), so it can be non-extensible and strongly connected at once. Three
candidate meanings are therefore in play — extensible, upward-connected, or a
conjunction — and which one "global" is meant to carry has not been decided.

**Terminology.** "Local" and "global" already carry three senses in this
document (§3, §5, §6). A fourth sense is best given its own labels once the
definition is fixed. The tier structure in `boundary_board_states.md` §14
item 1 is a separate axis (how much of the position's history is assumed),
and should not be confused with this one.

**Open (also listed in §9):** which of extensibility, connectivity, or both is
intended; whether extensibility can be defined operationally (for example, as
the set of single-piece additions under which the exactly-one property and the
value survive); and whether either notion correlates with the load-bearing
relation of §4.

### 5a-iii. A later refinement, as stated, then assessed

*Correction (`boundary_board_states.md` §14 item 8): this section uses "exactly one non-losing move" as if it were the definition of a boundary state. It is a measure of width (w = 1 is the strict case), not the definition. Where a conclusion below depends on w = 1 specifically, it holds only for that strict case.*

**As stated.** A true, global boundary board state is one in which each player
has sufficient material and the best outcome for both players is a draw: an
equilibrium of race conditions. It is not the case of one side holding at best
a draw against a side holding at best a win, which is what a bare king against
a material advantage is. Perfect play would not deliberately navigate into
such a position unless it were the only option. The rook-pawn KPvK draw is a
*localized* boundary state, emergent from the lack of composition: adding
further primitives would interfere with its locking race too closely for
piece-placement geometry alone to manage. What makes it a boundary state is
the edge-of-board condition. The pawn is topologically restricted to the board
edge, and king routing is affected by the asymmetric occupancy of the edge.
That restriction acts like a board edge inside the locking mechanism (see §8k).

**Assessment — separate from the record above.**

- **This appears to answer the open question of §5a-ii.** "Global" is defined
  by material and equilibrium (both sides materially live, mutual-draw
  equilibrium), not by closure contribution or extensibility. "Local" is a
  draw that exists only because the edge does structural work and nothing else
  is composed with it. If this supersedes the closure-contribution wording,
  §5a-i and §5a-ii stay as the record of the earlier attempt, and the
  objection in §5a-ii still applies to that earlier wording.
- **It is consistent with condition 1** of `boundary_board_states.md` §1 and
  with the resolution in that document's §14 item 2. "Sufficient material"
  needs an operational form; proposed: *materially live*, meaning some legal
  position of the material split is a win for the side under optimal play.
- **Under this vocabulary the rook-pawn KPvK draw fails the both-live test**
  (bare king) and is therefore local. It is best treated as a *primitive*: a
  source of structure to compose from, not a global state itself. This
  matches the way §8 already uses it.
- **Open: do hanging-piece positions in total-win materials count as local
  boundary states?** They are drawn and exactly-one, and one side is not live,
  like KPvK. They differ in mechanism: the draw comes from reduction to
  insufficient material, and there is no locking race or edge condition. If a
  local boundary state must involve a locking mechanism, they fall outside.
  If "local" is only "exactly-one with one side not live," they fall inside.
- **The "only option" remark** corresponds to the Tier 2 requirement in
  `boundary_board_states.md` §14 item 1: reachable only through positions where
  neither side had a better alternative. Recorded as a candidate additional
  requirement for the global class, not adopted.

---

## 6. Coverage — scoped and local, not global, and why the difference matters

A related, separate question: can a growing catalog ever let you claim
certainty about a *region* of positions, without storing that region's full
landscape the way `boundary_board_states.md` §3 already explains isn't
necessary for individual verdicts?

**The circularity in the global version, stated precisely**: claiming "this
catalog now covers every boundary state in this entire material" requires
already knowing the complete set of decisive positions to check that claim
against — which is exactly the full-landscape knowledge this whole approach
exists to avoid needing. That circularity doesn't go away with a bigger
catalog; it's structural.

**The scoped version, which has no such problem**: "every position in this
specific, named set — these particular seeds, this particular starting
position of interest — has had its validated retrograde ancestry fully
resolved" is a claim with no circularity at all, because the set being claimed
complete is explicit and finite from the start, not "everything." This is the
same shape as §10's scoping mechanism in `boundary_board_states.md` — targeted
certainty about what was actually asked for, not a claim about the whole
material.

**Explicit bookkeeping, not an emergent property of catalog size**: getting
scoped coverage right requires tracking, deliberately, which specific positions
have had their full ancestry resolved and which haven't — a real data
structure decision, not something that falls out automatically from having
more entries in the catalog. This is the honest version of "we don't need to
carry the full landscape forward" — coverage of what was actually scoped, kept
explicit, not an implicit claim about everything the catalog happens to touch.

---

## 7. Scaffolding outward — an incremental strategy that doesn't require crossing the whole gap at once

§5 established that checking whether any *specific* boundary state is
start-reachable requires a chain all the way back through the true starting
position — currently out of reach in one leap. This section is a different,
practical question: not "can we prove start-reachability directly," but "can
we make real, incremental, prioritized progress toward it, using exactly the
machinery already built, without needing the whole chain at once."

*Note (survey §4c and §4f): the route to the start position is now also the forward trace of a pruned structure from the opening position, scoped to a chosen opening (`boundary_board_states_opening_survey.md` §4c). The backward scaffolding of this section remains, and survey §4f states reachability from certified positions as a backward direction beside the forward trace. The two meet at the boundary line of the catalog (survey §4f). The statement here that the chain cannot be traced in one leap is unchanged.*

**The strategy, stated precisely.** Take the full set of currently-confirmed
boundary states — everything already validated via Syzygy or this project's
own solved materials (`boundary_board_states.md` §2). From each one, extend
Method 2's validated retrograde walk further upstream than strictly necessary
for that state's own confirmation — into its own ancestors, and, via
uncapture/unpromotion, into the next tier of larger materials. Two things can
happen at each newly-reached ancestor, and both are informative:

- It connects — through some further validated chain — to *another*
  already-confirmed boundary state reached independently, from a different
  starting seed. This reveals that two states previously treated as separate
  are part of the same connected region of resolved territory — real,
  new structural information, not assumed going in.
- It runs into unresolved territory — a material too large for current Syzygy
  coverage or this project's own solved set. This is not a dead end to record
  and forget; it's a specific, named material now flagged as exactly where the
  current scaffold's frontier sits.

**What this produces**: not proof of start-reachability for anything, but an
explicit, growing map — which confirmed boundary states are connected to each
other, forming larger unified regions of understood territory, and which
specific materials sit exactly on the current frontier where that scaffold
runs out. That frontier list is a real, structurally-derived prioritization
target: not "which large material seems interesting to tackle next," but
"which large materials are exactly where the current connected scaffold stops,
such that resolving one of them would extend it." This is the ordering-only
use of structure §1 and §4 already sanction — it decides *what to work on
next*, never *what's already true*.

**The honest limit on the strategy itself, not just on what it proves.**
Confirming that two boundary states are genuinely connected is subject to the
same scale constraint as everything else here — it only succeeds as far as a
validated chain can currently be traced, which may itself run into unresolved
territory before reaching the other state. That's not a flaw; it's why the
strategy is incremental by construction rather than a single computation:
extend as far as current reach allows, record exactly where it stops, work on
that specific frontier, then extend again. Each cycle either grows the
connected region or narrows the confirmed frontier — real, cumulative,
never-recomputed progress (`boundary_board_states.md` §3), whether or not
start-reachability is ever fully settled for any given material.

**Why this is the concrete, technical form of the Lakatosian framing from
§2, not just a restatement of it.** The protective belt doesn't grow by
guessing at the whole structure in advance — it grows by testing specific,
current-reach extensions against what's already confirmed, keeping what
connects and flagging what doesn't yet resolve. This section is that process
applied directly to the boundary-state scaffold: build outward from certainty,
tier by tier, letting each extension inform the next rather than attempting
the whole span in one step.

**Precedent, checked directly rather than assumed — and a major one, not an
adjacent one.** This exact strategy — an exhaustively-solved floor, extended
outward by forward search that terminates the instant it reconnects to that
floor, each successful connection permanently expanding what's proven — is
how checkers was weakly solved: retrograde-built endgame databases combined
with proof-number search from the starting position, described directly in
that literature as using the database "to raise the water level of where
game play lines terminate." Eighteen years, a 3.9×10¹³-position database, the
single most successful result in this entire research area — arrived at here
independently, before this precedent was checked, not after.

**A precise correction worth stating plainly rather than leaving the
stronger claim standing uncorrected.** It is not that applying this strategy
to chess has gone unconsidered. The checkers literature addresses it
directly: even with the exact technique that solved checkers, chess's own
scale — roughly 10^44 legal positions, against checkers' roughly 5×10^20 — is
explicitly described as needing "a breakthrough such as quantum computing" to
attempt the same way. The field has already looked at this strategy for
chess, in its full, whole-game-solving form, and already identified the scale
gap (§5) as the blocking obstacle.

**Where the genuine contribution actually sits — and why this scoping is the
sophisticated move, not a repeat of an already-rejected attempt.** What this
section describes is not "run checkers' strategy on chess and hope it
scales" — that is precisely the version already known to fail. It is
narrower in a way that specifically avoids the objection that closed off the
naive version: anchoring to narrow-margin, condition-2 transition points
specifically rather than an undifferentiated exhaustive floor; explicitly
tracking connectivity between multiple independently-discovered anchors as
its own derived structure, rather than one forward line toward one database;
and, critically, the goal is incremental topology-mapping (§5), not
full-game solving — scoped explicitly to never require the whole chain in one
step. That reframing is exactly what avoids needing the scale breakthrough
the full version requires. As far as this document's own research effort can
determine, that specific combination — narrow-margin anchoring, explicit
cross-anchor connectivity tracking, and topology-mapping in place of full
solving — has not been found published anywhere.

**A further extension of the scaffold once it's large enough to be data, not
just structure** (`boundary_board_states.md` §2, Method 3): each confirmed
boundary state's own upstream ancestry — its cone, as built by this section's
own scaffolding process — is itself a comparable object once enough of them
exist. Comparing the shape of many such cones to each other, not just static
features of the states at their root, is the natural extension of Method 3's
population analysis to the structure this section builds. Explicitly a
candidate-prioritization signal, not a certainty method, per the same
distinction held throughout §1 and Method 3 itself.

### 7a. Composition as the specific mechanism that extends the floor past the storage ceiling

**The idea, as stated.** The forward-trace-plus-retrograde strategy already
described in this section is not offered as analogy alone; it is being used
directly, matched with boundary-state composition, as a way to sidestep
Syzygy's own storage ceiling — a genuine attempt at a pathway beyond current
methods, related to checkers' weak solution in direct application. What
distinguishes it from any prior consideration of applying the checkers
strategy to chess is that no prior attempt, in the checkers literature or in
chess's own, has used forward tracing from boundary-state topological
composition, with joint compatibility as described in §8g-§8m, as the
mechanism for extending the floor.

**Assessment — separate from the record above.**

- **What was already established, restated so the new point is visible
  against it.** §7 already argues that anchoring to condition-2 transition
  points and topology-mapping, rather than an undifferentiated exhaustive
  floor, is what avoids the scale breakthrough the naive checkers-on-chess
  attempt needs. §10 already establishes that a validated boundary catalog
  acts as an early-termination oracle for forward search. §12 already
  establishes that the boundary-catalog approach is a different shape of
  problem than Syzygy's own — resolved one candidate at a time and cached
  permanently, not stored all at once for a whole material.
- **What is newly made explicit here, and is a real, load-bearing point, not
  a restatement.** Composition (§8d-§8m) is specifically what lets the floor
  in §7's scaffolding strategy include certified positions whose material has
  more pieces than Syzygy, or this project's own search-based methods, could
  ever exhaustively store or search — because a composed boundary state's
  certificate needs only its primitives' own small, already-solved
  boundaries plus the composition argument, not storage or search over its
  own full material. §8d already makes a narrower version of this point about
  a single composition (two locked pairs already exceeds Syzygy's practical
  reach). Connecting it to §7's scaffolding strategy directly — that composed
  states are the specific mechanism for growing the floor past the
  piece-count wall the checkers literature already identified as chess's
  blocking obstacle (§7, quoted above: "a breakthrough such as quantum
  computing") — had not been stated together in one place before this turn.
  This is the precise technical content behind "relates to checkers directly
  in application, differently from what's been tried": the field's own
  assessment concerns a floor built by exhaustive storage; a floor built
  partly from composed certificates is not subject to that argument.
- **A real, related precedent found on checking this session, recorded
  honestly rather than left unsearched.** Bleicher's "Freezer" tool builds
  new tablebases from existing (Nalimov) tables using a priori information —
  specifically, blocked pawns — producing results for material with seven or
  more pieces even though full, unconstrained tablebases at that piece count
  are not generally available. This is independent, existing validation of
  the underlying intuition that locked or frozen structure is a real lever
  against the storage ceiling. It is a different mechanism from Method 5,
  worth being precise about rather than treating as the same idea already
  found elsewhere: Freezer narrows the search space using known-frozen
  structure and then computes exhaustively within that narrowed space; Method
  5 proves a specific position by argument, with no exhaustive computation
  even within the narrowed space. Recorded as the closest related prior art
  found so far, not as prior art for Method 5 itself.
- **The novelty claim otherwise stands at the same hedge §7 already uses.**
  This project's own effort — including one search this session, not an
  exhaustive review — has not found the specific combination described here
  (checkers-style scaffolding extended by certified composition rather than
  exhaustive storage) published elsewhere. That is the same epistemic status
  §7 already claims for its own narrower point, carried forward, not
  strengthened.
- **What this does not yet establish, at the same standard §7 and §13 already
  hold elsewhere.** That a composed floor of real, useful size can actually
  be built — §8m's own open question about whether enough primitives recur.
  That forward search from a position of real interest actually reaches such
  a floor in practice — §7's own honest limit on the scaffolding strategy,
  unchanged by this section. That any of this bears on chess's starting
  position specifically — §13's standing caveat. This section states *why*
  the storage ceiling specifically does not bind composed states; it is not a
  claim that the ceiling has been overcome in practice, and should not be
  read as one.

---

## 8. Composition from race-condition primitives — constructing boundary states by argument

This section records a construction idea in full, in the terms it was developed,
before assessing it. The assessment (§8d) is deliberately kept separate, so the
idea is preserved as stated and the analysis can be read, or disagreed with,
on its own.

### 8a. The idea, as stated

Boundary states need not only be *found* by search over a material. Some can be
*constructed* from small structures whose values are already established, so
that the constructed position's value follows from those structures' own
arguments. Three parts:

- **Frozen pieces as landscape.** A piece with zero legal moves, whose
  disturbance provably costs the game under perfect play, is not an active
  piece. It occupies space and constrains other pieces exactly as a board edge
  does. If the setup's race condition is composed correctly, making the frozen
  structure "unstuck" costs a loss, so it can be relied on as fixed terrain.
- **Race conditions as the certificate.** Breaking the structure requires a
  piece to leave a post that was covering something else. Whether that costs
  the game is a race — tempo and distance — of the same kind as Uri's rule and
  the rule of the square. If the inequality favors the defender at setup, the
  structure holds.
- **The target outcome.** A position with sufficient material on both sides in
  which the only non-losing move is a bishop move that neither breaks the
  structure nor abandons a post; repetition follows; every attempted
  breakthrough loses. That is condition 2 of §1 in `boundary_board_states.md`,
  reached by construction instead of by search.

**The worked example this has been developed on.** Locked pawn pairs — each
pawn faced by an enemy pawn, so zero legal moves each — with each side's bishop
on the same color as its own pawns and therefore opposite in color to the
enemy's. A king is posted at the base of the chain, because the base pawn is
the weak point; optionally a knight aims at the base, itself placed so it
cannot move without loss. A diagonal chain is automatically single-colored, so
for example white pawns a3–b4–c5–d6 (all dark squares) against black pawns
a4–b5–c6–d7 (all light squares) satisfy the color assignment across the whole
structure; the assignment can equally be flipped.

### 8b. The catalog, as a class rather than a single example

The construction needs primitives, and what is wanted is a *catalog of them as
a class* — known race-condition and zugzwang structures, each recorded as: a
small ending, the inequality that decides it, and the pieces and squares it
locks. Entry families:

- **Edge and corner draws:** rook-pawn king-and-pawn against king (the defender
  reaches the corner); wrong-colored bishop with a rook pawn (the Uri's-rule
  material); knight plus a rook pawn on the seventh rank (a known exception to
  the general knight-and-pawn win).
- **Pure race:** the rule of the square.
- **The opposition family:** direct, distant and diagonal opposition; key
  squares; triangulation; corresponding squares.
- **Spare-tempo parity:** a side with one more spare pawn move wins the tempo
  battle; the standard illustration turns on a single pawn's square (h6
  versus h7).
- **Base-of-chain attack:** the weak base pawn, and the king posted to guard it.
- **Deflection and breakthrough:** the outside passed pawn; a pawn sacrifice
  that creates a passer.
- **Blockade fortresses:** opposite-colored bishops and similar (wide-margin;
  see §8d).
- **Existing catalogs:** Haworth's 2001 enumeration of mutual zugzwangs across
  essentially all endgames up to five pieces is already-built raw material,
  with the caveat that it uses the toggled definition, so its entries need
  re-reading against the single-position boundary-state definition
  (`boundary_board_states.md` §9).

An expectation, not a finding: the edge, opposition and spare-tempo entries are
narrow-margin and the likeliest sources of exactly-one positions; the fortress
entries are wide.

### 8c. Process notes worth keeping

- The construction was first attempted on king-bishop-rook-pawn against king
  (the Uri's-rule material) and then self-corrected: that material leaves the
  attacker too little force to punish a detour into a pawn structure, so the
  tempo-margin argument cannot be made there. A base for composition needs
  enough attacking material to enforce the inequality.
- Approximations such as Uri's rule are legitimate for informing candidates and
  topological relations. The standing rule (§1 of this document) is only that
  they never supply certainty.
- Real opening theory was proposed as a source of "pseudo-game-start"
  positions: known lines that reach locked structures with reduced rook
  mobility, used as a plausible-corridor heuristic, not a claim of
  start-reachability. Whether perfect play reaches such a position is
  irrelevant to the construction's validity, since all drawing routes are
  equally valid.

### 8d. Assessment — separate from the record above

**What is sound.** Certifying a value by argument instead of search is how
fortress and opposition theory works, so the idea is not exotic. The color
arguments are real: a bishop is confined to one color permanently; a pawn's
capture step preserves color, so a pawn attacks only squares of its own color;
two pawns locked head-to-head are necessarily on opposite colors. With each
bishop on its own pawns' color, neither a bishop nor any pawn can ever attack an
enemy pawn or bishop, and the bishops cannot touch each other; in a position of
kings, bishops and pawns, only kings can win pawns. Zero legal moves is a pure
inspection fact.

**The proof obligation, stated explicitly.** "Breaking costs a loss by
definition" holds relative to a complete case list: for the side to move, every
legal move loses except the waiting move; the waiting move preserves the
situation; and the opponent has no winning try, sacrifices included. The
conditions that make that list short and checkable:

1. **Bishop color-lock**, as above.
2. **Tempo margin:** for each route that breaks the structure, the breaker's
   cost exceeds the race margin, with the +1 tempo adjustment for whoever is not
   on move — the same inequality shape as Uri's rule.
3. **Mutual blocking:** every pawn faced on its file, with no capture available,
   so each has zero legal moves.
4. **Capture-diagonal integrity:** the diagonal squares in front of locked pawns
   must stay free of enemy pieces. Enemy kings are excluded for free (stepping
   there is moving into check); an enemy bishop's route must be checked; a
   knight changes color every move, so it is not covered by the color argument
   at all — its inertness has to come from a pin or self-blocking.
5. **Waiting moves:** a side that runs out of safe waiting moves is in
   zugzwang and the position is a win, not a draw. Free bishops supply them.
6. **Stalemate:** freezing pieces removes legal moves and raises stalemate risk.
   A draw reached because the defending king is stalemated is a different
   mechanism from the race draw being generalized.

**Wide versus narrow is a design decision, not a detail.** A free bishop that
neither enemy pawns nor the enemy bishop can touch has many safe squares, so
many non-losing moves — a proven *wide* fortress. Exactly-one (condition 2)
requires that each alternative bishop move abandons something, which has to
come from the kings' placement restricting the bishops. Provable and narrow pull
in opposite directions, so the target has to be chosen.

**Spare tempo closes only when pawns are truly frozen.** A pawn with an
available push carries spare tempi and can decide a zugzwang battle without
being touched; the standard illustration turns on one pawn's square. Zero-legal-
move locking removes that channel, which is why the frozen structure is the
right choice.

**Uri's rule is a heuristic with documented gaps, not an exact boundary.** The
source giving the Chebyshev formula (chessprogramming.org, *Wrong Color Bishop
and Rook Pawn*) says this ending must be encoded in an evaluation function
because search does not evaluate it properly; lists a separate "too far" case
where the stronger side wins; notes that the heuristic requires the program to
recognize the draw even after the stronger side sacrifices its bishop; and its
footnotes cite a 2013 engine game in which one engine apparently lacked the
knowledge (as reported in that source; not independently verified here). The
formula is a good informer of candidates. For composition by proof, each
primitive used should come from its exact tablebase boundary for the small
material, which is trivially within reach, not from the formula.

**Opposite-colored bishops are a different concept from a wrong-colored bishop,
and not the natural source of narrow boundary states.** Wrong-colored: one
bishop, wrong for a rook pawn's promotion square. Opposite-colored: each side
has a bishop on opposite colors. The latter is the well-known drawish structure
in opening and middlegame theory, but its draws are robust blockades, and the
endgame literature states that zugzwang is less common there than with
same-colored bishops. A search of opening theory for such structures should
therefore be expected to find wide fortresses. The transition *into* the
blockade, before it is established, is the more plausible place for a
tempo-critical moment — a hypothesis to test, not an assumption.

**Verification stance.** For frozen-obstacle constructions the case-list
argument is the proof; a direct check is insurance against a case-list error,
not a gate on the idea. The check is small: a solver over kings, bishops and at
most a knight, on a board with blocked squares, needs no Syzygy. The minimal
instance — one locked pair, two bishops, two kings — is six pieces and inside
Syzygy; two locked pairs is eight and beyond it, which is where argument or
this project's own engine has to take over.

### 8e. Where this sits among the methods

It is the fifth generation method (`boundary_board_states.md` §2, Method 5),
unusual in one respect: its certificate is a constructive argument rather than a
table lookup or a search, provided the case list is complete. It feeds the same
condition-2 requirement as the other four.

The sections below (§8f-§8j) extend §8 with ideas developed after it was first
written. Each records the idea in the terms it was stated, then assesses it
separately, following §8a/§8d. They are a design record for an idea that is
not yet fully crystallized, not a set of results.

### 8f. Locks generalized: duties, not only frozen walls

*Correction (`boundary_board_states.md` §14 item 8): this section uses "exactly one non-losing move" as if it were the definition of a boundary state. It is a measure of width (w = 1 is the strict case), not the definition. Where a conclusion below depends on w = 1 specifically, it holds only for that strict case.*

**The idea, as stated.** Frozen pawn walls are one instance of a wider
locking mechanism. In advanced play a king with a fortress has coverage; the
opponent, holding an extra piece, places it undefended near the fortress so
that it looks like a plain sacrifice. Taking the piece removes a defender from
the post that was protecting against the mating threat, breaking the race
condition that held the fortress. Accepting the "hung piece" therefore loses,
and the offer is a deliberate move to disturb one precise race-bound
topological relation. This is the same phenomenon as the pawn walls,
generalized: locking mechanisms bound by race conditions, whose disturbance
is the source of forced results, and whose primitives are what the
boundary-state catalog records. Forced mate can be derived from calculation on
these relations when the primitives are identified and composed.

**Assessment — separate from the record above.**

- **A lock is a duty.** A piece committed to a post is locked when every move
  that abandons the post loses, whether or not it has legal moves. A frozen
  piece is the degenerate case, enforced by having none. Chess names the
  attacker's tools against duties (deflection, decoy, overloading, removing
  the guard); the poisoned offer is a decoy.
- **It connects to the exactly-one count directly.** The single non-losing
  move of a boundary state is the one move that preserves every duty. Each
  alternative loses because it abandons something; tempting captures of
  poisoned pieces are a natural source of those losing alternatives. This is
  the dynamic form of §8d's "each alternative bishop move abandons something."
- **It sharpens the geometry/topology distinction.** The same geometry, a
  piece en prise to the king, has opposite values in two cases: when the
  capture is the only non-losing move (an undefended queen next to a king in
  KQvK, capture reaches a draw), and when the capture breaks a lock (the
  grandmaster case, capture loses). Nothing readable from placement separates
  them; only the value of the children does. This is a clean illustration of
  the point made in `boundary_board_states.md` §2 Method 4.
- **Overload is the composition failure mode.** Two locks are compatible only
  if no single piece has to serve both duties. Breaking a composition by
  overloading a defender is the attacker's method, so "what does this
  primitive share" is a mandatory catalog field (§8i).
- **Certification differs from the frozen case.** A frozen wall is certified
  by inspection plus a short case list. A duty-based lock's case list must
  cover sacrifices, checks and other forcing lines that extend the footprint.
  §8h states the conditions under which that can be done by construction
  instead of by search over the composed position. The grandmaster's read is
  a strong way to *find* candidates; it is a prioritization aid under §1, and
  refuted sacrifices show it can fail.
- **Footprints grow with the number of duties.** A fortress depends on every
  piece covering its squares, so the scope around it is larger than around a
  pawn wall, and the saving from scoped verification (§8h) shrinks accordingly.

### 8g. Composition as part-wise joint compatibility and complete interdependence

**The idea, as stated.** Composition is part-wise: primitives a, b and c are
composed so that a with b, a with c, and b with c are each properly compatible,
with their race conditions inside the threshold for a draw. Breaking any of
these relations unlocks the structure and gives the opponent a tempo to
exploit it. The joint compatibility is complete-graph in nature: the
compatibility of a with c *emerges with* b; that of a with b, with c; that of
b with c, with a. Introducing a fourth primitive d extends this: a/b/c emerges
with d, a/b/d with c, and so on. All of them are interdependent in a single
locking mechanism in which the race conditions form a tie, and the tie
resolves only if a side chooses to break a lock, at which point perfect play
loses. Any imbalance outside this equilibrium would, by definition, mean the
position is not a boundary state. From this view, game-theoretic race
conditions become visible that piece-placement geometry alone would hide.

**Assessment — separate from the record above.**

- **The picture, restated.** A boundary state is a self-supporting set of
  locks, like an arch. Each primitive's guarantee depends on assumptions that
  the other primitives discharge; the interface between any subset determines
  what the remaining primitive must supply; the exactly-one move is the only
  one leaving every lock in place [Read with `boundary_board_states.md` §14 item 8: the equilibrium moves keep every lock in place; there need not be exactly one.].
- **Two testable consequences.**
  1. *Every primitive is load-bearing.* This is the dominator relation of §4
     applied inside a single position rather than across a catalog. Removing
     any one primitive (the ablation of §8a) should change the value or the
     non-losing count. A primitive whose removal changes neither is decoration,
     and the composition is not tight.
  2. *The set is closed.* Every piece is either inside a primitive or frozen
     terrain, and no interaction path leaves the set. This is the case-list
     completeness requirement of §8d stated as a structural property.
- **Pairwise compatibility does not imply joint compatibility.** Example:
  a shared resource of capacity 2 (say, the defending king's tempo margin)
  that each of a, b and c draws 1 unit from. Every pair fits (1 + 1 ≤ 2), all
  three do not (3 > 2), and a deviation that draws on the resource for all
  three breaks the composition though no pair was ever in conflict. Joint
  compatibility therefore has to be checked as a constraint system over all
  shared resources at once, not as a set of pairwise checks. The
  "emergent" complete-graph structure is real, but it is higher-order, not
  pairwise.
- **Mutual dependence needs a base case.** If a is safe because b holds and b
  holds because a is safe, nothing has been established unless the circle is
  anchored. Circular assume-guarantee reasoning is sound only with a
  well-founded ordering. Tempo supplies one, since each deviation costs moves
  and the argument can be by induction on the tempo remaining; the anchors are
  the exact tablebase-derived primitives and the zero-legal-move terrain.
  Everything else must be discharged relative to those.
- **"Anything outside the equilibrium would make it not a boundary state" is
  the definition, not the certificate.** It says what a boundary state is. It
  does not show that a particular composition is one. The certificate is the
  closure argument: showing that no outside imbalance exists.
- **Complete interdependence is the signature of a tight, irredundant state,
  not a requirement of the definition.** A boundary state carrying extra
  redundant structure can also exist; the load-bearing test distinguishes the
  two.
- **Sharpened after a direct challenge to this exact point, because the two
  questions it conflates need to stay separated by name.** *Question one:*
  given a position already verified to be a boundary state (by tablebase or
  by a completed §8d certificate), are its composing primitives jointly
  compatible? Yes, trivially, by definition — a verified boundary state's
  exactly-one-non-losing-move property *is* what joint compatibility means,
  so there is no empirical question left to ask here, and there never was.
  This was already the content of the bullet above. *Question two:* given a
  set of primitives individually verified only in their own small material,
  each not yet tested together, will composing them produce a boundary
  state? This is not answered by the definition at all, because the
  definition specifies the target property, not a procedure for predicting
  whether an untested combination reaches it. This is the question every
  open item in §8g-§8p (the resource-accounting constraint, the base-case
  requirement, the population scan, the storage-ceiling prediction) is
  actually about, and it is the engineering problem of the programme, not
  because joint compatibility is ill-defined, but because *verifying that a
  specific proposal satisfies a well-defined target* is a different task
  from *defining the target*, and only the second one is settled.

### 8h. Scope of verification: footprints and interference closure

*Correction (`boundary_board_states.md` §14 item 8): this section uses "exactly one non-losing move" as if it were the definition of a boundary state. It is a measure of width (w = 1 is the strict case), not the definition. Where a conclusion below depends on w = 1 specifically, it holds only for that strict case.*

**The idea, as stated.** Deviations are scoped by the topological relations of
the composition, so exhaustive coverage is bounded by topological closure.
Dynamic locks need not be checked exhaustively if the primitive topological
relation is understood and the position is built from a composition of race
conditions. The problem is not a search problem; it is a composition of
discrete topological relations with calculable, interdependent race
conditions. Given a correct composition, any deviation fails the race by
construction, so the equilibrium holds.

**Assessment — separate from the record above.**

- **This has the shape of established compositional reasoning.** It is
  assume-guarantee verification: each part is verified against explicit
  assumptions about its environment, and the composition is justified by
  showing the assumptions are discharged. The related decomposition idea in
  chess is Noam Elkies' application of combinatorial game theory to endgames
  (arXiv:math/9905198), which splits positions into nearly independent
  subgames with values under optimal play; the tempo counting across
  components there is close to the race interdependence described here.
  Elkies notes that CGT does not apply to chess directly, because the winner
  is not determined by who moves last and games can be drawn by infinite
  play, and the decomposition works best in endgames, where the absence of
  long-range pieces reduces interdependence between regions. Both limits bear
  on this idea (see the fourth condition below).
- **Footprint.** Each primitive has a footprint: the squares and pieces whose
  moves can affect its race. A deviation matters only if it reaches the
  footprint, so the exhaustive check runs over the footprint and its overlaps
  with other primitives, not over the board. The tempo margin of the race
  bounds how far a deviating piece can matter (a light-cone effect), which
  is what makes the footprint finite.
- **The footprint must be defined conservatively.** Four conditions:
  1. *Forcing moves.* Checks, captures and threats gain tempo, so the
     footprint includes every forcing sequence that extends reach, sacrifices
     included.
  2. *Line closure.* Sliding pieces and discovered lines act at a distance
     independent of tempo; the footprint is closed under lines of attack, not
     only king-step distance.
  3. *Kings are the global resource.* Their reach is bounded only by the
     tempo margin, so they are where footprints overlap most. If every
     primitive needs the king, the composed scope may become the whole board
     and the saving disappears; composition pays off when overlaps are small.
  4. *Inertness outside the scope must be proven.* A piece outside the
     footprint must have zero legal moves and be unable to gain any as the
     in-scope pieces move, or its influence must be shown negligible.
- **The exactly-one count is not scoped.** It needs every legal move of the
  side to move, across the whole board. It is cheap when the out-of-scope
  pieces are frozen, and the proof that each alternative loses is what the
  scoped argument supplies.
- **Where the exhaustive work now sits.** Not in the composed position, but in
  (a) deriving each primitive exactly, from the exact tablebase boundary of a
  small material (cheap, done once), and (b) validating the composition rules
  on many small compositions before trusting them on large ones. The obstacle-
  board solver of §8d is the instrument for (b). It tests the composition
  theory; it is not a per-position gate.
- **Two senses of "closure" are in use and should be named apart.** The
  *interference closure* here (the smallest set of squares and pieces
  containing a primitive's footprint and closed under the ways moves can
  affect it) is a different object from the *transition-point closure* of §3
  (the set of all first-nonzero ancestors of decisive positions). The first is
  a scoping device inside a position; the second is a property of a
  material's whole landscape.

### 8i. The catalog schema for primitives (proposed)

*Correction (`boundary_board_states.md` §14 item 8): this section uses "exactly one non-losing move" as if it were the definition of a boundary state. It is a measure of width (w = 1 is the strict case), not the definition. Where a conclusion below depends on w = 1 specifically, it holds only for that strict case.*

*Withdrawn: the composition constraint that exactly-one needs one live primitive and every other piece frozen. It followed from the strict reading. With width allowed, two or more live primitives are compatible with an equilibrium, and composing one live primitive with frozen terrain is one option, not a requirement.*

*Implemented and extended in `boundary_board_state_primitive_catalog.md` §2; this section is kept as the original proposal.*

Extending the entry format of §8b (a small ending, its deciding inequality,
and what it locks) with the fields that composition requires. Proposed, not
yet exercised on any real entry.

- **Identity:** name; material; source of the deciding inequality (tablebase
  boundary, not a heuristic formula, per §8d).
- **Duties:** for each committed piece, the post, the squares or lines it
  covers, the tempo cost of rerouting, and what abandoning the duty costs.
- **Resources used:** each shared resource with the number of units drawn
  (king tempo margin, specific pieces, squares, spare pawn tempi).
- **Assumptions on the environment:** what other pieces must not do, stated
  as conditions (no tempo-gaining forcing move into the footprint, and so on).
- **Guarantee:** what holds if the assumptions hold, including the loss that
  follows a break.
- **Mobility signature:** frozen (zero legal moves) or live, and the count of
  safe moves. Exactly-one composition needs one live primitive and the rest
  frozen, since two live primitives with a safe move each yield a wide
  position (§8d).
- **Footprint:** the interference closure (§8h) of the primitive.
- **Forcing-move entry points:** the ways a forcing sequence can cross into
  or out of the footprint.
- **Extensibility and connectivity notes:** what additions it tolerates and
  what larger structures trade down into it (§5a).

Composition rules to be validated on small instances before use:
1. Per-resource accounting: for every shared resource, the sum of demands over
   the whole set of primitives does not exceed capacity.
2. Non-interference: no forcing move in one footprint gains tempo inside
   another.
3. Closure: every piece is inside a primitive or is proven frozen terrain.
4. Well-founded ordering: any circular dependence is discharged by induction on
   tempo from exact anchors.
5. Load-bearing test: ablating any primitive changes the value or the
   non-losing count.

### 8j. Path to implementation (provisional)

*Correction (`boundary_board_states.md` §14 item 8): this section uses "exactly one non-losing move" as if it were the definition of a boundary state. It is a measure of width (w = 1 is the strict case), not the definition. Where a conclusion below depends on w = 1 specifically, it holds only for that strict case.*

*Withdrawn: the composition constraint that exactly-one needs one live primitive and every other piece frozen. It followed from the strict reading. With width allowed, two or more live primitives are compatible with an equilibrium, and composing one live primitive with frozen terrain is one option, not a requirement.*

The idea is not yet fully crystallized, so this is an ordering of likely
steps, not a commitment. **Superseded in part by §8m**, which argues steps 2-3
as originally drafted here (starting from one hand-analyzed primitive) are
premature — the corrected order puts a population scan first. The steps below
carry that correction; the struck-through framing is kept for the record.

1. Decide the open definitional questions (`boundary_board_states.md` §14,
   and §5a here): which tier of boundary state is the target, and what
   "local versus global" is to distinguish.
2. **(Population-first, per §8m.)** Scan every material already covered by
   Syzygy or this project for *local* boundary states (§5a-iii) — drawn,
   exactly-one, one side failing the materially-live test — across as many
   distinct materials and entry families (§8b) as the scan reaches, not only
   rook-pawn KPvK.
3. Run minimization-by-ablation (§8l) across the resulting population and
   record each result's interface (§8i); look for what recurs.
4. Propose catalog primitives from what recurred in step 3. Where nothing
   recurs beyond rook-pawn KPvK, treat that as the answer to §8m's own open
   question about population size, and fall back to the single-example path
   (original steps 2-3 below) as the only available one.
5. *(Original single-example path, kept for the case step 4 falls back to it.)*
   Write one primitive as a contract in the §8i schema — the rook-pawn corner
   race is the natural candidate, and is a *live* component (the pawn has
   pushes), so composing it needs every other piece to be frozen terrain that
   leaves the race untouched, or a variant with the pawn itself blocked and
   the race played over something else. Derive its deciding inequality from
   the exact tablebase boundary of the small material.
6. Compose a proposed primitive (from step 4, or step 5's fallback) with one
   frozen structure (the locked pawn chains of §8a) on paper: write the
   resource accounting, the footprint, and the case list.
7. Run the small obstacle-board solver, and Syzygy where the piece count
   allows, as insurance; then run the load-bearing ablation on the result.
8. Only after several small compositions behave as the rules predict, attempt
   a duty-based (non-frozen) lock, which needs the full §8h conditions.
9. Feed each verified composition, and each failure, back into the catalog with
   the rule it confirmed or broke.

Step 1 and the scan in step 2 need no new code beyond the guarded scan in
`boundary_board_states.md` §6, restricted to the liveness test of §5a-iii;
step 3 and steps 6-7 need the obstacle-board solver already called for in
§8d, which §8m notes is required infrastructure here, not insurance only.

### 8k. The board edge as a structural requirement: role substitution

**The idea, as stated.** The rook-pawn KPvK draw is a boundary state because
of the edge-of-board condition: the pawn is topologically confined to the edge
and king routing is affected by asymmetric edge occupancy. The primitive
therefore carries a structural requirement, a restriction on piece routing
that acts like a board edge inside the locking mechanism. If nothing in a
composition supplies a compatible topological relation that fulfills what the
edge fulfills, the primitive is not viable there. If another topological
relation in the composition does account for that requirement, restricting
the routing of the same pieces, that opens an ablation: solve for an emergent
board geometry that meets both the race-condition requirements and the
topological requirements. Race-condition compatibility is the joint
compatibility of the composition, and when it fits, it fits perfectly.

**Assessment — separate from the record above.**

- **In assume-guarantee terms the edge is an assumption the board discharges
  for free** ("squares beyond this line do not exist"). Composition then asks
  for another structure whose guarantee discharges the same assumption inside
  the board. In the §8i schema the primitive lists the edge restriction under
  assumptions on the environment; a candidate substitute lists it under
  guarantees. This is the formal content of role substitution.
- **What the edge supplies should be itemized before substituting** (to be
  verified against the tablebase, not assumed): unavailable squares for both
  kings; an asymmetry in which king benefits, depending on placement relative
  to the edge and the corner; a constraint on the pawn (one capture direction,
  one neighboring file); and inertness, since the edge cannot be attacked,
  captured, or used to give check.
- **A substitute is weaker than an edge in a known way.** Pieces used as walls
  can be captured unless defended or unless capturing loses, can give or block
  checks, can create stalemate possibilities, open capture diagonals, and draw
  on the same king-tempo resources. A substitute is therefore itself a lock,
  either frozen terrain with its own certificate (§8d) or a duty (§8f). An
  edge is also symmetric in geometry, while a wall of pieces is oriented, so
  the asymmetry must be matched deliberately.
- **Hypothesis (assessment's own, untested): the edge also funnels.** On an
  open board the king's Chebyshev metric gives many equally short routes, so
  a tempo-preserving move is usually not unique. Exactly-one requires a unique
  tempo-preserving move, and uniqueness generally comes from constraints that
  remove alternatives. If so, a substitute must funnel as well as restrict,
  and these are two separable requirements. A direct test: in tablebase KPvK,
  check whether exactly-one positions concentrate near the edge. **Tested in
  §8t with an exact solver of king and pawn against king: not supported.** The
  fraction of drawn positions with exactly one saving move is lowest on the
  rook-pawn file and highest on the d-file, and the defending king is on an
  edge square about as often in those positions as in drawn positions
  generally. This hypothesis is withdrawn.
- **"Fits perfectly" corresponds to zero slack.** A race at the boundary is a
  tie: the mover's non-losing move preserves the margin exactly and each
  alternative costs at least one tempo. A wide fortress has slack. Zero slack
  is close to necessary for exactly-one but not sufficient, because the
  tempo-preserving move must also be unique (previous point).
- **The ablation as a solve.** With substitutes, the search for an emergent
  geometry is a constraint problem: variables are piece placements; constraints
  are each primitive's assumptions discharged by another's guarantees, the
  resource accounting of §8g, and the exactly-one count. The first cases are
  small enough for enumeration, consistent with steps 2-4 of §8j.

### 8l. Discovering primitives from tablebases, and generating candidates by subtraction

*Correction (`boundary_board_states.md` §14 item 8): this section uses "exactly one non-losing move" as if it were the definition of a boundary state. It is a measure of width (w = 1 is the strict case), not the definition. Where a conclusion below depends on w = 1 specifically, it holds only for that strict case.*

**The idea, as stated.** The rook-pawn KPvK primitive consists of: the
defending king's legal-move exclusion under the Chebyshev metric (supplied by
the board edge); the white pawn's influence at that edge; and the attacking
king's race to navigate the edge "pseudo-structure," which is too long to
prevent the draw, whether by race condition or by board constraints. On a real
board the same role can be played by inert, locked structure that provides the
same exclusion criteria, so the primitive is identifiable as a composition when
observed there. In any other boundary state built this way, what happens is
choice on non-edge pawn routing and exclusion of the other side, with the
emergent property arising from interdependencies among the other
compositional topological relations. Take the n compositions in such a state,
remove the KPvK primitive (n-1), and compose what remains: the closure
question is then the open topological requirement the removal leaves behind.
Any addition to a relation must preserve equilibrium across the complete graph
of relations. The joint compatibility must hold on the stated race conditions,
and if it does not, the position is not viable even as a boundary-state
candidate.

The program that follows: discover primitives from known tablebases such as
Syzygy; build a catalog of boundary states; classify each by the composition
of its topological relations; determine the number of compositions; deduct
relations from the primitive and boundary-state catalog; and solve the
resulting compositional problem to generate a new candidate.

**Assessment — separate from the record above.**

- **The program, stated as operations.** (1) *Scan* a solved material for
  boundary states (`boundary_board_states.md` §6, with the §14 guards).
  (2) *Decompose* each into primitives with explicit interfaces. (3) *Catalog*
  primitives and decompositions. (4) *Subtract* one primitive from a catalog
  entry and compute the residual: the assumptions it discharged that are now
  unmet. (5) *Solve* for a structure whose guarantees discharge the residual,
  within resource capacity, and verify the composed candidate by the
  certificate of §8d.
- **The decomposition operator is the crux, and it is not yet defined.** One
  candidate definition: minimization by ablation. Replace each piece in turn
  with nothing or with an inert wall, and test whether the property of interest
  survives (value, the exactly-one count, and the drawing mechanism). A
  primitive is a minimal fragment that keeps the property when everything else
  is removed or walled off. Syzygy cannot evaluate positions with walls, so
  this needs the obstacle-board solver of §8d; it is required infrastructure,
  not insurance only.
- **Decompositions will not be canonical.** A state can admit several minimal
  fragments, and removing a piece changes the value through interactions as
  well as through its own role. The catalog should therefore store the
  *interface* (assumptions and guarantees, §8i) of each fragment, not just its
  piece set, and record when two decompositions of one state disagree.
- **The generation operator can be named: lifting.** Convert an assumption
  that the *board* discharges (an edge, a corner, a color complex) into one
  that *pieces* discharge. For KPvK, this replaces the edge by locked
  structure. Which assumptions of a primitive are board-discharged and which
  are piece-discharged is then a required catalog field, because the
  board-discharged ones are the liftable ones.
- **The KPvK contract needs one more item.** The corner draw is realized
  through blockade and stalemate, and stalemate requires that the defender
  have no legal move other than the restricted king moves. The contract must
  therefore include *mobility exclusivity*: every other defending piece is
  frozen. This is a requirement, not a convenience, and it is why stalemate
  risk (§8d, item 6) is a first-order issue for any composition containing
  this primitive. It also suggests that locked pawn pairs (§8a), whose members
  have zero legal moves for both sides, are natural candidates for the
  substitute structure. A caveat to check: an enemy locked pawn next to the
  king is capturable unless defended or unless capturing loses.
  (Hypothesis, untested.)
- **The interface match is a necessary filter, not a certificate.** "If the
  joint compatibility does not hold, it is not even a candidate" is sound as
  pruning, in the same role as Method 1's structural filter. Passing it
  proves nothing; the candidate still needs the §8d certificate. This keeps
  the standing rule of §1.
- **Counting compositions.** The number of candidate compositions is finite
  for a fixed piece budget and can be bounded by a compatibility relation on
  the catalog. That relation is higher-order (§8g), so the enumeration runs
  over subsets with a global resource check, not over pairs.
- **What decides whether the program has leverage.** Whether boundary states
  in real tablebases decompose into a modest, recurring set of primitives. If
  they do, the catalog compresses and subtraction generates useful candidates.
  If nearly every state is its own primitive, the catalog does not compress and
  this approach reduces to enumerating states. The same data tests the
  coverage hypothesis raised about the start position (whether compensation
  for material imbalance in drawn positions reduces to composable locks).
  Nothing here has been run.

### 8m. Population-first extraction — the methodology's own justification

*Correction (`boundary_board_states.md` §14 item 8): this section uses "exactly one non-losing move" as if it were the definition of a boundary state. It is a measure of width (w = 1 is the strict case), not the definition. Where a conclusion below depends on w = 1 specifically, it holds only for that strict case.*

**The idea, as stated.** The decomposition process cannot be understood in
advance, independent of the structure of the relationships being composed;
joint compatibility is an emergent phenomenon, visible only once the
race-condition calculations behind it are understood in real cases. What is
needed first is a population of proven, localized boundary states, and the
primitives are to be identified *from* that population, not theorized ahead
of it. Once this novel structure — the choice/force boundary under perfect
play — is understood from first principles in real examples, the principles
become extractable primitives, and what composition looks like can be
characterized beyond the abstract complete-graph description already given.
This is offered as a point of certainty for the methodology, beyond
tautology, because it leverages a structure already tied to real precedent:
Uri's rule, as an existing geometric proxy for a topological fact; Možina,
Sadikov and Bratko's mining of the KBNK tablebase for structure beyond raw
values; and the checkers weak solution, as a retrograde floor extended by
forward search reconnecting to it.

**Assessment — separate from the record above.**

- **The core methodological point is sound, and it resolves a gap.** §8l left
  the decomposition operator undefined. This supplies the missing piece:
  don't define it a priori, extract it inductively from a corpus of already-
  confirmed cases. That is good empirical practice, and it is a real extension
  of Method 3 (`boundary_board_states.md` §2) — Method 3 as originally stated
  mines a confirmed population for features that order which *existing*
  unverified candidates to try first; this uses the same population to
  propose decomposition structure for *constructing new* candidates. The
  extension changes what the population is used for, not what it is allowed
  to certify: Method 3's own scope statement — a prioritization method, never
  a certainty method — still governs. A composition assembled from extracted
  primitives is a proposed candidate; it still requires the full §8d
  case-list certificate before it counts as confirmed. Nothing here removes
  that requirement.
- **Two different kinds of certainty are in play here, and keeping them
  separate is what the "beyond tautological" claim needs to rest on.** What is
  deductive and already established, independent of any population: once
  primitives are proposed, joint compatibility requires resource accounting,
  closure, and a well-founded ordering for any circular dependence (§8g,
  §8h). This follows from the minimax framework already adopted and needs no
  data. What is inductive and still open: which primitives actually recur
  across real materials, and what their interfaces look like. Only a
  population supplies that, and it remains a discovery process — real,
  useful, and still bounded by §1's standing rule that heuristics order and
  propose, never certify.
- **The three-precedent synthesis is a design rationale, not new evidence.**
  Combining checkers, Možina et al., and Uri's rule into one narrative does
  not raise any one of their individually-assessed strengths as recorded in
  `boundary_board_states.md` §9 — checkers stays a tight structural match,
  Možina et al. stay a loose family resemblance, Uri's rule stays a heuristic
  with documented gaps. What the synthesis adds is that three independent
  lines point the same direction: that solved chess data contains structure
  beyond its raw values, extractable and usable. That is a real, if modest,
  addition — an argument for trying the extraction, not a result of it.
- **The Možina et al. comparison, refined twice now — first corrected, then
  sharpened into the form that actually holds.** The first correction (above)
  addressed the claim that they attempted boundary-state composition and fell
  short from lacking this project's framework, and rejected it: KBNK has zero
  boundary states, so there was nothing in their material for a composition
  method to find, and that specific claim doesn't hold. The sharper claim
  made afterward is different, and does hold: not that they attempted this
  and fell short, but that their own project's structure made the connection
  to the checkers-style forward-search methodology (§7, §7a) invisible to
  them, for two compounding reasons, neither of which is a deficiency in their
  work — both follow from what their project was.
  - **Their goal never called for it.** They were explaining an already-
    exhaustively-solved material to humans. Nothing in that task requires
    extending search past a known floor into unsolved territory, because
    there is no unsolved territory adjacent to a material that is already
    completely known. The forward-search-to-floor idea has no purchase where
    everything is already on the floor.
  - **Their material couldn't have shown it even if the question had been
    posed.** KBNK has zero boundary states (as established in the correction
    above), so no version of their own project, run on their own case study,
    could have surfaced a boundary-state-based floor to search toward.
  Put together: the specific object that connects "mine a solved tablebase
  for recurring structure" (which they did) to "use that structure as a
  checkers-style floor for search beyond the storage ceiling" (§7a) is a
  sparse, individually-certified set of anchor points that doesn't require
  full-material storage — a boundary-state catalog. That object was outside
  the scope their project could have engaged with at all, given both its goal
  and its choice of material. This is the relation worth recording: not a
  failed attempt at the same task, but a real, explainable structural
  blindness — the piece needed to see the checkers connection as applicable
  beyond an already-solved material was not present in, or reachable from,
  their project as it was actually constituted.
  - **A note on this project's own history, recorded as stated and not
    independently verifiable here.** The account given is that the
    forward-trace approach matching the checkers methodology was derived
    directly from the abstract structure of the topological compositions,
    independent of and prior to encountering Uri's rule, Možina et al.'s
    work, or the checkers precedent itself — each found afterward, on
    checking what already existed. Recorded in the same spirit this
    document's own preamble already records real revision as part of the
    project's history, not as something verified here.
- **The concrete correction to build order.** §8j's step 2 and
  `boundary_board_states.md` step 13, as drafted, start from a single
  hand-analyzed example (rook-pawn KPvK). Under this idea that is premature:
  a single example cannot reveal what recurs. The corrected order —
  1. Scan every material Syzygy already covers (and any this project has
     solved) for *local* boundary states specifically — drawn, exactly-one,
     and failing the materially-live test on one side (§5a-iii) — not just
     rook-pawn KPvK, but the other entry families in §8b that exist as real
     tablebase positions: wrong-colored bishop with a rook pawn, knight with a
     rook pawn on the seventh rank, and any others the scan turns up.
  2. Run the minimization-by-ablation operator (§8l) across that whole
     population, not on one case, and record each result's interface (§8i).
  3. Look for what recurs across interfaces — shared exclusion mechanisms,
     shared race-margin shapes, shared mobility-exclusivity requirements — and
     only then propose these as catalog primitives.
  4. Attempt composition (§8j steps 4-7) using primitives that recurred in
     step 3, not primitives defined from a single worked example.
  This reordering is reflected in the updated build steps below and in
  `boundary_board_states.md` step 13.
- **What decides whether this changes anything, concretely.** Whether the
  scan in step 1 above turns up enough local boundary states, across enough
  distinct materials, for step 3's recurrence-finding to have anything to
  work with. If the true population of tablebase-confirmed local boundary
  states is small (rook-pawn edge draws being close to the only clean
  instance at very small material), the population-first approach still
  applies in spirit but has little to extract from until larger materials are
  reachable. Untested, and worth checking early, since it bears directly on
  whether step 1 is worth running before step 2's machinery is built.

### 8n. Uri's rule as a worked case: checked against its own published formula

**The idea, as stated.** Uri's rule may itself already be a composition — the
KPvK-style topological relationship joined with another, unknown one — an
instance of this project's own compositional operation, discovered by others
without the formal apparatus and simply given a name for heuristic engine
purposes, with that primitive understanding still in real use in engines
today.

**Assessment — checked against the published rule itself, not left as
conjecture.** Uri Blass's rule (as recorded on chessprogramming.org, "Wrong
Color Bishop and Rook Pawn") compares Chebyshev distances to the pawn's
promotion square: the defending king's distance, tempo-adjusted for whose
move it is; the attacking king's distance; and the pawn's own distance to
that square, capped to account for its double-step option. The bishop does
not enter the distance calculation at all. Its only role is established
separately: because it is the wrong color for the promotion square, it can
never help control or contest that square, and so it drops out of the race
entirely.

- **This supports a specific, sharper decomposition than "some unknown
  primitive."** The race itself is not the KPvK edge/corner draw as such —
  it is the already-separately-catalogued "rule of the square" family (§8b,
  listed under "pure race," distinct from the "edge and corner draws" family
  that wrong-colored-bishop material itself sits under), generalized to
  include the attacking king's own distance to the target square alongside
  the pawn's, rather than the pawn alone. The wrong-colored bishop composes
  with this race not as a second race but as a color-lock static fact — the
  same kind of argument already established in §8d, that a bishop is
  confined to one color permanently — which removes it from the race rather
  than contributing a second race condition of its own.
- **This is a genuine, checkable instance of the pattern §8f-§8k describe.**
  A live primitive (the generalized king-vs-king-vs-pawn race) composed with
  an inert one (the bishop's color-lock), where the inert component occupies
  space and cannot affect the outcome, exactly as §8f's frozen-terrain
  framing requires. Unlike the general population-first program of §8m,
  which is still speculative pending a real scan, this specific decomposition
  is directly checkable against a rule that already exists in the published
  literature, which is what was just done.
- **What this does and does not establish.** It supports the claim that Uri's
  rule is explainable as a composition of an existing catalog entry (rule of
  the square, generalized to include the attacking king) with a color-lock
  fact, arrived at independently by chess analysts and engine developers
  without being framed as boundary-state composition. It does not, by itself,
  establish that composition is a general method for discovering *new*
  boundary states beyond what's already known. One confirmed decomposition of
  an already-known heuristic is evidence the compositional structure is real
  and has already left a footprint in existing chess knowledge. It is not yet
  evidence that the same operation, run forward rather than backward,
  generates anything new — that forward direction is exactly what §8m's
  population scan and §8j's construction attempts are for.
- **Open, and directly checkable without new tooling.** Whether the margin in
  the generalized race (king versus king versus pawn, to the corner) is
  nothing more than the already-published "rule of the square" (king versus
  pawn) with a second racer's distance folded in symmetrically — that is,
  whether Uri's rule requires no structure beyond the pawn-and-two-kings race
  plus the bishop's color-lock. This is a comparison between two already-
  published formulas and needs no new computation, only the comparison itself
  being carried out.
- **This does not conflict with §8c's self-correction, and the two should be
  kept apart.** §8c records that this project's own Method 5 first attempted
  wrong-colored-bishop-and-rook-pawn material as a *construction base* — a
  starting point to compose something new on top of — and found that the
  attacker has too little material there to enforce a tempo-margin argument
  against a detour into a pawn structure. That is a question about whether
  the material can *host further composition*. This section asks a different
  question: whether the material's own already-known draw *itself*
  decomposes into recognizable primitives. §8c's self-correction does not
  bear on the decomposition proposed here and should not be read as
  contradicting it.
- **The "already used in engines" point was already on record and needed no
  revision.** `boundary_board_states.md` §2 Method 4 and §9 already state,
  from the same chessprogramming.org source, that this ending must be
  encoded in the evaluation function because search does not evaluate it
  properly — consistent with, and unaffected by, what's added here.

### 8o. Perturbation search as a candidate-generation method, and a worked case study of what it actually finds

**The idea, as stated, across two turns.** Boundary board states can be
searched for by reverse engineering from known decisive positions: take a
position, and ask whether some parameter — centrally, whose move it is, but
in general a tempo or a defensive resource added or removed — can be varied
until the result crosses into equilibrium (a draw). A growing catalog of
confirmed boundary states would make this calculation easier and help scope
exhaustive coverage. Applied to a specific chess.com daily puzzle: does giving
the side under threat the move instead, so that it can cover the mating
threat and step back from the race, locate a position where the race
conditions balance? When such a perturbation instead flips the result from
one decisive outcome to a different decisive outcome, rather than to a draw,
that indicates the position is not joint compatible at some level; decisive
flipping is itself an outcome of the topological primitive catalog, and once
a genuine equilibrium is found nearby, relativistic analysis between it and
the decisive-flipping positions could inform other such flipping mechanisms
at an attractor-landscape level.

**The worked case study.** The position given (White to move):

`r4rk1/6p1/2p4p/1q1p3Q/p1nP4/P7/KP2B1R1/2b3R1 w - - 0 1`

No mate-in-one exists for White in this position — checked directly, not by
engine score. Stockfish finds a forced mate in four: `Rxg7+ Kh8 Rh7+ Kxh7
Qg6+ Kh8 Qg7#`, confirmed as genuine checkmate by python-chess's own rules
engine after playing the full line out, not taken on the engine's score label
alone. Every other reasonable try for White (bishop moves — Bd1, Bxc4, Bf3)
loses on the order of five to six hundred centipawns, apparently because
Black's queen, knight and bishop are already positioned against White's king
on a2.

The perturbation proposed — give the other side the move — was run directly
on this exact position, nothing else changed. White to move: mate in four for
White, as above. **Black to move: mate in one for Black** (`Qxb2#`), with
`Qb3+ Kb1 Qxb2#` next best. The position's result inverts completely with
whose move it is: White-to-move mates in four, Black-to-move mates in one.

**Assessment — separate from the record above.**

- **This is not a boundary state in either move-order, and it is worth being
  precise about why, because the reason is informative.** A boundary state
  requires the position to be drawn. This position is never drawn — it is
  decisive both ways, with the result flipping between two different forced
  mates rather than settling toward a draw. The perturbation here did not
  approach an equilibrium; it crossed a seam between two wins, with nothing
  in between. That is also not zugzwang, a term it would be easy to reach for
  here: zugzwang is when moving is bad for the side to move, and it would
  rather pass. Here moving is good for whoever does it — the problem is
  purely whose turn it is. It is an ordinary, if unusually sharp, mutual
  attacking race, and the size of the swing (four plies to mate one way, one
  ply the other) is what makes it a striking example rather than an unusual
  category of position.
- **The method itself is sound, and is a real, useful addition — recorded as
  Method 6 in `boundary_board_states.md` §2.** It is Method 2's
  retrograde-seeding instinct (follow value across a change and find the
  seam) applied as a single parametric probe rather than a full backward
  walk, which is exactly what makes it usable outside tablebase-reachable
  material: nothing about running it required the position to be small. Its
  confirmation is whatever confirms the reference position's value — engine
  search above tablebase range — so it holds the same status as Methods 3 and
  4 under the standing rule of §1: it proposes and orders candidates, and
  certifies nothing. A candidate it turns up at any piece count still needs
  the catalog's own certification machinery, or the tablebase, before it
  counts as a confirmed boundary state.
- **On "not joint compatible at some level" — the qualitative reading holds;
  the formal one does not, and the difference matters.** Read loosely, as a
  description rather than a technical claim: nothing here contradicts the
  observation that neither king has a stable defensive structure holding
  against the other, regardless of tempo — both sides have a live, immediate
  mating resource, and nothing locks either one down. That is a reasonable
  informal gloss on what the data shows. But §8g's joint compatibility is a
  specific, technical claim about a *set of individually identified and
  certified small primitives* composed together, checked by a resource
  accounting across that set. Nothing has been decomposed into primitives
  here — this is one holistic engine evaluation of a nineteen-piece position,
  with no primitive catalog entries identified, let alone certified, for
  either king's situation. Calling this an instance of joint incompatibility
  in the formal sense would be the same category move already corrected
  twice in this conversation — reusing a precisely defined term for something
  that resembles it in spirit but hasn't met its conditions. The informal
  reading is fine and worth keeping; the formal one isn't yet earned.
- **On "eventually finding an equilibrium position... at an attractor
  landscape level" — recorded as a genuine, open, currently unactionable
  research question, not a result.** The proposal is to compare the basin
  structure around confirmed boundary states to the basin structure around
  decisive-to-decisive flip points such as this one, using the first to
  interpret the second. That is a coherent thing to ask. It is also several
  steps past anything this project currently has: no population of confirmed
  boundary states exists yet to compare against (§8m's own open question),
  and nothing about "basin structure" has been defined precisely enough to
  measure. Recorded as a direction, not a method.
- **What the worked case study is actually worth.** Not a boundary state, and
  not evidence of joint incompatibility in the technical sense. It is a
  clean, verified demonstration that Method 6 can land on a different real
  structure than the one being searched for — a decisive-to-decisive tempo
  seam — and that this is itself worth recording as a known outcome of
  running the method, not a failure of it. The next genuine test of Method 6
  is to run it where it's more likely to land on an actual draw: starting
  from a position already close to known drawing material (the small
  materials targeted in `boundary_board_states.md` step 13) rather than from
  an attacking middlegame position where both kings are already exposed.

### 8p. Relativized classification — using a catalog of boundary states as a network of reference points, not only as leaf-node termination targets

**The idea, as stated.** Beyond using a single known boundary state as a
termination point for forward search, consider a catalog of several. Compare
them to each other relativistically, and to an unknown position, under
perturbation: can the unknown position, under optimal play, be navigated
toward one of the known equilibria? If it can, its value is informed by that
equilibrium's own certainty. If it cannot — if every attempt to perturb the
unknown position toward any known boundary state instead resolves to one
side's advantage, no matter how the perturbation is tried — this itself is a
way of establishing that the unknown position carries a definite advantage,
derived from its relativistic position with respect to a population of known,
certain equilibrium points. This would let the catalog scope exhaustive
coverage: not by classifying every position, but by using known certainty as
reference points from which nearby unknown positions can be classified.

**Assessment — separate from the record above.**

- **This already exists, worked out exactly, as classical endgame theory —
  checked this session, not assumed.** "Key squares" are exactly this
  project's certainty anchors: reach one, and the result is known. "Corresponding
  squares" is the generalization this idea is describing, already named and
  already in standard use: a network of paired squares such that if one
  king occupies one member of a pair, the opposing king must occupy the
  other to hold the position, extended well past simple opposition
  specifically because a single opposition calculation stops working once a
  position gets "cluttered" with more structure. That is a precise structural
  match: key squares are the catalog of certain anchors; corresponding
  squares are the relativistic network of relationships between them and
  everything nearby; the parity/tempo calculation that determines whether two
  squares correspond is exactly the "perturbation under optimal play"
  described here, already made exact for king-and-pawn structures.
- **The one place this needs a guard rail, and classical theory already states
  it explicitly, independent of anything in this project.** A directly
  relevant source checked this session puts it plainly: reaching a key square
  is what wins; if the key squares are already unreachable, "you can hold the
  opposition all day and achieve nothing." That is the precise failure mode
  to guard against here: **failing to find a path from the unknown position to
  any known boundary state in the catalog does not, by itself, establish that
  the position is decisive.** It establishes that this catalog, as it
  currently stands, doesn't classify the position — which could mean the
  position is genuinely decisive, or could mean the position is itself a
  different kind of drawn equilibrium the catalog doesn't yet contain, or
  could mean the perturbations tried weren't the right ones. Absence of a
  reachable known anchor is absence of information, not presence of a
  result, and treating it as the latter is exactly the inference the standing
  rule of §1 has guarded against throughout this document: a heuristic
  search's failure to find something is not a certificate that the something
  doesn't exist.
- **The defensible version, and it's a strong one, is the version corresponding-
  squares theory actually realized.** The catalog becomes certainty-generating
  for a nearby position only when the correspondence relationship itself is
  exact and provably complete for the structure in question — not when a
  search merely fails to find a connecting path. Corresponding-squares theory
  achieves this by deriving the full network for a given blocked pawn
  structure exhaustively, the same bar every certified method in this
  catalog already has to clear (Methods 5 and 6 both explicitly scoped this
  way). The source checked this session states the actual use pattern
  precisely: key squares tell you where to go, the opposition tells you how
  to get there — scoping *where to focus the exact calculation*, not
  replacing the calculation with the absence of a found path.
- **Stated as affirmatively as the guard rail above was stated cautiously,
  because both are true at once and neither should read as contradicting the
  other: yes, this determines decisiveness, and it is not aspirational — it
  is exactly how real K+P endgame theory already works, proven, in daily
  use.** A king at some arbitrary square, not itself a key square, is
  routinely and correctly classified as won, drawn, or lost purely by parity
  counting against the known key-square network — no search of that specific
  position's full game tree required, because the offset itself (whose move
  it is, and the exact tempo distance to the nearest relevant key square) is
  what the classification runs on. That is a real, provable instance of
  exactly the mechanism described here: the race-condition primitive in an
  unequal position is the same primitive that sits in equilibrium at the
  boundary, and the measured offset between them is what determines the
  result. This is not in dispute, and the guard rail above was never meant to
  deny it — it was narrower than that, and worth restating precisely now that
  the positive case has been laid out in full: the guard rail applies when
  the offset is *estimated* (an engine's centipawn score, or a search that
  simply failed to find a connecting line) rather than *exact* (a closed-form
  parity count against a provably complete network), and when the catalog
  being checked against is not itself proven complete for the structure in
  question. K+P key squares clear both bars — the Chebyshev-distance
  calculation is exact and closed-form, and the key-square network for that
  material is exhaustively derived and complete. Neither bar is yet cleared
  by anything this project has built for the broader primitive types it is
  trying to catalog: every offset measured so far in this conversation
  (the puzzle's centipawn swings, the opening's centipawn field) has been
  Stockfish's estimate, not a closed-form count, and §8m's own population
  scan — whether primitives even recur at usable scale — is still unrun. The
  method is correct and already proven correct, by existing chess theory, for
  the materials where it has been made exact. Whether this project's own
  catalog can be brought to the same standard for its broader primitives is
  the actual open question underneath everything from §8f onward, not
  whether the underlying logic works — it does.
- **This is a real, different operation from Method 6, and is worth recording
  as its own direction rather than folded into it.** Method 6 perturbs a
  single position's parameters and watches where the value crosses a
  boundary. What's proposed here uses a whole catalog as a network of
  reference points against which other positions are classified by their
  relationship to the network, the way corresponding squares are used against
  an actual position on the board — structurally closer to the scoping
  mechanism already described in `boundary_board_states.md` §10 (a validated
  boundary catalog as an early-termination oracle) than to Method 6, but
  generalized from "the search literally reaches a cataloged state" to "the
  search's relationship to several cataloged states can be computed even
  short of reaching one." Cross-referenced there.
- **What corresponding-squares theory does not resolve, and this project's
  own scope question returns exactly as before.** It's exact and complete for
  king-and-pawn structures, a small and highly constrained class. Whether an
  analogous exact, complete correspondence network can be derived for the
  broader, more varied primitives this project's catalog is trying to build
  (duty-based locks, not just frozen walls; multiple piece types) is exactly
  the open question §8m already poses about whether primitives recur at
  useful scale — this section sharpens what "useful scale" would need to mean
  if the catalog is to be used as a classification network rather than only
  as termination points.
- **Named synthesis, proposed afterward: boundary-state families scoped by
  opening line.** `boundary_board_states.md` §10 frames the catalog as a
  termination oracle for forward search from one specific starting position.
  This section frames it as a network usable for relativized classification.
  Put together, as proposed: group the catalog's entries into families by
  which opening lines can reach them under optimal perturbation (§8o, §8p),
  so that forward search from a given opening isn't checked against the
  catalog one entry at a time, but against the family its line is already
  known to be able to reach — narrowing the relevant part of the network
  before the search even starts. This is a real tightening of what the two
  sections already said separately, not new machinery on top of them, and it
  inherits every open question both sections already carry: §8m's
  recurrence question for what populates a family at all, and this section's
  own guard rail that failing to reach a family doesn't establish a
  position's value. It also inherits the scale wall already named throughout
  this document (§8j, §13's standing caveat in `boundary_board_states.md`):
  classifying which families are reachable from which openings is, at full
  scale, the same size of problem as classifying large parts of chess
  itself, and naming the synthesis doesn't shrink it.

### 8q. Two corrections from one challenge: composition needs start-reachability, and a refutation proof is not a boundary-state certificate

*Correction (`boundary_board_states.md` §14 item 8): this section uses "exactly one non-losing move" as if it were the definition of a boundary state. It is a measure of width (w = 1 is the strict case), not the definition. Where a conclusion below depends on w = 1 specifically, it holds only for that strict case.*

*Withdrawn: the argument, made in discussion, that the opening position is not a boundary state because it has many non-losing first moves. That rested on the strict reading. The opening position is not ruled out by width, and whether it is an equilibrium of the kind defined in item 8 is not established either way. The distinction between a refutation proof and a certificate stands, and so does the finding that the position after Qf3 in the OTB game is decisive and not drawn.*

**The idea, as stated.** Two primitives are either jointly compatible or they
are not — a binary fact about a position, not a question of degree. But
composition as discussed in §§8g-8p has been treated as if primitives and
their geometry could emerge independent of where they come from, when in
reality primitives do not arise from pieces placed on a board with no regard
for how they got there — they stem from the initial board position's own
equilibrium, through the race toward development that every opening engages
in. The asymmetry that produces a result like early queen development's cost
emerges specifically because the position's race-condition requirements were
not accounted for relative to this origin. If a forcing sequence can be found
after such a move, with no response available that both saves the piece and
keeps pace with the opponent's own development, that single counter-example
is sufficient to establish the move's suboptimality — not because the whole
tree was exhaustively searched, but because the counter-example itself proves
no compensating line exists.

**Assessment — separate from the record above, in two parts, because the
message makes two different claims that need different answers.**

- **First claim: composition needs to be checked for start-reachability, not
  only internal validity. This is correct, and it closes a real gap.**
  Nothing in §§8g-8p has required that a composed primitive set actually be
  reachable via legal play from the game's own starting position — only that
  it be internally valid once placed. For the catalog's stated purpose
  (scoping real forward search from real games, §7a, §10, §8p), that gap
  matters: an unreachable but internally valid composition is not useless —
  it is still a legitimate catalog entry showing a structure *can* exist as a
  boundary state, which is real information about the primitive itself — but
  it is useless for the specific application of terminating search from an
  actual opening. This was not missing from the project's own framework,
  only from how Method 5 and §§8g-8p had been connecting to it:
  `boundary_board_states.md` §14 item 1 already defines exactly this
  distinction as the tier structure, with Tier 3 named "start-reachable" for
  precisely this reason. The correction is to make explicit, here, that any
  composition work aimed at the forward-search-scoping purpose needs its
  targets checked against Tier 3 specifically, not left at Tier 1's
  internal-validity bar. Recorded as a standing requirement on Method 5 and
  on §8j's build order from this point forward.
- **Second claim: one forcing counter-example suffices to prove a move
  suboptimal, without full-tree search. Also correct, and not new — it is
  the same proof shape §8h already built for a different purpose, which is
  worth crediting as a genuine, independently-arrived-at match.** Refuting a
  move by exhibiting one line the opponent cannot escape, rather than
  searching every continuation, is standard minimax refutation analysis —
  and it is structurally identical to §8h's interference closure: check
  everything in the relevant footprint, not everything on the board. The
  footprint there is defined by tempo-bounded reachability, line closure, and
  forcing-move entry points; the same three conditions govern what a
  refutation counter-example has to rule out to be genuinely inescapable.
  Arriving at this independently, from a completely different direction
  (opening tempo rather than endgame primitives), is real confirming evidence
  that the underlying proof technique is doing honest work rather than being
  an artifact of the endgame framing it was first built in.
- **But refutation and boundary-state certification prove different
  properties, and the difference is not a technicality — it is the same
  distinction that has come up twice already in this conversation (the
  opening position itself, and the puzzle's mutual-mate-flip), and is worth
  settling here rather than re-deriving a fourth time.** A refutation proof
  shows a move is dominated — some alternative beats it — and this is true of
  moves in decisive positions exactly as often as in drawn ones; dominance
  has nothing to do with whether the position's own value is a draw. A
  boundary-state certificate shows something narrower and more specific: the
  position's value *is* drawn, and exactly one move preserves that draw while
  every other move loses. These are different claims about different
  objects, provable by overlapping but not identical means, and a clean
  forcing line proves the first regardless of whether the second holds
  anywhere nearby.
- **Checked against the actual game, not argued in the abstract.** Before
  White's Qf3, Stockfish's top five candidates span 4 to −60 centipawns — an
  ordinary field, nothing forced. After Qf3, Black's top five span +7 to
  +166 — clearly better for Black, confirmed and real — but that is a field
  of improving options, not one forced line from which White has zero
  escape. White has continuations from there; they start from a worse
  position, which is exactly what a dominance/refutation claim establishes,
  and exactly what "no moves downstream, no compensation possible" would need
  to be stronger than to match the language used. The position only
  approaches the "no escape, forced" texture described much later, after
  several further moves compound the disadvantage (consistent with the
  earlier finding, §8o, that this kind of razor-thin margin is a real but
  distinct phenomenon, found at specific later junctures — Kd1 vs. Kf1, for
  instance — not as a property of the queen move's immediate consequence by
  itself).
- **What this leaves standing, precisely.** The queen-development read is
  real chess judgment, correctly identified, and provable by ordinary
  refutation technique without needing a whole game tree — all of that
  holds. It establishes that Qf3 was dominated by better alternatives. It
  does not establish, and the data above shows it did not in fact produce,
  an inescapable forced sequence at that specific move, and it was never a
  claim about the position being a drawn equilibrium in the first place.

### 8r. Forward anchoring, the discovery/certification split, and cooperative-mate relaxation

*Correction (§8y): the statement below that discovery only proposes applies to discovery by population features and by relaxations. Discovery by residual, meaning a proven verdict that the equalized catalog does not account for, is a deduction under the definition, and its output is exact.*

*Correction (`boundary_board_states.md` §14 item 8): this section uses "exactly one non-losing move" as if it were the definition of a boundary state. It is a measure of width (w = 1 is the strict case), not the definition. Where a conclusion below depends on w = 1 specifically, it holds only for that strict case.*

**The idea, as stated.** The initial position informs every downstream
mechanism up to the mates Syzygy encapsulates, so cataloging primitives and
boundary states should begin from the forward trace and treat it as more
important than backward tracing; a family of reachable boundary states is
properly anchored only when tied to the start position. Within that forward
search, exhaustive Syzygy-style coverage is not the focal mechanism. The focus
is perturbation of race conditions to mate, and the balancing and mitigation
of those races in a way that is joint compatible with the others and that
explains why the primitives emerge at all. Stated precisely: *confirming* a
boundary state (a draw under perfect play) depends on exhaustive coverage,
while *discovering* race-condition primitives does not, and can proceed by
perturbation and by finding the quickest routes to mate and asking which race
condition is driven toward zero. The thought experiment: take the fastest
route to mate with no optimality considered, as if both sides cooperated to
end the game as soon as possible. The structure that emerges (the loser's
pieces locked out, a pawn move opening a route for the queen, the delivery
beating any "too early" criterion) is a specific race-condition primitive
under joint compatibility. With a chosen opening line, each piece move's
effect on those primitives, and its direction relative to a boundary state,
can be tracked; enumerating downstream from sub-optimal play yields
candidate, proven and disproven catalog entries, paired with backward search
and a growing catalog, and the segmented race-condition understanding makes
perturbation more efficient. The structure reinforces itself.

**Assessment — separate from the record above.**

- **The discovery/certification split is correct, and it states cleanly
  what §1 has held from the start.** Primitive discovery may use any
  heuristic, perturbation, or relaxation, because discovery only proposes.
  Certification of a boundary state needs the complete case list of §8d,
  which is exhaustive within the position's footprint (§8h), not necessarily
  over the whole tree. Much of the earlier friction in this record was about
  which of those two a given argument was doing; with the split named, most
  of it dissolves.
- **Forward versus backward is a question of what each direction supplies,
  not which is more important.** The forward trace supplies *relevance*:
  which structures are reachable from the start position (Tier 3,
  `boundary_board_states.md` §14) and therefore worth cataloging first. The
  backward, exact side supplies *certified anchors*: the tablebase, closed-form
  offsets, and completed case lists. Forward trace alone cannot certify
  anything, and backward tracing alone cannot say which anchors matter. The
  checkers weak solution used both, a forward search reconnecting to a
  retrograde database. As a priority ordering for *selecting what to
  catalog*, forward-first is defensible and is adopted here; as an ordering
  for *certifying*, nothing changes. This should be reflected in the build
  order of `boundary_board_states.md` §11: its current steps are
  retrograde-first because exact anchors are cheapest there, and a
  forward-driven selection layer belongs on top of them. One honest limit: a
  forward trace from the start position over the whole tree is the huge
  search, and perturbation exists to avoid it. The only perturbation
  procedure defined so far is Method 6's single-parameter probe, and nothing
  yet describes a forward perturbation search over many plies.
- **The cooperative-mate thought experiment is a known and well-behaved
  relaxation, and it is useful as a primitive-discovery tool.** It is the
  helpmate / shortest-game idea from chess composition, and it is a standard
  relaxed-problem device in heuristic search. Let CMD be the cooperative
  distance to mate (both sides choose moves to minimize it) and DTM the
  adversarial distance (the defender maximizes delay). Then CMD ≤ DTM
  whenever a forced mate exists, so CMD is a lower bound, and the gap DTM −
  CMD is what the defender's best resistance adds. That gap framing is this
  assessment's own observation, not something from the record above, and is
  untested.
- **A worked case, enumerated against the real start position rather than
  argued.** Exhaustive enumeration from the initial position finds no mate in
  1-3 plies by either side. Black's fastest cooperative mate takes 4 plies, in
  exactly 8 sequences, every one a fool's-mate variant (White plays f3 or f4
  and g4, Black plays e6 or e5 and then Qh4#). White's fastest takes 5 plies,
  in 347 sequences, and the sample listed all contain one spare move (g1h3 or
  b1c3, say) plus the same three ingredients mirrored: Black weakens the king
  diagonal with f- and g-pawn moves, White moves the e-pawn, and Qh5 delivers.
  The structure is the same for both sides: two moves by the loser that
  weaken the king's diagonal, one line-opening move, and one delivery, four
  plies of structure. The asymmetry in distance (4 versus 5) appears to be
  ply parity: White's mating move must land on an odd ply, so the first
  mover's tempo shows up as one spare ply and not as a faster mate. (The
  parity explanation is this assessment's reading of the enumeration, not
  separately verified.) This matches the stated structure of a lock (the
  loser's weakened diagonal), a route (the opened line), and a delivery, and
  it gives the "locked out pieces" and "not early because it mates first"
  language a concrete minimal instance.
- **What the relaxation cannot supply.** Cooperation removes the defender's
  choice, and the defender's choice is exactly where boundary states live: a
  boundary state is one where the side to move has exactly one non-losing
  move. The cooperative game therefore yields the attacker-side primitives
  (what has to be routed, locked and opened) and says nothing about
  defender-side primitives (what the defender does to break the route). The
  composition needs both, and the adversarial half still comes from
  perturbation and verification, not from the relaxation. A possible formal
  handle, also this assessment's own and untested: at a boundary state CMD
  is small and finite for the side that threatens, while the adversarial
  value is a draw, so the position sits where a near mate is held off by one
  move.
- **"Early" is a relation to a race, not a property of a move.** This agrees
  with §8q: whether a queen move is early is decided case by case by whether
  it wins or loses the race, which makes it a dominance judgment (§8q) and not
  a universal law. That does not weaken it as a heuristic.
- **"Distance to a boundary state" needs a definition before it can be a
  per-move measure.** As stated, each move's effect on the primitives and its
  direction toward equilibrium would be tracked along an opening line. No
  computable definition exists yet for positions outside the catalog.
  Proposed operationalization, with its prediction stated before it is run so
  that it counts as a test in the sense of `boundary_board_states.md` §15: for
  each ply of the OTB game in `boundary_board_states_case_study_otb.md`,
  compute a small set of measurable race quantities (a depth-limited CMD for
  each side, the number of undeveloped pieces, the queen's count of safe
  squares), and compare their per-ply change against Stockfish's evaluation
  change. The prediction is that at the identified junctures (Qf3, Nb4, f5)
  the quantities move in the direction the evaluation does. If they do not,
  "race progress" as framed does not track value. Not run. The cost is real:
  exact CMD at mid-game positions is a deep search, so it would be depth-limited.
- **The loop is real, and it is coherence, not corroboration.** The pieces do
  support one another: forward reachability selects which primitives to
  catalog, the relaxation and perturbation propose them, exact methods certify
  them, the catalog scopes later forward search, and the families tie it back
  to the start position. A structure whose parts justify each other is
  internally consistent; it is not thereby externally confirmed. What keeps it
  from being circular is the set of outside anchors it can be checked against:
  tablebases, closed-form offsets, published endgame theory, and disagreement
  tests against an engine.
- **Rediscovery of known theory as a validation suite, by design.** Known
  results are ground truth for the architecture: a primitive decomposition
  that fails to reproduce a known result exactly is wrong somewhere. Candidate
  list, each to be given a decomposition and a pass criterion before it is
  run: the rule of the square, opposition and key squares, corresponding
  squares, Uri's rule (partly done, §8n), the Philidor and Lucena positions,
  the Vančura position, and triangulation. Passing is meaningful only where
  the offset is exact (§8p); opening principles such as tempo count as
  heuristic rediscovery until an exact version is derived.

### 8s. The fool's mate as a minimal primitive: ablation results, the four-term template, and what is new in it

**The idea, as stated.** The fool's mate is informative because it is the
earliest cooperative mate and needs the least forward trace to confirm. It
is not optimal play, but it shows which race-condition primitives arise from
the start position and why, and a population of such results anchored at the
start position would separate real primitives from noise. The primitive is
not the specific moves but the piece-locking structure: the mate depends on
undeveloped, locked pieces and on a diagonal opened by the loser's own pawn
moves, and it would still be valid with the knights and rooks removed; the
structure is not piece dependent apart from the king. The same structure
recurs elsewhere: a castled king's coverage, a diagonal opened by an advanced
pawn, a pinned piece that cannot block. In the fool's mate nothing needed
pinning because the defender had cooperated; where a blocker exists, pinning
it removes its coverage and does the same work, and so does a decoy, a piece
placed so that capturing it costs the coverage that was holding the position.
A boundary state can then emerge at the most upstream point of a commitment
that forces the opponent into a position where they can only lose or draw,
the race being equalized on the branch taken.

**The data, run rather than argued.** For each of three real mates, every
subset of the surrounding pieces was tested for whether the mate survives:

| Position | Necessary | Free (never changes the verdict) | Conditional |
|---|---|---|---|
| Fool's mate, White's pieces (15 varied) | Qd1, Bf1, Pd2, Pe2 | Ra1, Nb1, Bc1, Ng1, Pa2, Pb2, Pc2, Pf3, Pg4 | Rh1 with Ph2 |
| OTB game, Qxh7# threat, Black's pieces (13 varied) | Rf8 | Nb4, Bc5, Pe5, Pd6, Pa7, Pb7, Pc7, Ra8, Bc8, Qd8 | Pf7 with Pg7 |
| Puzzle, final Qg7#, both sides (14 varied) | Rg1 | every other piece | none |

Further checks on the fool's mate: with all four knights and rooks removed it
is still checkmate; with the queen on h4 replaced by a bishop it is still
checkmate; replaced by a rook, knight or pawn it is not; a white pawn put
back on g3 or f2 blocks the ray and there is no check. The h-file pair
behaves as a unit: with Ph2 removed and Rh1 present White answers Rxh4, but
removing both, or the rook alone, leaves the mate intact. In the OTB position
the f7/g7 pair behaves the same way, since g7 shields f7 from the queen's
rank attack, so f7 is needed as an own-piece block only while g7 stands. And
a pin version of the pattern was built: White knight on e4 (able to interpose
on g3 or f2) with a black rook on e8 pinning it to the king and a black bishop
on c4 covering e2 is checkmate; the same position without the rook is not (Ng3
and Nf2 answer), and with the knight absent it is checkmate again.

**Assessment — separate from the record above.**

- **The stated predictions held.** The mate survives removing the knights and
  rooks, and the pieces it cannot do without are exactly the king's four
  start-square neighbors (Qd1, Bf1, Pd2, Pe2), so "dependent on the initial
  start locking" is literally true: none of them can have moved. This was
  stated before it was checked, and it could have failed. It is a small
  corroboration in the sense of `boundary_board_states.md` §15, modest because
  the claim is an easy one.
- **"Not piece dependent apart from the king" holds, with the exact form of
  the dependence stated.** What is king-relative is the set of conditions;
  what fills each condition varies. The attacker on the ray can be a queen or
  a bishop and must be a piece that attacks along that line. The conditions
  are exactly the rules' own definition of checkmate: a check is delivered, no
  flight square is available, the checking piece cannot be captured, and the
  check cannot be blocked. So the template of the fool's mate is the
  definition of checkmate decomposed, and it is not a discovery in itself.
  The content of the primitive library is the *mechanisms* that satisfy each
  term and what they cost in tempo. Observed so far:

  | Term | Fool's mate | OTB Qxh7# threat | Puzzle Qg7# |
  |---|---|---|---|
  | Check | open diagonal h4–e1 (g3, f2 emptied by the loser's pawns) | open diagonal e4–h7 (f5, g6 empty) | adjacent queen, no ray |
  | Flights denied | own pieces on d1, d2, e2, f1, plus the ray covers f2 | own rook on f8, the f7/g7 pair, queen covers h8 | queen covers g8 and h7, no own pieces needed |
  | Checker not capturable | no white piece reaches h4 (the h1 rook masked by its own h2 pawn) | Bd3 supports h7 (held fixed in the test) | Rg1 supports g7 (necessary) |
  | Block | none possible | not applicable | not applicable |

  Mechanisms in that table: own-piece lock, attacker coverage, support,
  masking of a defender by its own piece, and defenders that are simply out of
  play (the "free" column: Black's queen, bishop and second rook in the OTB
  position were irrelevant to the mate).
- **The pin generalization is correct for absolute pins and needs a cost
  comparison for others.** An absolutely pinned piece cannot legally leave the
  pin line, so its coverage is removed exactly as if it were absent, which the
  pin construction above shows. A relative pin (to a piece other than the king)
  does not remove the coverage: the pinned piece can still move, at a price in
  material, so it counts as absent only if using it costs more than the
  position is worth. This is the duty picture of §8f: a pin is one way to make
  a coverage duty unperformable. A decoy is another; the puzzle's mate is built
  from them.
- **Forced enabling moves are the bridge to adversarial play, and are this
  assessment's own synthesis.** In the fool's mate the loser's enabling moves
  (the f- and g-pawn advances) are voluntary, which is why it is a cooperative
  mate. In an adversarial game each such move has to be forced, and the
  mechanisms for forcing it are the ones named: a check, a decoy sacrifice, a
  pin, a capture threat that cannot be answered while holding the coverage.
  The puzzle's mate is the forced version of the same skeleton: Rxg7+, Rh7+
  and Qg6+ leave Black only Kh8, Kxh7 and Kh8, and Qg7# follows. The OTB
  game's Qxh7# threat is the same skeleton waiting, with f5 closing the ray.
- **Single-piece ablation is not enough; the dependencies are joint.** The
  Rh1/Ph2 and Pf7/Pg7 pairs are small real examples of the point of §8g: a
  piece can be harmless alone, harmful with a partner, and harmless again when
  both are gone. The load-bearing test has to be run over sets of pieces, as it
  was here (all 32,768 subsets in the largest case), not one piece at a time.
- **The classification into necessary, free and conditional is the noise
  filter for a population.** The 347 White fool's-mate-length sequences found
  in §8r differ mostly in a spare move; the ablation procedure is what collapses
  them to one skeleton. It has been run on three final positions, and not yet
  on the 347 sequences.
- **Development discharges terms, by geometry and not run.** A knight on f3
  attacks h4 and also occupies f3, which blocks the f-pawn from opening the
  ray. One developing move therefore discharges two terms of the fool's mate.
  This is a concrete form of the claim that piece activation unlocks
  defenders, and it is stated here from the knight's move pattern, not
  tested.
- **A terminology point for the catalog.** The pawns doing the work are the
  king's shield pawns advancing and weakening the diagonal. A passed pawn is a
  different thing (no opposing pawns ahead or beside it), and using the term
  here would blur a standard concept.
- **On boundary states emerging at the most upstream commitment.** Read from the committing side, a commitment is a move after which the opponent is in an equilibrium: some of its moves lose, and the draw holds only through the equilibrium moves. "Can only lose or draw" is true of every drawn position for the side to move, so it is not what distinguishes a boundary state; what does is that the opponent has a losing deviation available while the draw is sustained by a cycle (`boundary_board_states.md` §14 item 8). An earlier version of this bullet said that the count of drawing replies being exactly one distinguishes it. That was the strict reading and is withdrawn.
- **What none of this touches.** All three positions are decisive mates. The
  template gives the attacker-side half. A boundary state needs the defender-
  side half: what the single saving move preserves, and which term of the mate
  it keeps unsatisfied. That is where the draw lives, and it is the next
  experiment.

### 8t. The mating-net predicate, irreversibility, and the first exact scan

**The idea, as stated.** The fool's mate should not be read piece by piece but
as a composition of occupancy, coverage and the threat of placement. If a pin
is added to lock a piece, it is the same primitive: pieces locked, by
development or by a pin, and a king in a fortress facing a mating threat whose
delivery route is uncovered because the race to cover it is lost. It is not
piece specific; only the king was named. Pawns matter because a pawn cannot
move backwards, so an irreversible pawn move can open the opportunity, and the
side that chose the enabling move is the one exploited; where no pawns are
involved the fortress mechanisms are non-pawn ones and more trivial.
Perturbing a pawn move can reveal when the race has been altered, and an
equilibrium has to be jointly compatible with this primitive. Removing a piece
that was covering the route is another way to break it. Reversible moves that
balance the equilibrium are a requirement of boundary states, and an
irreversible move that breaks the coverage with no way to re-cover is a
decisive break. Finally, enumerating the population that is meant to be
cataloged is useless before the structure exists.

**Assessment — separate from the record above.**

- **Correction, added afterwards (`boundary_board_states.md` §14 item 8).** This
  section counted drawn positions with exactly one saving move and called them
  boundary-like. Exactly-one is a measure of width, not the definition, and the
  positions counted are not boundary states: the stronger side in this material cannot
  lose, and the family includes captures of a hanging pawn into a dead draw, which are
  blunders by the pawn's side and not an equilibrium. What the scan did measure is the
  width of drawn positions by pawn file, and that finding stands: the rook-pawn draw is
  the widest, which counts against exactly-one as the definition, since the rook-pawn
  draw is the project's paradigm local primitive. The figures below are retained as a
  width measurement.
- **Scope of the scan in this section.** The king-and-pawn-against-king scan
  below answers only the side remarks about irreversibility and about reversible
  moves, and the project's own earlier edge guess. It says nothing about the
  fortress and mating-net primitive, and a forced choice of material is the wrong
  instrument for finding that primitive. It is not evidence for or against it.
  The test of the primitive is §8u.
- **The primitive is best stated as a predicate over occupancy and coverage,
  with no piece names, and this is the form to adopt.** Let k be the king's
  square, D the defender and A the attacker. Let Occ_D be the squares held by
  D's pieces, Cov_A the squares A attacks, and Cov_D the squares D's pieces
  can *legally* move or capture to, where a pinned piece's coverage is
  restricted to its pin line. A mating net on k holds when four conditions are
  all true: (C) A gives check along an open line or from an adjacent square;
  (F) every flight square of k lies in Occ_D or Cov_A; (X) every checking piece
  lies outside Cov_D; (B) every interposition square on the check line lies
  outside Cov_D. The locking mechanisms are then exactly the ways of removing
  squares from Cov_D (an own-piece block, i.e. undeveloped pieces; a pin; a
  deflection or decoy; overload) or adding squares to Occ_D or Cov_A. The race
  is the tempo needed to realize each. This is the §8s table restated without
  piece names, and piece names belong only in the catalog of instances. It is
  the definition of checkmate, as §8s noted, so the content lies in the
  mechanisms and their tempo costs.
- **Structure before population is the right order, and the §8s item about
  running the ablation filter on the 347 sequences is demoted.** A population
  is for validating a stated predicate, not for finding it. The ablation in
  §8s was over sets of pieces, and it found joint dependencies, but its
  results are most useful as instances of the predicate above.
- **Irreversibility is a real structural fact, and it is wider than pawns.**
  Pawn moves, captures, and the loss of castling rights cannot be undone. This
  is why tablebase construction is organized in layers: a position after such a
  move cannot recur, so the value of a position depends only on positions in the
  same layer and in layers reached by an irreversible move. Within a layer, a
  draw is sustained by reversible moves.
- **Whether an irreversible break is decisive depends on the tempo budget, not
  on irreversibility itself.** Checked in the fool's mate: after 1.f3 e5 White
  has 19 legal moves and exactly one of them (g4) allows mate in one; after
  1.g4 e5 two of 21 do (f3 and f4); after 1.f3 e5 2.Nh3 none does, because the
  knight on h3 covers f2. The first irreversible weakening is not decisive; the
  second is, and a developing move that restores coverage undoes the damage.
  This fits "no viable way to cover the change" as a statement about how many
  tempi remain, and it shows the hazard is narrow: one losing move among 19.
- **A first exact scan, run on king and pawn against king.** The position
  graph was solved by retrograde iteration (the pawn never changes file, so
  files a to d suffice by symmetry), with promotions to all four pieces
  evaluated. It gives 62,480 won positions with White to move and 48,802 with
  Black to move, 111,282 in total, which matches the count I recall for the
  standard bitbase (from memory, not checked against a source). It matched
  Stockfish on all 600 randomly sampled legal positions, 405 won and 195 drawn.
  The script is `kpk_boundary_scan.py`. These are Tier 1, local positions only:
  Black's bare king fails the materially-live test of §5a-iii, and
  reachability from the start was not checked.
- **What the scan found.** Of 35,201 drawn Black-to-move positions with at
  least one legal move, 5,125 (14.6%) have exactly one non-losing move. Of
  those, 26 are positions where Black had only one legal move, so no choice
  existed. The other 5,099 are genuine choices, and in 3,960 of them the unique
  saving move is a king move, while in 1,139 it is capturing the pawn.
- **The edge hypothesis of §8k is falsified, and this changes the picture of
  the rook-pawn draw.** The fraction of drawn positions with exactly one saving
  king move is 7.1% on the a-file, 9.4% on b, 13.3% on c and 15.0% on d, and the
  mean number of non-losing moves falls from 4.15 on the a-file to 3.43 on the
  d-file. The defending king stands on an edge square in 36.3% of the
  exactly-one positions against 33.6% of drawn positions generally. So the
  rook-pawn draw is the *widest* draw, the most forgiving, and the narrow
  draws are the central-file opposition positions. This does not make the
  rook-pawn draw less real as a primitive; it makes it a wide edge primitive,
  useful as terrain (the lifting idea of §8k), and it means the first
  population of narrow boundary states is the central-file one, not the
  rook-pawn one. A hypothesis stated in advance was tested and failed, and it
  is recorded as a failure in the sense of `boundary_board_states.md` §15.
- **Whether reversible moves are a requirement of boundary states.** Not as a
  logical requirement of the definition: 1,139 of the 5,099 genuine-choice
  positions (22%) have an irreversible capture as the unique saving move, though
  these are the hanging-pawn family, not a locking mechanism. For the 3,960
  king-move positions the claim holds, but trivially, because Black has no other
  piece to move, so king and pawn against king cannot test it. As a design
  principle for composed positions it is sound and matches §8a, where frozen
  pawn chains remove the irreversible options and leave reversible moves to hold
  the balance. Testing it properly needs a material in which the defender also
  has pawn moves.
- **Whether irreversible pawn moves are the typical way a position is lost
  is not supported in king and pawn against king.** In won positions with
  White to move, a king move throws the win away 21.6% of the time it is
  available (82,298 of 381,094), a pawn push 15.2% (8,662 of 57,095). Pushes
  are 13.0% of the available moves and 9.5% of the spoiling ones. A push can
  still be the decisive error, as g4 is in the fool's mate, and two instances
  do not settle a general rule in either direction.
- **The decisiveness of any break is a property of the position's value, not
  of the move's irreversibility.** A spoiling king move is not physically
  irreversible, yet the win is gone, because the resulting position is drawn.
  Irreversibility shapes the structure of the game graph; value is what makes
  a break decisive.
- **What now exists that did not.** A population of 3,960 narrow, local,
  Tier 1 positions in one small material, which is the thing §8m's open
  question asked for. They are not independent (many are neighbors of one
  another), so the number overstates how many distinct primitives there are.
  Decomposing them by ablation to see what recurs is now possible and has not
  been done. **Correction:** under the corrected definition (`boundary_board_states.md`
  §14 item 8) these are a width benchmark of drawn positions and not a population of
  boundary states, so §8m's question is still unanswered.

### 8u. A second fortress mate in an unrelated position: the predicate checked, and the net assembled by force

**The idea, as stated.** The primitive does not depend on the pieces on the
board. It is found by starting from the initial position, taking the first
primitive that arises, and checking it against other positions where the
fortress occupancy and pinning structure emerges from race conditions. Forcing
specific piece combinations, as in a tablebase scan, gives mostly noise for this
purpose. As a test, a daily puzzle from the chess.com app, Black to move, with
White: Kh2, Qd5, Rg2, Bf2, Nd4, pawns a2, b3, e3; Black: Kh7, Qb1, Rf6, Nh5,
pawns g7, g5, f5. The solution as given is Nf4, exd4, Rh6+, Kg3, Qd3, Qf3, g4,
Qxd3, Rh3#. Black sacrifices the queen to position into the fortress, and the
point is the primitive at work in a completely different position, where the
pieces matter less than the emergent structure of the fortress built by force.

**The position and the line, checked.** FEN
`8/6pk/5r2/3Q1ppn/3N4/1P2P3/P4BRK/1q6 b - - 0 1`. As typed, `exd4` is illegal,
since White's own knight stands on d4; the reading that fits is `exf4` (the e3
pawn takes the knight that arrived on f4), and with it the whole line is legal
and ends in checkmate, confirmed by python-chess. Stockfish gives Nf4 as Black's
only winning first move: of 30 legal moves it scores +3.67, and the next best
(Qc1) scores −4.00. Along the line White's replies are heavily restricted: 2
legal replies after Rh6+ (Kg3, Bh4), 3 after Qd3+ (Qf3, Nf3, Be3, the last losing
to mate in two). The mate is not forced at the end: after ...g4 White has 16
legal moves, of which Rh2, Rg1, Be3 and Bg1 avoid mate (the engine prefers Rh2,
with Black still clearly winning, about +3.7), and Qxd3 allows Rh3#.

**The net report, computed by `mating_net_report.py` with no reference to piece types.**

| Term | Fool's mate | OTB Qxh7# | Puzzle 1, Qg7# | Puzzle 2, Rh3# |
|---|---|---|---|---|
| Check | open diagonal h4–e1 | open diagonal e4–h7 | adjacent queen | adjacent rook along the 3rd rank |
| Flights denied | own pieces on d1, d2, e2, f1; ray covers f2 | own rook f8, own pawn f7, queen covers h8, g7 | queen covers g8, h7 | own bishop f2, own rook g2, own pawn f4; rook covers h2, h4; x-ray through the king covers f3; Black pawn g4 (protected by f5) |
| Checker not capturable | no white piece reaches h4; the h1 rook is masked by the h2 pawn | king cannot take: Bd3 supports h7 | king cannot take: Rg1 supports g7 | king cannot take: pawn g4 supports h3; the white queen on d3 is masked by its own king on g3 |
| Block | none possible | not applicable | not applicable | not applicable |

Roles found by removing each piece in turn (which defence reappears): in puzzle
2, Bf2, Rg2 and the pawn on f4 uphold the flight term, the Black pawn on g4
upholds the capture term, and the Black pawn on f5 is needed because without it
Black's own king would be in check from the queen on d3. The white queen, the
d4 knight and the a2 and b3 pawns are irrelevant to the final mate.

**How the net is assembled, move by move.**

| Black move | What it does to the net | White's options |
|---|---|---|
| ...Nf4 | attacks Qd5 and Rg2; after exf4 a white pawn sits on f4, filling a future flight square, and e3 is empty, opening the 3rd rank | 38 legal, 12 avoid a forced mate, all clearly losing |
| ...Rh6+ | check along the h-file, opened by the knight's departure; Kg3 puts the king among Bf2, Rg2 and the f4 pawn | 2 legal |
| ...Qd3+ | check along the 3rd rank; the queen cannot be taken, because White's knight on d4 blocks the d-file; Qf3 interposes and is then pinned | 3 legal |
| ...g4 | attacks the pinned queen (its only legal moves are Qe3 and Qxd3), protects h3, fills g4 behind f5 | 16 legal, 4 avoid mate |
| ...Rh3# after Qxd3 | the queen arrives on d3, behind its own king, so it cannot cover h3 | none |

**Assessment — separate from the record above.**

- **The predicate carries over.** Four mates with different mating pieces
  (queen, queen, queen, rook), different kings and different surroundings
  decompose into the same four terms, and the roles are generated without any
  piece names. The lock here is the king's own pieces, Bf2 and Rg2, plus a
  pawn lured onto f4; in the fool's mate it was the king's four start-square
  neighbors.
- **New mechanisms appear, and each fits the §8t predicate as a restriction of
  the defender's coverage or an addition to occupancy or attacker coverage.**
  A decoy that fills a flight square (Nf4 and exf4 put the defender's own pawn
  on f4). Pinning the interposing piece (Qf3 pinned by the queen on d3). A
  defender masked by its own king (the queen on d3 behind the king on g3). A
  checker protected by a pawn (g4 protecting h3). The attacker's own shield
  can be part of the net (f5).
- **The assembly is forcing, and that is where the "race" shows.** Every Black
  move is a check or an attack on a queen, which is why White has 2 and 3 legal
  replies at the key points. The enabling moves that were voluntary in the
  fool's mate are here compelled by threats, which is the bridge to adversarial
  play noted in §8s.
- **Stated precisely, because it was described as the queen being sacrificed to
  position into the fortress:** the queen is given up only if White takes it,
  and White can decline with Rh2 or three other moves; Black still wins by a
  wide margin. The fortress in question is White's own king enclosure, built
  mainly from White's own pieces, and the sacrifices serve to fill it and to
  mask the defenders.
- **Independence from piece types holds for the predicate, not for the
  mechanisms.** A pin needs a slider on a line, a decoy needs a capture, a
  shield needs a line to shield. The geometry lives in the mechanism catalog,
  which is where it should be.
- **This is a weak test of the claim that the primitive can be found in any
  position where it emerges, because every position checked is a mate already
  on the board.** A stronger test is to recognize the net before the mate, for
  example from the start of this puzzle, by finding that the terms can be
  completed with forcing moves. That has not been done.
- **The g4 move does three things at once** (attacks the pinned queen, protects
  the mating square, fills a flight square). This is a second instance of the
  structure flagged in `boundary_board_states_case_study_otb.md` §7, a single
  move discharging several duties without conflict, with one difference: there
  the alignment was a coincidence and here it is chosen. Whether this adds
  anything beyond an ordinary tempo move that builds the net is open, and the
  decision recorded there is still pending.
- **Black has exactly one winning move at the start of this puzzle** (30 legal
  moves, Nf4 wins, the next best loses by about four pawns), so the position is
  a decisive only-move position. It is not a draw and is not a boundary state;
  the winning-move count there is a puzzle property, unrelated to the definition of a
  boundary state (`boundary_board_states.md` §14 item 8).

### 8v. Cataloging the mating net as a portable primitive

*Correction (§8ab).* The names P-002, P-004, P-005, P-007 and SUB-1 below are kept as recorded. They are not primitives: the catalog has three, the mating and capture net (P-001), chained dependencies (P-003, which absorbs what this section calls P-007, SUB-1 and the frozen structure) and the coverage race (P-008). The classification of each name is in §8ab and in catalog §1. Catalog section numbers cited here (§7, §7a, §7b, §7c, §8, §9) are those of the catalog before that reduction; the map to the present numbers is catalog §10.5.*

**The idea, as stated.** The discovery should be recorded in a scalable,
generalized way: as a modular primitive that can be extracted and deployed in a
catalog, depending on nothing but the base primitive and the way it expands to the
fool's mate and to any other position it applies to, with no further scripts to
validate it. A primitive applied to a valid boundary state will either not apply or
show, in explicit terms, an equilibrium in the composition of primitives. A valid
boundary state requires a composition satisfying narrow constraints for both sides,
with moves that sustain a loop of repetition, each returning to the same
equilibrium structure. Primitives can be discovered without boundary states, and no
position has yet been evaluated with one.

**Assessment — separate from the record above.**

- **The catalog is `boundary_board_state_primitive_catalog.md`.** It has four levels
  (primitive, mechanism, derived form, instance), an entry template that implements
  and extends §8i, an index of five entries of which one is written, and P-001 in
  full: the four-term definition, the mechanisms for each term with the instances
  that exhibit them, a defence mirror, derived forms, a procedure that evaluates any
  position without naming pieces, a composition interface, its role in
  boundary-state evaluation, and an evidence table. The fool's mate is the derived
  form DF-A and the instance I-1.
- **How it expands, which was the requirement.** The four terms are fixed. A new
  position adds an instance row; a position needing a new way to satisfy a term adds
  a mechanism; a new combination of known mechanisms adds a derived form; a new route
  to a decisive result adds a primitive. Nothing already written is rewritten.
- **Why the base primitive applies to any position.** The four terms are the rules'
  definition of checkmate, decomposed, so applying to any position is true by
  construction. The non-trivial content is the mechanism catalog and the cost model,
  and the cost model is not defined.
- **"Either does not apply or shows an equilibrium" holds for the mating route by
  the way the status vocabulary is built:** DORMANT is "does not apply", and ENGAGED
  with held terms is the explicit equilibrium. It has not been tested as a
  classifier, because no boundary-state candidate has been evaluated with it, and it
  says nothing about the other routes to a decisive result. That is why P-002 to
  P-005 are indexed.
- **The cyclic requirement is recorded as a definition in
  `boundary_board_states.md` §14 item 7,** and it makes P-005 (repetition loop) the
  next entry the catalog needs, since the loop moves must preserve the held terms.
- **No script was added.** `mating_net_report.py` is a reference implementation of
  the terminal case only; the catalog entry does not depend on it. The procedure for
  non-terminal positions (catalog §4.7) is a specification and has not been run.

### 8w. Applying the net to pieces other than the king

*Correction (§8ab).* The "race primitive P-002" named below is the tempo comparison, which is a term (TP of P-003, TM of P-008) and not a primitive. See §8ab.

**The idea, as stated.** The primitive should also apply to a piece other than the
king, for example an early-developed queen: a net that captures the queen, joint
compatible with the king net while making progress on other race conditions. In the
early game it lets the attacker unlock pieces from the initial position (development)
while putting pressure on the queen, forcing queen moves in place of the defender's
development. The same mating-net understanding can be joint compatible with the
queen-capturing net while progress is made on other primitives. The question is whether
this is the same primitive or a different one that can be jointly compatible.

**Assessment — separate from the record above.**

- **The game record supports the description.** In the first 12 moves White made 5
  queen moves (6.Qf3, 7.Qxd5, 9.Qe4, 10.Qf3, 12.Qg2). Three of them (9, 10, 12) were made
  with the queen under attack, answering ...Nb4, ...f5 and ...e4, which were played by an
  already-developed knight and two pawns. The first (Qf3) was voluntary and the second a
  recapture. In the same 12 moves White made two first-time minor-piece developments
  (Nc3, Bd3). This is a count from the game record, not a verdict on causation.
- **The structure is the same.** A net on a target is a conjunction of escape routes:
  the attack on the target, its flight squares, the capture of the attacker, and
  interposition. Mechanisms attach to terms, so the mechanisms of the catalog carry over
  unchanged: own occupancy of flight squares, attacker coverage, support of the attacker,
  masking, pin, decoy. The second puzzle contains a queen net in this form: after
  ...Qd3+ Qf3 the queen on f3 is pinned to the king, then ...g4 attacks it, and its only
  legal moves are along the pin line (Qe3 or Qxd3).
- **Two terms did not carry over, which is why a sibling was first proposed (superseded by the decision below).**
  (1) *Protection and exchange (P).* Capturing the king ends the game, so the king net
  has no exchange accounting. A queen net completes only if the capture wins material
  (the static exchange evaluation familiar from engines); protecting a queen against a
  pawn attack does not save it. P depends on a value scale, which the king case does not
  need. (2) *Counter-threat (K).* A king in check must answer. A queen under attack may be
  left alone if the defender has a check or a larger threat, so a counter-threat is an
  extra escape (revised in the Limits bullet below: counter-threats are a composition rule and not a
  term). The outcomes differ too: the king net ends in mate; the piece net ends in
  material gain, or, short of that, in tempo yield, where each residual escape the target
  uses costs the defender a move.
- **Decision and classification (the designer's).** The mating net and the capture net on any
  other piece are one primitive race condition. The catalog records them as P-001, version 2, the
  capture net, with the kind of target as a parameter (royal or non-royal). The former sibling
  specification P-006 is retired and its ID is not reused. The two properties that distinguish
  the kinds are the legal-move rule for a check and the outcome clause; neither is a separate
  primitive. The earlier proposal in this section, a sibling with its own terms, is superseded.
- **Joint compatibility has a precise form for the attacker.** For an attacking move m,
  let Res_K(m) and Res_Q(m) be the sets of defender replies that resolve the threat on
  each net. The nets are jointly compatible for the attacker when both are engaged and no
  single reply lies in both sets, so that every reply leaves at least one net
  unresolved. This is the double threat. The ...Nb4 fork in the OTB game (queen, bishop
  and pawn at once) is an instance with three nets, and ...g4 in the second puzzle is a
  move that advances a queen net (attacks the pinned queen) and a king net (protects h3,
  fills g4) together. This restates the fork in the catalog's vocabulary and is not new
  chess; the gain is that forks and mating nets share one schema and one compatibility
  condition.
- **Interference works in both directions.** A defender piece that serves both nets is
  overloaded. A residual escape of one net can also be steered into a mechanism of the
  other: in the second puzzle the queen's escape Qxd3 puts it behind its own king, which
  is the masking mechanism (X3) of the king net. Compatibility is therefore also about
  whether the residual escapes of one net fall into the mechanisms of the other. The
  second point is this assessment's own observation from one instance.
- **A tempo ledger gives the "race" a measurable content.** For each move, record
  whether it unlocks an own piece from its start square, which terms of which nets it
  discharges, and whether it forces a reply that cannot also be a development move. A
  move that is positive in several columns is the multi-duty move already noted in §8u.
  The count above is a hand-sized instance of the ledger. It has not been built or
  tested further.
- **Limits, as first written, and the correction.** As first written this bullet said that
  the exchange term needs a value scale and the counter-threat term needs a notion of a
  threat large enough to matter, so that the piece trap brings outside evaluation into a
  primitive the mating net kept free of it. That was wrong in kind: it treated the absence
  of entries from the catalog as an inherent need for evaluation. Corrected reading:
  *Counter-threats* are the same primitives applied with the colors swapped. A check is a
  king net the defender has engaged, and a larger threat is another net. Whether a
  counter-threat beats an attacker's net is a comparison of the tempo of two nets in
  opposite directions, which is the race. It is therefore a composition rule (catalog
  §4.10) and not a term of the piece trap. In the king net no such term was needed because
  the legal-move rule already forces a check to be answered. *The exchange term* is a
  chain of nested nets, since a recapture is a net on the capturing piece. What remains is
  an outcome measure, which in this framework would come from the effect of removing the
  piece on the composed race map. Static exchange evaluation is a stand-in until the
  catalog can supply that.
- **What has to exist for these to be derivable, stated as dependencies and not as a
  limitation of the primitive:** the cost model for tempo (unfinished, catalog §4.8);
  the race primitive P-002, to compare nets by tempo; the other primitives, so that the
  composed race map exists; and a well-founded ordering for the mutual recursion between a
  net and its counter-nets, for which tempo is the natural candidate (§8g). That a
  complete catalog derives the exchange outcome and the counter-net comparison without
  outside evaluation is what the definition of a boundary state requires of it (§8x). It is
  recorded entry by entry through the residual interface (§8x); the catalog is being built to
  that size.
  The two instances of the piece trap so far are the queen in the second puzzle and the
  queen hunt in the game record. Nothing here bears on whether any position is a boundary
  state, and it uses neither the width criterion nor any claim about the start position.

### 8x. Recording the composition of a proven boundary state, and the architecture it sits in

**The idea, as stated.** Once all primitives are cataloged and boundary states are proven in
the catalog, the composition of primitives confirmed for each proven state can be denoted.
The complexity most likely comes from how the primitives are composed into jointly compatible
parts. This is not a hypothesis. It is what the definition of a boundary state requires and has
been understood throughout the project. The engineering problem is how primitives compose to
emerge a boundary state, and how they completely enumerate the relations of the full
equilibrium position. The architecture, stated by its designer: there is a primitives catalog;
there is a candidate boundary state catalog, a proven one and a disproven one; the composition
of the primitives is understood in the proven boundary state catalog; and joint composition
implies that lacking primitives are reducible. If a, b and c are composed into a boundary state,
then given the composition of a and b the missing primitive c is emergent from joint
compatibility. This had been iterated in the past.

**Assessment — separate from the record above.**

- **The architecture is the designer's and is unchanged.** It matches §3 of
  `boundary_board_states.md`: a candidate pool built from necessary-not-sufficient filters,
  exhaustive verification of each candidate, the verdict cached permanently whether boundary or
  not, and a pool that shrinks as verification proceeds. In the terms used here: a primitives
  catalog (`boundary_board_state_primitive_catalog.md`); a boundary state catalog with
  candidate, proven and disproven entries; and the composition of primitives recorded on proven
  entries.
- **A correction to this record's earlier handling.** An earlier version of this section, and a
  section of the catalog file, put a record specification for proven entries inside the
  primitives catalog; gave it status labels of its own (composed, partial, uncomposed) in place
  of candidate, proven and disproven; and described missing primitives as things to be written
  before entries could be complete. All three departed from the architecture. The specification
  was removed from the primitives catalog, the status labels were dropped, and the third point
  is replaced by the reduction rule below.
- **The reduction rule, as designed.** If a proven entry is composed of primitives a, b and c,
  and the catalog holds a and b, the missing c is emergent from joint compatibility: its
  interface is fixed by the residual interface that a and b leave in the proven entry, in the
  same way that a with c fixes b in the complete-graph relation of §8g. The proven catalog is
  therefore a source of primitives as well as a consumer of them, and the subtraction of §8l run
  on a known composition generates candidate entries. The form of the constraint system that
  fixes c, and whether it fixes one primitive or a family that the proven entry and its
  neighbors then select among, belongs to the composition-rule work.
- **The composition record carried by a proven entry** (a proposed form; the entry's existence
  and its proven status are unchanged). At each node of the equilibrium loop, every legal move of
  the side to move is an equilibrium move or a losing deviation. Each losing deviation is mapped
  to the primitive whose net or race completes against it, with the term, mechanism and depth.
  This places no requirement on the number of equilibrium moves. The deviations with no
  refutation among the known primitives form the **residual interface**: the specification of
  the primitives not yet in the catalog, which the reduction rule solves for. The record also
  holds how the primitives relate: shared resources, double-net relations, escape steering, and
  the ordering of the recursion.
- **The proof and its insurance.** A completed composition is itself the case-list proof of §8d:
  every deviation of each side is mapped to a primitive race condition that refutes it, and the
  equilibrium moves keep the invariants of every primitive. Exhaustive search is insurance. It is
  recommended on the first proven entry, to confirm the composition rules are right, and used
  afterward as needed. A rule that an exhaustive search contradicts is revised.
- **Where the complexity sits is agreed:** in the composition rules, meaning resource accounting
  across primitives and its higher-order constraints (§8g), the double-net condition (§8w),
  escape steering, and the ordering of the recursion, and not in the individual primitives.
- **Earlier wording is replaced.** §8w and the catalog had called the derivability of the
  exchange outcome and the counter-net comparison the programme's hypothesis. It is what a
  complete catalog is required to achieve, checked entry by entry through the residual interface.
- **No proven entry exists yet.** The project has not proven a global boundary state, and the
  earlier king-and-pawn data are local width measurements. Under the reduction rule the first
  entries do not have to wait for the other primitives to be written: their residual interfaces
  are where those primitives are extracted from.

### 8y. Activation, the primitive state, and perturbation

*Correction (§8ab).* The names P-002, P-004, P-005, P-007 and SUB-1 below are kept as recorded. They are not primitives: the catalog has three, the mating and capture net (P-001), chained dependencies (P-003, which absorbs what this section calls P-007, SUB-1 and the frozen structure) and the coverage race (P-008). The classification of each name is in §8ab and in catalog §1. Catalog section numbers cited here (§7, §7a, §7b, §7c, §8, §9) are those of the catalog before that reduction; the map to the present numbers is catalog §10.5.*

**The idea, as stated.** Another primitive, identified from endgame logic, is pawn promotion. In
king and pawn against king and pawn it is not a mating or capture net but a segmented process in
which the net depends on a race condition for promotion. It is also directly related to piece
activation and to a locking mechanism: there is a race condition that unlocks a piece or puts one
into play, whether by development from the start position or by promotion, which requires
restructuring pieces according to the opposition's structure and the composition of primitives on the
board. An embedded race condition, dependent on adversarial play, defines the viability of these moves
and depends on the primitive structure of the board as a whole. The dependency is about the symmetry
and asymmetry of the viable mitigation and obstruction of primitives, which is more informative than
pieces. The primitives that matter here include pawn promotion, development from the start position,
and the capture and mating net relation. A complete relation of these race-condition dependencies,
including the adversarial ones, together with how a piece move perturbs them, may be where chess
strategy emerges, as a search on primitive calculations over a complete graph. Gaps in it inform which
primitives are not yet considered. The capture and mating net is the foundational primitive, the reason
anything is done. What never changes, the geometry of files a to h and ranks 1 to 8 and the rules of
piece movement, is where the space for perturbation has to be filled. Moving towards an equilibrium in
the composition of the primitives is moving closer to a boundary state.

**Assessment — separate from the record above.**

- **Followed, and integrated.** It fits the architecture without changing it: a new entry in the
  primitives catalog (P-007, with a dependency map in catalog §8), a new method for candidate
  generation (Method 7 in `boundary_board_states.md` §2), and no change to the boundary-state catalogs.
- **Classification, offered as the reading the entries use and open to correction.** Promotion and
  development are the same kind of primitive, *activation*: a piece's legal coverage is enlarged from a
  locked state to an active one. Its route is a race, which is P-002 (the tempo comparison, unchanged); its
  obstruction is a capture net on the activator (P-001); its lock states come from frozen structure and the
  piece's own blockers (P-003). A promotion is an activation whose route is a P-002 race. P-002 keeps the
  race calculus, and P-007 is what the race is a race for.
- **The segmentation is concrete in the endgame.** In king and pawn against king the downstream net, a
  queen against a bare king, is completable whenever the new queen is not immediately lost and the position
  is not a stalemate (standard endgame knowledge), so the activation race alone decides the result. In
  richer positions the strength of the net depends on which pieces have been activated, and the
  opponent's net on the activator is the obstruction. The dependency is mutual, and tempo makes it
  well-founded.
- **Activation changes both sides' nets.** A developing move supplies mechanisms to one's own nets and
  removes one's own locks, which are the own-occupancy and masking mechanisms the opponent's net would use
  (F1 and X3 in P-001). In the fool's mate the lock is the defender's own undeveloped pieces, and Nf3
  discharges terms of the net by geometry (not run).
- **Symmetry and asymmetry as a state and not a piece count.** The primitive state S(P) is, for each
  side, the status of every catalog primitive that engages (dormant, engaged with held terms and tempo, or
  completed). The asymmetry is the componentwise difference between the two sides' states. It is
  piece-independent: the same asymmetry can arise from different piece sets. The game record is one
  instance. In the first 12 moves White has 3 activation events (Nc3, Qf3, Bd3) and 5 queen moves, 3 of them
  made under attack, while Black has 5 activation events (Nc6, Bc5, Nf6, O-O, Qg5).
- **Perturbation is the change in that state under a move,** ΔS(P, m) = S(P·m) − S(P). The invariants
  (the board geometry and the movement rules) are what every catalog mechanism derives from, so the
  effects of perturbation can be enumerated over a finite geometry. The mechanism tables of the catalog are
  that enumeration, and filling them in is the work the statement describes.
- **Moving towards equilibrium is the order of the reduction, not an approximation.** A measure of
  distance from a position's state to an equilibrium composition is not yet defined, and it needs the
  cost model (catalog §4.8). Whatever its form, it only chooses which residual to reduce next and never
  changes what is true. A boundary state is proven by its composition, as a complete case list, with
  exhaustive insurance on the first entry (§8x).
- **Gap detection is a deduction under the definition.** If the catalog's primitives are equalized at a
  position whose verdict is proven, by tablebase or exhaustive search, not to be an equilibrium, then some
  primitive race condition decides it, and it is one of three: a primitive not yet in the catalog, whose
  joint compatibility with the known ones carries the imbalance; a known primitive that is not jointly
  compatible with another in this position and is decisive; or a catalog entry that was evaluated wrongly.
  The first two are the designer's statement of it; the third is a correction to an entry. All three
  develop the protective belt, the primitives, and none touches the hard core, the definition of a
  boundary state. The residual interface fixes what the missing or incompatible primitive must discharge.
  An engine verdict can point to where to look; the gap is established by the proof.
- **What this does not change.** The definition of a boundary state, the architecture of the catalogs and
  the proof-and-insurance rule are unchanged. None of P-002 to P-005 is written, and the cost model is not
  defined, so Method 7 has not been run.

### 8z. Segmented dependencies generalized: chains, the movement relation, and the substrate

*Correction (§8ab).* The names P-002, P-004, P-005, P-007 and SUB-1 below are kept as recorded. They are not primitives: the catalog has three, the mating and capture net (P-001), chained dependencies (P-003, which absorbs what this section calls P-007, SUB-1 and the frozen structure) and the coverage race (P-008). The classification of each name is in §8ab and in catalog §1. Catalog section numbers cited here (§7, §7a, §7b, §7c, §8, §9) are those of the catalog before that reduction; the map to the present numbers is catalog §10.5.*

**The idea, as stated.** Chain structures are a key structure of the locking mechanism. A passed pawn rests on a chain
structure, the pawn wall, with the pawn passing unopposed. A knight outpost is a chain structure: taking the knight
allows the capture by a pawn, and the weak point is the chained end of the structure. This informs the segmented
dependency primitive directly, and the intuition is that it should be integrated into it. A segmented dependency chain
informs race conditions directly through the capture mechanisms inside the chain, and the capture potential of a chain
completely alters the viability of the race condition that depends on it.

A second family is the movement rules of the pieces. A dark-squared bishop navigates only dark squares. The position of a
knight changes the routing, and without the routing move the routing is predictable from the move rules, which affects the
joint compatibility of the primitives with the composition of chained structures. The rook navigates files and ranks, and
the opposition of occupancy and coverage shapes the viability of its race through chain structure and coverage nets. The
queen navigates the chained structures with the fewest restrictions other than the knight, and its routing moves shorten a
race relative to rook-only or bishop-only movement. All of this generalizes into a primitive race-condition structure at
the granular piece level while staying piece-agnostic: the hierarchy of the primitive down to the granular understanding
of segmented-dependency races, by piece perturbation and the rule sets of the pieces relative to specific structures.

The constraint on every addition is stated by its designer: it must be modular, piece- and board-agnostic, and define a
primitive race condition. It need not work identically for every piece. There are specialized rules for certain pieces (the
passed pawn, the movement difference of the knight), but they are specializations, and the goal is the most modular and
general form. A note on one piece's color or parity alone does not meet that constraint.

**Assessment — separate from the record above.**

- **Followed, and integrated, in the belt only.** The additions are: the substrate SUB-1 (catalog §7a), which derives from
  the movement relation of any unit its orbit, route-length parity, the dependency set of each edge, the permanence of its
  coverage and route distance under inclusion; the entry P-003, specified with the chain as its core (catalog §7b); new mechanisms for
  P-007 (RT3 to RT5, LK3, OB4, OB5) and a derived pair of forms (DF-P1, DF-P2); the mechanism O3 for P-001; and rows in the
  dependency map (catalog §8). The definition of a boundary state, the architecture of the catalogs and the proof-and-insurance
  rule are unchanged. Nothing in the earlier documents was deleted. The one existing cell that was replaced is the state of
  P-003 in the catalog index ("not written"), because its state changed; the other existing lines were extended by appended
  sentences.
- **The two families are one construction at two layers, and they meet at the dependency set.** The chain is a dependency
  graph among units, with two edge types: ENABLE (an element's removal or occupation is required before a stage opens, which is
  the lock removal of P-007) and SUPPORT (an element's presence makes the capture of another unfavorable, which is P-001's X2
  and outcome clause). The movement rule of a unit is a navigation graph of its own. A route edge has a dependency set, the
  squares whose occupancy decides whether it exists, and a wall or chain that stands on that set is an ENABLE edge to the
  route. That is the mechanism by which a chain "informs" a race: it cuts the edges of the relations that travel along a line
  and leaves the relations that leap unaffected except at the destination.
- **A correction to this record's earlier handling, stated plainly.** The first version of this assessment, given in the
  discussion, treated the color symmetry of one kind of unit and the parity of another as separate piece notes, and then made
  piece-specific route claims. The designer's correction is that the substrate must be one construction on the movement
  relation and that the primitive catalog must state no rule of any piece. Piece-specific special cases are inputs supplied
  with each unit's relation when a primitive is applied to a position, and they are accounted for during composition and
  race calculation on that board structure. Recording them in a primitive turns it into a subspace of piece composition and
  defeats its purpose. The catalog now states the construction and no piece rule. On route distance, the generic fact is
  inclusion: a relation that includes another's edges gives a distance no greater, by definition, and the stalemate and
  zugzwang considerations enter through composition. A register of piece-specific special cases kept apart from the
  primitives, as an input they consult, is a possible modular addition (catalog §9.3) and is not started.
- **Classification is open.** The designer's reading is that the content belongs inside the segmented-dependency primitive.
  The catalog classifies SUB-1 as a substrate and puts the chain in P-003, because by its growth rule (catalog §1) these supply
  mechanisms and costs to existing terms and decide no result of their own. P-003 itself may turn out to be a substrate on the
  same grounds (catalog §7b.11). The definitions serve either placement, and the choice is recorded as the designer's.
- **The chain closes an open item.** The outcome clause O of P-001 uses static exchange evaluation as a stand-in "until the
  catalog supplies" a composed race map (§8w, catalog §4.3). A SUPPORT chain with its unsupported end and its critical nodes is a
  candidate for supplying it: a capture at a node is answered by the supporter, and the cascade ends at an unsupported end or
  at a node whose loss completes another net. SUPPORT edges are read from coverage and the cascade removes a node at each step,
  so the reading is well-founded and does not circle back to O. This is a hypothesis and has not been run.
- **Permanence is this record's inference, offered for correction.** A passed pawn is a route whose squares no opposing pawn
  can ever cover, and a pawn cannot move backward, so its coverage can only be lost. A unit whose orbit excludes a square can never cover it. Both facts are read from the relation: the first from irreversibility (an acyclic navigation
  graph), the second from orbit confinement. Coverage by a unit with a reversible relation and one orbit is transient, since
  the unit can leave and return. That gives a piece-agnostic form of "passed" (a route with no permanent denial) and of "outpost"
  (a square permanently uncovered by irreversible units, with a SUPPORT edge on the occupant, whose weak end is the unsupported
  end of the supporting chain), and it matches the base-of-chain attack already in §8b's list. It must be extended to the
  region where a relation changes, since a promotion adds coverage.
- **A mirror of overload, in the graph.** Overload (§8f) is a node that serves two edges it cannot serve at once. The single
  move that discharges several duties at no cost (the OTB case study §7, and §8u) reads as a move that adds or satisfies two
  edges at once. This is noted as a reading. Whether it is a named pattern is still the open decision recorded in the case study.
- **Two predictions, stated in advance (catalog §7a.10).** (1) In the ablation of §8s, the minimal sets whose removal changes
  the verdict coincide with the cut sets of the graph read from the same position by P-001's mechanisms. (2) Orbit exclusion
  (RT3) reproduces the known color-lock verdicts across the whole class, including the color-lock fact found inside Uri's rule
  (§8n), which this record reads as an instance of orbit exclusion and which has to be checked against §8n. Neither is run.
- **Where this stands in the programme's terms.** It is a development of the protective belt that leaves the hard core
  untouched. By the stricter standard of `boundary_board_states.md` §15 it is not yet progressive: nothing in it has been run, and definitions alone are not novel facts. The two predictions above are the novel-fact tests, and they can fail.
- **What this does not change.** The definition of a boundary state, the architecture of the catalogs, the proof-and-insurance
  rule, and the status of P-002, P-004 and P-005 as not written. SUB-1 and P-003 are SPEC.


### 8aa. The coverage race: a third primitive derived from two, adversarial composition, and scoped coverage

*Correction (§8ab).* The names P-002, P-004, P-005, P-007 and SUB-1 below are kept as recorded. They are not primitives: the catalog has three, the mating and capture net (P-001), chained dependencies (P-003, which absorbs what this section calls P-007, SUB-1 and the frozen structure) and the coverage race (P-008). The classification of each name is in §8ab and in catalog §1. Catalog section numbers cited here (§7, §7a, §7b, §7c, §8, §9) are those of the catalog before that reduction; the map to the present numbers is catalog §10.5.*

**The idea, as stated.** With two primitives, the capture and mating net and the chain structures with their segmented
dependencies, a third is derived under joint compatibility, and it demonstrates the complete graph of the primitives in a concrete
case. The composition, together with the formal chess understanding of removing the defender, gives a new race condition: the
coverage of a chained structure, and how the structure can be manipulated by removing its key defenders. The coverage race after
captures and trades, based on the segmented dependencies of the chain, translates the board into a new position. It affects the
emergent piece locking and the piece activation, and it matters under optimal play as a stress test of how effective a chained
structure is. If a chained structure is covered, but trading one of the pieces leaves the coverage of another piece no longer
sufficient, the race condition on the trade is won. So the structure requires additional coverage, and the question is whether the
chain holds if a key defender is removed. This is not a combinatorics problem. It is a graph-related coverage compositional
problem, viewed through race-condition calculation and not through exhaustive combinatorial coverage.

The architecture requires modularity. If two primitives overlap, then packaging several primitive dependencies into one excludes
emergent possibilities, so the overlap compromises the architecture and must be fixed.

The primitives are composed adversarially. The calculations are not performed on one side only: the same calculations are performed
from the other player's end, with the perturbation of that side's candidate moves, so that the structure is informed by the
relativistic structure played against it. The removal of a defender depends on the opponent's options as those primitive
calculations give them. It is also not exclusive to captures. An opponent can place a piece that threatens the structure, and after
the piece is captured the occupancy of the capturing piece is incompatible with the chained structure's race, so the coverage race
is broken without a defender being captured explicitly, by forcing a new piece into the structure's dynamic.

Where combinatorics reaches a limit because of the perturbation, a hybrid approach is viable, which is the scoping of exhaustive
coverage already stated. If the composition of primitives requires minimal combinatorial calculation that mitigates the search
effectively, it remains viable, and it is kept open that additional primitives may mitigate it further. With that understood,
adding it to the documents is a sufficient modification of the protective belt.

The approach is race-condition primitives and not a heuristic approximation. What matters is making moves under perturbation and
understanding the effects on the primitive race conditions under perturbation. Certainty is derived from the growing catalog's
account of how each primitive race condition is affected by a perturbation, so recaptures and static conditions are not a separate
uncertainty. The hybrid is exhaustive coverage formed by the primitive race-condition understanding: the calculation of the
primitives under perturbation is integrated into the process, and that is what mitigates the search.

**Assessment — separate from the record above.**

- **Followed, and integrated, in the belt only.** The additions are: the entry P-008 (catalog §7c), with its register of held terms and their support, its effect vector, its four terms, its mechanisms and its two-sided evaluation; the narrowing of P-003 to version 2 (catalog §7b); new
  rows in the dependency map and the index; a convention on adversarial evaluation (catalog §6); and the provenance of §9.4. The
  definition of a boundary state, the architecture of the catalogs and the proof-and-insurance rule are unchanged. Nothing was
  deleted: the material moved out of P-003 is mapped item by item to its new place (catalog §7b.5 and §7b.6).
- **The reduction rule, exercised on two written entries.** The composition of P-001 and P-003 leaves a residual interface, and
  the residual fixes what the third entry must discharge. P-001 reads coverage as sets at one position. P-003 reads support as
  present or absent. A net on a supporter is an exchange, and an exchange changes the position. The interface that remains is: does
  the coverage a structure needs still suffice afterwards, can the disturbing side realize that before it is repaired, and what do
  both sides' options say. That is the designer's derivation. It is a derivation by the definitions and is not an empirical fact
  by itself (§8x).
- **The overlap was a defect in this record's earlier specification, and it is fixed.** Version 1 of P-003, specified in §8z, held
  the dissolution race of a structure. By the modularity rule that overlapped the coverage race. The race moved to P-008, and P-003
  keeps the static structure: the nodes, the roots, the critical nodes, the cascade closure. Consequently P-003 decides no result of
  its own, and whether it stays an indexed primitive or is classified with SUB-1 as a substrate is the designer's decision, noted as
  such.
- **The primitive, as now specified.** A register of held terms with their residual replies (the held terms of P-001 at §4.4, the
  edges of a P-003 structure, the terms of P-007); the support of each reply, which is the elements its existence and efficacy read,
  taken from the edges and dependency sets of SUB-1; the effect vector of a perturbation over the registers of both sides, each
  entry UNCHANGED, HELD or OPENED by the owning primitive's definition; opening sets and key elements; joint support, which is
  "coverage is used once" and overload with no special case; and the four terms OP, CR, FS and TM. A perturbation whose change
  misses the support of a held term leaves it held, and a held term is opened only by a perturbation whose change meets the support
  of each of its residual replies. Both follow from the definition of support.
- **Adversarial composition, and its coupling.** The same primitive is computed for both sides with the colors swapped, and the
  disturbances available to one side depend on whether its own structure still holds after it spends the units the disturbance uses.
  That feasibility term couples the two evaluations. Each disturbance spends a unit or a tempo, so the coupling is well-founded by
  induction on units and tempo, as in §8y. That is argued here and not proven. It restates for all primitives the counter-threat
  composition already recorded in §8w.
- **Beyond capture.** The primitive is defined by the effect of a perturbation on the held terms: a removed or displaced element,
  a changed occupancy on a dependency set, a changed exposure, an opposing entry. The carrier that compels it (a check, an attack, a
  sacrifice, a pin, or an offered piece) is read from P-001's carriers. The case of a piece placed to threaten the structure, whose
  capture leaves the capturer incompatible with the structure's race, is then a forced-in arrival whose change combines the removal
  of an element from its post with a change of occupancy and exposure (catalog OP6).
- **A second correction to this record's handling, stated plainly.** The first version of this assessment and of the catalog
  entry described the supply and demand condition as a first-level condition that a scoped search then confirms or overturns, and
  proposed the rate of exceptions as a measure. That presented the primitive as an approximation to a search verdict. The designer's
  point is that primitive race conditions are calculated under perturbation: the effect of a move on each primitive race condition
  is what matters, the condition is evaluated exactly at every position, and certainty about a line comes from those evaluations
  and grows with the catalog. The entry and this record were restated, and the withdrawn wording is noted in catalog §7c.13 and
  §9.4.
- **A third correction to this record's handling, stated plainly.** Version 1 of the catalog entry defined the coverage condition
  through a supply-and-demand formalization, a Hall-form condition with a slack, a robustness and a tight set. That was imported from
  outside the race semantics of the primitives, and it is not the primitive architecture. The designer's point is that the
  architecture is built from the primitives' own race conditions under perturbation. The entry is restated in terms of held terms,
  support and effect vectors, and the withdrawn formalism is recorded in catalog §7c.13 and §9.4.
- **The hybrid, as the designer states it.** The primitive race conditions are calculated under perturbation, and exhaustive
  coverage is formed by that understanding. Concretely: the mechanism tables enumerate the effects a perturbation can have on a
  structure, so a perturbation with a null effect on every engaged primitive in the footprint is covered without branching (a null
  effect is a change that misses the support of the condition, catalog §7c.9); a held term is opened only by a perturbation whose
  change meets the support of each of its residual replies, which decides whole classes of lines DORMANT by derivation; and over the
  perturbations whose change meets the support closure, the coverage is complete, each case decided by primitive calculation. That
  is the case-list proof form of §8x, with exhaustive search as insurance on the first entries. What bounds it is the footprint of
  §8h, with its four conditions, and the completeness of the effect enumeration, which is a claim of the catalog. When the closure
  grows to the board because kings are involved or lines reach across it, the number of perturbations to branch grows with it, which
  §8h condition 3 already warns of. Further primitives may reduce it, and that is kept open as the designer has it.
- **What the catalog does not yet supply, as gaps and not as uncertainty in the calculation.** The outcome clause of P-001 still has
  only a static exchange count as its stand-in, because its composed race map needs the effect of nested nets under perturbation,
  which P-008 supplies, together with P-002 and the cost model, which are not written. The effect enumeration of the coverage race
  (the five perturbation classes of catalog §7a.7) is a claim of the catalog's completeness. Whether the carrier is a term or an
  input is open. Each of these is a residual interface for the reduction rule, and discharging them is how the catalog grows (§8y).
- **Three predictions, stated in advance** (catalog §7c.13): the key elements of the held terms coincide with the load-bearing
  elements found by ablation (§8s); the evaluation of the primitive under perturbation, composed over a line, reproduces the verdict
  of exhaustive search for exchange-dependent structures, a disagreement being a residual interface that develops the catalog; and
  every verdict-changing perturbation in the ablation and deviation data has a change that meets the support closure, a failure
  specifying by the reduction rule what the catalog lacks. None is run.
- **Where this stands in the programme's terms.** It is a development of the protective belt and leaves the hard core untouched.
  It is not yet progressive by the stricter standard of `boundary_board_states.md` §15, because nothing in it has been run and the
  derivation of a primitive from two others is a derivation under the definitions. The predictions are the novel-fact tests and
  they can fail.
- **What this does not change.** The definition of a boundary state, the architecture of the catalogs, the proof-and-insurance rule,
  and the status of P-002, P-004 and P-005 as not written. P-008 is SPEC.


### 8ab. The reduction to three primitives

**The idea, as stated.** There are three primitives: the mating and capture net, chained dependencies (the segmented dependencies,
already detailed), and the coverage primitive. The repetition loop is the definition of a boundary state, the emergent effect of
the primitives being in equilibrium, and it is not a primitive. Activation describes the result of unlocking a locked piece from a
locked structure, which is the result of a primitive and not a primitive. Edge restriction has no race condition. The hard core is
the boundary board state, the protective belt is the primitive catalog, and the requirements on an entry are that it be a race
condition and that it be modular, generalizable and piece- and board-agnostic. Entries that are descriptors of primitives, that are not
modular or that are surface-level are not primitives and spoil the programme.

**Followed, and integrated.** The catalog is rebuilt on those three. The test an entry has to pass is written into the catalog as the
admission rule (catalog §1): a conjunction of terms with mechanisms that decides a race under perturbation; irreducible by the
reduction rule; not a description of a result or of an emergent state; piece- and board-agnostic; and modular, with no overlap. The
classification of what the earlier index carried:

| Earlier entry | Classification | Where its content lives |
|---|---|---|
| P-001 Capture net | primitive | catalog §4 |
| P-003, with P-007 and SUB-1 | one primitive, chained dependencies: the activation entry, the substrate and the static structure were three pieces of it | catalog §7 (version 3) |
| P-008 Coverage race | primitive | catalog §8 |
| P-002 Race to a square | the tempo comparison: a term (TP, TM), since race is what every primitive is | catalog §7.5, §8.5 |
| P-004 Edge restriction | an input: missing edges of the movement relation | catalog §7.3 |
| P-005 Repetition loop | the definition of a boundary state; emergent, so the hard core | `boundary_board_states.md` §14 item 7; catalog §5 item 2 |
| P-006 | retired; merged into P-001 version 2 | catalog §4 |
| P-007 Activation | the result when a chain of enabling dependencies resolves | catalog §7 (terms RT, LK, OB, TP) |

**What went wrong, recorded.** The first catalog indexed P-002 to P-005 as primitives to be written, from the families of §8b and the
cyclic definition, and later sections added P-007, SUB-1 and the narrowed P-003 to that index without testing each against the
requirement that a primitive is a race condition with terms. The chained-dependencies primitive, which the designer had described as
one, was split into three entries by that growth. The index was then taken as given by later sessions. The correction restores the
architecture the designer had stated at the start.

**What is unchanged.** The hard core and the definition of a boundary state; the architecture of the catalogs; the reduction rule and
the proof-and-insurance rule; every term, mechanism, derived form, instance, definition and prediction of the merged entries, which
were moved and not discarded (catalog §10.5 maps each earlier section to its new place); P-008 (renumbered, no definition changed);
the historical sections of this record, which keep the names as recorded and carry correction notes. The earlier index is kept in
`backup_before_three_primitive_restructure/`.

**Left open, and recorded as a residual and not as a fourth primitive.** Zugzwang without capture is not read by any of the three.
By the reduction rule it comes back as a residual interface if a proven boundary-state entry needs it, and it is admitted then only if it
passes the admission rule.

### 8ac. The anchored opening survey: local and global by path, and the position-keyed cache

**The idea, as stated.** The methods so far looked from the far end, retrograde from mating positions and the mitigation of
complexity. They ignored the anchoring point of every game, the initial position, and the navigation from it down a chosen
opening line. The anchored survey takes a well-understood opening, chosen for a reason, and uses the primitive catalog and the
perturbation process along it to find a boundary-state candidate. It is a survey for candidates and is not yet a tool.

How a state is reached matters for classifying it. A global boundary state is a member of a reachable family that is reachable from
the initial position under the optimal backtrace constraint (survey §5a), so it depends on how it is arrived at; it requires proof that it is a boundary state
and that its family reaches the initial position. A global state is sound holistically: reached by non-losing
play and holding against all deviation replies. A deviation may lead to a different global boundary state, which still closes the
phase between choice and force, so there are many global states and avoidability is not a ground for calling a state non-global. A
local boundary state is one where the equilibrium is verified (the primitives show it, and exhaustive verification shows they did their job) without requiring how it was reached or
that its family reaches the initial position. The survey's target is a local candidate anchored from the initial position by a chosen line, optimal or
sub-optimal.

Forcing sequences are read as race conditions. A forcing sequence is a decisive advantage to the side that forces it, and it may
force a win, a draw or a loss. A forced loss is uninformative. A forced win marks a winning race condition, and adversarial play
takes it, which prunes the lines that allow it. A forced draw is a possible boundary state, and a draw forced at each step to
prevent a forced loss is boundary-state behavior. The theory tree and perturbation up it locate where choice becomes force.

A cache keyed on the position records work done. If force is already deployable from a position, or a boundary state has been
passed, so that every branch below is a forced sequential win and the only decisive moves are forced out of zugzwang-type
situations, then any path that reaches the position inherits that, and the work below it need not be repeated. A found local
candidate is cached, and it may or may not be global. If it is concluded exhaustively that no local candidate can exist below a
position, that work is not repeated either. Where the candidates at a depth form a wall of boundary-state candidates in a cycle,
with interdependencies that branch into a large combinatorial space, the choice of path is bad, because it is too far upstream to
verify. Choosing a line and a depth is a calibration problem. The survey of lines is to be taken from a book of opening theory.

**Followed, and integrated.** The method is written as its own document, `boundary_board_states_opening_survey.md`, with the
definitions above, the procedure, the three cache entries (forced region, boundary found, none can exist), the calibration, and the
use of the theory tree as material (deviation map, width, evidence tier). Pointers were added to `boundary_board_states.md` (the
thirty-first addition, the resume box, and §14 item 9). The primitives, the hard core and the other definitions are unchanged.

**Assessment — separate from the record above.**

- **Part of this is already on record.** Forward anchoring and the discovery/certification split are in §8r, and start-reachability
  is in §5 (Tier 3 of `boundary_board_states.md` §14). The survey elaborates §8r's "with a chosen opening line, each piece move's
  effect on the primitives can be tracked" into a procedure. What is new is the path-dependent classification, the cache and the
  calibration.
- **This is consistent with §3.** A state's own status is local and unconditional there. The classification here is of how a
  verified state is reached, and it does not touch the status.
- **A reading offered for correction.** *Global* is read as start-reachable in the sense of §5, as membership of a family reachable from the initial position under the optimal backtrace constraint (survey §5a). `boundary_board_states.md` §14 item 4 asked for a firm definition of local and global, and this is offered as one. The
  record's "local" in §3 (a status that is local and unconditional) and this survey's "local boundary state" (verified, with no
  claim about arrival) agree.
- **Avoidability is not the test (the designer's correction).** A reading proposed in discussion, that a candidate which a side can route around is for that reason non-global, was withdrawn. A detour leads to another global boundary state, which still closes the phase between choice and force; global states are many, and a detour does not make a candidate non-global.
- **A scope point in the cache.** The forced-region entry is a claim under optimal play. A drawn position that the winning side
  could reach by an error is not on an optimal continuation and remains admissible as a local candidate. The entry says nothing
  about it. Whether the survey wants such candidates is left open in the survey document.
- **A guard rail from §8p applies.** A cache claim about everything below a position needs a certificate. A search that finds no
  boundary state below a node has found absence of information, and an assertion in the theory tree is an assertion.
- **A correction recorded.** In discussion, the global and local classification was first read as a scope tag on a verdict, on the
  ground that a position's status is independent of how it is reached. The designer corrected it: the status is the position's, and
  the classification is path-dependent by definition. The first reading was withdrawn.
- **Not run.** No line has been surveyed. The calibration signals of the survey document are criteria and not results.

### 8ad. Forced terminal draws as boundary states

**The idea, as stated.** By the definition, a forced stalemate or a forced dead position, that is, a forced draw, is the same as a
boundary state when it arises where there is no choice to deviate. The assumption that they are not boundary states comes from
classifying them as sub-optimal relative to other potential branches. Under exhaustive confirmation of all branches, some branches may
force stalemate or a dead position, and that is not evidence against the classification. The definition as stated should not exclude
them. What it informs is that the boundary-state catalog is not to be spammed with stalemate positions treated as discoveries.

**Followed, and integrated.** Recorded in `boundary_board_states.md` §14 item 10, with correction notes at the places that exclude
these draws (§1, item 7, item 8), a sentence in catalog §5 item 2, and the survey document (the forced-draw node, §3, and §10). The
equilibrium condition is unchanged: the position is drawn, each side has a losing deviation, both sides are live. What changes is
the realization: the cycle of repetition is the paradigm and the forced terminal draw is a second.

**Assessment — separate from the record above.**

- **This closes the realization question** recorded when the closure claim was discussed: a route into a forced liquidation or a
  forced stalemate with no avoiding move is covered by the closure, and the definition does not need a cycle on every route.
- **The filter is by choice and not by terminal type.** A terminal draw that could have been avoided is not an equilibrium, so it is
  not a boundary state. A terminal draw that cannot be avoided, with a losing deviation on each side, is. The 1,139 king-and-pawn
  positions in which the saving move captures the pawn into a dead draw remain excluded, by the stronger side being unable to lose.
- **A condition to check on each such candidate:** that the equilibrium holds upstream of the terminal position, that is, that both
  sides have a losing deviation and no non-losing move avoids the termination. The terminal position itself is the end of the draw
  and not the equilibrium.
- **Not run.**

### 8ae. The closure claim and the segmented progression

**The idea, as stated.** The closure of the global boundary board states must be ensured, otherwise the programme fails. A viable
deviation in the choice phase always reaches a boundary state: a route cannot pass from choice to force without hitting one, and if
it is not the same boundary state it is a different global one, by definition. Otherwise a win could be forced from the initial
position, by an opening in the closure that decides chess with no blunder possible. No such line is known to exist in chess theory,
so the closure is predicated on its nonexistence, and enumerating the closure would be the proof. If the programme fails in this way,
the result is a forced opening line, a discovery that reshapes chess. If it succeeds, the forward and backward trace is complete and
at least a weak solution is viable. A missing primitive or a failed belt prediction is not a failure of the programme, because the
catalog is incomplete and the reduction rule develops it. A failure of method is the admission of primitives that do not fit the
architecture. The closure is expected to be achieved through the cyclical structure, or through zugzwang that forces eventual losing
moves and classifies the branch upstream once found.

On the segmented progression: the space between the initial position and the global boundary states is vast, and it is mitigated by
choosing a depth that scopes the primitives. A segment may show that a local candidate was already passed, which is information (no
boundary state exists below it on those paths), and it is cached as a checkpoint. Whether a local candidate is global, or provably
not global, is settled upstream later by the progression of segments, the catalog and the cache. A file for the segmented progression
is to be started.

**Followed, and integrated.** The closure claim is §2a of `boundary_board_states_opening_survey.md` and §14 item 11 of
`boundary_board_states.md`. The segmented progression is recorded in the survey as §7a with only what is stated; the process is not
yet specified and nothing has been filled in beyond that.

**Assessment — separate from the record above.**

- **The relation to §3 and §5.** The existence of a transition on every route to a decisive position follows from a drawn start
  (§5; the collective version is §3). The closure claim adds that each transition is a boundary state in the sense of item 10 and the
  cyclic definition. That is the claim the enumeration tests.
- **The premise is recorded as a premise.** The closure is conditional on the initial position being a draw, which is unproven.
- **Not run.**

### 8af. Pruning as information

**The idea, as stated.** Pruning is informative and does not only remove data. A forcing sequence that wins is a won race condition and an identifiable primitive. Perturbing upstream shows whether it is baked into the structure or emerged from a chained dependency or an ordering of moves. Tracing it shows whether an equilibrium can be reached. The pruned structure carries why each edge exists (the primitive that refutes the alternative), and the balance of the position's equilibrium comes out of the pruning. If no equilibrium can be established upstream, the branch has no boundary state. If one is found, it is a candidate anchored from the initial position.

**Followed, and integrated.** Survey §4a. Each pruned edge is recorded with its refuting primitive, its certificate type and its upstream classification.

**Assessment — separate from the record above.**

- **"No equilibrium can be established" has to be a deduction.** It follows from the joint incompatibility of the primitives (§8x). A failure to compose with the present catalog is a residual interface and may be a missing primitive (§8p).
- **The pruning is the certificate when it is closed.** A found equilibrium is a candidate until the pruned structure is a closed system with a certificate on every refutation. The complete deviation map over that closed region is the complete case list (§8d, §8h), and it certifies the local state. A first reading in discussion treated the case list as a further step after the deviation map. The designer corrected it: they are the same thing, and the reading was withdrawn.
- **Not run.**

### 8ag. Equilibrium maintained until the move pool narrows

**The idea, as stated.** The initial position resides in the choice phase, with the primitive race conditions in an equilibrium from
which a forced draw can emerge and with no forced sequence to a decisive outcome. That it is a draw is assumed because no forced win
has been found, and because no primitive race has committed moves to produce an opening: with the symmetry of the starting position,
any forced attempt can be mitigated by the moves available. Imbalance emerges from decisive moves that create the opportunity for a
race to be won, so a primitive is reached only through the moves that enable it. In Bird's Opening the fool's-mate primitive depends
on White making the enabling moves, and Black cannot force it open, so under optimal play it is not reached. Optimal play exhausts
the moves that mitigate the primitives until the move pool forces an imbalance, as in zugzwang. The move pool is enumerated in the
pruned structure, and its values come from the pruning, from the losing conditions of the primitives, and not from heuristics. The
global boundary states are found by maintaining the equilibrium until then: take the longest line along which the pruned structure
shows the equilibrium maintained, isolate where any move of the pool enters decisive territory, and where the position can be forced
to a draw it is a boundary state if the deviating moves of both sides decisively lose. A decisive forcing line in a side's favor
forces the opponent to avoid the branch, and this interdependent composition narrows the pruned structure.

**Followed, and integrated.** Survey §4b.

**Assessment — separate from the record above.**

- **The drawn start stays a premise** (survey §2a). The enumeration is what would prove it.
- **Not run.**


### 8ah. The process as a whole; temporal analytical relativistic measuring

**The idea, as stated.** Build a pruned structure tree from the opening position onward, focused on a specific opening. Understand
the primitives as they evolve through segmented exhaustive coverage, to map the branch structure and prune it, aiming to map the
choice phase of the choice-force dichotomy, where the boundary board state is the boundary. Perturbation and mining map the pruned
structure; the primitive race conditions are discovered by segmented exhaustive coverage; the structure is pruned on their
certainties. A further pruning mechanism is the designer's coined *temporal analytical relativistic measuring*, which means
something specific: a position A with White to move, known to be to Black's advantage, reached again down the line with Black to
move, has its difference determined from that position and the downstream primitive calculations, including the impact on zugzwang
situations; the same position with the same tempo is cached, as a routing to the structure already mapped, symmetrically the same
because the downstream position is the same. A race condition closes either by a decisive move, which can be pruned, or because it
is forced by the move pool decision; whichever it is informs how to proceed. Pruning a position because it is not a viable path to
the equilibrium is a route. Where the closing rests on move pool reduction and forced calculations, relativistic comparison finds
whether this is the actual equilibrium, or whether none exists, exhaustively. The designer's last sentence was left trailing, and it
is recorded as given.

**Followed, and integrated.** Survey §4c.

**Assessment — separate from the record above.**

- **Nothing in it is new to the classification.** The same-position, same-tempo cache is survey §6 (convergence). The cross-tempo
  comparison is the new element, and how it is recorded is unspecified.
- **Not run.**

### 8ai. Symmetry emerges from the pruned structure

**The idea, as stated.** The symmetry emerges from the structure when it is pruned sufficiently. The pruned portions are decisive
forced positions, in favor or not, with the primitives understood from the pruning. What remains are the positions that are not
decisive, and within forced sequences those where draws are still possible, forced wins are not reachable, and there are only
forced losses on deviation from the drawing choice. That is where the symmetry emerges. If the decisive and forced positions are all
pruned, so that the entire structure is pruned, there is no equilibrium boundary board state position. If something remains, it can
be analyzed.

**Followed, and integrated.** Survey §4d.

**The designer's clarification.** The symmetry is a measurable emergent property of existing at all after the pruning. It can only
exist when there are no forced win moves, a draw is available, and a deviation from the draw moves opens the door to a forced
sequence in the losing favor.

**Assessment — separate from the record above.**

- **Not run.**

### 8aj. The roles of the two sides and the tempo of the initial position

**The idea, as stated.** The initial position has a tempo, and White has it: a bounded structural fact. Assuming no forced win from
the first move, the two phases exist. Under optimal play White is the constructor and Black the reactor. White chooses the opening
line and so the direction toward the boundary board state; Black chooses from the pool White's choices have set. White also chooses
whether to maintain the tempo or give it to Black, a choice Black cannot make; optimal play to a boundary state may not depend on
maintaining it. For discovery, branches that surrender the tempo can be pruned. Black's observed fianchetto schematic (g6, Bg7, d6,
and the b6 and Bb7 variant) is read as the primitives showing in play: development race progress that does not impact the other
primitives negatively and commits to no branching path, with deviating threats in a small pool. The start is therefore the richest
vein of the branch network, not Bird's Opening, constructing against the known Black schematic with the pruning methodology. The
designer regards this as a missing piece of the protective belt.

**Tempo, as defined.** Tempo is an explicit, first-principles calculable metric, not a sense or a heuristic, and an emergent property
of the whole structure. White, first to move from the initial position, has a temporal advantage relative to Black. The balance of
constructor and reactor is anchored at the equilibrium where both sides have the same tempo relative to the primitive race condition
calculations: the boundary board state positions, whose move pool is only draw moves or loss moves for each side. In a forced phase
the side with the tempo forces the outcome by winning on a race condition primitive. In the choice phase it constructs the path toward
the boundary board state with the tempo, and the other side reacts toward a limited pool of boundary board state positions. Tempo is
who is in the driver's seat under optimal play, and whole positions differ as win, draw or loss by whose turn it is, a temporal
differentiation marker. It is determined exhaustively, or from the board position, the side to move and the exhaustive joint
compatible primitive calculations. In the choice phase the calculations are not too large; any primitive race condition for a forcing
sequence can be mitigated through a forced draw under optimal conditions.

**Followed, and integrated.** Survey §4e.

**Assessment — separate from the record above.**

- **Scope.** A prune by the tempo criterion, and the fixed Black schematic, find candidates and do not close the structure (survey §4e).
- **Belt.** Not assessed against the admission rule (catalog §1). Not a catalog entry.
- **Not run.**

### 8ak. Composition, what is cheap, and reachability from certified positions

**The idea, as stated.** Primitive calculation includes identifying properly how the primitives compose. They are piece agnostic and
define race conditions for a goal; how the goals depend on each other in a segmented way, and how the primitives are composed into
the calculation, is another matter entirely. The only cheap things are the board position, whose turn it is, the catalog, and the
certified boundary board state positions. Closure from the initial position to the boundary board state positions gives the weak
solution. Proving those positions is exhaustive coverage from mate, hopefully within the 7-piece Syzygy tablesets. From a certified
position, a move before it must reside within the closure, so the work is determining where the closure is and the boundary line from
the initial position to whatever boundary exists. As the catalog grows the reachability can be exploited, along with other emergent
relationships, to optimize the pruning.

**Followed, and integrated.** Survey §4f.

**Assessment — separate from the record above.**

- **Tablebase scope.** Direct certification applies to positions of seven pieces or fewer; larger states rest on the case list, with
  tablebase leaves (survey §4f).
- **Composition** beyond the three primitives is answered by the designer in survey §4f: joint compatibility equilibrium analysis after
  a boundary board state position is confirmed, adding the primitives that do not cover.
- **Not run.**

### 8al. The order of work

**The idea, as stated.** Scope the combinatorial explosion by finding the furthest downstream branch and certifying its boundary
board state status first, hopefully within Syzygy; then use the table with reachability to enumerate, from the forward trace, the other
reachable positions that are enumerated in Syzygy; scope the difference, do exhaustive coverage from the known Syzygy solutions, and
build up from there. Certified local positions are cataloged first, and global status is certified afterwards. It is a constructive
process for the protective belt and cannot all be planned in advance. The uncertainty is whether downstream solutions can be reached
and cataloged within the Syzygy limitations, whether only local states can be certified, and whether the work can return to global
status. Everything except White's constructor path is, to the designer, resolved.

**Followed, and integrated.** Survey §4h.

**Assessment — separate from the record above.**

- **Consistent with §7.** The scaffolding there starts from confirmed boundary states; §4h adds the forward trace as the source of the
  first one.
- **Not run.**

### 8am. Pruning by primitive race loss: which moves force a deviation

**The idea, as stated.** Take a planned opening framed as an optimal opening for one side, and ask which moves by the other side force a deviation from it by making a move of the plan lose a primitive race condition. Look at what the plan's first move changes, ask which of the opponent's pieces can reach the contested line within the race's move count, and compare the result against the sum of the deviating responses. In the example (Black's g6, Bg7, d6), only pawn moves of the first move remain. The calculation runs from the anchor, which is a constant in every game, and one counter-example beating a race of a plan invalidates it as optimal play.

**Followed, and integrated.** Survey §4i.

**Assessment — separate from the record above.**

- **Not run.**

### 8an. The discovery procedure on any seed position: backtrace, race, perturbation, joint compatibility

**The idea, as stated.** There are four distinct, generalizable properties that can be made into a protocol to survey and discover
local boundary board states, and to work upon them to see whether they can be worked retrograde upstream to fill out the catalog.

1. **Backtracing and reachability under optimality.** If the position is a local boundary board state and is being backtraced, that
   primitive is the focus of the reduction in the first place.
2. **Perturbation.** A calibration change to a measurable primitive race, to see whether an equilibrium is established from an
   imbalanced position, and whether the position arrives closer to an equilibrium position. Perturbation is required to navigate toward
   discovering an equilibrium local boundary board state position.
3. **The nature of the primitive race.** Its forced, sequential character and its segmented dependency, which show how obstruction can
   be possible, and why a race whose delivery is forced beats one whose delivery is not.
4. **Whether joint compatibility of the primitives is possible** (perturbation-related; the designer says it may be the same thing).
   The compatibility must come from the obstruction of the other side's mating race. Maintaining that obstruction is then the only
   winning condition, and committing to the mating race means losing the race that is being obstructed. A race that cannot be
   obstructed without changing the race itself is composed so that it is interdependent with obstructing the other.
   A race that is not jointly compatible with obstructing the other is set aside, and the direction of the perturbation emerges from
   this: toward the obstructions the position needs to be viable (the designer's illustration, start-selection §2).

**It works on any position.** The designer states that the process is generalized and can be done on any position:

- **On a decisively won position** (a position that is won through already decisive primitives), it is used to discover new local
  boundary board state positions: understand the race, perturb to change the race calculations, and work toward an equilibrium.
- **On a certified local boundary board state**, the same mechanism determines the local boundary board states reachable from it. The
  state is retrograde analyzed under optimality, to determine whether it navigates to another boundary state, or whether an upstream
  position forces navigation to this equilibrium, and so on. This grows the catalog.

**Dead pieces.** Pieces that cannot affect the primitive race to mate are dead to the mating net: if they did not exist, the mate is
still as inevitable, and if they exist, the mate is still inevitable. The same can hold in a local boundary board state position,
for pieces that do not impact the primitive race conditions that enforce the draw equilibrium, and positions can be reduced by them.

**Dead pieces, further (the designer's statement).** Adding dead pieces to a position is like adding pieces to a forced sequential win
that do not impact that primitive race: no additional calculation is needed unless they can obstruct the race. This is a feature in
two directions. It expands the catalog, since a dead piece added to a certified local boundary board state gives another position. It
also allows a local boundary board state to be reduced by taking out dead pieces, to find the true primitive of the state and to isolate
the impactful primitive joint composition positions. The deviation lines need coverage only if they beat the race conditions. Where an
obstruction is present, exhaustive coverage of it is required. Where the cost of not following the forced draw cycle is a forced loss
on the primitive under analysis, the exhaustive coverage is scoped to that distance. The designer states that this is deliberately
woven into the program, which is general and constructive.

**With the tablebase.** With fewer than 8 pieces, Syzygy determines win, draw or loss, and the tablebase can be navigated using that
understanding for any of the steps. What the designer's machine can store (perhaps 5 pieces at most) is to be considered later.
Syzygy or opening theory may be searched to identify candidates, or positions close to a candidate. Endgame positions may be less
complex and more calculable.

**The bottleneck and the automation (the designer's statement).** The tablebase is a basis on which the catalog is built, and it only
helps the program. It is limited (strict storage limits on exhaustive coverage, and no forward trace at all) and it does not address
what this program does. The designer has used tablebases to navigate the full landscape before. The automated process depends on
understanding and building the protective belt, and on understanding the process, with the structures being built out as the
automation nodes. Race length (the primitive calculations), obstruction (something that impacts the race condition calculations by
existing, understood through its change under perturbation) and dependency (the interdependent race primitives that compose into a
goal, as in f7 then f8 with a discovered-check forced sequence beating a one-move mate by force, through the segmented dependency
primitive) are understood. The real bottleneck is the classification and structure of how to compose the primitives properly, and
understanding the generalizable principles well enough to construct a script that handles the structure properly, exhaustively and in
general. It is difficult without a program with explicit considerations for it, and the designer states that it is not that it cannot
be done.

**What it does for the program.** The process benefits from a more extensive catalog: local boundary board states can be certified
from local positions upstream, reliably, in a scoped manner. With a full enough catalog, reachability upstream to the initial board
state position can eventually be confirmed. The designer regards this as an important addition to the program that changes how
practical the entire program is.

**Followed, and integrated.** The routes it opens for the first certified local state, and the designer's illustration of the
procedure on a decisive seed, are in `boundary_board_states_start_selection.md` §2. The pointer from the survey is at §4f.

**Assessment — separate from the record above.**

- **Relation to §7.** The scaffolding there starts from confirmed boundary states and extends upstream. This procedure is read as the
  operation that does the extending, and it takes a decisive position as a starting point as well as a confirmed state.
- **Relation to the anchored survey.** The survey walks forward from the initial position. This procedure is seeded by a position,
  not anchored at the initial position. A local boundary board state is decoupled from the initial position and the forward trace by
  definition (survey §2), so a state found this way is local; a state certified as local is not thereby globally certified. The two
  directions meet at the boundary line of the catalog (survey §4f).
- **Dead pieces.** The reading that stood here, that the deviation lines are covered as well as the main line, was replaced by the
  designer's rule above. How a reduced state is entered in the catalog, as its own entry or as a link to the original, is not specified.
- **Composition.** The bottleneck stated above is the part survey §4f records as not specified, and the designer says it is found
  through a first example.
- **Relation to the closure.** The scope of the programme is unchanged: the closure of the global boundary board states from the
  initial position, and the proof of that outcome (survey §2a). This procedure builds the certified catalog and its reachability that
  the closure is proved against, and the population it produces can inform and optimize the rest of the programme. It is a means to the
  closure and does not replace it.
- **The catalog.** Nothing here is added to it. Perturbation and dead pieces are operations on a position and not primitives (catalog §1).
- **Not run.**

### 8ao. The proof object of a local state, and the proof of an upstream position

**The idea, as stated.** The cataloging of a locally certified boundary board state position is fairly simple and straightforward. The
pruning structure already detailed is applied to the boundary board state position being certified, and that structure is the proof
object. Any upstream boundary board state then only has to be the exhaustive coverage between the two proof objects.

**Followed, and integrated.** Survey §4f.

**Assessment — separate from the record above.**

- **Where the pieces were.** The proof object is the closed pruned structure (survey §4 step 6, §4a). The reuse is the position-keyed
  cache (survey §6) and the early-termination role of a certified entry (§8p). The upstream extension is §7 and survey §4f. The
  statement that an upstream position's proof is the coverage between proof objects was not in one place, and the survey now has it.
- **The condition.** For both sides, every branch of the upstream structure ends at a certified forced loss, a certified entry, or a
  non-losing move that stays inside the region (survey §4a), and not only along the path.
- **Not run.**

### 8ap. Two catalogs, perturbation in both settings, and tempo and its anchor

**The idea, as stated.** *Two catalogs.* There are two catalogs: the primitive catalog and a boundary board state certification catalog. Both are important. The primitive
catalog is important for understanding primitive composition, which dictates the calculations. The catalog of boundary board state
positions and reachability families is important to eventually enumerate the closure and to progress from local boundary board states
to families of reachable local boundaries, building up from there.

*Perturbation, in both settings.* Perturbation is the analysis of a board state via move changes, whether by changing piece position entirely or by movement up a
branch through legal move requirements. The change in primitive state caused by a move applies to moving up the tree to find where a
race was acted upon, through the constraint of legal moves. When surveying to discover local boundary board states through
construction, the legal move constraint is not required, since how to arrive at the position is not a concern and there is no
consideration of the initial game start anchor. Perturbation in that context allows a change in board state in piece position without
legal move restriction.

*The scope of a perturbation.* The perturbation scope is dictated by the intended purpose of the perturbation and by whether it is anchored or not in the first
place. Perturbation is to be understood for what it is, a change of state, whether in discovery or in following a branch.

*Tempo and its anchor.* Tempo on a given position is definitionally coupled to the initial anchoring measurement point. If the tempo relation can be
understood from arriving at the position from any route, that tempo is the same calculation regardless of route, by the definition of
tempo. To discover tempo, the calculations of tempo to arrive from the initial position must be discovered. It is definitionally
prohibited to understand tempo at all without this anchor, or without the anchor of the terminal game condition (the equilibrium, or
the forced sequence, anchoring position), to understand the temporal relative value in the first place. The designer adds that the
confusion about this is understandable and should be informed by the program itself.

**Followed, and integrated.** Survey §4j (the catalogs and perturbation) and §4e (tempo); the perturbation convention in the primitive
catalog §6; `boundary_board_states_start_selection.md` §2 (the construction route).

**Assessment — separate from the record above.**

- **What the record held.** The catalog status names the separate catalog of candidate, proven and disproven entries. Perturbation was
  defined by a move (catalog §1, §7.7, §8.4; Method 7, §8y), which is the legal-move setting. The initial-position tempo and the
  tempo terms of the catalog were left unrelated (survey §4e, "The belt").
- **Reachability families** are added to the survey's objects (§3), in the designer's phrase.
- **Dead pieces** are the pieces outside the load-bearing set and the footprint (§8s, §8h); the reading is in survey §4j.
- **Tempo** terms are read as measured to a terminal anchor (survey §4e).
- **Not run.**

### 8aq. Composition beyond the three primitives, and reachability families

**The idea, as stated.** *Composition beyond the three primitives.* Composition beyond the three primitives is understood through joint compatibility equilibrium analysis after the confirmation of a
boundary board state position. Either the primitives cover all the primitives or they do not, and the ones that do not are primitives
that need to be added. There is no uncertainty here.

*Reachability families.* Reachability is simple: if a boundary board state can reach another one, then it is in the same reachable family. If that family is
confirmed joint to the initial start position, that is an important certification, and it essentially classifies all of the reachable
positions in that family as certified global boundary board state positions, when reaching the initial position under optimal
backtrace constraints. This is what makes reachable families important: they are a classification methodology to resolve a family of
reachable positions when that family is eventually anchored with the initial start position.

*From discovery to the initial position.* When discovering a local boundary board state from construction, the constraint of backtrace does not apply: discovery has no
constraint to the initial position anchor. Eventually a reachable family is composed that either concludes it cannot reach the initial
position under the optimal reachability constraint, or concludes that it does. This covers the whole process from discovery to
linking to the initial position. A family must reach the initial game start to prove global, if that family is composed under the
optimality backtrace constraint, which is the requirement for any family created.

*Global, restated. Any upstream position that under optimal conditions leads to a boundary lies in the choice region, and the boundary actually being realized does not change the equilibrium that the boundary position and the upstream positions reside within, if they are part of a reachable family under the optimal backtrace constraint. A family that is reachable from the initial start position under that constraint is an optimally reachable boundary board state family, and every member of it is a global boundary board state position. The initial game start position is itself a boundary board state position: it resides in the choice region, and if it has a move that can lose, it is a boundary. The designer adds that a distinction between the choice region and the boundary drawn by how many drawing moves a position has is a meaningless classification criterion.*

**Followed, and integrated.** Survey §4f (composition), §3 and §5a (reachability families).

**Assessment — separate from the record above.**

- **Composition** is the reduction rule (§8x) run on a confirmed state: the residual interface of the composed primitives is the
  specification of a missing one, admitted only under the admission rule (catalog §1).
- **Families** are read with survey §5 and §7 here (survey §5a); the readings are offered for correction.
- **Not run.**

### 8ar. Mining the Syzygy tablebases for the certification catalog

**The idea, as stated.** *The catalog from proven results.* The designer's aim is to begin the boundary board state certification catalog by creating it from
the proven results of the Syzygy tablebase solution set, enumerating the boundary board states from it. The population of draws from
Syzygy is important. A branch can be enumerated by taking the tablebase's win, draw or loss value and distance and mapping it onto the actual branch
with its legal move pool.

*Optimal play, and what is left to later.* Branches are created under optimal constraints. The primitive calculations are not part of the script; they are done when a
drawn position is verified later, and primitive compositional analysis is done after the fact, from the population of certified local
boundary board state positions obtained by exhaustive coverage from Syzygy. The 50-move rule is ignored for now: it matters only in a
reachable family, when the anchor is reached, as a cutoff.

*How a draw is realized.* The draws are categorized by stalemate, insufficient material, or cyclical. Where the verified best move is to reduce to insufficient
material, that is a provable equilibrium of reduction to insufficient material rather than a stalemate choice, since optimality is the concern.

*The boundary and the two spaces.* A position in the choice region that has a move that can lose is a boundary board state position, and the initial game start position is one.
A distinction of the choice region from the boundary by how many drawing moves a position has is a meaningless classification
criterion. There is a literal boundary between force and choice, and a choice space and a forced space that together encapsulate the
entire space of chess, with sub-optimal areas pruned; the pruned structure is the full solution of chess. The designer wants to
confirm that the boundary board state understanding covers a majority of the entire space.

*The scale and the order.* The tablebase size does not matter eventually, for reachability, with the program's methodology of backtracing and exhaustive enumeration
from the cached database. The work is to mine all the Syzygy data the machine can handle, and then to automate off that basis, from
the catalog, with the pruning mechanism: mining first, then the catalog and reachable families from it, then the automation and the
matching of primitives. The designer has said that another angle is under consideration which may change this; none is recorded.

**Followed, and integrated.** The process is written as its own document, `boundary_board_states_syzygy_mining.md`. The statements are summarized in survey §4k, the
classification of a family and the global definition are in survey §5a, and the later specification of the process is in
`boundary_board_states_start_selection.md` §2a. Changelog entries 62 and 65 of `boundary_board_states.md`.

**Assessment — separate from the record above.**

- **The tablebase certifies the closure condition of survey §4a** for a drawn position of a held material: every legal move keeps the
  draw or loses, with the value of each from the table. That is a reading, offered for correction. The test for a losing deviation of
  each side at the same position is not fixed by the statements.
- **The term *interior*** (drawn positions with no losing move) is this file's, used in the process document to name the remainder of the
  choice space. The designer did not name it.
- **Not run.** The process is a specification. Nothing in it is carried over from any earlier test, and the designer will run and verify the data on return. Primitives are placed after the fact, as the designer states.
- **Scale.** The sizes of larger materials are estimates in the process document (§8), not measurements.

---

## 9. What remains open




- **Whether any placement exists that satisfies all of §8d's conditions at once**
  — bishop color-lock, a tempo margin for every route that breaks the structure,
  mutual blocking, capture-diagonal integrity, waiting moves for both sides, and
  no stalemate. This is the discovery itself, not a formality.
- **Whether a construction of this kind can yield an exactly-one boundary state
  or only a proven wide fortress (§8d)** — the two pull in opposite directions
  and the target has to be chosen.
- **Whether the catalog's primitives (§8b) keep their deciding inequalities
  when composed on one board** — untested; each composition has to be re-argued,
  not assumed.
- **Whether the transition into an opposite-colored-bishop blockade contains a
  genuinely narrow, tempo-critical moment (§8d)** — a hypothesis, no data.
- **Whether Method 3's population and cone-shape analysis
  (`boundary_board_states.md` §2, and this document's own §7) finds any real,
  usable correlation at all, once enough confirmed boundary states exist to
  test it against** — a genuine empirical question with no answer yet either
  way; the method is well-scoped as a prioritization tool regardless of the
  outcome, but whether it actually outperforms undirected candidate generation
  is untested.
- **Whether the scoped, topology-mapping version of the checkers precedent
  (§7) genuinely avoids chess's own scale wall, or merely hits the same wall
  in a different guise once pursued far enough** — the reframing away from
  full-game solving is real and precisely what sidesteps the objection already
  on record against the naive version, but whether it stays tractable deep
  into the game, rather than just near the already-solved floor, is untested.
- **Whether real, currently-confirmed boundary states actually connect to each
  other within reach of the scaffolding strategy (§7), or whether the
  connected regions it produces stay small and isolated in practice** — the
  strategy is sound regardless of the answer, but how much real structure it
  reveals versus how quickly it hits unresolved territory is untested.
- **Whether any confirmed boundary state will ever be checkable for
  start-reachability (§5) given the scale gap**, or whether this remains a
  well-defined but practically unanswerable question for the foreseeable
  future — the concept doesn't depend on this being resolved, but its
  usefulness as a research target does.
- Whether the collective closure argument (§3) is worth building as an actual
  running check — systematically retrograde-tracing decisive positions and
  flagging any chain that doesn't touch the current catalog — or whether it's
  better used informally, spot-checked by hand as the catalog grows. Not yet
  decided, and the practical cost of running it systematically at real scale
  is untested.
- Whether tracking load-bearing/redundant relationships (§4) at real scale
  changes verification order enough to matter in practice, or whether the
  effect is real but marginal — untested.
- Whether scoped coverage certificates (§6) are worth building as an explicit,
  separate data structure alongside the boundary-state cache in
  `boundary_board_states.md` §3, or whether they're better derived on demand
  from what's already cached — an engineering choice, not yet made.
- Whether dominator-style analysis (§4) has a more efficient algorithmic form
  specific to this problem's structure, beyond directly adapting standard
  directed-graph dominator algorithms — unexplored.
- **What "local versus global" is meant to distinguish (§5a)** — extensibility,
  upward connectivity, or both; and, if closure contribution is the intended
  criterion, what narrower object "the closure of the choice/force boundary"
  denotes, since every genuine choice-point boundary state is already a
  transition point in the §3 sense.
- **Whether extensibility can be defined operationally (§5a)** — for example as
  the set of single-piece additions under which the value and the exactly-one
  property survive — and whether it correlates with the load-bearing relation
  of §4.
- **Whether duty-based (non-frozen) locks can be certified by construction
  (§8f, §8h)** — that is, whether the four footprint conditions can be met for
  any real fortress-style structure without the footprint growing to the whole
  board.
- **Whether higher-order joint compatibility (§8g) has a compact
  representation** — a constraint system over shared resources, beyond pairwise
  checks — and whether the tempo-induction argument suffices to ground circular
  dependence in real compositions.
- **Whether the composition rules of §8i survive validation** on small
  compositions, including the load-bearing ablation test, before any large
  construction is trusted. Untested.
- **Whether complete interdependence (§8g) characterizes tight boundary states
  in practice**, or whether redundant-structure boundary states are common
  enough that it is only one family among several.
- **Which primitive to write first as a contract (§8j)** — the rook-pawn corner
  race is the proposed candidate; not yet attempted.
- **Whether "materially live" (§5a-iii) is the right operational form of
  "sufficient material,"** and whether the local/global distinction by
  liveness supersedes the closure-contribution wording of §5a-i.
- **Whether hanging-piece positions in total-win materials count as local
  boundary states (§5a-iii)** or fall outside because they involve no locking
  mechanism.
- **Whether the board edge's role decomposes into restriction plus funneling
  (§8k)**, whether a substitute must supply both, and whether exactly-one
  positions in tablebase KPvK correlate with edge adjacency.
- **Whether wall pieces can replace the edge (§8k)** without introducing
  capture, check or stalemate interference that breaks the primitive's own
  race.
- **The start-position compensation hypothesis, stated and unrecorded until
  now.** As stated: from the symmetric start, any material or positional
  asymmetry along optimal play is compensated by topological conditions
  (locking mechanisms); otherwise it would be a plain advantage. Assessment:
  if the start is a draw, some compensation must exist on every optimal line
  (§5), so the weak form is a near-necessity; the strong form, that every
  compensation reduces to composable race-certified locks, is an empirical
  hypothesis, because dynamic compensation (initiative, perpetual check) and
  temporary compensation (trading down before any lock forms) are other
  candidate mechanisms. Test: take drawn positions with a material imbalance in
  the tablebase range (for example KRPvKR, KQvKR fortresses, KNvKP), classify
  the drawing mechanism of each, and see whether nearly all fall into
  primitive types. Not run. The connection to the start position itself is the
  start-reachability problem of §5 and is not checkable at scale.
- **Whether minimization by ablation (§8l) yields a usable decomposition
  operator,** given that fragments need not be canonical, and whether the
  obstacle-board solver it requires can be built to handle wall squares.
- **Whether boundary states in real tablebases decompose into a modest,
  recurring set of primitives (§8l),** which decides whether the catalog
  compresses.
- **Whether lifting (§8l) works for the rook-pawn primitive:** whether locked
  structure can reproduce the edge's exclusion and its mobility exclusivity
  without introducing capture, check or stalemate interference.
- **Whether a population scan turns up enough confirmed local boundary states,
  across enough distinct materials, for recurrence-finding to have real
  content (§8m)** — the population-first approach applies in spirit
  regardless, but if rook-pawn edge draws are close to the only clean
  small-material instance, step 3 of §8m's corrected order has little to work
  with until larger materials are reachable, and the single-example fallback
  in §8j becomes the practical path. Untested, and worth checking early since
  it decides whether the scan is worth building before the ablation machinery
  is.
- **Whether a composed floor of real, useful size can actually be built
  (§7a)** — restates §8m's own open question in the specific terms of the
  checkers-style scaffolding strategy: growing the floor past the piece-count
  wall only matters if enough compositions can actually be certified.
- **Whether the generalized race in Uri's rule reduces exactly to the
  published "rule of the square" with a second racer's distance added
  (§8n)** — a direct comparison of two existing formulas, not yet carried
  out.
- **Whether other existing named endgame heuristics decompose the same way
  Uri's rule appears to (§8n)** — knight-with-a-rook-pawn-on-the-seventh as
  the next candidate, per §8b's own entry list.
- **Whether Bleicher's "Freezer" tool (§7a) or similar a-priori-constrained
  tablebase constructions have already been applied to any of this project's
  own primitive candidates**, which would be closer prior art than this
  session's brief check found.
- **Whether Method 6 (perturbation search, §8o), run starting from positions
  already close to known drawing material rather than from exposed
  middlegame attacking positions, lands on genuine drawn boundary-state
  candidates more often than it lands on decisive-to-decisive seams.** Not
  yet tried in that form.
- **Whether "basin structure" around a boundary state versus around a
  decisive-to-decisive flip point (§8o) can be defined precisely enough to
  compare at all**, and whether any such comparison would say something an
  individual case-by-case certification doesn't already say.
- **Whether an exact, complete correspondence network — analogous to
  corresponding-squares theory, but for this project's broader primitive
  types (§8p) — can be derived at all**, given that the classical version is
  exact only for a narrow, highly constrained class of king-and-pawn
  structures. This is §8m's recurrence question, asked again in the specific
  terms §8p needs: not just whether primitives recur, but whether their
  relationships to each other can be made exact rather than merely
  engine-estimated.
- **Whether failure-to-reach-a-known-anchor can ever be upgraded from "the
  catalog doesn't classify this" to a genuine decisiveness proof (§8p)**
  without first proving the catalog complete for the relevant structure —
  and whether "complete for the relevant structure" is itself a checkable
  condition short of exhaustive tablebase coverage.
- **Whether a formal "refutation certificate" (§8q) — reusing the footprint
  and interference-closure machinery of §8h, but targeting dominance of one
  move over another rather than a position's drawn value — is worth
  formalizing as its own documented proof type**, distinct from but
  structurally related to a boundary-state certificate. Not yet attempted;
  flagged because the match found this session suggests it would reuse
  existing machinery rather than need new.
- **Whether any Method 5 composition proposed so far (§8j's rook-pawn
  contract, the locked-pawn-pair substitutes of §8k) actually clears the
  Tier 3 start-reachability bar now required of it (§8q)**, or whether this
  requirement, once checked, rules out some of what's already been proposed
  as a composition target.
- **Whether a forward perturbation search over many plies (§8r) can be
  defined at all**, since the only perturbation procedure in the record is
  Method 6's single-parameter probe, and forward trace over the full tree is
  the search perturbation is meant to avoid.
- **Whether the gap DTM − CMD (§8r) says anything about defender resources or
  about boundary states**, and whether CMD is small and finite at confirmed
  boundary states as the §8r handle predicts. Both are untested and are this
  assessment's own proposals.
- **Whether the per-ply race quantities proposed in §8r (depth-limited CMD,
  undeveloped-piece count, queen safe squares) track Stockfish's evaluation
  changes in the OTB game**, with the prediction stated in advance. Not run.
- **Whether the validation suite in §8r reproduces each known result exactly**
  from a primitive decomposition, and which of them (rule of the square,
  opposition, Philidor, Lucena, Vančura, triangulation) can be given an exact
  rather than engine-estimated offset.
- **Whether the ablation filter (§8s) collapses the 347 White fool's-mate-length
  sequences to one skeleton, and what it does to cooperative mates at greater
  depth.** Run on three final positions only. Depth is the constraint: full
  enumeration to 6 plies is on the order of 10^8 positions.
- **Whether the mechanism catalog of §8s (own-piece lock, attacker coverage,
  support, masking, out-of-play defenders, absolute pin, decoy, forced
  enabling move) accounts for the named mating patterns** (back-rank mate,
  smothered mate, the Greek-gift mate on h7, Anastasia's mate), as a
  validation-suite extension. Not run.
- **The defender-side counterpart (§8s):** for a confirmed boundary state, which
  term of the opponent's possible mate or promotion does the single saving
  move keep unsatisfied, and does the answer recur across boundary states.
  Untested, and it is the step that connects the mate templates to draws.
- **Whether the 3,960 narrow king-move positions of §8t decompose, by
  ablation, into a small recurring set of primitives**, and whether the
  decomposition is the opposition and key-square structure of classical theory.
  This is the first real test of §8m's recurrence question. The positions are
  highly correlated, so recurrence must be measured over distinct structures,
  not over the raw count.
- **Whether the rook-pawn draw behaves as a wide edge primitive in composition
  (§8t, §8k)**, given that it has the widest draws of any file. If it does, it is
  the right terrain to lift the edge condition into pieces, and the wrong
  first candidate for a narrow boundary state.
- **Net recognition before the mate (§8u).** Whether the four-term predicate, with
  tempo costs attached to each mechanism, can detect from a position that has not
  yet been mated that the terms can be completed by forcing moves, and how that
  compares with the engine on a set of puzzles. Not attempted; the checks so far
  are on final mating positions only.
- **Whether the mechanism list of §8u (decoy fill, pin of the interposer, masking
  by the defender's own piece, pawn-protected checker, the attacker's shield)
  accounts for a larger corpus of mating positions**, using
  `mating_net_report.py`. Four mates have been run.
- **The decision flagged in the OTB case study §7** on whether a single move that
  discharges several terms without conflict is a named pattern; there are now
  two instances (f5, g4), of different kinds.
- **The cost model of P-001** (catalog §4.8): how tempo cost attaches to each
  mechanism, which is the unfinished part of the entry.
- **Whether the evaluation procedure of catalog §4.7, applied to non-terminal
  positions, gives a usable DORMANT or ENGAGED signal.** A specification only; not run.
- **The writing order of the primitives (§8v, §8ab).** P-001 is written; P-003 and P-008 are SPEC and not run. P-002, P-004, P-005 and P-007 are not primitives (§8ab), so they are not to be written. What the three leave unaccounted for, such as zugzwang without capture, is returned by the reduction rule as a residual if a proven entry needs it.
- **Resolved:** exactly-one is not required at any node of the equilibrium loop
  (`boundary_board_states.md` §14 item 8).
- **Which measure of an equilibrium is informative:** the width (non-losing moves at a
  node), the number of tied primitives, or the number of losing deviations available to
  each side. Width is cheap to compute; the others need the catalog. Not decided.
- **What the catalog needs in order to derive the outcome clause and the counter-net comparison for
  non-royal targets (§8w, catalog §4.3 and §4.10):** the cost model, the tempo comparison (terms TP and TM), the other
  primitives, and a well-founded ordering.
- **Whether the tempo ledger of §8w can be built over a whole game** and whether its
  columns track the engine's evaluation changes, stated in advance as in §8r. Not run.
- **Whether escape steering (§8w) recurs**: a residual escape of one net falling into a
  mechanism of another. One instance so far.
- **Which proven position to record first as a proven boundary-state entry (§8x)**, and
  whether the deviation map is the right form of the composition record where a refutation is
  a sequence of races and not a single net.
- **Exhaustive insurance on the first composed entry (§8x)**, to confirm the composition rules.
  It is run when the first entry is composed.
- **A material in which the defender also has pawn moves (§8t)**, so that the
  claim that reversible moves are required to hold an equilibrium can be tested.
  King and pawn against king cannot discriminate. A self-built solver for
  king and pawn against king and pawn is the next material up and is much
  larger.
- **Whether the finding that draws are narrowest on central files and widest on
  the rook-pawn file holds in other materials**, or is a feature of king and
  pawn against king alone.
- **The form of the constraint system of the reduction rule (§8x)**: whether the joint-compatibility
  constraints fix the missing primitive uniquely or a family, and how the proven entry and its
  neighbors select.
- **The unit of the activation gain, and the definition of the primitive state S and the asymmetry
  (§8y, catalog §7.6).** The measure of distance from a state to an equilibrium composition needs the
  cost model.
- **Resolved (§8ab):** there is no boundary between P-002 and P-007 to draw. The tempo comparison is a term of the chained-dependencies primitive and of the coverage race, and a promotion is the resolution of a chain by the route RT2 (catalog §7).
- **Whether castling is an activation or a separate lock removal (catalog §7.14).**
- **Whether the dependency graph of §8z, read from a position by P-001's mechanisms, reproduces as its cut sets the load-bearing
  sets that ablation finds (§8s).** Prediction P3-1 of catalog §7.13, stated in advance. Not run.
- **Whether orbit exclusion (catalog RT3) reproduces the known color-lock verdicts across the whole class, and whether it is the
  color-lock fact inside Uri's rule (§8n).** Prediction P3-2 of catalog §7.13. Not run, and the identification with §8n is a reading
  to be checked against that section.
- **Resolved (§8ab):** the movement relation and the dependency graph are part of the chained-dependencies primitive P-003 (catalog §7), not a substrate, and P-003 decides a result of its own: the race to resolve a chain.
- **How permanence of coverage (§8z, catalog §7.4) extends to the region where a relation changes**, since promotion adds
  coverage, and how pins and x-rays enter the dependency sets of route edges.
- **Whether the chain reading of the outcome clause (§8z, catalog O3) reproduces static exchange evaluation where that is exact
  and corrects it where a chain cuts across the exchange.** Not run.
- **Whether the single move that discharges several duties (case study §7, §8u) is an identifiable pattern in the dependency
  graph (§8z):** a move that adds or satisfies two edges at once. Still undecided, and noted as a reading only.
- **How the effect of nested nets under perturbation (§8aa, catalog §8) combines with the tempo comparison and the cost model to give the composed
  race map of the outcome clause (catalog §4.3, §4.13).** Until then P-001 evaluates the outcome clause as it defines it.
- **Whether the key elements of the held terms (the singleton opening sets of catalog §8.3) coincide with the load-bearing elements
  found by ablation (§8s).** Prediction P8-1, stated in advance. Not run.
- **Whether the evaluation of the coverage race under perturbation, composed over a line, reproduces the verdict of exhaustive
  search (§8aa).** Prediction P8-2, stated in advance. A disagreement is a residual interface. Not run.
- **Whether every verdict-changing perturbation in the ablation and deviation data (§8s) meets the support closure (catalog §8.9).**
  Prediction P8-3, stated in advance. A failure specifies by the reduction rule what the catalog lacks. Not run.
- **How many perturbations must be branched (§8aa, catalog §8.9), especially where kings are involved (§8h condition 3), and whether
  additional primitives reduce it.** Kept open by design.
- **Whether the carrier is a term of the coverage race or an input it reads from P-001 (§8aa).** Open.
- **Resolved (§8ab):** P-003 remains a primitive, as the single chained-dependencies primitive (catalog §7, §10.5).
- **How the repetition cycle that defines a boundary state keeps every entry of both registers held across the cycle (catalog §8.14, §5 item 2).**
- **The anchored survey (§8ac, `boundary_board_states_opening_survey.md` §10).** Whether the forced-region cache claims anything about deviations from a forced win; how a web of interdependent candidates is verified; and how the cache treats the repetition cycle of a candidate. Not run.
