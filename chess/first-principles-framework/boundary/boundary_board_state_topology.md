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
  one leaving every lock in place.
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

### 8h. Scope of verification: footprints and interference closure

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
  check whether exactly-one positions concentrate near the edge.
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
