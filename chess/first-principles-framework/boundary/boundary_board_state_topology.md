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

**One further precision, tied to an already-open question (§8, and
`boundary_board_states.md` §12).** The unconditional guarantee concerns
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

---

## 8. What remains open

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
