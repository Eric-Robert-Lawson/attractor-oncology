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
section, §12), and the corrected, complete version below. Kept separate from
`syzygy_improvement_proven_results.md` and the other Claim A/B/C documents because
this is a genuinely different kind of object: those are about distance and escape
count; this is about win/draw/loss pattern-matching at the single precise point
where a losing side's last real choice exists.

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

## 2. Candidate generation — two methods, complementary, both feeding the same verification step

This is the pragmatic, buildable core of the proposal. Two distinct ways to
generate boundary-state candidates are recorded here — a structural filter over
a material's whole space, and a targeted, seeded search starting from positions
already known to be decisive. Both produce the same thing: a candidate needing
the same exactly-one-non-losing-move check (§1, condition 2) before it counts as
a confirmed boundary state, not a replacement for that check.

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
3. **Extend both methods to this project's own solved materials** (KQvK, KBNvK,
   KRPvK once complete) the same way `analyze_longest_line.py` already scans for
   forced plies.
4. **Compare the resulting KBNvK boundary-state count against Haworth's
   already-published zero-mutual-zugzwang result** (§9) — the first real,
   falsifiable sanity check, checkable the moment steps 1-3 produce real numbers.
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
