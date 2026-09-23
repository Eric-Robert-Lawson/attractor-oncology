# Boundary Board States: The Choice-Force Boundary of Perfect Play

**This is recorded as the primary objective of this project, not a secondary
analysis — stated explicitly and preserved as a fact about the project's own
priorities, not something this document should soften.** A separate, honestly
reasoned assessment of this idea's actual scope and limits is kept in §12,
clearly labeled as a distinct, independently-argued position rather than folded
into or diluting the objective stated here.

A record of a concept developed across this project's own conversation history,
preserved here precisely because it went through real revision — an initial
framing that didn't hold up, a correction, a genuine gap in an earlier draft of
this very document (the candidate-list/incremental-verification workflow in §3,
first flattened into "ruled out" when only one of its component claims actually
was), a second gap (the scoping mechanism in §10 was missing entirely from the
prior draft), and the corrected, complete version below. Kept separate from
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

## 2. The candidate-list construction and pruning strategy, in full

This is the pragmatic, buildable core of the proposal.

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
2. **Extend to this project's own solved materials** (KQvK, KBNvK, KRPvK once
   complete) the same way `analyze_longest_line.py` already scans for forced
   plies.
3. **Compare the resulting KBNvK boundary-state count against Haworth's
   already-published zero-mutual-zugzwang result** (§9) — the first real,
   falsifiable sanity check, checkable the moment step 2 produces real numbers.
4. **Build the repurposed forward-search tool** (§10) — a new, dedicated
   pipeline, not the current `compositional_trajectory_solver_modular_shapes.cpp`
   or `run_full_sweep.py` run unmodified against a different output target.
   Reuses their memoization pattern and delta-checkpoint approach (§3), extended
   with boundary-catalog early termination during forward search from a
   specified starting position.
5. **Cache proven boundary-state verdicts to a small, dedicated file**, following
   §3's point 3 — persisting only the proven-verdict subset, not the whole
   landscape, using the same compute-fully/persist-selectively pattern the
   engine's own `pending_export` mechanism already demonstrates is sound.
6. **Test the scoping mechanism (§10) on a material this project can already
   fully classify**, comparing forward-search-with-early-termination against
   full classification directly — the first real, measured data point on whether
   the approach actually reaches further than full classification does, before
   trusting it on a material too large to check the answer against.
7. **For a material neither Syzygy nor this project can fully classify**: §10's
   mechanism, not full classification (§4b), is the candidate path forward — with
   step 6's own result as the evidence for whether it's real.

---

## 12. What remains open

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
- Whether `board.legal_moves` count alone soundly separates "genuine choice-point
  boundary state" from "only one legal move existed, and it happens to draw"
  (§8) — proposed, not yet tested against real data.
- Whether a boundary-state catalog, once built for a real material, reveals
  further structure of its own (families, clustering, anything analogous to the
  shape graph) — genuinely unexplored.
- Whether the insufficient-material ceiling check (§4d) is worth adding as a
  pre-generation material-category skip, rather than only the per-position
  terminal check it already is.
