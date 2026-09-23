# Boundary Board States: A Second Analytical Target Alongside Full-Landscape Classification

A record of a concept developed across this project's own conversation history,
preserved here precisely because it went through real revision — an initial framing
that didn't hold up, a correction, and a final, sound version worth building toward
with certainty rather than re-deriving from memory later. Kept separate from
`syzygy_improvement_proven_results.md` and the other Claim A/B/C documents because
this is a genuinely different kind of object: those are about distance and escape
count; this is about win/draw/loss pattern-matching at the single precise point
where a losing side's last real choice exists.

---

## 1. The concept, in its final, correct form

A **boundary board state** is a position satisfying three conditions simultaneously:

1. **Both sides have sufficient material to theoretically force a win** — this is a
   material-only ceiling check, not a position-specific one (see §4 for why this
   distinction matters enormously).
2. **The side to move has exactly one non-losing legal move.** Every other legal
   move loses. The one non-losing move preserves a draw — however that draw is
   ultimately realized (see §6 for the taxonomy of realization mechanisms, and why
   this project has chosen not to care which one applies).
3. **This is a genuine consequence of correct upstream play, not an artifact.**
   If a forced win were available anywhere upstream, it would have been taken
   instead of drifting toward this boundary — the boundary position is what
   perfect play from an already-non-winning position necessarily produces, not a
   special case requiring separate justification.

The set of all such positions, for a given material, is what this document calls
the **boundary board state set** — proven members, proven non-members, and
(for materials not yet fully resolved) a remaining candidate pool.

**What a boundary board state is not, stated plainly because this was a real point
of confusion worth recording**: it is not a property of *how* a position was
reached, and it does not require the position to have arisen from optimal play
upstream in some real game. A boundary board state is a pure, self-contained
property of the position itself — exactly one non-losing move exists — the same
way this project has always insisted a full landscape must be correct for every
legal position, not just ones a perfect game would produce (see
`general_workflow_reference.md`'s own scope statement). Condition 3 above describes
*why* boundary states are the natural output of exhaustive backward induction, not
a requirement that a boundary state only "counts" if reached via perfect play.

---

## 2. Relationship to this project's existing architecture — same mechanism, different organizing question

This is not a new algorithm. `classify()`'s own fixed-point backward induction
already computes, for every position, the exact value of every legal child. A
boundary board state is simply a specific, checkable *pattern* over values already
being computed — "exactly one child is non-losing" — the same shape of thing as
`analyze_longest_line.py`'s `NumTiedAlternativesThisPly == 1` check for forced
plies, just applied to the win/draw/loss dimension instead of the distance/escape
dimension. Nothing here requires a new solving engine. It requires a new *scan*
over data that either already exists (Syzygy, or a material this project has
already solved) or will exist once full classification finishes (a material this
project is solving for the first time).

---

## 3. What was tried and ruled out — kept in the record, not smoothed over

Two ideas were proposed and tested against real argument before the version in §1
was reached. Both are worth keeping, because the reasoning that ruled them out is
reusable.

### 3a. "Prune during construction using the boundary-state definition itself" — ruled out

The first framing treated boundary-state identification as a way to *avoid* full
exhaustive construction — skip branches that provably can't contain a boundary
state, using only local or material-based checks. This does not hold, for a
precise reason: whether a *specific arrangement* satisfies condition 2 depends on
the true minimax value of every one of its children, which is exactly the
expensive thing full classification computes. There is no local, position-level
shortcut to that value — the same reason there's no way to verify a position is
drawn without first knowing every reply doesn't lose, which was the original,
correct objection raised early in this discussion and never actually overturned,
only refined in what counts as "local."

### 3b. "No check potential" as a sound local filter — ruled out, with a concrete counterexample

A refined version proposed filtering out positions with no check activity as
incapable of being boundary states. This fails on real chess theory, not just in
principle: **king-and-pawn opposition and zugzwang** are exactly this shape —
one side's move loses (forced to cede a key square), the other's doesn't — occurring
in positions with zero check activity anywhere nearby. Opposition is a central,
well-known organizing concept of endgame theory precisely because the deciding
factor is tempo and king geometry, not tactics. A "no check potential" filter
would discard some of the most textbook examples of boundary states that exist.
The sound version of this filter — "no check is ever possible from this position
under any continuation" — is a true global structural fact, but establishing it
requires knowing the reachable future of the position, which costs exactly what
the filter was meant to avoid.

### 3c. "Cache less (boundary states only), recompute more" — ruled out by this project's own prior, directly relevant experiment

A later framing accepted the real cost of full computation but proposed storing
only the boundary states themselves, recomputing subtrees on demand rather than
caching the full landscape, to keep the persistent file small. This runs directly
into an experiment this project already ran and measured: caching *more*
(`reachable_shapes` by canonical form, and a shared cross-pass cache of move
tuples) to avoid recomputation was tested and came back **~6.5% slower, twice**,
because `explore()` still walks the full subtree regardless of cache hits —
**traversal, not value storage, is the dominant cost** in this codebase. Caching
*less* is the same tradeoff run in the opposite direction, and hits the identical,
already-diagnosed wall: a scheme that increases how much re-traversal is needed
makes the already-dominant cost worse, not better.

**What does survive from this whole line of thinking, real and narrow**: material
that is provably insufficient to force a win *regardless of position* — KvK,
KNvK, KBvK — is already recognized by this project's own engine
(`compositional_trajectory_solver_modular_shapes.cpp`, confirmed directly by its
own comments on unconditional-draw terminals). This is a genuine ceiling argument,
not a floor argument about specific arrangements within a material that can
support real wins — which is exactly the distinction that makes it sound where
§3a and §3b were not.

---

## 4. Where this actually gets easy — and it's easier than anyone in this conversation initially gave it credit for

The corrected, sound version of this idea doesn't try to avoid computation. It
recognizes that **for any material Syzygy already covers (up to 7 pieces), the
computation has already been done, by someone else, in full.** Finding boundary
states in Syzygy-covered material requires no tablebase construction at all —
just legal move generation plus WDL probes, both standard, already-solved
problems:

```python
import chess
import chess.syzygy

def is_boundary_state(board, tablebase):
    """
    Returns True if the side to move has exactly one non-losing legal move
    and every other legal move loses -- the pattern this whole document is
    about, checked directly against already-computed Syzygy WDL values.

    See §5 for the collapse-cursed/blessed choice this makes, and §6 for
    what "non-losing" is deliberately left agnostic about (mechanism).
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

The same logic applies, with zero conceptual change, to any material *this*
project has already fully solved — a scan over the completed database, structurally
identical to how `analyze_longest_line.py` already finds forced plies, just keyed
on W/D/L pattern rather than tied-move count. The only place real, heavy
computation is still required is a material neither Syzygy nor this project has
solved yet — and that's not a new cost this idea introduces, it's the same cost
already being paid for KRPvK right now, for entirely separate reasons.

---

## 5. The 50-move rule — verified directly against Syzygy's own documented WDL semantics, not assumed

This needed checking rather than guessing, and the real answer is more precise
than either possibility considered before checking.

**Syzygy's `probe_wdl` returns five values, not three**: `-2` (unconditional
loss), `-1` (**blessed loss** — mate can be forced against the side to move, but
the 50-move rule saves them), `0` (draw), `1` (**cursed win** — mate can be
forced, but the 50-move rule would turn it into a draw), `2` (unconditional win).
So Syzygy *does* account for the 50-move rule, structurally — the cursed-win and
blessed-loss categories exist specifically to flag "this position's true DTM value
exceeds what the 50-move rule allows."

**But it does so only relative to a hypothetically fresh clock at the probed
position, confirmed directly from the documentation's own wording**: the
unconditional win/loss values are stated as holding "assuming 50-move counter is
zero." `probe_wdl` does not read the real game's actual current halfmove clock —
it always evaluates as if the position were reached immediately after a capture
or pawn push. So a "cursed win" reported by Syzygy means "if the clock reset here,
forcing progress would still take too long" — not "given this exact real game's
current clock value, is a draw actually guaranteed." Those are different
questions, and Syzygy only answers the first one directly.

**The reconciling fact, and it resolves the open question cleanly**: this
project's own engine tracks no move history and no halfmove clock at all — it
computes pure distance-to-mate, exactly matching Syzygy's "assume a fresh clock"
semantics, not the real-game-clock-aware version. So for consistency with this
project's own established classification (which has never modeled the 50-move
rule, by explicit prior choice — see `general_workflow_reference.md`), the correct
practice when bootstrapping from Syzygy is to **collapse cursed-win into
unconditional win, and blessed-loss into unconditional loss** — treating `1` as
`2` and `-1` as `-2` for the purposes of `is_boundary_state()` above. This is not
a workaround; it's the direct, correct alignment between two systems that already
agree on the underlying question and differ only in whether they separately flag
the 50-move edge case. The `<= -1` / `>= 1` thresholds in the code above already
implement this collapse.

Whether the 50-move rule is "irrelevant" or "geometrically avoided" was the
original open question — the precise answer is neither: it's a real, structurally
tracked distinction in Syzygy's own data, and this project's deliberate choice not
to model move history at all is what makes it safe to collapse away rather than
something to reason further about.

---

## 6. Draw-mechanism taxonomy — what's in scope, what needs care, and why

Standard chess rules recognize exactly four ways a position resolves as a draw:
**stalemate**, **insufficient material**, **repetition** (threefold, claimable, or
fivefold, automatic), and the **50/75-move rule**. Nothing else exists under FIDE
rules as a position-level drawing mechanism (draw by agreement is a human choice,
not a property of any position, and irrelevant to this analysis).

- **Repetition and the 50-move rule**: explicitly not distinguished in this
  project's own boundary-state work, by deliberate choice — both represent "the
  position is drawn," and which specific rule a real game would invoke to claim it
  doesn't change whether the position qualifies as a boundary state.
- **Insufficient material**: already the real, narrow ceiling exclusion from §3c —
  a material-only fact, already recognized by this project's engine, not a
  position-specific mechanism requiring case-by-case handling.
- **Stalemate**: the one mechanism worth handling with actual care rather than
  folding in automatically. Stalemate as the *realization* of a boundary state's
  drawing branch — the non-losing move leaves the *opponent* stalemated on their
  reply — is a real, legitimate defensive resource with genuine precedent in
  endgame theory (well-known "stalemate tricks" in king-and-pawn endings), not
  something to exclude by default. What *should* be flagged, not silently
  conflated with genuine boundary states: a "boundary state" whose only
  distinguishing feature is that literally no other legal move existed at all
  (a heavily constrained position where the single legal move happens to be
  drawing, rather than a genuine choice among multiple real alternatives where
  all-but-one specifically lose). `is_boundary_state()` above doesn't yet
  distinguish these two cases — a natural refinement is checking `board.legal_moves`
  count independent of the non-losing count, and flagging (not necessarily
  excluding) the degenerate case separately.

---

## 7. Relationship to existing chess-computing literature — checked directly, not assumed

This needed the same treatment as the KBNvK findings document's own literature
check: verified against real sources before writing anything down, not asserted
from a general sense that "zugzwang is probably related."

**This sits inside a real, decades-old research tradition, and it already covers
this project's own material.** Althöfer and Walter, *"Weak Zugzwang: Statistics on
some Chess Endgames"*, ICCA Journal 17(2), 1994 — an exhaustive, tablebase-based
statistical study across four endgames, one of which is **KBNK**, the exact
material this project has real, verified data for. Guy Haworth's *"Mutual
Zugzwangs in Chess"* (Computer Olympiad Workshop 6, 2001) goes further: a
comprehensive, exhaustive enumeration of mutual zugzwang across essentially every
2-to-5-man endgame, explicitly listing which materials have **zero** such
positions — **KBNK is on that list, by name, confirmed exhaustively: zero mutual
zugzwangs.**

**The precise, important distinction — this project's own formalization is
related to, but not identical to, any of the three established variants found:**

- **Mutual zugzwang (mZZ)**: fix an arrangement of pieces, toggle which side is to
  move, and ask whether *neither* side wants to be the one on move. Requires
  comparing two distinct positions (same pieces, opposite mover).
- **Weak zugzwang** (Althöfer & Walter's own contribution): the same toggle, but
  compares *distance-to-mate* rather than W/D/L category — a softer measure, still
  a two-position comparison.
- **Traditional zugzwang** (the older concept both papers build on, per Roycroft):
  the same toggle, but asks whether the *value category itself* changes — closest
  in spirit to this document's own concept, but still fundamentally a comparison
  across a toggled pair, not a property of one position alone.
- **Boundary board state (this document, §1)**: no toggle at all. Fix the side to
  move as given, and count how many of *that side's own* legal moves stay
  non-losing out of everything available to them. A single-position property, not
  a comparison against a counterfactual mover-swap.

**Honest verdict, calibrated rather than rounded up or down**: this is not
operating in a vacuum, and it is not simply "zugzwang, already done" either. It's
a genuinely adjacent question to a well-established one, built on the same real
infrastructure (exhaustive tablebase enumeration) that already proved workable at
this exact material's scale. The overlap between which positions get flagged by
either lens is presumably large but unmeasured — a real, open, directly checkable
question once §8's build exists.

**The concrete, grounded action this adds to §8's plan**: Haworth's KBNK-zero-mZZ
result is a precise, already-published, external number to check this project's
own boundary-state scan against, the moment it's run on the real KBNvK data
already in hand. Since a boundary board state doesn't require the mutual property
mZZ does, the expectation is that boundary-state count should come back
*larger* than zero even where mutual-zugzwang count is exactly zero — a real,
falsifiable prediction, not just a plausible-sounding one, and a natural first
sanity check for the tool once built.

## 8. Concrete next build, once this document's understanding is confirmed as correctly captured

1. **Bootstrap from Syzygy** for any material already covered (up to 7 pieces),
   using `is_boundary_state()` above (with the §5 collapse applied) — a pure scan,
   no new tablebase construction, real boundary states out immediately.
2. **Extend to this project's own solved materials** (KQvK, KBNvK, and KRPvK once
   complete) the same way `analyze_longest_line.py` already scans a completed
   database for forced plies — same shape, W/D/L pattern instead of distance/escape.
3. **Compare the resulting KBNvK boundary-state count against Haworth's
   already-published zero-mutual-zugzwang result for KBNK** (§7) as the first real
   sanity check — a real, falsifiable prediction (boundary-state count should
   come back larger than zero, since the property doesn't require mutuality),
   checkable the moment step 2 produces real numbers, not just plausible-sounding.
4. **Cache proven boundary states to a small, dedicated file**, separate from the
   full landscape database — this is sound now, unlike §3c's version, because it's
   caching an *output* of already-completed computation, not attempting to replace
   the computation itself.
5. **For a material neither Syzygy nor this project has solved**: no shortcut
   exists past full classification — stated plainly, not glossed over, per §3a/§3b.

---

## 9. What remains open

- Whether `board.legal_moves` count alone is a sufficient, sound way to separate
  "genuine choice-point boundary state" from "only one legal move existed, and it
  happens to draw" (§6) — proposed, not yet tested against real data.
- Whether Syzygy's 7-piece ceiling itself becomes the limiting factor before this
  project's own engine does, for materials richer than KRPvK.
- Whether a boundary-state catalog, once built for a real material, reveals its
  own further structure (families, clustering, anything analogous to the shape
  graph) — genuinely unexplored, and a natural question once real data exists to
  ask it of.
