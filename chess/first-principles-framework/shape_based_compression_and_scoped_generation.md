# Shape-Based Compression and Scoped Landscape Generation

A design document for a methodology derived from this project's own two-stage
perfect-play criterion (minimize distance-to-mate, then minimize/maximize
cumulative defender escape count among distance-ties) — developed
independently, from first principles, not adapted from existing tablebase
literature. This document separates what has been **built and verified**
from what is **proposed and open**, and keeps those two categories visually
distinct throughout, because conflating them is exactly the failure mode
this project has hit and corrected multiple times already.

---

## 1. The foundational distinction: shapes vs. origins

**Shape** — the identity `(mate, cumulative escape count)`. Two positions
share a shape if and only if perfect play from each eventually reaches the
identical mate (up to board symmetry) having accumulated the identical total
escape count to get there. This is *not* just "which mate" — the same mate
is routinely reachable via genuinely different escape totals from genuinely
different starting points (confirmed directly: 88,410 real cases on KQvK
alone where the identical mate, from the identical position, was reachable
via two different totals — e.g. 21 vs 17 escapes).

**Origin** — a position where a shape genuinely *begins*: either a tie
(where perfect play branches to two or more distinct possibilities, some of
which may resolve to different shapes — this is where shapes *intersect*,
not where a new one is created by the tie itself), or a position with no
forced-in parent at all (nothing forces play into it — either a true
starting point, or a position only reachable via a sub-optimal deviation
from perfect play elsewhere; both are the same underlying phenomenon, since
neither is the forced target of any optimal move).

Every position that is *not* an origin is simply further along an
already-established shape's own forced flow — not a new decision point,
and not a new catalog entry. This distinction was arrived at after two
earlier, incorrect versions of the mechanism, both confirmed wrong with
concrete evidence before being replaced (see §6 for the history) — it is
not a refinement of an earlier idea, it is the corrected foundation.

**Verified on KQvK** (345,404 positions):
- 264,532 origins (76.6% of the entire landscape)
- 46 distinct mates
- 6,436 distinct (mate, escape) shapes

---

## 2. What "shape count" does and does not tell you

This is the single most important distinction for everything that follows,
and it was not obvious at first — it took a wrong turn to surface clearly.

The **shape count** (6,436) measures *outcome diversity* — how many
genuinely different final (mate, escape) destinations the whole landscape
resolves into. It does **not** measure how many distinct *move sequences*
exist. Two origins can share the identical shape (same eventual mate, same
escape total) while having completely unrelated specific paths to get
there — different squares, different pieces moving, no relationship beyond
landing on the same labeled destination.

The number that actually matters for **storage compression via shared
sequences** is the **origin count** (264,532), not the shape count. Only
1.3x smaller than the total position count on KQvK — most of the dramatic
reduction from 345,404 down to 6,436 comes from many *different* origins
converging on the *same label*, not from many positions sharing the *same
sequence*. This is a real, load-bearing distinction: it means the shape
catalog by itself is not a compression scheme, even though it looks like
one at first glance.

---

## 3. The compression mechanism that *is* real and available now

**Verified, built, real:** any non-origin position is fully reconstructible
from exactly two integers — `(origin_id, ply_offset)` — plus that origin's
own stored move sequence. Its own distance, escape count, and best move
never need to be stored independently; replaying the origin's sequence
forward from the given offset reconstructs them exactly. This scales with
*chain length*, not chain count: a position 40 plies into a long forced run
costs the same two integers as a position 2 plies in.

On KQvK, this directly compresses the 23.4% of the landscape that isn't an
origin. The origins themselves — the 76.6% — are the genuine hard core this
mechanism does not yet reduce.

**Proposed, not yet built:** whether the hard core itself compresses
further depends on whether origins sharing a shape *also* share a sequence
via symmetry — i.e., extending the same consistent-transform check
`analyze_reducibility` already applies within one finding's own tied
branches (Phase 1 of `reduce`) across the *entire* landscape's origins, not
just within single ties. If origin A's whole sequence is a D4-transform of
origin B's whole sequence, that pair is genuinely compressible: store one
sequence plus a transform index for the other. This has real precedent —
Syzygy tables already exploit D4 symmetry to avoid storing 8 copies of
mirror-image positions — but has not been measured at the *origin* level
here, and the actual compression ratio this would deliver on KBPvK is
currently unknown, not merely unstated.

---

## 4. Scoped generation — four genuinely different things worth naming separately

The phrase "generate only what you need" covers four distinct mechanisms
with very different levels of readiness. Conflating them would overstate
what's currently possible.

### 4a. Scoped forward construction (already available, not shape-dependent)
If only specific starting position(s) matter, `discover()` can already be
seeded with just those position(s) instead of the full exhaustive seed
list — this naturally scopes discovery to only what's reachable forward
from the positions of interest. This requires no shape catalog at all; it's
an existing capability of the current pipeline (`generate_general_seed_positions.py`
for a small, targeted seed set, vs. `generate_exhaustive_positions.py` for
full combinatorial coverage).

### 4b. Scoped querying of an already-solved landscape (proposed, buildable, not built — but ordinary engineering, not an open research question)
Once a landscape has been *fully* solved once — on hardware that can afford
it — and its origins' sequences have been extracted and stored (§3), a
resource-constrained machine could answer "what's the perfect-play
continuation from this specific position" by looking up which origin it
belongs to and replaying only that origin's stored sequence, without ever
loading the full per-position database. This is analogous to how real
tablebase implementations load and decompress only the relevant on-disk
chunk for a query rather than the whole table.

Worth being precise about what's trivial here versus what actually needs
building, since these were conflated in an earlier version of this
document. Determining *which material* an arbitrary position belongs to is
not a problem at all — it's read directly off the piece list on the board,
no computation required. What does need building is a `position →
(origin_id, ply_offset)` index for that material, constructed once, directly
from its own already-computed classify() output — `resolve_reachable_shapes`
already produces almost exactly this data during its own pass; recording it
keyed for direct lookup rather than iterated during aggregation is a
straightforward extension, not a new capability. One real wrinkle: a
non-origin position can in principle have more than one forced-in parent,
so it may belong to more than one origin's chain simultaneously — the index
needs to account for a set of (origin, offset) pairs per position, not
assume a clean one-to-one mapping.

The actual open question, once this index exists, is empirical, not
"how": is a `(origin_id, offset)` pair per position, plus each origin's
sequence stored once, meaningfully smaller in total than the raw
per-position table it replaces? That's measurable directly and hasn't been
measured yet — it isn't an unsolved design problem.

### 4c. Reusing already-solved, reduced materials during a larger material's construction (built and verified, not just grounded in theory)

This is different from 4a and 4b in kind, not just degree — it is about a
*larger* material's own construction reusing *already-independently-solved,
smaller* materials that its own tree provably reduces into via capture or
promotion. It does **not** involve predicting anything about an unrelated
material, and it is **not** circular — the smaller material is solved
completely on its own, with no dependency on the larger one, before the
larger one ever needs it.

**This is not an open question — the dependency is already known.** This
project already identified, when discussing Syzygy validation, exactly
which reduced materials KBPvK's own exhaustive tree touches: KBQvK and
KBNvK (and KBRvK under `--full-promotion`) when the pawn promotes with the
bishop still on the board; KQvK, KNvK, and KRvK when the bishop is captured
either before or after promotion; KBvK and KPvK when one attacking piece is
captured before the other does anything. Every one of these is a strictly
smaller, independently solvable material — solving KBPvK never requires
knowing anything about KBPvK to solve any of them first.

**No longer just a theoretical mechanism — a real, working, tested CLI
flag.** `general_solver` now takes `--preload-from <file>` (repeatable),
seeding a run with an already-completed, separately-solved material's
`.db` file, in addition to (not instead of) the normal same-material
resume. Built directly on `pack_general_state`'s confirmed material-
agnostic packing and `seed_from_preloaded`'s existing soundness — a proven
fact's truth never depends on which run discovered it.

Verified with a real, concrete, dramatic result, not just a design
argument: fed a single already-solved KQvK position back in as a root for
a fresh run. Without `--preload-from`, the run rediscovers the entire
372,064-position graph and runs all 21 classification passes again. With
`--preload-from` pointing at that same material's own completed database,
the root is sealed immediately, discovery finds **zero** new positions,
and classification completes in a single trivial pass, reporting exactly
how many facts came pre-proven and how many re-verification passes were
spent on them (zero). Same starting position, same engine, only the flag
differs.

Also built and verified: a hard safety check. `seed_from_preloaded`
unconditionally overwrites on a matching key, which is fine if multiple
preload sources agree and dangerous if they silently don't — so before
trusting any cross-material file, its facts are checked against whatever
this run already knows, and a genuine disagreement between two sources
aborts the run with the exact mismatched values shown, rather than
silently picking one. Tested both ways: a single, uncontested (deliberately
corrupted) source was accepted, since there is nothing to disagree with — a
real, honest limit of what this check can catch, not a bug — while feeding
that same corrupted file *alongside* the correct one triggered the fatal
error exactly as intended.

**Not yet built: this flag isn't wired into `run_full_sweep.py`.** Checked
directly — the sweep wrapper constructs its solver command explicitly
(`--full-dag`, `--positions`, `--db`, optionally `--max-nodes`/`--fresh`)
with no passthrough for arbitrary extra flags, so `--preload-from` is
currently only usable by invoking `general_solver` directly (§2 of the
workflow reference), not through a full, multi-batch sweep. Adding
explicit `--preload-from` support to the sweep script is a small, separate
piece of work, not done here.

**Where shapes specifically would add value on top of what's built now,
and where the real open work remains:** feeding the *entire* exported
database of a reduced material works today, with nothing further needed,
but the reduced material's own full database can itself be large — that's
the real memory cost `--preload-from` pays today, holding each preloaded
material's complete table in memory during the larger material's own run.
The open, genuinely useful extension is using that reduced material's own
*shape/origin catalog* — not its raw per-position table — to import only
the origins actually reachable given the larger material's specific
context (e.g., only the KBQvK shapes whose governing bishop square is
consistent with wherever KBPvK's own bishop can actually be when a given
promotion occurs), reconstructing exactly those via §3's origin-sequence
mechanism rather than loading the reduced material's full table into
memory at all. That filtering-and-reconstruction step is not built yet;
feeding the full reduced-material database directly, which now works, is
the version available today.

### 4d. The recursive, whole-hierarchy version of 4c — the actual scale of what's being proposed

4c described reusing *one* already-solved reduced material during *one*
larger material's construction. The fuller claim is recursive, and worth
stating at its actual scope rather than one level at a time: **no material
anywhere in a dependency chain needs its full, raw per-position table
stored, ever — only each material's own *complementary* set, the shapes and
origins that are genuinely unique to it and not already covered by
something simpler it reduces into.**

Concretely for KBPvK: its own landscape splits into (a) positions where the
pawn is still a pawn and the bishop is still present — genuinely new
territory no simpler material has — and (b) positions that are *already,
literally* KQvK, KBNvK, KBQvK, KNvK, KRvK, KBvK, or KPvK positions, reached
via a capture or promotion. Category (b) is not KBPvK data at all — it's a
reference into a simpler material's own, separately-stored catalog. The
genuinely-new-to-KBPvK complementary set is necessarily smaller than the
origin fraction measured on KQvK alone (§1), because a real portion of what
a naive, undifferentiated computation would count as "KBPvK's own origins"
actually belong to one of those simpler materials' catalogs instead.

Applied recursively down the entire dependency graph — KBPvK's own
complement, plus KBQvK's own complement (everything in KBQvK not already
covered by KQvK), plus KQvK's own complement (its base case), and so on for
every branch — the **total** storage across the whole hierarchy is the sum
of each material's own non-reducible remainder, not the sum of everyone's
complete table. This is a materially different shape of solution than
storing each material's table whole and cross-referencing only at probe
time, which is closer to how Syzygy appears to operate, to the best
understanding here — worth restating as a real, acknowledged uncertainty
about Syzygy's exact internals rather than a confirmed comparison.

**Why this should matter more, not less, at larger scale, and this is a
real, checkable prediction rather than a hope:** as material complexity
grows toward 7-piece territory, there are more distinct ways to capture or
promote into some simpler, already-solved material, and the "still
combining every piece at once" territory should shrink as a *proportion* of
each material's own exponentially-growing position count. If that holds,
the complementary-set approach doesn't just save a fixed amount of
storage — its relative advantage over storing full tables at every level
should *grow* with scale, which is exactly the regime where Syzygy's own
storage costs already become the binding constraint. This is stated as a
prediction to verify, not an established result — nothing here has measured
whether the complementary fraction actually shrinks this way in practice.

**The piece this entire recursive scheme depends on, restated precisely
after an earlier version of this document overstated it:** every step
above requires a `position → (origin, offset)` index for each simpler
material, built once from that material's own already-computed classify()
output (§4b). This is ordinary engineering, not an unsolved problem —
identifying which material a position belongs to is trivial from its piece
list, and building the index itself is a straightforward extension of data
this pipeline already produces during `resolve_reachable_shapes`. What
actually remains open is measurement, not method: is the resulting index,
plus each origin's sequence stored once, meaningfully smaller in total than
the raw per-position table it replaces — for KQvK, for KBQvK, and at
whatever scale this gets pushed toward. That's the number nothing here has
measured yet, and it's the one the whole recursive storage argument in this
section actually rests on.

---



Stated honestly, since the distinction matters for evaluating this work
correctly:

- Compressing a solved endgame table by exploiting structural redundancy
  rather than storing every position independently: **has precedent**
  (Syzygy's own D4 symmetry handling; other structural-compression work in
  computer chess).
- Using **escape count as a second, independent optimality dimension
  alongside distance-to-mate**, and deriving shape/origin identity from
  both together: **no known precedent**. Standard tablebases use DTZ/DTM
  alone. This is the lens that surfaced real, non-obvious structure a
  distance-only analysis would never have found (the 88,410 escape-count
  violations; the shape/origin distinction itself).
- Reusing an already-solved, smaller material's results when constructing a
  larger material that reduces into it (§4c): **has strong precedent** —
  this is exactly how tablebase construction is normally done, building up
  from smaller materials since captures only ever reduce material. Not a
  novel idea on its own, though now actually built and verified here (the
  `--preload-from` flag), not just argued to be sound by analogy.
- Using that smaller material's *shape/origin catalog* specifically, rather
  than its full per-position table, to import only the reachable-in-context
  subset during that reuse (§4c's open extension): **appears genuinely
  open** — not confirmed as prior art, and not yet built here either.
- The recursive, whole-hierarchy version of this — no material anywhere in
  a dependency chain needing its full table stored, only each one's own
  non-reducible complement (§4d): **appears genuinely open, and is a
  materially different shape of solution than storing each material's
  table whole**, to the best understanding here of how Syzygy's own storage
  operates — stated with real, acknowledged uncertainty about Syzygy's
  exact internals, not as a confirmed comparison. The index each level
  depends on (§4b) is ordinary engineering, not unsolved — what's actually
  untested is whether it nets out smaller than the raw tables it would
  replace, at any real scale.

---

## 6. Corrections along the way, kept as a record rather than erased

This methodology went through real, confirmed wrong turns before reaching
the current definition, worth keeping visible rather than smoothing over,
since each correction is what actually established the current definition
is right:

1. **Canonicalize each position independently** (position + its own
   required move, D4-deduped). Rejected: only recovers ordinary 2-to-8-fold
   board symmetry, barely reducing anything (KQvK: 345,040 → 43,195).
2. **Walk forward to the nearest terminal, treating ties and mates as
   equivalent stopping points.** Wrong on two counts, both confirmed: a tie
   is an intersection point between shapes, not a terminal itself; and
   canonicalizing by position alone drops escape count entirely (the
   88,410-violation finding).
3. **Escape-count-aware, but counting every position's own resolved set.**
   Still wrong: a single, real, zero-branching 14-ply chain on KQvK reported
   a different escape value at every ply (36, then 6, then 0) and was being
   counted as 14 separate shapes, when it is unambiguously one. Fixed by
   introducing genuine origins (§1) — the current, verified definition.

---

## 7. Open questions, stated as questions

- **The central one (§4d): once a `position → (origin, offset)` index is
  built for a material — an ordinary engineering task, not an open one —
  is it, plus each origin's sequence stored once, actually meaningfully
  smaller than the raw per-position table it replaces?** This is
  measurable directly and hasn't been measured yet for any material here;
  every claim in §4b/4c/4d about real storage savings depends on this
  coming back favorable, not on solving some unsolved mechanism.
- Does the fraction of a material's landscape that's genuinely
  non-reducible (§4d) actually shrink as a proportion of total size at
  larger piece counts, as predicted, or does it stay roughly constant or
  even grow? This determines whether the recursive storage scheme's
  advantage genuinely compounds toward 7-piece scale or stays fixed.
- What fraction of KBPvK's landscape is origins vs. non-origin? (KQvK's
  76.6% may not generalize — richer piece mobility could push this higher
  or lower; no confident prediction exists yet.)
- Of KBPvK's origins that share a shape, how many also share a sequence via
  symmetry? This is the number that determines whether §3's proposed
  extension delivers meaningful compression on the hard core, or only on
  the non-origin fraction.
- For §4c's open extension specifically: given a reduced material's shape
  catalog, what's the actual mechanism for determining which of its origins
  are reachable-in-context from the larger material (e.g., which KBQvK
  bishop-square origins are even possible given KBPvK's own bishop's
  reachable squares at the moment of promotion), and does filtering to that
  subset before reconstruction meaningfully reduce memory versus `--preload-from`'s
  now-built and verified full-database approach?
- Does the escape-count-as-second-dimension lens generalize usefully to
  materials with more than one non-king attacking piece in ways not yet
  explored here, or does its value concentrate specifically in materials
  with rich promotion/deviation structure like KBPvK?
