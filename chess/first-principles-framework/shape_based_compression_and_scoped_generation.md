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

## 4. Scoped generation — three genuinely different things worth naming separately

The phrase "generate only what you need" covers three distinct mechanisms
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

### 4b. Scoped querying of an already-solved landscape (proposed, buildable, not built)
Once a landscape has been *fully* solved once — on hardware that can afford
it — and its origins' sequences have been extracted and stored (§3), a
resource-constrained machine could answer "what's the perfect-play
continuation from this specific position" by looking up which origin it
belongs to and replaying only that origin's stored sequence, without ever
loading the full per-position database. This is analogous to how real
tablebase implementations load and decompress only the relevant on-disk
chunk for a query rather than the whole table. This requires building: the
origin-sequence storage format, and a lookup mechanism mapping an arbitrary
query position to its origin. Neither exists yet.

### 4c. Using shape/origin knowledge to prune a *new* landscape's construction (proposed, genuinely unresolved)
The most ambitious version: using an already-known catalog to decide, during
a *fresh* construction, which branches don't need to be discovered or
classified at all because they're already known not to matter. This is
meaningfully different from 4a and 4b, and has a real, unresolved tension at
its core worth stating plainly rather than glossing over: `classify()`'s own
fixed-point process is what *determines* a position's true value and hence
which shape it belongs to — knowing "which shape a position will turn out
to belong to" *before* running that process is close to circular. Making
this work would require some independent way to predict shape-relevance
ahead of full classification, which is not solved by anything built so far.
This is the part of the idea that is genuinely novel and genuinely
unresolved, not merely unimplemented.

---

## 5. What's genuinely novel here, and what has real precedent

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
- Using a solved landscape's own derived structure to scope or prune the
  *construction* of a related landscape (§4c): **appears genuinely
  open** — not confirmed as prior art, but also not yet demonstrated
  working here.

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

- What fraction of KBPvK's landscape is origins vs. non-origin? (KQvK's
  76.6% may not generalize — richer piece mobility could push this higher
  or lower; no confident prediction exists yet.)
- Of KBPvK's origins that share a shape, how many also share a sequence via
  symmetry? This is the number that determines whether §3's proposed
  extension delivers meaningful compression on the hard core, or only on
  the non-origin fraction.
- Is there a principled way to predict shape-relevance ahead of full
  classification (§4c), or does scoped pruning of a *new* construction
  remain fundamentally dependent on having solved the related landscape
  first?
- Does the escape-count-as-second-dimension lens generalize usefully to
  materials with more than one non-king attacking piece in ways not yet
  explored here, or does its value concentrate specifically in materials
  with rich promotion/deviation structure like KBPvK?
