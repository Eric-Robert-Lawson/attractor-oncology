# Piece-wise Construction, Trajectory Caching, and the Memory/Storage Tradeoff

This document records an architectural understanding that emerged from
directly testing, rather than assuming, how the two existing resources in
this project — the full attractor landscape and the shape trajectory
graph built on top of it — can be constructed incrementally, verified
locally, and cached across pieces. It is deliberately precise about what
has been tested and confirmed versus what is a sound, motivated idea that
has not yet been built or verified — both categories matter, and blurring
them would undersell the real, tested foundation while overselling the
parts that still need engineering work.

This document assumes familiarity with `general_workflow_reference.md`
(the operational pipeline — sweep, classify, reduce) and
`shape_based_compression_and_scoped_generation.md` (what a shape and an
origin are, and the compression vision built on that definition). It does
not repeat either; it sits above both, describing how the landscape those
documents work with can be built and reused piece by piece rather than
only as one monolithic sweep.

---

## 1. The two resources, restated precisely

**Resource A — the attractor landscape.** The retrograde-solved database:
every position mapped to its distance-to-mate and its tied-optimal
move(s), each with an escape contribution. Built by `general_solver`'s
`discover()` + `classify()` fixed point. Local in the sense that any one
entry only encodes the immediate next step, not a full path — but each
entry's value is exact and complete on its own.

**Resource B — the shape trajectory graph.** Built entirely by reading
Resource A (`trace_shape_graph`, in `general_analyzer.py`): a deduplicated
graph of which positions belong to which (mate, escape-count) shapes, with
convergence (multiple origins sharing a downstream position) and
divergence (a tie's branches serving different shapes) both represented
as first-class graph structure rather than left implicit. B has no
independent existence — every edge in it comes from reading A's own
`tied` field at that position. This was confirmed directly earlier in
this project's work, not assumed: B's construction code has no other
source of truth.

---

## 2. Piece-wise construction of Resource A — validated

**Claim tested**: a position's correct value can be established by
exploring only its own forward-reachable set (every position reachable
from it via legal play, not just optimal play), without needing the full
material's exhaustive landscape.

**Why this has to be true, mechanically, not just empirically**: a
position's value depends only on its own descendants — distance and tied
moves are computed from children's values, never from anything upstream.
The forward-reachable set from any starting position is closed under
taking children by construction (every legal child of a reachable
position is, trivially, also reachable) — so nothing outside that set is
ever needed to correctly value anything inside it.

**Verified directly, exhaustively**: seeded `general_solver` with exactly
one KQvK position (`Q:c4 K:h8 k:f8`, a genuine mid-distance case, not a
trivial near-mate one) and compared every resulting position against a
full, all-seeds exhaustive run. Zero missing positions, zero mismatched
distances, zero mismatched tied moves, across all 345,404 positions
compared.

**This is not a new mechanism — it is the same one the existing pipeline
already depends on, just isolated at its smallest possible unit.**
`run_full_sweep.py`'s own chunking, and `discover()`'s own sealing of
already-known positions, are already piece-wise construction in
production use — batches of seeds processed incrementally, with
already-proven positions never re-expanded. This test confirms the
foundation that mechanism has always relied on, at the n=1 limit, rather
than introducing something new.

**The honest, important caveat this same test surfaced**: on KQvK, that
single seed's forward closure reached the *entire* 372,064-position
discovery set — identical to the full exhaustive seed list's own count.
There was no smaller piece to find from this starting position, very
plausibly because a queen's mobility means White and Black can reach
almost anywhere from almost anywhere under merely legal play. Whether
this differs on material with real structural disconnection (the
bishop-color split is the proven example already in this project) is
untested. Correctness of the mechanism and a practical size reduction
from using it are two separate claims — only the first is confirmed.

---

## 3. Full-comparison exclusion — already built and running; what's new is that it decomposes across independent pieces

**This section originally described exclusion itself as something being tested here for the first time. That framing was wrong, and worth correcting directly rather than leaving stand.** Comparing every candidate's fully-determined value and keeping only the best is not a new mechanism — it is `classify()`'s own core loop, already running in production, confirmed directly by reading the actual code rather than assumed:

```cpp
optional<int> best_v, best_esc; uint64_t best_child = 0;
for (auto& ckey : node.children) {
    auto cit = classified.find(ckey);
    if (cit == classified.end()) continue;
    int v = cit->second.distance + 1;
    ...
    if (!best_v || v < *best_v || (v == *best_v && contribution < *best_esc)) {
        best_v = v; best_esc = contribution; best_child = ckey;
    }
}
```

Every candidate not matching the best value found is excluded from ever being recorded as tied-optimal, by construction. This is exact, and it is how every distance and tied-move set in this project's entire database has always been computed.

**What's actually new, tested here, and worth keeping distinct**: whether this same comparison can be reproduced using candidates built *independently*, as separate pieces, rather than requiring them all to emerge from one shared `discover()` + `classify()` run. Took a real tie (`Q:d5 K:h8 k:f8`, distance 7, tied moves `Qe5`/`Qb7`, plus a findable non-tied legal alternative `Qa1`). Built three fully independent, single-seed forward constructions — one per candidate child — computing each child's own distance in isolation, with no shared state between the three runs. The comparison across these three independently-built results exactly reproduced what `classify()`'s own loop would have produced: `Qe5` and `Qb7` correctly identified as tied-optimal (implied parent distance 7, matching exactly), `Qa1` correctly excluded (implied distance 13). This confirms the comparison step doesn't need to happen inside one unified run — it can be reconstructed from pieces built separately and merged afterward, which is the genuinely relevant finding for piece-wise construction, not the existence of exclusion itself.

**What remains actually unbuilt — a narrower, specific claim, not "exclusion" broadly**: stopping a candidate's evaluation *before* its exact value is known, using only a partial bound (true alpha-beta-style pruning). The `classify()` loop above requires a child to be fully, exactly resolved (`classified.find(ckey)` succeeding) before it's compared at all — there is no partial-information path. This is a real, sound, motivated idea, structurally analogous to a well-established game-tree search technique — but it was not tested here, and testing it properly would require new logic in the C++ search itself (a bounded, "prove worse than threshold and stop" traversal, not the current all-or-nothing `--max-nodes` cap), which was explicitly out of scope for this session's work on the engine file.

---

## 4. Caching across independently-built pieces — validated, with one important limit

**Claim tested**: if two pieces are built independently and their forward
closures overlap, the shared portion should be computable once and
reused, not redone.

**Why it's sound, mechanically**: since a position's value depends only
on its descendants, two independently-explored trees that converge on a
shared position compute the identical value for it either way — there is
no scenario where reuse changes the answer.

**Verified directly, via the existing `--preload-from` mechanism**: after
building the first candidate's forward closure (§3), fed it as
`--preload-from` into the second and third candidates' own runs. Both
completed with zero new discovery — every position already sealed from
the first run. The reuse mechanism itself is real, already built, and
confirmed correct.

**The important limit, checked directly in the source rather than
assumed, and explicitly not yet resolved**: `--preload-from` does not
currently stay memory-bounded. `load_general_db_for_resume` loads the
entire preload source into an in-memory map before integrating it — so
the caching just verified achieves correctness by holding the full prior
result in RAM, not by reading only what a new piece's boundary actually
needs from disk. This is the same shape of constraint that caused this
project's own KBPvK memory crisis (db + origins + parent_counts +
resolved, all simultaneously resident) — reusing prior work currently
means re-loading all of it, not selectively reading the relevant slice of
it.

---

## 5. What this means for scale, stated honestly

The value of this architecture is not, on the evidence gathered so far,
"less total computation." On a densely-connected material like KQvK, the
first piece built from any starting point costs the whole landscape
regardless of where you start — there was no smaller sub-problem to find.
What is real and confirmed: correctness is preserved at every granularity
tested, down to n=1, and reuse of already-proven work across pieces is
sound and already mechanically available.

What separates the current, real capability from the larger vision is
specifically an I/O-engineering problem, not a correctness or research
question: making cross-piece reuse read only the boundary a new piece
actually touches, rather than the entirety of whatever it's reusing.
Solving that would turn "hold the whole landscape in RAM at once" into
"hold one piece plus a thin, on-demand boundary layer at a time" — the
actual mechanism that would let piece-wise construction reduce *peak*
memory rather than only redistribute *total* work over time. This has not
been built or tested; it's recorded here as the concrete, well-scoped
next engineering step, not as something already delivered.

Whether piece-wise construction also reduces *total* work (not just peak
memory) on material with genuine structural disconnection remains open
too (§2) — a material where different starting regions provably don't
need to know about each other (unlike KQvK, where they apparently all
converge) is the natural next test.

---

## 6. Storage-vs-access tradeoff, and a real indexed lookup — tested, not just discussed

**The reframe that actually matters, and it's a standard, well-established engineering pattern, not a new insight specific to this project**: total storage size being larger is not itself a problem if it buys genuinely better access. Search engines' inverted indices are larger than the text they index; database indices add size specifically to buy fast lookup. §5's framing of "bigger storage" as a cost needed this correction.

**Where this reframe actually lands, once Syzygy is accounted for accurately**: Syzygy's own format is already access-efficient — compressed, block-addressable, decompressed on demand (§ discussed earlier in this conversation, confirmed via direct source). So "better access than a flat table" doesn't by itself differentiate this project's approach from what already exists. The real edge stays exactly where it was: the escape-count dimension, which Syzygy's WDL/DTZ format has no equivalent of at all.

**Real, confirmed precedent for federated/distributed hosting of exactly this kind of data**: Syzygy tables today are already served through multiple independent, cooperating entities simultaneously — HTTP mirrors (`tablebase.sesse.net`), a BitTorrent tracker, and separate API servers (`lichess.org`'s public tablebase API, `syzygy-tables.info`) all serving the same underlying dataset, none of them exclusive. "Many centralized entities cooperating, not one replacing the other" is not a design risk to introduce here — it's the proven, existing pattern. Separately, generating the first 7-piece Syzygy set reportedly needed over 1TB of RAM on specialized hardware, later distributed across multiple machines — real, historical confirmation that construction-time memory was a genuine bottleneck for exactly this class of problem, at real scale, for the system this project is being compared against.

**A real indexed lookup tool was built and measured against real KQvK data, not sketched**: `build_shape_index.py` takes `--trace-shape-graph`'s own output and builds a normalized SQLite database — `nodes`, `node_shapes` (indexed on both `(position, turn)` and `(shape_id, escape_count)`, so both query directions are indexed, not just one), and `edges`. Measured directly, same 345,404-position KQvK graph both ways:

| Query | Linear CSV scan | Indexed SQLite | Speedup | Correctness |
|---|---|---|---|---|
| Position → shape membership (20 random positions) | 168.0 ms/query | 0.04 ms/query | **4,696x** | Exact match, all 20 |
| Shape → full membership (one real shape, 4 positions) | 611 ms | 0.1 ms | **4,404x** | Exact match |

**The honest storage cost, measured on the same data, illustrating the tradeoff directly rather than asserting it**: the indexed SQLite file is 285.8 MB, against 46.5 MB combined for the two source CSVs it was built from — roughly 6.1x larger. This is precisely the tradeoff §6's opening paragraph describes, shown with real numbers on real data rather than argued abstractly: substantially more storage, in exchange for a genuinely massive, measured, correctness-verified access improvement.

**Still open, stated plainly**: this is tested only at KQvK's 345,404-position scale. Whether the index-building time, the index-to-source size ratio, and the query speedup all hold up at 6-7 piece scale — where the underlying graph itself could be orders of magnitude larger, per §2's and §5's own caveats about what's actually been measured — is untested.

---



## 7. Relationship to the multi-piece-Black refactor — a distant, contingent next step

`loser_count_and_multipiece_black_refactor.md` describes a separate,
larger change: letting Black hold real material and potentially win,
moving from a binary (White-wins-or-draws) outcome model to a genuine
three-way one. That document already notes, from earlier work, that the
shapes/origin machinery is mostly turn-agnostic and would survive that
transition with only its terminal-detection step needing revision (see
its own §Part 4, item 6) — a solved `(position, turn)` already implies
who won, so no new field is needed in a shape's own identity.

This document's own findings compound on top of that rather than
requiring rework: piece-wise construction's correctness argument (§2) —
a position's value depends only on its own descendants — does not depend
on which side is which or how many pieces Black holds; it is a property
of retrograde valuation itself. If and when the ternary-outcome engine
exists, the same single-seed and local-verification tests run here should
be re-run against it directly, not assumed to carry over — but there is
no structural reason to expect them not to. This is recorded explicitly
as a *distant* dependency: the multi-piece-Black work has its own
substantial, separately-scoped open items (per its own Part 4) and isn't
blocked on anything in this document, nor does it block anything here.

---

## 8. Summary table, for fast re-entry

| Claim | Status | Evidence |
|---|---|---|
| Single-seed forward construction is exactly correct | **Validated** | 0/345,404 mismatches vs. full exhaustive run |
| Single-seed construction produces a *smaller* sub-landscape | **Open — negative result on KQvK** | Same seed reached all 372,064 discovered positions |
| Full-comparison exclusion (keep best, discard rest) is correct | **Already built and running** — this is `classify()`'s own existing loop, not newly tested | Read directly from the engine source |
| This same comparison decomposes correctly across independently-built pieces | **Validated** | Exact match to ground truth on a real tie + non-tied alternative, built as 3 separate single-seed runs |
| Bound-based early exclusion (stop before a candidate's exact value is known) is cheap | **Untested, unbuilt** | Would require new C++ search logic, out of scope this session |
| Convergent pieces can be safely merged/reused | **Validated, mechanism sound** | `--preload-from` gave zero new discovery on reused pieces |
| Reuse is memory-bounded (reads only what's needed, not the whole prior piece) | **Not validated — confirmed false as currently built** | `load_general_db_for_resume` loads the full source into RAM |
| Peak-memory reduction at real scale | **Open, well-scoped next step** | Depends on solving the row above |
| Indexed lookup gives real, large access-speed gains over flat CSVs | **Validated on KQvK** | 4,696x (position→shape) and 4,404x (shape→membership), both correctness-verified; index itself is ~6.1x larger than its source CSVs |
| Indexed lookup approach holds up at 6-7 piece scale | **Open, untested** | Only exercised at 345,404 positions so far |
| Total-work reduction on structurally disconnected material | **Open** | Untested; bishop-color split is the natural next test case |
