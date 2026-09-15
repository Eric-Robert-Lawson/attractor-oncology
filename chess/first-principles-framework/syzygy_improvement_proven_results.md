# Improving on Syzygy: What's Proven, What's Fixed, What's Still Open

A record of the specific, verifiable claims this project has established toward
its goal of improving on Syzygy — kept separate from the broader design
documents (`shape_based_compression_and_scoped_generation.md`,
`piecewise_construction_and_trajectory_caching.md`) so there's one place to
come back to that states plainly: this was tested, here's the real number,
here's what it does and doesn't show yet. Written after finding and fixing a
real, serious bug in the mechanism Claim B depends on — the bug and its fix
are part of the record, not smoothed over, because the fact that it was
caught and verified is itself part of why the current result can be trusted.

**Revision note:** this version adds Claim C, which the original version of
this document didn't separate out on its own — the trajectory graph as a
materialized structure in its own right, distinct from both the escape-count
dimension it's built from (Claim A) and the cross-material storage mechanism
(Claim B). It also adds an explicit section on distribution/hosting,
scoped honestly as precedent rather than a result this project has produced.
Nothing in Claims A or B below has changed; §§6–7 are new.

---

## 1. The three claims this project has made about improving on Syzygy

Stated separately because they're genuinely different kinds of claims, with
different evidence behind them:

**Claim A — escape count as a second, independent optimality dimension
alongside distance-to-mate.** Standard tablebases (Syzygy included) use
DTZ/DTM alone. This project's two-stage criterion — minimize distance first,
then minimize/maximize cumulative defender escape count among distance-ties —
has no known precedent. Established earlier in this project, before this
document: 88,410 real cases on KQvK alone where the identical mate, from the
identical position, was reachable via two different escape totals. Every
shape's identity in this project's own output is `(mate, escape count)`, not
`mate` alone — this claim is load-bearing throughout everything else here,
including Claims B and C.

**Claim B — recursive, complement-only storage across a material dependency
hierarchy.** No material needs its own full table stored; only its
non-reducible complement — the positions that aren't already, literally, some
simpler material's own position, reached via a capture or promotion. §§2–5
below are the real, verified test of this claim's mechanism, including a bug
that was found and fixed along the way.

**Claim C — the shape/trajectory graph as a materialized structure that
Syzygy has no equivalent of, and that avoids repeated re-derivation of
family-level relationships.** Claim A establishes that ties at a given
distance aren't all equivalent — some reach the mate having forced more
escapes than others. That fact, by itself, is a per-position property. Claim
C is about what happens once you go one level up: grouping tied positions
into families by `(mate, escape count)` and mapping the transitions between
those families — which shape leads to which other shape, and how — is
*relational* information that doesn't exist anywhere in a flat distance/
escape/best-move table. Deriving it requires walking the landscape's own
structure, and doing that walk is a real, measured cost, not a hypothetical
one (§6). Materializing it once, so that family- and trajectory-level
questions become a lookup instead of a fresh derivation every time, is what
this project calls Resource B (the shape trajectory graph), built on top of
Resource A (the full attractor landscape Claim A's values live in).

---

## 2. Claim B: the mechanism, and the real test built to check it

`--preload-from <path>` lets a material's own solve consult an
already-solved, simpler material's database. A position reached via a capture
or promotion that's already known from a preloaded source should be excluded
from the *current* material's own output entirely — not recomputed, not
duplicated — living only in the simpler material's own file.

**The test:** built KQvK and KNvK (both real, independently-solved
materials — KQvK winning, 345,404 positions; KNvK always-draw, insufficient
material), then solved KPvK twice from the same exhaustive seed set:

- **Baseline** — no preloading. KPvK's own output necessarily contains its
  own copy of every KQvK position reachable via promotion (a pawn promoting
  to a queen produces a literal KQvK position), redundantly, in addition to
  KQvK's own separate file.
- **Preloaded** — `--preload-from` pointing at both KQvK's and KNvK's
  databases simultaneously. This is a genuine two-dependency test, not a
  single-material toy case — matching the actual shape of the claim (a
  material reducing into *several* simpler ones, as the project's own
  KBPvK example already described).

This specific pair (KPvK depending on both KQvK and KNvK, via promotion) was
chosen because it's small enough to fully complete — not aborted, not
scoped down — while still being a real, non-trivial, multi-dependency case.

---

## 3. Claim B: what's confirmed — the exclusion mechanism itself

**Exact match, not approximate.** Baseline KPvK's own output contains exactly
345,404 rows that are already-promoted (non-pawn) positions — an exact match
to KQvK's own full classified count. The preloaded version contains **zero**
such rows. Every single one is correctly excluded, confirmed by direct
row-count search (`grep -c "^[QN]:"`), not inferred.

**Real, measured storage reduction.** Baseline output: 22,821,324 bytes.
Preloaded output: 10,159,596 bytes (after the fix in §4 — the number before
the fix was similar in size, but built on wrong values, see below). A genuine
55.5% reduction, on real data, for this specific pair.

**This is a multi-dependency result, not a single-level one.** Both KQvK-
and KNvK-derived positions were simultaneously, correctly excluded in the
same run, directly demonstrating the *recursive* framing the claim depends
on — a material reducing into more than one simpler dependency at once, not
just a single parent-child pair.

---

## 4. Claim B: the bug found while verifying further — and why it mattered

Verifying the exclusion mechanism wasn't the end of the check. A further,
harder question: do the *genuinely new* KPvK positions — the ones that still
have a pawn, never excluded at all — get the *same values* whether or not
preloading is used? They didn't.

**The bug, precisely.** `classify()` is a fixed-point induction. For
White-to-move positions, it locks in a value the moment *any* child is
known, picking the best among currently-known children — safe only because,
without preloading, children become known strictly in increasing-distance
order (a natively-discovered position of true distance D only ever enters
`classified` during the pass matching D). Preloading broke that guarantee: a
preloaded position's distance was injected into `classified` immediately,
regardless of its actual size — "pass 0," not the pass where it would have
naturally arrived. A White-to-move position with one preloaded child
(available instantly, distance 12) and one native child that would
eventually resolve better (distance 8, several passes later) locked in on
the worse, preloaded option first. Once locked in, `classify()`'s own
`if (classified.count(key)) continue;` never reconsidered it.

**Confirmed directly, traced to one specific position**, not inferred from
aggregate statistics: `P:a7 K:d3 k:g3` (White to move). Baseline correctly
found an 11-ply mate via `Ke4`. The unfixed preloaded run missed this
entirely, settling for a 13-ply mate via immediate promotion — because
`Ke4`'s own subtree (verified, separately, to be identical and correct in
both runs) hadn't been given the chance to resolve before the preloaded
alternative became available.

**The fix.** Preloaded values are now staged by their own distance and
released into `classified` only at the pass matching that distance —
restoring the exact timing a natively-discovered position of that distance
would have had. `discover()`'s own sealing and enqueue checks needed a
matching update, since they read `classified` directly and would otherwise
have incorrectly re-expanded staged-but-not-yet-released positions.

**Verified after the fix, not just argued to be correct:**
- Rebuilt the no-preload baseline with the fixed binary — **byte-for-byte
  identical** to the pre-fix baseline, confirming the fix changes nothing
  when preloading isn't used.
- Re-ran the preloaded KPvK test with the fix: all 222,540 genuinely-new
  positions are now **byte-for-byte identical** to baseline. Zero
  differences, checked via full diff, not sampling.
- Re-confirmed the exclusion mechanism itself still holds after the fix:
  zero promoted rows leak into the output, same as before.

This bug is now fixed in `compositional_trajectory_solver_modular_shapes.cpp`.
Any database built with a version of this engine predating the fix, using
`--preload-from`, should be treated as suspect for its *values* specifically —
the exclusion/storage behavior was always correct; the *distance and
escape-count values* for genuinely-new positions were not, whenever a faster
native alternative existed that hadn't yet resolved by the time a preloaded
alternative became available.

---

## 5. Claim B: what this proves, stated precisely

The mechanism Claim B depends on is real, does what it's supposed to, and —
as of the fix in §4 — produces correct values while doing it. This was
checked, not assumed: exact exclusion counts, byte-for-byte value matches,
and a real, measured file-size reduction, all on genuine multi-dependency
data. The recursive storage-compression idea is no longer just a mechanism
that *should* work by construction — it's one that's been watched working
correctly on real data, including catching and fixing the one real way it
was silently wrong.

---

## 6. Claim C: the trajectory graph as materialized derived data

**What Resource A and Resource B actually are, precisely, since Claim C
depends on keeping them distinct.** Resource A is the full attractor
landscape — every position mapped to its distance, best move, and escape
contribution. It's a *pointwise* structure: ask it about one position, get
one answer. Resource B is the shape trajectory graph — positions grouped
into shapes by `(mate, escape count)`, with edges recording which shapes'
moves lead into which other shapes. It's a *relational* structure built on
top of Resource A, not a reformatting of it: no single row of Resource A
tells you which families a position's shape connects to, because that's a
property of the graph, not of any one position.

**The claim, stated as an argument, not yet as a controlled measurement
(see the honest gap below).** Multi-way ties at a given distance are common
— Claim A's own 88,410-case figure on KQvK alone is evidence of that. Once
escape count is taken seriously as a second dimension, a natural further
question is family-level, not position-level: which tied lines, as a group,
lead to which other groups, and how do those groups relate across a game's
course. Answering that from Resource A alone means walking the landscape's
own structure to reconstruct it — for every such question, unless the
relational structure is computed once and kept. Resource B is that
computed-once structure.

**Evidence that this derivation is a real, non-trivial cost — not
hypothetical, not assumed.** The clearest direct evidence is operational,
not rhetorical: checkpointing had to be added to both the cluster pre-pass
and the segmented write pass specifically because an uncheckpointed run, at
real scale (KBNvK, ~7.4 million origins), lost real, non-recoverable
progress with nothing to resume from. That's not a design choice made out
of caution — it's a fix made *after* the cost of not having it was already
paid once. If reconstructing this relational structure were cheap, losing
an in-progress run wouldn't have mattered enough to build and verify
checkpoint/resume (including genuine `SIGKILL` tests, not just graceful
shutdowns) for it specifically.

**What this is not, to keep it from blurring into the indexed-storage
result in the sibling document.** `piecewise_construction_and_trajectory
_caching.md` §6 measures how fast Resource B can be *queried* once it
exists (4,696x / 4,404x over a linear CSV scan, in exchange for ~6.1x more
storage) — and that document is explicit that this alone doesn't
differentiate from Syzygy, because Syzygy's own format is already
access-efficient. Claim C is a different, prior claim: not that the graph
is fast to query, but that it encodes something — family membership and
inter-family trajectories over escape-count ties — that a DTZ/DTM-only
format has no dimension to build such a graph from in the first place,
regardless of how fast either format's own lookups are.

**The honest gap.** What has *not* been directly measured is the specific
before/after comparison that would make this claim as airtight as Claim B's
byte-for-byte checks or the sibling document's 168ms-vs-0.04ms numbers: the
cost of deriving a given family's trajectory relationships on demand,
freshly, from Resource A alone, timed against reading the same answer from
an already-built Resource B. The checkpointing evidence above shows the
*construction* of Resource B is expensive enough to matter; it does not by
itself measure how expensive re-deriving one specific family-level answer
on the fly would be, each time, without Resource B. That controlled test —
analogous in spirit to the CSV-vs-index timing already run for raw access
speed — has not been run for this specific question and should be treated
as open, not assumed in this claim's favor.

---

## 7. Distribution and hosting — precedent, not a result

This section exists to close off one specific, reasonable objection to
Claims A/B/C together — "wouldn't a richer, relational dataset like this be
unwieldy to distribute?" — without overstating what's actually been done
about it, which is nothing yet.

**What's real:** Syzygy's own tables today are already served through
multiple independent, cooperating channels simultaneously — HTTP mirrors
(`tablebase.sesse.net`), a BitTorrent tracker, and separate API servers
(`lichess.org`'s public tablebase API, `syzygy-tables.info`), none of them
exclusive of the others. That's an existing, proven pattern for
distributing exactly this class of data — large, static, widely-reused
chess tablebase output — at real scale, today.

**What that does and doesn't mean for this project.** It means that *if*
Resource A and Resource B are ever distributed the way Syzygy's own tables
are, there's a known playbook to follow rather than a novel distribution
problem to solve from scratch — a de-risking point, not a differentiator.
It emphatically does not mean this project has built, tested, measured, or
demonstrated anything about hosting, mirroring, seeding, or serving its own
output. No file from this project has been distributed via any channel.
This section is a reason not to treat future distribution as a blocking
unknown — nothing more, and it should never be cited alongside Claims A, B,
or C's own verified numbers as though it were evidence of the same kind.

---

## 8. What remains open, stated as honestly as the rest of this project's documents

**On Claim B specifically:**
- **Scale.** This was one small pair (KPvK, KQvK, KNvK — hundreds of
  thousands of positions), not a real 6-7 piece hierarchy. Whether the same
  correctness and the same proportional storage savings hold at real scale
  is untested.
- **Systematic application.** This was one deliberately-constructed test,
  not a real, full dependency graph solved end-to-end with every material
  preloading from every simpler one it reduces into. `shape_based
  _compression_and_scoped_generation.md`'s own §4d recursive claim — the
  *whole-hierarchy* version, where no material anywhere in a chain stores
  its full table — has not been built or measured this way yet. This test
  is the first real, verified building block toward that, not the whole
  thing.
- **Does the origin/complement fraction actually shrink with complexity?**
  The project's own KQvK (76.6%) vs. KBPvK (69.0%) origin-fraction data
  point is suggestive, not confirmed as a trend — two materials, not enough
  to call a scaling law.
- **Whether this specific bug had already corrupted any existing output.**
  Any database built with `--preload-from` before the §4 fix should be
  treated as needing re-verification for its distance/escape values on
  positions that could plausibly have had a faster, not-yet-resolved
  native alternative at the time a preloaded value became available. This
  document does not audit any specific existing database for that.

**On Claim C specifically:**
- **The controlled re-derivation-cost comparison described in §6 hasn't
  been run.** This is the single biggest gap in Claim C as stated — the
  checkpointing evidence shows construction is expensive, but a direct,
  timed "recompute on demand vs. read from Resource B" test for a specific
  family-level question has not been done.
- **KBNvK's own Resource B is still being built as of this writing**, not
  yet complete or independently verified at that scale the way KQvK's was
  (3,601,748 node-shape rows, 412,152 edges, confirmed byte-for-byte
  identical across two genuine `SIGKILL`-and-resume tests).
- **No downstream consumer has actually queried Resource B for a real
  family/trajectory question yet.** The architectural claim — that this
  kind of question becomes a lookup instead of a fresh derivation — hasn't
  been exercised end-to-end by an actual use case, only by the
  infrastructure that would serve one.

**On distribution/hosting:**
- **Entirely unbuilt and untested for this project's own data.** §7's
  content is precedent from Syzygy's own infrastructure, not a result of
  any kind produced by this project. Nothing here should be read as
  evidence toward "improving on Syzygy" the way Claims A, B, or C are.

---

## 9. How to reproduce or extend the Claim B test

```bash
# Build the two dependencies
python3 generate_exhaustive_positions.py --white Q --out kqvk_seeds.txt
./general_solver --full-dag --positions kqvk_seeds.txt --db kqvk.db --fresh

python3 generate_exhaustive_positions.py --white N --out knvk_seeds.txt
./general_solver --full-dag --positions knvk_seeds.txt --db knvk.db --fresh

# Baseline: no preloading
python3 generate_exhaustive_positions.py --white P --out kpvk_seeds.txt
./general_solver --full-dag --positions kpvk_seeds.txt --db kpvk_baseline.db --fresh

# Preloaded: complement-only storage
./general_solver --full-dag --positions kpvk_seeds.txt --db kpvk_preloaded.db \
    --fresh --preload-from kqvk.db --preload-from knvk.db

# Verify exclusion (should be 0)
grep -c "^[QN]:" kpvk_preloaded.db

# Verify genuinely-new positions match exactly (should be empty)
grep "^P:" kpvk_baseline.db | sort > /tmp/base.txt
grep "^P:" kpvk_preloaded.db | sort > /tmp/pre.txt
diff /tmp/base.txt /tmp/pre.txt
```

To push Claim B further: pick a real three-level chain (e.g., something that
reduces into KPvK itself, matching the project's own KBPvK example) and
confirm the same exact-match and file-size results hold one level deeper —
this is the natural next test, not yet done.

To close Claim C's open gap (§8): pick one real, specific family-level
question answerable from Resource B (e.g. "which shapes does mate-group G's
worst-escape-count shape transition into"), time answering it by querying
the built Resource B, then separately time deriving the same answer by
walking Resource A directly with no precomputed graph available, and record
both numbers — the same shape of test as the CSV-vs-indexed-SQLite
comparison already run for raw access speed, but for the relational
question Claim C is actually about.

---

## 10. Relationship to this project's other documents

`shape_based_compression_and_scoped_generation.md` §4c/§4d describes the
architecture Claim B tests a real slice of — read that document for the
full design and its own, separately-tracked open questions.
`piecewise_construction_and_trajectory_caching.md` covers two related but
distinct things: constructing one material's *own* landscape piece by piece
(not reusing *other* materials' results, so complementary to Claim B rather
than overlapping with it), and, in its §6, the indexed-storage access-speed
result that this document's §6 explicitly distinguishes Claim C from. Both
sibling documents are where the fuller, less claim-focused design reasoning
lives; this document stays narrowly focused on what's been tested and what
that test did and didn't show.
