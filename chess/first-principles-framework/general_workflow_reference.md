# Chess Endgame Perfect-Play Workflow — Full Reference

Covers the complete pipeline from constructing a material's full landscape through analyzing it for multi-perfect-play findings. Current scope: **Black restricted to a bare king** (KXvK / KXYZ...vK, up to five non-king White pieces). This document is the reference to use until the multi-piece-Black refactor lands — see `loser_count_and_multipiece_black_refactor.md` for that separate track, which this document's §7 now links into directly.

---

## 0. One-time setup

Build the solver once per material-family variant you need. All variants come from the same source file.

```bash
# Single-attacker engines (only needed if working with these specific, older, single-piece tools)
g++ -std=c++17 -O2 -o kqvk_solver compositional_trajectory_solver_modular.cpp
g++ -std=c++17 -O2 -DPIECE_ROOK -o krvk_solver compositional_trajectory_solver_modular.cpp
g++ -std=c++17 -O2 -DPIECE_PAWN -o kpvk_solver compositional_trajectory_solver_modular.cpp

# General/multi-piece engine -- this is the one the rest of this document uses
g++ -std=c++17 -O2 -DPIECE_GENERAL -o general_solver compositional_trajectory_solver_modular.cpp
```

`general_solver` handles any White material combination — one king plus **up to five** other pieces (at most one pawn) — against a bare Black king. That's **7 total pieces on the board**, matching Syzygy's own tablebase convention, and the ceiling is a real bit-budget limit of the packed-state representation, not an arbitrary round number (see §5). Use it for everything below unless you specifically need one of the older single-piece tools.

Materials with **0-2** non-king White pieces are the most battle-tested (KQvK, KRvK, KBNvK, KBBvK have all been run end-to-end and cross-checked). **3-5** non-king pieces are supported by every tool in this pipeline and have been directly tested for correctness, but expect the reachable position graph to be dramatically larger — there's no disk-backed storage yet, so a genuinely large material may hit real memory limits regardless of how carefully it's swept. Test with a small `--max-nodes` cap before committing to a full run on anything past 2 pieces.

---

## 1. Full workflow for a new material (e.g. KBNvK)

### Step 1 — Generate the exhaustive seed list

```bash
python3 generate_exhaustive_positions.py --white B,N --out material_seeds.txt
```

- `--white` — comma-separated White piece letters, e.g. `Q`, `R`, `B,N`, `Q,P`, or any combination of up to 5 (e.g. `Q,R,B,N,P`).
- True combinatorial enumeration: every legal placement of every piece is tried, not just one placement per king pairing — for one non-king piece this is ~370-500K seeds; for two, ~24-30 million; for 3+, expect it to keep growing steeply (see §0's caution on scale).
- Each seed is checked for legality, including that Black's king isn't already attacked by a White piece before White's first move (a real bug this project found and fixed — see §6).

**If the material includes a bishop and no pawn, run this twice — once per bishop color, into the same database:**

```bash
python3 generate_exhaustive_positions.py --white B,N --out material_light.txt --bishop-color light
python3 generate_exhaustive_positions.py --white B,N --out material_dark.txt --bishop-color dark
```

A bishop's square color is conserved by every legal move — light-squared and dark-squared bishop positions are genuinely disconnected components of the state graph, proven directly (zero overlap across millions of explored positions). For **two** bishops, there are **three** distinct, disconnected cases, not two: `--bishop-color opposite` (one of each — the practically common configuration), or `light`/`dark` (a same-colored pair, requiring an underpromotion to reach in a real game).

**Running without `--bishop-color` at all is also a legitimate, deliberate choice now**, not just an oversight to catch — confirmed directly: discovery correctly finds each color's full component regardless of what other seeds share the file, and `--max-nodes` protects memory identically either way, since what's capped is discovered node count, not seed count. The tradeoff is practical, not correctness: an unsplit seed file for a two-piece material is roughly 24.5 million lines, taking ~14s and ~1.5GB just to load into the sweep script, versus doing two smaller, explicit passes. Either approach is safe; pick based on whether you want that load cost paid once (combined) or twice (split).

If the material has a pawn, no extra flag is needed — the generator automatically places it on rank 2 for forward-reachability. At most one pawn is supported.

**A pawn alongside a bishop breaks the clean light/dark disjointness the split above depends on.** The invariant only holds because a *non-promoting* bishop's square color is conserved by every move it makes — but a pawn can promote to a bishop, and a promoted bishop's color depends on which file it lands on, independent of the original bishop's color. So a "light" seed can, via promotion, reach a position with the original bishop on light and a newly-promoted one on dark — exactly the kind of position the split was never designed to reason about for a single bishop. This isn't a correctness risk: `discover()` seals already-known positions regardless of which sweep found them first, so any overlap between the "light" and "dark" runs just means the second one correctly recognizes some positions as already-proven rather than re-discovering them. It does mean the two runs for a bishop+pawn material likely aren't as cleanly disjoint as the naming suggests, so don't be surprised if the second half's discovery count looks smaller than a bishop-only material's would.

### Step 2 — Run the sweep

```bash
python3 run_full_sweep.py material_light.txt --db material_perfect_play.db --solver ./general_solver --fresh
python3 run_full_sweep.py material_dark.txt --db material_perfect_play.db --solver ./general_solver
```

(Omit `--bishop-color` entirely and just run once if the material has no bishop, or if you've deliberately chosen the combined, unsplit approach above.)

**This script was rewritten during this project and now works differently than earlier versions — if a command you've used before stops being recognized (e.g. `--batch-size` errors as unrecognized), you're looking at stale documentation or a stale copy of the file, not a bug.** What changed and why:

- The seed file is fed to the solver in large, **contiguous chunks** (default 200,000 lines, `--chunk-size` to change it) — known and unknown seeds mixed together freely. There's no separate Python-side pruning pass anymore; the C++ solver's own `discover()` seals already-known positions internally, in compiled code, far faster than a Python scan ever did.
- This exists because the actual dominant cost at scale wasn't discovery or classification — it was reload overhead. Measured directly on a real, ~22M-row in-progress KBNvK sweep: the C++ solver's own database reload took 68.1s, and the *old* script's Python-side prune scan added another 33.7s on top, **every single dispatch**, regardless of how little of that dispatch's work was actually new. Collapsing thousands of small dispatches into far fewer large ones is what actually fixes this — diluting a batch with known seeds was never the real problem.
- `--max-nodes` (or auto-detection from system memory if omitted) is the actual, direct memory-safety mechanism, and it's completely independent of chunk size — confirmed empirically by feeding a chunk that mixed known seeds with a seed guaranteed to trigger massive new discovery, capped low, and watching it stop cleanly at exactly the cap regardless of how many other seeds rode along. Raising `--chunk-size` well past the default is safe from a memory standpoint; the only real cost is a bigger one-time seed-file load.
- **Resumability changed too**: since there's no more per-seed prune scan to naturally pick up where a stopped run left off, progress through the seed *file* (not the database) is tracked in `<db>.sweep_progress`, a small JSON sidecar recording how many lines have been fed to the solver so far. Re-running the identical command resumes from that offset directly. If the positions file's name or line count doesn't match what's recorded, resuming is refused with an explicit prompt rather than silently sweeping from the wrong offset.
- `--fresh` only on the **very first** invocation for a database — wipes any existing file (and the `.sweep_progress` sidecar) after an explicit confirmation prompt if real proven positions already exist. Every subsequent invocation (dark half, or resuming after a stop) should **omit** `--fresh`.
- Progress output now reads `CHUNK N (... seed lines, known and unknown mixed) -- through X/TOTAL of the full list (Y%)`, not the older `BATCH N (... genuinely new seeds)` format.

**An older version of this script, using `--batch-size` (default 12) with Python-side pruning, is still functionally correct and safe** — it was used successfully to complete a real KBNvK sweep during this project by scaling `--batch-size` up progressively (1000, then 10000, then 100000), which works because `--max-nodes` protects memory independently of batch size there too. It's simply less efficient at scale, since it pays *two* reload costs per dispatch (its own Python prune-scan plus the C++ solver's reload) rather than the newer script's one. Either script produces an identical, correct database — this is a performance difference, not a correctness one.

**Discovery failing on `--max-nodes` is all-or-nothing for the whole invocation, not per-seed — confirmed directly.** `discover()` runs one shared frontier across every seed fed to it at once; if that combined discovery exceeds the cap, *nothing* from that invocation gets written, including any number of genuinely easy seeds that happened to be bundled alongside one hard one. Tested directly: two trivially-resolvable seeds bundled with one seed engineered to blow the cap produced no database file at all — not a partial one, none. This is the concrete reason batch/chunk size choice matters beyond just efficiency, and it generalizes to every material, not just a specific one:

- A **small** batch/chunk size (down to 1) isolates which specific seeds are genuinely hard, and protects every other seed's progress from being discarded alongside a hard one it happened to share an invocation with. Worth doing at the start of any material without an established track record for how well-behaved its state space is.
- A **large** batch/chunk size is more efficient (far fewer database reloads) but more exposed to this same failure — one hard seed anywhere in a 200,000-line default chunk discards the whole chunk's progress, easy seeds included.
- The practical approach that's worked in this project: start small, and **escalate gradually** (e.g. 1 → 20 → 200 → 2000 → 20000, not straight to a huge number) rather than jumping immediately to a large batch/chunk size on unfamiliar material. The signal to watch for when it's safe to jump larger isn't elapsed time or a round number of steps — it's the "already resolved" fraction of what's been scanned so far climbing consistently high (comfortably 90%+). That's what actually indicates the landscape is well-mapped enough that a big batch is unlikely to contain an undiscovered hard seed. This is exactly the situation that made KBNvK's 13.3%-scanned-but-21.9M-already-found jump straight to `--batch-size 100000` safe — the ratio, not the batch-size history, is what to trust.
- If a large batch/chunk ever does fail, the right recovery is to drop back down, not retry the same size — a failure means something specific in that range is the obstacle, and stepping back down is what lets you actually corner which seeds they are rather than repeatedly discarding the same batch's progress against an unknown cause.
- Materials combining a pawn with another piece are worth treating with extra caution here specifically: this project hit a real OOM testing Bishop+Pawn directly (21M+ nodes before running out), and an even earlier version of the engine was OOM-killed by the same material's state space (noted directly in the C++ source's own comments) — there's no reason to expect this is unique to that one combination, so any pawn-containing material without its own track record is worth the same gradual-escalation caution.

**This can be safely stopped (Ctrl+C) and resumed at any time** by re-running the exact same command (minus `--fresh`) — already-proven work (both wins and confirmed draws) is never redone.

### Step 3 — (Optional but recommended) Validate against Syzygy

```bash
python3 validate_general_syzygy.py material_perfect_play.db --syzygy-dir /path/to/syzygy --sample 20000
```

- Requires **both** `.rtbw` and `.rtbz` files for the material itself, **and** for every material reachable by a single capture (e.g. for KBNvK: also KBvK and KNvK) — DTZ probing structurally requires WDL data at every step, confirmed directly from python-chess's own source.
- A difference of exactly 1 ply between engine and Syzygy is expected and reported separately as "within tolerance" — this is a documented property of DTZ table compression, not a disagreement. Only differences of 2+ plies are flagged as genuine mismatches worth investigating.
- Only checks `Distance`, never move choice — Syzygy has no way to verify this project's tied-move/escape-count logic; that's the actual scope of what this project adds beyond what tablebases already provide.
- Omit `--sample` to check every row (slow on a large database) or use a smaller sample first to catch systemic setup problems quickly before committing to a full run.

### Step 4 — Classify

```bash
python3 general_analyzer.py classify material_perfect_play.db --out-dir material_analysis --max-memory-gb 10
```

- `--min-plies N` — skip positions shorter than this (rarely needed).
- `--max-positions N` — cap how many positions get classified (for a quick partial check).
- `--memo-cache-size N` (default 3000) — bounds the internal reachable-set cache; lower it if memory/swap becomes a problem on a large run, raise it if you have RAM to spare.
- **`--max-memory-gb N`** — new: checks the process's *actual* measured memory (not an estimate) every ~10 seconds, and if it reaches this many GB, checkpoints immediately and exits cleanly (exit code 2) rather than letting the OS start swapping. Re-running the identical command resumes from the checkpoint. Unset by default; worth setting explicitly on anything approaching KBNvK's scale or larger. Handles the real Linux-vs-macOS `ru_maxrss` unit difference (KB vs. bytes) internally, so the number you pass means the same thing on either platform.
- The underlying memory footprint per position was substantially reduced during this project without changing any output: the database's in-memory representation switched from a dict-per-position to a lighter namedtuple (~32% less memory for that structure alone), and classification results now drop fields (`children` always; `tied`/`components` for non-independent categories) that are never actually read back out for those categories — an ~80% reduction on the results list specifically. Both changes were verified to produce byte-for-byte identical output to the prior version on real data; this is a memory optimization, not a behavior change.
- Automatically checkpoints progress every 60 seconds to `<db_path>.classify_checkpoint.pkl`, and resumes from it automatically if the same command is re-run after an interruption. The checkpoint is deleted only after a fully successful run.
- Prints live progress every ~10 seconds (including current memory usage if `--max-memory-gb` is set), and separately warns when a single position has 10+ tied moves (these can legitimately take a while on their own — up to `n` reachable-set searches and `n²/2` pairwise checks for `n` tied moves).
- Produces `material_analysis/independent_findings_deduped.csv` — the input to the next step.

### Step 5 — Families

```bash
python3 general_analyzer.py families material_analysis/independent_findings_deduped.csv --out-dir material_families --render-trees material_perfect_play.db
```

- This is the stage that actually produces the "unique findings" count comparable to this project's earlier established KQvK figures — `classify`'s own count is a different, larger, pre-family-collapse number (see §6).
- `--render-trees DB_PATH` — also renders the full move tree for each family/unique finding, for inspection. This reloads the full database, so it benefits from the same memory optimization described in Step 4.
- `--max-lines N` (default 500) — caps how much tree output gets rendered per finding.
- **`--checkpoint-path FILE`** — new: this stage now checkpoints and resumes automatically too (previously it didn't, and an interruption meant starting over). Defaults to `<findings_csv>.families_checkpoint.pkl` if not given. Checkpointing happens at the granularity of one fully-completed group within a role-sweep, since a single large group's pairwise symmetry check — not a whole role — is where a genuinely slow stretch actually shows up. Just re-run the identical command to resume after an interruption.

### Step 6 — Reduce: finding genuinely irreducible perfect-play structure

```bash
python3 general_analyzer.py reduce material_analysis/independent_findings_deduped.csv material_perfect_play.db --out-dir material_reduction
```

This asks a different question than families' role-sweep, and doesn't depend on families having run at all — it reads classify's output directly. Where families asks "which findings recur under a fixed-role sweep," reduce asks whether a finding's own tied branches are actually as independent as they look, and whether the exhaustive landscape as a whole reduces to a compact catalog of genuinely distinct decisions.

**What this catalog does and doesn't claim, stated precisely because the distinction matters:** it is complete and exact for every position reachable via legal play from this material's own starting configurations — every sub-optimal branch, every promotion outcome, every capture-simplified descendant, not just the optimal-play skeleton. It is not an approximation of that reachability closure; it *is* that closure, deduplicated. It does **not** generalize beyond that closure to an unrelated material, and it does not yet recognize the same relative arrangement of pieces shifted to a different part of the board as the same shape — the symmetry check is literal board rotation/reflection (the same 8-fold group used everywhere else in this pipeline), not translation-invariant pattern matching. A genuinely portable catalog of principles — one that would recognize "this is the same kind of situation" in a position never literally searched — is a real, substantial next step past what's built here, not something that falls out automatically.

Two phases, run in sequence by a single invocation:

**Phase 1 — tie-branch reduction** (over classify's tied findings only):
- For every finding's own tied root branches, checks whether ONE consistent transform maps one branch's entire forced continuation onto another's, ply for ply — not just whether they eventually land on the same square configuration. This single test correctly catches three things that look different on the surface: literal downstream convergence, a root position whose own self-symmetry pairs up two of its tied moves, and two branches that never touch again but run in exact mirrored parallel toward two different mates. A genuine match gets collapsed into one branch-group — a real duplicate, safe to discard.
- Separately, checks whether each branch's forced continuation eventually reaches (up to symmetry) *another* finding's own root elsewhere in the same findings file. This is recorded as a subsumption edge, not collapsed away — the root-level tie that led there is still a real, distinct decision even if it isn't the final word; deleting it would erase exactly the thing worth preserving.
- `--max-seconds-per-finding N` (default 60, `0` disables) — a hard, per-finding wall-clock budget enforced via SIGALRM. Cost here is not uniform the way classify's is: a 2-way tie sitting at the front of a 70+-ply forced line can cost far more than a wide tie that resolves in a few plies, and this compounds badly under real memory pressure on a large database. A finding that hits the budget is recorded as genuinely unanalyzed, not silently guessed at — see `skipped_findings.csv`.
- `--retry-skipped SKIPPED_CSV` — re-analyze only the findings listed in a prior run's `skipped_findings.csv`, against the *same* findings_csv they came from. Pairs naturally with a much larger (or disabled) time budget, once a handful of genuinely hard findings aren't also blocking hundreds of thousands of ordinary ones.

**Phase 2 — full-landscape decision-shape catalog** (over the *entire* database, tied or not):
- Scans every position, canonicalizing each one's required move(s) under board symmetry. This is the part that covers the vast majority of any real landscape classify/families never touch at all, since only genuine multi-way ties ever made it into their output — a forced, single-answer position is just as much a real principle of perfect play as a tied one, it's the identical question with a tie-count of one.
- Produces `decision_shapes.csv` — one row per genuinely distinct shape, with its recurrence count and one example position. `--full-position-map` additionally writes the full `Position,Turn -> ShapeId` mapping for every position (real size at scale — off by default).
- `--skip-full-landscape` — run Phase 1 only.

**Checkpoint/resume, and a real bug this fixed:** an earlier version of this command deleted Phase 1's checkpoint immediately after the analysis finished but *before* its output files were actually written — an interruption anytime afterward, including mid-Phase-2, left no record Phase 1 had ever completed, so a re-run redid it from scratch every time. Confirmed directly and fixed: `reduce` now checks, before doing any work, whether Phase 1's three output files already exist with no lingering checkpoint — if so, it loads that data back instead of re-running potentially hours of analysis. `--force-rerun-phase1` overrides this if the underlying findings_csv or database genuinely changed since.

Phase 1 checkpoints after every finding's processing *fully* completes (never mid-finding, never before), on whichever of three triggers fires first:
- 60 seconds elapsed (`--checkpoint-path`, default `<findings_csv>.reduce_checkpoint.pkl`)
- `--checkpoint-every-findings` findings completed since the last save (default 200 — guarantees a predictable cadence even if every single finding finishes just under the time threshold, which the time trigger alone could never catch)
- the finding just completed alone took at least `--slow-finding-checkpoint-seconds` (default 5.0) — securing a hard-won result immediately rather than risking it to a kill shortly afterward

Live progress printing (every ~10s) names the actual finding currently in progress — its position and tie-width — specifically so "is this stuck, or grinding through one genuinely hard item" is answerable from the terminal, not something to guess at from a bare count. Phase 2 checkpoints the same way (time-based only, since it has no natural per-item boundary to trigger on) to `--checkpoint-path-landscape`, default `<db_path>.landscape_checkpoint.pkl`.

As with classify/families, checkpointing here only helps once the process is far enough along to be inside the relevant loop — it does nothing for an interruption during the initial database load itself, which is paid in full on every single resume regardless of any checkpoint. At real scale (tens of millions of positions), that load alone can be the dominant memory cost — see §5.

---

## 2. Quick reference — direct solver invocation (no sweep wrapper)

Useful for small materials, quick checks, or a handful of hand-picked seeds rather than the full exhaustive list.

```bash
python3 generate_general_seed_positions.py --white Q --out quick_seeds.txt --num-seeds 8
./general_solver --full-dag --positions quick_seeds.txt --db quick.db --fresh
```

- `generate_general_seed_positions.py` produces a small, deliberately-spread seed set rather than exhaustive king-pair coverage — relies on (verified, but not universally proven) connectivity, so prefer the exhaustive generator for anything you intend to trust as complete.
- Same `--bishop-color {light,dark,opposite}` flag is available here too, with the same three-case caveat for two bishops.
- `./general_solver` flags: `--full-dag` (required), `--positions FILE`, `--db FILE`, `--fresh`, `--max-nodes N`.

---

## 3. File outputs, per material

| File | What it is |
|---|---|
| `material_perfect_play.db` | Every proven win: `Position\|Turn\|BestMove\|Distance\|EscapeCum\|TiedMoves` |
| `material_perfect_play.db.draws` | Every proven draw (position + turn only) — used to skip re-proving draws on resumed runs |
| `material_perfect_play.db.sweep_progress` | Line-offset into the seed file the current sweep has reached — lets `run_full_sweep.py` resume without rescanning |
| `material_perfect_play.db.classify_checkpoint.pkl` | In-progress `classify` results — auto-deleted on successful completion |
| `<findings_csv>.families_checkpoint.pkl` | In-progress `families` role-sweep state — auto-deleted on successful completion |
| `material_analysis/independent_findings_deduped.csv` | Classify-stage independent findings, D4-orbit deduplicated |
| `material_analysis/independent_unverified.csv` | Positions where verification was incomplete (not findings — gaps, worth investigating if nonzero) |
| `material_analysis/all_classifications.csv` | Every classified position, all categories |
| `material_families/` | Families vs. genuinely unique findings, post-`classify` |
| `<findings_csv>.reduce_checkpoint.pkl` | In-progress `reduce` Phase 1 (tie-branch) state — auto-deleted only after its output files below are fully written |
| `<db_path>.landscape_checkpoint.pkl` | In-progress `reduce` Phase 2 (full-landscape) state — auto-deleted on successful completion |
| `material_reduction/branch_reduction.csv` | Per-finding raw vs. irreducible tied-branch counts, and which moves collapsed together |
| `material_reduction/skipped_findings.csv` | Findings that hit `--max-seconds-per-finding` — genuinely unanalyzed, not a zero result |
| `material_reduction/subsumption_edges.csv` | Which findings' branches lead into which other findings' roots, up to symmetry |
| `material_reduction/decision_shapes.csv` | The full-landscape catalog: every genuinely distinct decision-shape, tied or forced, with recurrence counts |
| `material_reduction/position_to_shape.csv` | Full per-position shape mapping — only written with `--full-position-map` (real size at scale) |

---

## 4. Command flag reference (quick lookup)

**`generate_exhaustive_positions.py`**
`--white LETTERS` (required, up to 5 non-king pieces) · `--out FILE` (required) · `--bishop-color {light,dark,opposite}`

**`generate_general_seed_positions.py`**
`--white LETTERS` (required, up to 5 non-king pieces) · `--out FILE` (required) · `--num-seeds N` (default 8) · `--bishop-color {light,dark,opposite}`

**`run_full_sweep.py`** *(rewritten — see Step 2 above if this looks unfamiliar)*
`positions_file` (positional) · `--db FILE` (required) · `--solver PATH` (required) · `--chunk-size N` (default 200,000) · `--max-nodes N` · `--tmp-dir DIR` (default `.sweep_batches`) · `--fresh`

**`general_solver`** (the compiled C++ binary)
`--full-dag` · `--positions FILE` · `--db FILE` · `--fresh` · `--max-nodes N`

**`validate_general_syzygy.py`**
`db_path` (positional) · `--syzygy-dir DIR` (required) · `--sample N` · `--out-mismatches FILE` (default `syzygy_mismatches.csv`)

**`general_analyzer.py classify`**
`db_path` (positional) · `--min-plies N` (default 0) · `--max-positions N` · `--memo-cache-size N` (default 3000) · `--max-memory-gb N` (new) · `--out-dir DIR` (default `general_tie_analysis`)

**`general_analyzer.py families`**
`findings_csv` (positional) · `--out-dir DIR` (default `general_families`) · `--render-trees DB_PATH` · `--max-lines N` (default 500) · `--checkpoint-path FILE` (new)

**`general_analyzer.py reduce`**
`findings_csv` (positional) · `db_path` (positional) · `--out-dir DIR` (default `general_reduction`) · `--force-rerun-phase1` · `--retry-skipped SKIPPED_CSV` · `--max-seconds-per-finding N` (default 60) · `--checkpoint-path FILE` · `--checkpoint-every-findings N` (default 200) · `--slow-finding-checkpoint-seconds N` (default 5.0) · `--checkpoint-path-landscape FILE` · `--skip-full-landscape` · `--full-position-map`

**`general_analyzer.py tree`**
`db_path` (positional) · `--position POS` · `--turn {W,B}` · `--out FILE` · `--findings-csv FILE` · `--top N` (default 10) · `--out-dir DIR` (default `general_trees`) · `--max-lines N` (default 500) · `--classifications FILE`

---

## 5. Scale: how far this actually goes

The packed-state representation the C++ engine uses internally has a hard, bit-budget-derived ceiling of **5 non-king White pieces** — confirmed by direct arithmetic (WK + BK + turn = 19 fixed bits, leaving exactly room for 5 more 9-bit piece slots in a 64-bit key, with nothing left over). That's 7 pieces total on the board, matching Syzygy's own convention — not a coincidence so much as two independent designs converging on the same practical bit-packing tradeoff.

That ceiling is about what the *data structure* can represent, not what's practically solvable on a given machine. 0-2 non-king pieces is thoroughly proven territory (KQvK, KRvK, KBNvK, KBBvK all completed and cross-validated). 3+ pieces work correctly through every tool in this pipeline — generation, sweeping, classification — but the reachable position graph grows steeply, and there's no disk-backed storage yet, so a genuinely large 3+ piece material may simply exceed available RAM regardless of how carefully it's swept. Treat anything past 2 pieces as "should work, confirm the scale directly with a bounded test first," not "already proven at this scale."

---

## 6. Known gotchas — worth remembering

- **`--bishop-color` is easy to forget.** Unlike before, forgetting it is no longer automatically a mistake to catch and rerun — running without it now correctly covers both colors in one pass, verified safe under the same `--max-nodes` memory protection. But it does mean a much larger seed file and load time than splitting deliberately, so know which tradeoff you're choosing rather than getting it by default.
- **`--fresh` will prompt before destroying real progress** (on both the seed-sweep database and, in the current `run_full_sweep.py`, its `.sweep_progress` sidecar too), but only if real content already exists — a genuine first run never prompts.
- **Draws live in a separate `.draws` file**, not the main database — don't be surprised the main `.db` file only ever contains wins.
- **The old exhaustive seed generator had a real, now-fixed bug**: about 16% of seeds could leave Black's king already in check before White's first move (an illegal, unreachable position). Fixed going forward; a database built before the fix may still carry a small number of tainted *seed* rows specifically (never anything discovered downstream from them — proven safe by how move generation works). Syzygy validation will flag these as "invalid board constructions" with a note distinguishing this known cause from a genuinely new problem.
- **The KBNvK "more than 2 disconnected components" question from earlier in this project is effectively resolved** by the move to true combinatorial exhaustive enumeration (§1 Step 1) — that concern came from the older, king-pairs-only seed generator relying on a connectivity assumption; true exhaustive enumeration tries every legal placement directly and doesn't depend on connectivity at all, so it can't miss a component regardless of how many there turn out to be.
- **`classify`'s raw count is not the number comparable to earlier project figures.** The "396" KQvK figure (and its currently-reproduced equivalent) is a post-`families` number. Comparing a `classify`-stage count against it directly is comparing two different pipeline stages. **Open, unresolved discrepancy worth flagging**: a from-scratch KQvK run through this current pipeline reproduces 1524 raw classify-stage findings (matching the previously-established figure exactly) but only 528 unique post-`families` findings, not the previously-reported 659 — confirmed this isn't caused by any change made during this project's memory-optimization work (the prior, unmodified code reproduces the identical 528 on identical input), so the cause is somewhere earlier and hasn't been tracked down yet.
- **Memory on a long `classify` run** — meaningfully reduced by default now (namedtuple positions, stripped unused result fields), and `--max-memory-gb` gives an explicit, measured safety net on top for large materials. `--memo-cache-size` remains available to lower further if still tight.
- **`run_full_sweep.py` changed its interface** (§1 Step 2) — if a command that used to work now errors on `--batch-size`, that's stale documentation or a stale file copy, not a regression. The older interface still works correctly if that's what you have; it's just less efficient at scale.
- **A single hard seed can silently discard an entire batch/chunk's worth of otherwise-easy progress** — confirmed directly, not assumed (§1 Step 2). Start small on any unfamiliar material and escalate batch/chunk size only once the "already resolved" ratio is consistently high, rather than guessing a large size up front.
- **`reduce` used to redo Phase 1 from scratch if interrupted anywhere during Phase 2** — a real, now-fixed bug: its checkpoint was deleted the moment Phase 1 finished, before Phase 1's own output files were written, so an interruption afterward left no record it had ever completed. Fixed directly: `reduce` now checks for Phase 1's already-complete output before doing any work, and loads it back instead of re-running. If a version predating this fix is still in use, an interruption during Phase 2 means Phase 1 silently reruns in full on the next invocation.

---

## 7. What this document does *not* cover

The multi-piece-Black refactor (letting Black have real material, and potentially win) is a separate, not-yet-built track — see `loser_count_and_multipiece_black_refactor.md` for the full plan, open questions, and what specifically needs to change in the engine before this workflow extends to that case. One concrete piece of groundwork from this project's White-side piece-count expansion (§0, §5) is now available to that track directly: the packed-state bit budget was confirmed to hold exactly 5 non-king pieces total *regardless* of how they end up distributed between White and Black in a future flexible-color redesign — so the piece-count ceiling itself won't need to move again once Black-side material is added, only how those 5 slots get allocated.
