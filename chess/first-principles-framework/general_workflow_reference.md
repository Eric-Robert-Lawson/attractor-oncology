# Chess Endgame Perfect-Play Workflow — Full Reference

Covers the complete pipeline from constructing a material's full landscape through analyzing it for multi-perfect-play findings. Current scope: **Black restricted to a bare king** (KXvK / KXYvK). This document is the reference to use until the multi-piece-Black refactor lands — see `loser_count_and_multipiece_black_refactor.md` for that separate track.

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

`general_solver` handles any White material combination (one king plus up to two other pieces) against a bare Black king. Use it for everything below unless you specifically need one of the older single-piece tools.

---

## 1. Full workflow for a new material (e.g. KBNvK)

### Step 1 — Generate the exhaustive seed list

```bash
python3 generate_exhaustive_positions.py --white B,N --out material_seeds.txt
```

- `--white` — comma-separated White piece letters, e.g. `Q`, `R`, `B,N`, `Q,P`.
- Produces one seed per legal White-king/Black-king pairing (3,612 of them) — an explicit, exhaustive coverage guarantee over king configurations, not a connectivity assumption.
- Each seed is checked for legality, including that Black's king isn't already attacked by a White piece before White's first move (a real bug this project found and fixed — see §5).

**If the material includes a bishop and no pawn, run this twice — once per bishop color, into the same database:**

```bash
python3 generate_exhaustive_positions.py --white B,N --out material_light.txt --bishop-color light
python3 generate_exhaustive_positions.py --white B,N --out material_dark.txt --bishop-color dark
```

A bishop's square color is conserved by every legal move — light-squared and dark-squared bishop positions are genuinely disconnected components of the state graph, proven directly (zero overlap across millions of explored positions). **Forgetting `--bishop-color` doesn't error — it silently produces a 50/50 mix of both colors in one file**, defeating the purpose of running them separately. Always double-check the flag is actually there before running a real sweep.

If the material has a pawn, no extra flag is needed — the generator automatically places it on rank 2 for forward-reachability.

### Step 2 — Run the sweep

```bash
python3 run_full_sweep.py material_light.txt --db material_perfect_play.db --solver ./general_solver --fresh
python3 run_full_sweep.py material_dark.txt --db material_perfect_play.db --solver ./general_solver
```

(Omit `--bishop-color` entirely and just run once if the material has no bishop.)

- `--fresh` only on the **very first** invocation for a database — wipes any existing file. If the file already has real proven positions in it, `--fresh` will prompt for explicit confirmation (`type 'delete'`) before wiping anything — this is deliberate, to prevent accidentally destroying completed work by re-running an old command from shell history.
- Every subsequent invocation (dark half, or resuming after a stop) should **omit** `--fresh`.
- The sweep automatically prunes any seed already proven on disk before ever invoking the solver, and only spends a real invocation once a full batch of genuinely new seeds has accumulated — so re-running the same command is always safe and cheap once the material is fully covered.
- Progress output shows exactly how far through the full seed list the scan has gotten (`scanned X/3612`), not just how many batches have run.
- `--batch-size N` (default 12) — seeds per batch.
- `--max-nodes N` — passed through to the solver if you need to override its auto-detected memory cap.

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
python3 general_analyzer.py classify material_perfect_play.db --out-dir material_analysis
```

- `--min-plies N` — skip positions shorter than this (rarely needed).
- `--max-positions N` — cap how many positions get classified (for a quick partial check).
- `--memo-cache-size N` (default 3000) — bounds the internal reachable-set cache; lower it if memory/swap becomes a problem on a large run, raise it if you have RAM to spare.
- Automatically checkpoints progress every 60 seconds to `<db_path>.classify_checkpoint.pkl`, and resumes from it automatically if the same command is re-run after an interruption. The checkpoint is deleted only after a fully successful run.
- Prints live progress every ~10 seconds, and separately warns when a single position has 10+ tied moves (these can legitimately take a while on their own — up to `n` reachable-set searches and `n²/2` pairwise checks for `n` tied moves).
- Produces `material_analysis/independent_findings_deduped.csv` — the input to the next step.

### Step 5 — Families

```bash
python3 general_analyzer.py families material_analysis/independent_findings_deduped.csv --out-dir material_families --render-trees material_perfect_play.db
```

- This is the stage that actually produces the "unique findings" count comparable to this project's earlier established KQvK figures — `classify`'s own count is a different, larger, pre-family-collapse number (see §5).
- `--render-trees DB_PATH` — also renders the full move tree for each family/unique finding, for inspection.
- `--max-lines N` (default 500) — caps how much tree output gets rendered per finding.

---

## 2. Quick reference — direct solver invocation (no sweep wrapper)

Useful for small materials, quick checks, or a handful of hand-picked seeds rather than the full exhaustive list.

```bash
python3 generate_general_seed_positions.py --white Q --out quick_seeds.txt --num-seeds 8
./general_solver --full-dag --positions quick_seeds.txt --db quick.db --fresh
```

- `generate_general_seed_positions.py` produces a small, deliberately-spread seed set rather than exhaustive king-pair coverage — relies on (verified, but not universally proven) connectivity, so prefer the exhaustive generator for anything you intend to trust as complete.
- Same `--bishop-color {light,dark}` flag is available here too.
- `./general_solver` flags: `--full-dag` (required), `--positions FILE`, `--db FILE`, `--fresh`, `--max-nodes N`.

---

## 3. File outputs, per material

| File | What it is |
|---|---|
| `material_perfect_play.db` | Every proven win: `Position\|Turn\|BestMove\|Distance\|EscapeCum\|TiedMoves` |
| `material_perfect_play.db.draws` | Every proven draw (position + turn only) — used to skip re-proving draws on resumed runs |
| `material_analysis/independent_findings_deduped.csv` | Classify-stage independent findings, D4-orbit deduplicated |
| `material_analysis/independent_unverified.csv` | Positions where verification was incomplete (not findings — gaps, worth investigating if nonzero) |
| `material_analysis/all_classifications.csv` | Every classified position, all categories |
| `material_families/` | Families vs. genuinely unique findings, post-`classify` |

---

## 4. Command flag reference (quick lookup)

**`generate_exhaustive_positions.py`**
`--white LETTERS` (required) · `--out FILE` (required) · `--bishop-color {light,dark}`

**`generate_general_seed_positions.py`**
`--white LETTERS` (required) · `--out FILE` (required) · `--num-seeds N` (default 8) · `--bishop-color {light,dark}`

**`run_full_sweep.py`**
`positions_file` (positional) · `--db FILE` (required) · `--solver PATH` (required) · `--batch-size N` (default 12) · `--max-nodes N` · `--tmp-dir DIR` (default `.sweep_batches`) · `--fresh`

**`general_solver`** (the compiled C++ binary)
`--full-dag` · `--positions FILE` · `--db FILE` · `--fresh` · `--max-nodes N`

**`validate_general_syzygy.py`**
`db_path` (positional) · `--syzygy-dir DIR` (required) · `--sample N` · `--out-mismatches FILE` (default `syzygy_mismatches.csv`)

**`general_analyzer.py classify`**
`db_path` (positional) · `--min-plies N` (default 0) · `--max-positions N` · `--memo-cache-size N` (default 3000) · `--out-dir DIR` (default `general_tie_analysis`)

**`general_analyzer.py families`**
`findings_csv` (positional) · `--out-dir DIR` (default `general_families`) · `--render-trees DB_PATH` · `--max-lines N` (default 500)

**`general_analyzer.py tree`**
`db_path` (positional) · `--position POS` · `--turn {W,B}` · `--out FILE` · `--findings-csv FILE` · `--top N` (default 10) · `--out-dir DIR` (default `general_trees`) · `--max-lines N` (default 500) · `--classifications FILE`

---

## 5. Known gotchas — worth remembering

- **`--bishop-color` is easy to forget and fails silently.** Always verify it's actually in the command before running a real sweep on bishop-containing material.
- **`--fresh` will prompt before destroying real progress**, but only if the database already exists with content — a genuine first run never prompts.
- **Draws live in a separate `.draws` file**, not the main database — don't be surprised the main `.db` file only ever contains wins.
- **The old exhaustive seed generator had a real, now-fixed bug**: about 16% of seeds could leave Black's king already in check before White's first move (an illegal, unreachable position). Fixed going forward; a database built before the fix may still carry a small number of tainted *seed* rows specifically (never anything discovered downstream from them — proven safe by how move generation works). Syzygy validation will flag these as "invalid board constructions" with a note distinguishing this known cause from a genuinely new problem.
- **KBNvK showed evidence of more than 2 disconnected components** (beyond the known light/dark bishop split) during the exhaustive sweep — a real, still-open structural question, not yet fully explained. Worth watching for on any material with a bishop and a knight together.
- **`classify`'s raw count is not the number comparable to earlier project figures.** The "396" KQvK figure (and its correct current-data equivalent, 659) is a post-`families` number. Comparing a `classify`-stage count against it directly is comparing two different pipeline stages.
- **Memory/swap on a long `classify` run** — bounded now via `--memo-cache-size`, but lower it further if still tight on a particular machine.

---

## 6. What this document does *not* cover

The multi-piece-Black refactor (letting Black have real material, and potentially win) is a separate, not-yet-built track — see `loser_count_and_multipiece_black_refactor.md` for the full plan, open questions, and what specifically needs to change in the engine before this workflow extends to that case.
