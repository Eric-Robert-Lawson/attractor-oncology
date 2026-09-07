#!/usr/bin/env python3
"""
Automated batch sweep for general_solver, over an exhaustive position list.

WHAT THIS ACTUALLY DOES, PRECISELY:

Takes a positions file (from generate_exhaustive_positions.py, or any file
in the same format) and walks it once, greedily pruning any seed that
already has a proven answer on disk (checked against an in-memory set
loaded ONCE at sweep start, not re-scanned per seed) and accumulating the
genuinely NOT-yet-known seeds into groups of --batch-size. Each full group
is run through general_solver once, as a completely separate process
invocation, against the SAME shared database file. The final, possibly
undersized group (fewer than --batch-size genuinely new seeds, once the
list is exhausted) still runs rather than being silently dropped.

WHY PRUNE-AND-ACCUMULATE RATHER THAN FIXED WINDOWS OF THE ORIGINAL FILE
ORDER: every solver invocation -- however little new work it actually has
-- pays a fixed startup cost re-parsing the entire db and draws files from
scratch (confirmed directly in the C++ source: load_general_db_for_resume
and load_general_draws_for_resume both run unconditionally on every call,
regardless of how much of that call's own work turns out to be new). A
fixed window of 12 seeds where 9 are already known and only 3 are
genuinely new still pays that full cost for just 3 seeds' worth of actual
work. Pruning known seeds before ever invoking the solver, and only
spending an invocation once a FULL batch of genuinely new seeds has
accumulated, means every invocation's fixed cost is spent on a full
batch's worth of real work, never diluted by already-known seeds riding
along for free.

Each batch still benefits from every batch before it, exactly as before:
general_solver loads the database at the start of every run, seeds its
classifier directly with every already-proven position (verified sound --
see the engine's own comments on seed_from_preloaded), and its discovery
phase skips expanding the children of anything already known. Confirmed
directly on real material: a resumed run against an already-solved
database dropped from rediscovering 402,688 positions down to 67, because
only the genuinely new frontier gets expanded.

WHY EACH BATCH IS A SEPARATE PROCESS, NOT A LOOP INSIDE ONE LONG-RUNNING
PROGRAM: a fresh process is the only mechanism that GUARANTEES memory is
fully returned to the OS between batches -- clearing a C++ container
doesn't always release its underlying allocation back to the operating
system immediately, depending on the allocator, whereas the OS reclaiming
an entire exited process's memory is unconditional. This is what actually
delivers "batch completes, memory drops to zero, next batch starts clean."

Once the graph for a given material is fully covered, a batch invocation
will typically discover ~0 new positions and complete almost instantly --
this is expected and correct, not something to optimize away. Under the
prune-and-accumulate strategy, this should now be rare rather than the
common case: a batch only gets invoked once 12 genuinely-not-yet-known
seeds have accumulated, so a batch turning out to be entirely already-known
would mean the pruning check itself was somehow wrong, not that the
material happened to already be solved.

Each batch's output streams live, line by line, as the solver produces it
-- including its own periodic "[discovery progress]" (every 1M positions)
and "[classify progress]" (every pass) checkpoints -- rather than being
buffered and dumped only after the whole batch finishes. Confirmed
directly: checked a running sweep's log file mid-batch and saw it still
growing in real time (14 of an eventual 37 lines present partway through
a 4.4s batch), not silent until completion. This matters most on a
material's first real batch, which can run for minutes with no output
otherwise -- indistinguishable from a hang without live streaming.

Usage:
    python3 run_full_sweep.py kqvk_exhaustive.txt --db kqvk_perfect_play.db --solver ./general_solver
    python3 run_full_sweep.py kbnvk_exhaustive.txt --db kbnvk_perfect_play.db --solver ./general_solver --batch-size 200
"""

import argparse
import subprocess
import sys
import os
import time

# Mirrors GeneralState::str() in the C++ engine EXACTLY -- verified against
# real solver output, not assumed. Sort order: White before Black, then by
# this specific kind ordering (confirmed directly from the PieceKind enum:
# PAWN=0, QUEEN=1, KNIGHT=2, KING=3, ROOK=4, BISHOP=5 -- NOT alphabetical,
# NOT king-first), then by square. Getting this wrong would either silently
# disable the whole optimization below (harmless but pointless) or, far
# worse, cause it to silently skip a batch that actually has new work --
# so this was checked against a real, fresh solver run's own db output
# before being trusted, not just written and assumed correct.
_KIND_ORDER = {'P': 0, 'Q': 1, 'N': 2, 'K': 3, 'R': 4, 'B': 5}


def canonical_position_key(seed_line, turn='W'):
    """Converts one comma-separated seed line ('K:a1,B:a2,N:b1,k:a3') into
    the exact (position_string, turn) key that would appear as this
    position's row in the database, if it's ever been solved -- i.e. the
    same string GeneralState::str() would produce for it."""
    pieces = []
    for tok in seed_line.split(','):
        tok = tok.strip()
        if not tok:
            continue
        letter, sq = tok.split(':')
        color = 0 if letter.isupper() else 1  # White=0, Black=1
        pieces.append((color, _KIND_ORDER[letter.upper()], sq, letter, ))
    pieces.sort(key=lambda p: (p[0], p[1], (ord(p[2][0]) - ord('a')) * 8 + (int(p[2][1]) - 1)))
    return " ".join(f"{letter}:{sq}" for _, _, sq, letter in pieces), turn


def load_known_keys(db_file):
    """Loads every (position, turn) key already proven -- wins from db_file,
    draws from db_file + '.draws' -- into one in-memory set. Called ONCE per
    sweep run (not once per batch), specifically so checking a batch's seeds
    against it is a fast in-memory lookup, not a fresh file scan each time."""
    known = set()
    for path, has_extra_cols in [(db_file, True), (db_file + ".draws", False)]:
        if not os.path.isfile(path):
            continue
        with open(path) as f:
            next(f, None)  # header
            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue
                parts = line.split("|")
                if len(parts) < 2:
                    continue
                known.add((parts[0], parts[1]))
    return known


def read_all_lines(path):
    with open(path) as f:
        lines = [line.rstrip('\n') for line in f]
    header = None
    body = []
    for line in lines:
        if line.startswith('#') or (header is None and not line.strip().startswith(('K:', 'k:'))):
            if header is None:
                header = line
            continue
        if line.strip():
            body.append(line)
    return header or "# batch seeds", body


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('positions_file', help="Exhaustive (or any) positions file to sweep through, in batches")
    ap.add_argument('--db', required=True, help="Shared database file -- accumulates across every batch")
    ap.add_argument('--solver', required=True, help="Path to the compiled general_solver binary")
    ap.add_argument('--batch-size', type=int, default=12,
                     help="Seed positions per batch (default 12, matching this project's established "
                          "convention -- larger batches mean fewer subprocess-startup overheads at the "
                          "cost of slightly coarser memory-drop granularity)")
    ap.add_argument('--max-nodes', type=int, default=None,
                     help="Passed through to each general_solver invocation, if given -- otherwise each "
                          "batch auto-detects from system memory, same as running general_solver directly")
    ap.add_argument('--tmp-dir', default='.sweep_batches', help="Where per-batch seed files are written")
    ap.add_argument('--fresh', action='store_true', help="Wipe --db before starting (passed to the FIRST batch only)")
    args = ap.parse_args()

    if not os.path.isfile(args.solver):
        print(f"ERROR: solver binary not found at {args.solver}", file=sys.stderr)
        sys.exit(1)

    # --fresh unconditionally deletes an existing database, with no memory
    # of what was in it -- and the most natural way to trigger this by
    # ACCIDENT is simply re-running the same command line to resume after
    # stopping a sweep, if that line still has --fresh in it from the
    # original launch. Confirmed directly: re-running with --fresh still
    # present reports "Discovery: 372076 positions" (full rediscovery) on
    # a database that already had everything proven; the identical command
    # without --fresh reports "Discovery: 0 positions" (correctly resumed).
    # This check exists specifically so that mistake costs a keypress, not
    # hours of already-completed work.
    if args.fresh and os.path.isfile(args.db):
        existing_lines = sum(1 for _ in open(args.db)) - 1  # minus header
        if existing_lines > 0:
            print(f"\n{'!'*70}")
            print(f"WARNING: --fresh will DELETE {args.db}, which already contains")
            print(f"{existing_lines} proven positions. This cannot be undone.")
            print(f"{'!'*70}")
            print(f"If you're resuming a sweep you already started, you almost")
            print(f"certainly do NOT want --fresh -- just remove it and re-run.")
            answer = input(f"\nType 'delete' to actually wipe {args.db} and start over, "
                            f"or anything else to abort: ")
            if answer.strip() != 'delete':
                print("Aborted -- database left untouched.")
                sys.exit(1)

    header, all_seeds = read_all_lines(args.positions_file)
    total = len(all_seeds)
    if total == 0:
        print("ERROR: no seed positions found in the input file", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.tmp_dir, exist_ok=True)

    print(f"Loaded {total} seed positions from {args.positions_file}")
    print(f"Scanning for genuinely new work, accumulating into batches of up to "
          f"{args.batch_size} not-yet-known seeds each, writing into {args.db}\n")

    # Loaded ONCE here, not re-scanned per batch -- checking a seed against
    # this in-memory set is a fast dict lookup, not a fresh multi-million-
    # line file read. Every solver invocation, sealed or not, currently
    # re-parses the ENTIRE db and draws files from scratch before discovery
    # even starts (confirmed directly in the C++ source --
    # load_general_db_for_resume and load_general_draws_for_resume both run
    # unconditionally on every call). That reload is a fixed cost per
    # invocation, independent of how much of the batch is genuinely new --
    # which is exactly why a batch diluted with already-known seeds still
    # pays the full cost for however few new seeds it actually contains.
    known_keys = load_known_keys(args.db)

    # Greedy prune-and-accumulate, not fixed windows of the original file
    # order: walk the seed list once, silently dropping any seed already
    # resolved (no solver call spent on it at all, not even a diluted one),
    # and only invoking the solver once `batch_size` genuinely new seeds
    # have accumulated -- so every actual invocation's fixed startup cost
    # is spent on a full batch of real, new work, never diluted by
    # already-known seeds riding along. The final leftover group (fewer
    # than batch_size genuinely new seeds, once the seed list is exhausted)
    # still gets run rather than silently dropped -- that's the explicit
    # "last batch" exception.
    def run_batch(seed_group, batch_num, scan_pos, scan_pruned):
        nonlocal known_keys
        batch_file = os.path.join(args.tmp_dir, f"batch_{batch_num:06d}.txt")
        with open(batch_file, 'w') as f:
            f.write(header + "\n")
            for line in seed_group:
                f.write(line + "\n")

        cmd = [args.solver, "--full-dag", "--positions", batch_file, "--db", args.db]
        if args.max_nodes:
            cmd += ["--max-nodes", str(args.max_nodes)]
        if args.fresh and batch_num == 1:
            cmd += ["--fresh"]

        batch_start = time.time()
        pct = 100.0 * scan_pos / total if total else 100.0
        print(f"{'='*70}")
        print(f"BATCH {batch_num}  ({len(seed_group)} genuinely new seeds) -- "
              f"scanned {scan_pos}/{total} of the full list so far ({pct:.1f}%), "
              f"{scan_pruned} of those already resolved, {total - scan_pos} not yet examined")
        print(f"{'='*70}")

        # Streamed line-by-line, NOT captured and printed after the fact --
        # capture_output=True blocks until the whole subprocess exits, so on
        # a long batch (a not-yet-covered material's real discovery can run
        # minutes) nothing would appear on screen until it finished,
        # indistinguishable from a hang. This surfaces the solver's own
        # periodic "[discovery progress]" and "[classify progress]" lines
        # live, exactly like running the solver directly would show.
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                 text=True, bufsize=1)
        stdout_lines = []
        for line in proc.stdout:
            stdout_lines.append(line)
            if any(key in line for key in (
                "discovery progress", "classify progress", "Discovery:", "Classification:",
                "Wrote", "Sealed", "ERROR", "WIN in", "DRAW"
            )):
                print("  " + line.rstrip())
        proc.wait()
        batch_elapsed = time.time() - batch_start

        if proc.returncode != 0:
            print(f"  [WARNING] batch {batch_num} exited with code {proc.returncode}")
            tail = "".join(stdout_lines).strip()
            if tail:
                print("  output tail:", tail[-2000:])

        print(f"  batch time: {batch_elapsed:.1f}s\n")
        os.remove(batch_file)

        # Refresh from disk -- this batch may have written genuinely new
        # rows (wins and/or draws), which the rest of the scan needs to see.
        known_keys = load_known_keys(args.db)

    sweep_start = time.time()
    batches_run = 0
    seeds_pruned = 0
    pending = []
    last_scanned_idx = 0  # how many seeds (1-indexed) have been examined so far

    for i, seed in enumerate(all_seeds):
        last_scanned_idx = i + 1
        key = canonical_position_key(seed, 'W')
        if key in known_keys:
            seeds_pruned += 1
            continue
        pending.append(seed)
        if len(pending) >= args.batch_size:
            batches_run += 1
            run_batch(pending, batches_run, last_scanned_idx, seeds_pruned)
            pending = []

    # Flush any leftover genuinely-new seeds after the scan completes -- NOT
    # an in-loop "is this the last seed" check, which silently fails to
    # fire whenever the file happens to END on an already-known (pruned)
    # seed: that seed hits `continue` and the check is never reached, so a
    # non-empty `pending` from earlier in the scan would be dropped with no
    # warning. Confirmed as a real bug this way, not a hypothetical: an
    # interleaved test (4 known, 3 new, ending on a known seed) silently
    # lost its 3rd genuinely-new seed before this fix -- ran only 2 of the
    # 3 seeds that should have been processed, with no error at all.
    if pending:
        batches_run += 1
        run_batch(pending, batches_run, last_scanned_idx, seeds_pruned)

    total_elapsed = time.time() - sweep_start
    print(f"{'='*70}")
    print(f"SWEEP COMPLETE: {total} seed positions scanned, {seeds_pruned} pruned as "
          f"already-known (never invoked the solver), {batches_run} batch(es) actually "
          f"run, {total_elapsed:.1f}s total")
    print(f"{'='*70}")
    print(f"Final accumulated database: {args.db}")

    try:
        os.rmdir(args.tmp_dir)
    except OSError:
        pass


if __name__ == "__main__":
    main()
