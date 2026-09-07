#!/usr/bin/env python3
"""
Automated batch sweep for general_solver, over an exhaustive position list.

WHAT THIS ACTUALLY DOES, PRECISELY:

Takes a positions file (from generate_exhaustive_positions.py, or any file
in the same format), slices it into batches of --batch-size seed positions,
and runs general_solver once PER BATCH, as a completely separate process
invocation, against the SAME shared database file.

Each batch benefits from every batch before it: general_solver loads the
database at the start of every run, seeds its classifier directly with
every already-proven position (verified sound -- see the engine's own
comments on seed_from_preloaded), and its discovery phase skips expanding
the children of anything already known, rather than rediscovering the
whole graph each time. Confirmed directly on real material: a resumed run
against an already-solved database dropped from rediscovering 402,688
positions down to 67, because only the genuinely new frontier gets
expanded.

WHY EACH BATCH IS A SEPARATE PROCESS, NOT A LOOP INSIDE ONE LONG-RUNNING
PROGRAM: a fresh process is the only mechanism that GUARANTEES memory is
fully returned to the OS between batches -- clearing a C++ container
doesn't always release its underlying allocation back to the operating
system immediately, depending on the allocator, whereas the OS reclaiming
an entire exited process's memory is unconditional. This is what actually
delivers "batch completes, memory drops to zero, next batch starts clean."

Once the graph for a given material is fully covered, later batches will
typically discover ~0 new positions and complete almost instantly -- this
is expected and correct, not something to optimize away: it's the direct
result of everything already being proven and sealed from disk.

Usage:
    python3 run_full_sweep.py kqvk_exhaustive.txt --db kqvk_perfect_play.db --solver ./general_solver
    python3 run_full_sweep.py kbnvk_exhaustive.txt --db kbnvk_perfect_play.db --solver ./general_solver --batch-size 200
"""

import argparse
import subprocess
import sys
import os
import time


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

    header, all_seeds = read_all_lines(args.positions_file)
    total = len(all_seeds)
    if total == 0:
        print("ERROR: no seed positions found in the input file", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.tmp_dir, exist_ok=True)
    num_batches = (total + args.batch_size - 1) // args.batch_size

    print(f"Loaded {total} seed positions from {args.positions_file}")
    print(f"Running {num_batches} batches of up to {args.batch_size} seeds each, "
          f"accumulating into {args.db}\n")

    sweep_start = time.time()
    for batch_idx in range(num_batches):
        start = batch_idx * args.batch_size
        end = min(start + args.batch_size, total)
        batch_seeds = all_seeds[start:end]

        batch_file = os.path.join(args.tmp_dir, f"batch_{batch_idx:06d}.txt")
        with open(batch_file, 'w') as f:
            f.write(header + "\n")
            for line in batch_seeds:
                f.write(line + "\n")

        cmd = [args.solver, "--full-dag", "--positions", batch_file, "--db", args.db]
        if args.max_nodes:
            cmd += ["--max-nodes", str(args.max_nodes)]
        if args.fresh and batch_idx == 0:
            cmd += ["--fresh"]

        batch_start = time.time()
        print(f"{'='*70}")
        print(f"BATCH {batch_idx + 1}/{num_batches}  (seeds {start}-{end - 1} of {total})")
        print(f"{'='*70}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        batch_elapsed = time.time() - batch_start

        # Surface the parts of the solver's own output that actually matter,
        # rather than a wall of text per batch -- most batches after the
        # first few will be nearly silent (0 new positions), which is the
        # expected, correct outcome once the graph is fully sealed.
        for line in result.stdout.splitlines():
            if any(key in line for key in ("Discovery:", "Classification:", "Wrote", "Sealed", "ERROR", "WIN in", "DRAW")):
                print("  " + line)
        if result.returncode != 0:
            print(f"  [WARNING] batch {batch_idx + 1} exited with code {result.returncode}")
            if result.stderr.strip():
                print("  stderr:", result.stderr.strip()[:2000])

        print(f"  batch time: {batch_elapsed:.1f}s\n")

        os.remove(batch_file)

    total_elapsed = time.time() - sweep_start
    print(f"{'='*70}")
    print(f"SWEEP COMPLETE: {num_batches} batches, {total} total seed positions, "
          f"{total_elapsed:.1f}s total")
    print(f"{'='*70}")
    print(f"Final accumulated database: {args.db}")

    try:
        os.rmdir(args.tmp_dir)
    except OSError:
        pass


if __name__ == "__main__":
    main()
