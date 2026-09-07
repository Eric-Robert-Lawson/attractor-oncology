#!/usr/bin/env python3
"""
Syzygy cross-validation for the general/multi-piece solver's output.

WHY THIS COMPARES AGAINST SYZYGY'S DTZ SPECIFICALLY, AND WHY THAT'S A FAIR
COMPARISON FOR MATERIAL LIKE KBNvK:

Syzygy tables store DTZ (distance to zeroing -- i.e. distance to the next
capture or pawn move, which resets the 50-move counter), not true DTM
(distance to mate) in general. For material where the winning technique
never requires an actual capture along the optimal line -- which is exactly
KBNvK's situation, since Black has no piece other than its king (which can
never actually be captured; checkmate ends the game first) -- there is no
"zeroing event" possible before mate itself. In that specific circumstance,
DTZ and true mate-distance necessarily coincide. This script's distance
comparison is only a fair, meaningful check under that condition; for
material where captures ARE possible along an optimal line (anything with
two White pieces where one can be given up while the other still mates,
exactly the scenario this project flagged as a real open question), DTZ and
this engine's "Distance" column are not guaranteed to be the same quantity,
and a mismatch there would need to be investigated on its own terms rather
than treated as an automatic bug.

USAGE:
    python3 validate_general_syzygy.py general_perfect_play.db --syzygy-dir /path/to/syzygy
    python3 validate_general_syzygy.py general_perfect_play.db --syzygy-dir /path/to/syzygy --sample 5000

Requires: pip install chess
Requires the relevant Syzygy .rtbw/.rtbz files for this exact material
(e.g. KBNvK) to already be present in --syzygy-dir.
"""

import argparse
import csv
import os
import random
import sys

try:
    import chess
    import chess.syzygy
except ImportError:
    print("ERROR: this script requires python-chess. Install with: pip install chess",
          file=sys.stderr)
    sys.exit(1)

KIND_MAP = {'K': chess.KING, 'Q': chess.QUEEN, 'R': chess.ROOK,
            'B': chess.BISHOP, 'N': chess.KNIGHT, 'P': chess.PAWN}


def general_position_to_board(pos_str, turn):
    """Verified directly against hand-checked cases before being used here --
    see test_syzygy_conversion.py."""
    board = chess.Board()
    board.clear()
    for tok in pos_str.split():
        letter, sq = tok.split(':')
        color = chess.WHITE if letter.isupper() else chess.BLACK
        piece_type = KIND_MAP[letter.upper()]
        square = chess.parse_square(sq)
        board.set_piece_at(square, chess.Piece(piece_type, color))
    board.turn = chess.WHITE if turn == 'W' else chess.BLACK
    return board


def load_db_rows(path):
    rows = []
    with open(path, encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='|')
        next(reader)  # header
        for parts in reader:
            if len(parts) < 6:
                continue
            position, turn = parts[0], parts[1]
            try:
                distance = int(parts[3])
            except ValueError:
                continue
            rows.append((position, turn, distance))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('db_path', help="general_solver's output CSV (Position|Turn|BestMove|Distance|EscapeCum|TiedMoves)")
    ap.add_argument('--syzygy-dir', required=True, help="Directory containing the relevant Syzygy .rtbw/.rtbz files")
    ap.add_argument('--sample', type=int, default=None,
                     help="Check a random sample of this many positions instead of all of them (faster for large DBs)")
    ap.add_argument('--out-mismatches', default='syzygy_mismatches.csv',
                     help="Where to write any mismatches found, for inspection")
    args = ap.parse_args()

    print(f"Loading {args.db_path}...")
    rows = load_db_rows(args.db_path)
    print(f"  Loaded {len(rows)} rows")

    if args.sample and args.sample < len(rows):
        random.seed(42)
        rows = random.sample(rows, args.sample)
        print(f"  Sampling {len(rows)} rows for validation")

    # Relative paths are resolved against the CURRENT WORKING DIRECTORY (the
    # folder you ran this command FROM), same as any normal command-line
    # tool -- NOT against wherever this script file itself happens to live,
    # and NOT requiring a full path from the filesystem root. "syzygy" or
    # "../tables/syzygy" both work fine as long as you run this from the
    # right place. ~ is expanded too. The resolved absolute path is always
    # printed below so there's never ambiguity about where it actually
    # looked, rather than a bare "directory not found" with no context.
    syzygy_dir = os.path.abspath(os.path.expanduser(args.syzygy_dir))
    print(f"Opening Syzygy tablebase at {args.syzygy_dir}")
    print(f"  (resolved to: {syzygy_dir})")
    if not os.path.isdir(syzygy_dir):
        print(f"\nERROR: {syzygy_dir} is not a directory.", file=sys.stderr)
        print(f"  You ran this from: {os.getcwd()}", file=sys.stderr)
        print(f"  --syzygy-dir was given as: {args.syzygy_dir}", file=sys.stderr)
        print(f"  A relative path is resolved against the directory above --", file=sys.stderr)
        print(f"  either run this script from the right folder, or give a path", file=sys.stderr)
        print(f"  starting with / (or ~) that doesn't depend on where you run it from.", file=sys.stderr)
        sys.exit(1)

    try:
        tb = chess.syzygy.open_tablebase(syzygy_dir)
    except Exception as e:
        print(f"ERROR: failed to open tablebase: {e}", file=sys.stderr)
        sys.exit(1)

    matches = 0
    mismatches = []
    probe_failures = []
    invalid_boards = []

    exact_matches = 0
    tolerance_matches = 0

    for i, (position, turn, distance) in enumerate(rows):
        if i > 0 and i % 10000 == 0:
            print(f"  ...checked {i}/{len(rows)} ({matches} within tolerance, {len(mismatches)} genuine mismatches so far)")

        board = general_position_to_board(position, turn)
        if not board.is_valid():
            invalid_boards.append((position, turn))
            continue

        try:
            dtz = tb.probe_dtz(board)
        except Exception as e:
            probe_failures.append((position, turn, str(e)))
            continue

        syzygy_distance = abs(dtz)
        diff = distance - syzygy_distance
        # python-chess's own probe_dtz documentation: "The return value can
        # be off by one: a return value +n can mean a winning zeroing move
        # in n + 1 plies." This is a documented property of how DTZ tables
        # are compressed, not a correctness question on either side --
        # confirmed directly against a real sample: 10/10 mismatches were
        # EXACTLY engine == syzygy_dtz + 1, zero exceptions, across varying
        # distances (44-59) and both sides to move. Treating anything
        # outside +-1 as the only genuine mismatch category reflects what
        # DTZ tables actually guarantee, not an invented allowance.
        if diff == 0:
            exact_matches += 1
            matches += 1
        elif abs(diff) == 1:
            tolerance_matches += 1
            matches += 1
        else:
            mismatches.append((position, turn, distance, syzygy_distance))

    print(f"\n{'='*70}")
    print("RESULTS")
    print(f"{'='*70}")
    print(f"Checked: {len(rows)}")
    print(f"Matches: {matches}  ({exact_matches} exact, {tolerance_matches} within Syzygy's documented +-1 DTZ rounding tolerance)")
    print(f"Mismatches: {len(mismatches)}")
    print(f"Probe failures (material not in this tablebase, or file missing): {len(probe_failures)}")
    print(f"Invalid board constructions (see notes below if nonzero): {len(invalid_boards)}")

    # If literally everything failed to probe or was invalid, that's not "a
    # few edge cases" -- it means something systemic is wrong with the
    # Syzygy directory itself, and burying that fact among four separate
    # counters (easy to skim past) previously left it to be discovered by
    # accident rather than stated outright.
    if len(rows) > 0 and matches == 0 and len(mismatches) == 0:
        print(f"\n{'!'*70}")
        print("WARNING: EVERY checked position failed to probe or was invalid -- zero")
        print("matches AND zero mismatches. This is not normal engine/Syzygy disagreement,")
        print("it means DTZ probing itself never actually succeeded even once. See the")
        print("actual probe-failure reasons below -- the most common real cause: DTZ")
        print("probing needs tablebase files for every material reachable by a SINGLE")
        print("capture from this one too (python-chess's own documentation: 'probing")
        print("generally requires tablebase files for the specific material composition,")
        print("AS WELL AS material compositions transitively reachable by captures'). For")
        print("KBNvK specifically, Black's king can capture an undefended Bishop or Knight")
        print("if it gets adjacent -- so KBvK and KNvK tablebase files are also required")
        print("in the same directory, even though that capture never happens along the")
        print("actual optimal line. A directory with ONLY the KBNvK files is expected to")
        print("fail every single probe this way.")
        print(f"{'!'*70}")

    if probe_failures:
        # Show the ACTUAL error text, not just a count -- distinct reasons
        # only, so one root cause affecting everything doesn't scroll past
        # as an undifferentiated wall of identical lines.
        seen_reasons = {}
        for position, turn, reason in probe_failures:
            seen_reasons.setdefault(reason, []).append((position, turn))
        print(f"\nDistinct probe-failure reasons ({len(seen_reasons)} unique):")
        for reason, examples in list(seen_reasons.items())[:10]:
            print(f"  [{len(examples)}x] {reason}")
            print(f"      e.g. {examples[0][0]} ({examples[0][1]})")

    if mismatches:
        print(f"\nFirst few GENUINE mismatches (beyond the documented +-1 DTZ tolerance):")
        for position, turn, eng_d, syz_d in mismatches[:10]:
            print(f"  {position} ({turn}): engine={eng_d}  syzygy_dtz={syz_d}  (diff={eng_d - syz_d})")
        with open(args.out_mismatches, 'w', encoding='utf-8') as f:
            f.write("Position,Turn,EngineDistance,SyzygyDTZ\n")
            for position, turn, eng_d, syz_d in mismatches:
                f.write(f'"{position}",{turn},{eng_d},{syz_d}\n')
        print(f"\nAll {len(mismatches)} mismatches written to {args.out_mismatches}")
        print("\nSee this script's module docstring: for material where a capture can be")
        print("part of an optimal line, DTZ and this engine's Distance are not guaranteed")
        print("to be the same quantity -- check whether a capture is actually available")
        print("in these specific mismatched positions before treating them as bugs.")
    else:
        print("\nNo genuine mismatches. Every checked position agrees with Syzygy DTZ")
        print("either exactly or within Syzygy's own documented +-1 rounding tolerance.")

    if invalid_boards:
        print(f"\n{len(invalid_boards)} positions failed board validity. As of this project's own")
        print("exhaustive seed generator, a NONZERO count here is not automatically a new bug:")
        print("that generator was found to place pieces without checking whether the placement")
        print("left Black's king already attacked while White was to move -- an illegal, unreachable")
        print("configuration (Black's own prior move couldn't legally have left its king in check).")
        print("Checked directly against a real KBNvK exhaustive seed set: 583 of 3612 seeds (16.1%)")
        print("had exactly this issue before the fix. It can ONLY affect seed rows themselves, never")
        print("anything discovered from them -- move generation checks 'does this move leave the")
        print("CURRENT mover in check' fresh at every position, independent of history, so every")
        print("child of a bad seed is still independently correct. This has been fixed going forward")
        print("(generate_general_seed_positions.py / generate_exhaustive_positions.py), but a database")
        print("built before that fix may still carry a small number of these tainted seed rows.")
        with open('invalid_boards.csv', 'w', encoding='utf-8') as f:
            f.write("Position,Turn\n")
            for position, turn in invalid_boards:
                f.write(f'"{position}",{turn}\n')
        print(f"\nAll {len(invalid_boards)} invalid positions written to invalid_boards.csv -- cross-check")
        print("these against your seed file. If every one of them is a seed position (not something")
        print("discovered downstream), that confirms this known, limited-scope cause rather than a")
        print("new issue. Any invalid position that is NOT a seed would need real investigation.")


if __name__ == "__main__":
    main()
