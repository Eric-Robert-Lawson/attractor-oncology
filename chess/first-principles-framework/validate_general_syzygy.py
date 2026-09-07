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

    print(f"Opening Syzygy tablebase at {args.syzygy_dir}...")
    try:
        tb = chess.syzygy.open_tablebase(args.syzygy_dir)
    except Exception as e:
        print(f"ERROR: failed to open tablebase: {e}", file=sys.stderr)
        sys.exit(1)

    matches = 0
    mismatches = []
    probe_failures = []
    invalid_boards = []

    for i, (position, turn, distance) in enumerate(rows):
        if i > 0 and i % 10000 == 0:
            print(f"  ...checked {i}/{len(rows)} ({matches} matches, {len(mismatches)} mismatches so far)")

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
        if syzygy_distance == distance:
            matches += 1
        else:
            mismatches.append((position, turn, distance, syzygy_distance))

    print(f"\n{'='*70}")
    print("RESULTS")
    print(f"{'='*70}")
    print(f"Checked: {len(rows)}")
    print(f"Matches: {matches}")
    print(f"Mismatches: {len(mismatches)}")
    print(f"Probe failures (material not in this tablebase, or file missing): {len(probe_failures)}")
    print(f"Invalid board constructions (should be 0 -- a real bug if not): {len(invalid_boards)}")

    if mismatches:
        print(f"\nFirst few mismatches (engine_distance vs syzygy_dtz):")
        for position, turn, eng_d, syz_d in mismatches[:10]:
            print(f"  {position} ({turn}): engine={eng_d}  syzygy_dtz={syz_d}")
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
        print("\nNo mismatches. All checked positions agree with Syzygy DTZ exactly.")

    if invalid_boards:
        print(f"\nWARNING: {len(invalid_boards)} positions failed board validity -- this needs")
        print("investigation, it means the position string itself may be malformed.")

    if probe_failures:
        print(f"\n{len(probe_failures)} positions could not be probed -- most likely cause: the")
        print("Syzygy directory doesn't have the tablebase file for this exact material,")
        print("or the material doesn't match what's in --syzygy-dir at all.")


if __name__ == "__main__":
    main()
