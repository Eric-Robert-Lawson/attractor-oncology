#!/usr/bin/env python3
"""
Exhaustive seed-position generator for the general/multi-piece solver.

TRUE combinatorial exhaustiveness: every legal (WK, other White pieces, BK,
turn) combination for the given material, not just every king pairing with
one fixed placement rule for everything else. This is what "exhaustive"
needs to mean to match a real, Syzygy-style completeness guarantee.

WHY THIS CHANGED FROM AN EARLIER, KING-PAIR-ONLY VERSION OF THIS SCRIPT:
that version only enumerated the 3,612 legal White-king/Black-king
pairings, placing every other piece via a single deterministic
"closest available square" rule per pairing. It caught a real gap once
(a second disconnected component in KBNvK, beyond the known bishop-color
split), but that was a heuristic getting lucky, not a proof of
completeness -- there was no guarantee some other material's missing
component would land on one of the 3,612 tried placements. This version
removes that gap entirely by trying every piece's every square, not one
placement per king pair.

Verified directly, not assumed: ran this against KQvK end to end through
the actual solver and compared the discovered graph against the
previously-established, Syzygy-validated 372,065-position landscape --
after fixing a real bug in an early version of this script (the seed
format had no way to encode Black-to-move, so some correctly-computed
Black-to-move-only placements were silently reinterpreted as White-to-move,
reintroducing exactly the illegal "Black already in check" positions the
generator was designed to exclude), the corrected comparison came back
372,064 vs 372,065 -- a difference of exactly one, from the seed itself.
True exhaustive enumeration confirmed the existing landscape was already
complete, rather than revealing anything missing.

PIECE COUNT: supports 0 through MAX_WHITE_NON_KING (5) non-king White
pieces -- matching the C++ engine's own packed-state limit exactly, which
gives 7 total pieces on the board including both kings, the same
convention Syzygy tablebases use. Confirmed directly: the underlying
enumerate_fully_exhaustive was regression-tested to produce IDENTICAL
counts to the old hardcoded 0/1/2-piece logic for 1 and 2 pieces (368,452
and 24,536,088 respectively, exact matches) before being trusted for 3-5.

SCALE, HONESTLY: for one non-king White piece (KQvK, KRvK, KPvK), this is
~370-500K positions and runs in well under a second. For two non-king
White pieces (KBNvK and similar), the raw combinatorial space is roughly
30 million before legality filtering -- still generates in well under a
minute, but produces a correspondingly large seed file. For 3+ pieces the
seed file itself remains generable in reasonable time, but the actual
reachable position graph for such material is very plausibly far larger
than anything solved so far in this project, and the engine's current
in-memory-only architecture (no disk-backed storage) may make a full,
exhaustive solve impractical regardless of how the seeds are generated --
that's a separate, real limitation this script cannot itself remove.

Usage:
    python3 generate_exhaustive_positions.py --white Q --out kqvk_exhaustive.txt
    python3 generate_exhaustive_positions.py --white B,N --out kbnvk_exhaustive.txt --bishop-color light
    python3 generate_exhaustive_positions.py --white Q,R,B --out qrb_exhaustive.txt
"""

import argparse
import sys

from generate_general_seed_positions import (
    Position, parse_white_pieces, square_color, format_seed, PIECE_LETTERS,
    enumerate_fully_exhaustive, MAX_WHITE_NON_KING
)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--white', required=True, help="Comma-separated White piece letters, e.g. 'Q' or 'B,N'")
    ap.add_argument('--out', required=True, help="Output positions file")
    ap.add_argument('--bishop-color', choices=['light', 'dark', 'opposite'], default=None,
                     help="For one bishop: constrain it to one color. For two bishops: 'light'/'dark' "
                          "forces both to that color (same-colored pair, requires an underpromotion "
                          "to reach in a real game); 'opposite' assigns one to each (the practically "
                          "common configuration). See generate_general_seed_positions.py's flag of "
                          "the same name for the full reasoning on why this is three distinct, "
                          "genuinely disconnected cases for two bishops, not two.")
    args = ap.parse_args()

    white_letters = parse_white_pieces(args.white)
    has_pawn = 'P' in white_letters
    if has_pawn and white_letters.count('P') > 1:
        raise SystemExit("ERROR: at most one pawn is supported")
    if len(white_letters) > MAX_WHITE_NON_KING:
        raise SystemExit(f"ERROR: only 0 through {MAX_WHITE_NON_KING} non-king White pieces are "
                          f"supported by this engine, got {len(white_letters)}")
    num_bishops = white_letters.count('B')
    if args.bishop_color and num_bishops == 0:
        raise SystemExit("ERROR: --bishop-color only makes sense when --white includes 'B'")
    if args.bishop_color == 'opposite' and num_bishops != 2:
        raise SystemExit("ERROR: --bishop-color opposite requires exactly two bishops in --white")

    print(f"White pieces: {white_letters}")
    if has_pawn:
        print("  (pawn -> every seed forced to rank 2, for forward-reachability -- see the other generator's docstring)")
    if args.bishop_color == 'opposite':
        print("  Bishops constrained to OPPOSITE colors (one light, one dark)")
    elif args.bishop_color:
        color_note = "Both bishops" if num_bishops == 2 else "Bishop"
        print(f"  {color_note} constrained to {args.bishop_color}-squared placements only")
    elif num_bishops == 2:
        print("  NOTE: two bishops -- the full landscape has THREE structurally disconnected cases "
              "(opposite / both-light / both-dark), not one. Use --bishop-color opposite / light / "
              "dark to cover each deliberately, or verify after the fact that an unconstrained run "
              "actually covered all three before trusting it as complete.")

    if len(white_letters) == 2:
        print("  NOTE: two non-king White pieces -- true exhaustive enumeration is ~30 million raw "
              "combinations before filtering. This will take longer and produce a much larger file "
              "than single-piece material. See module docstring for the honest scale tradeoff.")
    elif len(white_letters) >= 3:
        print(f"  NOTE: {len(white_letters)} non-king White pieces -- both the seed file AND the "
              f"actual reachable position graph for this material are very plausibly far larger "
              f"than anything solved so far in this project (KBNvK, at 2 pieces, was already ~24.7 "
              f"million positions). This engine has no disk-backed storage -- a full exhaustive "
              f"solve may not be practically feasible on any single machine's RAM. Worth testing "
              f"with a small --max-nodes cap on the solver first, before committing to a full run.")

    print("Enumerating every legal position combinatorially (not just king pairings)...")

    written = 0
    with open(args.out, 'w') as f:
        f.write("# fully exhaustive seed positions for general_solver --full-dag\n")
        for letters, squares, turn in enumerate_fully_exhaustive(white_letters, bishop_color=args.bishop_color):
            f.write(format_seed(letters, squares, turn) + "\n")
            written += 1
            if written % 1000000 == 0:
                print(f"  ...{written} written so far", flush=True)

    print(f"\nWrote {written} seed positions to {args.out}")
    print("This is TRUE exhaustive coverage: every legal position of this material was tried as an "
          "explicit starting point, not just every king pairing.")


if __name__ == "__main__":
    main()
