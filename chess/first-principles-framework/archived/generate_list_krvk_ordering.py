#!/usr/bin/env python3
"""
Generate all legal KRvK positions and rank by local Syzygy DTZ files.
Direct adaptation of the KQvK generator: same structure and same output
format (for compatibility with the C++ solver's --full-dag position loader),
rook geometry instead of queen geometry.
"""

import chess
import chess.syzygy
from pathlib import Path
from collections import defaultdict
import time
from tqdm import tqdm

class Position:
    def __init__(self, file, rank):
        self.file = file
        self.rank = rank

    def __eq__(self, other):
        return self.file == other.file and self.rank == other.rank

    def __hash__(self):
        return hash((self.file, self.rank))

    def distance_to(self, other):
        return max(abs(self.file - other.file), abs(self.rank - other.rank))

    def to_str(self):
        return chr(ord('a') + self.file) + str(self.rank + 1)

    @staticmethod
    def from_str(s):
        return Position(ord(s[0]) - ord('a'), int(s[1]) - 1)


class KRvKPositionGenerator:
    def generate_all_king_moves(self, pos):
        """Mirror C++ generate_all_king_moves -- unchanged, piece-agnostic already."""
        moves = []
        for df in [-1, 0, 1]:
            for dr in [-1, 0, 1]:
                if df == 0 and dr == 0:
                    continue
                nf = pos.file + df
                nr = pos.rank + dr
                if 0 <= nf <= 7 and 0 <= nr <= 7:
                    moves.append(Position(nf, nr))
        return moves

    def generate_all_rook_moves(self, pos):
        """Mirror C++ RookRules::generate_moves -- 4 orthogonal directions only,
        no diagonals (the only geometric difference from the queen version)."""
        moves = []
        dirs = [(1, 0), (-1, 0), (0, 1), (0, -1)]

        for d in dirs:
            for dist in range(1, 8):
                nf = pos.file + d[0] * dist
                nr = pos.rank + d[1] * dist
                if 0 <= nf <= 7 and 0 <= nr <= 7:
                    moves.append(Position(nf, nr))
                else:
                    break
        return moves

    def is_attacked_by_rook(self, pos, rp, wk, wr):
        """Mirror C++ RookRules::attacks -- file/rank blocking checks only,
        no diagonal branch (a rook simply doesn't threaten off-file/off-rank
        squares, unlike the queen version this was adapted from)."""
        if pos == rp:
            return False

        # Check file
        if pos.file == rp.file:
            start = min(pos.rank, rp.rank) + 1
            end = max(pos.rank, rp.rank)
            for r in range(start, end):
                if Position(pos.file, r) == wk or Position(pos.file, r) == wr:
                    return False
            return True

        # Check rank
        if pos.rank == rp.rank:
            start = min(pos.file, rp.file) + 1
            end = max(pos.file, rp.file)
            for f in range(start, end):
                if Position(f, pos.rank) == wk or Position(f, pos.rank) == wr:
                    return False
            return True

        # No diagonal branch here -- this is the entire geometric difference
        # from the KQvK version's is_attacked_by_queen.
        return False

    def is_legal_state(self, wk, wr, bk):
        """Mirror C++ is_legal_state -- unchanged, already piece-agnostic."""
        positions = {(wk.file, wk.rank), (wr.file, wr.rank), (bk.file, bk.rank)}
        if len(positions) != 3:
            return False
        if wk.distance_to(bk) < 2:
            return False
        return True

    def generate_positions(self):
        """Generate all legal KRvK positions (White to move only)."""
        positions = []

        for wkf in range(8):
            for wkr in range(8):
                for wrf in range(8):
                    for wrr in range(8):
                        for bkf in range(8):
                            for bkr in range(8):
                                wk = Position(wkf, wkr)
                                wr = Position(wrf, wrr)
                                bk = Position(bkf, bkr)

                                if wk == wr or wk == bk or wr == bk:
                                    continue

                                if wk.distance_to(bk) < 2:
                                    continue

                                if not self.is_legal_state(wk, wr, bk):
                                    continue

                                if self.is_attacked_by_rook(bk, wr, wk, wr):
                                    continue

                                positions.append((wk, wr, bk))

        return positions


def query_syzygy_local(wk, wr, bk, tb):
    """Query local Syzygy KRK tablebase."""
    try:
        board = chess.Board()
        board.clear()

        wk_square = chess.square(wk.file, wk.rank)
        wr_square = chess.square(wr.file, wr.rank)
        bk_square = chess.square(bk.file, bk.rank)

        board.set_piece_at(wk_square, chess.Piece(chess.KING, chess.WHITE))
        board.set_piece_at(wr_square, chess.Piece(chess.ROOK, chess.WHITE))
        board.set_piece_at(bk_square, chess.Piece(chess.KING, chess.BLACK))

        board.turn = chess.WHITE

        dtz = tb.probe_dtz(board)

        if dtz is not None:
            return {
                'wk': wk.to_str(),
                'wr': wr.to_str(),
                'bk': bk.to_str(),
                # Labeled "WQ:" -- NOT "WR:" -- deliberately. The C++ solver's
                # internal position-string format always uses "WQ:" for White's
                # second piece regardless of which piece it actually is (a
                # pragmatic choice to avoid renaming that field across the whole
                # engine). load_positions_from_file() hardcodes searching for the
                # literal substring "WQ:", so this file must use that label to be
                # read correctly by the --full-dag KRvK binary -- it's holding a
                # rook's square here, not a queen's.
                'position_str': f"WK:{wk.to_str()} WQ:{wr.to_str()} BK:{bk.to_str()}",
                'dtz': abs(dtz)
            }
        return None

    except Exception:
        return None


def main():
    print("=" * 80)
    print("KRvK Position Generator and Syzygy Ranker (Local Tablebase)")
    print("=" * 80)

    syzygy_dir = Path("syzygy")
    if not syzygy_dir.exists():
        print(f"ERROR: Syzygy directory not found at {syzygy_dir}")
        print("Please create a 'syzygy' directory and place KRK.rtbw/KRK.rtbz inside")
        return

    print(f"\nUsing Syzygy tablebase from: {syzygy_dir.absolute()}")
    print("(needs KRK.rtbw / KRK.rtbz specifically -- if your syzygy/ directory")
    print(" already holds the full 3-4-5-piece Syzygy set, these ship alongside")
    print(" KQK.rtbw/KQK.rtbz and you may already have them)")

    print("\nOpening Syzygy tablebase...")
    try:
        tb = chess.syzygy.open_tablebase(str(syzygy_dir))
    except Exception as e:
        print(f"ERROR: Failed to open tablebase: {e}")
        return

    print("Generating all legal KRvK positions...")
    gen = KRvKPositionGenerator()
    positions = gen.generate_positions()
    print(f"Generated {len(positions)} legal positions\n")

    print("Querying local Syzygy tablebase...\n")

    position_dtz = []
    failed_count = 0

    with tqdm(total=len(positions), desc="Querying", unit="pos",
              bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]') as pbar:

        for wk, wr, bk in positions:
            result = query_syzygy_local(wk, wr, bk, tb)

            if result is not None:
                position_dtz.append(result)
            else:
                failed_count += 1

            pbar.update(1)

    print(f"\n✓ Successfully retrieved: {len(position_dtz)}")
    print(f"✗ Failed: {failed_count}\n")

    if len(position_dtz) == 0:
        print("ERROR: No positions retrieved from Syzygy!")
        print("Most likely cause: KRK.rtbw/KRK.rtbz are missing from the syzygy/ directory.")
        return

    print("Sorting by DTZ (moves to mate)...")
    position_dtz.sort(key=lambda x: x['dtz'])

    output_file = "krvk_positions_by_dtz.txt"
    print(f"Writing to {output_file}...")

    with open(output_file, 'w') as f:
        f.write("DTZ,Position\n")
        for pos in position_dtz:
            f.write(f"{pos['dtz']},{pos['position_str']}\n")

    print("\n" + "=" * 80)
    print("STATISTICS")
    print("=" * 80)
    print(f"Total positions: {len(position_dtz)}")

    dtz_groups = defaultdict(int)
    for pos in position_dtz:
        dtz_groups[pos['dtz']] += 1

    print("\nPositions by DTZ:")
    for dtz in sorted(dtz_groups.keys()):
        count = dtz_groups[dtz]
        print(f"  DTZ {dtz:2d}: {count:6d} positions")

    dtz_values = [p['dtz'] for p in position_dtz]
    print(f"\nLowest DTZ: {min(dtz_values)}")
    print(f"Highest DTZ: {max(dtz_values)}")
    print(f"Average DTZ: {sum(dtz_values) / len(dtz_values):.1f}")
    print(f"\nOutput file: {output_file}")
    print("=" * 80)


if __name__ == "__main__":
    main()
