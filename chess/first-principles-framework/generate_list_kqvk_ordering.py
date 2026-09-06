#!/usr/bin/env python3
"""
Generate all legal KQvK positions and rank by local Syzygy DTZ files.
Reads from local tablebase directory instead of querying online.
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

class KQvKPositionGenerator:
    def generate_all_king_moves(self, pos):
        """Mirror C++ generate_all_king_moves"""
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
    
    def generate_all_queen_moves(self, pos):
        """Mirror C++ generate_all_queen_moves"""
        moves = []
        dirs = [(1,0), (-1,0), (0,1), (0,-1), (1,1), (1,-1), (-1,1), (-1,-1)]
        
        for d in dirs:
            for dist in range(1, 8):
                nf = pos.file + d[0] * dist
                nr = pos.rank + d[1] * dist
                if 0 <= nf <= 7 and 0 <= nr <= 7:
                    moves.append(Position(nf, nr))
                else:
                    break
        return moves
    
    def is_attacked_by_queen(self, pos, qp, wk, wq):
        """Mirror C++ is_attacked_by_queen with blocking"""
        if pos == qp:
            return False
        
        # Check file
        if pos.file == qp.file:
            start = min(pos.rank, qp.rank) + 1
            end = max(pos.rank, qp.rank)
            for r in range(start, end):
                if Position(pos.file, r) == wk or Position(pos.file, r) == wq:
                    return False
            return True
        
        # Check rank
        if pos.rank == qp.rank:
            start = min(pos.file, qp.file) + 1
            end = max(pos.file, qp.file)
            for f in range(start, end):
                if Position(f, pos.rank) == wk or Position(f, pos.rank) == wq:
                    return False
            return True
        
        # Check diagonals
        if abs(pos.file - qp.file) == abs(pos.rank - qp.rank):
            df = 1 if pos.file > qp.file else -1
            dr = 1 if pos.rank > qp.rank else -1
            f = qp.file + df
            r = qp.rank + dr
            while f != pos.file:
                if Position(f, r) == wk or Position(f, r) == wq:
                    return False
                f += df
                r += dr
            return True
        
        return False
    
    def is_legal_state(self, wk, wq, bk):
        """Mirror C++ is_legal_state"""
        positions = {(wk.file, wk.rank), (wq.file, wq.rank), (bk.file, bk.rank)}
        if len(positions) != 3:
            return False
        if wk.distance_to(bk) < 2:
            return False
        return True
    
    def generate_positions(self):
        """Generate all legal KQvK positions (White to move only)"""
        positions = []
        
        for wkf in range(8):
            for wkr in range(8):
                for wqf in range(8):
                    for wqr in range(8):
                        for bkf in range(8):
                            for bkr in range(8):
                                wk = Position(wkf, wkr)
                                wq = Position(wqf, wqr)
                                bk = Position(bkf, bkr)
                                
                                if wk == wq or wk == bk or wq == bk:
                                    continue
                                
                                if wk.distance_to(bk) < 2:
                                    continue
                                
                                if not self.is_legal_state(wk, wq, bk):
                                    continue
                                
                                if self.is_attacked_by_queen(bk, wq, wk, wq):
                                    continue
                                
                                positions.append((wk, wq, bk))
        
        return positions


def query_syzygy_local(wk, wq, bk, tb):
    """Query local Syzygy tablebase"""
    try:
        board = chess.Board()
        board.clear()
        
        wk_square = chess.square(wk.file, wk.rank)
        wq_square = chess.square(wq.file, wq.rank)
        bk_square = chess.square(bk.file, bk.rank)
        
        board.set_piece_at(wk_square, chess.Piece(chess.KING, chess.WHITE))
        board.set_piece_at(wq_square, chess.Piece(chess.QUEEN, chess.WHITE))
        board.set_piece_at(bk_square, chess.Piece(chess.KING, chess.BLACK))
        
        board.turn = chess.WHITE
        
        # Query local tablebase
        dtz = tb.probe_dtz(board)
        
        if dtz is not None:
            return {
                'wk': wk.to_str(),
                'wq': wq.to_str(),
                'bk': bk.to_str(),
                'position_str': f"WK:{wk.to_str()} WQ:{wq.to_str()} BK:{bk.to_str()}",
                'dtz': abs(dtz)
            }
        return None
    
    except Exception as e:
        return None


def main():
    print("=" * 80)
    print("KQvK Position Generator and Syzygy Ranker (Local Tablebase)")
    print("=" * 80)
    
    # Find Syzygy tablebase directory
    syzygy_dir = Path("syzygy")
    if not syzygy_dir.exists():
        print(f"ERROR: Syzygy directory not found at {syzygy_dir}")
        print("Please create a 'syzygy' directory and place kqvk.dtz inside")
        return
    
    print(f"\nUsing Syzygy tablebase from: {syzygy_dir.absolute()}")
    
    # Open tablebase
    print("Opening Syzygy tablebase...")
    try:
        tb = chess.syzygy.open_tablebase(str(syzygy_dir))
    except Exception as e:
        print(f"ERROR: Failed to open tablebase: {e}")
        return
    
    # Generate positions
    print("Generating all legal KQvK positions...")
    gen = KQvKPositionGenerator()
    positions = gen.generate_positions()
    print(f"Generated {len(positions)} legal positions\n")
    
    # Query each position
    print("Querying local Syzygy tablebase...\n")
    
    position_dtz = []
    failed_count = 0
    
    with tqdm(total=len(positions), desc="Querying", unit="pos",
              bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]') as pbar:
        
        for wk, wq, bk in positions:
            result = query_syzygy_local(wk, wq, bk, tb)
            
            if result is not None:
                position_dtz.append(result)
            else:
                failed_count += 1
            
            pbar.update(1)
    
    print(f"\n✓ Successfully retrieved: {len(position_dtz)}")
    print(f"✗ Failed: {failed_count}\n")
    
    if len(position_dtz) == 0:
        print("ERROR: No positions retrieved from Syzygy!")
        return
    
    # Sort by DTZ
    print("Sorting by DTZ (moves to mate)...")
    position_dtz.sort(key=lambda x: x['dtz'])
    
    # Write output
    output_file = "kqvk_positions_by_dtz.txt"
    print(f"Writing to {output_file}...")
    
    with open(output_file, 'w') as f:
        f.write("DTZ,Position\n")
        for pos in position_dtz:
            f.write(f"{pos['dtz']},{pos['position_str']}\n")
    
    # Print statistics
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