#!/usr/bin/env python3
"""
Validate KQvK landscape against Syzygy DTZ.
Compare moves-to-mate values and show FULL ENDGAME from landscape.
"""

import csv
from pathlib import Path
from collections import defaultdict

class ValidationAnalyzer:
    def __init__(self, dtz_file, landscape_file):
        self.dtz_data = {}
        self.landscape_data = {}
        self.mismatches = []
        self.matches = []
        
        self.load_dtz(dtz_file)
        self.load_landscape(landscape_file)
    
    def load_dtz(self, filename):
        """Load Syzygy DTZ data: DTZ,Position (DTZ positions are always White to move).
        NOTE: DTZ is expressed in plies (half-moves), matching the landscape's
        total_plies field -- NOT the landscape's M field, which is White's move
        count only."""
        print(f"Loading Syzygy DTZ from {filename}...")
        
        try:
            with open(filename, 'r') as f:
                reader = csv.reader(f)
                header = next(reader)
                
                for row in reader:
                    if len(row) < 2:
                        continue
                    
                    dtz_str = row[0].strip()
                    position = row[1].strip()
                    
                    try:
                        dtz = int(dtz_str)
                        self.dtz_data[(position, 'W')] = dtz
                    except ValueError:
                        pass
        except Exception as e:
            print(f"ERROR loading DTZ: {e}")
            return
        
        print(f"  Loaded {len(self.dtz_data)} positions from Syzygy DTZ\n")
    
    def load_landscape(self, filename):
        """Load landscape database with full data, keyed by (position, turn)"""
        print(f"Loading landscape database from {filename}...")
        
        try:
            with open(filename, 'r') as f:
                reader = csv.reader(f, delimiter='|')
                header = next(reader)
                
                for row in reader:
                    if len(row) < 4:
                        continue
                    
                    position = row[0].strip()
                    turn = row[1].strip() if len(row) > 1 else "?"
                    best_move = row[2].strip() if len(row) > 2 else "?"
                    m_value = int(row[3]) if len(row) > 3 else 0
                    total_plies = int(row[4]) if len(row) > 4 else 0
                    white_moves = int(row[5]) if len(row) > 5 else 0
                    black_moves = int(row[6]) if len(row) > 6 else 0
                    nodes_eval = int(row[7]) if len(row) > 7 else 0
                    comp_time = float(row[8]) if len(row) > 8 else 0.0
                    bn_trajectory = row[9].strip() if len(row) > 9 else ""
                    
                    self.landscape_data[(position, turn)] = {
                        'turn': turn,
                        'best_move': best_move,
                        'M': m_value,
                        'total_plies': total_plies,
                        'white_moves': white_moves,
                        'black_moves': black_moves,
                        'nodes_eval': nodes_eval,
                        'comp_time': comp_time,
                        'bn_trajectory': bn_trajectory
                    }
        except Exception as e:
            print(f"ERROR loading landscape: {e}")
            return
        
        print(f"  Loaded {len(self.landscape_data)} positions from landscape\n")
    
    def parse_position(self, pos_str):
        """Parse position string 'WK:a1 WQ:b2 BK:c3' into components"""
        parts = {}
        for part in pos_str.split():
            if ':' in part:
                key, val = part.split(':')
                parts[key] = val
        return parts
    
    def apply_move(self, pos_str, move_str):
        """Apply a move to a position and return the resulting position string"""
        parts = self.parse_position(pos_str)
        
        if not move_str or move_str == "mate":
            return None
        
        # Extract piece and destination from move (e.g., "Ka1" or "Qa5" or "kb2")
        piece = move_str[0]  # K, Q, or k
        dest = move_str[1:]  # destination square
        
        if piece == 'K':
            parts['WK'] = dest
        elif piece == 'Q':
            parts['WQ'] = dest
        elif piece == 'k':
            parts['BK'] = dest
        else:
            return None
        
        # Reconstruct position string
        result = f"WK:{parts.get('WK', '?')} WQ:{parts.get('WQ', '?')} BK:{parts.get('BK', '?')}"
        return result
    
    def reconstruct_full_game(self, start_position, start_turn='W'):
        """Reconstruct the full endgame line from a starting position/turn"""
        moves = []
        positions = [start_position]
        current_pos = start_position
        current_turn = start_turn
        
        for _ in range(100):  # Max 100 plies to prevent infinite loops
            key = (current_pos, current_turn)
            if key not in self.landscape_data:
                break
            
            data = self.landscape_data[key]
            best_move = data['best_move']
            
            if not best_move or best_move == "mate":
                moves.append(best_move)
                break
            
            moves.append(best_move)
            
            # Apply move to get next position
            next_pos = self.apply_move(current_pos, best_move)
            if next_pos is None:
                break
            
            current_turn = 'B' if best_move[0] in ('K', 'Q') else 'W'
            positions.append(next_pos)
            current_pos = next_pos
        
        return {
            'moves': moves,
            'positions': positions,
            'game_line': ' '.join(moves)
        }
    
    def validate(self):
        """Compare DTZ vs total_plies. DTZ is measured in plies, so it must be
        compared against the landscape's total_plies field, not M (which is
        White's move count only -- a different unit that will never line up
        with DTZ directly)."""
        print("Validating positions...\n")
        
        common_keys = set(self.dtz_data.keys()) & set(self.landscape_data.keys())
        print(f"Positions in both datasets: {len(common_keys)}")
        
        dtz_only = set(self.dtz_data.keys()) - set(self.landscape_data.keys())
        landscape_only = set(self.landscape_data.keys()) - set(self.dtz_data.keys())
        
        print(f"Positions only in Syzygy: {len(dtz_only)}")
        print(f"Positions only in landscape: {len(landscape_only)}")
        print(f"  (landscape entries with turn='B' will always fall here -- "
              f"Syzygy DTZ only covers White-to-move positions)\n")
        
        for key in common_keys:
            position, turn = key
            dtz = self.dtz_data[key]
            landscape = self.landscape_data[key]
            plies = landscape['total_plies']
            
            if dtz == plies:
                self.matches.append({
                    'position': position,
                    'turn': turn,
                    'DTZ': dtz,
                    'landscape': landscape
                })
            else:
                self.mismatches.append({
                    'position': position,
                    'turn': turn,
                    'DTZ': dtz,
                    'landscape': landscape,
                    'delta': plies - dtz
                })
        
        print(f"Matches: {len(self.matches)}")
        print(f"Mismatches: {len(self.mismatches)}\n")
    
    def analyze_mismatches(self):
        """Analyze mismatch patterns with full game lines"""
        if not self.mismatches:
            print("No mismatches found! ✓\n")
            return
        
        print("MISMATCH ANALYSIS WITH FULL ENDGAME LINES")
        print("=" * 120)
        
        deltas = defaultdict(list)
        for mismatch in self.mismatches:
            delta = mismatch['delta']
            deltas[delta].append(mismatch)
        
        print("\nMismatches by delta (total_plies - DTZ):")
        for delta in sorted(deltas.keys()):
            count = len(deltas[delta])
            pct = 100 * count / len(self.mismatches)
            print(f"  Delta {delta:+3d}: {count:6d} positions ({pct:5.1f}%)")
        
        all_deltas = [m['delta'] for m in self.mismatches]
        avg_delta = sum(all_deltas) / len(all_deltas)
        
        print(f"\nDelta statistics:")
        print(f"  Average delta: {avg_delta:.2f}")
        print(f"  Min delta: {min(all_deltas)}")
        print(f"  Max delta: {max(all_deltas)}")
        
        # Critical mismatches -- landscape claims mate FASTER than Syzygy proves possible
        critical = [m for m in self.mismatches if m['landscape']['total_plies'] < m['DTZ']]
        if critical:
            print(f"\n⚠ CRITICAL MISMATCHES (plies < DTZ): {len(critical)}")
            print("=" * 120)
            
            for mismatch in critical[:5]:
                self._print_mismatch_detail(mismatch)
        
        # Positive mismatches -- landscape takes longer than true optimal
        positive = [m for m in self.mismatches if m['landscape']['total_plies'] > m['DTZ']]
        if positive:
            print(f"\n⚠ SUBOPTIMAL MISMATCHES (plies > DTZ): {len(positive)}")
            print("=" * 120)
            
            for mismatch in positive[:5]:
                self._print_mismatch_detail(mismatch)
    
    def _print_mismatch_detail(self, mismatch):
        """Print detailed mismatch with FULL ENDGAME LINE"""
        pos = mismatch['position']
        dtz = mismatch['DTZ']
        landscape = mismatch['landscape']
        plies = landscape['total_plies']
        
        print(f"\nPosition: {pos}")
        print(f"  Syzygy DTZ: {dtz}")
        print(f"  Landscape total_plies: {plies} (delta: {plies - dtz:+d})")
        print(f"  Landscape M (White moves only, informational -- not comparable to DTZ): {landscape['M']}")
        
        # Reconstruct full game
        game_data = self.reconstruct_full_game(pos, mismatch['turn'])
        
        print(f"\n  LANDSCAPE FULL ENDGAME:")
        print(f"    Move sequence: {game_data['game_line']}")
        print(f"    Number of moves: {len([m for m in game_data['moves'] if m != 'mate'])}")
        print(f"    Positions in sequence: {len(game_data['positions'])}")
        
        print(f"\n  Landscape Solution Details:")
        print(f"    Best move from start: {landscape['best_move']}")
        print(f"    Total plies: {landscape['total_plies']}")
        print(f"    White moves: {landscape['white_moves']}")
        print(f"    Black moves: {landscape['black_moves']}")
        print(f"    Nodes evaluated: {landscape['nodes_eval']}")
        print(f"    Computation time: {landscape['comp_time']:.3f}s")
        
        if landscape['bn_trajectory']:
            bn_parts = landscape['bn_trajectory'].split(',')
            print(f"\n    Black node trajectory ({len(bn_parts)} moves):")
            print(f"      {landscape['bn_trajectory']}")
            
            try:
                bn_vals = [int(x.strip()) for x in bn_parts if x.strip()]
                if bn_vals:
                    print(f"      Min BN: {min(bn_vals)}, Max BN: {max(bn_vals)}, Avg BN: {sum(bn_vals)/len(bn_vals):.1f}")
            except:
                pass
    
    def write_mismatches(self, output_file):
        """Write mismatches with FULL ENDGAME LINES"""
        print(f"\nWriting mismatches to {output_file}...")
        
        with open(output_file, 'w') as f:
            f.write("Position,DTZ,Plies,Delta,Status,M_WhiteMovesOnly,BestMove,FullGameLine,NumMoves,WhiteMoves,BlackMoves,NodesEval,CompTime,BN_Trajectory\n")
            
            for mismatch in sorted(self.mismatches, key=lambda x: abs(x['delta']), reverse=True):
                pos = mismatch['position']
                dtz = mismatch['DTZ']
                landscape = mismatch['landscape']
                plies = landscape['total_plies']
                delta = plies - dtz
                
                if plies < dtz:
                    status = "CRITICAL_SHORTER"
                elif plies > dtz:
                    status = "SUBOPTIMAL_LONGER"
                else:
                    status = "MATCH"
                
                # Reconstruct full game
                game_data = self.reconstruct_full_game(pos, mismatch['turn'])
                game_line = game_data['game_line']
                num_moves = len([m for m in game_data['moves'] if m != 'mate'])
                
                best_move = landscape['best_move']
                m_value = landscape['M']
                white_moves = landscape['white_moves']
                black_moves = landscape['black_moves']
                nodes_eval = landscape['nodes_eval']
                comp_time = landscape['comp_time']
                bn_traj = landscape['bn_trajectory']
                
                f.write(f'"{pos}",{dtz},{plies},{delta:+d},{status},{m_value},"{best_move}","{game_line}",{num_moves},{white_moves},{black_moves},{nodes_eval},{comp_time:.3f},"{bn_traj}"\n')
        
        print(f"✓ Wrote {len(self.mismatches)} mismatches to {output_file}\n")
    
    def write_matches(self, output_file):
        """Write matches with FULL ENDGAME LINES"""
        print(f"Writing matches to {output_file}...")
        
        with open(output_file, 'w') as f:
            f.write("Position,DTZ,Plies,M_WhiteMovesOnly,BestMove,FullGameLine,NumMoves,WhiteMoves,BlackMoves,NodesEval,CompTime,BN_Trajectory\n")
            
            for match in self.matches:
                pos = match['position']
                dtz = match['DTZ']
                landscape = match['landscape']
                plies = landscape['total_plies']
                
                # Reconstruct full game
                game_data = self.reconstruct_full_game(pos, match['turn'])
                game_line = game_data['game_line']
                num_moves = len([m for m in game_data['moves'] if m != 'mate'])
                
                best_move = landscape['best_move']
                m_value = landscape['M']
                white_moves = landscape['white_moves']
                black_moves = landscape['black_moves']
                nodes_eval = landscape['nodes_eval']
                comp_time = landscape['comp_time']
                bn_traj = landscape['bn_trajectory']
                
                f.write(f'"{pos}",{dtz},{plies},{m_value},"{best_move}","{game_line}",{num_moves},{white_moves},{black_moves},{nodes_eval},{comp_time:.3f},"{bn_traj}"\n')
        
        print(f"✓ Wrote {len(self.matches)} matches to {output_file}\n")


def main():
    print("=" * 120)
    print("KQvK Landscape Validation Against Syzygy (with Full Endgame Lines)")
    print("=" * 120 + "\n")
    
    dtz_file = "kqvk_positions_by_dtz.txt"
    landscape_file = "kqvk_perfect_play.db"
    
    if not Path(dtz_file).exists():
        print(f"ERROR: {dtz_file} not found")
        return
    
    if not Path(landscape_file).exists():
        print(f"ERROR: {landscape_file} not found")
        return
    
    validator = ValidationAnalyzer(dtz_file, landscape_file)
    validator.validate()
    validator.analyze_mismatches()
    
    validator.write_mismatches("mismatches_full_games.csv")
    validator.write_matches("matches_full_games.csv")
    
    print("=" * 120)
    print("VALIDATION COMPLETE")
    print("=" * 120)
    
    total = len(validator.matches) + len(validator.mismatches)
    if total > 0:
        accuracy = 100 * len(validator.matches) / total
        print(f"\nAccuracy: {accuracy:.2f}%")
        print(f"Matches: {len(validator.matches)}/{total}")
        print(f"Mismatches: {len(validator.mismatches)}/{total}\n")

if __name__ == "__main__":
    main()
