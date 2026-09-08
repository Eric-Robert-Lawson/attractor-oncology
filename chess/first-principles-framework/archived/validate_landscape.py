#!/usr/bin/env python3
"""
Validate a generated landscape database against Syzygy DTZ ground truth.
Piece-agnostic: works for KQvK, KRvK, or any future endgame that shares the
same file formats -- pass the right --syzygy and --landscape files and it
just works, no code changes needed per endgame.

Compares moves-to-mate values and shows the full endgame line for every
mismatch, so a discrepancy is immediately followed by the exact game that
produced it rather than just a bare number.

Usage:
    python3 validate_landscape.py \\
        --syzygy kqvk_positions_by_dtz.txt --landscape kqvk_perfect_play.db \\
        --out-dir validation_kqvk

    python3 validate_landscape.py \\
        --syzygy krvk_positions_by_dtz.txt --landscape krvk_perfect_play.db \\
        --out-dir validation_krvk

--out-dir defaults to "validation" but should be given explicitly when
validating more than one endgame from the same working directory, so one
run's mismatches.csv/matches.csv don't overwrite another's.
"""

import argparse
import csv
import os
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
        duplicate_conflicts = []

        try:
            with open(filename, 'r') as f:
                reader = csv.reader(f)
                header = next(reader)

                for row_num, row in enumerate(reader, start=2):
                    if len(row) < 2:
                        continue

                    dtz_str = row[0].strip()
                    position = row[1].strip()

                    try:
                        dtz = int(dtz_str)
                        key = (position, 'W')
                        if key in self.dtz_data and self.dtz_data[key] != dtz:
                            duplicate_conflicts.append((row_num, position, self.dtz_data[key], dtz))
                        self.dtz_data[key] = dtz
                    except ValueError:
                        pass
        except Exception as e:
            print(f"ERROR loading DTZ: {e}")
            return

        print(f"  Loaded {len(self.dtz_data)} positions from Syzygy DTZ\n")

        if duplicate_conflicts:
            print(f"  \u26a0 WARNING: {len(duplicate_conflicts)} positions appeared more than once "
                  f"in {filename} with DIFFERENT DTZ values. Whichever row came last silently won "
                  f"-- this file should not be trusted as ground truth until these are resolved:")
            for row_num, pos, old_dtz, new_dtz in duplicate_conflicts[:20]:
                print(f"    line {row_num}: {pos}  previous_dtz={old_dtz}  new_dtz={new_dtz}")
            if len(duplicate_conflicts) > 20:
                print(f"    ... and {len(duplicate_conflicts) - 20} more")
            print()

    def load_landscape(self, filename):
        """Load landscape database with full data, keyed by (position, turn)"""
        print(f"Loading landscape database from {filename}...")
        duplicate_conflicts = []

        try:
            with open(filename, 'r') as f:
                reader = csv.reader(f, delimiter='|')
                header = next(reader)

                for row_num, row in enumerate(reader, start=2):
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

                    key = (position, turn)
                    if key in self.landscape_data and self.landscape_data[key]['total_plies'] != total_plies:
                        duplicate_conflicts.append(
                            (row_num, position, turn, self.landscape_data[key]['total_plies'], total_plies))

                    self.landscape_data[key] = {
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

        if duplicate_conflicts:
            print(f"  \u26a0 WARNING: {len(duplicate_conflicts)} (position, turn) pairs appeared more "
                  f"than once in {filename} with DIFFERENT total_plies values. Whichever row came "
                  f"last silently won -- this normally shouldn't happen given add_position's "
                  f"first-write-wins guard, so this points at something worth investigating:")
            for row_num, pos, turn, old_plies, new_plies in duplicate_conflicts[:20]:
                print(f"    line {row_num}: {pos} turn={turn}  previous_plies={old_plies}  new_plies={new_plies}")
            if len(duplicate_conflicts) > 20:
                print(f"    ... and {len(duplicate_conflicts) - 20} more")
            print()

    def parse_position(self, pos_str):
        """Parse position string 'WK:a1 WQ:b2 BK:c3' into components. The 'WQ'
        key holds White's second piece regardless of what that piece actually
        is (queen, rook, ...) -- see apply_move below."""
        parts = {}
        for part in pos_str.split():
            if ':' in part:
                key, val = part.split(':')
                parts[key] = val
        return parts

    def apply_move(self, pos_str, move_str):
        """Apply a move to a position and return the resulting position string.
        Piece-agnostic via .isupper() rather than hardcoded to 'K'/'Q'/'k': a
        KRvK move looks like "Rb8", not "Qb8", so a hardcoded 'Q' check would
        silently fail to recognize every rook move and leave the position
        unchanged -- the exact bug class already found and fixed in the C++
        solver and in kqvk_analyzer.py/krvk_analyzer.py. Checking "is this
        letter uppercase" instead correctly recognizes ANY white second-piece
        letter, so this carries over unchanged to a bishop, knight, etc. too.
        """
        parts = self.parse_position(pos_str)

        if not move_str or move_str == "mate":
            return None

        piece = move_str[0]  # 'K', some uppercase White second-piece letter, or 'k'
        dest = move_str[1:]

        if piece == 'K':
            parts['WK'] = dest
        elif piece.isupper():
            parts['WQ'] = dest  # position-string label stays "WQ" regardless of piece
        elif piece == 'k':
            parts['BK'] = dest
        else:
            return None

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

            next_pos = self.apply_move(current_pos, best_move)
            if next_pos is None:
                break

            # Piece-agnostic, same reasoning as apply_move above: whichever
            # uppercase letter moved (K or White's second piece), it was
            # White's move, so the turn flips to Black; a lowercase letter
            # ('k') means Black moved, flipping back to White.
            current_turn = 'B' if best_move[0].isupper() else 'W'
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
            print("No mismatches found! \u2713\n")
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
            print(f"\n\u26a0 CRITICAL MISMATCHES (plies < DTZ): {len(critical)}")
            print("=" * 120)

            for mismatch in critical[:5]:
                self._print_mismatch_detail(mismatch)

        # Positive mismatches -- landscape takes longer than true optimal
        positive = [m for m in self.mismatches if m['landscape']['total_plies'] > m['DTZ']]
        if positive:
            print(f"\n\u26a0 SUBOPTIMAL MISMATCHES (plies > DTZ): {len(positive)}")
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
            except Exception:
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

        print(f"\u2713 Wrote {len(self.mismatches)} mismatches to {output_file}\n")

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

        print(f"\u2713 Wrote {len(self.matches)} matches to {output_file}\n")


def main():
    ap = argparse.ArgumentParser(
        description="Validate a generated landscape database against Syzygy DTZ ground truth. "
                     "Piece-agnostic -- works for KQvK, KRvK, or any future endgame.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument('--syzygy', required=True,
                     help="Syzygy DTZ positions file, e.g. kqvk_positions_by_dtz.txt or "
                          "krvk_positions_by_dtz.txt")
    ap.add_argument('--landscape', required=True,
                     help="Generated landscape database, e.g. kqvk_perfect_play.db or "
                          "krvk_perfect_play.db")
    ap.add_argument('--out-dir', default='validation',
                     help="Directory to write mismatches.csv/matches.csv into (default: "
                          "'validation'). Give this explicitly when validating more than one "
                          "endgame from the same working directory.")
    args = ap.parse_args()

    print("=" * 120)
    print("Landscape Validation Against Syzygy (with Full Endgame Lines)")
    print("=" * 120 + "\n")

    if not Path(args.syzygy).exists():
        print(f"ERROR: {args.syzygy} not found")
        return

    if not Path(args.landscape).exists():
        print(f"ERROR: {args.landscape} not found")
        return

    os.makedirs(args.out_dir, exist_ok=True)

    validator = ValidationAnalyzer(args.syzygy, args.landscape)
    validator.validate()
    validator.analyze_mismatches()

    validator.write_mismatches(os.path.join(args.out_dir, "mismatches_full_games.csv"))
    validator.write_matches(os.path.join(args.out_dir, "matches_full_games.csv"))

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
