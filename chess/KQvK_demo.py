"""
Topological Chess Engine - Full Endgame Comparison
Compare topological search vs traditional minimax on complete games
"""

from dataclasses import dataclass
from typing import List, Tuple, Optional
import time
import random

@dataclass(frozen=True)
class Position:
    file: int  # 0-7 (a-h)
    rank: int  # 0-7 (1-8)
    
    def __str__(self):
        return f"{chr(ord('a') + self.file)}{self.rank + 1}"
    
    @staticmethod
    def from_str(s: str) -> 'Position':
        return Position(ord(s[0]) - ord('a'), int(s[1]) - 1)
    
    def distance_to(self, other: 'Position') -> int:
        """Chebyshev distance (king distance)"""
        return max(abs(self.file - other.file), abs(self.rank - other.rank))

@dataclass(frozen=True)
class GameState:
    wk: Position
    wq: Position
    bk: Position
    to_move: str
    
    def __str__(self):
        return f"WK:{self.wk} WQ:{self.wq} BK:{self.bk} {self.to_move}"

class BaseEngine:
    """Base class with common functionality"""
    
    def __init__(self):
        self.nodes_searched = 0
        self.moves_generated = 0
        
    def distance_to_nearest_edge(self, pos: Position) -> int:
        return min(pos.file, 7 - pos.file, pos.rank, 7 - pos.rank)
    
    def is_on_edge(self, pos: Position) -> bool:
        return pos.file in [0, 7] or pos.rank in [0, 7]
    
    def generate_all_king_moves(self, pos: Position) -> List[Position]:
        moves = []
        for df in [-1, 0, 1]:
            for dr in [-1, 0, 1]:
                if df == 0 and dr == 0:
                    continue
                new_file = pos.file + df
                new_rank = pos.rank + dr
                if 0 <= new_file <= 7 and 0 <= new_rank <= 7:
                    moves.append(Position(new_file, new_rank))
        return moves
    
    def generate_all_queen_moves(self, pos: Position) -> List[Position]:
        moves = []
        directions = [(1,0), (-1,0), (0,1), (0,-1), (1,1), (1,-1), (-1,1), (-1,-1)]
        for df, dr in directions:
            for dist in range(1, 8):
                new_file = pos.file + df * dist
                new_rank = pos.rank + dr * dist
                if 0 <= new_file <= 7 and 0 <= new_rank <= 7:
                    moves.append(Position(new_file, new_rank))
                else:
                    break
        return moves
    
    def is_attacked_by_queen(self, pos: Position, queen_pos: Position) -> bool:
        if pos == queen_pos:
            return False
        if pos.file == queen_pos.file or pos.rank == queen_pos.rank:
            return True
        if abs(pos.file - queen_pos.file) == abs(pos.rank - queen_pos.rank):
            return True
        return False
    
    def is_legal_state(self, state: GameState) -> bool:
        positions = [state.wk, state.wq, state.bk]
        if len(set(positions)) != 3:
            return False
        if state.wk.distance_to(state.bk) < 2:
            return False
        return True
    
    def is_checkmate(self, state: GameState) -> bool:
        if state.to_move != "BLACK":
            return False
        
        if not self.is_attacked_by_queen(state.bk, state.wq):
            return False
        
        for move in self.generate_all_king_moves(state.bk):
            if self.is_attacked_by_queen(move, state.wq):
                continue
            if move.distance_to(state.wk) <= 1:
                continue
            return False
        
        return True
    
    def is_stalemate(self, state: GameState) -> bool:
        if state.to_move != "BLACK":
            return False
        
        if self.is_attacked_by_queen(state.bk, state.wq):
            return False
        
        for move in self.generate_all_king_moves(state.bk):
            if self.is_attacked_by_queen(move, state.wq):
                continue
            if move.distance_to(state.wk) <= 1:
                continue
            return False
        
        return True
    
    def has_legal_black_moves(self, state: GameState) -> bool:
        """Check if black has any legal moves"""
        if state.to_move != "BLACK":
            return True
        
        for move in self.generate_all_king_moves(state.bk):
            if self.is_attacked_by_queen(move, state.wq):
                continue
            if move.distance_to(state.wk) <= 1:
                continue
            return True
        
        return False

class TraditionalEngine(BaseEngine):
    """Traditional minimax without topological guidance"""
    
    def minimax(self, state: GameState, depth: int, alpha: float = float('-inf'), beta: float = float('inf')) -> Tuple[Optional[int], Optional[GameState]]:
        self.nodes_searched += 1
        
        if self.is_checkmate(state):
            return (0, None)
        
        if self.is_stalemate(state):
            return (None, None)
        
        if depth == 0:
            return (None, None)
        
        best_mate = None
        best_state = None
        
        if state.to_move == "WHITE":
            candidates = []
            
            for wk_new in self.generate_all_king_moves(state.wk):
                new_state = GameState(wk_new, state.wq, state.bk, "BLACK")
                if not self.is_legal_state(new_state):
                    continue
                if self.is_stalemate(new_state):
                    continue
                candidates.append(new_state)
            
            for wq_new in self.generate_all_queen_moves(state.wq):
                new_state = GameState(state.wk, wq_new, state.bk, "BLACK")
                if not self.is_legal_state(new_state):
                    continue
                if self.is_stalemate(new_state):
                    continue
                candidates.append(new_state)
            
            self.moves_generated += len(candidates)
            
            for new_state in candidates:
                mate_dist, _ = self.minimax(new_state, depth - 1, alpha, beta)
                
                if mate_dist is not None:
                    actual_dist = mate_dist + 1
                    if best_mate is None or actual_dist < best_mate:
                        best_mate = actual_dist
                        best_state = new_state
                        alpha = max(alpha, -best_mate)
                        if alpha >= beta:
                            break
        
        else:  # BLACK
            candidates = []
            
            for bk_new in self.generate_all_king_moves(state.bk):
                if self.is_attacked_by_queen(bk_new, state.wq):
                    continue
                if bk_new.distance_to(state.wk) <= 1:
                    continue
                
                new_state = GameState(state.wk, state.wq, bk_new, "WHITE")
                if not self.is_legal_state(new_state):
                    continue
                candidates.append(new_state)
            
            self.moves_generated += len(candidates)
            
            for new_state in candidates:
                mate_dist, _ = self.minimax(new_state, depth - 1, alpha, beta)
                
                if mate_dist is not None:
                    actual_dist = mate_dist + 1
                    if best_mate is None or actual_dist > best_mate:
                        best_mate = actual_dist
                        best_state = new_state
                        beta = min(beta, best_mate)
                        if alpha >= beta:
                            break
        
        return (best_mate, best_state)
    
    def find_best_move(self, state: GameState, max_depth: int = 10) -> Tuple[Optional[GameState], Optional[int]]:
        self.nodes_searched = 0
        self.moves_generated = 0
        
        best_state = None
        best_mate = None
        
        for depth in range(2, max_depth + 1, 2):
            mate_dist, state_at_depth = self.minimax(state, depth)
            
            if mate_dist is not None:
                best_mate = mate_dist
                best_state = state_at_depth
                break
            
            if state_at_depth is not None:
                best_state = state_at_depth
        
        return (best_state, best_mate)

class TopologicalEngine(BaseEngine):
    """Topological distance-guided engine"""
    
    def __init__(self):
        super().__init__()
        self.moves_pruned = 0
    
    def get_target_edge(self, bk: Position) -> dict:
        dist_to_rank_1 = bk.rank
        dist_to_rank_8 = 7 - bk.rank
        dist_to_file_a = bk.file
        dist_to_file_h = 7 - bk.file
        
        min_dist = min(dist_to_rank_1, dist_to_rank_8, dist_to_file_a, dist_to_file_h)
        
        if min_dist == dist_to_rank_1:
            return {'type': 'rank', 'value': 0, 'cutoff': 1, 'distance': dist_to_rank_1}
        elif min_dist == dist_to_rank_8:
            return {'type': 'rank', 'value': 7, 'cutoff': 6, 'distance': dist_to_rank_8}
        elif min_dist == dist_to_file_a:
            return {'type': 'file', 'value': 0, 'cutoff': 1, 'distance': dist_to_file_a}
        else:
            return {'type': 'file', 'value': 7, 'cutoff': 6, 'distance': dist_to_file_h}
    
    def compute_topological_distance(self, state: GameState) -> float:
        distance = 0.0
        
        bk = state.bk
        wk = state.wk
        wq = state.wq
        
        edge_distance = self.distance_to_nearest_edge(bk)
        distance += edge_distance * 10.0
        
        wk_bk_distance = wk.distance_to(bk)
        optimal_wk_distance = 2
        wk_positioning_error = abs(wk_bk_distance - optimal_wk_distance)
        distance += wk_positioning_error * 5.0
        
        target_edge = self.get_target_edge(bk)
        
        if edge_distance == 0:
            if target_edge['type'] == 'rank':
                cutoff_rank = target_edge['cutoff']
                wq_cutoff_distance = abs(wq.rank - cutoff_rank)
                distance += wq_cutoff_distance * 3.0
                file_alignment = abs(wq.file - bk.file)
                distance += min(file_alignment, 3) * 1.0
            else:
                cutoff_file = target_edge['cutoff']
                wq_cutoff_distance = abs(wq.file - cutoff_file)
                distance += wq_cutoff_distance * 3.0
                rank_alignment = abs(wq.rank - bk.rank)
                distance += min(rank_alignment, 3) * 1.0
        
        elif edge_distance == 1:
            if target_edge['type'] == 'rank':
                cutoff_rank = target_edge['cutoff']
                wq_cutoff_distance = abs(wq.rank - cutoff_rank)
                
                if target_edge['value'] == 0:
                    if wq.rank > bk.rank:
                        wq_cutoff_distance = min(wq_cutoff_distance, 0.5)
                else:
                    if wq.rank < bk.rank:
                        wq_cutoff_distance = min(wq_cutoff_distance, 0.5)
                
                distance += wq_cutoff_distance * 2.0
            else:
                cutoff_file = target_edge['cutoff']
                wq_cutoff_distance = abs(wq.file - cutoff_file)
                
                if target_edge['value'] == 0:
                    if wq.file > bk.file:
                        wq_cutoff_distance = min(wq_cutoff_distance, 0.5)
                else:
                    if wq.file < bk.file:
                        wq_cutoff_distance = min(wq_cutoff_distance, 0.5)
                
                distance += wq_cutoff_distance * 2.0
        
        else:
            if target_edge['type'] == 'rank':
                if target_edge['value'] == 0:
                    if wq.rank <= bk.rank:
                        distance += (bk.rank - wq.rank + 1) * 2.0
                    else:
                        separation = wq.rank - bk.rank
                        distance += min(separation - 1, 2) * 0.5
                else:
                    if wq.rank >= bk.rank:
                        distance += (wq.rank - bk.rank + 1) * 2.0
                    else:
                        separation = bk.rank - wq.rank
                        distance += min(separation - 1, 2) * 0.5
            else:
                if target_edge['value'] == 0:
                    if wq.file <= bk.file:
                        distance += (bk.file - wq.file + 1) * 2.0
                    else:
                        separation = wq.file - bk.file
                        distance += min(separation - 1, 2) * 0.5
                else:
                    if wq.file >= bk.file:
                        distance += (wq.file - bk.file + 1) * 2.0
                    else:
                        separation = bk.file - wq.file
                        distance += min(separation - 1, 2) * 0.5
        
        return distance
    
    def minimax_topological(self, state: GameState, depth: int, current_distance: float) -> Tuple[Optional[int], Optional[GameState], float]:
        self.nodes_searched += 1
        
        if self.is_checkmate(state):
            return (0, None, 0.0)
        
        if self.is_stalemate(state):
            return (None, None, float('inf'))
        
        if depth == 0:
            return (None, None, self.compute_topological_distance(state))
        
        best_mate = None
        best_state = None
        best_distance = float('inf')
        
        if state.to_move == "WHITE":
            candidates = []
            
            for wk_new in self.generate_all_king_moves(state.wk):
                new_state = GameState(wk_new, state.wq, state.bk, "BLACK")
                if not self.is_legal_state(new_state):
                    continue
                if self.is_stalemate(new_state):
                    continue
                topo_dist = self.compute_topological_distance(new_state)
                candidates.append((new_state, topo_dist))
            
            for wq_new in self.generate_all_queen_moves(state.wq):
                new_state = GameState(state.wk, wq_new, state.bk, "BLACK")
                if not self.is_legal_state(new_state):
                    continue
                if self.is_stalemate(new_state):
                    continue
                topo_dist = self.compute_topological_distance(new_state)
                candidates.append((new_state, topo_dist))
            
            self.moves_generated += len(candidates)
            
            if candidates:
                min_candidate_distance = min(c[1] for c in candidates)
                pruning_threshold = 5.0
                candidates_filtered = [c for c in candidates if c[1] <= min_candidate_distance + pruning_threshold]
                self.moves_pruned += len(candidates) - len(candidates_filtered)
                candidates = candidates_filtered
            
            candidates.sort(key=lambda x: x[1])
            
            for new_state, topo_dist in candidates:
                mate_dist, _, final_dist = self.minimax_topological(new_state, depth - 1, topo_dist)
                
                if mate_dist is not None:
                    actual_dist = mate_dist + 1
                    if best_mate is None or actual_dist < best_mate:
                        best_mate = actual_dist
                        best_state = new_state
                        best_distance = topo_dist
                elif best_mate is None:
                    if topo_dist < best_distance:
                        best_distance = topo_dist
                        best_state = new_state
        
        else:  # BLACK
            candidates = []
            
            for bk_new in self.generate_all_king_moves(state.bk):
                if self.is_attacked_by_queen(bk_new, state.wq):
                    continue
                if bk_new.distance_to(state.wk) <= 1:
                    continue
                
                new_state = GameState(state.wk, state.wq, bk_new, "WHITE")
                if not self.is_legal_state(new_state):
                    continue
                
                topo_dist = self.compute_topological_distance(new_state)
                candidates.append((new_state, topo_dist))
            
            self.moves_generated += len(candidates)
            
            if not candidates:
                return (None, None, float('inf'))
            
            max_candidate_distance = max(c[1] for c in candidates)
            pruning_threshold = 5.0
            candidates_filtered = [c for c in candidates if c[1] >= max_candidate_distance - pruning_threshold]
            self.moves_pruned += len(candidates) - len(candidates_filtered)
            candidates = candidates_filtered
            
            candidates.sort(key=lambda x: x[1], reverse=True)
            
            for new_state, topo_dist in candidates:
                mate_dist, _, final_dist = self.minimax_topological(new_state, depth - 1, topo_dist)
                
                if mate_dist is not None:
                    actual_dist = mate_dist + 1
                    if best_mate is None or actual_dist > best_mate:
                        best_mate = actual_dist
                        best_state = new_state
                        best_distance = topo_dist
                elif best_mate is None:
                    if topo_dist > best_distance:
                        best_distance = topo_dist
                        best_state = new_state
        
        return (best_mate, best_state, best_distance)
    
    def find_best_move(self, state: GameState, max_depth: int = 10) -> Tuple[Optional[GameState], Optional[int]]:
        self.nodes_searched = 0
        self.moves_generated = 0
        self.moves_pruned = 0
        
        current_distance = self.compute_topological_distance(state)
        best_state = None
        best_mate = None
        
        for depth in range(2, max_depth + 1, 2):
            mate_dist, state_at_depth, final_dist = self.minimax_topological(state, depth, current_distance)
            
            if mate_dist is not None:
                best_mate = mate_dist
                best_state = state_at_depth
                break
            
            if state_at_depth is not None:
                best_state = state_at_depth
        
        return (best_state, best_mate)

# ==================== GAME PLAYER ====================

def play_full_game(engine_white, engine_black, initial_state: GameState, max_moves: int = 50, timeout_per_move: float = 900.0) -> dict:
    """
    Play a complete game and record all statistics with timeout protection.
    """
    state = initial_state
    moves = []
    total_time_white = 0
    total_time_black = 0
    total_nodes_white = 0
    total_nodes_black = 0
    total_moves_generated_white = 0
    total_moves_generated_black = 0
    
    for move_num in range(max_moves):
        if engine_white.is_checkmate(state):
            return {
                'moves': moves,
                'total_time_white': total_time_white,
                'total_time_black': total_time_black,
                'total_nodes_white': total_nodes_white,
                'total_nodes_black': total_nodes_black,
                'total_moves_generated_white': total_moves_generated_white,
                'total_moves_generated_black': total_moves_generated_black,
                'winner': 'WHITE',
                'move_count': len(moves),
                'outcome': 'checkmate'
            }
        
        if engine_white.is_stalemate(state):
            return {
                'moves': moves,
                'total_time_white': total_time_white,
                'total_time_black': total_time_black,
                'total_nodes_white': total_nodes_white,
                'total_nodes_black': total_nodes_black,
                'total_moves_generated_white': total_moves_generated_white,
                'total_moves_generated_black': total_moves_generated_black,
                'winner': 'DRAW',
                'move_count': len(moves),
                'outcome': 'stalemate'
            }
        
        # White's turn
        if state.to_move == "WHITE":
            start = time.time()
            
            try:
                new_state, mate = engine_white.find_best_move(state, max_depth=10)
                elapsed = time.time() - start
                
                if elapsed > timeout_per_move:
                    print(f"  WARNING: White move took {elapsed:.1f}s (timeout: {timeout_per_move}s)")
                    return {
                        'moves': moves,
                        'total_time_white': total_time_white + elapsed,
                        'total_time_black': total_time_black,
                        'total_nodes_white': total_nodes_white,
                        'total_nodes_black': total_nodes_black,
                        'total_moves_generated_white': total_moves_generated_white,
                        'total_moves_generated_black': total_moves_generated_black,
                        'winner': 'TIMEOUT',
                        'move_count': len(moves),
                        'outcome': 'timeout_white'
                    }
                
            except Exception as e:
                print(f"  ERROR in White move: {e}")
                return {
                    'moves': moves,
                    'total_time_white': total_time_white,
                    'total_time_black': total_time_black,
                    'total_nodes_white': total_nodes_white,
                    'total_nodes_black': total_nodes_black,
                    'total_moves_generated_white': total_moves_generated_white,
                    'total_moves_generated_black': total_moves_generated_black,
                    'winner': 'ERROR',
                    'move_count': len(moves),
                    'outcome': 'error_white'
                }
            
            if new_state is None:
                return {
                    'moves': moves,
                    'total_time_white': total_time_white,
                    'total_time_black': total_time_black,
                    'total_nodes_white': total_nodes_white,
                    'total_nodes_black': total_nodes_black,
                    'total_moves_generated_white': total_moves_generated_white,
                    'total_moves_generated_black': total_moves_generated_black,
                    'winner': 'BLACK',
                    'move_count': len(moves),
                    'outcome': 'no_legal_moves'
                }
            
            total_time_white += elapsed
            total_nodes_white += engine_white.nodes_searched
            total_moves_generated_white += engine_white.moves_generated
            
            moves.append({
                'move_number': len(moves) + 1,
                'player': 'WHITE',
                'from_state': str(state),
                'to_state': str(new_state),
                'time': elapsed,
                'nodes': engine_white.nodes_searched,
                'moves_generated': engine_white.moves_generated
            })
            
            state = new_state
        
        # Black's turn
        else:
            start = time.time()
            
            try:
                new_state, mate = engine_black.find_best_move(state, max_depth=10)
                elapsed = time.time() - start
                
                if elapsed > timeout_per_move:
                    print(f"  WARNING: Black move took {elapsed:.1f}s")
                    return {
                        'moves': moves,
                        'total_time_white': total_time_white,
                        'total_time_black': total_time_black + elapsed,
                        'total_nodes_white': total_nodes_white,
                        'total_nodes_black': total_nodes_black,
                        'total_moves_generated_white': total_moves_generated_white,
                        'total_moves_generated_black': total_moves_generated_black,
                        'winner': 'TIMEOUT',
                        'move_count': len(moves),
                        'outcome': 'timeout_black'
                    }
            
            except Exception as e:
                print(f"  ERROR in Black move: {e}")
                return {
                    'moves': moves,
                    'total_time_white': total_time_white,
                    'total_time_black': total_time_black,
                    'total_nodes_white': total_nodes_white,
                    'total_nodes_black': total_nodes_black,
                    'total_moves_generated_white': total_moves_generated_white,
                    'total_moves_generated_black': total_moves_generated_black,
                    'winner': 'ERROR',
                    'move_count': len(moves),
                    'outcome': 'error_black'
                }
            
            if new_state is None:
                return {
                    'moves': moves,
                    'total_time_white': total_time_white,
                    'total_time_black': total_time_black,
                    'total_nodes_white': total_nodes_white,
                    'total_nodes_black': total_nodes_black,
                    'total_moves_generated_white': total_moves_generated_white,
                    'total_moves_generated_black': total_moves_generated_black,
                    'winner': 'WHITE',
                    'move_count': len(moves),
                    'outcome': 'no_legal_moves'
                }
            
            total_time_black += elapsed
            total_nodes_black += engine_black.nodes_searched
            total_moves_generated_black += engine_black.moves_generated
            
            moves.append({
                'move_number': len(moves) + 1,
                'player': 'BLACK',
                'from_state': str(state),
                'to_state': str(new_state),
                'time': elapsed,
                'nodes': engine_black.nodes_searched,
                'moves_generated': engine_black.moves_generated
            })
            
            state = new_state
    
    return {
        'moves': moves,
        'total_time_white': total_time_white,
        'total_time_black': total_time_black,
        'total_nodes_white': total_nodes_white,
        'total_nodes_black': total_nodes_black,
        'total_moves_generated_white': total_moves_generated_white,
        'total_moves_generated_black': total_moves_generated_black,
        'winner': 'DRAW',
        'move_count': len(moves),
        'outcome': 'max_moves_reached'
    }

def generate_random_position() -> GameState:
    """Generate a random legal KQvK position with good separation"""
    max_attempts = 1000
    
    for _ in range(max_attempts):
        wk = Position(random.randint(0, 7), random.randint(0, 7))
        wq = Position(random.randint(0, 7), random.randint(0, 7))
        bk = Position(random.randint(0, 7), random.randint(0, 7))
        
        # Check basic legality
        if len(set([wk, wq, bk])) != 3:
            continue
        if wk.distance_to(bk) < 2:
            continue
        
        # Ensure reasonable separation (avoid trivial positions)
        # WK should be 2-5 squares from BK
        if wk.distance_to(bk) > 5:
            continue
        
        # WQ should not be immediately adjacent to BK (avoid near-mate positions)
        if wq.distance_to(bk) < 2:
            continue
        
        state = GameState(wk, wq, bk, "WHITE")
        
        # Check not already checkmate or stalemate
        engine = BaseEngine()
        if engine.is_checkmate(state) or engine.is_stalemate(state):
            continue
        
        # Ensure black has legal moves
        if not engine.has_legal_black_moves(GameState(wk, wq, bk, "BLACK")):
            continue
        
        return state
    
    # Fallback to a known good position
    return GameState(Position.from_str("e4"), Position.from_str("d1"), Position.from_str("h8"), "WHITE")

def print_game_summary(result: dict, title: str):
    """Print a summary of a game"""
    print(f"\n{'='*80}")
    print(f"{title}")
    print(f"{'='*80}")
    
    print(f"\nOutcome: {result['outcome'].upper()} - Winner: {result['winner']}")
    print(f"Total moves: {result['move_count']}")
    
    print(f"\n[WHITE Statistics]")
    print(f"  Total time: {result['total_time_white']:.3f}s")
    print(f"  Total nodes: {result['total_nodes_white']:,}")
    print(f"  Total moves generated: {result['total_moves_generated_white']:,}")
    if result['move_count'] > 0:
        white_moves = result['move_count'] // 2 + result['move_count'] % 2
        print(f"  Avg time/move: {result['total_time_white']/white_moves:.3f}s")
        print(f"  Avg nodes/move: {result['total_nodes_white']//white_moves:,}")
    
    print(f"\n[BLACK Statistics]")
    print(f"  Total time: {result['total_time_black']:.3f}s")
    print(f"  Total nodes: {result['total_nodes_black']:,}")
    print(f"  Total moves generated: {result['total_moves_generated_black']:,}")
    if result['move_count'] > 1:
        black_moves = result['move_count'] // 2
        if black_moves > 0:
            print(f"  Avg time/move: {result['total_time_black']/black_moves:.3f}s")
            print(f"  Avg nodes/move: {result['total_nodes_black']//black_moves:,}")
    
    print(f"\n[Move Sequence]")
    for move in result['moves'][:10]:
        print(f"  {move['move_number']}. {move['player']}: {move['to_state']} "
              f"(time: {move['time']:.3f}s, nodes: {move['nodes']:,})")
    
    if len(result['moves']) > 10:
        print(f"  ... ({len(result['moves']) - 10} more moves)")

# ==================== MAIN COMPARISON ====================

if __name__ == "__main__":
    print("KQvK Endgame Engine Comparison")
    print("="*80)
    
    # Test 1: Known positions
    print("\n\nTEST 1: Traditional White vs Topological Black")
    print("Starting position: WK:c3 WQ:b8 BK:f1")
    
    state1 = GameState(Position.from_str("c3"), Position.from_str("b8"), Position.from_str("f1"), "WHITE")
    state2 = GameState(Position.from_str("c3"), Position.from_str("b8"), Position.from_str("f2"), "WHITE")
    
    traditional_white = TraditionalEngine()
    topological_black = TopologicalEngine()
    
    result1 = play_full_game(traditional_white, topological_black, state1)
    print_game_summary(result1, "Traditional White vs Topological Black (f1)")
    
    traditional_white2 = TraditionalEngine()
    topological_black2 = TopologicalEngine()
    result2 = play_full_game(traditional_white2, topological_black2, state2)
    print_game_summary(result2, "Traditional White vs Topological Black (f2)")
    
    # Test 2: Topological vs Topological
    print("\n\nTEST 2: Topological White vs Topological Black")
    
    topological_white = TopologicalEngine()
    topological_black3 = TopologicalEngine()
    result3 = play_full_game(topological_white, topological_black3, state1)
    print_game_summary(result3, "Topological vs Topological (f1)")
    
    topological_white2 = TopologicalEngine()
    topological_black4 = TopologicalEngine()
    result4 = play_full_game(topological_white2, topological_black4, state2)
    print_game_summary(result4, "Topological vs Topological (f2)")
    
    # Test 3: THE HARD POSITION
    print("\n\nTEST 3: THE HARD POSITION - WK:g4 WQ:b6 BK:c4")
    hard_state = GameState(Position.from_str("g4"), Position.from_str("b6"), Position.from_str("c4"), "WHITE")
    hard_state2 = GameState(Position.from_str("c8"), Position.from_str("d5"), Position.from_str("b6"), "WHITE")
    
    print("\n[Topological Engine on HARD position]")
    topo_hard_w = TopologicalEngine()
    topo_hard_b = TopologicalEngine()
    result_hard_topo = play_full_game(topo_hard_w, topo_hard_b, hard_state, timeout_per_move=900.0)
    print_game_summary(result_hard_topo, "Topological vs Topological (HARD: BK in center)")
    tradhard1 = TraditionalEngine()
    tradhard2 = TraditionalEngine()
    result_hard_topo = play_full_game(tradhard1, tradhard2, hard_state, timeout_per_move=3600.0)
    print_game_summary(result_hard_topo, "Topological vs Topological (HARD: BK in center)")
    
    print("\n\nTEST 3: THE HARD POSITION - WK:c8 WQ:d5 BK:b6")

    result_hard_topo = play_full_game(topo_hard_w, topo_hard_b, hard_state2, timeout_per_move=900.0)
    print_game_summary(result_hard_topo, "Topological vs Topological (WK:c8 WQ:d5 BK:b6)")
    
    result_hard_topo = play_full_game(tradhard1, tradhard2, hard_state2, timeout_per_move=3600.0)
    print_game_summary(result_hard_topo, "Traditional vs Topological (WK:c8 WQ:d5 BK:b6)")
    
    # Test 4: Random positions
    print("\n\nTEST 4: Random Position Tests (5 games each)")
    
    traditional_results = []
    topological_results = []
    
    for i in range(5):
        random_state = generate_random_position()
        
        print(f"\n--- Random Game {i+1} ---")
        print(f"Position: {random_state}")
        
        # Traditional
        trad_w = TraditionalEngine()
        trad_b = TopologicalEngine()
        result_trad = play_full_game(trad_w, trad_b, random_state, timeout_per_move=900.0)
        traditional_results.append(result_trad)
        print(f"Traditional: {result_trad['move_count']} moves, "
              f"{result_trad['total_time_white']:.2f}s, "
              f"{result_trad['total_nodes_white']:,} nodes, "
              f"outcome: {result_trad['outcome']}")
        
        # Topological
        topo_w = TopologicalEngine()
        topo_b = TopologicalEngine()
        result_topo = play_full_game(topo_w, topo_b, random_state, timeout_per_move=900.0)
        topological_results.append(result_topo)
        print(f"Topological: {result_topo['move_count']} moves, "
              f"{result_topo['total_time_white']:.2f}s, "
              f"{result_topo['total_nodes_white']:,} nodes, "
              f"outcome: {result_topo['outcome']}")
    
    # Summary statistics (only count successful checkmates)
    print(f"\n\n{'='*80}")
    print("SUMMARY STATISTICS (5 Random Games)")
    print(f"{'='*80}")
    
    trad_successful = [r for r in traditional_results if r['outcome'] == 'checkmate']
    topo_successful = [r for r in topological_results if r['outcome'] == 'checkmate']
    
    print(f"\n[TRADITIONAL ENGINE] ({len(trad_successful)}/5 successful)")
    if trad_successful:
        avg_moves_trad = sum(r['move_count'] for r in trad_successful) / len(trad_successful)
        avg_time_trad = sum(r['total_time_white'] + r['total_time_black'] for r in trad_successful) / len(trad_successful)
        avg_nodes_trad = sum(r['total_nodes_white'] + r['total_nodes_black'] for r in trad_successful) / len(trad_successful)
        
        print(f"Average moves to mate: {avg_moves_trad:.1f}")
        print(f"Average total time: {avg_time_trad:.3f}s")
        print(f"Average total nodes: {avg_nodes_trad:,.0f}")
    
    print(f"\n[TOPOLOGICAL ENGINE] ({len(topo_successful)}/5 successful)")
    if topo_successful:
        avg_moves_topo = sum(r['move_count'] for r in topo_successful) / len(topo_successful)
        avg_time_topo = sum(r['total_time_white'] + r['total_time_black'] for r in topo_successful) / len(topo_successful)
        avg_nodes_topo = sum(r['total_nodes_white'] + r['total_nodes_black'] for r in topo_successful) / len(topo_successful)
        
        print(f"Average moves to mate: {avg_moves_topo:.1f}")
        print(f"Average total time: {avg_time_topo:.3f}s")
        print(f"Average total nodes: {avg_nodes_topo:,.0f}")
        
        if trad_successful and topo_successful:
            print("\n[COMPARISON]")
            print(f"Topological speedup: {avg_time_trad/avg_time_topo:.2f}x faster")
            print(f"Topological node reduction: {avg_nodes_trad/avg_nodes_topo:.2f}x fewer nodes")
            print(f"Move efficiency: {avg_moves_trad/avg_moves_topo:.2f}x (lower = fewer moves)")
