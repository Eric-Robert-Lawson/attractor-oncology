#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <string>
#include <cmath>
#include <algorithm>
#include <tuple>
#include <optional>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <limits>
#include <fstream>

using namespace std;

const int INF_MAX = 2147483647;
const int INF_MIN = -2147483648;

// ============================================================================
// SolvedPosition - Data structure for storing perfect play
// ============================================================================

struct SolvedPosition {
    string position_key;        // "WK:a1 WQ:b2 BK:d5"
    char turn;                  // white or black
    string best_move;           // ""
    int M_value;                // Moves to mate
    vector<int> BN_trajectory;  // Black node counts per ply
    int total_plies;            // How many plies until mate
    int white_moves;            // White's move count
    int black_moves;            // Black's move count
    int nodes_evaluated;        // Nodes searched
    double computation_time;    // Time in seconds
    
    // Convert to CSV line for export
    string to_csv() const {
        stringstream ss;
        ss << position_key << "|" 
           << turn << "|"
           << best_move << "|" 
           << M_value << "|"
           << total_plies << "|"
           << white_moves << "|"
           << black_moves << "|"
           << nodes_evaluated << "|"
           << fixed << setprecision(3) << computation_time << "|";
        
        // BN trajectory as comma-separated
        for (size_t i = 0; i < BN_trajectory.size(); i++) {
            if (i > 0) ss << ",";
            ss << BN_trajectory[i];
        }
        ss << "\n";
        return ss.str();
    }
    
    static string csv_header() {
        return "Position|Turn|BestMove|M|Plies|WhiteMoves|BlackMoves|NodesEval|Time|BN_Trajectory\n";
    }
};

// ============================================================================
// SolvedPositionDatabase - Manages the tablebase
// ============================================================================

class SolvedPositionDatabase {
    private:
        map<pair<string, char>, SolvedPosition> solved;
        string filename;
        
    public:
        SolvedPositionDatabase(const string& db_file = "kqvk_perfect_play.db")
            : filename(db_file) {
            load_from_file();
        }
        
        // Add a solved position to the database
        void add_position(const SolvedPosition& pos) {
            if (is_solved(pos.position_key, pos.turn)) {
                return;
            }
            solved[{pos.position_key, pos.turn}] = pos;
            
            // Helper lambdas for transformations
            auto rotate_square = [](const string& sq) {
                if (sq.length() < 2) return sq;
                int file = sq[0] - 'a';
                int rank = sq[1] - '1';
                int new_file = 7 - rank;
                int new_rank = file;
                return string(1, char('a' + new_file)) + char('1' + new_rank);
            };
            
            auto flip_vertical = [](const string& sq) {
                if (sq.length() < 2) return sq;
                int file = sq[0] - 'a';
                int rank = sq[1] - '1';
                int new_rank = 7 - rank;
                return string(1, char('a' + file)) + char('1' + new_rank);
            };
            
            auto flip_horizontal = [](const string& sq) {
                if (sq.length() < 2) return sq;
                int file = sq[0] - 'a';
                int rank = sq[1] - '1';
                int new_file = 7 - file;
                return string(1, char('a' + new_file)) + char('1' + rank);
            };
            
            auto flip_diagonal_a1h8 = [](const string& sq) {
                if (sq.length() < 2) return sq;
                int file = sq[0] - 'a';
                int rank = sq[1] - '1';
                return string(1, char('a' + rank)) + char('1' + file);
            };
            
            auto flip_diagonal_a8h1 = [](const string& sq) {
                if (sq.length() < 2) return sq;
                int file = sq[0] - 'a';
                int rank = sq[1] - '1';
                int new_file = 7 - rank;
                int new_rank = 7 - file;
                return string(1, char('a' + new_file)) + char('1' + new_rank);
            };
            
            // Helper macro to add transformed position
            auto add_transformed = [&](auto transform_func, const string& base_pos, const string& base_move) {
                size_t wk_pos = base_pos.find("WK:") + 3;
                size_t wq_pos = base_pos.find("WQ:") + 3;
                size_t bk_pos = base_pos.find("BK:") + 3;
                
                string wk_sq = base_pos.substr(wk_pos, 2);
                string wq_sq = base_pos.substr(wq_pos, 2);
                string bk_sq = base_pos.substr(bk_pos, 2);
                
                string new_wk = transform_func(wk_sq);
                string new_wq = transform_func(wq_sq);
                string new_bk = transform_func(bk_sq);
                
                string new_pos = "WK:" + new_wk + " WQ:" + new_wq + " BK:" + new_bk;
                
                string new_move = base_move;
                if (!new_move.empty() && new_move != "mate" && new_move != "checkmate") {
                    if (new_move.length() >= 3) {
                        char piece = new_move[0];
                        string dest = new_move.substr(1);
                        new_move = piece + transform_func(dest);
                    }
                }
                
                if (!is_solved(new_pos, pos.turn)) {
                    SolvedPosition transformed = pos;
                    transformed.position_key = new_pos;
                    transformed.best_move = new_move;
                    solved[{new_pos, pos.turn}] = transformed;
                }
            };
            
            // Generate 4 rotations
            string current_pos = pos.position_key;
            string current_move = pos.best_move;
            
            for (int rot = 1; rot < 4; rot++) {
                add_transformed(rotate_square, current_pos, current_move);
                
                // Update current for next rotation
                size_t wk_pos = current_pos.find("WK:") + 3;
                size_t wq_pos = current_pos.find("WQ:") + 3;
                size_t bk_pos = current_pos.find("BK:") + 3;
                
                string wk_sq = current_pos.substr(wk_pos, 2);
                string wq_sq = current_pos.substr(wq_pos, 2);
                string bk_sq = current_pos.substr(bk_pos, 2);
                
                string new_wk = rotate_square(wk_sq);
                string new_wq = rotate_square(wq_sq);
                string new_bk = rotate_square(bk_sq);
                
                current_pos = "WK:" + new_wk + " WQ:" + new_wq + " BK:" + new_bk;
                
                if (!current_move.empty() && current_move != "mate" && current_move != "checkmate") {
                    if (current_move.length() >= 3) {
                        char piece = current_move[0];
                        string dest = current_move.substr(1);
                        current_move = piece + rotate_square(dest);
                    }
                }
            }
            
            // Generate 4 reflections from original
            add_transformed(flip_vertical, pos.position_key, pos.best_move);
            add_transformed(flip_horizontal, pos.position_key, pos.best_move);
            add_transformed(flip_diagonal_a1h8, pos.position_key, pos.best_move);
            add_transformed(flip_diagonal_a8h1, pos.position_key, pos.best_move);
        }
    
    // Check if position is already solved
    bool is_solved(const string& position_key, char turn) const {
    return solved.count({position_key, turn}) > 0;
}
    
    // Get optimal move for a position
    optional<SolvedPosition> get_solution(const string& position_key, char turn) const {
        auto it = solved.find({position_key, turn});
        if (it != solved.end()) {
            return it->second;
        }
        return nullopt;
    }
    
    // Export all solved positions to file
    void export_to_file() {
        ofstream file(filename, ios::trunc);
        
        if (!file.is_open()) {
            cerr << "ERROR: Cannot open database file: " << filename << "\n";
            return;
        }
        
        file << SolvedPosition::csv_header();
        
        for (const auto& [key, pos] : solved) {
            file << pos.to_csv();
        }
        
        file.close();
    }
    
    // Load from file for future runs
    void load_from_file() {
        ifstream file(filename);
        
        if (!file.is_open()) {
            return;
        }
        
        string line;
        bool is_header = true;
        int line_num = 0;
        int loaded_count = 0;
        
        while (getline(file, line)) {
            line_num++;
            
            if (is_header) {
                is_header = false;
                continue;
            }
            
            if (line.empty()) {
                continue;
            }
            
            vector<string> parts;
            stringstream ss(line);
            string part;
            
            while (getline(ss, part, '|')) {
                parts.push_back(part);
            }
            
            if (parts.size() >= 8) {
                try {
                    SolvedPosition pos;
                    pos.position_key = parts[0];
                    pos.turn = parts[1].empty() ? 'W' : parts[1][0];
                    pos.best_move = parts[2];
                    pos.M_value = stoi(parts[3]);
                    pos.total_plies = stoi(parts[4]);
                    pos.white_moves = stoi(parts[5]);
                    pos.black_moves = stoi(parts[6]);
                    pos.nodes_evaluated = stoi(parts[7]);
                    pos.computation_time = stod(parts[8]);
                    if (parts.size() > 9) {
                        if (!parts[9].empty()) {
                            stringstream bn_stream(parts[9]);
                            string bn_part;
                            int bn_count = 0;
                            
                            while (getline(bn_stream, bn_part, ',')) {
                                try {
                                    int bn_val = stoi(bn_part);
                                    pos.BN_trajectory.push_back(bn_val);
                                    bn_count++;
                                } catch (const exception& e) {
                                    cout << "    -> ERROR parsing BN value: " << e.what() << "\n";
                                }
                            }
                        }
                    }
                    
                    solved[{pos.position_key, pos.turn}] = pos;
                    loaded_count++;
                    // cout << "  ✓ Successfully loaded position\n";
                    
                } catch (const exception& e) {
                    cout << "  ✗ ERROR parsing line: " << e.what() << "\n";
                }
            } else {
                cout << "  ✗ Skipping: parts.size() (" << parts.size() << ") < 8\n";
            }
        }
        
        file.close();
        cout << "\n[Database] Finished loading. Loaded " << loaded_count << " positions from " << filename << "\n";
        cout << "[Database] Total in map: " << solved.size() << " positions\n";
    }
    
    // Print summary statistics
    void print_summary() const {
        cout << "\n" << string(80, '=') << "\n";
        cout << "TABLEBASE SUMMARY\n";
        cout << string(80, '=') << "\n";
        cout << "Total positions solved: " << solved.size() << "\n";
        cout << "Database file: " << filename << "\n";
        
        if (solved.empty()) return;
        
        int total_M = 0;
        int total_plies = 0;
        int total_nodes = 0;
        double total_time = 0;
        
        for (const auto& [key, pos] : solved) {
            total_M += pos.M_value;
            total_plies += pos.total_plies;
            total_nodes += pos.nodes_evaluated;
            total_time += pos.computation_time;
        }
        
        cout << "\nStatistics:\n";
        cout << "  Average M: " << fixed << setprecision(1) << (double)total_M / solved.size() << "\n";
        cout << "  Average plies to mate: " << (double)total_plies / solved.size() << "\n";
        cout << "  Total nodes evaluated: " << total_nodes << "\n";
        cout << "  Total computation time: " << setprecision(1) << total_time << "s\n";
        cout << "\n";
    }
};


// ============================================================================
// Position class
// ============================================================================

class Position {
public:
    int file, rank;
    Position() : file(0), rank(0) {}
    Position(int f, int r) : file(f), rank(r) {}
    
    static Position from_str(const string& s) {
        return Position(s[0] - 'a', stoi(s.substr(1)) - 1);
    }
    
    string str() const {
        string result;
        result += char('a' + file);
        result += char('1' + rank);
        return result;
    }
    
    int distance_to(const Position& other) const {
        return max(abs(file - other.file), abs(rank - other.rank));
    }
    
    bool operator==(const Position& other) const {
        return file == other.file && rank == other.rank;
    }
    
    bool operator!=(const Position& other) const {
        return !(*this == other);
    }
};

// ============================================================================
// GameState class
// ============================================================================

class GameState {
public:
    Position wk, wq, bk;
    char to_move;
    
    GameState() : to_move('W') {}
    GameState(Position wk_, Position wq_, Position bk_, char to_move_)
        : wk(wk_), wq(wq_), bk(bk_), to_move(to_move_) {}
    
    string str() const {
        return "WK:" + wk.str() + " WQ:" + wq.str() + " BK:" + bk.str();
    }
    
    bool operator==(const GameState& other) const {
        return wk == other.wk && wq == other.wq && bk == other.bk && to_move == other.to_move;
    }
};

// ============================================================================
// BaseEngine
// ============================================================================

class BaseEngine {
public:
    int distance_to_nearest_edge(const Position& pos) const {
        return min({pos.file, 7 - pos.file, pos.rank, 7 - pos.rank});
    }
    
    bool is_on_edge(const Position& pos) const {
        return pos.file == 0 || pos.file == 7 || pos.rank == 0 || pos.rank == 7;
    }
    
    inline vector<Position> generate_all_king_moves(const Position& pos) const {
        vector<Position> moves;
        for (int df = -1; df <= 1; df++) {
            for (int dr = -1; dr <= 1; dr++) {
                if (df == 0 && dr == 0) continue;
                int nf = pos.file + df;
                int nr = pos.rank + dr;
                if (nf >= 0 && nf <= 7 && nr >= 0 && nr <= 7) {
                    moves.push_back(Position(nf, nr));
                }
            }
        }
        return moves;
    }
    
    inline vector<Position> generate_all_queen_moves(const Position& pos) const {
        vector<Position> moves;
        int dirs[][2] = {{1,0},{-1,0},{0,1},{0,-1},{1,1},{1,-1},{-1,1},{-1,-1}};
        
        for (auto& d : dirs) {
            for (int dist = 1; dist <= 7; dist++) {
                int nf = pos.file + d[0] * dist;
                int nr = pos.rank + d[1] * dist;
                if (nf >= 0 && nf <= 7 && nr >= 0 && nr <= 7) {
                    moves.push_back(Position(nf, nr));
                } else break;
            }
        }
        return moves;
    }
    
    bool is_attacked_by_queen(const Position& pos, const Position& qp, const Position& wk, const Position& wq) const {
        if (pos == qp) return false;
        
        // Check file
        if (pos.file == qp.file) {
            int start = min(pos.rank, qp.rank) + 1;
            int end = max(pos.rank, qp.rank);
            for (int r = start; r < end; r++) {
                if (Position(pos.file, r) == wk || Position(pos.file, r) == wq) return false;  // Blocked
            }
            return true;
        }
        
        // Check rank
        if (pos.rank == qp.rank) {
            int start = min(pos.file, qp.file) + 1;
            int end = max(pos.file, qp.file);
            for (int f = start; f < end; f++) {
                if (Position(f, pos.rank) == wk || Position(f, pos.rank) == wq) return false;  // Blocked
            }
            return true;
        }
        
        // Check diagonals
        if (abs(pos.file - qp.file) == abs(pos.rank - qp.rank)) {
            int df = (pos.file > qp.file) ? 1 : -1;
            int dr = (pos.rank > qp.rank) ? 1 : -1;
            int f = qp.file + df;
            int r = qp.rank + dr;
            while (f != pos.file) {
                if (Position(f, r) == wk || Position(f, r) == wq) return false;  // Blocked
                f += df;
                r += dr;
            }
            return true;
        }
        
        return false;
    }
    
    bool is_legal_state(const GameState& st) const {
        set<pair<int,int>> pos;
        pos.insert({st.wk.file, st.wk.rank});
        pos.insert({st.wq.file, st.wq.rank});
        pos.insert({st.bk.file, st.bk.rank});
        if (pos.size() != 3) return false;
        if (st.wk.distance_to(st.bk) < 2) return false;
        return true;
    }
    
    bool is_checkmate(const GameState& st) const {
        if (st.to_move != 'B') return false;
        if (!is_attacked_by_queen(st.bk, st.wq, st.wk, st.wq)) return false;
        for (auto& m : generate_all_king_moves(st.bk)) {
            if (is_attacked_by_queen(m, st.wq, st.wk, st.wq)) continue;
            if (m.distance_to(st.wk) <= 1) continue;
            return false;
        }
        return true;
    }
    
    bool is_stalemate(const GameState& st) const {
        if (st.to_move != 'B') return false;
        if (is_attacked_by_queen(st.bk, st.wq, st.wk, st.wq)) return false;
        for (auto& m : generate_all_king_moves(st.bk)) {
            if (is_attacked_by_queen(m, st.wq, st.wk, st.wq)) continue;
            if (m.distance_to(st.wk) <= 1) continue;
            return false;
        }
        return true;
    }
    
    string get_move_notation(const GameState& from, const GameState& to) const {
        if (from.wk != to.wk) return "K" + to.wk.str();
        if (from.wq != to.wq) return "Q" + to.wq.str();
        if (from.bk != to.bk) return "k" + to.bk.str();
        return "??";
    }
    
    int count_legal_moves(const GameState& st) const {
        if (st.to_move == 'W') {
            int cnt = 0;
            for (auto& wk_n : generate_all_king_moves(st.wk)) {
                GameState ns(wk_n, st.wq, st.bk, 'B');
                if (is_legal_state(ns) && !is_stalemate(ns)) cnt++;
            }
            for (auto& wq_n : generate_all_queen_moves(st.wq)) {
                if (wq_n.distance_to(st.bk) < 2 && wq_n.distance_to(st.wk) > 1) continue;
                
                // Check if White King blocks the Queen's path
                bool blocked = false;
                if (wq_n.file == st.wq.file) {
                    int start = min(st.wq.rank, wq_n.rank) + 1;
                    int end = max(st.wq.rank, wq_n.rank);
                    for (int r = start; r < end; r++) {
                        if (Position(wq_n.file, r) == st.wk) { blocked = true; break; }
                    }
                } else if (wq_n.rank == st.wq.rank) {
                    int start = min(st.wq.file, wq_n.file) + 1;
                    int end = max(st.wq.file, wq_n.file);
                    for (int f = start; f < end; f++) {
                        if (Position(f, wq_n.rank) == st.wk) { blocked = true; break; }
                    }
                } else if (abs(wq_n.file - st.wq.file) == abs(wq_n.rank - st.wq.rank)) {
                    int df = (wq_n.file > st.wq.file) ? 1 : -1;
                    int dr = (wq_n.rank > st.wq.rank) ? 1 : -1;
                    int f = st.wq.file + df;
                    int r = st.wq.rank + dr;
                    while (f != wq_n.file) {
                        if (Position(f, r) == st.wk) { blocked = true; break; }
                        f += df;
                        r += dr;
                    }
                }
                if (blocked) continue;

                GameState ns(st.wk, wq_n, st.bk, 'B');
                if (is_legal_state(ns) && !is_stalemate(ns)) cnt++;
            }
            return cnt;
        } else {
            int cnt = 0;
            for (auto& bk_n : generate_all_king_moves(st.bk)) {
                if (is_attacked_by_queen(bk_n, st.wq, st.wk, st.wq)) continue;
                if (bk_n.distance_to(st.wk) <= 1) continue;
                GameState ns(st.wk, st.wq, bk_n, 'W');
                if (is_legal_state(ns)) cnt++;
            }
            return cnt;
        }
    }
};

// ============================================================================
// CompositionalEngine
// ============================================================================

class CompositionalEngine : public BaseEngine {
public:
    int nodes_evaluated = 0;
    int candidates_measured = 0;
    unordered_map<uint64_t, int> M_cache;
// Hash function for memoization - converts position + depth to uint64_t
    uint64_t make_cache_key(const GameState& st, int depth) const {
        uint64_t key = 0;
        key |= ((uint64_t)st.wk.file << 56);
        key |= ((uint64_t)st.wk.rank << 52);
        key |= ((uint64_t)st.wq.file << 48);
        key |= ((uint64_t)st.wq.rank << 44);
        key |= ((uint64_t)st.bk.file << 40);
        key |= ((uint64_t)st.bk.rank << 36);
        key |= ((uint64_t)depth << 32);
        return key;
    }

    int compute_M_topological(const GameState& st) const {
        int edge_dist = distance_to_nearest_edge(st.bk);
        int target_rank = (st.bk.rank <= 3) ? 2 : 5;
        int wk_dist = max(abs(st.wk.rank - target_rank), abs(st.wk.file - st.bk.file));
        int q_moves = is_on_edge(st.bk) ? 1 : 2;
        return max(edge_dist, wk_dist) + q_moves + 1;
    }
    
    int compute_M_shallow(const GameState& st, int depth = 3) const {
        if (is_checkmate(st)) return 0;
        if (depth == 0) return compute_M_topological(st);
        
        if (st.to_move == 'W') {
            vector<GameState> cands;
            cands.reserve(64);
            for (auto& wk_n : generate_all_king_moves(st.wk)) {
                GameState ns(wk_n, st.wq, st.bk, 'B');
                if (is_legal_state(ns) && !is_stalemate(ns)) cands.push_back(ns);
            }
            for (auto& wq_n : generate_all_queen_moves(st.wq)) {
                if (wq_n.distance_to(st.bk) < 2 && wq_n.distance_to(st.wk) > 1) continue;
                
                // Check if White King blocks the Queen's path
                bool blocked = false;
                if (wq_n.file == st.wq.file) {
                    int start = min(st.wq.rank, wq_n.rank) + 1;
                    int end = max(st.wq.rank, wq_n.rank);
                    for (int r = start; r < end; r++) {
                        if (Position(wq_n.file, r) == st.wk) { blocked = true; break; }
                    }
                } else if (wq_n.rank == st.wq.rank) {
                    int start = min(st.wq.file, wq_n.file) + 1;
                    int end = max(st.wq.file, wq_n.file);
                    for (int f = start; f < end; f++) {
                        if (Position(f, wq_n.rank) == st.wk) { blocked = true; break; }
                    }
                } else if (abs(wq_n.file - st.wq.file) == abs(wq_n.rank - st.wq.rank)) {
                    int df = (wq_n.file > st.wq.file) ? 1 : -1;
                    int dr = (wq_n.rank > st.wq.rank) ? 1 : -1;
                    int f = st.wq.file + df;
                    int r = st.wq.rank + dr;
                    while (f != wq_n.file) {
                        if (Position(f, r) == st.wk) { blocked = true; break; }
                        f += df;
                        r += dr;
                    }
                }
                if (blocked) continue;

                GameState ns(st.wk, wq_n, st.bk, 'B');
                if (is_legal_state(ns) && !is_stalemate(ns)) cands.push_back(ns);
            }
            if (cands.empty()) return 99;
            
            int best_M = INF_MAX;
            int res = INF_MAX;

            vector<int> M_values;  // Cache the M values
            for (auto& c : cands) {
                int M = 1 + compute_M_shallow(c, depth - 1);
                M_values.push_back(M);
                best_M = min(best_M, M);
            }

            int thresh = max(2, best_M / 3);
            for (int M : M_values) {
                if (M <= best_M + thresh) res = min(res, M);
            }
            return res;
        } else {
            vector<GameState> cands;
            cands.reserve(64);
            for (auto& bk_n : generate_all_king_moves(st.bk)) {
                if (is_attacked_by_queen(bk_n, st.wq, st.wk, st.wq)) continue;
                if (bk_n.distance_to(st.wk) <= 1) continue;
                GameState ns(st.wk, st.wq, bk_n, 'W');
                if (is_legal_state(ns)) cands.push_back(ns);
            }
            if (cands.empty()) return 1;
            
            int best_M = INF_MIN;
            for (auto& c : cands) {
                best_M = max(best_M, 1 + compute_M_shallow(c, depth - 1));
            }
            int thresh = max(2, best_M / 3);
            
            int res = INF_MIN;
            for (auto& c : cands) {
                int M = 1 + compute_M_shallow(c, depth - 1);
                if (M >= best_M - thresh) res = max(res, M);
            }
            return res;
        }
    }
    
    int measure_black_nodes_after_trajectory(const GameState& st, int depth = 3) const {
        GameState curr = st;
        for (int ply = 0; ply < depth; ply++) {
            if (is_checkmate(curr)) break;
            
            if (curr.to_move == 'W') {
                vector<GameState> cands;
                cands.reserve(64);
                for (auto& wk_n : generate_all_king_moves(curr.wk)) {
                    GameState ns(wk_n, curr.wq, curr.bk, 'B');
                    if (is_legal_state(ns) && !is_stalemate(ns)) cands.push_back(ns);
                }
                for (auto& wq_n : generate_all_queen_moves(curr.wq)) {
                    if (wq_n.distance_to(curr.bk) < 2 && wq_n.distance_to(curr.wk) > 1) continue;
                    
                    // Check if White King blocks the Queen's path
                    bool blocked = false;
                    if (wq_n.file == curr.wq.file) {
                        int start = min(curr.wq.rank, wq_n.rank) + 1;
                        int end = max(curr.wq.rank, wq_n.rank);
                        for (int r = start; r < end; r++) {
                            if (Position(wq_n.file, r) == curr.wk) { blocked = true; break; }
                        }
                    } else if (wq_n.rank == curr.wq.rank) {
                        int start = min(curr.wq.file, wq_n.file) + 1;
                        int end = max(curr.wq.file, wq_n.file);
                        for (int f = start; f < end; f++) {
                            if (Position(f, wq_n.rank) == curr.wk) { blocked = true; break; }
                        }
                    } else if (abs(wq_n.file - curr.wq.file) == abs(wq_n.rank - curr.wq.rank)) {
                        int df = (wq_n.file > curr.wq.file) ? 1 : -1;
                        int dr = (wq_n.rank > curr.wq.rank) ? 1 : -1;
                        int f = curr.wq.file + df;
                        int r = curr.wq.rank + dr;
                        while (f != wq_n.file) {
                            if (Position(f, r) == curr.wk) { blocked = true; break; }
                            f += df;
                            r += dr;
                        }
                    }
                    if (blocked) continue;

                    GameState ns(curr.wk, wq_n, curr.bk, 'B');
                    if (is_legal_state(ns) && !is_stalemate(ns)) cands.push_back(ns);
                }
                if (cands.empty()) break;
                
                GameState best = cands[0];
                // HERE A COUPLE BAD ESTIMATES
                int best_M = compute_M_topological(cands[0]);
                for (size_t i = 1; i < cands.size(); i++) {
                    int M = compute_M_topological(cands[i]);
                    if (M < best_M) { best_M = M; best = cands[i]; }
                }
                curr = best;
            } else {
                vector<GameState> cands;
                cands.reserve(64);
                for (auto& bk_n : generate_all_king_moves(curr.bk)) {
                    if (is_attacked_by_queen(bk_n, curr.wq, curr.wk, curr.wk)) continue;
                    if (bk_n.distance_to(curr.wk) <= 1) continue;
                    GameState ns(curr.wk, curr.wq, bk_n, 'W');
                    if (is_legal_state(ns)) cands.push_back(ns);
                }
                if (cands.empty()) break;
                
                GameState best = cands[0];
                // HERE A COUPLE OF BAD ESTIMATES
                int best_M = compute_M_topological(cands[0]);
                for (size_t i = 1; i < cands.size(); i++) {
                    int M = compute_M_topological(cands[i]);
                    if (M > best_M) { best_M = M; best = cands[i]; }
                }
                curr = best;
            }
        }
        return count_legal_moves(curr);
    }
    
    pair<optional<int>, optional<GameState>> compositional_search_impl(
        const GameState& st, int depth, int ply,
        bool debug, unordered_map<uint64_t, pair<optional<int>, optional<GameState>>>& memo
    ) {
        uint64_t cache_key = make_cache_key(st, depth);
        if (memo.count(cache_key)) {
            return memo[cache_key];
        }
        
        if (is_checkmate(st)) {
            auto res = make_pair(optional<int>(0), optional<GameState>());
            memo[cache_key] = res;
            return res;
        }
        
        if (depth == 0) {
            // HERE A COUPLE OF BAD ESTIMATES
            int estimated_M = compute_M_topological(st);
            auto res = make_pair(optional<int>(estimated_M), optional<GameState>());
            memo[cache_key] = res;
            return res;
        }
        
        vector<GameState> cands;
        cands.reserve(64);
        if (st.to_move == 'W') {
            for (auto& wk_n : generate_all_king_moves(st.wk)) {
                if (wk_n.distance_to(st.bk) > st.wk.distance_to(st.bk)) continue;
                GameState ns(wk_n, st.wq, st.bk, 'B');
                if (is_legal_state(ns) && !is_stalemate(ns)) cands.push_back(ns);
            }
            for (auto& wq_n : generate_all_queen_moves(st.wq)) {
                if (wq_n.distance_to(st.bk) < 2 && wq_n.distance_to(st.wk) > 1) continue;
                // Check if White King blocks the Queen's path
                bool blocked = false;
                if (wq_n.file == st.wq.file) {
                    int start = min(st.wq.rank, wq_n.rank) + 1;
                    int end = max(st.wq.rank, wq_n.rank);
                    for (int r = start; r < end; r++) {
                        if (Position(wq_n.file, r) == st.wk) { blocked = true; break; }
                    }
                } else if (wq_n.rank == st.wq.rank) {
                    int start = min(st.wq.file, wq_n.file) + 1;
                    int end = max(st.wq.file, wq_n.file);
                    for (int f = start; f < end; f++) {
                        if (Position(f, wq_n.rank) == st.wk) { blocked = true; break; }
                    }
                } else if (abs(wq_n.file - st.wq.file) == abs(wq_n.rank - st.wq.rank)) {
                    int df = (wq_n.file > st.wq.file) ? 1 : -1;
                    int dr = (wq_n.rank > st.wq.rank) ? 1 : -1;
                    int f = st.wq.file + df;
                    int r = st.wq.rank + dr;
                    while (f != wq_n.file) {
                        if (Position(f, r) == st.wk) { blocked = true; break; }
                        f += df;
                        r += dr;
                    }
                }
                if (blocked) continue;
                GameState ns(st.wk, wq_n, st.bk, 'B');
                if (is_legal_state(ns) && !is_stalemate(ns)) cands.push_back(ns);
            }
        } else {
            vector<Position> all_bk_moves = generate_all_king_moves(st.bk);
            
            
            
            for (auto& bk_n : all_bk_moves) {
                if (is_attacked_by_queen(bk_n, st.wq, st.wk, st.wq)) {
                    continue;
                }
                
                if (bk_n.distance_to(st.wk) <= 1) {
                    continue;
                }
                
                GameState ns(st.wk, st.wq, bk_n, 'W');
                if (!is_legal_state(ns)) {
                    continue;
                }
                // if (st.to_move == 'B' && ply < 2) {
                //     cout << string(ply*2, ' ') << "[BLACK CAND ply=" << ply << "] bk -> "
                //          << bk_n.str() << " survived filters, state=" << ns.str() << "\n";
                // }
                cands.push_back(ns);
            }
        }
        
        if (cands.empty()) {
            auto res = make_pair(optional<int>(), optional<GameState>());
            memo[cache_key] = res;
            return res;
        }
        // cout << "[ply " << ply << "] cands: ";
        // for (auto& c : cands) cout << c.str() << " | ";
        // cout << "\n";
        candidates_measured += cands.size();
        candidates_measured += cands.size();
        vector<tuple<GameState,int,int>> meas;
        for (auto& c : cands) {
            uint64_t m_key = make_cache_key(c, 3);
            int M;
            
            auto it = M_cache.find(m_key);
            if (it != M_cache.end()) {
                M = it->second;
            } else {
                // HERE A BAD RECURSIVE CALL
                M = compute_M_shallow(c, 3);
                M_cache[m_key] = M;
            }
            
            int BN = count_legal_moves(c);
            meas.push_back({c, M, BN});
        }
        
        int best_M;
        string dir;
        if (st.to_move == 'W') {
            best_M = INF_MAX;
            for (auto& [c,M,BN] : meas) best_M = min(best_M, M);
            dir = "minimize";
        } else {
            best_M = INF_MIN;
            for (auto& [c,M,BN] : meas) best_M = max(best_M, M);
            dir = "maximize";
        }
        
        int thresh = max(2, best_M / 3);
        vector<GameState> viable;
        for (auto& [c,M,BN] : meas) {
            if (dir == "minimize" && M <= best_M + thresh) viable.push_back(c);
            else if (dir == "maximize" && M >= best_M - thresh) viable.push_back(c);
        }
        
        if (viable.size() < 3) {
            sort(meas.begin(), meas.end(), [](auto& a, auto& b) { return get<1>(a) < get<1>(b); });
            viable.clear();
            for (size_t i = 0; i < min(size_t(3), meas.size()); i++) {
                viable.push_back(get<0>(meas[i]));
            }
        }
        
        if (debug && ply < 4) {
            string ind(ply*2, ' ');
            double ratio = cands.empty() ? 0 : (double)viable.size()/cands.size();
            //cout << ind << "[Ply " << ply << "] " << st.to_move << ": " << cands.size()
                 //<< " -> " << viable.size() << " (" << fixed << setprecision(1) << (ratio*100)
                // << "%), best_M=" << best_M << ", thresh=" << thresh << "\n";
        }
        
        optional<int> best_val;
        optional<GameState> best_mv;
        
        for (size_t i = 0; i < viable.size(); i++) {
            auto [val, mv] = compositional_search_impl(viable[i], depth-1, ply+1, debug, memo);
            nodes_evaluated++;
            
                        if (val) {
                int v = *val + 1;
                if (val) {
                    int v = *val + 1;
                    if (!best_val) {
                        best_val = v;
                        best_mv = viable[i];
                    } else if (dir == "minimize" && v < *best_val) {
                        best_val = v;
                        best_mv = viable[i];
                    } else if (dir == "minimize" && v == *best_val) {
                        // ← TIE FOR WHITE: prefer the move that minimizes Black's resulting options
                        int curr_bn = count_legal_moves(viable[i]);
                        int best_bn = count_legal_moves(*best_mv);
                        if (curr_bn < best_bn) {
                            best_mv = viable[i];
                        }
                    } else if (dir == "maximize" && v > *best_val) {
                        best_val = v;
                        best_mv = viable[i];
                    } else if (dir == "maximize" && v == *best_val) {
                        // ← TIE FOR BLACK: use topological M as tiebreaker
                        // HERE A COUPLE OF BAD ESTIMATES
                        int curr_topo = compute_M_topological(viable[i]);
                        int best_topo = compute_M_topological(*best_mv);
                        if (curr_topo > best_topo) {
                            best_mv = viable[i];
                        }
                    }
                }
            }
        }
        
        auto res = make_pair(best_val, best_mv);
        memo[cache_key] = res;
        return res;
    }
    
    pair<optional<GameState>, optional<int>> find_best_move(
        const GameState& st, int max_depth = 10, bool debug = false
    ) {
        nodes_evaluated = 0;
        candidates_measured = 0;
        
        unordered_map<uint64_t, pair<optional<int>, optional<GameState>>> memo;
        
        optional<GameState> best_move;
        optional<int> best_value;
        
        for (int depth = 2; depth <= max_depth + 1; depth += 2) {
            
            auto [md, bm] = compositional_search_impl(st, depth, 0, debug, memo);
            
            if (md) {
                if (!best_value) {
                    best_value = md;
                    best_move = bm;
                } else if (st.to_move == 'W' && md < *best_value) {
                    best_value = md;
                    best_move = bm;
                } else if (st.to_move == 'B' && md > *best_value) {
                    best_value = md;
                    best_move = bm;
                }
            }}
        return make_pair(best_move, best_value);
    }
    
    tuple<vector<string>, int, vector<int>> play_complete_game(
        const GameState& first, SolvedPositionDatabase& db,
        int max_moves = 50, bool debug = false, int game_search_depth = 8,
        vector<string> initial_moves = {}, vector<int> initial_bnc = {}
    ) {
        vector<string> mvs = initial_moves;
        GameState curr = first;
        vector<int> bnc = initial_bnc;
        
        for (int move_num = 0; move_num < max_moves; move_num++) {
            if (is_checkmate(curr)) {
                return make_tuple(mvs, (int)mvs.size(), bnc);
            }
            
            auto [ns, md] = find_best_move(curr, 2*game_search_depth, debug);
            
            if (!ns) {
                return make_tuple(mvs, (int)mvs.size(), bnc);
            }
            
            // CHECK DB HERE - if next position already solved, stop
            if (db.is_solved(ns->str(), ns->to_move)) {
                auto cached = db.get_solution(ns->str(), ns->to_move);
                if (cached) {
                    string mv_str = get_move_notation(curr, *ns);
                    mvs.push_back(mv_str);
                    int tp = mvs.size() + cached->total_plies;  // Add remaining plies from cache
                    cout << "DB FINISHED COMPLETE GAME FROM SOLVED POINT\n";
                    return make_tuple(mvs, tp, bnc);
                }
            }

            string mv_str = get_move_notation(curr, *ns);
            mvs.push_back(mv_str);
            
            int bn = count_legal_moves(*ns);
            bnc.push_back(bn);
            
            curr = *ns;
        }
        
        return make_tuple(mvs, (int)mvs.size(), bnc);
    }
};

// ============================================================================ 
// batch mode
// ============================================================================
 
// Helper function to parse position from string
Position parse_position_from_string(const string& pos_str) {
    // Input: "a1" or "h8"
    char file_char = pos_str[0];
    char rank_char = pos_str[1];
    int file = file_char - 'a';
    int rank = rank_char - '1';
    return Position(file, rank);
}

// Load positions from the Python-generated file
vector<GameState> load_positions_from_file(const string& filename) {
    vector<GameState> positions;
    ifstream file(filename);
    
    if (!file.is_open()) {
        cerr << "ERROR: Cannot open file: " << filename << "\n";
        return positions;
    }
    
    string line;
    bool is_header = true;
    int line_num = 0;
    
    while (getline(file, line)) {
        line_num++;
        
        // Skip header
        if (is_header) {
            is_header = false;
            continue;
        }
        
        if (line.empty()) continue;
        
        // Parse line: "DTZ,WK:a1 WQ:b2 BK:c3"
        size_t comma_pos = line.find(',');
        if (comma_pos == string::npos) {
            cerr << "  ✗ Line " << line_num << " has no comma\n";
            continue;
        }
        
        string position_str = line.substr(comma_pos + 1);
        
        // Parse position string: "WK:a1 WQ:b2 BK:c3"
        size_t wk_pos = position_str.find("WK:") + 3;
        size_t wq_pos = position_str.find("WQ:") + 3;
        size_t bk_pos = position_str.find("BK:") + 3;
        
        try {
            Position wk = parse_position_from_string(position_str.substr(wk_pos, 2));
            Position wq = parse_position_from_string(position_str.substr(wq_pos, 2));
            Position bk = parse_position_from_string(position_str.substr(bk_pos, 2));
            
            GameState st(wk, wq, bk, 'W');
            positions.push_back(st);
        } catch (const exception& e) {
            cerr << "  ✗ Line " << line_num << " parse error: " << e.what() << "\n";
        }
    }
    
    file.close();
    cout << "  ✓ Loaded " << positions.size() << " positions from " << filename << "\n";
    return positions;
}

void batch_solve_all_kqvk_positions(CompositionalEngine& eng, SolvedPositionDatabase& db, int max_depth = 16) {
    cout << "\n" << string(80, '=') << "\n";
    cout << "BATCH SOLVER: ALL KQvK POSITIONS\n";
    cout << string(80, '=') << "\n\n";
    
    // Load positions from file (COMMENT OUT THE GENERATION BELOW)
    cout << "Loading positions from file...\n";
    vector<GameState> positions = load_positions_from_file("kqvk_positions_by_dtz.txt");
    
    if (positions.empty()) {
        cerr << "ERROR: No positions loaded!\n";
        return;
    }
    
    cout << "Loaded " << positions.size() << " legal positions from file\n\n";
    
    // Batch solve
    int solved_count = 0;
    int cache_hit_count = 0;
    auto batch_start = chrono::high_resolution_clock::now();

    for (size_t idx = 0; idx < positions.size(); idx++) {
        // BELOW IS ROOT CAUSE FOR FIRST HORRIBLE MOVE -- replaced with non-batch's exhaustive methodology
        const GameState& pos = positions[idx];
        // HERE A BAD ESTIMATE
        int m_est = eng.compute_M_topological(pos);
        
        // Check database first
        if (db.is_solved(pos.str(), pos.to_move)) {
            cache_hit_count++;
            if (idx % 100 == 0) {
                cout << "[" << idx << "/" << positions.size() << "] [CACHE] " 
                     << pos.str() << " M~" << m_est << "\n";
            }
            continue;
        }

        if (pos.to_move != 'W') {
            cout << "[" << idx << "/" << positions.size() << "] [FAILED] " 
                 << pos.str() << " (batch currently only handles White-to-move root positions)\n";
            continue;
        }

        int opt_d = -1;
        vector<tuple<GameState, int, int>> opt_meas;

        for (int td = 1; td <= 10; td++) {
            vector<GameState> root_cands;
            for (auto& wk_n : eng.generate_all_king_moves(pos.wk)) {
                if (wk_n.distance_to(pos.bk) > pos.wk.distance_to(pos.bk)) continue;
                GameState ns(wk_n, pos.wq, pos.bk, 'B');
                if (eng.is_legal_state(ns) && !eng.is_stalemate(ns)) root_cands.push_back(ns);
            }
            for (auto& wq_n : eng.generate_all_queen_moves(pos.wq)) {
                if (wq_n.distance_to(pos.bk) < 2 && wq_n.distance_to(pos.wk) > 1) continue;

                bool blocked = false;
                if (wq_n.file == pos.wq.file) {
                    int start = min(pos.wq.rank, wq_n.rank) + 1;
                    int end = max(pos.wq.rank, wq_n.rank);
                    for (int r = start; r < end; r++) {
                        if (Position(wq_n.file, r) == pos.wk) { blocked = true; break; }
                    }
                } else if (wq_n.rank == pos.wq.rank) {
                    int start = min(pos.wq.file, wq_n.file) + 1;
                    int end = max(pos.wq.file, wq_n.file);
                    for (int f = start; f < end; f++) {
                        if (Position(f, wq_n.rank) == pos.wk) { blocked = true; break; }
                    }
                } else if (abs(wq_n.file - pos.wq.file) == abs(wq_n.rank - pos.wq.rank)) {
                    int df = (wq_n.file > pos.wq.file) ? 1 : -1;
                    int dr = (wq_n.rank > pos.wq.rank) ? 1 : -1;
                    int f = pos.wq.file + df;
                    int r = pos.wq.rank + dr;
                    while (f != wq_n.file) {
                        if (Position(f, r) == pos.wk) { blocked = true; break; }
                        f += df;
                        r += dr;
                    }
                }
                if (blocked) continue;

                GameState ns(pos.wk, wq_n, pos.bk, 'B');
                if (eng.is_legal_state(ns) && !eng.is_stalemate(ns)) root_cands.push_back(ns);
            }

            if (root_cands.empty()) break;

            vector<tuple<GameState,int,int>> meas;
            for (auto& c : root_cands) {
                // HERE we have the recursive call.
                int M = eng.compute_M_shallow(c, td);
                int BN = eng.measure_black_nodes_after_trajectory(c, td);
                meas.push_back({c, M, BN});
            }

            int min_M = INF_MAX;
            for (auto& [c,M,BN] : meas) min_M = min(min_M, M);

            if (td >= min_M) {
                opt_d = td;
                opt_meas = meas;
                break;
            }
        }

        if (opt_d == -1) {
            opt_d = 7;
            vector<GameState> root_cands;
            for (auto& wk_n : eng.generate_all_king_moves(pos.wk)) {
                if (wk_n.distance_to(pos.bk) > pos.wk.distance_to(pos.bk)) continue;
                GameState ns(wk_n, pos.wq, pos.bk, 'B');
                if (eng.is_legal_state(ns) && !eng.is_stalemate(ns)) root_cands.push_back(ns);
            }
            for (auto& wq_n : eng.generate_all_queen_moves(pos.wq)) {
                if (wq_n.distance_to(pos.bk) < 2 && wq_n.distance_to(pos.wk) > 1) continue;

                bool blocked = false;
                if (wq_n.file == pos.wq.file) {
                    int start = min(pos.wq.rank, wq_n.rank) + 1;
                    int end = max(pos.wq.rank, wq_n.rank);
                    for (int r = start; r < end; r++) {
                        if (Position(wq_n.file, r) == pos.wk) { blocked = true; break; }
                    }
                } else if (wq_n.rank == pos.wq.rank) {
                    int start = min(pos.wq.file, wq_n.file) + 1;
                    int end = max(pos.wq.file, wq_n.file);
                    for (int f = start; f < end; f++) {
                        if (Position(f, wq_n.rank) == pos.wk) { blocked = true; break; }
                    }
                } else if (abs(wq_n.file - pos.wq.file) == abs(wq_n.rank - pos.wq.rank)) {
                    int df = (wq_n.file > pos.wq.file) ? 1 : -1;
                    int dr = (wq_n.rank > pos.wq.rank) ? 1 : -1;
                    int f = pos.wq.file + df;
                    int r = pos.wq.rank + dr;
                    while (f != wq_n.file) {
                        if (Position(f, r) == pos.wk) { blocked = true; break; }
                        f += df;
                        r += dr;
                    }
                }
                if (blocked) continue;

                GameState ns(pos.wk, wq_n, pos.bk, 'B');
                if (eng.is_legal_state(ns) && !eng.is_stalemate(ns)) root_cands.push_back(ns);
            }
            for (auto& c : root_cands) {
                // HERE WE HAVE A RECURSIVE CALL
                int M = eng.compute_M_shallow(c, opt_d);
                int BN = eng.measure_black_nodes_after_trajectory(c, opt_d);
                opt_meas.push_back({c, M, BN});
            }
        }

        if (opt_meas.empty()) {
            cout << "[" << idx << "/" << positions.size() << "] [FAILED] " 
                 << pos.str() << " (no legal moves)\n";
            continue;
        }

        sort(opt_meas.begin(), opt_meas.end(), [](auto& a, auto& b) {
            if (get<1>(a) != get<1>(b)) return get<1>(a) < get<1>(b);
            return get<2>(a) < get<2>(b);
        });

        int min_M_v = INF_MAX;
        for (auto& [c,M,BN] : opt_meas) min_M_v = min(min_M_v, M);

        vector<tuple<GameState,int,int>> opt_cands;
        for (auto& [c,M,BN] : opt_meas) {
            if (M == min_M_v) opt_cands.push_back({c, M, BN});
        }

        // Solve and record EVERY tied-optimal candidate -- exhaustive, not just the first
        for (size_t cand_idx = 0; cand_idx < opt_cands.size(); cand_idx++) {
            auto [best_cand, best_M, best_BN] = opt_cands[cand_idx];

            if (best_M == 0 || eng.is_checkmate(best_cand)) {
                SolvedPosition solution;
                solution.position_key = pos.str();
                solution.turn = pos.to_move;
                solution.best_move = eng.get_move_notation(pos, best_cand);
                solution.M_value = 1;
                solution.total_plies = 1;
                solution.white_moves = 1;
                solution.black_moves = 0;
                solution.nodes_evaluated = 0;
                solution.computation_time = 0;
                solution.BN_trajectory = {};

                db.add_position(solution);
                solved_count++;

                cout << "[" << idx << "/" << positions.size() << "] [SOLVED-MATE] " 
                     << pos.str() << " (candidate " << (cand_idx+1) << "/" << opt_cands.size() << ")\n";
            } else {
                auto solve_start = chrono::high_resolution_clock::now();
                auto [mvs, tp, bnc] = eng.play_complete_game(best_cand, db, 50, false, best_M);
                auto solve_end = chrono::high_resolution_clock::now();
                double solve_time = chrono::duration<double>(solve_end - solve_start).count();

                if (!mvs.empty()) {
                    SolvedPosition solution;
                    solution.position_key = pos.str();
                    solution.turn = pos.to_move;
                    solution.best_move = eng.get_move_notation(pos, best_cand);
                    solution.M_value = best_M;
                    solution.total_plies = tp;
                    solution.white_moves = (tp + 1) / 2;
                    solution.black_moves = tp / 2;
                    solution.nodes_evaluated = 0;
                    solution.computation_time = solve_time;
                    solution.BN_trajectory = bnc;

                    db.add_position(solution);
                    solved_count++;

                    GameState curr = best_cand;
                    for (size_t i = 0; i < mvs.size(); i++) {
                        SolvedPosition ply_solution;
                        ply_solution.position_key = curr.str();
                        ply_solution.turn = curr.to_move;
                        ply_solution.best_move = mvs[i];
                        ply_solution.M_value = best_M;
                        ply_solution.total_plies = tp;
                        ply_solution.white_moves = (tp + 1) / 2;
                        ply_solution.black_moves = tp / 2;
                        ply_solution.nodes_evaluated = 0;
                        ply_solution.computation_time = solve_time;
                        ply_solution.BN_trajectory = bnc;

                        db.add_position(ply_solution);
                        solved_count++;

                        char piece = mvs[i][0];
                        string dest = mvs[i].substr(1);
                        Position dest_pos = Position::from_str(dest);
                        if (piece == 'K') { curr.wk = dest_pos; curr.to_move = 'B'; }
                        else if (piece == 'Q') { curr.wq = dest_pos; curr.to_move = 'B'; }
                        else if (piece == 'k') { curr.bk = dest_pos; curr.to_move = 'W'; }
                    }

                    cout << "[" << idx << "/" << positions.size() << "] [SOLVED] " 
                         << pos.str() << " M=" << best_M
                         << " (candidate " << (cand_idx+1) << "/" << opt_cands.size()
                         << ", " << fixed << setprecision(3) << solve_time << "s)\n";
                } else {
                    cout << "[" << idx << "/" << positions.size() << "] [FAILED] " 
                         << pos.str() << " (candidate " << (cand_idx+1) << "/" << opt_cands.size() << ")\n";
                }
            }
        }
        // end of root cause issue

        // Export every 5 positions
        if (solved_count % 5 == 0 && solved_count > 0) {
            auto export_start = chrono::high_resolution_clock::now();
            db.export_to_file();
            auto export_end = chrono::high_resolution_clock::now();
            double export_time = chrono::duration<double>(export_end - export_start).count();
            
            auto batch_current = chrono::high_resolution_clock::now();
            double elapsed = chrono::duration<double>(batch_current - batch_start).count();
            double rate = solved_count / elapsed;
            
            cout << "\n[CHECKPOINT] Exported " << solved_count << " solutions\n";
            cout << "  Cache hits: " << cache_hit_count << "\n";
            cout << "  Rate: " << fixed << setprecision(1) << rate << " pos/sec\n";
            cout << "  Export time: " << fixed << setprecision(2) << export_time << "s\n";
            cout << "  Elapsed: " << fixed << setprecision(1) << elapsed << "s\n";
            cout << "  ETA: " << fixed << setprecision(1) 
                 << (positions.size() - idx) / rate << "s remaining\n\n";
        }
    }

    // Final export
    db.export_to_file();

    auto batch_end = chrono::high_resolution_clock::now();
    double total_time = chrono::duration<double>(batch_end - batch_start).count();

    cout << "\n" << string(80, '=') << "\n";
    cout << "BATCH SOLVER COMPLETE\n";
    cout << string(80, '=') << "\n";
    cout << "Total positions: " << positions.size() << "\n";
    cout << "Solved: " << solved_count << "\n";
    cout << "Cache hits: " << cache_hit_count << "\n";
    cout << "Total time: " << fixed << setprecision(1) << total_time << "s\n";
    cout << "Rate: " << fixed << setprecision(1) << solved_count / total_time << " pos/sec\n\n";
    }

// ============================================================================
// MAIN
// ============================================================================

int main(int argc, char* argv[]) {
    bool debug_en = false;
    bool batch_mode = false;
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        if (arg == "--debug" || arg == "-d") debug_en = true;
        if (arg == "--batch" || arg == "-b") batch_mode = true;
    }
    
    double root_t = 0;
    
    cout << "\n" << string(80, '=') << "\n";
    cout << "COMPOSITIONAL TOPOLOGICAL SEARCH\n";
    cout << "Measurement with Black Node Count Accumulation\n";
    cout << string(80, '=') << "\n";

    auto root_t_start = chrono::high_resolution_clock::now();
    CompositionalEngine eng;
    SolvedPositionDatabase db("kqvk_perfect_play.db");
    batch_solve_all_kqvk_positions(eng, db, 16);
    return 0;
}
