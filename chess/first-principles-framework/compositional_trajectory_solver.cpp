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
    int cumulative_bn = 0;      // Sum of Black's own escape-count at every Black-to-move
                                 // position from here through to mate, under the recorded
                                 // best_move and every move after it
    
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
           << fixed << setprecision(3) << computation_time << "|"
           << cumulative_bn << "|";
        
        // BN trajectory as comma-separated
        for (size_t i = 0; i < BN_trajectory.size(); i++) {
            if (i > 0) ss << ",";
            ss << BN_trajectory[i];
        }
        ss << "\n";
        return ss.str();
    }
    
    static string csv_header() {
        return "Position|Turn|BestMove|M|Plies|WhiteMoves|BlackMoves|NodesEval|Time|CumulativeBN|BN_Trajectory\n";
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
                    pos.cumulative_bn = (parts.size() > 9 && !parts[9].empty()) ? stoi(parts[9]) : 0;
                    if (parts.size() > 10) {
                        if (!parts[10].empty()) {
                            stringstream bn_stream(parts[10]);
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

struct SearchResult {
    optional<int> val;      // plies to mate -- same meaning as before
    optional<int> bn_cum;   // sum of Black's own escape-square count at every
                             // Black-to-move position along this result's chosen
                             // line, from here through to mate
    optional<GameState> mv; // chosen next state -- same meaning as before
};

class CompositionalEngine : public BaseEngine {
public:
    int nodes_evaluated = 0;
    int candidates_measured = 0;
    unordered_map<uint64_t, int> M_cache;

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

    SearchResult compositional_search_impl(
        const GameState& st, int depth, int ply,
        bool debug, unordered_map<uint64_t, SearchResult>& memo,
        SolvedPositionDatabase& db
    ) {
        uint64_t cache_key = make_cache_key(st, depth);
        if (memo.count(cache_key)) {
            return memo[cache_key];
        }

        if (is_checkmate(st)) {
            SearchResult res{0, 0, nullopt};   // mate: no further black escapes to accumulate
            memo[cache_key] = res;
            return res;
        }

        // DB CHECK ON st ITSELF -- if this exact position is already proven, use it directly.
        // Gated on depth: only trust this shortcut if the cached value would ALSO have been
        // provable by real recursion within the current budget (cached->total_plies <= depth).
        // Otherwise a DB hit could report a "resolved" value at a shallower depth than a
        // fresh-recursion sibling candidate ever gets a fair chance at, letting a worse but
        // cache-lucky branch win a comparison it shouldn't -- see per-candidate gate below
        // for the concrete failure mode this prevents.
        if (db.is_solved(st.str(), st.to_move)) {
            auto cached = db.get_solution(st.str(), st.to_move);
            if (cached && !cached->best_move.empty() && cached->total_plies <= depth) {
                GameState next = st;
                char piece = cached->best_move[0];
                string dest = cached->best_move.substr(1);
                Position dest_pos = Position::from_str(dest);
                if (piece == 'K') { next.wk = dest_pos; next.to_move = 'B'; }
                else if (piece == 'Q') { next.wq = dest_pos; next.to_move = 'B'; }
                else if (piece == 'k') { next.bk = dest_pos; next.to_move = 'W'; }
                SearchResult res{cached->total_plies, cached->cumulative_bn, next};
                memo[cache_key] = res;
                return res;
            }
        }

        if (depth == 0) {
            SearchResult res{nullopt, nullopt, nullopt};
            memo[cache_key] = res;
            return res;
        }

        vector<GameState> cands;
        cands.reserve(64);
        if (st.to_move == 'W') {
            for (auto& wk_n : generate_all_king_moves(st.wk)) {
                GameState ns(wk_n, st.wq, st.bk, 'B');
                if (is_legal_state(ns) && !is_stalemate(ns)) cands.push_back(ns);
            }
            for (auto& wq_n : generate_all_queen_moves(st.wq)) {
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
            for (auto& bk_n : generate_all_king_moves(st.bk)) {
                if (is_attacked_by_queen(bk_n, st.wq, st.wk, st.wq)) continue;
                if (bk_n.distance_to(st.wk) <= 1) continue;
                GameState ns(st.wk, st.wq, bk_n, 'W');
                if (!is_legal_state(ns)) continue;
                cands.push_back(ns);
            }
        }

        if (cands.empty()) {
            SearchResult res{nullopt, nullopt, nullopt};
            memo[cache_key] = res;
            return res;
        }
        candidates_measured += cands.size();

        string dir = (st.to_move == 'W') ? "minimize" : "maximize";

        optional<int> best_val;
        optional<int> best_bn_cum;
        optional<GameState> best_mv;

        for (auto& c : cands) {
            optional<int> val;
            optional<int> child_bn_cum;

            // Gated the same way as the st-itself check above: a cached value is only
            // usable here if it fits within the budget a fresh recursive call would have
            // (cached->total_plies < depth, since v = cached->total_plies + 1 must be <= depth).
            // If it doesn't fit yet, fall through to real recursion below -- which, for
            // consistency, will hit the same gate on its own "st itself" check and correctly
            // decline to use it too, forcing genuine step-by-step exploration until a later,
            // deeper find_best_move iteration naturally catches up to what the cache already
            // knows. Without this gate, a cache hit can report a "resolved" value at a much
            // shallower depth than a genuinely-better sibling candidate needs to prove itself,
            // letting a worse but cache-lucky move win a comparison it has no right to win.
            if (db.is_solved(c.str(), c.to_move)) {
                auto cached = db.get_solution(c.str(), c.to_move);
                if (cached && cached->total_plies < depth) {
                    val = cached->total_plies;
                    child_bn_cum = cached->cumulative_bn;
                }
            }

            if (!val) {
                SearchResult rec = compositional_search_impl(c, depth-1, ply+1, debug, memo, db);
                nodes_evaluated++;
                val = rec.val;
                child_bn_cum = rec.bn_cum;
            }

            if (val) {
                int v = *val + 1;
                // c's own contribution to the black-node sum: Black's escape count AT c,
                // counted only if c is actually a Black-to-move decision point.
                int own_contribution = 0;
                if (c.to_move == 'B') {
                    own_contribution = count_legal_moves(c);
                }
                optional<int> this_bn_cum;
                if (child_bn_cum) {
                    this_bn_cum = own_contribution + *child_bn_cum;
                }
                // if child_bn_cum is nullopt (cache-shortcut gap), this_bn_cum stays nullopt too

                if (!best_val) {
                    best_val = v;
                    best_bn_cum = this_bn_cum;
                    best_mv = c;
                } else if (dir == "minimize" && v < *best_val) {
                    best_val = v;
                    best_bn_cum = this_bn_cum;
                    best_mv = c;
                } else if (dir == "minimize" && v == *best_val) {
                    // TIE FOR WHITE: prefer the move minimizing Black's TOTAL cumulative
                    // node count across the entire remaining game -- only decidable when
                    // both candidates actually have a known cumulative value.
                    if (this_bn_cum && best_bn_cum && *this_bn_cum < *best_bn_cum) {
                        best_val = v;
                        best_bn_cum = this_bn_cum;
                        best_mv = c;
                    }
                } else if (dir == "maximize" && v > *best_val) {
                    best_val = v;
                    best_bn_cum = this_bn_cum;
                    best_mv = c;
                } else if (dir == "maximize" && v == *best_val) {
                    // TIE FOR BLACK: prefer the move maximizing Black's OWN TOTAL cumulative
                    // node count across the entire remaining game.
                    if (this_bn_cum && best_bn_cum && *this_bn_cum > *best_bn_cum) {
                        best_val = v;
                        best_bn_cum = this_bn_cum;
                        best_mv = c;
                    }
                }
            }
        }

        SearchResult res{best_val, best_bn_cum, best_mv};
        memo[cache_key] = res;
        return res;
    }

    pair<optional<GameState>, optional<int>> find_best_move(
        const GameState& st, SolvedPositionDatabase& db, int max_depth = 10, bool debug = false
    ) {
        nodes_evaluated = 0;
        candidates_measured = 0;

        unordered_map<uint64_t, SearchResult> memo;

        optional<GameState> best_move;
        optional<int> best_value;

        for (int depth = 2; depth <= max_depth + 1; depth += 2) {
            SearchResult r = compositional_search_impl(st, depth, 0, debug, memo, db);
            if (r.val) {
                if (!best_value) {
                    best_value = r.val;
                    best_move = r.mv;
                } else if (st.to_move == 'W' && r.val < *best_value) {
                    best_value = r.val;
                    best_move = r.mv;
                } else if (st.to_move == 'B' && r.val > *best_value) {
                    best_value = r.val;
                    best_move = r.mv;
                }
                break;
            }
        }
        return make_pair(best_move, best_value);
    }
    
    tuple<vector<string>, int, vector<int>, bool> play_complete_game(
        const GameState& first, SolvedPositionDatabase& db,
        int max_moves = 50, bool debug = false, int game_search_depth = 8,
        vector<string> initial_moves = {}, vector<int> initial_bnc = {}
    ) {
        vector<string> mvs = initial_moves;
        GameState curr = first;
        vector<int> bnc = initial_bnc;
        
        for (int move_num = 0; move_num < max_moves; move_num++) {
            if (is_checkmate(curr)) {
                return make_tuple(mvs, (int)mvs.size(), bnc, true);
            }
            
            auto [ns, md] = find_best_move(curr, db, 2*game_search_depth, debug);
            
            if (!ns) {
                return make_tuple(mvs, (int)mvs.size(), bnc, false);
            }
            
            // CHECK DB HERE - if next position already solved, stop
            if (db.is_solved(ns->str(), ns->to_move)) {
                auto cached = db.get_solution(ns->str(), ns->to_move);
                if (cached) {
                    string mv_str = get_move_notation(curr, *ns);
                    mvs.push_back(mv_str);
                    int bn = (ns->to_move == 'B') ? count_legal_moves(*ns) : 0;
                    bnc.push_back(bn);
                    int tp = mvs.size() + cached->total_plies;  // Add remaining plies from cache
                    cout << "DB FINISHED COMPLETE GAME FROM SOLVED POINT\n";
                    return make_tuple(mvs, tp, bnc, true);
                }
            }

            string mv_str = get_move_notation(curr, *ns);
            mvs.push_back(mv_str);
            
            int bn = (ns->to_move == 'B') ? count_legal_moves(*ns) : 0;
            bnc.push_back(bn);
            
            curr = *ns;
        }
        
        return make_tuple(mvs, (int)mvs.size(), bnc, false);
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
        
        // Check database first
        if (db.is_solved(pos.str(), pos.to_move)) {
            cache_hit_count++;
            if (idx % 100 == 0) {
                cout << "[" << idx << "/" << positions.size() << "] [CACHE] " 
                     << pos.str() << "\n";
            }
            continue;
        }

        if (pos.to_move != 'W') {
            cout << "[" << idx << "/" << positions.size() << "] [FAILED] " 
                 << pos.str() << " (batch currently only handles White-to-move root positions)\n";
            continue;
        }

        eng.nodes_evaluated = 0;
        eng.candidates_measured = 0;

        auto solve_start = chrono::high_resolution_clock::now();
        auto [mvs, tp, bnc, reached_mate] = eng.play_complete_game(pos, db, 50, false, max_depth);
        auto solve_end = chrono::high_resolution_clock::now();
        double solve_time = chrono::duration<double>(solve_end - solve_start).count();

        if (mvs.empty() || !reached_mate) {
            cout << "[" << idx << "/" << positions.size() << "] [FAILED] " 
                 << pos.str() << " (no proven mate found, mvs=" << mvs.size() << ")\n";
            continue;
        }

        GameState curr = pos;
        int recorded_count = 0;

        for (size_t i = 0; i < mvs.size(); i++) {
            SolvedPosition solution;
            solution.position_key = curr.str();
            solution.turn = curr.to_move;
            solution.best_move = mvs[i];
            solution.white_moves = (tp - i + 1) / 2;
            solution.black_moves = (tp - i) / 2;
            solution.M_value = solution.white_moves;
            solution.total_plies = tp - i;
            solution.nodes_evaluated = eng.nodes_evaluated;
            solution.computation_time = solve_time / mvs.size();

            // bnc[i] is the black-escape-count measured AFTER mvs[i] is played -- i.e. it's
            // exactly the "own_contribution" of curr's own chosen next move, matching how
            // bn_cum is defined recursively in compositional_search_impl. So the trajectory
            // "from here to mate" for the position being recorded THIS iteration is bnc[i:],
            // inclusive of i -- not bnc[i+1:].
            solution.BN_trajectory.assign(bnc.begin() + i, bnc.end());
            solution.cumulative_bn = 0;
            for (int v : solution.BN_trajectory) solution.cumulative_bn += v;

            db.add_position(solution);
            solved_count++;
            recorded_count++;

            char piece = mvs[i][0];
            string dest = mvs[i].substr(1);
            Position dest_pos = Position::from_str(dest);
            if (piece == 'K') { curr.wk = dest_pos; curr.to_move = 'B'; }
            else if (piece == 'Q') { curr.wq = dest_pos; curr.to_move = 'B'; }
            else if (piece == 'k') { curr.bk = dest_pos; curr.to_move = 'W'; }
        }

        cout << "[" << idx << "/" << positions.size() << "] [SOLVED] " 
             << pos.str() << " plies=" << tp << " recorded=" << recorded_count
             << " (" << fixed << setprecision(3) << solve_time << "s)\n";
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
    batch_solve_all_kqvk_positions(eng, db, 25);
    return 0;
}
