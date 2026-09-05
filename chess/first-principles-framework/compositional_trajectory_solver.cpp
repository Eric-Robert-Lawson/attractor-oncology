#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <deque>
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

    // Every move at THIS position that achieves the exact same M_value (i.e. every move
    // provably tied for game-theoretically optimal), each paired with its own cumulative
    // black-escape count. best_move is always one entry of this list -- kept as a separate
    // field for backward compatibility with anything that only reads BestMove/M and ignores
    // this column. Populated by compositional_search_impl's tie-collection pass (see the
    // `tied` field on SearchResult); NOT a re-derivation from best_move, so it is exact,
    // not approximated after the fact.
    vector<pair<string, int>> tied_moves;

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
        ss << "|";

        // Tied moves as "move:bncum" pairs, semicolon-separated (comma is already used
        // above by BN_trajectory, and pipe is the field separator, so semicolon/colon
        // are the only delimiters left unused by the existing format).
        for (size_t i = 0; i < tied_moves.size(); i++) {
            if (i > 0) ss << ";";
            ss << tied_moves[i].first << ":" << tied_moves[i].second;
        }
        ss << "\n";
        return ss.str();
    }
    
    static string csv_header() {
        return "Position|Turn|BestMove|M|Plies|WhiteMoves|BlackMoves|NodesEval|Time|CumulativeBN|BN_Trajectory|TiedMoves\n";
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
            
            // Applies transform_func to a single "K"/"Q"/"k" + destination-square move
            // string. Shared by best_move and by every entry of tied_moves so both stay
            // in exact agreement after any symmetry transform.
            auto transform_move = [&](auto transform_func, const string& mv) {
                string new_mv = mv;
                if (!new_mv.empty() && new_mv != "mate" && new_mv != "checkmate" && new_mv.length() >= 3) {
                    char piece = new_mv[0];
                    string dest = new_mv.substr(1);
                    new_mv = piece + transform_func(dest);
                }
                return new_mv;
            };

            // Helper macro to add transformed position
            auto add_transformed = [&](auto transform_func, const string& base_pos, const string& base_move,
                                        const vector<pair<string,int>>& base_tied) {
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
                
                string new_move = transform_move(transform_func, base_move);

                vector<pair<string,int>> new_tied;
                new_tied.reserve(base_tied.size());
                for (auto& [mv, bncum] : base_tied) {
                    new_tied.emplace_back(transform_move(transform_func, mv), bncum);
                }
                
                if (!is_solved(new_pos, pos.turn)) {
                    SolvedPosition transformed = pos;
                    transformed.position_key = new_pos;
                    transformed.best_move = new_move;
                    transformed.tied_moves = new_tied;
                    solved[{new_pos, pos.turn}] = transformed;
                }
            };
            
            // Generate 4 rotations
            string current_pos = pos.position_key;
            string current_move = pos.best_move;
            vector<pair<string,int>> current_tied = pos.tied_moves;
            
            for (int rot = 1; rot < 4; rot++) {
                add_transformed(rotate_square, current_pos, current_move, current_tied);
                
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
                current_move = transform_move(rotate_square, current_move);

                vector<pair<string,int>> next_tied;
                next_tied.reserve(current_tied.size());
                for (auto& [mv, bncum] : current_tied) {
                    next_tied.emplace_back(transform_move(rotate_square, mv), bncum);
                }
                current_tied = next_tied;
            }
            
            // Generate 4 reflections from original
            add_transformed(flip_vertical, pos.position_key, pos.best_move, pos.tied_moves);
            add_transformed(flip_horizontal, pos.position_key, pos.best_move, pos.tied_moves);
            add_transformed(flip_diagonal_a1h8, pos.position_key, pos.best_move, pos.tied_moves);
            add_transformed(flip_diagonal_a8h1, pos.position_key, pos.best_move, pos.tied_moves);
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

    // Patch ONLY tied_moves onto an already-solved entry -- used when DAG exploration
    // encounters a position solved by an older run (before tied-move tracking existed),
    // where M_value/best_move/etc. are already trustworthy and don't need to be
    // recomputed, but tied_moves is empty and needs to be filled in from a fresh
    // find_best_move call. No-op if the position isn't already solved (callers should
    // use add_position for that case instead). Note: this does NOT propagate to that
    // entry's 8 symmetric siblings -- if they're also legacy rows missing tied_moves,
    // each will independently backfill itself the first time DAG exploration visits it.
    void set_tied_moves(const string& position_key, char turn, const vector<pair<string,int>>& tied) {
        auto it = solved.find({position_key, turn});
        if (it != solved.end()) {
            it->second.tied_moves = tied;
        }
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
                    // TiedMoves: "move:bncum" entries separated by ';'. Older database
                    // files (written before this field existed) simply won't have this
                    // 12th column -- parts.size() <= 11 leaves tied_moves empty, which is
                    // the correct, honest representation of "this row predates tie
                    // tracking" rather than guessing at a fabricated tie set.
                    if (parts.size() > 11 && !parts[11].empty()) {
                        stringstream tie_stream(parts[11]);
                        string tie_part;
                        while (getline(tie_stream, tie_part, ';')) {
                            size_t colon = tie_part.rfind(':');
                            if (colon == string::npos) continue;
                            try {
                                string mv = tie_part.substr(0, colon);
                                int bncum = stoi(tie_part.substr(colon + 1));
                                pos.tied_moves.emplace_back(mv, bncum);
                            } catch (const exception& e) {
                                cout << "    -> ERROR parsing TiedMoves entry: " << e.what() << "\n";
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

    // Inverse of get_move_notation: applies a "K"/"Q"/"k" + destination-square move
    // string to `from`, returning the resulting GameState (including the turn flip).
    // Used by DAG exploration to reconstruct a tied move's child state directly from
    // a stored (or freshly computed) move string, without re-running search on it.
    GameState apply_move_notation(const GameState& from, const string& mv) const {
        GameState result = from;
        if (mv.length() < 2) return result;
        char piece = mv[0];
        Position dest = Position::from_str(mv.substr(1));
        if (piece == 'K') { result.wk = dest; result.to_move = 'B'; }
        else if (piece == 'Q') { result.wq = dest; result.to_move = 'B'; }
        else if (piece == 'k') { result.bk = dest; result.to_move = 'W'; }
        return result;
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

    // Every RESOLVED candidate whose v exactly equals val (the game-theoretic optimum),
    // each with its own bn_cum computed by the same own_contribution + child_bn_cum
    // formula used everywhere else in this engine. This is a pure addition: val/bn_cum/mv
    // above keep exactly their previous single-choice meaning (mv == tied[k].second for
    // whichever k the existing tie-break rule selected), so nothing that already reads
    // val/bn_cum/mv changes behavior. Unresolved candidates are never in this list --
    // see the comment above the resolution loop in compositional_search_impl for why an
    // unresolved candidate can never tie val, for either side to move.
    // Empty for terminal (checkmate) and unresolved (nullopt val) results.
    vector<pair<GameState, int>> tied;
};

class CompositionalEngine : public BaseEngine {
public:
    int nodes_evaluated = 0;
    int candidates_measured = 0;
    unordered_map<uint64_t, int> M_cache;

    uint64_t make_cache_key(const GameState& st, int depth) const {
        // Each board field only ever holds 0-7 (3 bits), given 4 bits of margin here.
        // depth has real range up to ~50 in practice -- given the full remaining 40 bits
        // (shift 0) so it cannot overflow into bk_rank's field the way the previous layout
        // did. Confirmed by direct collision test: the old layout (depth at shift 32, only
        // 4 bits before bk_rank's field at shift 36) made (bk_rank=3, depth=0) produce the
        // exact same key as the unrelated (bk_rank=2, depth=16) -- any depth >= 16 corrupted
        // memo lookups into a different board position entirely.
        uint64_t key = 0;
        key |= ((uint64_t)st.wk.file << 60);
        key |= ((uint64_t)st.wk.rank << 56);
        key |= ((uint64_t)st.wq.file << 52);
        key |= ((uint64_t)st.wq.rank << 48);
        key |= ((uint64_t)st.bk.file << 44);
        key |= ((uint64_t)st.bk.rank << 40);
        key |= (uint64_t)depth;
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

        // DB shortcuts (both "check st itself" and per-candidate, further below) have been
        // REMOVED from this comparison path entirely -- not just re-gated. The gate this
        // engine used (cached->total_plies <= depth) checks whether the ANSWER fits in the
        // current budget, but proving that answer -- under the rule that every one of
        // Black's legal replies must also fully resolve -- can require far more depth than
        // the answer itself. A cached entry whose original proof needed, say, 12 plies of
        // budget can still slip through a "value <= depth" gate at depth 5, getting treated
        // as resolved when a fresh, honest recursion at that same depth would correctly say
        // unresolved. That's not a timing-fairness problem this gate can catch; it's a gap
        // in what the gate checks. Confirmed empirically: WK:d1 WQ:f7 BK:c6 resolves to an
        // impossible 13 plies with this shortcut enabled (proven_cache, in that instance),
        // and to the correct, Syzygy-verified 15 plies with it removed. The database is
        // still safely used elsewhere -- see play_complete_game's own lookup, and the
        // top-level cache check in batch_solve_all_kqvk_positions -- both of which consult
        // it only AFTER a decision is already made, never as a substitute inside one.

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
                // NOT a chess heuristic -- a structural necessity. GameState has no way
                // to represent a captured queen (wq is always some square; there is no
                // "queen is gone" state), so is_legal_state rejects a king-onto-queen
                // capture attempt as "two pieces on one square" rather than resolving it
                // as a capture. That means this engine cannot currently represent or
                // correctly evaluate the aftermath of losing the queen at all. Until that
                // representational gap is actually fixed, generating a queen move the
                // opponent's king could capture for free is not explorable, it's a silent
                // blind spot -- the search would wrongly treat it as safe purely because
                // the refutation is structurally unreachable, not because it's actually safe.
                if (wq_n.distance_to(st.bk) < 2 && wq_n.distance_to(st.wk) > 1) continue;

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
        // Only load-bearing for Black (maximize) -- see the gate at the end of this loop.
        // Left in place unconditionally rather than branching the whole loop on `dir`, since
        // computing it costs nothing extra: it's just tracking whether every candidate resolved.
        bool all_candidates_resolved = true;

        // Every candidate that resolves this ply, regardless of whether it ends up winning
        // the best_val/best_bn_cum comparison below. best_val/best_bn_cum/best_mv's selection
        // logic is completely untouched by this addition -- this is purely bookkeeping so that,
        // once best_val is known, every candidate tied on it (not just the tie-break winner)
        // can be recovered in one filtering pass instead of being discarded as the loop runs.
        vector<pair<int, pair<optional<int>, GameState>>> resolved_candidates;
        resolved_candidates.reserve(cands.size());

        for (auto& c : cands) {
            optional<int> val;
            optional<int> child_bn_cum;

            // No database shortcut here -- see the note above this function for why a
            // "cached value fits within the current depth budget" gate is not sufficient
            // to make this safe. Every candidate is honestly recursed into.
            SearchResult rec = compositional_search_impl(c, depth-1, ply+1, debug, memo, db);
            nodes_evaluated++;
            val = rec.val;
            child_bn_cum = rec.bn_cum;

            if (!val) {
                // An unresolved candidate's true value is provably larger than the current
                // depth budget (real recursion needed more than depth-1 plies to prove it).
                // For White (minimize), that bound alone means it can never beat an already
                // resolved candidate (whose value is provably <= depth) -- safe to ignore,
                // exactly as this code already did before this fix.
                // For Black (maximize), that SAME bound means the opposite: an unresolved
                // candidate's true value is provably AT LEAST AS GOOD as anything resolved
                // so far, and quite possibly better. Silently ignoring it here would let a
                // move that merely resolves quickly (often precisely because it's bad --
                // shorter, easier-to-prove lines are exactly the ones an evading king should
                // be avoiding) beat a genuinely better move that simply needs more depth.
                all_candidates_resolved = false;
                continue;
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

                resolved_candidates.push_back({v, {this_bn_cum, c}});

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

        if (dir == "maximize" && !all_candidates_resolved) {
            // Black's decision cannot be trusted yet: at least one legal candidate is
            // unresolved, and per the bound above, any such candidate is provably at least
            // as good as best_val found so far. Report unproven and let the caller retry at
            // greater depth, exactly like running out of budget with no candidates at all.
            SearchResult res{nullopt, nullopt, nullopt};
            memo[cache_key] = res;
            return res;
        }

        vector<pair<GameState, int>> tied;
        if (best_val) {
            for (auto& [v, rest] : resolved_candidates) {
                if (v != *best_val) continue;
                auto& [bncum_opt, state] = rest;
                // See the inductive argument in the SearchResult struct comment: bncum_opt
                // should always be populated whenever v is; -1 is a defensive fallback only,
                // never expected to actually trigger given the current (post-DB-shortcut-
                // removal) code path.
                int bncum_val = bncum_opt ? *bncum_opt : -1;
                tied.emplace_back(state, bncum_val);
            }
        }

        SearchResult res{best_val, best_bn_cum, best_mv, tied};
        memo[cache_key] = res;
        return res;
    }

    tuple<optional<GameState>, optional<int>, vector<pair<GameState, int>>> find_best_move(
        const GameState& st, SolvedPositionDatabase& db, int max_depth = 10, bool debug = false
    ) {
        nodes_evaluated = 0;
        candidates_measured = 0;

        unordered_map<uint64_t, SearchResult> memo;

        optional<GameState> best_move;
        optional<int> best_value;
        vector<pair<GameState, int>> tied;

        for (int depth = 2; depth <= max_depth + 1; depth += 2) {
            SearchResult r = compositional_search_impl(st, depth, 0, debug, memo, db);
            if (r.val) {
                if (!best_value) {
                    best_value = r.val;
                    best_move = r.mv;
                    tied = r.tied;
                } else if (st.to_move == 'W' && r.val < *best_value) {
                    best_value = r.val;
                    best_move = r.mv;
                    tied = r.tied;
                } else if (st.to_move == 'B' && r.val > *best_value) {
                    best_value = r.val;
                    best_move = r.mv;
                    tied = r.tied;
                }
                break;
            }
        }
        return make_tuple(best_move, best_value, tied);
    }
    
    tuple<vector<string>, int, vector<int>, bool, vector<vector<pair<string,int>>>> play_complete_game(
        const GameState& first, SolvedPositionDatabase& db,
        int max_moves = 50, bool debug = false, int game_search_depth = 10,
        vector<string> initial_moves = {}, vector<int> initial_bnc = {},
        vector<vector<pair<string,int>>> initial_tied = {}
    ) {
        vector<string> mvs = initial_moves;
        GameState curr = first;
        vector<int> bnc = initial_bnc;
        // tied_per_ply[i] is the FULL set of moves (as "notation:bncum") that tied for
        // optimal AT the position mvs[i] was played from -- i.e. every move provably as
        // good as mvs[i] itself, not just the one this game trajectory actually took.
        vector<vector<pair<string,int>>> tied_per_ply = initial_tied;
        
        for (int move_num = 0; move_num < max_moves; move_num++) {
            if (is_checkmate(curr)) {
                return make_tuple(mvs, (int)mvs.size(), bnc, true, tied_per_ply);
            }
            
            auto [ns, md, tied] = find_best_move(curr, db, 2*game_search_depth, debug);
            
            if (!ns) {
                return make_tuple(mvs, (int)mvs.size(), bnc, false, tied_per_ply);
            }

            vector<pair<string,int>> tied_notation;
            tied_notation.reserve(tied.size());
            for (auto& [tstate, tbncum] : tied) {
                tied_notation.emplace_back(get_move_notation(curr, tstate), tbncum);
            }
            
            // CHECK DB HERE - if next position already solved, stop
            if (db.is_solved(ns->str(), ns->to_move)) {
                auto cached = db.get_solution(ns->str(), ns->to_move);
                if (cached) {
                    string mv_str = get_move_notation(curr, *ns);
                    mvs.push_back(mv_str);
                    int bn = (ns->to_move == 'B') ? count_legal_moves(*ns) : 0;
                    bnc.push_back(bn);
                    tied_per_ply.push_back(tied_notation);
                    int tp = mvs.size() + cached->total_plies;  // Add remaining plies from cache
                    cout << "DB FINISHED COMPLETE GAME FROM SOLVED POINT\n";
                    return make_tuple(mvs, tp, bnc, true, tied_per_ply);
                }
            }

            string mv_str = get_move_notation(curr, *ns);
            mvs.push_back(mv_str);
            
            int bn = (ns->to_move == 'B') ? count_legal_moves(*ns) : 0;
            bnc.push_back(bn);
            tied_per_ply.push_back(tied_notation);
            
            curr = *ns;
        }
        
        return make_tuple(mvs, (int)mvs.size(), bnc, false, tied_per_ply);
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

// ============================================================================
// Full attractor-DAG exploration
// ============================================================================
//
// Packs (WK,WQ,BK,turn) into a node identity key -- deliberately NOT keyed on
// search depth (unlike CompositionalEngine::make_cache_key), because this is
// used for a single "have I already fully processed this NODE at all, ever,
// in this run" visited-set, not a per-depth memo. 3 bits per coordinate x 6
// coordinates + 1 bit for turn = 19 bits; comfortably fits uint64_t.
uint64_t pack_state_key(const GameState& st) {
    uint64_t key = 0;
    key |= ((uint64_t)st.wk.file << 20);
    key |= ((uint64_t)st.wk.rank << 17);
    key |= ((uint64_t)st.wq.file << 14);
    key |= ((uint64_t)st.wq.rank << 11);
    key |= ((uint64_t)st.bk.file << 8);
    key |= ((uint64_t)st.bk.rank << 5);
    key |= (st.to_move == 'W' ? 1ull : 0ull) << 4;
    return key;
}

// Explores every position reachable from `roots` by following ONLY tied-optimal
// moves at every step -- the actual perfect-play attractor DAG, not a single
// canonical line through it. This is the direct generalization of
// batch_solve_all_kqvk_positions (which records one game-length line per root)
// to record every branch that is exactly as good as that line.
//
// The combinatorial-explosion concern that applies to PATH enumeration (tied
// branches compounding multiplicatively ply over ply) does NOT apply here,
// because this is a graph/BFS traversal with a visited-node set, not a tree
// walk: any branch that reconverges onto an already-visited node is cut
// immediately via `visited_this_run` (in-memory, this run) and via
// `db.is_solved` (on disk, persists across runs/resumptions). Total work is
// therefore bounded by the number of DISTINCT REACHABLE NODES, which for
// KQvK is bounded by the full legal-position count (roughly 64*63*62*2 minus
// illegal-adjacency exclusions -- on the order of a few hundred thousand,
// matching the ~178,856 rows already in a completed single-path table), not
// by the number of paths between them.
//
// Known, deliberate simplifications versus batch_solve_all_kqvk_positions:
//   - BN_trajectory (the full per-ply array) is left empty. It was defined as
//     "the sequence of black-mobility counts along ONE specific forward line",
//     which doesn't have a single canonical meaning for a DAG node that may
//     have several equally-optimal continuations. cumulative_bn (the scalar,
//     taken directly from find_best_move's own recursive bn_cum) IS populated
//     correctly per node and per tied move, and is the field the tie-ranking
//     analysis actually needs.
//   - computation_time/nodes_evaluated are now genuinely PER-NODE (timed
//     around each individual find_best_move call), which is more precise than
//     the old single-path code's behavior of stamping one aggregate number
//     from the whole game onto every position it visited.
void explore_full_attractor_dag(
    CompositionalEngine& eng, SolvedPositionDatabase& db,
    const vector<GameState>& roots, int max_depth = 16,
    int checkpoint_every = 2000
) {
    cout << "\n" << string(80, '=') << "\n";
    cout << "FULL ATTRACTOR-DAG EXPLORATION\n";
    cout << string(80, '=') << "\n\n";

    unordered_set<uint64_t> visited_this_run;
    deque<GameState> frontier;

    for (auto& r : roots) {
        if (!eng.is_checkmate(r)) frontier.push_back(r);
    }
    cout << "Seeded frontier with " << frontier.size() << " root positions\n\n";

    long long processed = 0, newly_solved = 0, cache_reused = 0, backfilled = 0, failed = 0;
    auto run_start = chrono::high_resolution_clock::now();

    while (!frontier.empty()) {
        GameState st = frontier.front();
        frontier.pop_front();

        uint64_t key = pack_state_key(st);
        if (visited_this_run.count(key)) continue;
        visited_this_run.insert(key);

        if (eng.is_checkmate(st)) continue;  // terminal: nothing to expand

        vector<pair<string,int>> tied_notation;

        if (db.is_solved(st.str(), st.to_move)) {
            auto cached = db.get_solution(st.str(), st.to_move);
            if (cached && !cached->tied_moves.empty()) {
                tied_notation = cached->tied_moves;
                cache_reused++;
            } else {
                // Legacy row (solved before tied-move tracking existed): M_value/
                // best_move are already trustworthy, but we need a fresh search at
                // THIS node specifically to recover its tied set.
                auto search_start = chrono::high_resolution_clock::now();
                auto [ns, fval, tied] = eng.find_best_move(st, db, 2*max_depth, false);
                auto search_end = chrono::high_resolution_clock::now();
                (void)search_start; (void)search_end;
                if (!fval) { failed++; continue; }
                for (auto& [tstate, tbncum] : tied) {
                    tied_notation.emplace_back(eng.get_move_notation(st, tstate), tbncum);
                }
                db.set_tied_moves(st.str(), st.to_move, tied_notation);
                backfilled++;
            }
            processed++;
        } else {
            auto search_start = chrono::high_resolution_clock::now();
            auto [ns, fval, tied] = eng.find_best_move(st, db, 2*max_depth, false);
            auto search_end = chrono::high_resolution_clock::now();
            double search_time = chrono::duration<double>(search_end - search_start).count();

            if (!fval) {
                failed++;
                cout << "  [FAILED] " << st.str() << " (to_move=" << st.to_move
                     << ") -- unresolved within depth budget " << (2*max_depth) << "\n";
                continue;
            }

            for (auto& [tstate, tbncum] : tied) {
                tied_notation.emplace_back(eng.get_move_notation(st, tstate), tbncum);
            }

            int val = *fval;
            SolvedPosition solution;
            solution.position_key = st.str();
            solution.turn = st.to_move;
            solution.best_move = tied_notation.empty() ? "" : tied_notation[0].first;
            solution.white_moves = (val + 1) / 2;
            solution.black_moves = val / 2;
            solution.M_value = solution.white_moves;
            solution.total_plies = val;
            solution.nodes_evaluated = eng.nodes_evaluated;
            solution.computation_time = search_time;
            solution.cumulative_bn = tied_notation.empty() ? 0 : tied_notation[0].second;
            solution.tied_moves = tied_notation;

            db.add_position(solution);
            newly_solved++;
            processed++;
        }

        for (auto& [mv_str, bncum] : tied_notation) {
            (void)bncum;
            GameState child = eng.apply_move_notation(st, mv_str);
            uint64_t ckey = pack_state_key(child);
            if (!visited_this_run.count(ckey)) {
                frontier.push_back(child);
            }
        }

        if (processed % checkpoint_every == 0) {
            db.export_to_file();
            double elapsed = chrono::duration<double>(
                chrono::high_resolution_clock::now() - run_start).count();
            cout << "[CHECKPOINT] processed=" << processed
                 << " newly_solved=" << newly_solved
                 << " cache_reused=" << cache_reused
                 << " backfilled=" << backfilled
                 << " failed=" << failed
                 << " frontier=" << frontier.size()
                 << " visited=" << visited_this_run.size()
                 << " elapsed=" << fixed << setprecision(1) << elapsed << "s\n";
        }
    }

    db.export_to_file();
    double total_time = chrono::duration<double>(
        chrono::high_resolution_clock::now() - run_start).count();

    cout << "\n" << string(80, '=') << "\n";
    cout << "FULL ATTRACTOR-DAG EXPLORATION COMPLETE\n";
    cout << string(80, '=') << "\n";
    cout << "Distinct nodes visited: " << visited_this_run.size() << "\n";
    cout << "Newly solved this run:  " << newly_solved << "\n";
    cout << "Reused from DB cache:   " << cache_reused << "\n";
    cout << "Backfilled (legacy):    " << backfilled << "\n";
    cout << "Failed to resolve:      " << failed << "\n";
    cout << "Total time:             " << fixed << setprecision(1) << total_time << "s\n\n";
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
        auto [mvs, tp, bnc, reached_mate, tied_per_ply] = eng.play_complete_game(pos, db, 50, false, max_depth);
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

            if (i < tied_per_ply.size()) {
                solution.tied_moves = tied_per_ply[i];
            }

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
    bool full_dag_mode = false;
    string positions_file = "kqvk_positions_by_dtz.txt";
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        if (arg == "--debug" || arg == "-d") debug_en = true;
        if (arg == "--batch" || arg == "-b") batch_mode = true;
        if (arg == "--full-dag" || arg == "-g") full_dag_mode = true;
        if (arg == "--positions" && i + 1 < argc) positions_file = argv[++i];
    }
    
    double root_t = 0;
    
    cout << "\n" << string(80, '=') << "\n";
    cout << "COMPOSITIONAL TOPOLOGICAL SEARCH\n";
    cout << "Measurement with Black Node Count Accumulation\n";
    cout << string(80, '=') << "\n";

    auto root_t_start = chrono::high_resolution_clock::now();
    CompositionalEngine eng;
    SolvedPositionDatabase db("kqvk_perfect_play.db");

    if (full_dag_mode) {
        // Full perfect-play attractor DAG: every branch tied for optimal, not just
        // one canonical line per root. See explore_full_attractor_dag's comment
        // block for why this doesn't combinatorially explode.
        vector<GameState> roots = load_positions_from_file(positions_file);
        if (roots.empty()) {
            cerr << "ERROR: No root positions loaded from " << positions_file << "\n";
            return 1;
        }
        explore_full_attractor_dag(eng, db, roots, 25);
    } else {
        batch_solve_all_kqvk_positions(eng, db, 25);
    }
    return 0;
}
