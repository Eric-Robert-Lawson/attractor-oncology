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
    string best_move;           // "Q:e5" or "K:c2"
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
        return "Position|BestMove|M|Plies|WhiteMoves|BlackMoves|NodesEval|Time|BN_Trajectory\n";
    }
};

// ============================================================================
// SolvedPositionDatabase - Manages the tablebase
// ============================================================================

class SolvedPositionDatabase {
private:
    map<string, SolvedPosition> solved;
    string filename;
    
public:
    SolvedPositionDatabase(const string& db_file = "kqvk_perfect_play.db")
        : filename(db_file) {
        load_from_file();
    }
    
    // Add a solved position to the database
    void add_position(const SolvedPosition& pos) {
        if (is_solved(pos.position_key)) {
            cout << "[Database] Position already exists, skipping: " << pos.position_key << "\n";
            return;
        }
        solved[pos.position_key] = pos;
        cout << "[Database] Added new position: " << pos.position_key << "\n";
    }
    
    // Check if position is already solved
    bool is_solved(const string& position_key) const {
        return solved.count(position_key) > 0;
    }
    
    // Get optimal move for a position
    optional<SolvedPosition> get_solution(const string& position_key) const {
        auto it = solved.find(position_key);
        if (it != solved.end()) {
            return it->second;
        }
        return nullopt;
    }
    
    // Export all solved positions to file
    void export_to_file() {
        ofstream file(filename, ios::trunc);  // ← Truncate mode (overwrite)
        
        if (!file.is_open()) {
            cerr << "ERROR: Cannot open database file: " << filename << "\n";
            return;
        }
        
        // Always write header
        file << SolvedPosition::csv_header();
        
        // Write all positions
        for (const auto& [key, pos] : solved) {
            file << pos.to_csv();
        }
        
        file.close();
    }
    
    // Load from file for future runs
    void load_from_file() {
        ifstream file(filename);
        
        if (!file.is_open()) {
            cout << "[Database] No existing database found. Starting fresh.\n";
            return;
        }
        
        string line;
        bool is_header = true;
        
        while (getline(file, line)) {
            if (is_header) {
                is_header = false;
                continue;  // Skip header
            }
            
            if (line.empty()) continue;
            
            // Parse CSV line
            vector<string> parts;
            stringstream ss(line);
            string part;
            
            while (getline(ss, part, '|')) {
                parts.push_back(part);
            }
            
            if (parts.size() >= 8) {
                SolvedPosition pos;
                pos.position_key = parts[0];
                pos.best_move = parts[1];
                pos.M_value = stoi(parts[2]);
                pos.total_plies = stoi(parts[3]);
                pos.white_moves = stoi(parts[4]);
                pos.black_moves = stoi(parts[5]);
                pos.nodes_evaluated = stoi(parts[6]);
                pos.computation_time = stod(parts[7]);
                
                // Parse BN trajectory
                if (!parts[8].empty()) {
                    stringstream bn_stream(parts[8]);
                    string bn_part;
                    while (getline(bn_stream, bn_part, ',')) {
                        pos.BN_trajectory.push_back(stoi(bn_part));
                    }
                }
                
                solved[pos.position_key] = pos;
            }
        }
        
        file.close();
        cout << "[Database] Loaded " << solved.size() << " positions from " << filename << "\n";
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
    
    bool is_attacked_by_queen(const Position& pos, const Position& qp) const {
        if (pos == qp) return false;
        if (pos.file == qp.file || pos.rank == qp.rank) return true;
        if (abs(pos.file - qp.file) == abs(pos.rank - qp.rank)) return true;
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
        if (!is_attacked_by_queen(st.bk, st.wq)) return false;
        for (auto& m : generate_all_king_moves(st.bk)) {
            if (is_attacked_by_queen(m, st.wq)) continue;
            if (m.distance_to(st.wk) <= 1) continue;
            return false;
        }
        return true;
    }
    
    bool is_stalemate(const GameState& st) const {
        if (st.to_move != 'B') return false;
        if (is_attacked_by_queen(st.bk, st.wq)) return false;
        for (auto& m : generate_all_king_moves(st.bk)) {
            if (is_attacked_by_queen(m, st.wq)) continue;
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
                GameState ns(st.wk, wq_n, st.bk, 'B');
                if (is_legal_state(ns) && !is_stalemate(ns)) cnt++;
            }
            return cnt;
        } else {
            int cnt = 0;
            for (auto& bk_n : generate_all_king_moves(st.bk)) {
                if (is_attacked_by_queen(bk_n, st.wq)) continue;
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
                if (is_attacked_by_queen(bk_n, st.wq)) continue;
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
                    GameState ns(curr.wk, wq_n, curr.bk, 'B');
                    if (is_legal_state(ns) && !is_stalemate(ns)) cands.push_back(ns);
                }
                if (cands.empty()) break;
                
                GameState best = cands[0];
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
                    if (is_attacked_by_queen(bk_n, curr.wq)) continue;
                    if (bk_n.distance_to(curr.wk) <= 1) continue;
                    GameState ns(curr.wk, curr.wq, bk_n, 'W');
                    if (is_legal_state(ns)) cands.push_back(ns);
                }
                if (cands.empty()) break;
                
                GameState best = cands[0];
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
                if (wq_n.distance_to(st.bk) > st.wq.distance_to(st.bk)) continue;
                if (wq_n.distance_to(st.bk) < 2 && wq_n.distance_to(st.wk) > 1) continue;
                GameState ns(st.wk, wq_n, st.bk, 'B');
                if (is_legal_state(ns) && !is_stalemate(ns)) cands.push_back(ns);
            }
        } else {
            for (auto& bk_n : generate_all_king_moves(st.bk)) {
                //cout << "\n  Checking BK move: " << st.bk.str() << " -> " << bk_n.str();
                if (is_attacked_by_queen(bk_n, st.wq)) continue;
                if (bk_n.distance_to(st.wk) <= 1) continue;
                GameState ns(st.wk, st.wq, bk_n, 'W');
                if (is_legal_state(ns)) cands.push_back(ns);
            }
        }
        
        // DEBUG: Show Black candidates
        if (st.to_move == 'B' && ply < 2) {
            //cout << "\n[BLACK MOVE DEBUG] Position: " << st.str() << "\n";
            //cout << "Black candidates generated:\n";
            for (size_t i = 0; i < cands.size(); i++) {
                // cout << "  " << i << ". " << cands[i].str() << "\n";
            }
        }
        
        if (cands.empty()) {
            auto res = make_pair(optional<int>(), optional<GameState>());
            memo[cache_key] = res;
            return res;
        }
        
        candidates_measured += cands.size();
        vector<tuple<GameState,int,int>> meas;
        for (auto& c : cands) {
            uint64_t m_key = make_cache_key(c, 3);
            int M;
            
            auto it = M_cache.find(m_key);
            if (it != M_cache.end()) {
                M = it->second;
            } else {
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
        
        // DEBUG: Show measured M values and best_M (NOW DECLARED)
        if (st.to_move == 'B' && ply < 2) {
            //cout << "Measured M values:\n";
            for (auto& [c,M,BN] : meas) {
                //cout << "  " << c.str() << " M=" << M << "\n";
            }
            //cout << "best_M=" << best_M << " direction=" << dir << "\n";
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
            
            // DEBUG: Show what each viable move returns
            if (st.to_move == 'B' && ply < 2) {
                //cout << "  Viable[" << i << "] " << viable[i].str() << " returned val=";
                //if (val) cout << *val; else cout << "NONE";
                //cout << "\n";
            }
            
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
                    } else if (dir == "maximize" && v > *best_val) {
                        best_val = v;
                        best_mv = viable[i];
                    } else if (dir == "maximize" && v == *best_val) {
                        // ← TIE FOR BLACK: use topological M as tiebreaker
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
            
            if (st.to_move == 'W') {
                // WHITE: Return immediately on first value found
                if (md) {
                    return make_pair(bm, md);
                }
            } else {
                // BLACK: Continue searching until max_depth regardless
                // Only update if we found a value (or keep previous best)
                if (md) {
                    best_value = md;
                    best_move = bm;
                }
            }
        }
        
        return make_pair(best_move, best_value);
    }
    
    tuple<vector<string>, int, vector<int>> play_complete_game(
        const GameState& first, int max_moves = 50, bool debug = false, int game_search_depth = 8,
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
// MAIN
// ============================================================================

int main(int argc, char* argv[]) {
    bool debug_en = false;
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        if (arg == "--debug" || arg == "-d") debug_en = true;
    }
    
    double root_t = 0;
    
    cout << "\n" << string(80, '=') << "\n";
    cout << "COMPOSITIONAL TOPOLOGICAL SEARCH\n";
    cout << "Trajectory Measurement with Black Node Count Accumulation\n";
    cout << string(80, '=') << "\n";
    
    GameState init_st(Position::from_str("d3"), Position::from_str("b2"), 
                      Position::from_str("g1"), 'W');
    cout << "\nInitial position: " << init_st.str() << "\n";
    
    cout << "\n" << string(80, '=') << "\n";
    cout << "ROOT CANDIDATE EVALUATION\n";
    cout << string(80, '=') << "\n\n";
    
    auto root_t_start = chrono::high_resolution_clock::now();
    CompositionalEngine eng;
    
    cout << "Finding optimal root evaluation depth...\n";
    cout << "(Will iterate: depth >= min_M(s) is termination condition)\n\n";
    
    int opt_d = -1;
    vector<tuple<GameState, int, int>> opt_meas;
    SolvedPositionDatabase db("kqvk_perfect_play.db");
    auto t_start = chrono::high_resolution_clock::now();
    for (int td = 1; td <= 10; td++) {
        cout << "Testing depth " << td << "..."; cout.flush();
        
        vector<GameState> cands;
        cands.reserve(64);
        for (auto& wk_n : eng.generate_all_king_moves(init_st.wk)) {
            if (wk_n.distance_to(init_st.bk) > init_st.wk.distance_to(init_st.bk)) continue;
            GameState ns(wk_n, init_st.wq, init_st.bk, 'B');
            if (eng.is_legal_state(ns) && !eng.is_stalemate(ns)) cands.push_back(ns);
        }
        for (auto& wq_n : eng.generate_all_queen_moves(init_st.wq)) {
            if (wq_n.distance_to(init_st.bk) < 2 && wq_n.distance_to(init_st.wk) > 1) continue;
            GameState ns(init_st.wk, wq_n, init_st.bk, 'B');
            if (eng.is_legal_state(ns) && !eng.is_stalemate(ns)) cands.push_back(ns);
        }
        
        vector<tuple<GameState, int, int>> meas;
        for (auto& c : cands) {
            int M = eng.compute_M_shallow(c, td);
            int BN = eng.measure_black_nodes_after_trajectory(c, td);
            meas.push_back({c, M, BN});
        }
        
        int min_M = INF_MAX;
        for (auto& [c,M,BN] : meas) min_M = min(min_M, M);
        cout << " min_M(s) = " << min_M << "\n";
        
        if (td >= min_M) {
            cout << "\n✓ SUFFICIENT DEPTH: " << td << " >= " << min_M << "\n";
            cout << "  Optimal root evaluation depth: " << td << "\n";
            opt_d = td;
            opt_meas = meas;
            break;
        }
    }
    
    if (opt_d == -1) {

        opt_d = 7;
        vector<GameState> cands;
        cands.reserve(64);
        // White's turn - White King shouldn't move away from Black King
        for (auto& wk_n : eng.generate_all_king_moves(init_st.wk)) {
            if (wk_n.distance_to(init_st.bk) > init_st.wk.distance_to(init_st.bk)) continue;  // Prune if moving away
            GameState ns(wk_n, init_st.wq, init_st.bk, 'B');
            if (eng.is_legal_state(ns) && !eng.is_stalemate(ns)) cands.push_back(ns);
        }
        for (auto& wq_n : eng.generate_all_queen_moves(init_st.wq)) {
            if (wq_n.distance_to(init_st.bk) < 2 && wq_n.distance_to(init_st.wk) > 1) continue;
            GameState ns(init_st.wk, wq_n, init_st.bk, 'B');
            if (eng.is_legal_state(ns) && !eng.is_stalemate(ns)) cands.push_back(ns);
        }
        for (auto& c : cands) {
            int M = eng.compute_M_shallow(c, opt_d);
            int BN = eng.measure_black_nodes_after_trajectory(c, opt_d);
            opt_meas.push_back({c, M, BN});
        }
    }
    
    sort(opt_meas.begin(), opt_meas.end(), [](auto& a, auto& b) {
        if (get<1>(a) != get<1>(b)) return get<1>(a) < get<1>(b);
        return get<2>(a) < get<2>(b);
    });
    
    auto root_t_end = chrono::high_resolution_clock::now();
    root_t = chrono::duration<double>(root_t_end - root_t_start).count();
    
    cout << "Total root candidates: " << opt_meas.size() << "\n\n";
    cout << "Top 10 candidates by M(s):\n";
    for (size_t i = 0; i < min(size_t(10), opt_meas.size()); i++) {
        cout << "  " << (i+1) << ". " << get<0>(opt_meas[i]).str()
             << " M=" << get<1>(opt_meas[i])
             << ", Black nodes=" << get<2>(opt_meas[i]) << "\n";
    }
    
    int min_M_v = INF_MAX;
    for (auto& [c,M,BN] : opt_meas) min_M_v = min(min_M_v, M);
    
    vector<tuple<GameState, int, int>> opt_cands;
    for (auto& [c,M,BN] : opt_meas) {
        if (M == min_M_v) opt_cands.push_back({c, M, BN});
    }
    
    cout << "\nCandidates with M(s) = " << min_M_v << ": " << opt_cands.size() << "\n";
    for (size_t i = 0; i < opt_cands.size(); i++) {
        cout << "  " << (i+1) << ". " << get<0>(opt_cands[i]).str()
             << " M=" << get<1>(opt_cands[i])
             << ", Black nodes=" << get<2>(opt_cands[i]) << "\n";
    }
    
    cout << "\nFull endgame evaluation for all " << opt_cands.size() << " candidates\n";
    cout << "Root evaluation time: " << fixed << setprecision(2) << root_t << "s\n";
    
    cout << "\n" << string(80, '=') << "\n";
    cout << "COMPLETE GAME EVALUATION\n";
    cout << string(80, '=') << "\n";
    
    vector<double> c_times;
    
    for (size_t idx = 0; idx < opt_cands.size(); idx++) {
        auto [cand, M_r, BN_r] = opt_cands[idx];
        
        cout << "\n" << string(80, '=') << "\n";
        cout << "CANDIDATE " << (idx+1) << ": " << cand.str() << "\n";
        cout << "M(s) at root: " << M_r << "\n";
        cout << "Black nodes at root: " << BN_r << "\n";
        cout << string(80, '=') << "\n\n";
        
        eng.nodes_evaluated = 0;
        eng.candidates_measured = 0;
        
        auto c_start = chrono::high_resolution_clock::now();
        auto [mvs, tp, bnc] = eng.play_complete_game(cand, 50, debug_en, M_r);
        auto c_end = chrono::high_resolution_clock::now();
        double c_el = chrono::duration<double>(c_end - c_start).count();
        c_times.push_back(c_el);
        
                // Record this solution in the database
        SolvedPosition solution;
        solution.position_key = cand.str();
        solution.best_move = mvs.empty() ? "none" : mvs[0];
        solution.M_value = M_r;
        solution.BN_trajectory = bnc;
        solution.total_plies = tp;
        solution.white_moves = (tp + 1) / 2;
        solution.black_moves = tp / 2;
        solution.nodes_evaluated = eng.nodes_evaluated;
        solution.computation_time = c_el;
        
        db.add_position(solution);

        cout << "Complete game: ";
        for (size_t i = 0; i < mvs.size(); i++) {
            if (i > 0) cout << " ";
            cout << mvs[i];
        }
        cout << "\n";
        cout << "Total plies: " << tp << "\n";
        cout << "White moves: " << ((tp+1)/2) << "\n";
        cout << "Black moves: " << (tp/2) << "\n";
        cout << "Nodes evaluated: " << eng.nodes_evaluated << "\n";
        cout << "Candidates measured: " << eng.candidates_measured << "\n";
        cout << "Candidate " << (idx+1) << " computation time: " << fixed << setprecision(2) << c_el << "s\n";
    }
    
    auto t_end = chrono::high_resolution_clock::now();
    double tot_t = chrono::duration<double>(t_end - t_start).count();
    
    cout << "\n" << string(80, '=') << "\n";
    cout << "TIMING ANALYSIS\n";
    cout << string(80, '=') << "\n\n";
    
    cout << "Root candidate evaluation:    " << fixed << setprecision(2) << setw(10) << root_t << "s\n";
    for (size_t i = 0; i < c_times.size(); i++) {
        cout << "Candidate " << (i+1) << " computation:      " << setw(10) << c_times[i] << "s\n";
    }
    cout << string(40, '-') << "\n";
    cout << "Total computation time:       " << setw(10) << tot_t << "s\n";
    cout << "\n";
    db.export_to_file();
    db.print_summary();
    return 0;
}
