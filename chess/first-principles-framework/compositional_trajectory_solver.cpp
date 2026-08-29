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

using namespace std;

const int INF_MAX = 2147483647;
const int INF_MIN = -2147483648;

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
            if (ply == 0) cout << "[DEBUG] *** MEMO HIT at depth=" << depth << " ***\n";
            return memo[cache_key];
        }
        
        if (ply == 0) {
            cout << "\n[DEBUG compositional_search] START depth=" << depth << ", position=" << st.str()
                 << ", to_move=" << st.to_move << "\n";
        }
        
        if (is_checkmate(st)) {
            if (ply == 0) cout << "[DEBUG] CHECKMATE at start\n";
            auto res = make_pair(optional<int>(0), optional<GameState>());
            memo[cache_key] = res;
            return res;
        }
        
        if (depth == 0) {
            if (ply == 0) cout << "[DEBUG] Depth limit reached\n";
            auto res = make_pair(optional<int>(), optional<GameState>());
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
                GameState ns(st.wk, wq_n, st.bk, 'B');
                if (is_legal_state(ns) && !is_stalemate(ns)) cands.push_back(ns);
            }
        } else {
            for (auto& bk_n : generate_all_king_moves(st.bk)) {
                if (is_attacked_by_queen(bk_n, st.wq)) continue;
                if (bk_n.distance_to(st.wk) <= 1) continue;
                GameState ns(st.wk, st.wq, bk_n, 'W');
                if (is_legal_state(ns)) cands.push_back(ns);
            }
        }
        
        if (ply == 0) cout << "[DEBUG] Generated " << cands.size() << " candidates\n";
        
        if (cands.empty()) {
            if (ply == 0) {
                cout << "[DEBUG] NO CANDIDATES FOUND!\n";
                cout << "[DEBUG] Position: " << st.str() << "\n";
                cout << "[DEBUG] to_move: " << st.to_move << "\n";
            }
            auto res = make_pair(optional<int>(), optional<GameState>());
            memo[cache_key] = res;
            return res;
        }
        
        candidates_measured += cands.size();
        vector<tuple<GameState,int,int>> meas;
        for (auto& c : cands) {
            int M = compute_M_shallow(c, 3);
            int BN = count_legal_moves(c);
            meas.push_back({c, M, BN});
        }
        
        if (ply == 0) cout << "[DEBUG] Measured " << meas.size() << " candidates\n";
        
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
        
        if (ply == 0) cout << "[DEBUG] best_M=" << best_M << ", direction=" << dir << "\n";
        
        int thresh = max(2, best_M / 3);
        vector<GameState> viable;
        for (auto& [c,M,BN] : meas) {
            if (dir == "minimize" && M <= best_M + thresh) viable.push_back(c);
            else if (dir == "maximize" && M >= best_M - thresh) viable.push_back(c);
        }
        
        if (ply == 0) cout << "[DEBUG] threshold=" << thresh << ", viable=" << viable.size() << "/" << cands.size() << "\n";
        
        if (viable.size() < 3) {
            if (ply == 0) cout << "[DEBUG] Increasing viable from " << viable.size() << " to 3\n";
            sort(meas.begin(), meas.end(), [](auto& a, auto& b) { return get<1>(a) < get<1>(b); });
            viable.clear();
            for (size_t i = 0; i < min(size_t(3), meas.size()); i++) {
                viable.push_back(get<0>(meas[i]));
            }
        }
        
        if (debug && ply < 4) {
            string ind(ply*2, ' ');
            double ratio = cands.empty() ? 0 : (double)viable.size()/cands.size();
            cout << ind << "[Ply " << ply << "] " << st.to_move << ": " << cands.size()
                 << " -> " << viable.size() << " (" << fixed << setprecision(1) << (ratio*100)
                 << "%), best_M=" << best_M << ", thresh=" << thresh << "\n";
        }
        
        optional<int> best_val;
        optional<GameState> best_mv;
        
        if (ply == 0) cout << "[DEBUG] Recursing on " << viable.size() << " candidates\n";
        
        for (size_t i = 0; i < viable.size(); i++) {
            auto [val, mv] = compositional_search_impl(viable[i], depth-1, ply+1, debug, memo);
            nodes_evaluated++;
            
            if (ply == 0) {
                cout << "[DEBUG] Candidate " << (i+1) << "/" << viable.size() << ": value=";
                if (val) cout << *val; else cout << "None";
                cout << "\n";
            }
            
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
                }
            }
        }
        
        if (ply == 0) {
            cout << "[DEBUG] RESULT depth=" << depth << ": best_value=";
            if (best_val) cout << *best_val; else cout << "None";
            cout << ", best_move=" << (best_mv ? "True" : "False") << "\n";
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
        
        cout << "\n[DEBUG find_best_move] Starting search from: " << st.str() << "\n";
        cout << "[DEBUG find_best_move] to_move: " << st.to_move << "\n";
        cout << "[DEBUG find_best_move] max_depth: " << max_depth << "\n";
        
        unordered_map<uint64_t, pair<optional<int>, optional<GameState>>> memo;
        
        for (int depth = 2; depth <= max_depth + 1; depth += 2) {
            cout << "\n[DEBUG find_best_move] Testing depth " << depth << "...\n";
            auto [md, bm] = compositional_search_impl(st, depth, 0, debug, memo);
            
            cout << "[DEBUG find_best_move] depth=" << depth << ": mate_dist=";
            if (md) cout << *md; else cout << "None";
            cout << ", best_move=" << (bm ? "True" : "False") << "\n";
            
            if (md) {
                cout << "[DEBUG find_best_move] FOUND at depth " << depth << "\n";
                return make_pair(bm, md);
            } else {
                cout << "[DEBUG find_best_move] Continuing...\n";
            }
        }
        
        cout << "\n[DEBUG find_best_move] EXHAUSTED ALL DEPTHS\n";
        return make_pair(optional<GameState>(), optional<int>());
    }
    
    tuple<vector<string>, int, vector<int>> play_complete_game(
        const GameState& first, int max_moves = 50, bool debug = false, int game_search_depth = 8
    ) {
        vector<string> mvs;
        GameState curr = first;
        vector<int> bnc;
        
        cout << "\n[DEBUG play_complete_game] Starting game evaluation\n";
        cout << "[DEBUG] Game search depth: " << game_search_depth << "\n";
        cout << "[DEBUG] Starting position: " << curr.str() << "\n";
        cout << "[DEBUG] to_move: " << curr.to_move << "\n";
        
        for (int move_num = 0; move_num < max_moves; move_num++) {
            cout << "\n[DEBUG] ========== MOVE " << (move_num+1) << " ==========\n";
            cout << "[DEBUG] Current position: " << curr.str() << "\n";
            cout << "[DEBUG] to_move: " << curr.to_move << "\n";
            
            if (is_checkmate(curr)) {
                cout << "[DEBUG] CHECKMATE DETECTED - Game over\n";
                return make_tuple(mvs, (int)mvs.size(), bnc);
            }
            
            cout << "[DEBUG] Calling find_best_move with max_depth=" << game_search_depth << "\n";
            auto [ns, md] = find_best_move(curr, 2*game_search_depth, debug);
            
            if (!ns) {
                cout << "[DEBUG] NO NEXT STATE FOUND\n";
                return make_tuple(mvs, (int)mvs.size(), bnc);
            }
            
            string mv_str = get_move_notation(curr, *ns);
            cout << "[DEBUG] Move notation: " << mv_str << "\n";
            mvs.push_back(mv_str);
            
            int bn = count_legal_moves(*ns);
            bnc.push_back(bn);
            cout << "[DEBUG] Black nodes at position: " << bn << "\n";
            
            curr = *ns;
        }
        
        cout << "[DEBUG] Reached max_moves limit\n";
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
    
    auto t_start = chrono::high_resolution_clock::now();
    double root_t = 0;
    
    cout << "\n" << string(80, '=') << "\n";
    cout << "COMPOSITIONAL TOPOLOGICAL SEARCH\n";
    cout << "Trajectory Measurement with Black Node Count Accumulation\n";
    if (debug_en) cout << "[DEBUG MODE ENABLED]\n";
    cout << string(80, '=') << "\n";
    
    GameState init_st(Position::from_str("g6"), Position::from_str("f6"), 
                    Position::from_str("d2"), 'W');
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
    
    for (int td = 1; td <= 7; td++) {
        cout << "Testing depth " << td << "..."; cout.flush();
        
        vector<GameState> cands;
        cands.reserve(64);
        for (auto& wk_n : eng.generate_all_king_moves(init_st.wk)) {
            GameState ns(wk_n, init_st.wq, init_st.bk, 'B');
            if (eng.is_legal_state(ns) && !eng.is_stalemate(ns)) cands.push_back(ns);
        }
        for (auto& wq_n : eng.generate_all_queen_moves(init_st.wq)) {
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
        for (auto& wk_n : eng.generate_all_king_moves(init_st.wk)) {
            GameState ns(wk_n, init_st.wq, init_st.bk, 'B');
            if (eng.is_legal_state(ns) && !eng.is_stalemate(ns)) cands.push_back(ns);
        }
        for (auto& wq_n : eng.generate_all_queen_moves(init_st.wq)) {
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
    
    cout << "\nWith C++ optimization (50x speedup):\n";
    cout << "  Root eval:     " << fixed << setprecision(3) << (root_t/50) << "s\n";
    for (size_t i = 0; i < c_times.size(); i++) {
        cout << "  Candidate " << (i+1) << ":   " << (c_times[i]/50) << "s\n";
    }
    cout << "  Total:         " << (tot_t/50) << "s\n";
    
    return 0;
}
