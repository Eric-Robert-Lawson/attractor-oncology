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

    // For pawn endgames only: which piece currently occupies the "WQ:" square --
    // 'P' (still a pawn), 'Q' (promoted to queen), or 'N' (promoted to knight).
    // Necessary because the position STRING alone is ambiguous: "WK:e2 WQ:e4
    // BK:e8" could be either a pawn on e4 or an already-promoted queen that
    // later moved there -- these are genuinely different game states requiring
    // different move rules, and the string can't distinguish them on its own.
    // Defaults to 'Q' so every EXISTING KQvK/KRvK row (which never had this
    // column) loads with a value that's simply never consulted by non-pawn
    // engines -- see load_from_file's backward-compatible parsing below.
    char attacker_kind = 'Q';

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
        ss << "|" << attacker_kind;
        ss << "\n";
        return ss.str();
    }
    
    static string csv_header() {
        return "Position|Turn|BestMove|M|Plies|WhiteMoves|BlackMoves|NodesEval|Time|CumulativeBN|"
               "BN_Trajectory|TiedMoves|AttackerKind\n";
    }
};

// ============================================================================
// SolvedPositionDatabase - Manages the tablebase
// ============================================================================

class SolvedPositionDatabase {
    private:
        map<pair<string, char>, SolvedPosition> solved;
        string filename;
        // Keys added or mutated since the last append_new_to_file() call. Tracking
        // ONLY the delta is what makes frequent, cheap on-disk checkpointing possible --
        // see append_new_to_file() below for why this replaces repeated full rewrites.
        vector<pair<string, char>> pending_export;
        
    public:
        SolvedPositionDatabase(const string& db_file = "kqvk_perfect_play.db")
            : filename(db_file) {
            load_from_file();
        }
        
        // Add a solved position to the database.
        //
        // full_symmetry controls how many of the 8 D4 board symmetries get used
        // to generate derived positions. Default true preserves EXACT existing
        // behavior (all 8) for every piece used so far (King, Queen, Rook --
        // Bishop/Knight below too): their movement rules depend only on
        // relative file/rank displacement, which every one of the 8 D4
        // transforms preserves (rotations and reflections permute "same file",
        // "same rank", and "same diagonal" among themselves, but never turn
        // one into something that isn't a line-of-movement relation at all).
        //
        // Pass full_symmetry=false whenever a PAWN is on the board. A pawn's
        // legality depends on an ABSOLUTE direction (forward = toward higher
        // rank for White, lower rank for Black) and ABSOLUTE ranks (promotion
        // on rank 8, double-step only from rank 2) -- properties tied to a
        // fixed axis, not to relative displacement. Every D4 transform except
        // identity and flip_horizontal either swaps the file/rank axes
        // (turning "forward" into "sideways": both rotations of 90/270, and
        // both diagonal flips) or reverses the rank axis outright (turning a
        // forward pawn push into a backward one: rotate_180 and flip_vertical).
        // flip_horizontal is the sole nontrivial survivor because it only ever
        // touches file, leaving rank -- and therefore "forward" and every
        // absolute-rank rule -- completely untouched. Verified directly: a
        // White d2-d4 push transforms to a legal e2-e4 push under
        // flip_horizontal, but to illegal or backward moves under every other
        // one of the remaining 6 transforms.
        void add_position(const SolvedPosition& pos, bool full_symmetry = true) {
            if (is_solved(pos.position_key, pos.turn)) {
                return;
            }
            solved[{pos.position_key, pos.turn}] = pos;
            pending_export.push_back({pos.position_key, pos.turn});
            
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
                    pending_export.push_back({new_pos, pos.turn});
                }
            };
            
            // Generate 4 rotations -- ONLY when full_symmetry, since a rotation
            // swaps the file/rank axes and is never valid with a pawn present
            // (see the comment on add_position's signature above).
            string current_pos = pos.position_key;
            string current_move = pos.best_move;
            vector<pair<string,int>> current_tied = pos.tied_moves;

            if (full_symmetry) {
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
            }

            // flip_horizontal (mirror left-right, rank UNCHANGED) is valid
            // regardless of pawns -- always applied.
            add_transformed(flip_horizontal, pos.position_key, pos.best_move, pos.tied_moves);

            // flip_vertical (reverses rank -- turns forward into backward) and
            // both diagonal flips (swap file/rank axes, same problem as
            // rotation) are only valid without a pawn on the board.
            if (full_symmetry) {
                add_transformed(flip_vertical, pos.position_key, pos.best_move, pos.tied_moves);
                add_transformed(flip_diagonal_a1h8, pos.position_key, pos.best_move, pos.tied_moves);
                add_transformed(flip_diagonal_a8h1, pos.position_key, pos.best_move, pos.tied_moves);
            }
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
            pending_export.push_back({position_key, turn});
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
        pending_export.clear();  // this rewrite already covers everything pending
    }

    // Appends ONLY the rows added or mutated since the last flush (export_to_file
    // or append_new_to_file) -- O(size of the delta), not O(total database size).
    // This is what makes it cheap to checkpoint constantly instead of every few
    // thousand positions: export_to_file() rewrites the entire file every single
    // call (that's what the original single-path driver's "every 5 positions"
    // checkpoint was doing, and it compounds badly as the file grows), whereas
    // this only ever writes what's new.
    //
    // One consequence worth knowing: if a key was already written to disk in an
    // earlier flush and is later MUTATED (this happens for set_tied_moves backfills
    // on legacy rows -- see that method), this appends a SECOND, updated line for
    // the same key rather than editing the earlier line in place. That's safe, not
    // silent corruption: load_from_file() reads the file top-to-bottom and does an
    // unconditional map assignment per row, so whichever occurrence of a key comes
    // LAST in the file wins on reload -- appends always land after whatever was
    // already there, so the final in-memory state after reloading is always
    // correct. The only cost is a few harmless stale duplicate lines accumulating
    // for positions that got backfilled mid-run, which export_to_file()'s full
    // rewrite (deduplicated by construction, since it iterates the in-memory map)
    // cleans up the next time it's called -- e.g. once at the very end of a run.
    void append_new_to_file() {
        if (pending_export.empty()) return;

        ifstream check(filename);
        bool need_header = !check.good() || check.peek() == std::ifstream::traits_type::eof();
        check.close();

        ofstream file(filename, ios::app);
        if (!file.is_open()) {
            cerr << "ERROR: Cannot open database file for append: " << filename << "\n";
            return;
        }
        if (need_header) {
            file << SolvedPosition::csv_header();
        }
        for (auto& key : pending_export) {
            auto it = solved.find(key);
            if (it != solved.end()) {
                file << it->second.to_csv();
            }
        }
        file.close();
        pending_export.clear();
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

                    // AttackerKind: 13th column, added for pawn-endgame support. Older
                    // database files (written before this existed) simply won't have it --
                    // parts.size() <= 12 leaves the struct's default ('Q') in place, which
                    // is exactly correct for every pre-existing KQvK/KRvK row: it's never
                    // consulted by anything except the pawn engine.
                    if (parts.size() > 12 && !parts[12].empty()) {
                        pos.attacker_kind = parts[12][0];
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
// AttackerRules -- the piece-agnostic interface
// ============================================================================
//
// A compile-time (template parameter) interface, not a virtual base class --
// this is a search-heavy engine and template instantiation lets the compiler
// fully inline move generation and attack checks, with zero runtime dispatch
// overhead versus the original hardcoded-to-queen code. Each endgame
// (KQvK, KRvK, ...) compiles to its own binary from the same source, chosen
// by which Rules struct BaseEngine/CompositionalEngine are instantiated with.
//
// Required members of any Rules struct:
//   static constexpr char letter        -- algebraic notation letter (e.g. 'Q')
//   static vector<Position> generate_moves(const Position& pos)
//       Pseudo-legal move generation: every square this piece could reach
//       from `pos` on an otherwise-empty board. Occupancy/blocking is
//       checked separately by the caller (matches how generate_all_queen_
//       moves already worked -- it never looked at other pieces).
//   static bool attacks(const Position& target, const Position& piece_pos,
//                        const Position& blocker_a, const Position& blocker_b)
//       Does the piece at `piece_pos` attack `target`, given up to two other
//       occupied squares (blocker_a, blocker_b) that can block the line?
//       Mirrors is_attacked_by_queen's exact signature and semantics,
//       including piece_pos==target returning false (used contextually to
//       mean "you can't be attacked by the square you're moving onto", which
//       lets a king safely walk onto an undefended attacker's square).
//
// QueenRules below is a direct, behavior-preserving extraction of the
// pre-refactor generate_all_queen_moves/is_attacked_by_queen -- verified via
// a full regression harness comparing this binary's output against the
// pre-refactor original, not just by inspection.
// ============================================================================

struct QueenRules {
    static constexpr char letter = 'Q';
    // All 8 D4 board symmetries are valid for a queen -- its legality depends
    // only on relative file/rank displacement (same file, same rank, same
    // diagonal), and every D4 transform preserves that set of relations.
    static constexpr bool full_board_symmetry = true;

    // Exposes the raw direction set for reuse by the general multi-piece
    // engine's blocking-aware sliding logic -- purely additive, doesn't change
    // anything about how generate_moves/attacks already work.
    static vector<pair<int,int>> directions() {
        return {{1,0},{-1,0},{0,1},{0,-1},{1,1},{1,-1},{-1,1},{-1,-1}};
    }

    static vector<Position> generate_moves(const Position& pos) {
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

    static bool attacks(const Position& pos, const Position& qp,
                         const Position& wk, const Position& wq) {
        if (pos == qp) return false;

        // Check file
        if (pos.file == qp.file) {
            int start = min(pos.rank, qp.rank) + 1;
            int end = max(pos.rank, qp.rank);
            for (int r = start; r < end; r++) {
                if (Position(pos.file, r) == wk || Position(pos.file, r) == wq) return false;
            }
            return true;
        }

        // Check rank
        if (pos.rank == qp.rank) {
            int start = min(pos.file, qp.file) + 1;
            int end = max(pos.file, qp.file);
            for (int f = start; f < end; f++) {
                if (Position(f, pos.rank) == wk || Position(f, pos.rank) == wq) return false;
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
                if (Position(f, r) == wk || Position(f, r) == wq) return false;
                f += df;
                r += dr;
            }
            return true;
        }

        return false;
    }
};

// RookRules -- a direct sibling of QueenRules, restricted to the 4
// orthogonal directions (no diagonals). Everything else about how it plugs
// into BaseEngine/CompositionalEngine is identical to the queen case.
struct RookRules {
    static constexpr char letter = 'R';
    static constexpr bool full_board_symmetry = true;  // same reasoning as QueenRules

    // Exposes the raw direction set for reuse by the general multi-piece
    // engine's blocking-aware sliding logic -- purely additive, doesn't change
    // anything about how generate_moves/attacks already work.
    static vector<pair<int,int>> directions() {
        return {{1,0},{-1,0},{0,1},{0,-1}};
    }

    static vector<Position> generate_moves(const Position& pos) {
        vector<Position> moves;
        int dirs[][2] = {{1,0},{-1,0},{0,1},{0,-1}};  // orthogonal only

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

    static bool attacks(const Position& pos, const Position& rp,
                         const Position& wk, const Position& wr) {
        if (pos == rp) return false;

        // Check file
        if (pos.file == rp.file) {
            int start = min(pos.rank, rp.rank) + 1;
            int end = max(pos.rank, rp.rank);
            for (int r = start; r < end; r++) {
                if (Position(pos.file, r) == wk || Position(pos.file, r) == wr) return false;
            }
            return true;
        }

        // Check rank
        if (pos.rank == rp.rank) {
            int start = min(pos.file, rp.file) + 1;
            int end = max(pos.file, rp.file);
            for (int f = start; f < end; f++) {
                if (Position(f, pos.rank) == wk || Position(f, pos.rank) == wr) return false;
            }
            return true;
        }

        // No diagonal attack -- a rook simply doesn't threaten off-file/off-rank
        // squares, unlike the queen's third branch above.
        return false;
    }
};

// BishopRules and KnightRules -- genuinely modular in the same sense Rook was:
// both pieces' legality depends only on relative file/rank displacement, so
// all 8 D4 symmetries apply cleanly, and neither needs anything beyond the
// existing AttackerRules interface (pseudo-legal move generation + an attacks
// check). Cross-validated against hand-computed geometry before delivery.
//
// IMPORTANT, not a caveat to skip: a LONE bishop or LONE knight can NEVER
// force checkmate against a lone king, under any circumstance, from any
// starting position -- this is elementary, universal chess theory (the
// "insufficient material" rule), not a limitation of this engine. Running
// --full-dag or --batch with either of these instantiated as the sole
// attacker will find precisely nothing: every position is an unconditional
// draw, so there is no winning technique to compare and no tied-move
// structure to analyze at all. They're included here because the geometry
// itself is correct and reusable -- e.g. as one of several pieces in a
// future multi-attacker endgame (KBNvK, which unlike either piece alone CAN
// force mate, being the classic hard case; or as a supporting piece in
// KQBvK, KRNvK, etc.) -- not because KBvK or KNvK are meaningful targets to
// run standalone.

struct BishopRules {
    static constexpr char letter = 'B';
    static constexpr bool full_board_symmetry = true;

    // Exposes the raw direction set for reuse by the general multi-piece
    // engine's blocking-aware sliding logic -- purely additive, doesn't change
    // anything about how generate_moves/attacks already work.
    static vector<pair<int,int>> directions() {
        return {{1,1},{1,-1},{-1,1},{-1,-1}};
    }

    static vector<Position> generate_moves(const Position& pos) {
        vector<Position> moves;
        int dirs[][2] = {{1,1},{1,-1},{-1,1},{-1,-1}};  // diagonals only

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

    static bool attacks(const Position& pos, const Position& bp,
                         const Position& wk, const Position& wb) {
        if (pos == bp) return false;

        // Diagonal only -- no file/rank branches at all, unlike Queen/Rook.
        if (abs(pos.file - bp.file) == abs(pos.rank - bp.rank)) {
            int df = (pos.file > bp.file) ? 1 : -1;
            int dr = (pos.rank > bp.rank) ? 1 : -1;
            int f = bp.file + df;
            int r = bp.rank + dr;
            while (f != pos.file) {
                if (Position(f, r) == wk || Position(f, r) == wb) return false;
                f += df;
                r += dr;
            }
            return true;
        }
        return false;
    }
};

struct KnightRules {
    static constexpr char letter = 'N';
    static constexpr bool full_board_symmetry = true;

    static vector<Position> generate_moves(const Position& pos) {
        vector<Position> moves;
        int offsets[][2] = {{1,2},{2,1},{2,-1},{1,-2},{-1,-2},{-2,-1},{-2,1},{-1,2}};
        for (auto& o : offsets) {
            int nf = pos.file + o[0];
            int nr = pos.rank + o[1];
            if (nf >= 0 && nf <= 7 && nr >= 0 && nr <= 7) {
                moves.push_back(Position(nf, nr));
            }
        }
        return moves;
    }

    static bool attacks(const Position& pos, const Position& np,
                         const Position& /*wk*/, const Position& /*wn*/) {
        if (pos == np) return false;
        // A knight jumps over pieces -- no blocking check at all, unlike
        // every sliding piece above. The two unused blocker parameters are
        // still required to satisfy the shared AttackerRules signature.
        int df = abs(pos.file - np.file);
        int dr = abs(pos.rank - np.rank);
        return (df == 1 && dr == 2) || (df == 2 && dr == 1);
    }
};

// ============================================================================
// PawnState / PawnEngine -- a genuinely separate, parallel implementation
// for pawn endgames (KPvK). Deliberately NOT built on GameState/AttackerRules.
// ============================================================================
//
// Every piece above (King, Queen, Rook, Bishop, Knight) shares two properties
// GameState and CompositionalEngine<AttackerRules> were designed around:
// (1) a piece's move squares and attack squares are the SAME set, and (2) a
// piece's type never changes during a game. A pawn breaks both -- it moves
// straight but attacks diagonally, and it can promote into a different piece
// mid-game. Rather than retrofit those two assumptions into the
// already-verified KQvK/KRvK machinery (real risk of disturbing something
// carefully proven correct), this is a self-contained parallel engine. It
// reuses everything that genuinely IS piece-agnostic already: Position,
// SolvedPosition (including the new AttackerKind column), SolvedPositionDatabase,
// and the symmetry-restriction logic from add_position's full_symmetry flag.
// QueenRules::generate_moves/attacks and KnightRules::generate_moves/attacks
// are reused directly for the post-promotion case -- no need to reimplement
// them a second time.

// KING/ROOK/BISHOP added here (purely additive -- PawnState/PawnEngine never
// construct these values, so nothing about their existing behavior changes)
// specifically so the general multi-piece engine below can reuse this same
// enum and these same helper functions instead of duplicating them.
enum class PieceKind : uint8_t { PAWN = 0, QUEEN = 1, KNIGHT = 2, KING = 3, ROOK = 4, BISHOP = 5 };

inline char kind_letter(PieceKind k) {
    switch (k) {
        case PieceKind::PAWN:   return 'P';
        case PieceKind::QUEEN:  return 'Q';
        case PieceKind::KNIGHT: return 'N';
        case PieceKind::KING:   return 'K';
        case PieceKind::ROOK:   return 'R';
        case PieceKind::BISHOP: return 'B';
    }
    return '?';
}

inline PieceKind kind_from_letter(char c) {
    if (c == 'Q') return PieceKind::QUEEN;
    if (c == 'N') return PieceKind::KNIGHT;
    if (c == 'K') return PieceKind::KING;
    if (c == 'R') return PieceKind::ROOK;
    if (c == 'B') return PieceKind::BISHOP;
    return PieceKind::PAWN;
}

class PawnState {
public:
    Position wk, wp, bk;
    PieceKind wp_kind;
    char to_move;

    PawnState() : wp_kind(PieceKind::PAWN), to_move('W') {}
    PawnState(Position wk_, Position wp_, PieceKind kind_, Position bk_, char to_move_)
        : wk(wk_), wp(wp_), bk(bk_), wp_kind(kind_), to_move(to_move_) {}

    // Label stays "WQ:" regardless of current piece kind, matching every
    // other endgame's convention -- SolvedPosition's separate AttackerKind
    // column (not this string) is the source of truth for what's actually
    // there, since the string alone can't distinguish a pawn on e4 from an
    // already-promoted queen that later moved to e4.
    string str() const {
        return "WK:" + wk.str() + " WQ:" + wp.str() + " BK:" + bk.str();
    }

    bool operator==(const PawnState& other) const {
        return wk == other.wk && wp == other.wp && wp_kind == other.wp_kind
            && bk == other.bk && to_move == other.to_move;
    }
};

struct PawnMove {
    Position dest;
    bool is_promotion;
};

// White pawn only -- the only case KPvK needs. Attacks the two
// diagonally-forward squares regardless of occupancy, exactly like every
// other piece's attacks(): occupancy determines whether a MOVE there is
// legal, not whether the square is under attack.
inline vector<Position> pawn_attack_squares(const Position& p) {
    vector<Position> squares;
    if (p.rank + 1 > 7) return squares;  // already on the back rank -- shouldn't occur mid-game
    for (int df : {-1, 1}) {
        int nf = p.file + df;
        if (nf >= 0 && nf <= 7) squares.push_back(Position(nf, p.rank + 1));
    }
    return squares;
}

// Quiet (non-capturing) forward moves only. Diagonal CAPTURE moves are
// deliberately not generated: in KPvK, Black's only piece is its king, and a
// king can never legally be captured (checkmate ends the game before that
// could happen), so a pawn's diagonal capture has no legal target in this
// material and would be dead, untested code if added now. pawn_attack_squares
// above still matters independently -- it determines which squares near
// Black's king are unsafe to move into, which has nothing to do with whether
// a capture move exists to make.
inline vector<PawnMove> generate_pawn_quiet_moves(const Position& p, const Position& wk, const Position& bk) {
    vector<PawnMove> moves;
    if (p.rank + 1 > 7) return moves;
    Position one_step(p.file, p.rank + 1);
    if (one_step == wk || one_step == bk) return moves;  // blocked -- no double-step possible either

    bool promotes = (one_step.rank == 7);
    moves.push_back({one_step, promotes});

    if (!promotes && p.rank == 1) {  // starting rank (index 1 = rank 2) -- double-step allowed
        Position two_step(p.file, p.rank + 2);
        if (!(two_step == wk) && !(two_step == bk)) {
            moves.push_back({two_step, false});
        }
    }
    return moves;
}

// Dispatches to the correct attacks() check based on whichever piece
// currently occupies the wp slot. The one place genuine runtime branching
// happens in this engine -- a plain 3-way switch, not a general N-piece scan,
// so the cost is negligible even though every AttackerRules piece above
// resolves this at compile time instead. piece_square is passed as both the
// "attacker position" and the "second blocker" argument to Queen/KnightRules,
// matching the exact calling convention used throughout the rest of this file
// (e.g. is_attacked_by_queen(m, st.wq, st.wk, st.wq)).
inline bool wp_attacks(PieceKind kind, const Position& target, const Position& piece_square, const Position& wk) {
    switch (kind) {
        case PieceKind::PAWN: {
            for (auto& s : pawn_attack_squares(piece_square)) if (s == target) return true;
            return false;
        }
        case PieceKind::QUEEN:
            return QueenRules::attacks(target, piece_square, wk, piece_square);
        case PieceKind::KNIGHT:
            return KnightRules::attacks(target, piece_square, wk, piece_square);
        default:
            // KING/ROOK/BISHOP are unreachable here -- PawnEngine (this
            // function's only caller) never assigns them to wp_kind. They
            // exist on PieceKind only for the general multi-piece engine
            // further below to reuse.
            return false;
    }
}

struct PawnSearchResult {
    optional<int> val;
    optional<int> bn_cum;
    optional<PawnState> mv;
    vector<pair<PawnState, int>> tied;
};

class PawnEngine {
public:
    int nodes_evaluated = 0;

    vector<Position> generate_all_king_moves(const Position& pos) const {
        vector<Position> moves;
        for (int df = -1; df <= 1; df++) {
            for (int dr = -1; dr <= 1; dr++) {
                if (df == 0 && dr == 0) continue;
                int nf = pos.file + df, nr = pos.rank + dr;
                if (nf >= 0 && nf <= 7 && nr >= 0 && nr <= 7) moves.push_back(Position(nf, nr));
            }
        }
        return moves;
    }

    bool is_legal_state(const PawnState& st) const {
        set<pair<int,int>> pos;
        pos.insert({st.wk.file, st.wk.rank});
        pos.insert({st.wp.file, st.wp.rank});
        pos.insert({st.bk.file, st.bk.rank});
        if (pos.size() != 3) return false;
        if (st.wk.distance_to(st.bk) < 2) return false;
        return true;
    }

    bool is_checkmate(const PawnState& st) const {
        if (st.to_move != 'B') return false;
        if (!wp_attacks(st.wp_kind, st.bk, st.wp, st.wk)) return false;
        for (auto& m : generate_all_king_moves(st.bk)) {
            if (wp_attacks(st.wp_kind, m, st.wp, st.wk)) continue;
            if (m.distance_to(st.wk) <= 1) continue;
            return false;
        }
        return true;
    }

    bool is_stalemate(const PawnState& st) const {
        if (st.to_move != 'B') return false;
        if (wp_attacks(st.wp_kind, st.bk, st.wp, st.wk)) return false;
        for (auto& m : generate_all_king_moves(st.bk)) {
            if (wp_attacks(st.wp_kind, m, st.wp, st.wk)) continue;
            if (m.distance_to(st.wk) <= 1) continue;
            return false;
        }
        return true;
    }

    string get_move_notation(const PawnState& from, const PawnState& to) const {
        if (from.wk != to.wk) return "K" + to.wk.str();
        if (from.wp != to.wp || from.wp_kind != to.wp_kind) {
            string mv(1, kind_letter(from.wp_kind));
            mv += to.wp.str();
            if (from.wp_kind != to.wp_kind) {
                mv += "=";
                mv += kind_letter(to.wp_kind);
            }
            return mv;
        }
        if (from.bk != to.bk) return "k" + to.bk.str();
        return "??";
    }

    // Inverse of get_move_notation. Handles a trailing "=Q"/"=N" promotion
    // suffix by updating wp_kind in addition to the destination square.
    PawnState apply_move_notation(const PawnState& from, const string& mv) const {
        PawnState result = from;
        if (mv.length() < 2) return result;
        char piece = mv[0];
        size_t eq = mv.find('=');
        string dest_str = (eq == string::npos) ? mv.substr(1) : mv.substr(1, eq - 1);
        Position dest = Position::from_str(dest_str);
        if (piece == 'K') {
            result.wk = dest; result.to_move = 'B';
        } else if (piece == 'k') {
            result.bk = dest; result.to_move = 'W';
        } else {
            result.wp = dest;
            if (eq != string::npos && eq + 1 < mv.length()) {
                result.wp_kind = kind_from_letter(mv[eq + 1]);
            }
            result.to_move = 'B';
        }
        return result;
    }

    // Only ever called with to_move=='B' throughout this file (mirroring
    // exactly how the existing engine's count_legal_moves is only ever
    // invoked for Black) -- Black's own escape-square count at this position.
    int count_legal_moves(const PawnState& st) const {
        int cnt = 0;
        for (auto& bk_n : generate_all_king_moves(st.bk)) {
            if (wp_attacks(st.wp_kind, bk_n, st.wp, st.wk)) continue;
            if (bk_n.distance_to(st.wk) <= 1) continue;
            PawnState ns = st; ns.bk = bk_n; ns.to_move = 'W';
            if (is_legal_state(ns)) cnt++;
        }
        return cnt;
    }

    vector<PawnState> generate_candidates(const PawnState& st) const {
        vector<PawnState> cands;
        cands.reserve(16);
        if (st.to_move == 'W') {
            for (auto& wk_n : generate_all_king_moves(st.wk)) {
                PawnState ns = st; ns.wk = wk_n; ns.to_move = 'B';
                if (is_legal_state(ns) && !is_stalemate(ns)) cands.push_back(ns);
            }

            if (st.wp_kind == PieceKind::PAWN) {
                for (auto& pm : generate_pawn_quiet_moves(st.wp, st.wk, st.bk)) {
                    // "Don't hang the piece" gate -- same reasoning as the
                    // existing engine, and now doubly load-bearing: per the
                    // fix from two runs back, hanging this piece is a
                    // provable permanent draw, so generating it here would
                    // hand White a candidate that's always strictly worse
                    // than avoiding it, for zero benefit.
                    if (pm.dest.distance_to(st.bk) < 2 && pm.dest.distance_to(st.wk) > 1) continue;

                    if (pm.is_promotion) {
                        // Only Queen and Knight are ever correct promotion
                        // choices. Any line achievable by promoting to Rook
                        // or Bishop is also achievable by promoting to Queen
                        // and simply choosing to only ever play the moves
                        // the Rook/Bishop would have played -- Queen weakly
                        // dominates both, for every possible continuation,
                        // under both stages of this engine's optimality
                        // criterion. Knight is NOT dominated (it reaches
                        // squares in one move a Queen cannot), so it's the
                        // only other real candidate worth searching.
                        for (PieceKind k : {PieceKind::QUEEN, PieceKind::KNIGHT}) {
                            PawnState ns = st; ns.wp = pm.dest; ns.wp_kind = k; ns.to_move = 'B';
                            if (is_legal_state(ns) && !is_stalemate(ns)) cands.push_back(ns);
                        }
                    } else {
                        PawnState ns = st; ns.wp = pm.dest; ns.to_move = 'B';
                        if (is_legal_state(ns) && !is_stalemate(ns)) cands.push_back(ns);
                    }
                }
            } else {
                vector<Position> piece_moves = (st.wp_kind == PieceKind::QUEEN)
                    ? QueenRules::generate_moves(st.wp)
                    : KnightRules::generate_moves(st.wp);
                for (auto& p_n : piece_moves) {
                    if (p_n.distance_to(st.bk) < 2 && p_n.distance_to(st.wk) > 1) continue;

                    if (st.wp_kind == PieceKind::QUEEN) {
                        bool blocked = false;
                        if (p_n.file == st.wp.file) {
                            int lo = min(st.wp.rank, p_n.rank) + 1, hi = max(st.wp.rank, p_n.rank);
                            for (int r = lo; r < hi; r++) if (Position(p_n.file, r) == st.wk) { blocked = true; break; }
                        } else if (p_n.rank == st.wp.rank) {
                            int lo = min(st.wp.file, p_n.file) + 1, hi = max(st.wp.file, p_n.file);
                            for (int f = lo; f < hi; f++) if (Position(f, p_n.rank) == st.wk) { blocked = true; break; }
                        } else if (abs(p_n.file - st.wp.file) == abs(p_n.rank - st.wp.rank)) {
                            int df = (p_n.file > st.wp.file) ? 1 : -1, dr = (p_n.rank > st.wp.rank) ? 1 : -1;
                            int f = st.wp.file + df, r = st.wp.rank + dr;
                            while (f != p_n.file) { if (Position(f, r) == st.wk) { blocked = true; break; } f += df; r += dr; }
                        }
                        if (blocked) continue;
                    }
                    // Knight never needs a blocking check -- it jumps over pieces.

                    PawnState ns = st; ns.wp = p_n; ns.to_move = 'B';
                    if (is_legal_state(ns) && !is_stalemate(ns)) cands.push_back(ns);
                }
            }
        } else {
            for (auto& bk_n : generate_all_king_moves(st.bk)) {
                if (wp_attacks(st.wp_kind, bk_n, st.wp, st.wk)) continue;
                if (bk_n.distance_to(st.wk) <= 1) continue;
                PawnState ns = st; ns.bk = bk_n; ns.to_move = 'W';
                if (!is_legal_state(ns)) continue;
                cands.push_back(ns);
            }
        }
        return cands;
    }

    uint64_t make_cache_key(const PawnState& st, int depth) const {
        uint64_t key = 0;
        key |= ((uint64_t)st.wk.file << 60);
        key |= ((uint64_t)st.wk.rank << 56);
        key |= ((uint64_t)st.wp.file << 52);
        key |= ((uint64_t)st.wp.rank << 48);
        key |= ((uint64_t)st.bk.file << 44);
        key |= ((uint64_t)st.bk.rank << 40);
        // wp_kind needs its own bits now that the piece's square alone no
        // longer determines what it is (see PawnState's comment on str()).
        // Placed at bit 10: comfortably above depth's realistic range
        // (bits 0-9, i.e. up to 1023 -- depth never remotely approaches
        // that in practice) and comfortably below the coordinate fields
        // starting at bit 40, so neither can collide with it.
        key |= ((uint64_t)static_cast<uint8_t>(st.wp_kind) << 10);
        key |= (uint64_t)depth;
        return key;
    }

    PawnSearchResult compositional_search_impl(
        const PawnState& st, int depth, int ply,
        unordered_map<uint64_t, PawnSearchResult>& memo
    ) {
        uint64_t cache_key = make_cache_key(st, depth);
        auto it = memo.find(cache_key);
        if (it != memo.end()) return it->second;

        if (is_checkmate(st)) {
            PawnSearchResult res{0, 0, nullopt, {}};
            memo[cache_key] = res;
            return res;
        }

        // Same fix as the existing engine, generalized to whichever piece
        // currently occupies wp: if Black's king is adjacent to it and
        // White's king doesn't defend it, capturing it is a completely
        // legal move, and for a lone-attacker endgame that always means
        // "White now has a bare king" -- an unconditional draw. Reported
        // unresolved exactly like running out of search depth, so a parent
        // that has this available never gets to claim a false fast mate
        // through it.
        if (st.to_move == 'B' && st.bk.distance_to(st.wp) <= 1 && st.wk.distance_to(st.wp) > 1) {
            PawnSearchResult res{nullopt, nullopt, nullopt, {}};
            memo[cache_key] = res;
            return res;
        }

        if (depth == 0) {
            PawnSearchResult res{nullopt, nullopt, nullopt, {}};
            memo[cache_key] = res;
            return res;
        }

        vector<PawnState> cands = generate_candidates(st);
        if (cands.empty()) {
            PawnSearchResult res{nullopt, nullopt, nullopt, {}};
            memo[cache_key] = res;
            return res;
        }

        string dir = (st.to_move == 'W') ? "minimize" : "maximize";
        optional<int> best_val, best_bn_cum;
        optional<PawnState> best_mv;
        bool all_candidates_resolved = true;
        vector<pair<int, pair<optional<int>, PawnState>>> resolved_candidates;
        resolved_candidates.reserve(cands.size());

        for (auto& c : cands) {
            PawnSearchResult rec = compositional_search_impl(c, depth - 1, ply + 1, memo);
            nodes_evaluated++;
            optional<int> val = rec.val;
            optional<int> child_bn_cum = rec.bn_cum;

            if (!val) {
                all_candidates_resolved = false;
                continue;
            }

            int v = *val + 1;
            int own_contribution = 0;
            if (c.to_move == 'B') own_contribution = count_legal_moves(c);
            optional<int> this_bn_cum;
            if (child_bn_cum) this_bn_cum = own_contribution + *child_bn_cum;

            resolved_candidates.push_back({v, {this_bn_cum, c}});

            if (!best_val) {
                best_val = v; best_bn_cum = this_bn_cum; best_mv = c;
            } else if (dir == "minimize" && v < *best_val) {
                best_val = v; best_bn_cum = this_bn_cum; best_mv = c;
            } else if (dir == "minimize" && v == *best_val) {
                if (this_bn_cum && best_bn_cum && *this_bn_cum < *best_bn_cum) {
                    best_val = v; best_bn_cum = this_bn_cum; best_mv = c;
                }
            } else if (dir == "maximize" && v > *best_val) {
                best_val = v; best_bn_cum = this_bn_cum; best_mv = c;
            } else if (dir == "maximize" && v == *best_val) {
                if (this_bn_cum && best_bn_cum && *this_bn_cum > *best_bn_cum) {
                    best_val = v; best_bn_cum = this_bn_cum; best_mv = c;
                }
            }
        }

        if (dir == "maximize" && !all_candidates_resolved) {
            PawnSearchResult res{nullopt, nullopt, nullopt, {}};
            memo[cache_key] = res;
            return res;
        }

        vector<pair<PawnState, int>> m_tied;
        if (best_val) {
            for (auto& [v, rest] : resolved_candidates) {
                if (v != *best_val) continue;
                auto& [bncum_opt, state] = rest;
                int bncum_val = bncum_opt ? *bncum_opt : -1;
                m_tied.emplace_back(state, bncum_val);
            }
        }
        vector<pair<PawnState, int>> tied;
        if (!m_tied.empty()) {
            int extremal = m_tied[0].second;
            for (auto& [state, bn] : m_tied) {
                extremal = (dir == "minimize") ? min(extremal, bn) : max(extremal, bn);
            }
            for (auto& [state, bn] : m_tied) {
                if (bn == extremal) tied.emplace_back(state, bn);
            }
        }

        PawnSearchResult res{best_val, best_bn_cum, best_mv, tied};
        memo[cache_key] = res;
        return res;
    }

    tuple<optional<PawnState>, optional<int>, vector<pair<PawnState, int>>> find_best_move(
        const PawnState& st, int max_depth = 30
    ) {
        nodes_evaluated = 0;
        unordered_map<uint64_t, PawnSearchResult> memo;
        optional<PawnState> best_move;
        optional<int> best_value;
        vector<pair<PawnState, int>> tied;

        for (int depth = 2; depth <= max_depth + 1; depth += 2) {
            PawnSearchResult r = compositional_search_impl(st, depth, 0, memo);
            if (r.val) {
                best_value = r.val; best_move = r.mv; tied = r.tied;
                break;
            }
        }
        return make_tuple(best_move, best_value, tied);
    }

    tuple<vector<string>, int, bool> play_complete_game(
        const PawnState& first, int max_moves = 50, int game_search_depth = 20
    ) {
        vector<string> mvs;
        PawnState curr = first;
        for (int move_num = 0; move_num < max_moves; move_num++) {
            if (is_checkmate(curr)) return make_tuple(mvs, (int)mvs.size(), true);
            auto [ns, val, tied] = find_best_move(curr, 2 * game_search_depth);
            if (!ns) return make_tuple(mvs, (int)mvs.size(), false);
            mvs.push_back(get_move_notation(curr, *ns));
            curr = *ns;
        }
        return make_tuple(mvs, (int)mvs.size(), false);
    }
};

// ============================================================================
// GeneralState / GeneralEngine -- arbitrary pieces, either color.
// ============================================================================
//
// This is the actual generalization: KQvK, KRvK, KPvK, and multi-attacker
// endgames like KBPvK are all instances of the SAME representation and the
// SAME search here, differing only in which pieces populate GeneralState's
// piece list at the start. It does not replace CompositionalEngine<AttackerRules>
// or PawnEngine -- both stay exactly as they were, so every existing
// kqvk_solver/krvk_solver/kpvk_solver binary keeps working unchanged and at
// zero added risk. This is verified, not just asserted: configured to
// White=King+Queen, Black=King only, this engine is tested to reproduce the
// existing Syzygy-verified CompositionalEngine<QueenRules> results exactly
// on the same positions used throughout this project, and likewise for Rook.
//
// Reuses every already-verified piece rule as a building block:
// QueenRules/RookRules/BishopRules::directions() for sliding pieces,
// KnightRules::generate_moves for the knight jump table, and the pawn
// quiet-move logic from the PawnEngine section above (generalized here to
// both colors, since Black can now actually have a pawn to capture with a
// White pawn's diagonal attack, unlike in pure KPvK where that branch was
// provably unreachable).
//
// Deliberately NOT implemented, and not silently pretended to work:
//   - Castling: requires move-history (has this king/rook ever moved), which
//     GeneralState doesn't track. Essentially never relevant to constructed
//     endgame positions (by the time material is this reduced, castling
//     rights are already gone in any realistic game), so this is scoped out
//     rather than built untested.
//   - En passant: same reasoning as the pawn engine -- needs move history,
//     only matters with two pawns adjacent on the same rank in a very
//     specific configuration. Can be added later following the same
//     optional-target-square design already described, once actually needed.
//   - Persisting these positions into the existing SolvedPositionDatabase/CSV
//     format with correct multi-piece symmetry expansion. add_position's
//     transform logic hardcodes exactly three roles (WK:/WQ:/BK:) and was
//     never designed for a variable-length piece list; getting symmetry
//     right for arbitrary combinations (which reflections are valid depends
//     on which colors have pawns, not just whether ANY piece is a pawn) is
//     its own real design problem. Rather than rush an unverified scheme in
//     alongside everything else here, this engine currently runs standalone
//     (find_best_move_general / play_complete_game_general), proven correct
//     in isolation. Wiring it into --full-dag-style persistent exploration
//     is a distinct next step once this core is trusted.

enum class Color : uint8_t { WHITE, BLACK };

struct GPiece {
    PieceKind kind;
    Color color;
    Position square;

    bool operator==(const GPiece& o) const {
        return kind == o.kind && color == o.color && square == o.square;
    }
};

class GeneralState {
public:
    vector<GPiece> pieces;
    char to_move;

    GeneralState() : to_move('W') {}

    const GPiece* piece_at(const Position& sq) const {
        for (auto& p : pieces) if (p.square == sq) return &p;
        return nullptr;
    }

    const GPiece* find_king(Color c) const {
        for (auto& p : pieces) if (p.kind == PieceKind::KING && p.color == c) return &p;
        return nullptr;
    }

    // Canonical, order-independent serialization: pieces sorted by color then
    // kind then square, so two GeneralStates representing the identical board
    // always produce the identical string regardless of internal vector order.
    string str() const {
        vector<GPiece> sorted = pieces;
        sort(sorted.begin(), sorted.end(), [](const GPiece& a, const GPiece& b) {
            if (a.color != b.color) return a.color < b.color;
            if (a.kind != b.kind) return a.kind < b.kind;
            return (a.square.file * 8 + a.square.rank) < (b.square.file * 8 + b.square.rank);
        });
        string s;
        for (auto& p : sorted) {
            if (!s.empty()) s += " ";
            char letter = kind_letter(p.kind);
            if (p.color == Color::BLACK) letter = tolower(letter);
            s += letter;
            s += ":";
            s += p.square.str();
        }
        return s;
    }

    bool operator==(const GeneralState& o) const {
        return pieces == o.pieces && to_move == o.to_move;
    }
};

// Does a piece of `kind`/`color` sitting at `from` attack `to`, given the
// full board for blocking? Reuses QueenRules/RookRules/BishopRules::directions()
// for sliding pieces and KnightRules' jump logic for knights; king and pawn
// attacks are simple enough to inline directly.
inline bool general_attacks(const GeneralState& gs, const Position& from, PieceKind kind, Color color, const Position& to) {
    if (from == to) return false;
    switch (kind) {
        case PieceKind::KING:
            return from.distance_to(to) == 1;
        case PieceKind::KNIGHT: {
            int df = abs(from.file - to.file), dr = abs(from.rank - to.rank);
            return (df == 1 && dr == 2) || (df == 2 && dr == 1);
        }
        case PieceKind::PAWN: {
            int forward = (color == Color::WHITE) ? 1 : -1;
            if (to.rank != from.rank + forward) return false;
            return abs(to.file - from.file) == 1;
        }
        case PieceKind::QUEEN:
        case PieceKind::ROOK:
        case PieceKind::BISHOP: {
            vector<pair<int,int>> dirs = (kind == PieceKind::QUEEN) ? QueenRules::directions()
                                        : (kind == PieceKind::ROOK)  ? RookRules::directions()
                                                                     : BishopRules::directions();
            for (auto& [df, dr] : dirs) {
                int f = from.file, r = from.rank;
                while (true) {
                    f += df; r += dr;
                    if (f < 0 || f > 7 || r < 0 || r > 7) break;
                    if (Position(f, r) == to) return true;
                    if (gs.piece_at(Position(f, r))) break;  // blocked before reaching `to`
                }
            }
            return false;
        }
    }
    return false;
}

inline bool square_attacked_by(const GeneralState& gs, const Position& sq, Color by) {
    for (auto& p : gs.pieces) {
        if (p.color != by) continue;
        if (general_attacks(gs, p.square, p.kind, p.color, sq)) return true;
    }
    return false;
}

inline bool is_in_check(const GeneralState& gs, Color c) {
    const GPiece* king = gs.find_king(c);
    if (!king) return false;  // shouldn't happen in any legal state, but never crash over it
    Color enemy = (c == Color::WHITE) ? Color::BLACK : Color::WHITE;
    return square_attacked_by(gs, king->square, enemy);
}

// Pseudo-legal destinations for one piece (occupancy-aware: sliding pieces
// stop at the first blocker and may capture it if it's an enemy; King/Knight
// filter out only squares occupied by an own piece; Pawn handles quiet
// forward moves and diagonal captures separately, since -- unlike every
// other piece here -- its move squares and attack squares are genuinely
// different sets). Does NOT check whether the move leaves the mover's own
// king in check -- that filtering happens once, generically, in
// generate_legal_moves below, rather than being duplicated per piece type.
inline vector<pair<Position,bool>> piece_destinations(const GeneralState& gs, const GPiece& p) {
    // second element of the pair: true if this destination is a promotion square (pawn only)
    //
    // CRITICAL: no branch below may ever generate a destination that lands
    // ON the enemy king's square. In real chess this can never come up --
    // checkmate ends the game before any piece could actually capture a
    // king -- but this function has no way to know "the game should already
    // be over" for an arbitrary constructed GeneralState, so it must
    // enforce the rule directly: a king is never a valid capture target,
    // for any piece, ever. Confirmed as a real, live bug: a sliding piece
    // was generating "capture the enemy king" as a legal move whenever a
    // constructed test position happened to already have the enemy king
    // under attack while it was the OTHER side's turn (an artificial
    // precondition used throughout this project's regression positions,
    // which the single-attacker engine never triggers because it has no
    // capture mechanic at all -- it only ever generates geometric
    // destinations, with "capture" implemented purely through checkmate
    // detection, never as a literal remove-the-piece move).
    vector<pair<Position,bool>> dests;
    switch (p.kind) {
        case PieceKind::KING: {
            for (int df = -1; df <= 1; df++) for (int dr = -1; dr <= 1; dr++) {
                if (df == 0 && dr == 0) continue;
                int nf = p.square.file + df, nr = p.square.rank + dr;
                if (nf < 0 || nf > 7 || nr < 0 || nr > 7) continue;
                Position d(nf, nr);
                const GPiece* occ = gs.piece_at(d);
                if (occ && occ->color == p.color) continue;
                if (occ && occ->kind == PieceKind::KING) continue;  // never capture a king
                dests.push_back({d, false});
            }
            break;
        }
        case PieceKind::KNIGHT: {
            int offs[][2] = {{1,2},{2,1},{2,-1},{1,-2},{-1,-2},{-2,-1},{-2,1},{-1,2}};
            for (auto& o : offs) {
                int nf = p.square.file + o[0], nr = p.square.rank + o[1];
                if (nf < 0 || nf > 7 || nr < 0 || nr > 7) continue;
                Position d(nf, nr);
                const GPiece* occ = gs.piece_at(d);
                if (occ && occ->color == p.color) continue;
                if (occ && occ->kind == PieceKind::KING) continue;  // never capture a king
                dests.push_back({d, false});
            }
            break;
        }
        case PieceKind::QUEEN:
        case PieceKind::ROOK:
        case PieceKind::BISHOP: {
            vector<pair<int,int>> dirs = (p.kind == PieceKind::QUEEN) ? QueenRules::directions()
                                        : (p.kind == PieceKind::ROOK)  ? RookRules::directions()
                                                                       : BishopRules::directions();
            for (auto& [df, dr] : dirs) {
                int f = p.square.file, r = p.square.rank;
                while (true) {
                    f += df; r += dr;
                    if (f < 0 || f > 7 || r < 0 || r > 7) break;
                    Position d(f, r);
                    const GPiece* occ = gs.piece_at(d);
                    if (!occ) { dests.push_back({d, false}); continue; }
                    // A real piece attacking the enemy king gives check --
                    // it does not get to move onto or past the king's
                    // square. Stop sliding here either way (occupied), but
                    // only add the destination if it's a capturable
                    // non-king enemy piece.
                    if (occ->color != p.color && occ->kind != PieceKind::KING) dests.push_back({d, false});
                    break;
                }
            }
            break;
        }
        case PieceKind::PAWN: {
            int forward = (p.color == Color::WHITE) ? 1 : -1;
            int start_rank = (p.color == Color::WHITE) ? 1 : 6;
            int promo_rank = (p.color == Color::WHITE) ? 7 : 0;

            int nr1 = p.square.rank + forward;
            if (nr1 >= 0 && nr1 <= 7) {
                Position one(p.square.file, nr1);
                if (!gs.piece_at(one)) {
                    dests.push_back({one, nr1 == promo_rank});
                    if (p.square.rank == start_rank && nr1 != promo_rank) {
                        Position two(p.square.file, p.square.rank + 2 * forward);
                        if (!gs.piece_at(two)) dests.push_back({two, false});
                    }
                }
                for (int df : {-1, 1}) {
                    int nf = p.square.file + df;
                    if (nf < 0 || nf > 7) continue;
                    Position diag(nf, nr1);
                    const GPiece* occ = gs.piece_at(diag);
                    if (occ && occ->color != p.color && occ->kind != PieceKind::KING) {
                        dests.push_back({diag, nr1 == promo_rank});
                    }
                }
            }
            break;
        }
    }
    return dests;
}

// Every pseudo-legal move for the side to move -- a promotion produces TWO
// candidates (Queen and Knight only, per the dominance argument already
// established: Rook/Bishop promotions are never better than Queen under
// this engine's exact optimality criterion, for either color).
inline vector<GeneralState> generate_pseudo_legal_moves(const GeneralState& gs) {
    vector<GeneralState> out;
    Color mover = (gs.to_move == 'W') ? Color::WHITE : Color::BLACK;
    for (size_t i = 0; i < gs.pieces.size(); i++) {
        if (gs.pieces[i].color != mover) continue;
        Position original_square = gs.pieces[i].square;  // stable identity, captured before any mutation
        for (auto& [dest, is_promo] : piece_destinations(gs, gs.pieces[i])) {
            auto make_child = [&](PieceKind resulting_kind) {
                GeneralState child = gs;
                // Remove any captured enemy piece at dest FIRST. This shifts
                // every element after the removed one down by one position --
                // which is exactly why re-indexing by the ORIGINAL `i` below
                // would be wrong (and was the actual bug here originally):
                // if the captured piece sat at a lower index than the mover,
                // `i` goes stale, silently reading/writing the wrong element
                // or past the end of the (now shorter) vector -- undefined
                // behavior, confirmed to produce exactly this kind of
                // corruption (a captured position being mis-evaluated as a
                // safe king move, an "own-king-safety" check on real garbage
                // that happened to look like a queen turning into a king).
                // Re-locating the mover by its ORIGINAL SQUARE instead (a
                // stable, unique identifier no erasure can invalidate) fixes
                // this regardless of which piece got captured or where.
                child.pieces.erase(
                    remove_if(child.pieces.begin(), child.pieces.end(),
                              [&](const GPiece& q) { return q.square == dest && q.color != mover; }),
                    child.pieces.end());
                for (auto& p : child.pieces) {
                    if (p.square == original_square && p.color == mover) {
                        p.square = dest;
                        p.kind = resulting_kind;
                        break;
                    }
                }
                child.to_move = (mover == Color::WHITE) ? 'B' : 'W';
                out.push_back(child);
            };
            if (is_promo) {
                make_child(PieceKind::QUEEN);
                make_child(PieceKind::KNIGHT);
            } else {
                make_child(gs.pieces[i].kind);
            }
        }
    }
    return out;
}

// Filters pseudo-legal moves down to legal ones: does NOT leave the mover's
// own king in check afterward. Implemented by literally trying the move and
// checking -- simpler and unambiguously correct versus hand-rolled pin
// detection, at the cost of a check-scan per candidate. Given this project's
// consistent priority on correctness over performance in new territory
// (explicitly agreed for the pawn engine too), that's the right tradeoff here.
inline vector<GeneralState> generate_legal_moves(const GeneralState& gs) {
    Color mover = (gs.to_move == 'W') ? Color::WHITE : Color::BLACK;
    vector<GeneralState> legal;
    for (auto& child : generate_pseudo_legal_moves(gs)) {
        if (!is_in_check(child, mover)) legal.push_back(child);
    }
    return legal;
}

inline bool is_checkmate_general(const GeneralState& gs) {
    Color mover = (gs.to_move == 'W') ? Color::WHITE : Color::BLACK;
    if (!is_in_check(gs, mover)) return false;
    return generate_legal_moves(gs).empty();
}

// Can this SET of piece kinds (position-independent -- just which types are
// present) ever force checkmate against a lone king, under best defense?
// This is a lookup against established chess theory, not something the
// search needs to discover -- and it's deliberately conservative: anything
// not clearly settled defaults to "assume yes" (no shortcut applied),
// falling back to the same slower-but-correct behavior as if this function
// didn't exist at all. A false "no" here would be a silent correctness bug
// (skipping a genuinely winnable position); a false "yes" only costs some
// wasted search time before correctly concluding "unresolved" -- so erring
// toward "yes" is the safe direction on both ends of any real uncertainty.
//
// Established facts this actually relies on:
//   - A lone Queen or Rook always suffices (already proven throughout this
//     project); a Pawn always suffices given time (it can promote into one).
//   - A lone Bishop or lone Knight can never mate alone (basic insufficient-
//     material rule).
//   - Two Knights alone cannot force mate against best defense (well-
//     established: KNNvK is a draw with correct defense, even though a
//     mate CAN occur against a mistake -- this engine searches for forced
//     mate under best play, so the correct classification here is "no").
//   - Bishop + Knight together CAN force mate (the classic "hard case" in
//     endgame theory, distinct from either piece alone).
//   - Two Bishops on opposite-colored squares CAN force mate; two Bishops
//     on the SAME colored squares cannot (they never cover the other
//     color's mating net) -- a case that essentially never arises from
//     real promotion (restricted to Queen/Knight here) but is included for
//     completeness if two literal Bishops are given as starting material.
inline bool material_can_mate(const vector<GPiece>& white_pieces) {
    if (white_pieces.empty()) return false;
    for (auto& p : white_pieces) {
        if (p.kind == PieceKind::QUEEN || p.kind == PieceKind::ROOK || p.kind == PieceKind::PAWN) {
            return true;
        }
    }
    // Everything remaining is Bishop and/or Knight only.
    int bishops = 0, knights = 0;
    bool bishop_light = false, bishop_dark = false;
    for (auto& p : white_pieces) {
        if (p.kind == PieceKind::KNIGHT) knights++;
        else if (p.kind == PieceKind::BISHOP) {
            bishops++;
            if ((p.square.file + p.square.rank) % 2 == 0) bishop_dark = true;
            else bishop_light = true;
        }
    }
    if (bishops + knights <= 1) return false;      // single minor piece: no
    if (bishops == 0) return false;                // 2+ knights, no bishop: no (best-defense draw)
    if (knights >= 1) return true;                  // bishop + knight together: yes
    if (bishop_light && bishop_dark) return true;   // opposite-colored bishops: yes
    if (bishops == 2) return false;                 // exactly 2 same-colored bishops: no
    return true;  // 3+ same-colored bishops or anything else unaccounted for: default to "assume yes"
}

inline bool is_stalemate_general(const GeneralState& gs) {
    Color mover = (gs.to_move == 'W') ? Color::WHITE : Color::BLACK;
    if (is_in_check(gs, mover)) return false;
    return generate_legal_moves(gs).empty();
}

// Generalizes "cumulative escape count" from "Black's king's legal move
// count" (well-defined when the king is Black's only piece) to "total legal
// moves available to the side under attack, across all its pieces" --
// explicitly a design decision, not something implicit in the original
// single-piece definition. Preserves the same spirit (how confined was the
// losing side at this decision point) while being meaningful once that side
// can have more than a king.
inline int defender_mobility(const GeneralState& gs) {
    return (int)generate_legal_moves(gs).size();
}

// Identifies which piece moved between `from` and `to` via set-difference
// matching on (kind, square), NOT positional index comparison. This matters
// because a capture calls vector::erase on the piece list, which shifts every
// subsequent element down by one index -- comparing from.pieces[i] to
// to.pieces[i] directly breaks the instant the captured piece's original
// index is lower than the mover's, misattributing the move to whichever
// piece happens to now share that index. Confirmed directly: King captures
// Bishop, with the Bishop listed before the King, produced the notation
// "bb2=K" -- claiming a Black bishop promoted into a King -- before this fix.
string get_move_notation_general(const GeneralState& from, const GeneralState& to) {
    Color mover = (from.to_move == 'W') ? Color::WHITE : Color::BLACK;

    vector<GPiece> from_mover, to_mover;
    for (auto& p : from.pieces) if (p.color == mover) from_mover.push_back(p);
    for (auto& p : to.pieces) if (p.color == mover) to_mover.push_back(p);

    // Every mover-color piece that DIDN'T move has an exact (kind, square)
    // match on both sides -- pair those off first, tracking usage to handle
    // any coincidental duplicates correctly. Whatever's left unmatched on
    // each side is the same piece's old and new state.
    vector<bool> used(from_mover.size(), false);
    GPiece new_piece{}; bool found_new = false;
    for (auto& tp : to_mover) {
        bool matched = false;
        for (size_t i = 0; i < from_mover.size(); i++) {
            if (!used[i] && from_mover[i].kind == tp.kind && from_mover[i].square == tp.square) {
                used[i] = true; matched = true; break;
            }
        }
        if (!matched) { new_piece = tp; found_new = true; break; }
    }
    if (!found_new) return "??";

    GPiece old_piece{}; bool found_old = false;
    for (size_t i = 0; i < from_mover.size(); i++) {
        if (!used[i]) { old_piece = from_mover[i]; found_old = true; break; }
    }
    if (!found_old) return "??";

    char letter = kind_letter(old_piece.kind);
    if (old_piece.color == Color::BLACK) letter = tolower(letter);
    string mv(1, letter);
    mv += new_piece.square.str();
    if (old_piece.kind != new_piece.kind) {
        mv += "=";
        mv += kind_letter(new_piece.kind);
    }
    return mv;
}

GeneralState apply_move_notation_general(const GeneralState& from, const string& mv) {
    // Re-derives the move by generating legal moves and matching notation --
    // simpler and less error-prone than hand-parsing which piece-list index
    // moved, given the piece list's order is already stable across a single
    // search tree but this function may be called independently of one.
    for (auto& child : generate_pseudo_legal_moves(from)) {
        if (get_move_notation_general(from, child) == mv) return child;
    }
    return from;
}

struct GeneralSearchResult {
    optional<int> val;
    optional<int> mobility_cum;
    optional<GeneralState> mv;
    vector<pair<GeneralState,int>> tied;
};

class GeneralEngine {
public:
    int nodes_evaluated = 0;

    uint64_t make_cache_key(const GeneralState& st, int depth) const {
        // Piece lists are small (this engine targets endgames, not full
        // boards) -- a simple order-independent hash over the canonical
        // string plus depth is more than adequate and avoids having to
        // design a new fixed-width bit-packed key for a variable-length
        // piece list.
        std::hash<string> h;
        uint64_t key = h(st.str() + st.to_move);
        key ^= ((uint64_t)depth << 1);
        return key;
    }

    GeneralSearchResult compositional_search_impl(
        const GeneralState& st, int depth, unordered_map<uint64_t, GeneralSearchResult>& memo
    ) {
        uint64_t cache_key = make_cache_key(st, depth);
        auto it = memo.find(cache_key);
        if (it != memo.end()) return it->second;

        if (is_checkmate_general(st)) {
            GeneralSearchResult res{0, 0, nullopt, {}};
            memo[cache_key] = res;
            return res;
        }

        if (depth == 0) {
            GeneralSearchResult res{nullopt, nullopt, nullopt, {}};
            memo[cache_key] = res;
            return res;
        }

        // material_can_mate is checked here, not just once at the root:
        // material changes DURING search (captures remove pieces, pawns
        // promote), so a branch that starts winnable can become insufficient
        // partway through, and vice versa is impossible (material never
        // increases) -- meaning this check has to run at every node to
        // correctly short-circuit branches that just became hopeless via a
        // capture, not only ones that started that way. Directly measured to
        // matter: a real KBvK-shaped hopeless branch cost 17s and 30.6M
        // nodes to correctly resolve as unresolved WITHOUT this check.
        {
            vector<GPiece> white_pieces;
            for (auto& p : st.pieces) if (p.color == Color::WHITE && p.kind != PieceKind::KING) white_pieces.push_back(p);
            if (!material_can_mate(white_pieces)) {
                GeneralSearchResult res{nullopt, nullopt, nullopt, {}};
                memo[cache_key] = res;
                return res;
            }
        }

        vector<GeneralState> cands = generate_legal_moves(st);
        if (cands.empty()) {
            // Stalemate (not checkmate, since that's already handled above)
            GeneralSearchResult res{nullopt, nullopt, nullopt, {}};
            memo[cache_key] = res;
            return res;
        }

        Color mover = (st.to_move == 'W') ? Color::WHITE : Color::BLACK;
        bool white_to_move = (mover == Color::WHITE);
        optional<int> best_val, best_mob;
        optional<GeneralState> best_mv;
        bool all_resolved = true;
        vector<pair<int, pair<optional<int>, GeneralState>>> resolved;
        resolved.reserve(cands.size());

        for (auto& c : cands) {
            GeneralSearchResult rec = compositional_search_impl(c, depth - 1, memo);
            nodes_evaluated++;
            if (!rec.val) { all_resolved = false; continue; }

            int v = *rec.val + 1;
            // Only accumulate mobility on the DEFENDING side's own turn --
            // exactly the same gate the original single-attacker engines use
            // (own_contribution = count_legal_moves(c) only if c.to_move ==
            // 'B'). White is always the attacker/minimizer and Black is
            // always the defender/maximizer throughout this whole project,
            // including here -- this sum measures how confined the DEFENDER
            // was, not a running total of everyone's mobility at every node.
            int own_contribution = (c.to_move == 'B') ? defender_mobility(c) : 0;
            optional<int> this_mob;
            if (rec.mobility_cum) this_mob = own_contribution + *rec.mobility_cum;

            resolved.push_back({v, {this_mob, c}});

            if (!best_val) {
                best_val = v; best_mob = this_mob; best_mv = c;
            } else if (white_to_move && v < *best_val) {
                best_val = v; best_mob = this_mob; best_mv = c;
            } else if (white_to_move && v == *best_val) {
                if (this_mob && best_mob && *this_mob < *best_mob) { best_val = v; best_mob = this_mob; best_mv = c; }
            } else if (!white_to_move && v > *best_val) {
                best_val = v; best_mob = this_mob; best_mv = c;
            } else if (!white_to_move && v == *best_val) {
                if (this_mob && best_mob && *this_mob > *best_mob) { best_val = v; best_mob = this_mob; best_mv = c; }
            }
        }

        if (!white_to_move && !all_resolved) {
            GeneralSearchResult res{nullopt, nullopt, nullopt, {}};
            memo[cache_key] = res;
            return res;
        }

        vector<pair<GeneralState,int>> tied;
        if (best_val) {
            int extremal = 0;
            bool have_extremal = false;
            for (auto& [v, rest] : resolved) {
                if (v != *best_val) continue;
                int m = rest.first ? *rest.first : -1;
                if (!have_extremal) { extremal = m; have_extremal = true; }
                else extremal = white_to_move ? min(extremal, m) : max(extremal, m);
            }
            for (auto& [v, rest] : resolved) {
                if (v != *best_val) continue;
                int m = rest.first ? *rest.first : -1;
                if (m == extremal) tied.emplace_back(rest.second, m);
            }
        }

        GeneralSearchResult res{best_val, best_mob, best_mv, tied};
        memo[cache_key] = res;
        return res;
    }

    tuple<optional<GeneralState>, optional<int>, vector<pair<GeneralState,int>>> find_best_move(
        const GeneralState& st, int max_depth = 30
    ) {
        nodes_evaluated = 0;
        unordered_map<uint64_t, GeneralSearchResult> memo;
        optional<GeneralState> best_move;
        optional<int> best_value;
        vector<pair<GeneralState,int>> tied;
        for (int depth = 2; depth <= max_depth + 1; depth += 2) {
            GeneralSearchResult r = compositional_search_impl(st, depth, memo);
            if (r.val) { best_value = r.val; best_move = r.mv; tied = r.tied; break; }
        }
        return make_tuple(best_move, best_value, tied);
    }

    tuple<vector<string>, int, bool> play_complete_game(const GeneralState& first, int max_moves = 50, int depth_budget = 20) {
        vector<string> mvs;
        GeneralState curr = first;
        for (int i = 0; i < max_moves; i++) {
            if (is_checkmate_general(curr)) return make_tuple(mvs, (int)mvs.size(), true);
            auto [ns, val, tied] = find_best_move(curr, 2 * depth_budget);
            if (!ns) return make_tuple(mvs, (int)mvs.size(), false);
            mvs.push_back(get_move_notation_general(curr, *ns));
            curr = *ns;
        }
        return make_tuple(mvs, (int)mvs.size(), false);
    }
};

// ============================================================================
// RetrogradeClassifier -- genuine exhaustive fixed-point classification.
// ============================================================================
//
// This replaces the top-down recursive search's approach to draws entirely,
// for the reason worked out at length in this project's history: a top-down
// search with "on-stack" cycle detection can only ever prove a position
// UNRESOLVED relative to information that is itself still mid-computation --
// three separate attempts at patching that in built three separate real
// bugs (unbounded recursion exhausting the state space before finding a
// cycle; a "sound" version so conservative it lost all transposition reuse
// and became 22x slower; a hybrid version that silently allowed a
// suboptimal Black defense to be reported as forced mate). None of those
// bugs are possible here, because nothing is EVER marked classified until
// it is fully justified by facts that are THEMSELVES already fully
// justified -- exactly how real retrograde tablebase generation works.
//
// Algorithm: discover every position reachable from the given roots via
// EVERY legal move for both sides (not just optimal ones -- both sides'
// full alternative sets have to be accounted for, or a "forced win" claim
// isn't actually proven). Seed checkmate positions as distance 0. Then
// repeatedly pass over every unclassified position: a WHITE (minimizing)
// position becomes classified the moment ANY of its children is already
// classified, with distance = 1 + that child's minimum; a BLACK
// (maximizing) position becomes classified only once ALL of its children
// are classified, with distance = 1 + their maximum (Black's only real
// choice is which of several forced losses to pick, so the WORST for
// White is correct). Repeat to a fixed point. Anything left unclassified
// is a proven, permanent draw -- not because any single path looked
// drawish, but because the complete process considered every alternative
// for both sides and none of them ever found a way out.
//
// Two implementation details were each found to matter through direct
// testing, not by inspection alone:
//   1. Each pass's newly-classified results are staged separately and
//      merged into the shared classified map only AFTER the full pass
//      completes. Merging immediately (checked directly, and confirmed
//      with a real example) lets an unordered_map's arbitrary iteration
//      order cause a parent to miss a sibling that becomes classifiable
//      later in the SAME pass, silently dropping a genuine tied move even
//      though the distance itself stays correct.
//   2. Positions are packed into a single 8-byte integer (kind+color+square
//      per piece, up to 2 non-king White pieces plus both kings, 41 bits
//      total) rather than a heap-allocated string. A string-keyed version
//      of this exact algorithm was confirmed correct on KQvK/KRvK but was
//      OOM-killed partway through KBPvK's larger state space; the packed
//      version is dramatically more memory-efficient per node, though even
//      it has a real, measured ceiling -- see material_state_space_size_ok
//      below.
//
// A separate, genuinely pre-existing bug was found and fixed while
// building this: piece_destinations allowed a sliding piece, a knight, or
// a pawn to "capture" the enemy king outright. In real chess this can
// never come up (checkmate ends the game first), but this project's
// constructed regression positions -- used successfully throughout with
// the single-attacker engine -- include ones where the enemy king is
// already under attack while it happens to be the other side's turn (an
// artificial but harmless precondition for that engine, which has no
// capture mechanic at all). GeneralState's move generator, unlike the
// single-attacker engine, does model captures literally, and took that
// artificial precondition as license to generate "capture the king" as a
// legal move -- producing states missing a king entirely and, once packed,
// wrongly reconstructed with a phantom king on a1 during unpacking.
// Confirmed directly: this produced a false "forced mate in 7" and a false
// "forced mate in 1" on two of this project's own long-standing regression
// positions before the fix. Fixed at the root: no piece may ever generate
// the enemy king's square as a valid destination, matching the real chess
// rule that attacking a king gives check rather than permitting capture.

constexpr uint8_t GSTATE_NO_PIECE_SENTINEL = 7;

// Packs WK + up to 2 additional White pieces + BK + turn into a uint64_t:
// 10 bits/piece (3 kind, 6 square, 1 spare) x 4 slots + 1 turn bit = 41 bits.
// Reversible -- states are reconstructed on demand, never stored directly,
// which is most of where the memory savings over a string key come from.
// Currently scoped to White having 0-2 non-king pieces and Black having
// only a king, matching this project's current explicit scope (multiple
// White pieces first, multiple Black pieces as a later, separate step).
inline uint64_t pack_general_state(const GeneralState& st) {
    Position wk{}, bk{};
    vector<GPiece> white_others;
    for (auto& p : st.pieces) {
        if (p.kind == PieceKind::KING && p.color == Color::WHITE) wk = p.square;
        else if (p.kind == PieceKind::KING && p.color == Color::BLACK) bk = p.square;
        else white_others.push_back(p);
    }
    sort(white_others.begin(), white_others.end(), [](const GPiece& a, const GPiece& b) {
        if (a.kind != b.kind) return (int)a.kind < (int)b.kind;
        return (a.square.file * 8 + a.square.rank) < (b.square.file * 8 + b.square.rank);
    });
    if (white_others.size() > 2) {
        // This packing reserves exactly 2 slots for White's non-king pieces
        // -- silently packing only the first 2 and dropping the rest would
        // be a real, silent data-corruption bug (two different game states
        // collapsing onto the same key). Failing loudly here is deliberate:
        // this project's whole history says a silent wrong answer is far
        // worse than a crash pointing at exactly what's unsupported.
        cerr << "FATAL: pack_general_state received " << white_others.size()
             << " non-king White pieces, but only 2 are supported by this "
             << "packed representation. State: " << st.str() << "\n";
        abort();
    }
    while (white_others.size() < 2) {
        white_others.push_back({(PieceKind)GSTATE_NO_PIECE_SENTINEL, Color::WHITE, Position(0, 0)});
    }
    uint64_t key = 0;
    auto pack_piece = [&](int shift, uint8_t kind, const Position& sq) {
        key |= ((uint64_t)kind << shift);
        key |= ((uint64_t)(sq.file * 8 + sq.rank) << (shift + 3));
    };
    pack_piece(0, 0, wk);
    pack_piece(9, (uint8_t)white_others[0].kind, white_others[0].square);
    pack_piece(18, (uint8_t)white_others[1].kind, white_others[1].square);
    pack_piece(27, 0, bk);
    key |= ((uint64_t)(st.to_move == 'W' ? 1 : 0) << 36);
    return key;
}

inline GeneralState unpack_general_state(uint64_t key) {
    GeneralState st;
    auto unpack_piece = [&](int shift) -> pair<uint8_t, Position> {
        uint8_t kind = (key >> shift) & 0x7;
        uint8_t sq = (key >> (shift + 3)) & 0x3F;
        return {kind, Position(sq / 8, sq % 8)};
    };
    auto [wk_kind, wk_sq] = unpack_piece(0); (void)wk_kind;
    auto [p1_kind, p1_sq] = unpack_piece(9);
    auto [p2_kind, p2_sq] = unpack_piece(18);
    auto [bk_kind, bk_sq] = unpack_piece(27); (void)bk_kind;
    st.pieces.push_back({PieceKind::KING, Color::WHITE, wk_sq});
    if (p1_kind != GSTATE_NO_PIECE_SENTINEL) st.pieces.push_back({(PieceKind)p1_kind, Color::WHITE, p1_sq});
    if (p2_kind != GSTATE_NO_PIECE_SENTINEL) st.pieces.push_back({(PieceKind)p2_kind, Color::WHITE, p2_sq});
    st.pieces.push_back({PieceKind::KING, Color::BLACK, bk_sq});
    st.to_move = ((key >> 36) & 1) ? 'W' : 'B';
    return st;
}

struct RetroNode {
    vector<uint64_t> children;
    uint8_t flags = 0;  // bit0: checkmate, bit1: terminal draw (stalemate/insufficient material), bit2: white to move
};

struct RetroResult {
    int distance;
    int escape_cum;
    uint64_t best_child;
};

class RetrogradeClassifier {
public:
    unordered_map<uint64_t, RetroNode> nodes;
    unordered_map<uint64_t, RetroResult> classified;
    int passes_run = 0;

    // Measured directly, not estimated: KQvK and KRvK's full reachable
    // graphs (~370-400K states) use well under 200MB with this packed
    // representation. KBPvK's full graph (two non-king White pieces) was
    // still climbing past 21 million discovered states -- and past 3.9GB
    // of actual memory -- when this sandbox's OOM killer stopped it, with
    // no sign of being close to finished. This is a genuine memory-capacity
    // ceiling of this specific environment, not a flaw in the algorithm
    // (which is exact-match verified against Syzygy-backed ground truth on
    // every position it was able to run to completion on). Discovery below
    // checks against this cap and stops cleanly rather than letting the
    // OS kill the process outright.
    long long max_nodes_before_abort = 15000000;
    bool aborted_for_memory = false;

    void discover(const vector<GeneralState>& roots) {
        deque<uint64_t> queue;
        for (auto& r : roots) queue.push_back(pack_general_state(r));
        while (!queue.empty()) {
            if ((long long)nodes.size() >= max_nodes_before_abort) {
                aborted_for_memory = true;
                return;
            }
            uint64_t key = queue.front(); queue.pop_front();
            if (nodes.count(key)) continue;
            GeneralState st = unpack_general_state(key);
            RetroNode node;
            if (st.to_move == 'W') node.flags |= 4;
            if (is_checkmate_general(st)) { node.flags |= 1; nodes[key] = node; continue; }
            vector<GPiece> white_pieces;
            for (auto& p : st.pieces) if (p.color == Color::WHITE && p.kind != PieceKind::KING) white_pieces.push_back(p);
            if (!material_can_mate(white_pieces)) { node.flags |= 2; nodes[key] = node; continue; }
            vector<GeneralState> cands = generate_legal_moves(st);
            if (cands.empty()) { node.flags |= 2; nodes[key] = node; continue; }
            node.children.reserve(cands.size());
            for (auto& c : cands) {
                uint64_t ckey = pack_general_state(c);
                node.children.push_back(ckey);
                if (!nodes.count(ckey)) queue.push_back(ckey);
            }
            nodes[key] = node;
        }
    }

    void classify() {
        for (auto& [key, node] : nodes) if (node.flags & 1) classified[key] = {0, 0, 0};
        passes_run = 0;
        bool changed = true;
        while (changed) {
            changed = false;
            passes_run++;
            // Staged separately, merged only after the pass completes -- see
            // the large comment above this class for why this matters.
            unordered_map<uint64_t, RetroResult> newly_classified;
            for (auto& [key, node] : nodes) {
                if (classified.count(key)) continue;
                if (node.flags & 1 || node.flags & 2) continue;
                bool white_to_move = (node.flags & 4);
                if (white_to_move) {
                    optional<int> best_v, best_esc; uint64_t best_child = 0;
                    for (auto& ckey : node.children) {
                        auto cit = classified.find(ckey);
                        if (cit == classified.end()) continue;
                        int v = cit->second.distance + 1;
                        auto& cnode = nodes[ckey];
                        int contribution = ((cnode.flags & 4) == 0 ? (int)cnode.children.size() : 0) + cit->second.escape_cum;
                        if (!best_v || v < *best_v || (v == *best_v && contribution < *best_esc)) {
                            best_v = v; best_esc = contribution; best_child = ckey;
                        }
                    }
                    if (best_v) {
                        newly_classified[key] = {*best_v, *best_esc, best_child};
                        changed = true;
                    }
                } else {
                    if (node.children.empty()) continue;  // stalemate, already flagged above; defensive only
                    bool all_done = true;
                    for (auto& ckey : node.children) if (!classified.count(ckey)) { all_done = false; break; }
                    if (!all_done) continue;
                    int worst_v = -1;
                    for (auto& ckey : node.children) worst_v = max(worst_v, classified[ckey].distance + 1);
                    optional<int> best_esc; uint64_t best_child = 0;
                    for (auto& ckey : node.children) {
                        auto& cres = classified[ckey];
                        int v = cres.distance + 1;
                        if (v != worst_v) continue;
                        auto& cnode = nodes[ckey];
                        int contribution = ((cnode.flags & 4) == 0 ? (int)cnode.children.size() : 0) + cres.escape_cum;
                        if (!best_esc || contribution > *best_esc) { best_esc = contribution; best_child = ckey; }
                    }
                    newly_classified[key] = {worst_v, *best_esc, best_child};
                    changed = true;
                }
            }
            for (auto& [key, res] : newly_classified) classified[key] = res;
        }
    }

    optional<RetroResult> get(uint64_t key) const {
        auto it = classified.find(key);
        if (it == classified.end()) return nullopt;
        return it->second;
    }

    // Full tied-move set: every child matching both the recorded distance
    // AND the extremal escape count, matching this project's established
    // two-stage lexicographic optimality criterion exactly.
    vector<pair<GeneralState,int>> tied_for(uint64_t key) const {
        vector<pair<GeneralState,int>> out;
        auto it = classified.find(key);
        if (it == classified.end()) return out;
        auto& node = nodes.at(key);
        for (auto& ckey : node.children) {
            auto cit = classified.find(ckey);
            if (cit == classified.end()) continue;
            int v = cit->second.distance + 1;
            if (v != it->second.distance) continue;
            auto& cnode = nodes.at(ckey);
            int contribution = ((cnode.flags & 4) == 0 ? (int)cnode.children.size() : 0) + cit->second.escape_cum;
            if (contribution == it->second.escape_cum) out.push_back({unpack_general_state(ckey), contribution});
        }
        return out;
    }
};

// ============================================================================
// Sweep/persistence driver for RetrogradeClassifier, mirroring the
// --full-dag pattern used by the single-attacker engines and PawnEngine.
// ============================================================================
//
// Positions file format is necessarily different from the single-attacker
// engines' "DTZ,WK:.. WQ:.. BK:.." lines: material here is genuinely
// variable (which piece types White has is DATA, not a compile-time
// template parameter), so a root position has to spell out every piece.
// Format, one root per line after a header: a comma-separated list of
// "kind:square" tokens covering White's king, White's other piece(s), and
// Black's king -- e.g. "K:e6,B:f5,P:e5,k:h8" for a KBPvK start. Uppercase
// letters are White, lowercase Black, matching GeneralState::str()'s own
// convention. Root positions are always White to move.
inline optional<GeneralState> parse_general_position(const string& line) {
    GeneralState st;
    st.to_move = 'W';
    stringstream ss(line);
    string token;
    while (getline(ss, token, ',')) {
        // trim whitespace
        size_t start = token.find_first_not_of(" \t\r\n");
        size_t end = token.find_last_not_of(" \t\r\n");
        if (start == string::npos) continue;
        token = token.substr(start, end - start + 1);
        if (token.empty()) continue;
        size_t colon = token.find(':');
        if (colon == string::npos || colon != 1) return nullopt;
        char letter = token[0];
        string sq_str = token.substr(2);
        if (sq_str.length() != 2) return nullopt;
        Position sq;
        try { sq = Position::from_str(sq_str); } catch (...) { return nullopt; }
        Color color = isupper(letter) ? Color::WHITE : Color::BLACK;
        char upper = toupper(letter);
        PieceKind kind;
        if (upper == 'K') kind = PieceKind::KING;
        else if (upper == 'Q') kind = PieceKind::QUEEN;
        else if (upper == 'R') kind = PieceKind::ROOK;
        else if (upper == 'B') kind = PieceKind::BISHOP;
        else if (upper == 'N') kind = PieceKind::KNIGHT;
        else if (upper == 'P') kind = PieceKind::PAWN;
        else return nullopt;
        st.pieces.push_back({kind, color, sq});
    }
    if (!st.find_king(Color::WHITE) || !st.find_king(Color::BLACK)) return nullopt;
    int white_non_king = 0;
    for (auto& p : st.pieces) if (p.color == Color::WHITE && p.kind != PieceKind::KING) white_non_king++;
    if (white_non_king > 2) {
        // Caught here, at parse time, rather than only in pack_general_state's
        // hard abort() deep inside a run -- a malformed input file should
        // fail with a clear, immediate message pointing at the exact line,
        // not crash the process after minutes of discovery work.
        cerr << "  Rejected (only up to 2 non-king White pieces are supported): " << line << "\n";
        return nullopt;
    }
    return st;
}

inline vector<GeneralState> load_general_positions_from_file(const string& filename) {
    vector<GeneralState> positions;
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
        if (is_header) { is_header = false; continue; }
        if (line.empty()) continue;
        auto st = parse_general_position(line);
        if (!st) { cerr << "  Line " << line_num << " failed to parse: " << line << "\n"; continue; }
        positions.push_back(*st);
    }
    cout << "  Loaded " << positions.size() << " positions from " << filename << "\n";
    return positions;
}

// Simple, standalone CSV format for general/multi-piece results -- not
// SolvedPosition, since that format's add_position hardcodes exactly three
// roles (WK:/WQ:/BK:) and was never designed for a variable-length piece
// list. Symmetry expansion is deliberately not applied here yet, matching
// the scoping decision already noted where GeneralState was first defined:
// getting symmetry right for arbitrary material combinations (which
// reflections are valid depends on which colors have pawns, not just
// whether any piece is a pawn) is its own real design problem, and every
// row written here is still exact and independently correct without it --
// there just isn't yet a multiplier reducing how many positions have to be
// solved directly.
inline void export_general_classification(RetrogradeClassifier& rc, const string& db_file) {
    ofstream out(db_file);
    out << "Position|Turn|BestMove|Distance|EscapeCum|TiedMoves\n";
    long long written = 0;
    for (auto& [key, node] : rc.nodes) {
        auto res = rc.get(key);
        if (!res) continue;  // unclassified = proven draw, deliberately not written as a row
        GeneralState st = unpack_general_state(key);
        auto tied = rc.tied_for(key);
        string best_move;
        stringstream tied_ss;
        bool first = true;
        for (auto& [cstate, esc] : tied) {
            string mv = get_move_notation_general(st, cstate);
            if (best_move.empty()) best_move = mv;
            if (!first) tied_ss << ";";
            tied_ss << mv << ":" << esc;
            first = false;
        }
        out << st.str() << "|" << st.to_move << "|" << best_move << "|"
            << res->distance << "|" << res->escape_cum << "|" << tied_ss.str() << "\n";
        written++;
    }
    out.close();
    cout << "Wrote " << written << " classified positions to " << db_file << "\n";
}

inline void run_general_sweep(const vector<GeneralState>& roots, const string& db_file,
                               long long max_nodes = 15000000) {
    cout << "\n" << string(80, '=') << "\n";
    cout << "RETROGRADE CLASSIFICATION SWEEP (multi-piece)\n";
    cout << string(80, '=') << "\n\n";

    RetrogradeClassifier rc;
    rc.max_nodes_before_abort = max_nodes;

    auto start = chrono::high_resolution_clock::now();
    rc.discover(roots);
    auto disc_end = chrono::high_resolution_clock::now();
    cout << "Discovery: " << rc.nodes.size() << " positions in "
         << fixed << setprecision(1) << chrono::duration<double>(disc_end - start).count() << "s\n";

    if (rc.aborted_for_memory) {
        cerr << "\nERROR: discovery exceeded " << max_nodes << " positions and was stopped\n"
             << "before running out of actual memory. This material's full reachable\n"
             << "graph is larger than this run's safety cap -- raise max_nodes if you\n"
             << "have the RAM to back it, or reduce scope (fewer/simpler root positions,\n"
             << "material with fewer non-king pieces).\n";
        return;
    }

    rc.classify();
    auto class_end = chrono::high_resolution_clock::now();
    cout << "Classification: " << rc.passes_run << " passes, "
         << rc.classified.size() << "/" << rc.nodes.size() << " positions classified as forced wins in "
         << fixed << setprecision(1) << chrono::duration<double>(class_end - disc_end).count() << "s\n";
    cout << "(" << (rc.nodes.size() - rc.classified.size()) << " positions proven drawn)\n\n";

    for (auto& r : roots) {
        uint64_t rkey = pack_general_state(r);
        auto res = rc.get(rkey);
        cout << r.str() << " " << r.to_move << ": ";
        if (res) cout << "WIN in " << res->distance << " plies\n";
        else cout << "DRAW\n";
    }

    export_general_classification(rc, db_file);

    auto total_end = chrono::high_resolution_clock::now();
    cout << "\nTotal time: " << fixed << setprecision(1)
         << chrono::duration<double>(total_end - start).count() << "s\n";
}

template<typename AttackerRules>
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
    
    // Thin forwarding wrappers to AttackerRules -- kept under their original
    // names so every existing call site in this file (there are many) needs
    // no changes at all. The actual logic now lives in QueenRules (or
    // whichever Rules struct this engine is instantiated with) above.
    inline vector<Position> generate_all_queen_moves(const Position& pos) const {
        return AttackerRules::generate_moves(pos);
    }

    bool is_attacked_by_queen(const Position& pos, const Position& qp, const Position& wk, const Position& wq) const {
        return AttackerRules::attacks(pos, qp, wk, wq);
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
        if (from.wq != to.wq) return string(1, AttackerRules::letter) + to.wq.str();
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
        else if (piece == AttackerRules::letter) { result.wq = dest; result.to_move = 'B'; }
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

    // The full set of moves this engine considers PERFECT PLAY from this position,
    // under a two-stage definition: (1) every RESOLVED candidate whose v exactly
    // equals val (the game-theoretic optimum -- mate distance), THEN (2) among
    // those, only the ones ALSO achieving the extremal cumulative black-escape
    // count (bn_cum) -- minimum for White, maximum for Black. A move tied on M
    // alone but with a worse bn_cum than another M-tied move is NOT perfect play
    // by this definition and will not appear here, even though it does resolve
    // to the same mate distance. mv above is always one entry of this set (the
    // same one the pre-existing tie-break logic already selected), so nothing
    // that reads val/bn_cum/mv changes behavior -- this field is what to use for
    // "give me every move that's actually perfect from here", not val/mv alone.
    // Unresolved candidates are never in this list -- see the comment above the
    // resolution loop in compositional_search_impl for why an unresolved
    // candidate can never tie val, for either side to move.
    // Empty for terminal (checkmate) and unresolved (nullopt val) results.
    vector<pair<GameState, int>> tied;

    // NOTE: a legal_move_count/choice_spread pair was tried here and removed.
    // Reading them off resolved_candidates (the same early-stopping that's
    // provably safe for val/tied) is NOT safe for a count or a worst-case
    // value: find_best_move stops at the first depth where val resolves,
    // which can leave most OTHER legal candidates simply not yet resolved
    // even though they don't affect val at all. Measured directly on
    // WK:b1 WQ:d8 BK:c4: only 4 of White's 23 legal candidates had resolved
    // at the depth the search stopped at, giving legal_move_count=4 -- wrong
    // by a factor of nearly 6. Computing it correctly requires fully
    // resolving every legal candidate independently at the point of
    // recording a position, which is a fundamentally different (and much
    // more expensive) operation than anything this struct can provide as a
    // byproduct of the existing search.
};

template<typename AttackerRules>
class CompositionalEngine : public BaseEngine<AttackerRules> {
public:
    // Required because BaseEngine<AttackerRules> is a DEPENDENT base class:
    // the compiler's first lookup phase won't find unqualified calls to
    // inherited members without these. Every BaseEngine method actually
    // called anywhere below is listed here.
    using BaseEngine<AttackerRules>::generate_all_king_moves;
    using BaseEngine<AttackerRules>::generate_all_queen_moves;
    using BaseEngine<AttackerRules>::is_attacked_by_queen;
    using BaseEngine<AttackerRules>::is_legal_state;
    using BaseEngine<AttackerRules>::is_checkmate;
    using BaseEngine<AttackerRules>::is_stalemate;
    using BaseEngine<AttackerRules>::get_move_notation;
    using BaseEngine<AttackerRules>::count_legal_moves;

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

    // Extracted from compositional_search_impl's own candidate-generation logic
    // so there's exactly one implementation of "what counts as a legal
    // candidate from this position", not a risk of hand-maintained copies
    // drifting apart.
    vector<GameState> generate_candidates(const GameState& st) const {
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
        return cands;
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

        // PRE-EXISTING REPRESENTATIONAL GAP, now actually closed: GameState has no
        // way to represent "White's second piece is gone" (wq is always some real
        // square), which is why the original code could only guard against WHITE
        // ever voluntarily hanging it (the "hang the piece" comment in
        // generate_candidates). What that comment didn't cover: if it's Black's
        // move and Black's king is already adjacent to an undefended wq (whether
        // because White hung it or simply because Black's king walked there over
        // several prior moves -- the very common case in KRvK's standard box/ladder
        // technique, where the rook sits still on a file/rank for many plies),
        // capturing it is a completely legal chess move. generate_candidates'
        // Black branch builds GameState(wk, wq, bk_n=wq, 'W') for that candidate --
        // two pieces on one square -- and is_legal_state correctly rejects it as an
        // invalid STATE, for the wrong reason: it's not an illegal state, it's a
        // legal capture the data structure can't express, so the candidate silently
        // vanished instead of being resolved as "White now has a bare king, which
        // can never force mate." Confirmed empirically: WK:f3 WQ:c2 BK:d1 (Black to
        // move) has Black's king adjacent to an undefended rook on c2, and this gap
        // caused the position -- and everything above it in the search -- to be
        // scored as a sound mate-in-9 when Black actually had a permanent escape.
        // Once White has only a king, KvK is an unconditional draw, so this is
        // reported unresolved exactly like running out of search depth: Black
        // (maximizing) correctly treats an unresolved option as at least as good as
        // anything resolved, so a parent that has this available never gets to
        // claim a fast forced mate through it.
        if (st.to_move == 'B' && st.bk.distance_to(st.wq) <= 1 && st.wk.distance_to(st.wq) > 1) {
            SearchResult res{nullopt, nullopt, nullopt};
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

        vector<GameState> cands = generate_candidates(st);

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

        // TWO-STAGE definition of "tied for perfect play", per the actual
        // definition of optimal play in this engine: (1) exactly optimal mate
        // distance (M), AND (2) among those, exactly optimal cumulative
        // black-escape count over the rest of the game -- minimum for White
        // (most confining), maximum for Black (most resistant). A move that
        // ties on M alone but has a WORSE cumulative escape count than another
        // M-tied move is not perfect play; it merely matches the mate distance.
        // `tied` below is therefore the set of every move an engine playing
        // perfectly under this definition could choose from THIS position --
        // not the broader (and looser) set of everything merely tied on M.
        vector<pair<GameState, int>> m_tied;  // stage 1: tied on M only
        if (best_val) {
            for (auto& [v, rest] : resolved_candidates) {
                if (v != *best_val) continue;
                auto& [bncum_opt, state] = rest;
                // See the inductive argument in the SearchResult struct comment: bncum_opt
                // should always be populated whenever v is; -1 is a defensive fallback only,
                // never expected to actually trigger given the current (post-DB-shortcut-
                // removal) code path.
                int bncum_val = bncum_opt ? *bncum_opt : -1;
                m_tied.emplace_back(state, bncum_val);
            }
        }

        vector<pair<GameState, int>> tied;  // stage 2: also extremal on bn_cum
        if (!m_tied.empty()) {
            int extremal_bncum = m_tied[0].second;
            for (auto& [state, bn] : m_tied) {
                if (dir == "minimize") extremal_bncum = min(extremal_bncum, bn);
                else extremal_bncum = max(extremal_bncum, bn);
            }
            for (auto& [state, bn] : m_tied) {
                if (bn == extremal_bncum) tied.emplace_back(state, bn);
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
template<typename AttackerRules>
void explore_full_attractor_dag(
    CompositionalEngine<AttackerRules>& eng, SolvedPositionDatabase& db,
    const vector<GameState>& roots, int max_depth = 16,
    int checkpoint_every = 100
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
        // LIFO (stack), not FIFO: pop from the back, same end new children get
        // pushed onto. This makes traversal depth-first -- a root's children
        // are processed immediately after it, not after every other root in
        // the frontier. With a FIFO queue, a large root file (yours has
        // 144,508 White-to-move roots) means EVERY root gets dequeued before
        // ANY child ever does, since children only ever get pushed to the
        // back -- that's exactly why an in-progress run showed zero
        // Black-to-move rows no matter how long it had been running.
        GameState st = frontier.back();
        frontier.pop_back();

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

            // best_move/cumulative_bn must come from `ns` -- the SAME candidate
            // find_best_move already selected via the existing (unmodified)
            // minimize/maximize-bn_cum tie-break -- NOT from tied_notation[0],
            // which is just whichever resolved candidate happened to be generated
            // first and carries no such guarantee. tied_notation itself is a
            // complete, correctly-computed enumeration regardless; this only
            // affects which one gets labeled the single "best_move".
            string best_move_str = eng.get_move_notation(st, *ns);
            int best_move_bncum = 0;
            for (auto& [mv_str, bncum] : tied_notation) {
                if (mv_str == best_move_str) { best_move_bncum = bncum; break; }
            }

            int val = *fval;
            SolvedPosition solution;
            solution.position_key = st.str();
            solution.turn = st.to_move;
            solution.best_move = best_move_str;
            solution.white_moves = (val + 1) / 2;
            solution.black_moves = val / 2;
            solution.M_value = solution.white_moves;
            solution.total_plies = val;
            solution.nodes_evaluated = eng.nodes_evaluated;
            solution.computation_time = search_time;
            solution.cumulative_bn = best_move_bncum;
            solution.tied_moves = tied_notation;

            db.add_position(solution, AttackerRules::full_board_symmetry);
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
            db.append_new_to_file();
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

    // One full, deduplicated rewrite at the end -- cleans up any stale duplicate
    // lines left behind by mid-run backfill re-exports (see append_new_to_file's
    // comment). Everything in between was already durable on disk via the cheap
    // incremental appends above; this isn't "the only real save", just tidying.
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

// ============================================================================
// Pawn-specific driver functions -- mirror explore_full_attractor_dag's
// structure exactly, adapted for PawnState/PawnEngine. Kept as free
// (non-templated) functions since PawnEngine has no AttackerRules parameter
// to be generic over.
// ============================================================================

vector<PawnState> load_positions_from_file_pawn(const string& filename) {
    vector<PawnState> positions;
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
        if (is_header) { is_header = false; continue; }
        if (line.empty()) continue;
        size_t comma_pos = line.find(',');
        if (comma_pos == string::npos) {
            cerr << "  Line " << line_num << " has no comma\n";
            continue;
        }
        string position_str = line.substr(comma_pos + 1);
        size_t wk_pos = position_str.find("WK:") + 3;
        size_t wq_pos = position_str.find("WQ:") + 3;
        size_t bk_pos = position_str.find("BK:") + 3;
        try {
            Position wk = parse_position_from_string(position_str.substr(wk_pos, 2));
            Position wp = parse_position_from_string(position_str.substr(wq_pos, 2));
            Position bk = parse_position_from_string(position_str.substr(bk_pos, 2));
            // Root/seed positions for a pawn endgame always start as an actual
            // pawn -- promotion is something that happens DURING the game,
            // never something a starting position begins as.
            positions.push_back(PawnState(wk, wp, PieceKind::PAWN, bk, 'W'));
        } catch (const exception& e) {
            cerr << "  Line " << line_num << " parse error: " << e.what() << "\n";
        }
    }
    file.close();
    cout << "  Loaded " << positions.size() << " positions from " << filename << "\n";
    return positions;
}

// 3 bits per coordinate x 6 coordinates + 1 bit turn + 2 bits piece kind = 21
// bits, comfortably fits uint64_t. Kind needs its own bits here for the same
// reason as make_cache_key above: once a pawn can promote, (square, turn)
// alone no longer uniquely identifies a node -- a pawn and an already-
// promoted piece can occupy the identical square as genuinely different
// states.
uint64_t pack_pawn_state_key(const PawnState& st) {
    uint64_t key = 0;
    key |= ((uint64_t)st.wk.file << 20);
    key |= ((uint64_t)st.wk.rank << 17);
    key |= ((uint64_t)st.wp.file << 14);
    key |= ((uint64_t)st.wp.rank << 11);
    key |= ((uint64_t)st.bk.file << 8);
    key |= ((uint64_t)st.bk.rank << 5);
    key |= (st.to_move == 'W' ? 1ull : 0ull) << 4;
    key |= (uint64_t)static_cast<uint8_t>(st.wp_kind) << 2;
    return key;
}

void explore_full_attractor_dag_pawn(
    PawnEngine& eng, SolvedPositionDatabase& db,
    const vector<PawnState>& roots, int max_depth = 30,
    int checkpoint_every = 100
) {
    cout << "\n" << string(80, '=') << "\n";
    cout << "FULL ATTRACTOR-DAG EXPLORATION (PAWN)\n";
    cout << string(80, '=') << "\n\n";

    unordered_set<uint64_t> visited_this_run;
    deque<PawnState> frontier;

    for (auto& r : roots) {
        if (!eng.is_checkmate(r)) frontier.push_back(r);
    }
    cout << "Seeded frontier with " << frontier.size() << " root positions\n\n";

    long long processed = 0, newly_solved = 0, cache_reused = 0, failed = 0;
    auto run_start = chrono::high_resolution_clock::now();

    while (!frontier.empty()) {
        PawnState st = frontier.back();
        frontier.pop_back();

        uint64_t key = pack_pawn_state_key(st);
        if (visited_this_run.count(key)) continue;
        visited_this_run.insert(key);

        if (eng.is_checkmate(st)) continue;

        vector<pair<string,int>> tied_notation;

        // AttackerKind must be part of the DB lookup key too -- reuse the
        // existing (position_string, turn) lookup, but the resulting cached
        // row's own AttackerKind column must match st.wp_kind, or this is a
        // DIFFERENT state that happens to share a position string with an
        // already-solved one (a pawn and an already-promoted piece on the
        // same square). Mismatched kind is treated as "not actually cached"
        // and solved fresh, exactly like a genuinely new position.
        auto cached = db.get_solution(st.str(), st.to_move);
        bool cache_valid = cached && kind_from_letter(cached->attacker_kind) == st.wp_kind;

        if (cache_valid && !cached->tied_moves.empty()) {
            tied_notation = cached->tied_moves;
            cache_reused++;
            processed++;
        } else {
            auto search_start = chrono::high_resolution_clock::now();
            auto [ns, fval, tied] = eng.find_best_move(st, 2 * max_depth);
            auto search_end = chrono::high_resolution_clock::now();
            double search_time = chrono::duration<double>(search_end - search_start).count();

            if (!fval) {
                failed++;
                cout << "  [FAILED] " << st.str() << " kind=" << kind_letter(st.wp_kind)
                     << " (to_move=" << st.to_move << ") -- unresolved within depth budget "
                     << (2 * max_depth) << "\n";
                continue;
            }

            for (auto& [tstate, tbncum] : tied) {
                tied_notation.emplace_back(eng.get_move_notation(st, tstate), tbncum);
            }

            string best_move_str = eng.get_move_notation(st, *ns);
            int best_move_bncum = 0;
            for (auto& [mv_str, bncum] : tied_notation) {
                if (mv_str == best_move_str) { best_move_bncum = bncum; break; }
            }

            int val = *fval;
            SolvedPosition solution;
            solution.position_key = st.str();
            solution.turn = st.to_move;
            solution.best_move = best_move_str;
            solution.white_moves = (val + 1) / 2;
            solution.black_moves = val / 2;
            solution.M_value = solution.white_moves;
            solution.total_plies = val;
            solution.nodes_evaluated = eng.nodes_evaluated;
            solution.computation_time = search_time;
            solution.cumulative_bn = best_move_bncum;
            solution.tied_moves = tied_notation;
            solution.attacker_kind = kind_letter(st.wp_kind);

            // Pawn positions never use full board symmetry (see add_position's
            // comment) -- only identity + horizontal-mirror are valid.
            db.add_position(solution, false);
            newly_solved++;
            processed++;
        }

        for (auto& [mv_str, bncum] : tied_notation) {
            (void)bncum;
            PawnState child = eng.apply_move_notation(st, mv_str);
            uint64_t ckey = pack_pawn_state_key(child);
            if (!visited_this_run.count(ckey)) {
                frontier.push_back(child);
            }
        }

        if (processed % checkpoint_every == 0) {
            db.append_new_to_file();
            double elapsed = chrono::duration<double>(
                chrono::high_resolution_clock::now() - run_start).count();
            cout << "[CHECKPOINT] processed=" << processed
                 << " newly_solved=" << newly_solved
                 << " cache_reused=" << cache_reused
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
    cout << "FULL ATTRACTOR-DAG EXPLORATION (PAWN) COMPLETE\n";
    cout << string(80, '=') << "\n";
    cout << "Distinct nodes visited: " << visited_this_run.size() << "\n";
    cout << "Newly solved this run:  " << newly_solved << "\n";
    cout << "Reused from DB cache:   " << cache_reused << "\n";
    cout << "Failed to resolve:      " << failed << "\n";
    cout << "Total time:             " << fixed << setprecision(1) << total_time << "s\n\n";
}

template<typename AttackerRules>
void batch_solve_all_kqvk_positions(CompositionalEngine<AttackerRules>& eng, SolvedPositionDatabase& db, int max_depth = 16) {
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

            db.add_position(solution, AttackerRules::full_board_symmetry);
            solved_count++;
            recorded_count++;

            char piece = mvs[i][0];
            string dest = mvs[i].substr(1);
            Position dest_pos = Position::from_str(dest);
            if (piece == 'K') { curr.wk = dest_pos; curr.to_move = 'B'; }
            else if (piece == AttackerRules::letter) { curr.wq = dest_pos; curr.to_move = 'B'; }
            else if (piece == 'k') { curr.bk = dest_pos; curr.to_move = 'W'; }
        }

        cout << "[" << idx << "/" << positions.size() << "] [SOLVED] " 
             << pos.str() << " plies=" << tp << " recorded=" << recorded_count
             << " (" << fixed << setprecision(3) << solve_time << "s)\n";
        // end of root cause issue

        // Checkpoint every 5 positions. append_new_to_file() only writes what's new
        // since the last flush (O(delta)), not the previous export_to_file()'s full
        // rewrite of the whole file (O(total size)) -- watch export_time below stay
        // small and roughly flat as solved_count grows, instead of climbing with it.
        if (solved_count % 5 == 0 && solved_count > 0) {
            auto export_start = chrono::high_resolution_clock::now();
            db.append_new_to_file();
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

// ============================================================================
// Which endgame this binary is built for. Default (no flag) is KQvK. Compile
// with -DPIECE_ROOK to build KRvK instead -- same source file, same search
// and DAG-exploration logic, just a different AttackerRules instantiation.
//
// -DPIECE_PAWN builds KPvK, and takes a genuinely different path through
// main() below: PawnState/PawnEngine are a separate, parallel implementation
// (see the large comment where they're defined), not another AttackerRules
// instantiation of Engine/GameState, so they need their own driver logic
// rather than slotting into the existing Engine alias.
// ============================================================================
#ifndef PIECE_PAWN
#ifndef PIECE_GENERAL
#ifdef PIECE_ROOK
using Engine = CompositionalEngine<RookRules>;
#else
using Engine = CompositionalEngine<QueenRules>;
#endif
#endif
#endif

int main(int argc, char* argv[]) {
    bool debug_en = false;
    bool batch_mode = false;
    bool full_dag_mode = false;
    bool fresh_start = false;
#ifdef PIECE_GENERAL
    string positions_file = "general_positions.txt";
    string db_file = "general_perfect_play.db";
#elif defined(PIECE_PAWN)
    string positions_file = "kpvk_positions_by_dtz.txt";
    string db_file = "kpvk_perfect_play.db";
#elif defined(PIECE_ROOK)
    string positions_file = "krvk_positions_by_dtz.txt";
    string db_file = "krvk_perfect_play.db";
#else
    string positions_file = "kqvk_positions_by_dtz.txt";
    string db_file = "kqvk_perfect_play.db";
#endif
    vector<string> unrecognized;

    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        if (arg == "--debug" || arg == "-d") { debug_en = true; }
        else if (arg == "--batch" || arg == "-b") { batch_mode = true; }
        else if (arg == "--full-dag" || arg == "-g") { full_dag_mode = true; }
        else if (arg == "--positions" && i + 1 < argc) { positions_file = argv[++i]; }
        else if (arg == "--db" && i + 1 < argc) { db_file = argv[++i]; }
        else if (arg == "--fresh") { fresh_start = true; }
        else { unrecognized.push_back(arg); }
    }

    // A silently-ignored flag is exactly the bug class that caused --db to be a
    // no-op in the previous version of this file: the arg-parsing loop just skipped
    // anything it didn't recognize instead of complaining. Refusing to start on an
    // unrecognized argument (rather than plausibly running against the wrong file)
    // is the safer failure mode here.
    if (!unrecognized.empty()) {
        cerr << "ERROR: unrecognized argument(s):";
        for (auto& u : unrecognized) cerr << " " << u;
        cerr << "\nKnown flags: --debug/-d, --batch/-b, --full-dag/-g, "
                "--positions <file>, --db <file>, --fresh\n";
        return 1;
    }

    // --fresh exists because of a real, demonstrated failure mode: if db_file
    // already contains an entry for a position, is_solved() treats it as
    // already-decided and reuses it directly (see explore_full_attractor_dag's
    // cache_reused path) -- it is NEVER recomputed, even after a search-logic
    // fix that would have produced a different, correct answer for it. That
    // makes "delete the old database before rerunning with a fixed binary" a
    // manual step easy to forget and expensive to forget silently: the run
    // completes normally, checkpoints normally, and just quietly carries every
    // stale answer forward unchanged. --fresh removes the file first so a
    // clean start is guaranteed by the tool rather than by remembering to `rm`/
    // `del` it yourself.
    if (fresh_start) {
        if (remove(db_file.c_str()) == 0) {
            cout << "--fresh: removed existing " << db_file << " before starting\n";
        } else {
            cout << "--fresh: no existing " << db_file << " to remove (starting fresh anyway)\n";
        }
    }
    
    double root_t = 0;
    
    cout << "\n" << string(80, '=') << "\n";
    cout << "COMPOSITIONAL TOPOLOGICAL SEARCH\n";
    cout << "Measurement with Black Node Count Accumulation\n";
    cout << string(80, '=') << "\n";
    cout << "Database file: " << db_file << "\n";
    if (full_dag_mode) cout << "Mode: full attractor-DAG exploration\n";
    else cout << "Mode: single-path batch solve\n";

    auto root_t_start = chrono::high_resolution_clock::now();

#ifdef PIECE_GENERAL
    // Uses its own CSV format (export_general_classification), not
    // SolvedPositionDatabase -- that format's add_position hardcodes
    // exactly three roles (WK:/WQ:/BK:) and was never designed for a
    // variable-length piece list. --batch has no meaning here either: the
    // classifier always exhaustively solves everything reachable from the
    // given roots in one pass, there is no separate single-path mode.
    if (!full_dag_mode) {
        cerr << "ERROR: --batch is not meaningful for PIECE_GENERAL; use --full-dag.\n";
        return 1;
    }
    vector<GeneralState> roots = load_general_positions_from_file(positions_file);
    if (roots.empty()) {
        cerr << "ERROR: No root positions loaded from " << positions_file << "\n";
        return 1;
    }
    run_general_sweep(roots, db_file);
#else
    SolvedPositionDatabase db(db_file);

#ifdef PIECE_PAWN
    PawnEngine eng;
    if (!full_dag_mode) {
        // batch_solve_all_kqvk_positions has no pawn-aware equivalent yet --
        // refusing to silently do the wrong thing is safer than pretending
        // to support a mode that was never built for this piece.
        cerr << "ERROR: --batch is not implemented for KPvK yet; use --full-dag.\n";
        return 1;
    }
    vector<PawnState> roots = load_positions_from_file_pawn(positions_file);
    if (roots.empty()) {
        cerr << "ERROR: No root positions loaded from " << positions_file << "\n";
        return 1;
    }
    explore_full_attractor_dag_pawn(eng, db, roots, 30);
#else
    Engine eng;
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
#endif
#endif
    return 0;
}