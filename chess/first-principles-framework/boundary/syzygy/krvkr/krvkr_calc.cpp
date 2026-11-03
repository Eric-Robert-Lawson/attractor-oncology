// KRvKR primitive calculator. No tablebase is read by the rules. The harness (main) reads the table-derived
// arrays only to compare the rules against them.
//
// Rules (all in terms of the catalog's primitives, written for this material):
//  Net terms (P-001): C attack, F escape, X capture, B block, O outcome.
//   - Royal target:  mate in 1 = the mover has a move after which the defender is checkmated (C,F,X,B hold).
//   - Non-royal target (the defender's rook): the capture is available (C) and its outcome O is favourable:
//       the capturing unit cannot be netted by the defender (no recapture) and the defender is not stalemated.
//     A capture whose outcome is not favourable is a drawing exit (the exit condition); a favourable one is a win.
//  Composition (P-003 across the irreversible boundary, P-001 nested): the capture is the boundary. Before it the goal is the net on the
//   rook; after a favourable capture the goal is the mating net on the king, which is decided (the net cannot be refused).
//   W(h): the mover wins within h plies  = a favourable capture, or mate in 1, or a move from the pool after which the
//         defender loses within h-1 plies.
//   L(h): the mover loses within h plies = checkmated (h=0), or the mover has no capture of the rook and every move leads to a position
//         where the opponent wins within h-1 plies (every move in the pool is refuted by a net).
//  Status: win if W(H); loss if L(H-1); otherwise drawn (the complement).
//  Pool (the scope of the winner's candidate moves) is a parameter: 0 full, 1 checks, 2 checks + moves that attack the rook,
//   3 and 4 as set below by the counterexamples found.
#include <bits/stdc++.h>
using namespace std;
typedef unsigned char u8;
static const int NPOS = 64 * 64 * 64 * 64 * 2;
static int32_t* canon;
static inline int fl(int s) { return s & 7; }
static inline int rk(int s) { return s >> 3; }
static inline bool adj(int a, int b) { return a != b && max(abs(fl(a) - fl(b)), abs(rk(a) - rk(b))) == 1; }
static inline bool rookAtt(int from, int to, const int* bl, int n) {
    if (from == to) return false;
    int ff = fl(from), fr = rk(from), tf = fl(to), tr = rk(to);
    if (ff != tf && fr != tr) return false;
    int df = (tf > ff) - (tf < ff), dr = (tr > fr) - (tr < fr);
    int f = ff + df, r = fr + dr;
    while (!(f == tf && r == tr)) {
        int s = r * 8 + f;
        for (int i = 0; i < n; i++) if (bl[i] == s) return false;
        f += df; r += dr;
    }
    return true;
}
struct Pos { int k[2], r[2]; int t; };
struct Mv { int piece; int to; bool cap; };
static inline int rawIdx(const Pos& p) { return ((((p.k[0] * 64 + p.r[0]) * 64 + p.k[1]) * 64 + p.r[1]) * 2 + p.t); }
static Pos fromIdx(int i) {
    Pos p; p.t = i & 1; i >>= 1; p.r[1] = i & 63; i >>= 6; p.k[1] = i & 63; i >>= 6; p.r[0] = i & 63; i >>= 6; p.k[0] = i & 63; return p;
}
static const int DF[8] = {1, 1, 1, 0, 0, -1, -1, -1}, DR[8] = {1, 0, -1, 1, -1, 1, 0, -1};
static bool inCheck(const Pos& p) {
    int m = p.t, o = 1 - m; int bl[2] = {p.k[o], p.r[m]};
    return rookAtt(p.r[o], p.k[m], bl, 2);
}
static long GENCALLS = 0;
static int genLegal(const Pos& p, Mv* out) {
    GENCALLS++;
    int m = p.t, o = 1 - m; int K = p.k[m], R = p.r[m], OK = p.k[o], OR = p.r[o]; int n = 0;
    for (int d = 0; d < 8; d++) {
        int f = fl(K) + DF[d], r = rk(K) + DR[d];
        if (f < 0 || f > 7 || r < 0 || r > 7) continue;
        int to = r * 8 + f;
        if (to == R || to == OK || adj(to, OK)) continue;
        bool cap = (to == OR);
        if (!cap) { int bl[2] = {OK, R}; if (rookAtt(OR, to, bl, 2)) continue; }
        out[n++] = {0, to, cap};
    }
    static const int RD[4][2] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
    for (int d = 0; d < 4; d++) {
        int f = fl(R), r = rk(R);
        while (true) {
            f += RD[d][0]; r += RD[d][1];
            if (f < 0 || f > 7 || r < 0 || r > 7) break;
            int to = r * 8 + f;
            if (to == K || to == OK) break;
            if (to == OR) { out[n++] = {1, to, true}; break; }
            int bl[2] = {OK, to};
            if (!rookAtt(OR, K, bl, 2)) out[n++] = {1, to, false};
        }
    }
    return n;
}
static inline Pos applyMv(const Pos& p, const Mv& mv) {
    Pos c = p; int m = p.t;
    if (mv.piece == 0) c.k[m] = mv.to; else c.r[m] = mv.to;
    c.t = 1 - m; return c;
}
// O term of the net on the defender's rook: after the capture, can the defender net the capturing side's rook, and is the defender stalemated?
static bool favourableCapture(const Pos& p, const Mv& mv) {
    int m = p.t, o = 1 - m; int K = p.k[m], R = p.r[m], OK = p.k[o];
    if (mv.piece == 0) K = mv.to; else R = mv.to;
    int cnt = 0; bool rookCap = false;
    for (int d = 0; d < 8; d++) {
        int f = fl(OK) + DF[d], r = rk(OK) + DR[d];
        if (f < 0 || f > 7 || r < 0 || r > 7) continue;
        int to = r * 8 + f;
        if (to == K || adj(to, K)) continue;
        if (to == R) { rookCap = true; cnt++; continue; }
        int bl[1] = {K};
        if (rookAtt(R, to, bl, 1)) continue;
        cnt++;
    }
    int bl[1] = {K};
    bool chk = rookAtt(R, OK, bl, 1);
    return !rookCap && !(cnt == 0 && !chk);
}
static int MODE = 0;
// the winner's candidate pool
static bool threat(const Pos& p, const Mv& mv);
static bool inPool(const Pos& p, const Mv& mv) {
    if (MODE == 0) return true;
    int m = p.t, o = 1 - m;
    int K = p.k[m], R = p.r[m], OK = p.k[o], OR = p.r[o];
    if (mv.piece == 0) K = mv.to; else R = mv.to;
    int bl[2] = {K, OR};
    bool chk = rookAtt(R, OK, bl, 2);
    if (chk) return true;
    if (MODE == 1) return false;
    if (MODE == 3 || MODE == 5) return threat(p, mv);
    if (MODE == 4) { int bl1[1] = {K}; if (rookAtt(R, OR, bl1, 1) || adj(K, OR)) return true; return threat(p, mv); }
    // attacks the rook: the rook sees it, with the defender's king allowed as the only blocker (x-ray), or the king stands next to it
    int bl1[1] = {K};
    bool att = rookAtt(R, OR, bl1, 1) || adj(K, OR);
    if (att) return true;
    if (MODE == 2) return false;
    return true;
}
static bool Wf(int h, const Pos& p);
// threat: if the defender passed, the mover would have a favourable capture or a mate in one (a net completable next move)
static int CURH = 1;
static bool threat(const Pos& p, const Mv& mv) {
    Pos c = applyMv(p, mv); c.t = p.t;
    return Wf(CURH, c);
}
static u8 *Wm[6], *Lm[6];
static bool Lf(int h, const Pos& p);
static bool Wf(int h, const Pos& p) {   // h odd
    int key = canon[rawIdx(p)]; u8& memo = Wm[(h - 1) / 2][key];
    if (memo) return memo == 1;
    bool res = false;
    Mv mv[32]; int n = genLegal(p, mv);
    for (int i = 0; i < n && !res; i++) if (mv[i].cap && favourableCapture(p, mv[i])) res = true;
    if (!res && h >= 3) res = Wf(h - 2, p);
    if (!res) for (int i = 0; i < n && !res; i++) {
        if (mv[i].cap) continue;
        if (h > 1) { CURH = (MODE == 5 ? h - 2 : 1); if (!inPool(p, mv[i])) continue; }   // mate in 1 is read directly (a mating move is a check, so it is in every pool)
        Pos c = applyMv(p, mv[i]);
        if (Lf(h - 1, c)) res = true;
    }
    memo = res ? 1 : 2; return res;
}
static bool Lf(int h, const Pos& p) {   // h even
    int key = canon[rawIdx(p)]; u8& memo = Lm[h / 2][key];
    if (memo) return memo == 1;
    bool res;
    Mv mv[32]; int n = genLegal(p, mv);
    if (h == 0) res = (n == 0 && inCheck(p));
    else if (Lf(h - 2, p)) res = true;
    else if (n == 0) res = false;
    else {
        res = true;
        for (int i = 0; i < n && res; i++) {
            if (mv[i].cap) { res = false; break; }
            Pos c = applyMv(p, mv[i]);
            if (!Wf(h - 1, c)) res = false;
        }
    }
    memo = res ? 1 : 2; return res;
}
static const int HW = 7, HL = 6;
static inline int value(const Pos& p) { if (Wf(HW, p)) return 1; if (Lf(HL, p)) return -1; return 0; }

// number of legal replies of the defender (a bare king) after the mover's capture of the rook
static int defenderReplies(const Pos& p, const Mv& mv) {
    int m = p.t, o = 1 - m; int K = p.k[m], R = p.r[m], OK = p.k[o];
    if (mv.piece == 0) K = mv.to; else R = mv.to;
    int cnt = 0;
    for (int d = 0; d < 8; d++) {
        int f = fl(OK) + DF[d], r = rk(OK) + DR[d];
        if (f < 0 || f > 7 || r < 0 || r > 7) continue;
        int to = r * 8 + f;
        if (to == K || adj(to, K)) continue;
        if (to == R) { cnt++; continue; }
        int bl[1] = {K};
        if (rookAtt(R, to, bl, 1)) continue;
        cnt++;
    }
    return cnt;
}
static int originOf(const Pos& p, const Mv& mv) { return mv.piece == 0 ? p.k[p.t] : p.r[p.t]; }
static bool findMove(const Pos& p, int piece, int to, Mv& out) {
    Mv mv[32]; int k = genLegal(p, mv);
    for (int i = 0; i < k; i++) if (!mv[i].cap && mv[i].piece == piece && mv[i].to == to) { out = mv[i]; return true; }
    return false;
}
// a back-and-forth cycle of four plies through p: both sides move a piece and move it back, every position drawn
static bool shuffle4(const Pos& p0) {
    Mv a[32]; int na = genLegal(p0, a);
    for (int i = 0; i < na; i++) {
        if (a[i].cap) continue;
        Pos c1 = applyMv(p0, a[i]); if (value(c1) != 0) continue;
        int o1 = originOf(p0, a[i]);
        Mv b[32]; int nb = genLegal(c1, b);
        for (int j = 0; j < nb; j++) {
            if (b[j].cap) continue;
            Pos c2 = applyMv(c1, b[j]); if (value(c2) != 0) continue;
            int o2 = originOf(c1, b[j]);
            Mv r1; if (!findMove(c2, a[i].piece, o1, r1)) continue;
            Pos c3 = applyMv(c2, r1); if (value(c3) != 0) continue;
            Mv r2; if (!findMove(c3, b[j].piece, o2, r2)) continue;
            Pos c4 = applyMv(c3, r2);
            if (rawIdx(c4) == rawIdx(p0)) return true;
        }
    }
    return false;
}
static void runCycles(const vector<int32_t>& reps, const vector<u8>& valid) {
    int M = reps.size();
    vector<int> repid(NPOS, -1);
    for (int n = 0; n < M; n++) if (valid[n]) repid[reps[n]] = n;
    vector<int8_t> v(M, 9);
    for (int n = 0; n < M; n++) if (valid[n]) v[n] = value(fromIdx(reps[n]));
    vector<int> eoff(M + 1, 0), eto; eto.reserve(15000000);
    vector<char> exitDraw(M, 0), exitDev(M, 0), dev(M, 0), stale(M, 0);
    long nd = 0;
    for (int n = 0; n < M; n++) {
        eoff[n] = eto.size();
        if (!valid[n] || v[n] != 0) continue;
        nd++;
        Pos p = fromIdx(reps[n]); Mv mv[32]; int k = genLegal(p, mv); int nl = 0;
        if (k == 0 && !inCheck(p)) stale[n] = 1;
        for (int i = 0; i < k; i++) {
            if (mv[i].cap) { if (!favourableCapture(p, mv[i])) { exitDraw[n] = 1; if (defenderReplies(p, mv[i]) >= 2) exitDev[n] = 1; } continue; }
            Pos c = applyMv(p, mv[i]); int cv = value(c);
            if (cv == 1) nl++; else if (cv == 0) { int cid = repid[canon[rawIdx(c)]]; if (cid < 0) { printf("bad child\n"); exit(3); } eto.push_back(cid); }
        }
        dev[n] = nl > 0;
    }
    eoff[M] = eto.size();
    printf("drawn positions %ld, drawing edges within the material %zu\n", nd, eto.size());
    // strongly connected components (iterative Tarjan)
    vector<int> idx(M, -1), low(M, 0), comp(M, -1), stk, callst, itpos(M, 0); vector<char> onst(M, 0);
    int counter = 0, ncomp = 0;
    vector<int> csize;
    for (int s0 = 0; s0 < M; s0++) {
        if (!valid[s0] || v[s0] != 0 || idx[s0] != -1) continue;
        callst.push_back(s0); idx[s0] = low[s0] = counter++; stk.push_back(s0); onst[s0] = 1; itpos[s0] = eoff[s0];
        while (!callst.empty()) {
            int u = callst.back();
            if (itpos[u] < eoff[u + 1]) {
                int w = eto[itpos[u]++];
                if (idx[w] == -1) { idx[w] = low[w] = counter++; stk.push_back(w); onst[w] = 1; itpos[w] = eoff[w]; callst.push_back(w); }
                else if (onst[w]) low[u] = min(low[u], idx[w]);
            } else {
                if (low[u] == idx[u]) { int sz = 0; while (true) { int w = stk.back(); stk.pop_back(); onst[w] = 0; comp[w] = ncomp; sz++; if (w == u) break; } csize.push_back(sz); ncomp++; }
                callst.pop_back();
                if (!callst.empty()) { int par = callst.back(); low[par] = min(low[par], low[u]); }
            }
        }
    }
    vector<char> selfl(M, 0), onCycle(M, 0);
    for (int n = 0; n < M; n++) for (int e = eoff[n]; e < eoff[n + 1]; e++) if (eto[e] == n) selfl[n] = 1;
    long ncyc = 0, big = 0; int bigc = -1;
    for (int c = 0; c < ncomp; c++) if (csize[c] > big) { big = csize[c]; bigc = c; }
    for (int n = 0; n < M; n++) if (valid[n] && v[n] == 0) { onCycle[n] = (csize[comp[n]] > 1) || selfl[n]; ncyc += onCycle[n]; }
    long nontriv = 0; for (int c = 0; c < ncomp; c++) if (csize[c] > 1) nontriv++;
    printf("strongly connected components: %d, with more than one position: %ld, largest %ld positions\n", ncomp, nontriv, big);
    printf("drawn positions on a repetition cycle: %ld, on none: %ld\n", ncyc, nd - ncyc);
    // reachability of a cycle and of a terminal draw, and of a losing deviation for each side
    vector<char> rc(M, 0), rt(M, 0), xe(M, 0), xo(M, 0);
    for (int n = 0; n < M; n++) if (valid[n] && v[n] == 0) { rc[n] = onCycle[n]; rt[n] = exitDraw[n] || stale[n]; xe[n] = dev[n]; xo[n] = exitDev[n]; }
    bool ch = true; int it = 0;
    while (ch) { ch = false; it++;
        for (int n = 0; n < M; n++) { if (!valid[n] || v[n] != 0) continue;
            for (int e = eoff[n]; e < eoff[n + 1]; e++) { int c = eto[e];
                if (rc[c] && !rc[n]) { rc[n] = 1; ch = true; }
                if (rt[c] && !rt[n]) { rt[n] = 1; ch = true; }
                if (xo[c] && !xe[n]) { xe[n] = 1; ch = true; }
                if (xe[c] && !xo[n]) { xo[n] = 1; ch = true; } } } }
    printf("propagation rounds %d\n", it);
    long cls[5] = {}, two[5] = {}, mixedTwo = 0;
    const char* cn[5] = {"stalemate now", "terminal only (forced)", "cycle only", "mixed", "neither"};
    for (int n = 0; n < M; n++) { if (!valid[n] || v[n] != 0) continue;
        int c = stale[n] ? 0 : (rc[n] && rt[n]) ? 3 : rc[n] ? 2 : rt[n] ? 1 : 4;
        cls[c]++; if (xe[n] && xo[n]) two[c]++; }
    printf("realization of the draw, and how many of each have a losing deviation for each side somewhere in the region of drawing moves:\n");
    for (int c = 0; c < 5; c++) printf("  %-24s %9ld   two-sided %9ld\n", cn[c], cls[c], two[c]);
    long tw = 0, strict = 0; for (int c = 0; c < 5; c++) tw += two[c];
    printf("drawn positions with a losing deviation for each side in the region: %ld (%.2f%% of the drawn)\n", tw, 100.0 * tw / nd);
    // a four-ply back-and-forth cycle through the position itself
    vector<char> sh(M, 0); long nsh = 0;
    for (int n = 0; n < M; n++) { if (!valid[n] || v[n] != 0) continue; if (shuffle4(fromIdx(reps[n]))) { sh[n] = 1; nsh++; } }
    printf("drawn positions with a four-ply back-and-forth cycle through them: %ld (%.2f%%)\n", nsh, 100.0 * nsh / nd);
    // distance in plies, along drawing moves, to the nearest position that has such a cycle
    vector<int> roff(M + 2, 0), rfrom(eto.size());
    for (int n = 0; n < M; n++) for (int e = eoff[n]; e < eoff[n + 1]; e++) roff[eto[e] + 1]++;
    for (int i = 0; i < M; i++) roff[i + 1] += roff[i];
    { vector<int> pos(roff.begin(), roff.begin() + M);
      for (int n = 0; n < M; n++) for (int e = eoff[n]; e < eoff[n + 1]; e++) rfrom[pos[eto[e]]++] = n; }
    vector<int> dist(M, -1); vector<int> q;
    for (int n = 0; n < M; n++) if (sh[n]) { dist[n] = 0; q.push_back(n); }
    for (size_t h = 0; h < q.size(); h++) { int u = q[h]; for (int e = roff[u]; e < roff[u + 1]; e++) { int w = rfrom[e]; if (dist[w] < 0) { dist[w] = dist[u] + 1; q.push_back(w); } } }
    map<int, long> dd; long none = 0, noneOnCycle = 0;
    for (int n = 0; n < M; n++) { if (!valid[n] || v[n] != 0) continue; if (dist[n] < 0) { none++; if (onCycle[n]) noneOnCycle++; } else dd[dist[n]]++; }
    printf("plies along drawing moves to the nearest position with a four-ply cycle:\n");
    long cum = 0; for (auto& kv : dd) { cum += kv.second; printf("  %2d plies: %9ld   cumulative %.4f\n", kv.first, kv.second, (double)cum / nd); }
    printf("  none reachable: %ld (of which on a longer cycle: %ld)\n", none, noneOnCycle);
    // the drawn positions that reach no cycle: what they are
    long noCyc = 0, noCycTwo = 0; for (int n = 0; n < M; n++) if (valid[n] && v[n] == 0 && !rc[n]) { noCyc++; if (xe[n] && xo[n]) noCycTwo++; }
    printf("drawn positions from which no cycle is reachable: %ld, two-sided %ld\n", noCyc, noCycTwo);
}
static void describe(const Pos& p);
int main(int argc, char** argv) {
    bool evalMode = argc > 1 && string(argv[1]) == "eval";
    bool cold = argc > 1 && string(argv[1]) == "cold";
    bool tempo = argc > 1 && string(argv[1]) == "tempo";
    bool cyclesMode = argc > 1 && string(argv[1]) == "cycles";
    MODE = (evalMode || tempo || cyclesMode) ? 5 : cold ? atoi(argv[2]) : (argc > 1 ? atoi(argv[1]) : 0);
    int limit = argc > 2 ? atoi(argv[2]) : 0;
    canon = (int32_t*)malloc(sizeof(int32_t) * NPOS);
    for (int i = 0; i < NPOS; i++) {
        Pos p = fromIdx(i); int best = i;
        for (int s = 0; s < 8; s++) {
            int a[4] = {p.k[0], p.r[0], p.k[1], p.r[1]}, b[4];
            for (int j = 0; j < 4; j++) {
                int f = fl(a[j]), r = rk(a[j]), nf, nr;
                switch (s) { case 0: nf = f; nr = r; break; case 1: nf = 7 - f; nr = r; break; case 2: nf = f; nr = 7 - r; break; case 3: nf = 7 - f; nr = 7 - r; break;
                    case 4: nf = r; nr = f; break; case 5: nf = 7 - r; nr = f; break; case 6: nf = r; nr = 7 - f; break; default: nf = 7 - r; nr = 7 - f; }
                b[j] = nr * 8 + nf;
            }
            int i1 = (((b[0] * 64 + b[1]) * 64 + b[2]) * 64 + b[3]) * 2 + p.t;
            int i2 = (((b[2] * 64 + b[3]) * 64 + b[0]) * 64 + b[1]) * 2 + (1 - p.t);
            best = min(best, min(i1, i2));
        }
        canon[i] = best;
    }
    for (int i = 0; i < 6; i++) { Wm[i] = (u8*)calloc(NPOS, 1); Lm[i] = (u8*)calloc(NPOS, 1); }
    if (tempo) {
        // the same board with each side to move, read from White's point of view
        long tab[3][3] = {}; long boards = 0;
        for (int wk = 0; wk < 64; wk++) for (int wr = 0; wr < 64; wr++) { if (wr == wk) continue;
         for (int bk = 0; bk < 64; bk++) { if (bk == wk || bk == wr || adj(bk, wk)) continue;
          for (int br = 0; br < 64; br++) { if (br == wk || br == wr || br == bk) continue;
            Pos p; p.k[0] = wk; p.r[0] = wr; p.k[1] = bk; p.r[1] = br;
            int b1[2] = {wk, br}, b2[2] = {bk, wr};
            if (rookAtt(wr, bk, b1, 2)) continue;        // black king attacked by the white rook
            if (rookAtt(br, wk, b2, 2)) continue;        // white king attacked by the black rook
            p.t = 0; int vw = value(p); p.t = 1; int vb = -value(p);   // White's point of view
            tab[vw + 1][vb + 1]++; boards++;
          } } }
        const char* nm3[3] = {"White loses", "drawn", "White wins"};
        printf("boards where neither king is in check: %ld (not folded)\nrows: White to move; columns: Black to move; all read as White's result\n", boards);
        for (int a = 0; a < 3; a++) { printf("%12s:", nm3[a]); for (int c = 0; c < 3; c++) printf(" %10ld", tab[a][c]); printf("\n"); }
        long diff = boards - tab[0][0] - tab[1][1] - tab[2][2];
        printf("boards whose result depends on the side to move: %ld (%.2f%%)\n", diff, 100.0 * diff / boards);
        return 0;
    }
    if (evalMode) {
        Pos p; p.k[0] = atoi(argv[2]); p.r[0] = atoi(argv[3]); p.k[1] = atoi(argv[4]); p.r[1] = atoi(argv[5]); p.t = atoi(argv[6]);
        describe(p); return 0;
    }
    int M = 1906548;
    vector<int32_t> reps(M); vector<u8> valid(M); vector<int8_t> val(M); vector<int16_t> nlose(M), ndi(M), ndc(M), nm(M);
    auto rd = [&](const char* f, void* d, size_t sz) { FILE* fp = fopen(f, "rb"); if (!fp) { perror(f); exit(1); } if (fread(d, 1, sz, fp) != sz) exit(2); fclose(fp); };
    rd("reps.bin", reps.data(), M * 4); rd("valid.bin", valid.data(), M); rd("val.bin", val.data(), M);
    rd("nlose.bin", nlose.data(), M * 2); rd("ndi.bin", ndi.data(), M * 2); rd("ndc.bin", ndc.data(), M * 2); rd("nm.bin", nm.data(), M * 2);
    if (cold) {
        int N = atoi(argv[3]); mt19937 rng(7); vector<long> calls[3];
        for (int c = 0; c < N; c++) {
            int n; do { n = rng() % M; } while (!valid[n]);
            for (int i = 0; i < 6; i++) { memset(Wm[i], 0, NPOS); memset(Lm[i], 0, NPOS); }
            GENCALLS = 0; Pos p = fromIdx(reps[n]); int v = value(p);
            calls[v + 1].push_back(GENCALLS);
        }
        const char* nm3[3] = {"loss", "draw", "win"};
        for (int a = 0; a < 3; a++) { auto& v = calls[a]; if (v.empty()) continue; sort(v.begin(), v.end()); double mean = 0; for (long x : v) mean += x; mean /= v.size();
            printf("cold start, %s (%zu positions): move generations per position, mean %.0f, median %ld, max %ld\n", nm3[a], v.size(), mean, v[v.size() / 2], v.back()); }
        return 0;
    }
    if (cyclesMode) { runCycles(vector<int32_t>(reps), valid); return 0; }
    long conf[3][3] = {}; long total = 0; long moveBad = 0, nmBad = 0; long posBadMoves = 0;
    vector<int> bad;
    for (int n = 0; n < M; n++) {
        if (!valid[n]) continue;
        if (limit && total >= limit) break;
        Pos p = fromIdx(reps[n]); total++;
        int v = value(p);
        conf[val[n] + 1][v + 1]++;
        if (v != val[n] && bad.size() < 4000) bad.push_back(n);
        // move level: losing moves, drawing in-material moves, drawing exits
        Mv mv[32]; int k = genLegal(p, mv); int nl = 0, di = 0, dc = 0;
        if (k != nm[n]) nmBad++;
        for (int i = 0; i < k; i++) {
            if (mv[i].cap) { if (!favourableCapture(p, mv[i])) dc++; continue; }
            Pos c = applyMv(p, mv[i]); int cv = value(c);
            if (cv == 1) nl++; else if (cv == 0) di++;
        }
        if (v == 0 && (nl != nlose[n] || di != ndi[n] || dc != ndc[n])) posBadMoves++;
        if (v == 0 && nl != nlose[n]) moveBad++;
    }
    printf("mode %d, positions %ld, move generations %ld (%.1f per position)\n", MODE, total, GENCALLS, (double)GENCALLS / total);
    const char* nmv[3] = {"loss", "draw", "win"};
    printf("rows = table, columns = rules\n");
    for (int a = 0; a < 3; a++) { printf("%5s:", nmv[a]); for (int b = 0; b < 3; b++) printf(" %9ld", conf[a][b]); printf("\n"); }
    long agree = conf[0][0] + conf[1][1] + conf[2][2];
    printf("agreement %ld of %ld, disagreements %ld\n", agree, total, total - agree);
    printf("move count disagreements %ld; drawn positions with a differing (losing, drawing, exit) count %ld (losing count only: %ld)\n", nmBad, posBadMoves, moveBad);
    FILE* fo = fopen("bad.txt", "w");
    for (int n : bad) { Pos p = fromIdx(reps[n]); fprintf(fo, "%d %d %d %d %d table=%d rule=%d\n", p.k[0], p.r[0], p.k[1], p.r[1], p.t, val[n], value(p)); }
    fclose(fo);
    return 0;
}

static const char* sqn(int s) { static char b[8][3]; static int i = 0; i = (i + 1) & 7; b[i][0] = 'a' + fl(s); b[i][1] = '1' + rk(s); b[i][2] = 0; return b[i]; }
static void describe(const Pos& p) {
    int m = p.t, o = 1 - m;
    printf("position: mover=%s king %s rook %s ; other king %s rook %s\n", m == 0 ? "white" : "black", sqn(p.k[m]), sqn(p.r[m]), sqn(p.k[o]), sqn(p.r[o]));
    int v = value(p);
    int depth = 0;
    if (v == 1) { for (int h = 1; h <= HW; h += 2) if (Wf(h, p)) { depth = h; break; } }
    if (v == -1) { for (int h = 0; h <= HL; h += 2) if (Lf(h, p)) { depth = h; break; } }
    printf("status: %s", v == 1 ? "mover wins" : v == -1 ? "mover loses" : "drawn"); if (v) printf(" (decided within %d plies)", depth); printf("\n");
    Mv mv[32]; int k = genLegal(p, mv); int nl = 0, di = 0, dc = 0, dw = 0;
    for (int i = 0; i < k; i++) {
        const char* cls;
        if (mv[i].cap) { bool f = favourableCapture(p, mv[i]); cls = f ? "capture of the rook, favourable (win)" : "capture of the rook, not favourable (drawing exit)"; if (!f) dc++; else dw++; }
        else { Pos c = applyMv(p, mv[i]); int cv = value(c); if (cv == 1) { cls = "loses"; nl++; } else if (cv == 0) { cls = "draws"; di++; } else { cls = "wins"; dw++; } }
        printf("  %s%s-%s : %s\n", mv[i].piece == 0 ? "K" : "R", mv[i].piece == 0 ? sqn(p.k[m]) : sqn(p.r[m]), sqn(mv[i].to), cls);
    }
    printf("legal moves %d: losing %d, drawing within the material %d, drawing exits %d, winning %d\n", k, nl, di, dc, dw);
    if (v == 0) printf("boundary reading A (the mover has a losing move): %s\n", nl > 0 ? "yes" : "no");
}
