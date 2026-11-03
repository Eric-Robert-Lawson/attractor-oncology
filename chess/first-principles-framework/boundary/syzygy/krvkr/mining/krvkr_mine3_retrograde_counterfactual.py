#!/usr/bin/env python3
# Independent retrograde solve of KRvKR from move structure + KRvK exit values only,
# DTZ distribution of decisive positions, exit-dependence counterfactual, independent move generator check.
import numpy as np, time, random, chess, chess.syzygy
import os
TABLES = os.environ.get("KRVKR_TABLES", "/home/claude/syz")      # the folder with KRvK and KRvKR .rtbw/.rtbz
OUT = os.environ.get("KRVKR_OUT", TABLES + "/out")           # where the mined arrays are written and read
T0 = time.time()
def log(*a): print("[%7.1fs]" % (time.time() - T0), *a, flush=True)
tb = chess.syzygy.open_tablebase(TABLES)
p1 = np.load(OUT + "/pass1.npz"); p2 = np.load(OUT + "/pass2.npz")
reps = p1["reps"]; valid = p1["valid"]; val = p1["val"]; nmoves = p2["nmoves"]
M = len(reps)

# canonical map recomputed
N = 64 ** 4 * 2
idx = np.arange(N, dtype=np.int32)
turn = (idx & 1).astype(np.int32)
r = idx >> 1
br = (r & 63); r >>= 6
bk = (r & 63); r >>= 6
wr = (r & 63); r >>= 6
wk = (r & 63)
def sym(s, k):
    f = (s & 7); rk = (s >> 3)
    nf, nr = [(f, rk), (7 - f, rk), (f, 7 - rk), (7 - f, 7 - rk), (rk, f), (7 - rk, f), (rk, 7 - f), (7 - rk, 7 - f)][k]
    return (nr * 8 + nf).astype(np.int32)
def mk(a, b, c, d, t): return ((((a * 64 + b) * 64 + c) * 64 + d) * 2 + t).astype(np.int32)
canon = idx.copy()
for k in range(8):
    a, b, c, d = sym(wk, k), sym(wr, k), sym(bk, k), sym(br, k)
    canon = np.minimum(canon, mk(a, b, c, d, turn))
    canon = np.minimum(canon, mk(c, d, a, b, 1 - turn))
repid = np.full(N, -1, np.int32); repid[reps[valid]] = np.nonzero(valid)[0].astype(np.int32)
del idx, turn, r, br, bk, wr, wk, a, b, c, d
log("canon ready")

def decode(i):
    t = i & 1; q = i >> 1
    brs = q & 63; q >>= 6; bks = q & 63; q >>= 6; wrs = q & 63; q >>= 6; wks = q & 63
    return [wks, wrs, bks, brs], t
def board_from(sq, t):
    b = chess.Board(None)
    b.set_piece_at(sq[0], chess.Piece(chess.KING, True)); b.set_piece_at(sq[1], chess.Piece(chess.ROOK, True))
    b.set_piece_at(sq[2], chess.Piece(chess.KING, False)); b.set_piece_at(sq[3], chess.Piece(chess.ROOK, False))
    b.turn = (t == 0); return b

# ---- all edges
es = []; ed = []; ee = []   # src, dst (-1 for exit), exit value for mover
dtz_win = [];
for n in np.nonzero(valid)[0]:
    sq, t = decode(int(reps[n])); b = board_from(sq, t)
    for mv in b.legal_moves:
        if b.is_capture(mv):
            b.push(mv)
            w = tb.probe_wdl(b) if not b.is_insufficient_material() else 0
            b.pop()
            ev = -((w > 0) - (w < 0))
            es.append(n); ed.append(-1); ee.append(ev)
        else:
            s = list(sq); s[s.index(mv.from_square)] = mv.to_square
            cid = int(repid[int(canon[int(mk(np.int32(s[0]), np.int32(s[1]), np.int32(s[2]), np.int32(s[3]), np.int32(1 - t)))])])
            assert cid >= 0
            es.append(n); ed.append(cid); ee.append(9)
    if n % 300000 == 0: log("edges", n, len(es))
es = np.array(es, np.int32); ed = np.array(ed, np.int32); ee = np.array(ee, np.int8)
log("edges total", len(es), " exits", int((ed < 0).sum()))
np.savez_compressed(OUT + "/edges.npz", es=es, ed=ed, ee=ee)

incheck = p2["incheck"]
def solve(exit_draw_value):
    """exit_draw_value: value assigned to exits whose table value is a draw (0 baseline, -1 counterfactual)."""
    ev = ee.copy()
    ev[(ed < 0) & (ee == 0)] = exit_draw_value
    L = np.zeros(M, np.int8)
    deg = np.bincount(es, minlength=M)
    # terminal nodes: no legal moves
    term = valid & (deg == 0)
    L[term & incheck] = -1   # mate
    # stalemate: 0 stays but must be fixed: mark with 2
    L[term & ~incheck] = 2
    it = 0
    while True:
        it += 1
        isexit = ed < 0
        # mover value per edge: for in-material edges, -L[child] when child is win(1)/loss(-1); draw-fixed(2) -> 0 ; unknown -> None
        cl = np.zeros(len(es), np.int8)
        cl[~isexit] = L[ed[~isexit]]
        mv = np.full(len(es), 7, np.int8)   # 7 unknown
        mv[isexit] = ev[isexit]
        known = (~isexit) & ((cl == 1) | (cl == -1))
        mv[known] = -cl[known]
        mv[(~isexit) & (cl == 2)] = 0
        win_n = np.zeros(M, bool); win_n[es[mv == 1]] = True
        lose_cnt = np.bincount(es[mv == -1], minlength=M)
        newL = L.copy()
        upd_w = valid & (L == 0) & win_n
        upd_l = valid & (L == 0) & (~win_n) & (lose_cnt == deg) & (deg > 0)
        newL[upd_w] = 1; newL[upd_l] = -1
        ch = int(upd_w.sum() + upd_l.sum())
        L = newL
        if ch == 0: break
    L[L == 2] = 0
    return L, it

L0, it0 = solve(0)
agree = (L0[valid] == val[valid]).mean()
log("baseline retrograde: iterations", it0, " agreement with the table on all valid positions:", agree,
    " disagreements", int((L0[valid] != val[valid]).sum()))
# race depth: iteration at which each position resolved is not stored; recompute depth distribution
def depths():
    deg = np.bincount(es, minlength=M)
    L = np.zeros(M, np.int8); dep = np.zeros(M, np.int16)
    term = valid & (deg == 0)
    L[term & incheck] = -1; dep[term & incheck] = 0; L[term & ~incheck] = 2
    it = 0
    isexit = ed < 0
    while True:
        it += 1
        cl = np.zeros(len(es), np.int8); cl[~isexit] = L[ed[~isexit]]
        mv = np.full(len(es), 7, np.int8); mv[isexit] = ee[isexit]
        known = (~isexit) & ((cl == 1) | (cl == -1)); mv[known] = -cl[known]
        mv[(~isexit) & (cl == 2)] = 0
        win_n = np.zeros(M, bool); win_n[es[mv == 1]] = True
        lose_cnt = np.bincount(es[mv == -1], minlength=M)
        uw = valid & (L == 0) & win_n
        ul = valid & (L == 0) & (~win_n) & (lose_cnt == deg) & (deg > 0)
        if uw.sum() + ul.sum() == 0: break
        L[uw] = 1; dep[uw] = it; L[ul] = -1; dep[ul] = it
    return L, dep
Ld, dep = depths()
log("race-depth (rounds of retrograde to resolve, 1 = immediate: a winning exit/capture of the rook or a move to a lost child already known)")
for sgn, name in ((1, "mover wins"), (-1, "mover loses")):
    m = valid & (Ld == sgn)
    u, c = np.unique(dep[m], return_counts=True)
    print(name, "count", int(m.sum()), "max rounds", int(u.max()))
    cum = 0
    for a, bb in zip(u, c):
        cum += bb
        print("   rounds", int(a), ":", int(bb), " cumulative share", round(cum / m.sum(), 4))
np.save(OUT + "/dep.npy", dep)

# ---- counterfactual: draw-exits count as a loss for whoever takes them
L1, it1 = solve(-1)
v = valid
print("\n== counterfactual: exits that draw are treated as losing for the side that takes them")
for a in (-1, 0, 1):
    for bb in (-1, 0, 1):
        c = int(((L0 == a) & (L1 == bb) & v).sum())
        if c: print("   baseline", {-1: "loss", 0: "draw", 1: "win"}[a], "-> counterfactual", {-1: "loss", 0: "draw", 1: "win"}[bb], ":", c)
d0 = (L0 == 0) & v
print("drawn positions that stay drawn without the exit resource:", int(((L1 == 0) & d0).sum()), "of", int(d0.sum()))
print("drawn positions whose draw depends on an exit:", int(((L1 != 0) & d0).sum()), round(int(((L1 != 0) & d0).sum()) / int(d0.sum()), 4))
np.save(OUT + "/L1.npy", L1)

# ---- DTZ of decisive positions (plies to zeroing, which in KRvKR is a capture)
dz = {}
vals = []
for n in np.nonzero(valid & (val != 0))[0][::7]:   # every 7th decisive position for speed
    sq, t = decode(int(reps[n])); b = board_from(sq, t)
    try: vals.append((abs(tb.probe_dtz(b)), int(val[n]), int(dep[n])))
    except Exception as e: pass
arr = np.array(vals)
log("DTZ sampled positions:", len(arr))
for sgn, name in ((1, "mover wins"), (-1, "mover loses")):
    a = arr[arr[:, 1] == sgn]
    u, c = np.unique(a[:, 0], return_counts=True)
    print(name, "dtz distribution (value:count):", dict(zip(map(int, u[:25]), map(int, c[:25]))), "... max", int(u.max()))
w = arr[arr[:, 1] == 1]
print("correlation of table DTZ with retrograde rounds for won positions:", round(float(np.corrcoef(w[:, 0], w[:, 2])[0, 1]), 3))

# ---- independent move generator vs python-chess counts (step 8)
def indep_moves(sq, t):
    wkq, wrq, bkq, brq = sq
    me_k, me_r, op_k, op_r = (wkq, wrq, bkq, brq) if t == 0 else (bkq, brq, wkq, wrq)
    def fr(s): return s & 7, s >> 3
    def attacked_by_rook(rs, target, occ):
        rf, rr = fr(rs); tf, tr = fr(target)
        if rf != tf and rr != tr: return False
        df = (tf > rf) - (tf < rf); dr = (tr > rr) - (tr < rr)
        f, r_ = rf + df, rr + dr
        while (f, r_) != (tf, tr):
            if (r_ * 8 + f) in occ: return False
            f += df; r_ += dr
        return True
    def kadj(a, b_):
        af, ar = fr(a); bf, br_ = fr(b_); return max(abs(af - bf), abs(ar - br_)) <= 1
    cnt = 0; caps = 0
    # king moves
    kf, kr = fr(me_k)
    for df in (-1, 0, 1):
        for dr in (-1, 0, 1):
            if df == 0 and dr == 0: continue
            f, r_ = kf + df, kr + dr
            if not (0 <= f < 8 and 0 <= r_ < 8): continue
            to = r_ * 8 + f
            if to == me_r: continue
            if kadj(to, op_k): continue
            occ = {op_k, op_r, me_r} - {to}   # king has left me_k; captured piece removed
            if to == op_r:
                capd = True
                if attacked_by_rook(op_r, to, occ) and False: pass
            else: capd = False
            if (not capd) and attacked_by_rook(op_r, to, {op_k, me_r, op_r}): continue
            cnt += 1; caps += capd
    # rook moves
    rf, rr = fr(me_r)
    for df, dr in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        f, r_ = rf + df, rr + dr
        while 0 <= f < 8 and 0 <= r_ < 8:
            to = r_ * 8 + f
            if to == me_k: break
            if to == op_k: break
            capd = (to == op_r)
            occ = {me_k, op_k} | ({op_r} if not capd else set())
            if not capd or True:
                # own king must not be attacked by enemy rook after move (rook now at 'to')
                occ2 = {op_k, to} | {me_k}
                if capd: attacked = False
                else: attacked = attacked_by_rook(op_r, me_k, {op_k, to, me_k})
                if not attacked:
                    cnt += 1; caps += capd
            if capd: break
            f += df; r_ += dr
    return cnt, caps
rng = np.random.default_rng(1)
sample = rng.choice(np.nonzero(valid)[0], 20000, replace=False)
bad = 0
for n in sample:
    sq, t = decode(int(reps[n]))
    c, cp = indep_moves(sq, t)
    if c != nmoves[n] or cp != p2["ncap_total"][n]: bad += 1
log("independent move generator: sample 20000, count disagreements:", bad)
# symmetry fold check: raw random positions (unfolded) value vs representative
bad2 = 0; tested = 0
for _ in range(20000):
    sq = rng.choice(64, 4, replace=False).tolist(); t = int(rng.integers(0, 2))
    b = board_from(sq, t)
    if b.is_attacked_by(b.turn, b.king(not b.turn)): continue
    w = tb.probe_wdl(b); raw = (w > 0) - (w < 0)
    ci = int(canon[int(mk(np.int32(sq[0]), np.int32(sq[1]), np.int32(sq[2]), np.int32(sq[3]), np.int32(t)))])
    tested += 1
    if raw != val[int(repid[ci])]: bad2 += 1
log("fold check: raw positions tested", tested, "value differs from representative:", bad2)
