#!/usr/bin/env python3
# Mining process of syzygy_mining.md section 3, run on KRvKR (pawnless, colour-symmetric material).
# Steps: 1 material, 2 enumerate+fold (16 images), 3 probe, 4 classify moves, 5 minimax check,
# 6 drawing-move graph, 7 store. Step 8 is mine2.py.
import sys, time, numpy as np, chess, chess.syzygy
import os
TABLES = os.environ.get("KRVKR_TABLES", "/home/claude/syz")      # the folder with KRvK and KRvKR .rtbw/.rtbz
OUT = os.environ.get("KRVKR_OUT", TABLES + "/out")           # where the mined arrays are written and read
import os; os.makedirs(OUT, exist_ok=True)
T0 = time.time()
def log(*a):
    print("[%7.1fs]" % (time.time() - T0), *a, flush=True)

tb = chess.syzygy.open_tablebase(TABLES)

N = 64 ** 4 * 2
idx = np.arange(N, dtype=np.int32)
turn = (idx & 1).astype(np.int8)
r = idx >> 1
br = (r & 63).astype(np.int8); r >>= 6
bk = (r & 63).astype(np.int8); r >>= 6
wr = (r & 63).astype(np.int8); r >>= 6
wk = (r & 63).astype(np.int8)
del r

def sym(s, k):
    f = (s & 7).astype(np.int16); rk = (s >> 3).astype(np.int16)
    if k == 0: nf, nr = f, rk
    elif k == 1: nf, nr = 7 - f, rk
    elif k == 2: nf, nr = f, 7 - rk
    elif k == 3: nf, nr = 7 - f, 7 - rk
    elif k == 4: nf, nr = rk, f
    elif k == 5: nf, nr = 7 - rk, f
    elif k == 6: nf, nr = rk, 7 - f
    else: nf, nr = 7 - rk, 7 - f
    return (nr * 8 + nf).astype(np.int32)

def mk(a, b, c, d, t):
    return ((((a * 64 + b) * 64 + c) * 64 + d) * 2 + t).astype(np.int32)

log("computing canonical index over 16 images")
canon = idx.copy()
for k in range(8):
    a, b, c, d = sym(wk, k), sym(wr, k), sym(bk, k), sym(br, k)
    canon = np.minimum(canon, mk(a, b, c, d, turn.astype(np.int32)))             # no colour swap
    canon = np.minimum(canon, mk(c, d, a, b, (1 - turn).astype(np.int32)))       # colour swap
distinct = (wk != wr) & (wk != bk) & (wk != br) & (wr != bk) & (wr != br) & (bk != br)
isrep = (canon == idx) & distinct
reps = np.nonzero(isrep)[0].astype(np.int32)
log("raw distinct-square positions:", int(distinct.sum()), " folded representatives:", len(reps))
del a, b, c, d

def board_from(i):
    t = i & 1; q = i >> 1
    brs = q & 63; q >>= 6
    bks = q & 63; q >>= 6
    wrs = q & 63; q >>= 6
    wks = q & 63
    b = chess.Board(None)
    b.set_piece_at(int(wks), chess.Piece(chess.KING, True))
    b.set_piece_at(int(wrs), chess.Piece(chess.ROOK, True))
    b.set_piece_at(int(bks), chess.Piece(chess.KING, False))
    b.set_piece_at(int(brs), chess.Piece(chess.ROOK, False))
    b.turn = (t == 0)
    return b, [int(wks), int(wrs), int(bks), int(brs)], t

# pass 1: legality and value
repid = np.full(N, -1, dtype=np.int32)
val = np.zeros(len(reps), dtype=np.int8)
valid = np.zeros(len(reps), dtype=bool)
rawwdl = np.zeros(len(reps), dtype=np.int8)
for n, i in enumerate(reps):
    b, sq, t = board_from(int(i))
    opp = not b.turn
    if b.is_attacked_by(b.turn, b.king(opp)):
        continue   # side not to move in check, or kings adjacent
    valid[n] = True
    w = tb.probe_wdl(b)
    rawwdl[n] = w
    val[n] = (w > 0) - (w < 0)
    repid[i] = n
    if n % 200000 == 0: log("pass1", n, "of", len(reps))
log("valid folded positions:", int(valid.sum()))
np.savez_compressed(OUT + "/pass1.npz", reps=reps, valid=valid, val=val, rawwdl=rawwdl)

# pass 2: move classification
M = len(reps)
nmoves = np.zeros(M, np.int16); nlose = np.zeros(M, np.int16); ndraw_in = np.zeros(M, np.int16)
ndraw_cap = np.zeros(M, np.int16); ndraw_stale = np.zeros(M, np.int16); nwin = np.zeros(M, np.int16)
ncap_total = np.zeros(M, np.int16)
incheck = np.zeros(M, bool)
mismatch = 0; checked = 0
src = []; dst = []
for n in np.nonzero(valid)[0]:
    i = int(reps[n])
    b, sq, t = board_from(i)
    mover_white = (t == 0)
    incheck[n] = b.is_check()
    best = -2
    nm = 0
    for mv in b.legal_moves:
        nm += 1
        if b.is_capture(mv):
            b.push(mv)
            w = tb.probe_wdl(b) if not b.is_insufficient_material() else 0
            stalemate_child = b.is_stalemate()
            b.pop()
            cv = (w > 0) - (w < 0)
            mv_val = -cv
            ncap_total[n] += 1
            if mv_val == 0:
                ndraw_cap[n] += 1
                if stalemate_child: ndraw_stale[n] += 1
        else:
            # child raw index
            s = list(sq)
            pi = s.index(mv.from_square)
            s[pi] = mv.to_square
            ci = mk(np.int32(s[0]), np.int32(s[1]), np.int32(s[2]), np.int32(s[3]), np.int32(1 - t))
            ci = int(ci)
            cc = int(canon[ci])
            cid = int(repid[cc])
            assert cid >= 0, ("child invalid", i, mv)
            mv_val = -int(val[cid])
            if mv_val == 0:
                ndraw_in[n] += 1
                src.append(n); dst.append(cid)
        if mv_val < 0: nlose[n] += 1
        elif mv_val > 0: nwin[n] += 1
        best = max(best, mv_val)
    nmoves[n] = nm
    if nm == 0:
        exp = -1 if b.is_check() else 0
    else:
        exp = best
    checked += 1
    if exp != val[n]:
        mismatch += 1
    if n % 200000 == 0: log("pass2", n, "mism", mismatch)
log("step5: positions checked", checked, "minimax mismatches", mismatch)
np.savez_compressed(OUT + "/pass2.npz", nmoves=nmoves, nlose=nlose, ndraw_in=ndraw_in, ndraw_cap=ndraw_cap,
                    ndraw_stale=ndraw_stale, nwin=nwin, ncap_total=ncap_total, incheck=incheck,
                    src=np.array(src, np.int32), dst=np.array(dst, np.int32))
log("done pass2")
