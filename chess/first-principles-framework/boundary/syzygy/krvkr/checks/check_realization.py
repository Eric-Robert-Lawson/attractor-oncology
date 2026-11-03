"""Checks realization() on sampled drawn positions against the classes of the exhaustive drawing graph (mining step 6)."""
import sys, random
sys.path.insert(0, ".")
import numpy as np, chess
import material, kvk, krvk, krvkr
OUT = sys.argv[1]
p1 = np.load(OUT + "/pass1.npz"); an = np.load(OUT + "/analysis.npz")
reps = p1["reps"]; cls = an["cls"]       # 1 terminal only, 3 mixed
def board(i):
    t = i & 1; q = i >> 1; br = q & 63; q >>= 6; bk = q & 63; q >>= 6; wr = q & 63; q >>= 6; wk = q & 63
    b = chess.Board(None)
    for s, pc, c in ((wk, chess.KING, True), (wr, chess.ROOK, True), (bk, chess.KING, False), (br, chess.ROOK, False)):
        b.set_piece_at(int(s), chess.Piece(pc, c))
    b.turn = (t == 0); return b
rng = np.random.default_rng(2)
tests = [(1, rng.choice(np.nonzero(cls == 1)[0], int(sys.argv[2]), replace=False)), (3, rng.choice(np.nonzero(cls == 3)[0], int(sys.argv[2]), replace=False))]
for c, idxs in tests:
    ok = 0; nodes = []; plies = []
    for n in idxs:
        r = krvkr.realization(board(int(reps[n])))
        want = "terminal only" if c == 1 else "cycle"
        if r[0] == want:
            ok += 1
        if r[0] == "cycle": plies.append(len(r[1]))
        else: nodes.append(r[1])
    print("class", "terminal only" if c == 1 else "mixed", ": sampled", len(idxs), "realization agrees", ok,
          ("; closure sizes of the terminal-only ones: max %d" % max(nodes)) if nodes else "",
          ("; plies to a position with a four-ply cycle: max %d" % max(plies)) if plies else "", flush=True)
