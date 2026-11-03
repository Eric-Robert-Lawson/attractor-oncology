"""Checks each segment module against the designer's tables (python-chess probes). The modules read no table."""
import sys, random, itertools
sys.path.insert(0, ".")
import chess, chess.syzygy
import material, kvk, krvk, krvkr
TB = sys.argv[1]
tb = chess.syzygy.open_tablebase(TB)

# KRvK: every legal position, both sides to move
n = bad = 0
for wk, wr, bk in itertools.permutations(range(64), 3):
    if chess.square_distance(wk, bk) <= 1:
        continue
    for turn in (True, False):
        b = chess.Board(None)
        b.set_piece_at(wk, chess.Piece(chess.KING, True)); b.set_piece_at(wr, chess.Piece(chess.ROOK, True)); b.set_piece_at(bk, chess.Piece(chess.KING, False))
        b.turn = turn
        if b.is_attacked_by(b.turn, b.king(not b.turn)):
            continue
        w = tb.probe_wdl(b); tv = (w > 0) - (w < 0)
        n += 1
        if krvk.status(b) != tv:
            bad += 1
print("KRvK module against the KRvK table: positions", n, "disagreements", bad, flush=True)

# KRvKR: random positions, read through the interface
N = int(sys.argv[2]); random.seed(5)
n = bad = 0; comp = {1: 0, 0: 0, -1: 0}
while n < N:
    sq = random.sample(range(64), 4)
    b = chess.Board(None)
    for s, pc, c in ((sq[0], chess.KING, True), (sq[1], chess.ROOK, True), (sq[2], chess.KING, False), (sq[3], chess.ROOK, False)):
        b.set_piece_at(s, chess.Piece(pc, c))
    b.turn = random.random() < 0.5
    if b.is_attacked_by(b.turn, b.king(not b.turn)):
        continue
    w = tb.probe_wdl(b); tv = (w > 0) - (w < 0)
    if tv == 0 and comp[0] >= N * 0.5:
        continue
    n += 1; comp[tv] += 1
    if krvkr.status(b) != tv:
        bad += 1; print("MISMATCH", b.fen(), tv, flush=True)
print("KRvKR module (exits read through the KRvK and KvK modules) against the KRvKR table: positions", n, comp, "disagreements", bad, flush=True)
