"""Segment KRvK: a king and a rook against a king. A certified base: its rule is closed-form and is checked against the
table on every position (check_segments.py).

Rule (the mating net on the royal target is decided in this material):
  the side with the rook, to move: wins.
  the bare king, to move: draws if it can capture the rook (the exit to KvK) or is stalemated; otherwise loses.
"""
import chess
from material import register, signature

SIGNATURE = "KRvK"


def status(b):
    rook_side = chess.WHITE if b.pieces(chess.ROOK, chess.WHITE) else chess.BLACK
    if b.turn == rook_side:
        return 1
    moves = list(b.legal_moves)
    if any(b.is_capture(m) for m in moves):
        return 0                                    # the exit: KvK, a dead position
    if not moves and not b.is_check():
        return 0                                    # stalemate
    return -1


register(SIGNATURE, status)
