"""Segment KvK: insufficient material, a dead position. Status 0 for the side to move (the draw)."""
import chess
from material import register

SIGNATURE = "KvK"


def status(b):
    return 0


register(SIGNATURE, status)
