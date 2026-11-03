"""Material signatures and the registry of certified segments.

A segment is a material together with the irreversible boundary that leaves it (a capture or a promotion in these materials). It exposes
one function, status(board) -> +1 / 0 / -1 for the side to move. A segment above another reads the lower one only through this
function, at the positions its irreversible moves reach (the interface). Nothing else of a lower segment is visible to the one above.
"""
import chess

ORDER = "KQRBNP"
LETTER = {chess.KING: "K", chess.QUEEN: "Q", chess.ROOK: "R", chess.BISHOP: "B", chess.KNIGHT: "N", chess.PAWN: "P"}


def side_string(b, color):
    return "".join(LETTER[pt] * len(b.pieces(pt, color)) for pt in (chess.KING, chess.QUEEN, chess.ROOK, chess.BISHOP, chess.KNIGHT, chess.PAWN))


def signature(b):
    w, k = side_string(b, chess.WHITE), side_string(b, chess.BLACK)
    a, c = sorted((w, k), key=lambda s: (len(s), [ORDER.index(x) for x in s]), reverse=True)
    return a + "v" + c


REGISTRY = {}


def register(sig, status_fn):
    REGISTRY[sig] = status_fn


def status_below(b):
    """The interface: the status, for the side to move, of a position of a lower segment."""
    sig = signature(b)
    if sig not in REGISTRY:
        raise KeyError("no certified segment for " + sig)
    return REGISTRY[sig](b)
