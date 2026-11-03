"""Segment KRvKR: a king and a rook against a king and a rook. The irreversible boundary is a capture.

Depends on: KRvK and KvK (through material.status_below, at the position reached by a capture). Provides: status(board).

The terms and the composition are in krvkr_rule_set.md. No table is read.
"""
import sys
import chess
from functools import lru_cache
from material import register, status_below, signature

SIGNATURE = "KRvKR"
HW, HL = 7, 6          # horizon in plies for a win and for a loss (the scope of the calculation)


def key(b):
    return b.board_fen() + (" w" if b.turn else " b")


def board_of(k):
    b = chess.Board(None); b.set_board_fen(k[:-2]); b.turn = (k[-1] == "w")
    return b


# ---- the exit: the value, for the side that captures, of the position that the capture reaches --------------------
def exit_value(b, mv):
    """Read through the interface from the lower segment (KRvK or KvK). Nothing of the lower segment is used but its status."""
    b.push(mv)
    try:
        return -status_below(b)                     # the status is for the side now to move, the capturer's opponent
    finally:
        b.pop()


def captures(b):
    return [m for m in b.legal_moves if b.is_capture(m)]


def noncaptures(b):
    return [m for m in b.legal_moves if not b.is_capture(m)]


# ---- the winner's candidate pool: checks, and moves that create a threat of a win ---------------------------------
def in_pool(b, mv, h):
    if b.gives_check(mv):
        return True
    b.push(mv)
    b.turn = not b.turn            # the defender passes
    try:
        return W(h - 2, b)         # the mover would win within h-2 plies
    finally:
        b.turn = not b.turn
        b.pop()


# ---- composition ---------------------------------------------------------------------------------------------------
@lru_cache(maxsize=None)
def _W(k, h):
    return _Wb(board_of(k), h)


def W(h, b):
    return _W(key(b), h)


def _Wb(b, h):
    """The mover wins within h plies (h odd)."""
    for mv in captures(b):
        if exit_value(b, mv) == 1:
            return True
    if h >= 3 and W(h - 2, b):
        return True
    for mv in noncaptures(b):
        if h > 1 and not in_pool(b, mv, h):
            continue
        b.push(mv)
        try:
            if L(h - 1, b):
                return True
        finally:
            b.pop()
    return False


@lru_cache(maxsize=None)
def _L(k, h):
    return _Lb(board_of(k), h)


def L(h, b):
    return _L(key(b), h)


def _Lb(b, h):
    """The mover loses within h plies (h even)."""
    if h == 0:
        return b.is_checkmate()
    if L(h - 2, b):
        return True
    moves = list(b.legal_moves)
    if not moves:
        return False
    for mv in moves:
        if b.is_capture(mv):
            if exit_value(b, mv) != -1:
                return False                       # a capture that does not lose is not refuted
        else:
            b.push(mv)
            try:
                if not W(h - 1, b):
                    return False                   # a move that is not refuted by a net
            finally:
                b.pop()
    return True


def status(b):
    if W(HW, b):
        return 1
    if L(HL, b):
        return -1
    return 0


register(SIGNATURE, status)


def move_classes(b):
    """For each legal move: loses, draws, wins (exits are classed by their value)."""
    out = []
    for mv in b.legal_moves:
        san = b.san(mv)
        if b.is_capture(mv):
            out.append((san, {1: "exit: wins", 0: "exit: draws", -1: "exit: loses"}[exit_value(b, mv)]))
        else:
            b.push(mv)
            v = status(b)
            b.pop()
            out.append((san, {1: "loses", 0: "draws", -1: "wins"}[v]))
    return out


# ---- the draw: drawing moves, repetition cycles and the realization ----------------------------------------------------
def drawing_moves(b):
    """The non-capture moves that keep the draw (the child is drawn), and the capture exits that draw."""
    stay, exits = [], []
    for mv in b.legal_moves:
        if b.is_capture(mv):
            if exit_value(b, mv) == 0:
                exits.append(mv)
        else:
            b.push(mv)
            d = status(b) == 0
            b.pop()
            if d:
                stay.append(mv)
    return stay, exits


def shuffle4(b):
    """A four-ply cycle through the position: both sides move a piece and move it back, every position drawn. Returns the line or None."""
    p0 = key(b)
    for m1 in drawing_moves(b)[0]:
        o1 = m1.from_square
        b.push(m1)
        for m2 in drawing_moves(b)[0]:
            o2 = m2.from_square
            b.push(m2)
            r1 = chess.Move(m1.to_square, o1)
            if r1 in b.legal_moves and not b.is_capture(r1):
                b.push(r1)
                if status(b) == 0:
                    r2 = chess.Move(m2.to_square, o2)
                    if r2 in b.legal_moves and not b.is_capture(r2):
                        b.push(r2)
                        if key(b) == p0:
                            line = [m1, m2, r1, r2]
                            b.pop(); b.pop(); b.pop(); b.pop()
                            return line
                        b.pop()
                b.pop()
            b.pop()
        b.pop()
    return None


def realization(b, node_cap=20000):
    """How the draw is sustained, from the position:
       ('cycle', path, line): along the drawing moves 'path' (possibly empty) a position is reached that has the four-ply cycle 'line';
       ('terminal only', n): the closure of the drawing moves is exhausted (n positions) with no cycle in it, so every drawing line ends in an exit;
       ('undetermined', n): the search stopped at node_cap positions."""
    start = key(b)
    seen = {start: None}
    frontier = [start]
    while frontier:
        nxt = []
        for k in frontier:
            c = board_of(k)
            line = shuffle4(c)
            if line is not None:
                path = []
                while seen[k] is not None:
                    k, mv = seen[k]
                    path.append(mv)
                return ("cycle", list(reversed(path)), line)
            for mv in drawing_moves(c)[0]:
                c.push(mv)
                k2 = key(c)
                c.pop()
                if k2 not in seen:
                    seen[k2] = (key(c), mv)
                    nxt.append(k2)
                    if len(seen) >= node_cap:
                        return ("undetermined", len(seen))
        frontier = nxt
    return ("terminal only", len(seen))


if __name__ == "__main__":
    b = chess.Board(None)
    parts = sys.argv[1].split()
    b.set_board_fen(parts[0]); b.turn = (len(parts) < 2 or parts[1] == "w")
    print({1: "mover wins", 0: "drawn", -1: "mover loses"}[status(b)])
    for san, c in move_classes(b):
        print("  ", san, c)
    if status(b) == 0:
        r = realization(b)
        print("realization:", r[0], (r[1:] if r[0] != "cycle" else "reached by %d drawing moves; cycle %s" % (len(r[1]), " ".join(m.uci() for m in r[2]))))
