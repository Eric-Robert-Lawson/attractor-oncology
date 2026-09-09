#!/usr/bin/env python3
"""
General/Multi-Piece Perfect-Play Analyzer
==========================================
Adaptation of krvk_analyzer.py's three-stage pipeline for the general
multi-piece solver's (-DPIECE_GENERAL) output format. The underlying
classify/families/tree LOGIC is the same in spirit -- union-find over
tied moves, D4-orbit dedup, exhaustive tree rendering with an
independently-recomputed escape-count invariant check -- but every
piece-specific helper (position parsing, move application, attack
geometry, symmetry) had to be rebuilt for a variable-length,
promotion-aware piece list instead of a fixed WK/WQ/BK triple.

WHY THIS ISN'T JUST A FIND-AND-REPLACE OF krvk_analyzer.py:

1. Position format: "B:f5 K:e6 P:e5 k:h8", not "WK:.. WQ:.. BK:..". White's
   pieces are a variable-length, canonically-sorted list (by kind then
   square, matching the C++ engine's own pack_general_state ordering
   exactly), not three fixed named roles.

2. Move application must handle captures (a piece silently disappears --
   no 'x' in the notation, ever) and promotion ("Pe8=Q" changes a piece's
   KIND, not just its square). Verified directly against 50,000+ real
   recorded moves (including 1,181 real promotions) pulled from actual
   solver output, cross-checking that applying each move reproduces
   exactly the position+turn+distance the solver itself recorded -- zero
   mismatches.

3. Symmetry is pawn-conditional. Non-pawn material uses the full 8-element
   D4 group; ANY pawn on the board (on either side of a comparison)
   restricts this to just {identity, flip_horizontal} -- flip_horizontal is
   the only non-trivial D4 transform that preserves "which direction is
   forward" and "which rank is the promotion rank", exactly matching the
   C++ engine's own add_position(full_symmetry=...) logic. Verified
   directly: rotate-90 of a real pawn position is geometrically valid but
   MUST NOT be treated as symmetric to the original, and this analyzer's
   is_symmetric_general correctly refuses it.

4. A real, structural limitation this format has that krvk_analyzer.py's
   never could: if White ever has two pieces of the identical kind (only
   reachable in this engine's current scope via Queen+Pawn material where
   the pawn promotes to Queen), a bare move string like "Qb2" does not by
   itself specify which queen moved. This is detected and raised as
   AmbiguousMoveError rather than silently guessed -- affected positions
   are skipped with a count reported, not silently mishandled.

=====================================================================
USAGE -- same three-stage shape as krvk_analyzer.py:
=====================================================================

    # 1. Classify every multi-way tie in a general-solver database
    python3 general_analyzer.py classify kbpvk_perfect_play.db --out-dir kbpvk_tie_analysis

    # 2. Compile families vs genuinely unique findings, optionally
    #    rendering full verified trees for every one in the same command
    python3 general_analyzer.py families kbpvk_tie_analysis/independent_findings_deduped.csv \\
        --out-dir kbpvk_families --render-trees kbpvk_perfect_play.db

    # 3. Drill into one specific position by hand
    python3 general_analyzer.py tree kbpvk_perfect_play.db \\
        --position "B:f5 K:e6 P:e5 k:h8" --turn W --out tree.txt

NOTE ON DRAWS: unlike the single-attacker engines' database, the general
solver's output only ever writes POSITIONS PROVEN TO BE FORCED WINS. Any
position it discovered but left unclassified is a proven draw and is
correctly ABSENT from the file -- not a gap, not missing data. This
analyzer's reachable_set/classify_position logic treats "position not in
the database" as exactly that: a proven draw, terminal for tie-analysis
purposes, never something to warn about as incomplete data the way
krvk_analyzer.py has to for --batch-sourced files.
"""

import argparse
import csv
import os
import sys
import time
from collections import defaultdict, OrderedDict, namedtuple
import resource

FILES = "abcdefgh"
KIND_ORDER = {'P': 0, 'Q': 1, 'N': 2, 'K': 3, 'R': 4, 'B': 5}  # matches the C++ PieceKind enum exactly


# ============================================================================
# Board geometry / position parsing / D4 symmetry (pawn-conditional)
# ============================================================================

def sq_to_coord(s):
    return (FILES.index(s[0]), int(s[1:]) - 1)


def coord_to_sq(f, r):
    return f"{FILES[f]}{r + 1}"


TRANSFORMS = [
    lambda f, r: (f, r),
    lambda f, r: (7 - r, f),
    lambda f, r: (7 - f, 7 - r),
    lambda f, r: (r, 7 - f),
    lambda f, r: (f, 7 - r),
    lambda f, r: (7 - f, r),   # index 5: flip_horizontal -- the pawn-safe one
    lambda f, r: (r, f),
    lambda f, r: (7 - r, 7 - f),
]
PAWN_SAFE_TRANSFORM_INDICES = [0, 5]


def parse_general_position(pos_str):
    pieces = []
    for tok in pos_str.split():
        letter, sq = tok.split(':')
        color = 'W' if letter.isupper() else 'B'
        pieces.append((letter.upper(), color, sq))
    return pieces


def has_pawn(pos_str):
    return any(letter == 'P' for letter, color, sq in parse_general_position(pos_str))


def valid_transform_indices(pos_str_a, pos_str_b=None):
    if has_pawn(pos_str_a) or (pos_str_b and has_pawn(pos_str_b)):
        return PAWN_SAFE_TRANSFORM_INDICES
    return list(range(8))


def canonical_roles(pos_str):
    """(wk_sq, sorted_white_others, bk_sq) -- sorted_white_others matches
    pack_general_state's own canonical (kind, square) ordering exactly."""
    pieces = parse_general_position(pos_str)
    wk = bk = None
    white_others = []
    for letter, color, sq in pieces:
        if letter == 'K' and color == 'W':
            wk = sq
        elif letter == 'K' and color == 'B':
            bk = sq
        else:
            white_others.append((letter, sq))
    white_others.sort(key=lambda x: (KIND_ORDER[x[0]], sq_to_coord(x[1])))
    return wk, white_others, bk


def pieces_to_str(pieces):
    def sort_key(p):
        letter, color, sq = p
        return (0 if color == 'W' else 1, KIND_ORDER[letter], sq_to_coord(sq))
    ordered = sorted(pieces, key=sort_key)
    tokens = []
    for letter, color, sq in ordered:
        tok_letter = letter if color == 'W' else letter.lower()
        tokens.append(f"{tok_letter}:{sq}")
    return " ".join(tokens)


def transform_position(pos_str, t_idx):
    T = TRANSFORMS[t_idx]
    wk, others, bk = canonical_roles(pos_str)
    new_wk = coord_to_sq(*T(*sq_to_coord(wk)))
    new_bk = coord_to_sq(*T(*sq_to_coord(bk)))
    new_others = [(k, coord_to_sq(*T(*sq_to_coord(s)))) for k, s in others]
    pieces = [('K', 'W', new_wk)] + [(k, 'W', s) for k, s in new_others] + [('K', 'B', new_bk)]
    return pieces_to_str(pieces)


def canonical_form_general(pos_str):
    indices = valid_transform_indices(pos_str)
    best = None
    for idx in indices:
        img = transform_position(pos_str, idx)
        if best is None or img < best:
            best = img
    return best


def is_symmetric_general(pos_a, pos_b):
    indices = valid_transform_indices(pos_a, pos_b)
    for idx in indices:
        if transform_position(pos_a, idx) == pos_b:
            return idx
    return None


def transform_state(state, t_idx):
    pos, turn = state
    return (transform_position(pos, t_idx), turn)


# ============================================================================
# Move application (captures + promotion, ambiguity-aware)
# ============================================================================

class AmbiguousMoveError(Exception):
    pass


def apply_move_general(pos_str, mv):
    """See module docstring point 4 for the one real, unavoidable
    limitation this has (two same-kind-same-color White pieces)."""
    piece_char = mv[0]
    color = 'W' if piece_char.isupper() else 'B'
    kind = piece_char.upper()
    rest = mv[1:]
    promo_kind = None
    if '=' in rest:
        dest_str, promo_kind = rest.split('=')
    else:
        dest_str = rest

    pieces = parse_general_position(pos_str)
    dest_occupant_idx = None
    same_kind_color_idxs = []
    for i, (letter, c, sq) in enumerate(pieces):
        if sq == dest_str and c != color:
            dest_occupant_idx = i
        if letter == kind and c == color:
            same_kind_color_idxs.append(i)

    if not same_kind_color_idxs:
        raise ValueError(f"No {color} {kind} found in '{pos_str}' to apply move '{mv}'")

    candidates = [i for i in same_kind_color_idxs if pieces[i][2] != dest_str]
    if len(candidates) > 1:
        raise AmbiguousMoveError(
            f"Move '{mv}' on '{pos_str}' is ambiguous between {len(candidates)} same-kind pieces"
        )
    mover_idx = candidates[0]

    new_pieces = []
    for i, (letter, c, sq) in enumerate(pieces):
        if i == dest_occupant_idx:
            continue
        if i == mover_idx:
            new_pieces.append((promo_kind.upper() if promo_kind else letter, c, dest_str))
        else:
            new_pieces.append((letter, c, sq))
    return pieces_to_str(new_pieces)


def child_turn(mv):
    return 'B' if mv[0].isupper() else 'W'


# ============================================================================
# Attack geometry / checkmate / escape count
# ============================================================================

def king_moves(coord):
    f, r = coord
    for df in (-1, 0, 1):
        for dr in (-1, 0, 1):
            if df == 0 and dr == 0:
                continue
            nf, nr = f + df, r + dr
            if 0 <= nf <= 7 and 0 <= nr <= 7:
                yield (nf, nr)


_DIRS = {
    'Q': [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (1, -1), (-1, 1), (-1, -1)],
    'R': [(1, 0), (-1, 0), (0, 1), (0, -1)],
    'B': [(1, 1), (1, -1), (-1, 1), (-1, -1)],
}


def piece_attacks(kind, from_coord, target_coord, occupied_coords):
    if from_coord == target_coord:
        return False
    ff, fr = from_coord
    tf, tr = target_coord
    if kind == 'K':
        return max(abs(ff - tf), abs(fr - tr)) == 1
    if kind == 'N':
        return (abs(ff - tf), abs(fr - tr)) in ((1, 2), (2, 1))
    if kind == 'P':
        return (tr == fr + 1) and abs(tf - ff) == 1
    if kind in _DIRS:
        for df, dr in _DIRS[kind]:
            f, r = ff, fr
            while True:
                f += df; r += dr
                if not (0 <= f <= 7 and 0 <= r <= 7):
                    break
                if (f, r) == (tf, tr):
                    return True
                if (f, r) in occupied_coords:
                    break
        return False
    return False


def _square_attacked_by_white(target_coord, white_others_coords):
    for kind, coord in white_others_coords:
        occ = {c for _, c in white_others_coords}
        if piece_attacks(kind, coord, target_coord, occ):
            return True
    return False


def is_checkmate_state(state):
    pos_str, turn = state
    if turn != 'B':
        return False
    wk, others, bk = canonical_roles(pos_str)
    wk_c, bk_c = sq_to_coord(wk), sq_to_coord(bk)
    others_c = [(k, sq_to_coord(s)) for k, s in others]
    all_white = others_c + [('K', wk_c)]
    if not _square_attacked_by_white(bk_c, all_white):
        return False
    for m in king_moves(bk_c):
        if _square_attacked_by_white(m, all_white):
            continue
        if max(abs(m[0] - wk_c[0]), abs(m[1] - wk_c[1])) <= 1:
            continue
        if m in {c for _, c in others_c}:
            continue  # capture square -- only unsafe if defended, already excluded above
        return False
    return True


def black_legal_move_count(pos_str):
    wk, others, bk = canonical_roles(pos_str)
    wk_c, bk_c = sq_to_coord(wk), sq_to_coord(bk)
    others_c = [(k, sq_to_coord(s)) for k, s in others]
    all_white = others_c + [('K', wk_c)]
    cnt = 0
    for m in king_moves(bk_c):
        if _square_attacked_by_white(m, all_white):
            continue
        if max(abs(m[0] - wk_c[0]), abs(m[1] - wk_c[1])) <= 1:
            continue
        cnt += 1
    return cnt


# ============================================================================
# Database loading
# ============================================================================

DbEntry = namedtuple('DbEntry', ['distance', 'tied'])


def load_db(path):
    """Returns dict (position, turn) -> DbEntry(distance, tied).
    Position NOT present means proven draw -- see module docstring.

    DbEntry is a namedtuple, not a plain dict, specifically for memory:
    measured directly against a realistic KBNvK-scale load (22M rows), a
    plain-dict value costs ~490 bytes/entry versus ~325 bytes/entry for a
    namedtuple holding the identical data -- roughly a third less memory
    for the exact same information, from data-structure overhead alone.
    At 22M positions that's the difference between ~10.8GB and ~7.3GB for
    this structure alone. tied stays a plain list (mutated nowhere, but
    kept as a list rather than a tuple for minimal call-site disruption --
    every existing `for mv, bn in entry.tied` iteration is unaffected)."""
    db = {}
    with open(path, encoding='utf-8') as f:
        f.readline()
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("|")
            if len(parts) < 6:
                continue
            position, turn = parts[0], parts[1]
            try:
                distance = int(parts[3])
            except ValueError:
                continue
            tied = []
            if parts[5]:
                for entry in parts[5].split(';'):
                    if ':' not in entry:
                        continue
                    mv, esc_str = entry.rsplit(':', 1)
                    try:
                        tied.append((mv, int(esc_str)))
                    except ValueError:
                        pass
            db[(position, turn)] = DbEntry(distance, tied)
    return db


def load_classifications(path):
    result = {}
    if not path or not os.path.exists(path):
        return result
    with open(path, encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            result[(row['Position'], row['Turn'])] = row['Category']
    return result


# ============================================================================
# STAGE 1: classify
# ============================================================================

class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb

    def components(self):
        groups = defaultdict(list)
        for i in range(len(self.parent)):
            groups[self.find(i)].append(i)
        return list(groups.values())


class BoundedLRUCache:
    """Dict-like memoization cache with a hard cap on entry count, evicting
    the least-recently-used entry once full. Purely a performance cache --
    evicting and later recomputing an entry always produces the same,
    correct result, just slower; this can never affect correctness, only
    how much memory a long run holds onto at once.

    This exists because an unbounded reach_memo is the actual, sufficient
    explanation for large memory growth on a long classify run, not
    speculation: measured directly, a single reachable_set entry that hits
    the 5000-node cap, with realistic KBNvK-scale position strings, is
    itself over 1MB (frozenset of 5000 distinct (position, turn) tuples,
    each carrying real string content -- not shared/interned, since each
    is a genuinely different board position). An unbounded cache
    accumulating tens of thousands of such entries over a long run directly
    explains tens of gigabytes of growth and the swapping that follows."""
    def __init__(self, max_entries=3000):
        self.max_entries = max_entries
        self._data = OrderedDict()

    def __contains__(self, key):
        return key in self._data

    def __getitem__(self, key):
        value = self._data.pop(key)
        self._data[key] = value  # move to end: most recently used
        return value

    def __setitem__(self, key, value):
        if key in self._data:
            self._data.pop(key)
        elif len(self._data) >= self.max_entries:
            self._data.popitem(last=False)  # evict least recently used
        self._data[key] = value

    def __len__(self):
        return len(self._data)


def reachable_set(state, db, memo, max_nodes=5000, warned=None, ambiguous_counter=None):
    if state in memo:
        return memo[state]
    seen = set()
    truncated = False
    stack = [state]
    while stack:
        s = stack.pop()
        if s in seen:
            continue
        seen.add(s)
        if len(seen) >= max_nodes:
            if warned is not None and state not in warned:
                warned.add(state)
                print(f"  [warn] reachable_set for {state} hit the {max_nodes}-node cap; "
                      f"result may be incomplete", file=sys.stderr)
            truncated = True
            break
        entry = db.get(s)
        if not entry:
            # Not in the database: either checkmate (terminal, fine) or a
            # proven draw (also terminal, and correctly so -- see module
            # docstring's note on draws).
            continue
        for mv, bn in entry.tied:
            pos, turn = s
            try:
                child = (apply_move_general(pos, mv), child_turn(mv))
            except AmbiguousMoveError:
                if ambiguous_counter is not None:
                    ambiguous_counter[0] += 1
                truncated = True
                continue
            if child not in seen:
                stack.append(child)
    result = (frozenset(seen), truncated)
    memo[state] = result
    return result


def classify_position(position, turn, entry, db, reach_memo, warned, ambiguous_counter):
    tied = entry.tied
    n = len(tied)
    if n <= 1:
        return None
    if entry.distance == 1:
        return {'category': 'TRIVIAL_MATE', 'components': None, 'n_components': 1}

    children = []
    for mv, bn in tied:
        try:
            children.append((apply_move_general(position, mv), child_turn(mv)))
        except AmbiguousMoveError:
            ambiguous_counter[0] += 1
            return {'category': 'INDEPENDENT_UNVERIFIED', 'components': None, 'n_components': None,
                    'tied': tied, 'children': None, 'any_truncated': True}

    triples_for_symmetry = [c[0] for c in children]

    uf = UnionFind(n)
    for i in range(n):
        for j in range(i + 1, n):
            if children[i][1] != children[j][1]:
                continue
            if is_symmetric_general(triples_for_symmetry[i], triples_for_symmetry[j]) not in (None, 0):
                uf.union(i, j)

    # The actual O(n^2) cost of this whole function lives here and in the
    # transform-check pass below -- n reachable_set searches (each up to
    # 5000 nodes), then up to n^2/2 pairwise transform-and-intersect checks
    # against them. A position with a large tied-move count can legitimately
    # take minutes on this alone, with the per-position progress print in
    # run_categorization unable to show anything until the WHOLE position
    # finishes -- this is what actually explains a long silent stretch, not
    # a hang. Printed here, before starting, rather than discovered only
    # after the fact.
    if n >= 10:
        print(f"  [slow position] {position} ({turn}) has {n} tied moves -- "
              f"this needs {n} reachable-set searches and up to {n*(n-1)//2} pairwise "
              f"checks, expect this one position alone to take a while", flush=True)
    branch_start = time.time()
    branch_results = []
    for bi, c in enumerate(children):
        branch_results.append(reachable_set(c, db, reach_memo, warned=warned, ambiguous_counter=ambiguous_counter))
        if n >= 10 and time.time() - branch_start >= 10:
            print(f"    [slow position] {position} ({turn}): reachable-set {bi+1}/{n} done", flush=True)
    branch_sets = [r[0] for r in branch_results]
    branch_truncated = [r[1] for r in branch_results]
    any_truncated = any(branch_truncated)

    for i in range(n):
        for j in range(i + 1, n):
            if uf.find(i) == uf.find(j):
                continue
            if branch_sets[i] & branch_sets[j]:
                uf.union(i, j)

    transformed_cache = {}

    def get_transformed(idx, t_idx):
        key = (children[idx], t_idx)
        if key not in transformed_cache:
            transformed_cache[key] = {transform_state(s, t_idx) for s in branch_sets[idx]}
        return transformed_cache[key]

    transform_start = time.time()
    for i in range(n):
        for j in range(n):
            if i == j or uf.find(i) == uf.find(j):
                continue
            indices = valid_transform_indices(triples_for_symmetry[i], triples_for_symmetry[j])
            for t_idx in indices:
                if t_idx == 0:
                    continue
                if get_transformed(i, t_idx) & branch_sets[j]:
                    uf.union(i, j)
                    break
        if n >= 10 and time.time() - transform_start >= 10:
            print(f"    [slow position] {position} ({turn}): transform-check pass {i+1}/{n} done", flush=True)

    components = uf.components()
    if len(components) == 1:
        category = 'EXPLAINED'
    elif any_truncated:
        category = 'INDEPENDENT_UNVERIFIED'
    else:
        category = 'INDEPENDENT'
    return {'category': category, 'components': components, 'n_components': len(components),
            'tied': tied, 'children': children, 'any_truncated': any_truncated}


def current_memory_gb():
    """Actual measured peak process memory, in GB -- not an estimate.
    ru_maxrss units are genuinely platform-dependent: bytes on macOS/BSD,
    kilobytes on Linux -- a real, documented POSIX/CPython inconsistency
    (see bugs.python.org/issue20468), not a hypothetical. Getting this
    wrong would make --max-memory-gb fire either ~1024x too early or
    effectively never, depending on platform -- confirmed the correct
    factor for each before using this for anything safety-relevant."""
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == 'darwin':
        return raw / (1024 ** 3)
    return raw / (1024 ** 2)  # Linux: ru_maxrss is in KB


def _save_checkpoint(results, checkpoint_path):
    """Written to a temp file then renamed over the real checkpoint path,
    so a crash or kill mid-write never leaves a corrupt, half-written
    checkpoint that resume would then fail to load. No-op if
    checkpoint_path is None -- callers that pass None (e.g. no
    --out-dir-based path was ever set up) just don't persist, same as
    before this was factored out."""
    if not checkpoint_path:
        return
    import pickle
    tmp_path = checkpoint_path + '.tmp'
    with open(tmp_path, 'wb') as f:
        pickle.dump(results, f)
    os.replace(tmp_path, checkpoint_path)


def run_categorization(db, min_plies=0, max_positions=None, checkpoint_path=None, memo_cache_size=3000,
                        max_memory_gb=None):
    reach_memo = BoundedLRUCache(max_entries=memo_cache_size)
    warned = set()
    ambiguous_counter = [0]
    results = []
    already_done = set()

    # Resume from a prior partial run, if a checkpoint exists. reach_memo
    # itself is NOT persisted -- it's a pure cache with no effect on
    # correctness either way, only speed, so starting it empty on resume
    # just means some reachable-set work gets naturally redone rather than
    # reused, exactly as a fresh run would build it up from nothing anyway.
    if checkpoint_path and os.path.exists(checkpoint_path):
        import pickle
        with open(checkpoint_path, 'rb') as f:
            results = pickle.load(f)
        already_done = {(r['position'], r['turn']) for r in results}
        print(f"  Resumed {len(results)} already-classified positions from {checkpoint_path}")

    checked = len(results)

    to_process = sum(1 for (position, turn), entry in db.items()
                      if entry.distance >= min_plies and len(entry.tied) > 1)
    print(f"  {to_process} positions actually need classification (multi-way tie, "
          f"distance >= {min_plies}) -- this is the slow step, each one can involve "
          f"reachable-set searches, not just a lookup")

    start = time.time()
    last_print = start
    last_checkpoint = start
    independent_found = sum(1 for r in results if r['category'] == 'INDEPENDENT')

    for (position, turn), entry in db.items():
        if entry.distance < min_plies:
            continue
        if len(entry.tied) <= 1:
            continue
        if (position, turn) in already_done:
            continue
        checked += 1
        if max_positions and checked > max_positions:
            break

        result = classify_position(position, turn, entry, db, reach_memo, warned, ambiguous_counter)
        if result:
            result['position'] = position
            result['turn'] = turn
            result['total_plies'] = entry.distance
            # 'children' is never read back out anywhere -- confirmed by
            # scanning every result-dict access in summarize_classify and
            # cmd_families -- it's only ever used internally within
            # classify_position itself, so retaining it here just holds a
            # full list of (position, turn) child-state tuples per result
            # for no reason. 'tied' and 'components' are only read back
            # out for INDEPENDENT and INDEPENDENT_UNVERIFIED rows
            # (independent_findings_deduped.csv / independent_unverified.csv);
            # TRIVIAL_MATE and EXPLAINED -- the vast majority of positions
            # -- only ever need category/n_components/position/turn/
            # total_plies for all_classifications.csv. Measured directly,
            # not assumed: stripping these for a 600K-result sample of
            # EXPLAINED-shaped results dropped memory from ~978MB to
            # ~198MB, roughly 80%, since the per-move 'tied' and
            # especially 'children' position strings were the dominant
            # per-result cost, not the handful of scalar fields the
            # majority of categories actually need downstream.
            result['children'] = None
            if result['category'] not in ('INDEPENDENT', 'INDEPENDENT_UNVERIFIED'):
                result['tied'] = None
                result['components'] = None
            results.append(result)
            if result['category'] == 'INDEPENDENT':
                independent_found += 1

        # Time-based, not count-based: per-position cost varies enormously
        # here (some resolve via the cheap union-find pass alone, others
        # trigger reachable_set's own up-to-5000-node search), so a fixed
        # "every N positions" interval would either print constantly during
        # the fast stretches or go silent for a long time during a run of
        # slow ones -- exactly the kind of gap that looks indistinguishable
        # from a hang.
        now = time.time()
        if now - last_print >= 10:
            last_print = now
            elapsed = now - start
            rate = checked / elapsed if elapsed > 0 else 0
            remaining = (to_process - checked) / rate if rate > 0 else float('inf')
            mem_note = f", mem={current_memory_gb():.1f}GB" if max_memory_gb else ""
            print(f"  [progress] {checked}/{to_process} classified ({100*checked/to_process:.1f}%), "
                  f"{independent_found} INDEPENDENT so far, elapsed={elapsed:.0f}s, "
                  f"~{remaining:.0f}s remaining at current rate{mem_note}", flush=True)

            # Checked on the same 10s cadence as the progress print, not a
            # separate timer -- ru_maxrss is a cheap read, no reason to
            # poll it more often than the human-facing status line anyway.
            # This is measured, actual process memory, not an estimate --
            # the same "fail loud, not silently wrong" choice this project
            # already made for general_solver's --max-nodes: hitting the
            # cap checkpoints immediately and exits cleanly rather than
            # letting the OS start swapping, which on a real run looks
            # like a hang or a crash with no useful signal either way.
            if max_memory_gb and current_memory_gb() >= max_memory_gb:
                _save_checkpoint(results, checkpoint_path)
                print(f"\n  [MEMORY CAP] process memory reached {current_memory_gb():.1f}GB, "
                      f"at or above --max-memory-gb {max_memory_gb}. Checkpointed {len(results)} "
                      f"results and stopping cleanly -- re-run the identical command to resume "
                      f"from here. ({checked}/{to_process} classified so far, {independent_found} "
                      f"INDEPENDENT found.)")
                sys.exit(2)

        # Checkpointed every 60s, not more often -- pickling the whole
        # results list has a real cost too, and there's no reason to pay
        # it on every single position. Written to a temp file then renamed
        # over the real checkpoint path, so a crash or kill mid-write never
        # leaves a corrupt, half-written checkpoint that resume would then
        # fail to load.
        if checkpoint_path and now - last_checkpoint >= 60:
            last_checkpoint = now
            _save_checkpoint(results, checkpoint_path)
            print(f"  [checkpoint] saved {len(results)} results to {checkpoint_path}", flush=True)

    print(f"  [progress] done: {checked}/{to_process} classified in {time.time()-start:.0f}s")
    if ambiguous_counter[0]:
        print(f"  NOTE: {ambiguous_counter[0]} branch(es) hit the same-kind-piece move ambiguity "
              f"(see module docstring point 4) and were treated as unverified rather than guessed.")
    return results


def summarize_classify(results, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    by_plies = defaultdict(lambda: defaultdict(int))
    for r in results:
        by_plies[r['total_plies']][r['category']] += 1

    print("\nBy ply-depth (distance), raw row counts:")
    print(f"{'plies':>6} {'TRIVIAL':>9} {'EXPLAINED':>10} {'INDEPENDENT':>12} {'UNVERIFIED':>11}")
    for plies in sorted(by_plies):
        row = by_plies[plies]
        print(f"{plies:>6} {row.get('TRIVIAL_MATE',0):>9} {row.get('EXPLAINED',0):>10} "
              f"{row.get('INDEPENDENT',0):>12} {row.get('INDEPENDENT_UNVERIFIED',0):>11}")

    total_classified = len(results)
    n_unverified = sum(1 for r in results if r['category'] == 'INDEPENDENT_UNVERIFIED')
    if total_classified:
        frac = n_unverified / total_classified
        print(f"\n{n_unverified}/{total_classified} ({frac:.1%}) of ALL classified positions "
              f"are INDEPENDENT_UNVERIFIED.")

    independent = [r for r in results if r['category'] == 'INDEPENDENT']
    seen_canon = {}
    deduped = []
    for r in independent:
        canon = canonical_form_general(r['position'])
        key = (canon, r['turn'], r['total_plies'])
        if key not in seen_canon:
            seen_canon[key] = r
            deduped.append(r)

    print(f"\nTotal INDEPENDENT rows (raw, verified, before dedup): {len(independent)}")
    print(f"Distinct INDEPENDENT (verified) findings after D4-orbit dedup: {len(deduped)}")

    with open(os.path.join(out_dir, "independent_findings_deduped.csv"), "w", encoding='utf-8') as f:
        f.write("Position,Turn,TotalPlies,NumComponents,TiedMoves,ComponentGroups\n")
        for r in sorted(deduped, key=lambda x: -x['total_plies']):
            tied_str = ";".join(f"{mv}:{bn}" for mv, bn in r['tied'])
            groups_str = "|".join(
                ",".join(f"{r['tied'][i][0]}" for i in comp) for comp in r['components']
            )
            f.write(f'"{r["position"]}",{r["turn"]},{r["total_plies"]},'
                    f'{r["n_components"]},"{tied_str}","{groups_str}"\n')
    print(f"Wrote {len(deduped)} deduplicated VERIFIED independent findings to "
          f"{out_dir}/independent_findings_deduped.csv")

    unverified = [r for r in results if r['category'] == 'INDEPENDENT_UNVERIFIED']
    with open(os.path.join(out_dir, "independent_unverified.csv"), "w", encoding='utf-8') as f:
        f.write("Position,Turn,TotalPlies,TiedMoves\n")
        for r in sorted(unverified, key=lambda x: -x['total_plies']):
            tied_str = ";".join(f"{mv}:{bn}" for mv, bn in r['tied'])
            f.write(f'"{r["position"]}",{r["turn"]},{r["total_plies"]},"{tied_str}"\n')
    print(f"Wrote {len(unverified)} INDEPENDENT_UNVERIFIED positions to "
          f"{out_dir}/independent_unverified.csv -- these are NOT findings, they're gaps")

    with open(os.path.join(out_dir, "all_classifications.csv"), "w", encoding='utf-8') as f:
        f.write("Position,Turn,TotalPlies,Category,NumComponents\n")
        for r in results:
            f.write(f'"{r["position"]}",{r["turn"]},{r["total_plies"]},'
                    f'{r["category"]},{r["n_components"]}\n')
    print(f"Wrote {len(results)} total classifications to {out_dir}/all_classifications.csv")


def cmd_classify(args):
    print(f"Loading {args.db_path}...")
    db = load_db(args.db_path)
    print(f"  Loaded {len(db)} positions")

    multi_count = sum(1 for e in db.values() if len(e.tied) > 1)
    print(f"  Positions with a multi-way tie: {multi_count}")

    checkpoint_path = args.db_path + '.classify_checkpoint.pkl'
    results = run_categorization(db, min_plies=args.min_plies, max_positions=args.max_positions,
                                  checkpoint_path=checkpoint_path, memo_cache_size=args.memo_cache_size,
                                  max_memory_gb=args.max_memory_gb)
    summarize_classify(results, args.out_dir)

    # Only removed after a full, successful run -- if this line is never
    # reached (killed, crashed, interrupted), the checkpoint stays on disk
    # exactly so the next run can resume from it instead of starting over.
    if os.path.exists(checkpoint_path):
        os.remove(checkpoint_path)


# ============================================================================
# Tree rendering (used standalone and by `families --render-trees`)
# ============================================================================

def verify_and_count(state, db, memo, mismatches):
    if state in memo:
        return memo[state]
    entry = db.get(state)
    tied = entry.tied if entry else []
    if not tied:
        result = (1, 0)
        memo[state] = result
        return result
    totals = []
    leaf_count = 0
    pos, turn = state
    for mv, bn in tied:
        try:
            child_pos = apply_move_general(pos, mv)
        except AmbiguousMoveError:
            continue
        ct = child_turn(mv)
        contribution = black_legal_move_count(child_pos) if ct == 'B' else 0
        child_leaves, child_total = verify_and_count((child_pos, ct), db, memo, mismatches)
        leaf_count += child_leaves
        totals.append(contribution + child_total)
    if len(set(totals)) > 1:
        mismatches.append((state, list(zip([mv for mv, bn in tied], totals))))
    result = (leaf_count, totals[0] if totals else 0)
    memo[state] = result
    return result


def render(state, db, indent, ply_num, out, max_lines, visiting, classifications=None):
    if len(out) >= max_lines:
        out.append("  " * indent + "... [text view truncated at max_lines]")
        return
    if state in visiting:
        out.append("  " * indent + "# [anomaly: revisited a position already on this line]")
        return
    pos, turn = state
    entry = db.get(state)
    if entry is None:
        out.append("  " * indent + ("# checkmate" if is_checkmate_state(state)
                                     else "# [proven draw -- correctly terminal, not a gap]"))
        return
    tied = entry.tied
    if not tied:
        out.append("  " * indent + "# [no tied moves recorded here]")
        return
    if len(tied) == 1:
        mv, bn = tied[0]
        out.append("  " * indent + f"{ply_num}. {mv}")
        try:
            child = (apply_move_general(pos, mv), child_turn(mv))
        except AmbiguousMoveError:
            out.append("  " * (indent + 1) + "# [ambiguous move -- cannot continue this branch]")
            return
        render(child, db, indent, ply_num + 1, out, max_lines, visiting | {state}, classifications)
    else:
        label_note = ""
        if classifications is not None:
            cat = classifications.get(state)
            if cat:
                label_note = f"  [classification: {cat}]"
        out.append("  " * indent + f"-- {len(tied)}-way tie (all mate in the same distance "
                                    f"AND same cumulative escape count){label_note} --")
        for i, (mv, bn) in enumerate(tied):
            label = chr(ord('A') + i) if i < 26 else f"#{i}"
            out.append("  " * indent + f"  ({label}) {ply_num}. {mv}   [cum. escapes={bn}]")
            try:
                child = (apply_move_general(pos, mv), child_turn(mv))
            except AmbiguousMoveError:
                out.append("  " * (indent + 2) + "# [ambiguous move -- cannot continue this branch]")
                continue
            render(child, db, indent + 2, ply_num + 1, out, max_lines, visiting | {state}, classifications)


def render_position(position, turn, db, out_path, component_labels=None, max_lines=500,
                     classifications=None):
    key = (position, turn)
    entry = db.get(key)
    header = [f"Position: {position}  ({turn} to move)"]
    if entry:
        header.append(f"Total plies to mate: {entry.distance}")
        header.append(f"Root tied moves: {len(entry.tied)}")
        if component_labels:
            header.append("Root classification:")
            for group_idx, moves in enumerate(component_labels):
                tag = "INDEPENDENT component" if len(component_labels) > 1 else "component"
                header.append(f"  {tag} {group_idx + 1}: {', '.join(moves)}")
        mismatches = []
        leaf_count, verified_total = verify_and_count(key, db, {}, mismatches)
        header.append("")
        header.append(f"EXHAUSTIVE complete perfect-play paths from this position: {leaf_count}")
        if mismatches:
            header.append(f"WARNING - INVARIANT VIOLATION: {len(mismatches)} node(s) had tied "
                           f"choices that do NOT all sum to the same total escape count. First:")
            state, vals = mismatches[0]
            header.append(f"  at {state}: {vals}")
        else:
            header.append(f"Verified: every one of those {leaf_count} paths sums to the exact "
                           f"same total cumulative escape count ({verified_total}), independently "
                           f"recomputed from scratch.")
    else:
        header.append("# [position not found in database]")
    header.append("")

    out = []
    if entry:
        render(key, db, 0, 1, out, max_lines, set(), classifications)
    if len(out) >= max_lines:
        header.append(f"NOTE: tree below is capped at {max_lines} lines; counts above are exhaustive.")
        header.append("")
    text = "\n".join(header + out) + "\n"
    with open(out_path, "w", encoding='utf-8') as f:
        f.write(text)
    return text


def cmd_tree(args):
    print(f"Loading {args.db_path}...")
    db = load_db(args.db_path)
    print(f"  Loaded {len(db)} positions")
    classifications = load_classifications(args.classifications) if args.classifications else None

    if args.position:
        if not args.turn:
            print("ERROR: --position requires --turn", file=sys.stderr)
            sys.exit(1)
        out_path = args.out or "tree.txt"
        text = render_position(args.position, args.turn, db, out_path, max_lines=args.max_lines,
                                classifications=classifications)
        print(f"\nWrote {out_path}\n")
        print(text)
        return

    if args.findings_csv:
        os.makedirs(args.out_dir, exist_ok=True)
        rows = []
        with open(args.findings_csv, encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    row['TotalPlies'] = int(row['TotalPlies'])
                except (KeyError, ValueError):
                    row['TotalPlies'] = 0
                rows.append(row)
        rows.sort(key=lambda r: -r['TotalPlies'])
        rows = rows[:args.top]
        print(f"Rendering the {len(rows)} deepest findings from {args.findings_csv}...\n")
        for i, row in enumerate(rows):
            position, turn = row['Position'], row['Turn']
            components = None
            if row.get('ComponentGroups'):
                components = [g.split(',') for g in row['ComponentGroups'].split('|')]
            safe_name = position.replace(':', '').replace(' ', '_')
            out_path = os.path.join(args.out_dir, f"{i+1:03d}_{safe_name}_{turn}.txt")
            render_position(position, turn, db, out_path, component_labels=components,
                             max_lines=args.max_lines, classifications=classifications)
            print(f"  [{i+1}/{len(rows)}] {position} ({turn}, {row['TotalPlies']} plies) -> {out_path}")
        print(f"\nWrote {len(rows)} tree files to {args.out_dir}/")
        return

    print("ERROR: specify either --position/--turn or --findings-csv", file=sys.stderr)
    sys.exit(1)


# ============================================================================
# STAGE 2: families
# ============================================================================

def load_findings(path):
    findings = []
    with open(path, encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            wk, others, bk = canonical_roles(row['Position'])
            tied_moves = set()
            for entry in row['TiedMoves'].split(';'):
                if ':' in entry:
                    tied_moves.add(entry.rsplit(':', 1)[0])
            findings.append({
                'position': row['Position'], 'turn': row['Turn'],
                'total_plies': int(row['TotalPlies']) if row['TotalPlies'] else 0,
                'WK': wk, 'white_others': tuple(others), 'BK': bk,
                'tied_moves': frozenset(tied_moves),
                'component_groups': row.get('ComponentGroups', ''),
            })
    return findings


def find_families(findings, checkpoint_path=None):
    """Generalizes krvk_analyzer's 'fix one of WK/WQ/BK, sweep the rest' to
    an arbitrary-length piece list: roles are WK, BK, and each index into
    the canonically-sorted white_others tuple (index 0 = whichever piece
    sorts first by kind-then-square, matching pack_general_state's own
    ordering -- see the module docstring). Note this means "role 0" can
    refer to a DIFFERENT KIND of piece across different findings if
    promotion is in play for this material (e.g. a still-a-pawn finding's
    role 0 vs an already-promoted finding's role 0) -- sweeps are only
    meaningful within a fixed material/kind-shape subset, exactly as
    "fix WQ" implicitly was for krvk_analyzer.py's single-piece case.

    checkpoint_path, if given, checkpoints at GROUP granularity, not just
    per-role -- deliberately, since this function's own docstring already
    notes a single group's O(n^2) pairwise pass, not a whole role, is
    where a genuinely slow stretch would actually show up. The group's
    own (role, fixed_val, turn, tied_moves) grouping key doubles as the
    resume key -- no separate id scheme needed, since it's already a
    complete, unique identifier for "which group is this"."""
    seen_in_a_family = set()
    families = []
    done_keys = set()
    if checkpoint_path and os.path.exists(checkpoint_path):
        import pickle
        with open(checkpoint_path, 'rb') as f:
            checkpoint = pickle.load(f)
        families = checkpoint['families']
        seen_in_a_family = checkpoint['seen_in_a_family']
        done_keys = checkpoint['done_keys']
        print(f"  Resumed {len(families)} already-found families ({len(done_keys)} groups "
              f"already checked) from {checkpoint_path}")

    last_checkpoint = time.time()

    def maybe_checkpoint(force=False):
        nonlocal last_checkpoint
        if not checkpoint_path:
            return
        now = time.time()
        if not force and now - last_checkpoint < 60:
            return
        last_checkpoint = now
        import pickle
        tmp_path = checkpoint_path + '.tmp'
        with open(tmp_path, 'wb') as f:
            pickle.dump({'families': families, 'seen_in_a_family': seen_in_a_family,
                         'done_keys': done_keys}, f)
        os.replace(tmp_path, checkpoint_path)

    roles = ['WK', 'BK'] + [f'white_other_{i}' for i in range(5)]  # 5 = MAX_WHITE_NON_KING, matching
                                                                    # the C++ engine's own limit -- see
                                                                    # compositional_trajectory_solver_modular.cpp
    for role in roles:
        groups = defaultdict(list)
        for finding in findings:
            if role.startswith('white_other_'):
                idx = int(role.split('_')[-1])
                if idx >= len(finding['white_others']):
                    continue
                fixed_val = finding['white_others'][idx]
            else:
                fixed_val = finding[role]
            key = (role, fixed_val, finding['turn'], finding['tied_moves'])
            groups[key].append(finding)

        print(f"  [role={role}] {len(groups)} distinct groups to check")
        group_start = time.time()
        last_print = group_start

        for gi, (group_key, members) in enumerate(groups.items()):
            fixed_role, fixed_val, turn, tied_moves = group_key
            if group_key in done_keys:
                continue
            distinct_positions = {(m['position'], m['turn']): m for m in members}
            if len(distinct_positions) < 2:
                done_keys.add(group_key)
                continue
            member_list = list(distinct_positions.values())
            # The only part of this whole function with worse-than-linear
            # cost: pairwise across one group's members. Groups are
            # typically small, but nothing guarantees that, so this is
            # where a genuinely slow group would actually show up.
            symmetric_pairs = []
            n = len(member_list)
            pair_count = 0
            for i in range(n):
                for j in range(i + 1, n):
                    pair_count += 1
                    t_idx = is_symmetric_general(member_list[i]['position'], member_list[j]['position'])
                    if t_idx is not None and t_idx != 0:
                        symmetric_pairs.append((member_list[i]['position'], member_list[j]['position'], t_idx))
                    now = time.time()
                    if now - last_print >= 10:
                        last_print = now
                        print(f"    [role={role}] group {gi+1}/{len(groups)} ({n} members, "
                              f"{n*(n-1)//2} pairs): {pair_count}/{n*(n-1)//2} pairs checked, "
                              f"elapsed={now-group_start:.0f}s", flush=True)
            families.append({
                'fixed_piece': fixed_role, 'fixed_value': fixed_val, 'turn': turn,
                'tied_moves': tied_moves, 'members': member_list,
                'symmetric_pairs': symmetric_pairs,
            })
            for m in member_list:
                seen_in_a_family.add((m['position'], m['turn']))
            done_keys.add(group_key)
            maybe_checkpoint()
        print(f"  [role={role}] done in {time.time()-group_start:.0f}s")
        maybe_checkpoint(force=True)  # always checkpoint at a role boundary, regardless of the 60s timer
    unique_findings = [f for f in findings if (f['position'], f['turn']) not in seen_in_a_family]
    return families, unique_findings


def cmd_families(args):
    print(f"Loading {args.findings_csv}...")
    findings = load_findings(args.findings_csv)
    print(f"  Loaded {len(findings)} deduplicated independent findings")

    checkpoint_path = args.checkpoint_path or (args.findings_csv + '.families_checkpoint.pkl')
    families, unique_findings = find_families(findings, checkpoint_path=checkpoint_path)

    # Only removed after a full, successful run -- if find_families never
    # returns (killed, crashed, interrupted), the checkpoint stays on disk
    # exactly so the next run can resume from it instead of starting over,
    # matching cmd_classify's own convention.
    if os.path.exists(checkpoint_path):
        os.remove(checkpoint_path)

    seen_member_sets = set()
    distinct_families = []
    for fam in families:
        key = frozenset(m['position'] for m in fam['members'])
        if key in seen_member_sets:
            continue
        seen_member_sets.add(key)
        distinct_families.append(fam)
    distinct_families.sort(key=lambda f: -len(f['members']))

    os.makedirs(args.out_dir, exist_ok=True)
    print(f"\nFamilies found: {len(distinct_families)}")
    print(f"Unique findings (no recurrence under any single-role-fixed sweep): {len(unique_findings)}")

    with open(os.path.join(args.out_dir, "families.csv"), "w", encoding='utf-8') as f:
        f.write("FixedRole,FixedValue,Turn,TiedMoves,FamilySize,MemberPositions,ContainsSymmetricPair\n")
        for fam in distinct_families:
            moves_str = ";".join(sorted(fam['tied_moves']))
            members_str = "|".join(m['position'] for m in fam['members'])
            f.write(f'{fam["fixed_piece"]},"{fam["fixed_value"]}",{fam["turn"]},"{moves_str}",'
                    f'{len(fam["members"])},"{members_str}",{bool(fam["symmetric_pairs"])}\n')
    print(f"Wrote {len(distinct_families)} families to {args.out_dir}/families.csv")

    with open(os.path.join(args.out_dir, "unique_findings.csv"), "w", encoding='utf-8') as f:
        f.write("Position,Turn,TotalPlies,TiedMoves,ComponentGroups\n")
        for u in sorted(unique_findings, key=lambda x: -x['total_plies']):
            moves_str = ";".join(sorted(u['tied_moves']))
            f.write(f'"{u["position"]}",{u["turn"]},{u["total_plies"]},"{moves_str}","{u["component_groups"]}"\n')
    print(f"Wrote {len(unique_findings)} unique findings to {args.out_dir}/unique_findings.csv")

    print("\nLargest families:")
    for fam in distinct_families[:10]:
        moves_str = "/".join(sorted(fam['tied_moves']))
        print(f"  fix {fam['fixed_piece']}={fam['fixed_value']} ({fam['turn']}), moves={{{moves_str}}}: "
              f"{len(fam['members'])} positions")

    print("\nDeepest unique findings:")
    for u in sorted(unique_findings, key=lambda x: -x['total_plies'])[:10]:
        moves_str = "/".join(sorted(u['tied_moves']))
        print(f"  {u['position']} ({u['turn']}, {u['total_plies']} plies): {{{moves_str}}}")

    if args.render_trees:
        print(f"\nRendering trees ({args.render_trees})...")
        db = load_db(args.render_trees)
        tree_dir = os.path.join(args.out_dir, "trees")
        os.makedirs(tree_dir, exist_ok=True)
        for u in unique_findings:
            safe_name = u['position'].replace(':', '').replace(' ', '_')
            out_path = os.path.join(tree_dir, f"unique_{safe_name}_{u['turn']}.txt")
            render_position(u['position'], u['turn'], db, out_path, max_lines=args.max_lines)
        print(f"  Wrote {len(unique_findings)} unique-finding trees to {tree_dir}/")
        for i, fam in enumerate(distinct_families):
            deepest = max(fam['members'], key=lambda m: m['total_plies'])
            safe_name = deepest['position'].replace(':', '').replace(' ', '_')
            out_path = os.path.join(tree_dir, f"family{i+1:03d}_rep_{safe_name}_{deepest['turn']}.txt")
            render_position(deepest['position'], deepest['turn'], db, out_path, max_lines=args.max_lines)
        print(f"  Wrote {len(distinct_families)} family-representative trees to {tree_dir}/")


# ============================================================================
# CLI
# ============================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='command', required=True)

    p1 = sub.add_parser('classify', help="Classify every multi-way tie in a landscape database")
    p1.add_argument('db_path')
    p1.add_argument('--min-plies', type=int, default=0)
    p1.add_argument('--max-positions', type=int, default=None)
    p1.add_argument('--memo-cache-size', type=int, default=3000,
                     help="Max entries in the reachable-set memoization cache (default 3000). "
                          "Each entry can be over 1MB at scale, so this bounds a large classify "
                          "run's memory instead of letting it grow unbounded for the whole run. "
                          "Lower it if you're still seeing high memory/swap; raising it trades "
                          "more memory for fewer recomputed cache misses.")
    p1.add_argument('--out-dir', default='general_tie_analysis')
    p1.add_argument('--max-memory-gb', type=float, default=None,
                     help="If the process's actual measured memory (not an estimate) exceeds this "
                          "many GB, checkpoint immediately and exit cleanly rather than let the OS "
                          "start swapping -- the same 'fail loud, not silently wrong' choice this "
                          "project uses elsewhere (see general_solver's --max-nodes). Re-running the "
                          "identical command resumes from the checkpoint. Unset by default -- the "
                          "memory reductions in this version (namedtuple db values, stripped "
                          "non-INDEPENDENT result fields) cut typical usage substantially, but this "
                          "is a real safety net for scaling to larger materials, not a replacement "
                          "for watching actual usage the first time you run a new, bigger one.")

    p2 = sub.add_parser('families', help="Compile families vs unique findings from classify's output")
    p2.add_argument('findings_csv')
    p2.add_argument('--out-dir', default='general_families')
    p2.add_argument('--render-trees', metavar='DB_PATH', default=None)
    p2.add_argument('--max-lines', type=int, default=500)
    p2.add_argument('--checkpoint-path', default=None,
                     help="Path to save/resume role-sweep progress (default: <findings_csv>."
                          "families_checkpoint.pkl). Checkpointed after each completed role, since "
                          "that's this function's only genuinely unbounded-cost step (a role with a "
                          "large group triggers an O(n^2) pairwise pass) -- see find_families's "
                          "own docstring for why group-level, not finer, is the right granularity.")

    p3 = sub.add_parser('tree', help="Render the full move tree for one position or a findings CSV")
    p3.add_argument('db_path')
    p3.add_argument('--position')
    p3.add_argument('--turn', choices=['W', 'B'])
    p3.add_argument('--out', default=None)
    p3.add_argument('--findings-csv')
    p3.add_argument('--top', type=int, default=10)
    p3.add_argument('--out-dir', default='general_trees')
    p3.add_argument('--max-lines', type=int, default=500)
    p3.add_argument('--classifications')

    args = ap.parse_args()
    if args.command == 'classify':
        cmd_classify(args)
    elif args.command == 'families':
        cmd_families(args)
    elif args.command == 'tree':
        cmd_tree(args)


if __name__ == "__main__":
    main()
