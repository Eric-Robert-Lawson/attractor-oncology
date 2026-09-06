#!/usr/bin/env python3
"""
KQvK Perfect-Play Tie Categorizer
====================================

For every position in your landscape database with a genuine multi-way tie
(multiple moves achieving BOTH optimal mate-distance AND optimal cumulative
black-escape count -- i.e. every entry in TiedMoves), this classifies WHY
the tie exists, so you aren't sorting real findings out of thousands of
trivial duplicates by hand.

IMPORTANT -- --batch vs --full-dag data completeness:
This script's transposition checks work by walking forward from a tied
move's child, following THAT position's own TiedMoves, and so on -- which
only works if the database actually contains those descendant positions.
`--full-dag` mode explores every tied branch, so it does. `--batch` mode
walks exactly ONE canonical path per root and never explores the sibling
branches of a tie at all, even though TiedMoves is still correctly recorded
at every position it does visit (that computation happens inside
find_best_move regardless of which driver calls it). That means in a
--batch-generated file, most non-canonical tied-move children simply don't
exist in the database -- not because they're terminal, but because they
were never explored. Silently treating "child not found" the same as "child
is a real dead end" would manufacture false INDEPENDENT verdicts on a
massive scale for --batch data: "we found no reconvergence" and "we have no
data to check" are NOT the same thing, and conflating them was the same
class of error as the FIFO-queue bug from earlier in this project, just
structural this time instead of incidental.

This version fixes that by independently verifying, using its own from-
scratch checkmate check (never reading TiedMoves to make this determination
-- see is_checkmate_state), whether a missing position is a genuine mate
(a real, informative dead end) or an unexplored gap (uninformative). Every
branch's reachability walk now reports whether it hit a genuine gap:

  - INDEPENDENT            -- 2+ connected components, AND every branch
                                involved was walked to completion (only real
                                mates were ever "missing", never a gap).
                                This is a real, load-bearing finding.
  - INDEPENDENT_UNVERIFIED -- 2+ connected components, but at least one
                                branch's walk hit an unexplored gap. The
                                apparent independence might be real, or it
                                might just mean the missing branch was never
                                explored far enough to find a connection.
                                Treat this as "unknown", not as a finding.

On a --full-dag file that finished exploring, expect INDEPENDENT_UNVERIFIED
to be rare (only true where the depth budget or node cap was hit) or your
own results to hold up. On a --batch file, expect it to dominate --batch was
never designed to produce the data this analysis actually needs.

Every tied move is a node in a small graph; an edge is drawn between two
tied moves whenever they're related by ANY of:

  (a) DIRECT SYMMETRY   -- some D4 board symmetry maps one move's resulting
                            child position exactly onto the other's. Doesn't
                            depend on the database at all, so this edge is
                            always trustworthy regardless of data gaps.
  (b) TRANSPOSITION      -- the two branches' full descendant trees (every
                            position reachable by repeatedly following
                            TiedMoves) share an identical future position.
                            A FOUND connection is always trustworthy evidence
                            -- it's the absence of one that becomes ambiguous
                            when a gap was hit.
  (c) SYMMETRIC           -- one branch's ENTIRE descendant tree, reflected
      TRANSPOSITION         or rotated by some single D4 element, shares a
                            position with the other branch's descendant
                            tree. Catches a mirror relationship that only
                            appears a few plies deeper.

Connected components of this graph are computed via union-find, which is
what correctly handles 3+-way ties: not every pair needs to be directly
related for the whole set to collapse into one explained group, as long as
a CHAIN of relations connects them.

Positions where the parent is already mate-in-1 are tagged TRIVIAL_MATE and
excluded from the above (multiple mating moves from an already-won position
carry no information -- there's nothing left to compare).

Every INDEPENDENT (verified) finding is canonicalized under the board's full
D4 symmetry group and deduplicated, because the database's own automatic
8-fold symmetry expansion means a single genuine discovery can appear as up
to 8 near-identical rows. Statistics are broken down BY PLY-DEPTH
specifically, so you can see whether independence emerges more, less, or
not at all as mate-distance grows, instead of one aggregate number that
buries the trend.

Usage:
    python3 tie_analysis.py full_dag.db --out-dir tie_analysis
    python3 tie_analysis.py kqvk_perfect_play.db --out-dir tie_analysis_batch
    python3 tie_analysis.py full_dag.db --min-plies 8 --out-dir tie_analysis_deep
"""

import argparse
import sys
from collections import defaultdict

FILES = "abcdefgh"


# ============================================================================
# Board geometry / D4 symmetry group (same 8 elements as the C++ engine's
# SolvedPositionDatabase::add_position and the earlier exact-landscape tool)
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
    lambda f, r: (7 - f, r),
    lambda f, r: (r, f),
    lambda f, r: (7 - r, 7 - f),
]


def parse_position(pos_str):
    parts = {}
    for tok in pos_str.split():
        k, v = tok.split(':')
        parts[k] = sq_to_coord(v)
    return parts


def position_triple(pos_str):
    p = parse_position(pos_str)
    return (p['WK'], p['WQ'], p['BK'])


def canonical_form(triple):
    """Smallest image of (wk,wq,bk) under the 8-element D4 group -- used to
    deduplicate findings that are the same underlying discovery replicated
    by the database's own symmetry expansion."""
    best = None
    for T in TRANSFORMS:
        img = tuple(T(f, r) for f, r in triple)
        if best is None or img < best:
            best = img
    return best


def apply_move(pos_str, mv):
    parts = {}
    for tok in pos_str.split():
        k, v = tok.split(':')
        parts[k] = v
    piece, dest = mv[0], mv[1:]
    if piece == 'K':
        parts['WK'] = dest
    elif piece == 'Q':
        parts['WQ'] = dest
    elif piece == 'k':
        parts['BK'] = dest
    return f"WK:{parts['WK']} WQ:{parts['WQ']} BK:{parts['BK']}"


def child_turn(mv):
    return 'B' if mv[0] in ('K', 'Q') else 'W'


def transform_state(state, T):
    pos, turn = state
    triple = position_triple(pos)
    tt = tuple(T(f, r) for f, r in triple)
    new_pos = f"WK:{coord_to_sq(*tt[0])} WQ:{coord_to_sq(*tt[1])} BK:{coord_to_sq(*tt[2])}"
    return (new_pos, turn)


# ============================================================================
# Independent, from-scratch checkmate check (mirrors BaseEngine::is_checkmate
# in the C++ engine exactly). Used ONLY to tell a genuine dead end (mate)
# apart from an unexplored database gap -- never reads TiedMoves or anything
# else derived from the search, so it can't be fooled by the same gap it's
# meant to detect.
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


def is_attacked_by_queen(pos, qp, wk, wq):
    if pos == qp:
        return False
    pf, pr = pos
    qf, qr = qp
    if pf == qf:
        lo, hi = sorted([pr, qr])
        for r in range(lo + 1, hi):
            if (pf, r) == wk or (pf, r) == wq:
                return False
        return True
    if pr == qr:
        lo, hi = sorted([pf, qf])
        for f in range(lo + 1, hi):
            if (f, pr) == wk or (f, pr) == wq:
                return False
        return True
    if abs(pf - qf) == abs(pr - qr):
        df = 1 if pf > qf else -1
        dr = 1 if pr > qr else -1
        f, r = qf + df, qr + dr
        while f != pf:
            if (f, r) == wk or (f, r) == wq:
                return False
            f += df
            r += dr
        return True
    return False


def is_checkmate_state(state):
    """state = (position_str, turn). Mate is only ever a Black-to-move fact
    in this material, matching BaseEngine::is_checkmate's own precondition."""
    pos_str, turn = state
    if turn != 'B':
        return False
    try:
        p = parse_position(pos_str)
        wk, wq, bk = p['WK'], p['WQ'], p['BK']
    except (KeyError, ValueError, IndexError):
        return False
    if not is_attacked_by_queen(bk, wq, wk, wq):
        return False
    for m in king_moves(bk):
        if is_attacked_by_queen(m, wq, wk, wq):
            continue
        if max(abs(m[0] - wk[0]), abs(m[1] - wk[1])) <= 1:
            continue
        return False
    return True


# ============================================================================
# Loading
# ============================================================================

def load_db(path):
    """Returns dict (position, turn) -> {'total_plies': int, 'tied': [(mv,bn),...]}"""
    db = {}
    with open(path) as f:
        header = f.readline()
        for line_num, line in enumerate(f, start=2):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("|")
            if len(parts) < 5:
                continue
            position, turn = parts[0], parts[1]
            try:
                total_plies = int(parts[4])
            except ValueError:
                continue
            tied = []
            if len(parts) > 11 and parts[11]:
                for entry in parts[11].split(';'):
                    if ':' not in entry:
                        continue
                    mv, bn_str = entry.rsplit(':', 1)
                    try:
                        tied.append((mv, int(bn_str)))
                    except ValueError:
                        pass
            db[(position, turn)] = {'total_plies': total_plies, 'tied': tied}
    return db


# ============================================================================
# Union-Find
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


# ============================================================================
# Reachability (memoized globally -- shared subtrees are only ever walked
# once across the whole analysis, not once per tie that touches them)
# ============================================================================

def reachable_set(state, db, memo, max_nodes=5000, warned=None):
    """Returns (frozenset of reached states, truncated: bool). truncated is
    True iff the walk ever hit a position that is NOT in the database and is
    NOT independently verified as a genuine checkmate -- i.e. a real data
    gap, not a real dead end. Also True if the max_nodes safety cap is hit,
    since that's just as much "we don't actually know what's further out"."""
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
            if not is_checkmate_state(s):
                truncated = True  # genuine gap: not in the DB, and not a real mate either
            continue
        for mv, bn in entry['tied']:
            pos, turn = s
            child = (apply_move(pos, mv), child_turn(mv))
            if child not in seen:
                stack.append(child)
    result = (frozenset(seen), truncated)
    memo[state] = result
    return result


# ============================================================================
# Classification
# ============================================================================

def classify_position(position, turn, entry, db, reach_memo, warned):
    tied = entry['tied']
    n = len(tied)
    if n <= 1:
        return None  # nothing to classify
    if entry['total_plies'] == 1:
        return {'category': 'TRIVIAL_MATE', 'components': None, 'n_components': 1}

    children = [(apply_move(position, mv), child_turn(mv)) for mv, bn in tied]
    triples = [position_triple(c[0]) for c in children]

    uf = UnionFind(n)

    # (a) direct symmetry between children -- doesn't touch the database at
    # all, so this is trustworthy regardless of any data gaps.
    for i in range(n):
        for j in range(i + 1, n):
            if children[i][1] != children[j][1]:
                continue  # different side to move can't be a board symmetry of each other
            for T in TRANSFORMS[1:]:
                img = tuple(T(f, r) for f, r in triples[i])
                if img == triples[j]:
                    uf.union(i, j)
                    break

    # Precompute reachable sets once per branch (memoized globally in reach_memo)
    branch_results = [reachable_set(c, db, reach_memo, warned=warned) for c in children]
    branch_sets = [r[0] for r in branch_results]
    branch_truncated = [r[1] for r in branch_results]
    any_truncated = any(branch_truncated)

    # (b) direct transposition -- a FOUND connection is always real evidence,
    # truncation elsewhere doesn't undermine an actual match.
    for i in range(n):
        for j in range(i + 1, n):
            if uf.find(i) == uf.find(j):
                continue
            if branch_sets[i] & branch_sets[j]:
                uf.union(i, j)

    # (c) symmetric transposition -- reflect/rotate branch i's WHOLE reachable
    # set by each non-identity transform and check for overlap with branch j's.
    # Transformed sets are cached per (branch root, transform index) since 3+-way
    # ties re-check the same branch against multiple siblings.
    transformed_cache = {}

    def get_transformed(idx, t_idx):
        key = (children[idx], t_idx)
        if key not in transformed_cache:
            T = TRANSFORMS[t_idx]
            transformed_cache[key] = {transform_state(s, T) for s in branch_sets[idx]}
        return transformed_cache[key]

    for i in range(n):
        for j in range(n):
            if i == j or uf.find(i) == uf.find(j):
                continue
            for t_idx in range(1, len(TRANSFORMS)):
                if get_transformed(i, t_idx) & branch_sets[j]:
                    uf.union(i, j)
                    break

    components = uf.components()
    if len(components) == 1:
        category = 'EXPLAINED'
    elif any_truncated:
        # 2+ components remain, but at least one branch was never fully
        # walked -- the missing part of the tree could still hide a
        # connection. This is NOT a finding, it's "unknown".
        category = 'INDEPENDENT_UNVERIFIED'
    else:
        category = 'INDEPENDENT'
    return {'category': category, 'components': components, 'n_components': len(components),
            'tied': tied, 'children': children, 'any_truncated': any_truncated}


def run_categorization(db, min_plies=0, max_positions=None):
    reach_memo = {}
    warned = set()
    results = []
    checked = 0

    items = list(db.items())
    for (position, turn), entry in items:
        if entry['total_plies'] < min_plies:
            continue
        if len(entry['tied']) <= 1:
            continue
        checked += 1
        if max_positions and checked > max_positions:
            break
        result = classify_position(position, turn, entry, db, reach_memo, warned)
        if result:
            result['position'] = position
            result['turn'] = turn
            result['total_plies'] = entry['total_plies']
            results.append(result)

    return results


def summarize(results, out_dir):
    import os
    os.makedirs(out_dir, exist_ok=True)

    by_plies = defaultdict(lambda: defaultdict(int))
    for r in results:
        by_plies[r['total_plies']][r['category']] += 1

    print("\nBy ply-depth (total_plies), raw row counts:")
    print(f"{'plies':>6} {'TRIVIAL':>9} {'EXPLAINED':>10} {'INDEPENDENT':>12} {'UNVERIFIED':>11}")
    for plies in sorted(by_plies):
        row = by_plies[plies]
        print(f"{plies:>6} {row.get('TRIVIAL_MATE',0):>9} {row.get('EXPLAINED',0):>10} "
              f"{row.get('INDEPENDENT',0):>12} {row.get('INDEPENDENT_UNVERIFIED',0):>11}")

    total_classified = len(results)
    n_unverified = sum(1 for r in results if r['category'] == 'INDEPENDENT_UNVERIFIED')
    n_independent = sum(1 for r in results if r['category'] == 'INDEPENDENT')
    if total_classified:
        unverified_frac = n_unverified / total_classified
        print(f"\n{n_unverified}/{total_classified} ({unverified_frac:.1%}) of ALL classified positions "
              f"are INDEPENDENT_UNVERIFIED -- apparent independence that could not be distinguished "
              f"from a missing branch. This will be large for --batch-sourced data (expected -- see "
              f"the module docstring) and should be small for a fully-explored --full-dag file.")
        if unverified_frac > 0.5:
            print("  \u26a0 Over half of all classified positions are unverified. If this file came "
                  "from --batch mode, that's exactly the expected outcome, not a bug: --batch never "
                  "explores sibling branches, so most ties simply can't be checked. Treat every "
                  "INDEPENDENT count below as a LOWER BOUND at best, not a real measurement, until "
                  "this is run against a --full-dag file instead.")

    independent = [r for r in results if r['category'] == 'INDEPENDENT']
    print(f"\nTotal INDEPENDENT rows (raw, verified, before dedup): {len(independent)}")

    # Deduplicate by D4 orbit of the PARENT position (the row this tie lives at)
    seen_canon = {}
    deduped = []
    for r in independent:
        triple = position_triple(r['position'])
        canon = canonical_form(triple)
        key = (canon, r['turn'], r['total_plies'])
        if key not in seen_canon:
            seen_canon[key] = r
            deduped.append(r)

    print(f"Distinct INDEPENDENT (verified) findings after D4-orbit dedup: {len(deduped)}")

    by_plies_dedup = defaultdict(int)
    for r in deduped:
        by_plies_dedup[r['total_plies']] += 1
    print("\nDeduplicated INDEPENDENT (verified) findings by ply-depth:")
    for plies in sorted(by_plies_dedup):
        print(f"  {plies:>3} plies: {by_plies_dedup[plies]}")

    with open(f"{out_dir}/independent_findings_deduped.csv", "w") as f:
        f.write("Position,Turn,TotalPlies,NumComponents,TiedMoves,ComponentGroups\n")
        for r in sorted(deduped, key=lambda x: -x['total_plies']):
            tied_str = ";".join(f"{mv}:{bn}" for mv, bn in r['tied'])
            groups_str = "|".join(
                ",".join(f"{r['tied'][i][0]}" for i in comp) for comp in r['components']
            )
            f.write(f'"{r["position"]}",{r["turn"]},{r["total_plies"]},'
                    f'{r["n_components"]},"{tied_str}","{groups_str}"\n')
    print(f"\n\u2713 Wrote {len(deduped)} deduplicated VERIFIED independent findings to "
          f"{out_dir}/independent_findings_deduped.csv")

    unverified = [r for r in results if r['category'] == 'INDEPENDENT_UNVERIFIED']
    with open(f"{out_dir}/independent_unverified.csv", "w") as f:
        f.write("Position,Turn,TotalPlies,NumComponents,TiedMoves\n")
        for r in sorted(unverified, key=lambda x: -x['total_plies']):
            tied_str = ";".join(f"{mv}:{bn}" for mv, bn in r['tied'])
            f.write(f'"{r["position"]}",{r["turn"]},{r["total_plies"]},'
                    f'{r["n_components"]},"{tied_str}"\n')
    print(f"\u2713 Wrote {len(unverified)} INDEPENDENT_UNVERIFIED (data-limited) positions to "
          f"{out_dir}/independent_unverified.csv -- these are NOT findings, they're gaps")

    with open(f"{out_dir}/all_classifications.csv", "w") as f:
        f.write("Position,Turn,TotalPlies,Category,NumComponents\n")
        for r in results:
            f.write(f'"{r["position"]}",{r["turn"]},{r["total_plies"]},'
                    f'{r["category"]},{r["n_components"]}\n')
    print(f"\u2713 Wrote {len(results)} total classifications to {out_dir}/all_classifications.csv")

    print("\nSample deduplicated VERIFIED independent findings (deepest first):")
    for r in sorted(deduped, key=lambda x: -x['total_plies'])[:15]:
        groups = [[r['tied'][i][0] for i in comp] for comp in r['components']]
        group_str = "  vs  ".join("/".join(g) for g in groups)
        print(f"  {r['position']} ({r['turn']} to move, {r['total_plies']} plies): {group_str}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("db_path")
    ap.add_argument("--min-plies", type=int, default=0,
                     help="Only classify positions with at least this many plies remaining "
                          "(use this to skip the shallow, mate-in-1-dominated end and focus "
                          "on deeper, more complex positions)")
    ap.add_argument("--max-positions", type=int, default=None,
                     help="Cap the number of multi-way-tied positions classified, for a quick look")
    ap.add_argument("--out-dir", default="tie_analysis")
    args = ap.parse_args()

    print(f"Loading {args.db_path}...")
    db = load_db(args.db_path)
    print(f"  Loaded {len(db)} positions")

    multi_count = sum(1 for e in db.values() if len(e['tied']) > 1)
    print(f"  Positions with a multi-way tie (any ply-depth): {multi_count}")
    if args.min_plies:
        multi_count_filtered = sum(1 for e in db.values()
                                    if len(e['tied']) > 1 and e['total_plies'] >= args.min_plies)
        print(f"  ...of which >= {args.min_plies} plies: {multi_count_filtered}")

    # Quick, cheap heuristic check up front: how many multi-way-tied positions
    # have AT LEAST ONE tied move whose immediate child isn't even in the
    # database? This doesn't require the full reachability walk and gives an
    # early signal of whether this looks like --batch-style data before the
    # (much slower) full categorization runs.
    missing_immediate = 0
    multi_positions = [(k, e) for k, e in db.items() if len(e['tied']) > 1]
    for (position, turn), e in multi_positions:
        for mv, bn in e['tied']:
            child = (apply_move(position, mv), child_turn(mv))
            if child not in db and not is_checkmate_state(child):
                missing_immediate += 1
                break
    if multi_positions:
        frac = missing_immediate / len(multi_positions)
        print(f"  Quick check: {missing_immediate}/{len(multi_positions)} ({frac:.1%}) of multi-way-tied "
              f"positions have at least one tied move whose child isn't in this file at all.")
        if frac > 0.3:
            print("  \u26a0 This looks like it may be --batch-sourced (or a partial --full-dag run) -- "
                  "expect INDEPENDENT_UNVERIFIED to dominate the results below.")

    results = run_categorization(db, min_plies=args.min_plies, max_positions=args.max_positions)
    summarize(results, args.out_dir)


if __name__ == "__main__":
    main()
