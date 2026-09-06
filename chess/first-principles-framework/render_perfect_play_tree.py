#!/usr/bin/env python3
"""
KQvK Perfect-Play Tree Renderer
==================================

Given a position (typically one flagged by tie_analysis.py's
independent_findings_deduped.csv), walks the FULL perfect-play DAG rooted
there -- following every entry in TiedMoves recursively, branching wherever
a genuine tie exists -- and renders it as readable, indented game notation
all the way to mate.

Forced (single-choice) stretches are printed as one flowing line, exactly
like a normal PGN main line. A branch point is only ever printed where a
REAL choice exists (len(TiedMoves) > 1), with each alternative labeled
(A), (B), (C)... and indented under it. This mirrors how a human annotator
writes up a line with sub-variations, rather than dumping one row per ply.

The root position's OWN top-level branches are additionally labeled with
their classification component (from tie_analysis.py's ComponentGroups
column, if you point --findings-csv at it) -- so you can see directly which
branches are "the same underlying idea" (explained) and which are the
independently-tied ones that were actually the point of looking at this
position at all.

Usage:
    # Render one specific position
    python3 render_perfect_play_tree.py full_dag.db \\
        --position "WK:a4 WQ:e4 BK:h4" --turn W --out tree_a4e4h4.txt

    # Render the N deepest verified independent findings from a prior
    # tie_analysis.py run, one file per finding
    python3 render_perfect_play_tree.py full_dag.db \\
        --findings-csv tie_analysis/independent_findings_deduped.csv \\
        --top 10 --out-dir trees/
"""

import argparse
import csv
import os
import sys

FILES = "abcdefgh"


# ============================================================================
# Shared utilities (same definitions as tie_analysis.py, kept self-contained
# here so this script doesn't depend on tie_analysis.py's file location)
# ============================================================================

def sq_to_coord(s):
    return (FILES.index(s[0]), int(s[1:]) - 1)


def parse_position(pos_str):
    parts = {}
    for tok in pos_str.split():
        k, v = tok.split(':')
        parts[k] = sq_to_coord(v)
    return parts


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


def black_legal_move_count(wk, wq, bk):
    """Independently recomputes Black's TRUE legal move count from scratch --
    never reads own_contribution/bn_cum from the database -- so it can serve
    as a real check on the data, not just a repeat of what's already stored."""
    cnt = 0
    for m in king_moves(bk):
        if is_attacked_by_queen(m, wq, wk, wq):
            continue
        if max(abs(m[0] - wk[0]), abs(m[1] - wk[1])) <= 1:
            continue
        if m == wk or m == wq:
            continue
        cnt += 1
    return cnt


def verify_and_count(state, db, memo, mismatches):
    """Returns (leaf_count, verified_total_bn_cum) for `state`, computed by
    MEMOIZED recursion over distinct nodes -- not brute-force enumeration of
    every path -- so this stays cheap (O(nodes)) even when the true number of
    leaf paths is astronomically larger than the node count. Independently
    recomputes Black's legal move count at every Black-to-move node rather
    than trusting the stored bn_cum; if the "every tied choice sums to the
    same total" invariant is ever violated by the actual data, `mismatches`
    records exactly where, instead of silently producing a wrong number."""
    if state in memo:
        return memo[state]
    entry = db.get(state)
    tied = entry['tied'] if entry else []
    if not tied:
        result = (1, 0)  # one complete path ends here; nothing further to add
        memo[state] = result
        return result

    totals = []
    leaf_count = 0
    pos, turn = state
    for mv, bn in tied:
        child_pos = apply_move(pos, mv)
        ct = child_turn(mv)
        contribution = 0
        if ct == 'B':
            p = parse_position(child_pos)
            contribution = black_legal_move_count(p['WK'], p['WQ'], p['BK'])
        child_leaves, child_total = verify_and_count((child_pos, ct), db, memo, mismatches)
        leaf_count += child_leaves
        totals.append(contribution + child_total)

    if len(set(totals)) > 1:
        mismatches.append((state, list(zip([mv for mv, bn in tied], totals))))

    result = (leaf_count, totals[0] if totals else 0)
    memo[state] = result
    return result


def load_classifications(path):
    """Loads all_classifications.csv from tie_analysis.py -- keyed by
    (position, turn), giving every multi-way-tied position's OWN verdict,
    not just the ones that made it into the deduplicated findings list."""
    result = {}
    if not path or not os.path.exists(path):
        return result
    with open(path, encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            result[(row['Position'], row['Turn'])] = row['Category']
    return result


def load_db(path):
    db = {}
    with open(path, encoding='utf-8') as f:
        f.readline()
        for line in f:
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
# Tree rendering
# ============================================================================

def render(state, db, indent, ply_num, out, max_lines, visiting, classifications=None):
    """Recursively renders the perfect-play tree from `state`. Forced
    (single-move) stretches print as one line each; a real tie prints an
    explicit branch header with each alternative labeled and indented.
    If `classifications` is provided (from all_classifications.csv), every
    branch point -- not just the root -- gets its own EXPLAINED/INDEPENDENT/
    TRIVIAL_MATE label, since every multi-way tie in the database was
    already classified when tie_analysis.py ran over the whole file, not
    just the ones that made it into the deduplicated findings list.
    `visiting` guards against re-entering a state already on the current
    path (defensive only -- the game strictly decreases in remaining plies
    every move, so a true cycle should be structurally impossible; this
    just turns a data anomaly into a visible note instead of infinite
    recursion if one ever existed)."""
    if len(out) >= max_lines:
        out.append("  " * indent + "... [text view truncated at max_lines -- "
                                    "see the exhaustive path count above, which is NOT truncated]")
        return

    if state in visiting:
        out.append("  " * indent + "# [anomaly: revisited a position already on this line -- "
                                    "should be impossible, flag this data]")
        return

    pos, turn = state
    entry = db.get(state)
    if entry is None:
        if is_checkmate_state(state):
            out.append("  " * indent + "# checkmate")
        else:
            out.append("  " * indent + "# [GAP: position not in database -- unexplored, not a real dead end]")
        return

    tied = entry['tied']
    if not tied:
        out.append("  " * indent + "# [no tied moves recorded here -- terminal or gap]")
        return

    if len(tied) == 1:
        mv, bn = tied[0]
        out.append("  " * indent + f"{ply_num}. {mv}")
        child = (apply_move(pos, mv), child_turn(mv))
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
            child = (apply_move(pos, mv), child_turn(mv))
            render(child, db, indent + 2, ply_num + 1, out, max_lines, visiting | {state}, classifications)


def render_position(position, turn, db, out_path, component_labels=None, max_lines=500,
                     classifications=None):
    key = (position, turn)
    entry = db.get(key)
    header = [
        f"Position: {position}  ({turn} to move)",
    ]
    if entry:
        header.append(f"Total plies to mate: {entry['total_plies']}")
        header.append(f"Root tied moves: {len(entry['tied'])}")
        if component_labels:
            header.append("Root classification (from tie_analysis.py):")
            for group_idx, moves in enumerate(component_labels):
                tag = "INDEPENDENT component" if len(component_labels) > 1 else "component"
                header.append(f"  {tag} {group_idx + 1}: {', '.join(moves)}")

        # Exhaustive, memoized (NOT brute-force, NOT truncated) count of every
        # complete tied-optimal path from this position to mate, with an
        # independent recomputation of whether they really do all share the
        # same total cumulative escape count -- see verify_and_count's own
        # comment for why this is cheap even when the true path count is huge.
        mismatches = []
        leaf_count, verified_total = verify_and_count(key, db, {}, mismatches)
        header.append("")
        header.append(f"EXHAUSTIVE complete perfect-play paths from this position: {leaf_count}")
        if mismatches:
            header.append(f"WARNING - INVARIANT VIOLATION: {len(mismatches)} node(s) had tied choices "
                           f"that do NOT all sum to the same total escape count -- this should be "
                           f"mathematically impossible given how TiedMoves is constructed, and means "
                           f"something in the data or this check is wrong. First violation:")
            state, vals = mismatches[0]
            header.append(f"  at {state}: {vals}")
        else:
            header.append(f"Verified: every one of those {leaf_count} paths sums to the exact "
                           f"same total cumulative escape count ({verified_total}), independently "
                           f"recomputed from scratch -- not read from the stored bn_cum values.")
    else:
        header.append("# [position not found in database]")
    header.append("")

    out = []
    if entry:
        render(key, db, 0, 1, out, max_lines, set(), classifications)
    if len(out) >= max_lines:
        header.append(f"NOTE: the line-by-line tree below is capped at {max_lines} lines for "
                       f"readability and does NOT show all {leaf_count} paths -- the counts and "
                       f"verification above are exhaustive regardless.")
        header.append("")
    text = "\n".join(header + out) + "\n"

    with open(out_path, "w", encoding='utf-8') as f:
        f.write(text)
    return text


# ============================================================================
# CLI
# ============================================================================

def parse_component_groups(groups_str):
    """Parses the ComponentGroups column from independent_findings_deduped.csv:
    groups separated by '|', moves within a group separated by ','."""
    if not groups_str:
        return None
    return [g.split(',') for g in groups_str.split('|')]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("db_path")
    ap.add_argument("--position", help='e.g. "WK:a4 WQ:e4 BK:h4"')
    ap.add_argument("--turn", choices=["W", "B"])
    ap.add_argument("--out", default=None, help="Output file for a single --position render")
    ap.add_argument("--findings-csv", help="independent_findings_deduped.csv from tie_analysis.py")
    ap.add_argument("--top", type=int, default=10, help="Render this many deepest findings from --findings-csv")
    ap.add_argument("--out-dir", default="trees", help="Output directory when using --findings-csv")
    ap.add_argument("--max-lines", type=int, default=500, help="Safety cap on rendered lines per position "
                     "(the exhaustive path count is never capped, only the line-by-line text view)")
    ap.add_argument("--classifications", help="all_classifications.csv from tie_analysis.py -- "
                     "annotates EVERY branch point in the tree with its own verdict, not just the root")
    args = ap.parse_args()

    print(f"Loading {args.db_path}...")
    db = load_db(args.db_path)
    print(f"  Loaded {len(db)} positions")

    classifications = load_classifications(args.classifications) if args.classifications else None
    if args.classifications:
        print(f"  Loaded {len(classifications)} classifications for branch-point annotation")

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
            position = row['Position']
            turn = row['Turn']
            components = parse_component_groups(row.get('ComponentGroups', ''))
            safe_name = position.replace(':', '').replace(' ', '_')
            out_path = os.path.join(args.out_dir, f"{i+1:03d}_{safe_name}_{turn}.txt")
            render_position(position, turn, db, out_path, component_labels=components,
                             max_lines=args.max_lines, classifications=classifications)
            print(f"  [{i+1}/{len(rows)}] {position} ({turn}, {row['TotalPlies']} plies) -> {out_path}")
        print(f"\nWrote {len(rows)} tree files to {args.out_dir}/")
        return

    print("ERROR: specify either --position/--turn or --findings-csv", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
