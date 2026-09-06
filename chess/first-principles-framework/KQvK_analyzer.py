#!/usr/bin/env python3
"""
KQvK Perfect-Play Analyzer
=============================
One tool, three stages, one shared codebase.

  classify   Runs over the WHOLE landscape database (from --batch or
             --full-dag), classifies every multi-way tie as TRIVIAL_MATE /
             EXPLAINED / INDEPENDENT / INDEPENDENT_UNVERIFIED. Writes
             all_classifications.csv and independent_findings_deduped.csv.

  families   Runs over classify's independent_findings_deduped.csv, groups
             recurring non-symmetric patterns into FAMILY vs true UNIQUE
             findings. Pass --render-trees DB_PATH to also render the full,
             exhaustive move tree for every unique finding (and one
             representative per family) in the same run -- no separate
             invocation needed to go from "here's the list" to "here's what
             it actually looks like."

  tree       Renders the full, exhaustive, invariant-verified move tree for
             one specific position, or for the top N deepest entries from
             any findings CSV.

This replaces tie_analysis.py, find_perfect_play_families.py, and
render_perfect_play_tree.py as three separate files -- same logic, ported
directly (not rewritten from scratch), now sharing one set of board-geometry
and database-loading functions instead of three drifting copies of them.

Usage:
    python3 kqvk_analyzer.py classify full_dag.db --out-dir tie_analysis

    python3 kqvk_analyzer.py families tie_analysis/independent_findings_deduped.csv \\
        --out-dir families --render-trees full_dag.db

    python3 kqvk_analyzer.py tree full_dag.db \\
        --position "WK:a4 WQ:e4 BK:h4" --turn W --out tree.txt \\
        --classifications tie_analysis/all_classifications.csv
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

FILES = "abcdefgh"


# ============================================================================
# Shared board geometry / D4 symmetry group -- used by all three subcommands,
# defined exactly once.
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
        parts[k] = v
    return parts


def parse_position_coords(pos_str):
    parts = {}
    for tok in pos_str.split():
        k, v = tok.split(':')
        parts[k] = sq_to_coord(v)
    return parts


def position_triple_coords(pos_str):
    p = parse_position_coords(pos_str)
    return (p['WK'], p['WQ'], p['BK'])


def canonical_form(triple):
    """Smallest image of (wk,wq,bk) coords under the 8-element D4 group."""
    best = None
    for T in TRANSFORMS:
        img = tuple(T(f, r) for f, r in triple)
        if best is None or img < best:
            best = img
    return best


def apply_move(pos_str, mv):
    parts = parse_position(pos_str)
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
    triple = position_triple_coords(pos)
    tt = tuple(T(f, r) for f, r in triple)
    new_pos = f"WK:{coord_to_sq(*tt[0])} WQ:{coord_to_sq(*tt[1])} BK:{coord_to_sq(*tt[2])}"
    return (new_pos, turn)


def is_symmetric_triple(triple_a, triple_b):
    """triple_* = (wk, wq, bk) as squares. Returns the transform index (0 =
    identity, i.e. literally the same position) or None."""
    ca = tuple(sq_to_coord(s) for s in triple_a)
    cb = tuple(sq_to_coord(s) for s in triple_b)
    for t_idx, T in enumerate(TRANSFORMS):
        if tuple(T(f, r) for f, r in ca) == cb:
            return t_idx
    return None


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
        p = parse_position_coords(pos_str)
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


def load_db(path):
    """Returns dict (position, turn) -> {'total_plies': int, 'tied': [(mv,bn),...]}"""
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
# STAGE 1: classify  (ported from tie_analysis.py)
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


def reachable_set(state, db, memo, max_nodes=5000, warned=None):
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
                truncated = True
            continue
        for mv, bn in entry['tied']:
            pos, turn = s
            child = (apply_move(pos, mv), child_turn(mv))
            if child not in seen:
                stack.append(child)
    result = (frozenset(seen), truncated)
    memo[state] = result
    return result


def classify_position(position, turn, entry, db, reach_memo, warned):
    tied = entry['tied']
    n = len(tied)
    if n <= 1:
        return None
    if entry['total_plies'] == 1:
        return {'category': 'TRIVIAL_MATE', 'components': None, 'n_components': 1}

    children = [(apply_move(position, mv), child_turn(mv)) for mv, bn in tied]
    triples = [position_triple_coords(c[0]) for c in children]

    uf = UnionFind(n)
    for i in range(n):
        for j in range(i + 1, n):
            if children[i][1] != children[j][1]:
                continue
            for T in TRANSFORMS[1:]:
                img = tuple(T(f, r) for f, r in triples[i])
                if img == triples[j]:
                    uf.union(i, j)
                    break

    branch_results = [reachable_set(c, db, reach_memo, warned=warned) for c in children]
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
    for (position, turn), entry in db.items():
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


def summarize_classify(results, out_dir):
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
    if total_classified:
        frac = n_unverified / total_classified
        print(f"\n{n_unverified}/{total_classified} ({frac:.1%}) of ALL classified positions "
              f"are INDEPENDENT_UNVERIFIED -- apparent independence indistinguishable from a "
              f"missing branch. Large for --batch-sourced data (expected), should be small for "
              f"a fully-explored --full-dag file.")
        if frac > 0.5:
            print("  WARNING: over half of all classified positions are unverified. If this file "
                  "came from --batch mode, that's expected -- --batch never explores sibling "
                  "branches. Treat INDEPENDENT counts as a lower bound, not a measurement.")

    independent = [r for r in results if r['category'] == 'INDEPENDENT']
    seen_canon = {}
    deduped = []
    for r in independent:
        triple = position_triple_coords(r['position'])
        canon = canonical_form(triple)
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
        f.write("Position,Turn,TotalPlies,NumComponents,TiedMoves\n")
        for r in sorted(unverified, key=lambda x: -x['total_plies']):
            tied_str = ";".join(f"{mv}:{bn}" for mv, bn in r['tied'])
            f.write(f'"{r["position"]}",{r["turn"]},{r["total_plies"]},'
                    f'{r["n_components"]},"{tied_str}"\n')
    print(f"Wrote {len(unverified)} INDEPENDENT_UNVERIFIED (data-limited) positions to "
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

    multi_count = sum(1 for e in db.values() if len(e['tied']) > 1)
    print(f"  Positions with a multi-way tie (any ply-depth): {multi_count}")

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
            print("  WARNING: this looks like it may be --batch-sourced (or a partial --full-dag "
                  "run) -- expect INDEPENDENT_UNVERIFIED to dominate the results below.")

    results = run_categorization(db, min_plies=args.min_plies, max_positions=args.max_positions)
    summarize_classify(results, args.out_dir)


# ============================================================================
# STAGE 3 code needed by STAGE 2's --render-trees integration: tree rendering
# (ported from render_perfect_play_tree.py). Defined before `families` so
# cmd_families can call it directly.
# ============================================================================

def verify_and_count(state, db, memo, mismatches):
    if state in memo:
        return memo[state]
    entry = db.get(state)
    tied = entry['tied'] if entry else []
    if not tied:
        result = (1, 0)
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
            p = parse_position_coords(child_pos)
            contribution = black_legal_move_count(p['WK'], p['WQ'], p['BK'])
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
        out.append("  " * indent + ("# checkmate" if is_checkmate_state(state)
                                     else "# [GAP: position not in database -- unexplored, not a real dead end]"))
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
    header = [f"Position: {position}  ({turn} to move)"]
    if entry:
        header.append(f"Total plies to mate: {entry['total_plies']}")
        header.append(f"Root tied moves: {len(entry['tied'])}")
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
                           f"choices that do NOT all sum to the same total escape count -- this "
                           f"should be mathematically impossible. First violation:")
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
                       f"readability and does NOT show all paths -- the counts above are exhaustive.")
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
# STAGE 2: families  (ported from find_perfect_play_families.py), with the
# new --render-trees integration so "compile the list" and "look at one"
# don't require a separate invocation.
# ============================================================================

def load_findings(path):
    findings = []
    with open(path, encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            parts = parse_position(row['Position'])
            tied_moves = set()
            for entry in row['TiedMoves'].split(';'):
                if ':' in entry:
                    tied_moves.add(entry.rsplit(':', 1)[0])
            findings.append({
                'position': row['Position'], 'turn': row['Turn'],
                'total_plies': int(row['TotalPlies']) if row['TotalPlies'] else 0,
                'WK': parts.get('WK'), 'WQ': parts.get('WQ'), 'BK': parts.get('BK'),
                'tied_moves': frozenset(tied_moves),
                'component_groups': row.get('ComponentGroups', ''),
            })
    return findings


def find_families(findings):
    seen_in_a_family = set()
    families = []
    for fix_piece in ('WK', 'WQ', 'BK'):
        groups = defaultdict(list)
        for finding in findings:
            key = (finding[fix_piece], finding['turn'], finding['tied_moves'])
            groups[key].append(finding)
        for (fixed_val, turn, tied_moves), members in groups.items():
            distinct_positions = {(m['position'], m['turn']): m for m in members}
            if len(distinct_positions) < 2:
                continue
            member_list = list(distinct_positions.values())
            triples = [(m['WK'], m['WQ'], m['BK']) for m in member_list]
            symmetric_pairs = []
            for i in range(len(triples)):
                for j in range(i + 1, len(triples)):
                    t_idx = is_symmetric_triple(triples[i], triples[j])
                    if t_idx is not None and t_idx != 0:
                        symmetric_pairs.append((member_list[i]['position'], member_list[j]['position'], t_idx))
            families.append({
                'fixed_piece': fix_piece, 'fixed_value': fixed_val, 'turn': turn,
                'tied_moves': tied_moves, 'members': member_list,
                'symmetric_pairs': symmetric_pairs,
            })
            for m in member_list:
                seen_in_a_family.add((m['position'], m['turn']))
    unique_findings = [f for f in findings if (f['position'], f['turn']) not in seen_in_a_family]
    return families, unique_findings


def cmd_families(args):
    print(f"Loading {args.findings_csv}...")
    findings = load_findings(args.findings_csv)
    print(f"  Loaded {len(findings)} deduplicated independent findings")

    families, unique_findings = find_families(findings)

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
    print(f"\nFamilies found (recurring, non-symmetric, describable patterns): {len(distinct_families)}")
    print(f"Unique findings (no recurrence under any single-piece-fixed sweep): {len(unique_findings)}")

    unexpected_symmetric = [f for f in distinct_families if f['symmetric_pairs']]
    if unexpected_symmetric:
        print(f"\nWARNING: {len(unexpected_symmetric)} 'family' actually contained a nontrivial "
              f"symmetric pair (expected only when the fixed piece sits on a main diagonal):")
        for fam in unexpected_symmetric[:5]:
            print(f"  fixed {fam['fixed_piece']}={fam['fixed_value']}: {fam['symmetric_pairs']}")

    with open(os.path.join(args.out_dir, "families.csv"), "w", encoding='utf-8') as f:
        f.write("FixedPiece,FixedValue,Turn,TiedMoves,FamilySize,MemberPositions,ContainsSymmetricPair\n")
        for fam in distinct_families:
            moves_str = ";".join(sorted(fam['tied_moves']))
            members_str = "|".join(m['position'] for m in fam['members'])
            f.write(f'{fam["fixed_piece"]},{fam["fixed_value"]},{fam["turn"]},"{moves_str}",'
                    f'{len(fam["members"])},"{members_str}",{bool(fam["symmetric_pairs"])}\n')
    print(f"Wrote {len(distinct_families)} families to {args.out_dir}/families.csv")

    with open(os.path.join(args.out_dir, "unique_findings.csv"), "w", encoding='utf-8') as f:
        f.write("Position,Turn,TotalPlies,TiedMoves,ComponentGroups\n")
        for u in sorted(unique_findings, key=lambda x: -x['total_plies']):
            moves_str = ";".join(sorted(u['tied_moves']))
            f.write(f'"{u["position"]}",{u["turn"]},{u["total_plies"]},"{moves_str}","{u["component_groups"]}"\n')
    print(f"Wrote {len(unique_findings)} unique findings to {args.out_dir}/unique_findings.csv "
          f"-- THESE are the genuinely singular cases")

    print("\nLargest families:")
    for fam in distinct_families[:10]:
        moves_str = "/".join(sorted(fam['tied_moves']))
        print(f"  fix {fam['fixed_piece']}={fam['fixed_value']} ({fam['turn']}), moves={{{moves_str}}}: "
              f"{len(fam['members'])} positions")

    print("\nDeepest unique (truly one-off) findings:")
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
    p1.add_argument('--out-dir', default='tie_analysis')

    p2 = sub.add_parser('families', help="Compile families vs unique findings from classify's output")
    p2.add_argument('findings_csv')
    p2.add_argument('--out-dir', default='families')
    p2.add_argument('--render-trees', metavar='DB_PATH', default=None,
                     help="Also render full move trees for every unique finding and one "
                          "representative per family, using this database file")
    p2.add_argument('--max-lines', type=int, default=500)

    p3 = sub.add_parser('tree', help="Render the full move tree for one position or a findings CSV")
    p3.add_argument('db_path')
    p3.add_argument('--position')
    p3.add_argument('--turn', choices=['W', 'B'])
    p3.add_argument('--out', default=None)
    p3.add_argument('--findings-csv')
    p3.add_argument('--top', type=int, default=10)
    p3.add_argument('--out-dir', default='trees')
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
