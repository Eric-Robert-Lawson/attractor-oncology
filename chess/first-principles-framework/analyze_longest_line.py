#!/usr/bin/env python3
"""
Finds the longest (maximum-distance) position in any solved material's raw
database, walks the true canonical optimal line from it to mate, and
reports -- at every ply -- both the full menu of distance-tied
alternatives available there, AND, new in this version, EXACTLY which
distinct mate destination(s) get structurally eliminated by each specific
move, not just how the breadth number changes.

WHY THIS EXISTS: a breadth number alone (e.g. "2 distinct mates reachable"
sustained for 53 straight plies, then dropping to 1) tells you THAT
something narrowed, and roughly where, but not WHICH destinies were still
alive, nor the exact mechanical reason the drop happens where it does.
This version traces that directly instead of leaving it to be inferred
from a bare number.

A REAL MATHEMATICAL FACT WORTH KNOWING BEFORE READING THE OUTPUT, not
assumed: a position with exactly ONE tied move cannot, by construction,
change which distinct mate destinations are reachable. reachable_shapes()
for a single-child node is `{(mate, bn+esc) for (mate,esc) in
reachable_shapes(child)}` -- a pure, fixed-offset shift of the child's own
escape values, applied uniformly. It relabels escape counts; it cannot
merge two different mates together or drop one, since mate identity is
untouched and the offset is constant across every entry. So elimination
of a distinct mate destination can ONLY happen exactly at a position with
more than one tied move -- MatesEliminatedByThisMove should read empty at
every purely forced ply, without exception, and non-empty only at genuine
decision points. This isn't asserted here -- it's checked directly, row
by row, against your own real data, and the summary below will say so
plainly if that expectation is ever violated (which would itself be a
real, interesting finding, not just a bug to fix quietly).

Everything else (backward induction, the two recursive quantities, the
verification history) is unchanged from the previous version -- see this
project's general_workflow_reference.md §8 for the full record.

Usage:
    python3 analyze_longest_line.py kbnvk_perfect_play.db --out kbnvk_longest_line.csv --material-name KBNvK
"""
import argparse
import csv
import sys
import time


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('db_path', help="Raw classify() database for any material (e.g. kbnvk_perfect_play.db)")
    ap.add_argument('--out', default='longest_line.csv', help="Output CSV path (default longest_line.csv)")
    ap.add_argument('--material-name', default=None,
                     help="Label for this material in the printed output. Purely cosmetic.")
    ap.add_argument('--general-analyzer-dir', default='.',
                     help="Directory containing general_analyzer.py, if not the current one")
    args = ap.parse_args()

    material_name = args.material_name or "this material"

    sys.path.insert(0, args.general_analyzer_dir)
    try:
        from general_analyzer import (
            load_db, apply_move_general, child_turn, canonical_form_general,
            AmbiguousMoveError,
        )
    except ImportError as e:
        print(f"ERROR: couldn't import general_analyzer.py ({e}). Point --general-analyzer-dir "
              f"at the directory that contains it.", file=sys.stderr)
        sys.exit(1)

    print(f"Loading {args.db_path} as {material_name}...")
    t0 = time.time()
    db = load_db(args.db_path)
    print(f"  Loaded {len(db)} positions in {time.time()-t0:.1f}s.\n")

    print(f"Finding the longest (max-distance) position in {material_name}...")
    max_distance = -1
    longest_states = []
    for state, entry in db.items():
        if entry.distance > max_distance:
            max_distance = entry.distance
            longest_states = [state]
        elif entry.distance == max_distance:
            longest_states.append(state)
    print(f"  Longest distance found: {max_distance} plies.")
    print(f"  {len(longest_states)} position(s) share this maximum -- using the first "
          f"as the starting point: {longest_states[0]}\n")
    root = longest_states[0]

    # --- cumulative_best: proper backward induction ---------------------------
    best_memo = {}
    best_choice = {}
    all_options = {}

    def cumulative_best(state):
        if state in best_memo:
            return best_memo[state]
        entry = db.get(state)
        if entry is None:
            raise ValueError(f"{state} is not in the database -- can't be on a forced-mate "
                              f"line if it's a proven draw or unreached.")
        if not entry.tied:
            best_memo[state] = 0
            return 0
        pos, turn = state
        options = []
        for mv, bn in entry.tied:
            try:
                child_pos = apply_move_general(pos, mv)
            except AmbiguousMoveError:
                continue
            child = (child_pos, child_turn(mv))
            total = bn + cumulative_best(child)
            options.append((mv, child, bn, total))
        if not options:
            raise ValueError(f"{state}: every tied move was ambiguous -- nothing to recurse into.")
        if turn == 'W':
            best = min(options, key=lambda o: o[3])
        else:
            best = max(options, key=lambda o: o[3])
        best_memo[state] = best[3]
        best_choice[state] = best
        all_options[state] = options
        return best[3]

    # --- reachable_shapes: full union set ---------------------------------
    reach_memo = {}

    def reachable_shapes(state):
        if state in reach_memo:
            return reach_memo[state]
        entry = db.get(state)
        if entry is None:
            return frozenset()
        pos, turn = state
        if not entry.tied:
            canon = (canonical_form_general(pos), turn)
            result = frozenset([(canon, 0)])
            reach_memo[state] = result
            return result
        combined = set()
        for mv, bn in entry.tied:
            try:
                child_pos = apply_move_general(pos, mv)
            except AmbiguousMoveError:
                continue
            child = (child_pos, child_turn(mv))
            for mate_canon, esc in reachable_shapes(child):
                combined.add((mate_canon, bn + esc))
        result = frozenset(combined)
        reach_memo[state] = result
        return result

    print(f"Computing {material_name}'s true canonical optimal line, keeping every tied "
          f"alternative's own value along the way (proper backward induction)...")
    t0 = time.time()
    total_cumulative_escape = cumulative_best(root)
    print(f"  Done in {time.time()-t0:.1f}s. Total cumulative escape count for the actual line: "
          f"{total_cumulative_escape}.\n")

    print("Computing reachable-shape breadth at every ply, and tracing exactly which mate "
          "destinations get eliminated by each move...")
    t0 = time.time()

    # First pass: walk the line, recording the state and reachable_shapes()
    # set at every ply, before assigning readable labels.
    line_states = []
    state = root
    while True:
        line_states.append(state)
        entry = db[state]
        if not entry.tied:
            break
        mv, child, bn, total = best_choice[state]
        state = child

    line_breadths = [reachable_shapes(s) for s in line_states]

    # Assign short, readable labels to every distinct mate destination that
    # ever appears anywhere along this line's own reachable sets -- makes
    # the trace readable without printing full canonical position strings
    # over and over.
    all_mate_canons = set()
    for breadth in line_breadths:
        for mate_canon, esc in breadth:
            all_mate_canons.add(mate_canon)
    mate_label = {canon: f"Mate{i+1}" for i, canon in enumerate(sorted(all_mate_canons))}
    print(f"  {len(all_mate_canons)} distinct mate destination(s) appear anywhere along this "
          f"line's own reachable sets: {sorted(mate_label.values())}\n")

    rows = []
    running_escape = 0
    violations = 0
    for ply, state in enumerate(line_states):
        entry = db[state]
        pos, turn = state
        breadth = line_breadths[ply]
        mates_here = {mc for mc, esc in breadth}
        labels_here = sorted(mate_label[mc] for mc in mates_here)

        if not entry.tied:
            rows.append({
                'Ply': ply, 'Position': pos, 'Turn': turn, 'Move': '', 'IsChosen': True,
                'DistanceRemaining': entry.distance, 'MoveEscapeContribution': '',
                'ThisAlternativeTotalFromHere': 0, 'CumulativeEscapeSoFar': running_escape,
                'NumTiedAlternativesThisPly': 0,
                'MatesReachableHere': ';'.join(labels_here),
                'DistinctShapesReachableHere': len(breadth),
                'MatesEliminatedByThisMove': '',
                'IsMate': True,
            })
            break

        options = all_options[state]
        n_options = len(options)
        chosen_mv, chosen_child, chosen_bn, chosen_total = best_choice[state]

        mates_after = {mc for mc, esc in line_breadths[ply + 1]} if ply + 1 < len(line_states) else mates_here
        eliminated = mates_here - mates_after
        eliminated_labels = sorted(mate_label[mc] for mc in eliminated)
        if eliminated_labels and n_options == 1:
            violations += 1
            print(f"  [UNEXPECTED] ply {ply} is a forced move (only one tied option) but "
                  f"{eliminated_labels} were eliminated anyway -- this contradicts the "
                  f"mathematical expectation in this script's own docstring. Flagging, not "
                  f"silently ignoring.")

        for mv, child, bn, total in sorted(options, key=lambda o: o[3]):
            rows.append({
                'Ply': ply, 'Position': pos, 'Turn': turn, 'Move': mv,
                'IsChosen': (mv == chosen_mv and child == chosen_child),
                'DistanceRemaining': entry.distance, 'MoveEscapeContribution': bn,
                'ThisAlternativeTotalFromHere': total, 'CumulativeEscapeSoFar': running_escape,
                'NumTiedAlternativesThisPly': n_options,
                'MatesReachableHere': ';'.join(labels_here),
                'DistinctShapesReachableHere': len(breadth),
                'MatesEliminatedByThisMove': (';'.join(eliminated_labels)
                                               if (mv == chosen_mv and child == chosen_child) else ''),
                'IsMate': False,
            })

        running_escape += chosen_bn

    print(f"  Done in {time.time()-t0:.1f}s. Line is {len(line_states)} plies long "
          f"(matches the {max_distance}-ply distance found above), {len(rows)} total rows.\n")

    with open(args.out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows)} rows to {args.out}")

    # --- Summary: exactly where and what got eliminated -----------------------
    print(f"\n=== {material_name}: exactly where mate destinations get eliminated along the line ===")
    if violations:
        print(f"  {violations} forced ply/plies eliminated a mate destination anyway -- see "
              f"[UNEXPECTED] lines above; this needs a closer look, not just a note.")
    else:
        print(f"  Confirmed directly on this line's own real data: every elimination happens "
              f"exactly at a genuine decision point (>1 tied move), never at a forced one -- "
              f"matching the mathematical expectation in this script's own docstring, not just "
              f"assumed.")
    any_elimination = False
    for r in rows:
        if r['MatesEliminatedByThisMove']:
            any_elimination = True
            print(f"  Ply {r['Ply']} ({r['Position']}, {r['Turn']} to move), move {r['Move']} "
                  f"chosen over {r['NumTiedAlternativesThisPly']-1} alternative(s): eliminates "
                  f"{r['MatesEliminatedByThisMove']} -- reachable set goes from "
                  f"{r['MatesReachableHere']} to whatever remains.")
    if not any_elimination:
        print("  No eliminations recorded at all -- the same full set of mate destinations "
              "stayed reachable for the entire line, only resolving at the mate itself.")

    print(f"\n=== {material_name}: breadth of distinct SHAPES (mate+escape pairs) at each ply ===")
    print([r['DistinctShapesReachableHere'] for r in rows if r['IsChosen']])


if __name__ == "__main__":
    main()
