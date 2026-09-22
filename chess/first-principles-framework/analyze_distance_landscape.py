#!/usr/bin/env python3
"""
Analyzes the real shape of a material's distance-to-mate landscape: how
many positions exist at each distance, how populated the high-distance
region actually is (a rare, isolated outlier vs. a real, richly-connected
neighborhood), and how many of the positions sharing the material's own
maximum distance are truly distinct once D4/mirror symmetry is accounted
for, rather than counted as separate when they're really the same shape
seen from a different orientation.

WHY THIS EXISTS: "Black wants to stay in the high-distance region as long
as possible" is mathematically identical to "Black maximizes distance" --
already exactly what classify() computes, not a new claim. The genuinely
open, checkable question underneath it is different: is the material's
own maximum a rare, isolated spike, or the edge of a real, well-populated
region Black could plausibly navigate within? This script answers that
directly from the same raw database, in one pass, rather than requiring
two separate ~7GB loads for what would otherwise be two separate scripts.

SYMMETRY REDUCTION, VERIFIED BEFORE BEING TRUSTED, NOT ASSUMED: a known
geometric mirror pair (K:a1 k:h8 vs its left-right reflection K:h1 k:a8)
was confirmed, directly, to canonicalize to the identical string via this
project's own canonical_form_general -- so grouping the maximal-distance
set by (canonical_form, turn) is a real, checked way to find out how many
of them are actually the same underlying shape, not assumed to work.

Usage:
    python3 analyze_distance_landscape.py kbnvk_perfect_play.db --general-analyzer-dir .
    python3 analyze_distance_landscape.py kbnvk_perfect_play.db --band-width 20 --material-name KBNvK
"""
import argparse
import sys
import time
from collections import Counter


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('db_path', help="Raw classify() database for any material")
    ap.add_argument('--material-name', default=None, help="Label for output. Purely cosmetic.")
    ap.add_argument('--general-analyzer-dir', default='.',
                     help="Directory containing general_analyzer.py, if not the current one")
    ap.add_argument('--band-width', type=int, default=20,
                     help="Report population in [max_distance - band_width, max_distance] "
                          "(default 20). Overridden by --band-min/--band-max if both given.")
    ap.add_argument('--band-min', type=int, default=None)
    ap.add_argument('--band-max', type=int, default=None)
    args = ap.parse_args()

    material_name = args.material_name or "this material"

    sys.path.insert(0, args.general_analyzer_dir)
    try:
        from general_analyzer import load_db, canonical_form_general
    except ImportError as e:
        print(f"ERROR: couldn't import general_analyzer.py ({e}). Point --general-analyzer-dir "
              f"at the directory that contains it.", file=sys.stderr)
        sys.exit(1)

    print(f"Loading {args.db_path} as {material_name}...")
    t0 = time.time()
    db = load_db(args.db_path)
    print(f"  Loaded {len(db)} positions in {time.time()-t0:.1f}s.\n")

    print("Building distance histogram...")
    t0 = time.time()
    hist = Counter(entry.distance for entry in db.values())
    max_distance = max(hist)
    print(f"  Done in {time.time()-t0:.1f}s. Distances range from 0 to {max_distance}.\n")

    print(f"=== {material_name}: population at every distance from the top down ===")
    print(f"{'Distance':>8} | {'Positions':>12} | {'% of total':>10}")
    total = len(db)
    for d in range(max_distance, -1, -1):
        count = hist.get(d, 0)
        pct = 100 * count / total if total else 0
        print(f"{d:>8} | {count:>12} | {pct:>9.4f}%")

    band_max = args.band_max if args.band_max is not None else max_distance
    band_min = args.band_min if args.band_min is not None else max(0, max_distance - args.band_width)
    band_count = sum(c for d, c in hist.items() if band_min <= d <= band_max)
    band_pct = 100 * band_count / total if total else 0
    print(f"\n=== High-distance band [{band_min}, {band_max}] ===")
    print(f"{band_count} of {total} total positions ({band_pct:.4f}%) fall in this band.")

    # Compare against the band immediately below it, same width, for direct
    # context on whether the top of the distribution is unusually sparse.
    lower_min = max(0, band_min - (band_max - band_min + 1))
    lower_max = band_min - 1
    if lower_max >= 0:
        lower_count = sum(c for d, c in hist.items() if lower_min <= d <= lower_max)
        lower_pct = 100 * lower_count / total if total else 0
        print(f"For direct comparison, the same-width band immediately below it, "
              f"[{lower_min}, {lower_max}]: {lower_count} positions ({lower_pct:.4f}%).")
        if lower_count > 0:
            print(f"Ratio: the top band has {band_count/lower_count:.3f}x as many positions "
                  f"as the band just below it.")

    # --- Symmetry reduction of the maximal-distance set specifically -------
    print(f"\n=== Symmetry reduction of the {max_distance}-ply maximal-distance set ===")
    max_states = [state for state, entry in db.items() if entry.distance == max_distance]
    print(f"{len(max_states)} raw positions share the maximum distance.")

    t0 = time.time()
    canon_groups = Counter()
    for pos, turn in max_states:
        canon_groups[(canonical_form_general(pos), turn)] += 1
    print(f"  Reduced to {len(canon_groups)} truly distinct positions under D4/mirror symmetry "
          f"in {time.time()-t0:.1f}s (a raw-to-distinct ratio of "
          f"{len(max_states)/len(canon_groups):.2f}x, for reference: full D4 symmetry with no "
          f"pawn on the board has 8 elements, so a ratio near 8 means the set is generic with "
          f"little self-symmetry; a ratio well below 8 means many of these positions are their "
          f"own partial symmetric fixed points, not full 8-fold orbits).")

    group_sizes = sorted(canon_groups.values(), reverse=True)
    print(f"  Group sizes (how many raw positions collapsed into each distinct one), "
          f"largest first: {group_sizes[:20]}{' ...' if len(group_sizes) > 20 else ''}")


if __name__ == "__main__":
    main()
