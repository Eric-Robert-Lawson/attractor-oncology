#!/usr/bin/env python3
"""
Segmented, hybrid version of build_shape_index.py: a small, centralized
routing index plus one SQLite file per mate, instead of one monolithic
database holding every shape.

WHY THIS EXISTS, AND WHY IT'S NOT JUST "SMALLER FILES ARE BETTER":

Tested directly, both directions, on real KQvK data, before building this
for real: for random access across many different shapes, the monolithic,
single-file approach won by 3.0x (0.026ms/query vs 0.080ms/query) --
opening a fresh file per query has real, measurable OS-level overhead that
one shared, cached B-tree index inside an already-open file doesn't pay.
But for high-locality access -- many queries within ONE shape before
moving to another, which is what actually walking a trajectory looks like
-- segmentation won by 1.9x (0.007ms/query vs 0.013ms/query), since once
the smaller file is open, searching it is faster than searching the whole
shared index. This tool exists for that second access pattern
specifically. build_shape_index.py's original monolithic output is still
the right choice for random, scattered lookups across many shapes.

SEGMENTED BY MATE, NOT BY THE FULL (MATE, ESCAPE-COUNT) SHAPE IDENTITY:

A position with a single mate ahead of it can still have several distinct
escape-count variants (measured on KQvK: every one of 46 mates has more
than one, one mate has 133). Splitting at that finer granularity would
mean an edge's own file assignment depends on knowing the exact remaining
escape budget at that specific ply -- information the current
shape_graph_edges.csv doesn't carry (it has no escape-count column at
all), so a genuinely correct escape-exact split isn't derivable from this
tool's own inputs without re-deriving it from the original database. Mate-
level segmentation avoids that ambiguity entirely -- an edge unambiguously
belongs to every mate its CHILD position's ShapeIds names (not the FROM
position's -- see the real, confirmed bug this fixed, in pass 2's own
comment below) -- and matches
the granularity --trace-shape-ids already lets a person pick shapes at.
The real cost: a mate's own file bundles all of its escape-count variants
together, not just the one narrow (mate, escape) slice measured at ~4-10KB
in earlier testing -- still dramatically smaller than the full graph, just
not quite that small.

MEMORY-BOUNDED BY DESIGN, NOT JUST BY DEFAULT:

Two real, separate memory risks were designed around directly rather than
discovered the hard way (the way find_shape_origins and
resolve_reachable_shapes both were, earlier in this project):

1. Holding every distinct mate's own SQLite connection open at once,
   unbounded, doesn't scale -- a material with many thousands of mates
   would exhaust file descriptors or memory long before finishing.
   BoundedConnectionCache below mirrors general_analyzer.py's own
   BoundedLRUCache pattern (already used for reachable_set memoization),
   applied to file handles instead of computed values: at most
   --max-open-files connections stay open at once, least-recently-used
   evicted (committed and closed, not discarded) when a new one is
   needed beyond the cap, and correctly re-opened in append mode if that
   mate is needed again later.

2. Assigning each edge to its FROM position's mate(s) needs to know
   those mates -- but holding a full position -> mates dict in memory
   for that purpose would itself scale with the whole database's
   position count, the exact shape of cost this project already hit and
   fixed twice before (find_shape_origins' parent_counts,
   resolve_reachable_shapes' resolved dict). Avoided entirely here by
   using the routing index ITSELF as that lookup -- fully built and
   committed in pass 1, then queried (not held in memory) during pass 2
   -- trading a disk-backed indexed query per edge for zero unbounded
   in-memory growth.

Usage:
    python3 build_segmented_shape_index.py material_reduction/shape_graph_nodes.csv \\
        material_reduction/shape_graph_edges.csv --out-dir material_shapes_segmented \\
        --routing material_routing.sqlite
"""

import argparse
import csv
import os
import sqlite3
import time
from collections import OrderedDict


class BoundedConnectionCache:
    """At most max_open SQLite connections open at once. Mirrors
    general_analyzer.py's own BoundedLRUCache -- same reasoning:
    unbounded resource growth on a long run at real scale is a real,
    previously-confirmed cause of trouble in this project, not a
    hypothetical one being guarded against preemptively.

    A REAL, CONFIRMED FAILURE MODE THIS MUST BE SIZED AGAINST: if the
    number of distinct mates in the material exceeds max_open, and
    entries in the source CSVs are not grouped by mate (they aren't --
    both are in position order), this cache thrashes -- constantly
    evicting and reopening connections, since each new row's mate is
    rarely the one just evicted. Measured directly on KQvK (46 mates):
    max_open=20 dropped throughput to ~200 rows/second (severe,
    confirmed thrashing -- a full run would have taken ~29 minutes for
    node-loading alone); max_open=100 (comfortably above 46) processed
    the same data in under 2 seconds. This is not a tuning nicety, it's
    the difference between the tool being usable and not. eviction_warned
    below fires once, loudly, the first time evictions start happening
    at a rate consistent with this exact problem -- so it's caught
    immediately at real scale, not discovered 29 minutes in."""

    def __init__(self, out_dir, max_open=100):
        self.out_dir = out_dir
        self.max_open = max_open
        self._conns = OrderedDict()  # mate_id -> (conn, cursor)
        self.opened_files = set()
        self.eviction_count = 0
        self.eviction_warned = False

    def get(self, mate_id):
        if mate_id in self._conns:
            conn, cur = self._conns.pop(mate_id)
            self._conns[mate_id] = (conn, cur)  # move to end: most recently used
            return cur
        if len(self._conns) >= self.max_open:
            _, (evict_conn, _) = self._conns.popitem(last=False)  # evict least recently used
            evict_conn.commit()
            evict_conn.close()
            self.eviction_count += 1
            if not self.eviction_warned and self.eviction_count >= 500:
                print(f"  WARNING: {self.eviction_count} connection evictions so far -- this "
                      f"material likely has more distinct mates than --max-open-files "
                      f"({self.max_open}), which causes severe thrashing (confirmed directly on "
                      f"KQvK: ~200 rows/sec instead of ~14,000+/sec). Consider killing this run "
                      f"and restarting with a much larger --max-open-files, at least comfortably "
                      f"above this material's actual distinct-mate count (check shapes.csv's own "
                      f"row count, or the number of distinct ShapeIds seen, for an estimate).",
                      flush=True)
                self.eviction_warned = True
        path = os.path.join(self.out_dir, f"shape_{mate_id}.sqlite")
        is_new = not os.path.exists(path)
        conn = sqlite3.connect(path)
        cur = conn.cursor()
        if is_new:
            cur.execute("PRAGMA journal_mode=WAL")
            cur.execute("CREATE TABLE nodes (position TEXT, turn TEXT, is_origin INTEGER, "
                        "PRIMARY KEY (position, turn))")
            cur.execute("CREATE TABLE node_shapes (position TEXT, turn TEXT, "
                        "shape_id INTEGER, escape_count INTEGER)")
            cur.execute("CREATE TABLE edges (position TEXT, turn TEXT, move TEXT, "
                        "child_position TEXT, child_turn TEXT)")
        self._conns[mate_id] = (conn, cur)
        self.opened_files.add(mate_id)
        return cur

    def close_all(self):
        for conn, cur in self._conns.values():
            conn.commit()
            conn.close()
        self._conns.clear()


def build_segmented_index(nodes_csv, edges_csv, out_dir, routing_path, max_open_files=100,
                           progress_every=500000):
    os.makedirs(out_dir, exist_ok=True)

    print(f"Pre-scan: counting distinct mates in {nodes_csv} before any writing begins, "
          f"specifically to catch a real, confirmed failure mode early rather than partway "
          f"through a slow run -- see BoundedConnectionCache's own docstring for the measured "
          f"~70x slowdown this causes when --max-open-files is too small relative to the "
          f"material's actual mate count.")
    t0 = time.time()
    distinct_mates = set()
    with open(nodes_csv, encoding='utf-8-sig', newline='') as f:
        for row in csv.DictReader(f):
            for entry in row['ShapeIds'].split(';'):
                if entry and ':' in entry:
                    distinct_mates.add(int(entry.split(':')[0]))
    print(f"  {len(distinct_mates)} distinct mates found in {time.time()-t0:.1f}s.")
    if len(distinct_mates) > max_open_files * 0.8:
        print(f"  WARNING: --max-open-files is {max_open_files}, uncomfortably close to (or "
              f"below) the {len(distinct_mates)} distinct mates this material actually has. "
              f"This will very likely thrash -- consider re-running with --max-open-files set "
              f"well above {len(distinct_mates)} before continuing. Proceeding anyway in 5s "
              f"(Ctrl+C to stop)...", flush=True)
        time.sleep(5)

    print(f"\nPass 1/2: reading {nodes_csv}, building the routing index and each mate's own "
          f"nodes/node_shapes tables...")
    t0 = time.time()

    if os.path.exists(routing_path):
        os.remove(routing_path)
    rconn = sqlite3.connect(routing_path)
    rcur = rconn.cursor()
    rcur.execute("PRAGMA journal_mode=WAL")
    rcur.execute("CREATE TABLE routing (position TEXT, turn TEXT, mate_id INTEGER)")

    cache = BoundedConnectionCache(out_dir, max_open=max_open_files)
    routing_batch = []
    node_rows = 0
    with open(nodes_csv, encoding='utf-8-sig', newline='') as f:
        for row in csv.DictReader(f):
            node_rows += 1
            if progress_every and node_rows % progress_every == 0:
                print(f"  [pass 1] {node_rows} node rows processed, "
                      f"{len(cache.opened_files)} distinct mate files opened so far, "
                      f"elapsed={time.time()-t0:.0f}s", flush=True)
            is_origin = 1 if row['IsOrigin'] == 'True' else 0
            mates_seen_this_row = set()
            for entry in row['ShapeIds'].split(';'):
                if not entry:
                    continue
                if ':' not in entry:
                    raise SystemExit(
                        f"ERROR: malformed ShapeIds entry {entry!r} in {nodes_csv} -- this looks "
                        f"like a STALE file predating the escape-count fix to trace_shape_graph "
                        f"(bare '36;37' instead of the correct '36:74;37:74'). Regenerate via "
                        f"'reduce ... --trace-shape-graph' before building an index from it.")
                sid, esc = entry.split(':')
                sid, esc = int(sid), int(esc)
                cur = cache.get(sid)
                cur.execute("INSERT OR IGNORE INTO nodes VALUES (?,?,?)",
                            (row['Position'], row['Turn'], is_origin))
                cur.execute("INSERT INTO node_shapes VALUES (?,?,?,?)",
                            (row['Position'], row['Turn'], sid, esc))
                if sid not in mates_seen_this_row:
                    mates_seen_this_row.add(sid)
                    routing_batch.append((row['Position'], row['Turn'], sid))
            if len(routing_batch) >= 50000:
                rcur.executemany("INSERT INTO routing VALUES (?,?,?)", routing_batch)
                routing_batch = []
    if routing_batch:
        rcur.executemany("INSERT INTO routing VALUES (?,?,?)", routing_batch)

    cache.close_all()
    print(f"  Building routing index... ", end='', flush=True)
    ti = time.time()
    rcur.execute("CREATE INDEX idx_routing_pos ON routing(position, turn)")
    rconn.commit()
    print(f"done in {time.time()-ti:.1f}s")
    print(f"Pass 1 complete in {time.time()-t0:.1f}s: {node_rows} node rows, "
          f"{len(cache.opened_files)} distinct mate files written.")

    print(f"\nPass 2/2: reading {edges_csv}, routing each edge to its CHILD position's mate "
          f"file(s) via an indexed query against the now-complete routing index (not an "
          f"in-memory dict) ...")
    t0 = time.time()
    cache = BoundedConnectionCache(out_dir, max_open=max_open_files)
    edge_rows = 0
    with open(edges_csv, encoding='utf-8-sig', newline='') as f:
        for row in csv.DictReader(f):
            edge_rows += 1
            if progress_every and edge_rows % progress_every == 0:
                print(f"  [pass 2] {edge_rows} edge rows processed, elapsed={time.time()-t0:.0f}s",
                      flush=True)
            # Routed by the CHILD's own mates, not the FROM position's --
            # a real, confirmed bug in an earlier version routed by the
            # FROM position instead, which is wrong whenever a position
            # belongs to more mates than any single one of its own moves
            # actually leads toward (a tie whose different branches serve
            # different mates -- confirmed directly on real KQvK data:
            # Q:a1 K:a2 k:b7's own Kb3 move only leads toward 19 of that
            # position's 30 total mates; routing by the FROM position
            # would have written that edge into the other 11 mates' files
            # too, which it has nothing to do with). Since escape count
            # accumulates backward from the mate (bn + esc), a child
            # having mate M in its own reachable set is both necessary
            # and sufficient for this specific edge to genuinely be part
            # of mate M's trajectory -- matching exactly the criterion
            # the original, monolithic trace_shape_graph already uses
            # (in_scope_shape_ids(child)).
            rcur.execute("SELECT mate_id FROM routing WHERE position=? AND turn=?",
                         (row['ChildPosition'], row['ChildTurn']))
            for (sid,) in rcur.fetchall():
                cur = cache.get(sid)
                cur.execute("INSERT INTO edges VALUES (?,?,?,?,?)",
                            (row['Position'], row['Turn'], row['Move'],
                             row['ChildPosition'], row['ChildTurn']))
    cache.close_all()
    rconn.close()
    print(f"Pass 2 complete in {time.time()-t0:.1f}s: {edge_rows} edge rows routed.")

    print(f"\nBuilding per-file indices on each mate's own nodes/node_shapes/edges tables...")
    t0 = time.time()
    n_indexed = 0
    for fname in sorted(os.listdir(out_dir)):
        if not fname.endswith('.sqlite'):
            continue
        path = os.path.join(out_dir, fname)
        c = sqlite3.connect(path)
        cc = c.cursor()
        cc.execute("CREATE INDEX IF NOT EXISTS idx_node_shapes_pos ON node_shapes(position, turn)")
        cc.execute("CREATE INDEX IF NOT EXISTS idx_node_shapes_shape ON node_shapes(shape_id, escape_count)")
        cc.execute("CREATE INDEX IF NOT EXISTS idx_edges_pos ON edges(position, turn)")
        c.commit()
        c.close()
        n_indexed += 1
    print(f"Indexed {n_indexed} mate files in {time.time()-t0:.1f}s.")
    print(f"\nDone. Routing index: {routing_path}. Segmented mate files: {out_dir}/shape_*.sqlite")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('nodes_csv')
    ap.add_argument('edges_csv')
    ap.add_argument('--out-dir', required=True, help="Directory for per-mate shape_<id>.sqlite files")
    ap.add_argument('--routing', required=True, help="Output path for the routing index")
    ap.add_argument('--max-open-files', type=int, default=100,
                     help="Max simultaneously-open per-mate SQLite connections (default 100) -- "
                          "bounds memory/file-descriptor use regardless of how many distinct "
                          "mates the material has.")
    ap.add_argument('--progress-every', type=int, default=500000,
                     help="Print progress every N rows processed in each pass (default 500000, "
                          "0 disables).")
    args = ap.parse_args()
    build_segmented_index(args.nodes_csv, args.edges_csv, args.out_dir, args.routing,
                           max_open_files=args.max_open_files, progress_every=args.progress_every)


if __name__ == "__main__":
    main()
