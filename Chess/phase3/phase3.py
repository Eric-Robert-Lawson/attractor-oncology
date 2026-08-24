#!/usr/bin/env python3
"""
CHESS TOPOLOGICAL SOLUTION — PHASE 3
Complex Endgames (multi-piece, tablebase-backed) + Middlegame Principle
Mining from Engine-Confirmed Games

Author: Eric Robert Lawson / OrganismCore
Continuation of: chess_topo_phase1.py, chess_topo_phase2.py
Date: 2026-08-23

────────────────────────────────────────────────────────────────
WHAT CHANGED FROM PHASE 2 (READ THIS FIRST)
────────────────────────────────────────────────────────────────
Phase 2's rook-endgame oracle fell back to a hand-written heuristic
evaluator whenever Syzygy files were absent. That heuristic computed
things like `lucena_score` and `philidor_score` directly from named
chess concepts, and was then used as the "ground truth" label to
train a classifier — which of course reported `lucena_score` as
important, because the labels were partly a function of it. That is
circular: it shows the heuristic is internally consistent with
itself, not that Lucena/Philidor structure "emerged" from data.

Phase 3 removes that fallback for anything reported as a principle.
If a real Syzygy tablebase file for a material class is not present,
that class is skipped and reported as skipped — it is not filled in
with a heuristic pretending to be ground truth. Every endgame ranking
this script produces is backed by a genuine tablebase probe, and if
you run it with no tablebase files at all, the endgame stage will
correctly report zero results rather than manufacture some.

Phase 3 also adds a second, independent evidence source: real games
(a PGN you supply) cross-checked against a UCI engine (Stockfish or
similar) to flag positions where the played move was confirmed, or
contradicted, by engine search. That is the Tier-2 pipeline the
original build plan calls for. Its output is a correlation study of
which structural features predict engine agreement — useful, but
explicitly not proof of a causal "principle," and it inherits
whatever blind spots the chosen engine/depth has.

────────────────────────────────────────────────────────────────
SCOPE
────────────────────────────────────────────────────────────────
  TIER 1 / TIER 3 (endgame, tablebase-backed):
    A configurable set of material classes up to `--max-men` pieces
    (kings included). NOTE: "all 7-piece tablebase positions" as
    stated in the original plan is not something most machines can
    host — the full 7-man Syzygy set runs to tens of terabytes. This
    script defaults to max_men=6 (the 6-man WDL-only set is ~150GB
    and realistically downloadable) and will use whichever material
    classes' files you actually have on disk, up to whatever
    --max-men you set. Point --syzygy at a directory containing any
    mix of 3..7-man files and it will use what it finds.

  TIER 2 (middlegame, engine-confirmed):
    Positions sampled from a PGN file you supply, each checked
    against a UCI engine you supply, to determine whether the played
    move matched (or nearly matched, within --cp-margin) the
    engine's top choice at --engine-depth.

  PRINCIPLE HIERARCHY:
    Cross-references which structural features rank as important
    across *multiple, unrelated* material classes and also in the
    middlegame set. Those recurring features are the more
    interesting candidates for a general factor, as opposed to
    features that are only locally relevant to one endgame type.

  NAVIGATION FUNCTION:
    One-ply, principle-weighted move selection (no tree search),
    generalised to route to the right weight table by material
    class, validated against genuine before/after tablebase queries
    (never against the heuristic that produced the weights).

REQUIREMENTS:
    pip install python-chess numpy pandas scikit-learn matplotlib

Each stage is independently optional and is skipped with a clear
message if its inputs are missing:

    --syzygy PATH      directory of .rtbw/.rtbz tablebase files
    --pgn PATH         a PGN file of games to mine for Tier 2
    --stockfish PATH   a UCI-compatible engine binary

USAGE:
    python chess_topo_phase3.py \\
        --syzygy ./syzygy --max-men 6 \\
        --pgn games.pgn --stockfish /usr/local/bin/stockfish

    Outputs (written to --outdir, default ./phase3_output/):
      phase3_endgame_positions.csv
      phase3_middlegame_positions.csv
      phase3_importance_<CLASS>.csv        (one per resolved material class)
      phase3_importance_middlegame.csv
      phase3_principle_hierarchy.csv
      phase3_hierarchy.png
      phase3_validation_samples.csv
      phase3_results.json
"""

import argparse
import json
import os
import random
from collections import Counter

import numpy as np
import pandas as pd

import chess
import chess.pgn

try:
    import chess.engine
    ENGINE_AVAILABLE = True
except ImportError:
    ENGINE_AVAILABLE = False

try:
    import chess.syzygy
    SYZYGY_MODULE_AVAILABLE = True
except ImportError:
    SYZYGY_MODULE_AVAILABLE = False

try:
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import accuracy_score
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False
    print("[WARNING] scikit-learn not available. pip install scikit-learn\n")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    print("[WARNING] matplotlib not available. Visualisation skipped.\n")


# ══════════════════════════════════════════════════════════════
# SECTION 0 — CONFIG & MATERIAL CLASSES
# ══════════════════════════════════════════════════════════════

DEFAULT_CONFIG = {
    "syzygy_path":                   "./syzygy",
    "max_men":                       6,
    "pgn_path":                      None,
    "stockfish_path":                None,
    "engine_depth":                  16,
    "cp_agreement_margin":           25,
    "positions_per_material_class":  3000,
    "max_games":                     500,
    "max_positions_per_game":        10,
    "skip_opening_plies":            16,
    "random_seed":                   20260823,
    "outdir":                        "./phase3_output",
}

# (white_extra_pieces, black_extra_pieces), NOT including kings.
# total_men = 2 + len(white_extra) + len(black_extra)
MATERIAL_CLASSES = {
    # -- 4/5-man: cheap sanity checks, small tablebase files --
    "KQKR":    ([chess.QUEEN], [chess.ROOK]),
    "KRKB":    ([chess.ROOK], [chess.BISHOP]),
    "KRBKR":   ([chess.ROOK, chess.BISHOP], [chess.ROOK]),
    "KRNKR":   ([chess.ROOK, chess.KNIGHT], [chess.ROOK]),
    "KQKRR":   ([chess.QUEEN], [chess.ROOK, chess.ROOK]),
    "KBBKN":   ([chess.BISHOP, chess.BISHOP], [chess.KNIGHT]),
    # -- 6-man --
    "KRPKR":   ([chess.ROOK, chess.PAWN], [chess.ROOK]),
    "KQPKQ":   ([chess.QUEEN, chess.PAWN], [chess.QUEEN]),
    "KRRKQ":   ([chess.ROOK, chess.ROOK], [chess.QUEEN]),
    "KRBKRN":  ([chess.ROOK, chess.BISHOP], [chess.ROOK, chess.KNIGHT]),
    # -- 7-man: only attempted if --max-men 7 AND the files exist --
    "KQRKQRP": ([chess.QUEEN, chess.ROOK], [chess.QUEEN, chess.ROOK, chess.PAWN]),
    "KRRPKRR": ([chess.ROOK, chess.ROOK, chess.PAWN], [chess.ROOK, chess.ROOK]),
}

PIECE_ORDER = [chess.QUEEN, chess.ROOK, chess.BISHOP, chess.KNIGHT, chess.PAWN]


def build_material_lookup(material_classes):
    lookup = {}
    for name, (w_extra, b_extra) in material_classes.items():
        key = (tuple(sorted(w_extra)), tuple(sorted(b_extra)))
        lookup[key] = name
    return lookup


def material_class_of(board, lookup):
    w_types, b_types = [], []
    for pt in PIECE_ORDER:
        w_types += [pt] * len(board.pieces(pt, chess.WHITE))
        b_types += [pt] * len(board.pieces(pt, chess.BLACK))
    key = (tuple(sorted(w_types)), tuple(sorted(b_types)))
    return lookup.get(key, "UNCLASSIFIED")


# ══════════════════════════════════════════════════════════════
# SECTION 1 — TABLEBASE ORACLE (REAL GROUND TRUTH ONLY, NO FALLBACK)
# ══════════════════════════════════════════════════════════════

class SyzygyOracle:
    """
    Thin wrapper around chess.syzygy. Deliberately has NO heuristic
    fallback: if a position can't be probed, `query` returns None,
    and callers must treat that as "unknown" rather than inventing
    an answer.
    """

    def __init__(self, path):
        self.path = path
        self.tb = None
        if not SYZYGY_MODULE_AVAILABLE:
            print("[INFO] chess.syzygy not available in this "
                  "python-chess install.")
            return
        if path and os.path.isdir(path):
            try:
                self.tb = chess.syzygy.open_tablebase(path)
                print(f"[OK] Syzygy tablebase opened at {path}")
            except Exception as e:
                print(f"[INFO] Could not open tablebase at {path}: {e}")
        else:
            print(f"[INFO] No tablebase directory at '{path}'. "
                  f"Endgame (Tier 1/3) mining will be skipped for any "
                  f"material class that can't be resolved. Download "
                  f"files from https://syzygy-tables.info/ and pass "
                  f"--syzygy PATH to enable this stage.")

    def query(self, board):
        """Returns 'WIN' / 'DRAW' / 'LOSS' (side-to-move perspective), or None."""
        if self.tb is None:
            return None
        try:
            wdl = self.tb.probe_wdl(board)
        except Exception:
            return None
        if wdl > 0:
            return "WIN"
        if wdl == 0:
            return "DRAW"
        return "LOSS"

    def close(self):
        if self.tb is not None:
            try:
                self.tb.close()
            except Exception:
                pass


# ══════════════════════════════════════════════════════════════
# SECTION 2 — GENERIC RANDOM POSITION GENERATOR
# ══════════════════════════════════════════════════════════════

def random_position(white_extra, black_extra, rng, max_tries=300):
    """
    Build a random legal position with two kings plus the given extra
    piece lists for each side. Returns a chess.Board, or None if no
    legal placement was found within max_tries attempts.
    """
    squares = list(chess.SQUARES)
    n_needed = 2 + len(white_extra) + len(black_extra)

    for _ in range(max_tries):
        chosen = rng.sample(squares, n_needed)
        wk_sq, bk_sq = chosen[0], chosen[1]
        rest = chosen[2:]

        if chess.square_distance(wk_sq, bk_sq) <= 1:
            continue

        board = chess.Board(fen=None)
        board.clear()
        board.set_piece_at(wk_sq, chess.Piece(chess.KING, chess.WHITE))
        board.set_piece_at(bk_sq, chess.Piece(chess.KING, chess.BLACK))

        ok, idx = True, 0
        for pt in white_extra:
            sq = rest[idx]
            idx += 1
            if pt == chess.PAWN and chess.square_rank(sq) in (0, 7):
                ok = False
                break
            board.set_piece_at(sq, chess.Piece(pt, chess.WHITE))
        if not ok:
            continue

        for pt in black_extra:
            sq = rest[idx]
            idx += 1
            if pt == chess.PAWN and chess.square_rank(sq) in (0, 7):
                ok = False
                break
            board.set_piece_at(sq, chess.Piece(pt, chess.BLACK))
        if not ok:
            continue

        board.turn = rng.choice([chess.WHITE, chess.BLACK])

        if not board.is_valid():
            continue

        return board

    return None


# ══════════════════════════════════════════════════════════════
# SECTION 3 — GENERIC STRUCTURAL FEATURE EXTRACTOR
# ══════════════════════════════════════════════════════════════
# Works on ANY position, endgame or middlegame. Unlike Phase 1/2,
# which hard-coded a separate feature set per material class, this
# extractor is material-agnostic so the same pipeline covers KQKR,
# KRBKRN, and full middlegame positions from real games alike.

CENTER_SQUARES_FR = [(3, 3), (3, 4), (4, 3), (4, 4)]  # d4, d5, e4, e5


def _center_distance(sq):
    f, r = chess.square_file(sq), chess.square_rank(sq)
    return min(abs(f - cf) + abs(r - cr) for cf, cr in CENTER_SQUARES_FR)


def _count_passed_pawns(board, color):
    pawns = board.pieces(chess.PAWN, color)
    opp_pawns = board.pieces(chess.PAWN, not color)
    count = 0
    for sq in pawns:
        f, r = chess.square_file(sq), chess.square_rank(sq)
        passed = True
        for osq in opp_pawns:
            of, orr = chess.square_file(osq), chess.square_rank(osq)
            if abs(of - f) <= 1:
                if color == chess.WHITE and orr > r:
                    passed = False
                    break
                if color == chess.BLACK and orr < r:
                    passed = False
                    break
        if passed:
            count += 1
    return count


def _count_isolated_pawns(pawn_squares):
    files = set(chess.square_file(s) for s in pawn_squares)
    return sum(
        1 for s in pawn_squares
        if (chess.square_file(s) - 1) not in files
        and (chess.square_file(s) + 1) not in files
    )


def _count_doubled_pawns(pawn_squares):
    c = Counter(chess.square_file(s) for s in pawn_squares)
    return sum(v - 1 for v in c.values() if v > 1)


def _rooks_on_open_files(board, color):
    rooks = board.pieces(chess.ROOK, color)
    pawn_files = set(chess.square_file(s) for s in board.pieces(chess.PAWN, chess.WHITE))
    pawn_files |= set(chess.square_file(s) for s in board.pieces(chess.PAWN, chess.BLACK))
    return sum(1 for r in rooks if chess.square_file(r) not in pawn_files)


def _pawn_shield_count(board, color):
    k = board.king(color)
    if k is None:
        return 0
    kf, kr = chess.square_file(k), chess.square_rank(k)
    shield_rank = kr + 1 if color == chess.WHITE else kr - 1
    if not (0 <= shield_rank <= 7):
        return 0
    count = 0
    for df in (-1, 0, 1):
        f = kf + df
        if 0 <= f <= 7:
            p = board.piece_at(chess.square(f, shield_rank))
            if p and p.piece_type == chess.PAWN and p.color == color:
                count += 1
    return count


def _center_control(board, color):
    squares = [chess.D4, chess.D5, chess.E4, chess.E5]
    return sum(1 for sq in squares if board.is_attacked_by(color, sq))


def _attackers_near_king(board, attacking_color, king_sq):
    if king_sq is None:
        return 0
    count = 0
    for sq in chess.SQUARES:
        if chess.square_distance(sq, king_sq) <= 2:
            p = board.piece_at(sq)
            if p and p.color == attacking_color and p.piece_type != chess.KING:
                count += 1
    return count


PIECE_VALUES = {
    chess.PAWN: 1, chess.KNIGHT: 3, chess.BISHOP: 3,
    chess.ROOK: 5, chess.QUEEN: 9,
}


def extract_generic_features(board):
    f = {}
    f["white_to_move"] = int(board.turn == chess.WHITE)
    f["fullmove_number"] = board.fullmove_number
    f["halfmove_clock"] = board.halfmove_clock

    total_material = 0
    for pt, val in PIECE_VALUES.items():
        sym = chess.piece_symbol(pt)
        wc = len(board.pieces(pt, chess.WHITE))
        bc = len(board.pieces(pt, chess.BLACK))
        f[f"w_{sym}"] = wc
        f[f"b_{sym}"] = bc
        f[f"diff_{sym}"] = wc - bc
        total_material += (wc + bc) * val

    f["material_balance"] = sum(
        PIECE_VALUES[pt] * (
            len(board.pieces(pt, chess.WHITE)) - len(board.pieces(pt, chess.BLACK))
        )
        for pt in PIECE_VALUES
    )
    f["total_material"] = total_material
    f["bishop_pair_w"] = int(len(board.pieces(chess.BISHOP, chess.WHITE)) >= 2)
    f["bishop_pair_b"] = int(len(board.pieces(chess.BISHOP, chess.BLACK)) >= 2)

    wk, bk = board.king(chess.WHITE), board.king(chess.BLACK)
    if wk is not None and bk is not None:
        f["king_distance"] = chess.square_distance(wk, bk)
        f["wk_file"], f["wk_rank"] = chess.square_file(wk), chess.square_rank(wk)
        f["bk_file"], f["bk_rank"] = chess.square_file(bk), chess.square_rank(bk)
        f["wk_center_dist"] = _center_distance(wk)
        f["bk_center_dist"] = _center_distance(bk)
    else:
        f["king_distance"] = 0
        f["wk_file"] = f["wk_rank"] = f["bk_file"] = f["bk_rank"] = 0
        f["wk_center_dist"] = f["bk_center_dist"] = 0

    f["mobility_stm"] = len(list(board.legal_moves))
    f["mobility_other"] = 0
    if not board.is_check():
        try:
            b2 = board.copy(stack=False)
            b2.push(chess.Move.null())
            f["mobility_other"] = len(list(b2.legal_moves))
        except Exception:
            pass

    f["in_check"] = int(board.is_check())

    w_pawns = list(board.pieces(chess.PAWN, chess.WHITE))
    b_pawns = list(board.pieces(chess.PAWN, chess.BLACK))
    f["w_passed_pawns"] = _count_passed_pawns(board, chess.WHITE)
    f["b_passed_pawns"] = _count_passed_pawns(board, chess.BLACK)
    f["w_isolated_pawns"] = _count_isolated_pawns(w_pawns)
    f["b_isolated_pawns"] = _count_isolated_pawns(b_pawns)
    f["w_doubled_pawns"] = _count_doubled_pawns(w_pawns)
    f["b_doubled_pawns"] = _count_doubled_pawns(b_pawns)
    f["w_most_advanced_pawn"] = max((chess.square_rank(s) for s in w_pawns), default=0)
    f["b_most_advanced_pawn"] = min((chess.square_rank(s) for s in b_pawns), default=7)

    f["w_rooks_open_file"] = _rooks_on_open_files(board, chess.WHITE)
    f["b_rooks_open_file"] = _rooks_on_open_files(board, chess.BLACK)
    f["w_rooks_7th"] = sum(1 for s in board.pieces(chess.ROOK, chess.WHITE) if chess.square_rank(s) == 6)
    f["b_rooks_2nd"] = sum(1 for s in board.pieces(chess.ROOK, chess.BLACK) if chess.square_rank(s) == 1)

    f["wk_pawn_shield"] = _pawn_shield_count(board, chess.WHITE)
    f["bk_pawn_shield"] = _pawn_shield_count(board, chess.BLACK)

    f["w_center_control"] = _center_control(board, chess.WHITE)
    f["b_center_control"] = _center_control(board, chess.BLACK)

    f["w_attackers_near_bk"] = _attackers_near_king(board, chess.WHITE, bk)
    f["b_attackers_near_wk"] = _attackers_near_king(board, chess.BLACK, wk)

    return f


# ══════════════════════════════════════════════════════════════
# SECTION 4 — TIER 1/3: ENDGAME PRINCIPLE MINING (TABLEBASE-BACKED)
# ══════════════════════════════════════════════════════════════

def compute_feature_importance(df, label_col, drop_cols, seed):
    X = df.drop(columns=[c for c in drop_cols if c in df.columns])
    X = X.select_dtypes(include=[np.number]).fillna(0)
    y = df[label_col]

    if y.nunique() < 2:
        return None, None

    strat = y if y.value_counts().min() >= 2 else None
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=seed, stratify=strat
    )
    clf = RandomForestClassifier(n_estimators=300, random_state=seed, n_jobs=-1)
    clf.fit(X_train, y_train)
    acc = accuracy_score(y_test, clf.predict(X_test))
    importances = pd.Series(
        clf.feature_importances_, index=X.columns
    ).sort_values(ascending=False)
    return importances, acc


def mine_endgame_principles(material_classes, oracle, rng, cfg):
    print("=" * 64)
    print("STEP: TIER 1/3 ENDGAME MINING (tablebase-backed only)")
    print("=" * 64)

    if oracle.tb is None:
        print("  [SKIP] No tablebase loaded -- nothing to mine here.\n")
        return pd.DataFrame(), {}

    all_frames = []
    per_class_importance = {}

    for name, (w_extra, b_extra) in material_classes.items():
        total_men = 2 + len(w_extra) + len(b_extra)
        if total_men > cfg["max_men"]:
            print(f"  [SKIP] {name}: {total_men}-man class exceeds "
                  f"--max-men {cfg['max_men']}")
            continue

        rows = []
        attempts = 0
        target = cfg["positions_per_material_class"]
        probe_failures = 0

        while len(rows) < target and attempts < target * 25:
            attempts += 1
            board = random_position(w_extra, b_extra, rng)
            if board is None:
                continue
            outcome = oracle.query(board)
            if outcome is None:
                probe_failures += 1
                if probe_failures > 200 and not rows:
                    break  # this class's tablebase files are almost certainly absent
                continue
            feats = extract_generic_features(board)
            feats["material_class"] = name
            feats["outcome"] = outcome
            feats["fen"] = board.fen()
            rows.append(feats)

        if not rows:
            print(f"  [SKIP] {name}: 0 positions resolved "
                  f"(tablebase files for this class not found in "
                  f"'{cfg['syzygy_path']}')")
            continue

        df = pd.DataFrame(rows)
        counts = df["outcome"].value_counts().to_dict()
        print(f"  [OK] {name}: {len(df)} positions {counts}")
        all_frames.append(df)

        if SKLEARN_AVAILABLE:
            imp, acc = compute_feature_importance(
                df, "outcome",
                drop_cols=["material_class", "outcome", "fen"],
                seed=cfg["random_seed"],
            )
            if imp is not None:
                per_class_importance[name] = imp
                print(f"        held-out accuracy: {acc:.3f} | "
                      f"top feature: {imp.index[0]} ({imp.iloc[0]:.3f})")
            else:
                print("        only one outcome class present -- "
                      "skipping importance ranking")

    combined = pd.concat(all_frames, ignore_index=True) if all_frames else pd.DataFrame()
    print()
    return combined, per_class_importance


# ══════════════════════════════════════════════════════════════
# SECTION 5 — TIER 2: ENGINE-CONFIRMED MIDDLEGAME MINING
# ══════════════════════════════════════════════════════════════

def _score_to_cp(pov_score_obj, color):
    if pov_score_obj is None:
        return None
    try:
        return pov_score_obj.pov(color).score(mate_score=100000)
    except Exception:
        return None


def mine_middlegame_principles(cfg):
    print("=" * 64)
    print("STEP: TIER 2 MIDDLEGAME MINING (engine-confirmed)")
    print("=" * 64)

    if not cfg["pgn_path"] or not os.path.exists(cfg["pgn_path"]):
        print("  [SKIP] No PGN file supplied (--pgn PATH).\n")
        return pd.DataFrame(), None

    if not ENGINE_AVAILABLE:
        print("  [SKIP] chess.engine not available in this "
              "python-chess install.\n")
        return pd.DataFrame(), None

    if not cfg["stockfish_path"] or not os.path.exists(cfg["stockfish_path"]):
        print("  [SKIP] No UCI engine supplied (--stockfish PATH) -- "
              "cannot confirm move optimality without one.\n")
        return pd.DataFrame(), None

    print(f"  Starting engine: {cfg['stockfish_path']}")
    try:
        engine = chess.engine.SimpleEngine.popen_uci(cfg["stockfish_path"])
        print(f"  Engine started: {engine.id['name']}\n")
    except Exception as e:
        print(f"  [ERROR] Could not start engine: {e}\n")
        import traceback
        traceback.print_exc()
        return pd.DataFrame(), None

    rows = []
    games_used = 0
    limit = chess.engine.Limit(depth=cfg["engine_depth"])

    try:
        with open(cfg["pgn_path"], encoding="utf-8", errors="ignore") as pgn_file:
            while games_used < cfg["max_games"]:
                game = chess.pgn.read_game(pgn_file)
                if game is None:
                    break

                board = game.board()
                ply, sampled = 0, 0

                for move in game.mainline_moves():
                    ply += 1
                    if board.is_game_over():
                        break
                    if ply < cfg["skip_opening_plies"] or sampled >= cfg["max_positions_per_game"]:
                        board.push(move)
                        continue

                    mover = board.turn
                    try:
                        info = engine.analyse(board, limit)
                    except Exception as e:
                        print(f"    [WARNING] Engine analysis failed: {e}")
                        board.push(move)
                        continue

                    pv = info.get("pv")
                    best_move = pv[0] if pv else None
                    best_cp = _score_to_cp(info.get("score"), mover)

                    confirmed = None
                    if best_move is not None and best_cp is not None:
                        if move == best_move:
                            confirmed = True
                        else:
                            b2 = board.copy()
                            b2.push(move)
                            try:
                                info2 = engine.analyse(b2, limit)
                                played_cp = _score_to_cp(info2.get("score"), mover)
                            except Exception:
                                played_cp = None
                            if played_cp is not None:
                                confirmed = (best_cp - played_cp) <= cfg["cp_agreement_margin"]

                    if confirmed is not None:
                        feats = extract_generic_features(board)
                        feats["confirmed_optimal"] = int(confirmed)
                        feats["fen"] = board.fen()
                        feats["played_move"] = board.san(move)
                        rows.append(feats)
                        sampled += 1

                    board.push(move)

                games_used += 1
                if games_used % 10 == 0:
                    print(f"  ... {games_used} games processed, "
                          f"{len(rows)} labelled positions so far")

    finally:
        engine.quit()

    df = pd.DataFrame(rows)
    print(f"  [OK] {len(df)} labelled positions from {games_used} games\n")

    importance = None
    if SKLEARN_AVAILABLE and not df.empty:
        imp, acc = compute_feature_importance(
            df, "confirmed_optimal",
            drop_cols=["fen", "played_move", "confirmed_optimal"],
            seed=cfg["random_seed"],
        )
        if imp is not None:
            importance = imp
            print(f"  held-out accuracy predicting engine agreement: {acc:.3f}")
            print(f"  top correlate: {imp.index[0]} ({imp.iloc[0]:.3f})\n")

    return df, importance


# ══════════════════════════════════════════════════════════════
# SECTION 6 — PRINCIPLE HIERARCHY / CROSS-CLASS ANALYSIS
# ══════════════════════════════════════════════════════════════

def analyse_principle_hierarchy(per_class_importance, middlegame_importance, top_n=8):
    """
    A feature that ranks in the top-N for several *unrelated* material
    classes, and also in the middlegame set, is a candidate for a
    genuinely general structural factor rather than something local to
    one endgame type. This is exploratory pattern-matching, not a
    proof of a causal hierarchy -- treat the output as a shortlist to
    investigate further, e.g. by checking against chess literature or
    testing whether navigation accuracy drops when the feature is
    removed.
    """
    records = []
    for cls, imp in per_class_importance.items():
        for rank, (feat, val) in enumerate(imp.head(top_n).items(), start=1):
            records.append({"source": cls, "feature": feat, "rank": rank, "importance": val})
    if middlegame_importance is not None:
        for rank, (feat, val) in enumerate(middlegame_importance.head(top_n).items(), start=1):
            records.append({"source": "MIDDLEGAME", "feature": feat, "rank": rank, "importance": val})

    if not records:
        return pd.DataFrame(), pd.DataFrame()

    rank_df = pd.DataFrame(records)
    cross_class = (
        rank_df.groupby("feature")
        .agg(n_sources=("source", "nunique"),
             mean_rank=("rank", "mean"),
             mean_importance=("importance", "mean"))
        .sort_values(["n_sources", "mean_importance"], ascending=[False, False])
    )
    return rank_df, cross_class


# ══════════════════════════════════════════════════════════════
# SECTION 7 — NAVIGATION FUNCTION (ONE-PLY, PRINCIPLE-WEIGHTED)
# ══════════════════════════════════════════════════════════════

def navigate(board, weight_lookup, material_lookup, fallback_weights=None):
    """
    No tree search. For each legal move, score the resulting position
    as a weighted sum of its structural features, using the importance
    weights derived for that position's material class (falling back
    to the middlegame weights, then to plain material balance).
    """
    legal = list(board.legal_moves)
    if not legal:
        return None, "no legal moves", 0.0

    cls = material_class_of(board, material_lookup)
    weights = weight_lookup.get(cls)
    source = cls
    if weights is None:
        weights = fallback_weights
        source = "MIDDLEGAME/fallback"
    if weights is None:
        weights = pd.Series({"material_balance": 1.0})
        source = "material-only fallback"

    mover_is_white = board.turn == chess.WHITE
    best_move, best_score = None, -float("inf")

    for mv in legal:
        b2 = board.copy()
        b2.push(mv)
        feats = extract_generic_features(b2)
        score = 0.0
        for feat, w in weights.items():
            score += w * feats.get(feat, 0) * (1 if mover_is_white else -1)
        if score > best_score:
            best_score, best_move = score, mv

    return best_move, source, best_score


def judge_wdl_after_move(before, after):
    """
    `before`/`after` are WDL strings from the perspective of whoever
    was to move at that point. A move flips the side to move, so
    `after` is from the opponent's perspective -- flip it back to the
    mover's perspective before comparing.
    """
    flip = {"WIN": "LOSS", "DRAW": "DRAW", "LOSS": "WIN"}
    after_movers_pov = flip.get(after, after)
    order = {"LOSS": 0, "DRAW": 1, "WIN": 2}
    return order.get(after_movers_pov, 1) >= order.get(before, 1)


def validate_navigation(material_classes, per_class_importance, oracle,
                         material_lookup, rng, cfg, n_tests=300):
    print("=" * 64)
    print("STEP: NAVIGATION VALIDATION (against genuine tablebase probes)")
    print("=" * 64)

    if not per_class_importance:
        print("  [SKIP] No endgame principle library to validate.\n")
        return 0.0, pd.DataFrame()

    per_class_n = max(1, n_tests // max(1, len(per_class_importance)))
    correct, total, samples = 0, 0, []

    for name in per_class_importance:
        w_extra, b_extra = material_classes[name]
        tries, done = 0, 0
        while done < per_class_n and tries < per_class_n * 30:
            tries += 1
            board = random_position(w_extra, b_extra, rng)
            if board is None:
                continue
            before = oracle.query(board)
            if before is None:
                continue
            mv, _, _ = navigate(board, per_class_importance, material_lookup)
            if mv is None:
                continue
            san = board.san(mv)
            b2 = board.copy()
            b2.push(mv)
            after = oracle.query(b2)
            if after is None:
                continue
            ok = judge_wdl_after_move(before, after)
            correct += int(ok)
            total += 1
            done += 1
            samples.append({
                "material_class": name, "fen": board.fen(), "move": san,
                "oracle_before": before, "oracle_after": after, "ok": ok,
            })

    acc = correct / total if total else 0.0
    print(f"  {correct}/{total} moves preserved-or-improved the "
          f"tablebase-verified outcome ({acc:.1%})\n")
    return acc, pd.DataFrame(samples)


# ══════════════════════════════════════════════════════════════
# SECTION 8 — REPORTING
# ══════════════════════════════════════════════════════════════

def visualise_hierarchy(cross_class_df, output_path, top_n=20):
    if not MATPLOTLIB_AVAILABLE or cross_class_df.empty:
        return
    top = cross_class_df.head(top_n)
    fig, ax = plt.subplots(figsize=(9, max(4, 0.35 * len(top))))
    ax.barh(top.index[::-1], top["n_sources"][::-1])
    ax.set_xlabel("Number of distinct sources this feature ranked top-N in")
    ax.set_title("Cross-class feature recurrence (Phase 3)")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {output_path}")


def print_methodological_notes():
    print("""
────────────────────────────────────────────────────────────────
METHODOLOGICAL NOTES -- READ BEFORE CALLING ANYTHING A "PRINCIPLE"
────────────────────────────────────────────────────────────────
  1. Every endgame ranking above came from a REAL Syzygy probe.
     Classes with 0 positions mean the tablebase files for that
     material signature were not found on disk -- not that anything
     was proven or disproven about that endgame.

  2. Middlegame ("Tier 2") rankings measure correlation with engine
     agreement at a fixed search depth. A high-ranking feature is a
     candidate worth investigating, not a confirmed causal principle.
     Change the engine, the depth, or the game sample and rankings
     may shift -- that sensitivity is itself useful information, not
     a flaw to explain away.

  3. Endgame positions were sampled uniformly at random from legal
     placements, which is NOT the distribution of positions that
     occur in actual play. A feature can distinguish WIN from DRAW
     across random placements without being the feature a player
     most needs in realistic positions.

  4. The most interesting output here is the cross-class table:
     features that rank highly in several unrelated material classes
     AND in the middlegame set. Those are worth checking against
     existing chess literature, and worth testing for necessity (does
     navigation accuracy drop if you remove them?) before treating
     them as new discoveries.
────────────────────────────────────────────────────────────────
""")


# ══════════════════════════════════════════════════════════════
# SECTION 9 — MAIN
# ══════════════════════════════════════════════════════════════

def build_arg_parser():
    p = argparse.ArgumentParser(description="Chess Topological Solution — Phase 3")
    p.add_argument("--syzygy", default=DEFAULT_CONFIG["syzygy_path"])
    p.add_argument("--max-men", type=int, default=DEFAULT_CONFIG["max_men"])
    p.add_argument("--pgn", default=DEFAULT_CONFIG["pgn_path"])
    p.add_argument("--stockfish", default=DEFAULT_CONFIG["stockfish_path"])
    p.add_argument("--engine-depth", type=int, default=DEFAULT_CONFIG["engine_depth"])
    p.add_argument("--cp-margin", type=int, default=DEFAULT_CONFIG["cp_agreement_margin"])
    p.add_argument("--positions-per-class", type=int,
                    default=DEFAULT_CONFIG["positions_per_material_class"])
    p.add_argument("--max-games", type=int, default=DEFAULT_CONFIG["max_games"])
    p.add_argument("--max-positions-per-game", type=int,
                    default=DEFAULT_CONFIG["max_positions_per_game"])
    p.add_argument("--skip-opening-plies", type=int,
                    default=DEFAULT_CONFIG["skip_opening_plies"])
    p.add_argument("--seed", type=int, default=DEFAULT_CONFIG["random_seed"])
    p.add_argument("--outdir", default=DEFAULT_CONFIG["outdir"])
    return p


def main():
    args = build_arg_parser().parse_args()
    cfg = dict(DEFAULT_CONFIG)
    cfg.update({
        "syzygy_path":                  args.syzygy,
        "max_men":                      args.max_men,
        "pgn_path":                     args.pgn,
        "stockfish_path":               args.stockfish,
        "engine_depth":                 args.engine_depth,
        "cp_agreement_margin":          args.cp_margin,
        "positions_per_material_class": args.positions_per_class,
        "max_games":                    args.max_games,
        "max_positions_per_game":       args.max_positions_per_game,
        "skip_opening_plies":           args.skip_opening_plies,
        "random_seed":                  args.seed,
        "outdir":                       args.outdir,
    })

    print("=" * 64)
    print("CHESS TOPOLOGICAL SOLUTION — PHASE 3")
    print("Multi-piece endgames (tablebase) + engine-confirmed middlegame")
    print("OrganismCore | Eric Robert Lawson | 2026-08-23")
    print("=" * 64)
    print()

    rng = random.Random(cfg["random_seed"])
    os.makedirs(cfg["outdir"], exist_ok=True)
    material_lookup = build_material_lookup(MATERIAL_CLASSES)

    # -- TIER 1/3: endgame --
    oracle = SyzygyOracle(cfg["syzygy_path"])
    endgame_df, per_class_importance = mine_endgame_principles(
        MATERIAL_CLASSES, oracle, rng, cfg
    )

    # -- TIER 2: middlegame --
    middlegame_df, middlegame_importance = mine_middlegame_principles(cfg)

    # -- HIERARCHY --
    print("=" * 64)
    print("STEP: PRINCIPLE HIERARCHY / CROSS-CLASS ANALYSIS")
    print("=" * 64)
    rank_df, cross_class_df = analyse_principle_hierarchy(
        per_class_importance, middlegame_importance
    )
    if cross_class_df.empty:
        print("  [SKIP] Not enough importance tables to cross-reference "
              "(need at least one resolved endgame class or the "
              "middlegame set).\n")
    else:
        print(cross_class_df.head(15).to_string())
        print()

    # -- SAVE ARTIFACTS --
    if not endgame_df.empty:
        endgame_df.to_csv(os.path.join(cfg["outdir"], "phase3_endgame_positions.csv"), index=False)
    if not middlegame_df.empty:
        middlegame_df.to_csv(os.path.join(cfg["outdir"], "phase3_middlegame_positions.csv"), index=False)
    for cls, imp in per_class_importance.items():
        imp.to_csv(os.path.join(cfg["outdir"], f"phase3_importance_{cls}.csv"), header=["importance"])
    if middlegame_importance is not None:
        middlegame_importance.to_csv(
            os.path.join(cfg["outdir"], "phase3_importance_middlegame.csv"), header=["importance"]
        )
    if not cross_class_df.empty:
        cross_class_df.to_csv(os.path.join(cfg["outdir"], "phase3_principle_hierarchy.csv"))
        visualise_hierarchy(cross_class_df, os.path.join(cfg["outdir"], "phase3_hierarchy.png"))

    # -- VALIDATE --
    acc, val_samples = validate_navigation(
        MATERIAL_CLASSES, per_class_importance, oracle, material_lookup, rng, cfg
    )
    if not val_samples.empty:
        val_samples.to_csv(os.path.join(cfg["outdir"], "phase3_validation_samples.csv"), index=False)

    # -- RESULTS SUMMARY --
    results = {
        "phase": 3,
        "date": "2026-08-23",
        "endgame_classes_resolved": list(per_class_importance.keys()),
        "endgame_positions": int(len(endgame_df)),
        "middlegame_positions": int(len(middlegame_df)),
        "cross_class_features": (
            cross_class_df.reset_index().to_dict(orient="records")
            if not cross_class_df.empty else []
        ),
        "navigation_accuracy": float(acc),
    }
    with open(os.path.join(cfg["outdir"], "phase3_results.json"), "w") as fh:
        json.dump(results, fh, indent=2, default=str)

    oracle.close()

    print("=" * 64)
    print("PHASE 3 COMPLETE — SUMMARY")
    print("=" * 64)
    print(f"  Endgame classes resolved: {len(per_class_importance)} "
          f"of {len(MATERIAL_CLASSES)} configured")
    print(f"  Endgame positions:        {len(endgame_df):,}")
    print(f"  Middlegame positions:     {len(middlegame_df):,}")
    print(f"  Navigation accuracy:      {acc:.1%}")
    print(f"  Output directory:         {cfg['outdir']}")
    print_methodological_notes()


if __name__ == "__main__":
    main()
