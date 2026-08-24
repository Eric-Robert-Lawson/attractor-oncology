============================================================
CHESS TOPOLOGICAL SOLUTION — PHASE 2
Rook Endgames: KRK, KRPKR, KRKP
OrganismCore | Eric Robert Lawson | 2026-08-23
============================================================

PURPOSE:
  Derive Lucena and Philidor positions as
  topological principles from structural analysis
  of position data — without hand-coding,
  without tree search, without statistics.

PHASE 1 CARRIED FORWARD:
  Navigation accuracy: 93.0%
  Methodology: VALIDATED
  Key gaps: Key Square refinement,
            Zugzwang recognition

[INFO] No syzygy/ directory found.
[INFO] Using heuristic rook endgame evaluation.
[INFO] Download tablebases from https://syzygy-tables.info/ for full accuracy.

============================================================
STEP 1: GENERATING POSITIONS
============================================================

──────────────────────────────────────────────────
Generating KRK positions...
  Generated: 8,000 | Skipped: 2,616
  WIN: 8,000 | DRAW: 0

──────────────────────────────────────────────────
Generating KRPKR positions...
  Generated 5,000 KRPKR positions...
  Generated 10,000 KRPKR positions...
  Generated 15,000 KRPKR positions...
  Generated: 15,000 | Skipped: 5,057
  WIN: 9,881 | DRAW: 5,119 | LOSS: 0

──────────────────────────────────────────────────
Generating KRKP positions...
  Generated: 5,000 | Skipped: 1,552
  WIN: 4,744 | DRAW: 256

  Total positions: 28,000

============================================================
STEP 2: TOPOLOGICAL PRINCIPLE EXTRACTION
============================================================

=======================================================
PRINCIPLE EXTRACTION: KRK
=======================================================
  Positions: 8,000
  Features:  18
  WIN: 8,000 (100.0%)
  Only one outcome class — skipping RF.
  KRK dataset saved: krk_positions.csv

=======================================================
PRINCIPLE EXTRACTION: KRPKR
=======================================================
  Positions: 15,000
  Features:  36
  WIN: 9,881 (65.9%)
  DRAW: 5,119 (34.1%)

  Model validation:
              precision    recall  f1-score   support

     NON-WIN      1.000     0.972     0.986      1024
         WIN      0.986     1.000     0.993      1976

    accuracy                          0.990      3000
   macro avg      0.993     0.986     0.989      3000
weighted avg      0.990     0.990     0.990      3000


  Rank  Feature                            Imp      WIN  NON-WIN  Name
  ────────────────────────────────────────────────────────────────────────
  1     pawn_rank                       0.4127    4.491    1.570  T1[KNOWN] Pawn Advancement
  2     pawn_steps_to_promo             0.3621    2.509    5.430  T3[NEW?]  UNNAMED #2
  3     lucena_score                    0.0747    3.848    2.164  T1[KNOWN] Lucena Position
  4     philidor_score                  0.0372    0.983    1.682  T1[KNOWN] Philidor Position
  5     wk_rank_advantage               0.0360   -0.945    1.988  T3[NEW?]  UNNAMED #5
  6     wk_leading                      0.0173    0.327    0.691  T3[NEW?]  UNNAMED #6
  7     bk_cut_off_distance             0.0088    2.708    2.673  T3[NEW?]  UNNAMED #7
  8     pawn_on_7th                     0.0081    0.242    0.000  T3[NEW?]  UNNAMED #8
  9     br_cuts_wk_from_pawn            0.0061    0.048    0.148  T3[NEW?]  UNNAMED #9
  10    wk_to_promo                     0.0053    4.308    4.290  T3[NEW?]  UNNAMED #10
  11    wk_shelters_pawn                0.0052    0.148    0.264  T1[KNOWN] King Shelter
  12    wr_to_bridge                    0.0040    4.357    3.625  T3[NEW?]  UNNAMED #12
  13    bk_to_pawn                      0.0036    3.519    3.701  T3[NEW?]  UNNAMED #13
  14    bk_to_promo                     0.0035    4.405    4.243  T3[NEW?]  UNNAMED #14
  15    wk_to_pawn                      0.0034    3.457    3.636  T3[NEW?]  UNNAMED #15
  16    br_on_6th                       0.0024    0.113    0.148  T1[KNOWN] Philidor Rook Cut-off
  17    promo_race                      0.0022   -0.097    0.047  T3[NEW?]  UNNAMED #17
  18    wk_pawn_file_dist               0.0021    2.647    2.662  T3[NEW?]  UNNAMED #18
  19    wr_to_promo                     0.0018    4.357    4.317  T3[NEW?]  UNNAMED #19
  20    kings_distance                  0.0009    4.074    4.061  T3[NEW?]  UNNAMED #20

  KNOWN PRINCIPLE VALIDATION:
  ✓ CONFIRMED  Rank # 3  Imp=0.0747  Lucena Position
  ✓ CONFIRMED  Rank # 4  Imp=0.0372  Philidor Position
  ✓ CONFIRMED  Rank # 1  Imp=0.4127  Pawn Advancement
  ~ WEAK       Rank #11  Imp=0.0052  King Shelter
  ~ WEAK       Rank #16  Imp=0.0024  Philidor Rook Cut-off
  ✗ NOT FOUND  Rook on 7th Rank
  ✗ NOT FOUND  Rook Pawn Exception

  Confirmed: 3/7

  TIER 3 CANDIDATES:

    Feature:  pawn_steps_to_promo
    Imp:      0.3621
    WIN mean: 2.509
    NON-WIN:  5.430

    Feature:  wk_rank_advantage
    Imp:      0.0360
    WIN mean: -0.945
    NON-WIN:  1.988
  KRPKR dataset saved: krpkr_positions.csv

=======================================================
PRINCIPLE EXTRACTION: KRKP
=======================================================
  Positions: 5,000
  Features:  16
  WIN: 4,744 (94.9%)
  DRAW: 256 (5.1%)

  Model validation:
              precision    recall  f1-score   support

     NON-WIN      1.000     1.000     1.000        51
         WIN      1.000     1.000     1.000       949

    accuracy                          1.000      1000
   macro avg      1.000     1.000     1.000      1000
weighted avg      1.000     1.000     1.000      1000


  Rank  Feature                            Imp      WIN  NON-WIN  Name
  ────────────────────────────────────────────────────────────────────────
  1     is_rook_pawn                    0.3135    0.219    0.754  T3[NEW?]  UNNAMED #1
  2     bp_to_promo                     0.2347    3.615    1.000  T3[NEW?]  UNNAMED #2
  3     bp_rank                         0.2216    3.615    1.000  T3[NEW?]  UNNAMED #3
  4     bk_to_promo                     0.0769    4.389    3.754  T3[NEW?]  UNNAMED #4
  5     bk_to_pawn                      0.0542    3.517    3.492  T3[NEW?]  UNNAMED #5
  6     wk_can_catch                    0.0248    0.610    0.098  T3[NEW?]  UNNAMED #6
  7     bk_ahead_of_pawn                0.0222    0.431    0.266  T3[NEW?]  UNNAMED #7
  8     wr_to_pawn                      0.0215    3.493    4.492  T3[NEW?]  UNNAMED #8
  9     wk_to_pawn                      0.0154    3.494    4.434  T3[NEW?]  UNNAMED #9
  10    wk_to_promo                     0.0114    4.327    4.938  T3[NEW?]  UNNAMED #10
  11    wr_mobility                     0.0019   12.991   13.242  T3[NEW?]  UNNAMED #11
  12    wr_behind_pawn                  0.0007    0.055    0.105  T3[NEW?]  UNNAMED #12
  13    wr_cuts_bk                      0.0007    0.241    0.230  T3[NEW?]  UNNAMED #13
  14    wr_on_pawn_file                 0.0005    0.112    0.129  T3[NEW?]  UNNAMED #14
  15    bk_supports_pawn                0.0001    0.059    0.066  T3[NEW?]  UNNAMED #15
  16    white_to_move                   0.0000    1.000    1.000  T3[NEW?]  UNNAMED #16

  TIER 3 CANDIDATES:

    Feature:  is_rook_pawn
    Imp:      0.3135
    WIN mean: 0.219
    NON-WIN:  0.754

    Feature:  bp_to_promo
    Imp:      0.2347
    WIN mean: 3.615
    NON-WIN:  1.000

    Feature:  bp_rank
    Imp:      0.2216
    WIN mean: 3.615
    NON-WIN:  1.000

    Feature:  bk_to_promo
    Imp:      0.0769
    WIN mean: 4.389
    NON-WIN:  3.754

    Feature:  bk_to_pawn
    Imp:      0.0542
    WIN mean: 3.517
    NON-WIN:  3.492

    Feature:  wk_can_catch
    Imp:      0.0248
    WIN mean: 0.610
    NON-WIN:  0.098

    Feature:  bk_ahead_of_pawn
    Imp:      0.0222
    WIN mean: 0.431
    NON-WIN:  0.266

    Feature:  wr_to_pawn
    Imp:      0.0215
    WIN mean: 3.493
    NON-WIN:  4.492
  KRKP dataset saved: krkp_positions.csv

============================================================
STEP 3: BUILDING COMBINED PRINCIPLE LIBRARY
============================================================

  Principles in library: 21
  KRK: 0 principles
  KRPKR: 11 principles
  KRKP: 10 principles

  Library saved: phase2_principles.json

  Generating visualisation: phase2_principle_space.png
  Saved: phase2_principle_space.png

=======================================================
PHASE 2 NAVIGATION VALIDATION
=======================================================

  Testing 302 positions...

  Correct:   298/300 (99.3%)
  Incorrect: 2

  ✓ STRONG VALIDATION (≥90%)
    Phase 2 methodology confirmed.
    Lucena/Philidor principles operational.

  Sample decisions:
  Config   Desc                            Oracle   Move   After   OK
  ─────────────────────────────────────────────────────────────────
  KRPKR    Random                             WIN   c3b4     WIN    ✓
  KRPKR    Random                             WIN   g4g5     WIN    ✓
  KRPKR    Random                             WIN   g3h4     WIN    ✓
  KRPKR    Random                             WIN   b6b7     WIN    ✓
  KRPKR    Random                            DRAW   b3b4     WIN    ✓
  KRPKR    Random                            DRAW   c2a2    DRAW    ✓
  KRPKR    Random                            DRAW   d2d3    DRAW    ✓
  KRPKR    Random                             WIN   e6e7     WIN    ✓
  KRPKR    Random                            DRAW   c1c2    DRAW    ✓
  KRPKR    Random                            DRAW   f7f8    DRAW    ✓
  KRPKR    Random                             WIN   h7g7     WIN    ✓
  KRPKR    Random                             WIN   b6b7     WIN    ✓

  Results saved: phase2_results.json

============================================================
PHASE 2 COMPLETE — SUMMARY
============================================================

  Positions analysed:  28,000
  Principles derived:  21
  Navigation accuracy: 99.3%
  Known confirmed:     3
  Tier 3 candidates:   10

  Output files:
    krk_positions.csv
    krpkr_positions.csv
    krkp_positions.csv
    phase2_principles.json
    phase2_principle_space.png
    phase2_results.json

  PHASE 3 READINESS:
  ✓ PHASE 3 READY
    Lucena and Philidor emerged from data.
    Navigation accuracy sufficient.
    Proceed to Phase 3: Middlegame Principles.

  Run interactive demo? (y/n): y

=======================================================
PHASE 2 INTERACTIVE DEMO
Supports: KRK, KRPKR, KRKP positions
Enter FEN or 'quit'
=======================================================

  Running demo positions:


  ──────────────────────────────────────────────────
  Lucena position
  FEN: 1K6/1P6/8/8/8/8/r7/2k1R3 w - - 0 1
  Oracle: WIN

  Config: KRPKR
  FEN:    1K6/1P6/8/8/8/8/r7/2k1R3 w - - 0 1...
  Legal moves: 14

  Top 3 moves:
    1. b8c7     Score: +0.7150 [Oracle: WIN]
       Principles: Pawn Advancement, UNNAMED #2, Lucena Position
    2. b8c8     Score: +0.6429 [Oracle: WIN]
       Principles: Pawn Advancement, UNNAMED #2, Lucena Position
    3. e1e8     Score: +0.6429 [Oracle: WIN]
       Principles: Pawn Advancement, UNNAMED #2, Lucena Position

  → RECOMMENDATION: b8c7
    Score: +0.7150
    Active: Pawn Advancement, UNNAMED #2, Lucena Position

  ──────────────────────────────────────────────────
  KRK — basic
  FEN: 8/8/8/8/8/2K5/8/1R2k3 w - - 0 1
  Oracle: WIN

  Config: KRK
  FEN:    8/8/8/8/8/2K5/8/1R2k3 w - - 0 1...
  Legal moves: 18

  Top 3 moves:
    1. c3d4     Score: +0.0000 [Oracle: WIN]
    2. c3c4     Score: +0.0000 [Oracle: WIN]
    3. c3b4     Score: +0.0000 [Oracle: WIN]

  → RECOMMENDATION: c3d4
    Score: +0.0000
    Active: none

  ──────────────────────────────────────────────────
  Pawn on 7th, king sheltering
  FEN: 8/1PK5/8/8/8/8/8/1k1r4 w - - 0 1
  Oracle: DRAW

  Config: UNKNOWN
  FEN:    8/1PK5/8/8/8/8/8/1k1r4 w - - 0 1...
  Legal moves: 8

  Top 3 moves:
    1. c7c8     Score: +0.0000 [Oracle: DRAW]
    2. c7b8     Score: +0.0000 [Oracle: DRAW]
    3. c7c6     Score: +0.0000 [Oracle: DRAW]

  → RECOMMENDATION: c7c8
    Score: +0.0000
    Active: none

  ──────────────────────────────────────────────────
  KRKP — stopping enemy pawn
  FEN: 8/8/8/8/8/8/1p6/1K2R3 w - - 0 1
  Error: unsupported operand type(s) for &: 'NoneType' and 'int'

──────────────────────────────────────────────────
  Enter FEN (or 'quit'): 1K6/1P6/8/8/8/8/r7/2k1R3 w - - 0 1
  Config:  KRPKR
  Oracle:  WIN

  Config: KRPKR
  FEN:    1K6/1P6/8/8/8/8/r7/2k1R3 w - - 0 1...
  Legal moves: 14

  Top 3 moves:
    1. b8c7     Score: +0.7150 [Oracle: WIN]
       Principles: Pawn Advancement, UNNAMED #2, Lucena Position
    2. b8c8     Score: +0.6429 [Oracle: WIN]
       Principles: Pawn Advancement, UNNAMED #2, Lucena Position
    3. e1e8     Score: +0.6429 [Oracle: WIN]
       Principles: Pawn Advancement, UNNAMED #2, Lucena Position

  → RECOMMENDATION: b8c7
    Score: +0.7150
    Active: Pawn Advancement, UNNAMED #2, Lucena Position

──────────────────────────────────────────────────
  Enter FEN (or 'quit'): 4k3/4P3/4K3/8/8/8/8/4R2r w - - 0 1
  Config:  KRPKR
  Oracle:  WIN

  Config: KRPKR
  FEN:    4k3/4P3/4K3/8/8/8/8/4R2r w - - 0 1...
  Legal moves: 16

  Top 3 moves:
    1. e6f6     Score: +0.6913 [Oracle: WIN]
       Principles: Pawn Advancement, UNNAMED #2, Lucena Position
    2. e6d6     Score: +0.6913 [Oracle: WIN]
       Principles: Pawn Advancement, UNNAMED #2, Lucena Position
    3. e6f5     Score: +0.6913 [Oracle: WIN]
       Principles: Pawn Advancement, UNNAMED #2, Lucena Position

  → RECOMMENDATION: e6f6
    Score: +0.6913
    Active: Pawn Advancement, UNNAMED #2, Lucena Position

──────────────────────────────────────────────────
  Enter FEN (or 'quit'): 8/1PK5/8/8/8/8/8/1k1r4 w - - 0 1
  Config:  UNKNOWN
  Oracle:  DRAW

  Config: UNKNOWN
  FEN:    8/1PK5/8/8/8/8/8/1k1r4 w - - 0 1...
  Legal moves: 8

  Top 3 moves:
    1. c7c8     Score: +0.0000 [Oracle: DRAW]
    2. c7b8     Score: +0.0000 [Oracle: DRAW]
    3. c7c6     Score: +0.0000 [Oracle: DRAW]

  → RECOMMENDATION: c7c8
    Score: +0.0000
    Active: none

──────────────────────────────────────────────────
  Enter FEN (or 'quit'): 8/8/8/8/8/8/1P6/1K2r1k1 w - - 0 1
  Config:  UNKNOWN
  Oracle:  DRAW

  Config: UNKNOWN
  FEN:    8/8/8/8/8/8/1P6/1K2r1k1 w - - 0 1...
  Legal moves: 2

  Top 3 moves:
    1. b1c2     Score: +0.0000 [Oracle: DRAW]
    2. b1a2     Score: +0.0000 [Oracle: DRAW]

  → RECOMMENDATION: b1c2
    Score: +0.0000
    Active: none

──────────────────────────────────────────────────
  Enter FEN (or 'quit'): K7/P7/8/8/8/8/8/1k2r2R w - - 0 1
  Config:  KRPKR
  Oracle:  WIN

  Config: KRPKR
  FEN:    K7/P7/8/8/8/8/8/1k2r2R w - - 0 1...
  Legal moves: 12

  Top 3 moves:
    1. a8b7     Score: +0.7150 [Oracle: WIN]
       Principles: Pawn Advancement, UNNAMED #2, Lucena Position
    2. a8b8     Score: +0.6429 [Oracle: WIN]
       Principles: Pawn Advancement, UNNAMED #2, Lucena Position
    3. h1h8     Score: +0.6429 [Oracle: WIN]
       Principles: Pawn Advancement, UNNAMED #2, Lucena Position

  → RECOMMENDATION: a8b7
    Score: +0.7150
    Active: Pawn Advancement, UNNAMED #2, Lucena Position

──────────────────────────────────────────────────
  Enter FEN (or 'quit'): 8/1P6/8/8/8/8/1K6/1k2r2R w - - 0 1
  Config:  KRPKR
  Oracle:  WIN

  Config: KRPKR
  FEN:    8/1P6/8/8/8/8/1K6/1k2r2R w - - 0 1...
  Legal moves: 3

  Top 3 moves:
    1. b2c3     Score: +0.7224 [Oracle: WIN]
       Principles: Pawn Advancement, UNNAMED #2, Lucena Position
    2. b2b3     Score: +0.7224 [Oracle: WIN]
       Principles: Pawn Advancement, UNNAMED #2, Lucena Position
    3. b2a3     Score: +0.7224 [Oracle: WIN]
       Principles: Pawn Advancement, UNNAMED #2, Lucena Position

  → RECOMMENDATION: b2c3
    Score: +0.7224
    Active: Pawn Advancement, UNNAMED #2, Lucena Position

──────────────────────────────────────────────────
  Enter FEN (or 'quit'): 8/1P1K4/8/8/8/8/8/rk5R w - - 0 1
  Config:  KRPKR
  Oracle:  WIN

  Config: KRPKR
  FEN:    8/1P1K4/8/8/8/8/8/rk5R w - - 0 1...
  Legal moves: 25

  Top 3 moves:
    1. d7e7     Score: +0.7189 [Oracle: WIN]
       Principles: Pawn Advancement, UNNAMED #2, Lucena Position
    2. d7e6     Score: +0.7189 [Oracle: WIN]
       Principles: Pawn Advancement, UNNAMED #2, Lucena Position
    3. d7d6     Score: +0.7189 [Oracle: WIN]
       Principles: Pawn Advancement, UNNAMED #2, Lucena Position

  → RECOMMENDATION: d7e7
    Score: +0.7189
    Active: Pawn Advancement, UNNAMED #2, Lucena Position

──────────────────────────────────────────────────
  Enter FEN (or 'quit'): 8/8/8/8/8/8/1P6/1K2r1k1 w - - 0 1
  Config:  UNKNOWN
  Oracle:  DRAW

  Config: UNKNOWN
  FEN:    8/8/8/8/8/8/1P6/1K2r1k1 w - - 0 1...
  Legal moves: 2

  Top 3 moves:
    1. b1c2     Score: +0.0000 [Oracle: DRAW]
    2. b1a2     Score: +0.0000 [Oracle: DRAW]

  → RECOMMENDATION: b1c2
    Score: +0.0000
    Active: none

──────────────────────────────────────────────────
  Enter FEN (or 'quit'): 8/8/8/8/8/8/p7/k3K2R w - - 0 1
  Config:  KRKP
  Oracle:  DRAW

  Config: KRKP
  FEN:    8/8/8/8/8/8/p7/k3K2R w - - 0 1...
  Legal moves: 14

  Top 3 moves:
    1. h1h8     Score: +0.0000 [Oracle: DRAW]
    2. h1h7     Score: +0.0000 [Oracle: DRAW]
    3. h1h6     Score: +0.0000 [Oracle: DRAW]

  → RECOMMENDATION: h1h8
    Score: +0.0000
    Active: none

──────────────────────────────────────────────────
  Enter FEN (or 'quit'): 8/8/8/8/8/8/3p4/3k1K1R w - - 0 1
  Config:  KRKP
  Oracle:  DRAW

  Config: KRKP
  FEN:    8/8/8/8/8/8/3p4/3k1K1R w - - 0 1...
  Legal moves: 11

  Top 3 moves:
    1. h1g1     Score: +0.2565 [Oracle: DRAW]
       Principles: UNNAMED #1, UNNAMED #6, UNNAMED #7
    2. h1h8     Score: +0.2416 [Oracle: DRAW]
       Principles: UNNAMED #1, UNNAMED #6, UNNAMED #7
    3. h1h7     Score: +0.2416 [Oracle: DRAW]
       Principles: UNNAMED #1, UNNAMED #6, UNNAMED #7

  → RECOMMENDATION: h1g1
    Score: +0.2565
    Active: UNNAMED #1, UNNAMED #6, UNNAMED #7

──────────────────────────────────────────────────
  Enter FEN (or 'quit'): 8/8/8/8/4K3/8/3p4/3k3R w - - 0 1
  Config:  KRKP
  Oracle:  DRAW

  Config: KRKP
  FEN:    8/8/8/8/4K3/8/3p4/3k3R w - - 0 1...
  Legal moves: 19

  Top 3 moves:
    1. e4f3     Score: +0.2416 [Oracle: DRAW]
       Principles: UNNAMED #1, UNNAMED #6, UNNAMED #7
    2. e4e3     Score: +0.2416 [Oracle: DRAW]
       Principles: UNNAMED #1, UNNAMED #6, UNNAMED #7
    3. e4d3     Score: +0.2416 [Oracle: WIN]
       Principles: UNNAMED #1, UNNAMED #6, UNNAMED #7

  → RECOMMENDATION: e4f3
    Score: +0.2416
    Active: UNNAMED #1, UNNAMED #6, UNNAMED #7

──────────────────────────────────────────────────
  Enter FEN (or 'quit'): 8/8/8/8/8/8/3p4/K2k3R w - - 0 1
  Config:  KRKP
  Oracle:  DRAW

  Config: KRKP
  FEN:    8/8/8/8/8/8/3p4/K2k3R w - - 0 1...
  Legal moves: 14

  Top 3 moves:
    1. a1b2     Score: +0.2416 [Oracle: DRAW]
       Principles: UNNAMED #1, UNNAMED #6, UNNAMED #7
    2. a1b1     Score: +0.2416 [Oracle: DRAW]
       Principles: UNNAMED #1, UNNAMED #6, UNNAMED #7
    3. a1a2     Score: +0.2192 [Oracle: DRAW]
       Principles: UNNAMED #1, UNNAMED #7, UNNAMED #9

  → RECOMMENDATION: a1b2
    Score: +0.2416
    Active: UNNAMED #1, UNNAMED #6, UNNAMED #7

──────────────────────────────────────────────────
  Enter FEN (or 'quit'): 8/8/8/8/8/2K5/8/1R2k3 w - - 0 1
  Config:  KRK
  Oracle:  WIN

  Config: KRK
  FEN:    8/8/8/8/8/2K5/8/1R2k3 w - - 0 1...
  Legal moves: 18

  Top 3 moves:
    1. c3d4     Score: +0.0000 [Oracle: WIN]
    2. c3c4     Score: +0.0000 [Oracle: WIN]
    3. c3b4     Score: +0.0000 [Oracle: WIN]

  → RECOMMENDATION: c3d4
    Score: +0.0000
    Active: none

──────────────────────────────────────────────────
  Enter FEN (or 'quit'): k7/8/1R6/8/8/8/8/7K w - - 0 1
  Config:  KRK
  Oracle:  WIN

  Config: KRK
  FEN:    k7/8/1R6/8/8/8/8/7K w - - 0 1...
  Legal moves: 17

  Top 3 moves:
    1. b6b8     Score: +0.0000 [Oracle: WIN]
    2. b6b7     Score: +0.0000 [Oracle: WIN]
    3. b6h6     Score: +0.0000 [Oracle: WIN]

  → RECOMMENDATION: b6b8
    Score: +0.0000
    Active: none

──────────────────────────────────────────────────
  Enter FEN (or 'quit'): 8/1P6/1K6/8/8/8/8/rk5R w - - 0 1
  Config:  KRPKR
  Oracle:  WIN

  Config: KRPKR
  FEN:    8/1P6/1K6/8/8/8/8/rk5R w - - 0 1...
  Legal moves: 21

  Top 3 moves:
    1. b6c6     Score: +0.7189 [Oracle: WIN]
       Principles: Pawn Advancement, UNNAMED #2, Lucena Position
    2. b6c5     Score: +0.7189 [Oracle: WIN]
       Principles: Pawn Advancement, UNNAMED #2, Lucena Position
    3. b6b5     Score: +0.7189 [Oracle: WIN]
       Principles: Pawn Advancement, UNNAMED #2, Lucena Position

  → RECOMMENDATION: b6c6
    Score: +0.7189
    Active: Pawn Advancement, UNNAMED #2, Lucena Position

──────────────────────────────────────────────────
  Enter FEN (or 'quit'): 8/3P4/3K4/8/8/8/8/r2k3R w - - 0 1
  Config:  KRPKR
  Oracle:  WIN

  Config: KRPKR
  FEN:    8/3P4/3K4/8/8/8/8/r2k3R w - - 0 1...
  Legal moves: 22

  Top 3 moves:
    1. d6e6     Score: +0.7189 [Oracle: WIN]
       Principles: Pawn Advancement, UNNAMED #2, Lucena Position
    2. d6c6     Score: +0.7189 [Oracle: WIN]
       Principles: Pawn Advancement, UNNAMED #2, Lucena Position
    3. d6e5     Score: +0.7189 [Oracle: WIN]
       Principles: Pawn Advancement, UNNAMED #2, Lucena Position

  → RECOMMENDATION: d6e6
    Score: +0.7189
    Active: Pawn Advancement, UNNAMED #2, Lucena Position

──────────────────────────────────────────────────
