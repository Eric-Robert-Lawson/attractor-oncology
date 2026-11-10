# KRvKR workflow record: force region rule set and status set

Run of 9 November 2026. Counts only. Nothing here is a verdict on a rule. This record replaces the scripts and tables of the run as the carried-forward certification: the numbers in section 3 are what any reproduction must match. No code is carried.

## 1. Scope and inputs

- Segment: KRvKR, white king, rook and black king, rook, no castling, no pawns. Fold by the 16-image symmetry (8 board symmetries times colour exchange, KRvKR being its own colour exchange). 1,347,906 classes, 21,561,456 positions.
- Index convention: index = ((((wK*64+wR)*64+bK)*64+bR)*2+stm), square = 8*rank+file, stm 0 = White; class = minimum index over the 16 images.
- Syzygy tables (counted: cursed and blessed as win and loss, 50-move rule not followed) are comparison only. No table value enters any recursion.
- Definitions used, as revised 6 November 2026: boundary status = drawn, reached under optimality from an anchor inside the segment, the draw realized by a cycle or a terminal draw reached under optimality. The losing deviation is a recorded property and not a condition. Anchor rule: start-selection §2c. A drawn class reached only through a won side's non-optimal move (a thrown win), or with no predecessor inside the segment, is outside the status set. The exception, entry by an irreversible crossing from a larger material, is outside this run's scope.
- The move-again entry filter is withdrawn (start-selection §2c, revision of 9 November 2026). It plays no role in the force region and excluded equilibrium positions. No run step uses it.

## 2. Procedure

1. Enumerate classes and generate moves with an own generator. A capture of the opponent's rook leads to a KRvK leaf: a win by value if the weak side is mated or has a legal move and cannot capture the rook; a draw if it is stalemated or can recapture.
2. Win(P,d), d odd: the mover has a capture into a won leaf, or a quiet move after which the opponent has at least one reply and every reply is Win(.,d-2); a mate in one counts at d=1. Lost(P,d), d even: the mover is mated (d=0), or has a legal move and every move is quiet and leads to Win(.,d-1). A capture into a drawn leaf is neither.
3. Run the class-level recursion to a fixed point, recording the minimum depth per class. Stop after two consecutive iterations with no new assignment.
4. Compare class by class with the table (win, draw, loss) and with |dtz|, after the recursion.
5. Composition ladder: rerun step 2 with the winning side's candidate moves restricted to S0 capture or mate only, S1 plus checks, S2 plus rook attacks, S3 plus checks and rook attacks, S4 every move. The defender keeps every legal reply. Restricted sets are subsets of the exact set by construction.
6. Quiet layers (9 November 2026): rerun the ladder recursion with a quota q of quiet moves for the winning side along the line, q = 0 to 4. A quiet move is one that neither gives check nor leaves the opponent's rook attacked. Captures, mates, checks and rook attacks are unrestricted; the defender has every legal reply and no quota. q = 0 is S3 and a large q is the exact recursion. Then, for the quiet-first classes, list every quiet move that passes (the position after it, with the defender to move, is lost for the defender by forcing primitives alone), the defender's replies with how each is won, and the same board with the winner to move (the pass test).
7. Drawn classes: build the drawing-move graph, take strongly connected components independently of any stored file, then a predecessor pass by reverse move generation. Anchor kinds: 1 on the cycle object, 2 off the cycle and reached from it by drawing moves, 3 reached only from a thrown win, 4 no predecessor inside the segment.
8. Table-free evaluator (10 November 2026): evaluate a position by search to 7 plies with the rules of section 3 and no table, in three modes: 0 every move; 1 forcing primitives plus one quiet move; 2 forcing primitives plus the stored lists of section 3. Run on every class and compare with the table (value and depth); then on a fresh evaluator per position on random samples (seed 20261110) to count the nodes expanded.
9. Status set = kinds 1 and 2. Record the losing deviation per class as a property.

## 3. Results of record

**Force region (exact recursion, 7 plies for wins, 6 for losses)**
- Won classes 392,459, drawn 946,389 (15,138,064 positions), lost 9,058. Recursion equals the table: 0 mismatches.
- New won classes at d = 1, 3, 5, 7: 376,262 / 14,503 / 1,636 / 58. New lost classes at d = 0, 2, 4, 6: 603 / 7,882 / 465 / 108.
- Won d equals |dtz| for every won class; lost d equals |dtz|+1 for every lost class that is not mated; the 603 mated classes have d=0.
- Every position without an imbalance within that depth is drawn.

**Composition ladder (won / lost classes found)**
- S0 376,262 / 8,485. S1 391,656 / 9,045. S2 378,332 / 8,831. S3 391,937 / 9,050. S4 392,459 / 9,058, equal to the exact recursion class by class. All restricted levels have 0 false wins and 0 false losses.
- Smallest restriction that wins a won class: capture or mate only 376,262; needs checks 15,394; needs rook attacks 83; needs checks and attacks together 198; needs a quiet move 522 (466 at d=5, 56 at d=7). Composition by checks and attacks covers 391,937 of 392,459 (99.867%).
- Lost classes: capture or mate only 8,485; needs checks 560; attacks 0; both 5; quiet 8.
- The 522 won and 8 lost classes that no forcing line wins: see Quiet layer below. The 466 d=5 classes have quiet-only first winning moves: king moves 358, rook moves 108.

**Quiet layer (quota of quiet moves along the line)**
- q = 0 equals S3 (0 differences); q = 1 reproduces the exact recursion class by class (392,459 won, 9,058 lost, 0 false, same minimal depth for every class); q = 2 to 4 add nothing. Every class that needs a quiet move needs exactly one: won d=5 466, d=7 56; lost 8.
- The 56 d=7 classes have a forcing first move and the quiet move later in the line; none has a quiet first move.
- For the 466 quiet-first classes the passing quiet move is a move into a position, defender to move, that is lost for the defender by forcing primitives alone (466 of 466; first move 598 passing moves in all, 394 classes with one, 30 with two, 24 with three, 18 with four, out of 5,534 quiet moves).
- In that position, in all 466, the defender has exactly one legal king move and 7 to 14 rook moves. Of 4,326 replies, 3,800 are won at once (940 by the king taking the rook, 352 by the rook taking the rook, 2,508 by a mate in one) and 526 by a forcing win at d=3 (the single king move in all 466, and 60 rook moves).
- Pass test: the same board with the winner to move is won at once in 458 of the 466 and within three plies in the other 8. Every passing move creates a threat (0 passing moves without one). A threat does not suffice: 714 quiet moves create one, 598 pass, 116 do not. The two move sets are identical in 350 of the 466 classes.
- Reading in primitive terms: the quiet move is one capture or mating net threat created by a move, composed with a coverage net over the defender's reply set, every reply of which is covered by a capture or mate at d=1 or a forcing win at d=3. The term "threat" is a label of this run.

**Rule set at counted depth 7, table-free (labels of this run)**
- R0 capture net: the mover captures the opponent's rook; won by value unless the weak side is stalemated or can recapture the rook (the leaf rule, checked against the KRvK table, 223,944 positions, 0 disagreements).
- R1 mating net: a move that mates.
- R2 forcing composition: a move that gives check or attacks the opponent's rook, after which every defender reply is lost, the losses being those of R0 to R4.
- R3 quiet composition: a quiet move into one of the 74 seeds (won in 5 plies), or a forcing move into one of the 8 stored lost classes (won in 7 plies). Lost: the mover is mated, or is one of the 8, or every move is a non-capture that leads to a won position by R0 to R3. Drawn: none of these within the window. The fold is the 16-image symmetry.
- Stored data, 82 class indices in the convention of section 1, listed in the appendix: seeds 74, lost classes 8. They are the backward construction below and can be re-derived by it.
- Certification: all three modes reproduce the table on every one of the 1,347,906 classes, value and depth (won 392,459, lost 9,058, drawn 946,389; depth equal to |dtz| for wins and |dtz|+1 for lost classes that are not mated): 0 disagreements. Mode 0 needs no stored data.
- Cost of one query on a fresh evaluator, nodes expanded (mean, max) and time per query:
  - mode 0, 120 per outcome: won 3 (135), 0.01 ms; lost 4,790 (27,012), 17.8 ms; drawn 9,042 (19,734), 34.1 ms.
  - mode 1, 1,500 per outcome: won 5 (977), 0.02 ms; lost 4,068 (34,779), 15.2 ms; drawn 9,019 (30,070), 34.9 ms.
  - mode 2, 1,500 per outcome: won 2 (654), 0.01 ms; lost 2,780 (24,474), 10.3 ms; drawn 4,627 (23,729), 17.3 ms.
  - 0 wrong values in every sample.
- Storage. Syzygy KRvKR: 12,944 bytes of win/draw/loss table and 3,408 bytes of dtz table, 16,352 in all. Stored data of mode 2: 82 class indices, 328 bytes, plus the rules as code. Mode 0 stores nothing. The table is already small, so there is no storage gain to measure on this segment; the comparison is storage against compute. A drawn position costs about 4,600 nodes in mode 2 (9,000 in mode 0), against one table probe.

**Quiet layer built backwards (no trial of candidate moves)**
- Seeds. The forcing-only-lost classes (defender to move, lost by S3 alone) are 9,050: 603 mated, 7,882 at d=2, 465 at d=4, 100 at d=6. A position can be the end of a quiet move only if the defender is not in check and the defender's rook is not attacked. The 603 mated classes are in check. Of the other 8,447, 8,373 are not quiet (7,882 at d=2, 391 at d=4, 100 at d=6) and are entered by forcing moves, so their predecessors are already won by forcing lines. Only 74 classes, all at d=4, are quiet. In all 74 the defender has exactly one legal king move.
- Step 1. The quiet un-moves of the 74 seeds give exactly the 466 quiet-first won classes, all at d=5: none missing, none extra, each class's exact depth equal to its seed depth plus 1.
- Step 2. The defender's un-moves from the 466 give 3,216 candidate positions. Those in which every defender move is covered (won by forcing, or one of the 466) are exactly the 8 lost classes at d=6. Each has a single legal move, and it leads to one of the 466.
- Step 3. The forcing un-moves of the 8 give exactly the 56 won classes at d=7. No class lies beyond them, which matches q = 1 sufficing.
- Reachable: the special classes are the quiet un-move predecessors of the 74 seeds, the defender predecessors of those, and the forcing predecessors of the 8. Not reachable by a quiet layer: any forcing-lost class that is in check or has its rook attacked, which is every forcing-lost class at d=0, 2 and 6 and 391 of 465 at d=4.

**Status set**
- Drawn classes 946,389. Strongly connected components 12,724, one with more than one class (933,666 classes), no self-loops.
- Anchor kinds (classes / positions): 1 on the cycle object 933,666 / 14,934,624; 2 off cycle, reached from it 609 / 9,744; 3 thrown-win only 11,819 / 188,976; 4 no predecessor in the segment 295 / 4,720 (294 have the mover in check).
- Status set (kinds 1 and 2): 934,275 classes, 14,944,368 positions. Outside it, drawn: 12,114 classes (kinds 3 and 4).
- Recorded properties on the cycle object: White-to-move members with a losing move 412,342, Black-to-move members 397,749. Status-set classes whose mover has a losing move: 810,680.

## 4. Certification

| Result | Checked by | Outcome |
|---|---|---|
| Recursion, wins, draws, losses | Class-by-class against the table | 0 mismatches over 1,347,906 classes |
| Depth | Against |dtz| | Equal as stated in section 3 |
| Capture-leaf rule | python-chess against the KRvK table, every legal position, weak side to move (223,944) | 0 disagreements |
| Minimal depth and winning-move structure | python-chess, own recursion, no table, random samples (seed 20261109) across d = 1, 3, 5, 7 and lost d = 2, 4, 6 | 0 disagreements |
| Drawn classes | python-chess, 150 sampled | none won within 7 or lost within 8 |
| Move generator | Counts of legal, losing and staying moves for every drawn class and total edges (20,137,706) against the stored files | 0 mismatches |
| Components | Rebuilt without the stored file | match the stored file exactly |
| Anchors | Reverse move generation: all 295 kind-4 classes have no legal predecessor; 400 random drawn classes agree on thrown-win and drawn-predecessor existence | 0 disagreements |
| Ladder | python-chess, own restricted recursion, 1,044 membership tests over S0 to S3 on samples from every category | 0 disagreements |
| Table-free evaluator, three modes | Every class against the table: value and depth. Cost on fresh-evaluator samples | 0 disagreements over 1,347,906 classes in each mode; 0 wrong in every cost sample |
| Quiet layer built backwards: 74 seeds to the 466, 3,216 candidates to the 8, the 8 to the 56 | python-chess unmove generation and own legality, S3 sets from the ladder table | exact equality with the C++ sets at every step |
| Quiet layer, all 466 quiet-first classes | python-chess, own restricted recursion: a forcing-only line never wins; a quiet move into a forcing-only-lost position exists; defender king moves; pass test; threat versus passing move sets | 0 disagreements with the counts of section 3 |

A first run of the python-chess checker disagreed on 94 samples at depth 5 and beyond. The cause was in the checker (its Lost ignored a mated position once d>0, and its Win skipped mate in one at d>=3). It was fixed and rerun. The C++ program was not changed.

## 5. What is not covered

- The independent checks are samples except the table comparison, the move-generator comparison and the unanchored-class check, which are exhaustive. Both programs are mine and share my reading of the programme. An independent recount is outside this session.
- "Forcing" and "quiet" are labels of this run, not terms of the other documents.
- No tempo or cost term (§4.8, TP, TM) was used. The quiet layer is covered by a one-move composition of the primitives already derived (section 3), with the candidate quiet moves enumerated. The quiet layer is also built backwards from the seeds (section 3), so the candidate quiet moves are not tried forward. What is not derived is a static description of the 74 seeds: all have one defender king move, which is a shared property and is not shown to define them. The seeds are read from the forcing-only-lost set, which comes from the C++ run. The backward construction was run in python-chess (unmove generation, own legality) with the S3 win and loss sets taken from the C++ ladder table. The ply of the quiet move inside the 56 follows from the construction (the second ply after the forcing first move leads to the 8, whose single reply is one of the 466).
- The storage comparison covers one small segment. Whether the rule set and its stored lists grow more slowly than a table for larger segments is not tested. The cost figures are samples; the table comparison is exhaustive. The code of the rules is not counted in the 328 bytes.
- Pawn and castling folds were not run. Entry by crossing from a larger material was not run.
- The thrown-win limit (a drawn position reached only through a won side's sub-optimal move is outside the status set unless it is an anchor entered by an irreversible crossing) is the designer's provisional statement of 9 November 2026. The designer is not certain of it. It affects 11,819 of 946,389 drawn classes, next to the force region, and none of the force-region results.

## 6. Reproduction

Reproduce with the conventions of section 1, the steps of section 2 in order, and the stopping rule of step 3. A reproduction is accepted when it matches the counts of section 3, each against its check in section 4. Random samples use seed 20261109 (checkers) and 20261110 (cost samples).

## Appendix: stored class indices (convention of section 1)

Seeds (74, defender to move, quiet, lost in 4): 222532 222548 222564 222580 526488 526490 526492 526494 526504 526506 526508 526510 526520 526522 526524 526526 526536 526538 526540 526542 526552 526554 526556 526558 526568 526570 526572 526574 526584 526586 526588 526590 1059098 1059100 1059102 1059114 1059116 1059118 1059130 1059132 1059134 1059146 1059148 1059150 1059162 1059164 1059166 1059178 1059180 1059182 1059194 1059196 1059198 1591708 1591710 1591724 1591726 1591740 1591742 1591756 1591758 1591772 1591774 1591788 1591790 1591804 1591806 1608080 1608096 1608112 1608128 1608144 1608160 1608176.

Lost classes needing a quiet move (8, lost in 6): 524738 524754 524770 524786 525762 525778 525794 525810.
