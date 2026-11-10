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
8. Status set = kinds 1 and 2. Record the losing deviation per class as a property.

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
| Quiet layer, all 466 quiet-first classes | python-chess, own restricted recursion: a forcing-only line never wins; a quiet move into a forcing-only-lost position exists; defender king moves; pass test; threat versus passing move sets | 0 disagreements with the counts of section 3 |

A first run of the python-chess checker disagreed on 94 samples at depth 5 and beyond. The cause was in the checker (its Lost ignored a mated position once d>0, and its Win skipped mate in one at d>=3). It was fixed and rerun. The C++ program was not changed.

## 5. What is not covered

- The independent checks are samples except the table comparison, the move-generator comparison and the unanchored-class check, which are exhaustive. Both programs are mine and share my reading of the programme. An independent recount is outside this session.
- "Forcing" and "quiet" are labels of this run, not terms of the other documents.
- No tempo or cost term (§4.8, TP, TM) was used. The quiet layer is covered by a one-move composition of the primitives already derived (section 3), with the candidate quiet moves enumerated. A generator that proposes the passing move without trying the quiet moves is not derived: the threat condition narrows 5,534 quiet moves to 714 and keeps all 598 passing, and the coverage test decides the rest. The quiet layer was checked in python-chess for the 466 quiet-first classes and not for the 56 quiet-later or the 8 lost classes, which rest on the C++ run alone. The ply of the quiet move in the 56 and the 8 is not located.
- Pawn and castling folds were not run. Entry by crossing from a larger material was not run.
- The thrown-win limit (a drawn position reached only through a won side's sub-optimal move is outside the status set unless it is an anchor entered by an irreversible crossing) is the designer's provisional statement of 9 November 2026. The designer is not certain of it. It affects 11,819 of 946,389 drawn classes, next to the force region, and none of the force-region results.

## 6. Reproduction

Reproduce with the conventions of section 1, the steps of section 2 in order, and the stopping rule of step 3. A reproduction is accepted when it matches the counts of section 3, each against its check in section 4. Random samples use seed 20261109.
