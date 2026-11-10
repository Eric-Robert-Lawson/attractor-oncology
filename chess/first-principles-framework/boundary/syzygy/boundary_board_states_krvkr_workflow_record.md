# KRvKR workflow record: applying the programme to a segment

Run of 9 and 10 November 2026. Counts only. Nothing here is a verdict on a rule. This record carries the procedure of the programme as applied to KRvKR, in the order the programme specifies, with what was run, what stood in for a specified element, and what was not run. The exhaustive landscape is the calibration ground. No code is carried: the numbers in sections 3 and 4 are what any reproduction must match.

## 1. The procedure, applied

**A. Hard core: segment, fold, anchor, status.** Segment KRvKR (king and rook against king and rook, no castling, no pawns). Fold by the 16-image symmetry (8 board symmetries times colour exchange, KRvKR being its own colour exchange): 1,347,906 classes, 21,561,456 positions. Index convention: index = ((((wK*64+wR)*64+bK)*64+bR)*2+stm), square = 8*rank+file, stm 0 = White, class = minimum index over the 16 images. Status, as revised 6 November 2026 (`boundary_board_states.md` §14 item 8): drawn, reached under optimality from an anchor inside the segment, the draw realized by a cycle or a terminal draw reached under optimality; the losing deviation is a recorded property and not a condition. Anchor rule: start-selection §2c. A drawn class reached only through a won side's non-optimal move (a thrown win), or with no predecessor inside the segment, is outside the status set. Entry by an irreversible crossing from a larger material is outside this run. The move-again entry filter is withdrawn (§2c, 9 November 2026) and no step uses it.

**B. Compose the primitives for the material** (`boundary_board_state_primitive_catalog.md`). Targets: the defender's king (royal) and the defender's rook (non-royal). Each is a capture net, the royal one being the mating net (P-001, §4.3 to §4.7). The nets are composed by the double-threat rule of §4.10: nets on different targets are jointly compatible when no reply resolves both, and a check must be answered by the legal-move rule. Nets are read within a horizon h in plies (§4.2). The outcome of the capture net on the rook is read from the value of the position the capture reaches (§4.3, outcome clause): the KRvK leaf is a win by value unless the weak side is stalemated or can recapture the rook (checked against the KRvK table, 223,944 positions, 0 disagreements).

**C. The landscape as calibration.** The exact recursion over the segment, with no table in it, and the Syzygy tables (cursed and blessed counted as win and loss, 50-move rule not followed), compared class by class. Win(P,d), d odd: the mover captures the rook into a won leaf, or mates, or makes a move after which the opponent has at least one reply and every reply is Win(.,d-2). Lost(P,d), d even: the mover is mated (d=0), or has a legal move and every move is a non-capture leading to Win(.,d-1). A capture into a drawn leaf is neither. Run to a fixed point (two consecutive iterations with no new assignment). The landscape gives the value, and the exact depth, of every class.

**D. Compare the composition with the landscape, and locate the residual.** The composition is restricted to the primitive moves and the landscape shows where it falls short. The ladder levels are S0 capture or mate only, S1 plus checks, S2 plus rook attacks, S3 plus checks and rook attacks, S4 every move. Checks and rook attacks stand in for the net being engaged (the term C of §4.3); this is this record's reading.

**E. Read the residual, articulate it, run again.** The classes S3 does not win are the special cases. They were read in the landscape (section 3.3), articulated as a threat-creating move composed with a net on the other target, and the composition was run again as mode 3 (section 3.1).

**F. Certify and record.** Section 4. What is left is section 5.

## 2. Status of each element of the catalog in this run

| Element | Status in this run |
|---|---|
| P-001 nets on the king and the rook, terms C, F, X, B, O (§4.3) | C held by the check and attack flags; O by the capture leaf; F, X, B were not computed as terms |
| Outputs: residual reply sets R_F, R_X, R_B, R_O and status COMPLETED, ENGAGED, DORMANT (§4.4, §4.7) | not computed; completion within 3 plies by captures, mates, checks and attacks stands in for them |
| Double threat, check forcing the answer (§4.10) | run as the rule of mode 3, read as a double threat by this record |
| Mechanism tags (§4.5) | not run |
| Cost model (§4.8): cooperative distance CMD and the gap to the adversarial depth | CMD run at class level (section 3.3); the cost model itself is not calibrated |
| P-003 route distance, route-length set, tempo TP (§7) | not run; a residual reply set at horizon 1 (section 3.5) stands in for the reach |
| P-008 coverage race, TM (§8) | not run; the residual reply set of section 3.5 is a one-step form of the held replies |
| Hard core status and anchor (items 7 to 10) | run on the segment (section 3.4) |
| Closure claim (item 11) | not tested; one segment |

## 3. Results

### 3.1 The composition rule and its certification (mode 3)

Rule set at counted depth 7 (labels of this run):
- R0 capture net: the mover captures the opponent's rook; won by value unless the weak side is stalemated or can recapture the rook.
- R1 mating net: a move that mates.
- R2 forcing composition: a move that gives check or attacks the opponent's rook, after which every defender reply is lost, losses being those of R0 to R3.
- R3 quiet composition: a quiet move (neither a check nor a rook attack) that creates a threat, meaning that with the same board and the mover to move again, a capture net on the rook, a mating net on the king, or a forcing line of checks and rook attacks completes within 3 plies; after it every defender reply is lost by R0 to R3. This is the horizon h of §4.2 and §4.7 applied to the nets on the two targets, and the double-threat reading of §4.10.
- Lost: the mover is mated, or every move is a non-capture that leads to a won position by R0 to R3. Drawn: none of these within the window. The fold is the 16-image symmetry.
- Stored data: none.

Certification: mode 3 reproduces the table on all 1,347,906 classes, value and depth (won 392,459, lost 9,058, drawn 946,389; depth equal to |dtz| for wins and |dtz|+1 for lost classes that are not mated): 0 disagreements. Three other modes were run on the same classes with 0 disagreements: mode 0 every move, mode 1 forcing primitives plus one quiet move anywhere in the line, mode 2 forcing primitives plus 82 stored class indices (the 74 seeds and 8 lost classes of section 3.3, in the appendix).

Cost of one query on a fresh evaluator (new memo per query, seed 20261110, 600 queries per outcome, 474 for the special cases). Nodes expanded per query, mean; for mode 3 the threat tests are counted as nodes too (an earlier version of this record left them out and kept the memo between queries, which understated mode 3):

| Mode and move order | won | lost | drawn | special (466 + 8) |
|---|---|---|---|---|
| 0 every move, generator order | 3 | (earlier run: 4,790) | 9,232 | 435 |
| 0 every move, forcing first | 2 | 909 | 4,029 | 176 |
| 2 stored lists, forcing first | 2 | 102 | 417 | 40 |
| 3 threat rule, generator order | 10 | 5,713 | 11,500 | 1,048 |
| 3 threat rule, forcing first | 1 | 318 | 2,107 | 216 |
| 3 + forcing first + reach pruning by cooperative distance | 8 | 6,046 | 37,914 | 4,724 |

0 wrong values in every cell.

Readings. (1) Move order is the large effect, and it applies to proving a draw as well: a drawn position needs every move of the mover refuted (a "for all"), but each refutation is a "there exists" reply of the opponent, and forcing-first puts the rook capture and checks first, so the refuting reply is found early. Mode 0 alone falls from 9,232 to 4,029 nodes by order. (An earlier statement that order saves nothing for draws applied only to the top "for all" level and is withdrawn.) (2) Mode 3 restricts the mover to primitives but pays a threat test per quiet move, so in generator order it costs more nodes than mode 0; with forcing first it is about half of mode 0 (2,107 against 4,029 on drawn) and 318 against 909 on lost. Mode 2, which has the 82 stored classes in place of the threat test, is cheaper still. (3) Reach pruning by the cooperative distance is sound (0 wrong) but costs 18 times more than it saves: that distance is flat on this segment and it is computed by a cooperative search at every node. Per-net reach needs a cheap local distance from the segmented dependencies (P-003); that is not built.

What the stored lists are. Modes 0 and 3 store nothing and mode 2 stores 82 class indices (the 74 seeds and 8 lost classes, appendix). These are the places where the composition of the primitives did not yet derive the move, so they mark the interface the catalog still has to close. They are not compared with any table: a table lists outcomes, and the composition computes them from the position with the catalog's elements, so the two are different measurements and this record does not set one against the other. The cost figures above measure how much calculation a composition takes, in order to compare compositions and orderings with each other.

### 3.2 The landscape (calibration ground)

- Won classes 392,459, drawn 946,389 (15,138,064 positions), lost 9,058. The recursion equals the table: 0 mismatches.
- New won classes at d = 1, 3, 5, 7: 376,262 / 14,503 / 1,636 / 58. New lost classes at d = 0, 2, 4, 6: 603 / 7,882 / 465 / 108.
- Won d equals |dtz| for every won class; lost d equals |dtz|+1 for every lost class that is not mated; the 603 mated classes have d=0. Every position without an imbalance within that depth is drawn.
- The exact depth of every class is the cost in plies against best resistance, available for the calibration of §4.8.

### 3.3 What the landscape shows about the special cases (calibration findings, not rules)

**Ladder (won / lost classes found).** S0 376,262 / 8,485. S1 391,656 / 9,045. S2 378,332 / 8,831. S3 391,937 / 9,050. S4 392,459 / 9,058, equal to the exact recursion class by class. All restricted levels have 0 false wins and 0 false losses. Smallest restriction that wins a won class: capture or mate only 376,262; needs checks 15,394; needs rook attacks 83; needs checks and attacks together 198; needs a quiet move 522 (466 at d=5, 56 at d=7). Lost classes: capture or mate only 8,485; needs checks 560; attacks 0; both 5; quiet 8. Checks alone carry nearly the whole force region (S1 against S3).

**Quiet layer.** With a quota q of quiet moves for the winning side along the line, q = 0 equals S3 and q = 1 reproduces the exact recursion class by class, with the same minimal depth for every class; q = 2 to 4 add nothing. Every class that needs a quiet move needs exactly one: won d=5 466, d=7 56; lost 8. The 56 d=7 classes have a forcing first move and the quiet move later. The 466 d=5 classes have quiet-only first winning moves: king moves 358, rook moves 108.

**The position after the quiet move.** For the 466 quiet-first classes the passing quiet move leads to a position, defender to move, that is lost for the defender by forcing primitives alone (466 of 466; 598 passing moves out of 5,534 quiet moves: 394 classes with one, 30 with two, 24 with three, 18 with four). In all 466 the defender has exactly one legal king move and 7 to 14 rook moves. Of 4,326 replies, 3,800 are won at once (940 by the king taking the rook, 352 by the rook taking the rook, 2,508 by a mate in one) and 526 by a forcing win at d=3 (the single king move in all 466, and 60 rook moves). The same board with the winner to move is won at once in 458 of the 466 and within three plies in the other 8. Every passing move creates a threat (0 passing moves without one); a threat does not suffice (714 quiet moves create one, 598 pass, 116 do not; the two sets are identical in 350 of the 466 classes).

**Built backwards from the primitive-lost positions.** The forcing-only-lost classes (defender to move, lost by S3 alone) are 9,050: 603 mated, 7,882 at d=2, 465 at d=4, 100 at d=6. A quiet move can end only in a position where the defender is not in check and its rook is not attacked. The 603 mated are in check; of the other 8,447, 8,373 are not quiet (7,882 at d=2, 391 at d=4, 100 at d=6) and are entered by forcing moves, so their predecessors are already won by forcing lines. Only 74 classes, all at d=4, are quiet (the seeds), and in all 74 the defender has exactly one legal king move. The quiet un-moves of the 74 give exactly the 466. The defender's un-moves from the 466 give 3,216 candidates; those with every move covered (won by forcing, or one of the 466) are exactly the 8 lost classes at d=6, each with a single legal move that leads to one of the 466. The forcing un-moves of the 8 give exactly the 56 won classes at d=7. Nothing lies beyond them. Not reachable by a quiet layer: every forcing-lost class that is in check or has its rook attacked.

**Tempo reserve.** Reserve is the number of the mover's moves that are not lost at once. Take the classes where the mover is not in check, its rook is not attacked, it has no winning capture, and the opponent to move on the same board would win at once: with reserve 1 there are 109 classes, 70 lost in 4 plies and 39 drawn; with reserve 2 to 15 none is lost. The other 4 seeds have only a 3-ply threat and reserve 7 to 10 (one each at 7, 8, 9 and 10), every move covered at d=3. Together 74, all the seeds. Reserve 1 therefore holds for 70 of the 74 seeds and is not sufficient (39 drawn). The term that separates the 70 from the 39 is not computed.

**Cooperative distance and the slack of the winning moves.** Net completed means A has a winning capture or a mate (this record's reading of COMPLETED, §4.4). CMD is the number of plies for A to complete a net with both sides cooperating, the completing move included, computed over the class graph. It is 1 for 376,262 classes, 3 for 970,797 and 5 for 214 of the 1,347,273 classes where A, the mover, can complete one; 30 drawn classes have none. Of the 946,389 drawn classes, 946,293 have CMD 3 and 66 have CMD 5. The minimum-move pool is therefore flat on this material: cooperation reaches a net within 3 plies almost everywhere, drawn or won. What separates the outcomes is the defender's resistance, the gap between the adversarial depth d and CMD: for won classes the gap is 0 at d=1 (376,262 classes) and d=3 (14,503), 2 at d=5 (1,636) and 4 at d=7 (58). Slack of a move is CMD of the position it reaches (defender to move) minus CMD of the position it leaves plus 1; slack 0 means the move is on a cooperative shortest route. Of 4,876,998 non-capturing moves from won classes, 268,869 (5.5%) have slack 0. Of the 20,798 minimal-depth winning moves, 20,361 have slack 0 and 437 have slack 2; none has another slack. A minimal winning move with slack 0 exists in 392,022 of the 392,459 won classes. The 437 others are the classes whose every minimal winning move is a detour of 2 plies: 381 at d=5 and 56 at d=7, all with a forcing first move. All 598 passing quiet moves of the 466 quiet-first classes have slack 0. The reading of the 2-ply detour as the tempo adjustment of the route-length parity of §7 is this record's and is not tested. The cooperative minimum is a rank for the order of enumeration and not an answer: every move stays in the search, and a winning move may lie at any rank.

**Priority order (the cooperative minimum used to order the enumeration).** One rank per net, because the pool of nets is separate for each: cooperative distance to a capture net, A the mover (plies:classes): 1:372,719, 3:974,340, 5:214, none 633; to a mating net: 1:4,215, 3:64,424, 5:262,211, 6:1,015,500, none 1,556. Slack is read per net. Measured on the 16,197 won classes with d>=3, which have 17.86 legal moves on average: the expected number of moves examined until the first minimal-depth winning move, ties in random order (all d>=3 / the 466 quiet-first / the 437 slack-2 classes):
- generator order, no priority: 9.38 / 6.97 / 9.63.
- forcing first (mate, check, attack, quiet): 1.65 / 8.21 / 1.34.
- slack, either net: 8.87 / 6.23 / 17.82.
- slack per net (the lower), then the sum: 10.52 / 1.68 / 17.59.
- forcing first, then quiet moves by per-net slack (lower, then sum): 1.53 / 3.89 / 1.34.
The ranks are read from the class graph, that is, computed backward from the terminal conditions over every position. That is a table of distances to terminal points, which the programme does not use: the programme's rank is the local reach from the position by the progression of moves towards a mating condition, through the dependencies and the coverage, and it serves two purposes, the order of enumeration and the question whether a win is reachable within the remaining plies at all. The table here is an oracle that shows the ordering power an ideal local computation has, and the local computation is to be checked against it. An order saves nothing in proving a draw, where every move must be examined; a saving there needs a sound bound, and CMD is a lower bound on the adversarial depth (§4.8) that is flat here.

**Depth-indexed certification (partial certificates).** The composition rule of section 3.1 has the exact depth of every one of the 1,347,906 classes, so it is exact at every depth k. Won within k, lost within k, and the composition needed at that depth:
- k = 1 (lost within 0 and 2): won 376,262, lost 603 and 8,485. Capture net and mating net only.
- k = 3: won 390,765, lost 8,485. Plus the forcing composition (checks and rook attacks); every class first won at d=3 has a forcing winning move.
- k = 5 (lost within 4): won 392,401, lost 8,950. Plus the quiet threat move, needed by 466 of the 1,636 classes first won at d=5. The 465 classes first lost at d=4 need only the forcing composition.
- k = 7 (lost within 6): won 392,459, lost 9,058. The quiet move appears later in the line for 56 of the 58 classes first won at d=7, and the 8 lost classes at d=6 need it.
A certificate at depth k is the set won (or lost) within k together with "no forced win within k" for the remainder. Here nothing lies beyond depth 7, so the remainder at k=7 is the draws. The cooperative bound (cooperative reach cannot exceed the adversarial depth, §4.8) certifies "no win within 3" without search for 847 classes only; on this material the reachability bound is weak.

### 3.5 Local reach at each compositional moment (run on 10 November)

This runs the reading that the minimum reach is computed during the composition, at each step, to choose which route to continue, and that a fixed ranking is only one approach. Catalog elements run: the capture net on the rook and the mating net on the king as completion at horizon 1 (§4.3, §4.4), and a residual reply set at horizon 1, which is a one-step form of R_F, R_X, R_B, R_O. Left out: route distance by the segmented dependencies (§7), the coverage race and TM (§8), mechanism tags (§4.5), horizons above 1. No table and no class graph is read by it.

At a node, for each candidate move of the mover, R1 is the number of opponent replies after which the mover has no capture or mate at once; R1 = 0 means the net completes by force within 3 plies. The mover's candidates are ordered by R1 (capture or mate first). In the "for all" positions the replies are ordered by whether the mover would then have a completion at once (C1). Every probe is counted as a generated position, the same unit as the search; the first table below counts each call of the move generator, mate tests only on checking moves.

Selection of the first move. Positions with the true depth at least 3 (the exact depth is used only to mark which moves are winning, afterwards):

| Order | first move is a winning move, random won positions (74 classes) | same, quiet-first special classes (466) | mean position of the first winning move (74 / 466) |
|---|---|---|---|
| generator order | 0.0% | 55.6% | 14.15 / 3.18 |
| forcing moves first (fixed list) | 44.6% | 10.5% | 1.84 / 4.94 |
| local R1 over every candidate | 95.9% | 76.4% | 1.04 / 1.26 |
| local R1 over forcing candidates, quiet moves after | 94.6% | 10.5% | 1.34 / 4.94 |

Calculation, in generated positions per query (0 wrong values everywhere):

| Composition | won | lost | drawn | special (466 + 8) |
|---|---|---|---|---|
| mode 0, forcing first | 1 | 934 | 6,439 | 276 |
| mode 0, local R1 over every candidate | 581 | 4,957 | 159,168 | 7,418 |
| mode 0, local R1 over forcing candidates | 330 | 2,712 | 65,733 | 3,865 |
| mode 3, forcing first | 3 | 583 | 3,597 | 370 |
| mode 3, local R1 over every candidate | 562 | 2,498 | 8,551 | 4,041 |
| mode 3, local R1 over forcing candidates | 328 | 1,812 | 3,293 | 1,759 |

Findings.
1. The residual set at horizon 1, computed from the position with the two nets and nothing else, identifies the move that continues a minimal winning route as the first choice in 96% of random won positions, and in 76% of the quiet-first special classes, where a fixed list of move kinds puts a winning move first in 10%. So the quiet moves the earlier record had to derive backwards are mostly visible to a local reach measure. This is the first result in which the choice of the quiet move comes from the position's own nets.
2. With that order the search itself is small: in mode 3 the search nodes per query fall from 428 to 76 on drawn positions and from 104 to 16 on lost positions. The calculation moves from the search into the reach measure, which is computed for every candidate at every node.
3. At horizon 1 the measure is accurate but not economical to compute for every candidate; restricted to forcing candidates it is cheaper, and then it is no more accurate than the fixed list on the special classes (10.5%), because those moves are quiet. The 76% needs the quiet moves in the measure. Cost and accuracy trade off at this horizon, and the way out the programme names is the route distance of §7 and the coverage race of §8, which should hold the reach between nodes instead of recounting the replies at each one. That is not built.
4. About 24% of the quiet-first special classes still do not have the winning move first. Those are candidates for the residual sets R_F, R_X, R_B, R_O at a longer horizon, and are the next material for the catalog (the interface of the quiet layer in section 3.3).
5. Cooperative reach used as a cut (earlier row, "reach pruning") is not this: it asks whether any reply allows completion, and the cooperative distance is flat on this segment, so it removes few routes. The residual set asks about all replies, and it discriminates.
6. Samples: the random won positions with depth at least 3 are 74 of 1,500 drawn, a small sample; the 466 are all of them.

### 3.4 Status set

- Drawn classes 946,389. Strongly connected components 12,724, one with more than one class (933,666 classes), no self-loops.
- Anchor kinds (classes / positions): 1 on the cycle object 933,666 / 14,934,624; 2 off cycle, reached from it 609 / 9,744; 3 thrown-win only 11,819 / 188,976; 4 no predecessor in the segment 295 / 4,720 (294 have the mover in check).
- Status set (kinds 1 and 2): 934,275 classes, 14,944,368 positions. Outside it, drawn: 12,114 classes (kinds 3 and 4).
- Recorded properties: White-to-move cycle members with a losing move 412,342, Black-to-move 397,749; status-set classes whose mover has a losing move 810,680.
- Form of the certificate for an equilibrium region. A depth-bounded search cannot certify it. It is certified by closure: a finite set S of positions reached under optimality in which every position has a value-preserving move that stays in S, in which every move leaving S reaches a position certified won for the opponent by a depth-bounded certificate (the thrown wins, section 3.3), and which holds a cycle and no mating terminal. Then "no forced win within k" holds without limit on k. In KRvKR: drawn positions with a winning move 0; edges from a lost class into a drawn class 0; one cycle object of 933,666 of the 946,389 drawn classes; the 12,723 others are 609 reached from it by drawing moves, 11,819 reached only by a thrown win and 295 with no predecessor. The closure here was built by enumerating the region. That the primitives supply it as an invariant, the defence form of §4.7 (invariant-preserving moves: after which, for every reply of the attacker within h plies, no net is COMPLETED), is the clause of item 8 that no run has tested.
- 2,792 drawn classes have every drawing move a capture (44,672 positions), none on the cycle object; they stay drawn equilibrium states.

## 4. Certification

| Result | Checked by | Outcome |
|---|---|---|
| Recursion, wins, draws, losses | Class by class against the table | 0 mismatches over 1,347,906 classes |
| Depth | Against \|dtz\| | Equal as stated in section 3.2 |
| Capture-leaf rule | python-chess against the KRvK table, every legal position, weak side to move (223,944) | 0 disagreements |
| Minimal depth and winning-move structure | python-chess, own recursion, no table, random samples (seed 20261109): d = 1, 3, 5, 7 and lost d = 2, 4, 6 | 0 disagreements |
| Drawn classes | python-chess, 150 sampled | none won within 7 or lost within 8 |
| Move generator | Counts of legal, losing and staying moves for every drawn class, and total edges (20,137,706), against the stored files | 0 mismatches |
| Components | Rebuilt without the stored file | match exactly |
| Anchors | Reverse move generation: all 295 kind-4 classes have no legal predecessor; 400 random drawn classes agree on thrown-win and drawn-predecessor existence | 0 disagreements |
| Ladder | python-chess, own restricted recursion, 1,044 membership tests over S0 to S3, samples from every category | 0 disagreements |
| Composition rule (mode 3) and the three other modes | Every class against the table, value and depth; cost on fresh-evaluator samples (seed 20261110) | 0 disagreements over 1,347,906 classes in each of the four modes; 0 wrong in every cost sample |
| Quiet layer, all 466 quiet-first classes | python-chess, own restricted recursion: a forcing-only line never wins; a quiet move into a forcing-only-lost position exists; defender king moves; pass test; threat versus passing move sets | 0 disagreements with section 3.3 |
| Quiet layer built backwards | python-chess unmove generation and own legality, S3 sets from the ladder table: 74 seeds to the 466, 3,216 candidates to the 8, the 8 to the 56 | exact equality with the C++ sets at every step |
| Tempo reserve | C++, every class meeting the conditions of section 3.3 | counts as stated; not independently re-derived |
| Cooperative distance, slack and the priority orders | C++, relaxation over the class graph, per net and combined | counts as stated; not independently re-derived |
| Local reach by the residual reply set at horizon 1 (section 3.5) | C++, fresh evaluator per query, probes counted as generated positions, first winning move marked afterwards with the exact recursion | 0 wrong values in every cost sample; selection counts as stated; not independently re-derived |

A first run of the python-chess checker disagreed on 94 samples at depth 5 and beyond. The cause was in the checker (its Lost ignored a mated position once d>0, and its Win skipped mate in one at d>=3). It was fixed and rerun. The C++ program was not changed.

## 5. Residual progress (where this run stopped)

In the order the specification asks for them:
1. Compute the terms F, X, B and O for the nets on the king and on the rook, and the residual reply sets R_F, R_X, R_B, R_O of §4.4, and test whether the double-threat rule of §4.10 over them reproduces the rule of mode 3 and the ladder levels. The check and attack flags and the 3-ply completion are the stand-ins to replace.
2. Tag the residual replies with the mechanisms of §4.5 whose establishment would remove them, and test whether the tags select the passing quiet moves without enumerating the quiet moves. The threat condition alone leaves 116 of 714 quiet moves that do not pass.
3. Compute the per-net route distance (§7) from the movement relation, held between nodes, so that the reach is computed from the position and not recounted at each node (section 3.5 does it by a residual reply set at horizon 1, accurate, not economical), and test that it reproduces the ranks of section 3.3. Run the coverage race (§8) on the defender's replies: the replies covered by a completed net and those held open, and whether the attacker's move shrinks the held set, as a move generator in place of enumeration. The cooperative distance is flat on this material, so this is where the content lies. Compute the route distance and the parity of the route-length set (§7) for the attacker's unit, for the 109 reserve-1 classes and for the 437 slack-2 classes, and test whether it separates the 70 lost from the 39 drawn.
4. Calibrate the cost model (§4.8) against the exact depths of section 3.2.
5. Close the status and closure questions on more than one segment (item 11), including the pawn fold (left-right reflection only) and the castling fold (not run), and entry by crossing from a larger material.

Other limits of this record:
- The independent checks are samples except the table comparison, the move-generator comparison, the unanchored-class check and the mode comparisons, which are exhaustive. Both programs are the author's and share one reading of the programme; an independent recount is outside this session.
- "Forcing", "quiet", "threat", "seed" and the ladder levels are labels of this run, not terms of the other documents.
- The thrown-win limit (a drawn position reached only through a won side's non-optimal move is outside the status set unless it is an anchor entered by an irreversible crossing) is the designer's provisional statement of 9 November 2026, with the designer not certain of it. It affects 11,819 of 946,389 drawn classes, next to the force region, and none of the force-region results.
- The composition was run on one segment. Whether the stored interface stays as small on larger segments is not tested. Cost figures are samples of 600 per outcome and 474 for the special cases. The code of the rules is not counted.

## 6. Reproduction

Reproduce with the conventions of section 1.A, the procedure of section 1 in order, and the stopping rule of 1.C. A reproduction is accepted when it matches the counts of section 3, each against its check in section 4. Random samples use seed 20261109 (checkers) and 20261110 (cost samples).

## Appendix: stored class indices of mode 2 (convention of section 1.A)

Seeds (74, defender to move, quiet, lost in 4): 222532 222548 222564 222580 526488 526490 526492 526494 526504 526506 526508 526510 526520 526522 526524 526526 526536 526538 526540 526542 526552 526554 526556 526558 526568 526570 526572 526574 526584 526586 526588 526590 1059098 1059100 1059102 1059114 1059116 1059118 1059130 1059132 1059134 1059146 1059148 1059150 1059162 1059164 1059166 1059178 1059180 1059182 1059194 1059196 1059198 1591708 1591710 1591724 1591726 1591740 1591742 1591756 1591758 1591772 1591774 1591788 1591790 1591804 1591806 1608080 1608096 1608112 1608128 1608144 1608160 1608176.

Lost classes needing a quiet move (8, lost in 6): 524738 524754 524770 524786 525762 525778 525794 525810.
