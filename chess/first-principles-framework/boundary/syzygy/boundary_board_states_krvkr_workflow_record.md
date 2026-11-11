# KRvKR workflow record: applying the programme to a segment

Run of 9 and 10 November 2026. Counts only. Nothing here is a verdict on a rule. This record carries the procedure of the programme as applied to KRvKR, in the order the programme specifies, with what was run, what stood in for a specified element, and what was not run. The exhaustive landscape is the calibration ground. No code is carried: the numbers in sections 3 and 4 are what any reproduction must match.

## 0. Scope of the application to KRvKR (the designer's statement, 10 November 2026)

1. Use Syzygy to enumerate the boundary board states of KRvKR. The status is the one of `boundary_board_states.md` §14 items 7, 8 and 10 as revised on 6 November 2026: the position is drawn, it is reached under optimality from an anchoring point, and the draw is realized by a cycle or by a terminal draw reached under optimality. The losing deviation of each side, the width and the realization class are recorded properties and not conditions, so nothing is excluded for them.
2. Apply the primitives of the catalog to show that these states are equilibrium states: the registers and equilibrium moves of P-008 §8.11, and the defence form of P-001 §4.7.
3. From the compositions found against the Syzygy ground truth, build a rule set that certifies the boundary board states by primitive calculation, without using Syzygy, and agrees with Syzygy on 100% of the classes.

| Step | State at the end of 10 November |
|---|---|
| 1 | Done. Status set 934,275 classes, 14,944,368 positions (section 3.4). |
| 2 | In progress. The P-001 terms agree with the generator on 100% of classes where a net is complete (section 3B). The defence-form closure over all classes is run in sections 3C and 3D: 0 disagreements in value and depth over 1,347,906 classes, with no economy. |
| 3 | Not done. The composition rule of section 3.1 matches the table on 100% of classes, but it enumerates replies and recurses with the move generator. It is the calibration, not the rule set of step 3. |

## 1. The procedure, applied

**A. Hard core: segment, fold, anchor, status.** Segment KRvKR (king and rook against king and rook, no castling, no pawns). Fold by the 16-image symmetry (8 board symmetries times colour exchange, KRvKR being its own colour exchange): 1,347,906 classes, 21,561,456 positions. Index convention: index = ((((wK*64+wR)*64+bK)*64+bR)*2+stm), square = 8*rank+file, stm 0 = White, class = minimum index over the 16 images. Status, as revised 6 November 2026 (`boundary_board_states.md` §14 item 8): drawn, reached under optimality from an anchor inside the segment, the draw realized by a cycle or a terminal draw reached under optimality; the losing deviation is a recorded property and not a condition. Anchor rule: start-selection §2c. A drawn class reached only through a won side's non-optimal move (a thrown win), or with no predecessor inside the segment, is outside the status set. Entry by an irreversible crossing from a larger material is outside this run. The move-again entry filter is withdrawn (§2c, 9 November 2026) and no step uses it.

**B. Compose the primitives for the material** (`boundary_board_state_primitive_catalog.md`). Targets: the defender's king (royal) and the defender's rook (non-royal). Each is a capture net, the royal one being the mating net (P-001, §4.3 to §4.7). The nets are composed by the double-threat rule of §4.10: nets on different targets are jointly compatible when no reply resolves both, and a check must be answered by the legal-move rule. Nets are read within a horizon h in plies (§4.2). The outcome of the capture net on the rook is read from the value of the position the capture reaches (§4.3, outcome clause): the KRvK leaf is a win by value unless the weak side is stalemated or can recapture the rook (checked against the KRvK table, 223,944 positions, 0 disagreements).

**C. The landscape as calibration.** The exact recursion over the segment, with no table in it, and the Syzygy tables (cursed and blessed counted as win and loss, 50-move rule not followed), compared class by class. Win(P,d), d odd: the mover captures the rook into a won leaf, or mates, or makes a move after which the opponent has at least one reply and every reply is Win(.,d-2). Lost(P,d), d even: the mover is mated (d=0), or has a legal move and every move is a non-capture leading to Win(.,d-1). A capture into a drawn leaf is neither. Run to a fixed point (two consecutive iterations with no new assignment). The landscape gives the value, and the exact depth, of every class.

**D. Compare the composition with the landscape, and locate the residual.** The composition is restricted to the primitive moves and the landscape shows where it falls short. The ladder levels are S0 capture or mate only, S1 plus checks, S2 plus rook attacks, S3 plus checks and rook attacks, S4 every move. Checks and rook attacks stand in for the net being engaged (the term C of §4.3); this is this record's reading.

**E. Read the residual, articulate it, run again.** The classes S3 does not win are the special cases. They were read in the landscape (section 3.3), articulated as a threat-creating move composed with a net on the other target, and the composition was run again as mode 3 (section 3.1).

**F. Certify and record.** Section 4. What is left is section 5.

## 1b. The programme as the designer articulated it, checked against the documents (10 November)

The statements are the designer's, from this session, in the order given. Each is checked against the text of the catalog (`boundary_board_state_primitive_catalog.md`), the main document and `boundary_board_states_start_selection.md`; the sections cited in the table were read for this check, and the rest of those documents was not re-read.

| # | The designer's statement | Where the documents say it | Result |
|---|---|---|---|
| 1 | The hard core is the boundary board state; the belt is the three primitives, joint compatible and compositionally dependent: any two need the third for a calculation that is not trivial | main §14 items 7 to 11; catalog §11.1 (the three primitives), §9 | matches |
| 2 | The mating net is an outcome of reachability from the position, not a precomputed array of terminal points | catalog §4.7 (the deliveries are the attacker's legal moves that attack the target, at the position, to the horizon h); §7.3 and §7.5 (a route to a configuration with a gain); §9 (the roles of a move: completing the net, a capture net across boundaries, a move that restricts the opponent's optimal replies, which is a segmented dependency, navigation) | matches |
| 3 | The coverage net informs the mating net | catalog §9 (P-008 to P-001: an opened held term lets a net complete; the coverage race is a dependency of the optimality calculation of the nets, since a net cannot be evaluated for both sides without the occupancy and coverage of the other pieces); §8.3 | matches |
| 4 | The minimum number of moves to reach a net, even with cooperation, informs the priority of enumeration; it is a priority order and never the answer | catalog §4.8 (the cooperative cost is the helpmate distance, CMD <= DTM); §6 ("a heuristic may order candidates and may propose mechanisms. It certifies nothing"); §9 (cheaper approaches can prune, not definitive; the ordering is not a rule of the catalog, the workflow decides) | matches |
| 5 | The minimum reachability is computed at each compositional moment, to choose which route is more practical to continue | catalog §7.3: the route distance d_M(s, t, P) is the cheapest path in the navigation graph, an edge costing one tempo plus its exposure, with the route-length set; §7.5 terms RT and TP; §8.9 items 2 to 4: a perturbation whose change misses the support closure needs no branching | matches. The documents define the local calculation as a route distance on the navigation graph and the support closure. It is not a cooperative search and not a count of the opponent's replies. Section 3B runs the route distance with the defender held in place |
| 6 | Navigation towards the net is dependency based, chronological from the anchoring position; a forced sequence is a cadence that segmented dependency race conditions calculate | catalog §7 (ENABLE and SUPPORT edges, the segmented process); §11.2 (check as a segmented dependency that restricts the viable moves); start-selection §2c (anchor and trajectory) | matches |
| 7 | The search runs inside the primitive calculations and does not enumerate every branch | catalog §8.9 items 3 and 4 (exact pruning, exhaustive coverage over what remains, as case lists decided by primitive calculation); §9 last lines; start-selection §2c ("within a segment ... a forcing sequence exists ... or the position is drawn, and the draw needs exhaustive coverage") | matches. Exhaustive coverage stays the certification and the insurance on the first entries |
| 8 | Depth-indexed partial certification, with the uncertain region as the focus of community effort | start-selection §2c: the exclusion form with its depth bound d; the coverage depth recorded as a scoped parameter and not as certified; the scoped region of uncertainty | matches |
| 9 | An equilibrium region cycles under optimality and never reaches a mating terminal | main §14 items 7, 8 and 10 (the cycle, the terminal draw); catalog §4.7 defence form, §8.11 (the equilibrium moves are those whose effect vector has no OPENED entry on the mover's own register) | matches. The certificate by invariant is specified there and was not run; the run's certificate is an enumeration |
| 10 | The special cases arise where placement removes the moves that keep the equilibrium, for king and for rook, measurable from placement | catalog §4.7 defence form (the invariant-preserving moves); §4.5 mechanisms F1, F2 and the mirror D-F; §8.3 held terms, residual replies, support, opening sets, key elements, joint support; §8.7 OP2 (joint spend) and OP4 (denial); the board edge as an input (§1, the P-004 row; §7.3 missing edges) | matches the framework. The quantities the documents ask for are held terms and their supports. The run computed stand-ins: the moves not removed by the generator and by a capture or mate at once (sections 3.7 to 3.9) |
| 11 | The primitives are piece-agnostic and board-agnostic; the edge and the movement are inputs; a special case recorded inside a primitive would make it a subspace of piece composition | catalog §1 (admission rule 4), §6, §10 (the P-004 row) | an arrangement such as a rook beside its king on an edge names a piece and an edge. It is an instance, and at most a candidate derived form (the mechanisms F1, the board edge as missing edges, and F2, as in DF-A), with no definitional weight. It is not a rule and is not to be written into a primitive |
| 12 | One workflow record, no scripts, each computation says which catalog element it runs and what it leaves out | the designer's instruction | followed; the code is not delivered |

What the run did that follows the programme. Both targets, the king as the royal target and the rook as the non-royal target, with the outcome clause read from the value of the position the capture reaches, tested against the KRvK table as the sandbox (start-selection §2c, the exit condition); the double-threat composition of §4.10 at horizon 3; the defence form of §4.7 at levels 1 and 2; the hard-core status on the segment (section 3.4); depth-indexed certification; the landscape used as calibration and not as a rule.

What it did not run, in the documents' terms. P-003 on KRvKR: the navigation graphs of the king and the rook, the route distance with exposure, the route-length set, the term TP. P-008 on KRvKR: the register of each side (held terms with residual replies), the supports, opening sets, key elements and joint support, the effect vector of a move, the terms OP, CR, FS, TM. P-001 as written: the terms C, F, X, B and O as terms, the residual partition R_F, R_X, R_B, R_O, the mechanism tags. The predictions P3-1, P3-2, P8-1, P8-2 and P8-3, none of which was stated before a run. The invariant-based equilibrium certificate. A second segment.

What the programme asks to run next, in its own terms. (1) P-003: build the navigation graphs on the board for the king and the rook with their dependency sets, compute route distance with exposure from the position, and compare the reach it gives with the ranks of section 3.3; this is the local calculation of statement 5. (2) P-008: the defender's register at the node and at the seeds, with supports, opening sets and key elements, and the effect vector of the attacker's moves, naming the mechanism (OP1 to OP6) of each; read the quiet move as a perturbation by that vector. (3) P-001 as §4.7 writes it. Predictions will be stated in the record before each is run. The designer's open items are kept apart: the thrown-win anchor limit (provisional), and whether the operating-rule block goes into the main document.

## 2. Status of each element of the catalog in this run

| Element | Status in this run |
|---|---|
| P-001 nets on the king and the rook, terms C, F, X, B, O (§4.3) | computed from coverage sets and occupancy, royal target and non-royal target (section 3B): 0 disagreements with the generator, with F read with the outcome at the destination square for the rook |
| Outputs: residual reply sets R_F, R_X, R_B, R_O and status COMPLETED, ENGAGED, DORMANT (§4.4, §4.7) | not computed at every position. The defence form of §4.7 as an invariant over the status set has not been run (section 0, step 2) |
| Double threat, check forcing the answer (§4.10) | run as the rule of mode 3, read as a double threat by this record |
| Mechanism tags (§4.5) | not run |
| Cost model (§4.8): cooperative distance CMD and the gap to the adversarial depth | CMD run at class level (section 3.3); the cost model itself is not calibrated |
| P-003 route distance, route-length set, tempo TP (§7) | route distance of the attacker's king and rook with the defender held in place, run (section 3B); TP not run |
| P-008 coverage race, TM (§8) | register at the lost classes, run (section 3B); effect vectors, equilibrium moves and TM not run |
| Hard core status and anchor (items 7 to 10) | run on the segment (section 3.4) |
| Closure claim (item 11) | not tested; one segment |
| Predictions P3-1, P3-2, P8-1, P8-2, P8-3 (§7.13, §8.13) | not tested. PA to PF of section 3A were stated before the run of section 3B |

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

### 3.4 Status set

- Drawn classes 946,389. Strongly connected components 12,724, one with more than one class (933,666 classes), no self-loops.
- Anchor kinds (classes / positions): 1 on the cycle object 933,666 / 14,934,624; 2 off cycle, reached from it 609 / 9,744; 3 thrown-win only 11,819 / 188,976; 4 no predecessor in the segment 295 / 4,720 (294 have the mover in check).
- Status set (kinds 1 and 2): 934,275 classes, 14,944,368 positions. Outside it, drawn: 12,114 classes (kinds 3 and 4).
- Recorded properties: White-to-move cycle members with a losing move 412,342, Black-to-move 397,749; status-set classes whose mover has a losing move 810,680.
- Form of the certificate for an equilibrium region. A depth-bounded search cannot certify it. It is certified by closure: a finite set S of positions reached under optimality in which every position has a value-preserving move that stays in S, in which every move leaving S reaches a position certified won for the opponent by a depth-bounded certificate (the thrown wins, section 3.3), and which holds a cycle and no mating terminal. Then "no forced win within k" holds without limit on k. In KRvKR: drawn positions with a winning move 0; edges from a lost class into a drawn class 0; one cycle object of 933,666 of the 946,389 drawn classes; the 12,723 others are 609 reached from it by drawing moves, 11,819 reached only by a thrown win and 295 with no predecessor. The closure here was built by enumerating the region. That the primitives supply it as an invariant, the defence form of §4.7 (invariant-preserving moves: after which, for every reply of the attacker within h plies, no net is COMPLETED), is the clause of item 8 that no run has tested.
- 2,792 drawn classes have every drawing move a capture (44,672 positions), none on the cycle object; they stay drawn equilibrium states.

## 3A. Primitive run on KRvKR: predictions registered before the run (10 November, 13:40)

Catalog elements run: P-001 §4.3 terms C, F, X, B, O for the royal target (the defender's king) and the non-royal target (the defender's rook), with the residual partition of §4.4; P-003 §7.3 navigation graphs of the king and the rook with dependency sets and edge exposure, route distance d_M, and the reach of the attacker's units to a configuration at which a net completes, with the defender's pieces held in place; the double-threat composition of §4.10. P-008 §8.3 is read at the lost classes (held terms, residual replies, supports, key elements; results in section 3B). Left out: the tempo budget TP and TM as a comparison against an opposing race (KRvKR has no pawn race), the cost model §4.8 beyond the route distance, mechanism tags beyond F1, F2, X3.

Predictions, stated before any of the code below was run:
1. **PA.** For the royal target, C, F, X and B as §4.3 defines them hold exactly when the defender is checkmated: 0 disagreements with the move generator over every class position with the defender to move.
2. **PB.** For the non-royal target, the net holds by the five terms exactly when no reply of the defender resolves it by R_F, R_X, R_B or R_O, and the only disagreements with the move generator are replies that give check (a counter-net, §4.10): none other.
3. **PC.** The reach of the attacker's units to a net completion, with the defender held in place, in attacker moves k, gives 2k - 1 <= the exact depth for every won class (it is a cooperative lower bound, §4.8).
4. **PE.** The gap between the exact depth and 2k - 1 is larger for the 466 and the 56 than for the forcing-solved classes of the same depth.
5. **PF.** The support closure of §8.9, computed for the defender's register, is the whole board in KRvKR, so §8.9 item 3 gives no pruning here.

## 3B. Primitive run on KRvKR: results (10 November, after section 3A)

Every figure below is a count over the 1,347,906 classes (or over the ladder rows named). Nothing here is a verdict on a rule.

**P-001 terms, royal target (PA).** C, F, X and B as §4.3 defines them, computed from the coverage sets and occupancy (Cov_A with the target king lifted, Occ_D, legal attackers, between-squares), agree with the move generator on every class: 0 disagreements. The net holds on 603 classes, all of them the 603 mated classes.

**P-001 terms, non-royal target (PB).** The target is the defender's rook, attacked by the attacker's rook or king in 367,815 classes where the defender has a legal reply (a further 42 attacked-rook classes have no legal reply and are outside the comparison). Terms C, F, X, B and O were computed, with O read as §4.3 states it (the capture is favourable and its outcome is read by the value of the position the capture reaches; each defender reply is checked for a favourable capture left to the attacker).
- F read literally (a covered square closes the flight) disagrees with the generator on 888 classes.
- F read with the outcome at the destination (a covered square where the capture is unfavourable counts as a flight) agrees on all 1,347,906 classes (net holds on 5,818), 0 disagreements. Recorded as a residual interface of §4.3 for a non-royal target: F needs the outcome at the square, which is the clause O already carries. The catalog text is unchanged.
- An earlier version of the code that ignored the legality of the attacker's king capture after a defender king move disagreed on 250 classes. That was an error in the code, found and fixed here, not a finding about the catalog.
- No disagreement was a counter-net check, so the exception named in PB was not needed.

**P-003 reach with the defender held in place (PC, PE).** The movement relation of the attacker's king and rook (§7.3), the defender's pieces fixed, completing move = a mate or a capture the leaf values as a won position, k = least number of attacker moves. Two variants: movement relation only (R), and R with edges that leave the moved rook attacked and unprotected denied (exposure, E). Compared with the exact depth, bound 2k - 1 for a won class (attacker to move) and 2k for a lost class (defender to move).
- PC is not confirmed. The bound is violated by 946 of the 392,459 won classes: 938 of the 14,503 depth-3 classes and 8 of the 58 depth-7 classes have an exact depth below 2k - 1. For the 8,455 lost classes with d > 0 it is violated by 1,941. The cause has not been diagnosed. The hypothesis, not yet tested, is that holding the defender in place is not a relaxation, because the defender has to move and the forced move changes the defender's pieces.
- Gap (exact depth minus 2k - 1) by exact depth: depth 1, all 376,262 at gap 0; depth 3, 13,565 at 0 and 938 at -2; depth 5, 4 at 0 and 1,632 at 2; depth 7, 28 at 0, 20 at 2, 2 at 4, 8 at -2.
- PE holds in aggregate and fails once depth is matched. The 522 quiet-first classes: 482 at gap 2, 32 at 0, 8 at -2. The other won classes: 1,170 at gap 2, 2 at gap 4, 938 at -2, the rest at 0. At depth 5 alone, all 1,170 forcing-solved classes have gap 2, as do 462 of the 466 quiet-first classes. The gap follows the depth, and it does not separate the quiet-first classes from the forcing-solved classes of the same depth.
- Exposure variant: violated on 10,310 of 391,937 non-special won classes and on 8 of the 522, with no unreachable class. The denied edges are used by the shortest lines, so the denial is not a bound either.

**P-008 register at the lost classes (§8.3).** The defender's register where the attacker's rook gives check, replies partitioned as R_F, R_X, R_B:

| Group | Classes | Rook gives check | Flights per position | Captures | Interpositions | Held terms | Joint support (narrow and wide) |
|---|---|---|---|---|---|---|---|
| Seeds (lost d=4, quiet move needed) | 74 | 0 | none | none | none | none | none |
| Lost d=4, not seeds | 391 | 307 | 4.04 | 0 | 0 | F only, 307 | 0 |
| Lost d=6, forcing-solved | 100 | 99 | 5.09 | 0 | 0 | F only, 99 | 0 |
| Lost d=6, quiet move needed | 8 | 8 | 1.00 | 0 | 0 | F only, 8 | 0 |

- At the 74 seeds the defender is not in check, so C does not hold and the net on the king is not engaged. The register at a seed is empty by §8.3 as written. The attacker's win there needs a quiet move that establishes C, and the catalog route to that is the ENGAGED-NEW status of §8.4, whose engagement support is read through P-003 route distance (RT). How the engagement squares of a dormant net are defined is not written in §8.4.
- Where the rook gives check, only F is held, in every group. Narrow support = the unit, the squares it uses and the line it moves along; wide = narrow plus the coverage reads (attacker king, attacker rook and its ray). Under both readings the three families of replies never share support, and the one element common to all residual replies of the register is the defender's king, a key element of F.
- Not computed: effect vectors of defender moves (UNCHANGED, HELD, OPENED), the equilibrium moves of §8.11, and the opening mechanisms OP1 to OP6. They depend on the register at the positions the defender moves from, and the register is empty at the seeds .

**Not run.** PD (ordering by reach; withdrawn from this record as outside the scope of section 0) and PF (the support closure).

## 3C. Defence-form closure from the P-001 completion: predictions registered before the run (10 November, 20:50)

Sections 3C and 3D of an earlier version of this record (an attacker-candidate rule built on a support definition that the documents do not give) are withdrawn.

**What is run, from the documents only.** The defence form of §4.7: for the side to move, a move is invariant-preserving when, after it, the opponent has no move that completes a net at once (a move m with R(P, m) empty: for the king, the terms C, F, X, B of §4.3 hold, which is checkmate; for the rook, a capture the leaf of §4.3 O values as a win), and every other reply of the opponent leads to a position in which the side again has an invariant-preserving move. A capture of the opponent's rook by the mover is a move after which the opponent has no net. The set of positions in which the side to move has such a move at every turn, forever, is the greatest fixed point of that condition: the closure certificate of section 3.4 ("every position has a value-preserving move that stays in S"), computed with the primitive's completion test and no table. Both sides are treated alike: a class is lost for the mover when it leaves the fixed point, and the mover wins when some move leads to a position that left it (or completes at once). Positions with no legal move: in check, the mover is checkmated and is outside the set; not in check, a stalemate, inside.

**What it leaves out.** P-003 route distance, P-008 support and any pruning: every legal move and every reply is taken, so this is the exact defence form with no economy. It uses P-001 §4.3 and §4.7 and the closure certificate of section 3.4 only.

Predictions, stated before the run, all exact:
1. **G1.** The fixed point is reached after the removal rounds 0 to 3 and a round 4 that removes nothing. The classes removed in round j are exactly the classes with exact loss depth 2j (round 0: the 603 checkmated classes; round 1: the lost d=2 classes; round 2: lost d=4; round 3: lost d=6).
2. **G2.** The classes inside the fixed point from which the mover has no move to a position outside it (and no completing move) are exactly the drawn classes of the status set and the other drawn classes; the mover wins exactly in the 392,459 won classes. 0 disagreements with Syzygy over 1,347,906 classes in value.
3. **G3.** The depth of every win, read as one plus the round of the best move's child (a mate or a capture the leaf values as a win: depth 1), equals the exact depth in every won class.

## 3D. Defence-form closure: results (10 November, run after section 3C)

Computed over the 1,347,906 classes with the P-001 completion (terms C, F, X, B of §4.3 for the king, the capture leaf for the rook) and the defence form of §4.7, taking every legal move and every reply. No table value is read during the computation. The table is read afterwards, for the comparison only.

| Round | Classes removed from the closure |
|---|---|
| 0 (checkmated) | 603 |
| 1 | 7,882 |
| 2 | 465 |
| 3 | 108 |
| 4 | 0 |

- **G1** confirmed. The rounds remove 603, 7,882, 465 and 108 classes, which are the classes of exact loss depth 0, 2, 4 and 6 (the 465 are the 391 non-seeds and the 74 seeds; the 108 are the 100 forcing-solved and the 8 quiet-move classes). Round 4 removes nothing: the fixed point is reached.
- **G2** confirmed. Value agrees with Syzygy on 1,347,906 of 1,347,906 classes, 0 disagreements: won 392,459, lost 9,058, drawn 946,389. The drawn classes are the classes inside the fixed point from which the mover has no move to a position outside it.
- **G3** confirmed. Win depth equals the exact depth on all 392,459 won classes (0 differ), and the loss depth, twice the round, equals the exact depth on all 9,058 lost classes (0 differ).

What this establishes: the closure certificate of section 3.4 is supplied by the defence form with the P-001 completion test, over the whole segment, with no table and 100% agreement in value and depth. What it does not establish: any economy. Every move and reply was taken, so it is exact enumeration with the primitive's completion test, and P-003 and P-008 were not used. The classes the 74 seeds and the 8 are in are found only because the recursion looks two and three rounds past the immediate completion.

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

A first run of the python-chess checker disagreed on 94 samples at depth 5 and beyond. The cause was in the checker (its Lost ignored a mated position once d>0, and its Win skipped mate in one at d>=3). It was fixed and rerun. The C++ program was not changed.

## 5. Residual progress (where this run stopped)

Steps 2 and 3 of section 0, in order:
1. The defence form of §4.7 as an invariant over the status set: the closure certificate of section 3.4 supplied by the primitives (invariant-preserving moves, after which no net is COMPLETED for any attacker reply within h plies), tested against the 934,275 classes. Predictions first.
2. The residual reply sets R_F, R_X, R_B, R_O at every position, and the test whether the double-threat rule of §4.10 over them reproduces the rule of mode 3 and the ladder levels.
3. P-008 registers and equilibrium moves (§8.11) at the positions the defender moves from. At the 74 seeds no net is engaged and the register is empty (section 3B).

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
