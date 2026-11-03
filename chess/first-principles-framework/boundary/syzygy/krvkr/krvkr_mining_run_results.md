# KRvKR mining run: results (a record of one run, kept)

**Status.** Results of running the process of `boundary_board_states_syzygy_mining.md` §3 on the designer's KRvKR and KRvK tables
(`.rtbw` and `.rtbz`, read with python-chess 1.11.2). Everything below is a measurement. Readings are marked as readings. The scripts are
`krvkr_mine1_enumerate_classify.py`, `krvkr_mine2_graph_classes.py`, `krvkr_mine3_retrograde_counterfactual.py`, and a run takes about
ten minutes on two cores. No document is changed by this file. Whether any of it is entered in the documents is the designer's decision.

## 1. The process ran, and the checks hold (§3 steps 1 to 5 and 8)

| Item | Result |
|---|---|
| Raw positions with distinct squares | 30,498,048 |
| Folded representatives (16 images: 8 board symmetries x colour swap) | 1,906,548 |
| Legal by the file's rule (side not to move not in check, so kings not adjacent) | 1,347,906 folded |
| Mover wins / draws / loses | 392,459 / 946,389 (70.2%) / 9,058 |
| Step 5: value equals the best child, negated, over every legal position | 0 mismatches in 1,347,906 |
| Independent retrograde solve from the move structure and the KRvK exit values only, not reading the KRvKR values | agrees on all 1,347,906 positions |
| Second move generator (own king and rook code), 20,000 positions, move and capture counts | 0 disagreements |
| Fold check: 14,181 random unfolded positions against their representative | 0 disagreements |
| Cursed or blessed values (1 or -1) present | none: every value is 2, 0 or -2 |
| Drawn positions with a winning move | 0 |

## 2. The drawing graph (§3 step 6)

- 14,507,570 drawing edges that stay in the material. One strongly connected component holds 933,666 drawn positions, which is 98.66% of
  the drawn positions. 12,723 drawn positions are on no cycle. There are 31 weak components, and the largest has 946,359 positions.
- Realization classes of the drawn positions: stalemate now 0; terminal only 2,857; cycle only 0; mixed 943,532.
- Terminal-only positions: 2,791 have one drawing move and 66 have two. In 2,776 of them the mover is in check. 2,792 have no drawing
  move that stays in the material (the draw is a capture exit). A typical one is a position where the mover must take the checking rook.
- Drawing moves that capture: 35,946 drawn positions have at least one, and 99 of them reach a stalemate.
- Width (drawing moves per drawn position): mean 12.7, maximum 21. Exactly one drawing move: 8,060 positions (0.85% of the drawn).

## 3. The boundary criterion (main §14 items 7, 8), on three readings of "each side has a losing deviation"

| Reading | Drawn positions | Share of drawn | Share of all legal |
|---|---|---|---|
| A: the mover has a losing move | 822,408 | 86.9% | 61.0% |
| B: A, and some drawing in-material move reaches a position where the opponent has a losing move | 819,023 | 86.5% | 60.8% |
| C: A, and every drawing in-material move reaches such a position | 252,994 | 26.7% | 18.8% |

Drawn positions with no losing move at all (the interior of mining §2): 123,981 (13.1% of the drawn).

## 4. The force region is shallow (the complement, start-selection §2c)

Race depth is the number of retrograde rounds to resolve a position, which is plies. For the won positions it equals the table's DTZ
on a sample of every seventh decisive position (57,360 positions; correlation 1.0, values 1, 3, 5, 7). In KRvKR the zeroing move is the
capture of a rook, so for a win DTZ is the length in plies of the race to the capture. For lost positions DTZ is one less than the
retrograde depth, and the mated positions (depth 0) have no DTZ.

- Mover wins: 95.9% are won in 1 ply (almost all by a capture that the table values as won), 3.7% in 3 plies, 0.4% in 5, 0.015% in 7. The maximum is 7.
- Mover loses: 6.7% at depth 0 (mated or equivalent), 87.0% at 2 plies, 5.1% at 4, 1.2% at 6.
- So the whole not-drawn region, 29.8% of the legal positions, is decided by capture races of at most 7 plies.

## 5. The exit resource (start-selection §2c, the exit condition)

Counterfactual: every capture exit whose table value is a draw is scored as a loss for the side that would take it, and the game is
re-solved (iteration over the full move graph).

| Baseline | Counterfactual | Positions |
|---|---|---|
| win | win | 392,459 |
| loss | loss | 9,058 |
| draw | draw | 938,813 |
| draw | loss for the mover | 3,805 |
| draw | win for the mover | 3,771 |

So 7,576 drawn positions (0.80% of the drawn) are drawn only because of an exit: 3,771 would be wins for the mover without the
defender's exit, and 3,805 would be losses without the mover's own. No won or lost position changes status.

## 6. Readings (this file's, offered for correction)

1. The mining process of §3 is executable as written and is internally verified; the generator and the tables agree.
2. In this material the draw set is one connected web (98.66% in a single component). The calibration signal of survey §7, that no
   candidate's verification can be separated from the rest, is the case here. The certificate for a draw in this material is the whole
   drawn region, and any compression has to come from a rule set for the complement, which section 4 shows to be shallow and small.
3. The designer's expectation that the information lies in the complement (start-selection §2c) is borne out in size: 29.8% of the
   positions, all decided in at most 7 plies. Whether a primitive rule set reproduces it without search is not tested here.
4. The exit condition does work, and it is small in this material: it is the only thing between 3,771 positions and a win for the
   mover. That is 0.4% of the drawn positions, so the measured effect of the condition here is rare, and it is a necessary condition on
   a win, as the statement says. The programme's version of it (stated in primitive terms) is not tested here.
5. On the boundary share (mining §1: whether it is a majority of the space): under the weakest reading it is, 61% of all legal
   positions and 87% of the drawn. Under reading C it is not. The criterion alone does not separate much of the drawn region in this
   material, so the realization and the width carry the information.
6. The two empty realization classes (stalemate now, cycle only) are properties of this material.
7. One material. None of this says what happens for the same measurements elsewhere.
