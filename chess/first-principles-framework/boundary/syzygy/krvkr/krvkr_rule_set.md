# KRvKR: the primitive rule set for the segment (a data-area record)

**Status.** A rule set for one irreversible segment, written as executable conditions in the terms of the three primitives, with the segment's
interface to the segments below it, and checked against the exhaustive classification of the designer's tables on every position. It is a
data-area record and is not part of the programme's documentation. It uses the pieces of this material, which the programme's documents do not.
Sections 1 to 6 specify the calculation. Section 7 is the certification record, and every number in it is a measurement.

**What the calculation does.** Given any position of the segment with the side to move, it returns the status (the mover wins, the position
is drawn, the mover loses), the number of plies within which a decisive status is decided, the status of every legal move, the realization of
a draw (how it is sustained, with the line that shows it), and the boundary readings. It reads no table and enumerates nothing of the segment.
The tables are used only to check it (section 7).

**What is compared.** The comparison with the tables is of the value class of a position and of the class of each move: does the move keep the
value, or lose it. Moves of equal value are one class, and which of several equal moves is chosen is never compared. A difference would matter
only where a move or a position changes the value, which is where optimality is compromised.

**Files.** `segments/` holds one module per segment (`kvk.py`, `krvk.py`, `krvkr.py`), the registry that joins them (`material.py`) and the
checks (`check_segments.py`, `check_realization.py`). `krvkr_calc.cpp` is the compiled calculator with the harness that runs it on every
position. `export_arrays.py` writes the arrays the harness compares against. The mining scripts and the first results note are in the earlier
delivery. This record and `segments/krvkr.py` supersede the earlier `krvkr_rules.py`. The checks run from the `segments` folder, for example
`python3 check_segments.py <folder with the table files> 300` and `python3 check_realization.py <folder with pass1.npz and analysis.npz> 60`.

---

## 1. Scope

The material is a king and a rook against a king and a rook, no pawns and no castling. The only irreversible move in the material is a
capture, so the segment is the set of positions reached by reversible moves, and the capture is its boundary. The legal positions have
distinct squares, the kings not adjacent and the side not to move not in check. A position is read with the side to move.

## 2. The interface and the dependency chain

A segment reads a lower segment only through its status function, at the positions that its irreversible moves reach. Nothing else of the
lower segment is visible.

- **Requires:** KRvK, which every capture in this material reaches. KRvK is a segment of its own (`krvk.py`), with a closed-form rule that is
  certified by exhaustion (section 7, item 6). It in turn requires KvK, a dead position (`kvk.py`), status 0, which its own capture reaches.
- **Provides:** `status(position)` for the side to move, for every position of KRvKR. A segment above this one, whose capture reaches
  KRvKR, reads it by that function.
- **The chain:** KvK, KRvK, KRvKR, then the materials from which one capture reaches KRvKR. Each segment record has the same parts: the
  material and its boundary; the segments it requires; the net terms; the composition with its horizon and its pool; the derived statuses; the
  realization of a draw; the certification record; and what it provides.
- **What is material-specific and what is not in the code.** The composition (`W`, `L`, the pool, the drawing moves, the realization) names
  no piece. It reads the movement relation, the test for an irreversible move, and the status of a lower segment. A segment is added by
  registering its module under its material signature. The one place that changes for a material with pawns is the test for an irreversible
  move, which is a capture here and would include a pawn move there. That the composition carries over unchanged is tested by writing the
  next segment, and it is not claimed here.

## 3. The net terms (P-001), for the two kinds of target

**Royal target, mate in 1.** The mover has a move after which the defender is checkmated: the defender's king is attacked (C), no square
in its neighbourhood is free of occupancy and coverage (F), the checking rook cannot be captured by the defender's king or rook (X), and
nothing can be placed between the rook and the king (B). In legal-move terms: the defender is in check and has no legal move.

**Non-royal target, the defender's rook.** The mover has a legal capture of it (C). The outcome O is read through the interface: it is the
status of the lower segment at the position that the capture reaches. In KRvK that is a win for the capturer unless the defending king can
capture the capturing side's rook (the net on the capturing unit completes in return) or is stalemated, and then it is a draw.

## 4. The boundary and the exit

A capture is the irreversible move, and across it the goal changes. Before it the goal is the net on the rook. After a favourable capture the
goal is the mating net on the king in the remaining material, decided there. After an unfavourable capture the position is drawn. A drawing
exit is a resource for the side that has it, and a win has to exclude every such exit (the exit condition). In the composition below, a
capture is valued by the status that the interface returns and by nothing else.

## 5. The composition

Let h be a horizon in plies, odd for a win and even for a loss.

- **W(h), the mover wins within h plies.** The mover has a capture whose exit value is a win, or W(h-2), or a move m from the pool after
  which the defender loses within h-1 plies. For h = 1 the last alternative is a mate in 1.
- **L(h), the mover loses within h plies.** h = 0: the mover is checkmated. Otherwise L(h-2), or: the mover has legal moves and each is
  refuted, a capture by its exit value being a loss and any other move by the opponent winning within h-1 plies. When every move is refuted
  by a net, the zugzwang-type composition of the three primitives holds.
- **Status.** The mover wins if W(7). The mover loses if L(6). Otherwise the position is drawn: the draw is the complement.
- **The pool.** The candidate moves of the side that is trying to win are the checks, and the moves after which the mover would win
  within h-2 plies if the defender passed (a threat that a net is completable). The defender is not restricted: all its legal moves are examined.

The horizon is a parameter of the calculation (primitive catalog input h). 7 and 6 are the values at which this segment is exact: the
exhaustive classification contains no decisive position beyond them.

## 6. The derived statuses

**A move.** A move is classed by the status it leads to: it loses if after it the opponent wins within 7 plies, it wins if after it the
opponent loses within 6, and otherwise it draws. A capture is classed by its exit value.

**The boundary readings of a drawn position** (main file §14 items 7 and 8 read three ways at the position, and once over the region):
A, the mover has a losing move; B, A holds and some drawing move that stays in the material leads to a position where the opponent has a
losing move; C, A holds and every such move does. The region reading: the closure of the drawing moves of both sides contains a position where
the first side to move has a losing move and one where the other side has.

**The realization of a draw** (the repetition cycle, and the forced terminal draw of item 10).
- *The drawing graph:* the nodes are the drawn positions, the edges the drawing moves that stay in the material. A drawing exit is a terminal
  draw.
- *The cycle, as a local certificate:* a four-ply back-and-forth line through the position, in which each side moves a piece and moves it
  back, every position drawn. If the position has none, the nearest position along drawing moves that has one is found, and the path to it is
  the certificate. (The closure is searched forward.)
- *Forced terminal:* the closure of the drawing moves from the position is exhausted with no cycle in it, so every drawing line ends in an exit.
- *The classes:* terminal only, cycle only, mixed (both a cycle and an exit are reachable), stalemate now.

## 7. The certification record (measurements)

1. **The rules against the exhaustive classification, every position.** 1,347,906 folded legal positions (16 images). The calculation
   gives the same value class on all of them: 392,459 wins, 946,389 draws, 9,058 losses. It gives the same number of legal moves for all, and
   at the 946,389 drawn positions the same number of losing moves, drawing moves that stay in the material and drawing exits. The
   calculation with the pool of all moves gives the same result.
2. **The horizon.** The classification contains no win decided in more than 7 plies and no loss in more than 6.
3. **The refinement of the pool, with the differences in value class counted at each pass.**

   | Pool of the winning side | Positions of 1,347,906 whose value class differs |
   |---|---|
   | all moves | 0 |
   | checks | 816 (803 wins and 13 losses missed) |
   | checks and moves that attack the rook (the rook sees it, with the king as the only blocker, or the king stands next to it) | 530 (522 wins and 8 losses missed) |
   | checks and threats of a one-move capture or mate | 8 (wins missed) |
   | checks and threats of a win within h-2 plies | 0 |

   The 8 positions of the fourth row are one family: the winning side's first move is a quiet king move that makes no check and no threat of
   a capture or mate, and the win is decided in 5 plies. They are closed by taking the threat to the horizon, which is the last row. The pool
   is part of this segment's rule set, and it is certified by the exhaustion above.
4. **The cost of one position, from a cold start.** For a drawn position the number of move generations has a mean of 5,888, a median of
   5,700 and a maximum of 11,230 (the pool of all moves: mean 10,732). For a winning position the mean is 5 and the median 1. The whole segment,
   run once with a shared memo, takes 5.4 seconds compiled. The exhaustive run is the check and is not needed to calculate a position.
5. **The readable version.** The modular `segments/krvkr.py`, with its exits read through the interface, has the same value class as the
   table on 300 random positions (147 wins, 150 draws, 3 losses). The earlier single-file version had the same on 800.
6. **The segments below.** The `krvk.py` rule (the rook side to move wins; the bare king to move draws if it can capture the rook or is
   stalemated, and otherwise loses) has the same value class as the KRvK table on all 399,112 legal positions, both sides to move.
7. **The tempo of the board.** For each board in which neither king is in check (8,530,032 boards, not folded), the result read from White's
   point of view with each side to move:

   | White to move, Black to move | White loses | drawn | White wins |
   |---|---|---|---|
   | White loses | 1,560 | 0 | 0 |
   | drawn | 685,020 | 5,486,740 | 0 |
   | White wins | 1,670,132 | 685,020 | 1,560 |

   The result of 3,040,172 boards (35.64%) depends on the side to move. In 1,670,132 of them whoever has the move wins. In the other 1,370,040 one
   side wins if it has the move, and the board is drawn if the other side has it. In no board is the result better for the side without the move
   than for the side with it, so having the move is never a disadvantage in this material.
8. **The exit as a resource.** 7,576 drawn positions (0.8%) are drawn only because of an exit: 3,771 would be won and 3,805 lost without it.
9. **The drawing graph, from the rules.** 946,389 drawn positions and 11,988,840 drawing edges within the material. One strongly connected
   component has more than one position, with 933,666 positions (98.66% of the drawn). 12,723 drawn positions are on no cycle. The classes are
   the same as those of the exhaustive classification: terminal only 2,857, cycle only 0, mixed 943,532, stalemate now 0.
10. **The local certificate of the cycle.** 932,206 drawn positions (98.50%) have a four-ply back-and-forth cycle through them. Of the rest,
    11,287 reach a position that has one in 1 ply, 38 in 2 plies and 1 in 3. The 2,857 left are exactly the terminal-only positions, from which
    no cycle is reachable. `realization()` agrees with the classes on 60 sampled terminal-only positions (closure of at most 2 positions) and 60
    sampled mixed positions.
11. **The region reading of the boundary.** 946,272 drawn positions (99.99% of the drawn) have, in the closure of their drawing moves, a losing
    deviation for each side. All 943,532 mixed positions do. Of the 2,857 terminal-only positions 2,740 do, and 117 do not. The readings at the
    position itself (A, B, C) are in the first results note and give 86.9%, 86.5% and 26.7% of the drawn.

## 8. What the record does not cover

- The mating net of the remaining material is taken as the decided status of the KRvK segment (item 6). It is not derived in the primitive terms here.
- The pool and the threat test are stated as moves and null-move tests. Writing them in the held terms and the tempo terms of P-003 and P-008 is
  the designer's step. The record uses the coverage that the legality of the moves gives and no separate coverage calculation.
- Which of the boundary readings is the designer's definition is the designer's decision. The record gives the counts of each.
