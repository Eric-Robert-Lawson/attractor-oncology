# Boundary board states in practice — a worked OTB case study

A new, separate document, created at the author's request after presenting a
real game (Black, vs. a chess.com bot rated ~1600) as evidence that the
topological framework developed in `boundary_board_states.md` and
`boundary_board_state_topology.md` is not purely tautological — that it
produces real strength in actual play, independent of formal training.

**Why this is a separate document rather than folded into the other two.**
Those two documents are a formal research program: definitions, tablebase-
grounded claims, and a construction methodology aimed at eventual proof. This
document is a different kind of object — a single game, played under
imperfect conditions by a self-described informally-trained player, examined
both through the player's own topological narrative and through independent
engine verification. It belongs beside the other two as supporting material,
not inside them. Where it touches the formal program directly (one candidate
new concept, see the end of this document), that is flagged for a decision
about whether to fold it into `boundary_board_state_topology.md` §8g.

**Method, stated plainly.** The game was replayed move by move with
`python-chess` and independently analyzed with Stockfish 16 (depth 18-20,
some positions to depth 20 with a 5-line MultiPV search, run in this
project's own sandbox — not chess.com's engine, which is not accessible
here). Every specific, checkable claim below was checked against that
independent analysis rather than accepted from memory or from the narrative
alone, in keeping with this project's standing practice of separating record
from assessment. Claims that reproduced are marked confirmed with the actual
numbers. Claims that did not reproduce are reported as such, plainly, with
what the independent analysis found instead.

---

## 1. The game, as played

Black (the author) vs. chess.com bot (~1600), White.

```
1. c4 e5           2. Nc3 Nc6          3. g3 Bc5           4. e3 Nf6
5. Nd5 O-O          6. Qf3 Nxd5         7. Qxd5 d6          8. Bd3 Nb4
9. Qe4 f5          10. Qf3 Nxd3+       11. Ke2 e4          12. Qg2 Qg5
13. f3 f4          14. exf4 Nxf4+      15. Kf1 Nxg2        16. Kxg2 Bxg1
17. Rxg1 exf3      18. Kh1 b6          19. h4 Qg4          20. Re1 f2
21. b4 fxe1=Q+     22. Kh2 Rf2#
```

Verified: fully legal from the starting position through to an actual
checkmate on the board (`python-chess` confirms `board.is_checkmate() ==
True` at the final position, FEN
`r1b3k1/p1p3pp/1p1p4/8/1PP3qP/6P1/P2P1r1K/R1B1q3 w - - 2 23`). This was not
assumed from the narrative; it was checked.

**Author's own framing, stated up front and kept as the framing throughout.**
Not claimed as perfect play — the author is explicit that it is not. The
claim is narrower: that a topological, boundary-board-state-informed reading
of the position, developed without formal chess training beyond one opening
book, produced strong, coherent play, including at least one moment where
the topological read aligned with a concrete tactical necessity (a missed
mate threat) that the author was not consciously calculating at all.

---

## 2. The evaluation trajectory — what the engine actually shows

The single clearest, least interpretive fact in this whole game: White's
position collapses in large, discrete jumps, and the jumps line up with
exactly the moments the author's narrative identifies as significant. This
table is Stockfish's own centipawn evaluation (positive = better for White)
after each ply, not the author's assessment:

| Ply | Move | Eval (cp, +White) | Note |
|---|---|---|---|
| 8 | Nf6 | +47 | Roughly balanced, slightly better for White |
| 9 | Nd5 | −33 | First dip |
| 10 | O-O | −10 | |
| **11** | **Qf3** | **−178** | **Queen commits early — −168cp in one move** |
| 12 | Nxd5 | −197 | |
| 13 | Qxd5 | −199 | |
| 14 | d6 | −194 | |
| **15** | **Bd3** | **−562** | **−368cp in one move** |
| 16 | Nb4 | −578 | |
| 17 | Qe4 | −682 | |
| 18 | f5 | −683 | |
| 19 | Qf3 (retreat) | −756 | |
| 20 | Nxd3+ | −779 | |
| 21 | Ke2 | −764 | |
| 22 | e4 | −763 | |
| 23 | Qg2 | −749 | |
| 24 | Qg5 | −756 | |
| 25 | f3 | −761 | |
| 26 | f4 | −770 | |
| **27** | **exf4** | **−1472** | **−702cp in one move** |
| 28 | Nxf4+ | −1572 | |
| **29** | **Kf1** | **−2150** | **−578cp in one move** |
| 30 | Nxg2 | mate in 15 | From here the position is a forced mate |

**Assessment.** The single largest, cleanest fact here is the first jump:
White's position is roughly balanced through move 8, and the moment the
queen comes out on move 6 (ply 11 in this table, since ply counts every
half-move from 1), the evaluation drops by 168 centipawns in one move, before
any tactic has occurred. This is an unambiguous, independently-verified
confirmation of the author's core read at exactly the moment the author
identified it — ply 6 (White's move count), described in the narrative as
"the causal reason for the move ... costs against topological equilibrium
balance." The engine agrees, in its own vocabulary, before any of the
subsequent tactics happen. Two further large jumps follow at the two other
moments the narrative flags as decisive errors (Bd3 and, later, the exf4/Kf1
pair around the final combination). The overall shape of the game — White's
position eroding in large, discrete steps that track the author's own
segmentation of the game into "topologically significant" moments — is real
and confirmed, not a retrospective story fitted to the result.

---

## 3. The mate-in-1 threat — checked directly, and it is exactly as described

This is the single most important finding in this document, because it is
the most concretely falsifiable claim made, and it survived the check.

**The claim, as stated.** After White's 9th move (Qe4), the author's reply
(f5) was chosen for reasons unrelated to any mate threat — described as
continuing a plan to attack the queen and pursue development — and only in
hindsight did the author realize this move happened to also prevent an
unnoticed mate threat.

**What was checked.** The exact position after 9. Qe4 (Black to move) was
taken and examined for whether White had an immediate mate available, by
testing every legal White move from that same piece arrangement for
checkmate.

**Result: confirmed, precisely.** From the position after 9. Qe4, if it were
White's move, **Qxh7# is checkmate.** The queen on e4 sits on the same
diagonal as h7 (e4-f5-g6-h7), the diagonal is completely open, and the black
king on g8 has no flight square (f7 and g7 are occupied by Black's own pawns,
h8 is covered by the queen along the h-file) and no piece can capture or
block on h7. This is real, forced, and was on the board the moment White
played Qe4.

**Was f5 the only defense, or one of several?** Checked exhaustively — every
one of Black's 34 legal replies to Qe4 was tested. **8 moves genuinely
neutralize the threat**: Re8, Qh4, Bf5, Nxd3+, Nc2+, g6, f6, and **f5**. The
other 26 legal moves all still hang the mate. So this was not a
one-in-thirty-four needle — a real, if minority, share of moves work — but it
also was not a trivial or obvious defensive requirement, and the author was,
by the author's own account, not calculating it at all.

**What makes f5 specifically notable among the eight.** It is not simply one
of the defenses; it does two separate jobs in one move. It blocks the
e4-h7 diagonal (removing the mate), and it attacks the White queen directly
along the newly-created f5-e4 relationship (forcing the queen to move again,
which is exactly the "queen must move, race condition imposed" reading the
author gives the position). The move was chosen for the second reason, in the
author's own account, and the first reason held anyway.

**Assessment.** This is real, and it is the strongest single piece of
evidence produced in this whole project, across every session, that the
topological framework is tracking something outside itself rather than
constructing an internally-consistent story after the fact. It is not proof
that the framework can be formalized into the exhaustive, composable catalog
described in `boundary_board_state_topology.md` §8 — that is a different,
much harder claim, addressed in §6 below. It is proof that the intuitive
version of the framework, applied under real OTB conditions with a ten-second
or so decision window, correctly favored a move that satisfied a real
tactical necessity the player was not consciously tracking, for reasons that
were themselves independently sound (tempo, continued development).

---

## 4. Other tactical claims checked

**The Nb4-to-d3 capture ("free bishop," ply 20, Nxd3+).** Checked: the bishop
on d3 had zero defenders at the moment of capture. Confirmed as a genuinely
free capture, delivered with check.

**The Nxf4+ fork (ply 28).** Checked: at the moment this move is played, the
knight lands on f4 and simultaneously attacks e2 (where the White king
stood, confirmed check) and g2 (where the White queen stood, confirmed by
piece location in the position). A real, simultaneous fork on king and queen,
not simply a capture with a follow-up threat.

**The White knight on g1.** Checked across the FEN record: White's g1 knight
never moves for the entire game, until it is captured directly on move 16
(Bxg1). The author's read — "white knight virtually useless" — is confirmed
literally: it made zero moves in a 22-move game before being captured for
free.

---

## 5. Claims that did not reproduce — reported plainly

In the interest of the same honesty this project has applied to its own
formal claims throughout, three specific recollections of engine-suggested
best moves were checked directly and did not reproduce as stated. This does
not affect anything in §§2-4 above, which rest on independent findings, not
on cross-checking the author's memory of chess.com's output.

- **"Best move suggested was Nf3 instead" (White's 9th move, i.e. the
  alternative to Qe4).** Checked with a 5-line search at the position before
  White's 9th move: the top five candidates found were Bxh7+, Qg2, h4, Ne2,
  and b3. Nf3 did not appear among them. Worth noting: White's actual move,
  Qe4, also does not appear among the top five — so the broader point (White
  had better options than committing the queen further) still holds, even
  though the specific move recalled does not match this independent check.
- **"Better move was Kd1" (near White's move 13, f3).** Checked at the
  position before White's 13th move: the top five were f3 (the move actually
  played), Rb1, Kf1, h4, a4. Kd1 was not among them, and the move actually
  played was, in this independent check, the top choice — the opposite of
  what was recalled.
- **"Kd1 still best" (near White's exf4 capture on move 14).** Checked at the
  position before that capture: the top five were fxe4 (a different capture,
  taking the other Black pawn), Nh3, Qf1, Kd1, b3. Kd1 does appear here, in
  fourth place, a real candidate but not the top choice — the top choice was
  a different move entirely (fxe4), and the move White actually played
  (exf4) does not appear in the top five at all, confirmed as a large error
  independent of which alternative is named.
- **A closer match, found while checking the above.** Right after Nxd3+
  (White's 10th move, forced to respond to check), the top three replies
  were nearly tied: Kf1 (−765), Ke2 (−773), Kd1 (−773). White played Ke2.
  Kd1 and Ke2 are here essentially equivalent, both slightly behind Kf1. This
  is the closest match found to the "Kd1 was best" recollection, though even
  here it is tied with the move actually played rather than clearly superior
  to it.

**Assessment.** These discrepancies most likely come from chess.com using a
different engine, version, depth, or move-highlighting convention than the
Stockfish 16 / depth 18-20 search run here, or from the recollections being
approximate (the game was recalled and narrated some time after being
played, not read from a saved log). None of this affects the confirmed
findings in §§2-4, which were derived independently rather than by
cross-checking memory. It is recorded for the same reason this project
records every discrepancy it finds in its own formal claims: because a
record that only keeps the parts that confirmed the thesis is not a record
worth trusting.

---

## 6. What this does, and does not, establish

**What is confirmed, plainly.** The evaluation trajectory shows that early,
undeveloped queen commitment carries a real, large, immediately-measurable
cost in this position, exactly where the author's narrative places it. The
mate-in-1 finding shows that an intuitive read — arrived at through the
author's own topological vocabulary, not through conscious calculation of the
specific threat — landed on a real tactical necessity. Several other specific
tactical claims (the free bishop, the fork, the wasted knight) hold up
exactly as described. This is real evidence that the intuitive framework
has practical value as a heuristic for actual play, produced independent of,
and prior to, any of the formal apparatus in the other two documents.

**What is not established, and should not be overstated.** Punishing early
queen development with tempo-gaining developing moves is one of the oldest,
most widely taught principles in chess — not a new discovery. What is
distinctive here is not the underlying chess content, which is classical, but
the author's own framework for holding and applying it: arrived at
independently, described in a topological vocabulary the author built for
their own reasons, and applied with more consistency across the whole game
than the author's stated formal background would predict. That is a genuine
and interesting finding about the framework's practical value as a personal
heuristic. It is a different claim from the formal program in
`boundary_board_state_topology.md` §§8f-8x, which is that a boundary state is a composition of
primitive race conditions recorded in a catalog. That follows from the definition of a boundary
state (`boundary_board_states.md` §14 item 8) and does not rest on this game; one game neither
supports it nor is needed for it. What the game does show is that the intuition being formalized
is tracking something real, which bears on which primitives to catalog first. Building the
composition calculus is the engineering work of the primitive catalog
(`boundary_board_state_primitive_catalog.md`).

**On the single-game 2150 performance figure.** Worth naming plainly: a
performance rating from one game is a very noisy statistic. A single strong
positional-plus-tactical theme, well executed against a weaker opponent's
compounding errors, can easily produce an inflated single-game number that
would not hold up as an average across many games. This does not diminish
what was real in this game — the confirmed findings above stand regardless —
but the 2150 figure itself should be read as a single-game estimate, not a
claim about standing playing strength.

---

## 7. One candidate new concept, flagged for a decision

The f5 move in §3 does two jobs at once — defends against a mate and gains a
tempo on the queen — at zero additional cost. This is worth naming precisely
because it is the *mirror image* of overload, the composition-failure mode
already recorded in `boundary_board_state_topology.md` §8f: overload is a
single piece unable to serve two duties, which an attacker exploits to break
a composition. What happened here is a single move serving two duties for
free, with no conflict between them. If this pattern recurs — a move whose
effect on one race condition happens to also satisfy a separate one, without
being chosen for that reason — it would be a genuinely useful addition to the
composition vocabulary in §8g: a best-case resource-sharing event, as opposed
to the worst-case (overload) already recorded there. Not added to that
document yet, since this is one instance, and whether it is a repeatable,
identifiable pattern rather than a single lucky alignment is itself an open
question. Flagged here for a decision on whether to fold it into §8g as its
own subsection once (or if) more instances turn up, in a future game or in
the tablebase-derived population once §8m's scan is run.

**Second instance, added later (see `boundary_board_state_topology.md` §8u).**
In a second chess.com puzzle, the move ...g4 attacks a pinned queen, protects the
square h3 on which the mating rook lands, and fills the flight square g4, all at
once. It differs from f5 in one respect that matters for the question above: f5
aligned by coincidence, g4 is chosen for the combination. The structure, one move
discharging several duties with no conflict between them, has therefore now
appeared twice. Whether it is a named pattern that adds something to the
composition vocabulary, or an ordinary tempo move that builds the mating net,
is still undecided and is left for that decision.
