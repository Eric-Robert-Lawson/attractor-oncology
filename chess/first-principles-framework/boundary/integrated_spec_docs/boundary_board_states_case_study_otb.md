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
through the player's own topological narrative. It belongs beside the other two as supporting material,
not inside them. Where it touches the formal program directly (one candidate
new concept, see the end of this document), that is flagged for a decision
about whether to fold it into `boundary_board_state_topology.md` §8g.

**Method, stated plainly.** The game record and the narrative are the author's.
Nothing in this document is machine-checked: no engine evaluation, engine-suggested
move or enumeration is kept here. A reading that rests on inspecting the board is
marked as such, and the assistant's readings are separate from the author's.

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

The game ends in checkmate on the board, as recorded by the author.

**Author's own framing, stated up front and kept as the framing throughout.**
Not claimed as perfect play — the author is explicit that it is not. The
claim is narrower: that a topological, boundary-board-state-informed reading
of the position, developed without formal chess training beyond one opening
book, produced strong, coherent play, including at least one moment where
the topological read aligned with a concrete tactical necessity (a missed
mate threat) that the author was not consciously calculating at all.

---

## 2. The moments the narrative identifies

The author's narrative places the significant moments of the game at White's early
queen commitment (Qf3), at Bd3, and at the exf4 and Kf1 pair around the final
combination. The author describes the early queen move as costing against
topological equilibrium balance: it commits the queen before development and
imposes a race on the position. White's position then erodes in discrete steps that
track the author's own segmentation of the game into topologically significant
moments. This is the author's reading. No engine data is kept here to support or
contradict it.

---

## 3. The mate-in-1 threat

The claim is the most concretely falsifiable one made. After White's 9th move
(Qe4), the author's reply (f5) was chosen for reasons unrelated to any mate threat,
described as continuing a plan to attack the queen and pursue development. Only in
hindsight did the author realize this move also prevented an unnoticed mate threat.

**The reading, by inspection of the board.** From the position after 9. Qe4, if it
were White's move, Qxh7# would be checkmate. The queen on e4 sits on the same
diagonal as h7 (e4-f5-g6-h7), the diagonal is open, and the black king on g8 has no
flight square (f7 and g7 are occupied by Black's own pawns, h8 is covered by the
queen along the h-file), with no piece able to capture or block on h7. The author
was, by the author's own account, not calculating this at all.

**What makes f5 notable.** It does two separate jobs in one move. It blocks the
e4-h7 diagonal (removing the mate), and it attacks the White queen directly (forcing
the queen to move again, which is the "queen must move, race condition imposed"
reading the author gives the position). The move was chosen for the second reason,
in the author's own account, and the first reason held anyway. How many of Black's
other replies also met the threat is not recorded here.

**Assessment.** If the reading holds, it is evidence that the
topological framework is tracking something outside itself rather than
constructing an internally-consistent story after the fact. It is not proof
that the framework can be formalized into the exhaustive, composable catalog
described in `boundary_board_state_topology.md` §8 — that is a different,
much harder claim, addressed in §6 below. It is evidence that the intuitive
version of the framework, applied under real OTB conditions with a ten-second
or so decision window, correctly favored a move that satisfied a real
tactical necessity the player was not consciously tracking, for reasons that
were themselves independently sound (tempo, continued development).

---

## 4. Other tactical claims

The author's reads, from the game record:

- **Nxd3+ (ply 20).** The author's read is a free capture, delivered with check.
- **Nxf4+ (ply 28).** The knight lands on f4 and attacks e2, where the White king
  stood, and g2, where the White queen stood: a fork on king and queen, from the
  knight's move pattern and the recorded positions.
- **The White knight on g1.** In the game record it never moves until it is captured
  on move 16 (Bxg1). The author's read, "white knight virtually useless", matches
  the record.

---

## 5. What is not kept here

Engine evaluations of the game, the author's recollections of engine-suggested best
moves (chess.com's suggestions), and counts of defences to the Qxh7# threat are not
part of this record. They were produced by the assistant's runs and are not the
author's data.

---

## 6. What this does, and does not, establish

**What the record supports.** The author's narrative places the cost of early,
undeveloped queen commitment where the author says it falls. The mate-in-1 reading
shows that an intuitive read, arrived at through the author's own topological
vocabulary and not through conscious calculation of the specific threat, landed on
a real tactical necessity. Several other tactical claims (the capture, the fork, the
wasted knight) read as the author describes. This is evidence that the intuitive
framework has practical value as a heuristic for actual play, produced independent
of, and prior to, any of the formal apparatus in the other two documents.

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
what was real in this game,
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

**Second instance (see `boundary_board_state_topology.md` §8u).**
In a second chess.com puzzle, the move ...g4 attacks a pinned queen, protects the
square h3 on which the mating rook lands, and fills the flight square g4, all at
once. It differs from f5 in one respect that matters for the question above: f5
aligned by coincidence, g4 is chosen for the combination. The structure, one move
discharging several duties with no conflict between them, has therefore now
appeared twice. Whether it is a named pattern that adds something to the
composition vocabulary, or an ordinary tempo move that builds the mating net,
is still undecided and is left for that decision.
