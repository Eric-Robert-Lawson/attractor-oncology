# Boundary Board States: Where to Start (the discovery of the first line to survey)

**Status.** A record of the discovery phase, not a method and not a result. It holds the designer's statements about which line the
anchored opening survey (`boundary_board_states_opening_survey.md`) should start with, and what has and has not been examined. It
decides nothing. Nothing in it has been run on a pruned structure.

**Companion documents.** `boundary_board_states_opening_survey.md` (the method; §4e is the statement this file expands, §7 is the
calibration of a line and a depth), `boundary_board_state_topology.md` (§8aj), `boundary_board_states.md` (changelog entries 38, 39
and 42), `boundary_board_state_primitive_catalog.md` (the three primitives). The designer's statements are marked as such and kept
apart from this file's readings, as in the catalog's §10.

---

## 1. The state of the discovery

The designer is still in the discovery phase: understanding the program well enough to deduce where to start. The question of where
to start is therefore open on the White side and scoped on the Black side, as below. It is not a settled plan.

## 2. The designer's statements

**What was scoped.** From high-level observation, where play is closest to optimal, over the board and in timed variants, the Black
schematic described in the survey §4e has emerged, with three specific moves (g6, Bg7, d6). The designer identified these moves
uniquely with the primitives from the initial position as being positioned to be a rich vein for this process: if a downstream
boundary board state position can be found, the reachability should be rich. The designer has therefore scoped the Black moves
properly for where the survey would start.

**What was gone over in detail.** The Black schematic. The detail of why, from the primitives, is in the survey §4e ("Why, from the
primitives"): development race progress, no negative impact on the other race condition primitives, no commitment to a branching
path, the mobility unlocked by d6, and a small pool of deviations that threaten.

**What has not been looked at.** White's position as the constructor, other than Bird's Opening, which came from another plan.

**White's side, begun: the worked example narrowed to the central pawns (a statement).** The operation and its first steps, down to pawn
moves as the only first moves that force a deviation from the plan, are in the survey §4i. This continues that example.

Within two moves, pawn moves on the a, b, c, f, g and h files are not as impactful as the e and d pawns in a way that would obstruct the
development of the plan inside the window. The scope is therefore the combinations of e and d pawn moves: ee (e4, then e5), ed (e4,
then d4), de (d4, then e4) and dd (d4, then d5).

This is why d6 is in the plan: any e or d pawn push after g6 can be followed by d6. From there it is a combinatorial tree.

- **dd.** Black completes the plan; the d pawn is blocked by the d6 pawn, and d6 covers a knight development square as a primitive
  impact. White may build a wall of pawns, but the wall is on light squares, so the dark squares are vulnerable to the dark-squared
  bishop on g7.
- **ed and de.** After the bishop move to g7, White pushes either the e pawn or the d pawn. The d push is an extension of the dd branch.
  The e push offers a trade that would advantage Black, one pawn up with no White central pawns. The move sequence of this line was
  not given.
- **ee.** e4, g6, e5, d6. Black may lose this race. If so, against the first move e4, d6 as Black's first move is forced. That makes
  e5 an important square for Black to contest, as a coverage primitive to calculate in this instance.

Deviation into back-rank piece moves or the other pawn pushes is, temporally, a combination of a different opening in which those
move-order choices could have been switched with no impact on the plan. It is excluded in the same way as the back-rank moves and the
a, b, c, f, g and h pawn pushes.

If de or ed is played, with any deviation moves, Black can follow the Sicilian fianchetto variant: c5 opens the queen's diagonal, and
the c pawn is pushed past the knight development point, which puts pressure on White's d pawn push and prevents that primitive
development race. The sequence is d4, g6, e4, d6, a knight or bishop developing move (not Ba6 or Bb5), Bg7, another developing move
that does not involve the e or d pawns, then c5, d5. The development of the bishop and the knight is then a matter of compositional
management. It is a bottleneck that can be foreseen from the calculation alone of c5 followed by the d5 push with the bishop and knight
undeveloped, where the non-committed moves of both allow a choice of the cadence of development, based on the occupancy of d7.

The designer has already calculated this, and it is the primitive calculation that shows this plan to be optimal for Black once the
purpose is understood and the points where deviation is forced are known.

**How the certainty arises (a statement).** Much of what chess is at an optimal level is an invisible scoping mechanism, under the
calibration of the optimality of one's opponent. In population-level play at a high level, the rules that govern optimality do not
change with who plays; they converge, though perhaps with different perspectives and understandings of them. At the highest level,
the program will be able to aid in understanding and identifying why a move was made: whether by deliberate calculation of the
primitives, or because the position looks good and is memorable to the player from a branch with good outcomes. Whatever the
mechanism of arriving at optimal play, optimal play is optimal play. The designer's certainty of Black's schematic is not that experts
play it and so it must be optimal. It comes from a population-level observation of a pattern: a common invariant opening structure
for Black, commonly played in speed chess across many players, in which each of the three moves has a direct, explainable impact on
the primitive race conditions, especially relative to the initial start position and to development, and relative to other moves and
their potential developmental cost. The deduction was made from the principles of the program, and it is seen to align with what is
observed in optimal play.

*Correction recorded: an earlier remark by this assistant treated agreement with expert play as the designer's evidence and called it
weak. That misread the statement above. The evidence is the deduction from the primitives, and the observation is what directed it.*

**The comparison against alternatives (a statement).** The alternatives are alternative moves, measured on the primitive catalog. If the
primitive catalog is exhaustive, the three moves can be evaluated on it, relative to the relativistic perturbations of the
alternative moves, to measure explicitly the opportunity cost of each move relative to its alternatives. The designer's sentence
was left trailing, and it is recorded as given.

**The claim made of the deduction.** Even an incomplete program can inform optimal understanding and over-the-board decisions. Using
the primitives already in the catalog and the program built to this point, the deduction identifies and explains why this
schematic is more optimal, in a generalized sense, than many other Black opening plans. The designer offers the deduction as an
example of this, and asks that it be recorded.

**The first milestone (a statement).** The first step is to catalog certified local boundary board state positions. Finding the optimal
forward path to a candidate can be certified afterwards; the classification as global comes later. What is sought first is a path to
a position where the primitives can be enumerated compositionally well enough to certify that it is in equilibrium. A first
certified local boundary board state position would state how it is arrived at from the initial board position, the primitive
calculations and the composition, and the exhaustive coverage certification. The protective belt then discovers and catalogs the
emergent observable facts that arise from the population of boundary board state positions, if any exist. That is to be discovered.

**The progress to be denoted.** The designer asked that progress on the discovery of where to start be denoted, and whether a
document existed for it. This file is that document.

## 3. Progress

| Item | State |
|---|---|
| Black's schematic (g6, Bg7, d6) scoped as a rich vein from the initial position | Scoped by the designer; detailed in survey §4e |
| The deduction of why it is more optimal than many other Black plans | Stated by the designer; recorded here; not compared in the documents against any other plan |
| White's choices as constructor against it | Begun: the worked example narrows White's first move to the central pawns (§2); not complete |
| Bird's Opening | From another plan, per the designer; the other documents still mention it and are not changed |
| A pruned structure on any line | Not run |
| The first certified local boundary board state, with its path from the initial position, primitive calculations, composition and certification | Not found; this is the first milestone |
| Tempo as a scoping criterion on this line | Stated (survey §4e); not applied |

## 4. This file's reading, offered for correction

- **The deduction is a prediction the structure can test.** Its properties (progress on the development race primitives, no negative
  effect on the others, no commitment to a branch) can each be checked move by move in the catalog's terms. The same check on other
  Black plans would give the comparison the designer describes. Neither check is recorded as done.
- **The comparison in the record's terms.** The perturbation of a move is its change in the primitive state (topology §8y, ΔS). The
  comparison is of the perturbation of each of the three moves with that of each alternative move at the same point. Where one
  move's perturbation is at least as good in every component, no weighting is needed. Where components trade against each other, the
  record's own position is that values come from the pruning and the losing conditions of the primitives, not from a score, so the
  ordering would come from there. The catalog has three primitives at present, so a comparison now would carry a stated residual
  for what the missing primitives would add. This is a reading, and it is not yet specified how the comparison is recorded.
- **Over-the-board use is separate from the closure claim.** Using the catalog to explain or choose moves in play does not require
  the closure of survey §2a, and it does not test it either.
- **A fixed Black schematic fixes Black's pool.** What the survey finds against it is a candidate line. The closure still requires the
  other Black replies (survey §4e).
- **Bird's Opening.** The earlier plan (From's Gambit as the first real branch, the fool's mate subtree as the dry run) is described
  in the survey and the pipeline document. The designer says it came from another plan. Whether it is dropped or kept as a dry
  run is not stated, and nothing is changed here.

## 5. Not yet specified

White's candidate constructions against the schematic; the order in which Black's three moves are taken in the survey; the other
Black plans to compare; how the comparison is recorded; and how the schematic's variant with b6 and Bb7 is treated beside the g6
line.

## 6. Provenance

**The designer's:** the observation of the Black schematic in high-level play; the three moves and their uniqueness to the
primitives; that the Black moves are scoped and White's side has not been examined; that the discovery is still in progress; and the
claim made of the deduction. **This file's:** the table, the section order, and the readings in §4.
