# Boundary Board States: Where to Start (the discovery of the first line to survey)

**Status.** A record of the discovery phase, not a method and not a result. It holds the designer's statements about which line the
anchored opening survey (`boundary_board_states_opening_survey.md`) should start with, and what has and has not been examined. It
decides nothing. Nothing in it has been run on a pruned structure.

**Companion documents.** `boundary_board_states_opening_survey.md` (the method; §4e is the statement this file expands, §7 is the
calibration of a line and a depth), `boundary_board_state_topology.md` (§8aj, §8al, §8an, §8ap, §8aq), `boundary_board_states.md` (changelog entries 38, 39,
42, and 62 to 70), `boundary_board_state_primitive_catalog.md` (the three primitives). The designer's statements are marked as such and kept
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

**The calculations in the example are primitive calculations (a statement).** The calculations in the example are literally a primitive
calculation: the race of the chained dependencies primitive (catalog §7), focused on scoping Black's opening schematic for optimality
analysis from the anchored opening starting position. Development is a segmented chained dependency, anchored to the initial position.
It emerges as a primitive calculation and is not framed as a development primitive: it arises through the segmented chaining of
dependencies, with the perturbation of committed moves as the pruning mechanism. White's development is calculated in parity from the
same anchor.

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

**The emphasis of the workflow: a population of certified local boundary board states (a statement).** A local boundary board state does
not require optimal move choices on the forward trace. It requires a path understood to be reachable in a branch, with optimal
understanding from that point forward, and the state certified. The only difference between a local and a global boundary board state
is the exhaustive certification of the optimal move choices from the anchoring initial point that reach the state in equilibrium; a
local state is just in equilibrium. Populating a catalog of certified local boundary board states allows a population of equilibrium
positions to be analyzed, to understand the nature of these equilibria. If properties can be derived that are shared by a population of
local boundary board states and by all boundary board states, then the boundary board states can be identified and understood more
easily, and the modular certification of global versus local can be focused on afterwards. Therefore, when first deriving from the
start, focusing on global constraints for exhaustive coverage is a huge mistake, without first discovering more information from a
population of local boundary board states, which are more easily certified and may themselves lie within the Syzygy scope. If one
boundary board state within the scope is certified, others can be deduced from the standpoint of reachability. The emphasis of the
workflow hinges on the discoverability from certifying a population of local boundary board states: to see a true example of primitive
joint compatibility with the race condition calculations in equilibrium. That is the designer's real goal from this point, the first
certified local boundary board state.

**A second route to the first certified local state: construction, decoupled from the forward trace (a statement).** In addition to
the pruning mechanism recorded in the survey, the designer proposes populating the catalog of certified local boundary board states
by constructing positions from piece combinations, isolating piece combinations and primitives to see what emerges, and analyzing the
structures through logic from a constructed position. Endgames with pawns may trivialize the process enough to populate some of the
catalog and construct a first local boundary board state position. This decouples the route from the initial-position anchor and
from the forward trace. In this setting perturbation changes piece positions without the legal move restriction (topology §8ap).
The designer's view is that it may help populate the catalog of certified local boundary board state
positions, which may be beneficial in itself, and that keeping the pruning mechanism recorded as it is while populating the catalog
through construction can be a worthwhile use of time. Endgame positions may be less complex and more calculable, and many can be
reduced by dead pieces (dead pieces can also be added to expand the catalog, or removed to isolate the true primitive composition of a local state; topology §8an). With fewer than 8 pieces, Syzygy determines win, draw or loss and can be navigated for any of the steps;
what the designer's machine can store (perhaps 5 pieces at most) is to be considered later.

**A third path: opening theory branches that force draws (a statement).** The designer also has in mind analyzing opening theory
branches to identify branches that result in forced draws, if possible and identified as optimal. Interesting results may come from
this alone. Syzygy or opening theory may be searched to identify candidates for local boundary board state positions, or positions
that are close.

**The procedure on any seed (a statement).** Both routes, and the pruning route, use one procedure on any position: backtracing and
reachability under optimality, perturbation, the nature of the primitive race, and joint compatibility. It is used on a decisively
won position to discover local boundary board states, and on a certified one to determine the local boundary board states reachable
from it. It is recorded in `boundary_board_state_topology.md` §8an. The example below is an illustration of it and not the method.

**An illustration of the procedure on a decisive seed (the designer's example).** The position is from a daily puzzle: White: Ka1,
pawns a2, b2, Bd4, pawn f6. Black: pawn b3, Qa4, Nc5, pawn h7, Kh8. White to move. It is a seed, not a candidate: White's race is
decisive, and the point of the example is the perturbation that would be needed to move it toward an equilibrium.

- **White's race.** f7+ and then f8, promoting to a queen, to mate. The Black king has only g7 and g8, since h7 is occupied; f7+ is a
  discovered check from the bishop and leaves the king no viable move. No black move can obstruct, and no black piece prevents f7.
  The only viable reply is the bishop capture, and then f8 follows with no obstruction possible. The bishop's discovered check and
  the pawn position form the chain that enables the promotion. The only ways to obstruct would be a black piece on f7 or f8, or a
  piece occupying coverage threats without modifying the race under optimality. Such a piece would also need to contribute to the
  mating net to be placeable in an equilibrium that draws by a cycle of repeated moves, and not by stalemate or insufficient material.
- **Black's race.** A back-rank mate race, made possible by the white pawn occupancy and the black pawn on b3; the white king needs 3
  moves to escape that threat. As the designer counts it, White's forced mate and Black's back-rank mate take 2 moves each, so the
  race is equal at 2.
- **Dead pieces.** The black queen and knight are dead to White's mating net, except that the queen can capture the bishop. A knight,
  a second bishop or a rook there would not matter.
- **The one-move mate is set aside.** The mate in one in the corner where White's king stands is not considered, because it is not
  jointly compatible: a race that close to the king cannot obstruct the race on the other side of the board. White's two-move race
  beats it because of the forced nature of the race that delivers it. For the position to be viable it has to be navigated toward
  obstructions in the middle of the board. The back-rank mate in 2 with an obstruction meeting the bishop, or an obstruction on f7 or
  f8, is the only viable route.
- **The perturbation and its directionality.** The directionality of the perturbation emerges from this. The race is to be composed
  so that Black's back-rank mate is interdependent with obstructing White's race ("cannot be delivered due to dependency to black f7
  pawn advancement", as the designer wrote it), so that maintaining the obstruction is the only winning condition and committing to the
  mating race loses the race that is being obstructed. Making the moves to mate longer by perturbing the position is how the other
  route toward obstruction is reached. The position is not jointly compatible as it stands, and it is used as an example to run the
  four steps.

**The progress to be denoted.** The designer asked that progress on the discovery of where to start be denoted, and whether a
document existed for it. This file is that document.

## 2a. The Syzygy mining process (the designer's statements, and the procedure as worked through)

**Status of this section.** The process below is a specification, recorded so that the work can be resumed. Nothing from any earlier
test run is carried over: no data, results, scripts or database. The designer will run the process and verify the data on returning
to it. The boundary board state understanding is first to be confirmed as covering a majority of the entire space. The designer also
intends to explore a further angle that may change the process, and this section records it as it stood.

**Mining the tablebases to begin the catalog (a statement).** The Syzygy tablebase solution set can be used to enumerate the draws of a
material, giving a database of proven draws under optimal conditions. The catalog of boundary board state positions and its
reachable families can be created from those proven results. The designer will mine all the Syzygy data the machine can handle,
then automate off that basis, from the catalog, with the pruning mechanism. Eventually the size of the Syzygy set does not matter,
because reachability under the methodology of the programme, backtracing and exhaustive enumeration from the cached database file,
carries past it.

**How a draw ends (a statement).** The draws are categorized by whether the paths draw from stalemate, from insufficient material, or
from cyclical draws.

**Branches under optimal constraints (a statement).** A branch is enumerated by taking the tablebase's win, draw or loss value and
distance and mapping it to the actual branch with its legal move pool, to determine which branches align with the count. Branches
are created under optimal constraints: the moves chosen are the optimal ones.

**Reduction to insufficient material (a statement).** Where the tablebase verifies that the best move is a reduction to insufficient
material, that is a provable equilibrium of reduction to insufficient material and not a choice of stalemate, since optimality is
what matters here.

**The 50-move rule (a statement).** It is ignored for now. It matters only for a reachable family once the anchor is reached, where
it gives a cutoff from the number of moves.

**The primitives afterwards (a statement).** The primitive compositional analysis is not part of the enumeration. It is done post hoc,
from the population of certified local boundary board state positions obtained by exhaustive coverage of Syzygy, and by the designer
when a draw position is verified.

**The criterion (a statement).** A position in the choice region that has a move that can lose is a boundary board state position, and
the initial game start position is one. Upstream positions that under optimal conditions lead to the boundary share its equilibrium.
A distinction drawn by how many drawing moves a position has is a meaningless classification criterion. A family reachable from the
initial start position under the optimal backtrace constraint has every member global (survey §5a).

**The two spaces (a statement).** There is a choice space and a forced space that together encapsulate the entire space of chess, with
the areas of sub-optimal play pruned. The boundary between them is literal. The pruned structure that emerges from the process is the
full solution of chess.

**The order (a statement).** Mine Syzygy. The catalog of boundary board states and the reachable families come out of that mining.
Then the process is automated and the primitives are matched.

**The procedure as worked through (this file's organization, offered for correction).**

1. **Choose a material and find its leaves.** Every capture, and every promotion for a material with pawns, leads to another material,
   whose table is then needed as a leaf. For two kings with a rook each, the leaves are a king and rook against a king, for either
   colour, and king against king. A missing table makes the probes fail.
2. **Enumerate every legal position**, with either side to move: distinct squares, kings not adjacent, and the side not to move not in
   check.
3. **Fold by symmetry.** Keep one canonical representative per class and its orbit size, so that counts over all positions can be
   recovered. With no pawns the symmetries are the eight board symmetries and the swap of colours with the side to move swapped. With
   pawns the direction of play limits them to the left-right mirror.
4. **Take the tablebase value** of each representative by its sign: win, draw or loss for the side to move. Cursed and blessed values
   count as wins and losses, because the 50-move rule is ignored. Check coverage: every legal position has a value.
5. **Check the table against an independent move generator.** Derive each position's value from the values of its children and compare
   (checkmate and stalemate where there are no legal moves). An inconsistency means an error in the move
   generator, the symmetry fold or the reading of the table. A sample is also compared with a second tablebase reader.
6. **The drawn population** is the set of positions with value zero for the side to move.
7. **Classify the moves of each drawn position** by the value of the resulting position for the mover: losing, or drawing (a winning
   move cannot occur). Record counts and not a verdict: legal, losing and drawing moves; and of the drawing moves, those that stay in
   the material, those that leave it by a capture into a terminal draw by stalemate or insufficient material, and those that reach a
   stalemate position. Storing counts lets any criterion be applied later by query.
8. **Build the graph of drawing moves** among the drawn classes. Compute the strongly connected components (positions that can reach one
   another, which is where the cycles are) and the weakly connected components (the reachability families). Propagate flags: a cycle, a
   stalemate or insufficient material is reachable by drawing moves. Label the realization as stalemate now, terminal only, cycle only
   or mixed. A cycle in the symmetry quotient is a real cycle: a path from a position to its own symmetric image repeats after the
   order of the symmetry.
9. **Apply the criterion by query.** Drawn with a move that loses is a boundary board state position. Drawn with no losing move is
   interior. Not drawn is forced space.
10. **Verify the components** with a second implementation.
11. **Store the result in a database.** For the post hoc composition analysis, also keep the edges of the graph, the losing moves with the
    positions they lead to, and the refutation distance of each loss (the distance to zeroing from the distance tables, and true mate
    distance only where a solver supplies it). These three are not built yet.
12. **Link to the anchor.** A family reachable from the initial position under the optimal backtrace constraint makes every member global
    (survey §5a). The catalog entries are certified entries, and act as leaves for any structure upstream of them (survey §4f). For
    materials beyond what the machine can store, the segmented method builds upward from the cached entries, and the automated pruning is
    checked against the tables wherever a table exists.
13. **Order of materials.** By piece count, up to what the machine can hold (the designer estimates perhaps five pieces).

**What is carried over.** Nothing from any earlier test run. The designer verifies the data on returning to the process. Step 11's
additional tables and step 12's backtrace are not specified further than above.

## 2b. The staged workflow from an anchor (the designer's statements)

**Status of this section.** The designer describes this as a staged process of the pruning mechanism and as the sum of the programme so
far, and expects it to be more automatable. Some rules are not fully resolved, as stated below. Nothing in it has been run.

**From the initial position to a segment anchor (a statement).** The initial game start position is assumed to be in the equilibrium
choice phase, and its legal move pool is vast, so branching on all of it is combinatorially explosive. The structure is treated as a
sum. One move is chosen from the pool and worked down its branch, and the whole is composed from the branch segments. A preset opening
path takes the work from the initial anchor to a pruned structure anchor, where the pruned structure begins. Segments are certified up to
a specific point, composing them may be additive once they are certified, and the rules of composition are to be understood later,
because they depend on what certifying the pruned structure involves.

**Discovery at a segment anchor (a statement).** From a pruned structure anchor the work is discovery, not building from the catalog.
The goal is to discover a local boundary board state, in order to create a reachable family segmented progression certificate. Since
the target is local, optimality of the path is not the concern at this stage. The object is the full landscape from the pruned
structure anchor to the terminal local boundary board state positions.

**Enumeration (a statement).** All legal moves are covered exhaustively from the anchor, keeping only the pruned structure that
maintains equilibrium, and keeping track of the coverage, so that the network of legal moves is built by brute force under a rule set
that stops it trailing off into repeated moves. Because the structure is pruned for equilibrium and not for advantage, the designer
expects the result to be much smaller than the decisive branches.

**The pruning rules (a statement).** Decisive wins and losses are pruned. A position already in the pruned structure, arrived at again,
is cached as work already done. The other pruning rules concern the optimality of reachability and opportunity cost against the other
potential branches, using recurring board positions and the measurement of relativistic branches against the structure already
enumerated, so the protocol is compositionally and relativistically interdependent. A self-referential procedure is to rule out a
sub-optimal route: one that arrives at a position that could have been arrived at before by optimal pathing, which is evident from the
pruned structure. The designer says the rule set for repeated positions and the self-referential procedure is not fully resolved, and
that it may work itself out from running the process from the start.

**Reversible and irreversible moves (a statement).** There are moves that can be reversed and moves that cannot. An irreversible move is
a true moment of progression that reshapes everything. A reversible move can be reversed, and sometimes it is optimal to reverse a move
or to return eventually to a position, so pruning such moves is trivially against the programme. The dancing is mitigated by referencing
the keys, to know whether a move is within a cycle in a segmented boundary that is already resolved. Sometimes it is good to move within
a segmented boundary for zugzwang purposes, if the opponent is forced to. A segmented boundary that is repeated to draw conditions is an
important indication of a boundary board state position.

**Mapping is not certification (a statement).** The above maps out the landscape. Local boundary board states are concerned with
optimality, so the potential draw positions found in the structure must be certified through optimality, by referencing the other
branches pruned for decisive reasons. At each pruned end the composition of the win and loss potential of each player is understood.
This determines whether a draw position is reachable under optimal conditions: if there exists a branch that a player can choose that
deviates into a forced win, that branch can be chosen instead of the draw under optimal conditions.

**Which regions to enumerate (a statement).** This informs which regions are enumerated at all. If a point has been reached where both
sides could already have forced a win, the work is already in the forced phase region, and the boundary between choice and force has
been crossed. That boundary then resides behind the pruned structure's anchor point, or within the object being enumerated, or the
initial position can be forced, which the designer considers most unlikely.

**Certifying a draw by the path to the local boundary board state (a statement).** The optimality of pruning potential draws is analysed
through the path between the proposed draw point and the ply of the local boundary board state position, under the optimal scoping of
whatever equilibrium can be exhaustively verified, which needs to be within the Syzygy scope. This is where the programme aligns with
exhaustive coverage: the scope of the forward trace needs to overlap with the backward trace to enumerate the local boundary board
state catalog. This is how the local catalog is populated with Syzygy, and the purpose is to identify and certify the space between
the forward and the back trace through the construction of the pruning structure.

**The cost, and the worst case (a statement).** Some underlying pattern must be found in the local boundary board state catalog that
can be enumerated, or the programme reaches the same limitations as Syzygy. Even then, the distance to a weak solution can be scoped,
which the designer regards as massive progress. Perhaps the local boundary board states create local boundaries that are
compositionally defined with the same structure and principles as the global ones, and that form their own informative decisive force
and choice boundary, or display holes for decisive force sequences that make optimality obvious. In the worst case, too little can be
enumerated within the Syzygy limitation to construct a population from which to develop an understanding of these structures at all.

**Modularity and community effort (a statement).** Linking to initial anchored pruning positions and the game start position is a good
way to approach the work in segments, and makes it approachable as a community effort if it can be shown to have a viable path at all.

**On novelty (a statement, not checked by this assistant).** The designer's belief is that this is concrete and provable, and that from
the designer's understanding and an earlier literature check it does not exist in the literature.

**Where the regions are, and closure outside the anchored region (a statement).** A move that wins is decisive and lies outside the choice
region. Choice phases and regions can exist that are disjoint from the initial start position, and they retain equilibriums and spaces
that contain reachable families of local boundary board state positions. The process is mining for them, and the designer expects to
find them eventually. Whether enough of them will be found to populate the catalog well enough to inform the composition of the global
boundary board states is still uncertain. The aim is to populate as many of these sparse regions as possible, and, with luck, to find that
the global boundary between choice and force already overlaps the back-trace region the Syzygy tables provide. The closure property for
regions can also exist outside the global region anchored at the initial position. Not every local boundary board state region has it:
some are unclosed, with decisive routing through the reachable family. A region that closes the branch graph network of a reachable
family is to be remarked upon as significant, because it implies a lack of forcing sequences and a localized choice phase or region under
optimality. It is a significant finding if ever found. A localized or global closed region will share the properties of the global
choice phase or region if the programme's presupposition holds.

**Mining without an anchor, and the choice of anchor (a statement).** The Syzygy landscape can be enumerated into reachable families, with
the family constructed backwards and met forwards. The boundary board states can be mined locally and the pruned structure basis built with
no anchoring at the initial position, linking back to the initial position afterwards by surveying the Syzygy database. The work
reaches the limits of Syzygy when it reaches them, and what can be enumerated with Syzygy is enumerated. The forward anchor can be the
initial position or wherever is chosen. A decisive position can have a choice region below it, reached if the player who is winning
deviates from the forced region and plays sub-optimally, or through ignorance of optimality, which is what the programme maps out. A
localized choice region can therefore lie past the global choice region. There cannot be disjoint global regions that share the initial
start position, by the geometric definition of global.

**The structure of the whole (a statement).** The global choice region leads into the global forced region. Localized choice regions lie
within the forced region, disjoint from the global choice region, and are therefore unreachable from the initial position under optimal
play. They are reached only by sub-optimal play or through ignorance of optimality. This resolves how the survey's list of global
boundary states relates to the two spaces.

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
| A second route: construction of endgame local boundary board states, decoupled from the forward trace | Stated by the designer (§2); not run |
| A third path: opening theory branches that force draws | Stated by the designer (§2); not run |
| The procedure on any seed (decisive position or certified state) | Stated by the designer; recorded in topology §8an; illustrated in §2; not run |
| Syzygy mining of a material into a database of proven draws (§2a) | Process stated by the designer and recorded as a specification; nothing from any earlier test run is carried over; to be run and verified by the designer |
| The boundary board state understanding confirmed as covering a majority of the whole space | Pending, per the designer, before any results are integrated |
| The link from a mined family to the initial position | Stated (survey §5a); the backtrace procedure is not specified beyond the statements |
| Closed local regions and sparse disjoint regions (§2b) | Stated by the designer; none found; whether enough regions exist to inform the global composition is uncertain |
| The staged workflow from a segment anchor (§2b) | Stated by the designer; the self-referential prune and the repeated-position rules are not fully resolved; not run |

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
- **The illustration.** The position has ten pieces, outside Syzygy, so the tablebase cannot check it.
- **The anchored and the unanchored route.** The first-milestone and workflow-emphasis statements above describe the anchored route,
  where a state is reached by a path from the initial position. The construction route has no backtrace constraint at discovery; the
  family composed from a discovered state is what links it to the initial position or concludes that it cannot (survey §5a, topology
  §8aq). The two are read as covering different routes into the same catalog. Offered for correction.
- **The criterion is applied to the side to move.** The mining applies the designer's criterion to the moves of the side to move at each
  position. The opponent's deviations are seen at the positions reached by a drawing move. The survey's hard core asks for a losing
  deviation for each side (survey §2, §4a). Offered for correction.
- **Forced reduction.** The reduction to insufficient material is read as the positions where every drawing move is a capture that ends
  in insufficient material, or where no cycle is reachable and every drawing line ends there. A table of values does not rank moves of
  equal value, so "the best move" is read as the unique drawing move. Offered for correction.
- **The two spaces and the weak solution.** The statement is read as the weak solution of survey §2a and §4g: the region of optimal play
  with every deviation out of it certified as a loss. The boundary is read as the rim of the choice space, which includes the initial
  position. The survey's sentence about global states forming a list, with the space between the initial position and that list, is
  unchanged. Offered for correction.
- **Keys and segments.** Nodes are read as keyed by the position and the side to move, so the structure is a finite graph and a repeat is a
  cache hit. A cycle can only lie within a fixed material, pawn structure and castling rights, so the irreversible moves (captures,
  pawn moves, loss of castling rights) are read as the boundaries between segments. Offered for correction.
- **The region is a fixpoint.** A node is in the region only if no move wins for the side to move, and every move is either a certified
  loss or leads to a node in the region. A branch found to be a forced win for the mover therefore removes the node from the region, and
  that removal passes up to the parents whose moves led into it. So the region is final only when every move of every node has been
  classified, and not when the forward enumeration stops. This is the designer's "if a player can choose a branch that deviates into a
  forced win" put as a rule. The survey's certification condition (§4a) says each move "stays inside the region as a non-losing move";
  it is read here as "leads to a node in the region", since a move into a position that wins for the mover is non-losing and is not in
  the region. The survey is not changed. Offered for correction.
- **What the status of an anchor decides.** The earlier reading here, that a decisive anchor "has no region below it", was wrong and is
  withdrawn (the designer's correction). A decisive position can have a choice region below it, reached when the side that is winning
  plays sub-optimally or through ignorance of optimality, and mapping such regions is what the process is for. The status of an anchor
  therefore does not decide whether a region exists below it. It decides what claim the result carries: a region at the initial
  position carries the global closure claim, and a region below any other anchor, decisive or not, is a local finding. Offered for
  correction.
- **Where the traces meet.** The families inside the tables are exact. The forward structure meets them only where it descends, by
  captures and promotions, to a material the tables hold, so the meeting is at the table boundary. A predecessor of a family member from
  a larger material is outside the tables, and that step is the pruning mechanism of the larger material. The forward side, between the
  full material and the tables, is the large part. Offered for correction.
- **The tablebase end of the overlap.** The part of the forward structure that the tables reach is the part that descends, by captures
  or promotions, to the material the tables hold. A boundary state that keeps more pieces than the tables hold closes by a cycle at a
  fixed material and is not reached by the tables. The first population is therefore of low material, and whether its pattern carries
  to higher material is the cost the designer names. Offered for correction.
- **Precedent.** Forward search meeting backward databases is the structure of the checkers precedent already cited in topology §7. This
  assistant has not checked the literature on the programme's own claims (primitives, the boundary as the rim of the choice space,
  reachability families, a catalog of boundary states), and the novelty statement in §2b is recorded as the designer's.
- **A closed local region and the closure claim.** A reachable family is closed when every move of every member either leads to a member or
  is a certified loss for the side that plays it, and no move wins for the mover. This is the region definition applied to a family
  that need not lead to the initial position. The global closure claim (the programme's closure and its failure condition, a certified
  route out) is about the region anchored at the initial position. A closed local region is a finding of its own, with no failure
  condition attached, and an unclosed one is a family with a certified decisive route out. The two are kept apart in the record.
  Offered for correction.
- **What the 50-move rule can change, and what it cannot.** It can only turn a win into a draw, so a branch cached as decisive with the rule
  ignored may be a draw on a route with a high move counter. The draw claim of a certified region is not affected, since everything in it
  is a draw under the stricter reading. What can be affected is the extent of the region: a move labelled losing may be a drawing move on
  such a route, so the region found may be smaller than the real one, and those routes remain to be enumerated once the anchor is
  reached. The counter resets on a capture or a pawn move, which are read as the segment boundaries, so a verdict at a segment entry does
  not depend on the route, and inside a segment a cached decisive verdict is exposed only where the distance to zeroing plus the counter
  on arrival comes close to 100 plies. The certification of a region is valid up to where it is valid. Offered for correction.
- **The list of global states and the localized regions.** With the statement of the whole, the survey's list of global states is read as
  the rim of the global choice region, and the localized choice regions are the disjoint pockets in the forced space. A reachable family
  that links to the initial position under optimal backtracing is therefore a family on the rim of the global region, and a family that
  cannot link lies in a pocket. Offered for correction.
- **Closure is a test, applied family by family.** It follows from the fixpoint reading above that a family is tested for closure by
  classifying every move of its members, and an unclosed family is identified by finding one move into decisive territory. Offered for
  correction.

## 5. Not yet specified

White's candidate constructions against the schematic; the order in which Black's three moves are taken in the survey; the other
Black plans to compare; how the comparison is recorded; and how the schematic's variant with b6 and Bb7 is treated beside the g6
line.

The backtrace that links a mined family to the initial position, beyond the designer's statements; the order of materials; a tablebase
reader suited to the size of the material; which of the edges, losing moves and refutation distances to build first; and whether interior
positions are entered in the catalog.

The self-referential prune and the rules for repeated positions in §2b; how a segment anchor is chosen and checked; the preset path from the
initial position to the first segment anchor; and how certified segments are composed.

How the closed local regions found are to be recorded and flagged in the catalog, and how an unclosed family is recorded with its decisive
route; and the criterion for stopping the search for sparse regions.

## 6. Provenance

**The designer's:** the observation of the Black schematic in high-level play; the three moves and their uniqueness to the
primitives; that the Black moves are scoped and White's side has not been examined; that the discovery is still in progress; and the
claim made of the deduction; the second route, the third path, the procedure on any seed and the illustration (§2); the statements of §2a and §2b. **This file's:** the table, the section order, the organization of the procedure in §2a, and the readings in §4.
