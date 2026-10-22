# Boundary Board States: the Anchored Opening Survey

**Status.** A method document, not a program and not a tool. It specifies how the primitive catalog is used to look for a
boundary-state candidate along a chosen opening line, anchored at the initial game position. Nothing in it has been run. The
hard core (the boundary board state) and the protective belt (the three primitives) are unchanged by it. Evidence tags as in the
catalog: **SPEC** (defined, not tested), **INSTANCES** (checked on specific positions), **CERTIFIED** (backed by an exhaustive or
complete case-list certificate).

Companion documents: `boundary_board_states.md` (the concept and methods; §1 the definition, §3 the caching architecture, §14
items 4 and 7 to 11), `boundary_board_state_topology.md` (the record; §3 status is local and unconditional, §5 start-reachability,
§8h footprints, §8p the guard rail on absent paths, §8r forward anchoring, §8x the proven entry, §8ac, §8ad, §8ae, §8af and §8ag, and §8ah to §8ar, the designer's
statements behind this file), `boundary_board_state_primitive_catalog.md` (the three primitives, the admission rule),
`boundary_board_states_case_study_otb.md` (the game), `boundary_board_states_start_selection.md` (where to start: the
discovery of the first line to survey), `boundary_board_states_syzygy_mining.md` (the process of mining the Syzygy tablebases
for the certification catalog).

**How this file is organized.** What is the designer's statement is marked as such and kept apart from this file's organizational
choices (§11), in the manner of the catalog's §10. Where this file reads one of the designer's statements against an existing
definition, it says so and offers the reading for correction.

---

## 1. What the survey is for

The earlier methods approach the boundary state from the far end: retrograde from mating positions, population analysis of
confirmed states, construction from primitives. They leave out the anchoring point of every game, which is the initial position,
and the path from it. The survey starts there. It takes a specific opening line, chosen for a reason (it is well understood in
theory, and its branching is known well enough to work with), and uses the primitive catalog and the perturbation process along
it to look for a **local boundary board state candidate**: a position where the primitives show an equilibrium.

It is a scouting exercise, not a build. It asks for no script that deploys the primitives on any opening. Its product is
candidates with their compositions, the boundary between choice and force along the line as far as it can be located, and cache
entries (§6) recording what work has been done.

It does not aim to establish that a position is global, to find the optimal opening, or to decide whether a global boundary state
has already been passed. Those are consequences that may fall out of it (§5, §6). They are not its target.

## 2. Definitions

**The designer's, as stated.**

- **Choice and force.** Between the initial position and any boundary state there is a region of choice, where each side has
  non-losing alternatives, and a region of force, where the moves are decided by the race conditions. A boundary state is where the
  region of choice closes: the first place along the way where the choice-versus-force structure is reached.
- **Forcing sequence.** A sequence that gives a decisive advantage to the side that chooses to force it. A side may force a win, a
  draw or a loss. Forcing a loss is not informative. Forcing a win marks a winning race condition. Forcing a draw marks a possible
  boundary state. A position where each side's every alternative either loses or is answered by a forced sequence that keeps the
  draw, so that the draw is forced at each step to prevent a forced loss, shows boundary-state behavior.
- **Local boundary board state.** A position where the equilibrium is verified: the catalog's primitives show it (the composition
  of `boundary_board_state_primitive_catalog.md` §5), and exhaustive verification shows that the primitives did their job. A local
  boundary state is not required to say how it was reached, and it is not required to be the first deviation point between choice
  and force. It informs us of the viability of the primitives, and it is the point on which exhaustive verification is performed. The exhaustive pruning of §4a, closed and fully certified, is that
  verification and is the certificate.
- **Global boundary board state.** A boundary board state that is a member of a reachable family that is reachable from the initial game
  position under the optimal backtrace constraint (§5a). Every member of such a family is global, and the initial game position is itself
  a boundary board state position. Whether a state is global therefore depends on how it is reached, by optimal play from the initial
  position. A state not shown to be in such a family is local.
- **Global is sound against every deviation, and there are many.** A global boundary state must hold holistically: it is reached by
  non-losing play, and it holds against all deviation replies. A deviation from the path may lead to a different global boundary
  state, a detour to another position that closes the phase between choice and force. Each global
  state closes the phase for the routes that lead to it, so the global states form a list, and the space between the initial
  position and that list is the enclosed phasing of choice and force. That a position can be avoided is therefore not a ground
  for classifying it as non-global. A verified state not shown to be in a family reachable from the initial position under the optimal
  backtrace constraint (§5a) is local.
- **The whole, as the designer states it later.** The global choice region leads into the global forced region. Localized choice regions
  lie within the forced region, disjoint from the global choice region, and so are unreachable from the initial position under optimal
  play (`boundary_board_states_start_selection.md` §2b). The list of global states above is read as the rim of the global choice region.
  This is the designer's resolution of how the list relates to the two spaces.
- **The target of this survey.** A local boundary state candidate, anchored at the initial position through a chosen line with a
  viable path of arrival. The path may be optimal or sub-optimal.

**What stays as recorded.** A boundary state's own status is local and unconditional (topology §3): it is evaluated from the
position's legal moves and the values of its children, and nothing about its ancestors can change it. The concept is a property of
the position (`boundary_board_states.md` §1). This file's definitions concern the *classification* of a verified state by how it is
reached, so they leave §3 and §1 untouched.

**This file's reading, offered for the designer to correct.** The record already defines a start-reachable boundary state
(topology §5; Tier 3 of `boundary_board_states.md` §14 item 1) as one connected to the initial position by a chain in which both
sides always had a non-losing option and neither had already played a losing move. The reading here is that *global* is
start-reachable in that sense and also first, across all such chains, so that the choice region is the region reachable by
non-losing play and the boundary is where it closes. A local candidate anchored by a sub-optimal path is then a local state with
no start-reachability claim. `boundary_board_states.md` §14 item 4 noted that "local" and "global" carried four meanings and needed
a firm definition. The definitions above are offered as that definition.

## 2a. The closure claim

**The designer's statement.** The closure of the global boundary board states must be ensured, otherwise the programme fails. A
viable deviation, while the position is still in the choice phase, always reaches a boundary board state. A route cannot go from the
choice phase to the forced phase without hitting one, and if it does not reach the same boundary state it reaches a different global
one, by definition. If this were not true, a win could be forced from the initial position: an opening in the closure would be the
forcing opening that decides chess from the moment it is taken, with no blunder possible. No such line is known to exist in all of
chess theory, so the closure is predicated on the fact that one does not exist, and the enumeration of the closure would be the proof.

The closure is expected to be achieved through the cyclical structure, or through zugzwang that forces eventual losing moves, which
classifies the branch upstream once it is discovered. With the ruling of `boundary_board_states.md` §14 item 10, a forced terminal
draw with no choice to deviate is also a realization.

**The existence of the boundary (the designer's statement).** A decisive position is pruned. If the pruning goes all the way up to the
initial position, chess is solved by force from the initial position. If it does not, there is an explicit boundary, either the
initial position or positions downstream of it. This is trivial logic and does not wait on an enumeration: given that the initial
position is not won by force, the boundary exists and cannot not exist. Whether it is the initial position or a closure more distant
from it is yet to be calculated. The boundary board state is the equilibrium of the primitives under joint compatibility: the
ability to force a draw down a line that would otherwise be decisively losing, so that the position would be a decisive loss if not
for the ability to force the draw.

**What can happen, as the designer states it.**

| Outcome | Reading |
|---|---|
| The closure holds on every route and is enumerated | the forward and backward trace is complete, and at least a weak solution is viable |
| A route escapes the closure, with a certificate | a forced opening line from the initial position that decides the game; a discovery that reshapes chess |
| A prediction of the belt fails, or a primitive is missing | not a failure of the programme: the catalog is incomplete, and the reduction rule develops it |
| A primitive is admitted that does not fit the architecture of the solution | a methodological failure; the admission rule (catalog §1) guards against it |

The hard core therefore has one failure condition, which is also a result: a certified route out of the closure.

**Assessment — separate from the record above.**

- **The relation to the record.** That a transition from choice to force exists on every route to a decisive position follows from a
  drawn start (topology §5), and the collective version is topology §3. The designer states that the existence of the boundary
  follows with certainty from this, and the enumeration does not test it. What the enumeration does is locate the boundary and
  certify it, and it would show a certified forcing line if the premise below failed. The closure claim states that the boundary is
  a boundary state in the sense of `boundary_board_states.md` §14 items 7, 8 and 10, which the designer holds to follow from the
  definition. No argument for or against that is recorded here beyond the designer's statement.
- **The premise.** That the initial position is a draw is not proven. The closure is stated conditional on it, and enumerating the
  closure would prove it.

## 3. The survey's objects

| Object | Meaning |
|---|---|
| Anchor | the initial game position |
| Line | a sequence of moves from the anchor, chosen from opening theory, optimal or not |
| Node | a position on or reached from the line, with the side to move |
| Primitive status | for each side and each engaged primitive, COMPLETED, ENGAGED (with held terms) or DORMANT (catalog §4.4, §7.6, §8.6) |
| Decisive node | a node where a side has a forced win: a net completed or completable against every reply |
| Choice node | a node where both sides have non-losing alternatives, recorded with its width |
| Forced-draw node | a node where the draw is held by forced replies against a forced loss, with the equilibrium realized by a repetition cycle sustained by both sides or by a forced terminal draw with no choice to deviate, both sides live (`boundary_board_states.md` §14 items 7 and 10) |
| Candidate | a forced-draw node whose primitive accounting is an equilibrium (catalog §5), entered in the candidate catalog |
| Web | the set of candidates whose verification depends on one another (§7) |
| Reachability family | boundary board states that can reach one another: if a boundary board state can reach another one, it is in the same reachable family (§5a). The classification of a family when it is anchored to the initial position is in §5a. Recorded in the boundary board state certification catalog (§4j) |

## 4. The procedure

1. **Choose the line and the depth.** The line comes from the theory tree. The depth is chosen so that the primitives scope the
   region properly (§7). This is the calibration choice.
2. **Walk the line.** At each node, compute the status of the primitives for both sides, with the colors swapped (catalog §6). A
   forcing sequence appears as a net or chain that the status shows as completable.
3. **Perturb up the tree.** When a forcing sequence is found, perturb the earlier nodes by Method 7 (`boundary_board_states.md`
   §2): find where the race condition was acted on, and how the choice of a piece or move at an earlier node made the forced
   sequence available. This locates where choice becomes force along the line.
4. **Prune.** Adversarial play always takes a forced win when it is available. A move that gives the opponent a completable net is
   a losing deviation, and the lines that contain it are pruned. What remains is the set of non-losing alternatives of each side,
   which is the choice region, and its narrowing toward forced replies is the boundary between choice and force.
5. **Separate the forced outcomes.** A forced win ends the line as decisive and is recorded (§6, C1). A forced loss is not
   informative. A forced draw is a possible boundary state: check the cyclic definition and that both sides are live, then compose
   the primitives (catalog §5).
6. **Enter and verify.** A candidate is entered in the candidate catalog with its composition. Its verification is the pruned
   structure itself, when it is a closed system with a certificate on every refutation (§4a): that structure is the complete case
   list of topology §8d, scoped by the footprint of §8h, and no further step is needed to certify it. Exhaustive search is
   insurance on the first entries. A verified candidate is a local boundary state. Whether it is global is a separate question (§5).
7. **Cache** what was learned, whichever way the search ends (§6).

### 4a. Pruning as information (the designer's statement)

Pruning is not only the removal of data. A forcing sequence that wins is a race condition that is won, and it is identifiable as a
primitive. Perturbing upstream with that primitive in mind shows whether it is baked into the structure or emerged from a chained
dependency or from the order of the move choices that made the race work. Tracing it upstream shows whether an equilibrium can be
reached at all. Identifying the losing sequences and their primitives is therefore informative about the equilibrium of the
potential boundary state, which requires the primitives to be jointly compatible.

The pruned structure thus carries more than its layout. For each edge it carries why the edge is there, which is the race-condition
primitive that refutes the alternative, and which moves activate a race condition and what the calculation on the primitives
becomes under the perturbation. The balance of a position's equilibrium comes out of the pruning process.

If the composition upstream shows that no equilibrium can be established, then no boundary state can exist in that branch. If one is
found, with the pruning excluding branches on the basis of the primitives and the structure standing as an equilibrium of race
conditions, that is a boundary-state candidate anchored from the initial position, for example from 1.f4 and a chosen reply.

**This file's organization, offered for correction.** Each pruned edge is recorded with the primitive that refutes it, the
certificate type, and the upstream classification (structural, or emerged from a chain or an ordering of moves). "No equilibrium can
be established" is a deduction from the joint incompatibility of the primitives (topology §8x), and a failure to compose with the
current catalog is first read as a residual interface, which may be a missing primitive and not an impossibility (topology §8p).

**Joint compatibility and the end of the tree (the designer's statement).** The equilibrium measurement in the primitives comes from the
primitives that lose in the other branches, and that is what informs the pruning. Joint compatibility emerges with the equilibrium
measurement from the perturbation. When the pruning is complete, the furthest downstream positions of the tree are the boundary-state
candidates, with nothing beyond them and the equilibrium cycling or resulting in a forced draw. If the tree does not end there, a
branch from the position being evaluated leads directly to a forced favorable outcome.

**This file's reading, offered for correction.** A conflict between invariants over a shared resource would show in the tree as a move
that is neither inside the region nor a certified loss, so the closed structure tests joint compatibility, and the compatibility
items of catalog §5 (4 and 5) are the account of why the structure closes. A tree that does not close has either an open edge, an
unresolved move that is work remaining or a residual for the reduction rule, or a found branch to a forced favorable outcome, which
refutes the candidate. The two are recorded separately.

**Certification condition (this file's reading of the designer's statement, offered for correction).** The pruning certifies an
equilibrium when the pruned structure is a closed system. At every position of the equilibrium region, which is the cycle together
with each non-losing alternative that stays drawn, every legal move of both sides is accounted for: it either leads to a
position in the region, or it is a certified forced loss for the side that plays it. Then no side can improve by leaving the
region, and the value of the region is the draw. The complete deviation map is the proof, and the primitives supply the reasons for
the refutations. It needs: each forced loss certified (a verified forced mate, a tablebase result, or a complete case list, and not
an engine assertion); closure under the non-losing moves; a losing deviation for each side, with the realization shown (the cycle,
or a forced terminal draw with no choice to deviate); and, for a global classification, the separate proof of §5.

### 4b. Equilibrium maintained until the move pool narrows (the designer's statement)

The initial position resides in the choice phase. There the primitive race conditions are in an equilibrium from which a forced draw
can emerge, with no forced sequence to a decisive outcome. The assumption that it is a draw rests on two things: no forced win has
been found, and no primitive race has yet committed moves to produce an opening, so that, with the symmetry of the starting
position, any forced attempt can be mitigated by the moves available. The primitives available from the initial position have not
committed enough to narrow the branches. Imbalance emerges from decisive moves that create the opportunity for a race to be won, so a
primitive is reached only through the moves that enable it. The fool's-mate primitive in Bird's Opening is the example: it depends on
White making the enabling moves, the g-pawn above all, and Black cannot make moves that force it open, so under optimal play it is
not reached.

Optimal play is then a matter of exhausting the moves that mitigate these primitives, until the move pool and the lack of other
opportunities force an imbalance, in a mechanic similar to zugzwang. Finding the global boundary states is maintaining the primitive
equilibrium until the move pool forces the imbalance.

**The move pool.** The move pool is enumerated in the pruned structure. Its values are not given by heuristics: they come from the
pruning, from the losing conditions of the primitives in the race conditions. A pool with a winning move is decisive for the mover. A
pool of only losing moves is a zugzwang loss, decisive against the mover. A pool with drawing moves is drawn. The choice region is a
drawn pool with many drawing moves, and the boundary is where the drawing moves hold the draw against deviations that decisively lose.

The method follows. From the opening position and its theory, use the same perturbation and pruning to take the longest line along
which the equilibrium is maintained, as the pruned structure shows. Further down the line, isolate when any move from the move pool
enters decisive territory. Where the position can be forced to a draw, it is a boundary state if the moves in the pool that deviate
from the forced draws decisively lose, for both sides. That is visible in the pruned structure. The boundary condition on the pool is
that it offers only decisively losing options or forced draws. If a decisive forcing line exists in a side's favor, the opponent must
be forced to a move that avoids that branch entirely. This interdependent composition narrows the pruned structure toward what is
sought. How large the pruned structure is remains the question.

**Recorded with it.** The drawn start is a premise (§2a), and the enumeration is what would prove it. Whether a given net stays
uncompleted because no carrier compels the move that enables it is established for each position (catalog §8.7, D-CR, with the
carriers K0 to K4 of catalog §4).

### 4c. The process as a whole (the designer's statement)

The approach is to build a pruned structure tree from the opening position onward, focused on a specific opening. The primitives are
understood as they evolve, through segmented exhaustive coverage. That coverage maps the branch structure, and the structure is
pruned on what it shows. The aim is to map the choice phase of the choice-force dichotomy, where the boundary board state is the
boundary. This is determined through perturbation and mining, which map out the pruned structure, and through discovering the
primitive race conditions by segmented exhaustive coverage. The structure is pruned on the certainties of those race conditions.

**Temporal analytical relativistic measuring.** The designer's coined term. It is not jargon; it means something specific. Take a
position A, within the pruned structure, with White to move, and suppose it is known that from this position Black has the
advantage. If, down the line, the same position is reached with Black to move, then it can be determined from that position and the
downstream primitive calculations what the difference in tempo does, including its impact on zugzwang situations. If it is the same
position with the same tempo (the same side to move), caching can be applied: simply a routing to the structure that has been
mapped out, which is symmetrically the same because the downstream position is the same.

**Reading the structure as it is mapped.** The primitives are calculated and understood as the structure is mapped, and a race
condition closes in one of two ways: by a decisive move, which can be pruned, or because it is forced by the move pool decision.
Whichever it is, it is informative on how to move forward. Pruning a position because it is evidently not a viable path to the
equilibrium is a route. If the closing rests on move pool reduction and essentially forced calculations, then relativistic
comparison finds whether this is the actual equilibrium, or whether none exists, exhaustively.

**As stated, unfinished.** The designer's last sentence on the third part was left trailing ("exhaustively....."). It is recorded as
given, and no completion of it is supplied here.

**This file's reading, offered for correction.** Same position with the same side to move is already one cache entry (§6,
convergence). What §4c adds, as read here, is the comparison of the same position across a different side to move. How that
comparison is recorded (a new cache key, or a relation between two entries) is not specified, and nothing is invented for it.
Position A's "advantage" is read as a value from the pruning (§4b) and not as a heuristic.

### 4d. Symmetry emerges from the pruned structure (the designer's statement)

The symmetry emerges from the structure when it is pruned sufficiently. What is pruned are the decisive forced positions, whether
in favor of the side or not, with the primitives understood from the pruning. What remains are the positions that are not
decisive, and within forced sequences those where draws are still possible, where forced wins are not reachable, and where the
only forced losses come from deviating from the drawing choice. That is where the symmetry emerges. If the decisive and forced
positions are all pruned, so that the entire structure is pruned, there is no equilibrium boundary board state position. If
something remains, it can be analyzed.

**The designer's clarification.** The symmetry is a measurable emergent property of existing at all after the pruning. It can only
exist when the following holds: there are no forced win moves, a draw is available, and a deviation from the draw moves opens the
door to a forced sequence in the losing favor.

### 4e. The roles of the two sides and the tempo of the initial position (the designer's statement)

This is grounded in the two-phase dichotomy of choice and force and in the boundary board state barrier, and not in particular
players, styles or techniques. It comes from the designer's own observation, independent of chess theory, and whether it already
exists in theory under a name (a specific opening, say) is to be noted if it does.

**The tempo of the initial position.** The survey anchors at the initial position, as the forward trace of a weak solution. A
subtlety beneath this is the tempo of the initial position: White has it. It is a bounded structural fact about the initial position
of chess, and it is informative. Assuming the initial position cannot be won by force from the first move, the two phases exist.

**Tempo, as defined (the designer's statement).** Tempo is not a loosely understood term here. It has a direct link to the program in
explicit terms, and it is an emergent property of the whole structure. From the initial position White is first to move, and this is
taken as a temporal advantage to White relative to Black, from the anchoring position of the initial game start. The calculation of
tempo, without exhaustive coverage of the combinatorial branching, is unlikely to be well understood, if at all, and this is where
the program can bring an explicit advantage. The balance between constructor and reactor has its anchoring point of equilibrium
where both sides have the same tempo relative to the primitive race condition calculations. Those are the boundary board state
positions, where the move pool is only draw moves or loss moves for each of White and Black. In a forced phase, the constructor (the
side with the tempo) forces the outcome by winning on a race condition primitive. In the choice phase, the constructor is
constructing the path toward the boundary board state position with the tempo, and Black is reacting to the constructed structure
of White, toward a limited pool of boundary board state positions under optimal play. Tempo in this sense is who is in the
driver's seat under optimal play, for both White and Black. Entire positions can be differentiated as win, draw or loss by whose
turn it is, which is a temporal differentiation marker. Tempo can therefore be determined exhaustively, or from the board
position, whose turn it is, and the exhaustive joint compatible primitive calculations, which show who has the tempo under optimal
conditions. In the choice phase the reason is not that the joint compatible exhaustive primitive calculations are too
combinatorially large. It is that any primitive race condition for a forcing sequence can be mitigated through a forced draw under
optimal conditions. Tempo was always understood explicitly in this way. It is not a sense or a heuristic: it is a first-principles
calculable metric that can be defined precisely with the program.

**Tempo and its anchor (the designer's statement).** Tempo on a given position is definitionally coupled to the initial anchoring measurement point. If the tempo relation can be
understood from arriving at the position from any route, that tempo is the same calculation regardless of route, by the definition of
tempo. To discover tempo, the calculations of tempo to arrive from the initial position must be discovered. It is definitionally
prohibited to understand tempo at all without this anchor, or without the anchor of the terminal game condition (the equilibrium, or
the forced sequence, anchoring position), to understand the temporal relative value in the first place. The designer adds that the
confusion about this is understandable and should be informed by the program itself.

**The roles under optimal play.** With the global boundary board state structure, White and Black have distinct roles. White is the
constructor and Black is the reactor. White chooses which opening line to navigate, the direction in which to arrive at the
boundary board state position where the game draws. Black reacts to those choices, from a pool of choices that White's choices have
set. White makes all of the choice, and Black makes choices within the structural composition White has started. White also has a
choice, under optimal play, to maintain the tempo or to give Black the tempo. Optimal play to arrive at a boundary board state
position may not depend on maintaining the tempo, but it is a choice that Black cannot make and White makes.

**A criterion to scope the search.** Maintaining the tempo when one has it is, under most situations, likely optimal for winning
conditions such as forcing a sequential mate. So, for discovery, branches that surrender the tempo can be pruned. Boundary board
state positions that maintain White's momentum from the beginning of the game are then found, and positions that give the tempo are
entirely pruned out of the discovery mining.

**The observation of Black's schematic.** The fianchetto is seen more often played by Black, unless White chooses a line in which the
fianchetto was part of deep study. In speed chess, Black appears to choose a fianchetto schematic almost as an automatic move
scheme: a fianchetto on the king side with a pawn push (typically the e pawn), or a pawn push reactive to White's development of
that pawn. The general schematic of Black's opening is: White moves, b6; White moves, either Bb7 or the pawn push; White moves.
The designer reads this as informative because of the dichotomy: Black is reactive to White and White is the constructor, at the
highest level and under optimal settings.

How the designer's certainty arises from this observation is in `boundary_board_states_start_selection.md` §2.

**Why, from the primitives.** The combination g6, Bg7, d6 makes progress on the development race primitives in a particular way,
anchored on the initial position. g6 is a bishop development, Bg7 is the bishop development move, and d6 is a development move that
unlocks the mobility of several pieces at once: it opens the potential of the c pawn push, with the queen ray and Black's bishop
ray, and gives the king mobility and a fortress without committing to the castling side. These moves also do not impact the other
race condition primitives negatively, and do not commit to a branching path that would make any of them a blunder from the initial
position anchor. This allows White to construct while Black makes moves that directionally resolve race condition primitive
dependencies, in reaction to White's construction. Moves that deviate and threaten the position are in a small pool. The designer
reads this as the structural result of the primitives on the causality of optimal play showing in actual play, with the program
allowing us to understand why this reactive schema is optimal and observed often.

**What it informs.** The place to start evaluating is not Bird's Opening, which is obscure, but the richest vein of the branch
network. There the pruned structure is enumerated by making White's move choices to construct against the known Black schematic,
with the pruning methodology, and possibly with the tempo constraint above. The designer regards this as a missing piece of the
protective belt, something to use in the pruning process and to understand optimal chess when it is seen.

**This file's reading, offered for correction.**

- **A prune by tempo is a prune for discovery.** The closed pruned structure of §4a requires every legal move to be inside the region
  or a certified forced loss. A branch that surrenders the tempo and is pruned only by the criterion is therefore recorded as pruned
  by scope, uncertified, and not as closed. The boundary states found under the criterion are candidates anchored on White's
  momentum. Whether a boundary state exists that does not maintain it is left open by the criterion.
- **The same holds for the fixed Black schematic.** Playing against a fixed Black schematic fixes Black's pool, so what it finds is
  a candidate line, and the closure of §2a still requires the other Black replies.
- **The names.** The setup g6, Bg7, d6 appears in opening theory in the Modern Defence, the Pirc Defence, the King's Indian Defence
  and the Dragon variation of the Sicilian, among others; b6 with Bb7 appears in Owen's Defence and the Queen's Indian. This is
  noted from memory and not checked here.
- **The belt.** The designer's placing of this in the protective belt has not been assessed against the admission rule (catalog §1).
  It is not a catalog entry. How the tempo of the initial position relates to the tempo terms of the catalog (§4, §8) is not
  specified.
- **Tempo and the earlier statements.** The definition above is the same marker as in §4c, where the same position with a different
  side to move is compared. The tempo of the initial position (White to move) is the one tempo fact that is known before any
  calculation. For positions further down the line, the definition makes tempo a result of the primitive calculations, so the
  criterion that prunes branches surrendering it is applied where that result is available.
- **Tempo and its anchor.** The tempo terms of the catalog (TP of P-003, TM of P-008) compare the lengths of races. Read with the
  statement above, a race length is measured to a terminal condition, the completion of a net or a promotion or the equilibrium, so
  that is the second kind of anchor, and the initial position is the first. This partly answers the open relation noted under "The
  belt" above. It is a reading, offered for correction. For a seed reached by construction (§4j), the initial-position anchor is not
  used and the terminal anchor is.

### 4f. Composition, what is cheap, and reachability from certified positions (the designer's statement)

**Composition.** Primitive calculation is not just knowing the primitives and calculating them. It is how to identify properly how
the primitives compose. The primitives are piece agnostic and define race conditions for a goal. How those goals depend on each
other compositionally, in a segmented way, and how the primitives are composed to create the calculation, is another matter
entirely.

**Composition beyond the three primitives (the designer's statement).** Composition beyond the three primitives is understood through joint compatibility equilibrium analysis after the confirmation of a
boundary board state position. Either the primitives cover all the primitives or they do not, and the ones that do not are primitives
that need to be added. There is no uncertainty here.

**What is cheap.** The only cheap things are: the board position; whose turn it is; the primitive catalog; and the boundary board
state positions that have been validated and certified.

**Closure and certification.** If closure is shown from the initial game position to the boundary board state positions, the weak
solution is provided. Proving the boundary board state positions themselves is a matter of exhaustive coverage from mate, which
hopefully resides in the known 7-piece Syzygy tablesets.

**Reachability from certified positions.** Once a proven certified position, or set of positions, exists, work proceeds on optimal
principles of reachability: a move before a certified position must reside within the closure. It is then a matter of determining
where that closure is, and the boundary line that must exist from the initial position to whatever boundary exists. As the catalog
grows, the reachability can be exploited, and other emergent relationships in the structure of the program can be used to optimize
the pruning mechanism and the program as a whole.

**Where the two directions meet, and where the work stands (the designer's statement).** The problem is scoped as a compositional
application of known primitives, and then the calculations are made with the known certified boundary board state positions, to
scope from the boundary positions when the closure of the boundary catalog is complete. Until then the work is to enumerate that
catalog and certify its positions, to prove the closure; then the tempo can be calculated explicitly. From that point the
Syzygy forward and backward trace and scope are reduced, from the boundary line of the catalog to the mating positions, and forward
trace from that boundary may help solve the backtrace. The two directions therefore meet at the boundary line of the catalog.

**This file's reading, offered for correction.**

- **Tablebase scope.** A tablebase certifies a boundary state directly only if the position has seven pieces or fewer. For a state
  with more pieces, the certificate is the complete case list of §4a, and the tablebase can certify the leaves of the refutation
  lines where they reach seven pieces or fewer.
- **Two directions.** Reachability from certified positions is a backward direction alongside the forward trace from the opening. They
  meet at the boundary line of the catalog, as the designer states above. An earlier version of this note called the combination
  unspecified; that was a misreading.
- **Composition.** The composition beyond the three primitives is stated by the designer above. The documents hold the joint
  compatibility of §4a and the catalog's composition rules for the three primitives; a primitive that the equilibrium analysis shows
  to be missing is added by the reduction rule (topology §8x) under the admission rule (catalog §1).

**The procedure that applies to any seed.** The discovery of new local states, and the reachability from certified ones, use one
procedure that works on a decisive position or on a certified local boundary board state alike: topology §8an (the designer's
statement). The routes it opens for the first certified state are recorded in `boundary_board_states_start_selection.md` §2. It is
seeded by a position and is not anchored at the initial position. It meets the forward trace of this survey at the boundary line of
the catalog, as stated above.

**The proof object of a local state, and of an upstream position (the designer's statement).** The cataloging of a locally certified
boundary board state position is fairly simple and straightforward. The pruning structure already detailed is applied to the boundary
board state position being certified, and that structure is the proof object. Any upstream boundary board state then only has to be the
exhaustive coverage between the two proof objects.

**This file's reading, offered for correction.** The proof object of a position is the closed pruned structure on it (§4 step 6, §4a).
A certified entry in the catalog is a leaf for the structure of any position upstream of it, since a resolved position is not worked
again (§6, convergence; topology §8p). The proof object of an upstream position is therefore the exhaustive coverage between that
position and the certified entries its branches reach. The closure condition of §4a applies to it: for both sides, every branch ends
at a certified forced loss, at a certified entry, or at a move that leads to a position in the region, and not only along the path.
The pieces of this were in the record (§4 step 6, §4a, §6, topology §7 and §8p) and the statement itself was not; the gap was found in
discussion and is closed here.

### 4g. The end-state test in one place (a summary of §2a, §4b and §4d)

This is a summary and adds nothing to them; the statements are the designer's, in those sections.

- The start is a premise: drawn, because no forced win has been found and no primitive race has committed moves (§4b, §2a).
- The pruned structure is built from the opening position by segmented exhaustive coverage. Decisive forced positions are pruned,
  each with its primitive and certificate (§4a, §4c).
- What remains is tested (§4d): there are no forced win moves, a draw is available, and a deviation from the draw moves opens the door
  to a forced sequence in the losing favor. The symmetry exists only when that holds, and it is the boundary board state (§4b).
- If the entire structure is pruned, there is no equilibrium boundary board state. If something remains, it is analyzed (§4d).
- A closure of the global boundary states from the initial position, given a drawn start, is the weak solution (§2a, §4f).

### 4h. The order of work (the designer's statement)

The combinatorial explosion is scoped by finding the furthest downstream branch and certifying its boundary board state status
first, hopefully within Syzygy. Then the table is used with reachability to enumerate, from the forward trace, the other reachable
positions that are enumerated in Syzygy. Then the difference is scoped, exhaustive coverage is done from the known Syzygy
solutions, and the work builds up from there. Certified local boundary board state positions are cataloged first (§4f, and
`boundary_board_states_start_selection.md`), and the classification as global is certified afterwards.

This is a constructive process for the protective belt, so it cannot all be planned out in advance. The designer understands what is
required and the steps needed to build out, and the steps beyond the above are not written here. Leveraging the structure of the
program along the way to add to the protective belt is part of the process.

The uncertainty, as the designer states it: whether downstream solutions can be arrived at and cataloged from the forward trace
within the Syzygy limitations; whether only local boundary board state positions can be certified; and whether the work can go back
from there to verifying global status eventually. Everything other than the starting position for White (a constructor path, against
the set Black response) is, to the designer, resolved.

**This file's reading, offered for correction.**

- **It matches the scaffolding of topology §7,** which also starts from confirmed boundary states and extends upstream, recording
  where it stops. What §4h adds is that the first confirmed state is sought by the forward trace from a chosen opening.
- **"Scoping the difference"** is read as the positions the forward trace reaches that Syzygy does not resolve, where exhaustive
  coverage is then done from the known Syzygy solutions. If it means something else, this reading is wrong.

### 4i. Pruning by primitive race loss: which moves force a deviation (the designer's statement)

*"Schematic" is the designer's word for an opening framed as an optimal opening for one side. It is not defined elsewhere in the record,
and no definition is added here.*

The operation: take a planned opening (here Black's g6, then Bg7 or d6 depending on White's move, then the remaining move) and ask
which moves by the other side invalidate it and force a deviation, by making a move of the plan lose a primitive race condition. Why
such a move works is informative for the pruning and for how segmented race condition dependencies compose.

Stepwise: look at what the plan's first move changes. g6 opens the dark-squared diagonal across the board to the rook, the ray Bg7 is
meant to occupy, and the plan's race needs two moves for it. Ask which of the opponent's pieces can reach that line within the race's
move count: the queen, knights and rooks cannot in time; pawn structure can; the dark-squared bishop can after a b-pawn move. Two
against two, with White moving first, Black loses the race. This is also compared against the sum of Black's responses that deviate
from the plan: if g6 is not yet committed, Black can abandon it. So the optimality of the disrupting move must itself be compared
against the other moves. This eliminates the impact of any piece move besides pawn moves for the first move that deviates Black from
the plan.

**The protocol it follows.** The anchor is the initial position, whose piece structure and race conditions are a constant in every
game. Committed moves decisively change race conditions, and the calculation runs continuously from the anchor, picking up after each
perturbation. One counter-example beating a race of a plan invalidates it as optimal play. The designer states this generalizes far
beyond the example.

### 4j. Two catalogs, and perturbation in both settings (the designer's statement)

**The two catalogs.** There are two catalogs: the primitive catalog and a boundary board state certification catalog. Both are important. The primitive
catalog is important for understanding primitive composition, which dictates the calculations. The catalog of boundary board state
positions and reachability families is important to eventually enumerate the closure and to progress from local boundary board states
to families of reachable local boundaries, building up from there.

**Perturbation, in both settings.** Perturbation is the analysis of a board state via move changes, whether by changing piece position entirely or by movement up a
branch through legal move requirements. The change in primitive state caused by a move applies to moving up the tree to find where a
race was acted upon, through the constraint of legal moves. When surveying to discover local boundary board states through
construction, the legal move constraint is not required, since how to arrive at the position is not a concern and there is no
consideration of the initial game start anchor. Perturbation in that context allows a change in board state in piece position without
legal move restriction.

**The scope of a perturbation (the designer's statement).** The perturbation scope is dictated by the intended purpose of the perturbation and by whether it is anchored or not in the first
place. Perturbation is to be understood for what it is, a change of state, whether in discovery or in following a branch.

**This file's reading, offered for correction.**

- **The names.** The candidate catalog of §3, and the candidate, proven and disproven catalogs of `boundary_board_states.md` §3, are
  read as the stages of the boundary board state certification catalog. The reachability family is added to the objects of §3 in the
  designer's phrase, with the wording of the row offered for correction.
- **Dead pieces.** Pieces that can be removed without changing the race are the pieces outside the load-bearing set and the footprint
  that the record already has: the ablation of topology §8s (catalog §4.12), the footprint (catalog §4.9, topology §8h), and the
  predictions P3-1 and P8-1. "Unless they can obstruct" is the footprint's membership. This is a reading, offered for correction.

### 4k. Mining the Syzygy tablebases for the certification catalog (the designer's statements)

The full statements and the process are in `boundary_board_states_syzygy_mining.md` §1 and §3; they are summarized here.

- The aim is to begin the boundary board state certification catalog (§4j) by creating it from the proven results of the Syzygy
  tablebases, enumerating the boundary board states from the solution set. A branch can be enumerated by mapping the tablebase value
  (and distance) onto the branch with its legal move pool. Branches are created under optimal constraints.
- Primitive calculations are not part of the script. They are done when a drawn position is verified later, and primitive
  compositional analysis is done after the fact, from the population of certified local boundary board states obtained by
  exhaustive coverage from Syzygy.
- The 50-move rule is ignored for now. It matters only in a reachable family, as a cutoff, when the anchor is reached.
- Draws are categorized by stalemate, insufficient material, or cyclical. Where the best move verified from the tablebase is to
  reduce to insufficient material, that is a provable equilibrium of reduction rather than a stalemate choice, since optimality is
  the concern.
- A position in the choice region with a move that can lose is a boundary, and the initial position is one (§5a). The choice space
  and the forced space together encapsulate the entire space of chess, with sub-optimal areas pruned, and the pruned structure is
  the full solution. The designer wants to confirm that boundary board states are a majority of the space.
- The tablebase size does not matter eventually, for reachability, with the program's methodology of backtracing and exhaustive
  enumeration from the cached database. The work is to mine all the Syzygy data the machine can handle and then automate off that
  basis with the pruning mechanism. The order is mining, then the catalog and families from it, then automation and the matching of
  primitives.

**This file's reading, offered for correction.**

- **The tablebase as the certificate of §4a.** For a drawn position of a material held in the tablebase, every legal move either
  keeps the draw or loses, and the table gives the value of each, so the closure condition of §4a (every legal move accounted for)
  is met by the table. Whether each side has a losing deviation at the position, and what "both sides live" tests, is not settled by
  the statements (mining document §9).
- **The mined population as the entry stage of the catalog.** The proven draws of a mined material are read as entries of the
  certification catalog at its entry stage (§4j), for the materials mined. No fields of an entry are proposed.
- **The order of work.** This is §4h run on tablebase data: the furthest downstream certified entries first, and the families and the
  upward link afterwards (§5a).

## 5. Local and global along the line

A local candidate verified at depth d is local until its family is shown to be reachable from the initial position under the
optimal backtrace constraint (§5a). Showing that is a matter of the positions from which the candidate can be reached by optimal
play from the anchor, with transpositions included, so it is a web and not a tree. It can be large, and the choice of depth keeps
it small enough to examine (§7).

Avoidance is not the test. A side may have a non-losing detour that never reaches the candidate, and the detour leads to another
boundary state or to a position of the same kind. The candidate still closes the phase for the paths that arrive at it, and the
detour does not make it non-global. (A reading that treated avoidability as a ground for non-global was proposed in discussion and
withdrawn on the designer's correction.)

A consequence is worth recording, the designer's reading:

- If the survey concludes exhaustively that no local candidate can exist below a node, the work is done once and need not be
  repeated for any path that reaches that node (§6, C3).

## 5a. Reachability families and the classification of a family (the designer's statement)

Reachability is simple: if a boundary board state can reach another one, then it is in the same reachable family. If that family is
confirmed joint to the initial start position, that is an important certification, and it essentially classifies all of the reachable
positions in that family as certified global boundary board state positions, when reaching the initial position under optimal
backtrace constraints. This is what makes reachable families important: they are a classification methodology to resolve a family of
reachable positions when that family is eventually anchored with the initial start position.

**Global, restated (the designer's statement).** Any upstream position that under optimal conditions leads to a boundary lies in the choice region, and the boundary actually being realized does not change the equilibrium that the boundary position and the upstream positions reside within, if they are part of a reachable family under the optimal backtrace constraint. A family that is reachable from the initial start position under that constraint is an optimally reachable boundary board state family, and every member of it is a global boundary board state position. The initial game start position is itself a boundary board state position: it resides in the choice region, and if it has a move that can lose, it is a boundary. The designer adds that a distinction between the choice region and the boundary drawn by how many drawing moves a position has is a meaningless classification criterion.

**From discovery to the initial position (the designer's statement).** When discovering a local boundary board state from construction, the constraint of backtrace does not apply: discovery has no
constraint to the initial position anchor. Eventually a reachable family is composed that either concludes it cannot reach the initial
position under the optimal reachability constraint, or concludes that it does. This covers the whole process from discovery to
linking to the initial position. A family must reach the initial game start to prove global, if that family is composed under the
optimality backtrace constraint, which is the requirement for any family created.

**This file's reading, offered for correction.**

- **With §2 and §5.** The statement above is the definition of *global* used in §2 and §5. Offered for correction.
- **With topology §7.** The connected regions of resolved territory there, two confirmed states joined by a validated chain, are read
  as related to a reachability family. Whether they are the same is not stated.
- **The flow.** Three stages, as read from the statement above. Discovery, by construction or otherwise, finds a certified local state
  with no backtrace constraint. A family is then composed from it under the optimal backtrace constraint. The family either reaches the
  initial position, which classifies it as in §5a, or concludes that it cannot, in which case it is not global and its states stay
  certified local. The final clause of the statement was put as a question; it is read as saying that the optimal backtrace constraint
  applies to every family. Offered for correction.

## 6. The cache

The cache is keyed on the position. Its claims concern every path from the anchor that arrives at that position, so a position
reached by another route, or reached again, inherits the entry. It follows the architecture of `boundary_board_states.md` §3
(a verdict, once resolved, is permanent; only finished verdicts are persisted).

| Entry | Content | What it saves |
|---|---|---|
| **C1 Forced region** | From this position every branch concludes in a forced sequential win for one side, and the only decisive moves are forced out of zugzwang-type situations, so the position is classified as forced up to the point where a non-forced deviation could occur. No boundary state exists downstream under optimal play. | Once any path reaches the position, the work below it is already done |
| **C2 Boundary found** | A local candidate, or a verified local state, exists here. It may or may not be global, and the entry records which has been shown. | A found boundary is not rediscovered, and later classification reuses it |
| **C3 None can exist** | It was concluded exhaustively that no local candidate exists below this position, on the closure of the non-losing moves and to the depth reached (scope below). | No future work is spent proving it again |

**What an entry must carry.** The position key and side to move; the claim (C1, C2 or C3); the certificate that supports it (a
complete case list, exhaustive search, or an assertion taken from theory, with its evidence tag); the scope of the claim (below);
and the depth and line through which it was reached.

**Guard rail.** A claim of C1 or C3 is a claim about everything downstream. It needs a certificate. The failure of a search to find
a boundary state below a node is absence of information and not a result (topology §8p). An entry that rests on an assertion in
the theory tree is tagged as such, and it does not count as certified.

**Scope of C1.** The claim is made *under optimal play*, as the designer stated. A forced win has no drawn position on any optimal
continuation, so no boundary state can lie on one. A drawn position that the winning side could reach by an error is not on an
optimal continuation, and under the target of §2 it is an admissible local candidate anchored through this line. C1 therefore says
nothing about such positions. Whether the survey wants candidates reached by deviations from a forced win was open (§10); the designer's later statement is that such regions exist.

**Scope of C3 and the pruning of sub-optimal paths (the designer's statement).** Sub-optimal paths are excluded by the distinction
between choice and force: a choice that leads to a forced loss is pruned. In a drawn region a move either preserves the value or
loses, so the non-losing moves of both sides are all optimal in value, and a move that is only inferior in play stays in the net.
The net the survey exhausts is therefore the closure of the non-losing moves, the chain of topology §5. C3 claims that no
boundary state exists on that closure below the position, to the depth reached. It makes no claim about paths that leave the closure
through a losing move. Each prune rests on a certified forced loss, and a claim that rests on an uncertified prune is conditional
and is tagged as such.

**Convergence.** Positions that recur or converge, by transposition or by a repetition cycle, are one entry in the cache. The
region below a node is a graph and not a tree, and a position already resolved is not worked again. This reduces the work, and it is
the reason the cache is keyed on the position.

## 7. Calibration: choosing a line and a depth

The survey fails in a specific way: the candidates it finds are not independent. If the forced-draw nodes at a depth form a wall,
a set of boundary-state candidates that reach and sustain one another in a cycle, the verification of each depends on the others,
and the tails of the cycle expand the scope of proof, because each candidate has interdependencies that branch into the others'
spaces. When the web branches into a combinatorial space too large to verify exhaustively, the choice of path is a bad one. The
line is too far upstream for its candidacy to be examined.

This is a calibration problem. The depth and the line are chosen so that the web of candidates and the footprint closure (topology
§8h, four conditions) stay of a size that the complete case list can cover. Signals that the choice is wrong:

- The footprint of the candidate's primitives grows to the board (kings as the global resource, §8h condition 3).
- The positions from which a global classification would have to link the candidate to the initial position (§5) are not enumerable.
- The candidates at the depth chosen are interdependent, and no candidate's verification can be separated from the rest.

The designer's statement on where to start (§4e) bears on the choice of line. The progress of that discovery is recorded in
`boundary_board_states_start_selection.md`.

When a signal fires, the line is changed or the depth is moved, and the result is recorded as a cache entry for the abandoned
choice (a note that the web was too large at this depth, which is information and not a certificate).

## 7a. The segmented progression (what is understood so far)

This section records only what the designer has stated. The full process is not yet specified and is not written here.

**The designer's statements.** The space from the initial position to the global boundary states is combinatorially vast. It is
mitigated by choosing a depth from the initial position at which the primitives scope the pruning properly (§7). The survey
proceeds in segments. A segment may show that a local candidate was already passed, which is information: it means that a boundary
state exists on those paths. That is cached as a checkpoint, and the progression
continues from it. The classification of a local candidate as global, or as provably not global, is settled upstream later, by
the progress of the segments, the catalog and the cache. A file for the segmented progression, to enumerate the landscape in
this way, is to be started.

**Segments are compared relativistically (the designer's statement).** The work is to reduce the move pool in the pruned structures by
relativistic perturbation branching and by the primitive races, through how a move opens the door to a race imbalance against the side
that makes it, done in parity for the opponent. This is more practical than enumerating full coverage until mate, because it reduces a
branch structure. The branching structure is worked in segments. A segment added to the pruned structure either wins its race against
another segment, or is pruned entirely because what already exists in the structure shows that the branch would never be taken.

**This file's reading, offered for correction.** The comparison of segments is a comparison of values, and a branch is excluded when
the side to move has an alternative of value at least as good. Each such exclusion carries a certificate like any other prune (§4a).

The whole process as the designer has stated it is in §4c.

**Not yet specified.** The fields of the progression file; the order in which segments are taken; how the depth of a segment is
chosen beyond the calibration signals of §7; how the boundary conditions of a segment at its depth limit are recorded; and how
segments compose into the closure of §2a. These are not decided and are not filled in here.

## 8. The theory tree as material

Opening theory is a tree of non-losing moves, with the deviations of each side annotated as losing. For the survey this has three
uses, none of which is a certificate:

- **Deviation map.** The refutations in the book are the raw material for the deviation map of a proven entry (topology §8x): each
  losing deviation of each side is refuted by a named race condition, and the book's refutation is the first candidate for it.
- **Width.** The number of equal alternatives at a node is its width, recorded as in `boundary_board_states.md` §14 item 8.
  Recording the width at each node of the line gives a profile along it. It shows the narrowing from choice to force as data. This
  profile is this file's suggestion (§11).
- **Evidence tier.** A book's assessment, or an engine's, is an assertion. It places a candidate in the candidate catalog with the
  tag INSTANCES. It becomes CERTIFIED only through the complete case list or exhaustive verification (topology §8r).

## 9. What the survey produces

For each line surveyed: the nodes classified (decisive, choice, forced-draw); the located boundary between choice and force with
the primitives that carry it; candidate entries with their primitive accounting; the web, if any, and the calibration judgment; the
cache entries; and the primitives that the line shows to be missing, found by the reduction rule on the residual interface
(topology §8x). The last item feeds the protective belt, subject to the catalog's admission rule (catalog §1).

## 10. Open questions

- The meaning of "first" over all paths: retired. The designer's later statement (§5a) defines *global* by membership of a family reachable from the initial position under the optimal backtrace constraint, so the question no longer arises.
- Whether C1 should also claim anything about positions reached by deviations from a forced win (§6). Answered in part by the designer's later
  statement (`boundary_board_states_start_selection.md` §2b): a decisive position can have a choice region below it, reached when the side that
  is winning plays sub-optimally or through ignorance of optimality, and a localized choice region can lie past the global one. Such regions
  exist and are what the mining looks for. Whether C1 records them is not specified.
- How a cycle of interdependent candidates is verified when it is small enough to verify at all: whether the web is verified as
  one unit.
- How the cache treats the repetition cycle that realizes a candidate, since the draw is a property of a cycle of positions and not
  of a single one. (The question of whether a forced terminal draw also realizes a candidate was ruled on by the designer and is closed:
  `boundary_board_states.md` §14 item 10.)
- Whether the width profile along a line tracks the choice-to-force boundary that the primitives locate. A prediction to be stated
  before it is run.
- How a book's assertion enters the case list where the position is outside the book's lines.
- Which line to start with, given §4e: the designer's statement is to start with the richest vein and not with Bird's Opening. The
  earlier plan for Bird's is not withdrawn by the documents, and the choice between them is the designer's.
- The sense of tempo in §4e: answered by the designer (§4e, "Tempo, as defined"). It is the explicit, calculable metric stated there,
  and not a sense or a heuristic. How its use as a scoping criterion is recorded on the pruned edges is not specified.
- How primitives are composed to create the calculation (§4f): answered by the designer (§4f, "Composition beyond the three
  primitives"). It is understood through joint compatibility equilibrium analysis after a boundary board state position is confirmed.

**The designer on what waits (a statement).** The comparison across tempo does not apply until there is enough exhaustive
certainty from both the local and the global catalogs of boundary board state positions. The compositional question is a matter of
finding a first example. The cache key for a repetition cycle is not something the designer recognizes, and is assumed to be
discovered. Pruning relative to temporal decisions is a modular addition to the data, to be discovered once the prerequisites for
making those decisions exist. These are therefore not prerequisites for the first certified local boundary board state.

## 11. Provenance

**The designer's, and unchanged:** that the initial position is the anchoring point of every game; that a global boundary state
is defined by membership of a family reachable from it under the optimal backtrace constraint (§5a), and so depends on how it is reached, and is
sound against all deviation replies, there being many global states, each of which closes the phase for its own routes; that a local one
is a verified equilibrium with no requirement on the path or on its family reaching the initial position; that the survey's target is a local candidate
anchored from the initial position by a chosen opening, optimal or not; that forced sequences are read as race conditions, that a
forced win is pruned by adversarial play, that a forced loss is uninformative and a forced draw is the boundary-state behavior;
that the theory tree and the perturbation process locate where choice becomes force; that the cache records, by position, that
force is already deployable from there, that a boundary was found, or that none can exist, so the work is not repeated; that the
depth is a calibration problem, and that a web of interdependent candidates marks a choice that is too far upstream; that
the closure of the global boundary states must hold (§2a), with the outcomes and the premise as stated there; that the progression
is segmented (§7a, as far as it is understood); and that a
survey of Modern Chess Openings is the way to select the lines.

**This file's organizational choices, offered for the designer to change:** the separation into definitions, objects, procedure,
cache and calibration; the node classes of §3; the three cache entries and their fields (the designer described their content; the
names C1 to C3 are this file's); the reading of *global* against the start-reachability of topology §5 (§2); the scope note on C1
(§6); the width profile (§8); and the list of open questions.

**A correction recorded.** In discussion, the global and local classification was first read as a scope tag attached to a
verdict, on the ground that a position's status does not depend on how it is reached. The designer corrected this: the status of a
boundary state is the position's, and the classification as global or local is path-dependent by definition. The first reading was withdrawn and none of it remains here.

**A later statement recorded.** The designer's later statement in §5a defines *global* by membership of a family reachable from the initial position under the optimal backtrace constraint, with the initial position itself a boundary board state position. It replaces the earlier first-boundary definition (global as the first boundary across all paths, with a second proof that no boundary state lies between the initial position and it), which the designer did not recognize as their understanding. That definition and the material resting on it were removed from this file. The history is in the pipeline file's changelog (entries 31 and 60 and 61).
