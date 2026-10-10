# Boundary Board States: the Anchored Opening Survey

**Status.** A method document, not a program and not a tool. It specifies how the primitive catalog is used to look for a
boundary-state candidate along a chosen opening line, anchored at the initial game position. Nothing in it has been run. The
hard core (the boundary board state) and the protective belt (the three primitives) are unchanged by it. Evidence tags as in the
catalog: **SPEC** (defined, not tested), **INSTANCES** (checked on specific positions), **CERTIFIED** (backed by an exhaustive or
complete case-list certificate).

Companion documents: `boundary_board_states.md` (the concept and methods; §1 the definition, §3 the caching architecture, §14
items 4, 7, 8 and 9), `boundary_board_state_topology.md` (the record; §3 status is local and unconditional, §5 start-reachability,
§8h footprints, §8p the guard rail on absent paths, §8r forward anchoring, §8x the proven entry, §8ac and §8ae the designer's statements
behind this file), `boundary_board_state_primitive_catalog.md` (the three primitives, the admission rule),
`boundary_board_states_case_study_otb.md` (the game).

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
  and force. It informs us of the viability of the primitives, and it is the point on which exhaustive verification is performed.
- **Global boundary board state.** A boundary state that is the first boundary between choice and force from the initial game
  position. It is path-dependent by definition, because it closes the phase between choice and force from the start to itself.
  Classifying a state as global requires proof of two things: that it is a boundary state, and that no boundary state exists
  between the initial position and it. A state with the first proof and not the second is local.
- **The paths are all paths.** The web of dependencies among candidates shows that the second proof concerns every way of arriving
  at the position, not one line. "First" is across all paths that reach the position from the initial position.
- **Global is sound against every deviation, and there are many.** A global boundary state must hold holistically: it is reached by
  non-losing play, and it holds against all deviation replies. A deviation from the path may lead to a different global boundary
  state, a detour to another position that closes the phase between choice and force. That does not touch the first. Each global
  state closes the phase for the routes that lead to it, so the global states form a list, and the space between the initial
  position and that list is the enclosed phasing of choice and force. That a position can be avoided is therefore not a ground
  for classifying it as non-global. The only ground is a boundary state that precedes it on a non-losing path that reaches it.
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
  drawn start (topology §5), and the collective version is topology §3. The closure claim adds that each such transition is a boundary
  state in the sense of `boundary_board_states.md` §14 items 7, 8 and 10. That addition is the programme's claim, and the enumeration
  tests it.
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
6. **Enter and verify.** A candidate is entered in the candidate catalog with its composition. Its verification is the complete
   case list of topology §8d with exhaustive insurance, scoped by the footprint of §8h. A verified candidate is a local boundary
   state. Whether it is global is a separate question (§5).
7. **Cache** what was learned, whichever way the search ends (§6).

## 5. Local and global along the line

A local candidate verified at depth d says nothing about whether an earlier boundary exists. To classify it as global, the survey
must show that no boundary state exists in the positions that precede it on any path. That is the set of positions from which the
candidate can be reached by legal play from the anchor, with transpositions included, so it is a web and not a tree. It can be
large, and the choice of depth keeps it small enough to examine (§7).

Avoidance is not the test. A side may have a non-losing detour that never reaches the candidate, and the detour leads to another
boundary state or to a position of the same kind. The candidate still closes the phase for the paths that arrive at it. Only a
boundary state earlier on a non-losing path to the candidate makes it non-global. (A reading that treated avoidability as a second
ground for non-global was proposed in discussion and withdrawn on the designer's correction.)

Two consequences are worth recording, both the designer's reading:

- If a local state is found at some depth, then no later state on a path through it is global, since it is not first. That is a
  fact about every position that can only be reached through it. It is recorded as a cache entry (§6, C2).
- If the survey concludes exhaustively that no local candidate can exist below a node, the work is done once and need not be
  repeated for any path that reaches that node (§6, C3).

## 6. The cache

The cache is keyed on the position. Its claims concern every path from the anchor that arrives at that position, so a position
reached by another route, or reached again, inherits the entry. It follows the architecture of `boundary_board_states.md` §3
(a verdict, once resolved, is permanent; only finished verdicts are persisted).

| Entry | Content | What it saves |
|---|---|---|
| **C1 Forced region** | From this position every branch concludes in a forced sequential win for one side, and the only decisive moves are forced out of zugzwang-type situations, so the position is classified as forced up to the point where a non-forced deviation could occur. No boundary state exists downstream under optimal play. | Once any path reaches the position, the work below it is already done |
| **C2 Boundary found** | A local candidate, or a verified local state, exists here. It may or may not be global, and the entry records which has been shown. Positions that can only be reached through it are not global. | A found boundary is not rediscovered, and later classification reuses it |
| **C3 None can exist** | It was concluded exhaustively that no local candidate exists below this position. | No future work is spent proving it again |

**What an entry must carry.** The position key and side to move; the claim (C1, C2 or C3); the certificate that supports it (a
complete case list, exhaustive search, or an assertion taken from theory, with its evidence tag); the scope of the claim (below);
and the depth and line through which it was reached.

**Guard rail.** A claim of C1 or C3 is a claim about everything downstream. It needs a certificate. The failure of a search to find
a boundary state below a node is absence of information and not a result (topology §8p). An entry that rests on an assertion in
the theory tree is tagged as such, and it does not count as certified.

**Scope of C1.** The claim is made *under optimal play*, as the designer stated. A forced win has no drawn position on any optimal
continuation, so no boundary state can lie on one. A drawn position that the winning side could reach by an error is not on an
optimal continuation, and under the target of §2 it is an admissible local candidate anchored through this line. C1 therefore says
nothing about such positions. Whether the survey wants candidates reached by deviations from a forced win is open (§10).

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
- The ancestor set that a global classification would have to check is not enumerable.
- The candidates at the depth chosen are interdependent, and no candidate's verification can be separated from the rest.

When a signal fires, the line is changed or the depth is moved, and the result is recorded as a cache entry for the abandoned
choice (a note that the web was too large at this depth, which is information and not a certificate).

## 7a. The segmented progression (what is understood so far)

This section records only what the designer has stated. The full process is not yet specified and is not written here.

**The designer's statements.** The space from the initial position to the global boundary states is combinatorially vast. It is
mitigated by choosing a depth from the initial position at which the primitives scope the pruning properly (§7). The survey
proceeds in segments. A segment may show that a local candidate was already passed, which is information: it means that a boundary
state exists on those paths and that nothing below it on them can be global. That is cached as a checkpoint, and the progression
continues from it. The classification of a local candidate as global, or as provably not global, is settled upstream later, by
the progress of the segments, the catalog and the cache. A file for the segmented progression, to enumerate the landscape in
this way, is to be started.

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

- The exact meaning of "first" over all paths: the ancestor set to be checked, and whether it is the set of positions reachable by
  non-losing play only (the §5 chain) or by all legal play. The reading of §2 takes the former for *global* and leaves the latter
  to local candidates.
- Whether C1 should also claim anything about positions reached by deviations from a forced win (§6).
- How a cycle of interdependent candidates is verified when it is small enough to verify at all: whether the web is verified as
  one unit.
- How the cache treats the repetition cycle that realizes a candidate, since the draw is a property of a cycle of positions and not
  of a single one. (The question of whether a forced terminal draw also realizes a candidate was ruled on by the designer and is closed:
  `boundary_board_states.md` §14 item 10.)
- Whether the width profile along a line tracks the choice-to-force boundary that the primitives locate. A prediction to be stated
  before it is run.
- How a book's assertion enters the case list where the position is outside the book's lines.

## 11. Provenance

**The designer's, and unchanged:** that the initial position is the anchoring point of every game; that a global boundary state
is the first boundary between choice and force from it, and so depends on how it is reached, across all paths that reach it, and is
sound against all deviation replies, there being many global states, each of which closes the phase for its own routes; that a local one
is a verified equilibrium with no requirement on the path or on being first; that the survey's target is a local candidate
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
boundary state is the position's, and the classification as global or local is path-dependent by definition, because global means
first. The first reading was withdrawn and none of it remains here.
