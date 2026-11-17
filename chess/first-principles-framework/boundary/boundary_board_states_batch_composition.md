# Segmented Batch Composition within the Paradigm Window: A Proposed Direction

**Status.** A record of a direction the designer articulated on 15 November 2026. The designer states that the batch composition is one approach among others, a modular component of the programme that applies to this part of it (certifying a segment within a storage window), and not the lower-level architecture of the Lakatosian programme. Other approaches can be added. d-depth scoping is not part of what is one approach: the designer states that it is a structural requirement of the programme, invariant to it and integrated into its fetching mechanism (`boundary_board_states.md`, the statement on d-depth scoping; §1.15 here records the modes). It is not worked out, and no result of any run is recorded in it. It is
kept as its own file so that the other documents are not changed by it. If the direction fails, it is removed by deleting this file. One
tension with an existing statement is recorded in §3 and is not resolved here. The direction is stated in general terms: no piece, material or
position is part of it, and the fractions of the window in §1.7 are an illustration and not a measurement. It adds no primitive and does not
touch the hard core (`boundary_board_states.md` §14 items 7 to 11). Statements marked **the designer's** are the designer's. Statements marked
**this file's reading** are readings of the designer's words, offered for correction.

Companion documents: `boundary_board_states_start_selection.md` (§2c the exclusion form and the depth bound d; §2d the forward anchor and the
sources of forward paths; §2e fronts and contribution; §2f seeded pruning as community work; §2g depth-scoped certification and dependency
across segments), `boundary_board_states_opening_survey.md` (§4f what the tablebase certifies; §4j the two catalogs; §5a reachability
families; §7 calibration), `boundary_board_state_primitive_catalog.md` (§8.9 the support closure; §9 the dependency map),
`boundary_board_states_syzygy_validation_and_scoping.md` (§1 the scope of the enumeration), `boundary_board_states.md` (the definition, §14, and the accounting of §15), `boundary_board_states_calibration.md` (the stored band, what absence means, the cost model).

---

## 1. The direction (the designer's)

### 1.1 What the program is without the storage limit

1. With no storage limit of the kind that bounds Syzygy, the program converges into Syzygy with the program's added benefits. Those benefits
   of the catalog have not yet been emergent or integrated, because nothing has been run.
2. The program is then Syzygy with certainty, showing full depth. The question is how to use it and scale it up progressively.

### 1.2 The depth layer and the gap between segments

1. Increase the depth d iteratively and compare the result with exhaustive coverage. For a segment with exhaustive coverage, the depth at
   which coverage resolves is sufficient to find a forced sequence if one exists, and a position where none exists is a boundary board state.
2. A position of a larger material that navigates into a smaller one takes a number of plies to do so. There is a gap between the segments,
   and it is scoped by d to arrive at the next irreversible segmented boundary, on which there is certainty at its own d.
3. This is the segmented progression and accumulation of the program without exhaustive coverage but with d-depth certification. Certainty
   from d is measured: it is the d one must calculate to in order to resolve an approach. How much d is required to realize that a move is
   optimal can be resolved.
4. This can be a public repository for chess training to which the community adds, so as to make progress in specific segments. As the
   combinatorial vastness increases, such contributions matter more. Games can be fed in to see whether any d-depth optimization is made and
   added, as a constant loop. Engines can complete games, and the program is a basis that informs decisions to progress the program
   explicitly.
5. Where Syzygy does not cover a segment, d-depth certification is a certification checkpoint of progress, and it allows that progress to be
   carried forward.

### 1.3 Wins, draws and the edges

1. Certainty accumulates from wins. Wins gain certainty, and certifying draws requires exhaustive certification.
2. The certainty from wins informs the gaps where certainty cannot be obtained because of the size that certification requires. The work is
   guided search and understanding the d of certainty that can be derived, working from the edges as games accumulate.

### 1.4 Seeds and the anchoring point

1. A seed is a game. Seeds fed into the structure enumerate certainty on carefully selected games, within the paradigm window.
2. Each member of a community holds their own repository of seeds and composes them to work towards the certification of a segment, including
   certifications of draws. The program repository does not keep seeds. Providing them would bound discovery to the discoveries reachable
   within that seed basis. The more diverse the seeds the better. What matters is making progress and understanding the impact of each game
   integrated into the seeds.
3. The anchoring point informs the seeds. Seeds converge into a transposition of an anchor point. That anchor point is the same board state
   as the place of certainty at which progress is made within the irreversible segmented boundary being certified. The convergence point
   must be the anchor point.
4. The program under-specifies the importance of the integration of the anchoring point. The seed approach functions from an anchored
   basis, and this is the prime reason for it.

### 1.5 The repository of games, and the engines

1. From a large repository of games, the games that arrive at a convergent anchoring point can be isolated and fed into the program as
   seeds. This seems an automatable process.
2. The repository has constraints. The games are fully played out, with no agreed draws and no resignations. There is no move limit for a
   draw, since engine games average over ninety moves and optimality cannot be constrained by move count.
3. Engines are added with a calibrated d-depth scope of search. The engine additions are not for certainty. Engines can search the space
   smartly within a calibrated d, especially when paired with the program and with the seeding it works from. The stronger the engine, the
   better its choices.
4. In an engine game against itself, both sides are aware of the program's progress, of the shared anchoring point of all seeds, and of the
   outcomes the seeds have built up. A choice is considered relative to that progress, including whether it follows the branch of another
   seed. The engine can then choose, from the progress already made, the deviation that would relieve the most work. This is a progressive,
   cumulative loop.

### 1.6 The gap, partial certificates and the window

1. The work then turns to certifying the gaps for draws within the reachable segments, to see whether they fit the paradigm window. How much
   is stored and how much work is spent on certification can be calibrated from the loop, and the data can be analysed to certify the gaps
   without building more storage.
2. The gap to be calculated, plus the storage from construction, is the exhaustive gap requirement. If the gap is segmented into reachable
   understanding, there is a finite uncertain region, and the optimality routing within it can be analysed to see whether routing can force
   itself into known forced regions. Partial gap certificates plus the storage from construction are then smaller than the full gap plus the
   storage from construction.
3. If a draw certificate depends on full coverage of the draw states, the partial certificates are interdependent on the completion of the
   others. A partial certificate that yields a conflict adds to the force region catalog of seeds. If the seed storage limit is hit, the
   paradigm window is hit.
4. The smaller the gap and the more iterations of this kind, the more certainty is gained. The games that cannot be certified have a scoped
   amount of uncertainty. A combination of seeds can in itself provide full coverage.
5. Suppose the forced region does not fit the storage but the uncertain region does. The full list of uncertain positions can be composed,
   and one can navigate within that uncertainty to discover whether a way can be forced into a region not on the list, which is then a force
   region. The storage limit then becomes the limit on storing all draws, and not full coverage.

### 1.7 The batch procedure (the designer's illustration, in fractions of the window)

1. Take an irreversible segmented boundary region whose full exhaustive coverage needs more than one window, say 150%. The forced region is
   70% and the choice region 80% of the window. Partial coverage is placed in seed batches that fit the window: batch a at 100% and batch b
   at 50%. A region of the segmentation that was left out cannot be reached, which gives uncertainty inherently, since resolving some
   positions needs the coverage that was segmented away.
2. Suppose batch a resolves 40% as forced under certainty and leaves 60% uncertain, and batch b resolves 20% and leaves 30% uncertain.
3. **Batch c** takes the forced region certified by a and b (40% + 20% = 60%). **Batch d** takes the uncertainty left by a and b (60% + 30% =
   90%). Of the 70% forced, 10% remains inside d.
4. **Batch e** is a smart selection from c added to d. Work on e uses the forced certainty added from c to resolve a forced branch of d. That
   branch is removed from d and placed in c. This gives more room to take from c into e, which resolves more of d.
5. Within d the goal is to force a way out of the set while keeping the advantage, arriving at an anchor point favourable to the side
   forcing. If that can be done, the result is added to the certainty region under force.
6. The room that d leaves in the window (10% in the illustration) is the room for forced-region certainty used in routing. As more is
   resolved the uncertain region is smaller, the storage requirement for certifications rises, and the depth required to force into the
   certified forced region rises. That depth is a stricter requirement on the certainty fed into the seeds. If a certification can be made
   within that room, combinations of seeding can be enumerated to obtain it up to a certain depth.

### 1.8 The two stores and the flag

1. Batch c need not be discarded. It can be an added output, and it can be stored on a different machine and imported once a smart selection
   is made. This requires managing the storage limit of the paradigm window explicitly.
2. Pruning the forced region to create this dichotomy is useful where the window does not fit the forced region. The forced region need not
   be kept, unless it requires a lot of work to certify: a region needing a depth of 15, say, should be kept over one needing 5. This is not
   what is being worked on at present.
3. The dichotomy between c and d resolves itself as a flag in itself.
4. The inverse question is whether to keep d. Batch d holds boundary board states in which the two sides are equal. The forced region
   dictates who wins or loses, which matters. A position that is not in d is not known to be winning for either side. A position that is not
   in c is known to be drawn under perfect play.

### 1.9 The terminal reason of d, the legend, and clean-up certification (the designer's, 15 November 2026)

1. **d is calculated to a stated terminal reason.** The designer's initial instinct was for d to measure the depth to the terminal mating net
   condition. The designer set that aside: it is not a scalable metric and does not inform the program in a segmented way, because it would
   require downstream calculations. Stopping the d calculation at the determination of the forced region, from an understanding of downstream
   certification, is important if the program approaches the work this way.
2. **The fit with segmented dependency.** If d is calculated to a stated terminal reason, this fits the segmented dependency primitive and the
   compositional architecture. It makes the segmented dependency explicit, and it presents an opportunity to consider scaling.
3. **An example of a terminal reason (an illustration).** The terminal condition can be a capture that places the position in the irreversible
   segmented boundary of a smaller material, with force-region certainty. The material imbalance of the position is what allows the forced
   region to be categorized. That holds unless the capture is followed by a recapture that puts the position into a boundary board state
   classification again. That case must be addressed.
4. **The terminal reason is denoted by a symbol.** A flag denotes the explanation, to save storage within the d-depth metric. As the approach
   scales, it can be addressed with a similar schema or with other applications that optimize.
5. **What is not addressed.** How this scales, in a piece-agnostic fashion, to other irreversible segmented boundaries with varying terminal
   mating positions is, in the designer's words, a different story.
6. **Clean-up certification.** A partial certification can change other calculation certainties in other areas. Clean-up certification work
   is then performed on the program, relative to the dependencies that the process altered. This is a large part of the segmented progression
   and the community effort that the designer is deliberately architecting into the program.
7. **The legend.** The symbol is integrated into the same d-depth header and is not an entry of its own, so that it adds no storage. The
   legend can be removed from the file and placed in the README of the irreversible segmented boundary itself. That README is the specified
   point for the piece-agnostic considerations and for the determination of the terminal endpoints. The legend denotes the terminal mating
   positions that are reachable, and also those left uncertain but likely forced, as a dependency certification progress point.
   The designer's statement of 16 November 2026 refines items 4 and 7: the header is [d][index][symmetry condition], and the legend is a
   file of terminal states of its own (`boundary_board_states.md`, "The d-depth metric").
8. **The earlier statement, from the same day.** The d-depth metric bounds the number of moves available to complete a terminal net. A d of 5
   leaves the winning side at most 3 moves, so the pool of nets that can complete must be reachable within that, and this guides the choices
   that keep that pool reachable. The designer's refinement of items 1 to 7 followed this statement. §1.10, items 4 to 6, speak to the bound
   again, and §2, item 18, gives this file's reading of how the two fit.
9. **Clean-up certification can be entirely automated, upon contribution.** It is a specification that requires integration and further work.
   The designer sees it as inevitable as the program scales: the work of resolving these issues will require this solution. It generalizes to
   the conditions that the program deliberately architects, which are segmented dependencies and reachability.
10. **It changes only what needs to be kept and carried over.** Where the paradigm window cannot fit the forced region, the forced region is
    not kept, unless its d-depth certification requirement is more than the lowest on the list. There are conditions, based on the
    architecture of the program, that dictate what the clean-up does when it is applied to particular paradigm window considerations.

### 1.10 The landscape, the pool of terminals, and the rule set (the designer's, 15 November 2026)

1. **The landscape of full coverage is the choice region and the force region.** The shape that the force region takes relative to the choice
   region, and the primitive composition calculations that apply, are directly linked.
2. **What is known of each.** In the choice region the primitive composition is such that a terminal is never reachable. In the force region
   there are reachable, calculable primitives. The positions of the forced regions for which the primitive calculations have not yet been
   composed, relative to the choice region and to other force-region board states, are also known.
3. **The structure that this gives.** It gives a way to understand how primitive calculations of an explicit d-depth scale relative to the
   landscape. The designer's interest in the composition of the landscape from choice to force, in relation to d-depth, is because of the link
   between primitive composition and d-depth, and the structural requirements of the primitive composition that follow.
4. **A bounded pool of terminals.** If a terminal condition is reachable within d-depth, a limited pool of potential terminal conditions exists
   from the board state. This bounds the scope of the viable choices, compared with the pool of choices there would be without the guidance
   of the program. (§1.11, item 1, distinguishes the pools.)
5. **Why d is kept.** d is kept so that the primitive composition can be recalculated with the program's architectural considerations. The
   scope of resolution went from an exhaustive search until the forced mate is found, to knowing the d-depth and the mating terminal
   requirements, under optimality, that must exist because of the d-depth calculation.
6. **The legend makes the recalculation easier.** With the legend denoting the type, the recalculation is told the mating condition and does
   not have to consider the whole pool of options that could resolve the composition.
7. **The seed batches inform the rule set.** The two batches of a segment, the force region with its d and the choice region, inform the
   primitive compositional process and the creation of a rule set. Considering d-depth for the force region is what lets the rule set
   compress properly.

### 1.11 The pools, and how d informs reachability (the designer's, 15 November 2026)

1. **Three pools, which are different.** (a) The pool of legal moves. (b) The pool of terminal mating conditions that are possible within an
   irreversible segmented boundary. (c) The choice of moves that fit the primitive composition understanding of d-depth, which are the moves
   that allow the composition to exist in the first place. Where a position is already certified and only the optimality is being recalculated,
   the symbol denoting the terminal condition informs the choice directly.
2. **What d tells the calculation.** d directly informs how many plies must be looked ahead to find the mating condition, or otherwise the
   boundary board state condition. Within that many plies there is only a set number of ways to reach the mating conditions, and these directly
   inform the reachability of moves.
3. **Segmented dependency units.** Each irreversible segmented boundary, from the initial game position to the end of the game, is in itself a
   segmented dependency unit. Within it are composed more segmented dependency chains of moves, which navigate the segmented dependency units.
   The designer doubted that this had been formalized and said it might be a gap in the documentation. It is since stated as a structural fact
   in `boundary_board_state_primitive_catalog.md` (Status, "Segmented dependency units, and the segmented dependency within them"), which does not
   depend on this approach.
4. **Moves as measurements.** The number of moves required to reach the condition informs the viable placements of the terminal mating condition.
   Each move can then be a measurement towards a pool of viable mating conditions. There exists a composition of moves that forces the outcome
   and follows the optimality constraints, choosing the move most adversarial against it.
5. **Each ply informs the next.** As the calculation progresses through each branch, each ply informs the calculation requirements further from
   each board state moment. The viable pool of terminal mating positions, given the response at each ply, must fit within the d-depth
   requirement, and that viability changes after each moment from the anchoring point.
6. **Termination.** Progress can be terminated at the point in a branch where there is certainty that the terminal cannot be reached optimally.
7. **Optimizations in the workflow.** Caching from previous work: landing on a resolved board state with known downstream solutions stops the
   recalculation. A repeated board state is not optimal, because the loop can be cut out of the branch.
8. **The pool along an optimal line.** Every d-optimal move leaves at least one option in the pool of terminal mating conditions, until the
   requirement is met.
9. **The budget.** A branch is not the correct one when no terminal condition is reachable within the required d-depth ply count, given the
   number of moves made to reach the board state: when d-depth minus the number of moves made is below zero, or is zero while the board state
   is not the terminal. A move that is not optimal excludes the branch from ever meeting the d-depth requirement.
10. **The conditions are calculable, for both sides.** The pool of viable terminal outcomes has conditions for reaching them that can be
    calculated. At any moment there is an understanding of this calculation for both sides: the calculations in their favour and against their
    favour. The calculation progresses relative to the perturbation of those conditions by any legal move.
11. **A bounded calculation.** This bounds the work as a quick calculation for a bounded solution.
12. **Sub-optimal moves.** A sub-optimal move can help fulfil the requirement of the condition, but not when optimality is considered in choosing
    from a pool of moves. The optimal choice is the move that calculably keeps the side the least far from the terminal condition, and
    closest to the terminal condition that favours it.
13. **What makes a move sub-optimal.** The primitive calculations advance one side in a way that cannot be escaped, and that excludes the
    sub-optimal move, especially when another option with less cost is already calculated to be provably better.
14. **The zero case.** When d-depth equals the number of moves made, so that the difference is zero, the line is one of the optimal routes (there
    can be many), if it is proven to have been chosen from optimal move decisions.
15. **Coverage remains exhaustive.** This still requires covering exhaustively. The exclusion that follows from the dependence on d-depth helps to
    alleviate depth concerns.

### 1.12 What optimality within a segment does and does not cover (the designer's, 15 November 2026)

1. Optimality found within an irreversible segmented boundary, within a certain depth, does not address how optimal it is to arrive at the
   lower dependency if that dependency is not resolved.
2. d can be set deep and never find a forcing move.
3. Without a full grasp of the terminal mating conditions that inform the calculations as they continue down, the accuracy of optimality is
   affected, and this absolutely has an impact.
4. This is what the clean-up is for. As the program progresses, gaps in the certification at certain d-depths resolve with more certainty,
   within the chain of segmented dependency units, and that resolves a branch for an outcome.

### 1.13 Building from the edge of the certified boundary (the designer's, 15 November 2026)

1. A comparison of exits whose lower values are unresolved can be left pending exhaustive resolution. It is in the pool of uncertainty.
2. It would be wise to have the program built up from the edge of the Syzygy certainty boundary itself, after enumerating the results for the full
   Syzygy landscape.
3. This scopes the clean-up pass more easily.

### 1.14 The 50-move rule as a calibration, the one-way ordering, and fetching against the frontier (the designer's, 15 November 2026)

1. The 50-move rule dependency can be a foundational rule for reasoning, to create several fronts of the program that have global rules. An
   example: the global sum of the segmented dependency unit move chain cannot exceed 50, otherwise it is a draw.
2. The main program cannot have this global rule, for finding true optimality. The 50-move rule is good for optimizing for normal chess.
3. Importing certificates is a concern that can be addressed depending on whether it is chosen to matter in the first place.
4. The ordering runs one way. This is a deliberate architectural choice, made from the beginning of the definition of the use of irreversible
   segmented boundaries as a formality.
5. The reason is that irreversible moves directly affect the reachable terminal mating conditions in the irreversible segmented boundary region
   and in the other reachable regions downstream from it. The further the progression, the narrower the scope.
6. **The rule the designer intends for a front.** The designer took the 50-move rule to be a total over the game, and states that this is what
   was meant: the game cannot exceed 50 moves globally from the initial start position. A requirement follows that the backward connection to
   the initial position respects the 50-move limit.
7. **The designer's correction of items 1 and 6.** The articulation of items 1 and 6 was based on a misunderstanding of the rule itself. The
   rule is the sum of the move chain within an irreversible segmented boundary, and not the sum of the units from the initial position. The
   designer then sees that the rule is based on what was defined as the irreversible segmented boundary.
8. **The rule as a calibration.** The 50-move rule remains a thing to consider against true optimality. It may be a good way to avoid having to
   record cursed wins, such as a very long manoeuvre to win a piece before a mate is forced. The program can allow for no limit, a 50-move
   limit, a 75-move limit and so on. This is similar to the d-depth calibration, as another variable that can be calibrated to prune for
   storage. Lowering it prunes the lines that need an extreme depth to resolve.
9. **Whether it is the same as a d-depth.** The designer asks whether the rule is not the same as setting a d-depth of 50 or 75, whether this
   matters, and whether the program can include more and opt into the rule through the d-depth calibration chosen at the start.
10. **Fetching and the frontier.** There is a difference between fetching data up to a d-depth certainty that fits within the storage limit of
    one's machine, and working at the frontier of the program's progression.

### 1.15 The primitives, the two modes, d-depth scoping, the legend, and what is deferred (the designer's, 15 November 2026)

1. **Fetching and the frontier should be formal modes.** The designer states that they should be formal modes of the program.
2. **The frontier will exceed most machines.** The designer is sure that, at a certain point, if the program is taken seriously, the
   frontier will exceed the capacity of most people's machines.
3. **Learning as well as progression.** The designer states that the program can be very useful for chess study, and that it is essential to
   allow for learning as well as for the progression of the program, for scientific discovery.
4. **The primitives are being specified, and the segmented dependency primitive is specified.** There are three primitives: the mating or
   capture net, the coverage net, and the segmented dependency primitive. The designer states that the primitives are being fully specified,
   and objects to calling the segmented dependency primitive unspecified. The terminal condition of the mating or capture net is specified,
   and the segmented dependency is understood through the compositional differentiation between the chains of moves and the irreversible
   segmented boundary units (the catalog's statement of 15 November 2026).
5. **Integrated specifications follow the progression of the program.** There are integrated specifications that are not yet specified,
   and that depend on the progression of the Lakatosian programme, as is happening in this work. The designer states that this is fine.
6. **The two modes.** Frontier mode is the full fetch of given irreversible segmented boundaries. Fetch mode is the calibrated d-depth
   fetch of an irreversible segmented boundary. The user then chooses which segments to fetch.
7. **The legend.** If the legend can be specified in the documentation, it should be, now, while following the architectural constraints of
   the program. (The specification is the designer's statement of 16 November 2026 in `boundary_board_states.md`; this file's reading is §2, item 31.)
8. **d-depth scoping is a structural requirement of the programme.** The designer states that it is integrated into the fetching mechanism, is as
   structurally required as the primitives, is invariant to the programme and is not isolated to an approach variant. The statement is in
   `boundary_board_states.md` and in the catalog (§4.10).
9. **The clean-up and the interdependencies are deferred.** Understanding how the interdependencies resolve, and the clean-up, is too large
   a task to undertake at this time. The designer leaves it to be built later, as an integrated specification.
10. **The study use is the calibrated fetch (the designer's statement, 16 November 2026).** Use for study is the same event as the calibrated d-depth fetch of item 6. The d is the metric threshold for the depth capacity under which the student wants to study (`boundary_board_states_d_depth_applications.md` §1, item 8).

---

---

## 2. This file's readings, offered for correction

1. **Two kinds of certificate.** A win is a finite proof. A draw is a closed set: at every position in it the mover has a move that stays in
   the set, and every move of the mover that leaves the set goes to a position certified lost for the mover. Wins are therefore built up from
   the forced terminals, and draws are what remains when the set stabilises. This agrees with the statement that wins gain certainty and
   draws require exhaustive certification (§1.3, item 1).
2. **A set is a certificate only when it has stabilised.** A position found decisive that other members relied on as their staying move
   removes those members unless they have another, and the removal can cascade. Until it stabilises, a partial certificate is a conditional
   statement with a dependency record (start-selection §4, §2f), and the unresolved dependent set measures the progress.
3. **The membership test needs a completeness condition.** "Not in d" means decisive only if every position of the scoped reachable region
   was enumerated and classified. A list of the uncertain positions found so far cannot certify anything, because absence from it may mean
   unvisited. Likewise "not in c means drawn" (§1.8, item 4) holds once c is certified complete for the scoped region, which is what the
   stabilised remainder establishes. Neither reading is shown here.
4. **The shell.** The part of c worth importing into e is the verdicts of the positions one move outside the branch of d being worked on,
   since only these can decide whether a position of d resolves. They are also what the closure check of item 1 consumes. The room of §1.7,
   item 6, is then a budget for the shell, and the final feasibility condition would be that the draw closure plus its shell verdicts fit the
   window. How large the shell is relative to the set is not known.
5. **Time against storage.** A verdict for a pruned position can be re-derived by a search to its d from the seed or the line that certified
   it. Pruning trades storage for time, and the re-derivation is a certification re-run, not a new claim.
6. **Merging.** A position certified by any batch is certified overall, so the uncertain list after a merge is the intersection of the
   uncertain lists, and the certificates are the union. This is the closure reading of start-selection §2f applied to batches.
7. **The anchor as the interface.** Certificates are keyed by position, so a seed that transposes into any certified position inherits its
   verdict. In this reading the anchor is distinguished as the point that carries the other half of the status, that a position is reached
   under optimality from an anchoring point (main §14 items 7 to 10). The seed path joins the certified region at the anchor, which is why
   independent seeds need no coordination. This distinction is put to the designer in §4.
8. **A depth that suffices is a property of a segment.** That a given depth suffices to find a forced sequence if one exists is known for a
   segment from the exhaustive coverage that settled it. Where there is no exhaustive coverage the depth is not known beforehand, and a
   position with no forced result within d is a lower bound on its status and not the status. The boundary board state status also needs the
   position to be reached under optimality and the draw to be realized (main §14 items 7, 8 and 10).
9. **Engines propose.** An engine move or a played game is an assertion until certified (start-selection §2e, readings). If each engine picks
   the deviation that relieves the most work, the proposals can converge, which works against the diversity of §1.4, item 2. Which
   deviation a contributor works on is the community's work selection, and the choice among deviations is not specified.
10. **The frontier.** The boundary of the uncertain region, and the points at which scoped attempts reach the window, would let a contributor
    choose seeds without any seed being supplied. Whether the program publishes it is this file's suggestion and not the designer's.
11. **Re-certification runs in one direction.** An irreversible boundary scopes each segment, so a result in a lower segment does not depend on
    a result above it. A change to a lower certificate, whether a conditional entry confirmed or overturned, bears on the results that cite
    it through the dependency record (start-selection §2g). The scope of the clean-up of §1.9, item 6, would then be the results whose terminal
    flag marks an entry into the changed segment, and the draw certificates that rest on such an entry. A draw that rests on an entry loses its
    closure if the entry's status is overturned. A result whose terminal ends the game inside its own segment depends on nothing below.
12. **What a legend declares (read with the designer's statement of 16 November 2026).** For each type of terminal: whether the terminal ends the game inside the segment or is an entry into a lower segment;
    the segment entered; and the status that the entry carries (certified, conditional, or likely forced and not yet certified). An entry that
    returns the position to a boundary board state classification is an edge of the choice region and is recorded there, not as a terminal of
    the forced region. The table of terminal states is a file of its own, the README of the segmented boundary describes it and declares the
    types (§2, item 31), and a general document states only that it exists and what it declares.
13. **"Likely forced" is a progress point and never a certificate.** The symbol records where the dependency sits and does not stand in for the
    status. A heuristic orders work and does not certify (the standing rules at the head of `boundary_board_states.md`).
14. **The flag as a set (superseded).** This reading expected that the terminal flag might need to name a set of types. The designer's statement
    of 16 November 2026 records one predominant branch, with the alternates reachable as anchor points, and the index names a single terminal
    state. No set is recorded.
15. **A clean-up pass on contribution.** Read as a pass that each contribution triggers: find the results that cite the changed certificates,
    through the dependency record; re-certify them or demote them to conditional; then apply the keep rule of §1.9, item 10, to the window.
    Each step is a re-check against the record and not a judgment, which is what would make it automatable. A certified result changes only
    through a re-certification that is itself recorded. The pass is not specified.
16. **The keep rule as an ordering.** "The list" is read as the forced-region entries ordered by the d that their certification needs. When
    the window cannot hold the forced region, entries of the lowest d are dropped first, and an entry is kept only if its d is above the
    lowest on the list. This agrees with §1.8, item 2 (a depth of 15 kept over a depth of 5). A dropped verdict can be re-derived by a search
    to its d (item 5).
17. **A pruned entry must be told from an unvisited one.** The membership reading of item 3 needs the record to show that the region was
    covered and that the missing forced entries were dropped by the rule. Without that, a forced entry absent from the store cannot be told
    from a position never visited, and "not in c" would no longer mean drawn. So the keep rule would need a record of what was pruned and by
    which ordering and threshold, and not only of what remains.
18. **The bound on the pool of terminal conditions is a necessary condition.** This is pool (b) of §1.11, item 1. A terminal reachable within
    d plies must be of the stated type and reachable by the pieces in that many plies. If the favoured side is to move it has (d+1)/2 moves within d, and if the other side is to move the favoured side has d/2.
    This removes candidate terminals and certifies nothing. This file reads the bound of §1.9, item 8, as applying to the stated terminal and not to
    the completed mate, in line with §1.10, items 4 to 6.
19. **A label where coverage is exhaustive, a scope under test where it is not.** Where exhaustive coverage exists, d and the terminal type of a
    position are labels. Where it does not, d is the scope of an attempt (can a composition force a terminal of this type within d?), and the
    pool bound holds for that attempt. The result stays conditional until the coverage of the segment settles it (§2, item 8).
20. **A composition proposes until it certifies.** A rule set tested against the two batches meets two kinds of failure. A prediction of force
    for a position of the choice region refutes the composition, and a force position for which no composition finds the terminal within its
    d shows a gap in the primitives. Both are findings for the protective belt. A composition certifies only where the primitive calculation
    is itself a proof (the standing rules at the head of `boundary_board_states.md`).
21. **The depth at which coverage resolves.** Where exhaustive coverage exists, the largest d found is the calibration for the depth bound of
    the exclusion form of the rule set (start-selection §2c and §2g, and §1.2, item 1 above). It is a property of the segment and is recorded in
    the segment's own record, not here.

22. **Exclusion only on certainty.** In this reading, §1.11 items 6 and 13 give the grounds for removing a move from the work: a proof that the
    terminal cannot be reached optimally in the plies that remain, or another option already proven to cost less. Short of a proof, a move is not excluded. It is covered within the stated d. A filter that removes moves by their kind, without a proof of unreachability, is a heuristic exclusion
    and is not part of the direction. A test of the direction would derive the compositions from the labels of the two batches, level by level.
    A move that the composition does not account for is an anomaly that revises the composition. It is not an error rate, because the composed
    calculation is exact (catalog §9).
23. **Two grounds for exclusion, both proofs.** In this reading, from §1.11 items 9, 12 and 13, a move is removed either because the terminals
    that favour the side are provably out of reach in the plies that remain, or because another option is already proven to cost less. The
    second needs the cost of the better option to be certified, not estimated. The moves that remain are compared adversarially, each side taking
    the move that is best for it against the other's reply. This follows the two readings of catalog §4.7 (the cooperative one that measures
    reachability, the adversarial one that does the work) and the statement in start-selection §2g that a found forcing sequence beats the
    branches it dominates, so that they need not be enumerated. Coverage stays exhaustive (§1.11, item 15), and the exclusions reduce the depth
    that has to be explored and not the region that has to be covered.
24. **The zero case, and routes.** Item 14 confirms this file's reading of the budget formula: at zero plies remaining, the line is on an optimal
    route if it is at the terminal, and it is not the correct branch if it is not. There can be several optimal routes, so the exclusion by cost
    leaves every route of least cost and removes only the dearer ones.
25. **Optimality upstream is relative to the values of the exits.** In this reading, what is certified inside a segment is the choice among
    moves given the values that the exits carry. A comparison between two exits is certified only if the lower values it uses are certified.
    Otherwise it is conditional, and the record would say which lower values were used and their status (§1.9, item 7; start-selection §2g).
    Wins and draws differ on clean-up. A forced win found below is permanent, so the upstream results that rest on it stand. A draw certified
    to a depth below can be overturned by a later forced win, and that can change the upstream results that rested on it. A clean-up therefore
    has to re-check the dependents and cannot only add.
26. **The frontier as the scope of the clean-up.** In this reading, from §1.13, a layer built from the certified edge takes the lower values that
    its exits use from the enumerated landscape, so those values are certified. Conditional statuses then arise only for results that rest on a
    layer that is itself only partly certified, and a clean-up has to re-check the layers above the one that changed and not the whole chain.
    The one-way ordering of the dependencies, within a layer and between layers, is the designer's deliberate choice (§1.14, items 4 and 5).
27. **The Laws' rule is a limit on the chain inside a unit.** The Laws of Chess (FIDE Article 9.3 and 9.6.2) count the moves by each player since the
    last capture or pawn move, which are irreversible moves, and the count starts again at each of them. A draw may be claimed once 50 moves
    by each player have been completed without either, and the game is drawn without a claim at 75. In the program's terms that is a limit on
    the move chain inside one segmented dependency unit, restarted by each crossing of an irreversible boundary, as the designer's correction
    (§1.14, item 7) says. The game as a whole is not bounded by it. A front that adopts the rule certifies under it, so its certificates are
    certificates under that front's rule and not under the main program's, which has no such rule (`boundary_board_states.md` §7).
28. **The rule against a d cap: the unit, and the result.** In this file d counts plies (§1.11, item 9; a d of 5 leaves the winning side at most 3
    moves). The rule's 50 and 75 are moves by each player, which is 100 and 150 plies, so the matching cap is a d of 100 or 150 and not 50 or
    75. The two also differ in what they say of a line beyond the limit. With a d cap and no rule, a win that needs more plies than the cap is
    unknown, since nothing certifies it. With the rule, the same line is a certified draw, because the side that is losing can hold by running
    out the count. A d cap therefore prunes the work and leaves the status open, and the rule changes the status. Both keep the very long
    lines out of the store.
29. **Two ways of using the program (superseded in part by §1.15, item 6).** An earlier reading of §1.14, item 10, took fetching to be the
    use of certified data up to the d that the window holds, and the frontier to be work that extends certification where none exists. The
    designer's statement of §1.15, item 6, defines the frontier mode differently: as the full fetch of given segments, with fetch mode as the
    calibrated d-depth fetch. This file takes the later statement as the definition. Whether "working at the frontier of the program's
    progression" (§1.14, item 10) and the full fetch of a segment are the same thing is for the designer to say.
30. **What the two modes would say.** In this reading, both modes are selected by the user per segment, and differ in what is requested of
    that segment. Fetch mode takes a segment and a d, and returns the results certified within that d, each with its d and terminal flag. A
    position with no result within the d is unknown, and is not reported as a draw (item 28). Frontier mode takes a segment and returns
    everything the program holds for it, which is the largest quantity in the program and the one that §1.15, item 2, expects to exceed most
    machines. Both would state the status that each result carries (certified, conditional or unresolved), the rule front if one applies
    (§1.14), and the legend of the segment. What "everything the program holds" includes, such as the pending and unresolved lists, how
    large a fetch is estimated before it is made, and how a user's machine is matched to a request, are not specified. The study use of
    §1.15, item 3, is read as the same data used for another purpose: d, the favoured side and the terminal flag say how many plies a
    certified result needs and what it ends in. The study use is the calibrated fetch (§1.15, item 10), and the program records the learning use and the discovery use as equally required.
31. **The header, the legend file and what the README declares (this file's reading of the designer's statement of 16 November 2026).** The
    statement is in `boundary_board_states.md` ("The d-depth metric"). This reading replaces the draft of 15 November, which put a terminal
    symbol in the header and the legend in the README.
    1. **The header** is [d][index][symmetry condition]. The index points into the legend file, and the entry it points to gives the condition
       type (mate, capture, or pawn push). No separate symbol and no set of terminals is recorded.
    2. **The legend file** holds the terminal board states, once per symmetry class, and only those. Colour is not stored. The group of
       transformations in force is the one that keeps a position inside the segment.
    3. **The README of the segmented boundary** describes the legend file and the transformation symbols, and declares for each type: whether the
       terminal ends the game inside the segment or enters a lower segment; the segment entered; the status of the entry (certified to a stated
       scope and by a stated source, conditional, or likely forced and not yet certified); and what d counts to. A README declares and does not
       certify. "Likely forced" is a progress point and never a certificate (item 13).
    4. **An entry that returns the position to a boundary board state classification** is an edge of the choice region and is recorded there,
       and it is not a terminal of the forced region (item 12).
    5. **Which terminal is recorded when several are reachable** is the designer's predominant-branch statement. The rule used to pick it, such
       as the lowest index, belongs to the application of a segment and not to a general document.
    6. **Not specified here.** How the status entries of a README are updated when a lower segment changes belongs to the clean-up, which the
       designer defers (§1.15, item 9). The pawn condition is stated by the designer in `boundary_board_states.md` and is not exercised in a
       segment without a pawn.

32. **The stored band and what absence means (this file's reading of the designer's statement of 16 November 2026 on calibrating storage against
    search time).** The reading is held in `boundary_board_states_calibration.md` (§2, items 1 to 3), together with the designer's statements it reads. The keep rule of §1.9, item 10 and §2, item 16, is the lower threshold in that file's terms: entries of the lowest d are dropped and re-derived by search.

## 3. Relation to the existing record, and one tension

- **Agrees with** start-selection §2g (depth-scoped certification, the depth ending at the irreversible boundary, the dependency record
  across segments), §2d (the forward anchor, and games as seeds that are witnesses for conditional links), §2e (the fronts, and that a
  certificate can be checked by a contributor other than the one who produced it) and §2f (batches, the closure reading, the interface
  values at a crossing).
- **A possible tension, or a different approach.** The designer indicated that the start-selection sections may describe a different approach (and that part of them may detail pruning the forced region to create this dichotomy). If so, the statements below are two options and not a conflict. Start-selection §2f records the designer's statement that the object holds the seeds that inform the reachability family
  under consideration, and also reads the seed paths as provenance. The statement of §1.4, item 2, is that the program repository does not
  keep seeds, because a supplied seed basis would bound discovery. These may be compatible if the object keeps witnesses as provenance while
  the program supplies no canonical seed set, but that is this file's guess and the two statements are not reconciled. Whether they are two approaches is for the designer to say. §2f is not changed by
  this file.
- **What this file adds to the record** is the placement of batches in the window (§1.7), the c and d stores and the flag (§1.8), the
  constraints on a repository of games (§1.5), the role of the anchor as the convergence point of seeds (§1.4), and the terminal reason of d,
  its legend and clean-up certification (§1.9), the landscape, the pool of terminals and the rule set (§1.10), the three pools (§1.11), what optimality within a segment covers (§1.12), and the two modes and the draft of the legend (§1.15; §2, items 30 and 31). The statement of §1.12 refines start-selection §2g (the depth ends at the irreversible boundary; dependency
  across segments) for this approach, and does not change it.

## 4. What is open

- **The anchor and interior positions.** Whether the distinction of §2, item 7, is what the designer means by "the convergence point must be
  the anchor point".
- **The two stores.** Whether the program is to hold c as a certified record and d as the working list, or the reverse; and the format of
  each, including what is kept on another machine and what is imported (§1.8).
- **The shell.** Its size against the window is unmeasured. The designer's calibration principle would measure it on a material where the
  tables hold the answer, before any use where they do not.
- **The impact of a seed.** Its marginal contribution depends on the state of the store when it is run, so an impact measure is relative to a
  snapshot. Not specified.
- **The smart selection from c** that forms e, and the order in which branches of d are worked. Not specified beyond §2, item 4.
- **Calibrating an engine's search depth to the program's d.** Not specified.
- **The illustration's split.** The forced and choice fractions were stated as 70% choice and 80% forced, and later as 70% forced and 80%
  choice. §1.7 uses the later. The illustration is not a measurement and does not depend on the split.
- **Whether start-selection §2f and this file are two approaches or one** (§3).
- **The scaling of the legend** to other irreversible segmented boundaries with other terminals. The designer called it a different story. Not
  addressed. The legend and the header are specified by the designer's statement of 16 November 2026; §2, item 31, is this file's reading of what the README declares.
- **Detecting and owning clean-up certification.** How a changed certificate is noticed, who re-checks the results that cite it, and the format
  of the record are not specified (start-selection §4 and §5 say the record is not specified). The designer defers this, with the resolution of
  the interdependencies, to a later integrated specification (§1.15, item 9).
- **The specification of the clean-up pass.** The trigger, the format of the record, the ordering and its ties, and how a demoted result is
  marked. The designer says it requires integration and further work, and defers it (§1.15, item 9). Not specified.
- **What "the list" is in the keep rule** (§1.9, item 10), which this file reads as an ordering by the d that certification needs (§2, item 16).
- **The conditions that dictate the clean-up for each window case.** The designer states that the architecture dictates them. They are not
  enumerated, and none is written here beyond the keep rule.
- **The first test of compositions derived from the batches.** Level by level, which primitive compositions account for each d-optimal move of the
  labels, and how an anomaly revises a composition (§2, item 22). Not run and not specified here. A first attempt that restricted the moves by
  kind was a heuristic exclusion, and it is withdrawn.
- **The measure of how composition scales with d** relative to the landscape (§1.10, item 3). The designer's interest is stated, and no measure
  is specified.
- **How the enumerated landscape is held and imported** (§1.13): what is kept within the window, in what convention (which depends on the choice
  of §1.14, item 3), and with what status for a certificate that comes from the tables. Not specified.
- **The unit of the rule's calibration** (§1.14, item 8; §2, item 28): moves by each player, as the Laws state it, or plies as d does. Which
  fronts exist, and how a certificate carries the rule it was certified under, are not specified.
- **The cost model and its scaling beyond one segment**, and the other open items of calibration, are in `boundary_board_states_calibration.md` §3.
- **The two modes** (§1.15, items 1 and 6; §2, items 29 and 30). The designer states that fetching and the frontier are formal modes: frontier
  mode is the full fetch of given segments, and fetch mode is the calibrated d-depth fetch. The definition of each is not written beyond
  that. Open: what a full fetch includes (pending and unresolved lists), how the size of a request is estimated, what a fetch states of the status of its data (see `boundary_board_states_calibration.md` §2, item 6, on the stored band).
- **How a comparison of exits with unresolved lower values is recorded**, and how its conditional status is carried up the chain until the lower
  values are resolved (§1.12; §2, item 25). Not specified. It belongs with the resolution of the interdependencies that the designer defers
  (§1.15, item 9).
- **How the viable pool of terminal conditions is computed at each ply** (§1.11, items 4 and 5), and what counts as a proof of unreachability
  (§1.11, item 6). Not specified.
- **Whether the loop reaches the final state.** It stalls if the remaining forced positions need a depth beyond reach, or if the verdicts they
  need exceed the room. Nothing here shows which happens.

## 5. Status

| Item | State |
| --- | --- |
| The depth layer and the segment gap (§1.2) | Stated by the designer; the d-scoped certification is in start-selection §2g; no run recorded here |
| Seeds converging on the anchor (§1.4) | Stated by the designer; no run recorded here |
| The repository constraints and the engines (§1.5) | Stated by the designer; no run recorded here |
| Partial certificates and the window (§1.6) | Stated by the designer; no run recorded here |
| The batches c, d and e (§1.7) | Stated by the designer as an illustration; no run recorded here |
| The shell and the feasibility condition (§2, item 4) | This file's reading; no run recorded here |
| d to a stated terminal reason, the flag and the segmented dependency (§1.9, items 1 to 5) | Stated by the designer; no run recorded here; scaling to other boundaries not addressed |
| Clean-up certification and the legend (§1.9, items 6 and 7) | Stated by the designer; the legend's symbols and convention not specified |
| The landscape as choice plus force, the bounded pool of terminals, d kept for recomposition, the batches as input to the rule set (§1.10) | Stated by the designer; no run recorded here |
| The 50-move rule as a calibration, the one-way ordering, fetching against the frontier (§1.14) | Stated by the designer, with the designer's correction of items 1 and 6 (item 7); the fronts are not specified; no run recorded here |
| Building from the edge of the certified boundary after enumerating its landscape, and the scope of the clean-up (§1.13) | Stated by the designer; no run recorded here; the import of the landscape is not specified |
| What optimality within a segment covers, and what the clean-up resolves (§1.12) | Stated by the designer; the record of conditional exit comparisons is not specified; no run recorded here |
| The three pools, d informing reachability, the segmented dependency unit, termination and the optimizations (§1.11) | Stated by the designer; the unit and the chains inside it are stated as a structural fact in the catalog; no run recorded here |
| Automated clean-up on contribution and the keep rule (§1.9, items 9 and 10) | Stated by the designer; the specification and its integration are not done; no run recorded here |
| The two modes, the primitives as specified, the legend to be specified now, the clean-up and interdependencies deferred (§1.15) | Stated by the designer; the modes are named and not defined beyond that; the clean-up is deferred by the designer |
| Direction of re-certification, what a legend declares, the flag as a set, the clean-up pass, the keep rule as an ordering, pruned against unvisited, the modes, the reading of the header, the legend file and the README (§2, items 11 to 32) | This file's readings; no run recorded here |

## 6. Provenance

**The designer's:** §1 and the framing of the Status, in the order stated on 15 November 2026, including the illustration of §1.7, the statements on the two stores and the statements of §1.9, §1.10, §1.11, §1.12, §1.13, §1.14 and §1.15.
**This file's:** §2, §3, the list of §4 and the table of §5. No result of any run is recorded in this file.
