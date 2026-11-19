# Calibrating Storage against Search Time: the d-depth Metric as the Control

**Status.** A record of a direction the designer articulated on 16 November 2026. It states how the d-depth metric calibrates what a store holds against the search that rederives the rest, and it is the file where that calibration is developed. It is stated in general terms: no piece, material or position is part of it. It adds no primitive and does not touch the hard core (`boundary_board_states.md` §14 items 7 to 11). Statements marked **the designer's** are the designer's. Statements marked **this file's reading** are readings of the designer's words, offered for correction. No measurement is recorded in this file, and the designer states that the cost of the trade is not known.

Companion documents: `boundary_board_states.md` (the designer's statements at "Where this stands": d-depth scoping, the d-depth metric, calibrating storage against search time), `boundary_board_states_start_selection.md` (§2g depth-scoped certification and dependency across segments), `boundary_board_states_batch_composition.md` (the keep rule, §1.9 item 10; pruned against unvisited; the two modes, §1.15), `boundary_board_states_syzygy_mining.md` (the certified base that this builds on), `boundary_board_states_d_depth_applications.md` (the applications of partial d-depth certification).

---

## 1. The designer's statements (16 November 2026)

### 1.1 Exhaustive coverage and the program

1. **Small segments.** For a small segment, exhaustive coverage may be the better method, because it fits the paradigm window and delivers at speed. The designer states that this is most likely the case for most of the small piece combinations, although each also occupies paradigm window space cumulatively.
2. **Where the program takes over.** At some point in the opportunity cost calculation it can be calculated that it is better to leave some plies unresolved and to convert them to a search problem. As the coverage of pieces is added to pass the limits of the Syzygy tables, this can be calibrated to optimize the speed of delivering the optimal branch while fitting the most within the paradigm window.
3. **The trade.** This increases what the storage can reach, at the calculable cost of the time to solve a search problem whose difficulty increases, in a calculable way, with each ply removed.
4. **The edges of the Syzygy solution.** Storage can be removed from the edges of the Syzygy solution, so that more complicated regions become reachable.
5. **The overhead against the Syzygy schematic.** Each stored board state records only one of the most optimal routes, in the schematic of the program. That is an explicit opportunity cost, and the storage overhead of the header and its information, against the schematic of Syzygy, is a calculable amount.

### 1.2 What the d-depth metric does here

6. **The metric calibrates it directly.** The d-depth metric is the control. The designer states that d can be calibrated differently on each region, and that it goes both ways. It can be applied to both regions, the forced and the choice regions, and which it is applied to is a choice for the user. For chess study the user would likely want the drawn regions that fit within the d-depth that they can calculate.
7. **A chosen d is an admission.** A chosen d is also an admission that optimal play is out of reach, with certainty, within the delivery time that is wanted. This ceiling is reached whether or not the program provides the answer, and performing work while the study continues is not a barrier.
8. **Cursed wins.** The Syzygy solution space can be pruned to what d-depth reachability of certainty gives, so that perfect play under partial d-depth certification keeps the cursed wins, which blow up the Syzygy storage, from being stored. Instead of being exact within an infinite d, the program is calibrated to remove that blow-up and the positions that surpass the certainty that can be derived within the d-depth certification. For study, the same exclusion keeps the student out of rabbit holes that only occupy storage and distract from the depths of analysis that matter (`boundary_board_states_d_depth_applications.md` §1, item 8).

### 1.3 What is stored, and what is rederived

9. **The chain between the anchor and the terminal is not indexed.** The d-depth gives the furthest anchoring point and the terminal condition. The board states between them in the segmented dependency move chain do not need to be indexed, because the resolution of the move chain in the primitive composition resolves the board states within the composition.
10. **Symmetry and folds.** The storage can be compressed with the folds, as Syzygy does. Symmetry conditions can also be used to resolve alternative branches.
11. **Alternative equally optimal routes (open, the designer's [A]).** When a branch is resolved, the points where equally optimal alternatives branch off inform reachability. The designer states that this is an open front to leverage for certainty, that the difference between the anchor and an alternative equally optimal forcing sequence can be calculated, and that the designer will return to it. The integrated specifications are to be found by working with it.
12. **Only what is reached.** The program does not need to enumerate the whole drawn region. Only the positions that are reached in study need to be enumerated. The program does not require exhaustive coverage, and can accomplish it or, to manage the paradigm window, not accomplish it, and delete from storage the branches that do not matter.
13. **A player's own branches.** Specific branches that recur in a player's own play can be enumerated. This informs the starting anchor and reachability, and it is progress that accumulates for the program.

### 1.4 Chess study

14. **Study benefits most.** A person studying chess does not need to know of the cursed wins that burden the Syzygy storage and limit what it can deliver beyond its piece limits. The d-depth metric is directly related to the amount of play the student is trying to master. At a d of 40 or more, a human has to consider optimality within a move chain of that size. The designer states that chess study benefits from the d-depth architecture the most, more than the scientific endeavour does. The designer refines the limit for study on 19 November 2026: the d-depth limit is not a hard-coded limit but a quantity derived from the proof-size limitation, which calculates it, and it depends on the proof size within d (`boundary_board_states_engine_integration.md` §1, item 5).
15. **Learning and progression.** It is essential to allow for learning as well as for the progression of the program (`boundary_board_states_batch_composition.md` §1.15, item 3).

### 1.5 Syzygy, compression and confidence

16. **Compression.** The Syzygy compression is strong. It is to be optimized later, as the program is fortified, building on the work that Syzygy has done, to expand the limits that the storage problem sets.
17. **Convergence.** The approach converges well with the Syzygy solution. The designer is confident that it will work.
18. **What is not known.** The cost of the trade is not known. What has to be found is how much the search time increases for each ply removed, how much the storage limitation is mitigated by the d-depth pruning, and how the segmented dependency units are best chained and calibrated.

### 1.6 Calibration patterns, the choice region and drawn states (the designer's, 16 November 2026, later in the session)

19. **Shallow forced results can be omitted from storage.** Positions with a forced mate in 1 ply or in 3 plies can be omitted from storage where they are numerous, because a search finds them faster than resolving the calculations would, as an ordinary chess engine already does. This takes a class of solutions out of the storage space.
20. **The d-depth puts optimality within reach.** The d-depth places the depth that is required into an obtainable window for practical play, and for understanding optimality within human reach. It can be calibrated for each irreversible segmented boundary.
21. **The moves between the links of the move chain are a calibration metric.** The positions kept along a chain can reach as far as possible, or can leave a calculated number of gaps, so that the calculations are more reachable to resolve in the time. For a required d-depth of 7, with o for kept and x for discarded: oxxxxxo, oxoxoxo and oxxoxxo are all possible.
22. **The best pattern differs between segments.** The difficulty of finding the optimal moves at the end of the chain is easier in some irreversible segmented boundaries than in others, so the emergent pattern can look more like ooooxxx or oooxxxx.
23. **The pattern is a rule over the values of d.** The pattern of o and x applies to the values of d, and not to one particular chain. This holds where d is calibrated to the d that exhaustive coverage requires, and below it.
24. **A cap below the d required.** If d is calibrated to 6 where 7 is required to resolve, then 6 is the maximum, and the gaps that follow are emergent from the 7 required.
25. **The scope of a gap.** The primitive compositional problem across a gap such as the one in oxxo is 2 plies, and the depth to reach certainty of optimal play is scoped between those two plies. That scopes the difficulty explicitly, and mitigates, in a calculable way, the storage needed for the indices. The longest gap bounds the search and the difficulty of resolving it, especially where the anchoring points on either side of the gap are already certified as forced optimally, so that the work is measuring optimality between two bounded points.
26. **Engines as part of a hybrid.** Chess engines can be used to resolve board states of less difficulty. This is a heuristic approach, guided, and it does not ensure optimality, but it gives a good starting place. It makes a hybrid possible in which branches are worked and resolved in parallel. The designer left this out as integrated specification, and states it here explicitly.
27. **The choice region needs the exhaustive threshold.** Resolving the choice region requires the exhaustive coverage threshold, so the force region is always available before the choice is resolved. This is a structural certainty.
28. **Drawn states, an open front.** The program will work on approaches to certify the drawn states, and to leverage the architecture to see whether emergent relations can be constructed between the property of a drawn equilibrium and the certainty that can be derived. If that is not possible, it is an inherent barrier to the program, and truly certifying drawn states exhaustively requires an expansion of the paradigm window. Short of that, the program resorts to d-depth calibration, which certifies partially within the d-depth scope.
29. **The choice region can be recorded from the same search.** A position that is not stored, and whose search of the length of the longest gap fails, can be recorded as a draw in the choice region, so that it is not left to guesswork when storage permits recording it. It is a draw certified within the d-depth scope, not a certified draw, and the designer states that this is still very useful.
30. **Shallow capture exits are trivial.** A capture exit that lies one to three plies away under force is an easy problem to reach, because the terminal condition is in the legend. It is the same as the xxx of ooooxxx, where the discarded rows are trivial.
31. **A branch is abandoned on proof.** Not every forced line has to be exhausted. A branch is excluded when it calculably cannot reach a terminal condition with advantage from where it stands on the move chain, given the d-depth requirement. This abandons a branch with confidence, without the full verification that exhaustive verification would require, and the designer states that it influences speed.

---

## 2. This file's readings, offered for correction

Items 9 to 14 are confirmed by the designer as the content that was left implicit. They remain this file's readings until the designer specifies them.

1. **The stored band.** Two different truncations are possible. A lower threshold k means that rows below it are rederived by search, and a cap K means that wins beyond it are not stored. A row below k is recoverable by a search of bounded depth. A win beyond K is not certified, and it is not a draw. A contiguous band is the special case of a pattern over d (item 9).
2. **What absence means.** With both, a position that is outside the store and whose search to k-1 fails is a conditional draw: no forced result within K. A position outside the store whose search to k-1 succeeds is a certified win of d at most k-1, and a win is permanent. This is the reading of item 8 and item 12 of §1 together, and it keeps the program's rule that wins are permanent and draws are conditional.
3. **Deletion.** A deleted row is lossless only if the search can rederive it within the budget. A pruned row has to be distinguishable from an unvisited one (`boundary_board_states_batch_composition.md` §1.9, item 10), or a certified result would be indistinguishable from an unresolved one.
4. **One space, two limits.** The exhaustive store is the limit with nothing below k and nothing above K. The pure search problem is the limit with nothing stored. A hybrid is anywhere between, set per region by d (§1, item 6), on the forced region and on the choice region alike, as the user chooses.
5. **What the cost model would hold.** Two quantities per value of d and per region: the bits stored for a kept row (its header and its share of the legend of the segment), and the search needed for a position that is not stored. The first is calculable from the header (§1, item 5). The second is not known beyond one segment (§1, item 18). Whether certifying a whole region by search, and answering one query by search, are the same cost is not assumed.
6. **The modes.** Fetch mode is the calibrated d-depth fetch, and frontier mode is the full fetch of given segments (`boundary_board_states_batch_composition.md` §1.15, item 6). Read here, what a fetch can state about a position depends on the band that the fetched store was built under, so a fetch would state that band. How it does so is not specified.
7. **The rule front.** A front that adopts a move-chain limit makes the wins beyond it certified draws under that front, and without the rule they stay unknown (`boundary_board_states_batch_composition.md` §2, item 28). The cap K and the rule front are different things, and a cap leaves the beyond unknown.
8. **Study.** The data of a stored band (d, the terminal reached) say how many plies a certified result needs and what it ends in. The choice of K sets the size of the chains that are certified, and the program records the learning use and the discovery use as equally required (§1, item 15).

9. **The stored band is a special case of a pattern over d.** Read with §1, items 21 to 24, what is stored is a set S of values of d. The band of item 1 is the contiguous set from k to K. Reading ooooxxx as a keep of the highest values and a drop of the lowest gives the lower threshold, and a cap of 6 where 7 is required gives K. A set with gaps inside it, such as oxoxoxo, is the general case. Which set a segment uses is left to measurement.
10. **A gap bounds the search.** A discarded row is re-derived by a search whose leaves are stored rows, the terminal states of the legend, and values from the lower segment. Reading §1, item 25, the search across a gap of g plies has a cost that depends on g and on the branching, and the longest gap in the set bounds the search needed for any position whose result lies within K, since such a result is at most one more ply than that gap away from a stored row or a terminal. The cost of a gap is not measured here.
11. **Absence with a pattern.** A position that is not stored is not a draw. If a search of the length that the longest gap sets, from that position, reaches no stored row and no terminal that decides it, the position is a conditional draw, certified within the d-depth scope: no forced result within K. The designer confirms this reading. It is item 2 restated for a set with gaps, and the lower threshold of item 2 is the case of one gap at the bottom. A draw recorded in the choice region (§1, item 29) is of this kind.
12. **Exits and the lower segment.** Omitting a shallow row (§1, item 19) whose result ends in a capture is lossless only if the lower segment's value of the position after the capture can be reached, stored or re-derived. This is the dependency across segments of start-selection §2g.
13. **The role of an engine.** A win that an engine finds is a certificate once its forced line has been checked against every reply of the other side, and a win is permanent. A position where an engine finds nothing proves nothing, so it stays a conditional result. This keeps the standing rule of the programme that heuristics scope the work and never certify (topology §1). An engine is a source of candidates for the anchors to keep, and the check that makes them part of the store is exact.
14. **The cap and the choice region.** Reading §1, item 27, with a cap K below the d that the segment requires, the choice region is not certified exhaustively, so its drawn status is a conditional one: no forced result within K (item 2). Whether the choice region can be certified exhaustively without the full force region is the open front of §1, item 28.
15. **The exclusion of §1, item 31 is the exclusion on proof of the catalog.** A branch is dropped when the plies that remain on the chain cannot reach a terminal with advantage. That is a lower bound on the plies needed (catalog §4.10 and §9: a move is excluded only on proof). Applied at each ply of the calculation, it is what makes exhaustive verification of every forced line unnecessary within the d-depth scope. How the bound is computed at a ply is not specified here.

---

## 3. What is open

- **The cost model.** Its units, per value of d and per region: bits stored, and search done per query. Nothing is measured in this file. The designer states that the trade is to be found by test.
- **The scaling of both costs beyond one segment.** Not known.
- **The setting of the set S of kept values of d** (§2, item 9), with k and K as its special case. Who or what sets it for a segment, for a contributor's machine, and for a mode. Not specified.
- **How a mode states the band** of the store it serves (§2, item 6). Not specified.
- **The format of the store.** The designer states that compression is to be optimized later (§1, item 16). Not specified here.
- **The alternative equally optimal routes** (§1, item 11). Open, and the designer will return to it.
- **The clean-up** that re-checks results when a lower segment changes. Deferred by the designer (`boundary_board_states_batch_composition.md` §1.15, item 9).
- **The rederivation by primitive composition** of the chain between an anchor and its terminal (§1, item 9). The composition has not been run.
- **The pattern that emerges for a segment** (§1, items 21 and 22). Which pattern a segment settles on. Not specified.
- **The cost of a gap against its length** (§2, item 10). Not measured.
- **How an engine is brought into the program in parallel** (§1, item 26; §2, item 13): what it proposes, how its candidates are checked, and how branches resolved in parallel are recorded. Not specified.
- **The drawn states** (§1, item 28). Whether emergent relations between the drawn equilibrium and certainty can be constructed. Open, and widening the paradigm window is the exhaustive route if they cannot.
- **How a recorded draw in the choice region is stored and marked** as scoped to d (§1, item 29). Not specified.
- **A request beyond the stored cap.** The student's limit, which depends on the proof size within d (`boundary_board_states_d_depth_applications.md` §1, item 8, as refined), and the cap K of the store (§2, item 1) are two limits. What a fetch gives when the requested limit lies beyond what the store was built to, and how a limit on proof size maps onto the cap K, are not specified.

## 4. Status

| Element | Status |
|---|---|
| The statements of §1 | Stated by the designer. The designer is confident it will work, and states that the cost of the trade is not known |
| The statements of §1.6 (items 19 to 31) | Stated by the designer. Item 28 is an open front |
| Items 1 to 8 of §2 | This file's readings, offered for correction |
| Items 9 to 15 of §2 | This file's readings; items 9 to 14 confirmed by the designer as left implicit |
| The cost model and its scaling, the cost of a gap, and the pattern that emerges for a segment | Open. Nothing measured in this file |
| The alternative equally optimal routes | Open, the designer's [A] |
| The rederivation by composition | Not run |

## 5. Provenance

The statements of §1 are the designer's, from the session of 16 November 2026, and the statement at "Where this stands" in `boundary_board_states.md` is the same direction in a shorter form. §2 is this file's reading. This file adds no result of any run.

The statements of §1.6 are the designer's, from later in the same session of 16 November 2026.
