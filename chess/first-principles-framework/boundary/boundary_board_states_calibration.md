# Calibrating Storage against Search Time: the d-depth Metric as the Control

**Status.** A record of a direction the designer articulated on 16 November 2026. It states how the d-depth metric calibrates what a store holds against the search that rederives the rest, and it is the file where that calibration is developed. It is stated in general terms: no piece, material or position is part of it. It adds no primitive and does not touch the hard core (`boundary_board_states.md` §14 items 7 to 11). Statements marked **the designer's** are the designer's. Statements marked **this file's reading** are readings of the designer's words, offered for correction. No measurement is recorded in this file, and the designer states that the cost of the trade is not known.

Companion documents: `boundary_board_states.md` (the designer's statements at "Where this stands": d-depth scoping, the d-depth metric, calibrating storage against search time), `boundary_board_states_start_selection.md` (§2g depth-scoped certification and dependency across segments), `boundary_board_states_batch_composition.md` (the keep rule, §1.9 item 10; pruned against unvisited; the two modes, §1.15), `boundary_board_states_syzygy_mining.md` (the certified base that this builds on).

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
8. **Cursed wins.** The Syzygy solution space can be pruned to what d-depth reachability of certainty gives, so that perfect play under partial d-depth certification keeps the cursed wins, which blow up the Syzygy storage, from being stored. Instead of being exact within an infinite d, the program is calibrated to remove that blow-up and the positions that surpass the certainty that can be derived within the d-depth certification.

### 1.3 What is stored, and what is rederived

9. **The chain between the anchor and the terminal is not indexed.** The d-depth gives the furthest anchoring point and the terminal condition. The board states between them in the segmented dependency move chain do not need to be indexed, because the resolution of the move chain in the primitive composition resolves the board states within the composition.
10. **Symmetry and folds.** The storage can be compressed with the folds, as Syzygy does. Symmetry conditions can also be used to resolve alternative branches.
11. **Alternative equally optimal routes (open, the designer's [A]).** When a branch is resolved, the points where equally optimal alternatives branch off inform reachability. The designer states that this is an open front to leverage for certainty, that the difference between the anchor and an alternative equally optimal forcing sequence can be calculated, and that the designer will return to it. The integrated specifications are to be found by working with it.
12. **Only what is reached.** The program does not need to enumerate the whole drawn region. Only the positions that are reached in study need to be enumerated. The program does not require exhaustive coverage, and can accomplish it or, to manage the paradigm window, not accomplish it, and delete from storage the branches that do not matter.
13. **A player's own branches.** Specific branches that recur in a player's own play can be enumerated. This informs the starting anchor and reachability, and it is progress that accumulates for the program.

### 1.4 Chess study

14. **Study benefits most.** A person studying chess does not need to know of the cursed wins that burden the Syzygy storage and limit what it can deliver beyond its piece limits. The d-depth metric is directly related to the amount of play the student is trying to master. At a d of 40 or more, a human has to consider optimality within a move chain of that size. The designer states that chess study benefits from the d-depth architecture the most, more than the scientific endeavour does.
15. **Learning and progression.** It is essential to allow for learning as well as for the progression of the program (`boundary_board_states_batch_composition.md` §1.15, item 3).

### 1.5 Syzygy, compression and confidence

16. **Compression.** The Syzygy compression is strong. It is to be optimized later, as the program is fortified, building on the work that Syzygy has done, to expand the limits that the storage problem sets.
17. **Convergence.** The approach converges well with the Syzygy solution. The designer is confident that it will work.
18. **What is not known.** The cost of the trade is not known. What has to be found is how much the search time increases for each ply removed, how much the storage limitation is mitigated by the d-depth pruning, and how the segmented dependency units are best chained and calibrated.

---

## 2. This file's readings, offered for correction

1. **The stored band.** Two different truncations are possible. A lower threshold k means that rows below it are rederived by search, and a cap K means that wins beyond it are not stored. A row below k is recoverable by a search of bounded depth. A win beyond K is not certified, and it is not a draw.
2. **What absence means.** With both, a position that is outside the store and whose search to k-1 fails is a conditional draw: no forced result within K. A position outside the store whose search to k-1 succeeds is a certified win of d at most k-1, and a win is permanent. This is the reading of item 8 and item 12 of §1 together, and it keeps the program's rule that wins are permanent and draws are conditional.
3. **Deletion.** A deleted row is lossless only if the search can rederive it within the budget. A pruned row has to be distinguishable from an unvisited one (`boundary_board_states_batch_composition.md` §1.9, item 10), or a certified result would be indistinguishable from an unresolved one.
4. **One space, two limits.** The exhaustive store is the limit with nothing below k and nothing above K. The pure search problem is the limit with nothing stored. A hybrid is anywhere between, set per region by d (§1, item 6), on the forced region and on the choice region alike, as the user chooses.
5. **What the cost model would hold.** Two quantities per value of d and per region: the bits stored for a kept row (its header and its share of the legend of the segment), and the search needed for a position that is not stored. The first is calculable from the header (§1, item 5). The second is not known beyond one segment (§1, item 18). Whether certifying a whole region by search, and answering one query by search, are the same cost is not assumed.
6. **The modes.** Fetch mode is the calibrated d-depth fetch, and frontier mode is the full fetch of given segments (`boundary_board_states_batch_composition.md` §1.15, item 6). Read here, what a fetch can state about a position depends on the band that the fetched store was built under, so a fetch would state that band. How it does so is not specified.
7. **The rule front.** A front that adopts a move-chain limit makes the wins beyond it certified draws under that front, and without the rule they stay unknown (`boundary_board_states_batch_composition.md` §2, item 28). The cap K and the rule front are different things, and a cap leaves the beyond unknown.
8. **Study.** The data of a stored band (d, the terminal reached) say how many plies a certified result needs and what it ends in. The choice of K sets the size of the chains that are certified, and the program records the learning use and the discovery use as equally required (§1, item 15).

---

## 3. What is open

- **The cost model.** Its units, per value of d and per region: bits stored, and search done per query. Nothing is measured in this file. The designer states that the trade is to be found by test.
- **The scaling of both costs beyond one segment.** Not known.
- **The setting of k and K.** Who or what sets them for a segment, for a contributor's machine, and for a mode. Not specified.
- **How a mode states the band** of the store it serves (§2, item 6). Not specified.
- **The format of the store.** The designer states that compression is to be optimized later (§1, item 16). Not specified here.
- **The alternative equally optimal routes** (§1, item 11). Open, and the designer will return to it.
- **The clean-up** that re-checks results when a lower segment changes. Deferred by the designer (`boundary_board_states_batch_composition.md` §1.15, item 9).
- **The rederivation by primitive composition** of the chain between an anchor and its terminal (§1, item 9). The composition has not been run.

## 4. Status

| Element | Status |
|---|---|
| The statements of §1 | Stated by the designer. The designer is confident it will work, and states that the cost of the trade is not known |
| Items 1 to 8 of §2 | This file's readings, offered for correction |
| The cost model and its scaling | Open. Nothing measured in this file |
| The alternative equally optimal routes | Open, the designer's [A] |
| The rederivation by composition | Not run |

## 5. Provenance

The statements of §1 are the designer's, from the session of 16 November 2026, and the statement at "Where this stands" in `boundary_board_states.md` is the same direction in a shorter form. §2 is this file's reading. This file adds no result of any run.
