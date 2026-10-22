# Boundary Board States: Mining the Syzygy Tablebases (the process)

**Status.** A process record, not a result and not a tool. It holds the process the designer is using to build the boundary board
state certification catalog (survey §4j) from the Syzygy tablebases: enumerate the proven draws of a material, classify them, group
them into reachable families, and link the families upward toward the initial position. Nothing in it has been run or verified as part of
the record; the designer will run and verify the data on return. The same process is held, as a specification, in
`boundary_board_states_start_selection.md` §2a, which is the fuller and later record. No primitive calculation is part of the
process (§1). The designer has said that another angle is under consideration which may change this process; the staged workflow
from an anchor in start-selection §2b is the one recorded since.

**Companion documents.** `boundary_board_states_opening_survey.md` (§4j the two catalogs, §4k the designer's statements behind this
file, §5a reachability families and the classification of a family), `boundary_board_state_topology.md` (§8ar, the record of the
statements), `boundary_board_states_start_selection.md` (§2a, the later specification of this process, and §2b, the staged workflow), `boundary_board_states.md`
(changelog entries 62, 63 and 65), `boundary_board_state_primitive_catalog.md` (the primitives, used afterwards, §1).

**How this file is organized.** What is the designer's statement is marked as such and kept apart from this file's readings, in the
manner of the catalog's §10. The steps of §3 describe the procedure, and where a step is generalized to other materials the
generalization is marked as this file's reading and is offered for correction.

---

## 1. What the process is for (the designer's statements)

- **The catalog from proven results.** The designer's aim is to begin the boundary board state certification catalog by creating it
  from the proven results of the Syzygy tablebase solution set, enumerating the boundary board states from it. Getting the
  population of draws from Syzygy is important.
- **Enumeration by mapping the tablebase onto the branch.** A branch can be enumerated by taking the tablebase's win, draw or loss
  value (and distance) and mapping it onto the actual branch with its legal move pool, to determine which branches align with the
  count.
- **Optimal play is the emphasis.** The moves chosen for branch creation are chosen under optimal constraints.
- **Primitive calculation is not part of the script.** The primitive calculations are done when a drawn position is verified later,
  by the designer. Primitive compositional analysis is done after the fact, from the population of certified local boundary board
  state positions obtained by exhaustive coverage from Syzygy.
- **The 50-move rule is ignored for now.** It matters only inside a reachable family, when the anchor is reached, as a cutoff.
- **How a draw is realized.** The draws are categorized by whether the paths draw by stalemate, by insufficient material, or
  cyclically. Where the verified best move from the tablebase is to reduce to insufficient material, that is a provable equilibrium
  of reduction to insufficient material, rather than a stalemate choice, since optimality is the concern.
- **What a boundary board state position is, here.** A position in the choice region that has a move that can lose is a boundary. The
  initial game start position is itself a boundary board state position. A distinction between the choice region and the boundary
  by how many drawing moves a position has is a meaningless classification criterion (survey §5a).
- **The two spaces.** There is a literal boundary between force and choice, and understanding that boundary, where the forced draws
  are, is meaningful. There is a choice space and a forced space, which together encapsulate the entire space of chess, with areas
  of sub-optimal play pruned. The pruned structure that emerges from the process is the full solution of chess.
- **Whether the boundary is a majority of the space.** The designer wants to confirm that the boundary board state understanding
  covers a majority of the entire space. This is the reason a test on one small material was not integrated into the understanding on its own.
- **Syzygy size and the automation.** The size of the tablebase does not matter eventually, for reachability, with the program's
  methodology of backtracing and exhaustive enumeration from the cached database file. The work is to mine all the Syzygy data the
  designer's machine can handle, and then to automate off that basis, from the catalog, with the pruning mechanism.
- **The order.** Syzygy tablebase mining first; then a boundary board state catalog and reachable families already in hand from the
  mining; then the automation of the process, the matching of primitives, and so on.

## 2. The two spaces, in the terms of the process

This section is this file's reading of the statements above, offered for correction. For a material held in the tablebase, every
legal position is in one of two spaces. The forced space is the positions whose value is a win or a loss for the side to move. The
choice space is the drawn positions. The boundary is the set of drawn positions that have at least one move that loses, the
positions from which a deviation leaves the choice space for the forced space. The rest of the choice space, the drawn positions
with no losing move, is called the interior here. The designer did not name it, and the term is this file's. For a material
small enough to be held completely, the tablebase is the complete solution of the material and the two spaces are exact.

## 3. The process, step by step

1. **Choose a material.** A set of pieces, within what Syzygy provides and what the designer's machine can store. The tables needed
   are those of the material and of every material reachable from it by a capture; a bare king against a bare king needs no
   file, and is a draw by insufficient material.
2. **Enumerate the legal positions, both sides to move.** The squares are distinct, the kings are not adjacent, and the side not to
   move is not in check. Fold by symmetry. For a material with no pawns, the eight symmetries of the board together with the swap of
   colours (which also swaps the side to move) give sixteen images, and a position is represented by the image with the smallest
   index. For a material with pawns the group is smaller; that is this file's reading, and has not been run.
3. **Probe the value of every representative.** Win, draw or loss for the side to move, taken by sign. A cursed win or a blessed loss
   counts as a win or a loss, since the 50-move rule is ignored (§1).
4. **Classify every legal move of every drawn position, under optimal play.** A move that loses is a deviation, certified by the
   table, and is counted and not followed. A move that keeps the draw is an optimal move, and is either a move that stays in the
   material, which is an edge of the graph of §5, or a capture that leaves the material, whose value is read from the lower
   material's table. A capture is a terminal exit when it leaves the opponent stalemated, or leaves a bare king that can take the
   last piece and so reach a bare-kings draw. No drawn position has a winning move, and the process checks that.
5. **Check the tables against the move generator, on every position.** The value of a position must equal the best of its children's
   values, negated, or the mate or the stalemate where there are no moves. This is a check of
   the generator against the tables, and does not test anything about the programme.
6. **Build the graph of the drawing moves.** The nodes are the drawn positions, the edges are the drawing moves that stay in the
   material. The strongly connected components are the sets of positions that can reach one another by drawing moves, and the
   weakly connected components are the families when direction is ignored. A position lies on a cycle when its component has more
   than one position or it has a self-loop; a cycle in the folded graph is a cycle in the real graph, because the symmetry that
   maps a position to its image has finite order. For each position, record whether a cycle, a stalemate, or a bare-kings draw can
   be reached by drawing moves, and the realization class that follows: *stalemate now*; *terminal only* (no cycle can be reached,
   every drawing line ends in a terminal draw); *cycle only* (no terminal draw can be reached); or *mixed*.
7. **Store.** The population, with the fields of §5.
8. **Verify independently.** Recount the moves of a sample of positions with a second implementation (python-chess), check the
   legality and the value of a random sample of raw positions, and recompute the components and the reachability flags with a
   second implementation (scipy). The sample sizes are for the designer to choose.

## 4. The classification the process applies

The designer's criterion is that a drawn position with a move that can lose is a boundary board state position (§1). The process
applies it to the side to move, from the stored counts. The designer's criterion has two sides, and whether the opponent's losing
deviation is also tested at the same position is not settled by the statement (§9). The interior of §2 is not a boundary board
state position under the criterion.

The realization of a draw (§3, step 6) is recorded separately from the criterion, as the designer asks (§1). Reduction to a bare
kings draw is a forced class when every drawing move of the position is a capture to a bare kings draw.

## 5. What the process stores, and what it leaves out

**Stored** (proposed), one row per drawn class: the position, the side to move, the number of positions the class stands for, whether the side
to move is in check, whether the position is itself a stalemate, the number of legal moves, losing moves, drawing moves, drawing
moves that stay in the material, and drawing moves that end in a stalemate or a bare kings draw by capture, the strongly connected
component and its size, the weakly connected component, whether the position is on a cycle, and the reachability flags and the
realization class of §3. These are the fields of the process's database. They are not a proposed template for an entry of the
certification catalog, and no template is proposed here.

**Left out**, so far (this file's reading): the edges themselves, the identity of the losing moves with the child each leads to, and how fast each loss is
refuted (the `.rtbz` files give the distance to zeroing, which is not the distance to mate where captures are possible). Whether
to store them is the designer's decision.

## 6. From the mined catalog to families and to the anchor

**The designer's statements, recorded in survey §5a.** If a boundary board state can reach another, they are in the same reachable
family. A family confirmed joint to the initial start position under the optimal backtrace constraint is an optimally reachable
family, and every member is global. When a state is discovered by construction there is no backtrace constraint at discovery, and a
family is then composed that either reaches the initial position or concludes that it cannot. Positions upstream that lead to a
boundary under optimal play share its equilibrium.

**This file's reading, offered for correction.**

- **Families from the mining.** The components of step 6 are the families within one material. Whether the strongly connected
  component, the weakly connected component, or something else is the designer's *reachable family* is not stated; both are stored.
- **The upward link.** A mined family is linked toward the initial position by backtracing under optimal play into materials with
  more pieces: positions of a larger material from which optimal play reaches a position of the mined family by a capture. That
  step needs data for the larger material. Where the larger material is within the tables the same process applies. Beyond them it
  is the automation of the pruning mechanism, whose leaves are the certified entries (survey §4f, the proof object), and which is
  not specified here.
- **A test of the automation.** On every material where the tablebase has the answer, the pruning mechanism's verdicts can be
  compared with the tables position by position, so the automation can be checked before it is trusted on materials the tables do
  not reach.
- **The 50-move rule** enters at this stage and not before, as the designer states (§1).

## 7. The order of work

As the designer states it (§1): mine all the Syzygy data the machine can handle; have the boundary board state catalog and the
reachable families from that mining; then automate the process off that basis with the pruning mechanism, matching the primitives
afterwards. This agrees with survey §4h (the order of work) and §4f (what is cheap).

## 8. Scale (this file's readings, offered for correction)

These are estimates and have not been checked against a run. Each added piece multiplies the count of legal positions by a factor of
tens, so a five-piece material is likely to be in the tens of millions after folding, more where pawns reduce the symmetry, and a
table of edges would be larger by about the average number of drawing moves. The sizes of the Syzygy sets (a gigabyte or so for
three to five pieces, tens of gigabytes for six, terabytes for seven) are from memory and not checked. Probing in a scripting
language is likely to be the slow step and may need a faster tool for larger materials.

## 9. Not decided

- Whether the boundary criterion tests the opponent's losing deviation at the same position, or only the side to move's.
- Whether to store the edges, the losing moves and the refutation distances (§5).
- Which components are the reachable families (§6).
- The fields of an entry of the boundary board state certification catalog; none is proposed.
- How the upward link is carried out beyond the tables (§6).
- The other angle the designer has said is under consideration.

## 10. Files

No script, database or result from any earlier test is carried over into this record. The designer's own files for the first
material to be mined, `KRvKR.rtbw` and `KRvK.rtbw`, and the designer's solver and validator, are the designer's, and the process is to
be built and verified from them on return. The procedure of §3 is specific to the pieces named at each step (step 2 assumes a
pawnless material; other materials need their own move generation).

## 11. Provenance

**The designer's:** everything marked as a statement in §1 and §6. **This file's:** the steps of §3 as generalized, the terms
*interior* and the realization classes, the readings of §2, §6 and §8, and the list of §9.
