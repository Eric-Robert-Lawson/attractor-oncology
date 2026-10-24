# Syzygy as Validation and as the Scope of Certification: A Placeholder

**Status.** A placeholder that holds a direction the designer has articulated, so that it is not lost. Nothing in it has been run. It is a
separate file so that the other documents are not changed by it. If the direction fails, it is removed by deleting this file, and no
other document depends on it. Statements marked **the designer's** are the designer's. Statements marked **this file's reading** are
readings, offered for correction. Section 3 holds what the designer stated and has not settled.

Companion documents: `boundary_board_states_start_selection.md` (§2b the relativistic comparison between local choice regions and the
primitive calculations between two of them, §4 the readings on tests inside the tables and on where the traces meet),
`boundary_board_states_syzygy_mining.md` (the mining process; §6 the test of the automation), `boundary_board_states_opening_survey.md`
(§4f what the tablebase certifies, §5a reachability families, §7 calibration), `boundary_board_state_primitive_catalog.md` (the primitives
and their tempo terms), `boundary_board_states_piece_structure_composition.md` (a related direction).

---

## 1. The direction (the designer's)

1. The Syzygy tablebases can be used to show the coherence of the whole landscape, following invariants that are verifiable within the
   scope of Syzygy. The deterministic potential of the program is then tested from that basis.
2. If the exhaustive coverage that Syzygy needs can be mitigated, so that a basis of optimal-play certainty with certification is
   formed that sidesteps the limitations of Syzygy, the result would be of high impact. This is the aim and not a result.
3. In the forced region the table gives, for a decisive position, the branch that wins fastest, and so the primitive race that wins in
   the shortest primitive distance. The program's assumption about reachable navigation and about the distance to an equilibrium can
   be considered against it.
4. A draw is stated by the table as a draw, with no countable metric for the size of the primitive calculation. That metric is supplied
   by the structure built in this program, or by knowing the Syzygy solution from its exhaustive coverage.
5. **The scope of the enumeration.** Evaluate whether there are positions reachable to the Syzygy tables. The enumeration is then scoped
   to the distance from the predicted local boundary board state to the closest already-enumerated Syzygy table point, for each branch.
   Each branch is checked exhaustively between the predicted state and Syzygy's certainty.
6. The segmented progression can be understood from Syzygy through reachability over the full enumeration of the landscape.
7. The designer identifies this as potentially compatible with the protective belt.

## 2. This file's readings, offered for correction

1. **What is already in the record.** The primitive calculation between two disjoint choice regions, posa and posb, is in start-selection
   §2b. That the calculation can be compared with the table's win, draw or loss value and distance, and extended beyond the tables if it
   tracks them, is in start-selection §4. That a certified entry is a leaf for the structure upstream of it, and that the tablebase
   certifies leaves of refutation lines where they reach seven pieces or fewer, is in survey §4f.
2. **What the table's distance is.** Syzygy gives win, draw or loss and a distance to zeroing, the next capture or pawn move. It does
   not give distance to mate. The python-chess documentation states that the value can be off by one, and that it encodes the 50-move
   rule, with cursed wins and blessed losses above 100. The program does not follow that rule, so the meaning of the distance for the
   program needs checking. The mining file (§5) already notes that the distance to zeroing is not the distance to mate.
3. **Zeroing and the segments.** A zeroing move is an irreversible move, so the table's distance is cut at the same places as the
   program's segmented boundaries. This suggests that segmented progression can be read from the table, as item 6 of §1 says. It is a
   reading and is not worked out.
4. **The scope of item 5, with two precisions.** A branch reaches a table only through a capture or promotion into a held material, so the
   bound per branch is the plies to mate or to the first such capture, whichever comes first. This bounds the refutation side. The
   non-losing alternatives of a state that keeps more pieces than the tables hold stay at that material, and the closure of the region
   has to be enumerated there. The tables a machine can hold set the boundary in practice; start-selection §2a puts it at perhaps five
   pieces. The distance also reads as the cost measure that the calibration of survey §7 asks for.
5. **A validation list inside the tables**, collecting readings already in the record:
   - the sign of the primitive calculation against the table's win, draw or loss, on every position of a material;
   - the length of the winning race against the table's distance, with the care of item 2;
   - that the equilibrium value is zero exactly on the drawn positions that have a losing deviation for each side;
   - the verdicts of the pruning automation against the tables, position by position (mining §6);
   - the piece-addition test and the king-and-pawn example (start-selection §4).
6. **The step from agreement inside the tables to certainty outside them is a separate step.** Agreement on every material the tables
   hold shows that the primitive calculation is coherent. It does not by itself show it stays exact where the tables stop, which
   start-selection §2b names as the programme's worst case. Outside the tables the certificate is still the exhaustive coverage of the
   pruned equilibrium region, which avoids covering the whole landscape.

## 3. Stated by the designer and not settled

- **The profile between posa and posb.** The designer stated that the equilibrium calculation at posa and at posb is 0 or 1, bouncing
  between them with the management of the equilibrium and the tempo; that the calculation for the positions between them rises to a
  maximum and falls back to 0 or 1; that this would be self-evident in the table, because the regions where the equilibrium is not
  established have enumerated branches; and that the values subtracted to find the equilibrium would fit as a deterministic calculation
  of the true primitive values from the forced-win enumeration. The designer adds that it may not scale linearly and that the primitive
  values may jump, in which case the structure can be used to determine how the segmented boundaries and other optimizations of the
  program account for the discrepancy.
- **What this file cannot yet read from it.** What the value 0 or 1 is (for example the tempo margin); which quantity is subtracted,
  since the table gives only the winning side's result and the primitive calculation has a length for each side's race; and whether the
  shape is a gradual rise or a jump followed by a fall of one per move under optimal conversion. These are not settled here.

## 4. Provenance

**The designer's:** §1 and §3. **This file's:** §2, and the list of what is not read in §3.
