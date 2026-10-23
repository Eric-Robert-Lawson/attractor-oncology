# Piece-Structure Composition: A Proposed Direction

**Status.** A record of a direction the designer has articulated. It is not worked out and nothing in it has been run. It is kept as its
own file so that the other documents are not changed by it. If the direction fails, it is removed by deleting this file, and no other
document depends on it. Statements marked **the designer's** are the designer's. Statements marked **this file's reading** are the
reading of the designer's words, offered for correction.

Companion documents: `boundary_board_state_primitive_catalog.md` (the primitives; §4.10 the composition interface; §5 evaluating a
candidate; §8.9 the support closure; §10 the table "Edge and corner draws"; §11 the rule set),
`boundary_board_states_start_selection.md` (the ablation and inverse ablation of pieces),
`boundary_board_states_opening_survey.md` (§5a reachability families), `boundary_board_states_syzygy_mining.md` (the source of
certified local boundary states), `boundary_board_states.md` (the concept and the definition).

---

## 1. The direction (the designer's)

1. A local boundary board state position can be encapsulated in a region smaller than 8x8. The piece structure is the pieces together
   with the region they occupy.
2. The board edge is not removed. It exists, and the structure must exist within the 8x8. Its placement within the board region matters
   because the edge closes the region and so informs the calculation.
3. The primitive calculations of a piece structure can be examined as the structure is displaced within the board region, and as it is
   placed in composition with other piece structures.
4. The result is two catalogs: a catalog of local boundary board state positions, and an accompanying catalog of the primitive
   calculations of those piece structures with respect to their placement within the board region.
5. With known primitive structures that have proven local boundary board state positions, the work is to find out whether piece
   structures compose, decompose, and are jointly compatible, and to formulate a process in the program for the certification checks
   that navigate this. If the structure of that process can be understood with certainty and the direction proves useful, it is a
   potential addition to the protective belt.

## 2. The anchor and compatibility (the designer's)

1. The terminal anchoring point of the mating net is the king. The anchor informs the segmented dependency calculation and the coverage
   net.
2. Piece structures can share the king in the same place, and then overlap at the anchor.
3. Piece structures with a different king placement are not jointly compatible by definition, because the anchoring point of the mating
   net primitive differs.
4. Joint compatibility has two parts: the joint compatibility of the primitives, and the joint compatibility with the composition of the
   rule set from the initial position. The second concerns the piece set under the rules of chess. A board state with three knights
   is not compatible with the initial piece setup, and is reachable only through a pawn promotion to a knight. It does not mean the
   anchoring of a global boundary board state to the initial position.
5. Sharing the king's square as the anchoring point is the minimum requirement for royal structures to be compatible. It is not
   sufficient: the other checks still apply.

## 3. Perturbation and ablation of a structure (the designer's)

1. Where a side has more of a piece than the usual number, for instance more than two knights, the integration of the piece structure
   can be understood, and the extra piece can then be removed or added as a perturbation of the board state.
2. What is examined is how the removal changes the composition of the primitives and their compatibility. Single pieces are the granular
   case under the rule set of perturbation. The focus here is structures of pieces and how to decompose and compose them.
3. The whole is a system-state manipulation schematic: the composition and decomposition of the joint compatibility of the primitives
   and of the rule set, and of the composition and decomposition of the piece structure.

## 4. The protocol as stated (the designer's)

1. Gain a catalog of boundary board state positions (the designer runs the mining).
2. Calculate the primitive calculations on those positions.
3. See whether the approach is viable.
4. Keep the direction on record, to be returned to when there are sufficient requirements to try it empirically and discover the
   emergent relationships between compositions of the primitives.

---

## 5. This file's readings, offered for correction

1. **The edge as an input.** The catalog already treats the board edge as missing edges of the movement relation, an input of P-003 (§10,
   "Edge and corner draws"). I read the designer's statement as extending that to a placement: the same structure at another
   placement has another edge closure, so it can have another calculation.
2. **Placement classes from the support closure.** A verdict depends on the board through the support closure of the structure (§8.9).
   Where the closure at one placement does not meet the edge, a displaced placement with the same relative closure would share the
   verdict and the proof. Where it meets the edge, the placement would have its own entry. This is a reading, not shown.
3. **A sketch of the checks for a composition certificate.** The closures of the two structures do not interfere at their relative
   offset. Neither structure's terms change through the other's occupancy or coverage (§4.10). The edge is handled in each closure at its
   placement. The combined material is reachable under the rules from the initial position piece set (Law 3.10.3: a position that
   cannot have been reached by any series of legal moves is illegal), applied to the combination. The anchoring of a global boundary
   board state to the initial position (`boundary_board_states_opening_survey.md` §5a) is a separate matter and is not part of this
   check. By the ruling recorded in the primitive catalog §11.2, a combination whose material is not reachable is kept as information
   and tagged, and does not enter a global closure.
4. **Tempo.** The obligation to move means a move in one structure is a free move for the other. An equilibrium that holds because the
   side to move must move may not hold in a composition that supplies a spare move. The reading is that equilibria divide into those that
   depend on tempo and those that do not, and that the composition check needs to know which. It is not shown.
5. **Reachable material.** An extra piece of a kind arises by promotion, which consumes a pawn and needs captures for the pawn to change
   file. This bears on whether a composition with an extra knight is reachable from the initial position. It is not worked out here.
6. **Fields of a composition record.** The anchor; the region and the placement; the support closure; the primitive calculations; the
   material required and whether it is reachable; whether the equilibrium depends on tempo.
7. **The status of non-royal structures.** Structures with no king in them, such as a capture net on a non-royal target or a chain, would
   be compatible with a royal structure so long as the other checks pass.

## 6. What is open

- **The anchor in a repetition cycle.** The king may alternate between squares in a cycle. Whether the same anchor means the same king
  square at the same point of the cycle, or the king's orbit through the cycle, is not stated.
- **Viability.** Whether the direction is useful is to be decided on a catalog of proven local boundary board states and the primitive
  calculations run on them. The primitive calculations of P-003 and P-008 are specified and not run (primitive catalog §5), and the cost
  model is unfinished (§4.8). The direction waits on those.
- **The composition check for royal structures** where the kings coincide at the anchor but other pieces differ.
- **The tempo-dependence of equilibria** (§5 item 4).

## 7. Provenance

**The designer's:** §1 to §4, the anchor and its compatibility, the perturbation of a structure, the protocol, and the decision to keep
the direction as a separate file that can be removed. **This file's:** §5 and the open questions of §6.
