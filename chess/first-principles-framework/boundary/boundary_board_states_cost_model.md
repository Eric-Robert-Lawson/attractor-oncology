# Cost Model for the d-depth Conversion

**Status.** A record of the cost model as the designer stated it on 17 November 2026, with this file's reading of its structure. It is stated in general terms: no piece, material, position or figure from any run is part of it. It adds no primitive and does not touch the hard core (`boundary_board_states.md` §14 items 7 to 11). Statements marked **the designer's** are the designer's. Statements marked **this file's reading** are readings of the designer's words, offered for correction. No measurement is recorded in this file.

Companion documents: `boundary_board_states_calibration.md` (the stored band, the patterns over d, the gap, what absence means; its §3 lists the cost model as open), `boundary_board_states_batch_composition.md` (the legend and header, §1.14 and §1.15), `boundary_board_states_syzygy_mining.md` (the certified base that the comparison is made against).

---

## 1. The designer's statements (17 November 2026)

1. **Part of the cost model is emergent.** The savings of d-depth calibration relative to a calibration with no d-depth limit, which is exhaustive coverage, are emergent case by case. The saving is an emergent value relative to the pruning of the resolution that is unreachable under the calibration.
2. **Part of it can be derived with certainty.** The cost of converting the approach from exhaustive coverage to the primitive compositional d-depth schematic metric approach can be derived, in bits, including the cost of storing a position as the Syzygy solution indexes it.
3. **The conversion gives an allowance.** That derivation identifies the optimizations and calculations that give the structurally inherent allowance from which the opportunity cost of the approach can be calculated, the opportunity cost being what the approach allows to be covered in certainty.
4. **A starting basis.** This is the starting basis of the cost model. The d-depth pruning is the main calibration parameter and allows the most savings. The base calculation of the conversion at a structural level gives the overhead storage that is inherent in the conversion, and the savings that the d-depth pruning must supply for the approach to be viable at all.
5. **A fair comparison.** The program's own schematic is to be optimized for storage as far as it can be, so that the comparison is a proper one. The storage of the Syzygy solution is optimized as far as it can be, and the program is a fresh Lakatosian program.
6. **Three fronts of savings.** There is savings potential on multiple fronts: the d-depth exclusion, through the measurement of the distance to the condition, the move chain enumerated and the d-depth limit; the d-depth pruning itself; and the cost of the conversion itself, as it scales with the gap method that resolves narrow, scoped gaps optimally.
7. **The threshold.** Solving for when the opportunity cost presents an opportunity, as far as it can be solved, gives the threshold. How to calculate that explicitly is the work of this model.
8. **The pruning is from both ends, against exhaustive certification.** The d-depth approach sidesteps the exhaustive certification that the Syzygy solution rests on, by certifying partially within what a human can reach. It prunes the calculations that are trivial to reach and the results whose resolution lies above the d-depth limit, and the storage that this frees mitigates the storage limitation. The designer states that the Syzygy solution does not do this.

---

## 2. This file's readings, offered for correction

1. **The terms of a kept row.** For a segment, with a table of reference size T, the program's store holds R rows. Reading §1, items 2 and 4, three terms make up the cost of the schematic:
   - **The key**, which identifies the position that a row is for. In a dense table in canonical order the key is implicit. In a store that omits rows (calibration §1, items 21 to 23) it is not, and what it costs depends on the encoding.
   - **The header** of a row: d, the index into the legend and the symmetry condition (batch composition §1.15, item 7, and the header of the designer's statement of 16 November 2026).
   - **The legend**, once per segment, counted for the entries that the kept rows use.
2. **The base calculation.** With c bits for the key and header of a kept row and L bits for the legend, the schematic is at or below the reference when R·c + L is at most T. The pruning that d-depth must supply for viability is then the fraction of the rows of the exhaustive store that does not have to be kept. This is §1, item 4, written as a condition.
3. **The allowance.** When T − (R·c + L) is positive, it is storage that the conversion releases against the reference. Reading §1, item 3, the opportunity cost is what that storage can be spent on: a higher cap K, more segments, or recorded results that study reaches. This reading is offered for correction.
4. **Which terms are derivable and which are emergent.** Reading §1, items 1 and 2: c and L follow from the schematic and the encoding, and are derivable. T is the size of the reference. R follows from the distribution of d in a segment and from the set of kept values of d (calibration §2, item 9), and it is emergent for each segment. The search cost of a gap (calibration §2, item 10) is a separate term and is not part of the storage terms.
5. **The three fronts and the terms.** Reading §1, item 6: the d-depth exclusion reduces the work of reaching a result, the pruning reduces R, and the gap method trades c·R against the search that re-derives the dropped rows. Whether each front acts on one term or on several is not stated here.
6. **What a fair comparison needs.** Reading §1, item 5: the reference must be counted on the same basis as the schematic. That means stating what T includes (both sides to move, the results that the schematic also holds, and the lower segments that exits depend on), and using on the program's side the best encoding of the key and the header that is available at the time of the comparison.
7. **The two stores deliver different statements.** The reference gives an exact value for every position within its limits. The schematic gives a result certified within the stated d, a draw that is conditional within that scope, and no statement beyond the cap (calibration §2, items 2 and 11). The allowance of item 3 is therefore not a saving of the same thing. It is storage released by a deliberate change in what is certified: the exact value beyond the cap is not stated (§1, item 8), and the rows at the bottom of the band are re-derived by search. For the use of study (`boundary_board_states_d_depth_applications.md` §1, item 8) the values beyond the cap are not wanted. A use that requires the exact value of every position is a use of the reference.

---

## 3. What is open

- **The encoding of the key** for a store that omits rows, and how far it can be reduced. Not specified.
- **The encoding of the header**, including whether d, the index and the symmetry condition can be coded on one another. Not specified.
- **The basis of T:** what the reference includes, and how the two are compared like for like (§2, item 6).
- **The exchange between storage and search time.** Its units, and how the threshold of §1, item 7 is stated across the two. Not specified.
- **R for a segment.** Emergent for each segment (§1, item 1). How it is found for a segment that has not been run is not specified.
- **The search cost of a gap** (calibration §2, item 10). Not measured.

## 4. Status

| Element | Status |
|---|---|
| The statements of §1 | Stated by the designer |
| Items 1 to 7 of §2 | This file's readings, offered for correction |
| The base calculation as a number for any segment | Not in this file. Open |

## 5. Provenance

The statements of §1 are the designer's, from the session of 17 November 2026. §2 is this file's reading. This file adds no result of any run, and records no figure from one.
