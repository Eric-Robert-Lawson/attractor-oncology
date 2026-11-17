# Applications of Partial d-depth Certification: a Front in Itself

**Status.** A record of a direction the designer articulated on 16 November 2026. It records what partial d-depth certification offers, apart from the extension of the Syzygy limits, and the applications the designer sees for it. It is stated in general terms: no piece, material or position is part of it. It adds no primitive and does not touch the hard core (`boundary_board_states.md` §14 items 7 to 11). Statements marked **the designer's** are the designer's. Statements marked **this file's reading** are readings of the designer's words, offered for correction. No measurement is recorded in this file.

Companion documents: `boundary_board_states_calibration.md` (the stored band, the patterns over d, what absence means), `boundary_board_states_start_selection.md` (§2g depth-scoped certification and dependency across segments), `boundary_board_states_batch_composition.md` (the legend, the header and the two modes), `boundary_board_states_syzygy_mining.md` (the certified base).

---

## 1. The designer's statements (16 November 2026)

### 1.1 The front stands on partial d-depth certification

1. **The potential does not wait on the rest.** The potential of partial d-depth certification is carried by the certification itself. The extension of the limits of the Syzygy solution, the inverse ablation and the cataloging of local boundary board states are additions to it. They are not dependencies of what can be applied with the program as it is specified.
2. **A second objective, for the chess community.** The main focus of the program was extending the storage limitation. The opportunity that partial d-depth certification gives to the chess community is a secondary objective of the program, and it is reachable when d-depth is calibrated.
3. **The storage requirement decreases.** Calibrating d-depth decreases the storage requirement. Whether the calibration stays inside the paradigm window and is still useful, and whether it scales up the limits already inherent in the Syzygy solution, is yet to be tested. The designer leans heavily towards it holding.
4. **Segments outside the reach of the tables.** Calibrating d-depth makes the resolution of irreversible segmented boundaries reachable within a useful d-depth outside the reach of the Syzygy tables. The designer is almost certain of this, and it is to be tested.
5. **The program's own trajectory.** The d-depth calibration alters the trajectory of the Syzygy solution into the program's own.

### 1.2 Chess study

6. **A player's own catalog of games.** A player can feed the catalog of their games into the program, to understand the positions they frequently encounter and how to play them better directly. The program points out where they diverge from the certified play, and can show patterns and blind spots directly from their own play. This informs the direction of their study.
7. **One of the first individualizations.** This is one application the designer thought of first. Much more integrated specification can go into individualizing the program for chess study. It is work for another day.
8. **The depth limit is the student's capacity.** The depth that a student studies under is a limit on d, set by the capacity of the student to comprehend the analysis and carry it forward into over-the-board play. It is the same as the calibrated fetch mode: the metric threshold for the depth capacity under which the student wants to study. A table with an unlimited d would give the student rabbit holes, such as cursed wins, that occupy storage and distract from the depths of analysis that matter. Humans and machines both have an inherent capacity limit on the depth of calculation, and the storage limitation adds another, so the d-depth limit is natural and is introduced as a limit.

### 1.3 Play analysis and detecting cheating

9. **d-depth as a measure of reachable certainty.** The d-depth certainty is a calculable calibration metric for the certainty that can be reached to resolve a position. If the program were integrated into a platform that analyzes games, it would show whether a player studied to that depth, and whether a move is reasonable, or an incalculable leap, relative to the certainty that the d-depth provides. Measures that detect cheating can benefit from this.

### 1.4 The test

10. **Access to test at the frontier.** The designer does not currently have the equipment to work at the frontier of the Syzygy solution and to put this to the test. The test waits for access, or for the blockers to be resolved. The direction is recorded here so that it is held until then.

---

## 2. This file's readings, offered for correction

1. **The graded certificate.** Read with `boundary_board_states_calibration.md`, partial d-depth certification gives three kinds of statement about a position in a certified segment: a result certified within d, a conditional draw certified within the d-depth scope (no forced result within d), and no statement beyond the cap. The tables give the first kind only, inside their limits.
2. **What an application can state.** d is certified within an irreversible segmented boundary. What an application can say about a position is therefore what the program has certified for the segment that position is in. A position in a segment that has not been certified carries no statement.
3. **A player's own positions as anchors.** The reading of §1, item 6 uses the enumeration of a player's recurring branches, which the calibration file states (§1, item 13), as the starting anchors and reachability. The divergence is read against the certified result at the position.
4. **The depth limit and the measure of item 9.** The limit of item 8 is the d of the calibrated fetch (`boundary_board_states_batch_composition.md` §1.15, item 6). The d recorded for a position is plies to the first irreversible event of its segment, and a position or a move whose d exceeds the limit lies outside the depth that the user set. How the limit is recorded for a student or for a platform is not specified.
5. **The roles of d found so far.** d scopes a certificate (start-selection §2g), controls what is stored against what is searched (`boundary_board_states_calibration.md`), is the depth limit a student studies under (§1, item 8), and is the measure used to analyze play (§1, item 9). The designer confirms this reading, and states that it is not determined whether these roles are exhaustive.


---

## 3. What is open

- **How a player's games are fed in and mapped** to the segments that have been certified.
- **The form of the report:** divergences, recurring positions, patterns and blind spots. Not specified.
- **How the depth limit is set and recorded** for a student and for a platform (§1, item 8; §2, item 4). Not specified.
- **How the program is integrated into a platform that analyzes games**, and the measure that detecting cheating would use. Not specified.
- **The individualization of the program for chess study** beyond §1, item 6. The designer says it is for another day.
- **A request for a d beyond the stored cap** (`boundary_board_states_calibration.md` §3). Not specified.
- **The test at the frontier of the Syzygy solution** (§1, items 3, 4 and 10). Waits for access.

## 4. Status

| Element | Status |
|---|---|
| The statements of §1 | Stated by the designer. Items 3 and 4 are to be tested, and item 10 states that the test waits for access |
| Items 1 to 5 of §2 | This file's readings, offered for correction |
| The applications of §1.2 and §1.3 | Direction recorded. Not specified further, and not run |

## 5. Provenance

The statements of §1 are the designer's, from the session of 16 November 2026. §2 is this file's reading. This file adds no result of any run.
