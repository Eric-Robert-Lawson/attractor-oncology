# F and P: Empirical Record

**Status.** A developmental record of the empirical application of the F/P workflow, begun on 22 November 2026, kept in parity with the code. Its base is `boundary_board_states_f_p_construction.md`, whose record was brought up to date for a first attempt on 22 November 2026 (§1, items 73 and 74; §4). That file is not changed by this one. This file adds no primitive to the primitive catalog and does not touch the hard core. It is stated in general terms: no piece, material, position or engine is part of it.

It does not accumulate attempts. It records what the development of the code makes understood, and in particular what is not willing to be baked into the code, including as a comment. The code was first written on 22 November 2026.

## 1. The designer's statements

1. **The record is developmental.** What matters is the developmental progress made in creating the empirical application of what the construction file specifies, so that it is understood what is not willing to be baked into the code, including as a comment (construction file, item 74).
2. **The shape of the code.** There are two functions, F and P. Each takes two engine parameters and two priors, or leaves the priors empty if there are none (item 74).
3. **The content is the open front.** The contents of F and P are the programmable goal front and are left open. The integration of the priors and the engines, as the foundational basis, is what proceeds empirically from there (item 74).
4. **What is integrated now, and what is not.** The loop can start from P with a batch or from F, and the designer states that a determined schematic would be a mistake to hard code. Recording engine settings with a batch is not integrated, because anything optional added to a batch or to the priors scales storage. Such properties can be a separate optional file, and only the vital information is held with the priors. Marking engine moves as heuristic or certified is not integrated, and is left to the goal front, since it belongs to an approach that deviates from the proof tree approach and is not a constraint at the lowest level of the architecture. Measuring storage and time per cycle can wait. These were stated on 22 November 2026.
5. **Certification follows from the d-depth limit.** The field for the limit is the d-depth limit for certification. Certification follows from the proof tree once the operation is completed and what to keep is resolved, so a separate certified distinction should not exist. Whether an entry has resolved is shown by whether it carries the index, and an entry at the outer bound of a tree carries d alone (construction file, items 80 and 87). The source field and the engine name are garbage that accumulates (22 November 2026).
6. **The move held at an entry is variable.** The move makes sense, and it can change as the proof tree depth informs the best move. It is determined by the viable pool of moves that satisfy the certification, as the most optimal move within the scope of the d-depth limitation certification (22 November 2026).
7. **The entry holds d and no move is added.** The metric of what d holds is what matters, and d is not one quantity, as the documentation states (construction file, item 79; 22 November 2026).
8. **Symmetry is not held in an entry.** The notation is a redundant storage requirement. One entry covers the folds of its position, and a transformation that lands in another segment gives the same priors for both segments (construction file, items 85 and 86; 23 November 2026).
9. **The d field carries status.** Nothing means not resolved, a d of 1 is the boundary of the proof tree for instance, and a number with the information about the terminal status is certified at the limit (construction file, item 87; 23 November 2026).
10. **Rederivation can use both searches, and the pool of equals can change with the scope.** The search within d and the search by the condition are used together, and the condition is omitted when the pool is to be checked. The check of the pool across a change of scope is left open (construction file, items 88 and 89; 23 November 2026).

## 2. What the code fixes

The code is in `code/` beside this file, as a package of seven parts. Each line below gives what the code fixes and where it comes from in the construction file or the engine integration file.

1. **Two functions with the same plumbing.** F takes an engine and a set of priors for each side and a goal for each side, and gives a batch. P takes a set of priors for each side, a batch, the two engines and a goal, and gives a set of priors for each side (engine integration, item 20; construction file, items 74 and 75).
2. **The goal is a parameter that is the body.** For F the goal is an object whose hooks, a move hook and a done hook, are both optional. For P the goal is the function that is called with the context and returns the new priors. A goal is loaded from stored code by path (construction file, items 76 and 77).
3. **The engine slot.** An engine has a name, settings, and a method that gives a move with an optional score and depth. An engine result carries no certificate (engine integration, readings 20 and 23).
4. **The priors container.** The priors are held in segments, one for each irreversible segmented boundary. A function gives the segment of a position, and it is a parameter. Its default is the key that the documents give: material, with bishops told apart by the colour of their squares, pawn structure and castling rights (construction file, item 83 and reading 72). This is about irreversible segmented boundaries, and not about the outer bound of a proof tree or a boundary board state, which are different boundaries (construction file, item 82). An entry holds the parts of the d-depth header that its state of resolution allows. d is always present. The index is present only where the position resolved to a destination. An entry with d alone marks the outer bound of a proof tree. The index is an opaque value until the legend is specified, and no symmetry condition is held (construction file, items 79 to 81, 85 and 87, readings 69 to 71 and 74 to 76). The key of a position is a function that can be replaced.
5. **The loop.** F, then P, then F on the priors that P gave (engine integration, item 20). It starts with P when a batch is given, and with F when none is (construction file, item 13). The order after the start is not claimed to be fixed, and the loop is a convenience around the two functions, which can each be called alone.
   - **Bringing games in.** A batch can be built from a list of moves, with each move checked as legal, or read from a PGN file or text, with a start position taken from its header where there is one.
6. **What F supplies and what a goal supplies.** F supplies the board, the legal moves, the end of a game, and the batch with duplicate games dropped. A goal supplies the choice of move, and when the batch is complete. The number of games is therefore a goal input (construction file, items 24 and 36). Both sets of priors are given to every goal, so visibility is the goal's choice (item 6).

7. **What P builds from a batch.** P merges the board states of the lines of the batch by position into one structure, with the moves played between them, the places where lines begin and end, and how many lines pass through each state. It gives this to the goal with the batch. It is working data for one call and is not stored with the priors (construction file, items 39, 44 and 67).
8. **The shape of what P gives back.** P requires a pair of priors, one for each side, and refuses anything else. Every entry must have a d, a whole number of zero or more, and an entry with no d is refused because what cannot be certified is not added. The check is of shape only and does not judge what a goal keeps or omits (construction file, items 2, 4, 28, 43 and 81).

**Baseline choices made by the code.** Each stands in for an open item and is not an answer to it.

- When a goal leaves the move to the default, the engine chooses. The priors are not consulted by the default, and how a goal uses them is the open item on the query surface.
- Duplicate games are dropped within a batch only, and not against the priors.
- A game ends at the rules' end of game with a draw claim, or at a cap on plies that is a parameter.
- The baseline P adds nothing, since no check has resolved anything and an entry with no d is not added (construction file, item 81). It replaces the baseline P that kept every line, which engine integration, reading 20, called the simplest form.
- The batch is a list of games, each a start, a list of moves, and a result. This stands in for the open item on the batch format. Nothing else is held per game.
- The key of a position is its text without clocks. This stands in for any definition of a boundary.
- The shape check refuses an entry with no d, and requires a d to be a whole number of zero or more. The segment function stands in for the assignment of a position to an irreversible segmented boundary, and uses the key of start-selection §4 and the topology file, which are readings of the designer's statements, with bishops told apart by the colour of their squares as the designer states (construction file, item 83).
- A first version held a limit and a source on every entry and the names of the two engines on every game. These were the code's own additions, with no statement behind them, and the limit was placed on each entry where item 2 of the construction file has it for each irreversible segmented boundary. They were removed on 22 November 2026.

## 3. What is kept out of the code

1. **The exact check.** Nothing in the framework performs a check. A check is what a goal for P supplies.
2. **A keep or omit rule.** None beyond the baseline that keeps everything.
3. **A unit for the window.** There is none, so nothing counts against a window yet.
4. **A boundary definition.** The key of a position is replaceable, and no boundary is computed.
5. **Any filtering of priors by side.** Both sets are always given to a goal.
6. **A fixed number of games.** The count is held only by a goal that supplies it.
7. **A hosted fetch of goals.** Goals are loaded from a local path. Fetching from a server is left for when it is wanted.
8. **A layout for storing the priors.** The container is held in memory and is not written to disk.
9. **Engine settings recorded with a batch.** Not recorded, because of the storage scaling of anything optional added to a batch or to the priors (§1, item 4).
10. **A heuristic or certified mark on an engine's move, and a certified field on an entry.** Not added. The first is left to the goal front, and the second is the presence of d (§1, items 4 and 5).
11. **Measurement of storage and time.** Not added yet (§1, item 4).
12. **The limit of the proof size that the priors emerged under, apart from d.** Not held. d is the d-depth limit for certification. Item 2 of the construction file also has the proof size limitation recorded for each irreversible segmented boundary, so where that is held, if at all, waits on the boundary key.
13. **A record of where an entry came from.** Not held. No statement asks for it, and anything optional added per entry scales storage (§1, item 4).
14. **A move per entry.** Not held (§1, item 7).
15. **The legend of a boundary.** Not in the code. The index is an opaque value with no meaning enforced, and a goal for P that resolves a position supplies it (construction file, items 81 and 87).
16. **Entries with no d.** Refused by the shape check (construction file, item 81).
17. **Transformations by symmetry.** Not in the code. An entry holds no symmetry condition (construction file, items 85 and 86). The folds within a segment and the equivalence between segments are a property of the key and are not implemented: the segment function and the key of a position do not put a position in a shared form (construction file, reading 75).
18. **A search to rederive a line.** Not in the code, which has no search of its own. Whether it goes by d, by the condition or by both is not fixed (construction file, item 88).
19. **A check of the priors for the pool of equals across a change of scope.** Not in the code and left open: to record the pool, to run a recurring check, or neither, and whether the check is the clean-up of batch composition, a function of P, or P with a goal. It is a mechanism that a goal may supply (construction file, item 89 and reading 77).

## 4. Readings

1. **Divergence came from the engine's state, and the code does not provide it.** This file's reading of a run with an engine that speaks UCI at fixed settings, offered for correction. The games of one batch did not repeat, so no game was dropped as a duplicate. This is read as an effect of state the engine carries between games, and not as a means of divergence that the code gives. Engine integration, reading 20, leaves the means unspecified, and the code leaves it to the goal. If an engine is run with no carried state, repeated games would be expected, and this has not been run.

## 5. Status

| Element | Status |
|---|---|
| Code | Written on 22 November 2026: the engine slot, the priors container, F, P, the loop and the loading of stored goals |
| The loop | Runs end to end at small scale with two trivial test engines, a null goal for F and the baseline P. Plumbing only: no check has been performed and no entry has resolved |
| The wrapper for an engine that speaks UCI | Run on 22 November 2026 against an engine installed from the workspace's system package source. It gave a move, a score and a depth, and the loop ran for three cycles with it at two different depths for the two sides, with the null goal and the baseline P that kept every line, which has since been withdrawn. Plumbing only, as above |
| Starting from a batch | Integrated on 22 November 2026, with a start from P on a batch or from F. Tested with a hand-made batch and with games read from PGN text. Plumbing only |
| What P builds and requires | Integrated on 22 November 2026: the merged structure of a batch, given to the goal, and a check of the shape of what a goal returns. Tested on lines that share a prefix and that transpose, and with outputs of the wrong shape. Plumbing only |
| The shape of an entry | Changed on 23 November 2026: the symmetry field and the rule that tied it to the index were removed after items 85 to 87. An entry holds d and, where resolved, the index. Plumbing only |
| What the code fixes | Recorded in §2 |
| What is kept out of the code | Recorded in §3 |
