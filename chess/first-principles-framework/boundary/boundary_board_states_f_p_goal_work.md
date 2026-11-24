# F and P Goal Work File

**Status.** A work file begun on 22 November 2026, for the goal of the exact check that P gives to F, and for making the F and P workflow executable. It takes `boundary_board_states_f_p_construction.md` as its base and does not change it. The work that carries a risk of failing is kept here, so that the construction record is not affected by attempts (construction file, item 23 and reading 19). It adds no primitive to the primitive catalog, does not touch the hard core (`boundary_board_states.md` §14 items 7 to 11), and belongs to the protective belt. It is stated in general terms: no piece, material, position or engine is part of it. Item and reading numbers refer to the construction file unless another file is named. Statements marked **the designer's** are the designer's, and statements marked **this file's reading** are readings offered for correction. Nothing in this file has been run, and no attempt has been made yet.

---

## 1. Carried forward from the construction file

The first two blocks and all the open items are copies of the construction file as it stood on 22 November 2026. The copies are unchanged, and the construction file keeps the originals. Where a copied text cites an item or a reading, or says "this file", it means the construction file.

### 1.1 The requirements on the goal for P (construction file, reading 43)

Note, not a copy. Item 72, made after reading 43 was written, states that the goal of P has a general preset and can be as variable as the goals of F. Where the copy below says that the goal is constant, it is read as the preset (construction file, reading 64).

The copy, unchanged:

43. **The requirements on the goal for P, collated.** This file's collation of statements that are already recorded, offered for correction. It adds no statement of the designer's, and each line gives its sources.
    - The goal is constant. P gives it to F when a proof tree is wanted. The goals that the agent programs into F are a separate matter (items 22 and 41).
    - The scope is the differential of the batch relative to the existing priors, from the anchoring points of the branches of the batch (items 30 and 41).
    - The batch is checked as a sum within the window of the machine, and whether it is one operation or one game per operation is a parameter (items 36 and 44).
    - The order considered is from the ends of the branches upward, beginning near the edge of exhaustive certainty, with positions inside the reach of the tables left to the tables and the engines (items 31, 32 and 33). It is recorded as a consideration and not as a procedure.
    - The check covers every reply of the side that does not choose, within d. The side that chooses follows the priors, and the engine is used where the priors hold none. The check ends at certified anchors, which reset the window (items 22, 26 and 29; readings 12 and 17; engine integration, reading 22). Reading 12 is the baseline form and does not fix the construction (item 40).
    - The check stops at closure or at the storage limit of the check, and it records the d reached and the proof size (items 28 and 29).
    - The output is a proof tree in which the refuting lines are the main addition. It is integrated cumulatively into the previous priors, and the lines of the batch need not be kept (item 30).
    - The proof size is the size of the certificate. The effort of producing it can be larger, and its efficiency matters because it recurs on every batch (items 41 and 42; reading 38).
    - At the nodes where the certified side chooses, the certificate holds one winning child, and the alternatives are not stored. Whether such cases are recorded depends on the balance (items 45 and 47; reading 42).
    - Gaps are re-derived within d from certified anchors, and that cost accumulates in each check that crosses them (item 42; engine integration, reading 22).
    - The priors keep certified entries. Heuristic entries are not kept after a check. Keeping an entry depends on the cost of re-deriving it within the balance of storage, coverage and speed, with the frequency of crossing as a consideration (items 34, 37, 43 and 47).
    - Left open: the construction of the goal, including how the evaluator is used (item 40), the unit of the window, the query on the priors, the batch format, and the reset logic (§3). How "already covered" is decided for the differential scope of item 41 is addressed in part by items 48 to 51 and reading 47, and the rule that compares what is stored with what a check seeks is not stated.

### 1.2 Where the evaluator enters each vector (construction file, reading 58)

58. **Where the evaluator enters each vector, collated.** This file's collation of what is recorded, offered for correction, after item 63. It adds no statement of the designer's, gives the sources of each line, and does not decide how the evaluator affects any vector.
    - **The engine slot.** The engine is a modular evaluator with the widest compatibility (items 3 and 11). The layer needs only an engine that runs locally, speaks UCI, evaluates to a limit and reports its depth (engine integration, reading 23). The exactness of a check does not depend on which engine produced a line, so engines can be exchanged (engine integration, reading 14). The optimal priors may depend on the engine in the slot (item 63).
    - **The goals.** The engine is used as an evaluator relative to the priors (item 5). In the layer the goals are hooks, which include which nodes are expanded, how a leaf is valued, how the opponent is modelled and when to stop (engine integration, reading 23).
    - **F.** The engine is a source of moves in the self-play (items 5 and 11), and which priors it sees is set by the goals (item 6).
    - **The exact check.** Its exactness does not depend on the evaluator (engine integration, reading 14). The baseline has the engine propose moves for the certified side where the priors hold none (reading 12), and how the evaluator is used inside the check is left open (item 40).
    - **The priors.** An entry carries a certified or heuristic mark, and an engine's finding is certified only once the line is checked against every reply (calibration, reading 13; engine integration, reading 20). The record of a certificate names the engine, its version and its settings (engine integration, reading 14). Whether certain priors are important at all may depend on the evaluator (item 63).
    - **P, keeping and omitting.** The cost of deriving a result decides keeping, with the example of a short mate that an engine finds easily (item 34), and the use application has an engine resolving the gap at play (items 38 and 43). That the keep decision may depend on the engine that crosses the gap is reading 56, a possibility that is not part of the definition of keeping (item 64), and the designer has not stated it in those words.
    - **The window and d.** An engine calibrated by d-depth turns heuristic search into partial certificates (engine integration, §1, item 2) and reports the depth that it reached (engine integration, reading 23). How the depth of the evaluator relates to the d of a check is not recorded.
    - **The balance.** How far apart the stored states can be depends on what an engine can resolve between them at an effective speed (items 38 and 43). Whether that makes the distance depend on the engine is open (§3).
    - **The two applications.** In progression the evaluator serves the construction of the priors, and in use it navigates the gaps (items 43 and 60 to 62).
    - **Not recorded.** How the evaluations should impact the priors and the exact check (item 59), and the role of the evaluator, which is not exactly clear throughout in the designer's mind (item 63). Open in §3: which engine the keep decision is made for, whether the priors are held per engine, how the depth of the evaluator relates to d, and what "provably better than the native line" means.

### 1.3 All open items (construction file, §3)

Note, not a copy. **The designer's:** these open items may or may not be important to establishing the goal front for P and F, and whether they are cannot be told now. They may also be important to how that front is implemented as a variable that can be modified for F, and possibly for P as well, since separating the optimization from the architecture (construction file, items 69 and 70, reading 63) may influence the calibration of P. The designer is not sure of this and keeps it open. **This file's reading, offered for correction:** no open item is read as unimportant, and none is excluded. Item 72, made after the first version of this note, states that the goal of P has a general preset and can be as variable as the goals of F. Parameters such as the working store and the lifetime of reading 63 are therefore read as settable for P as well, with the preset giving their defaults. Which defaults the preset has is not stated.

The open items, copied unchanged:

- **The query surface of the priors** (§2, item 3): which queries the goals can make, and what each returns.
- **How a goal is expressed** (§1, item 5; §2, item 5): what a goal is given, how a side's goal is passed to the framework, and what the framework may not interpret.
- **How visibility is specified** (§1, items 6 and 7): how it is set for each side and for each goal, and how its effect is tested.
- **The measurement of the optimality of a batch relative to a prior** (§1, item 10): what is counted, what counts as bias into bad positions, and how it is recorded.
- **Where differential testing sits** (§2, item 6): inside P or as iterations of F, and how much of it the cost of re-derivation replaces (§1, item 42; §2, reading 39).
- **Where the attempted d is set** (§2, item 2): in the iteration, in the goal, or elsewhere, and how the proof size is recorded with it (engine integration, §3, the measure of the space to resolve).
- **The batch format and what the batch records** (§2, item 5).
- **The baseline P** (§1, item 9): the form in which a batch of games becomes priors, as the simplest case.
- **How the re-derivation cost is recorded** (§2, item 9): the proof size, the time, or both, and at which d.
- **How the static measures of P combine into a rule for keeping or omitting** (§1, items 14, 31, 33 and 34; §2, readings 9, 34, 36 and 37): including how the speed tolerance is set or left to emerge, within the balance of storage, coverage and speed (item 34, reading 37).
- **How a promising position that is unresolved is recorded and prioritized** (§2, item 8).
- **When the limitation is raised, and by how much** (§1, item 15; §2, items 10 and 11).
- **The perimeter of a proof tree** (§1, items 49, 50, 52, 53 and 56; §2, readings 44 to 46, 48, 49 and 51): the designer views it as the lowest edge of a map of certification limits and leans towards storing it at d-depth 1, with the d-depth 0 states derived from it, and leaves the level open to the balance of storage, speed and coverage. Open: whether the states at which the tree ends at a certified anchor are outside it (readings 48 and 51 read them so, and item 57 does not address it); whether the certified move is held at the entries where the certified side chooses, which the derivation of the zero states needs; whether all of the edge or part of it is stored; and how it relates to the interior states that are kept.
- **How "already covered" is decided for the differential scope** (§1, items 41 and 48 to 51; §2, readings 43 and 47): the designer states that the understanding is closed and not yet fully articulated. What is stated is the perimeter of the tree and the role of batches. What is open is the rule that compares what is stored with what a check seeks, and what a check does at a perimeter state. Item 54 states that how a later check proceeds from the perimeter varies with the programmable goal, so the rule may be a baseline and not fixed.
- **How a goal that yields a proof tree relates to the two families of goals** (§1, items 24, 41 and 54; §2, reading 50): item 41 separates the variable goals of F, which do not result in a proof tree, from the constant goal of P, which does. Item 54 says that a goal could be a proof tree in itself and that F could be recursive with the properties of P. Whether that is the constant goal of P called from within F, a variable goal that yields a proof tree, or something else, is not stated. The designer states that it should not be excluded and can only be answered empirically (item 65). Item 72 adds that the goal of P has a general preset and can be as variable as the goals of F.
- **The check as a goal** (§1, items 19 and 29; §2, item 12): the designer states that this is unresolved and has not been attempted. It covers the form of a batch that is a proof tree, how the enumeration of every reply is requested in a goal, and how the proof size is counted. Item 56 states the objective as the coverage of the d-depth limit certification, which item 57 confirms is formally articulated as raising the lowest edge.
- **Whether work retained from a check counts against the window** (§1, items 29 and 42; §2, reading 39): or against a separate working allowance, and for how long it is held. Item 67 states that repeated work can be composed on transpositions through compression, which reading 61 reads as reducing repeats and not the number of distinct states examined. Item 68 raises whether sections of the working memory can be dumped and whether the evaluator can make re-derivation cheap enough to allow it. The designer does not know, and reading 62 gives a reading that is not tested. The designer classes the working memory as an open front for optimization more than architecture (item 69), and reading 63 reads the goal as stated with this as a parameter.
- **How the batch is split into operations** (§1, item 44; §2, reading 41): one operation for the sum or one per game, how the merge takes the composition of the games into account, and how merging two sets of priors relates to merging a batch into existing priors. Item 51 recalls this item and adds the comparison of the priors of batches that resolve complementary regions.
- **The repeated failed attempts at a larger d** (§1, item 46): how the program avoids repeating a lost cause, and where any record of a failed attempt is held, whether in the priors, which adds storage, or outside them as bookkeeping of the progression. The designer leaves this to be discovered empirically and added to the protective belt. Reading 47 raises a possible link to the perimeter of item 49, which the designer has not stated.
- **The efficiency of the recurring check** (§1, item 22; §2, item 18): how the check configuration is made inexpensive per call, since it recurs on every cycle. Not specified here.
- **The stored form of the d recorded for a board state** (§1, items 28, 52 and 57; §2, readings 21, 29 and 52): what is recorded is given. The d is a value held at a position, it is the certification limit of that position, and it is counted in layers of legal plies from the anchoring points. What is not specified is the layout in which it is stored, and its relation to the header of the earlier d-depth metric statement, which item 28 says does not apply.
- **How the coverage of the priors affects the reset points** (§1, items 25 and 26; §2, readings 23, 25 and 26): item 26 gives the mechanism, in which the reset leaves an allowance that widens or deepens the search. What is open is how the shape of the coverage of the heuristic entries and of the certified entries sets how often and how far apart the anchors are, and whether the existing priors make a check easier or longer. The designer identifies the first as unresolved and is not certain about the second.
- **The measure of the benefit at a boundary** (§1, item 18; §2, item 13): the difference between the recorded proof size and the human parameter, which needs the human parameter (engine integration, §3).
- **The unit of the window of proof size** (§1, items 20, 21 and 29; §2, item 15): item 29 has the stop at the storage limit reached for the check. What is counted (nodes per claim, time, storage, or something else), whether it applies per claim or per iteration, and whether and how it differs from the storage window of the priors, given that item 21 describes it as a storage based limitation like Syzygy. Item 66 adds whether the window is applied to one boundary at a time, with the priors of a boundary loaded when required and replaced by those of the next, or to all the boundaries together. The designer does not claim which is best.
- **How engine evaluations affect the priors and the exact check** (§1, items 3, 40 and 59; §2, readings 54 and 55): left unintegrated by the designer on purpose, because of the uncertainty about the effect and the potential to optimize it. Whether item 3, which calls the engine the least impactful anchor, and item 59 are about different things, as reading 54 reads them. The designer states that the role of the evaluator is not exactly clear throughout (item 63).
- **Whether the optimal priors depend on the engine in the slot** (§1, items 60, 62, 63 and 64; §2, readings 56, 58 and 59): the designer states that this is likely and that it is not a definitional influence on the approach. Whether different engines are used against each other to optimize the priors is left open in general. Open: how any dependence is held, which engine the keep-or-gap decision is made for, whether the priors are held per engine, and how the depth that an evaluator reaches relates to the d of a check.
- **Whether the gap criterion of item 38 depends on the engine** (§1, items 38, 43 and 64; §2, readings 36 and 59): item 38 has the gaps reduced until an engine can parse them with certainty at an effective speed. Whether that makes the stored distance depend on the engine, and how that sits with item 64, is not stated. The designer states that this is an integrated specification concern, expected and understood (item 65).
- **What "provably better than the native line" means** (§1, item 61): whether it refers to the result or to the cost of resolving the gap, what proves it, given that an engine's finding is certified only once it is checked against every reply (calibration, reading 13), and what the native line is.
- **Whether the best result of an engine navigating a gap enters P** (§1, item 61): how the best result within the time limit relates to a batch and to the priors, and whether the time limit is the speed tolerance of item 34.
- **Master priors and the irreversible segmented boundaries** (§1, item 1): how the approach scales with them. This is scoped to P and is not a strict consideration at this point.

### 1.4 Statements made after reading 43 was written (construction file, items 48 to 71, titles only)

- Item 48: "Already covered" is closed in the designer's understanding.
- Item 49: The perimeter of the proof tree has zero d-depth.
- Item 50: Storing the perimeter.
- Item 51: Batches, the unknown region, and the order of composition.
- Item 52: The tree as a map of certification limits, and where its boundary is stored.
- Item 53: Keeping the boundary departs in principle from an earlier statement, and this is accepted.
- Item 54: How a later check proceeds is open to the programmable goal.
- Item 55: One basis is chosen, and the scope is not narrowed.
- Item 56: The root, the map, and the coverage of the d-depth limit as a goal.
- Item 57: Confirmations on position, objective and distance.
- Item 58: Reachability has a general sense, and a cost sense in the priors.
- Item 59: The engine as evaluator was deliberately left unintegrated.
- Item 60: The basis: priors that carry an engine of lower strength towards one of higher strength.
- Item 61: The evaluator navigates the gaps.
- Item 62: The priors can be iterated beyond the strongest engine.
- Item 63: The test cases are not required, and the general integration of the evaluator is what matters.
- Item 64: The engine dependency of the priors is not a definitional influence.
- Item 65: The goal that yields a proof tree, and the engine dependence of the gap criterion, are left to emerge.
- Item 66: The segmented dependency units and the storage window.
- Item 67: Repeated work can be composed on transpositions.
- Item 68: The evaluator may aid the working memory concern, and whether it works is not known.
- Item 69: The working memory is an open front for optimization more than program architecture.
- Item 70: The parameterization is read as modular, and exclusion is the concern.
- Item 71: The lifetime part is optional.

### 1.5 Items and readings cited by the copied blocks (construction file, titles only)

- Item 3: The engine is an evaluator, with the widest compatibility.
- Item 5: The goals are the programmable front of F.
- Item 6: The fronts of variance.
- Item 11: The approach to construction.
- Item 22: The exact check fits into F, and the loop it closes.
- Item 26: A reset adds to the scope, and the allowance that remains can be used.
- Item 28: The d-depth recorded for a board state is its certification limit, and the earlier specification does not apply.
- Item 29: The goal for the exact check is unresolved, and what depends on it.
- Item 30: What the check does to the priors, and why a batch is used.
- Item 31: Where the check starts: from the ends of the branches.
- Item 32: The benefit is at the edge of the exhaustive tables.
- Item 33: The layers from a certified region, and closing the gaps between them.
- Item 34: What decides keeping and omitting, and the balance of storage, coverage and speed.
- Item 36: The number of games is a goal parameter.
- Item 37: The balance in use: proof size is the distance, and storage is the limit.
- Item 38: Reachable as a proof tree is not the same as resolvable by the engine at effective speed.
- Item 40: The approach to the goal is left open on purpose.
- Item 41: Two kinds of goals, and the constant goal of P.
- Item 42: The cost of the check includes re-deriving the omitted board states, and this may reduce the need for differential testing.
- Item 43: Two applications, the frequency of crossing, and heuristic states.
- Item 44: The batch is checked as a sum within the window, and whether it is one operation is regulated.
- Item 45: Pruning at determined optimal branches, and what is not carried into the priors.
- Item 47: The inversion for the other side, and no hard rules.

- Reading 12: The check as a goal.
- Reading 17: The reset at a chaining point.
- Reading 38: Two kinds of goals, and two optimizations.
- Reading 42: Where the pruning of item 45 applies.
- Reading 47: How "already covered" may follow from items 48 to 51, and what is left.
- Reading 56: The engine and the keep-or-gap decision.

---

### 1.6 Statements made after item 71 (construction file, items 72 to 99, titles only)

- Item 72: The goal of P has a general preset, and it can be as variable as the goals of F.
- Item 73: The preset is found by running, from goals programmed by hand.
- Item 74: The empirical record is developmental, and the two functions are a signature with an open body.
- Item 75: For now the goal front is the body of the function, and the engines may differ by side.
- Item 76: The goal is a parameter that is the body, held as stored code.
- Item 77: The stored goal code is disjoint from the priors and the storage window.
- Item 78: The viable pool is enumerated in the check and not held in the priors, and d is the depth of verification to the forced solution.
- Item 79: The entry holds d and no move is added, and the metric of what d holds is what matters.
- Item 80: The parts of the d-depth metric an entry holds depend on its state of resolution.
- Item 81: What cannot be certified is not added, the index and the symmetry stay opaque until the legend and the transformations are specified, and the priors are segmented.
- Item 82: Three different boundaries, the legend as a prior, and what d denotes.
- Item 83: Bishop colour is part of a segment, and what the symmetry is for.
- Item 84: The d-depth architecture was chosen for storage, coverage and speed, and the header is compression.
- Item 85: Symmetry is compression, one entry covers its folds, and its notation is redundant storage.
- Item 86: A transformation into another segment is predictable and is exploited.
- Item 87: The d field carries status by itself.
- Item 88: Directed rederivation and the condition based search are used together.
- Item 89: The condition chosen from a pool of equals is relative to the scope, and the check of it is left open.
- Item 90: A segment has a finite set of conditions, and where the condition sits is not hard-coded.
- Item 91: Certainty is resolved back to front, and the search is bounded between the priors.
- Item 92: The engine evaluates and is not trusted first where priors exist.
- Item 93: Engine variation has to be accounted for, and it is open.
- Item 94: The priors judge the engine's picks afterwards, and putting them into the evaluator is not excluded.
- Item 95: A segment and its counterparts are a class.
- Item 96: The stop of the check is a variable of the goal front.
- Item 97: What the goal can ask of the priors cannot be answered before the goal exists.
- Item 98: The metric and the indices measure the certainty of optimal routing, and not win, draw or loss.
- Item 99: The goal is not built until the specification is understood.

## 2. Attempts

None have been made. This section is begun when attempts begin. An attempt is written in four parts: what is given, with its sources in the construction file; the interactions between statements that are not stated, put as questions to the designer and not decided; the goal logic, with a tag on each part naming what it depends on; and what could not be derived. The goal is stated as a specification and not as a procedure, with the open pieces as named parameters, since the construction file leaves the approach to the goal open on purpose (item 40) and the working memory as a parameter (reading 63). Failed attempts are recorded here and not in the construction file. If a result bears on a statement in the construction file, the designer decides whether the statement is revised.

## 3. Status

| Element | Status |
|---|---|
| §1 | Copies of the construction file of 22 November 2026, unchanged, and §1.6, added on 24 November 2026 with the titles of items 72 to 99 |
| §2 | No attempt made |
