# Boundary Board State Primitive Catalog

**Status.** A specification, not a program. The catalog has exactly three primitives, and every entry is a race condition: a
conjunction of terms with an explicit interface, written so that it can be applied to any position without naming pieces and so
that it does not depend on any implementation. The three are **P-001**, the capture net, whose royal-target case is
the mating net (§4, written); **P-003**, chained dependencies, the segmented process of a chain of dependencies and the race to
resolve it (§7, SPEC); and **P-008**, the coverage race under disturbance (§8, SPEC). Nothing else is a
primitive (§1, the admission rule). This file is the primitives catalog only. Boundary states belong to a separate catalog
(candidate, proven and disproven; `boundary_board_states.md` §3), with the composition of primitives recorded on proven entries
(topology §8x). Evidence tags used throughout: **SPEC** (defined, not tested), **INSTANCES** (checked on specific positions),
**CERTIFIED** (backed by an exhaustive or complete case-list certificate; none yet).

Companion documents: `boundary_board_states.md` (the concept and methods; the definition of a global boundary state is in its
§14 items 7 and 8), `boundary_board_state_topology.md` (the record; P-001 was derived in its §8s, §8t and §8u; the reduction to
three primitives is recorded in its §8ab), `boundary_board_states_case_study_otb.md` (the game).

---

## 1. What a catalog entry is

An entry is the smallest reusable statement of how a position can be won, lost or held, written without naming pieces. There are
four levels, and only the first three are catalog content.

1. **Primitive.** The base definition: its terms, its inputs, its outputs.
2. **Mechanism.** One way of satisfying one term (attack side) or keeping it open (defence side), with a precondition and an
   effect. Geometry and piece types enter here, and only here.
3. **Derived form.** A named binding of one mechanism to each term. Derived forms are what get recognized in practice (the fool's
   mate is one).
4. **Instance.** A concrete position or line. Instances are evidence rows. They never carry definitional weight.

**Growth rule.** A new position adds an instance row. A position that needs a new way to satisfy a term adds a mechanism. A
position that combines known mechanisms in a new way adds a derived form. A new route to a decisive result adds a primitive, and
only through the admission rule below. Adding any of these rewrites nothing already in the catalog.

**Mechanisms are keyed by term.** A primitive that shares a term with another shares the mechanisms attached to it. Primitives
differ in which terms they have and in what they output, so a mechanism found in one primitive is available to every primitive
with that term. The merge of the mating net and the piece trap into one primitive (§4) is the first case: they share
every term and mechanism and differ only in the target kind.

**Admission rule: what a primitive is.** The catalog is the protective belt of the programme whose hard core is the boundary
board state, the drawn equilibrium sustained by a repetition cycle (`boundary_board_states.md` §14 items 7 and 8). An entry is
admitted as a primitive only if all five hold:

1. **A race condition.** It is a conjunction of terms. Each term is a condition on a position that one side can satisfy and the
   other can hold open, and each term has mechanisms. The entry decides who completes the conjunction first, or who holds it,
   under perturbation of the position by a move.
2. **Irreducible.** It is not the composition of admitted primitives. A missing primitive is a deduction from the residual
   interface of composed primitives (the reduction rule, topology §8x), and it is admitted only as that deduction. P-008 is the
   worked case: it was derived from P-001 and P-003.
3. **Not a description of a result.** A name for what happens when a primitive resolves is not a primitive. Neither is a name for
   the state that emerges when primitives are in equilibrium.
4. **Piece-agnostic and board-agnostic.** It states no rule of any piece and no feature of a particular board. Both enter as
   inputs when the primitive is applied (§6).
5. **Modular.** It does not overlap the terms of another primitive. An overlap is a defect and is fixed by moving the content,
   not by packing several dependencies into one entry (topology §8aa).

Content that fails the rule is not discarded. It is classified as one of: a **definition** (of the hard core), an **input** (to a
primitive), a **term** or **mechanism** (of an existing primitive), a **derived form** or an **instance**. It is recorded in the
table below, and it is never added to the index.

**Candidates that are not primitives.** Entries proposed as primitives that do not pass the admission rule, and what each one is:

| ID and name | What it is | Where its content lives |
|---|---|---|
| P-002 Race to a square | the tempo comparison. Race is what every primitive is, so it is a term and not a primitive. A pawn race is the promotion route of P-003 against obstruction | term TP of P-003 (§7.5), term TM of P-008 (§8.5), mechanism RT2 |
| P-004 Edge restriction | the board edge as missing edges in the movement relation. It has no race condition: it is an input | the navigation graph of P-003 (§7.3) |
| P-005 Repetition loop | the definition of a boundary state: the cycle sustained by both sides when the primitives are in equilibrium (a forced terminal draw with no choice to deviate is a second realization, `boundary_board_states.md` §14 item 10). It emerges from the primitives, so it is the hard core and not a primitive | §5 item 2; `boundary_board_states.md` §14 item 7 |
| P-006 Piece trap | merged into P-001: the capture net on any target, of which the mating net is the royal case | §4 |
| P-007 Activation | the result when a chain of enabling dependencies resolves (development, unlocking, castling, promotion). The race condition in it is the chain | P-003 (§7): terms RT, LK, OB, TP |
| SUB-1 Substrate | the movement relation and the dependency graph, the structure on which the chain is read. A part of the chained-dependencies primitive, not a separate entry | P-003 (§7.3, §7.4, §7.7) |
| P-003 Frozen structure | the chain read statically, without a side or a horizon | P-003 (§7.6, the structure record) |

**IDs are stable labels and not a count.** IDs are never reused (§6). P-002, P-004, P-005, P-006, P-007 and SUB-1 are not reassigned. The catalog has three entries, P-001, P-003 and P-008.

## 2. Entry template

Every entry has these fields ("none" is allowed, blank is not). This implements and extends the proposal in topology §8i.

| Field | Content |
|---|---|
| ID, Name, Version | IDs are never reused; the version increments when a definition changes |
| Scope | the route to a decisive result it covers, and what it says nothing about |
| Inputs | what must be supplied to evaluate it |
| Terms | the conjunction of conditions that define it |
| Outputs | its status vocabulary |
| Mechanisms | ways to satisfy or hold each term, with evidence |
| Derived forms | named bindings, with instances |
| Cost model | how tempo cost is attached to each mechanism |
| Footprint | the squares, lines and pieces that can affect it |
| Resources | what it consumes (pieces' coverage, tempo, safety) |
| Composition interface | what it requires, what it excludes, what it shares |
| Boundary-state role | what it reports for a candidate boundary state |
| Evidence | claim, status tag, where recorded |
| Open questions | what is not yet defined |

## 3. Index

| ID | Name | Covers | State | Source in the record |
|---|---|---|---|---|
| P-001 | Capture net (royal target: the mating net) | capture of a target piece by an attacker: checkmate for the king, material or tempo yield for any other piece | **written** (§4) | topology §8s, §8t, §8u, §8w |
| P-003 | Chained dependencies | the segmented process through a chain of dependencies and the race to resolve it: a pawn chain, an outpost, a lock held by duty, and the bringing of a unit or pawn into play (development, unlocking, castling, promotion). The movement relation of each unit supplies the routes and the dependency sets | **specified, SPEC** (§7, version 3; was P-007, SUB-1 and P-003 version 2) | topology §8a, §8d, §8f, §8y, §8z, §8ab |
| P-008 | Coverage race under disturbance | the effect of a perturbation on the held terms of every engaged primitive at once, the support those terms share, and the race to open one before it is restored, computed for both sides | **specified, SPEC** (§8, version 2; was §7c) | topology §8aa |

Three primitives are specified. P-001 is written, and P-003 and P-008 are SPEC and have not been run. What the three do not
account for is recorded as residual and is not pre-registered as a primitive: zugzwang without capture is not read by any of them.
By the reduction rule it is returned as a residual interface if a proven boundary-state entry needs it (topology §8x), and it is
admitted only as the deduction of that residual under the admission rule of §1.

---

## 4. P-001 Capture net

*The mating net and the capture net on any other piece are one primitive race condition (the designer's decision). The kind of
target is a parameter; the royal case is the mating net. The ID P-006 is not reused.*

### 4.1 Scope

**Covers:** the capture of a target piece by an attacker, by any combination of pieces:
checkmate when the target is the king, and the win of material, or tempo yield short of it, when
the target is any other piece. **Says nothing about:** decisive results by promotion, or by
zugzwang without capture. It names no piece types.

**Target kinds.** The primitive has one parameter, the kind of target: *royal* (the king) or
*non-royal* (any other piece). The kinds differ in two properties and in nothing else. (1) For
the royal target a check must be answered, so the legal-move rule already handles counter-threats;
for a non-royal target the defender may leave the attack unanswered, and that case is handled by
composition with counter-nets (§4.10). (2) Completing the net on the king ends the game, so the
outcome clause O is vacuous; for a non-royal target O is live (§4.3).

### 4.2 Inputs

A position P; the target T (a piece of the defender D, on square t) and its kind; the attacker A;
the side to move; and a horizon h in plies for evaluating positions that are not yet complete.

### 4.3 Definitions

- t is the square of T, and N(t) the set of squares T can legally move to. For the king these
  are the adjacent squares.
- Occ_D and Occ_A are the squares held by each side's pieces.
- Cov_A is the set of squares A attacks, computed with T lifted off the board so that lines
  passing through t are seen (the x-ray).
- Cov_D is the set of squares D's pieces can **legally** move or capture to, T's own moves
  included. A piece pinned to its king contributes only the squares on its pin line. This is the
  single place where a pin enters.
- The attackers are the A pieces attacking t (for the king, the checkers).

A reply by D **resolves** the net when, after it, T cannot be captured (royal: the king is no
longer in check) or the capture is unfavorable (non-royal). The **net holds** exactly when no reply
resolves it, which is when all five terms hold:

| Term | Condition |
|---|---|
| **C** Attack | T is attacked by A (royal: check) |
| **F** Escape | every square in N(t) is in Occ_D, or in Cov_A, or holds an A piece that is protected, so that T cannot capture it |
| **X** Capture | every attacker lies outside Cov_D: no D piece, T included, can legally capture it |
| **B** Block | every square between a sliding attacker and t lies outside Cov_D (vacuously true for an adjacent or knight attack) |
| **O** Outcome | royal: vacuous, since capture ends the game. Non-royal: the capture is favorable, meaning every D answer that leaves T on the board but makes the capture unfavorable (protection, or an exchange that equalizes) is closed. A recapture is a net on the capturing piece, so O is a chain of nested nets whose outcome is read from the composed race map. Static exchange evaluation is the stand-in until the catalog supplies that |

**Side condition S.** The sequence that builds the net must not leave A's own king in check. A
piece that shields A's king, or protects a net square, can therefore be part of the net.

For the royal target, C, F, X and B are the rules' own definition of checkmate, decomposed, which is
why the primitive applies to any position. The same conditions describe any capture, and the
extension to other targets adds only the outcome clause. Double check needs no special case,
because Cov_D is defined by legal moves and only king moves are legal in double check. The content
of the primitive lies in the mechanisms (§4.5) and their costs (§4.8).

### 4.4 Outputs

For an attacking move m, let **R(P, m)** be the replies of D in the position after m that resolve
the net, partitioned into:

- **R_F:** T moves to a safe square;
- **R_X:** an attacker is captured, by any D piece, T included;
- **R_B:** an interposition that saves T;
- **R_O:** T stays and the capture becomes unfavorable (empty for the king, since a king in check
  cannot stay in check).

R(P, m) is empty exactly when the net holds after m. The status of P for the target is one of:

| Status | Meaning |
|---|---|
| **COMPLETED** | T is already lost to the net (royal: D is checkmated), or A has a move m with R(P, m) empty |
| **ENGAGED** | A can attack T now or within h plies, but every such m leaves R non-empty. The report lists, for each m, the held terms (those with non-empty residuals) and the residual replies. For a non-royal target it also reports the **tempo yield**: the number of defender moves each residual reply would cost, and whether the reply can be a developing move |
| **DORMANT** | no A piece can reach an attacking square on t within h plies. The primitive does not apply |

Counter-threats are not a term: they are handled in composition (§4.10).

### 4.5 Mechanisms

Each row is a mechanism: how a term is satisfied on the attack side, with the
instances where it was observed. The instance codes are those of §4.6.

**C: check**

| ID | Mechanism | Instances |
|---|---|---|
| C1 | Open ray: a line from an A slider to k with every intermediate square empty, opened by the loser's own moves, by forcing, or by the attacker's own moves vacating it | I-1 (the f- and g-pawns left f2 and g3), I-2, I-4 (...Nf4 vacated h5 and exf4 vacated e3) |
| C2 | Adjacent or knight check: no line is needed | I-3, I-4 (Rh3) |

**F: flight**

| ID | Mechanism | Instances |
|---|---|---|
| F1 | Own occupancy: D's own piece stands on the flight square, by the undeveloped initial lock or by F3 | I-1 (d1, d2, e2, f1), I-2 (f8, f7), I-4 (f2, g2) |
| F2 | Attacker coverage: the square is in Cov_A, including through the king | all four |
| F3 | Decoy fill: a sacrifice or capture lands a D piece on the flight square | I-4 (...Nf4 exf4 puts a pawn on f4) |
| F4 | Protected attacker occupancy: an A piece stands on the square and is protected | I-4 (pawn g4 protected by f5) |
| F5 | Pinned target: T is absolutely pinned, so its legal moves are restricted to the pin line | I-5 (the queen on f3 pinned by the queen on d3) |

**X: capture of the checker**

| ID | Mechanism | Instances |
|---|---|---|
| X1 | No defender: no D piece attacks the checker's square | I-1 |
| X2 | Support: the checker is protected by an A piece | I-2 (Bd3 supports h7), I-3 (Rg1 supports g7), I-4 (g4 supports h3) |
| X3 | Masking: a D piece that would capture is blocked by a piece of its own, the king included | I-1 (the h1 rook behind the h2 pawn), I-4 (the queen on d3 behind the king on g3) |
| X4 | Pin: a D piece that would capture is absolutely pinned off the line | I-4 (the interposing queen on f3 is pinned); by the legal-move rule for an absolute pin |

**B: block**

| ID | Mechanism | Instances |
|---|---|---|
| B1 | Adjacency or knight check: nothing to block | I-2, I-3, I-4 |
| B2 | No interposer: no D piece can reach a between-square | I-1 |
| B3 | Pinned interposer: the only D pieces that could interpose are pinned off the line | by the legal-move rule for an absolute pin |

**S: attacker safety**

| ID | Mechanism | Instances |
|---|---|---|
| S1 | Attacker shield: a piece keeps A's own king out of check and also serves the net | I-4 (the f5 pawn shields Black's king from the queen on d3 and protects g4) |

**O: outcome (non-royal targets)**

| ID | Mechanism | Instances |
|---|---|---|
| O1 | Unprotected target: no recapture is available to D | none recorded yet |
| O2 | The equalizing recapture fails: the attacker is of lower worth than T in the composed race map, so the exchange favors A even if T is protected | I-5 (the g4 pawn attacking the queen), I-6 (a knight and two pawns attacking the queen) |
| O3 | Support chain: T is a node of a SUPPORT chain (§7), so a capture at T is answered by the next supporter and the exchange is read along the chain. The capture is favorable to A only if the cascade reaches an unsupported end of the chain or a node whose loss completes another net. The candidate replacement for static exchange evaluation where the structure is a chain | none recorded yet (SPEC) |

**Carriers.** How an enabling step is compelled rather than volunteered. These are
not terms; they say how a mechanism is reached.

| ID | Carrier | Instances |
|---|---|---|
| K0 | Cooperation: the defender makes the enabling move voluntarily | I-1 |
| K1 | Check: the enabling move is the only legal reply | I-3, I-4 |
| K2 | Attack on a piece, typically a queen: the move saves material | I-4 (...Nf4) |
| K3 | Sacrifice: the defender takes material and lands on the wrong square | I-3, I-4 |
| K4 | Pin: coverage is removed without any move by the defender | I-4 |

**Defence mirror.** The net is a conjunction, so the defence needs only one term to
stay false. This follows from the definition, with no further evidence needed:

| Term kept open | What that requires |
|---|---|
| D-C | every checking square for A unreachable, or the ray held by a piece that cannot be removed or pinned away |
| D-F | at least one flight square available and outside Cov_A: a luft, or no own-piece lock |
| D-X | an unpinned, unmasked D piece attacking each possible checker square, or the checker left unsupported |
| D-B | an unpinned interposer for each between-square |
| D-O | non-royal only: keep the capture unfavorable, by protecting T against an attacker that is not cheaper or keeping T out of reach of cheaper attackers |

The chess-standard names for these (luft, guarding the delivery square, and so on)
are not verified against any project data and are used only as labels.

### 4.6 Derived forms and instances

A derived form is a binding of one mechanism to each term (for a royal target the outcome clause O is vacuous and is omitted). The first column of the
instance table names the evidence.

| Derived form | Binding (C; F; X; B; O) |
|---|---|
| **DF-A** self-locked king, open ray, unguarded delivery | C1; F1 + F2; X1 + X3; B2 |
| **DF-B** locked king, open ray, supported delivery | C1; F1 + F2; X2; B1 |
| **DF-C** coverage net, supported delivery | C2; F2; X2; B1 |
| **DF-D** decoyed lock, masked and pinned defenders, protected delivery | C1 then C2; F1 + F2 + F3 + F4; X2 + X3; B1; S1 |
| **DF-E** pinned target attacked by a cheaper piece (non-royal) | C2; F5; X4; B1; O2 |

| Instance | Position or line | Form | Carriers | Status |
|---|---|---|---|---|
| I-1 | the fool's mate: 1.f3 e5 2.g4 Qh4# | DF-A | K0 | INSTANCES |
| I-2 | the OTB game after 9.Qe4, the threat Qxh7# (parried by ...f5) | DF-B | threat only | INSTANCES |
| I-3 | the first puzzle: Rxg7+ Kh8 Rh7+ Kxh7 Qg6+ Kh8 Qg7# | DF-C | K1, K3 | INSTANCES |
| I-4 | the second puzzle: Nf4 exf4 Rh6+ Kg3 Qd3+ Qf3 g4 Qxd3 Rh3# (the mate is reached only if White takes the queen) | DF-D | K1, K2, K3, K4 | INSTANCES |
| I-5 | the second puzzle after ...g4: the white queen on f3, pinned to its king by the queen on d3 and attacked by the g4 pawn. Residual replies Qe3 and Qxd3, and Qxd3 is steered into I-4 | DF-E (engaged, not completed) | K2, K4 | INSTANCES |
| I-6 | the queen hunt in the game record: ...Nb4, ...f5 and ...e4 each attack the queen, and 3 of White's 5 queen moves in the first 12 were made under attack | attacks by cheaper pieces (C, O2); escapes not closed | K2 | INSTANCES (tempo yield) |

### 4.7 Evaluation procedure

The procedure is stated for any position and refers to no piece types. It is a
specification; no implementation is required to apply it.

```
evaluate(P, T, D, h):
    t <- square of T ;  A <- the other side
    if T is already lost to the net in P (royal: D is checkmated):  return COMPLETED
    deliveries <- every legal move of A in P that attacks T (a check, if T is the king),
                  or, if D is to move, in P after each legal move of D (up to h plies)
    if deliveries is empty and no A piece can reach an attacking square on t within h plies:
        return DORMANT
    for each m in deliveries:
        Q <- P after m
        partition the replies of D in Q that resolve the net into
            R_F : T moves to a safe square
            R_X : an attacker is captured, by any D piece, T included
            R_B : an interposition that saves T
            R_O : T stays and the capture becomes unfavorable   (empty if T is the king)
        if R_F, R_X, R_B and R_O are all empty:  return COMPLETED (by m)
        record (m, R_F, R_X, R_B, R_O)
    for each residual reply r:  tag r with the mechanism of 4.5 whose establishment would remove it;
                                if T is non-royal, record the tempo yield of r
    return ENGAGED, with the held terms per delivery, the tags and the tempo yields
```
Counter-nets are applied in composition (§4.10), not here.

**Defence form.** For D to move, the *invariant-preserving* moves are those after
which, for every reply of A within h plies, the status is not COMPLETED. When the capture of the target
is the only route to a loss, these are the non-losing moves. At a boundary state
that set is the set of equilibrium moves: one element in the strict case, possibly several
otherwise (`boundary_board_states.md` §14 item 8).

### 4.8 Cost model (unfinished)

The cost of establishing a mechanism is the number of A moves it takes, with D
replying adversarially. No closed form exists for general positions, so the entry
carries cost only as "one or more, to be supplied by search". What is known: the
cost with both sides cooperating is the helpmate distance, the cost against best
resistance is the mate distance, and CMD ≤ DTM, with the gap measuring the
defender's resistance (topology §8r; the gap reading is this record's own and
untested).

### 4.9 Footprint and resources

**Footprint:** N(t); the attack lines from every A slider to t; the delivery squares; the pin lines through D pieces and through T; and the squares from which D pieces cover any of these.
**Resources:** D's coverage sets Cov_D, A's attacking pieces, the tempo budget, and A's own king safety (S).

### 4.10 Composition interface

- **Requires:** nothing. It applies to any position.
- **Overload:** one D piece cannot discharge two defence mirrors at once. If it is the only
  unpinned interposer and also the only guard of a delivery square, the attacker wins by forcing
  it to choose (topology §8f).
- **Shared tempo:** a completion needs a tempo budget that the tempo comparison (the TP term of P-003 and the TM term of P-008) may also
  claim. They must be accounted for jointly.
- **Nets on different targets.** For an attacking move m, let Res_1(m) and Res_2(m) be the
  replies of D that resolve each net. The nets are jointly compatible for the attacker when both
  are engaged and Res_1(m) ∩ Res_2(m) is empty, so that every reply leaves at least one net
  unresolved. This is the double threat, and it applies to any pair or set of targets, kings
  included.
- **Counter-nets (the composition rule for counter-threats).** The defender's counter-threats
  are the nets of this catalog applied with the colors swapped: a check is an engaged net on the
  attacker's king, a larger threat is another net. A net is completed in composition only if no
  counter-net resolves or outruns it, which is a comparison of tempo between nets in opposite
  directions. The recursion (a net against its counter-nets, each with counter-nets of its own)
  needs a well-founded ordering; tempo is the candidate. This rule depends on the cost model
  (§4.8) and on the tempo comparison, which is a term of P-003 (TP) and of P-008 (TM). For the royal target the
  legal-move rule already forces a check to be answered, so the rule is built in.
- **Escape steering.** A residual escape of one net can fall into a mechanism of another. In the
  second puzzle the queen's escape Qxd3 puts it behind its own king, which is the masking
  mechanism X3 of the net on the king (instance I-5 into I-4).
- **Chained dependencies, as structure (P-003)** supply F1 at no cost to the defender, and cannot be undone. An
  irreversible pawn move removes such occupancy permanently, so it is a one-way change to the F
  term.
- **Chained dependencies, as resolution (P-003)**: when a chain resolves, the mechanisms that nets use become available (gain), and the chain receives obstruction from nets (§9).
- **Chains (P-003, §7).** A SUPPORT chain turns the outcome clause O from a nested exchange into a structural reading: a capture at a node of a chain is answered by the supporter, and the attacker's route into the chain is a net on a critical node (the unsupported end, or a node with several dependents). Status SPEC; the reading is offered as the candidate to replace static exchange evaluation (§4.3, term O), and it has not been run.
- **Coverage race (P-008, §8).** A disturbance of a structure that changes the coverage a net depends on (F1, X2 or the support that closes O) is the way a net completes against a structure that held. P-008 reads P-001's nets, carriers K0 to K4 and outcome clause, and returns the status of the structure after the disturbance. Status SPEC.
- **Coverage is used once.** Each D piece appears in Cov_D once, so two defence invariants that
  both rely on the same piece are compatible only if they are consistent about how it is used.
- **Reports.** For each target, P-001 returns a status, the held terms, their mechanism tags and,
  for non-royal targets, the tempo yield. These are the quantities that other primitives and the
  boundary-state check consume.

### 4.11 Role in boundary-state evaluation

For a candidate boundary state, and for each target against which A has an engaged net (always
the two kings, and any other piece A can attack), P-001 must return **DORMANT** (the primitive
does not apply) or **ENGAGED** with at least one held term for every delivery (an explicit
equilibrium: this target survives because term t is held open by mechanism m). **COMPLETED**
against the king of the side not to move cannot occur in a drawn position. A net that completes
against a non-royal target can occur in a drawn position when composition with counter-nets
compensates (a perpetual, for example). A threat against the side to move is an ENGAGED state, held
by the equilibrium move(s). The loop moves of the equilibrium must preserve these held terms.
P-001 covers capture nets only: a boundary state must also be accounted for by every other
primitive that engages, and until those entries exist the accounting is incomplete.

### 4.12 Evidence

| Claim | Status | Where |
|---|---|---|
| The terms read checkmate on four unrelated mates, by reading the positions | INSTANCES (not machine-checked) | topology §8s, §8u |
| Non-royal targets: the queen in the second puzzle (I-5) and the queen hunt in the game record (I-6) | INSTANCES | topology §8w |
| Outcome clause O and the tempo yield, as defined | SPEC | not run |
| The chain reading of the outcome clause (mechanism O3) | SPEC | §7, topology §8z; not run |
| Non-terminal evaluation, §4.7 | SPEC | not run |
| Defence mirror, §4.5 | SPEC | follows from the conjunction |
| Cost model, §4.8 | not defined | |
| Use on a boundary-state candidate, §4.11 | not tested | no candidate has been evaluated with it |

### 4.13 Open questions

The cost model. The choice of horizon h. Whether the mechanism list is complete on a larger corpus
of captures. Whether the DORMANT/ENGAGED distinction gives a useful scoping signal on non-terminal
positions. How the derived forms relate to the named mating patterns (back-rank, smothered, the
Greek-gift mate). How a complete catalog derives the outcome clause O and the counter-net
comparison without outside evaluation: this needs the cost model (§4.8), the tempo comparison of P-003 and P-008,
the other primitives, so that a composed race map exists, and a well-founded ordering for the
recursion between a net and its counter-nets. The tempo ledger over a whole game (topology §8w). Whether the chain reading of O (§7) reproduces static exchange evaluation where that is exact and corrects it where a chain cuts across the exchange (topology §8z).

---

## 5. Evaluating a boundary-state candidate with the catalog

1. **Drawn, and an equilibrium:** each side has a deviation that loses against perfect
   play. The number of non-losing moves at the node is recorded as a width and is not
   required to be one. Established by the composition (a complete case list, topology §8d),
   with exhaustive search as insurance on the first entries.
2. **Cyclic realization.** The drawing branch closes into a loop whose moves are
   available to both sides, and both sides are live (the definition is in
   `boundary_board_states.md` §14 item 7). This is the definition of a boundary state, the
   emergent effect of the primitives being in equilibrium. It is not a primitive (§1). A forced terminal draw with no choice to
   deviate is a second realization of the same equilibrium (`boundary_board_states.md` §14 item 10).
   Established by search.
3. **Primitive accounting.** For each catalog primitive that engages, the status is
   DORMANT or an explicit set of held terms (P-001, §4.11). The accounting is recorded on the
   proven entry (topology §8x).
4. **The loop preserves the invariants.** The equilibrium moves at the node
   are the ones that keep them.
5. **Joint compatibility.** No resource is claimed by two invariants in conflict
   (§4.10).

Items 1 and 2 are established by the composition, and exhaustive search is insurance on the first
entries. Items 3 to 5 are the composition itself: they say *why* the position holds and scope the
search for neighbors. At present P-001 is written and P-003 and P-008 are specified and not run, so the
accounting is complete only along the capture nets, and it is as complete for chains and for
disturbances as their SPEC entries are.

## 6. Conventions

- IDs are never reused. Versions increment when a definition changes. Evidence rows
  are append-only.
- An entry is not tagged CERTIFIED without an exhaustive or complete case-list
  certificate.
- Rediscovery of known theory is recorded as such and is not counted as new
  structure.
- Instances carry no definitional weight. If an instance disagrees with a
  definition, the definition is revised and the version incremented.
- A heuristic may order candidates and may propose mechanisms. It certifies
  nothing.
- Discovery by residual, the reduction of a proven verdict that the equalized catalog does not
  account for, is not a heuristic. It is a deduction under the definition (topology §8y), and its
  output is the exact specification of what the missing or incompatible primitive must discharge.
- Primitives and substrates state no rule of any piece. The attributes of a unit's movement relation are inputs, and
  piece-specific special cases are accounted for during composition and race calculation on the board structure to which a
  primitive is applied. A special case recorded inside a primitive would make that primitive a subspace of piece
  composition. (Added with the movement relation, §7 and §10.3.)
- Every primitive is evaluated for both sides, with the colors swapped, and the candidate moves of each side perturb the other side's evaluation. This restates the counter-threat composition of §4.10 for all primitives; P-008 formalizes the coupling (§8.8). The effect of a move on a side's structure depends on the other side's options and on the other side's own coverage. (Added with P-008, §8 and §10.4; designer's statement, topology §8aa.)

---
- A primitive is admitted only under the admission rule of §1. Something that names the result of a primitive, a state that
  emerges from primitives, an input, or a term shared by several primitives is classified as such and is not indexed. (Added
  with the reduction to three primitives, §10.5; designer's statement, topology §8ab.)
- Perturbation has two settings (designer's statement, topology §8ap). In the legal-move setting, the change in primitive state caused
  by a move is what moves up a branch to find where a race was acted upon. In the construction setting, where a position is surveyed
  without regard to how it is reached and without the initial-position anchor, perturbation changes piece positions without the
  legal move restriction. Its scope is set by its purpose and by whether it is anchored. The definitions above, which state
  perturbation by a move, are unchanged. The boundary-state catalogs are one
  separate catalog, the boundary board state certification catalog (survey §4j), kept apart from this one.

---

## 7. P-003 Chained dependencies (specification)

**Status: SPEC.** Not run. Specified from the designer's statements recorded in topology §8y (the segmented process of
promotion, development and unlocking), §8f (locks as duties), §8z (chain structures and the movement rules of units) and §8ab (the
three primitives).

**One primitive.** Chained dependencies is one primitive, whose race condition is the segmented process through the chain (the
designer's statement). Activation (development, unlocking, castling, promotion) names the result when a chain resolves and is not a
primitive. The movement relation, the dependency graph and the static reading of the structure are parts of this entry (§7.3, §7.4, §7.6).

### 7.1 Scope

**Covers:** structures whose parts depend on one another, in sequence or in support, and the race to resolve them: a pawn chain, an
outpost, an immobile unit used as terrain, a duty lock (topology §8f), and the segmented process by which a unit or pawn is brought
into play (development from the start position, unlocking by moving its own blockers, castling, promotion), which is the resolution
of a chain of enabling dependencies. The movement rule of each unit supplies the routes and the dependency sets that the chain is
read on. **Says nothing about:** what a resolved chain then captures (P-001), or whether the coverage that a chain provides survives
a perturbation (P-008). It names no piece type and states no piece rule (§6).

### 7.2 Inputs

A position P; a chain, which is a sub-graph of the dependency graph G (§7.3); for each unit its movement relation M; the side S
that resolves the chain and the opponent O; a horizon h in plies. For the static reading (the structure record, §7.6) no side and no
horizon are needed. The values of a relation's attributes for each unit kind are inputs and are not recorded here (§7.7).

### 7.3 The structure: movement relation, navigation graph and dependency graph

**The movement relation.**

- **M(s, P)** is the set of squares a unit on s can legally move to in P, ignoring whose turn it is and ignoring pins (a pin enters
  once, in Cov_D of P-001, §4.3). **C(s, P)** is the set it can capture on. For most relations C = M, and where they differ the
  difference is an attribute (§7.7).
- **Navigation graph N(M, P):** the squares are the nodes, and there is an edge s → t exactly when t ∈ M(s, P). **N0(M)** is the
  same graph on the empty board.
- **Dependency set Dep(s, t):** the squares other than s and t whose occupancy decides whether the edge exists. It is empty when only
  the destination matters, and it is the set of between-squares when the unit travels along a line. This one attribute is how a wall
  or a chain cuts a route.
- **Edge cost:** one tempo, plus **exposure**. The arrival at t is exposed when t is in the opponent's coverage and a P-001 net on
  the arriving unit completes, or the exchange there is unfavorable. An edge whose exposure is infinite is **denied**.
- **Route distance d_M(s, t, P)** is the cheapest path, and the **route-length set L_M(s, t, P)** is the set of lengths of paths
  from s to t. The distance feeds the tempo term TP; the set matters where the parity of a route matters.

**Lock state and gain.**

- **Lock state of X:** its legal coverage Cov(X) is restricted by (a) its own blockers on its lines, own pieces and the king
  included; (b) being a pawn that has not reached the last rank; (c) obstruction by O, meaning squares denied by O's coverage or by a
  pin.
- **Gain G(X, c):** the set of catalog mechanisms X can supply once it stands in configuration c, for example a ray that serves C1,
  coverage that serves a defence mirror, a contribution to Cov_A. Gain is stated in mechanisms, not in squares or piece names.
- A **route** for X is a sequence of legal moves of S that brings X to a configuration c with a gain G: a path in N(M_X, P). Its
  tempo cost is the path length with each edge weighted by its exposure.

**The dependency graph G.** Nodes are units, squares and duties (a unit's coverage obligation, topology §8f). Edges are typed:

- **ENABLE (a → b):** b's route stage or coverage requires the removal, displacement or occupation of a. A blocker on Dep(s, t) is
  the source of an ENABLE edge to the route edge s → t. This is the lock removal of term LK.
- **SUPPORT (a → b):** the presence of a makes the capture of b unfavorable, by recapture or by coverage of the arrival square. This
  is P-001's X2 for a checker and its outcome clause O for a target.

Node attributes: the permanence of its coverage (§7.4), its mobility class (frozen, live, or duty-locked, topology §8f) and its
relation M. Defined quantities:

- **Root:** a node with no in-edge of the type. For SUPPORT it is the **unsupported end**, the weak end: a capture there is not
  answered.
- **Critical node:** a node whose capture or displacement makes a dependent fail one of its own terms. The set of critical nodes is
  the cut set of the structure.
- **Cascade:** the sequence of failures that follows the loss of a node along its out-edges, and its depth.
- **Completion cost of an ENABLE chain:** the sum of the tempo costs of clearing its nodes in order, each weighted by its exposure.

G adds no input. It represents what the mechanisms of P-001 already say: support and masking are SUPPORT and ENABLE edges.
**SUPPORT edges are read from coverage** (a recapture exists). The cascade is the iteration of that reading, and it terminates
because each step removes a node, so no circularity with the outcome clause of P-001 arises.

### 7.4 Invariants computed from the relation, and permanence of coverage

These are outputs of a breadth-first search on N0(M). None is an input, and none names a piece.

- **Orbit:** a class of squares reachable from one another in N0(M). A target configuration outside the orbit of the unit's square
  is unreachable by that unit alone, whatever the occupancy.
- **Grading:** a map g from squares to {0, 1} that every edge flips. If one exists, every route between s and t has a length
  congruent to g(s) + g(t) modulo 2, so the unit cannot adjust its route length by a single tempo.
- **Reversibility:** s → t implies t → s. A relation that is not reversible and has an acyclic N0 is **irreversible**: no unit can
  return to a square it has left.
- **Distance table and diameter on the empty board.**
- **Wall sensitivity:** whether the relation has non-empty dependency sets. A relation with empty dependency sets is not cut by
  occupancy along a line.

**Inclusion.** If the edges of one relation include those of another, then in the same position the route distance under the first
is no greater than under the second, for every pair of squares. This holds by definition, and it is a statement about route
distance only. It does not say that more edges are better for the side that owns them. Where the number of available moves matters
(stalemate, zugzwang, the parity of tempo in the race terms, the repetition of a boundary state's cycle), the extra edges enter
through composition and not here. Which relation includes which is an input supplied with each unit's relation, and it is not
recorded in this catalog.

**Permanence of coverage.** A square q is **permanently uncovered by side X**, relative to X's current units and with no relation
change, when for every unit u of X, q lies outside the attack closure of u: no square in the reach of u in N0(M_u) attacks q. Two
structural facts, both read from the relation, make coverage permanent or permanently absent:

- **Irreversibility:** N0 is acyclic, so the reach of the unit only shrinks as it moves. Coverage once lost is lost for good.
- **Orbit confinement:** the reach of the unit is its orbit, so squares outside it are never covered by that unit.

A unit with a reversible relation and a single orbit has only transient coverage: it can leave and return. Permanence is relative to
the unit set. A capture only removes coverage, and a relation change (promotion) or the arrival of a new unit adds it, so the
closure of an irreversible unit must include the region where its relation changes. That point is open (§7.14). This is the
structural content of "passed" and "outpost" (DF-P1 and DF-P2, §7.9), and of the weakness of the squares that lie outside the orbit
of every remaining unit of a side.

### 7.5 Terms

The chain resolves (the unit is in play, the pawn promotes, the lock is removed) when four terms hold:

| Term | Condition |
|---|---|
| **RT** Route | a route exists: a path in N(M_X, P) from X's square to a configuration with a gain |
| **LK** Lock removal | every own blocker on the route can be removed or avoided, and its removal is itself available: the ENABLE chain into each route edge has no node that is locked in turn |
| **OB** Obstruction-free | O has no capture net on X along the route (P-001), no SUPPORT chain of O that answers X's attempt to clear the route, and the route squares are not occupied or denied by O, permanent denial included (§7.4) |
| **TP** Tempo | the route's weighted cost fits the budget left after the competing races (the tempo comparison, which is the race term of each primitive: TP here, TM in P-008) |

The recursion in LK runs over a finite graph, so it terminates. A route that fails at RT because the target configuration lies
outside the unit's orbit is reported DORMANT before any tempo question is asked (RT3).

### 7.6 Outputs

**Status.**

| Status | Meaning |
|---|---|
| **COMPLETED** | the chain is resolved: the gain is realized, and the tempo spent is recorded |
| **ENGAGED** | the route is under way: the remaining tempo, and the residual obstructions, each a term held by a mechanism |
| **DORMANT** | no route within h plies, or the chain is inert terrain. A frozen structure with no legal moves is the degenerate case, held at no cost to the defender |

Both sides' chain states are reported, and their difference is the **activation asymmetry** of topology §8y.

**Structure record (static, with no side and no horizon).** For a chain S, P-003 reports from the graph G: its nodes and edges with
their types; its roots and unsupported ends; its critical nodes and cut set; the cascade closure of each node (the static reach of
its loss); whether it is **intact** (every edge's condition holds in P); and its **irreversible nodes** (those whose removal or
advance is a one-way change, §7.11). The structure record is what P-008 reads.

### 7.7 Attributes and perturbation classes

**Attributes of a relation.** A relation may specify the following, and the primitive reads them as inputs. Their values for each
unit kind are not recorded in this catalog (§6):

- **Reversibility,** with the irreversible case (§7.4).
- **Move set and capture set,** which may differ.
- **Relation change at a region:** the relation is replaced by another from that node on.
- **Kind of dependency set** (§7.3): empty, or the between-squares of a line.
- **State-dependent edges:** an edge available only while a state flag holds. An edge that couples two units is left open (§7.14).

**Perturbation classes.** A move changes the graphs in one of five ways. Each has a defined effect on RT, LK, OB and TP, and the sum
over the classes is the ΔS of Method 7 (topology §8y).

| Class | What changes | Effect |
|---|---|---|
| unit move | the unit's node in N; N itself is unchanged | route distance from the new square; under a grading, the parity of the remaining route flips |
| blocker move or capture | occupancy | edges of line-travelling relations whose Dep sets contain the square appear or disappear; leap edges change only at the destination |
| coverage change or pin | edge exposure | an edge becomes denied or is restored |
| relation change | the unit's relation | N is replaced from that node on |
| capture of a node of G | the node and its edges | cascade along SUPPORT and ENABLE edges |

These are also the perturbation classes that P-008 reads (§8.4).

### 7.8 Mechanisms

| ID | Mechanism | Instances |
|---|---|---|
| RT1 | Direct development: one move to an active square | I-A1 |
| RT2 | Promotion route: a pawn's march to the last rank, whose race is the tempo comparison TP against the opponent's route and obstruction | I-A3 |
| RT3 | Orbit and parity exclusion: the target configuration lies outside the orbit of the unit's square, or the route-length set cannot match the tempo class available (§7.4). The route fails with no tempo calculation | none recorded yet |
| RT4 | Relation inclusion: a unit whose relation includes another's has a route distance no greater than the other's in the same position (§7.4). Distance only: extra edges bear on stalemate, zugzwang and tempo parity through composition, not here | none recorded yet |
| RT5 | Wall bypass: a route edge with an empty dependency set is unaffected by occupancy along a line, so a structure standing on the dependency sets of other route edges does not cut it | none recorded yet |
| LK1 | Own blocker moves: a pawn move opens lines for a line-travelling unit | I-A1 (3.g3 and 4.e3 in the game) |
| LK2 | Castling: the rook is activated and the king leaves the lock | I-A1 (5...O-O) |
| LK3 | Enable chain: the lock is a chain of ENABLE edges, each node removable in turn; the completion cost is the sum along the chain with each step weighted by its exposure | none recorded yet |
| OB1 | Coverage obstruction: the activating square is attacked by a cheaper piece, so activating there loses X to a net | none recorded yet |
| OB2 | Occupation or blocking of the route | none recorded yet |
| OB3 | Pin: X is pinned, so its route is restricted to the pin line | none recorded yet |
| OB4 | Supported obstruction: the denying unit stands on a SUPPORT chain, so capturing it is answered (O3 of P-001) and the unit cannot clear its route by exchange | none recorded yet |
| OB5 | Permanent denial: a route square is covered by a unit whose coverage of that square is permanent (§7.4), by irreversibility or by orbit confinement. Its negation is what leaves a route free of denial for good | none recorded yet |
| TP1 | Advance with tempo: the move that advances the chain also attacks (carrier K2 of P-001), forcing a reply that cannot be a developing move | I-A1 (the knight and pawn attacks on the queen) |

**Defence mirror.** To keep a chain from resolving, O needs one term to stay false: deny the route, keep a blocker locked, keep an
obstruction (a net on X), or win the race.

### 7.9 Derived forms and instances

A derived form binds one mechanism to each term. Both forms below are SPEC with no instance recorded.

- **DF-P1 passed route.** RT2, with every square of the route permanently uncovered by the opponent's irreversible relations (the
  negation of OB5). The only obstruction left is a P-001 net by a unit with a cyclic relation, which is a race (TP) and not a
  structural fact.
- **DF-P2 supported outpost.** A configuration reached by RT1 whose square is permanently uncovered by the opponent's irreversible
  relations, and whose occupant carries a SUPPORT edge, so that capturing the occupant is answered (O3 of P-001) and an irreversible
  unit cannot dislodge it. The weak point is the unsupported end of the chain that supports it.

| Instance | Evidence | Status |
|---|---|---|
| I-A1 | the game record, first 12 moves: White has 3 activation events (Nc3, Qf3, Bd3) and 5 queen moves, 3 of them under attack; Black has 5 (Nc6, Bc5, Nf6, O-O, Qg5) | INSTANCES |
| I-A2 | the fool's mate: the lock is the defender's own undeveloped pieces, and a developing move such as Nf3 discharges terms of the net (the knight attacks h4 and occupies f3). The lock is recorded in topology §8s. The discharge is by geometry and was not run | INSTANCES (lock); not run (discharge) |
| I-A3 | king and pawn against king: resolution by RT2. The downstream net, a queen against a bare king, is completable whenever the new queen is not immediately lost and the position is not a stalemate, so the race to resolve the chain alone decides the result (standard endgame knowledge) | standard knowledge |

### 7.10 Cost model, footprint, resources

**Cost.** The route cost is the path length in the navigation graph weighted by edge exposure (§7.3), and the cost of an ENABLE chain
is its completion cost. It is a specification and has not been run. **Footprint:** the nodes and edges of the chain, the dependency
sets of the route edges, the squares from which supporters cover, and the nets of P-001 on the unit and on the chain's critical
nodes. **Resources:** the supporters' coverage, used once per unit (§4.10), the tempo budget, and permanence.

### 7.11 Composition interface

- **Requires:** P-001 for the obstruction term OB and for the exposure of an edge (a net on the arriving unit), and the tempo
  comparison that the race terms make.
- **Supplies:** to P-001, gain at resolution (the mechanisms newly available: coverage, rays, supporters, defenders), own occupancy
  (F1), support (X2) and the chain reading of the outcome clause (O3), the last two at no cost to the defender; to P-008, the
  structure record (§7.6), the support of each residual reply (the edges of the relation and their dependency sets) and the
  perturbation classes (§7.7).
- **Receives:** from P-001, obstruction: a net on the unit along its route or on a critical node of a chain, and the outcome of a
  promotion as a target; from P-008, the duty locks (an element in the support of two held terms is locked to both, topology §8f)
  that lock removal reads, and the effects of a perturbation that open ENABLE chains (LK3).
- **Obstruction is a net on the unit,** and the outcome of a promotion is itself a target. The mutual dependence between P-001 and
  this primitive is well-founded by tempo.
- **Irreversibility:** resolving or dissolving an irreversible node (an advance or a capture by a unit with an irreversible
  relation) is a one-way change to the F term of P-001 (§4.10), so such a node is a cut point in time as well as in the graph.
- **Duty trade-off:** a unit that holds a duty in a defence mirror (a guard, an interposer) leaves it when it advances elsewhere.
  This is the overload of topology §8f, in the graph a node that serves two edges it cannot serve at once (coverage is used once,
  §4.10). The single move that discharges several duties at no cost, recorded in the OTB case study §7 and topology §8u, reads as a
  move that adds or satisfies two edges at once. This is a reading offered for the designer's decision on whether the pattern is a
  named one. It is not asserted.
- **A frozen structure is the chain read statically, and unlocking is the same chain read as a process.** They are one primitive.

### 7.12 Role in boundary-state evaluation

For a candidate boundary state, each chain on the board is reported for both sides with its status and held terms. A chain that is
DORMANT is inert terrain: the frozen-structure fortress certificate of topology §8d is the case list for it. A chain whose
resolution is available but loses (resolving it opens a held term of P-008, or lets a net of P-001 complete) is ENGAGED with the held
term named: this is the duty lock of topology §8f, and it is an explicit equilibrium. The loop moves of the equilibrium must keep
these held terms (§5). The activation asymmetry (§7.6) is reported with them.

### 7.13 Evidence and predictions, stated before they are run

| Claim | Status | Where |
|---|---|---|
| The four terms, the mechanism table and the structure record | SPEC | not run |
| Orbit, grading, dependency set, inclusion and permanence, as definitions | SPEC | §7.3, §7.4; not run |
| The joint dependencies that G is to reproduce | to be established by ablation | not run |
| The dependency graph reproduces those joint dependencies | not run | prediction P3-1 |
| Orbit exclusion (RT3) reproduces the known orbit-confined verdicts | not run | prediction P3-2 |
| The chain reading of the outcome clause (O3 of P-001) | SPEC | §4.5, topology §8z; not run |
| The game record's activation events (I-A1) | INSTANCES | §7.9 |
| Base-of-chain attack as P-001 with P-003 | mapping only | §10.2 |

1. **P3-1. Cut sets and load-bearing sets coincide.** In an ablation of a position (topology §8s), the sets of pieces whose removal changes the
   verdict are the cut sets of the graph G read from the same position by the mechanisms of P-001. A position where the minimal
   load-bearing set is not a cut set of G refutes the construction of G.
2. **P3-2. Orbit exclusion reproduces the known orbit-confined verdicts.** For the recorded structures whose draw rests on a unit
   confined to an orbit that excludes the square it needs (the color-lock entries of topology §8b, and the color-lock fact that the
   record found inside Uri's rule, topology §8n), RT3 reports the unit's route onto the needed square as DORMANT, and the verdict is
   a draw, across the whole class and not only the textbook instance. These structures are instances and carry no definitional
   weight (§6). A won position in the class whose win needs that route, with no other route, refutes the prediction. The
   identification of the color-lock fact in §8n with orbit exclusion is this file's reading and is to be checked against §8n.

### 7.14 Open questions

The cost model for routes, specified in §7.10 and not run. The unit of the gain measure. Whether castling is a state-dependent edge
of two units or a separate lock removal. How the chain states of the two sides combine with the states of the other primitives
(topology §8y). How the closure of an irreversible unit includes the region of a relation change, which permanence needs to be
sound. How pins and x-rays enter the dependency sets, since P-001 lifts the target to see through it. Whether the weights of edge
exposure can be defined without the composed race map that P-001's outcome clause also needs. Whether a chain whose critical node is
a cut node (several dependents) needs a separate term. How the cascade is bounded where SUPPORT edges form a cycle (two units
protecting each other). Whether the chain reading of O reproduces static exchange evaluation where that is exact and corrects it
where a chain cuts across the exchange (topology §8z).

---

## 8. P-008 Coverage race under disturbance (specification)

**Status: SPEC. Version 2.** Not run. (Renumbered from §7c with the reduction to three primitives, §10.5; no definition changed, and the entries it read as P-007 and SUB-1 are now P-003.) Derived by the reduction rule of topology §8x from the composition of P-001 and P-003, as the
designer set out (topology §8aa). It is a race-condition primitive built from the held terms, residual replies and mechanisms that
the written entries already define, and it adds no evaluation of its own: every effect it reports is determined by the definitions
of the primitive whose term is affected.

**Basis.** The primitive is built from the primitives' own race conditions under perturbation (the designer's statement), not from a
formalization imported from outside the race semantics.

**How the calculation is read.** P-008 is a primitive race condition and is evaluated as one. Its condition is calculated exactly on a position, and again on the position after each perturbation (Method 7, topology §8y). Certainty about a line comes from those evaluations composed over a complete case list. Certainty about the catalog grows as mechanisms are added and as residual interfaces are discharged (topology §8y). Where this entry is silent about an effect, that is a gap in the catalog to be developed by the reduction rule, and it is not an uncertainty in the calculation.

**The residual.** P-001 reads coverage at one position, as sets. P-003 reads support as present or absent. Composed, a net on a
supporter or coverer is an exchange, and an exchange, or any other forcing move, changes the position: it removes or relocates
coverers, changes occupancy and can change exposure. Neither primitive says whether the coverage a structure needs still
suffices afterwards, or whether the side that disturbs it can make the change before the other side repairs it. The elements that
their reports read can be read by more than one held term (§4.10: coverage is used once), so a perturbation changes several held
terms at once. That interface is what this entry fixes. It is computed for both sides (§8.8).

### 8.1 Scope

**Covers:** whether the coverage a structure requires still holds after a disturbance; whether a disturbance is compelled,
feasible and timely; for both sides. **Says nothing about:** the identity and static properties of the structure (P-003), the net
on any single target (P-001), or the route and tempo comparison itself, which it reads (P-003, terms RT and TP). It names no piece
type and states no piece rule (§6).

### 8.2 Inputs

A position P; for each side, the reports of the engaged primitives (their status, held terms, residual replies with their
mechanism tags, and the structure records of P-003); the movement relations (P-003, §7.3); a horizon h in plies; the side to move.

### 8.3 Held terms and support

- **Held term.** A term with a non-empty set of residual replies (§4.4): the replies of the holding side that keep the term from
  being satisfied against it. The edges of a P-003 structure and its terms RT, LK, OB and TP are held terms of the primitive in the same
  sense. The **register** of a side is its set of held terms with their residual replies.
- **Support of a residual reply.** The elements whose state the existence and efficacy of the reply depends on: the unit that
  makes it, the squares it uses, the edges of its movement relation with their dependency sets (P-003, §7.3), and the occupancy those
  edges read. The support of a held term is the union over its residual replies.
- Two statements follow from this definition and need no evidence. (a) A perturbation whose change misses the support of every
  residual reply of a held term leaves the term held. (b) A held term is opened only by a perturbation whose change meets the
  support of each of its residual replies.
- **Opening sets and key elements.** An opening set of a held term is a minimal set of elements whose change removes every
  residual reply. It is found from the supports of the replies, which are few (P-001 records 2 and 3 replies at the key points of
  the second puzzle, topology §8u), and its effect is read from the owning primitive's definition on the changed position. A
  singleton opening set is a **key element**, and a **key defender** when the element is a unit.
- **Joint support.** An element in the supports of residual replies of two or more held terms, of the same or of different
  primitives. It is the content of "coverage is used once" (§4.10) and of overload (topology §8f): a reply that spends the element
  for one term removes it from the support of the other, and no special case is needed.
- **The condition is evaluated exactly at a position, and again on the position after every perturbation (§8.9).** Recaptures,
  discovered coverage and changed occupancy are perturbations whose effects are recomputed at the position each one produces. They
  are not corrections applied afterwards. Whether a definition in this entry is right is a question about the entry. It is answered
  by comparison with proven verdicts and, where they disagree, by the residual interface and the reduction rule (topology §8y). An
  entry found wrong is corrected and its version incremented (§6).

### 8.4 Disturbances and effect vectors

- A **disturbance** is a perturbation of the position by a move or a compelled sequence. Its **delta** is the set of elements it
  changes, classified by the perturbation classes of §7.7: a unit move; a blocker move or capture; a coverage change or a pin; a
  relation change; the capture of a node of G.
- Its **effect vector** gives, for every entry of the registers of both sides, one of: **UNCHANGED** (the delta misses the support,
  exact by (a)); **HELD** (touched, and still held after the owning primitive re-evaluates the term, possibly with other residual
  replies); **OPENED** (touched, and the owning primitive's definition on the new position establishes one of the term's
  mechanisms). A primitive that was DORMANT and becomes engaged is reported as **ENGAGED-NEW**; it is reached when the delta meets
  its engagement support, the units whose route within h (P-003, RT) reaches an engagement square.
- The effect of a class on a mechanism is defined where the mechanism is defined: the dependency sets of P-003 say which classes
  remove an edge, and the mechanisms of P-001 and P-003 say what an established mechanism does to a term. Where a class's effect
  on a term is not determined by a catalog entry, the reduction rule returns the residual, which is the specification of the
  missing mechanism or effect.
- An exchange is a sequence of perturbations by alternating sides, each with its effect vector. A recapture is the other side's
  perturbation whose change re-establishes a residual reply, and it spends the recapturing element, which is joint support. No
  special case is needed.
- The **carrier** is how a perturbation is compelled, read from K0 to K4 of P-001. A capture is one carrier among a check, an
  attack, a sacrifice and a pin. The primitive is defined by the effect and not by the capture: a piece offered so that its capture
  breaks a structure is a K3 whose delta combines several classes. Whether the carrier is a term of this entry is open (§8.14).

### 8.5 Terms

The perturber Y opens a held term h of X when all four hold:

| Term | Condition |
|---|---|
| **OP** Open | after the perturbation, h has no residual reply: the owning primitive's definition establishes one of its mechanisms on the new position |
| **CR** Compelled | each step of the perturbation is compelled by a carrier: every reply of X that avoids the delta loses elsewhere (the partition of §4.7) |
| **FS** Feasible | the effect vector of the same perturbation over Y's own register contains no OPENED entry: Y does not open its own held terms in making it |
| **TM** Timely | the perturbation completes before X's perturbation re-establishes a residual reply of h or opens a held term of Y, compared by tempo (P-003, TP) |

### 8.6 Outputs

For each entry of each register: its opening sets and key elements, and the joint support it takes part in. Status: **COMPLETED**
(all four terms hold, with the sequence), **ENGAGED** (some perturbation meets an opening set and some term is held, with the held
terms, their mechanisms and the residual replies listed), or **DORMANT** (no perturbation within h changes any element of an
opening set). The report also carries the support closure of §8.9.

### 8.7 Mechanisms

Mechanisms are the ways a perturbation opens a held term, by the class of its delta.

| ID | Mechanism | Instances |
|---|---|---|
| OP1 | Removal of a supporting element: a unit in an opening set is captured, displaced or drawn from its post, and no other residual reply remains (classes: unit move, capture of a node). The classical removal of the defender; rediscovery of known theory is recorded as such (§6) | none recorded yet |
| OP2 | Joint spend: a perturbation, or a reply to it, changes an element in the support of two held terms. Resolving or restoring one term removes the element from the other, which opens. Overload (topology §8f, §4.10) is the case in which the element is the only member of a support of each | none recorded yet |
| OP3 | Cut: occupancy changes on the dependency set of an edge in the support, so the edge is removed (class: blocker move or capture) | none recorded yet |
| OP4 | Denial: the post or recapture square of a supporting element becomes exposed or denied, or the element is pinned, so its reply no longer exists (class: coverage change or pin; Cov_D of P-001) | none recorded yet |
| OP5 | Opposing entry: a unit comes into the condition on the opposing side (a new attacker, a new line), or a relation changes, so a residual reply no longer resolves (classes: unit move, relation change) | none recorded yet |
| OP6 | Forced-in arrival: a piece is compelled by a carrier to capture or to occupy a square. Its delta combines OP1 (it leaves its post), OP3 (its arrival occupancy) and OP4 (its arrival exposure). The designer's case of a capture whose capturer is incompatible with the structure's race | none recorded yet |
| CR | the carriers K0 to K4 of P-001, not repeated here | I-4 and I-5 of P-001 for the carriers; no P-008 instance |
| FS1 | Disjoint support: the elements the perturbation changes are not in the support of any held term of Y's register | none recorded yet |
| FS2 | Compensation: the perturbation opens a held term of Y's, and Y's perturbation that opens a held term of X wins on tempo (the counter-nets of §4.10) | none recorded yet |
| TM1 | No restoring element can reach a square of the support in the tempo remaining (the route distance of P-003, RT and TP) | none recorded yet |
| TM2 | The restoration costs a tempo that the carrier denies by forcing the reply | none recorded yet |
| TM3 | The restoring perturbation itself opens a held term of X's (it spends joint support) | none recorded yet |

**Defence mirror.** The held term is opened only if all four terms hold, so the defence needs one to stay false: D-OP (every opening
set has a member that no perturbation of the opponent can change), D-CR (no carrier compels the delta), D-FS (the perturbation
opens the perturber's own held term) or D-TM (the restoration is in time).

**Illustration, not an instance (§6).** A unit s makes a residual reply for the edge that supports a node a of a chain, and also
makes the residual reply by which a net on a target T is held (it can capture an attacker). Its support is in both entries, so it
is joint support. A perturbation that exchanges s for an attacker spends s, and the effect vector shows the edge into a OPENED even
though no unit of the chain was attacked. Whether that wins is then CR, FS and TM: whether the exchange is compelled, whether Y's
own register stays held, and whether X can restore a residual reply in time.

### 8.8 Both sides

The primitive is evaluated for both sides with the colors swapped. The effect vector of a perturbation is computed over the
registers of both sides at once, so the feasibility term FS is read off the same vector: the effect on Y's own entries. For X's
register the perturber is Y, whose candidate moves inside the support closure of X's register (§8.9) generate the perturbations;
for Y's register the roles reverse. So the removal of a defender is not a property of one side's structure alone. It depends on the
other side's options, and those depend on the other side's own held terms. The recursion between a perturbation and its
counter-perturbations uses the well-founded ordering that §4.10 names for nets against counter-nets (tempo, with each
perturbation spending a unit or a tempo). This entry adds no ordering of its own. The counter-net composition of §4.10 is the
special case in which the compared outcomes are nets.

### 8.9 Cost model: primitive evaluation under perturbation, with complete coverage

1. **The calculation is exact at every position.** For a position, the registers of both sides are read from the reports of the
   engaged primitives, and the supports, opening sets, key elements and joint support are read from those reports by the
   definitions of §8.3. A perturbation is a sequence of moves. After each move the engaged primitives are evaluated again on the
   resulting position, the effect vector is read, and the terms OP, CR, FS and TM are decided from those evaluations. This is
   Method 7: ΔS(P, m) is the change in the primitive state under a move. A recapture, a discovered line and a change of occupancy
   are perturbations of the classes of §7.7, recomputed at the position they produce. Nothing is estimated, and nothing is
   corrected afterwards.
2. **What has to be covered is set by the mechanism tables.** The perturbation classes enumerate the ways a move can change the
   state of a structure. A perturbation with a null effect on every engaged primitive in the footprint leaves the state unchanged,
   and that fact covers it without branching. A perturbation has a null effect on a condition exactly when its change misses the
   support of that condition, the elements the condition reads, as defined in §8.3. This follows from the definition of the
   condition and needs no evidence. The perturbations that must be branched are those whose change meets a support.
3. **An exact pruning result.** By (b) of §8.3, a held term is opened only by a perturbation whose change meets the support of every
   one of its residual replies. A line on which no reachable change can meet the support of some residual reply is therefore
   decided DORMANT by derivation under the definitions, without enumerating it.
4. **Exhaustive coverage over what remains.** The **support closure** is the least set of elements that contains the supports of
   all entries of both registers and the engagement supports of the DORMANT primitives, and is closed under the way a perturbation
   reaches an element: the units that can change an element of the closure (the nets of P-001 on it, and the routes of P-003 to a square of a dependency set of an edge in the closure) are added with their own supports. It is a fixed point over the
   definitions. A perturbation whose change misses it cannot change any entry. Over the perturbations whose change meets it, the
   coverage is complete: a case list in which each case is decided by primitive calculation. That is the proof form of the project
   (topology §8x), with exhaustive search as insurance on the first entries, and it is the scoping that topology §8h states.
5. **What bounds it.** The footprint of topology §8h, with its four conditions (forcing moves, closure under lines, kings as the
   global resource, inertness proven), and the completeness of the effect enumeration. Both are claims of the catalog, and a gap in
   either is a residual interface that the reduction rule develops (topology §8y). Neither is an error rate. Where the footprint
   grows to the board, the number of perturbations to branch grows with it. Additional primitives that decide a term without
   enumeration, or that shrink the footprint, reduce it further. That is kept open (§8.14).

### 8.10 Composition interface

- **Requires:** the reports of the engaged primitives (P-001, P-003), including the supports, dependency sets and perturbation
  classes that P-003 defines, and the tempo comparison, which is a term of P-003 (TP) and of this entry (TM).
- **Supplies:** to P-001, the effect of nested nets under perturbation, which is part of the composed race map that the outcome
  clause names as missing (§4.3, §4.13). That map also needs the tempo comparison (TP of P-003) and the cost model of §4.8, so O is not yet fully derived. To
  P-003, the duty locks (joint support) that lock removal reads and the effects that open ENABLE chains (LK3).
- **Resources:** coverage used once per unit (§4.10), which is the joint-support accounting, and tempo. **Footprint:** the support
  closure of §8.9.

### 8.11 Role in boundary-state evaluation

For a candidate boundary state, the registers of both sides are reported with their opening sets, key elements and joint support.
The equilibrium moves are those whose effect vector contains no OPENED entry on the mover's own register, so every held term stays
held (§4.11). Each losing deviation is a perturbation whose effect vector opens an entry of the deviator's register, by a named
mechanism OP1 to OP6, and that mapping is recorded on the proven entry as its deviation map (topology §8x). A deviation with no
refutation among the known mechanisms is a residual interface for the reduction rule.

### 8.12 Derived forms and instances

| Derived form | Binding (OP; CR; FS; TM) |
|---|---|
| **DF-H1** base-of-chain attack | OP1 on the support of an unsupported end; K1 to K4; FS1; TM1 |
| **DF-H2** poisoned offer (topology §8f) | OP2 or OP6; K3; FS1; TM1 |
| **DF-H3** forced-in arrival | OP6; K3; FS1; TM2 |

No instance is recorded. Removal of the defender is named in classical theory and is a rediscovery (§6). The poisoned offer is the
case recorded in topology §8f, and it has not been checked.

### 8.13 Evidence and predictions, stated before they are run

| Claim | Status | Where |
|---|---|---|
| The two-sided coupling is well-founded | argued, not proven; §8.8 adds no ordering of its own | §8.8 |
| The primitive's evaluation under perturbation, composed over a line, reproduces the exhaustive verdict | not run | prediction P8-2 |
| A perturbation whose change misses the support of a condition leaves that condition unchanged | by definition | §8.3, §8.9 |
| A held term is opened only by a perturbation whose change meets the support of every residual reply | by definition | §8.3 |
| Held terms, support, opening sets, key elements, joint support and the support closure, as defined | SPEC | not run |

1. **P8-1.** In an ablation of a position (topology §8s), the elements whose removal changes the verdict are the key elements of the held
   terms, that is, the singleton opening sets. Where a held term has no singleton opening set, the smallest removal set that changes
   the verdict is a minimal opening set. A smaller verdict-changing set refutes the support definition.
2. **P8-2.** The evaluation of P-008 on the position after each perturbation, composed over a line, reproduces the verdict of
   exhaustive search for exchange-dependent structures. A disagreement is a residual interface, which is a missing mechanism, an
   incompatibility with another primitive, or a wrong entry (topology §8y), and it develops the catalog. The number of
   perturbations that must be branched, against the number that are legal in the footprint, measures the saving.
3. **P8-3.** In ablation and deviation data on a composition, every perturbation that changes the verdict of a composition has a
   change that meets the support closure. A verdict-changing perturbation outside the closure refutes it: a support is missing, which
   means that a mechanism or an effect is absent from the catalog, and by the reduction rule it specifies what must be added.


### 8.14 Open questions

How the effect of nested nets under perturbation combines with the tempo comparison (TP of P-003) and the cost model to give the composed race map of the
outcome clause O (§4.3, §4.13); until then P-001 evaluates O as it defines it. Whether the carrier is a term of this entry or an
input it reads. Whether additional primitives reduce the support closure; the designer keeps this open as part of the design.
How the repetition cycle that defines a boundary state (§5 item 2) keeps every entry of both registers held across
the cycle. How large the support closure is on real compositions, which is recorded with each verdict and measured by build step 20
of `boundary_board_states.md` §11. Where a perturbation class acts on a P-003 term that P-003's mechanisms leave undetermined, the
reduction rule returns the residual.

**Recaptures and discovered coverage.** Recaptures and discovered coverage are perturbations of the classes of §7.7 and are recomputed at each position
(§8.9). What remains open about them is whether the five perturbation classes are complete for this primitive. That is a claim of
the catalog's completeness (topology §3, §8p), and a gap in it is a residual interface for the reduction rule.

---

## 9. Primitive dependency map (provisional)

The relations among the three primitives, stated at the level of primitives. Cycles are resolved by tempo.

| From | To | What flows |
|---|---|---|
| P-003 Chained dependencies | P-001 Capture net | gain at resolution: the mechanisms newly available (coverage, rays, supporters, defenders); own occupancy (F1) and support (X2) at no cost to the defender; the chain reading of the outcome clause (O3) |
| P-001 Capture net | P-003 Chained dependencies | obstruction: a net on the unit along its route, or on a critical node of a chain, which is the attacker's route into the cascade; the outcome of a promotion is a target |
| P-001 Capture net | P-008 Coverage race | the held terms and residual replies, the carriers K0 to K4, and the nets by which a perturbation reaches an element |
| P-003 Chained dependencies | P-008 Coverage race | the structure record: its edges as held terms, its critical nodes and mobility classes; the support of each residual reply (edges and dependency sets); the effect classes of a perturbation; the tempo comparison of the TM term |
| P-008 Coverage race | P-001 Capture net | an opened held term (F1, X2, or the support that closes O) lets a net complete; the effect of nested nets that the outcome clause's composed race map needs (§4.3) |
| P-008 Coverage race | P-003 Chained dependencies | duty locks (an element in the support of two held terms) and the effects that open ENABLE chains |
| P-008 Coverage race (side X) | P-008 Coverage race (side Y) | the feasibility term, read from the same effect vector over both sides' registers |

P-001 is the terminal primitive: every plan ends in a capture net. Race to a square, edge restriction and activation are not nodes
of this map: the first is a term (TP, TM), the second is an input (missing edges of the movement relation), and the third is the
result of a chain resolving (§1).

---

## 10. Provenance, and the mapping from the original catalog families

### 10.1 What is the designer's and what is this file's organization

The designer's: the architecture (a primitives catalog; candidate, proven and disproven
boundary-state catalogs with verdicts cached permanently; composition of primitives recorded on proven
entries; the reduction rule), the definition of a boundary state, the requirement that primitives be race
conditions stated without piece names and without board features, and the three primitives: the mating and
capture net on any target, chained dependencies (the segmented process through a chain, which includes
development, unlocking, castling and promotion, with locks and the movement rules of units), and the coverage race
derived from the two (topology §8ab).

**The origin of the names (the designer's statement).** The structure of the primitives was built by the designer first, and the formal chess
names ("mating net", "decoy", "opposition" and the like) were applied afterwards as labels that map the structure onto formal chess
understanding. The catalog was not constructed from a list of formal chess theory. That the structure and the formal theory agree is
read by the designer as convergence: two different descriptions of the same structure. The names in this file were mostly supplied by the
assistant when the designer described the structure, and so the vocabulary can look borrowed from formal theory when the structure was not.
The designer states this as fact about the origin, and not as a supplement to any proof.

This file's organizational choices, offered for the designer to change: the ID scheme (P-001 and on); the
four levels of primitive, mechanism, derived form and instance; the fields of the entry template; the term
codes (C, F, X, B, O for P-001; RT, LK, OB, TP for P-003; OP, CR, FS, TM for P-008); the status words COMPLETED,
ENGAGED and DORMANT; the defence mirror; the evaluation procedure; the evidence tags SPEC, INSTANCES and
CERTIFIED; and the wording of the admission rule of §1.

### 10.2 Mapping from the original families (`boundary_board_state_topology.md` §8b)

Section 8b describes the primitives catalog as known race-condition and zugzwang structures, each
recorded as a small ending, the inequality that decides it, and the pieces and squares it locks. In the
generalized catalog those entries are the concrete layer beneath the primitives: a small ending is a
derived form or an instance, the deciding inequality is the cost model of the primitive for that instance
(the tempo comparison), and the pieces and squares it locks are the footprint and the lock states. The
mapping below is this file's reading, and no entry has been written.

| Family in §8b | Generalized primitive or primitives |
|---|---|
| Edge and corner draws (rook-pawn king and pawn against king; wrong-colored bishop with a rook pawn; knight with a rook pawn on the seventh) | P-003: the board edge as missing edges of the movement relation (an input), the defender's king to the corner as a route and tempo comparison (RT, TP), and, for the wrong-colored bishop, orbit exclusion (RT3). The cycle that holds the draw is the repetition of the boundary state (§5 item 2), which is not a primitive. The decomposition of Uri's rule in topology §8n is a route race composed with orbit exclusion |
| Pure race: the rule of the square | P-003: the promotion route RT2 against the opponent's route and obstruction, compared by the tempo term TP |
| The opposition family: direct, distant and diagonal opposition; key squares; triangulation; corresponding squares | P-003, as the tempo and parity comparison between routes (TP, and the route-length set of §7.4), with the repetition cycle of §5 item 2 where the opposition is held by shuffling |
| Spare-tempo parity | P-003, the route-length set and the tempo term TP |
| Base-of-chain attack | P-001 with a non-royal target (the weak base pawn) and P-003 (the chain); the disturbance of the chain is P-008 (DF-H1) |
| Deflection and breakthrough: the outside passed pawn; a pawn sacrifice that creates a passer | P-003 (the passer's promotion route, DF-P1); the sacrifice is a decoy (mechanism F3 of P-001) |
| Blockade fortresses (opposite-colored bishops and similar) | P-003 (the locked structure, read statically as DORMANT) with the repetition cycle of §5 item 2 |
| Haworth's enumeration of mutual zugzwangs | not a primitive: raw material for candidate and proven entries, to be re-read against the definition of a boundary state (`boundary_board_states.md` §9) |

The expectation recorded in §8b that the edge, opposition and spare-tempo entries are the narrow ones, and
that fortress entries are wide, was stated in terms of the exactly-one reading. Under the
definition the width of a position is a recorded parameter, and the expectation is not a classification
(`boundary_board_states.md` §14 item 8).

### 10.3 Provenance of the movement relation and of P-003

**The designer's:** that chain structures are a key structure of the locking mechanism, as the dependency of
one element on another (the passed pawn as a chain of the pawn wall with the pawn passing unopposed; the outpost as a chain
in which taking the occupant allows a recapture, with the chained end as the weak point); that the capture potential inside a
chain alters the viability of the race condition it belongs to, so that it informs the segmented-dependency primitive
directly; that the movement rules of every unit shape routing and its joint compatibility with the chain structures (the
color symmetry of a diagonal-line unit, the parity and routing of the leaping unit, the line restriction of the rook-type
unit and its dependence on occupancy and coverage, the queen-type unit as the least restricted navigator other than the leaper,
and its routing moves that shorten a race); that all of this must generalize to a primitive race-condition structure at the
granular level while staying modular and piece- and board-agnostic; and that piece-specific rules are allowed only as
specializations of that general structure.

**This file's organizational choices, offered for the designer to change:** the formalism of the movement relation (navigation graph, dependency set, orbit,
grading, union, edge exposure); the two edge types ENABLE and SUPPORT; the notion of permanence and its two sources; the
three terms and the mechanism codes of P-003; the mechanism codes RT3 to RT5, LK3, OB4, OB5 and O3; the five perturbation
classes; and the two predictions of §7.13.

**Piece-specific rules (the designer's statement).** The movement relation is one construction, and the catalog states no rule of any
piece. The piece-specific special cases are inputs supplied with each unit's relation when a primitive is applied to a position, and
they are accounted for during composition and race calculation on that board structure. Recording them inside a primitive would turn
the primitive into a subspace of piece composition and defeat its purpose.

**A possible modular addition, not part of the protective belt.** A register of piece-specific special cases, kept apart from
the primitives and consulted by them as an input during composition, is viable by the designer's statement. It is not
started here, and it would hold no definition of a primitive.

### 10.4 Provenance of P-008

**The designer's:** that two written primitives, the capture and mating net and the chained structures with their
segmented dependencies, are enough to derive a third under joint compatibility, as a demonstration of the reduction rule and of
the complete graph of primitives; that the third is the coverage race of a chained structure after trades and the removal of key
defenders, where a structure that is covered can lose its coverage when a covering piece is traded and the race on the trade is
then won; that this is a graph-related coverage compositional problem viewed through race-condition calculation and not
exhaustive combinatorial coverage; that overlap between primitives compromises modularity and must be fixed, because packing
several dependencies into one primitive excludes emergent possibilities; that the primitives are composed adversarially, with
the same calculations performed from the other side and the perturbation of that side's candidate moves informing the structure,
so that the removal of a defender depends on the opponent's options; that a disturbance need not be a capture, since a piece
placed to threaten the structure whose capture leaves the capturer incompatible with the chained structure's race is the same
race; and that where combinatorics reaches a limit, a hybrid that scopes exhaustive coverage is viable provided the composition
needs only minimal combinatorial calculation, with further primitives possibly mitigating it more, which is kept open; and that the approach is race-condition primitives and not a heuristic approximation, so that what matters is making moves under perturbation and understanding the effects on the primitive race conditions under perturbation, with certainty derived from the growing catalog's account of those effects, and with the hybrid being exhaustive coverage formed by that understanding.

**This file's organizational choices, offered for the designer to change:** the name and ID of P-008 and the placement of its
section as §8; the register of held terms with the support of their residual replies; the effect vector over both sides'
registers; the opening sets, key elements and joint support; the four terms OP, CR, FS and TM and the mechanism codes; the
support closure and the evaluation under perturbation with complete coverage of §8.9; the treatment of the carrier as an input;
the placement of the dissolution race in P-008 and not in P-003, with the codes DF-H1 to DF-H3; the conventions bullet on
adversarial evaluation in §6; and the three predictions P8-1 to P8-3.

**Principles of construction (the designer's statements).** A primitive carries one behavior. Overlap between primitives is a defect,
so the dissolution race of a structure is in P-008 and not in P-003. A race condition is evaluated exactly at each position, and again
after each perturbation, and the coverage of perturbations is complete over the effects that the mechanism tables enumerate; a
disagreement with a proven verdict is a residual interface for the reduction rule (topology §8y). The primitive architecture is built from
the primitives' own race conditions under perturbation: the held terms and residual replies that the written entries define, the support
those replies read, the effect of a perturbation on them, and the race for them. No formalization is imported from outside the race
semantics. P-003 is the single chained-dependencies primitive.

### 10.5 The reduction to three primitives

Only P-001, chained dependencies and P-008 are race conditions under the definition of the programme. The designer's statement is that
there are three primitives: the mating and capture net, chained dependencies, and the coverage race. The other candidates are
descriptors: a name for the state that emerges when primitives are in equilibrium (P-005), a name for the result of a chain resolving
(P-007), an input with no race condition (P-004), the tempo comparison that every race contains (P-002), and an ID merged into P-001
(P-006). Their classification is in the table of §1. The admission rule of §1 tests every later addition.

Catalog sections cited in `boundary_board_state_topology.md` §8 before §8ab use the section numbers below. This is a lookup for those
citations.

| Cited section | Section now |
|---|---|
| §7 P-007 Activation | §7 P-003 (terms RT, LK, OB, TP; mechanisms RT1 to RT5, LK1 to LK3, OB1 to OB5, TP1; derived forms and instances) |
| §7a SUB-1 | §7.3 (movement relation, dependency graph), §7.4 (invariants, permanence), §7.7 (attributes, perturbation classes), §7.13 (predictions) |
| §7b P-003 | §7.6 (structure record) |
| §7c P-008 | §8 |
| §8 dependency map | §9 |
| §9 provenance | §10 |

**This file's reading, offered for the designer to change:** the choice of P-003 as the ID of the chained-dependencies entry; the reading that the tempo comparison is a term of P-003 and of P-008 and not an entry; the reading that zugzwang
without capture is a residual and not a fourth primitive.
