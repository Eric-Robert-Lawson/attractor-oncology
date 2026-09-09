# Loser Node Count & the Multi-Piece-Black Refactor

A reference document preserving the reasoning behind the "loser node count" metric and the concrete plan for extending the engine so Black can have real pieces (and can win), rather than being structurally forced into a bare king. Written to be re-entered cold, without needing to replay the conversation that produced it.

---

## Part 1 — Where this started: the existing two-stage optimality criterion

Already built, verified, and in production use (KQvK, KRvK, KBNvK):

1. **Primary criterion**: minimize moves-to-mate (matches Syzygy DTZ exactly — independently verified).
2. **Secondary criterion**, among distance-tied candidates: White minimizes / Black maximizes a *cumulative* count of the defender's own legal-move count, summed across every defender-to-move position along the forced line.

This distinguishes "dominant" mates (same speed, but more confining throughout) from "sloppy" ones (same speed, but leaving the defender more options along the way) — something Syzygy-style DTZ tables cannot express at all, since they only ever need to return *a* fastest move, never rank among ties.

**Empirically established, not just theorized:**
- KQvK: 396 unique positions where two or more genuinely independent (non-symmetric, non-transposing) paths tie on *both* distance and this metric simultaneously, out of ~182,000 positions checked.
- KRvK: independently verified against Syzygy DTZ.
- This is a real, non-trivial empirical result — not something guaranteed to exist by definition. Distance and escape-count are both smallish integers; the interesting open question (not yet investigated) is whether the rate of these ties exceeds what plain combinatorial density alone would predict, and whether they cluster geometrically (king distance, board region, ply-depth) rather than appearing uniformly at random.

---

## Part 2 — How the metric itself was refined: "escape count" → "loser node count"

### 2a. The problem that motivated the refinement

Original definition: sum, at every defender-to-move node along the line, the *total* legal-move count across all of the defender's pieces.

This breaks the moment the defender has more than a king. Concrete counterexample: a defender's king castled behind its own pawns has reduced mobility *because of its own pieces sitting in the way*, not because of any geometric pressure from the attacker. Counting that reduced number as if it reflected "closure" imposed by the attacker is measuring the wrong thing — it's an artifact of the defender's own piece placement, unrelated to how forcing the attack actually is.

### 2b. The resolution: scope mobility to the piece actually chosen

Don't sum legal moves across every piece the defender has. For a given defender-to-move node, each *candidate reply* already and unambiguously determines which specific piece it moves — that's known immediately, with no lookahead required. Compute the contribution *per candidate*, using only the legal-move count of the piece *that candidate* moves, not the defender's total mobility across everything on the board.

This isn't a new, separate calculation bolted on afterward — it slots directly into the mechanism that already exists: every child's contribution is already computed independently, per candidate, in both the White and Black branches of `classify()`. Nothing about this requires restructuring the existing per-child computation loop.

### 2c. Confirmed non-circular

Concern raised and resolved: does knowing "which piece the defender's perfect play moves" require already knowing the answer you're trying to compute? No — each candidate's own move already fixes which piece is involved, independent of which candidate ultimately wins the comparison. And by the time a parent node needs a child's value, that child is *already* fully classified — the fixed-point process always resolves children before parents use them. No circularity.

### 2d. The "fixed point" question — resolved, not by a new rule

Raised: at what point along a forced line does this mobility-counting actually start being meaningful — does it need a rule like "first position in check" or "first position with a unique defensive reply"?

**Resolved: no new detection rule is needed.** A defender-to-move position only ever enters `classified` (as a forced loss) if *every single one* of the defender's legal replies also leads to a forced loss — that's the literal definition `classify()`'s Black branch already enforces (`all_candidates_resolved` gate). Therefore every defender-to-move position that ever contributes to the cumulative count, from the first move of any classified forced-loss line, already satisfies "the defender is choosing among options that all lead to a forced loss." The existing classification boundary *is* the fixed point, by construction — nothing new needs to be built to detect it.

### 2e. Handling ties where candidates move different pieces

Raised: if two distance-tied candidates move *different* pieces (say, a queen move and a king move both achieve the same mate distance), is "the piece perfect play moves" ambiguous?

**Resolved: no contradiction, because the calculation is inherently per-candidate, not per-node.** Each candidate's contribution (mobility of the piece *it* moves, plus its own downstream cumulative count) is computed independently. The defender then picks whichever candidate produces the extremal (max, for the defender) total — exactly the existing maximize-among-ties mechanism, just with mobility now scoped per-piece instead of summed across the board. If multiple candidates that move *different* pieces genuinely tie on both distance and this refined total, that's a legitimate multi-path tie — consistent with how every other tie in this project has already been treated (a real finding to catalogue, not an error to resolve away).

### 2f. Rename: "loser node count"

Once the defending side can genuinely have material capable of winning (see Part 3), "Black escape count" stops being the right name for this quantity — it isn't about Black specifically, and it isn't really about "escaping" once the defending side might instead be actively trying to force its own win. The more accurate name: **loser node count** — the cumulative, per-piece-scoped legal-move count of whichever side ends up losing, counted only for the specific piece chosen at each of that side's decision points, from the first move of a classified forced-loss line onward.

### 2g. Why none of this touches anything already computed

With a bare-king defender (every material solved so far — KQvK, KRvK, KBNvK), "the piece the candidate moves" is *always* the king — there is no other piece it could be. The piece-scoped definition and the original definition produce identical numbers in every case where the defender has no material besides a king. **This refinement is a complete no-op on all existing results.** Nothing needs to be recomputed. It only becomes a real, non-trivial distinction the moment the defending side has a second piece.

---

## Part 3 — What multi-piece Black actually requires: verified against the real code, not assumed

Prompted by the question "isn't this just remapping White's logic to Black, since piece movement mechanics don't change?" — investigated directly against `compositional_trajectory_solver_modular.cpp`.

### 3a. What's already genuinely generic (confirmed by reading the code)

- **`generate_pseudo_legal_moves`**: already determines the mover generically (`Color mover = (gs.to_move == 'W') ? Color::WHITE : Color::BLACK;`) and generates moves for whichever side is actually to move. Never hardcoded to White.
- **`is_in_check(gs, Color c)`**: takes a color parameter directly — already symmetric.
- **`is_checkmate_general(gs)`**: checks "is `mover` in check with no legal moves," where `mover` is derived generically from `gs.to_move`. Already correctly symmetric — would correctly detect *either* side being mated, in isolation.
- **`GeneralState`**: just a piece list with per-piece colors. Nothing in the data representation itself assumes Black is restricted to a king.

On piece movement and check-legality mechanics specifically, the "just remapping" framing is correct.

### 3b. Where it actually breaks — exact locations, not a general concern

**Location 1 — `discover()`, checkmate-terminal handling:**
```cpp
if (is_checkmate_general(st)) { node.flags |= 1; nodes[key] = node; continue; }
```
This flag unconditionally means "distance 0, forced win." Safe today *only* because a bare-king Black can never physically deliver check — so every checkmate detected here is necessarily White mating Black. The instant Black has an attacking piece, this breaks: a position where Black's piece has just checkmated *White's* king would be flagged identically to White mating Black — both "distance 0, forced win" — silently misrepresenting a White loss as a White win.

**Location 2 — `discover()`, material-sufficiency check:**
```cpp
for (auto& p : st.pieces) if (p.color == Color::WHITE && p.kind != PieceKind::KING) white_pieces.push_back(p);
if (!material_can_mate(white_pieces)) { node.flags |= 2; nodes[key] = node; continue; }
```
Only White's material is ever evaluated for mating potential. Black's material, if any, is never asked the same question — "can Black's own pieces ever force mate against White" is currently not a question the engine can even pose.

### 3c. The deeper issue these two locations point to: the outcome space itself

The whole classification model is currently **binary**: every position is either a forced White win (some finite distance) or a proven draw (`classified` vs. `proven_draws`). `classify()`'s Black branch requires *all* children to resolve and takes the *max* distance — a structure that only makes sense if Black's sole objective is delaying White's win for as long as possible.

That assumption is true by necessity when Black is a bare king — there is no other option. It stops being true the moment Black could instead actively try to *win*. With real material on both sides, a position isn't just "White wins in N" or "drawn" — it's potentially **White wins**, **Black wins**, or **drawn**, a genuine three-outcome space. And Black's objective at a given node stops being uniformly "maximize White's distance" — it becomes "choose whichever is better: delay indefinitely, or force my own mate first." That is a different minimax shape, not a sign-flip on the existing one.

### 3d. Overall assessment of scope

The retrograde fixed-point *methodology* itself remains sound — this is confirmed, not in question. What needs to change is the **outcome data model** (binary → ternary) and the **two specific spots above** where "only Black can ever lose" was baked in as a structural assumption, not a deliberate design choice. That is a genuine, scoped refactor of the classification and terminal-detection layers — not a rebuild of the underlying search/verification approach, and not a sign that the architecture chosen so far was wrong for the problem it was originally solving.

### 3e. A concrete constraint now confirmed for the eventual design: the piece-count bit budget

Separately from this refactor, the engine's White-side piece limit was raised from 2 to 5 non-king pieces (7 total on the board, matching Syzygy's own convention) during this project. That work directly confirmed a fact worth having on hand when Part 4 item 1 actually gets designed: the packed 64-bit state key has room for **exactly 5 non-king pieces total**, and that number is stable whether those 5 slots all belong to White (today's scope) or end up distributed flexibly between both sides.

The arithmetic: today's packing uses 9 bits/slot (3 kind + 6 square) for every piece including both kings, with no per-piece color bit — color is implicit in which fixed "role" a slot represents (always White for the two non-king slots, since Black is assumed bare-king). That scheme has exactly enough room for 5 non-king slots and no more. A future scheme that lets *either* side occupy a given slot needs an explicit color bit per slot (10 bits instead of 9), but can recover the difference by dropping the currently-wasted kind bits on the two king slots (a king's kind is already implied by its fixed position in the key, so those bits are spent for no reason today) — 2 kings × 6 bits (square only) + 1 turn bit + 5 pieces × 10 bits = 63 bits, one bit to spare. So: **the piece-count ceiling itself won't need to move when this refactor lands** — only how the existing 5-slot budget gets allocated between sides changes, not the size of the budget. Worth confirming this arithmetic again once the actual redesign is underway, rather than assuming it still holds after other changes, but it's a solid, checked starting point rather than an open question.

---

## Part 4 — Concrete next steps (not yet resolved, in priority order)

1. **Design the three-outcome classification model.** Decide how `RetroResult` (or its replacement) represents "White wins in N" vs. "Black wins in N" vs. "proven draw," and how `classify()`'s White/Black branches need to change so each side is evaluating "which of my available outcomes is best for me," not just "minimize/maximize a single distance-to-White-win number."

2. **Fix the two identified locations** (checkmate-terminal handling, material-sufficiency check) to be genuinely bilateral — asking "who got mated" and "can *each* side's material ever force mate," not just White's.

3. **Consider a narrower first target** before the full three-outcome model: material where Black has a second piece, but that piece is provably incapable of ever delivering mate (e.g., a lone extra pawn with no realistic promotion/mating path against a king that can always contain it). This would let the piece-scoped loser-count metric (Part 2) get exercised on genuinely multi-piece Black, without yet requiring the harder three-outcome redesign — a real intermediate milestone, not a detour.

4. **Re-audit move generation and check detection under real load**, not just by inspection — confirm `generate_pseudo_legal_moves` and `is_in_check` behave correctly once Black's pieces are actually exercised in anger (pins, discovered checks, and interactions between Black's own pieces should already be handled by the existing "try the move, check if the mover's own king is left in check" mechanism, but this has never been tested against a Black side with real material before).

5. **Once (1)-(2) land**, the loser-node-count refinement from Part 2 becomes immediately meaningful (not just theoretically motivated) — worth re-running the escape-count computation on a genuinely multi-piece-Black material as the first real test of the piece-scoping logic.

---

## One-paragraph summary, for a very fast re-entry

The escape-count metric was refined to count only the legal moves of the specific piece the losing side's perfect play actually moves at each decision point (not all of that side's pieces), renamed "loser node count" to reflect that the losing side need not be Black once Black can have real material — this refinement is a verified no-op on everything already computed, since a bare king has no other piece to disambiguate. Extending the engine so Black can have real pieces (and can win) is a genuine, scoped refactor — not a simple remap — because two specific places in `discover()` currently hard-code "only Black can ever be mated" (the checkmate-terminal flag, and the White-only material-sufficiency check), and the whole classification model needs to expand from binary (White-wins-or-draw) to a genuine three-outcome space (White wins / Black wins / draw). The underlying retrograde fixed-point methodology itself is not in question — it's sound, and doesn't need to be replaced, only extended. One piece of concrete groundwork is already banked, separately from this refactor: the engine's packed-state bit budget was confirmed (while raising White's own piece limit to 5) to hold exactly 5 non-king pieces total regardless of how they're eventually distributed between sides, so this refactor's design work won't need to revisit the piece-count ceiling itself — see 3e.
