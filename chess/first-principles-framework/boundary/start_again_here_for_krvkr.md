**Opening task (before building anything).** Read the programme documents and evaluate them. Report what is missing, what is unclear, and whether you understand how to proceed on every front (the primitives, the mining process, the dependency chain between segments, tempo, the exclusion rules). Where you don't understand a front, say so and ask. Do not build until the designer has answered.

**Goal.** Enumerate and progress the programme in any certifiable way, so that there is more structure for the Lakatosian programme to build on. Every step has to be considered against the programme's architecture and against scaling on all fronts: more materials, more segments, more of the space.

**Working frame.** For an irreversible segment (KRvKR is the first example), every position is a boundary board state *unless* a rule excludes it. The rules are written in primitive terms (coverage net, capture/mating net, segmented dependencies) and bound the forced region through ply-bounded temporal calculations, with closure. Depth extends only when the structure calls for it.

**Syzygy's role.** A measuring instrument for the branch to the win and for the landscape. It is not a target to reproduce, and reproducing win/draw/loss classes is not progress.

**What to measure.**
1. **Forced-region count.** The moves to mate along the Syzygy branch. Only the minimum matters for optimality: won in x, the opponent needs x+1 to equalize, or the extra move changes the outcome.
2. **Tempo flag.** One board state with a turn flag and no duplicate storage. Record whether the turn changes the forced-region status or the position is forced regardless.
3. **Perturbation.** How the equilibrium shifts when a move from a boundary state enters a forced region, and the primitive race inside it.
4. **Landscape.** The boundary states, the forced regions reachable from each, the value changes, and tempo derived from the changes in tempo between reachable states.

**Constraints.**
- A few small files. No storage blow-up, and no generated arrays in the repo.
- Run results never enter the programme documentation. The catalog and docs stay piece-agnostic, and specifics stay in the data.
- No search programs presented as the rule set.
- Before building anything, state what you will produce and why, and check it against this brief. Ask when a choice is the designer's.
- Use the names coverage net, capture/mating net and segmented dependencies. Use no numbered labels in discussion.
- Coherence is not evidence. Say what has been measured, and say what has not.

---

One caution for the new session: its first answer will rest on whatever documents it reads, so give it the programme documents and not the old `syzygy/` outputs. That way its evaluation comes from the programme itself.
