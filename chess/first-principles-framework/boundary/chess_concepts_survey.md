# Formal Chess Concepts: A Survey for Articulation

**Purpose.** A reference list of established chess concepts and their standard meanings, for the user's own use in articulating ideas. It is a vocabulary, not part of the boundary-board-state programme. Nothing here is admitted into the catalog. The last section lists *touchpoints* by name only and is marked as the assistant's reading.

**Conventions.**
- Definitions are the standard textbook/FIDE/engine-community usage. Where a term has several established meanings, they are numbered.
- "Ply" = one move by one side. "Move" in the traditional sense = a pair of plies.
- Square names are algebraic. "Side" = White or Black.
- Brief by design: each entry says what the term names, not how to play it.

---

## 0. Index

1. Rules-level and formal definitions
2. Time: tempo, initiative, development
3. Forcing moves and move classification
4. Tactical motifs
5. Mating patterns
6. Material, exchange and sacrifice
7. Pawn structure
8. Piece placement and piece relationships
9. Squares, lines, and space
10. King safety and king play
11. Strategic concepts and principles
12. Opening concepts
13. Endgame concepts
14. Zugzwang and related
15. Draw mechanisms and defensive concepts
16. Game-theoretic and solving concepts
17. Engine and analysis concepts
18. Error and judgement vocabulary
19. Terms with multiple meanings (disambiguation)
20. Touchpoints with the programme (assistant's reading, not admitted)

---

## 1. Rules-level and formal definitions

- **Legal move.** A move that obeys the piece's movement and does not leave the mover's own king in check. *Pseudo-legal* = obeys movement but may leave own king in check.
- **Check.** King attacked. Three ways to answer: move the king, capture the checking piece, interpose (not against knight or pawn checks, nor against double check where only a king move works).
- **Double check.** Two pieces give check at once. Only a king move answers it.
- **Checkmate.** Check with no legal reply. Ends the game.
- **Stalemate.** Side to move has no legal move and is not in check. Draw.
- **Insufficient material.** No sequence of legal moves can produce checkmate (e.g. K v K, K+B v K, K+N v K). FIDE also has a *dead position* rule below.
- **Dead position.** A position from which no sequence of legal moves can lead to checkmate by any series of legal moves (FIDE Art. 5.2.2). Includes some blocked-pawn positions. The game is drawn immediately.
- **Threefold repetition.** The same position, same side to move, same castling and en-passant rights, occurring three times. A claimable draw (FIDE 9.2).
- **Fivefold repetition.** Five occurrences. Automatic draw (FIDE 9.6).
- **Fifty-move rule.** 50 moves by each side with no pawn move and no capture. Claimable draw. **Seventy-five-move rule** is the automatic version. (Standard rule; the programme does not follow it.)
- **Draw by agreement.** Players agree.
- **Castling.** King moves two squares toward a rook, rook jumps over. Conditions: neither piece has moved, squares between empty, king not in check, does not cross or land on an attacked square. *Castling rights* are part of position identity.
- **En passant.** A pawn that has just advanced two squares may be captured as if it had advanced one, on the immediately following move only.
- **Promotion.** A pawn reaching the last rank becomes Q, R, B or N of its own colour. **Underpromotion** = to anything but a queen.
- **Touch-move.** Over-the-board rule; irrelevant to position analysis.
- **Position identity.** Board + side to move + castling rights + en-passant availability. (Engines hash this; repetition counts it.)
- **Forced move.** Only one legal move (e.g. a single answer to check). Not the same as *only good move*.

---

## 2. Time: tempo, initiative, development

- **Tempo** (several meanings, see §19):
  1. A unit of time = one move of one side.
  2. *Gaining a tempo* = achieving a result while forcing the opponent to spend a move responding, or reaching a position in fewer moves than otherwise.
  3. *Losing a tempo* = reaching a position in more moves than needed (moving a piece twice, a pointless back-and-forth).
  4. *Tempo move* = a move whose purpose is to pass the obligation to move to the opponent (§14).
  5. *Reserve tempo / spare tempo* = a harmless move available that does not weaken one's position; matters in zugzwang fights.
  6. *Tempo count* (pawn endings, races): number of moves needed to promote or to reach a key square.
  7. *Tempo bonus* (engines): a small evaluation credit for the side to move.
- **Initiative.** The ability to make threats that the opponent must answer, so that one side dictates the course of play. Distinct from *attack*: initiative is about who is forced to respond; attack is about targeting the king or a weakness.
- **Development.** Bringing pieces from the starting squares into play. **Lead in development** = more pieces developed than the opponent at the same moment.
- **Gambit.** A material offer (usually a pawn) in return for tempo, development, space or initiative. *Accepted / declined.* **Counter-gambit** = the defender returns a gambit.
- **Compensation.** Non-material advantages (time, activity, king safety, structure) offsetting a material deficit.
- **Dynamic factors** (time-dependent: tempo, initiative, activity) vs **static factors** (slow-changing: pawn structure, material, weak squares).
- **Zeitnot.** Shortage of clock time. A practical, not positional, concept.
- **Tempo-up / tempo-down structures.** Same position reached with different side to move (see *triangulation*, *Reti manoeuvre*).
- **Transposition.** Reaching the same position by a different move order. Position-level identity; matters for caching/trees.

---

## 3. Forcing moves and move classification

- **Forcing move.** A move leaving the opponent few or one reasonable reply. The classic order of forcing moves: **checks, captures, threats** (CCT).
- **Check / capture / threat.** Check attacks the king. Capture removes material. Threat creates an unstated next-move gain.
- **Quiet move.** Neither check nor capture, and not an immediate threat; often prepares a later forcing line.
- **Threat.** A move (or position) that would win something or mate on the next move if unanswered. **Mate threat**, **material threat**, **positional threat**.
- **Prophylactic move.** Prevents the opponent's plan or threat before it appears (Nimzowitsch, Petrosian).
- **Waiting move.** A neutral move that keeps the position intact while the opponent has to commit.
- **Developing move.** Brings a piece into play.
- **Preparatory move.** Sets up a later move (e.g. clearing a square).
- **Retreat / regroup.** Moving a piece back to improve placement (not a loss of tempo if it serves a purpose).
- **Counter-threat / counter-attack.** Answering a threat with a larger threat.
- **Intermediate move (zwischenzug, in-between move).** A forcing move inserted before the expected recapture or reply, changing the evaluation. **Intermezzo** is the same idea.
- **Desperado.** A piece that is doomed captures or checks as much as it can before being lost (or a piece traded at a gain after being attacked).
- **Only move / only good move.** The single move that holds the position or keeps a result.
- **Best defence / critical line / refutation.** The strongest defensive reply; the line on which a plan or sacrifice is tested; the sequence showing an idea fails.
- **Novelty (TN, theoretical novelty).** First time a move is played in a known line.
- **Candidate moves.** Moves worth analysing seriously (Kotov: list candidates, analyse each line once).
- **Variation / main line / sideline / transposition.** Branches in analysis.

---

## 4. Tactical motifs

Each motif is a pattern by which material or mate is won.

- **Fork.** One piece attacks two or more targets at once. (Knight fork, pawn fork, queen fork, king fork.)
- **Pin.** A piece cannot move without exposing a more valuable piece behind it.
  - *Absolute pin:* the piece behind is the king; the pinned piece may not move off the line.
  - *Relative pin:* the piece behind is something valuable but not the king; moving it is legal but costly.
  - *Cross-pin:* a piece is pinned on two lines.
- **Skewer.** The reverse of a pin: the more valuable piece is attacked first and must move, exposing the lesser piece behind.
- **Discovered attack.** Moving one piece unblocks a line of another, which then attacks. **Discovered check** is the same with check. **Double check** is a discovery in which the moved piece also checks.
- **Battery.** Two or more pieces aligned on one line so the rear one adds force behind the front (e.g. queen behind bishop on a diagonal; two rooks on a file; queen + rook). Types:
  - *Direct battery:* both pieces attack along the line.
  - *Discovery battery:* the front piece moves away with tempo, revealing the rear.
  - *Alekhine's gun:* two rooks doubled on a file with the queen behind them.
  - *Pawn battery* (informal): pawns on a chain supporting an advance.
- **X-ray.** A piece "sees through" another piece along a line: attacks or defends a square or piece behind the intervening one.
- **Deflection.** Forcing a defending piece away from its post or duty.
- **Decoy / attraction.** Forcing a piece (often the king) onto a bad square.
- **Interference.** Placing a piece on a line between an enemy piece and what it protects, cutting the connection.
- **Overloading.** A defender has too many duties; satisfying one leaves another unprotected.
- **Removal of the guard (undermining).** Capturing or driving off the piece that defends a target.
- **Clearance (line-clearing, square-clearing).** Moving a piece away to free a square or line for another.
- **Zwischenzug.** See §3.
- **Trapped piece.** A piece with no safe squares that can be won.
- **Hanging piece / en prise.** A piece that is attacked and insufficiently defended / that can be captured. *Loose piece* = undefended piece.
- **Windmill.** A repeated discovered-check/capture mechanism that wins material every move.
- **Perpetual check.** Unending checks leading to a draw.
- **Sacrifice** (see §6). **Quiet sacrifice**: a sacrifice without immediate forcing follow-up.
- **Greek gift.** Bishop sacrifice on h7 (h2) against a castled king.
- **Back-rank weakness.** King on its first rank trapped by its own pawns without luft.
- **Luft.** A flight square made by a pawn move so the king escapes back-rank mate.
- **Smothering.** King immobilised by its own pieces (see mating patterns).
- **Promotion tactics.** Pawn races, promotion with check, underpromotion to avoid stalemate or deliver mate.
- **Stalemate trick (swindle).** Sacrificing material to force stalemate when losing.
- **Tactical shot / combination.** A sequence of forcing moves with a concrete gain. A *combination* typically involves a sacrifice.
- **Counting tactics.** Counting attackers and defenders on a square to decide a capture sequence (exchange sequence).
- **Capture sequence, exchange sequence.** The chain of captures and recaptures on one square.
- **Overprotection.** More defenders than needed on a key point (Nimzowitsch).
- **Counter-sacrifice.** Answering a sacrifice with one's own.
- **Rook lift.** Raising a rook to the third/fourth rank to swing along it for attack.

---

## 5. Mating patterns

Named final positions. Descriptions are minimal.

- **Fool's mate.** Fastest mate: 1.f3 e5 2.g4 Qh4# (and variants). Black can do it after 1.f4/f3 plus g4.
- **Scholar's mate.** Queen and bishop target f7/f2.
- **Back-rank mate.** Rook or queen mates a king hemmed in by its own pawns.
- **Smothered mate.** Knight mates a king surrounded by its own pieces.
- **Anastasia's mate.** Knight and rook; king on the edge, blocked by its own pawn.
- **Arabian mate.** Rook and knight; king in the corner.
- **Boden's mate.** Two bishops on crossing diagonals mate a king boxed in by own pieces.
- **Damiano's mate.** Queen supported by a pawn mates on the h-file/h7 area.
- **Epaulette mate.** King flanked by own rooks; queen mates on the file in front.
- **Greco's mate.** Bishop and rook/queen on the h-file against a castled king.
- **Hook mate.** Rook, knight and pawn.
- **Lolli's mate.** Queen supported by a pawn on g7 against a fianchetto.
- **Légal's mate.** After a queen sacrifice, two knights and a bishop mate.
- **Morphy's mate.** Bishop and rook.
- **Opera mate.** Bishop and rook (Morphy v. Duke of Brunswick and Count Isouard).
- **Pillsbury's mate.** Bishop and rook, bishop covers g7.
- **Swallowtail (Guéridon) mate.** Queen in front of the king, own pieces blocking the others.
- **Blackburne's mate.** Bishops and knight.
- **Dovetail (Cozio's) mate.** Queen mates diagonally adjacent, own pieces blocking the flight squares.
- **Kill Box mate.** Rook adjacent to king, supported by queen.
- **Vukovic mate.** Rook and knight.
- **Reti's mate.** Bishop and rook, king blocked by own pieces.
- **Suffocation mate.** Knight and bishop.
- **Triangle mate.** Queen with the king confined to a triangle of squares.
- **Corridor mate.** King confined to a corridor; rook mates.
- **Ladder / lawnmower mate.** Two rooks (or queen and rook) check along adjacent ranks/files, pushing the king to the edge.
- **Mayet's, Balestra, Max Lange, David and Goliath, Anderssen's** and others: named patterns of piece pairs; each is a specific arrangement of mating piece + guards.
- **Basic endgame mates.** K+Q v K, K+R v K, K+2B v K, K+B+N v K (W-manoeuvre, Delétang). K+2N v K does not force mate.

---

## 6. Material, exchange and sacrifice

- **Material.** Pieces on the board. Conventional values (pawns): P=1, N=3, B=3, R=5, Q=9; the king has no exchange value.
- **Relative value.** Value depends on position (open/closed, endgame/middlegame). Knights prefer closed positions; bishops prefer open.
- **Material imbalance.** Unequal piece types, not necessarily unequal totals (e.g. rook vs two minors).
- **Bishop pair.** Having both bishops, valuable in open positions.
- **Exchange.** Trading pieces. **Exchange up / down** = winning/losing the exchange, i.e. rook for minor piece. **Exchange sacrifice** = giving a rook for a minor piece for compensation.
- **Trade-down.** Simplifying when ahead.
- **Principle of two weaknesses.** One weakness can be defended; two stretch the defender.
- **Sacrifice types.** Positional (for structure/squares), tactical (concrete gain), speculative (unsound but practical), sham sacrifice (regains material by force), *decoy*, *deflection*, *clearance*, *desperado*, *Greek gift*.
- **Sound vs unsound.** A sound sacrifice is justified by best play; unsound is refutable.
- **Material count vs activity.** Activity can outweigh a small material deficit.
- **Qualitative advantage.** Better piece quality or placement with equal material.
- **Passed-pawn compensation.** A pawn majority/passed pawn can compensate for a piece in the endgame.

---

## 7. Pawn structure

- **Pawn chain.** Diagonally linked pawns; *base* = rearmost, *head* = front.
- **Pawn island.** A group of pawns separated from others by empty files.
- **Isolated pawn.** No friendly pawn on adjacent files. (**IQP** = isolated queen's pawn.)
- **Doubled pawns.** Two on the same file. **Tripled** = three.
- **Backward pawn.** Behind adjacent pawns, cannot safely advance, and its stop square is controlled by the enemy.
- **Passed pawn.** No enemy pawn ahead on its file or adjacent files. **Protected passed pawn** = defended by another pawn. **Candidate pawn** = could become passed through exchange. **Outside passed pawn** = on a wing far from the main action, a decoy for the enemy king.
- **Connected pawns.** On adjacent files, side by side or on adjacent ranks. **Connected passed pawns** are strong.
- **Hanging pawns.** Two adjacent pawns (typically c+d) with no neighbour pawns, controlling the centre but vulnerable.
- **Pawn majority / minority.** More pawns on one wing. **Minority attack** = advancing fewer pawns against more to create a weakness (e.g. b4-b5 vs c6).
- **Pawn storm.** Advancing pawns against the enemy king.
- **Pawn break / lever.** A pawn advance to challenge an enemy pawn and open lines. **Tension** = two pawns that can capture each other.
- **Fixed pawns / locked pawns.** Pawns blocked head-to-head. **Blocked / closed structure.**
- **Pawn wedge.** A protected advanced pawn splitting the opponent's position.
- **Space advantage.** More territory behind advanced pawns.
- **Holes / weak squares.** Squares that can no longer be controlled by pawns. **Outpost** = a piece placed on a hole that cannot be driven off by pawns.
- **Colour complex.** Weak squares of one colour, often when bishops are traded.
- **Named structures** (used as shorthand): Carlsbad, Maroczy Bind, Hedgehog, Stonewall, Sicilian Dragon / Scheveningen / Najdorf setups, Isolani, Hanging pawns, Caro-Kann / French chains, Benoni structure.
- **Pawn-structure imbalance.** Asymmetry that gives each side targets.
- **Rook pawn.** A pawn on the a- or h-file. Special endgame properties (§13).
- **Zugzwang in pawn endings.** See §14.

---

## 8. Piece placement and piece relationships

- **Attack / defence / protection.** A piece *attacks* a square it can capture on; *defends* a friendly piece standing there.
- **Control of a square.** Having more or better attackers than the opponent on it.
- **Mobility.** Number of safe squares available to a piece.
- **Activity.** How effectively a piece is placed.
- **Centralisation.** Placing pieces in or near the centre.
- **Coordination / harmony.** Pieces supporting each other's tasks.
- **Outpost.** See §7. Classic: a knight on a protected advanced square.
- **Open file / half-open (semi-open) file.** No pawns / only enemy pawns on a file. **Rook on the open file.**
- **Seventh rank (second rank).** A rook on the 7th attacks pawns and confines the king. **Pig** = rook on the 7th; two rooks there = *pigs*.
- **Doubling rooks.** Two rooks on one file.
- **Good bishop / bad bishop.** A bad bishop is hemmed by its own pawns on its colour. **Opposite-coloured bishops.** Each side has one bishop and they travel on different colours.
- **Fianchetto.** Developing a bishop on b2/g2/b7/g7.
- **Knight vs bishop.** Contest of piece types depending on structure.
- **Knight on the rim.** Poorly placed.
- **Rook behind a passed pawn.** Tarrasch rule: rooks belong behind passed pawns.
- **Active vs passive.** Piece placement in terms of options.
- **Piece quality.** Combination of mobility, targets and safety.
- **Domination.** Restricting a piece so every move loses (e.g. trapping a knight with a bishop or king; *Dvoretsky*). Variant of trapped piece.
- **Restraint / blockade.** Preventing an enemy pawn or piece from advancing by occupying/controlling its stop square (Nimzowitsch).
- **Overprotection.** See §4.
- **Mined / vulnerable squares.** Squares where a piece placed would be captured at a loss.
- **Interposition / block.** Placing a piece between an attacker and its target.
- **Line, ray, diagonal, file, rank.** Lines of movement for sliding pieces.
- **Weak piece.** Loose, tied to defence, or low mobility.
- **Piece tied to defence.** Immobilised by a duty.

---

## 9. Squares, lines, and space

- **Centre.** d4, d5, e4, e5. *Extended centre* includes c3–f6.
- **Control vs occupation of the centre** (hypermodern vs classical).
- **Weak squares / strong squares / key squares.**
- **Critical square / key square** (endgame): a square which, if reached by the king, guarantees a result (§13).
- **Entry squares.** Squares through which pieces penetrate (e.g. on the 7th).
- **Open lines.** Files, ranks, diagonals clear of pawns.
- **Long diagonal.** a1–h8 / h1–a8.
- **Colour complex.** All squares of one colour in a region.
- **Space.** Territory behind one's pawn front.
- **Flank vs centre.** Where play occurs.
- **Kingside / queenside; opposite-side castling.**
- **Back rank / first rank / luft.** See §4.
- **Corner / edge.** King and knight limitations.
- **Diagonal / file opening** by pawn breaks.

---

## 10. King safety and king play

- **King safety.** Degree of exposure to attack.
- **Castling** (short / long) as a safety device; **castling rights** as a resource.
- **Pawn shelter / pawn shield.** Pawns in front of the castled king; *weakened shelter* after moves such as g-pawn advance.
- **King in the centre.** Uncastled, often vulnerable.
- **King hunt.** Chasing an exposed king across the board.
- **King walk.** King travelling voluntarily, often in endgames.
- **Active king.** In endgames the king is a fighting piece.
- **Opposition.** See §13.
- **Escape squares / flight squares.** Squares the king can move to.
- **Luft.** See §4.
- **Opposite-side castling.** Mutual pawn storms; races.
- **Fortress around the king.** Defensive structures.
- **Zugzwang for the king.** See §14.
- **King as a blocker.** Blockading a pawn.
- **Shouldering (body check).** King uses the board geometry to keep the opposing king away from a pawn.
- **Cutting off the king.** Restricting the enemy king along a file or rank.

---

## 11. Strategic concepts and principles

- **Plan.** A goal-directed sequence of moves chosen from the position's features.
- **Evaluation.** Judgement of who is better, by how much, and why.
- **Steinitz's principles.** Accumulate small advantages; do not attack without a prior advantage; the defender must neutralise; the centre.
- **Capablanca's style.** Simplify, endgame technique, liquidation.
- **Nimzowitsch's ideas.** Prophylaxis, blockade, overprotection, restraint.
- **Hypermodern vs classical.** Control the centre from afar vs occupy it with pawns.
- **Kotov's analysis method.** List candidate moves; examine each tree branch once.
- **Dvoretsky's ideas.** Prophylactic thinking, critical positions, trading rooks, domination.
- **Principle of the two weaknesses.** See §6.
- **Principle of activity.** Increase the activity of your pieces and reduce your opponent's.
- **Weakness.** Something that can be attacked: a pawn, square, king, piece.
- **Targets.** Weaknesses selected for attack.
- **Initiative vs attack.** See §2.
- **Counterplay.** Active chances for the defender.
- **Attack on the king.** Aggression directed at the king.
- **Consolidation.** Securing an advantage before pressing on.
- **Conversion.** Turning a won position into a win.
- **Technique.** Reliable method of converting an advantage.
- **Simplification.** Trade pieces to reduce counterplay.
- **Liquidation.** Trading into an endgame.
- **Overextension.** Advancing too far without support.
- **Dynamic equality.** Balance maintained through active play.
- **Static weakness vs dynamic compensation.**
- **Strategic sacrifice, positional sacrifice.**
- **Prophylaxis.** See §3.
- **Tension.** See §7; also maintaining unresolved contact between pieces.
- **Closed / open / semi-open position.** Character determined by pawn structure.
- **Phases.** Opening, middlegame, endgame.
- **Transition.** Moving between phases, e.g. exchanging into an endgame.
- **Imbalance.** Differences in material, structure, king safety, activity.
- **Harmony.** Pieces cooperate toward one plan.
- **Zwischenzug in strategy.** Intermediate manoeuvres.

---

## 12. Opening concepts

- **Principles.** Control the centre, develop pieces, castle, avoid moving a piece twice or the queen early, connect rooks.
- **Open game (1.e4 e5), semi-open (1.e4 other), closed (1.d4 d5), semi-closed (1.d4 other), flank openings (1.c4, 1.Nf3, 1.f4 ...).**
- **Opening names.** Named by ECO code (A00–E99). Bird's Opening 1.f4 is A02–A03. 1.f3 is the Barnes / Gedult Opening (A00).
- **Gambit; gambit accepted/declined; counter-gambit.** E.g. From's Gambit: 1.f4 e5 2.fxe5 d6.
- **Theory.** Established best-known lines. **Tabiya** = a standard position reached in a line.
- **Novelty.** See §3.
- **Repertoire.** A player's chosen set of openings.
- **Transposition.** Entering a known line by another move order.
- **Opening trap.** A plausible move that loses to a tactic (e.g. fool's mate patterns, Légal's trap, Fishing pole, Elephant trap).
- **Symmetrical positions.** Both sides mirror; the side to move has a tempo, which can be either an asset or a liability.
- **Copycat / mirror strategy.** Replying symmetrically; fails when a first-mover threat appears.
- **Development race.** Who completes development first.
- **Tempo in openings.** Counting moves to reach a setup.
- **Pawn structure from the opening.** Determines middlegame plans.
- **Opening book** (software) / **opening explorer** (database statistics).
- **Mainline / sideline / deviation.**
- **Hypermodern setups:** fianchetto systems, King's Indian, Nimzo-Indian.
- **Colour reversal.** Same structure with colours exchanged (reversed openings, e.g. Bird's as a reversed Dutch).

---

## 13. Endgame concepts

- **Opposition.** Kings face each other with one square between; the side not to move has the opposition. *Direct*, *distant* (more squares, same file/rank), *diagonal*.
- **Key (critical) squares.** If the attacking king reaches them, the pawn promotes regardless of who moves.
- **Rule of the square.** A king can catch a pawn if it can step into the square drawn from the pawn to its promotion rank.
- **Triangulation.** The king loses a tempo by travelling a triangle of squares, changing who is to move while the position stays.
- **Corresponding squares.** Pairs of squares for the two kings such that if one king stands on one, the other must stand on its counterpart to hold the position.
- **Zugzwang / mutual zugzwang.** See §14.
- **Outside passed pawn.** Draws the defending king away; then the attacking king takes the other wing.
- **Pawn majority; creating a passed pawn.**
- **Breakthrough.** A pawn sacrifice that creates a passed pawn.
- **Lucena position.** Winning technique with rook + pawn on the 7th: *building a bridge*.
- **Philidor position.** Drawing technique: rook on the 6th rank to restrict, then check from behind.
- **Vancura position.** Drawing defence against a rook pawn.
- **Tarrasch rule.** Rooks behind passed pawns.
- **Saavedra position.** Underpromotion to rook to win (avoiding stalemate).
- **Wrong-coloured rook pawn / wrong bishop.** Bishop and rook pawn draw if the bishop does not control the promotion square colour.
- **Fortress.** A setup the superior side cannot break (§15).
- **Réti manoeuvre.** King moves along a diagonal threatening two goals at once.
- **Troitsky line.** Two knights v pawn: winning zone limits.
- **Kling and Horwitz** etc. Named studies.
- **Opposite-coloured bishop endings.** Strong drawing tendencies.
- **Same-coloured bishop endings; good/bad bishop.**
- **Knight endings;** akin to pawn endings.
- **Rook endings.** Most common; many standard positions.
- **Queen vs pawn on the 7th.** Draws with rook or bishop pawns.
- **Queen vs rook.** Theoretically often a win, practically hard.
- **Minor piece vs pawns.**
- **Mating technique.** K+Q, K+R, K+2B, K+B+N.
- **Tablebases.** Exact solved endings (§16).
- **Active king, passive king.**
- **Shouldering; cutting off; outflanking.** King manoeuvres.
- **Pawn race.** Both sides promote; counting tempi.
- **Underpromotion** to dodge stalemate or give mate.
- **Zugzwang in endings.** The dominant theme.
- **Simplification into a won ending.**

---

## 14. Zugzwang and related

- **Zugzwang.** Any legal move worsens the position of the side to move (the obligation to move is a disadvantage). *Whoever must move loses something.*
  - *Mutual zugzwang (reciprocal zugzwang).* Whoever is to move loses.
  - *Trébuchet.* Mutual zugzwang in a pawn ending where each king attacks the other's pawn.
  - *Zugzwang vs forced move.* Zugzwang needs legal moves that all harm; stalemate is the case where there are no legal moves (draw).
  - *Partial zugzwang / positional zugzwang* (some moves restricted).
  - *Tempo move.* See §2.
  - *Reserve tempi.* The count of spare moves each side has before being forced to commit.
- **Domination.** See §8.
- **Opposition.** See §13.
- **Triangulation.** Method to transfer the move.
- **Waiting strategy.** Keeping a harmless move available.
- **Zugzwang as a goal.** Reducing the opponent's safe moves.
- **Mutual zugzwang in fortresses.**

---

## 15. Draw mechanisms and defensive concepts

- **Fortress.** A defensive setup which cannot be broken though the opponent is materially ahead.
- **Perpetual check.** Repeated checks.
- **Stalemate (resource).** Forcing stalemate.
- **Repetition.** See §1.
- **Blockade.** A position where nothing can move.
- **Dead position.** See §1.
- **Insufficient material.** See §1.
- **Fifty-move draw.** See §1. Not followed by the programme.
- **Opposite-coloured bishops.** See §13.
- **Wrong-rook-pawn draws.**
- **Perpetual attack on a piece.** Repeating threats.
- **Draw by agreement; "drawish".**
- **Counter-attack as defence.**
- **Prophylactic defence; active defence vs passive.**
- **Hold with the only move.**
- **Fortress-building;** pawn, piece, or king fortresses.
- **Swindle / trick.** A hope-based resource, e.g. a stalemate trap.
- **Simplification into a draw.**
- **Equilibrium.** A position that remains balanced under best play.

---

## 16. Game-theoretic and solving concepts

- **Chess as a game.** Two-player, zero-sum, finite, perfect information, no chance.
- **Game-theoretic value.** Win / draw / loss under optimal play from a position (WDL).
- **Minimax.** Maximise own worst-case outcome.
- **Game tree.** Nodes = positions, edges = legal moves. **Branching factor** ~ 30–35 on average. **Depth** in plies.
- **Perfect play, optimal play.**
- **Solving levels (Allis).**
  - *Ultra-weakly solved:* the value of the start position is known.
  - *Weakly solved:* value known and a strategy achieving it from the start is known.
  - *Strongly solved:* value known for every legal position.
  - Chess is unsolved at the start; endings of up to 7 pieces are strongly solved by tablebases.
- **Tablebase.** Database of exact values for all positions of a material set. **Syzygy** (WDL and DTZ), **Nalimov** (DTM), **Lomonosov** (7-piece). **DTM** = depth to mate; **DTZ** = distance to zeroing (pawn move or capture); **WDL** = win/draw/loss.
- **Retrograde analysis.** Building a tablebase backwards from terminal positions.
- **Forced win / forced line / forced sequence.**
- **Opening tree / theory tree.** The branching of known opening variations.
- **Transposition table.** Cache of positions already evaluated.
- **Zobrist hashing.** Position hash used for transposition tables.
- **Position repeats and path dependence.** Repetition draws depend on history, not just the board.
- **History heuristic / graph-history interaction (GHI).** The issue that a position's value depends on its path because of repetition rules.
- **Strategy-stealing, pairing strategies** (as arguments in other games; chess has zugzwang so no general strategy-stealing).
- **Symmetry.** Colour symmetry of the start; side to move breaks it.
- **Draw death.** At top levels, draws common.
- **Complexity.** State-space estimates (about 10^44 legal positions; Shannon number ~10^120 for game-tree size).
- **Proof-number search, retrograde, endgame databases.**

---

## 17. Engine and analysis concepts

- **Evaluation function.** Static score of a position. Units: **centipawns** (1/100 of a pawn). Neural evaluations (NNUE) in modern engines.
- **Search.** Minimax with **alpha-beta pruning**, **iterative deepening**, **quiescence search** (extend through captures/checks), **null-move pruning**, **late move reductions**, **aspiration windows**, **transposition table**, **killer moves**, **history heuristic**.
- **Horizon effect.** Search cut-off hides an unavoidable loss beyond the depth.
- **Principal variation (PV).** Best line found. **MultiPV** = several best lines at once.
- **Depth, seldepth, nodes, nps.** Search measures.
- **WDL model.** Evaluation expressed as win/draw/loss probabilities.
- **Mate score vs eval score.** Mate in N distinct from centipawn values.
- **Engine assertion vs proof.** An engine evaluation is a heuristic claim, not a proof, except in tablebase and forced-mate cases. (The programme classifies certificates as MATE-VERIFIED, TB-EXACT, ENGINE-ASSERTION.)
- **Forced mate.** Mate that the opponent cannot avoid.
- **Blunder check.** Detecting a tactical oversight.
- **Move quality classes:** best, excellent, good, inaccuracy, mistake, blunder (by centipawn loss thresholds).
- **Accuracy score / average centipawn loss (ACPL).**
- **Opening database / explorer / book.** Statistics from games.
- **Cloud analysis.** Precomputed evaluations.
- **Fixed-node / fixed-depth analysis.** For reproducibility.

---

## 18. Error and judgement vocabulary

- **Inaccuracy, mistake, blunder.** Increasing severity of error.
- **Oversight / miscalculation.**
- **Missed win / missed tactic.**
- **Positional error.** Strategic misjudgement.
- **Time trouble blunder.**
- **Good / bad / dubious / interesting move** (annotation symbols: !, ?, !!, ??, !?, ?!, =, ±, ∓, +-, -+, ∞ unclear).
- **Equal, slightly better, clearly better, winning.** Verbal evaluation scale.
- **Unclear.** Dynamic imbalance, evaluation uncertain.
- **Compensation.** See §2.
- **Practical chances.** Difficulty for the opponent.
- **Sound / unsound.** See §6.

---

## 19. Terms with multiple meanings (disambiguation)

- **Tempo.** One move; a gained/lost move; a move that passes the obligation; a spare move; a pawn-race count; an engine bonus (§2).
- **Initiative.** Dictating play with threats; not synonymous with attack or advantage.
- **Attack.** (1) Piece attacks a square/piece. (2) A campaign against the king. (3) A threat to win material.
- **Threat.** Immediate (next move) vs strategic (longer-term).
- **Control.** Attack on a square; domination of an area; occupation.
- **Weakness.** Pawn weakness, square weakness, king weakness, colour-complex weakness.
- **Outpost.** A square not attackable by enemy pawns; the piece on it.
- **Exchange.** A trade; also specifically rook vs minor piece.
- **Open.** An open file, an open position, an open game, an open king.
- **Closed.** Pawn structure blocked; closed opening.
- **Pin, skewer, x-ray.** Differ by which piece is front/back and whether the front is more valuable.
- **Battery.** Pieces aligned on a line; discovered attack. Different from the "Alekhine's gun" specifically.
- **Fortress.** A defence structure; sometimes a drawn structure in a materially unbalanced position.
- **Zugzwang.** See §14.
- **Prophylaxis.** Prevention of an opponent's idea; not passive.
- **Sacrifice.** Voluntary material loss; includes gambit, exchange sacrifice.
- **Compensation.** Non-material offsets.
- **Dynamic vs static.** Time-dependent vs structural factors.
- **Forced.** Only legal move; only good move; forced by threats; forced in analysis (opponent has no choice without losing).
- **Position / structure / formation.** Position = full state; structure = pawn skeleton; formation = piece setup.
- **Line.** A board line (file/rank/diagonal) or a sequence of moves (variation).
- **Theory.** Known analysis of an opening or ending.
- **Equal / equilibrium / balance.** Evaluation of 0; stable balance maintained by best play.

---

## 20. Touchpoints with the programme (assistant's reading, not admitted)

These are *names that sit near* the existing catalog terms, listed so the user can decide whether any is useful for articulation. **None is proposed as a primitive**, and nothing here is a mapping the user has stated. Under the admission rule a candidate must be a race condition with terms, irreducible, not a description of a result, and piece/board-agnostic; the survey above mostly lists results and descriptions, which the rule would classify as terms, inputs, definitions or results.

- Near **P-001 (mating/capture net):** fork, pin, skewer, trapped piece, hanging piece / en prise, loose piece, double check, smothering, mating patterns (§5), back-rank weakness, flight squares / luft.
- Near **P-003 (chained dependencies):** overloading, removal of the guard, deflection, decoy, interference, relative pin, discovered attack, clearance, overprotection, pawn chain base/head, capture sequence counting.
- Near **P-008 (coverage race under disturbance):** control of a square, key squares, outpost/hole, domination, mined squares, flight squares, colour complex, fortress, opposition and corresponding squares, restraint/blockade.
- Time-race vocabulary that names a *race*: tempo, tempo count, pawn race, rule of the square, lead in development, initiative, triangulation, reserve tempi.
- Recorded residual in the catalog: **zugzwang without capture** (zugzwang, mutual zugzwang, trébuchet, corresponding squares, §14). It remains a residual, not a primitive.
- Equilibrium vocabulary: dynamic equality, fortress, perpetual check, repetition, stalemate, dead position, insufficient material, blockade, opposite-coloured bishops (§15). These name draw *results*; the programme's boundary-state definition requires a cycle realization or a forced terminal draw.
- Solving vocabulary: WDL, DTZ/DTM, tablebase, weak solution, transposition table, GHI (§16–17). Path-dependence through repetition (GHI) is relevant to how a proof search caches a verdict that rests on a cycle. For a closed region the draw value does not depend on it, because the repetition rule can only end a game as a draw (the programme does not follow the 50-move rule). It bears on the certification of exits, which must not rest on a cycle argument, and on the cache key for a cycle, which the designer expects to discover (survey §10).

**Caution on use.** If a term here is felt to be a candidate primitive, the check is the admission rule in the catalog §1, not its presence in this list.
