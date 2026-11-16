# KRvKR seed batches: column guide

Two files, one per region. Together they partition every legal KRvKR position (21,561,456 raw positions, 2,695,812 classes).

- `krvkr_choice_region.csv.gz`: drawn positions (1,892,778 classes, 15,138,064 raw positions).
- `krvkr_forced_region.csv.gz`: forced positions with d and favoured side (803,034 classes, 6,423,392 raw positions).

## Shared columns

| column | meaning |
|---|---|
| `id` | index of the class representative: `(((wK*64+wR)*64+bK)*64+bR)*2+stm`, square = 8*rank+file with a1 = 0, stm 0 = white to move |
| `class16` | id of the class under the 16-image fold (8 board symmetries x colour swap). Rows with equal `class16` are colour mirrors of each other |
| `wK wR bK bR` | squares of the representative |
| `stm` | side to move, `w` or `b` |
| `orbit` | number of raw positions the class stands for (distinct images under the 8 board symmetries; colour is not folded) |

## Forced region

`d` is the number of plies from the position to the first irreversible event of the segment, which is the capture of a rook or mate (mate gives d = 0). Certification is scoped to the KRvKR segment and ends at that terminal, whose reason is the symbol in `tag`. It is not DTM.

| column | meaning |
|---|---|
| `d` | plies to the terminal reason named by the symbol in `tag`, optimal play from both sides |
| `side` | favoured colour, `w` or `b`. Odd d means the side to move is favoured, even d means the other side is |
| `tag` | `d`, then `side`, then one terminal symbol, for example `5wC` (legend below) |
| `mover_favoured` | 1 if the side to move is the favoured side |
| `winner_moves` | moves the favoured side makes within d: (d+1)/2 if it is to move, d/2 otherwise |
| `moves` | legal moves for the side to move |
| `to_mover_win`, `to_draw`, `to_opp_win` | those moves split by outcome: still winning for the mover, into the drawn region (a capture whose lower-segment result is a draw counts here), into a position the opponent wins |
| `opt` | number of d-optimal moves (win in d-1 for the favoured mover, maximum resistance for the other) |

### Terminal symbol (legend)

The symbol states the terminal reason at which d ends, over the d-optimal lines (the favoured side plays only d-optimal moves, the other side only maximum-resistance moves). Each entry gives the terminal kind, the segment entered, the status of the entry, and what d counts to.

| symbol | terminal kind | segment entered | status of the entry | what d counts to |
|---|---|---|---|---|
| `C` | capture of a rook: an irreversible entry into a lower segment, as a win for the capturer | KRvK | certified, valued from `KRvK.rtbw` (the capture rule matched that table on 223,944 positions). These results depend on the KRvK segment | plies to the capture |
| `M` | mate inside this segment, which ends the game | none | certified within this segment. Depends on nothing below | plies to the mate. A mated position has d = 0 |
| `B` | both terminals are reachable on d-optimal lines | KRvK, for the capture lines | as for `C` and `M`, line by line | plies to the first terminal on a d-optimal line, which is the same d for both |

The symbols are labels chosen for this file only. No other symbol set occurs, so no sets other than `B` are listed.

A capture whose KRvK result is a draw (recapture or stalemate) is not a terminal of this region. It is an edge back to the choice region (`to_draw` here, `stay_leaf` in the other file).

## Choice region

| column | meaning |
|---|---|
| `moves` | legal moves for the side to move |
| `stay` | moves that keep the position drawn, including captures whose lower-segment result is a draw |
| `stay_leaf` | how many of the `stay` moves are such captures |
| `dev` | losing deviations of the mover: moves into a position the opponent wins |
| `dev_d_min`, `dev_d_max` | after the cheapest and the dearest losing deviation, the plies in which the opponent wins (blank when `dev` = 0) |
| `stay_with_dev` | drawing moves into a (non-capture) drawn position where the other side has a losing deviation |
| `rests_on_leaf` | 1 if every `stay` move is a capture whose draw is decided by the lower segment |

## Provenance and checks

- d, favoured side and all move counts come from a table-free retrograde search (C++), with captures as leaves. Syzygy tables were used only as a check.
- The capture-leaf rule agrees with `KRvK.rtbw` on all 223,944 KRvK positions it was tested on.
- WDL agrees with `KRvKR.rtbw` on all 2,695,812 classes. There are no cursed wins or blessed losses.
- For positions where the mover wins, |DTZ| equals d in all 784,918 classes. For positions where the mover loses, |DTZ| is d-1 in all 17,116 classes, and 1 for mated positions. This is a counting convention of the stored DTZ, not a disagreement about the result.
- A separate python-chess minimax confirmed d on 2,086 stratified samples with 0 failures.
- The result and d are identical across all 16 images of every class. Move counts match python-chess on 30,000 sampled classes.
- Draw closure holds for every drawn class: no winning move, and at least one move that stays drawn.
- Terminal types were measured with a second pass of the same search. It reproduced d and the results exactly.
- The 50-move rule is not applied.
