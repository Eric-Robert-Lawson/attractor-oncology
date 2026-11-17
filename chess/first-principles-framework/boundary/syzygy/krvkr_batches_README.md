# KRvKR seed batches: README of the irreversible segmented boundary

Three files. The two regions together partition every legal KRvKR position (21,561,456 raw positions, 2,695,812 classes under the 8 board symmetries, colour kept).

- `krvkr_forced_region.csv.gz`: forced positions (803,034 classes, 6,423,392 raw positions), each with its d-depth header.
- `krvkr_choice_region.csv.gz`: drawn positions (1,892,778 classes, 15,138,064 raw positions), the position only. Either region determines the other, given the enumeration of legal positions.
- `krvkr_terminal_states_legend.csv.gz`: the legend file, the table of terminal board states (25,862 entries).

## Position columns (both region files)

| column | meaning |
|---|---|
| `id` | index of the class representative: `(((wK*64+wR)*64+bK)*64+bR)*2+stm`, square = 8*rank+file with a1 = 0, stm 0 = white to move |
| `wK wR bK bR` | squares of the representative |
| `stm` | side to move, `w` or `b` |

## The d-depth header (forced region file)

The header is `[d][index][symmetry condition]`, in the columns `d`, `index`, `sym`.

| column | meaning |
|---|---|
| `d` | plies from the position to the terminal condition, to the first irreversible event of the segment (the capture of a rook, or mate; a mated position has d = 0), under optimal play from both sides. Certification is scoped to this segment and ends at that terminal. It is not DTM |
| `index` | row of the legend file: the terminal board state at the end of the recorded move chain |
| `sym` | the transformation applied to the stored state to give the terminal state reached (symbols below) |

The favoured side is not stored. Odd d means the side to move is favoured, and even d means the other side is. The condition type is given by the entry that `index` points to.

When several equally viable terminals can be reached on d-optimal lines, the recorded one is the lowest `(index, sym)` among all d-optimal lines from the position. Mate entries come before capture entries in the legend file, so a position that can end by mate or by capture is recorded as mate. The other routes are reachable as anchor points and connect by the segmented dependency move chain.

## The legend file

One row per terminal state, stored once per symmetry class. The 8 board symmetries and the white/black exchange are folded, so colour is not stored.

| rows | type | columns `s1 s2 s3 s4` |
|---|---|---|
| 0 to 602 | `M` | the mated position: king and rook of the mated side, then king and rook of the mating side |
| 603 to 25,861 | `C` | the position after the capture, in the lower segment: king and rook of the capturing side, then the king of the captured side. `s4` is empty |

What each type declares:

| type | terminal kind | segment entered | status of the entry | what d counts to |
|---|---|---|---|---|
| `M` | mate inside this segment, which ends the game | none | certified within this segment. Depends on nothing below | plies to the mate. A mated position has d = 0 |
| `C` | capture of a rook: an irreversible entry into a lower segment, as a win for the capturer | KRvK, with the weak side to move | certified, valued from `KRvK.rtbw` (the capture rule matched that table on 223,944 positions). Results of this type depend on the KRvK segment | plies to the capture |

A capture whose KRvK result is a draw (recapture or stalemate) is not a terminal. It is an edge back to the choice region.

## Symmetry symbols

The group in force for this segment is the 8 board symmetries (no pawn on the board). A transformation is built from a file flip (`a` to `h`), a rank flip (`1` to `8`) and a transposition, in that order. A fixed point of several transformations takes the lowest symbol.

| symbol | transformation |
|---|---|
| `a` | identity |
| `b` | file flip |
| `c` | rank flip |
| `d` | file flip and rank flip (half turn) |
| `e` | transposition over the a1-h8 diagonal |
| `f` | file flip, then transposition (quarter turn) |
| `g` | rank flip, then transposition (quarter turn the other way) |
| `h` | file flip and rank flip, then transposition (reflection over the other diagonal) |

## Provenance and checks

- d, the results and the terminals come from a table-free retrograde search (C++), with captures as leaves. Syzygy tables were used only as a check.
- The capture-leaf rule agrees with `KRvK.rtbw` on all 223,944 KRvK positions it was tested on.
- WDL agrees with `KRvKR.rtbw` on all 2,695,812 classes. There are no cursed wins or blessed losses.
- For positions where the mover wins, |DTZ| equals d in all 784,918 classes. Where the mover loses, |DTZ| is d-1 in all 17,116 classes, and 1 for mated positions. This is a counting convention of the stored DTZ.
- A separate python-chess minimax confirmed d on 2,086 stratified samples with 0 failures. The result and d are identical across the images of every class.
- Draw closure holds for every drawn class: no winning move, and at least one move that stays drawn.
- The forced file was rebuilt with the header and then read back as a user would. For all 803,034 rows: d equals the retrograde label; the terminal rebuilt from `index` and `sym` is reached by a d-optimal line of d plies; and the lowest `(index, sym)` over all d-optimal lines, found by a direct walk, equals the recorded one (0 differences). The class set, d and favoured side are unchanged from the earlier forced file. Terminal types: 793,136 classes record a capture, 8,292 a mate where only mate is reachable, and 1,606 a mate where mate and capture are both reachable.
- The 50-move rule is not applied.

## Changes from the earlier files

The earlier files carried `class16`, `orbit`, `tag`, `side`, `mover_favoured`, `winner_moves`, `moves`, `to_mover_win`, `to_draw`, `to_opp_win` and `opt` (forced), and `moves`, `stay`, `stay_leaf`, `dev`, `dev_d_min`, `dev_d_max`, `stay_with_dev` and `rests_on_leaf` (choice). None is kept. The move composition is the work of the d-depth metric and not a stored record. The number of raw positions in a class is derivable from the symmetry.
