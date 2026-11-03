import numpy as np, scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
import os
TABLES = os.environ.get("KRVKR_TABLES", "/home/claude/syz")      # the folder with KRvK and KRvKR .rtbw/.rtbz
OUT = os.environ.get("KRVKR_OUT", TABLES + "/out")           # where the mined arrays are written and read
p1 = np.load(OUT + "/pass1.npz"); p2 = np.load(OUT + "/pass2.npz")
valid = p1["valid"]; val = p1["val"]; rawwdl = p1["rawwdl"]
nm, nl, ndi, ndc, nds, nw = [p2[k] for k in ("nmoves", "nlose", "ndraw_in", "ndraw_cap", "ndraw_stale", "nwin")]
inc = p2["incheck"]; src = p2["src"]; dst = p2["dst"]
M = len(valid)
print("valid folded positions", valid.sum())
for v, name in ((1, "mover wins"), (0, "drawn"), (-1, "mover loses")):
    print(name, int(((val == v) & valid).sum()))
print("raw wdl values 2/-2/1/-1 (cursed/blessed):", {int(k): int(((rawwdl == k) & valid).sum()) for k in (-2, -1, 0, 1, 2)})

D = valid & (val == 0)
nD = int(D.sum())
print("\n== drawn positions:", nD, " share of all valid:", round(nD / valid.sum(), 4))
print("drawn positions with a winning move (must be 0):", int((nw[D] > 0).sum()))
stal = D & (nm == 0)
print("stalemates:", int(stal.sum()))
print("drawn, side to move in check:", int((D & inc).sum()))

# step 4/6 counts
dr_total = ndi + ndc
print("\ndrawing moves per drawn position (in-material + capture): width distribution")
w = dr_total[D & (nm > 0)]
u, c = np.unique(w, return_counts=True)
for a, b in list(zip(u, c))[:12]: print("  width", int(a), ":", int(b))
print("  max width", int(w.max()), " mean", round(float(w.mean()), 2))
print("drawn positions with >=1 losing move (side-to-move deviation exists):", int((D & (nl > 0)).sum()),
      round(int((D & (nl > 0)).sum()) / nD, 4))
print("drawn with NO losing move (interior):", int((D & (nl == 0) & (nm > 0)).sum()))
print("drawn with exactly 1 drawing move:", int((D & (dr_total == 1)).sum()))

# graph
S = sp.csr_matrix((np.ones(len(src), np.int8), (src, dst)), shape=(M, M))
S.data[:] = 1
ncs, lab = connected_components(S, directed=True, connection="strong")
ncw, labw = connected_components(S, directed=True, connection="weak")
size = np.bincount(lab, minlength=ncs)
selfloop = np.zeros(M, bool); selfloop[src[src == dst]] = True
oncycle = (size[lab] > 1) | selfloop
print("\n== graph: edges", len(src))
print("SCC count among drawn nodes:", len(np.unique(lab[D])))
print("drawn on a cycle:", int((D & oncycle).sum()), " not on a cycle:", int((D & ~oncycle).sum()))
sc = size[lab[D]]
print("largest SCC size:", int(size[lab[D]].max()))
big = np.unique(lab[D & oncycle], return_counts=True)
print("number of non-trivial SCCs:", len(np.unique(lab[D & oncycle])))
wc = np.unique(labw[D], return_counts=True)
print("weak components among drawn:", len(wc[0]), " largest:", int(wc[1].max()))

# reachability flags by condensation DAG
# terminal draw: stalemate, or has a drawing capture exit
term = (D & (nm == 0)) | (D & (ndc > 0))
# condensation
cs = lab[src]; cd = lab[dst]
mask = cs != cd
cs, cd = cs[mask], cd[mask]
# component-level flags
cyc_c = np.zeros(ncs, bool); cyc_c[lab[D & oncycle]] = True
term_c = np.zeros(ncs, bool)
np.logical_or.at(term_c, lab[term], True)
# topological processing: Tarjan labels in scipy are not guaranteed topological, so iterate to a fixed point
reach_cyc = cyc_c.copy(); reach_term = term_c.copy()
changed = True; it = 0
while changed:
    it += 1
    old = (reach_cyc.copy(), reach_term.copy())
    np.logical_or.at(reach_cyc, cs, reach_cyc[cd])
    np.logical_or.at(reach_term, cs, reach_term[cd])
    changed = not (np.array_equal(old[0], reach_cyc) and np.array_equal(old[1], reach_term))
print("fixed-point iterations:", it)
rc = reach_cyc[lab]; rt = reach_term[lab]
cls = np.full(M, -1)
cls[D & stal] = 0   # stalemate now
cls[D & ~stal & rt & ~rc] = 1  # terminal only
cls[D & ~stal & rc & ~rt] = 2  # cycle only
cls[D & ~stal & rc & rt] = 3   # mixed
cls[D & ~stal & ~rc & ~rt] = 4 # neither (should not occur)
names = ["stalemate now", "terminal only", "cycle only", "mixed", "neither"]
print("\n== realization classes of drawn positions")
for k in range(5):
    print("  ", names[k], int((cls == k).sum()))
# among drawn that are not in-check etc.
print("\n== boundary criteria (main section 14 items 7, 8), reading A: mover has a losing move")
bA = D & (nl > 0)
print("reading A (side to move has a losing deviation):", int(bA.sum()), round(int(bA.sum()) / nD, 4), "of drawn")
# reading B: mover has a losing deviation AND for some drawing move, the opponent at the child has a losing deviation
opp_dev = np.zeros(M, bool)
# edge-based: child has nlose>0
e_ok = nl[dst] > 0
np.logical_or.at(opp_dev, src, e_ok)
# captures lead out of material; count them by existence of drawing capture exit
bB = D & (nl > 0) & opp_dev
print("reading B (mover has one, and some drawing in-material move reaches a position where the opponent has one):",
      int(bB.sum()), round(int(bB.sum()) / nD, 4), "of drawn")
# reading C: all drawing in-material moves reach positions where the opponent has a losing deviation
cnt_ok = np.zeros(M, np.int32); np.add.at(cnt_ok, src, e_ok.astype(np.int32))
bC = D & (nl > 0) & (ndi > 0) & (cnt_ok == ndi)
print("reading C (as B, for every drawing in-material move):", int(bC.sum()), round(int(bC.sum()) / nD, 4))
print("share of ALL valid positions that are boundary, reading A:", round(int(bA.sum()) / valid.sum(), 4))
print("share of ALL valid positions that are boundary, reading B:", round(int(bB.sum()) / valid.sum(), 4))
np.savez_compressed(OUT + "/analysis.npz", lab=lab, cls=cls, oncycle=oncycle, bA=bA, bB=bB)
