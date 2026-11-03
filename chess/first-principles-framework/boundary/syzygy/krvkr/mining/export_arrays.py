# Writes the table-derived arrays that the verification harness of krvkr_calc.cpp reads.
# Run krvkr_mine1_enumerate_classify.py first (it writes pass1.npz and pass2.npz into OUT).
import numpy as np
OUT = "/home/claude/syz/out"          # change to your own folder
p1 = np.load(OUT + "/pass1.npz"); p2 = np.load(OUT + "/pass2.npz")
p1["reps"].astype(np.int32).tofile("reps.bin")
p1["valid"].astype(np.uint8).tofile("valid.bin")
p1["val"].astype(np.int8).tofile("val.bin")
for k, f in (("nlose", "nlose.bin"), ("ndraw_in", "ndi.bin"), ("ndraw_cap", "ndc.bin"), ("nmoves", "nm.bin"), ("nwin", "nwin.bin")):
    p2[k].astype(np.int16).tofile(f)
