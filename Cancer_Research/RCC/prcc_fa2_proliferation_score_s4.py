"""
PRCC FA-2 — Combined proliferation-depth score
Hypothesis: TI_FA2 × MKI67 outperforms either alone in Type 2 OS.
This is a new prediction, locked here before running.

Locked prediction (2026-03-07):
  P1: TI_FA2 * MKI67_norm predicts OS in Type 2 better than
      MKI67 alone (C-index or log-rank p).
  P2: CDK4 * MKI67_norm predicts Type 2 OS (HR > 2, p < 0.01).
  P3: EZH2 survives adjustment for MKI67 in multivariate Cox
      (p < 0.05), indicating independent chromatin-layer signal.

OrganismCore | 2026-03-07
"""

import os, gzip, warnings
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")

try:
    from lifelines import KaplanMeierFitter, CoxPHFitter
    from lifelines.statistics import logrank_test
    from lifelines.utils import concordance_index
except ImportError:
    import subprocess, sys
    subprocess.run([sys.executable,"-m","pip","install",
                    "lifelines","--quiet"])
    from lifelines import KaplanMeierFitter, CoxPHFitter
    from lifelines.statistics import logrank_test
    from lifelines.utils import concordance_index

try:
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    PLOT_OK = True
except ImportError:
    PLOT_OK = False

RCC_DIR   = os.path.dirname(os.path.abspath(__file__))
PRCC_BASE = os.path.join(RCC_DIR,"PRCC","prcc_false_attractor")
PRCC_S2   = os.path.join(PRCC_BASE,"results_s2")
PRCC_S5   = os.path.join(PRCC_BASE,"results_s5")
S4_OUT    = os.path.join(RCC_DIR,"CCRCC",
                          "ccrcc_false_attractor","results_s4")

KIRP_EXPR    = os.path.join(PRCC_BASE,"TCGA_KIRP_HiSeqV2.gz")
KIRP_SURV    = os.path.join(PRCC_BASE,"KIRP_survival.txt")
GDC_SUBTYPES = os.path.join(PRCC_BASE,"KIRP_GDC_subtypes.tsv")
TI_FA2_FILE  = os.path.join(PRCC_S5,"TI_FA2.csv")

OUT_CSV = os.path.join(S4_OUT,"prcc_fa2_prolif_score.csv")
OUT_FIG = os.path.join(S4_OUT,"prcc_fa2_prolif_figure.png")
OUT_LOG = os.path.join(S4_OUT,"prcc_fa2_prolif_log.txt")

_log = []
def log(m=""): _log.append(m); print(m)
def flush(): open(OUT_LOG,"w").write("\n".join(_log))
def b12(s): return "-".join(str(s).split("-")[:3])
def fmt(x,d=4):
    return "nan" if (x is None or
                     (isinstance(x,float) and np.isnan(x))) \
           else f"{x:.{d}f}"
def fmt_p(p):
    if p is None or np.isnan(p): return "nan"
    return f"{p:.2e}" if p < 0.0001 else f"{p:.4f}"
def median_s(kmf):
    m = kmf.median_survival_time_
    return "NR" if (m is None or np.isnan(float(m))) \
           else f"{float(m):.0f}d"

def cox1(df, dur, evt, cov):
    sub = df[[dur,evt,cov]].dropna()
    if len(sub)<10 or sub[evt].sum()<3: return (np.nan,)*5
    try:
        cph = CoxPHFitter()
        cph.fit(sub,duration_col=dur,event_col=evt,
                show_progress=False)
        r   = cph.summary.loc[cov]
        hr  = float(np.exp(r["coef"]))
        cil = float(np.exp(r["coef lower 95%"]))
        ciu = float(np.exp(r["coef upper 95%"]))
        p   = float(r["p"])
        c   = concordance_index(sub[dur],-sub[cov],sub[evt])
        return hr,cil,ciu,p,c
    except Exception as e:
        log(f"  Cox error: {e}")
        return (np.nan,)*5

def cox_multi(df, dur, evt, covs):
    sub = df[[dur,evt]+covs].dropna()
    if len(sub)<20 or sub[evt].sum()<5: return None
    try:
        cph = CoxPHFitter()
        cph.fit(sub,duration_col=dur,event_col=evt,
                show_progress=False)
        return cph.summary, len(sub)
    except Exception as e:
        log(f"  Multi-Cox error: {e}")
        return None

def norm01(s):
    mn,mx = s.min(),s.max()
    return (s-mn)/(mx-mn) if mx>mn else pd.Series(0.0,index=s.index)

# ═════════════════���════════════════════════════════════════════════
log("PRCC FA-2 PROLIFERATION-DEPTH COMBINED SCORE")
log("OrganismCore | 2026-03-07")
log("Predictions locked before run (see header)")
log("="*60)

# ── Load subtypes ─────────────────────────────────────────────────
gdc = pd.read_csv(GDC_SUBTYPES,sep="\t")
gdc["b12"] = gdc["SAMPLE_ID"].apply(b12)
def nt(s):
    s = str(s).lower()
    if "1" in s: return "Type1"
    if "2" in s: return "Type2"
    return "Unknown"
gdc["type"] = gdc["SUBTYPE"].apply(nt)
gdc = gdc.set_index("b12")[["type"]]

# ── Load survival ─────────────────────────────────────────────────
surv = pd.read_csv(KIRP_SURV,sep="\t")
tc = next(c for c in surv.columns if "os.time" in c.lower())
ec = next(c for c in surv.columns if c.lower() in ("os","_os"))
sc = surv.columns[0]
surv = surv[[sc,tc,ec]].copy()
surv.columns = ["sample","os_time","os_event"]
surv["os_time"]  = pd.to_numeric(surv["os_time"],  errors="coerce")
surv["os_event"] = pd.to_numeric(surv["os_event"], errors="coerce")
surv = surv.dropna()
surv["b12"] = surv["sample"].apply(b12)
surv = surv.set_index("b12")

# ── Load TI_FA2 ───────────────────────────────────────────────────
ti2 = pd.read_csv(TI_FA2_FILE)[["sample_id","TI_FA2"]].copy()
ti2["b12"] = ti2["sample_id"].apply(b12)
ti2 = ti2.set_index("b12")

# ── Load expression ───────────────────────────────────────────────
GENES = ["MKI67","CDK4","EZH2","RUNX1","LAMC2","SLC7A9",
         "TOP2A","CCNE1","FH","OGDHL","CDKN2A","CDKN1A"]
with gzip.open(KIRP_EXPR,"rt") as f:
    raw = pd.read_csv(f,sep="\t",index_col=0)
avail = [g for g in GENES if g in raw.index]
miss  = [g for g in GENES if g not in raw.index]
if miss: log(f"Missing genes: {miss}")
tcols = [c for c in raw.columns
         if len(c.split("-"))>=4
         and c.split("-")[3][:2].isdigit()
         and 1<=int(c.split("-")[3][:2])<=9]
expr = raw.loc[avail,tcols].T
expr.index = [b12(s) for s in expr.index]

# ── Build master ──────────────────────────────────────────────────
master = surv.join(gdc,   how="left")
master = master.join(ti2, how="left")
master = master.join(expr, how="left")
master = master.dropna(subset=["os_time","os_event"])

t2 = master[master["type"]=="Type2"].copy()
log(f"\nType 2 frame: n={len(t2)}  "
    f"events={int(t2['os_event'].sum())}")

# ── Build composite scores ────────────────────────────────────────
log("\nBuilding composite scores...")

# Normalise all to [0,1] within Type 2
for g in ["MKI67","CDK4","EZH2","TI_FA2"]:
    if g in t2.columns:
        t2[f"{g}_n"] = norm01(t2[g].fillna(t2[g].median()))

# Score 1: TI_FA2 × MKI67 (multiplicative — both high = worst)
if "TI_FA2_n" in t2.columns and "MKI67_n" in t2.columns:
    t2["TI_MKI_product"] = t2["TI_FA2_n"] * t2["MKI67_n"]
    log("  TI_FA2 × MKI67 product: built")

# Score 2: TI_FA2 + MKI67 (additive sum)
if "TI_FA2_n" in t2.columns and "MKI67_n" in t2.columns:
    t2["TI_MKI_sum"] = t2["TI_FA2_n"] + t2["MKI67_n"]
    log("  TI_FA2 + MKI67 sum: built")

# Score 3: CDK4 × MKI67
if "CDK4_n" in t2.columns and "MKI67_n" in t2.columns:
    t2["CDK4_MKI_product"] = t2["CDK4_n"] * t2["MKI67_n"]
    log("  CDK4 × MKI67 product: built")

# Score 4: EZH2 + MKI67
if "EZH2_n" in t2.columns and "MKI67_n" in t2.columns:
    t2["EZH2_MKI_sum"] = t2["EZH2_n"] + t2["MKI67_n"]
    log("  EZH2 + MKI67 sum: built")

# ── Cox for each score vs individual genes ────────────────────────
log("\n" + "="*60)
log("COX — COMPOSITE SCORES vs INDIVIDUAL GENES")
log("Locked predictions: see header")
log("="*60)
log(f"\n  {'Predictor':<22} {'HR':>7} {'CI_lo':>7} {'CI_hi':>7} "
    f"{'p':>10}  {'C':>6}  n")
log("  " + "-"*65)

score_rows = []
test_vars = [
    ("MKI67",            "individual — benchmark"),
    ("CDK4",             "individual"),
    ("EZH2",             "individual"),
    ("TI_FA2",           "individual (not confirmed solo)"),
    ("TI_MKI_product",   "COMPOSITE: TI_FA2 × MKI67  [P1]"),
    ("TI_MKI_sum",       "COMPOSITE: TI_FA2 + MKI67"),
    ("CDK4_MKI_product", "COMPOSITE: CDK4 × MKI67    [P2]"),
    ("EZH2_MKI_sum",     "COMPOSITE: EZH2 + MKI67"),
]

for var, label in test_vars:
    if var not in t2.columns:
        log(f"  {var:<22} MISSING")
        continue
    hr,cil,ciu,p,c = cox1(t2,"os_time","os_event",var)
    log(f"  {var:<22} {fmt(hr):>7} {fmt(cil):>7} {fmt(ciu):>7} "
        f"{fmt_p(p):>10}  {fmt(c):>6}  {(~t2[[var,'os_time','os_event']].isnull().any(axis=1)).sum()}")
    score_rows.append({"predictor": var, "label": label,
                        "hr": hr, "ci_lo": cil, "ci_hi": ciu,
                        "p": p, "c_index": c})

# ── Verdict: P1 ───────────────────────────────────────────────────
log("\n" + "="*60)
log("P1 VERDICT: TI_FA2 × MKI67 vs MKI67 alone")
log("="*60)
c_mki = next((r["c_index"] for r in score_rows
              if r["predictor"]=="MKI67"), np.nan)
c_prod = next((r["c_index"] for r in score_rows
               if r["predictor"]=="TI_MKI_product"), np.nan)
if not np.isnan(c_mki) and not np.isnan(c_prod):
    v1 = "CONFIRMED ✓" if c_prod > c_mki else "NOT CONFIRMED ✗"
    log(f"  MKI67 alone C={fmt(c_mki)}")
    log(f"  TI×MKI67    C={fmt(c_prod)}")
    log(f"  P1: {v1}")

# ── Verdict: P3 — EZH2 survives MKI67 adjustment ─────────────────
log("\n" + "="*60)
log("P3 VERDICT: EZH2 independent of MKI67 (multivariate)")
log("="*60)
covs3 = [c for c in ["EZH2","MKI67","CDK4","TI_FA2"]
         if c in t2.columns]
res3 = cox_multi(t2,"os_time","os_event",covs3)
if res3 is not None:
    summ3, n3 = res3
    log(f"  Multivariate Cox  n={n3}")
    log(f"  {'Covariate':<12} {'HR':>7}  {'p':>10}")
    log("  " + "-"*35)
    for cov in covs3:
        if cov in summ3.index:
            hr_ = float(np.exp(summ3.loc[cov,"coef"]))
            p_  = float(summ3.loc[cov,"p"])
            log(f"  {cov:<12} {fmt(hr_):>7}  {fmt_p(p_):>10}")
            if cov == "EZH2":
                v3 = "CONFIRMED ✓" if p_ < 0.05 \
                     else "NOT CONFIRMED ✗"
                log(f"  P3 (EZH2 independent): {v3}")

# ── KM: TI×MKI67 tertile split ────────────────────────────────────
log("\n" + "="*60)
log("KM: TI×MKI67 PRODUCT — TERTILE SPLIT")
log("="*60)
if "TI_MKI_product" in t2.columns:
    t2v = t2.dropna(subset=["TI_MKI_product","os_time","os_event"]).copy()
    try:
        t2v["tertile"] = pd.qcut(t2v["TI_MKI_product"], 3,
                                  labels=["Low","Mid","High"])
        for tert in ["Low","Mid","High"]:
            g = t2v[t2v["tertile"]==tert]
            kmf = KaplanMeierFitter().fit(
                g["os_time"],g["os_event"],label=tert)
            log(f"  {tert}: n={len(g)}  "
                f"events={int(g['os_event'].sum())}  "
                f"median={median_s(kmf)}")
        lo = t2v[t2v["tertile"]=="Low"]
        hi = t2v[t2v["tertile"]=="High"]
        if len(lo)>=5 and len(hi)>=5:
            lr_ = logrank_test(
                lo["os_time"].values, hi["os_time"].values,
                event_observed_A=lo["os_event"].values,
                event_observed_B=hi["os_event"].values)
            log(f"  Low vs High log-rank p = {fmt_p(lr_.p_value)}")
    except Exception as e:
        log(f"  Tertile error: {e}")

# ── Figure ────────────────────────────────────────────────────────
if PLOT_OK:
    fig = plt.figure(figsize=(14,9))
    gs  = gridspec.GridSpec(2,3,hspace=0.38,wspace=0.32)

    def style(ax,title):
        ax.set_title(title,fontsize=9,fontweight="bold",pad=4)
        ax.set_xlabel("Time (days)",fontsize=8)
        ax.set_ylabel("Survival prob.",fontsize=8)
        ax.set_ylim(-0.03,1.07)
        ax.tick_params(labelsize=7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(fontsize=7,loc="upper right")

    PAL = {"Low":"#2ecc71","Mid":"#f39c12","High":"#e74c3c"}

    # ── Panels A–C: KM for top 3 composite scores ─────────────────
    combos = [
        ("TI_MKI_product", "A.  TI_FA2 × MKI67"),
        ("CDK4_MKI_product","B.  CDK4 × MKI67"),
        ("EZH2_MKI_sum",   "C.  EZH2 + MKI67"),
    ]
    for idx,(score,title) in enumerate(combos):
        ax = fig.add_subplot(gs[0,idx])
        if score not in t2.columns:
            ax.set_title(f"{score} missing",fontsize=8); continue
        sub = t2.dropna(
            subset=[score,"os_time","os_event"]).copy()
        if len(sub)<12:
            ax.set_title(f"{score}: n too small",fontsize=8)
            continue
        try:
            sub["tert"] = pd.qcut(sub[score],3,
                                   labels=["Low","Mid","High"])
        except Exception:
            ax.set_title(f"{score}: qcut failed",fontsize=8)
            continue
        for tert in ["Low","Mid","High"]:
            g = sub[sub["tert"]==tert]
            if len(g)>=3:
                KaplanMeierFitter().fit(
                    g["os_time"],g["os_event"],label=tert
                ).plot_survival_function(
                    ax=ax,ci_show=False,
                    color=PAL[tert],linewidth=1.5)
        lo = sub[sub["tert"]=="Low"]
        hi = sub[sub["tert"]=="High"]
        if len(lo)>=5 and len(hi)>=5:
            lr_ = logrank_test(
                lo["os_time"].values,hi["os_time"].values,
                event_observed_A=lo["os_event"].values,
                event_observed_B=hi["os_event"].values)
            ax.text(0.05,0.08,f"p={fmt_p(lr_.p_value)}",
                    transform=ax.transAxes,fontsize=8)
        style(ax,title)

    # ── Panel D: C-index comparison bar ──────────────────────────
    ax = fig.add_subplot(gs[1,0])
    c_vals = [(r["predictor"],r["c_index"])
              for r in score_rows
              if not np.isnan(r["c_index"])]
    c_vals.sort(key=lambda x: x[1])
    labels_ = [v[0] for v in c_vals]
    vals_   = [v[1] for v in c_vals]
    colors_ = ["#e74c3c" if "product" in l or "sum" in l
               else "#95a5a6" for l in labels_]
    ax.barh(range(len(labels_)),vals_,color=colors_,height=0.65)
    ax.set_yticks(range(len(labels_)))
    ax.set_yticklabels(labels_,fontsize=7)
    ax.axvline(0.5,color="k",linewidth=0.8,linestyle="--",alpha=0.5)
    ax.set_xlabel("C-index",fontsize=8)
    ax.set_title("D.  C-index comparison (Type 2)",
                 fontsize=9,fontweight="bold",pad=4)
    ax.tick_params(labelsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # ── Panel E: MKI67 median split (benchmark) ───────────────────
    ax = fig.add_subplot(gs[1,1])
    if "MKI67" in t2.columns:
        sub_m = t2[["MKI67","os_time","os_event"]].dropna().copy()
        med_m = sub_m["MKI67"].median()
        sub_m["grp"] = np.where(sub_m["MKI67"]>=med_m,"High","Low")
        for lbl,col in [("High","#e74c3c"),("Low","#2ecc71")]:
            g = sub_m[sub_m["grp"]==lbl]
            if len(g)>=3:
                KaplanMeierFitter().fit(
                    g["os_time"],g["os_event"],
                    label=f"MKI67 {lbl}"
                ).plot_survival_function(
                    ax=ax,ci_show=True,color=col,linewidth=1.5)
        lo_ = sub_m[sub_m["grp"]=="Low"]
        hi_ = sub_m[sub_m["grp"]=="High"]
        if len(lo_)>=5 and len(hi_)>=5:
            lr_m = logrank_test(
                lo_["os_time"].values,hi_["os_time"].values,
                event_observed_A=lo_["os_event"].values,
                event_observed_B=hi_["os_event"].values)
            ax.text(0.05,0.08,f"p={fmt_p(lr_m.p_value)}",
                    transform=ax.transAxes,fontsize=8)
    style(ax,"E.  PRCC Type 2: MKI67 (benchmark)")

    # ── Panel F: EZH2 median split in Type 2 ─────────────────────
    ax = fig.add_subplot(gs[1,2])
    if "EZH2" in t2.columns:
        sub_e = t2[["EZH2","os_time","os_event"]].dropna().copy()
        med_e = sub_e["EZH2"].median()
        sub_e["grp"] = np.where(sub_e["EZH2"]>=med_e,"High","Low")
        for lbl,col in [("High","#e74c3c"),("Low","#2ecc71")]:
            g = sub_e[sub_e["grp"]==lbl]
            if len(g)>=3:
                KaplanMeierFitter().fit(
                    g["os_time"],g["os_event"],
                    label=f"EZH2 {lbl}"
                ).plot_survival_function(
                    ax=ax,ci_show=True,color=col,linewidth=1.5)
        lo_ = sub_e[sub_e["grp"]=="Low"]
        hi_ = sub_e[sub_e["grp"]=="High"]
        if len(lo_)>=5 and len(hi_)>=5:
            lr_e = logrank_test(
                lo_["os_time"].values,hi_["os_time"].values,
                event_observed_A=lo_["os_event"].values,
                event_observed_B=hi_["os_event"].values)
            ax.text(0.05,0.08,f"p={fmt_p(lr_e.p_value)}",
                    transform=ax.transAxes,fontsize=8)
    style(ax,"F.  PRCC Type 2: EZH2")

    fig.suptitle(
        "PRCC Type 2 — Proliferation-Depth Combined Scores\n"
        "OrganismCore — Eric Robert Lawson | 2026-03-07"
        " | Document 95f-S4-addendum",
        fontsize=10,fontweight="bold",y=0.998)
    plt.savefig(OUT_FIG,dpi=150,bbox_inches="tight")
    plt.close()
    log(f"\nFigure: {OUT_FIG}")

# ── Save ──────────────────────────────────────────────────────────
pd.DataFrame(score_rows).to_csv(OUT_CSV,index=False)
log(f"CSV: {OUT_CSV}")
flush()
log("COMPLETE.")
