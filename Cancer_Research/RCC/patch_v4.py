"""
PRCC subtype OS patch — Script 4 addendum
Uses KIRP_GDC_subtypes.tsv (already cached) to split
Type 1 vs Type 2 and re-run TI analyses per subtype.

Run from: /Users/ericlawson/cancer/RCC/
OrganismCore | 2026-03-07
"""

import os, gzip, warnings
import numpy as np
import pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")

try:
    from lifelines import KaplanMeierFitter, CoxPHFitter
    from lifelines.statistics import logrank_test
    from lifelines.utils import concordance_index
except ImportError:
    import subprocess, sys
    subprocess.run([sys.executable, "-m", "pip", "install",
                    "lifelines", "--quiet"])
    from lifelines import KaplanMeierFitter, CoxPHFitter
    from lifelines.statistics import logrank_test
    from lifelines.utils import concordance_index

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    PLOT_OK = True
except ImportError:
    PLOT_OK = False

# ── Paths ─────────────────────────────────────────────────────────
RCC_DIR   = os.path.dirname(os.path.abspath(__file__))
PRCC_BASE = os.path.join(RCC_DIR, "PRCC", "prcc_false_attractor")
PRCC_S2   = os.path.join(PRCC_BASE, "results_s2")
PRCC_S5   = os.path.join(PRCC_BASE, "results_s5")
S4_OUT    = os.path.join(RCC_DIR, "CCRCC",
                          "ccrcc_false_attractor", "results_s4")

KIRP_EXPR    = os.path.join(PRCC_BASE, "TCGA_KIRP_HiSeqV2.gz")
KIRP_SURV    = os.path.join(PRCC_BASE, "KIRP_survival.txt")
GDC_SUBTYPES = os.path.join(PRCC_BASE, "KIRP_GDC_subtypes.tsv")
TI_FA1_FILE  = os.path.join(PRCC_S2, "transition_index.csv")
TI_FA2_FILE  = os.path.join(PRCC_S5, "TI_FA2.csv")

OUT_CSV  = os.path.join(S4_OUT, "prcc_subtype_os_patch.csv")
OUT_FIG  = os.path.join(S4_OUT, "prcc_subtype_figure.png")
OUT_LOG  = os.path.join(S4_OUT, "prcc_subtype_patch_log.txt")

_log = []
def log(m=""): _log.append(m); print(m)
def flush(): open(OUT_LOG,"w").write("\n".join(_log))

def b12(s): return "-".join(str(s).split("-")[:3])
def fmt(x, d=4):
    return "nan" if (x is None or
                     (isinstance(x,float) and np.isnan(x))) \
           else f"{x:.{d}f}"
def fmt_p(p):
    if p is None or np.isnan(p): return "nan"
    if p < 1e-100: return "<1e-100"
    if p < 0.0001: return f"{p:.2e}"
    return f"{p:.4f}"
def median_s(kmf):
    m = kmf.median_survival_time_
    return "NR" if (m is None or np.isnan(float(m))) \
           else f"{float(m):.0f}d"

def cox1(df, dur, evt, cov):
    sub = df[[dur,evt,cov]].dropna()
    if len(sub) < 10 or sub[evt].sum() < 3:
        return (np.nan,)*5
    try:
        cph = CoxPHFitter()
        cph.fit(sub, duration_col=dur, event_col=evt,
                show_progress=False)
        r   = cph.summary.loc[cov]
        hr  = float(np.exp(r["coef"]))
        cil = float(np.exp(r["coef lower 95%"]))
        ciu = float(np.exp(r["coef upper 95%"]))
        p   = float(r["p"])
        c   = concordance_index(sub[dur], -sub[cov], sub[evt])
        return hr, cil, ciu, p, c
    except Exception as e:
        log(f"  Cox error: {e}")
        return (np.nan,)*5

def km_split(df, dur, evt, cov_col, label_hi, label_lo):
    a = df[df[cov_col]==label_hi]
    b = df[df[cov_col]==label_lo]
    if len(a)<5 or len(b)<5: return None, None, np.nan
    kh = KaplanMeierFitter().fit(a[dur], a[evt], label=label_hi)
    kl = KaplanMeierFitter().fit(b[dur], b[evt], label=label_lo)
    lr = logrank_test(a[dur].values, b[dur].values,
                      event_observed_A=a[evt].values,
                      event_observed_B=b[evt].values)
    return kh, kl, lr.p_value

# ════════════════════════════════��═════════════════════════════════
log("PRCC SUBTYPE OS PATCH")
log("OrganismCore | 2026-03-07")
log("=" * 60)

# ── 1. Load GDC subtypes ──────────────────────────────────────────
log("\n1. GDC SUBTYPES")
log(f"   {GDC_SUBTYPES}")
if not os.path.exists(GDC_SUBTYPES):
    log("   FATAL: KIRP_GDC_subtypes.tsv not found.")
    log("   Expected at: " + GDC_SUBTYPES)
    flush(); raise SystemExit(1)

gdc = pd.read_csv(GDC_SUBTYPES, sep="\t")
log(f"   cols: {list(gdc.columns)}")
log(f"   rows: {len(gdc)}")

# Find sample + subtype columns
# GDC file typically has 'sample' or 'id' + 'Subtype' or 'Molecular_subtype'
sample_col  = gdc.columns[0]
subtype_col = next((c for c in gdc.columns
                    if any(k in c.lower()
                           for k in ["subtype","type","histol",
                                     "molecular","papillary"])),
                   None)
if subtype_col is None:
    log("   Columns available:")
    for c in gdc.columns: log(f"     {c}")
    log("   Cannot find subtype column — using all columns above.")
    subtype_col = gdc.columns[1]   # fallback: second column

log(f"   Sample col : {sample_col}")
log(f"   Subtype col: {subtype_col}")
log(f"   Values: {dict(gdc[subtype_col].value_counts())}")

gdc["b12"] = gdc[sample_col].apply(b12)
gdc = gdc[["b12", subtype_col]].rename(
    columns={subtype_col: "subtype"}).set_index("b12")

# Normalise subtype labels → "Type1" / "Type2" / "CIMP"
def normalise_type(s):
    s = str(s).lower().strip()
    if any(k in s for k in ["type1", "type 1", "papillary 1",
                             "prcc1", "prcc type 1", "1"]):
        return "Type1"
    if any(k in s for k in ["type2", "type 2", "papillary 2",
                             "prcc2", "prcc type 2", "2"]):
        return "Type2"
    if "cimp" in s or "hlrcc" in s:
        return "CIMP"
    return "Unknown"

gdc["type"] = gdc["subtype"].apply(normalise_type)
log(f"   Normalised: {dict(gdc['type'].value_counts())}")

# ── 2. Load survival ──────────────────────────────────────────────
log("\n2. SURVIVAL")
surv = pd.read_csv(KIRP_SURV, sep="\t")
tc   = next((c for c in surv.columns if "os.time" in c.lower()), None)
ec   = next((c for c in surv.columns if c.lower() in ("os","_os")), None)
sc   = surv.columns[0]
surv = surv[[sc, tc, ec]].copy()
surv.columns = ["sample", "os_time", "os_event"]
surv["os_time"]  = pd.to_numeric(surv["os_time"],  errors="coerce")
surv["os_event"] = pd.to_numeric(surv["os_event"], errors="coerce")
surv = surv.dropna()
surv["b12"] = surv["sample"].apply(b12)
surv = surv.set_index("b12")
log(f"   n={len(surv)}  events={int(surv['os_event'].sum())}")

# ── 3. Load TI files ──────────────────────────────────────────────
log("\n3. TI FILES")
ti1 = pd.read_csv(TI_FA1_FILE)[["sample_id","TI"]].copy()
ti1["b12"] = ti1["sample_id"].apply(b12)
ti1 = ti1.set_index("b12")[["TI"]].rename(columns={"TI":"TI_FA1"})
log(f"   TI_FA1: n={len(ti1)}  range={ti1['TI_FA1'].min():.4f}–"
    f"{ti1['TI_FA1'].max():.4f}")

ti2 = pd.read_csv(TI_FA2_FILE)[["sample_id","TI_FA2"]].copy()
ti2["b12"] = ti2["sample_id"].apply(b12)
ti2 = ti2.set_index("b12")
log(f"   TI_FA2: n={len(ti2)}  range={ti2['TI_FA2'].min():.4f}–"
    f"{ti2['TI_FA2'].max():.4f}")

# ── 4. Load expression ────────────────────────────────────────────
log("\n4. EXPRESSION (key genes)")
KEY_GENES = ["KRT19","SLC22A6","ERBB2","FH","OGDHL",
             "EZH2","CDK4","MKI67","LAMC2","SLC7A9","RUNX1"]
with gzip.open(KIRP_EXPR,"rt") as f:
    raw = pd.read_csv(f, sep="\t", index_col=0)
avail = [g for g in KEY_GENES if g in raw.index]
miss  = [g for g in KEY_GENES if g not in raw.index]
if miss: log(f"   Missing: {miss}")
tumour_cols = [c for c in raw.columns
               if len(c.split("-"))>=4
               and c.split("-")[3][:2].isdigit()
               and 1 <= int(c.split("-")[3][:2]) <= 9]
expr = raw.loc[avail, tumour_cols].T
expr.index = [b12(s) for s in expr.index]
log(f"   Tumour samples: {len(expr)}")

# ── 5. Build master frame ─────────────────────────────────────────
log("\n5. MASTER FRAME")
master = surv.join(gdc[["type"]], how="left")
master = master.join(ti1,  how="left")
master = master.join(ti2,  how="left")
master = master.join(expr, how="left")
master = master.dropna(subset=["os_time","os_event"])
log(f"   Total: n={len(master)}")
log(f"   Type distribution:")
for t, g in master.groupby("type"):
    log(f"     {t}: n={len(g)}  events={int(g['os_event'].sum())}")

# ── 6. Type 1 vs Type 2 OS ────────────────────────────────────────
log("\n" + "="*60)
log("ANALYSIS 1 — TYPE 1 vs TYPE 2 OS")
log("="*60)

t1 = master[master["type"]=="Type1"]
t2 = master[master["type"]=="Type2"]
log(f"  Type1 n={len(t1)}  events={int(t1['os_event'].sum())}")
log(f"  Type2 n={len(t2)}  events={int(t2['os_event'].sum())}")

results = []

if len(t1)>=5 and len(t2)>=5:
    lr = logrank_test(t1["os_time"].values, t2["os_time"].values,
                      event_observed_A=t1["os_event"].values,
                      event_observed_B=t2["os_event"].values)
    log(f"  Log-rank p = {fmt_p(lr.p_value)}")
    kmf1 = KaplanMeierFitter().fit(t1["os_time"], t1["os_event"],
                                    label="Type1")
    kmf2 = KaplanMeierFitter().fit(t2["os_time"], t2["os_event"],
                                    label="Type2")
    log(f"  Type1 median OS: {median_s(kmf1)}")
    log(f"  Type2 median OS: {median_s(kmf2)}")
    v = ("CONFIRMED ✓" if lr.p_value < 0.05
         else "NOT CONFIRMED ✗")
    log(f"  S4-P5b (Type2 < Type1): {v}")
    results.append({"analysis":"Type1_vs_Type2",
                    "logrank_p": lr.p_value,
                    "n_t1": len(t1), "n_t2": len(t2),
                    "med_t1": median_s(kmf1),
                    "med_t2": median_s(kmf2),
                    "verdict": v})
else:
    log("  Insufficient n for Type1/Type2 split.")
    log("  Check GDC subtype labels above.")

# ── 7. FA-1 TI vs OS within Type 1 ───────────────────────────────
log("\n" + "="*60)
log("ANALYSIS 2 — FA-1 TI vs OS (TYPE 1 only)")
log("  Prediction: TI_FA1 predicts OS within Type 1")
log("="*60)

t1_ti = t1.dropna(subset=["TI_FA1"])
log(f"  Type1 with TI_FA1: n={len(t1_ti)}  "
    f"events={int(t1_ti['os_event'].sum())}")

if len(t1_ti) >= 10:
    hr1, cil1, ciu1, p1, c1 = cox1(
        t1_ti, "os_time", "os_event", "TI_FA1")
    log(f"  Cox (TI_FA1 in Type1): HR={fmt(hr1)} "
        f"[{fmt(cil1)}–{fmt(ciu1)}]  p={fmt_p(p1)}  C={fmt(c1)}")

    med1 = t1_ti["TI_FA1"].median()
    t1_ti = t1_ti.copy()
    t1_ti["ti1_grp"] = np.where(
        t1_ti["TI_FA1"] >= med1,
        "TI1-High (deep FA1)", "TI1-Low (shallow)")
    kh, kl, pkm = km_split(
        t1_ti, "os_time", "os_event", "ti1_grp",
        "TI1-High (deep FA1)", "TI1-Low (shallow)")
    if kh is not None:
        log(f"  KM split: TI1-High median={median_s(kh)}  "
            f"TI1-Low median={median_s(kl)}")
        log(f"           log-rank p={fmt_p(pkm)}")
        v = ("CONFIRMED ✓" if pkm < 0.05
             else "NOT CONFIRMED ✗")
        log(f"  S4-P5a (FA1 TI in Type1): {v}")
        results.append({"analysis":"FA1_TI_Type1",
                         "hr": hr1, "p": p1, "c": c1,
                         "logrank_p": pkm, "verdict": v})

# ── 8. FA-2 TI vs OS within Type 2 ───────────────────────────────
log("\n" + "="*60)
log("ANALYSIS 3 — FA-2 TI vs OS (TYPE 2 only)")
log("  Prediction: TI_FA2 (LAMC2/SLC7A9) predicts OS in Type 2")
log("="*60)

t2_ti = t2.dropna(subset=["TI_FA2"])
log(f"  Type2 with TI_FA2: n={len(t2_ti)}  "
    f"events={int(t2_ti['os_event'].sum())}")

if len(t2_ti) >= 10:
    hr2, cil2, ciu2, p2, c2 = cox1(
        t2_ti, "os_time", "os_event", "TI_FA2")
    log(f"  Cox (TI_FA2 in Type2): HR={fmt(hr2)} "
        f"[{fmt(cil2)}–{fmt(ciu2)}]  p={fmt_p(p2)}  C={fmt(c2)}")

    med2 = t2_ti["TI_FA2"].median()
    t2_ti = t2_ti.copy()
    t2_ti["ti2_grp"] = np.where(
        t2_ti["TI_FA2"] >= med2,
        "TI2-High (deep FA2)", "TI2-Low (shallow)")
    kh2, kl2, pkm2 = km_split(
        t2_ti, "os_time", "os_event", "ti2_grp",
        "TI2-High (deep FA2)", "TI2-Low (shallow)")
    if kh2 is not None:
        log(f"  KM split: TI2-High median={median_s(kh2)}  "
            f"TI2-Low median={median_s(kl2)}")
        log(f"           log-rank p={fmt_p(pkm2)}")
        v2 = ("CONFIRMED ✓" if pkm2 < 0.05
              else "NOT CONFIRMED ✗")
        log(f"  FA-2 TI in Type2: {v2}")
        results.append({"analysis":"FA2_TI_Type2",
                         "hr": hr2, "p": p2, "c": c2,
                         "logrank_p": pkm2, "verdict": v2})

# ── 9. CDK4 vs OS within Type 2 (strongest PRCC finding) ─────────
log("\n" + "="*60)
log("ANALYSIS 4 — CDK4 vs OS (TYPE 2 only)")
log("  CDK4 was the strongest gene predictor (pooled HR=3.06)")
log("  Expected to be even stronger within Type 2")
log("="*60)

for gene in ["CDK4","MKI67","EZH2","RUNX1","FH","OGDHL"]:
    if gene not in t2.columns:
        continue
    sub = t2[[gene,"os_time","os_event"]].dropna()
    if len(sub)<10 or sub["os_event"].sum()<3:
        continue
    hr_g, cil_g, ciu_g, p_g, c_g = cox1(
        sub, "os_time", "os_event", gene)
    log(f"  {gene:<10} HR={fmt(hr_g)} [{fmt(cil_g)}–{fmt(ciu_g)}]  "
        f"p={fmt_p(p_g)}  C={fmt(c_g)}  n={len(sub)}")
    results.append({"analysis": f"gene_{gene}_Type2",
                     "hr": hr_g, "p": p_g, "c_index": c_g,
                     "n": len(sub),
                     "events": int(sub["os_event"].sum())})

# ── 10. Figure ────────────────────────────────────────────────────
if PLOT_OK and len(t1)>=5 and len(t2)>=5:
    fig = plt.figure(figsize=(14, 10))
    gs  = gridspec.GridSpec(2, 3, hspace=0.38, wspace=0.32)

    axes = [fig.add_subplot(gs[r,c])
            for r in range(2) for c in range(3)]

    def style(ax, title):
        ax.set_title(title, fontsize=9, fontweight="bold", pad=4)
        ax.set_xlabel("Time (days)", fontsize=8)
        ax.set_ylabel("Survival prob.", fontsize=8)
        ax.set_ylim(-0.03,1.07)
        ax.tick_params(labelsize=7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(fontsize=7, loc="upper right")

    # Panel A: Type 1 vs Type 2
    ax = axes[0]
    for t, col, lbl in [(t1,"#2ecc71","Type 1"),
                         (t2,"#e74c3c","Type 2")]:
        if len(t)>=5:
            KaplanMeierFitter().fit(
                t["os_time"],t["os_event"],label=lbl
            ).plot_survival_function(ax=ax, ci_show=True,
                                      color=col, linewidth=1.5)
    if len(t1)>=5 and len(t2)>=5:
        lr_ = logrank_test(
            t1["os_time"].values, t2["os_time"].values,
            event_observed_A=t1["os_event"].values,
            event_observed_B=t2["os_event"].values)
        ax.text(0.05,0.08,f"p={fmt_p(lr_.p_value)}",
                transform=ax.transAxes, fontsize=8)
    style(ax, "A.  PRCC: Type 1 vs Type 2")

    # Panel B: FA-1 TI split in Type 1
    ax = axes[1]
    t1_plot = t1.dropna(subset=["TI_FA1"]).copy()
    if len(t1_plot)>=10:
        m = t1_plot["TI_FA1"].median()
        t1_plot["grp"] = np.where(
            t1_plot["TI_FA1"]>=m, "High","Low")
        for lbl, col in [("High","#e74c3c"),("Low","#2ecc71")]:
            g = t1_plot[t1_plot["grp"]==lbl]
            if len(g)>=3:
                KaplanMeierFitter().fit(
                    g["os_time"],g["os_event"],
                    label=f"TI1-{lbl}"
                ).plot_survival_function(ax=ax,ci_show=True,
                                          color=col,linewidth=1.5)
        lr2_ = logrank_test(
            t1_plot[t1_plot["grp"]=="High"]["os_time"].values,
            t1_plot[t1_plot["grp"]=="Low"]["os_time"].values,
            event_observed_A=t1_plot[
                t1_plot["grp"]=="High"]["os_event"].values,
            event_observed_B=t1_plot[
                t1_plot["grp"]=="Low"]["os_event"].values)
        ax.text(0.05,0.08,f"p={fmt_p(lr2_.p_value)}",
                transform=ax.transAxes, fontsize=8)
    style(ax, "B.  PRCC Type 1: FA-1 TI split")

    # Panel C: FA-2 TI split in Type 2
    ax = axes[2]
    t2_plot = t2.dropna(subset=["TI_FA2"]).copy()
    if len(t2_plot)>=10:
        m2 = t2_plot["TI_FA2"].median()
        t2_plot["grp"] = np.where(
            t2_plot["TI_FA2"]>=m2, "High","Low")
        for lbl, col in [("High","#e74c3c"),("Low","#2ecc71")]:
            g = t2_plot[t2_plot["grp"]==lbl]
            if len(g)>=3:
                KaplanMeierFitter().fit(
                    g["os_time"],g["os_event"],
                    label=f"TI2-{lbl}"
                ).plot_survival_function(ax=ax,ci_show=True,
                                          color=col,linewidth=1.5)
        if len(t2_plot[t2_plot["grp"]=="High"])>=5 and \
           len(t2_plot[t2_plot["grp"]=="Low"])>=5:
            lr3_ = logrank_test(
                t2_plot[t2_plot["grp"]=="High"]["os_time"].values,
                t2_plot[t2_plot["grp"]=="Low"]["os_time"].values,
                event_observed_A=t2_plot[
                    t2_plot["grp"]=="High"]["os_event"].values,
                event_observed_B=t2_plot[
                    t2_plot["grp"]=="Low"]["os_event"].values)
            ax.text(0.05,0.08,f"p={fmt_p(lr3_.p_value)}",
                    transform=ax.transAxes, fontsize=8)
    style(ax, "C.  PRCC Type 2: FA-2 TI split")

    # Panels D-F: CDK4, MKI67, FH in Type 2 (median split KM)
    for idx, gene in enumerate(["CDK4","MKI67","FH"]):
        ax = axes[3+idx]
        if gene not in t2.columns:
            ax.set_title(f"{gene} not available", fontsize=8)
            continue
        sub = t2[[gene,"os_time","os_event"]].dropna().copy()
        if len(sub)<10:
            ax.set_title(f"{gene}: n too small", fontsize=8)
            continue
        med_g = sub[gene].median()
        sub["grp"] = np.where(sub[gene]>=med_g,"High","Low")
        for lbl, col in [("High","#e74c3c"),("Low","#2ecc71")]:
            g = sub[sub["grp"]==lbl]
            if len(g)>=3:
                KaplanMeierFitter().fit(
                    g["os_time"],g["os_event"],
                    label=f"{gene} {lbl}"
                ).plot_survival_function(ax=ax,ci_show=True,
                                          color=col,linewidth=1.5)
        if len(sub[sub["grp"]=="High"])>=5 and \
           len(sub[sub["grp"]=="Low"])>=5:
            lr_g = logrank_test(
                sub[sub["grp"]=="High"]["os_time"].values,
                sub[sub["grp"]=="Low"]["os_time"].values,
                event_observed_A=sub[
                    sub["grp"]=="High"]["os_event"].values,
                event_observed_B=sub[
                    sub["grp"]=="Low"]["os_event"].values)
            ax.text(0.05,0.08,f"p={fmt_p(lr_g.p_value)}",
                    transform=ax.transAxes,fontsize=8)
        style(ax, f"{chr(68+idx)}.  PRCC Type 2: {gene}")

    fig.suptitle(
        "PRCC — Type 1 / Type 2 Subtype OS Analysis\n"
        "OrganismCore — Eric Robert Lawson | 2026-03-07"
        " | Document 95f-S4-patch",
        fontsize=10, fontweight="bold", y=0.998)
    plt.savefig(OUT_FIG, dpi=150, bbox_inches="tight")
    plt.close()
    log(f"\n  Figure saved: {OUT_FIG}")

# ── 11. Save + summary ────────────────────────────────────────────
if results:
    pd.DataFrame(results).to_csv(OUT_CSV, index=False)
    log(f"\n  Results saved: {OUT_CSV}")

log("\n" + "="*60)
log("PATCH SCORECARD")
log("="*60)
for r in results:
    a = r.get("analysis","?")
    v = r.get("verdict","")
    p = r.get("logrank_p", r.get("p", np.nan))
    log(f"  {a:<25}  {v}  p={fmt_p(p)}")

flush()
log(f"\nLog saved: {OUT_LOG}")
log("PATCH COMPLETE.")
