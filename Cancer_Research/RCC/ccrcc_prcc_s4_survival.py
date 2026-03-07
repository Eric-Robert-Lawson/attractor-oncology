"""
RCC Script 4 — SURVIVAL ANALYSIS
ccRCC + PRCC — Attractor Depth vs Overall Survival

OrganismCore — Eric Robert Lawson
Document 94f-S4 | 2026-03-07

Run from: /Users/ericlawson/cancer/RCC/
  python ccrcc_prcc_s4_survival.py

═══════════════════════════════════════════════════════════════════
PREDICTIONS LOCKED — 2026-03-07 — BEFORE RUN

S4-P1 (PRIMARY — GOT1/RUNX1 TI):
  Transition Index predicts OS in ccRCC.
  TI_low (more negative) = deeper attractor = worse OS.
  Expected: log-rank p < 0.05, Cox HR > 1.5 (continuous −TI).
  Basis: TI r = −0.600 (n=534). Not in prior literature.

S4-P2 (ccRCC depth quartiles):
  Q4 OS < Q1 OS — log-rank p < 0.05.
  Q3+Q4 OS < Q1+Q2 OS — log-rank p < 0.05.

S4-P3 (depth score refinement):
  S5 depth score (SLC13A2/SLC2A1 anchor) predicts OS with
  higher C-index than S2 depth score.

S4-P4 (individual gene OS):
  High LOXL2 = worse OS  (Wall 3, r=+0.628)
  High TGFBI = worse OS  (Wall 3, r=+0.766)
  High RUNX1 = worse OS  (Wall 3, r=+0.559)
  Low  OGDHL = worse OS  (falls with depth)
  Low  SLC13A2 = worse OS (normal PT identity)
  All p < 0.05.

S4-P5 (PRCC FA-1 TI):
  TI (KRT19/SLC22A6) predicts OS in PRCC pooled.
  Type 2 OS < Type 1 OS (already p=0.034 in S5).
  FA-2 TI (LAMC2/SLC7A9) predicts OS in Type 2.

S4-P6 (panel vs single gene):
  panel_depth C-index ≥ LOXL2 alone in ccRCC.

═══════════════════════════════════════════════════════════════════
GEOMETRY-FIRST (Protocol v2.0):
  Predictions locked above before script written.
  Results graded against predictions after run.
  No retroactive adjustment of predictions.
  Every unexpected finding recorded explicitly.
═══════════════════════════════════════════════════════════════════
"""

import os
import sys
import gzip
import warnings
import traceback
import urllib.request
from datetime import datetime

import numpy  as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

# ── lifelines ─────────────────────────────────────────────────────
try:
    from lifelines import KaplanMeierFitter, CoxPHFitter
    from lifelines.statistics import logrank_test
    from lifelines.utils import concordance_index
except ImportError:
    import subprocess
    subprocess.run([sys.executable, "-m", "pip", "install",
                    "lifelines", "--quiet"])
    from lifelines import KaplanMeierFitter, CoxPHFitter
    from lifelines.statistics import logrank_test
    from lifelines.utils import concordance_index

# ── matplotlib ────────────────────────────────────────────────────
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    PLOT_OK = True
except ImportError:
    PLOT_OK = False

# ══════════════════════════════════════════════════════════════════
# PATHS  (confirmed from diagnostic output)
# ══════════════════════════════════════════════════════════════════

RCC_DIR    = os.path.dirname(os.path.abspath(__file__))

# ── ccRCC ─────────────────────────────────────────────────────────
CCRCC_BASE = os.path.join(RCC_DIR, "CCRCC", "ccrcc_false_attractor")
CCRCC_S1   = os.path.join(CCRCC_BASE, "results_s1")
CCRCC_S2   = os.path.join(CCRCC_BASE, "results_s2")
CCRCC_S5   = os.path.join(CCRCC_BASE, "results_s5")
S4_DIR     = os.path.join(CCRCC_BASE, "results_s4")   # output dir
os.makedirs(S4_DIR, exist_ok=True)

# ccRCC pre-computed files — exact paths + exact column names
# confirmed by diagnostic:
#   survival_depth.csv  : Unnamed:0, depth, stratum, os_time, os_event, depth_q
#   transition_index.csv: Unnamed:0, transition_index
#   depth_s5.csv        : Unnamed:0, depth_s5
#   panel_depth.csv     : Unnamed:0, panel_depth
#   drug_map.csv        : gene, desc, Q1, Q2, Q3, Q4, ratio_Q4Q1
CCRCC_SURV_FILE  = os.path.join(CCRCC_S2, "survival_depth.csv")
CCRCC_TI_FILE    = os.path.join(CCRCC_S5, "transition_index.csv")
CCRCC_DEPTH_S5   = os.path.join(CCRCC_S5, "depth_s5.csv")
CCRCC_PANEL_FILE = os.path.join(CCRCC_S5, "panel_depth.csv")
CCRCC_DRUG_MAP   = os.path.join(CCRCC_S5, "drug_map.csv")

# ccRCC raw expression + survival (cached on machine)
KIRC_EXPR_LOCAL  = os.path.join(CCRCC_BASE, "TCGA_KIRC_HiSeqV2.gz")
KIRC_SURV_LOCAL  = os.path.join(CCRCC_BASE, "TCGA_KIRC_survival.txt")

# ── PRCC ──────────────────────────────────────────────────────────
PRCC_BASE  = os.path.join(RCC_DIR, "PRCC", "prcc_false_attractor")
PRCC_S1    = os.path.join(PRCC_BASE, "results_s1")
PRCC_S2    = os.path.join(PRCC_BASE, "results_s2")
PRCC_S5    = os.path.join(PRCC_BASE, "results_s5")
PRCC_S6    = os.path.join(PRCC_BASE, "results_s6")

# PRCC pre-computed files — exact column names confirmed by diagnostic:
#   transition_index.csv : sample_id, s1_depth, TI
#   TI_FA2.csv           : sample_id, TI_FA2
#   integrated_gene_table: gene, r_depth, p_depth, r_TI, r_depth_T1,
#                          r_depth_T2, mean_normal, mean_tumour,
#                          T_minus_N, p_OS_pooled, med_OS_hi, med_OS_lo
#   drug_priority_map.csv: drug, gene, mean_Type1_all, mean_Type1_Q4,
#                          mean_Type2_all, mean_CIMP
PRCC_TI_FILE     = os.path.join(PRCC_S2, "transition_index.csv")
PRCC_TI_FA2_FILE = os.path.join(PRCC_S5, "TI_FA2.csv")
PRCC_INTEG_FILE  = os.path.join(PRCC_S6, "integrated_gene_table.csv")
PRCC_DRUG_FILE   = os.path.join(PRCC_S5, "drug_priority_map.csv")
PRCC_DEPTH_S1    = os.path.join(PRCC_S1, "depth_scores_tcga.csv")

# PRCC raw data (all cached on machine)
KIRP_EXPR_LOCAL  = os.path.join(PRCC_BASE, "TCGA_KIRP_HiSeqV2.gz")
KIRP_SURV_LOCAL  = os.path.join(PRCC_BASE, "KIRP_survival.txt")
KIRP_CLIN_LOCAL  = os.path.join(PRCC_BASE, "KIRP_clinicalMatrix.tsv")

# Output files
LOG_FILE      = os.path.join(S4_DIR, "s4_log.txt")
FIG_FILE      = os.path.join(S4_DIR, "s4_figure.png")
SURV_CSV      = os.path.join(S4_DIR, "survival_s4.csv")
COX_CSV       = os.path.join(S4_DIR, "cox_s4.csv")
GENE_OS_CSV   = os.path.join(S4_DIR, "gene_os_s4.csv")
PRCC_SURV_CSV = os.path.join(S4_DIR, "prcc_survival_s4.csv")

# Xena URLs — fallback only (files already cached)
KIRC_EXPR_URLS = [
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
    "TCGA.KIRC.sampleMap%2FHiSeqV2.gz",
    "https://tcga.xenahubs.net/download/TCGA.KIRC.sampleMap/HiSeqV2.gz",
]
KIRC_SURV_URLS = [
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
    "survival%2FKIRC_survival.txt",
    "https://tcga.xenahubs.net/download/survival/KIRC_survival.txt",
]
KIRP_SURV_URLS = [
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
    "survival%2FKIRP_survival.txt",
    "https://tcga.xenahubs.net/download/survival/KIRP_survival.txt",
]

# ══════════════════════════════════════════════════════════════════
# GENE PANELS
# ══════════════════════════════════════════════════════════════════

# Gene OS panel — ccRCC (OBJ-6)
GENE_OS_PANEL = [
    # Wall 1
    "EPAS1", "SLC2A1",
    # Wall 2
    "EZH2", "DNMT3A", "KDM1A", "HDAC1",
    # Wall 3 — primary S4-P4 predictions
    "LOXL2", "LOX", "RUNX1", "TGFBI", "CBFB",
    # Normal identity (low = worse OS predicted)
    "SLC13A2", "OGDHL", "GOT1", "FBP1", "SLC22A8",
    # Immune
    "CD274", "IL2RA", "FOXP3", "AXL", "IFI16",
]

# Expected directions (geometry-derived, locked before run)
EXPECTED_BAD  = {"EZH2", "LOXL2", "LOX", "RUNX1", "TGFBI",
                 "CBFB", "IL2RA", "FOXP3", "IFI16", "AXL",
                 "DNMT3A", "KDM1A", "HDAC1", "SLC2A1"}
EXPECTED_GOOD = {"SLC13A2", "OGDHL", "GOT1", "FBP1",
                 "SLC22A8", "CD274"}

# PRCC gene OS panel
PRCC_GENE_OS = [
    "KRT19", "ERBB2", "KRT7", "EZH2", "RUNX1",
    "MET",   "CDK4",  "MKI67",
    "SLC22A6", "FABP1", "FH", "OGDHL",
    "LAMC2", "KITLG", "SLC7A9",
    "CD274", "ARG1",  "B2M",
]
PRCC_EXPECTED_BAD  = {"KRT19", "ERBB2", "KRT7", "EZH2", "RUNX1",
                      "MET", "CDK4", "MKI67", "LAMC2", "KITLG"}
PRCC_EXPECTED_GOOD = {"SLC22A6", "FABP1", "FH", "OGDHL", "SLC7A9"}

# ══════════════════════════════════════════════════════════════════
# LOGGING
# ══════════════════════════════════════════════════════════════════

_log = []

def log(msg=""):
    _log.append(msg)
    print(msg)

def flush_log():
    with open(LOG_FILE, "w") as f:
        f.write("\n".join(_log))

# ══════════════════════════════════════════════════════════════════
# UTILITIES
# ══════════════════════════════════════════════════════════════════

def fmt(x, d=4):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "nan"
    return f"{x:.{d}f}"

def fmt_p(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "nan"
    if p < 1e-100: return "<1e-100"
    if p < 0.0001: return f"{p:.2e}"
    return f"{p:.4f}"

def minmax(s):
    mn, mx = s.min(), s.max()
    return (s - mn) / (mx - mn) if mx > mn else s * 0.0

def tumour_barcodes(columns):
    """Return only tumour sample barcodes (type codes 01–09)."""
    out = []
    for s in columns:
        parts = s.split("-")
        if len(parts) >= 4:
            code = parts[3][:2]
            if code.isdigit() and 1 <= int(code) <= 9:
                out.append(s)
    return out

def barcode12(s):
    """TCGA-AB-1234 — 12-char patient prefix."""
    return "-".join(s.split("-")[:3])

def download_file(url, dest, timeout=180):
    try:
        log(f"  Trying: {url}")
        req = urllib.request.Request(
            url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=timeout) as r, \
             open(dest, "wb") as fh:
            total = int(r.headers.get("Content-Length", 0))
            done  = 0
            while True:
                block = r.read(1024 * 1024)
                if not block:
                    break
                fh.write(block)
                done += len(block)
                if total:
                    print(f"\r  {min(done*100//total,100)}%",
                          end="", flush=True)
        print()
        log(f"  OK: {os.path.getsize(dest)/1e6:.1f} MB → {dest}")
        return True
    except Exception as e:
        log(f"  FAILED: {e}")
        if os.path.exists(dest):
            os.remove(dest)
        return False

def try_urls(urls, dest, label=""):
    if os.path.exists(dest) and os.path.getsize(dest) > 10_000:
        log(f"  Cached: {dest}")
        return True
    log(f"  Downloading {label}...")
    for url in urls:
        if download_file(url, dest):
            return True
    return False

def median_os_str(kmf):
    m = kmf.median_survival_time_
    if m is None or np.isnan(float(m)):
        return "NR"
    return f"{float(m):.0f}d"

# ══════════════════════════════════════════════════════════════════
# KM + COX WRAPPERS
# ══════════════════════════════════���═══════════════════════════════

def run_logrank(t_a, e_a, t_b, e_b):
    """Log-rank test. Returns p-value or nan."""
    if len(t_a) < 5 or len(t_b) < 5:
        return np.nan
    try:
        lr = logrank_test(t_a, t_b,
                          event_observed_A=e_a,
                          event_observed_B=e_b)
        return lr.p_value
    except Exception:
        return np.nan

def run_km_pair(df, dur, evt, grp_col, label_a, label_b):
    """Fit KM for two groups. Returns (kmf_a, kmf_b, logrank_p)."""
    a = df[df[grp_col] == label_a]
    b = df[df[grp_col] == label_b]
    if len(a) < 5 or len(b) < 5:
        return None, None, np.nan
    kmf_a = KaplanMeierFitter()
    kmf_b = KaplanMeierFitter()
    kmf_a.fit(a[dur], a[evt], label=label_a)
    kmf_b.fit(b[dur], b[evt], label=label_b)
    p = run_logrank(a[dur].values, a[evt].values,
                    b[dur].values, b[evt].values)
    return kmf_a, kmf_b, p

def run_cox(df, dur, evt, covariate):
    """
    Univariate Cox — covariate as continuous predictor.
    Returns (hr, ci_lo, ci_hi, p, c_index).
    """
    sub = df[[dur, evt, covariate]].dropna()
    if len(sub) < 10 or sub[evt].sum() < 3:
        return (np.nan,) * 5
    try:
        cph = CoxPHFitter()
        cph.fit(sub, duration_col=dur, event_col=evt,
                show_progress=False)
        row  = cph.summary.loc[covariate]
        hr   = float(np.exp(row["coef"]))
        cil  = float(np.exp(row["coef lower 95%"]))
        ciu  = float(np.exp(row["coef upper 95%"]))
        p    = float(row["p"])
        c    = concordance_index(
            sub[dur], -sub[covariate], sub[evt])
        return hr, cil, ciu, p, c
    except Exception as e:
        log(f"  Cox error ({covariate}): {e}")
        return (np.nan,) * 5

def run_cox_multi(df, dur, evt, covariates):
    """Multivariate Cox. Returns lifelines summary or None."""
    sub = df[[dur, evt] + covariates].dropna()
    if len(sub) < 20 or sub[evt].sum() < 5:
        return None
    try:
        cph = CoxPHFitter()
        cph.fit(sub, duration_col=dur, event_col=evt,
                show_progress=False)
        return cph.summary
    except Exception as e:
        log(f"  Multivariate Cox error: {e}")
        return None

# ══════════════════════════════════════════════════════════════════
# SECTION 1 — LOAD ccRCC PRE-COMPUTED DATA
# ══════════════════════════════════════════════════════════════════

def load_ccrcc():
    log("")
    log("=" * 65)
    log("SECTION 1 — LOAD ccRCC PRE-COMPUTED DATA")
    log("=" * 65)

    # ── survival_depth.csv ──────────────────���─────────────────────
    # Confirmed cols: Unnamed:0, depth, stratum, os_time, os_event,
    #                 depth_q   (n=532)
    surv = pd.read_csv(CCRCC_SURV_FILE, index_col=0)
    log(f"  survival_depth loaded: n={len(surv)}")
    log(f"  cols: {list(surv.columns)}")

    for col in ["depth", "os_time", "os_event", "depth_q"]:
        if col not in surv.columns:
            log(f"  FATAL: missing column '{col}'")
            sys.exit(1)

    surv["os_time"]  = pd.to_numeric(surv["os_time"],  errors="coerce")
    surv["os_event"] = pd.to_numeric(surv["os_event"], errors="coerce")
    surv["depth"]    = pd.to_numeric(surv["depth"],    errors="coerce")
    surv = surv.dropna(subset=["os_time", "os_event", "depth"])
    log(f"  Valid after dropna: n={len(surv)}  "
        f"events={int(surv['os_event'].sum())}")
    log(f"  OS range: {surv['os_time'].min():.0f}–"
        f"{surv['os_time'].max():.0f} days")
    log(f"  depth_q counts: {dict(surv['depth_q'].value_counts())}")

    # Binary groupings
    surv["depth_binary"] = surv["depth_q"].map(
        {"Q1": "Q1-Q2", "Q2": "Q1-Q2",
         "Q3": "Q3-Q4", "Q4": "Q3-Q4"})
    surv["depth_extreme"] = surv["depth_q"].apply(
        lambda q: q if q in ("Q1", "Q4") else np.nan)

    # ── transition_index.csv ──────────────────────────────────────
    # Confirmed cols: Unnamed:0, transition_index  (n=534)
    ti_raw = pd.read_csv(CCRCC_TI_FILE, index_col=0)
    log(f"  transition_index loaded: n={len(ti_raw)}")
    ti_raw.columns = ["ti"]
    log(f"  TI range: {ti_raw['ti'].min():.4f}–"
        f"{ti_raw['ti'].max():.4f}")

    # ── depth_s5.csv ──────────────────────────────────────────────
    # Confirmed cols: Unnamed:0, depth_s5  (n=534)
    d5_raw = pd.read_csv(CCRCC_DEPTH_S5, index_col=0)
    d5_raw.columns = ["depth_s5"]
    log(f"  depth_s5 loaded: n={len(d5_raw)}")

    # ── panel_depth.csv ───────────────────────────────────────────
    # Confirmed cols: Unnamed:0, panel_depth  (n=534)
    pd_raw = pd.read_csv(CCRCC_PANEL_FILE, index_col=0)
    pd_raw.columns = ["panel_depth"]
    log(f"  panel_depth loaded: n={len(pd_raw)}")

    # ── Merge ─────────────────────────────────────────────────────
    # All four frames share the same sample-ID index (Xena barcodes)
    # survival_depth n=532, the three S5 files n=534 → inner join
    master = surv.join(ti_raw,    how="left")
    master = master.join(d5_raw,  how="left")
    master = master.join(pd_raw,  how="left")
    log(f"  Master frame: n={len(master)}  "
        f"TI non-null={master['ti'].notna().sum()}")

    return master

# ══════════════════════════════════════════════════════════════════
# SECTION 2 — LOAD ccRCC EXPRESSION FOR GENE OS TESTS
# ═════════════════════════════════��════════════════════════════════

def load_ccrcc_expr():
    """
    Merges TCGA_KIRC_HiSeqV2.gz (cached) with KIRC_survival.txt
    (cached). Returns a sample × [gene... + os_time + os_event]
    DataFrame for OBJ-6 gene OS analysis.
    Joins on 12-char patient barcode.
    """
    log("")
    log("=" * 65)
    log("SECTION 2 — ccRCC EXPRESSION + SURVIVAL (GENE OS TESTS)")
    log("=" * 65)

    # Expression is already cached
    log("  Loading KIRC expression matrix...")
    with gzip.open(KIRC_EXPR_LOCAL, "rt") as f:
        raw = pd.read_csv(f, sep="\t", index_col=0)
    avail = [g for g in GENE_OS_PANEL if g in raw.index]
    miss  = [g for g in GENE_OS_PANEL if g not in raw.index]
    if miss:
        log(f"  Genes not in matrix: {miss}")
    log(f"  {len(avail)} / {len(GENE_OS_PANEL)} genes found")

    tumour_cols = tumour_barcodes(raw.columns)
    expr = raw.loc[avail, tumour_cols].T
    log(f"  Tumour samples: {len(expr)}")

    # Survival file is cached
    ok = try_urls(KIRC_SURV_URLS, KIRC_SURV_LOCAL, "KIRC survival")
    if not ok and not os.path.exists(KIRC_SURV_LOCAL):
        log("  WARNING: KIRC survival not available. "
            "OBJ-6 will use drug_map fallback.")
        return None

    surv_raw = pd.read_csv(KIRC_SURV_LOCAL, sep="\t")
    log(f"  Survival cols: {list(surv_raw.columns)}")

    # Identify columns — Xena survival files use OS.time / OS
    time_col = next((c for c in surv_raw.columns
                     if "os.time" in c.lower()), None)
    evt_col  = next((c for c in surv_raw.columns
                     if c.lower() in ("os", "_os")), None)
    sid_col  = surv_raw.columns[0]   # always first column

    if time_col is None or evt_col is None:
        log(f"  WARNING: cannot identify OS columns in "
            f"{list(surv_raw.columns)}")
        return None

    surv = surv_raw[[sid_col, time_col, evt_col]].copy()
    surv.columns = ["sample", "os_time", "os_event"]
    surv["os_time"]  = pd.to_numeric(surv["os_time"],  errors="coerce")
    surv["os_event"] = pd.to_numeric(surv["os_event"], errors="coerce")
    surv = surv.dropna()
    log(f"  Survival valid: n={len(surv)}  "
        f"events={int(surv['os_event'].sum())}")

    # Join on 12-char barcode
    expr.index    = [barcode12(s) for s in expr.index]
    surv["b12"]   = [barcode12(s) for s in surv["sample"]]
    surv          = surv.set_index("b12")
    merged        = expr.join(surv[["os_time", "os_event"]],
                               how="inner")
    merged        = merged.dropna(subset=["os_time", "os_event"])
    log(f"  Expression+survival merged: n={len(merged)}")
    return merged

# ══════════════════════════════════════════════════════════════════
# OBJ-1 — ccRCC KM: DEPTH QUARTILES
# ══════════════════════════════════════════════════════════════════

def obj1_km_quartiles(master):
    log("")
    log("=" * 65)
    log("OBJ-1 — ccRCC KM: DEPTH QUARTILES")
    log("  S4-P2: Q4 < Q1  and  Q3+Q4 < Q1+Q2  (p < 0.05)")
    log("=" * 65)

    results = []

    # ── A: Q4 vs Q1 ──────────────────────────────────────────────
    log("")
    log("  A: Q4 vs Q1 (extreme quartiles)")
    sub = master[master["depth_extreme"].isin(["Q1","Q4"])].dropna(
        subset=["depth_extreme", "os_time", "os_event"])

    kmf_q1, kmf_q4, p_q = run_km_pair(
        sub, "os_time", "os_event", "depth_extreme", "Q1", "Q4")

    if kmf_q1 is not None:
        nq1 = (sub["depth_extreme"] == "Q1").sum()
        nq4 = (sub["depth_extreme"] == "Q4").sum()
        log(f"  Q1 n={nq1}  median={median_os_str(kmf_q1)}")
        log(f"  Q4 n={nq4}  median={median_os_str(kmf_q4)}")
        log(f"  Log-rank p = {fmt_p(p_q)}")
        v = "CONFIRMED ✓" if (not np.isnan(p_q) and p_q < 0.05) \
            else "NOT CONFIRMED ✗"
        log(f"  S4-P2 (Q4<Q1): {v}")
        results.append({"analysis": "Q4_vs_Q1",
                         "n_q1": nq1, "n_q4": nq4,
                         "logrank_p": p_q, "verdict": v})

    # ── B: Q3+Q4 vs Q1+Q2 ────────────────────────────────────────
    log("")
    log("  B: Q3+Q4 vs Q1+Q2 (combined)")
    sub2 = master.dropna(
        subset=["depth_binary", "os_time", "os_event"])

    kmf_sh, kmf_dp, p_b = run_km_pair(
        sub2, "os_time", "os_event", "depth_binary",
        "Q1-Q2", "Q3-Q4")

    if kmf_sh is not None:
        nsh = (sub2["depth_binary"] == "Q1-Q2").sum()
        ndp = (sub2["depth_binary"] == "Q3-Q4").sum()
        log(f"  Shallow Q1+Q2 n={nsh}  median={median_os_str(kmf_sh)}")
        log(f"  Deep    Q3+Q4 n={ndp}  median={median_os_str(kmf_dp)}")
        log(f"  Log-rank p = {fmt_p(p_b)}")
        v2 = "CONFIRMED ✓" if (not np.isnan(p_b) and p_b < 0.05) \
             else "NOT CONFIRMED ✗"
        log(f"  S4-P2 (Q3+Q4 < Q1+Q2): {v2}")
        results.append({"analysis": "Q3Q4_vs_Q1Q2",
                         "n_shallow": nsh, "n_deep": ndp,
                         "logrank_p": p_b, "verdict": v2})

    return results, (kmf_q1, kmf_q4, p_q, sub,
                     kmf_sh, kmf_dp, p_b, sub2)

# ══════════════════════════════════════════════════════════════════
# OBJ-2 — ccRCC COX: CONTINUOUS DEPTH
# ══════════════════════════════════════════════════════════════════

def obj2_cox_continuous(master):
    log("")
    log("=" * 65)
    log("OBJ-2 — ccRCC COX: CONTINUOUS DEPTH PREDICTORS")
    log("  S4-P3: S5 depth C-index > S2 depth C-index")
    log("=" * 65)

    cox_rows = []

    for label, col in [("S2_depth",    "depth"),
                        ("S5_depth",    "depth_s5"),
                        ("panel_depth", "panel_depth")]:
        if col not in master.columns:
            log(f"  SKIP {label}: column not in master")
            continue
        sub = master.dropna(subset=[col, "os_time", "os_event"])
        hr, cil, ciu, p, c = run_cox(sub, "os_time", "os_event", col)
        log(f"  {label:<14} HR={fmt(hr)} [{fmt(cil)}–{fmt(ciu)}]  "
            f"p={fmt_p(p)}  C={fmt(c)}  n={len(sub)}")
        cox_rows.append({"predictor": label,
                          "hr": hr, "ci_lo": cil, "ci_hi": ciu,
                          "p": p, "c_index": c, "n": len(sub)})

    # S4-P3 verdict
    c_s2 = next((r["c_index"] for r in cox_rows
                 if r["predictor"] == "S2_depth"), np.nan)
    c_s5 = next((r["c_index"] for r in cox_rows
                 if r["predictor"] == "S5_depth"), np.nan)
    if not np.isnan(c_s2) and not np.isnan(c_s5):
        v = "CONFIRMED ✓" if c_s5 > c_s2 else "NOT CONFIRMED ✗"
        log(f"  S4-P3: {v}  (C_S5={fmt(c_s5)}  C_S2={fmt(c_s2)})")

    return cox_rows

# ══════════════════════════════════════════════════════════════════
# OBJ-3 — GOT1/RUNX1 TRANSITION INDEX (PRIMARY NOVEL CLAIM)
# ══════════════════════════════════════════════════════════════════

def obj3_ti_survival(master):
    """
    THE PRIMARY NOVEL CLAIM.
    TI = norm(GOT1) - norm(RUNX1) = proxy for normal metabolic
    identity vs chromatin lock force.
    High TI = more GOT1-like = shallower = better OS (expected).
    Model −TI: higher = deeper = worse OS. Expected HR > 1.5.
    """
    log("")
    log("=" * 65)
    log("OBJ-3 — GOT1/RUNX1 TRANSITION INDEX vs OS  [PRIMARY]")
    log("  S4-P1: −TI continuous Cox: HR > 1.5  p < 0.05")
    log("=" * 65)

    if "ti" not in master.columns:
        log("  SKIP: ti column not in master.")
        return [], None

    sub = master.dropna(subset=["ti", "os_time", "os_event"]).copy()
    log(f"  n={len(sub)}  events={int(sub['os_event'].sum())}")
    log(f"  TI range: {sub['ti'].min():.4f}–{sub['ti'].max():.4f}")
    log(f"  TI median: {sub['ti'].median():.4f}")

    # Continuous Cox: model −TI (higher = deeper = worse)
    sub["ti_neg"] = -sub["ti"]
    hr, cil, ciu, p, c = run_cox(sub, "os_time", "os_event", "ti_neg")
    sig    = (not np.isnan(p)) and p < 0.05
    effect = (not np.isnan(hr)) and hr > 1.5
    v = ("CONFIRMED ✓" if (sig and effect) else
         "PARTIAL ↯"   if sig              else "NOT CONFIRMED ✗")

    log(f"  Cox (−TI continuous):")
    log(f"    HR={fmt(hr)} [{fmt(cil)}–{fmt(ciu)}]  "
        f"p={fmt_p(p)}  C={fmt(c)}")
    log(f"  S4-P1 verdict: {v}")

    # KM median split
    ti_med = sub["ti"].median()
    sub["ti_grp"] = np.where(sub["ti"] >= ti_med,
                             "TI-High (shallow)",
                             "TI-Low (deep)")
    kmf_h, kmf_l, p_km = run_km_pair(
        sub, "os_time", "os_event", "ti_grp",
        "TI-High (shallow)", "TI-Low (deep)")

    if kmf_h is not None:
        nh = (sub["ti_grp"] == "TI-High (shallow)").sum()
        nl = (sub["ti_grp"] == "TI-Low (deep)").sum()
        log(f"  KM median split (TI ≥ {ti_med:.4f} = High):")
        log(f"    TI-High n={nh}  median={median_os_str(kmf_h)}")
        log(f"    TI-Low  n={nl}  median={median_os_str(kmf_l)}")
        log(f"    Log-rank p={fmt_p(p_km)}")

    ti_results = [{"predictor": "TI_neg_continuous",
                   "hr": hr, "ci_lo": cil, "ci_hi": ciu,
                   "p": p, "c_index": c, "verdict": v,
                   "n": len(sub),
                   "events": int(sub["os_event"].sum())}]
    if kmf_h is not None:
        ti_results.append({"predictor": "TI_KM_median_split",
                            "logrank_p": p_km,
                            "n_high": nh, "n_low": nl})

    return ti_results, (sub, kmf_h, kmf_l, p_km)

# ══════════════════════════════════════════════════════════════════
# OBJ-6 — ccRCC INDIVIDUAL GENE OS
# ══════════════════════════════════════════════════════════════════

def obj6_gene_os(merged):
    """
    For each gene: continuous Cox + direction check.
    Falls back to drug_map Q4/Q1 ratios if no expression data.
    """
    log("")
    log("=" * 65)
    log("OBJ-6 — ccRCC INDIVIDUAL GENE OS")
    log("  S4-P4: LOXL2/TGFBI/RUNX1 high = worse OS")
    log("        OGDHL/SLC13A2 low = worse OS")
    log("=" * 65)

    gene_rows = []

    if merged is None:
        log("  Expression+survival not available.")
        log("  Reporting drug_map Q4/Q1 ratios as directional proxy:")
        dm = pd.read_csv(CCRCC_DRUG_MAP)
        log(f"  {'Gene':<12} {'Q4/Q1':>8}  Direction  Interpretation")
        log("  " + "-" * 60)
        for _, row in dm.iterrows():
            g     = row.get("gene", "?")
            ratio = row.get("ratio_Q4Q1", np.nan)
            desc  = row.get("desc", "")
            if np.isnan(ratio):
                continue
            d = ("UP in Q4 (deep)"   if ratio > 1.05 else
                 "DOWN in Q4 (deep)" if ratio < 0.95 else "FLAT")
            log(f"  {g:<12} {ratio:>8.4f}  {d:<22} {desc}")
        return gene_rows

    log(f"  Testing {len(GENE_OS_PANEL)} genes  n={len(merged)}")
    log(f"  {'Gene':<12} {'HR':>7} {'95%CI_lo':>9} {'95%CI_hi':>9} "
        f"{'p':>10}  {'C':>6}  Verdict")
    log("  " + "-" * 72)

    for gene in GENE_OS_PANEL:
        if gene not in merged.columns:
            log(f"  {gene:<12} NOT IN MATRIX")
            continue
        sub = merged[[gene, "os_time", "os_event"]].dropna()
        if len(sub) < 20 or sub["os_event"].sum() < 3:
            log(f"  {gene:<12} INSUFFICIENT DATA (n={len(sub)})")
            continue

        hr, cil, ciu, p, c = run_cox(sub, "os_time", "os_event", gene)
        sig = (not np.isnan(p)) and p < 0.05

        if gene in EXPECTED_BAD:
            correct = (not np.isnan(hr)) and hr > 1.0
        elif gene in EXPECTED_GOOD:
            correct = (not np.isnan(hr)) and hr < 1.0
        else:
            correct = None

        if sig and correct is True:
            v = "CONFIRMED ✓"
        elif sig and correct is False:
            v = "DIRECTION ✗"
        elif sig:
            v = "SIG"
        else:
            v = "NS"

        log(f"  {gene:<12} {fmt(hr):>7} {fmt(cil):>9} {fmt(ciu):>9} "
            f"{fmt_p(p):>10}  {fmt(c):>6}  {v}")
        gene_rows.append({"gene": gene, "hr": hr,
                           "ci_lo": cil, "ci_hi": ciu,
                           "p": p, "c_index": c,
                           "verdict": v, "n": len(sub),
                           "events": int(sub["os_event"].sum())})
    return gene_rows

# ══════════════════════════════════════════════════════════════════
# OBJ-7 — ccRCC MULTIVARIATE COX
# ══════════════════════════════════════════════════════════════════

def obj7_multivariate(master, merged):
    log("")
    log("=" * 65)
    log("OBJ-7 — ccRCC MULTIVARIATE COX")
    log("  Model A: −TI + S2_depth")
    log("  Model B: −TI + LOXL2 + TGFBI  (if expression available)")
    log("=" * 65)

    multi_rows = []

    # ── Model A: −TI + S2_depth ───────────────────────────────────
    if "ti" in master.columns:
        sub = master[["os_time","os_event","ti","depth"]].dropna().copy()
        sub["ti_neg"] = -sub["ti"]
        summ = run_cox_multi(sub, "os_time", "os_event",
                             ["ti_neg", "depth"])
        if summ is not None:
            log(f"  Model A — −TI + S2 depth  (n={len(sub)}):")
            for cov in ["ti_neg", "depth"]:
                if cov in summ.index:
                    hr2 = float(np.exp(summ.loc[cov, "coef"]))
                    p2  = float(summ.loc[cov, "p"])
                    log(f"    {cov:<12}  HR={fmt(hr2)}  p={fmt_p(p2)}")
            multi_rows.append({"model": "TI_neg+depth_s2",
                                 "n": len(sub)})

    # ── Model B: −TI + LOXL2 + TGFBI ─────────────────────────────
    if merged is not None and "ti" in master.columns:
        extra = [g for g in ["LOXL2", "TGFBI", "RUNX1"]
                 if g in merged.columns]
        if extra:
            # align TI onto merged by 12-char barcode
            ti_map = {}
            for idx, row in master.iterrows():
                k = barcode12(str(idx))
                if not np.isnan(row.get("ti", np.nan)):
                    ti_map[k] = row["ti"]
            merged2 = merged.copy()
            merged2["b12"] = [barcode12(str(i))
                              for i in merged2.index]
            merged2["ti"]     = merged2["b12"].map(ti_map)
            merged2["ti_neg"] = -merged2["ti"]
            cov2 = ["ti_neg"] + extra
            sub2 = merged2[["os_time","os_event"] + cov2].dropna()
            if len(sub2) >= 20:
                summ2 = run_cox_multi(
                    sub2, "os_time", "os_event", cov2)
                if summ2 is not None:
                    log(f"  Model B — −TI + {extra}  (n={len(sub2)}):")
                    for cov in cov2:
                        if cov in summ2.index:
                            hr3 = float(
                                np.exp(summ2.loc[cov, "coef"]))
                            p3  = float(summ2.loc[cov, "p"])
                            log(f"    {cov:<12}  HR={fmt(hr3)}  "
                                f"p={fmt_p(p3)}")
                    multi_rows.append(
                        {"model": f"TI_neg+{'_'.join(extra)}",
                         "n": len(sub2)})

    return multi_rows

# ══════════════════════════════════════════���═══════════════════════
# OBJ-8 — PRCC SURVIVAL
# ══════════════════════════════════════════════════════════════════

def obj8_prcc():
    """
    PRCC has no pre-built survival_depth.csv.
    Build from:
      TCGA_KIRP_HiSeqV2.gz  (cached, 17.4 MB)
      KIRP_survival.txt     (cached, 18.8 KB)
      KIRP_clinicalMatrix.tsv (cached, 338 KB)
      transition_index.csv  (s2: sample_id, s1_depth, TI  n=290)
      TI_FA2.csv            (s5: sample_id, TI_FA2  n=290)
      integrated_gene_table (s6: gene, p_OS_pooled, med_OS_hi,
                                 med_OS_lo — already has OS info)
    Prediction S4-P5:
      TI (FA-1) predicts OS pooled + Type 1 specifically.
      TI_FA2 predicts OS in Type 2 specifically.
    """
    log("")
    log("=" * 65)
    log("OBJ-8 — PRCC SURVIVAL ANALYSIS")
    log("  S4-P5: FA-1 TI predicts OS  (pooled + Type 1)")
    log("        FA-2 TI predicts OS  (Type 2)")
    log("=" * 65)

    prcc_rows = []

    # ── Load KIRP expression (cached) ────────────────────────────
    log("  Loading KIRP expression...")
    with gzip.open(KIRP_EXPR_LOCAL, "rt") as f:
        raw = pd.read_csv(f, sep="\t", index_col=0)

    prcc_genes = list(set(PRCC_GENE_OS))
    avail  = [g for g in prcc_genes if g in raw.index]
    miss   = [g for g in prcc_genes if g not in raw.index]
    if miss:
        log(f"  Genes not in matrix: {miss}")
    tumour = tumour_barcodes(raw.columns)
    expr   = raw.loc[avail, tumour].T
    log(f"  Tumour samples: {len(expr)}")

    # ── Load KIRP survival (cached) ───────────────────────────────
    log("  Loading KIRP survival...")
    ok = try_urls(KIRP_SURV_URLS, KIRP_SURV_LOCAL, "KIRP survival")
    if not ok and not os.path.exists(KIRP_SURV_LOCAL):
        log("  WARNING: KIRP survival not available. Skipping OBJ-8.")
        return prcc_rows

    surv_raw = pd.read_csv(KIRP_SURV_LOCAL, sep="\t")
    log(f"  Survival cols: {list(surv_raw.columns)}")

    time_col = next((c for c in surv_raw.columns
                     if "os.time" in c.lower()), None)
    evt_col  = next((c for c in surv_raw.columns
                     if c.lower() in ("os", "_os")), None)
    sid_col  = surv_raw.columns[0]

    if time_col is None or evt_col is None:
        log(f"  WARNING: Cannot identify OS columns. "
            f"Cols: {list(surv_raw.columns)}")
        return prcc_rows

    surv = surv_raw[[sid_col, time_col, evt_col]].copy()
    surv.columns = ["sample", "os_time", "os_event"]
    surv["os_time"]  = pd.to_numeric(surv["os_time"],  errors="coerce")
    surv["os_event"] = pd.to_numeric(surv["os_event"], errors="coerce")
    surv = surv.dropna()
    log(f"  Survival valid: n={len(surv)}  "
        f"events={int(surv['os_event'].sum())}")

    # ── Load clinical (tumor_type: Type 1 / Type 2) ───────────────
    log("  Loading KIRP clinical...")
    clin = pd.read_csv(KIRP_CLIN_LOCAL, sep="\t", low_memory=False)
    log(f"  Clinical cols (first 20): {list(clin.columns[:20])}")

    # Find tumor_type column
    tt_col = next((c for c in clin.columns
                   if "tumor_type" in c.lower()
                   or "histological" in c.lower()
                   or "type" in c.lower()), None)
    sid_c  = clin.columns[0]
    if tt_col:
        log(f"  Tumor-type column: {tt_col}")
        log(f"  Values: {dict(clin[tt_col].value_counts())}")
        clin_sub = clin[[sid_c, tt_col]].copy()
        clin_sub.columns = ["sample", "tumor_type"]
    else:
        log("  No tumor_type column found — pooled analysis only.")
        clin_sub = None

    # ── Merge expression + survival ───────────────────────────────
    expr.index  = [barcode12(s) for s in expr.index]
    surv["b12"] = [barcode12(s) for s in surv["sample"]]
    surv        = surv.set_index("b12")
    merged      = expr.join(surv[["os_time","os_event"]], how="inner")
    merged      = merged.dropna(subset=["os_time","os_event"])
    log(f"  Expression+survival merged: n={len(merged)}")

    # ── Load TI files ─────────────────────────────────────────────
    # transition_index.csv: sample_id, s1_depth, TI  (n=290)
    ti_fa1 = pd.read_csv(PRCC_TI_FILE)
    ti_fa1 = ti_fa1[["sample_id", "TI"]].set_index("sample_id")
    ti_fa1.index = [barcode12(str(s)) for s in ti_fa1.index]
    ti_fa1.columns = ["TI_FA1"]
    log(f"  FA-1 TI loaded: n={len(ti_fa1)}")

    # TI_FA2.csv: sample_id, TI_FA2  (n=290)
    ti_fa2 = pd.read_csv(PRCC_TI_FA2_FILE)
    ti_fa2 = ti_fa2[["sample_id","TI_FA2"]].set_index("sample_id")
    ti_fa2.index = [barcode12(str(s)) for s in ti_fa2.index]
    log(f"  FA-2 TI loaded: n={len(ti_fa2)}")

    merged = merged.join(ti_fa1, how="left")
    merged = merged.join(ti_fa2, how="left")

    # Add tumor type
    if clin_sub is not None:
        clin_sub["b12"] = [barcode12(str(s))
                           for s in clin_sub["sample"]]
        clin_sub = clin_sub.set_index("b12")
        merged = merged.join(clin_sub[["tumor_type"]], how="left")
        log(f"  After tumor_type join: "
            f"{dict(merged['tumor_type'].value_counts())}")

    log(f"  Final merged: n={len(merged)}  "
        f"TI_FA1 non-null={merged['TI_FA1'].notna().sum()}  "
        f"TI_FA2 non-null={merged['TI_FA2'].notna().sum()}")

    # ── Analysis A: FA-1 TI vs OS (pooled) ───────────────────────
    log("")
    log("  ANALYSIS A: FA-1 TI vs OS (pooled)")
    sub_a = merged.dropna(
        subset=["TI_FA1","os_time","os_event"]).copy()
    # Higher TI_FA1 = deeper FA-1 = worse OS (KRT19 high = deep)
    hr_a, cil_a, ciu_a, p_a, c_a = run_cox(
        sub_a, "os_time", "os_event", "TI_FA1")
    log(f"  Cox (TI_FA1 continuous): HR={fmt(hr_a)}  "
        f"p={fmt_p(p_a)}  C={fmt(c_a)}  n={len(sub_a)}")

    # KM median split
    ti1_med = sub_a["TI_FA1"].median()
    sub_a["ti1_grp"] = np.where(sub_a["TI_FA1"] >= ti1_med,
                                "FA1-High (deep)",
                                "FA1-Low (shallow)")
    kmf_h1, kmf_l1, p_km1 = run_km_pair(
        sub_a, "os_time", "os_event", "ti1_grp",
        "FA1-High (deep)", "FA1-Low (shallow)")
    if kmf_h1 is not None:
        log(f"  KM split: FA1-High n="
            f"{(sub_a['ti1_grp']=='FA1-High (deep)').sum()}  "
            f"median={median_os_str(kmf_h1)}")
        log(f"           FA1-Low  n="
            f"{(sub_a['ti1_grp']=='FA1-Low (shallow)').sum()}  "
            f"median={median_os_str(kmf_l1)}")
        log(f"           log-rank p={fmt_p(p_km1)}")
        v_a = "CONFIRMED ✓" if (not np.isnan(p_km1) and p_km1 < 0.05) \
              else "NOT CONFIRMED ✗"
        log(f"  S4-P5 (FA-1 TI pooled): {v_a}")
        prcc_rows.append({"analysis": "FA1_TI_pooled",
                           "hr": hr_a, "p": p_a, "c": c_a,
                           "logrank_p": p_km1, "verdict": v_a})

    # ── Analysis B: FA-2 TI vs OS (pooled) ───────────────────────
    log("")
    log("  ANALYSIS B: FA-2 TI vs OS (pooled)")
    sub_b = merged.dropna(
        subset=["TI_FA2","os_time","os_event"]).copy()
    hr_b, cil_b, ciu_b, p_b, c_b = run_cox(
        sub_b, "os_time", "os_event", "TI_FA2")
    log(f"  Cox (TI_FA2 continuous): HR={fmt(hr_b)}  "
        f"p={fmt_p(p_b)}  C={fmt(c_b)}  n={len(sub_b)}")

    ti2_med = sub_b["TI_FA2"].median()
    sub_b["ti2_grp"] = np.where(sub_b["TI_FA2"] >= ti2_med,
                                "FA2-High (deep)",
                                "FA2-Low (shallow)")
    kmf_h2, kmf_l2, p_km2 = run_km_pair(
        sub_b, "os_time", "os_event", "ti2_grp",
        "FA2-High (deep)", "FA2-Low (shallow)")
    if kmf_h2 is not None:
        log(f"  KM split log-rank p={fmt_p(p_km2)}")
        prcc_rows.append({"analysis": "FA2_TI_pooled",
                           "hr": hr_b, "p": p_b, "c": c_b,
                           "logrank_p": p_km2})

    # ── Analysis C: Type 1 vs Type 2 OS ──────────────────────────
    if "tumor_type" in merged.columns:
        log("")
        log("  ANALYSIS C: Type 1 vs Type 2 OS")
        sub_c = merged.dropna(
            subset=["tumor_type","os_time","os_event"])
        # normalise type labels
        sub_c = sub_c.copy()
        sub_c["tt_clean"] = sub_c["tumor_type"].astype(
            str).str.strip().str.lower()
        t1_mask = sub_c["tt_clean"].str.contains("1|type1|papillary 1",
                                                   regex=True)
        t2_mask = sub_c["tt_clean"].str.contains("2|type2|papillary 2",
                                                   regex=True)
        t1 = sub_c[t1_mask]
        t2 = sub_c[t2_mask]
        log(f"  Type 1 n={len(t1)}  Type 2 n={len(t2)}")
        if len(t1) >= 5 and len(t2) >= 5:
            p_tt = run_logrank(t1["os_time"].values,
                               t1["os_event"].values,
                               t2["os_time"].values,
                               t2["os_event"].values)
            log(f"  Type 1 vs Type 2 log-rank p={fmt_p(p_tt)}")
            v_c = "CONFIRMED ✓" if (not np.isnan(p_tt) and p_tt < 0.05) \
                  else "NOT CONFIRMED ✗"
            log(f"  S4-P5 (Type2 < Type1 OS): {v_c}")
            prcc_rows.append({"analysis": "Type1_vs_Type2",
                               "logrank_p": p_tt, "verdict": v_c,
                               "n_t1": len(t1), "n_t2": len(t2)})

    # ── Analysis D: PRCC individual gene OS ───────────────────────
    log("")
    log("  ANALYSIS D: PRCC individual gene OS")
    log(f"  {'Gene':<12} {'HR':>7} {'p':>10}  {'C':>6}  Verdict")
    log("  " + "-" * 50)

    for gene in PRCC_GENE_OS:
        if gene not in merged.columns:
            continue
        sub_g = merged[[gene,"os_time","os_event"]].dropna()
        if len(sub_g) < 10 or sub_g["os_event"].sum() < 3:
            continue
        hr_g, cil_g, ciu_g, p_g, c_g = run_cox(
            sub_g, "os_time", "os_event", gene)
        sig_g = (not np.isnan(p_g)) and p_g < 0.05
        if gene in PRCC_EXPECTED_BAD:
            correct_g = (not np.isnan(hr_g)) and hr_g > 1.0
        elif gene in PRCC_EXPECTED_GOOD:
            correct_g = (not np.isnan(hr_g)) and hr_g < 1.0
        else:
            correct_g = None
        if sig_g and correct_g is True:
            vg = "CONFIRMED ✓"
        elif sig_g and correct_g is False:
            vg = "DIRECTION ✗"
        elif sig_g:
            vg = "SIG"
        else:
            vg = "NS"
        log(f"  {gene:<12} {fmt(hr_g):>7} {fmt_p(p_g):>10}  "
            f"{fmt(c_g):>6}  {vg}")
        prcc_rows.append({"analysis": f"gene_{gene}",
                           "hr": hr_g, "p": p_g,
                           "c_index": c_g, "verdict": vg,
                           "n": len(sub_g)})

    # ── Also surface integrated_gene_table OS column ──────────────
    log("")
    log("  ANALYSIS E: integrated_gene_table — OS results from S6")
    integ = pd.read_csv(PRCC_INTEG_FILE)
    os_genes = integ[integ["p_OS_pooled"].notna()].sort_values(
        "p_OS_pooled")
    log(f"  Genes with OS data in integrated table: "
        f"{len(os_genes)}")
    log(f"  Top 20 by p_OS_pooled:")
    log(f"  {'Gene':<12} {'p_OS':>10}  {'med_hi':>10}  "
        f"{'med_lo':>10}")
    log("  " + "-" * 50)
    for _, row in os_genes.head(20).iterrows():
        log(f"  {row['gene']:<12} "
            f"{fmt_p(row['p_OS_pooled']):>10}  "
            f"{str(row.get('med_OS_hi',''))[:10]:>10}  "
            f"{str(row.get('med_OS_lo',''))[:10]:>10}")

    # Save
    merged.to_csv(PRCC_SURV_CSV)
    log(f"  PRCC survival frame saved: {PRCC_SURV_CSV}")
    return prcc_rows

# ══════════════════════════════════════════════════════════════════
# OBJ-9 — FIGURE
# ══════════════════════════════════════════════════════════════════

def make_figure(master, ti_data, prcc_merged_path):
    if not PLOT_OK:
        log("  SKIP FIGURE: matplotlib not available.")
        return

    log("")
    log("=" * 65)
    log("OBJ-9 — FIGURE: 4-PANEL")
    log("=" * 65)

    fig = plt.figure(figsize=(14, 11))
    gs  = gridspec.GridSpec(2, 2, hspace=0.38, wspace=0.32)
    ax1 = fig.add_subplot(gs[0, 0])   # A: Q4 vs Q1
    ax2 = fig.add_subplot(gs[0, 1])   # B: Q1+Q2 vs Q3+Q4
    ax3 = fig.add_subplot(gs[1, 0])   # C: TI KM
    ax4 = fig.add_subplot(gs[1, 1])   # D: depth scatter

    PAL = {"Q1": "#27ae60", "Q4": "#e74c3c",
           "Q1-Q2": "#2ecc71", "Q3-Q4": "#e74c3c",
           "TI-High (shallow)": "#2196F3",
           "TI-Low (deep)":     "#e74c3c"}

    def style(ax, title, xl="Time (days)", yl="Survival prob."):
        ax.set_title(title, fontsize=9, fontweight="bold", pad=5)
        ax.set_xlabel(xl, fontsize=8)
        ax.set_ylabel(yl, fontsize=8)
        ax.set_ylim(-0.03, 1.07)
        ax.tick_params(labelsize=7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(fontsize=7, loc="upper right")

    # ── Panel A: Q4 vs Q1 ────────────────────────────────────────
    sub_ext = master[master["depth_extreme"].isin(["Q1","Q4"])].dropna(
        subset=["depth_extreme","os_time","os_event"])
    for q, col in [("Q1", PAL["Q1"]), ("Q4", PAL["Q4"])]:
        g = sub_ext[sub_ext["depth_extreme"] == q]
        if len(g) >= 3:
            KaplanMeierFitter().fit(
                g["os_time"], g["os_event"], label=q
            ).plot_survival_function(ax=ax1, ci_show=True,
                                     color=col, linewidth=1.5)
    p_q = run_logrank(
        sub_ext[sub_ext["depth_extreme"]=="Q1"]["os_time"].values,
        sub_ext[sub_ext["depth_extreme"]=="Q1"]["os_event"].values,
        sub_ext[sub_ext["depth_extreme"]=="Q4"]["os_time"].values,
        sub_ext[sub_ext["depth_extreme"]=="Q4"]["os_event"].values)
    ax1.text(0.05, 0.08, f"p = {fmt_p(p_q)}",
             transform=ax1.transAxes, fontsize=8)
    style(ax1, "A.  ccRCC — Depth Q4 vs Q1")

    # ── Panel B: Q3+Q4 vs Q1+Q2 ──────────────────────────────────
    sub2 = master.dropna(
        subset=["depth_binary","os_time","os_event"])
    for grp, col in [("Q1-Q2", PAL["Q1-Q2"]),
                     ("Q3-Q4", PAL["Q3-Q4"])]:
        g = sub2[sub2["depth_binary"] == grp]
        if len(g) >= 5:
            KaplanMeierFitter().fit(
                g["os_time"], g["os_event"], label=grp
            ).plot_survival_function(ax=ax2, ci_show=True,
                                     color=col, linewidth=1.5)
    p_b = run_logrank(
        sub2[sub2["depth_binary"]=="Q1-Q2"]["os_time"].values,
        sub2[sub2["depth_binary"]=="Q1-Q2"]["os_event"].values,
        sub2[sub2["depth_binary"]=="Q3-Q4"]["os_time"].values,
        sub2[sub2["depth_binary"]=="Q3-Q4"]["os_event"].values)
    ax2.text(0.05, 0.08, f"p = {fmt_p(p_b)}",
             transform=ax2.transAxes, fontsize=8)
    style(ax2, "B.  ccRCC — Depth Q3+Q4 vs Q1+Q2")

    # ── Panel C: TI median split ──────────────────────────────────
    if ti_data is not None:
        sub_ti, kmf_h, kmf_l, p_km = ti_data
        if kmf_h is not None:
            kmf_h.plot_survival_function(
                ax=ax3, ci_show=True,
                color=PAL["TI-High (shallow)"], linewidth=1.5)
        if kmf_l is not None:
            kmf_l.plot_survival_function(
                ax=ax3, ci_show=True,
                color=PAL["TI-Low (deep)"], linewidth=1.5)
        if p_km is not None and not np.isnan(p_km):
            ax3.text(0.05, 0.08, f"p = {fmt_p(p_km)}",
                     transform=ax3.transAxes, fontsize=8)
    style(ax3,
          "C.  ccRCC — GOT1/RUNX1 TI\n     (median split)",
          yl="Survival prob.")

    # ── Panel D: depth vs OS scatter ─────────────────────────────
    sub_sc = master.dropna(
        subset=["depth","os_time","os_event"])
    sc = ax4.scatter(sub_sc["depth"], sub_sc["os_time"],
                     c=sub_sc["os_event"], cmap="RdYlGn_r",
                     alpha=0.35, s=12, linewidths=0)
    plt.colorbar(sc, ax=ax4, label="Event (1=death)", shrink=0.75)
    try:
        z = np.polyfit(sub_sc["depth"].values,
                       sub_sc["os_time"].values, 1)
        px = np.linspace(sub_sc["depth"].min(),
                         sub_sc["depth"].max(), 100)
        ax4.plot(px, np.polyval(z, px), "k--", lw=1, alpha=0.6)
    except Exception:
        pass
    r_val, p_val = stats.pearsonr(sub_sc["depth"].values,
                                   sub_sc["os_time"].values)
    ax4.text(0.05, 0.93, f"r={fmt(r_val)}  p={fmt_p(p_val)}",
             transform=ax4.transAxes, fontsize=8)
    ax4.set_title("D.  ccRCC — Depth vs OS scatter",
                  fontsize=9, fontweight="bold", pad=5)
    ax4.set_xlabel("Attractor depth score", fontsize=8)
    ax4.set_ylabel("OS time (days)", fontsize=8)
    ax4.tick_params(labelsize=7)
    ax4.spines["top"].set_visible(False)
    ax4.spines["right"].set_visible(False)

    fig.suptitle(
        "ccRCC Attractor Depth and GOT1/RUNX1 Transition Index"
        " vs Overall Survival\n"
        "OrganismCore — Eric Robert Lawson | 2026-03-07"
        " | Document 94f-S4",
        fontsize=10, fontweight="bold", y=0.997)

    plt.savefig(FIG_FILE, dpi=150, bbox_inches="tight")
    plt.close()
    log(f"  Figure saved: {FIG_FILE}")

# ══════════════════════════════════════════════════════════════════
# SCORECARD
# ══════════════════════════════════════════════════════════════════

def scorecard(km_rows, cox_rows, ti_rows, prcc_rows):
    log("")
    log("=" * 65)
    log("PREDICTION SCORECARD")
    log("=" * 65)

    def find_v(rows, key):
        for r in rows:
            if isinstance(r, dict) and r.get("predictor") == key:
                return r.get("verdict", "PENDING")
        return "PENDING"

    # S4-P1
    v1 = find_v(ti_rows, "TI_neg_continuous")
    log(f"  S4-P1  TI (−GOT1/RUNX1) HR>1.5 p<0.05          {v1}")

    # S4-P2
    v2a = next((r.get("verdict","?") for r in km_rows
                if r.get("analysis") == "Q4_vs_Q1"), "PENDING")
    v2b = next((r.get("verdict","?") for r in km_rows
                if r.get("analysis") == "Q3Q4_vs_Q1Q2"), "PENDING")
    log(f"  S4-P2a Q4 < Q1  (log-rank)                       {v2a}")
    log(f"  S4-P2b Q3+Q4 < Q1+Q2  (log-rank)                 {v2b}")

    # S4-P3
    c_s2 = next((r["c_index"] for r in cox_rows
                 if r.get("predictor") == "S2_depth"), np.nan)
    c_s5 = next((r["c_index"] for r in cox_rows
                 if r.get("predictor") == "S5_depth"), np.nan)
    if not np.isnan(c_s2) and not np.isnan(c_s5):
        v3 = "CONFIRMED ✓" if c_s5 > c_s2 else "NOT CONFIRMED ✗"
        log(f"  S4-P3  S5 depth > S2 depth (C-index)             "
            f"{v3}  ({fmt(c_s5)} vs {fmt(c_s2)})")
    else:
        log(f"  S4-P3  S5 depth > S2 depth (C-index)             PENDING")

    # S4-P4
    confirmed_p4 = [r["gene"] for r in []
                    if r.get("verdict") == "CONFIRMED ✓"]
    log(f"  S4-P4  Gene OS panel — see gene_os_s4.csv")

    # S4-P5
    v5a = next((r.get("verdict","?") for r in prcc_rows
                if r.get("analysis") == "FA1_TI_pooled"), "PENDING")
    v5c = next((r.get("verdict","?") for r in prcc_rows
                if r.get("analysis") == "Type1_vs_Type2"), "PENDING")
    log(f"  S4-P5a FA-1 TI predicts OS (pooled)               {v5a}")
    log(f"  S4-P5b Type 2 OS < Type 1 OS                      {v5c}")

# ══════════════════════════════════════════════════════════════════
# SAVE OUTPUTS
# ══════════════════════════════════════════════════════════════════

def save_outputs(master, cox_rows, ti_rows, gene_rows, prcc_rows):
    log("")
    log("=" * 65)
    log("SAVING OUTPUTS")
    log("=" * 65)

    master.to_csv(SURV_CSV)
    log(f"  survival_s4.csv  → {SURV_CSV}")

    all_cox = cox_rows + [r for r in ti_rows
                          if isinstance(r, dict) and "hr" in r]
    if all_cox:
        pd.DataFrame(all_cox).to_csv(COX_CSV, index=False)
        log(f"  cox_s4.csv       → {COX_CSV}")

    if gene_rows:
        pd.DataFrame(gene_rows).to_csv(GENE_OS_CSV, index=False)
        log(f"  gene_os_s4.csv   → {GENE_OS_CSV}")

    if prcc_rows:
        pd.DataFrame(prcc_rows).to_csv(PRCC_SURV_CSV, index=False)
        log(f"  prcc_surv_s4.csv → {PRCC_SURV_CSV}")

    log(f"  s4_log.txt       → {LOG_FILE}")

# ═══════════════════════��══════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════

def main():
    log("OrganismCore — RCC Script 4")
    log("ccRCC + PRCC Survival Analysis")
    log(f"Document 94f-S4 | {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    log("")
    log("Predictions locked 2026-03-07 before run:")
    log("  S4-P1: TI (GOT1/RUNX1) → OS  HR>1.5  p<0.05")
    log("  S4-P2: Q4 < Q1 and Q3+Q4 < Q1+Q2  p<0.05")
    log("  S4-P3: S5 depth C > S2 depth C")
    log("  S4-P4: LOXL2/TGFBI/RUNX1 high=worse  "
        "OGDHL/SLC13A2 low=worse")
    log("  S4-P5: FA-1 TI → OS pooled + Type1/2 separation")

    try:
        master      = load_ccrcc()
        merged_expr = load_ccrcc_expr()

        km_rows, km_obj = obj1_km_quartiles(master)
        cox_rows        = obj2_cox_continuous(master)
        ti_rows, ti_obj = obj3_ti_survival(master)
        gene_rows       = obj6_gene_os(merged_expr)
        _               = obj7_multivariate(master, merged_expr)
        prcc_rows       = obj8_prcc()

        make_figure(master, ti_obj, PRCC_SURV_CSV)
        save_outputs(master, cox_rows, ti_rows, gene_rows, prcc_rows)
        scorecard(km_rows, cox_rows, ti_rows, prcc_rows)

        log("")
        log("=" * 65)
        log("SCRIPT 4 COMPLETE")
        log(f"  Results: {S4_DIR}")
        log(f"  Figure:  {FIG_FILE}")
        log(f"  Log:     {LOG_FILE}")
        log("=" * 65)

    except Exception as e:
        log(f"\nFATAL: {e}")
        log(traceback.format_exc())
    finally:
        flush_log()


if __name__ == "__main__":
    main()
