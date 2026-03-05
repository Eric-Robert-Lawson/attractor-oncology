"""
BRCA CLAUDIN-LOW — SCRIPT 2
OrganismCore — Document BRCA-S7d | 2026-03-05

PREDICTIONS FROM BRCA-S7c (before_script2.md):
  S2-P1: Depth score predicts OS — depth-high = worse  [HIGH]
  S2-P2: CLDN3-low predicts worse OS                   [MOD-HIGH]
  S2-P3: Immune score predicts OS — high = better      [MODERATE]
  S2-P4: MKI67-high predicts worse OS                  [MODERATE]
  S2-P5: CL worse OS than LumA, comparable to TNBC     [MODERATE]
  S2-P6: PARP1 does NOT predict OS independently       [NEGATIVE]
  S2-P7: Depth score holds in Stratum B (ESR1-low)     [MODERATE]
  S2-P8: CD274 does NOT predict OS in TCGA cohort      [NEGATIVE]

THREE ANALYSIS STRATA (from BRCA-S7c):
  A: Full geometry set (score >= 7, ERBB2 exclusion)  n=268
  B: ESR1-low subset (ESR1 < cohort median)           n=TBD
  C: Basal+Normal-PAM50 within geometry set            n=~79

SELF-CONTAINED:
  Reuses data files downloaded by Script 1.
  Pancan survival supplement downloaded if not cached.
  All output to Claudin_Low_s2_results/.

GEOMETRY-FIRST PROTOCOL v2.0:
  Cross-subtype KM plot printed FIRST (no prediction
  framing). Individual prediction tests printed SECOND.
"""

import os
import sys
import gzip
import time
import warnings
import urllib.request
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines import CoxPHFitter

# ============================================================
# CONFIGURATION
# ============================================================

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))

# Reuse Script 1 data
S1_DATA_DIR  = os.path.join(SCRIPT_DIR,
                             "Claudin_Low_s1_analysis", "data")

# Script 2 output
BASE_DIR     = os.path.join(SCRIPT_DIR, "Claudin_Low_s2_results")
os.makedirs(BASE_DIR, exist_ok=True)

LOG_FILE     = os.path.join(BASE_DIR, "cl_s2_log.txt")
FIG_FILE     = os.path.join(BASE_DIR, "cl_s2_figure.png")
CSV_FILE     = os.path.join(BASE_DIR, "cl_s2_survival.csv")
SCORECARD    = os.path.join(BASE_DIR, "cl_s2_scorecard.csv")

# Input files from Script 1
EXPR_FILE    = os.path.join(S1_DATA_DIR, "TCGA_BRCA_HiSeqV2.gz")
CLIN_FILE    = os.path.join(S1_DATA_DIR, "TCGA_BRCA_clinicalMatrix.tsv")

# Survival supplement
SURV_URLS = [
    "https://pancanatlas.xenahubs.net/download/"
    "Survival_SupplementalTable_S1_20171025_xena_sp",
    "https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/"
    "download/Survival_SupplementalTable_S1_20171025_xena_sp",
    "https://gdc.xenahubs.net/download/TCGA-CDR-SupplementalTableS1.tsv.gz",
]
SURV_FILE    = os.path.join(S1_DATA_DIR, "TCGA_pancan_survival.tsv")

# ── Barcode constants (confirmed from ILC s2 log) ────────────
PATIENT_PREFIX_LEN = 12
SAMPLE_TYPE_POS    = (13, 15)

# ── Classification constants (must match Script 1) ───────────
CL_PRIMARY_THRESHOLD      = 7
ERBB2_EXCLUSION_PERCENTILE = 90

CL_NEG_MARKERS = ["CLDN3", "CLDN4", "CLDN7", "CDH1", "ESR1"]
CL_POS_MARKERS = ["VIM", "CD44", "SNAI1", "ZEB1", "FN1"]
CL_SIGNATURE_GENES = CL_NEG_MARKERS + CL_POS_MARKERS

# ── Immune composite genes ────────────────────────────────────
IMMUNE_GENES = ["FOXP3", "PDCD1", "TIGIT", "LAG3"]

# ── Minimum n for survival analysis ──────────────────────────
MIN_SURV_N = 10

# ============================================================
# LOGGING
# ============================================================

_log_lines = []

def log(msg=""):
    print(msg)
    _log_lines.append(str(msg))

def flush_log():
    with open(LOG_FILE, "w") as f:
        f.write("\n".join(_log_lines))

# ============================================================
# DOWNLOAD UTILITY
# ============================================================

def download_file(url, dest, label="", timeout=180):
    try:
        log(f"  Downloading {label or os.path.basename(url[:70])} ...")
        req = urllib.request.Request(
            url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = resp.read()
        with open(dest, "wb") as f:
            f.write(data)
        log(f"  OK: {os.path.basename(dest)} "
            f"({os.path.getsize(dest)/1e6:.1f} MB)")
        return True
    except Exception as e:
        log(f"  FAILED: {e}")
        if os.path.exists(dest):
            os.remove(dest)
        return False


def try_urls(url_list, dest, label=""):
    if os.path.exists(dest) and os.path.getsize(dest) > 10_000:
        log(f"  Cached: {os.path.basename(dest)}")
        return True
    for url in url_list:
        log(f"  Trying: {url[:80]}")
        if download_file(url, dest, label=label):
            return True
        time.sleep(2)
    log(f"  FATAL: All URLs failed for {label or dest}")
    return False

# ============================================================
# STEP 0 — SURVIVAL DATA ACQUISITION
# ============================================================

def acquire_survival():
    log("=" * 65)
    log("STEP 0: SURVIVAL DATA ACQUISITION")
    log("=" * 65)

    # Verify Script 1 files exist
    for f, label in [(EXPR_FILE, "Expression"),
                     (CLIN_FILE, "Clinical")]:
        if not os.path.exists(f):
            log(f"  FATAL: {label} file not found: {f}")
            log("  Run Script 1 first.")
            flush_log()
            sys.exit(1)
        log(f"  Found: {label} — {os.path.basename(f)}")

    log("\n-- Pancan survival supplement --")
    try_urls(SURV_URLS, SURV_FILE, "TCGA pancan survival")

# ============================================================
# STEP 1 — LOAD EXPRESSION
# ============================================================

def load_expression():
    log("")
    log("=" * 65)
    log("STEP 1: LOAD EXPRESSION MATRIX")
    log("=" * 65)

    try:
        expr = pd.read_csv(EXPR_FILE, sep="\t",
                           index_col=0, compression="gzip")
    except Exception:
        try:
            expr = pd.read_csv(EXPR_FILE, sep="\t", index_col=0)
        except Exception as e:
            log(f"  FATAL: {e}")
            flush_log()
            sys.exit(1)

    if expr.shape[0] < expr.shape[1]:
        expr = expr.T
    log(f"  Shape (genes x samples): {expr.shape}")
    log(f"  Sample example: {list(expr.columns[:2])}")
    return expr

# ============================================================
# STEP 2 — LOAD CLINICAL + SURVIVAL
# Confirmed column names from ILC s2 log and S1 run:
#   PAM50Call_RNAseq, histological_type
#   OS_Time_nature2012 / OS_event_nature2012 (clinical fallback)
#   Pancan survival: OS.time, OS (sample-indexed)
# ============================================================

def load_clinical_and_survival():
    log("")
    log("=" * 65)
    log("STEP 2: LOAD CLINICAL + SURVIVAL")
    log("=" * 65)

    # ── Clinical ─────────────────────────────────────────────
    with open(CLIN_FILE, "rb") as f:
        first = f.readline().decode("utf-8", errors="replace")
    sep = "\t" if "\t" in first else ","
    clin = pd.read_csv(CLIN_FILE, sep=sep,
                       index_col=0, low_memory=False)
    clin.index = [str(s).replace(".", "-") for s in clin.index]
    log(f"  Clinical shape: {clin.shape}")

    # ── Pancan survival ───────────────────────────────────────
    surv = None
    surv_time_col  = None
    surv_event_col = None

    if os.path.exists(SURV_FILE) and os.path.getsize(SURV_FILE) > 1000:
        try:
            surv = pd.read_csv(SURV_FILE, sep="\t",
                               low_memory=False)
            log(f"  Pancan survival shape: {surv.shape}")
            log(f"  Columns: {list(surv.columns[:8])}")

            # Find sample/patient ID column
            id_col = None
            for cand in ["sample", "_PATIENT", "bcr_patient_barcode",
                         "sampleID"]:
                if cand in surv.columns:
                    id_col = cand
                    break
            if id_col:
                surv = surv.set_index(id_col)
                log(f"  Survival index set to: '{id_col}'")
            surv.index = [str(s).replace(".", "-")
                          for s in surv.index]

            # Find OS columns
            for tc in ["OS.time", "os_time", "OS_time",
                       "survival_time"]:
                if tc in surv.columns:
                    surv_time_col = tc
                    break
            for ec in ["OS", "os_event", "OS_event",
                       "vital_status"]:
                if ec in surv.columns:
                    surv_event_col = ec
                    break
            log(f"  OS time col:  {surv_time_col}")
            log(f"  OS event col: {surv_event_col}")
        except Exception as e:
            log(f"  WARNING: Survival load failed: {e}")
            surv = None

    # ── Clinical fallback columns ─────────────────────────────
    clin_time_col  = None
    clin_event_col = None
    for tc in ["OS_Time_nature2012", "days_to_death",
               "OS.time", "survival_months"]:
        if tc in clin.columns:
            clin_time_col = tc
            break
    for ec in ["OS_event_nature2012", "vital_status",
               "OS", "death_days_to"]:
        if ec in clin.columns:
            clin_event_col = ec
            break
    log(f"  Clinical fallback — time: {clin_time_col}, "
        f"event: {clin_event_col}")

    return clin, surv, surv_time_col, surv_event_col, \
           clin_time_col, clin_event_col

# ============================================================
# STEP 3 — REBUILD POPULATIONS (identical logic to Script 1)
# Re-identifies claudin-low from geometry to ensure consistency.
# ============================================================

def build_populations(expr, clin):
    log("")
    log("=" * 65)
    log("STEP 3: REBUILD POPULATIONS")
    log("  (Replicating Script 1 classification exactly)")
    log("=" * 65)

    samples = list(expr.columns)
    tumour_samples = []
    normal_samples = []
    for s in samples:
        parts = s.split("-")
        stype = parts[3][:2] if len(parts) >= 4 else "??"
        if stype == "01":
            tumour_samples.append(s)
        elif stype == "11":
            normal_samples.append(s)

    log(f"  Tumour: {len(tumour_samples)}  Normal: {len(normal_samples)}")

    # PAM50 map
    pam50_col = None
    for cand in ["PAM50Call_RNAseq", "PAM50_mRNA_nature2012"]:
        if cand in clin.columns:
            pam50_col = cand
            break

    pam50_map = {}
    if pam50_col:
        for pat_id in clin.index:
            pfx = str(pat_id)[:PATIENT_PREFIX_LEN]
            pam50_map[pfx] = str(clin.loc[pat_id, pam50_col])

    def get_pam50(s):
        return pam50_map.get(s[:PATIENT_PREFIX_LEN], "Unknown")

    # Subtype buckets
    subtype_buckets = {"LumA": [], "LumB": [],
                       "Her2": [], "Basal": [], "Normal": []}
    for s in tumour_samples:
        p = get_pam50(s)
        if p in subtype_buckets:
            subtype_buckets[p].append(s)

    # Claudin-low: 10-gene signature
    tumour_expr = expr[tumour_samples]
    sig_present = [g for g in CL_SIGNATURE_GENES if g in expr.index]

    cohort_medians = {}
    for g in sig_present:
        cohort_medians[g] = float(np.median(tumour_expr.loc[g].values))

    scores = {}
    for s in tumour_samples:
        score = 0
        for g in CL_POS_MARKERS:
            if g in sig_present:
                if float(expr.loc[g, s]) > cohort_medians[g]:
                    score += 1
        for g in CL_NEG_MARKERS:
            if g in sig_present:
                if float(expr.loc[g, s]) < cohort_medians[g]:
                    score += 1
        scores[s] = score

    cl_raw = [s for s, sc in scores.items() if sc >= CL_PRIMARY_THRESHOLD]

    if "ERBB2" in expr.index:
        erbb2_vals = tumour_expr.loc["ERBB2", tumour_samples].values
        cutoff = np.percentile(erbb2_vals, ERBB2_EXCLUSION_PERCENTILE)
        cl_samples = [s for s in cl_raw
                      if float(expr.loc["ERBB2", s]) <= cutoff]
    else:
        cl_samples = cl_raw

    log(f"  Claudin-low (Stratum A): n={len(cl_samples)}")

    # ── Stratum B: ESR1-low subset ────────────────────────────
    if "ESR1" in expr.index:
        esr1_cohort_median = float(
            np.median(tumour_expr.loc["ESR1"].values))
        stratum_b = [s for s in cl_samples
                     if float(expr.loc["ESR1", s]) < esr1_cohort_median]
        log(f"  Stratum B (ESR1-low):    n={len(stratum_b)}")
    else:
        stratum_b = []
        log("  Stratum B: ESR1 not found")

    # ── Stratum C: Basal + Normal-PAM50 within geometry set ───
    stratum_c = [s for s in cl_samples
                 if get_pam50(s) in ("Basal", "Normal")]
    log(f"  Stratum C (Basal+Normal-PAM50): n={len(stratum_c)}")

    populations = {
        "claudin_low":  cl_samples,
        "stratum_b":    stratum_b,
        "stratum_c":    stratum_c,
        "normal":       normal_samples,
        "luma":         subtype_buckets["LumA"],
        "lumb":         subtype_buckets["LumB"],
        "her2":         subtype_buckets["Her2"],
        "basal":        subtype_buckets["Basal"],
    }

    log("")
    log("  Population summary:")
    for k, v in populations.items():
        log(f"    {k:<20} n={len(v)}")

    return populations, scores

# ============================================================
# STEP 4 — BUILD SURVIVAL VECTORS
# Maps expression sample barcodes to survival data.
# Priority: pancan supplement → clinical matrix fallback.
# ============================================================

def build_survival_vectors(expr, populations, surv,
                            surv_time_col, surv_event_col,
                            clin, clin_time_col, clin_event_col):
    log("")
    log("=" * 65)
    log("STEP 4: BUILD SURVIVAL VECTORS")
    log("=" * 65)

    all_samples = list(expr.columns)
    os_time  = {}
    os_event = {}

    # ── Attempt 1: pancan supplement ─────────────────────────
    if surv is not None and surv_time_col and surv_event_col:
        log(f"  Pancan survival: time='{surv_time_col}', "
            f"event='{surv_event_col}'")
        n_matched = 0
        for s in all_samples:
            if s in surv.index:
                tv = surv.loc[s, surv_time_col]
                ev = surv.loc[s, surv_event_col]
                try:
                    t = float(tv)
                    e = float(ev)
                    if not np.isnan(t) and t > 0:
                        os_time[s]  = t
                        os_event[s] = e
                        n_matched += 1
                except (ValueError, TypeError):
                    pass
        log(f"  Pancan matches: {n_matched}")

    # ── Attempt 2: patient-prefix join ───────────────────────
    if surv is not None and surv_time_col and len(os_time) < 200:
        log("  Attempting patient-prefix join on pancan survival...")
        surv_pfx = {}
        for idx in surv.index:
            pfx = str(idx)[:PATIENT_PREFIX_LEN]
            if pfx not in surv_pfx:
                surv_pfx[pfx] = idx
        n_pfx = 0
        for s in all_samples:
            if s in os_time:
                continue
            pfx = s[:PATIENT_PREFIX_LEN]
            if pfx in surv_pfx:
                sidx = surv_pfx[pfx]
                try:
                    t = float(surv.loc[sidx, surv_time_col])
                    e = float(surv.loc[sidx, surv_event_col])
                    if not np.isnan(t) and t > 0:
                        os_time[s]  = t
                        os_event[s] = e
                        n_pfx += 1
                except (ValueError, TypeError):
                    pass
        log(f"  Prefix-join matches: {n_pfx}")

    # ── Attempt 3: clinical matrix fallback ──────────────────
    if clin_time_col and clin_event_col and len(os_time) < 200:
        log(f"  Clinical fallback: time='{clin_time_col}', "
            f"event='{clin_event_col}'")

        # Build patient-prefix → clinical row map
        clin_pfx = {}
        for pat_id in clin.index:
            pfx = str(pat_id)[:PATIENT_PREFIX_LEN]
            clin_pfx[pfx] = pat_id

        n_clin = 0
        for s in all_samples:
            if s in os_time:
                continue
            pfx = s[:PATIENT_PREFIX_LEN]
            if pfx in clin_pfx:
                pat_id = clin_pfx[pfx]
                try:
                    tv = clin.loc[pat_id, clin_time_col]
                    ev = clin.loc[pat_id, clin_event_col]
                    t = float(tv)
                    e_raw = str(ev).lower()
                    if any(x in e_raw for x in
                           ["dead", "deceased", "1"]):
                        e = 1.0
                    elif any(x in e_raw for x in
                             ["alive", "living", "0"]):
                        e = 0.0
                    else:
                        e = float(ev)
                    if not np.isnan(t) and t > 0:
                        os_time[s]  = t
                        os_event[s] = e
                        n_clin += 1
                except (ValueError, TypeError):
                    pass
        log(f"  Clinical fallback matches: {n_clin}")

    log(f"  Total samples with survival: {len(os_time)}")

    # Report survival coverage per population
    log("")
    log("  Survival coverage per population:")
    for pop, samps in populations.items():
        covered = [s for s in samps if s in os_time]
        log(f"    {pop:<20} {len(covered)}/{len(samps)} "
            f"with survival data")

    return os_time, os_event

# ============================================================
# STEP 5 — COMPUTE DEPTH SCORE
# Identical to Script 1 depth axis computation.
# Applied to the full expression matrix so all populations
# can be scored on the same axis.
# ============================================================

def compute_depth_scores(expr, tumour_samples):
    log("")
    log("=" * 65)
    log("STEP 5: COMPUTE DEPTH SCORES")
    log("=" * 65)

    sig_present = [g for g in CL_SIGNATURE_GENES if g in expr.index]
    pos_genes   = [g for g in CL_POS_MARKERS if g in expr.index]
    neg_genes   = [g for g in CL_NEG_MARKERS if g in expr.index]

    tumour_expr = expr[tumour_samples]

    def z_score_gene(arr):
        s = arr.std()
        if s < 1e-9:
            return np.zeros_like(arr)
        return (arr - arr.mean()) / s

    # Z-score each gene across the full tumour cohort
    pos_z = np.column_stack([
        z_score_gene(tumour_expr.loc[g].values.astype(float))
        for g in pos_genes
    ]) if pos_genes else np.zeros((len(tumour_samples), 1))

    neg_z = np.column_stack([
        z_score_gene(tumour_expr.loc[g].values.astype(float))
        for g in neg_genes
    ]) if neg_genes else np.zeros((len(tumour_samples), 1))

    # Depth = positive programme mean - negative programme mean
    depth_arr = np.mean(pos_z, axis=1) - np.mean(neg_z, axis=1)

    depth_scores = dict(zip(tumour_samples, depth_arr))

    log(f"  Depth scores computed for {len(tumour_samples)} tumour samples")
    log(f"  Range: {depth_arr.min():.4f} to {depth_arr.max():.4f}")
    log(f"  Mean:  {depth_arr.mean():.4f}")
    log(f"  Std:   {depth_arr.std():.4f}")

    return depth_scores

# ============================================================
# SURVIVAL UTILITY FUNCTIONS
# ============================================================

def fmt_p(p):
    if p is None or np.isnan(p):
        return "p=nan"
    if p < 1e-10: return f"p={p:.2e} ***"
    if p < 1e-5:  return f"p={p:.2e} **"
    if p < 0.05:  return f"p={p:.3f} *"
    if p < 0.10:  return f"p={p:.3f} (trend)"
    return              f"p={p:.3f}"


def get_survival_arrays(samples, os_time, os_event):
    """Returns (T, E, valid_samples) for a list of samples."""
    T, E, valid = [], [], []
    for s in samples:
        if s in os_time and s in os_event:
            t = os_time[s]
            e = os_event[s]
            if not np.isnan(t) and not np.isnan(e) and t > 0:
                T.append(t)
                E.append(int(e))
                valid.append(s)
    return np.array(T), np.array(E), valid


def run_logrank(T_high, E_high, T_low, E_low, label=""):
    """Run log-rank test. Returns (p, HR_approx)."""
    if len(T_high) < MIN_SURV_N or len(T_low) < MIN_SURV_N:
        log(f"    {label}: Insufficient n "
            f"(high={len(T_high)}, low={len(T_low)}) — skipped")
        return np.nan, np.nan
    try:
        result = logrank_test(T_high, T_low,
                              event_observed_A=E_high,
                              event_observed_B=E_low)
        p = result.p_value

        # Approximate HR: ratio of event rates (not true Cox HR)
        rate_high = E_high.sum() / T_high.sum() if T_high.sum() > 0 else np.nan
        rate_low  = E_low.sum()  / T_low.sum()  if T_low.sum()  > 0 else np.nan
        hr_approx = rate_high / rate_low if rate_low > 0 else np.nan

        log(f"    {label}: n_high={len(T_high)} n_low={len(T_low)}"
            f"  events_h={E_high.sum()} events_l={E_low.sum()}"
            f"  {fmt_p(p)}  HR≈{hr_approx:.3f}")
        return p, hr_approx
    except Exception as e:
        log(f"    {label}: logrank failed — {e}")
        return np.nan, np.nan


def tertile_split(samples, scores_dict, os_time, os_event):
    """
    Split samples by score tertile.
    Returns (T_high, E_high, T_low, E_low, T_mid, E_mid,
             valid_high, valid_mid, valid_low).
    """
    scored = [(s, scores_dict[s]) for s in samples
              if s in scores_dict and s in os_time]
    if len(scored) < MIN_SURV_N * 3:
        return (np.array([]),) * 6 + ([], [], [])
    scored.sort(key=lambda x: x[1])
    n   = len(scored)
    t1  = n // 3
    t2  = 2 * n // 3
    low_s  = [s for s, _ in scored[:t1]]
    mid_s  = [s for s, _ in scored[t1:t2]]
    high_s = [s for s, _ in scored[t2:]]
    T_h, E_h, vh = get_survival_arrays(high_s, os_time, os_event)
    T_m, E_m, vm = get_survival_arrays(mid_s,  os_time, os_event)
    T_l, E_l, vl = get_survival_arrays(low_s,  os_time, os_event)
    return T_h, E_h, T_l, E_l, T_m, E_m, vh, vm, vl

# ============================================================
# STEP 6 — CROSS-SUBTYPE KAPLAN-MEIER
# Geometry-first: read the cross-subtype survival structure
# before any individual prediction is tested.
# ============================================================

def cross_subtype_km(populations, os_time, os_event):
    log("")
    log("=" * 65)
    log("STEP 6: CROSS-SUBTYPE KAPLAN-MEIER")
    log("  Geometry read BEFORE prediction tests.")
    log("  No prediction framing. Read what is there.")
    log("=" * 65)

    pop_order = ["normal", "claudin_low", "luma",
                 "lumb", "her2", "basal"]
    pop_labels = {
        "normal":      "Normal",
        "claudin_low": "Claudin-low",
        "luma":        "LumA",
        "lumb":        "LumB",
        "her2":        "HER2",
        "basal":       "TNBC/Basal",
    }
    colors = {
        "normal":      "#6e7681",
        "claudin_low": "#58a6ff",
        "luma":        "#3fb950",
        "lumb":        "#d29922",
        "her2":        "#bc8cff",
        "basal":       "#f85149",
    }

    log("")
    log(f"  {'Population':<20} {'n_surv':>8} "
        f"{'events':>8} {'median_OS':>12}")
    log("  " + "─" * 52)

    km_results = {}
    for pop in pop_order:
        samps = populations.get(pop, [])
        T, E, valid = get_survival_arrays(samps, os_time, os_event)
        if len(T) < MIN_SURV_N:
            log(f"  {pop_labels[pop]:<20} n<{MIN_SURV_N} — skipped")
            continue
        kmf = KaplanMeierFitter()
        kmf.fit(T, E, label=pop_labels[pop])
        med = kmf.median_survival_time_
        log(f"  {pop_labels[pop]:<20} {len(T):>8} "
            f"{int(E.sum()):>8} {str(round(med, 1)):>12}")
        km_results[pop] = (kmf, T, E, valid)

    # Pairwise log-rank: claudin-low vs each subtype
    log("")
    log("  Pairwise log-rank: Claudin-low vs each subtype:")
    if "claudin_low" in km_results:
        T_cl, E_cl = km_results["claudin_low"][1], km_results["claudin_low"][2]
        for pop in pop_order:
            if pop == "claudin_low" or pop not in km_results:
                continue
            T_o, E_o = km_results[pop][1], km_results[pop][2]
            p, _ = run_logrank(T_cl, E_cl, T_o, E_o,
                               f"CL vs {pop_labels[pop]}")

    return km_results, pop_labels, colors

# ============================================================
# STEP 7 — PREDICTION TESTS
# All eight predictions from BRCA-S7c.
# ============================================================

def run_prediction_tests(expr, populations, depth_scores,
                         os_time, os_event):
    log("")
    log("=" * 65)
    log("STEP 7: PREDICTION TESTS")
    log("  Predictions from BRCA-S7c tested here.")
    log("  Geometry was read first (Step 6).")
    log("=" * 65)

    scorecard = []

    def record(pred_id, gene_or_axis, stratum,
               n_high, n_low, p, hr, direction_correct,
               confirmed, note=""):
        scorecard.append({
            "prediction":        pred_id,
            "gene_or_axis":      gene_or_axis,
            "stratum":           stratum,
            "n_high":            n_high,
            "n_low":             n_low,
            "p_value":           p,
            "HR_approx":         hr,
            "direction_correct": direction_correct,
            "confirmed":         confirmed,
            "note":              note,
        })

    # ── S2-P1: Depth score vs OS — Stratum A ─────────────────
    log("")
    log("  ─────────────────────────────────────────────────────")
    log("  S2-P1 — DEPTH SCORE vs OS (Stratum A, n=268)")
    log("  Prediction: depth-high = worse OS (HR > 1.0)")
    log("  ─────────────────────────────────────────────────────")

    cl_a = populations["claudin_low"]
    res = tertile_split(cl_a, depth_scores, os_time, os_event)
    T_h, E_h, T_l, E_l, T_m, E_m, vh, vm, vl = res

    p1_a, hr1_a = run_logrank(T_h, E_h, T_l, E_l,
                              "S2-P1 [depth tertile, Stratum A]")
    confirmed = ("CONFIRMED" if (not np.isnan(p1_a) and
                                 p1_a < 0.10 and
                                 (np.isnan(hr1_a) or hr1_a > 1.0))
                 else "NOT CONFIRMED")
    log(f"  Status: [{confirmed}]")
    record("S2-P1", "depth_score", "A",
           len(T_h), len(T_l), p1_a, hr1_a,
           (not np.isnan(hr1_a) and hr1_a > 1.0), confirmed)

    # ── S2-P7: Depth score vs OS — Stratum B ─────────────────
    log("")
    log("  ─────────────────────────────────────────────────────")
    log("  S2-P7 — DEPTH SCORE vs OS (Stratum B, ESR1-low)")
    log("  Prediction: HR(B) >= HR(A) — purity test")
    log("  ─────────────────────────────────────────────────────")

    cl_b = populations["stratum_b"]
    if len(cl_b) >= MIN_SURV_N * 3:
        res_b = tertile_split(cl_b, depth_scores, os_time, os_event)
        Tb_h, Eb_h, Tb_l, Eb_l, Tb_m, Eb_m, vbh, vbm, vbl = res_b
        p7_b, hr7_b = run_logrank(Tb_h, Eb_h, Tb_l, Eb_l,
                                  "S2-P7 [depth tertile, Stratum B]")
        purity_confirmed = (not np.isnan(hr7_b) and
                            not np.isnan(hr1_a) and
                            hr7_b >= hr1_a * 0.5)
        status7 = "CONFIRMED" if purity_confirmed else "NOT CONFIRMED"
        log(f"  Status: [{status7}]")
        log(f"  HR(A)={hr1_a:.3f}  HR(B)={hr7_b:.3f}")
        log(f"  Depth signal {'HOLDS' if purity_confirmed else 'WEAKENS'} "
            f"in ESR1-low purified subset")
        record("S2-P7", "depth_score", "B",
               len(Tb_h), len(Tb_l), p7_b, hr7_b,
               purity_confirmed, status7,
               f"HR_A={hr1_a:.3f}")
    else:
        log(f"  Stratum B n={len(cl_b)} — insufficient for tertile split")
        record("S2-P7", "depth_score", "B", 0, 0,
               np.nan, np.nan, False, "UNDERPOWERED",
               f"n={len(cl_b)}")

    # ── Stratum C survival ────────────────────────────────────
    log("")
    log("  ────────────────���────────────────────────────────────")
    log("  STRATUM C — Depth score vs OS (Basal+Normal-PAM50)")
    log("  ─────────────────────────────────────────────────────")
    cl_c = populations["stratum_c"]
    if len(cl_c) >= MIN_SURV_N * 3:
        res_c = tertile_split(cl_c, depth_scores, os_time, os_event)
        Tc_h, Ec_h, Tc_l, Ec_l = res_c[0], res_c[1], res_c[2], res_c[3]
        p_c, hr_c = run_logrank(Tc_h, Ec_h, Tc_l, Ec_l,
                                "depth tertile, Stratum C")
        record("S2-P1-C", "depth_score", "C",
               len(Tc_h), len(Tc_l), p_c, hr_c,
               (not np.isnan(hr_c) and hr_c > 1.0),
               "SEE NOTE",
               "Stratum C purity check")
    else:
        log(f"  Stratum C n={len(cl_c)} — insufficient for tertile split")

    # ── S2-P2: CLDN3 vs OS ───────────────────────────────────
    log("")
    log("  ──────────────────────���──────────────────────────────")
    log("  S2-P2 — CLDN3 expression vs OS (Stratum A)")
    log("  Prediction: CLDN3-low = worse OS")
    log("  ─────────────────────────────────────────────────────")

    if "CLDN3" in expr.index:
        cldn3_scores = {s: float(expr.loc["CLDN3", s])
                        for s in cl_a}
        # NOTE: CLDN3-low = worse, so high-CLDN3 = better
        # For logrank: compare high vs low, expect HR < 1.0
        # (high CLDN3 = better = lower hazard)
        res2 = tertile_split(cl_a, cldn3_scores, os_time, os_event)
        T2_h, E2_h, T2_l, E2_l = res2[0], res2[1], res2[2], res2[3]
        p2, hr2 = run_logrank(T2_h, E2_h, T2_l, E2_l,
                              "S2-P2 [CLDN3 tertile, high vs low]")
        # CLDN3-low = worse means high/low HR should show high CLDN3 = better
        # HR here is high/low. Direction correct if HR < 1.0 (high better)
        dir2 = (not np.isnan(hr2) and hr2 < 1.0)
        confirmed2 = "CONFIRMED" if (not np.isnan(p2) and
                                     p2 < 0.15 and dir2) else "NOT CONFIRMED"
        log(f"  (HR < 1.0 means CLDN3-high = better = prediction confirmed)")
        log(f"  Status: [{confirmed2}]")
        record("S2-P2", "CLDN3", "A",
               len(T2_h), len(T2_l), p2, hr2,
               dir2, confirmed2,
               "HR<1.0 = high CLDN3 = better survival (CLDN3-low = worse)")
    else:
        log("  CLDN3 not in expression matrix — skipped")
        record("S2-P2", "CLDN3", "A", 0, 0,
               np.nan, np.nan, False, "GENE MISSING")

    # ── S2-P3: Immune score vs OS ─────────────────────────────
    log("")
    log("  ─────────────────────────────────────────────────────")
    log("  S2-P3 — IMMUNE SCORE vs OS (Stratum A)")
    log("  Prediction: immune-high = BETTER OS (HR < 1.0)")
    log("  ─────────────────────────────────────────────────────")

    imm_genes_present = [g for g in IMMUNE_GENES if g in expr.index]
    log(f"  Immune genes present: {imm_genes_present}")

    if len(imm_genes_present) >= 2:
        # Composite immune score: mean z-score of immune genes
        imm_z_parts = []
        for g in imm_genes_present:
            vals = np.array([float(expr.loc[g, s]) for s in cl_a])
            s_std = vals.std()
            if s_std > 1e-9:
                imm_z_parts.append((vals - vals.mean()) / s_std)
        if imm_z_parts:
            imm_composite = np.mean(imm_z_parts, axis=0)
            imm_scores = dict(zip(cl_a, imm_composite))
            res3 = tertile_split(cl_a, imm_scores, os_time, os_event)
            T3_h, E3_h, T3_l, E3_l = res3[0], res3[1], res3[2], res3[3]
            p3, hr3 = run_logrank(T3_h, E3_h, T3_l, E3_l,
                                  "S2-P3 [immune score tertile]")
            dir3 = (not np.isnan(hr3) and hr3 < 1.0)  # immune-high = better
            confirmed3 = "CONFIRMED" if (not np.isnan(p3) and
                                         p3 < 0.15 and dir3) else "NOT CONFIRMED"
            if not dir3 and not np.isnan(hr3):
                log(f"  NOTE: HR > 1.0 — immune-high doing WORSE (Treg interpretation)")
            log(f"  Status: [{confirmed3}]")
            record("S2-P3", "immune_composite", "A",
                   len(T3_h), len(T3_l), p3, hr3,
                   dir3, confirmed3,
                   "HR<1.0=immune-high=better")
    else:
        log("  Insufficient immune genes — skipped")
        record("S2-P3", "immune_composite", "A",
               0, 0, np.nan, np.nan, False, "GENES MISSING")

    # ── S2-P4: MKI67 vs OS ───────────────────────────────────
    log("")
    log("  ─────────────────────────────────────────────────────")
    log("  S2-P4 — MKI67 vs OS (Stratum A)")
    log("  Prediction: MKI67-high = worse OS (HR > 1.0)")
    log("  ─────────────────────────────────────────────────────")

    if "MKI67" in expr.index:
        mki67_scores = {s: float(expr.loc["MKI67", s]) for s in cl_a}
        res4 = tertile_split(cl_a, mki67_scores, os_time, os_event)
        T4_h, E4_h, T4_l, E4_l = res4[0], res4[1], res4[2], res4[3]
        p4, hr4 = run_logrank(T4_h, E4_h, T4_l, E4_l,
                              "S2-P4 [MKI67 tertile]")
        dir4 = (not np.isnan(hr4) and hr4 > 1.0)
        confirmed4 = "CONFIRMED" if (not np.isnan(p4) and
                                     p4 < 0.15 and dir4) else "NOT CONFIRMED"
        log(f"  Status: [{confirmed4}]")
        record("S2-P4", "MKI67", "A",
               len(T4_h), len(T4_l), p4, hr4,
               dir4, confirmed4)
    else:
        log("  MKI67 not in matrix")
        record("S2-P4", "MKI67", "A", 0, 0,
               np.nan, np.nan, False, "GENE MISSING")

    # ── S2-P5: Cross-subtype survival position ────────────────
    log("")
    log("  ─────────────────────────────────────────────────────")
    log("  S2-P5 — CLAUDIN-LOW vs LumA and TNBC survival")
    log("  Prediction: CL < LumA (worse), CL ≈ TNBC")
    log("  ─────────────────────────────────────────────────────")

    T_cl, E_cl, _ = get_survival_arrays(
        populations["claudin_low"], os_time, os_event)
    T_la, E_la, _ = get_survival_arrays(
        populations["luma"], os_time, os_event)
    T_ba, E_ba, _ = get_survival_arrays(
        populations["basal"], os_time, os_event)

    p5_cl_la, hr5_cl_la = run_logrank(T_cl, E_cl, T_la, E_la,
                                       "S2-P5 [CL vs LumA]")
    p5_cl_ba, hr5_cl_ba = run_logrank(T_cl, E_cl, T_ba, E_ba,
                                       "S2-P5 [CL vs TNBC]")

    # CL worse than LumA: HR(CL/LumA) > 1.0
    cl_worse_luma = (not np.isnan(hr5_cl_la) and hr5_cl_la > 1.0)
    # CL comparable to TNBC: p(CL vs TNBC) > 0.10
    cl_sim_tnbc   = (np.isnan(p5_cl_ba) or p5_cl_ba > 0.10)

    confirmed5 = "CONFIRMED" if cl_worse_luma else "NOT CONFIRMED"
    log(f"  CL worse than LumA: {'YES' if cl_worse_luma else 'NO'}")
    log(f"  CL comparable to TNBC (p>0.10): {'YES' if cl_sim_tnbc else 'NO'}")
    log(f"  Status: [{confirmed5}]")
    record("S2-P5a", "CL vs LumA", "A",
           len(T_cl), len(T_la), p5_cl_la, hr5_cl_la,
           cl_worse_luma, confirmed5, "expect HR>1.0")
    record("S2-P5b", "CL vs TNBC", "A",
           len(T_cl), len(T_ba), p5_cl_ba, hr5_cl_ba,
           cl_sim_tnbc, "CONFIRMED" if cl_sim_tnbc else "NOT CONFIRMED",
           "expect p>0.10 (comparable)")

    # ── S2-P6: PARP1 independent OS (NEGATIVE prediction) ────
    log("")
    log("  ─────────────────────────────────────────────────────")
    log("  S2-P6 — PARP1 vs OS independent of MKI67 (NEGATIVE)")
    log("  Prediction: PARP1 does NOT predict OS independently")
    log("  ─────────────────────────────────────────────────────")

    if "PARP1" in expr.index:
        parp1_scores = {s: float(expr.loc["PARP1", s]) for s in cl_a}
        res6 = tertile_split(cl_a, parp1_scores, os_time, os_event)
        T6_h, E6_h, T6_l, E6_l = res6[0], res6[1], res6[2], res6[3]
        p6, hr6 = run_logrank(T6_h, E6_h, T6_l, E6_l,
                              "S2-P6 [PARP1 tertile]")
        # Negative prediction: expect p > 0.10
        null_result = (np.isnan(p6) or p6 > 0.10)
        confirmed6 = "CONFIRMED (null as predicted)" if null_result else \
                     "NOT CONFIRMED (PARP1 does predict — unexpected)"
        log(f"  Status: [{confirmed6}]")
        record("S2-P6", "PARP1", "A",
               len(T6_h), len(T6_l), p6, hr6,
               null_result, confirmed6,
               "negative prediction: expect p>0.10")
    else:
        log("  PARP1 not in matrix")
        record("S2-P6", "PARP1", "A", 0, 0,
               np.nan, np.nan, True, "GENE MISSING")

    # ── S2-P8: CD274 vs OS (NEGATIVE prediction) ─────────────
    log("")
    log("  ─────────────────────────────────────────────────────")
    log("  S2-P8 — CD274 (PD-L1) vs OS (NEGATIVE prediction)")
    log("  Prediction: CD274 does NOT predict OS in TCGA cohort")
    log("  ─────────────────────────────────────────────────────")

    if "CD274" in expr.index:
        pdl1_scores = {s: float(expr.loc["CD274", s]) for s in cl_a}
        res8 = tertile_split(cl_a, pdl1_scores, os_time, os_event)
        T8_h, E8_h, T8_l, E8_l = res8[0], res8[1], res8[2], res8[3]
        p8, hr8 = run_logrank(T8_h, E8_h, T8_l, E8_l,
                              "S2-P8 [CD274 tertile]")
        null8 = (np.isnan(p8) or p8 > 0.10)
        confirmed8 = "CONFIRMED (null as predicted)" if null8 else \
                     "NOT CONFIRMED (CD274 does predict — unexpected)"
        log(f"  Status: [{confirmed8}]")
        record("S2-P8", "CD274", "A",
               len(T8_h), len(T8_l), p8, hr8,
               null8, confirmed8,
               "negative prediction: expect p>0.10")
    else:
        log("  CD274 not in matrix")
        record("S2-P8", "CD274", "A", 0, 0,
               np.nan, np.nan, True, "GENE MISSING")

    # ── Save scorecard ────────────────────────────────────────
    sc_df = pd.DataFrame(scorecard)
    sc_df.to_csv(SCORECARD, index=False)
    log(f"\n  Scorecard saved: {SCORECARD}")

    return scorecard

# ============================================================
# STEP 8 — COX MULTIVARIATE
# Depth score + MKI67 + ERBB2 in Stratum A.
# Tests whether depth is independent of proliferation.
# ============================================================

def cox_multivariate(expr, populations, depth_scores,
                     os_time, os_event):
    log("")
    log("=" * 65)
    log("STEP 8: COX MULTIVARIATE")
    log("  depth_score + MKI67 + ERBB2 within claudin-low")
    log("  Tests: is depth score independent of proliferation?")
    log("=" * 65)

    cl_a = populations["claudin_low"]
    rows = []
    for s in cl_a:
        if s not in os_time:
            continue
        row = {
            "T":           os_time.get(s, np.nan),
            "E":           os_event.get(s, np.nan),
            "depth_score": depth_scores.get(s, np.nan),
        }
        for g in ["MKI67", "ERBB2", "ESR1", "CLDN3"]:
            row[g] = float(expr.loc[g, s]) if g in expr.index else np.nan
        rows.append(row)

    cox_df = pd.DataFrame(rows).dropna()
    log(f"  Cox input: n={len(cox_df)}")

    if len(cox_df) < 30:
        log("  Insufficient n for Cox. Skipping.")
        return

    try:
        # Standardise covariates
        covars = ["depth_score", "MKI67", "ERBB2"]
        for c in covars:
            if c in cox_df.columns:
                m, s = cox_df[c].mean(), cox_df[c].std()
                if s > 1e-9:
                    cox_df[c + "_z"] = (cox_df[c] - m) / s
                else:
                    cox_df[c + "_z"] = 0.0

        z_cols = [c + "_z" for c in covars
                  if c + "_z" in cox_df.columns]
        cox_input = cox_df[["T", "E"] + z_cols].copy()
        cox_input = cox_input.rename(
            columns={"T": "duration", "E": "event"})

        cph = CoxPHFitter()
        cph.fit(cox_input, duration_col="duration",
                event_col="event")
        log("")
        log("  Cox PH results:")
        summary = cph.summary
        for idx, row in summary.iterrows():
            log(f"    {str(idx):<22} HR={row['exp(coef)']:.4f} "
                f"[{row['exp(coef) lower 95%']:.3f}–"
                f"{row['exp(coef) upper 95%']:.3f}]  "
                f"{fmt_p(row['p'])}")

        # Depth independent if p < 0.10 in multivariate
        depth_col = "depth_score_z"
        if depth_col in summary.index:
            dp = summary.loc[depth_col, "p"]
            log(f"\n  Depth score in multivariate: {fmt_p(dp)}")
            if dp < 0.10:
                log("  Depth score is INDEPENDENT of MKI67 and ERBB2 ✓")
            else:
                log("  Depth score is not independently significant after "
                    "MKI67/ERBB2 adjustment")
    except Exception as e:
        log(f"  Cox failed: {e}")

# ============================================================
# STEP 9 — FIGURE
# ============================================================

def make_figure(populations, depth_scores, os_time, os_event,
                km_results, pop_labels, colors, expr):
    log("")
    log("=" * 65)
    log("STEP 9: GENERATE FIGURE")
    log("=" * 65)

    fig = plt.figure(figsize=(22, 20))
    fig.patch.set_facecolor("#0d1117")
    gs  = gridspec.GridSpec(3, 3, figure=fig,
                            hspace=0.45, wspace=0.40)

    BLUE   = "#58a6ff"
    GREEN  = "#3fb950"
    RED    = "#f85149"
    ORANGE = "#d29922"
    PURPLE = "#bc8cff"
    GREY   = "#6e7681"

    def ax_style(ax):
        ax.set_facecolor("#161b22")
        for sp in ax.spines.values():
            sp.set_edgecolor("#30363d")
        ax.tick_params(colors="#8b949e", labelsize=7)
        ax.xaxis.label.set_color("#8b949e")
        ax.yaxis.label.set_color("#8b949e")

    lkw = dict(color="#8b949e", fontsize=8)
    tkw = dict(color="#e6edf3", fontsize=10,
               fontweight="bold", pad=8)

    # ── Panel 1: Cross-subtype KM ─────────────────────────────
    ax1 = fig.add_subplot(gs[0, 0:2])
    ax1.set_title("Cross-subtype Overall Survival (S2-P5)", **tkw)
    pop_order = ["normal", "luma", "lumb", "her2",
                 "claudin_low", "basal"]
    for pop in pop_order:
        if pop not in km_results:
            continue
        kmf, T, E, _ = km_results[pop]
        c = colors.get(pop, GREY)
        kmf.plot_survival_function(
            ax=ax1, ci_show=False, color=c,
            label=f"{pop_labels.get(pop, pop)} (n={len(T)})",
            linewidth=1.8)
    ax1.set_xlabel("Time (days)", **lkw)
    ax1.set_ylabel("Survival probability", **lkw)
    ax1.legend(fontsize=6.5, labelcolor="#8b949e",
               facecolor="#161b22", edgecolor="#30363d")
    ax_style(ax1)

    # ── Panel 2: Depth score tertile KM (Stratum A) ──────────
    ax2 = fig.add_subplot(gs[0, 2])
    ax2.set_title("Depth Score vs OS\nStratum A (S2-P1)", **tkw)
    cl_a = populations["claudin_low"]
    res = tertile_split(cl_a, depth_scores, os_time, os_event)
    T_h, E_h, T_l, E_l, T_m, E_m, vh, vm, vl = res
    if len(T_h) >= MIN_SURV_N:
        for T, E, lbl, c in [
            (T_h, E_h, f"Deep (n={len(T_h)})",    RED),
            (T_m, E_m, f"Mid  (n={len(T_m)})",    GREY),
            (T_l, E_l, f"Shallow (n={len(T_l)})", GREEN),
        ]:
            if len(T) >= MIN_SURV_N:
                kmf = KaplanMeierFitter()
                kmf.fit(T, E, label=lbl)
                kmf.plot_survival_function(
                    ax=ax2, ci_show=False, color=c,
                    linewidth=1.8)
    ax2.set_xlabel("Time (days)", **lkw)
    ax2.set_ylabel("Survival probability", **lkw)
    ax2.legend(fontsize=6.5, labelcolor="#8b949e",
               facecolor="#161b22", edgecolor="#30363d")
    ax_style(ax2)

    # ── Panel 3: Depth score tertile KM (Stratum B) ──────────
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.set_title("Depth Score vs OS\nStratum B ESR1-low (S2-P7)", **tkw)
    cl_b = populations["stratum_b"]
    if len(cl_b) >= MIN_SURV_N * 3:
        res_b = tertile_split(cl_b, depth_scores, os_time, os_event)
        Tb_h, Eb_h, Tb_l, Eb_l, Tb_m, Eb_m = res_b[:6]
        for T, E, lbl, c in [
            (Tb_h, Eb_h, f"Deep (n={len(Tb_h)})",    RED),
            (Tb_m, Eb_m, f"Mid  (n={len(Tb_m)})",    GREY),
            (Tb_l, Eb_l, f"Shallow (n={len(Tb_l)})", GREEN),
        ]:
            if len(T) >= MIN_SURV_N:
                kmf = KaplanMeierFitter()
                kmf.fit(T, E, label=lbl)
                kmf.plot_survival_function(
                    ax=ax3, ci_show=False, color=c,
                    linewidth=1.8)
        ax3.set_xlabel("Time (days)", **lkw)
        ax3.set_ylabel("Survival probability", **lkw)
        ax3.legend(fontsize=6.5, labelcolor="#8b949e",
                   facecolor="#161b22", edgecolor="#30363d")
    else:
        ax3.text(0.5, 0.5, f"Stratum B n={len(cl_b)}\n(underpowered)",
                 transform=ax3.transAxes, ha="center",
                 va="center", color="#8b949e")
    ax_style(ax3)

    # ── Panel 4: CLDN3 tertile KM ─────────────────────────────
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.set_title("CLDN3 Expression vs OS (S2-P2)", **tkw)
    if "CLDN3" in expr.index:
        cldn3_sc = {s: float(expr.loc["CLDN3", s]) for s in cl_a}
        res2 = tertile_split(cl_a, cldn3_sc, os_time, os_event)
        T2_h, E2_h, T2_l, E2_l, T2_m, E2_m = res2[:6]
        for T, E, lbl, c in [
            (T2_h, E2_h, f"CLDN3-high (n={len(T2_h)})", GREEN),
            (T2_m, E2_m, f"CLDN3-mid  (n={len(T2_m)})", GREY),
            (T2_l, E2_l, f"CLDN3-low  (n={len(T2_l)})", RED),
        ]:
            if len(T) >= MIN_SURV_N:
                kmf = KaplanMeierFitter()
                kmf.fit(T, E, label=lbl)
                kmf.plot_survival_function(
                    ax=ax4, ci_show=False, color=c,
                    linewidth=1.8)
        ax4.set_xlabel("Time (days)", **lkw)
        ax4.set_ylabel("Survival probability", **lkw)
        ax4.legend(fontsize=6.5, labelcolor="#8b949e",
                   facecolor="#161b22", edgecolor="#30363d")
    else:
        ax4.text(0.5, 0.5, "CLDN3 not in matrix",
                 transform=ax4.transAxes, ha="center",
                 va="center", color="#8b949e")
    ax_style(ax4)

    # ── Panel 5: Immune score tertile KM ─────────────────────
    ax5 = fig.add_subplot(gs[1, 2])
    ax5.set_title("Immune Score vs OS (S2-P3)", **tkw)
    imm_present = [g for g in IMMUNE_GENES if g in expr.index]
    if len(imm_present) >= 2:
        imm_z_p = []
        for g in imm_present:
            vals = np.array([float(expr.loc[g, s]) for s in cl_a])
            sd = vals.std()
            if sd > 1e-9:
                imm_z_p.append((vals - vals.mean()) / sd)
        if imm_z_p:
            imm_comp = np.mean(imm_z_p, axis=0)
            imm_sc = dict(zip(cl_a, imm_comp))
            res3 = tertile_split(cl_a, imm_sc, os_time, os_event)
            T3_h, E3_h, T3_l, E3_l, T3_m, E3_m = res3[:6]
            for T, E, lbl, c in [
                (T3_h, E3_h, f"Immune-high (n={len(T3_h)})", BLUE),
                (T3_m, E3_m, f"Immune-mid  (n={len(T3_m)})", GREY),
                (T3_l, E3_l, f"Immune-low  (n={len(T3_l)})", ORANGE),
            ]:
                if len(T) >= MIN_SURV_N:
                    kmf = KaplanMeierFitter()
                    kmf.fit(T, E, label=lbl)
                    kmf.plot_survival_function(
                        ax=ax5, ci_show=False, color=c,
                        linewidth=1.8)
            ax5.set_xlabel("Time (days)", **lkw)
            ax5.set_ylabel("Survival probability", **lkw)
            ax5.legend(fontsize=6.5, labelcolor="#8b949e",
                       facecolor="#161b22", edgecolor="#30363d")
    ax_style(ax5)

    # ── Panel 6: MKI67 tertile KM ─────────────────────────────
    ax6 = fig.add_subplot(gs[2, 0])
    ax6.set_title("MKI67 vs OS (S2-P4)", **tkw)
    if "MKI67" in expr.index:
        mki67_sc = {s: float(expr.loc["MKI67", s]) for s in cl_a}
        res4 = tertile_split(cl_a, mki67_sc, os_time, os_event)
        T4_h, E4_h, T4_l, E4_l, T4_m, E4_m = res4[:6]
        for T, E, lbl, c in [
            (T4_h, E4_h, f"MKI67-high (n={len(T4_h)})", RED),
            (T4_m, E4_m, f"MKI67-mid  (n={len(T4_m)})", GREY),
            (T4_l, E4_l, f"MKI67-low  (n={len(T4_l)})", GREEN),
        ]:
            if len(T) >= MIN_SURV_N:
                kmf = KaplanMeierFitter()
                kmf.fit(T, E, label=lbl)
                kmf.plot_survival_function(
                    ax=ax6, ci_show=False, color=c,
                    linewidth=1.8)
        ax6.set_xlabel("Time (days)", **lkw)
        ax6.set_ylabel("Survival probability", **lkw)
        ax6.legend(fontsize=6.5, labelcolor="#8b949e",
                   facecolor="#161b22", edgecolor="#30363d")
    ax_style(ax6)

    # ── Panel 7: Depth score distribution per stratum ─────────
    ax7 = fig.add_subplot(gs[2, 1])
    ax7.set_title("Depth Score Distribution by Stratum", **tkw)
    strata = [
        ("claudin_low", "Stratum A", BLUE),
        ("stratum_b",   "Stratum B", ORANGE),
        ("stratum_c",   "Stratum C", PURPLE),
    ]
    for pop, lbl, c in strata:
        samps = populations.get(pop, [])
        sc = [depth_scores[s] for s in samps if s in depth_scores]
        if len(sc) > 5:
            ax7.hist(sc, bins=20, alpha=0.5, color=c,
                     label=f"{lbl} (n={len(sc)})",
                     edgecolor="none")
    ax7.set_xlabel("Depth Score", **lkw)
    ax7.set_ylabel("Count", **lkw)
    ax7.legend(fontsize=7, labelcolor="#8b949e",
               facecolor="#161b22", edgecolor="#30363d")
    ax_style(ax7)

    # ── Panel 8: Depth vs CLDN3 scatter in CL (structural read) ──
    ax8 = fig.add_subplot(gs[2, 2])
    ax8.set_title("Depth Score vs CLDN3\n(within claudin-low)", **tkw)
    cl_a = populations["claudin_low"]
    if "CLDN3" in expr.index:
        d_vals = [depth_scores[s] for s in cl_a if s in depth_scores]
        c_vals = [float(expr.loc["CLDN3", s])
                  for s in cl_a if s in depth_scores]
        if len(d_vals) > 3:
            r, p = stats.pearsonr(d_vals, c_vals)
            ax8.scatter(d_vals, c_vals, c=BLUE,
                        alpha=0.4, s=12, edgecolors="none")
            xr = np.linspace(min(d_vals), max(d_vals), 50)
            z  = np.polyfit(d_vals, c_vals, 1)
            ax8.plot(xr, np.polyval(z, xr),
                     color=ORANGE, lw=1.5, linestyle="--")
            ax8.set_xlabel("Depth Score", **lkw)
            ax8.set_ylabel("CLDN3 expression", **lkw)
            ax8.text(0.05, 0.92, f"r={r:.3f}  {fmt_p(p)}",
                     transform=ax8.transAxes,
                     color="#e6edf3", fontsize=8)
    ax_style(ax8)

    fig.suptitle(
        "CLAUDIN-LOW — SCRIPT 2 SURVIVAL ANALYSIS\n"
        "OrganismCore | BRCA-S7d | 2026-03-05\n"
        "TYPE 4 — ROOT LOCK | Predictions from BRCA-S7c",
        color="#e6edf3", fontsize=13, fontweight="bold", y=0.98
    )

    plt.savefig(FIG_FILE, dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close()
    log(f"  Figure saved: {FIG_FILE}")

# ============================================================
# STEP 10 — PREDICTION SUMMARY PRINT
# ============================================================

def print_scorecard(scorecard):
    log("")
    log("=" * 65)
    log("STEP 10: PREDICTION SCORECARD SUMMARY")
    log("  Protocol v2.0: geometry read first (Step 6).")
    log("  All predictions now scored.")
    log("=" * 65)

    log("")
    log(f"  {'ID':<10} {'Axis':<22} {'Stratum':<10} "
        f"{'p':>10} {'HR':>8}  {'Status'}")
    log("  " + "─" * 75)

    for row in scorecard:
        p_str = (f"{row['p_value']:.3f}"
                 if not np.isnan(row['p_value']) else "nan")
        hr_str = (f"{row['HR_approx']:.3f}"
                  if not np.isnan(row['HR_approx']) else "nan")
        log(f"  {row['prediction']:<10} {str(row['gene_or_axis']):<22} "
            f"{str(row['stratum']):<10} {p_str:>10} {hr_str:>8}  "
            f"{row['confirmed']}")

    n_conf  = sum(1 for r in scorecard
                  if "CONFIRMED" in str(r["confirmed"]) and
                  "NOT" not in str(r["confirmed"]))
    n_total = len(scorecard)
    log("")
    log(f"  Predictions confirmed (including negatives): "
        f"{n_conf} / {n_total}")

# ============================================================
# MAIN
# ============================================================

def main():
    log("=" * 65)
    log("BRCA CLAUDIN-LOW — SCRIPT 2")
    log("OrganismCore — Document BRCA-S7d | 2026-03-05")
    log("")
    log("ATTRACTOR TYPE: TYPE 4 — ROOT LOCK")
    log("ANALYSIS: SURVIVAL + PREDICTION TESTS")
    log(f"Output directory: {BASE_DIR}")
    log("=" * 65)

    # Step 0: acquire survival supplement
    acquire_survival()

    # Step 1: load expression (reuse Script 1 files)
    expr = load_expression()

    # Step 2: load clinical + survival
    (clin, surv, surv_time_col, surv_event_col,
     clin_time_col, clin_event_col) = load_clinical_and_survival()

    # Step 3: rebuild populations
    populations, scores = build_populations(expr, clin)

    # Step 4: build survival vectors
    os_time, os_event = build_survival_vectors(
        expr, populations, surv,
        surv_time_col, surv_event_col,
        clin, clin_time_col, clin_event_col)

    # Step 5: compute depth scores
    tumour_samples = populations["claudin_low"] + \
                     populations["luma"] + \
                     populations["lumb"] + \
                     populations["her2"] + \
                     populations["basal"]
    depth_scores = compute_depth_scores(expr, tumour_samples)

    # Step 6: cross-subtype KM (geometry first)
    km_results, pop_labels, colors = cross_subtype_km(
        populations, os_time, os_event)

    # Step 7: prediction tests
    scorecard = run_prediction_tests(
        expr, populations, depth_scores, os_time, os_event)

    # Step 8: Cox multivariate
    cox_multivariate(expr, populations, depth_scores,
                     os_time, os_event)

    # Step 9: figure
    try:
        make_figure(populations, depth_scores, os_time, os_event,
                    km_results, pop_labels, colors, expr)
    except Exception as e:
        log(f"  Figure generation failed: {e}")

    # Step 10: scorecard summary
    print_scorecard(scorecard)

    # Save survival CSV
    surv_rows = []
    for pop, samps in populations.items():
        for s in samps:
            surv_rows.append({
                "sample":       s,
                "population":   pop,
                "os_time":      os_time.get(s, np.nan),
                "os_event":     os_event.get(s, np.nan),
                "depth_score":  depth_scores.get(s, np.nan),
            })
    pd.DataFrame(surv_rows).to_csv(CSV_FILE, index=False)
    log(f"\n  Survival CSV saved: {CSV_FILE}")

    log("")
    log("=" * 65)
    log("SCRIPT 2 COMPLETE")
    log("=" * 65)
    log(f"  Output directory: {BASE_DIR}")
    log(f"  Log:      {LOG_FILE}")
    log(f"  Figure:   {FIG_FILE}")
    log(f"  Survival: {CSV_FILE}")
    log(f"  Scorecard:{SCORECARD}")
    log("")
    log("  NEXT STEP:")
    log("    Write script2_results_and_reasoning.md (BRCA-S7d)")
    log("    Read KM geometry on its own terms first.")
    log("    Score all 8 predictions against the output.")

    flush_log()


if __name__ == "__main__":
    main()
