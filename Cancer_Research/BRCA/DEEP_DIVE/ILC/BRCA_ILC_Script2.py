"""
BRCA ILC — SCRIPT 2 (v1)
OrganismCore — Document BRCA-S6c/d | 2026-03-05

PURPOSE:
  Survival analysis and treatment-context validation within ILC.
  Uses TCGA-BRCA clinical + expression data already downloaded by Script 1.
  No new large downloads required (only survival supplement if not present).

ANALYSES:
  1. OS survival: ESR1-high vs ESR1-low within ILC           (S2-P1)
  2. OS survival: CDH1-depth (low vs high) within ILC        (S2-P2)
  3. OS survival: ILC vs LumA                                (S2-P3)
  4. OS survival: SPDEF-high vs SPDEF-low within ILC         (S2-P4)
  5. OS survival: PTEN-low vs PTEN-high within ILC           (S2-P5)
  6. OS survival: EZH2-high vs EZH2-low within ILC           (S2-P6)
  7. OS survival: CCND1-high vs CCND1-low within ILC         (S2-P7)
  8. OS survival: MKI67-high vs MKI67-low within ILC         (S2-P8)
  9. OS survival: all subtypes compared (ILC/LumA/TNBC/HER2/LumB)
 10. BONUS: PIK3CA mRNA tertile survival within ILC

DATA REUSE:
  Expression and clinical files from ILC_s1_analysis/data/
  Survival supplement downloaded separately if not cached.

SELF-CONTAINED:
  All paths derived from SCRIPT_DIR.
  Results go to ILC_s2_results/ (separate from Script 1 output).
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

# ============================================================
# CONFIGURATION
# ============================================================

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))

# Reuse data already downloaded by Script 1
S1_DATA_DIR  = os.path.join(SCRIPT_DIR, "ILC_s1_analysis", "data")

# Script 2 output directory
BASE_DIR     = os.path.join(SCRIPT_DIR, "ILC_s2_results")
os.makedirs(BASE_DIR, exist_ok=True)

LOG_FILE     = os.path.join(BASE_DIR, "ilc_s2_log.txt")
FIG_FILE     = os.path.join(BASE_DIR, "ilc_s2_figure.png")
CSV_FILE     = os.path.join(BASE_DIR, "ilc_s2_survival.csv")

# ── Input files from Script 1 ────────────────────────────────
EXPR_FILE    = os.path.join(S1_DATA_DIR, "TCGA_BRCA_HiSeqV2.gz")
CLIN_FILE    = os.path.join(S1_DATA_DIR, "TCGA_BRCA_clinicalMatrix.tsv")

# ── Survival supplement (Xena pan-can) ───────────────────────
SURV_URLS = [
    "https://pancanatlas.xenahubs.net/download/Survival_SupplementalTable_S1_20171025_xena_sp",
    "https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/"
    "Survival_SupplementalTable_S1_20171025_xena_sp",
    "https://gdc.xenahubs.net/download/TCGA-CDR-SupplementalTableS1.tsv.gz",
]
SURV_FILE    = os.path.join(S1_DATA_DIR, "TCGA_pancan_survival.tsv")

# ── Cohort thresholds ────────────────────────────────────────
MIN_ILC      = 30
MIN_NORMAL   = 15
PATIENT_PREFIX_LEN = 12
SAMPLE_TYPE_POS    = (13, 15)

ILC_KEYWORDS  = ["infiltrating lobular", "invasive lobular",
                 "lobular carcinoma", "lobular"]
LUMA_LABELS   = {"LumA"}
TNBC_LABELS   = {"Basal"}
HER2_LABELS   = {"Her2"}
LUMB_LABELS   = {"LumB"}

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

def download_file(url, dest, label="", timeout=120, chunk=1024 * 1024):
    try:
        log(f"  Downloading {label or os.path.basename(url[:60])} ...")
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            total = int(resp.headers.get("Content-Length", 0))
            done  = 0
            with open(dest, "wb") as fh:
                while True:
                    buf = resp.read(chunk)
                    if not buf:
                        break
                    fh.write(buf)
                    done += len(buf)
                    if total > 0:
                        print(f"    {done*100//total}%...", end="\r", flush=True)
        print()
        log(f"  OK: {os.path.basename(dest)} ({os.path.getsize(dest)/1e6:.1f} MB)")
        return True
    except Exception as e:
        log(f"  Failed: {type(e).__name__}: {e}")
        if os.path.exists(dest):
            os.remove(dest)
        return False


def try_urls(urls, dest, label=""):
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        log(f"  Already present: {os.path.basename(dest)} "
            f"({os.path.getsize(dest)/1e6:.1f} MB)")
        return True
    for i, url in enumerate(urls, 1):
        log(f"  Attempt {i}/{len(urls)}: {url[:72]}...")
        if download_file(url, dest, label):
            return True
        time.sleep(2)
    return False


def open_gz(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

# ============================================================
# LOG-RANK TEST (self-contained — no lifelines dependency)
# ============================================================

def logrank_test(t1, e1, t2, e2):
    """
    Manual log-rank test between two groups.
    t = time array, e = event array (1=event, 0=censored)
    Returns: (test_statistic, p_value)
    """
    t1, e1 = np.array(t1, float), np.array(e1, float)
    t2, e2 = np.array(t2, float), np.array(e2, float)

    # Pool all unique event times
    all_t = np.unique(np.concatenate([t1[e1 == 1], t2[e2 == 1]]))
    if len(all_t) == 0:
        return 0.0, 1.0

    O1, E1 = 0.0, 0.0
    V = 0.0

    for t in all_t:
        n1 = np.sum(t1 >= t)
        n2 = np.sum(t2 >= t)
        n  = n1 + n2
        if n == 0:
            continue
        o1 = np.sum((t1 == t) & (e1 == 1))
        o2 = np.sum((t2 == t) & (e2 == 1))
        o  = o1 + o2
        e1_exp = o * n1 / n
        O1 += o1
        E1 += e1_exp
        if n > 1:
            V += o * n1 * n2 * (n - o) / (n ** 2 * (n - 1))

    if V <= 0:
        return 0.0, 1.0

    chi2 = (O1 - E1) ** 2 / V
    p = 1 - stats.chi2.cdf(chi2, df=1)
    return chi2, p


def kaplan_meier(times, events):
    """
    Compute KM survival curve.
    Returns (time_points, survival_probabilities)
    """
    times  = np.array(times, float)
    events = np.array(events, float)
    order  = np.argsort(times)
    times  = times[order]
    events = events[order]

    t_km = [0.0]
    s_km = [1.0]
    n    = len(times)
    S    = 1.0

    for i, t in enumerate(np.unique(times)):
        at_risk = np.sum(times >= t)
        died    = np.sum((times == t) & (events == 1))
        if died == 0:
            continue
        S = S * (1 - died / at_risk)
        t_km.append(t)
        s_km.append(S)

    return np.array(t_km), np.array(s_km)


def median_os(times, events):
    t, s = kaplan_meier(times, events)
    idx  = np.where(s <= 0.5)[0]
    return float(t[idx[0]]) if len(idx) > 0 else float("inf")


def hr_estimate(t1, e1, t2, e2):
    """
    Crude HR estimate: (O1/E1) / (O2/E2) from log-rank observed/expected.
    """
    t1, e1 = np.array(t1, float), np.array(e1, float)
    t2, e2 = np.array(t2, float), np.array(e2, float)
    all_t  = np.unique(np.concatenate([t1[e1 == 1], t2[e2 == 1]]))

    O1, E1, O2, E2 = 0.0, 0.0, 0.0, 0.0
    for t in all_t:
        n1 = np.sum(t1 >= t)
        n2 = np.sum(t2 >= t)
        n  = n1 + n2
        if n == 0:
            continue
        o1 = np.sum((t1 == t) & (e1 == 1))
        o2 = np.sum((t2 == t) & (e2 == 1))
        o  = o1 + o2
        O1 += o1; O2 += o2
        E1 += o * n1 / n if n > 0 else 0
        E2 += o * n2 / n if n > 0 else 0

    hr = (O1 / E1) / (O2 / E2) if E1 > 0 and E2 > 0 and O2 > 0 and E2 > 0 else float("nan")
    return hr

# ============================================================
# LOAD DATA (reuse from Script 1 downloads)
# ============================================================

def load_expression():
    log("")
    log("=" * 65)
    log("LOADING EXPRESSION MATRIX (from Script 1 data)")
    log("=" * 65)

    if not os.path.exists(EXPR_FILE):
        log(f"  FATAL: Expression file not found: {EXPR_FILE}")
        log("  Run Script 1 first to download data.")
        flush_log(); sys.exit(1)

    with open_gz(EXPR_FILE, "rt") as fh:
        first = fh.readline()
    sep  = "\t" if "\t" in first else ","
    expr = pd.read_csv(EXPR_FILE, sep=sep, index_col=0, low_memory=False)

    if expr.shape[0] < 5000 and expr.shape[1] > 10000:
        expr = expr.T

    expr = expr.apply(pd.to_numeric, errors="coerce")
    expr.dropna(how="all", inplace=True)
    log(f"  Expression shape: {expr.shape}")
    log(f"  Sample example: {str(expr.columns[0])}")
    return expr


def load_clinical():
    log("")
    log("=" * 65)
    log("LOADING CLINICAL MATRIX")
    log("=" * 65)

    if not os.path.exists(CLIN_FILE):
        log(f"  FATAL: Clinical file not found: {CLIN_FILE}")
        flush_log(); sys.exit(1)

    with open_gz(CLIN_FILE, "rt") as fh:
        first = fh.readline()
    sep  = "\t" if first.count("\t") > first.count(",") else ","
    clin = pd.read_csv(CLIN_FILE, sep=sep, low_memory=False)

    # Identify barcode column
    barcode_col = None
    if clin.iloc[:, 0].astype(str).str.startswith("TCGA").any():
        barcode_col = clin.columns[0]
    else:
        for col in ["sampleID", "submitter_id", "bcr_patient_barcode",
                    "sample", "barcode"]:
            if col in clin.columns and clin[col].astype(str).str.startswith("TCGA").any():
                barcode_col = col
                break

    if barcode_col:
        clin = clin.set_index(barcode_col)
    clin.index = clin.index.astype(str).str.strip()

    # Identify columns
    hist_col = None
    for col in ["histological_type", "primary_diagnosis",
                "histologic_diagnosis", "icd_o_3_histology"]:
        if col in clin.columns:
            hist_col = col
            break

    pam50_col = None
    for col in ["PAM50Call_RNAseq", "PAM50_mRNA_nature2012",
                "Integrated_Clusters_with_PAM50__nature2012"]:
        if col in clin.columns:
            pam50_col = col
            break

    log(f"  Clinical shape: {clin.shape}")
    log(f"  Histology column: '{hist_col}'")
    log(f"  PAM50 column: '{pam50_col}'")
    return clin, hist_col, pam50_col


def load_survival():
    log("")
    log("=" * 65)
    log("LOADING SURVIVAL DATA")
    log("=" * 65)

    if not try_urls(SURV_URLS, SURV_FILE, "TCGA pancan survival"):
        log("  WARNING: Survival supplement unavailable.")
        log("  Falling back to OS fields in clinical matrix.")
        return None

    with open_gz(SURV_FILE, "rt") as fh:
        first = fh.readline()
    sep  = "\t" if "\t" in first else ","
    surv = pd.read_csv(SURV_FILE, sep=sep, low_memory=False)
    log(f"  Survival shape: {surv.shape}")
    log(f"  Columns: {list(surv.columns[:12])}")

    # Try to find sample/barcode column
    for col in surv.columns:
        if surv[col].astype(str).str.startswith("TCGA").any():
            surv = surv.set_index(col)
            log(f"  Survival index set to: '{col}'")
            break
    surv.index = surv.index.astype(str).str.strip()
    return surv

# ============================================================
# CLASSIFY SAMPLES (same logic as Script 1)
# ============================================================

def build_patient_map(expr_cols):
    patient_map = {}
    all_tumour  = []
    all_normal  = []
    for col in expr_cols:
        s   = str(col).strip()
        pfx = s[:PATIENT_PREFIX_LEN] if len(s) >= PATIENT_PREFIX_LEN else s
        patient_map.setdefault(pfx, []).append(s)
        st  = s[SAMPLE_TYPE_POS[0]:SAMPLE_TYPE_POS[1]] if len(s) >= SAMPLE_TYPE_POS[1] else "??"
        if st == "01": all_tumour.append(s)
        elif st == "11": all_normal.append(s)
    return patient_map, all_tumour, all_normal


def classify_samples(expr, clin, hist_col, pam50_col):
    log("")
    log("=" * 65)
    log("CLASSIFYING SAMPLES")
    log("=" * 65)

    expr_cols   = [str(c).strip() for c in expr.columns]
    patient_map, all_tumour, all_normal = build_patient_map(expr_cols)

    ilc_s = []; luma_s = []; tnbc_s = []; her2_s = []; lumb_s = []
    norm_s = list(all_normal)
    clin.index = clin.index.astype(str).str.strip()

    for clin_id in clin.index:
        pfx     = str(clin_id)[:PATIENT_PREFIX_LEN]
        e_samps = patient_map.get(pfx, [])
        if not e_samps:
            continue
        try:
            row = clin.loc[clin_id]
        except KeyError:
            continue
        hist_val  = str(row[hist_col]).lower()  if hist_col  and hist_col  in clin.columns else ""
        pam50_val = str(row[pam50_col]).strip() if pam50_col and pam50_col in clin.columns else ""
        is_ilc    = any(kw in hist_val for kw in ILC_KEYWORDS)

        for s in e_samps:
            st = s[SAMPLE_TYPE_POS[0]:SAMPLE_TYPE_POS[1]] if len(s) >= SAMPLE_TYPE_POS[1] else "??"
            if st == "11":
                if s not in norm_s: norm_s.append(s)
            elif st == "01":
                if is_ilc:            ilc_s.append(s)
                elif pam50_val in LUMA_LABELS:  luma_s.append(s)
                elif pam50_val in TNBC_LABELS:  tnbc_s.append(s)
                elif pam50_val in HER2_LABELS:  her2_s.append(s)
                elif pam50_val in LUMB_LABELS:  lumb_s.append(s)

    log(f"  ILC: {len(ilc_s)}  LumA: {len(luma_s)}  TNBC: {len(tnbc_s)}  "
        f"HER2: {len(her2_s)}  LumB: {len(lumb_s)}  Normal: {len(norm_s)}")

    if len(ilc_s) < MIN_ILC:
        log(f"  FATAL: Only {len(ilc_s)} ILC samples. Run Script 1 first.")
        flush_log(); sys.exit(1)

    return {"ilc": ilc_s, "luma": luma_s, "tnbc": tnbc_s,
            "her2": her2_s, "lumb": lumb_s, "normal": norm_s}

# ============================================================
# BUILD SURVIVAL VECTORS
# Maps sample barcodes → (OS_time, OS_event)
# Tries: pancan supplement first, then clinical matrix fallback
# ============================================================

def build_survival_vectors(pops, clin, surv):
    log("")
    log("=" * 65)
    log("BUILDING SURVIVAL VECTORS")
    log("=" * 65)

    # ── Try pancan supplement ────────────────────────────────
    os_time_col  = None
    os_event_col = None

    if surv is not None:
        time_candidates  = ["OS.time", "OS_time", "overall_survival",
                            "days_to_death", "OS.time.cr"]
        event_candidates = ["OS", "OS_event", "vital_status",
                           "overall_survival_status", "OS.cr"]
        for col in time_candidates:
            if col in surv.columns:
                os_time_col = col
                break
        for col in event_candidates:
            if col in surv.columns:
                os_event_col = col
                break
        log(f"  Pancan survival: time='{os_time_col}', event='{os_event_col}'")

    # ── Clinical fallback columns ────────────────────────────
    clin_time_cols  = ["OS_Time_nature2012", "days_to_last_followup",
                       "days_to_death", "days_to_last_known_alive"]
    clin_event_cols = ["OS_event_nature2012", "vital_status",
                       "Vital_Status_nature2012"]
    clin_time_col   = next((c for c in clin_time_cols  if c in clin.columns), None)
    clin_event_col  = next((c for c in clin_event_cols if c in clin.columns), None)
    log(f"  Clinical fallback: time='{clin_time_col}', event='{clin_event_col}'")

    def get_os(sample_barcode):
        """Return (time_days, event_01) or (nan, nan) for a sample barcode."""
        # sample barcode e.g. TCGA-BH-A0BZ-01
        # pancan supplement indexed by TCGA-BH-A0BZ (12 chars) or sample barcode
        pfx15 = sample_barcode[:15]   # TCGA-XX-XXXX-01 (15 chars as in Xena)
        pfx12 = sample_barcode[:12]   # TCGA-XX-XXXX (patient)

        # Try pancan supplement
        if surv is not None and os_time_col and os_event_col:
            for key in [sample_barcode, pfx15, pfx12]:
                if key in surv.index:
                    t_raw = surv.loc[key, os_time_col]
                    e_raw = surv.loc[key, os_event_col]
                    t = pd.to_numeric(t_raw, errors="coerce")
                    # event: 1 = dead, 0 = alive/censored
                    if isinstance(e_raw, str):
                        e = 1.0 if e_raw.lower() in ["1", "dead", "deceased",
                                                       "1:deceased"] else 0.0
                    else:
                        e = float(pd.to_numeric(e_raw, errors="coerce") or 0)
                    if not np.isnan(t) and t > 0:
                        return t, e

        # Try clinical matrix
        if clin_time_col and clin_event_col:
            for key in [sample_barcode, pfx15, pfx12]:
                if key in clin.index:
                    t_raw = clin.loc[key, clin_time_col]
                    e_raw = clin.loc[key, clin_event_col]
                    t = pd.to_numeric(t_raw, errors="coerce")
                    if isinstance(e_raw, str):
                        e = 1.0 if e_raw.lower() in ["1", "dead", "deceased",
                                                       "dead with tumor"] else 0.0
                    else:
                        e = float(pd.to_numeric(e_raw, errors="coerce") or 0)
                    if not np.isnan(t) and t > 0:
                        return t, e

        return float("nan"), float("nan")

    # Build dict: sample_barcode → (time, event)
    surv_map = {}
    for grp, samps in pops.items():
        if grp == "normal":
            continue
        for s in samps:
            t, e = get_os(s)
            if not np.isnan(t):
                surv_map[s] = (t, e)

    total_with_surv = len(surv_map)
    ilc_with_surv   = sum(1 for s in pops["ilc"] if s in surv_map)
    log(f"  Total samples with survival data: {total_with_surv}")
    log(f"  ILC samples with survival data:   {ilc_with_surv}")

    if ilc_with_surv < 20:
        log("  WARNING: Fewer than 20 ILC samples with survival data.")
        log("  Survival analyses will have low power.")

    return surv_map


def get_group_survival(samps, surv_map):
    """Return (times, events) arrays for a list of sample barcodes."""
    times, events = [], []
    for s in samps:
        if s in surv_map:
            t, e = surv_map[s]
            times.append(t)
            events.append(e)
    return np.array(times), np.array(events)

# ============================================================
# TERTILE SPLIT UTILITY
# ============================================================

def tertile_split(expr, gene, sample_list):
    """
    Split sample_list into high (top tertile) and low (bottom tertile)
    by expression of gene. Returns (high_samples, low_samples).
    """
    if gene not in expr.index:
        return [], []
    vals   = expr.loc[gene, sample_list].values.astype(float)
    t33    = np.nanpercentile(vals, 33)
    t67    = np.nanpercentile(vals, 67)
    high_s = [s for s, v in zip(sample_list, vals) if v >= t67]
    low_s  = [s for s, v in zip(sample_list, vals) if v <= t33]
    return high_s, low_s


def median_split(expr, gene, sample_list):
    """Split at median. Returns (high_samples, low_samples)."""
    if gene not in expr.index:
        return [], []
    vals  = expr.loc[gene, sample_list].values.astype(float)
    med   = np.nanmedian(vals)
    high_s = [s for s, v in zip(sample_list, vals) if v > med]
    low_s  = [s for s, v in zip(sample_list, vals) if v <= med]
    return high_s, low_s

# ============================================================
# SURVIVAL ANALYSIS — WITHIN ILC
# ============================================================

def fmt_p(p):
    if p < 1e-10: return f"p={p:.2e} ***"
    if p < 1e-5:  return f"p={p:.2e} **"
    if p < 0.05:  return f"p={p:.4f} *"
    if p < 0.15:  return f"p={p:.4f} (trend)"
    return              f"p={p:.4f} (ns)"


def survival_test(label, pred_id, high_s, low_s, surv_map,
                  high_label="High", low_label="Low",
                  predicted_direction="high_better"):
    """
    Run log-rank test and report result.
    Returns dict with result metadata.
    """
    t_hi, e_hi = get_group_survival(high_s, surv_map)
    t_lo, e_lo = get_group_survival(low_s,  surv_map)

    n_hi = len(t_hi)
    n_lo = len(t_lo)
    ev_hi = int(e_hi.sum()) if len(e_hi) > 0 else 0
    ev_lo = int(e_lo.sum()) if len(e_lo) > 0 else 0

    if n_hi < 5 or n_lo < 5:
        log(f"  {pred_id} [{label}]: Insufficient n (high={n_hi}, low={n_lo}). Skipped.")
        return {"pred_id": pred_id, "label": label, "status": "INSUFFICIENT_N",
                "n_high": n_hi, "n_low": n_lo, "p": float("nan"),
                "hr": float("nan")}

    chi2, p = logrank_test(t_hi, e_hi, t_lo, e_lo)
    hr      = hr_estimate(t_hi, e_hi, t_lo, e_lo)
    med_hi  = median_os(t_hi, e_hi)
    med_lo  = median_os(t_lo, e_lo)

    # Determine if direction matches prediction
    if predicted_direction == "high_better":
        direction_correct = hr < 1.0
        direction_str = "HIGH better (predicted)" if hr < 1.0 else "LOW better (opposite)"
    elif predicted_direction == "low_better":
        direction_correct = hr > 1.0
        direction_str = "LOW better (predicted)" if hr > 1.0 else "HIGH better (opposite)"
    else:  # null prediction
        direction_correct = p > 0.10
        direction_str = "NULL confirmed" if p > 0.10 else f"SIGNIFICANT (unexpected)"

    if p < 0.05 and direction_correct:
        status = "CONFIRMED"
    elif p < 0.15 and direction_correct:
        status = "TREND"
    elif p >= 0.10 and predicted_direction == "null":
        status = "NULL CONFIRMED"
    elif p < 0.05 and not direction_correct:
        status = "OPPOSITE"
    else:
        status = "NOT CONFIRMED"

    log(f"  {pred_id} [{label}]")
    log(f"    {high_label}: n={n_hi}, events={ev_hi}, median OS={med_hi:.0f}d")
    log(f"    {low_label}:  n={n_lo}, events={ev_lo}, median OS={med_lo:.0f}d")
    log(f"    Log-rank: {fmt_p(p)}  HR={hr:.3f}  [{direction_str}]")
    log(f"    Status: [{status}]")
    log("")

    return {"pred_id": pred_id, "label": label, "status": status,
            "n_high": n_hi, "n_low": n_lo,
            "events_high": ev_hi, "events_low": ev_lo,
            "p": p, "hr": hr,
            "med_os_high": med_hi, "med_os_low": med_lo,
            "high_label": high_label, "low_label": low_label,
            "t_high": t_hi, "e_high": e_hi,
            "t_low": t_lo, "e_low": e_lo}


def survival_analyses(expr, pops, surv_map):
    log("")
    log("=" * 65)
    log("SURVIVAL ANALYSES")
    log("=" * 65)

    ilc_s  = pops["ilc"]
    results = []

    # ── S2-P1: ESR1 high vs low within ILC ───────────────────
    log("")
    log("  S2-P1 — ESR1 PREDICTS OS WITHIN ILC")
    log("  Prediction: ESR1-high = better survival")
    log("-" * 65)
    hi, lo = tertile_split(expr, "ESR1", ilc_s)
    res = survival_test("ESR1 tertile within ILC", "S2-P1",
                        hi, lo, surv_map,
                        high_label="ESR1-high", low_label="ESR1-low",
                        predicted_direction="high_better")
    results.append(res)

    # ── S2-P2: CDH1 depth (low CDH1 = deep = worse) ──────────
    log("  S2-P2 — CDH1 DEPTH PREDICTS OS WITHIN ILC")
    log("  Prediction: CDH1-low (deep) = worse survival")
    log("-" * 65)
    # For CDH1: low CDH1 = deeper. Split: high CDH1 = shallow, low CDH1 = deep.
    hi_cdh1, lo_cdh1 = tertile_split(expr, "CDH1", ilc_s)
    # Test: CDH1-low (deep) vs CDH1-high (shallow) — CDH1-low predicted worse
    res = survival_test("CDH1 tertile within ILC", "S2-P2",
                        lo_cdh1, hi_cdh1, surv_map,
                        high_label="CDH1-low (deep)", low_label="CDH1-high (shallow)",
                        predicted_direction="high_better")   # CDH1-high (shallow) predicted better
    results.append(res)

    # ── S2-P3: ILC vs LumA ───────────────────────────────────
    log("  S2-P3 — ILC vs LumA OVERALL SURVIVAL")
    log("  Prediction: ILC worse than LumA")
    log("-" * 65)
    luma_s = pops.get("luma", [])
    if len(luma_s) >= 10:
        res = survival_test("ILC vs LumA", "S2-P3",
                            luma_s, ilc_s, surv_map,
                            high_label="LumA", low_label="ILC",
                            predicted_direction="high_better")  # LumA predicted better
        results.append(res)
    else:
        log("  S2-P3: Insufficient LumA samples. Skipped.")

    # ── S2-P4: SPDEF high vs low within ILC ──────────────────
    log("  S2-P4 — SPDEF NOVEL PROGNOSTIC MARKER WITHIN ILC")
    log("  Prediction: SPDEF-high = better survival (novel)")
    log("-" * 65)
    hi, lo = tertile_split(expr, "SPDEF", ilc_s)
    res = survival_test("SPDEF tertile within ILC", "S2-P4",
                        hi, lo, surv_map,
                        high_label="SPDEF-high", low_label="SPDEF-low",
                        predicted_direction="high_better")
    results.append(res)

    # ── S2-P5: PTEN low vs high within ILC ───────────────────
    log("  S2-P5 — PTEN-LOW PREDICTS WORSE OS WITHIN ILC")
    log("  Prediction: PTEN-low = worse survival")
    log("-" * 65)
    hi_pten, lo_pten = tertile_split(expr, "PTEN", ilc_s)
    res = survival_test("PTEN tertile within ILC", "S2-P5",
                        lo_pten, hi_pten, surv_map,
                        high_label="PTEN-low", low_label="PTEN-high",
                        predicted_direction="high_better")  # PTEN-high predicted better
    results.append(res)

    # ── S2-P6: EZH2 (null prediction) ────────────────────────
    log("  S2-P6 — EZH2 DOES NOT PREDICT OS WITHIN ILC (null)")
    log("  Prediction: no significant difference")
    log("-" * 65)
    hi, lo = tertile_split(expr, "EZH2", ilc_s)
    res = survival_test("EZH2 tertile within ILC", "S2-P6",
                        hi, lo, surv_map,
                        high_label="EZH2-high", low_label="EZH2-low",
                        predicted_direction="null")
    results.append(res)

    # ── S2-P7: CCND1 (weak positive prediction) ──────────────
    log("  S2-P7 — CCND1-HIGH TRENDS TOWARD BETTER OS WITHIN ILC")
    log("  Prediction: trend p < 0.15, CCND1-high better (treatment confound)")
    log("-" * 65)
    hi, lo = tertile_split(expr, "CCND1", ilc_s)
    res = survival_test("CCND1 tertile within ILC", "S2-P7",
                        hi, lo, surv_map,
                        high_label="CCND1-high", low_label="CCND1-low",
                        predicted_direction="high_better")
    results.append(res)

    # ── S2-P8: MKI67 (conditional null) ──────────────────────
    log("  S2-P8 — MKI67 DOES NOT STRONGLY PREDICT OS IN ILC (conditional null)")
    log("  Prediction: p > 0.05 within ILC, or weak HR")
    log("-" * 65)
    hi, lo = tertile_split(expr, "MKI67", ilc_s)
    res = survival_test("MKI67 tertile within ILC", "S2-P8",
                        hi, lo, surv_map,
                        high_label="MKI67-high", low_label="MKI67-low",
                        predicted_direction="null")
    results.append(res)

    # ── BONUS: PIK3CA mRNA tertile within ILC ────────────────
    log("  BONUS — PIK3CA mRNA TERTILE WITHIN ILC")
    log("  Not a formal prediction. Exploratory.")
    log("-" * 65)
    hi, lo = tertile_split(expr, "PIK3CA", ilc_s)
    res = survival_test("PIK3CA mRNA tertile within ILC", "BONUS-PIK3CA",
                        hi, lo, surv_map,
                        high_label="PIK3CA-high", low_label="PIK3CA-low",
                        predicted_direction="null")
    results.append(res)

    # ── BONUS: FOXA1 within ILC ───────────────────────────────
    log("  BONUS — FOXA1 TERTILE WITHIN ILC")
    log("  Exploratory. FOXA1 was highest luminal TF by absolute level.")
    log("-" * 65)
    hi, lo = tertile_split(expr, "FOXA1", ilc_s)
    res = survival_test("FOXA1 tertile within ILC", "BONUS-FOXA1",
                        hi, lo, surv_map,
                        high_label="FOXA1-high", low_label="FOXA1-low",
                        predicted_direction="high_better")
    results.append(res)

    # ── Cross-subtype survival ────────────────────────────────
    log("")
    log("=" * 65)
    log("  CROSS-SUBTYPE SURVIVAL COMPARISON")
    log("  ILC vs LumA vs TNBC vs HER2 vs LumB")
    log("=" * 65)

    subtype_surv = {}
    for grp in ["ilc", "luma", "tnbc", "her2", "lumb"]:
        samps = pops.get(grp, [])
        t, e  = get_group_survival(samps, surv_map)
        if len(t) >= 5:
            med = median_os(t, e)
            n_ev = int(e.sum())
            subtype_surv[grp] = {"t": t, "e": e, "n": len(t),
                                  "events": n_ev, "median_os": med}
            log(f"  {grp.upper():<12} n={len(t):<5} events={n_ev:<5} "
                f"median OS={med:.0f}d")

    # ILC vs each other subtype
    log("")
    log("  ILC pairwise log-rank:")
    ilc_data = subtype_surv.get("ilc")
    if ilc_data:
        for grp in ["luma", "tnbc", "her2", "lumb"]:
            gdata = subtype_surv.get(grp)
            if gdata is None:
                continue
            chi2, p = logrank_test(ilc_data["t"], ilc_data["e"],
                                   gdata["t"],    gdata["e"])
            hr = hr_estimate(ilc_data["t"], ilc_data["e"],
                             gdata["t"],    gdata["e"])
            log(f"  ILC vs {grp.upper():<8} {fmt_p(p)}  HR(ILC/other)={hr:.3f}")

    return results, subtype_surv

# ============================================================
# FIGURE
# ============================================================

def generate_figure(results, subtype_surv, pops, expr, surv_map):
    log("")
    log("=" * 65)
    log("GENERATING FIGURE")
    log("=" * 65)

    COLORS = {"ilc":  "#8B0000", "luma": "#4169E1", "tnbc": "#FF4500",
              "her2": "#9400D3",  "lumb": "#FF8C00",
              "high": "#2E8B57",  "low":  "#DC143C"}

    fig = plt.figure(figsize=(22, 18))
    fig.suptitle(
        "BRCA ILC — Script 2 Survival  |  OrganismCore BRCA-S6d  |  2026-03-05\n"
        "Attractor Type: TYPE 3 VARIANT — ADHESION LOCK DISSOLUTION",
        fontsize=12, fontweight="bold", y=0.99
    )
    gs_fig = gridspec.GridSpec(3, 3, figure=fig, hspace=0.50, wspace=0.42)

    def km_panel(ax, res, title, color_hi=None, color_lo=None):
        if res is None or res.get("status") == "INSUFFICIENT_N":
            ax.text(0.5, 0.5, "Insufficient data", ha="center",
                    va="center", transform=ax.transAxes, fontsize=10)
            ax.set_title(title, fontsize=9)
            return
        t_hi, e_hi = res["t_high"], res["e_high"]
        t_lo, e_lo = res["t_low"],  res["e_low"]
        if len(t_hi) == 0 or len(t_lo) == 0:
            ax.text(0.5, 0.5, "No data", ha="center",
                    va="center", transform=ax.transAxes, fontsize=10)
            ax.set_title(title, fontsize=9)
            return
        ch = color_hi or COLORS["high"]
        cl = color_lo or COLORS["low"]
        t1_km, s1_km = kaplan_meier(t_hi, e_hi)
        t2_km, s2_km = kaplan_meier(t_lo, e_lo)
        ax.step(t1_km / 365, s1_km, where="post", color=ch, lw=1.8,
                label=f"{res['high_label']} (n={res['n_high']})")
        ax.step(t2_km / 365, s2_km, where="post", color=cl, lw=1.8,
                label=f"{res['low_label']} (n={res['n_low']})")
        ax.set_xlabel("Years", fontsize=8)
        ax.set_ylabel("Survival probability", fontsize=8)
        ax.set_ylim(0, 1.05)
        p = res["p"]
        hr = res["hr"]
        status = res["status"]
        ax.set_title(f"{title}\n{fmt_p(p)}  HR={hr:.2f}  [{status}]", fontsize=8)
        ax.legend(fontsize=7)

    # Get results by pred_id
    def get_res(pid):
        for r in results:
            if r.get("pred_id") == pid:
                return r
        return None

    # A — S2-P1: ESR1
    ax = fig.add_subplot(gs_fig[0, 0])
    km_panel(ax, get_res("S2-P1"), "A — S2-P1: ESR1 within ILC\n"
             "Predicted: High ESR1 = better OS")

    # B — S2-P2: CDH1 depth
    ax = fig.add_subplot(gs_fig[0, 1])
    km_panel(ax, get_res("S2-P2"), "B — S2-P2: CDH1 depth within ILC\n"
             "Predicted: Low CDH1 (deep) = worse OS",
             color_hi=COLORS["low"], color_lo=COLORS["high"])

    # C — S2-P3: ILC vs LumA
    ax = fig.add_subplot(gs_fig[0, 2])
    km_panel(ax, get_res("S2-P3"), "C — S2-P3: ILC vs LumA\n"
             "Predicted: ILC worse than LumA",
             color_hi=COLORS["luma"], color_lo=COLORS["ilc"])

    # D — S2-P4: SPDEF (novel)
    ax = fig.add_subplot(gs_fig[1, 0])
    km_panel(ax, get_res("S2-P4"), "D — S2-P4: SPDEF within ILC (NOVEL)\n"
             "Predicted: SPDEF-high = better OS")

    # E — S2-P5: PTEN
    ax = fig.add_subplot(gs_fig[1, 1])
    km_panel(ax, get_res("S2-P5"), "E — S2-P5: PTEN within ILC\n"
             "Predicted: PTEN-low = worse OS",
             color_hi=COLORS["low"], color_lo=COLORS["high"])

    # F — Cross-subtype KM
    ax = fig.add_subplot(gs_fig[1, 2])
    if subtype_surv:
        for grp, data in subtype_surv.items():
            t_km, s_km = kaplan_meier(data["t"], data["e"])
            ax.step(t_km / 365, s_km, where="post",
                    color=COLORS.get(grp, "gray"), lw=1.5,
                    label=f"{grp.upper()} (n={data['n']})")
        ax.set_xlabel("Years", fontsize=8)
        ax.set_ylabel("Survival probability", fontsize=8)
        ax.set_ylim(0, 1.05)
        ax.set_title("F — Cross-Subtype Survival\nS2-P3: ILC vs LumA comparison", fontsize=9)
        ax.legend(fontsize=7)

    # G — S2-P6: EZH2 (null)
    ax = fig.add_subplot(gs_fig[2, 0])
    km_panel(ax, get_res("S2-P6"), "G — S2-P6: EZH2 within ILC\n"
             "Predicted: NULL (no significant difference)")

    # H — S2-P7 + P8: CCND1 and MKI67 comparison bar
    ax = fig.add_subplot(gs_fig[2, 1])
    pred_ids  = ["S2-P7", "S2-P8", "BONUS-FOXA1"]
    pred_labs = ["CCND1", "MKI67", "FOXA1"]
    p_vals    = []
    hrs       = []
    for pid in pred_ids:
        r = get_res(pid)
        if r and not np.isnan(r.get("p", float("nan"))):
            p_vals.append(-np.log10(r["p"] + 1e-10))
            hrs.append(r.get("hr", 1.0))
        else:
            p_vals.append(0)
            hrs.append(1.0)
    x = np.arange(len(pred_labs))
    bar_colors = [COLORS["high"] if h < 1 else COLORS["low"] for h in hrs]
    ax.bar(x, p_vals, color=bar_colors, alpha=0.85)
    ax.axhline(-np.log10(0.05), color="black", lw=1, ls="--", label="p=0.05")
    ax.axhline(-np.log10(0.15), color="gray", lw=1, ls=":", label="p=0.15")
    ax.set_xticks(x)
    ax.set_xticklabels(pred_labs, fontsize=9)
    ax.set_ylabel("-log10(p)", fontsize=8)
    ax.set_title("H — P7/P8/FOXA1 Summary\n-log10(p) | green=high better, red=low better",
                 fontsize=9)
    ax.legend(fontsize=8)

    # I — Prediction scorecard
    ax = fig.add_subplot(gs_fig[2, 2])
    ax.axis("off")
    scorecard_lines = [
        "PREDICTION SCORECARD",
        "─" * 32,
    ]
    status_color = {
        "CONFIRMED":        "darkgreen",
        "TREND":            "olivedrab",
        "NULL CONFIRMED":   "darkgreen",
        "NOT CONFIRMED":    "firebrick",
        "OPPOSITE":         "darkred",
        "INSUFFICIENT_N":   "gray",
    }
    for r in results:
        pid    = r.get("pred_id", "?")
        status = r.get("status", "?")
        p      = r.get("p", float("nan"))
        p_str  = fmt_p(p) if not np.isnan(p) else "n/a"
        line   = f"{pid:<14} {status:<18} {p_str}"
        scorecard_lines.append(line)
    ax.text(0.02, 0.98, "\n".join(scorecard_lines),
            transform=ax.transAxes, fontsize=7.5,
            verticalalignment="top", fontfamily="monospace",
            bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))

    plt.savefig(FIG_FILE, dpi=150, bbox_inches="tight")
    plt.close()
    log(f"  Figure saved: {FIG_FILE}")

# ============================================================
# SAVE RESULTS
# ============================================================

def save_results(results):
    rows = []
    for r in results:
        rows.append({
            "pred_id":     r.get("pred_id"),
            "label":       r.get("label"),
            "status":      r.get("status"),
            "n_high":      r.get("n_high"),
            "n_low":       r.get("n_low"),
            "events_high": r.get("events_high"),
            "events_low":  r.get("events_low"),
            "p":           r.get("p"),
            "hr":          r.get("hr"),
            "med_os_high": r.get("med_os_high"),
            "med_os_low":  r.get("med_os_low"),
        })
    df = pd.DataFrame(rows)
    df.to_csv(CSV_FILE, index=False)
    log(f"  Results saved: {CSV_FILE}")


def print_final_summary(results):
    log("")
    log("=" * 65)
    log("FINAL PREDICTION SCORECARD")
    log("=" * 65)
    log(f"  {'ID':<14} {'Status':<20} {'p':>10}  {'HR':>7}  Label")
    log("  " + "-" * 65)
    for r in results:
        p_str = f"{r['p']:.4f}" if not np.isnan(r.get("p", float("nan"))) else "n/a"
        hr_str = f"{r['hr']:.3f}" if not np.isnan(r.get("hr", float("nan"))) else "n/a"
        log(f"  {r.get('pred_id',''):<14} {r.get('status',''):<20} "
            f"{p_str:>10}  {hr_str:>7}  {r.get('label','')}")

# ============================================================
# MAIN
# ============================================================

def main():
    log("=" * 65)
    log("BRCA ILC — SCRIPT 2: SURVIVAL ANALYSIS")
    log("OrganismCore — Document BRCA-S6c/d | 2026-03-05")
    log("")
    log("ATTRACTOR TYPE: TYPE 3 VARIANT — ADHESION LOCK DISSOLUTION")
    log(f"Output directory: {BASE_DIR}")
    log("=" * 65)

    # Load data
    expr                   = load_expression()
    clin, hist_col, pam50_col = load_clinical()
    surv                   = load_survival()

    # Classify
    pops     = classify_samples(expr, clin, hist_col, pam50_col)
    surv_map = build_survival_vectors(pops, clin, surv)

    # Normalise expression (same as Script 1)
    vals       = expr.values.ravel()
    max_val    = float(np.nanmax(vals))
    if max_val >= 40:
        log("  Applying log2(x+1) normalisation.")
        expr = np.log2(expr + 1)
    else:
        log("  Expression already log2-transformed. Using as-is.")

    # Analyses
    results, subtype_surv = survival_analyses(expr, pops, surv_map)

    # Output
    generate_figure(results, subtype_surv, pops, expr, surv_map)
    save_results(results)
    print_final_summary(results)

    log("")
    log("=" * 65)
    log("SCRIPT 2 COMPLETE")
    log("=" * 65)
    log(f"  Log:     {LOG_FILE}")
    log(f"  Figure:  {FIG_FILE}")
    log(f"  CSV:     {CSV_FILE}")
    log("")
    log("NEXT: Write BRCA-S6d (script2_results_and_reasoning.md)")
    log("      then run Literature Check (Phase 4 of Protocol).")

    flush_log()


if __name__ == "__main__":
    main()
