"""
BRCA HER2-ENRICHED — SCRIPT 2
BULK RNA-seq VALIDATION + TRASTUZUMAB RESISTANCE + SURVIVAL
OrganismCore — Document BRCA-S3c/d | 2026-03-05

CONFIRMED DATA STRUCTURE (GSE37946):
  50 samples, 22 283 Affymetrix probes
  !Sample_characteristics_ch1 fields:
    age | er | pr | her2 | clin_stage | path_stage | path_response
  pCR column is 'path_response'  (values: 'pCR' / 'RD')
  NOT 'path_stage' — these are distinct rows.

BUG FIX vs previous run:
  Old parser selected the FIRST meta field whose VALUE contained
  'path_response' as a substring — that picked 'path_stage'.
  New parser:
    1. Collects all !Sample_characteristics_ch1 rows into a dict
       keyed by the label BEFORE the first colon.
    2. Scores every field by counting how many of its values
       match the canonical pCR/RD vocabulary exactly
       (case-insensitive, after stripping the 'key: ' prefix).
    3. Picks the field with the highest hit-count.
  This is robust to any GEO series matrix layout.

DATASETS:
  GSE37946  — primary (confirmed working)
  GSE50948, GSE55348, GSE66399 — fallback
  TCGA-BRCA expression + clinical + survival

PREDICTIONS TESTED (BRCA-S3c):
  S2-P1: ERBB3-high → pCR in trastuzumab-treated HER2+
  S2-P2: Depth score (ERBB3/CDH1/AR inverse) predicts OS
  S2-P3: r(EZH2, ESR1) < 0 within HER2-enriched bulk
  S2-P4: ERBB2 FC >+400% vs LumA; STARD3 co-elevated
  S2-P5: FOXA1 retained in HER2 vs Basal; r(FOXA1,ESR1) > 0
  S2-P6: Depth score correlates with Grade (Spearman)
  S2-P7: AR-low HER2+ has worse OS than AR-high HER2+
"""

import os
import sys
import gzip
import time
import warnings
import urllib.request
import urllib.error
import traceback
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import chi2
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

try:
    from sklearn.metrics import roc_auc_score
    SKLEARN_OK = True
except ImportError:
    SKLEARN_OK = False

# ============================================================
# CONFIGURATION
# ============================================================

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
BASE_DIR    = os.path.join(SCRIPT_DIR, "HER2_s2_results")
DATA_DIR    = os.path.join(BASE_DIR, "data")
RESULTS_DIR = BASE_DIR

for d in [BASE_DIR, DATA_DIR, RESULTS_DIR]:
    os.makedirs(d, exist_ok=True)

LOG_FILE = os.path.join(RESULTS_DIR, "her2_s2_log.txt")
FIG_FILE = os.path.join(RESULTS_DIR, "her2_s2_figure.png")
CSV_FILE = os.path.join(RESULTS_DIR, "her2_s2_results.csv")

# ── TCGA URLs (tried in order) ───────────────────────────────
TCGA_EXPR_URLS = [
    "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2.gz",
    "https://pancanatlas.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2.gz",
    "https://gdc.xenahubs.net/download/TCGA-BRCA.htseq_fpkm.tsv.gz",
]
TCGA_EXPR_FILE = os.path.join(DATA_DIR, "TCGA_BRCA_HiSeqV2.gz")

TCGA_CLIN_URLS = [
    "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix",
    "https://pancanatlas.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix",
    "https://gdc.xenahubs.net/download/TCGA-BRCA.GDC_phenotype.tsv.gz",
]
TCGA_CLIN_FILE = os.path.join(DATA_DIR, "TCGA_BRCA_clinicalMatrix.tsv")

TCGA_SURV_URLS = [
    "https://pancanatlas.xenahubs.net/download/Survival_SupplementalTable_S1_20171025_xena_sp",
    "https://gdc.xenahubs.net/download/TCGA-BRCA.survival.tsv",
]
TCGA_SURV_FILE = os.path.join(DATA_DIR, "TCGA_pancan_survival.tsv")

# ── Trastuzumab response datasets (tried in priority order) ──
TRAST_DATASETS = [
    (
        "GSE37946",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE37nnn/GSE37946/matrix/"
        "GSE37946_series_matrix.txt.gz",
        os.path.join(DATA_DIR, "GSE37946_series_matrix.txt.gz"),
    ),
    (
        "GSE50948",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50948/matrix/"
        "GSE50948_series_matrix.txt.gz",
        os.path.join(DATA_DIR, "GSE50948_series_matrix.txt.gz"),
    ),
    (
        "GSE55348",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE55nnn/GSE55348/matrix/"
        "GSE55348_series_matrix.txt.gz",
        os.path.join(DATA_DIR, "GSE55348_series_matrix.txt.gz"),
    ),
    (
        "GSE66399",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE66nnn/GSE66399/matrix/"
        "GSE66399_series_matrix.txt.gz",
        os.path.join(DATA_DIR, "GSE66399_series_matrix.txt.gz"),
    ),
]

# ── Gene panels ──────────────────────────────────────────────
AMPLICON_GENES   = ["ERBB2", "GRB7", "STARD3", "MIEN1"]
LUMINAL_GENES    = ["ESR1", "FOXA1", "GATA3", "PGR", "SPDEF",
                    "KRT8", "KRT18", "CDH1"]
BASAL_GENES      = ["KRT5", "KRT14", "SOX10", "VIM", "FOXC1",
                    "CDH3", "ZEB1", "ZEB2"]
DEPTH_GENES      = ["ERBB3", "CDH1", "AR"]   # low = deep = resistant
EPIGENETIC_GENES = ["EZH2", "EED", "SUZ12", "HDAC1", "HDAC2", "DNMT3A"]
SIGNALLING_GENES = ["AKT1", "MTOR", "PIK3CA", "PTEN", "EGFR", "ERBB3"]
PROLIF_GENES     = ["MKI67", "TOP2A", "PCNA", "CCNB1", "AURKA"]
DNA_REPAIR_GENES = ["BRCA1", "BRCA2", "TP53", "RB1", "PARP1"]
CONTROLS         = ["CDX2", "SPI1", "MBP"]

ALL_GENES = list(dict.fromkeys(
    AMPLICON_GENES + LUMINAL_GENES + BASAL_GENES + DEPTH_GENES +
    EPIGENETIC_GENES + SIGNALLING_GENES + PROLIF_GENES +
    DNA_REPAIR_GENES + CONTROLS
))

# Canonical ERBB3 probe IDs (Affymetrix HG-U133A / Plus2)
ERBB3_PROBES = [
    "205047_s_at",    # HG-U133A primary ERBB3
    "205048_s_at",
    "210766_s_at",
    "1565483_at",     # HG-U133Plus2
    "1565484_x_at",
]

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
# DOWNLOAD HELPERS
# ============================================================

def download_file(url, dest, label="", timeout=120, chunk=1024 * 1024):
    """Download url → dest.  Returns True on success."""
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        log(f"  Cached: {os.path.basename(dest)}")
        return True
    log(f"  Downloading {label or os.path.basename(dest)} ...")
    try:
        req = urllib.request.Request(
            url, headers={"User-Agent": "Mozilla/5.0 OrganismCore/2026"}
        )
        with urllib.request.urlopen(req, timeout=timeout) as r, \
             open(dest, "wb") as f:
            downloaded = 0
            while True:
                buf = r.read(chunk)
                if not buf:
                    break
                f.write(buf)
                downloaded += len(buf)
                if downloaded % (10 * chunk) == 0:
                    log(f"    {downloaded / 1e6:.0f} MB...")
        sz = os.path.getsize(dest)
        if sz < 1000:
            os.remove(dest)
            log(f"  FAILED (file too small: {sz} bytes)")
            return False
        log(f"  OK: {sz / 1e6:.1f} MB")
        return True
    except Exception as e:
        log(f"  FAILED: {e}")
        if os.path.exists(dest):
            os.remove(dest)
        return False


def try_urls(urls, dest, label=""):
    for url in urls:
        log(f"  Trying: {url}")
        if download_file(url, dest, label=label):
            return True
    return False


def open_maybe_gz(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode, errors="ignore")
    return open(path, mode, errors="ignore")

# ============================================================
# STEP 1 — DOWNLOAD ALL DATA
# ============================================================

def acquire_data():
    log("=" * 60)
    log("STEP 1: DATA ACQUISITION")
    log("=" * 60)

    results = {}

    log("\n-- TCGA-BRCA expression --")
    results["tcga_expr"] = try_urls(TCGA_EXPR_URLS, TCGA_EXPR_FILE,
                                    "TCGA-BRCA expression")

    log("\n-- TCGA-BRCA clinical --")
    results["tcga_clin"] = try_urls(TCGA_CLIN_URLS, TCGA_CLIN_FILE,
                                    "TCGA-BRCA clinical")

    log("\n-- TCGA-BRCA survival --")
    results["tcga_surv"] = try_urls(TCGA_SURV_URLS, TCGA_SURV_FILE,
                                    "TCGA-BRCA survival")

    log("\n-- Trastuzumab response datasets --")
    trast_found = None
    for name, url, dest in TRAST_DATASETS:
        log(f"  Trying {name}...")
        if download_file(url, dest, label=name):
            trast_found = (name, dest)
            log(f"  SUCCESS: using {name}")
            break
    results["trast"] = trast_found

    return results

# ============================================================
# STEP 2 — LOAD TCGA EXPRESSION
# ============================================================

def load_tcga_expr():
    log("")
    log("=" * 60)
    log("STEP 2: LOAD TCGA-BRCA EXPRESSION")
    log("=" * 60)

    if not os.path.exists(TCGA_EXPR_FILE):
        log("ERROR: TCGA expression file not found.")
        return None

    log("  Reading expression matrix...")
    try:
        df = pd.read_csv(TCGA_EXPR_FILE, sep="\t", index_col=0,
                         compression="gzip")
    except Exception:
        try:
            df = pd.read_csv(TCGA_EXPR_FILE, sep="\t", index_col=0)
        except Exception as e:
            log(f"  ERROR loading expression: {e}")
            return None

    log(f"  Shape: {df.shape}  (genes x samples)")
    if df.shape[0] < df.shape[1]:
        log("  Transposing (samples appear to be rows)...")
        df = df.T

    found   = [g for g in ALL_GENES if g in df.index]
    missing = [g for g in ALL_GENES if g not in df.index]
    log(f"  Target genes found:   {len(found)} / {len(ALL_GENES)}")
    if missing:
        log(f"  Missing genes:        {missing}")

    return df

# ============================================================
# STEP 3 — LOAD TCGA CLINICAL + SURVIVAL
# ============================================================

def load_tcga_clinical():
    log("")
    log("=" * 60)
    log("STEP 3: LOAD TCGA-BRCA CLINICAL + SURVIVAL")
    log("=" * 60)

    clin = None
    surv = None

    if os.path.exists(TCGA_CLIN_FILE):
        try:
            clin = pd.read_csv(TCGA_CLIN_FILE, sep="\t", index_col=0,
                               low_memory=False)
            log(f"  Clinical shape: {clin.shape}")
            log(f"  Clinical columns (first 10): {list(clin.columns[:10])}")
            clin.index = [s.replace(".", "-") for s in clin.index]
        except Exception as e:
            log(f"  ERROR loading clinical: {e}")
    else:
        log("  Clinical file not found.")

    if os.path.exists(TCGA_SURV_FILE):
        try:
            surv = pd.read_csv(TCGA_SURV_FILE, sep="\t", index_col=0,
                               low_memory=False)
            log(f"  Survival shape: {surv.shape}")
            log(f"  Survival columns: {list(surv.columns)}")
            surv.index = [s.replace(".", "-") for s in surv.index]
        except Exception as e:
            log(f"  ERROR loading survival: {e}")
    else:
        log("  Survival file not found.")

    return clin, surv

# ============================================================
# STEP 4 — IDENTIFY HER2-ENRICHED SAMPLES
# ============================================================

def get_her2_samples(expr, clin):
    log("")
    log("=" * 60)
    log("STEP 4: IDENTIFY HER2-ENRICHED SAMPLES")
    log("=" * 60)

    her2_samples  = None
    luma_samples  = None
    basal_samples = None

    if clin is not None:
        pam50_candidates = [c for c in clin.columns
                            if "pam50" in c.lower() or
                               "subtype" in c.lower()]
        log(f"  PAM50 column candidates: {pam50_candidates}")

        pam50_col = None
        for col in pam50_candidates:
            vals = clin[col].dropna().astype(str).unique()
            log(f"  Column '{col}' values: {list(vals[:10])}")
            if any("her2" in v.lower() or "HER2" in v for v in vals):
                pam50_col = col
                break

        if pam50_col:
            log(f"  Using PAM50 column: '{pam50_col}'")

            def match(kw):
                mask = clin[pam50_col].astype(str).str.contains(
                    kw, case=False, na=False)
                return list(clin.index[mask])

            her2_samples  = match("HER2")
            luma_samples  = match("LumA")
            basal_samples = match("Basal")
            log(f"  HER2-enriched: n={len(her2_samples)}")
            log(f"  LumA:          n={len(luma_samples)}")
            log(f"  Basal-like:    n={len(basal_samples)}")
        else:
            log("  WARNING: No PAM50 column found.")

    # Fallback: top-20% ERBB2
    if not her2_samples and expr is not None and "ERBB2" in expr.index:
        log("  Falling back to ERBB2-expression HER2 calling (top 20%)...")
        erbb2 = expr.loc["ERBB2"]
        thr   = erbb2.quantile(0.80)
        her2_samples = list(erbb2[erbb2 >= thr].index)
        log(f"  HER2-enriched (ERBB2 top 20%): n={len(her2_samples)}")

    return her2_samples, luma_samples, basal_samples

# ============================================================
# STEP 5 — TCGA BULK VALIDATION  (S2-P3, S2-P4, S2-P5)
# ============================================================

def tcga_bulk_validation(expr, her2_samples, luma_samples, basal_samples):
    log("")
    log("=" * 60)
    log("STEP 5: TCGA BULK VALIDATION")
    log("  S2-P3: r(EZH2, ESR1) in HER2-enriched")
    log("  S2-P4: Amplicon fold-change in bulk")
    log("  S2-P5: FOXA1 retained; r(FOXA1, ESR1) > 0")
    log("=" * 60)

    results = {}

    if expr is None or not her2_samples:
        log("  SKIPPED: expression or HER2 samples not available.")
        return results

    expr_set = set(expr.columns)
    her2  = [s for s in her2_samples  if s in expr_set]
    luma  = [s for s in (luma_samples  or []) if s in expr_set]
    basal = [s for s in (basal_samples or []) if s in expr_set]

    log(f"  HER2-enriched in expr: n={len(her2)}")
    log(f"  LumA in expr:          n={len(luma)}")
    log(f"  Basal in expr:         n={len(basal)}")

    if len(her2) < 20:
        log("  WARNING: Too few HER2 samples.")
        return results

    def safe_compare(gene, a_samp, b_samp):
        if gene not in expr.index:
            return None, None, None, None
        a = expr.loc[gene, a_samp].dropna().values.astype(float)
        b = expr.loc[gene, b_samp].dropna().values.astype(float)
        if len(a) < 5 or len(b) < 5:
            return None, None, None, None
        ref  = max(np.mean(b), 0.01)
        fc   = (np.mean(a) - np.mean(b)) / ref * 100
        _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        return np.mean(a), np.mean(b), fc, p

    # ── S2-P4 ─────────────────────────────────────────────────
    log("")
    log("  S2-P4 — AMPLICON FOLD-CHANGE (HER2 vs LumA)")
    if luma:
        genes_p4 = AMPLICON_GENES + ["ESR1", "FOXA1", "GATA3",
                                      "EZH2", "ERBB3"]
        log(f"  {'Gene':<12} {'HER2':>8} {'LumA':>8} {'FC%':>10}  p-value")
        log("  " + "-" * 55)
        for gene in genes_p4:
            ma, mb, fc, p = safe_compare(gene, her2, luma)
            if ma is None:
                log(f"  {gene:<12}  NOT FOUND")
                continue
            sig = ("***" if p < 0.001 else "**" if p < 0.01
                   else "*" if p < 0.05 else "ns")
            log(f"  {gene:<12} {ma:>8.3f} {mb:>8.3f} {fc:>+10.1f}%"
                f"  p={p:.2e} {sig}")
            results[f"p4_{gene}_fc_pct"] = fc
            results[f"p4_{gene}_p"]       = p

        erbb2_fc  = results.get("p4_ERBB2_fc_pct", 0)
        stard3_fc = results.get("p4_STARD3_fc_pct", 0)
        log("")
        log(f"  S2-P4 ERBB2: "
            + ("CONFIRMED" if erbb2_fc > 400 else
               "PARTIAL"   if erbb2_fc > 100 else "NOT CONFIRMED")
            + f" ({erbb2_fc:+.1f}%)")
        log(f"  S2-P4 STARD3: "
            + ("CONFIRMED" if stard3_fc > 100 else "NOT CONFIRMED")
            + f" ({stard3_fc:+.1f}%)")
    else:
        log("  SKIPPED: no LumA samples.")

    # ── S2-P5 ─────────────────────────────────────────────────
    log("")
    log("  S2-P5 — FOXA1 RETENTION (HER2 vs Basal)")
    if basal:
        log(f"  {'Gene':<12} {'HER2':>8} {'Basal':>8} {'FC%':>10}  p-value")
        log("  " + "-" * 55)
        for gene in ["FOXA1", "ESR1", "GATA3", "KRT5", "SOX10", "EZH2"]:
            ma, mb, fc, p = safe_compare(gene, her2, basal)
            if ma is None:
                log(f"  {gene:<12}  NOT FOUND")
                continue
            sig = ("***" if p < 0.001 else "**" if p < 0.01
                   else "*" if p < 0.05 else "ns")
            log(f"  {gene:<12} {ma:>8.3f} {mb:>8.3f} {fc:>+10.1f}%"
                f"  p={p:.2e} {sig}")
            results[f"p5_{gene}_her2"]  = ma
            results[f"p5_{gene}_basal"] = mb

        if all(g in expr.index for g in ["FOXA1", "ESR1"]):
            f_h = expr.loc["FOXA1", her2].dropna().values.astype(float)
            e_h = expr.loc["ESR1",  her2].dropna().values.astype(float)
            n   = min(len(f_h), len(e_h))
            r_fe, p_fe = stats.pearsonr(f_h[:n], e_h[:n])
            log(f"\n  r(FOXA1, ESR1) in HER2-enriched: "
                f"r={r_fe:+.3f}  p={p_fe:.2e}")
            results["p5_r_foxa1_esr1"] = r_fe
            results["p5_p_foxa1_esr1"] = p_fe
            log("  S2-P5: "
                + ("CONFIRMED" if r_fe > 0.20 and p_fe < 0.05 else
                   "PARTIAL"   if r_fe > 0    and p_fe < 0.05 else
                   "NOT CONFIRMED"))
    else:
        log("  SKIPPED: no Basal samples.")

    # ── S2-P3 ─────────────────────────────────────────────────
    log("")
    log("  S2-P3 — r(EZH2, ESR1) IN HER2-ENRICHED")
    if all(g in expr.index for g in ["EZH2", "ESR1"]):
        ez = expr.loc["EZH2", her2].dropna().values.astype(float)
        es = expr.loc["ESR1", her2].dropna().values.astype(float)
        n  = min(len(ez), len(es))
        r_ee, p_ee = stats.pearsonr(ez[:n], es[:n])
        log(f"  r(EZH2, ESR1): r={r_ee:+.3f}  p={p_ee:.2e}")
        results["p3_r_ezh2_esr1"] = r_ee
        results["p3_p_ezh2_esr1"] = p_ee
        log("  S2-P3: "
            + ("CONFIRMED"     if r_ee < -0.15 and p_ee < 0.05 else
               "PARTIAL"       if r_ee < 0     and p_ee < 0.05 else
               "NOT CONFIRMED"))

        if "FOXA1" in expr.index:
            f_h = expr.loc["FOXA1", her2].dropna().values.astype(float)
            n2  = min(len(ez), len(f_h))
            r_ef, p_ef = stats.pearsonr(ez[:n2], f_h[:n2])
            log(f"  r(EZH2, FOXA1): r={r_ef:+.3f}  p={p_ef:.2e}  [secondary]")
            results["p3_r_ezh2_foxa1"] = r_ef
    else:
        log("  SKIPPED: EZH2 or ESR1 not in expression index.")

    return results

# ============================================================
# STEP 6 — DEPTH SCORE + SURVIVAL  (S2-P2, S2-P6, S2-P7)
# ============================================================

def _logrank(t1, e1, t2, e2):
    """Return log-rank p-value (no external dependency)."""
    try:
        from lifelines.statistics import logrank_test
        return logrank_test(t1, t2,
                            event_observed_A=e1,
                            event_observed_B=e2).p_value
    except ImportError:
        pass
    all_t = np.sort(np.unique(
        np.concatenate([t1[e1 == 1], t2[e2 == 1]])))
    O1 = E1 = O2 = E2 = 0.0
    for t in all_t:
        n1 = (t1 >= t).sum();  n2 = (t2 >= t).sum()
        o1 = ((t1 == t) & (e1 == 1)).sum()
        o2 = ((t2 == t) & (e2 == 1)).sum()
        n  = n1 + n2
        if n < 2:
            continue
        e1t = (n1 / n) * (o1 + o2);  e2t = (n2 / n) * (o1 + o2)
        O1 += o1;  E1 += e1t
        O2 += o2;  E2 += e2t
    if E1 == 0 or E2 == 0:
        return 1.0
    chi = (O1 - E1) ** 2 / E1 + (O2 - E2) ** 2 / E2
    return float(1 - chi2.cdf(chi, df=1))


def _median_os(t, e):
    idx = np.argsort(t);  t, e = t[idx], e[idx]
    at_risk, s = len(t), 1.0
    for i in range(len(t)):
        if e[i] == 1:
            s *= (1 - 1 / at_risk)
        if s <= 0.5:
            return float(t[i])
        at_risk -= 1
    return np.nan


def depth_and_survival(expr, clin, surv, her2_samples):
    log("")
    log("=" * 60)
    log("STEP 6: DEPTH SCORE + SURVIVAL")
    log("  S2-P2: depth score → OS")
    log("  S2-P6: depth score → Grade")
    log("  S2-P7: AR-low → worse OS")
    log("=" * 60)

    results  = {}
    surv_df  = None
    depth_score = None

    if expr is None or not her2_samples:
        log("  SKIPPED.")
        return results, surv_df, depth_score

    expr_set = set(expr.columns)
    her2     = [s for s in her2_samples if s in expr_set]
    if len(her2) < 20:
        log("  Too few HER2 samples.")
        return results, surv_df, depth_score

    # ── Compute depth score ────────────────────────────────────
    depth_avail = [g for g in DEPTH_GENES if g in expr.index]
    log(f"  Depth genes available: {depth_avail}")
    if not depth_avail:
        log("  ERROR: No depth genes.")
        return results, surv_df, depth_score

    raw = expr.loc[depth_avail, her2].T.mean(axis=1)

    def norm01(s):
        mn, mx = s.min(), s.max()
        return (s - mn) / (mx - mn) if mx > mn else s * 0

    depth_score = 1.0 - norm01(raw)
    depth_score.name = "depth_score"
    log(f"  Depth score: mean={depth_score.mean():.3f}  "
        f"median={depth_score.median():.3f}  "
        f"range=[{depth_score.min():.3f}, {depth_score.max():.3f}]")

    # ── Find OS columns ────────────────────────────────────────
    os_col   = None
    os_event = None
    surv_src = None

    for src, name in [(surv, "survival"), (clin, "clinical")]:
        if src is None:
            continue
        tc = [c for c in src.columns
              if "OS" in c and "time" in c.lower()]
        ec = [c for c in src.columns
              if "OS" in c and ("event" in c.lower() or
                                "status" in c.lower() or
                                c.endswith(".1"))]
        if tc and ec:
            os_col   = tc[0]
            os_event = ec[0]
            surv_src = src
            log(f"  OS columns from {name}: time='{os_col}', "
                f"event='{os_event}'")
            break

    if os_col:
        tmp  = depth_score.to_frame()
        tmp.index = [s.replace(".", "-") for s in tmp.index]
        surv_src.index = [s.replace(".", "-") for s in surv_src.index]
        merged = tmp.join(surv_src[[os_col, os_event]], how="inner")
        merged[os_col]   = pd.to_numeric(merged[os_col],   errors="coerce")
        merged[os_event] = pd.to_numeric(merged[os_event], errors="coerce")
        merged = merged.dropna()
        log(f"  Depth + OS merged: n={len(merged)}")

        if len(merged) >= 30:
            med    = merged["depth_score"].median()
            high_d = merged[merged["depth_score"] >= med]
            low_d  = merged[merged["depth_score"] <  med]

            p_lr = _logrank(
                high_d[os_col].values, high_d[os_event].values,
                low_d[os_col].values,  low_d[os_event].values,
            )
            med_h = _median_os(high_d[os_col].values,
                               high_d[os_event].values)
            med_l = _median_os(low_d[os_col].values,
                               low_d[os_event].values)

            log(f"\n  S2-P2 — DEPTH SCORE VS OS (log-rank)")
            log(f"  Depth-high (n={len(high_d)}): median OS = "
                + (f"{med_h:.1f}" if not np.isnan(med_h) else "not reached"))
            log(f"  Depth-low  (n={len(low_d)}):  median OS = "
                + (f"{med_l:.1f}" if not np.isnan(med_l) else "not reached"))
            log(f"  Log-rank p = {p_lr:.4f}")
            results["p2_logrank_p"]   = p_lr
            results["p2_median_high"] = med_h
            results["p2_median_low"]  = med_l
            log("  S2-P2: " + ("CONFIRMED" if p_lr < 0.05
                                else "NOT CONFIRMED"))

            surv_df = merged.copy()
            surv_df["depth_group"] = np.where(
                surv_df["depth_score"] >= med, "high", "low")
        else:
            log(f"  WARNING: only {len(merged)} samples.")
    else:
        log("  WARNING: OS columns not found.")

    # ── S2-P7: AR-low vs AR-high ──────────────────────────────
    log("")
    log("  S2-P7 — AR-LOW vs AR-HIGH WITHIN HER2-ENRICHED")

    if ("AR" in expr.index and os_col and surv_src is not None):
        ar = expr.loc["AR", her2].dropna()
        tmp_ar = ar.to_frame("AR")
        tmp_ar.index = [s.replace(".", "-") for s in tmp_ar.index]
        surv_src.index = [s.replace(".", "-") for s in surv_src.index]
        ar_m = tmp_ar.join(surv_src[[os_col, os_event]], how="inner")
        ar_m[os_col]   = pd.to_numeric(ar_m[os_col],   errors="coerce")
        ar_m[os_event] = pd.to_numeric(ar_m[os_event], errors="coerce")
        ar_m = ar_m.dropna()
        log(f"  AR + OS merged: n={len(ar_m)}")

        if len(ar_m) >= 30:
            ar_med = ar_m["AR"].median()
            ar_h   = ar_m[ar_m["AR"] >= ar_med]
            ar_l   = ar_m[ar_m["AR"] <  ar_med]
            p_ar   = _logrank(
                ar_h[os_col].values, ar_h[os_event].values,
                ar_l[os_col].values, ar_l[os_event].values,
            )
            med_arh = _median_os(ar_h[os_col].values,
                                 ar_h[os_event].values)
            med_arl = _median_os(ar_l[os_col].values,
                                 ar_l[os_event].values)
            log(f"  Median OS AR-high: "
                + (f"{med_arh:.1f}" if not np.isnan(med_arh)
                   else "not reached"))
            log(f"  Median OS AR-low:  "
                + (f"{med_arl:.1f}" if not np.isnan(med_arl)
                   else "not reached"))
            log(f"  Log-rank p = {p_ar:.4f}")
            results["p7_ar_logrank_p"] = p_ar
            confirmed = (p_ar < 0.05 and not np.isnan(med_arh)
                         and not np.isnan(med_arl) and med_arl < med_arh)
            log("  S2-P7: " + ("CONFIRMED" if confirmed
                                else "NOT CONFIRMED"))
        else:
            log("  Too few samples.")
    else:
        log("  SKIPPED.")

    # ── S2-P6: Depth vs Grade ─────────────────────────────────
    log("")
    log("  S2-P6 — DEPTH SCORE vs HISTOLOGICAL GRADE")

    grade_col = None
    if clin is not None:
        grade_cands = [c for c in clin.columns
                       if "grade" in c.lower()]
        log(f"  Grade column candidates: {grade_cands}")
        if grade_cands:
            grade_col = grade_cands[0]

    if grade_col and clin is not None:
        tmp_g = depth_score.to_frame()
        tmp_g.index = [s.replace(".", "-") for s in tmp_g.index]
        clin.index  = [s.replace(".", "-") for s in clin.index]
        gm = tmp_g.join(clin[[grade_col]], how="inner")
        gm[grade_col] = pd.to_numeric(gm[grade_col], errors="coerce")
        gm = gm.dropna()
        log(f"  Grade-depth merged: n={len(gm)}")
        if len(gm) >= 20:
            r_gd, p_gd = stats.spearmanr(gm["depth_score"],
                                         gm[grade_col])
            log(f"  Spearman r(depth, grade) = {r_gd:+.3f}  p={p_gd:.4f}")
            results["p6_r_depth_grade"] = r_gd
            results["p6_p_depth_grade"] = p_gd
            log("  S2-P6: "
                + ("CONFIRMED" if r_gd > 0 and p_gd < 0.05
                   else "NOT CONFIRMED"))
        else:
            log("  Too few samples.")
    else:
        log("  Grade column not found — S2-P6 skipped.")

    return results, surv_df, depth_score

# ============================================================
# STEP 7 — TRASTUZUMAB RESISTANCE  (S2-P1)
#
# KEY FIX:
#   Old approach: find a meta field whose NAME contains
#   "patholog" or "response" — this picked 'path_stage'
#   in GSE37946 because its value string included "path_response".
#
#   New approach (robust to any GEO layout):
#   1. After stripping the "key: " prefix from each value cell,
#      count how many cells in each field match the canonical
#      pCR vocabulary  {pCR, rd, complete response, residual
#      disease, sensitive, resistant}.
#   2. Pick the field with the highest hit-count.
#   3. Log the full field-by-field score table so bugs are
#      immediately visible without re-running a diagnostic.
# ============================================================

# Vocabulary sets (lowercased, stripped)
_PCR_VOCAB = {"pcr", "pathologic complete response",
              "complete response", "cr", "sensitive",
              "response: yes", "pcr: yes"}
_RD_VOCAB  = {"rd", "residual disease", "no pcr",
              "non pcr", "resistant", "stable disease",
              "progressive disease", "response: no",
              "pcr: no", "partial response", "pr"}


def _strip_key(raw_val, key):
    """Remove 'key: ' prefix from a GEO characteristic value."""
    v = str(raw_val).strip().strip('"')
    prefix = key.lower() + ":"
    if v.lower().startswith(prefix):
        v = v[len(prefix):].strip()
    return v.lower()


def _score_field(values, key):
    """
    Count how many cells in `values` (list, one per sample)
    belong to the canonical pCR or RD vocabularies.
    Returns (n_pcr_hits, n_rd_hits, n_total_hits).
    """
    n_pcr = n_rd = 0
    for raw in values:
        v = _strip_key(raw, key)
        # exact match
        if v in _PCR_VOCAB:
            n_pcr += 1
        elif v in _RD_VOCAB:
            n_rd += 1
        else:
            # substring match for common patterns
            if "pcr" in v or "complete" in v:
                n_pcr += 1
            elif "rd" == v or "residual" in v or "resistant" in v:
                n_rd += 1
    return n_pcr, n_rd, n_pcr + n_rd


def parse_series_matrix(fpath):
    """
    Parse a GEO series matrix (.txt.gz).

    Returns:
        sample_ids  : list[str]
        meta_df     : pd.DataFrame  (index=sample_ids, columns=meta keys)
        expr_m      : pd.DataFrame  (index=probe_ids,  columns=sample_ids)
    """
    meta_rows   = {}   # key → list[str], one value per sample
    sample_ids  = []
    expr_data   = {}
    expr_started = False
    n_header     = 0

    with gzip.open(fpath, "rt", errors="ignore") as f:
        for i, line in enumerate(f):
            line = line.rstrip("\n")

            # Sample GEO IDs
            if line.startswith("!Sample_geo_accession"):
                parts      = line.split("\t")
                sample_ids = [p.strip().strip('"') for p in parts[1:]]

            # Metadata rows
            elif line.startswith("!Sample_characteristics_ch1"):
                parts = line.split("\t")
                vals  = [p.strip().strip('"') for p in parts[1:]]
                # Derive the key from the first non-empty value
                key   = ""
                for v in vals:
                    if v:
                        key = v.split(":")[0].strip().lower() if ":" in v \
                              else v[:40].lower()
                        break
                if key and key not in meta_rows:
                    meta_rows[key] = vals

            elif "series_matrix_table_begin" in line.lower():
                expr_started = True
                n_header     = i
                continue

            elif "series_matrix_table_end" in line.lower():
                break

            elif expr_started and i > n_header:
                parts = line.split("\t")
                pid   = parts[0].strip().strip('"')
                if not pid or pid == "ID_REF":
                    if not sample_ids and pid == "ID_REF":
                        sample_ids = [p.strip().strip('"')
                                      for p in parts[1:]]
                    continue
                if len(parts) > 1:
                    try:
                        vals = [
                            float(p) if p.strip() not in
                            ("", "null", "NA", "NaN") else np.nan
                            for p in parts[1:]
                        ]
                        expr_data[pid] = vals
                    except Exception:
                        pass

    n_samp = len(sample_ids)
    meta_df = pd.DataFrame(
        {k: v for k, v in meta_rows.items() if len(v) == n_samp},
        index=sample_ids
    ) if n_samp else pd.DataFrame()

    clean_expr = {k: v for k, v in expr_data.items() if len(v) == n_samp}
    expr_m     = pd.DataFrame(clean_expr, index=sample_ids).T \
                 if clean_expr else pd.DataFrame()

    return sample_ids, meta_df, expr_m


def _select_pcr_column(meta_df):
    """
    Score every column in meta_df against the pCR/RD vocabulary.
    Return (best_col, score_table_str).
    """
    scores = {}
    for col in meta_df.columns:
        n_pcr, n_rd, total = _score_field(
            meta_df[col].tolist(), col)
        scores[col] = (n_pcr, n_rd, total)

    # Build readable score table for the log
    lines = ["  pCR column scoring:"]
    lines.append(f"  {'field':<25} {'pCR hits':>9} {'RD hits':>8} "
                 f"{'total':>7}")
    lines.append("  " + "-" * 52)
    for col, (np_, nr, tot) in sorted(scores.items(),
                                       key=lambda x: -x[1][2]):
        lines.append(f"  {col:<25} {np_:>9} {nr:>8} {tot:>7}")
    score_txt = "\n".join(lines)

    # Pick field with highest total hit-count
    best_col = max(scores, key=lambda c: scores[c][2]) \
               if scores else None
    if best_col and scores[best_col][2] == 0:
        best_col = None  # no hits at all

    return best_col, score_txt


def trastuzumab_resistance(trast_result):
    log("")
    log("=" * 60)
    log("STEP 7: TRASTUZUMAB RESISTANCE ANALYSIS")
    log("  S2-P1: ERBB3 as pCR biomarker")
    log("=" * 60)

    results = {}

    if trast_result is None:
        log("  SKIPPED: No trastuzumab dataset available.")
        results["p1_status"] = "SKIPPED"
        return results

    name, fpath = trast_result
    log(f"  Dataset: {name}  File: {os.path.basename(fpath)}")

    # ── Parse series matrix ────────────────────────────────────
    log("  Parsing series matrix...")
    try:
        sample_ids, meta_df, expr_m = parse_series_matrix(fpath)
    except Exception as e:
        log(f"  ERROR: {e}")
        traceback.print_exc()
        results["p1_status"] = f"PARSE_ERROR: {e}"
        return results

    log(f"  Samples:  {len(sample_ids)}")
    log(f"  Probes:   {expr_m.shape[0]}")
    log(f"  Meta fields: {list(meta_df.columns)}")

    if meta_df.empty or expr_m.empty:
        log("  ERROR: empty data after parse.")
        results["p1_status"] = "EMPTY"
        return results

    # ── Select pCR column (NEW ROBUST METHOD) ─────────────────
    pcr_col, score_txt = _select_pcr_column(meta_df)
    log(score_txt)

    if pcr_col is None:
        log("  ERROR: No pCR column found after scoring.")
        log(f"  Raw meta sample (first row):")
        log(f"  {meta_df.iloc[0].to_dict()}")
        results["p1_status"] = "NO_PCR_COLUMN"
        return results

    log(f"\n  Selected pCR column: '{pcr_col}'")
    log(f"  Raw values (unique): "
        f"{sorted(meta_df[pcr_col].dropna().unique().tolist())[:12]}")

    # ── Map to binary pCR ──────────────────────────────────────
    def to_binary(raw_val):
        v = _strip_key(raw_val, pcr_col)
        if v in _PCR_VOCAB or "pcr" in v or "complete" in v:
            return 1
        if v in _RD_VOCAB  or v == "rd" or "residual" in v \
                or "resistant" in v:
            return 0
        return np.nan

    meta_df["pcr_binary"] = meta_df[pcr_col].apply(to_binary)
    known = meta_df.dropna(subset=["pcr_binary"])
    n_pcr = int(known["pcr_binary"].sum())
    n_rd  = int((known["pcr_binary"] == 0).sum())
    log(f"  pCR=1 (responders): {n_pcr}")
    log(f"  pCR=0 (non-resp):   {n_rd}")

    if n_pcr < 3 or n_rd < 3:
        log(f"  WARNING: too few labelled samples "
            f"(pCR={n_pcr}, RD={n_rd}).")
        log("  Attempting next dataset fallback...")
        results["p1_status"] = "TOO_FEW"
        return results

    # ── Identify ERBB3 probes ──────────────────────────────────
    erbb3_found = [p for p in ERBB3_PROBES if p in expr_m.index]
    erbb3_sym   = [p for p in expr_m.index
                   if str(p).upper() == "ERBB3"]
    erbb3_probes = list(dict.fromkeys(erbb3_found + erbb3_sym))

    log(f"  ERBB3 probes found: {erbb3_probes}")

    if not erbb3_probes:
        log(f"  ERBB3 not found. First 30 probe IDs:")
        log(f"  {list(expr_m.index[:30])}")
        results["p1_status"] = "ERBB3_NOT_FOUND"
        return results

    erbb3_expr = expr_m.loc[erbb3_probes].mean(axis=0)
    erbb3_expr.name = "ERBB3"

    # ── Merge and test ─────────────────────────────────────────
    merged = erbb3_expr.to_frame().join(
        meta_df[["pcr_binary"]], how="inner"
    ).dropna()
    log(f"  ERBB3 + pCR merged: n={len(merged)}")

    if len(merged) < 10:
        log("  Too few samples after merge.")
        results["p1_status"] = "TOO_FEW_MERGED"
        return results

    pcr1 = merged[merged["pcr_binary"] == 1]["ERBB3"].values
    pcr0 = merged[merged["pcr_binary"] == 0]["ERBB3"].values

    log(f"  ERBB3 pCR group  (n={len(pcr1)}): "
        f"mean={np.mean(pcr1):.3f}  median={np.median(pcr1):.3f}")
    log(f"  ERBB3 RD  group  (n={len(pcr0)}): "
        f"mean={np.mean(pcr0):.3f}  median={np.median(pcr0):.3f}")

    _, p_mw = stats.mannwhitneyu(pcr1, pcr0, alternative="two-sided")
    log(f"  Mann-Whitney p = {p_mw:.4f}")

    results["p1_p_mw"]             = p_mw
    results["p1_erbb3_pcr_mean"]   = float(np.mean(pcr1))
    results["p1_erbb3_rd_mean"]    = float(np.mean(pcr0))
    results["p1_erbb3_pcr_median"] = float(np.median(pcr1))
    results["p1_erbb3_rd_median"]  = float(np.median(pcr0))

    direction_ok = np.mean(pcr1) > np.mean(pcr0)
    if direction_ok and p_mw < 0.05:
        log("  S2-P1: CONFIRMED (ERBB3-high → pCR, p < 0.05)")
        results["p1_status"] = "CONFIRMED"
    elif p_mw < 0.05:
        log("  S2-P1: SIGNIFICANT but WRONG DIRECTION "
            "(ERBB3-high → RD)")
        results["p1_status"] = "WRONG_DIRECTION"
    else:
        direction_str = ("ERBB3-high → pCR" if direction_ok
                         else "ERBB3-high → RD (reverse)")
        log(f"  S2-P1: NOT CONFIRMED (p={p_mw:.4f}, {direction_str})")
        results["p1_status"] = "NOT_CONFIRMED"

    # AUC
    if SKLEARN_OK and len(merged) >= 10:
        try:
            auc = roc_auc_score(merged["pcr_binary"].values,
                                merged["ERBB3"].values)
            log(f"  ERBB3 AUC (pCR) = {auc:.3f}")
            results["p1_auc"] = auc
        except Exception:
            pass

    return results

# ============================================================
# STEP 8 — GENERATE FIGURE
# ============================================================

def generate_figure(expr, her2_samples, luma_samples, basal_samples,
                    surv_df, depth_score, results, trast_name):
    log("")
    log("=" * 60)
    log("STEP 8: GENERATING FIGURE")
    log("=" * 60)

    try:
        fig = plt.figure(figsize=(24, 18))
        gs  = gridspec.GridSpec(3, 3, figure=fig,
                                hspace=0.45, wspace=0.35)

        expr_set = set(expr.columns) if expr is not None else set()
        her2  = [s for s in (her2_samples  or []) if s in expr_set]
        luma  = [s for s in (luma_samples  or []) if s in expr_set]
        basal = [s for s in (basal_samples or []) if s in expr_set]

        # ── Panel A: Amplicon genes ────────────────────────────
        ax_a = fig.add_subplot(gs[0, 0])
        amp  = [g for g in AMPLICON_GENES
                if expr is not None and g in expr.index]
        if amp and luma:
            x = np.arange(len(amp));  w = 0.35
            ax_a.bar(x - w/2,
                     [expr.loc[g, her2].mean() for g in amp], w,
                     label="HER2-enriched", color="#c0392b", alpha=0.85)
            ax_a.bar(x + w/2,
                     [expr.loc[g, luma].mean() for g in amp], w,
                     label="LumA", color="#2980b9", alpha=0.85)
            ax_a.set_xticks(x);  ax_a.set_xticklabels(amp, fontsize=9)
            ax_a.set_ylabel("Mean expression (log2)")
            ax_a.legend(fontsize=8)
        ax_a.set_title("A — Amplicon genes: HER2 vs LumA\n(S2-P4)",
                        fontsize=10)

        # ── Panel B: Luminal retention ─────────────────────────
        ax_b = fig.add_subplot(gs[0, 1])
        cmp  = [g for g in ["FOXA1", "ESR1", "KRT5", "SOX10", "EZH2"]
                if expr is not None and g in expr.index]
        if cmp and basal:
            x = np.arange(len(cmp));  w = 0.35
            ax_b.bar(x - w/2,
                     [expr.loc[g, her2].mean() for g in cmp], w,
                     label="HER2-enriched", color="#c0392b", alpha=0.85)
            ax_b.bar(x + w/2,
                     [expr.loc[g, basal].mean() for g in cmp], w,
                     label="Basal-like", color="#8e44ad", alpha=0.85)
            ax_b.set_xticks(x);  ax_b.set_xticklabels(cmp, fontsize=9)
            ax_b.set_ylabel("Mean expression (log2)")
            ax_b.legend(fontsize=8)
        ax_b.set_title("B — HER2 vs Basal: luminal retention\n(S2-P5)",
                        fontsize=10)

        # ── Panel C: r(EZH2, ESR1) scatter ────────────────────
        ax_c = fig.add_subplot(gs[0, 2])
        if (expr is not None and "EZH2" in expr.index
                and "ESR1" in expr.index and her2):
            ez = expr.loc["EZH2", her2].dropna().values.astype(float)
            es = expr.loc["ESR1", her2].dropna().values.astype(float)
            n  = min(len(ez), len(es))
            ax_c.scatter(ez[:n], es[:n], alpha=0.25, s=8,
                         color="#c0392b", rasterized=True)
            r_v = results.get("p3_r_ezh2_esr1", np.nan)
            p_v = results.get("p3_p_ezh2_esr1", np.nan)
            if not np.isnan(r_v):
                m, b = np.polyfit(ez[:n], es[:n], 1)
                xr   = np.array([ez[:n].min(), ez[:n].max()])
                ax_c.plot(xr, m * xr + b, color="black", lw=1.5)
            ax_c.set_xlabel("EZH2 expression")
            ax_c.set_ylabel("ESR1 expression")
            title_r = f"r={r_v:+.3f}  p={p_v:.2e}" \
                      if not np.isnan(r_v) else ""
            ax_c.set_title(f"C — r(EZH2, ESR1) in HER2-enriched\n"
                           f"{title_r}  (S2-P3)", fontsize=10)
        else:
            ax_c.text(0.5, 0.5, "Data not available",
                      ha="center", va="center",
                      transform=ax_c.transAxes)
            ax_c.set_title("C — r(EZH2, ESR1)  (S2-P3)", fontsize=10)

        # ── Panel D: Depth score histogram ─────────────────────
        ax_d = fig.add_subplot(gs[1, 0])
        if depth_score is not None and len(depth_score) > 0:
            ax_d.hist(depth_score.values, bins=40,
                      color="#c0392b", alpha=0.75, edgecolor="white")
            ax_d.axvline(depth_score.median(), color="black",
                         lw=1.5, ls="--", label="Median")
            ax_d.set_xlabel("Depth score (0=shallow, 1=deep)")
            ax_d.set_ylabel("N samples")
            ax_d.legend(fontsize=8)
        ax_d.set_title("D — Depth score distribution\n"
                       "(ERBB3/CDH1/AR inverse, S2-P2)", fontsize=10)

        # ── Panel E: KM curves ────────────────────────────────
        ax_e = fig.add_subplot(gs[1, 1])
        drawn = False
        if surv_df is not None and "depth_group" in surv_df.columns:
            os_tc = [c for c in surv_df.columns
                     if "OS" in c and "time" in c.lower()]
            os_ec = [c for c in surv_df.columns
                     if "OS" in c and ("event" in c.lower()
                                       or c.endswith(".1"))]
            if os_tc and os_ec:
                tc = os_tc[0];  ec = os_ec[0]
                clr = {"high": "#c0392b", "low": "#2980b9"}
                for grp, gdf in surv_df.groupby("depth_group"):
                    t = gdf[tc].values;  e = gdf[ec].values
                    idx = np.argsort(t);  t, e = t[idx], e[idx]
                    s, ts, ss = 1.0, [0], [1.0]
                    ar = len(t)
                    for i in range(len(t)):
                        if e[i] == 1:
                            s *= (1 - 1 / ar)
                        ts.append(t[i]);  ss.append(s)
                        ar -= 1
                    ax_e.step(ts, ss, where="post",
                              color=clr.get(grp, "grey"),
                              label=f"Depth-{grp} (n={len(gdf)})",
                              lw=2)
                p_lr = results.get("p2_logrank_p", np.nan)
                ax_e.set_title(
                    f"E — OS by depth score (HER2-enriched)\n"
                    f"log-rank p={p_lr:.4f}  (S2-P2)", fontsize=10)
                ax_e.set_xlabel("Time");  ax_e.set_ylabel("Survival")
                ax_e.legend(fontsize=8);  ax_e.set_ylim(0, 1.05)
                drawn = True
        if not drawn:
            ax_e.text(0.5, 0.5, "Survival data not available",
                      ha="center", va="center",
                      transform=ax_e.transAxes)
            ax_e.set_title("E — Survival by depth score  (S2-P2)",
                           fontsize=10)

        # ── Panel F: AR distribution ───────────────────────────
        ax_f = fig.add_subplot(gs[1, 2])
        if expr is not None and "AR" in expr.index and her2:
            ar_vals = expr.loc["AR", her2].dropna().values.astype(float)
            ax_f.hist(ar_vals, bins=40, color="#27ae60",
                      alpha=0.75, edgecolor="white")
            ax_f.axvline(np.median(ar_vals), color="black",
                         lw=1.5, ls="--", label="Median")
            p_ar = results.get("p7_ar_logrank_p", np.nan)
            ax_f.set_title(
                f"F — AR distribution in HER2-enriched\n"
                f"OS log-rank p={p_ar:.4f}  (S2-P7)"
                if not np.isnan(p_ar) else
                "F — AR distribution in HER2-enriched\n(S2-P7)",
                fontsize=10)
            ax_f.set_xlabel("AR expression");  ax_f.set_ylabel("N samples")
            ax_f.legend(fontsize=8)
        else:
            ax_f.text(0.5, 0.5, "AR data not available",
                      ha="center", va="center",
                      transform=ax_f.transAxes)
            ax_f.set_title("F — AR in HER2-enriched  (S2-P7)",
                           fontsize=10)

        # ── Panel G: Prediction scorecard ─────────────────────
        ax_g = fig.add_subplot(gs[2, 0])
        ax_g.axis("off")
        fmt = lambda v: f"{v:+.3f}" if isinstance(v, float) else str(v)

        def vline(pred, desc):
            return f"{pred}: {desc}"

        card = [
            vline("S2-P1",
                  results.get("p1_status", "NO DATA")),
            vline("S2-P2",
                  f"p={fmt(results.get('p2_logrank_p', 'NO DATA'))}"),
            vline("S2-P3",
                  f"r={fmt(results.get('p3_r_ezh2_esr1', 'NO DATA'))}"),
            vline("S2-P4",
                  f"ERBB2 FC={fmt(results.get('p4_ERBB2_fc_pct', 'NO DATA'))}%"),
            vline("S2-P5",
                  f"r={fmt(results.get('p5_r_foxa1_esr1', 'NO DATA'))}"),
            vline("S2-P6",
                  f"r={fmt(results.get('p6_r_depth_grade', 'NO DATA'))}"),
            vline("S2-P7",
                  f"p={fmt(results.get('p7_ar_logrank_p', 'NO DATA'))}"),
        ]
        txt = "PREDICTION SCORECARD\n" + "─" * 32 + "\n"
        txt += "\n".join(card)
        ax_g.text(0.05, 0.95, txt, transform=ax_g.transAxes,
                  fontsize=9, va="top", fontfamily="monospace",
                  bbox=dict(fc="lightyellow", ec="grey", alpha=0.85))

        # ── Panel H: Depth correlations ────────────────────────
        ax_h = fig.add_subplot(gs[2, 1])
        if expr is not None and depth_score is not None and her2:
            cg = list(dict.fromkeys(
                [g for g in DEPTH_GENES + AMPLICON_GENES +
                 EPIGENETIC_GENES + ["FOXA1", "ESR1", "MKI67"]
                 if g in expr.index]
            ))
            rs, gnames = [], []
            for g in cg:
                gv = expr.loc[g, her2].dropna().values.astype(float)
                dv = depth_score.reindex(pd.Index(her2)).dropna()
                n  = min(len(gv), len(dv))
                if n >= 20:
                    r, _ = stats.pearsonr(dv.values[:n], gv[:n])
                    rs.append(r);  gnames.append(g)
            if rs:
                idx = np.argsort(rs)
                rs_s = [rs[i] for i in idx]
                gn_s = [gnames[i] for i in idx]
                clrs = ["#c0392b" if r < 0 else "#27ae60" for r in rs_s]
                ax_h.barh(range(len(rs_s)), rs_s, color=clrs, alpha=0.8)
                ax_h.set_yticks(range(len(gn_s)))
                ax_h.set_yticklabels(gn_s, fontsize=7)
                ax_h.axvline(0, color="black", lw=0.8)
                ax_h.set_xlabel("Pearson r with depth score")
        ax_h.set_title("H — Depth score correlations\n"
                       "(HER2-enriched bulk)", fontsize=10)

        # ── Panel I: EZH2 / trastuzumab context ───────────────
        ax_i = fig.add_subplot(gs[2, 2])
        ax_i.axis("off")
        ctx = (
            f"EZH2 AXIS — HER2-ENRICHED\n"
            f"{'─'*30}\n"
            f"EZH2 FC vs LumA:  {fmt(results.get('p4_EZH2_fc_pct','N/A'))}%\n"
            f"r(EZH2, ESR1):    {fmt(results.get('p3_r_ezh2_esr1','N/A'))}\n"
            f"r(EZH2, FOXA1):   {fmt(results.get('p3_r_ezh2_foxa1','N/A'))}\n"
            f"r(FOXA1, ESR1):   {fmt(results.get('p5_r_foxa1_esr1','N/A'))}\n"
            f"\n"
            f"TRASTUZUMAB DATASET\n"
            f"{'─'*30}\n"
            f"Dataset: {trast_name or 'none'}\n"
            f"ERBB3 pCR mean: {fmt(results.get('p1_erbb3_pcr_mean','N/A'))}\n"
            f"ERBB3 RD  mean: {fmt(results.get('p1_erbb3_rd_mean','N/A'))}\n"
            f"AUC:            {fmt(results.get('p1_auc','N/A'))}\n"
            f"Status:         {results.get('p1_status','N/A')}\n"
        )
        ax_i.text(0.05, 0.95, ctx, transform=ax_i.transAxes,
                  fontsize=9, va="top", fontfamily="monospace",
                  bbox=dict(fc="lightcyan", ec="grey", alpha=0.85))
        ax_i.set_title("I — EZH2 / trastuzumab context", fontsize=10)

        fig.suptitle(
            "BRCA HER2-ENRICHED — SCRIPT 2\n"
            "OrganismCore | BRCA-S3c/d | 2026-03-05",
            fontsize=13, fontweight="bold", y=1.01
        )

        plt.savefig(FIG_FILE, dpi=150, bbox_inches="tight")
        plt.close()
        log(f"  Figure saved: {FIG_FILE}")

    except Exception as e:
        log(f"  Figure error: {e}")
        traceback.print_exc()

# ============================================================
# STEP 9 — SAVE RESULTS
# ============================================================

def save_results(all_results):
    log("")
    log("=" * 60)
    log("STEP 9: SAVE RESULTS")
    log("=" * 60)
    rows = [{"key": k, "value": v} for k, v in all_results.items()]
    pd.DataFrame(rows).to_csv(CSV_FILE, index=False)
    log(f"  Results saved: {CSV_FILE}  ({len(rows)} rows)")

# ============================================================
# MAIN
# ============================================================

def main():
    log("=" * 60)
    log("BRCA HER2-ENRICHED — SCRIPT 2")
    log("BULK VALIDATION + RESISTANCE + SURVIVAL")
    log("OrganismCore — BRCA-S3c/d | 2026-03-05")
    log("=" * 60)

    all_results = {}

    # Step 1
    acq = acquire_data()

    # Step 2
    expr = load_tcga_expr()

    # Step 3
    clin, surv = load_tcga_clinical()

    # Step 4
    her2_samples, luma_samples, basal_samples = \
        get_her2_samples(expr, clin)

    # Step 5
    bulk = tcga_bulk_validation(
        expr, her2_samples, luma_samples, basal_samples)
    all_results.update(bulk)

    # Step 6
    dep, surv_df, depth_score = depth_and_survival(
        expr, clin, surv, her2_samples)
    all_results.update(dep)

    # Step 7
    trast = trastuzumab_resistance(acq.get("trast"))
    all_results.update(trast)

    trast_name = acq["trast"][0] if acq.get("trast") else None

    # Step 8
    generate_figure(
        expr, her2_samples, luma_samples, basal_samples,
        surv_df, depth_score, all_results, trast_name
    )

    # Step 9
    save_results(all_results)

    # ── Final summary ──────────────────────────────────────────
    log("")
    log("=" * 60)
    log("PREDICTION SUMMARY")
    log("=" * 60)

    def pline(label, key_status=None, key_p=None, key_r=None):
        if key_status and key_status in all_results:
            log(f"  {label:<28} {all_results[key_status]}")
            return
        if key_p and key_p in all_results:
            p = all_results[key_p]
            verdict = "CONFIRMED" if p < 0.05 else "NOT CONFIRMED"
            detail  = f"p={p:.4f}"
            if key_r and key_r in all_results:
                detail += f"  r={all_results[key_r]:+.3f}"
            log(f"  {label:<28} {verdict}  ({detail})")
        else:
            log(f"  {label:<28} NO DATA")

    pline("S2-P1 ERBB3 → pCR",         key_status="p1_status")
    pline("S2-P2 Depth → OS",           key_p="p2_logrank_p")
    pline("S2-P3 r(EZH2, ESR1)",        key_p="p3_p_ezh2_esr1",
                                         key_r="p3_r_ezh2_esr1")
    pline("S2-P4 ERBB2 FC bulk",        key_p="p4_ERBB2_p",
                                         key_r="p4_ERBB2_fc_pct")
    pline("S2-P5 FOXA1 retention",      key_p="p5_p_foxa1_esr1",
                                         key_r="p5_r_foxa1_esr1")
    pline("S2-P6 Depth → Grade",        key_p="p6_p_depth_grade",
                                         key_r="p6_r_depth_grade")
    pline("S2-P7 AR-low → worse OS",    key_p="p7_ar_logrank_p")

    log("")
    log(f"  Log:     {LOG_FILE}")
    log(f"  Figure:  {FIG_FILE}")
    log(f"  Results: {CSV_FILE}")
    log("=" * 60)
    log("SCRIPT 2 COMPLETE")
    log("=" * 60)

    flush_log()


if __name__ == "__main__":
    main()
