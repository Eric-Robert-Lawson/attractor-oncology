"""
BRCA HER2-ENRICHED — SCRIPT 2
BULK RNA-seq VALIDATION + TRASTUZUMAB RESISTANCE + SURVIVAL
OrganismCore — Document BRCA-S3c/d | 2026-03-05

DATASETS:

TCGA-BRCA:
  Expression: https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2.gz
  Clinical:   https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix
  Survival:   https://pancanatlas.xenahubs.net/download/Survival_SupplementalTable_S1_20171025_xena_sp
  GDC expr:   https://gdc.xenahubs.net/download/TCGA-BRCA.htseq_fpkm.tsv.gz
  GDC surv:   https://gdc.xenahubs.net/download/TCGA-BRCA.survival.tsv
  GDC clin:   https://gdc.xenahubs.net/download/TCGA-BRCA.GDC_phenotype.tsv.gz

TRASTUZUMAB RESPONSE DATASETS (tried in priority order):
  GSE37946:   https://ftp.ncbi.nlm.nih.gov/geo/series/GSE37nnn/GSE37946/matrix/GSE37946_series_matrix.txt.gz
  GSE50948:   https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50948/matrix/GSE50948_series_matrix.txt.gz
  GSE55348:   https://ftp.ncbi.nlm.nih.gov/geo/series/GSE55nnn/GSE55348/matrix/GSE55348_series_matrix.txt.gz
  GSE66399:   https://ftp.ncbi.nlm.nih.gov/geo/series/GSE66nnn/GSE66399/matrix/GSE66399_series_matrix.txt.gz

PREDICTIONS TESTED (BRCA-S3c):
  S2-P1:  ERBB3-high predicts trastuzumab pCR
  S2-P2:  Depth score (ERBB3/CDH1/AR inverse) predicts OS
  S2-P3:  r(EZH2, ESR1) negative in HER2-enriched bulk
  S2-P4:  ERBB2 fold-change visible in bulk (>+400%), STARD3 co-elevated
  S2-P5:  FOXA1 retained in HER2 vs Basal; r(FOXA1, ESR1) positive
  S2-P6:  Depth score correlates with Grade 3
  S2-P7:  AR-low HER2+ has worse OS than AR-high HER2+
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

def download_file(url, dest, label="", timeout=120, chunk=1024*1024):
    """Download url → dest. Returns True on success."""
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        log(f"  Cached: {os.path.basename(dest)}")
        return True
    log(f"  Downloading {label or os.path.basename(dest)} ...")
    try:
        req = urllib.request.Request(
            url,
            headers={"User-Agent": "Mozilla/5.0 OrganismCore/2026"}
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
    """Try a list of URLs for the same file. Returns True if any succeed."""
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

    # TCGA expression
    log("\n-- TCGA-BRCA expression --")
    results["tcga_expr"] = try_urls(TCGA_EXPR_URLS, TCGA_EXPR_FILE,
                                    "TCGA-BRCA expression")

    # TCGA clinical
    log("\n-- TCGA-BRCA clinical --")
    results["tcga_clin"] = try_urls(TCGA_CLIN_URLS, TCGA_CLIN_FILE,
                                    "TCGA-BRCA clinical")

    # TCGA survival
    log("\n-- TCGA-BRCA survival --")
    results["tcga_surv"] = try_urls(TCGA_SURV_URLS, TCGA_SURV_FILE,
                                    "TCGA-BRCA survival")

    # Trastuzumab datasets
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
    log(f"  Sample IDs (first 3): {list(df.columns[:3])}")
    log(f"  Gene IDs (first 3):   {list(df.index[:3])}")

    # Ensure genes are rows
    if df.shape[0] < df.shape[1]:
        log("  Transposing (samples appear to be rows)...")
        df = df.T

    # Check gene coverage
    found = [g for g in ALL_GENES if g in df.index]
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

    # Clinical
    if os.path.exists(TCGA_CLIN_FILE):
        try:
            clin = pd.read_csv(TCGA_CLIN_FILE, sep="\t", index_col=0,
                               low_memory=False)
            log(f"  Clinical shape: {clin.shape}")
            log(f"  Clinical columns (first 10): {list(clin.columns[:10])}")
            # Normalise index format
            clin.index = [s.replace(".", "-") for s in clin.index]
        except Exception as e:
            log(f"  ERROR loading clinical: {e}")
    else:
        log("  Clinical file not found.")

    # Survival
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

    her2_samples = None
    luma_samples = None
    basal_samples = None
    pam50_col = None

    if clin is not None:
        # Find the PAM50 column
        pam50_candidates = [c for c in clin.columns
                            if "pam50" in c.lower() or
                            "subtype" in c.lower() or
                            "PAM50" in c]
        log(f"  PAM50 column candidates: {pam50_candidates}")

        for col in pam50_candidates:
            vals = clin[col].dropna().unique()
            log(f"  Column '{col}' values: {list(vals[:10])}")
            if any("HER2" in str(v) or "her2" in str(v).lower()
                   for v in vals):
                pam50_col = col
                break
            if any(str(v) in ["HER2-enriched", "Her2", "HER2",
                               "HER2_enriched"]
                   for v in vals):
                pam50_col = col
                break

        if pam50_col:
            log(f"  Using PAM50 column: '{pam50_col}'")
            vals = clin[pam50_col].dropna().unique()
            log(f"  All PAM50 values: {sorted([str(v) for v in vals])}")

            def match_subtype(keyword):
                mask = clin[pam50_col].astype(str).str.contains(
                    keyword, case=False, na=False)
                return list(clin.index[mask])

            her2_samples  = match_subtype("HER2")
            luma_samples  = match_subtype("LumA")
            basal_samples = match_subtype("Basal")

            log(f"  HER2-enriched: n={len(her2_samples)}")
            log(f"  LumA:          n={len(luma_samples)}")
            log(f"  Basal-like:    n={len(basal_samples)}")
        else:
            log("  WARNING: No PAM50 column found in clinical data.")

    # Fall back: use ERBB2 expression to call HER2-enriched
    if not her2_samples and expr is not None and "ERBB2" in expr.index:
        log("  Falling back to ERBB2-expression-based HER2 calling...")
        erbb2 = expr.loc["ERBB2"]
        threshold = erbb2.quantile(0.80)
        her2_samples = list(erbb2[erbb2 >= threshold].index)
        log(f"  HER2-enriched (ERBB2 top 20%): n={len(her2_samples)}")

    return her2_samples, luma_samples, basal_samples, pam50_col


# ============================================================
# STEP 5 — TCGA BULK VALIDATION
# S2-P3: r(EZH2, ESR1) negative
# S2-P4: ERBB2 fold-change visible in bulk
# S2-P5: FOXA1 retained vs Basal
# ============================================================

def tcga_bulk_validation(expr, her2_samples, luma_samples, basal_samples):
    log("")
    log("=" * 60)
    log("STEP 5: TCGA BULK VALIDATION")
    log(f"  S2-P3: r(EZH2, ESR1) in HER2-enriched")
    log(f"  S2-P4: Amplicon fold-change in bulk")
    log(f"  S2-P5: FOXA1 retention vs Basal")
    log("=" * 60)

    results = {}

    if expr is None or not her2_samples:
        log("  SKIPPED: expression or HER2 samples not available.")
        return results

    # Align sample sets to expression columns
    expr_samples = set(expr.columns)
    her2  = [s for s in her2_samples  if s in expr_samples]
    luma  = [s for s in (luma_samples  or []) if s in expr_samples]
    basal = [s for s in (basal_samples or []) if s in expr_samples]

    log(f"  HER2-enriched in expr: n={len(her2)}")
    log(f"  LumA in expr:          n={len(luma)}")
    log(f"  Basal in expr:         n={len(basal)}")

    if len(her2) < 20:
        log("  WARNING: Too few HER2 samples for reliable analysis.")
        return results

    def safe_fc(gene, a_samples, b_samples, label_a="A", label_b="B"):
        """Fold change and Mann-Whitney p for gene in group A vs B."""
        if gene not in expr.index:
            return None, None, None, None
        a_vals = expr.loc[gene, a_samples].dropna().values.astype(float)
        b_vals = expr.loc[gene, b_samples].dropna().values.astype(float)
        if len(a_vals) < 5 or len(b_vals) < 5:
            return None, None, None, None
        mean_a = np.mean(a_vals)
        mean_b = np.mean(b_vals)
        ref = max(mean_b, 0.01)
        fc_pct = (mean_a - mean_b) / ref * 100
        _, p = stats.mannwhitneyu(a_vals, b_vals, alternative="two-sided")
        return mean_a, mean_b, fc_pct, p

    # ── S2-P4: Amplicon fold-change (HER2 vs LumA) ────────────
    log("")
    log("  S2-P4 — AMPLICON FOLD-CHANGE (HER2 vs LumA)")
    if luma:
        log(f"  {'Gene':<12} {'HER2':>8} {'LumA':>8} {'FC%':>10}  p-value")
        log("  " + "-" * 55)
        for gene in AMPLICON_GENES + ["ESR1", "FOXA1", "GATA3",
                                       "EZH2", "ERBB3"]:
            ma, mb, fc, p = safe_fc(gene, her2, luma, "HER2", "LumA")
            if ma is None:
                log(f"  {gene:<12}  NOT FOUND")
                continue
            sig = "***" if p < 0.001 else "**" if p < 0.01 else \
                  "*" if p < 0.05 else "ns"
            log(f"  {gene:<12} {ma:>8.3f} {mb:>8.3f} {fc:>+10.1f}%"
                f"  p={p:.2e} {sig}")
            results[f"p4_{gene}_fc_pct"] = fc
            results[f"p4_{gene}_p"]       = p

        # Verdict
        erbb2_fc = results.get("p4_ERBB2_fc_pct", 0)
        stard3_fc = results.get("p4_STARD3_fc_pct", 0)
        log("")
        if erbb2_fc > 400:
            log(f"  S2-P4 ERBB2: CONFIRMED ({erbb2_fc:+.1f}% > +400%)")
        elif erbb2_fc > 100:
            log(f"  S2-P4 ERBB2: PARTIAL ({erbb2_fc:+.1f}%, threshold +400%)")
        else:
            log(f"  S2-P4 ERBB2: NOT CONFIRMED ({erbb2_fc:+.1f}%)")
        if stard3_fc > 100:
            log(f"  S2-P4 STARD3: CONFIRMED ({stard3_fc:+.1f}% > +100%)")
        else:
            log(f"  S2-P4 STARD3: NOT CONFIRMED ({stard3_fc:+.1f}%)")
    else:
        log("  SKIPPED: no LumA samples available.")

    # ── S2-P5: FOXA1 retention (HER2 vs Basal) ────────────────
    log("")
    log("  S2-P5 — FOXA1 RETENTION (HER2 vs Basal)")
    if basal:
        log(f"  {'Gene':<12} {'HER2':>8} {'Basal':>8} {'FC%':>10}  p-value")
        log("  " + "-" * 55)
        for gene in ["FOXA1", "ESR1", "GATA3", "KRT5", "SOX10", "EZH2"]:
            ma, mb, fc, p = safe_fc(gene, her2, basal, "HER2", "Basal")
            if ma is None:
                log(f"  {gene:<12}  NOT FOUND")
                continue
            sig = "***" if p < 0.001 else "**" if p < 0.01 else \
                  "*" if p < 0.05 else "ns"
            log(f"  {gene:<12} {ma:>8.3f} {mb:>8.3f} {fc:>+10.1f}%"
                f"  p={p:.2e} {sig}")
            results[f"p5_{gene}_her2_mean"] = ma
            results[f"p5_{gene}_basal_mean"] = mb

        # r(FOXA1, ESR1) within HER2-enriched
        if all(g in expr.index for g in ["FOXA1", "ESR1"]):
            foxa1_her2 = expr.loc["FOXA1", her2].dropna().values.astype(float)
            esr1_her2  = expr.loc["ESR1",  her2].dropna().values.astype(float)
            r_fe, p_fe = stats.pearsonr(foxa1_her2, esr1_her2)
            log(f"\n  r(FOXA1, ESR1) in HER2-enriched: r={r_fe:+.3f}  "
                f"p={p_fe:.2e}")
            results["p5_r_foxa1_esr1_her2"] = r_fe
            results["p5_p_foxa1_esr1_her2"] = p_fe
            if r_fe > 0.20 and p_fe < 0.05:
                log("  S2-P5 r(FOXA1,ESR1): CONFIRMED")
            elif r_fe > 0 and p_fe < 0.05:
                log("  S2-P5 r(FOXA1,ESR1): PARTIAL (positive but < 0.20)")
            else:
                log("  S2-P5 r(FOXA1,ESR1): NOT CONFIRMED")
    else:
        log("  SKIPPED: no Basal samples available.")

    # ── S2-P3: r(EZH2, ESR1) in HER2-enriched ────────────────
    log("")
    log("  S2-P3 — r(EZH2, ESR1) IN HER2-ENRICHED")
    if all(g in expr.index for g in ["EZH2", "ESR1"]):
        ezh2_her2 = expr.loc["EZH2", her2].dropna().values.astype(float)
        esr1_her2 = expr.loc["ESR1", her2].dropna().values.astype(float)
        r_ee, p_ee = stats.pearsonr(ezh2_her2, esr1_her2)
        log(f"  r(EZH2, ESR1) in HER2-enriched: r={r_ee:+.3f}  "
            f"p={p_ee:.2e}")
        results["p3_r_ezh2_esr1"] = r_ee
        results["p3_p_ezh2_esr1"] = p_ee
        if r_ee < -0.15 and p_ee < 0.05:
            log("  S2-P3: CONFIRMED (r < -0.15, p < 0.05)")
        elif r_ee < 0 and p_ee < 0.05:
            log("  S2-P3: PARTIAL (negative but > -0.15)")
        else:
            log("  S2-P3: NOT CONFIRMED")

        # Secondary: r(EZH2, FOXA1)
        if "FOXA1" in expr.index:
            foxa1_h = expr.loc["FOXA1", her2].dropna().values.astype(float)
            ezh2_h  = ezh2_her2[:len(foxa1_h)]
            min_n   = min(len(foxa1_h), len(ezh2_h))
            r_ef, p_ef = stats.pearsonr(ezh2_h[:min_n], foxa1_h[:min_n])
            log(f"  r(EZH2, FOXA1) in HER2-enriched: r={r_ef:+.3f}  "
                f"p={p_ef:.2e}  [secondary]")
            results["p3_r_ezh2_foxa1"] = r_ef
    else:
        log("  SKIPPED: EZH2 or ESR1 not in expression index.")

    return results


# ============================================================
# STEP 6 — DEPTH SCORE + SURVIVAL
# S2-P2: Depth score predicts OS
# S2-P6: Depth score correlates with Grade
# S2-P7: AR-low HER2+ has worse OS
# ============================================================

def depth_and_survival(expr, clin, surv, her2_samples):
    log("")
    log("=" * 60)
    log("STEP 6: DEPTH SCORE + SURVIVAL ANALYSIS")
    log(f"  S2-P2: depth score → OS")
    log(f"  S2-P6: depth score → Grade")
    log(f"  S2-P7: AR-low → worse OS")
    log("=" * 60)

    results = {}

    if expr is None or not her2_samples:
        log("  SKIPPED.")
        return results

    expr_samples = set(expr.columns)
    her2 = [s for s in her2_samples if s in expr_samples]

    if len(her2) < 20:
        log("  Too few HER2 samples.")
        return results

    # ── Compute depth score ────────────────────────────────────
    # Depth = 1 - norm(mean[ERBB3, CDH1, AR])
    # Low ERBB3 + low CDH1 + low AR → high depth → pre-resistant
    depth_avail = [g for g in DEPTH_GENES if g in expr.index]
    log(f"  Depth genes available: {depth_avail}")

    if not depth_avail:
        log("  ERROR: No depth genes found in expression data.")
        return results

    depth_vals = expr.loc[depth_avail, her2].T.mean(axis=1)

    def norm01(s):
        mn, mx = s.min(), s.max()
        if mx == mn:
            return s * 0
        return (s - mn) / (mx - mn)

    depth_score = 1.0 - norm01(depth_vals)
    depth_score.name = "depth_score"
    log(f"  Depth score range: {depth_score.min():.3f} – "
        f"{depth_score.max():.3f}")
    log(f"  Mean: {depth_score.mean():.3f}  "
        f"Median: {depth_score.median():.3f}")

    # ── Merge with survival ────────────────────────────────────
    os_col    = None
    os_event  = None
    surv_df   = None

    if surv is not None:
        # Look for OS columns in survival table
        surv_sub = surv.loc[surv.index.isin(her2)]
        os_time_candidates  = [c for c in surv.columns
                                if "OS" in c and "time" in c.lower()]
        os_event_candidates = [c for c in surv.columns
                                if "OS" in c and
                                ("event" in c.lower() or
                                 c.endswith(".1") or "status" in c.lower())]
        log(f"  OS time columns found:  {os_time_candidates}")
        log(f"  OS event columns found: {os_event_candidates}")

        if os_time_candidates and os_event_candidates:
            os_col   = os_time_candidates[0]
            os_event = os_event_candidates[0]

    if os_col is None and clin is not None:
        # Try clinical file for survival
        os_time_cands  = [c for c in clin.columns
                          if "OS" in c and
                          ("time" in c.lower() or "days" in c.lower())]
        os_ev_cands    = [c for c in clin.columns
                          if "OS" in c and
                          ("event" in c.lower() or "status" in c.lower()
                           or c.endswith(".1"))]
        if os_time_cands and os_ev_cands:
            os_col   = os_time_cands[0]
            os_event = os_ev_cands[0]
            surv     = clin

    if os_col and surv is not None:
        log(f"  Using OS columns: time='{os_col}', event='{os_event}'")
        tmp = depth_score.to_frame()
        tmp.index = [s.replace(".", "-") for s in tmp.index]
        surv.index = [s.replace(".", "-") for s in surv.index]

        merged = tmp.join(surv[[os_col, os_event]], how="inner")
        merged = merged.dropna()
        merged[os_col]   = pd.to_numeric(merged[os_col],   errors="coerce")
        merged[os_event] = pd.to_numeric(merged[os_event], errors="coerce")
        merged = merged.dropna()
        log(f"  Merged depth + survival: n={len(merged)}")

        if len(merged) >= 30:
            # Median-split log-rank
            median_d = merged["depth_score"].median()
            high_d   = merged[merged["depth_score"] >= median_d]
            low_d    = merged[merged["depth_score"] <  median_d]

            def logrank_p(t1, e1, t2, e2):
                """Simple log-rank test."""
                try:
                    from lifelines.statistics import logrank_test
                    r = logrank_test(t1, t2, event_observed_A=e1,
                                     event_observed_B=e2)
                    return r.p_value
                except ImportError:
                    pass
                # Manual Mantel-Cox
                all_t   = np.sort(np.unique(
                    np.concatenate([t1[e1 == 1], t2[e2 == 1]])))
                O1 = E1 = O2 = E2 = 0.0
                for t in all_t:
                    n1 = (t1 >= t).sum()
                    n2 = (t2 >= t).sum()
                    o1 = ((t1 == t) & (e1 == 1)).sum()
                    o2 = ((t2 == t) & (e2 == 1)).sum()
                    n  = n1 + n2
                    if n < 2:
                        continue
                    e1t = (n1 / n) * (o1 + o2)
                    e2t = (n2 / n) * (o1 + o2)
                    O1 += o1;  E1 += e1t
                    O2 += o2;  E2 += e2t
                if E1 == 0 or E2 == 0:
                    return 1.0
                chi = (O1 - E1) ** 2 / E1 + (O2 - E2) ** 2 / E2
                return 1 - chi2.cdf(chi, df=1)

            t1 = high_d[os_col].values
            e1 = high_d[os_event].values
            t2 = low_d[os_col].values
            e2 = low_d[os_event].values
            p_lr = logrank_p(t1, e1, t2, e2)

            # Median OS
            def median_os(t, e):
                idx = np.argsort(t)
                t, e = t[idx], e[idx]
                n = len(t)
                at_risk = n
                surv_est = 1.0
                for i in range(n):
                    if e[i] == 1:
                        surv_est *= (1 - 1 / at_risk)
                    if surv_est <= 0.5:
                        return t[i]
                    at_risk -= 1
                return np.nan

            med_high = median_os(t1, e1)
            med_low  = median_os(t2, e2)

            log(f"\n  S2-P2 — DEPTH SCORE VS OS (log-rank)")
            log(f"  Depth-high (n={len(high_d)}): median OS = "
                f"{med_high:.1f}" if not np.isnan(med_high) else
                f"  Depth-high (n={len(high_d)}): median OS = not reached")
            log(f"  Depth-low  (n={len(low_d)}):  median OS = "
                f"{med_low:.1f}" if not np.isnan(med_low) else
                f"  Depth-low  (n={len(low_d)}):  median OS = not reached")
            log(f"  Log-rank p = {p_lr:.4f}")

            results["p2_logrank_p"]    = p_lr
            results["p2_median_high"]  = med_high
            results["p2_median_low"]   = med_low

            if p_lr < 0.05:
                log("  S2-P2: CONFIRMED (p < 0.05)")
            else:
                log("  S2-P2: NOT CONFIRMED (p >= 0.05)")

            # Cox PH — simple univariate
            try:
                from lifelines import CoxPHFitter
                cph_data = merged[["depth_score", os_col, os_event]].copy()
                cph_data.columns = ["depth_score", "T", "E"]
                cph = CoxPHFitter()
                cph.fit(cph_data, duration_col="T", event_col="E")
                hr   = np.exp(cph.params_["depth_score"])
                p_cx = cph.summary.loc["depth_score", "p"]
                log(f"  Cox PH: HR={hr:.3f}  p={p_cx:.4f}")
                results["p2_cox_hr"] = hr
                results["p2_cox_p"]  = p_cx
            except Exception:
                log("  Cox PH: lifelines not available — skipped.")

            surv_df = merged.copy()
            surv_df["depth_group"] = np.where(
                surv_df["depth_score"] >= median_d, "high", "low")

        else:
            log(f"  WARNING: only {len(merged)} samples — "
                "survival test unreliable.")
    else:
        log("  WARNING: OS columns not found — survival test skipped.")

    # ── S2-P7: AR-low vs AR-high OS ───────────────────────────
    log("")
    log("  S2-P7 — AR-LOW vs AR-HIGH WITHIN HER2-ENRICHED")

    if "AR" in expr.index and surv is not None and os_col:
        ar_vals = expr.loc["AR", her2].dropna()
        ar_median = ar_vals.median()
        ar_high   = list(ar_vals[ar_vals >= ar_median].index)
        ar_low    = list(ar_vals[ar_vals <  ar_median].index)
        log(f"  AR-high: n={len(ar_high)}  AR-low: n={len(ar_low)}")

        tmp = ar_vals.to_frame("AR")
        tmp.index = [s.replace(".", "-") for s in tmp.index]
        surv.index = [s.replace(".", "-") for s in surv.index]
        ar_merged = tmp.join(surv[[os_col, os_event]], how="inner").dropna()
        ar_merged[os_col]   = pd.to_numeric(ar_merged[os_col],   errors="coerce")
        ar_merged[os_event] = pd.to_numeric(ar_merged[os_event], errors="coerce")
        ar_merged = ar_merged.dropna()

        if len(ar_merged) >= 30:
            ar_med_split = ar_merged["AR"].median()
            ar_h = ar_merged[ar_merged["AR"] >= ar_med_split]
            ar_l = ar_merged[ar_merged["AR"] <  ar_med_split]
            p_ar = logrank_p(
                ar_h[os_col].values, ar_h[os_event].values,
                ar_l[os_col].values, ar_l[os_event].values
            )
            log(f"  AR-high n={len(ar_h)}, AR-low n={len(ar_l)}")
            log(f"  Log-rank p = {p_ar:.4f}")
            results["p7_ar_logrank_p"] = p_ar
            if p_ar < 0.05:
                # Check direction
                med_ar_h = median_os(ar_h[os_col].values,
                                     ar_h[os_event].values)
                med_ar_l = median_os(ar_l[os_col].values,
                                     ar_l[os_event].values)
                log(f"  Median OS AR-high: {med_ar_h:.1f}")
                log(f"  Median OS AR-low:  {med_ar_l:.1f}")
                if (not np.isnan(med_ar_h) and not np.isnan(med_ar_l)
                        and med_ar_l < med_ar_h):
                    log("  S2-P7: CONFIRMED (AR-low → worse OS, p < 0.05)")
                else:
                    log("  S2-P7: SIGNIFICANT but direction check failed")
            else:
                log("  S2-P7: NOT CONFIRMED (p >= 0.05)")
        else:
            log("  WARNING: too few samples for AR survival test.")
    else:
        log("  SKIPPED: AR not in expression or no survival data.")

    # ── S2-P6: Depth score vs Grade ───────────────────────────
    log("")
    log("  S2-P6 — DEPTH SCORE vs HISTOLOGICAL GRADE")

    grade_col = None
    if clin is not None:
        grade_cands = [c for c in clin.columns
                       if "grade" in c.lower() or "Grade" in c]
        log(f"  Grade column candidates: {grade_cands}")
        if grade_cands:
            grade_col = grade_cands[0]

    if grade_col and clin is not None:
        tmp_d = depth_score.to_frame()
        tmp_d.index = [s.replace(".", "-") for s in tmp_d.index]
        clin.index  = [s.replace(".", "-") for s in clin.index]
        gm = tmp_d.join(clin[[grade_col]], how="inner").dropna()
        gm[grade_col] = pd.to_numeric(gm[grade_col], errors="coerce")
        gm = gm.dropna()
        log(f"  Grade-depth merged: n={len(gm)}")
        if len(gm) >= 20:
            r_gd, p_gd = stats.spearmanr(gm["depth_score"], gm[grade_col])
            log(f"  Spearman r(depth, grade) = {r_gd:+.3f}  p={p_gd:.4f}")
            results["p6_r_depth_grade"] = r_gd
            results["p6_p_depth_grade"] = p_gd
            if r_gd > 0 and p_gd < 0.05:
                log("  S2-P6: CONFIRMED (higher depth → higher grade)")
            else:
                log("  S2-P6: NOT CONFIRMED")
        else:
            log("  Too few samples for grade analysis.")
    else:
        log("  Grade column not found — S2-P6 skipped.")

    return results, surv_df, depth_score


# ============================================================
# STEP 7 — TRASTUZUMAB RESISTANCE
# S2-P1: ERBB3-high → pCR
# ============================================================

def trastuzumab_resistance(trast_result):
    log("")
    log("=" * 60)
    log("STEP 7: TRASTUZUMAB RESISTANCE ANALYSIS")
    log(f"  S2-P1: ERBB3 as pCR biomarker")
    log("=" * 60)

    results = {}

    if trast_result is None:
        log("  SKIPPED: No trastuzumab dataset available.")
        results["p1_status"] = "SKIPPED"
        return results

    name, fpath = trast_result
    log(f"  Dataset: {name}  File: {os.path.basename(fpath)}")

    # ── Parse series matrix ────────────────────────────────────
    meta_rows = {}
    probe_ids = []
    expr_data = {}
    sample_ids = []
    expr_started = False
    n_header = 0

    try:
        with gzip.open(fpath, "rt", errors="ignore") as f:
            for i, line in enumerate(f):
                line = line.rstrip("\n")

                if line.startswith("!Sample_geo_accession"):
                    parts = line.split("\t")
                    sample_ids = [p.strip().strip('"')
                                  for p in parts[1:]]

                if line.startswith("!Sample_characteristics_ch1"):
                    parts = line.split("\t")
                    val0  = parts[1].strip().strip('"') if len(parts) > 1 else ""
                    key   = (val0.split(":")[0].strip()
                             if ":" in val0 else val0[:40])
                    row   = [p.strip().strip('"') for p in parts[1:]]
                    if key not in meta_rows:
                        meta_rows[key] = row

                if "series_matrix_table_begin" in line.lower():
                    expr_started = True
                    n_header = i
                    continue

                if "series_matrix_table_end" in line.lower():
                    break

                if expr_started and i > n_header:
                    parts = line.split("\t")
                    if not parts[0].strip().strip('"'):
                        continue
                    pid = parts[0].strip().strip('"')
                    if pid == "ID_REF":
                        if not sample_ids:
                            sample_ids = [p.strip().strip('"')
                                          for p in parts[1:]]
                        continue
                    if len(parts) > 1:
                        try:
                            vals = [float(p) if p.strip() not in ("", "null",
                                    "NA", "NaN") else np.nan
                                    for p in parts[1:]]
                            expr_data[pid] = vals
                        except Exception:
                            pass

    except Exception as e:
        log(f"  ERROR parsing {name}: {e}")
        results["p1_status"] = f"PARSE_ERROR: {e}"
        return results

    log(f"  Samples parsed: {len(sample_ids)}")
    log(f"  Meta fields:    {list(meta_rows.keys())[:10]}")
    log(f"  Probes parsed:  {len(expr_data)}")

    if not expr_data or not sample_ids:
        log("  ERROR: No expression data parsed.")
        results["p1_status"] = "NO_DATA"
        return results

    # ── Build expression DataFrame ─────────────────────────────
    n_samp = len(sample_ids)
    clean  = {k: v for k, v in expr_data.items() if len(v) == n_samp}
    expr_m = pd.DataFrame(clean, index=sample_ids).T   # genes × samples

    log(f"  Expression matrix: {expr_m.shape}")

    # ── Identify pCR column ────────────────────────────────────
    pcr_col   = None
    pcr_vals  = None
    for key, row in meta_rows.items():
        vals_lower = [str(v).lower() for v in row]
        if any("pcr" in v or "pathologic" in v or
               "response" in v or "pCR" in v or
               "rCR" in v.lower()
               for v in vals_lower):
            pcr_col  = key
            pcr_vals = row
            break

    # Generic response search
    if pcr_col is None:
        for key, row in meta_rows.items():
            vals_lower = [str(v).lower() for v in row]
            if any("response" in v or "resistant" in v or
                   "sensitive" in v for v in vals_lower):
                pcr_col  = key
                pcr_vals = row
                break

    if pcr_col is None:
        log(f"  WARNING: pCR/response column not found in {name}.")
        log(f"  Available meta fields: {list(meta_rows.keys())}")
        results["p1_status"] = "NO_PCR_COLUMN"
        return results

    log(f"  pCR column: '{pcr_col}'")
    log(f"  pCR values (unique): "
        f"{list(set(str(v) for v in pcr_vals))[:10]}")

    # Build metadata DataFrame
    meta_df = pd.DataFrame(meta_rows, index=sample_ids)
    meta_df.columns = [c.strip() for c in meta_df.columns]

    # ── Map pCR labels → 0/1 ──────────────────────────────────
    def parse_pcr(v):
        v = str(v).lower().strip()
        if any(x in v for x in ["pcr", "complete", "cr ", "cr,",
                                 "sensitive", "response: yes",
                                 "pCR".lower()]):
            return 1
        if any(x in v for x in ["rd", "residual", "non", "no response",
                                 "resistant", "stable", "progressive"]):
            return 0
        return np.nan

    meta_df["pcr_binary"] = meta_df[pcr_col].apply(parse_pcr)
    pcr_known = meta_df.dropna(subset=["pcr_binary"])
    log(f"  pCR=1: n={int(pcr_known['pcr_binary'].sum())}  "
        f"pCR=0: n={int((pcr_known['pcr_binary'] == 0).sum())}")

    if len(pcr_known) < 10:
        log("  WARNING: too few labelled samples.")
        results["p1_status"] = "TOO_FEW"
        return results

    # ── Identify ERBB3 probes ──────────────────────────────────
    # For Affymetrix arrays, probe IDs matching ERBB3
    # Common HG-U133A/Plus2 ERBB3 probes:
    ERBB3_PROBES = [
        "205047_s_at",   # HG-U133A classic ERBB3
        "205048_s_at",
        "210766_s_at",
        "1565483_at",    # HG-U133Plus2
        "1565484_x_at",
    ]
    # Also try: probes whose IDs contain "ERBB3" (some normalised datasets)
    erbb3_found = [p for p in ERBB3_PROBES if p in expr_m.index]
    erbb3_sym   = [p for p in expr_m.index
                   if str(p).upper() == "ERBB3"]
    erbb3_probes = erbb3_found + erbb3_sym

    if not erbb3_probes:
        log("  ERBB3 probes not found — trying gene-symbol index...")
        # Some processed series matrices use gene symbols as row index
        if "ERBB3" in expr_m.index:
            erbb3_probes = ["ERBB3"]

    if not erbb3_probes:
        log(f"  WARNING: ERBB3 not found in {name}.")
        log(f"  First 20 probe IDs: {list(expr_m.index[:20])}")
        results["p1_status"] = "ERBB3_NOT_FOUND"
        return results

    log(f"  ERBB3 probes found: {erbb3_probes}")

    # Mean across probes
    erbb3_expr = expr_m.loc[erbb3_probes].mean(axis=0)
    erbb3_expr.name = "ERBB3"

    # Merge with pCR
    merged = erbb3_expr.to_frame().join(
        meta_df[["pcr_binary"]], how="inner"
    ).dropna()
    log(f"  ERBB3 + pCR merged: n={len(merged)}")

    if len(merged) < 10:
        log("  WARNING: too few samples after merge.")
        results["p1_status"] = "TOO_FEW_MERGED"
        return results

    # Mann-Whitney
    pcr1   = merged[merged["pcr_binary"] == 1]["ERBB3"].values
    pcr0   = merged[merged["pcr_binary"] == 0]["ERBB3"].values
    log(f"  ERBB3 in pCR  (n={len(pcr1)}): "
        f"mean={np.mean(pcr1):.3f}  median={np.median(pcr1):.3f}")
    log(f"  ERBB3 in RD   (n={len(pcr0)}): "
        f"mean={np.mean(pcr0):.3f}  median={np.median(pcr0):.3f}")

    if len(pcr1) >= 3 and len(pcr0) >= 3:
        _, p_mw = stats.mannwhitneyu(pcr1, pcr0, alternative="two-sided")
        log(f"  Mann-Whitney p = {p_mw:.4f}")
        results["p1_p_mw"] = p_mw
        results["p1_erbb3_pcr_mean"] = float(np.mean(pcr1))
        results["p1_erbb3_rd_mean"]  = float(np.mean(pcr0))

        direction = "ERBB3-high → pCR" if np.mean(pcr1) > np.mean(pcr0) \
                    else "ERBB3-low → pCR (OPPOSITE)"

        if np.mean(pcr1) > np.mean(pcr0) and p_mw < 0.05:
            log(f"  S2-P1: CONFIRMED ({direction}, p={p_mw:.4f})")
            results["p1_status"] = "CONFIRMED"
        elif p_mw < 0.05:
            log(f"  S2-P1: SIGNIFICANT but WRONG DIRECTION ({direction})")
            results["p1_status"] = "WRONG_DIRECTION"
        else:
            log(f"  S2-P1: NOT CONFIRMED (p={p_mw:.4f}, {direction})")
            results["p1_status"] = "NOT_CONFIRMED"

        # AUC
        if SKLEARN_OK and len(merged) >= 10:
            try:
                auc = roc_auc_score(merged["pcr_binary"],
                                    merged["ERBB3"])
                log(f"  ERBB3 AUC = {auc:.3f}")
                results["p1_auc"] = auc
            except Exception:
                pass
    else:
        log("  Too few samples in one group.")
        results["p1_status"] = "TOO_FEW_GROUP"

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

        expr_samples = set(expr.columns) if expr is not None else set()
        her2  = [s for s in (her2_samples  or []) if s in expr_samples]
        luma  = [s for s in (luma_samples  or []) if s in expr_samples]
        basal = [s for s in (basal_samples or []) if s in expr_samples]

        # ── Panel A: Amplicon genes HER2 vs LumA ──────────────
        ax_a = fig.add_subplot(gs[0, 0])
        amp_genes = [g for g in AMPLICON_GENES if
                     expr is not None and g in expr.index]
        if amp_genes and luma:
            her2_means  = [expr.loc[g, her2].mean() for g in amp_genes]
            luma_means  = [expr.loc[g, luma].mean() for g in amp_genes]
            x = np.arange(len(amp_genes))
            w = 0.35
            ax_a.bar(x - w/2, her2_means, w,
                     label="HER2-enriched", color="#c0392b", alpha=0.85)
            ax_a.bar(x + w/2, luma_means,  w,
                     label="LumA",          color="#2980b9", alpha=0.85)
            ax_a.set_xticks(x)
            ax_a.set_xticklabels(amp_genes, fontsize=9)
            ax_a.set_ylabel("Mean expression (log2)")
            ax_a.set_title("A — Amplicon genes: HER2 vs LumA\n(S2-P4)",
                           fontsize=10)
            ax_a.legend(fontsize=8)
        else:
            ax_a.text(0.5, 0.5, "Data not available",
                      ha="center", va="center", transform=ax_a.transAxes)
            ax_a.set_title("A — Amplicon genes")

        # ── Panel B: FOXA1/KRT5/SOX10 HER2 vs Basal ──────────
        ax_b = fig.add_subplot(gs[0, 1])
        comp_genes = [g for g in ["FOXA1", "ESR1", "KRT5", "SOX10", "EZH2"]
                      if expr is not None and g in expr.index]
        if comp_genes and basal:
            her2_m  = [expr.loc[g, her2].mean()  for g in comp_genes]
            basal_m = [expr.loc[g, basal].mean() for g in comp_genes]
            x = np.arange(len(comp_genes))
            w = 0.35
            ax_b.bar(x - w/2, her2_m,  w,
                     label="HER2-enriched", color="#c0392b", alpha=0.85)
            ax_b.bar(x + w/2, basal_m, w,
                     label="Basal-like",    color="#8e44ad", alpha=0.85)
            ax_b.set_xticks(x)
            ax_b.set_xticklabels(comp_genes, fontsize=9)
            ax_b.set_ylabel("Mean expression (log2)")
            ax_b.set_title("B — HER2 vs Basal: luminal retention\n(S2-P5)",
                           fontsize=10)
            ax_b.legend(fontsize=8)
        else:
            ax_b.text(0.5, 0.5, "Data not available",
                      ha="center", va="center", transform=ax_b.transAxes)
            ax_b.set_title("B — HER2 vs Basal")

        # ── Panel C: r(EZH2, ESR1) scatter ────────────────────
        ax_c = fig.add_subplot(gs[0, 2])
        if (expr is not None and "EZH2" in expr.index and
                "ESR1" in expr.index and her2):
            ez = expr.loc["EZH2", her2].dropna().values.astype(float)
            es = expr.loc["ESR1", her2].dropna().values.astype(float)
            n  = min(len(ez), len(es))
            ax_c.scatter(ez[:n], es[:n], alpha=0.25, s=8,
                         color="#c0392b", rasterized=True)
            r_val = results.get("p3_r_ezh2_esr1", np.nan)
            p_val = results.get("p3_p_ezh2_esr1", np.nan)
            if not np.isnan(r_val):
                m, b = np.polyfit(ez[:n], es[:n], 1)
                xr = np.array([ez[:n].min(), ez[:n].max()])
                ax_c.plot(xr, m * xr + b, color="black", lw=1.5)
                ax_c.set_title(
                    f"C — r(EZH2, ESR1) in HER2-enriched\n"
                    f"r={r_val:+.3f}  p={p_val:.2e}  (S2-P3)",
                    fontsize=10)
            else:
                ax_c.set_title("C — r(EZH2, ESR1)\n(S2-P3)", fontsize=10)
            ax_c.set_xlabel("EZH2 expression")
            ax_c.set_ylabel("ESR1 expression")
        else:
            ax_c.text(0.5, 0.5, "Data not available",
                      ha="center", va="center", transform=ax_c.transAxes)
            ax_c.set_title("C — r(EZH2, ESR1)")

        # ── Panel D: Depth score distribution ─────────────────
        ax_d = fig.add_subplot(gs[1, 0])
        if depth_score is not None and len(depth_score) > 0:
            ax_d.hist(depth_score.values, bins=40,
                      color="#c0392b", alpha=0.75, edgecolor="white")
            ax_d.axvline(depth_score.median(), color="black",
                         lw=1.5, ls="--", label="Median")
            ax_d.set_xlabel("Depth score (0=shallow, 1=deep)")
            ax_d.set_ylabel("Number of samples")
            ax_d.set_title(
                "D — Depth score distribution\n"
                "(ERBB3/CDH1/AR inverse, S2-P2)",
                fontsize=10)
            ax_d.legend(fontsize=8)
        else:
            ax_d.text(0.5, 0.5, "Data not available",
                      ha="center", va="center", transform=ax_d.transAxes)
            ax_d.set_title("D — Depth score")

        # ── Panel E: Survival curves ───────────────────────────
        ax_e = fig.add_subplot(gs[1, 1])
        if surv_df is not None and "depth_group" in surv_df.columns:
            os_cols = [c for c in surv_df.columns
                       if "OS" in c and "time" in c.lower()]
            os_ev   = [c for c in surv_df.columns
                       if "OS" in c and
                       ("event" in c.lower() or c.endswith(".1"))]
            if os_cols and os_ev:
                tc = os_cols[0]
                ec = os_ev[0]
                colors = {"high": "#c0392b", "low": "#2980b9"}
                for grp, grp_df in surv_df.groupby("depth_group"):
                    t = grp_df[tc].values
                    e = grp_df[ec].values
                    idx = np.argsort(t)
                    t, e = t[idx], e[idx]
                    n = len(t)
                    s = 1.0
                    times  = [0]
                    survs  = [1.0]
                    at_risk = n
                    for i in range(n):
                        if e[i] == 1:
                            s *= (1 - 1 / at_risk)
                        times.append(t[i])
                        survs.append(s)
                        at_risk -= 1
                    ax_e.step(times, survs, where="post",
                              color=colors.get(grp, "grey"),
                              label=f"Depth-{grp} (n={len(grp_df)})",
                              lw=2)
                p_val = results.get("p2_logrank_p", np.nan)
                ax_e.set_title(
                    f"E — OS by depth score (HER2-enriched)\n"
                    f"log-rank p={p_val:.4f}  (S2-P2)",
                    fontsize=10)
                ax_e.set_xlabel("Time")
                ax_e.set_ylabel("Survival probability")
                ax_e.legend(fontsize=8)
                ax_e.set_ylim(0, 1.05)
            else:
                ax_e.text(0.5, 0.5, "OS columns not found",
                          ha="center", va="center",
                          transform=ax_e.transAxes)
        else:
            ax_e.text(0.5, 0.5, "Survival data not available",
                      ha="center", va="center", transform=ax_e.transAxes)
        ax_e.set_title(
            ax_e.get_title() if ax_e.get_title() else
            "E — Survival by depth score", fontsize=10)

        # ── Panel F: AR expression in HER2 ────────────────────
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
            ax_f.set_xlabel("AR expression (log2)")
            ax_f.set_ylabel("Number of samples")
            ax_f.legend(fontsize=8)
        else:
            ax_f.text(0.5, 0.5, "AR data not available",
                      ha="center", va="center", transform=ax_f.transAxes)
            ax_f.set_title("F — AR in HER2-enriched")

        # ── Panel G: Prediction scorecard ─────────────────────
        ax_g = fig.add_subplot(gs[2, 0])
        ax_g.axis("off")

        def verdict(key_status=None, key_p=None, key_r=None,
                    direction=None, threshold_p=0.05,
                    r_threshold=None, r_direction=None):
            if key_status and key_status in results:
                s = results[key_status]
                return s
            if key_p and key_p in results:
                p = results[key_p]
                if p < threshold_p:
                    if key_r and key_r in results and r_threshold:
                        r = results[key_r]
                        if r_direction == "negative" and r < -r_threshold:
                            return "CONFIRMED"
                        if r_direction == "positive" and r > r_threshold:
                            return "CONFIRMED"
                        return "PARTIAL"
                    return "CONFIRMED"
                return "NOT CONFIRMED"
            return "NO DATA"

        scorecard = [
            ("S2-P1",
             f"ERBB3→pCR: {results.get('p1_status','NO DATA')}"),
            ("S2-P2",
             f"Depth→OS:  p={results.get('p2_logrank_p', np.nan):.4f}"
             if "p2_logrank_p" in results else "Depth→OS: NO DATA"),
            ("S2-P3",
             f"r(EZH2,ESR1)={results.get('p3_r_ezh2_esr1', np.nan):.3f}"
             if "p3_r_ezh2_esr1" in results else "r(EZH2,ESR1): NO DATA"),
            ("S2-P4",
             f"ERBB2 FC={results.get('p4_ERBB2_fc_pct', np.nan):.1f}%"
             if "p4_ERBB2_fc_pct" in results else "ERBB2 FC: NO DATA"),
            ("S2-P5",
             f"r(FOXA1,ESR1)={results.get('p5_r_foxa1_esr1_her2', np.nan):.3f}"
             if "p5_r_foxa1_esr1_her2" in results else "r(FOXA1,ESR1): NO DATA"),
            ("S2-P6",
             f"r(depth,grade)={results.get('p6_r_depth_grade', np.nan):.3f}"
             if "p6_r_depth_grade" in results else "r(depth,grade): NO DATA"),
            ("S2-P7",
             f"AR OS p={results.get('p7_ar_logrank_p', np.nan):.4f}"
             if "p7_ar_logrank_p" in results else "AR OS: NO DATA"),
        ]
        txt = "PREDICTION SCORECARD\n" + "─" * 34 + "\n"
        for pred, detail in scorecard:
            txt += f"{pred}:  {detail}\n"
        ax_g.text(0.05, 0.95, txt, transform=ax_g.transAxes,
                  fontsize=9, va="top", ha="left",
                  fontfamily="monospace",
                  bbox=dict(fc="lightyellow", ec="grey", alpha=0.8))

        # ── Panel H: Depth correlations bar chart ──────────────
        ax_h = fig.add_subplot(gs[2, 1])
        if expr is not None and depth_score is not None and her2:
            corr_genes = (DEPTH_GENES + AMPLICON_GENES + EPIGENETIC_GENES
                          + ["FOXA1", "ESR1", "MKI67"])
            corr_genes = list(dict.fromkeys(
                [g for g in corr_genes if g in expr.index]))
            rs, gs_names = [], []
            for g in corr_genes:
                gvals = expr.loc[g, her2].dropna().values.astype(float)
                common = depth_score.index.intersection(
                    pd.Index(her2))
                dvals = depth_score.loc[common].values.astype(float)
                n = min(len(gvals), len(dvals))
                if n >= 20:
                    r, _ = stats.pearsonr(dvals[:n], gvals[:n])
                    rs.append(r)
                    gs_names.append(g)
            if rs:
                idx = np.argsort(rs)
                rs_sorted = [rs[i] for i in idx]
                gn_sorted = [gs_names[i] for i in idx]
                colors_bar = ["#c0392b" if r < 0 else "#27ae60"
                              for r in rs_sorted]
                ax_h.barh(range(len(rs_sorted)), rs_sorted,
                          color=colors_bar, alpha=0.8)
                ax_h.set_yticks(range(len(gn_sorted)))
                ax_h.set_yticklabels(gn_sorted, fontsize=7)
                ax_h.axvline(0, color="black", lw=0.8)
                ax_h.set_xlabel("Pearson r with depth score")
                ax_h.set_title(
                    "H — Depth score correlations\n(HER2-enriched bulk)",
                    fontsize=10)
        else:
            ax_h.text(0.5, 0.5, "Data not available",
                      ha="center", va="center", transform=ax_h.transAxes)
            ax_h.set_title("H — Depth correlations")

        # ── Panel I: EZH2 context ──────────────────────────────
        ax_i = fig.add_subplot(gs[2, 2])
        ax_i.axis("off")
        ezh2_fc    = results.get("p4_EZH2_fc_pct",   "N/A")
        r_ee       = results.get("p3_r_ezh2_esr1",   "N/A")
        r_ef       = results.get("p3_r_ezh2_foxa1",  "N/A")
        r_fe       = results.get("p5_r_foxa1_esr1_her2", "N/A")
        trast_txt  = f"Dataset: {trast_name}" if trast_name else "No dataset"

        def fmt(v):
            return f"{v:+.3f}" if isinstance(v, float) else str(v)

        ctx_txt = (
            f"EZH2 AXIS — HER2-ENRICHED\n"
            f"{'─'*32}\n"
            f"EZH2 FC vs LumA:     {fmt(ezh2_fc)}%\n"
            f"r(EZH2, ESR1):       {fmt(r_ee)}\n"
            f"r(EZH2, FOXA1):      {fmt(r_ef)}\n"
            f"r(FOXA1, ESR1):      {fmt(r_fe)}\n"
            f"\n"
            f"TRASTUZUMAB DATASET\n"
            f"{'─'*32}\n"
            f"{trast_txt}\n"
            f"ERBB3 pCR mean: "
            f"{fmt(results.get('p1_erbb3_pcr_mean', 'N/A'))}\n"
            f"ERBB3 RD mean:  "
            f"{fmt(results.get('p1_erbb3_rd_mean', 'N/A'))}\n"
            f"AUC:            "
            f"{fmt(results.get('p1_auc', 'N/A'))}\n"
        )
        ax_i.text(0.05, 0.95, ctx_txt, transform=ax_i.transAxes,
                  fontsize=9, va="top", ha="left",
                  fontfamily="monospace",
                  bbox=dict(fc="lightcyan", ec="grey", alpha=0.8))
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
# STEP 9 — SAVE RESULTS CSV
# ============================================================

def save_results(all_results):
    log("")
    log("=" * 60)
    log("STEP 9: SAVE RESULTS")
    log("=" * 60)

    rows = []
    for k, v in all_results.items():
        rows.append({"key": k, "value": v})
    df = pd.DataFrame(rows)
    df.to_csv(CSV_FILE, index=False)
    log(f"  Results saved: {CSV_FILE}")
    log(f"  Rows: {len(df)}")


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

    # Step 1 — Acquire data
    acq = acquire_data()

    # Step 2 — Load TCGA expression
    expr = load_tcga_expr()

    # Step 3 — Load clinical + survival
    clin, surv = load_tcga_clinical()

    # Step 4 — Identify subtype samples
    her2_samples, luma_samples, basal_samples, pam50_col = \
        get_her2_samples(expr, clin)

    # Step 5 — Bulk validation (P3, P4, P5)
    bulk_results = tcga_bulk_validation(
        expr, her2_samples, luma_samples, basal_samples)
    all_results.update(bulk_results)

    # Step 6 — Depth + survival (P2, P6, P7)
    depth_results, surv_df, depth_score = depth_and_survival(
        expr, clin, surv, her2_samples)
    all_results.update(depth_results)

    # Step 7 — Trastuzumab resistance (P1)
    trast_results = trastuzumab_resistance(acq.get("trast"))
    all_results.update(trast_results)

    trast_name = acq["trast"][0] if acq.get("trast") else None

    # Step 8 — Figure
    generate_figure(
        expr, her2_samples, luma_samples, basal_samples,
        surv_df, depth_score, all_results, trast_name
    )

    # Step 9 — Save results
    save_results(all_results)

    # ── Final summary ──────────────────────────────────────────
    log("")
    log("=" * 60)
    log("PREDICTION SUMMARY")
    log("=" * 60)

    def p(label, key_status=None, key_p=None, key_r=None,
          threshold_p=0.05, expected_direction=None):
        if key_status and key_status in all_results:
            log(f"  {label}: {all_results[key_status]}")
            return
        if key_p in all_results:
            p_val = all_results[key_p]
            detail = f"p={p_val:.4f}"
            if key_r and key_r in all_results:
                detail += f"  r={all_results[key_r]:+.3f}"
            verdict = "CONFIRMED" if p_val < threshold_p else "NOT CONFIRMED"
            log(f"  {label}: {verdict}  ({detail})")
        else:
            log(f"  {label}: NO DATA")

    p("S2-P1 ERBB3→pCR",        key_status="p1_status")
    p("S2-P2 Depth→OS",         key_p="p2_logrank_p")
    p("S2-P3 r(EZH2,ESR1)",     key_p="p3_p_ezh2_esr1",
                                 key_r="p3_r_ezh2_esr1")
    p("S2-P4 ERBB2 FC bulk",    key_p="p4_ERBB2_p",
                                 key_r="p4_ERBB2_fc_pct")
    p("S2-P5 FOXA1 retention",  key_p="p5_p_foxa1_esr1_her2",
                                 key_r="p5_r_foxa1_esr1_her2")
    p("S2-P6 Depth→Grade",      key_p="p6_p_depth_grade",
                                 key_r="p6_r_depth_grade")
    p("S2-P7 AR-low→worse OS",  key_p="p7_ar_logrank_p")

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
