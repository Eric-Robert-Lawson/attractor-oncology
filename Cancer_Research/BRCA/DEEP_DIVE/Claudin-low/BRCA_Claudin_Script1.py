"""
BRCA CLAUDIN-LOW — SCRIPT 1
OrganismCore — Document BRCA-S7a/b | 2026-03-05

GEOMETRY-FIRST (Protocol v2.0):
  This script is a DISCOVERY script.
  It reveals what the data contains before testing predictions.
  The unfiltered top-mover scan is printed FIRST.
  Prediction panel tests are printed SECOND.

CLAUDIN-LOW IDENTIFICATION — FIRST PRINCIPLES:
  Claudin-low identity is derived from the 10-gene
  geometry-first signature. No external classifier label
  is imported. Samples are identified from TCGA-BRCA
  expression geometry:

    NEGATIVE markers (must be LOW vs. cohort median):
      CLDN3, CLDN4, CLDN7   — claudin junction proteins
      CDH1                   — E-cadherin
      ESR1                   — luminal identity TF

    POSITIVE markers (must be HIGH vs. cohort median):
      VIM                    — mesenchymal/stem marker
      CD44                   — breast stem cell marker
      SNAI1                  — EMT/stem TF
      ZEB1                   — EMT TF
      FN1                    — fibronectin (mesenchymal ECM)

  Classification rule: sample scores +1 for each positive
  marker above its cohort median, -1 for each negative
  marker below its cohort median. Max score = +10.
  Claudin-low = score >= 7 (at least 7 of 10 criteria met).
  Threshold is logged and tested at 6, 7, and 8 for sensitivity.

  ERBB2 check: any sample with ERBB2 expression in the
  top 10% of the cohort is excluded from claudin-low
  (claudin-low is ERBB2-non-amplified by definition).

  PAM50 cross-check: claudin-low samples are expected to
  call as Basal or Normal-like in PAM50. This is logged
  but not used as a filter — the geometry leads.

SELF-CONTAINED:
  Downloads TCGA-BRCA bulk RNA-seq (HiSeqV2) and clinical
  matrix to Claudin_Low_s1_analysis/data/.
  Does NOT require or reuse data from any prior subtype folder.
  All confirmed working URLs from prior scripts in this series.

CROSS-SUBTYPE COMPARISON:
  Claudin-low compared against all five completed subtypes
  (LumA, LumB, HER2-enriched, TNBC/Basal, ILC) and normal
  using the same TCGA-BRCA dataset.
  PAM50Call_RNAseq used for subtype classification
  (confirmed present in this dataset — see ILC s2 log).

PREDICTIONS TESTED (from BRCA-S7a — predictions.md):
  P1: ESR1/FOXA1/GATA3 lowest of all BRCA subtypes
  P2: VIM/ZEB1/SNAI1 elevated — partial EMT programme
  P3: CLDN3/CLDN4/CLDN7 lowest of all BRCA subtypes
  P4: CD44 high / CD24 low — stem phenotype
  P5: TP53 target signature less disrupted than TNBC
  P6: Immune programme present — PD-L1/CD8A elevated
  P7: Depth axis = stem programme depth, not ER axis
  P8: Most distant from LumA in cross-subtype PCA
  P9: Drug targets: immune + stem, NOT ER or HER2
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

# ============================================================
# CONFIGURATION — SELF-CONTAINED
# ============================================================

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
BASE_DIR    = os.path.join(SCRIPT_DIR, "Claudin_Low_s1_analysis")
DATA_DIR    = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")

for d in [BASE_DIR, DATA_DIR, RESULTS_DIR]:
    os.makedirs(d, exist_ok=True)

LOG_FILE        = os.path.join(RESULTS_DIR, "cl_s1_log.txt")
FIG_FILE        = os.path.join(RESULTS_DIR, "cl_s1_figure.png")
CSV_FILE        = os.path.join(RESULTS_DIR, "cl_s1_panel.csv")
TOP_MOVERS_FILE = os.path.join(RESULTS_DIR, "cl_s1_top_movers.csv")
CROSS_FILE      = os.path.join(RESULTS_DIR, "cl_s1_cross_subtype.csv")
DEPTH_FILE      = os.path.join(RESULTS_DIR, "cl_s1_depth_axis.csv")

# ── Confirmed working URLs (validated across prior BRCA scripts) ──
TCGA_EXPR_URLS = [
    "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2.gz",
    "https://pancanatlas.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2.gz",
    "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2_PANCAN.gz",
]
TCGA_CLIN_URLS = [
    "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix",
    "https://pancanatlas.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix",
]
TCGA_EXPR_FILE = os.path.join(DATA_DIR, "TCGA_BRCA_HiSeqV2.gz")
TCGA_CLIN_FILE = os.path.join(DATA_DIR, "TCGA_BRCA_clinicalMatrix.tsv")

# ── TCGA barcode constants (confirmed from ILC s2 log) ──────────
PATIENT_PREFIX_LEN = 12
SAMPLE_TYPE_POS    = (13, 15)  # 0-indexed; "01"=tumour, "11"=normal

# ── Claudin-low classification thresholds ───────────────────────
CL_SCORE_THRESHOLDS = [6, 7, 8]   # tested; 7 is primary
CL_PRIMARY_THRESHOLD = 7
ERBB2_EXCLUSION_PERCENTILE = 90   # top 10% ERBB2 → exclude

# ── Minimum sample sizes ─────────────────────────────────────────
MIN_CL_SAMPLES  = 10   # claudin-low is rare — accept smaller n
MIN_NORMAL      = 30

# ============================================================
# THE 10-GENE CLAUDIN-LOW GEOMETRY SIGNATURE
# Derived from first principles per BRCA-S7a
# ============================================================

# Negative markers: must be BELOW cohort median
CL_NEG_MARKERS = ["CLDN3", "CLDN4", "CLDN7", "CDH1", "ESR1"]

# Positive markers: must be ABOVE cohort median
CL_POS_MARKERS = ["VIM", "CD44", "SNAI1", "ZEB1", "FN1"]

CL_SIGNATURE_GENES = CL_NEG_MARKERS + CL_POS_MARKERS  # 10 total

# ============================================================
# GENE PANELS
# ============================================================

# Luminal identity axis
LUMINAL_AXIS = ["ESR1", "FOXA1", "GATA3", "PGR", "SPDEF", "AR"]

# Claudin / tight junction axis
CLAUDIN_AXIS = ["CLDN3", "CLDN4", "CLDN7", "CLDN8", "CDH1"]

# EMT / partial mesenchymal axis
EMT_AXIS = ["VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1",
            "FN1", "CDH2", "MMP2", "MMP9"]

# Stem programme axis
STEM_AXIS = ["CD44", "ALDH1A1", "KRT14", "KRT5", "TP63",
             "ITGA6", "SOX2", "CD24", "PROM1"]

# Immune programme axis
IMMUNE_AXIS = ["CD274", "CD8A", "FOXP3", "PDCD1", "LAG3",
               "TIGIT", "CD68", "CD163", "IFNG"]

# TP53 target genes (proxy for TP53 pathway activity)
TP53_TARGETS = ["CDKN1A", "MDM2", "BAX", "PUMA", "NOXA",
                "GADD45A", "TP53"]

# Proliferation
PROLIFERATION = ["MKI67", "TOP2A", "PCNA", "CCNB1", "CDK2"]

# HER2 check
HER2_CHECK = ["ERBB2", "ERBB3", "GRB7", "STARD3"]

# Drug target geometry
DRUG_TARGETS = ["TACSTD2", "BRCA1", "PARP1", "PIK3CA", "PTEN",
                "AKT1", "MTOR", "CDK4", "CDK6", "CCND1"]

# Controls (should be flat across all breast subtypes)
CONTROLS = ["CDX2", "SPI1", "MBP", "NKX2-1"]

ALL_GENES = list(dict.fromkeys(
    CL_SIGNATURE_GENES +
    LUMINAL_AXIS + CLAUDIN_AXIS + EMT_AXIS + STEM_AXIS +
    IMMUNE_AXIS + TP53_TARGETS + PROLIFERATION +
    HER2_CHECK + DRUG_TARGETS + CONTROLS
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
# DOWNLOAD UTILITY
# ============================================================

def download_file(url, dest, label="", timeout=180, chunk=1024 * 1024):
    try:
        log(f"  Downloading {label or os.path.basename(url[:70])} ...")
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            total = int(resp.headers.get("Content-Length", 0))
            downloaded = 0
            with open(dest, "wb") as f:
                while True:
                    buf = resp.read(chunk)
                    if not buf:
                        break
                    f.write(buf)
                    downloaded += len(buf)
                    if total > 0:
                        pct = downloaded * 100 // total
                        print(f"    {pct}%  ({downloaded/1e6:.1f} MB)...",
                              end="\r", flush=True)
        print()
        size_mb = os.path.getsize(dest) / 1e6
        log(f"  OK: {os.path.basename(dest)} ({size_mb:.1f} MB)")
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
# STEP 0 — DATA ACQUISITION
# ============================================================

def acquire_data():
    log("=" * 65)
    log("STEP 0: DATA ACQUISITION")
    log(f"  Output: {DATA_DIR}")
    log("=" * 65)

    log("\n-- TCGA-BRCA expression (HiSeqV2) --")
    ok_expr = try_urls(TCGA_EXPR_URLS, TCGA_EXPR_FILE,
                       "TCGA-BRCA expression")
    if not ok_expr:
        log("FATAL: Cannot proceed without expression data.")
        flush_log()
        sys.exit(1)

    log("\n-- TCGA-BRCA clinical matrix --")
    ok_clin = try_urls(TCGA_CLIN_URLS, TCGA_CLIN_FILE,
                       "TCGA-BRCA clinical")
    if not ok_clin:
        log("FATAL: Cannot proceed without clinical data.")
        flush_log()
        sys.exit(1)

# ============================================================
# STEP 1 — LOAD EXPRESSION MATRIX
# Confirmed format from ILC + HER2 + TNBC scripts:
#   Rows = genes (~20530), Columns = sample barcodes (~1218)
#   Values = log2(RSEM+1) — already normalised
# ============================================================

def load_expression():
    log("")
    log("=" * 65)
    log("STEP 1: LOAD EXPRESSION MATRIX")
    log("=" * 65)

    log("  Reading TCGA-BRCA HiSeqV2...")
    try:
        expr = pd.read_csv(TCGA_EXPR_FILE, sep="\t",
                           index_col=0, compression="gzip")
    except Exception:
        try:
            expr = pd.read_csv(TCGA_EXPR_FILE, sep="\t", index_col=0)
        except Exception as e:
            log(f"  FATAL: Cannot load expression: {e}")
            flush_log()
            sys.exit(1)

    log(f"  Raw shape: {expr.shape}")

    # Ensure orientation: rows=genes, cols=samples
    if expr.shape[0] < expr.shape[1]:
        log("  Transposing (samples appear to be rows)...")
        expr = expr.T

    log(f"  Shape (genes x samples): {expr.shape}")
    log(f"  Sample IDs (first 3): {list(expr.columns[:3])}")
    log(f"  Gene IDs (first 3): {list(expr.index[:3])}")

    # Check required genes
    found   = [g for g in ALL_GENES if g in expr.index]
    missing = [g for g in ALL_GENES if g not in expr.index]
    log(f"  Target genes found: {len(found)} / {len(ALL_GENES)}")
    if missing:
        log(f"  Missing genes: {missing}")

    # Check signature genes specifically
    sig_found   = [g for g in CL_SIGNATURE_GENES if g in expr.index]
    sig_missing = [g for g in CL_SIGNATURE_GENES if g not in expr.index]
    log(f"  Signature genes found: {len(sig_found)} / {len(CL_SIGNATURE_GENES)}")
    if sig_missing:
        log(f"  Missing signature genes: {sig_missing}")
        if len(sig_missing) >= 4:
            log("  WARNING: >3 signature genes missing. Classification will be degraded.")

    return expr

# ============================================================
# STEP 2 — LOAD CLINICAL MATRIX
# Confirmed columns from ILC s2 log:
#   PAM50Call_RNAseq — PAM50 subtype calls
#   histological_type — histology
#   ER_Status_nature2012, HER2_Final_Status_nature2012
# ============================================================

def load_clinical():
    log("")
    log("=" * 65)
    log("STEP 2: LOAD CLINICAL MATRIX")
    log("=" * 65)

    try:
        # Detect separator
        with open(TCGA_CLIN_FILE, "rb") as f:
            first_line = f.readline().decode("utf-8", errors="replace")
        sep = "\t" if "\t" in first_line else ","
        log(f"  Separator detected: {'TAB' if sep == chr(9) else 'COMMA'}")

        clin = pd.read_csv(TCGA_CLIN_FILE, sep=sep,
                           index_col=0, low_memory=False)
    except Exception as e:
        log(f"  FATAL: Cannot load clinical matrix: {e}")
        flush_log()
        sys.exit(1)

    # Normalise index: TCGA-XX-XXXX format
    clin.index = [str(s).replace(".", "-") for s in clin.index]
    log(f"  Clinical shape: {clin.shape}")
    log(f"  Clinical index sample: {list(clin.index[:3])}")

    # Report key columns
    key_cols = ["PAM50Call_RNAseq", "PAM50_mRNA_nature2012",
                "histological_type", "ER_Status_nature2012",
                "HER2_Final_Status_nature2012",
                "Integrated_Clusters_with_PAM50__nature2012"]
    for c in key_cols:
        if c in clin.columns:
            vals = clin[c].value_counts().head(8)
            log(f"  {c}:")
            for v, n in vals.items():
                log(f"    {str(v):<45} n={n}")

    return clin

# ============================================================
# STEP 3 — BUILD SAMPLE POPULATIONS
#
# APPROACH:
#   1. Assign tumour / normal from TCGA barcode position 13-14
#      (confirmed method from ILC Script 1 v3)
#   2. Among tumour samples, assign subtype from PAM50Call_RNAseq
#   3. Identify claudin-low from 10-gene geometry signature
#      (does not use PAM50 label — geometry leads)
#
# BARCODE ALIGNMENT:
#   Expression columns: TCGA-XX-XXXX-01A-... (sample barcode, 28 chars)
#   Clinical index:     TCGA-XX-XXXX (patient barcode, 12 chars)
#   Use patient prefix (first 12 chars) to join.
# ============================================================

def build_populations(expr, clin):
    log("")
    log("=" * 65)
    log("STEP 3: BUILD SAMPLE POPULATIONS")
    log("=" * 65)

    samples = list(expr.columns)
    log(f"  Total expression samples: {len(samples)}")
    log(f"  Sample barcode example: {samples[0]}")

    # ── Classify tumour vs normal from barcode ───────────────
    tumour_samples = []
    normal_samples = []
    for s in samples:
        parts = s.split("-")
        # TCGA barcodes: TCGA-XX-XXXX-{sample_type}...
        # sample_type position: field index 3 (0-based)
        if len(parts) >= 4:
            stype = parts[3][:2]  # "01", "11", "06" etc.
        else:
            stype = "??"
        if stype == "01":
            tumour_samples.append(s)
        elif stype == "11":
            normal_samples.append(s)
        # "06" = metastatic: excluded

    log(f"  Tumour samples (01):  {len(tumour_samples)}")
    log(f"  Normal samples (11):  {len(normal_samples)}")

    if len(normal_samples) < MIN_NORMAL:
        log(f"  WARNING: Only {len(normal_samples)} normal samples found.")
        log("  Attempting fallback: clinical sample_type column...")
        # Fallback: build patient prefix map
        patient_prefix_map = {}
        for s in samples:
            pfx = s[:PATIENT_PREFIX_LEN]
            patient_prefix_map.setdefault(pfx, []).append(s)

        # Look for sample_type in clinical columns
        type_col = None
        for cand in ["sample_type", "tissue_status", "_sample_type"]:
            if cand in clin.columns:
                type_col = cand
                break
        if type_col:
            for pat_id in clin.index:
                pat_pfx = str(pat_id)[:PATIENT_PREFIX_LEN]
                stype_val = str(clin.loc[pat_id, type_col]).lower()
                if "normal" in stype_val or "adjacent" in stype_val:
                    for s in patient_prefix_map.get(pat_pfx, []):
                        if s not in normal_samples:
                            normal_samples.append(s)
        log(f"  Normal samples after fallback: {len(normal_samples)}")

    if len(normal_samples) < MIN_NORMAL:
        log(f"  FATAL: Insufficient normal samples ({len(normal_samples)} < {MIN_NORMAL})")
        flush_log()
        sys.exit(1)

    # ── Build patient-prefix → PAM50 map ────────────────────
    # Uses confirmed PAM50Call_RNAseq column from ILC s2 log
    pam50_col = None
    for cand in ["PAM50Call_RNAseq", "PAM50_mRNA_nature2012",
                 "Integrated_Clusters_with_PAM50__nature2012"]:
        if cand in clin.columns:
            pam50_col = cand
            break

    if pam50_col is None:
        log("  WARNING: No PAM50 column found. Cross-check will be skipped.")
        pam50_map = {}
    else:
        log(f"  PAM50 column: '{pam50_col}'")
        pam50_map = {}
        for pat_id, row in clin.iterrows():
            pfx = str(pat_id)[:PATIENT_PREFIX_LEN]
            pam50_map[pfx] = str(row[pam50_col])

    def get_pam50(sample_id):
        return pam50_map.get(sample_id[:PATIENT_PREFIX_LEN], "Unknown")

    # ── Assign conventional subtypes from PAM50 ─────────────
    subtype_buckets = {
        "LumA":   [],
        "LumB":   [],
        "Her2":   [],
        "Basal":  [],
        "Normal": [],
    }
    for s in tumour_samples:
        p = get_pam50(s)
        if p in subtype_buckets:
            subtype_buckets[p].append(s)

    log("")
    log("  PAM50 subtype distribution (tumour samples):")
    for k, v in subtype_buckets.items():
        log(f"    {k:<12} n={len(v)}")
    unassigned = len(tumour_samples) - sum(len(v) for v in subtype_buckets.values())
    log(f"    Unassigned  n={unassigned}")

    # ── CLAUDIN-LOW IDENTIFICATION FROM GEOMETRY ─────────────
    log("")
    log("=" * 65)
    log("STEP 3b: CLAUDIN-LOW IDENTIFICATION FROM 10-GENE SIGNATURE")
    log("  Geometry-first protocol. No PAM50 label used.")
    log("=" * 65)

    # Work on tumour samples only
    tumour_expr = expr[tumour_samples]

    # Compute cohort medians across all tumour samples
    log("  Computing cohort medians across tumour samples...")
    sig_genes_present = [g for g in CL_SIGNATURE_GENES if g in expr.index]
    log(f"  Signature genes present: {sig_genes_present}")

    cohort_medians = {}
    for g in sig_genes_present:
        cohort_medians[g] = float(np.median(tumour_expr.loc[g].values))
        log(f"    {g:<10} cohort median = {cohort_medians[g]:.4f}")

    # Score each tumour sample
    log("")
    log("  Scoring tumour samples against 10-gene signature...")
    scores = {}
    for s in tumour_samples:
        score = 0
        for g in CL_POS_MARKERS:
            if g in sig_genes_present:
                if float(expr.loc[g, s]) > cohort_medians.get(g, 0):
                    score += 1
        for g in CL_NEG_MARKERS:
            if g in sig_genes_present:
                if float(expr.loc[g, s]) < cohort_medians.get(g, 0):
                    score += 1
        scores[s] = score

    score_arr = np.array(list(scores.values()))
    log(f"  Score distribution:")
    for thresh in range(0, 11):
        n = int((score_arr >= thresh).sum())
        log(f"    score >= {thresh}: n={n}")

    # Test thresholds and report
    log("")
    log("  Claudin-low calls at each threshold:")
    cl_samples_by_threshold = {}
    for thresh in CL_SCORE_THRESHOLDS:
        cl_cands = [s for s, sc in scores.items() if sc >= thresh]
        cl_samples_by_threshold[thresh] = cl_cands
        log(f"    Threshold {thresh}: n={len(cl_cands)}")

    # Apply primary threshold
    cl_samples_raw = cl_samples_by_threshold[CL_PRIMARY_THRESHOLD]
    log(f"\n  PRIMARY threshold ({CL_PRIMARY_THRESHOLD}): n={len(cl_samples_raw)}")

    # Exclude ERBB2-high samples (HER2 amplicon overlap)
    if "ERBB2" in expr.index:
        erbb2_vals = tumour_expr.loc["ERBB2", tumour_samples].values
        erbb2_cutoff = np.percentile(erbb2_vals, ERBB2_EXCLUSION_PERCENTILE)
        log(f"  ERBB2 exclusion cutoff (p{ERBB2_EXCLUSION_PERCENTILE}): {erbb2_cutoff:.4f}")
        cl_samples = [s for s in cl_samples_raw
                      if float(expr.loc["ERBB2", s]) <= erbb2_cutoff]
        excluded = len(cl_samples_raw) - len(cl_samples)
        log(f"  ERBB2-high excluded: {excluded}")
        log(f"  Claudin-low final n = {len(cl_samples)}")
    else:
        log("  WARNING: ERBB2 not found. Skipping ERBB2 exclusion.")
        cl_samples = cl_samples_raw

    if len(cl_samples) < MIN_CL_SAMPLES:
        log(f"  WARNING: Only {len(cl_samples)} claudin-low samples identified.")
        log("  Lowering threshold to 6 to recover samples...")
        cl_samples = cl_samples_by_threshold.get(6, cl_samples_raw)
        if "ERBB2" in expr.index:
            erbb2_cutoff = np.percentile(erbb2_vals, ERBB2_EXCLUSION_PERCENTILE)
            cl_samples = [s for s in cl_samples
                          if float(expr.loc["ERBB2", s]) <= erbb2_cutoff]
        log(f"  Claudin-low at threshold 6: n={len(cl_samples)}")

    log(f"\n  FINAL claudin-low n = {len(cl_samples)}")

    # ── PAM50 cross-check (informational — does not filter) ──
    if pam50_col:
        log("")
        log("  PAM50 cross-check for claudin-low samples (informational):")
        pam50_counts = {}
        for s in cl_samples:
            p = get_pam50(s)
            pam50_counts[p] = pam50_counts.get(p, 0) + 1
        for p, n in sorted(pam50_counts.items(), key=lambda x: -x[1]):
            log(f"    {p:<20} n={n}")
        log("  (Geometry-first: PAM50 label is a cross-check, not the classifier)")

    populations = {
        "claudin_low": cl_samples,
        "normal":      normal_samples,
        "luma":        subtype_buckets["LumA"],
        "lumb":        subtype_buckets["LumB"],
        "her2":        subtype_buckets["Her2"],
        "basal":       subtype_buckets["Basal"],  # TNBC
        "ilc":         [],  # ILC not classified by PAM50; set empty
    }

    log("")
    log("  Final population summary:")
    for pop, samps in populations.items():
        log(f"    {pop:<15} n={len(samps)}")

    return populations, scores

# ============================================================
# STEP 4 — UNFILTERED DISCOVERY SCAN
# Claudin-low vs. Normal — all expressed genes
# Geometry is read BEFORE prediction panel is applied
# ============================================================

def unfiltered_discovery(expr, populations):
    log("")
    log("=" * 65)
    log("STEP 4: UNFILTERED DISCOVERY SCAN")
    log("  Geometry first. No panel imposed.")
    log("  Claudin-low vs. Adjacent Normal — all expressed genes.")
    log("=" * 65)

    cl_s  = populations["claudin_low"]
    nor_s = populations["normal"]

    log(f"  Claudin-low samples: {len(cl_s)}")
    log(f"  Normal samples:      {len(nor_s)}")

    if len(cl_s) == 0:
        log("  FATAL: No claudin-low samples for discovery scan.")
        flush_log()
        sys.exit(1)

    cl_expr  = expr[cl_s]
    nor_expr = expr[nor_s]

    results = []
    n_genes  = len(expr.index)
    log(f"  Scanning {n_genes} genes...")

    for i, g in enumerate(expr.index):
        if i % 2000 == 0:
            print(f"    {i}/{n_genes}...", end="\r", flush=True)
        cl_vals  = cl_expr.loc[g].values.astype(float)
        nor_vals = nor_expr.loc[g].values.astype(float)
        mean_cl  = float(np.mean(cl_vals))
        mean_nor = float(np.mean(nor_vals))
        # Skip near-zero in both
        if mean_cl < 0.05 and mean_nor < 0.05:
            continue
        if mean_nor > 0.01:
            pct = (mean_cl - mean_nor) / mean_nor * 100
        elif mean_cl > 0.01:
            pct = float("inf")
        else:
            pct = 0.0
        try:
            _, pval = stats.mannwhitneyu(cl_vals, nor_vals,
                                         alternative="two-sided")
        except Exception:
            pval = 1.0
        results.append({
            "gene": g,
            "cl_mean": mean_cl,
            "nor_mean": mean_nor,
            "pct_change": pct,
            "pval": pval,
        })
    print()

    rdf = pd.DataFrame(results)
    rdf["abs_pct"] = rdf["pct_change"].abs()
    rdf = rdf.sort_values("abs_pct", ascending=False).reset_index(drop=True)
    rdf.to_csv(TOP_MOVERS_FILE, index=False)
    log(f"  Top movers saved: {TOP_MOVERS_FILE}")

    def fmt_p(p):
        if p < 1e-100: return "p<1e-100 ***"
        if p < 1e-10:  return f"p={p:.2e} ***"
        if p < 1e-5:   return f"p={p:.2e} **"
        if p < 0.05:   return f"p={p:.3f} *"
        return         f"p={p:.3f}"

    gained = rdf[rdf["pct_change"] > 0].head(30)
    lost   = rdf[rdf["pct_change"] < 0].sort_values("pct_change").head(30)

    log("")
    log("  TOP 30 GAINED IN CLAUDIN-LOW vs NORMAL:")
    log(f"  {'Gene':<12} {'CL Mean':>9} {'Nor Mean':>9} {'Change':>10}  p-value")
    log("  " + "-" * 62)
    for _, row in gained.iterrows():
        log(f"  {row['gene']:<12} {row['cl_mean']:>9.4f} "
            f"{row['nor_mean']:>9.4f} {row['pct_change']:>+9.1f}%  {fmt_p(row['pval'])}")

    log("")
    log("  TOP 30 LOST IN CLAUDIN-LOW vs NORMAL:")
    log(f"  {'Gene':<12} {'CL Mean':>9} {'Nor Mean':>9} {'Change':>10}  p-value")
    log("  " + "-" * 62)
    for _, row in lost.iterrows():
        log(f"  {row['gene']:<12} {row['cl_mean']:>9.4f} "
            f"{row['nor_mean']:>9.4f} {row['pct_change']:>+9.1f}%  {fmt_p(row['pval'])}")

    # ── Geometry read from unfiltered scan ───────────────────
    log("")
    log("  GEOMETRY READ (unfiltered scan — read before prediction check):")
    log("  ─────────────────────────────────────────────────────────────")

    key_reads = {
        "ESR1":    "Luminal identity TF — expect near-zero",
        "VIM":     "Mesenchymal/stem — expect strongly elevated",
        "CLDN3":   "Claudin-3 — expect strongly lost",
        "CLDN4":   "Claudin-4 — expect strongly lost",
        "CLDN7":   "Claudin-7 — expect strongly lost",
        "CD44":    "Breast stem cell marker — expect elevated",
        "CD24":    "Stem exclusion marker — expect low",
        "ZEB1":    "EMT TF — expect elevated",
        "SNAI1":   "EMT/stem TF — expect elevated",
        "ERBB2":   "HER2 amplicon — expect NOT elevated",
        "CD274":   "PD-L1 — immune programme",
        "FOXA1":   "Luminal pioneer TF — expect lost",
        "GATA3":   "Luminal differentiation TF — expect lost",
        "MKI67":   "Proliferation — relative to TNBC",
    }

    for gene, note in key_reads.items():
        row = rdf[rdf["gene"] == gene]
        if row.empty:
            log(f"  {gene:<10} NOT IN MATRIX — {note}")
            continue
        pct = float(row["pct_change"].values[0])
        pv  = float(row["pval"].values[0])
        rank = int(row.index[0]) + 1
        log(f"  {gene:<10} rank #{rank:<5} {pct:>+8.1f}%  {fmt_p(pv)}")
        log(f"             {note}")

    return rdf

# ============================================================
# STEP 5 — NAMED PANEL TESTS (PREDICTION CHECK)
# ============================================================

def panel_tests(expr, populations, rdf):
    log("")
    log("=" * 65)
    log("STEP 5: NAMED PANEL TESTS")
    log("  Prediction tests run AFTER geometry read.")
    log("=" * 65)

    cl_s  = populations["claudin_low"]
    nor_s = populations["normal"]
    cl_expr  = expr[cl_s]
    nor_expr = expr[nor_s]

    def fmt_p(p):
        if p < 1e-100: return "p<1e-100 ***"
        if p < 1e-10:  return f"p={p:.2e} ***"
        if p < 1e-5:   return f"p={p:.2e} **"
        if p < 0.05:   return f"p={p:.3f} *"
        return         f"p={p:.3f}"

    def test_gene(g):
        """Returns (cl_mean, nor_mean, pct_change, pval) for gene g."""
        if g not in expr.index:
            return None, None, None, None
        cl_v  = cl_expr.loc[g].values.astype(float)
        nor_v = nor_expr.loc[g].values.astype(float)
        mc = float(np.mean(cl_v))
        mn = float(np.mean(nor_v))
        pct = (mc - mn) / mn * 100 if mn > 0.01 else (
            float("inf") if mc > 0.01 else 0.0)
        try:
            _, pval = stats.mannwhitneyu(cl_v, nor_v, alternative="two-sided")
        except Exception:
            pval = 1.0
        return mc, mn, pct, pval

    panel_rows = []

    def section(title):
        log("")
        log(f"  {'─' * 60}")
        log(f"  {title}")
        log(f"  {'─' * 60}")

    def report(pred_id, gene, mc, mn, pct, pval, note=""):
        status = "?"
        if mc is None:
            status = "GENE MISSING"
            log(f"  [{pred_id}] {gene:<10}  MISSING FROM MATRIX  {note}")
        else:
            log(f"  [{pred_id}] {gene:<10}  CL={mc:.4f}  Nor={mn:.4f}  "
                f"{pct:>+8.1f}%  {fmt_p(pval)}  {note}")
        panel_rows.append({
            "prediction": pred_id,
            "gene": gene,
            "cl_mean": mc,
            "nor_mean": mn,
            "pct_change": pct,
            "pval": pval,
            "note": note,
        })

    # ── P1: ESR1/FOXA1/GATA3 LOWEST ─────────────────────────
    section("P1 — ESR1 / FOXA1 / GATA3 / PGR: lowest of all subtypes")
    for g in ["ESR1", "FOXA1", "GATA3", "PGR", "SPDEF", "AR"]:
        mc, mn, pct, pval = test_gene(g)
        report("P1", g, mc, mn, pct, pval)

    # ── P2: VIM / ZEB1 / SNAI1 ELEVATED ─────────────────────
    section("P2 — VIM / ZEB1 / SNAI1 / FN1: partial EMT programme")
    for g in ["VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "FN1", "CDH2"]:
        mc, mn, pct, pval = test_gene(g)
        report("P2", g, mc, mn, pct, pval)

    # ── P3: CLDN3 / CLDN4 / CLDN7 LOWEST ────────────────────
    section("P3 — CLDN3 / CLDN4 / CLDN7: lowest of all subtypes")
    for g in ["CLDN3", "CLDN4", "CLDN7", "CLDN8", "CDH1"]:
        mc, mn, pct, pval = test_gene(g)
        report("P3", g, mc, mn, pct, pval)

    # ── P4: CD44 HIGH / CD24 LOW ─────────────────────────────
    section("P4 — CD44 high / CD24 low: canonical breast stem phenotype")
    for g in ["CD44", "CD24", "ALDH1A1", "KRT14", "KRT5", "TP63", "ITGA6", "PROM1"]:
        mc, mn, pct, pval = test_gene(g)
        report("P4", g, mc, mn, pct, pval)

    # ── P5: TP53 TARGET SIGNATURE ────────────────────────────
    section("P5 — TP53 target signature: less disrupted than TNBC")
    for g in ["CDKN1A", "MDM2", "BAX", "GADD45A", "TP53", "PUMA"]:
        mc, mn, pct, pval = test_gene(g)
        report("P5", g, mc, mn, pct, pval)

    # ── P6: IMMUNE PROGRAMME ─────────────────────────────────
    section("P6 — Immune programme: PD-L1 / CD8A elevated")
    for g in ["CD274", "CD8A", "FOXP3", "PDCD1", "LAG3", "TIGIT",
              "CD68", "CD163", "IFNG"]:
        mc, mn, pct, pval = test_gene(g)
        report("P6", g, mc, mn, pct, pval)

    # ── P7: DEPTH AXIS STRUCTURE ─────────────────────────────
    section("P7 — Depth axis: CLDN3 orthogonal to VIM within claudin-low")
    if "CLDN3" in expr.index and "VIM" in expr.index:
        cldn3_vals = cl_expr.loc["CLDN3"].values.astype(float)
        vim_vals   = cl_expr.loc["VIM"].values.astype(float)
        r_cldn3_vim, p_cldn3_vim = stats.pearsonr(cldn3_vals, vim_vals)
        log(f"  Within claudin-low: r(CLDN3, VIM) = {r_cldn3_vim:.4f}  "
            f"{fmt_p(p_cldn3_vim)}")
        log(f"  Prediction: CLDN3 and VIM should be orthogonal (r near 0)")
        log(f"  Prediction: Negative r expected (as VIM rises, CLDN3 falls)")

        # CLDN3 vs luminal markers within claudin-low
        log("")
        log("  CLDN3 vs luminal identity markers (within claudin-low):")
        for g in ["ESR1", "FOXA1", "GATA3"]:
            if g in expr.index:
                g_vals = cl_expr.loc[g].values.astype(float)
                r, p = stats.pearsonr(cldn3_vals, g_vals)
                log(f"    r(CLDN3, {g}) = {r:.4f}  {fmt_p(p)}")
                log(f"    Prediction: positive r (residual claudin retention "
                    f"correlates with partial luminal identity)")

        # VIM vs stem markers within claudin-low
        log("")
        log("  VIM vs stem markers (within claudin-low):")
        for g in ["CD44", "ALDH1A1", "KRT14", "SNAI1", "ZEB1"]:
            if g in expr.index:
                g_vals = cl_expr.loc[g].values.astype(float)
                r, p = stats.pearsonr(vim_vals, g_vals)
                log(f"    r(VIM, {g}) = {r:.4f}  {fmt_p(p)}")
    else:
        log("  Cannot test P7: CLDN3 or VIM missing from expression matrix.")

    # ── P8: CROSS-SUBTYPE PCA POSITION ───────────────────────
    section("P8 — Cross-subtype PCA: claudin-low most distant from LumA")
    log("  (Computed in STEP 6 — cross-subtype comparison)")

    # ── P9: DRUG TARGET GEOMETRY ─────────────────────────────
    section("P9 — Drug target geometry: immune + stem, NOT ER or HER2")
    for g in ["TACSTD2", "BRCA1", "PARP1", "ERBB2", "ESR1",
              "CD274", "ALDH1A1", "PIK3CA", "PTEN"]:
        mc, mn, pct, pval = test_gene(g)
        report("P9", g, mc, mn, pct, pval)

    # ── Proliferation panel ───────────────────────────────────
    section("PROLIFERATION — context for cross-subtype comparison")
    for g in ["MKI67", "TOP2A", "PCNA", "CCNB1"]:
        mc, mn, pct, pval = test_gene(g)
        report("PROLIF", g, mc, mn, pct, pval)

    # ── Save panel results ────────────────────────────────────
    pdf = pd.DataFrame(panel_rows)
    pdf.to_csv(CSV_FILE, index=False)
    log("")
    log(f"  Panel results saved: {CSV_FILE}")

    return pdf

# ============================================================
# STEP 6 — CROSS-SUBTYPE COMPARISON
# Claudin-low vs. Normal vs. all five completed subtypes
# on the same gene panel
# ============================================================

def cross_subtype_comparison(expr, populations):
    log("")
    log("=" * 65)
    log("STEP 6: CROSS-SUBTYPE COMPARISON")
    log("  Claudin-low positioned relative to all BRCA subtypes.")
    log("  Same dataset, same reference, same gene panel.")
    log("=" * 65)

    pop_order = ["normal", "claudin_low", "luma", "lumb", "her2", "basal"]
    pop_labels = {
        "normal":      "Normal",
        "claudin_low": "Claudin-low",
        "luma":        "LumA",
        "lumb":        "LumB",
        "her2":        "HER2-enriched",
        "basal":       "TNBC/Basal",
    }

    # Key genes for cross-subtype table
    cross_genes = (
        ["ESR1", "FOXA1", "GATA3", "PGR"] +
        ["CLDN3", "CLDN4", "CLDN7"] +
        ["VIM", "ZEB1", "SNAI1", "CDH1"] +
        ["CD44", "CD24", "ALDH1A1", "KRT14", "KRT5"] +
        ["CD274", "CD8A"] +
        ["MKI67", "ERBB2"] +
        ["CDKN1A", "MDM2"]
    )
    cross_genes = [g for g in cross_genes if g in expr.index]

    rows = []
    for g in cross_genes:
        row = {"gene": g}
        for pop in pop_order:
            samps = populations.get(pop, [])
            if len(samps) == 0:
                row[pop_labels[pop]] = None
                continue
            vals = expr.loc[g, samps].values.astype(float)
            row[pop_labels[pop]] = float(np.mean(vals))
        rows.append(row)

    cdf = pd.DataFrame(rows)
    cdf.to_csv(CROSS_FILE, index=False)

    # Print table
    header = f"  {'Gene':<12}" + "".join(
        f"{pop_labels[p]:>14}" for p in pop_order
    )
    log(header)
    log("  " + "─" * (12 + 14 * len(pop_order)))

    for _, row in cdf.iterrows():
        line = f"  {row['gene']:<12}"
        for p in pop_order:
            lbl = pop_labels[p]
            val = row[lbl]
            line += f"{val:>14.4f}" if val is not None else f"{'N/A':>14}"
        log(line)

    log(f"\n  Cross-subtype table saved: {CROSS_FILE}")

    # ── PCA cross-subtype (P8 test) ──────────────────────────
    log("")
    log("  PCA — cross-subtype position (P8 test):")

    pca_pops = {k: v for k, v in populations.items()
                if len(v) > 0 and k != "ilc"}
    pop_centroids = {}
    for pop, samps in pca_pops.items():
        if len(samps) == 0:
            continue
        vals = expr.loc[cross_genes, samps].values.astype(float)
        pop_centroids[pop] = np.mean(vals, axis=1)

    if len(pop_centroids) < 2:
        log("  Insufficient populations for PCA. Skipping.")
        return cdf

    centroid_df = pd.DataFrame(pop_centroids, index=cross_genes).T

    try:
        scaler = StandardScaler()
        scaled = scaler.fit_transform(centroid_df.values)
        pca = PCA(n_components=min(3, scaled.shape[0] - 1))
        pcs = pca.fit_transform(scaled)
        log(f"  PCA variance explained: PC1={pca.explained_variance_ratio_[0]*100:.1f}%  "
            f"PC2={pca.explained_variance_ratio_[1]*100:.1f}%")

        log("")
        log(f"  {'Population':<20} {'PC1':>8} {'PC2':>8}")
        log("  " + "─" * 38)
        for i, pop in enumerate(pop_centroids.keys()):
            log(f"  {pop_labels.get(pop, pop):<20} {pcs[i, 0]:>8.3f} {pcs[i, 1]:>8.3f}")

        # Pairwise distances
        log("")
        log("  Pairwise centroid distances (Euclidean on PCA space):")
        pop_list = list(pop_centroids.keys())
        cl_idx = pop_list.index("claudin_low") if "claudin_low" in pop_list else None

        dist_rows = []
        for i, p1 in enumerate(pop_list):
            for j, p2 in enumerate(pop_list):
                if j <= i:
                    continue
                d = float(np.linalg.norm(pcs[i] - pcs[j]))
                dist_rows.append((pop_labels.get(p1, p1),
                                  pop_labels.get(p2, p2), d))
                log(f"    {pop_labels.get(p1, p1):<20} ↔ "
                    f"{pop_labels.get(p2, p2):<20} d={d:.4f}")

        if cl_idx is not None:
            log("")
            log("  Distances FROM claudin-low to each subtype (P8 test):")
            log("  Prediction: max distance to LumA")
            for j, p2 in enumerate(pop_list):
                if j == cl_idx:
                    continue
                d = float(np.linalg.norm(pcs[cl_idx] - pcs[j]))
                log(f"    Claudin-low ↔ {pop_labels.get(p2, p2):<20} d={d:.4f}")

    except Exception as e:
        log(f"  PCA failed: {e}")

    return cdf

# ============================================================
# STEP 7 — STEM DEPTH AXIS (P7 — within claudin-low)
# ============================================================

def stem_depth_axis(expr, populations):
    log("")
    log("=" * 65)
    log("STEP 7: STEM DEPTH AXIS (P7 — within claudin-low)")
    log("  How deeply locked in the false attractor is each sample?")
    log("=" * 65)

    cl_s = populations["claudin_low"]
    if len(cl_s) < 5:
        log("  Insufficient claudin-low samples for depth axis.")
        return

    cl_expr = expr[cl_s]

    # Depth score = mean of positive markers (normalised) minus
    # mean of negative markers (normalised) — all within CL
    # Normalise each gene to its range across CL samples
    pos_genes = [g for g in CL_POS_MARKERS if g in expr.index]
    neg_genes = [g for g in CL_NEG_MARKERS if g in expr.index]

    log(f"  Positive markers used: {pos_genes}")
    log(f"  Negative markers used: {neg_genes}")

    def z_score(arr):
        s = arr.std()
        if s < 1e-9:
            return np.zeros_like(arr)
        return (arr - arr.mean()) / s

    pos_matrix = np.column_stack([
        z_score(cl_expr.loc[g].values.astype(float))
        for g in pos_genes
    ]) if pos_genes else np.zeros((len(cl_s), 1))

    neg_matrix = np.column_stack([
        z_score(cl_expr.loc[g].values.astype(float))
        for g in neg_genes
    ]) if neg_genes else np.zeros((len(cl_s), 1))

    # Depth = high positive + low negative (invert neg markers)
    depth_scores = (np.mean(pos_matrix, axis=1) -
                    np.mean(neg_matrix, axis=1))

    depth_df = pd.DataFrame({
        "sample":      cl_s,
        "depth_score": depth_scores,
    }).sort_values("depth_score", ascending=False)

    depth_df.to_csv(DEPTH_FILE, index=False)
    log(f"  Depth scores computed for {len(cl_s)} claudin-low samples.")
    log(f"  Score range: {depth_scores.min():.4f} to {depth_scores.max():.4f}")
    log(f"  Score mean:  {depth_scores.mean():.4f}")
    log(f"  Score std:   {depth_scores.std():.4f}")
    log(f"  Depth scores saved: {DEPTH_FILE}")

    # Correlate depth with individual genes to validate
    log("")
    log("  Gene correlations with depth score (within claudin-low):")
    log("  (High positive r = gene marks deep attractor state)")
    log(f"  {'Gene':<12} {'r':>8} {'p':>14}")
    log("  " + "─" * 38)

    test_genes = pos_genes + neg_genes + [
        "FOXA1", "ESR1", "GATA3", "MKI67", "CD274", "ALDH1A1", "PROM1"
    ]
    corr_rows = []
    for g in test_genes:
        if g not in expr.index:
            continue
        g_vals = cl_expr.loc[g].values.astype(float)
        if g_vals.std() < 1e-9:
            continue
        r, p = stats.pearsonr(depth_scores, g_vals)
        log(f"  {g:<12} {r:>8.4f} {'p=' + f'{p:.2e}':>14}")
        corr_rows.append({"gene": g, "r_with_depth": r, "pval": p})

    corr_df = pd.DataFrame(corr_rows).sort_values("r_with_depth", ascending=False)
    corr_path = os.path.join(RESULTS_DIR, "cl_s1_depth_correlations.csv")
    corr_df.to_csv(corr_path, index=False)
    log(f"\n  Depth correlations saved: {corr_path}")

    return depth_df

# ============================================================
# STEP 8 — FIGURE
# ============================================================

def make_figure(expr, populations, rdf, cdf, depth_df):
    log("")
    log("=" * 65)
    log("STEP 8: GENERATE FIGURE")
    log("=" * 65)

    cl_s  = populations["claudin_low"]
    nor_s = populations["normal"]

    if len(cl_s) == 0:
        log("  No claudin-low samples — skipping figure.")
        return

    fig = plt.figure(figsize=(22, 20))
    fig.patch.set_facecolor("#0d1117")
    gs = gridspec.GridSpec(3, 3, figure=fig,
                           hspace=0.45, wspace=0.40)

    panel_style = dict(facecolor="#161b22", edgecolor="#30363d")
    title_kw    = dict(color="#e6edf3", fontsize=11, fontweight="bold", pad=8)
    label_kw    = dict(color="#8b949e", fontsize=8)
    tick_kw     = dict(colors="#8b949e", labelsize=7)

    BLUE   = "#58a6ff"
    GREEN  = "#3fb950"
    RED    = "#f85149"
    ORANGE = "#d29922"
    PURPLE = "#bc8cff"
    GREY   = "#6e7681"

    subtype_colors = {
        "normal":      GREY,
        "claudin_low": BLUE,
        "luma":        GREEN,
        "lumb":        ORANGE,
        "her2":        PURPLE,
        "basal":       RED,
    }
    subtype_labels = {
        "normal":      "Normal",
        "claudin_low": "Claudin-low",
        "luma":        "LumA",
        "lumb":        "LumB",
        "her2":        "HER2",
        "basal":       "TNBC",
    }

    def ax_style(ax):
        ax.set_facecolor("#161b22")
        for spine in ax.spines.values():
            spine.set_edgecolor("#30363d")
        ax.tick_params(colors="#8b949e", labelsize=7)
        ax.xaxis.label.set_color("#8b949e")
        ax.yaxis.label.set_color("#8b949e")

    # ── Panel 1: 10-gene signature expression — CL vs Normal ─
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title("10-Gene Claudin-low Signature", **title_kw)
    sig_present = [g for g in CL_SIGNATURE_GENES if g in expr.index]
    if sig_present and len(cl_s) > 0:
        cl_means  = [float(np.mean(expr.loc[g, cl_s].values))  for g in sig_present]
        nor_means = [float(np.mean(expr.loc[g, nor_s].values)) for g in sig_present]
        x = np.arange(len(sig_present))
        w = 0.35
        ax1.bar(x - w/2, nor_means, w, color=GREY,  label="Normal",  alpha=0.85)
        ax1.bar(x + w/2, cl_means,  w, color=BLUE,  label="CL",      alpha=0.85)
        ax1.set_xticks(x)
        ax1.set_xticklabels(sig_present, rotation=45, ha="right", fontsize=7)
        ax1.set_ylabel("log2(RSEM+1)", **label_kw)
        ax1.legend(fontsize=7, labelcolor="#8b949e",
                   facecolor="#161b22", edgecolor="#30363d")
        # Mark positive vs negative with colour
        for i, g in enumerate(sig_present):
            colour = GREEN if g in CL_POS_MARKERS else RED
            ax1.axvspan(i - 0.5, i + 0.5, alpha=0.05, color=colour)
    ax_style(ax1)

    # ── Panel 2: Luminal vs EMT axis — cross-subtype ─────────
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_title("Luminal vs EMT/Stem Axis (cross-subtype)", **title_kw)
    axis_genes = ["ESR1", "FOXA1", "GATA3", "VIM", "ZEB1", "SNAI1"]
    axis_genes = [g for g in axis_genes if g in expr.index]
    pop_order  = ["normal", "luma", "her2", "basal", "claudin_low"]
    if axis_genes:
        n_pops  = len(pop_order)
        n_genes = len(axis_genes)
        bar_w   = 0.8 / n_pops
        x = np.arange(n_genes)
        for pi, pop in enumerate(pop_order):
            samps = populations.get(pop, [])
            if not samps:
                continue
            means = [float(np.mean(expr.loc[g, samps].values))
                     if g in expr.index else 0 for g in axis_genes]
            offset = (pi - n_pops / 2 + 0.5) * bar_w
            ax2.bar(x + offset, means, bar_w,
                    label=subtype_labels[pop],
                    color=subtype_colors[pop], alpha=0.85)
        ax2.set_xticks(x)
        ax2.set_xticklabels(axis_genes, rotation=45, ha="right", fontsize=7)
        ax2.set_ylabel("log2(RSEM+1)", **label_kw)
        ax2.legend(fontsize=6, labelcolor="#8b949e",
                   facecolor="#161b22", edgecolor="#30363d")
    ax_style(ax2)

    # ── Panel 3: Claudin junction axis — cross-subtype ───────
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.set_title("Claudin Junction Axis (cross-subtype)", **title_kw)
    cl_genes_plot = ["CLDN3", "CLDN4", "CLDN7", "CDH1"]
    cl_genes_plot = [g for g in cl_genes_plot if g in expr.index]
    if cl_genes_plot:
        n_pops  = len(pop_order)
        n_genes = len(cl_genes_plot)
        bar_w   = 0.8 / n_pops
        x = np.arange(n_genes)
        for pi, pop in enumerate(pop_order):
            samps = populations.get(pop, [])
            if not samps:
                continue
            means = [float(np.mean(expr.loc[g, samps].values))
                     if g in expr.index else 0 for g in cl_genes_plot]
            offset = (pi - n_pops / 2 + 0.5) * bar_w
            ax3.bar(x + offset, means, bar_w,
                    label=subtype_labels[pop],
                    color=subtype_colors[pop], alpha=0.85)
        ax3.set_xticks(x)
        ax3.set_xticklabels(cl_genes_plot, rotation=45, ha="right", fontsize=7)
        ax3.set_ylabel("log2(RSEM+1)", **label_kw)
        ax3.legend(fontsize=6, labelcolor="#8b949e",
                   facecolor="#161b22", edgecolor="#30363d")
    ax_style(ax3)

    # ── Panel 4: Stem marker axis — CL vs TNBC vs Normal ─────
    ax4 = fig.add_subplot(gs[1, 0])
    ax4.set_title("Stem Marker Axis: CL vs TNBC vs Normal", **title_kw)
    stem_genes_plot = ["CD44", "CD24", "ALDH1A1", "KRT14", "KRT5", "ITGA6"]
    stem_genes_plot = [g for g in stem_genes_plot if g in expr.index]
    compare_pops = ["normal", "basal", "claudin_low"]
    if stem_genes_plot:
        n_pops  = len(compare_pops)
        bar_w   = 0.8 / n_pops
        x = np.arange(len(stem_genes_plot))
        for pi, pop in enumerate(compare_pops):
            samps = populations.get(pop, [])
            if not samps:
                continue
            means = [float(np.mean(expr.loc[g, samps].values))
                     if g in expr.index else 0 for g in stem_genes_plot]
            offset = (pi - n_pops / 2 + 0.5) * bar_w
            ax4.bar(x + offset, means, bar_w,
                    label=subtype_labels[pop],
                    color=subtype_colors[pop], alpha=0.85)
        ax4.set_xticks(x)
        ax4.set_xticklabels(stem_genes_plot, rotation=45, ha="right", fontsize=7)
        ax4.set_ylabel("log2(RSEM+1)", **label_kw)
        ax4.legend(fontsize=7, labelcolor="#8b949e",
                   facecolor="#161b22", edgecolor="#30363d")
    ax_style(ax4)

    # ── Panel 5: Immune programme — CL vs TNBC vs LumA ───────
    ax5 = fig.add_subplot(gs[1, 1])
    ax5.set_title("Immune Programme (P6)", **title_kw)
    immune_plot = ["CD274", "CD8A", "FOXP3", "PDCD1", "LAG3", "CD68"]
    immune_plot = [g for g in immune_plot if g in expr.index]
    compare_pops2 = ["normal", "luma", "basal", "claudin_low"]
    if immune_plot:
        n_pops = len(compare_pops2)
        bar_w  = 0.8 / n_pops
        x = np.arange(len(immune_plot))
        for pi, pop in enumerate(compare_pops2):
            samps = populations.get(pop, [])
            if not samps:
                continue
            means = [float(np.mean(expr.loc[g, samps].values))
                     if g in expr.index else 0 for g in immune_plot]
            offset = (pi - n_pops / 2 + 0.5) * bar_w
            ax5.bar(x + offset, means, bar_w,
                    label=subtype_labels[pop],
                    color=subtype_colors[pop], alpha=0.85)
        ax5.set_xticks(x)
        ax5.set_xticklabels(immune_plot, rotation=45, ha="right", fontsize=7)
        ax5.set_ylabel("log2(RSEM+1)", **label_kw)
        ax5.legend(fontsize=7, labelcolor="#8b949e",
                   facecolor="#161b22", edgecolor="#30363d")
    ax_style(ax5)

    # ── Panel 6: Depth score distribution ────────────────────
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.set_title("Stem Depth Score Distribution (within CL)", **title_kw)
    if depth_df is not None and len(depth_df) > 0:
        scores = depth_df["depth_score"].values
        ax6.hist(scores, bins=min(20, len(scores)//2 + 1),
                 color=BLUE, alpha=0.8, edgecolor="#30363d")
        ax6.axvline(np.median(scores), color=ORANGE, lw=1.5,
                    linestyle="--", label=f"Median={np.median(scores):.2f}")
        ax6.set_xlabel("Stem Lock Depth Score", **label_kw)
        ax6.set_ylabel("Count", **label_kw)
        ax6.legend(fontsize=7, labelcolor="#8b949e",
                   facecolor="#161b22", edgecolor="#30363d")
    else:
        ax6.text(0.5, 0.5, "No depth data", transform=ax6.transAxes,
                 ha="center", va="center", color="#8b949e")
    ax_style(ax6)

    # ── Panel 7: Cross-subtype PCA ────────────────────────────
    ax7 = fig.add_subplot(gs[2, 0:2])
    ax7.set_title("Cross-subtype PCA (P8 — claudin-low vs. all subtypes)", **title_kw)

    pca_pops = {k: v for k, v in populations.items() if len(v) > 0}
    pca_genes = [g for g in (LUMINAL_AXIS + CLAUDIN_AXIS + EMT_AXIS[:4] +
                              STEM_AXIS[:4] + IMMUNE_AXIS[:3])
                 if g in expr.index]

    if len(pca_genes) > 3 and len(pca_pops) > 1:
        try:
            all_samps = []
            samp_labels_pca = []
            for pop, samps in pca_pops.items():
                all_samps.extend(samps)
                samp_labels_pca.extend([pop] * len(samps))

            mat = expr.loc[pca_genes, all_samps].values.T.astype(float)
            scaler = StandardScaler()
            mat_sc = scaler.fit_transform(mat)
            pca = PCA(n_components=2)
            pcs  = pca.fit_transform(mat_sc)

            for pop in pca_pops.keys():
                idx = [i for i, l in enumerate(samp_labels_pca) if l == pop]
                if not idx:
                    continue
                ax7.scatter(pcs[idx, 0], pcs[idx, 1],
                            c=subtype_colors.get(pop, GREY),
                            label=subtype_labels.get(pop, pop),
                            alpha=0.5, s=15, edgecolors="none")

            ax7.set_xlabel(
                f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)",
                **label_kw)
            ax7.set_ylabel(
                f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)",
                **label_kw)
            ax7.legend(fontsize=7, labelcolor="#8b949e",
                       facecolor="#161b22", edgecolor="#30363d",
                       markerscale=1.5)
        except Exception as e:
            ax7.text(0.5, 0.5, f"PCA failed: {e}",
                     transform=ax7.transAxes,
                     ha="center", va="center", color="#f85149")
    ax_style(ax7)

    # ── Panel 8: CLDN3 vs VIM scatter within CL (P7 depth) ───
    ax8 = fig.add_subplot(gs[2, 2])
    ax8.set_title("CLDN3 vs VIM within Claudin-low (P7)", **title_kw)
    if ("CLDN3" in expr.index and "VIM" in expr.index and len(cl_s) > 3):
        cldn3_v = expr.loc["CLDN3", cl_s].values.astype(float)
        vim_v   = expr.loc["VIM",   cl_s].values.astype(float)
        r, p    = stats.pearsonr(cldn3_v, vim_v)
        ax8.scatter(cldn3_v, vim_v, c=BLUE, alpha=0.6, s=20, edgecolors="none")
        # Trend line
        if len(cldn3_v) > 2:
            z = np.polyfit(cldn3_v, vim_v, 1)
            xr = np.linspace(cldn3_v.min(), cldn3_v.max(), 50)
            ax8.plot(xr, np.polyval(z, xr), color=ORANGE,
                     lw=1.5, linestyle="--")
        ax8.set_xlabel("CLDN3 expression", **label_kw)
        ax8.set_ylabel("VIM expression",   **label_kw)
        ax8.text(0.05, 0.92, f"r={r:.3f}  p={p:.2e}",
                 transform=ax8.transAxes, color="#e6edf3", fontsize=8)
        ax8.text(0.05, 0.84,
                 "Prediction: near-zero or negative r",
                 transform=ax8.transAxes, color="#8b949e", fontsize=7)
    else:
        ax8.text(0.5, 0.5, "CLDN3 or VIM missing",
                 transform=ax8.transAxes,
                 ha="center", va="center", color="#8b949e")
    ax_style(ax8)

    # ── Title ─────────────────────────────────────────────────
    fig.suptitle(
        "CLAUDIN-LOW — SCRIPT 1 GEOMETRY\n"
        "OrganismCore | BRCA-S7b | 2026-03-05\n"
        "TYPE 4 — STEM LOCK: Differentiation Arrest at Progenitor Root",
        color="#e6edf3", fontsize=13, fontweight="bold", y=0.98
    )

    plt.savefig(FIG_FILE, dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close()
    log(f"  Figure saved: {FIG_FILE}")

# ============================================================
# STEP 9 — PREDICTION SCORECARD
# ============================================================

def prediction_scorecard(expr, populations, rdf, pdf):
    log("")
    log("=" * 65)
    log("STEP 9: PREDICTION SCORECARD")
    log("  Protocol v2.0: geometry was read first (Steps 4–6).")
    log("  Predictions now scored against that geometry.")
    log("=" * 65)

    cl_s    = populations["claudin_low"]
    nor_s   = populations["normal"]
    luma_s  = populations["luma"]
    basal_s = populations["basal"]

    def mean(samps, gene):
        if not samps or gene not in expr.index:
            return None
        return float(np.mean(expr.loc[gene, samps].values))

    log("")
    log("  P1 — ESR1/FOXA1/GATA3: lowest of all BRCA subtypes")
    for g in ["ESR1", "FOXA1", "GATA3"]:
        cl_m    = mean(cl_s, g)
        luma_m  = mean(luma_s, g)
        basal_m = mean(basal_s, g)
        if cl_m is None:
            log(f"    {g}: MISSING")
            continue
        cl_vs_basal = "LOWER ✓" if (basal_m and cl_m < basal_m) else "NOT LOWER ✗"
        cl_vs_luma  = "LOWER ✓" if (luma_m  and cl_m < luma_m)  else "NOT LOWER ✗"
        log(f"    {g}:  CL={cl_m:.4f}  Basal={basal_m:.4f}  LumA={luma_m:.4f}  "
            f"vs-Basal:{cl_vs_basal}  vs-LumA:{cl_vs_luma}")

    log("")
    log("  P2 — VIM/ZEB1/SNAI1: partial EMT elevated")
    for g in ["VIM", "ZEB1", "SNAI1"]:
        cl_m   = mean(cl_s, g)
        nor_m  = mean(nor_s, g)
        luma_m = mean(luma_s, g)
        if cl_m is None:
            log(f"    {g}: MISSING")
            continue
        vs_nor  = "UP ✓"  if (nor_m  and cl_m > nor_m)  else "NOT UP ✗"
        vs_luma = "UP ✓"  if (luma_m and cl_m > luma_m) else "NOT UP ✗"
        log(f"    {g}:  CL={cl_m:.4f}  Normal={nor_m:.4f}  LumA={luma_m:.4f}  "
            f"vs-Normal:{vs_nor}  vs-LumA:{vs_luma}")

    log("")
    log("  P3 — CLDN3/CLDN4/CLDN7: lowest of all subtypes")
    for g in ["CLDN3", "CLDN4", "CLDN7"]:
        cl_m    = mean(cl_s, g)
        basal_m = mean(basal_s, g)
        luma_m  = mean(luma_s, g)
        nor_m   = mean(nor_s, g)
        if cl_m is None:
            log(f"    {g}: MISSING")
            continue
        vs_basal = "LOWER ✓" if (basal_m and cl_m < basal_m) else "NOT LOWER ✗"
        vs_nor   = "LOWER ✓" if (nor_m   and cl_m < nor_m)   else "NOT LOWER ✗"
        log(f"    {g}:  CL={cl_m:.4f}  Basal={basal_m:.4f}  Nor={nor_m:.4f}  "
            f"vs-Basal:{vs_basal}  vs-Normal:{vs_nor}")

    log("")
    log("  P4 — CD44 high / CD24 low")
    cd44_cl  = mean(cl_s, "CD44")
    cd44_nor = mean(nor_s, "CD44")
    cd24_cl  = mean(cl_s, "CD24")
    cd24_nor = mean(nor_s, "CD24")
    if cd44_cl and cd44_nor:
        log(f"    CD44:  CL={cd44_cl:.4f}  Normal={cd44_nor:.4f}  "
            f"{'UP ✓' if cd44_cl > cd44_nor else 'NOT UP ✗'}")
    if cd24_cl and cd24_nor:
        log(f"    CD24:  CL={cd24_cl:.4f}  Normal={cd24_nor:.4f}  "
            f"{'DOWN ✓' if cd24_cl < cd24_nor else 'NOT DOWN ✗'}")

    log("")
    log("  P6 — Immune programme: PD-L1 elevated")
    pdl1_cl   = mean(cl_s, "CD274")
    pdl1_nor  = mean(nor_s, "CD274")
    pdl1_luma = mean(luma_s, "CD274")
    if pdl1_cl and pdl1_nor:
        log(f"    CD274: CL={pdl1_cl:.4f}  Normal={pdl1_nor:.4f}  "
            f"LumA={pdl1_luma:.4f}  "
            f"{'UP ✓' if pdl1_cl > pdl1_nor else 'NOT UP ✗'}")

    log("")
    log("  P9 — HER2 (ERBB2) confirmation: expect NOT elevated")
    erbb2_cl  = mean(cl_s, "ERBB2")
    erbb2_nor = mean(nor_s, "ERBB2")
    if erbb2_cl and erbb2_nor:
        log(f"    ERBB2: CL={erbb2_cl:.4f}  Normal={erbb2_nor:.4f}  "
            f"{'NOT amplified ✓' if erbb2_cl < erbb2_nor * 1.5 else 'ELEVATED — check exclusion ✗'}")

# ============================================================
# MAIN
# ============================================================

def main():
    log("=" * 65)
    log("BRCA CLAUDIN-LOW — SCRIPT 1")
    log("OrganismCore — Document BRCA-S7a/b | 2026-03-05")
    log("")
    log("ATTRACTOR TYPE: TYPE 4 — STEM LOCK")
    log("DIFFERENTIATION ARREST AT PROGENITOR ROOT")
    log("Cell of origin: Mammary stem cell")
    log("Depth axis: Stem programme depth")
    log(f"Output directory: {RESULTS_DIR}")
    log("=" * 65)

    # Step 0: acquire data
    acquire_data()

    # Step 1: load expression
    expr = load_expression()

    # Step 2: load clinical
    clin = load_clinical()

    # Step 3: build populations + identify claudin-low from geometry
    populations, scores = build_populations(expr, clin)

    # Step 4: unfiltered discovery scan (geometry first)
    rdf = unfiltered_discovery(expr, populations)

    # Step 5: named panel tests (prediction tests, run after geometry)
    pdf = panel_tests(expr, populations, rdf)

    # Step 6: cross-subtype comparison
    cdf = cross_subtype_comparison(expr, populations)

    # Step 7: stem depth axis
    depth_df = stem_depth_axis(expr, populations)

    # Step 8: figure
    try:
        make_figure(expr, populations, rdf, cdf, depth_df)
    except Exception as e:
        log(f"  Figure generation failed: {e}")

    # Step 9: prediction scorecard
    prediction_scorecard(expr, populations, rdf, pdf)

    # Final summary
    log("")
    log("=" * 65)
    log("SCRIPT 1 COMPLETE")
    log("=" * 65)
    log(f"  Claudin-low samples identified: {len(populations['claudin_low'])}")
    log(f"  Normal reference samples:       {len(populations['normal'])}")
    log(f"  Output directory: {RESULTS_DIR}")
    log(f"  Log:       {LOG_FILE}")
    log(f"  Figure:    {FIG_FILE}")
    log(f"  Panel CSV: {CSV_FILE}")
    log(f"  Top movers CSV: {TOP_MOVERS_FILE}")
    log(f"  Cross-subtype CSV: {CROSS_FILE}")
    log(f"  Depth axis CSV: {DEPTH_FILE}")
    log("")
    log("  NEXT STEP:")
    log("    Write script1_results_and_reasoning.md (BRCA-S7b)")
    log("    Read geometry on its own terms BEFORE scoring predictions.")
    log("    Then write before_script2.md (BRCA-S7c) with locked predictions.")

    flush_log()


if __name__ == "__main__":
    main()
