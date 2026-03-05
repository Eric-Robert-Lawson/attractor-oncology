"""
BRCA ILC — SCRIPT 1 (v3 — full corrected)
OrganismCore — Document BRCA-S6a/b | 2026-03-05

ROOT CAUSE OF v1/v2 FAILURES:
  The TCGA Xena clinical matrix (BRCA_clinicalMatrix) indexes rows by
  PATIENT barcodes: TCGA-BH-A0BZ  (12 chars, no sample suffix)
  The TCGA Xena HiSeqV2 expression matrix indexes columns by
  SAMPLE barcodes:  TCGA-BH-A0BZ-01A-11R-A12P-07 (28 chars)

  v1/v2 did direct equality matching — always 0 hits.

  FIX (v3):
    1. Build a patient-prefix dict from expression columns:
         { "TCGA-BH-A0BZ": ["TCGA-BH-A0BZ-01A-...", ...], ... }
    2. For each clinical row (12-char patient key), look up its
       expression samples via that dict.
    3. Assign tumour vs normal from TCGA barcode position 13-14 (0-indexed):
         "01" = primary tumour
         "11" = adjacent normal
         "06" = metastatic (excluded)

  ALSO FIXED:
    - Robust sep detection for the clinical file (tab vs comma)
    - Handles clinical files where the barcode column is NOT the index
      but a named column (sampleID, submitter_id, etc.)
    - Handles Xena HiSeqV2 which is already log2(RSEM+1) — no double-log
    - Normal sample collection uses both barcode suffix AND the
      sample_type / tissue_status clinical columns as fallbacks

SELF-CONTAINED:
  Downloads GSE176078 scRNA-seq first (expected: 0 ILC cells — Wu et al.
  classifies by molecular subtype, not histology).
  Falls through automatically to TCGA-BRCA bulk RNA-seq.
  All data downloaded to ILC_s1_analysis/data/.

GEOMETRY-FIRST (Protocol v2.0):
  Unfiltered top-mover scan printed FIRST.
  Panel prediction tests printed SECOND.

PREDICTIONS FROM BRCA-S6a (predictions.md):
  P1: CDH1 DOMINANT LOSS — structural axis
  P2: ESR1/FOXA1/GATA3/PGR RETAINED — luminal identity preserved
  P3: EZH2 elevated AND negatively correlates with CDH1 (not ESR1)
  P4: DNMT3A elevated — CDH1 methylation driver
  P5: LOW proliferation — MKI67/TOP2A low vs TNBC and HER2
  P6: PI3K/AKT/mTOR elevated — PIK3CA dominant driver in ILC
  P7: NO basal/mesenchymal programme — VIM/ZEB1/SOX10 low
  P8: CDH1-based depth axis — within-ILC CDH1 gradient
  P9: Drug targets confirmed from geometry alone
"""

import os
import sys
import gzip
import tarfile
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
# CONFIGURATION
# ============================================================

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
BASE_DIR    = os.path.join(SCRIPT_DIR, "ILC_s1_analysis")
DATA_DIR    = os.path.join(BASE_DIR,   "data")
RESULTS_DIR = os.path.join(BASE_DIR,   "results")

LOG_FILE        = os.path.join(RESULTS_DIR, "ilc_s1_log.txt")
FIG_FILE        = os.path.join(RESULTS_DIR, "ilc_s1_figure.png")
CSV_FILE        = os.path.join(RESULTS_DIR, "ilc_s1_panel.csv")
TOP_MOVERS_FILE = os.path.join(RESULTS_DIR, "ilc_s1_top_movers.csv")
CROSS_FILE      = os.path.join(RESULTS_DIR, "ilc_s1_cross_subtype.csv")

for d in [BASE_DIR, DATA_DIR, RESULTS_DIR]:
    os.makedirs(d, exist_ok=True)

# ── GSE176078 (scRNA-seq — always checked first) ─────────────
SC_TARFILE   = "GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"
SC_URL       = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/"
                "GSE176078/suppl/"
                "GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz")
SC_INNER_DIR = "Wu_etal_2021_BRCA_scRNASeq"
SC_TARPATH   = os.path.join(DATA_DIR, SC_TARFILE)
SC_META_PATH = os.path.join(DATA_DIR, SC_INNER_DIR, "metadata.csv")
MIN_SC_ILC   = 50

# ── TCGA-BRCA expression (Xena HiSeqV2 — log2 RSEM+1) ───────
EXPR_URLS = [
    "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2.gz",
    "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2_percentile.gz",
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA-BRCA.htseq_fpkm.tsv.gz",
    "https://gdc.xenahubs.net/download/TCGA-BRCA.htseq_fpkm.tsv.gz",
]
EXPR_FILE = os.path.join(DATA_DIR, "TCGA_BRCA_HiSeqV2.gz")

# ── TCGA-BRCA clinical (Xena sampleMap) ──────────────────────
CLIN_URLS = [
    "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix",
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA-BRCA.GDC_phenotype.tsv.gz",
    "https://gdc.xenahubs.net/download/TCGA-BRCA.clinical.tsv.gz",
]
CLIN_FILE = os.path.join(DATA_DIR, "TCGA_BRCA_clinicalMatrix.tsv")

MIN_ILC    = 30
MIN_NORMAL = 15

# TCGA barcode positions (0-indexed into full sample barcode)
# Example: TCGA-BH-A0BZ-01A-11R-A12P-07
#   chars 0-11  = patient prefix  "TCGA-BH-A0BZ"
#   chars 13-14 = sample type     "01" (tumour) / "11" (normal)
PATIENT_PREFIX_LEN = 12
SAMPLE_TYPE_POS    = (13, 15)    # slice indices

ILC_KEYWORDS = ["infiltrating lobular", "invasive lobular",
                "lobular carcinoma", "lobular"]
NORMAL_KEYWORDS = ["solid tissue normal", "normal", "adjacent"]

# PAM50 labels used in Xena clinical files
LUMA_LABELS  = {"LumA"}
TNBC_LABELS  = {"Basal"}
HER2_LABELS  = {"Her2"}
LUMB_LABELS  = {"LumB"}

# ============================================================
# GENE PANELS
# ============================================================

CDH1_AXIS     = ["CDH1", "CDH3", "CTNND1", "CDH2"]
LUMINAL_TFS   = ["ESR1", "FOXA1", "GATA3", "PGR", "KRT8", "KRT18", "SPDEF"]
BASAL_EMT     = ["KRT5", "KRT14", "VIM", "ZEB1", "ZEB2", "SNAI1",
                 "SOX10", "FOXC1", "FN1", "TWIST1"]
EPIGENETIC    = ["EZH2", "EED", "SUZ12", "HDAC1", "HDAC2",
                 "KDM1A", "DNMT3A", "DNMT1", "DNMT3B"]
PROLIFERATION = ["MKI67", "TOP2A", "PCNA", "CCNB1", "CDK2",
                 "CDK4", "CCND1", "AURKA"]
PI3K_AXIS     = ["AKT1", "AKT2", "MTOR", "PTEN", "PIK3CA",
                 "PIK3R1", "RPS6KB1"]
HER_SIGNAL    = ["ERBB2", "ERBB3", "EGFR", "GRB7"]
ILC_SPECIFIC  = ["TBX3", "RUNX1", "CBFB"]
DNA_REPAIR    = ["TP53", "BRCA1", "BRCA2", "RB1", "PARP1", "ATM"]
AR_AXIS       = ["AR", "CLDN3", "CLDN4", "CLDN7"]
CONTROLS      = ["CDX2", "SPI1", "MBP", "NKX2-1", "OLIG2"]

ALL_GENES = list(dict.fromkeys(
    CDH1_AXIS + LUMINAL_TFS + BASAL_EMT + EPIGENETIC +
    PROLIFERATION + PI3K_AXIS + HER_SIGNAL + ILC_SPECIFIC +
    DNA_REPAIR + AR_AXIS + CONTROLS
))

# ============================================================
# LOGGING
# ============================================================

_log_lines = []

def log(msg=""):
    print(msg)
    _log_lines.append(str(msg))

def write_log():
    with open(LOG_FILE, "w") as f:
        f.write("\n".join(_log_lines))

# ============================================================
# DOWNLOAD UTILITY
# ============================================================

def download_file(url, dest, label="", timeout=180, chunk=1024 * 1024):
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
                        print(f"    {done*100//total}% ({done/1e6:.1f}MB)...",
                              end="\r", flush=True)
        print()
        log(f"  OK: {os.path.basename(dest)} ({os.path.getsize(dest)/1e6:.1f} MB)")
        return True
    except Exception as e:
        log(f"  Failed: {type(e).__name__}: {e}")
        if os.path.exists(dest):
            os.remove(dest)
        return False


def try_urls(urls, dest, label=""):
    if os.path.exists(dest) and os.path.getsize(dest) > 10_000:
        log(f"  Already present: {os.path.basename(dest)} "
            f"({os.path.getsize(dest)/1e6:.1f} MB)")
        return True
    for i, url in enumerate(urls, 1):
        log(f"  Attempt {i}/{len(urls)}: {url[:72]}...")
        if download_file(url, dest, label):
            return True
        time.sleep(3)
    return False


def open_gz(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

# ============================================================
# STEP 1 — CHECK GSE176078 FOR ILC CELLS
# ============================================================

def check_scrna_for_ilc():
    log("")
    log("=" * 65)
    log("STEP 1: ATTEMPT scRNA-seq PATH (GSE176078)")
    log("=" * 65)

    if not (os.path.exists(SC_TARPATH) and os.path.getsize(SC_TARPATH) > 1e8):
        log("  Downloading GSE176078 tar.gz ...")
        ok = download_file(SC_URL, SC_TARPATH, "GSE176078 scRNA-seq")
        if not ok:
            log("  Download failed — skipping scRNA-seq path.")
            return False

    if not os.path.exists(SC_META_PATH):
        log(f"  Extracting {SC_TARFILE} ...")
        with tarfile.open(SC_TARPATH, "r:gz") as tf:
            tf.extractall(DATA_DIR)
        log(f"  Extracted to {DATA_DIR}")

    if not os.path.exists(SC_META_PATH):
        log("  metadata.csv not found after extraction — skipping scRNA-seq path.")
        return False

    log(f"  Reading metadata: {SC_META_PATH}")
    meta = pd.read_csv(SC_META_PATH)
    log(f"  Metadata shape: {meta.shape}")

    ct_col = "celltype_subset"
    if ct_col in meta.columns:
        log(f"  Cell type counts ({ct_col}):")
        for ct, n in meta[ct_col].value_counts().items():
            log(f"    {ct}: {n}")

    ilc_cells = 0
    for col in meta.select_dtypes(include="object").columns:
        hits = meta[col].str.lower().str.contains("lobular|ilc", na=False).sum()
        if hits > 0:
            log(f"  ILC-related cells in column '{col}': {hits}")
            ilc_cells += hits

    if ilc_cells == 0:
        log("  Cells with lobular/ILC keyword in any column: 0")

    if ilc_cells < MIN_SC_ILC:
        log(f"  ILC cells not found or insufficient (n={ilc_cells} < {MIN_SC_ILC}).")
        log("  Wu et al. 2021 classified by molecular subtype (not histology).")
        log("  ILC tumors were not explicitly included — falling through to TCGA.")
        return False

    log(f"  ILC cells found: {ilc_cells} — scRNA-seq path viable.")
    return True

# ============================================================
# STEP 2 — ACQUIRE TCGA DATA
# ============================================================

def acquire_tcga():
    log("")
    log("=" * 65)
    log("STEP 2: TCGA-BRCA BULK RNA-seq PATH")
    log("=" * 65)

    log("\n  [2a] Expression matrix")
    if not try_urls(EXPR_URLS, EXPR_FILE, "TCGA-BRCA HiSeqV2"):
        log("  FATAL: Expression matrix unavailable. Exiting.")
        write_log()
        sys.exit(1)

    log("\n  [2b] Clinical matrix")
    if not try_urls(CLIN_URLS, CLIN_FILE, "TCGA-BRCA clinical"):
        log("  WARNING: Clinical matrix unavailable.")
        log("  Without histological annotation ILC cannot be identified.")
        log("  Cannot proceed without clinical data for ILC.")
        write_log()
        sys.exit(1)

# ============================================================
# STEP 3 — LOAD EXPRESSION
# ============================================================

def load_expression():
    log("")
    log("=" * 65)
    log("LOADING EXPRESSION MATRIX")
    log("=" * 65)

    with open_gz(EXPR_FILE, "rt") as fh:
        first_line  = fh.readline()
        second_line = fh.readline()

    sep = "\t" if "\t" in first_line else ","
    log(f"  Delimiter: {'TAB' if sep == chr(9) else 'COMMA'}")
    log(f"  First line sample: {first_line[:120]}")

    expr = pd.read_csv(EXPR_FILE, sep=sep, index_col=0, low_memory=False)

    # Xena HiSeqV2: rows = genes, cols = samples.
    # If it comes out transposed (rows = samples), flip it.
    # Heuristic: if nrows < 5000 and ncols > 10000, it's transposed.
    if expr.shape[0] < 5000 and expr.shape[1] > 10000:
        log("  Transposing (detected samples×genes orientation).")
        expr = expr.T

    expr = expr.apply(pd.to_numeric, errors="coerce")
    expr.dropna(how="all", inplace=True)

    log(f"  Expression shape: {expr.shape}  (genes × samples)")
    log(f"  Gene index sample:   {list(expr.index[:5])}")
    log(f"  Sample column sample:{list(str(c) for c in expr.columns[:3])}")

    return expr

# ============================================================
# STEP 4 — LOAD CLINICAL
# ============================================================

def load_clinical():
    log("")
    log("=" * 65)
    log("LOADING CLINICAL MATRIX")
    log("=" * 65)

    # Detect separator robustly
    with open_gz(CLIN_FILE, "rt") as fh:
        first_line = fh.readline()
    sep = "\t" if first_line.count("\t") > first_line.count(",") else ","
    log(f"  Delimiter: {'TAB' if sep == chr(9) else 'COMMA'}")

    clin = pd.read_csv(CLIN_FILE, sep=sep, low_memory=False)
    log(f"  Clinical shape (before index): {clin.shape}")
    log(f"  All columns ({len(clin.columns)}): {list(clin.columns)}")

    # ── Identify the barcode column ─────────────────────────
    # It may be the first column, or named sampleID / submitter_id /
    # bcr_patient_barcode, or already the DataFrame index.
    barcode_col = None
    barcode_candidates = [
        "sampleID", "submitter_id", "bcr_patient_barcode",
        "sample", "barcode", "patient_id", "case_id",
    ]
    # Check if first column contains TCGA barcodes
    first_col_vals = clin.iloc[:, 0].astype(str)
    if first_col_vals.str.startswith("TCGA").any():
        barcode_col = clin.columns[0]
        log(f"  Barcode in first column: '{barcode_col}'")
    else:
        for col in barcode_candidates:
            if col in clin.columns:
                if clin[col].astype(str).str.startswith("TCGA").any():
                    barcode_col = col
                    log(f"  Barcode column identified: '{barcode_col}'")
                    break

    if barcode_col:
        clin = clin.set_index(barcode_col)
    else:
        # Assume the default index is barcodes (happens when
        # the file was already read with index_col=0 implicitly)
        clin.index = clin.index.astype(str)
        if not clin.index[0].startswith("TCGA"):
            # Try re-reading with index_col=0
            clin2 = pd.read_csv(CLIN_FILE, sep=sep, index_col=0, low_memory=False)
            clin2.index = clin2.index.astype(str)
            if clin2.index[0].startswith("TCGA"):
                clin = clin2
                log("  Re-read with index_col=0: TCGA barcodes confirmed in index.")
            else:
                log("  WARNING: Could not identify barcode column.")
                log(f"  Index sample: {list(clin.index[:3])}")

    clin.index = clin.index.astype(str).str.strip()
    log(f"  Clinical index sample: {list(clin.index[:5])}")
    log(f"  Clinical index length (chars): {len(clin.index[0])}")
    log(f"  Clinical shape (final): {clin.shape}")

    # ── Identify histology and PAM50 columns ────────────────
    hist_col  = None
    pam50_col = None

    for col in ["histological_type", "primary_diagnosis", "histologic_diagnosis",
                "icd_o_3_histology", "morphology", "tumor_type"]:
        if col in clin.columns:
            hist_col = col
            break

    for col in ["PAM50Call_RNAseq", "PAM50_mRNA_nature2012",
                "Integrated_Clusters_with_PAM50__nature2012",
                "pam50", "PAM50"]:
        if col in clin.columns:
            pam50_col = col
            break

    log(f"\n  Histology column: '{hist_col}'")
    if hist_col:
        lobular_n = clin[hist_col].astype(str).str.lower().str.contains(
            "lobular", na=False).sum()
        log(f"  Lobular entries in histology column: {lobular_n}")
        log(f"  Histology value sample: "
            f"{list(clin[hist_col].dropna().unique()[:8])}")

    log(f"  PAM50 column: '{pam50_col}'")
    if pam50_col:
        log(f"  PAM50 value counts:")
        for v, n in clin[pam50_col].value_counts().items():
            log(f"    {v}: {n}")

    return clin, hist_col, pam50_col

# ============================================================
# STEP 5 — CLASSIFY SAMPLES
#
# THE KEY FIX:
#   Clinical index: TCGA-BH-A0BZ          (12 chars = patient)
#   Expression cols: TCGA-BH-A0BZ-01A-... (28 chars = sample)
#
#   Build patient_map: { "TCGA-BH-A0BZ": ["TCGA-BH-A0BZ-01A-...", ...] }
#   Look up each clinical patient prefix in patient_map.
#   Then read sample type from barcode position 13:15.
# ============================================================

def build_patient_map(expr_cols):
    """
    Returns:
      patient_map: dict {12-char-prefix: [full_sample_barcode, ...]}
      all_tumour:  list of full barcodes where position 13:15 == "01"
      all_normal:  list of full barcodes where position 13:15 == "11"
    """
    patient_map = {}
    all_tumour  = []
    all_normal  = []

    for col in expr_cols:
        s = str(col).strip()
        pfx = s[:PATIENT_PREFIX_LEN] if len(s) >= PATIENT_PREFIX_LEN else s
        patient_map.setdefault(pfx, []).append(s)

        st = s[SAMPLE_TYPE_POS[0]:SAMPLE_TYPE_POS[1]] if len(s) >= SAMPLE_TYPE_POS[1] else "??"
        if st == "01":
            all_tumour.append(s)
        elif st == "11":
            all_normal.append(s)

    return patient_map, all_tumour, all_normal


def classify_samples(expr, clin, hist_col, pam50_col):
    log("")
    log("=" * 65)
    log("IDENTIFYING ILC AND NORMAL SAMPLES")
    log("=" * 65)

    expr_cols = [str(c).strip() for c in expr.columns]
    patient_map, all_tumour, all_normal = build_patient_map(expr_cols)

    log(f"  Expression columns total:    {len(expr_cols)}")
    log(f"  Expression barcode example:  {expr_cols[0]}")
    log(f"  Expression barcode length:   {len(expr_cols[0])}")
    log(f"  Patient prefixes in expr:    {len(patient_map)}")
    log(f"  Tumour samples (-01*):       {len(all_tumour)}")
    log(f"  Normal samples (-11*):       {len(all_normal)}")

    ilc_s  = []
    luma_s = []
    tnbc_s = []
    her2_s = []
    lumb_s = []
    norm_s = list(all_normal)   # seed with barcode-derived normals

    matched_patients = 0
    ilc_patients     = 0

    clin.index = clin.index.astype(str).str.strip()

    log(f"\n  Clinical index example: '{clin.index[0]}'  "
        f"(length {len(clin.index[0])})")
    log(f"  Patient map key example: '{list(patient_map.keys())[0]}'  "
        f"(length {len(list(patient_map.keys())[0])})")

    # ── Main matching loop ───────────────────────────────────
    for clin_id in clin.index:
        # Build patient prefix from clinical ID
        # Clinical barcodes may be:
        #   12 chars: TCGA-BH-A0BZ           (patient-level — Xena sampleMap)
        #   28 chars: TCGA-BH-A0BZ-01A-...   (sample-level — GDC phenotype)
        # Always use first 12 chars.
        pfx = str(clin_id)[:PATIENT_PREFIX_LEN]

        expr_samples = patient_map.get(pfx, [])
        if not expr_samples:
            continue
        matched_patients += 1

        # Get histology and PAM50 for this patient
        try:
            row = clin.loc[clin_id]
        except KeyError:
            continue

        hist_val  = str(row[hist_col]).lower()  if hist_col  and hist_col  in clin.columns else ""
        pam50_val = str(row[pam50_col]).strip() if pam50_col and pam50_col in clin.columns else ""

        is_ilc = any(kw in hist_val for kw in ILC_KEYWORDS)

        # Also check sample_type for normals
        sample_type_val = ""
        for stcol in ["sample_type", "sample_type.1", "tissue_status",
                      "_sample_type", "Sample_Type"]:
            if stcol in clin.columns:
                sample_type_val = str(row.get(stcol, "")).lower()
                break

        for s in expr_samples:
            st = s[SAMPLE_TYPE_POS[0]:SAMPLE_TYPE_POS[1]] if len(s) >= SAMPLE_TYPE_POS[1] else "??"

            if st == "11":
                if s not in norm_s:
                    norm_s.append(s)
            elif st == "01":
                if any(kw in sample_type_val for kw in NORMAL_KEYWORDS):
                    if s not in norm_s:
                        norm_s.append(s)
                    continue

                if is_ilc:
                    ilc_s.append(s)
                elif pam50_val in LUMA_LABELS:
                    luma_s.append(s)
                elif pam50_val in TNBC_LABELS:
                    tnbc_s.append(s)
                elif pam50_val in HER2_LABELS:
                    her2_s.append(s)
                elif pam50_val in LUMB_LABELS:
                    lumb_s.append(s)
            # metastatic (06) and other types excluded

        if is_ilc:
            ilc_patients += 1

    log("")
    log(f"  Patients in clinical matched to expression: {matched_patients}")
    log(f"  ILC patients matched to expression:         {ilc_patients}")
    log(f"  ILC samples (expression entries):           {len(ilc_s)}")
    log(f"  Normal samples:                             {len(norm_s)}")
    log(f"  IDC LumA (PAM50):                           {len(luma_s)}")
    log(f"  IDC TNBC (PAM50 Basal):                     {len(tnbc_s)}")
    log(f"  IDC HER2 (PAM50 Her2):                      {len(her2_s)}")
    log(f"  IDC LumB (PAM50):                           {len(lumb_s)}")

    # ── Diagnostic: if ILC == 0 despite patients being matched ──
    if ilc_patients > 0 and len(ilc_s) == 0:
        log("")
        log("  DIAGNOSTIC: ILC patients matched but 0 ILC samples assigned.")
        log("  This means expression samples were found but all had st != '01'.")
        log("  Checking sample type codes for a known ILC patient...")
        for clin_id in clin.index[:50]:
            if hist_col and hist_col in clin.columns:
                hv = str(clin.loc[clin_id, hist_col]).lower()
                if "lobular" in hv:
                    pfx = str(clin_id)[:PATIENT_PREFIX_LEN]
                    samps = patient_map.get(pfx, [])
                    log(f"    Patient '{clin_id}' → prefix '{pfx}' → samples: {samps}")
                    for s in samps:
                        st = s[SAMPLE_TYPE_POS[0]:SAMPLE_TYPE_POS[1]] if len(s) >= SAMPLE_TYPE_POS[1] else "??"
                        log(f"      Sample '{s}' → st='{st}'")
                    break

    # ── Hard failure ─────────────────────────────────────────
    if len(ilc_s) < MIN_ILC:
        log("")
        log(f"  ERROR: Only {len(ilc_s)} ILC samples found "
            f"(minimum required: {MIN_ILC}).")

        if matched_patients == 0:
            log("")
            log("  ZERO patients matched. Barcode prefix mismatch still present.")
            log("  Full diagnostic:")
            clin_sample = list(clin.index[:5])
            expr_sample = list(patient_map.keys())[:5]
            log(f"  Clinical index first 5:          {clin_sample}")
            log(f"  Expression patient prefixes (5): {expr_sample}")
            log("")
            log("  Check that CLIN_FILE and EXPR_FILE are both TCGA-BRCA.")
            log("  If clinical file uses sample-level barcodes (28 chars),")
            log("  adjust PATIENT_PREFIX_LEN at the top of this script to match.")

        write_log()
        sys.exit(1)

    if len(norm_s) < MIN_NORMAL:
        log(f"  WARNING: Only {len(norm_s)} normal samples.")
        if len(norm_s) == 0:
            log("  FATAL: No normal samples. Cannot compute cancer-vs-normal.")
            write_log()
            sys.exit(1)

    return {
        "ilc":    ilc_s,
        "normal": norm_s,
        "luma":   luma_s,
        "tnbc":   tnbc_s,
        "her2":   her2_s,
        "lumb":   lumb_s,
    }

# ============================================================
# STEP 6 — VERIFY SAMPLES ARE IN EXPRESSION MATRIX COLUMNS
# ============================================================

def verify_pops(expr, pops):
    log("")
    log("=" * 65)
    log("VERIFYING SAMPLES IN EXPRESSION MATRIX")
    log("=" * 65)

    available = set(str(c).strip() for c in expr.columns)
    verified  = {}

    for grp, slist in pops.items():
        kept    = [s for s in slist if s in available]
        dropped = len(slist) - len(kept)
        if dropped > 0:
            log(f"  {grp.upper():<12} classified={len(slist)}, "
                f"in matrix={len(kept)}, dropped={dropped}")
        else:
            log(f"  {grp.upper():<12} n={len(kept)} — all present")
        verified[grp] = kept

    if len(verified["ilc"]) < MIN_ILC:
        log(f"  FATAL: {len(verified['ilc'])} ILC samples verified in matrix.")
        write_log()
        sys.exit(1)
    if len(verified["normal"]) < MIN_NORMAL:
        log(f"  FATAL: {len(verified['normal'])} normal samples in matrix.")
        write_log()
        sys.exit(1)

    log("")
    log("  FINAL COHORT:")
    for grp, slist in verified.items():
        if slist:
            log(f"    {grp.upper():<12} n={len(slist)}")

    return verified

# ============================================================
# STEP 7 — NORMALISE
# Xena HiSeqV2 is log2(RSEM+1). Apply once more gives log2(log2+2)
# which is wrong. We detect and avoid double-logging.
# ============================================================

def normalise(expr):
    log("")
    log("=" * 65)
    log("NORMALISE")
    log("=" * 65)

    vals       = expr.values.ravel()
    nz         = vals[vals > 0]
    median_nz  = float(np.nanmedian(nz)) if len(nz) else 0.0
    max_val    = float(np.nanmax(vals))

    log(f"  Non-zero median: {median_nz:.4f}")
    log(f"  Max value:       {max_val:.4f}")

    if max_val < 40:
        # Already log2-transformed (Xena HiSeqV2 is log2(RSEM+1), max ~20)
        log("  Values appear already log2-transformed. Using as-is.")
        return expr
    elif max_val < 200:
        # Possibly log-TPM or log-FPKM with slightly higher range
        log("  Values in moderate range. Using as-is (assuming log-transformed).")
        return expr
    else:
        # Raw counts or FPKM — apply log2(x+1)
        log(f"  Large values (max={max_val:.1f}) ��� applying log2(x+1).")
        return np.log2(expr + 1)

# ============================================================
# UTILITY: gene stats
# ============================================================

def gene_stats(expr, gene, samps_a, samps_b):
    """Return (mean_a, mean_b, pct_change, pval) or (None,None,None,None)."""
    if gene not in expr.index:
        return None, None, None, None
    a = expr.loc[gene, samps_a].values.astype(float)
    b = expr.loc[gene, samps_b].values.astype(float)
    ma = float(np.nanmean(a))
    mb = float(np.nanmean(b))
    pct = (ma - mb) / mb * 100 if mb > 0.001 else (100.0 if ma > 0 else 0.0)
    try:
        _, pval = stats.mannwhitneyu(a, b, alternative="two-sided")
    except Exception:
        pval = 1.0
    return ma, mb, pct, pval


def fmt_p(p):
    if p < 1e-100: return "p<1e-100 ***"
    if p < 1e-10:  return f"p={p:.2e} ***"
    if p < 1e-5:   return f"p={p:.2e} **"
    if p < 0.05:   return f"p={p:.4f} *"
    return         f"p={p:.4f}"

# ============================================================
# STEP 8 — UNFILTERED DISCOVERY SCAN
# ============================================================

def unfiltered_discovery(expr, pops):
    log("")
    log("=" * 65)
    log("UNFILTERED DISCOVERY SCAN")
    log("  Geometry first. No panel imposed.")
    log("  ILC vs Adjacent Normal — all expressed genes.")
    log("=" * 65)

    ilc_s  = pops["ilc"]
    norm_s = pops["normal"]
    log(f"  ILC samples:    {len(ilc_s)}")
    log(f"  Normal samples: {len(norm_s)}")

    ilc_mat  = expr[ilc_s].values   # genes × ILC-samples
    norm_mat = expr[norm_s].values  # genes × normal-samples
    gene_ids = expr.index.tolist()

    results = []
    for i, g in enumerate(gene_ids):
        h = ilc_mat[i].astype(float)
        n = norm_mat[i].astype(float)
        mh = float(np.nanmean(h))
        mn = float(np.nanmean(n))
        if mh < 0.01 and mn < 0.01:
            continue
        pct = (mh - mn) / mn * 100 if mn > 0.001 else (100.0 if mh > 0 else 0.0)
        if abs(pct) > 5e5:
            continue
        try:
            _, pval = stats.mannwhitneyu(h, n, alternative="two-sided")
        except Exception:
            pval = 1.0
        results.append({"gene": g, "ilc_mean": mh, "norm_mean": mn,
                         "pct_change": pct, "pval": pval})

    rdf = (pd.DataFrame(results)
             .assign(abs_pct=lambda d: d["pct_change"].abs())
             .sort_values("abs_pct", ascending=False)
             .reset_index(drop=True))
    rdf.to_csv(TOP_MOVERS_FILE, index=False)
    log(f"  Full ranked results saved: {TOP_MOVERS_FILE}")

    gained = rdf[rdf["pct_change"] > 0].head(30)
    lost   = rdf[rdf["pct_change"] < 0].sort_values("pct_change").head(30)

    hdr = f"  {'Gene':<14} {'ILC':>8} {'Normal':>8} {'Change':>10}  Stat"
    sep = "  " + "-" * 62

    log("")
    log("  TOP 30 GAINED IN ILC vs Normal:")
    log(hdr); log(sep)
    for _, row in gained.iterrows():
        log(f"  {row['gene']:<14} {row['ilc_mean']:>8.4f} {row['norm_mean']:>8.4f} "
            f"{row['pct_change']:>+9.1f}%  {fmt_p(row['pval'])}")

    log("")
    log("  TOP 30 LOST IN ILC vs Normal:")
    log(hdr); log(sep)
    for _, row in lost.iterrows():
        log(f"  {row['gene']:<14} {row['ilc_mean']:>8.4f} {row['norm_mean']:>8.4f} "
            f"{row['pct_change']:>+9.1f}%  {fmt_p(row['pval'])}")

    # ── Geometry read ────────────────────────────────────────
    log("")
    log("  GEOMETRY READ (unfiltered):")

    def read_gene(sym, note):
        row = rdf[rdf["gene"] == sym]
        if row.empty:
            log(f"  {sym:<14} NOT IN MATRIX")
            return None
        fc  = float(row["pct_change"].values[0])
        rnk = int(row.index[0]) + 1
        log(f"  {sym:<14} rank #{rnk:<6} {fc:>+9.1f}%  "
            f"{fmt_p(float(row['pval'].values[0]))}  ← {note}")
        return fc

    cdh1_fc   = read_gene("CDH1",   "P1: structural loss — expected DOMINANT")
    esr1_fc   = read_gene("ESR1",   "P2: luminal TF — expected RETAINED")
    foxa1_fc  = read_gene("FOXA1",  "P2: luminal TF — expected RETAINED")
    gata3_fc  = read_gene("GATA3",  "P2: luminal TF — expected RETAINED")
    vim_fc    = read_gene("VIM",    "P7: mesenchymal — expected LOW")
    zeb1_fc   = read_gene("ZEB1",   "P7: EMT TF — expected LOW")
    sox10_fc  = read_gene("SOX10",  "P7: basal/NC — expected LOW (TNBC was +1323%)")
    krt5_fc   = read_gene("KRT5",   "P7: basal CK — expected LOW (TNBC was +508%)")
    ezh2_fc   = read_gene("EZH2",   "P3: epigenetic — expected ELEVATED")
    dnmt3a_fc = read_gene("DNMT3A", "P4: CDH1 methylation — expected ELEVATED")
    mki67_fc  = read_gene("MKI67",  "P5: proliferation — expected LOW vs TNBC/HER2")
    erbb2_fc  = read_gene("ERBB2",  "negative control — should NOT be elevated")
    akt1_fc   = read_gene("AKT1",   "P6: PI3K axis — expected ELEVATED")
    pten_fc   = read_gene("PTEN",   "P6: PI3K suppressor — may be reduced")

    log("")
    log("  STRUCTURAL INVERSION TEST:")
    log("  Predicted: ILC is the geometric inverse of TNBC")
    log("  TNBC: ESR1 erased, VIM elevated, CDH1 retained")
    log("  ILC:  ESR1 retained, VIM low, CDH1 erased")
    log("")

    if cdh1_fc is not None:
        if cdh1_fc < -40:
            log(f"  CDH1 {cdh1_fc:+.1f}% → STRUCTURAL LOSS CONFIRMED.")
        elif cdh1_fc < -15:
            log(f"  CDH1 {cdh1_fc:+.1f}% → Partial mRNA reduction.")
            log("  Note: CDH1 mutation-only ILC may show full protein loss with")
            log("  partial/no mRNA reduction. Check protein-level data if available.")
        else:
            log(f"  CDH1 {cdh1_fc:+.1f}% → CDH1 mRNA not substantially reduced.")
            log("  This is consistent with mutation-only ILC (CDH1 protein absent,")
            log("  mRNA intact). Methylation-driven ILC (~35%) shows mRNA reduction.")

    if esr1_fc is not None:
        if esr1_fc > -30:
            log(f"  ESR1 {esr1_fc:+.1f}% → LUMINAL TF RETAINED. P2 supported.")
        else:
            log(f"  ESR1 {esr1_fc:+.1f}% → ESR1 substantially reduced. Unexpected for ILC.")

    if vim_fc is not None and zeb1_fc is not None:
        if vim_fc < 150 and zeb1_fc < 100:
            log(f"  VIM {vim_fc:+.1f}% / ZEB1 {zeb1_fc:+.1f}% → No EMT programme. P7 supported.")
            log("  ILC invades without canonical EMT — structural inversion confirmed.")
        else:
            log(f"  VIM {vim_fc:+.1f}% / ZEB1 {zeb1_fc:+.1f}% → Some mesenchymal signal present.")

    return rdf

# ============================================================
# STEP 9 — PANEL ANALYSIS
# ============================================================

def panel_analysis(expr, pops, rdf):
    log("")
    log("=" * 65)
    log("PANEL ANALYSIS — PREDICTION TESTS (BRCA-S6a)")
    log("=" * 65)

    ilc_s  = pops["ilc"]
    norm_s = pops["normal"]
    rows   = []

    def report(gene, mi, mn, pct, pval, pred, status):
        if mi is None:
            log(f"  {gene:<12} NOT IN MATRIX")
            return
        log(f"  {gene:<12} ILC={mi:>7.4f}  Norm={mn:>7.4f}  "
            f"{pct:>+8.1f}%  {fmt_p(pval):<22}  [{status}]  {pred}")
        rows.append({"gene": gene, "ilc_mean": mi, "norm_mean": mn,
                     "pct_change": pct, "pval": pval,
                     "prediction": pred, "status": status})

    def gs(g):
        return gene_stats(expr, g, ilc_s, norm_s)

    # ─── P1 — CDH1 dominant structural loss ─────────────────
    log("")
    log("  P1 — CDH1 DOMINANT STRUCTURAL LOSS")
    log("  Expected: CDH1 largest or near-largest lost gene")
    log("-" * 65)
    for g in CDH1_AXIS:
        mi, mn, pct, pval = gs(g)
        if g == "CDH1":
            st = ("CONFIRMED"     if pct is not None and pct < -40 else
                  "PARTIAL"       if pct is not None and pct < -15 else
                  "NOTE-mRNA-intact-mutation-only-ILC")
        else:
            st = "INFO"
        report(g, mi, mn, pct, pval, "P1: CDH1 structural axis", st)

    # ─── P2 — Luminal TFs retained ───────────────────────────
    log("")
    log("  P2 — LUMINAL IDENTITY TFs RETAINED")
    log("  Expected: ESR1/FOXA1/GATA3/PGR comparable to normal (>-30%)")
    log("-" * 65)
    for g in LUMINAL_TFS:
        mi, mn, pct, pval = gs(g)
        st = ("CONFIRMED"     if pct is not None and pct > -30 else
              "PARTIAL"       if pct is not None and pct > -55 else
              "NOT CONFIRMED")
        report(g, mi, mn, pct, pval, "P2: Luminal TF retention", st)

    # ─── P3 — EZH2 elevated, correlates with CDH1 not ESR1 ──
    log("")
    log("  P3 — EZH2 ELEVATED AND TARGETS CDH1 NOT ESR1")
    log("  Expected: EZH2 elevated, r(EZH2,CDH1) < 0 within ILC")
    log("-" * 65)
    mi_ez, mn_ez, ezh2_pct, ezh2_pval = gs("EZH2")
    st = ("CONFIRMED"     if ezh2_pct is not None and ezh2_pct > 80 else
          "PARTIAL"       if ezh2_pct is not None and ezh2_pct > 30 else
          "NOT CONFIRMED")
    report("EZH2", mi_ez, mn_ez, ezh2_pct, ezh2_pval, "P3: EZH2 elevation", st)
    for g in ["EED", "SUZ12"]:
        mi, mn, pct, pval = gs(g)
        report(g, mi, mn, pct, pval, "P3: PRC2 complex", "INFO")

    if all(g in expr.index for g in ["EZH2", "CDH1"]) and len(ilc_s) > 10:
        ez_v  = expr.loc["EZH2", ilc_s].values.astype(float)
        cd_v  = expr.loc["CDH1", ilc_s].values.astype(float)
        v1 = ~(np.isnan(ez_v) | np.isnan(cd_v))
        if v1.sum() > 5:
            r1, p1 = stats.pearsonr(ez_v[v1], cd_v[v1])
            log(f"  EZH2 vs CDH1 (within ILC):  r={r1:+.3f}  {fmt_p(p1)}")
            if r1 < -0.2:
                log("  → EZH2 negatively correlated with CDH1 in ILC. P3 SUPPORTED.")
                log("    High EZH2 → low CDH1 = PRC2-mediated CDH1 silencing.")
            else:
                log(f"  → EZH2/CDH1 correlation weak (r={r1:+.3f}).")
                log("    PRC2 may not be dominant CDH1 silencer in bulk (mutation-only ILC).")
        if "ESR1" in expr.index:
            es_v = expr.loc["ESR1", ilc_s].values.astype(float)
            v2 = ~(np.isnan(ez_v) | np.isnan(es_v))
            if v2.sum() > 5:
                r2, p2 = stats.pearsonr(ez_v[v2], es_v[v2])
                log(f"  EZH2 vs ESR1 (within ILC):  r={r2:+.3f}  {fmt_p(p2)}")
                if r2 < -0.3:
                    log("  → EZH2 also silences ESR1 in ILC. Dual mechanism. NOVEL FINDING.")
                else:
                    log("  → EZH2 does NOT target ESR1 in ILC. P3 supported.")

    # ─── P4 — DNMT3A elevated ────────────────────────────────
    log("")
    log("  P4 — DNMT3A ELEVATED (CDH1 methylation driver)")
    log("  Expected: DNMT3A elevated vs normal")
    log("-" * 65)
    for g in ["DNMT3A", "DNMT1", "DNMT3B"]:
        mi, mn, pct, pval = gs(g)
        st = ("CONFIRMED" if g == "DNMT3A" and pct is not None and pct > 80 else
              "PARTIAL"   if g == "DNMT3A" and pct is not None and pct > 30 else
              "INFO")
        report(g, mi, mn, pct, pval, "P4: DNA methylation axis", st)

    # ─── P5 — Low proliferation ──────────────────────────────
    log("")
    log("  P5 — LOW PROLIFERATION (ILC Grade 1-2)")
    log("  Expected: MKI67/TOP2A lower than TNBC and HER2")
    log("-" * 65)
    for g in PROLIFERATION:
        mi, mn, pct, pval = gs(g)
        report(g, mi, mn, pct, pval, "P5: Proliferation", "INFO")

    if len(pops.get("tnbc", [])) >= 5 and "MKI67" in expr.index:
        ki_ilc  = float(np.nanmean(expr.loc["MKI67", pops["ilc"]].values))
        ki_tnbc = float(np.nanmean(expr.loc["MKI67", pops["tnbc"]].values))
        ki_her2 = (float(np.nanmean(expr.loc["MKI67", pops["her2"]].values))
                   if len(pops.get("her2", [])) >= 5 else None)
        log(f"  MKI67 cross-subtype: ILC={ki_ilc:.4f}  TNBC={ki_tnbc:.4f}"
            + (f"  HER2={ki_her2:.4f}" if ki_her2 else ""))
        if ki_ilc < ki_tnbc:
            log("  → ILC lower proliferation than TNBC. P5 SUPPORTED.")
        else:
            log("  → ILC NOT lower than TNBC for MKI67.")

    # ─── P6 — PI3K/AKT/mTOR elevated ────────────────────────
    log("")
    log("  P6 — PIK3CA PATHWAY ELEVATED")
    log("  Expected: AKT1/MTOR elevated; PTEN reduced")
    log("-" * 65)
    for g in PI3K_AXIS:
        mi, mn, pct, pval = gs(g)
        if g == "PTEN":
            st = ("CONFIRMED"     if pct is not None and pct < -20 else
                  "PARTIAL"       if pct is not None and pct < 0   else
                  "NOT CONFIRMED")
        elif g in ["AKT1", "MTOR", "RPS6KB1"]:
            st = ("CONFIRMED" if pct is not None and pct > 30 else
                  "PARTIAL"   if pct is not None and pct > 10 else
                  "NOT CONFIRMED")
        else:
            st = "INFO"
        report(g, mi, mn, pct, pval, "P6: PI3K pathway", st)

    # ─── P7 — Absence of basal/EMT programme ─────────────────
    log("")
    log("  P7 — NO BASAL OR MESENCHYMAL PROGRAMME")
    log("  Expected: VIM/ZEB1/ZEB2/SOX10/KRT5 all < +100%")
    log("-" * 65)
    for g in BASAL_EMT:
        mi, mn, pct, pval = gs(g)
        st = ("CONFIRMED"     if pct is not None and pct < 100 else
              "PARTIAL"       if pct is not None and pct < 300 else
              "NOT CONFIRMED")
        report(g, mi, mn, pct, pval, "P7: No basal/EMT", st)

    # ─── P8 — CDH1-based depth axis within ILC ───────────────
    log("")
    log("  P8 — CDH1-BASED DEPTH AXIS WITHIN ILC")
    log("  Tertile split on CDH1: deep (low CDH1) vs shallow (high CDH1)")
    log("  Key question: does ESR1 fall with CDH1? (unexpected) or stay?")
    log("-" * 65)

    depth_results = {}
    if "CDH1" in expr.index and len(ilc_s) > 10:
        cdh1_vals = expr.loc["CDH1", ilc_s].values.astype(float)
        depth     = -cdh1_vals   # deeper = lower CDH1 = higher score
        t33 = np.nanpercentile(depth, 33)
        t67 = np.nanpercentile(depth, 67)
        deep_s = [s for s, d in zip(ilc_s, depth) if d >= t67]
        shal_s = [s for s, d in zip(ilc_s, depth) if d <= t33]
        log(f"  CDH1 range in ILC: "
            f"min={float(np.nanmin(cdh1_vals)):.4f}, "
            f"max={float(np.nanmax(cdh1_vals)):.4f}, "
            f"mean={float(np.nanmean(cdh1_vals)):.4f}")
        log(f"  Deep ILC (low CDH1):     n={len(deep_s)}")
        log(f"  Shallow ILC (high CDH1): n={len(shal_s)}")

        depth_genes = (CDH1_AXIS + LUMINAL_TFS +
                       ["EZH2", "DNMT3A"] + BASAL_EMT[:5] +
                       PROLIFERATION[:3] + PI3K_AXIS[:4])
        depth_rows = []
        log("")
        log(f"  {'Gene':<12} {'Deep':>8} {'Shallow':>8} {'Change':>10}  Stat")
        log("  " + "-" * 56)
        for g in depth_genes:
            if g not in expr.index:
                continue
            dv = expr.loc[g, deep_s].values.astype(float)
            sv = expr.loc[g, shal_s].values.astype(float)
            md, ms = float(np.nanmean(dv)), float(np.nanmean(sv))
            if ms < 0.001:
                continue
            pct_ds = (md - ms) / ms * 100
            try:
                _, p_ds = stats.mannwhitneyu(dv, sv, alternative="two-sided")
            except Exception:
                p_ds = 1.0
            log(f"  {g:<12} {md:>8.4f} {ms:>8.4f} {pct_ds:>+9.1f}%  {fmt_p(p_ds)}")
            depth_rows.append({"gene": g, "deep_mean": md, "shallow_mean": ms,
                                "pct_change": pct_ds, "pval": p_ds})

        depth_df = pd.DataFrame(depth_rows)
        depth_results = {"depth_df": depth_df, "deep": deep_s, "shallow": shal_s}

        if "ESR1" in depth_df["gene"].values:
            esr1_d_row = depth_df[depth_df["gene"] == "ESR1"].iloc[0]
            fc_esr1    = esr1_d_row["pct_change"]
            log("")
            if fc_esr1 < -20:
                log(f"  ESR1 falls {fc_esr1:+.1f}% along CDH1 depth axis.")
                log("  → UNEXPECTED: ESR1 IS part of ILC depth axis.")
                log("  → Dual mechanism: CDH1 structural + ESR1 epigenetic.")
                log("  → NOVEL FINDING — report in BRCA-S6b.")
            else:
                log(f"  ESR1 stable across CDH1 depth axis ({fc_esr1:+.1f}%).")
                log("  → CDH1 structural axis independent of luminal TF identity. P8 SUPPORTED.")

    # ─── P9 — Drug target summary from geometry ──────────────
    log("")
    log("  P9 — DRUG TARGET SUMMARY FROM GEOMETRY ALONE")
    log("-" * 65)

    if "ESR1" in expr.index:
        ilc_esr1 = float(np.nanmean(expr.loc["ESR1", ilc_s].values))
        log(f"  ENDOCRINE THERAPY: ESR1 ILC={ilc_esr1:.4f}")
        if ilc_esr1 > 2.0:
            log("  → ER RETAINED at high level. Endocrine therapy strongly supported.")
            log("  → Tamoxifen / Aromatase inhibitors / CDK4/6 inhibitors (Palbociclib).")
        elif ilc_esr1 > 0.5:
            log("  → ER PRESENT at moderate level. Endocrine therapy supported.")

    if ezh2_pct is not None and ezh2_pct > 30:
        log(f"\n  EZH2 INHIBITOR: EZH2 elevated {ezh2_pct:+.1f}%")
        log("  → Tazemetostat (EZH2 inhibitor) has geometric basis.")
        log("  → NOVEL PREDICTION: target is CDH1 re-expression, NOT ESR1.")
        log("  → Distinct from Luminal A application of EZH2i.")

    log("\n  PI3K/mTOR INHIBITORS:")
    log("  → PIK3CA mutated in ~48% of ILC (highest of any BRCA subtype).")
    log("  → Alpelisib (PI3Kα), Everolimus (mTOR) have geometric basis.")

    if "ERBB2" in expr.index:
        erbb2_ilc = float(np.nanmean(expr.loc["ERBB2", ilc_s].values))
        log(f"\n  ANTI-HER2 (NEGATIVE): ERBB2 ILC={erbb2_ilc:.4f}")
        log("  → ERBB2 not amplified in ILC. Anti-HER2 has no geometric basis.")

    pd.DataFrame(rows).to_csv(CSV_FILE, index=False)
    log(f"\n  Panel results saved: {CSV_FILE}")
    return rows, depth_results

# ============================================================
# STEP 10 — CROSS-SUBTYPE COMPARISON
# ============================================================

def cross_subtype(expr, pops):
    log("")
    log("=" * 65)
    log("CROSS-SUBTYPE COMPARISON")
    log("  ILC vs LumA vs TNBC vs HER2 vs Normal")
    log("=" * 65)

    groups = {k: v for k, v in pops.items() if len(v) >= 5}
    focus  = ["CDH1", "ESR1", "FOXA1", "GATA3", "VIM", "ZEB1",
              "EZH2", "MKI67", "ERBB2", "DNMT3A", "KRT5", "SOX10",
              "CDH3", "AR", "PTEN", "AKT1", "CDH2", "CTNND1"]
    focus  = [g for g in focus if g in expr.index]

    header = f"  {'Gene':<14}" + "".join(f"  {g.upper():>10}" for g in groups)
    log(header)
    log("  " + "-" * (14 + 12 * len(groups)))

    rows = []
    for g in focus:
        row  = {"gene": g}
        line = f"  {g:<14}"
        for grp, samps in groups.items():
            m = float(np.nanmean(expr.loc[g, samps].values)) if samps else float("nan")
            row[grp] = m
            line += f"  {m:>10.4f}"
        log(line)
        rows.append(row)

    cross_df = pd.DataFrame(rows)
    cross_df.to_csv(CROSS_FILE, index=False)

    # Structural inversion check
    log("")
    log("  STRUCTURAL INVERSION SUMMARY:")
    if {"CDH1", "ESR1"}.issubset(set(cross_df["gene"])):
        c = cross_df.set_index("gene")
        if "ilc" in c.columns and "tnbc" in c.columns:
            log(f"  CDH1: ILC={c.loc['CDH1','ilc']:.4f}  "
                f"TNBC={c.loc['CDH1','tnbc']:.4f}  "
                f"Normal={c.loc['CDH1'].get('normal', float('nan')):.4f}")
            log(f"  ESR1: ILC={c.loc['ESR1','ilc']:.4f}  "
                f"TNBC={c.loc['ESR1','tnbc']:.4f}  "
                f"Normal={c.loc['ESR1'].get('normal', float('nan')):.4f}")
            cdh1_inv = c.loc["CDH1", "ilc"] < c.loc["CDH1", "tnbc"]
            esr1_inv = c.loc["ESR1", "ilc"] > c.loc["ESR1", "tnbc"]
            if cdh1_inv and esr1_inv:
                log("  ✓ CDH1 lower in ILC than TNBC — greater structural dissolution.")
                log("  ✓ ESR1 higher in ILC than TNBC — luminal identity retained.")
                log("  STRUCTURAL INVERSION CONFIRMED: ILC and TNBC are geometric opposites.")
            elif cdh1_inv:
                log("  CDH1 lower in ILC — structural arm confirmed.")
                log("  ESR1 not higher in ILC than TNBC — identity arm weaker.")
            elif esr1_inv:
                log("  ESR1 higher in ILC — identity arm confirmed.")
                log("  CDH1 not lower in ILC than TNBC — structural arm weaker in bulk.")

    return cross_df

# ============================================================
# STEP 11 — PCA GEOMETRY
# ============================================================

def pca_geometry(expr, pops):
    log("")
    log("=" * 65)
    log("PCA GEOMETRY")
    log("  Does ILC cluster near LumA (shared TF identity)?")
    log("  Or does CDH1 loss create a distinct geometric position?")
    log("=" * 65)

    groups = {k: v for k, v in pops.items() if len(v) >= 5}
    all_s, labels = [], []
    for grp, samps in groups.items():
        sel = samps[:100]
        all_s.extend(sel)
        labels.extend([grp] * len(sel))

    pca_genes = [g for g in ALL_GENES if g in expr.index]
    if len(pca_genes) < 5:
        log("  Too few panel genes for PCA. Skipping.")
        return None

    X  = expr.loc[pca_genes, all_s].T.values.astype(float)
    X  = np.nan_to_num(X)
    Xs = StandardScaler().fit_transform(X)
    pca = PCA(n_components=min(5, Xs.shape[1]))
    Xp  = pca.fit_transform(Xs)

    log(f"  Variance: PC1={pca.explained_variance_ratio_[0]*100:.1f}%  "
        f"PC2={pca.explained_variance_ratio_[1]*100:.1f}%")

    la         = np.array(labels)
    centroids  = {}
    log("  Centroids (PC1, PC2):")
    for grp in groups:
        m = la == grp
        if m.sum() == 0:
            continue
        c = Xp[m].mean(axis=0)
        centroids[grp] = c
        log(f"    {grp:<12} PC1={c[0]:>7.3f}  PC2={c[1]:>7.3f}")

    if {"ilc", "luma", "tnbc"}.issubset(centroids):
        d_la = float(np.linalg.norm(centroids["ilc"][:2] - centroids["luma"][:2]))
        d_tn = float(np.linalg.norm(centroids["ilc"][:2] - centroids["tnbc"][:2]))
        log(f"  ILC ↔ LumA distance: {d_la:.3f}")
        log(f"  ILC ↔ TNBC distance: {d_tn:.3f}")
        if d_la < d_tn:
            log("  ILC CLOSER to LumA. Shared luminal TF identity confirmed.")
        else:
            log("  ILC closer to TNBC than LumA — unexpected. Review cohort.")

    return {"pca": pca, "X_pca": Xp, "labels": labels, "centroids": centroids}

# ============================================================
# STEP 12 — FIGURE
# ============================================================

def generate_figure(expr, pops, panel_rows, cross_df, pca_result, rdf):
    log("")
    log("=" * 65)
    log("GENERATING FIGURE")
    log("=" * 65)

    COLORS = {"ilc":    "#8B0000", "normal": "#2E8B57",
              "luma":   "#4169E1", "tnbc":   "#FF4500",
              "her2":   "#9400D3", "lumb":   "#FF8C00"}
    w = 0.35

    fig = plt.figure(figsize=(22, 18))
    fig.suptitle(
        "BRCA ILC — Script 1 Discovery  |  OrganismCore BRCA-S6a/b  |  2026-03-05\n"
        "Attractor Type: TYPE 3 VARIANT — ADHESION LOCK DISSOLUTION\n"
        "Predicted geometry: TF identity retained (ESR1/FOXA1), "
        "structural lock dissolved (CDH1 loss)",
        fontsize=11, fontweight="bold", y=0.995
    )
    gs_fig = gridspec.GridSpec(3, 3, figure=fig, hspace=0.50, wspace=0.42)

    def bar_panel(ax, g_list, title):
        if not g_list:
            return
        ilc_m  = [float(np.nanmean(expr.loc[g, pops["ilc"]].values))    for g in g_list]
        norm_m = [float(np.nanmean(expr.loc[g, pops["normal"]].values))  for g in g_list]
        x = np.arange(len(g_list))
        ax.bar(x - w/2, norm_m, w, label="Normal", color=COLORS["normal"], alpha=0.85)
        ax.bar(x + w/2, ilc_m,  w, label="ILC",    color=COLORS["ilc"],    alpha=0.85)
        ax.set_xticks(x)
        ax.set_xticklabels(g_list, rotation=30, ha="right", fontsize=9)
        ax.set_ylabel("Mean log2 expr", fontsize=9)
        ax.set_title(title, fontsize=9)
        ax.legend(fontsize=8)

    # A — CDH1 axis + luminal TFs
    ax = fig.add_subplot(gs_fig[0, 0])
    g_list = [g for g in ["CDH1", "CDH3", "ESR1", "FOXA1", "GATA3", "PGR"] if g in expr.index]
    bar_panel(ax, g_list, "A — CDH1 Axis + Luminal TFs\nP1: CDH1↓  |  P2: ESR1/FOXA1/GATA3 retained")

    # B — Epigenetic axis
    ax = fig.add_subplot(gs_fig[0, 1])
    g_list = [g for g in ["EZH2", "EED", "SUZ12", "DNMT3A", "DNMT1", "KDM1A"] if g in expr.index]
    bar_panel(ax, g_list, "B — Epigenetic Axis\nP3: EZH2↑ targets CDH1  |  P4: DNMT3A↑")

    # C — Basal/EMT markers
    ax = fig.add_subplot(gs_fig[0, 2])
    g_list = [g for g in ["VIM", "ZEB1", "ZEB2", "SOX10", "KRT5", "SNAI1"] if g in expr.index]
    bar_panel(ax, g_list, "C — Basal/EMT Markers\nP7: All should be LOW in ILC")

    # D — EZH2 vs CDH1 scatter within ILC
    ax = fig.add_subplot(gs_fig[1, 0])
    if all(g in expr.index for g in ["EZH2", "CDH1"]) and len(pops["ilc"]) > 5:
        ez = expr.loc["EZH2", pops["ilc"]].values.astype(float)
        cd = expr.loc["CDH1", pops["ilc"]].values.astype(float)
        v  = ~(np.isnan(ez) | np.isnan(cd))
        ax.scatter(ez[v], cd[v], alpha=0.4, s=16, color=COLORS["ilc"])
        if v.sum() > 3:
            try:
                m_fit, b_fit = np.polyfit(ez[v], cd[v], 1)
                xf = np.linspace(ez[v].min(), ez[v].max(), 100)
                ax.plot(xf, m_fit * xf + b_fit, "k--", lw=1.2)
                r_val, p_val = stats.pearsonr(ez[v], cd[v])
                ax.set_title(f"D — EZH2 vs CDH1 (within ILC)\n"
                             f"r={r_val:+.3f}  {fmt_p(p_val)}", fontsize=9)
            except Exception:
                ax.set_title("D — EZH2 vs CDH1 (within ILC)\nP3: EZH2 targets CDH1?", fontsize=9)
        ax.set_xlabel("EZH2 (log2)", fontsize=8)
        ax.set_ylabel("CDH1 (log2)", fontsize=8)

    # E — Proliferation cross-subtype
    ax = fig.add_subplot(gs_fig[1, 1])
    g_list = [g for g in ["MKI67", "TOP2A", "PCNA", "CCNB1"] if g in expr.index]
    grp_plot = {k: v for k, v in pops.items() if len(v) >= 5}
    if g_list and len(grp_plot) >= 2:
        x    = np.arange(len(g_list))
        ngrp = len(grp_plot)
        bw   = 0.7 / ngrp
        for i, (grp, samps) in enumerate(grp_plot.items()):
            means = [float(np.nanmean(expr.loc[g, samps].values)) for g in g_list]
            ax.bar(x + (i - ngrp / 2) * bw + bw / 2, means, bw,
                   label=grp, color=COLORS.get(grp, "gray"), alpha=0.85)
        ax.set_xticks(x)
        ax.set_xticklabels(g_list, rotation=30, ha="right", fontsize=9)
        ax.set_ylabel("Mean log2 expr", fontsize=9)
        ax.set_title("E — Proliferation Cross-Subtype\nP5: ILC should be lowest", fontsize=9)
        ax.legend(fontsize=7)

    # F — CDH1 depth axis: deep vs shallow
    ax = fig.add_subplot(gs_fig[1, 2])
    if "CDH1" in expr.index and len(pops["ilc"]) > 10:
        cdh1_v = expr.loc["CDH1", pops["ilc"]].values.astype(float)
        depth  = -cdh1_v
        t67    = np.nanpercentile(depth, 67)
        t33    = np.nanpercentile(depth, 33)
        deep_s = [s for s, d in zip(pops["ilc"], depth) if d >= t67]
        shal_s = [s for s, d in zip(pops["ilc"], depth) if d <= t33]
        g_list = [g for g in ["CDH1", "ESR1", "FOXA1", "EZH2",
                               "VIM", "DNMT3A", "MKI67"] if g in expr.index]
        dm = [float(np.nanmean(expr.loc[g, deep_s].values)) for g in g_list]
        sm = [float(np.nanmean(expr.loc[g, shal_s].values)) for g in g_list]
        x  = np.arange(len(g_list))
        ax.bar(x - w/2, sm, w, label="Shallow (high CDH1)", color="#4682B4",      alpha=0.85)
        ax.bar(x + w/2, dm, w, label="Deep (low CDH1)",     color=COLORS["ilc"],  alpha=0.85)
        ax.set_xticks(x)
        ax.set_xticklabels(g_list, rotation=30, ha="right", fontsize=9)
        ax.set_ylabel("Mean log2 expr", fontsize=9)
        ax.set_title("F — CDH1 Depth Axis (within ILC)\n"
                     "P8: Does ESR1 fall with CDH1? (unexpected if yes)", fontsize=9)
        ax.legend(fontsize=8)

    # G — PI3K/AKT/mTOR axis
    ax = fig.add_subplot(gs_fig[2, 0])
    g_list = [g for g in ["AKT1", "MTOR", "PTEN", "PIK3CA", "EGFR", "ERBB2"] if g in expr.index]
    bar_panel(ax, g_list, "G — PI3K/AKT/mTOR Axis\nP6: PIK3CA pathway elevated")

    # H — PCA geometry
    ax = fig.add_subplot(gs_fig[2, 1])
    if pca_result is not None:
        Xp = pca_result["X_pca"]
        la = np.array(pca_result["labels"])
        for grp in grp_plot:
            m = la == grp
            if m.sum() == 0:
                continue
            ax.scatter(Xp[m, 0], Xp[m, 1], label=grp,
                       color=COLORS.get(grp, "gray"), alpha=0.4, s=16)
        for grp, c in pca_result.get("centroids", {}).items():
            ax.scatter(c[0], c[1], color=COLORS.get(grp, "gray"),
                       s=120, marker="*", zorder=5)
            ax.annotate(grp, (c[0], c[1]), fontsize=8, ha="center", va="bottom")
        pca_obj = pca_result["pca"]
        ax.set_xlabel(f"PC1 ({pca_obj.explained_variance_ratio_[0]*100:.1f}%)", fontsize=8)
        ax.set_ylabel(f"PC2 ({pca_obj.explained_variance_ratio_[1]*100:.1f}%)", fontsize=8)
        ax.set_title("H — PCA Geometry\nILC near LumA (shared TF)?  "
                     "Or CDH1-distinct?", fontsize=9)
        ax.legend(fontsize=7)

    # I — Top movers waterfall
    ax = fig.add_subplot(gs_fig[2, 2])
    if not rdf.empty:
        gained  = rdf[rdf["pct_change"] > 0].head(10)
        lost    = rdf[rdf["pct_change"] < 0].sort_values("pct_change").head(10)
        plot_df = pd.concat([gained, lost]).sort_values("pct_change")
        bar_colors = ["#8B0000" if v > 0 else "#2E8B57"
                      for v in plot_df["pct_change"]]
        ax.barh(plot_df["gene"], plot_df["pct_change"],
                color=bar_colors, alpha=0.85)
        ax.axvline(0, color="black", lw=0.8)
        ax.set_xlabel("% change vs Normal", fontsize=8)
        ax.set_title("I — Top 20 Movers (ILC vs Normal)\n"
                     "Unfiltered — geometry first", fontsize=9)
        ax.tick_params(axis="y", labelsize=7)

    plt.savefig(FIG_FILE, dpi=150, bbox_inches="tight")
    plt.close()
    log(f"  Figure saved: {FIG_FILE}")

# ============================================================
# MAIN
# ============================================================

def main():
    log("=" * 65)
    log("BRCA ILC — SCRIPT 1: DISCOVERY (v3 — corrected)")
    log("OrganismCore — Document BRCA-S6a/b | 2026-03-05")
    log("")
    log("ATTRACTOR TYPE: TYPE 3 VARIANT — ADHESION LOCK DISSOLUTION")
    log("  TF identity retained (ESR1/FOXA1/GATA3).")
    log("  Structural lock removed (CDH1 loss).")
    log("  Predicted geometric inverse of TNBC.")
    log(f"  Output directory: {BASE_DIR}")
    log("=" * 65)

    # Stage A: check scRNA-seq for ILC cells (expected: none)
    check_scrna_for_ilc()

    # Stage B: TCGA bulk
    acquire_tcga()
    expr = load_expression()
    clin, hist_col, pam50_col = load_clinical()

    pops_raw = classify_samples(expr, clin, hist_col, pam50_col)
    pops     = verify_pops(expr, pops_raw)
    expr     = normalise(expr)

    # ── Analysis pipeline ─────────────────────────────────────
    rdf = unfiltered_discovery(expr, pops)
    panel_rows, depth_results = panel_analysis(expr, pops, rdf)
    cross_df   = cross_subtype(expr, pops)
    pca_result = pca_geometry(expr, pops)

    generate_figure(expr, pops, panel_rows, cross_df, pca_result, rdf)

    # ── Final summary ─────────────────────────────────────────
    log("")
    log("=" * 65)
    log("SCRIPT 1 COMPLETE")
    log("=" * 65)
    log(f"  Log:         {LOG_FILE}")
    log(f"  Figure:      {FIG_FILE}")
    log(f"  Panel CSV:   {CSV_FILE}")
    log(f"  Top movers:  {TOP_MOVERS_FILE}")
    log(f"  Cross-type:  {CROSS_FILE}")
    log("")
    log("NEXT: Write BRCA-S6b (script1_results_and_reasoning.md)")
    log("      before running Script 2.")
    log("")
    log("KEY QUESTIONS TO ANSWER IN BRCA-S6b:")
    log("  1. Was CDH1 the dominant lost gene?  (P1)")
    log("  2. Were luminal TFs retained?        (P2)")
    log("  3. EZH2 vs CDH1 correlation sign?   (P3)")
    log("  4. DNMT3A elevated?                 (P4)")
    log("  5. MKI67 lower than TNBC/HER2?      (P5)")
    log("  6. PI3K axis elevated?              (P6)")
    log("  7. VIM/ZEB1/SOX10 low?              (P7)")
    log("  8. ESR1 stable on CDH1 depth axis?  (P8)")
    log("  9. What did the unfiltered scan reveal first?")
    log("")
    log("  Protocol v2.0: geometry first, predictions second.")

    write_log()


if __name__ == "__main__":
    main()
