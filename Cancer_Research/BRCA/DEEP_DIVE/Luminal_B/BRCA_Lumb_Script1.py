"""
BRCA LUMINAL B — SCRIPT 1
OrganismCore — Document BRCA-S5b | 2026-03-05

SELF-CONTAINED: Downloads GSE176078 to its OWN local
directory (lumb_s1_analysis/data/). Does NOT search for,
reference, or reuse data from any other subtype directory.
Running this script from a clean folder produces identical
results to running it anywhere else.

DATASET:
  GSE176078 — Wu et al. 2021, Nature Genetics
  PMID: 34493872
  100,064 cells — 26 primary breast tumors
  scRNA-seq 10X Chromium

  Pre-merged matrix TAR (NOT _RAW.tar):
  GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz
    → count_matrix_sparse.mtx   (genes × cells)
    → count_matrix_genes.tsv    (gene names)
    → count_matrix_barcodes.tsv (cell barcodes)
    → metadata.csv              (cell annotations)

GEOMETRY-FIRST (Protocol v2.0):
  Discovery script. Unfiltered top movers are printed
  BEFORE any prediction panel is applied.
  Read geometry output (Sections 1-3) before Section 4.

CRITICAL DATA NOTE:
  It is UNCERTAIN whether GSE176078 contains a
  "Cancer LumB SC" label in metadata.csv.
  Confirmed present: "Cancer LumA SC" n=7,742
                     "Cancer Basal SC" n=4,312
  LumB label: UNKNOWN — script detects at load time.
  Fallback: high-MKI67 ER+ Cancer proxy (lower rigor).
  Strategy is printed at Step 1 before any analysis runs.

PREDICTIONS LOCKED IN BRCA-S5a (2026-03-05):
  Attractor type: Type 3 (dominant) + partial Type 1
  P1.  FOXA1, GATA3 LOWER in LumB than LumA
       ESR1 EQUAL between LumB and LumA
       PGR LOWER in LumB than LumA
  P2.  MKI67, CCND1, CDK4, MYC, TOP2A ELEVATED vs LumA
       CDKN1A DOWN or equal vs LumA
  P3.  EZH2 INTERMEDIATE: above LumA (+18.5% flat),
       below TNBC. Predicted +30-60% vs Mature Luminal.
  P4.  Dual depth score:
         w1*(1 - norm(CDKN1A)) [proliferative axis]
       + w2*(1 - norm(mean(FOXA1,GATA3))) [identity axis]
       + w3*norm(MKI67)
       w1=0.4, w2=0.4, w3=0.2
  P5.  ER circuit more broken in LumB than LumA:
       PGR more suppressed, r(ESR1,PGR) lower in LumB
  P6.  Controls flat: KRT5, KRT14, SOX10 absent
  P9.  Cross-subtype gradient:
       ESR1/FOXA1/GATA3: LumA > LumB > TNBC
       EZH2:             LumA < LumB < TNBC

FIX vs v1:
  Removed 'global SC_DIR' from inside acquire_data().
  acquire_data() now returns sc_dir as a local variable.
  SC_DIR module-level constant is never mutated.
  main() receives sc_dir from acquire_data() return value
  and passes it to all downstream functions.
  Pattern matches BRCA_HER2_script1.py exactly.

Author:    Eric Robert Lawson
Framework: OrganismCore
Protocol:  Workflow_Protocol.md v2.0
Document:  BRCA-S5b
Date:      2026-03-05
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
from scipy.io import mmread
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# ============================================================
# CONFIGURATION — SELF-CONTAINED
# All paths are LOCAL to this script's directory.
# No path crawling. No cross-script cache reuse.
# SC_DIR is the EXPECTED extraction path — never mutated.
# acquire_data() returns the ACTUAL sc_dir as a local.
# ============================================================

GEO_ACCESSION  = "GSE176078"
SC_TARFILE     = "GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"
SC_INNER_DIR   = "Wu_etal_2021_BRCA_scRNASeq"
SC_URL         = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/"
    "GSE176078/suppl/"
    "GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"
)
BULK_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/"
    "GSE176078/suppl/"
    "GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts.txt.gz"
)
BULK_FILENAME = "GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts.txt.gz"

SC_FILES_INSIDE = [
    "count_matrix_sparse.mtx",
    "count_matrix_genes.tsv",
    "count_matrix_barcodes.tsv",
    "metadata.csv",
]

# All output is local to THIS script's folder
SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
BASE_DIR    = os.path.join(SCRIPT_DIR, "lumb_s1_analysis")
DATA_DIR    = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
LOG_FILE    = os.path.join(RESULTS_DIR, "lumb_s1_log.txt")
FIG_FILE    = os.path.join(RESULTS_DIR, "lumb_s1_figure.png")
CSV_FILE    = os.path.join(RESULTS_DIR, "lumb_s1_results.csv")
CACHE_FILE  = os.path.join(RESULTS_DIR, "expr_cache_lumb_s1.csv")

# Expected extraction path — used as default in acquire_data().
# acquire_data() returns the ACTUAL path (may differ if TAR
# extracts flat). Never assigned via 'global'.
SC_DIR_DEFAULT = os.path.join(DATA_DIR, SC_INNER_DIR)

# Local paths for downloaded files
SC_TARPATH  = os.path.join(DATA_DIR, SC_TARFILE)
BULK_PATH   = os.path.join(DATA_DIR, BULK_FILENAME)

for d in [BASE_DIR, DATA_DIR, RESULTS_DIR]:
    os.makedirs(d, exist_ok=True)

# ============================================================
# CELL POPULATION LABELS — Wu et al. 2021 celltype_subset
# ============================================================

CT_COL       = "celltype_subset"
CANCER_LUMB  = "Cancer LumB SC"
CANCER_LUMA  = "Cancer LumA SC"
CANCER_BASAL = "Cancer Basal SC"
MATURE_LUM   = "Mature Luminal"
LUMINAL_PROG = "Luminal Progenitors"
CANCER_CYCL  = "Cancer Cycling"

# ============================================================
# GENE PANELS — Locked in BRCA-S5a (2026-03-05)
# ============================================================

IDENTITY_TFS = ["FOXA1", "GATA3", "ESR1", "PGR", "SPDEF"]
PROLIF_GENES = ["MKI67", "CCND1", "CDK4", "MYC", "TOP2A",
                "CCNE1", "PCNA", "CDK2"]
EPIGENETIC   = ["EZH2", "HDAC1", "HDAC2", "KDM1A",
                "DNMT3A", "DNMT3B", "EED", "TET2"]
DEPTH_GENES  = ["CDKN1A", "CDKN2A", "RB1", "TP53", "TGFBR2"]
ER_CIRCUIT   = ["ESR1", "PGR", "NCOA1", "NCOA2", "TGFBR2",
                "CCND1", "ERBB2"]
DRUG_TARGETS = ["CDK4", "CDK6", "EZH2", "AR", "ERBB2",
                "PARP1", "CD274"]
GRADIENT_GENES = ["ESR1", "FOXA1", "GATA3", "EZH2",
                  "MKI67", "KRT5", "KRT14"]
CONTROLS = ["KRT5", "KRT14", "SOX10", "CDX2", "NKX2-1",
            "NKX3-1", "OLIG2", "SPI1"]

LUMA_REF = {
    "FOXA1":  {"pct_vs_mature": +37.3},
    "GATA3":  {"pct_vs_mature": +33.8},
    "ESR1":   {"pct_vs_mature": -30.1},
    "PGR":    {"pct_vs_mature": -54.8},
    "CDKN1A": {"pct_vs_mature": -74.3},
    "EZH2":   {"pct_vs_mature": +18.5},
    "MYC":    {"pct_vs_mature": -44.5},
}

ALL_GENES = list(dict.fromkeys(
    IDENTITY_TFS + PROLIF_GENES + EPIGENETIC +
    DEPTH_GENES + ER_CIRCUIT + DRUG_TARGETS +
    GRADIENT_GENES + CONTROLS
))

# ============================================================
# LOGGING
# ============================================================

log_lines = []

def log(msg=""):
    print(msg)
    log_lines.append(str(msg))

def write_log():
    with open(LOG_FILE, "w") as f:
        f.write("\n".join(log_lines))

# ============================================================
# DATA ACQUISITION — SELF-CONTAINED
# Returns: sc_dir (str) — actual path to extracted files
#          bulk_path (str or None)
# No global mutations. sc_dir is a local variable throughout.
# ============================================================

def _sc_dir_is_complete(path):
    """True if all four sc files are present and non-empty."""
    for fname in SC_FILES_INSIDE:
        fp = os.path.join(path, fname)
        if not os.path.exists(fp) or os.path.getsize(fp) < 1000:
            return False
    return True


def fetch_url(url, dest, retries=3):
    """Download url to dest with retry. Returns dest or None."""
    for attempt in range(1, retries + 1):
        try:
            log(f"  Attempt {attempt}: {url[-70:]}")
            req = urllib.request.Request(
                url, headers={"User-Agent": "Mozilla/5.0"}
            )
            with urllib.request.urlopen(req, timeout=300) as r:
                data = r.read()
            with open(dest, "wb") as f:
                f.write(data)
            log(f"  Saved: {dest} ({len(data)/1e6:.1f} MB)")
            return dest
        except Exception as e:
            log(f"  Attempt {attempt} failed: {e}")
            if attempt < retries:
                time.sleep(5)
    return None


def acquire_data():
    """
    Download and extract GSE176078 into DATA_DIR.
    Returns (sc_dir, bulk_path) where sc_dir is the ACTUAL
    path to the extracted matrix files — either
    DATA_DIR/Wu_etal_2021_BRCA_scRNASeq/ (normal TAR layout)
    or DATA_DIR/ (flat extraction). Never mutates SC_DIR_DEFAULT.
    """
    log("=" * 65)
    log("STEP 0: DATA ACQUISITION — SELF-CONTAINED")
    log(f"  All files downloaded to: {DATA_DIR}")
    log(f"  Dataset: {GEO_ACCESSION} — Wu et al. 2021")
    log("=" * 65)

    # ── scRNA-seq ──────────────────────────────────────────
    log("\n--- scRNA-seq (GSE176078 merged matrix TAR) ---")

    # Start with the expected subdirectory path
    sc_dir = SC_DIR_DEFAULT

    if _sc_dir_is_complete(sc_dir):
        log(f"  Already extracted: {sc_dir}")
        for fname in SC_FILES_INSIDE:
            fp = os.path.join(sc_dir, fname)
            log(f"    OK: {fname} ({os.path.getsize(fp)/1e6:.1f} MB)")
    else:
        # Check flat layout (files directly in DATA_DIR)
        if _sc_dir_is_complete(DATA_DIR):
            log(f"  Files found flat in DATA_DIR: {DATA_DIR}")
            sc_dir = DATA_DIR
        else:
            # Need to download and/or extract
            if (not os.path.exists(SC_TARPATH) or
                    os.path.getsize(SC_TARPATH) < 1_000_000):
                log("  TAR not found locally. Downloading (~2 GB)...")
                result = fetch_url(SC_URL, SC_TARPATH)
                if result is None:
                    log("  FATAL: TAR download failed.")
                    return None, None
            else:
                log(f"  TAR already present: {SC_TARPATH} "
                    f"({os.path.getsize(SC_TARPATH)/1e6:.0f} MB)")

            log(f"  Extracting to {DATA_DIR} ...")
            try:
                with tarfile.open(SC_TARPATH, "r:gz") as tf:
                    members = tf.getnames()
                    log(f"  TAR contains {len(members)} items:")
                    for m in members[:8]:
                        log(f"    {m}")
                    tf.extractall(DATA_DIR)
                log("  Extraction complete.")
            except Exception as e:
                log(f"  FATAL: TAR extraction failed: {e}")
                return None, None

            # Determine actual extraction layout
            if _sc_dir_is_complete(SC_DIR_DEFAULT):
                sc_dir = SC_DIR_DEFAULT
                log(f"  Extracted into subdirectory: {sc_dir}")
            elif _sc_dir_is_complete(DATA_DIR):
                sc_dir = DATA_DIR
                log(f"  TAR extracted flat into DATA_DIR: {sc_dir}")
            else:
                log("  FATAL: Expected files not found after extraction.")
                log(f"  DATA_DIR contents:")
                for f in os.listdir(DATA_DIR)[:20]:
                    log(f"    {f}")
                return None, None

        for fname in SC_FILES_INSIDE:
            fp = os.path.join(sc_dir, fname)
            log(f"    OK: {fname} ({os.path.getsize(fp)/1e6:.1f} MB)")

    # ── Bulk RNA-seq (optional) ────────────────────────────
    log("\n--- Bulk RNA-seq (GSE176078, optional) ---")
    bulk_path = None
    if os.path.exists(BULK_PATH) and os.path.getsize(BULK_PATH) > 10_000:
        log(f"  Already present: {BULK_PATH}")
        bulk_path = BULK_PATH
    else:
        log("  Downloading bulk RNA-seq...")
        result = fetch_url(BULK_URL, BULK_PATH)
        if result is None:
            log("  Bulk download failed — bulk step will be skipped.")
        else:
            bulk_path = result

    return sc_dir, bulk_path

# ============================================================
# STEP 1: LOAD METADATA + LumB LABEL DETECTION
# Takes sc_dir as argument — no globals read.
# ============================================================

def load_metadata(sc_dir):
    log("=" * 65)
    log("STEP 1: LOAD METADATA + LumB LABEL DETECTION")
    log("=" * 65)

    meta_path = os.path.join(sc_dir, "metadata.csv")
    meta = pd.read_csv(meta_path, index_col=0)
    log(f"  Cells: {len(meta)}")
    log(f"  Columns: {list(meta.columns[:12])}")

    if CT_COL not in meta.columns:
        log(f"  WARNING: '{CT_COL}' not found. Searching alternatives...")
        for col in meta.columns:
            if "celltype" in col.lower() or "subset" in col.lower():
                log(f"  Using '{col}' as cell type column.")
                meta[CT_COL] = meta[col]
                break
    if CT_COL not in meta.columns:
        log("  FATAL: Cannot find cell type column.")
        return None, "FATAL"

    log(f"\n  All cell types ({CT_COL}):")
    ct_counts = meta[CT_COL].value_counts()
    for ct, n in ct_counts.items():
        log(f"    {ct:<42} n={n}")

    log("\n  LABEL DETECTION:")
    has_lumb  = CANCER_LUMB  in ct_counts.index
    n_lumb    = ct_counts.get(CANCER_LUMB,  0)
    n_luma    = ct_counts.get(CANCER_LUMA,  0)
    n_basal   = ct_counts.get(CANCER_BASAL, 0)
    n_mature  = ct_counts.get(MATURE_LUM,   0)
    n_prog    = ct_counts.get(LUMINAL_PROG, 0)

    log(f"    'Cancer LumB SC' present: {has_lumb} (n={n_lumb})")
    log(f"    'Cancer LumA SC' present: "
        f"{CANCER_LUMA in ct_counts.index} (n={n_luma})")
    log(f"    'Cancer Basal SC' present: "
        f"{CANCER_BASAL in ct_counts.index} (n={n_basal})")
    log(f"    'Mature Luminal' present: "
        f"{MATURE_LUM in ct_counts.index} (n={n_mature})")
    log(f"    'Luminal Progenitors' present: "
        f"{LUMINAL_PROG in ct_counts.index} (n={n_prog})")

    if has_lumb and n_lumb >= 100:
        strategy = "DIRECT_LABEL"
        log("\n  STRATEGY: DIRECT_LABEL — using 'Cancer LumB SC' directly.")
    elif has_lumb and n_lumb > 0:
        strategy = "DIRECT_LABEL_LOW_N"
        log(f"\n  STRATEGY: DIRECT_LABEL_LOW_N — "
            f"n={n_lumb} cells (low power).")
    else:
        strategy = "FALLBACK_PROXY"
        log("\n  STRATEGY: FALLBACK_PROXY")
        log("  'Cancer LumB SC' NOT in metadata.")
        log("  Constructing LumB proxy: Cancer cells that are")
        log("  not LumA / Basal / Cycling, high MKI67 (>p50).")
        log("  NOTE: lower methodological rigor — flag in artifact.")

    meta["is_luma"]    = meta[CT_COL] == CANCER_LUMA
    meta["is_basal"]   = meta[CT_COL] == CANCER_BASAL
    meta["is_mature"]  = meta[CT_COL] == MATURE_LUM
    meta["is_prog"]    = meta[CT_COL] == LUMINAL_PROG
    meta["is_cycling"] = meta[CT_COL] == CANCER_CYCL

    if has_lumb:
        meta["is_lumb"]         = meta[CT_COL] == CANCER_LUMB
        meta["is_cancer_other"] = False
    else:
        meta["is_lumb"] = False
        meta["is_cancer_other"] = (
            meta[CT_COL].str.contains("Cancer", case=False, na=False) &
            ~meta["is_luma"] &
            ~meta["is_basal"] &
            ~meta["is_cycling"]
        )
        log(f"  'Cancer other' pool (for proxy): "
            f"{meta['is_cancer_other'].sum()} cells")

    log(f"\n  FINAL POPULATION SIZES:")
    log(f"    LumB (strategy={strategy}): {meta['is_lumb'].sum()}")
    log(f"    LumA SC:                    {n_luma}")
    log(f"    Basal SC (TNBC ref):        {n_basal}")
    log(f"    Mature Luminal (normal):    {n_mature}")
    log(f"    Luminal Progenitors:        {n_prog}")

    if n_mature < 50:
        log(f"  WARNING: Only {n_mature} Mature Luminal cells.")

    return meta, strategy

# ============================================================
# STEP 2: LOAD EXPRESSION
# Takes sc_dir as argument — no globals read.
# ============================================================

def load_sc_expression(sc_dir):
    log("=" * 65)
    log("STEP 2: LOAD EXPRESSION (MTX -> log1p)")
    log("=" * 65)

    if (os.path.exists(CACHE_FILE) and
            os.path.getsize(CACHE_FILE) > 100_000):
        sz = os.path.getsize(CACHE_FILE) / 1e6
        log(f"  Loading local cache: {CACHE_FILE} ({sz:.1f} MB)")
        expr = pd.read_csv(CACHE_FILE, index_col=0)
        log(f"  Cached shape: {expr.shape}")
        return expr

    mtx_f  = os.path.join(sc_dir, "count_matrix_sparse.mtx")
    gene_f = os.path.join(sc_dir, "count_matrix_genes.tsv")
    bc_f   = os.path.join(sc_dir, "count_matrix_barcodes.tsv")

    genes_df   = pd.read_csv(gene_f, sep="\t", header=None)
    gene_names = genes_df.iloc[:, 0].tolist()
    log(f"  Total genes in matrix: {len(gene_names)}")
    log(f"  First 5: {gene_names[:5]}")

    barcodes = pd.read_csv(bc_f, sep="\t", header=None
                           ).iloc[:, 0].tolist()
    log(f"  Total barcodes: {len(barcodes)}")

    log("  Reading MTX (may take 30-90 sec)...")
    mat = mmread(mtx_f).tocsr()
    log(f"  Matrix shape: {mat.shape}  (genes x cells)")

    gene_upper   = [g.upper() for g in gene_names]
    target_idx   = []
    target_names = []
    missing      = []

    for gene in ALL_GENES:
        gu = gene.upper()
        if gu in gene_upper:
            target_idx.append(gene_upper.index(gu))
            target_names.append(gene)
        else:
            missing.append(gene)

    log(f"  Target genes found:  {len(target_names)} / {len(ALL_GENES)}")
    if missing:
        log(f"  Not found: {missing}")

    log("  Extracting target gene rows...")
    expr_data = {}
    for i, (row_idx, gene_name) in enumerate(
        zip(target_idx, target_names)
    ):
        expr_data[gene_name] = mat[row_idx, :].toarray().flatten()
        if (i + 1) % 10 == 0:
            log(f"    {i+1}/{len(target_names)}")

    expr = pd.DataFrame(expr_data, index=barcodes)
    log(f"  Raw expression shape: {expr.shape}")

    lib = expr.sum(axis=1)
    log(f"\n  Library sizes (target genes, raw):")
    log(f"    Median={lib.median():.1f}  Mean={lib.mean():.1f}  "
        f"Min={lib.min():.1f}  Max={lib.max():.1f}")

    expr = np.log1p(expr)
    log("  log1p normalisation applied.")

    log("\n  Spot-check means (log1p):")
    for gene in (IDENTITY_TFS + PROLIF_GENES[:3]):
        if gene in expr.columns:
            log(f"    {gene:<14} mean={expr[gene].mean():.4f}  "
                f"nonzero={(expr[gene] > 0).sum()}")

    expr.to_csv(CACHE_FILE)
    log(f"\n  Cache saved: {CACHE_FILE}")
    return expr

# ============================================================
# STEP 3: CLASSIFY POPULATIONS
# ============================================================

def classify_populations(expr, meta, strategy):
    log("=" * 65)
    log("STEP 3: CLASSIFY POPULATIONS")
    log(f"  Strategy: {strategy}")
    log("=" * 65)

    common = expr.index.intersection(meta.index)
    log(f"  Expr/meta overlap: {len(common)} cells")
    expr_c = expr.loc[common]
    meta_c = meta.loc[common]

    pops = {}
    for name, flag in [
        ("mature", "is_mature"),
        ("prog",   "is_prog"),
        ("luma",   "is_luma"),
        ("basal",  "is_basal"),
    ]:
        idx = meta_c.index[
            meta_c.get(flag, pd.Series(False, index=meta_c.index))
        ]
        pops[name] = expr_c.loc[idx]
        log(f"  {name:<12}: {len(pops[name])} cells")

    if strategy in ("DIRECT_LABEL", "DIRECT_LABEL_LOW_N"):
        idx_lb = meta_c.index[meta_c["is_lumb"]]
        pops["lumb"] = expr_c.loc[idx_lb]
        log(f"  lumb (direct):  {len(pops['lumb'])} cells")
        if len(pops["lumb"]) < 50:
            log("  WARNING: n<50 LumB cells — results unreliable.")

    elif strategy == "FALLBACK_PROXY":
        log("\n  FALLBACK: Building LumB proxy from expression.")
        co_mask = meta_c.get(
            "is_cancer_other",
            pd.Series(False, index=meta_c.index)
        )
        co_idx  = meta_c.index[co_mask]
        co_expr = expr_c.loc[co_idx]
        log(f"  'Cancer other' pool: {len(co_expr)} cells")

        if len(co_expr) < 20:
            log("  FATAL: Insufficient cancer_other cells.")
            pops["lumb"] = pd.DataFrame()
        elif "MKI67" in co_expr.columns:
            thresh = co_expr["MKI67"].quantile(0.50)
            mask   = co_expr["MKI67"] > thresh
            pops["lumb"] = co_expr[mask]
            log(f"  MKI67 p50 threshold: {thresh:.4f}")
            log(f"  High-MKI67 LumB proxy: {len(pops['lumb'])} cells")
        else:
            log("  MKI67 not available — using all cancer_other.")
            pops["lumb"] = co_expr

    log(f"\n  FINAL POPULATION SIZES:")
    for name, df in pops.items():
        log(f"    {name:<12}: {len(df)} cells")

    if len(pops.get("lumb", pd.DataFrame())) < 20:
        log("  ERROR: Insufficient LumB cells.")
        return None

    return pops

# ============================================================
# SECTION 1: TOP MOVERS — UNFILTERED
# ============================================================

def top_movers_unfiltered(pops):
    log("")
    log("=" * 65)
    log("SECTION 1: TOP MOVERS — UNFILTERED")
    log("Read ALL THREE comparisons before Section 4.")
    log("=" * 65)

    lumb   = pops.get("lumb",   pd.DataFrame())
    mature = pops.get("mature", pd.DataFrame())
    luma   = pops.get("luma",   pd.DataFrame())
    prog   = pops.get("prog",   pd.DataFrame())

    def fmt_p(p):
        if p < 1e-100: return "p<1e-100 ***"
        if p < 0.001:  return f"p={p:.2e} ***"
        if p < 0.01:   return f"p={p:.2e}  **"
        if p < 0.05:   return f"p={p:.4f}   *"
        return              f"p={p:.4f}  ns"

    def compare(cancer_df, ref_df, cancer_lbl, ref_lbl):
        if len(cancer_df) < 10 or len(ref_df) < 3:
            log(f"  Skipped: {cancer_lbl}={len(cancer_df)} "
                f"{ref_lbl}={len(ref_df)}")
            return pd.DataFrame()

        rows = []
        for gene in cancer_df.columns:
            if gene not in ref_df.columns:
                continue
            tv = cancer_df[gene].values
            nv = ref_df[gene].values
            nm, tm = nv.mean(), tv.mean()
            if abs(nm) < 1e-6 and abs(tm) < 1e-6:
                continue
            chg = ((tm - nm) / abs(nm) * 100
                   if abs(nm) > 1e-6 else np.nan)
            try:
                _, p = stats.mannwhitneyu(
                    tv, nv, alternative="two-sided"
                )
            except Exception:
                p = 1.0
            rows.append({
                "gene":        gene,
                "ref_mean":    nm,
                "cancer_mean": tm,
                "change_pct":  chg,
                "abs_change":  abs(chg) if not np.isnan(chg) else 0,
                "p_value":     p,
                "direction":   "UP" if (chg or 0) > 0 else "DOWN",
            })

        rdf = pd.DataFrame(rows).sort_values(
            "abs_change", ascending=False
        )

        log(f"\n  {cancer_lbl} vs {ref_lbl}  "
            f"(n={len(cancer_df)} vs n={len(ref_df)})")

        for direction, label in [("DOWN", "LOST"), ("UP", "GAINED")]:
            subset = rdf[rdf["direction"] == direction].head(15)
            log(f"\n  TOP 15 {label}:")
            log(f"  {'Gene':<14} {'Ref':>9} {'Cancer':>9} "
                f"{'Change%':>9}  {'p-value':>16}")
            log(f"  {'-'*62}")
            for _, row in subset.iterrows():
                log(f"  {row['gene']:<14} {row['ref_mean']:>9.4f} "
                    f"{row['cancer_mean']:>9.4f} "
                    f"{row['change_pct']:>+8.1f}%  "
                    f"{fmt_p(row['p_value']):>16}")

        return rdf

    log("\n--- 1a: LumB vs Mature Luminal ---")
    rdf_mature = compare(lumb, mature, "LumB", "Mature_Luminal")

    log("\n--- 1b: LumB vs LumA (CRITICAL TEST) ---")
    log("  FOXA1 direction in this comparison is the")
    log("  single most important result of this script.")
    rdf_luma = compare(lumb, luma, "LumB", "LumA")

    log("\n--- 1c: LumB vs Luminal Progenitor ---")
    rdf_prog = compare(lumb, prog, "LumB", "Luminal_Progenitor")

    if not rdf_mature.empty:
        rdf_mature.to_csv(CSV_FILE, index=False)
        log(f"\n  Primary results (vs Mature) saved: {CSV_FILE}")

    return rdf_mature, rdf_luma, rdf_prog

# ============================================================
# SECTION 2: PCA GEOMETRY
# ============================================================

def pca_geometry(pops):
    log("")
    log("=" * 65)
    log("SECTION 2: PCA GEOMETRY")
    log("Landscape position — read before Section 4.")
    log("=" * 65)

    frames, labels = [], []
    for name, lbl in [
        ("lumb",   "LumB"),
        ("luma",   "LumA"),
        ("mature", "Mature"),
        ("prog",   "Progenitor"),
        ("basal",  "TNBC"),
    ]:
        df = pops.get(name, pd.DataFrame())
        if len(df) >= 5:
            frames.append(df)
            labels.extend([lbl] * len(df))

    if not frames:
        log("  No populations available.")
        return None, None

    combined = pd.concat(frames).fillna(0)
    labels   = np.array(labels[:len(combined)])

    lib  = combined.sum(axis=1)
    med  = lib.median()
    mask = (lib < med / 5) | (lib > med * 5)
    if mask.sum() > 0:
        log(f"  Removed {mask.sum()} library-size outliers.")
        combined = combined[~mask]
        labels   = labels[~mask]

    scaled = StandardScaler().fit_transform(combined)
    n_comp = min(5, combined.shape[1], combined.shape[0] - 1)
    pca    = PCA(n_components=n_comp)
    pcs    = pca.fit_transform(scaled)

    log("\n  Variance explained:")
    for i, v in enumerate(pca.explained_variance_ratio_):
        log(f"    PC{i+1}: {v*100:.2f}%")

    loadings = pd.Series(
        pca.components_[0], index=combined.columns
    ).sort_values(ascending=False)

    log("\n  PC1 top positive loadings:")
    for g, v in loadings.head(10).items():
        log(f"    {g:<18} {v:>+8.4f}")
    log("\n  PC1 top negative loadings:")
    for g, v in loadings.tail(10).items():
        log(f"    {g:<18} {v:>+8.4f}")

    log("\n  PC1 mean per population:")
    for lbl in ["Mature", "Progenitor", "LumA", "LumB", "TNBC"]:
        m = labels == lbl
        if m.sum() > 0:
            log(f"    {lbl:<15} mean={pcs[m, 0].mean():>+8.4f}  "
                f"n={m.sum()}")

    pc1_lb = pcs[labels == "LumB", 0]
    pc1_la = pcs[labels == "LumA", 0]
    if len(pc1_lb) > 5 and len(pc1_la) > 5:
        _, p_sep = stats.mannwhitneyu(pc1_lb, pc1_la,
                                      alternative="two-sided")
        q = ("SEPARATED" if p_sep < 0.001 else
             "PARTIAL"   if p_sep < 0.05  else "OVERLAP")
        log(f"\n  LumB vs LumA PC1 separation: {q}  p={p_sep:.2e}")
        log(f"    LumB mean={pc1_lb.mean():.4f}  "
            f"LumA mean={pc1_la.mean():.4f}")

    return pcs, labels

# ============================================================
# SECTION 3: DEPTH SCORE
# ============================================================

def compute_depth_score(pops):
    log("")
    log("=" * 65)
    log("SECTION 3: DEPTH SCORE (BRCA-S5a P4)")
    log("  w1=0.4*(1-norm(CDKN1A))")
    log("  w2=0.4*(1-norm(mean(FOXA1,GATA3)))")
    log("  w3=0.2*norm(MKI67)")
    log("=" * 65)

    lumb = pops.get("lumb", pd.DataFrame())
    luma = pops.get("luma", pd.DataFrame())
    if len(lumb) < 10:
        log("  Insufficient LumB cells.")
        return None, []

    def norm01(s):
        mn, mx = s.min(), s.max()
        if mx - mn < 1e-8:
            return pd.Series(0.5, index=s.index)
        return (s - mn) / (mx - mn)

    parts, total_w = [], 0.0
    w1, w2, w3 = 0.4, 0.4, 0.2

    if "CDKN1A" in lumb.columns:
        c1 = 1 - norm01(lumb["CDKN1A"])
        parts.append((w1, c1, "CDKN1A_axis"))
        total_w += w1
        log(f"  C1 (1-norm(CDKN1A)): mean={c1.mean():.4f}")
    elif "CDK4" in lumb.columns:
        c1 = norm01(lumb["CDK4"])
        parts.append((w1, c1, "CDK4_fallback"))
        total_w += w1
        log(f"  C1 fallback (norm(CDK4)): mean={c1.mean():.4f}")

    id_avail = [g for g in ["FOXA1", "GATA3"] if g in lumb.columns]
    if id_avail:
        c2 = 1 - norm01(lumb[id_avail].mean(axis=1))
        parts.append((w2, c2, f"identity({id_avail})"))
        total_w += w2
        log(f"  C2 (1-norm({id_avail})): mean={c2.mean():.4f}")

    if "MKI67" in lumb.columns:
        c3 = norm01(lumb["MKI67"])
        parts.append((w3, c3, "MKI67"))
        total_w += w3
        log(f"  C3 (norm(MKI67)): mean={c3.mean():.4f}")

    if not parts:
        log("  FATAL: No components available.")
        return None, []

    depth = sum(w * c for w, c, _ in parts) / total_w
    formula = " + ".join(f"w={w}*{lbl}" for w, _, lbl in parts)
    log(f"\n  Formula: {formula}")
    log(f"  Depth score ({len(lumb)} LumB cells):")
    log(f"    Mean={depth.mean():.4f}  Median={depth.median():.4f}  "
        f"Std={depth.std():.4f}")
    log(f"    Min={depth.min():.4f}  Max={depth.max():.4f}  "
        f"Q25={depth.quantile(0.25):.4f}  Q75={depth.quantile(0.75):.4f}")

    if len(luma) > 10 and "CDKN1A" in luma.columns:
        log("\n  LumB vs LumA depth (joint normalisation, CDKN1A axis):")
        joint   = pd.concat([lumb, luma])
        lumb_d  = 1 - norm01(joint["CDKN1A"].loc[lumb.index])
        luma_d  = 1 - norm01(joint["CDKN1A"].loc[luma.index])
        log(f"    LumB mean={lumb_d.mean():.4f}  "
            f"std={lumb_d.std():.4f}")
        log(f"    LumA mean={luma_d.mean():.4f}  "
            f"std={luma_d.std():.4f}")
        try:
            _, p = stats.mannwhitneyu(lumb_d, luma_d,
                                      alternative="two-sided")
            log(f"    p(LumB vs LumA depth): {p:.4e}")
        except Exception as e:
            log(f"    Comparison error: {e}")

    log(f"\n  Top depth correlations (LumB, |r| ranked):")
    log(f"  {'Gene':<18} {'r':>8}  p-value")
    log(f"  {'-'*36}")
    corrs = []
    for gene in lumb.columns:
        try:
            r, p = stats.pearsonr(depth.values, lumb[gene].values)
            corrs.append((gene, r, p))
        except Exception:
            pass
    corrs.sort(key=lambda x: abs(x[1]), reverse=True)
    for g, r, p in corrs[:25]:
        sig = ("***" if p < 0.001 else
               "**"  if p < 0.01  else
               "*"   if p < 0.05  else "ns")
        log(f"  {g:<18} {r:>+8.4f}  p={p:.2e} {sig}")

    return depth, corrs

# ============================================================
# SECTION 4: PREDICTION PANEL CHECK
# ============================================================

def prediction_check(pops, depth, corrs, rdf_luma):
    log("")
    log("=" * 65)
    log("SECTION 4: PREDICTION PANEL CHECK")
    log("Predictions locked 2026-03-05 in BRCA-S5a.")
    log("=" * 65)

    lumb   = pops.get("lumb",   pd.DataFrame())
    mature = pops.get("mature", pd.DataFrame())
    luma   = pops.get("luma",   pd.DataFrame())
    basal  = pops.get("basal",  pd.DataFrame())

    def fmt_p(p):
        if p is None:   return "  ns (pre-stated)"
        if p < 1e-100:  return "p<1e-100 ***"
        if p < 0.001:   return f"p={p:.2e} ***"
        if p < 0.01:    return f"p={p:.2e}  **"
        if p < 0.05:    return f"p={p:.4f}   *"
        return          f"p={p:.4f}  ns"

    def test(gene, ref, cancer):
        if (gene not in ref.columns or
                gene not in cancer.columns or
                len(ref) < 3 or len(cancer) < 3):
            return None, None, None, None
        nv, tv = ref[gene].values, cancer[gene].values
        nm, tm = nv.mean(), tv.mean()
        chg = (tm - nm) / abs(nm) * 100 if abs(nm) > 1e-6 else np.nan
        try:
            _, p = stats.mannwhitneyu(tv, nv, alternative="two-sided")
        except Exception:
            p = 1.0
        return nm, tm, chg, p

    def dr(gene):
        if depth is None or gene not in lumb.columns:
            return np.nan, np.nan
        try:
            return stats.pearsonr(depth.values, lumb[gene].values)
        except Exception:
            return np.nan, np.nan

    log(f"\n  P1 — IDENTITY TFs vs MATURE (predicted: RETAINED):")
    log(f"  {'Gene':<12} {'Mature':>8} {'LumB':>8} {'%chg':>8}  "
        f"{'p':>14}  LumA_ref  Verdict")
    log(f"  {'-'*75}")
    for gene in IDENTITY_TFS:
        nm, tm, chg, p = test(gene, mature, lumb)
        if nm is None:
            log(f"  {gene:<12} NOT IN DATA")
            continue
        rv, _ = dr(gene)
        luma_ref = LUMA_REF.get(gene, {}).get("pct_vs_mature", "?")
        actual = ("SUPPRESSED" if (chg or 0) < -20 and p < 0.05 else
                  "ELEVATED"   if (chg or 0) >  20 and p < 0.05 else
                  "FLAT")
        sym = "✓" if actual != "SUPPRESSED" else "✗"
        log(f"  {sym} {gene:<10} {nm:>8.4f} {tm:>8.4f} "
            f"{chg:>+7.1f}%  {fmt_p(p):>14}  "
            f"ref={str(luma_ref):>6}%  {actual}")

    log(f"\n  P1b — IDENTITY TFs: LumB vs LumA (CRITICAL TEST):")
    log(f"  Predicted: FOXA1 < LumA, GATA3 < LumA, ESR1 ~= LumA")
    log(f"  {'Gene':<12} {'LumA':>8} {'LumB':>8} {'%diff':>8}  "
        f"{'p':>14}  Predicted  -> Result")
    log(f"  {'-'*72}")
    foxa1_result = None
    for gene in IDENTITY_TFS:
        nm, tm, chg, p = test(gene, luma, lumb)
        if nm is None:
            log(f"  {gene:<12} NOT IN DATA")
            continue
        predicted = ("EQUAL"  if gene == "ESR1" else "LOWER")
        actual    = ("LOWER"  if (chg or 0) < -10 and p < 0.05 else
                     "HIGHER" if (chg or 0) >  10 and p < 0.05 else
                     "EQUAL")
        sym = "✓" if actual == predicted else "✗"
        crit = "  <- CRITICAL" if gene == "FOXA1" else ""
        log(f"  {sym} {gene:<10} {nm:>8.4f} {tm:>8.4f} "
            f"{chg:>+7.1f}%  {fmt_p(p):>14}  "
            f"pred={predicted:<7} -> {actual}{crit}")
        if gene == "FOXA1":
            foxa1_result = actual

    log(f"\n  FOXA1 verdict (P1b): {foxa1_result}")
    if foxa1_result == "LOWER":
        log("  CONFIRMED: LumB is less differentiated than LumA.")
        log("    Intermediate position model HOLDS.")
    elif foxa1_result == "EQUAL":
        log("  NOT CONFIRMED: LumB and LumA share FOXA1 level.")
        log("    Difference is proliferative, not differentiation-state.")
    else:
        log(f"  UNEXPECTED: FOXA1 {foxa1_result} in LumB vs LumA.")

    log(f"\n  P2 — PROLIFERATION vs LumA:")
    log(f"  {'Gene':<12} {'LumA':>8} {'LumB':>8} {'%diff':>8}  "
        f"{'p':>14}  r_depth  Verdict")
    log(f"  {'-'*68}")
    for gene in ["MKI67", "CCND1", "CDK4", "MYC", "TOP2A", "CDKN1A"]:
        nm, tm, chg, p = test(gene, luma, lumb)
        if nm is None:
            log(f"  {gene:<12} NOT IN DATA")
            continue
        rv, _ = dr(gene)
        if gene == "CDKN1A":
            verdict = ("DOWN"  if (chg or 0) < -5 else
                       "EQUAL" if abs(chg or 0) <= 15 else "UP-UNEXPECTED")
        else:
            verdict = ("UP"   if (chg or 0) > 10 and p < 0.05 else
                       "FLAT" if abs(chg or 0) <= 10 else "DOWN-UNEXPECTED")
        log(f"  {gene:<12} {nm:>8.4f} {tm:>8.4f} "
            f"{chg:>+7.1f}%  {fmt_p(p):>14}  "
            f"{rv:>+6.3f}  {verdict}")

    log(f"\n  P3 — EZH2 GRADIENT (CRITICAL TEST):")
    log(f"  Predicted: Mature < LumA (+18.5%) < LumB < TNBC")
    _, _, chg_la, p_la = test("EZH2", mature, luma)
    _, _, chg_lb, p_lb = test("EZH2", mature, lumb)
    _, _, chg_tn, _    = test("EZH2", mature, basal)
    rv_e, _ = dr("EZH2")
    log(f"    LumA vs Mature: {chg_la:>+8.1f}%  (ref: +18.5% ns)")
    log(f"    LumB vs Mature: {chg_lb:>+8.1f}%  {fmt_p(p_lb)}")
    log(f"    TNBC vs Mature: {chg_tn:>+8.1f}%")
    log(f"    r(EZH2, depth) in LumB: {rv_e:>+.4f}")
    if (chg_lb or 0) > (chg_la or 0) and (chg_lb or 0) < (chg_tn or 0) and (p_lb or 1) < 0.05:
        log("  P3 CONFIRMED — EZH2 intermediate. EZH2i flag active.")
    elif abs((chg_lb or 0) - (chg_la or 0)) < 15:
        log("  P3 FLAT — EZH2 equal to LumA. No EZH2i differentiation.")
    else:
        log("  P3 PARTIAL/UNEXPECTED — document in BRCA-S5b.")

    log(f"\n  P5 — ER CIRCUIT (LumA vs LumB vs Mature):")
    log(f"  Predicted: PGR more suppressed in LumB than LumA")
    log(f"  {'Gene':<12} {'LumA%':>10} {'LumB%':>10}  Result")
    log(f"  {'-'*46}")
    for gene in ["ESR1", "PGR", "NCOA1", "NCOA2", "TGFBR2"]:
        _, _, chg_la, _ = test(gene, mature, luma)
        _, _, chg_lb, _ = test(gene, mature, lumb)
        if chg_la is None or chg_lb is None:
            continue
        more_broken = ((chg_lb or 0) < (chg_la or 0) and gene != "ESR1")
        sym = "+" if more_broken else "="
        log(f"  {sym} {gene:<10} "
            f"LumA={chg_la:>+8.1f}%  LumB={chg_lb:>+8.1f}%  "
            f"{'MORE BROKEN' if more_broken else 'EQUAL/HIGHER'}")
    for lbl, df in [("LumB", lumb), ("LumA", luma)]:
        if ("ESR1" in df.columns and "PGR" in df.columns
                and len(df) > 20):
            try:
                r, p = stats.pearsonr(df["ESR1"].values, df["PGR"].values)
                log(f"  r(ESR1, PGR) within {lbl}: "
                    f"r={r:>+.4f}  p={p:.4e}")
            except Exception as e:
                log(f"  r(ESR1,PGR) error ({lbl}): {e}")

    log(f"\n  P6 — CONTROLS (expected absent):")
    for gene in ["KRT5", "KRT14", "SOX10", "CDX2", "NKX2-1"]:
        nm, tm, chg, p = test(gene, mature, lumb)
        if nm is None:
            log(f"  {gene:<12} not in data")
            continue
        flag = "ABSENT" if abs(chg or 0) < 30 else "ELEVATED-UNEXPECTED"
        log(f"  {flag} {gene:<10} {chg:>+7.1f}%  {fmt_p(p)}")

    log(f"\n  P7 — EPIGENETIC PANEL (LumB vs Mature):")
    log(f"  {'Gene':<12} {'Mature':>8} {'LumB':>8} {'%chg':>8}  "
        f"{'p':>14}  r_depth")
    log(f"  {'-'*58}")
    for gene in EPIGENETIC:
        nm, tm, chg, p = test(gene, mature, lumb)
        if nm is None:
            continue
        rv, _ = dr(gene)
        log(f"  {gene:<12} {nm:>8.4f} {tm:>8.4f} "
            f"{chg:>+7.1f}%  {fmt_p(p):>14}  {rv:>+.3f}")

    log(f"\n  P9 — CROSS-SUBTYPE GRADIENT:")
    log(f"  Predicted ESR1/FOXA1/GATA3: LumA > LumB > TNBC (DESC)")
    log(f"  Predicted EZH2:             LumA < LumB < TNBC (ASC)")
    log(f"  {'Gene':<12} {'LumA':>8} {'LumB':>8} {'TNBC':>8}  Gradient?")
    log(f"  {'-'*54}")
    for gene in GRADIENT_GENES:
        la = luma[gene].mean()  if gene in luma.columns  else np.nan
        lb = lumb[gene].mean()  if gene in lumb.columns  else np.nan
        tn = basal[gene].mean() if gene in basal.columns else np.nan
        if any(np.isnan(v) for v in [la, lb, tn]):
            log(f"  {gene:<12} MISSING DATA")
            continue
        if gene in ("ESR1", "FOXA1", "GATA3"):
            expected = "DESC"
            actual = ("DESC" if la > lb > tn else
                      "PART" if la > tn       else "FAIL")
        elif gene == "EZH2":
            expected = "ASC"
            actual = ("ASC"  if la < lb < tn else
                      "PART" if la < tn       else "FAIL")
        else:
            expected, actual = "?", "?"
        sym = "+" if actual in (expected, "PART") else "X"
        log(f"  {sym} {gene:<10} "
            f"LumA={la:.4f}  LumB={lb:.4f}  TNBC={tn:.4f}  -> {actual}")

# ============================================================
# SECTION 5: NOVEL SIGNALS
# ============================================================

def novel_signals(rdf_mature, corrs):
    log("")
    log("=" * 65)
    log("SECTION 5: NOVEL SIGNALS")
    log("Signals not in the locked panel (p<0.05).")
    log("=" * 65)

    if rdf_mature is None or rdf_mature.empty:
        log("  No top-movers data.")
        return

    panel = set(ALL_GENES)
    for label, direction in [
        ("Novel SUPPRESSED in LumB vs Mature", "DOWN"),
        ("Novel ELEVATED in LumB vs Mature",   "UP"),
    ]:
        subset = rdf_mature[
            (rdf_mature["direction"] == direction) &
            (~rdf_mature["gene"].isin(panel)) &
            (rdf_mature["p_value"] < 0.05)
        ].head(10)
        log(f"\n  {label}:")
        if subset.empty:
            log("    (none beyond panel)")
        for _, row in subset.iterrows():
            log(f"    {row['gene']:<20} {row['change_pct']:>+8.1f}%  "
                f"p={row['p_value']:.2e}")

    if corrs:
        novel_c = [(g, r, p) for g, r, p in corrs
                   if g not in panel and abs(r) > 0.15]
        if novel_c:
            log(f"\n  Novel depth correlates (|r|>0.15, not in panel):")
            for g, r, p in novel_c[:10]:
                log(f"    {g:<20} r={r:>+.4f}  p={p:.2e}")

# ============================================================
# STEP 7: BULK RNA-SEQ (optional)
# ============================================================

def bulk_analysis(bulk_path):
    log("")
    log("=" * 65)
    log("STEP 7: BULK RNA-SEQ ANALYSIS (OPTIONAL)")
    log("=" * 65)

    if bulk_path is None or not os.path.exists(bulk_path):
        log("  Bulk file not found. Step skipped.")
        return None

    try:
        with gzip.open(bulk_path, "rt") as f:
            bulk = pd.read_csv(f, sep="\t", index_col=0)
        log(f"  Bulk shape: {bulk.shape}")

        available = [g for g in ALL_GENES if g in bulk.index]
        log(f"  Target genes in bulk: {len(available)}")
        bulk_sub = bulk.loc[available]

        depth_gns = [g for g in ["CDKN1A", "FOXA1", "GATA3",
                                  "CDK4", "MKI67"] if g in available]
        if depth_gns:
            sub    = bulk_sub.loc[depth_gns]
            normed = sub.apply(
                lambda row: (row - row.min()) /
                            (row.max() - row.min() + 1e-8), axis=1)
            for g in ["CDKN1A", "FOXA1", "GATA3"]:
                if g in normed.index:
                    normed.loc[g] = 1 - normed.loc[g]
            bulk_depth = normed.mean(axis=0)
            log(f"\n  Bulk depth scores ({len(bulk_depth)} tumors):")
            log(f"    Mean={bulk_depth.mean():.4f}  "
                f"Std={bulk_depth.std():.4f}  "
                f"Min={bulk_depth.min():.4f}  "
                f"Max={bulk_depth.max():.4f}")
        return bulk_sub
    except Exception as e:
        log(f"  Bulk error: {e}")
        return None

# ============================================================
# FIGURE — 9-panel discovery output
# ============================================================

def generate_figure(pops, rdf_mature, rdf_luma,
                    depth, corrs, pca_result):
    log("")
    log("=" * 65)
    log("GENERATING FIGURE")
    log("=" * 65)

    lumb   = pops.get("lumb",   pd.DataFrame())
    luma   = pops.get("luma",   pd.DataFrame())
    mature = pops.get("mature", pd.DataFrame())
    basal  = pops.get("basal",  pd.DataFrame())
    pcs, labels = pca_result if pca_result else (None, None)

    CLR = {
        "lumb":   "#8e44ad",
        "luma":   "#2980b9",
        "mature": "#95a5a6",
        "basal":  "#c0392b",
        "prog":   "#27ae60",
        "up":     "#c0392b",
        "down":   "#2980b9",
        "flat":   "#27ae60",
    }

    def bar_clr(pct):
        if pct > 20:  return CLR["up"]
        if pct < -20: return CLR["down"]
        return CLR["flat"]

    fig = plt.figure(figsize=(24, 18))
    fig.suptitle(
        "BRCA Luminal B — Script 1 Discovery\n"
        "OrganismCore | GSE176078 | 2026-03-05\n"
        "Document BRCA-S5b | Type 3 + partial Type 1\n"
        f"LumB={len(lumb)}  LumA={len(luma)}  "
        f"Mature={len(mature)}  TNBC={len(basal)}",
        fontsize=9, fontweight="bold", y=1.005
    )
    gs = gridspec.GridSpec(3, 3, figure=fig,
                           hspace=0.55, wspace=0.44)

    # A — Identity TFs
    ax = fig.add_subplot(gs[0, 0])
    avail = [g for g in IDENTITY_TFS if g in lumb.columns]
    if avail:
        x, w = np.arange(len(avail)), 0.28
        ax.bar(x-w, [mature[g].mean() if g in mature.columns else 0
                     for g in avail], w,
               color=CLR["mature"], label="Mature", alpha=0.85)
        ax.bar(x,   [luma[g].mean()   if g in luma.columns   else 0
                     for g in avail], w,
               color=CLR["luma"],   label="LumA",   alpha=0.85)
        ax.bar(x+w, [lumb[g].mean() for g in avail], w,
               color=CLR["lumb"],   label="LumB",   alpha=0.85)
        ax.set_xticks(x)
        ax.set_xticklabels(avail, rotation=45, ha="right", fontsize=7)
        ax.set_title("A — Identity TFs\n(P1b: LumB < LumA?)",
                     fontsize=8, fontweight="bold")
        ax.legend(fontsize=6)
        ax.set_ylabel("log1p expression", fontsize=7)

    # B — Proliferation
    ax = fig.add_subplot(gs[0, 1])
    avail = [g for g in ["MKI67","CCND1","CDK4","MYC","TOP2A","CDKN1A"]
             if g in lumb.columns]
    if avail:
        x, w = np.arange(len(avail)), 0.28
        ax.bar(x-w, [mature[g].mean() if g in mature.columns else 0
                     for g in avail], w,
               color=CLR["mature"], label="Mature", alpha=0.85)
        ax.bar(x,   [luma[g].mean()   if g in luma.columns   else 0
                     for g in avail], w,
               color=CLR["luma"],   label="LumA",   alpha=0.85)
        ax.bar(x+w, [lumb[g].mean() for g in avail], w,
               color=CLR["lumb"],   label="LumB",   alpha=0.85)
        ax.set_xticks(x)
        ax.set_xticklabels(avail, rotation=45, ha="right", fontsize=7)
        ax.set_title("B — Proliferation\n(P2: LumB > LumA?)",
                     fontsize=8, fontweight="bold")
        ax.legend(fontsize=6)
        ax.set_ylabel("log1p expression", fontsize=7)

    # C — Waterfall: LumB vs Mature
    ax = fig.add_subplot(gs[0, 2])
    if rdf_mature is not None and len(rdf_mature) > 0:
        sig  = rdf_mature[rdf_mature["p_value"] < 0.05].sort_values(
            "change_pct")
        show = pd.concat([sig.head(12), sig.tail(12)]).drop_duplicates()
        ax.barh(show["gene"], show["change_pct"],
                color=[bar_clr(v) for v in show["change_pct"]])
        ax.axvline(0, color="black", lw=0.8)
        ax.set_title("C — Waterfall LumB vs Mature\n(unfiltered, p<0.05)",
                     fontsize=8, fontweight="bold")
        ax.set_xlabel("% change", fontsize=7)
        ax.tick_params(axis="y", labelsize=5)

    # D — Waterfall: LumB vs LumA (CRITICAL)
    ax = fig.add_subplot(gs[1, 0])
    if rdf_luma is not None and len(rdf_luma) > 0:
        sig  = rdf_luma[rdf_luma["p_value"] < 0.05].sort_values(
            "change_pct")
        show = pd.concat([sig.head(12), sig.tail(12)]).drop_duplicates()
        ax.barh(show["gene"], show["change_pct"],
                color=[bar_clr(v) for v in show["change_pct"]])
        ax.axvline(0, color="black", lw=0.8)
        ax.set_title("D — Waterfall LumB vs LumA\n(CRITICAL test)",
                     fontsize=8, fontweight="bold")
        ax.set_xlabel("% diff (LumB - LumA)", fontsize=7)
        ax.tick_params(axis="y", labelsize=5)

    # E — EZH2 gradient
    ax = fig.add_subplot(gs[1, 1])
    ezh2_pops = [
        ("Mature", mature, CLR["mature"]),
        ("LumA",   luma,   CLR["luma"]),
        ("LumB",   lumb,   CLR["lumb"]),
        ("TNBC",   basal,  CLR["basal"]),
    ]
    ezh2_vals = {lbl: df["EZH2"].mean()
                 for lbl, df, _ in ezh2_pops
                 if "EZH2" in df.columns and len(df) > 0}
    if ezh2_vals:
        lbls = list(ezh2_vals.keys())
        vals = list(ezh2_vals.values())
        cols = [c for lbl, _, c in ezh2_pops if lbl in ezh2_vals]
        bars = ax.bar(lbls, vals, color=cols, alpha=0.85)
        ax.set_title("E — EZH2 Gradient (P3)\n"
                     "Predicted: Mature < LumA < LumB < TNBC",
                     fontsize=8, fontweight="bold")
        ax.set_ylabel("log1p EZH2", fontsize=7)
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width()/2,
                    v + 0.002, f"{v:.4f}",
                    ha="center", fontsize=6)

    # F — Depth score distribution
    ax = fig.add_subplot(gs[1, 2])
    if depth is not None and len(depth) > 2:
        ax.hist(depth.values, bins=40,
                color=CLR["lumb"], alpha=0.75, edgecolor="white")
        ax.axvline(depth.mean(), color="black", lw=1.5, ls="--",
                   label=f"mean={depth.mean():.3f}")
        ax.axvline(depth.quantile(0.25), color="green", lw=1, ls=":",
                   label=f"Q25={depth.quantile(0.25):.3f}")
        ax.axvline(depth.quantile(0.75), color="red",   lw=1, ls=":",
                   label=f"Q75={depth.quantile(0.75):.3f}")
        ax.set_title("F — Depth Score Distribution\n(within LumB)",
                     fontsize=8, fontweight="bold")
        ax.set_xlabel("Depth Score", fontsize=7)
        ax.set_ylabel("n cells", fontsize=7)
        ax.legend(fontsize=6)

    # G — PCA landscape
    ax = fig.add_subplot(gs[2, 0])
    if pcs is not None and labels is not None:
        for lbl, ck, al in [
            ("LumB",       "lumb",   0.20),
            ("LumA",       "luma",   0.15),
            ("Mature",     "mature", 0.25),
            ("Progenitor", "prog",   0.25),
            ("TNBC",       "basal",  0.15),
        ]:
            m = labels == lbl
            if m.sum() > 0:
                ax.scatter(pcs[m, 0], pcs[m, 1],
                           c=CLR.get(ck, "#888888"),
                           alpha=al, s=3, label=lbl)
        ax.set_title("G — PCA Landscape\n(LumB intermediate?)",
                     fontsize=8, fontweight="bold")
        ax.set_xlabel("PC1", fontsize=7)
        ax.set_ylabel("PC2", fontsize=7)
        ax.legend(fontsize=6, markerscale=3)

    # H — Depth correlations
    ax = fig.add_subplot(gs[2, 1])
    if corrs:
        top20 = corrs[:20]
        gc    = [g for g, r, p in top20]
        vc    = [r for g, r, p in top20]
        ax.barh(gc, vc,
                color=[CLR["up"] if v > 0 else CLR["down"] for v in vc])
        ax.axvline(0, color="black", lw=0.8)
        ax.set_title("H — Depth Correlations\n(top 20 |r|)",
                     fontsize=8, fontweight="bold")
        ax.set_xlabel("r", fontsize=7)
        ax.tick_params(axis="y", labelsize=6)

    # I — Summary text
    ax = fig.add_subplot(gs[2, 2])
    ax.axis("off")

    ezh2_lumb_pct, ezh2_luma_pct = "?", "?"
    if ("EZH2" in lumb.columns and "EZH2" in mature.columns
            and mature["EZH2"].mean() > 1e-6):
        mm = mature["EZH2"].mean()
        ezh2_lumb_pct = f"{(lumb['EZH2'].mean()-mm)/mm*100:>+.1f}%"
        if "EZH2" in luma.columns:
            ezh2_luma_pct = f"{(luma['EZH2'].mean()-mm)/mm*100:>+.1f}%"

    summary = (
        f"I — BRCA-S5b SUMMARY\n{'─'*28}\n"
        f"GSE176078 | Wu et al. 2021\n"
        f"SELF-CONTAINED download\n"
        f"LumB={len(lumb)}  LumA={len(luma)}\n"
        f"Mature={len(mature)}  TNBC={len(basal)}\n\n"
        f"ATTRACTOR TYPE: T3 + partial T1\n\n"
        f"KEY TESTS:\n"
        f"  FOXA1 in LumB vs LumA?\n"
        f"  EZH2 gradient confirmed?\n\n"
        f"EZH2:\n"
        f"  LumA vs mature: {ezh2_luma_pct}\n"
        f"  LumB vs mature: {ezh2_lumb_pct}\n\n"
        f"DRUGS (BRCA-S5a):\n"
        f"  CDK4/6i (palbociclib)\n"
        f"  EZH2i if P3 confirmed\n"
        f"  Endocrine + CDK4/6i\n\n"
        f"READ ORDER:\n"
        f"  1. Sec 1b  LumB vs LumA\n"
        f"  2. Sec 1a  vs Mature\n"
        f"  3. Sec 2   PCA\n"
        f"  4. Sec 3   depth score\n"
        f"  5. Sec 4   predictions\n"
        f"  6. Panels D, E, G\n\n"
        f"OrganismCore | BRCA-S5b\n"
        f"2026-03-05"
    )
    ax.text(0.03, 0.97, summary,
            transform=ax.transAxes,
            fontsize=7, verticalalignment="top",
            fontfamily="monospace",
            bbox=dict(boxstyle="round",
                      facecolor="#f8f8f8",
                      edgecolor="#cccccc"))

    plt.savefig(FIG_FILE, dpi=150, bbox_inches="tight")
    log(f"  Figure saved: {FIG_FILE}")
    plt.close()

# ============================================================
# MAIN
# sc_dir is passed as a return value from acquire_data()
# through every downstream function. No global mutations.
# ============================================================

def main():
    log("=" * 65)
    log("BRCA LUMINAL B — SCRIPT 1")
    log("OrganismCore — Document BRCA-S5b")
    log("Date: 2026-03-05")
    log("=" * 65)
    log("")
    log("SELF-CONTAINED: All data downloaded to:")
    log(f"  {DATA_DIR}")
    log("No path crawling. No cache from other scripts.")
    log("")
    log("FIX vs v1: global SC_DIR removed from acquire_data().")
    log("  sc_dir is now a local variable returned to main().")
    log("  load_metadata() and load_sc_expression() take sc_dir")
    log("  as an argument — no module-level state mutation.")
    log("")

    # Step 0 — acquire data, get ACTUAL sc_dir back
    sc_dir, bulk_path = acquire_data()
    if sc_dir is None:
        log("FATAL: Data acquisition failed.")
        write_log()
        sys.exit(1)

    # Step 1 — metadata + label detection
    meta, strategy = load_metadata(sc_dir)
    if meta is None:
        log("FATAL: Metadata load failed.")
        write_log()
        sys.exit(1)

    # Step 2 — expression (sc_dir passed explicitly)
    expr = load_sc_expression(sc_dir)
    if expr is None:
        log("FATAL: Expression load failed.")
        write_log()
        sys.exit(1)

    # Step 3 — populations
    pops = classify_populations(expr, meta, strategy)
    if pops is None:
        log("FATAL: Population classification failed.")
        write_log()
        sys.exit(1)

    # Section 1 — top movers (geometry first)
    rdf_mature, rdf_luma, rdf_prog = top_movers_unfiltered(pops)

    # Section 2 — PCA
    pca_result = pca_geometry(pops)

    # Section 3 — depth score
    depth, corrs = compute_depth_score(pops)

    # Section 4 — prediction check
    prediction_check(pops, depth, corrs, rdf_luma)

    # Section 5 — novel signals
    novel_signals(rdf_mature, corrs)

    # Step 7 — bulk (optional, bulk_path from acquire_data)
    bulk_analysis(bulk_path)

    # Figure
    generate_figure(pops, rdf_mature, rdf_luma,
                    depth, corrs, pca_result)

    write_log()

    log("")
    log("=" * 65)
    log("SCRIPT 1 COMPLETE")
    log(f"  Data dir : {DATA_DIR}")
    log(f"  Log      : {LOG_FILE}")
    log(f"  Figure   : {FIG_FILE}")
    log(f"  CSV      : {CSV_FILE}")
    log("")
    log("READ ORDER (Protocol v2.0):")
    log("  1. Step 1     — strategy used (DIRECT vs PROXY)")
    log("  2. Section 1b — LumB vs LumA (FOXA1 result first)")
    log("  3. Section 1a — LumB vs Mature (absolute state)")
    log("  4. Section 2  — PCA position")
    log("  5. Section 3  — depth score distribution")
    log("  6. Section 4  — prediction check")
    log("  7. Section 5  — novel signals")
    log("  8. Figure E   �� EZH2 gradient bar chart")
    log("     Figure D   — LumB vs LumA waterfall")
    log("     Figure G   — PCA landscape")
    log("")
    log("THEN write BRCA-S5b reasoning artifact.")
    log("  First sentence: FOXA1 result from Section 1b.")
    log("  Second sentence: EZH2 position from panel E.")
    log("=" * 65)


if __name__ == "__main__":
    main()
