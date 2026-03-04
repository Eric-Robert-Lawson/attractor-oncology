"""
BRCA TNBC — SCRIPT 1
OrganismCore — Document BRCA-S4a/b | 2026-03-04

KEY FIX vs previous version:
  GSE176078 has TWO different archives on GEO:
    1. GSE176078_RAW.tar          ← 26 per-sample tar.gz files
                                     DO NOT USE THIS ONE
    2. GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz
                                  ← pre-merged matrix (4 files)
                                     THIS IS THE CORRECT ONE
  The LumA script used #2. This script uses #2.
  The correct URL is on the GEO FTP suppl directory.
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
# CONFIGURATION
# ============================================================

GEO_ACCESSION = "GSE176078"

# --- Path discovery: look for data in common locations ---
# The original BRCA analysis may have already downloaded
# the files. We search up the directory tree for them.
SEARCH_PATHS = [
    "./",
    "../",                              # DEEP_DIVE/Luminal_A/ → DEEP_DIVE/
    "../../",                           # → BRCA/
    "../../../",                        # → Cancer_Research/
    "../../../../",
    os.path.expanduser("~/cancer/BRCA/"),
    os.path.expanduser("~/cancer/"),
    f"./{GEO_ACCESSION}/",
    f"../BRCA/{GEO_ACCESSION}/",
    f"../../{GEO_ACCESSION}/",
]

SC_TARFILE   = "GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"
BULK_FILE    = "GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts.txt.gz"

SC_URL   = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/"
    "GSE176078/suppl/"
    "GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"
)
BULK_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/"
    "GSE176078/suppl/"
    "GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts.txt.gz"
)

# Output directory — always local to this script
SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
BASE_DIR    = "./TNBC_s1_analysis/"
DATA_DIR    = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
LOG_FILE    = os.path.join(RESULTS_DIR, "tnbc_s1_log.txt")
FIG_FILE    = os.path.join(RESULTS_DIR, "tnbc_s1_figure.png")
CSV_FILE    = os.path.join(RESULTS_DIR, "tnbc_s1_saddle.csv")
CACHE_FILE  = os.path.join(RESULTS_DIR, "expr_cache_tnbc_s1.csv")

for d in [BASE_DIR, DATA_DIR, RESULTS_DIR]:
    os.makedirs(d, exist_ok=True)

SC_INNER_DIR = "Wu_etal_2021_BRCA_scRNASeq"

# Known filenames INSIDE the TAR (discovered from LumA run)
SC_FILES_INSIDE = [
    "count_matrix_sparse.mtx",
    "count_matrix_genes.tsv",
    "count_matrix_barcodes.tsv",
    "metadata.csv",
]

# GSE25066 — bulk pCR dataset (separate from GSE176078)
GEO_BULK    = "GSE25066"
BULK_URL    = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/"
    "GSE25nnn/GSE25066/matrix/"
    "GSE25066_series_matrix.txt.gz"
)
BULK_FILE   = "GSE25066_series_matrix.txt.gz"

# ============================================================
# CELL TYPE LABELS (exact strings from Wu et al. metadata)
# Confirmed from LumA log: CT_COL = "celltype_major"
# ============================================================

CT_COL       = "celltype_major"
CANCER_BASAL = "Cancer Basal SC"
MATURE_LUM   = "Normal Luminal Mature"
LUMINAL_PROG = "Normal Luminal Progenitors"
MYOEPITH     = "Normal Myoepithelial"
CANCER_LUMA  = "Cancer LumA SC"

# ============================================================
# GENE PANELS — Locked in BRCA-S4a (2026-03-04)
# ============================================================

SWITCH_GENES = [
    "ESR1", "FOXA1", "GATA3", "SPDEF", "PGR",
    "KRT8", "KRT18", "BRCA1",
]
FA_MARKERS = [
    "KRT5", "KRT14", "SOX10", "FOXC1", "EGFR",
    "VIM", "CDH3",
]
EPIGENETIC = [
    "EZH2", "HDAC1", "HDAC2", "KDM1A",
    "TET2", "DNMT3A", "ASXL1", "EED",
]
PROLIFERATION = [
    "MKI67", "TOP2A", "PCNA", "CCNB1", "CDK2",
]
COMPOSITE_TYPE = [
    "BRCA1", "BRCA2", "TP53", "PTEN", "RB1",
    "PIK3CA", "AKT1",
]
DRUG_TARGETS = [
    "EZH2", "PARP1", "AURKA", "CD274",
]
HETEROGENEITY = [
    "AR", "CDH1",
    "ZEB1", "ZEB2", "SNAI1",
    "CLDN3", "CLDN4",
]
CONTROLS = [
    "CDX2", "SPI1", "NKX2-1", "NKX3-1", "OLIG2",
]

ALL_GENES = list(dict.fromkeys(
    SWITCH_GENES + FA_MARKERS + EPIGENETIC +
    PROLIFERATION + COMPOSITE_TYPE + DRUG_TARGETS +
    HETEROGENEITY + CONTROLS
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
# FILE UTILITIES
# ============================================================

def find_file(filename):
    """Search SEARCH_PATHS for filename. Returns path or None."""
    for base in SEARCH_PATHS:
        candidate = os.path.join(base, filename)
        if os.path.exists(candidate) and os.path.getsize(candidate) > 10000:
            log(f"  Found: {candidate} "
                f"({os.path.getsize(candidate)/1e6:.1f} MB)")
            return os.path.abspath(candidate)
    return None


def find_extracted_sc_dir():
    """
    Search for the already-extracted Wu_etal_2021_BRCA_scRNASeq/
    directory containing the four MTX files.
    """
    for base in SEARCH_PATHS:
        # Standard subdirectory
        candidate = os.path.join(base, SC_INNER_DIR)
        if (os.path.isdir(candidate) and
                os.path.exists(os.path.join(
                    candidate, "count_matrix_sparse.mtx"))):
            log(f"  Found extracted dir: {candidate}")
            return os.path.abspath(candidate)

        # Files extracted directly into base
        if os.path.exists(os.path.join(
                base, "count_matrix_sparse.mtx")):
            log(f"  Found sc files directly in: {base}")
            return os.path.abspath(base)

        # One level deeper (sub-extraction)
        if os.path.isdir(base):
            try:
                for sub in os.listdir(base):
                    sub_path = os.path.join(base, sub)
                    if (os.path.isdir(sub_path) and
                            os.path.exists(os.path.join(
                                sub_path,
                                "count_matrix_sparse.mtx"))):
                        log(f"  Found sc files in subdir: {sub_path}")
                        return os.path.abspath(sub_path)
            except PermissionError:
                pass
    return None


def fetch_url(url, dest, retries=3):
    for attempt in range(retries):
        try:
            log(f"  Attempt {attempt+1}: {url[-70:]}")
            req = urllib.request.Request(
                url, headers={"User-Agent": "Mozilla/5.0"}
            )
            with urllib.request.urlopen(req, timeout=300) as r:
                data = r.read()
            with open(dest, "wb") as f:
                f.write(data)
            log(f"  Downloaded: {dest} ({len(data)/1e6:.1f} MB)")
            return dest
        except Exception as e:
            log(f"  Attempt {attempt+1} failed: {e}")
            if attempt < retries - 1:
                time.sleep(4)
    return None

# ============================================================
# STEP 0: DATA ACQUISITION
# Uses GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz
# NOT GSE176078_RAW.tar (which contains per-sample tarballs)
# ============================================================

def acquire_data():
    log("=" * 65)
    log("STEP 0: DATA ACQUISITION")
    log(f"Dataset: {GEO_ACCESSION} — Wu et al. 2021")
    log(f"File: {SC_TARFILE}")
    log("NOTE: Using pre-merged matrix TAR, NOT _RAW.tar")
    log("=" * 65)

    # --- scRNA-seq ---
    log("\n--- scRNA-seq (GSE176078 merged matrix) ---")
    sc_dir = find_extracted_sc_dir()

    if sc_dir is None:
        log("  Extracted files not found. Looking for TAR...")
        tar_path = find_file(SC_TARFILE)

        if tar_path is None:
            log(f"  TAR not found. Downloading (~2 GB)...")
            tar_dest = os.path.join(DATA_DIR, SC_TARFILE)
            tar_path = fetch_url(SC_URL, tar_dest)
            if tar_path is None:
                log("  FATAL: Download failed.")
                return None, None

        log(f"  Extracting: {os.path.basename(tar_path)}")
        extract_dir = os.path.join(DATA_DIR, "sc_extracted")
        os.makedirs(extract_dir, exist_ok=True)

        try:
            with tarfile.open(tar_path, "r:gz") as tf:
                members = tf.getnames()
                log(f"  TAR contents ({len(members)} items):")
                for m in members[:8]:
                    log(f"    {m}")
                tf.extractall(extract_dir)
            log("  Extraction complete.")
        except Exception as e:
            log(f"  TAR extraction error: {e}")
            return None, None

        sc_dir = find_extracted_sc_dir()
        if sc_dir is None:
            log("  FATAL: Cannot locate extracted sc files.")
            log(f"  Contents of {extract_dir}:")
            for root, dirs, files in os.walk(extract_dir):
                for fname in files[:5]:
                    log(f"    {os.path.join(root, fname)}")
            return None, None

    log(f"  scRNA-seq dir: {sc_dir}")

    # Verify required files
    for fname in SC_FILES_INSIDE:
        fp = os.path.join(sc_dir, fname)
        if not os.path.exists(fp):
            log(f"  MISSING: {fname}")
            return None, None
        log(f"  OK: {fname} ({os.path.getsize(fp)/1e6:.1f} MB)")

    # --- Bulk pCR data ---
    log("\n--- Bulk pCR data (GSE25066) ---")
    bulk_dest = os.path.join(DATA_DIR, BULK_FILE)
    bulk_path = find_file(BULK_FILE)
    if bulk_path is None:
        log("  Downloading GSE25066 series matrix...")
        bulk_path = fetch_url(BULK_URL, bulk_dest)
        if bulk_path is None:
            log("  Download failed. P6 (pCR test) will be skipped.")

    return sc_dir, bulk_path

# ============================================================
# STEP 1: LOAD METADATA
# ============================================================

def load_metadata(sc_dir):
    log("=" * 65)
    log("STEP 1: LOAD METADATA")
    log("=" * 65)

    meta = pd.read_csv(os.path.join(sc_dir, "metadata.csv"),
                       index_col=0)
    log(f"  Cells: {len(meta)}")
    log(f"  Columns: {list(meta.columns[:10])}")

    # Print value distributions for all three celltype columns
    for col in ["celltype_major", "celltype_minor", "celltype_subset"]:
        if col in meta.columns:
            log(f"\n  Cell type distribution ({col}):")
            for ct, n in meta[col].value_counts().items():
                log(f"    {ct:<40} n={n}")

    # Auto-detect the correct fine-grained column by checking which one
    # contains the expected fine-grained labels. Try in priority order.
    _fine_labels = {CANCER_BASAL, MATURE_LUM}
    _candidate_cols = ["celltype_minor", "celltype_subset", "celltype_major"]

    def _build_combined(df, cols):
        combined = pd.Series("", index=df.index)
        for c in cols:
            if c in df.columns:
                combined = combined + " " + df[c].fillna("")
        return combined

    ct_col = None
    for col in _candidate_cols:
        if col in meta.columns:
            vals = set(meta[col].unique())
            if _fine_labels & vals:
                ct_col = col
                log(f"\n  Auto-detected fine-grained column: {col}")
                break

    if ct_col is None:
        # Partial-match fallback: combine all three celltype columns into one
        # search column and look for cells matching "Cancer" AND "Basal"
        log(f"\n  WARNING: No column contains exact fine-grained labels.")
        log(f"  Falling back to partial match across all celltype columns...")
        meta["_ct_combined"] = _build_combined(meta, _candidate_cols)
        ct_col = "_ct_combined"

    if ct_col not in meta.columns:
        log("  FATAL: Cannot find cell type column.")
        return None

    # Tag populations
    meta["is_tnbc"]   = meta[ct_col] == CANCER_BASAL
    meta["is_luma"]   = meta[ct_col] == CANCER_LUMA
    meta["is_mature"] = meta[ct_col] == MATURE_LUM
    meta["is_prog"]   = meta[ct_col] == LUMINAL_PROG
    meta["is_myo"]    = meta[ct_col] == MYOEPITH

    n_tnbc = meta["is_tnbc"].sum()
    log(f"\n  KEY POPULATIONS:")
    log(f"    TNBC / Basal SC:       {n_tnbc}")
    log(f"    LumA SC:               {meta['is_luma'].sum()}")
    log(f"    Mature Luminal:        {meta['is_mature'].sum()}")
    log(f"    Luminal Progenitors:   {meta['is_prog'].sum()}")
    log(f"    Myoepithelial:         {meta['is_myo'].sum()}")

    if n_tnbc < 50:
        log(f"  WARNING: Only {n_tnbc} TNBC cells with exact label.")
        log(f"  Trying partial match 'Basal' + 'Cancer'...")
        combined = _build_combined(meta, _candidate_cols)
        meta["is_tnbc"] = (
            combined.str.contains("Basal", case=False, na=False) &
            combined.str.contains("Cancer", case=False, na=False)
        )
        log(f"  After partial match: {meta['is_tnbc'].sum()}")

    return meta

# ============================================================
# STEP 2: LOAD EXPRESSION
# Reads uncompressed MTX directly (no gzip — files inside
# the tar.gz are stored uncompressed)
# ============================================================

def load_sc_expression(sc_dir, meta):
    log("=" * 65)
    log("STEP 2: LOAD EXPRESSION")
    log("Reading MTX via scipy.io.mmread (same as LumA)")
    log("=" * 65)

    # Use TNBC-specific cache — never touches LumA cache
    if os.path.exists(CACHE_FILE) and os.path.getsize(CACHE_FILE) > 100:
        sz = os.path.getsize(CACHE_FILE) / 1e6
        log(f"  Loading TNBC cache: {CACHE_FILE} ({sz:.1f} MB)")
        expr = pd.read_csv(CACHE_FILE, index_col=0)
        log(f"  Cached shape: {expr.shape}")
        return expr

    mtx_file  = os.path.join(sc_dir, "count_matrix_sparse.mtx")
    gene_file = os.path.join(sc_dir, "count_matrix_genes.tsv")
    bc_file   = os.path.join(sc_dir, "count_matrix_barcodes.tsv")

    # Gene names
    genes_df   = pd.read_csv(gene_file, sep="\t", header=None)
    gene_names = genes_df.iloc[:, 0].tolist()
    log(f"  Genes in matrix: {len(gene_names)}")
    log(f"  Sample gene names: {gene_names[:5]}")

    # Barcodes
    barcodes = pd.read_csv(
        bc_file, sep="\t", header=None
    ).iloc[:, 0].tolist()
    log(f"  Barcodes: {len(barcodes)}")

    # Load full sparse matrix via scipy (same as LumA)
    log("  Reading MTX via scipy.io.mmread (~30-90 sec)...")
    mat = mmread(mtx_file).tocsr()
    log(f"  Matrix shape: {mat.shape}  (genes × cells)")

    # Map target genes to row indices
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

    log(f"  Target genes found: {len(target_names)} / {len(ALL_GENES)}")
    if missing:
        log(f"  Missing: {missing}")

    # Extract target gene rows
    log("  Extracting target genes from sparse matrix...")
    expr_data = {}
    for i, (row_idx, gene_name) in enumerate(
        zip(target_idx, target_names)
    ):
        expr_data[gene_name] = mat[row_idx, :].toarray().flatten()
        if (i + 1) % 10 == 0:
            log(f"    {i+1}/{len(target_names)} genes extracted")

    expr = pd.DataFrame(expr_data, index=barcodes)
    log(f"  Raw expression matrix: {expr.shape}")

    # Library size check (Protocol v2.0)
    lib_sizes = expr.sum(axis=1)
    med       = lib_sizes.median()
    log(f"\n  Library sizes (raw counts, target genes):")
    log(f"    Median: {med:.1f}  Mean: {lib_sizes.mean():.1f}  "
        f"Min: {lib_sizes.min():.1f}  Max: {lib_sizes.max():.1f}")
    outliers = lib_sizes[
        (lib_sizes < med / 5) | (lib_sizes > med * 5)
    ]
    if len(outliers) > 0:
        log(f"  Library size outliers (outside 5x median): {len(outliers)}")

    # Log1p normalize
    expr = np.log1p(expr)
    log("  Log1p normalization applied")

    # Spot check
    for gene in (SWITCH_GENES + FA_MARKERS)[:6]:
        if gene in expr.columns:
            log(f"    {gene:<12} mean={expr[gene].mean():.4f}  "
                f"nonzero={(expr[gene]>0).sum()}")

    # Save TNBC-specific cache
    expr.to_csv(CACHE_FILE)
    log(f"\n  TNBC cache saved: {CACHE_FILE}")
    return expr

# ============================================================
# STEP 3: CLASSIFY POPULATIONS
# ============================================================

def classify_populations(expr, meta):
    log("=" * 65)
    log("STEP 3: CLASSIFY POPULATIONS")
    log("=" * 65)

    common = expr.index.intersection(meta.index)
    log(f"  Expr/meta overlap: {len(common)}")

    expr_c = expr.loc[common]
    meta_c = meta.loc[common]

    pops = {}
    for name, flag in [
        ("tnbc",   "is_tnbc"),
        ("mature", "is_mature"),
        ("prog",   "is_prog"),
        ("myo",    "is_myo"),
        ("luma",   "is_luma"),
    ]:
        idx = meta_c.index[meta_c[flag]]
        pops[name] = expr_c.loc[idx]
        log(f"  {name:<10}: {len(pops[name])} cells")

    if len(pops["tnbc"]) < 50:
        log("  ERROR: Insufficient TNBC cells. Check cell type labels.")
        return None

    return pops

# ============================================================
# STEP 4: LOAD BULK pCR DATA (GSE25066)
# ============================================================

def load_bulk_pcr(bulk_path):
    log("=" * 65)
    log("STEP 4: LOAD BULK pCR DATA (GSE25066)")
    log("=" * 65)

    if bulk_path is None or not os.path.exists(bulk_path):
        log("  Not available. P6 will be skipped.")
        return None, None

    sample_ids, pcr_v, er_v, pr_v, her2_v = [], {}, {}, {}, {}
    expr_rows, col_ids = {}, []

    try:
        open_fn = gzip.open if bulk_path.endswith(".gz") else open
        with open_fn(bulk_path, "rt", errors="ignore") as f:
            in_table = False
            for line in f:
                line = line.rstrip("\n")
                if line.startswith("!Sample_geo_accession"):
                    sample_ids = [
                        s.strip().strip('"')
                        for s in line.split("\t")[1:]
                    ]
                elif line.startswith("!Sample_characteristics_ch1"):
                    parts = line.split("\t")
                    v0 = parts[1].strip().strip('"').lower() if len(parts) > 1 else ""
                    if any(k in v0 for k in ["pcr", "prc", "pathological", "residual"]):
                        for i, p in enumerate(parts[1:]):
                            v = p.strip().strip('"').lower()
                            if i < len(sample_ids):
                                pcr_v[sample_ids[i]] = (
                                    1 if ("complete" in v or "pcr" in v)
                                    else 0
                                )
                    elif "er status" in v0 or v0.startswith("er:"):
                        for i, p in enumerate(parts[1:]):
                            if i < len(sample_ids):
                                er_v[sample_ids[i]] = (
                                    1 if "pos" in p.lower() else 0
                                )
                    elif "pr status" in v0 or v0.startswith("pr:"):
                        for i, p in enumerate(parts[1:]):
                            if i < len(sample_ids):
                                pr_v[sample_ids[i]] = (
                                    1 if "pos" in p.lower() else 0
                                )
                    elif "her2" in v0:
                        for i, p in enumerate(parts[1:]):
                            if i < len(sample_ids):
                                her2_v[sample_ids[i]] = (
                                    1 if "pos" in p.lower() else 0
                                )
                elif "series_matrix_table_begin" in line:
                    in_table = True
                    hdr = f.readline().rstrip("\n")
                    col_ids = [c.strip().strip('"') for c in hdr.split("\t")[1:]]
                elif "series_matrix_table_end" in line:
                    break
                elif in_table and line:
                    parts = line.split("\t")
                    if len(parts) < 2:
                        continue
                    probe = parts[0].strip().strip('"')
                    try:
                        vals = [float(v.strip().strip('"')) for v in parts[1:]]
                        expr_rows[probe] = vals
                    except ValueError:
                        pass

        if not expr_rows:
            log("  No expression data parsed.")
            return None, None

        use_cols  = col_ids if col_ids else sample_ids
        n_vals    = len(list(expr_rows.values())[0])
        expr_bulk = pd.DataFrame.from_dict(
            expr_rows, orient="index",
            columns=use_cols[:n_vals]
        )
        log(f"  Bulk matrix: {expr_bulk.shape}")

        meta_bulk = pd.DataFrame({
            "pCR": pd.Series(pcr_v),
            "ER":  pd.Series(er_v),
            "PR":  pd.Series(pr_v),
            "HER2":pd.Series(her2_v),
        })

        if meta_bulk["ER"].notna().sum() > 10:
            tnbc_m = (
                (meta_bulk["ER"]   == 0) &
                (meta_bulk["PR"]   == 0) &
                (meta_bulk["HER2"] == 0)
            )
            meta_tnbc = meta_bulk[tnbc_m]
            log(f"  TNBC subset: {len(meta_tnbc)}")
            log(f"  pCR=1: {meta_tnbc['pCR'].sum():.0f}  "
                f"pCR=0: {(meta_tnbc['pCR']==0).sum()}")
        else:
            log("  ER/PR/HER2 not parseable. Using all samples.")
            meta_tnbc = meta_bulk

        return expr_bulk, meta_tnbc

    except Exception as e:
        log(f"  Bulk parse error: {e}")
        import traceback
        log(traceback.format_exc())
        return None, None

# ============================================================
# SECTION 1: TOP MOVERS (Protocol v2.0 — geometry first)
# ============================================================

def top_movers_unfiltered(pops):
    log("")
    log("=" * 65)
    log("SECTION 1: TOP MOVERS — UNFILTERED")
    log("Protocol v2.0: geometry reveals itself first")
    log("No prediction panel imposed. Read before Section 4.")
    log("=" * 65)

    tnbc   = pops.get("tnbc",   pd.DataFrame())
    mature = pops.get("mature", pd.DataFrame())

    if len(tnbc) < 10 or len(mature) < 3:
        log("  Insufficient cells.")
        return pd.DataFrame()

    def fmt_p(p):
        if p < 1e-100: return "p<1e-100 ***"
        if p < 0.001:  return f"p={p:.2e} ***"
        if p < 0.01:   return f"p={p:.2e}  **"
        if p < 0.05:   return f"p={p:.4f}   *"
        return             f"p={p:.4f}  ns"

    results = []
    for gene in tnbc.columns:
        nv = mature[gene].values if gene in mature.columns else np.array([])
        tv = tnbc[gene].values
        if len(nv) < 3 or len(tv) < 5:
            continue
        nm = nv.mean()
        tm = tv.mean()
        if abs(nm) < 1e-6 and abs(tm) < 1e-6:
            continue
        chg = (tm - nm) / abs(nm) * 100 if abs(nm) > 1e-6 else np.nan
        try:
            _, p = stats.mannwhitneyu(tv, nv, alternative="two-sided")
        except Exception:
            p = 1.0
        results.append({
            "gene":        gene,
            "normal_mean": nm,
            "tumor_mean":  tm,
            "change_pct":  chg,
            "abs_change":  abs(chg) if not np.isnan(chg) else 0,
            "p_value":     p,
            "direction":   "UP" if (chg or 0) > 0 else "DOWN",
        })

    rdf = pd.DataFrame(results).sort_values("abs_change", ascending=False)

    for label, direction in [
        ("LOST IN TNBC vs mature luminal", "DOWN"),
        ("GAINED IN TNBC vs mature luminal", "UP"),
    ]:
        subset = rdf[rdf["direction"] == direction].head(20)
        log(f"\n  TOP 20 GENES {label}")
        log(f"  {'Gene':<14} {'Normal':>9} {'TNBC':>9} "
            f"{'Change%':>9}  {'p-value':>16}")
        log(f"  {'-'*62}")
        for _, row in subset.iterrows():
            log(f"  {row['gene']:<14} {row['normal_mean']:>9.4f} "
                f"{row['tumor_mean']:>9.4f} {row['change_pct']:>+8.1f}%  "
                f"{fmt_p(row['p_value']):>16}")

    rdf.to_csv(CSV_FILE, index=False)
    log(f"\n  Saved: {CSV_FILE}")
    return rdf

# ============================================================
# SECTION 2: PCA GEOMETRY
# ============================================================

def pca_geometry(pops):
    log("")
    log("=" * 65)
    log("SECTION 2: PCA GEOMETRY")
    log("Read PC1 loadings before Section 4 (prediction check)")
    log("=" * 65)

    frames, labels = [], []
    for name, label in [
        ("tnbc",   "TNBC"),
        ("mature", "Mature"),
        ("prog",   "Progenitor"),
        ("myo",    "Myoepithelial"),
    ]:
        df = pops.get(name, pd.DataFrame())
        if len(df) >= 5:
            frames.append(df)
            labels.extend([label] * len(df))

    if not frames:
        log("  Insufficient populations.")
        return None, None

    combined = pd.concat(frames).fillna(0)
    labels   = np.array(labels[:len(combined)])

    # Outlier removal (Protocol v2.0)
    lib   = combined.sum(axis=1)
    med   = lib.median()
    omask = (lib < med / 5) | (lib > med * 5)
    if omask.sum() > 0:
        log(f"  Library size outliers removed: {omask.sum()}")
        combined = combined[~omask]
        labels   = labels[~omask]

    scaler = StandardScaler()
    scaled = scaler.fit_transform(combined)
    n_comp = min(5, combined.shape[1], combined.shape[0] - 1)
    pca    = PCA(n_components=n_comp)
    pcs    = pca.fit_transform(scaled)

    log(f"\n  Variance explained:")
    for i, v in enumerate(pca.explained_variance_ratio_):
        log(f"    PC{i+1}: {v*100:.2f}%")

    loadings = pd.Series(
        pca.components_[0], index=combined.columns
    ).sort_values(ascending=False)

    log(f"\n  PC1 top positive loadings:")
    for g, v in loadings.head(10).items():
        log(f"    {g:<18} {v:>+8.4f}")

    log(f"\n  PC1 top negative loadings:")
    for g, v in loadings.tail(10).items():
        log(f"    {g:<18} {v:>+8.4f}")

    pc1_t = pcs[labels == "TNBC",   0]
    pc1_m = pcs[labels == "Mature", 0]
    if len(pc1_t) > 5 and len(pc1_m) > 5:
        _, p_sep = stats.mannwhitneyu(pc1_t, pc1_m, alternative="two-sided")
        q = "STRONG" if p_sep < 0.001 else ("MODERATE" if p_sep < 0.05 else "WEAK")
        log(f"\n  PC1 separation TNBC vs Mature Luminal: {q}  p={p_sep:.2e}")
        log(f"    TNBC mean={pc1_t.mean():.4f}  Mature mean={pc1_m.mean():.4f}")

    return pcs, labels

# ============================================================
# SECTION 3: DEPTH SCORE
# ============================================================

def compute_depth(pops):
    log("")
    log("=" * 65)
    log("SECTION 3: DEPTH SCORE")
    log("Predicted = norm(basal_FAs) + (1-norm(luminal_SWs)) / 2")
    log("=" * 65)

    tnbc = pops.get("tnbc", pd.DataFrame())
    if len(tnbc) < 10:
        log("  Insufficient cells.")
        return None, []

    basal_g   = [g for g in ["KRT5", "KRT14", "SOX10", "FOXC1", "EGFR"]
                 if g in tnbc.columns]
    luminal_g = [g for g in ["ESR1", "FOXA1", "GATA3", "SPDEF"]
                 if g in tnbc.columns]

    log(f"  Basal component:   {basal_g}")
    log(f"  Luminal component: {luminal_g}")

    def norm01(s):
        mn, mx = s.min(), s.max()
        return (s - mn) / (mx - mn) if mx > mn else pd.Series(0.5, index=s.index)

    parts = []
    if basal_g:
        parts.append(norm01(tnbc[basal_g].mean(axis=1)))
    if luminal_g:
        parts.append(1 - norm01(tnbc[luminal_g].mean(axis=1)))
    if not parts:
        log("  No genes available.")
        return None, []

    depth = sum(parts) / len(parts)

    log(f"\n  Depth score ({len(tnbc)} TNBC cells):")
    log(f"    Mean={depth.mean():.4f}  Median={depth.median():.4f}  "
        f"Std={depth.std():.4f}  Min={depth.min():.4f}  Max={depth.max():.4f}")

    log(f"\n  Top depth correlations:")
    log(f"  {'Gene':<18} {'r':>8}  p-value")
    log(f"  {'-'*38}")
    corrs = []
    for gene in tnbc.columns:
        try:
            rv, pv = stats.pearsonr(depth.values, tnbc[gene].values)
            corrs.append((gene, rv, pv))
        except Exception:
            pass
    corrs.sort(key=lambda x: abs(x[1]), reverse=True)
    for g, r, p in corrs[:25]:
        log(f"  {g:<18} {r:>+8.4f}  p={p:.2e}")

    return depth, corrs

# ============================================================
# SECTION 4: PREDICTION PANEL CHECK
# ============================================================

def prediction_check(pops, depth, corrs):
    log("")
    log("=" * 65)
    log("SECTION 4: PREDICTION PANEL CHECK")
    log("Locked in BRCA-S4a 2026-03-04. Checked here.")
    log("=" * 65)

    tnbc   = pops.get("tnbc",   pd.DataFrame())
    mature = pops.get("mature", pd.DataFrame())
    prog   = pops.get("prog",   pd.DataFrame())

    def fmt_p(p):
        if p < 1e-100: return "p<1e-100 ***"
        if p < 0.001:  return f"p={p:.2e} ***"
        if p < 0.01:   return f"p={p:.2e}  **"
        if p < 0.05:   return f"p={p:.4f}   *"
        return             f"p={p:.4f}  ns"

    def test(gene, ref, cancer):
        if gene not in ref.columns or gene not in cancer.columns:
            return None, None, None, None
        nv = ref[gene].values
        tv = cancer[gene].values
        nm, tm = nv.mean(), tv.mean()
        chg = (tm - nm) / abs(nm) * 100 if abs(nm) > 1e-6 else np.nan
        try:
            _, p = stats.mannwhitneyu(tv, nv, alternative="two-sided")
        except Exception:
            p = 1.0
        return nm, tm, chg, p

    def dr(gene):
        if depth is None or gene not in tnbc.columns:
            return np.nan, np.nan
        try:
            return stats.pearsonr(depth.values, tnbc[gene].values)
        except Exception:
            return np.nan, np.nan

    for label, panel, pred_dir in [
        ("P1 — SWITCH GENES (predicted DOWN)", SWITCH_GENES, "DOWN"),
        ("P2 — FA MARKERS   (predicted UP)",   FA_MARKERS,   "UP"),
    ]:
        log(f"\n  {label}:")
        log(f"  {'Gene':<12} {'Mature':>8} {'TNBC':>8} {'%chg':>8}  "
            f"{'p':>14}  {'r_depth':>8}  Verdict")
        log(f"  {'-'*72}")
        for gene in panel:
            nm, tm, chg, p = test(gene, mature, tnbc)
            if nm is None:
                log(f"  {gene:<12} NOT IN DATA")
                continue
            rv, _ = dr(gene)
            actual = ("DOWN" if (chg or 0) < -10
                      else ("UP" if (chg or 0) > 10 else "FLAT"))
            if actual == pred_dir and p < 0.05:
                verdict = "✓ CONFIRMED"
            elif actual != pred_dir and actual != "FLAT":
                verdict = "✗ INVERTED"
            else:
                verdict = "? WEAK"
            log(f"  {gene:<12} {nm:>8.4f} {tm:>8.4f} "
                f"{chg:>+7.1f}%  {fmt_p(p):>14}  "
                f"{rv:>+7.3f}   {verdict}")

    log(f"\n  P3 — EZH2 CONVERGENCE NODE (UP, r>+0.30):")
    nm, tm, chg, p = test("EZH2", mature, tnbc)
    rv, _ = dr("EZH2")
    if nm is not None:
        v = ("✓ CONFIRMED" if (chg or 0) > 20 and p < 0.05 and rv > 0.30
             else "? PARTIAL")
        log(f"  EZH2: mature={nm:.4f} TNBC={tm:.4f} "
            f"{chg:>+.1f}%  {fmt_p(p)}  r={rv:>+.3f}  {v}")

    log(f"\n  P4 — COMPOSITE TYPE TEST (r(BRCA1, ESR1) > +0.15):")
    if "BRCA1" in tnbc.columns and "ESR1" in tnbc.columns:
        try:
            rv, pv = stats.pearsonr(tnbc["BRCA1"].values, tnbc["ESR1"].values)
            v = "✓ CONFIRMED" if rv > 0.15 and pv < 0.05 else "✗ NOT CONFIRMED"
            log(f"  r(BRCA1, ESR1) = {rv:+.4f}  p={pv:.2e}  {v}")
        except Exception as e:
            log(f"  Error: {e}")

    log(f"\n  P7 — EPIGENETIC PANEL:")
    log(f"  {'Gene':<12} {'Mature':>8} {'TNBC':>8} {'%chg':>8}  "
        f"{'p':>14}  r_depth")
    log(f"  {'-'*60}")
    for gene in EPIGENETIC:
        nm, tm, chg, p = test(gene, mature, tnbc)
        if nm is None:
            continue
        rv, _ = dr(gene)
        log(f"  {gene:<12} {nm:>8.4f} {tm:>8.4f} "
            f"{chg:>+7.1f}%  {fmt_p(p):>14}  {rv:>+.3f}")

    log(f"\n  P9 — HETEROGENEITY:")
    for gene in ["AR", "VIM", "CDH1", "ZEB1", "ZEB2", "SNAI1",
                 "CLDN3", "CLDN4"]:
        nm, tm, chg, p = test(gene, mature, tnbc)
        if nm is None:
            continue
        rv, _ = dr(gene)
        log(f"  {gene:<12} {chg:>+7.1f}% vs mature  "
            f"r_depth={rv:>+.3f}  {fmt_p(p)}")

    log(f"\n  CONTROLS (should be flat):")
    for gene in CONTROLS:
        nm, tm, chg, p = test(gene, mature, tnbc)
        if nm is None:
            log(f"  {gene:<12} not in data")
            continue
        flag = "✓ FLAT" if abs(chg or 0) < 30 else "⚠ ELEVATED"
        log(f"  {gene:<12} {chg:>+7.1f}%  {fmt_p(p):>14}  {flag}")

    if len(prog) > 5:
        log(f"\n  VS LUMINAL PROGENITOR:")
        log(f"  {'Gene':<12} {'Prog':>8} {'TNBC':>8} {'%chg':>8}  p")
        log(f"  {'-'*50}")
        for gene in (SWITCH_GENES + FA_MARKERS)[:10]:
            nm, tm, chg, p = test(gene, prog, tnbc)
            if nm is None:
                continue
            log(f"  {gene:<12} {nm:>8.4f} {tm:>8.4f} "
                f"{chg:>+7.1f}%  {fmt_p(p)}")

# ============================================================
# SECTION 5: NOVEL SIGNALS
# ============================================================

def novel_signals(rdf, corrs):
    log("")
    log("=" * 65)
    log("SECTION 5: NOVEL SIGNALS")
    log("Signals not in prediction panel — named without literature")
    log("=" * 65)

    panel = set(ALL_GENES)

    for label, direction in [
        ("Novel suppressed", "DOWN"),
        ("Novel elevated",   "UP"),
    ]:
        subset = rdf[
            (rdf["direction"] == direction) &
            (~rdf["gene"].isin(panel))
        ].head(10)
        log(f"\n  {label} (not in prediction panel):")
        for _, row in subset.iterrows():
            log(f"    {row['gene']:<20} {row['change_pct']:>+8.1f}%  "
                f"p={row['p_value']:.2e}")

    if corrs:
        novel_c = [
            (g, r, p) for g, r, p in corrs
            if g not in panel and abs(r) > 0.15
        ]
        if novel_c:
            log(f"\n  Novel depth correlates (|r|>0.15):")
            for g, r, p in novel_c[:10]:
                log(f"    {g:<20} r={r:>+.4f}  p={p:.2e}")

# ============================================================
# P6: pCR DEPTH TEST
# ============================================================

def pcr_depth_test(expr_bulk, meta_tnbc_bulk):
    log("")
    log("=" * 65)
    log("P6: pCR DEPTH TEST (GSE25066)")
    log("Predicted: r(depth_proxy, pCR_binary) < 0")
    log("=" * 65)

    if expr_bulk is None or meta_tnbc_bulk is None:
        log("  Bulk data unavailable. P6 skipped.")
        return

    if "pCR" not in meta_tnbc_bulk.columns:
        log("  pCR column not found. P6 skipped.")
        return

    pcr_known = meta_tnbc_bulk["pCR"].notna()
    if pcr_known.sum() < 10:
        log(f"  Only {pcr_known.sum()} pCR-annotated samples. P6 skipped.")
        return

    bulk_samples = [
        s for s in meta_tnbc_bulk.index
        if s in expr_bulk.columns
    ]
    log(f"  Matched bulk samples: {len(bulk_samples)}")

    if len(bulk_samples) < 10:
        log("  Cannot match samples to bulk matrix.")
        log(f"  Meta index sample: {list(meta_tnbc_bulk.index[:3])}")
        log(f"  Bulk cols sample:  {list(expr_bulk.columns[:3])}")
        return

    bulk_sub  = expr_bulk[bulk_samples].T
    meta_sub  = meta_tnbc_bulk.loc[bulk_samples].copy()

    def norm01_s(s):
        mn, mx = s.min(), s.max()
        return (s - mn) / (mx - mn) if mx > mn else pd.Series(0.5, index=s.index)

    basal_p   = [c for c in bulk_sub.columns
                 if any(g.lower() in str(c).lower()
                        for g in ["KRT5", "KRT14", "SOX10"])]
    luminal_p = [c for c in bulk_sub.columns
                 if any(g.lower() in str(c).lower()
                        for g in ["ESR1", "FOXA1", "GATA3"])]
    log(f"  Basal probes: {len(basal_p)}  Luminal probes: {len(luminal_p)}")

    parts = []
    if basal_p:
        parts.append(norm01_s(bulk_sub[basal_p].mean(axis=1)))
    if luminal_p:
        parts.append(1 - norm01_s(bulk_sub[luminal_p].mean(axis=1)))
    if not parts:
        log("  No target probes found. Using top-variance proxy.")
        top_var = bulk_sub.var(axis=0).nlargest(50).index
        parts.append(norm01_s(bulk_sub[top_var].mean(axis=1)))

    meta_sub["depth_proxy"] = (sum(parts) / len(parts)).values

    pcr_col   = meta_sub["pCR"].dropna()
    depth_col = meta_sub.loc[pcr_col.index, "depth_proxy"]

    try:
        rv, pv = stats.pointbiserialr(pcr_col.values, depth_col.values)
        log(f"\n  r(depth_proxy, pCR) = {rv:>+.4f}   p = {pv:.4f}")
        if rv < 0 and pv < 0.05:
            log("  P6: ✓ CONFIRMED — deeper TNBC → lower pCR probability")
        elif rv < 0:
            log(f"  P6: ✓ DIRECTIONALLY CONFIRMED (ns, p={pv:.3f})")
        else:
            log(f"  P6: ✗ NOT CONFIRMED — consider alternative hypothesis")
            log("  (High proliferation in deep tumors may increase chemo sensitivity)")
    except Exception as e:
        log(f"  P6 error: {e}")

    pcr1 = meta_sub.loc[meta_sub["pCR"] == 1, "depth_proxy"]
    pcr0 = meta_sub.loc[meta_sub["pCR"] == 0, "depth_proxy"]
    log(f"  pCR=1 mean depth: {pcr1.mean():.4f}  n={len(pcr1)}")
    log(f"  pCR=0 mean depth: {pcr0.mean():.4f}  n={len(pcr0)}")

# ============================================================
# FIGURE
# ============================================================

def generate_figure(pops, rdf, depth, corrs, pca_result):
    log("")
    log("--- Generating figure ---")

    tnbc   = pops.get("tnbc",   pd.DataFrame())
    mature = pops.get("mature", pd.DataFrame())
    prog   = pops.get("prog",   pd.DataFrame())
    pcs, labels = pca_result if pca_result else (None, None)

    clr = {
        "tnbc":   "#c0392b",
        "mature": "#2980b9",
        "prog":   "#27ae60",
        "myo":    "#8e44ad",
        "up":     "#c0392b",
        "down":   "#2980b9",
    }

    fig = plt.figure(figsize=(24, 18))
    fig.suptitle(
        f"TNBC — Script 1 | OrganismCore BRCA-S4a/b | 2026-03-04\n"
        f"Composite Type 1→2 | EZH2 Convergence Node | GSE176078\n"
        f"TNBC={len(tnbc)}  Mature={len(mature)}  Progenitor={len(prog)}",
        fontsize=10, fontweight="bold", y=1.005
    )
    gs = gridspec.GridSpec(3, 3, figure=fig,
                           hspace=0.52, wspace=0.42)

    # A: Switch genes
    ax_a = fig.add_subplot(gs[0, 0])
    sw = [g for g in SWITCH_GENES if g in tnbc.columns]
    if sw and len(mature) > 0:
        x = np.arange(len(sw)); w = 0.26
        mv = [mature[g].mean() if g in mature.columns else 0 for g in sw]
        pv = [prog[g].mean()   if g in prog.columns   else 0 for g in sw]
        tv = [tnbc[g].mean()                                  for g in sw]
        ax_a.bar(x-w, mv, w, color=clr["mature"], label="Mature",   alpha=0.85)
        ax_a.bar(x,   pv, w, color=clr["prog"],   label="Progenitor", alpha=0.85)
        ax_a.bar(x+w, tv, w, color=clr["tnbc"],   label="TNBC",     alpha=0.85)
        ax_a.set_xticks(x)
        ax_a.set_xticklabels(sw, rotation=45, ha="right", fontsize=6)
        ax_a.set_title("A — Switch Genes (P1)\nPredicted DOWN", fontsize=8, fontweight="bold")
        ax_a.legend(fontsize=6)
        ax_a.set_ylabel("log1p expression", fontsize=7)

    # B: FA markers
    ax_b = fig.add_subplot(gs[0, 1])
    fa = [g for g in FA_MARKERS if g in tnbc.columns]
    if fa and len(mature) > 0:
        x = np.arange(len(fa)); w = 0.26
        mv = [mature[g].mean() if g in mature.columns else 0 for g in fa]
        pv = [prog[g].mean()   if g in prog.columns   else 0 for g in fa]
        tv = [tnbc[g].mean()                                  for g in fa]
        ax_b.bar(x-w, mv, w, color=clr["mature"], label="Mature",   alpha=0.85)
        ax_b.bar(x,   pv, w, color=clr["prog"],   label="Progenitor", alpha=0.85)
        ax_b.bar(x+w, tv, w, color=clr["tnbc"],   label="TNBC",     alpha=0.85)
        ax_b.set_xticks(x)
        ax_b.set_xticklabels(fa, rotation=45, ha="right", fontsize=6)
        ax_b.set_title("B — FA Markers (P2)\nPredicted UP", fontsize=8, fontweight="bold")
        ax_b.legend(fontsize=6)
        ax_b.set_ylabel("log1p expression", fontsize=7)

    # C: Waterfall
    ax_c = fig.add_subplot(gs[0, 2])
    if rdf is not None and len(rdf) > 0:
        sig  = rdf[rdf["p_value"] < 0.05].sort_values("change_pct")
        show = pd.concat([sig.head(15), sig.tail(15)]).drop_duplicates()
        bc   = [clr["down"] if v < 0 else clr["up"] for v in show["change_pct"]]
        ax_c.barh(show["gene"], show["change_pct"], color=bc)
        ax_c.axvline(0, color="black", lw=0.8)
        ax_c.set_title(f"C — Waterfall\nTop {len(show)} movers", fontsize=8, fontweight="bold")
        ax_c.set_xlabel("% change vs mature luminal", fontsize=7)
        ax_c.tick_params(axis="y", labelsize=5)

    # D: Depth distribution
    ax_d = fig.add_subplot(gs[1, 0])
    if depth is not None:
        ax_d.hist(depth.values, bins=40, color=clr["tnbc"], alpha=0.75, edgecolor="white")
        ax_d.axvline(depth.mean(), color="black", lw=1.5, ls="--",
                     label=f"mean={depth.mean():.3f}")
        ax_d.set_title("D — Depth Distribution", fontsize=8, fontweight="bold")
        ax_d.set_xlabel("Depth Score", fontsize=7)
        ax_d.set_ylabel("n cells", fontsize=7)
        ax_d.legend(fontsize=6)

    # E: EZH2 vs depth
    ax_e = fig.add_subplot(gs[1, 1])
    if depth is not None and "EZH2" in tnbc.columns:
        ax_e.scatter(tnbc["EZH2"].values, depth.values,
                     c=clr["tnbc"], alpha=0.1, s=3)
        try:
            rv, pv = stats.pearsonr(tnbc["EZH2"].values, depth.values)
            check  = "✓ r>0.30" if rv > 0.30 else f"? r={rv:+.3f}"
            ax_e.set_title(f"E — EZH2 vs Depth (P3)\nr={rv:+.3f}  {check}",
                           fontsize=8, fontweight="bold")
        except Exception:
            ax_e.set_title("E — EZH2 vs Depth", fontsize=8)
        ax_e.set_xlabel("EZH2", fontsize=7)
        ax_e.set_ylabel("Depth", fontsize=7)
        if "EZH2" in mature.columns:
            ax_e.axvline(mature["EZH2"].mean(), color=clr["mature"],
                         lw=1.5, ls="--", label="Mature mean")
            ax_e.legend(fontsize=6)

    # F: Depth correlations
    ax_f = fig.add_subplot(gs[1, 2])
    if corrs:
        top20 = corrs[:20]
        gc = [g for g, r, p in top20]
        vc = [r for g, r, p in top20]
        cc = [clr["up"] if v > 0 else clr["down"] for v in vc]
        ax_f.barh(gc, vc, color=cc)
        ax_f.axvline(0, color="black", lw=0.8)
        ax_f.set_title("F — Depth Correlations\nTop 20", fontsize=8, fontweight="bold")
        ax_f.set_xlabel("r", fontsize=7)
        ax_f.tick_params(axis="y", labelsize=6)

    # G: PCA
    ax_g = fig.add_subplot(gs[2, 0])
    if pcs is not None and labels is not None:
        for lbl, ck in [("TNBC","tnbc"),("Mature","mature"),
                         ("Progenitor","prog"),("Myoepithelial","myo")]:
            mask = labels == lbl
            if mask.sum() > 0:
                ax_g.scatter(pcs[mask,0], pcs[mask,1],
                             c=clr[ck], alpha=0.15, s=3, label=lbl)
        ax_g.set_title("G — PCA", fontsize=8, fontweight="bold")
        ax_g.set_xlabel("PC1", fontsize=7)
        ax_g.set_ylabel("PC2", fontsize=7)
        ax_g.legend(fontsize=6, markerscale=3)

    # H: Composite type test
    ax_h = fig.add_subplot(gs[2, 1])
    if "BRCA1" in tnbc.columns and "ESR1" in tnbc.columns:
        ax_h.scatter(tnbc["BRCA1"].values, tnbc["ESR1"].values,
                     c=clr["tnbc"], alpha=0.1, s=3)
        try:
            rv, pv = stats.pearsonr(tnbc["BRCA1"].values, tnbc["ESR1"].values)
            check = "✓ r>0.15" if rv > 0.15 else "✗ r<0.15"
            ax_h.set_title(f"H — Composite Type (P4)\nr(BRCA1,ESR1)={rv:+.3f}  {check}",
                           fontsize=8, fontweight="bold")
        except Exception:
            ax_h.set_title("H — BRCA1 vs ESR1", fontsize=8)
        ax_h.set_xlabel("BRCA1", fontsize=7)
        ax_h.set_ylabel("ESR1", fontsize=7)

    # I: Summary
    ax_i = fig.add_subplot(gs[2, 2])
    ax_i.axis("off")
    ax_i.text(0.03, 0.97,
        f"I — SCRIPT 1 SUMMARY\n{'─'*30}\n"
        f"GSE176078 | Wu et al. 2021\n"
        f"TNBC={len(tnbc)}  Mature={len(mature)}\n\n"
        f"TYPE: Composite T1→T2\n"
        f"  T1: BRCA1 loss\n"
        f"  T2: EZH2 basal attractor\n\n"
        f"P1 Switch: panel A\n"
        f"P2 FA:     panel B\n"
        f"P3 EZH2:   panel E\n"
        f"P4 Comp:   panel H\n\n"
        f"DRUGS:\n"
        f"1. EZH2i (tazemetostat)\n"
        f"2. PARPi (olaparib)\n"
        f"3. EZH2+PARP (novel)\n"
        f"4. PD-L1 depth-stratified\n\n"
        f"OrganismCore | BRCA-S4a/b\n"
        f"2026-03-04",
        transform=ax_i.transAxes, fontsize=7.5,
        verticalalignment="top", fontfamily="monospace",
        bbox=dict(boxstyle="round", facecolor="#f8f8f8",
                  edgecolor="#cccccc"))

    plt.savefig(FIG_FILE, dpi=150, bbox_inches="tight")
    log(f"  Figure saved: {FIG_FILE}")
    plt.close()

# ============================================================
# MAIN
# ============================================================

def main():
    log("=" * 65)
    log("BRCA TNBC — SCRIPT 1")
    log("OrganismCore — BRCA-S4a/b | 2026-03-04")
    log("=" * 65)
    log("")
    log("CRITICAL: Uses GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz")
    log("          NOT GSE176078_RAW.tar (that is per-sample files)")
    log("")

    sc_dir, bulk_path = acquire_data()
    if sc_dir is None:
        log("FATAL: Data acquisition failed.")
        write_log()
        return

    meta = load_metadata(sc_dir)
    if meta is None:
        log("FATAL: Metadata load failed.")
        write_log()
        return

    expr = load_sc_expression(sc_dir, meta)
    if expr is None:
        log("FATAL: Expression load failed.")
        write_log()
        return

    pops = classify_populations(expr, meta)
    if pops is None:
        log("FATAL: Population classification failed.")
        write_log()
        return

    expr_bulk, meta_bulk = load_bulk_pcr(bulk_path)

    rdf        = top_movers_unfiltered(pops)        # Section 1
    pca_result = pca_geometry(pops)                  # Section 2
    depth, corrs = compute_depth(pops)               # Section 3
    prediction_check(pops, depth, corrs)             # Section 4
    if rdf is not None and len(rdf) > 0:
        novel_signals(rdf, corrs)                    # Section 5
    pcr_depth_test(expr_bulk, meta_bulk)             # P6

    generate_figure(pops, rdf, depth, corrs, pca_result)

    write_log()

    log("")
    log("=" * 65)
    log("SCRIPT 1 COMPLETE")
    log(f"  Log    : {LOG_FILE}")
    log(f"  Figure : {FIG_FILE}")
    log(f"  CSV    : {CSV_FILE}")
    log("=" * 65)


if __name__ == "__main__":
    main()
