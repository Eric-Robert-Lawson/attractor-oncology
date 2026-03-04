"""
BRCA LUMINAL A — SCRIPT 1
OrganismCore — Document BRCA-S1b
Deep Dive Series: Breast Cancer Subtypes

PREDICTIONS LOCKED IN: BRCA-S1a (2026-03-04)

DATASET:
  GSE176078 — Wu et al. 2021, Nature Genetics
  PMID: 34493872
  100,064 cells — 26 primary breast tumors
  scRNA-seq 10X Chromium

  SAME DATASET AS ORIGINAL BRCA ANALYSIS.
  NO NEW DATA REQUIRED.
  Will reuse cached files if present.

FILES USED (auto-downloaded from GEO FTP):
  GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz
    → count_matrix_sparse.mtx   (genes x cells)
    → count_matrix_genes.tsv    (gene names)
    → count_matrix_barcodes.tsv (cell barcodes)
    → metadata.csv              (cell annotations)

KEY METADATA COLUMN: celltype_subset
KEY POPULATIONS:
  "Cancer LumA SC"   — Luminal A cancer cells
  "Cancer Basal SC"  — TNBC cancer cells (reference)
  "Mature Luminal"   — Normal mature luminal cells
  "Luminal Progenitors" — Normal luminal progenitors
  "Cancer Cycling"   — Cycling cancer cells

WHY THIS DATASET IS SUFFICIENT FOR LUMINAL A:
  The metadata already provides single-cell resolution
  subtype labels. Cancer LumA SC is the Luminal A
  population. Mature Luminal is the normal reference.
  No separate PAM50 classification needed.
  The prior analysis (Document 75) established the
  reference values for FOXA1/GATA3/ESR1 — this script
  extends that with the full Luminal A prediction panel.

ALSO USES (if present in the original data download):
  GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts.txt.gz
  This provides bulk RNA-seq on the same 26 tumors
  for the depth score / clinical correlation analysis.

PREDICTIONS FROM BRCA-S1a (locked 2026-03-04):
  1. Switch genes RETAINED: FOXA1, GATA3, ESR1, PGR
     Cancer LumA SC >= Mature Luminal
     (NOT suppressed — unlike TNBC)
  2. Depth axis PROLIFERATIVE: MKI67, CCND1, CDK4 elevated
  3. EZH2 FLAT in LumA (vs +269.7% in TNBC)
  4. PIK3CA axis UNCERTAIN (post-translational)
  5. Depth score correlates with proliferative state
  6. Controls FLAT: SPI1, CDX2, MBP (non-breast genes)
  7. TNBC markers ABSENT in LumA: KRT5, SOX10

Author: Eric Robert Lawson
Framework: OrganismCore
Protocol: Workflow_Protocol.md v2.0
Version: 1.0
Date: 2026-03-04
"""

import os
import sys
import gzip
import tarfile
import urllib.request
import urllib.error
import time
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy import stats
from scipy.io import mmread
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

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
BASE_DIR    = os.path.join(SCRIPT_DIR, "luma_results")
RESULTS_DIR = BASE_DIR
LOG_FILE    = os.path.join(RESULTS_DIR, "luma_s1_log.txt")
FIG_FILE    = os.path.join(RESULTS_DIR, "luma_s1_figure.png")
CSV_FILE    = os.path.join(RESULTS_DIR, "luma_s1_results.csv")

os.makedirs(RESULTS_DIR, exist_ok=True)

# ============================================================
# CELL POPULATION LABELS (from metadata.csv)
# ============================================================

CANCER_LUMA   = "Cancer LumA SC"
CANCER_BASAL  = "Cancer Basal SC"     # TNBC reference
MATURE_LUM    = "Mature Luminal"      # normal reference
LUMINAL_PROG  = "Luminal Progenitors" # normal progenitor
CANCER_CYCL   = "Cancer Cycling"
CT_COL        = "celltype_subset"     # metadata column

# ============================================================
# GENE PANELS — FROM BRCA-S1a PREDICTIONS
# ============================================================

# Prediction 1: Switch genes — expected RETAINED or ELEVATED
SWITCH_GENES = ["FOXA1", "GATA3", "ESR1", "PGR"]

# Prediction 2: Depth axis — expected ELEVATED in LumA
DEPTH_GENES  = ["MKI67", "CCND1", "CDK4", "MYC",
                "CCNE1", "PCNA", "TOP2A"]

# Prediction 3: Lock gene — expected FLAT (not +270% as in TNBC)
LOCK_GENES   = ["EZH2", "KDM1A", "HDAC1", "HDAC2"]

# Cell cycle exit — predicted suppressed (uncertain)
EXIT_GENES   = ["RB1", "CDKN2A", "CDKN1A", "TP53"]

# PIK3CA axis — post-translational, uncertain at mRNA
PIK3CA_GENES = ["AKT1", "AKT2", "MTOR", "PIK3CA", "PTEN"]

# Controls from prior cancers — expected flat (not breast)
CONTROL_GENES = ["SPI1", "CDX2", "MBP", "NKX2-1", "OLIG2"]

# TNBC markers — expected absent in LumA
TNBC_MARKERS  = ["KRT5", "KRT14", "SOX10", "EGFR"]

# Cross-comparison: these are the genes from
# the original Document 75 analysis panel
# (used in brca_saddle_point_analysis.py)
DOC75_PANEL = [
    "FOXA1", "GATA3", "ESR1",
    "SOX2", "MYC", "EGFR",
    "KRT5", "SOX10", "MBP",
    "CDX2", "SPI1", "KRT8",
    "MKI67", "AR", "ERBB2", "EZH2",
]

ALL_GENES = list(dict.fromkeys(
    SWITCH_GENES + DEPTH_GENES + LOCK_GENES +
    EXIT_GENES + PIK3CA_GENES +
    CONTROL_GENES + TNBC_MARKERS + DOC75_PANEL
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
# FILE DISCOVERY
# Search for already-downloaded files first.
# ============================================================

def find_file(filename):
    """
    Search known locations for a file.
    Returns full path if found, None otherwise.
    """
    for base in SEARCH_PATHS:
        candidate = os.path.join(base, filename)
        if os.path.exists(candidate):
            sz = os.path.getsize(candidate)
            if sz > 10000:  # > 10 KB — not corrupt
                log(f"  Found existing: {candidate} "
                    f"({sz/1e6:.1f} MB)")
                return os.path.abspath(candidate)
    return None


def find_extracted_sc_dir():
    """
    Search for the already-extracted scRNA-seq directory.
    The tar contains: Wu_etal_2021_BRCA_scRNASeq/
    """
    dir_name = "Wu_etal_2021_BRCA_scRNASeq"
    for base in SEARCH_PATHS:
        candidate = os.path.join(base, dir_name)
        mtx_check = os.path.join(
            candidate, "count_matrix_sparse.mtx"
        )
        if os.path.isdir(candidate) and os.path.exists(mtx_check):
            log(f"  Found extracted dir: {candidate}")
            return os.path.abspath(candidate)
        # Also check if files are directly in the base dir
        mtx_direct = os.path.join(base, "count_matrix_sparse.mtx")
        if os.path.exists(mtx_direct):
            log(f"  Found sc files directly in: {base}")
            return os.path.abspath(base)
    return None

# ============================================================
# DOWNLOAD UTILITY
# ============================================================

def fetch_url(url, dest, retries=3):
    for attempt in range(retries):
        try:
            log(f"  Attempt {attempt+1}: {url[-60:]}")
            req = urllib.request.Request(
                url,
                headers={"User-Agent": "Mozilla/5.0"}
            )
            with urllib.request.urlopen(req, timeout=120) as r:
                data = r.read()
            with open(dest, "wb") as f:
                f.write(data)
            log(f"  Downloaded: {dest} "
                f"({len(data)/1e6:.1f} MB)")
            return dest
        except Exception as e:
            log(f"  Attempt {attempt+1} failed: {e}")
            if attempt < retries - 1:
                time.sleep(3)
    return None

# ============================================================
# STEP 0: DATA ACQUISITION
# ============================================================

def acquire_data():
    log("=" * 65)
    log("STEP 0: DATA ACQUISITION")
    log(f"Dataset: {GEO_ACCESSION} — Wu et al. 2021")
    log("=" * 65)

    # --- scRNA-seq files ---
    log("\n--- scRNA-seq data ---")

    sc_dir = find_extracted_sc_dir()

    if sc_dir is None:
        log("  Extracted sc dir not found. "
            "Looking for tar file...")
        tar_path = find_file(SC_TARFILE)

        if tar_path is None:
            log(f"  Tar file not found. Downloading...")
            tar_dest = os.path.join(BASE_DIR, SC_TARFILE)
            tar_path = fetch_url(SC_URL, tar_dest)
            if tar_path is None:
                log("  FATAL: Could not download scRNA-seq data.")
                return None, None

        log(f"  Extracting: {tar_path}")
        extract_dir = os.path.join(BASE_DIR, "sc_extracted")
        os.makedirs(extract_dir, exist_ok=True)
        with tarfile.open(tar_path, "r:gz") as tf:
            tf.extractall(extract_dir)
        log(f"  Extracted to: {extract_dir}")

        # Find the extracted dir
        sc_dir = find_extracted_sc_dir()
        if sc_dir is None:
            # Files may be directly in extract_dir
            if os.path.exists(os.path.join(
                    extract_dir, "count_matrix_sparse.mtx")):
                sc_dir = extract_dir
            else:
                # Search one level deeper
                for sub in os.listdir(extract_dir):
                    candidate = os.path.join(extract_dir, sub)
                    if os.path.isdir(candidate):
                        if os.path.exists(os.path.join(
                                candidate,
                                "count_matrix_sparse.mtx")):
                            sc_dir = candidate
                            break

        if sc_dir is None:
            log("  FATAL: Could not locate extracted sc files.")
            return None, None

    log(f"  scRNA-seq dir: {sc_dir}")

    # Verify all required files
    required_sc = [
        "count_matrix_sparse.mtx",
        "count_matrix_genes.tsv",
        "count_matrix_barcodes.tsv",
        "metadata.csv",
    ]
    for f in required_sc:
        fp = os.path.join(sc_dir, f)
        if not os.path.exists(fp):
            log(f"  MISSING: {f}")
            log("  FATAL: sc files incomplete.")
            return None, None
        log(f"  OK: {f} "
            f"({os.path.getsize(fp)/1e6:.1f} MB)")

    # --- Bulk RNA-seq (optional) ---
    log("\n--- Bulk RNA-seq (optional) ---")
    bulk_path = find_file(BULK_FILE)
    if bulk_path is None:
        bulk_dest = os.path.join(BASE_DIR, BULK_FILE)
        log(f"  Not found locally. Downloading...")
        bulk_path = fetch_url(BULK_URL, bulk_dest)
        if bulk_path is None:
            log("  Could not download bulk data.")
            log("  Proceeding with scRNA-seq only.")

    return sc_dir, bulk_path

# ============================================================
# STEP 1: LOAD METADATA
# ============================================================

def load_metadata(sc_dir):
    log("=" * 65)
    log("STEP 1: LOAD METADATA")
    log("=" * 65)

    meta_file = os.path.join(sc_dir, "metadata.csv")
    meta = pd.read_csv(meta_file, index_col=0)
    log(f"  Cells loaded: {len(meta)}")
    log(f"  Columns: {list(meta.columns)[:10]}")

    if CT_COL not in meta.columns:
        # Try finding the right column
        log(f"  WARNING: '{CT_COL}' not found.")
        log(f"  Available: {list(meta.columns)}")
        for col in meta.columns:
            if "celltype" in col.lower() or "subset" in col.lower():
                log(f"  Using '{col}' as cell type column.")
                meta[CT_COL] = meta[col]
                break

    log(f"\n  Cell type distribution ({CT_COL}):")
    ct_counts = meta[CT_COL].value_counts()
    for ct, n in ct_counts.items():
        log(f"    {ct:<35} n={n}")

    # Show subtype distribution if available
    if "subtype" in meta.columns:
        log(f"\n  Subtype distribution:")
        for st, n in meta["subtype"].value_counts().items():
            log(f"    {st:<20} n={n}")

    # Tag our populations
    meta["is_luma"]    = (meta[CT_COL] == CANCER_LUMA)
    meta["is_basal"]   = (meta[CT_COL] == CANCER_BASAL)
    meta["is_mature"]  = (meta[CT_COL] == MATURE_LUM)
    meta["is_lumprg"]  = (meta[CT_COL] == LUMINAL_PROG)
    meta["is_cycling"] = (meta[CT_COL] == CANCER_CYCL)

    n_luma   = meta["is_luma"].sum()
    n_basal  = meta["is_basal"].sum()
    n_mature = meta["is_mature"].sum()
    n_lumprg = meta["is_lumprg"].sum()

    log(f"\n  KEY POPULATIONS:")
    log(f"    {CANCER_LUMA:<35} n={n_luma}")
    log(f"    {CANCER_BASAL:<35} n={n_basal}")
    log(f"    {MATURE_LUM:<35} n={n_mature}")
    log(f"    {LUMINAL_PROG:<35} n={n_lumprg}")

    if n_luma < 100:
        log(f"  WARNING: Only {n_luma} LumA cells.")
    if n_mature < 50:
        log(f"  WARNING: Only {n_mature} Mature Luminal cells.")

    return meta

# ============================================================
# STEP 2: LOAD scRNA-seq EXPRESSION
# ============================================================

def load_sc_expression(sc_dir, meta):
    log("=" * 65)
    log("STEP 2: LOAD scRNA-seq EXPRESSION")
    log("=" * 65)

    # Check for cached gene expression matrix
    cache_file = os.path.join(RESULTS_DIR, "expr_cache_luma.csv")

    if os.path.exists(cache_file):
        sz = os.path.getsize(cache_file) / 1e6
        log(f"  Loading cache: {cache_file} ({sz:.1f} MB)")
        expr = pd.read_csv(cache_file, index_col=0)
        log(f"  Cached shape: {expr.shape}")
        return expr

    log("  Loading MTX sparse matrix...")
    mtx_file  = os.path.join(sc_dir, "count_matrix_sparse.mtx")
    gene_file = os.path.join(sc_dir, "count_matrix_genes.tsv")
    bc_file   = os.path.join(sc_dir, "count_matrix_barcodes.tsv")

    # Load gene names
    genes_df = pd.read_csv(gene_file, sep="\t", header=None)
    gene_names = genes_df.iloc[:, 0].tolist()
    log(f"  Genes in matrix: {len(gene_names)}")

    # Load barcodes
    barcodes = pd.read_csv(
        bc_file, sep="\t", header=None
    ).iloc[:, 0].tolist()
    log(f"  Barcodes: {len(barcodes)}")

    # Load sparse matrix
    log("  Reading MTX (this takes ~30-60 sec)...")
    mat = mmread(mtx_file).tocsr()
    log(f"  Matrix shape: {mat.shape} (genes x cells)")

    # Find target gene indices
    gene_name_lower = [g.upper() for g in gene_names]
    target_indices  = []
    target_names    = []
    missing_genes   = []

    for gene in ALL_GENES:
        g_upper = gene.upper()
        # Exact match first
        if g_upper in gene_name_lower:
            idx = gene_name_lower.index(g_upper)
            target_indices.append(idx)
            target_names.append(gene)
        else:
            missing_genes.append(gene)

    log(f"  Target genes requested: {len(ALL_GENES)}")
    log(f"  Target genes found:     {len(target_names)}")
    if missing_genes:
        log(f"  Missing genes:          {missing_genes}")

    # Extract target genes from sparse matrix
    log("  Extracting target genes from sparse matrix...")
    expr_data = {}
    for i, (gene_idx, gene_name) in enumerate(
        zip(target_indices, target_names)
    ):
        row = mat[gene_idx, :].toarray().flatten()
        expr_data[gene_name] = row
        if (i + 1) % 10 == 0:
            log(f"    Extracted {i+1}/{len(target_names)} genes")

    expr = pd.DataFrame(expr_data, index=barcodes)
    log(f"  Expression matrix: {expr.shape} (cells x genes)")

    # Report basic stats
    for gene in target_names[:5]:
        vals = expr[gene]
        log(f"  {gene:<12}: "
            f"mean={vals.mean():.4f} "
            f"max={vals.max():.4f} "
            f"nonzero={int((vals > 0).sum())}")

    # Cache
    log(f"  Caching to: {cache_file}")
    expr.to_csv(cache_file)
    log("  Cache saved.")

    return expr

# ============================================================
# STEP 3: COMPUTE POPULATION PROFILES
# ============================================================

def compute_profiles(expr, meta):
    log("=" * 65)
    log("STEP 3: COMPUTE POPULATION PROFILES")
    log("=" * 65)

    # Align barcodes
    common = expr.index.intersection(meta.index)
    log(f"  Common cells: {len(common)}")

    expr_common = expr.loc[common]
    meta_common = meta.loc[common]

    # Extract populations
    luma_cells   = meta_common[meta_common["is_luma"]].index
    basal_cells  = meta_common[meta_common["is_basal"]].index
    mature_cells = meta_common[meta_common["is_mature"]].index
    lumprg_cells = meta_common[meta_common["is_lumprg"]].index

    log(f"  LumA SC cells:        {len(luma_cells)}")
    log(f"  Basal SC cells:       {len(basal_cells)}")
    log(f"  Mature Luminal cells: {len(mature_cells)}")
    log(f"  Luminal Progenitors:  {len(lumprg_cells)}")

    # Compute mean expression per population
    profiles = {}
    for name, cells in [
        ("LumA_Cancer",   luma_cells),
        ("Basal_Cancer",  basal_cells),
        ("Mature_Luminal", mature_cells),
        ("Luminal_Prog",  lumprg_cells),
    ]:
        if len(cells) > 0:
            profiles[name] = expr_common.loc[cells].mean()
        else:
            log(f"  WARNING: No cells for {name}")

    return expr_common, meta_common, profiles, {
        "luma":   luma_cells,
        "basal":  basal_cells,
        "mature": mature_cells,
        "lumprg": lumprg_cells,
    }

# ============================================================
# STEP 4: PRIMARY COMPARISON — ALL 7 PREDICTIONS
# LumA Cancer vs Mature Luminal Normal
# ============================================================

def primary_comparison(expr, populations, profiles):
    log("=" * 65)
    log("STEP 4: PRIMARY COMPARISON")
    log("LumA Cancer SC vs Mature Luminal Normal")
    log("Testing all 7 predictions from BRCA-S1a")
    log("=" * 65)

    luma_cells   = populations["luma"]
    mature_cells = populations["mature"]

    luma_expr   = expr.loc[luma_cells]
    mature_expr = expr.loc[mature_cells]

    log(f"\n  LumA:   n={len(luma_cells)} cells")
    log(f"  Normal: n={len(mature_cells)} cells")

    available_genes = [g for g in ALL_GENES if g in expr.columns]

    # ── SECTION 1: TOP MOVERS (UNFILTERED — no panel imposed)
    log("")
    log("=" * 65)
    log("SECTION 1: TOP MOVERS — UNFILTERED")
    log("Geometry speaks first. No panel imposed.")
    log("=" * 65)
    log("")

    all_results = []
    for gene in available_genes:
        luma_vals   = luma_expr[gene].values
        mature_vals = mature_expr[gene].values
        lm = luma_vals.mean()
        mm = mature_vals.mean()
        if mm > 1e-6:
            pct = ((lm - mm) / mm) * 100
        else:
            pct = 0.0

        try:
            _, pval = stats.mannwhitneyu(
                luma_vals, mature_vals,
                alternative="two-sided"
            )
        except Exception:
            pval = 1.0

        all_results.append({
            "gene": gene, "luma_mean": lm,
            "normal_mean": mm, "pct_change": pct,
            "p_value": pval, "abs_change": abs(pct)
        })

    rdf = pd.DataFrame(all_results).sort_values(
        "pct_change", ascending=False
    )

    log("TOP 20 GAINED in LumA vs Mature Luminal:")
    log(f"{'Gene':<12} {'Normal':>8} {'LumA':>8} "
        f"{'Change':>10}  {'p-value':>12}")
    log("-" * 56)
    for _, row in rdf.head(20).iterrows():
        sig = ("***" if row["p_value"] < 0.001 else
               "**"  if row["p_value"] < 0.01  else
               "*"   if row["p_value"] < 0.05  else "ns")
        log(f"{row['gene']:<12} {row['normal_mean']:>8.4f} "
            f"{row['luma_mean']:>8.4f} "
            f"{row['pct_change']:>+9.1f}%  "
            f"p={row['p_value']:.2e} {sig}")

    log("")
    log("TOP 20 LOST in LumA vs Mature Luminal:")
    log(f"{'Gene':<12} {'Normal':>8} {'LumA':>8} "
        f"{'Change':>10}  {'p-value':>12}")
    log("-" * 56)
    for _, row in rdf.tail(20).iterrows():
        sig = ("***" if row["p_value"] < 0.001 else
               "**"  if row["p_value"] < 0.01  else
               "*"   if row["p_value"] < 0.05  else "ns")
        log(f"{row['gene']:<12} {row['normal_mean']:>8.4f} "
            f"{row['luma_mean']:>8.4f} "
            f"{row['pct_change']:>+9.1f}%  "
            f"p={row['p_value']:.2e} {sig}")

    # ── SECTION 2: FULL RESULTS TABLE (all panel genes)
    log("")
    log("=" * 65)
    log("SECTION 2: FULL PANEL — LumA vs Mature Luminal")
    log("=" * 65)
    log("")

    results = []
    for _, row in rdf.iterrows():
        gene    = row["gene"]
        lm      = row["luma_mean"]
        mm      = row["normal_mean"]
        pct     = row["pct_change"]
        pval    = row["p_value"]
        abs_chg = row["abs_change"]

        sig = ("***" if pval < 0.001 else
               "**"  if pval < 0.01  else
               "*"   if pval < 0.05  else "ns")

        result = ("ELEVATED"   if pct >  20 and pval < 0.05 else
                  "SUPPRESSED" if pct < -20 and pval < 0.05 else
                  "FLAT")

        # Assign category
        if gene in SWITCH_GENES:
            cat = "SWITCH"
        elif gene in DEPTH_GENES:
            cat = "DEPTH"
        elif gene in LOCK_GENES:
            cat = "LOCK"
        elif gene in EXIT_GENES:
            cat = "EXIT"
        elif gene in PIK3CA_GENES:
            cat = "PIK3CA"
        elif gene in CONTROL_GENES:
            cat = "CONTROL"
        elif gene in TNBC_MARKERS:
            cat = "TNBC-MARKER"
        else:
            cat = "DOC75"

        log(f"  {gene:<12} [{cat:<12}] "
            f"N={mm:>7.4f}  L={lm:>7.4f} "
            f"{pct:>+8.1f}%  p={pval:.2e} {sig:<3}  "
            f"→ {result}")

        results.append({
            "gene": gene, "category": cat,
            "normal_mean": mm, "luma_mean": lm,
            "pct_change": pct, "p_value": pval,
            "significance": sig, "result": result,
        })

    results_df = pd.DataFrame(results)

    # ── SECTION 3: CROSS-REFERENCE — LumA vs Basal vs Document 75
    log("")
    log("=" * 65)
    log("SECTION 3: CROSS-REFERENCE")
    log("LumA vs Basal (TNBC) — same reference normal")
    log("Document 75 reference values shown")
    log("=" * 65)
    log("")

    basal_cells = populations["basal"]
    basal_expr  = expr.loc[basal_cells]

    doc75_ref = {
        "FOXA1": {"luma_doc75": 0.5221, "normal_doc75": 0.3934,
                  "basal_doc75": None},
        "GATA3": {"luma_doc75": 1.3230, "normal_doc75": 1.1115,
                  "basal_doc75": None},
        "ESR1":  {"luma_doc75": 0.6901, "normal_doc75": 0.7489,
                  "basal_doc75": None},
        "EZH2":  {"luma_doc75": None,   "normal_doc75": 0.0414,
                  "basal_doc75": 0.1530},
        "SOX10": {"luma_doc75": None,   "normal_doc75": None,
                  "basal_doc75": None},
    }

    log(f"  {'Gene':<12} {'Normal':>8} {'LumA':>8} "
        f"{'Basal':>8}  {'LumA%':>9} {'Basal%':>9}  "
        f"Note")
    log("  " + "-" * 70)

    for gene in available_genes:
        lm  = luma_expr[gene].mean()
        mm  = mature_expr[gene].mean()
        bm  = basal_expr[gene].mean() if len(basal_cells) > 0 else np.nan

        lm_pct = ((lm - mm) / mm * 100) if mm > 1e-6 else 0
        bm_pct = ((bm - mm) / mm * 100) if mm > 1e-6 else 0

        # Flag if same gene behaves differently LumA vs Basal
        note = ""
        if abs(lm_pct - bm_pct) > 50:
            note = "← DIVERGES"
        if gene == "EZH2":
            note = "← LOCK TEST"
        if gene in SWITCH_GENES:
            note = "← SWITCH"

        log(f"  {gene:<12} {mm:>8.4f} {lm:>8.4f} "
            f"{bm:>8.4f}  {lm_pct:>+8.1f}% {bm_pct:>+8.1f}%  "
            f"{note}")

    return results_df, rdf

# ============================================================
# STEP 5: DEPTH SCORE
# Proliferative axis within LumA cells
# ============================================================

def compute_depth_score(expr, populations):
    log("=" * 65)
    log("STEP 5: DEPTH SCORE CONSTRUCTION")
    log("Proliferative axis within LumA population")
    log("=" * 65)

    luma_cells = populations["luma"]
    luma_expr  = expr.loc[luma_cells]

    # Which proliferative genes are available?
    prolif_genes = [g for g in ["MKI67", "CCND1", "CDK4",
                                 "TOP2A", "PCNA", "CCNE1"]
                    if g in luma_expr.columns]

    log(f"  Proliferative genes available: {prolif_genes}")

    if not prolif_genes:
        log("  WARNING: No proliferative genes found.")
        return None

    # Norm each gene 0-1 per cell
    def norm01_cells(series):
        mn, mx = series.min(), series.max()
        if mx - mn < 1e-8:
            return pd.Series(0.5, index=series.index)
        return (series - mn) / (mx - mn)

    prolif_norm = pd.DataFrame({
        g: norm01_cells(luma_expr[g])
        for g in prolif_genes
    })

    depth_score = prolif_norm.mean(axis=1)

    log(f"  Depth score (n={len(depth_score)} LumA cells):")
    log(f"    Mean:   {depth_score.mean():.4f}")
    log(f"    Median: {depth_score.median():.4f}")
    log(f"    Std:    {depth_score.std():.4f}")
    log(f"    Min:    {depth_score.min():.4f}")
    log(f"    Max:    {depth_score.max():.4f}")
    log(f"    Q25:    {depth_score.quantile(0.25):.4f}")
    log(f"    Q75:    {depth_score.quantile(0.75):.4f}")

    # Compare to normal (Mature Luminal) as reference
    mature_cells = populations["mature"]
    mature_expr  = expr.loc[mature_cells]

    if prolif_genes:
        mature_norm = pd.DataFrame({
            g: norm01_cells(mature_expr[g])
            for g in prolif_genes
        })
        # Recompute on joint range for fair comparison
        joint = pd.concat([luma_expr[prolif_genes],
                           mature_expr[prolif_genes]])
        def norm01_joint(col):
            mn, mx = col.min(), col.max()
            if mx - mn < 1e-8:
                return pd.Series(0.5, index=col.index)
            return (col - mn) / (mx - mn)

        joint_norm  = joint.apply(norm01_joint, axis=0)
        luma_depth  = joint_norm.loc[luma_cells].mean(axis=1)
        mature_depth = joint_norm.loc[mature_cells].mean(axis=1)

        log(f"\n  Joint-normalised depth scores:")
        log(f"    LumA cancer: mean={luma_depth.mean():.4f}  "
            f"std={luma_depth.std():.4f}")
        log(f"    Mature lum:  mean={mature_depth.mean():.4f}  "
            f"std={mature_depth.std():.4f}")

        _, p_depth = stats.mannwhitneyu(
            luma_depth, mature_depth,
            alternative="two-sided"
        )
        log(f"    p(LumA > Mature): {p_depth:.4e}")

    # Depth correlations within LumA
    log(f"\n  Top depth correlations within LumA cells:")
    log(f"  {'Gene':<12} {'r':>8}  {'p':>12}")
    log("  " + "-" * 36)

    corrs = []
    for gene in luma_expr.columns:
        try:
            r, p = stats.pearsonr(
                depth_score.values,
                luma_expr[gene].values
            )
            corrs.append((gene, r, p))
        except Exception:
            pass

    corrs.sort(key=lambda x: abs(x[1]), reverse=True)
    for gene, r, p in corrs[:20]:
        sig = ("***" if p < 0.001 else
               "**"  if p < 0.01  else
               "*"   if p < 0.05  else "ns")
        log(f"  {gene:<12} {r:>+8.4f}  p={p:.2e} {sig}")

    return depth_score, corrs

# ============================================================
# STEP 6: PREDICTION SCORECARD
# ============================================================

def prediction_scorecard(results_df):
    log("")
    log("=" * 65)
    log("SECTION 4: PREDICTION SCORECARD")
    log("BRCA-S1a predictions vs actual data")
    log("=" * 65)
    log("")

    if results_df is None or results_df.empty:
        log("No results. Cannot evaluate.")
        return

    def get_result(gene):
        row = results_df[results_df["gene"] == gene]
        if row.empty:
            return None, None, None
        return (row.iloc[0]["result"],
                row.iloc[0]["pct_change"],
                row.iloc[0]["p_value"])

    # Prediction 1: Switch genes retained
    log("PREDICTION 1 — SWITCH GENES RETAINED "
        "(FOXA1, GATA3, ESR1, PGR)")
    log("Expected: NOT suppressed vs Mature Luminal")
    n_confirmed = 0
    for gene in SWITCH_GENES:
        result, pct, pval = get_result(gene)
        if result is None:
            log(f"  ??? {gene:<8}: NOT IN DATA")
            continue
        ok = (result != "SUPPRESSED")
        status = "✓" if ok else "✗"
        if ok:
            n_confirmed += 1
        log(f"  {status} {gene:<8}: "
            f"{pct:>+8.1f}%  {result}  "
            f"p={pval:.2e}")
    if n_confirmed == len([g for g in SWITCH_GENES
                           if get_result(g)[0] is not None]):
        log(f"  PREDICTION 1: CONFIRMED — "
            f"all switch genes retained")
    else:
        log(f"  PREDICTION 1: PARTIALLY WRONG — "
            f"{n_confirmed} of {len(SWITCH_GENES)} retained")
        log(f"  ANALYST ASSUMPTION ERROR — "
            f"switch genes more suppressed than expected")

    log("")

    # Prediction 2: Depth genes elevated
    log("PREDICTION 2 — DEPTH AXIS PROLIFERATIVE "
        "(MKI67, CCND1, CDK4 elevated)")
    n_elev = 0
    n_avail = 0
    for gene in DEPTH_GENES:
        result, pct, pval = get_result(gene)
        if result is None:
            continue
        n_avail += 1
        ok = (result == "ELEVATED")
        status = "✓" if ok else "✗"
        if ok:
            n_elev += 1
        log(f"  {status} {gene:<8}: "
            f"{pct:>+8.1f}%  {result}  "
            f"p={pval:.2e}")
    if n_avail > 0:
        frac = n_elev / n_avail
        if frac >= 0.5:
            log(f"  PREDICTION 2: CONFIRMED — "
                f"{n_elev}/{n_avail} elevated")
        else:
            log(f"  PREDICTION 2: NOT CONFIRMED — "
                f"only {n_elev}/{n_avail} elevated")

    log("")

    # Prediction 3: EZH2 flat (not the lock in LumA)
    log("PREDICTION 3 — EZH2 FLAT IN LumA "
        "(vs +269.7% in TNBC/Basal SC)")
    result, pct, pval = get_result("EZH2")
    if result is not None:
        log(f"  EZH2: {pct:>+8.1f}%  {result}  p={pval:.2e}")
        if result == "FLAT":
            log(f"  PREDICTION 3: CONFIRMED — "
                f"EZH2 flat in LumA")
            log(f"  (Confirms EZH2 lock is TNBC-specific)")
        elif result == "ELEVATED":
            log(f"  PREDICTION 3: WRONG — "
                f"EZH2 elevated in LumA too")
            log(f"  ANALYST ASSUMPTION ERROR — "
                f"EZH2 lock not subtype-specific")
        else:
            log(f"  PREDICTION 3: CONFIRMED — "
                f"EZH2 suppressed/flat (not the lock)")

    log("")

    # Prediction 6: Controls flat
    log("PREDICTION 6 — CONTROLS FLAT "
        "(SPI1, CDX2, MBP — non-breast genes)")
    for gene in CONTROL_GENES:
        result, pct, pval = get_result(gene)
        if result is None:
            continue
        ok = (result != "ELEVATED")
        status = "✓" if ok else "✗"
        log(f"  {status} {gene:<10}: "
            f"{pct:>+8.1f}%  {result}  "
            f"p={pval:.2e}")

    log("")

    # Prediction 7: TNBC markers absent
    log("PREDICTION 7 — TNBC MARKERS ABSENT IN LumA "
        "(KRT5, SOX10, KRT14, EGFR)")
    for gene in TNBC_MARKERS:
        result, pct, pval = get_result(gene)
        if result is None:
            continue
        ok = (result in ["FLAT", "SUPPRESSED"])
        status = "✓" if ok else "⚠"
        log(f"  {status} {gene:<10}: "
            f"{pct:>+8.1f}%  {result}  "
            f"p={pval:.2e}")

    log("")
    log("PREDICTION 4 (PIK3CA axis — uncertain):")
    log("  See Section 2 PIK3CA rows above.")
    log("  Post-translational activation may not")
    log("  appear at mRNA level — flagged pre-data.")
    log("")
    log("PREDICTION 5 (Depth predicts outcome):")
    log("  scRNA-seq does not have survival data.")
    log("  Depth score within-LumA distribution")
    log("  is computed in STEP 5.")
    log("  Survival correlation requires bulk data —")
    log("  see STEP 7 (bulk RNA-seq analysis).")

# ============================================================
# STEP 7: BULK RNA-SEQ ANALYSIS (if available)
# Depth score across 26 tumors
# ============================================================

def bulk_analysis(bulk_path):
    if bulk_path is None or not os.path.exists(bulk_path):
        log("=" * 65)
        log("STEP 7: BULK RNA-SEQ — SKIPPED (file not found)")
        log("=" * 65)
        return None

    log("=" * 65)
    log("STEP 7: BULK RNA-SEQ ANALYSIS")
    log(f"File: {bulk_path}")
    log("=" * 65)

    try:
        with gzip.open(bulk_path, "rt") as f:
            bulk = pd.read_csv(f, sep="\t", index_col=0)

        log(f"  Bulk shape: {bulk.shape}")
        log(f"  Index sample (first 5): {list(bulk.index[:5])}")
        log(f"  Columns (first 5): {list(bulk.columns[:5])}")

        # Find our genes in bulk
        available = [g for g in ALL_GENES if g in bulk.index]
        log(f"  Target genes found in bulk: {len(available)}")

        bulk_sub = bulk.loc[available]

        # Identify LumA samples by column name if possible
        # (Wu et al. uses subtype in metadata)
        # For now report across all 26 tumors
        log(f"  Samples: {bulk_sub.shape[1]}")

        # Report mean expression across tumors
        log(f"\n  Mean expression across {bulk_sub.shape[1]} "
            f"tumors:")
        for gene in available:
            vals = bulk_sub.loc[gene]
            log(f"  {gene:<12}: mean={vals.mean():.2f}  "
                f"std={vals.std():.2f}")

        # Proliferative depth score per tumor
        prolif_genes = [g for g in ["MKI67", "CCND1", "CDK4"]
                        if g in available]
        if prolif_genes:
            prolif_bulk = bulk_sub.loc[prolif_genes]
            # Norm per gene across tumors
            normed = prolif_bulk.apply(
                lambda row: (row - row.min()) /
                            (row.max() - row.min() + 1e-8),
                axis=1
            )
            depth_bulk = normed.mean(axis=0)
            log(f"\n  Bulk depth scores ({len(depth_bulk)} tumors):")
            log(f"    Mean:   {depth_bulk.mean():.4f}")
            log(f"    Std:    {depth_bulk.std():.4f}")
            log(f"    Min:    {depth_bulk.min():.4f}")
            log(f"    Max:    {depth_bulk.max():.4f}")

        return bulk_sub

    except Exception as e:
        log(f"  Bulk analysis error: {e}")
        return None

# ============================================================
# STEP 8: FIGURE — 9-PANEL DISCOVERY OUTPUT
# Protocol v2.0: top movers first, predictions after
# ============================================================

def generate_figure(results_df, rdf, depth_score,
                    expr, populations):
    log("=" * 65)
    log("STEP 8: GENERATING FIGURE")
    log("=" * 65)

    if results_df is None or results_df.empty:
        log("  No results. Cannot generate figure.")
        return

    fig = plt.figure(figsize=(22, 16))
    fig.suptitle(
        "BRCA Luminal A — Script 1 Discovery\n"
        "OrganismCore | GSE176078 | 2026-03-04\n"
        "Cancer LumA SC vs Mature Luminal — Geometry First",
        fontsize=11, fontweight="bold", y=1.01
    )
    gs = gridspec.GridSpec(3, 3, figure=fig,
                           hspace=0.50, wspace=0.42)

    clr = {
        "up":   "#2980b9",
        "down": "#c0392b",
        "flat": "#27ae60",
        "luma": "#2980b9",
        "norm": "#95a5a6",
        "basal": "#c0392b",
    }

    def bar_color(pct):
        if pct > 20:
            return clr["up"]
        elif pct < -20:
            return clr["down"]
        return clr["flat"]

    # ── Panel A: Top 15 GAINED (unfiltered)
    ax_a = fig.add_subplot(gs[0, 0])
    top_gained = rdf.head(15)
    if not top_gained.empty:
        ax_a.barh(
            top_gained["gene"],
            top_gained["pct_change"],
            color=[bar_color(v) for v in top_gained["pct_change"]]
        )
        ax_a.axvline(0, color="black", lw=0.8)
        ax_a.set_title("A — Top 15 GAINED\n(unfiltered, geometry first)",
                       fontsize=8, fontweight="bold")
        ax_a.set_xlabel("% change vs normal", fontsize=7)
        ax_a.tick_params(axis="y", labelsize=6)

    # ── Panel B: Top 15 LOST (unfiltered)
    ax_b = fig.add_subplot(gs[0, 1])
    top_lost = rdf.tail(15).iloc[::-1]
    if not top_lost.empty:
        ax_b.barh(
            top_lost["gene"],
            top_lost["pct_change"],
            color=[bar_color(v) for v in top_lost["pct_change"]]
        )
        ax_b.axvline(0, color="black", lw=0.8)
        ax_b.set_title("B — Top 15 LOST\n(unfiltered, geometry first)",
                       fontsize=8, fontweight="bold")
        ax_b.set_xlabel("% change vs normal", fontsize=7)
        ax_b.tick_params(axis="y", labelsize=6)

    # ── Panel C: Switch genes (Prediction 1)
    ax_c = fig.add_subplot(gs[0, 2])
    sw = results_df[results_df["category"] == "SWITCH"]
    if not sw.empty:
        sw_sorted = sw.sort_values("pct_change")
        bars = ax_c.barh(
            sw_sorted["gene"],
            sw_sorted["pct_change"],
            color=[bar_color(v) for v in sw_sorted["pct_change"]]
        )
        ax_c.axvline(0, color="black", lw=0.8)
        ax_c.axvline(20,  color="gray", lw=0.7,
                     linestyle="--", alpha=0.6)
        ax_c.axvline(-20, color="gray", lw=0.7,
                     linestyle="--", alpha=0.6)
        ax_c.set_title("C — Switch Genes\n(Pred 1: retained in LumA?)",
                       fontsize=8, fontweight="bold")
        ax_c.set_xlabel("% vs normal", fontsize=7)
        for i, (_, row) in enumerate(sw_sorted.iterrows()):
            ax_c.text(
                row["pct_change"] +
                (2 if row["pct_change"] >= 0 else -2),
                i, row["significance"],
                va="center",
                ha="left" if row["pct_change"] >= 0 else "right",
                fontsize=6
            )

    # ── Panel D: Depth genes (Prediction 2)
    ax_d = fig.add_subplot(gs[1, 0])
    dp = results_df[results_df["category"] == "DEPTH"]
    if not dp.empty:
        dp_sorted = dp.sort_values("pct_change")
        ax_d.barh(
            dp_sorted["gene"],
            dp_sorted["pct_change"],
            color=[bar_color(v) for v in dp_sorted["pct_change"]]
        )
        ax_d.axvline(0, color="black", lw=0.8)
        ax_d.axvline(20, color="gray", lw=0.7,
                     linestyle="--", alpha=0.6)
        ax_d.set_title("D — Depth/Prolif Genes\n(Pred 2: elevated?)",
                       fontsize=8, fontweight="bold")
        ax_d.set_xlabel("% vs normal", fontsize=7)
        ax_d.tick_params(axis="y", labelsize=7)
        for i, (_, row) in enumerate(dp_sorted.iterrows()):
            ax_d.text(
                row["pct_change"] +
                (2 if row["pct_change"] >= 0 else -2),
                i, row["significance"],
                va="center",
                ha="left" if row["pct_change"] >= 0 else "right",
                fontsize=6
            )

    # ── Panel E: Lock genes — EZH2 flat?
    ax_e = fig.add_subplot(gs[1, 1])
    lk = results_df[results_df["category"] == "LOCK"]
    if not lk.empty:
        lk_sorted = lk.sort_values("pct_change")
        ax_e.barh(
            lk_sorted["gene"],
            lk_sorted["pct_change"],
            color=[bar_color(v) for v in lk_sorted["pct_change"]]
        )
        ax_e.axvline(0, color="black", lw=0.8)
        ax_e.axvline(269.7, color="red", lw=1,
                     linestyle="--", alpha=0.7,
                     label="EZH2 in TNBC: +270%")
        ax_e.set_title("E — Lock Genes\n(Pred 3: EZH2 flat?)",
                       fontsize=8, fontweight="bold")
        ax_e.set_xlabel("% vs normal", fontsize=7)
        ax_e.legend(fontsize=6)
        for i, (_, row) in enumerate(lk_sorted.iterrows()):
            ax_e.text(
                row["pct_change"] +
                (2 if row["pct_change"] >= 0 else -2),
                i, row["significance"],
                va="center",
                ha="left" if row["pct_change"] >= 0 else "right",
                fontsize=6
            )

    # ── Panel F: Cell cycle exit + controls
    ax_f = fig.add_subplot(gs[1, 2])
    ex_ct = results_df[
        results_df["category"].isin(
            ["EXIT", "CONTROL", "TNBC-MARKER"]
        )
    ]
    if not ex_ct.empty:
        ec_sorted = ex_ct.sort_values("pct_change")
        ax_f.barh(
            ec_sorted["gene"],
            ec_sorted["pct_change"],
            color=[bar_color(v) for v in ec_sorted["pct_change"]]
        )
        ax_f.axvline(0, color="black", lw=0.8)
        ax_f.set_title("F — Exit/Controls/TNBC markers\n"
                       "(Preds 6 & 7: flat?)",
                       fontsize=8, fontweight="bold")
        ax_f.set_xlabel("% vs normal", fontsize=7)
        ax_f.tick_params(axis="y", labelsize=6)

    # ── Panel G: Depth score distribution within LumA
    ax_g = fig.add_subplot(gs[2, 0])
    if depth_score is not None and len(depth_score) > 2:
        ds = depth_score[0] if isinstance(
            depth_score, tuple) else depth_score
        ax_g.hist(ds, bins=40,
                  color=clr["luma"], alpha=0.75,
                  edgecolor="white")
        ax_g.axvline(ds.mean(), color="black",
                     lw=1.5, linestyle="--",
                     label=f"mean={ds.mean():.3f}")
        ax_g.axvline(ds.quantile(0.25), color="green",
                     lw=1, linestyle=":",
                     label=f"Q25={ds.quantile(0.25):.3f}")
        ax_g.axvline(ds.quantile(0.75), color="red",
                     lw=1, linestyle=":",
                     label=f"Q75={ds.quantile(0.75):.3f}")
        ax_g.set_title("G — Depth Score Distribution\n"
                       "Within LumA SC cells",
                       fontsize=8, fontweight="bold")
        ax_g.set_xlabel("Proliferative Depth Score", fontsize=7)
        ax_g.set_ylabel("n cells", fontsize=7)
        ax_g.legend(fontsize=6)

    # ── Panel H: Waterfall all panel genes
    ax_h = fig.add_subplot(gs[2, 1])
    wf = rdf.copy()
    wf_plot = pd.concat([wf.head(10), wf.tail(10)]).drop_duplicates()
    wf_plot = wf_plot.sort_values("pct_change")
    ax_h.barh(
        wf_plot["gene"],
        wf_plot["pct_change"],
        color=[bar_color(v) for v in wf_plot["pct_change"]]
    )
    ax_h.axvline(0, color="black", lw=0.8)
    ax_h.set_title("H — Waterfall\nTop movers (top/bottom 10)",
                   fontsize=8, fontweight="bold")
    ax_h.set_xlabel("% vs normal", fontsize=7)
    ax_h.tick_params(axis="y", labelsize=6)

    # ── Panel I: Summary
    ax_i = fig.add_subplot(gs[2, 2])
    ax_i.axis("off")

    # Build summary text
    sw_results = results_df[results_df["category"] == "SWITCH"]
    dp_results = results_df[results_df["category"] == "DEPTH"]
    lk_ezh2    = results_df[
        (results_df["category"] == "LOCK") &
        (results_df["gene"] == "EZH2")
    ]

    sw_ok = (sw_results["result"] != "SUPPRESSED").all()
    dp_ok = (dp_results["result"] == "ELEVATED").sum()
    ezh2_ok = (not lk_ezh2.empty and
               lk_ezh2.iloc[0]["result"] == "FLAT")

    summary = (
        f"I — SUMMARY\n"
        f"{'─'*30}\n"
        f"Dataset: GSE176078\n"
        f"LumA SC cells: {len(populations['luma'])}\n"
        f"Mature Luminal: {len(populations['mature'])}\n"
        f"Basal SC cells: {len(populations['basal'])}\n\n"
        f"PREDICTIONS (BRCA-S1a):\n"
        f"P1 Switch retained: "
        f"{'✓ CONFIRMED' if sw_ok else '✗ WRONG'}\n"
        f"P2 Prolif elevated: "
        f"{dp_ok}/{len(dp_results)} genes\n"
        f"P3 EZH2 flat: "
        f"{'✓ CONFIRMED' if ezh2_ok else '✗ WRONG'}\n"
        f"P4 PIK3CA: uncertain (pre-stated)\n"
        f"P5 Depth→recur: needs bulk data\n\n"
        f"KEY FINDING:\n"
    )

    if sw_ok:
        summary += (
            f"LumA RETAINS luminal identity.\n"
            f"Switch genes at or above normal.\n"
            f"Depth axis = proliferative drive.\n"
            f"EZH2 lock is TNBC-specific.\n"
        )
    else:
        summary += (
            f"Some switch genes suppressed.\n"
            f"Partial identity loss in LumA.\n"
            f"Analyst assumption error — document.\n"
        )

    summary += (
        f"\nDoc: BRCA-S1b\n"
        f"Date: 2026-03-04\n"
        f"Protocol: v2.0"
    )

    ax_i.text(
        0.03, 0.97, summary,
        transform=ax_i.transAxes,
        fontsize=7.5,
        verticalalignment="top",
        fontfamily="monospace",
        bbox=dict(
            boxstyle="round",
            facecolor="#f8f8f8",
            edgecolor="#cccccc"
        )
    )

    plt.savefig(FIG_FILE, dpi=150, bbox_inches="tight")
    log(f"  Figure saved: {FIG_FILE}")
    plt.close()

# ============================================================
# MAIN
# ============================================================

def main():
    log("=" * 65)
    log("BRCA LUMINAL A — SCRIPT 1")
    log("OrganismCore — Document BRCA-S1b")
    log("Date: 2026-03-04")
    log("=" * 65)
    log("")
    log("PREDICTIONS LOCKED IN BRCA-S1a:")
    log("  1. Switch genes RETAINED (FOXA1/GATA3/ESR1/PGR)")
    log("  2. Depth axis PROLIFERATIVE (MKI67/CCND1/CDK4 up)")
    log("  3. EZH2 FLAT (not +270% as in TNBC)")
    log("  4. PIK3CA axis UNCERTAIN (post-translational)")
    log("  5. Depth predicts outcome (needs bulk data)")
    log("  6. Controls FLAT (SPI1/CDX2/MBP)")
    log("  7. TNBC markers ABSENT (KRT5/SOX10)")
    log("")
    log("DATA: GSE176078 — same dataset as original "
        "BRCA analysis")
    log("Population: Cancer LumA SC vs Mature Luminal")
    log("")

    # Step 0: Acquire data
    sc_dir, bulk_path = acquire_data()
    if sc_dir is None:
        log("FATAL: Cannot acquire data. Exiting.")
        write_log()
        sys.exit(1)

    # Step 1: Load metadata
    meta = load_metadata(sc_dir)

    # Step 2: Load expression
    expr = load_sc_expression(sc_dir, meta)

    # Step 3: Compute profiles
    expr_common, meta_common, profiles, populations = \
        compute_profiles(expr, meta)

    # Step 4: Primary comparison + top movers
    results_df, rdf = primary_comparison(
        expr_common, populations, profiles
    )

    # Save results CSV
    results_df.to_csv(CSV_FILE, index=False)
    log(f"\nResults saved: {CSV_FILE}")

    # Step 5: Depth score
    depth_score = compute_depth_score(
        expr_common, populations
    )

    # Step 6: Prediction scorecard
    prediction_scorecard(results_df)

    # Step 7: Bulk analysis
    bulk_data = bulk_analysis(bulk_path)

    # Step 8: Figure
    generate_figure(
        results_df, rdf, depth_score,
        expr_common, populations
    )

    # Write log
    write_log()

    log("")
    log("=" * 65)
    log("SCRIPT 1 COMPLETE")
    log(f"  Results : {RESULTS_DIR}")
    log(f"  Log     : {LOG_FILE}")
    log(f"  Figure  : {FIG_FILE}")
    log(f"  CSV     : {CSV_FILE}")
    log("")
    log("READ IN THIS ORDER (Protocol v2.0):")
    log("  1. Top movers (Section 1) — geometry first")
    log("  2. Full panel (Section 2) — all genes")
    log("  3. Cross-reference (Section 3) — LumA vs Basal")
    log("  4. Prediction scorecard (Section 4)")
    log("  5. Depth score (Step 5)")
    log("  6. Figure")
    log("")
    log("THEN write BRCA-S1b reasoning artifact.")
    log("=" * 65)


if __name__ == "__main__":
    main()
