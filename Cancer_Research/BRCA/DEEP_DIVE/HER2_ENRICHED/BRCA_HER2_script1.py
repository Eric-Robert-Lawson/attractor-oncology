"""
BRCA HER2-ENRICHED — SCRIPT 1
OrganismCore — Document BRCA-S3a/b | 2026-03-04

SELF-CONTAINED: Downloads GSE176078 to its own local
directory. Does NOT search for or reuse data from
LumA or TNBC directories.

GEOMETRY-FIRST (Protocol v2.0):
  This script is a DISCOVERY script.
  It reveals what the data contains before
  testing predictions. The unfiltered top movers
  are printed first, before any panel is applied.
  The panel then provides prediction tests.

PREDICTIONS TO TEST (BRCA-S3a):
  P1:  FOXA1 suppressed but NOT near-zero (-30% to -60%)
  P2:  GATA3 suppressed but partially retained
  P3:  ESR1 near-zero (ER-negative, >-85%)
  P4:  PGR near-zero (PR-negative)
  P5:  ERBB2 is the LARGEST elevated gene (>+500%)
  P6:  GRB7 co-elevated (17q12 amplicon confirmation)
  P7:  SOX10 NOT elevated (<+50%) — not a Type 2 geometry
  P8:  KRT5 NOT elevated (<+100%)
  P9:  EZH2 elevated +100% to +200% (less than TNBC +270%)
  P10: MKI67 elevated (Grade 3 biology)
  P11: HER2 geometry occupies intermediate PCA space
       between LumA and TNBC

NOVEL GEOMETRY QUESTION:
  Is HER2-enriched a "slope arrest" Type 1 variant?
  Or does it share Type 2 features with TNBC?
  The unfiltered scan answers this before predictions
  are checked.

CROSS-SUBTYPE COMPARISON:
  This is the first script to place all three
  analysed subtypes (LumA, HER2, TNBC) in the
  same geometric space using the same dataset
  and reference population.

FIX vs v1:
  define_populations() rebuilt to use the
  metadata 'celltype_subset' column directly
  by iterating over all barcodes and matching
  against the metadata DataFrame by index
  position. Avoids ct_map key mismatch when
  the MTX barcodes and metadata row order
  differ.
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
# ============================================================

GEO_ACCESSION    = "GSE176078"
SC_TARFILE       = "GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"
SC_URL           = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/"
    "GSE176078/suppl/"
    "GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"
)
SC_INNER_DIR     = "Wu_etal_2021_BRCA_scRNASeq"
SC_FILES_INSIDE  = [
    "count_matrix_sparse.mtx",
    "count_matrix_genes.tsv",
    "count_matrix_barcodes.tsv",
    "metadata.csv",
]

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
BASE_DIR    = os.path.join(SCRIPT_DIR, "HER2_s1_analysis")
DATA_DIR    = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
LOG_FILE    = os.path.join(RESULTS_DIR, "her2_s1_log.txt")
FIG_FILE    = os.path.join(RESULTS_DIR, "her2_s1_figure.png")
CSV_FILE    = os.path.join(RESULTS_DIR, "her2_s1_saddle.csv")
CACHE_FILE  = os.path.join(RESULTS_DIR, "expr_cache_her2_s1.csv")
CROSS_FILE  = os.path.join(RESULTS_DIR, "her2_s1_cross_subtype.csv")
TOP_MOVERS_FILE = os.path.join(RESULTS_DIR, "her2_s1_top_movers.csv")

for d in [BASE_DIR, DATA_DIR, RESULTS_DIR]:
    os.makedirs(d, exist_ok=True)

# ============================================================
# CELL TYPE LABELS
# ============================================================

CT_COL        = "celltype_subset"
CANCER_HER2   = "Cancer Her2 SC"
MATURE_LUM    = "Mature Luminal"
LUMINAL_PROG  = "Luminal Progenitors"
MYOEPITH      = "Myoepithelial"
CANCER_LUMA   = "Cancer LumA SC"
CANCER_BASAL  = "Cancer Basal SC"

# ============================================================
# GENE PANELS
# ============================================================

SWITCH_GENES  = ["ESR1", "FOXA1", "GATA3", "SPDEF", "PGR", "KRT8", "KRT18", "CDH1"]
HER2_AMPLICON = ["ERBB2", "GRB7", "STARD3", "MIEN1"]
FA_MARKERS    = ["KRT5", "KRT14", "SOX10", "FOXC1", "VIM", "CDH3", "ZEB1", "ZEB2", "SNAI1"]
EPIGENETIC    = ["EZH2", "EED", "SUZ12", "HDAC1", "HDAC2", "KDM1A", "DNMT3A"]
PROLIFERATION = ["MKI67", "TOP2A", "PCNA", "CCNB1", "CDK2"]
SIGNALING     = ["PIK3CA", "AKT1", "MTOR", "PTEN", "EGFR", "ERBB3"]
DNA_REPAIR    = ["BRCA1", "BRCA2", "TP53", "RB1", "PARP1"]
HETEROGENEITY = ["AR", "CLDN3", "CLDN4", "AURKA"]
CONTROLS      = ["CDX2", "SPI1", "MBP", "NKX2-1", "OLIG2"]

ALL_GENES = list(dict.fromkeys(
    SWITCH_GENES + HER2_AMPLICON + FA_MARKERS +
    EPIGENETIC + PROLIFERATION + SIGNALING +
    DNA_REPAIR + HETEROGENEITY + CONTROLS
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
# DOWNLOAD UTILITY
# ============================================================

def fetch_url(url, dest_path, retries=3):
    for attempt in range(1, retries + 1):
        try:
            log(f"  Attempt {attempt}: {os.path.basename(url)}")
            def reporthook(count, block, total):
                if total > 0 and count % 500 == 0:
                    pct = min(100, count * block * 100 // total)
                    print(f"    {pct}%...", end="\r", flush=True)
            urllib.request.urlretrieve(url, dest_path, reporthook)
            print()
            log(f"  Downloaded: {dest_path} ({os.path.getsize(dest_path)/1e6:.1f} MB)")
            return True
        except Exception as e:
            log(f"  Attempt {attempt} failed: {e}")
            if os.path.exists(dest_path):
                os.remove(dest_path)
            if attempt < retries:
                time.sleep(5)
    return False

# ============================================================
# STEP 0 — DATA ACQUISITION
# ============================================================

def acquire_data():
    log("=" * 65)
    log("STEP 0: DATA ACQUISITION — SELF-CONTAINED")
    log(f"  All data downloaded to: {DATA_DIR}")
    log("=" * 65)

    sc_dir = os.path.join(DATA_DIR, SC_INNER_DIR)
    files_present = all(
        os.path.exists(os.path.join(sc_dir, f))
        for f in SC_FILES_INSIDE
    )

    if files_present:
        log("  scRNA-seq files already present.")
        for f in SC_FILES_INSIDE:
            p = os.path.join(sc_dir, f)
            log(f"  OK: {f} ({os.path.getsize(p)/1e6:.1f} MB)")
        return sc_dir

    tar_dest = os.path.join(DATA_DIR, SC_TARFILE)
    if not os.path.exists(tar_dest):
        log(f"  Downloading from GEO FTP...")
        success = fetch_url(SC_URL, tar_dest)
        if not success:
            log("  FATAL: Download failed.")
            sys.exit(1)

    log(f"  Extracting {SC_TARFILE} ...")
    with tarfile.open(tar_dest, "r:gz") as tf:
        members = tf.getmembers()
        log(f"  Archive members: {len(members)}")
        for m in members:
            log(f"    {m.name}  ({m.size/1e6:.1f} MB)")
        tf.extractall(DATA_DIR)

    if not all(os.path.exists(os.path.join(sc_dir, f)) for f in SC_FILES_INSIDE):
        # Try root of DATA_DIR
        if all(os.path.exists(os.path.join(DATA_DIR, f)) for f in SC_FILES_INSIDE):
            sc_dir = DATA_DIR
        else:
            log("  FATAL: Expected files not found after extraction.")
            sys.exit(1)

    for f in SC_FILES_INSIDE:
        p = os.path.join(sc_dir, f)
        log(f"  OK: {f} ({os.path.getsize(p)/1e6:.1f} MB)")

    return sc_dir

# ============================================================
# STEP 1 — LOAD METADATA
# Returns metadata DataFrame with barcodes as the index.
# The index will be used for all downstream alignment.
# ============================================================

def load_metadata(sc_dir):
    log("")
    log("=" * 65)
    log("STEP 1: LOAD METADATA")
    log("=" * 65)

    meta_path = os.path.join(sc_dir, "metadata.csv")
    meta = pd.read_csv(meta_path)
    log(f"  Rows loaded: {len(meta)}")
    log(f"  Columns: {list(meta.columns)}")

    # ── Identify the barcode column ─────────────────────────
    # Wu et al. metadata has barcodes as the first unnamed column
    # (reads in as 'Unnamed: 0') or as the index itself.
    # We normalise: set barcodes as the DataFrame index.
    first_col = meta.columns[0]
    log(f"  First column (barcode candidate): '{first_col}'")
    log(f"  Sample values: {list(meta[first_col].values[:3])}")

    meta = meta.set_index(first_col)
    meta.index.name = "barcode"
    log(f"  Barcodes set as index. Index sample: {list(meta.index[:3])}")

    # ── Confirm CT_COL present ──────────────────────────────
    if CT_COL not in meta.columns:
        log(f"  WARNING: '{CT_COL}' not found. Available: {list(meta.columns)}")
        # Try common alternatives
        for candidate in ["celltype_subset", "celltype_minor", "celltype_major", "cell_type"]:
            if candidate in meta.columns:
                log(f"  Using '{candidate}' instead.")
                meta[CT_COL] = meta[candidate]
                break
        else:
            log("  FATAL: No cell type column found.")
            sys.exit(1)

    log(f"\n  Cell type distribution ({CT_COL}):")
    dist = meta[CT_COL].value_counts()
    for ct, n in dist.items():
        log(f"    {ct:<45} n={n}")

    her2_n = (meta[CT_COL] == CANCER_HER2).sum()
    log(f"\n  {CANCER_HER2}: n={her2_n}")
    if her2_n == 0:
        log("  FATAL: No Cancer Her2 SC cells found.")
        sys.exit(1)

    return meta

# ============================================================
# STEP 2 — LOAD EXPRESSION MATRIX
# The MTX barcodes file defines the row order of the matrix.
# We read those barcodes explicitly and use them as the
# expression DataFrame index — then align with metadata by
# index intersection. This avoids any column-order mismatch.
# ============================================================

def load_expression(sc_dir, meta):
    log("")
    log("=" * 65)
    log("STEP 2: LOAD EXPRESSION MATRIX")
    log("=" * 65)

    if os.path.exists(CACHE_FILE):
        log(f"  Loading cache: {CACHE_FILE}")
        expr = pd.read_csv(CACHE_FILE, index_col=0)
        log(f"  Cache shape: {expr.shape}")
        log(f"  Cache index sample: {list(expr.index[:3])}")
        missing = [g for g in ALL_GENES if g not in expr.columns]
        if missing:
            log(f"  Missing genes in cache: {missing}")
            log("  Re-extracting from MTX...")
            expr = extract_genes_from_mtx(sc_dir, ALL_GENES)
            expr.to_csv(CACHE_FILE)
        else:
            log("  All required genes in cache.")
        return expr

    log("  No cache found. Extracting from MTX (~2 min)...")
    expr = extract_genes_from_mtx(sc_dir, ALL_GENES)
    expr.to_csv(CACHE_FILE)
    log(f"  Cache saved: {CACHE_FILE}")
    return expr


def extract_genes_from_mtx(sc_dir, gene_list):
    mtx_path      = os.path.join(sc_dir, "count_matrix_sparse.mtx")
    genes_path    = os.path.join(sc_dir, "count_matrix_genes.tsv")
    barcodes_path = os.path.join(sc_dir, "count_matrix_barcodes.tsv")

    # ── Gene names ──────────────────────────────────────────
    genes_df   = pd.read_csv(genes_path, sep="\t", header=None)
    gene_names = genes_df.iloc[:, 1].values if genes_df.shape[1] >= 2 \
                 else genes_df.iloc[:, 0].values
    log(f"  Genes in matrix: {len(gene_names)}")

    # ── Barcodes — read from the barcodes file ───────────────
    # These are the EXACT index values for the expression matrix.
    barcodes_df = pd.read_csv(barcodes_path, sep="\t", header=None)
    barcodes    = barcodes_df.iloc[:, 0].values.astype(str)
    log(f"  Barcodes in matrix: {len(barcodes)}")
    log(f"  Barcode sample: {list(barcodes[:3])}")

    # ── Gene index lookup ────────────────────────────────────
    target_idx = {}
    for g in gene_list:
        matches = np.where(gene_names == g)[0]
        if len(matches) > 0:
            target_idx[g] = int(matches[0])
        else:
            log(f"  WARNING: {g} not found in gene list")

    found_genes = list(target_idx.keys())
    log(f"  Genes found: {len(found_genes)} / {len(gene_list)}")

    # ── Parse MTX ───────────────────────────────────────────
    n_cells  = len(barcodes)
    n_found  = len(found_genes)
    gene_pos = {target_idx[g]: i for i, g in enumerate(found_genes)}
    data_arr = np.zeros((n_cells, n_found), dtype=np.float32)

    log("  Parsing MTX sparse matrix (streaming)...")
    with open(mtx_path, "r") as f:
        line = f.readline()
        while line.startswith("%"):
            line = f.readline()
        parts = line.strip().split()
        n_genes_mtx  = int(parts[0])
        n_cells_mtx  = int(parts[1])
        n_entries    = int(parts[2])
        log(f"  MTX: {n_genes_mtx} genes × {n_cells_mtx} cells, {n_entries} entries")

        count = 0
        for line in f:
            parts = line.split()
            row = int(parts[0]) - 1   # gene index (1-based in MTX)
            col = int(parts[1]) - 1   # cell index (1-based in MTX)
            val = float(parts[2])
            if row in gene_pos:
                data_arr[col, gene_pos[row]] = val
            count += 1
            if count % 5_000_000 == 0:
                print(f"    {count:,} / {n_entries:,} entries...",
                      end="\r", flush=True)
    print()

    expr = pd.DataFrame(data_arr, index=barcodes, columns=found_genes)
    log(f"  Expression matrix: {expr.shape}")
    return expr

# ============================================================
# STEP 3 — NORMALISE
# ============================================================

def normalise(expr):
    log("")
    log("=" * 65)
    log("STEP 3: NORMALISE (log1p CPM)")
    log("=" * 65)
    totals = expr.sum(axis=1).replace(0, 1)
    expr_norm = expr.div(totals, axis=0) * 1e4
    expr_norm = np.log1p(expr_norm)
    log(f"  Shape: {expr_norm.shape}")
    return expr_norm

# ============================================================
# STEP 4 — DEFINE POPULATIONS
#
# KEY FIX: align expression and metadata by index intersection.
# Both expr.index (from barcodes file) and meta.index (set
# from first column in load_metadata) must use the same
# barcode strings. We intersect them explicitly and diagnose
# any mismatch before failing.
# ============================================================

def define_populations(expr, meta):
    log("")
    log("=" * 65)
    log("STEP 4: DEFINE POPULATIONS")
    log("=" * 65)

    # ── Diagnose index overlap ───────────────────────────────
    expr_idx = set(expr.index.astype(str))
    meta_idx = set(meta.index.astype(str))
    overlap  = expr_idx & meta_idx
    log(f"  Expression barcodes: {len(expr_idx)}")
    log(f"  Metadata barcodes:   {len(meta_idx)}")
    log(f"  Overlap:             {len(overlap)}")

    if len(overlap) == 0:
        # Diagnose the mismatch
        log("")
        log("  OVERLAP IS ZERO — diagnosing barcode format mismatch...")
        e_sample = sorted(expr_idx)[:5]
        m_sample = sorted(meta_idx)[:5]
        log(f"  Expression index samples: {e_sample}")
        log(f"  Metadata index samples:   {m_sample}")

        # Common fix 1: metadata barcodes have a suffix like "-1"
        # Expression barcodes are plain. Try stripping suffix.
        meta_stripped = {b.rsplit("-", 1)[0]: b for b in meta_idx}
        overlap_stripped = expr_idx & set(meta_stripped.keys())
        if len(overlap_stripped) > 1000:
            log(f"  Fix found: stripping '-N' suffix from metadata barcodes")
            log(f"  Overlap after strip: {len(overlap_stripped)}")
            meta.index = meta.index.astype(str).str.rsplit("-", n=1).str[0]
            meta_idx   = set(meta.index.astype(str))
            overlap    = expr_idx & meta_idx
            log(f"  Overlap after fix: {len(overlap)}")
        else:
            # Common fix 2: expression barcodes have suffix, metadata plain
            expr_stripped = {b.rsplit("-", 1)[0]: b for b in expr_idx}
            overlap_stripped2 = meta_idx & set(expr_stripped.keys())
            if len(overlap_stripped2) > 1000:
                log(f"  Fix found: stripping suffix from expression index")
                expr.index = expr.index.astype(str).str.rsplit("-", n=1).str[0]
                expr_idx   = set(expr.index.astype(str))
                overlap    = expr_idx & meta_idx
                log(f"  Overlap after fix: {len(overlap)}")
            else:
                log("  FATAL: Cannot reconcile barcode formats.")
                log(f"  Expression sample: {sorted(expr_idx)[:10]}")
                log(f"  Metadata sample:   {sorted(meta_idx)[:10]}")
                write_log()
                sys.exit(1)

    # ── Build cell type map from overlapping barcodes ────────
    # Use only barcodes present in both expression and metadata.
    common = sorted(overlap)
    meta_aligned = meta.loc[meta.index.astype(str).isin(overlap)]

    # Convert to dict: barcode -> cell type
    ct_map = dict(zip(
        meta_aligned.index.astype(str),
        meta_aligned[CT_COL].astype(str)
    ))

    log(f"  Barcodes used for population assignment: {len(ct_map)}")

    def get_cells(label):
        return [b for b in expr.index.astype(str) if ct_map.get(b) == label]

    her2  = get_cells(CANCER_HER2)
    lum   = get_cells(MATURE_LUM)
    prog  = get_cells(LUMINAL_PROG)
    myo   = get_cells(MYOEPITH)
    luma  = get_cells(CANCER_LUMA)
    tnbc  = get_cells(CANCER_BASAL)

    log("")
    log(f"  {CANCER_HER2:<45} n={len(her2)}")
    log(f"  {MATURE_LUM:<45} n={len(lum)}")
    log(f"  {LUMINAL_PROG:<45} n={len(prog)}")
    log(f"  {MYOEPITH:<45} n={len(myo)}")
    log(f"  {CANCER_LUMA:<45} n={len(luma)}")
    log(f"  {CANCER_BASAL:<45} n={len(tnbc)}")

    # ── Hard failures ────────────────────────────────────────
    if len(her2) < 100:
        log("")
        log("  FATAL: Fewer than 100 HER2 cells found.")
        log("  CT_COL cell type values found in metadata:")
        log(f"  {sorted(set(ct_map.values()))}")
        write_log()
        sys.exit(1)
    if len(lum) < 50:
        log("  FATAL: Fewer than 50 Mature Luminal cells found.")
        write_log()
        sys.exit(1)

    return expr, pops_dict(her2, lum, prog, myo, luma, tnbc)


def pops_dict(her2, lum, prog, myo, luma, tnbc):
    return {
        "her2": her2,
        "lum":  lum,
        "prog": prog,
        "myo":  myo,
        "luma": luma,
        "tnbc": tnbc,
    }

# ============================================================
# STEP 5 — UNFILTERED DISCOVERY SCAN
# ============================================================

def unfiltered_discovery(expr, pops):
    log("")
    log("=" * 65)
    log("STEP 5: UNFILTERED DISCOVERY SCAN")
    log("  Geometry first. No panel imposed.")
    log("  HER2 vs Mature Luminal — all expressed genes.")
    log("=" * 65)

    her2_expr = expr.loc[pops["her2"]]
    lum_expr  = expr.loc[pops["lum"]]

    log(f"  HER2 cells: {len(pops['her2'])}")
    log(f"  Lum cells:  {len(pops['lum'])}")
    log(f"  Genes:      {expr.shape[1]}")

    results = []
    for g in expr.columns:
        h = her2_expr[g].values
        l = lum_expr[g].values
        mean_h = float(np.mean(h))
        mean_l = float(np.mean(l))
        if mean_l < 0.01 and mean_h < 0.01:
            continue
        pct = (mean_h - mean_l) / mean_l * 100 if mean_l > 0 else (
            float("inf") if mean_h > 0 else 0.0
        )
        try:
            _, pval = stats.mannwhitneyu(h, l, alternative="two-sided")
        except Exception:
            pval = 1.0
        results.append({"gene": g, "her2_mean": mean_h,
                         "lum_mean": mean_l, "pct_change": pct, "pval": pval})

    rdf = pd.DataFrame(results)
    rdf["abs_pct"] = rdf["pct_change"].abs()
    rdf = rdf.sort_values("abs_pct", ascending=False)
    rdf.to_csv(TOP_MOVERS_FILE, index=False)

    def fmt_p(p):
        if p < 1e-100: return "p<1e-100 ***"
        if p < 1e-10:  return f"p={p:.2e} ***"
        if p < 1e-5:   return f"p={p:.2e} **"
        if p < 0.05:   return f"p={p:.2e} *"
        return f"p={p:.3f}"

    gained = rdf[rdf["pct_change"] > 0].head(30)
    lost   = rdf[rdf["pct_change"] < 0].sort_values("pct_change").head(30)

    log("")
    log("  TOP 30 GAINED IN HER2 vs Mature Luminal:")
    log(f"  {'Gene':<12} {'HER2':>8} {'Lum':>8} {'Change':>10}  p-value")
    log("  " + "-" * 58)
    for _, row in gained.iterrows():
        log(f"  {row['gene']:<12} {row['her2_mean']:>8.4f} "
            f"{row['lum_mean']:>8.4f} {row['pct_change']:>+9.1f}%  {fmt_p(row['pval'])}")

    log("")
    log("  TOP 30 LOST IN HER2 vs Mature Luminal:")
    log(f"  {'Gene':<12} {'HER2':>8} {'Lum':>8} {'Change':>10}  p-value")
    log("  " + "-" * 58)
    for _, row in lost.iterrows():
        log(f"  {row['gene']:<12} {row['her2_mean']:>8.4f} "
            f"{row['lum_mean']:>8.4f} {row['pct_change']:>+9.1f}%  {fmt_p(row['pval'])}")

    # Geometry read
    log("")
    log("  GEOMETRY READ (unfiltered):")
    top_gained = gained["gene"].tolist()

    erbb2_row  = rdf[rdf["gene"] == "ERBB2"]
    foxa1_row  = rdf[rdf["gene"] == "FOXA1"]
    sox10_row  = rdf[rdf["gene"] == "SOX10"]
    krt5_row   = rdf[rdf["gene"] == "KRT5"]
    esr1_row   = rdf[rdf["gene"] == "ESR1"]

    if not erbb2_row.empty:
        rank = rdf.index.get_loc(erbb2_row.index[0]) + 1
        log(f"  ERBB2: rank #{rank}  change={erbb2_row['pct_change'].values[0]:+.1f}%")
    if not foxa1_row.empty:
        fc = foxa1_row["pct_change"].values[0]
        log(f"  FOXA1: {fc:+.1f}%  (TNBC was -81% | LumA was partial)")
        if fc > -70:
            log("  → Partial luminal retention — slope arrest geometry supported")
        else:
            log("  → Deep suppression — reassess geometry type")
    if not esr1_row.empty:
        log(f"  ESR1:  {esr1_row['pct_change'].values[0]:+.1f}%  (expected < -85%)")
    if not sox10_row.empty:
        fc = sox10_row["pct_change"].values[0]
        log(f"  SOX10: {fc:+.1f}%  (TNBC was +1323% — neural crest signal)")
        if fc < 50:
            log("  → No neural crest signal. Type 2 geometry NOT present.")
        elif fc < 300:
            log("  → Moderate SOX10 elevation. Partial Type 2 component possible.")
        else:
            log("  → SOX10 dominates. Type 2 component present — composite type evaluation needed.")
    if not krt5_row.empty:
        log(f"  KRT5:  {krt5_row['pct_change'].values[0]:+.1f}%  (TNBC was +508%)")

    return rdf

# ============================================================
# STEP 6 — PANEL ANALYSIS AND PREDICTION TESTS
# ============================================================

def panel_analysis(expr, pops):
    log("")
    log("=" * 65)
    log("STEP 6: PANEL ANALYSIS — PREDICTION TESTS")
    log("  HER2 vs Mature Luminal")
    log("=" * 65)

    her2_expr = expr.loc[pops["her2"]]
    lum_expr  = expr.loc[pops["lum"]]

    def fmt_p(p):
        if p < 1e-100: return "p<1e-100 ***"
        if p < 1e-10:  return f"p={p:.2e} ***"
        if p < 1e-5:   return f"p={p:.2e} **"
        if p < 0.05:   return f"p={p:.2e} *"
        return f"p={p:.3f}"

    panel_groups = [
        ("LUMINAL SWITCH GENES",        SWITCH_GENES),
        ("HER2 AMPLICON",               HER2_AMPLICON),
        ("FALSE ATTRACTOR / BASAL",     FA_MARKERS),
        ("EPIGENETIC AXIS",             EPIGENETIC),
        ("PROLIFERATION",               PROLIFERATION),
        ("PI3K / SIGNALING",            SIGNALING),
        ("DNA REPAIR / COMPOSITE TYPE", DNA_REPAIR),
        ("HETEROGENEITY PROBES",        HETEROGENEITY),
        ("CROSS-CANCER CONTROLS",       CONTROLS),
    ]

    all_results = []
    for group_name, gene_list in panel_groups:
        log(f"\n  --- {group_name} ---")
        log(f"  {'Gene':<12} {'HER2':>8} {'Lum':>8} {'Change':>10}  p-value")
        log("  " + "-" * 58)
        avail = [g for g in gene_list if g in expr.columns]
        if not avail:
            log("  (none found in matrix)")
            continue
        for g in avail:
            h = her2_expr[g].values
            l = lum_expr[g].values
            mean_h = float(np.mean(h))
            mean_l = float(np.mean(l))
            pct = (mean_h - mean_l) / mean_l * 100 if mean_l > 0.001 else (
                float("inf") if mean_h > 0.001 else 0.0
            )
            try:
                _, pval = stats.mannwhitneyu(h, l, alternative="two-sided")
            except Exception:
                pval = 1.0
            arrow = "↑" if pct > 0 else "↓"
            log(f"  {g:<12} {mean_h:>8.4f} {mean_l:>8.4f} "
                f"  {arrow}{abs(pct):>7.1f}%  {fmt_p(pval)}")
            all_results.append({
                "gene": g, "group": group_name,
                "her2_mean": mean_h, "lum_mean": mean_l,
                "pct_change": pct, "pval": pval,
            })

    results_df = pd.DataFrame(all_results)

    # ── PREDICTION SCORECARD ──────────────────────────────────
    log("")
    log("=" * 65)
    log("  PREDICTION SCORECARD (BRCA-S3a)")
    log("=" * 65)

    def get_pct(gene):
        row = results_df[results_df["gene"] == gene]
        return row["pct_change"].values[0] if not row.empty else None

    checks = [
        ("P1",  "FOXA1", lambda x: -70 < x < 0,   "partial suppression -30% to -70%"),
        ("P2",  "GATA3", lambda x: -60 < x < 0,   "partial suppression"),
        ("P3",  "ESR1",  lambda x: x < -85,        "near-zero >-85%"),
        ("P4",  "PGR",   lambda x: x < -80,        "near-zero >-80%"),
        ("P5",  "ERBB2", lambda x: x > 500,        "ERBB2 >+500%"),
        ("P6",  "GRB7",  lambda x: x > 200,        "GRB7 >+200% (17q12 amplicon)"),
        ("P7",  "SOX10", lambda x: x < 50,         "SOX10 <+50% (not Type 2)"),
        ("P8",  "KRT5",  lambda x: x < 100,        "KRT5 <+100% (not basal)"),
        ("P9",  "EZH2",  lambda x: 100 < x < 270,  "EZH2 +100% to +270%"),
        ("P10", "MKI67", lambda x: x > 100,        "MKI67 >+100% (Grade 3)"),
    ]

    for pred_id, gene, cond, desc in checks:
        pct = get_pct(gene)
        if pct is None:
            log(f"  {pred_id}: {gene:<8} NOT IN MATRIX (inconclusive)")
        else:
            status = "CONFIRMED" if cond(pct) else "NOT CONFIRMED"
            log(f"  {pred_id}: {gene:<8} {pct:>+8.1f}%  → {status}  [{desc}]")

    return results_df

# ============================================================
# STEP 7 — WITHIN-HER2 DEPTH ANALYSIS
# ============================================================

def depth_analysis(expr, pops):
    log("")
    log("=" * 65)
    log("STEP 7: WITHIN-HER2 DEPTH ANALYSIS")
    log("  What drives heterogeneity inside the HER2 attractor?")
    log("=" * 65)

    her2_expr = expr.loc[pops["her2"]]
    depth_genes = [g for g in ["FOXA1", "GATA3", "ESR1", "PGR", "SPDEF"]
                   if g in expr.columns]

    if not depth_genes:
        log("  No depth genes available.")
        return None, None

    log(f"  Depth score genes (low luminal = deep): {depth_genes}")

    def norm01(s):
        mn, mx = s.min(), s.max()
        return (s - mn) / (mx - mn + 1e-9)

    depth_score = np.mean(
        [1 - norm01(her2_expr[g]).values for g in depth_genes], axis=0
    )

    log(f"  Mean: {depth_score.mean():.4f}  Std: {depth_score.std():.4f}")

    corr_genes = [g for g in ALL_GENES
                  if g in expr.columns and g not in depth_genes]
    corr_results = []
    for g in corr_genes:
        vals = her2_expr[g].values
        if np.std(vals) < 1e-6:
            continue
        r, p = stats.spearmanr(depth_score, vals)
        corr_results.append({"gene": g, "r": r, "p": p})

    corr_df = pd.DataFrame(corr_results).sort_values("r", ascending=False)

    log("")
    log("  Top 15 POSITIVE r(gene, depth) — deeper = higher:")
    log(f"  {'Gene':<12} {'r':>8}  p-value")
    log("  " + "-" * 35)
    for _, row in corr_df.head(15).iterrows():
        log(f"  {row['gene']:<12} {row['r']:>+8.3f}  p={row['p']:.2e}")

    log("")
    log("  Top 15 NEGATIVE r(gene, depth) — deeper = lower:")
    log(f"  {'Gene':<12} {'r':>8}  p-value")
    log("  " + "-" * 35)
    for _, row in corr_df.tail(15).iterrows():
        log(f"  {row['gene']:<12} {row['r']:>+8.3f}  p={row['p']:.2e}")

    for g in ["ERBB2", "EZH2", "EED", "MKI67", "GRB7"]:
        row = corr_df[corr_df["gene"] == g]
        if not row.empty:
            log(f"  {g} r(depth): {row['r'].values[0]:+.3f}  p={row['p'].values[0]:.2e}")

    return pd.Series(depth_score, index=pops["her2"]), corr_df

# ============================================================
# STEP 8 — CROSS-SUBTYPE COMPARISON
# ============================================================

def cross_subtype_comparison(expr, pops):
    log("")
    log("=" * 65)
    log("STEP 8: CROSS-SUBTYPE COMPARISON")
    log("  LumA / HER2 / TNBC — same dataset, same reference")
    log("=" * 65)

    lum_expr = expr.loc[pops["lum"]]
    sub_pops  = {"LumA": pops["luma"], "HER2": pops["her2"], "TNBC": pops["tnbc"]}

    key_genes = [
        "FOXA1", "GATA3", "ESR1", "PGR",
        "ERBB2", "GRB7",
        "SOX10", "KRT5", "ZEB1", "VIM",
        "EZH2",  "EED",
        "MKI67", "TOP2A",
        "AR",    "EGFR",
        "SPI1",  "CDX2", "MBP",
    ]
    avail = [g for g in key_genes if g in expr.columns]

    log(f"\n  Reference: Mature Luminal (n={len(pops['lum'])})")
    header = f"  {'Gene':<12}" + "".join(f"  {s:>12}" for s in sub_pops)
    log(header)
    log("  " + "-" * (12 + len(sub_pops) * 14 + 4))

    cross_rows = []
    for g in avail:
        lum_mean = lum_expr[g].mean()
        row_str  = f"  {g:<12}"
        row_data = {"gene": g}
        for sub, cells in sub_pops.items():
            sub_mean = expr.loc[cells][g].mean()
            pct = (sub_mean - lum_mean) / lum_mean * 100 if lum_mean > 0.001 else 0.0
            row_str += f"  {pct:>+10.1f}%"
            row_data[sub] = pct
        log(row_str)
        cross_rows.append(row_data)

    cross_df = pd.DataFrame(cross_rows)
    cross_df.to_csv(CROSS_FILE, index=False)
    log(f"\n  Saved: {CROSS_FILE}")

    log("")
    log("  GEOMETRIC ORDERING (Predictions H15-H17):")
    for g, ascending in [("FOXA1", False), ("GATA3", False), ("EZH2", True)]:
        row = cross_df[cross_df["gene"] == g]
        if row.empty:
            continue
        luma_v = row["LumA"].values[0]
        her2_v = row["HER2"].values[0]
        tnbc_v = row["TNBC"].values[0]
        if ascending:
            ordered = luma_v < her2_v < tnbc_v
        else:
            ordered = luma_v > her2_v > tnbc_v
        status = "ORDER CONFIRMED" if ordered else "ORDER NOT AS PREDICTED"
        log(f"  {g:<8}: LumA={luma_v:+.1f}%  HER2={her2_v:+.1f}%  "
            f"TNBC={tnbc_v:+.1f}%  → {status}")

    return cross_df

# ============================================================
# STEP 9 — PCA GEOMETRY
# ============================================================

def pca_geometry(expr, pops):
    log("")
    log("=" * 65)
    log("STEP 9: PCA GEOMETRY")
    log("=" * 65)

    pca_genes = [g for g in ALL_GENES if g in expr.columns]
    if len(pca_genes) < 5:
        log("  Not enough genes. Skipping.")
        return None, None, None, None

    pop_cfg = {
        "her2": (CANCER_HER2,  "red"),
        "lum":  (MATURE_LUM,   "steelblue"),
        "prog": (LUMINAL_PROG, "cornflowerblue"),
        "luma": (CANCER_LUMA,  "goldenrod"),
        "tnbc": (CANCER_BASAL, "dimgray"),
    }

    all_cells, cell_labels, cell_colors = [], [], []
    rng = np.random.default_rng(42)
    for key, (label, color) in pop_cfg.items():
        cells = pops[key]
        if not cells:
            continue
        sample = list(rng.choice(cells, min(2000, len(cells)), replace=False))
        all_cells   += sample
        cell_labels += [label] * len(sample)
        cell_colors += [color] * len(sample)

    sub_expr = expr.loc[all_cells, pca_genes].fillna(0)
    X_scaled = StandardScaler().fit_transform(sub_expr.values)
    pca      = PCA(n_components=3, random_state=42)
    X_pca    = pca.fit_transform(X_scaled)
    var_exp  = pca.explained_variance_ratio_

    log(f"  Genes: {len(pca_genes)}   Cells: {len(all_cells)}")
    log(f"  PC1: {var_exp[0]*100:.1f}%   PC2: {var_exp[1]*100:.1f}%   PC3: {var_exp[2]*100:.1f}%")

    pca_df = pd.DataFrame(X_pca, columns=["PC1", "PC2", "PC3"])
    pca_df["label"] = cell_labels

    log("\n  Centroids:")
    for label, _ in pop_cfg.values():
        sub = pca_df[pca_df["label"] == label]
        if sub.empty:
            continue
        c = sub[["PC1", "PC2", "PC3"]].mean().values
        log(f"    {label:<45} PC1={c[0]:+.3f}  PC2={c[1]:+.3f}  PC3={c[2]:+.3f}")

    # Distance tests (H19-H20)
    def centroid(lbl):
        return pca_df[pca_df["label"] == lbl][["PC1","PC2","PC3"]].mean().values

    her2_c = centroid(CANCER_HER2)
    lum_c  = centroid(MATURE_LUM)
    tnbc_c = centroid(CANCER_BASAL)
    luma_c = centroid(CANCER_LUMA)

    d_her2 = np.linalg.norm(her2_c - lum_c)
    d_tnbc = np.linalg.norm(tnbc_c - lum_c)
    d_luma = np.linalg.norm(luma_c - lum_c)

    log(f"\n  Centroid distances from Mature Luminal:")
    log(f"    LumA → Lum: {d_luma:.4f}")
    log(f"    HER2 → Lum: {d_her2:.4f}")
    log(f"    TNBC → Lum: {d_tnbc:.4f}")

    if d_tnbc > 0:
        rel = d_her2 / d_tnbc
        log(f"\n  HER2 relative distance (TNBC=1.0): {rel:.3f}")
        if rel < 0.70:
            log("  → P11 CONFIRMED: HER2 geometrically between LumA and TNBC")
        elif rel < 1.0:
            log("  → P11 PARTIAL: HER2 closer to Luminal than TNBC but above threshold")
        else:
            log("  → P11 NOT CONFIRMED: HER2 as far or further from Luminal than TNBC")

    return X_pca, pca_df, var_exp, cell_colors

# ============================================================
# STEP 10 — FIGURE
# ============================================================

def generate_figure(expr, pops, panel_results, cross_df, pca_result):
    log("")
    log("=" * 65)
    log("STEP 10: GENERATING FIGURE")
    log("=" * 65)

    fig = plt.figure(figsize=(22, 18))
    fig.suptitle(
        "HER2-Enriched Breast Cancer — Script 1\n"
        "OrganismCore | False Attractor Framework | 2026-03-04",
        fontsize=13, fontweight="bold", y=0.98
    )
    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.38)

    # Panel A — bar chart
    ax_a = fig.add_subplot(gs[0, :2])
    if panel_results is not None and not panel_results.empty:
        key_plot = ["ERBB2","GRB7","EZH2","EED","MKI67","TOP2A","EGFR",
                    "FOXA1","GATA3","ESR1","PGR","SPDEF",
                    "SOX10","KRT5","VIM","ZEB1","AR"]
        sub = panel_results[panel_results["gene"].isin(key_plot)].copy()
        sub = sub.set_index("gene").reindex(key_plot).dropna()
        colors = ["tomato" if v > 0 else "steelblue" for v in sub["pct_change"]]
        ax_a.barh(sub.index, sub["pct_change"], color=colors)
        ax_a.axvline(0, color="black", linewidth=0.8)
        ax_a.set_xlabel("% change vs Mature Luminal")
        ax_a.set_title("HER2 vs Mature Luminal — Key Gene Panel", fontweight="bold")
        ax_a.invert_yaxis()

    # Panel B — attractor annotation
    ax_b = fig.add_subplot(gs[0, 2])
    ax_b.axis("off")
    ann = [
        "ATTRACTOR TYPE",
        "Type 1 Variant — Slope Arrest",
        "",
        "Cell of origin:",
        "  Luminal Progenitor",
        "",
        "Founding event:",
        "  ERBB2 amplicon 17q12",
        "  Copy number, not TF loss",
        "",
        "Key tests:",
        "  ERBB2 dominant signal?",
        "  SOX10 absent (not Type 2)?",
        "  FOXA1 partial (not TNBC)?",
    ]
    ax_b.text(0.05, 0.95, "\n".join(ann), transform=ax_b.transAxes,
              fontsize=8, verticalalignment="top", fontfamily="monospace",
              bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))

    # Panel C — cross-subtype heatmap
    ax_c = fig.add_subplot(gs[1, :])
    if cross_df is not None and not cross_df.empty:
        heat_genes = ["FOXA1","GATA3","ESR1","PGR","ERBB2","GRB7",
                      "SOX10","KRT5","ZEB1","VIM","EZH2","EED",
                      "MKI67","TOP2A","AR","EGFR","SPI1","CDX2"]
        plot_heat  = [g for g in heat_genes if g in cross_df["gene"].values]
        hm  = cross_df[cross_df["gene"].isin(plot_heat)].set_index("gene").reindex(plot_heat)
        cols = [c for c in ["LumA","HER2","TNBC"] if c in hm.columns]
        hm_vals = hm[cols].apply(pd.to_numeric, errors="coerce").fillna(0)
        vmax = max(abs(hm_vals.values.min()), abs(hm_vals.values.max()), 300)
        im = ax_c.imshow(hm_vals.values.T, aspect="auto",
                         cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax_c.set_yticks(range(len(cols)))
        ax_c.set_yticklabels(cols, fontsize=9)
        ax_c.set_xticks(range(len(plot_heat)))
        ax_c.set_xticklabels(plot_heat, rotation=45, ha="right", fontsize=8)
        ax_c.set_title("Cross-Subtype: LumA / HER2 / TNBC vs Mature Luminal",
                       fontweight="bold")
        plt.colorbar(im, ax=ax_c, label="% change vs Mature Luminal", shrink=0.6)

    # Panel D — PCA
    ax_d = fig.add_subplot(gs[2, :2])
    if pca_result[0] is not None:
        X_pca, pca_df, var_exp, _ = pca_result
        cmap = {CANCER_HER2:"red", MATURE_LUM:"steelblue",
                LUMINAL_PROG:"cornflowerblue", CANCER_LUMA:"goldenrod",
                CANCER_BASAL:"dimgray"}
        for label, color in cmap.items():
            sub = pca_df[pca_df["label"] == label]
            if sub.empty:
                continue
            ax_d.scatter(sub["PC1"], sub["PC2"], c=color, s=1.5,
                         alpha=0.4, rasterized=True,
                         label=label.replace("Cancer ","").replace(" SC",""))
        ax_d.set_xlabel(f"PC1 ({var_exp[0]*100:.1f}%)")
        ax_d.set_ylabel(f"PC2 ({var_exp[1]*100:.1f}%)")
        ax_d.set_title("PCA: LumA / HER2 / TNBC / Reference", fontweight="bold")
        ax_d.legend(markerscale=5, fontsize=7)

    # Panel E — scorecard
    ax_e = fig.add_subplot(gs[2, 2])
    ax_e.axis("off")
    score_lines = ["PREDICTION SCORECARD", ""]
    if panel_results is not None:
        chk = [("P1 FOXA1","FOXA1",lambda x:-70<x<0),
               ("P2 GATA3","GATA3",lambda x:-60<x<0),
               ("P3 ESR1", "ESR1", lambda x:x<-85),
               ("P4 PGR",  "PGR",  lambda x:x<-80),
               ("P5 ERBB2","ERBB2",lambda x:x>500),
               ("P6 GRB7", "GRB7", lambda x:x>200),
               ("P7 SOX10","SOX10",lambda x:x<50),
               ("P8 KRT5", "KRT5", lambda x:x<100),
               ("P9 EZH2", "EZH2", lambda x:100<x<270),
               ("P10 MKI67","MKI67",lambda x:x>100)]
        for label, gene, fn in chk:
            row = panel_results[panel_results["gene"] == gene]
            if row.empty:
                score_lines.append(f"? {label}: N/A")
            else:
                pct = row["pct_change"].values[0]
                score_lines.append(
                    f"{'✓' if fn(pct) else '✗'} {label}: {pct:+.0f}%")
    ax_e.text(0.05, 0.95, "\n".join(score_lines), transform=ax_e.transAxes,
              fontsize=8, verticalalignment="top", fontfamily="monospace",
              bbox=dict(boxstyle="round", facecolor="lightcyan", alpha=0.8))

    plt.savefig(FIG_FILE, dpi=150, bbox_inches="tight")
    plt.close()
    log(f"  Figure saved: {FIG_FILE}")

# ============================================================
# MAIN
# ============================================================

def main():
    log("=" * 65)
    log("BRCA HER2-ENRICHED — SCRIPT 1")
    log("OrganismCore — Document BRCA-S3a/b | 2026-03-04")
    log("Self-contained: HER2_s1_analysis/data/")
    log("=" * 65)

    sc_dir       = acquire_data()
    meta         = load_metadata(sc_dir)
    expr_raw     = load_expression(sc_dir, meta)
    expr         = normalise(expr_raw)
    expr, pops   = define_populations(expr, meta)
    top_movers   = unfiltered_discovery(expr, pops)
    panel_results = panel_analysis(expr, pops)
    depth_result  = depth_analysis(expr, pops)
    cross_df      = cross_subtype_comparison(expr, pops)
    pca_result    = pca_geometry(expr, pops)

    generate_figure(expr, pops, panel_results, cross_df,
                    pca_result if pca_result[0] is not None
                    else (None, None, None, None))

    if panel_results is not None:
        panel_results.to_csv(CSV_FILE, index=False)
        log(f"\n  Panel results: {CSV_FILE}")

    write_log()
    log("")
    log("=" * 65)
    log("SCRIPT 1 COMPLETE")
    log(f"  Log:        {LOG_FILE}")
    log(f"  Figure:     {FIG_FILE}")
    log(f"  Results:    {CSV_FILE}")
    log(f"  Cross:      {CROSS_FILE}")
    log(f"  Top movers: {TOP_MOVERS_FILE}")
    log("=" * 65)
    write_log()


if __name__ == "__main__":
    main()
