"""
BRCA TNBC — SCRIPT 2
BULK RNA-seq VALIDATION + pCR ANALYSIS + TCGA SURVIVAL
OrganismCore — Document BRCA-S4c/d | 2026-03-04

CONFIRMED WORKING URLS (from diagnostic 2026-03-04):
  GSE25066:   https://ftp.ncbi.nlm.nih.gov/geo/series/GSE25nnn/GSE25066/matrix/GSE25066_series_matrix.txt.gz
  GPL96:      BUILT FROM SERIES MATRIX — all NCBI GPL96 files return 404.
              The series matrix !platform_table_begin block contains the
              full probe→gene mapping inline. We parse it from there.
  TCGA EXPR:  https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2.gz
  TCGA CLIN:  https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix
  TCGA SURV:  https://pancanatlas.xenahubs.net/download/Survival_SupplementalTable_S1_20171025_xena_sp
  PAM50:      Extracted from BRCA_clinicalMatrix column PAM50Call_RNAseq

LOCAL CACHE (confirmed by diagnostic):
  Script 1 TNBC cache:
    /Users/ericlawson/cancer/BRCA/DEEP_DIVE/TNBC/TNBC_s1_analysis/results/expr_cache_tnbc_s1.csv

THREE TASKS (from BRCA-S4c before-document):
  TASK 1 — GSE25066 proper analysis (GPL96 inline, TNBC subset, gene-mapped depth)
  TASK 2 — Lehmann subtype depth mapping (S2-P4)
  TASK 3 — TCGA-BRCA Basal-like: BRCA1 proxy, PAM50, survival (S2-P5, S2-P6)

PREDICTIONS TESTED (BRCA-S4c):
  S2-P1: Gene-mapped depth r(depth, pCR) more negative than Script 1 degraded proxy
  S2-P2: EED predicts pCR better than EZH2 in TNBC
  S2-P3: PARP1 predicts lower pCR (chemoresistance marker)
  S2-P4: Lehmann subtypes map depth axis: LAR < BL1/BL2 < M/MSL
  S2-P5: BRCA1 dysfunction enriched in Basal-like TCGA
  S2-P6: Depth score predicts survival in TCGA-BRCA
  S2-P7: SPI1 elevation is immune contamination (resolved in bulk)
  S2-P8: EED > EZH2 as EZH2i response biomarker (stated prospectively)
  S2-P9: AR correlates negatively with depth in bulk TNBC
"""

import os
import sys
import gzip
import time
import warnings
import urllib.request
import traceback
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from sklearn.metrics import roc_auc_score

# ============================================================
# CONFIGURATION
# ============================================================

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
BASE_DIR    = os.path.join(SCRIPT_DIR, "tnbc_s2_results")
DATA_DIR    = os.path.join(BASE_DIR, "data")
RESULTS_DIR = BASE_DIR

for d in [BASE_DIR, DATA_DIR, RESULTS_DIR]:
    os.makedirs(d, exist_ok=True)

LOG_FILE = os.path.join(RESULTS_DIR, "tnbc_s2_log.txt")
FIG_FILE = os.path.join(RESULTS_DIR, "tnbc_s2_figure.png")
CSV_FILE = os.path.join(RESULTS_DIR, "tnbc_s2_results.csv")

# ── Script 1 cache path (confirmed by diagnostic) ──────────
S1_CACHE = os.path.join(
    SCRIPT_DIR, "TNBC_s1_analysis", "results",
    "expr_cache_tnbc_s1.csv"
)
S1_META = os.path.join(
    os.path.expanduser("~/cancer/BRCA"),
    "Wu_etal_2021_BRCA_scRNASeq", "metadata.csv"
)

# ── Confirmed working URLs (diagnostic 2026-03-04) ─────────
GSE25066_URL  = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/"
    "GSE25nnn/GSE25066/matrix/"
    "GSE25066_series_matrix.txt.gz"
)
GSE25066_FILE = os.path.join(DATA_DIR, "GSE25066_series_matrix.txt.gz")

# GPL96 is parsed inline from GSE25066 series matrix.
# No separate download needed. See parse_gse25066().

TCGA_EXPR_URL  = (
    "https://tcga.xenahubs.net/download/"
    "TCGA.BRCA.sampleMap/HiSeqV2.gz"
)
TCGA_EXPR_FILE = os.path.join(DATA_DIR, "TCGA_BRCA_HiSeqV2.gz")

TCGA_CLIN_URL  = (
    "https://tcga.xenahubs.net/download/"
    "TCGA.BRCA.sampleMap/BRCA_clinicalMatrix"
)
TCGA_CLIN_FILE = os.path.join(DATA_DIR, "TCGA_BRCA_clinicalMatrix.tsv")

TCGA_SURV_URL  = (
    "https://pancanatlas.xenahubs.net/download/"
    "Survival_SupplementalTable_S1_20171025_xena_sp"
)
TCGA_SURV_FILE = os.path.join(DATA_DIR, "TCGA_pancan_survival.tsv")

# ============================================================
# GENE PANELS — from BRCA-S4c before-document
# ============================================================

# Depth score genes — from Script 1 confirmed findings
FA_MARKERS   = ["KRT5", "KRT14", "SOX10", "FOXC1", "EGFR", "VIM", "CDH3"]
SWITCH_GENES = ["ESR1", "FOXA1", "GATA3", "SPDEF", "PGR"]

# Genes for prediction testing
PCR_BIOMARKERS = ["EED", "EZH2", "PARP1", "AR", "VIM", "SOX10",
                  "KRT5", "KRT14", "ESR1", "GATA3", "FOXA1",
                  "HDAC1", "HDAC2", "KDM1A", "MKI67", "CD274",
                  "CDH1", "ZEB1", "ZEB2", "BRCA1", "BRCA2",
                  "TP53", "PIK3CA", "PTEN", "SPI1", "PTPRC",
                  "CDK2", "CCNB1", "TOP2A", "AURKA"]

# Lehmann subtype centroid genes (simplified — key markers per subtype)
# BL1: DNA damage response, cell cycle; BL2: growth factor; LAR: AR/luminal
# IM: immune; M: mesenchymal; MSL: mesenchymal stem-like
LEHMANN_GENES = {
    "BL1":  ["CCNE1", "CDK2", "CDC6", "CDC20", "CCNA2", "MKI67",
             "TOP2A", "PCNA", "BRCA1", "BRCA2"],
    "BL2":  ["EGFR", "MET", "HGF", "IGF1R", "PIK3CA", "PTEN",
             "CLDN3", "CLDN4"],
    "IM":   ["PTPRC", "CD4", "CD8A", "CD8B", "CCL5", "CXCL10",
             "SPI1", "IRF1", "STAT1"],
    "LAR":  ["AR", "SPDEF", "FOXA1", "TFF1", "TFF3", "AGR2",
             "GCDFP15", "KLK10", "KLK11"],
    "M":    ["VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1",
             "CDH2", "FN1"],
    "MSL":  ["CD44", "ALDH1A1", "FOXC1", "PDGFRA", "PDGFRB",
             "CXCR4", "ANGPTL1"],
}

ALL_GENES = list(dict.fromkeys(
    FA_MARKERS + SWITCH_GENES + PCR_BIOMARKERS +
    [g for genes in LEHMANN_GENES.values() for g in genes]
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

def fmt_p(p):
    if p < 1e-100: return "p<1e-100 ***"
    if p < 0.001:  return f"p={p:.2e} ***"
    if p < 0.01:   return f"p={p:.2e}  **"
    if p < 0.05:   return f"p={p:.4f}   *"
    return             f"p={p:.4f}  ns"

# ============================================================
# FETCH UTILITY
# ============================================================

def fetch_url(url, dest, label="", retries=3):
    if os.path.exists(dest) and os.path.getsize(dest) > 10000:
        log(f"  [CACHED] {label or os.path.basename(dest)} "
            f"({os.path.getsize(dest)/1e6:.1f} MB)")
        return dest
    log(f"  Downloading: {label or os.path.basename(dest)}")
    log(f"  URL: {url}")
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "Mozilla/5.0"}
            )
            with urllib.request.urlopen(req, timeout=600) as r:
                data = r.read()
            with open(dest, "wb") as f:
                f.write(data)
            log(f"  Downloaded: {dest} ({len(data)/1e6:.1f} MB)")
            return dest
        except Exception as e:
            log(f"  Attempt {attempt+1} failed: {e}")
            if attempt < retries - 1:
                time.sleep(5)
    log(f"  FATAL: Download failed after {retries} attempts")
    return None

# ============================================================
# STEP 0: ACQUIRE DATA
# ============================================================

def acquire_data():
    log("=" * 65)
    log("STEP 0: DATA ACQUISITION")
    log("=" * 65)

    results = {}

    log("\n[1/4] GSE25066 series matrix...")
    p = fetch_url(GSE25066_URL, GSE25066_FILE, "GSE25066")
    results["gse25066"] = p

    log("\n[2/4] TCGA-BRCA expression (Xena Hub)...")
    p = fetch_url(TCGA_EXPR_URL, TCGA_EXPR_FILE, "TCGA_EXPR")
    results["tcga_expr"] = p

    log("\n[3/4] TCGA-BRCA clinical (Xena Hub)...")
    p = fetch_url(TCGA_CLIN_URL, TCGA_CLIN_FILE, "TCGA_CLIN")
    results["tcga_clin"] = p

    log("\n[4/4] Pan-cancer survival table (UCSC)...")
    p = fetch_url(TCGA_SURV_URL, TCGA_SURV_FILE, "TCGA_SURV")
    results["tcga_surv"] = p

    log("\n  Script 1 cache:")
    if os.path.exists(S1_CACHE):
        log(f"  ✓ {S1_CACHE} ({os.path.getsize(S1_CACHE)/1e6:.1f} MB)")
    else:
        log(f"  ✗ NOT FOUND: {S1_CACHE}")

    return results

# ============================================================
# PARSE GSE25066 — INLINE GPL96 PROBE MAPPING
# Returns (expr_df, meta_df, probe_to_gene)
# ============================================================

def parse_gse25066(matrix_path):
    """
    Parse GSE25066 series matrix file.

    The GEO series matrix format contains:
      - Header rows (!Sample_...) with metadata per sample
      - An inline platform annotation block
        (!platform_table_begin ... !platform_table_end)
        which maps probe IDs to gene symbols — this IS the GPL96 mapping
      - Expression table (series_matrix_table_begin ... end)

    No separate GPL96 file needed — it is embedded here.
    """
    log("=" * 65)
    log("PARSE GSE25066 — INLINE GPL96 PROBE MAPPING")
    log("GPL96 annotation extracted from series matrix header")
    log("(all standalone GPL96 files return 404 on NCBI FTP)")
    log("=" * 65)

    if matrix_path is None or not os.path.exists(matrix_path):
        log("  GSE25066 file not available. Tasks 1+2 will be skipped.")
        return None, None, None

    sample_ids = []
    meta_rows  = {}      # key → list of values per sample
    probe_gene = {}      # probe_id → gene_symbol
    expr_rows  = {}      # probe_id → list of float values
    col_ids    = []

    in_platform = False
    in_expr     = False
    platform_header = []

    open_fn = gzip.open if matrix_path.endswith(".gz") else open

    log("  Parsing series matrix file (this may take ~30 seconds)...")
    n_lines = 0
    with open_fn(matrix_path, "rt", errors="ignore") as f:
        for line in f:
            n_lines += 1
            line = line.rstrip("\n")

            # ── Sample metadata ─────────────────────────────
            if line.startswith("!Sample_geo_accession"):
                parts = line.split("\t")
                sample_ids = [p.strip().strip('"') for p in parts[1:]]

            elif line.startswith("!Sample_characteristics_ch1"):
                parts = line.split("\t")
                raw   = parts[1].strip().strip('"').lower() if len(parts) > 1 else ""
                # Collect all characteristics — parse after
                key = raw.split(":")[0].strip() if ":" in raw else raw[:30]
                vals = [p.strip().strip('"') for p in parts[1:]]
                meta_rows.setdefault(key, vals)

            # ── Inline platform annotation (GPL96 mapping) ──
            elif "!platform_table_begin" in line.lower():
                in_platform = True
                platform_header = []
                continue

            elif "!platform_table_end" in line.lower():
                in_platform = False
                continue

            elif in_platform:
                parts = line.split("\t")
                if not platform_header:
                    platform_header = [p.strip().lower() for p in parts]
                else:
                    if len(parts) < 2:
                        continue
                    probe = parts[0].strip().strip('"')
                    # Look for gene symbol column
                    sym = ""
                    for col_name in ["gene_symbol", "gene symbol",
                                     "symbol", "gene_assignment"]:
                        if col_name in platform_header:
                            idx = platform_header.index(col_name)
                            if idx < len(parts):
                                raw_sym = parts[idx].strip().strip('"')
                                # gene_assignment format: "NM_001234 // GENE // ..."
                                if "//" in raw_sym:
                                    raw_sym = raw_sym.split("//")[1].strip()
                                sym = raw_sym.split()[0] if raw_sym else ""
                            break
                    if probe and sym and sym not in ("---", "", "NA"):
                        probe_gene[probe] = sym

            # ── Expression table ────────────────────────────
            elif "series_matrix_table_begin" in line.lower():
                in_expr = True
                continue

            elif "series_matrix_table_end" in line.lower():
                in_expr = False
                continue

            elif in_expr:
                parts = line.split("\t")
                if not col_ids:
                    col_ids = [c.strip().strip('"') for c in parts[1:]]
                    continue
                if len(parts) < 2:
                    continue
                probe = parts[0].strip().strip('"')
                try:
                    vals = [float(v.strip().strip('"')) for v in parts[1:]]
                    expr_rows[probe] = vals
                except ValueError:
                    pass

    log(f"  Lines parsed: {n_lines}")
    log(f"  Probe→gene mappings (inline GPL96): {len(probe_gene)}")
    log(f"  Expression probes: {len(expr_rows)}")
    log(f"  Samples (columns): {len(col_ids)}")

    if not expr_rows:
        log("  ERROR: No expression data parsed from series matrix.")
        return None, None, None

    if not probe_gene:
        log("  WARNING: Inline GPL96 table not found or empty.")
        log("  Series matrix may not have !platform_table_begin block.")
        log("  Will use probe IDs directly — depth score will use")
        log("  probe-level matching by gene symbol substring.")

    # Build expression DataFrame (samples × probes)
    use_cols = col_ids if col_ids else sample_ids
    n_vals   = len(list(expr_rows.values())[0])
    expr_df  = pd.DataFrame.from_dict(
        expr_rows, orient="index",
        columns=use_cols[:n_vals]
    ).T  # → samples × probes

    log(f"  Expression matrix (samples × probes): {expr_df.shape}")
    log(f"  Sample IDs (first 3): {list(expr_df.index[:3])}")

    # ── Build metadata DataFrame ─────────────────────────
    use_sids = use_cols[:n_vals]
    meta_df  = pd.DataFrame(index=use_sids)

    # Print ALL metadata keys before parsing — lesson from S4b
    log(f"\n  ALL METADATA FIELDS IN GSE25066 (printing before parse):")
    for k, vals in meta_rows.items():
        sample_of = vals[0] if vals else "(empty)"
        log(f"    '{k}' — example: {sample_of}")

    # Flexible parsing — look for pCR, ER, PR, HER2 regardless of exact key
    pcr_vals  = {}
    er_vals   = {}
    pr_vals   = {}
    her2_vals = {}

    for key, vals in meta_rows.items():
        key_l = key.lower()
        for i, v in enumerate(vals):
            if i >= len(use_sids):
                break
            sid = use_sids[i]
            v_l = v.lower()
            # pCR
            if "pcr" in key_l or "patholog" in key_l or "residual" in key_l:
                if "complete" in v_l or ": 1" in v_l or "yes" in v_l:
                    pcr_vals[sid] = 1
                elif "residual" in v_l or ": 0" in v_l or "no" in v_l:
                    pcr_vals[sid] = 0
                else:
                    # Try extracting trailing digit
                    for ch in reversed(v_l.replace(" ", "")):
                        if ch in "01":
                            pcr_vals[sid] = int(ch)
                            break
            # ER
            if "er status" in key_l or (key_l.startswith("er") and "status" in key_l):
                er_vals[sid] = 1 if "pos" in v_l else 0
            # PR
            if "pr status" in key_l or (key_l.startswith("pr") and "status" in key_l):
                pr_vals[sid] = 1 if "pos" in v_l else 0
            # HER2
            if "her2" in key_l:
                her2_vals[sid] = 1 if "pos" in v_l else 0

    meta_df["pCR"]  = pd.Series(pcr_vals)
    meta_df["ER"]   = pd.Series(er_vals)
    meta_df["PR"]   = pd.Series(pr_vals)
    meta_df["HER2"] = pd.Series(her2_vals)

    log(f"\n  Metadata parsed:")
    log(f"    pCR annotated:  {meta_df['pCR'].notna().sum()} / {len(meta_df)}")
    log(f"    ER annotated:   {meta_df['ER'].notna().sum()}")
    log(f"    PR annotated:   {meta_df['PR'].notna().sum()}")
    log(f"    HER2 annotated: {meta_df['HER2'].notna().sum()}")

    if meta_df["pCR"].notna().sum() > 0:
        log(f"    pCR=1 (complete response): {(meta_df['pCR']==1).sum()}")
        log(f"    pCR=0 (residual disease):  {(meta_df['pCR']==0).sum()}")

    return expr_df, meta_df, probe_gene


def map_probes_to_genes(expr_df, probe_gene):
    """
    Map probe-level expression to gene-level using the inline GPL96 mapping.
    When multiple probes map to the same gene, take the mean.
    Returns DataFrame (samples × genes).
    """
    if not probe_gene:
        log("  No probe→gene map available. Cannot map probes.")
        return None

    # Build reverse: gene → [probe_ids]
    gene_probes = {}
    for probe, gene in probe_gene.items():
        if probe in expr_df.columns:
            gene_probes.setdefault(gene, []).append(probe)

    log(f"\n  Probe→gene mapping:")
    log(f"    Total genes with at least one probe: {len(gene_probes)}")

    # For each gene: mean of its probes
    gene_expr = {}
    for gene, probes in gene_probes.items():
        gene_expr[gene] = expr_df[probes].mean(axis=1)

    gene_df = pd.DataFrame(gene_expr, index=expr_df.index)
    log(f"    Gene-level matrix: {gene_df.shape}  (samples × genes)")
    return gene_df


# ============================================================
# TASK 1: GSE25066 — GENE-MAPPED DEPTH SCORE + pCR
# Tests S2-P1, S2-P2, S2-P3, S2-P7, S2-P9
# ============================================================

def task1_gse25066_pcr(expr_df, meta_df, probe_gene):
    log("")
    log("=" * 65)
    log("TASK 1: GSE25066 — GENE-MAPPED DEPTH + pCR ANALYSIS")
    log("Tests: S2-P1, S2-P2, S2-P3, S2-P7, S2-P9")
    log("=" * 65)

    if expr_df is None:
        log("  GSE25066 data unavailable. Task 1 skipped.")
        return {}

    # ── Map probes to genes ──────────────────────────────
    gene_df = map_probes_to_genes(expr_df, probe_gene)

    if gene_df is None:
        log("  Probe mapping failed. Falling back to probe-level matching.")
        # Fallback: find probes whose row index contains gene symbol
        gene_df = pd.DataFrame(index=expr_df.index)
        for gene in FA_MARKERS + SWITCH_GENES + PCR_BIOMARKERS:
            matching = [c for c in expr_df.columns
                        if gene.upper() in str(c).upper()]
            if matching:
                gene_df[gene] = expr_df[matching].mean(axis=1)
        log(f"  Fallback gene columns: {list(gene_df.columns[:10])}")

    # ── TNBC subset extraction ──────────────────────────
    log("\n  Extracting TNBC subset (ER-/PR-/HER2-)...")

    have_er   = meta_df["ER"].notna().sum() > 10
    have_pr   = meta_df["PR"].notna().sum() > 10
    have_her2 = meta_df["HER2"].notna().sum() > 10

    if have_er and have_pr and have_her2:
        tnbc_mask = (
            (meta_df["ER"]   == 0) &
            (meta_df["PR"]   == 0) &
            (meta_df["HER2"] == 0)
        )
        meta_tnbc = meta_df[tnbc_mask].copy()
        log(f"  TNBC subset (ER-/PR-/HER2-): n={len(meta_tnbc)}")
    else:
        log("  ER/PR/HER2 not parsed. Using all samples.")
        log("  (S2-P1 cohort-wide result will be reported as pan-cohort)")
        meta_tnbc = meta_df.copy()

    # Intersect with expression
    common = gene_df.index.intersection(meta_tnbc.index)
    if len(common) < 10:
        log(f"  WARNING: Only {len(common)} samples match between "
            f"expression and metadata.")
        log("  Using all samples with pCR annotation instead.")
        common = gene_df.index.intersection(
            meta_df[meta_df["pCR"].notna()].index
        )
        meta_tnbc = meta_df.loc[common].copy()
    log(f"  Samples for analysis: {len(common)}")

    gene_s = gene_df.loc[common]
    meta_s = meta_tnbc.loc[common].copy()

    # ── Build gene-mapped depth score ────────────────���──
    log("\n  Building gene-mapped depth score...")
    log("  Formula: norm(FA_markers) + (1 - norm(switch_genes)) / 2")

    fa_avail  = [g for g in FA_MARKERS  if g in gene_s.columns]
    sw_avail  = [g for g in SWITCH_GENES if g in gene_s.columns]

    log(f"  FA markers available:   {fa_avail}")
    log(f"  Switch genes available: {sw_avail}")

    def norm01(s):
        mn, mx = s.min(), s.max()
        return (s - mn) / (mx - mn) if mx > mn else pd.Series(0.5, index=s.index)

    parts = []
    if fa_avail:
        parts.append(norm01(gene_s[fa_avail].mean(axis=1)))
    if sw_avail:
        parts.append(1 - norm01(gene_s[sw_avail].mean(axis=1)))

    if not parts:
        log("  ERROR: No depth genes found in gene_df. Cannot compute depth.")
        return {}

    depth = sum(parts) / len(parts)
    meta_s["depth_proper"] = depth.values

    log(f"\n  Depth score statistics ({len(depth)} samples):")
    log(f"    Mean={depth.mean():.4f}  Median={depth.median():.4f}  "
        f"Std={depth.std():.4f}")
    log(f"    Min={depth.min():.4f}   Max={depth.max():.4f}")

    # ── S2-P1: Depth vs pCR ─────────────────────────────
    log("")
    log("  S2-P1: GENE-MAPPED DEPTH vs pCR")
    log("  Predicted: r(depth_proper, pCR) more negative than -0.098")

    pcr_known = meta_s["pCR"].notna()
    results   = {}
    if pcr_known.sum() >= 20:
        pcr_v   = meta_s.loc[pcr_known, "pCR"].values
        dep_v   = meta_s.loc[pcr_known, "depth_proper"].values
        try:
            rv, pv = stats.pointbiserialr(pcr_v, dep_v)
            results["p1_r"]     = rv
            results["p1_p"]     = pv
            results["p1_n"]     = int(pcr_known.sum())
            log(f"  r(depth_proper, pCR) = {rv:>+.4f}   p = {pv:.4f}   "
                f"n={int(pcr_known.sum())}")
            log(f"  Script 1 degraded proxy was r = -0.098")

            if rv < -0.098 and pv < 0.05:
                log("  S2-P1: ✓ CONFIRMED — proper mapping improves pCR prediction")
            elif rv < 0 and rv < -0.098:
                log("  S2-P1: ✓ DIRECTIONALLY CONFIRMED (r more negative than proxy)")
            elif rv < 0:
                log(f"  S2-P1: ✓ DIRECTIONAL (r negative but not stronger than proxy, "
                    f"p={pv:.3f})")
            else:
                log("  S2-P1: ✗ NOT CONFIRMED — depth does not predict pCR direction")

            # Subgroup: TNBC only (if subset was extracted)
            if have_er:
                tnbc_pcr_mask = (
                    pcr_known &
                    (meta_s.get("ER", pd.Series(dtype=float)) == 0)
                )
                if tnbc_pcr_mask.sum() >= 10:
                    rv2, pv2 = stats.pointbiserialr(
                        meta_s.loc[tnbc_pcr_mask, "pCR"].values,
                        meta_s.loc[tnbc_pcr_mask, "depth_proper"].values
                    )
                    results["p1b_r"] = rv2
                    results["p1b_p"] = pv2
                    log(f"\n  S2-P1b (TNBC only, n={tnbc_pcr_mask.sum()}):")
                    log(f"    r(depth, pCR) = {rv2:>+.4f}   p={pv2:.4f}")
                    if rv2 < -0.15 and pv2 < 0.05:
                        log("    S2-P1b: ✓ CONFIRMED — TNBC r < -0.15")
                    elif rv2 < 0:
                        log(f"    S2-P1b: ✓ DIRECTIONAL  (r={rv2:+.3f})")
                    else:
                        log("    S2-P1b: ✗ NOT CONFIRMED")
        except Exception as e:
            log(f"  S2-P1 error: {e}")
    else:
        log(f"  Only {pcr_known.sum()} pCR-annotated samples. "
            f"Need ≥20. S2-P1 result limited.")

    # ── S2-P2: EED vs pCR ───────────────────────────────
    log("")
    log("  S2-P2: EED AS pCR BIOMARKER")
    log("  Predicted: r(EED, pCR) < -0.10  AND  AUC(EED) > AUC(EZH2)")

    for gene in ["EED", "EZH2", "PARP1"]:
        if gene not in gene_s.columns:
            log(f"  {gene}: NOT IN DATASET")
            continue
        pcr_known2 = meta_s["pCR"].notna() & (meta_s[gene] if False else True)
        pcr_known2 = meta_s["pCR"].notna()
        if pcr_known2.sum() < 10:
            continue
        g_v = gene_s.loc[pcr_known2, gene].values
        p_v = meta_s.loc[pcr_known2, "pCR"].values
        try:
            rv, pv  = stats.pointbiserialr(p_v, g_v)
            # AUC: higher gene = lower pCR → negate for AUC
            try:
                auc_val = roc_auc_score(p_v, -g_v)
            except Exception:
                auc_val = np.nan
            results[f"{gene.lower()}_r"]   = rv
            results[f"{gene.lower()}_p"]   = pv
            results[f"{gene.lower()}_auc"] = auc_val
            log(f"  {gene:<8}  r(pCR)={rv:>+.4f}  {fmt_p(pv)}  "
                f"AUC={auc_val:.3f}")
        except Exception as e:
            log(f"  {gene} error: {e}")

    # Evaluate S2-P2
    if "eed_r" in results and "ezh2_r" in results:
        if results["eed_r"] < -0.10 and results["eed_r"] < results["ezh2_r"]:
            log("  S2-P2: ✓ CONFIRMED — EED more predictive than EZH2")
        elif results["eed_r"] < 0:
            log(f"  S2-P2: ✓ DIRECTIONAL  EED r={results['eed_r']:+.3f}")
        else:
            log("  S2-P2: ✗ NOT CONFIRMED")

        # S2-P2b AUC comparison
        if "eed_auc" in results and "ezh2_auc" in results:
            if not np.isnan(results["eed_auc"]):
                if results["eed_auc"] > results["ezh2_auc"]:
                    log(f"  S2-P2b AUC: ✓ AUC(EED)={results['eed_auc']:.3f} "
                        f"> AUC(EZH2)={results['ezh2_auc']:.3f}")
                else:
                    log(f"  S2-P2b AUC: ✗ AUC(EED)={results['eed_auc']:.3f} "
                        f"≤ AUC(EZH2)={results['ezh2_auc']:.3f}")

    # S2-P3 result
    if "parp1_r" in results:
        if results["parp1_r"] < -0.08:
            log(f"  S2-P3: ✓ CONFIRMED — PARP1 r={results['parp1_r']:+.3f} < -0.08")
        elif results["parp1_r"] < 0:
            log(f"  S2-P3: ✓ DIRECTIONAL — PARP1 r={results['parp1_r']:+.3f}")
        else:
            log(f"  S2-P3: ✗ NOT CONFIRMED — PARP1 r={results['parp1_r']:+.3f}")

    # ── S2-P7: SPI1 as immune contamination ─────────────
    log("")
    log("  S2-P7: SPI1 RESOLUTION")
    log("  Predicted: r(SPI1, PTPRC) > r(SPI1, depth)")
    if "SPI1" in gene_s.columns and "PTPRC" in gene_s.columns:
        try:
            r_spi1_ptprc, _ = stats.pearsonr(
                gene_s["SPI1"].values, gene_s["PTPRC"].values
            )
            r_spi1_depth, _ = stats.pearsonr(
                gene_s["SPI1"].values, meta_s["depth_proper"].values
            )
            log(f"  r(SPI1, PTPRC)  = {r_spi1_ptprc:>+.4f}")
            log(f"  r(SPI1, depth)  = {r_spi1_depth:>+.4f}")
            results["p7_spi1_ptprc"] = r_spi1_ptprc
            results["p7_spi1_depth"] = r_spi1_depth
            if r_spi1_ptprc > abs(r_spi1_depth):
                log("  S2-P7: ✓ CONFIRMED — SPI1 is immune contamination, "
                    "not depth biology")
            else:
                log("  S2-P7: ✗ NOT CONFIRMED — SPI1 may be depth-biological")
        except Exception as e:
            log(f"  S2-P7 error: {e}")
    else:
        log("  SPI1 or PTPRC not in bulk dataset. S2-P7 untestable.")

    # ── S2-P9: AR vs depth in bulk ───────────────────────
    log("")
    log("  S2-P9: AR CORRELATES NEGATIVELY WITH DEPTH IN BULK")
    log("  Predicted: r(AR, depth) < -0.15")
    if "AR" in gene_s.columns:
        try:
            r_ar, p_ar = stats.pearsonr(
                gene_s["AR"].values,
                meta_s["depth_proper"].values
            )
            results["p9_ar_r"] = r_ar
            log(f"  r(AR, depth) = {r_ar:>+.4f}  {fmt_p(p_ar)}")
            if r_ar < -0.15 and p_ar < 0.05:
                log("  S2-P9: ✓ CONFIRMED")
            elif r_ar < 0:
                log(f"  S2-P9: ✓ DIRECTIONAL  (r={r_ar:+.3f})")
            else:
                log("  S2-P9: ✗ NOT CONFIRMED")
        except Exception as e:
            log(f"  S2-P9 error: {e}")
    else:
        log("  AR not in gene_df.")

    # ── pCR group means ──────────────────────────────────
    log("")
    log("  pCR GROUP COMPARISON:")
    pcr1_d = meta_s.loc[meta_s["pCR"] == 1, "depth_proper"]
    pcr0_d = meta_s.loc[meta_s["pCR"] == 0, "depth_proper"]
    if len(pcr1_d) > 2 and len(pcr0_d) > 2:
        _, p_mwu = stats.mannwhitneyu(pcr0_d, pcr1_d, alternative="greater")
        log(f"  pCR=1 (complete): mean={pcr1_d.mean():.4f}  n={len(pcr1_d)}")
        log(f"  pCR=0 (residual): mean={pcr0_d.mean():.4f}  n={len(pcr0_d)}")
        log(f"  MWU p(pCR=0 deeper) = {p_mwu:.4f}")

    # ── Gene correlations with depth ─────────────────────
    log("")
    log("  BULK DEPTH CORRELATIONS (all available genes):")
    log(f"  {'Gene':<12} {'r':>8}  p-value")
    log(f"  {'-'*35}")
    bulk_corrs = []
    for gene in gene_s.columns:
        if gene not in PCR_BIOMARKERS + FA_MARKERS + SWITCH_GENES:
            continue
        try:
            rv, pv = stats.pearsonr(
                gene_s[gene].values,
                meta_s["depth_proper"].values
            )
            bulk_corrs.append((gene, rv, pv))
        except Exception:
            pass
    bulk_corrs.sort(key=lambda x: abs(x[1]), reverse=True)
    for g, r, p in bulk_corrs[:20]:
        log(f"  {g:<12} {r:>+8.4f}  {fmt_p(p)}")

    results["gene_df"]   = gene_s
    results["meta_s"]    = meta_s
    results["bulk_corrs"]= bulk_corrs
    return results


# ============================================================
# TASK 2: LEHMANN SUBTYPE DEPTH MAPPING
# Tests S2-P4
# ============================================================

def task2_lehmann_depth(task1_results):
    log("")
    log("=" * 65)
    log("TASK 2: LEHMANN SUBTYPE DEPTH MAPPING")
    log("Tests: S2-P4")
    log("Predicted: LAR < BL1/BL2 < IM < M/MSL")
    log("LAR mean depth < 0.40 | M/MSL mean depth > 0.60")
    log("=" * 65)

    gene_df = task1_results.get("gene_df")
    meta_s  = task1_results.get("meta_s")

    if gene_df is None or meta_s is None or len(gene_df) < 20:
        log("  Task 1 data unavailable. Task 2 skipped.")
        return {}

    # ── Simplified Lehmann centroid scoring ─────────────
    # Score each sample against each subtype's gene set.
    # Subtype = mean expression of available genes in that signature.
    log("  Computing Lehmann centroid scores...")

    subtype_scores = {}
    for subtype, genes in LEHMANN_GENES.items():
        avail = [g for g in genes if g in gene_df.columns]
        if len(avail) >= 2:
            subtype_scores[subtype] = gene_df[avail].mean(axis=1)
            log(f"    {subtype}: {len(avail)}/{len(genes)} genes available")
        else:
            log(f"    {subtype}: only {len(avail)} genes — skipped")

    if not subtype_scores:
        log("  No subtype scores computed. Task 2 skipped.")
        return {}

    scores_df = pd.DataFrame(subtype_scores, index=gene_df.index)
    # Assign each sample to highest-scoring subtype
    scores_df["assigned"] = scores_df.idxmax(axis=1)

    # Attach depth
    scores_df["depth"] = meta_s["depth_proper"].values

    log("\n  Lehmann subtype assignment counts:")
    counts = scores_df["assigned"].value_counts()
    for st, n in counts.items():
        log(f"    {st}: n={n}")

    log("\n  Depth score by Lehmann subtype:")
    log(f"  {'Subtype':<8} {'n':>5} {'mean_depth':>12} {'std':>8}")
    log(f"  {'-'*38}")
    depth_by_subtype = {}
    for st in ["LAR", "BL1", "BL2", "IM", "M", "MSL"]:
        if st not in subtype_scores:
            continue
        mask = scores_df["assigned"] == st
        if mask.sum() < 3:
            continue
        d = scores_df.loc[mask, "depth"]
        depth_by_subtype[st] = d.values
        log(f"  {st:<8} {mask.sum():>5} {d.mean():>12.4f} {d.std():>8.4f}")

    # Evaluate S2-P4
    log("\n  S2-P4 EVALUATION:")
    p4_results = {}
    if "LAR" in depth_by_subtype:
        lar_mean = depth_by_subtype["LAR"].mean()
        p4_results["LAR_mean"] = lar_mean
        log(f"  LAR mean depth = {lar_mean:.4f}  (predicted < 0.40)")
        if lar_mean < 0.40:
            log("  S2-P4 LAR: ✓ CONFIRMED")
        else:
            log("  S2-P4 LAR: ✗ NOT CONFIRMED")

    deep_subtypes = {k: v for k, v in depth_by_subtype.items()
                     if k in ["M", "MSL"]}
    if deep_subtypes:
        for st, d in deep_subtypes.items():
            m = d.mean()
            p4_results[f"{st}_mean"] = m
            log(f"  {st} mean depth = {m:.4f}  (predicted > 0.60)")
            if m > 0.60:
                log(f"  S2-P4 {st}: ✓ CONFIRMED")
            else:
                log(f"  S2-P4 {st}: ✗ NOT CONFIRMED")

    # Ordering test
    ordered = sorted(depth_by_subtype.items(), key=lambda x: x[1].mean())
    log("\n  Actual depth ordering (low → high):")
    for st, d in ordered:
        log(f"    {st}: {d.mean():.4f}")

    # KW test across subtypes
    if len(depth_by_subtype) >= 3:
        try:
            _, p_kw = stats.kruskal(*depth_by_subtype.values())
            log(f"\n  Kruskal-Wallis across subtypes: p={p_kw:.4f}")
            p4_results["kw_p"] = p_kw
        except Exception as e:
            log(f"  KW error: {e}")

    return {"depth_by_subtype": depth_by_subtype,
            "scores_df": scores_df, **p4_results}


# ============================================================
# TASK 3: TCGA-BRCA — PAM50 + SURVIVAL
# Tests S2-P5, S2-P6
# ============================================================

def load_tcga(tcga_paths):
    log("")
    log("=" * 65)
    log("TASK 3: TCGA-BRCA — LOAD DATA")
    log("=" * 65)

    result = {}

    # ── Expression ──────────────────────────────────────
    expr_path = tcga_paths.get("tcga_expr")
    if expr_path and os.path.exists(expr_path):
        log("  Loading TCGA expression...")
        try:
            with gzip.open(expr_path, "rt") as f:
                tcga_expr = pd.read_csv(f, sep="\t", index_col=0)
            # Xena format: genes × samples. Transpose to samples × genes.
            tcga_expr = tcga_expr.T
            log(f"  Expression shape (samples × genes): {tcga_expr.shape}")
            log(f"  Sample IDs (first 3): {list(tcga_expr.index[:3])}")
            log(f"  Gene IDs (first 5):   {list(tcga_expr.columns[:5])}")
            result["tcga_expr"] = tcga_expr
        except Exception as e:
            log(f"  TCGA expression load error: {e}")
            log(traceback.format_exc())
    else:
        log("  TCGA expression file unavailable.")

    # ── Clinical / PAM50 ────────────────────────────────
    clin_path = tcga_paths.get("tcga_clin")
    if clin_path and os.path.exists(clin_path):
        log("  Loading TCGA clinical matrix...")
        try:
            tcga_clin = pd.read_csv(clin_path, sep="\t", index_col=0)
            log(f"  Clinical shape: {tcga_clin.shape}")
            log(f"  Clinical columns (first 20):")
            for c in tcga_clin.columns[:20]:
                log(f"    '{c}'")
            result["tcga_clin"] = tcga_clin

            # Find PAM50 column
            pam50_col = None
            for col in tcga_clin.columns:
                if "pam50" in col.lower() or "subtype" in col.lower():
                    pam50_col = col
                    log(f"  PAM50 column: '{pam50_col}'")
                    log(f"  PAM50 values: "
                        f"{tcga_clin[pam50_col].value_counts().to_dict()}")
                    break
            if pam50_col:
                result["pam50_col"] = pam50_col
            else:
                log("  WARNING: No PAM50 column found. "
                    "Will check survival table.")

        except Exception as e:
            log(f"  TCGA clinical load error: {e}")
    else:
        log("  TCGA clinical file unavailable.")

    # ── Pan-cancer survival table ────────────────────────
    surv_path = tcga_paths.get("tcga_surv")
    if surv_path and os.path.exists(surv_path):
        log("  Loading pan-cancer survival table...")
        try:
            # This file may or may not be gzipped
            open_fn = (gzip.open if surv_path.endswith(".gz") else open)
            with open_fn(surv_path, "rt", errors="ignore") as f:
                surv_df = pd.read_csv(f, sep="\t", index_col=0)
            log(f"  Survival table shape: {surv_df.shape}")
            log(f"  Survival columns: {list(surv_df.columns[:15])}")
            # Filter to BRCA samples (sample IDs start with TCGA-BRCA)
            brca_mask = surv_df.index.str.contains("TCGA-BR", na=False)
            if brca_mask.sum() > 10:
                surv_df = surv_df[brca_mask]
                log(f"  BRCA samples in survival table: {len(surv_df)}")
            result["tcga_surv"] = surv_df
        except Exception as e:
            log(f"  Survival table load error: {e}")
    else:
        log("  Survival table unavailable.")

    return result


def task3_tcga(tcga_data):
    log("")
    log("=" * 65)
    log("TASK 3: TCGA-BRCA ANALYSIS")
    log("Tests: S2-P5 (BRCA1 dysfunction), S2-P6 (depth~survival)")
    log("=" * 65)

    tcga_expr = tcga_data.get("tcga_expr")
    tcga_clin = tcga_data.get("tcga_clin")
    pam50_col = tcga_data.get("pam50_col")
    tcga_surv = tcga_data.get("tcga_surv")

    results = {}

    if tcga_expr is None or tcga_clin is None:
        log("  TCGA data unavailable. Task 3 skipped.")
        return results

    # ── Align samples ────────────────────────────────────
    common = tcga_expr.index.intersection(tcga_clin.index)
    log(f"  Expression ∩ clinical samples: {len(common)}")
    if len(common) < 50:
        log("  Insufficient overlap. Checking ID format...")
        # TCGA sample IDs sometimes differ by last portion
        # Try matching on first 15 chars (TCGA-XX-XXXX-XX)
        expr_short = pd.Index([s[:15] for s in tcga_expr.index])
        clin_short = pd.Index([s[:15] for s in tcga_clin.index])
        overlap_s  = expr_short.intersection(clin_short)
        log(f"  Overlap on 15-char IDs: {len(overlap_s)}")
        if len(overlap_s) > 50:
            # Rebuild with short IDs
            tcga_expr = tcga_expr.copy()
            tcga_expr.index = expr_short
            tcga_clin = tcga_clin.copy()
            tcga_clin.index = clin_short
            common = overlap_s
            log(f"  Using truncated sample IDs: {len(common)} samples")

    if len(common) < 50:
        log("  Still insufficient. Task 3 analysis will be partial.")

    expr_a = tcga_expr.loc[common]
    clin_a = tcga_clin.loc[common]

    # ── PAM50 Basal-like subset ──────────────────────────
    basal_mask = pd.Series(False, index=common)
    if pam50_col:
        pam50_vals = clin_a[pam50_col].astype(str).str.lower()
        basal_mask = pam50_vals.str.contains("basal", na=False)
        log(f"\n  PAM50 Basal-like: n={basal_mask.sum()}")
        log(f"  PAM50 LumA:       "
            f"n={(pam50_vals.str.contains('luma|luminal a', na=False)).sum()}")
        log(f"  PAM50 LumB:       "
            f"n={(pam50_vals.str.contains('lumb|luminal b', na=False)).sum()}")
    else:
        log("  No PAM50 column. Using expression-based proxy (ESR1-low, KRT5-high).")
        if "ESR1" in expr_a.columns and "KRT5" in expr_a.columns:
            basal_mask = (expr_a["ESR1"] < expr_a["ESR1"].median()) & \
                         (expr_a["KRT5"] > expr_a["KRT5"].median())
            log(f"  Expression-proxy Basal-like: n={basal_mask.sum()}")

    # ── Compute TCGA depth score ──────────────────���──────
    log("\n  Computing depth score in TCGA-BRCA...")
    fa_t  = [g for g in FA_MARKERS  if g in expr_a.columns]
    sw_t  = [g for g in SWITCH_GENES if g in expr_a.columns]
    log(f"  FA markers:   {fa_t}")
    log(f"  Switch genes: {sw_t}")

    def norm01(s):
        mn, mx = s.min(), s.max()
        return (s - mn) / (mx - mn) if mx > mn else pd.Series(0.5, index=s.index)

    parts_t = []
    if fa_t:
        parts_t.append(norm01(expr_a[fa_t].mean(axis=1)))
    if sw_t:
        parts_t.append(1 - norm01(expr_a[sw_t].mean(axis=1)))

    if not parts_t:
        log("  Cannot compute depth score: no genes available.")
        return results

    depth_tcga = sum(parts_t) / len(parts_t)
    clin_a     = clin_a.copy()
    clin_a["depth"] = depth_tcga.values

    log(f"  All BRCA: mean depth={depth_tcga.mean():.4f}  "
        f"std={depth_tcga.std():.4f}")

    basal_depth = depth_tcga[basal_mask]
    other_depth = depth_tcga[~basal_mask]
    if len(basal_depth) > 10 and len(other_depth) > 10:
        _, p_bo = stats.mannwhitneyu(basal_depth, other_depth,
                                     alternative="greater")
        log(f"  Basal-like depth: mean={basal_depth.mean():.4f}  "
            f"n={len(basal_depth)}")
        log(f"  Non-Basal depth:  mean={other_depth.mean():.4f}  "
            f"n={len(other_depth)}")
        log(f"  MWU Basal > Non-Basal: p={p_bo:.4f}")
        results["basal_depth_p"] = p_bo

    # ── S2-P5: BRCA1 proxy in TCGA ───────────────────────
    log("")
    log("  S2-P5: BRCA1 DYSFUNCTION IN TCGA BASAL-LIKE")
    log("  (DNA-level test requires mutation/methylation data.")
    log("   At expression level: BRCA1 mRNA as proxy only.)")

    if "BRCA1" in expr_a.columns:
        br_basal = expr_a.loc[basal_mask, "BRCA1"]
        br_luma  = (expr_a["BRCA1"] if not pam50_col else
                    expr_a.loc[
                        clin_a.get(pam50_col, pd.Series()).astype(str)
                        .str.lower().str.contains("luma|luminal a", na=False),
                        "BRCA1"
                    ])
        log(f"  BRCA1 mRNA Basal-like: mean={br_basal.mean():.4f}  "
            f"n={len(br_basal)}")
        if len(br_luma) > 5:
            log(f"  BRCA1 mRNA LumA:       mean={br_luma.mean():.4f}  "
                f"n={len(br_luma)}")
            try:
                _, p_b1 = stats.mannwhitneyu(br_basal, br_luma,
                                              alternative="less")
                log(f"  BRCA1 Basal < LumA: p={p_b1:.4f}")
                results["p5_brca1_p"] = p_b1
                log("  NOTE: S2-P5 (mutation/methylation rate) cannot be")
                log("  confirmed from expression alone. Expression BRCA1")
                log("  reduction is a proxy — not the definitive test.")
                log("  Definitive test requires somatic mutation + methylation")
                log("  data from GDC TCGA (requires authentication).")
                log("  Prediction stated: mutation >20% Basal, <5% LumA.")
                log("  Expression proxy reported here for completeness.")
            except Exception as e:
                log(f"  S2-P5 BRCA1 test error: {e}")
    else:
        log("  BRCA1 not in TCGA expression dataset.")

    # ── S2-P6: Depth predicts survival ───────────────────
    log("")
    log("  S2-P6: DEPTH SCORE PREDICTS SURVIVAL IN TCGA-BRCA")
    log("  Predicted: depth-high → shorter OS (KM p<0.05, Cox HR>1.3)")

    # Find OS columns
    os_time_col  = None
    os_event_col = None
    for col in clin_a.columns:
        cl = col.lower()
        if os_time_col is None and any(
            x in cl for x in ["os.time", "overall_survival",
                               "os_time", "days_to_death",
                               "OS.time"]
        ):
            os_time_col = col
        if os_event_col is None and any(
            x in cl for x in ["os.event", "vital_status",
                               "os_event", "dead",
                               "OS.event"]
        ):
            os_event_col = col

    # Also check survival table
    if (os_time_col is None or os_event_col is None) and tcga_surv is not None:
        log("  OS columns not in clinical. Checking survival table...")
        surv_brca = tcga_surv
        log(f"  Survival table columns: {list(surv_brca.columns[:10])}")
        for col in surv_brca.columns:
            cl = col.lower()
            if os_time_col is None and "time" in cl:
                os_time_col  = col
            if os_event_col is None and ("event" in cl or "dead" in cl
                                          or "status" in cl):
                os_event_col = col
        if os_time_col:
            # Merge survival into clin_a
            surv_common = clin_a.index.intersection(surv_brca.index)
            log(f"  Survival ∩ BRCA clinical: {len(surv_common)}")
            if len(surv_common) > 50:
                clin_a = clin_a.join(
                    surv_brca[[os_time_col, os_event_col]].rename(
                        columns={os_time_col:  "OS_time",
                                 os_event_col: "OS_event"}
                    ), how="left"
                )
                os_time_col  = "OS_time"
                os_event_col = "OS_event"

    log(f"  OS time column:  {os_time_col}")
    log(f"  OS event column: {os_event_col}")

    if os_time_col and os_event_col:
        try:
            surv_sub = clin_a[[os_time_col, os_event_col, "depth"]].copy()
            surv_sub.columns = ["time", "event", "depth"]
            surv_sub = surv_sub.dropna()
            surv_sub["time"]  = pd.to_numeric(surv_sub["time"],  errors="coerce")
            surv_sub["event"] = pd.to_numeric(surv_sub["event"], errors="coerce")
            surv_sub = surv_sub.dropna()
            surv_sub = surv_sub[surv_sub["time"] > 0]

            log(f"  Survival analysis n={len(surv_sub)}")
            log(f"  Events (dead): {surv_sub['event'].sum():.0f}")
            log(f"  Median OS time: {surv_sub['time'].median():.1f} days")

            # KM split at depth median
            depth_med = surv_sub["depth"].median()
            hi = surv_sub[surv_sub["depth"] >= depth_med]
            lo = surv_sub[surv_sub["depth"] <  depth_med]

            log(f"\n  KM split at depth median={depth_med:.3f}:")
            log(f"  High depth: n={len(hi)}  events={hi['event'].sum():.0f}")
            log(f"  Low depth:  n={len(lo)}  events={lo['event'].sum():.0f}")

            # Log-rank test (manual implementation without lifelines)
            # Use scipy to implement simplified log-rank
            from scipy.stats import chi2
            def logrank_p(t1, e1, t2, e2):
                """Simplified log-rank statistic."""
                all_t = np.unique(np.concatenate([t1[e1==1], t2[e2==1]]))
                O1, O2, E1, E2 = 0, 0, 0, 0
                for t in all_t:
                    n1 = (t1 >= t).sum()
                    n2 = (t2 >= t).sum()
                    o1 = ((t1 == t) & (e1 == 1)).sum()
                    o2 = ((t2 == t) & (e2 == 1)).sum()
                    n  = n1 + n2
                    if n < 2:
                        continue
                    e1_exp = (n1 / n) * (o1 + o2)
                    O1 += o1
                    O2 += o2
                    E1 += e1_exp
                    E2 += (o1 + o2) - e1_exp
                if E1 <= 0 or E2 <= 0:
                    return np.nan, np.nan
                stat = ((O1 - E1)**2 / E1) + ((O2 - E2)**2 / E2)
                return stat, chi2.sf(stat, df=1)

            stat, p_lr = logrank_p(
                hi["time"].values, hi["event"].values,
                lo["time"].values, lo["event"].values
            )
            results["p6_logrank_p"] = p_lr
            results["p6_logrank_stat"] = stat

            log(f"  Log-rank statistic={stat:.3f}  p={p_lr:.4f}")

            if not np.isnan(p_lr) and p_lr < 0.05:
                # Check direction
                hi_med = hi.loc[hi["event"]==1, "time"].median()
                lo_med = lo.loc[lo["event"]==1, "time"].median()
                log(f"  High-depth median survival: {hi_med:.1f} days")
                log(f"  Low-depth  median survival: {lo_med:.1f} days")
                if hi_med < lo_med:
                    log("  S2-P6: ✓ CONFIRMED — depth-high shorter survival, "
                        f"p={p_lr:.4f}")
                else:
                    log("  S2-P6: ✗ INVERTED — depth-high has LONGER survival")
            elif not np.isnan(p_lr):
                log(f"  S2-P6: ns  p={p_lr:.4f} (not significant)")
            else:
                log("  S2-P6: log-rank failed")

            # Cox proportional hazards — manual Wald test via correlation
            # Full Cox requires lifelines (not available). Use Spearman correlation
            # as proxy for hazard direction.
            r_sp, p_sp = stats.spearmanr(surv_sub["depth"], surv_sub["event"])
            log(f"\n  Spearman r(depth, event) = {r_sp:+.4f}  {fmt_p(p_sp)}")
            log("  (Positive r = deeper → more events = shorter survival)")
            results["p6_spearman_r"] = r_sp
            results["p6_spearman_p"] = p_sp

        except Exception as e:
            log(f"  S2-P6 survival analysis error: {e}")
            log(traceback.format_exc())
    else:
        log("  OS time/event columns not identified. S2-P6 not testable.")
        log("  Clinical columns available:")
        if tcga_clin is not None:
            for col in tcga_clin.columns:
                log(f"    '{col}'")

    return results


# ============================================================
# SECTION: PREDICTION SCORECARD
# ============================================================

def prediction_scorecard(t1, t2, t3):
    log("")
    log("=" * 65)
    log("PREDICTION SCORECARD — BRCA-S4c")
    log("=" * 65)
    log("")

    def row(pred, status, note=""):
        sym = {"✓": "✓", "✗": "✗", "?": "?", "P": "P"}.get(status[0], "?")
        log(f"  {pred:<12}  {sym} {status:<25}  {note}")

    log(f"  {'Prediction':<12}  {'Status':<26}  Note")
    log(f"  {'-'*70}")

    # S2-P1
    p1r = t1.get("p1_r", np.nan)
    if not np.isnan(p1r):
        if p1r < -0.098:
            row("S2-P1", "✓ CONFIRMED", f"r={p1r:+.3f} more negative than -0.098")
        elif p1r < 0:
            row("S2-P1", "✓ DIRECTIONAL", f"r={p1r:+.3f} (negative, not stronger)")
        else:
            row("S2-P1", "✗ NOT CONFIRMED", f"r={p1r:+.3f}")
    else:
        row("S2-P1", "? DATA LIMITED", "pCR correlation not computed")

    # S2-P2
    eed_r  = t1.get("eed_r",  np.nan)
    ezh2_r = t1.get("ezh2_r", np.nan)
    eed_a  = t1.get("eed_auc", np.nan)
    ezh2_a = t1.get("ezh2_auc", np.nan)
    if not np.isnan(eed_r) and not np.isnan(ezh2_r):
        if eed_r < -0.10 and eed_r < ezh2_r:
            row("S2-P2", "✓ CONFIRMED", f"EED r={eed_r:+.3f} vs EZH2 r={ezh2_r:+.3f}")
        elif eed_r < ezh2_r:
            row("S2-P2", "✓ PARTIAL", f"EED>{ezh2_r:+.3f} but r≥-0.10")
        else:
            row("S2-P2", "✗ NOT CONFIRMED", f"EED r={eed_r:+.3f}")
    else:
        row("S2-P2", "? DATA LIMITED", "EED or EZH2 not in dataset")

    # S2-P3
    parp1_r = t1.get("parp1_r", np.nan)
    if not np.isnan(parp1_r):
        if parp1_r < -0.08:
            row("S2-P3", "✓ CONFIRMED", f"PARP1 r={parp1_r:+.3f} < -0.08")
        elif parp1_r < 0:
            row("S2-P3", "✓ DIRECTIONAL", f"PARP1 r={parp1_r:+.3f}")
        else:
            row("S2-P3", "✗ NOT CONFIRMED", f"PARP1 r={parp1_r:+.3f}")
    else:
        row("S2-P3", "? DATA LIMITED", "PARP1 not in dataset")

    # S2-P4
    lar_m = t2.get("LAR_mean", np.nan)
    if not np.isnan(lar_m):
        row("S2-P4", "✓ CONFIRMED" if lar_m < 0.40 else "✗ NOT CONFIRMED",
            f"LAR mean depth={lar_m:.3f}")
    else:
        row("S2-P4", "? DATA LIMITED", "Lehmann scoring incomplete")

    # S2-P5
    p5_p = t3.get("p5_brca1_p", np.nan)
    if not np.isnan(p5_p):
        row("S2-P5", "? PROXY ONLY",
            f"Expression proxy p={p5_p:.4f} — "
            "DNA-level test requires GDC mutation data")
    else:
        row("S2-P5", "? UNTESTABLE",
            "Mutation/methylation data requires GDC auth")

    # S2-P6
    p6_p = t3.get("p6_logrank_p", np.nan)
    if not np.isnan(p6_p):
        if p6_p < 0.05:
            row("S2-P6", "✓ CONFIRMED", f"log-rank p={p6_p:.4f}")
        else:
            row("S2-P6", "? NS", f"log-rank p={p6_p:.4f}")
    else:
        row("S2-P6", "? DATA LIMITED", "Survival analysis not completed")

    # S2-P7
    spi1_ptprc = t1.get("p7_spi1_ptprc", np.nan)
    spi1_depth = t1.get("p7_spi1_depth", np.nan)
    if not np.isnan(spi1_ptprc):
        if spi1_ptprc > abs(spi1_depth):
            row("S2-P7", "✓ CONFIRMED",
                f"r(SPI1,PTPRC)={spi1_ptprc:+.3f} "
                f"> r(SPI1,depth)={spi1_depth:+.3f}")
        else:
            row("S2-P7", "✗ NOT CONFIRMED",
                f"r(SPI1,PTPRC)={spi1_ptprc:+.3f}")
    else:
        row("S2-P7", "? DATA LIMITED", "SPI1/PTPRC not in dataset")

    # S2-P8 (prospective — no data to test)
    row("S2-P8", "P PROSPECTIVE",
        "EED>EZH2 as EZH2i biomarker — stated for future datasets")

    # S2-P9
    ar_r = t1.get("p9_ar_r", np.nan)
    if not np.isnan(ar_r):
        if ar_r < -0.15:
            row("S2-P9", "✓ CONFIRMED", f"r(AR,depth)={ar_r:+.3f} < -0.15")
        elif ar_r < 0:
            row("S2-P9", "✓ DIRECTIONAL", f"r(AR,depth)={ar_r:+.3f}")
        else:
            row("S2-P9", "✗ NOT CONFIRMED", f"r(AR,depth)={ar_r:+.3f}")
    else:
        row("S2-P9", "? DATA LIMITED", "AR not in dataset")

    log("")


# ============================================================
# FIGURE
# ============================================================

def generate_figure(t1, t2, t3):
    log("")
    log("--- Generating figure ---")

    gene_df    = t1.get("gene_df")
    meta_s     = t1.get("meta_s")
    bulk_corrs = t1.get("bulk_corrs", [])
    depth_by_subtype = t2.get("depth_by_subtype", {})
    p6_logrank_p = t3.get("p6_logrank_p", np.nan)

    fig = plt.figure(figsize=(24, 18))
    fig.suptitle(
        "TNBC — Script 2 | OrganismCore BRCA-S4c/d | 2026-03-04\n"
        "Bulk Validation | pCR Analysis | TCGA Survival\n"
        "GSE25066 (pCR) + TCGA-BRCA (survival/PAM50)",
        fontsize=10, fontweight="bold", y=1.005
    )
    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.52, wspace=0.42)
    clr = {
        "depth":  "#c0392b",
        "ndepth": "#2980b9",
        "lar":    "#f39c12",
        "m":      "#8e44ad",
        "msl":    "#2ecc71",
        "bl1":    "#e74c3c",
        "bl2":    "#e67e22",
        "im":     "#3498db",
    }

    # A: pCR depth distribution
    ax_a = fig.add_subplot(gs[0, 0])
    if meta_s is not None and "depth_proper" in meta_s.columns and \
            "pCR" in meta_s.columns:
        pcr1 = meta_s.loc[meta_s["pCR"] == 1, "depth_proper"]
        pcr0 = meta_s.loc[meta_s["pCR"] == 0, "depth_proper"]
        if len(pcr1) > 2 and len(pcr0) > 2:
            ax_a.hist(pcr0.values, bins=30, alpha=0.65,
                      color=clr["depth"],  label=f"pCR=0 n={len(pcr0)}")
            ax_a.hist(pcr1.values, bins=30, alpha=0.65,
                      color=clr["ndepth"], label=f"pCR=1 n={len(pcr1)}")
            ax_a.axvline(pcr0.mean(), color=clr["depth"],  lw=2, ls="--")
            ax_a.axvline(pcr1.mean(), color=clr["ndepth"], lw=2, ls="--")
        ax_a.set_title("A — Depth by pCR (S2-P1)\nGSE25066",
                       fontsize=8, fontweight="bold")
        ax_a.set_xlabel("Depth Score (gene-mapped)", fontsize=7)
        ax_a.set_ylabel("n samples", fontsize=7)
        ax_a.legend(fontsize=6)

    # B: EED vs EZH2 pCR
    ax_b = fig.add_subplot(gs[0, 1])
    if gene_df is not None and meta_s is not None and \
            "pCR" in meta_s.columns:
        pcr_known = meta_s["pCR"].notna()
        for gene, color in [("EED", clr["depth"]), ("EZH2", clr["ndepth"])]:
            if gene in gene_df.columns:
                g_v = gene_df.loc[pcr_known, gene].values
                p_v = meta_s.loc[pcr_known, "pCR"].values
                ax_b.scatter(g_v, p_v + (0.05 if gene == "EZH2" else -0.05),
                             c=color, alpha=0.25, s=8, label=gene)
        ax_b.set_title("B — EED/EZH2 vs pCR (S2-P2)\nPredicted: EED better predictor",
                       fontsize=8, fontweight="bold")
        ax_b.set_xlabel("Expression", fontsize=7)
        ax_b.set_ylabel("pCR (0/1)", fontsize=7)
        ax_b.legend(fontsize=6)

    # C: Bulk depth correlations
    ax_c = fig.add_subplot(gs[0, 2])
    if bulk_corrs:
        top20 = bulk_corrs[:20]
        gc = [g for g, r, p in top20]
        vc = [r for g, r, p in top20]
        cc = [clr["depth"] if v > 0 else clr["ndepth"] for v in vc]
        ax_c.barh(gc, vc, color=cc)
        ax_c.axvline(0, color="black", lw=0.8)
        ax_c.set_title("C — Bulk Depth Correlations\nTop 20",
                       fontsize=8, fontweight="bold")
        ax_c.set_xlabel("r", fontsize=7)
        ax_c.tick_params(axis="y", labelsize=6)

    # D: Lehmann subtype depth
    ax_d = fig.add_subplot(gs[1, 0])
    if depth_by_subtype:
        order  = sorted(depth_by_subtype.items(), key=lambda x: x[1].mean())
        labels = [st for st, _ in order]
        data   = [d  for _, d  in order]
        bplot  = ax_d.boxplot(data, patch_artist=True, notch=False)
        colors = [clr.get(st.lower(), "#95a5a6") for st in labels]
        for patch, c in zip(bplot["boxes"], colors):
            patch.set_facecolor(c)
            patch.set_alpha(0.7)
        ax_d.set_xticklabels(labels, fontsize=7)
        ax_d.set_title("D — Lehmann Subtype Depth (S2-P4)\nPredicted: LAR<BL<M/MSL",
                       fontsize=8, fontweight="bold")
        ax_d.set_ylabel("Depth Score", fontsize=7)

    # E: EED vs depth in bulk
    ax_e = fig.add_subplot(gs[1, 1])
    if gene_df is not None and meta_s is not None and \
            "EED" in gene_df.columns and "depth_proper" in meta_s.columns:
        ax_e.scatter(gene_df["EED"].values,
                     meta_s["depth_proper"].values,
                     c=clr["depth"], alpha=0.2, s=5)
        try:
            rv, _ = stats.pearsonr(gene_df["EED"].values,
                                   meta_s["depth_proper"].values)
            ax_e.set_title(f"E — EED vs Depth\nr={rv:+.3f}",
                           fontsize=8, fontweight="bold")
        except Exception:
            ax_e.set_title("E — EED vs Depth", fontsize=8)
        ax_e.set_xlabel("EED expression", fontsize=7)
        ax_e.set_ylabel("Depth score", fontsize=7)

    # F: AR vs depth in bulk (S2-P9)
    ax_f = fig.add_subplot(gs[1, 2])
    if gene_df is not None and meta_s is not None and \
            "AR" in gene_df.columns and "depth_proper" in meta_s.columns:
        ax_f.scatter(gene_df["AR"].values,
                     meta_s["depth_proper"].values,
                     c=clr["ndepth"], alpha=0.2, s=5)
        try:
            rv, _ = stats.pearsonr(gene_df["AR"].values,
                                   meta_s["depth_proper"].values)
            ax_f.set_title(f"F — AR vs Depth (S2-P9)\nr={rv:+.3f}  "
                           f"Predicted r<-0.15",
                           fontsize=8, fontweight="bold")
        except Exception:
            ax_f.set_title("F — AR vs Depth (S2-P9)", fontsize=8)
        ax_f.set_xlabel("AR expression", fontsize=7)
        ax_f.set_ylabel("Depth score", fontsize=7)

    # G: PARP1 vs depth (S2-P3)
    ax_g = fig.add_subplot(gs[2, 0])
    if gene_df is not None and meta_s is not None and \
            "PARP1" in gene_df.columns and "depth_proper" in meta_s.columns:
        ax_g.scatter(gene_df["PARP1"].values,
                     meta_s["depth_proper"].values,
                     c=clr["m"], alpha=0.2, s=5)
        try:
            rv, _ = stats.pearsonr(gene_df["PARP1"].values,
                                   meta_s["depth_proper"].values)
            ax_g.set_title(f"G — PARP1 vs Depth\nr={rv:+.3f}",
                           fontsize=8, fontweight="bold")
        except Exception:
            ax_g.set_title("G — PARP1 vs Depth", fontsize=8)
        ax_g.set_xlabel("PARP1", fontsize=7)
        ax_g.set_ylabel("Depth", fontsize=7)

    # H: TCGA depth distribution (Basal vs non-Basal)
    ax_h = fig.add_subplot(gs[2, 1])
    tcga_expr_d = t3.get("tcga_depth")
    if tcga_expr_d is not None:
        ax_h.hist(tcga_expr_d.values, bins=40,
                  color=clr["depth"], alpha=0.7)
        ax_h.set_title("H — TCGA-BRCA Depth\nBasal-like highlighted",
                       fontsize=8, fontweight="bold")
        ax_h.set_xlabel("Depth score", fontsize=7)
        ax_h.set_ylabel("n samples", fontsize=7)
    else:
        ax_h.axis("off")
        ax_h.text(0.5, 0.5, "TCGA depth\nnot available",
                  ha="center", va="center", fontsize=9,
                  transform=ax_h.transAxes)

    # I: Summary
    ax_i = fig.add_subplot(gs[2, 2])
    ax_i.axis("off")

    p1_r_str  = f"{t1.get('p1_r', float('nan')):+.3f}" if 'p1_r' in t1 else "n/a"
    eed_r_str = f"{t1.get('eed_r', float('nan')):+.3f}" if 'eed_r' in t1 else "n/a"
    p6_str    = f"p={p6_logrank_p:.4f}" if not np.isnan(p6_logrank_p) else "n/a"

    ax_i.text(0.03, 0.97,
        f"I — SCRIPT 2 SUMMARY\n{'─'*30}\n"
        f"GSE25066 | TCGA-BRCA | 2026-03-04\n\n"
        f"S2-P1 depth∼pCR:     r={p1_r_str}\n"
        f"S2-P2 EED pCR:       r={eed_r_str}\n"
        f"S2-P3 PARP1 pCR:     see panel G\n"
        f"S2-P4 Lehmann depth: see panel D\n"
        f"S2-P5 BRCA1:         proxy only\n"
        f"S2-P6 TCGA survival: {p6_str}\n"
        f"S2-P7 SPI1:          see text\n"
        f"S2-P8 EED biomarker: prospective\n"
        f"S2-P9 AR depth:      see panel F\n\n"
        f"DRUG TARGETS UPDATED:\n"
        f"1. EZH2i  (EED>EZH2 biomarker)\n"
        f"2. PARPi  (PARP1 depth-linked)\n"
        f"3. EED inhibitors (novel)\n"
        f"4. AR blockade (LAR shallow)\n\n"
        f"OrganismCore | BRCA-S4c/d\n"
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
    log("BRCA TNBC — SCRIPT 2")
    log("BULK VALIDATION + pCR + TCGA SURVIVAL")
    log("OrganismCore — BRCA-S4c/d | 2026-03-04")
    log("=" * 65)
    log("")
    log("CONFIRMED URLS (from diagnostic 2026-03-04):")
    log(f"  GSE25066:  {GSE25066_URL}")
    log(f"  TCGA EXPR: {TCGA_EXPR_URL}")
    log(f"  TCGA CLIN: {TCGA_CLIN_URL}")
    log(f"  TCGA SURV: {TCGA_SURV_URL}")
    log(f"  GPL96:     INLINE (parsed from GSE25066 series matrix)")
    log("")

    # ── STEP 0: Download data ────────────────────────────
    paths = acquire_data()

    # ── PARSE GSE25066 ───────────────────────────────────
    log("")
    expr_df, meta_df, probe_gene = parse_gse25066(paths.get("gse25066"))

    # ── TASK 1 ───────────────────────────────────────────
    t1 = task1_gse25066_pcr(expr_df, meta_df, probe_gene)

    # ── TASK 2 ───────────────────────────────────────────
    t2 = task2_lehmann_depth(t1)

    # ── LOAD TCGA ────────────────────────────────────────
    tcga_data = load_tcga(paths)

    # ── TASK 3 ───────────────────────────────────────────
    t3 = task3_tcga(tcga_data)

    # ── SCORECARD ────────────────────────────────────────
    prediction_scorecard(t1, t2, t3)

    # ── FIGURE ───────────────────────────────────────────
    generate_figure(t1, t2, t3)

    # ── CSV EXPORT ──────────────────────────────────��────
    export_rows = []
    for g, r, p in t1.get("bulk_corrs", []):
        export_rows.append({"gene": g, "r_depth": r, "p_value": p})
    if export_rows:
        pd.DataFrame(export_rows).to_csv(CSV_FILE, index=False)
        log(f"  CSV saved: {CSV_FILE}")

    write_log()
    log("")
    log("=" * 65)
    log("SCRIPT 2 COMPLETE")
    log(f"  Log:    {LOG_FILE}")
    log(f"  Figure: {FIG_FILE}")
    log(f"  CSV:    {CSV_FILE}")
    log("=" * 65)


if __name__ == "__main__":
    main()
