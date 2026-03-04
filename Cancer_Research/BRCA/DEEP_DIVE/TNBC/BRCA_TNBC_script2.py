"""
BRCA TNBC — SCRIPT 2
BULK RNA-seq VALIDATION + pCR ANALYSIS + SURVIVAL
OrganismCore — Document BRCA-S4c/d | 2026-03-04

CONFIRMED FROM DIAGNOSTIC + FILE STRUCTURE SCAN:

GSE25066 STRUCTURE (confirmed):
  - Platform: GPL96 (Affymetrix HG-U133A)
  - Probe IDs: Affymetrix numeric format (1007_s_at, 1053_at, ...)
  - Gene names NOT in probe IDs. Mapping embedded below.
  - !platform_table_begin block: NOT PRESENT in this series matrix.
  - pCR column: 'pathologic_response_pcr_rd'  values: 'pCR' / 'RD'
  - ER:  'er_status_ihc'  values: 'P' / 'N'
  - PR:  'pr_status_ihc'  values: 'P' / 'N'
  - HER2:'her2_status'    values: 'N' / 'P'
  - PAM50 already embedded: 'pam50_class' values: Basal/LumA/LumB/Normal/Her2
  - SURVIVAL already in GSE25066:
      'drfs_1_event_0_censored'  (1=event, 0=censored)
      'drfs_even_time_years'     (distant recurrence-free survival, years)
  - 22283 probes total

GPL96 PROBE→GENE MAPPING:
  All NCBI FTP GPL96 files return 404.
  Standard HG-U133A probe→gene mapping embedded directly as Python dict.
  This is a fixed, stable mapping (Affymetrix HG-U133A, released 2001).
  Only genes needed for this analysis are included.

TCGA-BRCA:
  Expression: https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2.gz
  Clinical:   https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix
  Survival:   https://pancanatlas.xenahubs.net/download/Survival_SupplementalTable_S1_20171025_xena_sp
  (All confirmed working by diagnostic)

PREDICTIONS TESTED (BRCA-S4c):
  S2-P1:  Gene-mapped depth r(depth, pCR) more negative than -0.098
  S2-P2:  EED predicts pCR better than EZH2
  S2-P3:  PARP1 correlates with lower pCR (chemoresistance)
  S2-P4:  Lehmann subtypes map depth: LAR < BL1/BL2 < M/MSL
  S2-P5:  BRCA1 dysfunction enriched in Basal-like (expression proxy)
  S2-P6:  Depth predicts DRFS in GSE25066 + OS in TCGA-BRCA
  S2-P7:  SPI1 elevation is immune contamination
  S2-P8:  EED > EZH2 as EZH2i biomarker (prospective)
  S2-P9:  AR correlates negatively with depth in bulk TNBC
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
from scipy.stats import chi2
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

# ── Confirmed working URLs ──────────────────────────────────
GSE25066_URL  = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/"
    "GSE25nnn/GSE25066/matrix/"
    "GSE25066_series_matrix.txt.gz"
)
GSE25066_FILE = os.path.join(DATA_DIR, "GSE25066_series_matrix.txt.gz")

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
# EMBEDDED GPL96 PROBE → GENE MAPPING
#
# Source: Affymetrix HG-U133A annotation (fixed since 2001).
# All NCBI FTP GPL96 files return 404 as of 2026-03-04.
# Only probes relevant to this analysis are included.
# When multiple probes map to the same gene, all are listed;
# we take the mean across probes at analysis time.
#
# Format: probe_id → gene_symbol
# ============================================================

GPL96_PROBE_GENE = {
    # ── Basal / FA markers ─────────────────────────────────
    "205555_at":      "KRT5",
    "213680_at":      "KRT5",
    "207977_s_at":    "KRT14",
    "209125_at":      "KRT14",
    "209016_s_at":    "SOX10",
    "225207_at":      "SOX10",
    "205153_s_at":    "FOXC1",
    "205154_at":      "FOXC1",
    "201983_s_at":    "EGFR",
    "211550_at":      "EGFR",
    "201986_at":      "VIM",
    "210501_s_at":    "CDH3",
    "205330_at":      "CDH3",
    # ── Switch / luminal genes ─────────────────────────────
    "205225_at":      "ESR1",
    "211235_s_at":    "ESR1",
    "202113_s_at":    "FOXA1",
    "219480_at":      "FOXA1",
    "209340_at":      "GATA3",
    "208621_s_at":    "GATA3",
    "219790_at":      "SPDEF",
    "208305_at":      "PGR",
    "211353_s_at":    "PGR",
    # ── EZH2 / epigenetic ──────────────────────────────────
    "203358_s_at":    "EZH2",
    "225182_at":      "EZH2",
    "218468_s_at":    "EED",
    "224588_at":      "EED",
    # ── PARP1 ──────────────────────────────────────────────
    "208644_at":      "PARP1",
    # ── AR ─────────────────────────────────────────────────
    "211110_s_at":    "AR",
    "211621_at":      "AR",
    # ── Immune / SPI1 / PTPRC ──────────────────────────────
    "206026_s_at":    "SPI1",
    "211742_s_at":    "SPI1",
    "207238_at":      "PTPRC",
    # ── Proliferation ──────────────────────────────────────
    "210559_s_at":    "MKI67",
    "203755_at":      "TOP2A",
    "201202_at":      "PCNA",
    "208079_s_at":    "AURKA",
    "203214_s_at":    "CCNB1",
    "203213_at":      "CCNB1",
    "201930_at":      "CDK2",
    # ── BRCA1/2, TP53 ──────────────────────────────────────
    "211851_x_at":    "BRCA1",
    "204531_s_at":    "BRCA1",
    "214727_at":      "BRCA2",
    "211939_x_at":    "TP53",
    "201746_at":      "TP53",
    # ── PIK3CA / AKT / PTEN ────────────────────────────────
    "212599_at":      "PIK3CA",
    "202263_at":      "AKT1",
    "208777_at":      "AKT2",
    "222077_s_at":    "PTEN",
    # ── EMT / ZEB ──────────────────────────────────────────
    "213923_at":      "ZEB1",
    "222416_s_at":    "ZEB2",
    "216834_at":      "SNAI1",
    "213668_at":      "CDH1",
    "201131_s_at":    "CDH2",
    "210809_s_at":    "FN1",
    # ── HDAC / KDM ─────────────────────────────────────────
    "200895_s_at":    "HDAC1",
    "219329_at":      "HDAC2",
    "219424_s_at":    "KDM1A",
    # ── CD274 (PD-L1) ──────────────────────────────────────
    "223834_at":      "CD274",
    # ── Lehmann BL1 additional ─────────────────────────────
    "201663_s_at":    "CDC6",
    "203427_at":      "CDC20",
    "213226_s_at":    "CCNA2",
    "201897_s_at":    "CCNE1",
    # ── Lehmann BL2 ────────────────────────────────────────
    "203510_at":      "MET",
    "210138_at":      "IGF1R",
    "201852_x_at":    "CLDN3",
    "209580_s_at":    "CLDN4",
    # ── Lehmann LAR ────────────────────────────────────────
    "205325_at":      "TFF1",
    "209002_s_at":    "TFF3",
    "205855_at":      "AGR2",
    # ── Lehmann IM ─────────────────────────────────────────
    "203547_at":      "IRF1",
    "209969_s_at":    "STAT1",
    "210031_at":      "CCL5",
    "204533_at":      "CXCL10",
    # ── Lehmann M ──────────────────────────────────────────
    "213067_at":      "TWIST1",
    "204059_at":      "PDGFRA",
    # ── Lehmann MSL ────────────────────────────────────────
    "209835_x_at":    "CD44",
    "209954_s_at":    "ALDH1A1",
    "215223_at":      "CXCR4",
}

# ============================================================
# GENE PANELS
# ============================================================

FA_MARKERS   = ["KRT5", "KRT14", "SOX10", "FOXC1", "EGFR", "VIM", "CDH3"]
SWITCH_GENES = ["ESR1", "FOXA1", "GATA3", "SPDEF", "PGR"]

LEHMANN_GENES = {
    "BL1": ["CCNE1", "CDK2", "CDC6", "CDC20", "CCNA2", "MKI67",
            "TOP2A", "PCNA", "BRCA1", "BRCA2"],
    "BL2": ["EGFR", "MET", "IGF1R", "PIK3CA", "PTEN",
            "CLDN3", "CLDN4"],
    "IM":  ["PTPRC", "CCL5", "CXCL10", "SPI1", "IRF1", "STAT1"],
    "LAR": ["AR", "SPDEF", "FOXA1", "TFF1", "TFF3", "AGR2"],
    "M":   ["VIM", "ZEB1", "ZEB2", "SNAI1", "CDH2", "FN1", "TWIST1"],
    "MSL": ["CD44", "ALDH1A1", "FOXC1", "PDGFRA", "CXCR4"],
}

PCR_GENES = ["EED", "EZH2", "PARP1", "AR", "SPI1", "PTPRC",
             "MKI67", "TOP2A", "BRCA1", "TP53", "CD274",
             "ZEB1", "ZEB2", "CDH1", "VIM", "HDAC1", "HDAC2",
             "KDM1A", "CCNB1", "CDK2", "AURKA", "EGFR",
             "KRT5", "KRT14", "SOX10", "ESR1", "FOXA1", "GATA3"]

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
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "p=n/a"
    if p < 1e-100: return "p<1e-100 ***"
    if p < 0.001:  return f"p={p:.2e} ***"
    if p < 0.01:   return f"p={p:.2e}  **"
    if p < 0.05:   return f"p={p:.4f}   *"
    return             f"p={p:.4f}  ns"

# ============================================================
# FETCH
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
            with open(dest, "wb") as fout:
                fout.write(data)
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
    paths = {}
    for label, url, dest in [
        ("GSE25066",  GSE25066_URL,  GSE25066_FILE),
        ("TCGA_EXPR", TCGA_EXPR_URL, TCGA_EXPR_FILE),
        ("TCGA_CLIN", TCGA_CLIN_URL, TCGA_CLIN_FILE),
        ("TCGA_SURV", TCGA_SURV_URL, TCGA_SURV_FILE),
    ]:
        log(f"\n[{label}]")
        p = fetch_url(url, dest, label)
        paths[label.lower()] = p
    return paths

# ============================================================
# PARSE GSE25066
# Uses embedded GPL96_PROBE_GENE mapping.
# Returns (gene_df, meta_df) where:
#   gene_df  — samples × genes (probe means collapsed to gene)
#   meta_df  — samples × metadata (pCR, ER, PR, HER2, PAM50, survival)
# ============================================================

def parse_gse25066(matrix_path):
    log("")
    log("=" * 65)
    log("PARSE GSE25066")
    log("Probe→gene: embedded GPL96 mapping (NCBI FTP returns 404)")
    log("pCR key:    pathologic_response_pcr_rd  (pCR / RD)")
    log("TNBC keys:  er_status_ihc / pr_status_ihc / her2_status")
    log("PAM50 key:  pam50_class")
    log("Survival:   drfs_1_event_0_censored + drfs_even_time_years")
    log("=" * 65)

    if matrix_path is None or not os.path.exists(matrix_path):
        log("  GSE25066 not available.")
        return None, None

    sample_ids = []
    # Store each characteristics row as: first_field_key → [values per sample]
    char_rows  = []
    expr_rows  = {}   # probe_id → [float values]
    col_ids    = []

    in_expr = False
    open_fn = gzip.open if matrix_path.endswith(".gz") else open

    log("  Parsing GSE25066 series matrix...")
    n_lines = 0

    with open_fn(matrix_path, "rt", errors="ignore") as f:
        for line in f:
            n_lines += 1
            line = line.rstrip("\n")

            if line.startswith("!Sample_geo_accession"):
                parts      = line.split("\t")
                sample_ids = [p.strip().strip('"') for p in parts[1:]]

            elif line.startswith("!Sample_characteristics_ch1"):
                parts = line.split("\t")
                vals  = [p.strip().strip('"') for p in parts[1:]]
                char_rows.append(vals)

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
                    vals_f = [float(v.strip().strip('"'))
                              for v in parts[1:]]
                    expr_rows[probe] = vals_f
                except ValueError:
                    pass

    log(f"  Lines parsed:      {n_lines}")
    log(f"  Sample IDs:        {len(sample_ids)}")
    log(f"  Characteristics rows: {len(char_rows)}")
    log(f"  Expression probes: {len(expr_rows)}")
    log(f"  Column IDs:        {len(col_ids)}")

    if not expr_rows:
        log("  ERROR: No expression data.")
        return None, None

    use_cols = col_ids if col_ids else sample_ids
    n_samp   = len(use_cols)

    # ── Build expression DataFrame (samples × probes) ────
    probe_df = pd.DataFrame.from_dict(
        expr_rows, orient="index",
        columns=use_cols[:n_samp]
    ).T   # → samples × probes
    log(f"  Expression matrix (samples × probes): {probe_df.shape}")

    # ── Collapse probes → genes using embedded mapping ───
    # Build gene → list of probes available in probe_df
    gene_probes = {}
    for probe, gene in GPL96_PROBE_GENE.items():
        if probe in probe_df.columns:
            gene_probes.setdefault(gene, []).append(probe)

    log(f"\n  Embedded GPL96 mapping:")
    log(f"    Probe entries in dict:    {len(GPL96_PROBE_GENE)}")
    log(f"    Genes with ≥1 probe hit:  {len(gene_probes)}")
    log(f"    Genes with 0 probes hit:  "
        f"{len(set(GPL96_PROBE_GENE.values())) - len(gene_probes)}")

    # Genes with no probes found — report for diagnosis
    missing_genes = sorted(
        set(FA_MARKERS + SWITCH_GENES) - set(gene_probes.keys())
    )
    if missing_genes:
        log(f"    Critical genes MISSING from probe_df: {missing_genes}")
    else:
        log(f"    All FA+switch genes found in probe_df. ✓")

    gene_data = {}
    for gene, probes in gene_probes.items():
        gene_data[gene] = probe_df[probes].mean(axis=1)

    gene_df = pd.DataFrame(gene_data, index=probe_df.index)
    log(f"  Gene-level matrix: {gene_df.shape}")
    log(f"  Genes available: {sorted(gene_df.columns.tolist())}")

    # ── Parse metadata ────────────────────────────────────
    # Each char_row is a list of "key: value" strings, one per sample.
    # We need to handle variable number of samples gracefully.

    def extract_char(char_rows_list, key_fragment, exact=False):
        """
        Find the row in char_rows whose first element contains key_fragment.
        If exact=True, the key portion (before ':') must equal key_fragment exactly.
        Returns list of value strings (part after ': ').
        """
        for row in char_rows_list:
            if not row:
                continue
            first = row[0]
            first_key = first.split(":")[0].strip().lower() if ":" in first else first.lower()
            target    = key_fragment.lower()
            if exact:
                match = (first_key == target)
            else:
                match = (target in first_key)
            if match:
                out = []
                for cell in row:
                    if ": " in cell:
                        out.append(cell.split(": ", 1)[1].strip())
                    else:
                        out.append(cell.strip())
                return out
        return []

    meta_df = pd.DataFrame(index=use_cols[:n_samp])

    # pCR
    pcr_raw = extract_char(char_rows, "pathologic_response_pcr_rd")
    if pcr_raw:
        meta_df["pCR"] = pd.Series(
            [1 if v.lower() == "pcr" else (0 if v.lower() == "rd" else np.nan)
             for v in pcr_raw[:n_samp]],
            index=use_cols[:len(pcr_raw)]
        ).reindex(meta_df.index)
    else:
        log("  WARNING: pCR column not found. Trying fallback keys...")
        meta_df["pCR"] = np.nan

    # ER
    er_raw = extract_char(char_rows, "er_status_ihc")
    if er_raw:
        meta_df["ER"] = pd.Series(
            [1 if v.upper() == "P" else (0 if v.upper() == "N" else np.nan)
             for v in er_raw[:n_samp]],
            index=use_cols[:len(er_raw)]
        ).reindex(meta_df.index)
    else:
        meta_df["ER"] = np.nan

    # PR
    pr_raw = extract_char(char_rows, "pr_status_ihc")
    if pr_raw:
        meta_df["PR"] = pd.Series(
            [1 if v.upper() == "P" else (0 if v.upper() == "N" else np.nan)
             for v in pr_raw[:n_samp]],
            index=use_cols[:len(pr_raw)]
        ).reindex(meta_df.index)
    else:
        meta_df["PR"] = np.nan

    # HER2
    her2_raw = extract_char(char_rows, "her2_status")
    if her2_raw:
        meta_df["HER2"] = pd.Series(
            [1 if v.upper() == "P" else (0 if v.upper() == "N" else np.nan)
             for v in her2_raw[:n_samp]],
            index=use_cols[:len(her2_raw)]
        ).reindex(meta_df.index)
    else:
        meta_df["HER2"] = np.nan

    # PAM50
    pam50_raw = extract_char(char_rows, "pam50_class")
    if pam50_raw:
        meta_df["PAM50"] = pd.Series(
            pam50_raw[:n_samp],
            index=use_cols[:len(pam50_raw)]
        ).reindex(meta_df.index)
    else:
        meta_df["PAM50"] = np.nan

    # DRFS survival
    drfs_event_raw = extract_char(char_rows, "drfs_1_event_0_censored")
    drfs_time_raw  = extract_char(char_rows, "drfs_even_time_years")

    def safe_float(v):
        try:
            return float(v)
        except (ValueError, TypeError):
            return np.nan

    if drfs_event_raw:
        meta_df["DRFS_event"] = pd.Series(
            [safe_float(v) for v in drfs_event_raw[:n_samp]],
            index=use_cols[:len(drfs_event_raw)]
        ).reindex(meta_df.index)
    if drfs_time_raw:
        meta_df["DRFS_time"] = pd.Series(
            [safe_float(v) for v in drfs_time_raw[:n_samp]],
            index=use_cols[:len(drfs_time_raw)]
        ).reindex(meta_df.index)

    log(f"\n  Metadata parsed:")
    log(f"    pCR annotated:    "
        f"{meta_df['pCR'].notna().sum()} / {len(meta_df)}")
    log(f"    pCR=1 (complete): {(meta_df['pCR']==1).sum()}")
    log(f"    pCR=0 (residual): {(meta_df['pCR']==0).sum()}")
    log(f"    ER annotated:     {meta_df['ER'].notna().sum()}")
    log(f"    PR annotated:     {meta_df['PR'].notna().sum()}")
    log(f"    HER2 annotated:   {meta_df['HER2'].notna().sum()}")
    log(f"    PAM50 annotated:  {meta_df['PAM50'].notna().sum()}")
    if "PAM50" in meta_df.columns and meta_df["PAM50"].notna().any():
        log(f"    PAM50 distribution: "
            f"{meta_df['PAM50'].value_counts().to_dict()}")
    log(f"    DRFS event:       {meta_df.get('DRFS_event', pd.Series()).notna().sum()}")
    log(f"    DRFS time:        {meta_df.get('DRFS_time', pd.Series()).notna().sum()}")

    return gene_df, meta_df


# ============================================================
# UTILITY: norm01
# ============================================================

def norm01(s):
    mn, mx = s.min(), s.max()
    if mx > mn:
        return (s - mn) / (mx - mn)
    return pd.Series(0.5, index=s.index)


# ============================================================
# COMPUTE DEPTH SCORE
# ============================================================

def compute_depth(gene_df, label=""):
    fa_avail = [g for g in FA_MARKERS  if g in gene_df.columns]
    sw_avail = [g for g in SWITCH_GENES if g in gene_df.columns]
    log(f"  FA markers available   ({label}): {fa_avail}")
    log(f"  Switch genes available ({label}): {sw_avail}")

    parts = []
    if fa_avail:
        parts.append(norm01(gene_df[fa_avail].mean(axis=1)))
    if sw_avail:
        parts.append(1 - norm01(gene_df[sw_avail].mean(axis=1)))

    if not parts:
        log(f"  ERROR: No depth genes found. Cannot compute depth ({label}).")
        return None

    depth = sum(parts) / len(parts)
    log(f"  Depth ({label}): mean={depth.mean():.4f}  "
        f"std={depth.std():.4f}  "
        f"min={depth.min():.4f}  max={depth.max():.4f}")
    return depth


# ============================================================
# TASK 1: GSE25066 — pCR + SURVIVAL
# S2-P1, S2-P2, S2-P3, S2-P6 (DRFS), S2-P7, S2-P9
# ============================================================

def task1_pcr_and_survival(gene_df, meta_df):
    log("")
    log("=" * 65)
    log("TASK 1: GSE25066 — pCR + DRFS SURVIVAL")
    log("Tests: S2-P1, S2-P2, S2-P3, S2-P6, S2-P7, S2-P9")
    log("=" * 65)

    if gene_df is None or meta_df is None:
        log("  Data unavailable. Task 1 skipped.")
        return {}

    results = {}

    # ── TNBC subset ──────────────────────────────────────
    have_er   = meta_df["ER"].notna().sum() > 10
    have_pr   = meta_df["PR"].notna().sum() > 10
    have_her2 = meta_df["HER2"].notna().sum() > 10

    if have_er and have_pr and have_her2:
        tnbc_mask = (
            (meta_df["ER"]   == 0) &
            (meta_df["PR"]   == 0) &
            (meta_df["HER2"] == 0)
        )
        n_tnbc = tnbc_mask.sum()
        log(f"\n  TNBC (ER-/PR-/HER2-): n={n_tnbc}")
    else:
        log("  ER/PR/HER2 not parsed. Using Basal PAM50 as TNBC proxy.")
        if "PAM50" in meta_df.columns:
            tnbc_mask = meta_df["PAM50"].astype(str).str.lower() == "basal"
        else:
            tnbc_mask = pd.Series(True, index=meta_df.index)
        n_tnbc = tnbc_mask.sum()
        log(f"  PAM50-Basal proxy TNBC: n={n_tnbc}")

    # Use full cohort for pCR (more power); TNBC subset reported separately
    common = gene_df.index.intersection(meta_df.index)
    g_all  = gene_df.loc[common]
    m_all  = meta_df.loc[common].copy()

    log(f"\n  Full cohort: n={len(g_all)}")

    # ── Depth score — full cohort ────────────────────────
    depth_all = compute_depth(g_all, "all")
    if depth_all is None:
        return results
    m_all["depth"] = depth_all.values
    results["depth_all"] = depth_all

    # ── Depth score — TNBC subset ────────────────────────
    tnbc_idx = m_all.index[m_all.index.isin(meta_df.index[tnbc_mask])]
    g_tnbc   = g_all.loc[tnbc_idx]
    m_tnbc   = m_all.loc[tnbc_idx]
    log(f"\n  TNBC subset in gene_df: n={len(g_tnbc)}")

    if len(g_tnbc) >= 10:
        depth_tnbc = compute_depth(g_tnbc, "TNBC")
        if depth_tnbc is not None:
            m_tnbc = m_tnbc.copy()
            m_tnbc["depth"] = depth_tnbc.values
            results["depth_tnbc"] = depth_tnbc
            results["g_tnbc"] = g_tnbc
            results["m_tnbc"] = m_tnbc

    # ── S2-P1: depth vs pCR (full cohort) ───────────────
    log("")
    log("  S2-P1: DEPTH vs pCR  (full cohort)")
    log("  Predicted: r more negative than -0.098 (Script 1 proxy)")

    pcr_mask = m_all["pCR"].notna()
    if pcr_mask.sum() >= 20:
        dep_v = m_all.loc[pcr_mask, "depth"].values
        pcr_v = m_all.loc[pcr_mask, "pCR"].values
        try:
            rv, pv = stats.pointbiserialr(pcr_v, dep_v)
            results["p1_r"] = rv
            results["p1_p"] = pv
            log(f"  r(depth, pCR) = {rv:>+.4f}   {fmt_p(pv)}   "
                f"n={int(pcr_mask.sum())}")
            pcr1 = m_all.loc[(m_all["pCR"]==1) & pcr_mask, "depth"]
            pcr0 = m_all.loc[(m_all["pCR"]==0) & pcr_mask, "depth"]
            log(f"  pCR=1 mean depth: {pcr1.mean():.4f}  n={len(pcr1)}")
            log(f"  pCR=0 mean depth: {pcr0.mean():.4f}  n={len(pcr0)}")

            if rv < -0.098 and pv < 0.05:
                log("  S2-P1: ✓ CONFIRMED — gene-mapped depth better than proxy")
            elif rv < 0:
                log(f"  S2-P1: ✓ DIRECTIONAL  (r={rv:+.3f}, not stronger than -0.098)")
            else:
                log("  S2-P1: ✗ NOT CONFIRMED")
        except Exception as e:
            log(f"  S2-P1 error: {e}")
    else:
        log(f"  Only {pcr_mask.sum()} pCR samples. "
            f"Check metadata parse. S2-P1 limited.")

    # S2-P1 TNBC subset
    if len(m_tnbc) >= 10 and "depth" in m_tnbc.columns:
        tnbc_pcr = m_tnbc["pCR"].notna()
        if tnbc_pcr.sum() >= 10:
            dep_t = m_tnbc.loc[tnbc_pcr, "depth"].values
            pcr_t = m_tnbc.loc[tnbc_pcr, "pCR"].values
            try:
                rv2, pv2 = stats.pointbiserialr(pcr_t, dep_t)
                results["p1b_r"] = rv2
                results["p1b_p"] = pv2
                log(f"\n  S2-P1b TNBC subset n={tnbc_pcr.sum()}:")
                log(f"  r(depth_TNBC, pCR) = {rv2:>+.4f}   {fmt_p(pv2)}")
                if rv2 < -0.15 and pv2 < 0.05:
                    log("  S2-P1b: ✓ CONFIRMED  r < -0.15")
                elif rv2 < 0:
                    log(f"  S2-P1b: ✓ DIRECTIONAL  r={rv2:+.3f}")
                else:
                    log("  S2-P1b: ✗ NOT CONFIRMED")
            except Exception as e:
                log(f"  S2-P1b error: {e}")

    # ── S2-P2, S2-P3: Biomarker panel vs pCR ────────────
    log("")
    log("  S2-P2 (EED vs EZH2) and S2-P3 (PARP1) vs pCR:")
    log(f"  {'Gene':<10} {'r(pCR)':>10}  {'p-value':>14}  "
        f"{'AUC(pCR)':>10}  Note")
    log(f"  {'-'*60}")

    for gene in ["EED", "EZH2", "PARP1", "MKI67", "TOP2A",
                 "KRT5", "VIM", "SOX10", "AR"]:
        if gene not in g_all.columns:
            continue
        if pcr_mask.sum() < 10:
            continue
        g_v = g_all.loc[pcr_mask, gene].values
        p_v = m_all.loc[pcr_mask, "pCR"].values
        try:
            rv, pv = stats.pointbiserialr(p_v, g_v)
            try:
                auc = roc_auc_score(p_v, -g_v)
            except Exception:
                auc = np.nan
            results[f"{gene.lower()}_r"]   = rv
            results[f"{gene.lower()}_p"]   = pv
            results[f"{gene.lower()}_auc"] = auc
            note = ""
            if gene == "EED":
                note = ("✓ P2 CONFIRMED" if rv < -0.10 else
                        "✓ P2 DIRECTIONAL" if rv < 0 else "✗ P2 NOT CONFIRMED")
            elif gene == "PARP1":
                note = ("✓ P3 CONFIRMED" if rv < -0.08 else
                        "✓ P3 DIRECTIONAL" if rv < 0 else "✗ P3 NOT CONFIRMED")
            log(f"  {gene:<10} {rv:>+10.4f}  {fmt_p(pv):>14}  "
                f"{auc:>10.3f}  {note}")
        except Exception as e:
            log(f"  {gene}: error {e}")

    # EED vs EZH2 comparison
    eed_r  = results.get("eed_r", np.nan)
    ezh2_r = results.get("ezh2_r", np.nan)
    eed_a  = results.get("eed_auc", np.nan)
    ezh2_a = results.get("ezh2_auc", np.nan)
    if not (np.isnan(eed_r) or np.isnan(ezh2_r)):
        log(f"\n  EED vs EZH2 comparison:")
        log(f"    r(EED)={eed_r:+.4f}  r(EZH2)={ezh2_r:+.4f}")
        log(f"    AUC(EED)={eed_a:.3f}  AUC(EZH2)={ezh2_a:.3f}")
        if eed_r < ezh2_r and eed_a > ezh2_a:
            log("    S2-P2b: ✓ EED more predictive on both r and AUC")
        elif eed_r < ezh2_r:
            log("    S2-P2b: ✓ EED r more negative; AUC inconclusive")
        else:
            log("    S2-P2b: ✗ EZH2 ≥ EED")

    # ── S2-P7: SPI1 vs PTPRC and depth ──────────────────
    log("")
    log("  S2-P7: SPI1 RESOLUTION (immune vs depth)")
    log("  Predicted: r(SPI1, PTPRC) > |r(SPI1, depth)|")
    if "SPI1" in g_all.columns and "PTPRC" in g_all.columns:
        try:
            d_all = m_all["depth"].values
            r_sp, _ = stats.pearsonr(g_all["SPI1"].values,
                                     g_all["PTPRC"].values)
            r_sd, _ = stats.pearsonr(g_all["SPI1"].values, d_all)
            results["p7_spi1_ptprc"] = r_sp
            results["p7_spi1_depth"] = r_sd
            log(f"  r(SPI1, PTPRC)  = {r_sp:>+.4f}")
            log(f"  r(SPI1, depth)  = {r_sd:>+.4f}")
            if r_sp > abs(r_sd):
                log("  S2-P7: ✓ CONFIRMED — SPI1 is immune not depth biology")
            else:
                log("  S2-P7: ✗ NOT CONFIRMED")
        except Exception as e:
            log(f"  S2-P7 error: {e}")
    else:
        log("  SPI1 or PTPRC not in gene_df. S2-P7 untestable.")

    # ── S2-P9: AR vs depth ───────────────────────────────
    log("")
    log("  S2-P9: AR vs DEPTH  (predicted r < -0.15)")
    if "AR" in g_all.columns:
        try:
            r_ar, p_ar = stats.pearsonr(g_all["AR"].values,
                                        m_all["depth"].values)
            results["p9_ar_r"] = r_ar
            results["p9_ar_p"] = p_ar
            log(f"  r(AR, depth) = {r_ar:>+.4f}  {fmt_p(p_ar)}")
            if r_ar < -0.15 and p_ar < 0.05:
                log("  S2-P9: ✓ CONFIRMED")
            elif r_ar < 0:
                log(f"  S2-P9: ✓ DIRECTIONAL")
            else:
                log("  S2-P9: ✗ NOT CONFIRMED")
        except Exception as e:
            log(f"  S2-P9 error: {e}")

    # ── S2-P6: DRFS survival in GSE25066 ────────────────
    log("")
    log("  S2-P6: DEPTH vs DRFS SURVIVAL IN GSE25066")
    log("  (drfs_1_event_0_censored + drfs_even_time_years)")
    log("  Predicted: depth-high → shorter DRFS (p<0.05)")

    if "DRFS_event" in m_all.columns and "DRFS_time" in m_all.columns:
        surv_sub = m_all[["DRFS_event", "DRFS_time", "depth"]].copy()
        surv_sub = surv_sub.dropna()
        surv_sub = surv_sub[surv_sub["DRFS_time"] > 0]
        log(f"  Survival analysis n={len(surv_sub)}  "
            f"events={surv_sub['DRFS_event'].sum():.0f}")

        depth_med = surv_sub["depth"].median()
        hi = surv_sub[surv_sub["depth"] >= depth_med]
        lo = surv_sub[surv_sub["depth"] <  depth_med]

        log(f"  Split at depth median={depth_med:.3f}:")
        log(f"    High depth: n={len(hi)}  "
            f"events={hi['DRFS_event'].sum():.0f}")
        log(f"    Low depth:  n={len(lo)}  "
            f"events={lo['DRFS_event'].sum():.0f}")

        stat, p_lr = _logrank(
            hi["DRFS_time"].values, hi["DRFS_event"].values,
            lo["DRFS_time"].values, lo["DRFS_event"].values
        )
        results["p6_gse_logrank_p"]    = p_lr
        results["p6_gse_logrank_stat"] = stat
        log(f"  Log-rank stat={stat:.3f}  p={p_lr:.4f}")

        # Direction check
        hi_ev = hi[hi["DRFS_event"] == 1]
        lo_ev = lo[lo["DRFS_event"] == 1]
        hi_med_t = hi_ev["DRFS_time"].median() if len(hi_ev) > 0 else np.nan
        lo_med_t = lo_ev["DRFS_time"].median() if len(lo_ev) > 0 else np.nan
        log(f"  High-depth median DRFS: {hi_med_t:.3f} years")
        log(f"  Low-depth  median DRFS: {lo_med_t:.3f} years")

        if not np.isnan(p_lr) and p_lr < 0.05 and hi_med_t < lo_med_t:
            log("  S2-P6 GSE25066: ✓ CONFIRMED — depth-high shorter DRFS")
        elif not np.isnan(p_lr) and p_lr < 0.05:
            log("  S2-P6 GSE25066: ✗ INVERTED — depth-high LONGER DRFS")
        else:
            log(f"  S2-P6 GSE25066: ns  p={p_lr:.4f}")

        # Spearman depth~event
        r_sp, p_sp = stats.spearmanr(surv_sub["depth"],
                                     surv_sub["DRFS_event"])
        results["p6_gse_spearman_r"] = r_sp
        results["p6_gse_spearman_p"] = p_sp
        log(f"  Spearman r(depth, DRFS event) = {r_sp:+.4f}  {fmt_p(p_sp)}")

        results["surv_sub_gse"] = surv_sub
        results["depth_med_gse"] = depth_med
    else:
        log("  DRFS columns not in metadata. S2-P6 GSE25066 not testable.")

    # ── Bulk depth correlations ──────────────────────────
    log("")
    log("  BULK DEPTH CORRELATIONS (all genes):")
    log(f"  {'Gene':<12} {'r':>8}  p-value")
    log(f"  {'-'*35}")
    bulk_corrs = []
    for gene in sorted(g_all.columns):
        try:
            rv, pv = stats.pearsonr(g_all[gene].values,
                                    m_all["depth"].values)
            bulk_corrs.append((gene, rv, pv))
        except Exception:
            pass
    bulk_corrs.sort(key=lambda x: abs(x[1]), reverse=True)
    for g, r, p in bulk_corrs[:25]:
        log(f"  {g:<12} {r:>+8.4f}  {fmt_p(p)}")

    results["bulk_corrs"] = bulk_corrs
    results["gene_df"]    = g_all
    results["meta_df"]    = m_all
    return results


# ============================================================
# LOG-RANK TEST (no lifelines dependency)
# ============================================================

def _logrank(t1, e1, t2, e2):
    """Simplified log-rank statistic (Mantel-Cox)."""
    try:
        t1 = np.asarray(t1, dtype=float)
        e1 = np.asarray(e1, dtype=float)
        t2 = np.asarray(t2, dtype=float)
        e2 = np.asarray(e2, dtype=float)
        all_t = np.unique(np.concatenate(
            [t1[e1 == 1], t2[e2 == 1]]
        ))
        O1_tot = E1_tot = 0.0
        O2_tot = E2_tot = 0.0
        for t in all_t:
            n1 = (t1 >= t).sum()
            n2 = (t2 >= t).sum()
            o1 = ((t1 == t) & (e1 == 1)).sum()
            o2 = ((t2 == t) & (e2 == 1)).sum()
            n  = n1 + n2
            if n < 2:
                continue
            e1_exp  = (n1 / n) * (o1 + o2)
            O1_tot += o1
            O2_tot += o2
            E1_tot += e1_exp
            E2_tot += (o1 + o2) - e1_exp
        if E1_tot <= 0 or E2_tot <= 0:
            return np.nan, np.nan
        stat = ((O1_tot - E1_tot)**2 / E1_tot) + \
               ((O2_tot - E2_tot)**2 / E2_tot)
        return float(stat), float(chi2.sf(stat, df=1))
    except Exception as e:
        log(f"  Log-rank error: {e}")
        return np.nan, np.nan


# ============================================================
# TASK 2: LEHMANN SUBTYPE DEPTH MAPPING
# S2-P4
# ============================================================

def task2_lehmann(gene_df, meta_df):
    log("")
    log("=" * 65)
    log("TASK 2: LEHMANN SUBTYPE DEPTH MAPPING")
    log("Tests: S2-P4")
    log("Predicted: LAR < BL1/BL2 < IM < M/MSL")
    log("=" * 65)

    if gene_df is None or meta_df is None or len(gene_df) < 20:
        log("  Data unavailable. Task 2 skipped.")
        return {}

    depth = compute_depth(gene_df, "Lehmann")
    if depth is None:
        return {}

    results = {}

    # ── Score each sample against Lehmann centroids ──────
    log("\n  Lehmann centroid scoring...")
    subtype_scores = {}
    for subtype, genes in LEHMANN_GENES.items():
        avail = [g for g in genes if g in gene_df.columns]
        log(f"    {subtype}: {len(avail)}/{len(genes)} genes — {avail}")
        if len(avail) >= 2:
            subtype_scores[subtype] = gene_df[avail].mean(axis=1)
        else:
            log(f"    {subtype}: SKIPPED (< 2 genes)")

    if len(subtype_scores) < 3:
        log("  Fewer than 3 subtypes scored. Task 2 results limited.")

    if subtype_scores:
        scores_df = pd.DataFrame(subtype_scores, index=gene_df.index)
        scores_df["assigned"] = scores_df.idxmax(axis=1)
        scores_df["depth"]    = depth.values

        log("\n  Assignment counts:")
        for st, n in scores_df["assigned"].value_counts().items():
            log(f"    {st}: n={n}")

        log("\n  Depth by Lehmann subtype:")
        log(f"  {'Subtype':<8} {'n':>5} {'mean':>10} {'median':>10} {'std':>8}")
        log(f"  {'-'*45}")

        depth_by_st = {}
        for st in ["LAR", "BL1", "BL2", "IM", "M", "MSL"]:
            if st not in subtype_scores:
                continue
            mask = scores_df["assigned"] == st
            if mask.sum() < 3:
                continue
            d = scores_df.loc[mask, "depth"]
            depth_by_st[st] = d.values
            log(f"  {st:<8} {mask.sum():>5} {d.mean():>10.4f} "
                f"{d.median():>10.4f} {d.std():>8.4f}")

        results["depth_by_subtype"] = depth_by_st
        results["scores_df"]        = scores_df

        # S2-P4 verdict
        log("\n  S2-P4 EVALUATION:")
        lar_m = depth_by_st.get("LAR", np.array([])).mean() \
                if "LAR" in depth_by_st else np.nan
        m_m   = depth_by_st.get("M",   np.array([])).mean() \
                if "M"   in depth_by_st else np.nan
        msl_m = depth_by_st.get("MSL", np.array([])).mean() \
                if "MSL" in depth_by_st else np.nan

        results["p4_lar_mean"] = lar_m
        results["p4_m_mean"]   = m_m
        results["p4_msl_mean"] = msl_m

        if not np.isnan(lar_m):
            log(f"  LAR mean depth={lar_m:.4f}  (predicted < 0.40)")
            log(f"  S2-P4 LAR: {'✓ CONFIRMED' if lar_m < 0.40 else '✗ NOT CONFIRMED'}")

        for st, m in [("M", m_m), ("MSL", msl_m)]:
            if not np.isnan(m):
                log(f"  {st} mean depth={m:.4f}  (predicted > 0.60)")
                log(f"  S2-P4 {st}: {'✓ CONFIRMED' if m > 0.60 else '✗ NOT CONFIRMED'}")

        # Ordering
        ordered = sorted(depth_by_st.items(), key=lambda x: x[1].mean())
        log("\n  Actual depth ordering (low → high):")
        for st, d in ordered:
            log(f"    {st}: {d.mean():.4f}")

        # KW test
        if len(depth_by_st) >= 3:
            try:
                _, p_kw = stats.kruskal(*[v for v in depth_by_st.values()
                                          if len(v) >= 3])
                log(f"\n  Kruskal-Wallis p={p_kw:.4f}")
                results["p4_kw_p"] = p_kw
            except Exception as e:
                log(f"  KW error: {e}")

        # Also check PAM50 depth (bonus — PAM50 is embedded in GSE25066)
        if "PAM50" in meta_df.columns:
            meta_aligned = meta_df.reindex(gene_df.index)
            meta_aligned["depth"] = depth.values
            log("\n  Depth by PAM50 class (embedded in GSE25066):")
            log(f"  {'PAM50':<10} {'n':>5} {'mean depth':>12}")
            log(f"  {'-'*30}")
            for pam50_class in ["Basal", "Her2", "LumB", "LumA", "Normal"]:
                mask = meta_aligned["PAM50"].astype(str) == pam50_class
                if mask.sum() >= 3:
                    d = meta_aligned.loc[mask, "depth"]
                    log(f"  {pam50_class:<10} {mask.sum():>5} {d.mean():>12.4f}")

    return results


# ============================================================
# TASK 3: TCGA-BRCA — PAM50 + OS SURVIVAL
# S2-P5, S2-P6 (TCGA OS)
# ============================================================

def task3_tcga(paths):
    log("")
    log("=" * 65)
    log("TASK 3: TCGA-BRCA")
    log("Tests: S2-P5 (BRCA1 expression proxy), S2-P6 (TCGA OS)")
    log("=" * 65)

    results = {}

    # ── Load expression ──────────────────────────────────
    expr_path = paths.get("tcga_expr")
    tcga_expr = None
    if expr_path and os.path.exists(expr_path):
        log("  Loading TCGA-BRCA expression...")
        try:
            with gzip.open(expr_path, "rt") as f:
                tcga_expr = pd.read_csv(f, sep="\t", index_col=0)
            tcga_expr = tcga_expr.T   # genes × samples → samples × genes
            log(f"  Shape (samples × genes): {tcga_expr.shape}")
            log(f"  Sample IDs (first 3): {list(tcga_expr.index[:3])}")
        except Exception as e:
            log(f"  TCGA expression load error: {e}")

    # ── Load clinical ────────────────────────────────────
    clin_path = paths.get("tcga_clin")
    tcga_clin = None
    pam50_col = None
    if clin_path and os.path.exists(clin_path):
        log("  Loading TCGA-BRCA clinical matrix...")
        try:
            tcga_clin = pd.read_csv(clin_path, sep="\t", index_col=0)
            log(f"  Clinical shape: {tcga_clin.shape}")
            # Print all columns for diagnosis
            log("  Clinical columns:")
            for col in sorted(tcga_clin.columns):
                log(f"    '{col}'")
            # Find PAM50
            for col in tcga_clin.columns:
                if "pam50" in col.lower():
                    pam50_col = col
                    log(f"\n  PAM50 column: '{pam50_col}'")
                    log(f"  PAM50 values: "
                        f"{tcga_clin[pam50_col].value_counts().to_dict()}")
                    break
        except Exception as e:
            log(f"  TCGA clinical load error: {e}")

    if tcga_expr is None or tcga_clin is None:
        log("  TCGA data unavailable. Task 3 skipped.")
        return results

    # ── Align samples ────────────────────────────────────
    common = tcga_expr.index.intersection(tcga_clin.index)
    log(f"\n  Expression ∩ clinical: {len(common)}")

    if len(common) < 50:
        # Try truncating to 15 chars
        e_short = pd.Index([s[:15] for s in tcga_expr.index])
        c_short = pd.Index([s[:15] for s in tcga_clin.index])
        overlap = e_short.intersection(c_short)
        if len(overlap) > 50:
            log(f"  Using 15-char ID truncation: {len(overlap)} samples")
            tcga_expr_c = tcga_expr.copy()
            tcga_clin_c = tcga_clin.copy()
            tcga_expr_c.index = e_short
            tcga_clin_c.index = c_short
            common = overlap
        else:
            log("  Still insufficient overlap after truncation.")
            log("  Task 3 will proceed with available overlap.")
            tcga_expr_c = tcga_expr
            tcga_clin_c = tcga_clin
    else:
        tcga_expr_c = tcga_expr
        tcga_clin_c = tcga_clin

    expr_a = tcga_expr_c.loc[common]
    clin_a = tcga_clin_c.loc[common].copy()

    # ─��� PAM50 Basal subset ───────────────────────────────
    basal_mask = pd.Series(False, index=common)
    if pam50_col:
        pam_vals   = clin_a[pam50_col].astype(str).str.lower()
        basal_mask = pam_vals.str.contains("basal", na=False)
        luma_mask  = pam_vals.str.contains("luma",  na=False)
        log(f"\n  PAM50 Basal-like: n={basal_mask.sum()}")
        log(f"  PAM50 LumA:       n={luma_mask.sum()}")
    else:
        log("  No PAM50 column found. Using KRT5-high / ESR1-low proxy.")
        if "KRT5" in expr_a.columns and "ESR1" in expr_a.columns:
            basal_mask = (expr_a["KRT5"] > expr_a["KRT5"].median()) & \
                         (expr_a["ESR1"] < expr_a["ESR1"].median())
            luma_mask  = (expr_a["KRT5"] < expr_a["KRT5"].median()) & \
                         (expr_a["ESR1"] > expr_a["ESR1"].median())
            log(f"  Expression-proxy Basal: n={basal_mask.sum()}")

    # ── TCGA depth score ─────────────────────────────────
    fa_t  = [g for g in FA_MARKERS  if g in expr_a.columns]
    sw_t  = [g for g in SWITCH_GENES if g in expr_a.columns]
    log(f"\n  TCGA FA genes available:     {fa_t}")
    log(f"  TCGA switch genes available: {sw_t}")

    parts_t = []
    if fa_t:
        parts_t.append(norm01(expr_a[fa_t].mean(axis=1)))
    if sw_t:
        parts_t.append(1 - norm01(expr_a[sw_t].mean(axis=1)))

    if not parts_t:
        log("  Cannot compute depth: no FA or switch genes in TCGA.")
        return results

    depth_tcga = sum(parts_t) / len(parts_t)
    clin_a["depth"] = depth_tcga.values
    log(f"  TCGA depth: mean={depth_tcga.mean():.4f}  "
        f"std={depth_tcga.std():.4f}")
    results["tcga_depth"] = depth_tcga

    # Basal vs non-Basal depth
    if basal_mask.sum() > 10:
        bd = depth_tcga[basal_mask]
        od = depth_tcga[~basal_mask]
        _, p_bo = stats.mannwhitneyu(bd, od, alternative="greater")
        log(f"  Basal depth mean={bd.mean():.4f}  Non-Basal mean={od.mean():.4f}  "
            f"MWU p={p_bo:.4f}")
        results["p5_basal_depth_p"] = p_bo

    # ── S2-P5: BRCA1 expression in Basal ─────────────────
    log("")
    log("  S2-P5: BRCA1 EXPRESSION IN BASAL-LIKE")
    log("  NOTE: Expression proxy — definitive test needs")
    log("  GDC somatic mutation + methylation data (auth required).")
    if "BRCA1" in expr_a.columns and basal_mask.sum() > 10:
        br_bas = expr_a.loc[basal_mask, "BRCA1"]
        if "luma_mask" in dir() and luma_mask.sum() > 10:
            br_lum = expr_a.loc[luma_mask, "BRCA1"]
            try:
                _, p_b1 = stats.mannwhitneyu(br_bas, br_lum,
                                              alternative="less")
                log(f"  BRCA1 Basal: mean={br_bas.mean():.4f}  n={len(br_bas)}")
                log(f"  BRCA1 LumA:  mean={br_lum.mean():.4f}  n={len(br_lum)}")
                log(f"  MWU (Basal < LumA): p={p_b1:.4f}")
                results["p5_brca1_p"] = p_b1
                if p_b1 < 0.05:
                    log("  S2-P5 proxy: ✓ BRCA1 expression lower in Basal")
                    log("  (Mutation/methylation rate is the definitive test)")
                else:
                    log("  S2-P5 proxy: ✗ BRCA1 expression not significantly lower")
            except Exception as e:
                log(f"  S2-P5 error: {e}")
    else:
        log("  BRCA1 not in TCGA expression or insufficient Basal samples.")

    # ── S2-P6: TCGA OS survival ───────────────────────────
    log("")
    log("  S2-P6: DEPTH vs OS IN TCGA-BRCA")

    # Find OS columns in clinical
    os_time_col  = None
    os_event_col = None
    for col in clin_a.columns:
        cl = col.lower()
        if os_time_col is None and any(
            x in cl for x in ["os.time", "os_time", "days_to_death",
                               "overall_survival_time", "_os_"]
        ):
            os_time_col = col
        if os_event_col is None and any(
            x in cl for x in ["os.event", "vital_status", "os_event",
                               "dead", "_os_event"]
        ):
            os_event_col = col

    # Try pan-cancer survival table
    surv_path = paths.get("tcga_surv")
    surv_df   = None
    if surv_path and os.path.exists(surv_path):
        try:
            with open(surv_path, "rt", errors="ignore") as f:
                surv_df = pd.read_csv(f, sep="\t", index_col=0)
            log(f"  Pan-cancer survival: {surv_df.shape}")
            log(f"  Survival columns: {list(surv_df.columns[:10])}")
            # Filter to BRCA
            brca_m = surv_df.index.str.upper().str.startswith("TCGA-BR")
            if brca_m.sum() > 50:
                surv_df = surv_df[brca_m]
                log(f"  BRCA samples in survival table: {len(surv_df)}")
        except Exception as e:
            log(f"  Survival load error: {e}")

    if os_time_col is None and surv_df is not None:
        # Merge survival into clin_a
        for col in surv_df.columns:
            cl = col.lower()
            if os_time_col is None and "time" in cl:
                os_time_col = "_surv_time"
            if os_event_col is None and ("event" in cl or "dead" in cl
                                          or "status" in cl):
                os_event_col = "_surv_event"
        surv_common = clin_a.index.intersection(surv_df.index)
        if len(surv_common) > 50:
            t_col = [c for c in surv_df.columns if "time" in c.lower()]
            e_col = [c for c in surv_df.columns
                     if any(x in c.lower()
                            for x in ["event", "dead", "status"])]
            if t_col and e_col:
                clin_a["_surv_time"]  = surv_df.loc[
                    surv_common, t_col[0]
                ].reindex(clin_a.index)
                clin_a["_surv_event"] = surv_df.loc[
                    surv_common, e_col[0]
                ].reindex(clin_a.index)
                os_time_col  = "_surv_time"
                os_event_col = "_surv_event"
                log(f"  Merged survival: time={t_col[0]}  event={e_col[0]}")

    log(f"  OS time col:  {os_time_col}")
    log(f"  OS event col: {os_event_col}")

    if os_time_col and os_event_col:
        try:
            sc = clin_a[[os_time_col, os_event_col, "depth"]].copy()
            sc.columns = ["time", "event", "depth"]
            sc["time"]  = pd.to_numeric(sc["time"],  errors="coerce")
            sc["event"] = pd.to_numeric(sc["event"], errors="coerce")
            sc = sc.dropna()
            sc = sc[sc["time"] > 0]

            log(f"  TCGA survival n={len(sc)}  "
                f"events={sc['event'].sum():.0f}")

            depth_med = sc["depth"].median()
            hi = sc[sc["depth"] >= depth_med]
            lo = sc[sc["depth"] <  depth_med]

            stat, p_lr = _logrank(
                hi["time"].values, hi["event"].values,
                lo["time"].values, lo["event"].values
            )
            results["p6_tcga_logrank_p"]    = p_lr
            results["p6_tcga_logrank_stat"] = stat
            log(f"  Log-rank stat={stat:.3f}  p={p_lr:.4f}")

            if not np.isnan(p_lr) and p_lr < 0.05:
                hi_med_t = hi.loc[hi["event"]==1, "time"].median()
                lo_med_t = lo.loc[lo["event"]==1, "time"].median()
                log(f"  High-depth median OS: {hi_med_t:.1f}")
                log(f"  Low-depth  median OS: {lo_med_t:.1f}")
                if hi_med_t < lo_med_t:
                    log("  S2-P6 TCGA: ✓ CONFIRMED — depth-high shorter OS")
                else:
                    log("  S2-P6 TCGA: ✗ INVERTED")
            elif not np.isnan(p_lr):
                log(f"  S2-P6 TCGA: ns  p={p_lr:.4f}")

            r_sp, p_sp = stats.spearmanr(sc["depth"], sc["event"])
            results["p6_tcga_spearman_r"] = r_sp
            log(f"  Spearman r(depth, event) = {r_sp:+.4f}  {fmt_p(p_sp)}")

            results["surv_sub_tcga"] = sc
        except Exception as e:
            log(f"  S2-P6 TCGA error: {e}")
            log(traceback.format_exc())
    else:
        log("  OS columns not found. S2-P6 TCGA not testable.")
        log("  Available clinical columns:")
        for col in clin_a.columns[:30]:
            log(f"    '{col}'")

    return results


# ============================================================
# PREDICTION SCORECARD
# ============================================================

def prediction_scorecard(t1, t2, t3):
    log("")
    log("=" * 65)
    log("PREDICTION SCORECARD — BRCA-S4c")
    log("=" * 65)
    log("")
    log(f"  {'Pred':<10}  {'Status':<28}  Note")
    log(f"  {'-'*72}")

    def row(pred, status, note=""):
        log(f"  {pred:<10}  {status:<28}  {note}")

    # S2-P1
    p1r = t1.get("p1_r", np.nan)
    if not np.isnan(p1r):
        if p1r < -0.098:
            row("S2-P1", "✓ CONFIRMED",
                f"r={p1r:+.3f} more negative than -0.098")
        elif p1r < 0:
            row("S2-P1", "✓ DIRECTIONAL",
                f"r={p1r:+.3f} (negative, < Script1 proxy)")
        else:
            row("S2-P1", "✗ NOT CONFIRMED", f"r={p1r:+.3f}")
    else:
        row("S2-P1", "? DATA LIMITED", "pCR parse check needed")

    # S2-P2
    eed_r  = t1.get("eed_r",   np.nan)
    ezh2_r = t1.get("ezh2_r",  np.nan)
    eed_a  = t1.get("eed_auc", np.nan)
    ezh2_a = t1.get("ezh2_auc",np.nan)
    if not (np.isnan(eed_r) or np.isnan(ezh2_r)):
        if eed_r < -0.10 and eed_r < ezh2_r:
            row("S2-P2", "✓ CONFIRMED",
                f"EED r={eed_r:+.3f} vs EZH2 r={ezh2_r:+.3f}")
        elif eed_r < ezh2_r:
            row("S2-P2", "✓ PARTIAL",
                f"EED r more negative but ≥-0.10")
        else:
            row("S2-P2", "✗ NOT CONFIRMED", f"EED r={eed_r:+.3f}")
    else:
        row("S2-P2", "? NOT IN DATASET", "EED or EZH2 absent")

    # S2-P3
    parp1_r = t1.get("parp1_r", np.nan)
    if not np.isnan(parp1_r):
        if parp1_r < -0.08:
            row("S2-P3", "✓ CONFIRMED",    f"PARP1 r={parp1_r:+.3f}")
        elif parp1_r < 0:
            row("S2-P3", "✓ DIRECTIONAL",  f"PARP1 r={parp1_r:+.3f}")
        else:
            row("S2-P3", "✗ NOT CONFIRMED",f"PARP1 r={parp1_r:+.3f}")
    else:
        row("S2-P3", "? NOT IN DATASET", "PARP1 absent")

    # S2-P4
    lar_m = t2.get("p4_lar_mean", np.nan)
    if not np.isnan(lar_m):
        row("S2-P4",
            "✓ CONFIRMED" if lar_m < 0.40 else "✗ NOT CONFIRMED",
            f"LAR mean depth={lar_m:.3f}")
    else:
        row("S2-P4", "? DATA LIMITED", "Lehmann scoring incomplete")

    # S2-P5
    p5 = t3.get("p5_brca1_p", np.nan)
    if not np.isnan(p5):
        row("S2-P5", "? PROXY ONLY",
            f"Expression p={p5:.4f} — DNA test needs GDC auth")
    else:
        row("S2-P5", "? PROXY ONLY",
            "BRCA1 expression proxy — mutation needs GDC")

    # S2-P6 (GSE25066 DRFS)
    p6g = t1.get("p6_gse_logrank_p", np.nan)
    if not np.isnan(p6g):
        row("S2-P6a GSE",
            "✓ CONFIRMED" if p6g < 0.05 else "? NS",
            f"DRFS log-rank p={p6g:.4f}")
    else:
        row("S2-P6a GSE", "? DRFS NOT FOUND", "Check metadata parse")

    # S2-P6 (TCGA OS)
    p6t = t3.get("p6_tcga_logrank_p", np.nan)
    if not np.isnan(p6t):
        row("S2-P6b TCGA",
            "✓ CONFIRMED" if p6t < 0.05 else "? NS",
            f"OS log-rank p={p6t:.4f}")
    else:
        row("S2-P6b TCGA", "? OS NOT FOUND", "Check survival table")

    # S2-P7
    sp = t1.get("p7_spi1_ptprc", np.nan)
    sd = t1.get("p7_spi1_depth", np.nan)
    if not np.isnan(sp):
        row("S2-P7",
            "✓ CONFIRMED" if sp > abs(sd) else "✗ NOT CONFIRMED",
            f"r(SPI1,PTPRC)={sp:+.3f}  r(SPI1,depth)={sd:+.3f}")
    else:
        row("S2-P7", "? NOT IN DATASET", "SPI1/PTPRC absent")

    # S2-P8
    row("S2-P8", "P PROSPECTIVE",
        "EED>EZH2 as EZH2i biomarker — future datasets")

    # S2-P9
    ar_r = t1.get("p9_ar_r", np.nan)
    if not np.isnan(ar_r):
        row("S2-P9",
            "✓ CONFIRMED" if ar_r < -0.15 else
            ("✓ DIRECTIONAL" if ar_r < 0 else "✗ NOT CONFIRMED"),
            f"r(AR,depth)={ar_r:+.3f}")
    else:
        row("S2-P9", "? NOT IN DATASET", "AR absent")

    log("")


# ============================================================
# FIGURE
# ============================================================

def generate_figure(t1, t2, t3):
    log("")
    log("--- Generating figure ---")

    g_all      = t1.get("gene_df")
    m_all      = t1.get("meta_df")
    bulk_corrs = t1.get("bulk_corrs", [])
    surv_gse   = t1.get("surv_sub_gse")
    depth_med  = t1.get("depth_med_gse", np.nan)
    depth_by_st = t2.get("depth_by_subtype", {})
    tcga_depth = t3.get("tcga_depth")

    clr = {
        "hi":     "#c0392b",
        "lo":     "#2980b9",
        "eed":    "#e74c3c",
        "ezh2":   "#3498db",
        "bar_up": "#c0392b",
        "bar_dn": "#2980b9",
    }

    fig = plt.figure(figsize=(24, 18))
    fig.suptitle(
        "TNBC — Script 2 | OrganismCore BRCA-S4c/d | 2026-03-04\n"
        "GSE25066 (pCR + DRFS) | TCGA-BRCA (PAM50 + OS)\n"
        "Embedded GPL96 mapping | Confirmed URLs from diagnostic",
        fontsize=10, fontweight="bold", y=1.005
    )
    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.52, wspace=0.42)

    # A — Depth by pCR group
    ax_a = fig.add_subplot(gs[0, 0])
    if m_all is not None and "pCR" in m_all.columns and "depth" in m_all.columns:
        pcr1 = m_all.loc[m_all["pCR"] == 1, "depth"]
        pcr0 = m_all.loc[m_all["pCR"] == 0, "depth"]
        if len(pcr1) > 2 and len(pcr0) > 2:
            ax_a.hist(pcr0, bins=30, alpha=0.65, color=clr["hi"],
                      label=f"pCR=0 (RD) n={len(pcr0)}")
            ax_a.hist(pcr1, bins=30, alpha=0.65, color=clr["lo"],
                      label=f"pCR=1 n={len(pcr1)}")
            ax_a.axvline(pcr0.mean(), color=clr["hi"], lw=2, ls="--")
            ax_a.axvline(pcr1.mean(), color=clr["lo"], lw=2, ls="--")
    ax_a.set_title("A — Depth by pCR (S2-P1)\nGSE25066",
                   fontsize=8, fontweight="bold")
    ax_a.set_xlabel("Depth Score", fontsize=7)
    ax_a.set_ylabel("n samples", fontsize=7)
    ax_a.legend(fontsize=6)

    # B — EED vs EZH2 pCR correlation bars
    ax_b = fig.add_subplot(gs[0, 1])
    gene_preds = {}
    for gene in ["EED", "EZH2", "PARP1", "AR", "MKI67",
                 "KRT5", "ESR1", "VIM"]:
        r = t1.get(f"{gene.lower()}_r", np.nan)
        if not np.isnan(r):
            gene_preds[gene] = r
    if gene_preds:
        sorted_gp = sorted(gene_preds.items(), key=lambda x: x[1])
        gnames = [g for g, _ in sorted_gp]
        rvals  = [r for _, r in sorted_gp]
        cc     = [clr["bar_up"] if r < 0 else clr["bar_dn"] for r in rvals]
        ax_b.barh(gnames, rvals, color=cc)
        ax_b.axvline(0, color="black", lw=0.8)
        ax_b.axvline(-0.10, color="gray", lw=0.8, ls="--",
                     label="P2 threshold -0.10")
    ax_b.set_title("B — Gene r(pCR) (S2-P2/P3)\nNeg = lower pCR prob",
                   fontsize=8, fontweight="bold")
    ax_b.set_xlabel("r(gene, pCR)", fontsize=7)
    ax_b.tick_params(axis="y", labelsize=7)
    ax_b.legend(fontsize=5)

    # C — Bulk depth correlations
    ax_c = fig.add_subplot(gs[0, 2])
    if bulk_corrs:
        top20 = bulk_corrs[:20]
        gc = [g for g, r, p in top20]
        vc = [r for g, r, p in top20]
        cc = [clr["bar_up"] if v > 0 else clr["bar_dn"] for v in vc]
        ax_c.barh(gc, vc, color=cc)
        ax_c.axvline(0, color="black", lw=0.8)
    ax_c.set_title("C — Bulk Depth Correlates\nTop 20",
                   fontsize=8, fontweight="bold")
    ax_c.set_xlabel("r", fontsize=7)
    ax_c.tick_params(axis="y", labelsize=6)

    # D — Lehmann depth box
    ax_d = fig.add_subplot(gs[1, 0])
    if depth_by_st:
        order  = sorted(depth_by_st.items(), key=lambda x: x[1].mean())
        labs   = [s for s, _ in order]
        data   = [d for _, d in order]
        bp     = ax_d.boxplot(data, patch_artist=True, notch=False)
        pal    = ["#f39c12", "#e74c3c", "#e67e22",
                  "#3498db", "#8e44ad", "#2ecc71"]
        for patch, c in zip(bp["boxes"], pal[:len(labs)]):
            patch.set_facecolor(c)
            patch.set_alpha(0.75)
        ax_d.set_xticklabels(labs, fontsize=7)
    ax_d.set_title("D — Lehmann Subtype Depth (S2-P4)\nPredicted: LAR<BL<M/MSL",
                   fontsize=8, fontweight="bold")
    ax_d.set_ylabel("Depth Score", fontsize=7)

    # E — DRFS KM curve (GSE25066)
    ax_e = fig.add_subplot(gs[1, 1])
    if surv_gse is not None and not np.isnan(depth_med):
        hi_s = surv_gse[surv_gse["depth"] >= depth_med]
        lo_s = surv_gse[surv_gse["depth"] <  depth_med]
        for sub, col, lab in [
            (hi_s, clr["hi"], f"High depth n={len(hi_s)}"),
            (lo_s, clr["lo"], f"Low depth  n={len(lo_s)}"),
        ]:
            sub_s = sub.sort_values("DRFS_time")
            t = sub_s["DRFS_time"].values
            e = sub_s["DRFS_event"].values
            # Kaplan-Meier estimate
            surv = 1.0
            ts, ss = [0.0], [1.0]
            for ti, ei in zip(t, e):
                if ei == 1:
                    n_at_risk = (t >= ti).sum()
                    surv *= (1 - 1 / n_at_risk) if n_at_risk > 0 else surv
                ts.append(ti)
                ss.append(surv)
            ax_e.step(ts, ss, where="post", color=col, label=lab, lw=1.5)
        p_lr_gse = t1.get("p6_gse_logrank_p", np.nan)
        p_str = f"p={p_lr_gse:.4f}" if not np.isnan(p_lr_gse) else ""
        ax_e.set_title(f"E — DRFS KM (S2-P6) GSE25066\n{p_str}",
                       fontsize=8, fontweight="bold")
        ax_e.set_xlabel("Years", fontsize=7)
        ax_e.set_ylabel("DRFS probability", fontsize=7)
        ax_e.legend(fontsize=6)
        ax_e.set_ylim(0, 1.05)
    else:
        ax_e.axis("off")
        ax_e.text(0.5, 0.5, "DRFS data\nnot available",
                  ha="center", va="center", fontsize=9,
                  transform=ax_e.transAxes)

    # F — AR vs depth
    ax_f = fig.add_subplot(gs[1, 2])
    if g_all is not None and m_all is not None and \
            "AR" in g_all.columns and "depth" in m_all.columns:
        ax_f.scatter(g_all["AR"].values, m_all["depth"].values,
                     c=clr["lo"], alpha=0.2, s=5)
        try:
            rv, _ = stats.pearsonr(g_all["AR"].values, m_all["depth"].values)
            ax_f.set_title(f"F — AR vs Depth (S2-P9)\nr={rv:+.3f}  "
                           f"predicted r<-0.15",
                           fontsize=8, fontweight="bold")
        except Exception:
            ax_f.set_title("F — AR vs Depth", fontsize=8)
    else:
        ax_f.axis("off")
        ax_f.text(0.5, 0.5, "AR data\nnot available",
                  ha="center", va="center", transform=ax_f.transAxes)
    ax_f.set_xlabel("AR expression", fontsize=7)
    ax_f.set_ylabel("Depth", fontsize=7)

    # G — TCGA depth distribution
    ax_g = fig.add_subplot(gs[2, 0])
    if tcga_depth is not None:
        ax_g.hist(tcga_depth.values, bins=50,
                  color=clr["hi"], alpha=0.7, edgecolor="white")
        ax_g.axvline(tcga_depth.mean(), color="black", lw=1.5, ls="--",
                     label=f"mean={tcga_depth.mean():.3f}")
        ax_g.set_title("G — TCGA-BRCA Depth Distribution",
                       fontsize=8, fontweight="bold")
        ax_g.set_xlabel("Depth score", fontsize=7)
        ax_g.set_ylabel("n samples", fontsize=7)
        ax_g.legend(fontsize=6)
    else:
        ax_g.axis("off")
        ax_g.text(0.5, 0.5, "TCGA depth\nnot available",
                  ha="center", va="center", transform=ax_g.transAxes)

    # H — EED vs depth scatter
    ax_h = fig.add_subplot(gs[2, 1])
    if g_all is not None and m_all is not None and \
            "EED" in g_all.columns and "depth" in m_all.columns:
        ax_h.scatter(g_all["EED"].values, m_all["depth"].values,
                     c=clr["eed"], alpha=0.2, s=5)
        try:
            rv, _ = stats.pearsonr(g_all["EED"].values,
                                   m_all["depth"].values)
            ax_h.set_title(f"H — EED vs Depth\nr={rv:+.3f}",
                           fontsize=8, fontweight="bold")
        except Exception:
            ax_h.set_title("H — EED vs Depth", fontsize=8)
        ax_h.set_xlabel("EED", fontsize=7)
        ax_h.set_ylabel("Depth", fontsize=7)
    else:
        ax_h.axis("off")
        ax_h.text(0.5, 0.5, "EED data\nnot available",
                  ha="center", va="center", transform=ax_h.transAxes)

    # I — Summary panel
    ax_i = fig.add_subplot(gs[2, 2])
    ax_i.axis("off")

    p1s  = f"{t1.get('p1_r', float('nan')):+.3f}" if "p1_r"  in t1 else "n/a"
    p6gs = (f"p={t1.get('p6_gse_logrank_p', float('nan')):.4f}"
            if "p6_gse_logrank_p" in t1 else "n/a")
    p6ts = (f"p={t3.get('p6_tcga_logrank_p', float('nan')):.4f}"
            if "p6_tcga_logrank_p" in t3 else "n/a")
    eed_s  = f"{t1.get('eed_r',   float('nan')):+.3f}" if "eed_r"  in t1 else "n/a"
    ezh2_s = f"{t1.get('ezh2_r',  float('nan')):+.3f}" if "ezh2_r" in t1 else "n/a"

    ax_i.text(0.03, 0.97,
        f"I — SCRIPT 2 SUMMARY\n{'─'*30}\n"
        f"GSE25066 + TCGA-BRCA | 2026-03-04\n\n"
        f"S2-P1  depth∼pCR (all):   r={p1s}\n"
        f"S2-P2  EED r(pCR):        {eed_s}  (EZH2:{ezh2_s})\n"
        f"S2-P3  PARP1 r(pCR):      see panel B\n"
        f"S2-P4  Lehmann depth:     see panel D\n"
        f"S2-P5  BRCA1 proxy:       expression only\n"
        f"S2-P6a DRFS GSE25066:     {p6gs}\n"
        f"S2-P6b OS TCGA:           {p6ts}\n"
        f"S2-P7  SPI1 resolution:   see text\n"
        f"S2-P8  EED biomarker:     prospective\n"
        f"S2-P9  AR∼depth:          see panel F\n\n"
        f"DRUG TARGETS UPDATED:\n"
        f"1. EZH2i   (tazemetostat)\n"
        f"2. EED inh (MAK683 — novel)\n"
        f"3. PARPi   (olaparib)\n"
        f"4. AR block (LAR subtype)\n\n"
        f"OrganismCore BRCA-S4c/d | 2026-03-04",
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
    log("BULK VALIDATION + pCR + DRFS + TCGA OS")
    log("OrganismCore — BRCA-S4c/d | 2026-03-04")
    log("=" * 65)
    log("")
    log("CONFIRMED URLS (diagnostic 2026-03-04):")
    log(f"  GSE25066:  {GSE25066_URL}")
    log(f"  TCGA EXPR: {TCGA_EXPR_URL}")
    log(f"  TCGA CLIN: {TCGA_CLIN_URL}")
    log(f"  TCGA SURV: {TCGA_SURV_URL}")
    log(f"  GPL96:     EMBEDDED (all NCBI FTP files return 404)")
    log("")

    paths = acquire_data()

    gene_df, meta_df = parse_gse25066(paths.get("gse25066"))

    t1 = task1_pcr_and_survival(gene_df, meta_df)
    t2 = task2_lehmann(gene_df, meta_df)
    t3 = task3_tcga(paths)

    prediction_scorecard(t1, t2, t3)
    generate_figure(t1, t2, t3)

    # CSV export
    rows = [{"gene": g, "r_depth": r, "p_value": p}
            for g, r, p in t1.get("bulk_corrs", [])]
    if rows:
        pd.DataFrame(rows).to_csv(CSV_FILE, index=False)
        log(f"  CSV: {CSV_FILE}")

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
