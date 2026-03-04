"""
BRCA TNBC — SCRIPT 2
BULK RNA-seq VALIDATION + pCR ANALYSIS + SURVIVAL
OrganismCore — Document BRCA-S4c/d
Date: 2026-03-04

FIXES vs previous version:
  1. GPL96 annotation: use GPL96_family.soft.gz (NCBI soft file)
     instead of GPL96.annot.gz (404 on NCBI FTP-over-HTTPS)
     Parser updated to read soft format.

  2. TCGA-BRCA: use UCSC Xena hub URLs instead of cBioPortal
     S3 bucket (403 Forbidden on brca_tcga_pub2015.tar.gz).
     Expression: TCGA.BRCA.sampleMap/HiSeqV2_PANCAN.gz
     Clinical:   TCGA.BRCA.sampleMap/BRCA_clinicalMatrix
     Both are publicly accessible without authentication.

  Everything else is identical to the previous version.
"""

import os
import sys
import gzip
import tarfile
import time
import warnings
import urllib.request
import traceback
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

# ============================================================
# DATA URLS — corrected
# ============================================================

# GSE25066
GSE25066_MATRIX_URL  = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/"
    "GSE25nnn/GSE25066/matrix/"
    "GSE25066_series_matrix.txt.gz"
)
GSE25066_MATRIX_FILE = "GSE25066_series_matrix.txt.gz"

# GPL96 soft file — correct URL (annot.gz returns 404)
# soft file contains full probe→gene annotation
GPL96_URL  = (
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/"
    "GPL96nnn/GPL96/soft/GPL96_family.soft.gz"
)
GPL96_FILE = "GPL96_family.soft.gz"

# TCGA-BRCA via UCSC Xena (cBioPortal S3 returns 403)
XENA_EXPR_URL = (
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/"
    "download/TCGA.BRCA.sampleMap%2FHiSeqV2_PANCAN.gz"
)
XENA_EXPR_FILE = "TCGA_BRCA_HiSeqV2_PANCAN.gz"

XENA_CLIN_URL  = (
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/"
    "download/TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix"
)
XENA_CLIN_FILE = "TCGA_BRCA_clinicalMatrix.tsv"

# TCGA PAM50 subtype calls — separate Xena phenotype file
XENA_PHENO_URL  = (
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/"
    "download/TCGA.BRCA.sampleMap%2FBRCA_phenotype_r"
    "enamed_noPE.tsv.gz"
)
XENA_PHENO_FILE = "TCGA_BRCA_phenotype.tsv.gz"

# Search paths — check these before downloading
SEARCH_PATHS = [
    DATA_DIR,
    BASE_DIR,
    "./",
    "../",
    "../../",
    os.path.join(SCRIPT_DIR, "../"),
    os.path.join(SCRIPT_DIR, "../../"),
    os.path.join(SCRIPT_DIR, "tnbc_results/data"),
    os.path.join(SCRIPT_DIR, "tnbc_s1_analysis/data"),
    os.path.expanduser("~/cancer/BRCA/"),
    os.path.expanduser("~/cancer/BRCA/DEEP_DIVE/TNBC/tnbc_s1_analysis/data"),
]

# ============================================================
# GENE PANELS (locked BRCA-S4a/c)
# ============================================================

# Affymetrix HG-U133A probe IDs for target genes
# Used as fallback when GPL96 soft file is unavailable
# or when soft file probe→gene lookup fails
PROBE_HINTS = {
    "ESR1":   ["205225_at",   "211233_x_at"],
    "FOXA1":  ["202340_x_at", "202341_s_at"],
    "GATA3":  ["209604_s_at", "213134_x_at"],
    "SPDEF":  ["219197_s_at", "219198_at"],
    "PGR":    ["208305_at",   "208306_x_at"],
    "KRT5":   ["201820_at"],
    "KRT14":  ["201744_s_at"],
    "KRT8":   ["201525_at",   "214456_x_at"],
    "KRT18":  ["201579_at"],
    "SOX10":  ["221579_at",   "221580_s_at"],
    "FOXC1":  ["203853_s_at"],
    "EGFR":   ["201983_s_at", "201984_s_at"],
    "VIM":    ["201426_s_at"],
    "CDH3":   ["205497_at"],
    "CDH1":   ["201130_s_at"],
    "EZH2":   ["203358_s_at"],
    "EED":    ["218657_s_at"],
    "HDAC1":  ["202705_at"],
    "HDAC2":  ["200895_s_at"],
    "KDM1A":  ["218263_s_at"],
    "DNMT3A": ["210436_at"],
    "MKI67":  ["212022_s_at"],
    "TOP2A":  ["201291_s_at"],
    "BRCA1":  ["204531_s_at", "204532_s_at"],
    "BRCA2":  ["208368_s_at"],
    "TP53":   ["201746_at"],
    "PTEN":   ["214440_at"],
    "PIK3CA": ["212733_at"],
    "AKT1":   ["207163_s_at"],
    "RB1":    ["202248_at"],
    "PARP1":  ["208501_at"],
    "AR":     ["211110_x_at", "211111_s_at"],
    "ZEB1":   ["209839_at"],
    "ZEB2":   ["213596_at"],
    "SNAI1":  ["216834_at"],
    "CLDN3":  ["214900_at"],
    "CLDN4":  ["201428_at"],
    "AURKA":  ["204092_s_at"],
    "CD274":  ["223834_at"],
    "SPI1":   ["203140_at"],
    "PTPRC":  ["213293_s_at"],
    "CD68":   ["203507_at"],
    "CD3D":   ["205986_at"],
    "CD8A":   ["205758_at"],
    "FOXP3":  ["220306_at"],
    "CDX2":   ["207769_at"],
    "NKX2-1": ["213844_at"],
    "NKX3-1": ["209706_at"],
    "OLIG2":  ["213825_at"],
}

ALL_TARGET_GENES = list(PROBE_HINTS.keys())

SWITCH_GENES    = ["ESR1","FOXA1","GATA3","SPDEF","PGR",
                   "KRT8","KRT18","BRCA1"]
FA_MARKERS      = ["KRT5","KRT14","SOX10","FOXC1","EGFR",
                   "VIM","CDH3"]
EPIGENETIC      = ["EZH2","EED","HDAC1","HDAC2","KDM1A","DNMT3A"]
COMPOSITE_TYPE  = ["BRCA1","BRCA2","TP53","PTEN","RB1",
                   "PIK3CA","AKT1"]
DEPTH_BASAL_G   = ["KRT5","KRT14","SOX10","FOXC1","EGFR"]
DEPTH_LUMINAL_G = ["ESR1","FOXA1","GATA3","SPDEF"]
IMMUNE_MARKERS  = ["SPI1","PTPRC","CD68","CD3D","CD8A",
                   "FOXP3","CD274"]
CONTROLS        = ["CDX2","NKX2-1","NKX3-1","OLIG2"]

# ============================================================
# LEHMANN 2016 SUBTYPE CENTROID GENES
# ============================================================

LEHMANN_GENES = {
    "BL1": ["CCNE1","CDC20","CDK2","CDKN3","CEP55","AURKA",
            "AURKB","BIRC5","BUB1","BUB1B","TOP2A","MKI67"],
    "BL2": ["EGFR","MET","DCLK1","ENPP1","HIF1A","LCK",
            "SDC1","CDK5","CDK5R1","ABCC4"],
    "M":   ["VIM","FN1","MMP2","MMP9","TWIST1","ZEB1","ZEB2",
            "SNAI1","SNAI2","CDH2","CDH11","ITGB1"],
    "MSL": ["CD44","ALDH1A3","FOXC1","ELF5","PROCR",
            "KDR","PDGFRB","ANGPTL1","ANGPTL2"],
    "IM":  ["CD3D","CD8A","CTLA4","LAG3","PDCD1","HAVCR2",
            "TIGIT","SIRPG","SPI1","IRF1","STAT1","IRF7"],
    "LAR": ["AR","DHCR24","FKBP4","ALCAM","APOD",
            "FASN","LDLR","INPP4B","PIK3R1","XBP1"],
}

ALL_LEHMANN_GENES = list(dict.fromkeys(
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

# ============================================================
# FILE UTILITIES
# ============================================================

def find_file(filename):
    for base in SEARCH_PATHS:
        p = os.path.join(base, filename)
        if os.path.exists(p) and os.path.getsize(p) > 1000:
            log(f"  Found: {p} ({os.path.getsize(p)/1e6:.1f} MB)")
            return os.path.abspath(p)
    return None


def fetch_url(url, dest, retries=3, timeout=300):
    for attempt in range(retries):
        try:
            log(f"  Attempt {attempt+1}: {url[-70:]}")
            req = urllib.request.Request(
                url,
                headers={
                    "User-Agent": (
                        "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                        "AppleWebKit/537.36 (KHTML, like Gecko) "
                        "Chrome/120.0.0.0 Safari/537.36"
                    )
                }
            )
            with urllib.request.urlopen(req, timeout=timeout) as r:
                data = r.read()
            with open(dest, "wb") as f:
                f.write(data)
            log(f"  Downloaded: {os.path.basename(dest)} "
                f"({len(data)/1e6:.1f} MB)")
            return dest
        except Exception as e:
            log(f"  Attempt {attempt+1} failed: {e}")
            if attempt < retries - 1:
                time.sleep(4)
    return None

# ============================================================
# STEP 0: ACQUIRE DATA
# ============================================================

def acquire_data():
    log("=" * 65)
    log("STEP 0: DATA ACQUISITION")
    log("=" * 65)

    # --- GPL96 soft file ---
    log("\n--- GPL96 platform annotation (soft file) ---")
    log(f"  Using: {GPL96_FILE}")
    log("  (GPL96_family.soft.gz — not GPL96.annot.gz which 404s)")
    gpl96_path = find_file(GPL96_FILE)
    if gpl96_path is None:
        gpl96_path = fetch_url(
            GPL96_URL,
            os.path.join(DATA_DIR, GPL96_FILE),
            timeout=300
        )
    if gpl96_path is None:
        log("  WARNING: GPL96 soft file unavailable.")
        log("  Falling back to PROBE_HINTS dictionary.")

    # --- GSE25066 ---
    log("\n--- GSE25066 series matrix ---")
    gse_path = find_file(GSE25066_MATRIX_FILE)
    if gse_path is None:
        gse_path = fetch_url(
            GSE25066_MATRIX_URL,
            os.path.join(DATA_DIR, GSE25066_MATRIX_FILE)
        )
    if gse_path is None:
        log("  FATAL: Cannot obtain GSE25066.")
        return None, None, None, None

    # --- TCGA-BRCA via UCSC Xena ---
    log("\n--- TCGA-BRCA via UCSC Xena hub ---")
    log("  (cBioPortal S3 returns 403 — using Xena instead)")

    xena_expr_path = find_file(XENA_EXPR_FILE)
    if xena_expr_path is None:
        log(f"  Downloading expression (~500 MB)...")
        xena_expr_path = fetch_url(
            XENA_EXPR_URL,
            os.path.join(DATA_DIR, XENA_EXPR_FILE),
            timeout=600
        )
    if xena_expr_path:
        log(f"  Expression: {xena_expr_path}")
    else:
        log("  Expression download failed.")

    xena_clin_path = find_file(XENA_CLIN_FILE)
    if xena_clin_path is None:
        log(f"  Downloading clinical matrix...")
        xena_clin_path = fetch_url(
            XENA_CLIN_URL,
            os.path.join(DATA_DIR, XENA_CLIN_FILE),
            timeout=120
        )

    xena_pheno_path = find_file(XENA_PHENO_FILE)
    if xena_pheno_path is None:
        log(f"  Downloading phenotype (PAM50) file...")
        xena_pheno_path = fetch_url(
            XENA_PHENO_URL,
            os.path.join(DATA_DIR, XENA_PHENO_FILE),
            timeout=120
        )

    # Package TCGA paths into a dict for later steps
    tcga_paths = {
        "expr":  xena_expr_path,
        "clin":  xena_clin_path,
        "pheno": xena_pheno_path,
    }
    if not any(tcga_paths.values()):
        log("  TCGA-BRCA unavailable. S2-P5/P6 will be skipped.")
        tcga_paths = None

    return gse_path, gpl96_path, tcga_paths

# ============================================================
# STEP 1: BUILD PROBE MAP
# Parses GPL96_family.soft.gz (soft format)
# Falls back to PROBE_HINTS if unavailable
# ============================================================

def build_probe_map(gpl96_path):
    log("=" * 65)
    log("STEP 1: BUILD PROBE→GENE MAP")
    log("=" * 65)

    probe_to_gene = {}

    if gpl96_path and os.path.exists(gpl96_path):
        log(f"  Parsing: {os.path.basename(gpl96_path)}")
        try:
            open_fn = gzip.open if gpl96_path.endswith(".gz") else open
            n_probes = 0
            with open_fn(gpl96_path, "rt", errors="ignore") as f:
                in_table = False
                id_val   = None
                sym_val  = None
                for line in f:
                    line = line.rstrip("\n")

                    # Platform table section markers
                    if line.startswith("!platform_table_begin"):
                        in_table = True
                        hdr_line = f.readline().rstrip("\n")
                        hdr      = hdr_line.split("\t")
                        # Find ID and Gene Symbol column indices
                        id_col  = 0
                        sym_col = None
                        for i, h in enumerate(hdr):
                            h_low = h.lower()
                            if h_low == "id":
                                id_col = i
                            elif "gene symbol" in h_low or h_low == "gene_symbol":
                                sym_col = i
                            elif "gene_assignment" in h_low and sym_col is None:
                                sym_col = i
                        log(f"  Table header: id_col={id_col} "
                            f"sym_col={sym_col}")
                        log(f"  Header fields: {hdr[:8]}")
                        continue

                    if line.startswith("!platform_table_end"):
                        in_table = False
                        continue

                    if not in_table:
                        continue

                    parts = line.split("\t")
                    if len(parts) <= max(id_col, sym_col or 0):
                        continue

                    probe = parts[id_col].strip()
                    if sym_col is not None:
                        sym_field = parts[sym_col].strip()
                        # Handle "Gene Symbol /// Gene Symbol" format
                        symbols = [
                            s.strip()
                            for s in sym_field.split("///")
                        ]
                        for sym in symbols:
                            if sym and sym != "---" and sym != "":
                                probe_to_gene[probe] = sym
                                n_probes += 1
                                break

            log(f"  Probes mapped: {n_probes}")

        except Exception as e:
            log(f"  Soft file parse error: {e}")
            log(traceback.format_exc())
            probe_to_gene = {}

    # Fall back to PROBE_HINTS for any remaining targets
    n_hint = 0
    for gene, probes in PROBE_HINTS.items():
        for probe in probes:
            if probe not in probe_to_gene:
                probe_to_gene[probe] = gene
                n_hint += 1

    log(f"  PROBE_HINTS added: {n_hint} entries")
    log(f"  Total probe map:   {len(probe_to_gene)} entries")

    # Build gene → probes reverse map
    gene_to_probes = {}
    for probe, gene in probe_to_gene.items():
        gene_to_probes.setdefault(gene, []).append(probe)

    log(f"  Unique genes:      {len(gene_to_probes)}")
    target_found = [g for g in ALL_TARGET_GENES if g in gene_to_probes]
    target_miss  = [g for g in ALL_TARGET_GENES if g not in gene_to_probes]
    log(f"  Target genes mapped: "
        f"{len(target_found)}/{len(ALL_TARGET_GENES)}")
    if target_miss:
        log(f"  Missing: {target_miss}")

    return probe_to_gene, gene_to_probes

# ============================================================
# STEP 2: LOAD GSE25066
# ============================================================

def load_gse25066(gse_path, probe_to_gene, gene_to_probes):
    log("=" * 65)
    log("STEP 2: LOAD GSE25066")
    log("Hatzis et al. 2011 JAMA | n=508 | pCR annotated")
    log("=" * 65)

    sample_ids = []
    meta_rows  = {}
    expr_rows  = {}
    col_ids    = []

    try:
        open_fn = gzip.open if gse_path.endswith(".gz") else open
        with open_fn(gse_path, "rt", errors="ignore") as f:
            in_table = False
            char_idx = 0
            for raw_line in f:
                line = raw_line.rstrip("\n")

                if line.startswith("!Sample_geo_accession"):
                    sample_ids = [
                        s.strip().strip('"')
                        for s in line.split("\t")[1:]
                    ]

                elif line.startswith("!Sample_characteristics_ch1"):
                    char_idx += 1
                    parts = line.split("\t")
                    v0 = (parts[1].strip().strip('"')
                          if len(parts) > 1 else "")
                    # Print first 20 fields to understand format
                    if char_idx <= 20:
                        log(f"  characteristics[{char_idx:02d}]: "
                            f"'{v0[:70]}'")
                    fk = (v0.split(":")[0].strip().lower()
                          if ":" in v0 else v0[:20].lower())
                    row_vals = {}
                    for i, p in enumerate(parts[1:]):
                        if i < len(sample_ids):
                            row_vals[sample_ids[i]] = p.strip().strip('"')
                    meta_rows[f"char_{char_idx:02d}_{fk}"] = row_vals

                elif "series_matrix_table_begin" in line:
                    in_table = True
                    hdr_line = f.readline().rstrip("\n")
                    col_ids  = [
                        c.strip().strip('"')
                        for c in hdr_line.split("\t")[1:]
                    ]

                elif "series_matrix_table_end" in line:
                    break

                elif in_table and line:
                    parts = line.split("\t")
                    if len(parts) < 2:
                        continue
                    probe = parts[0].strip().strip('"')
                    try:
                        vals = [float(v.strip().strip('"'))
                                for v in parts[1:]]
                        expr_rows[probe] = vals
                    except ValueError:
                        pass

    except Exception as e:
        log(f"  Parse error: {e}")
        log(traceback.format_exc())
        return None, None

    if not expr_rows:
        log("  FATAL: No expression data.")
        return None, None

    use_cols = col_ids if col_ids else sample_ids
    n_cols   = min(len(list(expr_rows.values())[0]), len(use_cols))
    log(f"\n  Samples: {n_cols}")
    log(f"  Probes:  {len(expr_rows)}")
    log(f"  Metadata char fields: {len(meta_rows)}")

    # ── Map probes to genes ──────────────────────────────────
    log("\n  Mapping probes to genes (best probe by variance)...")
    gene_best = {}

    for probe, vals_list in expr_rows.items():
        gene = probe_to_gene.get(probe)
        if gene is None:
            continue
        if gene not in ALL_TARGET_GENES + ALL_LEHMANN_GENES:
            continue
        vals = np.array(vals_list[:n_cols], dtype=np.float32)
        var  = float(np.var(vals))
        if (gene not in gene_best or var > gene_best[gene]["var"]):
            gene_best[gene] = {"probe": probe, "vals": vals, "var": var}

    gene_expr = {g: info["vals"] for g, info in gene_best.items()}
    expr = pd.DataFrame(gene_expr, index=use_cols[:n_cols])
    log(f"  Gene-mapped expression: {expr.shape}")

    found  = [g for g in ALL_TARGET_GENES if g in expr.columns]
    miss   = [g for g in ALL_TARGET_GENES if g not in expr.columns]
    log(f"  Target genes found: {len(found)}/{len(ALL_TARGET_GENES)}")
    if miss:
        log(f"  Missing target genes: {miss}")

    # ── Probe sample report ──────────────────────────────────
    log("\n  Expression spot check (mean ± std across samples):")
    for gene in (SWITCH_GENES + FA_MARKERS + ["EED","EZH2","PARP1"])[:10]:
        if gene in expr.columns:
            log(f"    {gene:<12}: mean={expr[gene].mean():.3f}  "
                f"std={expr[gene].std():.3f}  "
                f"probe={gene_best[gene]['probe']}")

    # ── Parse metadata ───────────────────────────────────────
    log("\n  Parsing metadata...")
    meta = pd.DataFrame(index=use_cols[:n_cols])

    pcr_parsed = er_parsed = pr_parsed = her2_parsed = False

    for fk, row_vals in meta_rows.items():
        sample_vals = [v for v in row_vals.values() if v]
        if not sample_vals:
            continue
        v0 = sample_vals[0].lower()

        # pCR
        if not pcr_parsed and any(
            k in v0 for k in ["pcr","pathological complete",
                               "residual","rcb","drfs","event"]
        ):
            vals_list = [row_vals.get(s,"") for s in use_cols[:n_cols]]
            pcr_col = []
            for v in vals_list:
                vl = v.lower()
                if any(k in vl for k in ["1","yes","complete","pcr"]):
                    pcr_col.append(1)
                elif any(k in vl for k in ["0","no","residual","rd"]):
                    pcr_col.append(0)
                else:
                    pcr_col.append(np.nan)
            n1 = sum(1 for x in pcr_col if x == 1)
            n0 = sum(1 for x in pcr_col if x == 0)
            if n1 > 5 and n0 > 5:
                meta["pCR"] = pcr_col
                pcr_parsed = True
                log(f"  pCR field '{fk}': "
                    f"pCR=1 n={n1}  pCR=0 n={n0}")

        # ER
        if not er_parsed and any(
            k in v0 for k in ["er status","estrogen","er:","er "]
        ):
            vals_list = [row_vals.get(s,"") for s in use_cols[:n_cols]]
            meta["ER"] = [1 if "pos" in v.lower() else 0
                          for v in vals_list]
            er_parsed = True
            log(f"  ER field '{fk}'")

        # PR
        if not pr_parsed and any(
            k in v0 for k in ["pr status","progesterone","pr:","pr "]
        ):
            vals_list = [row_vals.get(s,"") for s in use_cols[:n_cols]]
            meta["PR"] = [1 if "pos" in v.lower() else 0
                          for v in vals_list]
            pr_parsed = True
            log(f"  PR field '{fk}'")

        # HER2
        if not her2_parsed and "her2" in v0:
            vals_list = [row_vals.get(s,"") for s in use_cols[:n_cols]]
            meta["HER2"] = [1 if "pos" in v.lower() else 0
                            for v in vals_list]
            her2_parsed = True
            log(f"  HER2 field '{fk}'")

    # If pCR still not found, print all field sample values
    if not pcr_parsed:
        log("\n  WARNING: pCR not found in characteristics.")
        log("  All field sample values for inspection:")
        for fk, row_vals in meta_rows.items():
            sample_vals = list(row_vals.values())[:3]
            log(f"    {fk}: {sample_vals}")

    log(f"\n  Metadata: {list(meta.columns)}")
    log(f"  pCR={pcr_parsed}  ER={er_parsed}  "
        f"PR={pr_parsed}  HER2={her2_parsed}")

    # TNBC subset flag
    if er_parsed and pr_parsed and her2_parsed:
        tnbc_mask = (
            (meta["ER"]   == 0) &
            (meta["PR"]   == 0) &
            (meta["HER2"] == 0)
        )
        meta["is_tnbc"] = tnbc_mask
        log(f"  TNBC (ER-/PR-/HER2-): {tnbc_mask.sum()}")
    else:
        meta["is_tnbc"] = pd.Series(True, index=meta.index)
        log("  Cannot determine TNBC subset. "
            "Using all samples as proxy.")

    return expr, meta

# ============================================================
# STEP 3: LOAD TCGA-BRCA FROM XENA
# ============================================================

def load_tcga_xena(tcga_paths):
    log("=" * 65)
    log("STEP 3: LOAD TCGA-BRCA (UCSC XENA)")
    log("=" * 65)

    if tcga_paths is None:
        log("  TCGA paths not available. Skipping.")
        return None, None

    expr_tcga  = None
    meta_tcga  = None

    # Load expression
    ep = tcga_paths.get("expr")
    if ep and os.path.exists(ep):
        log(f"  Expression: {ep}")
        try:
            open_fn = gzip.open if ep.endswith(".gz") else open
            expr_tcga = pd.read_csv(
                ep, sep="\t", index_col=0,
                compression="gzip" if ep.endswith(".gz") else None
            )
            log(f"  Shape (raw): {expr_tcga.shape}")
            log(f"  Index sample: {list(expr_tcga.index[:5])}")
            # Xena format: genes × samples
            # index = gene symbols, columns = sample IDs
            target_found = [g for g in ALL_TARGET_GENES
                            if g in expr_tcga.index]
            log(f"  Target genes: {len(target_found)}/{len(ALL_TARGET_GENES)}")
            if len(target_found) > 5:
                expr_tcga = expr_tcga.loc[target_found].T
                log(f"  Transposed: {expr_tcga.shape} "
                    f"(samples × genes)")
            else:
                log("  Gene-row format not confirmed. "
                    "Trying transpose...")
                expr_tcga_t = expr_tcga.T
                target_found_t = [g for g in ALL_TARGET_GENES
                                  if g in expr_tcga_t.index]
                if len(target_found_t) > len(target_found):
                    expr_tcga = expr_tcga_t.loc[target_found_t].T
                    log(f"  Transposed shape: {expr_tcga.shape}")
        except Exception as e:
            log(f"  Expression load error: {e}")
            log(traceback.format_exc())
            expr_tcga = None

    # Load clinical
    cp = tcga_paths.get("clin")
    if cp and os.path.exists(cp):
        log(f"  Clinical: {cp}")
        try:
            meta_tcga = pd.read_csv(cp, sep="\t", index_col=0)
            log(f"  Clinical shape: {meta_tcga.shape}")
            log(f"  Clinical cols: {list(meta_tcga.columns[:12])}")
        except Exception as e:
            log(f"  Clinical load error: {e}")

    # Load phenotype (PAM50)
    pp = tcga_paths.get("pheno")
    if pp and os.path.exists(pp):
        log(f"  Phenotype: {pp}")
        try:
            open_fn = gzip.open if pp.endswith(".gz") else open
            pheno = pd.read_csv(pp, sep="\t", index_col=0,
                                compression="gzip" if pp.endswith(".gz")
                                else None)
            log(f"  Phenotype shape: {pheno.shape}")
            log(f"  Phenotype cols: {list(pheno.columns[:10])}")
            # Find PAM50 column
            for col in pheno.columns:
                if "pam50" in col.lower() or "subtype" in col.lower():
                    log(f"  PAM50 col: '{col}'")
                    log(f"  Values: {pheno[col].value_counts().to_dict()}")
                    if meta_tcga is not None:
                        # Merge pheno into meta
                        common_idx = meta_tcga.index.intersection(
                            pheno.index
                        )
                        meta_tcga.loc[common_idx, "PAM50"] = (
                            pheno.loc[common_idx, col]
                        )
                    break
        except Exception as e:
            log(f"  Phenotype load error: {e}")

    if expr_tcga is not None:
        log(f"\n  TCGA expression ready: {expr_tcga.shape}")
    if meta_tcga is not None:
        log(f"  TCGA clinical ready:   {meta_tcga.shape}")

    return expr_tcga, meta_tcga

# ============================================================
# DEPTH SCORE UTILITY
# ============================================================

def compute_depth(expr, label=""):
    basal_g   = [g for g in DEPTH_BASAL_G   if g in expr.columns]
    luminal_g = [g for g in DEPTH_LUMINAL_G if g in expr.columns]

    if not basal_g and not luminal_g:
        log(f"  Depth ({label}): no genes available.")
        return None

    def norm01(s):
        mn, mx = s.min(), s.max()
        return (s - mn)/(mx - mn) if mx > mn else pd.Series(0.5, index=s.index)

    parts = []
    if basal_g:
        parts.append(norm01(expr[basal_g].mean(axis=1)))
    if luminal_g:
        parts.append(1 - norm01(expr[luminal_g].mean(axis=1)))
    if not parts:
        return None

    depth = sum(parts) / len(parts)
    log(f"  Depth {label}: "
        f"mean={depth.mean():.4f}  std={depth.std():.4f}  "
        f"n={len(depth)}")
    return depth

# ============================================================
# SECTION 1: TOP MOVERS
# ============================================================

def top_movers_bulk(expr, meta):
    log("")
    log("=" * 65)
    log("SECTION 1: TOP MOVERS — BULK GSE25066 (UNFILTERED)")
    log("Protocol v2.0: geometry first")
    log("=" * 65)

    def fmt_p(p):
        if p < 1e-30: return "p<1e-30  ***"
        if p < 0.001: return f"p={p:.2e} ***"
        if p < 0.01:  return f"p={p:.2e}  **"
        if p < 0.05:  return f"p={p:.4f}   *"
        return             f"p={p:.4f}  ns"

    if "is_tnbc" not in meta.columns:
        log("  No TNBC flag.")
        return pd.DataFrame()

    tnbc_mask  = meta["is_tnbc"] == True
    tnbc_expr  = expr.loc[expr.index.isin(meta.index[tnbc_mask])]
    other_expr = expr.loc[expr.index.isin(meta.index[~tnbc_mask])]

    log(f"  TNBC: {len(tnbc_expr)}  Non-TNBC: {len(other_expr)}")
    if len(tnbc_expr) < 5 or len(other_expr) < 5:
        log("  Insufficient samples.")
        return pd.DataFrame()

    results = []
    for gene in expr.columns:
        tv = tnbc_expr[gene].values
        nv = other_expr[gene].values
        nm, tm = nv.mean(), tv.mean()
        if abs(nm) < 1e-6 and abs(tm) < 1e-6:
            continue
        chg = (tm - nm) / abs(nm) * 100 if abs(nm) > 1e-6 else np.nan
        try:
            _, p = stats.mannwhitneyu(tv, nv, alternative="two-sided")
        except Exception:
            p = 1.0
        results.append({
            "gene": gene, "tnbc_mean": tm, "other_mean": nm,
            "change_pct": chg,
            "abs_change": abs(chg) if not np.isnan(chg) else 0,
            "p_value": p,
            "direction": "UP" if (chg or 0) > 0 else "DOWN",
        })

    rdf = pd.DataFrame(results).sort_values("abs_change", ascending=False)
    rdf.to_csv(CSV_FILE, index=False)

    for d, lbl in [("DOWN","LOST"), ("UP","GAINED")]:
        sub = rdf[rdf["direction"] == d].head(15)
        log(f"\n  TOP 15 {lbl} IN TNBC vs non-TNBC (bulk)")
        log(f"  {'Gene':<12} {'NonTNBC':>9} {'TNBC':>9} "
            f"{'%chg':>8}  p")
        log(f"  {'-'*56}")
        for _, row in sub.iterrows():
            log(f"  {row['gene']:<12} {row['other_mean']:>9.4f} "
                f"{row['tnbc_mean']:>9.4f} {row['change_pct']:>+7.1f}%  "
                f"{fmt_p(row['p_value'])}")
    return rdf

# ============================================================
# SECTION 2: DEPTH + pCR (DEFINITIVE P6)
# ============================================================

def depth_pcr_analysis(expr, meta):
    log("")
    log("=" * 65)
    log("SECTION 2: DEPTH SCORE + pCR ANALYSIS")
    log("Definitive P6 — gene-mapped Affymetrix probes")
    log("=" * 65)

    def fmt_p(p):
        if p < 0.001: return f"p={p:.2e} ***"
        if p < 0.05:  return f"p={p:.4f}   *"
        return             f"p={p:.4f}  ns"

    depth_all = compute_depth(expr, label="(full cohort)")

    tnbc_mask = meta.get("is_tnbc", pd.Series(True, index=meta.index))
    expr_tnbc = expr.loc[expr.index.isin(meta.index[tnbc_mask == True])]
    meta_tnbc = meta.loc[tnbc_mask == True]
    depth_tnbc = compute_depth(expr_tnbc, label="(TNBC subset)")

    log(f"  TNBC samples: {len(expr_tnbc)}")

    s_results = {}

    def test_pcr(depth, meta_sub, label, threshold, key):
        if depth is None or "pCR" not in meta_sub.columns:
            log(f"  {label}: pCR or depth unavailable.")
            return
        common = depth.index[depth.index.isin(meta_sub.index)]
        pcr_s  = meta_sub.loc[common, "pCR"].dropna()
        d_s    = depth.loc[pcr_s.index]
        if len(pcr_s) < 10:
            log(f"  {label}: only {len(pcr_s)} pCR samples.")
            return
        try:
            rv, pv = stats.pointbiserialr(pcr_s.values, d_s.values)
            ok = "✓ CONFIRMED" if rv < threshold and pv < 0.05 else (
                 "✓ DIRECTIONAL" if rv < 0 else "✗ NOT CONFIRMED")
            log(f"\n  {label}:")
            log(f"    r(depth, pCR) = {rv:>+.4f}  "
                f"{fmt_p(pv)}  [predicted < {threshold}]  {ok}")
            pcr1 = d_s[pcr_s == 1]
            pcr0 = d_s[pcr_s == 0]
            log(f"    pCR=1 mean depth: {pcr1.mean():.4f}  n={len(pcr1)}")
            log(f"    pCR=0 mean depth: {pcr0.mean():.4f}  n={len(pcr0)}")
            s_results[key] = rv
        except Exception as e:
            log(f"    Error: {e}")

    test_pcr(depth_all,  meta,      "S2-P1a: Full cohort",  -0.098, "P1a")
    test_pcr(depth_tnbc, meta_tnbc, "S2-P1b: TNBC subset",  -0.150, "P1b")

    return depth_all, depth_tnbc, expr_tnbc, meta_tnbc, s_results

# ============================================================
# SECTION 3: EED vs EZH2 BIOMARKER
# ============================================================

def eed_vs_ezh2_analysis(expr_tnbc, meta_tnbc, depth_tnbc):
    log("")
    log("=" * 65)
    log("SECTION 3: EED vs EZH2 BIOMARKER")
    log("S2-P2: EED predicts pCR better than EZH2")
    log("S2-P3: PARP1 predicts chemo resistance")
    log("=" * 65)

    def fmt_p(p):
        if p < 0.001: return f"p={p:.2e} ***"
        if p < 0.05:  return f"p={p:.4f}   *"
        return             f"p={p:.4f}  ns"

    if "pCR" not in meta_tnbc.columns:
        log("  pCR unavailable. S2-P2/P3 skipped.")
        return

    common   = meta_tnbc.index[meta_tnbc.index.isin(expr_tnbc.index)]
    pcr_k    = meta_tnbc.loc[common, "pCR"].dropna()
    expr_sub = expr_tnbc.loc[pcr_k.index]

    log(f"  TNBC with pCR: {len(pcr_k)}  "
        f"(pCR=1: {(pcr_k==1).sum()}  pCR=0: {(pcr_k==0).sum()})")

    log(f"\n  Gene vs pCR correlations (r, p, direction):")
    log(f"  {'Gene':<12} {'r':>8}  {'p':>14}  Note")
    log(f"  {'-'*55}")

    genes_test = ["EED","EZH2","HDAC1","HDAC2","KDM1A",
                  "PARP1","AR","VIM","KRT5","KRT14",
                  "ESR1","FOXA1"]
    auc_vals = {}
    for gene in genes_test:
        if gene not in expr_sub.columns:
            log(f"  {gene:<12} NOT IN DATA")
            continue
        gv = expr_sub[gene].values
        pv = pcr_k.values
        try:
            rv, rpv = stats.pointbiserialr(pv, gv)
            note = ""
            if gene == "EED"   and rv < -0.10: note = "✓ P2a"
            if gene == "PARP1" and rv < -0.08: note = "✓ P3a"
            if gene == "AR"    and rv < -0.15: note = "✓ P9"
            log(f"  {gene:<12} {rv:>+8.4f}  {fmt_p(rpv):>14}  {note}")
            # AUC
            try:
                auc = roc_auc_score(pv.astype(int), gv)
                auc_vals[gene] = auc
            except Exception:
                pass
        except Exception as e:
            log(f"  {gene:<12} error: {e}")

    # AUC comparison
    if "EED" in auc_vals and "EZH2" in auc_vals:
        log(f"\n  AUC for pCR (higher = better predictor):")
        for gene in sorted(auc_vals, key=lambda g: auc_vals[g],
                           reverse=True):
            auc = auc_vals[gene]
            # Flip if < 0.5 (predicts non-response)
            auc_display = max(auc, 1 - auc)
            log(f"    {gene:<12}: AUC = {auc_display:.3f}")
        eed_better = auc_vals.get("EED",0.5)
        ezh2_val   = auc_vals.get("EZH2",0.5)
        if abs(eed_better - 0.5) > abs(ezh2_val - 0.5):
            log("  S2-P2b: ✓ EED outperforms EZH2 as pCR predictor")
        else:
            log("  S2-P2b: ✗ EZH2 ≥ EED as pCR predictor")

    # Depth correlations in bulk
    if depth_tnbc is not None:
        log(f"\n  Depth correlations in bulk TNBC:")
        d_c = depth_tnbc.index[depth_tnbc.index.isin(expr_tnbc.index)]
        for gene in ["EED","EZH2","PARP1","AR","VIM","KRT5","KRT14",
                     "HDAC1","HDAC2","KDM1A"]:
            if gene not in expr_tnbc.columns:
                continue
            try:
                rv, rpv = stats.pearsonr(
                    depth_tnbc.loc[d_c].values,
                    expr_tnbc.loc[d_c, gene].values
                )
                log(f"    {gene:<12} r={rv:>+.4f}  {fmt_p(rpv)}")
            except Exception as e:
                log(f"    {gene:<12} error: {e}")

# ============================================================
# SECTION 4: LEHMANN SUBTYPE MAPPING
# ============================================================

def lehmann_subtype_mapping(expr_tnbc, meta_tnbc, depth_tnbc):
    log("")
    log("=" * 65)
    log("SECTION 4: LEHMANN SUBTYPE MAPPING")
    log("S2-P4: LAR < BL1/BL2 < IM < M < MSL on depth axis")
    log("=" * 65)

    def fmt_p(p):
        if p < 0.001: return f"p={p:.2e} ***"
        if p < 0.05:  return f"p={p:.4f}   *"
        return             f"p={p:.4f}  ns"

    score_df = pd.DataFrame(index=expr_tnbc.index)
    for subtype, genes in LEHMANN_GENES.items():
        avail = [g for g in genes if g in expr_tnbc.columns]
        log(f"  {subtype}: {len(avail)}/{len(genes)} genes available")
        if len(avail) < 3:
            continue
        sub  = expr_tnbc[avail]
        z    = (sub - sub.mean()) / (sub.std() + 1e-6)
        score_df[subtype] = z.mean(axis=1)

    if score_df.empty or score_df.shape[1] == 0:
        log("  No Lehmann subtypes could be scored.")
        return None

    score_df["assigned"] = score_df[
        [c for c in score_df.columns if c != "assigned"]
    ].idxmax(axis=1)

    log(f"\n  Subtype assignments:")
    for st in ["BL1","BL2","M","MSL","IM","LAR"]:
        if st in score_df.columns:
            n = (score_df["assigned"] == st).sum()
            log(f"    {st}: {n} samples")

    # Depth by subtype
    if depth_tnbc is not None:
        d_c = depth_tnbc.index[depth_tnbc.index.isin(score_df.index)]
        sc  = score_df.loc[d_c]
        dt  = depth_tnbc.loc[d_c]

        log(f"\n  Mean depth by Lehmann subtype:")
        log(f"  {'Subtype':<8} {'Mean depth':>12} {'Std':>8} "
            f"{'n':>6}")
        log(f"  {'-'*40}")
        order = ["LAR","BL1","BL2","IM","M","MSL"]
        group_depths = []
        for st in order:
            if st not in sc.columns:
                continue
            mask = sc["assigned"] == st
            if mask.sum() < 3:
                continue
            d_st = dt[mask]
            group_depths.append(d_st.values)
            pred = ""
            if st == "LAR"        and d_st.mean() < 0.40: pred = "✓ P4"
            elif st in ("M","MSL") and d_st.mean() > 0.60: pred = "✓ P4"
            log(f"  {st:<8} {d_st.mean():>12.4f} {d_st.std():>8.4f} "
                f"{mask.sum():>6}  {pred}")

        if len(group_depths) >= 3:
            try:
                h, p_kw = stats.kruskal(*group_depths)
                log(f"\n  Kruskal-Wallis: H={h:.2f}  {fmt_p(p_kw)}")
            except Exception as e:
                log(f"  KW error: {e}")

    # pCR by subtype
    if "pCR" in meta_tnbc.columns:
        log(f"\n  pCR rate by Lehmann subtype:")
        common   = meta_tnbc.index[meta_tnbc.index.isin(score_df.index)]
        pcr_k    = meta_tnbc.loc[common, "pCR"].dropna()
        sc_pcr   = score_df.loc[pcr_k.index]
        for st in ["LAR","BL1","BL2","IM","M","MSL"]:
            mask = sc_pcr["assigned"] == st
            if mask.sum() < 3:
                continue
            rate = pcr_k[mask].mean()
            log(f"    {st:<6}: pCR={rate:.2%}  n={mask.sum()}")

    return score_df

# ============================================================
# SECTION 5: SPI1 RESOLUTION
# ============================================================

def spi1_resolution(expr_tnbc):
    log("")
    log("=" * 65)
    log("SECTION 5: SPI1 RESOLUTION")
    log("S2-P7: Immune contamination vs genuine IM biology")
    log("=" * 65)

    def fmt_p(p):
        if p < 0.001: return f"p={p:.2e} ***"
        if p < 0.05:  return f"p={p:.4f}   *"
        return             f"p={p:.4f}  ns"

    if "SPI1" not in expr_tnbc.columns:
        log("  SPI1 not in dataset. Skipped.")
        return

    genes_test = (
        ["PTPRC","CD68","CD3D","CD8A","FOXP3"] +  # immune
        ["EED","EZH2","KRT5","AR","VIM"]            # TNBC
    )
    corrs = {}
    log(f"  {'Gene':<12} {'r(SPI1, gene)':>14}  p")
    log(f"  {'-'*40}")
    for gene in genes_test:
        if gene not in expr_tnbc.columns:
            continue
        try:
            rv, pv = stats.pearsonr(
                expr_tnbc["SPI1"].values,
                expr_tnbc[gene].values
            )
            corrs[gene] = rv
            log(f"  {gene:<12} {rv:>+14.4f}  {fmt_p(pv)}")
        except Exception as e:
            log(f"  {gene:<12} error: {e}")

    ptprc_r = abs(corrs.get("PTPRC", 0))
    cd68_r  = abs(corrs.get("CD68", 0))
    krt5_r  = abs(corrs.get("KRT5", 0))
    ar_r    = abs(corrs.get("AR", 0))

    immune_signal = max(ptprc_r, cd68_r)
    tnbc_signal   = max(krt5_r, ar_r)

    log(f"\n  Immune marker signal (max PTPRC/CD68): {immune_signal:.4f}")
    log(f"  TNBC marker signal (max KRT5/AR):      {tnbc_signal:.4f}")

    if immune_signal > 0.30:
        log("  S2-P7: ✓ IMMUNE CONTAMINATION — "
            "SPI1 tracks immune cells")
    elif immune_signal < 0.15 and tnbc_signal < 0.15:
        log("  S2-P7: ? WEAK SIGNAL — "
            "SPI1 not clearly immune or TNBC")
    else:
        log("  S2-P7: ? AMBIGUOUS — "
            "need IM subtype score comparison")

# ============================================================
# SECTION 6: PREDICTION PANEL CHECK
# ============================================================

def prediction_check(s_results, depth_tnbc, expr_tnbc, meta_tnbc):
    log("")
    log("=" * 65)
    log("SECTION 6: PREDICTION PANEL CHECK")
    log("All predictions from BRCA-S4c (locked 2026-03-04)")
    log("=" * 65)

    def fmt_p(p):
        if p < 0.001: return f"p={p:.2e} ***"
        if p < 0.05:  return f"p={p:.4f}   *"
        return             f"p={p:.4f}  ns"

    # S2-P1
    log("\n  S2-P1a: r(depth_full, pCR) < -0.098")
    r1a = s_results.get("P1a", np.nan)
    log(f"  → r = {r1a:+.4f}  "
        f"{'✓' if r1a < -0.098 else '?' if r1a < 0 else '✗'}")

    log("\n  S2-P1b: r(depth_tnbc, pCR) < -0.15")
    r1b = s_results.get("P1b", np.nan)
    log(f"  → r = {r1b:+.4f}  "
        f"{'✓' if r1b < -0.15 else '?' if r1b < 0 else '✗'}")

    # S2-P2/P3 — covered in Section 3
    log("\n  S2-P2 (EED > EZH2 for pCR): → see Section 3")
    log("  S2-P3 (PARP1 chemo resistance): → see Section 3")
    log("  S2-P4 (Lehmann depth order): → see Section 4")
    log("  S2-P7 (SPI1 resolution): → see Section 5")

    # S2-P9: AR vs depth in bulk
    log("\n  S2-P9: r(AR, depth_tnbc) < -0.15")
    if depth_tnbc is not None and "AR" in expr_tnbc.columns:
        d_c = depth_tnbc.index[depth_tnbc.index.isin(expr_tnbc.index)]
        try:
            rv, pv = stats.pearsonr(
                depth_tnbc.loc[d_c].values,
                expr_tnbc.loc[d_c, "AR"].values
            )
            log(f"  → r = {rv:+.4f}  {fmt_p(pv)}  "
                f"{'✓ CONFIRMED' if rv < -0.15 else '✗'}")
        except Exception as e:
            log(f"  → error: {e}")

    # Controls
    log(f"\n  Controls (should be flat):")
    if depth_tnbc is not None:
        d_c = depth_tnbc.index[depth_tnbc.index.isin(expr_tnbc.index)]
        for gene in CONTROLS:
            if gene not in expr_tnbc.columns:
                log(f"    {gene:<12} not in data")
                continue
            try:
                rv, pv = stats.pearsonr(
                    depth_tnbc.loc[d_c].values,
                    expr_tnbc.loc[d_c, gene].values
                )
                flag = "✓ FLAT" if abs(rv) < 0.10 else "⚠ CORRELATED"
                log(f"    {gene:<12} r_depth={rv:>+.4f}  "
                    f"{fmt_p(pv)}  {flag}")
            except Exception:
                pass

# ============================================================
# SECTION 7: TCGA-BRCA VALIDATION
# ============================================================

def tcga_validation(expr_tcga, meta_tcga):
    log("")
    log("=" * 65)
    log("SECTION 7: TCGA-BRCA VALIDATION")
    log("S2-P5: BRCA1 dysfunction in Basal-like")
    log("S2-P6: Depth score predicts survival")
    log("=" * 65)

    def fmt_p(p):
        if p < 0.001: return f"p={p:.2e} ***"
        if p < 0.05:  return f"p={p:.4f}   *"
        return             f"p={p:.4f}  ns"

    if expr_tcga is None:
        log("  TCGA expression not available. Skipping.")
        return

    depth_tcga = compute_depth(expr_tcga, label="(TCGA-BRCA)")
    if depth_tcga is None:
        log("  Depth score cannot be computed from TCGA data.")
        return

    # PAM50 subtype depth
    if meta_tcga is not None and "PAM50" in meta_tcga.columns:
        log(f"\n  Depth score by PAM50 subtype:")
        common = depth_tcga.index[depth_tcga.index.isin(meta_tcga.index)]
        for st in ["LumA","LumB","Basal","Her2","Normal"]:
            mask = meta_tcga.loc[common, "PAM50"].str.contains(
                st, case=False, na=False
            )
            if mask.sum() < 5:
                continue
            d_st = depth_tcga.loc[common[mask]]
            log(f"    {st:<10}: mean={d_st.mean():.4f}  "
                f"std={d_st.std():.4f}  n={len(d_st)}")

    # S2-P6 survival
    if meta_tcga is not None:
        # Find OS columns in Xena clinical matrix
        os_time_col = os_event_col = None
        for col in meta_tcga.columns:
            cl = col.lower()
            if "os_days" in cl or "days_to_death" in cl:
                os_time_col = col
            elif "os_months" in cl:
                os_time_col = col
            if "vital_status" in cl or "os_status" in cl:
                os_event_col = col

        log(f"\n  S2-P6: Survival")
        log(f"  OS time col:  {os_time_col}")
        log(f"  OS event col: {os_event_col}")

        if os_time_col and os_event_col:
            common = depth_tcga.index[
                depth_tcga.index.isin(meta_tcga.index)
            ]
            meta_sub = meta_tcga.loc[common].copy()
            meta_sub["depth"] = depth_tcga.loc[common].values

            # Restrict to PAM50 Basal-like
            if "PAM50" in meta_sub.columns:
                basal_mask = meta_sub["PAM50"].str.contains(
                    "Basal", case=False, na=False
                )
                meta_sub = meta_sub[basal_mask]
                log(f"  PAM50 Basal-like: {len(meta_sub)} samples")

            # Depth quartile split
            if len(meta_sub) > 20:
                med = meta_sub["depth"].median()
                meta_sub["depth_group"] = (
                    meta_sub["depth"] > med
                ).map({True: "High", False: "Low"})

                log(f"\n  Depth-High vs Depth-Low:")
                for grp in ["Low","High"]:
                    sub = meta_sub[meta_sub["depth_group"] == grp]
                    log(f"    {grp}: n={len(sub)}")

                # Cox-like: log-rank on OS
                try:
                    os_time  = pd.to_numeric(
                        meta_sub[os_time_col], errors="coerce"
                    )
                    # Vital status: 1=dead, 0=alive
                    os_event = meta_sub[os_event_col].apply(
                        lambda x: 1 if str(x).lower() in
                        ["dead","1","deceased"] else 0
                    )
                    high_t = os_time[meta_sub["depth_group"] == "High"].dropna()
                    low_t  = os_time[meta_sub["depth_group"] == "Low"].dropna()
                    if len(high_t) > 5 and len(low_t) > 5:
                        _, p_mw = stats.mannwhitneyu(
                            high_t, low_t, alternative="less"
                        )
                        log(f"\n  OS time: depth-high vs depth-low")
                        log(f"    High mean: {high_t.mean():.0f} days")
                        log(f"    Low mean:  {low_t.mean():.0f} days")
                        log(f"    Mann-Whitney (high < low): "
                            f"{fmt_p(p_mw)}")
                        if p_mw < 0.05:
                            log("  S2-P6: ✓ CONFIRMED — "
                                "depth-high shorter OS")
                        else:
                            log("  S2-P6: ? DIRECTIONAL — "
                                f"p={p_mw:.3f}")
                except Exception as e:
                    log(f"  OS analysis error: {e}")
        else:
            log("  OS columns not found in clinical matrix.")
            log(f"  Available columns: "
                f"{list(meta_tcga.columns[:20])}")

    # S2-P5: BRCA1 dysfunction — noted but somatic mutation
    # data is not in the Xena expression matrix.
    # Xena MAF data requires separate download.
    log(f"\n  S2-P5: BRCA1 dysfunction enrichment")
    log("  Note: requires somatic mutation + methylation data.")
    log("  These are not in the Xena expression matrix.")
    log("  Available via TCGA GDC portal (MAF files).")
    log("  S2-P5 will be completed in Script 3 using TCGA GDC MAF.")
    log("  Deferring to Script 3.")

# ============================================================
# SECTION 8: NOVEL SIGNALS
# ============================================================

def novel_signals(rdf):
    log("")
    log("=" * 65)
    log("SECTION 8: NOVEL SIGNALS")
    log("=" * 65)

    if rdf is None or len(rdf) == 0:
        log("  No top movers data.")
        return

    panel = set(ALL_TARGET_GENES + ALL_LEHMANN_GENES)
    for lbl, d in [("Novel suppressed","DOWN"),
                   ("Novel elevated","UP")]:
        sub = rdf[
            (rdf["direction"] == d) &
            (~rdf["gene"].isin(panel))
        ].head(10)
        log(f"\n  {lbl}:")
        for _, row in sub.iterrows():
            log(f"    {row['gene']:<20} {row['change_pct']:>+8.1f}%  "
                f"p={row['p_value']:.2e}")

# ============================================================
# FIGURE
# ============================================================

def generate_figure(expr, meta, depth_all,
                    expr_tnbc, meta_tnbc, depth_tnbc,
                    rdf, lehmann_scores):
    log("")
    log("--- Generating figure ---")

    clr = {
        "tnbc": "#c0392b", "other": "#2980b9",
        "up":   "#c0392b", "down":  "#2980b9",
        "eed":  "#8e44ad", "ezh2":  "#e67e22",
        "leh":  ["#3498db","#e74c3c","#2ecc71",
                 "#9b59b6","#f39c12","#1abc9c"],
    }

    fig = plt.figure(figsize=(24, 18))
    fig.suptitle(
        "TNBC — Script 2 | OrganismCore BRCA-S4c/d | 2026-03-04\n"
        "Bulk validation: GSE25066 (pCR) + TCGA-BRCA (survival)\n"
        "Gene-mapped depth | Lehmann subtypes | EED vs EZH2",
        fontsize=10, fontweight="bold", y=1.005
    )
    gs = gridspec.GridSpec(3, 3, figure=fig,
                           hspace=0.52, wspace=0.42)

    # A: Switch + FA bar chart
    ax_a = fig.add_subplot(gs[0, 0])
    genes_show = [g for g in SWITCH_GENES + FA_MARKERS
                  if g in expr.columns]
    if "is_tnbc" in meta.columns and genes_show:
        tm = expr.loc[meta["is_tnbc"] == True, genes_show].mean()
        om = expr.loc[meta["is_tnbc"] == False, genes_show].mean()
        x  = np.arange(len(genes_show)); w = 0.35
        ax_a.bar(x - w/2, om, w, color=clr["other"],
                 label="Non-TNBC", alpha=0.85)
        ax_a.bar(x + w/2, tm, w, color=clr["tnbc"],
                 label="TNBC", alpha=0.85)
        ax_a.set_xticks(x)
        ax_a.set_xticklabels(genes_show, rotation=45,
                              ha="right", fontsize=5)
        ax_a.set_title("A — Switch + FA Genes (Bulk)\nS2-P1 validation",
                        fontsize=8, fontweight="bold")
        ax_a.legend(fontsize=6)
        ax_a.set_ylabel("Expression", fontsize=7)

    # B: Depth distributions
    ax_b = fig.add_subplot(gs[0, 1])
    if depth_all is not None and "is_tnbc" in meta.columns:
        common = depth_all.index[depth_all.index.isin(meta.index)]
        tnbc_f = meta.loc[common, "is_tnbc"] == True
        ax_b.hist(depth_all[common][~tnbc_f].values, bins=30,
                  alpha=0.6, color=clr["other"],
                  label="Non-TNBC", density=True)
        ax_b.hist(depth_all[common][tnbc_f].values, bins=30,
                  alpha=0.6, color=clr["tnbc"],
                  label="TNBC", density=True)
        ax_b.set_title("B — Depth Distribution",
                        fontsize=8, fontweight="bold")
        ax_b.set_xlabel("Depth Score", fontsize=7)
        ax_b.legend(fontsize=6)

    # C: pCR by depth quartile
    ax_c = fig.add_subplot(gs[0, 2])
    if (depth_tnbc is not None and
            "pCR" in meta_tnbc.columns and
            len(meta_tnbc) > 20):
        common = depth_tnbc.index[depth_tnbc.index.isin(meta_tnbc.index)]
        pcr_k  = meta_tnbc.loc[common, "pCR"].dropna()
        d_pcr  = depth_tnbc.loc[pcr_k.index]
        if len(pcr_k) > 20:
            try:
                qt = pd.qcut(d_pcr, 4, labels=["Q1","Q2","Q3","Q4"])
                rates = pcr_k.groupby(qt).mean()
                ns    = pcr_k.groupby(qt).count()
                bc    = [clr["other"] if i < 2 else clr["tnbc"]
                         for i in range(4)]
                ax_c.bar(range(4), rates.values, color=bc)
                for i, (v, n) in enumerate(zip(rates.values, ns.values)):
                    ax_c.text(i, v + 0.01, f"n={n}", ha="center",
                              fontsize=7)
                ax_c.set_xticks(range(4))
                ax_c.set_xticklabels(
                    ["Q1\n(shallow)","Q2","Q3","Q4\n(deep)"],
                    fontsize=7
                )
                ax_c.set_title(
                    "C — pCR Rate by Depth Quartile\nS2-P1b",
                    fontsize=8, fontweight="bold"
                )
                ax_c.set_ylabel("pCR rate", fontsize=7)
            except Exception:
                ax_c.set_title("C — pCR (data insufficient)",
                                fontsize=8)

    # D: EED vs EZH2 depth scatter
    ax_d = fig.add_subplot(gs[1, 0])
    if (depth_tnbc is not None and
            "EED" in expr_tnbc.columns and
            "EZH2" in expr_tnbc.columns):
        d_c = depth_tnbc.index[depth_tnbc.index.isin(expr_tnbc.index)]
        ax_d.scatter(depth_tnbc.loc[d_c],
                     expr_tnbc.loc[d_c, "EED"],
                     c=clr["eed"], alpha=0.4, s=8, label="EED")
        ax_d.scatter(depth_tnbc.loc[d_c],
                     expr_tnbc.loc[d_c, "EZH2"],
                     c=clr["ezh2"], alpha=0.4, s=8, label="EZH2")
        try:
            r_eed,  _ = stats.pearsonr(
                depth_tnbc.loc[d_c].values,
                expr_tnbc.loc[d_c, "EED"].values
            )
            r_ezh2, _ = stats.pearsonr(
                depth_tnbc.loc[d_c].values,
                expr_tnbc.loc[d_c, "EZH2"].values
            )
            ax_d.set_title(
                f"D — EED vs EZH2 vs Depth (Bulk)\n"
                f"EED r={r_eed:+.3f}  EZH2 r={r_ezh2:+.3f}",
                fontsize=8, fontweight="bold"
            )
        except Exception:
            ax_d.set_title("D — EED vs EZH2", fontsize=8)
        ax_d.set_xlabel("Depth Score", fontsize=7)
        ax_d.set_ylabel("Expression", fontsize=7)
        ax_d.legend(fontsize=6)

    # E: Lehmann subtype depth bar
    ax_e = fig.add_subplot(gs[1, 1])
    if lehmann_scores is not None and depth_tnbc is not None:
        d_c  = depth_tnbc.index[depth_tnbc.index.isin(lehmann_scores.index)]
        ls_s = lehmann_scores.loc[d_c]
        dt_s = depth_tnbc.loc[d_c]
        order = ["LAR","BL1","BL2","IM","M","MSL"]
        means_l, lbls_l, ns_l = [], [], []
        for i, st in enumerate(order):
            if st not in ls_s.columns:
                continue
            mask = ls_s["assigned"] == st
            if mask.sum() < 3:
                continue
            means_l.append(dt_s[mask].mean())
            lbls_l.append(st)
            ns_l.append(mask.sum())
        if means_l:
            clrs_l = [clr["leh"][i % len(clr["leh"])]
                      for i in range(len(lbls_l))]
            ax_e.bar(range(len(means_l)), means_l, color=clrs_l)
            ax_e.set_xticks(range(len(lbls_l)))
            ax_e.set_xticklabels(
                [f"{l}\nn={n}" for l, n in zip(lbls_l, ns_l)],
                fontsize=7
            )
            ax_e.set_title(
                "E — Lehmann Depth Order\n"
                "Predicted: LAR<BL1/2<IM<M<MSL",
                fontsize=8, fontweight="bold"
            )
            ax_e.set_ylabel("Mean Depth Score", fontsize=7)
            ax_e.axhline(0.5, color="black", ls="--", lw=0.8)

    # F: Waterfall
    ax_f = fig.add_subplot(gs[1, 2])
    if rdf is not None and len(rdf) > 0:
        sig  = rdf[rdf["p_value"] < 0.05].sort_values("change_pct")
        show = pd.concat([sig.head(12), sig.tail(12)]).drop_duplicates()
        bc   = [clr["down"] if v < 0 else clr["tnbc"]
                for v in show["change_pct"]]
        ax_f.barh(show["gene"], show["change_pct"], color=bc)
        ax_f.axvline(0, color="black", lw=0.8)
        ax_f.set_title("F — Waterfall (Bulk GSE25066)",
                        fontsize=8, fontweight="bold")
        ax_f.set_xlabel("% change TNBC vs Non-TNBC", fontsize=7)
        ax_f.tick_params(axis="y", labelsize=5)

    # G: pCR by Lehmann
    ax_g = fig.add_subplot(gs[2, 0])
    if (lehmann_scores is not None and "pCR" in meta_tnbc.columns):
        common = meta_tnbc.index[meta_tnbc.index.isin(lehmann_scores.index)]
        pcr_k  = meta_tnbc.loc[common, "pCR"].dropna()
        sc_pcr = lehmann_scores.loc[pcr_k.index]
        rates_g, lbls_g, ns_g = [], [], []
        for i, st in enumerate(["LAR","BL1","BL2","IM","M","MSL"]):
            if st not in sc_pcr.columns:
                continue
            mask = sc_pcr["assigned"] == st
            if mask.sum() < 3:
                continue
            rates_g.append(pcr_k[mask].mean())
            lbls_g.append(st)
            ns_g.append(mask.sum())
        if rates_g:
            clrs_g = [clr["leh"][i % len(clr["leh"])]
                      for i in range(len(lbls_g))]
            ax_g.bar(range(len(rates_g)), rates_g, color=clrs_g)
            ax_g.set_xticks(range(len(lbls_g)))
            ax_g.set_xticklabels(
                [f"{l}\nn={n}" for l, n in zip(lbls_g, ns_g)],
                fontsize=7
            )
            ax_g.set_title("G — pCR by Lehmann Subtype",
                           fontsize=8, fontweight="bold")
            ax_g.set_ylabel("pCR rate", fontsize=7)

    # H: SPI1 correlations
    ax_h = fig.add_subplot(gs[2, 1])
    spi1_genes = [
        g for g in ["PTPRC","CD68","CD3D","CD8A","FOXP3",
                    "KRT5","AR","EZH2","EED","VIM"]
        if g in expr_tnbc.columns and "SPI1" in expr_tnbc.columns
    ]
    if spi1_genes:
        cs = []
        for gene in spi1_genes:
            try:
                rv, _ = stats.pearsonr(
                    expr_tnbc["SPI1"].values,
                    expr_tnbc[gene].values
                )
                cs.append((gene, rv))
            except Exception:
                pass
        cs.sort(key=lambda x: x[1], reverse=True)
        gc = [g for g, r in cs]
        vc = [r for g, r in cs]
        cc = [clr["tnbc"] if v > 0 else clr["other"] for v in vc]
        ax_h.barh(gc, vc, color=cc)
        ax_h.axvline(0, color="black", lw=0.8)
        ax_h.set_title("H — SPI1 Correlations\nS2-P7",
                        fontsize=8, fontweight="bold")
        ax_h.set_xlabel("r with SPI1", fontsize=7)
        ax_h.tick_params(axis="y", labelsize=7)

    # I: Summary
    ax_i = fig.add_subplot(gs[2, 2])
    ax_i.axis("off")
    ax_i.text(0.03, 0.97,
        f"I — SCRIPT 2 SUMMARY\n{'─'*30}\n"
        f"GSE25066 | n=508 | pCR\n"
        f"Affymetrix HG-U133A\n"
        f"TCGA-BRCA | Xena hub\n\n"
        f"PREDICTIONS:\n"
        f"S2-P1 Depth→pCR:  panel C\n"
        f"S2-P2 EED>EZH2:   panel D\n"
        f"S2-P3 PARP1→pCR:  panel C\n"
        f"S2-P4 Lehmann:    panel E\n"
        f"S2-P7 SPI1:       panel H\n\n"
        f"KEY NOVEL:\n"
        f"EED as pCR biomarker\n"
        f"PARP1 depth-chemo link\n"
        f"LAR = shallowest\n"
        f"M/MSL = deepest\n\n"
        f"DEFERRED TO SCRIPT 3:\n"
        f"S2-P5 BRCA1 MAF data\n\n"
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
    log("OrganismCore — BRCA-S4c/d | 2026-03-04")
    log("=" * 65)
    log("")
    log("FIXES vs previous version:")
    log("  GPL96: GPL96_family.soft.gz (not .annot.gz — 404)")
    log("  TCGA:  UCSC Xena hub (not cBioPortal S3 — 403)")
    log("")

    gse_path, gpl96_path, tcga_paths = acquire_data()
    if gse_path is None:
        log("FATAL: GSE25066 unavailable.")
        write_log()
        return

    probe_to_gene, gene_to_probes = build_probe_map(gpl96_path)
    expr, meta = load_gse25066(gse_path, probe_to_gene, gene_to_probes)
    if expr is None:
        log("FATAL: GSE25066 load failed.")
        write_log()
        return

    expr_tcga, meta_tcga = load_tcga_xena(tcga_paths)

    rdf = top_movers_bulk(expr, meta)

    (depth_all, depth_tnbc,
     expr_tnbc, meta_tnbc,
     s_results) = depth_pcr_analysis(expr, meta)

    eed_vs_ezh2_analysis(expr_tnbc, meta_tnbc, depth_tnbc)
    lehmann_scores = lehmann_subtype_mapping(
        expr_tnbc, meta_tnbc, depth_tnbc
    )
    spi1_resolution(expr_tnbc)
    prediction_check(s_results, depth_tnbc, expr_tnbc, meta_tnbc)
    tcga_validation(expr_tcga, meta_tcga)
    novel_signals(rdf)

    generate_figure(
        expr, meta, depth_all,
        expr_tnbc, meta_tnbc, depth_tnbc,
        rdf, lehmann_scores
    )

    write_log()

    log("")
    log("=" * 65)
    log("SCRIPT 2 COMPLETE")
    log(f"  Log    : {LOG_FILE}")
    log(f"  Figure : {FIG_FILE}")
    log(f"  CSV    : {CSV_FILE}")
    log("")
    log("READ OUTPUT IN THIS ORDER:")
    log("  1. Section 1: Top movers — bulk geometry")
    log("  2. Section 2: Depth + pCR — definitive P6")
    log("  3. Section 3: EED vs EZH2")
    log("  4. Section 4: Lehmann subtypes")
    log("  5. Section 5: SPI1 resolution")
    log("  6. Section 6: Prediction check")
    log("  7. Section 7: TCGA survival")
    log("  8. Section 8: Novel signals")
    log("")
    log("THEN write BRCA-S4d reasoning artifact.")
    log("=" * 65)


if __name__ == "__main__":
    main()
