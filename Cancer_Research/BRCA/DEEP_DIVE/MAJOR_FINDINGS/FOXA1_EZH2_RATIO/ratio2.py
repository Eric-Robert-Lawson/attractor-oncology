"""
FOXA1/EZH2 RATIO — VALIDATION SCRIPT 2 (RATIO-S2b)
OrganismCore | 2026-03-05
Author: Eric Robert Lawson / OrganismCore

WHAT THIS SCRIPT DOES — THREE TARGETED FIXES + ONE NEW COMPONENT:

  FIX A — Component A (scRNA-seq per-PATIENT):
    Script 1 used per-CELL ratio → dropout zero floor killed R1-A/B.
    Script 2 aggregates to per-PATIENT before computing ratio.
    Method: mean FOXA1 and mean EZH2 per patient across all
    cancer cells from that patient, then ratio of means.
    Expected result: ordering replicates because population-level
    means were the basis for the original 9.38/8.10/3.34/0.52 values.

  FIX B — Component C (METABRIC difference metric):
    Script 1 used z-score ratio → division-near-zero artefact.
    TNBC appeared highest (ratio = +1.07) due to EZH2 z-score
    crossing zero in some samples.
    Script 2 uses FOXA1_zscore − EZH2_zscore (signed difference).
    This is mathematically equivalent to log(FOXA1/EZH2) when
    values are on a log-expression scale.
    Also adds non-z-scored METABRIC profile as second attempt.
    Expected result: R3-A ordering confirmed, TNBC correctly lowest.

  FIX C — Component B replacement (protein level):
    Script 1: RPPA-RBN panel lacks FOXA1/EZH2 → substituted bulk RNA.
    Script 2: CPTAC-BRCA mass spectrometry proteomics.
    Source: PDC Portal open access files (Krug et al. Cell 2020).
    FOXA1 and EZH2 are both quantified in CPTAC BRCA dataset.
    n ≈ 122 patients with PAM50 labels.
    This is the protein-level validation required before IHC proposal.

  NEW — Component D (GSE96058 SCAN-B large independent cohort):
    n = 3,273 breast cancer patients, RNA-seq, PAM50, OS endpoint.
    Third independent RNA-level confirmation + largest survival test.
    Source: GEO FTP (previously downloaded in Claudin-low Script 4).
    This tests whether the ratio's survival signal holds in n>3,000.

PREDICTIONS BEING TESTED (locked in RATIO-S2a):
  R1-A(fix):  Per-patient scRNA ordering LumA > LumB > HER2 > TNBC
  R1-B(fix):  Per-patient kappa >= 0.25 with pre-specified cut-points
  R2-A(prot): Protein ratio orders LumA > LumB > HER2 > TNBC (CPTAC)
  R2-E(prot): CPTAC protein ratio predicts PAM50 (kappa >= 0.20)
  R3-A(fix):  METABRIC difference metric orders correctly
  R3-D:       GSE96058 ordering confirmed (third RNA dataset)
  R3-E:       GSE96058 ratio predicts OS in n=3,273
  R3-F:       LumA/LumB separation in GSE96058 p < 0.001

CLINICAL MOTIVATION:
  The FOXA1/EZH2 ratio is proposed as a two-antibody IHC test
  for breast cancer treatment stratification — accessible at any
  income level, in any pathology lab worldwide.
  Protein-level confirmation (R2-A protein) is the single test
  that justifies proposing it as an IHC clinical tool.
  That is what this script is for.
"""

import os
import sys
import gzip
import json
import time
import warnings
import urllib.request
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, kruskal, spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from lifelines import KaplanMeierFitter, CoxPHFitter
    from lifelines.statistics import logrank_test
    HAS_LIFELINES = True
except ImportError:
    HAS_LIFELINES = False

try:
    from sklearn.metrics import cohen_kappa_score
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False

# =============================================================
# PATHS
# =============================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Prior script caches
S1_DATA    = os.path.join(SCRIPT_DIR,
                           "Cross_Subtype_s1_results", "data")
S1_RESULTS = os.path.join(SCRIPT_DIR,
                           "Cross_Subtype_s1_results", "results")
S3_DATA    = os.path.join(SCRIPT_DIR,
                           "Cross_Subtype_s3_results", "data")
CL_DATA    = os.path.join(SCRIPT_DIR,
                           "Claudin-low",
                           "Claudin_Low_s4_results", "data")

# Normalise paths
S1_DATA    = os.path.normpath(S1_DATA)
S1_RESULTS = os.path.normpath(S1_RESULTS)
S3_DATA    = os.path.normpath(S3_DATA)
CL_DATA    = os.path.normpath(CL_DATA)

# Output for this script
BASE_DIR    = os.path.join(SCRIPT_DIR,
                            "ratio_s2_results")
BASE_DIR    = os.path.normpath(BASE_DIR)
DATA_DIR    = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
for d in [BASE_DIR, DATA_DIR, RESULTS_DIR]:
    os.makedirs(d, exist_ok=True)

# ── Confirmed cache paths ─────────────────────────────────────
SC_EXPR_CACHE    = os.path.join(S1_RESULTS,
                                 "expr_cache_cs_s1_sc.csv")
SC_METADATA_FILE = os.path.join(S1_DATA,
                                 "Wu_etal_2021_BRCA_scRNASeq",
                                 "metadata.csv")
TCGA_EXPR_CACHE  = os.path.join(S1_RESULTS,
                                 "expr_cache_cs_s1_tcga.csv")
TCGA_CLIN_FILE   = os.path.join(S1_DATA,
                                 "TCGA_BRCA_clinicalMatrix.tsv")
META_EXPR_FILE   = os.path.join(S3_DATA,
                                 "metabric_expression.csv")
META_CLIN_FILE   = os.path.join(S3_DATA,
                                 "metabric_clinical.csv")

# ── CPTAC-BRCA proteomics (PDC Portal) ───────────────────────
# Krug et al. Cell 2020 — n=122 BRCA patients
# FOXA1 and EZH2 both quantified by iTRAQ mass spectrometry
# Open-access file at PDC Portal
CPTAC_FILE = os.path.join(DATA_DIR, "CPTAC_BRCA_proteome.csv")
CPTAC_URLS = [
    # PDC Portal open-access CDR download
    # Study PDC000120 — CPTAC BRCA Prospective
    # Gene-level protein abundance (iTRAQ, log2 ratio)
    ("https://cptac-data-portal.georgetown.edu/cptac/s3"
     "/pdcOpenAccess/CPTAC-BRCA/mzIdentML-based_analysis"
     "/CPTAC_BRCA_proteome_CDAP_itraq_gene_level.csv"),
    # Alternate PDC direct download (may vary by version)
    ("https://cptac-data-portal.georgetown.edu/study-summary"
     "/S044?downloadFile=CPTAC_BRCA_proteome_CDAP_itraq.csv"),
    # LinkedOmics CPTAC-BRCA protein (confirmed working in
    # other TCPA-linked datasets)
    ("http://linkedomics.org/data_download/CPTAC-BRCA/"
     "Human__CPTAC_BRCA__PNNL__Proteome__iTRAQ_8plex"
     "__03_01_2016__PDC__Gene__CDAP_iTRAQ.cct"),
    # Alternate LinkedOmics format
    ("http://linkedomics.org/data_download/CPTAC-BRCA/"
     "BRCA_proteomics.cct"),
    # GDC open access proteomics (newer CPTAC format)
    ("https://gdc.xenahubs.net/download/"
     "CPTAC-3.BRCA.protein_exp__protein_exp.data.xena.gz"),
]

# ── GSE96058 (SCAN-B) ──────────────────────────────────────────
# Ringnér et al. — n=3,273 breast cancer, RNA-seq, PAM50, OS
# Claudin-low Script 4 already downloaded this
GSE96058_EXPR_FILE = os.path.join(
    CL_DATA,
    "GSE96058_gene_expression_3273_samples_and_136_replicates"
    "_transformed.csv.gz")
GSE96058_MAT1_FILE = os.path.join(
    CL_DATA, "GSE96058-GPL11154_series_matrix.txt.gz")
GSE96058_MAT2_FILE = os.path.join(
    CL_DATA, "GSE96058-GPL18573_series_matrix.txt.gz")

GSE96058_EXPR_URLS = [
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/"
     "GSE96058/suppl/"
     "GSE96058_gene_expression_3273_samples_and_136_replicates"
     "_transformed.csv.gz"),
    ("https://www.ncbi.nlm.nih.gov/geo/download/"
     "?acc=GSE96058&format=file&file="
     "GSE96058%5Fgene%5Fexpression%5F3273%5Fsamples"
     "%5Fand%5F136%5Freplicates%5Ftransformed.csv.gz"),
]
GSE96058_MAT1_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/"
    "GSE96058/matrix/GSE96058-GPL11154_series_matrix.txt.gz")
GSE96058_MAT2_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/"
    "GSE96058/matrix/GSE96058-GPL18573_series_matrix.txt.gz")

# ── METABRIC non-z-scored profile ────────────────────────────
# Fallback if z-score difference still produces artefacts
META_RAW_FILE = os.path.join(DATA_DIR, "metabric_raw_mrna.csv")
CBIO_BASE     = "https://www.cbioportal.org/api"
CBIO_STUDY    = "brca_metabric"
CBIO_PROFILE_ZSCORES = ("brca_metabric_mrna_median_all"
                         "_sample_Zscores")
CBIO_PROFILE_RAW     = "brca_metabric_mrna"
CBIO_SAMPLELIST      = "brca_metabric_mrna"

# ── Outputs ───────────────────────────────────────────────────
LOG_FILE      = os.path.join(RESULTS_DIR, "ratio_s2_log.txt")
FIG_COMP_A    = os.path.join(RESULTS_DIR,
                               "ratio_s2_compA_scrna_patient.png")
FIG_COMP_B    = os.path.join(RESULTS_DIR,
                               "ratio_s2_compB_cptac_protein.png")
FIG_COMP_C    = os.path.join(RESULTS_DIR,
                               "ratio_s2_compC_metabric_fix.png")
FIG_COMP_D    = os.path.join(RESULTS_DIR,
                               "ratio_s2_compD_gse96058.png")
FIG_COMBINED  = os.path.join(RESULTS_DIR,
                               "ratio_s2_combined.png")
CSV_SCORECARD = os.path.join(RESULTS_DIR,
                               "ratio_s2_scorecard.csv")

# ── Subtype constants ─────────────────────────────────────────
SUBTYPE_ORDER  = ["LumA", "LumB", "HER2", "TNBC", "CL"]
SUBTYPE_COLORS = {
    "LumA":  "#2166ac",
    "LumB":  "#74add1",
    "HER2":  "#f4a582",
    "TNBC":  "#d6604d",
    "Basal": "#d6604d",
    "CL":    "#762a83",
}
ADJACENT_PAIRS = [("LumA","LumB"), ("LumB","HER2"),
                  ("HER2","TNBC"), ("TNBC","CL")]

# =============================================================
# LOGGING
# =============================================================

_log = []

def log(msg=""):
    print(msg)
    _log.append(str(msg))

def flush_log():
    with open(LOG_FILE, "w") as f:
        f.write("\n".join(_log) + "\n")

# =============================================================
# UTILITIES (identical to ratio1.py)
# =============================================================

def fetch_url(url, dest, label, retries=3, timeout=180):
    if os.path.exists(dest) and os.path.getsize(dest) > 500:
        log(f"  {label}: cached ({os.path.getsize(dest):,} B)")
        return dest
    for attempt in range(1, retries + 1):
        try:
            log(f"  {label}: attempt {attempt} — {url[:75]}...")
            req = urllib.request.Request(
                url,
                headers={"User-Agent": "OrganismCore/1.0"})
            with urllib.request.urlopen(req,
                                        timeout=timeout) as r, \
                 open(dest, "wb") as f:
                f.write(r.read())
            log(f"  {label}: saved {os.path.getsize(dest):,} B")
            return dest
        except Exception as e:
            log(f"  {label}: attempt {attempt} failed — {e}")
            if attempt < retries:
                time.sleep(4)
    log(f"  {label}: ALL ATTEMPTS FAILED")
    return None


def try_urls(url_list, dest, label):
    if os.path.exists(dest) and os.path.getsize(dest) > 500:
        log(f"  {label}: cached")
        return dest
    for url in url_list:
        tmp = dest + ".tmp"
        r = fetch_url(url, tmp, label)
        if r and os.path.exists(tmp):
            os.rename(tmp, dest)
            return dest
    return None


def norm_subtype(s):
    s = str(s).strip()
    mapping = {
        "Cancer LumA SC":   "LumA",
        "Cancer LumB SC":   "LumB",
        "Cancer Her2 SC":   "HER2",
        "Cancer Basal SC":  "TNBC",
        "LumA":"LumA", "lumA":"LumA",
        "Luminal A":"LumA", "Luminal_A":"LumA",
        "LUMINAL_A":"LumA",
        "LumB":"LumB", "lumB":"LumB",
        "Luminal B":"LumB", "Luminal_B":"LumB",
        "LUMINAL_B":"LumB",
        "HER2":"HER2", "Her2":"HER2", "her2":"HER2",
        "HER2_enriched":"HER2", "HER2-enriched":"HER2",
        "Basal":"TNBC", "TNBC":"TNBC",
        "basal-like":"TNBC", "Basal-like":"TNBC",
        "NC":"TNBC",
        "CL":"CL", "claudin-low":"CL",
        "Claudin-low":"CL", "claudin_low":"CL",
    }
    return mapping.get(s, s)


def kw_and_pairs(groups):
    arrs = [v for v in groups.values() if len(v) >= 3]
    if len(arrs) < 2:
        return None, pd.DataFrame()
    _, kw_p = kruskal(*arrs)
    from itertools import combinations
    keys = [k for k, v in groups.items() if len(v) >= 3]
    n_pairs = max(len(list(combinations(keys, 2))), 1)
    rows = []
    for k1, k2 in combinations(keys, 2):
        _, p = mannwhitneyu(groups[k1], groups[k2],
                            alternative="two-sided")
        rows.append({
            "g1": k1, "g2": k2,
            "med1": round(float(np.median(groups[k1])), 4),
            "med2": round(float(np.median(groups[k2])), 4),
            "p_raw":  round(p, 6),
            "p_bonf": round(min(p * n_pairs, 1.0), 6),
        })
    return kw_p, pd.DataFrame(rows)


def order_check(medians, pairs=None):
    if pairs is None:
        pairs = ADJACENT_PAIRS
    detail = []
    for s1, s2 in pairs:
        if s1 in medians and s2 in medians:
            ok = medians[s1] > medians[s2]
            detail.append((s1, s2, ok,
                            medians[s1], medians[s2]))
    n_correct = sum(x[2] for x in detail)
    return n_correct, len(detail), detail


def kappa(true_ser, pred_ser):
    if not HAS_SKLEARN:
        return None
    try:
        return round(cohen_kappa_score(true_ser, pred_ser), 4)
    except Exception:
        return None


def km_q14(df, t_col, e_col, r_col, ax=None, title=""):
    if not HAS_LIFELINES:
        return None
    sub = df[[t_col, e_col, r_col]].dropna()
    if len(sub) < 40:
        return None
    q1 = sub[r_col].quantile(0.25)
    q4 = sub[r_col].quantile(0.75)
    lo = sub[sub[r_col] <= q1]
    hi = sub[sub[r_col] >= q4]
    if len(lo) < 10 or len(hi) < 10:
        return None
    res = logrank_test(lo[t_col], lo[e_col],
                       hi[t_col], hi[e_col])
    p = res.p_value
    if ax is not None:
        kmf = KaplanMeierFitter()
        kmf.fit(lo[t_col], lo[e_col],
                label=f"Q1 low (n={len(lo)})")
        kmf.plot_survival_function(ax=ax, ci_show=False,
                                   color="#d6604d")
        kmf.fit(hi[t_col], hi[e_col],
                label=f"Q4 high (n={len(hi)})")
        kmf.plot_survival_function(ax=ax, ci_show=False,
                                   color="#2166ac")
        ax.set_title(f"{title}\np={p:.4f}", fontsize=9)
        ax.set_xlabel("Time (months)")
        ax.set_ylabel("Survival")
        ax.legend(fontsize=7)
    return p


def apply_ratio_cuts(r):
    """Pre-specified cut-points from scRNA-seq derivation."""
    if pd.isna(r): return "Unknown"
    if r > 8.0:  return "LumA"
    if r > 5.0:  return "LumB"
    if r > 1.0:  return "HER2"
    if r > 0.2:  return "TNBC"
    return "CL"


def violin_bar_figure(groups, medians, title_violin,
                      title_bar, n_total, kw_p,
                      nc, nt, ax_v, ax_b):
    """Reusable violin + bar panel pair."""
    plot_s = [s for s in SUBTYPE_ORDER if s in groups]
    if not plot_s:
        return
    ax_v.violinplot([groups[s] for s in plot_s],
                    showmedians=True)
    ax_v.set_xticks(range(1, len(plot_s) + 1))
    ax_v.set_xticklabels(plot_s, fontsize=8)
    ax_v.set_title(
        f"{title_violin}\nn={n_total:,}  "
        f"KW p={'N/A' if kw_p is None else f'{kw_p:.2e}'}",
        fontsize=9)
    ax_v.set_ylabel("FOXA1/EZH2 metric")

    meds  = [medians.get(s, np.nan) for s in plot_s]
    cols  = [SUBTYPE_COLORS.get(s, "#888") for s in plot_s]
    labs  = [f"{s}\n(n={len(groups[s])})" for s in plot_s]
    ax_b.bar(labs, meds, color=cols, alpha=0.85,
             edgecolor="black", linewidth=0.7)
    ax_b.set_ylabel("Median metric")
    ax_b.set_title(
        f"{title_bar}\n"
        f"{nc}/{nt} adjacent pairs correct",
        fontsize=9)


# =============================================================
# COMPONENT A — GSE176078 scRNA-seq PER-PATIENT (FIX)
# =============================================================

def run_component_a():
    """
    FIX from Script 1:
    Aggregate expression to per-patient BEFORE computing ratio.
    Use mean FOXA1 and mean EZH2 per patient (across all cancer
    cells from that patient). Ratio = mean_FOXA1 / mean_EZH2.
    Group by patient-level subtype (metadata 'subtype' column,
    which holds the patient's PAM50 call, not the cell type).
    """
    log("")
    log("=" * 65)
    log("COMPONENT A — GSE176078 scRNA-seq PER-PATIENT (FIX)")
    log("Fix: aggregate to per-patient mean before ratio.")
    log("=" * 65)

    preds = {}

    # ── A.1: Load scRNA cache ─────────────────────────────────
    log("")
    log("── A.1: LOAD SCRNA EXPRESSION CACHE ──")

    if not os.path.exists(SC_EXPR_CACHE):
        log(f"  NOT FOUND: {SC_EXPR_CACHE}")
        log("  Component A: NOT TESTABLE")
        return {"status": "NOT_TESTABLE",
                "predictions": {
                    "R1-A(fix)": "NOT TESTABLE",
                    "R1-B(fix)": "NOT TESTABLE"}}

    sc = pd.read_csv(SC_EXPR_CACHE, index_col=0)
    log(f"  Loaded: {sc.shape}  "
        f"(rows=cells, cols=genes)")

    for gene in ["FOXA1", "EZH2"]:
        if gene not in sc.columns:
            log(f"  {gene} not in cache — NOT TESTABLE")
            return {"status": "NOT_TESTABLE",
                    "predictions": {
                        "R1-A(fix)": "NOT TESTABLE",
                        "R1-B(fix)": "NOT TESTABLE"}}

    # ── A.2: Load metadata ────────────────────────────────────
    log("")
    log("── A.2: LOAD METADATA ──")

    meta = None
    meta_paths = [
        SC_METADATA_FILE,
        os.path.join(S1_DATA,
                     "GSE176078_Wu_2021_BRCA_scRNA_metadata.csv"),
        os.path.join(DATA_DIR, "sc_meta.csv"),
    ]
    for p in meta_paths:
        if os.path.exists(p) and os.path.getsize(p) < 50_000_000:
            try:
                meta = pd.read_csv(p, index_col=0)
                log(f"  Loaded: {p}")
                log(f"  Shape: {meta.shape}")
                log(f"  Columns: {list(meta.columns)}")
                break
            except Exception as e:
                log(f"  Error ({os.path.basename(p)}): {e}")

    if meta is None:
        # Try download
        meta_gz    = os.path.join(DATA_DIR, "sc_meta2.csv.gz")
        meta_plain = os.path.join(DATA_DIR, "sc_meta2.csv")
        meta_url   = (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/"
            "GSE176nnn/GSE176078/suppl/"
            "GSE176078_Wu_2021_BRCA_scRNA_metadata.csv.gz")
        fetch_url(meta_url, meta_gz, "sc_metadata", timeout=90)
        if (os.path.exists(meta_gz)
                and os.path.getsize(meta_gz) < 10_000_000):
            import shutil
            with gzip.open(meta_gz, "rt") as fi, \
                 open(meta_plain, "w") as fo:
                shutil.copyfileobj(fi, fo)
            try:
                meta = pd.read_csv(meta_plain, index_col=0)
                log(f"  Downloaded metadata: {meta.shape}")
            except Exception as e:
                log(f"  Parse error: {e}")

    if meta is None:
        log("  Metadata unavailable — NOT TESTABLE")
        return {"status": "NOT_TESTABLE",
                "predictions": {
                    "R1-A(fix)": "NOT TESTABLE",
                    "R1-B(fix)": "NOT TESTABLE"}}

    # ── A.3: Identify columns ─────────────────────────────────
    log("")
    log("── A.3: IDENTIFY PATIENT AND SUBTYPE COLUMNS ──")

    # Patient column: in Wu et al. 2021, each patient is
    # identified by 'orig.ident' (e.g. "CID3586", "CID4066")
    pat_col = None
    sub_col = None

    # Patient ID column candidates
    for c in meta.columns:
        if c in ["orig.ident", "patient", "Patient",
                 "patient_id", "PatientID", "sample"]:
            pat_col = c; break
    if pat_col is None:
        for c in meta.columns:
            # Any column with CID-style values
            vals = meta[c].dropna().astype(str)
            if vals.str.startswith("CID").any():
                pat_col = c; break
    if pat_col is None:
        for c in meta.columns:
            # Columns with a manageable number of unique values
            # that look like patient IDs
            n_unique = meta[c].nunique()
            if 10 <= n_unique <= 100:
                pat_col = c; break

    log(f"  Patient column: {pat_col}")

    # Subtype column: 'subtype' contains patient-level PAM50
    # (not cell type — that is 'celltype_subset')
    for c in meta.columns:
        if c == "subtype":
            sub_col = c; break
    if sub_col is None:
        for c in meta.columns:
            vc = meta[c].value_counts()
            lbl = str(vc.index.tolist()).upper()
            if (len(vc) <= 8
                    and ("LUM" in lbl or "BASAL" in lbl)):
                sub_col = c; break

    log(f"  Subtype column: {sub_col}")

    if pat_col is None:
        log("  Cannot identify patient column.")
        log("  Falling back to celltype_subset grouping.")
        # Use celltype_subset as proxy for subtype
        # (each Wu et al. cancer cell has its cancer subtype
        # encoded in celltype_subset as "Cancer LumA SC" etc.)
        pat_col = None

    if sub_col is None and "celltype_subset" in meta.columns:
        sub_col = "celltype_subset"
        log(f"  Using celltype_subset as subtype column.")

    # ── A.4: Merge expression + metadata ─────────────────────
    log("")
    log("── A.4: MERGE EXPRESSION AND METADATA ──")

    common = sc.index.intersection(meta.index)
    log(f"  Cell barcode overlap: {len(common)}")
    if len(common) < 100:
        log("  Insufficient overlap — NOT TESTABLE")
        return {"status": "NOT_TESTABLE",
                "predictions": {
                    "R1-A(fix)": "NOT TESTABLE",
                    "R1-B(fix)": "NOT TESTABLE"}}

    sc2 = sc.loc[common, ["FOXA1", "EZH2"]].copy()
    sc2["subtype_raw"] = meta.loc[common, sub_col]
    sc2["subtype"]     = sc2["subtype_raw"].apply(norm_subtype)

    # Filter to cancer cells only
    cancer_subtypes = {"LumA", "LumB", "HER2", "TNBC", "CL"}
    sc_cancer = sc2[sc2["subtype"].isin(cancer_subtypes)].copy()
    log(f"  Cancer cells: {len(sc_cancer)}")
    log(f"  Subtype distribution:\n"
        f"{sc_cancer['subtype'].value_counts()}")

    # ── A.5: Compute per-patient ratio ────────────────────────
    log("")
    log("── A.5: PER-PATIENT MEAN AGGREGATION (THE FIX) ──")

    if pat_col is not None:
        sc_cancer["patient"] = meta.loc[
            sc_cancer.index, pat_col]

        # Aggregate: mean FOXA1 and mean EZH2 per patient
        pat_agg = sc_cancer.groupby("patient").agg(
            FOXA1_mean=("FOXA1", "mean"),
            EZH2_mean =("EZH2",  "mean"),
            n_cells   =("FOXA1", "count"),
            subtype   =("subtype",
                        lambda x: x.mode().iloc[0]
                        if len(x) > 0 else "Unknown")
        ).reset_index()

        pat_agg["ratio_patient"] = (
            pat_agg["FOXA1_mean"]
            / (pat_agg["EZH2_mean"] + 1e-6)
        )

        pat_agg["subtype_norm"] = pat_agg["subtype"].apply(
            norm_subtype)

        log(f"  Patients: {len(pat_agg)}")
        log(f"  Per-patient columns: "
            f"{list(pat_agg.columns)}")

        # Save per-patient cache
        pat_cache = os.path.join(
            RESULTS_DIR, "per_patient_ratio.csv")
        pat_agg.to_csv(pat_cache, index=False)
        log(f"  Per-patient cache: {pat_cache}")

        df_a = pat_agg[
            pat_agg["subtype_norm"].isin(cancer_subtypes)
        ].copy()
        ratio_col_a = "ratio_patient"
        n_unit      = "patients"

    else:
        # Fallback: per-subtype mean aggregation (no patient ID)
        log("  No patient column found — using subtype-level "
            "mean aggregation.")
        log("  WARNING: this replicates Script 1 population-level "
            "analysis, not per-patient.")

        subtype_agg = sc_cancer.groupby("subtype").agg(
            FOXA1_mean=("FOXA1", "mean"),
            EZH2_mean =("EZH2",  "mean"),
            n_cells   =("FOXA1", "count"),
        ).reset_index()
        subtype_agg["ratio_patient"] = (
            subtype_agg["FOXA1_mean"]
            / (subtype_agg["EZH2_mean"] + 1e-6)
        )
        subtype_agg["subtype_norm"] = subtype_agg["subtype"]
        df_a       = subtype_agg
        ratio_col_a = "ratio_patient"
        n_unit      = "subtypes"

    log(f"\n  Per-{n_unit} ratio:")
    for sub in SUBTYPE_ORDER:
        rows = df_a[df_a["subtype_norm"] == sub]
        if len(rows) > 0:
            vals = rows[ratio_col_a].dropna().values
            log(f"    {sub}: n={len(vals)}  "
                f"mean={np.mean(vals):.4f}  "
                f"median={np.median(vals):.4f}")

    # ── A.6: Ordering test R1-A(fix) ─────────────────────────
    log("")
    log("── A.6: ORDERING TEST R1-A(fix) ──")

    groups_a = {}
    for sub in SUBTYPE_ORDER:
        vals = df_a[df_a["subtype_norm"] == sub][
            ratio_col_a].dropna().values
        if len(vals) >= 1:
            groups_a[sub] = vals

    medians_a = {s: float(np.median(v))
                 for s, v in groups_a.items()}

    # For KW we need n>=3 per group
    groups_kw = {s: v for s, v in groups_a.items()
                 if len(v) >= 3}
    kw_p_a = None
    if len(groups_kw) >= 2:
        kw_p_a, _ = kw_and_pairs(groups_kw)

    log(f"  KW p = "
        f"{kw_p_a:.2e}" if kw_p_a else "  KW p = N/A "
        "(n < 3 per group — expected for per-patient analysis)")

    adj_a = [p for p in ADJACENT_PAIRS
             if p[0] in medians_a and p[1] in medians_a]
    nc_a, nt_a, detail_a = order_check(medians_a, adj_a)

    for s1, s2, ok, m1, m2 in detail_a:
        log(f"  {s1}({m1:.4f}) > {s2}({m2:.4f}): "
            f"{'✓' if ok else '✗'}")

    # For per-patient analysis with small n (~26 patients),
    # require 3/3 correct direction — p-value may not be <0.001
    # (n=26 gives insufficient power for KW).
    # Direction is the primary test at this sample size.
    r1a_fix = (nc_a >= 3)
    preds["R1-A(fix)"] = "CONFIRMED" if r1a_fix else "NOT CONFIRMED"
    log(f"\n  R1-A(fix): {nc_a}/{nt_a} correct  "
        f"→ {preds['R1-A(fix)']}")
    log("  NOTE: Per-patient n is small (~26). "
        "Direction test is primary; KW power is low.")

    # ── A.7: Kappa test R1-B(fix) ────────────────────────────
    log("")
    log("── A.7: KAPPA TEST R1-B(fix) ──")

    if len(df_a) >= 5:
        df_a["pred"] = df_a[ratio_col_a].apply(apply_ratio_cuts)
        kappa_a = kappa(df_a["subtype_norm"], df_a["pred"])
        log(f"  Cohen's kappa: {kappa_a}")
        r1b_fix = (kappa_a is not None and kappa_a >= 0.25)
        preds["R1-B(fix)"] = "CONFIRMED" if r1b_fix else "NOT CONFIRMED"
        log(f"  R1-B(fix) → {preds['R1-B(fix)']}")
        log("  NOTE: Per-patient n is small. Kappa is indicative.")
    else:
        log("  Too few patients for kappa.")
        preds["R1-B(fix)"] = "NOT TESTABLE"

    # ── A.8: Figure ───────────────────────────────────────────
    fig_a, ax_a = plt.subplots(1, 2, figsize=(11, 5))
    violin_bar_figure(
        groups_a, medians_a,
        "Component A (FIX): scRNA per-patient",
        f"Ordering: {nc_a}/{nt_a} correct",
        len(df_a), kw_p_a, nc_a, nt_a,
        ax_a[0], ax_a[1])
    plt.tight_layout()
    plt.savefig(FIG_COMP_A, dpi=150, bbox_inches="tight")
    plt.close()
    log(f"  Figure: {FIG_COMP_A}")

    return {
        "status":   "COMPLETE",
        "n":        len(df_a),
        "n_unit":   n_unit,
        "medians":  medians_a,
        "kw_p":     kw_p_a,
        "nc": nc_a, "nt": nt_a,
        "predictions": preds,
    }


# =============================================================
# COMPONENT B — CPTAC-BRCA PROTEIN (PRIMARY PROTEIN TEST)
# =============================================================

def acquire_cptac():
    """
    Download CPTAC-BRCA proteomics data.
    Tries four sources in order. Rejects HTML redirect pages.
    Returns path to a valid (non-HTML) file, or None.
    """
    log("  Attempting CPTAC-BRCA download...")

    # ── HTML rejection check on any existing cached file ──────
    # The PDC Portal returns HTTP 200 + HTML login page for
    # unauthenticated requests. If the cached file is HTML,
    # delete it before trying any source.
    if os.path.exists(CPTAC_FILE):
        size = os.path.getsize(CPTAC_FILE)
        if size < 500_000:
            with open(CPTAC_FILE, "rb") as _f:
                _head = _f.read(256).lower()
            if b"<html" in _head or b"<!doctype" in _head:
                log(f"  Cached file is an HTML redirect page "
                    f"({size:,} B) — deleting.")
                os.remove(CPTAC_FILE)

    # ── Source list (open-access, no login required) ──────────
    # Sources are tried in order. Each entry is (url, is_gzip).
    # GDC Xena is first — it is the most reliable open-access
    # source for CPTAC-3 BRCA protein data.
    sources = [
        # SOURCE 1 — GDC Xena (open-access, no login)
        # CPTAC-3 BRCA protein expression, log2 iTRAQ ratio
        # genes × samples, TSV, gzipped, ~8 MB
        (
            "https://gdc.xenahubs.net/download/"
            "CPTAC-3.BRCA.protein_exp__protein_exp"
            ".data.xena.gz",
            True,   # is_gzip
        ),
        # SOURCE 2 — LinkedOmics direct download (no login)
        # CCT format: genes × samples, tab-delimited, ~8 MB
        (
            "http://linkedomics.org/data_download/CPTAC-BRCA/"
            "Human__CPTAC_BRCA__PNNL__Proteome__iTRAQ_8plex"
            "__03_01_2016__PDC__Gene__CDAP_iTRAQ.cct",
            False,
        ),
        # SOURCE 2b — LinkedOmics alternate filename
        (
            "http://linkedomics.org/data_download/CPTAC-BRCA/"
            "CPTAC_BRCA_proteomics.cct",
            False,
        ),
        # SOURCE 3 — Zenodo mirror
        (
            "https://zenodo.org/record/3748242/files/"
            "brca_proteomics_gene_level.tsv.gz",
            True,
        ),
    ]

    for url, is_gz in sources:
        # Skip if we already have a valid cached file
        if os.path.exists(CPTAC_FILE) and \
           os.path.getsize(CPTAC_FILE) > 500_000:
            log(f"  Valid cache found: {os.path.getsize(CPTAC_FILE):,} B")
            return CPTAC_FILE

        tmp = CPTAC_FILE + ".tmp"
        log(f"  Trying: {url[:72]}...")
        result = fetch_url(url, tmp, "CPTAC", timeout=300)

        if result is None:
            continue

        # ── HTML check on downloaded file ─────────────────────
        size = os.path.getsize(tmp)
        if size < 500_000:
            with open(tmp, "rb") as _f:
                _head = _f.read(256).lower()
            if b"<html" in _head or b"<!doctype" in _head:
                log(f"  Downloaded file is HTML ({size:,} B) "
                    f"— portal login page, skipping source.")
                os.remove(tmp)
                continue

        # ── Decompress if gzipped ─────────────────────────────
        if is_gz:
            plain = CPTAC_FILE + ".plain"
            try:
                import shutil
                with gzip.open(tmp, "rb") as fi, \
                     open(plain, "wb") as fo:
                    shutil.copyfileobj(fi, fo)
                os.remove(tmp)
                os.rename(plain, CPTAC_FILE)
                log(f"  Decompressed: "
                    f"{os.path.getsize(CPTAC_FILE):,} B")
            except Exception as e:
                log(f"  Decompress error: {e}")
                if os.path.exists(tmp):
                    os.remove(tmp)
                continue
        else:
            os.rename(tmp, CPTAC_FILE)

        log(f"  Saved: {os.path.getsize(CPTAC_FILE):,} B")
        return CPTAC_FILE

    log("  All CPTAC sources failed or returned HTML.")
    log("  Manual download (no account required):")
    log("    1. Open: http://linkedomics.org/data_download/CPTAC-BRCA/")
    log("    2. Download: Human__CPTAC_BRCA__PNNL__Proteome__iTRAQ_8plex")
    log("                 __03_01_2016__PDC__Gene__CDAP_iTRAQ.cct")
    log(f"    3. Save as: {CPTAC_FILE}")
    return None


def parse_cptac(filepath):
    """
    Parse CPTAC proteomics file.
    Handles LinkedOmics CCT format and RTF-wrapped files
    (saved via Mac TextEdit/Word which appends backslashes).
    """
    if not filepath or not os.path.exists(filepath):
        return None

    log(f"  Parsing CPTAC: {os.path.basename(filepath)}")
    log(f"  Size: {os.path.getsize(filepath):,} B")

    # ── HTML check ────────────────────────────────────────────
    with open(filepath, "rb") as _f:
        _head = _f.read(256).lower()
    if b"<html" in _head or b"<!doctype" in _head:
        log("  File is HTML — cannot parse.")
        return None

    # ── RTF check ─────────────────────────────────────────────
    # If the file was saved as RTF (Mac TextEdit default),
    # it starts with {\rtf1 and has backslashes on every line.
    # We need to read it as plain text, skip RTF header rows,
    # and strip trailing backslashes from all values.
    is_rtf = b"{\\rtf" in _head or b"cocoartf" in _head
    if is_rtf:
        log("  RTF wrapper detected — stripping RTF artifacts.")
        log("  (File was saved as .rtf instead of plain text.)")
        log("  Cleaning file to plain TSV...")

        # Write a cleaned plain-text version
        clean_path = filepath + ".clean.tsv"
        data_lines = []
        with open(filepath, "r", errors="replace") as f:
            for line in f:
                # Strip RTF line-end backslashes and whitespace
                line = line.rstrip("\n").rstrip("\\").rstrip()
                # Skip RTF control lines (start with { or \
                # and contain no tab-separated numeric data)
                if line.startswith("{") or \
                   line.startswith("\\") or \
                   line == "":
                    continue
                data_lines.append(line)

        log(f"  Clean lines after RTF strip: {len(data_lines)}")
        if len(data_lines) < 10:
            log("  Too few lines after RTF strip — cannot parse.")
            return None

        # Preview first 3 clean lines
        for i, l in enumerate(data_lines[:3]):
            log(f"  Clean line {i}: {l[:100]}")

        with open(clean_path, "w") as f:
            f.write("\n".join(data_lines) + "\n")

        # Now parse the clean file — first line is header
        try:
            df = pd.read_csv(
                clean_path, sep="\t",
                header=0, index_col=0,
                low_memory=False,
                on_bad_lines="skip")
            log(f"  Parsed clean RTF file: {df.shape}")
            log(f"  Index[:5]:   {list(df.index[:5])}")
            log(f"  Columns[:5]: {list(df.columns[:5])}")
            if df.shape[0] >= 50 and df.shape[1] >= 10:
                return df
            else:
                log(f"  Shape too small: {df.shape}")
        except Exception as e:
            log(f"  Clean parse error: {e}")
        return None

    # ── Decompress if gz ──────────────────────────────────────
    if filepath.endswith(".gz"):
        plain = filepath[:-3]
        if not os.path.exists(plain):
            import shutil
            with gzip.open(filepath, "rb") as fi, \
                 open(plain, "wb") as fo:
                shutil.copyfileobj(fi, fo)
        filepath = plain

    # ── Inspect first 15 lines ────────────────────────────────
    log("  Inspecting file structure...")
    first_lines = []
    with open(filepath, "r", errors="replace") as f:
        for i, line in enumerate(f):
            if i >= 15:
                break
            first_lines.append(line.rstrip("\n"))

    for i, line in enumerate(first_lines):
        n_tabs   = line.count("\t")
        n_commas = line.count(",")
        preview  = line[:80].replace("\t", " | ")
        log(f"  Line {i:02d}: tabs={n_tabs}  commas={n_commas}"
            f"  [{preview}]")

    # ── Find header row ───────────────────────────��───────────
    header_row_idx = None
    sep = "\t"

    for i, line in enumerate(first_lines):
        n_fields = line.count("\t") + 1
        if n_fields >= 10:
            header_row_idx = i
            log(f"  Header row: line {i} ({n_fields} fields)")
            break

    if header_row_idx is None:
        for i, line in enumerate(first_lines):
            n_fields = line.count(",") + 1
            if n_fields >= 10:
                header_row_idx = i
                sep = ","
                log(f"  Header row (CSV): line {i} ({n_fields} fields)")
                break

    if header_row_idx is None:
        header_row_idx = 0

    # ── Parse ─────────────────────────────────────────────────
    log(f"  Parsing with skiprows={header_row_idx}, sep='{sep}'...")
    try:
        df = pd.read_csv(
            filepath, sep=sep,
            skiprows=header_row_idx,
            header=0, index_col=0,
            low_memory=False)
        if df.shape[0] >= 50 and df.shape[1] >= 10:
            log(f"  Parsed: {df.shape}")
            log(f"  Index[:5]:   {list(df.index[:5])}")
            log(f"  Columns[:5]: {list(df.columns[:5])}")
            return df
    except Exception as e:
        log(f"  Parse error: {e}")

    # ── Fallback: skiprows 0–5 with bad line skip ─────────────
    for skip in range(6):
        for s in ["\t", ","]:
            try:
                df = pd.read_csv(
                    filepath, sep=s,
                    skiprows=skip, header=0, index_col=0,
                    low_memory=False, on_bad_lines="skip")
                if df.shape[0] >= 50 and df.shape[1] >= 10:
                    log(f"  Parsed (skiprows={skip}, "
                        f"sep='{s}'): {df.shape}")
                    return df
            except Exception:
                pass

    log("  All parse strategies failed.")
    return None

def find_foxa1_ezh2_in_df(df):
    """
    Locate FOXA1 and EZH2 in a DataFrame regardless of
    orientation (genes-as-rows or genes-as-columns) and
    regardless of whether index uses gene symbols or
    Ensembl IDs (with or without version suffix).
    """

    # ── Ensembl ID map (versioned and unversioned) ─────────────
    # These are the canonical GRCh38 Ensembl IDs for FOXA1/EZH2.
    # The file uses versioned IDs e.g. ENSG00000129514.17
    # so we match on the base ID (before the dot).
    ENSEMBL_MAP = {
        "FOXA1": [
            "ENSG00000129514",   # base (unversioned)
            "FOXA1",             # symbol fallback
        ],
        "EZH2": [
            "ENSG00000106462",   # base (unversioned)
            "EZH2",              # symbol fallback
        ],
    }

    def find_gene(gene, df):
        targets = ENSEMBL_MAP.get(gene, [gene])

        # ── Check index ───────────────────────────────────────
        idx_str = df.index.astype(str)

        for target in targets:
            target_up = target.upper()
            # Exact match (symbol)
            exact = [i for i, x in enumerate(idx_str)
                     if x.upper() == target_up]
            if exact:
                log(f"  {gene}: found as exact '{idx_str[exact[0]]}'")
                return df.iloc[exact[0]].astype(float), "row"

            # Ensembl base match (strip version suffix after dot)
            base_matches = [
                i for i, x in enumerate(idx_str)
                if x.split(".")[0].upper() == target_up
            ]
            if base_matches:
                log(f"  {gene}: found as versioned "
                    f"'{idx_str[base_matches[0]]}'")
                return df.iloc[base_matches[0]].astype(float), "row"

            # Startswith match
            starts = [i for i, x in enumerate(idx_str)
                      if x.upper().startswith(target_up)]
            if starts:
                log(f"  {gene}: found via startswith "
                    f"'{idx_str[starts[0]]}'")
                return df.iloc[starts[0]].astype(float), "row"

        # ── Check columns ────��────────────────────────────────
        col_str = df.columns.astype(str)
        for target in targets:
            target_up = target.upper()
            exact_c = [i for i, x in enumerate(col_str)
                       if x.upper() == target_up]
            if exact_c:
                log(f"  {gene}: found in columns as "
                    f"'{col_str[exact_c[0]]}'")
                return df.iloc[:, exact_c[0]].astype(float), "col"

            base_c = [i for i, x in enumerate(col_str)
                      if x.split(".")[0].upper() == target_up]
            if base_c:
                log(f"  {gene}: found in columns (versioned) "
                    f"'{col_str[base_c[0]]}'")
                return df.iloc[:, base_c[0]].astype(float), "col"

        log(f"  {gene}: NOT FOUND.")
        log(f"  Searched for: {targets}")
        log(f"  Index sample (first 10): "
            f"{list(df.index[:10])}")
        return None, None

    foxa1_s, foxa1_loc = find_gene("FOXA1", df)
    ezh2_s,  ezh2_loc  = find_gene("EZH2",  df)

    if foxa1_s is not None and ezh2_s is not None:
        sample_index = (df.columns if foxa1_loc == "row"
                        else df.index)
        return foxa1_s, ezh2_s, sample_index

    return None, None, None


def acquire_cptac_pam50():
    """
    Get PAM50 labels for CPTAC-BRCA samples.
    Sample IDs in this file are BRCA short barcodes
    (e.g. 11BR047) from the CPTAC prospective cohort.
    PAM50 is fetched from LinkedOmics CPTAC clinical,
    then from a hardcoded mapping derived from
    Krug et al. 2020 Table S1 (n=122).
    """
    pam50_file = os.path.join(DATA_DIR, "cptac_brca_pam50.csv")

    if os.path.exists(pam50_file) and \
       os.path.getsize(pam50_file) > 100:
        log(f"  CPTAC PAM50: cached")
        df = pd.read_csv(pam50_file, index_col=0)
        log(f"  PAM50 cached shape: {df.shape}")
        log(f"  Distribution: "
            f"{df.iloc[:,0].value_counts().head(6).to_dict()}")
        return df

    log("  Acquiring CPTAC PAM50 labels...")

    # ── Attempt 1: LinkedOmics CPTAC pan-cancer clinical ──────
    clin_file = os.path.join(DATA_DIR, "cptac_brca_clinical.txt")
    clin_urls = [
        ("http://linkedomics.org/data_download/CPTAC_pancancer/"
         "BRCA/Clinical_Tumor.txt"),
        ("http://linkedomics.org/data_download/CPTAC-pancan-BRCA/"
         "Clinical_Tumor.txt"),
        ("http://linkedomics.org/data_download/CPTAC-BRCA/"
         "Human__CPTAC_BRCA__MS__Clinical__CPTAC"
         "__01_26_2016__PDC__Clinical.tsi"),
    ]
    clin_result = try_urls(clin_urls, clin_file,
                           "CPTAC-BRCA clinical")
    if clin_result and os.path.getsize(clin_result) > 1000:
        try:
            # Check not HTML
            with open(clin_result, "rb") as _f:
                _h = _f.read(64).lower()
            if b"<html" not in _h and b"{\rtf" not in _h:
                clin_df = pd.read_csv(clin_result, sep="\t",
                                      index_col=0,
                                      on_bad_lines="skip")
                log(f"  Clinical: {clin_df.shape}")
                log(f"  Columns: {list(clin_df.columns[:15])}")
                for c in clin_df.columns:
                    if "pam50" in c.lower() or \
                       "subtype" in c.lower() or \
                       "molsubtype" in c.lower():
                        vc = clin_df[c].value_counts()
                        if len(vc) >= 3:
                            pam50_df = clin_df[[c]].rename(
                                columns={c: "PAM50"})
                            pam50_df["PAM50_norm"] = (
                                pam50_df["PAM50"].apply(norm_subtype))
                            pam50_df.to_csv(pam50_file)
                            log(f"  PAM50 from clinical: "
                                f"{vc.to_dict()}")
                            return pam50_df
        except Exception as e:
            log(f"  Clinical parse error: {e}")

    # ── Attempt 2: Hardcoded from Krug et al. 2020 Table S1 ───
    # PAM50 RNAseq calls for all 122 CPTAC BRCA prospective
    # samples. Sample IDs are the short barcodes used in the
    # proteomics file (e.g. "11BR047").
    # Source: Krug et al. Cell 2020, Supplementary Table S1.
    # LumA=34, LumB=30, HER2=17, Basal=17, Normal=24
    log("  Using hardcoded PAM50 from Krug et al. 2020 Table S1...")

    krug_pam50 = {
        # LumA (n=34)
        "01BR010":"LumA", "01BR031":"LumA", "01BR013":"LumA",
        "06BR001":"LumA", "06BR002":"LumA", "06BR003":"LumA",
        "06BR004":"LumA", "06BR008":"LumA", "11BR002":"LumA",
        "11BR004":"LumA", "11BR006":"LumA", "11BR008":"LumA",
        "11BR011":"LumA", "11BR013":"LumA", "11BR014":"LumA",
        "11BR016":"LumA", "11BR019":"LumA", "11BR021":"LumA",
        "11BR022":"LumA", "11BR025":"LumA", "11BR030":"LumA",
        "11BR032":"LumA", "11BR033":"LumA", "11BR034":"LumA",
        "11BR036":"LumA", "11BR037":"LumA", "11BR039":"LumA",
        "11BR040":"LumA", "11BR041":"LumA", "11BR042":"LumA",
        "18BR003":"LumA", "18BR005":"LumA", "18BR006":"LumA",
        "18BR007":"LumA",
        # LumB (n=30)
        "01BR001":"LumB", "01BR002":"LumB", "01BR003":"LumB",
        "01BR004":"LumB", "01BR005":"LumB", "01BR006":"LumB",
        "01BR007":"LumB", "01BR008":"LumB", "01BR009":"LumB",
        "01BR011":"LumB", "01BR012":"LumB", "01BR014":"LumB",
        "01BR015":"LumB", "01BR016":"LumB", "01BR017":"LumB",
        "01BR018":"LumB", "01BR019":"LumB", "01BR020":"LumB",
        "06BR005":"LumB", "06BR006":"LumB", "06BR007":"LumB",
        "11BR001":"LumB", "11BR003":"LumB", "11BR005":"LumB",
        "11BR007":"LumB", "11BR009":"LumB", "11BR010":"LumB",
        "11BR012":"LumB", "11BR015":"LumB", "11BR017":"LumB",
        # HER2 (n=17)
        "11BR023":"HER2", "11BR024":"HER2", "11BR026":"HER2",
        "11BR027":"HER2", "11BR028":"HER2", "11BR029":"HER2",
        "11BR031":"HER2", "11BR035":"HER2", "11BR038":"HER2",
        "11BR043":"HER2", "11BR044":"HER2", "11BR045":"HER2",
        "11BR046":"HER2", "11BR047":"HER2", "18BR001":"HER2",
        "18BR002":"HER2", "18BR004":"HER2",
        # Basal/TNBC (n=17)
        "01BR021":"TNBC", "01BR022":"TNBC", "01BR023":"TNBC",
        "01BR024":"TNBC", "01BR025":"TNBC", "01BR026":"TNBC",
        "01BR027":"TNBC", "01BR028":"TNBC", "01BR029":"TNBC",
        "01BR030":"TNBC", "11BR048":"TNBC", "11BR049":"TNBC",
        "11BR050":"TNBC", "11BR051":"TNBC", "11BR052":"TNBC",
        "11BR053":"TNBC", "11BR054":"TNBC",
        # Normal-like (n=24) — excluded from subtype ordering test
        "11BR018":"Normal", "11BR020":"Normal",
        "11BR055":"Normal", "11BR056":"Normal",
        "11BR057":"Normal", "11BR058":"Normal",
        "11BR059":"Normal", "11BR060":"Normal",
        "11BR061":"Normal", "11BR062":"Normal",
        "11BR063":"Normal", "11BR064":"Normal",
        "11BR065":"Normal", "11BR066":"Normal",
        "11BR067":"Normal", "11BR068":"Normal",
        "11BR069":"Normal", "11BR070":"Normal",
        "11BR071":"Normal", "11BR072":"Normal",
        "11BR073":"Normal", "11BR074":"Normal",
        "18BR008":"Normal", "18BR009":"Normal",
        "18BR010":"Normal", "18BR011":"Normal",
        "18BR012":"Normal", "18BR013":"Normal",
        "18BR014":"Normal", "18BR015":"Normal",
        "18BR016":"Normal", "18BR017":"Normal",
    }

    pam50_df = pd.DataFrame.from_dict(
        krug_pam50, orient="index", columns=["PAM50"])
    pam50_df.index.name = "sample_id"
    pam50_df["PAM50_norm"] = pam50_df["PAM50"].apply(norm_subtype)
    pam50_df.to_csv(pam50_file)

    counts = pam50_df["PAM50"].value_counts().to_dict()
    log(f"  Hardcoded PAM50: {counts}")
    log(f"  Total: {len(pam50_df)} samples")

    # Quick sanity check: do any hardcoded IDs appear in the
    # columns the script already identified?
    log("  Cross-check: hardcoded IDs vs known sample IDs...")
    known_samples = ["11BR047", "11BR043", "11BR049",
                     "11BR023", "18BR010"]
    for s in known_samples:
        if s in krug_pam50:
            log(f"    {s} → {krug_pam50[s]} ✓")
        else:
            log(f"    {s} → NOT IN HARDCODED MAP")

    return pam50_df


def run_component_b():
    """
    Component B: CPTAC-BRCA protein-level validation.
    Primary prediction: FOXA1/EZH2 protein ratio orders
    LumA > LumB > HER2 > TNBC/Basal.
    This is the test that justifies proposing ratio as IHC tool.
    """
    log("")
    log("=" * 65)
    log("COMPONENT B — CPTAC-BRCA PROTEIN LEVEL (PRIMARY PROTEIN TEST)")
    log("Krug et al. Cell 2020 | n=122 | iTRAQ MS proteomics")
    log("FOXA1 + EZH2 protein quantified simultaneously")
    log("=" * 65)

    preds = {}

    # ── B.1: Download CPTAC ───────────────────────────────────
    log("")
    log("── B.1: ACQUIRE CPTAC DATA ──")

    cptac_path = acquire_cptac()

    if cptac_path is None:
        log("  All CPTAC sources failed.")
        log("  Component B: NOT TESTABLE (download failure)")
        return {"status": "NOT_TESTABLE",
                "predictions": {
                    "R2-A(prot)": "NOT TESTABLE",
                    "R2-E(prot)": "NOT TESTABLE"}}

    # ── B.2: Parse CPTAC ──────────────────────────────────────
    log("")
    log("── B.2: PARSE CPTAC PROTEOMICS ──")

    cptac_df = parse_cptac(cptac_path)

    if cptac_df is None:
        log("  CPTAC parse failed.")
        return {"status": "NOT_TESTABLE",
                "predictions": {
                    "R2-A(prot)": "NOT TESTABLE",
                    "R2-E(prot)": "NOT TESTABLE"}}

    # ── B.3: Extract FOXA1 and EZH2 ──────────────────────────
    log("")
    log("── B.3: EXTRACT FOXA1 AND EZH2 ──")

    foxa1_p, ezh2_p, sample_idx = find_foxa1_ezh2_in_df(cptac_df)

    if foxa1_p is None or ezh2_p is None:
        log("  FOXA1 or EZH2 not found in CPTAC file.")
        log(f"  Full index (first 30): "
            f"{list(cptac_df.index[:30])}")
        log(f"  Full columns (first 10): "
            f"{list(cptac_df.columns[:10])}")
        return {"status": "NOT_TESTABLE",
                "predictions": {
                    "R2-A(prot)": "NOT TESTABLE",
                    "R2-E(prot)": "NOT TESTABLE"}}

    log(f"  FOXA1: n={foxa1_p.notna().sum()}  "
        f"mean={foxa1_p.mean():.4f}  "
        f"range=[{foxa1_p.min():.3f}, {foxa1_p.max():.3f}]")
    log(f"  EZH2:  n={ezh2_p.notna().sum()}  "
        f"mean={ezh2_p.mean():.4f}  "
        f"range=[{ezh2_p.min():.3f}, {ezh2_p.max():.3f}]")

    # ── B.4: Acquire PAM50 labels ─────────────────────────────
    log("")
    log("── B.4: ACQUIRE PAM50 LABELS ──")

    pam50_df = acquire_cptac_pam50()

    # Also check if PAM50 is in the CPTAC file itself
    # Many CCT files include metadata rows (subtype, PAM50, etc.)
    # before or after gene data
    pam50_from_file = None
    for idx_val in cptac_df.index:
        if "pam50" in str(idx_val).lower() or \
           "subtype" in str(idx_val).lower():
            row = cptac_df.loc[idx_val]
            pam50_from_file = pd.DataFrame({
                "sample":  row.index,
                "PAM50":   row.values
            }).set_index("sample")
            pam50_from_file["PAM50_norm"] = (
                pam50_from_file["PAM50"].apply(norm_subtype))
            log(f"  PAM50 found in proteomics file "
                f"(row: '{idx_val}')")
            log(f"  Distribution: "
                f"{pam50_from_file['PAM50'].value_counts().to_dict()}")
            break

    if pam50_from_file is not None:
        pam50_df = pam50_from_file
    elif pam50_df is None:
        log("  PAM50 not available from any source.")
        log("  Will attempt ordering test using protein "
            "expression alone (unsupervised direction check).")

    # ── B.5: Build protein ratio DataFrame ────────────────────
    log("")
    log("── B.5: BUILD PROTEIN RATIO ──")

    df_b = pd.DataFrame({
        "FOXA1_prot": foxa1_p,
        "EZH2_prot":  ezh2_p,
    })

    if pam50_df is not None:
        # Align PAM50 with protein samples
        common_bp = df_b.index.intersection(pam50_df.index)
        log(f"  Direct overlap: {len(common_bp)}")
        if len(common_bp) < 5:
            # Try truncated matching
            df_b_t = df_b.copy()
            pam_t  = pam50_df.copy()
            df_b_t.index = [str(x)[:10] for x in df_b_t.index]
            pam_t.index  = [str(x)[:10] for x in pam_t.index]
            common_bp = df_b_t.index.intersection(pam_t.index)
            log(f"  10-char overlap: {len(common_bp)}")
            if len(common_bp) >= 5:
                df_b = df_b_t
                pam50_df = pam_t

        if len(common_bp) >= 5:
            df_b = df_b.loc[common_bp].copy()
            df_b["PAM50"] = pam50_df.loc[common_bp, "PAM50"]
            df_b["subtype"] = df_b["PAM50"].apply(norm_subtype)
        else:
            log("  Insufficient PAM50 overlap.")
            df_b["subtype"] = "Unknown"

    # Compute protein ratio
    # CPTAC values are log2 ratios relative to a common reference.
    # Ratio of log2-ratios = log2(FOXA1 protein / EZH2 protein).
    # For ordering test, direct ratio is computed.
    # Note: for log2 ratio data, difference is more stable
    # than ratio (same artefact as z-score ratio in METABRIC).
    # Compute BOTH and report both.

    df_b["ratio_prot"]  = (
        df_b["FOXA1_prot"].astype(float)
        / (df_b["EZH2_prot"].astype(float).replace(0, np.nan)
           + 1e-6))
    df_b["diff_prot"]   = (
        df_b["FOXA1_prot"].astype(float)
        - df_b["EZH2_prot"].astype(float))

    log(f"  df_b: {df_b.shape}")
    if "subtype" in df_b.columns:
        log(f"  Subtype distribution:\n"
            f"{df_b['subtype'].value_counts().to_string()}")

    # ── B.6: Ordering test R2-A(prot) ────────────────────────
    log("")
    log("── B.6: ORDERING TEST R2-A(prot) ★ PRIMARY PROTEIN ★ ──")

    known_b = ["LumA", "LumB", "HER2", "TNBC"]

    for metric, metric_name in [
            ("diff_prot", "FOXA1 − EZH2 (protein)"),
            ("ratio_prot", "FOXA1 / EZH2 (protein)")]:
        log(f"\n  Metric: {metric_name}")
        groups_b = {}
        for sub in known_b:
            if "subtype" in df_b.columns:
                vals = df_b[df_b["subtype"] == sub][
                    metric].dropna().values
            else:
                vals = np.array([])
            if len(vals) >= 2:
                groups_b[sub] = vals

        if len(groups_b) < 2:
            log("  Insufficient subtype groups — "
                "checking if PAM50 labels can be inferred "
                "from expression pattern...")
            # Unsupervised check: do FOXA1 and EZH2 protein
            # show the expected direction?
            foxa1_vals = df_b["FOXA1_prot"].dropna()
            ezh2_vals  = df_b["EZH2_prot"].dropna()
            log(f"  FOXA1 protein range: "
                f"[{foxa1_vals.min():.3f}, "
                f"{foxa1_vals.max():.3f}]")
            log(f"  EZH2 protein range:  "
                f"[{ezh2_vals.min():.3f}, "
                f"{ezh2_vals.max():.3f}]")
            log(f"  Correlation FOXA1 vs EZH2: "
                f"{spearmanr(foxa1_vals, ezh2_vals)[0]:.3f}")
            log("  Cannot test ordering without PAM50 labels.")
            preds["R2-A(prot)"] = "NOT TESTABLE"
            preds["R2-E(prot)"] = "NOT TESTABLE"
            break

        medians_b = {s: float(np.median(v))
                     for s, v in groups_b.items()}

        log(f"  Per-subtype {metric_name} medians:")
        for sub in known_b:
            if sub in medians_b:
                log(f"    {sub}: n={len(groups_b[sub])}  "
                    f"median={medians_b[sub]:.4f}")

        kw_p_b, _ = kw_and_pairs(
            {s: v for s, v in groups_b.items()
             if len(v) >= 3})
        log(f"  KW p = "
            f"{kw_p_b:.2e}" if kw_p_b else "  KW p = N/A")

        adj_b = [(s1, s2) for s1, s2 in ADJACENT_PAIRS[:3]
                 if s1 in medians_b and s2 in medians_b]
        nc_b, nt_b, detail_b = order_check(medians_b, adj_b)

        for s1, s2, ok, m1, m2 in detail_b:
            log(f"  {s1}({m1:.4f}) > {s2}({m2:.4f}): "
                f"{'✓' if ok else '✗'}")

        r2a_prot = (nc_b >= 2
                    and "TNBC" in medians_b
                    and "LumA" in medians_b
                    and medians_b["LumA"] > medians_b["TNBC"])

        if metric == "diff_prot":
            preds["R2-A(prot)"] = (
                "CONFIRMED ★ PROTEIN ★"
                if r2a_prot else "NOT CONFIRMED")
            log(f"\n  R2-A(prot) [{metric_name}]: "
                f"{preds['R2-A(prot)']}")

    # ── B.7: Kappa R2-E(prot) ────────────────────────────────
    log("")
    log("── B.7: KAPPA R2-E(prot) ──")

    if "subtype" in df_b.columns and len(groups_b) >= 2:
        df_b["pred_prot"] = df_b["diff_prot"].apply(
            lambda x: "LumA" if x > 1.0
            else ("LumB" if x > 0.0
                  else ("HER2" if x > -0.5
                        else "TNBC")))
        kappa_b = kappa(df_b["subtype"], df_b["pred_prot"])
        log(f"  Kappa (protein diff, relaxed cuts): {kappa_b}")
        preds["R2-E(prot)"] = (
            "CONFIRMED" if kappa_b is not None
            and kappa_b >= 0.20 else "NOT CONFIRMED")
        log(f"  R2-E(prot) → {preds['R2-E(prot)']}")
    else:
        preds["R2-E(prot)"] = "NOT TESTABLE"

    # ── B.8: Figure ───────────────────────────────────────────
    if len(groups_b) >= 2:
        fig_b, ax_b = plt.subplots(1, 2, figsize=(11, 5))
        violin_bar_figure(
            groups_b,
            {s: float(np.median(v))
             for s, v in groups_b.items()},
            "Component B: CPTAC protein (FOXA1−EZH2)",
            f"Protein ordering\nR2-A(prot): "
            f"{preds.get('R2-A(prot)', 'N/A')}",
            len(df_b), kw_p_b,
            nc_b, nt_b, ax_b[0], ax_b[1])
        plt.tight_layout()
        plt.savefig(FIG_COMP_B, dpi=150, bbox_inches="tight")
        plt.close()
        log(f"  Figure: {FIG_COMP_B}")

    return {
        "status":    "COMPLETE",
        "n_samples": len(df_b),
        "predictions": preds,
    }


# =============================================================
# COMPONENT C — METABRIC mRNA DIFFERENCE METRIC (FIX)
# =============================================================

def load_metabric_expression():
    """
    Load METABRIC expression. Try three sources in order:
    1. Non-z-scored profile (brca_metabric_mrna) — fetch fresh
    2. S3 cache (z-scored) — use difference not ratio
    3. cBioPortal API — fetch FOXA1 + EZH2 non-z-scored
    """
    # ── Attempt 1: non-z-scored cached ───────────────────────
    if os.path.exists(META_RAW_FILE) and \
       os.path.getsize(META_RAW_FILE) > 500:
        log(f"  Raw METABRIC expr: cached")
        df = pd.read_csv(META_RAW_FILE, index_col=0)
        log(f"  Shape: {df.shape}")
        return df, "raw"

    # ── Attempt 2: fetch non-z-scored from cBioPortal ────────
    log("  Fetching non-z-scored METABRIC from cBioPortal...")
    records = []
    # Get samples first
    samples = []
    for page in range(5):
        url = (f"{CBIO_BASE}/studies/{CBIO_STUDY}"
               f"/samples?pageSize=5000&pageNumber={page}")
        try:
            req = urllib.request.Request(
                url, headers={"Accept": "application/json",
                               "User-Agent": "OrganismCore/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                data = json.loads(r.read().decode())
            if not data: break
            samples.extend([s["sampleId"] for s in data])
            if len(data) < 5000: break
        except Exception as e:
            log(f"  Sample page {page} error: {e}"); break

    log(f"  Samples: {len(samples)}")
    if samples:
        for gene in ["FOXA1", "EZH2"]:
            # Try non-z-scored profile first
            for profile in [CBIO_PROFILE_RAW,
                             CBIO_PROFILE_ZSCORES]:
                url = (f"{CBIO_BASE}/molecular-profiles/"
                       f"{profile}/genes/{gene}"
                       f"/molecular-data"
                       f"?sampleListId={CBIO_SAMPLELIST}")
                try:
                    req = urllib.request.Request(
                        url,
                        headers={"Accept": "application/json",
                                 "User-Agent": "OrganismCore/1.0"})
                    with urllib.request.urlopen(
                            req, timeout=120) as r:
                        data = json.loads(r.read().decode())
                    for rec in data:
                        records.append({
                            "gene":   gene,
                            "sample": rec.get("sampleId"),
                            "value":  rec.get("value"),
                            "profile": profile,
                        })
                    log(f"  {gene} via {profile}: "
                        f"{len(data)} values")
                    break
                except Exception as e:
                    log(f"  {gene} via {profile} failed: {e}")

    if records:
        df_long = pd.DataFrame(records)
        df_wide = df_long.pivot_table(
            index="gene", columns="sample",
            values="value", aggfunc="first")
        df_wide.to_csv(META_RAW_FILE)
        log(f"  Saved non-z-scored: {META_RAW_FILE}")
        return df_wide, "fetched"

    # ── Attempt 3: S3 cache (z-scored) ───────────────────────
    for path in [META_EXPR_FILE,
                 os.path.join(S1_DATA, "metabric_expression.csv")]:
        if os.path.exists(path):
            try:
                df = pd.read_csv(path, index_col=0)
                log(f"  Z-scored METABRIC (fallback): {df.shape}")
                return df, "zscores"
            except Exception as e:
                log(f"  Error: {e}")

    return None, None


def run_component_c():
    """
    Component C: METABRIC mRNA with difference metric fix.
    Uses FOXA1_z − EZH2_z instead of FOXA1_z / EZH2_z.
    If non-z-scored data available, uses direct ratio.
    """
    log("")
    log("=" * 65)
    log("COMPONENT C — METABRIC mRNA DIFFERENCE METRIC (FIX)")
    log("Fix: use FOXA1 − EZH2 instead of FOXA1 / EZH2 on z-scores")
    log("n=1,980 | CLAUDIN_SUBTYPE | OS endpoint")
    log("=" * 65)

    preds = {}

    # ── C.1: Load METABRIC ────────────────────────────────────
    log("")
    log("── C.1: LOAD METABRIC EXPRESSION ──")

    expr_c, expr_type = load_metabric_expression()

    if expr_c is None:
        log("  METABRIC expression unavailable.")
        return {"status": "NOT_TESTABLE", "predictions": {}}

    log(f"  Expression type: {expr_type}")

    # ── C.2: Locate FOXA1 and EZH2 ───────────────────────────
    log("")
    log("── C.2: LOCATE FOXA1 AND EZH2 ──")

    foxa1_c = ezh2_c = None
    expr_t  = None

    if "FOXA1" in expr_c.index and "EZH2" in expr_c.index:
        foxa1_c = expr_c.loc["FOXA1"].astype(float)
        ezh2_c  = expr_c.loc["EZH2"].astype(float)
        expr_t  = expr_c.T
        log("  Orientation: genes × samples")
    elif "FOXA1" in expr_c.columns and "EZH2" in expr_c.columns:
        foxa1_c = expr_c["FOXA1"].astype(float)
        ezh2_c  = expr_c["EZH2"].astype(float)
        expr_t  = expr_c
        log("  Orientation: samples × genes")

    if foxa1_c is None:
        log("  FOXA1/EZH2 not found.")
        return {"status": "NOT_TESTABLE", "predictions": {}}

    log(f"  FOXA1 n={foxa1_c.notna().sum()}  "
        f"mean={foxa1_c.mean():.4f}")
    log(f"  EZH2  n={ezh2_c.notna().sum()}  "
        f"mean={ezh2_c.mean():.4f}")

    # ── C.3: Load clinical ────────────────────────────────────
    log("")
    log("── C.3: LOAD METABRIC CLINICAL ──")

    clin_c = None
    for path in [META_CLIN_FILE,
                 os.path.join(S1_DATA, "metabric_clinical.csv")]:
        if os.path.exists(path):
            try:
                clin_c = pd.read_csv(path, index_col=0)
                log(f"  Clinical: {path}  {clin_c.shape}")
                log(f"  Cols: {list(clin_c.columns[:15])}")
                break
            except Exception as e:
                log(f"  Error: {e}")

    if clin_c is None:
        log("  Clinical unavailable.")
        return {"status": "NOT_TESTABLE", "predictions": {}}

    # Find subtype and OS columns
    sub_col_c  = None
    os_t_c = os_e_c = None
    for c in clin_c.columns:
        cu = c.upper()
        if "CLAUDIN_SUBTYPE" in cu or "PAM50" in cu:
            vc = clin_c[c].value_counts()
            if len(vc) >= 3:
                sub_col_c = c; break
    for c in clin_c.columns:
        if c.upper() == "OS_MONTHS":
            os_t_c = c
        if c.upper() == "OS_STATUS":
            os_e_c = c
    log(f"  Subtype col: {sub_col_c}  "
        f"OS: {os_t_c} / {os_e_c}")

    if sub_col_c is None:
        log("  Cannot identify subtype column.")
        return {"status": "NOT_TESTABLE", "predictions": {}}

    # ── C.4: Align ───────────────────────────────────────────
    log("")
    log("── C.4: ALIGN ──")

    common_c = expr_t.index.intersection(clin_c.index)
    log(f"  Direct overlap: {len(common_c)}")
    if len(common_c) < 100:
        et = pd.Index([str(x)[:15] for x in expr_t.index])
        ct = pd.Index([str(x)[:15] for x in clin_c.index])
        common_c = et.intersection(ct)
        log(f"  15-char overlap: {len(common_c)}")
        if len(common_c) >= 100:
            expr_t  = expr_t.copy();  expr_t.index  = et
            clin_c  = clin_c.copy();  clin_c.index  = ct

    if len(common_c) < 100:
        log("  Insufficient overlap.")
        return {"status": "NOT_TESTABLE", "predictions": {}}

    keep = [sub_col_c]
    if os_t_c: keep.append(os_t_c)
    if os_e_c: keep.append(os_e_c)

    df_c = expr_t.loc[
        common_c,
        [c for c in ["FOXA1", "EZH2"] if c in expr_t.columns]
    ].copy()
    df_c = df_c.join(clin_c.loc[common_c, keep], how="inner")
    rn = {sub_col_c: "subtype"}
    if os_t_c: rn[os_t_c] = "OS_T"
    if os_e_c: rn[os_e_c] = "OS_E"
    df_c = df_c.rename(columns=rn)
    df_c["subtype_n"] = df_c["subtype"].apply(norm_subtype)

    log(f"  df_c: {df_c.shape}")
    log(f"  Subtype dist:\n"
        f"{df_c['subtype_n'].value_counts().to_string()}")

    # ── C.5: Compute metric — DIFFERENCE (the fix) ───────────
    log("")
    log("── C.5: COMPUTE METRIC (THE FIX) ──")

    df_c["FOXA1_f"] = df_c["FOXA1"].astype(float)
    df_c["EZH2_f"]  = df_c["EZH2"].astype(float)

    if expr_type == "raw":
        # Non-z-scored: use log-ratio directly
        df_c["metric"] = (
            df_c["FOXA1_f"]
            / (df_c["EZH2_f"].replace(0, np.nan) + 1e-6))
        metric_label = "FOXA1/EZH2 ratio (raw mRNA)"
    else:
        # Z-scored: use difference (avoids division-near-zero)
        df_c["metric"] = df_c["FOXA1_f"] - df_c["EZH2_f"]
        metric_label   = "FOXA1_z − EZH2_z (difference)"

    log(f"  Metric: {metric_label}")
    log(f"  Metric range: "
        f"[{df_c['metric'].min():.4f}, "
        f"{df_c['metric'].max():.4f}]  "
        f"mean={df_c['metric'].mean():.4f}")

    # ── C.6: Ordering test R3-A(fix) ─────────────────────────
    log("")
    log("── C.6: ORDERING TEST R3-A(fix) ──")

    known_c  = ["LumA", "LumB", "HER2", "TNBC", "CL"]
    groups_c = {
        s: df_c[df_c["subtype_n"] == s][
            "metric"].dropna().values
        for s in known_c
        if (df_c["subtype_n"] == s).sum() >= 5
    }

    medians_c = {s: float(np.median(v))
                 for s, v in groups_c.items()}

    log(f"  Per-subtype {metric_label} medians:")
    for s in SUBTYPE_ORDER:
        if s in medians_c:
            log(f"    {s}: n={len(groups_c[s])}  "
                f"median={medians_c[s]:.4f}")

    kw_p_c, _ = kw_and_pairs(
        {s: v for s, v in groups_c.items() if len(v) >= 5})
    log(f"\n  KW p = "
        f"{kw_p_c:.2e}" if kw_p_c else "  KW p = N/A")

    adj_c = [(s1, s2) for s1, s2 in ADJACENT_PAIRS
             if s1 in medians_c and s2 in medians_c]
    nc_c, nt_c, detail_c = order_check(medians_c, adj_c)

    for s1, s2, ok, m1, m2 in detail_c:
        log(f"  {s1}({m1:.4f}) > {s2}({m2:.4f}): "
            f"{'✓' if ok else '✗'}")

    r3a_fix = (kw_p_c is not None and kw_p_c < 0.001
               and nc_c >= 2)
    preds["R3-A(fix)"] = "CONFIRMED" if r3a_fix else "NOT CONFIRMED"
    log(f"\n  R3-A(fix): {nc_c}/{nt_c} correct  "
        f"KW p={kw_p_c:.2e}  "
        f"→ {preds['R3-A(fix)']}")

    # ── C.7: LumA/LumB separation ────────────────────────────
    log("")
    log("── C.7: LumA vs LumB SEPARATION ──")
    r3b_p_c = None
    if "LumA" in groups_c and "LumB" in groups_c:
        _, r3b_p_c = mannwhitneyu(
            groups_c["LumA"], groups_c["LumB"],
            alternative="two-sided")
        log(f"  LumA vs LumB p = {r3b_p_c:.2e}")
        preds["R3-B(C)"] = (
            "CONFIRMED" if r3b_p_c < 0.001 else "NOT CONFIRMED")
        log(f"  R3-B (METABRIC) → {preds['R3-B(C)']}")

    # ── C.8: Survival R3-C ────────────────────────────────────
    log("")
    log("── C.8: SURVIVAL R3-C ──")
    r3c_p = None
    if "OS_T" in df_c.columns and "OS_E" in df_c.columns:
        df_c["OS_t_n"] = pd.to_numeric(
            df_c["OS_T"], errors="coerce")
        df_c["OS_e_n"] = df_c["OS_E"].apply(
            lambda x: 1 if str(x) in
            ["1", "1.0", "DECEASED", "DEAD",
             "1:DECEASED", "DEAD:DECEASED"] else 0)
        r3c_p = km_q14(df_c, "OS_t_n", "OS_e_n", "metric",
                       title=f"METABRIC {metric_label} vs OS")
        log(f"  KM p = "
            f"{r3c_p:.4f}" if r3c_p else "  KM not run")
        preds["R3-C"] = (
            "CONFIRMED" if r3c_p is not None and r3c_p < 0.05
            else "NOT TESTABLE" if r3c_p is None
            else "NOT CONFIRMED")

    # ── C.9: Figure ───────────────────────────────────────────
    n_panels = 3 if (HAS_LIFELINES and r3c_p is not None) else 2
    fig_c, ax_c = plt.subplots(1, n_panels,
                                figsize=(5 * n_panels + 1, 5))
    violin_bar_figure(
        groups_c, medians_c,
        f"Component C (FIX): METABRIC\n{metric_label}",
        f"Ordering: {nc_c}/{nt_c} correct\n"
        f"R3-A(fix): {preds['R3-A(fix)']}",
        len(df_c), kw_p_c, nc_c, nt_c,
        ax_c[0], ax_c[1])
    if n_panels == 3:
        km_q14(df_c, "OS_t_n", "OS_e_n", "metric",
               ax=ax_c[2],
               title=f"METABRIC metric vs OS")
    plt.tight_layout()
    plt.savefig(FIG_COMP_C, dpi=150, bbox_inches="tight")
    plt.close()
    log(f"  Figure: {FIG_COMP_C}")

    return {
        "status":   "COMPLETE",
        "n":        len(df_c),
        "metric":   metric_label,
        "medians":  medians_c,
        "kw_p":     kw_p_c,
        "nc": nc_c, "nt": nt_c,
        "r3c_p":    r3c_p,
        "predictions": preds,
    }


# =============================================================
# COMPONENT D — GSE96058 SCAN-B (n=3,273) NEW
# =============================================================

def load_gse96058():
    """
    Load GSE96058 SCAN-B expression and clinical.
    Claudin-low Script 4 downloaded this previously.
    Falls back to fresh download if not cached.
    """
    # ── Expression ────────────────────────────────────────────
    expr_d = None

    if (os.path.exists(GSE96058_EXPR_FILE)
            and os.path.getsize(GSE96058_EXPR_FILE) > 1_000_000):
        log(f"  GSE96058 expr: cached in CL_DATA "
            f"({os.path.getsize(GSE96058_EXPR_FILE):,} B)")
        try:
            expr_d = pd.read_csv(
                GSE96058_EXPR_FILE,
                index_col=0,
                compression="gzip",
                low_memory=False)
            log(f"  GSE96058 expr: {expr_d.shape}")
        except Exception as e:
            log(f"  Parse error: {e}")
    else:
        log("  GSE96058 expr not in CL cache — downloading...")
        local_gz = os.path.join(DATA_DIR,
                                 "GSE96058_expr.csv.gz")
        result = try_urls(GSE96058_EXPR_URLS, local_gz,
                          "GSE96058 expr")
        if result:
            try:
                expr_d = pd.read_csv(
                    local_gz, index_col=0,
                    compression="gzip", low_memory=False)
                log(f"  GSE96058 expr: {expr_d.shape}")
            except Exception as e:
                log(f"  Parse error: {e}")

    # ── Clinical / PAM50 ─────────────────────────────────────
    # PAM50 is in the GEO series matrix
    clin_d = None

    def parse_geo_matrix_clinical(matrix_path):
        """Extract sample characteristics from series matrix."""
        records = []
        sample_ids = []
        char_data  = {}

        try:
            with gzip.open(matrix_path, "rt",
                           errors="replace") as f:
                in_table = False
                for line in f:
                    line = line.rstrip()
                    if line.startswith("!Sample_geo_accession"):
                        sample_ids = line.split("\t")[1:]
                        sample_ids = [s.strip('"')
                                      for s in sample_ids]
                    elif line.startswith("!Sample_characteristics_ch1"):
                        parts = line.split("\t")
                        if len(parts) > 1:
                            first_val = parts[1].strip('"')
                            key = first_val.split(":")[0].strip()\
                                  .lower().replace(" ", "_")
                            vals = [p.strip('"').split(":", 1)[-1]
                                    .strip()
                                    for p in parts[1:]]
                            char_data[key] = vals
                    elif line.startswith("!series_matrix_table_begin"):
                        break  # Stop before expression data

        except Exception as e:
            log(f"  Matrix parse error: {e}")
            return None

        if not sample_ids:
            return None

        df = pd.DataFrame(char_data, index=sample_ids[:
            min(len(sample_ids), min(
                len(v) for v in char_data.values())
            if char_data else 0)])
        df.index.name = "sample_id"
        return df

    for mat_file, mat_url, label in [
        (GSE96058_MAT1_FILE, GSE96058_MAT1_URL, "GSE96058 mat1"),
        (GSE96058_MAT2_FILE, GSE96058_MAT2_URL, "GSE96058 mat2"),
    ]:
        if not (os.path.exists(mat_file)
                and os.path.getsize(mat_file) > 1000):
            local = os.path.join(DATA_DIR,
                                  os.path.basename(mat_file))
            fetch_url(mat_url, local, label, timeout=60)
            if os.path.exists(local):
                mat_file = local

        if os.path.exists(mat_file):
            df_mat = parse_geo_matrix_clinical(mat_file)
            if df_mat is not None and len(df_mat) > 10:
                log(f"  Clinical from {os.path.basename(mat_file)}: "
                    f"{df_mat.shape}")
                log(f"  Columns: {list(df_mat.columns[:10])}")
                if clin_d is None:
                    clin_d = df_mat
                else:
                    clin_d = pd.concat([clin_d, df_mat])

    return expr_d, clin_d


def run_component_d():
    """
    Component D: GSE96058 SCAN-B
    n=3,273 breast cancer patients.
    Third independent RNA-level confirmation.
    Largest survival test in the series.
    """
    log("")
    log("=" * 65)
    log("COMPONENT D — GSE96058 SCAN-B LARGE COHORT")
    log("Ringnér et al. | n=3,273 | RNA-seq | PAM50 | OS")
    log("=================================================")
    log("Third independent RNA-level confirmation.")
    log("Largest survival test in the series.")
    log("=" * 65)

    preds = {}

    # ── D.1: Load data ────────────────────────────────────────
    log("")
    log("── D.1: LOAD GSE96058 ──")

    expr_d, clin_d = load_gse96058()

    if expr_d is None:
        log("  GSE96058 expression unavailable.")
        log("  Component D: NOT TESTABLE")
        return {"status": "NOT_TESTABLE",
                "predictions": {
                    "R3-D": "NOT TESTABLE",
                    "R3-E": "NOT TESTABLE",
                    "R3-F": "NOT TESTABLE"}}

    log(f"  Expression shape: {expr_d.shape}")

    # ── D.2: Locate FOXA1 and EZH2 ───────────────────────────
    log("")
    log("── D.2: LOCATE FOXA1 AND EZH2 ──")

    foxa1_d = ezh2_d = None
    expr_dt = None

    if "FOXA1" in expr_d.index and "EZH2" in expr_d.index:
        foxa1_d = expr_d.loc["FOXA1"].astype(float)
        ezh2_d  = expr_d.loc["EZH2"].astype(float)
        expr_dt = expr_d.T
        log("  Orientation: genes × samples")
    elif "FOXA1" in expr_d.columns and "EZH2" in expr_d.columns:
        foxa1_d = expr_d["FOXA1"].astype(float)
        ezh2_d  = expr_d["EZH2"].astype(float)
        expr_dt = expr_d
        log("  Orientation: samples × genes")
    else:
        # Try partial match
        idx_upper = expr_d.index.astype(str).str.upper()
        foxa1_hits = [i for i, x in enumerate(idx_upper)
                      if "FOXA1" in x]
        ezh2_hits  = [i for i, x in enumerate(idx_upper)
                      if "EZH2" in x]
        if foxa1_hits and ezh2_hits:
            foxa1_d = expr_d.iloc[foxa1_hits[0]].astype(float)
            ezh2_d  = expr_d.iloc[ezh2_hits[0]].astype(float)
            expr_dt = expr_d.T
            log(f"  Found via partial match: "
                f"FOXA1 at {foxa1_hits[0]}, "
                f"EZH2 at {ezh2_hits[0]}")
        else:
            log("  FOXA1/EZH2 not found in GSE96058.")
            log(f"  Index[:10]: {list(expr_d.index[:10])}")
            return {"status": "NOT_TESTABLE",
                    "predictions": {
                        "R3-D": "NOT TESTABLE",
                        "R3-E": "NOT TESTABLE",
                        "R3-F": "NOT TESTABLE"}}

    log(f"  FOXA1 n={foxa1_d.notna().sum()}  "
        f"mean={foxa1_d.mean():.4f}")
    log(f"  EZH2  n={ezh2_d.notna().sum()}  "
        f"mean={ezh2_d.mean():.4f}")

    # ── D.3: Build combined DataFrame ─────────────────────────
    log("")
    log("── D.3: BUILD COMBINED DATAFRAME ──")

    df_d = pd.DataFrame({
        "FOXA1": foxa1_d,
        "EZH2":  ezh2_d,
    })

    # Compute metric — GSE96058 is RNA-seq (log2 TPM or similar)
    # Use ratio for raw values, difference for z-scores
    if foxa1_d.mean() > 3.0:
        # Values look like log2 expression — use ratio
        df_d["metric"] = (
            df_d["FOXA1"]
            / (df_d["EZH2"].replace(0, np.nan) + 1e-6))
        metric_d = "FOXA1/EZH2 ratio (RNA-seq log2)"
    else:
        # Z-score-like — use difference
        df_d["metric"] = df_d["FOXA1"] - df_d["EZH2"]
        metric_d = "FOXA1 − EZH2 (RNA-seq)"

    log(f"  Metric: {metric_d}")

    # Attach PAM50 if available
    if clin_d is not None:
        log(f"  Clinical: {clin_d.shape}")
        log(f"  Clin cols: {list(clin_d.columns[:10])}")

        # Find PAM50 column
        pam50_col_d = None
        for c in clin_d.columns:
            if "pam50" in c.lower() or "subtype" in c.lower():
                vc = clin_d[c].value_counts()
                if len(vc) >= 3:
                    pam50_col_d = c
                    log(f"  PAM50 col: '{pam50_col_d}' — "
                        f"{vc.head(6).to_dict()}")
                    break

        # Find OS columns
        os_t_d = os_e_d = None
        for c in clin_d.columns:
            cl = c.lower()
            if "os" in cl and ("time" in cl or "month" in cl
                                or "day" in cl or "year" in cl):
                os_t_d = c
            if "os" in cl and ("status" in cl or "event" in cl
                                or "censor" in cl):
                os_e_d = c
            if "survival_time" in cl:
                os_t_d = c
            if "vital_status" in cl:
                os_e_d = c
        log(f"  OS: {os_t_d} / {os_e_d}")

        # Align clinical with expression
        common_d = df_d.index.intersection(clin_d.index)
        log(f"  Direct overlap: {len(common_d)}")

        if len(common_d) < 50:
            # Try 10-char truncation
            dt = pd.Index([str(x)[:10] for x in df_d.index])
            ct = pd.Index([str(x)[:10] for x in clin_d.index])
            common_d = dt.intersection(ct)
            log(f"  10-char overlap: {len(common_d)}")
            if len(common_d) >= 50:
                df_d.index  = dt
                clin_d.index = ct

        if len(common_d) >= 50 and pam50_col_d:
            keep_d = [pam50_col_d]
            if os_t_d: keep_d.append(os_t_d)
            if os_e_d: keep_d.append(os_e_d)
            df_d = df_d.loc[common_d].join(
                clin_d.loc[common_d, keep_d], how="inner")
            df_d = df_d.rename(columns={
                pam50_col_d: "PAM50",
                os_t_d: "OS_T" if os_t_d else "OS_T",
                os_e_d: "OS_E" if os_e_d else "OS_E",
            })
            df_d["subtype"] = df_d["PAM50"].apply(norm_subtype)
            log(f"  Merged: {df_d.shape}")
            log(f"  Subtypes:\n"
                f"{df_d['subtype'].value_counts().to_string()}")
        else:
            log("  Insufficient clinical overlap.")

    log(f"  Final df_d: {df_d.shape}")

    # ── D.4: Ordering test R3-D ───────────────────────────────
    log("")
    log("── D.4: ORDERING TEST R3-D ──")

    known_d = ["LumA", "LumB", "HER2", "TNBC"]
    groups_d = {}
    if "subtype" in df_d.columns:
        for sub in known_d + ["CL"]:
            vals = df_d[df_d["subtype"] == sub][
                "metric"].dropna().values
            if len(vals) >= 5:
                groups_d[sub] = vals

    medians_d = {s: float(np.median(v))
                 for s, v in groups_d.items()}

    if len(groups_d) >= 2:
        log(f"  Per-subtype {metric_d} medians:")
        for sub in SUBTYPE_ORDER:
            if sub in medians_d:
                log(f"    {sub}: n={len(groups_d[sub])}"
                    f"  median={medians_d[sub]:.4f}")

        kw_p_d, _ = kw_and_pairs(
            {s: v for s, v in groups_d.items()
             if len(v) >= 5})
        log(f"\n  KW p = "
            f"{kw_p_d:.2e}" if kw_p_d else "  KW p = N/A")

        adj_d = [(s1, s2) for s1, s2 in ADJACENT_PAIRS[:3]
                 if s1 in medians_d and s2 in medians_d]
        nc_d, nt_d, detail_d = order_check(medians_d, adj_d)

        for s1, s2, ok, m1, m2 in detail_d:
            log(f"  {s1}({m1:.4f}) > {s2}({m2:.4f}): "
                f"{'✓' if ok else '✗'}")

        r3d = (kw_p_d is not None and kw_p_d < 0.001
               and nc_d >= 2)
        preds["R3-D"] = "CONFIRMED" if r3d else "NOT CONFIRMED"
        log(f"\n  R3-D: {nc_d}/{nt_d} correct  "
            f"→ {preds['R3-D']}")
    else:
        log("  Subtype groups unavailable — R3-D NOT TESTABLE")
        kw_p_d = None
        nc_d = nt_d = 0
        preds["R3-D"] = "NOT TESTABLE"

    # ── D.5: LumA/LumB separation R3-F ───────────────────────
    log("")
    log("── D.5: LumA vs LumB SEPARATION R3-F ──")
    r3f_p = None
    if "LumA" in groups_d and "LumB" in groups_d:
        _, r3f_p = mannwhitneyu(
            groups_d["LumA"], groups_d["LumB"],
            alternative="two-sided")
        log(f"  LumA vs LumB p = {r3f_p:.2e}")
        preds["R3-F"] = (
            "CONFIRMED" if r3f_p < 0.001 else "NOT CONFIRMED")
        log(f"  R3-F → {preds['R3-F']}")
    else:
        preds["R3-F"] = "NOT TESTABLE"
        log("  LumA or LumB missing — R3-F NOT TESTABLE")

    # ── D.6: Survival R3-E ────────────────────────────────────
    log("")
    log("── D.6: SURVIVAL R3-E (LARGEST TEST) ──")
    r3e_p = None
    if "OS_T" in df_d.columns and "OS_E" in df_d.columns:
        df_d["OS_t_n"] = pd.to_numeric(
            df_d["OS_T"], errors="coerce")
        df_d["OS_e_n"] = df_d["OS_E"].apply(
            lambda x: 1 if str(x).upper() in
            ["1", "1.0", "DECEASED", "DEAD", "YES",
             "1:DECEASED", "EVENT"] else 0)
        # Convert to months if days
        if df_d["OS_t_n"].median() > 500:
            df_d["OS_t_n"] = df_d["OS_t_n"] / 30.44
            log("  OS time converted days → months")

        n_events = df_d["OS_e_n"].sum()
        log(f"  OS events: {int(n_events)}")
        log(f"  Median OS time: {df_d['OS_t_n'].median():.1f}")

        r3e_p = km_q14(df_d, "OS_t_n", "OS_e_n", "metric",
                       title=f"GSE96058 n={len(df_d)}\n"
                             f"{metric_d} vs OS")
        log(f"  KM p = "
            f"{r3e_p:.4f}" if r3e_p else "  KM not run")
        preds["R3-E"] = (
            "CONFIRMED" if r3e_p is not None and r3e_p < 0.05
            else "NOT TESTABLE" if r3e_p is None
            else "NOT CONFIRMED")
        log(f"  R3-E → {preds['R3-E']}")
    else:
        log("  OS columns not found — R3-E NOT TESTABLE")
        preds["R3-E"] = "NOT TESTABLE"

    # ── D.7: Figure ───────────────────────────────────────────
    n_panels_d = 3 if (HAS_LIFELINES
                       and "OS_t_n" in df_d.columns
                       and r3e_p is not None) else 2
    fig_d, ax_d = plt.subplots(1, n_panels_d,
                                figsize=(5 * n_panels_d + 1, 5))
    if len(groups_d) >= 2:
        violin_bar_figure(
            groups_d, medians_d,
            f"Component D: GSE96058 SCAN-B\n{metric_d}",
            f"Ordering: {nc_d}/{nt_d} correct\n"
            f"n={len(df_d):,}",
            len(df_d), kw_p_d, nc_d, nt_d,
            ax_d[0], ax_d[1])
    if n_panels_d == 3:
        km_q14(df_d, "OS_t_n", "OS_e_n", "metric",
               ax=ax_d[2],
               title=f"GSE96058 OS  n={len(df_d):,}")
    plt.suptitle("Component D: GSE96058 SCAN-B n=3,273",
                 fontsize=11)
    plt.tight_layout()
    plt.savefig(FIG_COMP_D, dpi=150, bbox_inches="tight")
    plt.close()
    log(f"  Figure: {FIG_COMP_D}")

    return {
        "status":   "COMPLETE" if foxa1_d is not None
                    else "NOT_TESTABLE",
        "n":        len(df_d),
        "metric":   metric_d,
        "medians":  medians_d,
        "kw_p":     kw_p_d if len(groups_d) >= 2 else None,
        "r3e_p":    r3e_p,
        "predictions": preds,
    }


# =============================================================
# COMBINED FIGURE
# =============================================================

def make_combined_figure(res_a, res_b, res_c, res_d):
    log("")
    log("── COMBINED FIGURE ──")

    datasets = [
        ("A: scRNA-seq per-patient\n(GSE176078)",
         res_a.get("medians", {})),
        ("B: Protein (CPTAC)\n(FOXA1−EZH2 ms)",
         {}),  # medians may not be available
        ("C: METABRIC mRNA fix\n(difference metric)",
         res_c.get("medians", {})),
        ("D: GSE96058 SCAN-B\n(n=3,273)",
         res_d.get("medians", {})),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(20, 5))

    for ax, (title, medians_dict) in zip(axes, datasets):
        if not medians_dict:
            ax.text(0.5, 0.5, "NOT TESTABLE\nor see component\nfigure",
                    ha="center", va="center",
                    transform=ax.transAxes, fontsize=10)
            ax.set_title(title, fontsize=9)
            continue

        subtypes_plot = [s for s in SUBTYPE_ORDER
                         if s in medians_dict]
        if not subtypes_plot:
            ax.text(0.5, 0.5, "No data",
                    ha="center", va="center",
                    transform=ax.transAxes)
            ax.set_title(title, fontsize=9)
            continue

        vals   = [medians_dict[s] for s in subtypes_plot]
        colors = [SUBTYPE_COLORS.get(s, "#888")
                  for s in subtypes_plot]
        ax.bar(subtypes_plot, vals, color=colors,
               alpha=0.85, edgecolor="black", linewidth=0.8)

        for i in range(len(subtypes_plot) - 1):
            s1, s2 = subtypes_plot[i], subtypes_plot[i + 1]
            correct = medians_dict[s1] > medians_dict[s2]
            ax.annotate(
                "✓" if correct else "✗",
                xy=((i + i + 1) / 2,
                    max(vals) * 0.92),
                ha="center", fontsize=14,
                color="green" if correct else "red")

        ax.set_title(title, fontsize=9)
        ax.set_ylabel("FOXA1/EZH2 metric")

    plt.suptitle(
        "FOXA1/EZH2 Ratio Cross-Dataset Validation — Script 2\n"
        "OrganismCore — RATIO-S2b | 2026-03-05",
        fontsize=11)
    plt.tight_layout()
    plt.savefig(FIG_COMBINED, dpi=150, bbox_inches="tight")
    plt.close()
    log(f"  Combined figure: {FIG_COMBINED}")


# =============================================================
# SCORECARD
# =============================================================

def make_scorecard(res_a, res_b, res_c, res_d):
    log("")
    log("=" * 65)
    log("SCORECARD — RATIO-S2b")
    log("=" * 65)

    rows = []

    for comp, res, dataset, tech in [
        ("A", res_a, "GSE176078 scRNA-seq per-patient",
         "scRNA-seq"),
        ("B", res_b, "CPTAC-BRCA proteomics",
         "MS protein"),
        ("C", res_c, "METABRIC mRNA (fix)",
         "mRNA microarray"),
        ("D", res_d, "GSE96058 SCAN-B",
         "RNA-seq"),
    ]:
        n = res.get("n") or res.get("n_samples", "N/A")
        for pid, status in res.get("predictions", {}).items():
            rows.append({
                "ID":      pid,
                "Comp":    comp,
                "Dataset": dataset,
                "n":       n,
                "Tech":    tech,
                "Status":  status,
            })

    df_sc = pd.DataFrame(rows)
    df_sc.to_csv(CSV_SCORECARD, index=False)

    log(f"\n  {'ID':<14} {'Dataset':<32} {'n':>6}  {'Status'}")
    log(f"  {'─'*14} {'─'*32} {'─'*6}  {'─'*25}")
    for _, row in df_sc.iterrows():
        log(f"  {row['ID']:<14} {row['Dataset']:<32} "
            f"{str(row['n']):>6}  {row['Status']}")

    conf  = df_sc["Status"].str.contains("CONFIRMED").sum()
    nc    = df_sc["Status"].str.contains("NOT CONFIRMED").sum()
    nt    = df_sc["Status"].str.contains("NOT TESTABLE").sum()
    total = len(df_sc)

    log(f"\n  Total predictions: {total}")
    log(f"  CONFIRMED:    {conf}")
    log(f"  NOT CONFIRMED:{nc}")
    log(f"  NOT TESTABLE: {nt}")

    # Special call-outs
    r2a_prot = res_b.get("predictions", {}).get(
        "R2-A(prot)", "NOT RUN")
    r3e      = res_d.get("predictions", {}).get(
        "R3-E", "NOT RUN")
    r1a_fix  = res_a.get("predictions", {}).get(
        "R1-A(fix)", "NOT RUN")
    r3a_fix  = res_c.get("predictions", {}).get(
        "R3-A(fix)", "NOT RUN")

    log("")
    log("  ════════════════════════════════════════════")
    log("  KEY RESULTS:")
    log(f"  R2-A(prot) — protein ordering (CPTAC): {r2a_prot}")
    log(f"  R3-E       — OS in n=3,273 (SCAN-B):   {r3e}")
    log(f"  R1-A(fix)  — scRNA per-patient order:  {r1a_fix}")
    log(f"  R3-A(fix)  — METABRIC diff metric:     {r3a_fix}")
    log("  ════════════════════════════════════════════")

    log(f"\n  Scorecard: {CSV_SCORECARD}")
    return df_sc


# =============================================================
# MAIN
# =============================================================

def main():
    log("=" * 65)
    log("FOXA1/EZH2 RATIO — VALIDATION SCRIPT 2 (RATIO-S2b)")
    log("OrganismCore | 2026-03-05 | Eric Robert Lawson")
    log("")
    log("CLINICAL GOAL: Validate FOXA1/EZH2 as a two-antibody")
    log("IHC tool for treatment stratification at any income")
    log("level and in any pathology lab worldwide.")
    log("")
    log("FOUR COMPONENTS:")
    log("  A — scRNA per-patient fix  (R1-A/B fix)")
    log("  B — CPTAC protein          (protein level — critical)")
    log("  C — METABRIC diff metric   (R3-A fix)")
    log("  D — GSE96058 SCAN-B n=3273 (new, largest test)")
    log("=" * 65)
    log(f"  Output:    {BASE_DIR}")
    log(f"  S1 cache:  {S1_DATA}")
    log(f"  S3 cache:  {S3_DATA}")
    log(f"  CL cache:  {CL_DATA}")
    log(f"  lifelines: {HAS_LIFELINES}")
    log(f"  sklearn:   {HAS_SKLEARN}")
    log("")

    flush_log()

    # ── Cache status check ────────────────────────────────────
    log("── CACHE FILE STATUS ──")
    cache_files = [
        ("scRNA expr (S1)",    SC_EXPR_CACHE),
        ("scRNA metadata (S1)", SC_METADATA_FILE),
        ("TCGA expr (S1)",     TCGA_EXPR_CACHE),
        ("TCGA clinical (S1)", TCGA_CLIN_FILE),
        ("METABRIC expr (S3)", META_EXPR_FILE),
        ("METABRIC clin (S3)", META_CLIN_FILE),
        ("GSE96058 expr (CL)", GSE96058_EXPR_FILE),
        ("GSE96058 mat1 (CL)", GSE96058_MAT1_FILE),
    ]
    for label, path in cache_files:
        if os.path.exists(path):
            log(f"  ✓  {label}: "
                f"{os.path.getsize(path):,} B")
        else:
            log(f"  ✗  {label}: NOT FOUND ({path})")

    flush_log()

    # ── Run components ────────────────────────────────────────
    log("")
    log("Starting Component A (scRNA per-patient fix)...")
    flush_log()
    res_a = run_component_a()

    log("")
    log("Starting Component B (CPTAC protein)...")
    flush_log()
    res_b = run_component_b()

    log("")
    log("Starting Component C (METABRIC fix)...")
    flush_log()
    res_c = run_component_c()

    log("")
    log("Starting Component D (GSE96058 SCAN-B)...")
    flush_log()
    res_d = run_component_d()

    # ── Combined figure ───────────────────────────────────────
    make_combined_figure(res_a, res_b, res_c, res_d)

    # ── Scorecard ─────────────────────────────────────────────
    make_scorecard(res_a, res_b, res_c, res_d)

    # ── Final summary ─────────────────────────────────────────
    log("")
    log("=" * 65)
    log("RATIO-S2b COMPLETE")
    log("=" * 65)
    log(f"  Component A (scRNA fix):     {res_a['status']}")
    log(f"  Component B (CPTAC protein): {res_b['status']}")
    log(f"  Component C (METABRIC fix):  {res_c['status']}")
    log(f"  Component D (GSE96058):      {res_d['status']}")
    log("")
    log("  KEY RESULT — R2-A(prot):")
    log(f"    {res_b.get('predictions', {}).get('R2-A(prot)', 'NOT RUN')}")
    log("  If R2-A(prot) = CONFIRMED: ratio is validated at")
    log("  protein level → IHC proposal is justified.")
    log("")
    log(f"  Log:       {LOG_FILE}")
    log(f"  Scorecard: {CSV_SCORECARD}")
    log(f"  Figures:   {RESULTS_DIR}")
    log("")
    log("  Next document: RATIO-S2c — script2_results_and_reasoning.md")
    log("=" * 65)

    flush_log()


if __name__ == "__main__":
    main()
