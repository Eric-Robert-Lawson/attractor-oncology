"""
FOXA1/EZH2 RATIO — VALIDATION SCRIPT 1 (RATIO-S1b) v2
OrganismCore — Document RATIO-S1b
Three-Component Public Data Validation

PREDICTIONS LOCKED IN: RATIO-S1a (before_script_ratio.md)
Date: 2026-03-05
Author: Eric Robert Lawson / OrganismCore

CHANGES FROM v1:
  FIX 1 — PATH CORRECTION
    Script lives in DEEP_DIVE/ (not a subdirectory).
    S1_DATA and S3_DATA must NOT use "../".
    Correct:
      S1_DATA = SCRIPT_DIR/Cross_Subtype_s1_results/data
      S3_DATA = SCRIPT_DIR/Cross_Subtype_s3_results/data

  FIX 2 — RPPA ANTIBODY NAME DETECTION
    Xena RPPA_RBN uses MD Anderson display names.
    FOXA1 is stored as "FOXA1-R-V", "FOXA1_V", or similar.
    EZH2  is stored as "EZH2-R-V", "EZH2_V", or similar.
    Updated detection: substring search on index AND columns,
    case-insensitive, for any token starting with FOXA1 / EZH2.
    Also added TCPA Level 4 TSV as second download candidate.

  FIX 3 — fetch_metabric_clinical() StopIteration
    rows_dict is keyed by sampleId (the outer key).
    pd.DataFrame.from_dict(rows_dict, orient="index") correctly
    places sampleId as the index. Previous code tried to find
    "SAMPLE_ID"/"sampleId" as a COLUMN — they were never columns.

COMPONENT A — GSE176078 scRNA-seq (per-patient concordance)
  Reuses cached data from Cross_Subtype_s1_results/data/
  Formal ordering test: LumA > LumB > HER2 > TNBC > CL
  Cohen's kappa vs pre-specified cut-points

COMPONENT B — TCGA RPPA protein level
  Downloads FOXA1 + EZH2 protein from TCPA/Xena
  PAM50 labels from TCGA clinical cache (S1) or cBioPortal
  Primary prediction: protein ratio ordering confirmed

COMPONENT C — METABRIC mRNA (large-cohort)
  Reuses or re-fetches metabric_expression.csv +
  metabric_clinical.csv from Cross_Subtype_s3_results/data/
  n=1,980 | OS endpoint | CLAUDIN_SUBTYPE column
"""

import os
import sys
import gzip
import json
import time
import warnings
import urllib.request
import urllib.error
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

# ============================================================
# PATHS
# ============================================================
# Script lives directly in DEEP_DIVE/ — same level as the
# Cross_Subtype_s*_results/ directories.
# No "../" needed.
# ============================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Prior script output directories (siblings of this script)
S1_DATA    = os.path.join(SCRIPT_DIR,
                           "Cross_Subtype_s1_results", "data")
S1_RESULTS = os.path.join(SCRIPT_DIR,
                           "Cross_Subtype_s1_results", "results")
S2_DATA    = os.path.join(SCRIPT_DIR,
                           "Cross_Subtype_s2_results", "data")
S2_RESULTS = os.path.join(SCRIPT_DIR,
                           "Cross_Subtype_s2_results", "results")
S3_DATA    = os.path.join(SCRIPT_DIR,
                           "Cross_Subtype_s3_results", "data")
S3_RESULTS = os.path.join(SCRIPT_DIR,
                           "Cross_Subtype_s3_results", "results")

# Output for this script
BASE_DIR    = os.path.join(SCRIPT_DIR, "ratio_s1_results")
DATA_DIR    = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
for d in [BASE_DIR, DATA_DIR, RESULTS_DIR]:
    os.makedirs(d, exist_ok=True)

# ── Cache files from prior scripts ────────────────────────────
# Named explicitly in Script 2 docstring:
TCGA_EXPR_CACHE = os.path.join(S1_RESULTS,
                                "expr_cache_cs_s1_tcga.csv")
TCGA_CLIN_FILE  = os.path.join(S1_DATA,
                                "TCGA_BRCA_clinicalMatrix.tsv")
SC_EXPR_CACHE   = os.path.join(S1_RESULTS,
                                "expr_cache_cs_s1_sc.csv")

# Named in Script 3 PATH block:
META_EXPR_FILE  = os.path.join(S3_DATA,
                                "metabric_expression.csv")
META_CLIN_FILE  = os.path.join(S3_DATA,
                                "metabric_clinical.csv")
GSE25066_EXPR   = os.path.join(S3_DATA,
                                "gse25066_expression.csv")
GSE25066_CLIN   = os.path.join(S3_DATA,
                                "gse25066_clinical.csv")

# ── RPPA download targets ──────────────────────────────────────
# Priority 1: Xena RPPA_RBN (131 antibodies × 747 samples)
# Priority 2: TCPA Level 4 (cleaner gene symbols as row labels)
# Priority 3: UCSC Xena pancanatlas
RPPA_GZ   = os.path.join(DATA_DIR, "RPPA_RBN.gz")
RPPA_FILE = os.path.join(DATA_DIR, "RPPA_RBN.tsv")

RPPA_URLS = [
    "https://tcga.xenahubs.net/download/"
    "TCGA.BRCA.sampleMap/RPPA_RBN.gz",
    "https://pancanatlas.xenahubs.net/download/"
    "TCGA.BRCA.sampleMap/RPPA_RBN.gz",
    # TCPA Level4 — gene symbols as rows, samples as columns
    "https://tcga.xenahubs.net/download/"
    "TCGA.BRCA.sampleMap/"
    "RPPA_protein_level_20160128.gz",
]

# cBioPortal API (same config as Script 3)
CBIO_BASE       = "https://www.cbioportal.org/api"
CBIO_STUDY      = "brca_metabric"
CBIO_PROFILE    = ("brca_metabric_mrna_median_all"
                   "_sample_Zscores")
CBIO_SAMPLELIST = "brca_metabric_mrna"

# TCGA PAM50 / clinical cache
TCGA_CBIO_STUDY  = "brca_tcga_pan_can_atlas_2018"
TCGA_PAM50_CACHE = os.path.join(DATA_DIR,
                                 "tcga_pam50_cbio.tsv")

# ── Outputs ────────────────────────────────────────────────────
LOG_FILE      = os.path.join(RESULTS_DIR, "ratio_s1_log.txt")
FIG_COMP_A    = os.path.join(RESULTS_DIR,
                              "ratio_s1_compA_scrna.png")
FIG_COMP_B    = os.path.join(RESULTS_DIR,
                              "ratio_s1_compB_rppa.png")
FIG_COMP_C    = os.path.join(RESULTS_DIR,
                              "ratio_s1_compC_metabric.png")
FIG_COMBINED  = os.path.join(RESULTS_DIR,
                              "ratio_s1_combined.png")
CSV_SCORECARD = os.path.join(RESULTS_DIR,
                              "ratio_s1_scorecard.csv")

# ── Pre-specified ordering / cut-points ───────────────────────
SUBTYPE_ORDER  = ["LumA", "LumB", "HER2", "TNBC", "CL"]
SUBTYPE_COLORS = {
    "LumA":  "#2166ac",
    "LumB":  "#74add1",
    "HER2":  "#f4a582",
    "TNBC":  "#d6604d",
    "Basal": "#d6604d",
    "CL":    "#762a83",
}

# ============================================================
# LOGGING
# ============================================================

_log = []

def log(msg=""):
    print(msg)
    _log.append(str(msg))

def flush_log():
    with open(LOG_FILE, "w") as f:
        f.write("\n".join(_log) + "\n")

# ============================================================
# UTILITIES
# ============================================================

def fetch_url(url, dest, label, retries=3, timeout=120):
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


def try_urls(url_list, dest_gz, label):
    if os.path.exists(dest_gz) and os.path.getsize(dest_gz) > 500:
        log(f"  {label}: cached")
        return dest_gz
    for url in url_list:
        tmp = dest_gz + ".tmp"
        r = fetch_url(url, tmp, label)
        if r and os.path.exists(tmp):
            os.rename(tmp, dest_gz)
            return dest_gz
    return None


def norm_subtype(s):
    s = str(s).strip()
    luma  = ["LumA", "lumA", "Luminal_A", "Luminal A",
             "LUMINAL_A", "luminal_a"]
    lumb  = ["LumB", "lumB", "Luminal_B", "Luminal B",
             "LUMINAL_B", "luminal_b"]
    her2  = ["HER2", "Her2", "her2", "HER2_enriched",
             "HER2-enriched", "HER2 ENRICHED", "Her2 SC"]
    tnbc  = ["Basal", "TNBC", "basal-like", "Basal-like",
             "Basal SC", "Cancer Basal SC", "TNBC/Basal",
             "Basal-Like", "BASAL"]
    cl    = ["CL", "claudin-low", "Claudin-low",
             "claudin_low", "Claudin_low", "CLAUDIN_LOW"]
    for label, aliases in [("LumA", luma), ("LumB", lumb),
                            ("HER2", her2), ("TNBC", tnbc),
                            ("CL", cl)]:
        if s in aliases:
            return label
    return s


def kw_and_pairs(groups):
    arrs = [v for k, v in groups.items() if len(v) >= 3]
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


def order_check(medians, pairs):
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
                label=f"Q1 low ratio (n={len(lo)})")
        kmf.plot_survival_function(ax=ax, ci_show=False,
                                   color="#d6604d")
        kmf.fit(hi[t_col], hi[e_col],
                label=f"Q4 high ratio (n={len(hi)})")
        kmf.plot_survival_function(ax=ax, ci_show=False,
                                   color="#2166ac")
        ax.set_title(f"{title}\np={p:.4f}", fontsize=9)
        ax.set_xlabel("Time")
        ax.set_ylabel("Survival")
        ax.legend(fontsize=7)
    return p


# ============================================================
# METABRIC FETCH — mirrors Script 3 exactly
# ============================================================

def fetch_metabric_expression(genes):
    """Bulk expression from cBioPortal — same endpoint as
    Script 3 (molecular-profiles / query / bulk)."""
    log("  Fetching METABRIC samples from cBioPortal...")
    samples = []
    for page in range(5):
        url = (f"{CBIO_BASE}/studies/{CBIO_STUDY}"
               f"/samples?pageSize=5000&pageNumber={page}")
        try:
            req = urllib.request.Request(
                url,
                headers={"Accept": "application/json",
                         "User-Agent": "OrganismCore/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                data = json.loads(r.read().decode())
            if not data:
                break
            samples.extend([s["sampleId"] for s in data])
            if len(data) < 5000:
                break
        except Exception as e:
            log(f"  Sample page {page} error: {e}")
            break

    log(f"  METABRIC samples fetched: {len(samples)}")
    if not samples:
        return None

    # ── Try per-gene molecular data endpoint ──────────────────
    # Script 3 used the v2 single-gene endpoint that worked.
    # The /query/bulk 404 is a v2 API change — fall back to
    # single-gene /molecular-data endpoint.
    records = []
    for gene in genes:
        url = (f"{CBIO_BASE}/molecular-profiles/"
               f"{CBIO_PROFILE}/genes/{gene}/molecular-data"
               f"?sampleListId={CBIO_SAMPLELIST}")
        try:
            req = urllib.request.Request(
                url,
                headers={"Accept": "application/json",
                         "User-Agent": "OrganismCore/1.0"})
            with urllib.request.urlopen(req, timeout=120) as r:
                data = json.loads(r.read().decode())
            for rec in data:
                records.append({
                    "gene":   gene,
                    "sample": rec.get("sampleId"),
                    "value":  rec.get("value"),
                })
            log(f"  {gene}: {len(data)} values via "
                "/molecular-data")
        except Exception as e1:
            log(f"  {gene}: /molecular-data failed — {e1}")
            # Second fallback: /molecular-data with POST body
            url2 = (f"{CBIO_BASE}/molecular-profiles/"
                    f"{CBIO_PROFILE}/molecular-data/fetch")
            body = json.dumps({
                "sampleIds":      samples[:5000],
                "entrezGeneIds":  [],
                "hugoGeneSymbols": [gene],
            }).encode()
            try:
                req2 = urllib.request.Request(
                    url2, data=body,
                    headers={
                        "Content-Type": "application/json",
                        "Accept":       "application/json",
                        "User-Agent":   "OrganismCore/1.0",
                    })
                with urllib.request.urlopen(req2,
                                            timeout=120) as r2:
                    data2 = json.loads(r2.read().decode())
                for rec in data2:
                    records.append({
                        "gene":   gene,
                        "sample": rec.get("sampleId"),
                        "value":  rec.get("value"),
                    })
                log(f"  {gene}: {len(data2)} values via "
                    "/molecular-data/fetch")
            except Exception as e2:
                log(f"  {gene}: /molecular-data/fetch "
                    f"also failed — {e2}")

    if not records:
        return None

    df = (pd.DataFrame(records)
          .pivot(index="gene", columns="sample",
                 values="value"))
    log(f"  Expression pivot: {df.shape}")
    return df


def fetch_metabric_clinical():
    """
    Paginated clinical fetch — mirrors Script 3.

    FIX 3: rows_dict is keyed by sampleId (the outer key).
    Use pd.DataFrame.from_dict(orient='index') which
    naturally sets the outer key as the DataFrame index.
    The old code tried to find 'sampleId' as a column name —
    it was never a column, it was always the row key.
    """
    log("  Fetching METABRIC clinical from cBioPortal...")
    rows_dict = {}
    for page in range(10):
        url = (f"{CBIO_BASE}/studies/{CBIO_STUDY}"
               f"/clinical-data"
               f"?clinicalDataType=SAMPLE"
               f"&pageSize=5000&pageNumber={page}")
        try:
            req = urllib.request.Request(
                url,
                headers={"Accept": "application/json",
                         "User-Agent": "OrganismCore/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                data = json.loads(r.read().decode())
            if not data:
                break
            for rec in data:
                sid = rec.get("sampleId", "")
                key = rec.get("clinicalAttributeId", "")
                val = rec.get("value", "")
                if sid not in rows_dict:
                    rows_dict[sid] = {}
                rows_dict[sid][key] = val
            log(f"  Clinical page {page}: {len(data)} records")
            if len(data) < 5000:
                break
        except Exception as e:
            log(f"  Clinical page {page} error: {e}")
            break

    if not rows_dict:
        log("  Clinical fetch returned no records.")
        return None

    # ── FIX 3: orient='index' — outer key = sampleId = index ──
    clin = pd.DataFrame.from_dict(rows_dict, orient="index")
    clin.index.name = "SAMPLE_ID"
    log(f"  Clinical DataFrame: {clin.shape}")
    log(f"  Clinical columns: {list(clin.columns[:10])}")
    return clin


# ============================================================
# COMPONENT A — GSE176078 scRNA-seq
# ============================================================

def run_component_a():
    log("")
    log("=" * 65)
    log("COMPONENT A — GSE176078 scRNA-seq PER-PATIENT")
    log("Reuses expr_cache_cs_s1_sc.csv from S1 cache")
    log("=" * 65)

    preds = {}

    # ── A.1: Load scRNA-seq cache ─────────────────────────────
    log("")
    log("── A.1: LOAD scRNA-seq EXPRESSION CACHE ──")
    log(f"  Looking for: {SC_EXPR_CACHE}")

    if not os.path.exists(SC_EXPR_CACHE):
        log("  Cache not found — Component A NOT TESTABLE")
        log("  (Run BRCA_Cross_Subtype_Script1.py first to "
            "generate expr_cache_cs_s1_sc.csv)")
        return {"status": "NOT_TESTABLE",
                "predictions": {
                    "R1-A": "NOT TESTABLE",
                    "R1-B": "NOT TESTABLE",
                    "R1-C": "NOT TESTABLE",
                }}

    log(f"  Loading: {SC_EXPR_CACHE}")
    sc = pd.read_csv(SC_EXPR_CACHE, index_col=0)
    log(f"  scRNA cache shape: {sc.shape}")
    log(f"  Columns sample: {list(sc.columns[:12])}")

    for gene in ["FOXA1", "EZH2"]:
        if gene not in sc.columns:
            log(f"  {gene} missing from scRNA cache — "
                "NOT TESTABLE")
            return {"status": "NOT_TESTABLE",
                    "predictions": {"R1-A": "NOT TESTABLE"}}

    # ── A.2: Attach cell type labels ─────────────────────────
    log("")
    log("── A.2: ATTACH CELL TYPE LABELS ──")

    meta_paths = [
        os.path.join(S1_DATA, "sc_metadata.csv"),
        os.path.join(S1_DATA,
                     "GSE176078_Wu_2021_BRCA_scRNA_metadata.csv"),
        os.path.join(S1_DATA, "metadata.csv"),
        os.path.join(S1_DATA, "cell_metadata.csv"),
    ]
    meta = None
    for p in meta_paths:
        if os.path.exists(p):
            meta = pd.read_csv(p, index_col=0)
            log(f"  Metadata: {p}  shape={meta.shape}")
            break

    if meta is None:
        meta_url = (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/"
            "GSE176078/suppl/"
            "GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"
        )
        meta_gz    = os.path.join(DATA_DIR, "sc_meta.csv.gz")
        meta_plain = os.path.join(DATA_DIR, "sc_meta.csv")
        fetch_url(meta_url, meta_gz, "sc_metadata")
        if os.path.exists(meta_gz):
            import shutil
            with gzip.open(meta_gz, "rb") as fi, \
                 open(meta_plain, "wb") as fo:
                shutil.copyfileobj(fi, fo)
            meta = pd.read_csv(meta_plain, index_col=0)
            log(f"  Downloaded metadata: shape={meta.shape}")

    if meta is None:
        log("  No metadata found — Component A NOT TESTABLE")
        return {"status": "NOT_TESTABLE",
                "predictions": {"R1-A": "NOT TESTABLE"}}

    # ── A.3: Find subtype column ──────────────────────────────
    sub_col = None
    for c in meta.columns:
        vc = meta[c].value_counts()
        vals = str(vc.index.tolist()).upper()
        if (("LUMA" in vals or "LUM_A" in vals
             or "CANCER LUMA" in vals
             or "CANCER LUMB" in vals)
                and len(vc) >= 3):
            sub_col = c
            break

    if sub_col is None:
        for c in meta.columns:
            if ("subtype" in c.lower() or "pam50" in c.lower()
                    or "celltype" in c.lower()):
                vc = meta[c].value_counts()
                if len(vc) >= 3:
                    sub_col = c
                    break

    if sub_col is None:
        log("  Cannot identify subtype column in metadata.")
        log(f"  Columns: {list(meta.columns[:15])}")
        return {"status": "NOT_TESTABLE",
                "predictions": {"R1-A": "NOT TESTABLE"}}

    log(f"  Subtype column: '{sub_col}'")
    log(f"  Values:\n{meta[sub_col].value_counts().head(10)}")

    # ── A.4: Per-cell ratio ───────────────────────────────────
    log("")
    log("── A.4: PER-CELL FOXA1/EZH2 RATIO ──")

    common = sc.index.intersection(meta.index)
    log(f"  Cell barcode overlap: {len(common)}")
    if len(common) < 50:
        log("  Not enough overlap — NOT TESTABLE")
        return {"status": "NOT_TESTABLE",
                "predictions": {"R1-A": "NOT TESTABLE"}}

    sc2 = sc.loc[common].copy()
    sc2["subtype_raw"] = meta.loc[common, sub_col]
    sc2["subtype"]     = sc2["subtype_raw"].apply(norm_subtype)
    sc2["ratio"]       = (sc2["FOXA1"]
                          / (sc2["EZH2"] + 1e-6))

    known  = ["LumA", "LumB", "HER2", "TNBC", "CL"]
    df_a   = sc2[sc2["subtype"].isin(known)].copy()
    log(f"  Cells with known subtype: {len(df_a)}")
    log(f"  Per subtype:\n{df_a['subtype'].value_counts()}")

    if len(df_a) < 50:
        log("  Insufficient labelled cells — NOT TESTABLE")
        return {"status": "NOT_TESTABLE",
                "predictions": {"R1-A": "NOT TESTABLE"}}

    # ── A.5: Ordering test R1-A ───────────────────────────────
    log("")
    log("── A.5: ORDERING TEST (R1-A) ──")

    groups_a = {
        s: df_a[df_a["subtype"] == s]["ratio"].dropna().values
        for s in known
        if (df_a["subtype"] == s).sum() >= 3
    }
    medians_a = {s: float(np.median(v))
                 for s, v in groups_a.items()}

    log("  Ratio medians per subtype:")
    for s in SUBTYPE_ORDER:
        if s in medians_a:
            log(f"    {s:6s}: n={len(groups_a[s]):5d}  "
                f"med={medians_a[s]:.4f}")

    kw_p_a, pw_a = kw_and_pairs(groups_a)
    log(f"\n  KW p = {kw_p_a:.2e}")

    adj = [("LumA","LumB"), ("LumB","HER2"),
           ("HER2","TNBC"), ("TNBC","CL")]
    nc_a, nt_a, detail_a = order_check(medians_a, adj)
    for s1, s2, ok, m1, m2 in detail_a:
        log(f"  {s1}({m1:.4f}) > {s2}({m2:.4f}): "
            f"{'✓' if ok else '✗'}")

    r1a = (kw_p_a is not None and kw_p_a < 0.001
           and nc_a >= 3)
    preds["R1-A"] = "CONFIRMED" if r1a else "NOT CONFIRMED"
    log(f"  R1-A: {nc_a}/{nt_a} correct  "
        f"KW p={kw_p_a:.2e}  → {preds['R1-A']}")

    # ── A.6: Concordance R1-B ────────────────────────────────
    log("")
    log("── A.6: CONCORDANCE WITH CUTS (R1-B) ──")

    cuts = {"LumA": 8.0, "LumB": 5.0, "HER2": 1.0,
            "TNBC": 0.2}

    def apply_cuts(r):
        if pd.isna(r):
            return "Unknown"
        if r > 8.0:  return "LumA"
        if r > 5.0:  return "LumB"
        if r > 1.0:  return "HER2"
        if r > 0.2:  return "TNBC"
        return "CL"

    df_a["pred"] = df_a["ratio"].apply(apply_cuts)
    kappa_a = kappa(df_a["subtype"], df_a["pred"])
    log(f"  Cohen's kappa: {kappa_a}")
    r1b = (kappa_a is not None and kappa_a >= 0.25)
    preds["R1-B"] = "CONFIRMED" if r1b else "NOT CONFIRMED"
    log(f"  R1-B → {preds['R1-B']}")

    # ── A.7: Ratio vs single gene R1-C ───────────────────────
    log("")
    log("── A.7: RATIO OUTPERFORMS SINGLE GENE (R1-C) ──")

    f_med  = df_a["FOXA1"].median()
    e_med  = df_a["EZH2"].median()
    binary = df_a["subtype"].map(
        lambda x: "LumA" if x in ["LumA","LumB"] else "TNBC")
    pred_f = (df_a["FOXA1"] >= f_med).map(
        {True:"LumA", False:"TNBC"})
    pred_e = (df_a["EZH2"]  >= e_med).map(
        {True:"TNBC", False:"LumA"})
    kf = kappa(binary, pred_f)
    ke = kappa(binary, pred_e)
    log(f"  kappa ratio={kappa_a}  FOXA1={kf}  EZH2={ke}")
    r1c = (kappa_a is not None and kf is not None
           and ke is not None
           and kappa_a > kf and kappa_a > ke)
    preds["R1-C"] = "CONFIRMED" if r1c else "NOT CONFIRMED"
    log(f"  R1-C → {preds['R1-C']}")

    # ── A.8: Figure ───────────────────────────────────────────
    fig_a, axes_a = plt.subplots(1, 2, figsize=(12, 5))
    plot_s = [s for s in SUBTYPE_ORDER if s in groups_a]

    axes_a[0].violinplot(
        [groups_a[s] for s in plot_s], showmedians=True)
    axes_a[0].set_xticks(range(1, len(plot_s)+1))
    axes_a[0].set_xticklabels(plot_s)
    axes_a[0].set_ylabel("FOXA1/EZH2 (log1p scRNA)")
    axes_a[0].set_title(
        f"GSE176078  n={len(df_a):,} cells\n"
        f"KW p={kw_p_a:.2e}")

    meds = [medians_a.get(s, np.nan) for s in plot_s]
    cols = [SUBTYPE_COLORS.get(s, "#888888") for s in plot_s]
    axes_a[1].bar(plot_s, meds, color=cols, alpha=0.85,
                  edgecolor="black", linewidth=0.7)
    axes_a[1].set_ylabel("Median ratio")
    axes_a[1].set_title(
        f"Ordering: {nc_a}/{nt_a} correct\n"
        f"R1-A: {preds['R1-A']}")

    plt.tight_layout()
    plt.savefig(FIG_COMP_A, dpi=150, bbox_inches="tight")
    plt.close()
    log(f"  Figure: {FIG_COMP_A}")

    return {
        "status":   "COMPLETE",
        "n_cells":  len(df_a),
        "medians":  medians_a,
        "kw_p":     kw_p_a,
        "nc": nc_a, "nt": nt_a,
        "kappa":    kappa_a,
        "predictions": preds,
    }


# ============================================================
# COMPONENT B — TCGA RPPA protein-level validation
# ============================================================

def load_rppa():
    """
    Load TCGA BRCA RPPA.
    Returns (foxa1_series, ezh2_series, rppa_df) or
            (None, None, None).

    FIX 2: Xena RPPA_RBN has proteins as ROWS
    (e.g., '4EBP1', 'FOXA1-R-V', 'EZH2-C-V') and
    samples as COLUMNS.
    Search the index (rows) with substring matching,
    case-insensitive, for tokens that START with FOXA1 / EZH2.
    """
    got = try_urls(RPPA_URLS, RPPA_GZ, "RPPA")
    rppa = None

    if got and os.path.exists(got):
        try:
            with gzip.open(got, "rt") as f:
                rppa = pd.read_csv(f, sep="\t", index_col=0)
            log(f"  RPPA loaded from gz: {rppa.shape}")
        except Exception as e1:
            try:
                rppa = pd.read_csv(got, sep="\t", index_col=0)
                log(f"  RPPA loaded as plain: {rppa.shape}")
            except Exception as e2:
                log(f"  RPPA parse error: {e1} / {e2}")

    if rppa is None and os.path.exists(RPPA_FILE):
        rppa = pd.read_csv(RPPA_FILE, sep="\t", index_col=0)
        log(f"  RPPA loaded from local TSV: {rppa.shape}")

    if rppa is None:
        return None, None, None

    log(f"  RPPA index[:5]:   {list(rppa.index[:5])}")
    log(f"  RPPA columns[:3]: {list(rppa.columns[:3])}")

    # ── Search ROWS (proteins) for FOXA1 and EZH2 ────────────
    # Xena RPPA_RBN format:
    #   rows = antibodies (e.g. 'FOXA1-R-V', 'EZH2-C-V')
    #   columns = TCGA barcodes
    foxa1_row = ezh2_row = None

    all_rows_upper = {str(r).upper(): str(r)
                      for r in rppa.index}
    for token, orig in all_rows_upper.items():
        if token.startswith("FOXA1"):
            foxa1_row = orig
            break
    for token, orig in all_rows_upper.items():
        if token.startswith("EZH2"):
            ezh2_row = orig
            break

    if foxa1_row and ezh2_row:
        log(f"  Found FOXA1 row: '{foxa1_row}'")
        log(f"  Found EZH2  row: '{ezh2_row}'")
        return (rppa.loc[foxa1_row].astype(float),
                rppa.loc[ezh2_row].astype(float),
                rppa)

    # ── Search COLUMNS (antibodies) as fallback ───────────────
    all_cols_upper = {str(c).upper(): str(c)
                      for c in rppa.columns}
    for token, orig in all_cols_upper.items():
        if token.startswith("FOXA1"):
            foxa1_row = orig
            break
    for token, orig in all_cols_upper.items():
        if token.startswith("EZH2"):
            ezh2_row = orig
            break

    if foxa1_row and ezh2_row:
        log(f"  Cols orientation: FOXA1='{foxa1_row}' "
            f"EZH2='{ezh2_row}'")
        rppa2 = rppa.T
        return (rppa2.loc[foxa1_row].astype(float),
                rppa2.loc[ezh2_row].astype(float),
                rppa2)

    # ── Report what IS available ──────────────────────────────
    log("  FOXA1/EZH2 not found in RPPA.")
    fox_rows = [x for x in rppa.index
                if "FOX" in str(x).upper()][:8]
    ezh_rows = [x for x in rppa.index
                if "EZH" in str(x).upper()][:8]
    log(f"  FOX* rows: {fox_rows}")
    log(f"  EZH* rows: {ezh_rows}")
    fox_cols = [x for x in rppa.columns
                if "FOX" in str(x).upper()][:8]
    log(f"  FOX* cols: {fox_cols}")

    # ── Save full index for inspection ───────────────────────
    idx_file = os.path.join(DATA_DIR, "rppa_row_index.txt")
    with open(idx_file, "w") as f:
        f.write("\n".join(str(x) for x in rppa.index))
    log(f"  Full RPPA row index saved to: {idx_file}")

    return None, None, rppa


def load_tcga_pam50():
    """
    Load TCGA PAM50 labels.
    1. Try TCGA_BRCA_clinicalMatrix.tsv from S1_DATA
    2. Try TCGA_PAM50_CACHE (prior cBioPortal fetch)
    3. Fetch from cBioPortal
    Returns (clin_df, pam50_col) or (None, None).
    """
    for src, path in [("S1 cache", TCGA_CLIN_FILE),
                      ("local cache", TCGA_PAM50_CACHE)]:
        if os.path.exists(path):
            sep = "\t" if path.endswith(".tsv") else ","
            try:
                clin = pd.read_csv(path, sep=sep, index_col=0,
                                   low_memory=False)
                log(f"  TCGA clinical from {src}: {clin.shape}")
                pam50_col = None
                for c in clin.columns:
                    cu = c.upper()
                    if "PAM50" in cu or "SUBTYPE" in cu:
                        vc = clin[c].value_counts()
                        labels = str(vc.index.tolist()).upper()
                        if (len(vc) >= 3
                                and ("LUMA" in labels
                                     or "BASAL" in labels)):
                            pam50_col = c
                            break
                if pam50_col:
                    log(f"  PAM50 col: '{pam50_col}'")
                    return clin, pam50_col
                else:
                    log(f"  No PAM50 col found in {src}.")
                    log(f"  Cols: {list(clin.columns[:20])}")
            except Exception as e:
                log(f"  {src} load error: {e}")

    # ── cBioPortal fetch ──────────────────────────────────────
    log("  Fetching TCGA clinical from cBioPortal...")
    rows_dict = {}
    for page in range(10):
        url = (f"{CBIO_BASE}/studies/{TCGA_CBIO_STUDY}"
               f"/clinical-data"
               f"?clinicalDataType=SAMPLE"
               f"&pageSize=5000&pageNumber={page}")
        try:
            req = urllib.request.Request(
                url,
                headers={"Accept": "application/json",
                         "User-Agent": "OrganismCore/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                data = json.loads(r.read().decode())
            if not data:
                break
            for rec in data:
                sid = rec.get("sampleId","")
                key = rec.get("clinicalAttributeId","")
                val = rec.get("value","")
                if sid not in rows_dict:
                    rows_dict[sid] = {}
                rows_dict[sid][key] = val
            log(f"  TCGA page {page}: {len(data)} records")
            if len(data) < 5000:
                break
        except Exception as e:
            log(f"  TCGA page {page} error: {e}")
            break

    if not rows_dict:
        return None, None

    # FIX 3 applied here as well
    clin = pd.DataFrame.from_dict(rows_dict, orient="index")
    clin.index.name = "SAMPLE_ID"
    clin.to_csv(TCGA_PAM50_CACHE, sep="\t")
    log(f"  TCGA clinical fetched: {clin.shape}")
    log(f"  Cols: {list(clin.columns[:15])}")

    pam50_col = None
    for c in clin.columns:
        cu = c.upper()
        if "PAM50" in cu or "SUBTYPE" in cu:
            vc = clin[c].value_counts()
            labels = str(vc.index.tolist()).upper()
            if len(vc) >= 3 and "LUMA" in labels:
                pam50_col = c
                break
    if pam50_col:
        log(f"  PAM50 col: '{pam50_col}'")
    return clin, pam50_col


def run_component_b():
    log("")
    log("=" * 65)
    log("COMPONENT B — TCGA RPPA PROTEIN-LEVEL VALIDATION")
    log("FOXA1 + EZH2 protein | TCGA BRCA | n≈750")
    log("=" * 65)

    preds = {}

    # ── B.1: Load RPPA ────────────────────────────────────────
    log("")
    log("── B.1: LOAD RPPA ──")
    foxa1_p, ezh2_p, rppa_full = load_rppa()

    if foxa1_p is None:
        log("  RPPA FOXA1/EZH2 not accessible.")
        log("  See rppa_row_index.txt for available proteins.")
        log("  Component B NOT TESTABLE")
        return {"status":      "NOT_TESTABLE",
                "predictions": {
                    "R2-A": "NOT TESTABLE",
                    "R2-B": "NOT TESTABLE",
                    "R2-C": "NOT TESTABLE",
                    "R2-D": "NOT TESTABLE",
                }}

    log(f"  FOXA1 protein n={foxa1_p.notna().sum()}")
    log(f"  EZH2  protein n={ezh2_p.notna().sum()}")

    # ── B.2: Load PAM50 ───────────────────────────────────────
    log("")
    log("── B.2: LOAD TCGA PAM50 ──")
    clin_b, pam50_col_b = load_tcga_pam50()

    if clin_b is None or pam50_col_b is None:
        log("  PAM50 unavailable — Component B NOT TESTABLE")
        return {"status":      "NOT_TESTABLE",
                "predictions": {"R2-A": "NOT TESTABLE"}}

    log(f"  PAM50 distribution:\n"
        f"{clin_b[pam50_col_b].value_counts().head(8)}")

    # OS columns
    os_t = os_e = None
    for c in clin_b.columns:
        cu = c.upper()
        if cu == "OS_MONTHS":    os_t = c
        if cu == "OS_STATUS":    os_e = c
    log(f"  OS time={os_t}  OS event={os_e}")

    # ── B.3: Align RPPA ↔ PAM50 ──────────────────────────────
    log("")
    log("── B.3: ALIGN RPPA ↔ PAM50 ──")

    # Try direct overlap
    common_b = foxa1_p.index.intersection(clin_b.index)
    log(f"  Full barcode overlap: {len(common_b)}")

    if len(common_b) < 20:
        # 15-char truncation (TCGA-XX-XXXX-01 → same)
        rppa_idx  = pd.Index(
            [str(x)[:15] for x in foxa1_p.index])
        clin_idx  = pd.Index(
            [str(x)[:15] for x in clin_b.index])
        common_b  = rppa_idx.intersection(clin_idx)
        log(f"  15-char overlap: {len(common_b)}")
        if len(common_b) >= 20:
            foxa1_p.index = rppa_idx
            ezh2_p.index  = rppa_idx
            clin_b.index  = clin_idx

    if len(common_b) < 20:
        log("  Insufficient overlap (<20) — NOT TESTABLE")
        return {"status":      "NOT_TESTABLE",
                "predictions": {"R2-A": "NOT TESTABLE"}}

    log(f"  Using {len(common_b)} matched samples")

    keep_cols = [pam50_col_b]
    if os_t: keep_cols.append(os_t)
    if os_e: keep_cols.append(os_e)

    df_b = pd.DataFrame({
        "FOXA1_p": foxa1_p[foxa1_p.index.isin(common_b)],
        "EZH2_p":  ezh2_p[ezh2_p.index.isin(common_b)],
    })
    df_b = df_b.join(
        clin_b.loc[clin_b.index.isin(common_b),
                   keep_cols],
        how="inner")
    rename_map = {pam50_col_b: "PAM50"}
    if os_t: rename_map[os_t] = "OS_MONTHS"
    if os_e: rename_map[os_e] = "OS_STATUS"
    df_b = df_b.rename(columns=rename_map)
    df_b["subtype"] = df_b["PAM50"].apply(norm_subtype)
    df_b["ratio_p"] = (
        df_b["FOXA1_p"]
        / (df_b["EZH2_p"].replace(0, np.nan) + 1e-6))

    log(f"  df_b: {df_b.shape}")
    log(f"  Subtype dist:\n"
        f"{df_b['subtype'].value_counts()}")

    # ── B.4: Ordering test R2-A (PRIMARY) ─────────────────────
    log("")
    log("── B.4: R2-A ★ PRIMARY ★ ──")

    known_b  = ["LumA", "LumB", "HER2", "TNBC"]
    groups_b = {
        s: df_b[df_b["subtype"] == s]["ratio_p"].dropna().values
        for s in known_b
        if (df_b["subtype"] == s).sum() >= 5
    }
    medians_b = {s: float(np.median(v))
                 for s, v in groups_b.items()}
    log("  Protein ratio medians:")
    for s in known_b:
        if s in medians_b:
            log(f"    {s:6s}: n={len(groups_b[s]):3d}  "
                f"med={medians_b[s]:.4f}")

    kw_p_b, pw_b = kw_and_pairs(groups_b)
    log(f"\n  KW p = {kw_p_b:.2e}" if kw_p_b else
        "  KW p = None")

    adj_b = [("LumA","LumB"), ("LumB","HER2"), ("HER2","TNBC")]
    nc_b, nt_b, detail_b = order_check(medians_b, adj_b)
    for s1, s2, ok, m1, m2 in detail_b:
        log(f"  {s1}({m1:.4f}) > {s2}({m2:.4f}): "
            f"{'✓' if ok else '✗'}")

    r2a = (kw_p_b is not None and kw_p_b < 0.001
           and nc_b >= 2)
    preds["R2-A"] = ("CONFIRMED ★ PRIMARY ★"
                     if r2a else "NOT CONFIRMED")
    log(f"  R2-A: {nc_b}/{nt_b}  KW p={kw_p_b}  "
        f"→ {preds['R2-A']}")

    # ── B.5: EZH2 highest in TNBC R2-B ───────────────────────
    ezh2_meds = {
        s: float(np.median(
            df_b[df_b["subtype"]==s]["EZH2_p"].dropna()))
        for s in known_b
        if (df_b["subtype"]==s).sum() >= 5
    }
    log(f"  EZH2 protein meds: {ezh2_meds}")
    r2b = ("TNBC" in ezh2_meds and "LumA" in ezh2_meds
           and ezh2_meds["TNBC"] > ezh2_meds["LumA"])
    preds["R2-B"] = "CONFIRMED" if r2b else "NOT CONFIRMED"
    log(f"  R2-B → {preds['R2-B']}")

    # ── B.6: FOXA1 highest in LumA R2-C ──────────────────────
    foxa1_meds = {
        s: float(np.median(
            df_b[df_b["subtype"]==s]["FOXA1_p"].dropna()))
        for s in known_b
        if (df_b["subtype"]==s).sum() >= 5
    }
    log(f"  FOXA1 protein meds: {foxa1_meds}")
    r2c = ("LumA" in foxa1_meds and "TNBC" in foxa1_meds
           and foxa1_meds["LumA"] > foxa1_meds["TNBC"])
    preds["R2-C"] = "CONFIRMED" if r2c else "NOT CONFIRMED"
    log(f"  R2-C → {preds['R2-C']}")

    # ── B.7: Survival R2-D ────────────────────────────────────
    log("")
    log("── B.7: SURVIVAL (R2-D) ──")
    r2d_p = None
    if "OS_MONTHS" in df_b.columns and \
       "OS_STATUS" in df_b.columns:
        df_b["OS_t"] = pd.to_numeric(
            df_b["OS_MONTHS"], errors="coerce")
        df_b["OS_e"] = df_b["OS_STATUS"].apply(
            lambda x: 1 if str(x).upper() in
            ["DECEASED","1","1.0","DEAD",
             "1:DECEASED","DEAD:DECEASED"]
            else 0)
        r2d_p = km_q14(df_b, "OS_t", "OS_e", "ratio_p",
                       title="RPPA ratio vs OS")
        log(f"  KM p = {r2d_p:.4f}" if r2d_p is not None
            else "  KM: not testable")
    preds["R2-D"] = (
        "CONFIRMED" if r2d_p is not None and r2d_p < 0.05
        else ("NOT TESTABLE" if r2d_p is None
              else "NOT CONFIRMED"))
    log(f"  R2-D → {preds['R2-D']}")

    # ── B.8: Figure ───────────────────────────────────────────
    plot_s_b = [s for s in known_b if s in groups_b]
    fig_b, ax_b = plt.subplots(1, 3, figsize=(15, 5))

    ax_b[0].violinplot(
        [df_b[df_b["subtype"]==s]["FOXA1_p"].dropna().values
         for s in plot_s_b], showmedians=True)
    ax_b[0].set_xticks(range(1, len(plot_s_b)+1))
    ax_b[0].set_xticklabels(plot_s_b)
    ax_b[0].set_title("FOXA1 protein (RPPA)\n"
                      "Expected: LumA>LumB>HER2>TNBC")
    ax_b[0].set_ylabel("RPPA level")

    ax_b[1].violinplot(
        [df_b[df_b["subtype"]==s]["EZH2_p"].dropna().values
         for s in plot_s_b], showmedians=True)
    ax_b[1].set_xticks(range(1, len(plot_s_b)+1))
    ax_b[1].set_xticklabels(plot_s_b)
    ax_b[1].set_title("EZH2 protein (RPPA)\n"
                      "Expected: TNBC>HER2>LumB>LumA")
    ax_b[1].set_ylabel("RPPA level")

    ax_b[2].violinplot(
        [groups_b[s] for s in plot_s_b], showmedians=True)
    ax_b[2].set_xticks(range(1, len(plot_s_b)+1))
    ax_b[2].set_xticklabels(plot_s_b)
    ax_b[2].set_title(
        f"FOXA1/EZH2 protein ratio\n"
        f"KW p={'N/A' if kw_p_b is None else f'{kw_p_b:.2e}'}"
        f"  R2-A: {preds['R2-A']}")
    ax_b[2].set_ylabel("FOXA1/EZH2 protein ratio")

    plt.suptitle(f"Component B: TCGA RPPA  n={len(df_b)}",
                 y=1.01)
    plt.tight_layout()
    plt.savefig(FIG_COMP_B, dpi=150, bbox_inches="tight")
    plt.close()
    log(f"  Figure: {FIG_COMP_B}")

    return {
        "status":    "COMPLETE",
        "n_samples": len(df_b),
        "medians":   medians_b,
        "kw_p":      kw_p_b,
        "nc": nc_b, "nt": nt_b,
        "r2d_p":     r2d_p,
        "predictions": preds,
    }


# ============================================================
# COMPONENT C — METABRIC mRNA large-cohort
# ============================================================

def load_metabric():
    """
    Load METABRIC expression and clinical.
    Expression: prefer Script 3 cache (genes × samples or
    samples × genes — both orientations handled).
    Clinical: prefer Script 3 cache, else fetch via cBioPortal
    using the corrected fetch_metabric_clinical().
    """
    expr = None
    clin = None

    # ── Expression ────────────────────────────────────────────
    for label, path in [
        ("S3 cache",  META_EXPR_FILE),
        ("S1 alt",    os.path.join(
                          S1_DATA, "metabric_expression.csv")),
    ]:
        if os.path.exists(path):
            try:
                expr = pd.read_csv(path, index_col=0)
                log(f"  METABRIC expr from {label}: "
                    f"{expr.shape}")
                break
            except Exception as e:
                log(f"  {label} expr error: {e}")

    if expr is None:
        log("  METABRIC expression not cached — "
            "fetching FOXA1 + EZH2 from cBioPortal...")
        expr = fetch_metabric_expression(["FOXA1", "EZH2"])
        if expr is not None:
            os.makedirs(os.path.dirname(META_EXPR_FILE),
                        exist_ok=True)
            expr.to_csv(META_EXPR_FILE)
            log(f"  Saved: {META_EXPR_FILE}")

    # ── Clinical ────────────────────────────────��─────────────
    for label, path in [
        ("S3 cache",  META_CLIN_FILE),
        ("S1 alt",    os.path.join(
                          S1_DATA, "metabric_clinical.csv")),
    ]:
        if os.path.exists(path):
            try:
                clin = pd.read_csv(path, index_col=0)
                log(f"  METABRIC clin from {label}: "
                    f"{clin.shape}")
                break
            except Exception as e:
                log(f"  {label} clin error: {e}")

    if clin is None:
        log("  METABRIC clinical not cached — fetching...")
        clin = fetch_metabric_clinical()   # FIX 3 applied here
        if clin is not None:
            os.makedirs(os.path.dirname(META_CLIN_FILE),
                        exist_ok=True)
            clin.to_csv(META_CLIN_FILE)
            log(f"  Saved: {META_CLIN_FILE}")

    return expr, clin


def run_component_c():
    log("")
    log("=" * 65)
    log("COMPONENT C — METABRIC mRNA LARGE-COHORT")
    log("n=1,980 | CLAUDIN_SUBTYPE | OS endpoint")
    log("=" * 65)

    preds = {}

    # ── C.1: Load METABRIC ────────────────────────────────────
    log("")
    log("── C.1: LOAD METABRIC DATA ──")
    expr_c, clin_c = load_metabric()

    if expr_c is None:
        log("  Expression unavailable — NOT TESTABLE")
        return {"status": "NOT_TESTABLE", "predictions": {}}
    if clin_c is None:
        log("  Clinical unavailable — NOT TESTABLE")
        return {"status": "NOT_TESTABLE", "predictions": {}}

    log(f"  Expression: {expr_c.shape}")
    log(f"  Clinical:   {clin_c.shape}")
    log(f"  Clin cols:  {list(clin_c.columns[:15])}")

    # ── C.2: Locate FOXA1 and EZH2 ───────────────────────────
    log("")
    log("── C.2: LOCATE FOXA1 AND EZH2 ──")

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
    else:
        log("  FOXA1/EZH2 absent — NOT TESTABLE")
        log(f"  Index[:5]:   {list(expr_c.index[:5])}")
        log(f"  Columns[:5]: {list(expr_c.columns[:5])}")
        return {"status": "NOT_TESTABLE", "predictions": {}}

    # ── C.3: Subtype and OS columns ───────────────────────────
    log("")
    log("── C.3: FIND SUBTYPE + OS COLUMNS ──")

    sub_col_c = None
    for c in clin_c.columns:
        if "CLAUDIN_SUBTYPE" in c.upper():
            sub_col_c = c
            break
    if sub_col_c is None:
        for c in clin_c.columns:
            cu = c.upper()
            if "SUBTYPE" in cu or "PAM50" in cu:
                vc = clin_c[c].value_counts()
                if len(vc) >= 3:
                    sub_col_c = c
                    break

    if sub_col_c is None:
        log("  No subtype column found — NOT TESTABLE")
        log(f"  Cols: {list(clin_c.columns)}")
        return {"status": "NOT_TESTABLE", "predictions": {}}

    log(f"  Subtype col: '{sub_col_c}'")
    log(f"  Distribution:\n"
        f"{clin_c[sub_col_c].value_counts().head(10)}")

    os_t_c = os_e_c = None
    for c in clin_c.columns:
        cu = c.upper()
        if cu == "OS_MONTHS": os_t_c = c
        if cu == "OS_STATUS":  os_e_c = c
    log(f"  OS time={os_t_c}  OS event={os_e_c}")

    # ── C.4: Align ────────────────────────────────────────────
    log("")
    log("── C.4: ALIGN ──")

    common_c = expr_t.index.intersection(clin_c.index)
    log(f"  Direct overlap: {len(common_c)}")

    if len(common_c) < 100:
        ei = pd.Index([str(x)[:15] for x in expr_t.index])
        ci = pd.Index([str(x)[:15] for x in clin_c.index])
        common_c = ei.intersection(ci)
        log(f"  15-char overlap: {len(common_c)}")
        if len(common_c) >= 100:
            expr_t  = expr_t.copy()
            expr_t.index  = ei
            clin_c  = clin_c.copy()
            clin_c.index  = ci

    if len(common_c) < 50:
        log("  Insufficient overlap — NOT TESTABLE")
        return {"status": "NOT_TESTABLE", "predictions": {}}

    gene_cols = [c for c in ["FOXA1","EZH2"]
                 if c in expr_t.columns]
    df_c = expr_t.loc[common_c, gene_cols].copy()

    extra = [sub_col_c]
    if os_t_c: extra.append(os_t_c)
    if os_e_c: extra.append(os_e_c)
    df_c = df_c.join(clin_c.loc[common_c, extra])

    rn = {sub_col_c: "subtype"}
    if os_t_c: rn[os_t_c] = "OS_MONTHS"
    if os_e_c: rn[os_e_c] = "OS_STATUS"
    df_c = df_c.rename(columns=rn)
    df_c["subtype_n"] = df_c["subtype"].apply(norm_subtype)
    df_c["ratio_c"]   = (
        df_c["FOXA1"].astype(float)
        / (df_c["EZH2"].astype(float).replace(0, np.nan)
           + 1e-6))

    log(f"  df_c: {df_c.shape}")
    log(f"  Subtype dist:\n"
        f"{df_c['subtype_n'].value_counts().head(8)}")

    # ── C.5: Ordering R3-A ───────────────────────────────────
    log("")
    log("── C.5: ORDERING TEST (R3-A) ──")

    all_s = ["LumA", "LumB", "HER2", "TNBC", "CL"]
    groups_c = {
        s: df_c[df_c["subtype_n"]==s]["ratio_c"].dropna().values
        for s in all_s
        if (df_c["subtype_n"]==s).sum() >= 5
    }
    medians_c = {s: float(np.median(v))
                 for s, v in groups_c.items()}

    log("  Ratio medians (mRNA):")
    for s in SUBTYPE_ORDER:
        if s in medians_c:
            log(f"    {s:6s}: n={len(groups_c[s]):4d}  "
                f"med={medians_c[s]:.4f}")

    main_s = {s: groups_c[s] for s in ["LumA","LumB",
                                         "HER2","TNBC"]
              if s in groups_c}
    kw_p_c, pw_c = kw_and_pairs(main_s)
    log(f"\n  KW p = {kw_p_c:.2e}"
        if kw_p_c is not None else "  KW p = None")

    adj_c = [("LumA","LumB"), ("LumB","HER2"), ("HER2","TNBC")]
    nc_c, nt_c, detail_c = order_check(medians_c, adj_c)
    for s1, s2, ok, m1, m2 in detail_c:
        log(f"  {s1}({m1:.4f}) > {s2}({m2:.4f}): "
            f"{'✓' if ok else '✗'}")

    r3a = (kw_p_c is not None and kw_p_c < 0.001
           and nc_c >= 2)
    preds["R3-A"] = "CONFIRMED" if r3a else "NOT CONFIRMED"
    log(f"  R3-A → {preds['R3-A']}")

    # ── C.6: LumA vs LumB R3-B ───────────────────────────────
    r3b_p = None
    if "LumA" in groups_c and "LumB" in groups_c:
        _, r3b_p = mannwhitneyu(
            groups_c["LumA"], groups_c["LumB"],
            alternative="two-sided")
        log(f"  LumA vs LumB p = {r3b_p:.2e}")
    r3b = (r3b_p is not None and r3b_p < 0.001)
    preds["R3-B"] = "CONFIRMED" if r3b else "NOT CONFIRMED"
    log(f"  R3-B → {preds['R3-B']}")

    # ── C.7: Survival R3-C ────────────────────────────────────
    log("")
    log("── C.7: SURVIVAL (R3-C) ──")
    r3c_p = None
    if ("OS_MONTHS" in df_c.columns
            and "OS_STATUS" in df_c.columns):
        df_c["OS_t"] = pd.to_numeric(
            df_c["OS_MONTHS"], errors="coerce")
        df_c["OS_e"] = df_c["OS_STATUS"].apply(
            lambda x: 1 if str(x) in
            ["1","1.0","DECEASED","1:DECEASED",
             "DEAD:DECEASED"]
            else 0)
        r3c_p = km_q14(df_c, "OS_t", "OS_e", "ratio_c",
                       title="METABRIC ratio vs OS")
        log(f"  KM p = {r3c_p:.4f}"
            if r3c_p is not None else "  KM: not run")
    preds["R3-C"] = (
        "CONFIRMED" if r3c_p is not None and r3c_p < 0.05
        else ("NOT TESTABLE" if r3c_p is None
              else "NOT CONFIRMED"))
    log(f"  R3-C → {preds['R3-C']}")

    # ── C.8: Figure ───────────────────────────────────────────
    n_p = 3 if (HAS_LIFELINES and r3c_p is not None) else 2
    fig_c, ax_c = plt.subplots(1, n_p,
                                figsize=(5*n_p + 1, 5))
    if n_p == 2:
        ax_c = list(ax_c)

    plot_s_c = [s for s in SUBTYPE_ORDER if s in groups_c]

    ax_c[0].violinplot(
        [groups_c[s] for s in plot_s_c],
        showmedians=True)
    ax_c[0].set_xticks(range(1, len(plot_s_c)+1))
    ax_c[0].set_xticklabels(plot_s_c, fontsize=8)
    ax_c[0].set_ylabel("FOXA1/EZH2 (z-score mRNA)")
    ax_c[0].set_title(
        f"METABRIC  n={len(df_c)}\n"
        f"KW p={'N/A' if kw_p_c is None else f'{kw_p_c:.2e}'}"
        f"  R3-A: {preds['R3-A']}")

    meds_c = [medians_c.get(s, np.nan) for s in plot_s_c]
    cols_c = [SUBTYPE_COLORS.get(s, "#888888")
              for s in plot_s_c]
    labs_c = [f"{s}\n(n={len(groups_c[s])})"
              for s in plot_s_c]
    ax_c[1].bar(labs_c, meds_c, color=cols_c, alpha=0.85,
                edgecolor="black", linewidth=0.7)
    ax_c[1].set_ylabel("Median FOXA1/EZH2 (mRNA)")
    ax_c[1].set_title(
        f"{nc_c}/{nt_c} correct\n"
        f"LumA>LumB p={r3b_p:.2e}" if r3b_p else
        f"{nc_c}/{nt_c} correct")

    if n_p == 3:
        km_q14(df_c, "OS_t", "OS_e", "ratio_c",
               ax=ax_c[2],
               title="METABRIC ratio vs OS")

    plt.tight_layout()
    plt.savefig(FIG_COMP_C, dpi=150, bbox_inches="tight")
    plt.close()
    log(f"  Figure: {FIG_COMP_C}")

    return {
        "status":    "COMPLETE",
        "n_samples": len(df_c),
        "medians":   medians_c,
        "kw_p":      kw_p_c,
        "nc": nc_c, "nt": nt_c,
        "r3b_p":     r3b_p,
        "r3c_p":     r3c_p,
        "predictions": preds,
    }


# ============================================================
# COMBINED FIGURE
# ============================================================

def make_combined(res_a, res_b, res_c):
    log("")
    log("── COMBINED FIGURE ──")
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    panels = [
        ("A: scRNA-seq\n(GSE176078)",
         res_a.get("medians",{}),
         res_a.get("status","—"),
         res_a.get("nc","—"), res_a.get("nt","—")),
        ("B: RPPA protein\n(TCGA)",
         res_b.get("medians",{}),
         res_b.get("status","—"),
         res_b.get("nc","—"), res_b.get("nt","—")),
        ("C: mRNA microarray\n(METABRIC)",
         res_c.get("medians",{}),
         res_c.get("status","—"),
         res_c.get("nc","—"), res_c.get("nt","—")),
    ]

    for ax, (title, meds, status, nc, nt) in \
            zip(axes, panels):
        if not meds or status == "NOT_TESTABLE":
            ax.text(0.5, 0.5, "NOT TESTABLE\n" + title,
                    ha="center", va="center",
                    transform=ax.transAxes, fontsize=10)
            ax.set_title(title)
            continue

        subs  = [s for s in SUBTYPE_ORDER if s in meds]
        vals  = [meds[s] for s in subs]
        cols  = [SUBTYPE_COLORS.get(s,"#888") for s in subs]
        ax.bar(subs, vals, color=cols, alpha=0.85,
               edgecolor="black", linewidth=0.8)

        for i in range(len(subs)-1):
            s1, s2 = subs[i], subs[i+1]
            ok = meds[s1] > meds[s2]
            ax.annotate(
                "✓" if ok else "✗",
                xy=(i+0.5, max(vals)*0.95),
                ha="center", fontsize=14,
                color="green" if ok else "red")

        ax.set_title(
            f"{title}\n{nc}/{nt} adjacent pairs", fontsize=9)
        ax.set_ylabel("Median FOXA1/EZH2 ratio")

    plt.suptitle(
        "FOXA1/EZH2 — Three-Dataset Validation\n"
        "OrganismCore RATIO-S1b | 2026-03-05",
        fontsize=11)
    plt.tight_layout()
    plt.savefig(FIG_COMBINED, dpi=150, bbox_inches="tight")
    plt.close()
    log(f"  Combined figure: {FIG_COMBINED}")


# ============================================================
# SCORECARD
# ============================================================

def make_scorecard(res_a, res_b, res_c):
    log("")
    log("=" * 65)
    log("SCORECARD — RATIO-S1b")
    log("=" * 65)

    rows = []
    for pred_id, status in \
            res_a.get("predictions", {}).items():
        rows.append({
            "id":         pred_id,
            "component":  "A",
            "dataset":    "GSE176078 scRNA-seq",
            "n":          res_a.get("n_cells","N/A"),
            "technology": "scRNA-seq log1p",
            "status":     status,
        })
    for pred_id, status in \
            res_b.get("predictions", {}).items():
        rows.append({
            "id":         pred_id,
            "component":  "B",
            "dataset":    "TCGA RPPA protein",
            "n":          res_b.get("n_samples","N/A"),
            "technology": "RPPA protein",
            "status":     status,
        })
    for pred_id, status in \
            res_c.get("predictions", {}).items():
        rows.append({
            "id":         pred_id,
            "component":  "C",
            "dataset":    "METABRIC mRNA",
            "n":          res_c.get("n_samples","N/A"),
            "technology": "mRNA microarray z-score",
            "status":     status,
        })

    df_sc = pd.DataFrame(rows)
    if not df_sc.empty:
        df_sc.to_csv(CSV_SCORECARD, index=False)
        log(f"  {'ID':<8}  {'Dataset':<22}  {'n':>6}  "
            f"Status")
        log(f"  {'─'*8}  {'─'*22}  {'─'*6}  {'─'*22}")
        for _, r in df_sc.iterrows():
            log(f"  {r['id']:<8}  {r['dataset']:<22}  "
                f"{str(r['n']):>6}  {r['status']}")

        n_conf = df_sc["status"].str.startswith(
            "CONFIRMED").sum()
        n_nc   = (df_sc["status"] == "NOT CONFIRMED").sum()
        n_nt   = df_sc["status"].str.contains(
            "NOT TESTABLE").sum()
        log(f"\n  CONFIRMED:    {n_conf}")
        log(f"  NOT CONFIRMED:{n_nc}")
        log(f"  NOT TESTABLE: {n_nt}")

    r2a = res_b.get("predictions",{}).get("R2-A","NOT RUN")
    log("")
    log("  ══════════════════════════════════════")
    log(f"  PRIMARY: R2-A = {r2a}")
    log("  ══════════════════════════════════════")
    log(f"\n  Scorecard: {CSV_SCORECARD}")
    return df_sc


# ============================================================
# MAIN
# ============================================================

def main():
    log("=" * 65)
    log("FOXA1/EZH2 RATIO — VALIDATION SCRIPT 1 v2 (RATIO-S1b)")
    log("OrganismCore | 2026-03-05 | FIXES: paths, RPPA, "
        "clinical index")
    log("=" * 65)
    log(f"Script dir:  {SCRIPT_DIR}")
    log(f"S1 data:     {S1_DATA}")
    log(f"S3 data:     {S3_DATA}")
    log(f"Output:      {BASE_DIR}")
    log(f"lifelines:   {HAS_LIFELINES}")
    log(f"sklearn:     {HAS_SKLEARN}")

    log("")
    log("── CACHE FILE STATUS ──")
    for label, path in [
        ("TCGA expr (S1)",     TCGA_EXPR_CACHE),
        ("TCGA clinical (S1)", TCGA_CLIN_FILE),
        ("scRNA cache (S1)",   SC_EXPR_CACHE),
        ("METABRIC expr (S3)", META_EXPR_FILE),
        ("METABRIC clin (S3)", META_CLIN_FILE),
        ("GSE25066 expr (S3)", GSE25066_EXPR),
        ("GSE25066 clin (S3)", GSE25066_CLIN),
    ]:
        exists = os.path.exists(path)
        size   = (f"{os.path.getsize(path):,} B"
                  if exists else "—")
        log(f"  {'✓' if exists else '✗'}  {label}: {size}")

    flush_log()

    log("")
    log("Starting Component A (scRNA-seq)...")
    flush_log()
    res_a = run_component_a()

    log("")
    log("Starting Component B (RPPA protein)...")
    flush_log()
    res_b = run_component_b()

    log("")
    log("Starting Component C (METABRIC mRNA)...")
    flush_log()
    res_c = run_component_c()

    make_combined(res_a, res_b, res_c)
    make_scorecard(res_a, res_b, res_c)

    log("")
    log("=" * 65)
    log("RATIO-S1b COMPLETE")
    log("=" * 65)
    log(f"  Component A status: {res_a['status']}")
    log(f"  Component B status: {res_b['status']}")
    log(f"  Component C status: {res_c['status']}")
    log(f"\n  Log:       {LOG_FILE}")
    log(f"  Scorecard: {CSV_SCORECARD}")
    log(f"  Figures:   {RESULTS_DIR}")
    log("")
    log("  Next: RATIO-S1c — "
        "ratio_script1_results_and_reasoning.md")
    log("=" * 65)
    flush_log()


if __name__ == "__main__":
    main()
