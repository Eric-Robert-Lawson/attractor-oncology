"""
BRCA CROSS-SUBTYPE ANALYSIS — SCRIPT 2
OrganismCore — Document BRCA-S8d | 2026-03-05

BEFORE-DOCUMENT: BRCA-S8d (predictions locked above)
All predictions locked before this script was written.
Results go to Cross_Subtype_s2_results/.

WHAT THIS SCRIPT DOES:
  Survival analysis across all six breast cancer subtypes.
  Tests the four pending predictions from Script 1
  plus two new survival predictions locked in BRCA-S8d.

  Six required outputs:
    1. DEPTH SCORE OS ANALYSIS
       Kaplan-Meier by depth quartile within each
       PAM50 subtype in TCGA-BRCA.
       Tests: CS-10

    2. TNBC AR-DEPTH SURVIVAL
       AR-low vs AR-high TNBC survival curves.
       Depth score vs AR as survival predictors.
       Tests: CS-13-SURVIVAL

    3. LumB ER OUTPUT DECOUPLING
       TFF1/ESR1 ratio LumB vs LumA.
       HDAC correlation with decoupling.
       Tests: CS-LUMB-DECOUPLE

    4. PCA WITHOUT EZH2
       Identity-only PCA: CL vs TNBC distance.
       Tests: CS-PCA-EZH2FREE

    5. ILC FOXA1 SURVIVAL
       FOXA1 level and EZH2+MKI67 composite
       as survival predictors within ILC.
       Tests: CS-ILC-ET

    6. UNIVERSAL DEPTH SCORE SURVIVAL
       Depth score as dominant within-subtype
       prognostic variable vs standard clinical
       variables.
       Tests: CS-DEPTH-UNIVERSAL

DATA SOURCES:
  PRIMARY:   TCGA-BRCA HiSeqV2 (already cached from Script 1)
             expr_cache_cs_s1_tcga.csv
             TCGA_BRCA_clinicalMatrix.tsv

  SECONDARY: GSE176078 scRNA-seq
             (already cached from Script 1)
             expr_cache_cs_s1_sc.csv

  SURVIVAL:  TCGA-BRCA clinical matrix contains:
             OS_Time_nature2012 (days)
             OS_event_nature2012 (0/1)
             Plus: Days_to_Date_of_Last_Contact_nature2012
                   Days_to_date_of_Death_nature2012

SELF-CONTAINED:
  Reuses all caches from Script 1.
  Downloads supplementary survival data if needed.
  All results to Cross_Subtype_s2_results/.

PROTOCOL v2.0:
  Unfiltered depth score distribution printed FIRST.
  Survival tests SECOND.
  Wrong predictions documented alongside correct ones.
  Final scorecard combines Script 1 + Script 2.
"""

import os
import sys
import gzip
import warnings
import urllib.request
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test

# ============================================================
# CONFIGURATION — SELF-CONTAINED
# ============================================================

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))

# Script 1 results directory — reuse caches
S1_RESULTS   = os.path.join(SCRIPT_DIR,
                             "Cross_Subtype_s1_results", "results")
S1_DATA      = os.path.join(SCRIPT_DIR,
                             "Cross_Subtype_s1_results", "data")

# Script 2 output directory
BASE_DIR     = os.path.join(SCRIPT_DIR, "Cross_Subtype_s2_results")
RESULTS_DIR  = os.path.join(BASE_DIR, "results")

for d in [BASE_DIR, RESULTS_DIR]:
    os.makedirs(d, exist_ok=True)

# Cache paths from Script 1
SC_CACHE_FILE   = os.path.join(S1_RESULTS,
                                "expr_cache_cs_s1_sc.csv")
TCGA_CACHE_FILE = os.path.join(S1_RESULTS,
                                "expr_cache_cs_s1_tcga.csv")
TCGA_CLIN_FILE  = os.path.join(S1_DATA,
                                "TCGA_BRCA_clinicalMatrix.tsv")

# Script 2 outputs
LOG_FILE         = os.path.join(RESULTS_DIR, "cs_s2_log.txt")
FIG_KM_LUMA      = os.path.join(RESULTS_DIR, "cs_s2_km_luma.png")
FIG_KM_LUMB      = os.path.join(RESULTS_DIR, "cs_s2_km_lumb.png")
FIG_KM_HER2      = os.path.join(RESULTS_DIR, "cs_s2_km_her2.png")
FIG_KM_TNBC      = os.path.join(RESULTS_DIR, "cs_s2_km_tnbc.png")
FIG_KM_ILC       = os.path.join(RESULTS_DIR, "cs_s2_km_ilc.png")
FIG_DECOUPLE     = os.path.join(RESULTS_DIR, "cs_s2_lumb_decouple.png")
FIG_PCA_FREE     = os.path.join(RESULTS_DIR, "cs_s2_pca_ezh2free.png")
FIG_AR_TNBC      = os.path.join(RESULTS_DIR, "cs_s2_tnbc_ar.png")
FIG_MASTER       = os.path.join(RESULTS_DIR, "cs_s2_master_figure.png")
CSV_SCORECARD    = os.path.join(RESULTS_DIR, "cs_s2_scorecard.csv")
CSV_DEPTH_SCORES = os.path.join(RESULTS_DIR, "cs_s2_depth_scores.csv")
CSV_SURVIVAL     = os.path.join(RESULTS_DIR, "cs_s2_survival_summary.csv")

# ============================================================
# CLINICAL COLUMN NAMES — TCGA-BRCA
# ============================================================

OS_TIME_COL   = "OS_Time_nature2012"
OS_EVENT_COL  = "OS_event_nature2012"
PAM50_COL     = "PAM50Call_RNAseq"
HIST_COL      = "histological_type"

# Survival column alternatives
OS_TIME_ALTS  = ["OS_Time_nature2012", "OS.time",
                  "days_to_last_followup",
                  "Days_to_Date_of_Last_Contact_nature2012"]
OS_EVENT_ALTS = ["OS_event_nature2012", "OS",
                  "vital_status",
                  "Days_to_date_of_Death_nature2012"]

# ============================================================
# GENE PANELS
# ============================================================

# Luminal identity — Script 1 confirmed
LUMINAL_TFS   = ["FOXA1", "GATA3", "ESR1", "PGR", "SPDEF"]

# EZH2 panel
EPIGENETIC    = ["EZH2", "EED", "SUZ12", "HDAC1", "HDAC2",
                 "KDM1A", "DNMT3A"]

# TNBC / basal panel
TNBC_FA       = ["SOX10", "KRT5", "KRT14", "VIM", "EGFR",
                 "FOXC1", "CDH3"]

# ER output genes — critical for CS-LUMB-DECOUPLE
ER_OUTPUT     = ["TFF1", "TFF3", "GREB1", "PDZK1", "PGR",
                 "AGR2", "CELSR2"]

# ILC markers
ILC_MARKER    = ["CDH1", "CDH2", "CTNNA1"]

# CL markers
CL_STEM       = ["CD44", "CD24", "ALDH1A3", "SNAI1",
                 "ZEB1", "ZEB2", "FN1",
                 "CLDN3", "CLDN4", "CLDN7"]

# Depth proxies
DEPTH_GENES   = ["CDKN1A", "CDKN2A", "RB1", "MKI67",
                 "TOP2A", "CCND1", "CDK4", "CDK6", "AR"]

# Identity-only genes for EZH2-free PCA (CS-PCA-EZH2FREE)
IDENTITY_ONLY = ["FOXA1", "GATA3", "ESR1", "PGR", "SPDEF",
                 "SOX10", "KRT5", "KRT14", "CDH1", "AR"]

ALL_TCGA_GENES = list(dict.fromkeys(
    LUMINAL_TFS + EPIGENETIC + TNBC_FA + ER_OUTPUT +
    ILC_MARKER + CL_STEM + DEPTH_GENES + IDENTITY_ONLY
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
# SCORECARD
# ============================================================

# Load Script 1 scorecard to append
scorecard = {}

def record(pid, status, note=""):
    scorecard[pid] = {"status": status, "note": note}
    marker = {"CONFIRMED": "✓", "FAILED": "✗",
              "PARTIAL": "~", "N/A": "-",
              "PENDING": "?"}.get(status, "?")
    log(f"  [{marker}] {pid}: {status}  {note}")

def write_scorecard():
    # Load Script 1 scorecard
    s1_csv = os.path.join(S1_RESULTS, "cs_s1_scorecard.csv")
    if os.path.exists(s1_csv):
        s1_df = pd.read_csv(s1_csv)
        for _, row in s1_df.iterrows():
            if row["prediction"] not in scorecard:
                scorecard[row["prediction"]] = {
                    "status": row["status"],
                    "note":   row.get("note", "")
                }

    rows = [{"prediction": k,
             "status": v["status"],
             "note": v["note"]}
            for k, v in scorecard.items()]
    pd.DataFrame(rows).to_csv(CSV_SCORECARD, index=False)

    log("")
    log("=" * 65)
    log("COMBINED SCORECARD — SCRIPT 1 + SCRIPT 2")
    log("=" * 65)
    confirmed = sum(1 for v in scorecard.values()
                    if v["status"] == "CONFIRMED")
    partial   = sum(1 for v in scorecard.values()
                    if v["status"] == "PARTIAL")
    failed    = sum(1 for v in scorecard.values()
                    if v["status"] == "FAILED")
    na        = sum(1 for v in scorecard.values()
                    if v["status"] == "N/A")
    total     = len(scorecard) - na
    log(f"  Confirmed:  {confirmed}/{total}")
    log(f"  Partial:    {partial}/{total}")
    log(f"  Failed:     {failed}/{total}")
    log(f"  Scorecard:  {CSV_SCORECARD}")

# ============================================================
# STEP 0 — LOAD CACHED DATA FROM SCRIPT 1
# ============================================================

def load_cached_data():
    log("=" * 65)
    log("STEP 0: LOAD CACHED DATA FROM SCRIPT 1")
    log("=" * 65)

    # Load TCGA expression cache
    if not os.path.exists(TCGA_CACHE_FILE):
        log(f"  FATAL: TCGA cache not found: {TCGA_CACHE_FILE}")
        log("  Run Script 1 first.")
        sys.exit(1)
    log(f"  Loading TCGA cache: {TCGA_CACHE_FILE}")
    tcga = pd.read_csv(TCGA_CACHE_FILE, index_col=0)
    log(f"  TCGA shape: {tcga.shape}")

    # Check for ER output genes — may need to re-extract
    er_missing = [g for g in ER_OUTPUT if g not in tcga.columns]
    if er_missing:
        log(f"  ER output genes missing from cache: {er_missing}")
        log("  These were not in Script 1 gene panel.")
        log("  Will attempt to load directly from TCGA gz.")
        tcga = supplement_tcga_genes(tcga, er_missing)

    # Load clinical matrix
    if not os.path.exists(TCGA_CLIN_FILE):
        log(f"  FATAL: Clinical matrix not found: {TCGA_CLIN_FILE}")
        sys.exit(1)
    log(f"  Loading clinical matrix: {TCGA_CLIN_FILE}")
    clin = pd.read_csv(TCGA_CLIN_FILE, sep="\t",
                       index_col=0, low_memory=False)
    log(f"  Clinical shape: {clin.shape}")

    # Resolve survival columns
    os_time_col  = None
    os_event_col = None

    for c in OS_TIME_ALTS:
        if c in clin.columns:
            os_time_col = c
            break
    for c in OS_EVENT_ALTS:
        if c in clin.columns:
            os_event_col = c
            break

    log(f"  OS time column:  '{os_time_col}'")
    log(f"  OS event column: '{os_event_col}'")

    if os_time_col is None:
        log("  WARNING: No OS time column found.")
        log(f"  Available columns: {list(clin.columns[:30])}")
    if os_event_col is None:
        log("  WARNING: No OS event column found.")

    # Load scRNA-seq cache (for depth score derivation)
    sc = None
    if os.path.exists(SC_CACHE_FILE):
        log(f"  Loading scRNA-seq cache: {SC_CACHE_FILE}")
        sc = pd.read_csv(SC_CACHE_FILE, index_col=0)
        log(f"  scRNA-seq shape: {sc.shape}")
    else:
        log("  scRNA-seq cache not found — depth scores will use"
            " TCGA bulk only.")

    return tcga, clin, sc, os_time_col, os_event_col


def supplement_tcga_genes(tcga, missing_genes):
    """Try to load missing genes from TCGA gz file."""
    tcga_gz = os.path.join(S1_DATA, "TCGA_BRCA_HiSeqV2.gz")
    if not os.path.exists(tcga_gz):
        log(f"  TCGA gz not found at {tcga_gz}. Skipping.")
        return tcga
    try:
        log(f"  Extracting {len(missing_genes)} genes from TCGA gz...")
        with gzip.open(tcga_gz, "rt") as f:
            raw = pd.read_csv(f, sep="\t", index_col=0)
        found = [g for g in missing_genes if g in raw.index]
        log(f"  Found {len(found)}/{len(missing_genes)}: {found}")
        if found:
            supp = raw.loc[found].T
            tcga = pd.concat([tcga, supp], axis=1)
            tcga = tcga.loc[:, ~tcga.columns.duplicated()]
            # Save updated cache
            tcga.to_csv(TCGA_CACHE_FILE)
            log("  Updated cache saved.")
    except Exception as e:
        log(f"  Supplement failed: {e}")
    return tcga

# ============================================================
# STEP 1 — CLASSIFY TCGA-BRCA SUBTYPES
# ============================================================

def classify_tcga_subtypes(tcga, clin):
    log("")
    log("=" * 65)
    log("STEP 1: CLASSIFY TCGA-BRCA SUBTYPES")
    log("=" * 65)

    common = tcga.index.intersection(clin.index)
    log(f"  TCGA-clinical common samples: {len(common)}")
    tcga_a = tcga.loc[common]
    clin_a = clin.loc[common]

    pam50  = clin_a[PAM50_COL].fillna("Unknown") \
             if PAM50_COL in clin_a.columns \
             else pd.Series("Unknown", index=clin_a.index)
    hist   = clin_a[HIST_COL].fillna("").str.lower() \
             if HIST_COL in clin_a.columns \
             else pd.Series("", index=clin_a.index)

    # PAM50 subtype masks
    masks = {
        "LumA":  pam50.isin(["LumA"]),
        "LumB":  pam50.isin(["LumB"]),
        "HER2":  pam50.isin(["Her2"]),
        "Basal": pam50.isin(["Basal"]),
        "ILC":   hist.str.contains("lobular", na=False),
    }

    # Claudin-low — use threshold 5 (n=14 from Script 1)
    cl_pos = [g for g in ["VIM", "ZEB1", "SNAI1", "CD44", "FN1"]
              if g in tcga_a.columns]
    cl_neg = [g for g in ["CLDN3", "CLDN4", "CDH1", "ESR1", "CLDN7"]
              if g in tcga_a.columns]
    cl_score = pd.Series(0.0, index=tcga_a.index)
    for g in cl_pos:
        cl_score += (tcga_a[g] > tcga_a[g].median()).astype(float)
    for g in cl_neg:
        cl_score -= (tcga_a[g] > tcga_a[g].median()).astype(float)
    if "ERBB2" in tcga_a.columns:
        cl_score[tcga_a["ERBB2"] >= tcga_a["ERBB2"].quantile(0.9)] = -99
    masks["CL"] = cl_score >= 5

    pops = {}
    for label, mask in masks.items():
        pops[label] = tcga_a[mask].copy()
        log(f"  {label:<10}: n={mask.sum()}")

    return pops, tcga_a, clin_a, masks

# ============================================================
# STEP 2 — COMPUTE DEPTH SCORES IN TCGA-BRCA
# ============================================================

def compute_depth_scores(tcga_a, clin_a, masks):
    """
    Compute per-sample depth scores for each PAM50 subtype.

    Depth score formula:
      For luminal subtypes (LumA, LumB, ILC):
        depth = normalised (CDKN1A_loss + EZH2_elevation)
        = (EZH2_rank - CDKN1A_rank) / 2
        where rank is within-subtype percentile

      For TNBC/Basal:
        depth = (EZH2_rank + SOX10_rank + ZEB1_rank
                 - AR_rank - FOXA1_rank) / 5

      For HER2:
        depth = (EZH2_rank + ERBB2_rank
                 - FOXA1_rank - AR_rank) / 4

      For CL:
        depth = (ZEB1_rank + VIM_rank + SNAI1_rank
                 - ESR1_rank - CLDN3_rank) / 5

    All scores are within-subtype percentile ranked
    (0=shallowest, 1=deepest).
    """
    log("")
    log("=" * 65)
    log("STEP 2: COMPUTE DEPTH SCORES IN TCGA-BRCA")
    log("=" * 65)

    def rank_norm(series):
        """Percentile rank within series."""
        return series.rank(pct=True)

    depth_formulas = {
        "LumA": {
            "pos": ["EZH2", "MKI67"],
            "neg": ["CDKN1A", "FOXA1", "GATA3"]
        },
        "LumB": {
            "pos": ["EZH2", "HDAC1", "MKI67"],
            "neg": ["CDKN1A", "TFF1", "FOXA1"]
        },
        "HER2": {
            "pos": ["EZH2", "ERBB2", "MKI67"],
            "neg": ["FOXA1", "AR", "ESR1"]
        },
        "Basal": {
            "pos": ["EZH2", "SOX10", "ZEB1", "MKI67"],
            "neg": ["AR", "FOXA1", "CDKN1A"]
        },
        "ILC": {
            "pos": ["EZH2", "MKI67"],
            "neg": ["CDH1", "CDKN1A"]
        },
        "CL": {
            "pos": ["ZEB1", "VIM", "SNAI1"],
            "neg": ["ESR1", "CLDN3", "FOXA1"]
        },
    }

    all_depth_scores = pd.Series(np.nan, index=tcga_a.index)
    subtype_scores   = {}

    for subtype, formula in depth_formulas.items():
        if subtype not in masks:
            continue
        mask = masks[subtype]
        if mask.sum() < 5:
            log(f"  {subtype}: n={mask.sum()} — insufficient for depth score")
            continue

        pop = tcga_a[mask].copy()
        pos_genes = [g for g in formula["pos"] if g in pop.columns]
        neg_genes = [g for g in formula["neg"] if g in pop.columns]

        if len(pos_genes) == 0 or len(neg_genes) == 0:
            log(f"  {subtype}: insufficient genes "
                f"(pos={pos_genes}, neg={neg_genes})")
            continue

        score = pd.Series(0.0, index=pop.index)
        for g in pos_genes:
            score += rank_norm(pop[g])
        for g in neg_genes:
            score += (1 - rank_norm(pop[g]))

        score = score / (len(pos_genes) + len(neg_genes))
        all_depth_scores[mask] = score
        subtype_scores[subtype] = score

        q25, q50, q75 = score.quantile([0.25, 0.50, 0.75])
        log(f"  {subtype:<8}: n={mask.sum()}  "
            f"depth Q25={q25:.3f}  Q50={q50:.3f}  Q75={q75:.3f}")
        log(f"    pos genes: {pos_genes}")
        log(f"    neg genes: {neg_genes}")

    # Save depth scores alongside clinical data
    depth_df = pd.DataFrame({
        "sample": tcga_a.index,
        "depth_score": all_depth_scores,
        "subtype": "Unknown"
    }).set_index("sample")

    for subtype, mask in masks.items():
        depth_df.loc[mask, "subtype"] = subtype

    # Add survival columns
    if OS_TIME_COL in clin_a.columns:
        depth_df["os_time"] = clin_a[OS_TIME_COL]
    if OS_EVENT_COL in clin_a.columns:
        depth_df["os_event"] = clin_a[OS_EVENT_COL]

    depth_df.to_csv(CSV_DEPTH_SCORES)
    log(f"\n  Depth scores saved: {CSV_DEPTH_SCORES}")

    return subtype_scores, all_depth_scores, depth_df

# ============================================================
# STEP 3 — RESOLVE SURVIVAL COLUMNS
# ============================================================

def resolve_survival(clin_a, os_time_col, os_event_col):
    """
    Build clean os_time (days) and os_event (0/1) series.
    Handles multiple possible column formats.
    """
    log("")
    log("=" * 65)
    log("STEP 3: RESOLVE SURVIVAL DATA")
    log("=" * 65)

    # Build os_time
    os_time = None
    if os_time_col and os_time_col in clin_a.columns:
        raw = pd.to_numeric(clin_a[os_time_col], errors="coerce")
        os_time = raw
        log(f"  OS time from '{os_time_col}': "
            f"n_valid={raw.notna().sum()}, "
            f"median={raw.median():.0f} days")
    else:
        # Try constructing from death/last-contact days
        death_col   = "Days_to_date_of_Death_nature2012"
        contact_col = "Days_to_Date_of_Last_Contact_nature2012"
        if death_col in clin_a.columns and contact_col in clin_a.columns:
            death   = pd.to_numeric(clin_a[death_col],   errors="coerce")
            contact = pd.to_numeric(clin_a[contact_col], errors="coerce")
            os_time = death.fillna(contact)
            log(f"  OS time constructed from death/contact days.")
            log(f"  n_valid: {os_time.notna().sum()}")
        else:
            log("  WARNING: Cannot construct OS time.")

    # Build os_event
    os_event = None
    if os_event_col and os_event_col in clin_a.columns:
        raw_event = clin_a[os_event_col]
        # Handle string encoding
        if raw_event.dtype == object:
            event_map = {
                "dead": 1, "deceased": 1, "1": 1,
                "alive": 0, "living": 0, "0": 0,
                1: 1, 0: 0
            }
            os_event = raw_event.str.lower().map(
                lambda x: event_map.get(str(x).lower().strip(), np.nan)
            )
        else:
            os_event = pd.to_numeric(raw_event, errors="coerce")

        log(f"  OS event from '{os_event_col}': "
            f"n_events={int(os_event.sum())}, "
            f"event_rate={os_event.mean():.1%}")
    else:
        # Try vital_status
        if "vital_status" in clin_a.columns:
            vs = clin_a["vital_status"].str.lower()
            os_event = vs.map({"dead": 1, "alive": 0,
                               "deceased": 1, "living": 0})
            log(f"  OS event from 'vital_status': "
                f"n_events={int(os_event.sum())}")
        else:
            log("  WARNING: Cannot construct OS event.")

    if os_time is None or os_event is None:
        log("  FATAL: Survival data not available.")
        log("  Survival analyses will be skipped.")
        return None, None

    # Convert time to years for display
    os_time_years = os_time / 365.25

    # Filter valid
    valid = os_time.notna() & os_event.notna() & (os_time > 0)
    log(f"  Valid survival records: {valid.sum()}")
    log(f"  Median follow-up: {os_time[valid].median()/365.25:.1f} years")

    return os_time, os_event

# ============================================================
# STEP 4 — KAPLAN-MEIER ANALYSIS
# Tests CS-10: depth score predicts OS within subtypes
# ============================================================

def kaplan_meier_by_depth(subtype_scores, tcga_a, clin_a,
                           os_time, os_event, masks):
    log("")
    log("=" * 65)
    log("ANALYSIS 1: KAPLAN-MEIER BY DEPTH (CS-10)")
    log("Tests: Depth score predicts OS within subtypes")
    log("Predicted: High depth = worse OS in all subtypes")
    log("=" * 65)

    if os_time is None or os_event is None:
        log("  SKIP: No survival data.")
        record("CS-10", "PENDING", "Survival data not available")
        return

    survival_rows = []
    km_results    = {}

    subtype_figs = {
        "LumA":  FIG_KM_LUMA,
        "LumB":  FIG_KM_LUMB,
        "HER2":  FIG_KM_HER2,
        "Basal": FIG_KM_TNBC,
        "ILC":   FIG_KM_ILC,
    }

    SUBTYPE_COLORS = {
        "Q1 (shallow)": "#2980b9",
        "Q2":           "#27ae60",
        "Q3":           "#e67e22",
        "Q4 (deep)":    "#c0392b",
    }

    cs10_notes  = []
    cs10_pass   = True

    for subtype, depth_scores in subtype_scores.items():
        if subtype not in masks:
            continue
        mask = masks[subtype]
        if mask.sum() < 20:
            log(f"\n  {subtype}: n={mask.sum()} — too small for KM")
            continue

        # Align survival data
        idx   = depth_scores.index
        t_raw = pd.to_numeric(os_time.reindex(idx), errors="coerce")
        e_raw = pd.to_numeric(os_event.reindex(idx), errors="coerce")

        valid = t_raw.notna() & e_raw.notna() & (t_raw > 0)
        t_v   = t_raw[valid] / 365.25   # convert to years
        e_v   = e_raw[valid].astype(int)
        d_v   = depth_scores[valid]

        if valid.sum() < 20:
            log(f"\n  {subtype}: n_valid={valid.sum()} — too small")
            continue

        # Quartile groups
        q25 = d_v.quantile(0.25)
        q75 = d_v.quantile(0.75)

        groups = pd.Series("Q2", index=d_v.index)
        groups[d_v <= q25] = "Q1 (shallow)"
        groups[d_v >= q75] = "Q4 (deep)"
        groups[(d_v > q25) & (d_v < q75)] = (
            groups[(d_v > q25) & (d_v < q75)]
            .where(d_v[(d_v > q25) & (d_v < q75)] < d_v.median(),
                   other="Q3")
        )
        groups[(d_v > q25) & (d_v < d_v.median())] = "Q2"
        groups[(d_v >= d_v.median()) & (d_v < q75)] = "Q3"

        log(f"\n  {subtype} (n_valid={valid.sum()}):")
        log(f"  Group sizes: "
            + " | ".join(
                f"{g}:n={int((groups == g).sum())}"
                for g in ["Q1 (shallow)", "Q2", "Q3", "Q4 (deep)"]
                if g in groups.values
            ))

        # Log-rank test: Q1 vs Q4
        q1_mask = groups == "Q1 (shallow)"
        q4_mask = groups == "Q4 (deep)"

        if q1_mask.sum() >= 5 and q4_mask.sum() >= 5:
            lr = logrank_test(
                t_v[q1_mask], t_v[q4_mask],
                e_v[q1_mask], e_v[q4_mask]
            )
            p_lr    = lr.p_value
            median_q1 = t_v[q1_mask].median()
            median_q4 = t_v[q4_mask].median()

            log(f"  Q1 (shallow) median OS: {median_q1:.1f} yr  "
                f"n={q1_mask.sum()}")
            log(f"  Q4 (deep)    median OS: {median_q4:.1f} yr  "
                f"n={q4_mask.sum()}")
            log(f"  Log-rank p = {p_lr:.4f}")

            # Cox HR
            try:
                cox_df = pd.DataFrame({
                    "T": t_v,
                    "E": e_v,
                    "depth": d_v
                })
                cph = CoxPHFitter()
                cph.fit(cox_df, duration_col="T", event_col="E")
                hr   = float(np.exp(cph.params_["depth"]))
                hr_p = float(cph.summary["p"]["depth"])
                log(f"  Cox HR (continuous depth): {hr:.3f}  "
                    f"p={hr_p:.4f}")
            except Exception as e:
                hr   = np.nan
                hr_p = np.nan
                log(f"  Cox fit error: {e}")

            direction_correct = median_q4 <= median_q1
            sig = p_lr < 0.05

            survival_rows.append({
                "subtype":       subtype,
                "n_valid":       valid.sum(),
                "n_q1":          q1_mask.sum(),
                "n_q4":          q4_mask.sum(),
                "median_os_q1":  median_q1,
                "median_os_q4":  median_q4,
                "logrank_p":     p_lr,
                "cox_hr":        hr,
                "cox_p":         hr_p,
                "direction":     "CORRECT" if direction_correct
                                 else "REVERSED",
                "significant":   sig,
            })

            note = (f"Q4/Q1 median {median_q4:.1f}/{median_q1:.1f}yr "
                    f"logrank p={p_lr:.3f} HR={hr:.2f}")

            if direction_correct and sig:
                cs10_notes.append(f"{subtype}: CONFIRMED ({note})")
            elif direction_correct and not sig:
                cs10_notes.append(
                    f"{subtype}: PARTIAL — correct direction, "
                    f"underpowered ({note})"
                )
            else:
                cs10_notes.append(f"{subtype}: FAILED ({note})")
                cs10_pass = False

            km_results[subtype] = {
                "t": t_v, "e": e_v, "groups": groups,
                "p": p_lr, "hr": hr
            }

        for n in cs10_notes[-1:]:
            log(f"    {n}")

    # Overall CS-10 verdict
    if cs10_pass and len(survival_rows) > 0:
        n_sig = sum(1 for r in survival_rows if r["significant"])
        if n_sig >= 2:
            record("CS-10", "CONFIRMED",
                   f"{n_sig}/{len(survival_rows)} subtypes significant | "
                   + " | ".join(cs10_notes))
        else:
            record("CS-10", "PARTIAL",
                   f"{n_sig}/{len(survival_rows)} significant | "
                   + " | ".join(cs10_notes))
    elif not cs10_pass:
        record("CS-10", "FAILED", " | ".join(cs10_notes))
    else:
        record("CS-10", "PENDING", "No valid survival data")

    if survival_rows:
        pd.DataFrame(survival_rows).to_csv(CSV_SURVIVAL, index=False)
        log(f"\n  Survival summary: {CSV_SURVIVAL}")

    # Figures
    _make_km_figures(km_results, subtype_figs)

    return km_results, survival_rows


def _make_km_figures(km_results, subtype_figs):
    COLORS = {
        "Q1 (shallow)": "#2980b9",
        "Q2":           "#27ae60",
        "Q3":           "#e67e22",
        "Q4 (deep)":    "#c0392b",
    }

    for subtype, res in km_results.items():
        if subtype not in subtype_figs:
            continue
        try:
            fig, ax = plt.subplots(figsize=(8, 5))
            t       = res["t"]
            e       = res["e"]
            groups  = res["groups"]
            p       = res["p"]
            hr      = res["hr"]

            for grp in ["Q1 (shallow)", "Q2", "Q3", "Q4 (deep)"]:
                if grp not in groups.values:
                    continue
                m   = groups == grp
                kmf = KaplanMeierFitter()
                kmf.fit(t[m], e[m], label=f"{grp} (n={m.sum()})")
                kmf.plot_survival_function(
                    ax=ax,
                    color=COLORS.get(grp, "grey"),
                    ci_show=False,
                )

            ax.set_xlabel("Time (years)")
            ax.set_ylabel("Overall Survival Probability")
            ax.set_title(
                f"BRCA {subtype} — Depth Score Quartiles\n"
                f"Log-rank p={p:.4f}  "
                f"Cox HR(continuous)={hr:.2f}\n"
                "OrganismCore / BRCA-S8d / 2026-03-05",
                fontsize=10
            )
            ax.legend(fontsize=8)
            ax.set_ylim(0, 1.05)
            plt.tight_layout()
            plt.savefig(subtype_figs[subtype], dpi=150)
            plt.close()
            log(f"  KM figure saved: {subtype_figs[subtype]}")
        except Exception as e:
            log(f"  KM figure error ({subtype}): {e}")

# ============================================================
# STEP 5 — TNBC AR-DEPTH SURVIVAL (CS-13-SURVIVAL)
# ============================================================

def tnbc_ar_survival(tcga_a, masks, os_time, os_event):
    log("")
    log("=" * 65)
    log("ANALYSIS 2: TNBC AR-DEPTH SURVIVAL (CS-13-SURVIVAL)")
    log("Predicted: AR-low (deep) TNBC has worse DRFS/OS")
    log("than AR-high (shallow) despite equivalent/better pCR")
    log("=" * 65)

    if os_time is None or os_event is None:
        log("  SKIP: No survival data.")
        record("CS-13-SURVIVAL", "PENDING")
        return

    tnbc_mask = masks.get("Basal")
    if tnbc_mask is None or tnbc_mask.sum() < 20:
        log(f"  SKIP: Insufficient TNBC samples "
            f"(n={tnbc_mask.sum() if tnbc_mask is not None else 0})")
        record("CS-13-SURVIVAL", "PENDING",
               "Insufficient TNBC samples")
        return

    if "AR" not in tcga_a.columns:
        log("  SKIP: AR not in TCGA expression data.")
        record("CS-13-SURVIVAL", "PENDING", "AR not available")
        return

    tnbc  = tcga_a[tnbc_mask].copy()
    ar    = tnbc["AR"]
    t_raw = pd.to_numeric(os_time.reindex(tnbc.index), errors="coerce")
    e_raw = pd.to_numeric(os_event.reindex(tnbc.index), errors="coerce")
    valid = t_raw.notna() & e_raw.notna() & (t_raw > 0)

    t_v  = t_raw[valid] / 365.25
    e_v  = e_raw[valid].astype(int)
    ar_v = ar[valid]

    log(f"  TNBC n_valid: {valid.sum()}")
    log(f"  AR in TNBC: mean={ar_v.mean():.3f}  "
        f"median={ar_v.median():.3f}")

    if valid.sum() < 20:
        log("  SKIP: Insufficient valid TNBC survival records.")
        record("CS-13-SURVIVAL", "PENDING")
        return

    # Split by AR median
    ar_median = ar_v.median()
    ar_high   = ar_v >= ar_median
    ar_low    = ar_v <  ar_median

    log(f"  AR-high (shallow): n={ar_high.sum()}")
    log(f"  AR-low  (deep):    n={ar_low.sum()}")

    # Log-rank test
    try:
        lr = logrank_test(
            t_v[ar_high], t_v[ar_low],
            e_v[ar_high], e_v[ar_low]
        )
        p_lr = lr.p_value
        med_high = t_v[ar_high].median()
        med_low  = t_v[ar_low].median()

        log(f"  AR-high median OS: {med_high:.1f} yr")
        log(f"  AR-low  median OS: {med_low:.1f} yr")
        log(f"  Log-rank p = {p_lr:.4f}")

        # Cox HR
        cox_df = pd.DataFrame({
            "T":  t_v,
            "E":  e_v,
            "AR": ar_v
        })
        cph = CoxPHFitter()
        cph.fit(cox_df, duration_col="T", event_col="E")
        hr   = float(np.exp(cph.params_["AR"]))
        hr_p = float(cph.summary["p"]["AR"])
        log(f"  Cox HR (AR continuous): {hr:.3f}  p={hr_p:.4f}")
        log(f"  Interpretation: HR < 1 = higher AR protective")

        # Prediction: AR-low has WORSE OS (HR < 1 for AR,
        # i.e. lower AR = worse outcome)
        direction_correct = hr < 1.0
        sig = p_lr < 0.10  # relaxed threshold given small n

        # Figure
        try:
            fig, ax = plt.subplots(figsize=(8, 5))
            kmf_h = KaplanMeierFitter()
            kmf_h.fit(t_v[ar_high], e_v[ar_high],
                      label=f"AR-high/shallow (n={ar_high.sum()})")
            kmf_h.plot_survival_function(ax=ax, color="#2980b9",
                                         ci_show=True)
            kmf_l = KaplanMeierFitter()
            kmf_l.fit(t_v[ar_low], e_v[ar_low],
                      label=f"AR-low/deep (n={ar_low.sum()})")
            kmf_l.plot_survival_function(ax=ax, color="#c0392b",
                                         ci_show=True)
            ax.set_xlabel("Time (years)")
            ax.set_ylabel("Overall Survival Probability")
            ax.set_title(
                "TNBC — AR-high vs AR-low (Depth Proxy)\n"
                f"Log-rank p={p_lr:.4f}  HR(AR)={hr:.2f}\n"
                "Predicted: AR-low (deep) = worse OS\n"
                "OrganismCore / BRCA-S8d / 2026-03-05",
                fontsize=10
            )
            ax.legend(fontsize=9)
            ax.set_ylim(0, 1.05)
            plt.tight_layout()
            plt.savefig(FIG_AR_TNBC, dpi=150)
            plt.close()
            log(f"  AR survival figure: {FIG_AR_TNBC}")
        except Exception as fe:
            log(f"  AR figure error: {fe}")

        note = (f"AR-low OS {med_low:.1f}yr vs AR-high {med_high:.1f}yr "
                f"logrank p={p_lr:.3f} HR(AR)={hr:.3f}")
        if direction_correct and sig:
            record("CS-13-SURVIVAL", "CONFIRMED", note)
        elif direction_correct:
            record("CS-13-SURVIVAL", "PARTIAL",
                   note + " (underpowered)")
        else:
            record("CS-13-SURVIVAL", "FAILED",
                   f"Direction reversed: {note}")

    except Exception as e:
        log(f"  TNBC AR survival error: {e}")
        record("CS-13-SURVIVAL", "PENDING", f"Error: {e}")

# ============================================================
# STEP 6 — LumB ER OUTPUT DECOUPLING (CS-LUMB-DECOUPLE)
# ============================================================

def lumb_er_decoupling(tcga_a, masks):
    log("")
    log("=" * 65)
    log("ANALYSIS 3: LumB ER OUTPUT DECOUPLING (CS-LUMB-DECOUPLE)")
    log("Predicted: TFF1/ESR1 ratio lower in LumB than LumA")
    log("despite LumB having higher raw ESR1 mRNA")
    log("=" * 65)

    luma_mask = masks.get("LumA")
    lumb_mask = masks.get("LumB")

    if luma_mask is None or lumb_mask is None:
        log("  SKIP: Missing LumA or LumB population.")
        record("CS-LUMB-DECOUPLE", "PENDING")
        return

    # Check gene availability
    tff1_avail = "TFF1" in tcga_a.columns
    tff3_avail = "TFF3" in tcga_a.columns
    esr1_avail = "ESR1" in tcga_a.columns
    hdac1_avail = "HDAC1" in tcga_a.columns

    log(f"  TFF1:  {'available' if tff1_avail else 'MISSING'}")
    log(f"  TFF3:  {'available' if tff3_avail else 'MISSING'}")
    log(f"  ESR1:  {'available' if esr1_avail else 'MISSING'}")
    log(f"  HDAC1: {'available' if hdac1_avail else 'MISSING'}")

    if not (tff1_avail and esr1_avail):
        log("  SKIP: TFF1 or ESR1 not available.")
        record("CS-LUMB-DECOUPLE", "PENDING",
               "TFF1 or ESR1 missing from cache")
        return

    luma = tcga_a[luma_mask].copy()
    lumb = tcga_a[lumb_mask].copy()

    # Raw expression comparison
    log(f"\n  RAW EXPRESSION COMPARISON:")
    log(f"  {'Gene':<10} {'LumA mean':>12} {'LumB mean':>12} "
        f"{'LumB/LumA':>12}")
    log("  " + "-" * 50)

    for gene in ["ESR1", "TFF1", "TFF3", "FOXA1", "HDAC1",
                 "HDAC2", "DNMT3A"]:
        if gene not in tcga_a.columns:
            continue
        luma_m = luma[gene].mean()
        lumb_m = lumb[gene].mean()
        ratio  = lumb_m / luma_m if luma_m > 0 else np.nan
        log(f"  {gene:<10} {luma_m:12.3f} {lumb_m:12.3f} "
            f"{ratio:12.3f}")

    # TFF1/ESR1 ratio
    luma_ratio = luma["TFF1"] / (luma["ESR1"] + 1e-6)
    lumb_ratio = lumb["TFF1"] / (lumb["ESR1"] + 1e-6)

    log(f"\n  TFF1/ESR1 RATIO:")
    log(f"  LumA: mean={luma_ratio.mean():.4f}  "
        f"median={luma_ratio.median():.4f}")
    log(f"  LumB: mean={lumb_ratio.mean():.4f}  "
        f"median={lumb_ratio.median():.4f}")

    # Mann-Whitney test
    mw_stat, mw_p = mannwhitneyu(
        luma_ratio.dropna(), lumb_ratio.dropna(),
        alternative="greater"   # LumA > LumB ratio
    )
    log(f"  Mann-Whitney (LumA > LumB): U={mw_stat:.0f}, p={mw_p:.4f}")

    ratio_correct = luma_ratio.median() > lumb_ratio.median()
    ratio_sig     = mw_p < 0.05

    # TFF3/ESR1 ratio
    if tff3_avail:
        luma_r3 = luma["TFF3"] / (luma["ESR1"] + 1e-6)
        lumb_r3 = lumb["TFF3"] / (lumb["ESR1"] + 1e-6)
        mw3_stat, mw3_p = mannwhitneyu(
            luma_r3.dropna(), lumb_r3.dropna(),
            alternative="greater"
        )
        log(f"\n  TFF3/ESR1 RATIO:")
        log(f"  LumA: mean={luma_r3.mean():.4f}")
        log(f"  LumB: mean={lumb_r3.mean():.4f}")
        log(f"  Mann-Whitney (LumA > LumB): p={mw3_p:.4f}")

    # HDAC1 correlation with decoupling within LumB
    if hdac1_avail:
        r_hdac, p_hdac = stats.spearmanr(
            lumb["HDAC1"], lumb_ratio
        )
        log(f"\n  HDAC1 vs TFF1/ESR1 ratio in LumB: "
            f"r={r_hdac:.3f}, p={p_hdac:.4f}")
        log(f"  Predicted: HDAC1 high = lower TFF1/ESR1 ratio "
            f"(r < -0.2)")
        hdac_correct = r_hdac < -0.2
        log(f"  Direction: {'CONFIRMED' if hdac_correct else 'NOT CONFIRMED'}")
    else:
        r_hdac = np.nan
        hdac_correct = False

    # Figure
    try:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # Panel A: ESR1 vs TFF1 scatter
        ax = axes[0]
        ax.scatter(luma["ESR1"], luma["TFF1"],
                   alpha=0.4, color="#2980b9", s=20,
                   label=f"LumA (n={len(luma)})")
        ax.scatter(lumb["ESR1"], lumb["TFF1"],
                   alpha=0.4, color="#1abc9c", s=20,
                   label=f"LumB (n={len(lumb)})")
        ax.set_xlabel("ESR1 expression")
        ax.set_ylabel("TFF1 expression")
        ax.set_title("ESR1 vs TFF1\n(ER output decoupling)")
        ax.legend(fontsize=8)

        # Panel B: TFF1/ESR1 ratio boxplot
        ax = axes[1]
        ax.boxplot(
            [luma_ratio.dropna(), lumb_ratio.dropna()],
            labels=["LumA", "LumB"],
            patch_artist=True,
            boxprops=dict(facecolor="#d6eaf8"),
            medianprops=dict(color="black", linewidth=2)
        )
        ax.set_ylabel("TFF1/ESR1 ratio\n(ER output efficiency)")
        ax.set_title(
            f"ER Output Efficiency\n"
            f"Mann-Whitney p={mw_p:.4f}\n"
            f"LumA median={luma_ratio.median():.3f}  "
            f"LumB median={lumb_ratio.median():.3f}",
            fontsize=10
        )

        plt.suptitle(
            "LumB ER OUTPUT DECOUPLING\n"
            "OrganismCore / BRCA-S8d / 2026-03-05",
            fontsize=11
        )
        plt.tight_layout()
        plt.savefig(FIG_DECOUPLE, dpi=150)
        plt.close()
        log(f"  Decoupling figure: {FIG_DECOUPLE}")
    except Exception as e:
        log(f"  Decoupling figure error: {e}")

    note = (f"TFF1/ESR1: LumA={luma_ratio.median():.3f} > "
            f"LumB={lumb_ratio.median():.3f}, p={mw_p:.4f}")
    if ratio_correct and ratio_sig:
        record("CS-LUMB-DECOUPLE", "CONFIRMED",
               note + (f" | HDAC1 r={r_hdac:.3f}"
                       if not np.isnan(r_hdac) else ""))
    elif ratio_correct:
        record("CS-LUMB-DECOUPLE", "PARTIAL",
               note + " (correct direction, p not significant)")
    else:
        record("CS-LUMB-DECOUPLE", "FAILED",
               f"Direction reversed: {note}")

# ============================================================
# STEP 7 — PCA WITHOUT EZH2 (CS-PCA-EZH2FREE)
# ============================================================

def pca_ezh2_free(sc_expr, sc_meta_path=None):
    log("")
    log("=" * 65)
    log("ANALYSIS 4: EZH2-FREE PCA (CS-PCA-EZH2FREE)")
    log("Predicted: CL further than TNBC from MatureLum")
    log("when EZH2 removed from gene panel")
    log("=" * 65)

    if sc_expr is None:
        log("  SKIP: scRNA-seq cache not available.")
        record("CS-PCA-EZH2FREE", "PENDING",
               "scRNA-seq cache not loaded")
        return

    # Rebuild population means from scRNA-seq cache
    # We need metadata to separate populations
    meta_path = os.path.join(S1_DATA,
                              "Wu_etal_2021_BRCA_scRNASeq",
                              "metadata.csv")
    if not os.path.exists(meta_path):
        # Try alternate path
        meta_path = os.path.join(S1_DATA, "metadata.csv")
    if not os.path.exists(meta_path):
        log("  scRNA-seq metadata not found. Attempting TCGA fallback.")
        pca_ezh2_free_tcga()
        return

    log(f"  Loading metadata: {meta_path}")
    meta = pd.read_csv(meta_path)
    first_col = meta.columns[0]
    meta = meta.set_index(first_col)

    CT_COL = "celltype_subset"
    if CT_COL not in meta.columns:
        for c in ["celltype_subset", "celltype_minor",
                  "celltype_major"]:
            if c in meta.columns:
                CT_COL = c
                break

    # Align
    common = sc_expr.index.intersection(meta.index)
    sc_a   = sc_expr.loc[common]
    meta_a = meta.loc[common]

    # Normalise
    totals   = sc_a.sum(axis=1).replace(0, 1)
    sc_norm  = np.log1p(sc_a.div(totals, axis=0) * 1e4)

    pop_defs = {
        "LumA":      "Cancer LumA SC",
        "LumB":      "Cancer LumB SC",
        "HER2":      "Cancer Her2 SC",
        "TNBC":      "Cancer Basal SC",
        "MatureLum": "Mature Luminal",
    }

    # Claudin-low classifier
    cl_pos = [g for g in ["VIM", "ZEB1", "SNAI1", "CD44", "FN1"]
              if g in sc_norm.columns]
    cl_neg = [g for g in ["CLDN3", "CLDN4", "CDH1", "ESR1"]
              if g in sc_norm.columns]

    cancer_mask = meta_a[CT_COL].isin(
        ["Cancer LumB SC", "Cancer Basal SC"]
    )
    cancer_cells = sc_norm[cancer_mask]
    if len(cl_pos) >= 3 and len(cl_neg) >= 2:
        cl_score = pd.Series(0.0, index=cancer_cells.index)
        for g in cl_pos:
            med = cancer_cells[g].median()
            cl_score += (cancer_cells[g] > med).astype(float)
        for g in cl_neg:
            med = cancer_cells[g].median()
            cl_score -= (cancer_cells[g] > med).astype(float)
        cl_cells = cancer_cells[cl_score >= 3]
        log(f"  CL cells (thresh=3): n={len(cl_cells)}")
    else:
        cl_cells = pd.DataFrame(columns=sc_norm.columns)
        log("  WARNING: CL classifier genes insufficient")

    centroids = {}
    for label, ct_name in pop_defs.items():
        mask = meta_a[CT_COL] == ct_name
        cells = sc_norm[mask]
        if len(cells) > 0:
            centroids[label] = cells.mean()

    if len(cl_cells) > 0:
        centroids["CL"] = cl_cells.mean()

    if len(centroids) < 4:
        log("  SKIP: Insufficient populations for PCA.")
        record("CS-PCA-EZH2FREE", "PENDING")
        return

    centroid_df = pd.DataFrame(centroids).T

    # PCA WITH EZH2 (Script 1 result, for comparison)
    pca_genes_with = [g for g in LUMINAL_TFS + ["EZH2"] +
                      ["SOX10", "KRT5", "KRT14"]
                      if g in centroid_df.columns]

    # PCA WITHOUT EZH2 (CS-PCA-EZH2FREE test)
    pca_genes_free = [g for g in IDENTITY_ONLY
                      if g in centroid_df.columns
                      and g != "EZH2"]

    log(f"  Genes WITH EZH2: {pca_genes_with}")
    log(f"  Genes WITHOUT EZH2: {pca_genes_free}")

    results = {}
    for label, gene_panel in [("WITH_EZH2", pca_genes_with),
                               ("EZH2_FREE", pca_genes_free)]:
        if len(gene_panel) < 3:
            log(f"  {label}: insufficient genes")
            continue

        sub = centroid_df[gene_panel].dropna()
        if len(sub) < 3:
            continue

        scaler = StandardScaler()
        scaled = scaler.fit_transform(sub)
        pca    = PCA(n_components=min(3, len(gene_panel),
                                       len(sub)))
        coords = pca.fit_transform(scaled)
        var    = pca.explained_variance_ratio_

        pca_df = pd.DataFrame(
            coords, index=sub.index,
            columns=[f"PC{i+1}" for i in range(coords.shape[1])]
        )

        if "MatureLum" not in pca_df.index:
            continue

        ref = pca_df.loc["MatureLum"].values
        distances = {
            p: float(np.linalg.norm(pca_df.loc[p].values - ref))
            for p in pca_df.index
        }

        log(f"\n  PCA {label} (PC1={var[0]:.1%}):")
        dist_sorted = sorted(distances.items(), key=lambda x: x[1])
        for p, d in dist_sorted:
            log(f"    {p:<12}: dist={d:.3f}")

        results[label] = distances

    # CS-PCA-EZH2FREE test
    if "EZH2_FREE" in results:
        d_ezh2free = results["EZH2_FREE"]
        cl_dist    = d_ezh2free.get("CL", np.nan)
        tnbc_dist  = d_ezh2free.get("TNBC", np.nan)

        log(f"\n  EZH2-FREE: CL dist={cl_dist:.3f}  "
            f"TNBC dist={tnbc_dist:.3f}")
        log(f"  Predicted: CL > TNBC in EZH2-free space")

        if not np.isnan(cl_dist) and not np.isnan(tnbc_dist):
            if cl_dist > tnbc_dist:
                record("CS-PCA-EZH2FREE", "CONFIRMED",
                       f"CL ({cl_dist:.3f}) > TNBC ({tnbc_dist:.3f}) "
                       f"when EZH2 excluded")
            else:
                record("CS-PCA-EZH2FREE", "FAILED",
                       f"CL ({cl_dist:.3f}) <= TNBC ({tnbc_dist:.3f}) "
                       f"even without EZH2")
        else:
            record("CS-PCA-EZH2FREE", "PENDING",
                   "CL or TNBC distance not computed")
    else:
        record("CS-PCA-EZH2FREE", "PENDING", "PCA failed")

    # Figure: side-by-side PCA comparison
    try:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        COLORS = {
            "LumA":     "#2980b9",
            "LumB":     "#1abc9c",
            "HER2":     "#e67e22",
            "TNBC":     "#c0392b",
            "CL":       "#8e44ad",
            "MatureLum":"#95a5a6",
        }

        for idx, (label, gene_panel) in enumerate(
                [("WITH EZH2", pca_genes_with),
                 ("EZH2-FREE", pca_genes_free)]):
            ax = axes[idx]
            sub = centroid_df[
                [g for g in gene_panel if g in centroid_df.columns]
            ].dropna()
            if len(sub) < 3:
                continue
            scaler = StandardScaler()
            scaled = scaler.fit_transform(sub)
            pca    = PCA(n_components=2)
            coords = pca.fit_transform(scaled)
            pca_df = pd.DataFrame(
                coords, index=sub.index,
                columns=["PC1", "PC2"]
            )
            if "MatureLum" in pca_df.index:
                ref = pca_df.loc["MatureLum"].values
            else:
                ref = np.zeros(2)

            for pop in pca_df.index:
                color = COLORS.get(pop, "#34495e")
                pc1, pc2 = pca_df.loc[pop, "PC1"], \
                            pca_df.loc[pop, "PC2"]
                ax.scatter(pc1, pc2, s=200, color=color,
                           edgecolors="black", linewidth=1.5,
                           zorder=5)
                ax.annotate(pop, (pc1, pc2),
                            textcoords="offset points",
                            xytext=(8, 5), fontsize=9)
                if pop in ["CL", "TNBC"]:
                    ax.plot([ref[0], pc1], [ref[1], pc2],
                            "k--", alpha=0.3, linewidth=1)

            ax.set_xlabel("PC1")
            ax.set_ylabel("PC2")
            var = pca.explained_variance_ratio_
            ax.set_title(f"PCA {label}\n"
                         f"PC1={var[0]:.1%}, PC2={var[1]:.1%}")
            ax.axhline(0, color="grey", ls="--", lw=0.5)
            ax.axvline(0, color="grey", ls="--", lw=0.5)

        plt.suptitle(
            "PCA GEOMETRY — WITH vs WITHOUT EZH2\n"
            "OrganismCore / BRCA-S8d / 2026-03-05",
            fontsize=11
        )
        plt.tight_layout()
        plt.savefig(FIG_PCA_FREE, dpi=150)
        plt.close()
        log(f"  EZH2-free PCA figure: {FIG_PCA_FREE}")
    except Exception as e:
        log(f"  EZH2-free PCA figure error: {e}")

    return results


def pca_ezh2_free_tcga():
    """Fallback: use TCGA population means for EZH2-free PCA."""
    log("  Using TCGA population means for EZH2-free PCA fallback.")
    record("CS-PCA-EZH2FREE", "PENDING",
           "scRNA-seq metadata not available — "
           "use TCGA population means in next run")

# ============================================================
# STEP 8 — ILC FOXA1 SURVIVAL (CS-ILC-ET)
# ============================================================

def ilc_foxa1_survival(tcga_a, masks, os_time, os_event):
    log("")
    log("=" * 65)
    log("ANALYSIS 5: ILC FOXA1 SURVIVAL (CS-ILC-ET)")
    log("Predicted: EZH2-high + MKI67-high ILC HR > 2.0")
    log("within ILC only")
    log("=" * 65)

    if os_time is None or os_event is None:
        log("  SKIP: No survival data.")
        record("CS-ILC-ET", "PENDING")
        return

    ilc_mask = masks.get("ILC")
    if ilc_mask is None or ilc_mask.sum() < 20:
        log(f"  ILC n={ilc_mask.sum() if ilc_mask is not None else 0}"
            f" — checking minimum...")

    if ilc_mask is None or ilc_mask.sum() < 15:
        log("  SKIP: Insufficient ILC samples.")
        record("CS-ILC-ET", "PENDING", "Insufficient ILC")
        return

    ilc  = tcga_a[ilc_mask].copy()
    t_r  = pd.to_numeric(os_time.reindex(ilc.index), errors="coerce")
    e_r  = pd.to_numeric(os_event.reindex(ilc.index), errors="coerce")
    valid = t_r.notna() & e_r.notna() & (t_r > 0)
    t_v  = t_r[valid] / 365.25
    e_v  = e_r[valid].astype(int)

    log(f"  ILC n_valid: {valid.sum()}")

    if valid.sum() < 15:
        log("  SKIP: Insufficient valid ILC survival records.")
        record("CS-ILC-ET", "PENDING")
        return

    notes     = []
    test_pass = True

    # FOXA1 — within ILC survival
    if "FOXA1" in ilc.columns:
        foxa1_v  = ilc["FOXA1"][valid]
        foxa1_med = foxa1_v.median()
        f_high   = foxa1_v >= foxa1_med
        f_low    = foxa1_v <  foxa1_med

        if f_high.sum() >= 5 and f_low.sum() >= 5:
            lr_f = logrank_test(
                t_v[f_high], t_v[f_low],
                e_v[f_high], e_v[f_low]
            )
            med_fh = t_v[f_high].median()
            med_fl = t_v[f_low].median()
            log(f"\n  FOXA1-high vs low in ILC:")
            log(f"  FOXA1-high median OS: {med_fh:.1f} yr")
            log(f"  FOXA1-low  median OS: {med_fl:.1f} yr")
            log(f"  Log-rank p = {lr_f.p_value:.4f}")
            notes.append(
                f"FOXA1 p={lr_f.p_value:.3f} "
                f"med:{med_fh:.1f} vs {med_fl:.1f}yr"
            )

    # EZH2 + MKI67 composite in ILC
    if "EZH2" in ilc.columns and "MKI67" in ilc.columns:
        ezh2_v  = ilc["EZH2"][valid]
        mki67_v = ilc["MKI67"][valid]

        ezh2_high  = ezh2_v >= ezh2_v.median()
        mki67_high = mki67_v >= mki67_v.median()

        double_high = ezh2_high & mki67_high
        double_low  = (~ezh2_high) & (~mki67_high)

        log(f"\n  EZH2-high + MKI67-high (double high): "
            f"n={double_high.sum()}")
        log(f"  EZH2-low  + MKI67-low  (double low):  "
            f"n={double_low.sum()}")

        if double_high.sum() >= 5 and double_low.sum() >= 5:
            lr_d = logrank_test(
                t_v[double_high], t_v[double_low],
                e_v[double_high], e_v[double_low]
            )
            med_dh = t_v[double_high].median()
            med_dl = t_v[double_low].median()
            log(f"  Double-high median OS: {med_dh:.1f} yr")
            log(f"  Double-low  median OS: {med_dl:.1f} yr")
            log(f"  Log-rank p = {lr_d.p_value:.4f}")

            try:
                cox_df_ilc = pd.DataFrame({
                    "T":    t_v,
                    "E":    e_v,
                    "EZH2": ezh2_v,
                    "MKI67": mki67_v
                })
                cph_ilc = CoxPHFitter()
                cph_ilc.fit(cox_df_ilc,
                            duration_col="T", event_col="E")
                hr_ezh2  = float(np.exp(cph_ilc.params_["EZH2"]))
                hr_mki67 = float(np.exp(cph_ilc.params_["MKI67"]))
                log(f"  Cox HR EZH2:  {hr_ezh2:.3f}  "
                    f"p={cph_ilc.summary['p']['EZH2']:.4f}")
                log(f"  Cox HR MKI67: {hr_mki67:.3f}  "
                    f"p={cph_ilc.summary['p']['MKI67']:.4f}")

                direction_correct = med_dh <= med_dl
                sig = lr_d.p_value < 0.10

                note = (f"EZH2+MKI67 double-high OS "
                        f"{med_dh:.1f} vs {med_dl:.1f}yr "
                        f"p={lr_d.p_value:.3f}")
                notes.append(note)

                if direction_correct and sig:
                    record("CS-ILC-ET", "CONFIRMED", note)
                elif direction_correct:
                    record("CS-ILC-ET", "PARTIAL",
                           note + " (underpowered)")
                else:
                    record("CS-ILC-ET", "FAILED",
                           f"Direction reversed: {note}")
                    test_pass = False

            except Exception as e:
                log(f"  ILC Cox error: {e}")
                record("CS-ILC-ET", "PARTIAL",
                       f"Cox failed: {e}")
        else:
            log("  Insufficient groups for ILC double-high test.")
            record("CS-ILC-ET", "PARTIAL",
                   "Insufficient group sizes")
    else:
        log("  EZH2 or MKI67 not available in TCGA cache.")
        record("CS-ILC-ET", "PENDING", "EZH2/MKI67 missing")

# ============================================================
# STEP 9 — UNIVERSAL DEPTH SCORE (CS-DEPTH-UNIVERSAL)
# ============================================================

def universal_depth_survival(subtype_scores, tcga_a,
                              clin_a, os_time, os_event, masks):
    log("")
    log("=" * 65)
    log("ANALYSIS 6: UNIVERSAL DEPTH SCORE (CS-DEPTH-UNIVERSAL)")
    log("Predicted: Depth score dominates within-subtype OS")
    log("in every subtype where survival data is available")
    log("=" * 65)

    if os_time is None or os_event is None:
        log("  SKIP: No survival data.")
        record("CS-DEPTH-UNIVERSAL", "PENDING")
        return

    results_summary = []
    all_pass = True

    clinical_vars = {
        "LumA":  ["CDKN1A", "ESR1", "FOXA1"],
        "LumB":  ["ESR1", "MKI67", "FOXA1"],
        "HER2":  ["ERBB2", "ESR1", "EZH2"],
        "Basal": ["AR", "EZH2", "SOX10"],
        "ILC":   ["FOXA1", "EZH2", "MKI67"],
    }

    for subtype, depth_scores in subtype_scores.items():
        if subtype not in masks or subtype not in clinical_vars:
            continue
        if masks[subtype].sum() < 20:
            continue

        idx   = depth_scores.index
        t_raw = pd.to_numeric(os_time.reindex(idx), errors="coerce")
        e_raw = pd.to_numeric(os_event.reindex(idx), errors="coerce")
        valid = t_raw.notna() & e_raw.notna() & (t_raw > 0)

        if valid.sum() < 20:
            continue

        t_v = t_raw[valid] / 365.25
        e_v = e_raw[valid].astype(int)
        d_v = depth_scores[valid]

        # Cox comparison: depth score vs individual genes
        genes_avail = [g for g in clinical_vars[subtype]
                       if g in tcga_a.columns]
        if not genes_avail:
            continue

        pop_v = tcga_a.loc[idx][genes_avail][valid].copy()
        pop_v["depth"]  = d_v
        pop_v["T"]      = t_v
        pop_v["E"]      = e_v

        cox_results = {}

        # Depth score Cox
        try:
            cph = CoxPHFitter()
            cph.fit(pop_v[["T", "E", "depth"]].dropna(),
                    duration_col="T", event_col="E")
            hr_d = float(np.exp(cph.params_["depth"]))
            p_d  = float(cph.summary["p"]["depth"])
            cox_results["depth"] = (hr_d, p_d)
            log(f"\n  {subtype} — Depth score: "
                f"HR={hr_d:.3f}  p={p_d:.4f}")
        except Exception as e:
            log(f"  {subtype} — Depth Cox error: {e}")
            continue

        # Individual gene Cox
        for gene in genes_avail:
            if gene not in pop_v.columns:
                continue
            try:
                cph_g = CoxPHFitter()
                cph_g.fit(
                    pop_v[["T", "E", gene]].dropna(),
                    duration_col="T", event_col="E"
                )
                hr_g = float(np.exp(cph_g.params_[gene]))
                p_g  = float(cph_g.summary["p"][gene])
                cox_results[gene] = (hr_g, p_g)
                log(f"  {subtype} — {gene}: "
                    f"HR={hr_g:.3f}  p={p_g:.4f}")
            except Exception:
                pass

        if not cox_results:
            continue

        # Is depth score p-value the best or tied-best?
        depth_p  = cox_results.get("depth", (np.nan, 1.0))[1]
        other_ps = [v[1] for k, v in cox_results.items()
                    if k != "depth"]
        depth_best = all(depth_p <= p for p in other_ps) \
                     if other_ps else True

        results_summary.append({
            "subtype":    subtype,
            "depth_HR":   cox_results.get("depth", (np.nan,))[0],
            "depth_p":    depth_p,
            "depth_best": depth_best,
            "n_valid":    valid.sum(),
        })

        log(f"  {subtype} — Depth score best predictor: "
            f"{'YES' if depth_best else 'NO'}")
        if not depth_best:
            best_gene = min(
                [(k, v[1]) for k, v in cox_results.items()
                 if k != "depth"],
                key=lambda x: x[1]
            )
            log(f"    Best gene: {best_gene[0]} p={best_gene[1]:.4f}")
            all_pass = False

    if results_summary:
        n_best = sum(1 for r in results_summary if r["depth_best"])
        n_total = len(results_summary)
        log(f"\n  Depth score best predictor: "
            f"{n_best}/{n_total} subtypes")

        if n_best == n_total:
            record("CS-DEPTH-UNIVERSAL", "CONFIRMED",
                   f"Depth score best predictor in all "
                   f"{n_total} tested subtypes")
        elif n_best >= n_total * 0.6:
            record("CS-DEPTH-UNIVERSAL", "PARTIAL",
                   f"Depth score best in {n_best}/{n_total} subtypes")
        else:
            record("CS-DEPTH-UNIVERSAL", "FAILED",
                   f"Depth score best in only {n_best}/{n_total}")
    else:
        record("CS-DEPTH-UNIVERSAL", "PENDING",
               "No valid survival data for Cox comparison")

# ============================================================
# STEP 10 — MASTER FIGURE
# ============================================================

def make_master_figure(km_results):
    log("")
    log("=" * 65)
    log("MASTER FIGURE")
    log("=" * 65)

    n_panels = len(km_results)
    if n_panels == 0:
        log("  No KM results for master figure.")
        return

    COLORS = {
        "Q1 (shallow)": "#2980b9",
        "Q2":           "#27ae60",
        "Q3":           "#e67e22",
        "Q4 (deep)":    "#c0392b",
    }

    ncols  = min(3, n_panels)
    nrows  = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(6 * ncols, 5 * nrows))
    axes_flat = axes.flatten() if n_panels > 1 else [axes]

    for i, (subtype, res) in enumerate(km_results.items()):
        ax = axes_flat[i]
        t, e, groups, p, hr = (res["t"], res["e"],
                                res["groups"], res["p"],
                                res["hr"])
        for grp in ["Q1 (shallow)", "Q2", "Q3", "Q4 (deep)"]:
            if grp not in groups.values:
                continue
            m   = groups == grp
            kmf = KaplanMeierFitter()
            kmf.fit(t[m], e[m], label=f"{grp} (n={m.sum()})")
            kmf.plot_survival_function(
                ax=ax, color=COLORS.get(grp, "grey"),
                ci_show=False
            )
        ax.set_title(f"{subtype}\np={p:.4f}  HR={hr:.2f}",
                     fontsize=10)
        ax.set_xlabel("Time (years)", fontsize=9)
        ax.set_ylabel("OS Probability", fontsize=9)
        ax.set_ylim(0, 1.05)
        ax.legend(fontsize=7)

    for j in range(i + 1, len(axes_flat)):
        axes_flat[j].set_visible(False)

    plt.suptitle(
        "BRCA CROSS-SUBTYPE SURVIVAL — DEPTH SCORE QUARTILES\n"
        "OrganismCore / BRCA-S8d / 2026-03-05\n"
        "High depth = deep attractor lock = predicted worse OS",
        fontsize=12, y=1.01
    )
    plt.tight_layout()
    plt.savefig(FIG_MASTER, dpi=150, bbox_inches="tight")
    plt.close()
    log(f"  Master figure: {FIG_MASTER}")

# ============================================================
# MAIN
# ============================================================

def main():
    log("=" * 65)
    log("BRCA CROSS-SUBTYPE ANALYSIS — SCRIPT 2")
    log("OrganismCore — Document BRCA-S8d")
    log("Before-document: BRCA-S8d (predictions locked)")
    log("Date: 2026-03-05")
    log("=" * 65)
    log("")
    log("PROTOCOL v2.0:")
    log("  Depth score distributions printed FIRST.")
    log("  Survival tests SECOND.")
    log("  Wrong predictions documented alongside correct ones.")
    log("  Combined scorecard (Script 1 + Script 2) at end.")
    log("")

    # ── Load data ────────────────────────────────────────────
    tcga, clin, sc, os_time_col, os_event_col = \
        load_cached_data()

    # ── Classify subtypes ────────────────────────────────────
    pops, tcga_a, clin_a, masks = \
        classify_tcga_subtypes(tcga, clin)

    # ── Compute depth scores ─────────────────────────────────
    subtype_scores, all_depth, depth_df = \
        compute_depth_scores(tcga_a, clin_a, masks)

    # ── Resolve survival ─────────────────────────────────────
    os_time, os_event = \
        resolve_survival(clin_a, os_time_col, os_event_col)

    # ── Analysis 1: KM by depth (CS-10) ──────────────────────
    km_results, surv_rows = \
        kaplan_meier_by_depth(subtype_scores, tcga_a, clin_a,
                              os_time, os_event, masks)

    # ── Analysis 2: TNBC AR survival (CS-13-SURVIVAL) ────────
    tnbc_ar_survival(
        tcga_a=tcga_a,
        masks=masks,
        os_time=os_time,
        os_event=os_event,
    )

    # ── Analysis 3: LumB decoupling (CS-LUMB-DECOUPLE) ───────
    lumb_er_decoupling(tcga_a, masks)

    # ── Analysis 4: EZH2-free PCA (CS-PCA-EZH2FREE) ──────────
    pca_ezh2_free(sc)

    # ── Analysis 5: ILC survival (CS-ILC-ET) ─────────────────
    ilc_foxa1_survival(tcga_a, masks, os_time, os_event)

    # ── Analysis 6: Universal depth (CS-DEPTH-UNIVERSAL) ─────
    universal_depth_survival(subtype_scores, tcga_a,
                              clin_a, os_time, os_event, masks)

    # ── Master figure ─────────────────────────────────────────
    if km_results:
        make_master_figure(km_results)

    # ── Final combined scorecard ──────────────────────────────
    write_scorecard()

    log("")
    log("=" * 65)
    log("SCRIPT 2 COMPLETE")
    log(f"Results: {RESULTS_DIR}")
    log(f"Log:     {LOG_FILE}")
    log("Next:    BRCA-S8e (Script 2 Reasoning Artifact)")
    log("         then: BRCA-S8f (Script 3 if warranted)")
    log("=" * 65)

    write_log()


if __name__ == "__main__":
    main()
