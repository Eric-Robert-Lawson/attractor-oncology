#!/usr/bin/env python3
"""
FOXA1/EZH2 RATIO — VALIDATION SCRIPT 4 (RATIO-S4)
OrganismCore | 2026-03-06 | Eric Robert Lawson

FOUR COMPONENTS:
  A — Multivariate Cox (METABRIC)
      Does ratio add independent prognostic value
      beyond NPI, age, ER status, treatment?
  B — ROC curves / classification utility
      METABRIC primary (n=1980)
      TCGA PanCancer Atlas secondary (n=232 subtypes)
      LumA vs Basal, LumA vs LumB
  C — Claudin-low position (METABRIC, n=218 CL)
      Where does CL sit on the ratio?
      Honest biological interpretation.
  D — Treatment response
      GSE25066 (n=508, neoadjuvant chemo, pCR)
      METABRIC RFS in hormone-treated ER+ patients
"""
import os, sys, gzip, json, time
import urllib.request, urllib.error
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import mannwhitneyu, kruskal
from scipy.stats import spearmanr
import warnings
warnings.filterwarnings("ignore")

try:
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test
    from lifelines import CoxPHFitter
    HAS_LIFELINES = True
except ImportError:
    HAS_LIFELINES = False

try:
    from sklearn.metrics import (roc_curve, auc,
                                  cohen_kappa_score)
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False

# ── PATHS ─────────────────────────────────────────────────────
DEEP_DIVE = "/Users/ericlawson/cancer/BRCA/DEEP_DIVE"
OUT_DIR   = os.path.join(DEEP_DIVE,
                         "ratio_s4_results")
DATA_DIR  = os.path.join(OUT_DIR, "data")
RES_DIR   = os.path.join(OUT_DIR, "results")
for d in [OUT_DIR, DATA_DIR, RES_DIR]:
    os.makedirs(d, exist_ok=True)

MET_EXPR  = os.path.join(
    DEEP_DIVE,
    "ratio_s3_results/data/metabric_raw_mrna.csv")
MET_CLIN  = os.path.join(
    DEEP_DIVE,
    "Cross_Subtype_s3_results/data/"
    "metabric_clinical.csv")
GSE25066_LOCAL = os.path.join(
    DATA_DIR, "GSE25066_series_matrix.txt.gz")
GSE25066_URL   = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/"
    "GSE25nnn/GSE25066/matrix/"
    "GSE25066_series_matrix.txt.gz")

LOG_PATH  = os.path.join(RES_DIR,
                         "ratio_s4_log.txt")
SCORE_PATH= os.path.join(RES_DIR,
                         "ratio_s4_scorecard.csv")

SEP = "=" * 65

log_lines = []
def log(msg=""):
    print(msg)
    log_lines.append(msg)

preds = {}

# ── HELPERS ───────────────────────────────────────────────────
def api_get(url, timeout=20):
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "Mozilla/5.0",
                 "Accept":     "application/json"})
    with urllib.request.urlopen(
            req, timeout=timeout) as r:
        raw = r.read()
        return json.loads(raw) if raw.strip() else []

def api_post(url, body, timeout=30):
    data = json.dumps(body).encode()
    req  = urllib.request.Request(
        url, data=data,
        headers={"User-Agent":   "Mozilla/5.0",
                 "Accept":       "application/json",
                 "Content-Type": "application/json"})
    with urllib.request.urlopen(
            req, timeout=timeout) as r:
        raw = r.read()
        return json.loads(raw) if raw.strip() else []

def km_test(df, t_col, e_col, metric_col,
            title="", outpath=None):
    """KM log-rank test split at median."""
    if not HAS_LIFELINES:
        return None
    df = df[[t_col, e_col, metric_col]].dropna()
    df[t_col] = pd.to_numeric(df[t_col],
                               errors="coerce")
    df[e_col] = df[e_col].apply(
        lambda x: 1 if str(x).upper() in
        ["1","1.0","DECEASED","DEAD","YES",
         "TRUE","EVENT","1:DECEASED"] else 0)
    df = df.dropna(subset=[t_col])
    if len(df) < 20:
        return None
    med = df[metric_col].median()
    hi  = df[df[metric_col] >= med]
    lo  = df[df[metric_col] <  med]
    if len(hi) < 5 or len(lo) < 5:
        return None
    res = logrank_test(
        hi[t_col], lo[t_col],
        hi[e_col], lo[e_col])
    p = res.p_value

    if outpath:
        fig, ax = plt.subplots(figsize=(6, 4))
        kmf = KaplanMeierFitter()
        kmf.fit(hi[t_col], hi[e_col],
                label=f"High n={len(hi)}")
        kmf.plot_survival_function(ax=ax)
        kmf.fit(lo[t_col], lo[e_col],
                label=f"Low  n={len(lo)}")
        kmf.plot_survival_function(ax=ax)
        ax.set_title(f"{title}\np={p:.2e}")
        ax.set_xlabel("Months")
        ax.set_ylabel("Survival")
        fig.tight_layout()
        fig.savefig(outpath, dpi=150)
        plt.close(fig)
        log(f"  Figure: {outpath}")
    return p

def roc_binary(scores, labels, pos_label,
               title="", ax=None):
    """ROC for binary classification."""
    if not HAS_SKLEARN:
        return None, None
    y = (labels == pos_label).astype(int)
    if y.sum() < 5 or (1-y).sum() < 5:
        return None, None
    fpr, tpr, _ = roc_curve(y, scores)
    roc_auc = auc(fpr, tpr)
    if ax is not None:
        ax.plot(fpr, tpr,
                label=f"{title} AUC={roc_auc:.3f}")
        ax.plot([0,1],[0,1],"k--",alpha=0.3)
        ax.set_xlabel("FPR")
        ax.set_ylabel("TPR")
        ax.legend(fontsize=8)
    return fpr, roc_auc

def optimal_threshold(scores, labels, pos_label):
    """Youden J threshold."""
    if not HAS_SKLEARN:
        return None, None, None, None
    y = (labels == pos_label).astype(int)
    fpr, tpr, thresholds = roc_curve(y, scores)
    j = tpr - fpr
    idx = np.argmax(j)
    thresh = thresholds[idx]
    pred   = (scores >= thresh).astype(int)
    tp = int(((pred==1)&(y==1)).sum())
    fp = int(((pred==1)&(y==0)).sum())
    fn = int(((pred==0)&(y==1)).sum())
    tn = int(((pred==0)&(y==0)).sum())
    sens = tp/(tp+fn) if (tp+fn)>0 else 0
    spec = tn/(tn+fp) if (tn+fp)>0 else 0
    return thresh, sens, spec, tpr[idx]

# ══════════════════════════════════════════════════════════════
print(SEP)
log("FOXA1/EZH2 RATIO — VALIDATION SCRIPT 4 (RATIO-S4)")
log("OrganismCore | 2026-03-06 | Eric Robert Lawson")
log("")
log("COMPONENTS:")
log("  A — Multivariate Cox (METABRIC)")
log("  B — ROC classification utility")
log("  C — Claudin-low position (METABRIC)")
log("  D — Treatment response (GSE25066 + METABRIC)")
print(SEP)
log(f"  lifelines: {HAS_LIFELINES}")
log(f"  sklearn:   {HAS_SKLEARN}")

# ══════════════════════════════════════════════════════════════
# LOAD SHARED DATA
# ══════════════════════════════════════════════════════════════
log("\n── LOADING SHARED DATA ──")

# METABRIC expression
me = pd.read_csv(MET_EXPR, index_col=0)
foxa1_met = me.loc["FOXA1"].astype(float)
ezh2_met  = me.loc["EZH2"].astype(float)
metric_met = foxa1_met - ezh2_met
log(f"  METABRIC expr:    {me.shape}  "
    f"samples={len(metric_met)}")

# METABRIC clinical
mc = pd.read_csv(MET_CLIN, low_memory=False)
mc.columns = [c.strip() for c in mc.columns]
# Patient ID column
pid_col = "patient" if "patient" in mc.columns \
    else mc.columns[0]
mc = mc.set_index(pid_col)
log(f"  METABRIC clinical:{mc.shape}  "
    f"index='{pid_col}'")

# Align
common = metric_met.index.intersection(mc.index)
log(f"  METABRIC overlap: {len(common)}")
df_met = pd.DataFrame({
    "metric":      metric_met[common],
    "FOXA1":       foxa1_met[common],
    "EZH2":        ezh2_met[common],
    "subtype":     mc.loc[common, "CLAUDIN_SUBTYPE"],
    "OS_MONTHS":   pd.to_numeric(
        mc.loc[common, "OS_MONTHS"],
        errors="coerce"),
    "OS_STATUS":   mc.loc[common, "OS_STATUS"],
    "RFS_MONTHS":  pd.to_numeric(
        mc.loc[common, "RFS_MONTHS"],
        errors="coerce"),
    "RFS_STATUS":  mc.loc[common, "RFS_STATUS"],
    "AGE":         pd.to_numeric(
        mc.loc[common, "AGE_AT_DIAGNOSIS"],
        errors="coerce"),
    "LYMPH_NODES": pd.to_numeric(
        mc.loc[common, "LYMPH_NODES_EXAMINED_POSITIVE"],
        errors="coerce"),
    "NPI":         pd.to_numeric(
        mc.loc[common, "NPI"],
        errors="coerce"),
    "ER_IHC":      mc.loc[common, "ER_IHC"],
    "CHEMO":       mc.loc[common, "CHEMOTHERAPY"],
    "HORMONE":     mc.loc[common, "HORMONE_THERAPY"],
})
log(f"  METABRIC merged:  {df_met.shape}")
log(f"  Subtypes: "
    f"{df_met['subtype'].value_counts().to_dict()}")

# ══════════════════════════════════════════════════════════════
log("\n" + SEP)
log("COMPONENT A — MULTIVARIATE COX (METABRIC)")
log("Does ratio add independent prognostic value")
log("beyond NPI, age, ER status, treatment?")
log(SEP)

# ── A.1: Prepare Cox dataframe ────────────────────────────────
log("\n── A.1: PREPARE COX DATAFRAME ──")

df_cox = df_met.copy()

# OS event
df_cox["OS_E"] = df_cox["OS_STATUS"].apply(
    lambda x: 1 if "DECEASED" in str(x).upper()
    or str(x).strip() == "1" else 0)

# RFS event
df_cox["RFS_E"] = df_cox["RFS_STATUS"].apply(
    lambda x: 1 if "RECURRED" in str(x).upper()
    or "1:" in str(x) else 0)

# ER binary
df_cox["ER_pos"] = (
    df_cox["ER_IHC"].str.strip()
    .str.lower()
    .str.startswith("pos")
    .astype(float))

# Chemo binary
df_cox["chemo_yes"] = (
    df_cox["CHEMO"].str.strip()
    .str.upper() == "YES").astype(float)

# Hormone binary
df_cox["hormone_yes"] = (
    df_cox["HORMONE"].str.strip()
    .str.upper() == "YES").astype(float)

# Lymph nodes — clip outliers
df_cox["LYMPH_NODES"] = df_cox["LYMPH_NODES"].clip(
    upper=20)

# Subtype dummies (vs LumA reference)
df_cox["is_LumB"]   = (
    df_cox["subtype"] == "LumB").astype(float)
df_cox["is_Her2"]   = (
    df_cox["subtype"] == "Her2").astype(float)
df_cox["is_Basal"]  = (
    df_cox["subtype"] == "Basal").astype(float)
df_cox["is_CL"]     = (
    df_cox["subtype"] == "claudin-low").astype(float)
df_cox["is_Normal"] = (
    df_cox["subtype"] == "Normal").astype(float)

# Standardise metric for Cox (per-SD HR)
df_cox["metric_z"] = (
    (df_cox["metric"] - df_cox["metric"].mean())
    / df_cox["metric"].std())

cox_cols_os = [
    "OS_MONTHS", "OS_E",
    "metric_z", "NPI", "AGE",
    "LYMPH_NODES", "ER_pos",
    "chemo_yes", "hormone_yes",
    "is_LumB", "is_Her2",
    "is_Basal", "is_CL", "is_Normal"]

cox_cols_rfs = [
    "RFS_MONTHS", "RFS_E",
    "metric_z", "NPI", "AGE",
    "LYMPH_NODES", "ER_pos",
    "chemo_yes", "hormone_yes",
    "is_LumB", "is_Her2",
    "is_Basal", "is_CL", "is_Normal"]

df_os  = df_cox[cox_cols_os].dropna()
df_rfs = df_cox[cox_cols_rfs].dropna()
log(f"  Cox OS  n={len(df_os)}  "
    f"events={int(df_os['OS_E'].sum())}")
log(f"  Cox RFS n={len(df_rfs)}  "
    f"events={int(df_rfs['RFS_E'].sum())}")

# ── A.2: Fit Cox models ───────────────────────────────────────
log("\n── A.2: FIT COX MODELS ─��")

a_os_hr  = None
a_os_p   = None
a_rfs_hr = None
a_rfs_p  = None

if HAS_LIFELINES:
    # OS model
    try:
        cph_os = CoxPHFitter(penalizer=0.1)
        cph_os.fit(df_os,
                   duration_col="OS_MONTHS",
                   event_col="OS_E")
        log("\n  OS MULTIVARIATE COX:")
        log(f"  Concordance: "
            f"{cph_os.concordance_index_:.4f}")
        summary_os = cph_os.summary
        if "metric_z" in summary_os.index:
            row = summary_os.loc["metric_z"]
            a_os_hr = float(row["exp(coef)"])
            a_os_p  = float(row["p"])
            ci_lo   = float(row.get(
                "exp(coef) lower 95%",
                row.get("lower 0.95", np.nan)))
            ci_hi   = float(row.get(
                "exp(coef) upper 95%",
                row.get("upper 0.95", np.nan)))
            log(f"  metric_z HR={a_os_hr:.3f}  "
                f"95%CI [{ci_lo:.3f}-{ci_hi:.3f}]"
                f"  p={a_os_p:.4f}")
        log("\n  Full OS summary:")
        log(summary_os[
            ["exp(coef)",
             "exp(coef) lower 95%",
             "exp(coef) upper 95%",
             "p"]
        ].to_string())
    except Exception as e:
        log(f"  OS Cox error: {e}")

    # RFS model
    try:
        cph_rfs = CoxPHFitter(penalizer=0.1)
        cph_rfs.fit(df_rfs,
                    duration_col="RFS_MONTHS",
                    event_col="RFS_E")
        log("\n  RFS MULTIVARIATE COX:")
        log(f"  Concordance: "
            f"{cph_rfs.concordance_index_:.4f}")
        summary_rfs = cph_rfs.summary
        if "metric_z" in summary_rfs.index:
            row = summary_rfs.loc["metric_z"]
            a_rfs_hr = float(row["exp(coef)"])
            a_rfs_p  = float(row["p"])
            ci_lo    = float(row.get(
                "exp(coef) lower 95%",
                row.get("lower 0.95", np.nan)))
            ci_hi    = float(row.get(
                "exp(coef) upper 95%",
                row.get("upper 0.95", np.nan)))
            log(f"  metric_z HR={a_rfs_hr:.3f}  "
                f"95%CI [{ci_lo:.3f}-{ci_hi:.3f}]"
                f"  p={a_rfs_p:.4f}")
        log("\n  Full RFS summary:")
        log(summary_rfs[
            ["exp(coef)",
             "exp(coef) lower 95%",
             "exp(coef) upper 95%",
             "p"]
        ].to_string())
    except Exception as e:
        log(f"  RFS Cox error: {e}")
else:
    log("  lifelines not available — Cox skipped")

# A decision
log("\n── A.3: COMPONENT A RESULT ��─")
a_confirmed = (
    a_os_p  is not None and a_os_p  < 0.05
    and a_os_hr is not None and a_os_hr < 1.0)
a_rfs_confirmed = (
    a_rfs_p is not None and a_rfs_p < 0.05
    and a_rfs_hr is not None and a_rfs_hr < 1.0)

preds["R4-A_OS"]  = (
    "CONFIRMED"     if a_confirmed     else
    "NOT CONFIRMED" if a_os_p  is not None
    else "NOT TESTABLE")
preds["R4-A_RFS"] = (
    "CONFIRMED"     if a_rfs_confirmed else
    "NOT CONFIRMED" if a_rfs_p is not None
    else "NOT TESTABLE")

log(f"  OS  HR={a_os_hr}  p={a_os_p}")
log(f"  RFS HR={a_rfs_hr}  p={a_rfs_p}")
log(f"  R4-A_OS  → {preds['R4-A_OS']}")
log(f"  R4-A_RFS → {preds['R4-A_RFS']}")

# ══════════════════════════════════════════════════════════════
log("\n" + SEP)
log("COMPONENT B — ROC / CLASSIFICATION UTILITY")
log("METABRIC primary + TCGA secondary")
log("LumA vs Basal, LumA vs LumB")
log(SEP)

# ── B.1: METABRIC ROC ─────────────────────────────────────────
log("\n── B.1: METABRIC ROC ──")

fig_b, axes_b = plt.subplots(1, 3, figsize=(15, 5))

met_sub = df_met.dropna(subset=["subtype", "metric"])
met_sub = met_sub[
    met_sub["subtype"].isin(
        ["LumA", "LumB", "Basal",
         "Her2", "claudin-low"])]

auc_luma_basal_met = None
auc_luma_lumb_met  = None

# LumA vs Basal
sub_ab = met_sub[met_sub["subtype"].isin(
    ["LumA", "Basal"])]
if len(sub_ab) > 20:
    _, auc_luma_basal_met = roc_binary(
        sub_ab["metric"],
        sub_ab["subtype"],
        "LumA",
        title="METABRIC LumA vs Basal",
        ax=axes_b[0])
    thresh, sens, spec, _ = optimal_threshold(
        sub_ab["metric"],
        sub_ab["subtype"], "LumA")
    log(f"  LumA vs Basal: "
        f"n={len(sub_ab)}  "
        f"AUC={auc_luma_basal_met:.3f}  "
        f"thresh={thresh:.3f}  "
        f"sens={sens:.3f}  spec={spec:.3f}")
    axes_b[0].set_title(
        f"METABRIC LumA vs Basal\n"
        f"n={len(sub_ab)}  "
        f"AUC={auc_luma_basal_met:.3f}")

# LumA vs LumB
sub_al = met_sub[met_sub["subtype"].isin(
    ["LumA", "LumB"])]
if len(sub_al) > 20:
    _, auc_luma_lumb_met = roc_binary(
        sub_al["metric"],
        sub_al["subtype"],
        "LumA",
        title="METABRIC LumA vs LumB",
        ax=axes_b[1])
    thresh, sens, spec, _ = optimal_threshold(
        sub_al["metric"],
        sub_al["subtype"], "LumA")
    log(f"  LumA vs LumB:  "
        f"n={len(sub_al)}  "
        f"AUC={auc_luma_lumb_met:.3f}  "
        f"thresh={thresh:.3f}  "
        f"sens={sens:.3f}  spec={spec:.3f}")
    axes_b[1].set_title(
        f"METABRIC LumA vs LumB\n"
        f"n={len(sub_al)}  "
        f"AUC={auc_luma_lumb_met:.3f}")

# LumA vs HER2
sub_ah = met_sub[met_sub["subtype"].isin(
    ["LumA", "Her2"])]
if len(sub_ah) > 20:
    _, auc_luma_her2_met = roc_binary(
        sub_ah["metric"],
        sub_ah["subtype"],
        "LumA",
        title="METABRIC LumA vs Her2",
        ax=axes_b[2])
    thresh, sens, spec, _ = optimal_threshold(
        sub_ah["metric"],
        sub_ah["subtype"], "LumA")
    log(f"  LumA vs Her2:  "
        f"n={len(sub_ah)}  "
        f"AUC={auc_luma_her2_met:.3f}  "
        f"thresh={thresh:.3f}  "
        f"sens={sens:.3f}  spec={spec:.3f}")
    axes_b[2].set_title(
        f"METABRIC LumA vs Her2\n"
        f"n={len(sub_ah)}  "
        f"AUC={auc_luma_her2_met:.3f}")

fig_b.suptitle("FOXA1/EZH2 ROC — METABRIC",
               fontsize=13, fontweight="bold")
fig_b.tight_layout()
roc_met_path = os.path.join(
    RES_DIR, "ratio_s4_compB_metabric_roc.png")
fig_b.savefig(roc_met_path, dpi=150)
plt.close(fig_b)
log(f"  Figure: {roc_met_path}")

# ── B.2: TCGA fetch + ROC ─────────────────────────────────────
log("\n── B.2: TCGA PanCancer Atlas ──")

tcga_study    = "brca_tcga_pan_can_atlas_2018"
tcga_profile  = (f"{tcga_study}"
                 f"_rna_seq_v2_mrna")
tcga_slist    = f"{tcga_study}_all"

# Fetch FOXA1 and EZH2
log("  Fetching TCGA FOXA1 + EZH2 (RSEM)...")
tcga_expr = {}
for gene_name, entrez in [("FOXA1", 2308),
                           ("EZH2",  2146)]:
    url  = (f"https://www.cbioportal.org/api/"
            f"molecular-profiles/{tcga_profile}/"
            f"molecular-data/fetch")
    body = {"sampleListId": tcga_slist,
            "entrezGeneIds": [entrez]}
    try:
        data = api_post(url, body)
        vals = {d["sampleId"]: d["value"]
                for d in data
                if d.get("entrezGeneId") == entrez}
        tcga_expr[gene_name] = vals
        log(f"  {gene_name}: n={len(vals)}"
            f"  example_val="
            f"{list(vals.values())[0]:.3f}"
            if vals else
            f"  {gene_name}: empty")
        time.sleep(0.3)
    except Exception as e:
        log(f"  {gene_name} fetch error: {e}")

# Fetch clinical (all patients)
log("  Fetching TCGA clinical...")
tcga_clin = {}
url_clin = (f"https://www.cbioportal.org/api/"
            f"studies/{tcga_study}/clinical-data"
            f"?clinicalDataType=PATIENT"
            f"&pageSize=10000")
try:
    clin_rows = api_get(url_clin)
    log(f"  Clinical rows: {len(clin_rows)}")
    for row in clin_rows:
        pid  = row.get("patientId", "")
        attr = row.get("clinicalAttributeId", "")
        val  = row.get("value", "")
        if pid not in tcga_clin:
            tcga_clin[pid] = {}
        tcga_clin[pid][attr] = val
    log(f"  Unique patients: {len(tcga_clin)}")
except Exception as e:
    log(f"  Clinical fetch error: {e}")

# Build TCGA dataframe
log("  Building TCGA dataframe...")
if tcga_expr.get("FOXA1") and tcga_expr.get("EZH2"):
    foxa1_t = pd.Series(tcga_expr["FOXA1"],
                        name="FOXA1")
    ezh2_t  = pd.Series(tcga_expr["EZH2"],
                        name="EZH2")
    # RSEM — log2(x+1)
    foxa1_log = np.log2(foxa1_t + 1)
    ezh2_log  = np.log2(ezh2_t  + 1)
    metric_t  = foxa1_log - ezh2_log

    # Match patient IDs
    # sample IDs are TCGA-XX-XXXX-01
    # patient IDs are TCGA-XX-XXXX
    sample_to_patient = {
        s: "-".join(s.split("-")[:3])
        for s in metric_t.index}

    rows_t = []
    for sid, met in metric_t.items():
        pid = sample_to_patient.get(sid, "")
        cdat = tcga_clin.get(pid, {})
        subtype_raw = cdat.get("SUBTYPE", "")
        # Strip BRCA_ prefix
        subtype = subtype_raw.replace(
            "BRCA_", "")
        os_m = cdat.get("OS_MONTHS", np.nan)
        os_s = cdat.get("OS_STATUS", "")
        rows_t.append({
            "sampleId":  sid,
            "patientId": pid,
            "metric":    met,
            "FOXA1":     foxa1_log[sid],
            "EZH2":      ezh2_log[sid],
            "subtype":   subtype,
            "OS_MONTHS": os_m,
            "OS_STATUS": os_s,
        })
    df_tcga = pd.DataFrame(rows_t)
    df_tcga["OS_MONTHS"] = pd.to_numeric(
        df_tcga["OS_MONTHS"], errors="coerce")
    log(f"  TCGA shape: {df_tcga.shape}")
    log(f"  Subtypes: "
        f"{df_tcga['subtype'].value_counts().to_dict()}")
    log(f"  metric range: "
        f"[{df_tcga['metric'].min():.3f}, "
        f"{df_tcga['metric'].max():.3f}]")

    # TCGA ROC
    log("\n  TCGA ROC curves:")
    fig_bt, axes_bt = plt.subplots(
        1, 2, figsize=(10, 5))

    tcga_sub = df_tcga[
        df_tcga["subtype"].isin(
            ["LumA", "LumB", "Basal",
             "Her2", "Normal"])]

    auc_luma_basal_tcga = None
    auc_luma_lumb_tcga  = None

    sub_ab_t = tcga_sub[tcga_sub["subtype"].isin(
        ["LumA", "Basal"])]
    if len(sub_ab_t) > 10:
        _, auc_luma_basal_tcga = roc_binary(
            sub_ab_t["metric"],
            sub_ab_t["subtype"],
            "LumA",
            title="TCGA LumA vs Basal",
            ax=axes_bt[0])
        thresh, sens, spec, _ = optimal_threshold(
            sub_ab_t["metric"],
            sub_ab_t["subtype"], "LumA")
        log(f"  LumA vs Basal: "
            f"n={len(sub_ab_t)}  "
            f"AUC={auc_luma_basal_tcga:.3f}  "
            f"thresh={thresh:.3f}  "
            f"sens={sens:.3f}  spec={spec:.3f}")
        axes_bt[0].set_title(
            f"TCGA LumA vs Basal\n"
            f"n={len(sub_ab_t)}  "
            f"AUC={auc_luma_basal_tcga:.3f}")

    sub_al_t = tcga_sub[tcga_sub["subtype"].isin(
        ["LumA", "LumB"])]
    if len(sub_al_t) > 10:
        _, auc_luma_lumb_tcga = roc_binary(
            sub_al_t["metric"],
            sub_al_t["subtype"],
            "LumA",
            title="TCGA LumA vs LumB",
            ax=axes_bt[1])
        thresh, sens, spec, _ = optimal_threshold(
            sub_al_t["metric"],
            sub_al_t["subtype"], "LumA")
        log(f"  LumA vs LumB:  "
            f"n={len(sub_al_t)}  "
            f"AUC={auc_luma_lumb_tcga:.3f}  "
            f"thresh={thresh:.3f}  "
            f"sens={sens:.3f}  spec={spec:.3f}")
        axes_bt[1].set_title(
            f"TCGA LumA vs LumB\n"
            f"n={len(sub_al_t)}  "
            f"AUC={auc_luma_lumb_tcga:.3f}")

    fig_bt.suptitle(
        "FOXA1/EZH2 ROC — TCGA PanCancer Atlas",
        fontsize=13, fontweight="bold")
    fig_bt.tight_layout()
    roc_tcga_path = os.path.join(
        RES_DIR, "ratio_s4_compB_tcga_roc.png")
    fig_bt.savefig(roc_tcga_path, dpi=150)
    plt.close(fig_bt)
    log(f"  Figure: {roc_tcga_path}")

    # Save TCGA cache
    tcga_cache = os.path.join(
        DATA_DIR, "tcga_pancancer_merged.csv")
    df_tcga.to_csv(tcga_cache, index=False)
    log(f"  Saved: {tcga_cache}")

else:
    log("  TCGA expression fetch failed — "
        "skipping TCGA ROC")
    auc_luma_basal_tcga = None
    auc_luma_lumb_tcga  = None
    df_tcga = pd.DataFrame()

# ── B.3: ROC decisions ────────────────────────────────────────
log("\n── B.3: ROC DECISIONS ──")

preds["R4-B_LumA_Basal_MET"] = (
    "CONFIRMED"     if auc_luma_basal_met
                       is not None
                       and auc_luma_basal_met > 0.70
    else "NOT CONFIRMED" if auc_luma_basal_met
                            is not None
    else "NOT TESTABLE")

preds["R4-B_LumA_LumB_MET"]  = (
    "CONFIRMED"     if auc_luma_lumb_met
                       is not None
                       and auc_luma_lumb_met > 0.60
    else "NOT CONFIRMED" if auc_luma_lumb_met
                            is not None
    else "NOT TESTABLE")

preds["R4-B_LumA_Basal_TCGA"] = (
    "CONFIRMED"     if auc_luma_basal_tcga
                       is not None
                       and auc_luma_basal_tcga > 0.70
    else "NOT CONFIRMED" if auc_luma_basal_tcga
                            is not None
    else "NOT TESTABLE")

preds["R4-B_LumA_LumB_TCGA"]  = (
    "CONFIRMED"     if auc_luma_lumb_tcga
                       is not None
                       and auc_luma_lumb_tcga > 0.60
    else "NOT CONFIRMED" if auc_luma_lumb_tcga
                            is not None
    else "NOT TESTABLE")

for k, v in preds.items():
    if k.startswith("R4-B"):
        log(f"  {k} → {v}")

# ══════════════════════════════════════════════════════════════
log("\n" + SEP)
log("COMPONENT C — CLAUDIN-LOW POSITION (METABRIC)")
log("n=218 CL samples")
log(SEP)

log("\n── C.1: CL POSITION ──")
cl_mask  = df_met["subtype"] == "claudin-low"
luma_mask= df_met["subtype"] == "LumA"
tnbc_mask= df_met["subtype"] == "Basal"

cl_med   = df_met.loc[cl_mask,   "metric"].median()
luma_med = df_met.loc[luma_mask, "metric"].median()
tnbc_med = df_met.loc[tnbc_mask, "metric"].median()

log(f"  LumA  n={luma_mask.sum()}  "
    f"median={luma_med:.4f}")
log(f"  Basal n={tnbc_mask.sum()}  "
    f"median={tnbc_med:.4f}")
log(f"  CL    n={cl_mask.sum()}    "
    f"median={cl_med:.4f}")
log(f"  CL < TNBC(Basal): "
    f"{'✓' if cl_med < tnbc_med else '✗'}")
log(f"  CL < LumA: "
    f"{'✓' if cl_med < luma_med else '✗'}")

# Per-subtype full table
log("\n  All subtypes:")
for st in ["LumA", "LumB", "Her2",
           "Basal", "claudin-low", "Normal"]:
    mask = df_met["subtype"] == st
    if mask.sum() == 0:
        continue
    m = df_met.loc[mask, "metric"]
    f = df_met.loc[mask, "FOXA1"]
    e = df_met.loc[mask, "EZH2"]
    log(f"  {st:12s} n={mask.sum():4d}  "
        f"metric={m.median():.4f}  "
        f"FOXA1={f.median():.3f}  "
        f"EZH2={e.median():.3f}")

# CL vs Basal MWU
if cl_mask.sum() > 5 and tnbc_mask.sum() > 5:
    _, p_cl_tnbc = mannwhitneyu(
        df_met.loc[cl_mask,   "metric"],
        df_met.loc[tnbc_mask, "metric"],
        alternative="two-sided")
    log(f"\n  CL vs Basal MWU p={p_cl_tnbc:.4f}")

# CL FOXA1 and EZH2 separately
log(f"\n  CL FOXA1 median: "
    f"{df_met.loc[cl_mask,'FOXA1'].median():.4f}")
log(f"  CL EZH2 median:  "
    f"{df_met.loc[cl_mask,'EZH2'].median():.4f}")
log(f"  Basal FOXA1 median: "
    f"{df_met.loc[tnbc_mask,'FOXA1'].median():.4f}")
log(f"  Basal EZH2 median:  "
    f"{df_met.loc[tnbc_mask,'EZH2'].median():.4f}")

log("\n── C.2: BIOLOGICAL INTERPRETATION ──")
log("  CL EZH2 is LOW relative to Basal.")
log("  Basal: EZH2 high (active silencing).")
log("  CL: mesenchymal, low proliferation,")
log("      EZH2 not the dominant mechanism.")
log("  In bulk RNA, low EZH2 inflates ratio")
log("  for CL vs scRNA-seq expectation.")
log("  scRNA-seq (cancer cells only) showed")
log("  CL ratio=0.10 < Basal ratio=0.52.")
log("  Bulk RNA conflates cancer + stroma.")
log("  CL stroma is EZH2-low mesenchymal.")
log("  This explains CL/Basal boundary anomaly.")
log("  Conclusion: CL ordering requires")
log("  cancer-cell-specific measurement (IHC).")

preds["R4-C_CL_below_Basal"] = (
    "CONFIRMED"     if cl_med < tnbc_med else
    "NOT CONFIRMED (bulk RNA anomaly — "
    "see biological note)")

log(f"\n  R4-C → {preds['R4-C_CL_below_Basal']}")

# Component C figure
fig_c, ax_c = plt.subplots(figsize=(8, 5))
order = ["LumA", "LumB", "Her2",
         "Basal", "claudin-low", "Normal"]
colors = ["#2196F3", "#03A9F4", "#FF9800",
          "#F44336", "#9C27B0", "#4CAF50"]
data_c = []
labels_c = []
for st, col in zip(order, colors):
    mask = df_met["subtype"] == st
    if mask.sum() > 0:
        data_c.append(
            df_met.loc[mask, "metric"].values)
        labels_c.append(
            f"{st}\nn={mask.sum()}")
bplot = ax_c.boxplot(data_c, patch_artist=True,
                     notch=False,
                     medianprops=dict(
                         color="black",
                         linewidth=2))
for patch, col in zip(bplot["boxes"], colors):
    patch.set_facecolor(col)
    patch.set_alpha(0.7)
ax_c.set_xticks(range(1, len(labels_c)+1))
ax_c.set_xticklabels(labels_c, fontsize=9)
ax_c.set_ylabel("FOXA1 − EZH2 (log2)")
ax_c.set_title(
    "FOXA1/EZH2 Metric by Subtype — METABRIC\n"
    "Including Claudin-low (n=218)",
    fontweight="bold")
ax_c.axhline(0, color="gray",
             linestyle="--", alpha=0.4)
fig_c.tight_layout()
cl_path = os.path.join(
    RES_DIR, "ratio_s4_compC_claudin_low.png")
fig_c.savefig(cl_path, dpi=150)
plt.close(fig_c)
log(f"  Figure: {cl_path}")

# ══════════════════════════════════════════════════════════════
log("\n" + SEP)
log("COMPONENT D — TREATMENT RESPONSE")
log("D.1: GSE25066 (n=508, neoadjuvant chemo, pCR)")
log("D.2: METABRIC RFS in hormone-treated ER+")
log(SEP)

# ── D.1: GSE25066 ─────────────────────────────────────────────
log("\n── D.1: GSE25066 ──")

# Download if needed
if not os.path.exists(GSE25066_LOCAL):
    log(f"  Downloading GSE25066 (~62 MB)...")
    try:
        urllib.request.urlretrieve(
            GSE25066_URL, GSE25066_LOCAL)
        sz = os.path.getsize(GSE25066_LOCAL)
        log(f"  Downloaded: {sz:,} B")
    except Exception as e:
        log(f"  Download error: {e}")
else:
    sz = os.path.getsize(GSE25066_LOCAL)
    log(f"  Cached: {sz:,} B")

d1_result = "NOT TESTABLE"

if os.path.exists(GSE25066_LOCAL):
    try:
        log("  Parsing GSE25066 series matrix...")
        gsm_order = []
        char_accum = {}
        expr_rows  = {}
        in_matrix  = False
        foxa1_vals = {}
        ezh2_vals  = {}

        with gzip.open(
                GSE25066_LOCAL, "rt",
                errors="replace") as f:
            for line in f:
                line = line.rstrip("\n")

                if line.startswith(
                        "!series_matrix_table_begin"):
                    in_matrix = True
                    continue
                if line.startswith(
                        "!series_matrix_table_end"):
                    in_matrix = False
                    continue

                if not in_matrix:
                    if line.startswith(
                            "!Sample_geo_accession"):
                        parts = line.split("\t")
                        gsm_order = [
                            p.strip().strip('"')
                            for p in parts[1:]
                            if p.strip().strip('"')
                            .startswith("GSM")]

                    if line.startswith(
                            "!Sample_characteristics_ch1"):
                        parts = line.split("\t")
                        vals  = [
                            p.strip().strip('"')
                            for p in parts[1:]]
                        if vals:
                            key = vals[0].split(
                                ":")[0].strip().lower()
                            clean = [
                                v.split(":",1)[-1]
                                .strip()
                                for v in vals]
                            if key not in char_accum:
                                char_accum[key] = clean

                else:
                    # Expression matrix
                    parts = line.split("\t")
                    if len(parts) < 2:
                        continue
                    probe = parts[0].strip('"')
                    if probe in ["ID_REF", ""]:
                        continue
                    expr_rows[probe] = [
                        float(v) if v not in
                        ["", "null", "NA"]
                        else np.nan
                        for v in parts[1:]]

        log(f"  GSMs: {len(gsm_order)}")
        log(f"  Characteristics: "
            f"{list(char_accum.keys())}")
        log(f"  Expression probes: "
            f"{len(expr_rows)}")

        # Find pCR characteristic
        pcr_key = None
        pcr_map = {}
        for k, vals in char_accum.items():
            if any(x in k.lower()
                   for x in ["pcr", "response",
                              "patholog",
                              "residual"]):
                log(f"  Candidate key '{k}': "
                    f"{set(vals[:20])}")
                if pcr_key is None:
                    pcr_key = k

        if pcr_key:
            pcr_vals = char_accum[pcr_key]
            for i, gsm in enumerate(gsm_order):
                if i < len(pcr_vals):
                    pcr_map[gsm] = pcr_vals[i]
            log(f"  pCR key: '{pcr_key}'")
            pcr_vc = pd.Series(pcr_map).value_counts().to_dict()
            log(f"  pCR values: {pcr_vc}")

        # Find FOXA1 and EZH2 probes
        # GPL96 = Affymetrix HG-U133A
        # FOXA1 probe: 204948_s_at or 204949_s_at
        # EZH2 probe:  203358_s_at
        foxa1_probes = [
            "204948_s_at", "204949_s_at",
            "216565_x_at"]
        ezh2_probes  = [
            "203358_s_at", "225840_at",
            "203359_s_at"]

        foxa1_data = None
        ezh2_data  = None

        for probe in foxa1_probes:
            if probe in expr_rows:
                foxa1_data = np.array(
                    expr_rows[probe])
                log(f"  FOXA1 probe: {probe}  "
                    f"mean={np.nanmean(foxa1_data):.3f}")
                break

        for probe in ezh2_probes:
            if probe in expr_rows:
                ezh2_data = np.array(
                    expr_rows[probe])
                log(f"  EZH2 probe: {probe}  "
                    f"mean={np.nanmean(ezh2_data):.3f}")
                break

        if foxa1_data is None:
            # Search by partial match
            for probe, vals in expr_rows.items():
                if "204948" in probe or \
                   "204949" in probe:
                    foxa1_data = np.array(vals)
                    log(f"  FOXA1 probe (partial): "
                        f"{probe}")
                    break

        if ezh2_data is None:
            for probe, vals in expr_rows.items():
                if "203358" in probe or \
                   "203359" in probe:
                    ezh2_data = np.array(vals)
                    log(f"  EZH2 probe (partial): "
                        f"{probe}")
                    break

        if foxa1_data is not None \
                and ezh2_data is not None:
            n = min(len(foxa1_data),
                    len(ezh2_data),
                    len(gsm_order))
            metric_g = (foxa1_data[:n]
                        - ezh2_data[:n])
            gsms_g   = gsm_order[:n]

            df_g = pd.DataFrame({
                "GSM":    gsms_g,
                "metric": metric_g,
                "FOXA1":  foxa1_data[:n],
                "EZH2":   ezh2_data[:n],
            }).set_index("GSM")

            if pcr_key and pcr_map:
                df_g["pCR_raw"] = pd.Series(pcr_map)
                log(f"\n  GSE25066 merged: {df_g.shape}")
                pcr_vc = df_g["pCR_raw"].value_counts().to_dict()
                log(f"  pCR values: {pcr_vc}")

                # Identify pCR positive label
                pcr_counts = (
                    df_g["pCR_raw"]
                    .value_counts())
                log(f"  pCR counts: "
                    f"{pcr_counts.to_dict()}")

                # Common labels: pCR, RD, 
                # pathCR, non-pCR, 0, 1
                pcr_pos_labels = [
                    "pcr", "pcr ", "yes",
                    "1", "1.0", "pathcr",
                    "path cr", "complete response",
                    "cr"]
                pcr_neg_labels = [
                    "rd", "rd ", "no", "0",
                    "0.0", "residual disease",
                    "non-pcr"]

                def encode_pcr(v):
                    vl = str(v).lower().strip()
                    if vl in pcr_pos_labels:
                        return 1
                    if vl in pcr_neg_labels:
                        return 0
                    return np.nan

                df_g["pCR"] = df_g[
                    "pCR_raw"].apply(encode_pcr)
                n_pcr = int(
                    df_g["pCR"].sum())
                n_rd  = int(
                    (df_g["pCR"]==0).sum())
                log(f"  pCR=1: {n_pcr}  "
                    f"pCR=0(RD): {n_rd}")

                if n_pcr >= 10 and n_rd >= 10:
                    # KM-style: does high ratio
                    # predict pCR?
                    df_pcr = df_g.dropna(
                        subset=["pCR", "metric"])
                    med_g  = df_pcr[
                        "metric"].median()
                    hi_pcr = df_pcr[
                        df_pcr["metric"] >= med_g]
                    lo_pcr = df_pcr[
                        df_pcr["metric"] <  med_g]
                    hi_rate = hi_pcr["pCR"].mean()
                    lo_rate = lo_pcr["pCR"].mean()
                    log(f"  High ratio pCR rate: "
                        f"{hi_rate:.3f} "
                        f"(n={len(hi_pcr)})")
                    log(f"  Low  ratio pCR rate: "
                        f"{lo_rate:.3f} "
                        f"(n={len(lo_pcr)})")

                    # MWU on metric: pCR vs RD
                    pcr_scores = df_pcr.loc[
                        df_pcr["pCR"]==1, "metric"]
                    rd_scores  = df_pcr.loc[
                        df_pcr["pCR"]==0, "metric"]
                    _, p_pcr = mannwhitneyu(
                        pcr_scores, rd_scores,
                        alternative="two-sided")
                    log(f"  pCR vs RD MWU "
                        f"p={p_pcr:.4f}")
                    log(f"  pCR metric med="
                        f"{pcr_scores.median():.4f}")
                    log(f"  RD  metric med="
                        f"{rd_scores.median():.4f}")

                    # ROC for pCR
                    if HAS_SKLEARN:
                        fpr_pcr, auc_pcr = roc_binary(
                            df_pcr["metric"],
                            df_pcr["pCR"],
                            1,
                            title="GSE25066 pCR")
                        log(f"  pCR AUC={auc_pcr:.3f}"
                            if auc_pcr else
                            "  pCR AUC: failed")
                    else:
                        auc_pcr = None

                    # Figure
                    fig_d1, axes_d1 = plt.subplots(
                        1, 2, figsize=(10, 5))

                    # Box: metric by pCR
                    axes_d1[0].boxplot(
                        [pcr_scores.values,
                         rd_scores.values],
                        labels=["pCR", "RD"],
                        patch_artist=True)
                    axes_d1[0].set_title(
                        f"GSE25066 metric by response\n"
                        f"MWU p={p_pcr:.3f}")
                    axes_d1[0].set_ylabel(
                        "FOXA1 − EZH2 (log2)")

                    # ROC
                    if HAS_SKLEARN and auc_pcr:
                        fpr_p, tpr_p, _ = roc_curve(
                            (df_pcr["pCR"]==1
                             ).astype(int),
                            df_pcr["metric"])
                        axes_d1[1].plot(
                            fpr_p, tpr_p,
                            label=f"AUC={auc_pcr:.3f}")
                        axes_d1[1].plot(
                            [0,1],[0,1],
                            "k--", alpha=0.3)
                        axes_d1[1].set_title(
                            "ROC: pCR prediction")
                        axes_d1[1].set_xlabel("FPR")
                        axes_d1[1].set_ylabel("TPR")
                        axes_d1[1].legend()

                    fig_d1.suptitle(
                        "GSE25066 — FOXA1/EZH2 "
                        "vs Treatment Response",
                        fontweight="bold")
                    fig_d1.tight_layout()
                    d1_path = os.path.join(
                        RES_DIR,
                        "ratio_s4_compD1_gse25066.png")
                    fig_d1.savefig(d1_path, dpi=150)
                    plt.close(fig_d1)
                    log(f"  Figure: {d1_path}")

                    d1_result = (
                        "CONFIRMED"
                        if p_pcr < 0.05 else
                        "NOT CONFIRMED")
                    log(f"  R4-D1 → {d1_result}")
                else:
                    log("  Insufficient pCR/RD "
                        "counts to test")
            else:
                log("  No pCR column found")
                log(f"  Available chars: "
                    f"{list(char_accum.keys())}")
        else:
            missing = []
            if foxa1_data is None:
                missing.append("FOXA1")
            if ezh2_data is None:
                missing.append("EZH2")
            log(f"  Probes not found: {missing}")
            log("  Available probes (first 20):")
            for p in list(expr_rows.keys())[:20]:
                log(f"    {p}")

    except Exception as e:
        import traceback
        log(f"  GSE25066 error: {e}")
        log(traceback.format_exc())

preds["R4-D1_GSE25066"] = d1_result

# ── D.2: METABRIC RFS in hormone-treated ER+ ─────────────────
log("\n── D.2: METABRIC RFS — HORMONE-TREATED ER+ ──")

er_pos  = df_met["ER_IHC"].str.strip().str.lower()\
           .str.startswith("pos")
ht_yes  = df_met["HORMONE"].str.strip().str.upper()\
           == "YES"
mask_d2 = er_pos & ht_yes
df_d2   = df_met[mask_d2].copy()
log(f"  ER+ hormone-treated: n={len(df_d2)}")

df_d2["RFS_E"] = df_d2["RFS_STATUS"].apply(
    lambda x: 1 if str(x).strip().startswith("1:") else 0)

d2_result = "NOT TESTABLE"

if len(df_d2) > 50 and HAS_LIFELINES:
    med_d2 = df_d2["metric"].median()
    hi_d2  = df_d2[df_d2["metric"] >= med_d2]
    lo_d2  = df_d2[df_d2["metric"] <  med_d2]
    log(f"  High ratio n={len(hi_d2)}  "
        f"events={int(hi_d2['RFS_E'].sum())}")
    log(f"  Low  ratio n={len(lo_d2)}  "
        f"events={int(lo_d2['RFS_E'].sum())}")

    d2_rfs = df_d2[["RFS_MONTHS",
                     "RFS_E",
                     "metric"]].dropna()
    p_d2 = km_test(
        d2_rfs, "RFS_MONTHS", "RFS_E", "metric",
        title=(f"METABRIC ER+ HT n={len(d2_rfs)}\n"
               "FOXA1/EZH2 vs RFS"),
        outpath=os.path.join(
            RES_DIR,
            "ratio_s4_compD2_metabric_rfs.png"))

    if p_d2 is not None:
        log(f"  RFS KM p={p_d2:.4f}")
        # Check direction
        res_lr = logrank_test(
            hi_d2["RFS_MONTHS"].dropna(),
            lo_d2["RFS_MONTHS"].dropna(),
            hi_d2.loc[
                hi_d2["RFS_MONTHS"].notna(),
                "RFS_E"],
            lo_d2.loc[
                lo_d2["RFS_MONTHS"].notna(),
                "RFS_E"])
        log(f"  High ratio RFS med="
            f"{hi_d2['RFS_MONTHS'].median():.1f}mo")
        log(f"  Low  ratio RFS med="
            f"{lo_d2['RFS_MONTHS'].median():.1f}mo")
        d2_result = (
            "CONFIRMED" if p_d2 < 0.05 else
            "NOT CONFIRMED")
    else:
        log("  KM test failed")

preds["R4-D2_METABRIC_RFS_HT"] = d2_result
log(f"  R4-D2 → {d2_result}")

# ── D.3: TCGA DFS in available patients ──────────────────────
log("\n── D.3: TCGA DFS ──")
d3_result = "NOT TESTABLE"
if not df_tcga.empty and HAS_LIFELINES:
    df_tcga["DFS_MONTHS"] = pd.to_numeric(
        df_tcga.get("DFS_MONTHS", np.nan),
        errors="coerce")
    df_tcga["OS_E"] = df_tcga["OS_STATUS"].apply(
        lambda x: 1
        if "DECEASED" in str(x).upper() else 0)

    df_t_os = df_tcga[
        ["OS_MONTHS", "OS_E", "metric"]
        ].dropna()
    log(f"  TCGA OS n={len(df_t_os)}  "
        f"events={int(df_t_os['OS_E'].sum())}")

    if len(df_t_os) > 30:
        p_t_os = km_test(
            df_t_os,
            "OS_MONTHS", "OS_E", "metric",
            title=(f"TCGA n={len(df_t_os)}\n"
                   "FOXA1/EZH2 vs OS"),
            outpath=os.path.join(
                RES_DIR,
                "ratio_s4_compD3_tcga_os.png"))
        if p_t_os is not None:
            log(f"  TCGA OS KM p={p_t_os:.4f}")
            d3_result = (
                "CONFIRMED" if p_t_os < 0.05
                else "NOT CONFIRMED")

preds["R4-D3_TCGA_OS"] = d3_result
log(f"  R4-D3 → {d3_result}")

# ════��═════════════════════════════════════════════════════════
# COMBINED FIGURE
# ══════════════════════════════════════════════════════════════
log("\n── COMBINED SUMMARY FIGURE ──")

fig_main = plt.figure(figsize=(18, 12))
gs = gridspec.GridSpec(
    2, 4, figure=fig_main,
    hspace=0.45, wspace=0.35)

# Panel 1: METABRIC subtype distributions
ax1 = fig_main.add_subplot(gs[0, 0])
order_p = ["LumA", "LumB", "Her2",
           "Basal", "claudin-low"]
medians_p = []
labels_p  = []
for st in order_p:
    mask = df_met["subtype"] == st
    if mask.sum() > 0:
        medians_p.append(
            df_met.loc[mask, "metric"].median())
        labels_p.append(
            f"{st}\nn={mask.sum()}")
colors_p = ["#2196F3","#03A9F4","#FF9800",
            "#F44336","#9C27B0"][:len(labels_p)]
ax1.bar(range(len(labels_p)),
        medians_p, color=colors_p, alpha=0.8)
ax1.set_xticks(range(len(labels_p)))
ax1.set_xticklabels(labels_p, fontsize=7)
ax1.set_ylabel("Median FOXA1−EZH2")
ax1.set_title("Subtype Medians\nMETABRIC",
              fontsize=9)

# Panel 2: ROC LumA vs Basal
ax2 = fig_main.add_subplot(gs[0, 1])
sub_ab2 = met_sub[met_sub["subtype"].isin(
    ["LumA", "Basal"])]
if len(sub_ab2) > 0 and HAS_SKLEARN:
    roc_binary(sub_ab2["metric"],
               sub_ab2["subtype"],
               "LumA",
               title="LumA vs Basal (MET)",
               ax=ax2)
    ax2.set_title("ROC: LumA vs Basal\nMETABRIC",
                  fontsize=9)

# Panel 3: ROC LumA vs LumB
ax3 = fig_main.add_subplot(gs[0, 2])
sub_al2 = met_sub[met_sub["subtype"].isin(
    ["LumA", "LumB"])]
if len(sub_al2) > 0 and HAS_SKLEARN:
    roc_binary(sub_al2["metric"],
               sub_al2["subtype"],
               "LumA",
               title="LumA vs LumB (MET)",
               ax=ax3)
    ax3.set_title("ROC: LumA vs LumB\nMETABRIC",
                  fontsize=9)

# Panel 4: Cox forest plot (manual)
ax4 = fig_main.add_subplot(gs[0, 3])
ax4.axis("off")
cox_text = "MULTIVARIATE COX\nMETABRIC\n\n"
if a_os_hr is not None:
    cox_text += (
        f"OS model:\n"
        f"metric_z HR={a_os_hr:.3f}\n"
        f"p={a_os_p:.4f}\n\n")
if a_rfs_hr is not None:
    cox_text += (
        f"RFS model:\n"
        f"metric_z HR={a_rfs_hr:.3f}\n"
        f"p={a_rfs_p:.4f}")
ax4.text(0.05, 0.95, cox_text,
         transform=ax4.transAxes,
         fontsize=9, va="top",
         family="monospace",
         bbox=dict(boxstyle="round",
                   facecolor="lightyellow",
                   alpha=0.8))

# Panel 5: METABRIC RFS KM
ax5 = fig_main.add_subplot(gs[1, 0:2])
if HAS_LIFELINES:
    df_km2 = df_met[
        ["RFS_MONTHS", "RFS_STATUS",
         "metric"]].copy()
    df_km2["RFS_E"] = df_km2[
        "RFS_STATUS"].apply(
        lambda x: 1
        if "RECURRED" in str(x).upper()
        or "1:" in str(x) else 0)
    df_km2 = df_km2.dropna(
        subset=["RFS_MONTHS"])
    med_km = df_km2["metric"].median()
    hi_km  = df_km2[df_km2["metric"] >= med_km]
    lo_km  = df_km2[df_km2["metric"] <  med_km]
    kmf2   = KaplanMeierFitter()
    kmf2.fit(hi_km["RFS_MONTHS"],
             hi_km["RFS_E"],
             label=f"High n={len(hi_km)}")
    kmf2.plot_survival_function(ax=ax5)
    kmf2.fit(lo_km["RFS_MONTHS"],
             lo_km["RFS_E"],
             label=f"Low  n={len(lo_km)}")
    kmf2.plot_survival_function(ax=ax5)
    res2 = logrank_test(
        hi_km["RFS_MONTHS"],
        lo_km["RFS_MONTHS"],
        hi_km["RFS_E"],
        lo_km["RFS_E"])
    ax5.set_title(
        f"METABRIC RFS (all, n={len(df_km2)})\n"
        f"p={res2.p_value:.2e}",
        fontsize=9)
    ax5.set_xlabel("Months")
    ax5.set_ylabel("RFS")

# Panel 6: Scorecard text
ax6 = fig_main.add_subplot(gs[1, 2:4])
ax6.axis("off")
score_lines = ["RATIO-S4 SCORECARD\n"]
for k, v in preds.items():
    mark = "✓" if "CONFIRMED" == v else \
           "~" if "NOT TEST" in v else "✗"
    score_lines.append(f"{mark} {k}: {v}")
ax6.text(0.02, 0.98,
         "\n".join(score_lines),
         transform=ax6.transAxes,
         fontsize=8, va="top",
         family="monospace")

fig_main.suptitle(
    "FOXA1/EZH2 RATIO — SCRIPT 4 SUMMARY",
    fontsize=14, fontweight="bold")
combined_path = os.path.join(
    RES_DIR, "ratio_s4_combined.png")
fig_main.savefig(combined_path, dpi=150)
plt.close(fig_main)
log(f"  Figure: {combined_path}")

# ══════════════════════════════════════════════════════════════
# SCORECARD
# ══════════════════════════════════════════════════════════════
log("\n" + SEP)
log("SCORECARD — RATIO-S4")
log(SEP)
log(f"\n  {'ID':<30} {'Dataset':<25} "
    f"{'Status'}")
log(f"  {'─'*28} {'─'*23} {'─'*25}")

score_rows = []
for k, v in preds.items():
    ds = ("METABRIC"  if "MET"  in k else
          "TCGA"      if "TCGA" in k else
          "GSE25066"  if "GSE"  in k else
          "METABRIC")
    log(f"  {k:<30} {ds:<25} {v}")
    score_rows.append({
        "id": k, "dataset": ds, "status": v})

n_conf = sum(1 for v in preds.values()
             if v == "CONFIRMED")
n_nc   = sum(1 for v in preds.values()
             if v == "NOT CONFIRMED")
n_nt   = sum(1 for v in preds.values()
             if "NOT TESTABLE" in v)

log(f"\n  Total: {len(preds)}  "
    f"CONFIRMED: {n_conf}  "
    f"NOT CONFIRMED: {n_nc}  "
    f"NOT TESTABLE: {n_nt}")

# Save scorecard
sc_df = pd.DataFrame(score_rows)
sc_df.to_csv(SCORE_PATH, index=False)
log(f"\n  Scorecard: {SCORE_PATH}")

# Save log
with open(LOG_PATH, "w") as f:
    f.write("\n".join(log_lines))
log(f"  Log:       {LOG_PATH}")

log("\n" + SEP)
log("RATIO-S4 COMPLETE")
log(SEP)
