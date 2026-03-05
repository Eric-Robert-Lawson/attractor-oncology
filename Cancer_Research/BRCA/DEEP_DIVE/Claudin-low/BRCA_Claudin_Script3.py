"""
BRCA CLAUDIN-LOW — SCRIPT 3
OrganismCore — Document BRCA-S7g | 2026-03-05

ANALYSES AND PREDICTIONS FROM BRCA-S7f (before_script3.md):

ANALYSIS A — TREG/EFFECTOR IMMUNE RATIO
  S3-P1: FOXP3/CD8A ratio correlates with depth  r>0.20  [MOD]
  S3-P2: FOXP3/CD8A ratio predicts OS Stratum B  HR>1.0  [MOD-HIGH]
  S3-P3: Composite ratio strongest immune predictor       [MOD]

ANALYSIS B — CT ANTIGEN DEPTH CORRELATION AND SURVIVAL
  S3-P4: CT antigen composite correlates with depth r>0.25 [MOD]
  S3-P5: CT antigen score predicts OS in Stratum B HR>1.0  [MOD]
  S3-P6: CT antigen score correlates with Treg ratio r>0.15[MOD]

ANALYSIS C — LINEAGE MEMORY SUBGROUP ANALYSIS (Pommier 2020 proxy)
  S3-P7: Memory-low worse OS than memory-high (Stratum B)  [MOD-HIGH]
  S3-P8: Memory-low has higher CT antigen score             [MOD]
  S3-P9: Memory-low has deeper depth scores                 [HIGH]
  S3-P10: Memory-low has higher FOXP3/CD8A ratio            [MOD]

DERIVED FROM:
  Morel JCI 2017 — anti-PD1 paradox in claudin-low (Analysis A)
  Pommier Nature Comm 2020 — three subgroups (Analysis C)
  Axioms Prediction G — CT antigen as TYPE 4 marker (Analysis B)

SELF-CONTAINED:
  Reuses TCGA-BRCA expression + survival from Scripts 1–2.
  No new downloads required.
  Output to Claudin_Low_s3_results/.
"""

import os
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines import CoxPHFitter

# ============================================================
# CONFIGURATION
# ============================================================

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
S1_DATA_DIR  = os.path.join(SCRIPT_DIR,
                             "Claudin_Low_s1_analysis", "data")

BASE_DIR     = os.path.join(SCRIPT_DIR, "Claudin_Low_s3_results")
os.makedirs(BASE_DIR, exist_ok=True)

LOG_FILE         = os.path.join(BASE_DIR, "cl_s3_log.txt")
FIG_FILE         = os.path.join(BASE_DIR, "cl_s3_figure.png")
SCORECARD_FILE   = os.path.join(BASE_DIR, "cl_s3_scorecard.csv")
IMMUNE_FILE      = os.path.join(BASE_DIR, "cl_s3_immune_ratios.csv")
CT_FILE          = os.path.join(BASE_DIR, "cl_s3_ct_antigen.csv")
MEMORY_FILE      = os.path.join(BASE_DIR, "cl_s3_lineage_memory.csv")

EXPR_FILE        = os.path.join(S1_DATA_DIR, "TCGA_BRCA_HiSeqV2.gz")
CLIN_FILE        = os.path.join(S1_DATA_DIR, "TCGA_BRCA_clinicalMatrix.tsv")
SURV_FILE        = os.path.join(S1_DATA_DIR, "TCGA_pancan_survival.tsv")

# ── Classification constants (must match Scripts 1 and 2) ────
PATIENT_PREFIX_LEN         = 12
CL_PRIMARY_THRESHOLD       = 7
ERBB2_EXCLUSION_PERCENTILE = 90

CL_NEG_MARKERS = ["CLDN3", "CLDN4", "CLDN7", "CDH1", "ESR1"]
CL_POS_MARKERS = ["VIM", "CD44", "SNAI1", "ZEB1", "FN1"]
CL_SIGNATURE_GENES = CL_NEG_MARKERS + CL_POS_MARKERS

# ── Gene lists ────────────────────────────────────────────────
TREG_GENES     = ["FOXP3", "PDCD1", "TIGIT"]
EFFECTOR_GENES = ["CD8A", "GZMB", "PRF1"]
ALL_IMMUNE     = TREG_GENES + EFFECTOR_GENES + ["LAG3", "CD4",
                  "IFNG", "CD274", "CTLA4", "HAVCR2"]

CT_ANTIGEN_GENES = ["GAGE1", "GAGE2D", "GAGE4", "GAGE12D",
                    "GAGE12J", "CT45A3", "CT45A4",
                    "STRA8", "DPPA2"]

LINEAGE_MEMORY_GENES = ["FOXA1", "SPDEF", "GATA3"]

MIN_SURV_N = 10

# ============================================================
# LOGGING
# ============================================================

_log_lines = []

def log(msg=""):
    print(msg)
    _log_lines.append(str(msg))

def flush_log():
    with open(LOG_FILE, "w") as f:
        f.write("\n".join(_log_lines))

# ============================================================
# STEP 1 — LOAD DATA (reuse Script 1/2 files)
# ============================================================

def load_data():
    log("=" * 65)
    log("STEP 1: LOAD DATA")
    log("=" * 65)

    for f, label in [(EXPR_FILE, "Expression"),
                     (CLIN_FILE, "Clinical"),
                     (SURV_FILE, "Survival")]:
        if not os.path.exists(f):
            log(f"  FATAL: {label} missing: {f}")
            flush_log()
            sys.exit(1)
        log(f"  Found: {label}")

    # Expression
    try:
        expr = pd.read_csv(EXPR_FILE, sep="\t",
                           index_col=0, compression="gzip")
    except Exception:
        expr = pd.read_csv(EXPR_FILE, sep="\t", index_col=0)
    if expr.shape[0] < expr.shape[1]:
        expr = expr.T
    log(f"  Expression shape: {expr.shape}")

    # Clinical
    with open(CLIN_FILE, "rb") as f:
        first = f.readline().decode("utf-8", errors="replace")
    sep = "\t" if "\t" in first else ","
    clin = pd.read_csv(CLIN_FILE, sep=sep,
                       index_col=0, low_memory=False)
    clin.index = [str(s).replace(".", "-") for s in clin.index]

    # Survival
    surv = pd.read_csv(SURV_FILE, sep="\t", low_memory=False)
    id_col = None
    for cand in ["sample", "_PATIENT", "bcr_patient_barcode"]:
        if cand in surv.columns:
            id_col = cand
            break
    if id_col:
        surv = surv.set_index(id_col)
    surv.index = [str(s).replace(".", "-") for s in surv.index]

    return expr, clin, surv

# ============================================================
# STEP 2 — REBUILD POPULATIONS (identical to Scripts 1 and 2)
# ============================================================

def build_populations(expr, clin):
    log("")
    log("=" * 65)
    log("STEP 2: REBUILD POPULATIONS")
    log("=" * 65)

    samples = list(expr.columns)
    tumour_samples = []
    normal_samples = []
    for s in samples:
        parts = s.split("-")
        stype = parts[3][:2] if len(parts) >= 4 else "??"
        if stype == "01":
            tumour_samples.append(s)
        elif stype == "11":
            normal_samples.append(s)

    pam50_col = None
    for cand in ["PAM50Call_RNAseq", "PAM50_mRNA_nature2012"]:
        if cand in clin.columns:
            pam50_col = cand
            break
    pam50_map = {}
    if pam50_col:
        for pat_id in clin.index:
            pfx = str(pat_id)[:PATIENT_PREFIX_LEN]
            pam50_map[pfx] = str(clin.loc[pat_id, pam50_col])

    def get_pam50(s):
        return pam50_map.get(s[:PATIENT_PREFIX_LEN], "Unknown")

    subtype_buckets = {"LumA": [], "LumB": [],
                       "Her2": [], "Basal": [], "Normal": []}
    for s in tumour_samples:
        p = get_pam50(s)
        if p in subtype_buckets:
            subtype_buckets[p].append(s)

    tumour_expr = expr[tumour_samples]
    sig_present = [g for g in CL_SIGNATURE_GENES if g in expr.index]
    cohort_medians = {g: float(np.median(tumour_expr.loc[g].values))
                      for g in sig_present}

    scores = {}
    for s in tumour_samples:
        score = 0
        for g in CL_POS_MARKERS:
            if g in sig_present:
                if float(expr.loc[g, s]) > cohort_medians[g]:
                    score += 1
        for g in CL_NEG_MARKERS:
            if g in sig_present:
                if float(expr.loc[g, s]) < cohort_medians[g]:
                    score += 1
        scores[s] = score

    cl_raw = [s for s, sc in scores.items() if sc >= CL_PRIMARY_THRESHOLD]
    if "ERBB2" in expr.index:
        cutoff = np.percentile(
            tumour_expr.loc["ERBB2", tumour_samples].values,
            ERBB2_EXCLUSION_PERCENTILE)
        cl_samples = [s for s in cl_raw
                      if float(expr.loc["ERBB2", s]) <= cutoff]
    else:
        cl_samples = cl_raw

    esr1_med = float(np.median(tumour_expr.loc["ESR1"].values)) \
               if "ESR1" in expr.index else 0.0
    stratum_b = [s for s in cl_samples
                 if "ESR1" in expr.index and
                 float(expr.loc["ESR1", s]) < esr1_med]
    stratum_c = [s for s in cl_samples
                 if get_pam50(s) in ("Basal", "Normal")]

    populations = {
        "claudin_low": cl_samples,
        "stratum_b":   stratum_b,
        "stratum_c":   stratum_c,
        "normal":      normal_samples,
        "luma":        subtype_buckets["LumA"],
        "basal":       subtype_buckets["Basal"],
    }

    log(f"  claudin_low (Stratum A): n={len(cl_samples)}")
    log(f"  Stratum B (ESR1-low):    n={len(stratum_b)}")
    log(f"  Stratum C (Basal+Norm):  n={len(stratum_c)}")

    return populations

# ============================================================
# STEP 3 — COMPUTE DEPTH SCORES
# ============================================================

def compute_depth_scores(expr, populations):
    log("")
    log("=" * 65)
    log("STEP 3: COMPUTE DEPTH SCORES")
    log("=" * 65)

    all_tumour = (populations["claudin_low"] +
                  populations["luma"] +
                  populations["basal"])
    tumour_expr = expr[all_tumour]
    pos_genes = [g for g in CL_POS_MARKERS if g in expr.index]
    neg_genes = [g for g in CL_NEG_MARKERS if g in expr.index]

    def z_score(arr):
        s = arr.std()
        return (arr - arr.mean()) / s if s > 1e-9 else np.zeros_like(arr)

    pos_z = np.column_stack([
        z_score(tumour_expr.loc[g].values.astype(float))
        for g in pos_genes]) if pos_genes else np.zeros((len(all_tumour), 1))
    neg_z = np.column_stack([
        z_score(tumour_expr.loc[g].values.astype(float))
        for g in neg_genes]) if neg_genes else np.zeros((len(all_tumour), 1))

    depth_arr = np.mean(pos_z, axis=1) - np.mean(neg_z, axis=1)
    depth_scores = dict(zip(all_tumour, depth_arr))
    log(f"  Depth scores: n={len(depth_scores)}")
    return depth_scores

# ============================================================
# STEP 4 — BUILD SURVIVAL VECTORS
# ============================================================

def build_survival_vectors(expr, populations, surv):
    log("")
    log("=" * 65)
    log("STEP 4: BUILD SURVIVAL VECTORS")
    log("=" * 65)

    os_time, os_event = {}, {}
    surv_time_col, surv_event_col = None, None
    for tc in ["OS.time", "os_time", "OS_time"]:
        if tc in surv.columns:
            surv_time_col = tc
            break
    for ec in ["OS", "os_event", "OS_event"]:
        if ec in surv.columns:
            surv_event_col = ec
            break

    if surv_time_col and surv_event_col:
        all_samples = list(expr.columns)
        surv_pfx = {str(idx)[:PATIENT_PREFIX_LEN]: idx
                    for idx in surv.index}
        for s in all_samples:
            pfx = s[:PATIENT_PREFIX_LEN]
            if s in surv.index:
                idx = s
            elif pfx in surv_pfx:
                idx = surv_pfx[pfx]
            else:
                continue
            try:
                t = float(surv.loc[idx, surv_time_col])
                e = float(surv.loc[idx, surv_event_col])
                if not np.isnan(t) and t > 0:
                    os_time[s]  = t
                    os_event[s] = e
            except (ValueError, TypeError):
                pass

    log(f"  Samples with survival: {len(os_time)}")
    cl = populations["claudin_low"]
    log(f"  Claudin-low covered: "
        f"{sum(1 for s in cl if s in os_time)}/{len(cl)}")
    return os_time, os_event

# ============================================================
# SURVIVAL UTILITIES
# ============================================================

def fmt_p(p):
    if p is None or np.isnan(p):
        return "p=nan"
    if p < 1e-5: return f"p={p:.2e} **"
    if p < 0.05: return f"p={p:.3f} *"
    if p < 0.10: return f"p={p:.3f} (trend)"
    return              f"p={p:.3f}"

def get_surv_arrays(samples, os_time, os_event):
    T, E, valid = [], [], []
    for s in samples:
        if s in os_time and s in os_event:
            t, e = os_time[s], os_event[s]
            if not np.isnan(t) and not np.isnan(e) and t > 0:
                T.append(t); E.append(int(e)); valid.append(s)
    return np.array(T), np.array(E), valid

def run_logrank(T_h, E_h, T_l, E_l, label=""):
    if len(T_h) < MIN_SURV_N or len(T_l) < MIN_SURV_N:
        log(f"    {label}: n insufficient — skipped")
        return np.nan, np.nan
    try:
        res = logrank_test(T_h, T_l,
                           event_observed_A=E_h,
                           event_observed_B=E_l)
        p = res.p_value
        rh = E_h.sum() / T_h.sum() if T_h.sum() > 0 else np.nan
        rl = E_l.sum() / T_l.sum() if T_l.sum() > 0 else np.nan
        hr = rh / rl if rl and rl > 0 else np.nan
        log(f"    {label}: n_h={len(T_h)} n_l={len(T_l)} "
            f"ev_h={E_h.sum()} ev_l={E_l.sum()}  "
            f"{fmt_p(p)}  HR≈{hr:.3f}")
        return p, hr
    except Exception as e:
        log(f"    {label}: failed — {e}")
        return np.nan, np.nan

def median_split(samples, score_dict, os_time, os_event):
    """Median split on score_dict within samples."""
    scored = [(s, score_dict[s]) for s in samples if s in score_dict]
    if len(scored) < MIN_SURV_N * 2:
        return (np.array([]),) * 4, [], []
    med = np.median([v for _, v in scored])
    high_s = [s for s, v in scored if v >= med]
    low_s  = [s for s, v in scored if v <  med]
    T_h, E_h, vh = get_surv_arrays(high_s, os_time, os_event)
    T_l, E_l, vl = get_surv_arrays(low_s,  os_time, os_event)
    return T_h, E_h, T_l, E_l, vh, vl

def tertile_split(samples, score_dict, os_time, os_event):
    scored = [(s, score_dict[s]) for s in samples if s in score_dict]
    if len(scored) < MIN_SURV_N * 3:
        return (np.array([]),) * 6, [], [], []
    scored.sort(key=lambda x: x[1])
    n = len(scored); t1 = n // 3; t2 = 2 * n // 3
    low_s  = [s for s, _ in scored[:t1]]
    mid_s  = [s for s, _ in scored[t1:t2]]
    high_s = [s for s, _ in scored[t2:]]
    T_h, E_h, vh = get_surv_arrays(high_s, os_time, os_event)
    T_m, E_m, vm = get_surv_arrays(mid_s,  os_time, os_event)
    T_l, E_l, vl = get_surv_arrays(low_s,  os_time, os_event)
    return T_h, E_h, T_l, E_l, T_m, E_m, vh, vm, vl

def safe_corr(x, y, label=""):
    valid = [(xi, yi) for xi, yi in zip(x, y)
             if not np.isnan(xi) and not np.isnan(yi)]
    if len(valid) < 10:
        log(f"    {label}: n<10 for correlation")
        return np.nan, np.nan
    xv, yv = zip(*valid)
    r, p = stats.pearsonr(xv, yv)
    rs, ps = stats.spearmanr(xv, yv)
    log(f"    {label}: r={r:+.4f} {fmt_p(p)}  "
        f"(Spearman rs={rs:+.4f} {fmt_p(ps)})  n={len(valid)}")
    return r, p

# ============================================================
# STEP 5 — COMPUTE SCORE DICTIONARIES
# ============================================================

def compute_all_scores(expr, populations, depth_scores):
    log("")
    log("=" * 65)
    log("STEP 5: COMPUTE SCORE DICTIONARIES")
    log("=" * 65)

    cl_a = populations["claudin_low"]

    # Gene availability
    log("\n  Gene availability check:")
    for g in (TREG_GENES + EFFECTOR_GENES + CT_ANTIGEN_GENES +
              LINEAGE_MEMORY_GENES):
        log(f"    {g:<12} {'FOUND' if g in expr.index else 'MISSING'}")

    # ── Treg and effector raw values ──────────────────────────
    def safe_val(g, s):
        return float(expr.loc[g, s]) if g in expr.index else np.nan

    foxp3  = {s: safe_val("FOXP3", s)  for s in cl_a}
    cd8a   = {s: safe_val("CD8A", s)   for s in cl_a}
    gzmb   = {s: safe_val("GZMB", s)   for s in cl_a}
    prf1   = {s: safe_val("PRF1", s)   for s in cl_a}
    tigit  = {s: safe_val("TIGIT", s)  for s in cl_a}
    pdcd1  = {s: safe_val("PDCD1", s)  for s in cl_a}

    # ── FOXP3/CD8A ratio ──────────────────────────────────────
    eps = 1e-3
    foxp3_cd8a = {
        s: foxp3[s] / (cd8a[s] + eps) if not np.isnan(foxp3[s])
        and not np.isnan(cd8a[s]) else np.nan
        for s in cl_a}

    # ── FOXP3/GZMB ratio ─────────────────────────────────────
    foxp3_gzmb = {
        s: foxp3[s] / (gzmb[s] + eps) if not np.isnan(foxp3[s])
        and not np.isnan(gzmb[s]) else np.nan
        for s in cl_a}

    # ── TIGIT/CD8A ratio ──────────────────────────────────────
    tigit_cd8a = {
        s: tigit[s] / (cd8a[s] + eps) if not np.isnan(tigit[s])
        and not np.isnan(cd8a[s]) else np.nan
        for s in cl_a}

    # ── Composite Treg:effector ratio ─────────────────────────
    # (FOXP3 + TIGIT) / (CD8A + GZMB + PRF1)
    composite_treg = {}
    for s in cl_a:
        treg_num = sum([v for v in [foxp3.get(s, np.nan),
                                     tigit.get(s, np.nan)]
                        if not np.isnan(v)])
        eff_den  = sum([v for v in [cd8a.get(s, np.nan),
                                     gzmb.get(s, np.nan),
                                     prf1.get(s, np.nan)]
                        if not np.isnan(v)])
        composite_treg[s] = treg_num / (eff_den + eps) \
                             if eff_den > 0 else np.nan

    log(f"\n  Treg/effector scores computed: n={len(foxp3_cd8a)}")

    # ── CT antigen composite score ────────────────────────────
    ct_avail = [g for g in CT_ANTIGEN_GENES if g in expr.index]
    log(f"\n  CT antigen genes available: {ct_avail}")

    ct_z_parts = []
    ct_raw = {}
    for g in ct_avail:
        vals = np.array([safe_val(g, s) for s in cl_a])
        valid = vals[~np.isnan(vals)]
        if len(valid) < 10:
            continue
        sd = valid.std()
        if sd < 1e-9:
            continue
        z = (vals - valid.mean()) / sd
        ct_z_parts.append(z)
        for s, v in zip(cl_a, vals):
            ct_raw.setdefault(s, {})[g] = v

    if ct_z_parts:
        ct_composite_arr = np.mean(ct_z_parts, axis=0)
        ct_composite = dict(zip(cl_a, ct_composite_arr))
        log(f"  CT composite score: n={len(ct_composite)} genes={len(ct_z_parts)}")
    else:
        ct_composite = {}
        log("  WARNING: No CT antigen genes available for composite score")

    # ── Lineage memory score ───────────────────────────────────
    mem_avail = [g for g in LINEAGE_MEMORY_GENES if g in expr.index]
    log(f"\n  Lineage memory genes available: {mem_avail}")

    mem_z_parts = []
    for g in mem_avail:
        vals = np.array([safe_val(g, s) for s in cl_a])
        valid = vals[~np.isnan(vals)]
        if len(valid) < 10:
            continue
        sd = valid.std()
        if sd < 1e-9:
            continue
        mem_z_parts.append((vals - valid.mean()) / sd)

    if mem_z_parts:
        mem_arr = np.mean(mem_z_parts, axis=0)
        lineage_memory = dict(zip(cl_a, mem_arr))
        log(f"  Lineage memory scores: n={len(lineage_memory)}")
    else:
        lineage_memory = {}
        log("  WARNING: No lineage memory genes available")

    all_scores = {
        "foxp3_cd8a":    foxp3_cd8a,
        "foxp3_gzmb":    foxp3_gzmb,
        "tigit_cd8a":    tigit_cd8a,
        "composite_treg": composite_treg,
        "ct_composite":  ct_composite,
        "lineage_memory": lineage_memory,
        # Raw single genes for comparison
        "FOXP3":  foxp3,
        "CD8A":   cd8a,
        "GZMB":   gzmb,
        "PRF1":   prf1,
        "TIGIT":  tigit,
        "PDCD1":  pdcd1,
    }

    # Per-CT-antigen raw scores for correlation
    for g in ct_avail:
        all_scores[g] = {s: safe_val(g, s) for s in cl_a}

    return all_scores

# ============================================================
# STEP 6 — ANALYSIS A: TREG/EFFECTOR RATIO
# ============================================================

def analysis_a(populations, depth_scores, all_scores,
               os_time, os_event, scorecard):
    log("")
    log("=" * 65)
    log("ANALYSIS A — TREG/EFFECTOR IMMUNE RATIO")
    log("  Source: Morel JCI 2017")
    log("=" * 65)

    cl_a = populations["claudin_low"]
    cl_b = populations["stratum_b"]

    def record(pid, gene, stratum, n_h, n_l, p, hr,
               dir_ok, status, note=""):
        scorecard.append({
            "analysis": "A",
            "prediction": pid, "gene": gene,
            "stratum": stratum, "n_high": n_h,
            "n_low": n_l, "p_value": p,
            "HR_approx": hr, "direction_correct": dir_ok,
            "confirmed": status, "note": note,
        })

    # S3-P1: FOXP3/CD8A ratio vs depth score
    log("")
    log("  ── S3-P1: FOXP3/CD8A ratio vs depth score (Stratum A) ──")
    d_vals = [depth_scores[s] for s in cl_a if s in depth_scores
              and s in all_scores["foxp3_cd8a"]
              and not np.isnan(all_scores["foxp3_cd8a"][s])]
    r_vals = [all_scores["foxp3_cd8a"][s] for s in cl_a
              if s in depth_scores
              and s in all_scores["foxp3_cd8a"]
              and not np.isnan(all_scores["foxp3_cd8a"][s])]
    r1, p1 = safe_corr(d_vals, r_vals, "FOXP3/CD8A vs depth")
    p1_conf = ("CONFIRMED" if not np.isnan(r1) and r1 > 0.20
               else "NOT CONFIRMED")
    log(f"  Status: [{p1_conf}]")
    record("S3-P1", "FOXP3/CD8A", "A",
           len(d_vals), 0, p1, r1,
           (not np.isnan(r1) and r1 > 0.20),
           p1_conf, "correlation not survival")

    # Additional correlations for context
    log("\n  Additional Treg/effector correlations with depth:")
    for name in ["foxp3_gzmb", "tigit_cd8a", "composite_treg",
                 "FOXP3", "CD8A", "GZMB", "TIGIT"]:
        if name not in all_scores:
            continue
        dv = [depth_scores[s] for s in cl_a
              if s in depth_scores and s in all_scores[name]
              and not np.isnan(all_scores[name][s])]
        rv = [all_scores[name][s] for s in cl_a
              if s in depth_scores and s in all_scores[name]
              and not np.isnan(all_scores[name][s])]
        safe_corr(dv, rv, name)

    # S3-P2: FOXP3/CD8A ratio vs OS, Stratum B
    log("")
    log("  ── S3-P2: FOXP3/CD8A ratio vs OS (Stratum B) ──")
    res2 = tertile_split(cl_b, all_scores["foxp3_cd8a"],
                         os_time, os_event)
    T2_h, E2_h, T2_l, E2_l = res2[0], res2[1], res2[2], res2[3]
    p2, hr2 = run_logrank(T2_h, E2_h, T2_l, E2_l,
                          "S3-P2 [FOXP3/CD8A tertile, Stratum B]")
    dir2 = not np.isnan(hr2) and hr2 > 1.0
    conf2 = ("CONFIRMED" if not np.isnan(p2) and p2 < 0.15 and dir2
             else "NOT CONFIRMED")
    log(f"  Status: [{conf2}]")
    record("S3-P2", "FOXP3/CD8A", "B",
           len(T2_h), len(T2_l), p2, hr2, dir2, conf2)

    # S3-P3: Composite ratio vs single genes
    log("")
    log("  ── S3-P3: Composite Treg:effector ratio vs OS (Stratum B) ──")
    log("  Testing all immune metrics for comparison:")

    best_p  = 1.0
    best_hr = np.nan
    best_lbl = ""
    for name, label in [
        ("foxp3_cd8a",     "FOXP3/CD8A"),
        ("foxp3_gzmb",     "FOXP3/GZMB"),
        ("tigit_cd8a",     "TIGIT/CD8A"),
        ("composite_treg", "Composite Treg:Eff"),
        ("FOXP3",          "FOXP3"),
        ("CD8A",           "CD8A"),
        ("GZMB",           "GZMB"),
        ("TIGIT",          "TIGIT"),
    ]:
        if name not in all_scores:
            continue
        res = tertile_split(cl_b, all_scores[name], os_time, os_event)
        Tx_h, Ex_h, Tx_l, Ex_l = res[0], res[1], res[2], res[3]
        px, hrx = run_logrank(Tx_h, Ex_h, Tx_l, Ex_l, label)
        if not np.isnan(px) and px < best_p:
            best_p = px; best_hr = hrx; best_lbl = label

    log(f"\n  Best single immune predictor: {best_lbl} p={best_p:.3f} HR={best_hr:.3f}")

    comp_res = tertile_split(cl_b, all_scores["composite_treg"],
                             os_time, os_event)
    Tc_h, Ec_h, Tc_l, Ec_l = comp_res[0], comp_res[1], comp_res[2], comp_res[3]
    pc, hrc = run_logrank(Tc_h, Ec_h, Tc_l, Ec_l,
                          "S3-P3 [Composite tertile, Stratum B]")
    is_best = (not np.isnan(pc) and
               not np.isnan(best_p) and
               pc <= best_p + 0.02)
    conf3 = "CONFIRMED" if is_best else "NOT CONFIRMED"
    log(f"  Composite is best predictor: {is_best}")
    log(f"  Status: [{conf3}]")
    record("S3-P3", "composite_treg", "B",
           len(Tc_h), len(Tc_l), pc, hrc,
           is_best, conf3,
           f"best_other={best_lbl} p={best_p:.3f}")

    return scorecard

# ============================================================
# STEP 7 — ANALYSIS B: CT ANTIGEN
# ============================================================

def analysis_b(populations, depth_scores, all_scores,
               os_time, os_event, scorecard):
    log("")
    log("=" * 65)
    log("ANALYSIS B — CT ANTIGEN DEPTH CORRELATION AND SURVIVAL")
    log("  Source: Axioms Prediction G")
    log("=" * 65)

    cl_a = populations["claudin_low"]
    cl_b = populations["stratum_b"]

    def record(pid, gene, stratum, n_h, n_l, p, hr,
               dir_ok, status, note=""):
        scorecard.append({
            "analysis": "B",
            "prediction": pid, "gene": gene,
            "stratum": stratum, "n_high": n_h,
            "n_low": n_l, "p_value": p,
            "HR_approx": hr, "direction_correct": dir_ok,
            "confirmed": status, "note": note,
        })

    ct_avail = [g for g in CT_ANTIGEN_GENES if g in all_scores]

    # S3-P4: CT antigen genes vs depth
    log("")
    log("  ── S3-P4: CT antigen genes vs depth score ──")
    log("  Individual CT antigen correlations with depth:")
    ct_depth_rs = {}
    for g in ct_avail:
        dv = [depth_scores[s] for s in cl_a
              if s in depth_scores and s in all_scores[g]
              and not np.isnan(all_scores[g][s])]
        rv = [all_scores[g][s] for s in cl_a
              if s in depth_scores and s in all_scores[g]
              and not np.isnan(all_scores[g][s])]
        r, p = safe_corr(dv, rv, g)
        ct_depth_rs[g] = r

    # Composite CT vs depth
    if "ct_composite" in all_scores and all_scores["ct_composite"]:
        dv = [depth_scores[s] for s in cl_a
              if s in depth_scores and s in all_scores["ct_composite"]
              and not np.isnan(all_scores["ct_composite"][s])]
        rv = [all_scores["ct_composite"][s] for s in cl_a
              if s in depth_scores and s in all_scores["ct_composite"]
              and not np.isnan(all_scores["ct_composite"][s])]
        r4, p4 = safe_corr(dv, rv, "CT composite vs depth")
        n_above_020 = sum(1 for r in ct_depth_rs.values()
                          if not np.isnan(r) and r > 0.20)
        conf4 = ("CONFIRMED"
                 if not np.isnan(r4) and r4 > 0.25 and n_above_020 >= 3
                 else "NOT CONFIRMED")
        log(f"  CT genes r>0.20: {n_above_020}/{len(ct_avail)}")
        log(f"  Status: [{conf4}]")
        record("S3-P4", "ct_composite", "A",
               len(dv), 0, p4, r4,
               (not np.isnan(r4) and r4 > 0.25),
               conf4,
               f"genes_above_0.20={n_above_020}")
    else:
        log("  CT composite not available — skipped")
        record("S3-P4", "ct_composite", "A",
               0, 0, np.nan, np.nan, False, "DATA MISSING")

    # S3-P5: CT composite vs OS, Stratum B
    log("")
    log("  ── S3-P5: CT composite score vs OS (Stratum B) ──")
    if "ct_composite" in all_scores and all_scores["ct_composite"]:
        res5 = tertile_split(cl_b, all_scores["ct_composite"],
                             os_time, os_event)
        T5_h, E5_h, T5_l, E5_l = res5[0], res5[1], res5[2], res5[3]
        p5, hr5 = run_logrank(T5_h, E5_h, T5_l, E5_l,
                              "S3-P5 [CT composite tertile, Stratum B]")
        dir5 = not np.isnan(hr5) and hr5 > 1.0
        conf5 = ("CONFIRMED" if not np.isnan(p5) and p5 < 0.15 and dir5
                 else "NOT CONFIRMED")
        log(f"  Status: [{conf5}]")
        record("S3-P5", "ct_composite", "B",
               len(T5_h), len(T5_l), p5, hr5, dir5, conf5)
    else:
        log("  CT composite not available")
        record("S3-P5", "ct_composite", "B",
               0, 0, np.nan, np.nan, False, "DATA MISSING")

    # S3-P6: CT composite vs FOXP3/CD8A ratio
    log("")
    log("  ── S3-P6: CT composite vs FOXP3/CD8A ratio ──")
    if ("ct_composite" in all_scores and all_scores["ct_composite"]
            and "foxp3_cd8a" in all_scores):
        ct_v  = [all_scores["ct_composite"][s] for s in cl_a
                 if s in all_scores["ct_composite"]
                 and s in all_scores["foxp3_cd8a"]
                 and not np.isnan(all_scores["ct_composite"][s])
                 and not np.isnan(all_scores["foxp3_cd8a"][s])]
        rat_v = [all_scores["foxp3_cd8a"][s] for s in cl_a
                 if s in all_scores["ct_composite"]
                 and s in all_scores["foxp3_cd8a"]
                 and not np.isnan(all_scores["ct_composite"][s])
                 and not np.isnan(all_scores["foxp3_cd8a"][s])]
        r6, p6 = safe_corr(ct_v, rat_v, "CT composite vs FOXP3/CD8A")
        dir6 = not np.isnan(r6) and r6 > 0.15
        conf6 = "CONFIRMED" if dir6 else "NOT CONFIRMED"
        log(f"  Status: [{conf6}]")
        record("S3-P6", "ct_vs_treg_ratio", "A",
               len(ct_v), 0, p6, r6, dir6, conf6)
    else:
        record("S3-P6", "ct_vs_treg_ratio", "A",
               0, 0, np.nan, np.nan, False, "DATA MISSING")

    return scorecard

# ============================================================
# STEP 8 — ANALYSIS C: LINEAGE MEMORY
# ============================================================

def analysis_c(populations, depth_scores, all_scores,
               os_time, os_event, scorecard):
    log("")
    log("=" * 65)
    log("ANALYSIS C — LINEAGE MEMORY SUBGROUP ANALYSIS")
    log("  Source: Pommier Nature Comm 2020")
    log("=" * 65)

    cl_a = populations["claudin_low"]
    cl_b = populations["stratum_b"]

    def record(pid, gene, stratum, n_h, n_l, p, hr,
               dir_ok, status, note=""):
        scorecard.append({
            "analysis": "C",
            "prediction": pid, "gene": gene,
            "stratum": stratum, "n_high": n_h,
            "n_low": n_l, "p_value": p,
            "HR_approx": hr, "direction_correct": dir_ok,
            "confirmed": status, "note": note,
        })

    if not all_scores.get("lineage_memory"):
        log("  Lineage memory scores not available — Analysis C skipped")
        for pid in ["S3-P7","S3-P8","S3-P9","S3-P10"]:
            record(pid, "lineage_memory", "B",
                   0, 0, np.nan, np.nan, False, "DATA MISSING")
        return scorecard

    mem = all_scores["lineage_memory"]

    # Memory-high vs memory-low split (median split)
    mem_vals = [mem[s] for s in cl_b if s in mem and not np.isnan(mem[s])]
    if len(mem_vals) < MIN_SURV_N * 2:
        log(f"  Insufficient Stratum B samples with memory score "
            f"(n={len(mem_vals)})")
        return scorecard

    mem_med = np.median(mem_vals)
    mem_high_b = [s for s in cl_b if s in mem and
                  not np.isnan(mem[s]) and mem[s] >= mem_med]
    mem_low_b  = [s for s in cl_b if s in mem and
                  not np.isnan(mem[s]) and mem[s] <  mem_med]

    log(f"\n  Lineage memory split (Stratum B):")
    log(f"    Memory-high (≥ median): n={len(mem_high_b)}")
    log(f"    Memory-low  (< median): n={len(mem_low_b)}")
    log(f"    Memory score median = {mem_med:.4f}")

    # S3-P7: Memory-low vs memory-high OS
    log("")
    log("  ── S3-P7: Memory-low vs memory-high OS (Stratum B) ──")
    T_mh, E_mh, _ = get_surv_arrays(mem_high_b, os_time, os_event)
    T_ml, E_ml, _ = get_surv_arrays(mem_low_b,  os_time, os_event)
    p7, hr7 = run_logrank(T_ml, E_ml, T_mh, E_mh,
                          "S3-P7 [memory-low vs memory-high]")
    # HR here = memory-low / memory-high. Prediction: memory-low worse → HR > 1.0
    dir7 = not np.isnan(hr7) and hr7 > 1.0
    conf7 = ("CONFIRMED" if not np.isnan(p7) and p7 < 0.15 and dir7
             else "NOT CONFIRMED")
    log(f"  Status: [{conf7}]")
    record("S3-P7", "lineage_memory_OS", "B",
           len(T_ml), len(T_mh), p7, hr7, dir7, conf7,
           "memory-low in n_high position")

    # S3-P8: CT antigen in memory-high vs memory-low
    log("")
    log("  ── S3-P8: CT antigen score in memory groups ──")
    if all_scores.get("ct_composite"):
        ct = all_scores["ct_composite"]
        ct_mh = [ct[s] for s in mem_high_b if s in ct and
                 not np.isnan(ct[s])]
        ct_ml = [ct[s] for s in mem_low_b  if s in ct and
                 not np.isnan(ct[s])]
        if len(ct_mh) > 5 and len(ct_ml) > 5:
            st8, pt8 = stats.mannwhitneyu(ct_ml, ct_mh,
                                           alternative="greater")
            mean_ml = np.mean(ct_ml)
            mean_mh = np.mean(ct_mh)
            log(f"    CT score — memory-low mean={mean_ml:.4f}  "
                f"memory-high mean={mean_mh:.4f}")
            log(f"    Mann-Whitney (memory-low > memory-high): {fmt_p(pt8)}")
            dir8 = mean_ml > mean_mh
            conf8 = ("CONFIRMED" if pt8 < 0.10 and dir8
                     else "NOT CONFIRMED")
            log(f"  Status: [{conf8}]")
            record("S3-P8", "ct_vs_memory", "B",
                   len(ct_ml), len(ct_mh), pt8,
                   mean_ml / mean_mh if mean_mh > 0 else np.nan,
                   dir8, conf8)
        else:
            log("  Insufficient CT data in memory groups")
            record("S3-P8", "ct_vs_memory", "B",
                   0, 0, np.nan, np.nan, False, "DATA MISSING")
    else:
        record("S3-P8", "ct_vs_memory", "B",
               0, 0, np.nan, np.nan, False, "CT MISSING")

    # S3-P9: Depth score in memory groups
    log("")
    log("  ── S3-P9: Depth score in memory groups ──")
    depth_mh = [depth_scores[s] for s in mem_high_b
                if s in depth_scores]
    depth_ml = [depth_scores[s] for s in mem_low_b
                if s in depth_scores]
    if len(depth_mh) > 5 and len(depth_ml) > 5:
        st9, pt9 = stats.mannwhitneyu(depth_ml, depth_mh,
                                       alternative="greater")
        mean_dml = np.mean(depth_ml)
        mean_dmh = np.mean(depth_mh)
        log(f"    Depth — memory-low mean={mean_dml:.4f}  "
            f"memory-high mean={mean_dmh:.4f}")
        log(f"    Mann-Whitney (memory-low deeper): {fmt_p(pt9)}")
        dir9 = mean_dml > mean_dmh
        conf9 = ("CONFIRMED" if pt9 < 0.05 and dir9
                 else "NOT CONFIRMED")
        log(f"  Status: [{conf9}]")
        record("S3-P9", "depth_vs_memory", "B",
               len(depth_ml), len(depth_mh), pt9,
               mean_dml / mean_dmh if mean_dmh != 0 else np.nan,
               dir9, conf9)
    else:
        record("S3-P9", "depth_vs_memory", "B",
               0, 0, np.nan, np.nan, False, "DATA MISSING")

    # S3-P10: FOXP3/CD8A ratio in memory groups
    log("")
    log("  ── S3-P10: FOXP3/CD8A ratio in memory groups ──")
    ratio = all_scores["foxp3_cd8a"]
    rat_mh = [ratio[s] for s in mem_high_b if s in ratio
              and not np.isnan(ratio[s])]
    rat_ml = [ratio[s] for s in mem_low_b  if s in ratio
              and not np.isnan(ratio[s])]
    if len(rat_mh) > 5 and len(rat_ml) > 5:
        st10, pt10 = stats.mannwhitneyu(rat_ml, rat_mh,
                                         alternative="greater")
        mean_rml = np.mean(rat_ml)
        mean_rmh = np.mean(rat_mh)
        log(f"    FOXP3/CD8A — memory-low mean={mean_rml:.4f}  "
            f"memory-high mean={mean_rmh:.4f}")
        log(f"    Mann-Whitney (memory-low higher): {fmt_p(pt10)}")
        dir10 = mean_rml > mean_rmh
        conf10 = ("CONFIRMED" if pt10 < 0.10 and dir10
                  else "NOT CONFIRMED")
        log(f"  Status: [{conf10}]")
        record("S3-P10", "treg_ratio_vs_memory", "B",
               len(rat_ml), len(rat_mh), pt10,
               mean_rml / mean_rmh if mean_rmh > 0 else np.nan,
               dir10, conf10)
    else:
        record("S3-P10", "treg_ratio_vs_memory", "B",
               0, 0, np.nan, np.nan, False, "DATA MISSING")

    return scorecard

# ============================================================
# STEP 9 — SAVE INTERMEDIATE CSV FILES
# ============================================================

def save_intermediate_files(populations, depth_scores,
                             all_scores, os_time, os_event):
    cl_a = populations["claudin_low"]

    # Immune ratios
    rows = []
    for s in cl_a:
        row = {"sample": s,
               "depth_score": depth_scores.get(s, np.nan),
               "os_time":  os_time.get(s, np.nan),
               "os_event": os_event.get(s, np.nan)}
        for k in ["foxp3_cd8a", "foxp3_gzmb", "tigit_cd8a",
                  "composite_treg", "FOXP3", "CD8A", "GZMB",
                  "PRF1", "TIGIT", "PDCD1"]:
            row[k] = all_scores.get(k, {}).get(s, np.nan)
        rows.append(row)
    pd.DataFrame(rows).to_csv(IMMUNE_FILE, index=False)

    # CT antigen
    ct_rows = []
    for s in cl_a:
        row = {"sample": s,
               "depth_score": depth_scores.get(s, np.nan),
               "ct_composite": all_scores.get("ct_composite",
                                               {}).get(s, np.nan)}
        for g in CT_ANTIGEN_GENES:
            row[g] = all_scores.get(g, {}).get(s, np.nan)
        ct_rows.append(row)
    pd.DataFrame(ct_rows).to_csv(CT_FILE, index=False)

    # Lineage memory
    mem_rows = []
    for s in cl_a:
        row = {"sample": s,
               "depth_score": depth_scores.get(s, np.nan),
               "lineage_memory": all_scores.get("lineage_memory",
                                                 {}).get(s, np.nan),
               "ct_composite": all_scores.get("ct_composite",
                                               {}).get(s, np.nan),
               "foxp3_cd8a": all_scores.get("foxp3_cd8a",
                                             {}).get(s, np.nan)}
        mem_rows.append(row)
    pd.DataFrame(mem_rows).to_csv(MEMORY_FILE, index=False)

    log(f"\n  Intermediate CSVs saved:")
    log(f"    {IMMUNE_FILE}")
    log(f"    {CT_FILE}")
    log(f"    {MEMORY_FILE}")

# ============================================================
# STEP 10 — FIGURE
# ============================================================

def make_figure(populations, depth_scores, all_scores,
                os_time, os_event, scorecard):
    log("")
    log("=" * 65)
    log("STEP 10: GENERATE FIGURE")
    log("=" * 65)

    fig = plt.figure(figsize=(22, 20))
    fig.patch.set_facecolor("#0d1117")
    gs  = gridspec.GridSpec(3, 4, figure=fig,
                            hspace=0.50, wspace=0.40)

    BLUE   = "#58a6ff"
    GREEN  = "#3fb950"
    RED    = "#f85149"
    ORANGE = "#d29922"
    PURPLE = "#bc8cff"
    GREY   = "#6e7681"

    def ax_style(ax):
        ax.set_facecolor("#161b22")
        for sp in ax.spines.values():
            sp.set_edgecolor("#30363d")
        ax.tick_params(colors="#8b949e", labelsize=7)
        ax.xaxis.label.set_color("#8b949e")
        ax.yaxis.label.set_color("#8b949e")

    lkw = dict(color="#8b949e", fontsize=7)
    tkw = dict(color="#e6edf3", fontsize=9,
               fontweight="bold", pad=6)

    cl_a = populations["claudin_low"]
    cl_b = populations["stratum_b"]

    # ── Panel 1: FOXP3/CD8A ratio vs depth scatter ──────────
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title("FOXP3/CD8A vs Depth\n(S3-P1)", **tkw)
    if all_scores.get("foxp3_cd8a") and depth_scores:
        dv = [depth_scores[s] for s in cl_a
              if s in depth_scores and s in all_scores["foxp3_cd8a"]
              and not np.isnan(all_scores["foxp3_cd8a"][s])]
        rv = [all_scores["foxp3_cd8a"][s] for s in cl_a
              if s in depth_scores and s in all_scores["foxp3_cd8a"]
              and not np.isnan(all_scores["foxp3_cd8a"][s])]
        if len(dv) > 3:
            ax1.scatter(dv, rv, c=BLUE, alpha=0.35, s=10,
                        edgecolors="none")
            z = np.polyfit(dv, rv, 1)
            xr = np.linspace(min(dv), max(dv), 50)
            ax1.plot(xr, np.polyval(z, xr),
                     color=ORANGE, lw=1.5, linestyle="--")
            r, p = stats.pearsonr(dv, rv)
            ax1.text(0.05, 0.92, f"r={r:.3f}  {fmt_p(p)}",
                     transform=ax1.transAxes,
                     color="#e6edf3", fontsize=7)
    ax1.set_xlabel("Depth Score", **lkw)
    ax1.set_ylabel("FOXP3/CD8A ratio", **lkw)
    ax_style(ax1)

    # ── Panel 2: FOXP3/CD8A tertile KM (Stratum B) ──────────
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_title("FOXP3/CD8A vs OS\nStratum B (S3-P2)", **tkw)
    res2 = tertile_split(cl_b, all_scores.get("foxp3_cd8a", {}),
                         os_time, os_event)
    T2_h, E2_h, T2_l, E2_l, T2_m, E2_m = res2[:6]
    for T, E, lbl, c in [
        (T2_h, E2_h, f"High ratio (n={len(T2_h)})", RED),
        (T2_m, E2_m, f"Mid  (n={len(T2_m)})",       GREY),
        (T2_l, E2_l, f"Low ratio  (n={len(T2_l)})", GREEN),
    ]:
        if len(T) >= MIN_SURV_N:
            kmf = KaplanMeierFitter()
            kmf.fit(T, E, label=lbl)
            kmf.plot_survival_function(ax=ax2, ci_show=False,
                                       color=c, linewidth=1.8)
    ax2.set_xlabel("Time (days)", **lkw)
    ax2.set_ylabel("Survival probability", **lkw)
    ax2.legend(fontsize=5.5, labelcolor="#8b949e",
               facecolor="#161b22", edgecolor="#30363d")
    ax_style(ax2)

    # ── Panel 3: CT composite vs depth scatter ───────────────
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.set_title("CT Antigen vs Depth\n(S3-P4)", **tkw)
    if all_scores.get("ct_composite"):
        dv3 = [depth_scores[s] for s in cl_a
               if s in depth_scores and s in all_scores["ct_composite"]
               and not np.isnan(all_scores["ct_composite"][s])]
        cv3 = [all_scores["ct_composite"][s] for s in cl_a
               if s in depth_scores and s in all_scores["ct_composite"]
               and not np.isnan(all_scores["ct_composite"][s])]
        if len(dv3) > 3:
            ax3.scatter(dv3, cv3, c=PURPLE, alpha=0.35, s=10,
                        edgecolors="none")
            z3 = np.polyfit(dv3, cv3, 1)
            xr3 = np.linspace(min(dv3), max(dv3), 50)
            ax3.plot(xr3, np.polyval(z3, xr3),
                     color=ORANGE, lw=1.5, linestyle="--")
            r3, p3 = stats.pearsonr(dv3, cv3)
            ax3.text(0.05, 0.92, f"r={r3:.3f}  {fmt_p(p3)}",
                     transform=ax3.transAxes,
                     color="#e6edf3", fontsize=7)
    else:
        ax3.text(0.5, 0.5, "CT antigen data\nnot available",
                 transform=ax3.transAxes, ha="center",
                 va="center", color="#8b949e")
    ax3.set_xlabel("Depth Score", **lkw)
    ax3.set_ylabel("CT Antigen Composite", **lkw)
    ax_style(ax3)

    # ── Panel 4: CT composite KM (Stratum B) ────────────────
    ax4 = fig.add_subplot(gs[0, 3])
    ax4.set_title("CT Antigen vs OS\nStratum B (S3-P5)", **tkw)
    if all_scores.get("ct_composite"):
        res4 = tertile_split(cl_b, all_scores["ct_composite"],
                             os_time, os_event)
        T4_h, E4_h, T4_l, E4_l, T4_m, E4_m = res4[:6]
        for T, E, lbl, c in [
            (T4_h, E4_h, f"CT-high (n={len(T4_h)})", RED),
            (T4_m, E4_m, f"CT-mid  (n={len(T4_m)})", GREY),
            (T4_l, E4_l, f"CT-low  (n={len(T4_l)})", GREEN),
        ]:
            if len(T) >= MIN_SURV_N:
                kmf = KaplanMeierFitter()
                kmf.fit(T, E, label=lbl)
                kmf.plot_survival_function(ax=ax4, ci_show=False,
                                           color=c, linewidth=1.8)
        ax4.set_xlabel("Time (days)", **lkw)
        ax4.set_ylabel("Survival probability", **lkw)
        ax4.legend(fontsize=5.5, labelcolor="#8b949e",
                   facecolor="#161b22", edgecolor="#30363d")
    else:
        ax4.text(0.5, 0.5, "CT antigen data\nnot available",
                 transform=ax4.transAxes, ha="center",
                 va="center", color="#8b949e")
    ax_style(ax4)

    # ── Panel 5: Lineage memory KM (Stratum B) ──────────────
    ax5 = fig.add_subplot(gs[1, 0:2])
    ax5.set_title("Lineage Memory (FOXA1/SPDEF/GATA3) vs OS\n"
                  "Stratum B — Memory-low vs Memory-high (S3-P7)", **tkw)
    if all_scores.get("lineage_memory"):
        mem = all_scores["lineage_memory"]
        mem_b = [(s, mem[s]) for s in cl_b if s in mem
                 and not np.isnan(mem[s])]
        if len(mem_b) >= MIN_SURV_N * 2:
            med_b = np.median([v for _, v in mem_b])
            mh = [s for s, v in mem_b if v >= med_b]
            ml = [s for s, v in mem_b if v <  med_b]
            T_mh, E_mh, _ = get_surv_arrays(mh, os_time, os_event)
            T_ml, E_ml, _ = get_surv_arrays(ml, os_time, os_event)
            for T, E, lbl, c in [
                (T_mh, E_mh, f"Memory-high (n={len(T_mh)})", GREEN),
                (T_ml, E_ml, f"Memory-low  (n={len(T_ml)})", RED),
            ]:
                if len(T) >= MIN_SURV_N:
                    kmf = KaplanMeierFitter()
                    kmf.fit(T, E, label=lbl)
                    kmf.plot_survival_function(ax=ax5, ci_show=True,
                                               color=c, linewidth=2.0)
            ax5.set_xlabel("Time (days)", **lkw)
            ax5.set_ylabel("Survival probability", **lkw)
            ax5.legend(fontsize=7, labelcolor="#8b949e",
                       facecolor="#161b22", edgecolor="#30363d")
    else:
        ax5.text(0.5, 0.5, "Lineage memory data\nnot available",
                 transform=ax5.transAxes, ha="center",
                 va="center", color="#8b949e")
    ax_style(ax5)

    # ── Panel 6: Memory score vs depth scatter ──────────────
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.set_title("Memory Score vs Depth\n(S3-P9 check)", **tkw)
    if all_scores.get("lineage_memory"):
        mem = all_scores["lineage_memory"]
        dv6 = [depth_scores[s] for s in cl_a
               if s in depth_scores and s in mem
               and not np.isnan(mem[s])]
        mv6 = [mem[s] for s in cl_a
               if s in depth_scores and s in mem
               and not np.isnan(mem[s])]
        if len(dv6) > 3:
            ax6.scatter(dv6, mv6, c=GREEN, alpha=0.35, s=10,
                        edgecolors="none")
            z6 = np.polyfit(dv6, mv6, 1)
            xr6 = np.linspace(min(dv6), max(dv6), 50)
            ax6.plot(xr6, np.polyval(z6, xr6),
                     color=ORANGE, lw=1.5, linestyle="--")
            r6, p6 = stats.pearsonr(dv6, mv6)
            ax6.text(0.05, 0.92, f"r={r6:.3f}  {fmt_p(p6)}",
                     transform=ax6.transAxes,
                     color="#e6edf3", fontsize=7)
    ax6.set_xlabel("Depth Score", **lkw)
    ax6.set_ylabel("Memory Score (FOXA1/SPDEF/GATA3)", **lkw)
    ax_style(ax6)

    # ── Panel 7: Memory group boxplots ───────────────────────
    ax7 = fig.add_subplot(gs[1, 3])
    ax7.set_title("Memory Groups: Depth,\nCT Antigen, Treg Ratio", **tkw)
    if all_scores.get("lineage_memory"):
        mem = all_scores["lineage_memory"]
        mem_b = [(s, mem[s]) for s in cl_b if s in mem
                 and not np.isnan(mem[s])]
        if len(mem_b) >= MIN_SURV_N * 2:
            med_b = np.median([v for _, v in mem_b])
            mh_s = [s for s, v in mem_b if v >= med_b]
            ml_s = [s for s, v in mem_b if v <  med_b]
            metrics = []
            for label, slist, c in [("Mem-high", mh_s, GREEN),
                                     ("Mem-low",  ml_s, RED)]:
                d_v = [depth_scores[s] for s in slist
                       if s in depth_scores]
                metrics.append((label, d_v, c))
            for i, (lbl, vals, c) in enumerate(metrics):
                bp = ax7.boxplot(vals, positions=[i], widths=0.6,
                                 patch_artist=True,
                                 boxprops=dict(facecolor=c, alpha=0.6),
                                 medianprops=dict(color="#e6edf3", lw=2),
                                 whiskerprops=dict(color="#8b949e"),
                                 capprops=dict(color="#8b949e"),
                                 flierprops=dict(marker="o", alpha=0.3,
                                                 markerfacecolor=c,
                                                 markeredgewidth=0))
            ax7.set_xticks([0, 1])
            ax7.set_xticklabels(["Mem-high", "Mem-low"],
                                 color="#8b949e", fontsize=7)
            ax7.set_ylabel("Depth Score", **lkw)
    ax_style(ax7)

    # ── Panels 8–10: Per-CT gene depth correlations ─────────
    ct_avail = [g for g in CT_ANTIGEN_GENES if g in all_scores]
    ct_rs = {}
    for g in ct_avail:
        dv = [depth_scores[s] for s in cl_a
              if s in depth_scores and s in all_scores[g]
              and not np.isnan(all_scores[g][s])]
        gv = [all_scores[g][s] for s in cl_a
              if s in depth_scores and s in all_scores[g]
              and not np.isnan(all_scores[g][s])]
        if len(dv) > 3:
            r, _ = stats.pearsonr(dv, gv)
            ct_rs[g] = r

    ax8 = fig.add_subplot(gs[2, 0:2])
    ax8.set_title("CT Antigen Genes: Pearson r with Depth Score\n"
                  "(Axioms Prediction G test)", **tkw)
    if ct_rs:
        genes = list(ct_rs.keys())
        rs    = [ct_rs[g] for g in genes]
        colors = [GREEN if r > 0.20 else
                  (ORANGE if r > 0 else RED) for r in rs]
        bars = ax8.barh(genes, rs, color=colors, alpha=0.7,
                        edgecolor="#30363d")
        ax8.axvline(0, color="#8b949e", lw=0.8, linestyle="--")
        ax8.axvline(0.20, color=GREEN, lw=0.8, linestyle=":")
        ax8.set_xlabel("Pearson r (with depth score)", **lkw)
        ax8.text(0.21, 0.95, "r=0.20\nthreshold",
                 transform=ax8.get_yaxis_transform(),
                 color=GREEN, fontsize=6, va="top")
    else:
        ax8.text(0.5, 0.5, "CT antigen genes\nnot available",
                 transform=ax8.transAxes, ha="center",
                 va="center", color="#8b949e")
    ax_style(ax8)

    # ── Panel 9: Composite Treg:Eff ratio KM (Stratum B) ────
    ax9 = fig.add_subplot(gs[2, 2])
    ax9.set_title("Composite Treg:Effector vs OS\n"
                  "Stratum B (S3-P3)", **tkw)
    if all_scores.get("composite_treg"):
        res9 = tertile_split(cl_b, all_scores["composite_treg"],
                             os_time, os_event)
        T9_h, E9_h, T9_l, E9_l, T9_m, E9_m = res9[:6]
        for T, E, lbl, c in [
            (T9_h, E9_h, f"High (n={len(T9_h)})", RED),
            (T9_m, E9_m, f"Mid  (n={len(T9_m)})", GREY),
            (T9_l, E9_l, f"Low  (n={len(T9_l)})", GREEN),
        ]:
            if len(T) >= MIN_SURV_N:
                kmf = KaplanMeierFitter()
                kmf.fit(T, E, label=lbl)
                kmf.plot_survival_function(ax=ax9, ci_show=False,
                                           color=c, linewidth=1.8)
        ax9.set_xlabel("Time (days)", **lkw)
        ax9.set_ylabel("Survival probability", **lkw)
        ax9.legend(fontsize=5.5, labelcolor="#8b949e",
                   facecolor="#161b22", edgecolor="#30363d")
    ax_style(ax9)

    # ── Panel 10: CT vs FOXP3/CD8A scatter ──────────────────
    ax10 = fig.add_subplot(gs[2, 3])
    ax10.set_title("CT Antigen vs\nFOXP3/CD8A Ratio (S3-P6)", **tkw)
    if (all_scores.get("ct_composite") and
            all_scores.get("foxp3_cd8a")):
        ct_v = [all_scores["ct_composite"][s] for s in cl_a
                if s in all_scores["ct_composite"]
                and s in all_scores["foxp3_cd8a"]
                and not np.isnan(all_scores["ct_composite"][s])
                and not np.isnan(all_scores["foxp3_cd8a"][s])]
        rat_v = [all_scores["foxp3_cd8a"][s] for s in cl_a
                 if s in all_scores["ct_composite"]
                 and s in all_scores["foxp3_cd8a"]
                 and not np.isnan(all_scores["ct_composite"][s])
                 and not np.isnan(all_scores["foxp3_cd8a"][s])]
        if len(ct_v) > 3:
            ax10.scatter(ct_v, rat_v, c=PURPLE, alpha=0.35, s=10,
                         edgecolors="none")
            z10 = np.polyfit(ct_v, rat_v, 1)
            xr10 = np.linspace(min(ct_v), max(ct_v), 50)
            ax10.plot(xr10, np.polyval(z10, xr10),
                      color=ORANGE, lw=1.5, linestyle="--")
            r10, p10 = stats.pearsonr(ct_v, rat_v)
            ax10.text(0.05, 0.92, f"r={r10:.3f}  {fmt_p(p10)}",
                      transform=ax10.transAxes,
                      color="#e6edf3", fontsize=7)
    ax10.set_xlabel("CT Antigen Composite", **lkw)
    ax10.set_ylabel("FOXP3/CD8A ratio", **lkw)
    ax_style(ax10)

    fig.suptitle(
        "CLAUDIN-LOW — SCRIPT 3\n"
        "Treg/Effector Ratio | CT Antigen Depth | Lineage Memory Subgroups\n"
        "OrganismCore | BRCA-S7g | 2026-03-05 | TYPE 4 ROOT LOCK",
        color="#e6edf3", fontsize=12, fontweight="bold", y=0.99)

    plt.savefig(FIG_FILE, dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close()
    log(f"  Figure saved: {FIG_FILE}")

# ============================================================
# STEP 11 — SCORECARD SUMMARY
# ============================================================

def print_scorecard(scorecard):
    log("")
    log("=" * 65)
    log("STEP 11: PREDICTION SCORECARD SUMMARY")
    log("=" * 65)
    log(f"\n  {'ID':<8} {'Analysis':<10} {'Gene':<22} "
        f"{'Str':<5} {'p':>8} {'metric':>8}  Status")
    log("  " + "─" * 75)
    for row in scorecard:
        p_str = (f"{row['p_value']:.3f}"
                 if not np.isnan(row['p_value']) else "nan")
        hr_str = (f"{row['HR_approx']:.3f}"
                  if not np.isnan(row['HR_approx']) else "nan")
        log(f"  {row['prediction']:<8} {row['analysis']:<10} "
            f"{str(row['gene']):<22} {str(row['stratum']):<5} "
            f"{p_str:>8} {hr_str:>8}  {row['confirmed']}")

    n_conf  = sum(1 for r in scorecard
                  if "CONFIRMED" in str(r["confirmed"]) and
                  "NOT" not in str(r["confirmed"]))
    log(f"\n  Confirmed: {n_conf} / {len(scorecard)}")
    pd.DataFrame(scorecard).to_csv(SCORECARD_FILE, index=False)
    log(f"  Scorecard saved: {SCORECARD_FILE}")

# ============================================================
# MAIN
# ============================================================

def main():
    log("=" * 65)
    log("BRCA CLAUDIN-LOW — SCRIPT 3")
    log("OrganismCore — Document BRCA-S7g | 2026-03-05")
    log("TYPE 4 — ROOT LOCK")
    log("")
    log("ANALYSES:")
    log("  A: Treg/effector ratio (Morel JCI 2017 test)")
    log("  B: CT antigen depth correlation (Axioms Pred G test)")
    log("  C: Lineage memory subgroups (Pommier 2020 test)")
    log(f"Output: {BASE_DIR}")
    log("=" * 65)

    scorecard = []

    # Steps 1–4: load and prepare
    expr, clin, surv    = load_data()
    populations          = build_populations(expr, clin)
    depth_scores         = compute_depth_scores(expr, populations)
    os_time, os_event    = build_survival_vectors(expr, populations, surv)

    # Step 5: compute all score dictionaries
    all_scores = compute_all_scores(expr, populations, depth_scores)

    # Steps 6–8: three analyses
    scorecard = analysis_a(populations, depth_scores, all_scores,
                           os_time, os_event, scorecard)
    scorecard = analysis_b(populations, depth_scores, all_scores,
                           os_time, os_event, scorecard)
    scorecard = analysis_c(populations, depth_scores, all_scores,
                           os_time, os_event, scorecard)

    # Step 9: save intermediate CSVs
    save_intermediate_files(populations, depth_scores,
                            all_scores, os_time, os_event)

    # Step 10: figure
    try:
        make_figure(populations, depth_scores, all_scores,
                    os_time, os_event, scorecard)
    except Exception as e:
        log(f"  Figure failed: {e}")

    # Step 11: scorecard
    print_scorecard(scorecard)

    log("")
    log("=" * 65)
    log("SCRIPT 3 COMPLETE")
    log("=" * 65)
    log(f"  Output: {BASE_DIR}")
    log(f"  Log:        {LOG_FILE}")
    log(f"  Figure:     {FIG_FILE}")
    log(f"  Scorecard:  {SCORECARD_FILE}")
    log(f"  Immune CSV: {IMMUNE_FILE}")
    log(f"  CT CSV:     {CT_FILE}")
    log(f"  Memory CSV: {MEMORY_FILE}")
    log("")
    log("  NEXT STEP:")
    log("    Write script3_results_and_reasoning.md (BRCA-S7h)")
    log("    Geometry read first — no prediction framing.")
    log("    Key questions:")
    log("      Does FOXP3/CD8A ratio predict OS better than")
    log("        total immune score (S2-P3 was null)?")
    log("      Does CT antigen track depth (first Pred G test)?")
    log("      Does lineage memory split separate OS curves?")
    log("      Do memory-low and memory-high differ in CT antigen")
    log("        AND Treg ratio (linking all three analyses)?")

    flush_log()


if __name__ == "__main__":
    main()
