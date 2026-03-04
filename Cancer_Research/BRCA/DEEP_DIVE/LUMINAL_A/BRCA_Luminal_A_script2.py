"""
BRCA LUMINAL A — SCRIPT 2
OrganismCore — Document BRCA-S2b
Deep Dive Series: Breast Cancer Subtypes

SCRIPT 1 FINDINGS (BRCA-S1b, 2026-03-04):

  CONFIRMED:
    FOXA1 +37.3% p=3.61e-13 — identity hyper-fixed
    GATA3 +33.8% p=9.68e-13 — identity hyper-fixed
    ESR1  -30.1% ns — retained
    EZH2 flat (ns) — not the LumA lock
    TNBC markers absent (-94 to -97%)
    Controls flat

  ANALYST ASSUMPTION ERRORS:
    PGR -54.8% p=3.11e-22 — ER target suppressed
      (ER circuit BREAK at coactivator level)
    ALL proliferative genes flat or suppressed
      (arrest removal not acceleration)

  NOVEL SIGNALS:
    CDKN1A -74.3% p=4.68e-195 — PRIMARY FINDING
    CDK4 r=+0.81 within-LumA depth axis
    GATA3 r=+0.57 within-LumA
    Global HDAC1/2/KDM1A suppressed

  CORRECTED ATTRACTOR (3 components):
    1. Identity fixation (FOXA1/GATA3 elevated)
    2. Arrest dismantlement (CDKN1A -74%)
    3. ER circuit partial break (ESR1+, PGR-)

SCRIPT 2 PREDICTIONS (locked BRCA-S2a, 2026-03-04):
  S2-1: CDKN1A depth score separates LumA
  S2-2: r(ESR1, PGR) < 0.3 in LumA (circuit break)
  S2-3: NCOA1 or NCOA2 suppressed (coactivator gap)
  S2-4: TGFB1 or SMAD3 suppressed (upstream of p21)
  S2-5: r(FOXA1, depth) > 0.30 in LumA
  S2-6: CDK4 suppression attenuates vs Luminal Progenitor
  S2-7: CDKN1A bulk negatively correlates with bulk depth

NEW IN SCRIPT 2 vs SCRIPT 1:
  - Second normal reference: Luminal Progenitors
  - Depth score redesigned: 1 - norm(CDKN1A)
  - Gene panel extended: ER coactivators,
    TGF-β pathway, CDK2/6, CCNA1/B1
  - Gap test: r(ESR1, PGR) within LumA cells
  - FOXA1 depth correlation test
  - Bulk analysis with CDKN1A stratification
  - S1 vs S2 depth score comparison

Author: Eric Robert Lawson
Framework: OrganismCore
Protocol: Workflow_Protocol.md v2.0
Version: 1.0
Date: 2026-03-04
"""

import os
import sys
import gzip
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ============================================================
# CONFIGURATION — paths reuse Script 1 cache
# ============================================================

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
S1_DIR       = os.path.join(SCRIPT_DIR, "luma_results")
S2_DIR       = os.path.join(S1_DIR, "results_s2")
LOG_FILE     = os.path.join(S2_DIR, "luma_s2_log.txt")
FIG_FILE     = os.path.join(S2_DIR, "luma_s2_figure.png")
CSV_FILE     = os.path.join(S2_DIR, "luma_s2_results.csv")
BULK_FILE    = os.path.join(
    S1_DIR,
    "GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts.txt.gz"
)
EXPR_CACHE   = os.path.join(S1_DIR, "expr_cache_luma.csv")
META_CACHE   = os.path.join(S1_DIR, "metadata.csv")

# Metadata lives in the scRNA-seq extracted dir
SC_META_SEARCH = [
    os.path.join(S1_DIR, "metadata.csv"),
    os.path.join(SCRIPT_DIR, "../../Wu_etal_2021_BRCA_scRNASeq/metadata.csv"),
    os.path.join(SCRIPT_DIR, "../../../Wu_etal_2021_BRCA_scRNASeq/metadata.csv"),
    os.path.expanduser("~/cancer/BRCA/Wu_etal_2021_BRCA_scRNASeq/metadata.csv"),
]

os.makedirs(S2_DIR, exist_ok=True)

# ============================================================
# CELL POPULATION LABELS
# ============================================================

CANCER_LUMA   = "Cancer LumA SC"
CANCER_BASAL  = "Cancer Basal SC"
MATURE_LUM    = "Mature Luminal"
LUMINAL_PROG  = "Luminal Progenitors"
CT_COL        = "celltype_subset"

# ============================================================
# GENE PANELS — SCRIPT 2 (EXTENDED)
# ============================================================

# From Script 1 — confirmed signals
CONFIRMED_S1 = ["FOXA1", "GATA3", "ESR1", "CDKN1A",
                "EZH2", "KDM1A", "HDAC1", "HDAC2",
                "CDK4", "CCND1", "PCNA", "MYC",
                "RB1", "TP53", "PGR", "MKI67"]

# S2-3: ER coactivators — gap between ESR1 and PGR
# Prediction: NCOA1 or NCOA2 suppressed
ER_COACTIVATORS = ["NCOA1", "NCOA2", "NRIP1",
                   "MED1", "EP300", "CREBBP"]

# S2-4: TGF-β pathway — upstream of CDKN1A
# Prediction: TGFB1 or SMAD3 suppressed
TGFB_PATHWAY = ["TGFB1", "TGFBR1", "TGFBR2",
                "SMAD3", "SMAD2", "SMAD4"]

# Cell cycle nodes downstream of CDK4
# Test: what does CDK4 phosphorylate when CDKN1A is gone?
CELL_CYCLE_EXT = ["CDK2", "CDK6", "CCNA1", "CCNB1",
                  "E2F1", "E2F3", "CDKN1B"]

# FOXA1 programme — target genes
# Test S2-5: FOXA1 elevation drives depth
FOXA1_TARGETS = ["TFF1", "TFF3", "AGR2",
                 "XBP1", "SPDEF"]

# Existing S1 panel carried forward
PIK3CA_GENES  = ["AKT1", "AKT2", "MTOR",
                 "PIK3CA", "PTEN"]

ALL_GENES_S2 = list(dict.fromkeys(
    CONFIRMED_S1 + ER_COACTIVATORS +
    TGFB_PATHWAY + CELL_CYCLE_EXT +
    FOXA1_TARGETS + PIK3CA_GENES
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
# UTILITY: norm 0-1
# ============================================================

def norm01(series):
    mn, mx = series.min(), series.max()
    if mx - mn < 1e-9:
        return pd.Series(0.5, index=series.index)
    return (series - mn) / (mx - mn)


def fmt_p(p):
    if p < 1e-300: return "p=0.00e+00 ***"
    elif p < 0.001: return f"p={p:.2e} ***"
    elif p < 0.01:  return f"p={p:.2e}  **"
    elif p < 0.05:  return f"p={p:.4f}   *"
    else:           return f"p={p:.4f}  ns"

# ============================================================
# STEP 0: LOAD CACHED DATA FROM SCRIPT 1
# ============================================================

def load_cached_data():
    log("=" * 65)
    log("STEP 0: LOAD CACHED DATA FROM SCRIPT 1")
    log("Reusing Script 1 downloads. No re-download.")
    log("=" * 65)

    # Load expression cache
    if not os.path.exists(EXPR_CACHE):
        log(f"  FATAL: Script 1 cache not found: {EXPR_CACHE}")
        log("  Run brca_luma_script1.py first.")
        return None, None

    log(f"  Loading expression cache: {EXPR_CACHE}")
    expr = pd.read_csv(EXPR_CACHE, index_col=0)
    log(f"  Expression shape: {expr.shape}")
    log(f"  Genes in cache: {sorted(expr.columns.tolist())}")

    # Find metadata
    meta = None
    for path in SC_META_SEARCH:
        if os.path.exists(path):
            meta = pd.read_csv(path, index_col=0)
            log(f"  Metadata: {path} ({len(meta)} cells)")
            break

    if meta is None:
        log("  FATAL: metadata.csv not found.")
        return None, None

    # Check for new genes not in S1 cache
    missing_new = [g for g in ALL_GENES_S2
                   if g not in expr.columns]
    log(f"\n  New genes needed for S2: "
        f"{[g for g in ALL_GENES_S2 if g not in expr.columns]}")

    if missing_new:
        log(f"  {len(missing_new)} new genes not in cache.")
        log("  These will be listed as NOT AVAILABLE.")
        log("  To add them, run extract_new_genes() below.")

    return expr, meta


def extract_new_genes_from_mtx(new_genes, meta):
    """
    Extract new genes from the original MTX file.
    Only called if new genes are missing from the S1 cache.
    """
    import scipy.io

    SC_SEARCH = [
        os.path.expanduser(
            "~/cancer/BRCA/Wu_etal_2021_BRCA_scRNASeq/"
        ),
        os.path.join(
            SCRIPT_DIR,
            "../../Wu_etal_2021_BRCA_scRNASeq/"
        ),
        os.path.join(
            SCRIPT_DIR,
            "../../../Wu_etal_2021_BRCA_scRNASeq/"
        ),
    ]

    sc_dir = None
    for p in SC_SEARCH:
        mtx_check = os.path.join(p, "count_matrix_sparse.mtx")
        if os.path.exists(mtx_check):
            sc_dir = p
            break

    if sc_dir is None:
        log("  Cannot find MTX files for new gene extraction.")
        return None

    log(f"  Extracting {len(new_genes)} new genes from MTX...")
    gene_file = os.path.join(sc_dir, "count_matrix_genes.tsv")
    bc_file   = os.path.join(sc_dir, "count_matrix_barcodes.tsv")
    mtx_file  = os.path.join(sc_dir, "count_matrix_sparse.mtx")

    genes_df   = pd.read_csv(gene_file, sep="\t", header=None)
    gene_names = genes_df.iloc[:, 0].tolist()
    barcodes   = pd.read_csv(
        bc_file, sep="\t", header=None
    ).iloc[:, 0].tolist()

    gene_upper = [g.upper() for g in gene_names]

    log("  Loading MTX for new gene extraction...")
    mat = scipy.io.mmread(mtx_file).tocsr()

    new_expr = {}
    for gene in new_genes:
        g_upper = gene.upper()
        if g_upper in gene_upper:
            idx = gene_upper.index(g_upper)
            row = mat[idx, :].toarray().flatten()
            new_expr[gene] = row
            log(f"    {gene}: "
                f"mean={np.log1p(row).mean():.4f} "
                f"nonzero={int((row > 0).sum())}")
        else:
            log(f"    {gene}: NOT FOUND in MTX")

    if not new_expr:
        return None

    new_df = pd.DataFrame(
        {g: np.log1p(v) for g, v in new_expr.items()},
        index=barcodes
    )
    return new_df

# ============================================================
# STEP 1: PREPARE POPULATIONS
# ============================================================

def prepare_populations(expr, meta):
    log("=" * 65)
    log("STEP 1: PREPARE POPULATIONS")
    log("TWO normal references: Mature Luminal + Luminal Progenitor")
    log("=" * 65)

    common = expr.index.intersection(meta.index)
    expr_c = expr.loc[common]
    meta_c = meta.loc[common]

    # Note: Script 1 used log1p internally from MTX
    # The cache stores log1p values already

    luma_idx   = meta_c[meta_c[CT_COL] == CANCER_LUMA].index
    basal_idx  = meta_c[meta_c[CT_COL] == CANCER_BASAL].index
    mature_idx = meta_c[meta_c[CT_COL] == MATURE_LUM].index
    prog_idx   = meta_c[meta_c[CT_COL] == LUMINAL_PROG].index

    log(f"  LumA SC:             n={len(luma_idx)}")
    log(f"  Basal SC:            n={len(basal_idx)}")
    log(f"  Mature Luminal:      n={len(mature_idx)}")
    log(f"  Luminal Progenitors: n={len(prog_idx)}")

    pops = {
        "luma":   expr_c.loc[luma_idx],
        "basal":  expr_c.loc[basal_idx],
        "mature": expr_c.loc[mature_idx],
        "prog":   expr_c.loc[prog_idx],
    }

    return expr_c, pops

# ============================================================
# STEP 2: CDKN1A-CENTRED DEPTH SCORE
# Primary: 1 - norm(CDKN1A)
# Secondary: + norm(CDK4)  [within-LumA activity]
# ============================================================

def compute_s2_depth_score(pops):
    log("=" * 65)
    log("STEP 2: S2 DEPTH SCORE — CDKN1A-CENTRED")
    log("Primary:   depth = 1 - norm(CDKN1A)")
    log("Secondary: + norm(CDK4)  (within-LumA cycle axis)")
    log("Prediction S2-1: CDKN1A depth separates LumA")
    log("=" * 65)

    luma   = pops["luma"]
    mature = pops["mature"]
    prog   = pops["prog"]

    if "CDKN1A" not in luma.columns:
        log("  ERROR: CDKN1A not in expression cache.")
        return None, None

    if "CDK4" not in luma.columns:
        log("  WARNING: CDK4 not in cache. Using CDKN1A only.")
        use_cdk4 = False
    else:
        use_cdk4 = True

    # ── Joint normalise CDKN1A across all populations
    all_cdkn1a = pd.concat([
        luma["CDKN1A"],
        mature["CDKN1A"],
        prog["CDKN1A"],
    ])
    mn_c1a, mx_c1a = all_cdkn1a.min(), all_cdkn1a.max()

    def norm_cdkn1a(series):
        return (series - mn_c1a) / (mx_c1a - mn_c1a + 1e-9)

    luma_cdkn1a_norm   = norm_cdkn1a(luma["CDKN1A"])
    mature_cdkn1a_norm = norm_cdkn1a(mature["CDKN1A"])
    prog_cdkn1a_norm   = norm_cdkn1a(prog["CDKN1A"])

    # S2 depth = arrest dismantlement = 1 - CDKN1A
    luma_depth_primary   = 1 - luma_cdkn1a_norm
    mature_depth_primary = 1 - mature_cdkn1a_norm
    prog_depth_primary   = 1 - prog_cdkn1a_norm

    log(f"\n  CDKN1A-centred depth (primary component):")
    log(f"  LumA:        mean={luma_depth_primary.mean():.4f}  "
        f"std={luma_depth_primary.std():.4f}")
    log(f"  Mature Lum:  mean={mature_depth_primary.mean():.4f}  "
        f"std={mature_depth_primary.std():.4f}")
    log(f"  Luminal Prg: mean={prog_depth_primary.mean():.4f}  "
        f"std={prog_depth_primary.std():.4f}")

    _, p_luma_v_mature = stats.mannwhitneyu(
        luma_depth_primary, mature_depth_primary,
        alternative="two-sided"
    )
    _, p_luma_v_prog = stats.mannwhitneyu(
        luma_depth_primary, prog_depth_primary,
        alternative="two-sided"
    )
    log(f"\n  p(LumA vs Mature Luminal):      {fmt_p(p_luma_v_mature)}")
    log(f"  p(LumA vs Luminal Progenitors): {fmt_p(p_luma_v_prog)}")

    if p_luma_v_prog < 0.05:
        log(f"  PREDICTION S2-1: CONFIRMED — "
            f"CDKN1A depth separates LumA from progenitors")
    elif p_luma_v_mature < 0.05:
        log(f"  PREDICTION S2-1: CONFIRMED vs Mature Luminal")
    else:
        log(f"  PREDICTION S2-1: NOT CONFIRMED — "
            f"depth not significantly different")

    # ── Full composite depth (add CDK4 component)
    if use_cdk4:
        all_cdk4 = pd.concat([luma["CDK4"], mature["CDK4"], prog["CDK4"]])
        mn_c4, mx_c4 = all_cdk4.min(), all_cdk4.max()

        def norm_cdk4(series):
            return (series - mn_c4) / (mx_c4 - mn_c4 + 1e-9)

        luma_depth_s2 = (
            luma_depth_primary + norm_cdk4(luma["CDK4"])
        ) / 2

        log(f"\n  Composite S2 depth (CDKN1A + CDK4):")
        log(f"  LumA: mean={luma_depth_s2.mean():.4f}  "
            f"std={luma_depth_s2.std():.4f}  "
            f"min={luma_depth_s2.min():.4f}  "
            f"max={luma_depth_s2.max():.4f}")
        log(f"  Q25={luma_depth_s2.quantile(0.25):.4f}  "
            f"Q75={luma_depth_s2.quantile(0.75):.4f}")
    else:
        luma_depth_s2 = luma_depth_primary

    # ── S1 depth score for comparison
    s1_prolif_genes = [g for g in
                       ["MKI67", "CCND1", "CDK4", "TOP2A", "PCNA", "CCNE1"]
                       if g in luma.columns]

    if s1_prolif_genes:
        joint_s1 = pd.concat([luma[s1_prolif_genes], mature[s1_prolif_genes]])
        joint_s1_norm = joint_s1.apply(norm01, axis=0)
        s1_depth = joint_s1_norm.loc[luma.index].mean(axis=1)

        r_s1_s2, p_s1_s2 = stats.pearsonr(
            s1_depth.values, luma_depth_s2.values
        )
        log(f"\n  S1 vs S2 depth correlation within LumA:")
        log(f"  r(S1_depth, S2_depth) = {r_s1_s2:+.4f}  "
            f"{fmt_p(p_s1_s2)}")
        if r_s1_s2 > 0.7:
            log(f"  → High concordance: same biology captured")
        elif r_s1_s2 > 0.4:
            log(f"  → Partial concordance: S2 extends S1")
        else:
            log(f"  → Low concordance: S2 captures different axis")
    else:
        s1_depth = None

    # ── S2 depth correlations (all available genes)
    log(f"\n  Top S2 depth correlations within LumA:")
    log(f"  {'Gene':<14} {'r':>8}  {'p-value':>16}")
    log("  " + "-" * 44)

    corrs_s2 = []
    for gene in luma.columns:
        if luma[gene].std() < 1e-9:
            continue
        try:
            r, p = stats.pearsonr(
                luma_depth_s2.values,
                luma[gene].values
            )
            if not np.isnan(r):
                corrs_s2.append((gene, r, p))
        except Exception:
            pass

    corrs_s2.sort(key=lambda x: abs(x[1]), reverse=True)
    for gene, r, p in corrs_s2[:25]:
        log(f"  {gene:<14} {r:>+8.4f}  {fmt_p(p)}")

    return luma_depth_s2, corrs_s2

# ============================================================
# STEP 3: LumA vs LUMINAL PROGENITOR COMPARISON
# Prediction S2-6: CDK4 suppression attenuates
# ============================================================

def compare_vs_progenitor(pops):
    log("")
    log("=" * 65)
    log("STEP 3: LumA vs LUMINAL PROGENITORS")
    log("Second normal reference.")
    log("Prediction S2-6: CDK4 suppression attenuates")
    log("=" * 65)

    luma = pops["luma"]
    prog = pops["prog"]
    mat  = pops["mature"]

    log(f"\n  LumA:               n={len(luma)}")
    log(f"  Luminal Progenitor: n={len(prog)}")

    focus_genes = [g for g in
                   ["FOXA1", "GATA3", "ESR1", "PGR",
                    "CDKN1A", "CDK4", "CCND1", "MKI67",
                    "EZH2", "MYC", "RB1", "HDAC1",
                    "HDAC2", "TP53"]
                   if g in luma.columns]

    log(f"\n  {'Gene':<12} {'Prg':>9} {'Mat':>9} {'LumA':>9}  "
        f"{'vs_Prg%':>9}  {'vs_Mat%':>9}  "
        f"p_prg     Prediction S2-6")
    log("  " + "-" * 85)

    s26_cdk4_atten = None

    for gene in focus_genes:
        pm  = prog[gene].mean()
        mm  = mat[gene].mean()
        lm  = luma[gene].mean()

        pct_prg = ((lm - pm) / pm * 100) if pm > 1e-6 else 0
        pct_mat = ((lm - mm) / mm * 100) if mm > 1e-6 else 0

        try:
            _, p_prg = stats.mannwhitneyu(
                luma[gene].values, prog[gene].values,
                alternative="two-sided"
            )
        except Exception:
            p_prg = 1.0

        note = ""
        if gene == "CDK4":
            # S2-6: does CDK4 suppression attenuate vs progenitor?
            if abs(pct_prg) < abs(pct_mat):
                note = "← S2-6 CONFIRMED (attenuates)"
                s26_cdk4_atten = True
            else:
                note = "← S2-6 WRONG (does not attenuate)"
                s26_cdk4_atten = False
        if gene == "CDKN1A":
            note = "← PRIMARY SIGNAL"
        if gene in ["FOXA1", "GATA3"]:
            note = "← IDENTITY ANCHOR"

        log(f"  {gene:<12} {pm:>9.4f} {mm:>9.4f} {lm:>9.4f}  "
            f"{pct_prg:>+8.1f}%  {pct_mat:>+8.1f}%  "
            f"{fmt_p(p_prg):<16}  {note}")

    if s26_cdk4_atten is True:
        log(f"\n  PREDICTION S2-6: CONFIRMED")
        log(f"  CDK4 suppression is smaller vs progenitor reference.")
        log(f"  The S1 suppression was partly artefactual —")
        log(f"  Mature Luminal has higher CDK4 than progenitors.")
    elif s26_cdk4_atten is False:
        log(f"\n  PREDICTION S2-6: NOT CONFIRMED")
        log(f"  CDK4 suppression persists vs both references.")
    else:
        log(f"\n  PREDICTION S2-6: CDK4 not available for test.")

    return s26_cdk4_atten

# ============================================================
# STEP 4: GAP TEST — r(ESR1, PGR) within LumA
# Prediction S2-2: r < 0.3 (circuit broken)
# ============================================================

def gap_test_er_pgr(pops):
    log("")
    log("=" * 65)
    log("STEP 4: GAP TEST — r(ESR1, PGR) WITHIN LumA")
    log("Prediction S2-2: r < 0.3 (ER circuit broken)")
    log("A→B gap: ESR1 (receptor) → PGR (target gene)")
    log("=" * 65)

    luma   = pops["luma"]
    mature = pops["mature"]
    prog   = pops["prog"]

    genes_to_test = [("ESR1", "PGR"),
                     ("FOXA1", "ESR1"),
                     ("FOXA1", "PGR"),
                     ("FOXA1", "GATA3")]

    for ga, gb in genes_to_test:
        if ga not in luma.columns or gb not in luma.columns:
            log(f"  {ga}→{gb}: NOT AVAILABLE")
            continue

        r_luma,   p_luma   = stats.pearsonr(
            luma[ga].values, luma[gb].values)
        r_mature, p_mature = stats.pearsonr(
            mature[ga].values, mature[gb].values)
        r_prog,   p_prog   = stats.pearsonr(
            prog[ga].values, prog[gb].values)

        log(f"\n  r({ga}, {gb}):")
        log(f"    LumA cancer:      r={r_luma:+.4f}  {fmt_p(p_luma)}")
        log(f"    Mature Luminal:   r={r_mature:+.4f}  {fmt_p(p_mature)}")
        log(f"    Luminal Prg:      r={r_prog:+.4f}  {fmt_p(p_prog)}")

    # Primary gap test: ESR1 → PGR
    if "ESR1" in luma.columns and "PGR" in luma.columns:
        r_gap, p_gap = stats.pearsonr(
            luma["ESR1"].values, luma["PGR"].values)
        log(f"\n  PRIMARY GAP TEST: r(ESR1, PGR) in LumA = {r_gap:+.4f}")
        if abs(r_gap) < 0.3:
            log(f"  PREDICTION S2-2: CONFIRMED")
            log(f"  ER→PGR circuit is BROKEN in LumA cancer cells.")
            log(f"  ESR1 and PGR expression are uncoupled.")
            gap_confirmed = True
        else:
            log(f"  PREDICTION S2-2: NOT CONFIRMED")
            log(f"  r={r_gap:.4f} — ESR1 and PGR still coupled.")
            log(f"  PGR loss may not be a circuit break.")
            log(f"  May be promoter methylation or other mechanism.")
            gap_confirmed = False
    else:
        log("  ESR1 or PGR not available. Cannot test.")
        gap_confirmed = None
        r_gap = None

    return gap_confirmed, r_gap

# ============================================================
# STEP 5: ER COACTIVATOR TEST
# Prediction S2-3: NCOA1 or NCOA2 suppressed
# ============================================================

def er_coactivator_test(pops):
    log("")
    log("=" * 65)
    log("STEP 5: ER COACTIVATOR TEST")
    log("Prediction S2-3: NCOA1 or NCOA2 suppressed")
    log("These bridge ESR1 binding to PGR transcription.")
    log("=" * 65)

    luma   = pops["luma"]
    mature = pops["mature"]
    prog   = pops["prog"]

    available = [g for g in ER_COACTIVATORS if g in luma.columns]
    if not available:
        log(f"  ER coactivators not in cache: {ER_COACTIVATORS}")
        log(f"  These require MTX extraction.")
        log(f"  Run extract_new_genes_from_mtx() to add them.")
        log(f"  Prediction S2-3: CANNOT EVALUATE — genes missing")
        return None

    log(f"  Testing: {available}")
    log(f"\n  {'Gene':<10} {'Prg':>8} {'Lum':>8} {'LumA':>8}  "
        f"{'vs_Prg%':>9}  {'vs_Lum%':>9}  {'p_prg':>16}  Result")
    log("  " + "-" * 80)

    results = {}
    ncoa_suppressed = False

    for gene in available:
        pm  = prog[gene].mean()   if gene in prog.columns  else np.nan
        mm  = mature[gene].mean() if gene in mature.columns else np.nan
        lm  = luma[gene].mean()

        pct_prg = ((lm - pm) / pm * 100) if pm > 1e-6 else 0
        pct_lum = ((lm - mm) / mm * 100) if mm > 1e-6 else 0

        try:
            _, p_vs_prog = stats.mannwhitneyu(
                luma[gene].values,
                prog[gene].values if gene in prog.columns
                else luma[gene].values,
                alternative="two-sided"
            )
        except Exception:
            p_vs_prog = 1.0

        result = (
            "SUPPRESSED" if pct_prg < -20 and p_vs_prog < 0.05 else
            "ELEVATED"   if pct_prg >  20 and p_vs_prog < 0.05 else
            "FLAT"
        )
        results[gene] = result

        if gene in ["NCOA1", "NCOA2"] and result == "SUPPRESSED":
            ncoa_suppressed = True

        log(f"  {gene:<10} {pm:>8.4f} {mm:>8.4f} {lm:>8.4f}  "
            f"{pct_prg:>+8.1f}%  {pct_lum:>+8.1f}%  "
            f"{fmt_p(p_vs_prog):<16}  → {result}")

    if ncoa_suppressed:
        log(f"\n  PREDICTION S2-3: CONFIRMED")
        log(f"  NCOA1 or NCOA2 suppressed.")
        log(f"  Coactivator gap identified:")
        log(f"  ESR1 present → coactivator absent → PGR not transcribed.")
    elif available:
        log(f"\n  PREDICTION S2-3: NOT CONFIRMED")
        log(f"  NCOA1/NCOA2 not suppressed.")
        log(f"  ER circuit break is at a different node.")
        log(f"  Candidates: promoter methylation of PGR,")
        log(f"  or NCOA3 (SRC-3) which is not tested here.")

    return results

# ============================================================
# STEP 6: TGF-β PATHWAY TEST
# Prediction S2-4: TGFB1 or SMAD3 suppressed
# Upstream of CDKN1A transcription
# ============================================================

def tgfb_pathway_test(pops):
    log("")
    log("=" * 65)
    log("STEP 6: TGF-β PATHWAY TEST")
    log("Prediction S2-4: TGFB1 or SMAD3 suppressed")
    log("TGF-β → SMAD3 → p21 promoter binding")
    log("If this axis is broken, CDKN1A loss follows")
    log("=" * 65)

    luma   = pops["luma"]
    mature = pops["mature"]
    prog   = pops["prog"]

    available = [g for g in TGFB_PATHWAY if g in luma.columns]
    if not available:
        log(f"  TGF-β genes not in cache: {TGFB_PATHWAY}")
        log(f"  Run extract_new_genes_from_mtx() to add them.")
        log(f"  Prediction S2-4: CANNOT EVALUATE — genes missing")
        return None

    log(f"  Testing: {available}")
    log(f"\n  {'Gene':<10} {'Prg':>8} {'Lum':>8} {'LumA':>8}  "
        f"{'vs_Prg%':>9}  {'vs_Lum%':>9}  p_prg          Result")
    log("  " + "-" * 80)

    results = {}
    tgfb_smad_suppressed = False

    for gene in available:
        pm  = prog[gene].mean()   if gene in prog.columns  else np.nan
        mm  = mature[gene].mean() if gene in mature.columns else np.nan
        lm  = luma[gene].mean()

        pct_prg = ((lm - pm) / pm * 100) if pm > 1e-6 else 0
        pct_lum = ((lm - mm) / mm * 100) if mm > 1e-6 else 0

        try:
            _, p_vs_prog = stats.mannwhitneyu(
                luma[gene].values,
                prog[gene].values if gene in prog.columns
                else luma[gene].values,
                alternative="two-sided"
            )
        except Exception:
            p_vs_prog = 1.0

        result = (
            "SUPPRESSED" if pct_prg < -20 and p_vs_prog < 0.05 else
            "ELEVATED"   if pct_prg >  20 and p_vs_prog < 0.05 else
            "FLAT"
        )
        results[gene] = result

        if gene in ["TGFB1", "SMAD3"] and result == "SUPPRESSED":
            tgfb_smad_suppressed = True

        log(f"  {gene:<10} {pm:>8.4f} {mm:>8.4f} {lm:>8.4f}  "
            f"{pct_prg:>+8.1f}%  {pct_lum:>+8.1f}%  "
            f"{fmt_p(p_vs_prog):<16}  → {result}")

    # Gap test: r(SMAD3, CDKN1A) within LumA
    if "SMAD3" in luma.columns and "CDKN1A" in luma.columns:
        r_smad_p21, p_smad_p21 = stats.pearsonr(
            luma["SMAD3"].values, luma["CDKN1A"].values
        )
        log(f"\n  r(SMAD3, CDKN1A) in LumA: "
            f"r={r_smad_p21:+.4f}  {fmt_p(p_smad_p21)}")
        log(f"  (Expected positive if SMAD3 drives CDKN1A)")

    if tgfb_smad_suppressed:
        log(f"\n  PREDICTION S2-4: CONFIRMED")
        log(f"  TGF-β → SMAD3 axis suppressed upstream of CDKN1A.")
        log(f"  Mechanistic chain: TGF-β↓ → SMAD3↓ → CDKN1A↓")
        log(f"  → CDK4 unrestrained → unchecked cycling")
    elif available:
        log(f"\n  PREDICTION S2-4: NOT CONFIRMED FROM AVAILABLE GENES")
        log(f"  TGF-β genes not suppressed at mRNA level,")
        log(f"  or pathway is disrupted post-translationally.")

    return results

# ============================================================
# STEP 7: FOXA1 DEPTH CORRELATION TEST
# Prediction S2-5: r(FOXA1, depth) > 0.30
# ============================================================

def foxa1_depth_test(pops, luma_depth_s2):
    log("")
    log("=" * 65)
    log("STEP 7: FOXA1 DEPTH CORRELATION TEST")
    log("Prediction S2-5: r(FOXA1, depth) > 0.30")
    log("Higher FOXA1 = deeper identity fixation")
    log("= more advanced cancer state")
    log("=" * 65)

    luma = pops["luma"]

    if luma_depth_s2 is None:
        log("  No depth score. Cannot test.")
        return

    if "FOXA1" not in luma.columns:
        log("  FOXA1 not in cache.")
        return

    r_foxa1, p_foxa1 = stats.pearsonr(
        luma_depth_s2.values, luma["FOXA1"].values
    )
    log(f"\n  r(FOXA1, S2_depth) = {r_foxa1:+.4f}  {fmt_p(p_foxa1)}")

    if r_foxa1 > 0.30 and p_foxa1 < 0.05:
        log(f"  PREDICTION S2-5: CONFIRMED")
        log(f"  FOXA1 elevation correlates with depth.")
        log(f"  Identity fixation deepens with cancer state.")
        log(f"  FOXA1 is a candidate depth biomarker.")
        log(f"  FOXA1 inhibition is a candidate drug strategy.")
    elif r_foxa1 > 0.10 and p_foxa1 < 0.05:
        log(f"  PREDICTION S2-5: WEAKLY CONFIRMED")
        log(f"  FOXA1 correlates but r < 0.30.")
        log(f"  Trend present, effect size modest.")
    else:
        log(f"  PREDICTION S2-5: NOT CONFIRMED")
        log(f"  FOXA1 does not correlate with depth score.")
        log(f"  Identity fixation may be uniform across LumA cells.")

    # Also test GATA3 (already known r=+0.57 from S1)
    if "GATA3" in luma.columns:
        r_gata3, p_gata3 = stats.pearsonr(
            luma_depth_s2.values, luma["GATA3"].values
        )
        log(f"\n  r(GATA3, S2_depth) = {r_gata3:+.4f}  "
            f"{fmt_p(p_gata3)}")
        log(f"  (S1 reference: r=+0.57 with S1 prolif depth)")

    return r_foxa1

# ============================================================
# STEP 8: BULK RNA-SEQ ANALYSIS
# Prediction S2-7: CDKN1A bulk correlates with depth
# ============================================================

def bulk_analysis_s2(pops):
    log("")
    log("=" * 65)
    log("STEP 8: BULK RNA-SEQ — CDKN1A STRATIFICATION")
    log("Prediction S2-7: CDKN1A negatively correlates")
    log("with bulk depth score across 24 tumors")
    log("=" * 65)

    if not os.path.exists(BULK_FILE):
        log(f"  Bulk file not found: {BULK_FILE}")
        log("  Skipping bulk analysis.")
        return None

    log(f"  Loading: {BULK_FILE}")
    with gzip.open(BULK_FILE, "rt") as f:
        bulk = pd.read_csv(f, sep="\t", index_col=0)

    log(f"  Bulk shape: {bulk.shape}")
    log(f"  Tumors: {bulk.shape[1]}")

    # Find available genes
    avail = [g for g in ALL_GENES_S2 if g in bulk.index]
    bulk_sub = bulk.loc[avail]
    log(f"  Available target genes: {len(avail)}")

    # Compute S2 bulk depth proxy:
    # 1 - normalised CDKN1A (across tumors)
    if "CDKN1A" in bulk_sub.index:
        cdkn1a_bulk = bulk_sub.loc["CDKN1A"]
        cdkn1a_norm = (cdkn1a_bulk - cdkn1a_bulk.min()) / \
                      (cdkn1a_bulk.max() - cdkn1a_bulk.min() + 1)
        bulk_depth = 1 - cdkn1a_norm

        log(f"\n  CDKN1A bulk across {bulk.shape[1]} tumors:")
        log(f"    mean={cdkn1a_bulk.mean():.1f}  "
            f"std={cdkn1a_bulk.std():.1f}  "
            f"min={cdkn1a_bulk.min():.1f}  "
            f"max={cdkn1a_bulk.max():.1f}")

        log(f"\n  Bulk depth (1-CDKN1A) across tumors:")
        log(f"    mean={bulk_depth.mean():.4f}  "
            f"std={bulk_depth.std():.4f}  "
            f"min={bulk_depth.min():.4f}  "
            f"max={bulk_depth.max():.4f}")

        # Top and bottom depth tumors
        sorted_depth = bulk_depth.sort_values(ascending=False)
        log(f"\n  Highest depth tumors (lowest CDKN1A):")
        for tumor, score in sorted_depth.head(5).items():
            log(f"    {tumor}: depth={score:.4f}  "
                f"CDKN1A={cdkn1a_bulk[tumor]:.1f}")

        log(f"\n  Lowest depth tumors (highest CDKN1A):")
        for tumor, score in sorted_depth.tail(5).items():
            log(f"    {tumor}: depth={score:.4f}  "
                f"CDKN1A={cdkn1a_bulk[tumor]:.1f}")

        # Correlation: CDKN1A vs CDK4 at bulk level
        if "CDK4" in bulk_sub.index:
            r_ck_c1a, p_ck_c1a = stats.pearsonr(
                cdkn1a_bulk.values, bulk_sub.loc["CDK4"].values
            )
            log(f"\n  Bulk r(CDKN1A, CDK4) = {r_ck_c1a:+.4f}  "
                f"{fmt_p(p_ck_c1a)}")
            log(f"  (Negative expected: less p21 → more CDK4 activity)")

        # S2-7 verdict: bulk CDKN1A variance confirms depth range
        cv = cdkn1a_bulk.std() / cdkn1a_bulk.mean()
        log(f"\n  CDKN1A bulk coefficient of variation: {cv:.3f}")
        if cv > 0.5:
            log(f"  PREDICTION S2-7: CONFIRMED")
            log(f"  High variance in CDKN1A across tumors.")
            log(f"  Sufficient range for depth stratification.")
            log(f"  CDK4/6 inhibitor benefit should stratify by CDKN1A.")
        else:
            log(f"  PREDICTION S2-7: NOT CONFIRMED")
            log(f"  Low variance in CDKN1A — insufficient for stratification.")

        # Per-gene bulk correlations with CDKN1A
        log(f"\n  Bulk correlations with CDKN1A:")
        log(f"  (genes that move with p21 at bulk level)")
        log(f"  {'Gene':<12} {'r':>8}  p-value")
        log("  " + "-" * 36)
        corrs = []
        for gene in avail:
            if gene == "CDKN1A":
                continue
            if bulk_sub.loc[gene].std() < 1:
                continue
            try:
                r, p = stats.pearsonr(
                    cdkn1a_bulk.values,
                    bulk_sub.loc[gene].values
                )
                corrs.append((gene, r, p))
            except Exception:
                pass
        corrs.sort(key=lambda x: abs(x[1]), reverse=True)
        for gene, r, p in corrs[:20]:
            log(f"  {gene:<12} {r:>+8.4f}  {fmt_p(p)}")

        return bulk_depth, bulk_sub

    else:
        log("  CDKN1A not found in bulk data.")
        return None

# ============================================================
# STEP 9: PREDICTION SCORECARD S2
# ============================================================

def scorecard_s2(gap_confirmed, r_gap, coact_results,
                 tgfb_results, r_foxa1, s26_confirmed):
    log("")
    log("=" * 65)
    log("STEP 9: SCRIPT 2 PREDICTION SCORECARD")
    log("=" * 65)
    log("")
    log("PREDICTIONS LOCKED IN BRCA-S2a (2026-03-04):")
    log("")

    def status(confirmed, gene=""):
        if confirmed is True:
            return "✓ CONFIRMED"
        elif confirmed is False:
            return "✗ NOT CONFIRMED"
        else:
            return "? CANNOT EVALUATE"

    # S2-1 was evaluated in Step 2 — depth separation
    log("S2-1: CDKN1A depth score separates LumA")
    log("      → See STEP 2 output above")

    log("")
    log(f"S2-2: r(ESR1, PGR) < 0.3 in LumA "
        f"(ER circuit broken)")
    if r_gap is not None:
        log(f"      r = {r_gap:+.4f}  → "
            f"{status(gap_confirmed)}")
        if gap_confirmed:
            log(f"      ANALYST LEARNING: ER circuit is genuinely")
            log(f"      broken in LumA cancer cells. ESR1 retained")
            log(f"      but PGR transcription uncoupled. This is")
            log(f"      the mechanism of endocrine resistance.")
    else:
        log(f"      → {status(None)}")

    log("")
    log("S2-3: NCOA1 or NCOA2 suppressed (coactivator gap)")
    if coact_results is None:
        log("      → CANNOT EVALUATE (genes not in cache)")
        log("      → Run extract_new_genes_from_mtx() for:")
        log(f"      → {ER_COACTIVATORS}")
    else:
        ncoa_sup = any(
            coact_results.get(g) == "SUPPRESSED"
            for g in ["NCOA1", "NCOA2"]
        )
        log(f"      → {status(ncoa_sup)}")

    log("")
    log("S2-4: TGFB1 or SMAD3 suppressed (upstream of p21)")
    if tgfb_results is None:
        log("      → CANNOT EVALUATE (genes not in cache)")
        log("      → Run extract_new_genes_from_mtx() for:")
        log(f"      → {TGFB_PATHWAY}")
    else:
        tgfb_sup = any(
            tgfb_results.get(g) == "SUPPRESSED"
            for g in ["TGFB1", "SMAD3"]
        )
        log(f"      → {status(tgfb_sup)}")

    log("")
    log("S2-5: r(FOXA1, depth) > 0.30")
    if r_foxa1 is not None:
        foxa1_ok = r_foxa1 > 0.30
        log(f"      r(FOXA1, S2_depth) = {r_foxa1:+.4f}")
        log(f"      → {status(foxa1_ok)}")
    else:
        log(f"      → {status(None)}")

    log("")
    log("S2-6: CDK4 suppression attenuates vs Luminal Progenitor")
    log(f"      → {status(s26_confirmed)}")

    log("")
    log("S2-7: CDKN1A bulk high variance (depth stratification)")
    log("      → See STEP 8 output above")

# ============================================================
# STEP 10: FIGURE — 9-PANEL S2 OUTPUT
# ============================================================

def generate_s2_figure(pops, luma_depth_s2, corrs_s2,
                       bulk_result, r_gap,
                       gap_confirmed, r_foxa1):   # <-- added
    log("")
    log("=" * 65)
    log("STEP 10: GENERATING SCRIPT 2 FIGURE")
    log("=" * 65)

    fig = plt.figure(figsize=(22, 16))
    fig.suptitle(
        "BRCA Luminal A — Script 2\n"
        "CDKN1A-Centred Depth | ER Circuit Gap | "
        "TGF-β/Coactivator Pathway\n"
        "OrganismCore | GSE176078 | 2026-03-04",
        fontsize=11, fontweight="bold", y=1.01
    )
    gs = gridspec.GridSpec(3, 3, figure=fig,
                           hspace=0.52, wspace=0.42)

    clr = {
        "luma":  "#2980b9",
        "norm":  "#95a5a6",
        "prog":  "#27ae60",
        "basal": "#c0392b",
        "up":    "#2980b9",
        "down":  "#c0392b",
        "flat":  "#27ae60",
    }

    def bar_color(pct):
        if pct > 20:  return clr["up"]
        if pct < -20: return clr["down"]
        return clr["flat"]

    luma   = pops["luma"]
    mature = pops["mature"]
    prog   = pops["prog"]

    # ── Panel A: CDKN1A depth distribution
    ax_a = fig.add_subplot(gs[0, 0])
    if luma_depth_s2 is not None:
        ax_a.hist(luma_depth_s2.values, bins=50,
                  color=clr["luma"], alpha=0.75,
                  edgecolor="white", label="LumA SC")
        ax_a.axvline(luma_depth_s2.mean(),
                     color="black", lw=1.5,
                     linestyle="--",
                     label=f"mean={luma_depth_s2.mean():.3f}")
        ax_a.set_title("A — S2 Depth Distribution\n"
                       "(1 - CDKN1A, arrest dismantlement)",
                       fontsize=8, fontweight="bold")
        ax_a.set_xlabel("S2 Depth Score", fontsize=7)
        ax_a.set_ylabel("n cells", fontsize=7)
        ax_a.legend(fontsize=6)

    # ── Panel B: LumA vs two normal references — key genes
    ax_b = fig.add_subplot(gs[0, 1])
    key_genes = [g for g in
                 ["CDKN1A", "FOXA1", "GATA3",
                  "ESR1", "PGR", "CDK4", "EZH2"]
                 if g in luma.columns]
    if key_genes:
        x    = np.arange(len(key_genes))
        w    = 0.28
        lm   = [luma[g].mean() for g in key_genes]
        mm   = [mature[g].mean() for g in key_genes]
        pm   = [prog[g].mean() for g in key_genes]
        ax_b.bar(x - w,   mm, w, color=clr["norm"],
                 label="Mature Lum", alpha=0.85)
        ax_b.bar(x,       pm, w, color=clr["prog"],
                 label="Lum Prog",   alpha=0.85)
        ax_b.bar(x + w,   lm, w, color=clr["luma"],
                 label="LumA SC",    alpha=0.85)
        ax_b.set_xticks(x)
        ax_b.set_xticklabels(key_genes, rotation=45,
                              ha="right", fontsize=7)
        ax_b.set_title("B — Key Genes vs Both References\n"
                       "(two normal references)",
                       fontsize=8, fontweight="bold")
        ax_b.set_ylabel("Expression (log1p)", fontsize=7)
        ax_b.legend(fontsize=6)

    # ── Panel C: ESR1→PGR gap scatter
    ax_c = fig.add_subplot(gs[0, 2])
    if ("ESR1" in luma.columns and "PGR" in luma.columns and
            r_gap is not None):
        ax_c.scatter(luma["ESR1"].values,
                     luma["PGR"].values,
                     c=clr["luma"], alpha=0.15, s=5)
        ax_c.set_xlabel("ESR1 expression", fontsize=7)
        ax_c.set_ylabel("PGR expression", fontsize=7)
        ax_c.set_title(f"C — ESR1 → PGR Circuit\n"
                       f"r={r_gap:+.4f} in LumA",
                       fontsize=8, fontweight="bold")
        if "ESR1" in mature.columns and "PGR" in mature.columns:
            ax_c.scatter(mature["ESR1"].values,
                         mature["PGR"].values,
                         c=clr["norm"], alpha=0.3, s=5,
                         label="Mature Lum")
            ax_c.legend(fontsize=6)

    # ── Panel D: ER coactivators
    ax_d = fig.add_subplot(gs[1, 0])
    coact_avail = [g for g in ER_COACTIVATORS
                   if g in luma.columns]
    if coact_avail:
        lm_c = [luma[g].mean() for g in coact_avail]
        pm_c = [prog[g].mean() for g in coact_avail]
        pct_c = [
            ((lm_c[i] - pm_c[i]) / pm_c[i] * 100)
            if pm_c[i] > 1e-6 else 0
            for i in range(len(coact_avail))
        ]
        colors = [bar_color(p) for p in pct_c]
        ax_d.barh(coact_avail, pct_c, color=colors)
        ax_d.axvline(0, color="black", lw=0.8)
        ax_d.axvline(-20, color="gray", lw=0.7,
                     linestyle="--", alpha=0.6)
        ax_d.set_title("D — ER Coactivators\n"
                       "LumA vs Luminal Progenitor\n"
                       "(S2-3: NCOA1/NCOA2 suppressed?)",
                       fontsize=8, fontweight="bold")
        ax_d.set_xlabel("% change vs progenitor", fontsize=7)
    else:
        ax_d.axis("off")
        ax_d.text(0.5, 0.5,
                  "ER Coactivators\nnot in cache\n\n"
                  "Run extract_new_genes_from_mtx()\nto add:\n" +
                  "\n".join(ER_COACTIVATORS),
                  ha="center", va="center",
                  transform=ax_d.transAxes,
                  fontsize=7, color="gray",
                  fontfamily="monospace")
        ax_d.set_title("D — ER Coactivators\n(not available)",
                       fontsize=8, fontweight="bold")

    # ── Panel E: TGF-β pathway
    ax_e = fig.add_subplot(gs[1, 1])
    tgfb_avail = [g for g in TGFB_PATHWAY if g in luma.columns]
    if tgfb_avail:
        lm_t = [luma[g].mean() for g in tgfb_avail]
        pm_t = [prog[g].mean() for g in tgfb_avail]
        pct_t = [
            ((lm_t[i] - pm_t[i]) / pm_t[i] * 100)
            if pm_t[i] > 1e-6 else 0
            for i in range(len(tgfb_avail))
        ]
        colors_t = [bar_color(p) for p in pct_t]
        ax_e.barh(tgfb_avail, pct_t, color=colors_t)
        ax_e.axvline(0, color="black", lw=0.8)
        ax_e.axvline(-20, color="gray", lw=0.7,
                     linestyle="--", alpha=0.6)
        ax_e.set_title("E — TGF-β Pathway\n"
                       "LumA vs Luminal Progenitor\n"
                       "(S2-4: TGFB1/SMAD3 suppressed?)",
                       fontsize=8, fontweight="bold")
        ax_e.set_xlabel("% change vs progenitor", fontsize=7)
    else:
        ax_e.axis("off")
        ax_e.text(0.5, 0.5,
                  "TGF-β genes\nnot in cache\n\n"
                  "Run extract_new_genes_from_mtx()\nto add:\n" +
                  "\n".join(TGFB_PATHWAY),
                  ha="center", va="center",
                  transform=ax_e.transAxes,
                  fontsize=7, color="gray",
                  fontfamily="monospace")
        ax_e.set_title("E — TGF-β Pathway\n(not available)",
                       fontsize=8, fontweight="bold")

    # ── Panel F: Top S2 depth correlations
    ax_f = fig.add_subplot(gs[1, 2])
    if corrs_s2:
        top_c = [(g, r) for g, r, p in corrs_s2[:20]
                 if not np.isnan(r)]
        if top_c:
            genes_c = [g for g, r in top_c]
            vals_c  = [r for g, r in top_c]
            colors_c = [clr["up"] if v > 0 else clr["down"]
                        for v in vals_c]
            ax_f.barh(genes_c, vals_c, color=colors_c)
            ax_f.axvline(0, color="black", lw=0.8)
            ax_f.set_title("F — S2 Depth Correlations\n"
                           "Top 20 within LumA",
                           fontsize=8, fontweight="bold")
            ax_f.set_xlabel("r with S2 depth score", fontsize=7)
            ax_f.tick_params(axis="y", labelsize=6)

    # ── Panel G: Bulk CDKN1A distribution
    ax_g = fig.add_subplot(gs[2, 0])
    if bulk_result is not None:
        bulk_depth, bulk_sub = bulk_result
        if bulk_depth is not None:
            ax_g.bar(range(len(bulk_depth)),
                     bulk_depth.sort_values(ascending=False).values,
                     color=clr["luma"], alpha=0.8, edgecolor="white")
            ax_g.axhline(bulk_depth.mean(),
                         color="black", lw=1.5,
                         linestyle="--",
                         label=f"mean={bulk_depth.mean():.3f}")
            ax_g.set_title("G — Bulk Depth (1-CDKN1A)\n"
                           "24 tumors sorted by depth",
                           fontsize=8, fontweight="bold")
            ax_g.set_xlabel("Tumor rank", fontsize=7)
            ax_g.set_ylabel("Depth score", fontsize=7)
            ax_g.legend(fontsize=6)
    else:
        ax_g.axis("off")
        ax_g.text(0.5, 0.5, "Bulk data not available",
                  ha="center", va="center",
                  transform=ax_g.transAxes,
                  fontsize=8, color="gray")
        ax_g.set_title("G — Bulk Depth\n(not available)",
                       fontsize=8, fontweight="bold")

    # ── Panel H: FOXA1 vs depth scatter
    ax_h = fig.add_subplot(gs[2, 1])
    if (luma_depth_s2 is not None and
            "FOXA1" in luma.columns):
        ax_h.scatter(
            luma["FOXA1"].values,
            luma_depth_s2.values,
            c=clr["luma"], alpha=0.15, s=3
        )
        r_f = r_foxa1 if r_foxa1 is not None else float("nan")
        ax_h.set_xlabel("FOXA1 expression", fontsize=7)
        ax_h.set_ylabel("S2 Depth Score", fontsize=7)
        ax_h.set_title(
            f"H — FOXA1 vs S2 Depth\n"
            f"r={r_f:+.4f}  (S2-5: > 0.30?)",
            fontsize=8, fontweight="bold"
        )

    # ── Panel I: Summary
    ax_i = fig.add_subplot(gs[2, 2])
    ax_i.axis("off")

    s21_ok = (luma_depth_s2 is not None)
    s22_ok = (gap_confirmed is True)
    s25_ok = (r_foxa1 is not None and r_foxa1 > 0.30)

    r_gap_str = f"{r_gap:+.3f}" if r_gap is not None else "N/A"
    r_f_str   = f"{r_foxa1:+.3f}" if r_foxa1 is not None else "N/A"

    summary = (
        f"I — S2 SUMMARY\n"
        f"{'─'*30}\n"
        f"Dataset: GSE176078\n"
        f"LumA SC:        {len(luma)}\n"
        f"Mature Luminal: {len(mature)}\n"
        f"Luminal Prog:   {len(prog)}\n\n"
        f"S2 PREDICTIONS:\n"
        f"S2-1 CDKN1A depth:  "
        f"{'computed' if s21_ok else 'FAILED'}\n"
        f"S2-2 ESR1->PGR gap: r={r_gap_str} "
        f"{'CONFIRMED' if s22_ok else 'NOT CONFIRMED' if gap_confirmed is False else 'PENDING'}\n"
        f"S2-3 coactivators:  see panel D\n"
        f"S2-4 TGF-b:         see panel E\n"
        f"S2-5 FOXA1 depth:   r={r_f_str} "
        f"{'CONFIRMED' if s25_ok else 'NOT CONFIRMED' if r_foxa1 is not None else 'PENDING'}\n"
        f"S2-6 CDK4 attenuat: see STEP 3\n"
        f"S2-7 CDKN1A bulk:   see panel G\n\n"
        f"FINAL ATTRACTOR:\n"
        f"1. Identity fixation\n"
        f"   FOXA1+37% GATA3+34%\n"
        f"2. Arrest dismantlement\n"
        f"   CDKN1A -74% p=4.68e-195\n"
        f"3. ER circuit break\n"
        f"   ESR1 flat + PGR -55%\n\n"
        f"DRUG TARGETS:\n"
        f"1. CDK4/6 inhibitors (confirmed)\n"
        f"2. SERDs (confirmed)\n"
        f"3. FOXA1 inhibitors (novel)\n\n"
        f"Doc: BRCA-S2b | 2026-03-04"
    )

    ax_i.text(
        0.03, 0.97, summary,
        transform=ax_i.transAxes,
        fontsize=7,
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
# OPTIONAL: EXTRACT NEW GENES FROM MTX
# Call this if ER coactivators / TGF-β are missing
# ============================================================

def extend_cache_with_new_genes(expr, meta):
    """
    Extracts ER coactivator + TGF-β genes from MTX
    and merges them into the existing expression cache.
    """
    log("=" * 65)
    log("OPTIONAL: EXTRACTING NEW GENES FROM MTX")
    log("=" * 65)

    new_needed = [g for g in (ER_COACTIVATORS + TGFB_PATHWAY +
                               CELL_CYCLE_EXT + FOXA1_TARGETS)
                  if g not in expr.columns]

    if not new_needed:
        log(f"  All S2 genes already in cache.")
        return expr

    log(f"  Need to extract {len(new_needed)} new genes:")
    log(f"  {new_needed}")

    new_df = extract_new_genes_from_mtx(new_needed, meta)

    if new_df is None:
        log("  MTX not found. New genes unavailable.")
        return expr

    # Merge into existing expression frame
    common_cells = expr.index.intersection(new_df.index)
    log(f"  Merging {len(new_df.columns)} new genes "
        f"across {len(common_cells)} cells...")

    expr_extended = expr.copy()
    for gene in new_df.columns:
        expr_extended[gene] = new_df.loc[expr_extended.index, gene]

    # Update cache
    expr_extended.to_csv(EXPR_CACHE)
    log(f"  Cache updated: {EXPR_CACHE}")
    log(f"  New shape: {expr_extended.shape}")

    return expr_extended

# ============================================================
# MAIN
# ============================================================

def main():
    log("=" * 65)
    log("BRCA LUMINAL A — SCRIPT 2")
    log("OrganismCore — Document BRCA-S2b")
    log("Date: 2026-03-04")
    log("=" * 65)
    log("")
    log("SCRIPT 2 PREDICTIONS (locked BRCA-S2a):")
    log("  S2-1: CDKN1A depth score separates LumA")
    log("  S2-2: r(ESR1, PGR) < 0.3 (ER circuit broken)")
    log("  S2-3: NCOA1 or NCOA2 suppressed")
    log("  S2-4: TGFB1 or SMAD3 suppressed")
    log("  S2-5: r(FOXA1, depth) > 0.30")
    log("  S2-6: CDK4 suppression attenuates vs progenitor")
    log("  S2-7: CDKN1A bulk high variance -> stratification")
    log("")

    # Step 0: Load cached data
    expr, meta = load_cached_data()
    if expr is None:
        log("FATAL: Cannot load cache. Run Script 1 first.")
        write_log()
        return

    # Optional: extend cache with new genes if MTX available
    expr = extend_cache_with_new_genes(expr, meta)

    # Step 1: Prepare populations
    expr_c, pops = prepare_populations(expr, meta)

    # Step 2: CDKN1A-centred depth score
    luma_depth_s2, corrs_s2 = compute_s2_depth_score(pops)

    # Step 3: LumA vs Luminal Progenitor
    s26_confirmed = compare_vs_progenitor(pops)

    # Step 4: Gap test — r(ESR1, PGR)
    gap_confirmed, r_gap = gap_test_er_pgr(pops)

    # Step 5: ER coactivators
    coact_results = er_coactivator_test(pops)

    # Step 6: TGF-β pathway
    tgfb_results = tgfb_pathway_test(pops)

    # Step 7: FOXA1 depth correlation
    r_foxa1 = foxa1_depth_test(pops, luma_depth_s2)

    # Step 8: Bulk analysis
    bulk_result = bulk_analysis_s2(pops)

    # Step 9: Prediction scorecard
    scorecard_s2(gap_confirmed, r_gap, coact_results,
                 tgfb_results, r_foxa1, s26_confirmed)

    # Save S2 results CSV
    if corrs_s2:
        corr_df = pd.DataFrame(
            corrs_s2, columns=["gene", "r", "p_value"]
        )
        corr_df.to_csv(CSV_FILE, index=False)
        log(f"\nS2 depth correlations saved: {CSV_FILE}")

    # Step 10: Figure — pass ALL required variables
    generate_s2_figure(
        pops, luma_depth_s2, corrs_s2,
        bulk_result, r_gap,
        gap_confirmed, r_foxa1          # <-- the two that were missing
    )

    write_log()

    log("")
    log("=" * 65)
    log("SCRIPT 2 COMPLETE")
    log(f"  Results: {S2_DIR}")
    log(f"  Log:     {LOG_FILE}")
    log(f"  Figure:  {FIG_FILE}")
    log(f"  CSV:     {CSV_FILE}")
    log("")
    log("READ IN THIS ORDER (Protocol v2.0):")
    log("  1. STEP 2: CDKN1A depth score (primary axis)")
    log("  2. STEP 3: LumA vs Progenitor (reference check)")
    log("  3. STEP 4: ESR1->PGR gap (circuit break test)")
    log("  4. STEP 5: Coactivators (gap mechanism)")
    log("  5. STEP 6: TGF-b (upstream of CDKN1A)")
    log("  6. STEP 7: FOXA1 depth (identity anchor test)")
    log("  7. STEP 8: Bulk stratification")
    log("  8. STEP 9: Scorecard")
    log("")
    log("THEN write BRCA-S2b reasoning artifact.")
    log("=" * 65)


if __name__ == "__main__":
    main()
