"""
BRCA LUMINAL B — SCRIPT 2
OrganismCore — Document BRCA-S5c
Deep Dive Series: Breast Cancer Subtypes

SCRIPT 1 FINDINGS (BRCA-S5b, 2026-03-05):

  ATTRACTOR REVISED:
    Type 3 dominant model REFUTED.
    Revised to: Type 1-L dominant
      (Locked Luminal — hyperproliferative with intact identity)

  CONFIRMED:
    FOXA1 EQUAL to LumA (-7.4%, EQUAL threshold)
    GATA3 +42.1% vs LumA  (HIGHER — not lower)
    ESR1  +64.4% vs LumA  (HIGHER — inverts gradient)
    CDKN1A -69.3% vs LumA (PRIMARY SIGNAL, p<1e-100)
    MKI67 +1487%, TOP2A +911%, CCND1 +126% vs LumA
    EZH2 FLAT vs LumA (ns) — EZH2i REFUTED
    DNMT3A +146.7% vs LumA (dominant epigenetic signal)
    HDAC2 +99.5% vs LumA, HDAC1 +74.6% vs LumA
    CD274 +209.6% vs LumA (novel immune signal)
    r(ESR1, PGR) = 0.234 in LumB vs 0.333 in LumA

  DEPTH SCORE MISSPECIFICATION (S5b):
    The identity component (1 - norm(FOXA1+GATA3)) was
    misspecified: LumB has HIGHER FOXA1/GATA3 than LumA,
    so the composite incorrectly penalises LumB identity.
    Script 2 uses CDKN1A-ONLY depth as primary.

  OPEN QUESTIONS (BRCA-S5b → BRCA-S5c):
    OQ-1: SPI1 +315% vs LumA — immune contamination?
    OQ-2: ESR1-high/PGR-low subpopulation within LumB
    OQ-3: CDKN1A-depleted subpopulation structure
    OQ-4: DNMT3A + HDAC1/2 convergence on CDKN1A
    OQ-5: ERBB2 quantification in LumB
    OQ-6: CD274/MYC co-expression
    OQ-7: CDKN2A/CDKN1A paradox (p16 up, p21 down)

SCRIPT 2 PREDICTIONS (locked BRCA-S5a/b → BRCA-S5c):
  SB2-1: CDKN1A depth separates LumB from all references
  SB2-2: r(ESR1, PGR) in LumB < LumA (confirmed partial,
          now test mechanistic node — NCOA1/SPDEF axis)
  SB2-3: DNMT3A co-expresses with HDAC2 within LumB
          (convergence hypothesis: r > 0.20)
  SB2-4: SPI1-expressing LumB cells co-express PTPRC
          (immune contamination, not epithelial ectopic)
  SB2-5: CDKN2A and CDKN1A are inversely correlated
          within LumB (senescence bypass signature)
  SB2-6: CD274 and MYC co-express within LumB (r > 0.25)
  SB2-7: ERBB2 in LumB exceeds Mature Luminal and LumA
  SB2-8: Bulk CDKN1A more variable than CCND1 across
          tumors (CDKN1A loss = primary event, not CCND1 gain)

NEW IN SCRIPT 2 vs SCRIPT 1:
  - Three normal references: Mature + Progenitor + LumA
  - Depth score: CDKN1A-only (identity misspec corrected)
  - Co-expression tests within LumB (OQ resolution)
  - SPI1/PTPRC immune test
  - CDKN2A/CDKN1A co-expression (p16/p21 paradox)
  - ER coactivators: NCOA1, NCOA2, SPDEF
  - ERBB2 formal quantification
  - Cross-subtype depth ordering (LumA, LumB, TNBC)
  - Bulk CDKN1A vs CCND1 variance comparison
  - Extended gene panel from S5b signals

Author: Eric Robert Lawson
Framework: OrganismCore
Protocol: Workflow_Protocol.md v2.0
Version: 1.0
Date: 2026-03-05
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
# CONFIGURATION — reuse Script 1 cache
# ============================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

S1_DIR     = os.path.join(SCRIPT_DIR, "lumb_s1_analysis")
S1_RESULTS = os.path.join(S1_DIR, "results")
S1_DATA    = os.path.join(S1_DIR, "data")
S2_DIR     = os.path.join(S1_RESULTS, "results_s2")

LOG_FILE   = os.path.join(S2_DIR, "lumb_s2_log.txt")
FIG_FILE   = os.path.join(S2_DIR, "lumb_s2_figure.png")
CSV_FILE   = os.path.join(S2_DIR, "lumb_s2_results.csv")

EXPR_CACHE = os.path.join(S1_RESULTS, "expr_cache_lumb_s1.csv")
BULK_FILE  = os.path.join(S1_DATA,
    "GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts.txt.gz")

SC_META_SEARCH = [
    os.path.join(S1_DATA, "Wu_etal_2021_BRCA_scRNASeq", "metadata.csv"),
    os.path.join(S1_DATA, "metadata.csv"),  # flat extraction fallback
]

os.makedirs(S2_DIR, exist_ok=True)

# ============================================================
# CELL POPULATION LABELS
# ============================================================

CANCER_LUMB   = "Cancer LumB SC"
CANCER_LUMA   = "Cancer LumA SC"
CANCER_BASAL  = "Cancer Basal SC"
MATURE_LUM    = "Mature Luminal"
LUMINAL_PROG  = "Luminal Progenitors"
CT_COL        = "celltype_subset"

# ============================================================
# GENE PANELS — SCRIPT 2 (EXTENDED FROM S5b)
# ============================================================

# Script 1 confirmed signals
CONFIRMED_S1 = [
    "FOXA1", "GATA3", "ESR1", "PGR", "SPDEF",
    "CDKN1A", "CDKN2A", "CDK4", "CCND1", "RB1",
    "MKI67", "TOP2A", "PCNA", "MYC",
    "EZH2", "HDAC1", "HDAC2", "KDM1A", "DNMT3A", "DNMT3B",
    "EED", "TET2", "TP53",
    "KRT5", "KRT14", "SOX10",
    "AR", "NCOA1", "NCOA2", "SPI1", "CD274",
    "ERBB2", "TGFBR2", "CDK6",
]

# SB2-2: ER coactivator / circuit gap panel
ER_CIRCUIT = ["NCOA1", "NCOA2", "NRIP1", "MED1", "EP300",
              "CREBBP", "SPDEF", "TFF1", "TFF3", "AGR2"]

# SB2-3: DNMT3A/HDAC convergence panel
EPIGENETIC_CONVERGENCE = [
    "DNMT3A", "DNMT3B", "HDAC1", "HDAC2",
    "KDM1A", "EZH2", "EED", "PARP1",
]

# SB2-4: SPI1/immune
IMMUNE_PANEL = ["SPI1", "PTPRC", "CD14", "LYZ", "CSF1R",
                "CD68", "ITGAM", "FCGR3A"]

# SB2-5/7: Cell cycle / p21/p16
CELL_CYCLE_S2 = [
    "CDKN1A", "CDKN2A", "CDKN1B",
    "CDK4", "CDK6", "CDK2",
    "CCND1", "CCNE1", "CCNA1", "CCNB1",
    "RB1", "E2F1", "E2F3",
]

# SB2-6: CD274/MYC immune escape
IMMUNE_ESCAPE = ["CD274", "MYC", "CD80", "PDCD1LG2"]

# SB2-7: ERBB2
ERBB2_PANEL = ["ERBB2", "ERBB3", "EGFR", "MUC4", "GRB7"]

# PIK3CA pathway
PIK3CA_GENES = ["PIK3CA", "AKT1", "AKT2", "MTOR", "PTEN"]

ALL_GENES_S2 = list(dict.fromkeys(
    CONFIRMED_S1 + ER_CIRCUIT +
    EPIGENETIC_CONVERGENCE + IMMUNE_PANEL +
    CELL_CYCLE_S2 + IMMUNE_ESCAPE +
    ERBB2_PANEL + PIK3CA_GENES
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
    if p < 1e-300: return "p=0.00e+00 ***"
    elif p < 0.001: return f"p={p:.2e} ***"
    elif p < 0.01:  return f"p={p:.2e}  **"
    elif p < 0.05:  return f"p={p:.4f}   *"
    else:           return f"p={p:.4f}  ns"


def norm01(series):
    mn, mx = series.min(), series.max()
    if mx - mn < 1e-9:
        return pd.Series(0.5, index=series.index)
    return (series - mn) / (mx - mn)


# ============================================================
# STEP 0: LOAD CACHED DATA FROM SCRIPT 1
# ============================================================

def load_cached_data():
    log("=" * 65)
    log("STEP 0: LOAD CACHED DATA FROM SCRIPT 1")
    log("Reusing Script 1 downloads. No re-download.")
    log("=" * 65)

    if not os.path.exists(EXPR_CACHE):
        log(f"  FATAL: Script 1 cache not found: {EXPR_CACHE}")
        log("  Run BRCA_LumB_script1.py first.")
        return None, None

    log(f"  Loading expression cache: {EXPR_CACHE}")
    expr = pd.read_csv(EXPR_CACHE, index_col=0)
    log(f"  Expression shape: {expr.shape}")
    log(f"  Genes in cache: {sorted(expr.columns.tolist())}")

    meta = None
    for path in SC_META_SEARCH:
        if os.path.exists(path):
            meta = pd.read_csv(path, index_col=0)
            log(f"  Metadata: {path} ({len(meta)} cells)")
            break

    if meta is None:
        log("  FATAL: metadata.csv not found.")
        return None, None

    missing = [g for g in ALL_GENES_S2 if g not in expr.columns]
    log(f"\n  New genes needed for S2 ({len(missing)}): {missing}")
    if missing:
        log("  These will be extracted from MTX if available.")

    return expr, meta


def extract_new_genes_from_mtx(new_genes, meta):
    """Extract additional genes from original MTX."""
    import scipy.io

    # Inside extract_new_genes_from_mtx() — replace SC_SEARCH:
    SC_SEARCH = [
        os.path.join(S1_DATA, "Wu_etal_2021_BRCA_scRNASeq"),
        os.path.join(S1_DATA),  # flat extraction fallback
    ]

    sc_dir = None
    for p in SC_SEARCH:
        if os.path.exists(os.path.join(p, "count_matrix_sparse.mtx")):
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
    barcodes   = pd.read_csv(bc_file, sep="\t",
                             header=None).iloc[:, 0].tolist()
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
            log(f"    {gene}: mean={np.log1p(row).mean():.4f}  "
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


def extend_cache_with_new_genes(expr, meta):
    """Extend cache with any S2 genes missing from S1."""
    log("=" * 65)
    log("OPTIONAL: EXTRACTING NEW GENES FROM MTX")
    log("=" * 65)

    new_needed = [g for g in ALL_GENES_S2 if g not in expr.columns]
    if not new_needed:
        log("  All S2 genes already in cache.")
        return expr

    log(f"  Need to extract {len(new_needed)} new genes:")
    log(f"  {new_needed}")

    new_df = extract_new_genes_from_mtx(new_needed, meta)
    if new_df is None:
        log("  MTX not found. New genes unavailable.")
        return expr

    common_cells = expr.index.intersection(new_df.index)
    log(f"  Merging {len(new_df.columns)} new genes "
        f"across {len(common_cells)} cells...")

    expr_extended = expr.copy()
    for gene in new_df.columns:
        expr_extended[gene] = new_df.loc[expr_extended.index, gene]

    expr_extended.to_csv(EXPR_CACHE)
    log(f"  Cache updated: {EXPR_CACHE}")
    log(f"  New shape: {expr_extended.shape}")
    return expr_extended


# ============================================================
# STEP 1: PREPARE POPULATIONS
# THREE normal references: Mature Luminal, Luminal Progenitor,
# and Cancer LumA SC (cross-subtype anchor)
# ============================================================

def prepare_populations(expr, meta):
    log("=" * 65)
    log("STEP 1: PREPARE POPULATIONS")
    log("THREE references: Mature Luminal + Luminal Progenitor"
        " + Cancer LumA SC")
    log("=" * 65)

    common = expr.index.intersection(meta.index)
    expr_c = expr.loc[common]
    meta_c = meta.loc[common]

    lumb_idx   = meta_c[meta_c[CT_COL] == CANCER_LUMB].index
    luma_idx   = meta_c[meta_c[CT_COL] == CANCER_LUMA].index
    basal_idx  = meta_c[meta_c[CT_COL] == CANCER_BASAL].index
    mature_idx = meta_c[meta_c[CT_COL] == MATURE_LUM].index
    prog_idx   = meta_c[meta_c[CT_COL] == LUMINAL_PROG].index

    log(f"  LumB SC:             n={len(lumb_idx)}")
    log(f"  LumA SC:             n={len(luma_idx)}")
    log(f"  Basal SC:            n={len(basal_idx)}")
    log(f"  Mature Luminal:      n={len(mature_idx)}")
    log(f"  Luminal Progenitors: n={len(prog_idx)}")

    pops = {
        "lumb":   expr_c.loc[lumb_idx],
        "luma":   expr_c.loc[luma_idx],
        "basal":  expr_c.loc[basal_idx],
        "mature": expr_c.loc[mature_idx],
        "prog":   expr_c.loc[prog_idx],
    }

    return expr_c, pops


# ============================================================
# STEP 2: CDKN1A-ONLY DEPTH SCORE
#
# S5b CORRECTION: The composite depth score that included
# 1 - norm(FOXA1+GATA3) was misspecified for LumB because
# FOXA1 and GATA3 are HIGHER in LumB than LumA. Using those
# genes in a 1-norm component caused them to pull the score
# toward LOW DEPTH (the opposite of biological reality for
# LumB cells with intact luminal identity).
#
# Script 2 uses CDKN1A-only depth as primary:
#   depth = 1 - norm(CDKN1A)   [joint normalisation]
#
# This is pure CDKN1A depletion = cell cycle exit failure
# = the defining biological axis of LumB vs LumA.
#
# Prediction SB2-1: CDKN1A depth separates LumB from all
# three references (Mature, Progenitor, LumA).
# ============================================================

def compute_s2_depth_score(pops):
    log("=" * 65)
    log("STEP 2: CDKN1A-ONLY DEPTH SCORE")
    log("Formula: depth = 1 - norm(CDKN1A)")
    log("S5b correction: identity component REMOVED")
    log("  (FOXA1/GATA3 are HIGHER in LumB — misspecified)")
    log("Joint normalisation across LumB + Mature + Progenitor + LumA")
    log("Prediction SB2-1: depth separates LumB from all references")
    log("=" * 65)

    lumb   = pops["lumb"]
    luma   = pops["luma"]
    mature = pops["mature"]
    prog   = pops["prog"]

    if "CDKN1A" not in lumb.columns:
        log("  ERROR: CDKN1A not in cache. Cannot compute depth.")
        return None, None

    # Joint normalise across all four populations
    all_cdkn1a = pd.concat([
        lumb["CDKN1A"],
        luma["CDKN1A"],
        mature["CDKN1A"],
        prog["CDKN1A"],
    ])
    mn_c1a = all_cdkn1a.min()
    mx_c1a = all_cdkn1a.max()

    def norm_cdkn1a(series):
        return (series - mn_c1a) / (mx_c1a - mn_c1a + 1e-9)

    lumb_depth   = 1 - norm_cdkn1a(lumb["CDKN1A"])
    luma_depth   = 1 - norm_cdkn1a(luma["CDKN1A"])
    mature_depth = 1 - norm_cdkn1a(mature["CDKN1A"])
    prog_depth   = 1 - norm_cdkn1a(prog["CDKN1A"])

    log(f"\n  CDKN1A-only depth scores:")
    log(f"  LumB:        mean={lumb_depth.mean():.4f}  "
        f"std={lumb_depth.std():.4f}")
    log(f"  LumA:        mean={luma_depth.mean():.4f}  "
        f"std={luma_depth.std():.4f}")
    log(f"  Mature Lum:  mean={mature_depth.mean():.4f}  "
        f"std={mature_depth.std():.4f}")
    log(f"  Luminal Prg: mean={prog_depth.mean():.4f}  "
        f"std={prog_depth.std():.4f}")

    # SB2-1 — separation tests
    for ref_name, ref_depth in [
        ("Mature Luminal", mature_depth),
        ("Luminal Progenitor", prog_depth),
        ("Cancer LumA", luma_depth),
    ]:
        try:
            _, p = stats.mannwhitneyu(
                lumb_depth, ref_depth,
                alternative="two-sided"
            )
        except Exception:
            p = 1.0
        direction = ("LumB > " if lumb_depth.mean() > ref_depth.mean()
                     else "LumB < ")
        log(f"\n  p(LumB vs {ref_name}): {fmt_p(p)}"
            f"  [{direction}{ref_name}]")

    # S5b S1 composite vs S2 CDKN1A-only depth
    s1_composite_genes = [g for g in
                          ["CDKN1A", "FOXA1", "GATA3", "MKI67"]
                          if g in lumb.columns]
    if len(s1_composite_genes) >= 2:
        # Reconstruct approximate S1 composite within LumB
        joint_tmp = pd.concat([lumb[s1_composite_genes],
                                mature[s1_composite_genes]])
        nrm = joint_tmp.apply(norm01, axis=0)
        # S1 composite used: 0.4*(1-CDKN1A) + 0.4*(1-FOXA1+GATA3 mean)
        #                     + 0.2*MKI67
        s1_approx = pd.Series(
            index=lumb.index, dtype=float
        )
        c1a_n = nrm.loc[lumb.index, "CDKN1A"]
        fg_n  = nrm.loc[lumb.index,
                        ["FOXA1", "GATA3"]].mean(axis=1) \
                if all(g in nrm.columns
                       for g in ["FOXA1", "GATA3"]) else None
        mki_n = nrm.loc[lumb.index, "MKI67"] \
                if "MKI67" in nrm.columns else None

        if fg_n is not None and mki_n is not None:
            s1_approx = (
                0.4 * (1 - c1a_n) +
                0.4 * (1 - fg_n) +
                0.2 * mki_n
            )
            r_s1s2, p_s1s2 = stats.pearsonr(
                s1_approx.values, lumb_depth.values
            )
            log(f"\n  S1-composite vs S2-CDKN1A depth correlation:")
            log(f"  r = {r_s1s2:+.4f}  {fmt_p(p_s1s2)}")
            if r_s1s2 > 0.7:
                log("  → High concordance: both capture same cell-cycle axis")
            elif r_s1s2 > 0.4:
                log("  → Partial concordance: S2 purer on CDKN1A axis")
            else:
                log("  → Low concordance: S2 is meaningfully different")

    # Depth correlations within LumB
    log(f"\n  Top S2 depth correlations within LumB:")
    log(f"  {'Gene':<14} {'r':>8}  {'p-value':>16}")
    log("  " + "-" * 44)

    corrs_s2 = []
    for gene in lumb.columns:
        if lumb[gene].std() < 1e-9:
            continue
        try:
            r, p = stats.pearsonr(
                lumb_depth.values, lumb[gene].values
            )
            if not np.isnan(r):
                corrs_s2.append((gene, r, p))
        except Exception:
            pass

    corrs_s2.sort(key=lambda x: abs(x[1]), reverse=True)
    for gene, r, p in corrs_s2[:25]:
        log(f"  {gene:<14} {r:>+8.4f}  {fmt_p(p)}")

    return lumb_depth, corrs_s2


# ============================================================
# STEP 3: THREE-REFERENCE COMPARISON
# LumB vs Mature / Progenitor / LumA on focused gene set
# ============================================================

def three_reference_comparison(pops):
    log("")
    log("=" * 65)
    log("STEP 3: THREE-REFERENCE COMPARISON")
    log("LumB vs Mature Luminal / Progenitor / Cancer LumA")
    log("=" * 65)

    lumb   = pops["lumb"]
    luma   = pops["luma"]
    mature = pops["mature"]
    prog   = pops["prog"]

    focus_genes = [g for g in [
        "FOXA1", "GATA3", "ESR1", "PGR", "SPDEF",
        "CDKN1A", "CDKN2A", "CDK4", "CCND1", "MKI67",
        "TOP2A", "MYC", "EZH2", "DNMT3A", "HDAC1", "HDAC2",
        "RB1", "TP53", "SPI1", "CD274", "ERBB2", "NCOA1",
    ] if g in lumb.columns]

    hdr = (f"\n  {'Gene':<12}  {'Mature':>8}  {'Prog':>8}  "
           f"{'LumA':>8}  {'LumB':>8}  "
           f"{'vsMat%':>8}  {'vsLumA%':>8}  p_luma")
    log(hdr)
    log("  " + "-" * 80)

    ref_results = {}
    for gene in focus_genes:
        mm   = mature[gene].mean() if gene in mature.columns else np.nan
        pm   = prog[gene].mean()   if gene in prog.columns   else np.nan
        lam  = luma[gene].mean()   if gene in luma.columns   else np.nan
        lbm  = lumb[gene].mean()

        pct_mat  = ((lbm - mm)  / mm  * 100) if mm  > 1e-6 else 0
        pct_luma = ((lbm - lam) / lam * 100) if lam > 1e-6 else 0

        try:
            _, p_luma = stats.mannwhitneyu(
                lumb[gene].values,
                luma[gene].values if gene in luma.columns
                else lumb[gene].values,
                alternative="two-sided"
            )
        except Exception:
            p_luma = 1.0

        ref_results[gene] = {
            "mature": mm, "prog": pm, "luma": lam, "lumb": lbm,
            "pct_mat": pct_mat, "pct_luma": pct_luma,
            "p_luma": p_luma,
        }

        log(f"  {gene:<12}  {mm:>8.4f}  {pm:>8.4f}  "
            f"{lam:>8.4f}  {lbm:>8.4f}  "
            f"{pct_mat:>+7.1f}%  {pct_luma:>+7.1f}%  "
            f"{fmt_p(p_luma)}")

    return ref_results


# ============================================================
# STEP 4: ER CIRCUIT GAP TEST
# SB2-2: r(ESR1, PGR) in LumB — mechanism of decoupling
# Test NCOA1/NCOA2 as the missing link
# ESR1 elevated +64% vs LumA but r(ESR1,PGR)=0.234 < LumA 0.333
# ============================================================

def er_circuit_gap_test(pops):
    log("")
    log("=" * 65)
    log("STEP 4: ER CIRCUIT GAP TEST — MECHANISM")
    log("SB2-2: ESR1 elevated but PGR output lower than expected")
    log("S5b: r(ESR1,PGR)=0.234 in LumB vs 0.333 in LumA")
    log("Test nodes: NCOA1, NCOA2, SPDEF, NRIP1")
    log("=" * 65)

    lumb   = pops["lumb"]
    luma   = pops["luma"]
    mature = pops["mature"]
    prog   = pops["prog"]

    # Reproduce S5b gap test result
    circuit_pairs = [
        ("ESR1", "PGR"),
        ("ESR1", "SPDEF"),
        ("FOXA1", "ESR1"),
        ("FOXA1", "PGR"),
        ("FOXA1", "GATA3"),
        ("ESR1", "NCOA1"),
        ("ESR1", "NCOA2"),
    ]

    log(f"\n  Within-population co-expression correlations:")
    log(f"  {'Pair':<22}  {'r_LumB':>8}  {'r_LumA':>8}  "
        f"{'r_Mature':>8}  Gap_LumB<LumA?")
    log("  " + "-" * 68)

    gap_results = {}
    for ga, gb in circuit_pairs:
        if ga not in lumb.columns or gb not in lumb.columns:
            log(f"  {ga}→{gb:<22} NOT AVAILABLE")
            continue

        r_lb = r_la = r_mt = np.nan
        try:
            r_lb, _ = stats.pearsonr(
                lumb[ga].values, lumb[gb].values)
        except Exception:
            pass
        try:
            r_la, _ = stats.pearsonr(
                luma[ga].values, luma[gb].values) \
                if ga in luma.columns and gb in luma.columns \
                else (np.nan, np.nan)
        except Exception:
            pass
        try:
            r_mt, _ = stats.pearsonr(
                mature[ga].values, mature[gb].values) \
                if ga in mature.columns and gb in mature.columns \
                else (np.nan, np.nan)
        except Exception:
            pass

        label = f"{ga}→{gb}"
        gap_direction = (
            "YES — LumB < LumA" if r_lb < r_la
            else "NO — LumB ≥ LumA"
        )
        log(f"  {label:<22}  {r_lb:>+8.4f}  {r_la:>+8.4f}  "
            f"{r_mt:>+8.4f}  {gap_direction}")
        gap_results[label] = {
            "r_lumb": r_lb, "r_luma": r_la, "r_mature": r_mt
        }

    # Coactivator availability test
    log(f"\n  ER coactivator levels in LumB vs LumA:")
    log(f"  {'Gene':<10}  {'LumA':>8}  {'LumB':>8}  "
        f"{'vs_LumA%':>10}  p_luma")
    log("  " + "-" * 55)

    ncoa_results = {}
    for gene in ER_CIRCUIT:
        if gene not in lumb.columns:
            continue
        lam = luma[gene].mean() if gene in luma.columns else np.nan
        lbm = lumb[gene].mean()
        pct = ((lbm - lam) / lam * 100) if lam > 1e-6 else 0
        try:
            _, p = stats.mannwhitneyu(
                lumb[gene].values,
                luma[gene].values if gene in luma.columns
                else lumb[gene].values,
                alternative="two-sided"
            )
        except Exception:
            p = 1.0
        result = (
            "SUPPRESSED" if pct < -20 and p < 0.05 else
            "ELEVATED"   if pct >  20 and p < 0.05 else
            "FLAT"
        )
        ncoa_results[gene] = result
        log(f"  {gene:<10}  {lam:>8.4f}  {lbm:>8.4f}  "
            f"{pct:>+9.1f}%  {fmt_p(p)}  → {result}")

    # SB2-2 verdict
    r_esr1_pgr_lb = gap_results.get("ESR1→PGR", {}).get("r_lumb", np.nan)
    r_esr1_pgr_la = gap_results.get("ESR1→PGR", {}).get("r_luma", np.nan)
    log(f"\n  SB2-2 VERDICT:")
    log(f"  r(ESR1, PGR) LumB={r_esr1_pgr_lb:+.4f} "
        f"vs LumA={r_esr1_pgr_la:+.4f}")
    ncoa1_s = ncoa_results.get("NCOA1", "UNKNOWN")
    ncoa2_s = ncoa_results.get("NCOA2", "UNKNOWN")
    spdef_s = ncoa_results.get("SPDEF", "UNKNOWN")
    log(f"  NCOA1: {ncoa1_s}  NCOA2: {ncoa2_s}  SPDEF: {spdef_s}")

    if r_esr1_pgr_lb < r_esr1_pgr_la:
        log("  ESR1→PGR decoupling confirmed in LumB vs LumA.")
        if ncoa1_s == "SUPPRESSED" or ncoa2_s == "SUPPRESSED":
            log("  Mechanism: NCOA1/2 coactivator gap confirmed. "
                "ESR1 cannot efficiently drive PGR without coactivator.")
            log("  SB2-2: ✓ CONFIRMED — coactivator gap identified")
        elif spdef_s == "SUPPRESSED":
            log("  NCOA1/2 not suppressed. Mechanism is post-FOXA1/GATA3: "
                "SPDEF as downstream decoupling node.")
            log("  SB2-2: ✓ PARTIAL — SPDEF suppressed; NCOA1/2 intact")
        else:
            log("  Mechanism not resolved at expression level. "
                "May be post-translational or epigenetic.")
            log("  SB2-2: ✓ DIRECTIONAL — decoupling confirmed, "
                "mechanism unresolved")
    else:
        log("  SB2-2: ✗ NOT CONFIRMED — r(ESR1,PGR) not lower in LumB")

    return gap_results, ncoa_results


# ============================================================
# STEP 5: EPIGENETIC CONVERGENCE TEST
# SB2-3: DNMT3A co-expresses with HDAC2 within LumB
# Hypothesis: DNMT3A + HDAC1/2 converge on CDKN1A suppression
# ============================================================

def epigenetic_convergence_test(pops):
    log("")
    log("=" * 65)
    log("STEP 5: EPIGENETIC CONVERGENCE TEST")
    log("SB2-3: DNMT3A co-expresses with HDAC2 within LumB (r > 0.20)")
    log("Tests: DNMT3A-HDAC1/2 co-expression")
    log("       DNMT3A-CDKN1A inverse correlation")
    log("       HDAC1/2-CDKN1A inverse correlation")
    log("=" * 65)

    lumb = pops["lumb"]
    luma = pops["luma"]

    epigen_pairs = [
        ("DNMT3A", "HDAC2"),
        ("DNMT3A", "HDAC1"),
        ("DNMT3A", "CDKN1A"),
        ("HDAC2",  "CDKN1A"),
        ("HDAC1",  "CDKN1A"),
        ("EZH2",   "CDKN1A"),
        ("KDM1A",  "CDKN1A"),
        ("DNMT3A", "EZH2"),
        ("HDAC1",  "HDAC2"),
        ("PARP1",  "CDKN1A"),
    ]

    log(f"\n  {'Pair':<22}  {'r_LumB':>8}  {'r_LumA':>8}  "
        f"p_LumB           SB2-3?")
    log("  " + "-" * 65)

    epigen_results = {}
    for ga, gb in epigen_pairs:
        if ga not in lumb.columns or gb not in lumb.columns:
            continue
        try:
            r_lb, p_lb = stats.pearsonr(
                lumb[ga].values, lumb[gb].values)
        except Exception:
            r_lb, p_lb = np.nan, np.nan
        try:
            r_la, _ = stats.pearsonr(
                luma[ga].values, luma[gb].values) \
                if ga in luma.columns and gb in luma.columns \
                else (np.nan, np.nan)
        except Exception:
            r_la = np.nan

        label = f"{ga}–{gb}"
        note = ""
        if ga == "DNMT3A" and gb == "HDAC2":
            note = ("SB2-3 ✓ CONFIRMED" if r_lb > 0.20 and p_lb < 0.05
                    else "SB2-3 ✗ NOT CONFIRMED")
        epigen_results[label] = {"r_lumb": r_lb, "r_luma": r_la}
        log(f"  {label:<22}  {r_lb:>+8.4f}  {r_la:>+8.4f}  "
            f"{fmt_p(p_lb):<16}  {note}")

    # Summary
    r_d3a_hdac2 = epigen_results.get(
        "DNMT3A–HDAC2", {}).get("r_lumb", np.nan)
    r_d3a_c1a   = epigen_results.get(
        "DNMT3A–CDKN1A", {}).get("r_lumb", np.nan)
    r_hd2_c1a   = epigen_results.get(
        "HDAC2–CDKN1A", {}).get("r_lumb", np.nan)

    log(f"\n  Convergence summary:")
    log(f"  r(DNMT3A, HDAC2)   = {r_d3a_hdac2:+.4f}  "
        f"(SB2-3 threshold > 0.20)")
    log(f"  r(DNMT3A, CDKN1A)  = {r_d3a_c1a:+.4f}  "
        f"(negative = DNMT3A suppresses CDKN1A)")
    log(f"  r(HDAC2,  CDKN1A)  = {r_hd2_c1a:+.4f}  "
        f"(negative = HDAC2 suppresses CDKN1A)")

    if (r_d3a_c1a < -0.10 and r_hd2_c1a < -0.10
            and r_d3a_hdac2 > 0.10):
        log("  CONVERGENCE HYPOTHESIS SUPPORTED:")
        log("  DNMT3A and HDAC2 both negatively correlate with CDKN1A")
        log("  and positively co-express with each other.")
        log("  Mechanistic interpretation: DNMT3A + HDAC2 co-operate")
        log("  to suppress CDKN1A in the most deeply CDKN1A-depleted LumB cells.")
    else:
        log("  Convergence hypothesis not clearly supported at expression level.")

    return epigen_results


# ============================================================
# STEP 6: SPI1 RESOLUTION TEST
# SB2-4: SPI1-expressing LumB cells co-express PTPRC/CD14/LYZ
# If SPI1 cells are myeloid contamination → high PTPRC/CD14
# If SPI1 is genuine epithelial → high FOXA1/GATA3 in same cell
# ============================================================

def spi1_resolution_test(pops):
    log("")
    log("=" * 65)
    log("STEP 6: SPI1 RESOLUTION TEST")
    log("SB2-4: SPI1 elevation — immune contamination or ectopic?")
    log("Test: r(SPI1, PTPRC) > |r(SPI1, FOXA1)| → contamination")
    log("      r(SPI1, FOXA1) > |r(SPI1, PTPRC)| → ectopic epithelial")
    log("=" * 65)

    lumb = pops["lumb"]
    luma = pops["luma"]

    spi1_available = "SPI1" in lumb.columns
    if not spi1_available:
        log("  SPI1 not in cache. SB2-4 cannot be evaluated.")
        log("  Extract SPI1 and PTPRC from MTX to resolve OQ-1.")
        return None

    log(f"\n  SPI1 in LumB: mean={lumb['SPI1'].mean():.4f}  "
        f"nonzero={(lumb['SPI1'] > 0).sum()} / {len(lumb)}")

    # Test co-expression pairs
    immune_genes = [g for g in IMMUNE_PANEL if g in lumb.columns]
    epithelial_genes = [g for g in
                        ["FOXA1", "GATA3", "ESR1", "CDKN1A", "KRT5"]
                        if g in lumb.columns]

    log(f"\n  r(SPI1, immune markers) in LumB:")
    r_spi1_immune = {}
    for gene in immune_genes:
        if gene == "SPI1":
            continue
        try:
            r, p = stats.pearsonr(
                lumb["SPI1"].values, lumb[gene].values)
            r_spi1_immune[gene] = r
            log(f"    SPI1–{gene:<12}  r={r:+.4f}  {fmt_p(p)}")
        except Exception:
            pass

    log(f"\n  r(SPI1, epithelial markers) in LumB:")
    r_spi1_epi = {}
    for gene in epithelial_genes:
        try:
            r, p = stats.pearsonr(
                lumb["SPI1"].values, lumb[gene].values)
            r_spi1_epi[gene] = r
            log(f"    SPI1–{gene:<12}  r={r:+.4f}  {fmt_p(p)}")
        except Exception:
            pass

    # Compare in LumA
    log(f"\n  r(SPI1, immune markers) in LumA (reference):")
    for gene in immune_genes:
        if gene == "SPI1" or gene not in luma.columns:
            continue
        if "SPI1" not in luma.columns:
            continue
        try:
            r, p = stats.pearsonr(
                luma["SPI1"].values, luma[gene].values)
            log(f"    SPI1–{gene:<12}  r={r:+.4f}  {fmt_p(p)}")
        except Exception:
            pass

    # Verdict
    ptprc_r = r_spi1_immune.get("PTPRC", np.nan)
    foxa1_r = r_spi1_epi.get("FOXA1", np.nan)
    cd14_r  = r_spi1_immune.get("CD14",  np.nan)

    log(f"\n  SB2-4 VERDICT:")
    log(f"  r(SPI1, PTPRC) = {ptprc_r:+.4f}")
    log(f"  r(SPI1, FOXA1) = {foxa1_r:+.4f}")
    log(f"  r(SPI1, CD14)  = {cd14_r:+.4f}")

    if not np.isnan(ptprc_r) and not np.isnan(foxa1_r):
        if ptprc_r > abs(foxa1_r):
            log("  SB2-4: ✓ CONFIRMED — SPI1 is immune contamination")
            log("  SPI1-expressing cells are myeloid, not epithelial LumB cells.")
            log("  The SPI1 elevation in S5b should be excluded from LumB biology.")
        elif foxa1_r > abs(ptprc_r):
            log("  SB2-4: ✗ NOT CONFIRMED — SPI1 co-expresses with FOXA1")
            log("  SPI1 elevation appears to be genuine epithelial expression.")
            log("  Mechanistic role of SPI1 in LumB requires investigation.")
        else:
            log("  SB2-4: INCONCLUSIVE — both correlations near zero.")
            log("  SPI1 expression is diffuse, not cell-type-associated.")
    else:
        log("  PTPRC or FOXA1 not available. Cannot resolve SPI1 biology.")

    return {
        "ptprc_r": ptprc_r,
        "foxa1_r": foxa1_r,
        "r_spi1_immune": r_spi1_immune,
        "r_spi1_epi": r_spi1_epi,
    }


# ============================================================
# STEP 7: CDKN2A / CDKN1A PARADOX
# SB2-5: CDKN2A and CDKN1A are INVERSELY correlated within LumB
# S5b finding: CDKN2A +93.9% vs LumA despite CDKN1A -69.3%
# Hypothesis: senescence bypass — p16 accumulates (failed
# senescence signal) while p21 is actively suppressed
# (senescence effector removed).
# Expected: r(CDKN2A, CDKN1A) < 0 within LumB
# ============================================================

def cdkn2a_cdkn1a_paradox(pops):
    log("")
    log("=" * 65)
    log("STEP 7: CDKN2A / CDKN1A PARADOX — SENESCENCE BYPASS")
    log("SB2-5: CDKN2A elevated +93.9% vs LumA, CDKN1A -69.3%")
    log("Prediction: r(CDKN2A, CDKN1A) < 0 within LumB")
    log("Interpretation: p16 accumulates as senescence signal;")
    log("p21 suppressed as the effector exit point.")
    log("=" * 65)

    lumb = pops["lumb"]
    luma = pops["luma"]

    if "CDKN2A" not in lumb.columns or "CDKN1A" not in lumb.columns:
        log("  CDKN2A or CDKN1A not in cache. SB2-5 cannot be evaluated.")
        return None

    try:
        r_lb, p_lb = stats.pearsonr(
            lumb["CDKN2A"].values, lumb["CDKN1A"].values)
        r_la, p_la = stats.pearsonr(
            luma["CDKN2A"].values, luma["CDKN1A"].values) \
            if "CDKN2A" in luma.columns and "CDKN1A" in luma.columns \
            else (np.nan, np.nan)
    except Exception as e:
        log(f"  Correlation error: {e}")
        return None

    log(f"\n  r(CDKN2A, CDKN1A) in LumB: {r_lb:+.4f}  {fmt_p(p_lb)}")
    log(f"  r(CDKN2A, CDKN1A) in LumA: {r_la:+.4f}  {fmt_p(p_la)}")

    # Distribution stats
    log(f"\n  CDKN2A in LumB: mean={lumb['CDKN2A'].mean():.4f}  "
        f"std={lumb['CDKN2A'].std():.4f}")
    log(f"  CDKN1A in LumB: mean={lumb['CDKN1A'].mean():.4f}  "
        f"std={lumb['CDKN1A'].std():.4f}")
    log(f"  CDKN2A in LumA: mean={luma['CDKN2A'].mean():.4f}  "
        f"std={luma['CDKN2A'].std():.4f}")
    log(f"  CDKN1A in LumA: mean={luma['CDKN1A'].mean():.4f}  "
        f"std={luma['CDKN1A'].std():.4f}")

    # Identify CDKN2A-high / CDKN1A-low cells
    cdkn2a_med = lumb["CDKN2A"].median()
    cdkn1a_med = lumb["CDKN1A"].median()

    bypass_cells = lumb[
        (lumb["CDKN2A"] > cdkn2a_med) &
        (lumb["CDKN1A"] < cdkn1a_med)
    ]
    n_bypass = len(bypass_cells)
    pct_bypass = n_bypass / len(lumb) * 100

    log(f"\n  CDKN2A-high / CDKN1A-low (senescence bypass) cells:")
    log(f"  n={n_bypass} / {len(lumb)} ({pct_bypass:.1f}%)")

    if bypass_cells.shape[0] > 10:
        for gene in ["CCND1", "MYC", "FOXA1", "GATA3", "ESR1",
                     "HDAC2", "DNMT3A"]:
            if gene in bypass_cells.columns:
                bp_m = bypass_cells[gene].mean()
                all_m = lumb[gene].mean()
                log(f"  Bypass vs all LumB: {gene} {bp_m:.4f} vs {all_m:.4f} "
                    f"({(bp_m - all_m) / all_m * 100:+.1f}%)")

    # SB2-5 verdict
    if r_lb < 0 and p_lb < 0.05:
        log(f"\n  SB2-5: ✓ CONFIRMED — r={r_lb:+.4f}")
        log("  CDKN2A and CDKN1A are inversely correlated within LumB.")
        log("  Senescence bypass signature confirmed:")
        log("  p16 (CDKN2A) accumulates; p21 (CDKN1A) is suppressed.")
        log("  LumB cells have activated the senescence program (p16+)")
        log("  but escaped it through CDKN1A depletion.")
    elif r_lb < 0:
        log(f"\n  SB2-5: ✓ DIRECTIONAL — r={r_lb:+.4f} (ns at p={p_lb:.4f})")
        log("  Trend toward inverse correlation present.")
    else:
        log(f"\n  SB2-5: ✗ NOT CONFIRMED — r={r_lb:+.4f}")
        log("  CDKN2A and CDKN1A do not inversely correlate in LumB.")
        log("  p16/p21 paradox is not a within-cell phenomenon.")
        log("  May reflect between-cell-state heterogeneity.")

    return {"r_lumb": r_lb, "p_lumb": p_lb,
            "r_luma": r_la, "n_bypass": n_bypass}


# ============================================================
# STEP 8: CD274 / MYC CO-EXPRESSION
# SB2-6: r(CD274, MYC) > 0.25 within LumB
# Hypothesis: MYC directly transactivates CD274 promoter
# in LumB. If confirmed, MYC-high LumB cells are the
# immune-evasive subpopulation.
# ============================================================

def cd274_myc_test(pops):
    log("")
    log("=" * 65)
    log("STEP 8: CD274 / MYC CO-EXPRESSION TEST")
    log("SB2-6: S5b — CD274 +209.6% vs LumA")
    log("Prediction: r(CD274, MYC) > 0.25 within LumB")
    log("Mechanism: MYC is known to directly activate CD274")
    log("           MYC-high LumB cells = immune-evasive state")
    log("=" * 65)

    lumb = pops["lumb"]
    luma = pops["luma"]

    if "CD274" not in lumb.columns or "MYC" not in lumb.columns:
        log("  CD274 or MYC not in cache. SB2-6 cannot be evaluated.")
        return None

    try:
        r_lb, p_lb = stats.pearsonr(
            lumb["CD274"].values, lumb["MYC"].values)
        r_la, p_la = stats.pearsonr(
            luma["CD274"].values, luma["MYC"].values) \
            if "CD274" in luma.columns and "MYC" in luma.columns \
            else (np.nan, np.nan)
    except Exception as e:
        log(f"  Correlation error: {e}")
        return None

    log(f"\n  r(CD274, MYC) in LumB: {r_lb:+.4f}  {fmt_p(p_lb)}")
    log(f"  r(CD274, MYC) in LumA: {r_la:+.4f}  {fmt_p(p_la)}")

    # Additional CD274 co-expression panel
    for gene in ["CCND1", "CDK4", "CDKN1A", "FOXA1",
                 "GATA3", "ESR1", "HDAC2", "DNMT3A"]:
        if gene not in lumb.columns:
            continue
        try:
            r, p = stats.pearsonr(
                lumb["CD274"].values, lumb[gene].values)
            log(f"  r(CD274, {gene:<10}) = {r:+.4f}  {fmt_p(p)}")
        except Exception:
            pass

    # Depth correlation
    try:
        _, p_cd274_cdkn1a = stats.pearsonr(
            lumb["CD274"].values, lumb["CDKN1A"].values)
        r_cd274_cdkn1a, _ = stats.pearsonr(
            lumb["CD274"].values, lumb["CDKN1A"].values)
        log(f"\n  r(CD274, CDKN1A) = {r_cd274_cdkn1a:+.4f}  "
            f"{fmt_p(p_cd274_cdkn1a)}")
        log("  (Negative: CD274 and CDKN1A suppress each other)")
    except Exception:
        pass

    # SB2-6 verdict
    if r_lb > 0.25 and p_lb < 0.05:
        log(f"\n  SB2-6: ✓ CONFIRMED — r={r_lb:+.4f}")
        log("  CD274 and MYC co-express within LumB.")
        log("  MYC-driven CD274 upregulation is the candidate mechanism.")
        log("  High-MYC LumB cells are the immune-evasive subpopulation.")
        log("  Rationale for CDK4/6i + anti-PD-L1 combination in LumB:")
        log("  CDK4/6 inhibition reduces MYC → may also reduce CD274.")
    elif r_lb > 0.10 and p_lb < 0.05:
        log(f"\n  SB2-6: ✓ DIRECTIONAL — r={r_lb:+.4f} (below 0.25 threshold)")
    else:
        log(f"\n  SB2-6: ✗ NOT CONFIRMED — r={r_lb:+.4f}")
        log("  CD274 elevation in LumB is not primarily MYC-driven.")

    return {"r_lumb": r_lb, "p_lumb": p_lb, "r_luma": r_la}


# ============================================================
# STEP 9: ERBB2 QUANTIFICATION
# SB2-7: ERBB2 in LumB > Mature Luminal and LumA
# S5b: ERBB2 +96.8% vs Progenitor (p=1.47e-56)
# This step formally quantifies ERBB2 vs all references and
# assesses fraction of LumB cells meeting HER2-low criteria.
# ============================================================

def erbb2_quantification(pops):
    log("")
    log("=" * 65)
    log("STEP 9: ERBB2 QUANTIFICATION")
    log("SB2-7: ERBB2 in LumB vs all references")
    log("S5b: ERBB2 +96.8% vs Progenitor (not in S1 panel)")
    log("Assess: HER2-low expression territory in LumB")
    log("=" * 65)

    lumb   = pops["lumb"]
    luma   = pops["luma"]
    mature = pops["mature"]
    prog   = pops["prog"]

    if "ERBB2" not in lumb.columns:
        log("  ERBB2 not in cache. SB2-7 cannot be evaluated.")
        log("  Extract ERBB2 from MTX to evaluate HER2-low territory.")
        return None

    log(f"\n  ERBB2 across populations:")
    erbb2_vals = {}
    for name, pop in [("Progenitor", prog), ("Mature", mature),
                      ("LumA", luma), ("LumB", lumb)]:
        if "ERBB2" in pop.columns:
            m = pop["ERBB2"].mean()
            s = pop["ERBB2"].std()
            erbb2_vals[name] = m
            log(f"  {name:<12}: mean={m:.4f}  std={s:.4f}")

    # Pairwise tests
    for ref_name, ref_pop in [("Mature", mature), ("LumA", luma)]:
        if "ERBB2" not in ref_pop.columns:
            continue
        try:
            _, p = stats.mannwhitneyu(
                lumb["ERBB2"].values, ref_pop["ERBB2"].values,
                alternative="two-sided"
            )
            lbm  = lumb["ERBB2"].mean()
            rfm  = ref_pop["ERBB2"].mean()
            pct  = (lbm - rfm) / rfm * 100 if rfm > 1e-6 else 0
            log(f"  LumB vs {ref_name}: {pct:+.1f}%  {fmt_p(p)}")
        except Exception:
            pass

    # HER2-low proxy: top quartile ERBB2 expressors within LumB
    erbb2_q75 = lumb["ERBB2"].quantile(0.75)
    erbb2_q50 = lumb["ERBB2"].quantile(0.50)
    erbb2_q90 = lumb["ERBB2"].quantile(0.90)
    log(f"\n  LumB ERBB2 distribution:")
    log(f"  Q50 (median) = {erbb2_q50:.4f}")
    log(f"  Q75          = {erbb2_q75:.4f}")
    log(f"  Q90          = {erbb2_q90:.4f}")

    hi_erbb2 = lumb[lumb["ERBB2"] > erbb2_q75]
    pct_hi = len(hi_erbb2) / len(lumb) * 100
    log(f"\n  ERBB2-Q75+ cells (HER2-low territory proxy): "
        f"n={len(hi_erbb2)} ({pct_hi:.1f}% of LumB)")

    if len(hi_erbb2) > 10:
        log(f"  Profile of ERBB2-high LumB cells:")
        for gene in ["CDKN1A", "MYC", "CCND1", "FOXA1", "GATA3",
                     "ESR1", "CD274"]:
            if gene in hi_erbb2.columns:
                hi_m = hi_erbb2[gene].mean()
                all_m = lumb[gene].mean()
                log(f"    {gene:<10}: hi={hi_m:.4f}  all={all_m:.4f}  "
                    f"({(hi_m - all_m) / all_m * 100:+.1f}%)")

    # SB2-7 verdict
    lumb_erbb2 = erbb2_vals.get("LumB", np.nan)
    luma_erbb2 = erbb2_vals.get("LumA", np.nan)
    mat_erbb2  = erbb2_vals.get("Mature", np.nan)

    log(f"\n  SB2-7 VERDICT:")
    if lumb_erbb2 > luma_erbb2 and lumb_erbb2 > mat_erbb2:
        log(f"  SB2-7: ✓ CONFIRMED — LumB ERBB2 > LumA and Mature")
        log(f"  LumB occupies HER2-intermediate expression territory.")
        log("  Clinically relevant: HER2-low criteria (IHC 1+ or 2+/ISH-)")
        log("  may be enriched in LumB relative to LumA.")
        log("  Supports ADC sensitivity (trastuzumab deruxtecan) in LumB.")
    elif lumb_erbb2 > luma_erbb2:
        log(f"  SB2-7: ✓ PARTIAL — LumB > LumA but not > Mature")
    else:
        log(f"  SB2-7: ✗ NOT CONFIRMED")

    return {"lumb_erbb2": lumb_erbb2, "luma_erbb2": luma_erbb2,
            "mat_erbb2": mat_erbb2}


# ============================================================
# STEP 10: BULK RNA-SEQ ANALYSIS
# SB2-8: CDKN1A more variable than CCND1 across tumors
# Hypothesis: CDKN1A loss is the primary initiating event;
# CCND1 gain is secondary/reactive.
# ============================================================

def bulk_analysis_s2(pops):
    log("")
    log("=" * 65)
    log("STEP 10: BULK RNA-SEQ — CDKN1A vs CCND1 VARIANCE")
    log("SB2-8: CDKN1A CV > CCND1 CV across tumors")
    log("Hypothesis: CDKN1A loss = primary LumB → LumA event")
    log("=" * 65)

    if not os.path.exists(BULK_FILE):
        log(f"  Bulk file not found: {BULK_FILE}")
        log("  Skipping bulk analysis.")
        return None

    log(f"  Loading: {BULK_FILE}")
    try:
        with gzip.open(BULK_FILE, "rt") as f:
            bulk = pd.read_csv(f, sep="\t", index_col=0)
    except Exception as e:
        log(f"  Bulk load error: {e}")
        return None

    log(f"  Bulk shape: {bulk.shape}")

    target_genes = [
        "CDKN1A", "CCND1", "CDK4", "MKI67", "TOP2A",
        "FOXA1", "GATA3", "ESR1", "PGR",
        "EZH2", "DNMT3A", "HDAC2", "MYC",
        "ERBB2", "RB1", "CDKN2A",
    ]
    avail = [g for g in target_genes if g in bulk.index]
    log(f"  Available target genes: {len(avail)} / {len(target_genes)}")

    bulk_sub = bulk.loc[avail]

    # Coefficient of variation for each gene
    log(f"\n  Gene CV across tumors (higher CV = more variable):")
    log(f"  {'Gene':<12}  {'Mean':>10}  {'Std':>10}  {'CV':>8}")
    log("  " + "-" * 44)

    cv_data = {}
    for gene in avail:
        vals = bulk_sub.loc[gene]
        m = vals.mean()
        s = vals.std()
        cv = s / m if m > 0 else 0
        cv_data[gene] = cv
        log(f"  {gene:<12}  {m:>10.1f}  {s:>10.1f}  {cv:>8.3f}")

    cv_cdkn1a = cv_data.get("CDKN1A", np.nan)
    cv_ccnd1  = cv_data.get("CCND1",  np.nan)

    log(f"\n  SB2-8 VERDICT:")
    log(f"  CV(CDKN1A) = {cv_cdkn1a:.3f}")
    log(f"  CV(CCND1)  = {cv_ccnd1:.3f}")

    if not np.isnan(cv_cdkn1a) and not np.isnan(cv_ccnd1):
        if cv_cdkn1a > cv_ccnd1:
            log("  SB2-8: ✓ CONFIRMED — CDKN1A more variable than CCND1")
            log("  Supports CDKN1A loss as primary event in LumB.")
            log("  CCND1 expression is more uniform — may be constitutive.")
            log("  Drug target implication: CDKN1A loss is the variable")
            log("  axis → tumors stratify by CDKN1A level, not CCND1.")
        else:
            log("  SB2-8: ✗ NOT CONFIRMED — CCND1 more variable")
            log("  CCND1 gain may be the primary axis, not CDKN1A loss.")

    # CDKN1A-depth proxy and correlations
    if "CDKN1A" in bulk_sub.index:
        cdkn1a_bulk = bulk_sub.loc["CDKN1A"]
        cdkn1a_norm = norm01(cdkn1a_bulk)
        bulk_depth = 1 - cdkn1a_norm

        log(f"\n  Bulk depth (1-CDKN1A) across {bulk.shape[1]} tumors:")
        log(f"  mean={bulk_depth.mean():.4f}  "
            f"std={bulk_depth.std():.4f}  "
            f"min={bulk_depth.min():.4f}  "
            f"max={bulk_depth.max():.4f}")

        sorted_d = bulk_depth.sort_values(ascending=False)
        log(f"\n  Highest depth tumors (lowest CDKN1A):")
        for tumor, score in sorted_d.head(5).items():
            log(f"    {tumor}: depth={score:.4f}  "
                f"CDKN1A={cdkn1a_bulk[tumor]:.1f}")

        log(f"\n  Lowest depth tumors (highest CDKN1A):")
        for tumor, score in sorted_d.tail(5).items():
            log(f"    {tumor}: depth={score:.4f}  "
                f"CDKN1A={cdkn1a_bulk[tumor]:.1f}")

        # Bulk correlations
        log(f"\n  Bulk correlations with CDKN1A:")
        log(f"  {'Gene':<12}  {'r':>8}  p-value")
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
        for gene, r, p in corrs:
            log(f"  {gene:<12}  {r:>+8.4f}  {fmt_p(p)}")

        return bulk_depth, bulk_sub, cv_data

    return None, None, cv_data


# ============================================================
# STEP 11: PREDICTION SCORECARD
# ============================================================

def scorecard_s2(depth_result, gap_results, ncoa_results,
                 epigen_results, spi1_results,
                 cdkn_paradox, cd274_results, erbb2_results,
                 cv_data):
    log("")
    log("=" * 65)
    log("STEP 11: SCRIPT 2 PREDICTION SCORECARD")
    log("=" * 65)
    log("")
    log("PREDICTIONS LOCKED IN BRCA-S5a/b → BRCA-S5c:")
    log("")

    def verdict(confirmed, label=""):
        if confirmed is True:    return "✓ CONFIRMED"
        if confirmed is False:   return "✗ NOT CONFIRMED"
        return "? CANNOT EVALUATE"

    # SB2-1
    log("SB2-1: CDKN1A depth separates LumB from all references")
    log("       → See STEP 2 output above")

    # SB2-2
    log("")
    log("SB2-2: r(ESR1, PGR) in LumB < LumA — NCOA/SPDEF mechanism")
    r_ep = gap_results.get("ESR1→PGR", {}) if gap_results else {}
    r_lb = r_ep.get("r_lumb", np.nan)
    r_la = r_ep.get("r_luma", np.nan)
    if not np.isnan(r_lb) and not np.isnan(r_la):
        sb22 = (r_lb < r_la)
        log(f"       r_LumB={r_lb:+.4f}  r_LumA={r_la:+.4f}  "
            f"→ {verdict(sb22)}")
        if ncoa_results:
            nc1 = ncoa_results.get("NCOA1", "?")
            nc2 = ncoa_results.get("NCOA2", "?")
            spd = ncoa_results.get("SPDEF", "?")
            log(f"       NCOA1={nc1}  NCOA2={nc2}  SPDEF={spd}")
    else:
        log(f"       → {verdict(None)}")

    # SB2-3
    log("")
    log("SB2-3: r(DNMT3A, HDAC2) > 0.20 within LumB")
    r_d3h2 = (epigen_results.get("DNMT3A–HDAC2", {}).get("r_lumb", np.nan)
              if epigen_results else np.nan)
    if not np.isnan(r_d3h2):
        sb23 = r_d3h2 > 0.20
        log(f"       r(DNMT3A, HDAC2) = {r_d3h2:+.4f}  "
            f"→ {verdict(sb23)}")
    else:
        log(f"       → {verdict(None)}")

    # SB2-4
    log("")
    log("SB2-4: SPI1 co-expresses with PTPRC (immune contamination)")
    if spi1_results:
        ptprc_r = spi1_results.get("ptprc_r", np.nan)
        foxa1_r = spi1_results.get("foxa1_r", np.nan)
        if not np.isnan(ptprc_r) and not np.isnan(foxa1_r):
            sb24 = ptprc_r > abs(foxa1_r)
            log(f"       r(SPI1,PTPRC)={ptprc_r:+.4f}  "
                f"r(SPI1,FOXA1)={foxa1_r:+.4f}  "
                f"→ {verdict(sb24)}")
        else:
            log(f"       → {verdict(None)}")
    else:
        log(f"       → {verdict(None)} (SPI1 not in cache)")

    # SB2-5
    log("")
    log("SB2-5: r(CDKN2A, CDKN1A) < 0 in LumB (senescence bypass)")
    if cdkn_paradox:
        r_c2a_c1a = cdkn_paradox.get("r_lumb", np.nan)
        p_c2a_c1a = cdkn_paradox.get("p_lumb", np.nan)
        if not np.isnan(r_c2a_c1a):
            sb25 = r_c2a_c1a < 0 and p_c2a_c1a < 0.05
            log(f"       r={r_c2a_c1a:+.4f}  {fmt_p(p_c2a_c1a)}  "
                f"→ {verdict(sb25)}")
    else:
        log(f"       → {verdict(None)}")

    # SB2-6
    log("")
    log("SB2-6: r(CD274, MYC) > 0.25 within LumB")
    if cd274_results:
        r_cd274_myc = cd274_results.get("r_lumb", np.nan)
        p_cd274_myc = cd274_results.get("p_lumb", np.nan)
        if not np.isnan(r_cd274_myc):
            sb26 = r_cd274_myc > 0.25 and p_cd274_myc < 0.05
            log(f"       r={r_cd274_myc:+.4f}  {fmt_p(p_cd274_myc)}  "
                f"→ {verdict(sb26)}")
    else:
        log(f"       → {verdict(None)}")

    # SB2-7
    log("")
    log("SB2-7: ERBB2 in LumB > Mature Luminal and LumA")
    if erbb2_results:
        lb = erbb2_results.get("lumb_erbb2", np.nan)
        la = erbb2_results.get("luma_erbb2", np.nan)
        mt = erbb2_results.get("mat_erbb2",  np.nan)
        if not np.isnan(lb):
            sb27 = lb > la and lb > mt
            log(f"       LumB={lb:.4f}  LumA={la:.4f}  Mature={mt:.4f}  "
                f"→ {verdict(sb27)}")
    else:
        log(f"       → {verdict(None)} (ERBB2 not in cache)")

    # SB2-8
    log("")
    log("SB2-8: CV(CDKN1A) > CV(CCND1) across bulk tumors")
    if cv_data:
        cv1 = cv_data.get("CDKN1A", np.nan)
        cv2 = cv_data.get("CCND1",  np.nan)
        if not np.isnan(cv1) and not np.isnan(cv2):
            sb28 = cv1 > cv2
            log(f"       CV(CDKN1A)={cv1:.3f}  CV(CCND1)={cv2:.3f}  "
                f"→ {verdict(sb28)}")
        else:
            log(f"       → {verdict(None)} (bulk data not available)")
    else:
        log(f"       → {verdict(None)}")


# ============================================================
# STEP 12: FIGURE — 9-PANEL S2 OUTPUT
# ============================================================

def generate_s2_figure(pops, lumb_depth, corrs_s2,
                       bulk_result, gap_results,
                       epigen_results, cdkn_paradox,
                       cd274_results, erbb2_results):
    log("")
    log("=" * 65)
    log("STEP 12: GENERATING SCRIPT 2 FIGURE")
    log("=" * 65)

    fig = plt.figure(figsize=(22, 16))
    fig.suptitle(
        "BRCA Luminal B — Script 2\n"
        "CDKN1A Depth | ER Circuit Gap | Epigenetic Convergence |"
        " Co-expression Tests\n"
        "OrganismCore | GSE176078 | 2026-03-05",
        fontsize=11, fontweight="bold", y=1.01
    )
    gs = gridspec.GridSpec(3, 3, figure=fig,
                           hspace=0.52, wspace=0.42)

    # Colour scheme: LumB = orange, LumA = blue, normal = grey/green
    clr = {
        "lumb":  "#e67e22",   # orange
        "luma":  "#2980b9",   # blue
        "norm":  "#95a5a6",   # grey (mature)
        "prog":  "#27ae60",   # green (progenitor)
        "up":    "#c0392b",   # red
        "down":  "#2980b9",   # blue
        "flat":  "#27ae60",   # green
    }

    def bar_color(pct):
        if pct > 20:  return clr["up"]
        if pct < -20: return clr["down"]
        return clr["flat"]

    lumb   = pops["lumb"]
    luma   = pops["luma"]
    mature = pops["mature"]
    prog   = pops["prog"]

    # ── Panel A: CDKN1A depth distribution (LumB vs LumA)
    ax_a = fig.add_subplot(gs[0, 0])
    if lumb_depth is not None:
        ax_a.hist(lumb_depth.values, bins=50,
                  color=clr["lumb"], alpha=0.65,
                  edgecolor="white", label="LumB SC")
        ax_a.axvline(lumb_depth.mean(),
                     color=clr["lumb"], lw=2, linestyle="--",
                     label=f"LumB mean={lumb_depth.mean():.3f}")

        # LumA depth on same axes
        if "CDKN1A" in luma.columns:
            all_c1a = pd.concat([lumb["CDKN1A"], luma["CDKN1A"],
                                  mature["CDKN1A"] if "CDKN1A" in mature.columns
                                  else pd.Series(),
                                  prog["CDKN1A"] if "CDKN1A" in prog.columns
                                  else pd.Series()])
            mn_c = all_c1a.min(); mx_c = all_c1a.max()
            luma_depth = 1 - (luma["CDKN1A"] - mn_c) / (mx_c - mn_c + 1e-9)
            ax_a.hist(luma_depth.values, bins=50,
                      color=clr["luma"], alpha=0.45,
                      edgecolor="white", label="LumA SC")
            ax_a.axvline(luma_depth.mean(),
                         color=clr["luma"], lw=2, linestyle="--",
                         label=f"LumA mean={luma_depth.mean():.3f}")

    ax_a.set_title("A — CDKN1A Depth: LumB vs LumA\n"
                   "(1 – CDKN1A, joint normalisation)",
                   fontsize=8, fontweight="bold")
    ax_a.set_xlabel("Depth Score", fontsize=7)
    ax_a.set_ylabel("n cells", fontsize=7)
    ax_a.legend(fontsize=6)

    # ── Panel B: Key genes across three references + LumB
    ax_b = fig.add_subplot(gs[0, 1])
    key_genes = [g for g in
                 ["CDKN1A", "CCND1", "FOXA1", "GATA3",
                  "ESR1", "PGR", "DNMT3A", "HDAC2"]
                 if g in lumb.columns]
    if key_genes:
        x = np.arange(len(key_genes))
        w = 0.22
        lbv = [lumb[g].mean() for g in key_genes]
        lav = [luma[g].mean() if g in luma.columns else 0
               for g in key_genes]
        mmv = [mature[g].mean() if g in mature.columns else 0
               for g in key_genes]
        pgv = [prog[g].mean() if g in prog.columns else 0
               for g in key_genes]
        ax_b.bar(x - 1.5*w, mmv, w, color=clr["norm"],
                 label="Mature Lum", alpha=0.85)
        ax_b.bar(x - 0.5*w, pgv, w, color=clr["prog"],
                 label="Lum Prog",   alpha=0.85)
        ax_b.bar(x + 0.5*w, lav, w, color=clr["luma"],
                 label="LumA SC",    alpha=0.85)
        ax_b.bar(x + 1.5*w, lbv, w, color=clr["lumb"],
                 label="LumB SC",    alpha=0.85)
        ax_b.set_xticks(x)
        ax_b.set_xticklabels(key_genes, rotation=45,
                              ha="right", fontsize=7)
        ax_b.set_title("B — Key Genes vs Three References\n"
                       "(Mature + Progenitor + LumA + LumB)",
                       fontsize=8, fontweight="bold")
        ax_b.set_ylabel("Expression (log1p)", fontsize=7)
        ax_b.legend(fontsize=5)

    # ── Panel C: ESR1→PGR scatter (LumB vs Mature)
    ax_c = fig.add_subplot(gs[0, 2])
    if ("ESR1" in lumb.columns and "PGR" in lumb.columns):
        ax_c.scatter(lumb["ESR1"].values, lumb["PGR"].values,
                     c=clr["lumb"], alpha=0.12, s=4,
                     label="LumB SC")
        if "ESR1" in luma.columns and "PGR" in luma.columns:
            ax_c.scatter(luma["ESR1"].values, luma["PGR"].values,
                         c=clr["luma"], alpha=0.10, s=3,
                         label="LumA SC")
        r_ep = np.nan
        if gap_results:
            r_ep = gap_results.get("ESR1→PGR", {}).get("r_lumb", np.nan)
        ax_c.set_xlabel("ESR1 expression", fontsize=7)
        ax_c.set_ylabel("PGR expression", fontsize=7)
        ax_c.set_title(f"C — ESR1 → PGR Circuit\n"
                       f"r_LumB={r_ep:+.4f}  (S5b: 0.234)",
                       fontsize=8, fontweight="bold")
        ax_c.legend(fontsize=6)

    # ── Panel D: Epigenetic convergence (DNMT3A vs HDAC2)
    ax_d = fig.add_subplot(gs[1, 0])
    if ("DNMT3A" in lumb.columns and "HDAC2" in lumb.columns):
        ax_d.scatter(lumb["DNMT3A"].values, lumb["HDAC2"].values,
                     c=clr["lumb"], alpha=0.15, s=4)
        r_d3h2 = np.nan
        if epigen_results:
            r_d3h2 = epigen_results.get(
                "DNMT3A–HDAC2", {}).get("r_lumb", np.nan)
        ax_d.set_xlabel("DNMT3A expression", fontsize=7)
        ax_d.set_ylabel("HDAC2 expression", fontsize=7)
        ax_d.set_title(f"D — DNMT3A vs HDAC2 in LumB\n"
                       f"r={r_d3h2:+.4f}  (SB2-3: > 0.20?)",
                       fontsize=8, fontweight="bold")
    else:
        ax_d.axis("off")
        ax_d.text(0.5, 0.5, "DNMT3A or HDAC2\nnot in cache",
                  ha="center", va="center",
                  transform=ax_d.transAxes,
                  fontsize=8, color="gray")
        ax_d.set_title("D — DNMT3A vs HDAC2\n(not available)",
                       fontsize=8, fontweight="bold")

    # ── Panel E: CDKN2A vs CDKN1A (senescence bypass)
    ax_e = fig.add_subplot(gs[1, 1])
    if ("CDKN2A" in lumb.columns and "CDKN1A" in lumb.columns):
        ax_e.scatter(lumb["CDKN2A"].values, lumb["CDKN1A"].values,
                     c=clr["lumb"], alpha=0.15, s=4)
        r_c2c1 = np.nan
        if cdkn_paradox:
            r_c2c1 = cdkn_paradox.get("r_lumb", np.nan)
        ax_e.set_xlabel("CDKN2A (p16) expression", fontsize=7)
        ax_e.set_ylabel("CDKN1A (p21) expression", fontsize=7)
        ax_e.set_title(f"E — CDKN2A vs CDKN1A\n"
                       f"r={r_c2c1:+.4f}  (SB2-5: < 0?)",
                       fontsize=8, fontweight="bold")
    else:
        ax_e.axis("off")
        ax_e.text(0.5, 0.5, "CDKN2A or CDKN1A\nnot in cache",
                  ha="center", va="center",
                  transform=ax_e.transAxes,
                  fontsize=8, color="gray")
        ax_e.set_title("E — CDKN2A vs CDKN1A\n(not available)",
                       fontsize=8, fontweight="bold")

    # ── Panel F: S2 depth correlations (top 20)
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
                           "Top 20 within LumB",
                           fontsize=8, fontweight="bold")
            ax_f.set_xlabel("r with CDKN1A depth", fontsize=7)
            ax_f.tick_params(axis="y", labelsize=6)

    # ── Panel G: Bulk CDKN1A depth
    ax_g = fig.add_subplot(gs[2, 0])
    if bulk_result is not None and bulk_result[0] is not None:
        bulk_depth, _, cv_data = bulk_result
        ax_g.bar(range(len(bulk_depth)),
                 bulk_depth.sort_values(ascending=False).values,
                 color=clr["lumb"], alpha=0.8, edgecolor="white")
        ax_g.axhline(bulk_depth.mean(),
                     color="black", lw=1.5, linestyle="--",
                     label=f"mean={bulk_depth.mean():.3f}")
        ax_g.set_title("G — Bulk Depth (1-CDKN1A)\n"
                       "Tumors sorted by depth",
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

    # ── Panel H: CD274 vs MYC (immune escape)
    ax_h = fig.add_subplot(gs[2, 1])
    if "CD274" in lumb.columns and "MYC" in lumb.columns:
        ax_h.scatter(lumb["MYC"].values, lumb["CD274"].values,
                     c=clr["lumb"], alpha=0.15, s=4)
        r_cm = np.nan
        if cd274_results:
            r_cm = cd274_results.get("r_lumb", np.nan)
        ax_h.set_xlabel("MYC expression", fontsize=7)
        ax_h.set_ylabel("CD274 (PD-L1) expression", fontsize=7)
        ax_h.set_title(f"H — MYC vs CD274 in LumB\n"
                       f"r={r_cm:+.4f}  (SB2-6: > 0.25?)",
                       fontsize=8, fontweight="bold")
    else:
        ax_h.axis("off")
        ax_h.text(0.5, 0.5, "CD274 or MYC\nnot in cache",
                  ha="center", va="center",
                  transform=ax_h.transAxes,
                  fontsize=8, color="gray")
        ax_h.set_title("H — MYC vs CD274\n(not available)",
                       fontsize=8, fontweight="bold")

    # ── Panel I: Summary
    ax_i = fig.add_subplot(gs[2, 2])
    ax_i.axis("off")

    r_d3h2_s = "N/A"
    r_c2c1_s = "N/A"
    r_cm_s   = "N/A"
    lb_erbb2 = "N/A"

    if epigen_results:
        v = epigen_results.get("DNMT3A–HDAC2", {}).get("r_lumb", np.nan)
        if not np.isnan(v):
            r_d3h2_s = f"{v:+.3f}"
    if cdkn_paradox:
        v = cdkn_paradox.get("r_lumb", np.nan)
        if not np.isnan(v):
            r_c2c1_s = f"{v:+.3f}"
    if cd274_results:
        v = cd274_results.get("r_lumb", np.nan)
        if not np.isnan(v):
            r_cm_s = f"{v:+.3f}"
    if erbb2_results:
        v = erbb2_results.get("lumb_erbb2", np.nan)
        if not np.isnan(v):
            lb_erbb2 = f"{v:.4f}"

    summary = (
        f"I — S2 SUMMARY\n"
        f"{'─'*30}\n"
        f"Dataset: GSE176078\n"
        f"LumB SC:        {len(lumb)}\n"
        f"LumA SC:        {len(luma)}\n"
        f"Mature Luminal: {len(mature)}\n"
        f"Luminal Prog:   {len(prog)}\n\n"
        f"SB2 PREDICTIONS:\n"
        f"SB2-1 CDKN1A depth:    see STEP 2\n"
        f"SB2-2 ESR1-PGR gap:    see STEP 4\n"
        f"SB2-3 DNMT3A-HDAC2:    r={r_d3h2_s}\n"
        f"SB2-4 SPI1 immune:     see STEP 6\n"
        f"SB2-5 CDKN2A-CDKN1A:  r={r_c2c1_s}\n"
        f"SB2-6 CD274-MYC:       r={r_cm_s}\n"
        f"SB2-7 ERBB2 elev.:     LumB={lb_erbb2}\n"
        f"SB2-8 CDKN1A CV:       see STEP 10\n\n"
        f"ATTRACTOR TYPE:\n"
        f"Type 1-L dominant\n"
        f"(Locked Luminal)\n\n"
        f"PRIMARY MECHANISM:\n"
        f"CDKN1A -69% vs LumA\n"
        f"CCND1 +126% vs LumA\n\n"
        f"DRUG TARGETS:\n"
        f"1. CDK4/6i (confirmed)\n"
        f"2. ET + CDK4/6i (backbone)\n"
        f"3. HDAC2 (investigational)\n"
        f"4. TOP2A / anthracyclines\n\n"
        f"Doc: BRCA-S5c | 2026-03-05"
    )

    ax_i.text(
        0.03, 0.97, summary,
        transform=ax_i.transAxes,
        fontsize=7,
        verticalalignment="top",
        fontfamily="monospace",
        bbox=dict(boxstyle="round",
                  facecolor="#fff7ee",
                  edgecolor="#e67e22")
    )

    plt.savefig(FIG_FILE, dpi=150, bbox_inches="tight")
    log(f"  Figure saved: {FIG_FILE}")
    plt.close()


# ============================================================
# MAIN
# ============================================================

def main():
    log("=" * 65)
    log("BRCA LUMINAL B — SCRIPT 2")
    log("OrganismCore — Document BRCA-S5c")
    log("Date: 2026-03-05")
    log("=" * 65)
    log("")
    log("SCRIPT 2 PREDICTIONS (locked BRCA-S5a/b → BRCA-S5c):")
    log("  SB2-1: CDKN1A depth separates LumB from all references")
    log("  SB2-2: r(ESR1, PGR) < r in LumA — NCOA/SPDEF mechanism")
    log("  SB2-3: r(DNMT3A, HDAC2) > 0.20 within LumB")
    log("  SB2-4: SPI1 co-expresses with PTPRC (immune, not ectopic)")
    log("  SB2-5: r(CDKN2A, CDKN1A) < 0 in LumB (senescence bypass)")
    log("  SB2-6: r(CD274, MYC) > 0.25 within LumB")
    log("  SB2-7: ERBB2 > Mature Luminal and LumA")
    log("  SB2-8: CV(CDKN1A) > CV(CCND1) — CDKN1A is primary variable")
    log("")

    # Step 0: Load cached data
    expr, meta = load_cached_data()
    if expr is None:
        log("FATAL: Cannot load cache. Run BRCA_LumB_script1.py first.")
        write_log()
        return

    # Extend cache with S2 genes if MTX available
    expr = extend_cache_with_new_genes(expr, meta)

    # Step 1: Prepare populations
    expr_c, pops = prepare_populations(expr, meta)

    # Step 2: CDKN1A-only depth score
    lumb_depth, corrs_s2 = compute_s2_depth_score(pops)

    # Step 3: Three-reference comparison
    ref_results = three_reference_comparison(pops)

    # Step 4: ER circuit gap test
    gap_results, ncoa_results = er_circuit_gap_test(pops)

    # Step 5: Epigenetic convergence
    epigen_results = epigenetic_convergence_test(pops)

    # Step 6: SPI1 resolution
    spi1_results = spi1_resolution_test(pops)

    # Step 7: CDKN2A / CDKN1A paradox
    cdkn_paradox = cdkn2a_cdkn1a_paradox(pops)

    # Step 8: CD274 / MYC co-expression
    cd274_results = cd274_myc_test(pops)

    # Step 9: ERBB2 quantification
    erbb2_results = erbb2_quantification(pops)

    # Step 10: Bulk analysis
    bulk_result = bulk_analysis_s2(pops)

    # Step 11: Prediction scorecard
    cv_data = bulk_result[2] if (bulk_result is not None and
                                  len(bulk_result) == 3) else None
    scorecard_s2(
        lumb_depth, gap_results, ncoa_results,
        epigen_results, spi1_results, cdkn_paradox,
        cd274_results, erbb2_results, cv_data
    )

    # Save results CSV
    if corrs_s2:
        corr_df = pd.DataFrame(
            corrs_s2, columns=["gene", "r", "p_value"]
        )
        corr_df.to_csv(CSV_FILE, index=False)
        log(f"\nS2 depth correlations saved: {CSV_FILE}")

    # Step 12: Generate figure
    generate_s2_figure(
        pops, lumb_depth, corrs_s2,
        bulk_result, gap_results,
        epigen_results, cdkn_paradox,
        cd274_results, erbb2_results
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
    log("  1. STEP 2:  CDKN1A depth (primary axis, corrected formula)")
    log("  2. STEP 3:  Three-reference comparison (geometry context)")
    log("  3. STEP 4:  ER circuit gap (ESR1/PGR/NCOA1 mechanism)")
    log("  4. STEP 5:  Epigenetic convergence (DNMT3A + HDAC2)")
    log("  5. STEP 6:  SPI1 resolution (immune vs ectopic)")
    log("  6. STEP 7:  CDKN2A/CDKN1A paradox (senescence bypass)")
    log("  7. STEP 8:  CD274/MYC co-expression (immune escape)")
    log("  8. STEP 9:  ERBB2 quantification (HER2-low territory)")
    log("  9. STEP 10: Bulk variance (CDKN1A vs CCND1 primary event)")
    log(" 10. STEP 11: Prediction scorecard")
    log("")
    log("THEN write BRCA-S5c reasoning artifact.")
    log("=" * 65)


if __name__ == "__main__":
    main()
