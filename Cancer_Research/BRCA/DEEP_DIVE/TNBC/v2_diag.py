"""
TNBC SCRIPT 2 — DATA DIAGNOSTIC
OrganismCore | 2026-03-04

Run this BEFORE Script 2.
It tests every URL and local path Script 2 depends on.
Produces a report of what is available and what must be downloaded.
No data is written except the report file.

Usage:
    python TNBC_s2_diagnostic.py
"""

import os
import sys
import gzip
import time
import urllib.request
import urllib.error

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
REPORT_FILE = os.path.join(SCRIPT_DIR, "tnbc_s2_diagnostic_report.txt")

lines = []
def log(msg=""):
    print(msg)
    lines.append(str(msg))

# ============================================================
# ALL CANDIDATE URLs TO TEST
# ============================================================

URLS_TO_TEST = [
    # ── GSE25066 series matrix ────────────────────────────
    (
        "GSE25066_SERIES_MATRIX",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/"
        "GSE25nnn/GSE25066/matrix/"
        "GSE25066_series_matrix.txt.gz",
    ),

    # ── GPL96 annotation ─────────────────────────────────
    # Candidate 1: soft file (large ~200MB but definitely works)
    (
        "GPL96_SOFT",
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/"
        "GPL96nnn/GPL96/soft/GPL96_family.soft.gz",
    ),
    # Candidate 2: annot file (small ~3MB but often 404)
    (
        "GPL96_ANNOT",
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/"
        "GPL96nnn/GPL96/annot/GPL96.annot.gz",
    ),
    # Candidate 3: miniml
    (
        "GPL96_MINIML",
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/"
        "GPL96nnn/GPL96/miniml/GPL96_family.xml.tgz",
    ),

    # ── TCGA-BRCA expression ─────────────────────────────
    # Candidate 1: pancanatlas hub (most reliable)
    (
        "TCGA_EXPR_PANCAN",
        "https://pancanatlas.xenahubs.net/download/"
        "TCGA.BRCA.sampleMap/HiSeqV2.gz",
    ),
    # Candidate 2: pancanatlas PANCAN normalized
    (
        "TCGA_EXPR_PANCAN_NORM",
        "https://pancanatlas.xenahubs.net/download/"
        "TCGA.BRCA.sampleMap/HiSeqV2_PANCAN.gz",
    ),
    # Candidate 3: Xena S3 (previously tried, often 403)
    (
        "TCGA_EXPR_XENA_S3",
        "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/"
        "download/TCGA.BRCA.sampleMap%2FHiSeqV2_PANCAN.gz",
    ),
    # Candidate 4: Xena public hub direct
    (
        "TCGA_EXPR_XENA_HUB",
        "https://tcga.xenahubs.net/download/"
        "TCGA.BRCA.sampleMap/HiSeqV2.gz",
    ),

    # ── TCGA-BRCA clinical ───────────────────────────────
    (
        "TCGA_CLIN_PANCAN",
        "https://pancanatlas.xenahubs.net/download/"
        "TCGA.BRCA.sampleMap/BRCA_clinicalMatrix",
    ),
    (
        "TCGA_CLIN_XENA_S3",
        "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/"
        "download/TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix",
    ),
    (
        "TCGA_CLIN_XENA_HUB",
        "https://tcga.xenahubs.net/download/"
        "TCGA.BRCA.sampleMap/BRCA_clinicalMatrix",
    ),

    # ── TCGA-BRCA phenotype / PAM50 ─────────────────────
    (
        "TCGA_PHENO_PANCAN",
        "https://pancanatlas.xenahubs.net/download/"
        "TCGA.BRCA.sampleMap/BRCA_phenotype_renamed_noPE.tsv.gz",
    ),
    (
        "TCGA_PHENO_XENA_S3",
        "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/"
        "download/TCGA.BRCA.sampleMap%2FBRCA_phenotype_renamed_noPE.tsv.gz",
    ),
    # Candidate: survival data from pancan
    (
        "TCGA_SURVIVAL_PANCAN",
        "https://pancanatlas.xenahubs.net/download/"
        "TCGA.BRCA.sampleMap/BRCA_survival_data",
    ),
    # Candidate: TCGA pan-cancer clinical from UCSC
    (
        "TCGA_PANCLIN_UCSC",
        "https://pancanatlas.xenahubs.net/download/"
        "Survival_SupplementalTable_S1_20171025_xena_sp",
    ),
    # Candidate: GDC TCGA BRCA survival from Xena GDC hub
    (
        "TCGA_GDC_SURVIVAL",
        "https://gdc.xenahubs.net/download/"
        "TCGA-BRCA.survival.tsv",
    ),
    (
        "TCGA_GDC_CLIN",
        "https://gdc.xenahubs.net/download/"
        "TCGA-BRCA.GDC_phenotype.tsv.gz",
    ),
    (
        "TCGA_GDC_EXPR",
        "https://gdc.xenahubs.net/download/"
        "TCGA-BRCA.htseq_fpkm.tsv.gz",
    ),
]

# ============================================================
# LOCAL PATH CHECKS
# ============================================================

# Paths to check for existing downloaded/cached data
LOCAL_PATHS_TO_CHECK = [
    # S1 cache (from TNBC Script 1 run)
    os.path.join(SCRIPT_DIR, "TNBC_s1_analysis", "results",
                 "expr_cache_tnbc_s1.csv"),
    os.path.join(SCRIPT_DIR, "TNBC_s1_analysis", "results",
                 "tnbc_s1_saddle.csv"),
    # LumA cache (may have been downloaded already)
    os.path.join(SCRIPT_DIR, "..", "LUMINAL_A", "luma_results",
                 "expr_cache_luma.csv"),
    # scRNA-seq source files
    os.path.join(SCRIPT_DIR, "..", "..", "Wu_etal_2021_BRCA_scRNASeq",
                 "count_matrix_sparse.mtx"),
    os.path.join(os.path.expanduser("~/cancer/BRCA"),
                 "Wu_etal_2021_BRCA_scRNASeq",
                 "count_matrix_sparse.mtx"),
    # Any previously downloaded GSE25066 or GPL96
    os.path.join(SCRIPT_DIR, "tnbc_s2_results", "data",
                 "GSE25066_series_matrix.txt.gz"),
    os.path.join(SCRIPT_DIR, "tnbc_s2_results", "data",
                 "GPL96_family.soft.gz"),
    os.path.join(SCRIPT_DIR, "tnbc_s2_results", "data",
                 "GPL96.annot.gz"),
]


def test_url(name, url, head_only=True, max_bytes=1024):
    """
    Test a URL. Returns (status, size_or_note).
    head_only=True: use HEAD request (fast, no download).
    If HEAD fails, falls back to a partial GET.
    """
    headers = {"User-Agent": "Mozilla/5.0"}
    try:
        req = urllib.request.Request(url, headers=headers, method="HEAD")
        with urllib.request.urlopen(req, timeout=20) as r:
            code    = r.getcode()
            clen    = r.headers.get("Content-Length", "unknown")
            ctype   = r.headers.get("Content-Type", "unknown")
            return "OK", f"HTTP {code}  size={clen}  type={ctype}"
    except urllib.error.HTTPError as e:
        # HEAD disallowed — try partial GET
        if e.code in (403, 405):
            try:
                req2 = urllib.request.Request(
                    url, headers={**headers, "Range": "bytes=0-1023"}
                )
                with urllib.request.urlopen(req2, timeout=20) as r2:
                    _ = r2.read(max_bytes)
                    return "OK_PARTIAL", f"Partial GET succeeded (HEAD blocked)"
            except urllib.error.HTTPError as e2:
                return "FAIL", f"HTTP {e2.code}: {e2.reason}"
            except Exception as e2:
                return "FAIL", f"GET error: {e2}"
        return "FAIL", f"HTTP {e.code}: {e.reason}"
    except urllib.error.URLError as e:
        return "FAIL", f"URL error: {e.reason}"
    except Exception as e:
        return "FAIL", f"Exception: {e}"


def check_local(path):
    abs_p = os.path.abspath(path)
    if os.path.exists(abs_p):
        sz = os.path.getsize(abs_p)
        return True, f"{sz/1e6:.1f} MB  [{abs_p}]"
    return False, f"NOT FOUND  [{abs_p}]"


def peek_gzip_lines(path, n=5):
    """Read first n lines of a gzip file to confirm content."""
    try:
        with gzip.open(path, "rt", errors="ignore") as f:
            lines_out = []
            for i, line in enumerate(f):
                if i >= n:
                    break
                lines_out.append(line.rstrip("\n")[:120])
        return lines_out
    except Exception as e:
        return [f"ERROR: {e}"]


# ============================================================
# MAIN DIAGNOSTIC
# ============================================================

def main():
    log("=" * 70)
    log("TNBC SCRIPT 2 — DATA DIAGNOSTIC")
    log("OrganismCore | 2026-03-04")
    log("=" * 70)
    log("")

    # ── SECTION 1: URL TESTS ─────────────────────────────
    log("=" * 70)
    log("SECTION 1: URL ACCESSIBILITY TESTS")
    log("Testing all candidate download URLs...")
    log("=" * 70)
    log("")

    url_results = {}
    for name, url in URLS_TO_TEST:
        log(f"  [{name}]")
        log(f"    {url}")
        status, note = test_url(name, url)
        url_results[name] = (status, url, note)
        symbol = "✓" if status.startswith("OK") else "✗"
        log(f"    {symbol} {status}: {note}")
        log("")
        time.sleep(0.3)  # be polite to servers

    # ── SECTION 2: LOCAL FILE CHECKS ─────────────────────
    log("=" * 70)
    log("SECTION 2: LOCAL FILE CHECKS")
    log("Checking for already-downloaded/cached data...")
    log("=" * 70)
    log("")

    local_results = {}
    for path in LOCAL_PATHS_TO_CHECK:
        found, note = check_local(path)
        symbol = "✓" if found else "✗"
        log(f"  {symbol} {note}")
        local_results[path] = (found, note)
        if found and path.endswith(".gz") and os.path.getsize(path) > 1000:
            lines_preview = peek_gzip_lines(path)
            log(f"    First lines:")
            for ln in lines_preview:
                log(f"      {ln}")
        log("")

    # ── SECTION 3: SUMMARY AND RECOMMENDATION ────────────
    log("=" * 70)
    log("SECTION 3: SUMMARY AND RECOMMENDATION")
    log("=" * 70)
    log("")

    ok_urls  = {k: v for k, v in url_results.items() if v[0].startswith("OK")}
    bad_urls = {k: v for k, v in url_results.items() if not v[0].startswith("OK")}

    log(f"  URLs accessible:     {len(ok_urls)} / {len(URLS_TO_TEST)}")
    log(f"  URLs failing:        {len(bad_urls)} / {len(URLS_TO_TEST)}")
    log("")

    log("  WORKING URLs:")
    for k, (status, url, note) in ok_urls.items():
        log(f"    ✓ {k}")
        log(f"      {url}")
    log("")
    log("  FAILING URLs:")
    for k, (status, url, note) in bad_urls.items():
        log(f"    ✗ {k}  — {note}")
    log("")

    # Derive recommendation
    log("  SCRIPT 2 CONFIGURATION RECOMMENDATION:")
    log("")

    # GSE25066
    if "GSE25066_SERIES_MATRIX" in ok_urls:
        log("  GSE25066:  ✓ Use NCBI GEO series matrix URL (confirmed working)")
    else:
        log("  GSE25066:  ✗ NCBI FTP blocked. Alternative: download manually from")
        log("               https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25066")

    # GPL96
    gpl96_ok = [k for k in ["GPL96_SOFT", "GPL96_ANNOT", "GPL96_MINIML"] if k in ok_urls]
    if gpl96_ok:
        log(f"  GPL96:     ✓ Use: {gpl96_ok[0]}")
        log(f"               {ok_urls[gpl96_ok[0]][1]}")
    else:
        log("  GPL96:     ✗ All GPL96 candidates failed.")
        log("               Manual download: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL96")

    # TCGA expression
    expr_ok = [k for k in [
        "TCGA_EXPR_PANCAN", "TCGA_EXPR_PANCAN_NORM",
        "TCGA_EXPR_XENA_HUB", "TCGA_EXPR_XENA_S3",
        "TCGA_GDC_EXPR"
    ] if k in ok_urls]
    if expr_ok:
        log(f"  TCGA EXPR: ✓ Use: {expr_ok[0]}")
        log(f"               {ok_urls[expr_ok[0]][1]}")
    else:
        log("  TCGA EXPR: ✗ All Xena candidates failed.")
        log("               Try GDC Data Portal: https://portal.gdc.cancer.gov/")

    # TCGA clinical/survival
    clin_ok = [k for k in [
        "TCGA_CLIN_PANCAN", "TCGA_CLIN_XENA_HUB", "TCGA_CLIN_XENA_S3",
        "TCGA_GDC_CLIN"
    ] if k in ok_urls]
    if clin_ok:
        log(f"  TCGA CLIN: ✓ Use: {clin_ok[0]}")
        log(f"               {ok_urls[clin_ok[0]][1]}")
    else:
        log("  TCGA CLIN: ✗ All clinical candidates failed.")

    surv_ok = [k for k in [
        "TCGA_GDC_SURVIVAL", "TCGA_SURVIVAL_PANCAN", "TCGA_PANCLIN_UCSC"
    ] if k in ok_urls]
    if surv_ok:
        log(f"  TCGA SURV: ✓ Use: {surv_ok[0]}")
        log(f"               {ok_urls[surv_ok[0]][1]}")
    else:
        log("  TCGA SURV: ✗ No survival URL confirmed. Will use clinical OS fields.")

    # PAM50
    pam50_ok = [k for k in [
        "TCGA_PHENO_PANCAN", "TCGA_PHENO_XENA_S3", "TCGA_GDC_CLIN"
    ] if k in ok_urls]
    if pam50_ok:
        log(f"  TCGA PAM50:✓ Use: {pam50_ok[0]}")
    else:
        log("  TCGA PAM50:✗ No PAM50 phenotype URL confirmed.")

    log("")
    log("  LOCAL CACHE STATUS:")
    s1_cache = os.path.join(
        SCRIPT_DIR, "TNBC_s1_analysis", "results", "expr_cache_tnbc_s1.csv"
    )
    found, note = check_local(s1_cache)
    if found:
        log(f"  ✓ Script 1 TNBC cache exists — Script 2 can reuse it")
        log(f"    {note}")
    else:
        log(f"  ✗ Script 1 TNBC cache NOT found at expected path")
        log(f"    Script 2 will need to rerun MTX loading OR")
        log(f"    set CACHE_FILE path manually in Script 2 config")
    log("")

    log("=" * 70)
    log("DIAGNOSTIC COMPLETE")
    log(f"Report saved: {REPORT_FILE}")
    log("=" * 70)

    with open(REPORT_FILE, "w") as f:
        f.write("\n".join(lines))


if __name__ == "__main__":
    main()
