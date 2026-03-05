"""
BRCA CLAUDIN-LOW — SCRIPT 4 (v6 — GSE96058 ID join fix)
External Dataset Validation — METABRIC + GSE96058
OrganismCore — Document BRCA-S7k | 2026-03-05

ROOT CAUSE OF "survival matched: 0":
  Clinical index = geoAcc (GSMxxxxxxx accession IDs)
  Expression columns = sample title strings
  These are TWO DIFFERENT ID namespaces.
  The join loop tried cidx[str(s)] where s is a title —
  it will never match a GSM accession.

  Fix: parse the "title" characteristic from the series
  matrix alongside geoAcc, build a title→geoAcc bridge,
  then join expression (by title) → bridge → clinical
  (by geoAcc).

  Confirmed from 12379Monty/GSE96058 getGSEData.Rmd:
    rownames(sampDesc) <- sampDesc$geoAcc
    colnames(geneExpression) ← sampDesc$title
  So: expression col == title; clinical row == geoAcc.

  The series matrix stores per-sample metadata as:
    !Sample_geo_accession   GSM001  GSM002  ...
    !Sample_title           title1  title2  ...
    !Sample_characteristics_ch1  "key: val" ...

  _parse_gse_matrix() already grabs geoAcc as the index
  and all characteristics as columns.
  We just also need to grab !Sample_title and add it as
  a column, then build the bridge dict:
    title_to_geoacc[title] = geoAcc

PREVIOUS FIXES RETAINED:
  v4: POST API parse fix (item["gene"]["hugoGeneSymbol"])
  v4: cached bad POST file deleted at startup
  v5: GSE96058 PAM50 has no claudin-low (N/A for P4/P8/P9)
  v5: _find_os() covers ovrallSurvDays / ovrallSurvEvent
  v5: build_sv() auto-converts days→months
"""

import os
import gzip
import json
import ssl
import time
import warnings
warnings.filterwarnings("ignore")

import urllib.request
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

# ── SSL ──────────────────────────────────────────────────────
try:
    import certifi
    SSL_CTX = ssl.create_default_context(cafile=certifi.where())
except ImportError:
    SSL_CTX = ssl.create_default_context()
SSL_CTX_NOVERIFY = ssl._create_unverified_context()

# ============================================================
# CONFIGURATION
# ============================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR   = os.path.join(SCRIPT_DIR, "Claudin_Low_s4_results")
DATA_DIR   = os.path.join(BASE_DIR, "data")
os.makedirs(BASE_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

LOG_FILE       = os.path.join(BASE_DIR, "cl_s4_log.txt")
FIG_FILE       = os.path.join(BASE_DIR, "cl_s4_figure.png")
SCORECARD_FILE = os.path.join(BASE_DIR, "cl_s4_scorecard.csv")

METABRIC_EXPR_FILE  = os.path.join(DATA_DIR, "metabric_expr.tsv")
METABRIC_CLIN_FILE  = os.path.join(DATA_DIR, "metabric_clinical.tsv")
METABRIC_POST_FILE  = os.path.join(DATA_DIR, "metabric_mrna_post.json")
METABRIC_COMBO_FILE = os.path.join(DATA_DIR, "METABRIC_combo.csv")
GSE96058_EXPR_FILE  = os.path.join(DATA_DIR, "GSE96058_expr.csv.gz")
GSE96058_MAT1_FILE  = os.path.join(DATA_DIR,
    "GSE96058-GPL11154_series_matrix.txt.gz")
GSE96058_MAT2_FILE  = os.path.join(DATA_DIR,
    "GSE96058-GPL18573_series_matrix.txt.gz")

CBIO_BASE = "https://www.cbioportal.org/api"

CL_NEG = ["CLDN3","CLDN4","CLDN7","CDH1","ESR1"]
CL_POS = ["VIM","CD44","SNAI1","ZEB1","FN1"]
CL_SIG = CL_NEG + CL_POS
IMMUNE_GENES = ["FOXP3","PDCD1","TIGIT","CD8A","GZMB",
                "PRF1","LAG3","CD4","IFNG"]
CT_GENES     = ["GAGE1","GAGE2D","GAGE4","GAGE12D",
                "GAGE12J","CT45A3","CT45A4","STRA8","DPPA2"]
MEMORY_GENES = ["FOXA1","SPDEF","GATA3"]
ALL_PANEL    = list(dict.fromkeys(
    CL_SIG + IMMUNE_GENES + CT_GENES + MEMORY_GENES + ["MKI67"]))

GENE_ENTREZ = {
    "CLDN3":   1365,  "CLDN4":   1364,  "CLDN7":   1366,
    "CDH1":     999,  "ESR1":    2099,  "VIM":     7431,
    "CD44":     960,  "SNAI1":   6615,  "ZEB1":    6935,
    "FN1":     2335,  "FOXP3":  50943,  "PDCD1":   5133,
    "TIGIT":  201633, "CD8A":     925,  "GZMB":    3002,
    "PRF1":    5551,  "LAG3":    3902,  "CD4":      920,
    "IFNG":    3458,  "GAGE1":   2543,  "GAGE2D": 729447,
    "GAGE4":   2576,  "CT45A3": 441521, "CT45A4": 728911,
    "STRA8":  90799,  "DPPA2":  151871, "FOXA1":   3169,
    "SPDEF":  25803,  "GATA3":   2625,  "MKI67":   4288,
    "GAGE12D": 26578, "GAGE12J": 26580,
}

ZENODO_METABRIC_URLS = [
    "https://zenodo.org/records/6320936/files/"
    "METABRIC_RNA_Mutation.csv?download=1",
    "https://zenodo.org/record/6320936/files/"
    "METABRIC_RNA_Mutation.csv?download=1",
    "https://zenodo.org/records/7689036/files/"
    "METABRIC_RNA_Mutation.csv?download=1",
    "https://raw.githubusercontent.com/malikyousef/"
    "METABRIC-RNA-Mutation/main/METABRIC_RNA_Mutation.csv",
    "https://raw.githubusercontent.com/aakarsh-2004/"
    "METABRIC-Breast-Cancer/main/METABRIC_RNA_Mutation.csv",
    "https://raw.githubusercontent.com/Subaru-Forester/"
    "METABRIC/master/METABRIC_RNA_Mutation.csv",
    "https://raw.githubusercontent.com/arnabghosh75/"
    "METABRIC/main/METABRIC_RNA_Mutation.csv",
]

GSE96058_EXPR_URLS = [
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/"
    "GSE96058/suppl/"
    "GSE96058_gene_expression_3273_samples_and_136_replicates"
    "_transformed.csv.gz",
    "https://www.ncbi.nlm.nih.gov/geo/download/"
    "?acc=GSE96058&format=file&file="
    "GSE96058%5Fgene%5Fexpression%5F3273%5Fsamples"
    "%5Fand%5F136%5Freplicates%5Ftransformed.csv.gz",
]
GSE96058_MAT1_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/"
    "GSE96058/matrix/GSE96058-GPL11154_series_matrix.txt.gz")
GSE96058_MAT2_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/"
    "GSE96058/matrix/GSE96058-GPL18573_series_matrix.txt.gz")

CL_THRESHOLD = 7
MIN_N   = 10
MIN_EVT = 5
TIMEOUT = 300
RETRIES = 3

# ============================================================
# LOGGING
# ============================================================

_log = []

def log(msg=""):
    print(msg)
    _log.append(str(msg))

def flush_log():
    with open(LOG_FILE, "w") as f:
        f.write("\n".join(_log))

# ============================================================
# HTTP UTILITIES
# ============================================================

def _urlopen(url, data=None, extra_headers=None):
    headers = {"User-Agent": "Mozilla/5.0", "Accept": "*/*"}
    if extra_headers:
        headers.update(extra_headers)
    req = urllib.request.Request(url, data=data,
                                 headers=headers)
    for ctx in [SSL_CTX, SSL_CTX_NOVERIFY]:
        try:
            return urllib.request.urlopen(req, context=ctx,
                                          timeout=TIMEOUT)
        except ssl.SSLError:
            continue
        except Exception:
            raise
    return urllib.request.urlopen(req, timeout=TIMEOUT)

def fetch_url(url, dest, label="", retries=RETRIES,
              force=False):
    if (not force and os.path.exists(dest) and
            os.path.getsize(dest) > 500):
        log(f"  {label}: cached ({os.path.getsize(dest):,} B)")
        return True
    for attempt in range(1, retries + 1):
        try:
            log(f"  {label}: attempt {attempt} — "
                f"{url[:70]}...")
            with _urlopen(url) as r:
                data = r.read()
            with open(dest, "wb") as f:
                f.write(data)
            log(f"  {label}: saved {len(data):,} B")
            return True
        except Exception as e:
            log(f"  {label}: attempt {attempt} failed — {e}")
            time.sleep(min(2 ** attempt, 20))
    log(f"  {label}: ALL attempts failed")
    return False

def fetch_any(urls, dest, label="", force=False):
    for url in urls:
        if fetch_url(url, dest, label, force=force):
            return True
    return False

def post_json(url, body, dest, label="",
              retries=RETRIES, force=False):
    if (not force and os.path.exists(dest) and
            os.path.getsize(dest) > 100):
        log(f"  {label}: cached")
        return True
    data = json.dumps(body).encode("utf-8")
    for attempt in range(1, retries + 1):
        try:
            log(f"  {label}: POST attempt {attempt}")
            with _urlopen(url, data=data,
                          extra_headers={
                              "Content-Type": "application/json",
                              "Accept": "application/json"
                          }) as r:
                resp = r.read()
            with open(dest, "wb") as f:
                f.write(resp)
            log(f"  {label}: saved {len(resp):,} B")
            return True
        except Exception as e:
            log(f"  {label}: attempt {attempt} failed — {e}")
            time.sleep(min(2 ** attempt, 20))
    log(f"  {label}: ALL POST attempts failed")
    return False

# ============================================================
# METABRIC ACQUISITION
# ============================================================

def discover_metabric_profile():
    log(""); log("=" * 65)
    log("STEP 1: DISCOVER METABRIC MOLECULAR PROFILE")
    log("=" * 65)
    pf  = os.path.join(DATA_DIR, "metabric_profiles.json")
    url = (f"{CBIO_BASE}/studies/brca_metabric/"
           f"molecular-profiles?pageSize=100")
    if fetch_url(url, pf, "METABRIC profiles"):
        try:
            with open(pf) as f:
                ps = json.load(f)
            ids = [p.get("molecularProfileId", "")
                   for p in ps if isinstance(p, dict)]
            log(f"  Profiles: {ids}")
            for pid in ids:
                pl = pid.lower()
                if "mrna" in pl and (
                        "zscore" in pl or "median" in pl or
                        "z_score" in pl):
                    log(f"  Selected: {pid}"); return pid
            for pid in ids:
                if "mrna" in pid.lower():
                    log(f"  Fallback: {pid}"); return pid
            if ids: return ids[0]
        except Exception as e:
            log(f"  Profile parse error: {e}")
    return "brca_metabric_mrna_median_Zscores"

def get_metabric_samples():
    log(""); log("=" * 65)
    log("STEP 2: GET METABRIC SAMPLE IDs")
    log("=" * 65)
    sf  = os.path.join(DATA_DIR, "metabric_samples.json")
    url = (f"{CBIO_BASE}/studies/brca_metabric/samples"
           f"?pageSize=5000")
    if fetch_url(url, sf, "METABRIC samples"):
        try:
            with open(sf) as f:
                samps = json.load(f)
            ids = [s["sampleId"] for s in samps
                   if isinstance(s, dict) and "sampleId" in s]
            log(f"  Sample count: {len(ids)}")
            return ids
        except Exception as e:
            log(f"  Sample parse error: {e}")
    return []

def acquire_metabric_expr_post(profile_id, sample_ids):
    log(""); log("=" * 65)
    log(f"STEP 3: FETCH METABRIC EXPRESSION VIA POST")
    log(f"  profile: {profile_id}")
    log("=" * 65)
    entrez_ids     = [GENE_ENTREZ[g] for g in ALL_PANEL
                      if g in GENE_ENTREZ]
    entrez_to_hugo = {v: k for k, v in GENE_ENTREZ.items()}
    log(f"  Requesting {len(entrez_ids)} genes, "
        f"{len(sample_ids)} samples")
    url = (f"{CBIO_BASE}/molecular-profiles/{profile_id}/"
           f"molecular-data/fetch"
           f"?projection=SUMMARY&pageSize=500000000"
           f"&pageNumber=0")
    variants = [
        ({"entrezGeneIds": entrez_ids,
          "sampleListId": "brca_metabric_all"},
         "sampleListId"),
        ({"entrezGeneIds": entrez_ids,
          "sampleIds": sample_ids[:2500]},
         "sampleIds_2500"),
        ({"entrezGeneIds": entrez_ids,
          "sampleIds": sample_ids},
         "sampleIds_all"),
    ]
    for body, tag in variants:
        dest = os.path.join(DATA_DIR,
                            f"metabric_post_{tag}.json")
        if not post_json(url, body, dest,
                         f"METABRIC POST {tag}"):
            continue
        try:
            with open(dest) as f:
                raw = json.load(f)
            if not raw:
                log(f"  {tag}: empty response"); continue
            log(f"  {tag}: {len(raw)} data points")
            records = []
            for item in raw:
                gene_obj = item.get("gene", {})
                symbol   = (gene_obj.get("hugoGeneSymbol") or
                             item.get("hugoGeneSymbol") or
                             entrez_to_hugo.get(
                                 item.get("entrezGeneId"), ""))
                if not symbol: continue
                records.append({
                    "gene":   symbol.upper(),
                    "sample": item.get("sampleId", ""),
                    "value":  item.get("value", np.nan),
                })
            if not records:
                log(f"  {tag}: 0 records after parse"); continue
            df    = pd.DataFrame(records)
            pivot = df.pivot_table(
                index="gene", columns="sample",
                values="value", aggfunc="first")
            log(f"  {tag}: matrix {pivot.shape}")
            sig_ok = len([g for g in CL_SIG
                          if g in pivot.index])
            log(f"  {tag}: sig genes = {sig_ok}/10")
            if sig_ok < 4:
                log(f"  {tag}: too few — skip"); continue
            pivot.to_csv(METABRIC_EXPR_FILE, sep="\t")
            log(f"  Saved: {METABRIC_EXPR_FILE}")
            return True
        except Exception as e:
            log(f"  {tag}: parse failed — {e}")
            import traceback; log(traceback.format_exc())
    return False

def acquire_metabric_combo():
    log(""); log("=" * 65)
    log("STEP 4: ACQUIRE METABRIC COMBO CSV (FALLBACK)")
    log("=" * 65)
    if (os.path.exists(METABRIC_COMBO_FILE) and
            os.path.getsize(METABRIC_COMBO_FILE) > 100000):
        log(f"  Cached: {METABRIC_COMBO_FILE}"); return True
    ok = fetch_any(ZENODO_METABRIC_URLS,
                   METABRIC_COMBO_FILE, "METABRIC combo CSV")
    if ok:
        try:
            peek = pd.read_csv(METABRIC_COMBO_FILE,
                               nrows=3, low_memory=False)
            log(f"  CSV peek: {peek.shape[1]} columns")
            if peek.shape[1] < 5:
                os.remove(METABRIC_COMBO_FILE); return False
            return True
        except Exception as e:
            log(f"  CSV sanity check failed: {e}")
            if os.path.exists(METABRIC_COMBO_FILE):
                os.remove(METABRIC_COMBO_FILE)
    return False

def acquire_metabric_clinical():
    log(""); log("=" * 65)
    log("STEP 5: ACQUIRE METABRIC CLINICAL")
    log("=" * 65)
    if (os.path.exists(METABRIC_CLIN_FILE) and
            os.path.getsize(METABRIC_CLIN_FILE) > 1000):
        log("  Cached"); return True
    url  = (f"{CBIO_BASE}/studies/brca_metabric/clinical-data"
            f"?clinicalDataType=PATIENT&pageSize=100000")
    dest = os.path.join(DATA_DIR, "metabric_clin_api.json")
    if fetch_url(url, dest, "METABRIC clinical API"):
        try:
            with open(dest) as f: raw = json.load(f)
            records = {}
            for item in raw:
                pid = item.get("patientId", "")
                key = item.get("clinicalAttributeId", "")
                val = item.get("value", "")
                records.setdefault(pid, {})[key] = val
            clin = pd.DataFrame.from_dict(records,
                                           orient="index")
            clin.to_csv(METABRIC_CLIN_FILE, sep="\t")
            log(f"  Clinical: {clin.shape}")
            return True
        except Exception as e:
            log(f"  Clinical parse error: {e}")
    return False

def acquire_gse96058():
    log(""); log("=" * 65)
    log("STEP 6: ACQUIRE GSE96058")
    log("=" * 65)
    expr_ok = fetch_any(GSE96058_EXPR_URLS,
                        GSE96058_EXPR_FILE, "GSE96058 expr")
    mat1_ok = fetch_url(GSE96058_MAT1_URL,
                        GSE96058_MAT1_FILE, "GSE96058 mat1")
    mat2_ok = fetch_url(GSE96058_MAT2_URL,
                        GSE96058_MAT2_FILE, "GSE96058 mat2")
    return expr_ok, mat1_ok or mat2_ok

# ============================================================
# LOAD FUNCTIONS
# ============================================================

def load_metabric_expr():
    log(""); log("=" * 65)
    log("LOAD: METABRIC EXPRESSION"); log("=" * 65)
    all_g_up = {g.upper() for g in ALL_PANEL}
    for cand, default_sep in [(METABRIC_EXPR_FILE,  "\t"),
                               (METABRIC_COMBO_FILE, ",")]:
        if not os.path.exists(cand): continue
        try:
            sep = "\t" if cand.endswith(".tsv") else ","
            tmp = pd.read_csv(cand, sep=sep, index_col=0,
                              low_memory=False, comment="#")
            if "Entrez_Gene_Id" in tmp.columns:
                tmp.drop(columns=["Entrez_Gene_Id"],
                         inplace=True)
            log(f"  Loaded {os.path.basename(cand)}: "
                f"{tmp.shape}")
            tmp.index = (tmp.index.astype(str)
                         .str.strip().str.upper())
            row_ov = len(all_g_up & set(tmp.index))
            col_ov = len(all_g_up & set(
                tmp.columns.astype(str).str.upper()))
            log(f"  Gene overlap — rows:{row_ov} cols:{col_ov}")
            if col_ov > row_ov:
                gene_map = {c: str(c).strip().upper()
                            for c in tmp.columns
                            if str(c).strip().upper() in all_g_up}
                tmp = tmp[list(gene_map.keys())].rename(
                    columns=gene_map).T
                tmp.index = tmp.index.str.upper()
                log(f"  Transposed: {tmp.shape}")
            sig_ok = len([g for g in CL_SIG if g in tmp.index])
            log(f"  Sig genes: {sig_ok}/10")
            if sig_ok < 4:
                log("  Too few — skipping"); continue
            tmp = tmp.apply(pd.to_numeric, errors="coerce")
            log(f"  Final shape: {tmp.shape}")
            for g in CL_SIG + MEMORY_GENES + ["MKI67"]:
                log(f"    {g:<10} "
                    f"{'OK' if g in tmp.index else '--'}")
            return tmp
        except Exception as e:
            log(f"  Failed: {e}")
            import traceback; log(traceback.format_exc())
    log("  FATAL: Could not load METABRIC expression")
    return None

def load_metabric_clinical():
    log(""); log("=" * 65)
    log("LOAD: METABRIC CLINICAL"); log("=" * 65)
    for cand, sep in [
        (METABRIC_CLIN_FILE, "\t"),
        (METABRIC_COMBO_FILE, ","),
        (os.path.join(DATA_DIR, "metabric_clin_api.json"), None),
    ]:
        if not os.path.exists(cand): continue
        try:
            if cand.endswith(".json"):
                with open(cand) as f: raw = json.load(f)
                records = {}
                for item in raw:
                    pid = item.get("patientId", "")
                    key = item.get("clinicalAttributeId", "")
                    val = item.get("value", "")
                    records.setdefault(pid, {})[key] = val
                clin = pd.DataFrame.from_dict(records,
                                               orient="index")
            else:
                clin = pd.read_csv(cand, sep=sep,
                                   index_col=0,
                                   low_memory=False,
                                   comment="#")
            log(f"  Loaded {os.path.basename(cand)}: "
                f"{clin.shape}")
            log(f"  Columns: {list(clin.columns)[:20]}")
            ot, oe = _find_os(clin)
            if ot and oe:
                log(f"  OS: time={ot}  event={oe}")
                return clin, ot, oe
            log("  OS columns not found in this file")
        except Exception as e:
            log(f"  Failed: {e}")
    return None, None, None

# ============================================================
# GSE96058 SERIES MATRIX PARSER  ← KEY CHANGE IN v6
#
# Parses BOTH the sample characteristics AND !Sample_title.
# Returns:
#   meta  — DataFrame indexed by geoAcc (GSMxxxxxx)
#           with characteristics as columns
#           PLUS a "title" column
#   title_to_geoacc — dict mapping title → geoAcc
#
# This bridge is required because:
#   expression matrix columns  = sample titles
#   clinical DataFrame index   = geoAcc
# ============================================================

def _parse_gse_matrix(mat_file):
    """
    Parse a GEO series matrix file.
    Returns (meta_df, title_to_geoacc_dict).
    meta_df is indexed by geoAcc with clinical columns.
    title_to_geoacc maps sample title → geoAcc.
    """
    sample_ids  = []   # geoAcc list (from Sample_geo_accession)
    titles      = []   # title list  (from Sample_title)
    char_rows   = []   # list of [val, val, ...] per characteristic

    opener = gzip.open if mat_file.endswith(".gz") else open
    try:
        with opener(mat_file, "rt", errors="ignore") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith("!Sample_geo_accession"):
                    parts = line.split("\t")
                    sample_ids = [p.strip().strip('"')
                                  for p in parts[1:]]
                elif line.startswith("!Sample_title"):
                    parts  = line.split("\t")
                    titles = [p.strip().strip('"')
                              for p in parts[1:]]
                elif line.startswith("!Sample_characteristics_ch1"):
                    parts = line.split("\t")
                    char_rows.append(
                        [p.strip().strip('"')
                         for p in parts[1:]])
                elif "series_matrix_table_begin" in line.lower():
                    break
    except Exception as e:
        log(f"  Matrix parse error: {e}")
        return None, {}

    if not sample_ids:
        log("  No sample IDs found in matrix")
        return None, {}

    meta = pd.DataFrame(index=sample_ids)

    # Add title as a column if we got them
    if titles and len(titles) == len(sample_ids):
        meta["title"] = titles
        log(f"  Parsed {len(titles)} sample titles")
    else:
        log(f"  WARNING: title count ({len(titles)}) != "
            f"sample count ({len(sample_ids)})")

    # Parse characteristics
    for row in char_rows:
        if not row or ":" not in row[0]: continue
        key  = (row[0].split(":", 1)[0].strip()
                .lower().replace(" ", "_"))
        vals = [(c.split(": ", 1)[1].strip()
                 if ": " in c else c.strip())
                for c in row]
        if len(vals) == len(sample_ids):
            meta[key] = vals

    # Build bridge: title → geoAcc
    title_to_geoacc = {}
    if "title" in meta.columns:
        for geoAcc, title in zip(meta.index,
                                  meta["title"]):
            t = str(title).strip()
            if t:
                title_to_geoacc[t] = str(geoAcc).strip()
        log(f"  title→geoAcc bridge: {len(title_to_geoacc)} "
            f"entries")
        # Show a few examples
        examples = list(title_to_geoacc.items())[:3]
        for t, g in examples:
            log(f"    '{t}' → '{g}'")

    return meta, title_to_geoacc

def load_gse96058():
    log(""); log("=" * 65)
    log("LOAD: GSE96058"); log("=" * 65)

    # Expression
    expr = None
    for cand in [GSE96058_EXPR_FILE,
                 GSE96058_EXPR_FILE[:-3]]:
        if not os.path.exists(cand): continue
        try:
            kw = {"index_col": 0, "low_memory": False}
            if cand.endswith(".gz"):
                kw["compression"] = "gzip"
            tmp = pd.read_csv(cand, **kw)
            if tmp.shape[0] < tmp.shape[1]: tmp = tmp.T
            tmp.index = (tmp.index.astype(str)
                         .str.strip().str.upper())
            expr = tmp
            log(f"  GSE96058 expr: {expr.shape}")
            # Show a few column name examples
            log(f"  Expr col examples: "
                f"{list(expr.columns[:3])}")
            break
        except Exception as e:
            log(f"  GSE96058 expr load failed: {e}")

    # Clinical + title bridge
    clin             = None
    title_to_geoacc  = {}
    for mf in [GSE96058_MAT1_FILE, GSE96058_MAT2_FILE]:
        if not os.path.exists(mf): continue
        try:
            parsed, bridge = _parse_gse_matrix(mf)
            if parsed is not None and len(parsed) > 10:
                if clin is None:
                    clin            = parsed
                    title_to_geoacc = bridge
                else:
                    # Merge both platform matrices
                    clin = pd.concat(
                        [clin, parsed], axis=0,
                        ignore_index=False)
                    title_to_geoacc.update(bridge)
                log(f"  GSE96058 clinical after "
                    f"{os.path.basename(mf)}: "
                    f"{clin.shape}")
        except Exception as e:
            log(f"  GSE96058 matrix parse: {e}")

    if clin is not None:
        log(f"  GSE96058 clinical total: {clin.shape}")
        log(f"  Cols: {list(clin.columns)}")
        log(f"  title→geoAcc entries total: "
            f"{len(title_to_geoacc)}")

        # Audit PAM50
        for col in clin.columns:
            if "pam50" in col.lower():
                vals = clin[col].value_counts(dropna=False)
                log(f"  PAM50 '{col}': "
                    + str(dict(vals)))
                break

        # Confirm OS columns
        ot, oe = _find_os(clin)
        if ot and oe:
            log(f"  Survival: time='{ot}' event='{oe}'")
            if "day" in ot.lower():
                log("  Units: DAYS → will convert to months")
        else:
            log("  WARNING: OS columns not found")
            log(f"  Available: {list(clin.columns)}")

    log("")
    log("  NOTE: GSE96058 PAM50 = 5-class (no claudin-low).")
    log("  EXT-P4/P8/P9 → N/A.  EXT-P1b → will run.")

    return expr, clin, title_to_geoacc

# ============================================================
# SURVIVAL UTILITIES
# ============================================================

def _find_os(clin):
    TIME_C = [
        "OS_MONTHS", "os_months", "Overall_Survival_months",
        "overall_survival_months", "SURVIVAL_MONTHS",
        "ovrallSurvDays", "overall_survival_days",
        "overall_survival_days_", "overall_survival_(days)",
        "overall_survival", "T", "t", "time",
    ]
    EVT_C = [
        "OS_STATUS", "os_status", "vital_status",
        "Overall_Survival_Status", "VITAL_STATUS",
        "death_from_cancer",
        "ovrallSurvEvent", "overall_survival_event",
        "overall_survival_event_", "overall_survival_status",
        "E", "e", "event",
    ]
    ot = next((c for c in TIME_C if c in clin.columns), None)
    oe = next((c for c in EVT_C  if c in clin.columns), None)
    if not ot:
        for c in clin.columns:
            cl = c.lower().replace(" ", "_")
            if (("surv" in cl or "ovr" in cl or
                     "os_" in cl[:3]) and
                    ("day" in cl or "month" in cl or
                     "time" in cl)):
                ot = c; break
    if not oe:
        for c in clin.columns:
            cl = c.lower().replace(" ", "_")
            if (("event" in cl or "status" in cl or
                     "vital" in cl or "death" in cl) and
                    ("surv" in cl or "ovr" in cl or
                     "os" in cl or "cancer" in cl)):
                oe = c; break
    return ot, oe

def build_sv(expr, clin, ot, oe, label="",
             title_to_geoacc=None):
    """
    Build survival dicts.
    title_to_geoacc: if provided (GSE96058), maps
      expression column name (title) → geoAcc (clin index).
    Without it (METABRIC), direct string match is used.
    Auto-converts days→months when 'day' in ot name.
    """
    Td, Ed = {}, {}
    if clin is None or not ot or not oe:
        return Td, Ed

    days_mode = "day" in ot.lower()

    def pe(v):
        sv = str(v).upper().strip()
        if sv in ("1","1:DECEASED","DEAD","DECEASED",
                  "TRUE","YES","DIED","1.0"):  return 1.0
        if sv in ("0","0:LIVING","ALIVE","LIVING",
                  "FALSE","NO","CENSORED","0.0"): return 0.0
        try:
            f = float(sv); return 1.0 if f >= 0.5 else 0.0
        except Exception: return np.nan

    # Build index lookup: clin_key → clin_index_value
    cidx = {}
    for i in clin.index:
        cidx[str(i).strip()]       = i
        cidx[str(i).strip()[:12]]  = i

    n = 0
    for s in expr.columns:
        s_str = str(s).strip()

        # Determine the clinical lookup key
        if title_to_geoacc:
            # GSE96058: expression col is title → look up geoAcc
            clin_key = title_to_geoacc.get(s_str, "")
        else:
            # METABRIC: expression col IS the clinical ID
            clin_key = s_str

        idx = cidx.get(clin_key) or cidx.get(clin_key[:12])
        if idx is None: continue

        try:
            t_raw = float(clin.loc[idx, ot])
            e     = pe(clin.loc[idx, oe])
            if np.isnan(t_raw) or np.isnan(e) or t_raw <= 0:
                continue
            t      = t_raw / 30.44 if days_mode else t_raw
            Td[s]  = t
            Ed[s]  = e
            n     += 1
        except Exception:
            pass

    log(f"  {label} survival matched: {n} "
        f"({'days→months' if days_mode else 'months'})")
    return Td, Ed

def fmt_p(p):
    if p is None or np.isnan(p): return "p=nan"
    if p < 1e-5: return f"p={p:.2e} **"
    if p < 0.01: return f"p={p:.4f} **"
    if p < 0.05: return f"p={p:.3f} *"
    if p < 0.10: return f"p={p:.3f} (trend)"
    return              f"p={p:.3f}"

def get_sv(ss, Td, Ed):
    T, E, ok = [], [], []
    for s in ss:
        if s in Td and s in Ed:
            t, e = Td[s], Ed[s]
            if not np.isnan(t) and not np.isnan(e) and t > 0:
                T.append(t); E.append(int(e)); ok.append(s)
    return np.array(T), np.array(E), ok

def lrt(Th, Eh, Tl, El, lbl=""):
    if len(Th) < MIN_N or len(Tl) < MIN_N:
        log(f"    {lbl}: n too small"); return np.nan, np.nan
    if int(Eh.sum()) + int(El.sum()) < MIN_EVT * 2:
        log(f"    {lbl}: events too few"); return np.nan, np.nan
    try:
        p  = logrank_test(Th, Tl, Eh, El).p_value
        rh = Eh.sum() / Th.sum() if Th.sum() else np.nan
        rl = El.sum() / Tl.sum() if Tl.sum() else np.nan
        hr = rh / rl if (rl and rl > 0) else np.nan
        log(f"    {lbl}: n={len(Th)}/{len(Tl)} "
            f"ev={int(Eh.sum())}/{int(El.sum())} "
            f"{fmt_p(p)} HR≈{hr:.3f}")
        return p, hr
    except Exception as e:
        log(f"    {lbl}: {e}"); return np.nan, np.nan

def tertile(ss, sd, Td, Ed):
    sc = [(s, sd[s]) for s in ss
          if s in sd and not np.isnan(sd[s])]
    if len(sc) < MIN_N * 3:
        return (np.array([]),) * 6, [], [], []
    sc.sort(key=lambda x: x[1])
    n  = len(sc); t1 = n // 3; t2 = 2 * n // 3
    lo = [s for s, _ in sc[:t1]]
    mi = [s for s, _ in sc[t1:t2]]
    hi = [s for s, _ in sc[t2:]]
    Th, Eh, vh = get_sv(hi, Td, Ed)
    Tm, Em, vm = get_sv(mi, Td, Ed)
    Tl, El, vl = get_sv(lo, Td, Ed)
    return (Th, Eh, Tl, El, Tm, Em), vh, vm, vl

def medsplit(ss, sd, Td, Ed):
    sc = [(s, sd[s]) for s in ss
          if s in sd and not np.isnan(sd[s])]
    if len(sc) < MIN_N * 2:
        return (np.array([]),) * 4, [], []
    med = np.median([v for _, v in sc])
    hi  = [s for s, v in sc if v >= med]
    lo  = [s for s, v in sc if v < med]
    Th, Eh, vh = get_sv(hi, Td, Ed)
    Tl, El, vl = get_sv(lo, Td, Ed)
    return Th, Eh, Tl, El, vh, vl

def add(sc, pid, ds, nh, nl, p, hr, dok, conf, note=""):
    sc.append({"prediction": pid, "dataset": ds,
               "n_high": nh, "n_low": nl,
               "p_value": p, "HR_approx": hr,
               "direction_correct": dok,
               "confirmed": conf, "note": note})

# ============================================================
# CLASSIFICATION & SCORING
# ============================================================

def classify_cl(expr, label=""):
    pos   = [g for g in CL_POS if g in expr.index]
    neg   = [g for g in CL_NEG if g in expr.index]
    found = pos + neg
    if len(found) < 4:
        log(f"  {label}: only {len(found)} sig genes")
        return [], {}
    med     = {g: float(np.nanmedian(
                   expr.loc[g].values.astype(float)))
               for g in found}
    samples = list(expr.columns)
    scores  = {}
    for s in samples:
        sc = 0
        for g in pos:
            v = float(expr.loc[g, s])
            if not np.isnan(v) and v > med[g]: sc += 1
        for g in neg:
            v = float(expr.loc[g, s])
            if not np.isnan(v) and v < med[g]: sc += 1
        scores[s] = sc
    cl = [s for s, sc in scores.items()
          if sc >= CL_THRESHOLD]
    log(f"  {label}: CL n={len(cl)}/{len(samples)}")
    return cl, scores

def depth_sc(expr, samples):
    pos = [g for g in CL_POS if g in expr.index]
    neg = [g for g in CL_NEG if g in expr.index]
    if not pos or not neg: return {}
    sub = expr[samples]
    def z(arr):
        v = arr[~np.isnan(arr)]
        if len(v) < 2: return np.zeros_like(arr)
        mu, sd = v.mean(), v.std()
        return (arr-mu)/sd if sd > 1e-9 else np.zeros_like(arr)
    pz = np.column_stack(
        [z(sub.loc[g].values.astype(float)) for g in pos])
    nz = np.column_stack(
        [z(sub.loc[g].values.astype(float)) for g in neg])
    return dict(zip(samples,
                    np.mean(pz,axis=1)-np.mean(nz,axis=1)))

def all_scores(expr, cl):
    eps = 1e-3
    def sv(g, s):
        return float(expr.loc[g,s]) if g in expr.index else np.nan
    foxp3 = {s: sv("FOXP3",s) for s in cl}
    cd8a  = {s: sv("CD8A", s) for s in cl}
    gzmb  = {s: sv("GZMB", s) for s in cl}
    prf1  = {s: sv("PRF1", s) for s in cl}
    tigit = {s: sv("TIGIT",s) for s in cl}
    treg  = {}
    for s in cl:
        num = sum(v for v in [foxp3[s],tigit[s]]
                  if not np.isnan(v))
        den = sum(v for v in [cd8a[s],gzmb[s],prf1[s]]
                  if not np.isnan(v))
        treg[s] = num/(den+eps) if den > 0 else np.nan
    def composite(genes):
        parts = []
        for g in genes:
            if g not in expr.index: continue
            v  = np.array([sv(g,s) for s in cl])
            vv = v[~np.isnan(v)]
            if len(vv) >= MIN_N and vv.std() > 1e-9:
                parts.append((v-vv.mean())/vv.std())
        return dict(zip(cl,np.mean(parts,axis=0))) if parts else {}
    im  = composite(["FOXP3","PDCD1","TIGIT","LAG3"])
    ct  = composite(CT_GENES)
    mem = composite(MEMORY_GENES)
    log(f"  CT genes active: "
        f"{sum(1 for g in CT_GENES if g in expr.index)}")
    log(f"  Memory genes: "
        f"{[g for g in MEMORY_GENES if g in expr.index]}")
    return {"treg":treg,"immune":im,"ct":ct,"memory":mem}

# ============================================================
# ANALYSIS BLOCKS
# ============================================================

def block1(ds, cl, dp, Td, Ed, sc_list, esr1=None):
    log(""); log("="*65)
    log(f"BLOCK 1 — DEPTH SURVIVAL: {ds}"); log("="*65)
    pop = cl
    if esr1:
        med = np.nanmedian(
            [esr1.get(s,np.nan) for s in cl
             if not np.isnan(esr1.get(s,np.nan))])
        pop = [s for s in cl if esr1.get(s,np.nan) < med]
        log(f"  ESR1-low stratum: n={len(pop)}")
    (Th,Eh,Tl,El,Tm,Em),_,_,_ = tertile(pop,dp,Td,Ed)
    p,hr = lrt(Th,Eh,Tl,El,f"{ds} depth tertile")
    pid  = "EXT-P1" if ds == "METABRIC" else "EXT-P1b"
    ok   = not np.isnan(p) and p < 0.05 and hr > 2.0
    conf = "CONFIRMED" if ok else "NOT CONFIRMED"
    log(f"  [{pid}] {conf}")
    add(sc_list, pid, ds, len(Th), len(Tl), p, hr, hr>2.0, conf)
    return [(Th,Eh,f"Depth-high n={len(Th)}","#f85149"),
            (Tm,Em,f"Mid n={len(Tm)}","#6e7681"),
            (Tl,El,f"Depth-low n={len(Tl)}","#3fb950")]

def block2(cl, expr, Td, Ed, sc_list):
    log(""); log("="*65)
    log("BLOCK 2 — INDIVIDUAL GENE SURVIVAL"); log("="*65)
    km = {}
    for pid, gene, thr, direction in [
        ("EXT-P2",  "CLDN3", 0.10, "down"),
        ("EXT-P2b", "MKI67", 0.10, "up"),
    ]:
        log(f"\n  ── {pid}: {gene} ──")
        if gene in expr.index:
            d = {s: float(expr.loc[gene,s]) for s in cl}
            Th,Eh,Tl,El,_,_ = medsplit(cl,d,Td,Ed)
            p,hr = lrt(Th,Eh,Tl,El,gene)
            dok  = not np.isnan(hr) and (
                hr < 0.80 if direction=="down" else hr > 1.30)
            conf = ("CONFIRMED"
                    if not np.isnan(p) and p<thr and dok
                    else "NOT CONFIRMED")
            log(f"  [{pid}] {conf}")
            add(sc_list,pid,"METABRIC",len(Th),len(Tl),
                p,hr,dok,conf)
            if gene == "CLDN3":
                km["CLDN3"] = [
                    (Th,Eh,f"CLDN3-high n={len(Th)}","#58a6ff"),
                    (Tl,El,f"CLDN3-low n={len(Tl)}","#f85149")]
        else:
            add(sc_list,pid,"METABRIC",0,0,
                np.nan,np.nan,False,"DATA MISSING")
    # EXT-P3
    log("\n  ── EXT-P3: immune null ──")
    ig = [g for g in ["FOXP3","PDCD1","TIGIT","LAG3"]
          if g in expr.index]
    if ig:
        parts = []
        for g in ig:
            v  = np.array([float(expr.loc[g,s]) for s in cl])
            vv = v[~np.isnan(v)]
            if len(vv)>=MIN_N and vv.std()>1e-9:
                parts.append((v-vv.mean())/vv.std())
        if parts:
            ic = dict(zip(cl,np.mean(parts,axis=0)))
            Th,Eh,Tl,El,_,_ = medsplit(cl,ic,Td,Ed)
            p,hr = lrt(Th,Eh,Tl,El,"immune composite")
            null = np.isnan(p) or p > 0.15
            conf = "CONFIRMED" if null else "NOT CONFIRMED"
            log(f"  [EXT-P3] {conf} (null=p>0.15)")
            add(sc_list,"EXT-P3","METABRIC",len(Th),len(Tl),
                p,hr,null,conf,"negative prediction")
        else:
            add(sc_list,"EXT-P3","METABRIC",0,0,
                np.nan,np.nan,False,"DATA MISSING")
    else:
        add(sc_list,"EXT-P3","METABRIC",0,0,
            np.nan,np.nan,False,"DATA MISSING")
    return km, sc_list

def block3(cl, scrs, Td, Ed, sc_list):
    log(""); log("="*65)
    log("BLOCK 3 — SCRIPT 3 REPLICATION"); log("="*65)
    km = {}
    # EXT-P5
    log("\n  ── EXT-P5 Treg:eff ──")
    tr = scrs.get("treg",{})
    if tr:
        (Th,Eh,Tl,El,Tm,Em),_,_,_ = tertile(cl,tr,Td,Ed)
        p,hr = lrt(Th,Eh,Tl,El,"EXT-P5 Treg:eff")
        dok  = not np.isnan(hr) and hr > 1.5
        conf = ("CONFIRMED"
                if not np.isnan(p) and p<0.05 and dok
                else "NOT CONFIRMED")
        log(f"  [EXT-P5] {conf}")
        add(sc_list,"EXT-P5","METABRIC",len(Th),len(Tl),
            p,hr,dok,conf)
        km["EXT-P5"]=[(Th,Eh,f"Treg-hi n={len(Th)}","#f85149"),
                      (Tl,El,f"Treg-lo n={len(Tl)}","#3fb950")]
    else:
        add(sc_list,"EXT-P5","METABRIC",0,0,
            np.nan,np.nan,False,"DATA MISSING")
    # EXT-P6
    log("\n  ── EXT-P6 CT×memory ──")
    mem = scrs.get("memory",{})
    ct  = scrs.get("ct",{})
    if mem and ct:
        mv = [mem[s] for s in cl if s in mem
              and not np.isnan(mem[s])]
        if len(mv) >= MIN_N*2:
            mm  = np.median(mv)
            mh  = [s for s in cl if s in mem and mem[s]>=mm]
            ml  = [s for s in cl if s in mem and mem[s]< mm]
            cth = [ct[s] for s in mh if s in ct
                   and not np.isnan(ct[s])]
            ctl = [ct[s] for s in ml if s in ct
                   and not np.isnan(ct[s])]
            if len(cth)>=MIN_N and len(ctl)>=MIN_N:
                _,p = stats.mannwhitneyu(
                    ctl,cth,alternative="greater")
                mml = np.mean(ctl); mmh = np.mean(cth)
                log(f"    mem-low={mml:.4f} "
                    f"mem-high={mmh:.4f} {fmt_p(p)}")
                dok  = mml > mmh
                conf = ("CONFIRMED"
                        if p<0.001 and dok else "NOT CONFIRMED")
                log(f"  [EXT-P6] {conf}")
                add(sc_list,"EXT-P6","METABRIC",
                    len(ctl),len(cth),p,
                    mml/mmh if mmh else np.nan,dok,conf)
                km["EXT-P6"] = (ctl,cth)
            else:
                add(sc_list,"EXT-P6","METABRIC",0,0,
                    np.nan,np.nan,False,"DATA MISSING")
        else:
            add(sc_list,"EXT-P6","METABRIC",0,0,
                np.nan,np.nan,False,"DATA MISSING")
    else:
        add(sc_list,"EXT-P6","METABRIC",0,0,
            np.nan,np.nan,False,"DATA MISSING")
    # EXT-P7
    log("\n  ── EXT-P7 memory OS null ──")
    if mem:
        mv = [mem[s] for s in cl if s in mem
              and not np.isnan(mem[s])]
        if len(mv) >= MIN_N*2:
            mm  = np.median(mv)
            mh  = [s for s in cl if s in mem and mem[s]>=mm]
            ml  = [s for s in cl if s in mem and mem[s]< mm]
            Tmh,Emh,_ = get_sv(mh,Td,Ed)
            Tml,Eml,_ = get_sv(ml,Td,Ed)
            p,hr = lrt(Tml,Eml,Tmh,Emh,"EXT-P7 mem OS")
            null = np.isnan(p) or p > 0.30
            conf = "CONFIRMED" if null else "NOT CONFIRMED"
            log(f"  [EXT-P7] {conf} (null=p>0.30)")
            add(sc_list,"EXT-P7","METABRIC",
                len(Tml),len(Tmh),p,hr,null,conf,
                "negative prediction")
            km["EXT-P7"]=[(Tmh,Emh,
                f"Mem-high n={len(Tmh)}","#3fb950"),
               (Tml,Eml,f"Mem-low n={len(Tml)}","#f85149")]
        else:
            add(sc_list,"EXT-P7","METABRIC",0,0,
                np.nan,np.nan,False,"DATA MISSING")
    else:
        add(sc_list,"EXT-P7","METABRIC",0,0,
            np.nan,np.nan,False,"DATA MISSING")
    return km, sc_list

def block4(expr_g, clin_g, title_to_geoacc, sc_list):
    log(""); log("="*65)
    log("BLOCK 4 — GSE96058"); log("="*65)
    km = {}
    if expr_g is None:
        log("  GSE96058 expression not available — skipped")
        for pid in ["EXT-P1b","EXT-P4","EXT-P8","EXT-P9"]:
            add(sc_list,pid,"GSE96058",0,0,
                np.nan,np.nan,False,"DATA MISSING")
        return km, sc_list

    geo_cl,_ = classify_cl(expr_g,"GSE96058-geo")
    all_g    = list(expr_g.columns)
    dp_g     = depth_sc(expr_g,all_g)

    Tg, Eg = {}, {}
    if clin_g is not None:
        ot, oe = _find_os(clin_g)
        if ot and oe:
            Tg, Eg = build_sv(
                expr_g, clin_g, ot, oe, "GSE96058",
                title_to_geoacc=title_to_geoacc)
        else:
            log("  GSE96058: no OS columns found")

    # EXT-P8
    log("\n  ── EXT-P8 overlap ──")
    log("  [EXT-P8] N/A — GSE96058 PAM50 has no "
        "claudin-low category")
    add(sc_list,"EXT-P8","GSE96058",len(geo_cl),0,
        np.nan,np.nan,False,
        "canonical CL N/A in SCAN-B PAM50")

    if Tg:
        esr1_g = ({s: float(expr_g.loc["ESR1",s])
                   for s in geo_cl}
                  if "ESR1" in expr_g.index else None)
        pop_b = geo_cl
        if esr1_g:
            m     = np.nanmedian(list(esr1_g.values()))
            pop_b = [s for s in geo_cl
                     if esr1_g.get(s,np.nan) < m]
            log(f"  ESR1-low geo CL: n={len(pop_b)}")

        # EXT-P1b
        log("\n  ── EXT-P1b: depth vs OS geo CL ──")
        (Th,Eh,Tl,El,Tm,Em),_,_,_ = tertile(pop_b,dp_g,Tg,Eg)
        p1b,hr1b = lrt(Th,Eh,Tl,El,"EXT-P1b geo")
        d1b  = not np.isnan(hr1b) and hr1b > 2.0
        c1b  = ("CONFIRMED"
                if not np.isnan(p1b) and p1b<0.05 and d1b
                else "NOT CONFIRMED")
        log(f"  [EXT-P1b] {c1b}")
        add(sc_list,"EXT-P1b","GSE96058",
            len(Th),len(Tl),p1b,hr1b,d1b,c1b)
        km["EXT-P1b"]=[(Th,Eh,f"Depth-high n={len(Th)}","#f85149"),
                       (Tl,El,f"Depth-low n={len(Tl)}","#3fb950")]

        # EXT-P9/P4 require canonical CL — unavailable
        for pid in ["EXT-P9","EXT-P4"]:
            add(sc_list,pid,"GSE96058",0,0,
                np.nan,np.nan,False,
                "canonical CL N/A in SCAN-B PAM50")
    else:
        log("  No GSE96058 survival matched")
        for pid in ["EXT-P1b","EXT-P9","EXT-P4"]:
            add(sc_list,pid,"GSE96058",0,0,
                np.nan,np.nan,False,"no survival matched")

    return km, sc_list

# ============================================================
# FIGURE
# ============================================================

def make_figure(km_m, km_g, scorecard):
    log(""); log("="*65); log("GENERATING FIGURE"); log("="*65)
    fig = plt.figure(figsize=(22,18))
    fig.patch.set_facecolor("#0d1117")
    gs  = gridspec.GridSpec(3,4,figure=fig,
                            hspace=0.48,wspace=0.40)
    RED, GREEN = "#f85149","#3fb950"
    def style(ax):
        ax.set_facecolor("#161b22")
        for sp in ax.spines.values():
            sp.set_edgecolor("#30363d")
        ax.tick_params(colors="#8b949e",labelsize=7)
        ax.xaxis.label.set_color("#8b949e")
        ax.yaxis.label.set_color("#8b949e")
    lkw = dict(color="#8b949e",fontsize=7)
    tkw = dict(color="#e6edf3",fontsize=8.5,
               fontweight="bold",pad=5)
    def km_p(ax,series,ci=True):
        for T,E,lbl,c in series:
            if len(T)>=MIN_N:
                KaplanMeierFitter().fit(T,E,label=lbl)\
                    .plot_survival_function(ax=ax,ci_show=ci,
                        color=c,linewidth=1.8)
    for col,key,title in [
        (0,"METABRIC_depthA",
         "Depth vs OS — METABRIC\nStratum A"),
        (1,"METABRIC_depthB",
         "Depth vs OS — METABRIC\nESR1-low (EXT-P1)"),
        (2,"EXT-P5",
         "Treg:Eff vs OS\nMETABRIC (EXT-P5)"),
        (3,"EXT-P7",
         "Memory OS — null\nMETABRIC (EXT-P7)"),
    ]:
        ax = fig.add_subplot(gs[0,col])
        ax.set_title(title,**tkw)
        d = km_m.get(key)
        if d: km_p(ax,d,ci=(col==1))
        ax.set_xlabel("Time (months)",**lkw)
        ax.set_ylabel("Survival",**lkw)
        ax.legend(fontsize=5,labelcolor="#8b949e",
                  facecolor="#161b22",edgecolor="#30363d")
        style(ax)
    for col,key,title,span in [
        (0,"EXT-P1b","Depth vs OS GSE96058\nGeo CL (EXT-P1b)",2),
        (2,"EXT-P9","Depth vs OS GSE96058\nCanonical (EXT-P9)",2),
    ]:
        ax = fig.add_subplot(gs[1,col:col+span])
        ax.set_title(title,**tkw)
        d = km_g.get(key)
        if d: km_p(ax,d)
        ax.set_xlabel("Time (months)",**lkw)
        ax.set_ylabel("Survival",**lkw)
        ax.legend(fontsize=6,labelcolor="#8b949e",
                  facecolor="#161b22",edgecolor="#30363d")
        style(ax)
    ax_ct = fig.add_subplot(gs[2,0])
    ax_ct.set_title("CT Antigen × Memory\nMETABRIC (EXT-P6)",**tkw)
    ctd = km_m.get("EXT-P6")
    if ctd:
        ctl,cth = ctd
        for i,(vals,lbl,c) in enumerate([
            (ctl,"Mem-low",RED),(cth,"Mem-high",GREEN)]):
            if vals:
                ax_ct.boxplot(vals,positions=[i],widths=0.6,
                    patch_artist=True,
                    boxprops=dict(facecolor=c,alpha=0.5),
                    medianprops=dict(color="#e6edf3",lw=2),
                    whiskerprops=dict(color="#8b949e"),
                    capprops=dict(color="#8b949e"),
                    flierprops=dict(marker="o",alpha=0.2,
                        markerfacecolor=c,markeredgewidth=0))
        ax_ct.set_xticks([0,1])
        ax_ct.set_xticklabels(["Mem-low","Mem-high"],
                              color="#8b949e",fontsize=7)
        ax_ct.set_ylabel("CT Z-score",**lkw)
    style(ax_ct)
    ax_sc = fig.add_subplot(gs[2,1:])
    ax_sc.set_title("Prediction Scorecard",**tkw)
    if scorecard:
        df   = pd.DataFrame(scorecard)
        pids = df["prediction"].tolist()
        bclr = [GREEN if ("CONFIRMED" in str(c) and
                          "NOT" not in str(c))
                else ("#d29922" if ("MISSING" in str(c) or
                                    "N/A" in str(c)) else RED)
                for c in df["confirmed"].tolist()]
        hrs  = [r["HR_approx"]
                if not np.isnan(r["HR_approx"]) else 0
                for r in scorecard]
        ax_sc.barh(pids,hrs,color=bclr,alpha=0.7,
                   edgecolor="#30363d")
        ax_sc.axvline(1.0,color="#8b949e",lw=0.8,ls="--")
        ax_sc.axvline(2.0,color=GREEN,lw=0.8,ls=":")
        ax_sc.set_xlabel("HR / overlap metric",**lkw)
    style(ax_sc)
    fig.suptitle(
        "CLAUDIN-LOW — SCRIPT 4 — EXTERNAL VALIDATION\n"
        "METABRIC + GSE96058  |  OrganismCore  |  "
        "BRCA-S7k  |  2026-03-05  |  TYPE 4 ROOT LOCK",
        color="#e6edf3",fontsize=10,fontweight="bold",y=0.998)
    plt.savefig(FIG_FILE,dpi=150,bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close()
    log(f"  Figure: {FIG_FILE}")

# ============================================================
# SCORECARD
# ============================================================

def print_scorecard(sc_list):
    log(""); log("="*65)
    log("PREDICTION SCORECARD"); log("="*65)
    log(f"  {'ID':<10} {'Dataset':<12} {'n_h':>5} {'n_l':>5} "
        f"{'p':>9} {'HR':>7}  Status")
    log("  "+"─"*68)
    for r in sc_list:
        ps = (f"{r['p_value']:.3f}"
              if not np.isnan(r['p_value']) else "nan")
        hs = (f"{r['HR_approx']:.3f}"
              if not np.isnan(r['HR_approx']) else "nan")
        log(f"  {r['prediction']:<10} "
            f"{str(r['dataset']):<12} "
            f"{r['n_high']:>5} {r['n_low']:>5} "
            f"{ps:>9} {hs:>7}  {r['confirmed']}")
    nc = sum(1 for r in sc_list
             if "CONFIRMED" in str(r["confirmed"])
             and "NOT" not in str(r["confirmed"]))
    log(f"\n  Confirmed: {nc} / {len(sc_list)}")
    pd.DataFrame(sc_list).to_csv(SCORECARD_FILE,index=False)
    log(f"  Scorecard: {SCORECARD_FILE}")

# ============================================================
# MAIN
# ============================================================

def main():
    log("="*65)
    log("BRCA CLAUDIN-LOW — SCRIPT 4 (v6)")
    log("OrganismCore — BRCA-S7k | 2026-03-05 | TYPE 4")
    log(f"Output: {BASE_DIR}")
    log("="*65)

    sc_list = []
    km_m    = {}
    km_g    = {}
    expr_m  = None

    # Delete cached bad POST file (v3 parse bug: empty symbols)
    if os.path.exists(METABRIC_POST_FILE):
        try:
            with open(METABRIC_POST_FILE) as f:
                test = json.load(f)
            if isinstance(test,list) and test:
                first = test[0]
                sym   = (first.get("gene",{})
                         .get("hugoGeneSymbol","") or
                         first.get("hugoGeneSymbol",""))
                if not sym:
                    log("  Deleting cached POST (empty symbols)")
                    os.remove(METABRIC_POST_FILE)
        except Exception:
            pass

    # ── METABRIC ─────────────────────────────────────────────
    log("\n"+"="*65+"\nMETABRIC ANALYSIS\n"+"="*65)
    profile_id = discover_metabric_profile()
    sample_ids = get_metabric_samples()

    if acquire_metabric_expr_post(profile_id, sample_ids):
        expr_m = load_metabric_expr()

    if expr_m is None:
        log("\n  POST route failed — trying Zenodo combo CSV")
        if acquire_metabric_combo():
            expr_m = load_metabric_expr()

    if expr_m is not None:
        acquire_metabric_clinical()
        clin_m, ot_m, oe_m = load_metabric_clinical()
        cl_m,_   = classify_cl(expr_m,"METABRIC")
        dp_m     = depth_sc(expr_m,list(expr_m.columns))
        scr_m    = all_scores(expr_m,cl_m)
        Tm, Em   = (build_sv(expr_m,clin_m,ot_m,oe_m,
                              "METABRIC")
                    if clin_m is not None else ({},{}))
        esr1_m   = ({s: float(expr_m.loc["ESR1",s])
                     for s in cl_m}
                    if "ESR1" in expr_m.index else None)
        km_b1    = block1("METABRIC",cl_m,dp_m,
                           Tm,Em,sc_list,esr1_m)
        km_m["METABRIC_depthA"] = km_b1
        km_m["METABRIC_depthB"] = km_b1
        km2,sc_list = block2(cl_m,expr_m,Tm,Em,sc_list)
        km3,sc_list = block3(cl_m,scr_m,Tm,Em,sc_list)
        km_m.update(km2); km_m.update(km3)
    else:
        log("  METABRIC expression unavailable — skipped")
        for pid in ["EXT-P1","EXT-P2","EXT-P2b",
                    "EXT-P3","EXT-P5","EXT-P6","EXT-P7"]:
            add(sc_list,pid,"METABRIC",0,0,
                np.nan,np.nan,False,"DATA MISSING")

    # ── GSE96058 ──────────────────────────────────────────────
    log("\n"+"="*65+"\nGSE96058 ANALYSIS\n"+"="*65)
    acquire_gse96058()
    expr_g, clin_g, title_to_geoacc = load_gse96058()
    km_g_out, sc_list = block4(
        expr_g, clin_g, title_to_geoacc, sc_list)
    km_g.update(km_g_out)

    # ── Output ────────────────────────────────────────────────
    try:
        make_figure(km_m, km_g, sc_list)
    except Exception as e:
        log(f"  Figure failed: {e}")
        import traceback; log(traceback.format_exc())

    print_scorecard(sc_list)
    log("")
    log("="*65); log("SCRIPT 4 COMPLETE"); log("="*65)
    log(f"  Log:       {LOG_FILE}")
    log(f"  Figure:    {FIG_FILE}")
    log(f"  Scorecard: {SCORECARD_FILE}")
    log("")
    log("  KEY QUESTIONS FOR BRCA-S7l:")
    log("  1. EXT-P1:  depth HR>2.0, p<0.05 METABRIC?")
    log("  2. EXT-P1b: depth HR>2.0, p<0.05 GSE96058?")
    log("  3. EXT-P5:  Treg:eff OS confirmed?")
    log("  4. EXT-P7:  memory null replicated?")
    log("  5. EXT-P6:  CT antigen memory split replicated?")
    log("  NOTE: EXT-P4/P8/P9 N/A — no canonical CL in "
        "GSE96058 PAM50.")
    flush_log()


if __name__ == "__main__":
    main()
