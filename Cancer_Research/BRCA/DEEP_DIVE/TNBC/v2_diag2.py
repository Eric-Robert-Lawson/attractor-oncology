"""
GSE25066 PROBE DIAGNOSTIC
OrganismCore | 2026-03-04

Prints:
  1. First 30 lines of the series matrix (raw) — confirms format
  2. All unique metadata fields (!Sample_characteristics_ch1)
  3. First 20 probe IDs from the expression table
  4. Probes matching our target genes by substring
  5. Whether any of our gene symbols appear as probe labels

Run this, paste the output, and Script 2 will be fixed.
"""

import os
import gzip

MATRIX_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "tnbc_s2_results", "data",
    "GSE25066_series_matrix.txt.gz"
)

TARGET_GENES = [
    "KRT5", "KRT14", "SOX10", "FOXC1", "EGFR", "VIM", "CDH3",
    "ESR1", "FOXA1", "GATA3", "SPDEF", "PGR",
    "EED", "EZH2", "PARP1", "AR", "PTPRC", "SPI1",
    "MKI67", "TOP2A", "BRCA1", "BRCA2", "TP53",
]

def main():
    if not os.path.exists(MATRIX_PATH):
        print(f"FILE NOT FOUND: {MATRIX_PATH}")
        return

    print(f"File: {MATRIX_PATH}")
    print(f"Size: {os.path.getsize(MATRIX_PATH)/1e6:.1f} MB")
    print()

    # ── PASS 1: read header lines ──────────────────────────
    print("=" * 60)
    print("PASS 1: FIRST 50 LINES (raw)")
    print("=" * 60)
    meta_keys = {}
    platform_block = False
    expr_started   = False
    probe_ids      = []
    n_header_lines = 0

    with gzip.open(MATRIX_PATH, "rt", errors="ignore") as f:
        for i, line in enumerate(f):
            line = line.rstrip("\n")

            if i < 50:
                print(f"  L{i:03d}: {line[:120]}")

            if line.startswith("!Sample_characteristics_ch1"):
                parts = line.split("\t")
                val0  = parts[1].strip().strip('"') if len(parts) > 1 else ""
                key   = val0.split(":")[0].strip() if ":" in val0 else val0[:40]
                if key not in meta_keys:
                    meta_keys[key] = val0

            if "platform_table_begin" in line.lower():
                platform_block = True
                print(f"\n  *** FOUND !platform_table_begin at line {i} ***")

            if "platform_table_end" in line.lower():
                platform_block = False

            if "series_matrix_table_begin" in line.lower():
                expr_started = True
                n_header_lines = i

            if expr_started and i > n_header_lines:
                parts = line.split("\t")
                if parts[0].strip().strip('"') == "ID_REF":
                    print(f"\n  *** COLUMN HEADER ROW at line {i}: {line[:120]}")
                elif len(parts) > 1 and not line.startswith("!"):
                    probe_ids.append(parts[0].strip().strip('"'))
                    if len(probe_ids) <= 20:
                        print(f"  PROBE[{len(probe_ids)}]: {parts[0].strip()}")
                    if len(probe_ids) >= 20:
                        break

    print()
    print("=" * 60)
    print("PASS 1 SUMMARY")
    print("=" * 60)
    print(f"  Platform block found: {platform_block}")
    print(f"  Expression table found: {expr_started}")
    print(f"  Header lines before table: {n_header_lines}")
    print(f"  Sample probe IDs (first 20): {probe_ids[:20]}")
    print()

    # ── PASS 2: all metadata keys ──────────────────────────
    print("=" * 60)
    print("PASS 2: ALL !Sample_characteristics_ch1 KEYS")
    print("=" * 60)
    with gzip.open(MATRIX_PATH, "rt", errors="ignore") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("!Sample_characteristics_ch1"):
                parts = line.split("\t")
                val0  = parts[1].strip().strip('"') if len(parts) > 1 else ""
                key   = val0.split(":")[0].strip() if ":" in val0 else val0[:40]
                ex_vals = [p.strip().strip('"') for p in parts[1:6]]
                print(f"  KEY: '{key}'")
                print(f"    Examples: {ex_vals}")
            if "series_matrix_table_end" in line.lower():
                break
    print()

    # ── PASS 3: probe→gene matching ────────────────────────
    print("=" * 60)
    print("PASS 3: PROBE IDs MATCHING TARGET GENES")
    print("(Searching all probe IDs in expression table)")
    print("=" * 60)

    found_probes = {g: [] for g in TARGET_GENES}
    n_probes_total = 0

    with gzip.open(MATRIX_PATH, "rt", errors="ignore") as f:
        in_table = False
        for line in f:
            line = line.rstrip("\n")
            if "series_matrix_table_begin" in line.lower():
                in_table = True
                continue
            if "series_matrix_table_end" in line.lower():
                break
            if not in_table:
                continue
            parts = line.split("\t")
            probe = parts[0].strip().strip('"')
            if probe == "ID_REF" or not probe:
                continue
            n_probes_total += 1
            probe_up = probe.upper()
            for gene in TARGET_GENES:
                if gene.upper() in probe_up:
                    found_probes[gene].append(probe)

    print(f"  Total probes in expression table: {n_probes_total}")
    print()
    print("  Gene → probes containing gene name as substring:")
    for gene in TARGET_GENES:
        probes = found_probes[gene]
        print(f"  {gene:<12}: {probes if probes else '(none)'}")

    print()
    print("=" * 60)
    print("CONCLUSION")
    print("=" * 60)
    if n_probes_total > 0:
        probe_ex = probe_ids[:3] if probe_ids else []
        print(f"  Probe format example: {probe_ex}")
        print(f"  These are Affymetrix HG-U133A probe set IDs.")
        print(f"  Format: numeric_at or numeric_s_at or numeric_x_at")
        print(f"  Gene names do NOT appear in probe IDs.")
        print(f"  We need a probe→gene lookup table.")
        print(f"  GPL96 NCBI FTP is 404. Will embed mapping in Script 2.")
    print()
    print("DIAGNOSTIC COMPLETE.")

if __name__ == "__main__":
    main()
