#!/usr/bin/env python3
"""
NECROMASS DIAGNOSTIC SCRIPT
============================
Document ID:  NECROMASS_DIAGNOSTIC_v1.0
Date:         2026-03-13
Author:       Eric Robert Lawson / OrganismCore
Pre-reg DOI:  10.5281/zenodo.18986790

PURPOSE
-------
This is the DIAGNOSTIC script, not the analysis script.
Its job is to:

  1. State exactly what data we are testing against.
  2. State exactly what the hypothesis predicts.
  3. Define the data structure we need to build.
  4. Confirm that the published values are correctly
     understood before we run any statistics.
  5. Identify every gap, ambiguity, and assumption.
  6. Print a complete diagnostic report so we know
     exactly what the analysis script needs to do.

This script contains NO external dependencies beyond
the Python standard library.

It encodes the published data directly from:
  Changela & Bridges (2010)
  Meteoritics & Planetary Science 45:1847-1867
  DOI: 10.1111/j.1945-5100.2010.01123.x

And depth estimates from:
  Mikouchi et al. (2006)
  Wang et al. (2021) Earth Planets Space
  DOI: 10.1186/s40623-021-01492-3

Run this script first.
Read every line of output before proceeding.
The analysis script will only be written after
this diagnostic confirms the data structure is
correctly understood.

WHAT WE ARE TESTING
-------------------
The Mars Iron Necromass Hypothesis predicts:

  Iron oxide content in Martian subsurface rock
  should INCREASE with DEPTH in the original
  Mars rock pile, because:

    A subsurface biological iron-cycling community
    operates deeper in the Martian crust.
    Organisms at depth mine iron from olivine
    and basalt and precipitate it as iron oxide.
    The deeper the rock, the more time it spent
    in the biological iron-cycling zone.
    Therefore: Fe oxide enrichment should be
    highest in the deepest material.

The nakhlite meteorite depth sequence gives us
a direct test:
  We have physical Mars rock from three known
  depths in the same original rock pile.
  We can measure Fe oxide content at each depth.
  We can ask: does it increase with depth?

The null hypothesis:
  Fe oxide content is uniform with depth,
  or decreases with depth.
  Consistent with abiotic surface weathering
  that operates top-down, not bottom-up.

The biological hypothesis:
  Fe oxide content increases with depth.
  Consistent with a deep biological iron-cycling
  community operating from the bottom up.
"""

import sys
import math

# ─────────────────────────────────────────────────────────────
# SECTION 1: PUBLISHED DATA — NAKHLITE DEPTH SEQUENCE
# ─────────────────────────────────────────────────────────────
# Source: Changela & Bridges (2010)
# DOI: 10.1111/j.1945-5100.2010.01123.x
#
# These are the alteration vein FeO weight percent values
# from secondary iron oxide phases in nakhlite fractures.
# This is NOT the bulk rock FeO.
# This IS the iron oxide deposited by aqueous processes
# in the fracture network — the thing most likely to
# record biological iron cycling if it occurred.
#
# The distinction matters:
#   BULK FeO = iron in primary igneous minerals
#              (pyroxene, olivine)
#              Reflects original magma composition.
#              NOT what we are testing.
#
#   VEIN FeO  = iron in secondary alteration products
#              (ferrihydrite, goethite, Fe-oxide gel,
#               iddingsite, siderite)
#              Deposited from aqueous solution.
#              This IS what we are testing.
#              This is the record of water-rock
#              iron chemistry on Mars.
#              This is where biological iron cycling
#              would leave its signature.
# ─────────────────────────────────────────────────────────────

# Nakhlite meteorite data
# Each entry: (name, depth_m, vein_FeO_wt_pct_mean,
#              vein_FeO_wt_pct_min, vein_FeO_wt_pct_max,
#              alteration_volume_pct, notes)
#
# Depth values from Mikouchi et al. 2006 and Changela &
# Bridges 2010. Using midpoint of estimated range.
#
# FeO values from Changela & Bridges 2010, Table 4
# (secondary alteration vein compositions).
# These are representative values — the analysis script
# will need to load the full published table.

NAKHLITE_DATA = [
    {
        "name": "Yamato 000593",
        "name_short": "Y-000593",
        "depth_m": 3.5,           # estimated ~<7m, using 3.5m midpoint
        "depth_m_min": 0.0,
        "depth_m_max": 7.0,
        "depth_note": "Shallowest. <7m estimated. "
                      "Mikouchi et al. 2006.",
        "vein_FeO_mean": 64.7,    # wt% from Changela & Bridges Table 4
        "vein_FeO_min": 60.0,
        "vein_FeO_max": 69.0,
        "alteration_vol_pct": 0.1,  # ~0.1% by volume (least altered)
        "alteration_note": "Least altered nakhlite. "
                           "Thin veinlets only.",
        "source": "Changela & Bridges 2010, Table 4; "
                  "Mikouchi et al. 2006"
    },
    {
        "name": "Nakhla",
        "name_short": "Nakhla",
        "depth_m": 8.5,           # estimated 7-10m, using 8.5m midpoint
        "depth_m_min": 7.0,
        "depth_m_max": 10.0,
        "depth_note": "Intermediate depth. 7-10m estimated. "
                      "Changela & Bridges 2010.",
        "vein_FeO_mean": 51.7,    # wt% from Changela & Bridges Table 4
        "vein_FeO_min": 45.0,
        "vein_FeO_max": 55.0,
        "alteration_vol_pct": 0.5,  # intermediate alteration
        "alteration_note": "Intermediate alteration. "
                           "Iddingsite veins.",
        "source": "Changela & Bridges 2010, Table 4; "
                  "Changela & Bridges 2010 depth estimate"
    },
    {
        "name": "Governador Valadares",
        "name_short": "Gov.Val.",
        "depth_m": 8.5,           # same as Nakhla, paired
        "depth_m_min": 7.0,
        "depth_m_max": 10.0,
        "depth_note": "Intermediate depth. Paired with Nakhla. "
                      "7-10m estimated.",
        "vein_FeO_mean": 50.0,    # approximate, similar to Nakhla
        "vein_FeO_min": 44.0,
        "vein_FeO_max": 55.0,
        "alteration_vol_pct": 0.4,
        "alteration_note": "Intermediate alteration. "
                           "Similar to Nakhla.",
        "source": "Changela & Bridges 2010; "
                  "Mikouchi et al. 2006"
    },
    {
        "name": "Lafayette",
        "name_short": "Lafayette",
        "depth_m": 30.0,          # estimated 20-40m, using 30m midpoint
        "depth_m_min": 20.0,
        "depth_m_max": 40.0,
        "depth_note": "Deepest. 20-40m estimated. "
                      "Wang et al. 2021.",
        "vein_FeO_mean": 57.5,    # wt% from Changela & Bridges Table 4
        "vein_FeO_min": 54.9,
        "vein_FeO_max": 60.08,
        "alteration_vol_pct": 1.5,  # most altered nakhlite
        "alteration_note": "Most altered nakhlite. "
                           "Fe-oxide gel fills fractures.",
        "source": "Changela & Bridges 2010, Table 4; "
                  "Wang et al. 2021"
    },
    {
        "name": "NWA 998",
        "name_short": "NWA998",
        "depth_m": 30.0,          # estimated ~same as Lafayette
        "depth_m_min": 20.0,
        "depth_m_max": 40.0,
        "depth_note": "Deepest. Paired with Lafayette. "
                      "20-40m estimated.",
        "vein_FeO_mean": 52.0,    # approximate from published data
        "vein_FeO_min": 47.0,
        "vein_FeO_max": 55.0,
        "alteration_vol_pct": 1.2,
        "alteration_note": "Highly altered. "
                           "Deep pile position confirmed.",
        "source": "Changela & Bridges 2010; "
                  "Mikouchi et al. 2006"
    },
]

# ─────────────────────────────────────────────────────────────
# SECTION 2: ADDITIONAL TEST DATASETS
# ─────────────────────────────────────────────────────────────
# Beyond FeO in veins, we have two additional datasets
# that test the same hypothesis from different angles.
# ─────────────────────────────────────────────────────────────

# ALH84001 magnetite — biological morphology criteria
# Source: Thomas-Keprta et al. 2001 (GCA)
#         Thomas-Keprta et al. 2009 (GCA)
# This is the biological crystal morphology test.
# Six criteria for biogenic magnetite.
# Status of each criterion in ALH84001:
ALH84001_MAGNETITE = {
    "meteorite": "ALH84001",
    "age_Ga": 4.1,
    "depth_m": "Unknown — ancient crust sample",
    "magnetite_population": "Truncated hexa-octahedral (THO) subset",
    "criteria": [
        {
            "criterion": "1_morphology",
            "description": "Truncated hexa-octahedral (THO) shape",
            "status_in_ALH84001": "PRESENT in specific population",
            "abiotic_explained": False,
            "note": "THO morphology not reproducible abiotically "
                    "in published experiments. "
                    "Thomas-Keprta et al. 2001."
        },
        {
            "criterion": "2_size_range",
            "description": "Narrow size range 35-120nm "
                           "single-domain magnetic",
            "status_in_ALH84001": "PRESENT",
            "abiotic_explained": True,
            "note": "Size range can be produced abiotically "
                    "by thermal decomposition. "
                    "Golden et al. 2004."
        },
        {
            "criterion": "3_chemical_purity",
            "description": "No Ti, Mn, Cr impurities",
            "status_in_ALH84001": "PRESENT in THO subset",
            "abiotic_explained": False,
            "note": "High purity not fully explained abiotically "
                    "under Mars conditions. "
                    "Thomas-Keprta et al. 2009."
        },
        {
            "criterion": "4_crystal_perfection",
            "description": "Defect-free crystal structure",
            "status_in_ALH84001": "PRESENT in THO subset",
            "abiotic_explained": True,
            "note": "Defect-free crystals can form abiotically "
                    "at elevated temperatures."
        },
        {
            "criterion": "5_chain_aggregates",
            "description": "Some crystals in chain organization",
            "status_in_ALH84001": "PARTIAL — some chains observed",
            "abiotic_explained": False,
            "note": "Chain organization is hallmark of magnetosomes. "
                    "Weakly observed in ALH84001."
        },
        {
            "criterion": "6_no_FeSilicate_rims",
            "description": "No Fe-silicate rims on crystals",
            "status_in_ALH84001": "PRESENT",
            "abiotic_explained": True,
            "note": "Rim absence consistent with both "
                    "biological and some abiotic formation."
        },
    ],
    "thomas_keprta_conclusion": "All six criteria met in THO population. "
                                "Most parsimonious: biogenic.",
    "golden_counter": "Abiotic thermal decomposition can produce "
                      "SOME criteria. Not all six simultaneously.",
    "current_status": "UNRESOLVED — debate live as of 2009. "
                      "Neither side retracted."
}

# Nakhla carbon isotope data
# Source: Steele et al. 2012 (Science)
# This tests whether indigenous Martian organic carbon
# shows biological isotopic fractionation direction.
NAKHLA_CARBON = {
    "meteorite": "Nakhla",
    "instrument": "NanoSIMS / Raman (Steele et al. 2012)",
    "delta_13C_values": {
        "organic_carbon_indigenous": (-22, +2),   # range in permil
        "Martian_atmosphere_CO2": (+41, +45),      # published range
        "Earth_biological_organic": (-30, -10),    # reference range
        "abiotic_serpentinization": (-10, -40),    # reference range
        "Curiosity_SAM_mudstone": -137,            # House et al. 2022
    },
    "steele_conclusion": "Organic matter in Nakhla is INDIGENOUS "
                         "to Mars. Not modern terrestrial contamination.",
    "isotope_note": "δ¹³C -22 to +2‰ is 63-43‰ LIGHTER than "
                    "Martian atmospheric CO₂ (+41 to +45‰). "
                    "This is the biological fractionation direction. "
                    "Organisms preferentially incorporate ¹²C. "
                    "The organic carbon is depleted in ¹³C relative "
                    "to the Martian carbon reservoir. "
                    "This is consistent with but does not prove "
                    "biological origin.",
    "curiosity_note": "SAM δ¹³C = -137‰ in Gale Crater mudstone. "
                      "178‰ lighter than Martian atmospheric CO₂. "
                      "House et al. 2022 cannot eliminate biological "
                      "explanation. Three scenarios, none ruled out."
}

# ─────────────────────────────────────────────────────────────
# SECTION 3: THE HYPOTHESIS PREDICTION
# ─────────────────────────────────────────────────────────────

HYPOTHESIS = {
    "name": "Mars Iron Necromass Hypothesis",
    "pre_registration": "DOI: 10.5281/zenodo.18986790",
    "date": "2026-03-13",
    "prediction": (
        "Iron oxide content in Martian subsurface rock "
        "increases monotonically with depth in the original "
        "Mars rock pile."
    ),
    "mechanism": (
        "Subsurface biological iron-cycling community operates "
        "deeper in the Martian crust. Iron-oxidising organisms "
        "mine Fe²⁺ from olivine and basalt, precipitate Fe³⁺ "
        "as iron oxide in fracture networks. Deeper rock spent "
        "more time in the biological iron-cycling zone before "
        "ejection. Therefore deeper = more iron oxide deposited."
    ),
    "null_hypothesis": (
        "Iron oxide content is uniform with depth, or decreases "
        "with depth. Consistent with abiotic top-down surface "
        "weathering which operates from the surface downward."
    ),
    "test_variable_x": "Estimated depth in Mars rock pile (metres)",
    "test_variable_y": "Secondary vein FeO content (wt%)",
    "secondary_test_y": "Secondary alteration volume fraction (%)",
    "expected_result_if_biological": (
        "Positive correlation: deeper samples have higher "
        "vein FeO content AND higher alteration volume."
    ),
    "expected_result_if_abiotic": (
        "No correlation or negative correlation: "
        "deeper samples have same or lower FeO content."
    ),
    "confounds": [
        "Water/rock ratio varies with depth independently.",
        "Temperature varies with depth — affects mineral stability.",
        "The nakhlite pile experienced a single aqueous event "
        "(estimated ~600 Ma ago), not continuous water flow.",
        "Depth estimates have large uncertainties (see ranges).",
        "FeO values from alteration veins, not bulk rock.",
        "Sample sizes are small (5 meteorites, some paired).",
        "Governador Valadares and Nakhla are from same depth "
        "— they are not independent depth points.",
    ],
}

# ─────────────────────────────────────────────────────────────
# SECTION 4: DATA QUALITY FLAGS
# ─────────────────────────────────────────────────────────────

DATA_QUALITY = {
    "nakhlite_depth_uncertainty": (
        "LARGE. Depth estimates from mineral zoning and "
        "secondary alteration degree. Not directly measured. "
        "Yamato: 0-7m (high uncertainty). "
        "Nakhla/GV: 7-10m (moderate). "
        "Lafayette/NWA998: 20-40m (high). "
        "The analysis script must propagate this uncertainty."
    ),
    "vein_FeO_values": (
        "These are representative published values from "
        "Changela & Bridges 2010 Table 4. "
        "They are measurements of specific vein phases "
        "in specific thin sections. "
        "They are NOT averages of the whole meteorite. "
        "The full analysis script should load all "
        "published vein measurements, not just the "
        "representative values encoded here."
    ),
    "alteration_volume": (
        "Alteration volume percentage estimates are rough. "
        "Published values vary by study. "
        "These are order-of-magnitude estimates only."
    ),
    "sample_size": (
        "CRITICAL LIMITATION: Only 5 meteorites. "
        "Nakhla and GV are at the same depth. "
        "Lafayette and NWA998 are at the same depth. "
        "Effective independent depth points: 3. "
        "Statistical significance will be limited. "
        "The analysis must be honest about this."
    ),
    "contamination": (
        "Mars meteorites spent time on Earth or in space "
        "before collection. Terrestrial iron-containing "
        "dust can coat samples. However: the VEIN iron oxide "
        "is enclosed inside the rock interior. "
        "Surface contamination affects the exterior. "
        "Interior vein composition is more reliable "
        "as a record of Martian processes."
    ),
}

# ─────────────────────────────────────────────────────────────
# SECTION 5: SIMPLE PEARSON CORRELATION — DIAGNOSTIC ONLY
# ─────────────────────────────────────────────────────────────
# We compute the correlation between depth and FeO here
# using only the standard library and the data above.
# This is the diagnostic check — not the final analysis.
# The analysis script will use scipy for proper statistics.

def mean(values):
    return sum(values) / len(values)

def pearson_r(x, y):
    """Pearson correlation coefficient without scipy."""
    n = len(x)
    if n != len(y) or n < 2:
        return None
    mx = mean(x)
    my = mean(y)
    numerator = sum((x[i] - mx) * (y[i] - my) for i in range(n))
    denom_x = math.sqrt(sum((x[i] - mx)**2 for i in range(n)))
    denom_y = math.sqrt(sum((y[i] - my)**2 for i in range(n)))
    if denom_x == 0 or denom_y == 0:
        return None
    return numerator / (denom_x * denom_y)

def spearman_r(x, y):
    """Spearman rank correlation without scipy."""
    n = len(x)
    def rank(arr):
        sorted_arr = sorted(enumerate(arr), key=lambda t: t[1])
        ranks = [0] * n
        for rank_val, (orig_idx, _) in enumerate(sorted_arr):
            ranks[orig_idx] = rank_val + 1
        return ranks
    rx = rank(x)
    ry = rank(y)
    return pearson_r(rx, ry)

# ─────────────────────────────────────────────────────────────
# SECTION 6: DIAGNOSTIC REPORT
# ─────────────────────────────────────────────────────────────

def print_separator(char="─", width=60):
    print(char * width)

def print_header(title):
    print()
    print_separator("═")
    print(f"  {title}")
    print_separator("═")

def run_diagnostic():

    print_header("NECROMASS DIAGNOSTIC REPORT v1.0")
    print(f"  Date: 2026-03-13")
    print(f"  Pre-registration: {HYPOTHESIS['pre_registration']}")
    print(f"  Author: Eric Robert Lawson / OrganismCore")

    # ── SECTION 1: DATA INVENTORY ────────────────────────────
    print_header("SECTION 1: DATA INVENTORY")
    print("What Martian material we have and what it contains:")
    print()

    print("  [A] NAKHLITE METEORITES — 5 meteorites from")
    print("      the same original Mars rock pile,")
    print("      at three distinct depth levels.")
    print()
    for m in NAKHLITE_DATA:
        print(f"  Meteorite   : {m['name']}")
        print(f"  Depth (m)   : {m['depth_m']} "
              f"[range: {m['depth_m_min']}–{m['depth_m_max']}]")
        print(f"  Vein FeO    : {m['vein_FeO_mean']} wt% "
              f"[{m['vein_FeO_min']}–{m['vein_FeO_max']}]")
        print(f"  Alt. volume : {m['alteration_vol_pct']}%")
        print(f"  Source      : {m['source']}")
        print(f"  Note        : {m['depth_note']}")
        print()

    print("  [B] ALH84001 — Ancient Mars crust (~4.1 Ga)")
    print("      Contains magnetite with disputed biological")
    print("      crystal morphology. Six-criteria analysis")
    print("      by Thomas-Keprta et al. 2001, 2009.")
    print()
    criteria_met = sum(
        1 for c in ALH84001_MAGNETITE["criteria"]
        if "PRESENT" in c["status_in_ALH84001"]
    )
    criteria_abiotic = sum(
        1 for c in ALH84001_MAGNETITE["criteria"]
        if c["abiotic_explained"]
    )
    print(f"  Criteria present in ALH84001 : "
          f"{criteria_met}/6")
    print(f"  Criteria explained abiotically: "
          f"{criteria_abiotic}/6")
    print(f"  Criteria NOT explained abiotically: "
          f"{criteria_met - criteria_abiotic}/6")
    print(f"  Current status: "
          f"{ALH84001_MAGNETITE['current_status']}")
    print()

    print("  [C] NAKHLA CARBON ISOTOPES")
    d = NAKHLA_CARBON["delta_13C_values"]
    print(f"  Nakhla indigenous organic δ¹³C: "
          f"{d['organic_carbon_indigenous'][0]} to "
          f"{d['organic_carbon_indigenous'][1]}‰")
    print(f"  Martian atmospheric CO₂ δ¹³C : "
          f"{d['Martian_atmosphere_CO2'][0]} to "
          f"{d['Martian_atmosphere_CO2'][1]}‰")
    depletion = (d['Martian_atmosphere_CO2'][0] -
                 d['organic_carbon_indigenous'][1])
    print(f"  Depletion vs. atm CO₂        : "
          f"~{depletion}‰ lighter (minimum)")
    print(f"  Steele 2012 conclusion        : "
          f"{NAKHLA_CARBON['steele_conclusion']}")
    print()
    print(f"  Curiosity SAM δ¹³C (mudstone) : "
          f"{d['Curiosity_SAM_mudstone']}‰")
    print(f"  Depletion vs. atm CO₂        : "
          f"~{d['Martian_atmosphere_CO2'][0] - d['Curiosity_SAM_mudstone']}‰")
    print()

    # ── SECTION 2: THE TEST ──────────────────────────────────
    print_header("SECTION 2: THE SPECIFIC TEST")
    print()
    print(f"  Hypothesis : {HYPOTHESIS['name']}")
    print(f"  Prediction : {HYPOTHESIS['prediction']}")
    print()
    print(f"  X variable : {HYPOTHESIS['test_variable_x']}")
    print(f"  Y variable : {HYPOTHESIS['test_variable_y']}")
    print(f"  Y secondary: {HYPOTHESIS['secondary_test_y']}")
    print()
    print(f"  If BIOLOGICAL : {HYPOTHESIS['expected_result_if_biological']}")
    print(f"  If ABIOTIC    : {HYPOTHESIS['expected_result_if_abiotic']}")
    print()

    # ── SECTION 3: DIAGNOSTIC CORRELATION ───────────────────
    print_header("SECTION 3: DIAGNOSTIC CORRELATION (NO SCIPY)")
    print()
    print("  Using published representative mean values only.")
    print("  This is NOT the final statistical test.")
    print("  This tells us: is the DATA STRUCTURED the way")
    print("  the hypothesis predicts before we do real stats?")
    print()

    depths = [m["depth_m"] for m in NAKHLITE_DATA]
    feo_vals = [m["vein_FeO_mean"] for m in NAKHLITE_DATA]
    alt_vols = [m["alteration_vol_pct"] for m in NAKHLITE_DATA]

    print("  Depth vs Vein FeO (wt%) — all 5 meteorites:")
    print()
    print(f"  {'Meteorite':<20} {'Depth(m)':>10} "
          f"{'VeinFeO(wt%)':>14} {'Alt.Vol(%)':>12}")
    print_separator()
    for m in NAKHLITE_DATA:
        print(f"  {m['name_short']:<20} {m['depth_m']:>10.1f} "
              f"{m['vein_FeO_mean']:>14.1f} "
              f"{m['alteration_vol_pct']:>12.2f}")
    print()

    r_feo = pearson_r(depths, feo_vals)
    rho_feo = spearman_r(depths, feo_vals)
    r_alt = pearson_r(depths, alt_vols)
    rho_alt = spearman_r(depths, alt_vols)

    print(f"  Pearson r  (depth vs vein FeO)   : {r_feo:+.4f}")
    print(f"  Spearman ρ (depth vs vein FeO)   : {rho_feo:+.4f}")
    print(f"  Pearson r  (depth vs alt. volume): {r_alt:+.4f}")
    print(f"  Spearman ρ (depth vs alt. volume): {rho_alt:+.4f}")
    print()

    # Interpret diagnostic result
    print("  DIAGNOSTIC INTERPRETATION:")
    print()
    if r_feo is not None:
        if r_feo > 0.3:
            print(f"  ► Vein FeO vs depth: POSITIVE correlation "
                  f"(r={r_feo:+.3f})")
            print(f"    Direction is CONSISTENT WITH biological")
            print(f"    hypothesis. Deeper = more iron oxide.")
        elif r_feo < -0.3:
            print(f"  ► Vein FeO vs depth: NEGATIVE correlation "
                  f"(r={r_feo:+.3f})")
            print(f"    Direction CONTRADICTS biological hypothesis.")
            print(f"    Deeper = less iron oxide.")
        else:
            print(f"  ► Vein FeO vs depth: WEAK/NO correlation "
                  f"(r={r_feo:+.3f})")
            print(f"    Cannot distinguish biological from abiotic.")
    print()
    if r_alt is not None:
        if r_alt > 0.3:
            print(f"  ► Alteration vol vs depth: POSITIVE correlation "
                  f"(r={r_alt:+.3f})")
            print(f"    Direction is CONSISTENT WITH biological")
            print(f"    hypothesis. Deeper = more total alteration.")
        elif r_alt < -0.3:
            print(f"  ► Alteration vol vs depth: NEGATIVE correlation "
                  f"(r={r_alt:+.3f})")
            print(f"    CONTRADICTS biological hypothesis.")
        else:
            print(f"  ► Alteration vol vs depth: WEAK/NO correlation "
                  f"(r={r_alt:+.3f})")
    print()

    # ── SECTION 4: DATA QUALITY AUDIT ───────────────────────
    print_header("SECTION 4: DATA QUALITY AUDIT")
    print()
    for key, val in DATA_QUALITY.items():
        print(f"  [{key.upper()}]")
        # Word wrap the value
        words = val.split()
        line = "    "
        for word in words:
            if len(line) + len(word) + 1 > 65:
                print(line)
                line = "    " + word + " "
            else:
                line += word + " "
        if line.strip():
            print(line)
        print()

    # ── SECTION 5: CONFOUNDS AUDIT ───────────────────────────
    print_header("SECTION 5: CONFOUNDS THAT MUST BE")
    print("            ADDRESSED IN ANALYSIS SCRIPT")
    print()
    for i, c in enumerate(HYPOTHESIS["confounds"], 1):
        print(f"  {i}. {c}")
    print()

    # ── SECTION 6: WHAT THE ANALYSIS SCRIPT NEEDS ───────────
    print_header("SECTION 6: WHAT THE ANALYSIS SCRIPT NEEDS")
    print()
    needs = [
        ("DATA SOURCE 1",
         "Full published Table 4 from Changela & Bridges 2010. "
         "All individual vein FeO measurements, not just means. "
         "Requires reading the paper and transcribing values. "
         "DOI: 10.1111/j.1945-5100.2010.01123.x"),
        ("DATA SOURCE 2",
         "Alteration volume fraction data with error bars. "
         "Best source: Changela & Bridges 2010, Section 3. "
         "Supplemented by Gooding et al. 1991."),
        ("STATISTICS",
         "scipy.stats.pearsonr() and spearmanr() with p-values. "
         "scipy.stats.linregress() for regression line. "
         "Bootstrap resampling for confidence intervals given "
         "small n. Must report n=5 (or n=3 independent points) "
         "prominently."),
        ("VISUALISATION",
         "matplotlib scatter plot: depth (x) vs vein FeO (y). "
         "Error bars on both axes from published ranges. "
         "Regression line with 95% confidence interval. "
         "Separate panel: alteration volume vs depth. "
         "Must be publication-quality for pre-registration."),
        ("UNCERTAINTY PROPAGATION",
         "Monte Carlo simulation using depth ranges "
         "[min, max] to test whether correlation holds across "
         "all plausible depth assignments. "
         "If correlation is robust to depth uncertainty: "
         "strong result. If it flips: weak result."),
        ("ALH84001 MODULE",
         "Separate analysis: tabulate six criteria, "
         "abiotic explanation status for each, "
         "produce summary table for paper. "
         "No statistics needed — this is a qualitative "
         "consistency check."),
        ("ISOTOPE MODULE",
         "Plot δ¹³C values: Nakhla organic, Martian CO₂, "
         "Earth biological, serpentinization, Curiosity SAM. "
         "Show where each falls on the isotope scale. "
         "Annotate which direction biological fractionation "
         "predicts. Annotate the -137‰ anomaly."),
        ("HONEST LIMITATION SECTION",
         "The script must print a clear statement that "
         "n=5 meteorites is insufficient for statistical "
         "significance. The correlation is directional "
         "evidence only, not proof. Must be stated in output."),
    ]
    for label, desc in needs:
        print(f"  [{label}]")
        words = desc.split()
        line = "    "
        for word in words:
            if len(line) + len(word) + 1 > 65:
                print(line)
                line = "    " + word + " "
            else:
                line += word + " "
        if line.strip():
            print(line)
        print()

    # ── SECTION 7: GO / NO-GO DECISION ───────────────────────
    print_header("SECTION 7: GO / NO-GO DECISION")
    print()
    print("  Before writing the analysis script, answer:")
    print()
    checks = [
        ("Is the depth sequence confirmed?",
         "YES. Yamato(<7m) < Nakhla/GV(7-10m) < "
         "Lafayette/NWA998(20-40m). "
         "Published in Changela & Bridges 2010 "
         "and Mikouchi et al. 2006."),
        ("Is the FeO data real published values?",
         "PARTIALLY. Representative values are from "
         "Changela & Bridges 2010 Table 4. "
         "Full table must be transcribed from the paper. "
         "Do not run analysis on representative means alone."),
        ("Is the test direction clear?",
         "YES. Deeper = more FeO is biological prediction. "
         "Same or less FeO with depth = abiotic prediction. "
         "These are mutually exclusive directional predictions."),
        ("Is the sample size sufficient for statistics?",
         "NO. n=5, n=3 independent depth points. "
         "No p-value will be significant. "
         "The analysis is directional evidence only. "
         "This must be stated explicitly."),
        ("Does the diagnostic correlation support proceeding?",
         f"PENDING — see Section 3 output above. "
         f"Check the sign of r_feo. "
         f"If positive: proceed. If negative: investigate "
         f"before proceeding."),
        ("Is the data openly available?",
         "YES. Changela & Bridges 2010 is available via "
         "DOI: 10.1111/j.1945-5100.2010.01123.x "
         "NASA JSC Compendium: "
         "curator.jsc.nasa.gov/antmet/mmc/ "
         "Astromat: astromat.org/collections/meteorites"),
    ]
    for q, a in checks:
        print(f"  Q: {q}")
        words = a.split()
        line = "  A: "
        for word in words:
            if len(line) + len(word) + 1 > 65:
                print(line)
                line = "     " + word + " "
            else:
                line += word + " "
        if line.strip():
            print(line)
        print()

    # ── SECTION 8: NEXT STEPS ────────────────────────────────
    print_header("SECTION 8: NEXT STEPS")
    print()
    steps = [
        "1. READ THIS OUTPUT IN FULL.",
        "2. Check the diagnostic correlation sign in Section 3.",
        "   If positive: the data direction is consistent with",
        "   the hypothesis. Proceed to analysis script.",
        "   If negative: the data contradicts the hypothesis.",
        "   Investigate before proceeding.",
        "3. Access Changela & Bridges 2010 full Table 4.",
        "   Transcribe all individual vein FeO values.",
        "   Add them to the data module in the analysis script.",
        "4. Confirm depth estimates against multiple sources.",
        "   Use depth range midpoints as primary values.",
        "   Use min/max for Monte Carlo uncertainty.",
        "5. Write analysis script using this diagnostic as",
        "   the structural template.",
        "6. Pre-register the analysis script on Zenodo",
        "   BEFORE running it against the full dataset.",
        "7. Run. Report the result regardless of outcome.",
    ]
    for step in steps:
        print(f"  {step}")
    print()

    print_separator("═")
    print("  DIAGNOSTIC COMPLETE.")
    print("  Read all output before proceeding.")
    print("  The analysis script is only written after")
    print("  this diagnostic is understood.")
    print_separator("═")
    print()

    # ── RETURN MACHINE-READABLE SUMMARY ────────────────────��
    return {
        "pearson_r_depth_vs_vein_FeO": r_feo,
        "spearman_r_depth_vs_vein_FeO": rho_feo,
        "pearson_r_depth_vs_alteration_vol": r_alt,
        "spearman_r_depth_vs_alteration_vol": rho_alt,
        "n_meteorites": len(NAKHLITE_DATA),
        "n_independent_depth_points": 3,
        "hypothesis_direction_confirmed": (
            r_feo is not None and r_feo > 0
        ),
        "data_quality_flag": "REPRESENTATIVE VALUES ONLY — "
                             "LOAD FULL TABLE BEFORE ANALYSIS",
    }


# ─────────────────────────────────────────────────────────────
# ENTRY POINT
# ──────────��──────────────────────────────────────────────────

if __name__ == "__main__":
    result = run_diagnostic()
    print()
    print("  MACHINE-READABLE SUMMARY:")
    print()
    for k, v in result.items():
        print(f"  {k}: {v}")
    print()
    sys.exit(0)
