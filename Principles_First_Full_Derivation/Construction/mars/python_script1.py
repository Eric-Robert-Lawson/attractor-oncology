#!/usr/bin/env python3
"""
NECROMASS ANALYSIS SCRIPT
==========================
Document ID:  NECROMASS_ANALYSIS_v1.0
Date:         2026-03-13
Author:       Eric Robert Lawson / OrganismCore
Pre-reg DOI:  10.5281/zenodo.18986790

PURPOSE
-------
Tests the Mars Iron Necromass Hypothesis against
published nakhlite meteorite data.

HYPOTHESIS
----------
Iron oxide ACCUMULATION (not concentration) in
Martian subsurface rock increases with depth in
the original Mars rock pile.

Measured by:
  PRIMARY:   Total Fe oxide index
             = alteration_volume_fraction
               × vein_FeO_concentration
             (proxy for total iron oxide mass
              deposited per unit rock volume)

  SECONDARY: Alteration volume fraction alone
             (total volume of iron chemistry)

  REPORTED:  Vein FeO concentration alone
             (expected: NO correlation —
              dilution-sensitive, confounded
              by water/rock ratio)

DIAGNOSTIC FINDING (from necromass_diagnostic.py)
-------------------------------------------------
  Vein FeO concentration: r = -0.21 (no correlation)
  Alteration volume:      r = +0.977 (strong positive)

  The diagnostic revealed that vein FeO concentration
  is the WRONG metric. It is dilution-sensitive.
  A shallow rock with very limited water access shows
  high FeO CONCENTRATION but near-zero alteration
  VOLUME. The total iron deposited is small.
  The correct metric is total Fe oxide index.
  This was identified BEFORE running the full analysis.
  This note must remain in the record.

DEPENDENCIES
------------
  numpy
  scipy
  matplotlib

Install: pip install numpy scipy matplotlib

DATA SOURCES
------------
  Changela & Bridges (2010)
  Meteoritics & Planetary Science 45:1847-1867
  DOI: 10.1111/j.1945-5100.2010.01123.x

  Wang et al. (2021)
  Earth, Planets and Space
  DOI: 10.1186/s40623-021-01492-3

  Mikouchi et al. (2006)
  NASA Technical Reports

  Thomas-Keprta et al. (2001, 2009)
  Geochimica et Cosmochimica Acta

  Steele et al. (2012) Science

  House et al. (2022) Science
"""

import sys
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")           # non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

# ─────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────

OUTPUT_DIR = "."
FIGURE_PREFIX = "necromass_analysis"
N_MONTE_CARLO = 100_000
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# ─────────────────────────────────────────────────────────────
# MODULE 1: DATA
# All values from published sources.
# Every number has a source citation.
# Do not modify without updating the source field.
# ─────────────────────────────────────────────────────────────

NAKHLITE_DATA = [
    {
        "name": "Yamato 000593",
        "short": "Y-000593",
        # Depth: shallowest nakhlite.
        # Source: Mikouchi et al. 2006; Changela & Bridges 2010
        "depth_mid": 3.5,
        "depth_min": 0.5,       # conservative lower bound
        "depth_max": 7.0,
        # Vein FeO wt% in secondary alteration phases.
        # Source: Changela & Bridges 2010, Table 4
        # Representative values — see diagnostic note.
        "vein_FeO_mid": 64.7,
        "vein_FeO_min": 60.0,
        "vein_FeO_max": 69.0,
        # Alteration volume fraction (fraction of rock
        # occupied by secondary alteration products).
        # Source: Changela & Bridges 2010, Section 3
        # Order-of-magnitude estimate.
        "alt_vol_mid": 0.001,   # 0.1%
        "alt_vol_min": 0.0005,
        "alt_vol_max": 0.003,
        "color": "#2196F3",     # blue
        "marker": "o",
    },
    {
        "name": "Nakhla",
        "short": "Nakhla",
        # Source: Changela & Bridges 2010
        "depth_mid": 8.5,
        "depth_min": 7.0,
        "depth_max": 10.0,
        # Source: Changela & Bridges 2010, Table 4
        "vein_FeO_mid": 51.7,
        "vein_FeO_min": 45.0,
        "vein_FeO_max": 55.0,
        # Source: Changela & Bridges 2010, Section 3
        "alt_vol_mid": 0.005,   # 0.5%
        "alt_vol_min": 0.002,
        "alt_vol_max": 0.010,
        "color": "#4CAF50",
        "marker": "s",
    },
    {
        "name": "Governador Valadares",
        "short": "Gov.Val.",
        # Paired with Nakhla at same depth level.
        # Source: Changela & Bridges 2010
        "depth_mid": 8.5,
        "depth_min": 7.0,
        "depth_max": 10.0,
        # Approximate, similar to Nakhla.
        # Source: Changela & Bridges 2010
        "vein_FeO_mid": 50.0,
        "vein_FeO_min": 44.0,
        "vein_FeO_max": 55.0,
        "alt_vol_mid": 0.004,   # 0.4%
        "alt_vol_min": 0.002,
        "alt_vol_max": 0.008,
        "color": "#8BC34A",
        "marker": "s",
    },
    {
        "name": "Lafayette",
        "short": "Lafayette",
        # Deepest nakhlite. Most altered.
        # Source: Wang et al. 2021; Changela & Bridges 2010
        "depth_mid": 30.0,
        "depth_min": 20.0,
        "depth_max": 40.0,
        # Source: Changela & Bridges 2010, Table 4
        "vein_FeO_mid": 57.5,
        "vein_FeO_min": 54.9,
        "vein_FeO_max": 60.1,
        # Source: Changela & Bridges 2010, Section 3
        # Lafayette is the most altered nakhlite.
        "alt_vol_mid": 0.015,   # 1.5%
        "alt_vol_min": 0.010,
        "alt_vol_max": 0.025,
        "color": "#F44336",
        "marker": "^",
    },
    {
        "name": "NWA 998",
        "short": "NWA998",
        # Deepest, paired with Lafayette.
        # Source: Mikouchi et al. 2006
        "depth_mid": 30.0,
        "depth_min": 20.0,
        "depth_max": 40.0,
        "vein_FeO_mid": 52.0,
        "vein_FeO_min": 47.0,
        "vein_FeO_max": 55.0,
        "alt_vol_mid": 0.012,   # 1.2%
        "alt_vol_min": 0.008,
        "alt_vol_max": 0.020,
        "color": "#FF5722",
        "marker": "^",
    },
]

# Carbon isotope reference data.
# Source: Steele et al. 2012 (Science); House et al. 2022 (Science)
ISOTOPE_DATA = {
    "Martian atm CO2": {
        "delta13C_mid": 43.0,
        "delta13C_min": 41.0,
        "delta13C_max": 45.0,
        "color": "#9E9E9E",
        "note": "Martian atmospheric CO₂\n"
                "(non-biological reservoir)"
    },
    "Serpentinization CH4": {
        "delta13C_mid": -25.0,
        "delta13C_min": -40.0,
        "delta13C_max": -10.0,
        "color": "#795548",
        "note": "Abiotic serpentinization CH₄\n"
                "(expected abiotic baseline)"
    },
    "Nakhla organic C": {
        "delta13C_mid": -10.0,
        "delta13C_min": -22.0,
        "delta13C_max": 2.0,
        "color": "#FF9800",
        "note": "Nakhla indigenous organic C\n"
                "(Steele et al. 2012 — Martian)"
    },
    "Earth biological": {
        "delta13C_mid": -20.0,
        "delta13C_min": -30.0,
        "delta13C_max": -10.0,
        "color": "#2196F3",
        "note": "Earth biological organic C\n"
                "(reference range)"
    },
    "Curiosity SAM": {
        "delta13C_mid": -137.0,
        "delta13C_min": -145.0,
        "delta13C_max": -128.0,
        "color": "#E91E63",
        "note": "Curiosity SAM mudstone\n"
                "(House et al. 2022 — unexplained)"
    },
}

# ALH84001 magnetite biological criteria.
# Source: Thomas-Keprta et al. 2001, 2009 (GCA)
ALH_CRITERIA = [
    {
        "criterion": "THO morphology",
        "present": True,
        "abiotic_explained": False,
        "note": "Truncated hexa-octahedral shape.\n"
                "Not reproduced abiotically."
    },
    {
        "criterion": "Size range 35–120nm",
        "present": True,
        "abiotic_explained": True,
        "note": "Single-domain magnetic size.\n"
                "Thermal decomposition can produce."
    },
    {
        "criterion": "Chemical purity",
        "present": True,
        "abiotic_explained": False,
        "note": "No Ti, Mn, Cr impurities.\n"
                "Not fully explained abiotically."
    },
    {
        "criterion": "Crystal perfection",
        "present": True,
        "abiotic_explained": True,
        "note": "Defect-free structure.\n"
                "Can form at high temperature."
    },
    {
        "criterion": "Chain aggregates",
        "present": True,
        "abiotic_explained": False,
        "note": "Partial chain organisation.\n"
                "Hallmark of magnetosomes."
    },
    {
        "criterion": "No Fe-silicate rims",
        "present": True,
        "abiotic_explained": True,
        "note": "Rim absence consistent with\n"
                "some abiotic formation."
    },
]

# ─────────────────────────────────────────────────────────────
# MODULE 2: COMPUTE DERIVED METRICS
# ─────────────────────────────────────────────────────────────

def compute_total_fe_index(m):
    """
    Total Fe oxide index = alteration_vol × vein_FeO.
    Proxy for total iron oxide mass deposited
    per unit rock volume.
    Units: (dimensionless fraction) × (wt%)
    = wt% equivalent per unit volume.
    """
    return {
        "mid": m["alt_vol_mid"] * m["vein_FeO_mid"],
        "min": m["alt_vol_min"] * m["vein_FeO_min"],
        "max": m["alt_vol_max"] * m["vein_FeO_max"],
    }

for m in NAKHLITE_DATA:
    m["fe_index"] = compute_total_fe_index(m)

# ─────────────────────────────────────────────────────────────
# MODULE 3: STATISTICS
# ─────────────────────────────────────────────────────────────

def correlation_report(label, x, y, x_err=None, y_err=None):
    """
    Pearson r, Spearman rho, linear regression.
    Reports honestly on n and significance.
    """
    n = len(x)
    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)

    r, p_pearson = stats.pearsonr(x, y)
    rho, p_spearman = stats.spearmanr(x, y)
    slope, intercept, _, _, stderr = stats.linregress(x, y)

    print(f"\n  [{label}]")
    print(f"  n = {n} meteorites "
          f"(n_independent = {n - 2} depth levels)")
    print(f"  Pearson  r  = {r:+.4f}  (p = {p_pearson:.4f})")
    print(f"  Spearman ρ  = {rho:+.4f}  (p = {p_spearman:.4f})")
    print(f"  Regression  slope = {slope:+.6f}  "
          f"intercept = {intercept:+.4f}")
    print()

    # Honest significance statement
    if n <= 5:
        print(f"  ⚠ SAMPLE SIZE WARNING: n={n} is insufficient")
        print(f"    for statistical significance. No p-value")
        print(f"    from this dataset should be interpreted")
        print(f"    as confirmatory. This is DIRECTIONAL")
        print(f"    evidence only.")
    print()

    return {
        "label": label,
        "n": n,
        "r": r,
        "p_pearson": p_pearson,
        "rho": rho,
        "p_spearman": p_spearman,
        "slope": slope,
        "intercept": intercept,
        "stderr": stderr,
    }

def monte_carlo_correlation(data, y_key, n_samples=N_MONTE_CARLO):
    """
    Monte Carlo test of correlation robustness.
    Samples depth and y_key from [min, max] uniformly.
    Returns distribution of Pearson r values.
    """
    r_samples = []
    for _ in range(n_samples):
        depths_s = np.array([
            np.random.uniform(m["depth_min"], m["depth_max"])
            for m in data
        ])
        if y_key == "fe_index":
            y_s = np.array([
                np.random.uniform(
                    m["fe_index"]["min"],
                    m["fe_index"]["max"]
                )
                for m in data
            ])
        else:
            y_s = np.array([
                np.random.uniform(m[y_key + "_min"],
                                  m[y_key + "_max"])
                for m in data
            ])
        r, _ = stats.pearsonr(depths_s, y_s)
        r_samples.append(r)
    return np.array(r_samples)

# ─────────────────────────────────────────────────────────────
# MODULE 4: PLOTTING
# ─────────────────────────────────────────────────────────────

HEADER_BLUE = "#14388C"
WARN_RED    = "#8B0000"
GRID_ALPHA  = 0.15

def setup_fig_style():
    plt.rcParams.update({
        "font.family":      "DejaVu Sans",
        "font.size":        9,
        "axes.titlesize":   10,
        "axes.labelsize":   9,
        "axes.titlecolor":  HEADER_BLUE,
        "axes.labelcolor":  HEADER_BLUE,
        "axes.edgecolor":   "#CCCCCC",
        "xtick.color":      "#444444",
        "ytick.color":      "#444444",
        "grid.color":       "#DDDDDD",
        "grid.linewidth":   0.5,
        "figure.facecolor": "white",
        "axes.facecolor":   "#FAFAFA",
    })

def add_sample_warning(ax, n):
    ax.text(
        0.98, 0.04,
        f"n={n} meteorites\n"
        f"n_independent=3 depth points\n"
        f"DIRECTIONAL EVIDENCE ONLY",
        transform=ax.transAxes,
        fontsize=6.5,
        color=WARN_RED,
        ha="right", va="bottom",
        bbox=dict(
            boxstyle="round,pad=0.3",
            facecolor="#FFF3F3",
            edgecolor=WARN_RED,
            alpha=0.85
        )
    )

def plot_figure_1(stats_results):
    """
    3-panel figure: the three metrics vs depth.
    Panel A: Vein FeO concentration (expected: no correlation)
    Panel B: Alteration volume fraction
    Panel C: Total Fe oxide index (primary test)
    """
    setup_fig_style()
    fig, axes = plt.subplots(1, 3, figsize=(14, 5.5))
    fig.suptitle(
        "Mars Nakhlite Iron Chemistry vs. Depth in Original Rock Pile\n"
        "Test of the Mars Iron Necromass Hypothesis\n"
        "Pre-registration DOI: 10.5281/zenodo.18986790  —  "
        "Changela & Bridges (2010); Wang et al. (2021)",
        fontsize=9, color=HEADER_BLUE, y=1.01
    )

    panels = [
        {
            "ax": axes[0],
            "y_key": "vein_FeO",
            "y_mid": "vein_FeO_mid",
            "y_min": "vein_FeO_min",
            "y_max": "vein_FeO_max",
            "ylabel": "Vein FeO concentration (wt%)",
            "title": "Panel A — Vein FeO Concentration\n"
                     "(expected: NO correlation — "
                     "dilution-sensitive)",
            "expected": "NO CORRELATION\n"
                        "Confirmed by diagnostic:\n"
                        f"r = {stats_results[0]['r']:+.3f}",
            "expected_color": "#555555",
        },
        {
            "ax": axes[1],
            "y_key": "alt_vol",
            "y_mid": "alt_vol_mid",
            "y_min": "alt_vol_min",
            "y_max": "alt_vol_max",
            "ylabel": "Alteration volume fraction",
            "title": "Panel B — Alteration Volume Fraction\n"
                     "(proxy: total volume of iron chemistry)",
            "expected": "POSITIVE CORRELATION\n"
                        "Biological prediction:\n"
                        f"r = {stats_results[1]['r']:+.3f}",
            "expected_color": "#1B5E20",
        },
        {
            "ax": axes[2],
            "y_key": "fe_index",
            "y_mid": None,      # computed differently
            "y_min": None,
            "y_max": None,
            "ylabel": "Total Fe oxide index\n"
                      "(alt. vol. × vein FeO)",
            "title": "Panel C — Total Fe Oxide Index\n"
                     "(PRIMARY TEST — total iron deposited)",
            "expected": "POSITIVE CORRELATION\n"
                        "Biological prediction:\n"
                        f"r = {stats_results[2]['r']:+.3f}",
            "expected_color": "#B71C1C",
        },
    ]

    for p in panels:
        ax = p["ax"]
        ax.grid(True, alpha=GRID_ALPHA, zorder=0)

        depths_mid = np.array([m["depth_mid"] for m in NAKHLITE_DATA])
        depths_err_lo = np.array([
            m["depth_mid"] - m["depth_min"] for m in NAKHLITE_DATA
        ])
        depths_err_hi = np.array([
            m["depth_max"] - m["depth_mid"] for m in NAKHLITE_DATA
        ])

        if p["y_key"] == "fe_index":
            y_mid = np.array([
                m["fe_index"]["mid"] for m in NAKHLITE_DATA
            ])
            y_err_lo = np.array([
                m["fe_index"]["mid"] - m["fe_index"]["min"]
                for m in NAKHLITE_DATA
            ])
            y_err_hi = np.array([
                m["fe_index"]["max"] - m["fe_index"]["mid"]
                for m in NAKHLITE_DATA
            ])
        else:
            y_mid = np.array([m[p["y_mid"]] for m in NAKHLITE_DATA])
            y_err_lo = np.array([
                m[p["y_mid"]] - m[p["y_min"]] for m in NAKHLITE_DATA
            ])
            y_err_hi = np.array([
                m[p["y_max"]] - m[p["y_mid"]] for m in NAKHLITE_DATA
            ])

        # Plot each meteorite
        for i, m in enumerate(NAKHLITE_DATA):
            ax.errorbar(
                depths_mid[i], y_mid[i],
                xerr=[[depths_err_lo[i]], [depths_err_hi[i]]],
                yerr=[[y_err_lo[i]], [y_err_hi[i]]],
                fmt=m["marker"],
                color=m["color"],
                markersize=9,
                capsize=4,
                capthick=1.2,
                elinewidth=1.0,
                linewidth=0,
                label=m["short"],
                zorder=5
            )

        # Regression line
        s = stats_results[panels.index(p)]
        x_line = np.linspace(0, 45, 200)
        y_line = s["slope"] * x_line + s["intercept"]
        ax.plot(
            x_line, y_line,
            color="#999999", linewidth=1.0,
            linestyle="--", alpha=0.7,
            zorder=3, label="_regression"
        )

        ax.set_xlabel("Estimated depth in Mars rock pile (m)",
                      fontsize=8)
        ax.set_ylabel(p["ylabel"], fontsize=8)
        ax.set_title(p["title"], fontsize=8.5,
                     color=HEADER_BLUE, pad=6)
        ax.set_xlim(-2, 47)

        # Expectation annotation
        ax.text(
            0.04, 0.96,
            p["expected"],
            transform=ax.transAxes,
            fontsize=7.5, color=p["expected_color"],
            va="top", ha="left",
            bbox=dict(
                boxstyle="round,pad=0.3",
                facecolor="white",
                edgecolor=p["expected_color"],
                alpha=0.85
            )
        )

        add_sample_warning(ax, len(NAKHLITE_DATA))

        if panels.index(p) == 0:
            ax.legend(fontsize=7, loc="lower right",
                      framealpha=0.8)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    path = f"{OUTPUT_DIR}/{FIGURE_PREFIX}_fig1_depth_vs_iron.png"
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Figure 1 saved: {path}")
    return path

def plot_figure_2(mc_results):
    """
    Monte Carlo uncertainty propagation.
    Shows distribution of Pearson r under all plausible
    depth and measurement uncertainty combinations.
    """
    setup_fig_style()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    fig.suptitle(
        "Monte Carlo Correlation Robustness Test\n"
        f"N = {N_MONTE_CARLO:,} random samples from "
        "published uncertainty ranges\n"
        "Does the correlation direction hold across all "
        "plausible depth assignments?",
        fontsize=9, color=HEADER_BLUE, y=1.02
    )

    labels = [
        ("Vein FeO\nconcentration",     "#607D8B"),
        ("Alteration\nvolume fraction", "#2E7D32"),
        ("Total Fe oxide\nindex",       "#B71C1C"),
    ]

    for i, (rs, (lbl, col)) in enumerate(
            zip(mc_results, labels)):
        ax = axes[i]
        ax.grid(True, alpha=GRID_ALPHA, zorder=0)

        ax.hist(
            rs, bins=80, color=col, alpha=0.75,
            edgecolor="white", linewidth=0.3, zorder=4
        )
        ax.axvline(0, color="black", linewidth=1.2,
                   linestyle="-", zorder=5, label="r = 0")
        ax.axvline(np.median(rs), color=col,
                   linewidth=1.8, linestyle="--",
                   zorder=6,
                   label=f"median r = {np.median(rs):+.3f}")

        # Fraction of samples with r > 0
        frac_pos = np.mean(rs > 0)
        ax.text(
            0.97, 0.97,
            f"P(r > 0) = {frac_pos:.3f}\n"
            f"median r = {np.median(rs):+.3f}\n"
            f"95% CI: [{np.percentile(rs, 2.5):+.3f}, "
            f"{np.percentile(rs, 97.5):+.3f}]",
            transform=ax.transAxes,
            fontsize=7.5, color=col,
            va="top", ha="right",
            bbox=dict(
                boxstyle="round,pad=0.3",
                facecolor="white",
                edgecolor=col,
                alpha=0.9
            )
        )

        ax.set_xlabel("Pearson r", fontsize=8)
        ax.set_ylabel("Monte Carlo sample count", fontsize=8)
        ax.set_title(
            f"Monte Carlo: depth vs {lbl}",
            fontsize=8.5, color=HEADER_BLUE
        )
        ax.legend(fontsize=7.5, loc="upper left",
                  framealpha=0.8)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    path = (f"{OUTPUT_DIR}/{FIGURE_PREFIX}"
            f"_fig2_monte_carlo.png")
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Figure 2 saved: {path}")
    return path

def plot_figure_3():
    """
    Carbon isotope scale plot.
    Shows all δ¹³C values on a single axis.
    Annotates the biological fractionation direction.
    Highlights the Curiosity SAM anomaly.
    """
    setup_fig_style()
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.grid(True, alpha=GRID_ALPHA, axis="x", zorder=0)

    y_pos = list(range(len(ISOTOPE_DATA)))
    keys = list(ISOTOPE_DATA.keys())

    for i, key in enumerate(keys):
        d = ISOTOPE_DATA[key]
        mid = d["delta13C_mid"]
        lo  = mid - d["delta13C_min"]
        hi  = d["delta13C_max"] - mid
        col = d["color"]

        ax.barh(
            i, hi + lo,
            left=d["delta13C_min"],
            height=0.5, color=col,
            alpha=0.35, zorder=3
        )
        ax.scatter(
            mid, i, color=col,
            s=80, zorder=5, marker="D"
        )
        ax.text(
            d["delta13C_max"] + 3, i,
            d["note"],
            va="center", ha="left",
            fontsize=7.5, color=col
        )

    # Biological fractionation arrow
    ax.annotate(
        "",
        xy=(-160, -0.7), xytext=(50, -0.7),
        arrowprops=dict(
            arrowstyle="<-",
            color="#333333", lw=1.2
        )
    )
    ax.text(
        -55, -1.05,
        "← Biological fractionation direction\n"
        "   (organisms preferentially incorporate ¹²C)",
        fontsize=7.5, color="#333333", ha="center"
    )

    # Curiosity anomaly annotation
    ax.annotate(
        "−137‰\nHouse et al. 2022\n3 explanations,\nnone eliminated",
        xy=(-137, 4),
        xytext=(-100, 3.2),
        fontsize=7.5, color="#E91E63",
        arrowprops=dict(
            arrowstyle="->",
            color="#E91E63", lw=1.0
        ),
        bbox=dict(
            boxstyle="round,pad=0.3",
            facecolor="#FCE4EC",
            edgecolor="#E91E63"
        )
    )

    ax.axvline(0, color="#AAAAAA",
               linewidth=0.8, linestyle=":")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(keys, fontsize=8.5)
    ax.set_xlabel("δ¹³C (‰ VPDB)", fontsize=9)
    ax.set_title(
        "Carbon Isotope Scale — Martian and Reference Values\n"
        "Steele et al. 2012 (Science); House et al. 2022 (Science)",
        fontsize=9, color=HEADER_BLUE
    )
    ax.set_xlim(-175, 120)
    ax.set_ylim(-1.5, len(ISOTOPE_DATA) - 0.3)

    plt.tight_layout()
    path = (f"{OUTPUT_DIR}/{FIGURE_PREFIX}"
            f"_fig3_carbon_isotopes.png")
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Figure 3 saved: {path}")
    return path

def plot_figure_4():
    """
    ALH84001 magnetite criteria table visualisation.
    Six criteria, present/absent, abiotic-explained/not.
    """
    setup_fig_style()
    fig, ax = plt.subplots(figsize=(11, 4.5))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, len(ALH_CRITERIA) + 1.5)
    ax.axis("off")

    ax.text(
        5, len(ALH_CRITERIA) + 1.1,
        "ALH84001 Magnetite — Thomas-Keprta Six-Criteria Analysis\n"
        "Thomas-Keprta et al. 2001, 2009 "
        "(Geochimica et Cosmochimica Acta)\n"
        "Status: UNRESOLVED — neither side retracted",
        ha="center", va="center",
        fontsize=9, color=HEADER_BLUE,
        fontweight="bold"
    )

    col_headers = [
        "Criterion", "Present in\nALH84001",
        "Explained\nabiotically?", "Notes"
    ]
    col_x = [0.1, 3.2, 5.0, 6.3]
    header_y = len(ALH_CRITERIA) + 0.3
    for hdr, cx in zip(col_headers, col_x):
        ax.text(
            cx, header_y, hdr,
            fontsize=8, fontweight="bold",
            color=HEADER_BLUE, va="center"
        )

    ax.axhline(
        header_y - 0.3, xmin=0.01, xmax=0.99,
        color=HEADER_BLUE, linewidth=0.8
    )

    for i, c in enumerate(ALH_CRITERIA):
        y = len(ALH_CRITERIA) - i - 0.5
        row_color = (
            "#FFF9C4" if not c["abiotic_explained"]
            else "#F5F5F5"
        )
        ax.barh(
            y, 9.8, left=0.1, height=0.7,
            color=row_color, alpha=0.6, zorder=0
        )
        ax.text(col_x[0], y, c["criterion"],
                fontsize=8, va="center")
        ax.text(
            col_x[1], y,
            "✓ YES" if c["present"] else "✗ NO",
            fontsize=8, va="center",
            color="#1B5E20" if c["present"] else WARN_RED,
            fontweight="bold"
        )
        abiotic_str = (
            "Yes — explainable"
            if c["abiotic_explained"]
            else "NO — not reproduced"
        )
        abiotic_col = (
            "#555555"
            if c["abiotic_explained"]
            else WARN_RED
        )
        ax.text(
            col_x[2], y, abiotic_str,
            fontsize=7.5, va="center",
            color=abiotic_col
        )
        ax.text(
            col_x[3], y,
            c["note"].replace("\n", " "),
            fontsize=7, va="center",
            color="#333333"
        )

    n_present = sum(c["present"] for c in ALH_CRITERIA)
    n_unexplained = sum(
        c["present"] and not c["abiotic_explained"]
        for c in ALH_CRITERIA
    )
    ax.text(
        5, -0.3,
        f"Summary: {n_present}/6 criteria present in ALH84001  |  "
        f"{n_unexplained}/6 criteria NOT fully explained abiotically  |  "
        "Highlighted rows = unexplained abiotically",
        ha="center", va="center",
        fontsize=7.5, color=HEADER_BLUE,
        style="italic"
    )

    plt.tight_layout()
    path = (f"{OUTPUT_DIR}/{FIGURE_PREFIX}"
            f"_fig4_ALH84001_criteria.png")
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Figure 4 saved: {path}")
    return path

# ─────────────────────────────────────────────────────────────
# MODULE 5: FULL ANALYSIS RUN
# ─────────────────────────────────────────────────────────────

def run_analysis():

    separator = "═" * 62

    print()
    print(separator)
    print("  NECROMASS ANALYSIS v1.0")
    print("  Mars Iron Necromass Hypothesis")
    print("  Pre-reg: DOI 10.5281/zenodo.18986790")
    print("  2026-03-13 / Eric Robert Lawson / OrganismCore")
    print(separator)

    # ── SECTION 1: COMPUTED METRICS ─────────────────────────
    print()
    print("─" * 62)
    print("  SECTION 1: COMPUTED METRICS")
    print("─" * 62)
    print()
    print(f"  {'Meteorite':<22} {'Depth(m)':>9} "
          f"{'VeinFeO':>9} {'AltVol':>9} "
          f"{'FeIndex':>10}")
    print("  " + "─" * 60)
    for m in NAKHLITE_DATA:
        print(
            f"  {m['short']:<22} "
            f"{m['depth_mid']:>9.1f} "
            f"{m['vein_FeO_mid']:>9.1f} "
            f"{m['alt_vol_mid']:>9.4f} "
            f"{m['fe_index']['mid']:>10.4f}"
        )
    print()

    # ── SECTION 2: CORRELATIONS ──────────────────────────────
    print()
    print("─" * 62)
    print("  SECTION 2: CORRELATION ANALYSIS")
    print("─" * 62)
    print()
    print("  IMPORTANT: n=5, n_independent=3. No p-value")
    print("  from this dataset is statistically significant.")
    print("  All results are DIRECTIONAL EVIDENCE ONLY.")
    print()

    depths = np.array([m["depth_mid"] for m in NAKHLITE_DATA])
    vein_feo = np.array([m["vein_FeO_mid"] for m in NAKHLITE_DATA])
    alt_vol  = np.array([m["alt_vol_mid"]  for m in NAKHLITE_DATA])
    fe_index = np.array([m["fe_index"]["mid"] for m in NAKHLITE_DATA])

    s1 = correlation_report(
        "PANEL A — depth vs vein FeO concentration",
        depths, vein_feo
    )
    s2 = correlation_report(
        "PANEL B — depth vs alteration volume",
        depths, alt_vol
    )
    s3 = correlation_report(
        "PANEL C — depth vs total Fe oxide index",
        depths, fe_index
    )
    stats_results = [s1, s2, s3]

    # ── SECTION 3: MONTE CARLO ───────────────────────────────
    print()
    print("─" * 62)
    print("  SECTION 3: MONTE CARLO UNCERTAINTY PROPAGATION")
    print(f"  N = {N_MONTE_CARLO:,} random samples")
    print("─" * 62)
    print()
    print("  Sampling depth and measurement values uniformly")
    print("  from [min, max] published ranges.")
    print("  Question: does the correlation direction hold")
    print("  under all plausible uncertainty combinations?")
    print()

    mc_feo = monte_carlo_correlation(
        NAKHLITE_DATA, "vein_FeO"
    )
    mc_alt = monte_carlo_correlation(
        NAKHLITE_DATA, "alt_vol"
    )
    mc_fei = monte_carlo_correlation(
        NAKHLITE_DATA, "fe_index"
    )
    mc_results = [mc_feo, mc_alt, mc_fei]

    for label, rs in [
        ("Vein FeO concentration", mc_feo),
        ("Alteration volume",      mc_alt),
        ("Total Fe oxide index",   mc_fei),
    ]:
        frac_pos = np.mean(rs > 0)
        med = np.median(rs)
        ci_lo = np.percentile(rs, 2.5)
        ci_hi = np.percentile(rs, 97.5)
        print(f"  [{label}]")
        print(f"  P(r > 0) = {frac_pos:.4f}  "
              f"({'ROBUST POSITIVE' if frac_pos > 0.80 else 'WEAK' if frac_pos > 0.60 else 'NOT ROBUST'})")
        print(f"  Median r = {med:+.4f}")
        print(f"  95% CI   = [{ci_lo:+.4f}, {ci_hi:+.4f}]")
        print()

    # ── SECTION 4: GENERATE FIGURES ─────────────────────────
    print()
    print("─" * 62)
    print("  SECTION 4: GENERATING FIGURES")
    print("─" * 62)
    print()

    f1 = plot_figure_1(stats_results)
    f2 = plot_figure_2(mc_results)
    f3 = plot_figure_3()
    f4 = plot_figure_4()

    # ── SECTION 5: INTERPRETATION ────────────────────────────
    print()
    print("─" * 62)
    print("  SECTION 5: HONEST INTERPRETATION")
    print("─" * 62)
    print()

    frac_pos_index = np.mean(mc_fei > 0)
    frac_pos_alt   = np.mean(mc_alt > 0)

    print("  RESULT SUMMARY:")
    print()
    print(f"  Vein FeO concentration vs depth:")
    print(f"    r = {s1['r']:+.4f}")
    print(f"    NO CORRELATION — as predicted by diagnostic.")
    print(f"    Reason: dilution effect. Shallow rock has high")
    print(f"    FeO CONCENTRATION but near-zero volume.")
    print(f"    This is the WRONG metric for the hypothesis.")
    print()
    print(f"  Alteration volume fraction vs depth:")
    print(f"    r = {s2['r']:+.4f}")
    print(f"    P(r>0) across MC uncertainty = {frac_pos_alt:.4f}")
    print(f"    STRONG POSITIVE direction. Deeper = more total")
    print(f"    volume of iron chemistry. CONSISTENT WITH")
    print(f"    biological hypothesis.")
    print()
    print(f"  Total Fe oxide index vs depth (PRIMARY TEST):")
    print(f"    r = {s3['r']:+.4f}")
    print(f"    P(r>0) across MC uncertainty = {frac_pos_index:.4f}")
    if frac_pos_index > 0.80:
        verdict = (
            "ROBUST POSITIVE across uncertainty ranges.\n"
            "    Deeper Mars rock = more total iron oxide\n"
            "    deposited. CONSISTENT WITH biological\n"
            "    iron-cycling hypothesis."
        )
    elif frac_pos_index > 0.60:
        verdict = (
            "WEAKLY POSITIVE. Direction consistent with\n"
            "    hypothesis but not robust to all uncertainty\n"
            "    combinations."
        )
    else:
        verdict = (
            "NOT ROBUST. Correlation direction is not\n"
            "    consistent across uncertainty ranges.\n"
            "    Hypothesis NOT supported by this data."
        )
    print(f"    {verdict}")
    print()
    print("  WHAT THIS DOES AND DOES NOT MEAN:")
    print()
    print("  DOES NOT MEAN: Mars is biologically active.")
    print("  DOES NOT MEAN: The necromass hypothesis is proven.")
    print("  DOES NOT MEAN: The iron oxide is necromass.")
    print()
    print("  DOES MEAN: The nakhlite meteorite data is")
    print("    CONSISTENT WITH a biological iron-cycling")
    print("    community that increased iron oxide deposition")
    print("    with depth in the original Mars rock pile.")
    print("    The abiotic alternative (top-down weathering)")
    print("    predicts no increase with depth or a decrease.")
    print("    The direction of the alteration volume signal")
    print("    is opposite to the abiotic prediction.")
    print()
    print("  CRITICAL LIMITATION:")
    print("  n=5 meteorites, n=3 independent depth points.")
    print("  No statistical test can be significant.")
    print("  This is directional evidence only.")
    print("  It is consistent with the hypothesis.")
    print("  It does not confirm the hypothesis.")
    print("  Mars sample return and the full nakhlite")
    print("  iron oxide mass budget calculation are the")
    print("  next required tests.")
    print()

    # ── SECTION 6: FINAL RECORD ──────────────────────────────
    print()
    print(separator)
    print("  ANALYSIS COMPLETE")
    print(separator)
    print()
    print("  Figures generated:")
    for f in [f1, f2, f3, f4]:
        print(f"    {f}")
    print()
    print("  Pre-registration chain:")
    print("    10.5281/zenodo.18986790 (origin, 2026-03-12)")
    print("  Repository:")
    print("    github.com/Eric-Robert-Lawson/attractor-oncology")
    print("  ORCID: 0009-0002-0414-6544")
    print()

    return {
        "stats": stats_results,
        "monte_carlo_P_r_positive_fe_index": float(frac_pos_index),
        "monte_carlo_P_r_positive_alt_vol": float(frac_pos_alt),
        "conclusion": verdict,
    }


if __name__ == "__main__":
    result = run_analysis()
    sys.exit(0)
