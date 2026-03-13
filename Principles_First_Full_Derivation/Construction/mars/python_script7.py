#!/usr/bin/env python3
"""
NECROMASS SAM ANOMALY — SCRIPT 7 v1.0
=======================================
Document ID:  NECROMASS_SAM_v1.0
Date:         2026-03-13
Author:       Eric Robert Lawson / OrganismCore
Pre-reg DOI:  10.5281/zenodo.18986790

PURPOSE
-------
The Curiosity SAM instrument measured the most
isotopically depleted carbon ever recorded on Mars
in-situ. This is the anchor point for Script 7.

  Observed: δ¹³C = -137‰ (most depleted)
  Range:    -137 to +22‰ across 24 samples
  Source:   House et al. 2022, PNAS 119:e2115651119
  Location: Gale Crater sedimentary rocks
            (Murray formation)
  Method:   Pyrolysis-GC-MS in SAM instrument

This is 180‰ lighter than Martian atmospheric CO₂
(+43‰). That depletion requires explanation.

THE THREE HOUSE 2022 EXPLANATIONS
-----------------------------------
House et al. 2022 proposed three explanations,
none eliminated:

  1. UV photolysis of biologically produced methane.
     Subsurface microbes produce CH₄ depleted ~60-90‰
     from DIC. CH₄ rises to surface, UV photolysis
     creates even-more-depleted organics.

  2. UV photoreduction of atmospheric CO₂.
     Atmospheric CO₂ photolytically reduced to CO,
     which is strongly ¹³C-depleted due to isotope
     selective absorption. CO polymerises to organics.

  3. Deposition of cosmic dust from galactic
     molecular cloud passage.
     Interstellar chemistry produces extremely depleted
     carbon compounds delivered as cosmic dust.

POST-HOUSE 2022 DEVELOPMENT — CRITICAL
----------------------------------------
Ueno et al. 2024 (Nature Geoscience,
DOI: 10.1038/s41561-024-01443-z):

  Experimental + photochemical model showing
  UV photolysis of CO₂ in early Mars reducing
  atmosphere (CO-rich) produces CO depleted
  to δ¹³C ≈ -135‰ to -170‰.

  Subsequent organic synthesis from this CO
  produces organics matching the SAM observations
  WITHOUT biological activity.

  This is the strongest abiotic competitor
  published to date for the -137‰ signal.
  It was published AFTER House 2022.
  It MUST be reported honestly.

THE BIOLOGICAL FRACTIONATION CHAIN
------------------------------------
This script tests whether the biological
explanation (Explanation 1) is quantitatively
reachable from the fractionation chain
established in Scripts 1-6.

The chain for biological explanation:

  Step 0: Mars atmospheric CO₂
            δ¹³C = +43‰
            Source: Niles et al. 2013

  Step 1: CO₂ dissolves to DIC in fracture water
            δ¹³C = +10‰ (Lafayette carbonate proxy)
            Fractionation: -33‰ (CO₂ → DIC → carbonate)
            Source: Scripts 4-5

  Step 2: FeOB/methanogen fixes DIC
            Fractionation from DIC: -25 to -90‰
            (methanogen WL, subsurface H₂-limited)
            Source: Okumura 2016; Hattori 2012
            → CH₄ δ¹³C: -15 to -80‰

  Step 3: CH₄ migrates to surface
            No fractionation during migration.
            δ¹³C unchanged: -15 to -80‰

  Step 4: UV photolysis of CH₄ at surface
            Produces organic residue depleted
            relative to source methane.
            House 2022: this step is the mechanism.
            Fractionation factor for UV photolysis
            of CH₄: -30 to -60‰ additional depletion
            → Final organic δ¹³C: -45 to -140‰

  Target: -137‰
  Question: Does the biological chain reach this?

THE ABIOTIC FRACTIONATION CHAIN
---------------------------------
  Step 0: Mars atmospheric CO₂
            δ¹³C = +43‰

  Step 1: CO₂ photolysis by UV
            Fractionation: -135 to -180‰
            Source: Ueno et al. 2024 Nat Geosci
            Yoshida et al. 2023 arxiv:2302.12457
            → CO δ¹³C: -92 to -137‰

  Step 2: CO polymerises to organic compounds
            Fractionation: near 0 to -10‰
            → Organic δ¹³C: -92 to -147‰

  Target: -137‰
  This chain reaches -137‰ without biology.

WHAT THIS SCRIPT COMPUTES
--------------------------
1. Monte Carlo fractionation chain for biological
   pathway. What range of final δ¹³C values can
   the biological chain produce?

2. Monte Carlo fractionation chain for abiotic
   photochemical pathway. What range?

3. Monte Carlo for galactic dust pathway.

4. Overlap test: What fraction of each pathway's
   predicted distribution reaches the observed
   SAM range (-137 to -70‰)?

5. HONEST VERDICT: Which pathway(s) can produce
   the SAM anomaly? Can biology be excluded?
   Can it be confirmed?

6. Connection to Lafayette signals:
   Are the Lafayette signals (Scripts 1-5)
   quantitatively on the same fractionation chain
   as the SAM anomaly, or are they independent?

ALL DATA SOURCES
----------------
  House et al. 2022, PNAS 119:e2115651119
    SAM observations: -137 to +22‰, 24 samples.

  Ueno et al. 2024, Nat Geosci
    DOI: 10.1038/s41561-024-01443-z
    UV photolysis produces -135‰ organics abiotically.

  Yoshida et al. 2023, arXiv:2302.12457
    CO from CO₂ photolysis: -170‰ minimum.

  Okumura et al. 2016, ResearchGate
    Methanogen ε (CO₂→CH₄): -65 to -90‰ (H₂-limited)

  Hattori et al. 2012, Geochemical Journal 46
    Methanogen ε: -52 to -74‰ (deep aquifer)

  Niles et al. 2013, SSR 174:301
    Mars atmospheric CO₂: +43‰

  Scripts 1-6 (this series)
    Lafayette DIC proxy: +10.2‰
    Lafayette organic: -31.5‰
    Lafayette fractionation: -41.3‰

DEPENDENCIES
------------
  numpy, scipy, matplotlib
"""

import sys
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

np.random.seed(42)
N_MC = 200_000

# ─────────────────────────────────────────────────────────────
# MODULE 1: PUBLISHED PARAMETERS
# ─────────────────────────────────────────────────────────────

# SAM observations
SAM = {
    "d13C_min":    -137.0,   # most depleted
    "d13C_max":    +22.0,    # most enriched
    "d13C_depleted_threshold": -70.0,  # House 2022 cutoff
    "n_samples":   24,
    "n_depleted":  10,       # <-70‰
    "source":      "House et al. 2022, PNAS 119:e2115651119",
}

# Mars atmospheric CO���
MARS_CO2 = {
    "d13C_mid":  +43.0,
    "d13C_lo":   +38.0,
    "d13C_hi":   +46.0,
    "source":    "Niles et al. 2013, SSR 174:301",
}

# Lafayette DIC proxy (from Script 4)
LAFAYETTE_DIC = {
    "d13C_mid":  +10.2,
    "d13C_lo":    +7.3,
    "d13C_hi":   +13.1,
    "source":    "Bridges & Grady 2000, Scripts 4-5",
}

# Lafayette organic carbon (from Scripts 3-5)
LAFAYETTE_ORGANIC = {
    "d13C_mid": -31.5,
    "d13C_lo":  -36.6,
    "d13C_hi":  -27.6,
    "source":   "Steele et al. 2012; Scripts 3-5",
}

# ─────────────────────────────────────────────────────────────
# MODULE 2: FRACTIONATION STEP PARAMETERS
# All from published sources
# ─────────────────────────────────────────────────────────────

STEPS = {

    # ── BIOLOGICAL PATHWAY ───────────────────────────────────

    "CO2_to_DIC": {
        "description": "CO₂ dissolution + water-rock → DIC",
        "frac_lo":  -35.0,
        "frac_mid": -33.0,
        "frac_hi":  -28.0,
        "source":   "Scripts 4-5, carbonate proxy",
        "notes":    (
            "CO₂(atm,+43‰) → DIC(+10‰). "
            "Observed depletion = 33‰."
        ),
    },

    "DIC_to_CH4_methanogen_subsurface": {
        "description": "Methanogen DIC→CH₄ (H₂-limited)",
        "frac_lo":  -90.0,
        "frac_mid": -70.0,
        "frac_hi":  -52.0,
        "source":   (
            "Okumura et al. 2016 (H₂-limited: -65 to -90‰); "
            "Hattori et al. 2012 (-52 to -74‰ deep aquifer)"
        ),
        "notes":    (
            "H₂-limited subsurface conditions "
            "produce maximum fractionation."
        ),
    },

    "DIC_to_biomass_FeOB": {
        "description": "FeOB WL carbon fixation (DIC→biomass)",
        "frac_lo":  -41.0,
        "frac_mid": -25.0,
        "frac_hi":  -15.0,
        "source":   (
            "Scripts 4-5; House 2003 GCA; "
            "Chivian 2008. Lafayette obs: -41.3‰"
        ),
        "notes":    (
            "Lafayette measured fractionation "
            "used as upper bound (-41‰). "
            "Standard WL -15 to -36‰."
        ),
    },

    "CH4_UV_photolysis_depletion": {
        "description": "UV photolysis of CH₄ → organic residue",
        "frac_lo":  -60.0,
        "frac_mid": -40.0,
        "frac_hi":  -20.0,
        "source":   (
            "House et al. 2022 PNAS (mechanism); "
            "Estimated from photolysis fractionation "
            "literature. Exact value not fully constrained. "
            "This is the most uncertain step."
        ),
        "notes":    (
            "UV photolysis of CH₄ at Mars surface "
            "produces organic residue depleted "
            "relative to source methane. "
            "Range estimated from House 2022 "
            "supplementary discussion."
        ),
    },

    # ── ABIOTIC PHOTOCHEMICAL PATHWAY ────────────────────────

    "CO2_UV_photolysis_to_CO": {
        "description": "CO₂ UV photolysis → CO (abiotic)",
        "frac_lo":  -180.0,
        "frac_mid": -150.0,
        "frac_hi":  -105.0,
        "source":   (
            "Yoshida et al. 2023 arXiv:2302.12457 "
            "(model minimum -170‰); "
            "Ueno et al. 2024 Nat Geosci "
            "(experimental -135‰ to -170‰)"
        ),
        "notes":    (
            "Preferential UV absorption by ¹²CO₂ "
            "produces strongly ¹³C-depleted CO. "
            "Early Mars reducing (CO-rich) atmosphere "
            "amplifies this fractionation."
        ),
    },

    "CO_to_organic": {
        "description": "CO polymerisation → organic compounds",
        "frac_lo":  -10.0,
        "frac_mid":  -3.0,
        "frac_hi":    0.0,
        "source":   "Ueno et al. 2024 Nat Geosci",
        "notes":    (
            "Organic synthesis from CO introduces "
            "small additional fractionation. "
            "Ueno 2024 measured near-zero to -10‰."
        ),
    },

    # ── GALACTIC DUST PATHWAY ────────────────────────────────

    "galactic_dust_delivery": {
        "description": "Galactic molecular cloud dust",
        "frac_lo":  -200.0,
        "frac_mid": -150.0,
        "frac_hi":  -90.0,
        "source":   (
            "House et al. 2022 (hypothesis); "
            "ISM carbon isotope chemistry "
            "can produce extreme depletion."
        ),
        "notes":    (
            "Not well constrained. "
            "Ion-molecule reactions in cold dense "
            "molecular clouds can fractionate "
            "carbon to extreme values. "
            "Flux to Mars surface uncertain."
        ),
    },
}

# ─────────────────────────────────────────────────────────────
# MODULE 3: FRACTIONATION CHAIN MONTE CARLO
# ─────────────────────────────────────────────────────────────

def sample_step(step_key, n=1):
    """Sample fractionation for a given step."""
    s = STEPS[step_key]
    return np.random.uniform(s["frac_lo"], s["frac_hi"], n)

def run_biological_chain(n=N_MC):
    """
    Full biological fractionation chain:
    Mars CO₂ → DIC → CH₄ (methanogen) → UV photolysis → organic

    Returns array of final organic δ¹³C values.
    Also returns intermediate values.
    """
    # Step 0: Mars atmospheric CO₂
    d_CO2 = np.random.uniform(
        MARS_CO2["d13C_lo"], MARS_CO2["d13C_hi"], n
    )

    # Step 1: CO₂ → DIC (dissolution + water-rock)
    frac_DIC = sample_step("CO2_to_DIC", n)
    d_DIC = d_CO2 + frac_DIC

    # Step 2: DIC → CH₄ (methanogen, H₂-limited)
    frac_CH4 = sample_step(
        "DIC_to_CH4_methanogen_subsurface", n
    )
    d_CH4 = d_DIC + frac_CH4

    # Step 3: CH₄ → UV photolysis → organic residue
    frac_UV = sample_step("CH4_UV_photolysis_depletion", n)
    d_organic_final = d_CH4 + frac_UV

    return {
        "d_CO2":          d_CO2,
        "d_DIC":          d_DIC,
        "d_CH4":          d_CH4,
        "d_organic_final":d_organic_final,
    }

def run_biological_FeOB_chain(n=N_MC):
    """
    FeOB biological chain (no methane step):
    Mars CO₂ → DIC → FeOB biomass (WL fixation)

    This is the Lafayette chain — no UV step.
    Returns final organic δ¹³C (biomass/necromass).
    """
    d_CO2 = np.random.uniform(
        MARS_CO2["d13C_lo"], MARS_CO2["d13C_hi"], n
    )
    frac_DIC = sample_step("CO2_to_DIC", n)
    d_DIC = d_CO2 + frac_DIC
    frac_fix = sample_step("DIC_to_biomass_FeOB", n)
    d_biomass = d_DIC + frac_fix
    return {
        "d_CO2":    d_CO2,
        "d_DIC":    d_DIC,
        "d_biomass":d_biomass,
    }

def run_abiotic_photochemical_chain(n=N_MC):
    """
    Abiotic photochemical chain:
    Mars CO₂ → UV photolysis → CO → organic compounds

    Source: Ueno et al. 2024; Yoshida et al. 2023
    """
    d_CO2 = np.random.uniform(
        MARS_CO2["d13C_lo"], MARS_CO2["d13C_hi"], n
    )
    frac_CO = sample_step("CO2_UV_photolysis_to_CO", n)
    d_CO = d_CO2 + frac_CO
    frac_org = sample_step("CO_to_organic", n)
    d_organic = d_CO + frac_org
    return {
        "d_CO2":    d_CO2,
        "d_CO":     d_CO,
        "d_organic":d_organic,
    }

def run_galactic_dust_chain(n=N_MC):
    """
    Galactic molecular cloud dust delivery.
    Most uncertain pathway.
    """
    d_CO2 = np.full(n, MARS_CO2["d13C_mid"])
    frac_gal = sample_step("galactic_dust_delivery", n)
    # Galactic dust δ¹³C: relative to solar system baseline
    d_organic = MARS_CO2["d13C_mid"] + frac_gal
    return {
        "d_CO2":    d_CO2,
        "d_organic":d_organic,
    }

# ─────────────────────────────────────────────────────────────
# MODULE 4: OVERLAP AND CONSISTENCY TESTS
# ─────────────────────────────────────────────────────────────

def P_reaches_target(d_array, target_lo, target_hi):
    """Fraction of MC samples reaching target range."""
    return float(np.mean(
        (d_array >= target_lo) & (d_array <= target_hi)
    ))

def P_more_depleted_than(d_array, threshold):
    """Fraction more depleted than threshold."""
    return float(np.mean(d_array <= threshold))

# ─────────────────────────────────────────────────────────────
# MODULE 5: LAFAYETTE CONNECTION TEST
# ─────────────────────────────────────────────────────────────

def lafayette_connection_test(bio_FeOB, bio_meth):
    """
    Test whether Lafayette signals (Scripts 1-5) sit
    on the same fractionation chain as the SAM anomaly.

    If the FeOB biomass prediction (Lafayette chain)
    AND the methanogen + UV chain (SAM chain)
    are both consistent with their respective observations:
    → Same biological system, different expressions.

    If Lafayette FeOB biomass matches Lafayette obs
    but methanogen + UV chain does not reach SAM:
    → Lafayette and SAM are independent signals.
    """
    # Lafayette FeOB biomass vs observed
    P_FeOB_matches_Lafayette = P_reaches_target(
        bio_FeOB["d_biomass"],
        LAFAYETTE_ORGANIC["d13C_lo"],
        LAFAYETTE_ORGANIC["d13C_hi"]
    )

    # Lafayette fractionation observed -41.3‰
    # Is this reachable from the FeOB chain?
    lafayette_DIC_obs = np.random.uniform(
        LAFAYETTE_DIC["d13C_lo"],
        LAFAYETTE_DIC["d13C_hi"],
        N_MC
    )
    frac_FeOB = sample_step("DIC_to_biomass_FeOB", N_MC)
    predicted_Lafayette_frac = frac_FeOB
    P_FeOB_frac = P_reaches_target(
        predicted_Lafayette_frac, -47.1, -35.5
    )

    return P_FeOB_matches_Lafayette, P_FeOB_frac

# ─────────────────────────────────────────────────────────────
# MODULE 6: PLOTTING
# ─────────────────────────────────────────────────────────────

HEADER_BLUE = "#14388C"
WARN_RED    = "#8B0000"

def setup_style():
    plt.rcParams.update({
        "font.family":      "DejaVu Sans",
        "font.size":        8.5,
        "axes.titlesize":   9,
        "axes.labelsize":   8,
        "axes.titlecolor":  HEADER_BLUE,
        "axes.labelcolor":  HEADER_BLUE,
        "axes.edgecolor":   "#CCCCCC",
        "figure.facecolor": "white",
        "axes.facecolor":   "#FAFAFA",
        "grid.color":       "#DDDDDD",
        "grid.linewidth":   0.5,
    })

def plot_figure_20(bio_meth, bio_FeOB, abiotic, galactic):
    """
    Figure 20: Four-pathway fractionation chain.
    Final δ¹³C distributions for each pathway.
    SAM observed range overlaid.
    """
    setup_style()
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        "Script 7 — SAM -137‰ Anomaly: "
        "Three Explanations + Lafayette Chain\n"
        "Final δ¹³C distributions for each pathway "
        "vs. SAM observations\n"
        "House 2022 PNAS  |  "
        "Ueno 2024 Nat Geosci  |  "
        "Yoshida 2023 arXiv  |  Scripts 1-6\n"
        "Pre-reg DOI: 10.5281/zenodo.18986790  |  "
        "2026-03-13  |  Eric Robert Lawson / OrganismCore",
        fontsize=8.5, color=HEADER_BLUE, y=1.01
    )

    panels = [
        (axes[0, 0],
         bio_meth["d_organic_final"],
         "Biological: methanogen + UV photolysis",
         "#1565C0", "Bio-methane + UV"),
        (axes[0, 1],
         bio_FeOB["d_biomass"],
         "Biological: FeOB WL (Lafayette chain)",
         "#2E7D32", "FeOB WL (Lafayette)"),
        (axes[1, 0],
         abiotic["d_organic"],
         "Abiotic: CO₂ UV photolysis → CO → organic\n"
         "(Ueno 2024 / Yoshida 2023)",
         "#BF360C", "Abiotic photochemical"),
        (axes[1, 1],
         galactic["d_organic"],
         "Galactic dust (molecular cloud)",
         "#795548", "Galactic dust"),
    ]

    for ax, data, title, color, label in panels:
        ax.hist(
            data, bins=100,
            color=color, alpha=0.65,
            edgecolor="white", linewidth=0.3,
            zorder=4, label=label
        )
        # SAM depleted zone
        ax.axvspan(
            SAM["d13C_min"],
            SAM["d13C_depleted_threshold"],
            color="#E91E63", alpha=0.15,
            zorder=2,
            label=f"SAM depleted zone\n({SAM['d13C_min']} "
                  f"to {SAM['d13C_depleted_threshold']}‰)"
        )
        ax.axvline(
            SAM["d13C_min"],
            color="#E91E63", linewidth=2.0,
            linestyle="--", zorder=5,
            label=f"SAM most depleted ({SAM['d13C_min']}‰)"
        )
        ax.axvline(
            LAFAYETTE_ORGANIC["d13C_mid"],
            color="#FF9800", linewidth=1.5,
            linestyle=":", zorder=5, alpha=0.8,
            label=f"Lafayette organic "
                  f"({LAFAYETTE_ORGANIC['d13C_mid']}‰)"
        )

        P_sam = P_reaches_target(
            data, SAM["d13C_min"],
            SAM["d13C_depleted_threshold"]
        )
        P_lt = P_more_depleted_than(
            data, SAM["d13C_min"]
        )
        med = np.median(data)

        ax.text(
            0.97, 0.97,
            f"Median: {med:.1f}‰\n"
            f"P(in SAM depleted zone): {P_sam:.4f}\n"
            f"P(< -137‰): {P_lt:.4f}",
            transform=ax.transAxes,
            fontsize=7.5, color=HEADER_BLUE,
            va="top", ha="right",
            bbox=dict(
                boxstyle="round,pad=0.3",
                facecolor="white",
                edgecolor=HEADER_BLUE,
                alpha=0.92
            )
        )
        ax.set_xlabel("Final organic δ¹³C (‰)", fontsize=8)
        ax.set_ylabel("MC count", fontsize=8)
        ax.set_title(title, fontsize=8, color=HEADER_BLUE)
        ax.legend(fontsize=6.5, framealpha=0.85)
        ax.grid(True, alpha=0.12, zorder=0)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    path = "./necromass_sam_fig20_pathways.png"
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Figure 20 saved: {path}")
    return path

def plot_figure_21(bio_meth, abiotic):
    """
    Figure 21: Step-by-step fractionation chains.
    Biological (methanogen + UV) vs abiotic photochemical.
    Shows where each pathway sits at each step.
    """
    setup_style()
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(
        "Script 7 — Fractionation Chain Anatomy\n"
        "Biological (methane + UV) vs "
        "Abiotic (CO₂ photolysis)\n"
        "Cumulative δ¹³C depletion at each step",
        fontsize=8.5, color=HEADER_BLUE
    )

    # Panel A: Biological chain steps
    ax = axes[0]
    steps_bio = [
        ("Martian CO₂",
         bio_meth["d_CO2"]),
        ("DIC (fracture water)",
         bio_meth["d_DIC"]),
        ("CH₄ (methanogen)",
         bio_meth["d_CH4"]),
        ("Organic (UV photolysis)",
         bio_meth["d_organic_final"]),
    ]
    y_pos = list(range(len(steps_bio)))
    for i, (label, data) in enumerate(steps_bio):
        med = np.median(data)
        lo  = np.percentile(data, 5)
        hi  = np.percentile(data, 95)
        ax.errorbar(
            med, i,
            xerr=[[med - lo], [hi - med]],
            fmt="o", color="#1565C0",
            markersize=8, capsize=5,
            elinewidth=1.2, capthick=1.2,
            zorder=5
        )
        ax.text(
            med + 3, i,
            f"{med:.0f}‰",
            va="center", fontsize=8,
            color="#1565C0"
        )

    ax.axvline(
        SAM["d13C_min"], color="#E91E63",
        linewidth=1.5, linestyle="--",
        zorder=4,
        label=f"SAM target ({SAM['d13C_min']}‰)"
    )
    ax.set_yticks(y_pos)
    ax.set_yticklabels(
        [s[0] for s in steps_bio], fontsize=8.5
    )
    ax.set_xlabel("δ¹³C (‰)", fontsize=8.5)
    ax.set_title(
        "Panel A — Biological chain\n"
        "(methanogen + UV photolysis)\n"
        "Median ± 90% CI at each step",
        fontsize=8.5, color=HEADER_BLUE
    )
    ax.invert_yaxis()
    ax.legend(fontsize=7.5, framealpha=0.85)
    ax.grid(True, alpha=0.12, axis="x", zorder=0)

    # Panel B: Abiotic chain steps
    ax2 = axes[1]
    steps_ab = [
        ("Martian CO₂",
         abiotic["d_CO2"]),
        ("CO (UV photolysis of CO₂)",
         abiotic["d_CO"]),
        ("Organic compound",
         abiotic["d_organic"]),
    ]
    y_pos2 = list(range(len(steps_ab)))
    for i, (label, data) in enumerate(steps_ab):
        med = np.median(data)
        lo  = np.percentile(data, 5)
        hi  = np.percentile(data, 95)
        ax2.errorbar(
            med, i,
            xerr=[[med - lo], [hi - med]],
            fmt="s", color="#BF360C",
            markersize=8, capsize=5,
            elinewidth=1.2, capthick=1.2,
            zorder=5
        )
        ax2.text(
            med + 3, i,
            f"{med:.0f}‰",
            va="center", fontsize=8,
            color="#BF360C"
        )

    ax2.axvline(
        SAM["d13C_min"], color="#E91E63",
        linewidth=1.5, linestyle="--",
        zorder=4,
        label=f"SAM target ({SAM['d13C_min']}‰)"
    )
    ax2.set_yticks(y_pos2)
    ax2.set_yticklabels(
        [s[0] for s in steps_ab], fontsize=8.5
    )
    ax2.set_xlabel("δ¹³C (‰)", fontsize=8.5)
    ax2.set_title(
        "Panel B — Abiotic photochemical chain\n"
        "(CO₂ UV photolysis, Ueno 2024)\n"
        "Median ± 90% CI at each step",
        fontsize=8.5, color=HEADER_BLUE
    )
    ax2.invert_yaxis()
    ax2.legend(fontsize=7.5, framealpha=0.85)
    ax2.grid(True, alpha=0.12, axis="x", zorder=0)

    plt.tight_layout()
    path = "./necromass_sam_fig21_chains.png"
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Figure 21 saved: {path}")
    return path

def plot_figure_22(bio_meth, bio_FeOB, abiotic,
                   galactic, P_dict):
    """
    Figure 22: Summary comparison bar chart.
    P(reaches SAM depleted zone) for each pathway.
    The key discriminating figure.
    """
    setup_style()
    fig, ax = plt.subplots(figsize=(11, 6))

    pathways = [
        ("Biological\nmethane + UV",    "#1565C0",
         P_dict["bio_meth_SAM"]),
        ("Biological\nFeOB WL\n(Lafayette chain)",
         "#2E7D32", P_dict["FeOB_SAM"]),
        ("Abiotic\nCO₂ photolysis\n(Ueno 2024)",
         "#BF360C", P_dict["abiotic_SAM"]),
        ("Galactic\ndust",              "#795548",
         P_dict["galactic_SAM"]),
    ]

    labels = [p[0] for p in pathways]
    colors = [p[1] for p in pathways]
    P_vals = [p[2] for p in pathways]

    bars = ax.bar(
        range(len(pathways)), P_vals,
        color=colors, alpha=0.78,
        edgecolor="#333333", linewidth=0.7,
        zorder=4
    )
    ax.axhline(
        0.05, color="#888888",
        linewidth=1.0, linestyle="--",
        zorder=5, label="P = 0.05 threshold"
    )
    ax.axhline(
        0.50, color="#388E3C",
        linewidth=1.0, linestyle=":",
        zorder=5, label="P = 0.50"
    )

    for bar, P in zip(bars, P_vals):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            P + 0.01,
            f"P={P:.4f}",
            ha="center", va="bottom",
            fontsize=9.5, fontweight="bold",
            color="#111111"
        )

    ax.set_xticks(range(len(pathways)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel(
        "P(predicted δ¹³C reaches SAM depleted zone\n"
        f"[{SAM['d13C_min']} to "
        f"{SAM['d13C_depleted_threshold']}‰])",
        fontsize=8.5
    )
    ax.set_title(
        "Script 7 — Which Pathway Reaches the SAM "
        "-137‰ Anomaly?\n"
        "P(consistent) for each explanation "
        "from Monte Carlo chains\n"
        "House 2022 PNAS  |  "
        "Ueno 2024 Nat Geosci  |  Scripts 1-6",
        fontsize=8.5, color=HEADER_BLUE
    )
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=7.5, framealpha=0.88)
    ax.grid(True, alpha=0.12, axis="y", zorder=0)

    plt.tight_layout()
    path = "./necromass_sam_fig22_comparison.png"
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Figure 22 saved: {path}")
    return path

# ───────────────────────────────────────────────────────��─────
# MODULE 7: MAIN RUN
# ─────────────────────────────────────────────────────────────

def run():

    sep = "═" * 62

    print()
    print(sep)
    print("  NECROMASS SAM ANOMALY — SCRIPT 7 v1.0")
    print("  Mars Iron Necromass Hypothesis")
    print("  Pre-reg: DOI 10.5281/zenodo.18986790")
    print("  2026-03-13 / Eric Robert Lawson / OrganismCore")
    print(sep)

    # ── SECTION 1: DATA ──────────────────────────────────────
    print()
    print("─" * 62)
    print("  SECTION 1: SAM OBSERVATIONS AND PARAMETERS")
    print("─" * 62)
    print()
    print(f"  SAM most depleted:    {SAM['d13C_min']}‰")
    print(f"  SAM range:            "
          f"{SAM['d13C_min']} to {SAM['d13C_max']}‰")
    print(f"  SAM depleted samples: "
          f"{SAM['n_depleted']} / {SAM['n_samples']} "
          f"< {SAM['d13C_depleted_threshold']}‰")
    print(f"  Source:               {SAM['source']}")
    print()
    print(f"  Mars CO₂:             "
          f"{MARS_CO2['d13C_mid']}‰ "
          f"[{MARS_CO2['d13C_lo']}, {MARS_CO2['d13C_hi']}]")
    print(f"  Total depletion req:  "
          f"{MARS_CO2['d13C_mid'] - SAM['d13C_min']:.0f}‰")
    print()
    print("  POST-HOUSE 2022 DEVELOPMENT:")
    print("  Ueno et al. 2024, Nat Geosci:")
    print("    UV photolysis CO₂ → CO → organic")
    print("    produces -135‰ to -170‰ ABIOTICALLY.")
    print("    This matches SAM without biology.")
    print("    Published after House 2022.")
    print("    Must be reported honestly.")
    print()

    # ── SECTION 2: FRACTIONATION CHAINS ─────────────────────
    print()
    print("─" * 62)
    print(f"  SECTION 2: MONTE CARLO CHAINS  N={N_MC:,}")
    print("─" * 62)
    print()
    print("  Running biological (methanogen + UV)...")
    bio_meth  = run_biological_chain(N_MC)
    print("  Running biological (FeOB WL, Lafayette)...")
    bio_FeOB  = run_biological_FeOB_chain(N_MC)
    print("  Running abiotic photochemical (Ueno 2024)...")
    abiotic   = run_abiotic_photochemical_chain(N_MC)
    print("  Running galactic dust...")
    galactic  = run_galactic_dust_chain(N_MC)
    print("  Complete.")
    print()

    # Chain summary statistics
    chains = {
        "Biological methane+UV":
            bio_meth["d_organic_final"],
        "Biological FeOB WL (Lafayette)":
            bio_FeOB["d_biomass"],
        "Abiotic CO₂ photolysis (Ueno 2024)":
            abiotic["d_organic"],
        "Galactic dust":
            galactic["d_organic"],
    }
    print(f"  {'Pathway':<40} {'Median':>8} "
          f"{'5th pc':>8} {'95th pc':>8}")
    print("  " + "─" * 66)
    for name, data in chains.items():
        med = np.median(data)
        lo  = np.percentile(data, 5)
        hi  = np.percentile(data, 95)
        print(
            f"  {name:<40} "
            f"{med:>8.1f}‰ "
            f"{lo:>8.1f}‰ "
            f"{hi:>8.1f}‰"
        )
    print()

    # ── SECTION 3: SAM OVERLAP TESTS ─────────────────────────
    print()
    print("─" * 62)
    print("  SECTION 3: SAM DEPLETED ZONE OVERLAP TESTS")
    print(f"  Target: {SAM['d13C_min']} to "
          f"{SAM['d13C_depleted_threshold']}‰")
    print("─" * 62)
    print()

    P_bio_meth = P_reaches_target(
        bio_meth["d_organic_final"],
        SAM["d13C_min"],
        SAM["d13C_depleted_threshold"]
    )
    P_FeOB = P_reaches_target(
        bio_FeOB["d_biomass"],
        SAM["d13C_min"],
        SAM["d13C_depleted_threshold"]
    )
    P_abiotic = P_reaches_target(
        abiotic["d_organic"],
        SAM["d13C_min"],
        SAM["d13C_depleted_threshold"]
    )
    P_galactic = P_reaches_target(
        galactic["d_organic"],
        SAM["d13C_min"],
        SAM["d13C_depleted_threshold"]
    )

    P_dict = {
        "bio_meth_SAM":  P_bio_meth,
        "FeOB_SAM":      P_FeOB,
        "abiotic_SAM":   P_abiotic,
        "galactic_SAM":  P_galactic,
    }

    print(f"  {'Pathway':<40} {'P(in SAM zone)':>16}")
    print("  " + "─" * 58)
    for name, P in [
        ("Biological methane + UV", P_bio_meth),
        ("Biological FeOB WL (Lafayette)", P_FeOB),
        ("Abiotic CO₂ photolysis", P_abiotic),
        ("Galactic dust", P_galactic),
    ]:
        flag = (
            "*** REACHES TARGET ***"
            if P > 0.05 else ""
        )
        print(
            f"  {name:<40} {P:>16.4f}  {flag}"
        )
    print()

    # ── SECTION 4: LAFAYETTE CONNECTION ──────────────────────
    print()
    print("─" * 62)
    print("  SECTION 4: LAFAYETTE CONNECTION TEST")
    print("─" * 62)
    print()
    P_FeOB_Laf, P_FeOB_frac = lafayette_connection_test(
        bio_FeOB, bio_meth
    )
    print(f"  P(FeOB chain matches Lafayette organic δ¹³C):")
    print(f"    {P_FeOB_Laf:.4f}")
    print()
    print(f"  P(FeOB chain matches Lafayette fractionation):")
    print(f"    {P_FeOB_frac:.4f}")
    print()
    print(f"  Lafayette organic mid: "
          f"{LAFAYETTE_ORGANIC['d13C_mid']}‰")
    print(f"  FeOB chain median:     "
          f"{np.median(bio_FeOB['d_biomass']):.1f}‰")
    print()

    # ── SECTION 5: FIGURES ───────────────────────────────────
    print()
    print("─" * 62)
    print("  SECTION 5: GENERATING FIGURES")
    print("─" * 62)
    print()
    f20 = plot_figure_20(bio_meth, bio_FeOB,
                         abiotic, galactic)
    f21 = plot_figure_21(bio_meth, abiotic)
    f22 = plot_figure_22(bio_meth, bio_FeOB,
                         abiotic, galactic, P_dict)

    # ── SECTION 6: INTERPRETATION ────────────────────────────
    print()
    print("─" * 62)
    print("  SECTION 6: HONEST INTERPRETATION")
    print("─" * 62)
    print()

    print("  WHICH PATHWAYS REACH THE SAM ANOMALY?")
    print()

    if P_bio_meth > 0.05:
        print(f"  Biological methane + UV:  P = {P_bio_meth:.4f}")
        print("    REACHES the SAM depleted zone.")
        print("    The biological fractionation chain")
        print("    (methanogen + UV photolysis) can")
        print("    produce -137‰ organics.")
        print("    This is quantitatively consistent with")
        print("    the biological explanation in House 2022.")
    else:
        print(f"  Biological methane + UV:  P = {P_bio_meth:.4f}")
        print("    Does NOT consistently reach SAM zone.")

    print()
    if P_abiotic > 0.05:
        print(f"  Abiotic CO₂ photolysis:   P = {P_abiotic:.4f}")
        print("    REACHES the SAM depleted zone.")
        print("    Ueno et al. 2024 (Nat Geosci):")
        print("    UV photolysis of CO₂ in reducing")
        print("    atmosphere produces -135‰ organics.")
        print("    This matches SAM WITHOUT biology.")
        print("    PUBLISHED AFTER HOUSE 2022.")
        print("    This is the current strongest")
        print("    abiotic explanation for the SAM signal.")
    else:
        print(f"  Abiotic CO₂ photolysis:   P = {P_abiotic:.4f}")
        print("    Does NOT consistently reach SAM zone.")

    print()
    print(f"  Galactic dust:            P = {P_galactic:.4f}")
    gal_status = ("REACHES" if P_galactic > 0.05
                  else "uncertain — parameter range poorly constrained")
    print(f"    {gal_status}.")
    print("    Least constrained pathway.")
    print("    Cannot be evaluated quantitatively.")
    print()

    print("  CRITICAL COMPARISON:")
    print()
    print("  If P_bio_meth AND P_abiotic both > 0.05:")
    print("    CANNOT DISTINGUISH biology from abiotic")
    print("    using the fractionation chain alone.")
    print("    The -137‰ SAM anomaly is explained by")
    print("    BOTH pathways. Not discriminating.")
    print()
    print("  If only one reaches the target:")
    print("    That pathway is the better explanation.")
    print()

    both_reach = (P_bio_meth > 0.05 and P_abiotic > 0.05)
    neither_reach = (P_bio_meth <= 0.05 and
                     P_abiotic <= 0.05)

    if both_reach:
        print("  RESULT: BOTH PATHWAYS REACH SAM TARGET.")
        print("  The SAM -137‰ anomaly is NOT a")
        print("  discriminating biosignature on its own.")
        print("  Both biological methane + UV AND")
        print("  abiotic CO₂ photolysis (Ueno 2024)")
        print("  predict overlapping δ¹³C distributions")
        print("  that include -137‰.")
        print()
        print("  The correct interpretation of House 2022")
        print("  in light of Ueno 2024 is:")
        print("    The SAM anomaly is CONSISTENT with")
        print("    biology but is NOT uniquely biological.")
        print("    Abiotic photochemistry can also produce")
        print("    it. The anomaly cannot be used as")
        print("    standalone evidence for biology.")
    elif neither_reach:
        print("  RESULT: NEITHER PATHWAY CONSISTENTLY")
        print("  REACHES SAM TARGET.")
        print("  The fractionation parameters require")
        print("  revision. Galactic dust or unknown")
        print("  mechanism may be required.")
    else:
        winner = ("biological" if P_bio_meth > P_abiotic
                  else "abiotic photochemical")
        print(f"  RESULT: {winner.upper()} PATHWAY")
        print("  is the better fit for the SAM anomaly.")

    print()
    print("  LAFAYETTE CONNECTION:")
    print()
    print(f"    P(FeOB matches Lafayette organic): "
          f"{P_FeOB_Laf:.4f}")
    print(f"    P(FeOB matches Lafayette frac):    "
          f"{P_FeOB_frac:.4f}")
    print()
    if P_FeOB_Laf > 0.05:
        print("    Lafayette organic carbon IS on the")
        print("    FeOB WL fractionation chain.")
        print("    Lafayette and SAM are DIFFERENT")
        print("    points on the Mars carbon cycle.")
        print("    Lafayette = subsurface biological")
        print("      necromass (if biological).")
        print("    SAM = surface expression via CH₄")
        print("      migration + UV (if biological),")
        print("      or CO₂ photolysis (if abiotic).")
        print("    The two signals are connected by")
        print("    a common Mars carbon fractionation")
        print("    framework, not the same process.")
    else:
        print("    Lafayette organic carbon does NOT")
        print("    consistently fall on the FeOB chain")
        print("    from this starting point.")
        print("    The chain parameterisation needs")
        print("    revision for the FeOB pathway.")

    # ── SECTION 7: FINAL CUMULATIVE RECORD ──────────────────
    print()
    print("─" * 62)
    print("  SECTION 7: COMPLETE SCIENTIFIC RECORD")
    print("  ALL SEVEN SCRIPTS")
    print("─" * 62)
    print()
    print("  SCRIPTS 1-5: Five signals in Lafayette.")
    print("    All consistent, no single abiotic")
    print("    explanation for combination.")
    print("    Best template: biological iron-cycling")
    print("    (match 0.833). P(bio leaning) = 0.990.")
    print()
    print("  SCRIPT 6: Desulforudis analogue failed.")
    print("    H₂ radiolysis energy too small 2000×.")
    print("    Fe(II) energy feasible 3×10⁶× larger.")
    print("    Hypothesis refined: FeOB, not Desulforudis.")
    print()
    print("  SCRIPT 7: SAM -137‰ anomaly.")
    print(f"    P(bio methane+UV reaches SAM): "
          f"{P_bio_meth:.4f}")
    print(f"    P(abiotic photochem reaches SAM): "
          f"{P_abiotic:.4f}")
    print()
    print("  THE COMPLETE HONEST POSITION:")
    print()
    print("  Lafayette (Scripts 1-5):")
    print("    Five independent signals. All consistent")
    print("    with biological iron-cycling necromass.")
    print("    No known single abiotic process explains")
    print("    all five simultaneously.")
    print("    Molecular identity (Script 5): biological")
    print("    leaning at P = 0.990.")
    print()
    print("  SAM anomaly (Script 7):")
    if both_reach:
        print("    CANNOT discriminate biology from")
        print("    abiotic CO₂ photolysis on this signal.")
        print("    Ueno 2024 provides a complete abiotic")
        print("    explanation published after House 2022.")
        print("    SAM anomaly is NOT a uniquely")
        print("    biological signal.")
    print()
    print("  What is and is not claimed:")
    print()
    print("  NOT CLAIMED:")
    print("    Mars has life.")
    print("    The SAM anomaly is biological.")
    print("    Any single signal is proof.")
    print()
    print("  WHAT IS CLAIMED:")
    print("    The Lafayette signals (five independent,")
    print("    real Martian material, published sources)")
    print("    have no known single abiotic explanation")
    print("    for their combination.")
    print("    The biological iron-cycling template")
    print("    matches the molecular data at 0.833.")
    print("    The Fe(II) energy budget is feasible.")
    print("    The SAM anomaly is consistent with both")
    print("    biological and abiotic photochemical")
    print("    explanations and cannot discriminate.")
    print()
    print("  The pattern is real.")
    print("  The record is complete.")
    print("  Seven scripts. Twenty-two figures.")
    print("  All real Martian material.")
    print("  All published sources.")
    print("  Pre-registered before analysis.")
    print("  Timestamp: 2026-03-13.")
    print()

    print(sep)
    print("  SCRIPT 7 COMPLETE")
    print("  ANALYSIS SERIES COMPLETE")
    print(sep)
    print()
    print("  Figures generated:")
    for f in [f20, f21, f22]:
        print(f"    {f}")
    print()
    print("  Full figure set (all scripts):")
    print("    Figs 1-7:   Script 1-2 (depth, rate)")
    print("    Figs 8-10:  Script 3 (morphology, carbon)")
    print("    Figs 11-13: Script 4 (WL fractionation)")
    print("    Figs 14-16: Script 5 (molecular identity)")
    print("    Figs 17-19: Script 6 (Daud analogue)")
    print("    Figs 20-22: Script 7 (SAM anomaly)")
    print("    Total: 22 figures")
    print()
    print("  Pre-registration chain:")
    print("    10.5281/zenodo.18986790")
    print("  Repository:")
    print("    github.com/Eric-Robert-Lawson/"
          "attractor-oncology")
    print("  ORCID: 0009-0002-0414-6544")
    print()

    return {
        "P_biological_methane_UV_SAM":   float(P_bio_meth),
        "P_biological_FeOB_SAM":         float(P_FeOB),
        "P_abiotic_photochem_SAM":       float(P_abiotic),
        "P_galactic_dust_SAM":           float(P_galactic),
        "P_FeOB_matches_Lafayette":      float(P_FeOB_Laf),
        "P_FeOB_matches_Lafayette_frac": float(P_FeOB_frac),
        "both_bio_and_abiotic_reach_SAM":bool(both_reach),
        "SAM_anomaly_discriminating":    bool(not both_reach),
        "bio_methane_chain_median_d13C":
            float(np.median(bio_meth["d_organic_final"])),
        "abiotic_chain_median_d13C":
            float(np.median(abiotic["d_organic"])),
        "Lafayette_FeOB_chain_median_d13C":
            float(np.median(bio_FeOB["d_biomass"])),
    }


if __name__ == "__main__":
    result = run()
    print("  MACHINE-READABLE SUMMARY:")
    print()
    for k, v in result.items():
        if isinstance(v, float):
            print(f"  {k}: {v:.6f}")
        else:
            print(f"  {k}: {v}")
    print()
    sys.exit(0)
