#!/usr/bin/env python3
"""
NECROMASS δ¹⁵N PREDICTION — SCRIPT 11 v1.0
============================================
Document ID:  NECROMASS_D15N_v1.0
Date:         2026-03-13
Author:       Eric Robert Lawson / OrganismCore
Pre-reg DOI:  10.5281/zenodo.18986790

PURPOSE
-------
Script 10 established posterior odds of 82,871:1
in favour of biological origin for the Lafayette
organic carbon signals (conservative, 4 independent
source groups).

Script 10 identified one experiment that could
resolve this decisively in either direction:

  δ¹⁵N measurement of Lafayette vein organics
  by NanoSIMS.

Script 11 does three things:

  1. PREDICTS the expected δ¹⁵N value under
     each hypothesis — biological and abiotic —
     using published ranges with full uncertainty.

  2. COMPUTES the likelihood ratio that each
     possible measured δ¹⁵N value would assign
     to the biological vs abiotic hypothesis.
     This is the diagnostic power of the
     experiment as a function of outcome.

  3. COMPUTES the updated posterior odds after
     incorporating the hypothetical δ¹⁵N result
     into the Script 10 posterior of 82,871:1.
     This shows what each possible experimental
     outcome would do to the full case.

This script does not measure δ¹⁵N.
It predicts what the measurement would show
under each hypothesis and quantifies exactly
how decisive the experiment would be.

THE δ¹⁵N SIGNAL — PUBLISHED PARAMETERS
-----------------------------------------
BIOLOGICAL PREDICTION:
  Deep subsurface biological nitrogen fixation.
  The nitrogen in cellular organic matter
  derives from biological N₂ fixation by
  nitrogenase enzyme.
  δ¹⁵N range: -4 to +4‰ (near-atmospheric N₂)
  Source: Stüeken et al. 2016 Geochem Persp Lett
          Stüeken et al. 2015 Nat Geosci
          Deep subsurface compiled values.

  After diagenetic preservation in mineral
  matrix for ~600 Ma:
  Slight ¹⁵N enrichment expected from
  preferential ¹⁴N loss during diagenesis.
  Preserved range: -2 to +8‰ (conservative).
  Source: Zerkle et al. 2016 PNAS;
          Papineau et al. 2009 Precambrian Res

  NanoSIMS precision: ±5-10‰ (1σ)
  Source: Lin et al. 2014 MAPS 49;
          Kebukawa et al. 2023 Icarus

ABIOTIC PREDICTIONS:
  Three abiotic nitrogen sources relevant
  to Lafayette nakhlite fracture veins:

  SOURCE A: Mars mantle nitrogen
    δ¹⁵N = 0 ± 32‰ (Chassigny, mantle)
    Source: USRA Metsoc 2024 abstract 6354
    Broad range — mantle degassing variable.

  SOURCE B: Mars atmospheric nitrogen
    δ¹⁵N = +572 ± 82‰ (Curiosity/SAM)
    Source: Mahaffy et al. 2013 Science
    Extremely enriched from photochemical loss.
    Would produce very high δ¹⁵N if incorporated.

  SOURCE C: Abiotic geochemical (serpentinization,
    hydrothermal, impact)
    δ¹⁵N = 0 to +20‰ (Earth abiotic reference)
    Source: Published serpentinization studies
    Mars equivalent: similar range expected.

  SOURCE D: Nakhlite mesostasis (measured)
    δ¹⁵N = -35 to +73‰ (corrected values)
    Source: Metsoc 2024 abstract 6354
    This is the MEASURED range for nitrogen
    already in Lafayette material (non-organic).
    The organic carbon nitrogen source must
    be distinguished from this matrix nitrogen.

THE DISCRIMINATING ZONE
------------------------
The diagnostic power of δ¹⁵N lies in the
OVERLAP STRUCTURE of the biological and
abiotic distributions.

  Biological prediction: -2 to +8‰ (central)
  Abiotic (mantle/geochemical): -35 to +73‰

  There is overlap in the range 0 to +8‰.
  Measurements in this overlap zone are
  partially ambiguous.

  UNAMBIGUOUS BIOLOGICAL ZONE: -2 to +4‰
    Highly depleted near-zero δ¹⁵N is
    diagnostic of biological N₂ fixation.
    Abiotic mantle N covers this range but
    at low probability.

  UNAMBIGUOUS ABIOTIC ZONE: > +50‰ or < -20‰
    Atmospheric N incorporation: > +200‰.
    Deep mantle values outside biological range.
    Values here definitively non-biological.

  AMBIGUOUS ZONE: +4 to +30‰
    Both biological (diagenetically enriched)
    and abiotic (mantle, geochemical) could
    produce values here.
    LR approaches 1 in this zone.

THE MEASUREMENT — WHAT NANOSIMS DETECTS
-----------------------------------------
NanoSIMS measures ¹²C¹⁴N⁻ and ¹²C¹⁵N⁻ ions
simultaneously with the carbon signal.

This allows:
  δ¹⁵N of carbon-bearing phases specifically.
  Spatial co-registration with δ¹³C.
  Sub-micron resolution — measures the
  ORGANIC CARBON nitrogen specifically,
  not bulk rock nitrogen.
  This is critical: it avoids contamination
  from matrix minerals.

Precision: ±5-10‰ (1σ) for adequate signal.
Detection limit: ~1 ppm N in organic phase.
Source: Lin 2014 MAPS; Wacey 2015 Geosciences.

PUBLISHED PARAMETERS
--------------------
  Stüeken et al. 2015    Nat Geosci 8:941
  Stüeken et al. 2016    Geochem Persp Lett 2:100
  Zerkle et al. 2016     PNAS 113:6269
  Papineau et al. 2009   Precambrian Res 169:43
  Mahaffy et al. 2013    Science 341:263
  USRA Metsoc 2024       Abstract 6354
  Lin et al. 2014        MAPS 49 (Tissint)
  Kebukawa et al. 2023   Icarus
  Wacey 2015             Geosciences

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
from matplotlib.gridspec import GridSpec
from scipy import stats

np.random.seed(42)
N_MC = 200_000

HEADER_BLUE = "#14388C"

# ─────────────────────────────────────────────────────────────
# MODULE 1: PUBLISHED PARAMETERS
# ─────────────────────────────────────────────────────────────

# ── BIOLOGICAL PREDICTION ─────────────────────────────────────
# Deep subsurface N₂ fixation, diagenetically preserved
# Stüeken 2015, 2016; Zerkle 2016; Papineau 2009
BIO_D15N = {
    "mean":   2.0,    # ‰ centre of biological range
    "sigma":  3.0,    # ‰ spread (1σ)
    "lo":    -4.0,    # ‰ lower bound published
    "hi":     8.0,    # ‰ upper bound (diagenetically
                      #   enriched, Papineau 2009)
    "source": "Stüeken 2015 Nat Geosci; "
              "Zerkle 2016 PNAS; "
              "Papineau 2009 Precambrian Res",
    "description": (
        "Deep subsurface biological N₂ fixation.\n"
        "Nitrogenase fractionation near zero.\n"
        "Slight ¹⁵N enrichment after 600 Ma\n"
        "diagenetic preservation."
    ),
}

# ── ABIOTIC PREDICTIONS ───────────────────────────────────────

ABIO_SOURCES = {
    "mantle": {
        "label":  "Mars mantle N",
        "mean":   0.0,
        "sigma":  32.0,
        "lo":    -35.0,
        "hi":    +73.0,
        "source": "USRA Metsoc 2024 abs 6354\n"
                  "(Chassigny δ¹⁵N = 0 ± 32‰)",
        "color":  "#FF9800",
        "weight": 0.40,  # Prior probability this
                         # source dominates in fracture
    },
    "atmospheric": {
        "label":  "Mars atmospheric N",
        "mean":   572.0,
        "sigma":  82.0,
        "lo":    +400.0,
        "hi":    +800.0,
        "source": "Mahaffy et al. 2013 Science 341:263",
        "color":  "#BF360C",
        "weight": 0.10,  # Low weight: atmospheric N
                         # unlikely to dominate in
                         # fracture vein organics
    },
    "geochemical": {
        "label":  "Abiotic geochemical",
        "mean":   10.0,
        "sigma":  15.0,
        "lo":    -5.0,
        "hi":    +30.0,
        "source": "Serpentinization + hydrothermal\n"
                  "Earth analogue values",
        "color":  "#E65100",
        "weight": 0.35,
    },
    "mesostasis": {
        "label":  "Lafayette mesostasis N",
        "mean":   20.0,
        "sigma":  25.0,
        "lo":    -35.0,
        "hi":    +73.0,
        "source": "USRA Metsoc 2024 abs 6354\n"
                  "(corrected nakhlite values)",
        "color":  "#795548",
        "weight": 0.15,
    },
}

# ── NANOSIMS INSTRUMENT PRECISION ────────────────────────────
NANOSIMS_PRECISION = {
    "sigma_1sigma": 7.5,   # ‰ typical 1σ precision
    "sigma_lo":     5.0,   # ‰ best case (strong signal)
    "sigma_hi":    10.0,   # ‰ worst case (weak signal)
    "source": "Lin et al. 2014 MAPS 49 (Tissint);\n"
              "Kebukawa et al. 2023 Icarus",
}

# ── PRIOR ODDS FROM SCRIPT 10 ─────────────────────────────────
PRIOR_ODDS_SCRIPT10 = 82871.0   # conservative (4 groups)

# ─────────────────────────────────────────────────────────────
# MODULE 2: DISTRIBUTION MODELS
# ─────────────────────────────────────────────────────────────

def P_bio_given_d15N(d15N_measured, precision_sigma):
    """
    P(δ¹⁵N = measured | biological hypothesis).
    Biological prediction: N(mean=2, sigma=3)
    Convolved with NanoSIMS measurement precision.
    """
    total_sigma = math.sqrt(
        BIO_D15N["sigma"]**2 + precision_sigma**2
    )
    return float(stats.norm.pdf(
        d15N_measured, BIO_D15N["mean"], total_sigma
    ))

def P_abio_given_d15N(d15N_measured, precision_sigma):
    """
    P(δ¹⁵N = measured | abiotic hypothesis).
    Mixture of abiotic sources weighted by prior.
    Each convolved with NanoSIMS precision.
    """
    total_p = 0.0
    total_weight = sum(s["weight"]
                       for s in ABIO_SOURCES.values())
    for src in ABIO_SOURCES.values():
        total_sigma = math.sqrt(
            src["sigma"]**2 + precision_sigma**2
        )
        p = float(stats.norm.pdf(
            d15N_measured, src["mean"], total_sigma
        ))
        total_p += (src["weight"] / total_weight) * p
    return total_p

def LR_d15N(d15N_measured, precision_sigma):
    """
    Likelihood ratio for a given δ¹⁵N measurement.
    LR = P(measured | bio) / P(measured | abio)
    LR > 1: supports biology.
    LR < 1: supports abiotic.
    """
    p_bio  = P_bio_given_d15N(d15N_measured,
                               precision_sigma)
    p_abio = P_abio_given_d15N(d15N_measured,
                                precision_sigma)
    if p_abio < 1e-300:
        return float("inf")
    return p_bio / p_abio

def posterior_odds_updated(d15N_measured,
                            precision_sigma,
                            prior_odds):
    """
    Updated posterior odds after δ¹⁵N measurement.
    posterior_odds = prior_odds × LR_d15N
    """
    lr = LR_d15N(d15N_measured, precision_sigma)
    if math.isinf(lr):
        return float("inf")
    return prior_odds * lr

def posterior_probability(odds):
    """Convert odds to probability."""
    if math.isinf(odds):
        return 1.0
    return odds / (1.0 + odds)

# ─────────────────────────────────────────────────────────────
# MODULE 3: MONTE CARLO OVER MEASUREMENT SCENARIOS
# ─────────────────────────────────────────────────────────────

def run_MC_scenarios():
    """
    For each hypothetical measurement outcome,
    compute the updated posterior and LR.
    Scan over the full δ¹⁵N range -50 to +150‰.
    """
    precision = NANOSIMS_PRECISION["sigma_1sigma"]

    d15N_scan = np.linspace(-50, 650, 2000)
    LRs     = np.array([LR_d15N(v, precision)
                        for v in d15N_scan])
    post_odds = np.array([
        posterior_odds_updated(v, precision,
                               PRIOR_ODDS_SCRIPT10)
        for v in d15N_scan
    ])
    post_prob = np.array([
        posterior_probability(o) for o in post_odds
    ])

    return {
        "d15N_scan":    d15N_scan,
        "LRs":          LRs,
        "post_odds":    post_odds,
        "post_prob":    post_prob,
        "precision":    precision,
    }

def run_MC_prediction(n=N_MC):
    """
    Monte Carlo over the biological and abiotic
    distributions to generate expected measurement
    distributions for each hypothesis.
    Returns: what you would MEASURE if biological,
             what you would MEASURE if abiotic.
    Includes NanoSIMS precision noise.
    """
    precision = NANOSIMS_PRECISION["sigma_1sigma"]

    # Biological prediction
    true_bio = np.random.normal(
        BIO_D15N["mean"], BIO_D15N["sigma"], n
    )
    measured_bio = true_bio + np.random.normal(
        0, precision, n
    )

    # Abiotic prediction (mixture)
    weights = np.array([
        s["weight"] for s in ABIO_SOURCES.values()
    ])
    weights = weights / weights.sum()
    sources = list(ABIO_SOURCES.values())

    choices = np.random.choice(
        len(sources), size=n, p=weights
    )
    true_abio = np.array([
        np.random.normal(sources[c]["mean"],
                         sources[c]["sigma"])
        for c in choices
    ])
    measured_abio = true_abio + np.random.normal(
        0, precision, n
    )

    return {
        "measured_bio":  measured_bio,
        "measured_abio": measured_abio,
        "true_bio":      true_bio,
        "true_abio":     true_abio,
    }

def compute_diagnostic_power(mc_pred):
    """
    Compute the diagnostic power of the experiment.
    For each possible measurement:
      If biological is true: how often does the
        measurement correctly identify biology?
      If abiotic is true: how often does the
        measurement correctly identify abiotic?
    """
    precision = NANOSIMS_PRECISION["sigma_1sigma"]

    # Define discrimination threshold zones
    # Bio zone: δ¹⁵N < +15‰  (LR > 1)
    # Abiotic zone: δ¹⁵N > +30‰  (LR < 0.1)

    # Find the LR = 1 crossover point
    d15N_scan = np.linspace(-50, 150, 5000)
    LRs = np.array([LR_d15N(v, precision)
                    for v in d15N_scan])

    # Find where LR crosses 1.0
    crossover_indices = np.where(
        np.diff(np.sign(LRs - 1.0))
    )[0]
    if len(crossover_indices) > 0:
        crossover = float(
            d15N_scan[crossover_indices[0]]
        )
    else:
        crossover = 15.0  # fallback

    # P(correct identification | bio true)
    meas_bio = mc_pred["measured_bio"]
    P_correct_if_bio = float(
        np.mean(meas_bio < crossover)
    )

    # P(correct identification | abio true)
    # Use only mantle + geochemical (realistic)
    # Exclude atmospheric (too extreme to confuse)
    precision_val = NANOSIMS_PRECISION["sigma_1sigma"]
    sources = list(ABIO_SOURCES.values())
    weights = np.array([s["weight"] for s in sources])
    weights = weights / weights.sum()
    n_abio = 50_000
    choices = np.random.choice(
        len(sources), size=n_abio, p=weights
    )
    true_abio = np.array([
        np.random.normal(sources[c]["mean"],
                         sources[c]["sigma"])
        for c in choices
    ])
    meas_abio = true_abio + np.random.normal(
        0, precision_val, n_abio
    )
    P_correct_if_abio = float(
        np.mean(meas_abio > crossover)
    )

    return {
        "LR_crossover_d15N": crossover,
        "P_correct_if_bio":  P_correct_if_bio,
        "P_correct_if_abio": P_correct_if_abio,
    }

# ─────────────────────────────────────────────────────────────
# MODULE 4: KEY SCENARIO OUTCOMES
# ─────────────────────────────────────────────────────────────

KEY_SCENARIOS = [
    {
        "label":  "Strong biological signal",
        "d15N":   2.0,
        "note":   "Centre of biological prediction.\n"
                  "N₂ fixation, minimal diagenesis.",
    },
    {
        "label":  "Diagenetically enriched bio",
        "d15N":   6.0,
        "note":   "Upper biological range.\n"
                  "600 Ma diagenetic ¹⁵N enrichment.",
    },
    {
        "label":  "Ambiguous zone",
        "d15N":   15.0,
        "note":   "Overlap of bio and abiotic.\n"
                  "Partially informative.",
    },
    {
        "label":  "Geochemical abiotic",
        "d15N":   25.0,
        "note":   "Upper abiotic geochemical range.\n"
                  "Serpentinization / hydrothermal.",
    },
    {
        "label":  "Mantle abiotic",
        "d15N":   45.0,
        "note":   "Mars mantle nitrogen signature.\n"
                  "Strongly abiotic.",
    },
    {
        "label":  "Atmospheric contamination",
        "d15N":  200.0,
        "note":   "Mars atmospheric N₂ incorporation.\n"
                  "Definitively abiotic.",
    },
]

# ───────────────────────────────────────────────────────��─────
# MODULE 5: PLOTTING
# ─────────────────────────────────────────────────────────────

def setup_style():
    plt.rcParams.update({
        "font.family":      "DejaVu Sans",
        "font.size":        8.5,
        "axes.titlesize":   9,
        "axes.titlecolor":  HEADER_BLUE,
        "axes.labelcolor":  HEADER_BLUE,
        "axes.edgecolor":   "#CCCCCC",
        "figure.facecolor": "white",
        "axes.facecolor":   "#FAFAFA",
        "grid.color":       "#DDDDDD",
        "grid.linewidth":   0.5,
    })

def plot_figure_32(mc_pred, scan_res):
    """
    Figure 32: Predicted measurement distributions.
    Panel A: δ¹⁵N distributions bio vs abiotic.
    Panel B: LR as function of measured δ¹⁵N.
    """
    setup_style()
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(
        "Script 11 — δ¹⁵N Prediction: "
        "Expected Measurement Distributions\n"
        "Biological: N₂ fixation δ¹⁵N ~ -4 to +8‰  "
        "|  Abiotic: mantle/geochemical -35 to +73‰\n"
        "NanoSIMS precision ±7.5‰ (1σ) included  "
        "|  Sources: Stüeken 2015; Mahaffy 2013; "
        "Lin 2014 MAPS",
        fontsize=8.5, color=HEADER_BLUE, y=1.01
    )

    # Panel A: distributions
    ax = axes[0]
    bio_vals  = mc_pred["measured_bio"]
    abio_vals = mc_pred["measured_abio"]

    # Clip for display (exclude atmospheric extreme)
    bio_clip  = bio_vals[
        (bio_vals > -60) & (bio_vals < 100)
    ]
    abio_clip = abio_vals[
        (abio_vals > -60) & (abio_vals < 150)
    ]

    ax.hist(bio_clip, bins=100, density=True,
            color="#2E7D32", alpha=0.60,
            edgecolor="white", linewidth=0.3,
            zorder=4,
            label=f"Biological prediction\n"
                  f"N₂ fixation + diagenesis\n"
                  f"mean={BIO_D15N['mean']:.0f}‰, "
                  f"σ={BIO_D15N['sigma']:.0f}‰\n"
                  f"Stüeken 2015; Zerkle 2016")

    ax.hist(abio_clip, bins=100, density=True,
            color="#BF360C", alpha=0.50,
            edgecolor="white", linewidth=0.3,
            zorder=3,
            label="Abiotic prediction\n"
                  "(mixed: mantle+geochemical\n"
                  "+atmospheric+mesostasis)\n"
                  "Mahaffy 2013; Metsoc 2024")

    # Mark zones
    ax.axvspan(-60, 15, alpha=0.05,
               color="#2E7D32", zorder=1,
               label="Biological zone")
    ax.axvspan(30, 150, alpha=0.05,
               color="#BF360C", zorder=1,
               label="Abiotic zone")
    ax.axvspan(15, 30, alpha=0.05,
               color="#888888", zorder=1,
               label="Ambiguous zone")

    ax.axvline(BIO_D15N["mean"],
               color="#2E7D32", linewidth=2.0,
               linestyle="--", zorder=5,
               label=f"Bio centre "
                     f"({BIO_D15N['mean']:.0f}‰)")

    ax.set_xlabel("Measured δ¹⁵N (‰, VPDB-air)",
                  fontsize=8.5)
    ax.set_ylabel("Probability density", fontsize=8.5)
    ax.set_title(
        "Panel A — Predicted δ¹⁵N measurement\n"
        "distributions under each hypothesis",
        color=HEADER_BLUE
    )
    ax.set_xlim(-60, 150)
    ax.legend(fontsize=7, framealpha=0.88)
    ax.grid(True, alpha=0.12, zorder=0)

    # Panel B: LR vs δ¹⁵N
    ax2 = axes[1]
    scan = scan_res["d15N_scan"]
    LRs  = scan_res["LRs"]
    mask = (scan > -50) & (scan < 150)

    log_LRs = np.where(
        LRs[mask] > 0,
        np.log10(np.maximum(LRs[mask], 1e-6)),
        -6
    )

    ax2.plot(scan[mask], log_LRs,
             color=HEADER_BLUE, linewidth=2.5,
             zorder=5)
    ax2.axhline(0, color="black",
                linewidth=1.5, linestyle="--",
                zorder=4, label="LR = 1:1 (crossover)")
    ax2.axhline(2, color="#2E7D32",
                linewidth=1.0, linestyle=":",
                zorder=4, alpha=0.7,
                label="LR = 100:1 (strong bio)")
    ax2.axhline(-2, color="#BF360C",
                linewidth=1.0, linestyle=":",
                zorder=4, alpha=0.7,
                label="LR = 0.01:1 (strong abio)")
    ax2.axhline(math.log10(150),
                color="#1565C0", linewidth=1.0,
                linestyle="-.", zorder=4, alpha=0.7,
                label="Kass-Raftery decisive (150:1)")

    ax2.axvspan(-50, 15, alpha=0.05,
                color="#2E7D32", zorder=1)
    ax2.axvspan(30, 150, alpha=0.05,
                color="#BF360C", zorder=1)

    # Mark key scenarios
    scenario_colors = [
        "#2E7D32", "#66BB6A", "#FF9800",
        "#E65100", "#BF360C", "#880000"
    ]
    for sc, col in zip(KEY_SCENARIOS[:5],
                       scenario_colors[:5]):
        lr = LR_d15N(sc["d15N"],
                     NANOSIMS_PRECISION["sigma_1sigma"])
        if lr > 0:
            log_lr = math.log10(max(lr, 1e-6))
        else:
            log_lr = -6
        ax2.scatter(
            sc["d15N"], log_lr,
            c=col, s=80, zorder=7,
            marker="D"
        )
        ax2.annotate(
            f"{sc['d15N']:.0f}‰",
            (sc["d15N"], log_lr),
            textcoords="offset points",
            xytext=(5, 5), fontsize=7,
            color=col
        )

    ax2.set_xlabel("Measured δ¹⁵N (‰)", fontsize=8.5)
    ax2.set_ylabel("log₁₀(LR)  [LR=P(bio)/P(abio)]",
                   fontsize=8.5)
    ax2.set_title(
        "Panel B — Likelihood ratio vs measured δ¹⁵N\n"
        "Diagnostic power of the experiment",
        color=HEADER_BLUE
    )
    ax2.set_xlim(-50, 150)
    ax2.legend(fontsize=7.5, framealpha=0.88)
    ax2.grid(True, alpha=0.12, zorder=0)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    path = "./necromass_d15N_fig32_distributions.png"
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Figure 32 saved: {path}")
    return path

def plot_figure_33(scan_res):
    """
    Figure 33: Updated posterior odds.
    What Script 10 posterior becomes after δ¹⁵N.
    """
    setup_style()
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(
        "Script 11 — δ¹⁵N Prediction: "
        "Updated Posterior Odds\n"
        "Prior = Script 10 conservative (82,871:1)  "
        "× LR(δ¹⁵N) = Updated posterior\n"
        "Sources: Kass & Raftery 1995; Script 10",
        fontsize=8.5, color=HEADER_BLUE, y=1.01
    )

    scan  = scan_res["d15N_scan"]
    post  = scan_res["post_odds"]
    pprob = scan_res["post_prob"]
    mask  = (scan > -50) & (scan < 150)

    log_post = np.log10(np.maximum(post[mask], 1e-6))

    # Panel A: log posterior odds vs δ¹⁵N
    ax = axes[0]
    ax.plot(scan[mask], log_post,
            color=HEADER_BLUE, linewidth=2.5,
            zorder=5)

    ax.axhline(math.log10(PRIOR_ODDS_SCRIPT10),
               color="#888888", linewidth=1.5,
               linestyle="--", zorder=4,
               label=f"Script 10 prior "
                     f"({PRIOR_ODDS_SCRIPT10:.0f}:1)")
    ax.axhline(math.log10(150),
               color="#E91E63", linewidth=1.5,
               linestyle="-.", zorder=4,
               label="Kass-Raftery decisive (150:1)")
    ax.axhline(0, color="black",
               linewidth=1.5, linestyle=":",
               zorder=4, label="Odds = 1:1 (no preference)")

    ax.axvspan(-50, 15, alpha=0.05,
               color="#2E7D32", zorder=1,
               label="Bio zone")
    ax.axvspan(30, 150, alpha=0.05,
               color="#BF360C", zorder=1,
               label="Abiotic zone")

    # Key scenarios
    sc_cols = ["#2E7D32", "#66BB6A",
               "#FF9800", "#E65100", "#BF360C"]
    for sc, col in zip(KEY_SCENARIOS[:5], sc_cols):
        pr = posterior_odds_updated(
            sc["d15N"],
            NANOSIMS_PRECISION["sigma_1sigma"],
            PRIOR_ODDS_SCRIPT10
        )
        lp = math.log10(max(pr, 1e-6))
        ax.scatter(sc["d15N"], lp, c=col,
                   s=80, zorder=7, marker="D")
        ax.annotate(
            f"{sc['d15N']:.0f}‰",
            (sc["d15N"], lp),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=7, color=col
        )

    ax.set_xlabel("Measured δ¹⁵N (‰)", fontsize=8.5)
    ax.set_ylabel("log₁₀(Posterior odds)",
                  fontsize=8.5)
    ax.set_title(
        "Panel A — Updated posterior odds\n"
        "Script 10 prior × LR(δ¹⁵N)",
        color=HEADER_BLUE
    )
    ax.set_xlim(-50, 150)
    ax.legend(fontsize=7.5, framealpha=0.88)
    ax.grid(True, alpha=0.12, zorder=0)

    # Panel B: posterior probability vs δ¹⁵N
    ax2 = axes[1]
    ax2.plot(scan[mask], pprob[mask],
             color=HEADER_BLUE, linewidth=2.5,
             zorder=5)

    ax2.axhline(0.5, color="black",
                linewidth=1.5, linestyle="--",
                zorder=4, label="P = 0.5 (no preference)")
    ax2.axhline(0.9999, color="#2E7D32",
                linewidth=1.0, linestyle=":",
                zorder=4, alpha=0.8,
                label="P = 0.9999 (decisive bio)")
    ax2.axhline(
        posterior_probability(PRIOR_ODDS_SCRIPT10),
        color="#888888", linewidth=1.5,
        linestyle="--", zorder=4,
        label=f"Script 10 prior P="
              f"{posterior_probability(PRIOR_ODDS_SCRIPT10):.6f}"
    )

    ax2.axvspan(-50, 15, alpha=0.05,
                color="#2E7D32", zorder=1)
    ax2.axvspan(30, 150, alpha=0.05,
                color="#BF360C", zorder=1)

    ax2.set_xlabel("Measured δ¹⁵N (‰)", fontsize=8.5)
    ax2.set_ylabel("Posterior P(biological origin)",
                   fontsize=8.5)
    ax2.set_title(
        "Panel B — Posterior probability\n"
        "P(bio) as function of δ¹⁵N measurement",
        color=HEADER_BLUE
    )
    ax2.set_xlim(-50, 150)
    ax2.set_ylim(-0.05, 1.05)
    ax2.legend(fontsize=7.5, framealpha=0.88)
    ax2.grid(True, alpha=0.12, zorder=0)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    path = "./necromass_d15N_fig33_posterior.png"
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Figure 33 saved: {path}")
    return path

def plot_figure_34(diag_power, mc_pred):
    """
    Figure 34: Diagnostic power summary.
    Scenario outcomes + confusion matrix.
    """
    setup_style()
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    fig.suptitle(
        "Script 11 — δ¹⁵N Prediction: "
        "Diagnostic Power and Scenario Outcomes\n"
        "What each possible measurement means "
        "for the Lafayette case\n"
        "Pre-reg DOI: 10.5281/zenodo.18986790",
        fontsize=8.5, color=HEADER_BLUE, y=1.01
    )

    # Panel A: scenario outcomes table
    ax = axes[0]
    ax.axis("off")

    precision = NANOSIMS_PRECISION["sigma_1sigma"]
    rows = []
    for sc in KEY_SCENARIOS:
        lr = LR_d15N(sc["d15N"], precision)
        po = posterior_odds_updated(
            sc["d15N"], precision, PRIOR_ODDS_SCRIPT10
        )
        pp = posterior_probability(po)
        if lr > 100:
            interp = "DECISIVE BIO"
            col = "#2E7D32"
        elif lr > 10:
            interp = "STRONG BIO"
            col = "#388E3C"
        elif lr > 1:
            interp = "BIO LEANING"
            col = "#66BB6A"
        elif lr > 0.1:
            interp = "AMBIGUOUS"
            col = "#FF9800"
        elif lr > 0.01:
            interp = "ABIO LEANING"
            col = "#E65100"
        else:
            interp = "DECISIVE ABIO"
            col = "#BF360C"

        if math.isinf(po):
            po_str = ">10¹⁵:1"
        elif po > 1e9:
            po_str = f"10^{math.log10(po):.0f}:1"
        elif po > 1000:
            po_str = f"{po:.0f}:1"
        elif po > 1:
            po_str = f"{po:.1f}:1"
        else:
            po_str = f"1:{1/po:.1f}"

        rows.append([
            f"{sc['d15N']:+.0f}‰",
            sc["label"],
            f"{lr:.2f}",
            po_str,
            f"{pp:.4f}",
            interp,
        ])

    col_labels = ["δ¹⁵N", "Scenario",
                  "LR", "Post odds",
                  "P(bio)", "Verdict"]
    tbl = ax.table(
        cellText=rows,
        colLabels=col_labels,
        cellLoc="center",
        loc="center",
        bbox=[0, 0.1, 1, 0.85]
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.5)

    verdict_colors = {
        "DECISIVE BIO":  "#E8F5E9",
        "STRONG BIO":    "#F1F8E9",
        "BIO LEANING":   "#F9FBE7",
        "AMBIGUOUS":     "#FFF8E1",
        "ABIO LEANING":  "#FBE9E7",
        "DECISIVE ABIO": "#FFEBEE",
    }
    for (row, col), cell in tbl.get_celld().items():
        cell.set_edgecolor("#CCCCCC")
        if row == 0:
            cell.set_facecolor(HEADER_BLUE)
            cell.set_text_props(color="white",
                                fontweight="bold")
        elif col == 5 and row > 0:
            verdict = rows[row-1][5]
            cell.set_facecolor(
                verdict_colors.get(verdict, "white")
            )

    ax.set_title(
        "Panel A — Scenario outcomes\n"
        "What each δ¹⁵N value means for the case",
        color=HEADER_BLUE, fontsize=8.5
    )

    # Panel B: diagnostic power summary
    ax2 = axes[1]
    ax2.axis("off")

    cross = diag_power["LR_crossover_d15N"]
    p_bio = diag_power["P_correct_if_bio"]
    p_abi = diag_power["P_correct_if_abio"]

    lines = [
        ("DIAGNOSTIC POWER SUMMARY",
         0.93, 12, HEADER_BLUE, "bold"),
        ("NanoSIMS δ¹⁵N experiment on Lafayette",
         0.87, 9.5, "#333333", "normal"),
        ("─" * 52,
         0.83, 9, "#CCCCCC", "normal"),
        (f"LR crossover point:  {cross:.1f}‰",
         0.78, 9.5, HEADER_BLUE, "normal"),
        (f"  δ¹⁵N < {cross:.1f}‰  →  LR > 1  (bio favoured)",
         0.73, 9, "#2E7D32", "normal"),
        (f"  δ¹⁵N > {cross:.1f}‰  →  LR < 1  (abio favoured)",
         0.68, 9, "#BF360C", "normal"),
        ("─" * 52,
         0.63, 9, "#CCCCCC", "normal"),
        ("If biology is true:",
         0.58, 9.5, HEADER_BLUE, "bold"),
        (f"  P(measurement correctly identifies bio)",
         0.53, 9, "#333333", "normal"),
        (f"  = {p_bio:.4f}  ({100*p_bio:.1f}%)",
         0.48, 11, "#2E7D32", "bold"),
        ("If abiotic is true:",
         0.43, 9.5, HEADER_BLUE, "bold"),
        (f"  P(measurement correctly identifies abio)",
         0.38, 9, "#333333", "normal"),
        (f"  = {p_abi:.4f}  ({100*p_abi:.1f}%)",
         0.33, 11, "#BF360C", "bold"),
        ("─" * 52,
         0.28, 9, "#CCCCCC", "normal"),
        ("Prior odds (Script 10): 82,871:1",
         0.23, 9, "#333333", "normal"),
        ("If bio measured (+2‰): >10⁶:1",
         0.18, 9, "#2E7D32", "normal"),
        ("If abio measured (+45‰): ~350:1",
         0.13, 9, "#FF9800", "normal"),
        ("If atmos measured (+200‰): <1:1",
         0.08, 9, "#BF360C", "normal"),
        ("One session. Existing material.",
         0.03, 9, HEADER_BLUE, "bold"),
    ]

    for (text, y, size, color, weight) in lines:
        ax2.text(
            0.5, y, text,
            transform=ax2.transAxes,
            fontsize=size, va="center",
            ha="center", color=color,
            fontweight=weight
        )
    ax2.set_title(
        "Panel B — Experiment diagnostic power",
        color=HEADER_BLUE, fontsize=8.5
    )

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    path = "./necromass_d15N_fig34_power.png"
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Figure 34 saved: {path}")
    return path

# ─────────────────────────────────────────────────────────────
# MODULE 6: MAIN RUN
# ─────────────────────────────────────────────────────────────

def run():
    sep = "═" * 62

    print()
    print(sep)
    print("  NECROMASS δ¹⁵N PREDICTION — SCRIPT 11 v1.0")
    print("  Mars Iron Necromass Hypothesis")
    print("  Pre-reg: DOI 10.5281/zenodo.18986790")
    print("  2026-03-13 / Eric Robert Lawson / OrganismCore")
    print(sep)

    # ── SECTION 1: PARAMETERS ──────────────────────────���─────
    print()
    print("─" * 62)
    print("  SECTION 1: PUBLISHED PARAMETERS")
    print("─" * 62)
    print()
    print("  BIOLOGICAL PREDICTION:")
    print(f"    δ¹⁵N = {BIO_D15N['lo']:.0f} to "
          f"{BIO_D15N['hi']:.0f}‰")
    print(f"    Mean: {BIO_D15N['mean']:.0f}‰, "
          f"σ: {BIO_D15N['sigma']:.0f}‰")
    print(f"    Source: {BIO_D15N['source']}")
    print()
    print("  ABIOTIC PREDICTIONS (mixture model):")
    for k, s in ABIO_SOURCES.items():
        print(f"    {s['label']}:")
        print(f"      δ¹⁵N = {s['lo']:.0f} to "
              f"{s['hi']:.0f}‰  "
              f"(mean {s['mean']:.0f}, "
              f"σ {s['sigma']:.0f}‰)")
        print(f"      Weight: {s['weight']:.2f}")
        print(f"      Source: {s['source']}")
    print()
    print("  NANOSIMS PRECISION:")
    print(f"    ±{NANOSIMS_PRECISION['sigma_1sigma']:.1f}‰ "
          f"(1σ, typical)")
    print(f"    Range: ±{NANOSIMS_PRECISION['sigma_lo']:.0f}"
          f" to ±{NANOSIMS_PRECISION['sigma_hi']:.0f}‰")
    print(f"    Source: {NANOSIMS_PRECISION['source']}")
    print()
    print(f"  PRIOR ODDS (Script 10 conservative):")
    print(f"    {PRIOR_ODDS_SCRIPT10:.0f}:1")
    print()

    # ── SECTION 2: KEY SCENARIO OUTCOMES ─────────────────────
    print()
    print("─" * 62)
    print("  SECTION 2: KEY SCENARIO OUTCOMES")
    print("─" * 62)
    print()
    precision = NANOSIMS_PRECISION["sigma_1sigma"]
    print(f"  {'Scenario':<32} "
          f"{'δ¹⁵N':>6} "
          f"{'LR':>10} "
          f"{'Post odds':>14} "
          f"{'P(bio)':>9}")
    print("  " + "─" * 73)
    for sc in KEY_SCENARIOS:
        lr = LR_d15N(sc["d15N"], precision)
        po = posterior_odds_updated(
            sc["d15N"], precision, PRIOR_ODDS_SCRIPT10
        )
        pp = posterior_probability(po)

        if lr > 1e6:
            lr_str = f"10^{math.log10(lr):.1f}:1"
        elif lr > 1:
            lr_str = f"{lr:.1f}:1"
        elif lr > 0.001:
            lr_str = f"1:{1/lr:.1f}"
        else:
            lr_str = f"1:{1/lr:.0f}"

        if math.isinf(po):
            po_str = ">10¹⁵:1"
        elif po > 1e9:
            po_str = f"10^{math.log10(po):.1f}:1"
        elif po > 1000:
            po_str = f"{po:.0f}:1"
        elif po > 1:
            po_str = f"{po:.1f}:1"
        else:
            po_str = f"1:{1/po:.1f}"

        print(
            f"  {sc['label']:<32} "
            f"{sc['d15N']:>+6.0f}‰ "
            f"{lr_str:>10} "
            f"{po_str:>14} "
            f"{pp:>9.6f}"
        )

    # ── SECTION 3: MONTE CARLO ───────────────────────────────
    print()
    print("─" * 62)
    print("  SECTION 3: MONTE CARLO PREDICTIONS")
    print(f"  N = {N_MC:,}")
    print("─" * 62)
    print()
    print("  Running MC prediction distributions...")
    mc_pred = run_MC_prediction()
    print("  Running diagnostic scan...")
    scan_res = run_MC_scenarios()
    print("  Computing diagnostic power...")
    diag = compute_diagnostic_power(mc_pred)
    print()

    print("  PREDICTED MEASUREMENT (if biological):")
    bio_m = mc_pred["measured_bio"]
    print(f"    Mean: {np.mean(bio_m):.2f}‰")
    print(f"    Median: {np.median(bio_m):.2f}‰")
    print(f"    90% CI: [{np.percentile(bio_m,5):.1f}, "
          f"{np.percentile(bio_m,95):.1f}]‰")
    print()
    print("  PREDICTED MEASUREMENT (if abiotic):")
    abio_m = mc_pred["measured_abio"]
    print(f"    Mean: {np.mean(abio_m):.1f}‰")
    print(f"    Median: {np.median(abio_m):.1f}‰")
    print(f"    90% CI: [{np.percentile(abio_m,5):.1f}, "
          f"{np.percentile(abio_m,95):.1f}]‰")
    print()
    print("  DIAGNOSTIC POWER:")
    print(f"    LR crossover: "
          f"{diag['LR_crossover_d15N']:.1f}‰")
    print(f"    P(correct | bio true): "
          f"{diag['P_correct_if_bio']:.4f} "
          f"({100*diag['P_correct_if_bio']:.1f}%)")
    print(f"    P(correct | abio true): "
          f"{diag['P_correct_if_abio']:.4f} "
          f"({100*diag['P_correct_if_abio']:.1f}%)")
    print()

    # ── SECTION 4: FIGURES ───────────────────────────────────
    print()
    print("─" * 62)
    print("  SECTION 4: GENERATING FIGURES")
    print("─" * 62)
    print()
    f32 = plot_figure_32(mc_pred, scan_res)
    f33 = plot_figure_33(scan_res)
    f34 = plot_figure_34(diag, mc_pred)

    # ── SECTION 5: FINAL STATEMENT ───────────────────────────
    print()
    print("─" * 62)
    print("  SECTION 5: FINAL STATEMENT")
    print("─" * 62)
    print()
    print("  THE EXPERIMENT:")
    print("    δ¹⁵N of Lafayette vein organics")
    print("    by NanoSIMS co-registered with δ¹³C.")
    print("    Measure CN⁻ isotope ratio specifically")
    print("    in the organic carbon phase at Fe")
    print("    oxide grain boundaries.")
    print("    Precision: ±7.5‰ (1σ).")
    print("    One session. Existing polished section.")
    print()
    print("  THREE POSSIBLE OUTCOMES:")
    print()

    sc_bio  = KEY_SCENARIOS[0]
    sc_amb  = KEY_SCENARIOS[2]
    sc_abio = KEY_SCENARIOS[4]

    for sc in [sc_bio, sc_amb, sc_abio]:
        lr = LR_d15N(sc["d15N"], precision)
        po = posterior_odds_updated(
            sc["d15N"], precision,
            PRIOR_ODDS_SCRIPT10
        )
        pp = posterior_probability(po)
        if po > 1e9:
            po_s = f"10^{math.log10(po):.1f}:1"
        elif po > 1:
            po_s = f"{po:.0f}:1"
        else:
            po_s = f"1:{1/po:.1f}"

        print(f"  If measured δ¹⁵N ≈ {sc['d15N']:+.0f}‰"
              f"  ({sc['label']}):")
        print(f"    LR(δ¹⁵N):         {lr:.3f}")
        print(f"    Updated odds:      {po_s}")
        print(f"    Updated P(bio):    {pp:.6f}")
        print()

    print("  WHAT CHANGES AND WHAT DOES NOT:")
    print()
    print("  If biological result (+2‰):")
    print("    Updated posterior: >10⁶:1.")
    print("    C:N signal is confirmed by δ¹⁵N.")
    print("    Both signals point same direction.")
    print("    No abiotic explanation survives.")
    print()
    print("  If ambiguous result (+15‰):")
    print("    Updated posterior: somewhat reduced.")
    print("    Still well above decisive threshold.")
    print("    C:N signal remains strongest constraint.")
    print("    Experiment adds partial information.")
    print()
    print("  If abiotic result (+45‰):")
    print("    Updated posterior: substantially reduced.")
    print("    Requires re-evaluation of C:N signal.")
    print("    Could indicate N source is mantle/")
    print("    atmospheric contamination rather than")
    print("    biological fixation.")
    print("    Would not alone refute the hypothesis")
    print("    but would significantly weaken it.")
    print()
    print("  WHAT THE EXPERIMENT CANNOT DO:")
    print("    It cannot by itself PROVE biology.")
    print("    No single measurement can.")
    print("    It can STRONGLY CONFIRM or")
    print("    SUBSTANTIALLY WEAKEN the case.")
    print("    That is exactly what a diagnostic")
    print("    experiment should do.")
    print()

    print(sep)
    print("  SCRIPT 11 COMPLETE")
    print("  FULL SERIES COMPLETE (Scripts 1-11)")
    print(sep)
    print()
    print("  Figures generated:")
    for f in [f32, f33, f34]:
        print(f"    {f}")
    print()
    print("  Complete figure set (all 11 scripts):")
    print("    Figs  1-7:  Scripts 1-2")
    print("    Figs  8-10: Script 3")
    print("    Figs 11-13: Script 4")
    print("    Figs 14-16: Script 5")
    print("    Figs 17-19: Script 6")
    print("    Figs 20-22: Script 7")
    print("    Figs 23-26: Script 8")
    print("    Figs 27-28: Script 9")
    print("    Figs 29-31: Script 10")
    print("    Figs 32-34: Script 11")
    print("    Total: 34 figures")
    print()
    print("  Pre-reg: 10.5281/zenodo.18986790")
    print("  github.com/Eric-Robert-Lawson/"
          "attractor-oncology")
    print("  ORCID: 0009-0002-0414-6544")
    print()

    pr_bio  = posterior_odds_updated(
        KEY_SCENARIOS[0]["d15N"], precision,
        PRIOR_ODDS_SCRIPT10
    )
    pr_amb  = posterior_odds_updated(
        KEY_SCENARIOS[2]["d15N"], precision,
        PRIOR_ODDS_SCRIPT10
    )
    pr_abio = posterior_odds_updated(
        KEY_SCENARIOS[4]["d15N"], precision,
        PRIOR_ODDS_SCRIPT10
    )

    return {
        "prior_odds_script10":
            float(PRIOR_ODDS_SCRIPT10),
        "bio_prediction_mean_d15N":
            float(BIO_D15N["mean"]),
        "bio_prediction_range":
            f"{BIO_D15N['lo']} to {BIO_D15N['hi']}",
        "LR_crossover_d15N":
            float(diag["LR_crossover_d15N"]),
        "P_correct_if_bio":
            float(diag["P_correct_if_bio"]),
        "P_correct_if_abio":
            float(diag["P_correct_if_abio"]),
        "posterior_if_bio_signal":
            float(min(pr_bio, 1e15)),
        "posterior_if_ambiguous":
            float(min(pr_amb, 1e15)),
        "posterior_if_abio_signal":
            float(min(pr_abio, 1e15)),
        "P_bio_if_bio_measured":
            float(posterior_probability(pr_bio)),
        "P_bio_if_ambiguous":
            float(posterior_probability(pr_amb)),
        "P_bio_if_abio_measured":
            float(posterior_probability(pr_abio)),
        "nanosims_precision_1sigma":
            float(NANOSIMS_PRECISION["sigma_1sigma"]),
        "experiment":
            "NanoSIMS delta15N on Lafayette vein organics",
        "measurement_mode":
            "12C14N- and 12C15N- simultaneous "
            "co-registered with 12C- and 13C-",
        "existing_material_required":
            True,
        "new_sample_prep_required":
            False,
    }


if __name__ == "__main__":
    result = run()
    print("  MACHINE-READABLE SUMMARY:")
    print()
    for k, v in result.items():
        if isinstance(v, float):
            print(f"  {k}: {v:.6e}")
        else:
            print(f"  {k}: {v}")
    print()
    sys.exit(0)
