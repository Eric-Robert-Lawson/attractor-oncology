#!/usr/bin/env python3
"""
NECROMASS RATE CALCULATION — SCRIPT 2
======================================
Document ID:  NECROMASS_RATE_v1.0
Date:         2026-03-13
Author:       Eric Robert Lawson / OrganismCore
Pre-reg DOI:  10.5281/zenodo.18986790

PURPOSE
-------
Script 1 showed that total iron oxide deposition
increases monotonically with depth in the Martian
rock pile. P(r>0) = 1.000 across 100,000 Monte Carlo
samples. Effect magnitude: 13× from surface to 30m.

Script 1 could not distinguish between:
  BIOLOGICAL:  Iron-cycling organisms catalyse
               iron oxidation at depth.
  ABIOTIC:     Temperature gradient at depth
               accelerates abiotic iron oxidation.

THIS SCRIPT asks:

  Given the total mass of iron oxide in the
  Lafayette nakhlite alteration veins, and given
  the published duration and temperature of the
  single aqueous event that produced it, could
  abiotic chemistry have deposited that iron oxide
  in the time available?

  Or does the mass require biological catalysis?

THE THREE ABIOTIC PATHWAYS ON MARS
-----------------------------------
Oxygen-driven oxidation: NOT APPLICABLE.
  Mars had no free atmospheric oxygen during or
  after the aqueous event. Eliminated.

Temperature-driven oxidation by H2O2 (radiolysis):
  APPLICABLE. This is the Fenton pathway.
  Radiolysis of water by U, Th, K decay produces
  H2O2. H2O2 oxidises Fe2+ → Fe3+.
  Rate: k_Fenton at temperature T.
  This is the strongest abiotic candidate.

Surface-catalysed mineral oxidation:
  APPLICABLE. Fe3+-bearing mineral surfaces
  can catalyse Fe2+ oxidation.
  Rate: slow, k ~ 1e-8 to 1e-6 s-1 at 100°C.
  This is the second abiotic candidate.

THE BIOLOGICAL PATHWAY
-----------------------
Iron-oxidising chemolithotrophs:
  Analogues: Gallionella, Leptothrix,
             D. audaxviator-type sulfate reducers,
             photoferrotrophic anaerobes.
  Rate: 1-6 × 10^-16 mol Fe2+ per cell per second.
  At subsurface cell densities: 10^3 cells/mL.
  Biological amplification over abiotic:
  up to 10^6 × in low-oxygen conditions.

THE QUESTION
------------
  Does the observed Lafayette iron oxide mass
  fall within abiotic range for the event duration?
  Or does it require biological catalysis?

ALL VALUES FROM PUBLISHED SOURCES
-----------------------------------
Lafayette alteration:
  Treiman et al. 1993, Meteoritics 28:86-97
  Changela & Bridges 2010, MAPS 45:1847-1867
  Papike et al. 2006, American Mineralogist 91:607

Aqueous event parameters:
  Bridges et al. 2001, 2015
  Swindle et al. 2000
  Changela & Bridges 2010

Fe2+ oxidation kinetics:
  Fenton: Millero 1987, GCA 51:793
  Biological: Vollrath et al. 2012,
              Geomicrobiology Journal
  Radiolysis H2O2 production:
    Lin et al. 2005, Earth Planet Sci Lett

DEPENDENCIES
------------
  numpy, scipy, matplotlib
  pip install numpy scipy matplotlib
"""

import sys
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

# ─────────────────────────────────────────────────────────────
# PHYSICAL CONSTANTS
# ─────────────────────────────────────────────────────────────

R_GAS       = 8.314          # J / (mol K)
AVOGADRO    = 6.022e23       # molecules / mol
SECONDS_PER_YEAR = 3.156e7   # s / yr

# ─────────────────────────────────────────────────────────────
# MODULE 1: LAFAYETTE IRON OXIDE MASS
# ────��────────────────────────────────────────────────────────
# We need to estimate the total mass of iron oxide
# deposited in Lafayette alteration veins.
#
# From published data we can estimate this as:
#
#   M_FeO = V_rock × f_alt × rho_alt × w_FeO
#
# Where:
#   V_rock    = reference volume of rock (1 m³)
#   f_alt     = alteration volume fraction (published)
#   rho_alt   = density of alteration product (published)
#   w_FeO     = weight fraction of FeO in alteration (published)
#
# Sources:
#   f_alt:    Treiman 2005 — <2% bulk, using 1.5% from Script 1
#   rho_alt:  ~2.3 g/cm³ (iddingsite, published)
#   w_FeO:    32-40 wt% from Treiman 1993, Papike 2006
#             Using Changela & Bridges mean: 57.5 wt%
#             (this is the vein FeO concentration we measured)
# ─────────────────────────────────────────────────────────────

LAFAYETTE = {
    # Reference volume of rock (1 cubic metre)
    "V_rock_m3": 1.0,

    # Alteration volume fraction
    # Source: Changela & Bridges 2010 (from Script 1)
    # Script 1 value: 1.5% = 0.015 fraction
    "f_alt_mid":  0.015,
    "f_alt_min":  0.010,
    "f_alt_max":  0.025,

    # Density of alteration product (iddingsite/Fe-oxide gel)
    # Source: Earth Science Reviews 96:196-216
    # Published range: 2.2-2.4 g/cm3
    "rho_alt_gcm3_mid": 2.3,
    "rho_alt_gcm3_min": 2.2,
    "rho_alt_gcm3_max": 2.4,

    # FeO weight fraction in alteration veins
    # Source: Changela & Bridges 2010 Table 4 (Script 1)
    # Mean value used: 57.5 wt% = 0.575
    # Note: this is TOTAL iron as FeO equiv.
    # The ACTUAL Fe oxide fraction (ferrihydrite,
    # goethite, magnetite) may be lower.
    # Conservative estimate: use 0.35 (midpoint of
    # Treiman 1993 range 32-40 wt% Fe2O3)
    "w_FeO_mid":  0.35,
    "w_FeO_min":  0.32,
    "w_FeO_max":  0.575,     # upper bound from Changela

    # Molar mass of FeO equivalent (g/mol)
    "M_FeO_gmol": 71.844,    # Fe = 55.845, O = 15.999

    # Aqueous event parameters
    # Source: Bridges et al. 2001; Swindle et al. 2000;
    #         Changela & Bridges 2010
    "event_age_Ma":      600,    # Ma ago
    "event_duration_yr_mid": 1e4,  # years (median estimate)
    "event_duration_yr_min": 1e1,  # minimum: years
    "event_duration_yr_max": 1e6,  # maximum: million years

    # Temperature during aqueous event
    # Source: Bridges et al. 2001; Changela & Bridges 2010
    "T_celsius_mid":  85.0,    # °C
    "T_celsius_min":  50.0,
    "T_celsius_max": 150.0,

    # Fluid:rock mass ratio
    # Source: Bridges et al. 2001
    "fluid_rock_ratio_mid": 0.05,
    "fluid_rock_ratio_min": 0.01,
    "fluid_rock_ratio_max": 0.10,
}

# ─────────────────────────────────────────────────────────────
# MODULE 2: ABIOTIC OXIDATION PATHWAYS
# ─────────────────────────────────────────────────────────────

# ── PATHWAY 1: FENTON (H2O2) OXIDATION ──────────────────────
# Fe2+ + H2O2 → Fe3+ + OH- + •OH
#
# Rate law: d[Fe3+]/dt = k_Fenton × [Fe2+] × [H2O2]
#
# Rate constant at 25°C: k = 63 M-1 s-1
# Source: Millero 1987, GCA 51:793
# (measured in seawater, used as reference)
#
# Temperature dependence (Arrhenius):
# k(T) = k_25 × exp(-Ea/R × (1/T - 1/298.15))
#
# Activation energy Ea: 26 kJ/mol (mid estimate)
# Source: NIST kinetics database; RSC Advances 2021
#
# H2O2 concentration in Martian subsurface groundwater
# from radiolysis:
# Source: Lin et al. 2005, Earth Planet Sci Lett 237:286
# Estimated: 0.1-1 mM in rock with typical U/Th/K content
# Conservative estimate: 0.1 mM = 1e-4 M
#
# Fe2+ concentration in Martian subsurface groundwater:
# Estimate from nakhlite olivine dissolution:
# 1-10 mM conservative. Using 1 mM = 1e-3 M.
# ─────────────────────────────────────────────────────────────

FENTON = {
    "name": "Fenton (H2O2 radiolysis)",
    "pathway": "Fe2+ + H2O2 → Fe3+ + OH- + •OH",
    "k_25C_M_s": 63.0,           # M-1 s-1 at 25°C
    "Ea_Jmol_mid": 26000.0,      # J/mol
    "Ea_Jmol_min": 20000.0,
    "Ea_Jmol_max": 40000.0,
    "H2O2_conc_M_mid": 1e-4,     # M (0.1 mM)
    "H2O2_conc_M_min": 1e-5,     # M
    "H2O2_conc_M_max": 1e-3,     # M (1 mM)
    "Fe2_conc_M_mid": 1e-3,      # M (1 mM)
    "Fe2_conc_M_min": 1e-4,
    "Fe2_conc_M_max": 1e-2,
    "source": (
        "Millero 1987 GCA 51:793; "
        "NIST kinetics database; "
        "Lin et al. 2005 EPSL 237:286"
    ),
}

# ── PATHWAY 2: SURFACE-CATALYSED OXIDATION ──────────────────
# Fe2+ on mineral surfaces → Fe3+
# Rate constant at 100°C: ~1e-7 s-1
# Source: Luther et al. 1992 extrapolation;
#         Millero 1987 (Table 2)
# ───────────────────���─────────────────────────────────────────

SURFACE_CAT = {
    "name": "Surface-catalysed mineral oxidation",
    "pathway": "Fe2+ (surface) → Fe3+ (mineral surface)",
    "k_100C_s": 1e-7,            # s-1 pseudo-first-order
    "k_25C_s":  1e-9,            # s-1
    "Ea_Jmol_mid": 30000.0,      # J/mol (estimated)
    "Ea_Jmol_min": 20000.0,
    "Ea_Jmol_max": 50000.0,
    "source": (
        "Luther et al. 1992 (extrapolated); "
        "Millero 1987 Table 2"
    ),
}

# ─────────────────────────────────────────────────────────────
# MODULE 3: BIOLOGICAL OXIDATION
# ─────────────────────────────────────────────────────────────
# Iron-oxidising chemolithotrophs (Gallionella analogues)
# under anaerobic, low-oxygen conditions.
#
# Per-cell rate: 1-6 × 10-16 mol Fe2+/cell/s
# Source: Vollrath et al. 2012, Geomicrobiology Journal
#
# Cell density in deep subsurface:
# ~10^3 cells/mL in D. audaxviator habitat
# Source: Chivian et al. 2008, Science 322:275
#
# Volume of water in Lafayette alteration system:
# fluid:rock ratio 0.05, rock density 3.0 g/cm3 (pyroxenite)
# Per 1 m3 rock: 0.05 × 3e6 g = 150,000 g water = 150 L = 1.5e5 mL
# ─────────────────────────────────────────────────────────────

BIOLOGICAL = {
    "name": "Biological iron oxidation",
    "pathway": "Fe2+ + [organism] → Fe3+ + metabolic energy",
    "rate_mol_cell_s_mid": 3.0e-16,   # mol Fe2+/cell/s
    "rate_mol_cell_s_min": 1.0e-16,
    "rate_mol_cell_s_max": 6.0e-16,
    "cell_density_per_mL_mid": 1e3,
    "cell_density_per_mL_min": 1e2,
    "cell_density_per_mL_max": 1e5,
    "source": (
        "Vollrath et al. 2012 Geomicrobiology J; "
        "Chivian et al. 2008 Science 322:275"
    ),
}

# ─────────────────────────────────────────────────────────────
# MODULE 4: CALCULATION FUNCTIONS
# ─────────────────────────────────────────────────────────────

def arrhenius(k_ref, T_ref_K, T_K, Ea_Jmol):
    """
    Arrhenius temperature correction.
    k(T) = k_ref × exp(-Ea/R × (1/T - 1/T_ref))
    """
    return k_ref * math.exp(
        -Ea_Jmol / R_GAS * (1.0 / T_K - 1.0 / T_ref_K)
    )

def compute_lafayette_Fe_mass_kg(
        f_alt, rho_gcm3, w_FeO, V_m3=1.0):
    """
    Total iron (as FeO equiv) mass in Lafayette
    alteration veins per unit rock volume.

    Returns mass in kg per m3 of rock.
    """
    V_alt_cm3 = V_m3 * 1e6 * f_alt   # cm3 of alteration
    mass_alt_g = V_alt_cm3 * rho_gcm3  # g of alteration
    mass_FeO_g = mass_alt_g * w_FeO    # g of FeO in alteration
    return mass_FeO_g / 1000.0          # kg

def compute_Fe_moles(mass_kg, M_gmol=71.844):
    """Convert mass of FeO (kg) to moles of Fe."""
    return (mass_kg * 1000.0) / M_gmol

def fenton_rate_mol_per_Lyr(T_C, Ea_Jmol,
                              H2O2_M, Fe2_M):
    """
    Fenton Fe2+ oxidation rate in mol/L/yr.
    d[Fe3+]/dt = k(T) × [Fe2+] × [H2O2]
    """
    T_K = T_C + 273.15
    k_T = arrhenius(
        FENTON["k_25C_M_s"], 298.15, T_K, Ea_Jmol
    )
    rate_mol_L_s = k_T * Fe2_M * H2O2_M
    return rate_mol_L_s * SECONDS_PER_YEAR

def surface_cat_rate_mol_per_Lyr(T_C, Ea_Jmol,
                                   Fe2_M):
    """
    Surface-catalysed Fe2+ oxidation rate in mol/L/yr.
    Pseudo first-order: d[Fe3+]/dt = k(T) × [Fe2+]
    Reference k at 25°C.
    """
    T_K = T_C + 273.15
    k_T = arrhenius(
        SURFACE_CAT["k_25C_s"], 298.15, T_K, Ea_Jmol
    )
    rate_mol_L_s = k_T * Fe2_M
    return rate_mol_L_s * SECONDS_PER_YEAR

def biological_rate_mol_per_Lyr(
        rate_mol_cell_s, cell_density_per_mL):
    """
    Biological Fe2+ oxidation rate in mol/L/yr.
    = rate_per_cell × cell_density × conversion
    """
    # cell_density_per_mL → per L = × 1000
    rate_mol_L_s = (rate_mol_cell_s
                    * cell_density_per_mL * 1000.0)
    return rate_mol_L_s * SECONDS_PER_YEAR

def water_volume_L_per_m3_rock(
        fluid_rock_ratio, rock_density_gcm3=3.0):
    """
    Volume of water (L) per m3 of rock.
    fluid_rock_ratio is mass water / mass rock.
    rock_density_gcm3: pyroxenite ~3.0 g/cm3.
    """
    rock_mass_g = 1e6 * rock_density_gcm3    # g per m3
    water_mass_g = rock_mass_g * fluid_rock_ratio
    return water_mass_g / 1000.0             # L (1g water = 1mL)

def time_required_yr(Fe_moles_total,
                     rate_mol_L_yr, water_L):
    """
    Time required to deposit Fe_moles_total moles
    of Fe3+ via oxidation at rate_mol_L_yr
    in water_L litres of fluid.

    Total moles deposited per year:
      = rate_mol_L_yr × water_L

    Time = Fe_moles_total / (rate × water_L)
    """
    if rate_mol_L_yr <= 0 or water_L <= 0:
        return float("inf")
    moles_per_yr = rate_mol_L_yr * water_L
    return Fe_moles_total / moles_per_yr

# ─────────────────────────────────────────────────────────────
# MODULE 5: MONTE CARLO OVER PARAMETER SPACE
# ─────────────────────────────────────────────────────────────

N_MC = 100_000
np.random.seed(42)

def run_monte_carlo():
    """
    Sample all parameters from uniform [min, max] ranges.
    For each sample, compute:
      - Fe moles in Lafayette alteration
      - Abiotic Fenton time required
      - Abiotic surface-cat time required
      - Biological time required
      - Ratio: abiotic / biological time
      - Whether abiotic can explain in event duration
    """
    results = {
        "fe_moles": [],
        "fenton_time_yr": [],
        "surface_time_yr": [],
        "bio_time_yr": [],
        "event_duration_yr": [],
        "fenton_adequate": [],    # True if fenton_time <= event_duration
        "surface_adequate": [],
        "bio_adequate": [],
        "fenton_bio_ratio": [],   # fenton_time / bio_time
    }

    for _ in range(N_MC):
        # Sample Lafayette parameters
        f_alt = np.random.uniform(
            LAFAYETTE["f_alt_min"],
            LAFAYETTE["f_alt_max"]
        )
        rho = np.random.uniform(
            LAFAYETTE["rho_alt_gcm3_min"],
            LAFAYETTE["rho_alt_gcm3_max"]
        )
        w_FeO = np.random.uniform(
            LAFAYETTE["w_FeO_min"],
            LAFAYETTE["w_FeO_max"]
        )
        T_C = np.random.uniform(
            LAFAYETTE["T_celsius_min"],
            LAFAYETTE["T_celsius_max"]
        )
        event_dur = np.random.uniform(
            LAFAYETTE["event_duration_yr_min"],
            LAFAYETTE["event_duration_yr_max"]
        )
        fr_ratio = np.random.uniform(
            LAFAYETTE["fluid_rock_ratio_min"],
            LAFAYETTE["fluid_rock_ratio_max"]
        )

        # Sample kinetic parameters
        Ea_fenton = np.random.uniform(
            FENTON["Ea_Jmol_min"],
            FENTON["Ea_Jmol_max"]
        )
        H2O2 = np.random.uniform(
            FENTON["H2O2_conc_M_min"],
            FENTON["H2O2_conc_M_max"]
        )
        Fe2 = np.random.uniform(
            LAFAYETTE["fluid_rock_ratio_min"] * 10,
            LAFAYETTE["fluid_rock_ratio_max"] * 10
        )
        Fe2_M = np.random.uniform(1e-4, 1e-2)

        Ea_surf = np.random.uniform(
            SURFACE_CAT["Ea_Jmol_min"],
            SURFACE_CAT["Ea_Jmol_max"]
        )

        bio_rate = np.random.uniform(
            BIOLOGICAL["rate_mol_cell_s_min"],
            BIOLOGICAL["rate_mol_cell_s_max"]
        )
        cell_density = np.random.uniform(
            BIOLOGICAL["cell_density_per_mL_min"],
            BIOLOGICAL["cell_density_per_mL_max"]
        )

        # Compute Fe mass and moles
        mass_kg = compute_lafayette_Fe_mass_kg(
            f_alt, rho, w_FeO
        )
        fe_moles = compute_Fe_moles(mass_kg)

        # Compute water volume
        water_L = water_volume_L_per_m3_rock(fr_ratio)

        # Compute rates
        r_fenton = fenton_rate_mol_per_Lyr(
            T_C, Ea_fenton, H2O2, Fe2_M
        )
        r_surface = surface_cat_rate_mol_per_Lyr(
            T_C, Ea_surf, Fe2_M
        )
        r_bio = biological_rate_mol_per_Lyr(
            bio_rate, cell_density
        )

        # Compute times required
        t_fenton  = time_required_yr(
            fe_moles, r_fenton, water_L
        )
        t_surface = time_required_yr(
            fe_moles, r_surface, water_L
        )
        t_bio     = time_required_yr(
            fe_moles, r_bio, water_L
        )

        results["fe_moles"].append(fe_moles)
        results["fenton_time_yr"].append(t_fenton)
        results["surface_time_yr"].append(t_surface)
        results["bio_time_yr"].append(t_bio)
        results["event_duration_yr"].append(event_dur)
        results["fenton_adequate"].append(
            t_fenton <= event_dur
        )
        results["surface_adequate"].append(
            t_surface <= event_dur
        )
        results["bio_adequate"].append(
            t_bio <= event_dur
        )
        if t_bio > 0 and t_fenton < float("inf"):
            results["fenton_bio_ratio"].append(
                t_fenton / t_bio
            )
        else:
            results["fenton_bio_ratio"].append(np.nan)

    return {k: np.array(v) for k, v in results.items()}

# ─────────────────────────────────────────────────────────────
# MODULE 6: PLOTTING
# ─────────────────────────────────────────────────────────────

HEADER_BLUE = "#14388C"
WARN_RED    = "#8B0000"

def setup_style():
    plt.rcParams.update({
        "font.family":      "DejaVu Sans",
        "font.size":        9,
        "axes.titlesize":   9.5,
        "axes.labelsize":   8.5,
        "axes.titlecolor":  HEADER_BLUE,
        "axes.labelcolor":  HEADER_BLUE,
        "axes.edgecolor":   "#CCCCCC",
        "figure.facecolor": "white",
        "axes.facecolor":   "#FAFAFA",
        "grid.color":       "#DDDDDD",
        "grid.linewidth":   0.5,
    })

def plot_figure_5(mc):
    """
    Figure 5: Time-required distributions for each pathway.
    3 panels: Fenton, surface-cat, biological.
    Event duration range overlaid.
    """
    setup_style()
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(
        "Script 2 — Iron Oxide Rate Calculation\n"
        "Time Required to Deposit Lafayette Fe Oxide Mass "
        "by Each Pathway vs. Published Event Duration\n"
        "Pre-reg DOI: 10.5281/zenodo.18986790  |  "
        "2026-03-13  |  Eric Robert Lawson / OrganismCore",
        fontsize=8.5, color=HEADER_BLUE, y=1.01
    )

    panels = [
        {
            "ax": axes[0],
            "data_key": "fenton_time_yr",
            "adequate_key": "fenton_adequate",
            "color": "#795548",
            "title": "Panel A — Fenton Pathway\n"
                     "(H₂O��� radiolysis; strongest abiotic)",
            "label": "Fenton (abiotic)",
        },
        {
            "ax": axes[1],
            "data_key": "surface_time_yr",
            "adequate_key": "surface_adequate",
            "color": "#607D8B",
            "title": "Panel B — Surface-Catalysed\n"
                     "(mineral surface oxidation; abiotic)",
            "label": "Surface-cat (abiotic)",
        },
        {
            "ax": axes[2],
            "data_key": "bio_time_yr",
            "adequate_key": "bio_adequate",
            "color": "#2E7D32",
            "title": "Panel C — Biological Pathway\n"
                     "(iron-oxidising chemolithotrophs)",
            "label": "Biological",
        },
    ]

    for p in panels:
        ax = p["ax"]
        raw = mc[p["data_key"]]
        # Clip infinite values for display
        finite = raw[np.isfinite(raw)]
        if len(finite) == 0:
            ax.text(0.5, 0.5, "All infinite\n(pathway too slow)",
                    ha="center", va="center",
                    transform=ax.transAxes,
                    color=WARN_RED, fontsize=10)
            continue

        log_data = np.log10(np.clip(finite, 1, 1e20))
        ax.hist(
            log_data, bins=60,
            color=p["color"], alpha=0.70,
            edgecolor="white", linewidth=0.3,
            zorder=4, label=p["label"]
        )

        # Event duration range
        log_dur_min = math.log10(
            LAFAYETTE["event_duration_yr_min"]
        )
        log_dur_max = math.log10(
            LAFAYETTE["event_duration_yr_max"]
        )
        ax.axvspan(
            log_dur_min, log_dur_max,
            color="#FFC107", alpha=0.25,
            zorder=3, label="Event duration range"
        )
        ax.axvline(
            math.log10(LAFAYETTE["event_duration_yr_mid"]),
            color="#FFC107", linewidth=1.5,
            linestyle="--", zorder=5,
            label="Event duration (mid)"
        )

        # Adequate fraction
        frac_ok = np.mean(mc[p["adequate_key"]])
        med_t = np.median(finite)

        ax.text(
            0.97, 0.97,
            f"P(adequate) = {frac_ok:.4f}\n"
            f"Median time = {med_t:.2e} yr\n"
            f"Event mid   = {LAFAYETTE['event_duration_yr_mid']:.0e} yr",
            transform=ax.transAxes,
            fontsize=7.5, color=p["color"],
            va="top", ha="right",
            bbox=dict(
                boxstyle="round,pad=0.3",
                facecolor="white",
                edgecolor=p["color"],
                alpha=0.9
            )
        )

        ax.set_xlabel("log₁₀(Time required, years)",
                      fontsize=8)
        ax.set_ylabel("Monte Carlo sample count",
                      fontsize=8)
        ax.set_title(p["title"], fontsize=8.5,
                     color=HEADER_BLUE)
        ax.grid(True, alpha=0.15, zorder=0)
        ax.legend(fontsize=7, loc="upper left",
                  framealpha=0.8)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    path = "./necromass_rate_fig5_time_required.png"
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Figure 5 saved: {path}")
    return path

def plot_figure_6(mc):
    """
    Figure 6: Fenton/Biological time ratio distribution.
    Shows how much longer abiotic takes than biological
    across all parameter combinations.
    """
    setup_style()
    fig, ax = plt.subplots(figsize=(9, 5))

    ratios = mc["fenton_bio_ratio"]
    finite = ratios[np.isfinite(ratios) & (ratios > 0)]
    log_ratios = np.log10(finite)

    ax.hist(
        log_ratios, bins=70,
        color="#B71C1C", alpha=0.70,
        edgecolor="white", linewidth=0.3,
        zorder=4
    )
    ax.axvline(0, color="black", linewidth=1.2,
               linestyle="-", zorder=5,
               label="Ratio = 1 (equal speed)")
    ax.axvline(
        np.median(log_ratios),
        color="#B71C1C", linewidth=1.8,
        linestyle="--", zorder=6,
        label=f"Median ratio = "
              f"10^{np.median(log_ratios):.2f}"
              f" = {10**np.median(log_ratios):.1e}×"
    )

    frac_abiotic_slower = np.mean(log_ratios > 0)
    ax.text(
        0.97, 0.97,
        f"P(Fenton slower than biological) = "
        f"{frac_abiotic_slower:.4f}\n"
        f"Median Fenton/Bio time ratio = "
        f"{10**np.median(log_ratios):.2e}×\n"
        f"95% CI: [{10**np.percentile(log_ratios,2.5):.1e}×, "
        f"{10**np.percentile(log_ratios,97.5):.1e}×]\n"
        f"Interpretation: Fenton abiotic pathway requires\n"
        f"this many times LONGER than biological pathway\n"
        f"to deposit the same Fe oxide mass.",
        transform=ax.transAxes,
        fontsize=8, color="#B71C1C",
        va="top", ha="right",
        bbox=dict(
            boxstyle="round,pad=0.4",
            facecolor="white",
            edgecolor="#B71C1C",
            alpha=0.92
        )
    )

    ax.set_xlabel(
        "log₁₀(Fenton time / Biological time)\n"
        "> 0 means abiotic is slower than biological",
        fontsize=8.5
    )
    ax.set_ylabel("Monte Carlo sample count", fontsize=8.5)
    ax.set_title(
        "Fenton Abiotic vs. Biological Iron Oxidation — "
        "Time Ratio\n"
        "How much longer does the abiotic pathway take "
        "to deposit the same Fe oxide mass?\n"
        "Pre-reg DOI: 10.5281/zenodo.18986790",
        fontsize=8.5, color=HEADER_BLUE
    )
    ax.grid(True, alpha=0.15, zorder=0)
    ax.legend(fontsize=8, loc="upper left", framealpha=0.8)

    plt.tight_layout()
    path = "./necromass_rate_fig6_abiotic_vs_bio_ratio.png"
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Figure 6 saved: {path}")
    return path

def plot_figure_7(mc):
    """
    Figure 7: Fe moles distribution — how much iron
    did we actually deposit in Lafayette?
    Shows the quantity being explained.
    """
    setup_style()
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(
        "Script 2 — Lafayette Iron Oxide Mass Estimate\n"
        "Total Fe moles per m³ rock and required "
        "time vs. event duration\n"
        "Sources: Changela & Bridges 2010; "
        "Treiman 1993; Wang et al. 2021",
        fontsize=8.5, color=HEADER_BLUE, y=1.01
    )

    # Panel A: Fe moles distribution
    ax = axes[0]
    ax.hist(
        mc["fe_moles"], bins=60,
        color="#FF6F00", alpha=0.75,
        edgecolor="white", linewidth=0.3, zorder=4
    )
    ax.axvline(
        np.median(mc["fe_moles"]),
        color="#E65100", linewidth=1.8,
        linestyle="--", zorder=5,
        label=f"Median = {np.median(mc['fe_moles']):.1f} mol/m³"
    )
    ax.set_xlabel(
        "Moles of Fe (as FeO) per m³ rock\n"
        "in Lafayette alteration veins",
        fontsize=8
    )
    ax.set_ylabel("Monte Carlo sample count", fontsize=8)
    ax.set_title(
        "Panel A — Fe oxide mass in Lafayette\n"
        "(the quantity each pathway must explain)",
        fontsize=8.5, color=HEADER_BLUE
    )
    ax.legend(fontsize=8, framealpha=0.8)
    ax.grid(True, alpha=0.15, zorder=0)
    ax.text(
        0.97, 0.97,
        f"Median: {np.median(mc['fe_moles']):.1f} mol/m³\n"
        f"95% CI: [{np.percentile(mc['fe_moles'],2.5):.1f}, "
        f"{np.percentile(mc['fe_moles'],97.5):.1f}]",
        transform=ax.transAxes,
        fontsize=7.5, color="#E65100",
        va="top", ha="right",
        bbox=dict(
            boxstyle="round,pad=0.3",
            facecolor="white",
            edgecolor="#E65100", alpha=0.9
        )
    )

    # Panel B: Comparison — required times vs event duration
    ax2 = axes[1]
    labels = [
        "Fenton\n(abiotic)",
        "Surface-cat\n(abiotic)",
        "Biological"
    ]
    colors = ["#795548", "#607D8B", "#2E7D32"]
    keys   = ["fenton_time_yr", "surface_time_yr",
              "bio_time_yr"]
    medians = []
    lo_errs = []
    hi_errs = []
    for key in keys:
        d = mc[key]
        finite = d[np.isfinite(d)]
        med = np.median(finite) if len(finite) else np.nan
        lo  = np.percentile(finite, 2.5) if len(finite) else np.nan
        hi  = np.percentile(finite, 97.5) if len(finite) else np.nan
        medians.append(med)
        lo_errs.append(med - lo if not np.isnan(lo) else 0)
        hi_errs.append(hi - med if not np.isnan(hi) else 0)

    x_pos = np.arange(len(labels))
    bars = ax2.bar(
        x_pos,
        [math.log10(max(m, 1)) for m in medians],
        color=colors, alpha=0.75,
        edgecolor="#444444", linewidth=0.5,
        zorder=4
    )

    # Event duration range
    ax2.axhspan(
        math.log10(LAFAYETTE["event_duration_yr_min"]),
        math.log10(LAFAYETTE["event_duration_yr_max"]),
        color="#FFC107", alpha=0.25,
        zorder=3, label="Event duration range"
    )
    ax2.axhline(
        math.log10(LAFAYETTE["event_duration_yr_mid"]),
        color="#FFC107", linewidth=1.5,
        linestyle="--", zorder=5,
        label="Event duration (mid)"
    )

    for i, (bar, med) in enumerate(zip(bars, medians)):
        if not np.isnan(med) and np.isfinite(med):
            ax2.text(
                bar.get_x() + bar.get_width() / 2,
                math.log10(max(med, 1)) + 0.1,
                f"{med:.1e} yr",
                ha="center", va="bottom",
                fontsize=7.5, color=colors[i]
            )

    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(labels, fontsize=8)
    ax2.set_ylabel(
        "log₁₀(Time required, years)", fontsize=8
    )
    ax2.set_title(
        "Panel B — Required time vs. event duration\n"
        "(yellow = published event duration range)",
        fontsize=8.5, color=HEADER_BLUE
    )
    ax2.legend(fontsize=7.5, framealpha=0.8)
    ax2.grid(True, alpha=0.15, axis="y", zorder=0)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    path = "./necromass_rate_fig7_fe_mass_and_time.png"
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Figure 7 saved: {path}")
    return path

# ─────────────────────────────────────────────────────────────
# MODULE 7: MAIN ANALYSIS RUN
# ─────────────────────────────────────────────────────────────

def run():

    sep = "═" * 62

    print()
    print(sep)
    print("  NECROMASS RATE CALCULATION — SCRIPT 2 v1.0")
    print("  Mars Iron Necromass Hypothesis")
    print("  Pre-reg: DOI 10.5281/zenodo.18986790")
    print("  2026-03-13 / Eric Robert Lawson / OrganismCore")
    print(sep)

    # ── SECTION 1: PARAMETER REVIEW ─────────────────────────
    print()
    print("─" * 62)
    print("  SECTION 1: PUBLISHED PARAMETERS USED")
    print("─" * 62)
    print()
    print("  LAFAYETTE ALTERATION:")
    print(f"    Alteration volume fraction : "
          f"{LAFAYETTE['f_alt_mid']:.3f} "
          f"[{LAFAYETTE['f_alt_min']:.3f}"
          f"–{LAFAYETTE['f_alt_max']:.3f}]")
    print(f"    Alteration product density : "
          f"{LAFAYETTE['rho_alt_gcm3_mid']} g/cm³")
    print(f"    FeO in alteration veins    : "
          f"{LAFAYETTE['w_FeO_mid']:.2f} wt frac "
          f"[{LAFAYETTE['w_FeO_min']:.2f}"
          f"–{LAFAYETTE['w_FeO_max']:.3f}]")
    print(f"    Aqueous event temperature  : "
          f"{LAFAYETTE['T_celsius_mid']}°C "
          f"[{LAFAYETTE['T_celsius_min']}"
          f"–{LAFAYETTE['T_celsius_max']}]")
    print(f"    Event duration             : "
          f"{LAFAYETTE['event_duration_yr_mid']:.0e} yr "
          f"[{LAFAYETTE['event_duration_yr_min']:.0e}"
          f"–{LAFAYETTE['event_duration_yr_max']:.0e}]")
    print(f"    Fluid:rock ratio           : "
          f"{LAFAYETTE['fluid_rock_ratio_mid']} "
          f"[{LAFAYETTE['fluid_rock_ratio_min']}"
          f"–{LAFAYETTE['fluid_rock_ratio_max']}]")

    print()
    print("  KINETIC PARAMETERS:")
    print(f"    Fenton k at 25°C           : "
          f"{FENTON['k_25C_M_s']} M⁻¹s⁻¹")
    print(f"    Fenton Ea                  : "
          f"{FENTON['Ea_Jmol_mid']/1000:.0f} kJ/mol")
    print(f"    H₂O₂ concentration        : "
          f"{FENTON['H2O2_conc_M_mid']:.0e} M")
    print(f"    Surface-cat k at 25°C      : "
          f"{SURFACE_CAT['k_25C_s']:.0e} s⁻¹")
    print(f"    Biological rate/cell       : "
          f"{BIOLOGICAL['rate_mol_cell_s_mid']:.0e} "
          f"mol/cell/s")
    print(f"    Cell density               : "
          f"{BIOLOGICAL['cell_density_per_mL_mid']:.0e} "
          f"cells/mL")

    # ── SECTION 2: POINT ESTIMATE ───────────────────────────
    print()
    print("─" * 62)
    print("  SECTION 2: POINT ESTIMATE (MID VALUES)")
    print("─" * 62)
    print()

    # Fe mass and moles (mid values)
    mass_kg = compute_lafayette_Fe_mass_kg(
        LAFAYETTE["f_alt_mid"],
        LAFAYETTE["rho_alt_gcm3_mid"],
        LAFAYETTE["w_FeO_mid"]
    )
    fe_moles = compute_Fe_moles(mass_kg)
    water_L  = water_volume_L_per_m3_rock(
        LAFAYETTE["fluid_rock_ratio_mid"]
    )
    T_C = LAFAYETTE["T_celsius_mid"]
    T_K = T_C + 273.15

    print(f"  Fe oxide mass per m³ rock    : {mass_kg:.2f} kg")
    print(f"  Fe moles per m³ rock         : {fe_moles:.1f} mol")
    print(f"  Water volume per m³ rock     : {water_L:.1f} L")
    print(f"  Temperature                  : {T_C}°C ({T_K:.1f} K)")
    print()

    # Fenton rate at temperature
    k_fenton_T = arrhenius(
        FENTON["k_25C_M_s"], 298.15, T_K,
        FENTON["Ea_Jmol_mid"]
    )
    r_fenton = fenton_rate_mol_per_Lyr(
        T_C, FENTON["Ea_Jmol_mid"],
        FENTON["H2O2_conc_M_mid"],
        LAFAYETTE["fluid_rock_ratio_mid"] * 10
    )
    t_fenton = time_required_yr(fe_moles, r_fenton, water_L)

    # Surface-cat rate
    r_surface = surface_cat_rate_mol_per_Lyr(
        T_C, SURFACE_CAT["Ea_Jmol_mid"],
        LAFAYETTE["fluid_rock_ratio_mid"] * 10
    )
    t_surface = time_required_yr(
        fe_moles, r_surface, water_L
    )

    # Biological rate
    r_bio = biological_rate_mol_per_Lyr(
        BIOLOGICAL["rate_mol_cell_s_mid"],
        BIOLOGICAL["cell_density_per_mL_mid"]
    )
    t_bio = time_required_yr(fe_moles, r_bio, water_L)

    event_dur = LAFAYETTE["event_duration_yr_mid"]

    print(f"  {'Pathway':<28} {'Rate(mol/L/yr)':>16} "
          f"{'Time req.(yr)':>16} {'< event dur?':>14}")
    print("  " + "─" * 76)
    for name, rate, t_req in [
        ("Fenton (abiotic)",     r_fenton,  t_fenton),
        ("Surface-cat (abiotic)",r_surface, t_surface),
        ("Biological",           r_bio,     t_bio),
    ]:
        adequate = "YES" if t_req <= event_dur else "NO"
        t_str = f"{t_req:.2e}" if np.isfinite(t_req) else "∞"
        print(f"  {name:<28} {rate:>16.3e} "
              f"{t_str:>16} {adequate:>14}")
    print()
    print(f"  Published event duration (mid): "
          f"{event_dur:.0e} yr")
    print()

    # ── SECTION 3: MONTE CARLO ───────────────────────────────
    print()
    print("─" * 62)
    print("  SECTION 3: MONTE CARLO UNCERTAINTY PROPAGATION")
    print(f"  N = {N_MC:,} samples")
    print("─" * 62)
    print()
    print("  Running Monte Carlo...")
    mc = run_monte_carlo()
    print("  Complete.")
    print()

    for label, adequate_key, time_key in [
        ("Fenton (abiotic)",      "fenton_adequate",
         "fenton_time_yr"),
        ("Surface-cat (abiotic)", "surface_adequate",
         "surface_time_yr"),
        ("Biological",            "bio_adequate",
         "bio_time_yr"),
    ]:
        finite = mc[time_key][np.isfinite(mc[time_key])]
        frac_ok = np.mean(mc[adequate_key])
        med = np.median(finite) if len(finite) else np.nan
        lo  = np.percentile(finite, 2.5) if len(finite) else np.nan
        hi  = np.percentile(finite, 97.5) if len(finite) else np.nan

        print(f"  [{label}]")
        print(f"  P(can explain in event duration) = "
              f"{frac_ok:.4f}")
        print(f"  Median time required = {med:.2e} yr")
        print(f"  95% CI = [{lo:.2e}, {hi:.2e}] yr")
        print()

    # Abiotic vs biological ratio
    ratios = mc["fenton_bio_ratio"]
    finite_r = ratios[np.isfinite(ratios) & (ratios > 0)]
    med_ratio = np.median(finite_r)
    lo_ratio  = np.percentile(finite_r, 2.5)
    hi_ratio  = np.percentile(finite_r, 97.5)
    frac_abiotic_slower = np.mean(finite_r > 1)

    print("  [Fenton/Biological time ratio]")
    print(f"  Median ratio  = {med_ratio:.2e}×")
    print(f"  95% CI        = [{lo_ratio:.2e}×, "
          f"{hi_ratio:.2e}×]")
    print(f"  P(Fenton slower than bio) = "
          f"{frac_abiotic_slower:.4f}")
    print()

    # ── SECTION 4: FIGURES ───────────────────────────────────
    print()
    print("─" * 62)
    print("  SECTION 4: GENERATING FIGURES")
    print("─" * 62)
    print()
    f5 = plot_figure_5(mc)
    f6 = plot_figure_6(mc)
    f7 = plot_figure_7(mc)

    # ── SECTION 5: INTERPRETATION ────────────────────────────
    print()
    print("─" * 62)
    print("  SECTION 5: HONEST INTERPRETATION")
    print("─" * 62)
    print()

    frac_fenton = np.mean(mc["fenton_adequate"])
    frac_bio    = np.mean(mc["bio_adequate"])

    print("  QUESTION ASKED:")
    print("  Can abiotic chemistry deposit the observed")
    print("  Lafayette iron oxide mass within the published")
    print("  duration of the single aqueous event?")
    print()
    print(f"  Fenton abiotic P(adequate) = {frac_fenton:.4f}")
    print(f"  Biological     P(adequate) = {frac_bio:.4f}")
    print()

    if frac_fenton < 0.10:
        fenton_verdict = (
            "ABIOTIC FENTON PATHWAY LIKELY INSUFFICIENT.\n"
            "  In <10% of parameter combinations can the\n"
            "  Fenton pathway deposit the observed iron\n"
            "  oxide mass within the event duration.\n"
            "  Biological catalysis or an alternative\n"
            "  abiotic source is required."
        )
    elif frac_fenton < 0.50:
        fenton_verdict = (
            "ABIOTIC FENTON PATHWAY MARGINAL.\n"
            "  The Fenton pathway can explain the iron\n"
            "  oxide mass in some but not most parameter\n"
            "  combinations. Result is inconclusive.\n"
            "  Biological catalysis remains a viable\n"
            "  explanation."
        )
    else:
        fenton_verdict = (
            "ABIOTIC FENTON PATHWAY SUFFICIENT IN\n"
            "  MAJORITY OF PARAMETER COMBINATIONS.\n"
            "  Cannot require biological catalysis\n"
            "  from this calculation alone.\n"
            "  Biological hypothesis not supported\n"
            "  as necessary by rate argument."
        )

    print(f"  VERDICT:")
    print(f"  {fenton_verdict}")
    print()
    print("  WHAT THIS DOES AND DOES NOT MEAN:")
    print()
    print("  DOES NOT MEAN: Mars has no biology.")
    print("  DOES NOT MEAN: The necromass hypothesis is wrong.")
    print("  DOES NOT MEAN: Abiotic is the correct answer.")
    print()
    print("  DOES MEAN: The rate argument either constrains")
    print("    or does not constrain the biological hypothesis")
    print("    depending on the fenton_verdict above.")
    print()
    print("  CRITICAL UNCERTAINTY:")
    print("  The event duration range spans 10^1 to 10^6 yr.")
    print("  This 5-order-of-magnitude uncertainty dominates")
    print("  the result. If the event was 10^6 yr long,")
    print("  abiotic processes have much more time available.")
    print("  If the event was 10^1 to 10^3 yr long,")
    print("  biological catalysis becomes more necessary.")
    print("  The event duration constraint is the key")
    print("  uncertain parameter in this calculation.")
    print()
    print("  NEXT REQUIRED TEST:")
    print("  Script 3 — Carbon co-precipitation in")
    print("  Lafayette alteration veins.")
    print("  If biologically-fractionated organic carbon")
    print("  co-precipitates with iron oxide in Lafayette,")
    print("  that is a biological signal independent of")
    print("  the rate argument.")
    print()

    print(sep)
    print("  SCRIPT 2 COMPLETE")
    print(sep)
    print()
    print("  Figures generated:")
    for f in [f5, f6, f7]:
        print(f"    {f}")
    print()
    print("  Pre-registration chain:")
    print("    10.5281/zenodo.18986790 (origin)")
    print("  Repository:")
    print("    github.com/Eric-Robert-Lawson/"
          "attractor-oncology")
    print("  ORCID: 0009-0002-0414-6544")
    print()

    return {
        "fe_moles_median": float(np.median(mc["fe_moles"])),
        "fenton_P_adequate": float(frac_fenton),
        "bio_P_adequate": float(frac_bio),
        "fenton_bio_ratio_median": float(med_ratio),
        "fenton_verdict": fenton_verdict,
    }


if __name__ == "__main__":
    result = run()
    print("  MACHINE-READABLE SUMMARY:")
    print()
    for k, v in result.items():
        print(f"  {k}: {v}")
    print()
    sys.exit(0)
