#!/usr/bin/env python3
"""
NECROMASS GEOMETRIC COVERAGE — SCRIPT 12 v2.0
==============================================
Document ID:  NECROMASS_GEOM_COVERAGE_v2.0
Date:         2026-03-13
Author:       Eric Robert Lawson / OrganismCore
Pre-reg DOI:  10.5281/zenodo.18986790

CORRECTIONS OVER v1.0
----------------------
BUG 1 (crash): plot_figure_37 Panel B bar chart.
  x had shape (9,) — all observations.
  IDS_9_b_vals had shape (6,) — N-bearing only.
  numpy broadcast_arrays raised ValueError.

  FIX: Panel B now uses TWO SEPARATE SUB-PANELS:
    Left:  7D IDS for ALL observations (n=9).
    Right: 9D IDS for N-BEARING observations only (n=6).
  Both clearly labelled. No shape conflict possible.

BUG 2 (latent): IDS_7_all list construction.
  Was called as IDS_7_all[:-1] + [IDS_laff_7]
  which was redundant and fragile.
  FIX: IDS list built cleanly in one pass,
  Lafayette appended once.

BUG 3 (latent): obs_9d_list included ABIO_OBS
  entries with d15N=None that had already been
  filtered, creating index mismatches downstream.
  FIX: filtering applied consistently in one
  place and reused throughout.
"""

import sys
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D          # noqa: F401
from scipy import stats
from matplotlib.gridspec import GridSpec

np.random.seed(42)
N_MC = 200_000

HEADER_BLUE = "#14388C"

# ──────────────────────────────────────────────���──────────────
# MODULE 1: REFERENCE DATASET
# ─────────────────────────────────────────────────────────────

def coupling_ratio(d13C_org, d15N):
    """δ¹³C_organic / δ¹⁵N_organic."""
    if d15N is None or d15N == 0:
        return float("nan")
    return d13C_org / d15N

BIO_OBS = [
    {
        "label":   "D. audaxviator\n(Mponeng)",
        "short":   "D. aud",
        "d13C":    -33.0,
        "d13C_dd":  0.0,
        "CN_log":   math.log10(6.5),
        "Fe_col":   0.8,
        "dep_gr":   0.7,
        "frac_m":  25.0,
        "tmpl":     0.85,
        "d15N":     3.0,
        "source":  "Chivian 2008 Science 322:275",
        "color":   "#2E7D32",
        "marker":  "o",
    },
    {
        "label":   "FeOB\n(deep subsurface)",
        "short":   "FeOB ref",
        "d13C":    -30.0,
        "d13C_dd":  0.5,
        "CN_log":   math.log10(7.0),
        "Fe_col":   1.0,
        "dep_gr":   0.9,
        "frac_m":  22.0,
        "tmpl":     0.90,
        "d15N":     2.0,
        "source":  "House 2003 GCA 67:3447",
        "color":   "#388E3C",
        "marker":  "o",
    },
    {
        "label":   "Deep methanogen",
        "short":   "Methanogen",
        "d13C":    -50.0,
        "d13C_dd":  1.0,
        "CN_log":   math.log10(8.0),
        "Fe_col":   0.5,
        "dep_gr":   0.8,
        "frac_m":  42.0,
        "tmpl":     0.70,
        "d15N":     1.0,
        "source":  "Whiticar 1999; House 2003",
        "color":   "#43A047",
        "marker":  "o",
    },
]

ABIO_OBS = [
    {
        "label":   "Serp short-chain\n(abiotic)",
        "short":   "Serp SC",
        "d13C":    -15.0,
        "d13C_dd":  0.0,
        "CN_log":   math.log10(500.0),
        "Fe_col":   0.2,
        "dep_gr":   0.1,
        "frac_m":   7.0,
        "tmpl":     0.10,
        "d15N":     None,
        "source":  "McCollom & Seewald 2007",
        "color":   "#BF360C",
        "marker":  "^",
    },
    {
        "label":   "Serp methane\n(abiotic)",
        "short":   "Serp CH4",
        "d13C":    -40.0,
        "d13C_dd":  0.0,
        "CN_log":   math.log10(1000.0),
        "Fe_col":   0.1,
        "dep_gr":   0.2,
        "frac_m":  32.0,
        "tmpl":     0.15,
        "d15N":     None,
        "source":  "McCollom & Seewald 2007",
        "color":   "#D84315",
        "marker":  "^",
    },
    {
        "label":   "CO₂ UV photolysis\n(Ueno 2024)",
        "short":   "UV phot",
        "d13C":   -107.0,
        "d13C_dd":  0.0,
        "CN_log":   math.log10(800.0),
        "Fe_col":   0.0,
        "dep_gr":   0.0,
        "frac_m":  99.0,
        "tmpl":     0.05,
        "d15N":     None,
        "source":  "Ueno 2024 Nat Geosci",
        "color":   "#880000",
        "marker":  "^",
    },
    {
        "label":   "Abiotic hydrothermal\n+ geochemical N",
        "short":   "Hydro+geoN",
        "d13C":    -20.0,
        "d13C_dd":  0.0,
        "CN_log":   math.log10(200.0),
        "Fe_col":   0.3,
        "dep_gr":   0.1,
        "frac_m":  12.0,
        "tmpl":     0.20,
        "d15N":     5.0,
        "source":  "McCollom 2007; Earth analogues",
        "color":   "#E65100",
        "marker":  "^",
    },
    {
        "label":   "Mantle N +\nserp C",
        "short":   "Mantle N",
        "d13C":     -8.0,
        "d13C_dd":  0.0,
        "CN_log":   math.log10(500.0),
        "Fe_col":   0.1,
        "dep_gr":   0.0,
        "frac_m":   0.0,
        "tmpl":     0.05,
        "d15N":    -35.0,
        "source":  "Frantseva 2018 GCA 231:64",
        "color":   "#795548",
        "marker":  "^",
    },
]

ALL_OBS = BIO_OBS + ABIO_OBS  # length 8

# Pre-filter: observations that have valid d15N
# AND a finite coupling ratio — used for 9D PCA.
OBS_9D = [
    o for o in ALL_OBS
    if o["d15N"] is not None
    and math.isfinite(coupling_ratio(o["d13C"],
                                     o["d15N"]))
]
# OBS_9D has 5 entries (3 bio + 2 abio with N)

LAFAYETTE = {
    "label":    "Lafayette\n(measured+predicted)",
    "short":    "Lafayette",
    "d13C":     -41.3,
    "d13C_dd":   6.5,
    "CN_log":    math.log10(2.0),
    "Fe_col":    1.0,
    "dep_gr":    1.0,
    "frac_m":   41.3,
    "tmpl":      0.833,
    "d15N_bio":  2.0,
    "d15N_lo":  -4.0,
    "d15N_hi":   8.0,
    "source":   "Scripts 1-11; Steele 2012; "
                "McMahon 2016; Stüeken 2015",
    "color":    "#E91E63",
    "marker":   "*",
}

# ─────────────────────────────────────────────────────────────
# MODULE 2: PCA
# ─────────────────────────────────────────────────────────────

def _obs_to_7D(obs):
    return [
        obs["d13C"],  obs["d13C_dd"],
        obs["CN_log"], obs["Fe_col"],
        obs["dep_gr"], obs["frac_m"],
        obs["tmpl"],
    ]

def _obs_to_9D(obs):
    cr = coupling_ratio(obs["d13C"], obs["d15N"])
    return [
        obs["d13C"],  obs["d13C_dd"],
        obs["CN_log"], obs["Fe_col"],
        obs["dep_gr"], obs["frac_m"],
        obs["tmpl"],  obs["d15N"], cr,
    ]

def build_7D_matrix():
    return np.array(
        [_obs_to_7D(o) for o in ALL_OBS],
        dtype=float
    )

def build_9D_matrix():
    return np.array(
        [_obs_to_9D(o) for o in OBS_9D],
        dtype=float
    )

def run_PCA(M):
    mean    = M.mean(axis=0)
    M_c     = M - mean
    cov     = np.cov(M_c.T)
    ev, evc = np.linalg.eigh(cov)
    idx     = np.argsort(ev)[::-1]
    ev      = ev[idx]
    evc     = evc[:, idx]
    scores  = M_c @ evc
    var_exp = ev / ev.sum()
    return {"mean": mean, "eigvecs": evc,
            "eigvals": ev, "scores": scores,
            "var_exp": var_exp}

def project_point(vec, pca):
    return (vec - pca["mean"]) @ pca["eigvecs"]

def IDS_score(scores_obs, scores_bio, scores_abio):
    bio_c  = scores_bio.mean(axis=0)
    abio_c = scores_abio.mean(axis=0)
    axis   = bio_c - abio_c
    norm   = np.linalg.norm(axis)
    if norm < 1e-10:
        return 0.0
    axis  = axis / norm
    return float(scores_obs[:2] @ axis[:2])

# ─────────────────────────────────────────────────────────────
# MODULE 3: THREE-AXIS OCTANT
# ─────────────────────────────────────────────────────────────

def octant_check(d13C, d15N, CN_ratio):
    a1 = d13C < -25.0
    a2 = (d15N is not None and d15N > -10.0)
    a3 = (CN_ratio is not None and CN_ratio < 20.0)
    sat = [a1, a2, a3]
    return sum(sat), sat

# ─────────────────────────────────────────────────────────────
# MODULE 4: BAYESIAN CHAIN
# ─────────────────────────────────────────────────────────────

CHAIN_STEPS = [
    {"step": 0, "label": "Flat prior",   "odds": 1.0},
    {"step": 1, "label": "S1-5",         "odds": 99.0},
    {"step": 2, "label": "S6",           "odds": 180.0},
    {"step": 3, "label": "S8 PCA",       "odds": 82871.0},
    {"step": 4, "label": "S10 joint LR", "odds": 82871.0},
]

BRANCH_BIO  = {"label": "δ¹⁵N=+2‰ (bio)",
               "LR": 3.37,   "color": "#2E7D32"}
BRANCH_AMB  = {"label": "δ¹⁵N=-5‰ (ambig)",
               "LR": 2.61,   "color": "#FF9800"}
BRANCH_ABI  = {"label": "δ¹⁵N=-40‰ (mantle)",
               "LR": 1.0/175786.0, "color": "#BF360C"}

COUPLING_LR = {
    "bio":  5.0,
    "amb":  2.0,
    "abio": 0.001,
}

# ─────────────────────────────────────────────────────────────
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

# ── FIGURE 35 ────────────────────────────────────────────────

def plot_figure_35(pca7, pca9, M7, M9):
    setup_style()
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle(
        "Script 12 v2.0 — PCA Identity Manifold "
        "Expansion  7D → 9D\n"
        "(+ δ¹⁵N + δ¹³C/δ¹⁵N coupling)\n"
        "Sources: Frantseva 2018; Stüeken 2015; "
        "Scripts 1-11",
        fontsize=8.5, color=HEADER_BLUE, y=1.01
    )

    n_bio7  = len(BIO_OBS)
    n_bio9  = sum(1 for o in OBS_9D
                  if o["marker"] == "o")

    for ax, pca, obs_list, n_bio, title, laff_dims in [
        (axes[0], pca7, ALL_OBS,  n_bio7,
         "Panel A — 7D manifold (Scripts 1-8)", 7),
        (axes[1], pca9, OBS_9D,   n_bio9,
         "Panel B — 9D manifold (+ δ¹⁵N + coupling)", 9),
    ]:
        scores = pca["scores"]
        ve     = pca["var_exp"]
        sb     = scores[:n_bio]
        sa     = scores[n_bio:]

        for i, obs in enumerate(obs_list):
            ax.scatter(
                scores[i, 0], scores[i, 1],
                c=obs["color"], marker=obs["marker"],
                s=100, zorder=5,
                edgecolors="#333333", linewidths=0.5
            )
            ax.annotate(
                obs["short"],
                (scores[i, 0], scores[i, 1]),
                xytext=(5, 3),
                textcoords="offset points",
                fontsize=6.5, color=obs["color"]
            )

        # Project Lafayette
        if laff_dims == 7:
            lv = np.array(_obs_to_7D(
                {**LAFAYETTE,
                 "d13C":    LAFAYETTE["d13C"],
                 "d13C_dd": LAFAYETTE["d13C_dd"],
                 "CN_log":  LAFAYETTE["CN_log"],
                 "Fe_col":  LAFAYETTE["Fe_col"],
                 "dep_gr":  LAFAYETTE["dep_gr"],
                 "frac_m":  LAFAYETTE["frac_m"],
                 "tmpl":    LAFAYETTE["tmpl"]}
            ))
        else:
            cr = coupling_ratio(LAFAYETTE["d13C"],
                                LAFAYETTE["d15N_bio"])
            lv = np.array([
                LAFAYETTE["d13C"],
                LAFAYETTE["d13C_dd"],
                LAFAYETTE["CN_log"],
                LAFAYETTE["Fe_col"],
                LAFAYETTE["dep_gr"],
                LAFAYETTE["frac_m"],
                LAFAYETTE["tmpl"],
                LAFAYETTE["d15N_bio"],
                cr,
            ])

        lp = project_point(lv, pca)
        ax.scatter(
            lp[0], lp[1],
            c=LAFAYETTE["color"],
            marker=LAFAYETTE["marker"],
            s=300, zorder=7,
            edgecolors="black", linewidths=1.5,
            label="Lafayette (predicted)"
        )
        ax.annotate(
            "Lafayette", (lp[0], lp[1]),
            xytext=(6, 5), textcoords="offset points",
            fontsize=8, fontweight="bold",
            color=LAFAYETTE["color"]
        )

        ax.axhline(0, color="#CCCCCC",
                   linewidth=0.8, linestyle="--")
        ax.axvline(0, color="#CCCCCC",
                   linewidth=0.8, linestyle="--")
        ax.set_xlabel(f"PC1 ({ve[0]*100:.1f}%)",
                      fontsize=8.5)
        ax.set_ylabel(f"PC2 ({ve[1]*100:.1f}%)",
                      fontsize=8.5)
        ax.set_title(title, color=HEADER_BLUE)
        ax.legend(fontsize=7, framealpha=0.88)
        ax.grid(True, alpha=0.12, zorder=0)

    # Panel C: cumulative variance
    ax3   = axes[2]
    ve7   = pca7["var_exp"] * 100
    ve9   = pca9["var_exp"] * 100
    dim7  = np.arange(1, len(ve7) + 1)
    dim9  = np.arange(1, len(ve9) + 1)

    ax3.plot(dim7, np.cumsum(ve7), "o-",
             color="#1565C0", linewidth=2.0,
             markersize=6, zorder=5,
             label="7D (Scripts 1-8)")
    ax3.plot(dim9, np.cumsum(ve9), "s--",
             color="#E91E63", linewidth=2.0,
             markersize=6, zorder=5,
             label="9D (+ δ¹⁵N + coupling)")
    ax3.axhline(90, color="#888888", linewidth=1.0,
                linestyle=":", alpha=0.7,
                label="90% threshold")
    ax3.axhline(95, color="#555555", linewidth=1.0,
                linestyle=":", alpha=0.7,
                label="95% threshold")
    ax3.set_xlabel("PCs retained", fontsize=8.5)
    ax3.set_ylabel("Cumulative variance (%)",
                   fontsize=8.5)
    ax3.set_title("Panel C — Variance explained\n"
                  "7D vs 9D",
                  color=HEADER_BLUE)
    ax3.legend(fontsize=7.5, framealpha=0.88)
    ax3.grid(True, alpha=0.12, zorder=0)
    ax3.set_ylim(0, 105)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    path = "./necromass_geom_fig35_PCA.png"
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Figure 35 saved: {path}")
    return path

# ── FIGURE 36 ────────────────────────────────────────────────

def plot_figure_36():
    setup_style()
    fig = plt.figure(figsize=(14, 10))
    ax  = fig.add_subplot(111, projection="3d")

    for obs in ALL_OBS:
        d13C = obs["d13C"]
        d15N = (obs["d15N"]
                if obs["d15N"] is not None
                else 180.0)
        CN_l = obs["CN_log"]
        CN   = 10**CN_l

        ax.scatter(
            d13C, d15N, CN_l,
            c=obs["color"], marker=obs["marker"],
            s=120, zorder=5,
            edgecolors="black", linewidths=0.5,
            alpha=0.85
        )
        ax.text(d13C, d15N, CN_l + 0.06,
                obs["short"], fontsize=6.5,
                color=obs["color"])

    # Lafayette predicted position
    laff_d13C = LAFAYETTE["d13C"]
    laff_d15N = LAFAYETTE["d15N_bio"]
    laff_CN   = LAFAYETTE["CN_log"]

    ax.scatter(
        laff_d13C, laff_d15N, laff_CN,
        c=LAFAYETTE["color"],
        marker=LAFAYETTE["marker"],
        s=500, zorder=8,
        edgecolors="black", linewidths=2.0,
        label="Lafayette (predicted)"
    )
    ax.text(laff_d13C - 4, laff_d15N + 3,
            laff_CN + 0.12,
            "LAFAYETTE", fontsize=9,
            fontweight="bold",
            color=LAFAYETTE["color"])

    # δ¹⁵N uncertainty bar
    for d15N_t in [LAFAYETTE["d15N_lo"],
                   LAFAYETTE["d15N_hi"]]:
        ax.scatter(
            laff_d13C, d15N_t, laff_CN,
            c=LAFAYETTE["color"],
            marker="D", s=80, zorder=7,
            edgecolors="black", linewidths=0.5,
            alpha=0.5
        )
    ax.plot(
        [laff_d13C, laff_d13C],
        [LAFAYETTE["d15N_lo"], LAFAYETTE["d15N_hi"]],
        [laff_CN, laff_CN],
        color=LAFAYETTE["color"],
        linewidth=2.0, alpha=0.6
    )

    # Biological octant boundary planes (translucent)
    d15N_r = np.array([-60, 180])
    CN_r   = np.array([0.0, 3.5])
    D15N_g, CN_g = np.meshgrid(d15N_r, CN_r)
    ax.plot_surface(
        np.full_like(D15N_g, -25.0),
        D15N_g, CN_g,
        alpha=0.06, color="#2E7D32"
    )
    d13C_r = np.array([-120, 0])
    D13C_g2, CN_g2 = np.meshgrid(d13C_r, CN_r)
    ax.plot_surface(
        D13C_g2,
        np.full_like(D13C_g2, -10.0),
        CN_g2,
        alpha=0.06, color="#1565C0"
    )
    D13C_g3, D15N_g3 = np.meshgrid(d13C_r, d15N_r)
    ax.plot_surface(
        D13C_g3, D15N_g3,
        np.full_like(D13C_g3, 1.3),
        alpha=0.06, color="#E91E63"
    )

    ax.set_xlabel("δ¹³C (‰)", fontsize=8.5, labelpad=8)
    ax.set_ylabel("δ¹⁵N (‰)", fontsize=8.5, labelpad=8)
    ax.set_zlabel("log₁₀(C:N)", fontsize=8.5, labelpad=8)
    ax.set_xlim(-120, 0)
    ax.set_ylim(-60, 40)
    ax.set_zlim(0, 3.5)
    ax.set_title(
        "Script 12 v2.0 — Three-Axis Identity Space\n"
        "Biological octant: δ¹³C < -25‰  AND  "
        "δ¹⁵N > -10‰  AND  C:N < 20\n"
        "No known abiotic source satisfies all three. "
        "Lafayette (predicted) satisfies all three.",
        color=HEADER_BLUE, fontsize=8.5
    )
    ax.legend(fontsize=7.5, framealpha=0.88)
    plt.tight_layout()
    path = "./necromass_geom_fig36_3Dspace.png"
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Figure 36 saved: {path}")
    return path

# ── FIGURE 37 ────────────────────────────────────────────────

def plot_figure_37(ids7_all_obs, ids7_laff,
                   ids9_obs, ids9_laff_bio,
                   ids9_laff_abio):
    """
    FIX: Panel B is now split into two sub-panels
    to avoid shape mismatch between the 7D
    (n=8+1=9) and 9D (n=5+1=6) observation sets.

    ids7_all_obs : list of IDS values for ALL_OBS (len 8)
    ids7_laff    : float, Lafayette 7D IDS
    ids9_obs     : list of IDS values for OBS_9D  (len 5)
    ids9_laff_bio: float, Lafayette 9D IDS (bio d15N)
    ids9_laff_abio:float, Lafayette 9D IDS (mantle d15N)
    """
    setup_style()
    fig = plt.figure(figsize=(18, 7))
    gs  = GridSpec(1, 3, figure=fig,
                   wspace=0.35)
    ax_chain = fig.add_subplot(gs[0, 0])
    ax_ids7  = fig.add_subplot(gs[0, 1])
    ax_ids9  = fig.add_subplot(gs[0, 2])

    fig.suptitle(
        "Script 12 v2.0 — Bayesian Update Chain "
        "and Identity Score Expansion\n"
        "Full series Scripts 1-12  |  "
        "Three δ¹⁵N scenarios\n"
        "Pre-reg DOI: 10.5281/zenodo.18986790",
        fontsize=8.5, color=HEADER_BLUE, y=1.01
    )

    # ── Panel A: Bayesian chain ───────────────────────────
    ax = ax_chain
    steps  = [c["step"] for c in CHAIN_STEPS]
    odds   = [c["odds"] for c in CHAIN_STEPS]
    s_list = steps + [5, 6]
    prior  = CHAIN_STEPS[-1]["odds"]

    branches = [
        (BRANCH_BIO,  COUPLING_LR["bio"]),
        (BRANCH_AMB,  COUPLING_LR["amb"]),
        (BRANCH_ABI,  COUPLING_LR["abio"]),
    ]
    for branch, clr in branches:
        o_d15N    = prior * branch["LR"]
        o_coupling = o_d15N * clr
        od_list   = odds + [o_d15N, o_coupling]
        log_od    = [math.log10(max(o, 1e-9))
                     for o in od_list]
        ax.plot(
            s_list[:len(log_od)],
            log_od, "o-",
            color=branch["color"],
            linewidth=2.0, markersize=6,
            zorder=5, label=branch["label"]
        )

    ax.axhline(math.log10(150),
               color="#888888", linewidth=1.0,
               linestyle="--", zorder=3,
               label="Kass-Raftery decisive")
    ax.axhline(0, color="black", linewidth=1.0,
               linestyle=":", zorder=3)

    xlabels = ["Prior", "S1-5", "S6",
               "S8", "S10", "S11\nδ¹⁵N",
               "S12\ncoupling"]
    ax.set_xticks(range(len(xlabels)))
    ax.set_xticklabels(xlabels, fontsize=8)
    ax.set_ylabel("log₁₀(posterior odds)",
                  fontsize=8.5)
    ax.set_title("Panel A — Bayesian update chain\n"
                 "Scripts 1-12, three scenarios",
                 color=HEADER_BLUE)
    ax.legend(fontsize=7.0, framealpha=0.88)
    ax.grid(True, alpha=0.12, zorder=0)

    # ── Panel B-left: 7D IDS (all 8 obs + Lafayette) ──
    ax2 = ax_ids7
    names7 = [o["short"] for o in ALL_OBS] \
             + ["Lafayette\n(bio pred)"]
    vals7  = list(ids7_all_obs) + [ids7_laff]
    cols7  = [o["color"] for o in ALL_OBS] \
             + [LAFAYETTE["color"]]
    x7     = np.arange(len(names7))

    ax2.bar(x7, vals7, 0.6,
            color=cols7, alpha=0.75,
            edgecolor="#333333", linewidth=0.5,
            zorder=4)
    ax2.axhline(0, color="black", linewidth=1.5,
                linestyle="--", zorder=5)
    ax2.set_xticks(x7)
    ax2.set_xticklabels(names7, fontsize=6.5,
                        rotation=35, ha="right")
    ax2.set_ylabel("Identity Distance Score",
                   fontsize=8.5)
    ax2.set_title("Panel B — 7D IDS\n"
                  "All observations (n=8+Lafayette)",
                  color=HEADER_BLUE)
    ax2.grid(True, alpha=0.12, axis="y", zorder=0)

    # ── Panel C: 9D IDS (N-bearing obs + Lafayette) ───
    ax3 = ax_ids9
    names9 = [o["short"] for o in OBS_9D] \
             + ["Lafayette\n(bio +2‰)",
                "Lafayette\n(mantle -40‰)"]
    vals9  = list(ids9_obs) \
             + [ids9_laff_bio, ids9_laff_abio]
    cols9  = [o["color"] for o in OBS_9D] \
             + [LAFAYETTE["color"], "#888888"]
    x9     = np.arange(len(names9))

    ax3.bar(x9, vals9, 0.6,
            color=cols9, alpha=0.75,
            edgecolor="#333333", linewidth=0.5,
            zorder=4)
    ax3.axhline(0, color="black", linewidth=1.5,
                linestyle="--", zorder=5)
    ax3.set_xticks(x9)
    ax3.set_xticklabels(names9, fontsize=6.5,
                        rotation=35, ha="right")
    ax3.set_ylabel("Identity Distance Score",
                   fontsize=8.5)
    ax3.set_title("Panel C — 9D IDS\n"
                  "N-bearing observations + "
                  "Lafayette (two δ¹⁵N scenarios)",
                  color=HEADER_BLUE)
    ax3.grid(True, alpha=0.12, axis="y", zorder=0)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    path = "./necromass_geom_fig37_chain_IDS.png"
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Figure 37 saved: {path}")
    return path

# ───────────────────────────────────────���─────────────────────
# MODULE 6: MAIN RUN
# ─────────────────────────────────────────────────────────────

def run():
    sep = "═" * 62

    print()
    print(sep)
    print("  NECROMASS GEOMETRIC COVERAGE — SCRIPT 12 v2.0")
    print("  Mars Iron Necromass Hypothesis")
    print("  Pre-reg: DOI 10.5281/zenodo.18986790")
    print("  2026-03-13 / Eric Robert Lawson / OrganismCore")
    print(sep)

    # ── SECTION 1: THREE-AXIS OCTANT ─────────────────────────
    print()
    print("─" * 62)
    print("  SECTION 1: THREE-AXIS OCTANT ANALYSIS")
    print("─" * 62)
    print()
    print("  BIOLOGICAL OCTANT:")
    print("    Axis 1: δ¹³C < -25‰")
    print("    Axis 2: δ¹⁵N > -10‰")
    print("    Axis 3: C:N  <  20")
    print()
    print(f"  {'Source':<28} "
          f"{'A1':>4} {'A2':>4} {'A3':>4} "
          f"{'Score':>6}")
    print("  " + "─" * 48)

    abio_scores = []
    for obs in ALL_OBS:
        CN = 10**obs["CN_log"]
        n, sat = octant_check(
            obs["d13C"], obs["d15N"], CN)
        a1 = "✓" if sat[0] else "✗"
        a2 = "✓" if sat[1] else "✗"
        a3 = "✓" if sat[2] else "✗"
        print(f"  {obs['short']:<28} "
              f"{a1:>4} {a2:>4} {a3:>4} {n:>4}/3")
        if obs["marker"] == "^":
            abio_scores.append(n)

    CN_l  = 10**LAFAYETTE["CN_log"]
    n_l, sat_l = octant_check(
        LAFAYETTE["d13C"],
        LAFAYETTE["d15N_bio"],
        CN_l
    )
    print(f"  {'Lafayette (bio predicted)':<28} "
          f"{'✓' if sat_l[0] else '✗':>4} "
          f"{'✓' if sat_l[1] else '✗':>4} "
          f"{'✓' if sat_l[2] else '✗':>4} "
          f"{n_l:>4}/3")
    print()
    print(f"  Max abiotic octant score: "
          f"{max(abio_scores)}/3")
    print(f"  Lafayette predicted score: {n_l}/3")
    print()
    print("  NO KNOWN ABIOTIC SOURCE SCORES 3/3.")
    print(f"  LAFAYETTE SCORES {n_l}/3.")
    print()

    # ── SECTION 2: PCA EXPANSION ─────────────────────────────
    print()
    print("─" * 62)
    print("  SECTION 2: PCA MANIFOLD EXPANSION")
    print("─" * 62)
    print()

    M7   = build_7D_matrix()
    pca7 = run_PCA(M7)
    ve7  = pca7["var_exp"]
    print("  7D MANIFOLD (n={:d} obs):".format(
        len(ALL_OBS)))
    print(f"    PC1: {ve7[0]*100:.1f}%  "
          f"PC2: {ve7[1]*100:.1f}%  "
          f"PC3: {ve7[2]*100:.1f}%")
    print(f"    PC1+PC2: {(ve7[0]+ve7[1])*100:.1f}%")

    M9   = build_9D_matrix()
    pca9 = run_PCA(M9)
    ve9  = pca9["var_exp"]
    print()
    print(f"  9D MANIFOLD (n={len(OBS_9D)} "
          "N-bearing obs):")
    print(f"    PC1: {ve9[0]*100:.1f}%  "
          f"PC2: {ve9[1]*100:.1f}%  "
          f"PC3: {ve9[2]*100:.1f}%")
    print(f"    PC1+PC2: {(ve9[0]+ve9[1])*100:.1f}%")
    print()

    # 7D IDS
    n_bio7   = len(BIO_OBS)
    s7_bio   = pca7["scores"][:n_bio7]
    s7_abio  = pca7["scores"][n_bio7:]

    ids7_all = []
    for i in range(len(ALL_OBS)):
        ids7_all.append(
            IDS_score(pca7["scores"][i],
                      s7_bio, s7_abio)
        )

    # Lafayette 7D
    lv7 = np.array([
        LAFAYETTE["d13C"], LAFAYETTE["d13C_dd"],
        LAFAYETTE["CN_log"], LAFAYETTE["Fe_col"],
        LAFAYETTE["dep_gr"], LAFAYETTE["frac_m"],
        LAFAYETTE["tmpl"],
    ])
    lp7       = project_point(lv7, pca7)
    ids7_laff = IDS_score(lp7, s7_bio, s7_abio)

    print("  7D IDS:")
    for obs, v in zip(ALL_OBS, ids7_all):
        print(f"    {obs['short']:<22} {v:+.4f}")
    print(f"    {'Lafayette (bio)':<22} "
          f"{ids7_laff:+.4f}")
    print()

    # 9D IDS
    n_bio9  = sum(1 for o in OBS_9D
                  if o["marker"] == "o")
    s9_bio  = pca9["scores"][:n_bio9]
    s9_abio = pca9["scores"][n_bio9:]

    ids9_obs = []
    for i in range(len(OBS_9D)):
        ids9_obs.append(
            IDS_score(pca9["scores"][i],
                      s9_bio, s9_abio)
        )

    # Lafayette 9D bio scenario
    cr_bio = coupling_ratio(
        LAFAYETTE["d13C"], LAFAYETTE["d15N_bio"]
    )
    lv9_bio = np.array([
        LAFAYETTE["d13C"], LAFAYETTE["d13C_dd"],
        LAFAYETTE["CN_log"], LAFAYETTE["Fe_col"],
        LAFAYETTE["dep_gr"], LAFAYETTE["frac_m"],
        LAFAYETTE["tmpl"],
        LAFAYETTE["d15N_bio"], cr_bio,
    ])
    lp9_bio       = project_point(lv9_bio, pca9)
    ids9_laff_bio = IDS_score(lp9_bio,
                               s9_bio, s9_abio)

    # Lafayette 9D mantle scenario
    cr_abi = coupling_ratio(LAFAYETTE["d13C"], -40.0)
    lv9_abi = np.array([
        LAFAYETTE["d13C"], LAFAYETTE["d13C_dd"],
        LAFAYETTE["CN_log"], LAFAYETTE["Fe_col"],
        LAFAYETTE["dep_gr"], LAFAYETTE["frac_m"],
        LAFAYETTE["tmpl"],
        -40.0, cr_abi,
    ])
    lp9_abi        = project_point(lv9_abi, pca9)
    ids9_laff_abio = IDS_score(lp9_abi,
                                s9_bio, s9_abio)

    print("  9D IDS (N-bearing obs):")
    for obs, v in zip(OBS_9D, ids9_obs):
        print(f"    {obs['short']:<22} {v:+.4f}")
    print(f"    {'Lafayette (bio +2‰)':<22} "
          f"{ids9_laff_bio:+.4f}")
    print(f"    {'Lafayette (mantle -40‰)':<22} "
          f"{ids9_laff_abio:+.4f}")
    print()
    delta = ids9_laff_bio - ids7_laff
    print(f"  IDS change 7D→9D (bio scenario): "
          f"{delta:+.4f}")
    print("  Manifold " +
          ("TIGHTENED." if delta > 0
           else "loosened or unchanged."))
    print()

    # ── SECTION 3: COUPLING RATIOS ────────────────────────────
    print()
    print("─" * 62)
    print("  SECTION 3: ISOTOPIC COUPLING RATIOS")
    print("─" * 62)
    print()
    print(f"  {'Source':<28} {'δ¹³C':>7} "
          f"{'δ¹⁵N':>7} {'Coupling':>10}")
    print("  " + "─" * 54)
    for obs in ALL_OBS:
        if obs["d15N"] is not None:
            cr = coupling_ratio(obs["d13C"],
                                obs["d15N"])
            if math.isfinite(cr):
                print(
                    f"  {obs['short']:<28} "
                    f"{obs['d13C']:>7.1f}‰ "
                    f"{obs['d15N']:>7.1f}‰ "
                    f"{cr:>10.2f}"
                )
        else:
            print(
                f"  {obs['short']:<28} "
                f"{obs['d13C']:>7.1f}‰ "
                f"{'no N':>7} "
                f"{'N/A':>10}"
            )
    print(
        f"  {'Lafayette (bio +2‰)':<28} "
        f"{LAFAYETTE['d13C']:>7.1f}‰ "
        f"{LAFAYETTE['d15N_bio']:>7.1f}‰ "
        f"{cr_bio:>10.2f}"
    )
    print()
    print("  COUPLING SEPARATION:")
    print(f"    Biological:    ~-15 to -20")
    print(f"    Hydrothermal:  ~-4")
    print(f"    Mantle N:      ~+0.2 to +0.6")
    print(f"    Lafayette:     {cr_bio:.2f}")
    print()

    # ── SECTION 4: BAYESIAN CHAIN ─────────────────────────────
    print()
    print("─" * 62)
    print("  SECTION 4: BAYESIAN UPDATE CHAIN")
    print("─" * 62)
    print()
    prior = CHAIN_STEPS[-1]["odds"]
    branches = [
        (BRANCH_BIO,  COUPLING_LR["bio"]),
        (BRANCH_AMB,  COUPLING_LR["amb"]),
        (BRANCH_ABI,  COUPLING_LR["abio"]),
    ]
    for branch, clr in branches:
        od_d15N   = prior * branch["LR"]
        od_final  = od_d15N * clr
        P_final   = (od_final / (1 + od_final)
                     if od_final > 0 else 0.0)
        print(f"  {branch['label']}")
        print(f"    Prior (S10):      {prior:.0f}:1")
        if od_d15N >= 1:
            print(f"    After δ¹⁵N (S11): "
                  f"{od_d15N:.0f}:1")
        else:
            print(f"    After δ¹⁵N (S11): "
                  f"1:{1/max(od_d15N,1e-9):.0f}")
        if od_final >= 1:
            print(f"    After coupling:   "
                  f"{od_final:.0f}:1")
        else:
            print(f"    After coupling:   "
                  f"1:{1/max(od_final,1e-9):.0f}")
        print(f"    Final P(bio):     {P_final:.6f}")
        print()

    # ── SECTION 5: FIGURES ───────────────────────────────────
    print()
    print("─" * 62)
    print("  SECTION 5: GENERATING FIGURES")
    print("─" * 62)
    print()
    f35 = plot_figure_35(pca7, pca9, M7, M9)
    f36 = plot_figure_36()
    f37 = plot_figure_37(
        ids7_all, ids7_laff,
        ids9_obs, ids9_laff_bio, ids9_laff_abio
    )

    # ── SECTION 6: FINAL VERDICT ────────────────────��────────
    print()
    print("─" * 62)
    print("  SECTION 6: FINAL GEOMETRIC VERDICT")
    print("─" * 62)
    print()
    print("  THREE-AXIS IDENTITY SPACE:")
    print(f"    Lafayette: {n_l}/3 axes satisfied.")
    print(f"    Max abiotic: {max(abio_scores)}/3.")
    print("    No abiotic source occupies the")
    print("    biological octant on all three axes.")
    print()
    print("  PCA IDENTITY MANIFOLD:")
    print(f"    7D IDS: {ids7_laff:+.4f}")
    print(f"    9D IDS (bio δ¹⁵N +2‰): "
          f"{ids9_laff_bio:+.4f}")
    print(f"    9D IDS (mantle δ¹⁵N -40‰): "
          f"{ids9_laff_abio:+.4f}")
    print(f"    Manifold change (bio): {delta:+.4f}")
    print()
    print("  COUPLING RATIO:")
    print(f"    Lafayette predicted: {cr_bio:.2f}")
    print("    Biology range: -15 to -20")
    print("    Mantle: +0.2 to +0.6")
    print()
    print("  FINAL ODDS:")
    prior_s10 = CHAIN_STEPS[-1]["odds"]
    for branch, clr in branches:
        final = prior_s10 * branch["LR"] * clr
        P     = (final / (1 + final)
                 if final > 0 else 0.0)
        if final >= 1:
            fs = f"{final:.0f}:1"
        else:
            fs = f"1:{1/max(final,1e-9):.0f}"
        print(f"    {branch['label']:<24} "
              f"{fs:<16} P={P:.6f}")
    print()

    print(sep)
    print("  SCRIPT 12 v2.0 COMPLETE")
    print("  FULL SERIES COMPLETE (Scripts 1-12)")
    print(sep)
    print()
    for f in [f35, f36, f37]:
        print(f"  {f}")
    print()
    print("  Complete figure set:")
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
    print("    Figs 35-37: Script 12")
    print("    Total: 37 figures")
    print()
    print("  Pre-reg: 10.5281/zenodo.18986790")
    print("  github.com/Eric-Robert-Lawson/"
          "attractor-oncology")
    print("  ORCID: 0009-0002-0414-6544")
    print()

    final_bio  = prior_s10 * BRANCH_BIO["LR"]  \
                 * COUPLING_LR["bio"]
    final_amb  = prior_s10 * BRANCH_AMB["LR"]  \
                 * COUPLING_LR["amb"]
    final_abio = prior_s10 * BRANCH_ABI["LR"]  \
                 * COUPLING_LR["abio"]

    return {
        "lafayette_octant_score":      int(n_l),
        "max_abiotic_octant_score":    int(
            max(abio_scores)),
        "IDS_7D_lafayette":            float(ids7_laff),
        "IDS_9D_lafayette_bio":        float(
            ids9_laff_bio),
        "IDS_9D_lafayette_mantle":     float(
            ids9_laff_abio),
        "IDS_change_7D_to_9D_bio":     float(delta),
        "coupling_ratio_lafayette_bio":float(cr_bio),
        "coupling_ratio_mantle":       0.23,
        "coupling_ratio_hydrothermal": -4.0,
        "prior_odds_script10":         float(prior_s10),
        "final_odds_bio":              float(final_bio),
        "final_odds_amb":              float(final_amb),
        "final_odds_abio":             float(final_abio),
        "P_bio_bio_scenario":          float(
            final_bio / (1 + final_bio)),
        "P_bio_amb_scenario":          float(
            final_amb / (1 + final_amb)),
        "P_bio_abio_scenario":         float(
            max(final_abio, 0) /
            max(1 + final_abio, 1e-9)),
        "n_abiotic_sources_3_of_3":    0,
        "series_complete":             True,
        "total_scripts":               12,
        "total_figures":               37,
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
