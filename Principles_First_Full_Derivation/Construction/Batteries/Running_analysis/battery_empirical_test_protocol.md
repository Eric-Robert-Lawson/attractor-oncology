# OC-BATTERY-EMPIRICAL-TEST-PROTOCOL.md
## Empirical Test Protocol: SCE as Unified Causal Variable
## Using Publicly Available Data and Python Analysis
## OrganismCore — Eric Robert Lawson
## Date: 2026-04-01

---

## WHAT WE ARE DOING

We are testing whether solvation
configuration entropy (SCE) predicts
lithium metal battery cycling
performance better than the variables
the field currently uses.

We are doing this using:
— Publicly available published data
— Python scripts running on a MacBook Air
— No new experiments
— No new synthesis
— No laboratory access required

If the test holds, we have the first
empirical confirmation of the unified
causal variable from existing data.

That is a publishable finding.

---

## THE CORE PREDICTION BEING TESTED

```
PREDICTION 5 (cheapest, first to test):

Across a dataset of characterised
electrolyte systems with known cycling
performance, the variance of the
solvation shell geometry distribution
— quantified as Shannon entropy (SCE)
— will correlate with cycling
performance more strongly than the
primary variable each study used to
explain its own results.

In plain terms:
The spread of bubble shapes the lithium
ion arrives in predicts how long the
battery lasts better than any single
molecular property currently measured.
```

---

## THE DATA SOURCES

These are publicly available.
No subscription required for the
key data we need.
All coordination distribution data
is reported in paper text, figures,
or supplementary tables.

### PRIMARY TARGET — USE THIS FIRST

```
PAPER:
"Statistical insights into solvation
sheaths of lithium-ions"

Journal: Journal of Power Sources
Year: 2025
DOI: 10.1016/j.jpowsour.2025.237088

WHY THIS PAPER:
This paper literally computed
coordination number distributions
for lithium ions in DEC/EC electrolytes.
They tracked the population fraction
of every coordination geometry
configuration across the Li⁺ population.

They reported:
— Population fractions for each
  coordination configuration
  (nDEC, nEC, nPF6⁻) as tuples
— Most common: five-coordinate
  (1DEC/4EC/0PF6⁻) and (2DEC/3EC/0PF6⁻)
— Most stable: four-coordinate CIP
  (1DEC/1EC/2PF6⁻)
— EC drives >55% of exchange events

This data IS the coordination geometry
distribution. We compute SCE directly
from it. One afternoon.

WHAT TO ACCESS:
The paper itself (try Sci-Hub,
unpaywall.org, or email the authors —
most will share for research).
The supplementary information.
Look for tables of coordination
configuration populations.
```

### SECONDARY SOURCES — USE AFTER PRIMARY

```
SOURCE 2:
arxiv.org/abs/2501.11932
"Computational Study of Li+ Solvation
Structures in Fluorinated Ether,
Non-Fluorinated Ether, and Organic
Carbonate-Based Electrolytes"
FREE on arXiv. No paywall.
Contains coordination number
distribution data for multiple
electrolyte types with different
performance characteristics.

SOURCE 3:
NOMAD Repository
nomad-lab.eu
Free, publicly accessible.
Search: "lithium electrolyte solvation
molecular dynamics"
Many groups deposit raw MD data here.
If raw trajectory files available,
we can extract distributions ourselves.

SOURCE 4:
Zenodo
zenodo.org
Search: "lithium solvation MD trajectory"
Authors increasingly deposit simulation
data here open access.

SOURCE 5:
Materials Project
materialsproject.org
Free account. Some electrolyte
solvation data accessible.

SOURCE 6:
The Hunan University / Joule 2025 paper
"Solvation-Configurational Entropy
Governs Interfacial Kinetics in
Low-Temperature Batteries"
This paper computed Ssc values for
their electrolytes. If we can access
their values, we can correlate against
their performance data directly.
Their data tests the band hypothesis.
```

---

## THE SCE CALCULATION

This is the computation we are running.
It is trivial. The bottleneck is data
extraction, not computation.

### What SCE Is

```
SCE = Shannon entropy of the Li⁺
coordination geometry distribution.

SCE = -Σ p(gᵢ) × log(p(gᵢ))

where:
gᵢ = a distinct solvation geometry
     configuration (e.g. 1DEC/1EC/2PF6⁻)
p(gᵢ) = the fraction of Li⁺ ions
         in that configuration

LOW SCE = tight distribution =
  few dominant configurations =
  coherent navigator = good performance

HIGH SCE = broad distribution =
  many competing configurations =
  incoherent navigator = dendrites
```

### The Python Script

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

# ============================================
# STEP 1: INPUT YOUR COORDINATION DATA
# ============================================
# Extract population fractions from the paper.
# Each entry is one electrolyte condition.
# p values must sum to 1.0 for each electrolyte.

# Example format — replace with real paper data:
# Each row = one electrolyte system
# Each column = one coordination geometry
# configuration population fraction

electrolyte_labels = [
    "EC/DEC standard",
    "EC/DEC high concentration",
    "LHCE fluorinated",
    "LiFSI/DME",
    # add more as extracted from papers
]

# Population fractions for each geometry
# configuration in each electrolyte.
# Rows = electrolytes, Columns = configurations
# These numbers come from the paper tables.
coordination_distributions = [
    [0.45, 0.30, 0.15, 0.10],  # EC/DEC standard
    [0.60, 0.25, 0.10, 0.05],  # EC/DEC high conc
    [0.75, 0.15, 0.07, 0.03],  # LHCE fluorinated
    [0.85, 0.10, 0.03, 0.02],  # LiFSI/DME
    # replace with actual values from paper
]

# Known cycling performance for each electrolyte
# Use whatever metric the paper reports:
# — Cycle number at 80% capacity retention
# — Coulombic efficiency (%)
# — Capacity retention at N cycles (%)
cycling_performance = [
    150,   # EC/DEC standard (cycles at 80%)
    300,   # EC/DEC high conc
    450,   # LHCE fluorinated
    600,   # LiFSI/DME
    # replace with actual values from paper
]

# ============================================
# STEP 2: COMPUTE SCE FOR EACH ELECTROLYTE
# ============================================

def compute_SCE(p_array):
    """
    Compute Shannon entropy of a probability
    distribution. This is the Solvation
    Configuration Entropy (SCE).

    Parameters:
        p_array: list or array of population
                 fractions. Must sum to 1.0.

    Returns:
        SCE value (float). Lower = less variance
        = more coherent navigator geometry.
    """
    p = np.array(p_array)
    # Remove zeros to avoid log(0)
    p = p[p > 0]
    # Normalise to ensure sum = 1
    p = p / p.sum()
    return -np.sum(p * np.log(p))

SCE_values = []
for dist in coordination_distributions:
    sce = compute_SCE(dist)
    SCE_values.append(sce)
    print(f"SCE = {sce:.4f}")

# ============================================
# STEP 3: CORRELATION ANALYSIS
# ============================================

SCE_arr = np.array(SCE_values)
perf_arr = np.array(cycling_performance)

# Pearson correlation (linear)
r_pearson, p_pearson = pearsonr(SCE_arr, perf_arr)

# Spearman correlation (rank-based, more robust)
r_spearman, p_spearman = spearmanr(SCE_arr, perf_arr)

print("\n=== CORRELATION RESULTS ===")
print(f"Pearson r:  {r_pearson:.4f}  (p = {p_pearson:.4f})")
print(f"Spearman r: {r_spearman:.4f}  (p = {p_spearman:.4f})")
print(f"\nNote: Negative correlation expected.")
print(f"Lower SCE → better performance → r should be negative.")

# ============================================
# STEP 4: COMPARE AGAINST PAPER'S OWN
# PRIMARY VARIABLE
# ============================================
# The paper will report its own explanatory
# variable (e.g. average coordination number,
# average solvation energy, donor number etc.)
# Enter those values here and compare R².

paper_primary_variable = [
    4.8,   # e.g. average coordination number
    4.2,
    3.1,
    1.8,
    # replace with actual values from paper
]

pv_arr = np.array(paper_primary_variable)
r_pv, p_pv = pearsonr(pv_arr, perf_arr)

print("\n=== COMPARISON ===")
print(f"SCE vs performance:              R² = {r_pearson**2:.4f}")
print(f"Paper primary variable vs perf:  R² = {r_pv**2:.4f}")
print()
if r_pearson**2 > r_pv**2:
    print("SCE is MORE predictive than paper's primary variable.")
    print("PREDICTION 5 SUPPORTED.")
else:
    print("Paper's primary variable is more predictive.")
    print("PREDICTION 5 NOT SUPPORTED. Revise geometry.")

# ============================================
# STEP 5: VISUALISATION
# ============================================

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Plot 1: SCE vs cycling performance
axes[0].scatter(SCE_arr, perf_arr, color='navy',
                s=100, zorder=5)
for i, label in enumerate(electrolyte_labels):
    axes[0].annotate(label,
                     (SCE_arr[i], perf_arr[i]),
                     textcoords="offset points",
                     xytext=(5, 5), fontsize=8)

# Trend line
z = np.polyfit(SCE_arr, perf_arr, 1)
p_line = np.poly1d(z)
x_line = np.linspace(SCE_arr.min(), SCE_arr.max(), 100)
axes[0].plot(x_line, p_line(x_line), "r--", alpha=0.7)

axes[0].set_xlabel("SCE (Solvation Configuration Entropy)",
                   fontsize=11)
axes[0].set_ylabel("Cycling Performance", fontsize=11)
axes[0].set_title(
    f"SCE vs Performance\nPearson r = {r_pearson:.3f}, "
    f"R² = {r_pearson**2:.3f}",
    fontsize=11)
axes[0].invert_xaxis()  # Low SCE on right = better
axes[0].grid(True, alpha=0.3)

# Plot 2: Paper's primary variable vs performance
axes[1].scatter(pv_arr, perf_arr, color='darkgreen',
                s=100, zorder=5)
z2 = np.polyfit(pv_arr, perf_arr, 1)
p_line2 = np.poly1d(z2)
x_line2 = np.linspace(pv_arr.min(), pv_arr.max(), 100)
axes[1].plot(x_line2, p_line2(x_line2), "r--", alpha=0.7)

axes[1].set_xlabel("Paper's Primary Variable", fontsize=11)
axes[1].set_ylabel("Cycling Performance", fontsize=11)
axes[1].set_title(
    f"Primary Variable vs Performance\nPearson r = {r_pv:.3f}, "
    f"R² = {r_pv**2:.3f}",
    fontsize=11)
axes[1].grid(True, alpha=0.3)

plt.suptitle(
    "SCE vs Paper's Primary Variable: Predictive Power Comparison\n"
    "OC Battery Framework — Pre-Registration Test — 2026-04-01",
    fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig("SCE_vs_performance_test.png", dpi=150,
            bbox_inches='tight')
plt.show()
print("\nPlot saved: SCE_vs_performance_test.png")
```

---

## THE BAND HYPOTHESIS TEST

This tests whether the Joule 2025
finding (high Ssc = better low-T)
and the framework finding (low SCE =
better room-T) together define a band.

```python
import numpy as np
import matplotlib.pyplot as plt

# ============================================
# BAND HYPOTHESIS TEST
# ============================================
# We need electrolytes with known SCE values
# and performance data at TWO temperatures:
# room temperature AND low temperature.
#
# If the band exists:
# — Low SCE electrolytes perform best at room-T
# — Intermediate SCE electrolytes perform best
#   across both temperatures
# — High SCE electrolytes perform best at low-T
#   but worse at room-T
# — A performance peak exists at intermediate
#   SCE at low temperature

# Enter SCE values and performance at each temp
# Data from multiple papers if needed

SCE_values = [0.1, 0.3, 0.6, 0.9, 1.2, 1.5]
# ^ replace with computed values from data

room_temp_performance = [95, 88, 75, 60, 45, 30]
# ^ cycles at 80% capacity, room temperature

low_temp_performance = [20, 45, 70, 80, 75, 65]
# ^ capacity retention at -20°C or similar

SCE_arr = np.array(SCE_values)
rt_arr = np.array(room_temp_performance)
lt_arr = np.array(low_temp_performance)

fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(SCE_arr, rt_arr, 'b-o', linewidth=2,
        markersize=8, label='Room Temperature')
ax.plot(SCE_arr, lt_arr, 'r-s', linewidth=2,
        markersize=8, label='Low Temperature')

# Mark the band
# Band = region where both curves are above
# acceptable threshold
threshold = 60  # define acceptable performance
ax.axhline(y=threshold, color='gray',
           linestyle='--', alpha=0.5,
           label=f'Acceptable threshold ({threshold})')

# Shade the band region
band_mask = (rt_arr >= threshold) & (lt_arr >= threshold)
if band_mask.any():
    band_SCE = SCE_arr[band_mask]
    ax.axvspan(band_SCE.min(), band_SCE.max(),
               alpha=0.2, color='green',
               label='Optimal SCE Band')

ax.set_xlabel("SCE (Solvation Configuration Entropy)",
              fontsize=12)
ax.set_ylabel("Performance Metric", fontsize=12)
ax.set_title(
    "The Band Hypothesis\n"
    "Optimal SCE Range for All-Temperature Operation\n"
    "OC Framework — 2026-04-01",
    fontsize=12, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("band_hypothesis_test.png", dpi=150,
            bbox_inches='tight')
plt.show()
print("Plot saved: band_hypothesis_test.png")

# Find the band
print("\n=== BAND ANALYSIS ===")
for i, sce in enumerate(SCE_values):
    rt = room_temp_performance[i]
    lt = low_temp_performance[i]
    in_band = rt >= threshold and lt >= threshold
    print(f"SCE={sce:.1f}: RT={rt}, LT={lt}, "
          f"In band: {in_band}")
```

---

## DATA EXTRACTION GUIDE

This is how to get the numbers out of
the papers and into the Python scripts.

### What to look for in each paper

```
TARGET DATA TYPE 1:
Coordination number distribution tables.

What they look like in papers:
A table with columns like:
(nDEC, nEC, nPF6⁻) | Population %
(0, 4, 0)           | 12%
(1, 3, 0)           | 28%
(2, 2, 0)           | 31%
(1, 2, 1)           | 18%
(0, 3, 1)           | 11%

These percentages ARE the p(gᵢ)
values for the SCE calculation.
Divide by 100 to get fractions.
Plug directly into compute_SCE().

TARGET DATA TYPE 2:
Radial distribution function (RDF)
peak heights or coordination numbers.

What they look like:
"The average coordination number
of Li⁺ with DME is 2.3, with
FSI⁻ is 1.8"

This is less precise than Type 1
but can be used as a proxy.
If only average coordination numbers
are reported, use them to estimate
the distribution shape.

TARGET DATA TYPE 3:
Cluster population data.

What they look like:
"Solvent-separated ion pairs (SSIP):
42%, Contact ion pairs (CIP): 38%,
Aggregates: 20%"

Three states. Three population
fractions. Compute SCE from these.
Coarser than Type 1 but usable.

PERFORMANCE DATA to record alongside:
— Cycle number at 80% capacity (best)
— Coulombic efficiency (%)
— Capacity retention at fixed cycles
— Any quantitative performance metric
  the paper reports
```

### Step by step data extraction

```
STEP 1:
Get the paper.
Try in order:
1. unpaywall.org (paste DOI)
2. arxiv.org (many are deposited here)
3. sci-hub.se (paste DOI)
4. Email corresponding author
   "I am researching solvation
   structure and would appreciate
   access to your 2025 paper..."
   Most will respond within days.

STEP 2:
Search the paper for these terms:
"coordination number distribution"
"cluster population"
"solvation structure population"
"contact ion pair"
"aggregate"
"SSIP" "CIP" "AGG"
These signal the data we need.

STEP 3:
Find the supplementary information.
Most coordination distribution data
is in supplementary tables or figures.
Download the SI separately.

STEP 4:
If data is in a figure (bar chart,
pie chart), use WebPlotDigitizer:
automeris.io/WebPlotDigitizer
Free, browser-based.
Upload the figure image.
Click on each bar to extract values.

STEP 5:
Enter values into the Python script.
Run. Record result.
```

---

## THE FULL TEST SEQUENCE

Run these in order. Stop if
a test fails and revise before
proceeding.

### TEST 1 — The Afternoon Test
**Time: 2-4 hours**
**Data: JPS 2025 paper**
**What we are computing:**
SCE from their coordination
distribution tables.
Correlate against their reported
performance differences.

**Success condition:**
SCE correlates more strongly with
performance than their reported
primary variable.
R²(SCE) > R²(primary variable)

**What it confirms:**
Prediction 5 partially confirmed.
Unified causal variable demonstrated
from a single paper's existing data.

---

### TEST 2 — The Multi-Paper Test
**Time: 1-2 weeks**
**Data: arXiv 2501.11932 + 2 more**
**What we are computing:**
SCE across multiple electrolyte
types from the arXiv paper.
Build a dataset of 8-15 data points
across different electrolyte classes.

**Success condition:**
Monotonic relationship between SCE
and performance across all electrolyte
classes. The gradient scale confirmed
computationally.

**What it confirms:**
The monotonic gradient is real.
The structural invariant holds across
electrolyte classes.

---

### TEST 3 — The Band Test
**Time: 1-2 weeks**
**Data: Joule 2025 Hunan University
paper + room-temperature data for
the same electrolytes**
**What we are computing:**
SCE for each electrolyte in their
dataset. Plot against both their
low-temperature performance and the
room-temperature performance of the
same electrolytes from other papers.

**Success condition:**
Room-temperature performance peaks
at low SCE. Low-temperature performance
peaks at intermediate SCE. The band
is visible in the data.

**What it confirms:**
The band hypothesis is real.
The temperature-responsive SCE
electrolyte is the correct target.

---

### TEST 4 — The Derived Candidates Test
**Time: Weeks (computational)**
**Data: MD simulation of derived
candidate molecules**
**What we are computing:**
SCE for Target 4 (fluorinated
orthocarbonate esters) from MD
simulation of the molecule in mixture.

This test requires more than a
MacBook Air for the simulation itself
— but the simulation can be run on
Google Colab with GPU (free tier)
using GROMACS or LAMMPS.
The SCE calculation from the
trajectory output runs on MacBook Air.

**Success condition:**
Target 4 SCE is lower than HFTHP/BTFMD.
If confirmed computationally,
proceed to synthesis.

**What it confirms:**
The derived candidates are real.
The framework found molecules the
field has not explored.

---

## SETUP INSTRUCTIONS

### Environment Setup

```bash
# Install required Python packages
pip install numpy scipy matplotlib pandas

# Optional but useful
pip install mdanalysis  # for MD trajectory analysis
pip install scikit-learn  # for clustering if needed

# Verify installation
python -c "import numpy, scipy, matplotlib; print('Ready')"
```

### File Structure

```
OC-Battery-Test/
│
├── data/
│   ├── JPS_2025_coordination_data.csv
│   ├── arXiv_2501_11932_data.csv
│   ├── hunan_joule_2025_data.csv
│   └── [add papers as extracted]
│
├── scripts/
│   ├── compute_SCE.py
│   ├── correlation_analysis.py
│   ├── band_hypothesis_test.py
│   └── visualise_gradient.py
│
├── results/
│   ├── SCE_vs_performance_test.png
│   ├── band_hypothesis_test.png
│   └── gradient_scale_confirmed.png
│
└── OC-BATTERY-EMPIRICAL-TEST-PROTOCOL.md
```

### Data CSV Format

```csv
electrolyte,config_1,config_2,config_3,config_4,config_5,performance_metric,metric_type,source
EC_DEC_standard,0.28,0.31,0.18,0.12,0.11,150,cycles_80pct,JPS_2025
LiFSI_DME_3M,0.72,0.18,0.07,0.03,0.00,600,cycles_80pct,multiple
HFTHP_diluent,0.91,0.07,0.02,0.00,0.00,900,cycles_80pct,NatEnergy_2024
```

---

## WHAT SUCCESS LOOKS LIKE

If the tests hold, you will have:

```
A scatter plot showing:
X axis: SCE (low on right, high on left)
Y axis: Cycling performance (cycles, CE%)

Points from 10-20 different published
electrolyte systems spanning the full
range from standard carbonate to HFTHP.

A clear monotonic relationship:
Low SCE → high performance
High SCE → low performance

R²(SCE) > R²(any other variable)
across the dataset.

And a second plot showing:
The band — where room-temperature
performance and low-temperature
performance both sit above threshold
only in an intermediate SCE window.
```

That is the result. Two plots.
Built from published data.
On a MacBook Air.
In one to two weeks.

If those plots look the way the
geometry predicts they will look,
you bring them to the expert.
You do not bring a theory.
You bring a result.

---

## WHAT TO DO IF IT DOES NOT HOLD

If SCE does not predict performance
better than the paper's primary variable:

```
DIAGNOSTIC QUESTIONS:

1. Was the data precise enough?
   Coordination NUMBER distributions
   (coarser) vs coordination GEOMETRY
   distributions (more precise).
   If only numbers were available,
   find a paper with geometry data.

2. Was the performance metric comparable?
   Mixing Coulombic efficiency with
   cycle life across papers adds noise.
   Restrict to one metric first.

3. Is the dataset too small?
   Fewer than 6 data points produces
   unreliable correlations.
   Expand the dataset before concluding.

4. Is the relationship nonlinear?
   Plot SCE vs performance as a scatter.
   If there is a curve rather than a
   line, use Spearman not Pearson.
   A nonlinear monotonic relationship
   still confirms the framework.

5. Is there a confounding variable?
   If two electrolytes have similar SCE
   but very different performance, there
   is a second variable operating.
   Record this. The geometry is then
   incomplete, not wrong. Add the
   confound to the model.

A failed test is not a dead end.
It is a geometric signal.
It tells you what the model is missing.
Record the result. Revise. Retest.
```

---

## THE STATEMENT WHEN IT HOLDS

When the tests hold, this is
what you bring to the expert
and eventually to a journal:

```
We tested whether solvation
configuration entropy (SCE) —
the Shannon entropy of the Li⁺
coordination geometry distribution
extracted from published molecular
dynamics datasets — predicts
lithium metal battery cycling
performance more precisely than
the primary variables reported by
each study.

Across [N] electrolyte systems
from [X] independent published
studies spanning standard carbonate
electrolytes to fluorinated cyclic
ether diluent systems, SCE showed
a monotonic negative correlation
with cycling performance
(Pearson r = [value], R² = [value])
that exceeded the predictive power
of [the field's primary variable]
(R² = [value]).

Additionally, consistent with
Luo et al. (Joule, 2025), who
independently confirmed that
solvation-configurational entropy
governs low-temperature interfacial
kinetics, we observed a performance
band at intermediate SCE values
that maximises performance across
both room-temperature and low-
temperature operating conditions.

These results support the hypothesis
that SCE is the unified causal
variable underlying five independent
experimental research programmes
in lithium metal battery electrolyte
design, and suggest that SCE should
be incorporated as a primary
screening criterion in molecular
discovery platforms.

All predictions were pre-registered
before analysis at:
[repository URL]
Timestamp: 2026-04-01
```

---

## DOCUMENT METADATA

```
document_id:
  OC-BATTERY-EMPIRICAL-TEST-PROTOCOL

version:    1.0
date:       2026-04-01
author:     Eric Robert Lawson
            OrganismCore

purpose:
  Complete operational protocol for
  testing the SCE unified causal
  variable hypothesis using publicly
  available data and Python analysis
  on a MacBook Air.

contains:
  — Core prediction being tested
  — Primary and secondary data sources
  — Complete Python scripts (copy-paste
    ready, placeholder data labelled)
  — Data extraction guide
  — Full test sequence (4 tests)
  — Environment setup instructions
  — File structure
  — CSV data format
  — Success criteria
  — Failure diagnostic protocol
  — Statement template for expert
    conversation and journal submission

primary_data_source:
  Journal of Power Sources 2025
  DOI: 10.1016/j.jpowsour.2025.237088
  "Statistical insights into solvation
  sheaths of lithium-ions"
  Already has coordination distribution
  data. SCE computable in one afternoon.

secondary_data_source:
  arXiv 2501.11932 — FREE, no paywall
  Multiple electrolyte types with
  coordination structure data.

band_hypothesis_data:
  Joule 2025 Hunan University paper
  "Solvation-Configurational Entropy
  Governs Interfacial Kinetics in
  Low-Temperature Batteries"
  Their Ssc values + performance data
  can be used directly.

related_documents:
  OC-BATTERY-MASTER-RECORD-2026-04-01.md
  OC-BATTERY-SOLVATION-VARIANCE-
    PREREGISTRATION.md
  OC-BATTERY-VARIANCE-DERIVATION-
    AND-CANDIDATES.md
  OC-BATTERY-LITERATURE-CHECK-
    DERIVED-CANDIDATES.md
  OC-BATTERY-SOLVATION-COHERENCE-
    DERIVATION.md

repository:
  OrganismCore
  attractor-oncology

ORCID: 0009-0002-0414-6544
```

---

*The data already exists.*
*The scripts are here.*
*The prediction is timestamped.*

*One afternoon.*
*One paper.*
*One result.*

*Then you bring the result*
*not the theory.*
