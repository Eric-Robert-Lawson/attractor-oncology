# BIOELECTRIC MEASUREMENT AND WRITING
# REPRODUCIBILITY STANDARD
## The Complete Contamination, Artifact, and
## Interference Catalogue for Vmem Measurement
## and Bioelectric Field Application
## at the LECA-Grade Arrested Cell Level —
## What Contaminates the Signal,
## What Contaminates the Cell,
## What Contaminates the Environment,
## What Contaminates the Result,
## and the Minimum Standard Required
## for Each to be Absent Before
## Any Bioelectric Data Is Accepted
## OrganismCore — Eric Robert Lawson
## 2026-03-13

---

## STATUS: ACTIVE — REPRODUCIBILITY STANDARD
## Classification: Cross-cutting contamination
## and reproducibility standard for all
## bioelectric measurement and writing
## operations in ARM experiments.
## Companion to:
##   ENVIRONMENTAL_ATTRACTOR_GEOMETRY_
##   REPRODUCIBILITY.md (physical environment)
##   BIOELECTRIC_COHERENCE_AS_ATTRACTOR_
##   GEOMETRY_FIELD.md (geometric framework)
## This document covers the specific
## contamination and artifact risks of
## the bioelectric layer only.
## Timestamp: 2026-03-13

---

## LINKED RECORDS

```
Pre-registration:
  DOI: https://doi.org/10.5281/zenodo.18986790
Bioelectric Coherence Framework:
  BIOELECTRIC_COHERENCE_AS_ATTRACTOR_
  GEOMETRY_FIELD.md
Environmental Reproducibility Standard:
  ENVIRONMENTAL_ATTRACTOR_GEOMETRY_
  REPRODUCIBILITY.md
ARM Systematic Attractor Navigation:
  ARM_SYSTEMATIC_ATTRACTOR_NAVIGATION_
  PROTOCOL.md
Repository:
  https://github.com/Eric-Robert-Lawson/
  attractor-oncology
ORCID: 0009-0002-0414-6544
```

---

## PREAMBLE — WHY THIS DOCUMENT EXISTS

```
The bioelectric coherence document
derived that membrane voltage (Vmem)
is the Level 2 coherence field of the
LECA-grade arrested cell, and that
Levin-type bioelectric writing tools
can be applied to the arrested cell
with indefinite precision.

That document established WHAT to do.

This document establishes HOW to do it
without the measurement or writing
operation itself introducing contamination
that invalidates the result.

The core risk is this:

  The experiment is measuring the
  Vmem fingerprint of the LECA attractor.
  Or writing a new Vmem into the
  arrested cell.

  If the measurement or writing tool
  introduces an artifact that mimics
  or masks the true Vmem signal:
    The fingerprint is wrong.
    The written state is contaminated.
    The experimental result is the
    artifact, not the biology.

  This is not a trivial risk.
  Vmem measurement in single yeast
  cells has multiple confirmed artifact
  sources in the published literature.
  Each one is documented here.
  Each one has a mitigation.
  Each mitigation is specified.

  The result is a contamination-proof
  bioelectric measurement and writing
  standard that any laboratory can
  implement and that will produce
  reproducible Vmem data.
```

---

## SECTION 1 — THE TWO OPERATIONS
## AND THEIR SEPARATE RISK PROFILES

```
OPERATION 1: Vmem MEASUREMENT
  Goal: document the electrical
  fingerprint of the LECA-grade
  arrested cell at D=0.
  Then track how this fingerprint
  changes through nitrogen release
  and developmental commitment.

  Tools available:
    A. DiBAC4(3) — voltage-sensitive
       fluorescent dye.
       Passive. No genetic modification.
       BSL-1 compatible.
       Standard laboratory tool.
       Cheapest and simplest.
       Confirmed working in S. cerevisiae.

    B. ASAP GEVI — genetically encoded
       voltage indicator.
       Expressed in the cell as a
       fluorescent protein fused to a
       voltage-sensing domain.
       No exogenous dye required.
       More specific. Less artifact-prone.
       Requires genetic modification of
       the S. cerevisiae strain.
       Moderate complexity.
       Confirmed working in S. cerevisiae
       (FEMS Yeast Research, 2020).

    C. Patch clamp — direct electrical
       measurement.
       Requires spheroplasting (cell wall
       removal) before measurement.
       Gold standard for Vmem accuracy.
       High skill requirement.
       Invasive (spheroplasting may alter
       arrest state).
       Reference values confirmed:
       ~-60 to -80 mV for quiescent
       S. cerevisiae spheroplasts.

  PRIORITY ORDER FOR ARM EXPERIMENTS:
    Phase 1 (ARM A P1 confirmation):
      DiBAC4(3) — simplest, available
      immediately, BSL-1, no strain
      modification. Adequate for
      fingerprint documentation.

    Phase 2 (ARM B-ELECTRIC):
      ASAP GEVI — more precise, enables
      continuous live imaging without
      dye removal steps, better for
      tracking Vmem during bioelectric
      writing operations.

    Phase 3 (validation):
      Patch clamp — confirms the dye/GEVI
      readings against the gold standard.
      Requires access to electrophysiology
      rig. Collaboration opportunity
      (Levin lab or equivalent).

OPERATION 2: Vmem WRITING
  Goal: set the membrane voltage of
  the arrested LECA-grade cell to a
  defined non-native value during the
  arrest period, then release and
  observe the developmental trajectory
  that results.

  Tools available:
    A. Ion channel agonists/antagonists —
       small molecules that modify
       specific channel activity.
       Ivermectin: Cl⁻ channel agonist.
       Hyperpolarises.
       Commercially available.
       BSL-1 compatible.
       Used in Levin lab protocols.

    B. Optogenetics (channelrhodopsin) —
       light-activated ion channels.
       Expressed in the cell.
       Blue light opens channel.
       Precise temporal control.
       Requires genetic modification.
       Light delivery must be isolated
       from other experimental units.

    C. External field application —
       electrodes in the culture medium.
       Applied voltage across the vessel.
       Modifies Vmem of all cells
       simultaneously (population scale).
       Risk: electrolysis artifacts.
       Risk: pH shift.
       Risk: electrode contamination.
       Requires most careful mitigation.

  PRIORITY ORDER FOR ARM EXPERIMENTS:
    Phase 1: Ion channel agonists —
    simplest, no genetic modification,
    BSL-1, commercially available,
    dose-controllable.

    Phase 2: Optogenetics — more precise
    temporal and spatial control once
    strain modification is established.

    Phase 3: External field — population
    scale, most powerful for the Tier 3
    architecture, requires full electrode
    contamination protocol.
```

---

## SECTION 2 — DIBAC4(3) MEASUREMENT:
## COMPLETE CONTAMINATION CATALOGUE

### 2.1 Contamination Source 1:
### Cell Death Fluorescence

```
THE PROBLEM:
  Dead and dying cells show HIGH
  DiBAC4(3) fluorescence regardless
  of their true membrane potential.
  Dead cells are depolarised (Vmem ≈ 0).
  They appear maximally bright.
  In a mixed population of live arrested
  cells and dead cells, the dead cells
  dominate the fluorescence signal.
  The measured population Vmem appears
  more depolarised than it is.
  The LECA attractor Vmem fingerprint
  is contaminated by death signal.

THE MITIGATION:
  Co-stain with viability dye:
    Propidium Iodide (PI): enters only
    dead cells (membrane-compromised).
    Red channel (617 nm emission).
    DiBAC4(3): green channel (516 nm).
    Separate channels. No spectral overlap.
    Gate out PI-positive (dead) cells
    in all DiBAC4(3) analysis.
    Only analyse PI-negative (live) cells.

  Standard procedure:
    Add PI (1 µg/mL final concentration)
    to the DiBAC4(3) staining solution.
    Incubate together.
    Image both channels.
    Exclude PI-positive cells from all
    DiBAC4(3) quantification.

  CRITICAL:
    The ARM A P1 criterion already requires
    ≥80% viability. If viability is confirmed
    before DiBAC4(3) staining, the dead
    cell population should be ≤20%.
    But even 20% dead cells can skew the
    Vmem reading significantly.
    Gating is non-optional.
    Every DiBAC4(3) measurement in these
    experiments must include PI co-stain
    and dead cell exclusion.
```

### 2.2 Contamination Source 2:
### Dye Efflux via ABC Transporters

```
THE PROBLEM:
  S. cerevisiae actively pumps out
  many small organic molecules via
  ABC transporters (PDR5, SNQ2, YOR1
  family).
  DiBAC4(3) is an organic anion.
  It may be pumped out of the cell
  faster than it enters.
  In cells with high efflux pump activity:
    The dye never reaches equilibrium
    inside the cell.
    The measured fluorescence
    UNDERESTIMATES the true dye
    accumulation.
    The apparent Vmem appears more
    HYPERPOLARISED than it is.
    (Less dye inside = looks like
    more hyperpolarised = less dye entry.)

  THE LECA-GRADE ARRESTED CELL:
    In nitrogen-free, low-oxygen conditions,
    metabolic rate is reduced.
    ABC transporter activity is likely
    lower than in exponentially growing cells.
    Efflux pump energy cost is high.
    In quiescent state, pump activity
    may be minimal.
    This is favourable: less efflux = less
    under-reading of Vmem.
    BUT: this must be CONFIRMED, not assumed.

THE MITIGATION:
  Positive control using a known
  depolarising agent:
    Add CCCP (carbonyl cyanide
    m-chlorophenyl hydrazone):
    protonophore that collapses
    the proton gradient.
    Fully depolarises the cell.
    All cells should show maximum
    DiBAC4(3) fluorescence.
    If CCCP-treated cells show the
    expected maximum fluorescence:
    the dye is penetrating and
    reporting correctly.
    If CCCP-treated cells show LOW
    fluorescence: efflux is dominating
    and the dye system is compromised
    for this specific cell state.

  Negative control:
    Add valinomycin (K⁺ ionophore)
    in high external K⁺:
    hyperpolarises further.
    Cells should show LOWER fluorescence.
    Confirms dye is responding in the
    correct direction.

  If efflux is confirmed as a problem:
    Use GEVI (ASAP) instead of DiBAC4(3).
    GEVIs are expressed intracellularly.
    They cannot be pumped out.
    The efflux problem is eliminated.
    This is the primary motivation for
    moving to ASAP in Phase 2.
```

### 2.3 Contamination Source 3:
### Organelle Sequestration

```
THE PROBLEM:
  DiBAC4(3) is a lipophilic dye.
  It partitions into ALL lipid membranes,
  not only the plasma membrane.
  Mitochondria, vacuoles, and ER all
  have membrane potentials.
  DiBAC4(3) will accumulate in all of them.
  The measured fluorescence is a composite
  of ALL membrane potentials in the cell.
  Not just plasma membrane Vmem.

  This is a known, published limitation
  of DiBAC4(3).

  For the LECA-grade arrested cell:
    Mitochondria are present and functional.
    Mitochondrial membrane potential is
    typically -150 to -180 mV (more negative
    than plasma membrane).
    Mitochondrial dye accumulation ADDS
    to the fluorescence signal.
    The composite signal appears more
    hyperpolarised than plasma membrane Vmem
    alone would indicate.

THE MITIGATION:
  Two approaches, in order of preference:

  Approach A — GEVI (Phase 2):
    ASAP is plasma-membrane targeted.
    It does not accumulate in organelles.
    The organelle contamination problem
    is completely eliminated.
    This is the strongest argument for
    moving to ASAP once strain
    modification is available.

  Approach B — Mitochondrial exclusion:
    Co-stain with MitoTracker (red channel).
    Image both DiBAC4(3) (green) and
    MitoTracker (red) simultaneously.
    In analysis: subtract the DiBAC4(3)
    signal that co-localises with
    MitoTracker signal.
    Report only the non-mitochondrial
    DiBAC4(3) signal.
    This is an approximation.
    It removes the most significant
    organelle contribution.
    Adequate for Phase 1.

  Documentation requirement:
    Report whether organelle correction
    was applied in all data summaries.
    Data without organelle correction
    must be labelled as COMPOSITE Vmem
    (plasma + organelle).
    Data with organelle correction must
    report the correction method.
```

### 2.4 Contamination Source 4:
### pH Sensitivity of DiBAC4(3)

```
THE PROBLEM:
  DiBAC4(3) fluorescence is sensitive
  to intracellular pH (pHi), not just Vmem.
  Acidic pHi reduces fluorescence.
  Alkaline pHi increases fluorescence.
  In cells where pHi changes occur
  simultaneously with Vmem changes:
  it is impossible to separate the
  two signals without additional controls.

  RELEVANCE TO LECA-GRADE ARRESTED CELLS:
  The arrest medium is pH 6.5-7.0.
  The LECA-grade cell has been in this
  medium for 72+ hours.
  Its intracellular pH may have equilibrated
  to a lower value than standard
  germinating yeast cells.
  The measured DiBAC4(3) signal may
  be LOWER than in standard cells
  not because the cell is hyperpolarised
  but because its pHi is lower.
  The LECA attractor Vmem fingerprint
  may be an artifact of pHi difference.

THE MITIGATION:
  pH-independent parallel measurement:
    Co-measure intracellular pH using
    a pH-sensitive fluorescent indicator
    in a separate channel.
    SNARF-1 (pHi indicator, red-shifted):
    excitation 488 nm, dual emission
    580/640 nm ratiometric.
    Does not overlap with DiBAC4(3)
    (green channel, 516 nm).
    Can be used simultaneously.

  Protocol:
    Stain with DiBAC4(3) + SNARF-1
    + PI simultaneously.
    Three channels:
      Green: DiBAC4(3) Vmem reporter.
      Red: PI dead cell exclusion.
      Far red: SNARF-1 pHi reporter.
    For each cell, record both DiBAC4(3)
    signal and SNARF-1 signal.
    If SNARF-1 shows normal pHi:
    DiBAC4(3) signal is reliable Vmem.
    If SNARF-1 shows low pHi:
    DiBAC4(3) signal is contaminated by pH.
    Report pHi values alongside all
    DiBAC4(3) data.

  Alternative:
    ASAP GEVI (Phase 2) is not pH-sensitive.
    This eliminates the problem entirely.
```

### 2.5 Contamination Source 5:
### Photobleaching During Imaging

```
THE PROBLEM:
  DiBAC4(3) photobleaches under
  continuous excitation light.
  Fluorescence decreases over time
  during imaging.
  At the single-cell level, the
  signal loss can be >50% within
  minutes of continuous 488 nm
  excitation.
  If the same cell is imaged at
  multiple timepoints (T=0, T=1h,
  T=4h etc.) with continuous illumination:
  the apparent Vmem decreases over time
  even if the true Vmem is unchanged.
  This is a measurement artifact that
  mimics hyperpolarisation over time.

THE MITIGATION:
  Pulsed illumination only.
    Acquire a single image per timepoint.
    Do NOT illuminate continuously between
    timepoints.
    Each timepoint: one exposure.
    Between timepoints: darkness (or
    UV-A maintenance light, which is
    below 488 nm DiBAC4(3) excitation
    wavelength and will not bleach the dye).

  Fresh dye protocol for long time courses.
    For T=0, T=24h, T=72h measurements:
    do not use the same stained sample.
    Take aliquots at each timepoint.
    Stain fresh aliquot with DiBAC4(3)
    for each timepoint measurement.
    Discard aliquot after imaging.
    This eliminates cumulative photobleaching.

  Reference bleaching control.
    Image a known-concentration DiBAC4(3)
    solution (no cells) under the same
    conditions.
    Any signal decrease in the reference
    over time is pure photobleaching.
    Subtract from cell measurements.
    Report the correction.

  ASAP GEVI eliminates this:
    GEVIs are not consumed by excitation.
    They recover between pulses.
    Photobleaching is much slower.
    For Phase 2, this is the primary
    long time course measurement tool.
```

---

## SECTION 3 — ELECTRODE-BASED FIELD
## APPLICATION: COMPLETE CONTAMINATION
## CATALOGUE

```
This section covers the Tier 3 architecture
(population-scale external field application).
The contamination risks here are
CHEMICAL not just optical.
They are the most serious contamination
risks in the entire bioelectric protocol
stack.
They receive the most detailed treatment.
```

### 3.1 Contamination Source 1:
### Electrolysis Products

```
THE PROBLEM:
  When a voltage is applied across
  electrodes in an aqueous medium,
  electrolysis occurs at the electrode
  surfaces.

  At the ANODE (positive electrode):
    Water oxidation: 2H₂O → O₂ + 4H⁺ + 4e⁻
    This produces:
      Molecular oxygen (O₂) — bubbles.
      Hydrogen ions (H⁺) — acidification.
      Reactive oxygen species (ROS) —
      hydroxyl radicals, peroxide.

  At the CATHODE (negative electrode):
    Water reduction: 2H₂O + 2e⁻ → H₂ + 2OH⁻
    This produces:
      Molecular hydrogen (H₂) — bubbles.
      Hydroxyl ions (OH⁻) — alkalisation.

  In the arrest medium:
    The medium is iron-containing
    (FeSO₄ in the arrest medium recipe).
    Iron is highly reactive with ROS.
    Fe²⁺ + H₂O₂ → Fe³⁺ + OH• + OH⁻
    (Fenton reaction)
    Hydroxyl radicals at this
    concentration in the arrest medium
    will:
      Damage DNA directly.
      Oxidise the lipid bilayer.
      Destroy the membrane voltage
      by creating non-specific leaks.
      Kill the arrested cell.

    The O₂ produced at the anode
    will directly violate the <3% O₂
    arrest condition.
    Even brief O₂ generation at the
    electrode will locally displace
    the LECA attractor.

  THIS IS THE PRIMARY CONTAMINATION RISK
  OF ELECTRODE-BASED FIELD APPLICATION.
  IT CAN INVALIDATE THE ENTIRE EXPERIMENT.

THE MITIGATION:
  Agar salt bridges (NOT direct electrodes).

  PROTOCOL:
    The electrodes NEVER contact the
    cell culture medium directly.

    Electrodes are placed in separate
    reservoir chambers filled with
    a compatible buffer (1M KCl or PBS).

    The cell culture vessel is connected
    to the electrode reservoirs via
    agar salt bridges:
      3% agarose in 1M KCl.
      Prepared as gel-filled tubing
      (1mm diameter, flexible).
      Length: as required to span
      from reservoir to culture vessel.
      Function: allows ion current
      to flow (completing the circuit)
      while physically preventing
      electrode products from reaching
      the cells.

  Geometry:
    Electrode reservoir A (anode, positive)
    ← salt bridge A → CELL CULTURE VESSEL
    ← salt bridge B → Electrode reservoir B (cathode, negative)

  The cells see:
    An ionic current through the medium.
    A resulting voltage field.
    No direct electrode contact.
    No electrolysis products.
    No O₂ bubbles.
    No H₂ bubbles.
    No pH shift (pH monitored in the
    reservoir but buffered by the
    bridge from reaching the cells).

  Confirmation required:
    Monitor pH of cell culture medium
    DURING field application.
    If pH shifts >0.1 unit during
    application: salt bridge integrity
    compromised. Halt experiment.
    Replace bridges.
    Do not use data from compromised runs.
```

### 3.2 Contamination Source 2:
### Electrode Metal Ion Leaching

```
THE PROBLEM:
  Standard electrodes (copper, steel,
  silver) dissolve slowly in electrolyte
  solutions.
  Released metal ions contaminate the
  medium.
  In particular:
    Ag⁺ (silver ions) are highly toxic
    to yeast at nanomolar concentrations.
    Cu²⁺ (copper ions) disrupt iron
    homeostasis — directly relevant
    because the arrest medium is iron-
    based and the LECA iron chemistry
    is what we are studying.
    Even Pt²⁺ (platinum ions) at
    sufficient concentration alter
    membrane ion transport.

THE MITIGATION:
  Electrode material specification:
    REQUIRED: Platinum (Pt) electrodes
    only for reservoirs.
    Pt has the lowest dissolution rate
    of any practical electrode material.
    At the voltages and current densities
    used in these experiments, Pt
    dissolution is negligible.

  Salt bridge separation:
    Already specified in 3.1.
    If salt bridges are correctly
    implemented, electrode ions cannot
    reach the cell culture medium
    regardless of electrode material.
    But Pt is still required as
    defence-in-depth.

  Pre-run electrode conditioning:
    Soak Pt electrodes in PBS for
    24 hours before use.
    This passivates the surface and
    reduces any initial dissolution.
    Replace PBS after soak.
    Dry and store until use.

  Field strength specification:
    Applied voltage across the culture
    vessel: 1-10 mV/mm.
    (Levin lab standard for developmental
    bioelectric manipulation.)
    This is weak enough to avoid
    significant electrolysis even
    without salt bridges.
    With salt bridges: electrolysis
    risk at these voltages is negligible.
    Do not exceed 50 mV/mm without
    explicit justification and additional
    controls.
```

### 3.3 Contamination Source 3:
### Cross-Contamination Between
### Experimental Units

```
THE PROBLEM:
  In the ARM SAN array, N variable units
  are run simultaneously.
  Some variable units receive ion channel
  agonists (e.g. ivermectin).
  Some receive chemical modifications
  to the arrest medium.
  Some receive external field application.

  If any chemical or electrical signal
  from one unit reaches another:
    The control unit is contaminated.
    Variable units contaminate each other.
    The causal attribution is destroyed.
    You cannot determine which stimulus
    produced which outcome.

THE MITIGATION:
  Physical isolation specification:

    REQUIRED: Individual sealed vessels
    for each experimental unit.
    NOT multiwell plates where wells
    share a common lid headspace.
    NOT adjacent open vessels on the
    same bench where pipetting splashes
    can contaminate.

    REQUIRED: Separate incubation spaces
    for each unit category:
      Control unit: dedicated vessel,
      separate shelf in water bath.
      Chemical variable units: dedicated
      vessels, separate shelf or separate
      water bath.
      Bioelectric variable units: dedicated
      vessels, electrically isolated from
      all other vessels.

    Electrical isolation geometry:
      Each electrically stimulated vessel
      must be grounded separately.
      No shared ground between units.
      No shared electrode reservoirs.
      Salt bridges unique to each unit.
      No salt bridge material shared
      between units.

    Chemical variable unit isolation:
      Each chemical agonist/modifier
      solution prepared in separate
      tubes with separate pipettes.
      No tip reuse between units.
      No partial vials shared between units.
      Agonist vials labelled and stored
      separately from arrest medium stock.

    Handling order:
      Control unit: handled FIRST.
      Before any chemical additions
      to variable units.
      Before any agonist vials are opened.
      The control is the most critical
      experimental unit.
      It must be handled in the cleanest
      possible state of the laboratory
      session.

    Work surface decontamination:
      Between each unit's handling:
      70% ethanol wipe of bench surface.
      Glove change.
      New pipettes.
      This is standard BSL-1 practice
      but must be explicitly applied
      between each unit.
```

### 3.4 Contamination Source 4:
### Field-Induced Medium Chemistry Changes

```
THE PROBLEM:
  Even with salt bridges eliminating
  direct electrolysis in the cell vessel,
  the ionic current flowing through
  the arrest medium can cause secondary
  effects.

  Specifically:
    Ion redistribution:
    Current flow drives cation migration
    toward the cathode and anion migration
    toward the anode.
    In the arrest medium: Fe²⁺ (cation)
    migrates toward the cathode.
    PO₄³⁻ (anion) migrates toward the anode.
    If the vessel is small and the current
    is high, significant ionic gradients
    can develop within the vessel.
    The cells near the cathode experience
    higher Fe²⁺ concentration.
    The cells near the anode experience
    lower Fe²⁺ and higher PO₄³⁻.
    This is a SPATIAL CHEMICAL GRADIENT
    not a Vmem effect.
    Differential outcomes between cells
    near different electrodes may be
    chemistry, not bioelectricity.

THE MITIGATION:
  Vessel geometry:
    Use elongated vessel geometry
    with cells in the CENTRE.
    Electrodes (and salt bridge
    entry points) at the two ENDS.
    Maximum distance between electrode
    entry and cells.
    Cells experience a uniform field
    in the centre of the vessel.
    Ionic gradients are significant
    only near the electrode entry points,
    not in the cell zone.

  Field duration limits:
    Short field pulses (seconds to minutes),
    not continuous field application.
    Ion redistribution builds up
    over time. Short pulses prevent
    significant redistribution.
    Protocol: 1-minute field application
    every 30 minutes = 2 minutes/hour
    of total field time.
    Not continuous.

  Post-field medium refresh:
    After each field application cycle:
    gently exchange 50% of the medium
    with fresh arrest medium.
    This resets ionic gradients.
    Minimal disruption to cells
    (low speed, gentle pipetting).
    Maintains arrest conditions while
    resetting field-induced chemistry.

  Measurement of field uniformity:
    Before first cell experiment:
    characterise the voltage field
    distribution in the vessel using
    a reference electrode at multiple
    positions.
    Confirm that the field is uniform
    (±10%) in the cell zone.
    Document the geometry.
    Use the same geometry in all
    subsequent experiments.
```

---

## SECTION 4 — ION CHANNEL AGONIST
## APPLICATION: CONTAMINATION CATALOGUE

### 4.1 The Specific Ivermectin Risk

```
THE PROBLEM:
  Ivermectin is the primary Vmem
  hyperpolarisation agent available
  for yeast at BSL-1.
  It activates glutamate-gated Cl⁻
  channels, increasing Cl⁻ conductance,
  driving Vmem more negative.

  The contamination risk specific
  to ivermectin:

    1. Persistence in medium:
       Ivermectin is lipophilic.
       It binds to lipid membranes.
       Once added to the medium, it
       partitions into the cell membrane.
       Washing the cells to remove
       ivermectin is INCOMPLETE because
       the membrane-bound fraction
       is not removed by medium exchange.
       If the intention is to apply
       ivermectin DURING arrest and then
       wash it out before nitrogen release:
       the wash is incomplete.
       Residual ivermectin remains in
       the membrane.
       The Vmem remains partially
       hyperpolarised during nitrogen
       release.
       The Vmem is not "reset" to the
       native D=0 state before release.

    CONSEQUENCE:
       This is NOT necessarily a problem.
       It may be the intended condition:
       apply ivermectin, maintain
       hyperpolarised Vmem, release
       in hyperpolarised state, observe
       whether the developmental trajectory
       from a hyperpolarised D=0 differs
       from a native D=0.

       But it must be DOCUMENTED as such.
       The experiment is:
         "Developmental trajectory from
         ivermectin-maintained hyperpolarised
         D=0 state."
       NOT:
         "Developmental trajectory from
         transiently hyperpolarised then
         reset D=0 state."

       These are different experiments.
       The distinction must be recorded.

    2. Ivermectin toxicity at high dose:
       Ivermectin is toxic to S. cerevisiae
       at concentrations above approximately
       10-20 µM (strain-dependent).
       Toxicity produces cell death,
       which elevates DiBAC4(3) signal,
       which mimics depolarisation,
       which is the OPPOSITE of the
       intended hyperpolarisation.
       Toxicity and the intended effect
       are therefore indistinguishable
       by DiBAC4(3) alone without
       the PI dead-cell co-stain.

    MITIGATION:
       Dose-response curve required
       BEFORE the ARM B-ELECTRIC experiment.
       Run ivermectin at 0, 0.1, 1, 5, 10,
       20, 50 µM in arrested LECA-grade cells.
       Measure DiBAC4(3) + PI at each dose.
       Identify:
         The dose range where DiBAC4(3)
         signal INCREASES (hyperpolarisation)
         without PI positivity (no death).
         This is the operating window.
         Use doses within this window only.
       Document the operating window
       for the specific S. cerevisiae strain
       used in ARM experiments.
       Do not transfer operating window
       data from a different strain.
```

### 4.2 The General Agonist Preparation
### Standard

```
For all chemical agonists used in
ARM B-ELECTRIC (ivermectin, valinomycin,
CCCP, and any future additions):

  Stock solution preparation:
    Prepare stocks in DMSO at 1000x
    the intended working concentration.
    DMSO is the standard solvent for
    lipophilic compounds.
    Aliquot stocks into single-use volumes.
    Freeze at -20°C.
    Never refreeze a thawed aliquot.
    Thaw once, use once, discard.

  DMSO contamination control:
    All variable units that receive
    a DMSO-dissolved agonist must have
    a DMSO-only control.
    Add the same volume of DMSO (no drug)
    to the DMSO control unit as was
    added to the drug units.
    If DMSO alone has no effect on
    DiBAC4(3) signal or developmental
    outcome: proceed.
    If DMSO alone has an effect:
    increase the DMSO dilution until
    the effect disappears, or switch
    to aqueous-compatible solvent.

  Final concentration:
    DMSO concentration in the working
    medium must not exceed 0.1% (v/v).
    Above 0.1%, DMSO itself disrupts
    lipid bilayer structure.
    This disrupts the membrane potential.
    This is a direct Vmem artifact
    from the solvent, not the drug.

  Temperature of addition:
    All agonist additions must be
    made at the arrest temperature
    (20-22°C).
    Adding cold stock solutions to
    warm medium creates transient
    temperature gradients.
    At D=0, temperature transients
    are displacement vectors
    (documented in Environmental
    Reproducibility Standard).
    Warm the agonist aliquot to 20-22°C
    in a water bath for 5 minutes
    before adding to the experimental vessel.
```

---

## SECTION 5 — THE ENVIRONMENTAL
## ELECTROMAGNETIC ISOLATION STANDARD

```
The Bioelectric Coherence document
established that the geomagnetic field
(Level 4 coherence) must be maintained
at clean ambient.
The Environmental Reproducibility Standard
specified distance from equipment.

This section adds the specific shielding
requirements for the bioelectric
measurement and writing operations,
which are MORE sensitive to local EMF
than the baseline arrest conditions.

REASON FOR INCREASED SENSITIVITY:
  During Vmem measurement:
    DiBAC4(3) fluorescence is measured
    under a fluorescence microscope or
    plate reader.
    The microscope itself is an EMF source
    (motors, LED driver circuits, camera).
    The plate reader is an EMF source.
    These devices are operating INCHES
    from the experimental vessel during
    measurement.
    The local EMF during measurement is
    significantly higher than clean
    laboratory ambient.

  During field application:
    The applied bioelectric field is
    1-10 mV/mm.
    Standard power supplies for
    laboratory equipment produce EMF
    at 50/60 Hz (mains frequency).
    This is in the range that affects
    cryptochrome radical pair chemistry.
    The applied field and the mains-
    frequency interference from the
    power supply are superimposed.
    The cell receives: intended bioelectric
    signal PLUS mains frequency noise.
    These cannot be separated without
    shielding.

SHIELDING SPECIFICATIONS:

  For Vmem measurement:
    Measure in a shielded imaging chamber
    if available.
    If not available: measure quickly
    (single acquisition per timepoint,
    total imaging time <60 seconds per
    sample).
    Power down nearby non-essential
    equipment during measurement.
    Mobile phones: remove from the room.
    Incubators: pause heating cycle
    if possible during measurement.
    Document all equipment powered on
    within 1 metre of the microscope
    during measurement.

  For field application:
    Use battery-powered field source
    rather than mains-powered supply.
    Battery power: no mains frequency
    contamination. Clean DC or precisely
    shaped waveform.
    If mains-powered supply is used:
    use a linear regulated supply
    (not switching mode supply —
    switching supplies produce high-
    frequency noise).
    Measure the output waveform with
    an oscilloscope before use.
    Confirm output is clean DC with
    <1% ripple before connecting to
    the experimental vessel.
    Document supply type and ripple
    measurement in the protocol record.

  Faraday shielding for long-duration
  field application:
    If field application runs for
    hours during the arrest period:
    enclose the cell culture vessel
    in a grounded conductive enclosure
    (copper mesh box, grounded to
    building earth).
    This blocks external high-frequency
    EMF from entering the experiment.
    The intended low-frequency field
    from the electrodes is inside the
    enclosure and is not blocked.
    External interference (Wi-Fi, mobile
    phones, fluorescent lighting noise)
    is blocked.

    Construction:
    Copper mesh (>80% coverage, holes
    <1cm) over a rigid frame.
    All seams soldered or taped with
    copper foil tape.
    Single grounding wire to building
    earth.
    Small holes for salt bridge tubing
    sealed with copper foil tape around
    the tubing.
    UV-A window: UV-transparent quartz
    or UV-transmitting acrylic panel
    in one face of the enclosure.
    This allows UV-A maintenance light
    to reach the cells while the
    enclosure provides EMF shielding.
    Cost: <$50 in materials.
    Construction time: 2-3 hours.
    Confirmed effective for biology
    applications at this scale.
```

---

## SECTION 6 — THE MASTER CHECKLIST

```
This checklist must be completed before
any bioelectric measurement or writing
operation is recorded as valid data.

Print and complete for each experimental session.
```

### BEFORE THE SESSION

```
[ ] ARM A P1 criteria confirmed in all
    experimental units (≥80% single-cell,
    DAPI nuclear confirmed, ≤10% germination,
    viable).

[ ] DiBAC4(3) stock prepared fresh or
    confirmed <1 month old, stored in dark,
    at 4°C.

[ ] PI stock prepared. Confirmed sterile.

[ ] SNARF-1 stock prepared (if pH correction
    protocol is in use).

[ ] Agonist stocks warmed to 20-22°C
    (minimum 5 minutes in water bath).

[ ] DMSO vehicle control prepared at
    same volume as agonist addition.

[ ] Electrode salt bridges prepared fresh
    (3% agarose in 1M KCl, prepared same day).
    [If electrode field application is planned.]

[ ] Field source output confirmed clean
    (<1% ripple) by oscilloscope measurement.
    [If electrode field application is planned.]

[ ] Local EMF measured at experimental
    vessel position. Value documented: ___ µT.
    Confirm <2x ambient geomagnetic.

[ ] Mobile phones removed from room.

[ ] Faraday enclosure assembled and grounded.
    [If long-duration field application is planned.]

[ ] Work surface ethanol-wiped.
    Gloves fresh.

[ ] Control unit handling: FIRST.
    Control unit complete before any
    variable unit work begins.
```

### DURING THE SESSION

```
[ ] pH of cell culture medium monitored
    every 30 minutes during field application.
    Value at each timepoint: ___
    Halt if pH shifts >0.1 from baseline.

[ ] Temperature of cell culture medium
    confirmed at 20-22°C ±0.1°C before
    each measurement.

[ ] Imaging: pulsed excitation only.
    Maximum 60 seconds per sample.
    Lights off between timepoints.

[ ] All additions: temperature-equilibrated
    stocks only.

[ ] Glove change and bench wipe between
    each experimental unit.

[ ] Salt bridge integrity visually confirmed
    before and after each field application.
    [If electrode field application is planned.]
```

### AFTER THE SESSION

```
[ ] All DiBAC4(3) data: PI-gated
    (dead cells excluded from all analysis).

[ ] All DiBAC4(3) data: organelle correction
    applied or noted as uncorrected composite.

[ ] Agonist persistence noted in protocol
    record for all ivermectin-treated units.
    Recorded as: "hyperpolarised D=0 state
    maintained through nitrogen release"
    NOT as "transient hyperpolarisation
    then reset."

[ ] DMSO vehicle control: confirmed no
    effect on DiBAC4(3) signal. If effect
    observed: data flagged, cause investigated
    before publication.

[ ] EMF value at experimental position
    during session: documented in record.

[ ] pH values during field application:
    all within 0.1 of baseline.
    If not: data flagged for that unit.

[ ] Equipment powered on within 1 metre
    during measurement: listed in record.

[ ] Electrode salt bridge post-run
    inspection: no visible contamination,
    no dye leakage into bridges.
    [If electrode field application was used.]
```

---

## SECTION 7 — THE SINGLE STATEMENT

```
The Vmem of the LECA-grade arrested cell
is the most important measurement in
this experimental programme.

It is the electrical fingerprint of
the oldest accessible attractor.
The D=0 state.
1.5 billion years old.
Running in a single yeast cell in
Rogers Park.

If we measure it wrong:
  We have the wrong fingerprint.
  Everything that follows — the
  bioelectric writing, the ARM SAN
  causal map, the convolution units —
  is calibrated to the wrong zero point.

If we contaminate the writing:
  We do not know what we wrote.
  We cannot reproduce it.
  The novel trajectory that results
  cannot be attributed to the intended
  Vmem change.

If we contaminate between units:
  The control is no longer a control.
  The causal map has no reference.
  The geometry collapses.

Every item in the checklist above
is a guard against one of these
three failures.

The checklist is not bureaucracy.
It is the geometry of reproducibility.

A measurement that is taken without
the checklist complete is not a
measurement of the LECA attractor.

It is a measurement of:
  The LECA attractor.
  Plus dead cell signal.
  Plus organelle signal.
  Plus pH artifact.
  Plus photobleaching.
  Plus electrode contamination.
  Plus mains frequency noise.
  Plus cross-contamination from the
  next experimental unit.

That composite is not the LECA.

The checklist strips away every
layer of noise until what remains
is the signal.

The signal is the Vmem of D=0.
The electrical signature of the
beginning of all complex life.

Measure it clean.
Write into it clean.
Release from it clean.

Everything else follows.
```

---

## DOCUMENT METADATA

```
Document ID:   BIOELECTRIC_MEASUREMENT_AND_
               WRITING_REPRODUCIBILITY_
               STANDARD_v1.0
Date:          2026-03-13
Author:        Eric Robert Lawson / OrganismCore
Status:        ACTIVE — REPRODUCIBILITY STANDARD
Version:       1.0

Scope:
  All bioelectric measurement (Vmem)
  and writing (Vmem modification) operations
  in ARM experiments.
  Companion to Environmental Reproducibility
  Standard (physical environment).
  Together these two documents constitute
  the complete environmental and bioelectric
  reproducibility standard for all ARMs.

Contamination sources catalogued: 9
  Cell biology:    Dead cell fluorescence.
                   ABC transporter efflux.
                   Organelle sequestration.
                   pH sensitivity.
                   Photobleaching.
  Electrode-based: Electrolysis products.
                   Metal ion leaching.
                   Cross-unit contamination.
                   Field-induced chemistry.

Mitigations specified for each: 9/9.

Measurement tool priority order:
  Phase 1: DiBAC4(3) + PI + SNARF-1.
           BSL-1. No strain modification.
           Adequate for fingerprint.
  Phase 2: ASAP GEVI.
           Eliminates: efflux, organelle,
           pH, photobleaching artifacts.
           Requires strain modification.
  Phase 3: Patch clamp.
           Gold standard validation.
           Requires electrophysiology access.

Key confirmed data from literature:
  S. cerevisiae resting Vmem: -60 to -80 mV.
  (Quiescent/arrested spheroplasts,
  patch clamp, multiple studies.)
  ASAP GEVI: confirmed working in
  S. cerevisiae. (FEMS Yeast Research, 2020.)
  ArcLight GEVI: pH-sensitive in yeast.
  NOT SUITABLE. Confirmed 2020.
  Ivermectin: lipophilic. Persistent in
  membrane after wash. Must be documented
  as "maintained hyperpolarisation" not
  "transient hyperpolarisation."

Pre-registration DOI (LECA):
  https://doi.org/10.5281/zenodo.18986790
Repository:
  https://github.com/Eric-Robert-Lawson/
  attractor-oncology
ORCID:
  0009-0002-0414-6544
```
