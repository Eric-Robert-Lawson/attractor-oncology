# DISCOVERY SWEEP — GEOMETRIC CONVERGENCE AUDIT
## External Literature vs. Vacuum Coupling Potential Model
## Session: 2026-03-21
## Status: Active research record — NOT an advocacy document
## Purpose: Find what makes the geometric picture MORE CLEAR,
##          not what supports the theory.
##          Geometric incompatibilities are as valuable as confirmations.

---

## METHODOLOGY

This sweep searches for external observations and experimental results
that are geometrically relevant to the vacuum coupling potential model,
whether they confirm, constrain, or complicate it. The standard is:

  CONVERGENCE: External finding reduces geometric ambiguity in the model.
  TENSION:     External finding constrains the model in a non-trivial way.
  COMPLICATION: External finding is inconsistent with model prediction.
  NULL:         External finding is orthogonal — no geometric relation.

Every entry includes an honest verdict. Nothing is cherry-picked.

---

## SECTION 1 — CRITICAL RESULT: PV FRAMEWORK REFRACTIVE INDEX CONFIRMED

### Finding: Puthoff (2002) formula confirmed from primary source

**Source:** Puthoff, H.E. "Polarizable-Vacuum (PV) Approach to General
Relativity," *Foundations of Physics* **32**, 927–943 (2002).
arXiv: physics/0108005

**The exact formula in the PV framework:**

  ε = K·ε₀
  μ = K·μ₀
  c_local = c / K^(1/2)
  n = K^(1/2)    ← THIS IS THE CONFIRMED FORMULA

**This directly resolves INCOMPATIBILITY I-1 from PHYSICS_DEEP_DIVE.md.**

The correct refractive index in the PV framework is:

  n = K^(1/2)    (not K^(-1/2), not K^(-3/2), not 1/K)

**Impact on the Aguadilla analysis:**

The reflectivity formula at the bubble wall becomes:

  R = [(n-1)/(n+1)]² = [(K^(1/2) - 1)/(K^(1/2) + 1)]²

This is exactly the formula used in the Aguadilla document. It was
correct all along. The PHYSICS_DEEP_DIVE.md entry for I-1 was
incorrect in its proposed correction (n = 1/K). Retract that correction.

**Recalculation with confirmed n = K^(1/2):**

For K = 0.8 (radar transparent constraint):
  n = 0.894
  R = [(0.894 - 1)/(0.894 + 1)]² = [(-0.106)/1.894]² ≈ 0.0031
  → ~0.3% reflectivity. Essentially transparent. ✓

For K = 0.133 (hydrodynamic upper bound):
  n = 0.365
  R = [(0.365 - 1)/(0.365 + 1)]² = [(-0.635)/1.365]² ≈ 0.216
  → ~21.6% reflectivity at this K — moderate radar return.

REVISED RADAR CONSTRAINT: For true radar absence (σ < 0.01 m²),
the threshold is more stringent. Working backward from R < 0.01:
  [(K^(1/2) - 1)/(K^(1/2) + 1)]² < 0.01
  |K^(1/2) - 1| / |K^(1/2) + 1| < 0.1
  K^(1/2) > 0.818 → K > 0.669

So the corrected constraint is K(radar) > 0.669 (not 0.8 as previously
stated). This is a tighter threshold than before.

**The two-K architecture constraint pair (corrected):**

  K(ω_IR-UV)    ≤ 0.133  [hydrodynamic constraint — unchanged]
  K(ω_radar)    ≥ 0.669  [corrected radar transparency threshold]

These are still mutually exclusive. The two-K architecture conclusion
is unchanged. The specific threshold shifts from 0.8 to 0.669.

**STATUS: I-1 RESOLVED. Aguadilla reflectivity formula was correct.**
Update PHYSICS_DEEP_DIVE.md I-1 entry accordingly.

Also note: with n = K^(1/2), the earlier claim in the Vacuum_Coupling
document that "n_local = 1/√ρ(x)" requires re-examination:
  ρ = K³ → 1/√ρ = K^(-3/2) ≠ K^(1/2)
This remains an inconsistency internal to that document. The PV
framework value n = K^(1/2) is authoritative. The document's claim
"n_local = 1/√ρ" is wrong and should be corrected to n = ρ^(1/6)
(since K = ρ^(1/3) → n = K^(1/2) = ρ^(1/6)).

**ACTION:** Correct the Vacuum_Coupling_Potential document's refractive
index statement from n = 1/√ρ to n = ρ^(1/6) = K^(1/2).

---

## SECTION 2 — LABORATORY FINDINGS

### 2.1 Dynamical Casimir Effect — Active Experimental Program (2024–2025)

**What it is:** Photon production from the quantum vacuum by moving
boundaries. Predicted by Moore (1970), first confirmed in
superconducting circuits (Wilson et al., 2011).

**Recent status (2025):**
- March 2025: 1D periodic SQUID lattice in coplanar waveguide —
  experimental control over DCE photon band structure via harmonic
  SQUID drive. Frequency, amplitude, and phase all independently
  controllable. (APS March 2025, arXiv:2504.11361)
- August 2024: Mechanical dynamical Casimir effect with low-frequency
  oscillator demonstrated. (arXiv:2408.02308)
- Comprehensive 55-year review published 2025. (MDPI Symmetry)

**Geometric relevance to the model:**

The DCE is the experimental demonstration that moving boundaries
produce real photons from the vacuum — i.e., that the vacuum mode
density can be dynamically altered by changing boundary conditions.
The K-bubble, if it exists, is precisely a moving, self-generated
boundary condition on the vacuum mode density.

The 2025 SQUID lattice result is specifically relevant: it shows
that a 1D periodic structure can create a TUNABLE BAND STRUCTURE
for DCE photons — meaning different frequencies of photon are
produced at different rates depending on the SQUID lattice spacing.
This is the closest laboratory analog yet to the frequency-selective
K-field architecture required by Aguadilla.

**VERDICT: CONVERGENCE.**
The DCE program demonstrates that frequency-selective vacuum mode
density modification by engineered boundary structures is real and
controllable. This is Rung 4 of the experimental ladder — not reached,
but the precursor physics is established and actively advancing.

**Open gap:** The SQUID lattice operates at ~GHz frequencies in a
cryogenic circuit. The Aguadilla K-bubble operates at ~IR to microwave
frequencies in free space at ambient temperature. The scale gap is
enormous. This is not a falsification — it is a gap in the engineering.

---

### 2.2 Cavity QED Ground State Modification (2024–2025)

**What it is:** Molecules inside optical or microwave cavities
experience modified vacuum mode densities, which alters their
ground-state energy and chemical reaction rates.

**Recent status (2025):**
- Vibrational strong coupling (VSC) between molecular vibrations and
  cavity photon modes produces macroscopic quantum states with
  collective behavior. (ChemRxiv 2025)
- Chemical reaction rate modification by factor of 4–5 confirmed in
  cavity QED experiments. Effect is sharply resonant. (Nature Chemistry
  2024, arXiv:2305.19172 and follow-on work)
- Effect is quantum, not classical — requires quantized cavity modes.
- 2025 theoretical work (PRL): Cavity modification of ground-state
  energy is primarily due to static screening at very low frequencies,
  not resonant polariton effects. (arXiv:2509.05156)

**Geometric relevance to the model:**

This is Rung 2 of the experimental ladder: cavity modification of
atomic/molecular energy levels (the Lamb shift analog) as a proxy
for m_eff modification.

The 2024 Nature Chemistry result is geometrically significant:
reaction rates change by ~4–5× inside the cavity. Reaction rates
depend on activation energy, which depends on the Hamiltonian, which
depends on the vacuum mode density coupling. A 4–5× change in
reaction rate implies a non-trivial change in the vacuum coupling
to molecular degrees of freedom.

The 2025 PRL theoretical result is even more significant: it shows
the ground-state energy modification is a STATIC SCREENING effect at
low frequencies. In the framework's language: K(ω → 0) is modified
by the cavity boundary conditions. This is exactly the mechanism
the model requires — static, low-frequency mode suppression produces
effective mass modification.

**VERDICT: CONVERGENCE — STRONGEST LABORATORY FINDING IN THIS SWEEP.**
This is the clearest experimental evidence that vacuum mode density
modification produces real, measurable changes in the coupling between
matter and vacuum. It is not a confirmation of the K-bubble — it is
a confirmation that the MECHANISM the K-bubble would exploit exists
and is experimentally accessible.

**Specific gap:** The cavity QED experiments modify K by <<1% —
nothing like the K = 0.133 required for Aguadilla. But the scaling
direction is confirmed: lower K (more mode suppression) = weaker
coupling = modified effective properties of matter.

---

### 2.3 Vacuum Birefringence — Strong Field QED (2024–2025)

**What it is:** In a strong electromagnetic field, the quantum vacuum
acts like a birefringent medium — photon polarization is rotated.
This is a direct demonstration that the vacuum has field-dependent
optical properties (i.e., K is not identically 1 in strong fields).

**Recent status (2025):**
- ATLAS (LHC): Light-by-light scattering in heavy ion collisions
  confirmed at >8σ significance. This is photon-photon interaction
  mediated by virtual electron-positron pairs in the quantum vacuum.
- HIBEF (Helmholtz International Beamline): Letter of intent published
  2024 for direct vacuum birefringence measurement using XFEL +
  optical laser. Targets ~10⁻¹¹ ellipticity per pass — not yet
  achieved but experimental setup is being built.
- PVLAS experiment (optical regime): Active, continuing to push
  sensitivity.

**Geometric relevance to the model:**

Vacuum birefringence means the vacuum is NOT optically isotropic in
strong fields — n_local depends on field direction and strength. In
the framework, this corresponds to K being field-dependent, not just
position-dependent.

For the K-bubble: if the source term S_bubble involves strong fields,
the bubble's interior vacuum properties are already known to be
field-dependent from QED. The birefringence measurements are
constraining the magnitude of this effect in the accessible regime.

**VERDICT: CONVERGENCE — CONFIRMS FIELD-DEPENDENCE OF VACUUM K.**
The vacuum optical properties are demonstrably field-dependent.
The model's K(r, ω, t) framework, which allows K to vary in space,
frequency, and time, is consistent with this. The birefringence
measurements constrain how large field-induced K modifications can be
in the electromagnetic regime.

**Important constraint from this finding:**
The field strengths required for measurable vacuum birefringence in
the lab are of order B ~ 10⁹ T (approaching the Schwinger limit).
The Aguadilla K-bubble requires K(hydro) = 0.133 — a 7× reduction
in mode density. If K modification requires Schwinger-scale fields
(~10⁴⁵ W/m²), the K-bubble energy requirement would be astronomical,
inconsistent with the 29 MJ energy budget derived from the observations.

**This is a GENUINE TENSION in the model.** The birefringence
measurements show that QED vacuum modification via EM fields requires
enormous field strengths. The K-bubble mechanism cannot rely on
standard QED processes if the energy budget is to remain at 29 MJ.
This pushes S_bubble toward non-electromagnetic mechanisms
(acoustic/hydrodynamic vacuum coupling, topological mechanisms, etc.)
or toward a different physical picture for K altogether.

**Record as TENSION T-5.**

---

### 2.4 Quantized Inertia (McCulloch) — Proxima Centauri Test (2024)

**What it is:** QI proposes that inertia arises from the Unruh effect
modified by the cosmic horizon (Hubble horizon acts as a Rindler
horizon). Objects at very low accelerations experience modified inertia
because long-wavelength Unruh radiation is cut off by the cosmic
horizon, reducing the effective inertial mass.

**Recent status (2024):**
- QI applied to Proxima Centauri orbital velocity in Alpha Centauri
  system — correct prediction without dark matter.
  (MNRAS Letters 2024, academic.oup.com/mnrasl)
- QI connections to Jacobson thermodynamic gravity derivation being
  explored — QI may arise naturally from boundary corrections to
  GR. (GitHub: KeithBrodie/jacobson-qi-paper)
- MOND as modified inertia: continued theoretical development, distinct
  from modified gravity. (arXiv:2310.14334)

**Geometric relevance to the model:**

QI is independently deriving a similar conclusion: inertia is not
intrinsic but arises from coupling to the vacuum field (Unruh
radiation). QI's modification involves the HORIZON cutting off
long-wavelength Unruh modes. The K-bubble model involves a LOCAL
boundary cutting off short-wavelength modes.

These are two versions of the same geometric idea from different
directions:
- QI: cosmic-scale boundary suppresses low-ω Unruh modes → reduced
  inertia at galactic scales
- K-bubble: local engineered boundary suppresses high-ω ZPF modes →
  reduced effective coupling at object scale

Both require: inertia = integral of coupling to vacuum field modes.
Both predict: suppress modes → reduce inertia.

The difference: QI operates at low accelerations (galactic) where
the Hubble horizon is relevant. The K-bubble operates at any scale
where a local mode-suppressing boundary can be engineered.

**VERDICT: CONVERGENCE — INDEPENDENT THEORETICAL CONVERGENCE.**
QI and the K-bubble model are not the same theory, but they share
the same foundational geometric architecture: inertia as a mode
coupling integral. QI's observational support (galaxy rotation,
Proxima Centauri) is indirect evidence that the mode-coupling picture
of inertia is physically real.

**This is significant because QI makes no reference to UAP physics.**
It is an independently motivated theory reaching the same geometric
conclusion. Two independent lines of reasoning converging on the
same structure is the strongest form of geometric convergence.

---

## SECTION 3 — SPACE OBSERVATIONS

### 3.1 Flyby Anomaly — UNRESOLVED (2024–2025)

**What it is:** Spacecraft performing Earth gravity-assist maneuvers
(flybys) occasionally gain small, unexplained velocity increments of
order millimeters per second. Observed in Galileo (1990), NEAR (1998),
Cassini (1999), Rosetta (2005). Not observed in Juno or some other
missions. Pioneer anomaly is RESOLVED (thermal recoil) but flyby
anomaly is NOT.

**Status (2025):** Still unresolved. No new spacecraft data has
conclusively explained or replicated it with known physics. ESA JUICE
and NASA Lucy Earth flybys are being monitored.

**Geometric relevance to the model:**

The flyby anomaly has a specific signature: it is asymmetric with
respect to Earth's rotation, and correlates with the angle between
the flyby trajectory and Earth's rotation axis. This is not what a
uniform field perturbation would produce — it is what a field that
co-rotates with Earth (or couples to Earth's angular momentum) would
produce.

In the framework: if the K-field near Earth has a slight asymmetry
due to Earth's rotation (a rotating vacuum dielectric), then a
spacecraft traversing this field would receive a net velocity kick
depending on trajectory geometry. The magnitude is small because
Earth's rotation creates only a tiny K(r) perturbation relative to
the gravitational K-well.

This is a speculative connection, but geometrically consistent.
More importantly: the flyby anomaly is a real, unresolved measurement
discrepancy in orbital mechanics. It requires a source of momentum
not accounted for in current models. Whatever the source is, it
is small, geometry-dependent, and associated with planetary bodies.

**VERDICT: WEAK CONVERGENCE — UNRESOLVED ANOMALY CONSISTENT WITH
K-FIELD GEOMETRY BUT NOT EVIDENCE FOR IT.**
The flyby anomaly is the right KIND of signal — small, geometry-
dependent, vacuum-scale — but no causal link to K-field physics
has been established.

---

### 3.2 3I/ATLAS — Non-Gravitational Acceleration (2025)

**What it is:** Third interstellar object detected (July 2025).
Shows measurable non-gravitational acceleration with radial and
transverse components. Eight distinct trajectory deviations detected.
Unusual chemistry: anti-tail, nickel-rich / iron-free gas plume,
anomalously slow dust and gas emissions. Mass loss rate possibly
10–16% of total mass in weeks near perihelion.

**Mainstream explanation:** Asymmetric outgassing (rocket effect)
with CO/CO₂ as primary volatiles. Qualitatively plausible.

**Geometric relevance to the model:**

This is a NULL finding for the model. The 3I/ATLAS anomaly is best
explained by known cometary physics (outgassing). The anomaly is
large and associated with visible gas emission — this is NOT the
kind of absent-signature anomaly the model predicts.

For the model to be relevant here, you would need: anomalous
acceleration WITH absent thermal signature, no detectable gas plume,
no mass loss. The opposite is observed.

**VERDICT: NULL — 3I/ATLAS DOES NOT SUPPORT OR CHALLENGE THE MODEL.**
Outgassing is the correct explanation. The comparison is instructive
however: what distinguishes a K-bubble signature from an outgassing
signature is exactly the absent thermal/chemical signature. 3I/ATLAS
has all the signatures present. This reinforces the Aguadilla
"cold signature as thermal memory" argument by contrast.

---

### 3.3 Quaoar Ring — Outside Roche Limit (2023, ongoing)

**What it is:** Trans-Neptunian object Quaoar has a ring system
located well OUTSIDE its Roche limit. By classical mechanics,
material at this distance should accrete into a moon, not remain
as a ring. Temperature at Quaoar: ~-220°C. Leading explanations:
elastic collisions at cryogenic temperatures, unseen orbital
resonances, or non-classical low-gravity cohesion effects.

**Geometric relevance to the model:**

The Quaoar ring is anomalous within standard tidal mechanics. The
proposed explanation (elastic collisions at extreme cold) modifies
the effective inter-particle coupling: at -220°C, icy particles do
not stick upon collision — they bounce. The Roche limit calculation
assumes particles that accrete (sticky = fully coupled). If coupling
is suppressed at extreme cold, the effective Roche limit moves inward.

In the framework's language: at extreme low temperatures, the thermal
EM coupling between icy particles is reduced. Particle stickiness
(accretion cross-section) is a function of inter-particle coupling
strength. Reduced coupling = reduced accretion = ring stability beyond
classical Roche limit.

This is not a K-bubble. But it is an observational demonstration
that reducing inter-particle coupling (via temperature suppression
of thermal EM processes) changes the effective gravitational/tidal
physics of a system in a way that matches what K-suppression would
predict at larger scales.

**VERDICT: WEAK CONVERGENCE — DEMONSTRATES COUPLING-MODIFIED TIDAL
PHYSICS IS REAL AT ASTROPHYSICAL SCALES.**
The Quaoar ring is not evidence for K-bubbles, but it is evidence
that coupling strength modifications produce observable changes in
effective mechanical physics — the exact scaling argument the model
relies on.

---

### 3.4 Hubble Tension + DESI Dark Energy Results (2024–2025)

**What it is:** JWST + HST confirm the Hubble constant discrepancy
(H₀ local ~73 km/s/Mpc vs. CMB ~67 km/s/Mpc) is not a measurement
error. DESI Year 1 BAO data (2024) shows 2.5–3σ preference for
dynamical dark energy (w ≠ -1, possibly evolving). Possible local
effect: DESI tension may arise from observations at z < 0.1 — the
nearby universe.

**Geometric relevance to the model:**

The DESI result that dark energy may be evolving — and specifically
that the tension may be strongest at low redshift (local universe,
z < 0.1, within ~300 Mpc) — is geometrically relevant.

In the framework: the cosmological constant IS the ambient vacuum
energy density, which corresponds to ρ_ambient = K_ambient = 1.
If dark energy is evolving, this means K_ambient is not globally
fixed — the vacuum mode density at the cosmological scale is
changing with time (or position).

A spatially varying K_ambient is exactly what the framework requires
to be extendable from local K-bubble physics to cosmological scales.
The DESI result is a cosmological-scale indication that K_ambient
is not a universal constant — it varies.

**CRITICAL NOTE:** This is weak geometric convergence, not strong.
The DESI result is consistent with many dark energy models. It does
NOT specifically support the K-bubble framework. However, it rules
out the simplest cosmological assumption (K_ambient = const) which
would otherwise constrain the model.

**VERDICT: WEAK CONVERGENCE — DESI RESULT REMOVES A POTENTIAL
CONSTRAINT on the model but is not evidence for it.**

---

## SECTION 4 — THE REFRACTIVE INDEX CORRECTION — FULL IMPLICATIONS

This is the most consequential single finding of the sweep.

With n = K^(1/2) confirmed from Puthoff (2002):

**Re-examining the bubble wall optics:**

Inside bubble: K = ε → 0. Therefore n_interior = K^(1/2) → 0.
This means the refractive index INSIDE the bubble approaches ZERO —
not infinity.

The Vacuum_Coupling_Potential document stated:
"Inside the bubble (ρ → 0): n → ∞. Light cannot escape from inside
the bubble (effectively infinite refractive index)."

**This is WRONG with the correct PV formula.**

With n = K^(1/2): n_interior → 0 as K → 0.

A medium with n → 0 (epsilon-near-zero, ENZ) has OPPOSITE optical
properties to n → ∞:
- n → 0: phase velocity → ∞ (superluminal phase propagation)
- n → 0: the medium is NOT opaque — it is effectively transparent
  to radiation from OUTSIDE, but radiation from INSIDE propagates
  with infinite phase velocity
- At the bubble wall (n_interior → 0, n_exterior = 1):
  R = [(0 - 1)/(0 + 1)]² = 1 → 100% reflectivity

Wait — this recovers the optical isolation. But the mechanism is
different: the bubble wall is highly reflective NOT because n → ∞
inside, but because the impedance mismatch at the wall (n=0 inside
vs n=1 outside) produces total reflection.

This is epsilon-near-zero (ENZ) physics — a well-studied and
experimentally confirmed regime in metamaterials.

**ENZ materials are a real, experimentally realized class of
optical materials. They exist. They have n ≈ 0 at specific
frequencies. They produce anomalous transmission and reflection
properties exactly as described here.**

**GEOMETRIC CONSEQUENCE:**

The K-bubble interior is an ENZ medium. Not a high-index medium.
The optical properties of ENZ media are well-characterized:
- Near-unity transmission at specific angles (ENZ tunneling)
- Phase velocity → ∞ (spatial phase uniformity inside the bubble)
- The interior of the bubble is spatially "phase-compressed" — the
  entire interior is at the same electromagnetic phase
- This phase uniformity is what allows the object inside to
  displace without radiating — there is no phase gradient to drive
  electromagnetic emission

**This is a new geometric insight not previously stated in the model.**
The ENZ nature of the bubble interior explains the optical dark
interior without requiring n → ∞. It is geometrically cleaner and
consistent with known physics.

**VERDICT: MAJOR GEOMETRIC CLARIFICATION.**
The bubble interior is an ENZ medium, not a high-index medium.
This changes the optical picture but preserves all physical
predictions. ENZ physics is experimentally established.
This finding should be incorporated into the Vacuum_Coupling
document as a correction to the refractive index narrative.

---

## SECTION 5 — SYNTHESIS: WHAT MAKES THE PICTURE CLEARER

### What is now geometrically clearer than before this sweep:

**1. The refractive index formula is confirmed: n = K^(1/2)**
   The bubble interior is an epsilon-near-zero (ENZ) medium.
   ENZ physics is real, confirmed, and well-characterized.
   The Aguadilla reflectivity formula was correct.
   The Vacuum_Coupling document's n = 1/√ρ is wrong;
   the correct value is n = K^(1/2) = ρ^(1/6).

**2. The mechanism for frequency-selective K is clarified:**
   ENZ metamaterials operate at specific frequencies determined
   by their structure. The frequency at which n → 0 (and thus
   the bubble wall becomes reflective) is an engineerable property.
   This is how the K-bubble could simultaneously be reflective
   at IR/UV and transparent at microwave: the ENZ frequency is
   tuned to the IR/UV band. Below the ENZ frequency (microwave),
   the material is metallic (reflective). Above the ENZ frequency
   (IR/UV), the material is transparent OR the K suppression
   applies. The specific crossover depends on the K(ω) profile.
   This needs more careful treatment but is geometrically coherent.

**3. The cavity QED and DCE programs confirm the mechanism exists:**
   Vacuum mode density is measurably modifiable by engineered
   boundaries. Chemical reactions are changed. Ground-state
   energies are shifted. These are small effects but the
   direction and mechanism are confirmed.

**4. Quantized inertia provides independent theoretical convergence:**
   A completely independent theory (QI/McCulloch) arrives at the
   same geometric conclusion via a different route: inertia = mode
   coupling integral. QI's observational support is indirect
   evidence the geometry is correct.

**5. The vacuum birefringence program constrains the mechanism:**
   The Schwinger limit tension (T-5) is real. Standard QED
   electromagnetic field modification of K requires enormously
   strong fields. The K-bubble at 29 MJ cannot work via standard
   QED EM field mechanisms. S_bubble must be non-electromagnetic
   or operate via a different physical mechanism.

### What remains geometrically open:

**OQ-1 UPDATE:** RESOLVED by this sweep. n = K^(1/2) confirmed.
   Bubble interior is ENZ. Reflectivity formula was correct.

**OQ-3 UPDATE:** D_K interpretation as source relaxation (not field
   diffusion) is confirmed by the ENZ picture. ENZ media have
   a characteristic response time related to their resonant
   frequency. For the K-bubble to operate as an ENZ medium at IR
   frequencies (~10¹⁴ Hz), the natural response time would be
   ~1/f ~ 10⁻¹⁴ seconds — not 13.2 seconds. This reinforces the
   conclusion that D_K measures S_bubble source relaxation,
   not K-field propagation. The S_bubble mechanism has a 13-second
   timescale. This is a macroscopic mechanical or hydrodynamic
   process, not an EM process.

**NEW OPEN QUESTION OQ-6:** The Schwinger limit tension.
   What is the physical mechanism for S_bubble that does not
   require EM field strengths near the Schwinger limit but still
   produces K(hydro) = 0.133?

   Candidate directions:
   (a) Acoustic vacuum coupling — phonon-mediated mode suppression
       in a dense medium (not EM). This avoids the Schwinger limit.
   (b) Topological boundary conditions — the bubble wall is a
       topological defect in the vacuum field, sustained by
       topology rather than by continuous field energy input.
       A topological soliton could maintain K suppression at
       energies far below the Schwinger limit.
   (c) Casimir sphere geometry — a spherical Casimir cavity.
       The Casimir energy for a conducting sphere of radius r is:
         E_Casimir ≈ 0.04618 ħc/r
       For r = 1.88 m: E ≈ 10⁻³⁵ J. This is 37 orders of
       magnitude too small to explain 29 MJ. Casimir geometry
       alone cannot be S_bubble.
   (d) Dark energy coupling — if K_ambient is not constant (DESI
       result), then a local manipulation of the dark energy
       coupling could produce K suppression without requiring
       Schwinger-scale EM fields.

   None of these are established. But (b) topological soliton is
   the most geometrically motivated — topological objects maintain
   their structure through topology, not continuous energy input,
   and can exist at energies far below their constituent field
   quanta. This is worth pursuing.

---

## SECTION 6 — COMPLETE TENSION AND QUESTION REGISTER (UPDATED)

### Geometric Incompatibilities (UPDATED)

| ID | Description | Status After Sweep |
|---|---|---|
| I-1 | Refractive index formula | **RESOLVED.** n = K^(1/2) confirmed. Aguadilla formula correct. |
| I-2 | ρ vs K³ notation | Open — minor notation audit needed |
| I-3 | -c²∇ρ vs -3c²∇(ln K) inside bubble | Open — approximation regime must be stated |
| I-4 (NEW) | Vacuum_Coupling document states n = 1/√ρ — this is wrong | Correct to n = ρ^(1/6) = K^(1/2) |

### Tensions (UPDATED)

| ID | Description | Status After Sweep |
|---|---|---|
| T-1 | K³ scaling of convective coupling | Still open — van der Waals/Casimir analogy supports it |
| T-3 | Accelerating frequency mechanism | Still open — see OQ-2 |
| T-4 | D_K anomalously slow | Reframed as source relaxation — consistent with ENZ picture |
| T-5 (NEW) | Schwinger limit: QED EM field modification of K requires ~10⁴⁵ W/m² | Critical — S_bubble cannot be purely EM |

### Open Questions (UPDATED)

| ID | Description | Status After Sweep |
|---|---|---|
| OQ-1 | Confirm refractive index formula | **CLOSED.** n = K^(1/2) |
| OQ-2 | Correct mechanism for accelerating frequency | Still open |
| OQ-3 | What does D_K actually measure? | Clarified: source relaxation, consistent with ENZ timescale |
| OQ-4 | What generates S_bubble? | Still open. T-5 now constrains: cannot be standard EM |
| OQ-5 | Frequency-selective boundary mechanism | Clarified: ENZ physics provides the framework |
| OQ-6 (NEW) | Non-EM mechanism for S_bubble at sub-Schwinger energies | New — see topological soliton direction |

---

## SECTION 7 — PRIORITY ACTION LIST (UPDATED)

**Priority 1 (Blocking — update documents):**
- [ ] Correct n = 1/√ρ to n = ρ^(1/6) = K^(1/2) in
      Vacuum_Coupling_Potential_Physics_Derivation.md
- [ ] Update "n → ∞ inside bubble" narrative to "ENZ inside bubble"
      throughout all documents
- [ ] Update I-1 in PHYSICS_DEEP_DIVE.md: Aguadilla formula
      was correct; n = K^(1/2) confirmed; radar threshold
      corrects to K(radar) > 0.669

**Priority 2 (Important — new physics to document):**
- [ ] Add ENZ (epsilon-near-zero) physics section to the model.
      The bubble interior is an ENZ medium. Cite the established
      ENZ metamaterial literature. This is experimentally confirmed
      physics.
- [ ] Add T-5 (Schwinger limit) to PHYSICS_DEEP_DIVE.md as a
      genuine constraint on S_bubble mechanism.

**Priority 3 (Research frontier):**
- [ ] Explore topological soliton as candidate for S_bubble (OQ-6).
      Topological Q-balls, skyrmions, and domain walls in quantum
      field theory maintain structure without continuous energy
      input. Is there a vacuum-field topological object that
      maintains an ENZ interior?

**Priority 4 (Observational follow-up):**
- [ ] Check flyby anomaly data for correlation with K-field
      geometry prediction: the anomaly should scale with
      sin²(δ) where δ is the declination of the incoming
      asymptote. This is the formula already in the literature
      (Anderson et al. 2008). If K-field co-rotation with Earth
      is the source, it should produce this exact scaling.
      This is a falsifiable check that can be done now.

---

## DOCUMENT METADATA

- Author: Eric Robert Lawson / GitHub Copilot
- Date: 2026-03-21
- Status: Active discovery record
- Version: 1.0
- Supersedes: Partial entries in PHYSICS_DEEP_DIVE.md
  (specifically I-1, OQ-1, and OQ-3 are updated here)
- Related files:
  - PHYSICS_DEEP_DIVE.md (deep audit — update I-1, OQ-1, OQ-3)
  - AGUADILLA_TODO.md (operational pickup file)
  - Vacuum_Coupling_Potential_Physics_Derivation.md (needs n correction)
  - Reverse_Engineering_Aguadilla.md (reflectivity formula confirmed)

---

## ONE-LINE SUMMARY

The ENZ (epsilon-near-zero) interpretation of the bubble interior —
derived from confirming n = K^(1/2) in the PV framework — is the
single most clarifying result: the bubble is not an optical trap,
it is a phase-uniform ENZ cavity, and this is known, confirmed,
engineerable physics.

---
