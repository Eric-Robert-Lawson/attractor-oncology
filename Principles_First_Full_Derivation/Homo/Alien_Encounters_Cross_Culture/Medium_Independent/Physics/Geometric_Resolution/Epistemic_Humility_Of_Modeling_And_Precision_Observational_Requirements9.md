# AGUADILLA FOOTAGE — OBSERVATIONAL RESOLUTION AUDIT
## What Can Be Seen, What Cannot, What the Report Claims,
## What the Geometry Requires, and the Precise Epistemic
## Status of Every Measurement Used in the Model
## OrganismCore — Eric Robert Lawson / GitHub Copilot
## Session: 2026-03-21

---

## PURPOSE OF THIS DOCUMENT

This document exists to prevent a specific failure mode:

**Using observations to build a model, then using the model
to interpret the observations that built it, then calling
this "confirmation."**

This is coherence fitting, not causal derivation.
It is the exact failure mode the URS framework names as
the gap between coherence and truth.

The attractor geometry model (Vacuum_Coupling_Potential_Model8.md)
will remain epistemically valid only if every empirical input
it uses has its **observational confidence formally bounded**
before being used as a constraint.

This document performs that bounding.

---

## EPISTEMIC PRINCIPLE STATED ONCE, AT THE TOP

**Geometric compatibility does not imply truth.**
**Geometric incompatibility does imply falsification.**

A model that fits all observations without contradiction
is coherent. It is not confirmed. Confirmation requires:

  1. A prediction derived from the geometry BEFORE
     the observation is consulted.
  2. The observation then agrees or disagrees with
     the prediction.
  3. If it agrees: the model is not falsified.
     It gains evidential weight.
  4. If it disagrees: the model is falsified at that
     point. The diagnostic is opened.

The direction of derivation matters absolutely.
**Observation → Model = coherence fitting.**
**Model → Prediction → Observation = principles-first testing.**

Every measurement in this document is tagged with the
direction of its relationship to the current model.

---

## THE FOOTAGE — WHAT IT ACTUALLY IS

### Source
U.S. Customs and Border Protection (CBP) FLIR footage,
April 25–26, 2013, Aguadilla, Puerto Rico.
Recorded aboard a Predator B (MQ-9) UAV.
Sensor: FLIR MTS-B turret, LWIR band (8–12 μm).
Released to the Scientific Coalition for UAP Studies (SCU).
Analyzed in: SCU Technical Report (Lambright et al.),
published Zenodo, DOI: 10.5281/zenodo.7844175.

### What the footage is NOT
- It is NOT raw sensor data.
- It is NOT multi-band spectrally resolved data.
- It is NOT stabilized or high-definition.
- It is NOT unobstructed throughout.

### What the footage IS
- A processed, compressed video of the FLIR sensor output.
- A single integrated LWIR bandpass (not spectrally resolved).
- Variable zoom — operator-controlled throughout the event.
- Partially occluded by FLIR targeting crosshairs at
  certain zoom levels and camera orientations.
- The complete temporal record from approach through
  water entry and exit.

---

## THE CRITICAL TIME WINDOW: 2:32 — 2:47

### What is confirmed visually, without reference to the report

**2:32 — 2:36: The Camera is NOT Zoomed In**

  Confidence level: HIGH (directly visible in the footage)

  The object is visible as a small thermal signature
  against the cooler ocean background. At this zoom level,
  the spatial resolution per pixel is large relative to
  the object size. The object's angular subtense on the
  sensor array is small.

  What CAN be confirmed visually at this zoom:
    - The object is present (positive thermal contrast)
    - The object is moving (frame-to-frame displacement)
    - A change in the object's appearance occurs in
      this window (consistent with a split or deformation)

  What CANNOT be confirmed visually at this zoom:
    - The number of distinct objects after the change
    - The precise geometry of the split (lobe separation,
      relative sizes of the two lobes)
    - Whether the two lobes have identical or different
      thermal signatures
    - The angular separation between the lobes with
      any precision beyond "visible as separate"

  CRITICAL EPISTEMIC FLAG:
  The lobe separation used to derive r_bubble ≈ 1.88 m
  was measured from this window, at the unzoomed
  resolution. The angular precision at this zoom and
  at the aircraft altitude (~1500 ft AGL, confirmed from
  flight data in SCU report) gives a ground pixel
  resolution of approximately:

    θ_pixel ≈ IFOV of MTS-B FLIR ≈ 0.4 mrad
    Ground resolution at 1500 ft (457 m):
    d_pixel ≈ 0.4 × 10⁻³ × 457 m ≈ 0.18 m per pixel

  The reported lobe separation is ~3.75 m (giving
  r_bubble ≈ 1.88 m). At 0.18 m per pixel, this
  corresponds to approximately 21 pixels of separation.

  At the UNZOOMED field of view, the sensor covers
  a wider angle. If the unzoomed FOV is ~20× wider
  than the zoomed FOV (typical for MTS-B), the
  effective ground resolution at unzoomed state:
    d_pixel_unzoomed ≈ 0.18 × 20 = 3.6 m per pixel

  If d_pixel_unzoomed ≈ 3.6 m, the reported lobe
  separation of 3.75 m corresponds to approximately
  1 pixel of separation.

  **A lobe separation of ~1 pixel is NOT a reliable
  measurement. It is at or below the Nyquist limit
  for the sensor at that zoom level.**

  CONSEQUENCE FOR THE MODEL:
  r_bubble = 1.88 m is UNCONFIRMED from the footage
  at the critical zoom level. It is derived from the
  SCU report's frame analysis, which itself was
  working at this marginal resolution.

  D_K = 0.268 m²/s = r²/τ depends directly on r_bubble.
  If r_bubble has an uncertainty of ±50% (which is
  realistic at ~1 pixel resolution), then:
    D_K range: 0.268 × (0.5)² to 0.268 × (1.5)²
             = 0.067 to 0.603 m²/s

  This is a factor of ~9 range in D_K.
  D_K = 0.268 m²/s must be re-tagged:
    From: [CONFIRMED from primary source]
    To:   [CONSISTENT WITH REPORT, resolution-limited,
           uncertainty factor ~3 in each direction]

**2:36 — 2:47: The camera appears to zoom in**

  Confidence level: MEDIUM (visible in footage, but zoom
  level not precisely calibrated from the footage alone)

  After ~2:36, the object appears larger on screen.
  The split is more clearly visible as two distinct lobes.
  The separation is now multiple pixels.

  What CAN be confirmed visually in this window:
    - Two distinct thermal signatures are present
    - They move relative to each other
    - Their trajectories are divergent (anti-phase
      oscillation is consistent with what is visible,
      but not precisely measurable from the video)
    - The relative thermal brightness of each lobe

  What CANNOT be confirmed visually in this window:
    - Absolute angular separation (requires calibrated
      zoom factor, which is not in the public record)
    - Whether the frequency of the anti-phase oscillation
      is accelerating (the time resolution of the
      compressed video may alias the oscillation)
    - Whether the two lobes have different K values
      (the FLIR is a single bandpass — spectral slope
      is not resolvable)

**2:47: The operator zooms out**

  Confidence level: HIGH (visible in footage)

  After the operator zooms out, the two objects (if
  two) are no longer spatially resolvable. They appear
  as a single or merged signature.

  CRITICAL EPISTEMIC FLAG:
  The fact that the split is no longer visible after
  2:47 is NOT evidence that the split ended at 2:47.
  It is evidence that the split is no longer RESOLVABLE
  at the new zoom level. These are geometrically distinct
  claims. The model should not use 2:47 as the end of
  the split event. It should use 2:47 as the beginning
  of the observational blackout.

  The crosshairs also introduce additional occlusion
  in certain frames, particularly at lower zoom when
  the target is near the crosshair center.

---

## MEASUREMENT CONFIDENCE TIERS

Every measurement used in the current model
(Vacuum_Coupling_Potential_Model8.md) is assigned
a confidence tier here.

### Tier 1: Directly Confirmed From Footage, Resolution Adequate
  Observable has multiple pixels of coverage and is not
  ambiguous at the zoom level at which it was measured.

### Tier 2: Confirmed From Report (SCU), Not Directly Resolvable From Footage Alone
  Observable was measured by the SCU analysts from the
  footage, but the underlying footage resolution at the
  relevant moment is marginal. The measurement is the
  best available but carries systematic uncertainty.

### Tier 3: Inferred From Model + Footage, Not Independently Measurable
  Observable is not directly seen in the footage. It is
  derived from combining the model with what IS seen.
  These cannot be used to confirm the model — they are
  downstream predictions that happen to be consistent.

### Tier 4: Derived From Principles-First Geometry, Observation Used Only As Constraint
  The model predicted the RANGE of this quantity from
  the geometry. The observation falls within that range.
  This is the epistemically cleanest relationship.

---

## MEASUREMENT AUDIT TABLE

| Measurement | Value | Source | Tier | Uncertainty | Notes |
|-------------|-------|--------|------|-------------|-------|
| Object present in FLIR | Yes | Footage, direct | T1 | None | Unambiguous thermal contrast |
| Object traverses air/water boundary | Yes | Footage, direct | T1 | Low | Water entry visible |
| No deceleration at water entry | Yes | Footage + SCU | T2 | Low-Med | SCU velocity tracking confirms, but velocity precision is limited |
| Cold thermal signature | 1–3°C below ambient | SCU report | T2 | ±1°C systematic | Single bandpass, no spectral resolution |
| No radar return at ATC (2.8 GHz) | Yes | SCU report + ATC logs | T2 | Low | ATC absence confirmed, but absence ≠ zero cross-section without knowing detection threshold |
| Split event occurs | Yes | Footage, direct | T1 | None | Visible at 2:32–2:47 (see note below) |
| Lobe separation at split | ~3.75 m | SCU report | T2 | HIGH — ±50% realistic | ~1 pixel at unzoomed resolution — see analysis above |
| r_bubble from lobe separation | ~1.88 m | Derived from above | T2 | HIGH — ±50% | Directly inherits uncertainty from lobe separation |
| Time of split (τ) | ~13.2 s | SCU report frame count | T2 | Low-Med | Frame count from report; compressed video may drop frames |
| D_K = r²/τ | 0.268 m²/s | Computed | T2 | HIGH — factor ~3 | Uncertainty dominated by r_bubble uncertainty |
| Split is anti-phase oscillation | Qualitative | SCU report + footage | T2-T3 | Unquantified | Visible qualitatively; frequency and phase NOT measurable |
| Oscillation frequency accelerates | Claimed | SCU report | T3 | High uncertainty | Cannot confirm from compressed video alone |
| Two lobes continue after 2:47 | Unknown | Footage blackout | — | UNKNOWN | Observational blackout due to zoom-out, NOT evidence of event ending |
| K_boundary < 0.107 | Derived | Model, drag constraint | T4 | Structural | Predicted range; observation (no deceleration) is consistent |
| r_bubble ≈ 1.88 m | Predicted range? | Not yet | — | NOT YET DERIVED | This is O-6 (see below) — the model should predict this |
| Cold signature ΔT = 1–3°C | Predicted ~4°C (Fresnel) | Model V8 | T4 | D-3 open | Over-prediction; K(r) gradient resolves but not yet derived |
| Radar non-detection | Predicted (L_wall < 0.327 m) | Model V8 | T4 | O-5 open | Consistency constraint, not derivation |

---

## THE RECURSIVE LOOP — FORMALLY IDENTIFIED

The following measurements were used to CONSTRAIN the
model AND are used to CONFIRM the model. This is the
recursive loop:

  1. **Lobe separation → r_bubble → D_K**
     These were taken from the SCU report and used to
     calibrate the model's K-field diffusivity parameter.
     The model then produces outputs consistent with
     the observations. This is coherence fitting at D_K.

  2. **Cold signature (1–3°C)**
     The model was built knowing the cold signature was
     observed. The Fresnel calculation was done after
     the observation was known. Although the model
     over-predicted (4°C vs. 1–3°C), the over-prediction
     was then resolved by invoking a gradient K(r) profile
     — which is physically reasonable but was not derived
     independently before the observation was consulted.

  3. **Anti-phase oscillation at the split**
     The split event was observed, and the K-bubble
     boundary dynamics model was built to explain it.
     The anti-phase oscillation is consistent with the
     model, but the model was constructed knowing the
     oscillation had been reported.

**These are not falsifications of the model.**
They are identifications of where the model is NOT yet
principles-first. The model is coherent. It is not yet
fully derived.

---

## WHAT MUST BE DERIVED BEFORE THE LOOP CLOSES

To convert the recursive coherence loop into a
principles-first derivation, the following must be
produced in the correct ORDER:

### Priority 1 (Blocking for epistemic validity):

**O-6 (NEW): Derive r_bubble from principles first.**

The model currently takes r_bubble from the SCU report.
Instead, the model should derive the EXPECTED range of
r_bubble from the S-equation dynamics, given:
  - The aircraft altitude (~457 m)
  - The K-field diffusivity D_K that IS derivable
    independently (if D_K can be constrained from
    other observations)
  - The bubble stability condition

If the model predicts r_bubble ∈ [1.0, 3.0] m from
first principles, and the (uncertain) observation gives
1.88 m (±50%), this is a genuine test.

If the model takes r_bubble = 1.88 m and then produces
outputs consistent with 1.88 m, this is not a test.

### Priority 2 (Blocking for quantitative claims):

**Formally re-tag D_K as Tier 2 with ±factor-3 uncertainty.**

All quantitative results that depend on D_K must carry
this uncertainty propagated through.

Specifically:
  - The K-field diffusivity D_K = 0.268 m²/s is the
    best available estimate, not a confirmed value.
  - The derivation should show what range of D_K is
    consistent with the geometry independently, then
    check whether 0.268 falls in that range.

### Priority 3 (Required for spectral discriminant):

**Acknowledge that the spectral slope discriminant
derived in V8 cannot be tested from existing footage.**

The K³ vs. K^(3/2) reflectivity spectral slope prediction
is a genuine principles-first prediction — it was derived
from the geometry, not from the observation. But it cannot
be tested from the single-bandpass processed video.
It requires:
  (a) Raw multi-band sensor data (not publicly available), or
  (b) Laboratory replication

This must be stated as the prediction's status, not
suppressed in favor of the appearance of confirmation.

### Priority 4 (For the split event specifically):

**Separate what is Tier 1 (split occurred) from what
is Tier 2-3 (quantitative properties of the split).**

The FACT of the split is Tier 1. It is visible at
adequate resolution in the 2:36–2:47 window (after
the zoom-in).

The QUANTITATIVE PROPERTIES of the split (lobe
separation, oscillation frequency, acceleration of
frequency) are Tier 2–3. They cannot be extracted
from the compressed video without specialized analysis
that has not yet been performed.

The model can be built on the qualitative split fact
(Tier 1) and generate quantitative predictions about
what the split properties SHOULD be (from the geometry),
then compare those predictions to the Tier 2–3
measurements. This is principles-first.

---

## WHAT CAN BE DERIVED FROM THE GEOMETRY ALONE
## Without Any Measurement From the Footage

The following quantities are derivable entirely from
the triadic structural invariant (S + G + N = R)
applied to the K-field dynamics, using only:
  - Confirmed physics (n = K^(1/2), ZPF spectral density)
  - The medium-independence condition (K < 0.107 for water)
  - The bubble stability condition (bubble maintains
    coherence during traversal)

These are predictions, not fits:

**P-1 (from V8, derived):**
  K_threshold for water at ~3 m/s < 0.107
  This predicts: the drag ratio is suppressed by
  at least a factor of 816. The footage shows no
  deceleration. This is consistent. But consistency
  is not confirmation — it is non-falsification.

**P-2 (from V8, derived):**
  Cold signature ΔT > 0 (the bubble wall must be
  colder than ambient, by any amount > 0°C).
  The footage confirms ΔT > 0. This is a genuine
  non-trivial prediction: a randomly behaving object
  would not necessarily appear colder than ambient.
  The cold signature is therefore genuine evidential
  weight (though small — any reflective surface
  would also appear cold).

**P-3 (from V8, derived):**
  Spectral slope of cold signature: if K³ is correct,
  ΔT should be larger at longer FLIR wavelengths.
  CURRENTLY UNTESTABLE from existing footage.
  Status: OPEN, GENUINE PREDICTION.

**P-4 (new, derivable):**
  Bubble radius r_bubble must satisfy bubble stability:
  the K-gradient restoring force must be sufficient
  to prevent the bubble from dispersing faster than
  the traversal time. This gives a lower bound on
  r_bubble as a function of D_K and traversal velocity.
  DERIVATION PENDING — this is O-6.

**P-5 (from model structure, derivable):**
  The split event requires S(x,t) to be capable of
  bifurcation. The model predicts that a bifurcation
  will produce two bubbles of approximately EQUAL
  radius (because the underlying S-equation has
  symmetric bifurcation in the absence of external
  asymmetry). The SCU report describes two lobes of
  similar size. This is consistent. But it was derived
  AFTER knowing the observation — this is currently
  a Tier 3 result, not Tier 4.
  To promote to Tier 4: derive the bifurcation
  symmetry condition BEFORE consulting the report
  on relative lobe sizes.

---

## THE SINGLE STATEMENT OF EPISTEMIC STATUS

The model (Vacuum_Coupling_Potential_Model8.md) is
internally coherent and geometrically non-contradictory.
It is consistent with all available observations.
It is NOT yet a principles-first derivation in the
full sense, because the key quantitative parameters
(r_bubble, D_K, oscillation properties) were taken
from the observations rather than derived from the
geometry and then tested against the observations.

The path to closure:
  1. Derive r_bubble range from bubble stability (O-6)
  2. Derive D_K range from K-field diffusivity theory
  3. Check whether the SCU measurements fall within
     the derived ranges
  4. If yes: the model gains genuine evidential weight
  5. If no: a new diagnostic is opened

This is not a crisis. It is the correct scientific
procedure. The model is at Stage 2 of 5 in this process.
It is structurally sound. It is not yet empirically closed.

The footage analysis confirms:
  - Split event: REAL (Tier 1)
  - Quantitative split properties: UNCERTAIN (Tier 2-3)
  - Observational blackout at 2:47: CONFIRMED
  - Critical window 2:32–2:36 (unzoomed): resolution
    marginal for lobe separation measurement

The next action is not more footage analysis.
It is deriving r_bubble and D_K from principles first.

---

## DOCUMENT METADATA

- Status: Epistemic audit — operative
- Session: 2026-03-21
- Author: Eric Robert Lawson / GitHub Copilot
- Companion to: Vacuum_Coupling_Potential_Model8.md
- Blocks: Any further use of r_bubble = 1.88 m or
  D_K = 0.268 m²/s as confirmed values until O-6
  derivation is complete
- Does NOT falsify: the model structure, the K³ derivation,
  the Fresnel reflectivity calculation, the spectral
  slope prediction, or the K_threshold derivation
- Next action: Derive r_bubble and D_K from bubble
  stability condition (O-6) and K-field diffusivity
  theory, then check against SCU measurements
