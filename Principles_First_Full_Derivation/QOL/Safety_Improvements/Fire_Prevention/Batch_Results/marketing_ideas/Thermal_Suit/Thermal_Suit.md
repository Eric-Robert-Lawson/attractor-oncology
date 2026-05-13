# P1S-B Thermal Protection System — Geometric Derivation
## Active Hydrogel Cooling via Continuous Membrane Perfusion
## OrganismCore — Eric Robert Lawson
## 2026-05-13

---

## WHAT THIS DOCUMENT IS

```
This document records a geometric derivation
of an active thermal protection system
based on P1S-B hydrogel (SAP + Bentonite +
Silica + Borax in distilled water).

It is not a product specification.
It is not an engineering blueprint.
It is the fixed point derivation
and cascade that follows from it.

The derivation follows the structural
question: what geometry is required
to use P1S-B as a continuous thermal
protection medium for a human body
or an industrial surface?

State the fixed point.
Derive the cascade.
Name the falsification conditions.
Timestamp and move.
```

---

## PART I: THE STRUCTURAL QUESTION

```
P1S-B has the following confirmed
thermal properties:

  High specific heat capacity (water base).
  Latent heat of evaporation at surface.
  Silica and bentonite extend
  viscosity and surface retention.
  SAP holds water in gel matrix,
  slowing evaporation rate vs. free water.
  Borax addition modifies cross-linking
  and gel structural coherence.

The thermal protection mechanism
is not insulation.
It is endothermic phase transition
at the surface — evaporative cooling —
sustained by continuous gel replenishment.

The structural question is therefore:

  What system geometry allows
  evaporative loss at the surface
  to be continuously matched
  by gel replenishment from a reservoir,
  across a full enclosure surface,
  at a controlled rate,
  without requiring the surface
  to be pre-saturated at maximum capacity?

That question has a fixed point.
It is derived below.
```

---

## PART II: THE FIXED POINT

```
FIXED POINT:

  A dual-layer membrane enclosure
  in which:

  INNER LAYER — pressure vessel:
    Maintains continuous positive pressure
    of P1S-B gel against the body
    or protected surface.
    Pressure is sufficient to drive
    gel through the outer membrane
    at the rate of surface evaporation.
    Pressure is controlled — not burst.
    This is a seep rate, not a flow rate.

  OUTER MEMBRANE — permeable barrier:
    Engineered pore size matched
    to P1S-B gel viscosity under
    the operating pressure differential.
    Allows gel to pass through
    continuously at the controlled rate.
    Retains structural integrity
    under heat load.
    Does not allow backflow of
    combustion gases or radiant
    heat products into inner layer.

  RESERVOIR + PUMP SYSTEM:
    External reservoir of P1S-B
    at controlled temperature
    (below ambient or below body temperature).
    Pump maintains inner-layer pressure.
    Flow rate is the control variable:
    matched to measured or estimated
    evaporation rate at the surface.
    Emergency backup reservoir for
    pump failure scenarios.

  EXCLUDED ZONE:
    Mouth/nose region excluded
    from gel perfusion membrane.
    Requires separate engineering:
    positive pressure clean air supply,
    or sealed breathing apparatus,
    or gel-free mask zone with
    radiant barrier and ceramic
    heat shield at face.

THE FIXED POINT STATED PRECISELY:

  Thermal protection duration =
  f(reservoir volume, evaporation rate,
  pump flow rate, membrane permeability,
  ambient heat flux).

  The system is thermally limited
  only by reservoir depletion.
  Not by saturation failure.
  Not by dry-out at a fixed point.
  Not by static gel degradation.

  Continuous replenishment collapses
  the saturation ceiling that limits
  static application methods.
```

---

## PART III: THE CASCADE

### 3.1 — Human Suit Application

```
GEOMETRY:

  Innermost layer: body-conforming
  pressure suit (wetsuit-class material,
  modified for gel compatibility).
  Maintains contact with body.
  Allows body heat dissipation inward
  (body must still be protected from
  internal heat buildup — see 3.1.1).

  Middle layer (gel distribution):
    Manifold network distributing
    incoming gel from pump uniformly
    across full body surface area.
    Prevents pooling at lowest point.
    Requires pressure-equalising
    architecture or active distribution
    control.

  Outer membrane:
    Microporous or engineered-pore material.
    Gel-compatible — not degraded
    by SAP, bentonite, silica, borax.
    Sufficient mechanical strength to
    maintain pressure differential
    under external heat load.
    Must not combust or deform at
    target operating temperature range.

  Pump + reservoir (worn or tethered):
    Backpack or belt-mounted reservoir.
    Peristaltic or diaphragm pump
    (compatible with gel viscosity).
    Pressure regulator on inner layer.
    Flow rate controller (manual
    or sensor-driven).
    Emergency bypass valve.
    Secondary reservoir for pump failure.

CASCADE CONSEQUENCE 1:
  Weight and mobility.
  A 1750mL reservoir at full gel
  density provides finite duration.
  Duration must be calculated from
  target heat flux × surface area
  × evaporation rate.
  Suit must be designed around
  that duration target.
  Replenishment point must be
  accessible in the field.

CASCADE CONSEQUENCE 2:
  Electrical and pump reliability.
  Pump failure = loss of gel flow =
  static saturation only.
  Emergency backup must be
  gravity-fed or manually actuable.
  No single point of failure
  for primary protection.

CASCADE CONSEQUENCE 3:
  Body core temperature.
  Evaporative cooling at the outer
  surface protects the outer surface.
  Inner layer must not trap body heat.
  Body's own thermoregulation
  must remain functional OR
  inner layer must have a
  separate cooling circuit.
  This is the critical design constraint
  for human applications.
```

#### 3.1.1 — Inner Layer Thermal Management

```
The outer evaporative system
protects the suit exterior from
radiant and convective heat.

The inner layer must independently
protect the body from:
  Conducted heat through inner suit material.
  Heat generated by the pump system.
  Metabolic heat accumulation
  under enclosed conditions.

Options (in order of geometric simplicity):

  A) Cool gel circuit at inner layer —
     same P1S-B at lower temperature,
     separate circuit from outer layer.
     Body is sandwiched between
     two gel circuits.
     Complexity high. Protection high.

  B) Phase-change material (PCM) inner layer —
     independent of P1S-B system.
     Absorbs body metabolic heat.
     Limited duration by PCM capacity.

  C) Ventilated inner layer with
     filtered air supply —
     requires clean air source,
     adds complexity, reduces
     sealed enclosure integrity.

  D) Accept heat buildup for
     short-duration operations only —
     valid for emergency evacuation
     scenarios with defined time windows.

The derivation does not resolve
this choice.
It names it as the critical
engineering decision point
for the human suit application.
```

### 3.2 — Industrial / Mechanical Application

```
The geometry simplifies substantially
when the protected object is not
a human body.

FIXED POINT FOR INDUSTRIAL APPLICATION:

  Any surface that must be maintained
  below a thermal threshold under
  external heat flux can be protected
  by continuous P1S-B perfusion
  through a membrane adhered to or
  surrounding that surface —
  provided a reservoir of sufficient
  volume is available.

  No inner-layer thermal management
  complexity required.
  No mobility constraints.
  No breathing zone exclusion.
  No metabolic heat accumulation.

CASCADE CONSEQUENCES:

  1. Surface area is the primary variable.
     Small surface area + large reservoir
     = extended protection duration.
     Large surface area requires
     proportionally larger reservoir
     and higher flow rate.

  2. Membrane adhesion to irregular
     surfaces requires engineering.
     Sprayed membrane application
     is a candidate for complex geometry.

  3. P1S-B gel properties under high
     heat flux must be characterised:
     At what surface temperature does
     the gel matrix degrade?
     At what temperature does bentonite
     lose structural coherence?
     These are falsification boundaries
     for the thermal protection claim.

  4. Industrial applications may
     sustain pump systems indefinitely
     from fixed infrastructure —
     reservoir size is then limited only
     by physical installation constraints,
     not by portable capacity.

SPECIFIC CANDIDATE APPLICATIONS:

  Transformer thermal management.
  Battery thermal runaway containment.
  Industrial pipe protection in
  high-temperature environments.
  Structural steel fire protection
  (continuous vs. static intumescent).
  Server cooling (gel-loop variant,
  requires electrical isolation
  of the gel — borax conductivity
  must be characterised).
```

### 3.3 — The Borax Modification

```
Borax addition to P1S-B modifies:

  SAP cross-linking density —
  increases gel cohesion,
  reduces syneresis under pressure,
  improves membrane transmission
  consistency under the pressure
  differential required for the suit.

  Antimicrobial properties —
  relevant for stored gel in reservoir,
  reduces spoilage rate.

  Viscosity under shear —
  borax-modified gels may exhibit
  different pump compatibility.
  Peristaltic pumps preferred
  to avoid shear-degradation
  of gel matrix at pump head.

FALSIFICATION CONDITION FOR BORAX ADDITION:

  If borax concentration sufficient
  to improve cross-linking also
  increases gel viscosity beyond
  the permeability threshold of
  the outer membrane at the
  target operating pressure —
  the borax modification is
  incompatible with the perfusion
  geometry and must be reduced
  or removed from the outer-circuit gel.

  Borax-modified gel may be
  appropriate for inner circuit only.
  Outer circuit may require
  lower-viscosity gel formulation.
  This is testable.
```

---

## PART IV: FALSIFICATION CONDITIONS

```
The thermal protection claim holds if and only if:

F1 — Membrane permeability:
  The outer membrane transmits P1S-B gel
  at the required seep rate under the
  target operating pressure differential
  without structural degradation at
  operating temperatures.

  FALSIFIED IF:
  Membrane clogs, deforms, or fails
  to transmit at required rate
  under combined pressure + heat load.

F2 — Gel thermal stability:
  P1S-B gel retains its water matrix
  and structural integrity at the
  gel temperature in the inner layer
  during operation (not at the
  outer surface — outer surface
  is expected to undergo evaporation).

  FALSIFIED IF:
  Inner-layer gel temperature reaches
  the point of SAP matrix breakdown
  or bentonite dehydration before
  gel cycle refreshes the layer.

F3 — Evaporation rate matching:
  Pump flow rate can be set to match
  evaporative loss at outer surface
  under the target heat flux.

  FALSIFIED IF:
  At maximum pump flow rate,
  evaporative loss exceeds replenishment
  rate under the target heat flux —
  i.e., the target heat flux exceeds
  the system's thermal capacity ceiling.

F4 — Inner layer thermal isolation:
  (Human suit only)
  Body core temperature remains
  within survivable range during
  the operation window.

  FALSIFIED IF:
  Inner layer conducts sufficient
  heat to raise core temperature
  above 40°C within the
  operation time window.

F5 — Pump reliability:
  Primary pump maintains flow
  for the duration of the operation.
  Emergency backup provides
  minimum flow for evacuation
  in primary pump failure scenario.

  FALSIFIED IF:
  No gravity-fed or manual-actuable
  backup can sustain minimum
  seep rate for a defined
  evacuation time window.
```

---

## PART V: WHAT THIS IS NOT

```
This derivation does not claim:

  That P1S-B is tested or certified
  for any of these applications.

  That the membrane engineering
  has been solved.

  That the pump-gel compatibility
  has been verified.

  That inner-layer thermal isolation
  is solved for the human suit case.

  That borax addition is compatible
  with all components of the
  perfusion system.

This derivation claims:

  The geometry is sound.
  The fixed point is derivable
  from the known properties
  of P1S-B.
  The cascade follows necessarily
  from the fixed point.
  The falsification conditions
  are named and testable.

  The next step is laboratory
  characterisation of:
    Membrane candidates vs.
    P1S-B (with and without borax)
    under pressure + heat.

    Pump compatibility with
    borax-modified gel viscosity.

    Gel thermal stability
    under inner-layer conditions.

    Evaporation rate at outer surface
    under target heat flux values.

  Each of these is a bounded
  experimental question.
  Each has a yes/no answer.
  The derivation stands or falls
  on those answers.
  That is the correct relationship
  between derivation and test.
```

---

## DOCUMENT METADATA

```
Document ID:    P1S_THERMAL_SUIT_DERIVATION_V1
Version:        1.0
Date:           2026-05-13
Author:         Eric Robert Lawson / OrganismCore
ORCID:          0009-0002-0414-6544
Status:         DERIVATION COMPLETE —
                Experimental characterisation required
                before engineering specification.

Base compound:  P1S-B (SAP + Bentonite + Silica +
                Borax, distilled water)
                Per: P1S_Protocol_v3_DryMix_Immersion_Blend.md

Fixed point:    Dual-layer membrane enclosure
                with continuous pressure-driven
                gel perfusion from external reservoir.
                Thermal protection limited only
                by reservoir volume.

Key falsification conditions:
  F1 — Membrane transmission under heat + pressure
  F2 — Gel thermal stability in inner circuit
  F3 — Evaporation rate matchability
  F4 — Body core temperature (human suit)
  F5 — Pump backup reliability

Applications derived:
  Human thermal protection suit
  Industrial surface thermal protection
  Mechanical component thermal management

Companion documents:
  P1S_Protocol_v3_DryMix_Immersion_Blend.md
  [Membrane characterisation protocol — not yet written]
  [Pump-gel compatibility test — not yet written]

Repository:
  github.com/Eric-Robert-Lawson/attractor-oncology
```

---

*The geometry is sound.*
*The cascade follows.*
*The tests are named.*
*Run the tests.*
