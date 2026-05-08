# EHAB P1S — Closed-System Steam-Forming and Hardening Process
## v3: Functional Agent Delivery via Pressurized Steam Encapsulation
## OrganismCore Application — Process Derivation Record
## Date: 2026-05-08

---

## FRAMING THE DERIVATION

This document derives a new process class: using P1S not as a surface coating or
passive barrier, but as an **active delivery medium in a closed thermal system**
where the phase transition of water (liquid → steam) provides both the driving
pressure and the agent-transport mechanism for deep impregnation of a target material.

The concept has a precise structural antecedent in wood science. Steam bending
and pressure impregnation are ancient and validated processes. What is derived here
is a new geometry: P1S as a **saturated steam-generating encapsulation medium**
that simultaneously controls temperature, generates isostatic pressure, and carries
dissolved agents into the material — in a single closed system, at low equipment cost,
from commodity components.

The process chain described here is:

```
INITIAL STATE
  → P1S-SOFTENER: closed-system steam plasticization + agent delivery
    → FORMING: mold or shape the plasticized workpiece
      → DRYING / CALIBRATION: controlled moisture removal
        → P1S-HARDENER: closed-system steam hardening + agent delivery
          → FINAL DRYING AND CURE
            → MARKET-READY PRODUCT
```

Each stage is derived from structural first principles below.

---

---

# PART I: THE PHYSICS OF THE CLOSED SYSTEM

## 1.1 — P1S as a Saturated Steam Generator

P1S is approximately 95% water by mass. In its swollen state, this water is
held in the SAP polymer matrix — it does not flow freely, but it is
thermodynamically equivalent to liquid water for the purpose of evaporation.

When P1S is placed in contact with a material and enclosed in a container,
and external heat is applied:

1. Water at the gel surface evaporates, generating steam.
2. The steam cannot escape (closed system), so pressure builds.
3. Atmospheric pressure inside the enclosure exceeds the external ambient.
4. The elevated pressure raises the boiling point of water (the same mechanism
   as a pressure cooker). At 120°C, saturated steam pressure is ~2 atm.
   At 130°C, ~2.7 atm. This is achievable with ordinary kitchen-grade
   containment (mason jar, sealed metal tin, simple clamp box).
5. The hot, pressurized steam penetrates the surface of the enclosed material
   and condenses inside it — transferring both heat and moisture.

This is not identical to a pressure cooker. The key geometric difference:

**In a pressure cooker, the food is surrounded by water vapor.**
**In the P1S system, the material is surrounded by a gel that is simultaneously
releasing steam AND pressing against the material surface.**

The gel contact pressure + the steam pressure operate together. This dual-mode
pressure delivery is the structural novelty.

## 1.2 — The Isostatic Pressure Geometry

P1S is thixotropic and transparent to compression (P1S Property 8 from the v1
derivation). It distributes mechanical pressure uniformly across any surface it
contacts, regardless of surface geometry. This means:

- Irregular surfaces receive uniform contact pressure
- The gel fills pores, cracks, and anatomical features at the wood surface
- There are no pressure concentration points that would cause surface damage

As steam builds inside the closed system, the pressure acts on all surfaces
simultaneously — including the outer surface of the gel, which transmits that
pressure inward to the material surface. This is geometrically equivalent to
isostatic pressing, without the need for a hydraulic press or inert gas system.

## 1.3 — Volume Difference as a Direct Measurement of Uptake

This is structurally elegant and practically important.

At the start of the process:
- Total system volume = container volume (constant, rigid container)
- Water distributed as: liquid in gel + any air in headspace

At the end of the process:
- Water distributed as: less liquid in gel + some water absorbed into material
  + water vapor filling headspace

If the container is sealed and rigid, and the initial and final states can be
weighed (total system mass is constant), then:

**Mass of water absorbed by material = initial gel water mass − final gel water mass**

This can be measured by weighing the gel before sealing, and weighing the drained
residual gel after opening. The difference is what entered the material.

This gives a **direct, instrument-free measurement of impregnation depth** without
destroying the sample. This measurement is currently only possible in research
settings using vacuum pressure equipment and gravimetric techniques. The P1S
closed system makes it accessible on any bench with a kitchen scale.

## 1.4 — The Critical Temperature Points

Wood plasticization physics is well-established and directly applicable:

```
CRITICAL TEMPERATURES FOR WOOD CELL WALL COMPONENTS:

Component           Glass Transition (Tg, wet wood)   Effect above Tg
─────────────────────────────────────────────────────────────────────────
Hemicellulose       ~40–60°C (wet conditions)          Viscous flow; first
                                                       to plasticize
Lignin              ~80–100°C (wet conditions)         Thermoplastic flow;
                                                       dominant bending enabler
Cellulose (amorphous) 60–80°C (wet conditions)        Partial softening;
                                                       contributes to ductility
Cellulose (crystalline) >200°C                        Does not soften;
                                                       structural backbone
                                                       of wood throughout process

TARGET PROCESSING TEMPERATURE:
  Steam bending / plasticization: 100°C at minimum, 100–130°C optimal
  Achievable in P1S closed system at: modest external heat (oven, heat gun,
  hot plate) + containment long enough for the internal temperature to
  equilibrate to steam saturation point at system pressure
```

The P1S gel acts as a thermal flywheel during the initial heat ramp. It absorbs
incoming heat through evaporation, which:
- Prevents the wood surface from being subjected to a heat spike above 100°C
  before the interior has reached plasticization temperature
- Ensures the internal temperature rise is coupled to steam condensation inside
  the wood, not to conducted heat from a dry hot surface

This is the same property that makes P1S effective for fire protection, here
applied to controlled material processing.

---

---

# PART II: P1S-SOFTENER FORMULATION

## 2.1 — What the Softener Must Do

The P1S-softener is a P1S base gel with dissolved plasticizing agents added
to the water phase. When steam from the gel condenses inside the wood, the
dissolved plasticizer travels with the water into the wood cell wall.

This is the agent delivery mechanism: **steam condensation carries dissolved
agent into the material via the same pathway that water travels.**

The wood science research literature identifies several water-soluble wood
plasticizers that are safe, accessible, and validated:

```
CANDIDATE PLASTICIZING AGENTS FOR P1S-SOFTENER:

Agent             Concentration   Mechanism                     Notes
────────────────────────────────────────────────────────────────────────────────
Urea              10–20% by mass  Disrupts hydrogen bonds       Food-grade.
                  in water phase  in cellulose/hemicellulose.   Dissolves directly.
                                  Raises wood moisture above     Decomposes ~135°C
                                  equilibrium. Research-         releasing ammonia —
                                  validated plasticizer.         stay below 130°C.

Glycerol          10–15%          Bulks cell wall, acts as      Food-grade.
                                  internal plasticizer by        Cosmetic-grade.
                                  reducing Tg of amorphous       Fully miscible
                                  regions. Well documented.      with water.

Lignosulfonate    5–10%           Water-soluble lignin           Sold as agricultural
(sodium or        (lower          derivative — functions as      dust suppressant.
ammonium)         concentration   humectant and dispersant.      Available from farm
                  than urea)      Penetrates wood cell wall      supply stores.
                                  and modifies lignin–
                                  cellulose bonding.

Ammonia water     1–3% NH₃        Swells cell walls, increases   Requires ventilation.
(dilute)          by volume       accessibility of water to      Use only if ventilation
                                  cell wall microfibrils.        is reliable. Validated
                                  Used historically in wood      in the literature.
                                  bending prep.

Propylene glycol  5–10%           Similar to glycerol. Cell      Food/cosmetic grade.
                                  wall bulking, Tg depression.
```

**Recommended starting formulation for bench development:**
```
P1S-SOFTENER BASE FORMULATION:

  P1S gel (standard Phase 1 formula)     1000g
  Urea                                     80g  (dissolved in water before
                                               adding to gel, or mixed
                                               into gel water phase at A1)
  Glycerol                                 60g  (added to gel water phase)

  Total water-phase plasticizer load:     ~14% by mass of the water in the gel
  
  Method: Dissolve urea and glycerol in the 1000mL water used for SAP hydration
  (Step A1 of the standard protocol). The rest of the P1S synthesis proceeds
  identically. The plasticizers distribute uniformly through the water phase
  and will travel with steam condensate into the wood during processing.
```

## 2.2 — The Softener Process Protocol

```
SOFTENER PROCESS SEQUENCE:

Pre-condition the wood:
  Preferred: air-dried wood at 10–20% moisture content
  Avoid: kiln-dried below 8% (lignin has been set by dry heat — partially
         irreversible without extended re-hydration first)
  Avoid: green wood above 25% MC (already plasticized — may not benefit from
         additional moisture, and the steam process may over-plasticize,
         reducing structural quality after drying)

Step 1 — Apply P1S-softener to all surfaces of the wood:
  Coat the workpiece uniformly with 3–5mm of P1S-softener gel.
  Ensure full surface coverage including end grain.
  End grain is the dominant entry path for steam — priority surface.

Step 2 — Seal in container:
  Place gel-coated workpiece in a rigid container with a loose-fitting
  lid or a lid with a small pressure relief (a small hole or gap that
  allows steam to escape above a threshold pressure, preventing unsafe
  pressure buildup, while maintaining elevated pressure below that threshold).
  
  Simple bench container options:
  - Metal biscuit tin with tight lid (vents at corners under pressure)
  - Wooden box with removable lid weighted by a heavy object
  - PVC pipe sealed at both ends with caps for long thin workpieces
  - Commercially available food steamer basket inside a covered pot
  
  The container does not need to be perfectly sealed. It needs to maintain
  elevated humidity and partial steam pressure. A small gap that allows
  steam to escape slowly is safer than a fully sealed pressure vessel.

Step 3 — Apply heat:
  External heat source: oven at 110–130°C, or heat gun directed at container,
  or place container in a pot of simmering water.
  
  Duration guidelines (from standard steam bending practice, adjusted for
  gel delivery method):
    Softwood (pine, cedar, spruce): 30–45 min per 25mm thickness
    Hardwood (oak, ash, beech, maple): 45–60 min per 25mm thickness
    Thin material (<10mm): 20–30 min
    
  These are starting points. The P1S system accelerates moisture delivery
  vs. open steam box, so actual softening may occur faster. Test and refine.

Step 4 — Check for plasticization:
  Open container. Using gloved hands, try to bend a thin test piece of the
  same wood species and dimension that was processed alongside the workpiece.
  
  Adequately plasticized: bends significantly without cracking (rubber-like)
  Under-plasticized: stiff, cracks at small bend radius
  Over-plasticized: very soft, surface feels spongy — reduce time in future
  
  Work quickly. The workpiece is malleable while hot and wet. It begins to
  set as it cools below lignin Tg (~80–100°C).

Step 5 — Form the workpiece:
  See Part III.
```

---

---

# PART III: FORMING

## 3.1 — The Geometry of Forming Plasticized Wood

When lignin is above its glass transition temperature (~80°C in wet wood), wood
behaves like a thermoplastic composite. The lignin flows, the hemicellulose
softens, and the cellulose fibers (which remain crystalline and do not flow)
act as a flexible fiber reinforcement in a softened matrix.

The key mechanical constraint: wood fibers in compression can support very
large strains (up to 30% compression without failure when plasticized).
Wood fibers in tension fail at much lower strains (~1–2% even plasticized).

This means all forming operations should **compress the outer fibers and
allow tension on the inner radius to be minimized or avoided**. A back-strap
(rigid support on the tension face) is standard steam-bending practice for
this reason.

```
FORMING EQUIPMENT CLASSIFICATION (by scale):

SMALL SCALE (utensils, handles, decorative pieces, boat ribs):
  - Bent wood forms: wood or metal template, clamps
  - Compression back strap: a steel strap or thick aluminum strip held
    against the outer tension face of the bend
  - Toggle clamps for rapid securing before workpiece cools
  - Time to secure to form: must happen within ~2–4 minutes of removal
    from the steam enclosure

MEDIUM SCALE (chairs, sleds, paddles, curved structural members):
  - Laminated bending forms (plywood jig)
  - Multiple clamp stations at regular intervals
  - Two-person operation: one bends, one clamps
  - Target: all clamps engaged within 3–5 minutes

LARGE SCALE (canoe gunwales, ribs, architectural elements):
  - Pre-built bending jig with indexed clamp positions
  - P1S-softener process done in sections (one end clamped while other
    end is still in the steam enclosure)
  - OR: large container allows full workpiece to be processed then transferred
    to bending jig completely

CANOE RIBS (specific case mentioned in the derivation):
  Canoe rib bending is a traditional and well-documented use case for
  steam bending. White cedar, white oak, and white ash are the preferred
  species. Standard thicknesses: 3–6mm for ribs.
  
  At 3–6mm thickness, P1S-softener processing time is 15–25 minutes.
  The workpiece is thin enough that steam penetration from all surfaces
  achieves full through-thickness plasticization rapidly.
  
  The P1S system advantage over a standard steam box for canoe ribs:
  - Urea/glycerol plasticizers in the steam condensate may allow tighter
    bend radii than water steam alone (cell wall is chemically modified,
    not just thermally)
  - The oxygen exclusion property of the gel prevents any surface oxidation
    or char during processing
  - The gel can be applied in the field with no steam box infrastructure —
    relevant for traditional boatbuilding contexts
```

## 3.2 — Drying and Dimension Setting

After forming on the mold, the workpiece must be held until it cools below
lignin Tg AND until sufficient moisture has been lost to set the new shape.

```
CRITICAL INSIGHT FROM WOOD SCIENCE:

Wood "sets" in the bent position through two distinct mechanisms:
1. THERMAL SETTING: Cooling lignin below its Tg locks the new configuration.
   This happens within minutes after removing from heat.
2. MOISTURE SETTING: Drying the wood below ~10% MC creates permanent hydrogen
   bonds across the cell wall in the new configuration.
   This takes hours to days depending on thickness and conditions.

If the workpiece is removed from the form before moisture setting is complete,
spring-back occurs — the wood partially returns to its original shape.

P1S-mediated drying (applying dry P1S base gel, or allowing natural drying
while held in the form) can accelerate moisture removal by:
  - The gel will absorb surface moisture from the wood as it dries
    (reverse osmotic gradient — the dry gel has lower water activity than
    the saturated wood surface)
  - The bentonite clay creates a diffusion path for moisture movement
    out of the wood surface and into the gel matrix
  - This is the same physics as using a desiccant against a wet surface,
    but conformally applied and self-removing

Standard steam bending practice: hold on form for 24–72 hours.
P1S-assisted drying: likely reduces this to 12–48 hours depending on
thickness. Quantify by measuring workpiece weight loss on form over time.
```

---

---

# PART IV: P1S-HARDENER FORMULATION

## 4.1 — What the Hardener Must Do

After the softener process, forming, and drying, the wood has a new shape and
modified cell wall chemistry. The cell wall has been plasticized — the lignin
has been redistributed slightly, the hemicellulose has been hydrolyzed
partially, and urea/glycerol are present in the cell walls.

The hardener process reverses the plasticization and adds structural reinforcement
agents into the cell wall and lumen space. The target end state is a wood product
that is:
- Harder than the original wood (surface and structural hardness)
- More dimensionally stable (less swelling/shrinking with humidity changes)
- More resistant to decay and biological attack
- Potentially fire-resistant at the surface
- Set in its new formed shape permanently

Sodium silicate (water glass) is the primary hardening agent. Its validated
effects on wood are:

```
VALIDATED EFFECTS OF SODIUM SILICATE IMPREGNATION IN WOOD:

Property              Effect                                  Source
──────────────────────────────────────────────────────────────────────────
Surface hardness      ~40% increase over untreated pine       Literature validated
Dimensional stability MEE (moisture excluding efficiency)      Requires controlled
                      depends critically on drying method      drying (see 4.3)
Fire resistance       High — glass matrix at surface           Literature validated
                      prevents ignition at low heat fluxes
Decay resistance      High — alkaline pH inhibits fungal       Literature validated
                      growth; glass matrix blocks hyphal
                      penetration
Water resistance      Moderate — silicate is hygroscopic       Drawback — address
                      if not fully cured; can re-absorb        with drying process
                      moisture unless properly dehydrated
Insect resistance     High — glass matrix is inedible          Literature validated
```

## 4.2 — P1S-Hardener Formulation

```
P1S-HARDENER BASE FORMULATION:

  P1S gel (standard Phase 1 formula)          1000g base

  PRIMARY HARDENING AGENTS (add to water phase before SAP hydration):
  
  Sodium silicate solution (40% concentration) 100–150g
    (Available: hardware store as concrete sealer, art supply as
     "water glass", pottery supply. Usually sold as 40% solution.
     Dilute to ~20% in the gel water phase for good penetration.)
    
    Chemistry note: In the P1S-hardener gel, sodium silicate is dissolved
    in the water phase. During steam condensation into the wood, the silicate
    travels with the water into the wood lumen and cell wall. On subsequent
    drying and curing at 44–60°C, the silicate polymerizes within the wood
    structure, forming an insoluble silica glass network bonded to the wood
    cell wall. This is petrification chemistry, time-compressed from
    geological scale to hours.

  SECONDARY AGENTS (optional, add to water phase):
  
  PVA emulsion (wood glue base)        20–40g   Cell wall and lumen reinforcement.
                                                 Water-based polymer that deposits
                                                 on cell wall surfaces during
                                                 condensation drying.
  
  Borax (sodium borate)                10–20g   Decay and fire resistance.
                                                 Water-soluble, penetrates cell wall.
                                                 Synergistic with silicate for
                                                 fire performance (boron + silica).

  FORMULATION CONSTRAINT:
  Do NOT add the sodium silicate directly to the P1S gel at full concentration.
  Sodium silicate is strongly alkaline (pH ~12). At high concentration, it will
  disrupt the SAP polymer network (SAP is sensitive to high ionic strength and
  pH extremes). Add to the gel water phase before SAP hydration at the dilute
  concentration (20% sodium silicate in water, then add SAP to that solution).
  
  Test gel consistency after compounding: if the gel is too thin (SAP disrupted
  by salt), increase SAP by 2g and re-test. Alternatively, use the hardener
  as a lower-viscosity liquid (thinner gel) — penetration is improved with
  lower viscosity.
```

## 4.3 — The Hardener Process Protocol

```
HARDENER PROCESS SEQUENCE:

Pre-condition the formed, dried workpiece:
  Target moisture content of workpiece before hardener process: 8–12% MC
  Workpiece must be fully formed and set — no further bending after this step
  
  If the workpiece is too wet (>15% MC): allow to air-dry further before
  hardener process. A wet cell wall has less capacity to absorb additional
  water-borne agents (cells already filled with water).

Step 1 — Apply P1S-hardener to all surfaces:
  Same geometry as softener process: 3–5mm gel coating on all surfaces.
  Ensure end grain coverage — end grain penetration is the primary pathway.

Step 2 — Seal in container and apply heat:
  Same procedure as softener process.
  Temperature target: 100–120°C (lower than softener process is acceptable —
  the goal is agent delivery, not plasticization; lower temperature is fine).
  
  Duration: 30–60 minutes depending on thickness.
  
  CAUTION: At temperatures above 130°C, sodium silicate begins to gel prematurely
  and may block penetration pathways before reaching depth. Stay below 130°C.

Step 3 — Controlled drying and cure:
  This step is critical for sodium silicate performance.
  
  Immediate post-process drying at low temperature (44–60°C, oven or in sun):
    - Drives water out of cell lumens
    - Triggers sodium silicate polymerization in place within the wood structure
    - Produces an insoluble silica glass network bonded to wood cells
  
  Avoid:
    - Rapid drying at high temperature (>90°C) immediately post-process:
      causes internal steam pressure that can check (crack) the wood
    - Drying in free air without controlled temperature:
      sodium silicate can absorb ambient moisture and fail to cure properly
  
  Recommended drying schedule:
    Phase 1: 44–60°C in oven for 2–4 hours (triggers silicate polymerization)
    Phase 2: Air dry at room temperature to equilibrium MC (24–48 hours)
    Phase 3 (optional): Final cure at 80–100°C for 1 hour
                        (maximizes silicate cross-linking and water resistance)

Step 4 — Surface finishing:
  After curing, the surface has a slight glassy sheen from silicate at the surface.
  This is aesthetically usable or can be:
    - Lightly sanded (240+ grit) for a matte surface
    - Oiled (linseed, tung) — silicate-treated wood accepts oil well,
      and oil further improves water resistance
    - Left as-is for maximum fire and water resistance at surface
```

---

---

# PART V: PROCESS CHAIN ANALYSIS AND FAILURE MODES

## 5.1 — Where the Process Works Well

```
OPTIMAL CANDIDATE MATERIALS:

Species               Suitability   Notes
────────────────────────────────────────────────────────────────────────
White oak             EXCELLENT     Classic steam bending species.
                                    Dense; benefits from urea plasticizer.
White ash             EXCELLENT     Flexible even dry; P1S accelerates.
White cedar           VERY GOOD     Traditional canoe rib material.
                                    Low density; fast penetration.
Hickory               VERY GOOD     High strength after hardening.
Beech                 VERY GOOD     European tradition; good results.
Cherry                GOOD          Moderate bending radius.
Walnut                MODERATE      Heartwood resists penetration.
Pine (air-dried)      MODERATE      Softwood; plasticizes easily but
                                    spring-back is higher than hardwoods.
Bamboo                POTENTIAL     High silica content already; the
                                    hardener process may be less effective
                                    (competing silica); softener process
                                    highly relevant.
```

## 5.2 — Failure Modes and Structural Mitigations

```
FAILURE MODE 1: SPRING-BACK (wood returns toward original shape after forming)

  Cause: Insufficient moisture setting. Wood removed from form while still
         wet and above ~10% MC.
  
  Mitigation:
    1. Hold on form longer (minimum 24 hours; 48 for thick sections)
    2. Increase urea concentration in softener (urea itself contributes to
       permanent set by modifying hydrogen bonding in cell wall on drying)
    3. Apply hardener process before removing from form — harden in place
       on the forming jig (the silicate network resists spring-back mechanically)

FAILURE MODE 2: SURFACE CHECKING (surface cracks during drying)

  Cause: Outer surface dries faster than interior, creating a tensile stress
         in the surface layer that exceeds modulus of rupture.
  
  Mitigation:
    1. Slow the drying rate by re-applying a thin layer of plain P1S gel
       to the outer surface during early drying (gel controls surface
       drying rate via the same evaporative buffer mechanism used in
       fire protection)
    2. Reduce temperature during drying phase
    3. Cover workpiece loosely with plastic wrap for first 4–6 hours
       of drying (reduces rate of surface moisture loss)

FAILURE MODE 3: INADEQUATE PENETRATION (agent does not reach interior)

  Cause: Cell lumens blocked by extractives or prior kiln treatment.
         Also: processing time too short.
  
  Mitigation:
    1. Pre-soak the workpiece in plain warm water for 12–24 hours before
       P1S process — opens lumens and saturates cell walls before agent delivery
    2. Increase processing time
    3. Reduce agent concentration (lower viscosity = better penetration)
    4. Use end-grain as the primary delivery surface — seal side grain with
       plain P1S and leave end grain open for maximum penetration

FAILURE MODE 4: SAP DISRUPTION BY SODIUM SILICATE (gel too thin in hardener formulation)

  Cause: Sodium ions from silicate compete with the SAP's osmotic uptake,
         reducing gel viscosity.
  
  Mitigation:
    1. Increase SAP by 2–3g in hardener formulation
    2. Reduce silicate concentration (20% vs 40% solution)
    3. Accept lower viscosity hardener — a thinner gel penetrates better anyway
    4. Use the hardener as a brush-applied liquid rather than a gel (less SAP,
       more water, same dissolved silicate)
```

## 5.3 — Volume Difference Measurement Protocol

```
QUANTIFYING IMPREGNATION:

Before process:
  1. Weigh workpiece (W₁)
  2. Weigh total gel to be applied (G₁)

After process, before drying:
  3. Weigh workpiece (W₂)
  4. Weigh recovered gel (remove from surface) (G₂)

Calculation:
  Agent uptake into wood = W₂ − W₁
  Gel water lost to steam = G₁ − G₂ − (agent uptake)
  
  Cross-check: G₁ − G₂ should equal agent uptake + steam lost to headspace
  
  If agent uptake / wood volume (cm³) > 0.05 g/cm³: meaningful penetration
  If agent uptake / wood volume (cm³) < 0.02 g/cm³: surface treatment only

This measurement directly calibrates process parameters (temperature, time,
concentration, wood species) and builds a dataset for process optimization
without any specialized equipment.
```

---

---

# PART VI: SCALABILITY AND PRODUCT PATHWAY

## 6.1 — Scale Spectrum

```
BENCH SCALE (current, this protocol):

  Container size: 1–5 liters
  Workpiece size: thin strips, small blocks, utensil blanks
  Equipment: kitchen oven, metal container, clamps
  Cost: under $5 per workpiece in materials
  
  Validate here: species response, agent concentration, timing.
  Produce: prototype product (spoon, small paddle blade, bent bracket, 
           decorative piece)

WORKSHOP SCALE:

  Container: custom-built insulated wooden box with metal interior lining,
             external heat source (propane burner, heat blanket)
  Workpiece size: chair components, canoe ribs, gunwales up to 3m
  Equipment: box, clamps, bending jig, oven for drying phase
  Cost: under $200 to build the box; under $20 per workpiece in materials

INDUSTRIAL SCALE:

  Container: pressurized impregnation cylinder (standard in pressure treatment
             industry — already exists commercially)
  P1S gel replaces or supplements conventional impregnation liquid
  Key industrial advantage: P1S gel's thixotropic behavior means it does
  NOT drain out of vertical or complex-geometry surfaces during processing —
  conventional liquid impregnants drain away from vertical surfaces
  
  Products: architectural bent timber, performance sporting goods
            (hockey sticks, archery bows, kayak paddles, skis),
            custom furniture, marine components
```

## 6.2 — The Canoe Derivation (Explicit)

```
CANOE RIB FORMING — DERIVATION FROM FIRST PRINCIPLES:

Target product: bent wood canoe rib, white cedar or white oak, 5mm × 20mm × ~1.2m

Standard process (without P1S):
  - Steam box: 30–45 minutes at 100°C atmospheric steam
  - Transfer to bending jig: 2–3 minutes window before cooling
  - Clamp and hold: 24–72 hours
  - Yield: ~70–80% (brittle failures are common; extractives in dry wood
    resist plasticization)
  - No chemical modification — same material properties as untreated wood
    after processing

P1S-softener process:
  - Apply P1S-softener (urea + glycerol in gel), seal in pipe section, 20–25 minutes
  - Transfer to bending jig: same 2–3 minute window
  - Clamp and hold: 24 hours (gel-assisted drying may reduce this)
  - Expected yield improvement: urea-assisted bending is documented to
    reduce bend failures and allow tighter radii than steam alone
  - Cell wall modification: urea present in cell wall on drying creates
    a slightly modified lignin–cellulose interaction that reduces spring-back

P1S-hardener process (optional second stage for canoe ribs):
  - Apply P1S-hardener (silicate + borax) to bent and dried ribs
  - 30 minutes in sealed enclosure at 110°C
  - Low-temperature cure at 50°C for 2 hours
  - Result: silicate + boron throughout the rib — increased hardness,
    fire resistance, decay resistance, insect resistance
  - The treated rib requires no additional wood preservative for exterior use

Final canoe product comparison to conventional:
  CONVENTIONAL CANOE RIB: bent cedar, no treatment, requires paint or varnish
                           for weather protection, susceptible to rot at
                           gunwale fastening points over time
  
  P1S-PROCESSED RIB: bent, hardened with silicate, borax-protected,
                     naturally fire-resistant, rot-resistant, dimensionally
                     stable — no paint required for structural protection
                     (finish can be applied for aesthetics only)
```

## 6.3 — Process Chain Summary

```
PROCESS CHAIN DECISION TREE:

Is the workpiece requiring SHAPING (bending, forming)?
  YES → P1S-softener first → Form → Dry → optionally P1S-hardener
  NO  → P1S-hardener only (impregnation without forming)

Does the product require DURABILITY (outdoor, marine, structural)?
  YES → P1S-hardener (sodium silicate + borax)
  NO  → P1S-softener forming only; standard drying

Does the product require FIRE RESISTANCE?
  YES → P1S-hardener (sodium silicate, optimized for fire)
        OR standard P1S as surface coating (wet, temporary protection)
  NO  → Process without fire additive

Does the product require MAXIMUM BEND RADIUS (tight curves)?
  YES → P1S-softener with urea + glycerol at upper concentration range
        + back-strap compression on tension face
        + maximum processing time before forming
  NO  → Standard softener formulation, standard time
```

---

---

# PART VII: WHAT IS NOVEL WITH RESPECT TO PRIOR ART

## 7.1 — What Exists

Steam bending: ancient practice, well documented, commercially practiced.
Urea wood plasticization: researched in academic literature.
Glycerol wood plasticization: researched.
Sodium silicate wood impregnation: researched and practiced commercially.
Pressure impregnation of wood: large industrial sector (pressure treatment).

## 7.2 — What Does Not Exist

```
NOVEL ELEMENTS IN THIS DERIVATION:

1. P1S AS SIMULTANEOUS STEAM GENERATOR + AGENT CARRIER + ISOSTATIC PRESSURE MEDIUM
   
   No prior system combines: (a) conformal encapsulation of the workpiece,
   (b) in-situ steam generation from the encapsulant, (c) agent delivery via
   steam condensation, and (d) isostatic pressure from gel contact — in a single
   material applied at room temperature from commodity components.
   
   The pressure treatment industry uses tanks, pumps, valves, and liquid
   impregnants that drain from complex surfaces. The P1S system uses a gel
   that clings to all surfaces, generates its own steam, and delivers agents
   conformally to every surface simultaneously.

2. DUAL-PHASE PROCESS (SOFTEN THEN HARDEN) AS INTEGRATED SYSTEM
   
   Steam bending and silicate hardening have never been described as a
   two-stage P1S-mediated integrated process. Each has been studied
   independently. Combining them via the same delivery mechanism (P1S closed
   system) with the volume-difference measurement protocol creates a
   verifiable, calibratable, bench-accessible process chain.

3. VOLUME DIFFERENCE AS BENCH-ACCESSIBLE IMPREGNATION MEASUREMENT
   
   Determining impregnation depth in wood normally requires destructive
   cross-section analysis or specialized equipment (NMR, X-ray). The
   gravimetric mass-balance approach enabled by the P1S closed system
   makes this measurement available without any specialized equipment.

4. P1S-MEDIATED DRYING AS MOISTURE CALIBRATION STEP
   
   Using the gel's osmotic differential to actively pull moisture from
   the workpiece surface during the drying phase (by applying dry P1S
   base gel to the formed workpiece on the mold) is not described in
   the steam bending or wood processing literature.
```

---

*End of EHAB P1S — Closed-System Steam-Forming and Hardening Process (v3)*
*This document extends the v1 application class 1.1 (wood modification under encapsulation)*
*and v2 application class 10.1 (anoxic hot-press processing medium) into a full process chain.*
*Agent-specific formulations, process protocols, failure modes, and scalability analysis*
*are new in v3 and not present in prior documents.*
