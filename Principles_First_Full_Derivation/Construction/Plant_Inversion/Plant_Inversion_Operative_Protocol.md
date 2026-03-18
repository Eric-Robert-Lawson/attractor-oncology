# PLANT ARCHITECTURE INVERSION
# OPERATIVE PROTOCOL v2.3
## Two-Chamber Cone Membrane System
## Complete Step-by-Step Execution Document
## OrganismCore — Eric Robert Lawson
## 2026-03-18

---

## STATUS: ACTIVE — OPERATIVE PROTOCOL v2.3
## Classification: Full execution document.
## Supersedes: Plant_Inversion_Operative_Protocol_v2.2.md
## Amendment notes v2.3:
## (1) Scenario 1 corrected to true normal control.
## S1 now contains a UV-A LED free in air space with
## zero penetration cost. This is the only scenario
## where UV is freely available without a foam barrier.
## Previous v2.2 Scenario 1 had no UV anywhere — this
## was not a valid control.
## (2) Substrate types revised and expanded to three
## distinct types plus one free UV condition:
## Type A — hydroponic foam, no UV
## Type B — UV-only foam, no hydroponic nutrients
## Type C — combined foam, hydroponic + UV co-located
## Free UV — no foam, UV-A LED free in S1 air space
## (3) Scenarios 2 and 3 now use Type B foam —
## UV signal separated from hydroponic signal into
## different chambers. Seedling faces competing
## attractors in opposite directions.
## (4) Scenarios 4 and 5 now use Type C foam —
## UV and hydroponic signals co-located in the same
## chamber. No competing signal in the opposite
## direction.
## (5) Cone orientation table updated for all 5
## scenarios under the revised substrate definitions.
## Timestamp: 2026-03-18

---

## LINKED RECORDS

- Pre-registration (v1.0): Plant_Inversion_Pre_Registry.tex
- Geometric framework: Plant_Inversion.md
- Triadic invariant:
  TRIADIC_INVARIANT_BIOLOGY_REASONING_ARTIFACT.md
- Confirmed cases:
  ATTRACTOR_GEOMETRY_CONFIRMED_CASES_REASONING_ARTIFACT.md
- One-pager: One_Pager_Plant_Inversion.tex
- Cost and failure budget:
  plant_inversion_cost_and_failure_budget.md
- Supersedes:
  Plant_Inversion_Operative_Protocol_v2.2.md
- Supersedes:
  Plant_Inversion_Operative_Protocol_v2.1.md
- Supersedes:
  Plant_Inversion_Operative_Protocol_v2.0.md
- Supersedes:
  Plant_Inversion_Operative_Protocol.md (v1.0)

---

## PREAMBLE — WHAT THIS DOCUMENT IS

This is the document you have open on the bench while
running the experiment. It is not a framework document.
It is not a reasoning artifact. It is a step-by-step
procedure with failure modes, recovery actions, and
success confirmation criteria at every step.

Read Plant_Inversion.md and
TRIADIC_INVARIANT_BIOLOGY_REASONING_ARTIFACT.md first
for the geometric framework and the predictions this
experiment is testing. This document assumes you have
read those documents.

The core prediction being tested:

> When the coherence gradient field experienced by a
> germinating seed is inverted before germination
> completes — the substrate carrying UV-A radiation
> is placed below the seed rather than above — the
> resulting organism will exhibit architectural
> inversion: root growing upward toward the UV
> substrate, shoot growing downward to penetrate it
> and access the embedded UV-A signal.
> Same genome. Different field. Different organism.

---

## THE CENTRAL DESIGN PRINCIPLE —
## SIGNAL DEPRIVATION AND SKOTOMORPHOGENESIS

This principle governs the entire experimental design
and must be understood before any component of the
protocol is executed.

**The experiment does not depend on UV-A being a
maximally efficient phototropic signal. It depends
on UV-A being the only photonic signal the seedling
can detect in all scenarios except Scenario 1.**

Phototropin receptors (phot1 and phot2) in plant
seedlings respond to the UV-A and blue light range
(390–500nm). UV-A alone activates these receptors
less efficiently than peak blue light at 450nm.
Under normal conditions with ambient light present,
a UV-A-only source would be a weak directional signal
competing against background photonic noise.

This experiment does not use normal conditions for
Scenarios 2 through 5.

When a germinating seedling develops in complete
darkness, it enters skotomorphogenesis — the
developmental program for light-seeking under total
light deprivation. In this state:

- Hypocotyl elongation is maximal — the seedling is
  physically reaching for any signal
- All photoreceptors are primed at maximum sensitivity
  — detection threshold is at its lowest possible value
- There is zero ambient photonic noise — no competing
  signal from any direction

Under these conditions the absolute efficiency of
UV-A as a phototropic driver is irrelevant. The
embedded UV-A source is not competing against
anything. It is the only signal in the coherence
field. A seedling in a completely dark chamber with
a single UV-A source embedded in a substrate is not
asking whether UV-A is optimal — it is responding
to the only photonic input that exists in its entire
environment.

This is the difference between a shout in a noisy
room and a book falling in a silent room. The book
falling is orders of magnitude more significant in
the silent room than the shout is in the noisy room
— not because the book is loud but because nothing
else is making any sound.

**Scenario 1 is the explicit exception to this
principle.** In Scenario 1, UV is freely available
in the S1 air space with no foam barrier. This is
the true normal control — it replicates the
conditions under which a seedling would germinate
in nature with light above and nutrients below.
Scenario 1 confirms the biological system and
apparatus are functioning correctly before any
inversion results are interpreted.

**The engineering consequence:**
Complete external light deprivation of all chambers
during the growth phase is the load-bearing condition
of Scenarios 2 through 5. If ambient light enters
these chambers from outside, it contaminates the
coherence field and undermines the signal isolation
the experiment depends on.

Complete darkness outside the embedded UV-A source
is what makes Scenarios 2 through 5 work.

---

## PART I — EXPERIMENTAL DESIGN OVERVIEW

### The Substrate States

All experimental variation is built from three
substrate types and one free UV condition.

---

**FREE UV CONDITION — S1 of Scenario 1 only:**
- No foam substrate in S1 chamber
- UV-A LED mounted free in the S1 chamber air space,
  facing downward toward the seed
- UV-A radiation freely available throughout the S1
  air space
- Zero penetration cost for the shoot to access UV-A
- No water, no minerals in S1
- Used in: S1 of Scenario 1 only

---

**SUBSTRATE TYPE A — Hydroponic foam, no UV:**
- Material: reticulated polyurethane foam 20ppi,
  saturated with hydroponic mineral solution
- Signals present: water, minerals
- UV present: NO — at any depth
- Penetration requirement: seedling must grow into
  foam to access water and minerals. No UV signal
  is present within it at any depth.
- Used in: S2 of Scenario 1, S2 of Scenario 2,
  S1 of Scenario 3

---

**SUBSTRATE TYPE B — UV-only foam, no hydroponic:**
- Material: reticulated polyurethane foam 20ppi,
  dry or minimally moist (not saturated with
  nutrient solution), with UV-A LED strip embedded
  within the foam body at 10–20mm depth from the
  seed-facing surface
- Signals present: UV-A radiation only
- Water and minerals present: NO
- UV present: YES — embedded inside foam body.
  Seedling must physically penetrate foam to reach
  UV-A signal.
- Penetration requirement: seedling must grow into
  and through the foam to reach the UV-A signal
- Used in: S1 of Scenario 2, S2 of Scenario 3

---

**SUBSTRATE TYPE C — Combined foam, hydroponic
+ UV co-located:**
- Material: reticulated polyurethane foam 20ppi,
  saturated with hydroponic mineral solution, with
  UV-A LED strip embedded within the foam body at
  10–20mm depth from the seed-facing surface
- Signals present: water, minerals, UV-A radiation
  — all co-located in the same substrate
- UV present: YES — embedded inside foam body.
  Seedling must physically penetrate foam to reach
  UV-A signal along with water and minerals.
- Penetration requirement: high — seedling must
  grow into and through the foam to access all
  signals simultaneously
- Used in: S1 of Scenario 4, S2 of Scenario 5

---

### The 5 Scenarios

| Scenario | S1 (Top Chamber) | S2 (Bottom Chamber) | Cone Orientation |
|---|---|---|---|
| **1** | FREE UV — UV-A LED free in air, no foam | TYPE A — hydroponic foam, water + minerals, no UV | Apex up |
| **2** | TYPE B — UV-only foam, UV embedded, no hydroponic | TYPE A — hydroponic foam, water + minerals, no UV | Apex up |
| **3** | TYPE A — hydroponic foam, water + minerals, no UV | TYPE B — UV-only foam, UV embedded, no hydroponic | Apex down |
| **4** | TYPE C — hydroponic foam, water + minerals + UV embedded | EMPTY — void | Apex up |
| **5** | EMPTY — void | TYPE C — hydroponic foam, water + minerals + UV embedded | Apex down |

---

### The Governing Cone Orientation Rule

> **The apex always points toward the void or least
> active substrate. The wide opening always faces
> the most active substrate.**

Empty void is always least active.
Free UV in air (Scenario 1 S1) is active but presents
no liquid threat — apex up is correct.
Type A is less active than Type B or Type C.
Type B and Type C are the most active substrates.

| Scenario | Most active location | Cone orientation | Apex points toward |
|---|---|---|---|
| 1 | Free UV in S1 above, Type A in S2 below | Apex up | Away from Type A liquid below |
| 2 | Type B UV foam in S1 above | Apex up | Away from Type B above |
| 3 | Type B UV foam in S2 below | Apex down | Away from Type B below |
| 4 | Type C combined foam in S1 above | Apex up | Toward empty S2 below |
| 5 | Type C combined foam in S2 below | Apex down | Toward empty S1 above |

---

### What Each Scenario Tests

**Scenario 1 — True normal control:**
UV-A is freely available in the S1 air space above
the seed. Zero penetration cost. The shoot emerges
into UV-lit air. The root must penetrate Type A foam
below for water and minerals. This replicates normal
natural conditions — light above, nutrients below,
no barrier to either signal.

This scenario confirms the biological system and
apparatus function correctly. It must record NORMAL
architecture in every experimental run. If it does
not, the run is invalid regardless of what other
scenarios show.

Expected outcome: root down into S2 Type A foam,
shoot up into free UV-lit air. Normal architecture.
No penetration required by shoot.

---

**Scenario 2 — UV signal above, separated from
hydroponic, penetration required:**
The UV-A signal is embedded in Type B foam above.
The shoot must penetrate the Type B foam to access
UV-A. There are no hydroponic nutrients in S1 — the
water and mineral signal is in S2 below. The seedling
faces two separated attractors: UV above requiring
upward penetration, hydroponic nutrients below
requiring downward penetration.

Expected outcome: shoot grows upward penetrating
Type B foam to reach UV-A. Root grows downward into
Type A foam to access water and minerals. Normal
architecture but with active shoot penetration of
S1 Type B required to reach UV.

---

**Scenario 3 — Full inversion of Scenario 2.
PRIMARY EXPERIMENTAL CLAIM:**
Everything is structurally identical to Scenario 2
except the positions of Type B and Type A are
swapped and the cone membrane is flipped. The UV-A
signal is now embedded in Type B foam below. The
hydroponic nutrients are now in Type A foam above.
The seedling faces the same two separated attractors
but in opposite directions: UV below requiring
downward penetration, hydroponic nutrients above
requiring upward penetration.

The prediction: the UV-A signal, as the only photonic
input in the entire chamber system, drives
architectural inversion. The root grows upward
toward nutrients in Type A foam above. The shoot
grows downward penetrating Type B foam below to
reach the only UV-A source available.

Expected outcome: root up penetrating S1 Type A
foam toward water and minerals, shoot down
penetrating S2 Type B foam to reach embedded UV-A.
Full architectural inversion.

---

**Scenario 4 — All signals co-located above,
nothing below:**
All signals — water, minerals, and UV-A — are
co-located in Type C foam in S1 above. S2 is void.
The seedling has every signal concentrated above and
nothing below. No competing attractor in any other
direction.

Expected outcome: shoot grows upward penetrating
S1 Type C foam to access all signals simultaneously.
Root grows downward into void or is indeterminate.
Normal architecture with maximum signal concentration
above and zero competition.

---

**Scenario 5 — All signals co-located below,
nothing above:**
All signals — water, minerals, and UV-A — are
co-located in Type C foam in S2 below. S1 is void.
Every signal is concentrated below. No competing
attractor above.

Expected outcome: shoot grows downward penetrating
S2 Type C foam to access all signals simultaneously.
Root grows upward into void or is indeterminate.
Full architectural inversion with maximum signal
concentration below and zero competition.

---

### The Critical Comparison Structure

**Primary inversion test:**
Scenario 2 vs Scenario 3. Structurally identical —
same substrate types, same signal strengths, same
penetration requirements. The only difference is
which direction holds Type B UV foam. Cone is
flipped. Prediction: architecture flips with it.

**Penetration cost for UV access:**
Scenario 1 (UV free in air, zero cost) vs
Scenario 2 (UV in Type B foam above, high cost).
Same direction, same hydroponic below, only
penetration cost differs.

**Signal separation vs co-location:**
Scenarios 2 and 3 (UV and hydroponic separated,
competing attractors) vs Scenarios 4 and 5 (UV and
hydroponic co-located, single coherent attractor).

**Maximum inversion signal:**
Scenario 3 (inversion with competing hydroponic
above) vs Scenario 5 (inversion with nothing above,
maximum downward signal, no competition).

---

### Prediction Summary Table

| Scenario | Root Direction | Shoot Direction | Shoot Penetrates | Root Penetrates | Architecture |
|---|---|---|---|---|---|
| **1** | Down into S2 Type A | Up into free UV air | No penetration — free air | Yes — Type A foam | NORMAL — free UV access |
| **2** | Down into S2 Type A | Up, penetrates S1 Type B to reach UV | Yes — Type B UV foam | Yes — Type A foam | NORMAL + ACTIVE PENETRATION |
| **3** | Up, penetrates S1 Type A toward nutrients | Down, penetrates S2 Type B to reach UV | Yes — Type B UV foam | Yes — Type A foam | **INVERTED + ACTIVE PENETRATION** |
| **4** | Down into void or indeterminate | Up, penetrates S1 Type C to reach all signals | Yes — Type C foam | Nothing to penetrate | NORMAL — MAXIMUM SIGNAL ABOVE |
| **5** | Up into void or indeterminate | Down, penetrates S2 Type C to reach all signals | Yes — Type C foam | Nothing to penetrate | INVERTED — MAXIMUM SIGNAL BELOW |

---

## PART II — MATERIALS

### 2.1 Biological Materials

| Item | Specification | Source | Cost |
|---|---|---|---|
| Radish seeds | *Raphanus sativus*, Cherry Belle variety | Garden centre / Amazon | $2–4 |
| Quantity | Minimum 25 seeds (5 per scenario plus spares) | — | — |

**Seed sourcing note:** Use seeds from a single batch
and single supplier for all 5 scenarios. Do not mix
seed batches. Record batch number and supplier name
in the experimental log before beginning.

### 2.2 Chamber Materials

| Item | Specification | Source | Cost |
|---|---|---|---|
| PVC pipe | 50mm internal diameter, 150mm length per section, opaque black preferred | Hardware / plumbing supply | $2–4 per section |
| Quantity | 10 sections (2 per scenario × 5 scenarios) | — | $10–20 total |
| End caps | Fitted to pipe diameter, removable | Hardware store | $1–2 each, 10 needed |
| Aquarium silicone sealant | Clear, aquarium-safe, non-toxic | Hardware / aquarium store | $5–8 |
| Black opaque tape or black plastic wrap | External light exclusion | Hardware store | $2–5 |

### 2.3 Cone Membrane Materials

| Item | Specification | Source | Cost |
|---|---|---|---|
| Metallized Mylar film | Emergency blanket, aluminized both sides | Camping supply / dollar store | $1–3 per blanket |
| Small plastic funnels | Rim diameter matching pipe inner diameter (50mm) | Kitchen / laboratory supply | $1–3 each, 5 needed |
| Hydroponic net cups | 16mm or 25mm diameter | Hydroponics supplier / Amazon | $3–8 for pack of 50 |
| Aquarium silicone sealant | As above | As above | Same tube |

**Membrane barrier functions:**
The cone membrane requires two barrier functions:
1. Liquid barrier — prevent water and nutrient
   solution migrating between S1 and S2
2. UV and light barrier — prevent UV-A from the
   substrate in one chamber crossing into the other
   chamber, and prevent external light entering
   through the membrane joint

### 2.4 Substrate and Signal Materials

| Item | Specification | Use | Source | Cost |
|---|---|---|---|---|
| Reticulated polyurethane foam | 20ppi open cell, sheet form | Type A, Type B, and Type C substrate carrier | Aquarium supply / Amazon | $5–10 per sheet |
| Hydroponic mineral solution | Balanced hydroponic nutrient, seedling dilution 1/4 to 1/2 strength | Nutrient signal in Type A and Type C | Garden centre / hydroponics supplier | $5–12 |
| Distilled water | 1L minimum | Nutrient solution base | Grocery store | $1–2 |
| UV-A LED strip | 315–400nm — confirm wavelength before purchase | UV signal in Free UV condition, Type B, and Type C | Amazon | $8–15 |

**UV-A wavelength critical note:**
Confirm the LED strip wavelength before purchase.
- UV-A (315–400nm): correct
- UV-B (280–315nm): do not use — DNA damage
- UV-C (100–280nm): do not use — germicidal
- Visible blue light (400–500nm): acceptable
  substitution — more efficient phototropin
  activator than UV-A, must be documented as
  protocol deviation if used

**Note on Type B foam moisture:**
Type B foam carries UV-A only — no hydroponic
nutrients. The foam should be dry or at most
minimally moist (dampened with plain distilled
water only, no nutrients). This ensures the UV-A
signal in S1 of Scenario 2 and S2 of Scenario 3
is cleanly separated from the hydroponic signal
in the opposite chamber.

### 2.5 Free UV Condition Materials
### (Scenario 1 S1 only)

| Item | Specification | Source | Cost |
|---|---|---|---|
| UV-A LED strip | Same as above | Same as above | Same unit — cut a short length |
| Small mounting bracket or adhesive clip | To mount LED strip to inner face of S1 end cap facing downward | Hardware store | $0–2 |

**Free UV mounting note:**
In Scenario 1, the UV-A LED is mounted to the inner
face of the S1 end cap, oriented to face downward
into the S1 chamber air space toward the seed. The
LED illuminates the full air space of S1 freely.
The shoot emerges into this UV-lit space without
any foam barrier. This is the only scenario where
UV is not embedded in foam.

### 2.6 Ventilation Port Materials

| Item | Specification | Source | Cost |
|---|---|---|---|
| Reticulated foam plugs | Small pieces of 20ppi foam cut to fit ventilation port diameter | Cut from substrate foam sheet — no additional cost | $0 |
| Small USB fan (optional) | 5V USB powered, 40–50mm diameter | Amazon / electronics supply | $3–8 |
| Drill bit | Sized to ventilation port diameter (10–15mm) | Hardware store | $2–5 |

**Ventilation port design:**
Each chamber requires at least one ventilation port
in the end cap for ambient CO₂ and O₂ exchange.
The port must allow gas movement while preventing
light from entering.

The light trap is a plug of open-cell reticulated
foam pressed into the ventilation port hole. Gas
and air diffuse freely through the open-cell
structure. Light photons cannot travel in a
straight line through the tortuous foam cell path
— they are scattered and absorbed before reaching
the interior.

Active ventilation via small USB fan is recommended
for all chambers containing Type B or Type C
substrate where the LED generates heat.

### 2.7 Tools and Equipment

| Item | Use | Source | Cost |
|---|---|---|---|
| Utility knife or scalpel | Cutting foam and Mylar | Hardware / craft store | $3–8 |
| Ruler and permanent marker | Measuring and marking | — | $0–2 |
| Lighter or candle | Heating pin for port formation | — | $0–1 |
| Sewing pin or needle | Forming seed port in cone apex | — | $0–1 |
| Tweezers | Handling germinated seeds | — | $2–5 |
| Damp paper towel | Seed germination substrate | Grocery store | $0–1 |
| Ziplock bags (small) | Germination chambers | Grocery store | $1–2 |
| Camera or smartphone | Photographic record | — | $0 |
| Nitrile gloves | Handling seedlings | Pharmacy / hardware | $3–5 |
| Millimeter ruler | Measuring penetration depth | — | $0–2 |
| Red-light torch or headlamp | Observation without disrupting dark condition | Photography / astronomy supply | $5–15 |
| Drill | Ventilation ports in end caps | Hardware store | $0 if borrowed |

**Red-light torch note:**
Red light (>650nm) does not activate phototropin
receptors and does not disrupt skotomorphogenic
state. Use red-light only for any interior
inspection during the growth phase. Never use
white light, blue light, or UV near open chambers
at any point during the growth phase.

### 2.8 Complete Cost Summary

| Category | Estimated Cost |
|---|---|
| Biological materials | $2–4 |
| Chamber materials | $15–30 |
| Cone membrane materials | $8–15 |
| Substrate and signal materials | $18–37 |
| Ventilation materials | $3–13 |
| Tools and equipment | $15–35 |
| **Total estimate** | **$61–134** |

---

## PART III — LIGHT EXCLUSION REQUIREMENTS

Light exclusion is the load-bearing engineering
requirement of Scenarios 2 through 5. All design
decisions are secondary to it.

### 3.1 External Light Exclusion — Scenarios 2, 3, 4, 5

All four chamber units for Scenarios 2 through 5
must be in complete external light deprivation
during the full growth phase from T=0 to T=168.

Methods in order of preference:

**Method A — Opaque pipe material:**
Use black opaque PVC pipe. The pipe material itself
blocks external light. Preferred because it is
passive and cannot be accidentally removed.

**Method B — External opaque wrap:**
If transparent pipe is used, wrap the entire
exterior of each chamber in black opaque tape,
black plastic wrap, or black fabric before T=0.
Must cover full length of both tube sections and
all joints with no gaps.

**Method C — Dark enclosure:**
Place all chambers inside a light-tight box,
cabinet, or dark room for the full growth phase.
Must be confirmed light-tight before T=0.

Any combination of methods is acceptable. The
requirement is the result — zero external light
reaching the chamber interior.

### 3.2 Scenario 1 Light Conditions

Scenario 1 is the exception. The S1 chamber of
Scenario 1 contains a freely available UV-A LED
that illuminates the chamber interior. External
light exclusion of the S1 section of Scenario 1
is not required and would defeat the purpose of
the control.

However, the S2 section of Scenario 1 (containing
Type A hydroponic foam) should be opaque to prevent
external light reaching the hydroponic substrate
and potentially creating an unintended signal in S2.

The cone membrane between S1 and S2 of Scenario 1
must block UV-A from the free LED in S1 from
passing into S2. The Mylar membrane handles this.

### 3.3 Ventilation Port Light Exclusion

Every ventilation port in every chamber must be
fitted with a foam plug light trap before T=0.

Confirm each foam plug is seated fully with no gap
around the circumference. Apply aquarium silicone
around the plug circumference if needed.

Test each port: shine a torch at the exterior foam
plug surface in darkness and confirm no light
visible from the interior side.

### 3.4 Membrane Joint Light Exclusion

The cone membrane joint must be fully light-sealed
by the aquarium silicone bead applied during
membrane assembly. Test after curing by placing
a torch inside the assembled chamber and checking
for light leaks at the joint from outside.

### 3.5 Observation Protocol for Light Exclusion

At each observation time point:
1. Prepare camera before opening or unwrapping
   any chamber.
2. Open or unwrap one chamber at a time. Replace
   wrap before opening the next.
3. Maximum 2 minutes of light exposure per chamber
   per observation time point.
4. Use red-light torch for any interior inspection.
5. Record duration of light exposure for each
   chamber at each time point in the experimental
   log.

---

## PART IV — SUBSTRATE CONSTRUCTION

### Type B Substrate — UV-Only Foam

Type B carries UV-A only. No hydroponic nutrients.
Construct all Type B units in the same session.

**Step B1: Cut foam disc**
Cut a disc of 20ppi reticulated foam to fit snugly
in the pipe inner diameter (50mm). Thickness: 30mm.

**Step B2: Mark seed-facing surface**
Mark the face that will be oriented toward the seed
with permanent marker on the outer edge. Write
scenario number and "SEED SIDE."

**Step B3: Cut LED channel**
Cut a channel across the full diameter of the foam
disc at 10–20mm depth from the seed-facing surface.
The exact depth within this range is not critical.
The LED must be interior to the foam body, not at
the surface. Record the depth chosen for each unit.

Channel width and depth: just sufficient to seat
the UV-A LED strip flush within the foam body
without protruding above the channel surface.

**Step B4: Seat LED strip in channel**
Press the UV-A LED strip into the channel flush
with the foam body. Route power cable out through
a small notch cut in the outer edge of the foam
disc at channel depth level.

**Step B5: Do not saturate with nutrient solution**
Type B foam is dry or minimally moist (plain
distilled water only, no nutrients). Do not apply
hydroponic solution to Type B foam.

If minimal moisture is applied to keep foam from
desiccating under LED heat: use plain distilled
water only. Apply sparingly.

**Step B6: Confirm and record**
For each Type B unit record:
- UV-A LED strip: product name, stated wavelength,
  date of purchase
- LED embedding depth: mm from seed-facing surface
- Foam moisture condition: dry / minimally moist
  with plain water only
- Scenario number assigned to this unit

---

### Type A Substrate — Hydroponic Foam, No UV

Type A carries water and minerals only. No UV.

**Step A1: Cut foam disc**
Cut a disc of 20ppi reticulated foam to fit snugly
in pipe inner diameter (50mm). Thickness: 20mm.
No LED channel required.

**Step A2: Mark seed-facing surface**
Mark seed-facing surface on outer edge with scenario
number and "SEED SIDE."

**Step A3: Saturate with nutrient solution**
Submerge foam disc in prepared hydroponic nutrient
solution. Saturate fully — approximately 5 minutes
submersion with gentle compression and release.
Remove and drain 2 minutes. Foam should be fully
saturated but not dripping.

**Step A4: Confirm and record**
For each Type A unit record:
- Nutrient solution: product name, batch, dilution
  ratio, date prepared
- Foam disc: ppi rating, thickness (mm), supplier
- Scenario number assigned to this unit

---

### Type C Substrate — Combined Hydroponic + UV Foam

Type C carries water, minerals, and UV-A all
co-located. Construction combines Type A and Type B
procedures.

**Step C1: Cut foam disc**
Cut a disc of 20ppi reticulated foam to fit snugly
in pipe inner diameter (50mm). Thickness: 30mm.

**Step C2: Mark seed-facing surface**
Mark seed-facing surface with scenario number and
"SEED SIDE."

**Step C3: Cut LED channel**
As per Step B3 — cut channel at 10–20mm depth from
seed-facing surface. Record depth.

**Step C4: Seat LED strip in channel**
As per Step B4 — seat LED strip flush in channel.
Route power cable out through notch in foam edge.

**Step C5: Saturate with nutrient solution**
As per Step A3 — submerge fully in hydroponic
nutrient solution including the LED-seated foam.
The LED strip is waterproof or water-resistant for
this purpose. Confirm LED strip is rated for moist
conditions before saturating. If LED strip is not
moisture-rated, apply a thin coat of clear silicone
to the strip surface before seating in channel and
allow to cure before saturation.

Saturate fully. Drain 2 minutes.

**Step C6: Confirm and record**
For each Type C unit record:
- UV-A LED strip: product name, stated wavelength,
  moisture rating, date of purchase
- LED embedding depth: mm from seed-facing surface
- Nutrient solution: product name, batch, dilution
  ratio, date prepared
- Scenario number assigned to this unit

---

### Free UV Condition — Scenario 1 S1 Only

**Step F1: Mount UV-A LED to S1 end cap**
Cut a short length of UV-A LED strip sufficient
to cover the inner face of the S1 end cap.
Adhere to the inner face of the S1 end cap using
the LED strip adhesive backing or a mounting
bracket. Orient the LED face downward into the
S1 chamber air space toward the seed.

Route the power cable through a small port in the
S1 end cap. Seal around the cable with aquarium
silicone.

**Step F2: Confirm and record**
- UV-A LED strip: product name, stated wavelength,
  date of purchase
- Mounting position: inner face of S1 end cap,
  facing downward

---

## PART V — CONE MEMBRANE CONSTRUCTION

Construct 5 cone membrane units — one per scenario.
All units are identical in construction. Orientation
is determined at assembly time.

### Step M1: Confirm funnel dimensions

Funnel rim diameter must match pipe inner diameter
(50mm) within ±1mm.

If rim is slightly too large: trim with scissors.
If rim is slightly too small: fill gap with
aquarium silicone at assembly.

### Step M2: Form the apex port

Apex port must match net cup outer diameter
(16mm or 25mm). Measure net cups before forming.

1. Heat pin or needle with lighter until glowing.
2. Melt small pilot hole at exact center of funnel
   apex.
3. Enlarge incrementally until net cup outer rim
   seats snugly in port.
4. Net cup rim rests on funnel surface at apex
   without falling through.
5. Test: insert net cup, invert funnel, confirm
   net cup does not fall out under own weight.

Failure mode: port too large — net cup falls through.
Recovery: discard funnel, begin again.

Failure mode: port off-center.
Recovery: discard funnel, begin again.

### Step M3: Wrap funnel in Mylar

Mylar provides two barrier functions:
1. UV-A and light barrier — aluminum layer reflects
   and absorbs UV-A and visible light, preventing
   any photonic signal crossing the membrane in
   either direction
2. Liquid barrier — prevents nutrient solution
   migration across the membrane

1. Cut Mylar to cover full outer cone surface with
   15mm overlap at all edges.
2. Apply thin layer of aquarium silicone to outer
   funnel surface.
3. Press Mylar onto silicone-coated surface from
   apex downward toward rim. Eliminate all air
   pockets and gaps.
4. At apex port: cut star pattern in Mylar (4–6
   radial cuts from center). Fold tabs through
   port opening and press flat against inner funnel
   surface. Secure with aquarium silicone ring
   on inside.
5. At rim: fold excess Mylar over rim edge and
   press flat against inner rim surface. Secure
   with aquarium silicone.
6. Cure minimum 12 hours.

Confirm after curing: no gaps in Mylar visible
when held to strong light source. Apex port
remains open. Net cup still seats correctly.

### Step M4: Light-seal test of membrane unit

Before installing in chamber, test for light
leakage:

1. In darkened room, hold assembled membrane unit
   over bright torch with wide opening facing
   light source.
2. Observe apex side for any light transmission
   through membrane surface or around rim.
3. Any visible light transmission: re-apply Mylar
   and silicone to gap. Cure again. Retest.

Do not install a membrane unit that fails this
test.

### Step M5: Determine orientation before sealing

| Scenario | Cone Orientation | Rationale |
|---|---|---|
| 1 | Apex up | Free UV above, Type A liquid below — apex up away from liquid in S2 |
| 2 | Apex up | Type B UV foam in S1 above — apex up, away from liquid threat above |
| 3 | Apex down | Type B UV foam in S2 below — apex down, away from liquid threat below |
| 4 | Apex up | Type C combined foam in S1 above — apex up, toward empty void in S2 |
| 5 | Apex down | Type C combined foam in S2 below — apex down, toward empty void in S1 |

Mark orientation on outside of chamber pipe with
permanent marker before sealing. Once silicone
cures, orientation cannot be changed.

### Step M6: Seal membrane into chamber

1. Insert cone membrane into pipe at membrane plane
   — junction between S1 and S2 pipe sections.
2. Apply continuous bead of aquarium silicone
   around full circumference of funnel rim on both
   S1 side and S2 side.
3. Cure minimum 12 hours.

Confirm: no gap between funnel rim and pipe wall
around full circumference. Perform light-seal test
on fully assembled joint — shine torch at joint
from outside, confirm no light visible from
interior.

---

## PART VI — VENTILATION PORT INSTALLATION

Install ventilation ports in all 5 chamber units
before loading substrates or seeds.

### Step V1: Drill ventilation ports

For each chamber unit, drill one ventilation port
in each end cap:
- Port diameter: 10–15mm
- Location: centered in end cap face

For active ventilation, drill a second port in the
same or opposite end cap.

### Step V2: Install foam light trap plugs

Cut cylindrical plugs of 20ppi reticulated foam:
- Diameter: equal to port hole diameter
- Length: 20–25mm

Press each plug firmly into its port hole. Apply
thin ring of aquarium silicone around plug
circumference at outer face of end cap to seal
any gap. Do not seal inner face — gas must exit
through plug.

Cure silicone minimum 12 hours.

### Step V3: Test each port for light exclusion

In darkened room, shine torch at outer face of
each foam plug. Observe inner side of end cap
for any visible light transmission.

Any light transmission: re-seat plug, add
silicone, cure, retest.

### Step V4: Install optional fan

Mount small USB fan at one end cap port using
aquarium silicone to seal between fan housing and
end cap. Fan draws ambient air through foam plug.

Confirm fan functioning before T=0. Record fan
model and voltage in experimental log.

Running fan 10–15 minutes at each observation
time point is sufficient for atmosphere refresh.

---

## PART VII — PRE-EXPERIMENT PREPARATION

### Step P1: Seed Viability Test
### (3–4 days before experiment day)

1. Place 10 seeds on damp paper towel folded in
   half.
2. Seal in ziplock bag with small air pocket.
3. Place in warm location (18–24°C) in complete
   darkness.
4. Check at 24, 48, and 72 hours under red-light
   only.
5. Record: number germinated, time to germination,
   radicle length at each check.

Proceed only if germination rate ≥ 80%.

Failure mode: germination rate below 80%.
Recovery: source new seeds from different batch.

### Step P2: Prepare Hydroponic Nutrient Solution

1. Measure 1L distilled water.
2. Add hydroponic nutrient concentrate at seedling
   dilution (1/4 to 1/2 strength per manufacturer
   instructions).
3. Mix thoroughly.
4. Record: product name, batch number, dilution
   ratio, date prepared.
5. Use within 48 hours of preparation.
6. Use same solution batch for all Type A and
   Type C units in a single experimental run.

### Step P3: Construct All Substrate Units

Following Part IV, construct:
- 1 Free UV condition unit (Scenario 1 S1)
- 2 Type B units (Scenarios 2 and 3)
- 2 Type C units (Scenarios 4 and 5)
- 3 Type A units (Scenarios 1, 2, and 3)

Record all construction parameters. Leave LED
cables unconnected until substrate loading.

### Step P4: Construct All 5 Cone Membrane Units

Following Part V, construct and cure all 5 cone
membrane units. Allow minimum 12 hours curing.
Perform light-seal test on all 5 units before
installation. Label each by scenario number and
cone orientation.

### Step P5: Assemble and Cure Chamber Units

For each scenario, join S1 and S2 pipe sections
at cured cone membrane. Seal with aquarium
silicone. Allow minimum 12 hours curing.

Install ventilation ports per Part VI. Test all
ports for light exclusion.

Confirm all 5 chamber units are:
- Sealed at all joints
- Ventilation ports installed and light-tested
- External light exclusion confirmed or enclosure
  ready
- Labelled by scenario number, cone orientation,
  S1 contents, S2 contents

### Step P6: Full System Light Exclusion Test

Before loading any substrates or seeds, in a
fully darkened room, allow eyes to adapt minimum
3 minutes. Examine each chamber unit for any
visible light leakage at joints, ventilation
ports, end cap edges, and pipe walls.

Note: Scenario 1 S1 section is exempt from this
test — it will contain a freely emitting UV-A LED.
Test only the S2 section of Scenario 1.

Correct any leakage before proceeding.

---

## PART VIII — SEED GERMINATION AND PLACEMENT

### Step G1: Germinate Seeds
### (48 hours before placement)

1. Place 10 seeds on damp paper towel.
2. Fold paper towel over seeds.
3. Seal in ziplock bag with small air pocket.
4. Place in warm location (20–24°C) in complete
   darkness.
5. Keep in complete darkness throughout germination
   — seeds must enter chambers in dark-adapted
   skotomorphogenic state.
6. Check at 24 and 48 hours under red-light only.

Target state for placement:
- Radicle visible and clearly directed
- Radicle length: 3–8mm
- Hypocotyl not yet significantly elongated
- Seed coat may still be partially attached

Do not wait for radicle to exceed 10mm.

### Step G2: Select Seeds for Placement

Select 5 seeds — one per scenario.

Selection criteria:
- Clear radicle direction visible
- Radicle intact, no kinking or mechanical damage
- Radicle length 3–8mm
- No signs of fungal contamination

Record for each selected seed:
- Radicle length at placement (mm)
- Hours since germination began

All selection and placement steps under red-light
only. Do not expose germinated seeds to white
light, blue light, or UV at any point.

### Step G3: Place Germinated Seed in Net Cup

Using tweezers and clean gloves, under red-light:

1. Handle seed by seed body only. Do not grasp
   radicle.
2. Lower seed into net cup with radicle pointing
   downward through net cup mesh bottom.
3. Seed body rests on net cup rim or upper mesh.
4. Radicle extends downward through mesh without
   compression or bending.
5. Stabilize with small moist foam fragment inside
   net cup if needed. Do not pack tightly.

Confirm: radicle exits downward through mesh
without bending. Seed body stable, does not rock.

### Step G4: Seat Net Cup in Cone Membrane

Under red-light only:

1. Lower net cup with seed into apex port of cone
   membrane.
2. Net cup rim seats on cone surface at apex.
3. Radicle extends through port into substrate
   space on apex side of membrane.
4. Seed body and emerging hypocotyl are in
   substrate space on wide-opening side of membrane.

**Apex-up scenarios (1, 2, 4):**
- Radicle extends downward into S2
- Seed body and hypocotyl in S1 space

**Apex-down scenarios (3, 5):**
- Cone inverted — wide opening faces upward
- Net cup seats at bottom apex
- Radicle extends downward through port into S2

Confirm under red-light: net cup seated snugly,
no rocking, radicle not compressed against port
edge.

### Step G5: Load Substrates and Activate Signals

This is T=0. Record exact time for each scenario.
All 5 scenarios loaded and activated within 30
minutes maximum. Perform all steps under red-light.

---

**Scenario 1:**
- S1: Mount Free UV unit — UV-A LED strip adhered
  to inner face of S1 end cap, facing downward.
  Connect LED to power. Confirm UV-A active.
  Seal S1 end cap with ventilation port. Route
  LED cable through sealed cable port.
  Note: S1 of Scenario 1 is NOT externally wrapped
  — free UV in air is the condition.
  S2 end cap and pipe body should be opaque or
  wrapped to prevent external light entering S2.
- S2: Insert Type A substrate unit, seed-facing
  surface oriented upward toward seed. Confirm
  foam contacts S2 inner surface. Seal S2 end cap
  with ventilation port.

---

**Scenario 2:**
- S1: Insert Type B substrate unit, seed-facing
  surface oriented downward toward seed. Connect
  UV-A LED strip to power. Confirm UV-A active.
  Seal S1 end cap with ventilation port. Route
  LED cable through sealed cable port.
- S2: Insert Type A substrate unit, seed-facing
  surface oriented upward toward seed. Seal S2
  end cap with ventilation port.
- Apply external light exclusion wrap to full
  chamber immediately after sealing.

---

**Scenario 3:**
- S1: Insert Type A substrate unit, seed-facing
  surface oriented downward toward seed. Seal S1
  end cap with ventilation port.
- S2: Insert Type B substrate unit, seed-facing
  surface oriented upward toward seed. Connect
  UV-A LED strip. Confirm active. Seal S2 end cap
  with ventilation port. Route LED cable.
- Apply external light exclusion wrap to full
  chamber immediately after sealing.

---

**Scenario 4:**
- S1: Insert Type C substrate unit, seed-facing
  surface oriented downward toward seed. Connect
  UV-A LED strip. Confirm active. Seal S1 end cap
  with ventilation port. Route LED cable.
- S2: Empty void. Seal S2 end cap with ventilation
  port.
- Apply external light exclusion wrap to full
  chamber immediately after sealing.

---

**Scenario 5:**
- S1: Empty void. Seal S1 end cap with ventilation
  port.
- S2: Insert Type C substrate unit, seed-facing
  surface oriented upward toward seed. Connect
  UV-A LED strip. Confirm active. Seal S2 end cap
  with ventilation port. Route LED cable.
- Apply external light exclusion wrap to full
  chamber immediately after sealing.

---

Record T=0 for each scenario in experimental log.
Record LED activation confirmed Y/N for all
scenarios containing UV substrate or Free UV.

Place all Scenario 2 through 5 chamber units in
dark enclosure or confirm external wrap fully
applied. No further light exposure to these
chambers until T=24 observation.

---

## PART IX — OBSERVATION PROTOCOL

All observations at specified time points. Between
time points, Scenarios 2 through 5 are in complete
darkness with no exceptions.

Equipment ready before opening any chamber:
- Camera charged and ready
- Written log open
- Millimeter ruler
- Red-light torch
- Fresh wrap material ready for re-sealing

---

### Substrate Penetration Depth Measurement

Definition: distance (mm) from seed-facing surface
of foam substrate to visible tip of root or shoot
growing into that substrate.

Measured by estimating position of growing tip
through pipe wall or via brief direct observation
under red-light when end cap is removed.

Record:
- Shoot penetration depth into S1 substrate (mm)
  — if S1 contains foam substrate
- Root penetration depth into S2 substrate (mm)
  — if S2 contains foam substrate
- N/A for empty chambers and Free UV condition

When penetration depth exceeds LED embedding depth
in Type B or Type C substrate, the growing tip has
reached the UV-A signal. Record this crossing event
with the time point it was first observed.

---

### T=24 Hours — First Germination Check

For each scenario record:
- Radicle visible and directed as at placement? Y/N
- Radicle elongated since T=0? Estimated length (mm)
- Root penetration depth into substrate (mm)
- Hypocotyl beginning to emerge? Y/N
- Shoot penetration depth into substrate (mm) if
  applicable
- Contamination visible? Y/N
- UV-A LED active? Y/N
- Photograph each chamber
- Re-seal and re-wrap Scenarios 2–5 immediately
  after photography

### T=48 Hours — Directional Growth Assessment

For each scenario record:
- Radicle length (mm)
- Radicle direction: downward / upward / horizontal
  / indeterminate
- Root penetration depth into substrate (mm)
- Root tip passed LED depth in Type B or C?
  Y/N/N/A — record time of crossing if first
  observed here
- Hypocotyl direction: upward / downward /
  horizontal / indeterminate
- Shoot penetration depth into substrate (mm)
- Shoot tip passed LED depth in Type B or C?
  Y/N/N/A — record time of crossing if first
  observed here
- Architectural direction determinable? Y/N
- UV-A LED active? Y/N
- Photograph from side and above
- Re-seal and re-wrap Scenarios 2–5 immediately

### T=72 Hours — Primary Outcome Point

For each scenario record:
- Radicle length (mm)
- Radicle direction. Angle from vertical (degrees).
  0° = straight down. 180° = straight up.
- Root penetration depth into substrate (mm)
- Root tip passed LED depth? Y/N/N/A
- Hypocotyl length (mm)
- Hypocotyl direction. Angle from vertical (degrees)
- Shoot penetration depth into substrate (mm)
- Shoot tip passed LED depth? Y/N/N/A
- Cotyledon emergence: Y/N. Direction if emerged.
- UV-A LED active? Y/N
- Contamination: Y/N. Location and description.
- Overall architectural assessment:
  NORMAL / INVERTED / INDETERMINATE / FAILED
- Photograph from side, from above, close-up of
  membrane port region
- Re-seal and re-wrap Scenarios 2–5 immediately

**Primary result determination at T=72:**

| Scenario | NORMAL | INVERTED | INDETERMINATE |
|---|---|---|---|
| 1 | Root down into S2 Type A, shoot up into free UV air | Root up, shoot down | Neither clear |
| 2 | Root down into S2 Type A, shoot up penetrating S1 Type B | Root up, shoot down | Neither clear |
| **3** | Root down, shoot up | **Root up penetrating S1 Type A, shoot down penetrating S2 Type B** | Neither clear |
| 4 | Root down into void, shoot up penetrating S1 Type C | Root down, shoot down | Neither clear |
| 5 | Root up into void, shoot down penetrating S2 Type C | Root down, shoot up | Neither clear |

**The experiment succeeds if Scenario 3 records
INVERTED at T=72 with measurable root penetration
into S1 Type A foam and measurable shoot penetration
into S2 Type B foam.**

Scenario 1 recording NORMAL is required as the
control confirmation in every run.

### T=120 Hours (Day 5) and T=168 Hours (Day 7)

Repeat full T=72 observation protocol at T=120
and T=168.

Purpose: confirm directional growth is maintained
and not corrected by gravitropic override over time.
Track continued substrate penetration depth.

If Scenario 3 shows INVERTED at T=72 but NORMAL
at T=120, gravitropic correction has occurred.
Record as qualified result — report full angle data
at all time points. This is data about relative
signal strength versus gravitropism, not framework
failure.

---

## PART X — DATA RECORDING TEMPLATE

Complete for each scenario at each time point:

```
SCENARIO: ___
TIME POINT: T= ___ hours
DATE/TIME OF OBSERVATION: ___
OBSERVER INITIALS: ___
DURATION OF LIGHT EXPOSURE THIS TIME POINT (minutes): ___

RADICLE / ROOT:
  Length (mm): ___
  Direction: downward / upward / horizontal /
             indeterminate
  Angle from vertical (degrees): ___
  Penetration depth into substrate (mm): ___
  Substrate type being penetrated:
    Type A / Type B / Type C / void / free UV air
  Root tip has passed LED embedding depth?
    Y / N / N/A
  If Y, first observed at this time point? Y / N
  Notes: ___

HYPOCOTYL / SHOOT:
  Length (mm): ___
  Direction: upward / downward / horizontal /
             indeterminate
  Angle from vertical (degrees): ___
  Penetration depth into substrate (mm): ___
  Substrate type being penetrated:
    Type A / Type B / Type C / void / free UV air
  Shoot tip has passed LED embedding depth?
    Y / N / N/A
  If Y, first observed at this time point? Y / N
  Notes: ___

COTYLEDONS:
  Emerged: Y / N
  Direction if emerged: ___

UV-A LED STATUS:
  Confirmed active: Y / N
  Method of confirmation: ___

CONTAMINATION:
  Observed: Y / N
  Location: seed / root / shoot / S1 substrate /
            S2 substrate / chamber wall
  Description: ___

OVERALL ARCHITECTURAL ASSESSMENT:
  NORMAL / INVERTED / INDETERMINATE / FAILED

LIGHT EXCLUSION (Scenarios 2–5):
  External wrap / enclosure intact before opening:
    Y / N
  Re-wrapped / re-enclosed immediately after: Y / N
  Any accidental light exposure: Y / N
  Description if yes: ___

PHOTOGRAPHS TAKEN: Y / N
  Filename(s): ___

INTERVENTIONS: Y / N
  Description if yes: ___
  Data flagged as potentially compromised: Y / N
```

---

## PART XI — REPRODUCIBILITY REQUIREMENTS

1. **Scenario 3 INVERTED result must be observed
   in minimum 3 independent seeds across minimum
   2 independent experimental runs.**

2. **Scenario 1 must record NORMAL in every run
   where Scenario 3 records INVERTED.** If Scenario
   1 fails to record NORMAL, the run is invalid.

3. **Measurable shoot penetration into S2 Type B
   foam must be recorded in Scenario 3 at T=72.**
   Directional growth alone without penetration
   does not confirm active navigation toward the
   UV signal.

4. **Complete light exclusion must be confirmed
   and logged for Scenarios 2 through 5 at every
   time point.** Runs with compromised light
   exclusion must be flagged and reported
   separately from primary results.

5. **All 5 scenarios must be run simultaneously
   in each experimental run** from the same seed
   batch, same substrate preparation batch, and
   same UV-A LED product batch.

6. **Photographic record at minimum T=0, T=24,
   T=48, T=72, T=120, T=168 for all 5 scenarios.**

7. **All materials, batch numbers, UV-A wavelength
   confirmation, LED embedding depth, preparation
   dates, light exposure durations, and all
   observation records must be logged** per Part X
   template.

8. **Any intervention must be recorded** and that
   chamber's data flagged as potentially
   compromised.

---

## PART XII — FAILURE MODES AND RECOVERY ACTIONS

| Failure Mode | Detection Point | Recovery Action |
|---|---|---|
| Seed does not germinate | T=24h, T=48h | Check temperature (18–24°C) and moisture. Replace with spare seed before T=48h if available. Record replacement. |
| Radicle compressed in port | T=24h | Attempt gentle repositioning under red-light with clean pin. Record intervention. Flag data. |
| Fungal contamination in substrate | Any time point | If substrate only: continue and note. If on seed or radicle: flag data, continue to T=72h, report contamination in results. |
| UV-A LED failure | Any time point | Replace LED under red-light. Record replacement and duration of outage. If outage exceeded 12 hours flag run for that scenario. |
| Light exclusion breach in Scenarios 2–5 | Any time point | Record duration and source of accidental exposure. Flag affected scenario data as potentially compromised. Continue run. Report in results. |
| Silicone seal failure — liquid migration between chambers | Any time point | Run compromised for that scenario. Complete observation to T=72h. Report seal failure. Discard that scenario's data from primary results. Rebuild membrane for next run. |
| Shoot does not penetrate Type B or C foam at T=72h | T=72h | Record penetration depth as zero or near-zero. Valid result — signal may be insufficient to drive penetration at current LED depth or intensity. Consider reducing embedding depth or increasing LED intensity for next run. Report as is. |
| Gravitropic correction — INVERTED at T=72h, NORMAL at T=120h | T=120h | Record as qualified result. Report full angle data at all time points. Data about competition between coherence gradient and gravitropism. Not a framework failure. |
| Scenario 1 records NORMAL but Scenario 3 records NORMAL also | T=72h | Primary claim not confirmed in this run. Check light exclusion integrity for Scenario 3. Check Type B LED was active. Assess whether UV signal was reaching seedling. Redesign or repair for next run. |
| All 5 scenarios record NORMAL | T=72h | Either light exclusion was compromised in Scenarios 2–5 or signal strength was insufficient. Verify LED wavelength and function. Verify light exclusion was intact. Increase signal strength or improve exclusion for next run. |
| Active ventilation fan failure | Any time point | Check USB connection. Replace fan if failed. If fan off for extended period in sealed chamber, open ventilation port briefly under red-light for gas exchange before resealing. Record. |

---

## DOCUMENT METADATA

- Author: Eric Robert Lawson
- Affiliation: OrganismCore (Independent Research)
- ORCID: 0009-0002-0414-6544
- Contact: OrganismCore@proton.me
- Repository:
  github.com/Eric-Robert-Lawson/attractor-oncology
- Document version: 2.3 (didn't upload other drafted version with this)
- Supersedes:
  Plant_Inversion_Operative_Protocol.md (v1.0)
- Timestamp: 2026-03-18
- Pre-registration DOI: 10.5281/zenodo.18986790
- License: CC BY-NC-SA 4.0
