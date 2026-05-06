# EHAB Prototype Protocol
## Encapsulating Hydrogel-Aerogel Barrier (EHAB) — From-Scratch Fabrication, Costing, Testing, and Iterative Improvement
## OrganismCore Application — Practical Prototyping Record
## Date: 2026-05-06
## Author: Eric Robert Lawson / Copilot Synthesis Agent

---

## DOCUMENT STATUS

```
Type:         Practical prototyping protocol — step 1 to finish
Audience:     Prototyper with access to basic lab or workshop equipment
              No specialized equipment required for Phase 1 prototype
Predecessor:  encapsulating_gel_suppression_causal_geometry_derivation.md
              gel_substance_deployment_gradient_analysis.md
Scope:        (1) Ingredient sourcing and pricing
              (2) Step-by-step synthesis protocol
              (3) Testing battery — surface adhesion, fire protection,
                  aerogel transition, deployability
              (4) Predictable failure modes and correction protocol
Status:       First-principles prototyping guide
              Coherence ≠ truth — all outputs require empirical validation
Version:      1.0
```

---

## PART I: INGREDIENT SOURCING AND PRICING

### 1.1 The Component List

The EHAB material is built from five component layers. Each can be sourced
independently. The Phase 1 prototype uses Layers 1-3 only (base system).
Layers 4 and 5 are tested in Phase 2 once the base system is validated.

---

### Layer 1 — Superabsorbent Polymer (SAP) — The Gel Matrix

**What it is:** Crosslinked sodium polyacrylate powder. Absorbs 300-500x its
weight in water. Forms the sticky, viscous gel that coats surfaces.

**Sourcing options (ordered by accessibility):**

| Source | Product | Cost | Notes |
|---|---|---|---|
| Amazon / craft supply | Sodium polyacrylate SAP crystals (sold as "water crystals", "fake snow", "soil crystals") | $8-15 / 100g | Immediate. Fine for Phase 1 prototype. Particle size variable — use fine powder grade if available. |
| Alibaba / industrial supplier (MOQ 1kg) | Industrial-grade SAP powder | $2-5/kg | Best for Phase 2 and beyond. 2-4 week lead time. |
| Sigma-Aldrich / lab supplier | Sodium polyacrylate, reagent grade | $50-80/100g | Unnecessary for this application — use industrial or consumer grade. |

**Phase 1 recommendation:** Buy consumer "water absorbing crystals" or "SAP
powder" on Amazon for initial testing. 100g is sufficient for 20+ test batches.
**Phase 1 cost: ~$10-15 for 100g**

---

### Layer 2 — Colloidal Silica / Fumed Silica — The Aerogel Transition Agent

**What it is:** Fine silica nanoparticles that self-assemble into aerogel structure
as gel water evaporates under heat. This is the Phase C protection agent.

**Critical note for Phase 1:** Two forms available:
- **Fumed silica (powder):** Aerosil 200, Cab-O-Sil — dry powder, disperses in water.
  Easier to source and handle. Fine for Phase 1.
- **Colloidal silica (suspension):** Pre-dispersed nanoparticles in water.
  More expensive but easier to mix uniformly. Use in Phase 2.

| Source | Product | Cost | Notes |
|---|---|---|---|
| Amazon / pottery supply | Fumed silica / "Cab-O-Sil" (sold as matte agent or epoxy thickener) | $15-30 / 200g | Immediate. Aerosil 200 equivalent. WEAR RESPIRATOR — fine silica dust is hazardous to inhale. |
| Aliexpress / industrial | Fumed silica powder, hydrophilic | $3-8/kg (bulk) | Best pricing. 2-4 week lead time. |
| Sigma-Aldrich | Fumed silica, 7nm BET | $45/100g | Fine for Phase 1 lab use if already ordering from them. |

**Phase 1 recommendation:** Source fumed silica from Amazon (sold as epoxy
thickener or matte agent) or pottery supply. 200g is more than sufficient.
**Phase 1 cost: ~$20-30 for 200g**

**⚠️ SAFETY:** Fumed silica is an extremely fine powder. Inhaling it causes
silicosis. Always handle under a fume hood or outdoors with an N95 respirator.
Never handle in an enclosed space without respiratory protection.

---

### Layer 3 — Clay Filler — The Rheology and Barrier Agent

**What it is:** Bentonite or nanoclay platelets that give the gel thixotropic
(shear-thinning) behavior and create a gas-barrier effect in the film.

**Phase 1 approach:** Two options with very different costs:

| Source | Product | Cost | Notes |
|---|---|---|---|
| Pottery/ceramics supply, Amazon | Sodium bentonite clay powder | $5-15/kg | Immediate. Coarser particle size but works for Phase 1 baseline. |
| Bulk apothecary / soap supply | Bentonite clay (cosmetic grade) | $8-15/kg | Fine particle size, clean, good for consistent results. |
| Sigma-Aldrich | Nanoclay (Laponite equivalent), Cloisite | $80-200/100g | Lab grade. Better performance. Consider for Phase 2 when optimizing. |
| Alibaba | Bentonite nanoclay, montmorillonite | $0.5-2/kg (bulk, MOQ 10kg) | Best option once scaling. |

**Phase 1 recommendation:** Cosmetic-grade bentonite clay from soap/craft supply.
500g is sufficient for all Phase 1 tests.
**Phase 1 cost: ~$10-15 for 500g**

---

### Layer 4 — CO₂ Foaming Agent (Phase 2) — The Thermal Activation Component

**What it is:** Two dry ingredients stored separately, encapsulated in wax
microcapsules that break at ~80-100°C, releasing CO₂ within the gel film
to create an in-place foam and O₂ displacement effect.

| Source | Product | Cost | Notes |
|---|---|---|---|
| Grocery / supermarket | Baking soda (sodium bicarbonate, food grade) | $1-3/kg | Immediately available everywhere. |
| Grocery / supplement store | Citric acid powder (food grade) | $3-6/kg | Immediately available — brewing supply or health food store. |
| Craft supply / candle supply | Paraffin wax pellets | $5-10/kg | For wax microencapsulation shell. Melting point ~55-65°C. |

**Phase 1:** Mix baking soda and citric acid UNENCAPSULATED first to confirm
CO₂ generation and foaming behavior. Encapsulation is a Phase 2 refinement.

**Phase 2:** Wax microencapsulation (see Section 3.4 below).
**Phase 1 cost: ~$10-15 for all three components**

---

### Layer 5 — PCM Microcapsules (Phase 3 Enhancement — Optional)

**Not required for Phase 1 or Phase 2 prototype. Skip unless specifically
testing extended protection window scenarios.**

| Source | Product | Cost | Notes |
|---|---|---|---|
| Alibaba / specialty chemical | Microencapsulated paraffin PCM, Tm ~50-60°C | $15-30/kg (bulk) | 2-4 week lead time. |
| Research supplier | Micronal (BASF), Microtek PCM | $30-80/kg | Better quality control for research. |

---

### 1.2 Phase 1 Prototype — Complete Shopping List and Budget

```
PHASE 1 PROTOTYPE — COMPLETE INGREDIENT LIST

Item                          Qty needed    Source              Cost (est.)
──────────────────────────────────────────────────────────────────────────
SAP powder (sodium            100g          Amazon / craft       $10-15
polyacrylate, fine grade)

Fumed silica (Aerosil         200g          Amazon / pottery/    $20-30
equivalent, hydrophilic)                    epoxy supply

Bentonite clay powder         500g          Soap/cosmetic        $10-15
(sodium bentonite,                          supply / Amazon
cosmetic or pottery grade)

Baking soda                   500g          Grocery              $2-3
(sodium bicarbonate)

Citric acid powder            250g          Grocery /            $3-5
(food grade)                                brewing supply

Paraffin wax pellets          200g          Candle supply /      $5-8
                                            Amazon

Distilled water               5 liters      Grocery /            $3-5
                                            auto parts store

──────────────────────────────────────────────────────────────────────────
TOTAL PHASE 1 INGREDIENT COST:                                   $53-81

EQUIPMENT (if not already owned):
Mixing bowls (glass or       2-3            Dollar store /       $5-10
stainless steel)                            kitchen supply

Kitchen scale (0.1g          1              Amazon               $15-25
resolution)

N95 respirator               1 pack         Hardware store       $10-15
(for fumed silica handling)

Nitrile gloves               1 box          Hardware store       $8-12

Heat-resistant sample        6-10 pieces    Hardware store /     $10-20
holders (metal cookie        (10x10cm)      kitchen supply
sheets or metal plate
offcuts)

Propane torch (small         1              Hardware store       $15-25
handheld "Bernzomatic"
style)

Stopwatch or phone           1              Already owned        $0

Ruler, marker, camera        standard       Already owned        $0

──────────────────────────────────────────────────────────────────────────
TOTAL PHASE 1 ALL-IN COST:                                       $116-188
(ingredients + basic equipment)
```

---

### 1.3 Per-Liter Production Cost at Phase 1 Scale

```
Phase 1 batch (1 liter gel):

SAP powder: 5g × ($0.15/g) = $0.75
Fumed silica: 25g × ($0.15/g) = $3.75
Bentonite: 30g × ($0.02/g) = $0.60
Water: 940mL (negligible)
──────────────────────────────
Phase 1 cost per liter: ~$5.10

Note: At industrial scale (bulk component pricing):
SAP: 5g × ($0.003/g) = $0.015
Fumed silica: 25g × ($0.006/g) = $0.15
Bentonite: 30g × ($0.001/g) = $0.03
Water: negligible
Processing: ~$0.30/liter estimate
──────────────────────────────
Industrial-scale cost per liter: ~$0.50
(Order-of-magnitude estimate only)
```

---

## PART II: STEP-BY-STEP SYNTHESIS PROTOCOL

### Overview of the Protocol

```
PHASE 1: BASE GEL (Layers 1-3)
  Step 1: SAP hydration — prepare the gel base
  Step 2: Fumed silica dispersion — prepare the aerogel agent
  Step 3: Bentonite dispersion — prepare the rheology/barrier agent
  Step 4: Component compounding — combine into final gel
  Step 5: Quality check — viscosity, adhesion, appearance

PHASE 2: FOAMING ACTIVATION (Layer 4)
  Step 6: Prepare unencapsulated activation test (baking soda + citric acid)
  Step 7: Wax microencapsulation protocol
  Step 8: Incorporate encapsulated agent into gel

PHASE 3: FORMULATION OPTIMIZATION
  Step 9: Systematic variation trials
  Step 10: Record, compare, converge on optimal formulation
```

---

### Step 1: SAP Hydration — Preparing the Gel Base

**Purpose:** Create the base hydrogel that will carry all other components
and provide surface adhesion.

**You will make:** A stock hydrated SAP gel at 1% w/v concentration.

```
MATERIALS:
  SAP powder (sodium polyacrylate): 10g
  Distilled water: 1000 mL (1 liter)
  Mixing bowl (large — the gel will expand significantly)
  Kitchen scale
  Spatula or whisk

PROCEDURE:

1. Weigh 10g of SAP powder into the large mixing bowl.
   (Note: 10g of SAP will absorb approximately 3-5 liters of water.
   Start with 1 liter — the gel will not be fully saturated, which is
   intentional. You want a flowable gel, not a solid block.)

2. Slowly add 1000 mL of distilled water to the SAP powder while stirring
   continuously. Add water in 100mL increments, stirring between additions.
   
   WHY: Adding water too fast creates lumps that are hard to disperse. Slow
   addition ensures uniform hydration.

3. Continue stirring for 5-10 minutes after all water is added.
   The gel will thicken progressively.

4. Allow to rest for 30 minutes. The SAP will fully hydrate and form a
   clear-to-slightly-hazy viscous gel.

5. Observe: The gel should hold its shape when scooped but flow slowly
   under its own weight. If it is a solid block, add more water (50mL at
   a time). If it is completely liquid, add more SAP (1g at a time).

TARGET APPEARANCE:
  Clear to hazy, viscous gel
  Consistency similar to: thick hair conditioner or aloe vera gel
  Pours slowly, does not drip freely
  
TROUBLESHOOTING:
  Too thick / lumpy: Add water in 50mL increments, mix thoroughly between additions
  Too thin / watery: Add SAP 1g at a time, wait 10 minutes between additions
  
RECORD:
  Final mass of SAP used: ____g
  Final volume of water used: ____mL
  Visual assessment of consistency: ____
```

---

### Step 2: Fumed Silica Dispersion — Preparing the Aerogel Agent

**Purpose:** Create a uniform dispersion of fumed silica nanoparticles in
water that can be incorporated into the gel. Fumed silica in dry form will
not disperse easily — it must be wetted properly.

**⚠️ CRITICAL SAFETY STEP:** Perform this step outdoors or under a fume hood.
Wear N95 respirator and nitrile gloves. Never open the fumed silica container
without respiratory protection. The particles are invisible and will lodge
in lung tissue.

```
MATERIALS:
  Fumed silica powder: 25g
  Distilled water: 200 mL
  Mixing bowl (small, stainless steel or glass)
  Whisk or immersion blender
  N95 respirator (MANDATORY)
  Nitrile gloves (MANDATORY)

PROCEDURE:

1. Put on N95 respirator before opening the fumed silica container.
   Keep respirator on until the silica is fully incorporated in water
   and no dry powder remains visible.

2. Weigh 25g of fumed silica powder.

3. Add 200 mL of distilled water to the mixing bowl.

4. SLOWLY sift the fumed silica into the water while stirring vigorously.
   Add it in four or five portions, stirring for 1 minute between additions.
   
   WHY: Fumed silica is hydrophobic in its dry state and will initially
   float on the water surface. High shear mixing is required to force it
   to wet. If added all at once it forms a floating raft that is very
   difficult to disperse.

5. Once all powder is added and roughly wetted, mix vigorously with a whisk
   or immersion blender for 3-5 minutes.
   
   TARGET: A milky white, uniformly cloudy suspension with no floating
   white powder. It will feel slightly slippery/thickened.

6. Allow to rest for 10 minutes, then check for dry clumps at the surface.
   If present, continue mixing.

7. The silica dispersion is now ready to be incorporated in Step 4.

TARGET APPEARANCE:
  Milky white suspension (like skim milk or diluted paint)
  No dry white powder floating
  Uniform cloudiness throughout

TROUBLESHOOTING:
  Powder won't wet (floating raft): Add a single drop of dish soap as wetting
  agent, then mix. This introduces a trace surfactant which aids initial wetting.
  The small amount will not significantly affect gel performance.
  
  Clumps persist: Immersion blender or higher-speed mixing for 5+ minutes.

RECORD:
  Final mass of silica used: ____g
  Visual assessment: ____
```

---

### Step 3: Bentonite Dispersion — Preparing the Rheology Agent

**Purpose:** Create a bentonite clay dispersion that will give the combined gel
thixotropic behavior (flows under stress, holds at rest).

```
MATERIALS:
  Bentonite clay powder: 30g
  Distilled water: 500 mL
  Mixing bowl (medium)
  Whisk

PROCEDURE:

1. Add 500 mL of distilled water to the mixing bowl.

2. Sift bentonite powder into the water slowly while whisking.
   Add in 5-6 portions, whisking for 1 minute between additions.

3. Once all bentonite is added, whisk vigorously for 5 minutes.

4. Allow to hydrate and swell for 30-60 minutes undisturbed.
   Bentonite forms a gel structure as the clay platelets hydrate
   and stack into networks.

5. After resting, the bentonite dispersion should be noticeably thicker
   and exhibit the following behavior: if you stir it quickly and then
   leave it undisturbed, it should thicken back up within 1-2 minutes.
   This is thixotropic behavior — confirmation you are on the right track.

TARGET APPEARANCE:
  Grey-tan cloudy suspension, medium viscosity
  Stirs easily, thickens when left undisturbed
  
TROUBLESHOOTING:
  No thickening observed after 60 min: Add another 10g of bentonite,
  mix thoroughly, wait another 30 min.
  
  Too thick to stir: Add water 50mL at a time.

RECORD:
  Final mass of bentonite used: ____g
  Thixotropy observed? Y / N
  Visual assessment: ____
```

---

### Step 4: Component Compounding — Assembling the Final Gel

**Purpose:** Combine the SAP gel, silica dispersion, and bentonite dispersion
into the final EHAB base gel (Phase 1 formulation).

```
MATERIALS:
  SAP gel from Step 1 (all of it, ~1000g total weight)
  Fumed silica dispersion from Step 2 (~225g total)
  Bentonite dispersion from Step 3 (~530g total)
  Large mixing bowl
  Spatula or immersion blender

PHASE 1 TARGET FORMULATION (1750g total gel):
  Component                     Weight (g)    % of total
  ─────────────────────────────────────────────────────
  SAP hydrated gel              1000g         57%
  Fumed silica dispersion        225g         13%
  Bentonite dispersion           530g         30%
  ─────────────────────────────────────────────────────
  TOTAL                         1755g         ~1.75 liters

PROCEDURE:

1. Pour the SAP gel into the large mixing bowl.

2. Add the bentonite dispersion to the SAP gel. Mix thoroughly with spatula
   or immersion blender for 3 minutes.
   
   WHY: Bentonite goes in before silica because it is a structural network
   component. The clay platelets will integrate with the SAP polymer chains.

3. Slowly add the fumed silica dispersion to the SAP+bentonite mixture
   while mixing continuously. Add in four portions, mixing 2 minutes
   between additions.
   
   WHY: Fumed silica added slowly to maintain uniform distribution.
   Localized high-concentration silica creates lumps that are difficult
   to re-disperse.

4. Once all components are combined, mix thoroughly for 5 minutes.

5. Final mixing: Use the immersion blender at medium speed for 2-3 minutes
   for maximum uniformity. If using hand mixing only, extend to 10 minutes
   of vigorous stirring.

6. Allow the gel to rest undisturbed for 30 minutes.

7. After resting, perform the quick quality checks in Step 5 before
   proceeding to testing.

TARGET APPEARANCE:
  Milky white to pale grey viscous gel
  Significantly thicker than any individual component
  Thixotropic: stirs with some effort, holds shape when left undisturbed
  Should coat a spatula and hold a peak for 5+ seconds before flowing

RECORD:
  Total final gel mass: ____g
  Total volume (approximate): ____mL
  Consistency description: ____
  Thixotropy: flows when stirred Y/N, holds when still Y/N
```

---

### Step 5: Quality Check — Base Gel Properties

Before testing fire performance, confirm the base gel has the required
physical properties. These checks require no specialized equipment.

```
CHECK 1 — VERTICAL ADHESION TEST

Procedure:
  (a) Apply a 5cm wide, ~3mm thick stripe of gel to a vertical surface
      (glass, tile, painted drywall, or plywood).
  (b) Start timer.
  (c) Observe at 1 min, 5 min, 15 min, 30 min intervals.
  (d) Photograph at each interval.

Pass criteria:
  PASS: Gel holds on vertical surface for >15 minutes with <20% runoff
  FAIL: Gel runs off within 5 minutes (too thin — increase SAP or bentonite)

Record:
  % runoff at 15 min: ____%
  % runoff at 30 min: ____%
  PASS / FAIL: ____
  Formulation notes: ____

──────────────────────────────────────────────────────────────────────

CHECK 2 — SHEAR-THINNING CONFIRMATION

Procedure:
  (a) Stir gel vigorously with spatula for 30 seconds.
  (b) Remove spatula and leave gel undisturbed.
  (c) Observe and time how long it takes to stop flowing and hold a peak.
  (d) Try to pour gel from bowl — observe flow rate vs. at rest state.

Pass criteria:
  PASS: Gel flows when stirred/poured, slows and holds within 2 minutes at rest
  FAIL A: Gel is too rigid — does not flow even under stirring (need to add more water)
  FAIL B: Gel has no recovery — stays liquid after stirring (need more SAP or clay)

Record:
  Time to recovery (stop flowing after stirring): ____min
  Flow behavior when poured: ____
  PASS / FAIL: ____

──────────────────────────────────────────────────────────────────────

CHECK 3 — FILM FORMATION ON HORIZONTAL SURFACE

Procedure:
  (a) Pour ~50mL of gel onto a horizontal metal plate or tile.
  (b) Spread to approximately 2mm thickness with spatula.
  (c) Allow to air dry for 30 minutes at room temperature.
  (d) Observe film integrity — does it crack, peel, or remain continuous?

Pass criteria:
  PASS: Film remains continuous (may shrink slightly as water evaporates).
        Does not crack into islands before 30 minutes.
  FAIL: Film cracks severely within 10 minutes (SAP concentration too high
        or bentonite too low)

Record:
  Film integrity at 30 min: ____
  Cracking observed? Y/N, if Y, time to first crack: ____
  PASS / FAIL: ____
```

---

### Step 6: Phase 2 Preparation — Unencapsulated CO₂ Foaming Test

**Purpose:** Before committing to wax microencapsulation, test that the
baking soda + citric acid mixture actually produces useful foaming behavior
when incorporated in the gel. This is the validation step before Step 7.

```
MATERIALS:
  Phase 1 gel (200g, from Steps 1-5)
  Baking soda: 4g (2% of gel weight)
  Citric acid: 3g (1.5% of gel weight)
  Hot water bath (80-100°C) — pot of hot water on stovetop
  Small glass jar or metal container

PROCEDURE:

1. Take 200g of the Phase 1 gel into a small container.

2. Add baking soda and citric acid to the gel and mix well. Keep them
   in the gel but do not add water — they will react slowly in the gel's
   water content. Note any immediate fizzing.
   
   EXPECTED: Some slow background fizzing (acid + bicarb in water).
   This is why encapsulation is needed for the final product — unencapsulated,
   the reaction begins immediately. This step is just to observe the foaming
   behavior.

3. Place the container into the hot water bath (80-100°C).

4. Observe the gel over the next 5-10 minutes as the temperature rises.

5. Record: does the gel expand/foam? How much? What is the consistency
   of the foam?

TARGET BEHAVIOR:
  Gel volume should increase 2-4x as CO₂ is generated.
  Resulting foam should be stable (not collapse immediately).
  Foam should still adhere to the container walls — not purely liquid.

RECORD:
  Foaming observed? Y/N
  Volume increase factor: ____x
  Foam stability at 5 min after heat removal: ____
  Notes: ____

NOTE ON ENCAPSULATION DECISION:
  If foaming is observed and is substantial, proceed to wax encapsulation (Step 7).
  If foaming is minimal, increase baking soda ratio to 4%, retry.
```

---

### Step 7: Phase 2 — Wax Microencapsulation of CO₂ Agent

**Purpose:** Encapsulate the baking soda + citric acid in wax that melts
at ~60-70°C, preventing premature reaction during gel storage and ensuring
the CO₂ generation is triggered by heat at the fire boundary.

**This creates a "latent activation" layer — the gel sits inertly until
fire temperature triggers the wax to melt and the reaction to begin.**

```
MATERIALS:
  Paraffin wax pellets: 50g (melting point ~55-65°C preferred — check product)
  Baking soda: 15g
  Citric acid: 12g
  Double boiler (or microwave + heat-resistant glass bowl)
  Ice bath or cold plate (a metal tray in the freezer)
  Parchment paper
  Gloves (the wax will be hot)

PROCEDURE:

1. Mix baking soda and citric acid powders together in a dry bowl.
   Mix thoroughly. Keep dry — any moisture will trigger premature reaction.

2. Melt paraffin wax in double boiler to ~70°C. Do not overheat.
   Molten wax should be fully liquid and clear/slightly cloudy.

3. Pull the cold plate from the freezer (or place a metal tray on ice).
   Line with parchment paper.

4. Using a teaspoon or dropper, drop a small amount of molten wax
   (~1mL) onto the cold parchment. Allow to partially solidify for 5 seconds
   until it is still soft but not liquid.

5. Quickly place a small amount (~0.3g) of the baking soda/citric acid
   mixture in the center of the soft wax disc.

6. Immediately drizzle another ~0.5mL of molten wax over the top to
   seal the powder inside.

7. Allow to fully solidify (1-2 minutes on cold plate).

8. Repeat to produce 20-30 capsules.

9. Optional refinement: Once solidified, briefly dip each capsule into
   fresh molten wax to ensure complete sealing (second coat).

10. Test a capsule: Drop one capsule into a cup of hot water (~80°C).
    Observe wax melting and CO₂ fizzing.
    
    PASS: Wax melts, capsule releases contents, fizzing observed.
    FAIL: No fizzing (moisture got in during encapsulation — remake dry)

TARGET CAPSULE PROPERTIES:
  Size: ~1-2 cm diameter (lab scale — will be much smaller at production scale)
  Shell: Complete wax coverage, no exposed powder
  Release temperature: 60-75°C (verify with hot water test)

RECORD:
  Capsule count produced: ____
  Release test result: PASS / FAIL
  Release temperature (observed): ____°C
  Notes: ____
```

---

### Step 8: Incorporating Encapsulated Agent into Phase 2 Gel

```
PROCEDURE:

1. Prepare fresh Phase 1 gel (Steps 1-5).

2. Coarsely crush the wax capsules into smaller pieces (5-10mm)
   using your fingers. Keep them as intact as possible — just
   reduce size for better distribution in the gel.

3. Fold the crushed capsule pieces into the gel gently.
   Use a folding motion, not vigorous stirring, to avoid
   mechanical rupture of the capsule shells.

4. The Phase 2 gel is ready for testing once capsules are distributed.

5. Visual check: capsule pieces should be visible throughout the gel.
   No white powder visible (if powder is visible, capsules ruptured —
   see failure modes in Part IV).

RECORD:
  Capsule-to-gel ratio (weight %): ____%
  Visible powder in gel? Y/N (if Y, remake)
  Notes: ____
```

---

## PART III: TESTING PROTOCOL

Four tests. Each is designed to isolate one property of the gel.
Tests are ordered from simplest to most fire-proximate.

---

### Test 1 — Vertical Adhesion and Coverage Test

**What it measures:** Can the gel stay on vertical surfaces without running off?
Does it coat uniformly? This is the deployment precondition test.

```
SETUP:
  Three vertical panels, 30cm × 30cm each:
    Panel A: Plywood (representative combustible surface)
    Panel B: Painted drywall (interior wall surface)
    Panel C: Glass (transparent — lets you observe film uniformity)
    
  Application method: Spatula/brush at 3mm thickness
  Alternative: Fill a spray bottle with gel (if thin enough to spray);
  record nozzle pressure and output rate

PROCEDURE:

1. Mount panels vertically (lean against wall at 90°).

2. Apply gel to each panel from top edge, spread to 3mm thickness
   using a spatula. Try to achieve even coverage.

3. Start timer.

4. Photograph each panel at 0 min, 5 min, 15 min, 30 min.

5. At 30 minutes, measure runoff: weigh any gel that has run off
   into a collection tray at the base of each panel.

6. Calculate retention percentage:
   Retention % = ((initial gel weight - runoff weight) / initial gel weight) × 100

PASS CRITERIA:
  >80% gel retention at 15 minutes: PASS (suitable for pre-flashover window)
  60-80% at 15 minutes: MARGINAL (increase bentonite by 5-10g per batch)
  <60% at 15 minutes: FAIL (see correction protocol in Part IV)

RECORD TABLE:
  Panel    | Applied (g) | Runoff at 5min (g) | Runoff at 15min (g) | Retention% | Result
  ────────────────────────────────────────────────────────────────────────────────────────
  A (Wood) |             |                    |                     |            |
  B (Drywall)|           |                    |                     |            |
  C (Glass)|             |                    |                     |            |
```

---

### Test 2 — Fire Protection Duration (Torch Test)

**What it measures:** How long does the gel coating protect the underlying
wood substrate from ignition? This is the core performance test.

**⚠️ SAFETY:** Conduct outdoors or with excellent ventilation. Have a bucket
of water and fire extinguisher immediately accessible. Wear face shield,
gloves, and non-synthetic clothing. Work on a non-combustible surface.
Never leave a burning test sample unattended.

```
SETUP:
  6 plywood squares, 10cm × 10cm × 12mm thick (standard dimensional lumber)
  Group A: 3 squares — UNCOATED (control)
  Group B: 3 squares — GEL COATED (3mm thickness, applied to top face)
  
  Propane torch (standard hardware store "Bernzomatic" type)
  Stopwatch
  Fireproof surface (concrete block, brick, or metal stand)
  Bucket of water (safety)
  Fire extinguisher nearby (safety)

GEL COATING PROCEDURE:
  Apply gel to the top face of Group B panels to 3mm thickness.
  Allow gel to sit undisturbed on the panel for 5 minutes before testing
  (simulates the gel having had time to adhere and partially set).

TEST PROCEDURE:

1. Place one panel (coated side up) on the fireproof surface.

2. Position propane torch so flame tip is ~25mm (1 inch) from panel surface.
   Light torch.

3. Start stopwatch.

4. Hold torch position fixed. Observe panel continuously.

5. Record the following events with timestamps:
   Event A: First visible steam from gel surface (water evaporation beginning)
   Event B: First visible change in gel appearance (color change, drying)
   Event C: First visible charring or darkening of substrate (if gel fails)
   Event D: First visible ignition (sustained flame from wood surface)
   Event E: Flame self-sustains after torch removal (torch off test:
             at 2 min, briefly remove torch for 5 seconds — does flame sustain?)
   
6. Maximum test duration: 10 minutes per panel.
   If no ignition in 10 minutes, record as ">600 seconds protection."
   Extinguish any sustained flame with water immediately after test.

7. Repeat for all 6 panels (3 control, 3 coated).

8. After cooling, photograph charring pattern on each panel.
   Measure maximum char depth with ruler.

RECORD TABLE:
  Panel ID | Coated? | Time to steam | Time to char | Time to ignition | Char depth | Notes
  ─────────────────────────────────────────────────────────────────────────────────────────
  A1       | No      |      n/a      |              |                  |            |
  A2       | No      |      n/a      |              |                  |            |
  A3       | No      |      n/a      |              |                  |            |
  B1       | Yes     |               |              |                  |            |
  B2       | Yes     |               |              |                  |            |
  B3       | Yes     |               |              |                  |            |

PASS CRITERIA:
  STRONG PASS: Coated panels show >3x longer time-to-ignition vs. control
  PASS: Coated panels show >2x longer time-to-ignition vs. control
  MARGINAL: Coated panels show 1.5-2x improvement
  FAIL: Coated panels show <1.5x improvement vs. control
```

---

### Test 3 — Aerogel Transition Observation

**What it measures:** When the water evaporates from the gel under heat,
does the silica component form a visible residual solid barrier (the aerogel
transition)? This is the Phase C protection validation.

**This is an observational test, not a fire test.**

```
SETUP:
  3 metal plates (cookie sheet sections or similar), 10cm × 10cm
  Phase 1 gel: 50g per plate, spread to 3mm
  Oven or heat gun
  Infrared thermometer (if available; otherwise use oven dial setting)
  Camera

PROCEDURE:

1. Apply 50g of gel to each metal plate, spread to 3mm thickness.
   Weigh each plate: total weight = plate + gel weight.

2. Place plates in oven at 200°C (or use heat gun directed at plate).
   OR: use the propane torch on LOW, held at 15-20cm distance.
   The goal is steady heat, not direct impingement flame.

3. Observe and photograph at 2 min intervals.

4. Record: how long until gel appears dry/cracked?
   What does the residue look like?

5. After 20 minutes or when gel appears completely dried:
   Remove from heat, allow to cool fully.

6. Weigh residue (plate + dried residue).
   Calculate: residue mass = (final total weight) - (plate weight alone)

7. Attempt to brush the residue off the plate with a dry finger.
   
   AEROGEL FORMATION INDICATOR:
   If the residue is a WHITE or OFF-WHITE fluffy/powdery coherent layer
   that resists brushing off and holds together as a porous structure
   → Silica aerogel transition has occurred. PASS.
   
   FAILURE INDICATOR:
   If the residue is a flat, dense, dark layer, or crumbles to
   nothing, or the plate is bare — transition did not occur or
   silica concentration is insufficient.

RECORD TABLE:
  Plate | Initial gel mass | Time to dry | Residue mass | Residue description | Aerogel? Y/N
  ────────────────────────────────────────────────────────────────────────────────────────────
  1     |                  |             |              |                     |
  2     |                  |             |              |                     |
  3     |                  |             |              |                     |

PASS CRITERIA:
  PASS: White/off-white coherent residue observed, resists brushing = aerogel formation
  MARGINAL: Partial residue present but fragile/patchy = increase silica concentration
  FAIL: No residue or bare plate after drying = silica not transitioning (see Part IV)
```

---

### Test 4 — Phase 2 Foaming Activation Test (If Step 7 completed)

**What it measures:** Do the wax capsules remain inert at room temperature
but activate under heat, producing the CO₂ foam layer?

```
PROCEDURE:

Phase 4a — Cold stability test (room temperature, 24 hours):
  1. Apply Phase 2 gel (with capsules) to a vertical panel.
  2. Leave undisturbed at room temperature for 24 hours.
  3. Observe: is there any visible CO₂ foaming, bubbling, or gel disruption?
  PASS: No foaming observed. Capsules remain intact at ambient temperature.
  FAIL: Foaming visible (capsules ruptured prematurely — see Part IV)

Phase 4b — Heat activation test:
  1. Apply Phase 2 gel (with capsules) to a horizontal metal plate.
  2. Apply heat gun or low torch at ~15cm until plate surface reaches
     ~70-90°C (estimated — gel will start steaming around this range).
  3. Observe: does the gel foam and increase in volume as it heats?
  PASS: Visible foam formation, volume increase >1.5x, foam holds structure.
  FAIL: No foaming (capsules not releasing — check encapsulation integrity)

RECORD:
  Phase 4a result: PASS / FAIL
  Phase 4b result: PASS / FAIL
  Foam expansion factor: ____x
  Foam stability (holds structure at 5 min after heat removal): Y/N
  Notes: ____
```

---

### Test 5 — Target Use Case Simulations

**What it measures:** Does the gel actually perform in the specific use case
you are designing for? This is the translation from lab properties to
real-world application.

**Choose the test most relevant to your primary use case:**

```
USE CASE A: Room surface pre-flashover suppression
  Setup: Apply gel to all surfaces of a small model enclosure
  (a cardboard box is NOT adequate — use a metal or gypsum board box,
  ~30cm cube). Place a small candle inside. Apply gel to all interior
  surfaces before igniting candle.
  Test: Does the candle's small flame spread to the gel-coated surfaces?
  Measure: How long before any visible char on coated surfaces?
  Compare to: Identical enclosure without gel coating.

USE CASE B: Handheld extinguisher after-treatment
  Setup: Ignite a small piece of kindling wood on a fireproof surface.
  Allow to burn for 30 seconds (small flame established).
  Extinguish with water (simulating conventional extinguisher knockdown).
  Immediately apply gel to the charred surface (simulating gel extinguisher use).
  Test: Does the previously ignited surface re-ignite after gel application
  when a match is held to it?
  Compare to: Same protocol without gel application — measure re-ignition time.

USE CASE C: Surface coating for structural element protection
  Setup: Pre-apply gel to a wood panel, then apply direct flame for
  a standardized 60-second burn (from Test 2).
  Compare panel to uncoated control.
  Measure char depth after 60 seconds. Is it meaningfully less?

RECORD all use case tests with:
  Test date, formulation batch number, procedure used,
  pass/fail criteria, observations, photographs.
```

---

## PART IV: PREDICTABLE FAILURE MODES AND CORRECTION PROTOCOL

This section states in advance what will most likely go wrong, why it goes
wrong, and exactly how to fix it. Most first-batch failures are predictable
from the material science.

---

### Failure Mode 1 — Gel Runs Off Vertical Surface (Too Thin)

```
SYMPTOMS:
  Gel slides down surface within 5 minutes.
  Retention < 60% at 15 minutes.
  Looks and flows like a liquid when applied.

ROOT CAUSE:
  Insufficient yield stress — the gel has no structural network strong enough
  to resist gravity on a vertical surface.
  Most common cause: Too little bentonite (clay) or too little SAP.

CORRECTION PROTOCOL:

  Step 1: Add bentonite in 10g increments to the batch.
  Mix thoroughly between additions. Retest adhesion after each addition.
  Maximum add: 50g additional bentonite (beyond baseline 30g) before
  considering a different approach.
  
  Step 2 (if Step 1 insufficient): Increase SAP by 2g increments.
  SAP increases both viscosity and elastic network strength.
  Retest after each addition.
  
  Step 3 (if Steps 1-2 insufficient): Add 5g of xanthan gum (available
  in grocery store baking section) per liter of gel. Xanthan gum is a
  powerful shear-thinning thickener. Dissolve in 50mL warm water before
  adding to gel.
  
  Step 4 (reformulate): Reduce water in Step 1. Use 800mL instead of 1000mL
  for the SAP hydration step. This creates a denser base gel.

RECORD: Which correction applied, new adhesion test result.
```

---

### Failure Mode 2 — Gel Too Thick to Spray or Pour (Too Viscous)

```
SYMPTOMS:
  Gel will not flow from a container.
  Cannot be sprayed through a nozzle.
  Requires significant force to stir.
  Forms a rubbery solid rather than a viscous gel.

ROOT CAUSE:
  Too much SAP or too much bentonite.
  The gel has crosslinked into a solid-like network beyond the target
  shear-thinning region.

CORRECTION PROTOCOL:

  Step 1: Add distilled water in 100mL increments.
  Mix thoroughly between additions.
  Test: Does the gel flow when poured? If yes, stop.
  
  Step 2: Confirm shear-thinning is maintained after dilution.
  (Run Check 2 from Step 5 after dilution.)
  
  Step 3 (if dilution doesn't restore flow): The SAP concentration is
  too high. Start fresh batch with reduced SAP (7g instead of 10g per liter).
  
  NOTE: There is a formulation sweet spot. Too thick = can't deploy.
  Too thin = can't hold. Iterating between these two failure modes is the
  primary Phase 1 calibration task.

RECORD: Water added, new flow assessment.
```

---

### Failure Mode 3 — No Aerogel Transition (Bare Plate After Drying)

```
SYMPTOMS:
  After Test 3, the dried residue is negligible or absent.
  The metal plate is bare or has only a faint white powder that brushes off.
  No coherent white layer visible.

ROOT CAUSE:
  Insufficient fumed silica concentration, or
  silica not uniformly dispersed in the gel (clumps that remain on the
  plate as discrete particles rather than an interconnected network).

CORRECTION PROTOCOL:

  Step 1: Double the fumed silica concentration.
  Use 50g silica in 200mL water (instead of 25g) for the dispersion step.
  Redo Test 3 with this higher-silica formulation.
  
  Step 2 (if Step 1 insufficient): The silica must be more finely dispersed.
  Extend the mixing step (Step 2 in synthesis) to 10+ minutes.
  If an immersion blender is available, use it for 5 minutes on the silica
  dispersion before incorporating into the gel.
  
  Step 3: Check the product specification of your fumed silica.
  If BET surface area is not specified, you may have a coarser grade.
  Aerosil 200 or Cab-O-Sil M5 are the target grades (BET ~200 m²/g).
  Lower surface area grades will not form aerogel structure effectively.
  Order a confirmed Aerosil 200-grade product if unknown grade was used.
  
  Step 4: Increase silica to 10% of gel dry weight (from 5%).
  Test 3 should now show visible white aerogel layer.

RECORD: Silica concentration in failed batch, correction applied, new Test 3 result.
```

---

### Failure Mode 4 — Gel Ignites or Burns (Combustion of Gel)

```
SYMPTOMS:
  During Test 2, the gel itself catches fire.
  Flame spreads along the gel surface rather than being suppressed.

ROOT CAUSE (most likely):
  SAP polymer component is burning. This can occur if:
  (a) SAP concentration is very high and water has evaporated from the surface
      before the test, leaving a dry polymer layer (which IS combustible)
  (b) The SAP product used contains combustible additives (some consumer
      "water crystals" contain fragrance oils or other organics)

CORRECTION PROTOCOL:

  Step 1: Immediately confirm this is not a product-quality issue.
  Test the SAP powder alone: hold a small pinch in metal tongs and expose
  to a lighter flame. Pure sodium polyacrylate will NOT ignite easily —
  it chars but does not sustain flame. If the powder ignites easily,
  the product has combustible additives. Source a purer grade.
  
  Step 2: Ensure gel is freshly applied for Test 2 (not dried out).
  Never test gel that has been sitting on a surface long enough to lose
  significant water content. Reapply fresh gel immediately before torch test.
  
  Step 3: If the gel is being tested after partial drying, increase the
  silica and clay concentrations — these non-combustible components will
  dominate the dried residue and prevent ignition of the polymer char.
  
  Step 4 (worst case): Switch from SAP base to cellulose-based thickener.
  HPMC (hydroxypropyl methylcellulose, available from brewing supply or
  pharmacy as "Methocel") is less combustible than SAP and USFS-approved.
  Substitute 1% HPMC in water for the SAP hydration step.

RECORD: Source of SAP used, whether pure SAP test was performed, correction applied.
```

---

### Failure Mode 5 — Wax Capsules Rupture Prematurely in Gel

```
SYMPTOMS:
  Phase 4a cold stability test fails — foaming observed at room temperature.
  Gel has visible bubbles when stored.
  Smell of CO₂ fizzing from gel container.

ROOT CAUSE:
  Wax capsule shells are too thin, were incompletely sealed, or the wax
  melting point is too low (capsules soften at room temperature).

CORRECTION PROTOCOL:

  Step 1: Remake capsules with a SECOND WAX COAT.
  After initial solidification, dip each capsule into fresh molten wax
  twice more to build up shell thickness.
  
  Step 2: Check the melting point of your paraffin wax.
  Target is 55-65°C. If your wax melts at <50°C, it may soften at warm
  ambient temperatures (or from hands during handling). Source a higher-
  melting-point paraffin (available from candle supply — specify 60-65°C
  melting point).
  
  Step 3: Keep capsules refrigerated until incorporation into gel.
  Cold storage prevents any thermally-induced shell softening.
  
  Step 4: Add capsules to gel only immediately before use if the application
  is same-day. Do not store Phase 2 gel with capsules for more than 24 hours.

RECORD: Capsule failure timeline, wax melt point used, correction applied.
```

---

### Failure Mode 6 — Fire Protection Duration Not Significantly Better Than Control

```
SYMPTOMS:
  Test 2 results: coated panels ignite at similar time to uncoated control.
  <1.5x improvement in time-to-ignition.

ROOT CAUSE ANALYSIS (work through in order):

  ROOT CAUSE A: Gel film thickness insufficient
  3mm may be too thin for the heat flux from the propane torch.
  The torch heat flux (~50-100 kW/m²) is very high — comparable to
  post-flashover conditions, which exceed this system's primary scope.
  (Note: This is a test design issue as much as a gel performance issue.)
  
  CORRECTION A: Increase gel film to 6-8mm. Retest.
  Also test at a LONGER torch distance (50mm instead of 25mm) to simulate
  lower heat flux (nearer to pre-flashover conditions rather than direct flame).

  ROOT CAUSE B: Aerogel transition not occurring (see Failure Mode 3)
  If the silica is insufficient, there is no Phase C protection.
  After water evaporates (which absorbs heat for a time), the substrate
  is exposed with no persistent barrier.
  
  CORRECTION B: Resolve Failure Mode 3 first, then retest.

  ROOT CAUSE C: Gel dried too quickly on the panel before the test
  If the gel was applied long before the test, it may have dried.
  
  CORRECTION C: Apply gel fresh (<2 minutes before torch application).
  Or add 5% glycerol (food-grade, pharmacy) to the gel — glycerol is a
  humectant that slows water evaporation without combusting.

  ROOT CAUSE D: Propane torch is testing the wrong domain
  The propane torch represents an extremely high heat flux — far beyond
  what a growing residential fire produces in the pre-flashover window.
  The system is designed for pre-flashover surface propagation suppression,
  not direct impingement flame resistance.
  
  CORRECTION D: Reframe the test. Test the gel against a small candle flame
  (much lower heat flux, representative of early fire growth). Or test
  flame SPREAD along a surface rather than direct impingement resistance.
  Use two panels in a V-configuration: ignite the bottom edge of the
  angled assembly and observe if flame travels up coated vs. uncoated panel.
  This tests surface propagation arm suppression directly.

RECORD: Root cause identified, correction applied, retest results.
```

---

### Failure Mode 7 — Fumed Silica Clumping — Visible White Lumps in Gel

```
SYMPTOMS:
  White lumpy or grainy texture in the finished gel.
  Silica visible as discrete particles rather than uniform dispersion.

ROOT CAUSE:
  Fumed silica was added too fast, or insufficient mixing energy.

CORRECTION PROTOCOL:

  Step 1: Use an immersion blender instead of a whisk for the silica
  dispersion step. High-shear mixing is essential for fumed silica.
  
  Step 2: Add silica in smaller portions (5-6 additions rather than 4).
  More passes through the mixing cycle = better dispersion.
  
  Step 3: Add 1-2 drops of dish soap to the water before adding silica.
  The surfactant dramatically improves wetting of the hydrophobic silica
  surface. This small amount will not compromise the gel.
  
  Step 4 (most effective): Prepare silica dispersion the day before.
  Allow the wetted silica to sit overnight in the water before mixing.
  The prolonged contact time allows more complete surface wetting.
  Then blend vigorously the following day before incorporation.

RECORD: Mixing method used, observation of dispersion quality.
```

---

## PART V: ITERATION ROADMAP — MOVING FROM PASS TO OPTIMIZED

Once Phase 1 tests show a PASS on all four tests, the gel is a validated
baseline. Optimization proceeds through systematic single-variable changes.

```
ITERATION PHASE 1 → PHASE 2 (Establish performance baselines)

Variable 1: Silica concentration
  Test: 3%, 5%, 8%, 12% silica (by dry weight of gel)
  Measure: Test 3 aerogel quality score + Test 2 protection duration
  Find: Minimum silica for consistent aerogel formation

Variable 2: Bentonite concentration
  Test: 2%, 4%, 6%, 8% bentonite
  Measure: Test 1 vertical retention %
  Find: Minimum bentonite for >80% retention at 15 min

Variable 3: SAP concentration
  Test: 0.5%, 1%, 1.5%, 2% SAP
  Measure: Test 1 retention + deployability (shear-thinning quality)
  Find: Sweet spot between deployability and adhesion

Variable 4: Gel film thickness
  Test: 1mm, 2mm, 3mm, 5mm
  Measure: Test 2 protection duration
  Find: Minimum useful film thickness (cost implications for at-scale)

ITERATION PHASE 2 → PHASE 3 (Enhance for specific use case)

Variable 5: CO₂ foaming agent loading
  Test: 1%, 2%, 4% encapsulated agent
  Measure: Test 4 foam expansion factor + Test 2 protection duration
  Find: Whether foam enhancement meaningfully extends protection duration

Variable 6: Glycerol addition (humectant)
  Test: 0%, 3%, 5%, 8% glycerol
  Measure: Test 2 protection duration + gel shelf life at ambient temperature
  Find: Whether glycerol extends protection window by slowing water loss

SYSTEMATIC RECORD FORMAT FOR ALL ITERATIONS:
  Batch ID:
  Date:
  Formulation (all components and concentrations):
  Test 1 result (retention %):
  Test 2 result (protection duration, coated vs. control ratio):
  Test 3 result (aerogel observed Y/N, quality):
  Test 4 result (if applicable):
  Notes (anomalies, observations):
  Next iteration planned:
```

---

## PART VI: MINIMUM VIABLE PROTOTYPE ACCEPTANCE CRITERIA

A prototype batch of EHAB is considered validated for target use case
exploration when it meets ALL of the following:

```
CRITERION 1 — Vertical Adhesion:
  >80% retention on a painted wood or drywall vertical surface at 15 minutes

CRITERION 2 — Fire Protection Duration:
  >2x improvement in time-to-ignition vs. uncoated control
  (at 25mm torch distance, 3mm gel film, fresh application)

CRITERION 3 — Aerogel Transition:
  Visible white/off-white coherent residue after 200°C oven drying test,
  residue resists brushing with dry finger

CRITERION 4 — Non-Ignition of Gel:
  Gel applied fresh to a surface does not sustain independent flame
  when torch is removed

CRITERION 5 — Deployability:
  Gel flows through a 5mm nozzle opening under hand-squeeze pressure
  (simulating basic sprayer/extinguisher deployment)

WHEN ALL 5 CRITERIA ARE MET:
  The formulation is recorded as the MVP formulation.
  Component concentrations and mixing procedure are locked.
  Further iterations are ENHANCEMENT only, not correction.
  The MVP formulation is the basis for all use-case application testing.
```

---

## DOCUMENT METADATA

```
Document ID:    EHAB_PROTOTYPE_PROTOCOL_v1.0
Date:           2026-05-06
Framework:      Causal Geometry / Triadic Invariant (OrganismCore)
Author:         Eric Robert Lawson / Copilot Synthesis Agent

TOTAL PHASE 1 ESTIMATED COST:     $116-188 (ingredients + basic equipment)
TOTAL INGREDIENT COST ONLY:       $53-81
ESTIMATED BATCHES POSSIBLE:       20+ test batches from Phase 1 inventory
BATCH SIZE (Phase 1):             ~1.75 liters per batch
ESTIMATED TIME TO MVP:            2-4 weeks of weekend prototyping sessions

KEY SOURCES FOR FORMULATION TARGETS:
  Stanford / Dong, Appel et al.: Cellulose-SAP + colloidal silica
  hydrogel-to-aerogel transition (Advanced Materials 2024; arXiv 2503.14923)
  MDPI Gels / Scientific Reports: Multi-synergistic flame-retardant hydrogel
  coatings for structural applications (2023-2024)
  UL 94 / ASTM D3801: Small-scale flammability test protocol basis
  Sigma-Aldrich / Alibaba / bulk chemical pricing: Phase 1 cost estimates

SAFETY SUMMARY:
  MANDATORY: N95 respirator for fumed silica handling — no exceptions.
  Nitrile gloves for all synthesis steps.
  Face shield for fire tests.
  Fire tests outdoors or with active ventilation only.
  Water bucket + fire extinguisher present during all fire tests.
  No synthetic clothing during fire tests (cotton/wool only).
  Never test indoors without active ventilation.
  Gel components (SAP, silica, bentonite, baking soda, citric acid,
  paraffin) are non-toxic at handling quantities. The fire tests
  are the only meaningful hazard in the protocol.

RELATED DOCUMENTS:
  encapsulating_gel_suppression_causal_geometry_derivation.md
  gel_substance_deployment_gradient_analysis.md
  Onboarding_3/causal_geometry_onboarding.md
  Onboarding_3/Coherence_Is_Not_Truth.md

NEXT DOCUMENT:
  EHAB_MVP_FORMULATION_RECORD.md
  (to be created after MVP criteria are met — locks the formulation
  and transitions from prototyping to use-case application testing)
```

---

*The protocol follows the geometry to where it can be physically tested.*
*The formulation is derivable from components on a shelf.*
*The tests are falsifiable against the predictions.*
*The failure modes are predictable before the first batch is made.*
*Coherence is not truth — the bench confirms or corrects.*
