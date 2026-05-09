# P1S Synthesis — Timing, Chronology, and Optimal Mixing Windows
## EHAB Phase 1 Standard Protocol — Time-State Analysis
## OrganismCore Technical Record
## Date: 2026-05-09

---

## PURPOSE OF THIS DOCUMENT

The standard Phase 1 protocol describes *what* to do and *why* each step matters
structurally. What it does not map explicitly is the time-state of each component —
the window during which each prepared component is in its optimal state for
compounding, and what happens to it if it sits too long before being used.

This document analyzes those time-states from first principles, maps the chronology
of all three components against each other, and produces an optimized parallel
preparation schedule that minimizes the total session time while ensuring every
component is compounded at its quality peak.

The conclusion is that the standard protocol's *sequential* approach is the primary
source of preventable quality loss. Two of the three components (SAP and bentonite)
have rest periods that can run in parallel. The silica suspension has a specific
degradation window that is the most time-sensitive constraint in the system.

---

---

# PART I: TIME-STATE ANALYSIS OF EACH COMPONENT

## 1.1 — SAP Hydrogel (A1): The Overcrosslinking Problem

### What is happening during the SAP rest period

When SAP powder contacts water, individual polymer chains absorb water
through osmotic uptake and begin to swell. The crosslinked network
expands, forming a gel. This is not instantaneous — individual granules
absorb water from the outside in, forming a gel shell around a still-dry
interior. The 30-minute rest period is the time needed for the interior
of each granule to fully hydrate and for the entire polymer network to
reach thermodynamic equilibrium.

### The optimal state window

```
SAP TIME-STATE PROFILE:

  0–12 min:   ACTIVE HYDRATION. Water being added incrementally.
              Not yet usable — polymer not fully swollen.

  12–30 min:  NETWORK FORMATION REST. Polymer chains interconnecting.
              Do not disturb. Do not compound yet.

  30–60 min:  ★ OPTIMAL STATE ★
              Network fully formed. Maximum water retention.
              Viscosity at peak. Surface silanol groups on SAP
              carboxylate chains are maximally available for
              hydrogen bonding with silica particles during compounding.
              This is the compounding window.

  60–120 min: ACCEPTABLE — slight over-relaxation begins.
              The network continues to equilibrate past its kinetic
              peak. Practically, the gel is still usable, but
              viscosity may have increased slightly (continued
              chain entanglement) or decreased slightly (chain
              relaxation) depending on granule size.

  >120 min:   DEGRADATION BEGINS — syneresis risk.
              In undisturbed SAP hydrogel at room temperature, the
              polymer network can begin to expel water (syneresis)
              as the crosslinked chains contract. Visually: a thin
              layer of clear water appears on top of the gel.
              The gel beneath becomes denser. This is not catastrophic
              immediately, but uniformity is compromised.
              If syneresis is observed: stir thoroughly before
              compounding. Accept that silica distribution will
              be slightly less uniform.

  >4 hours:   SIGNIFICANT DEGRADATION.
              Syneresis pronounced. Network cohesion reduced.
              Do not use for precision applications — acceptable
              only for fire suppression tests where uniformity
              is less critical.
```

**Key constraint:** The SAP gel should be compounded within 30–90 minutes
of the rest period beginning. 60 minutes is the sweet spot. Do not let it
sit unattended for more than 2 hours before compounding.

**What "sitting too long" looks like:** Clear water accumulating on the surface.
The gel pulling away slightly from the bowl edges. The top surface appearing
drier or more matte than expected.

**Recovery:** Stir vigorously for 60 seconds before compounding. The gel will
re-homogenize but will not fully recover the peak network cohesion. Record this
if it occurs.

---

## 1.2 — Fumed Silica Suspension (A2): The Settling and Re-Agglomeration Problem

### What is happening in the silica suspension

Fumed silica particles are 7–40 nm in diameter. At this scale, they are
technically colloidal — they do not sediment quickly by gravity in the way
that larger particles do. Instead, the instability mechanism is different:
**particle re-agglomeration**.

When fumed silica particles are dispersed in water by vigorous mixing,
the particle-particle silanol bonds (Si-OH...HO-Si) that hold the original
powder aggregates together are broken. The particles are individually
wetted and carry a negative surface charge (from deprotonated silanol groups
at neutral pH), which creates electrostatic repulsion that keeps them apart.

Over time at rest, two things happen:
1. The electrostatic repulsion weakens slightly as particles drift closer.
2. Hydrogen bonding between silanol groups on adjacent particles begins
   to re-form, creating new inter-particle bridges.

This re-agglomeration process produces a silica network — the suspension
begins to self-gel. This is the same chemistry that produces the aerogel
layer under heat in the final P1S gel, but here it is happening prematurely
in the suspension, before compounding.

### The critical issue with the P1S silica concentration

In the P1S formulation, silica is 40g in approximately 350mL water — a
concentration of roughly 10.3% by mass. At this concentration, the silica
suspension is above the threshold at which colloidal silica suspensions
begin to exhibit significant gel network formation over time. This is
the concentration regime where self-gelation is a real and measurable concern.

```
SILICA SUSPENSION TIME-STATE PROFILE:

  0–15 min:   ACTIVE DISPERSION. Silica being added and whisked.
              Not yet usable — dispersion not complete.

  15–25 min:  REST + INSPECTION. The 10-minute rest period in the protocol.
              Silica particles fully wetted. Suspension at maximum
              uniformity. Electrostatic repulsion at maximum.
              This is when the suspension is most homogeneous.

  25–45 min:  ★ OPTIMAL COMPOUNDING WINDOW ★
              Particles fully dispersed. Re-agglomeration not yet
              significant. Suspension flows freely when poured.
              This is the target window for compounding silica
              into the SAP+bentonite matrix.

  45–90 min:  EARLY RE-AGGLOMERATION.
              Particle clusters beginning to reform. Suspension
              viscosity increasing measurably. The milky appearance
              may become slightly more opaque at the bottom of the
              container as small aggregates form.
              Still usable but requires vigorous re-stirring
              immediately before addition to SAP+bentonite.
              Re-stir for 2 full minutes before use.

  >90 min:    SIGNIFICANT RE-AGGLOMERATION.
              Visible white sediment or thickened zones in the
              suspension. A surface film may appear (floating
              silica clusters). The suspension no longer disperses
              as uniformly into the SAP matrix — pre-formed
              aggregates resist breakup during compounding.
              This is the primary cause of white lumps and patchy
              aerogel formation in the final gel.
              
              Recovery: Use immersion blender at high speed for
              3 minutes immediately before compounding.
              This partially breaks re-agglomerates but cannot
              fully restore the original dispersion quality.

  >3 hours:   PARTIAL GELATION.
              The silica suspension may begin to hold its shape
              (cohesive network forming). This is pre-aerogel
              chemistry occurring at room temperature. The suspension
              is significantly compromised. If this occurs, add
              50mL additional distilled water, blend at high speed
              for 3 minutes, use immediately.
```

**Key constraint:** The silica suspension has the tightest quality window
of the three components. It should be compounded within 45 minutes of
completing the dispersion step. It is the time-critical bottleneck of the
entire synthesis.

**What "sitting too long" looks like:** The bottom of the container appears
whiter or thicker than the top. The suspension pours with slight resistance.
A thin white film on the surface that re-incorporates when stirred.

---

## 1.3 — Bentonite Dispersion (A3): The Over-Gelation Problem

### What is happening during the bentonite rest period

Sodium bentonite clay platelets carry a net negative charge on their faces
and a net positive charge on their edges. In water, these platelets
hydrate — water intercalates between the aluminum silicate layers, causing
the clay to swell to up to 15× its dry volume. As hydration proceeds,
the charged platelet surfaces repel each other while the edge-face charge
attractions create a loose "house of cards" network — this is the source of
thixotropy.

The 45-minute rest period is the minimum needed for the interlayer water
to fully intercalate and for the thixotropic network to establish. Longer
rest improves the network.

### The over-gelation ceiling

Unlike SAP (which begins to degrade past ~2 hours) and silica (which
re-agglomerates past ~45 minutes), bentonite is the most patient of the
three components. It improves with rest up to approximately 8–12 hours.
However, there is a ceiling beyond which more rest produces a stiffer,
more elastic gel that is harder to compound uniformly into the SAP matrix.

```
BENTONITE TIME-STATE PROFILE:

  0–10 min:   ACTIVE DISPERSION. Water being added. Platelets
              beginning to separate and hydrate.

  10–45 min:  HYDRATION REST. Interlayer water intercalating.
              Thixotropic network forming. Do not disturb.

  45–60 min:  ★ MINIMUM ACCEPTABLE — GOOD FOR MOST USES ★
              Thixotropy present. Edge-face network established.
              Compounding now is acceptable.

  60–90 min:  ★ OPTIMAL STATE ★
              Thixotropic network fully developed. Maximum gel
              strength. The bentonite at this stage integrates
              most smoothly with the SAP network during compounding
              — the developed platelet network distributes
              throughout the SAP matrix rather than forming
              isolated dense zones.

  90 min–4 hours: EXTENDED REST — ACCEPTABLE, SLIGHTLY STIFFER.
              The network continues developing. The gel becomes
              incrementally firmer. Compounding may require more
              vigorous stirring to fully integrate with SAP.
              Performance in the final gel is equivalent or slightly
              better than the 60-minute state.

  >6 hours:   OVER-GELLED — requires careful handling.
              The bentonite network has developed to the point
              where it resists integration. When added to SAP,
              it may not fully disperse — remaining as visible
              grey-tan streaks through the gel. If this occurs:
              stir the bentonite vigorously for 90 seconds before
              adding to SAP, and whisk the combined SAP+bentonite
              for 5 minutes (vs. the standard 3 minutes) to ensure
              full integration.
```

**Key constraint:** Bentonite is the most tolerant of the three components
to extended rest. A bentonite dispersion prepared 60–90 minutes before
compounding is ideal. Anything from 45 minutes to 4 hours is usable.
Beyond 6 hours, expect to compensate with extra mixing energy at compounding.

---

---

# PART II: THE PROBLEM WITH SEQUENTIAL PREPARATION

## 2.1 — What the Standard Protocol Sequential Order Produces

The standard protocol prepares components in this sequence:

```
STANDARD SEQUENTIAL PROTOCOL (as written):

  T+00:00  Begin A1 — SAP hydration (add water incrementally)
  T+00:12  Finish adding water to SAP
  T+00:12  SAP rest begins
  T+00:42  SAP rest complete (30 min) ★ SAP at optimal state ★
  T+00:42  Begin A2 — Silica dispersion
  T+00:57  Finish adding silica to water + 10-min rest begins
  T+01:07  Silica rest complete ★ Silica at optimal state ★
  T+01:07  Begin A3 — Bentonite dispersion
  T+01:17  Finish adding bentonite to water
  T+01:17  Bentonite rest begins (45 min minimum)
  T+02:02  Bentonite rest complete (45 min) ★ Bentonite at acceptable state ★
  T+02:02  Begin A4 — Compounding (bentonite into SAP first)
```

**The problem this creates:**

At the point compounding begins (T+02:02):
- SAP has been sitting for **80 minutes** since water addition completed.
  It is past its optimal 30–60 minute window.
  Syneresis risk is beginning. Network is over-relaxed.
- Silica suspension has been sitting for **55 minutes** since dispersion completed.
  It is past its 25–45 minute optimal window.
  Early re-agglomeration is occurring.
- Bentonite is at its minimum acceptable state (45 minutes).
  This is the only component that is exactly on time.

The sequential protocol inadvertently puts the two most time-sensitive
components (SAP and silica) into suboptimal states before compounding,
while the most time-tolerant component (bentonite) is compounded immediately
at its minimum — with no buffer built in.

**The total session time is also inflated:** 2 hours 2 minutes of preparation
before compounding even begins, with the subsequent 25-minute compounding
and 30-minute final rest taking total time to ~3 hours.

---

# PART III: THE OPTIMIZED PARALLEL SCHEDULE

## 3.1 — The Core Insight

The three rest periods are independent. Nothing prevents preparing bentonite
and SAP simultaneously. The only hard constraint is:

1. SAP must be at 30+ minutes of rest before compounding.
2. Bentonite must be at 45+ minutes of rest before compounding.
3. Silica must be compounded within 45 minutes of dispersion completing.
4. Silica is added LAST (after SAP+bentonite are combined) — this sequencing
   is non-negotiable for network structure reasons.
5. The silica dispersion step (A2) requires the N95 respirator and is
   the only step that cannot be easily paused.

The optimization: **start bentonite first, start SAP second, do silica last
immediately before compounding**. This brings all three components to
their optimal states simultaneously and eliminates the drift periods
that the sequential protocol creates.

## 3.2 — The Optimized Parallel Protocol

```
OPTIMIZED PARALLEL SCHEDULE:

SESSION START — PRE-WORK (5 min):
  Set out all equipment. Tare scale. Put on gloves.
  Measure and pre-set distilled water volumes:
    - 400mL in a container labeled BENTONITE
    - 1000mL in a container labeled SAP
    - 350mL in a container labeled SILICA (set aside for later)
  Do NOT open silica container yet. No N95 needed yet.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
T+00:00  ► START A3 — BENTONITE DISPERSION
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Weigh 30g bentonite into medium bowl.
  Add 400mL water in two 200mL additions with 90 sec whisking each.
  Whisk continuously for 3 minutes.
  Cover bowl with a plate or cling wrap.
  SET TIMER: 60 minutes (targets the optimal 60-90 min rest range).

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
T+00:13  ► START A1 — SAP HYDRATION (immediately after bentonite)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Weigh 10g SAP into large bowl.
  Add 1000mL water in 10 × 100mL increments.
    60 seconds vigorous stirring after each increment.
    Total addition time: ~12 minutes.
  After final addition, set bowl aside undisturbed.
  SET TIMER: 30 minutes for SAP rest.
  Do NOT stir during rest.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
T+00:25  ► CHECKPOINT — Both rests are now running simultaneously
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Current state:
    Bentonite: 25 min into 60-min rest → 35 min remaining
    SAP:       12 min into 30-min rest → 18 min remaining

  Use this time:
    - Set up the large compounding bowl
    - Set up recording materials (notebook, phone)
    - Set out the silica container (do NOT open)
    - Have the 350mL water for silica ready
    - Confirm N95 respirator is accessible
    - Rest / observe

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
T+00:55  SAP REST COMPLETE (42 min actual rest — within optimal window)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Check SAP:
    Expected: Clear to milky gel, holds a slow peak.
    Troubleshoot if needed (see protocol).
  Do NOT compound yet — bentonite needs 18 more minutes.
  SAP can hold in optimal state for another 30–45 min from here.
  Leave it undisturbed.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
T+01:00  ► START A2 — FUMED SILICA DISPERSION
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  ⚠ N95 RESPIRATOR ON BEFORE THIS STEP. Gloves on.

  Weigh 40g fumed silica.
  Seal silica container immediately after weighing.
  Disperse 40g silica into 350mL water in 8 × 5g portions.
    90 seconds vigorous whisking per portion.
    Final 4-minute continuous whisk (or 3-min immersion blender).
  Rest 10 minutes. Inspect for floating powder.
  Re-whisk 2 minutes if floating powder observed.
  N95 may be removed after suspension is confirmed uniform.

  SET TIMER: 15 minutes (this is the start of the silica quality clock).
  The silica suspension is now in its optimal window.
  It must be compounded before the 45-minute mark from NOW.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
T+01:13  BENTONITE REST COMPLETE (60 min actual rest — OPTIMAL STATE)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Check bentonite thixotropy:
    Stir vigorously 30 seconds. Observe recovery.
    PASS: surface stills and firms within 2 minutes.
    FAIL: remains liquid → add 10g bentonite, wait 30 min.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
T+01:13 to T+01:20  (after bentonite check clears)

  Confirm silica time status:
    Silica dispersion was completed at approximately T+01:15.
    At T+01:20, silica has been resting ~5 minutes. ✓ Still optimal.
    Compounding must begin by approximately T+02:00 (45-min silica window).

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
T+01:20  ► BEGIN A4 — COMPOUNDING
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  COMPONENT STATES AT COMPOUNDING:
    SAP:        ~67 min total rest — ★ OPTIMAL WINDOW ★
    Bentonite:  ~60 min rest       — ★ OPTIMAL STATE ★
    Silica:     ~5-10 min rest     — ★ OPTIMAL WINDOW ★

  All three components simultaneously at their best state.
  This is the target the optimized schedule is designed to achieve.

  STEP 1 — BENTONITE INTO SAP:
    Pour bentonite in thin stream into SAP while stirring.
    Whisk vigorously 3 minutes.
    Mixture should be uniform grey-tan.

  STEP 2 — SILICA INTO SAP+BENTONITE:
    Re-stir silica suspension vigorously for 30 seconds before use.
    Divide into 6 equal portions.
    Add each portion in a thin stream around the bowl.
    Whisk 2 minutes between portions.
    Final extended mix: 6 min by hand or 3 min immersion blender.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
T+01:55  COMPOUNDING COMPLETE. SET FINAL REST.
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Cover bowl. Do not stir. 30-minute final rest.
  During this rest: set up validation test surfaces.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
T+02:25  FINAL GEL READY
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Total active session time: ~2 hours 25 minutes
  vs. standard sequential protocol: ~3 hours
  
  Time saved: ~35 minutes
  Quality improvement: significant — all three components
  compounded at peak state simultaneously.
```

---

# PART IV: COMPONENT STATE AT COMPOUNDING — COMPARISON TABLE

```
COMPONENT       SEQUENTIAL PROTOCOL        OPTIMIZED PARALLEL PROTOCOL
                (as written)               (this document)
────────────────────────────────────────────────────────────────────────────────
SAP             80 min rest at compounding  67 min rest at compounding
                PAST optimal window         ★ WITHIN optimal window (30–90 min)
                Syneresis risk beginning    At viscosity and network peak
                
Bentonite       45 min rest at compounding  60 min rest at compounding
                Minimum acceptable          ★ OPTIMAL state
                Just passed threshold        Fully developed thixotropic network
                
Silica          55 min rest at compounding  5–10 min rest at compounding
                PAST optimal window         ★ WITHIN optimal window (10–45 min)
                Early re-agglomeration      Freshest possible dispersion
                occurring

Total session   ~3 hours to gel ready       ~2 hours 25 minutes to gel ready
────────────────────────────────────────────────────────────────────────────────
```

The most significant improvement is the silica state at compounding. In the
sequential protocol, the silica sits for 55 minutes — well into the re-agglomeration
window. In the optimized protocol, it is freshly prepared and incorporated
within its peak window. This directly improves aerogel layer continuity in the
final gel — the single most important quality variable for fire protection performance.

---

# PART V: DEGRADATION INDICATORS — WHAT TO WATCH FOR

## 5.1 — Before Compounding: Component-Level Warnings

```
COMPONENT   DEGRADATION INDICATOR           MEANING / RESPONSE
────────────────────────────────────────────────────────────────────────────────
SAP         Clear water layer on top of gel  Syneresis beginning.
                                             Stir vigorously 60 sec.
                                             Compound within 20 min.
                                             
            Gel pulls away from bowl edge    Network contraction. Same response.
            
            Surface appears dry/matte        Water redistributing internally.
            but gel feels normal below       Minor — stir and proceed.

────────────────────────────────────────────────────────────────────────────────
Silica      White sediment visible at        Settling / re-agglomeration.
            bottom of container              Re-stir 2 min before use.
            
            Suspension visibly thicker       Early network formation.
            at bottom than top               Use immersion blender 2 min.
                                             Compound immediately.
                                             
            Suspension holds its shape       Significant gelation.
            when poured (pours slowly)       Add 50mL water, blend 3 min,
                                             use within 15 min.

────────────────────────────────────────────────────────────────────────────────
Bentonite   Visible grey-tan gel zones       Slight over-gel. 
            separated from clear water       Normal if rested >4 hours.
            on top                           Stir vigorously 90 sec.
                                             Increase compounding whisk
                                             time to 5 minutes.
                                             
            No thixotropy (remains           Under-gelled. Add 10g bentonite,
            completely liquid)               wait 30 min, re-test.
```

## 5.2 — After Compounding: Final Gel Degradation

The final compounded gel has a different degradation profile from the
individual components. Once compounded, the SAP-bentonite-silica network
stabilizes each component relative to its isolated state.

```
FINAL GEL TIME-STATE PROFILE (post-compounding, post-rest):

  0–2 hours post-rest:   ★ OPTIMAL — use for all tests and applications.
                          Network fully integrated. All properties at peak.

  2–6 hours post-rest:   GOOD — minor viscosity drift (slight thickening
                          as bentonite network continues developing).
                          Still fully usable. Stir before application.

  6–12 hours post-rest:  ACCEPTABLE for fire suppression applications.
                          Silica may have partially settled within the gel
                          matrix (visible as slightly clearer layer on top,
                          whiter/denser layer below).
                          Stir thoroughly before use.
                          Aerogel formation may be slightly less uniform.
                          NOT ideal for forming or moisture modification
                          applications where uniformity is critical.

  12–24 hours:           USABLE WITH CAUTION.
                          Seal the bowl with plastic wrap to prevent
                          surface drying. Significant silica settling
                          may have occurred. Stir thoroughly — the settled
                          silica should re-disperse if the gel has not dried.
                          Test on a small sample before full application.

  >24 hours:             SURFACE DRYING PROBABLE.
                          The top layer of the gel will have begun to dry
                          and form a semi-rigid skin (partial aerogel
                          transition at the air interface). This skin
                          must be removed before use — scrape off the
                          top 3–5mm of gel, discard, stir the remaining gel.
                          Performance is significantly reduced.
                          Recommend making a fresh batch.
```

**Practical rule:** Make the batch, use it within 2 hours for critical tests.
If storing overnight for next-day use, seal with plastic wrap, stir before use,
and accept a small quality reduction.

---

# PART VI: CRITICAL TIMING RULES — SUMMARY

These are the derived rules that govern all timing decisions in P1S synthesis:

```
RULE 1 — BENTONITE FIRST, ALWAYS.
  Bentonite has the longest minimum rest requirement (45 min) and improves
  with time up to ~90 min. It must be started first to avoid being the
  bottleneck at compounding time.

RULE 2 — SAP SECOND, IMMEDIATELY AFTER BENTONITE STARTS.
  Starting SAP ~13 minutes after bentonite means SAP reaches its 30-min
  rest minimum at approximately the same time bentonite approaches optimal.
  The two rest periods align.

RULE 3 — SILICA LAST, TIMED TO COMPOUND WITHIN 45 MINUTES.
  The silica suspension quality degrades fastest of all three components.
  It should be prepared last, immediately before compounding is ready to begin.
  Never prepare silica first and let it sit while waiting for SAP or bentonite.

RULE 4 — COMPOUND SAP+BENTONITE BEFORE ADDING SILICA.
  This sequence is non-negotiable. The clay network must establish within
  the SAP matrix before silica is introduced. Silica added to a clay-free
  SAP gel forms different (less uniform) aggregates than silica added to
  a pre-established SAP-bentonite matrix.

RULE 5 — THE SILICA CLOCK STARTS AT END OF DISPERSION, NOT START.
  The quality window (45 min) begins when dispersion is complete and
  confirmed uniform — not when mixing started. Set the timer after the
  final whisk or inspection, not at the beginning of the dispersion step.

RULE 6 — RE-STIR SILICA IMMEDIATELY BEFORE EACH ADDITION.
  Even within the optimal window, stir the silica suspension vigorously
  for 30 seconds before each of the 6 addition portions during compounding.
  Re-agglomeration is continuous; keeping the suspension agitated minimizes
  aggregate size at the point of addition.

RULE 7 — DO NOT LET SAP REST EXCEED 90 MINUTES BEFORE COMPOUNDING.
  Past 90 minutes, syneresis risk becomes real. If the schedule slips and
  SAP has been resting for >90 minutes, stir it thoroughly before compounding
  and inspect for water separation. Record if this occurred.

RULE 8 — USE DISTILLED WATER THROUGHOUT.
  Calcium and magnesium ions from tap water compress the SAP swelling capacity
  by ion-exchange competition and reduce bentonite swelling by interlayer
  cation exchange. These effects are not recoverable by any protocol adjustment.
  Distilled water is not optional — it is what keeps the timing windows above
  valid.
```

---

# PART VII: ONE-PAGE QUICK REFERENCE

```
╔══════════════════════════════════════════════════════════════════════════════╗
║          P1S OPTIMIZED SYNTHESIS SCHEDULE — QUICK REFERENCE                 ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  PRE-WORK (T-05): Measure all water volumes. Set out all equipment.         ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  T+00:00  START BENTONITE (A3)                                              ║
║           30g bentonite + 400mL water. Whisk 3 min. Cover. Rest 60 min.    ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  T+00:13  START SAP (A1)                                                    ║
║           10g SAP + 1000mL water in 10 increments. Rest 30 min.            ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  T+00:55  SAP REST DONE — Check consistency. Do not compound yet.           ║
║           Leave undisturbed. Bentonite needs 18 more minutes.               ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  T+01:00  START SILICA (A2) — N95 ON FIRST                                  ║
║           40g silica + 350mL water. 8 portions, 90 sec whisk each.         ║
║           Final 4-min whisk. 10-min rest + inspect.                        ║
║           ⏱ SILICA QUALITY CLOCK STARTS NOW. Use by T+01:45 latest.       ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  T+01:13  BENTONITE REST DONE — Run thixotropy check.                       ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  T+01:20  ► BEGIN COMPOUNDING (A4)                                          ║
║                                                                              ║
║           SAP state:      67 min rest    ★ OPTIMAL                          ║
║           Bentonite state: 60 min rest    ★ OPTIMAL                          ║
║           Silica state:    ~5-10 min rest ★ OPTIMAL                          ║
║                                                                              ║
║           STEP 1: Bentonite → SAP. Whisk 3 min.                            ║
║           STEP 2: Silica → SAP+bentonite. 6 portions, 2 min each.          ║
║                   Re-stir silica 30 sec before each portion.               ║
║                   Final whisk 6 min (or blender 3 min).                    ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  T+01:55  COMPOUNDING DONE. Final rest 30 min. Do not stir.                 ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  T+02:25  ★ GEL READY — Use within 2 hours for optimal results ★            ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  DEGRADATION QUICK-REFERENCE:                                               ║
║    Water on top of SAP → syneresis → stir and compound within 20 min       ║
║    White sediment in silica → re-agglomeration → blend 2 min, use now      ║
║    No thixotropy in bentonite → add 10g more, wait 30 min                  ║
║    Silica past 45-min window → blend 3 min, compound immediately           ║
╚══════════════════════════════════════════════════════════════════════════════╝
```

---

*End of P1S Synthesis Timing and Chronology Optimization Document*
*Companion to: EHAB Phase 1 Protocol v2.0*
*Author: OrganismCore / Eric Robert Lawson*
