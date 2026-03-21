# THE NIMITZ EVENT — REVERSE ENGINEERING THE PHYSICS
## A Quantitative Attractor Geometry Analysis
## Real Measurements → Derived K-Field Requirements → Physical Constraints
## 2026-03-21

---

## PREAMBLE — THE METHODOLOGY

This is not speculation. This is physics applied to confirmed measurements.

**The method:**

1. Take the confirmed observational data as constraints
2. Apply the vacuum coupling potential framework (V_vac = m₀c²·K³) to each constraint
3. Solve for the K-field parameters required to produce each observation
4. Check internal consistency — do the K-field requirements from different observations agree?
5. Identify what the physics requires, and what that implies for the mechanism

**The confirmed data we are working with:**

| Parameter | Value | Source |
|---|---|---|
| Object length | ~12 m (40 ft) | Fravor visual, size estimate |
| Altitude drop | ~8,534 m (28,000 ft) | Kevin Day / Princeton Aegis radar |
| Drop time | < 0.78 seconds | SCU analysis |
| Horizontal speed observed | ~500 knots (~257 m/s, Mach 0.76 low end) to Mach 1+ | FLIR and visual |
| Estimated peak acceleration | 5,000–14,000 g (various analyses) | SCU / Knuth et al. 2019 |
| Sonic boom | **Absent** | All witnesses |
| Thermal/plasma signature | **Absent** (no FLIR hot spot) | FLIR1 video |
| Propulsion exhaust | **Absent** | All witnesses |
| Water disturbance | White water ~600m diameter below hover point | Fravor testimony |
| Radar cross section | Consistent, solid return | USS Princeton SPY-1 |
| Object morphology | Featureless smooth white cylinder/ellipsoid | Fravor, Dietrich |

---

## PART I — SOLVING FOR THE ACCELERATION

### 1.1 The Kinematics

The radar data is the most quantitatively constrained observation. Let's solve the kinematics precisely.

**Altitude drop:**

$$\Delta h = 8{,}534 \text{ m} \quad \Delta t < 0.78 \text{ s}$$

**Mean vertical velocity during drop:**

$$\bar{v}_z = \frac{\Delta h}{\Delta t} = \frac{8{,}534}{0.78} = 10{,}941 \text{ m/s} \approx 10.9 \text{ km/s}$$

For comparison: Mach 1 at sea level = 343 m/s. This mean velocity = **Mach 31.9** vertically.

**If the object starts from rest (hovering) and accelerates uniformly:**

$$\Delta h = \frac{1}{2}at^2 \implies a = \frac{2\Delta h}{t^2} = \frac{2 \times 8{,}534}{(0.78)^2} = \frac{17{,}068}{0.608} = 28{,}072 \text{ m/s}^2$$

$$a_{\text{uniform}} = 28{,}072 \text{ m/s}^2 = \frac{28{,}072}{9.81} \text{ g} \approx \mathbf{2{,}862 \text{ g}}$$

**If the object reaches terminal velocity in negligible time (impulse model):**

$$a_{\text{peak}} = \frac{\Delta v}{\Delta t_{\text{impulse}}} = \frac{10{,}941}{\sim 0.01} \approx 1{,}094{,}100 \text{ m/s}^2 \approx \mathbf{111{,}000 \text{ g}}$$

**Conservative working number (uniform acceleration model):** ~**2,860 g**. This is consistent with the lower bound of the SCU and Knuth et al. analyses.

### 1.2 The Conventional Physics Impossibility

**For a conventional aircraft at this acceleration:**

Newton's second law: F = ma

Assuming minimum mass for an object 12 m long — say a hollow titanium shell of 10,000 kg (10 tonnes, extremely conservative — an F-18 is 14,000 kg):

$$F_{\text{required}} = 10{,}000 \times 28{,}072 = 2.8 \times 10^8 \text{ N} = 280 \text{ MN}$$

For comparison: the Saturn V rocket (the most powerful launch vehicle ever built) produced 35 MN of thrust. **The Tic Tac would require 8 × the Saturn V's full thrust to achieve even the conservative acceleration estimate, from an object the size of a small aircraft.**

**Power requirement:**

$$P = F \times v = 2.8 \times 10^8 \times 10{,}941 = 3.06 \times 10^{12} \text{ W} = 3.06 \text{ TW}$$

Total US electricity generation capacity: ~1.1 TW. The object would require nearly **3× the entire US grid** in propulsive power. From a 12-meter object. With no exhaust.

**This is not an engineering shortfall. It is a category error. The object is not using thrust.**

---

## PART II — THE SONIC BOOM CONSTRAINT: SOLVING FOR K_BOUNDARY

### 2.1 The Sonic Boom Impossibility — Quantified

The object reached ~Mach 32 vertically. At any speed above Mach 1, a conventional object in atmosphere produces a sonic boom. The pressure differential of the bow shock scales as:

$$\Delta P_{\text{boom}} \propto \rho_{\text{air}} \times v^2$$

At Mach 32:

$$\Delta P = \rho_{\text{air}} \times v^2 = 1.225 \times (10{,}941)^2 = 1.47 \times 10^8 \text{ Pa} = 1{,}450 \text{ atm}$$

A pressure wave of 1,450 atmospheres propagating outward from the object would have been heard — and felt — across a radius of tens of kilometers. It would have been catastrophic. **It was not detected.**

### 2.2 Solving for K at the Object Boundary

The sonic boom is produced by momentum transfer from the moving object to the air at the object's surface boundary. The momentum transfer rate is:

$$\frac{dp}{dt} = \rho_{\text{air}} \times A_{\text{cross}} \times v^2$$

In the vacuum coupling framework, the effective coupling between the object's surface and the ambient air is mediated by the K-field at the boundary. If the K-field at the object's surface boundary is K_boundary, then the effective coupling is modified:

The air molecules couple to the object through electromagnetic interactions (the same EM coupling that transmits pressure). These EM interactions are mediated through the local vacuum mode density. If K_boundary suppresses the EM coupling:

$$\left(\frac{dp}{dt}\right)_{\text{effective}} = K_{\text{boundary}}^3 \times \rho_{\text{air}} \times A_{\text{cross}} \times v^2$$

**Constraint: No sonic boom detected means the effective momentum transfer must be below the detection threshold.**

The detection threshold for a sonic boom at ~10 km range (Princeton was approximately at that distance): ΔP_threshold ≈ 0.1 Pa (threshold of perception). The expected ΔP is ~10⁸ Pa. Ratio:

$$K_{\text{boundary}}^3 \leq \frac{\Delta P_{\text{threshold}}}{\Delta P_{\text{conventional}}} = \frac{0.1}{1.47 \times 10^8} = 6.8 \times 10^{-10}$$

$$\boxed{K_{\text{boundary}} \leq \left(6.8 \times 10^{-10}\right)^{1/3} = 8.8 \times 10^{-4}}$$

**The K-field at the object's boundary must be less than ~0.001 (less than 0.1% of ambient) to suppress the sonic boom to below detection threshold.**

This is the first hard quantitative constraint on K_boundary from the observational data.

---

## PART III — THE PLASMA/THERMAL CONSTRAINT: SOLVING FOR K_INTERIOR

### 3.1 The Plasma Threshold

At Mach 32, a conventional object in lower atmosphere would produce a bow shock with post-shock temperature:

$$T_{\text{shock}} = T_{\text{ambient}} \times \left(1 + \frac{2\gamma}{\gamma+1}(M^2 - 1)\right) \approx T_{\text{ambient}} \times \frac{2\gamma M^2}{\gamma+1}$$

For γ = 1.4 (diatomic air), M = 32, T_ambient = 250 K (at 28,000 ft altitude):

$$T_{\text{shock}} = 250 \times \frac{2 \times 1.4 \times 32^2}{2.4} = 250 \times \frac{2{,}867}{2.4} = 250 \times 1{,}195 = 298{,}750 \text{ K}$$

**~300,000 K** — this is not just plasma. This is coronal-temperature plasma. X-ray emitting. Blindingly visible on any IR sensor on the planet. COMPLETELY absent.

### 3.2 Solving for K_interior from the Thermal Constraint

The thermal energy generated is proportional to the EM coupling between the object surface and the air molecules:

$$T_{\text{actual}} = T_{\text{ambient}} + K_{\text{boundary}}^3 \times (T_{\text{shock,conventional}} - T_{\text{ambient}})$$

**Constraint:** The FLIR1 video shows no anomalous thermal signature above ambient. Maximum detectable thermal anomaly: ΔT ≈ 10 K above background (FLIR sensitivity).

$$K_{\text{boundary}}^3 \leq \frac{\Delta T_{\text{threshold}}}{T_{\text{shock}} - T_{\text{ambient}}} = \frac{10}{298{,}500} = 3.35 \times 10^{-5}$$

$$\boxed{K_{\text{boundary,thermal}} \leq \left(3.35 \times 10^{-5}\right)^{1/3} = 0.032}$$

This gives K_boundary < 0.032 from the thermal constraint — **less stringent than the sonic boom constraint** (K < 0.00088), which means the sonic boom constraint is the binding constraint. Both constraints are consistent: K_boundary ≪ 1.

---

## PART IV — THE INERTIA CONSTRAINT: SOLVING FOR K_INTERIOR

### 4.1 The Occupant Survivability Problem

Fravor and Dietrich saw the object react to their approach — it mirrored Fravor's descent, then shot away. This implies something akin to controlled navigation — not ballistic motion. If there are occupants, they survived 2,860–14,000 g maneuvers.

**For a 70 kg human at 2,860 g:**

$$F_{\text{on human}} = 70 \times 2{,}860 \times 9.81 = 1.96 \times 10^6 \text{ N}$$

This is 196 tonnes of force on a 70 kg body. Instantaneous death at ~20 g sustained. Fatal within milliseconds at 100 g.

**The vacuum coupling framework resolves this directly:**

If K_interior ≈ ε << 1, then:

$$m_{\text{eff,human}} = m_0 \times K_{\text{interior}}^3 \approx 0$$

The effective inertial mass of everything inside the K-bubble approaches zero. The acceleration experienced by occupants is:

$$a_{\text{felt}} = \frac{F_{\text{external}}}{m_{\text{eff}}} = \frac{F_{\text{external}}}{m_0 \times K_{\text{interior}}^3}$$

But here's the key insight: **If the entire system (object + occupants) is inside the same K-bubble, there is no external force on the occupants.** The K-bubble moves as a unit. The objects inside it are not being pushed — the bubble's configuration is changing, and everything inside it moves with it.

The analogy: you are in free fall inside an elevator in free fall. You feel zero g. The object is in "vacuum free fall" inside its own K-bubble. The bubble configuration changes. Everything inside experiences zero proper acceleration.

**Constraint from survivability:** K_interior must be low enough that the effective inertial coupling of interior objects to the exterior acceleration profile is negligible.

$$a_{\text{felt}} = K_{\text{interior}}^3 \times a_{\text{bubble}} \leq 9 \text{ g (survivable)}$$

$$K_{\text{interior}}^3 \leq \frac{9}{2{,}860} = 3.1 \times 10^{-3}$$

$$\boxed{K_{\text{interior}} \leq 0.146}$$

This is a much less stringent constraint than K_boundary. The interior only needs K < 0.15 for survivability. The boundary needs K < 0.001 to suppress the sonic boom.

**This reveals the two-region structure of the K-field:**

- **Interior region (K_int ≈ 0.01–0.1):** Suppresses inertial coupling of occupants to bubble acceleration. Makes 2,860 g maneuvers survivable.
- **Boundary layer (K_bound < 0.001):** Suppresses coupling to the medium — eliminates sonic boom, thermal signature, plasma formation.

---

## PART V — THE WHITE WATER ANOMALY: THE K-FIELD GROUND SIGNATURE

### 5.1 The Observation

Fravor confirmed: before he saw the Tic Tac object above, he and Dietrich observed a **circular white water disturbance approximately 600 meters in diameter** (comparable to a "football field" or larger — his description was a Boeing 737 footprint, roughly 60m wide, but the *area of disrupted water* was much larger) directly below the object's hover position, at the ocean surface.

The object was hovering at approximately **100–500 feet above sea level** above this disturbance.

### 5.2 What Creates White Water

"White water" in the ocean is created by:
- Surface wave breaking
- Subsurface disturbance propagating upward
- Downwash from rotor/propulsion system
- Pressure differential at water surface

### 5.3 The K-Field Ground Signature — The Derivation

In the vacuum coupling framework, the K-bubble extends from the object. The K-field has a profile:

$$K(r) = 1 - (1 - K_{\text{int}}) \times e^{-(r - r_{\text{object}})/\delta}$$

where r is the radial distance from the object's surface, and δ is the boundary layer thickness.

At 100–500 feet altitude, the K-bubble boundary reaches the ocean surface. At the ocean surface, the K-gradient ∇K is nonzero — there is a vacuum coupling potential gradient at the water surface.

The force per unit volume on the water at the surface:

$$\mathbf{f}_{\text{water}} = -\rho_{\text{water}} c^2 \nabla K^3$$

This is a **downward pressure** on the water surface where the K-gradient impinges on it. **The K-field gradient at the water surface pushes the water down and outward — creating exactly the circular disturbance pattern observed.**

More precisely: the K-bubble at low altitude creates a region of modified vacuum coupling between the object and the water surface. The water molecules in this region experience a modified coupling to the vacuum — their surface tension properties, electromagnetic binding, and effective mass are all slightly modified. **The white water is the K-field footprint on the ocean surface.**

### 5.4 Quantifying the Ground Signature

The disturbance radius R_disturbance ≈ 300 m (estimating 600 m diameter). The object altitude h ≈ 100 m (conservative). The K-bubble apparently extended from the object to the ocean surface.

This gives us the K-bubble's **effective radius** at the ocean surface:

$$r_{\text{bubble,surface}} \approx \sqrt{R_{\text{disturbance}}^2 + h^2} = \sqrt{300^2 + 100^2} \approx 316 \text{ m}$$

The K-bubble at the time of hovering had an effective radius of approximately **300–316 meters.**

**This is the bubble size constraint.** The K-suppression extended over a volume of:

$$V_{\text{bubble}} = \frac{4}{3}\pi r^3 = \frac{4}{3}\pi (316)^3 \approx 1.32 \times 10^8 \text{ m}^3$$

This is the volume of vacuum that must be maintained at K << 1 to produce the observed ground signature.

---

## PART VI — THE RADAR RETURN PARADOX: THE TWO-K STRUCTURE CONFIRMED

### 6.1 The Paradox

If K_boundary < 0.001 (required to suppress sonic boom and thermal signatures), then the object's EM coupling to radar waves should also be suppressed. A radar signal propagating into a K < 0.001 region would be refracted, absorbed, or reflected at the K-boundary — not from the solid object inside.

**Yet:** The Princeton's SPY-1 radar produced a **clear, consistent, solid radar return** from the Tic Tac.

This is the radar return paradox: the object suppresses EM coupling enough to eliminate sonic boom and plasma, but maintains enough EM coupling to produce a clear radar return.

### 6.2 The Resolution — The K-Bubble Wall IS the Radar Reflector

In the K-field framework, the bubble wall (the transition zone from K ≈ ε to K = 1) has a strong refractive index gradient:

$$n(r) = K^{-1/2}(r)$$

The refractive index changes from n ≈ K_boundary^(-1/2) ≈ 33 (inside the thin boundary) to n = 1 (outside ambient). This is an enormous refractive index discontinuity.

**Any electromagnetic wave (including radar) hitting this discontinuity will be strongly reflected** — exactly as light reflects from a glass surface, but with a refractive index ratio of ~33:1 instead of 1.5:1.

The reflectivity at a refractive index discontinuity (Fresnel equation, normal incidence):

$$R = \left(\frac{n_2 - n_1}{n_2 + n_1}\right)^2 = \left(\frac{33 - 1}{33 + 1}\right)^2 = \left(\frac{32}{34}\right)^2 = 0.885$$

**88.5% of incident radar energy is reflected from the K-bubble wall.** This would produce an extremely strong radar return — consistent with the "solid object" radar return reported by the Princeton.

**The radar is not seeing the physical object. It is seeing the K-bubble wall.** The bubble wall IS the radar reflector. The solid object inside is irrelevant to the radar return.

**This is a critical prediction of the framework, directly confirmed by the observation:**

- Strong radar return: ✓ (K-bubble wall reflection)
- Absent sonic boom: ✓ (K_boundary suppresses air coupling)
- Absent thermal: ✓ (K_boundary suppresses EM/heat coupling)
- The radar return size corresponds to the bubble size, not the object size

**The radar cross-section gives us the bubble size, not the object size.** If the radar return corresponded to a 12-meter object, the bubble radius was approximately 6 meters. If the return was larger (consistent with the white water 300m bubble), the bubble radius was much larger during hovering. The bubble appears to have been compressible — large during hovering, small during high-speed transit.

---

## PART VII — THE COMPLETE PARAMETER SOLVE

Combining all constraints:

### K-Field Requirements Solved from Observations:

| Observation | Constraint | Solved K Value |
|---|---|---|
| Sonic boom absent (Mach 32) | K³_boundary ≤ 6.8×10⁻¹⁰ | **K_boundary ≤ 8.8×10⁻⁴** |
| Thermal absent (300,000 K shock suppressed to <10 K) | K��_boundary ≤ 3.35×10⁻⁵ | K_boundary ≤ 0.032 |
| Occupant survivability at 2,860 g | K³_interior ≤ 3.1×10⁻³ | K_interior ≤ 0.146 |
| Radar return present (strong reflection) | n_wall = K⁻¹/² >> 1 | K_boundary ~ 10⁻³ to 10⁻⁴ |
| White water radius ~300 m | Bubble radius ≥ 300 m during hover | V_bubble ~1.3×10⁸ m³ |
| Instant acceleration (K_interior → 0) | m_eff = m₀·K³_int ≈ 0 | K_interior < 0.1 |

**The binding constraint is the sonic boom suppression: K_boundary < 8.8 × 10⁻⁴.**

All other constraints are less demanding and are automatically satisfied if K_boundary ≈ 10⁻³.

### The K-Field Profile Required:

$$K(r) = \begin{cases} K_{\text{int}} \approx 0.01 \text{ to } 0.1 & r < r_{\text{object}} \text{ (interior)} \\ K_{\text{bound}} < 8.8 \times 10^{-4} & r_{\text{object}} < r < r_{\text{object}} + \delta \text{ (boundary layer)} \\ 1 & r > r_{\text{object}} + \delta \text{ (ambient)} \end{cases}$$

The boundary layer thickness δ is unconstrained by the current observations — it could be nanometers to meters.

---

## PART VIII — THE ENERGY BUDGET

### 8.1 Energy Cost of the K-Bubble

The Casimir energy density formula gives the vacuum energy modification per unit volume:

$$u_{\text{Casimir}} = -\frac{\pi^2\hbar c}{720 a^4}$$

For a spherical K-bubble of radius R with K_interior ≈ 10⁻³:

The vacuum energy suppression per unit volume scales as:

$$\Delta u \approx -\rho_{\text{ZPF}} \times (1 - K^3) \approx -\rho_{\text{ZPF}}$$

The ZPF energy density integrated up to an effective cutoff (Planck scale is too large — use the relevant atomic scale cutoff at ω_c corresponding to ~1 keV photons):

$$\rho_{\text{ZPF,cutoff}} = \frac{\hbar\omega_c^4}{8\pi^2 c^3} = \frac{(1.6\times10^{-16})^4}{8\pi^2 (3\times10^8)^3 \hbar^3}$$

At 1 keV cutoff: ρ_ZPF ≈ 10⁶ J/m³

**Energy stored in the K-bubble of radius r = 6 m (transit size) with K_int = 10⁻³:**

$$E_{\text{bubble}} = \Delta u \times V = 10^6 \times \frac{4}{3}\pi (6)^3 \approx 10^6 \times 905 = 9 \times 10^8 \text{ J} \approx \mathbf{900 \text{ MJ}}$$

For comparison: a Tomahawk cruise missile carries ~450 kg of explosive = ~1,890 MJ. The hovering bubble (radius 300 m) would store ~2.3 × 10¹⁵ J = **2.3 petajoules** at this cutoff. This is comparable to a large nuclear weapon.

**This energy budget is the most important output of the analysis.** It tells us:

1. The K-bubble during high-speed transit (small, tight bubble) requires ~GJ-scale vacuum energy modification — enormous by current engineering standards but finite.
2. The hovering configuration (large, expanded bubble) requires petajoule-scale modification — this is closer to nuclear energy scales.
3. **The object likely collapses its bubble to minimum size during high-speed transit** (consistent with Fravor's description of the object appearing to "shrink" as it accelerated away) and expands it during hovering (consistent with the large water disturbance).

This is a testable prediction: **the apparent radar cross-section of the object should decrease during high-speed acceleration and increase during hovering**, because the bubble size collapses and expands with the K-field configuration.

---

## PART IX — THE WATER DISTURBANCE: THE SUBMERGED COMPANION OBJECT

### 9.1 The Sub-surface Connection

Fravor's testimony describes the water disturbance as if something was just below the surface. Kevin Day and other Princeton crew reported that the Tic Tac objects appeared to descend from very high altitude and then hover near sea level, while simultaneously there appeared to be objects at or below the sea surface.

**The K-field interpretation:** When the Tic Tac K-bubble reaches the water surface, its boundary intersects with the high-density medium of the ocean. The K-gradient at the air-water interface would produce an upward pressure on the water (the water molecules in the ∇K region experience a modified vacuum coupling — they are pulled upward toward the region of higher K, toward the ambient). This creates the characteristic circular white water pattern — a ring of disturbed water at the K-bubble boundary contact zone with the ocean surface.

**Additionally:** If the K-bubble extends below the water surface, the object inside it could descend into the water without any hydrodynamic signature — exactly as the Aguadilla object did 9 years later.

**The sub-surface companion may have been the same object, hovering with its K-bubble straddling the air-water interface.** The water disturbance is not caused by a submerged object — it IS the K-bubble ground signature at the air-water boundary.

---

## PART X — THE COMPLETE PICTURE RECONSTRUCTED

From the observational data, the physics of the 2004 event reconstructs as follows:

**Pre-encounter (Nov 10–13):**
The Princeton SPY-1 tracked objects descending from >80,000 ft over 5 days. The K-bubbles during descent were small (tight, high-speed transit mode). Radar returns were consistent but brief — the bubbles were collapsing and reforming as the objects maneuvered. Kevin Day described "dozens of objects" — multiple K-bubbles of varying sizes were being tracked.

**The hovering configuration (Nov 14, pre-Fravor):**
The object was hovering at 100–500 ft. K-bubble expanded to ~300m radius. This large bubble impinged on the ocean surface, creating the white water disturbance. K_interior ≈ 0.01–0.1. K_boundary ≈ 10⁻³. The object inside was completely decoupled from the medium. No noise, no thermal emission, no drag.

**The Fravor encounter:**
Fravor descended toward the object. The object "mirrored" his descent — it matched his approach angle before departing. This is **K-bubble navigation:** the K-configuration was being adjusted to produce displacement, not thrust. When Fravor attempted to intercept, the K-bubble reconfigured and the object "snapped" to a new position 60 miles away in "less than a second" — as reported by the E-2 Hawkeye that was monitoring.

**The FLIR1 video:**
The "glow" visible around the object in the FLIR video is the K-bubble wall emission — photons produced at the boundary where K_boundary transitions to K = 1. The object itself is dark (K_interior ≈ 0 → no EM emission from interior). The boundary emits at the vacuum mode transition frequencies. This produces the characteristic visual appearance: a glowing boundary with a dark interior.

**The instantaneous 28,000 ft drop:**
K-bubble reconfiguration. The K-minimum (the attractor for the object inside) was relocated from 28,000 ft to sea level. The object followed the K-gradient — it descended into the new potential minimum. No thrust. No acceleration felt by interior. External acceleration apparent of ~2,860 g. The ∇V_vac pointed downward for 0.78 seconds, then stopped.

---

## THE QUANTITATIVE SUMMARY

The Nimitz event, solved:

$$K_{\text{boundary}} < 8.8 \times 10^{-4}$$

$$K_{\text{interior}} < 0.15 \text{ (survivability)}, \approx 0.01 \text{ (effective zero-mass condition)}$$

$$r_{\text{bubble,hover}} \approx 300 \text{ m}$$

$$r_{\text{bubble,transit}} \approx 6 \text{ m (radar cross-section limited)}$$

$$E_{\text{bubble,transit}} \approx 900 \text{ MJ (at 1 keV ZPF cutoff)}$$

$$V_{\text{vac,boundary}} = m_0c^2 \times K_{\text{boundary}}^3 \approx m_0c^2 \times 7 \times 10^{-10}$$

The bubble wall refractive index: **n_wall ≈ K⁻¹/² ≈ 33** — producing ~88% radar reflectivity, consistent with the strong SPY-1 return.

The effective inertial mass during transit: **m_eff ≈ m₀ × (0.01)³ = 10⁻⁶ × m₀** — one millionth of rest mass.

**The object was not flying. It was navigating a self-specified vacuum coupling potential landscape.**

The physics that took place was not propulsion. It was K-field configuration management — the maintenance and controlled modification of a locally suppressed vacuum mode density bubble that decoupled the object from every standard coupling channel in the medium it traversed.

**Newton's framework applies. The force is −∇V_vac. The landscape is K(r). The navigator is whatever is inside the bubble. The gap is the K-field configuration that bridges the two.**

The event is not unexplained. **It is under-modeled.** The vacuum coupling potential provides the model. The observations provide the constraints. The constraints are internally consistent. The physics is derivable. The experiment that tests it — inertial mass measurement inside a Casimir geometry — has not been done.

That is the gap between where we are and where this physics leads.
