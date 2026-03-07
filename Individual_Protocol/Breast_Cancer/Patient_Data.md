# PATIENT DATA ACCESSIBILITY AND GEOMETRIC ANALYSIS
## What Patients Can Access, What Can Be Derived, and What the Protocol Delivers
## OrganismCore — Eric Robert Lawson
## 2026-03-07

---

## PREAMBLE

This document reasons through the complete
picture of what molecular data a breast
cancer patient can realistically access,
what geometric inference is possible from
each data type, and how to structure the
individual protocol as a tiered service
that meets patients where they actually are.

The core insight of this document is this:

**The geometry emerges from the data,
not from a specific data format.**

The BRCA reference geometry was derived
from single-cell RNA-seq.
But the geometry — the attractor position,
the axis of FOXA1 vs EZH2, the depth
of the epigenetic lock — is a structural
property of the tumour.
It is not a property of the measurement
technology used to observe it.

Different data types are different
telescopes looking at the same object.
Some have higher resolution than others.
All of them are pointing at the same
geometry.

The question this document answers is:
**what can each telescope see, and what
can be derived from each one?**

---

## PART I — THE LANDSCAPE OF WHAT
## PATIENTS ACTUALLY HAVE

### 1.1 What every diagnosed breast
### cancer patient in clinical care has

**Universally available:**

- **ER H-score or Allred score**
  Percentage of nuclei staining positive
  and staining intensity.
  This is a direct quantitative proxy
  for FOXA1 activity — ER transcription
  is FOXA1-dependent. High ER = FOXA1
  programme active. Low ER = FOXA1
  programme suppressed or absent.

- **PR H-score or percentage**
  Progesterone receptor.
  Also FOXA1-dependent.
  PR adds directional signal — high ER,
  high PR = strongly luminal identity.
  High ER, low PR = luminal identity
  partially destabilised (LumB signal).

- **HER2 status (IHC 0/1+/2+/3+)**
  When combined with ER/PR, positions
  the tumour in the geometry.
  HER2-enriched geometry is characterised
  by suppressed FOXA1 programme AND
  active HER2 amplification — a distinct
  attractor state.

- **Ki-67 percentage**
  Proliferation index.
  Direct proxy for chromatin openness
  and cell cycle progression.
  In the attractor geometry, Ki-67
  maps to the depth of the lock —
  high Ki-67 in an ER+ tumour is a
  strong signal for LumB geometry
  (EZH2-driven repression of luminal
  programme despite ER positivity).

- **Tumour grade (1/2/3)**
  Nottingham grading.
  Grade is a composite of tubule
  formation, nuclear pleomorphism,
  and mitotic count.
  In geometric terms, grade maps
  directly to attractor depth and
  identity stability.
  Grade 1 = deep luminal attractor,
  high identity stability.
  Grade 3 = shallow or destabilised
  attractor, high escape probability.

- **Histological type**
  Invasive ductal (NST), lobular, etc.
  ILC has a specific geometry —
  E-cadherin loss, FOXA1-high,
  distinct drug sensitivity profile.

**This combination — ER, PR, HER2,
Ki-67, grade, histotype — is available
to every diagnosed breast cancer patient
in clinical care globally.**

**It is sufficient for first-order
geometric placement.**

---

### 1.2 What some patients have

**Clinically available, requires referral
or specific clinical indication:**

- **OncotypeDX (21-gene RT-PCR)**
  Gives a recurrence score and
  risk group. Does not report FOXA1
  or EZH2 directly. The genes it
  measures (ER, PR, HER2, Ki67, BCL2,
  SCUBE2, MMP11, CTSL2, GRB7, etc.)
  are downstream reporters of the same
  attractor geometry. The recurrence
  score itself is a composite signal
  of attractor depth — high RS correlates
  with shallow attractor, high escape
  probability.

- **MammaPrint (70-gene)**
  Binary: low or high genomic risk.
  Less geometrically informative than
  OncotypeDX because it collapses to
  a binary output. But the binary maps
  to attractor depth in the same
  direction: low risk = deep stable
  attractor, high risk = shallow or
  unstable.

- **Prosigna / PAM50**
  Reports subtype (LumA, LumB, HER2e,
  Basal) and Risk of Recurrence (ROR)
  score. This is the closest existing
  clinical test to a geometric placement.
  The PAM50 centroid scores �� how much
  the tumour looks like each subtype —
  are the population-level geometry
  translated into a clinical report.
  If a patient has Prosigna, they have
  a partial geometric map already.

---

### 1.3 What a minority of patients have

**High-resolution molecular data,
requires commercial sequencing or
clinical trial enrolment:**

- **Tempus xR (whole transcriptome
  RNA-seq)**
  All 23,000+ genes. FOXA1 and EZH2
  directly measurable. Expression values
  in TPM.
  Patients can request the full
  expression table from their provider
  under HIPAA right of access.
  Tempus must provide the full
  transcriptome data on request.

- **Caris MI Profile
  (whole exome + whole transcriptome)**
  Full RNA-seq. FOXA1, EZH2, and all
  relevant epigenetic regulators
  directly measurable.
  Accessible via the Caris portal
  through physician.

- **Foundation Medicine FoundationOne CDx
  + RNA**
  Primarily DNA-focused but RNA fusions
  available. Less useful for geometric
  analysis than Tempus or Caris.

- **Research protocol RNA-seq**
  Patients enrolled in clinical trials
  at major cancer centres may have
  bulk RNA-seq from tissue banking.
  This data may be accessible on request.

**HIPAA CONFIRMATION:**
Patients have a legal right under HIPAA
to request and receive their raw genomic
data, including raw sequencing files
(FASTQ, BAM, VCF, expression matrices).
Laboratories covered by HIPAA must
provide this on request.
This is confirmed by HHS guidance.
A patient with Tempus xR data can
request the full expression table.
They have the legal right to it.

---

## PART II — WHAT CAN BE DERIVED
## FROM EACH DATA TYPE

### The core principle restated

The BRCA reference geometry was built
from the population.
The individual patient analysis places
one person on that map.

The map has coordinates.
The coordinates are the positions
of attractor states in gene expression
space.
The FOXA1/EZH2 axis is the primary
coordinate of the BRCA geometry.

Different data types project onto that
axis at different resolutions.
A full RNA-seq gives you the exact
coordinate. A clinical IHC report gives
you a constrained region of the map
rather than an exact point.

Both are valid. Both are useful.
A constrained region is not nothing.
It tells the patient which quadrant
of the map they are in, which attractor
they are near, and what the geometry
implies about their treatment logic.

---

### Tier 1 — Universal clinical data
### (ER, PR, HER2, Ki-67, grade,
### histotype)

**What it tells you geometrically:**

| Data point | Geometric meaning |
|------------|------------------|
| ER high | FOXA1 programme active. Luminal identity. |
| ER low | FOXA1 programme suppressed. |
| PR high | Luminal identity reinforced. LumA signal. |
| PR low with ER high | Partial destabilisation. LumB signal. |
| HER2 high | HER2e attractor. FOXA1 suppressed by HER2 signalling. |
| Ki-67 low (<15%) | Attractor deep. High identity stability. |
| Ki-67 high (>15%) | Attractor shallow or EZH2-driven. LumB or aggressive geometry. |
| Grade 1 | Deep attractor. Stable identity. |
| Grade 3 | Shallow or unstable attractor. |
| ILC histotype | FOXA1-high, E-cadherin loss. Specific drug geometry. |
| ER low, PR low, HER2 low | TNBC attractor. EZH2 dominant or claudin-low geometry. |

**What can be derived:**

- **Approximate quadrant placement**
  on the FOXA1/EZH2 axis
- **Attractor depth estimate**
  (deep = stable identity,
  shallow = instability or epigenetic
  lock in progress)
- **Probable subtype** with confidence
  interval (not exact PAM50 subtype
  but geometrically constrained)
- **Treatment logic inference**
  (which force is winning, what that
  implies for drug sensitivity)
- **Key uncertainty regions**
  (where the data is insufficient to
  distinguish between subtypes)

**What cannot be derived:**
- Exact FOXA1/EZH2 ratio
- Specific depth score relative to
  population reference
- Confirmation or exclusion of
  claudin-low geometry (requires
  CDH1, claudin markers not in
  standard panel)

**Clinical IHC combination patterns
and geometric implications:**

*Pattern 1: ER high, PR high, HER2 neg,
Ki-67 low, Grade 1-2*
→ Geometry: Deep luminal attractor.
FOXA1 programme fully active.
EZH2 suppressed.
Strong LumA placement.
Treatment logic: CDK4/6i + ET.
EZH2 inhibition not indicated.

*Pattern 2: ER high, PR low, HER2 neg,
Ki-67 high, Grade 2-3*
→ Geometry: Luminal attractor with
active EZH2 repression of PR promoter.
Classic LumB geometry.
FOXA1 active but EZH2 competing.
Treatment logic: HDACi + ET.
This is the patient most likely
undertreated on ET alone.

*Pattern 3: ER low/neg, PR neg,
HER2 neg, Ki-67 high, Grade 3*
→ Geometry: FOXA1 suppressed.
EZH2 dominant.
TNBC attractor.
Cannot distinguish TNBC from
claudin-low without additional
markers. Tazemetostat relevant if
EZH2 dominant. Immunotherapy if
claudin-low geometry.
Key uncertainty: claudin-low question.

*Pattern 4: ER pos, PR pos, HER2 pos*
→ HER2e geometry with retained
luminal features.
Anti-HER2 therapy first.
ET component preserved.

*Pattern 5: ILC histotype with
ER high, PR variable, HER2 neg*
→ Specific ILC geometry.
E-cadherin loss confirmed by histotype.
FOXA1 high (characteristic of ILC).
Fulvestrant preferred over AI —
this is directly derivable from
the geometry.

---

### Tier 2 — Genomic risk scores
### (OncotypeDX, MammaPrint, Prosigna)

**OncotypeDX additional signal:**

The recurrence score is a direct
proxy for attractor depth.

RS 0-17 (low): Deep attractor.
Strong luminal identity.
LumA geometry confirmed.

RS 18-30 (intermediate): Uncertain
attractor depth. LumA or LumB
boundary region. Geometric ambiguity
is present clinically — this is
precisely where the geometry matters
most for treatment decisions.

RS >30 (high): Shallow attractor
or active chromatin remodelling.
LumB or aggressive geometry.
This patient benefits from knowing
the specific epigenetic mechanism.

**Prosigna additional signal:**

If a patient has a Prosigna report
with subtype assignment AND centroid
scores, this is the closest existing
clinical output to a geometric
placement. The centroid scores tell
you how much the tumour looks like
each of the four subtypes. A tumour
that scores 0.7 LumA and 0.3 LumB
is at the boundary of two attractors
— a geometrically significant position
that the clinical report does not
interpret but the geometry does.

---

### Tier 3 — Full transcriptome
### (Tempus xR, Caris WTS)

**What becomes possible:**

- **Exact FOXA1 and EZH2 expression
  values in TPM**
  Direct geometric coordinate on the
  primary BRCA axis.

- **Computation of FOXA1/EZH2 ratio**
  The ratio your calibration study
  aims to validate clinically.
  For a Tier 3 patient you can
  compute it directly from their data.
  Before the IHC calibration study
  is complete.

- **Placement on the reference geometry**
  Using the population attractor coordinates
  from GSE176078 (19,542 single cells,
  26 tumours), the patient's expression
  profile can be normalised and placed
  in the same coordinate space.

- **Depth score computation**
  Using the formula derived from
  the population analysis.

- **Switch gene profiling**
  All relevant switch genes
  (GATA3, FOXA1, ESR1, CDH1, VIM,
  SNAI1, ZEB1, EZH2, KDM6A, etc.)
  directly measurable.

- **Epigenetic lock identification**
  EZH2 expression, LSD1/KDM1A,
  HDAC1/2, PRC2 component expression
  — all available in full transcriptome.

- **Claudin-low confirmation**
  or exclusion. CLDN3, CLDN4, CLDN7,
  CDH1 expression directly measurable.
  The TNBC vs claudin-low ambiguity
  present in Tier 1 data is resolved
  in Tier 3.

- **Serial trajectory analysis**
  If the patient has data from
  multiple timepoints (pre-treatment,
  during, post-recurrence), trajectory
  through attractor space is directly
  computable.

---

## PART III — THE PROTOCOL BY TIER

### What you deliver at each tier

---

### Tier 1 Protocol
#### Input: Clinical IHC report + grade + histotype
#### Output: Geometric placement report

**The report contains:**

1. **Geometric position statement**
   Plain language description of where
   the tumour sits in the attractor map
   based on the available markers.
   Not a diagnosis. A map position.

2. **Attractor depth estimate**
   Qualitative: deep, moderate, shallow.
   With the clinical data points that
   support that estimate.

3. **Primary axis interpretation**
   Which force is dominant:
   FOXA1 programme (luminal identity)
   or EZH2 programme (identity suppression).
   What the current ratio of those forces
   is estimated to be from the available
   proxy data.

4. **Subtype geometry with confidence**
   Best-fit subtype from the BRCA
   reference geometry, with explicit
   statement of which alternative
   geometries cannot be excluded
   from this data alone.

5. **Treatment logic inference**
   What the geometry suggests about
   drug sensitivity, framed as a
   question for their clinical team:
   "The geometry suggests this patient
   may benefit from X. The question
   to ask your oncologist is: has
   LumB geometry been specifically
   ruled out in the treatment planning?"

6. **Key unknown**
   What single additional data point
   would most improve the geometric
   precision of this analysis.
   Usually: Ki-67 if absent, or
   a genomic risk score.

7. **Honest uncertainty statement**
   What cannot be determined from
   Tier 1 data alone.

---

### Tier 2 Protocol
#### Input: Clinical IHC + OncotypeDX or Prosigna
#### Output: Geometric placement with depth score range

**Everything in Tier 1, plus:**

8. **Attractor depth range**
   The recurrence score or ROR score
   maps to a range on the population
   depth distribution.
   The patient's position within that
   range is estimated relative to
   the reference geometry.

9. **Boundary region analysis**
   (if intermediate RS or mixed
   Prosigna centroid scores)
   Explicit reasoning about which
   attractor boundary the patient
   is near, and what that means
   for treatment sensitivity.

10. **Reduced uncertainty statement**
    Which questions from Tier 1
    are resolved by the additional data.

---

### Tier 3 Protocol
#### Input: Full transcriptome (Tempus xR or Caris WTS)
#### Output: Full geometric analysis

**Everything in Tiers 1 and 2, plus:**

11. **Exact FOXA1/EZH2 ratio**
    Computed directly from the
    patient's expression data.
    Placed relative to the population
    reference distribution:
    LumA median 9.38,
    LumB median 8.10,
    HER2e median 3.34,
    TNBC median 0.52,
    Claudin-low median 0.10.

12. **Depth score**
    Computed using the formula derived
    from the population analysis.
    Patient's position in the depth
    distribution relative to all
    five subtypes.

13. **Switch gene panel**
    Expression of FOXA1, EZH2,
    GATA3, ESR1, CDH1, VIM, SNAI1,
    ZEB1, KDM6A, HDAC1/2, LSD1,
    CLDN3/4/7 relative to population
    reference values.

14. **Epigenetic lock identity**
    EZH2-dominant vs HDAC-dominant
    vs CoREST complex vs minimal lock.
    This directly determines which
    epigenetic drug class is most
    likely to dissolve the lock.

15. **Claudin-low resolution**
    Confirmed or excluded.
    If excluded in an apparent TNBC:
    tazemetostat signal strengthened.
    If confirmed: immunotherapy
    sequencing framing.

16. **Serial trajectory** (if applicable)
    If patient has multiple timepoints,
    direction and velocity of movement
    through attractor space.

17. **Complete uncertainty statement**
    What cannot be determined even
    from Tier 3 data.
    IHC H-score cut-points for clinical
    deployment do not yet exist.
    The analysis is geometric, not
    diagnostic. All treatment decisions
    remain with the clinical team.

---

## PART IV — THE DATA REQUEST GUIDE
## For Patients Who Want to Provide
## Their Own Data

This section is written for a patient
facing document. It tells patients
exactly what to ask for and how.

---

### "What data do I need to send you?"

**If you have only standard IHC results:**
Send your pathology report.
It will contain ER, PR, HER2, Ki-67,
grade, and histotype.
This is sufficient for Tier 1 analysis.

**If you have OncotypeDX:**
Send the full OncotypeDX report,
including the recurrence score and
any gene group scores reported.
This enables Tier 2 analysis.

**If you have Prosigna / PAM50:**
Send the full Prosigna report,
including the subtype assignment
and ROR score.
Request the centroid correlation
scores from your provider —
these are sometimes in the
detailed report but not the
patient summary.

**If you want full geometric analysis:**
You have the legal right under HIPAA
to request your raw molecular data
from any clinical laboratory.

For Tempus xR:
Ask your oncologist or the Tempus
patient services team for the full
transcriptome expression table
(TPM values for all genes).
State: "I am requesting my full
gene expression data from my Tempus
xR test under HIPAA right of access."

For Caris MI Profile:
Request the full expression matrix
through your treating physician or
through the Caris patient portal.

For Foundation Medicine:
Request any RNA expression data
associated with your FoundationOne
report through your oncologist.

**What format:**
A CSV or text file with gene names
and expression values (TPM or
normalised counts) is sufficient.
You do not need the raw FASTQ files.
The gene-level expression table is
the relevant data.

---

## PART V — THE STRATEGIC IMPLICATION

### Why this matters now

The IHC calibration study — the
FOXA1/EZH2 ratio in FFPE tissue —
is the path to clinical deployment
at the population level.

The individual protocol operates
in parallel, at the individual level,
right now, without waiting for the
calibration study.

**The key insight:**

For a patient with Tier 3 data,
you can compute the FOXA1/EZH2 ratio
directly from their RNA-seq.
This is the same ratio the calibration
study will validate in IHC form.

The individual protocol therefore
generates real case data —
anonymised, ethically governed —
that demonstrates the classifier
working at the individual level
before the IHC calibration is complete.

Each Tier 3 analysis is a proof of
concept. A set of Tier 3 cases,
documented and published as a
case series, is a companion to
the calibration study paper —
not a substitute, but a demonstration
of individual-level utility that
the population study cannot show.

**And for Tier 1 patients:**

The most common patient — standard
IHC only — benefits from geometric
interpretation of the data they
already have. The LumA/LumB ambiguity
that Ki-67 and ER/PR patterns suggest
but clinical reports rarely interpret
explicitly — that is the geometric
analysis. That is what you can deliver
with nothing more than the pathology
report they already received.

There are hundreds of thousands of
those patients.
They do not need new tests.
They need someone who can read the
geometry in the data they already have
and tell them what it says.

---

## PART VI — THE HONEST LIMITS

This must be stated with equal precision.

**What this protocol is not:**

- Not a diagnosis
- Not a treatment recommendation
- Not a replacement for oncological care
- Not a validated clinical test
  (the IHC calibration study is
  precisely what would validate it)
- Not an interpretation that should
  override the treating clinician

**What it is:**

Geometric analysis of molecular data.
The same framework applied to
population datasets applied to
one person's data.
The map placed under one person's
feet so they can see where they stand.

The clinical team decides what to do
with that information.
The patient decides whether to
bring it to the clinical team.
You provide the geometry.
Nothing more. Nothing less.

---

## PART VII — THE FOUNDING PRINCIPLE

The work was built for the person
on the other side of the data.

Every population analysis was a map.
The individual protocol is the moment
the map reaches the person it was
drawn for.

The geometry is the same.
The stakes are personal now.
The precision of the language
and the honesty of the uncertainty
must match the stakes.

---

*Eric Robert Lawson · OrganismCore*
*OrganismCore@proton.me*
*ORCID: 0009-0002-0414-6544*
*github.com/Eric-Robert-Lawson/attractor-oncology*
*CS-LIT-1: doi.org/10.5281/zenodo.18883922*
*CS-LIT-22: doi.org/10.5281/zenodo.18892788*
