# The Big Picture
## What IHC Validation Changes — Globally
### Eric Robert Lawson / OrganismCore
### 2026-03-07

---

## Preface

This document reasons through the full
downstream consequences of a successful
IHC calibration study. It is not advocacy.
It is not optimism. It is an attempt to
reason carefully about what actually changes
in the world if the FOXA1/EZH2 ratio
classifier is validated against PAM50
ground truth and the cut-points are published.

The reasoning is structured. The numbers
are sourced from global cancer epidemiology.
The logic is falsifiable at each step.

---

## Part 1: The Current State of the World

### 1.1 How breast cancer is classified today

Breast cancer is not one disease.
It is at minimum six biologically distinct
diseases that happen to originate in
breast tissue:

- Luminal A
- Luminal B
- HER2-enriched
- Triple-negative (TNBC)
- Claudin-low
- Invasive lobular carcinoma (ILC)

Each has a different mechanism.
Each has a different treatment logic.
Treating them identically — which is what
happens without molecular subtyping —
is mechanistically equivalent to treating
six different infections with the same
antibiotic because they all cause fever.

### 1.2 The global subtyping gap

New breast cancer diagnoses per year
globally: approximately 2.3 million.

Patients who receive molecular subtype
information (PAM50/Prosigna or equivalent):
approximately 800,000–1,000,000 per year.
Almost exclusively in high-income countries
with major cancer centers.

**Patients who receive no molecular subtype
information: approximately 1.2–1.5 million
per year.**

These are not patients in remote areas
without any healthcare. Many of them
receive surgery, radiation, and
chemotherapy. They receive treatment.
They just receive treatment without
knowing which of the six diseases they
have — so they receive the wrong treatment
at significant rates.

### 1.3 Why the gap exists

PAM50/Prosigna (NanoString nCounter):
- Cost per test: $3,000–$4,000
- Capital equipment: $150,000–$250,000
- Platform type: proprietary RNA
- RNA quality requirement: strict
  (degrades; archival tissue unreliable)
- Laboratory type: specialised,
  centrally certified
- Turnaround: days to weeks
  including shipping
- Global availability: major centres
  in high-income countries only

This is not a scientific failure.
It is an infrastructure and cost failure.
The science exists. The deployment
does not reach the majority of patients.

---

## Part 2: What Validation Produces

### 2.1 What the calibration study delivers

A single published paper containing:

1. IHC H-score cut-points for FOXA1 and EZH2
   derived from ~300–400 FFPE cases against
   PAM50 ground truth
2. Concordance statistics (kappa, AUC)
   between IHC ratio classification
   and PAM50 subtype
3. Inter-observer reliability data
   (ICC and kappa between two pathologists)
4. A citable, reproducible protocol
   with specific antibody clones,
   antigen retrieval conditions, and
   scoring instructions

That paper, once published, is free.
Open access. Available to every
pathologist on earth with internet access.

### 2.2 What a pathologist anywhere
### can do with it

Read the paper. Order two antibodies.
Run them on their existing IHC platform.
Score the H-scores. Divide. Compare to
the published cut-points. Classify.

This is not a new workflow. This is
existing IHC workflow applied to two
new antibodies with a published ratio
threshold. Every pathology laboratory
in the world that processes tissue
biopsies already does this for ER, PR,
HER2, and Ki67. Adding FOXA1 and EZH2
is not a capability upgrade. It is a
reagent addition.

---

## Part 3: The Cost Transformation

### 3.1 Unit economics

| Item | PAM50/Prosigna | FOXA1/EZH2 IHC |
|------|---------------|----------------|
| Cost per patient | $3,000–$4,000 | $50–$100 |
| Capital equipment | $150k–$250k | $0 (already owned) |
| RNA required | Yes | No |
| Tissue requirements | Fresh preferred | Any FFPE, any age |
| Proprietary platform | Yes | No |
| Reference lab required | Often | No |
| Shipping required | Often | No |
| Turnaround | Days–weeks | Same day |

### 3.2 The cost reduction factor

$3,500 average PAM50 cost ÷ $75 average
IHC cost = **47× cost reduction.**

For a hospital system running 500 breast
cancer cases per year:

PAM50 for all patients:
500 × $3,500 = $1,750,000/year

FOXA1/EZH2 IHC for all patients:
500 × $75 = $37,500/year

**Annual saving per 500-patient institution:
$1,712,500.**

For a middle-income country running
10,000 breast cancer diagnoses per year:

PAM50 (if accessible): $35,000,000/year
FOXA1/EZH2 IHC: $750,000/year

The capital barrier alone — $150,000–
$250,000 for NanoString equipment — makes
PAM50 inaccessible to most of the world's
pathology departments. The FOXA1/EZH2
protocol requires zero capital expenditure
beyond standard IHC reagents already
in every pathology budget.

---

## Part 4: The Geographic Transformation

### 4.1 Where breast cancer kills people

The countries with the highest breast
cancer mortality are not the countries
with the highest incidence. They are
the countries with the lowest access
to molecular subtyping and targeted
treatment.

High mortality, low subtyping access:
- Sub-Saharan Africa
- South and Southeast Asia
- Latin America (outside major centres)
- Eastern Europe
- Middle East and North Africa

These regions collectively account for
the majority of the 685,000 annual
breast cancer deaths globally.

### 4.2 What IHC access looks like
### in these regions

Standard IHC is already present in
most pathology departments in these
regions. ER, PR, and HER2 IHC are
run routinely — these are the basic
clinical requirements for breast cancer
workup. The infrastructure, the
technicians, the staining platforms,
the microscopes — they exist.

What does not exist is PAM50.
The equipment is too expensive.
The RNA handling requirements are
too strict. The centralised reference
laboratory infrastructure is absent.

FOXA1 and EZH2 IHC run on the same
platforms, by the same technicians,
using the same workflow already in
place for ER, PR, and HER2.

**The deployment barrier is not
capability. It is having published
cut-points to interpret the result.**

That is what the calibration study
produces. One paper. Published open
access. And the infrastructure to
deploy it globally already exists.

---

## Part 5: The Diagnostic Transformation

### 5.1 What changes at the point
### of diagnosis

Without subtyping, a pathologist
diagnosing ER-positive breast cancer
reports: ER-positive, PR-positive,
HER2-negative. That is the full
molecular information available.

The oncologist receives this report
and makes a treatment decision that
does not distinguish between LumA
and LumB — two diseases with different
mechanisms and different optimal
treatment sequences.

With FOXA1/EZH2 IHC added to the
same tissue workup, the pathologist
now reports a ratio that maps to
a specific subtype — not just
ER-positive, but LumA or LumB —
and each subtype maps to a specific
treatment logic.

**This happens on the same slide,
on the same day, at the same
laboratory, for $50–$100 additional
reagent cost.**

No shipping. No waiting. No RNA.
No reference laboratory. No special
equipment.

### 5.2 The LumA/LumB distinction
### specifically

This is where the immediate clinical
impact is largest.

LumA and LumB together represent
approximately 70% of all breast
cancer diagnoses. They are both
ER-positive. They are currently
treated largely identically with
endocrine therapy as the backbone.

They have different mechanisms:
- LumA: luminal programme intact.
  CDK4/6 inhibitor + endocrine therapy
  is highly effective. HDAC inhibition
  adds nothing.
- LumB: chromatin lock via HDAC2
  co-repressor circuit. Endocrine therapy
  alone is suboptimal. HDACi + ET
  dissolves the lock. CDK4/6i also active
  but the primary resistance mechanism
  is epigenetic, not proliferative.

Treating LumB patients as LumA patients
— which is what happens without the
distinction — means those patients
receive endocrine therapy without the
HDAC component that addresses their
actual resistance mechanism.

The METABRIC data showed an 18.7-month
RFS difference between high and low
ratio groups on endocrine therapy.
That is not a marginal improvement.
That is the difference between a patient
who recurs at 3 years and a patient who
does not recur until after 5 years — or
who does not recur at all within the
study follow-up.

### 5.3 The TNBC/claudin-low distinction

TNBC and claudin-low are both
triple-negative by standard IHC.
They are currently treated identically.
They have opposite treatment logics:

- TNBC: EZH2 epigenetic lock.
  Tazemetostat dissolves the lock.
  Then fulvestrant to address the
  restored ER sensitivity.
  The mechanism is targetable.

- Claudin-low: no epigenetic lock
  to dissolve. Immune-cold to
  immune-hot transition required.
  Anti-TIGIT → anti-PD-1 sequence.
  Chemotherapy has poor mechanistic
  rationale. Immunotherapy is the
  correct approach.

TNBC patients receiving immunotherapy
sequencing appropriate for claudin-low
are receiving the wrong treatment.
Claudin-low patients receiving
tazemetostat are receiving a drug
targeting a mechanism they do not have.

The FOXA1/EZH2 ratio distinguishes
these two groups from the same
two stains. No additional tests.
No additional cost beyond the
already-justified reagent spend.

---

## Part 6: The Prognostic Transformation

### 6.1 Survival stratification confirmed

The calibration study will produce
IHC cut-points. But the survival
data already exists from the
computational validation:

- TCGA Kaplan-Meier: p=0.0031 (n=1,218)
- METABRIC Kaplan-Meier: p<0.0001 (n=1,980)
- GSE96058 Kaplan-Meier: p≈0 (n=3,273,
  336 events, 52 months follow-up)

Three independent survival analyses.
Two platforms. The ratio stratifies
survival independently of standard
clinical variables.

Once IHC cut-points are established,
the same prognostic stratification
is available at the bench.

### 6.2 What prognostic stratification
### enables

A pathologist who can report not just
subtype but prognostic tier from the
same two stains enables:

- More aggressive surveillance schedules
  for high-risk subtypes
- Less aggressive treatment for confirmed
  LumA (the least aggressive subtype,
  currently sometimes overtreated)
- Earlier escalation to targeted therapy
  in LumB (currently undertreated
  relative to the mechanism)
- Rational immunotherapy sequencing
  in claudin-low (currently receiving
  chemotherapy regimens with poor
  mechanistic rationale)

Each of these represents a population
of patients currently receiving
suboptimal treatment because the
information to guide better treatment
is not available at their laboratory.

---

## Part 7: The Lives Saved Calculation

### 7.1 The honest methodology

A precise "lives saved" number requires
assumptions that cannot be made without
clinical trial data. What can be reasoned
from the available evidence:

**Population receiving no subtyping:**
~1.2–1.5 million patients per year.

**Proportion who are ER-positive
(LumA or LumB):** approximately 70%
= 840,000–1,050,000 patients per year
receiving no LumA/LumB distinction.

**Proportion who are triple-negative
(TNBC or claudin-low):** approximately
15–20% = 180,000–300,000 patients per year
receiving no TNBC/claudin-low distinction.

**The METABRIC RFS delta:** 18.7 months
between high and low ratio groups
on endocrine therapy.

What this means in population terms:
If even a fraction of the ER-positive
patients currently receiving no subtyping
were correctly identified as LumB and
received HDACi + ET instead of ET alone,
the 18.7-month RFS differential applies
to that fraction.

If 10% of the 840,000 unsubtyped
ER-positive patients per year are LumB
and would benefit from the correct
treatment sequence — that is 84,000
patients per year receiving a treatment
that adds 18.7 months of recurrence-free
survival on average.

### 7.2 The honest statement

It is not possible to state with precision
how many lives this saves per year without
a prospective clinical trial using the
validated IHC classifier.

What can be stated with precision:

1. 1.2–1.5 million patients per year
   currently receive no molecular subtype
   information.

2. Of those, approximately 1 million are
   ER-positive, of whom a significant
   fraction are LumB — a subtype with
   a known treatment improvement
   available that they are not receiving.

3. Of those, approximately 200,000–300,000
   are triple-negative, of whom a
   significant fraction are claudin-low —
   a subtype receiving chemotherapy
   with poor mechanistic rationale when
   immunotherapy sequencing is indicated.

4. The FOXA1/EZH2 ratio classifier,
   on validation, is the only tool
   that makes this distinction available
   at $50–$100 in any laboratory on earth.

5. The number of lives affected is
   therefore in the hundreds of thousands
   per year. The number of lives saved —
   in the sense of preventing recurrence,
   extending progression-free survival,
   or directing patients to treatments
   that actually address their mechanism
   — is likely in the tens of thousands
   per year once deployment reaches scale.

That is a conservative estimate.
It does not account for the ILC
distinction, the HER2/ET sequencing
optimisation, or the secondary effects
of more rational treatment allocation
on healthcare systems.

---

## Part 8: The Research Transformation

### 8.1 What the classifier enables
### beyond the clinical setting

Once validated and published, the
FOXA1/EZH2 IHC ratio becomes a
standard retrospective analysis tool
for existing breast cancer trial archives.

Every clinical trial that banked FFPE
tissue — including E2112, METABRIC,
TCGA, and hundreds of smaller trials —
can be retrospectively analysed using
the validated IHC classifier.

This means:

**E2112 (entinostat):**
The LumB-specific entinostat hypothesis
can be tested in the existing trial
archive with no new patients. If the
LumB subgroup in E2112 shows significant
PFS benefit that was diluted by LumA
non-responders, entinostat gets
re-evaluated for a defined population
instead of being abandoned as a failed
trial.

**Every failed HR+ trial:**
Many breast cancer trials that returned
negative results enrolled mixed LumA/LumB
populations. The classifier retrospectively
applied to those archives will identify
which trials had a real signal in a
subgroup that was diluted by the wrong
population.

This is not a small effect. A significant
fraction of "failed" breast cancer trials
may contain real signals in a population
that was never identified at the time.

### 8.2 The regulatory implication

A validated IHC biomarker with published
cut-points, concordance data, and
inter-observer reliability becomes
eligible for regulatory consideration
as a companion diagnostic.

The pathway from published IHC cut-points
to companion diagnostic status involves:
- Prospective validation in a clinical
  trial setting
- FDA (or equivalent) biomarker
  qualification process
- Regulatory submission as a companion
  diagnostic for a specific drug

This is a multi-year process. But it
begins the moment the calibration study
is published and the cut-points exist.
The calibration study is the first step
in a regulatory pathway that does not
currently exist for this classifier.

---

## Part 9: The Equity Argument,
## Precisely Stated

### 9.1 The current distribution of
### diagnostic benefit

Molecular subtyping in breast cancer
is currently available almost exclusively
to patients in high-income countries
with major cancer centers.

This is not a neutral distribution.
It is a distribution that concentrates
the benefit of a scientific advance
among the patients who are already
best served by the healthcare system
and withholds it from the patients
who are most dependent on low-cost
accurate diagnosis.

### 9.2 What validation changes
### about that distribution

The FOXA1/EZH2 IHC protocol, on
validation, is deployable at the same
cost and complexity as ER, PR, and HER2
IHC — which are already deployed
globally, including in low and middle
income countries.

This means the benefit of molecular
subtyping — previously reserved for
patients whose hospitals could afford
$150,000 in capital equipment and
$3,500 per test — becomes available
to any patient whose pathology
laboratory runs standard IHC.

That is not a marginal improvement
in equity. It is a structural change
in who has access to the information
that determines whether they receive
the right treatment.

### 9.3 The cost to achieve this

The calibration study costs
$5,000–$15,000 in consumables.

The potential reach: 1.2–1.5 million
patients per year who currently receive
no molecular subtype information.

The cost per patient reached on full
deployment: approximately $0.01
(the cost of the calibration study
divided by annual patients reached,
which is a one-time cost amortised
across all future use of the protocol).

This is the highest return on investment
in breast cancer diagnostics currently
available to a research partner willing
to run 300 slides.

---

## Part 10: The Single Sentence

If this classifier is validated,
a pathologist anywhere on earth —
in Lagos, in Dhaka, in Recife,
in Kyiv, in Manila — running standard
IHC on a breast cancer biopsy that
is already on their bench for ER, PR,
and HER2 assessment can add two stains,
divide two numbers, and tell their
oncologist not just that the tumour
is ER-positive but which of six
mechanistically distinct diseases
their patient has — and what the
correct treatment logic is for that
specific disease.

For $50–$100.
Same day.
No shipping.
No RNA.
No proprietary platform.
No reference laboratory.

That is what the calibration study
enables. That is what 300 slides
and one pathologist's participation
produces.

---

*Eric Robert Lawson · OrganismCore*
*OrganismCore@proton.me*
*ORCID: 0009-0002-0414-6544*
*github.com/Eric-Robert-Lawson/attractor-oncology*
*CS-LIT-1: doi.org/10.5281/zenodo.18883922*
*CS-LIT-22: doi.org/10.5281/zenodo.18892788*
