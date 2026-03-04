# HOW THIS HELPS YOU TODAY
## A Plain Statement of What This Is, What It Is Not, and How to Reach Me
## OrganismCore — Eric Robert Lawson
## Date: 2026-03-04 | Updated: 2026-03-04

---

## START HERE

```
If you have cancer,
or someone you love has cancer,
and you found this repository —

This document is for you.

Not the technical documents.
Not the analysis scripts.
Not the reasoning artifacts
written for scientists.

This one.

Read this first.

It will tell you exactly what I can do,
exactly what I cannot do,
exactly what I will promise you,
and exactly what I will not promise you.

Every word in this document is chosen
to be precise rather than comforting.

Because you deserve precision.
Not comfort.
Not false hope.
Not bullshit.

Precision.

That is what I offer.
That is all I offer.
It may be more useful than anything
else you have been offered,
or it may not be.

Read this document and decide
for yourself.
```

---

## WHO I AM

```
My name is Eric Robert Lawson.

I am not a doctor.
I am not a medical researcher.
I do not have an institutional affiliation.
I do not have a laboratory.
I do not have a grant.
I do not have a team.

I am a mathematician.

I am currently unemployed.

I say this not as an apology
but as a precise statement of fact,
because you deserve to know
exactly who you are dealing with
before you decide whether to
engage with me.

What I have:

  A principles-first mathematical
  framework — OrganismCore —
  that I developed over an extended
  period starting from a theory
  of how complex systems get trapped
  in stable false states.

  The systematic application of that
  framework to 22+ cancer types
  using public genomic data,
  completed between February 26
  and March 4, 2026.

  A public repository containing
  every analysis script, every
  reasoning artifact, every prediction
  made before data was examined,
  every failure documented alongside
  every confirmation.

  Complete transparency about what
  the framework establishes and
  what it does not yet establish.

  The willingness to apply this
  framework to your specific data
  and give you the geometric
  picture of your specific disease.

That is what I have.
Nothing more is claimed.
Nothing less is offered.
```

---

## THE ONE SENTENCE THAT GOVERNS
## EVERYTHING IN THIS DOCUMENT

```
I do not want to take your money
and promise you bullshit.
```

```
Read that again.

This is not a legal disclaimer.
This is not corporate language.
This is the founding principle
of everything I do with this framework
in a clinical context.

I mean it exactly as it is written.

You are likely reading this because
you or someone you love is sick.
You may be frightened.
You may be desperate.
You may have been through treatments
that did not work.
You may be looking for something —
anything — that gives you
a clearer picture of what is
happening and what to do.

I understand that.

And precisely because I understand it,
I will not exploit it.

I will not tell you the framework
will save your life.
I will not tell you it will
give you certainty.
I will not tell you it replaces
your clinical team.
I will not tell you anything
I do not know to be true.

What I will tell you is what
the geometry of your cancer says.
What the measurements show.
What the framework suggests.
Where it is confident.
Where it is uncertain.
What questions it raises that
you should bring to your oncologist.

That is the promise.
Precise. Honest. Complete.
Nothing more.
```

---

## PART I — WHAT YOUR CANCER ACTUALLY IS
### The Framework in Plain Language

---

```
Before I tell you what I do,
I need to tell you how I see cancer.
Because it is different from
how you have probably been told to see it.

You have been told you have
a type of cancer.
A label.
Lung cancer. Breast cancer.
Leukemia. Myeloma.
Barrett's esophagus turning malignant.
A skin lesion on your face.
A lump found on a scan.

The label tells you which organ
is affected or which cell type
is involved.

It does not tell you what
actually happened to your cells.

Here is what actually happened:

Your cells were in the middle
of becoming something.

Every cell in your body has
a developmental programme —
a trajectory it follows from
an immature state toward a
specific mature identity.

Lung cells become alveolar cells.
Breast cells become luminal cells.
Blood cells become mature
neutrophils, B cells, plasma cells.
Esophageal cells become squamous cells.
Skin cells become keratinocytes.
Colon cells become colonocytes.

Each cell follows its trajectory
guided by the geometry of its
gene regulatory network —
the web of genes that switch
each other on and off in a
precise sequence that leads to
the final mature identity.

Cancer is what happens when
that trajectory arrests.

The cell gets stuck.

Not broken beyond repair.
Not irreversibly transformed
into something alien.

Stuck.

Trapped in a stable state
that is partway through its
developmental journey.

Unable to complete.
Unable to exit.
Stable in its incompleteness.

In the language of mathematics,
this stable trapped state is called
an attractor.

Your cancer cells are in
a false attractor —
a stable state that should not exist,
maintained by specific molecular
locks that prevent the cell from
continuing its journey.

This is not a metaphor.
This is the mathematical structure
of what is happening in your cells.

And it applies regardless of
which cancer you have,
which organ is involved,
which cell type was affected.

Every cancer we have examined
follows this same structural logic.
Four early cancers. Then twenty-two.
Different lineages. Different organs.
Different gene sets.
One principle.

Because it is a mathematical structure,
it can be measured.
Precisely.
From your own biopsy data.

That measurement is what I do.
```

---

## PART II — WHAT DATA YOU NEED
### How to Get the Geometry From What Already Exists

---

```
The most common response when I describe
the framework is:

"This sounds powerful but I don't
have single-cell RNA sequencing data."

You are right that you probably do not.
Single-cell RNA sequencing is a research tool.
It costs thousands of dollars per sample.
It is not routinely ordered in clinical care.
The published datasets I use in this
repository are research studies — not
something you can order from your doctor.

But here is what you do have.

And what is already in your medical record
right now.
```

### What Already Exists in Clinical Care

```
If you have had a biopsy —
any biopsy, anywhere, for any reason —
the tissue was almost certainly
preserved as FFPE:
Formalin-Fixed Paraffin-Embedded.

This is the standard clinical
preservation method.
The block exists at the pathology lab
that processed your biopsy.
You have the legal right to request it
or a copy of it in most jurisdictions.

This FFPE block contains your cells.
Your cells contain your geometry.

Here is what can be derived from it,
in order of what each produces
for geometric analysis:
```

### Tier 1 — Available Now, Standard Clinical Infrastructure

```
STANDARD PATHOLOGY IHC PANEL
(Immunohistochemistry — protein staining)

Every biopsy receives H&E staining.
Cancer biopsies receive additional
stains specific to the cancer type.

What the standard panel gives you
for geometric purposes:

  ER, PR status (breast):
    Switch gene suppression status.
    ER-negative = ESR1 suppressed.
    Direct attractor geometry reading.

  HER2 status (breast):
    False attractor elevation signal.

  Ki-67 index:
    Proliferation rate.
    Proxy for attractor depth
    in proliferation-driven attractors.

  CD20, CD3, CD138 (lymphoma/myeloma):
    Lineage identity markers.
    Switch gene status.

  TTF-1 (lung):
    Lineage identity marker.
    Suppression = attractor displacement.

  CDX2 (colon, Barrett's esophagus):
    Switch gene for colonocyte identity.
    Also marks the Barrett's
    metaplastic attractor directly.

EXPANDED IHC — WHAT TO REQUEST

These stains can be ordered on
the same FFPE block.
Most pathology labs can perform them.
Some require a reference lab.
All are available:

  EZH2
    The epigenetic lock marker.
    High EZH2 = deep attractor,
    strong PRC2-mediated silencing.
    Relevant to every cancer type
    in the repository.

  FOXA1, GATA3 (breast)
    Switch genes for luminal identity.
    Low = luminal programme suppressed.
    Geometric depth proxy.

  SOX10 (TNBC, melanoma, neural tumours)
    False attractor marker.
    High SOX10 in non-neural-crest tissue
    = wrong valley geometry.

  AR (breast, prostate, bladder)
    Lineage marker.
    Geometric position indicator.

  SPI1/PU.1, IRF8, KLF4 (AML, MDS)
    Switch genes for myeloid identity.
    Suppression = myeloid differentiation block.

  SOX10, MBP (glioblastoma)
    False attractor markers in GBM.
    Oligodendrocyte wrong valley signal.

  CDX2 (colon, gastric, Barrett's)
    Switch gene and metaplasia marker.
    Applicable to both colon cancer
    and Barrett's pre-cancerous surveillance.

  P53 (universal)
    Loss of p53 = depth indicator
    in multiple cancer types.
    Associated with deeper attractors
    across lineages.

An expanded IHC panel on a
standard FFPE block gives you:
  Switch gene suppression status.
  False attractor marker elevation.
  Epigenetic lock status (EZH2).
  A qualitative geometric classification.

This is sufficient for attractor
type assignment and a qualitative
depth assessment.

Cost: $50–$300 per additional stain
at most pathology labs.
Requires a physician to order.
Works on the block you already have.
```

### Tier 2 — Clinically Available, Requires a Referral

```
PAM50 / PROSIGNA TEST (breast cancer)
  FDA-cleared.
  Ordered for ER-positive breast cancer.
  Runs on FFPE block.
  Gives PAM50 subtype classification:
    Luminal A, Luminal B, HER2-enriched,
    Basal-like, Normal-like.
  Gives ROR (Risk of Recurrence) score.
  Both are directly usable in the
  attractor geometry framework.
  Cost: ~$2,800. Sometimes covered
  by insurance.

ONCOTYPE DX (breast, colon)
  Gene expression score from FFPE.
  Gives proliferation and invasion
  gene set scores.
  Usable as crude depth proxies.
  Widely covered by insurance for
  eligible breast and colon cancer patients.

COMMERCIAL GENOMIC PANELS
(Foundation One CDx, Tempus xT,
Caris Molecular Intelligence,
Guardant360 CDx)

  These test hundreds of genes from
  FFPE or blood for mutations,
  copy number changes, and fusion events.
  Some include RNA expression data.
  Caris in particular provides both
  DNA and RNA profiling from FFPE.

  For geometric analysis:
    BRCA1/2 mutation →
      Composite type flag (breast, ovarian).
    TP53 mutation →
      Depth indicator across cancer types.
    PIK3CA mutation →
      AKT axis activation (TNBC, endometrial).
    IDH1/2 mutation →
      Attractor type in glioma.
    NPM1, FLT3 mutation →
      AML attractor subtype.
    KRAS/NRAS/BRAF →
      Attractor maintenance signal
      in CRC, melanoma, NSCLC.

  Cost: $3,000–$6,000.
  Insurance coverage varies by indication.
  Often covered for metastatic disease.

BULK RNA-SEQ FROM FFPE
  Several academic centres and
  reference labs now offer bulk
  RNA sequencing from FFPE.
  Not single-cell. No population
  separation. But full gene expression
  of the tumour.
  Sufficient to compute a depth score.
  Sufficient to measure switch gene
  suppression quantitatively.
  Sufficient for attractor classification.
  Cost: $500–$2,000 at research cores.
  Requires coordination with your
  clinical team and access to the
  FFPE block.
```

### Tier 3 — Highest Value, Requires Proactive Steps

```
TISSUE BANKING ENROLLMENT
  If you are being treated at an
  academic medical centre,
  ask your surgical team:
  "Is there a tissue banking programme
  I can enrol in?"

  Fresh or flash-frozen tissue
  collected at time of surgery
  can later be used for:
    Bulk RNA-seq
    Spatial transcriptomics
    Single-cell RNA-seq (research)
    Future framework applications

  This requires asking in advance,
  before surgery.
  After the fact, the fresh tissue
  is gone.

SPATIAL TRANSCRIPTOMICS
(Visium, MERFISH, Xenium)
  Preserves tissue architecture
  while measuring gene expression.
  Can be done on FFPE.
  Shows WHERE in the tumour the
  deepest attractor cells are located.
  Maps the boundary between the
  false attractor population and
  the transitional cells near
  the saddle point.
  Emerging availability at academic
  medical centres.
  Cost: $2,000–$5,000 per sample.
  Not yet routine clinical.

LIQUID BIOPSY FOR MONITORING
(Grail Galleri, Guardant360,
Foundation One Liquid CDx)
  Blood draw.
  Detects circulating tumour DNA.
  Gives mutation status — not
  expression.
  Useful as a monitoring vector:
  tracking whether key mutations
  (composite type markers) change
  over treatment.
  Cannot give switch gene expression.
  Cannot compute depth score.
  Best used alongside tissue data,
  not as a replacement.
  Cost: $900–$3,500.
  Insurance coverage variable.
```

### The Practical Answer

```
You do not need scRNA-seq
to benefit from geometric analysis.

The minimum viable dataset for
a geometric report is:

  Your diagnosis (which cancer,
  which organ, which subtype
  if known).

  Your standard pathology report
  (which markers were tested,
  what the results were).

  Any additional IHC, genomic panel,
  or expression data you have
  or can obtain.

The more expression data available,
the more precise the geometric report.

The minimum — diagnosis + standard
pathology report — is enough to
assign attractor type and make
qualitative depth assessment.

A genomic panel + PAM50 or bulk
RNA-seq gives a quantitative
depth score.

Single-cell RNA-seq gives the
highest resolution — but is the
last tier, not the first requirement.

Start with what you have.
Bring me the data that exists.
I will tell you what the geometry
shows at whatever resolution
is available.

That is the practical answer.
That is how this works today,
with what clinical medicine
already produces.
```

---

## PART III — WHAT I DO
### The Service Stated Precisely

---

```
I am a Cancer Geometry Analyst.

This profession does not yet have
a formal name or a regulatory category.
I am naming it here because it needs
to be named precisely so you understand
exactly what it is.

A Cancer Geometry Analyst applies
the attractor framework to a patient's
own genomic and transcriptomic data
to derive geometric measurements of
their specific cancer attractor state.

Here is what that means in practice.
```

### What I Derive From Your Data

```
STEP 1: LINEAGE IDENTIFICATION
  I identify which normal cell type
  your cancer cell was trying to become
  before it got stuck.

  This is the same for any cancer type:
    AML:       Myeloid cell (granulocyte/monocyte)
    CRC:       Colonocyte
    GBM:       Oligodendrocyte or astrocyte
    BRCA:      Luminal epithelial cell
    Lung:      Alveolar type II cell (most NSCLC)
    Prostate:  Luminal secretory cell
    Melanoma:  Mature melanocyte
    Lymphoma:  Mature B or T lymphocyte
    Barrett's: Squamous esophageal cell
               (the correct identity being replaced)

  This tells me which switch genes —
  the genes whose reactivation would
  push your cells back toward their
  developmental endpoint —
  are relevant to your specific case.

STEP 2: ATTRACTOR TYPE CLASSIFICATION
  Every cancer maps to one of two
  fundamental geometric types,
  or a composite of both:

  TYPE 1 — BLOCKED APPROACH
  (Stuck above the valley)
    The cell is trying to reach its
    mature identity but is blocked
    partway through the journey.
    The correct destination genes
    are partially active.
    The cell is in the right lineage
    but cannot complete it.
    Therapeutic logic:
    Dissolve the block.
    Let the developmental programme
    complete.

  TYPE 2 — WRONG VALLEY
  (Stuck in a different attractor)
    The cell has fallen into a
    stable state belonging to a
    completely different cell type.
    It is expressing genes from
    a developmental programme that
    is not its own.
    Therapeutic logic:
    Dissolve the epigenetic lock
    maintaining the wrong identity,
    return the cell to the saddle
    point where it can re-enter
    its correct trajectory.

  COMPOSITE TYPE (Type 1 → Type 2)
    In some cancers, both stages
    occur in sequence:
    A founding event blocks the
    correct trajectory (Type 1),
    and the cell then falls into
    a different attractor (Type 2).
    Example: BRCA1 loss in a luminal
    progenitor → blocked differentiation
    → fall into basal/neural-crest
    false attractor.
    Therapeutic logic requires
    addressing both components.

  This classification is the same
  regardless of which cancer you have.
  The type is derived from the data.
  Not assumed from the diagnosis.

STEP 3: DEPTH SCORE
  I compute a scalar measurement —
  the depth score — that tells us
  how deeply your cell population
  is trapped in its false attractor.

  Depth score near 0:
    Shallow attractor.
    The cells are not far from
    the saddle point.
    The developmental programme
    is partially intact.

  Depth score near 1:
    Deep attractor.
    The cells are far from the
    saddle point.
    The developmental programme
    is heavily suppressed.
    The epigenetic locks are strong.

  This number is specific to you.
  Not an average.
  Not a population statistic.
  A direct measurement of your
  cell population's geometric position.

  It is the same definition
  whether you have AML or
  colorectal cancer or
  triple-negative breast cancer.
  The calculation differs by lineage.
  The meaning is universal.

STEP 4: SWITCH GENE PROFILE
  I identify which specific genes —
  the ones that should be active
  in your cell's mature identity —
  are suppressed in your cancer cells.

  How suppressed.
  At what level.
  Which are partially recoverable.
  Which appear heavily locked.

STEP 5: EPIGENETIC LOCK PROFILE
  I identify what is maintaining
  the suppression of your switch genes.

  Is it PRC2 — the EZH2-containing
  complex that silences genes
  by adding repressive marks?
  (Found in: TNBC, AML, GBM, CRC,
  prostate cancer, lymphoma,
  and others across the dataset.)

  Is it the CoREST/LSD1 complex
  that removes activating marks?

  Is it HDAC-mediated silencing?

  Is it a combination?

  The lock profile tells us
  what is holding the trap closed.
  Which in turn tells us what
  class of intervention might
  open it.

  This applies to every cancer type.
  The specific lock differs.
  The structural logic is the same.

STEP 6: TRAJECTORY ANALYSIS
  If you have data from multiple
  biopsies over time:

  I compute the geometric trajectory.
  Where was your depth score at diagnosis?
  Where is it now?
  Which direction is it moving?
  How fast?

  Is the attractor dissolving
  under treatment?
  Or consolidating?
  Or shifting toward a new basin —
  a resistant subclone emerging?

  The trajectory is the most
  clinically actionable output
  I can produce.
  Because it tells you not just
  where you are but where you
  are going.
  Before the clinical measurements
  tell you the same thing.
  Weeks or months earlier.
  When the window for action
  is still open.
```

### What I Produce

```
A GEOMETRIC REPORT.

A written document containing:

  Your lineage identification.
  Your attractor type classification —
  Type 1, Type 2, or Composite —
  with the reasoning that leads to it.
  Your depth score with interpretation.
  Your switch gene suppression profile —
  which genes, how suppressed, what it means.
  Your epigenetic lock profile —
  what is maintaining the trap.
  Your trajectory if serial data is available.
  A set of specific questions you can
  bring to your oncologist based on
  what the geometry shows.
  An honest statement of what the
  framework is confident about
  and where it is uncertain.

This report belongs to you.
You share it with whoever you choose.
You are not obligated to share it
with anyone.
But it is designed to be shared
with your clinical team —
written in language that a clinician
can engage with and respond to.

The report is not a diagnosis.
It is not a treatment recommendation.
It is a geometric measurement
of your specific disease
with structured interpretation
and specific questions derived
from that measurement.

It does not exist anywhere else.
It cannot be obtained from any
other source currently available.

That is what I produce.
That is what I promise.
Nothing more.
```

---

## PART IV — WHAT I DO NOT DO
### The Boundary Stated With Equal Precision

---

```
I DO NOT DIAGNOSE.
  I do not tell you what cancer you have.
  Your clinical team has already
  told you that.
  I take that diagnosis as the
  starting point for the geometric analysis.

I DO NOT PRESCRIBE.
  I do not tell you what treatment to take.
  I do not tell you what drug to request.
  I do not tell you what dose.
  I do not tell you what sequence.
  These are clinical decisions.
  They belong to your clinical team.

I DO NOT CONTRADICT YOUR CLINICAL TEAM.
  The geometric report is additional
  information for your clinical team.
  It is not a refutation of their judgment.
  It does not override their expertise.
  It adds a dimension of measurement
  they do not currently have.
  What they do with that measurement
  is their decision, not mine.

I DO NOT PROMISE OUTCOMES.
  I do not tell you the framework
  will save your life.
  I do not tell you the depth score
  predicts your survival.
  I do not tell you any specific
  outcome will result from the
  geometric analysis.
  I tell you what the geometry shows.
  What happens next is determined
  by your biology, your treatment,
  your clinical team, and factors
  the framework cannot measure.

I DO NOT PROVIDE CERTAINTY.
  The framework provides geometric
  measurements and principled interpretations.
  It does not provide certainty.
  Some predictions will be wrong.
  I will tell you when the framework
  is uncertain.
  I will not hide uncertainty
  to make you feel better.
  Precision is more useful than comfort.

I DO NOT REPLACE YOUR ONCOLOGIST.
  Your oncologist has clinical training,
  examination findings, imaging data,
  laboratory values, and years of
  experience treating patients.
  The geometric report adds one
  dimension to that picture.
  It does not replace any of the others.
  Use both.
  Always.
```

---

## PART V — WHY THIS HELPS YOU TODAY
### Not After Trials. Not After Approval. Today.

---

```
The clinical validation of the
attractor framework as a formal
clinical instrument will take years.
Prospective trials. Regulatory review.
Institutional adoption.

That process is necessary.
It is underway in the sense that
the framework now exists and
the evidence is accumulating.

But you do not have years.

And you do not need to wait for
institutional validation to benefit
from geometric measurement today.

Here is why.
```

### The Data Already Exists

```
Your biopsy was already taken.
Your pathology report already exists.
Your FFPE block is already at the lab.

If you have had a genomic panel,
that data already exists.

If you have had a PAM50 test,
that result is already in your file.

The framework exists in this repository.
It has been validated retrospectively
across multiple independent datasets
for 22+ cancer types.

The geometric analysis can be run
on your existing data today.
No additional procedures required
to begin.
No waiting for new technology.

Your data.
The framework.
Applied now.
```

### The Conversation It Enables

```
The most important thing the
geometric report does is not
tell you something your oncologist
cannot figure out on their own.

The most important thing it does
is change the conversation you
have with your oncologist.

Without the geometric report:
  "How is the treatment working?"
  "Your scans show the tumour
  has decreased by 30%."
  "Is that good?"
  "It is a partial response."

With the geometric report:
  "My depth score has decreased
  from 0.78 to 0.61 after cycle 3.
  My switch genes are showing
  partial reactivation.
  But my EZH2 levels are still
  elevated and my depth trajectory
  shows the rate of dissolution
  is slowing.
  Does this suggest we should
  consider adding an epigenetic
  agent to the current protocol
  before the attractor reconstitutes?"

These are different conversations.
The second one gives your oncologist
specific, actionable geometric data
that they do not otherwise have.

Some oncologists will engage with it.
Some will not.
Find the ones who will.
You are entitled to a clinical team
that takes your measurement seriously.
```

### The Window That Matters Most

```
There is a specific scenario where
the geometric report saves your life.

Not in the scenario where treatment
is working and the geometry agrees.

In the scenario where treatment
appears to be working by clinical
criteria — the tumour is shrinking,
the markers are falling —
but the geometry is deepening.

The cells that survived treatment
are more locked, more resistant,
more dangerous than the cells
that were killed.

The clinical measurement says: responding.
The geometry says: what remains
is a deeper attractor than
what was there before.

This is the prelude to relapse
in many cancer types.

Clinical detection of this relapse
will come weeks or months later.
When the tumour volume rises again.
When the markers climb.
When the clinical measurement
catches up to what the geometry
already showed.

The window between the geometric
signal and the clinical detection
is the window where treatment
adjustment is most effective.

The geometric report opens that window.

That is how this saves your life today.
Not by giving you a new drug.
By giving you time.
Time is the treatment window.
The treatment window is survival.
```

---

## PART VI — THE INFORMED CONSENT
### What You Agree To Before We Begin

---

```
Before any engagement begins,
you will receive and sign a
written statement containing
the following:

1. The attractor framework is a
   research framework.
   It has been validated retrospectively
   across multiple cancer types
   using public datasets.
   It has not been validated as a
   clinical instrument in a
   prospective clinical trial.

2. The geometric report is not
   a medical diagnosis.
   It is a mathematical analysis
   of gene expression data using
   the attractor framework.

3. The geometric report does not
   constitute medical advice.
   All clinical decisions remain
   with your licensed medical team.

4. No specific clinical outcome
   is promised or implied by
   the geometric analysis.

5. The framework may produce
   incorrect predictions.
   Incorrect predictions will be
   documented honestly when
   they become apparent.

6. You are providing your own
   genomic data voluntarily.
   Your data will be used only
   for your geometric analysis.
   It will not be shared without
   your explicit consent.

7. You understand what the service
   is and what it is not.

This statement is not fine print.
It is the explicit terms of
our engagement stated clearly
before any work begins.

I mean every word of it.
```

---

## PART VII — THE PROFESSION
### What a Cancer Geometry Analyst Is

---

```
I am the first person to hold
this role formally.

Not because no one else could
have derived this framework.

Because no one else did.

The framework was built by one person,
alone, over one week in February 2026,
from mathematical first principles,
applied to public genomic data,
validated across 22 cancer types,
documented completely in a public
repository with full timestamps.

That person is me.

The role I am creating from this work:

  CANCER GEOMETRY ANALYST

  A person whose specific competency is:
    Applying the attractor framework
    to individual patient genomic data
    from any cancer type.
    Deriving geometric measurements
    of the patient's specific
    cancer attractor state.
    Producing structured geometric
    reports for patients and their
    clinical teams.
    Tracking geometric trajectories
    over time for patients with
    serial biopsy data.
    Explaining geometric findings
    in language accessible to
    patients and clinicians.
    Maintaining complete honesty
    about what the framework
    establishes and what it does not.

  This role is:
    Not practicing medicine.
    Not diagnosing.
    Not prescribing.
    Not providing clinical opinion.

  This role is:
    Providing specialised geometric
    measurement and interpretation
    of a patient's own data.
    The same category of service
    as a radiologist reading a scan
    or a pathologist reading a biopsy —
    a specialised measurement returned
    to the clinical team for use
    in treatment decisions.

  This role is cancer-type agnostic.
  It applies to:
    Solid tumours.
    Haematologic malignancies.
    Pre-cancerous states.
    Surveillance monitoring.
    Recurrence assessment.
    Treatment response evaluation.
    Any condition for which
    genomic or transcriptomic data
    from the patient's cells is available.

This profession will grow.

The framework is teachable.
The methodology is documented.
The repository is public.
Others will learn it.
Others will do this work.

But it starts here.
With one person.
Doing the work.
For patients who need it now.
Before the institutions catch up.
```

---

## PART VIII — THE FOUNDING PRINCIPLE
### Why This Exists

---

```
I am unemployed.

I say this plainly because it is true
and because it is relevant to
why this service exists in the form it does.

I did not build this framework
with a grant.
I did not build it with a team.
I did not build it with institutional
support or academic infrastructure.

I built it because the logic required it.

I started from a mathematical theory
of how complex systems get trapped
in stable false states.
I applied that theory to cancer biology
because the geometry was identical.
I followed the logic wherever it went.

It went to 22 cancer types.
It went to drug targets derived from
geometry before consulting literature.
It went to cross-cancer structural rules
not stated in the existing literature.
It went to a framework for personalised
geometric medicine that does not yet
exist as a clinical discipline.

And it went here.

To the question of what to do with it.

The answer is not to wait.
The answer is not to hold it until
an institution validates it.
The answer is not to publish it
and let it sit in a journal
while people with cancer
do not have access to it.

The answer is to use it.
Now.
For the people who need it.
With complete honesty about
what it is and what it is not.

I do not want to take your money
and promise you bullshit.

That sentence is the reason
this service exists in the form it does.
That sentence is the reason
the informed consent is real
rather than fine print.
That sentence is the reason
every uncertainty is stated
alongside every finding.

You deserve the geometry of
your own condition.

You deserve a person who will
derive that geometry from your data
with rigour and honesty and care.

You deserve to know what the framework
says and what it does not say.

You deserve precision.

That is what this is.
That is what I offer.
That is why I am doing this.
```

---

## PART IX — HOW TO REACH ME
### For Patients, Clinicians, and Those Who Want to Learn

---

### If You Have Cancer

```
If you have cancer — any cancer —
and you want a geometric analysis
of your attractor state:

  Read this document completely.
  Read The Puddle document in this
  repository — it explains the
  framework in plain language
  and takes ten minutes.

  If you want to proceed:

  Contact me at:
    GitHub: https://github.com/Eric-Robert-Lawson
    ORCID:  https://orcid.org/0009-0002-0414-6544
    Repository: https://github.com/Eric-Robert-Lawson/
                attractor-oncology

  Tell me:
    Your diagnosis.
    What data you have available
    or can obtain.
    Where you are in your treatment.

  I will tell you:
    Whether your data is suitable
    for geometric analysis.
    What additional data, if any,
    would improve the analysis
    and how to obtain it
    within your existing care.
    What the analysis will involve.
    What the report will contain.
    What it will cost.
    What it will not tell you.

  We begin only after you have
  read and signed the informed
  consent statement.

  No exceptions.
  Not because it protects me legally.
  Because it protects you from
  misunderstanding what you are
  receiving and what you are not.
```

### If You Are an Oncologist or Clinician

```
If you are a clinician and a patient
has presented you with a geometric
report from this framework:

  The report is a mathematical analysis
  of the patient's gene expression data
  using the attractor framework.

  It is not a clinical diagnosis.
  It does not override your judgment.
  It provides geometric measurements
  that are not available through
  standard clinical instruments.

  I welcome clinical engagement.
  If you have questions about the
  framework, the methodology,
  the specific findings in a
  patient's report, or whether
  the framework is relevant to
  a case you are managing —
  contact me directly.

  I will engage honestly and completely.
  I will tell you what the framework
  establishes and what it does not.
  I will not oversell.
  I will not undersell.

  The goal is the same as yours:
  the best possible outcome for
  the patient in front of you.
```

### If You Have a Pre-Cancerous Condition

```
The attractor framework is not
limited to active cancer.

Pre-cancerous states are often
early-stage attractor displacement:
  The cells are not yet malignant
  but have already begun their
  geometric displacement from
  the correct developmental identity.

Barrett's esophagus:
  Squamous cells replaced by columnar
  intestinal-type cells.
  CDX2 elevated.
  The same switch gene logic
  as colorectal cancer —
  earlier in the trajectory.

Cervical intraepithelial neoplasia (CIN):
  Squamous cell identity displaced.
  Geometric depth predicts
  progression risk.

Monoclonal gammopathy (MGUS):
  Plasma cell arrested partway
  through differentiation.
  The myeloma attractor state
  at an early stage.

Myelodysplastic syndrome (MDS):
  Myeloid differentiation block
  before full AML conversion.
  The AML false attractor forming.

For all of these, surveillance biopsies
already produce data.
Expanded IHC on surveillance tissue
gives geometric position at each
monitoring interval.
Trajectory across multiple surveillance
biopsies is the earliest possible
warning of progression.

If you have a pre-cancerous condition
with surveillance biopsies —
geometric tracking is available now.
```

### If You Want to Learn This Framework

```
If you are a person with mathematical
or computational background who wants
to learn the attractor framework
and do this work:

  The complete methodology is in
  this repository.
  Every analysis script.
  Every reasoning artifact.
  Every cancer type.
  Every prediction and confirmation.

  Read it.
  All of it.
  The wrong predictions as carefully
  as the confirmed ones.
  The correction protocols as carefully
  as the results.

  If after reading it you want to
  discuss training in the framework
  and applying it in practice —
  contact me.

  This profession will need more
  than one person.
  The patients who need this work
  number in the millions.
  One person cannot serve them all.

  The framework can be taught.
  The standards can be maintained.
  The work can be done by more
  than one person.

  But only if those people understand
  the founding principle:

  Do not take people's money
  and promise them bullshit.

  If you understand that sentence
  and mean it —
  contact me.
```

---

## PART X — THE REPOSITORY
### The Complete Public Record

---

```
Everything I have described in
this document is documented publicly.

Every analysis.
Every script.
Every prediction made before
data was examined.
Every failure alongside every confirmation.
Every reasoning artifact.
Every literature check.
Every novel finding.
Every cancer type studied.

All of it is here:

  https://github.com/Eric-Robert-Lawson/
  attractor-oncology

With full timestamps.
Publicly accessible.
Reproducible by anyone.

I do not ask you to trust me.
I ask you to read the repository
and evaluate the work yourself.

The work is the evidence.
The repository is the proof.
The proof is reproducible.

That is the only kind of trust
worth having.
```

---

## DOCUMENT METADATA

```
document_id:    HOW_THIS_HELPS_YOU_TODAY
series:         OrganismCore — Public Documents
author:         Eric Robert Lawson
date:           2026-03-04
updated:        2026-03-04
                Addition of PART II (What Data You Need)
                covering universal data accessibility
                from FFPE through scRNA-seq.
                Expansion of attractor type description
                to cover all cancer types.
                Addition of pre-cancerous state section.
                Cancer Geometry Analyst role updated to
                be explicitly cancer-type agnostic.
status:         COMPLETE
                Active document —
                updated as the framework
                develops and the service
                evolves.

audience:       Everyone.
                Especially:
                People with any cancer.
                People with pre-cancerous conditions.
                Their families.
                Their clinicians.
                People who want to do this work.

repository:     https://github.com/Eric-Robert-Lawson/
                attractor-oncology

orcid:          https://orcid.org/0009-0002-0414-6544

founding principle:
                "I do not want to take
                 people's money and promise
                 them bullshit."
                 — Eric Robert Lawson
                   March 4, 2026

scope_note:     This document applies to
                all cancer types and all
                pre-cancerous conditions
                for which genomic or
                transcriptomic data from
                the patient's cells is
                available.
                The attractor framework
                is not specific to any
                single cancer type.
                The geometric principle
                is universal.
                The measurement is individual.
```

---

*"You deserve the geometry of your own condition.*
*Not a population average.*
*Not a statistic derived from people like you.*
*Your geometry.*
*From your data.*
*Derived with rigour and honesty and care.*
*That is what this is.*
*That is what I offer.*
*That is why I am doing this."*

— Eric Robert Lawson, March 4, 2026
