# INDIVIDUAL PATIENT GEOMETRIC ANALYSIS PROTOCOL
## The Complete Methodology for Single-Patient Attractor Geometry Scoring
## OrganismCore — Eric Robert Lawson
## Date: 2026-03-04

---

## PREAMBLE

```
This document exists because a question
was asked precisely:

  "Suppose someone gives me data
   from their biopsy.
   How do I utilize it properly?"

The answer requires a complete protocol.

Not a summary.
Not a guide.
A protocol.

The same rigour that governs the
population analyses in this repository —
predictions before data, full output
preserved, honest uncertainty documented —
governs the individual patient analysis.

The difference is this:

The population analyses establish
the reference geometry.
The individual patient analysis
locates one person within it.

Both require the same epistemic discipline.
Both require the same honest accounting
of what is established and what is not.
Both require complete documentation.

This document specifies the complete
methodology for individual patient
geometric analysis from first receipt
of patient data to delivery of the
geometric report.

It is written so that:
  Any competent analyst trained in
  this framework can follow it
  without asking how.

  Any patient who receives a report
  from this protocol can understand
  exactly how it was produced.

  Any clinician who receives a report
  can evaluate the methodology that
  produced it.

  Any future analyst learning this
  work can reproduce it exactly.

The protocol is the ethics made operational.
The ethics are stated in
HOW_THIS_HELPS_YOU_TODAY.md.
The protocol is how those ethics
are preserved in practice.

Read this document completely
before accepting any patient data.
```

---

## PART I — THE TWO PROBLEMS
### Understanding What This Protocol Solves

---

```
THE POPULATION PROBLEM:
  Given a cancer type,
  what is the attractor geometry
  of that cancer in general?
  What are the switch genes?
  What is the depth axis?
  What are the drug targets?

  Input: GEO dataset with hundreds
         of tumour and normal samples.
  Output: Reference attractor geometry
          for that lineage.
  Protocol: Workflow_Protocol.md
  Status: Complete for 22+ cancer types.

─────────────────────────────────────────

THE INDIVIDUAL PROBLEM:
  Given one patient's biopsy data,
  where does THIS PATIENT sit in the
  attractor geometry already established?

  Input: One patient's expression data
         from their biopsy.
  Reference: The established attractor
             geometry for their cancer type
             from the population analysis.
  Output: That patient's geometric report.
  Protocol: This document.
  Status: The protocol being established here.

─────────────────────────────────────────

WHY THEY ARE DIFFERENT:

  The population analysis uses
  statistics across many samples.
  Mann-Whitney U tests.
  Spearman correlations.
  Mean expressions compared.
  Population-level patterns extracted.

  The individual analysis uses
  the population geometry as a
  coordinate system.
  One patient's values are placed
  into that coordinate system.
  Not compared to other patients.
  Measured against the reference geometry.

  This distinction is the difference
  between cartography and navigation.

  The population analysis draws the map.
  The individual analysis says:
  here is where you are on the map.

  Both are necessary.
  Neither replaces the other.
  The map must exist before
  navigation is possible.
  This is why the population analyses
  come first.
```

---

## PART II — PREREQUISITES
### What Must Exist Before Any Patient Analysis Begins

---

```
PREREQUISITE 1:
  The reference geometry for the
  patient's cancer type must already
  exist in the repository.

  This means:
    A complete population analysis
    has been run for their cancer type
    following Workflow_Protocol.md.
    Documents [N]a, [N]b, and [N]c
    are complete for that cancer type.
    The switch gene panel is established.
    The FA marker panel is established.
    The depth score formula is established.
    The reference normal expression
    means are recorded.
    The reference cancer expression
    distributions are recorded.
    The drug map is derived and
    literature-validated.

  If the reference geometry does not
  exist for a patient's cancer type:
    DO NOT PROCEED with individual
    patient analysis.
    Run the population analysis first.
    Follow Workflow_Protocol.md
    Phase 0 through Phase 5.
    Only then return to this protocol.

PREREQUISITE 2:
  The patient's informed consent
  must be obtained before any
  data is received.

  The informed consent document
  (specified in Part IV of this protocol)
  must be:
    Read by the patient completely.
    Signed by the patient.
    Dated.
    Filed before any data transfer.

  No data is accepted without
  completed informed consent.
  No exceptions.

PREREQUISITE 3:
  The analyst must have read and
  understood:
    HOW_THIS_HELPS_YOU_TODAY.md
    Patient_Geometric_Sovereignty_
    Reasoning_Artifact.md
    This document completely.

  The ethics and the methodology
  are inseparable.
  Knowing the methodology without
  understanding the ethics produces
  reports that may cause harm.
  Knowing the ethics without the
  methodology produces intentions
  without deliverables.
  Both are required.
```

---

## PART III — WHAT DATA TO REQUEST
### Precisely What the Patient Needs to Obtain

---

```
THE CORE REQUEST:
  The patient's gene expression data
  from their tumour biopsy.
  This data already exists if they
  have had a biopsy as part of
  standard clinical care.
  It was generated from their cells.
  It belongs to them.
  They have the right to request it.

HOW TO REQUEST IT:
  The patient contacts their oncologist
  or the pathology department where
  their biopsy was processed and asks:

  "I would like a copy of the gene
  expression data from my biopsy.
  Specifically I am requesting the
  RNA-seq data in a standard format
  such as a raw counts matrix or
  normalized expression values
  (CPM or TPM) with gene symbols.
  I understand this data belongs to
  me and I am requesting it for
  my own use."

  In most jurisdictions this right
  is legally protected.
  If the institution resists:
    Request it in writing formally.
    Reference patient data rights
    under applicable law.
    Persist.

DATA FORMATS — RANKED BY QUALITY:

  TIER 1 — IDEAL:
    Bulk RNA-seq from tumour biopsy.
    Format: raw counts matrix or
            normalized CPM/TPM table.
    Gene identifiers: HGNC symbols preferred.
            Ensembl IDs acceptable.
    Coverage: whole transcriptome
              (10,000+ genes).
    Source: clinical RNA-seq run on
            the biopsy tissue.

  TIER 2 — EXCELLENT:
    Single-cell RNA-seq from biopsy.
    Format: counts matrix (genes × cells).
    The tumour cell population can be
    extracted and averaged.
    Provides richer resolution than bulk.
    More computationally intensive.

  TIER 3 — ACCEPTABLE:
    Targeted gene expression panel.
    Format: expression values for
            a defined gene panel.
    Requires: the panel must include
    the switch genes, FA markers, and
    epigenetic lock genes for the
    patient's cancer type.
    Coverage: limited but sufficient
    if the right genes are present.
    Sources: NanoString, targeted
    RNA panels, OncoScan, etc.

  TIER 4 — MINIMUM:
    IHC staining values for the
    3-gene clinical panel.
    Each cancer analysis produces
    a 3-gene panel approximating
    the depth score at r > 0.85.
    If all that exists is three
    IHC H-scores or staining intensities —
    an approximate depth score
    is computable.
    Less precise. Directionally valid.
    Honest about the approximation.

  NOT USABLE:
    DNA mutation data alone.
    Copy number variation alone.
    Imaging data.
    Blood marker values without
    expression data.
    Clinical notes without
    molecular measurements.

    The framework is geometric.
    It requires expression measurements.
    Geometry cannot be derived from
    mutation profiles alone.

SERIAL BIOPSY DATA:
  If the patient has had multiple
  biopsies over time — at diagnosis,
  after treatment, at recurrence —
  ALL available expression datasets
  are requested.

  Serial data enables:
    Trajectory computation.
    Treatment response geometry.
    Early relapse detection.
    The most clinically valuable
    output the framework produces.

  Even two time points provide
  a direction vector.
  Direction is more informative
  than position alone.
```

---

## PART IV — INFORMED CONSENT
### The Complete Document

---

```
PATIENT GEOMETRIC ANALYSIS
INFORMED CONSENT DOCUMENT

Patient ID (anonymous):  ________________
Cancer Type:             ________________
Date:                    ________________
Analyst:                 Eric Robert Lawson

─────────────────────────────────────────
WHAT YOU ARE CONSENTING TO

You are consenting to provide your
gene expression data from your tumour
biopsy to Eric Robert Lawson for
the purpose of computing a geometric
analysis of your cancer attractor state
using the OrganismCore attractor framework.

─────────────────────────────────────────
WHAT THE ANALYSIS IS

The attractor framework applies a
mathematical model of cancer — derived
from the Waddington epigenetic landscape —
to your gene expression data.

It computes:
  A depth score — measuring how deeply
  your cells are trapped in a false
  attractor state.
  A switch gene profile — identifying
  which developmental genes are suppressed
  in your cancer cells.
  An epigenetic lock profile — identifying
  what is maintaining those suppressions.
  An attractor classification — whether
  your cancer is a differentiation
  attractor or a survival attractor.
  A drug map position — where your
  geometry falls on the treatment
  landscape derived from the framework.

─────────────────────────────────────────
WHAT THE ANALYSIS IS NOT

The analysis is NOT a medical diagnosis.
The analysis is NOT a treatment
recommendation or prescription.
The analysis is NOT a clinical prognosis.
The analysis is NOT a replacement for
your medical team's judgment.
The analysis is NOT validated as a
clinical instrument in a prospective
clinical trial.

─────────────────────────────────────────
THE FRAMEWORK STATUS

The OrganismCore attractor framework
has been validated retrospectively
across 22+ cancer types using public
gene expression datasets. In every
validated cancer type, the framework
correctly identified drug targets
confirmed by published pharmacology
and active clinical trials.

The framework has NOT been validated
in a prospective clinical trial.
It has NOT been approved as a clinical
diagnostic instrument by any regulatory
body.
It is a research framework with strong
retrospective validation and no
prospective clinical validation.

This status is stated honestly and
completely. It will not be overstated.

─────────────────────────────────────────
WHAT IS PROMISED

A geometric report will be produced
containing:
  Your depth score with interpretation.
  Your switch gene profile.
  Your epigenetic lock profile.
  Your attractor classification.
  Your drug map position.
  Specific questions you can bring
  to your clinical team.
  An honest statement of what the
  framework is confident about and
  where it is uncertain.

The report will be completed within
the timeline agreed at engagement.
The report will be explained to you
in language you can understand.
Questions about the report will be
answered honestly and completely.

─────────────────────────────────────────
WHAT IS NOT PROMISED

No specific clinical outcome is
promised or implied by this analysis.
No claim is made that the analysis
will save your life, extend your
survival, or improve your response
to treatment.
No specific treatment is recommended.

─────────────────────────────────────────
YOUR DATA

Your gene expression data will be:
  Used exclusively for your geometric
  analysis.
  Stored securely and identified only
  by an anonymous patient ID.
  Not shared with any third party
  without your explicit written consent.
  Not published in any form without
  your explicit written consent.
  Deleted upon your request at any time.

─────────────────────────────────────────
YOUR RIGHT TO WITHDRAW

You may withdraw from this analysis
at any time before the report is
delivered. Your data will be deleted
immediately upon request.

─────────────────────────────────────────
CONFIRMATION

By signing below I confirm that:

  I have read this document completely.
  I understand what the analysis is
  and what it is not.
  I understand the validation status
  of the framework.
  I understand that this is not
  medical advice.
  I am providing my data voluntarily.
  I have had the opportunity to ask
  questions and they have been answered.

Patient signature:  ________________
Date:              ________________

Analyst signature:  ________________
Date:              ________________
─────────────────────────────────────────
```

---

## PART V — DATA RECEIPT AND VERIFICATION
### What Happens When the Data Arrives

---

```
STEP 1: ACKNOWLEDGE RECEIPT
  Confirm receipt of data to the patient.
  Record:
    Date received.
    Format received.
    File names.
    Data type (RNA-seq / panel / IHC).
    Time points if serial.

STEP 2: VERIFY CANCER TYPE MATCH
  Confirm the patient's cancer type
  matches an existing reference geometry
  in the repository.

  If it does not match:
    Stop immediately.
    Notify the patient.
    Explain that the reference geometry
    for their cancer type must be
    established first.
    Run Workflow_Protocol.md Phase 0-5
    for their cancer type.
    Only then proceed.

STEP 3: VERIFY DATA QUALITY

  CHECK 3.1 — FORMAT VERIFICATION:
    Is the data in a readable format?
    Gene expression values present?
    Gene identifiers present?
    If IDs are Ensembl — map to symbols
    using the Ensembl map in the scripts.

  CHECK 3.2 — COVERAGE VERIFICATION:
    For each required gene in the
    reference geometry panel — check:
    Is this gene present in the data?

    REQUIRED GENES FOR ANY CANCER:
      All switch genes from [N]b
      All FA markers from [N]b
      All epigenetic lock genes:
        EZH2, KDM1A, HDAC1, HDAC2,
        DNMT3A, TET2, ASXL1
      All 3-gene clinical panel genes

    COVERAGE TIERS:
      Full coverage:
        All required genes present.
        Full geometric report possible.
        Proceed.

      Partial coverage (>70% present):
        Most required genes present.
        Geometric report possible with
        noted gaps.
        State which genes are missing
        and how this limits the analysis.
        Proceed with honest limitation statement.

      Minimal coverage (50-70% present):
        Significant gaps.
        Approximate report only.
        Explicitly state the approximation.
        Depth score will be less precise.
        Proceed only if patient accepts
        the limitation in writing.

      Insufficient coverage (<50%):
        Cannot produce a meaningful
        geometric report.
        Stop.
        Notify patient.
        Explain what data is needed.
        Do not produce a report from
        insufficient data.
        A wrong report is more harmful
        than no report.

  CHECK 3.3 — NORMALISATION VERIFICATION:
    Are the values in a standard
    normalised format?
      TPM: sum across genes ~1,000,000
      CPM: sum across genes ~1,000,000
      Log2 CPM: values typically 0-16
      Raw counts: integers, high variance
    Record the format.
    Apply appropriate normalisation
    before scoring.

  CHECK 3.4 — SAMPLE IDENTITY VERIFICATION:
    For bulk RNA-seq: one column of values.
    For scRNA-seq: matrix with cell barcodes.
    For serial data: multiple columns
    labelled by time point or biopsy date.
    Confirm the structure before analysis.

STEP 4: RECORD THE DATA RECEIPT
  Create a case record file:

    CASE_[anonymous_ID]/
      intake_record.md
        Patient ID (anonymous)
        Cancer type
        Date of consent
        Date of data receipt
        Data format
        Gene coverage check result
        Normalisation format
        Time points (if serial)
        Quality tier (Full / Partial /
                      Minimal / Insufficient)
      raw_data/
        [patient expression files —
         as received, unmodified]
      analysis/
        [all analysis outputs]
      report/
        [final geometric report]
```

---

## PART VI — THE REFERENCE GEOMETRY EXTRACTION
### Pulling What Is Needed From the Population Analysis

---

```
Before running the patient scoring,
extract the reference parameters
from the relevant cancer analysis
in the repository.

These parameters are already in the
Documents [N]b for that cancer type.
They need to be formalised for
programmatic use.

WHAT TO EXTRACT:

FROM DOCUMENT [N]b — SCRIPT 2:

  SWITCH GENE PANEL:
    Gene names
    Normal mean expression
      (from the population analysis)
    Cancer mean expression
      (from the population analysis)
    % change in population
    Direction: always DOWN
    Depth correlation r values
    Which genes are in the 3-gene
    clinical panel

  FA MARKER PANEL:
    Gene names
    Normal mean expression
    Cancer mean expression
    Direction: always UP
    Depth correlation r values

  EPIGENETIC LOCK PROFILE:
    EZH2 direction in this cancer
      (elevated = gain-of-function lock;
       suppressed = loss-of-function)
    KDM1A direction
    HDAC direction
    Which lock type dominates

  DEPTH SCORE FORMULA:
    Exact formula used in Script 2.
    Which genes are the components.
    How normalisation is applied.
    The 0-1 scaling method.

  REFERENCE POPULATION STATISTICS:
    Mean depth score in the population
    Std of depth score in population
    Q25, Q50, Q75 of depth score
    These are needed to compute
    the patient's percentile rank.

  ATTRACTOR TYPE CLASSIFICATION:
    Differentiation or survival.
    The distinguishing criteria.

  DRUG MAP:
    Depth strata (thresholds).
    Drug class per stratum.
    Contraindications per stratum.
    Literature confirmation for each.

STORE AS REFERENCE FILE:
  reference_geometries/
    [cancer_type]_reference.json

  This file is used by the
  patient scoring script.
  It is derived from [N]b.
  It is not changed between patients.
  If the population analysis is
  updated, the reference file
  is regenerated.
  The version used for each patient
  analysis is recorded in the
  case intake record.
```

---

## PART VII — THE PATIENT SCORING COMPUTATION
### The Analysis Step by Step

---

### Step 7.1 — Normalise the Patient Data

```
PURPOSE:
  Place the patient's expression values
  in the same coordinate system as
  the reference population.

METHOD:
  For each gene g in the patient data:

    patient_z[g] = (patient_expr[g] -
                    reference_normal_mean[g]) /
                   reference_normal_std[g]

  This expresses the patient's expression
  as standard deviations from the
  reference normal.

  Positive z-score: gene elevated
  relative to normal.
  Negative z-score: gene suppressed
  relative to normal.

  If reference_normal_std[g] is zero
  or very small (< 0.01):
    Use reference_cancer_std[g] instead.
    Flag this gene in the output as
    "low variance in normal reference —
     z-score less reliable."

RECORD:
  The normalisation method used.
  Any genes flagged for low variance.
  The z-score for every gene
  in the analysis panel.
```

### Step 7.2 — Compute the Depth Score

```
PURPOSE:
  A single scalar (0 to 1) measuring
  how deeply this patient's cancer
  cells are trapped in the false attractor.

METHOD:

  COMPONENT 1 — SWITCH GENE SUPPRESSION:
    For each switch gene s in the panel:
      suppression[s] = -patient_z[s]
      (negative z = suppressed = positive
       contribution to depth)
      Clip to 0 if patient_z[s] > 0
      (a switch gene that is not suppressed
       does not contribute positively to depth)

    switch_component = mean(suppression[s]
                            for s in switch_panel)

    Normalise to 0-1 using reference
    population parameters:
      switch_component_norm =
        (switch_component - ref_switch_min) /
        (ref_switch_max - ref_switch_min)
      Clip to [0, 1].

  COMPONENT 2 — FA MARKER ELEVATION:
    For each FA marker f in the panel:
      elevation[f] = patient_z[f]
      Clip to 0 if patient_z[f] < 0
      (an FA marker that is not elevated
       does not contribute to depth)

    fa_component = mean(elevation[f]
                        for f in fa_panel)

    Normalise to 0-1 using reference
    population parameters:
      fa_component_norm =
        (fa_component - ref_fa_min) /
        (ref_fa_max - ref_fa_min)
      Clip to [0, 1].

  DEPTH SCORE:
    depth_score = mean(switch_component_norm,
                       fa_component_norm)

  PERCENTILE RANK:
    Using the reference population
    depth score distribution:
      depth_percentile = percentile rank of
      depth_score in the reference
      cancer population distribution.

    Report both the scalar (0-1) and
    the percentile (0-100th).

  INTERPRETATION:
    0.00 - 0.30: Shallow attractor.
                 Cells not far from normal.
                 Minimal epigenetic lock.
                 Differentiation therapy
                 may work with minimal
                 epigenetic augmentation.

    0.30 - 0.55: Moderate depth.
                 Partial switch gene suppression.
                 Intermediate lock.
                 Combined differentiation +
                 epigenetic agent likely needed.

    0.55 - 0.75: Deep attractor.
                 Strong switch gene suppression.
                 Significant epigenetic lock.
                 Epigenetic agent is a priority.

    0.75 - 1.00: Very deep attractor.
                 Near-complete switch gene loss.
                 Strong epigenetic lock.
                 Epigenetic intervention may need
                 to precede differentiation therapy.

RECORD:
  switch_component_raw
  switch_component_norm
  fa_component_raw
  fa_component_norm
  depth_score (final scalar)
  depth_percentile
  interpretation tier
  Any genes missing from the panel
  and how this affects the score.
```

### Step 7.3 — Profile Each Switch Gene Individually

```
PURPOSE:
  Move from the scalar to the gene-level.
  Understand which switches are most
  suppressed, partially retained,
  or near-normal.

METHOD:
  For each switch gene s:
    patient_value[s] = patient_expr[s]
    reference_normal[s] = ref_normal_mean[s]
    reference_cancer_median[s] = ref_cancer_q50[s]
    patient_z[s] = normalised z-score

    CLASSIFICATION:
      patient_z[s] < -2.0:
        DEEPLY SUPPRESSED
        This gene is more suppressed than
        the majority of the reference
        cancer population.
        Epigenetic lock likely strong here.

      -2.0 ≤ patient_z[s] < -1.0:
        SIGNIFICANTLY SUPPRESSED
        Below normal, within typical
        cancer range.
        Lock present.

      -1.0 ≤ patient_z[s] < -0.3:
        MODERATELY SUPPRESSED
        Partial loss.
        May be partially recoverable.

      -0.3 ≤ patient_z[s] < 0.3:
        NEAR NORMAL
        This switch gene is largely
        retained in this patient.
        Not a dominant contributor
        to the false attractor in
        this specific case.

      patient_z[s] ≥ 0.3:
        UNEXPECTED — NOT SUPPRESSED
        Flag explicitly.
        State: "This gene was predicted
        to be suppressed in this cancer
        type but is near normal or
        elevated in this patient.
        This may indicate a different
        attractor subtype or a data
        quality issue with this gene."

RECORD:
  Full table: gene, patient value,
  reference normal, reference cancer
  median, z-score, classification.
  Flag any unexpected direction findings.
```

### Step 7.4 — Profile the Epigenetic Locks

```
PURPOSE:
  Identify what is maintaining
  the switch gene suppressions.
  This points toward the epigenetic
  agent class most likely to dissolve
  the lock.

METHOD:
  For each epigenetic lock gene:
    EZH2, KDM1A (LSD1), HDAC1, HDAC2,
    DNMT3A, TET2, ASXL1

    Compute patient_z[gene].

  Compare to reference direction
  for this cancer type
  (from the reference geometry file):
    If reference says EZH2 UP
    and patient EZH2 is strongly positive:
      CONSISTENT WITH EXPECTED LOCK.
      EZH2 inhibition is predicted.

    If reference says EZH2 UP
    and patient EZH2 is near zero or negative:
      DIVERGES FROM EXPECTED LOCK.
      EZH2 may not be the dominant
      lock in this patient.
      Flag explicitly.
      State: "This patient's epigenetic
      lock profile diverges from the
      typical profile for this cancer type.
      The dominant lock mechanism may
      differ. This warrants discussion
      with the clinical team."

CLASSIFICATION:
  PRC2-DOMINANT LOCK:
    EZH2 patient_z > +1.0.
    Predicted: EZH2/PRC2 inhibitor
    (tazemetostat class).

  COREST-DOMINANT LOCK:
    KDM1A patient_z > +1.0.
    Predicted: LSD1 inhibitor
    (iadademstat class).

  HDAC-DOMINANT LOCK:
    HDAC1 or HDAC2 patient_z > +1.0.
    Predicted: HDAC inhibitor class.

  MIXED LOCK:
    Multiple epigenetic genes elevated.
    More than one lock type present.
    Combination epigenetic therapy
    may be required.

  UNUSUAL PROFILE:
    Pattern does not match any
    standard classification.
    State: "Epigenetic lock profile
    is atypical for this cancer type.
    This finding should be flagged
    for clinical discussion."

RECORD:
  Full epigenetic gene profile table.
  Lock type classification.
  Any divergence from expected noted.
```

### Step 7.5 — Classify the Attractor Type

```
PURPOSE:
  Determine whether this patient's
  cancer is a differentiation attractor
  or a survival attractor.
  These require different therapeutic logic.

DIFFERENTIATION ATTRACTOR:
  Definition: Cells blocked before
  reaching terminal differentiation.
  They are immature, proliferating,
  expressing progenitor markers.
  The switch genes that are suppressed
  are terminal differentiation genes.

  Criteria:
    Terminal differentiation TFs suppressed.
    Proliferation markers elevated.
    Progenitor/stem markers retained.

  Therapeutic logic:
    Dissolve the epigenetic block.
    Restore the differentiation programme.
    Let the cells complete their journey.
    Differentiation agents +
    epigenetic lock dissolution.

SURVIVAL ATTRACTOR:
  Definition: Cells that have reached
  mature identity but cannot execute
  programmed death.
  They look mature but are immortal.
  The switched-off genes are the
  apoptotic exit genes.

  Criteria:
    Anti-apoptotic genes elevated
    (BCL2, MCL1, BCL-XL).
    Mature identity markers present.
    Apoptotic pathway components suppressed.

  Therapeutic logic:
    Restore the apoptotic exit.
    Remove the survival signal.
    BH3 mimetics (venetoclax class).
    Anti-survival agents.

CLASSIFICATION METHOD:
  Use the reference geometry file
  for this cancer type.
  The attractor type is determined
  at the cancer-type level.

  Confirm with patient-specific data:
    Does the patient's profile
    match the expected type?
    If yes: proceed with the
    classification from the reference.
    If no: flag the divergence.

RECORD:
  Attractor type classification.
  Confirmation that patient data
  matches the expected type.
  Any divergence noted.
```

### Step 7.6 — Locate on the Drug Map

```
PURPOSE:
  Translate the geometric position
  into specific drug predictions.

METHOD:
  Using the drug map from the
  reference geometry file:

  LOCATE DEPTH STRATUM:
    Compare patient depth_score to
    the stratum thresholds in the
    reference drug map.
    Assign to the appropriate stratum.

  READ THE DRUG MAP:
    What drug class does this stratum
    predict as primary?
    What does it predict as secondary?
    What is contraindicated at this depth?

  CHECK ATTRACTOR TYPE:
    Does the attractor type
    (differentiation vs survival)
    modify the drug prediction?
    Apply any attractor-type
    modifications from the reference.

  CHECK EPIGENETIC LOCK:
    Does the specific lock type
    (PRC2 vs CoREST vs HDAC)
    modify the drug prediction?
    Apply specific agent class
    from the lock profile.

  COMPILE PREDICTION:
    Primary drug class prediction.
    Rationale from geometry.
    Secondary prediction if applicable.
    Contraindications from geometry.
    Literature confirmation status
    for each prediction
    (from Document [N]c).

CRITICAL STATEMENT ON DRUG PREDICTIONS:
  The drug map predictions are
  geometric predictions derived from
  population-level analysis.
  They are NOT treatment recommendations.
  They are predictions about which
  molecular targets the geometry suggests
  are active at this patient's depth position.
  Whether these agents are appropriate
  for this specific patient is a
  clinical decision requiring:
    Full medical history.
    Performance status.
    Organ function.
    Prior treatment history.
    Comorbidities.
    Drug availability and access.
  None of these are available to
  the geometric analysis.
  All clinical decisions remain with
  the patient's medical team.

  The drug map prediction generates
  specific questions for the clinical team.
  It does not replace the clinical team.

RECORD:
  Depth stratum assigned.
  Drug map position.
  Primary prediction with rationale.
  Secondary prediction with rationale.
  Contraindications.
  Honest uncertainty about each prediction.
```

### Step 7.7 — Serial Trajectory Analysis (If Applicable)

```
PURPOSE:
  If the patient has data from multiple
  time points, compute the geometric
  trajectory and detect directional signals.

METHOD:
  Run Steps 7.1 through 7.6 independently
  for each time point.
  Record depth_score at each time point.
  Record switch gene profiles at each
  time point.

  COMPUTE VELOCITY:
    Velocity = depth_score[t2] - depth_score[t1]
    Negative velocity = depth decreasing =
    attractor dissolving = treatment working.
    Positive velocity = depth increasing =
    attractor reconstituting = concern.

  COMPUTE ACCELERATION:
    If three or more time points:
    Acceleration = velocity[t2-t3] - velocity[t1-t2]
    Negative acceleration of dissolution =
    rate of response is slowing.
    This is an early warning signal even
    if the trajectory is still negative.

  TRAJECTORY CLASSIFICATION:

  DISSOLVING TRAJECTORY:
    Depth score falling across time points.
    Switch genes reactivating.
    FA markers falling.
    Interpretation: treatment is working
    at the geometric level.

  STABLE HIGH TRAJECTORY:
    Depth score not changing despite treatment.
    Switch genes remain suppressed.
    FA markers remain elevated.
    Interpretation: the attractor is
    not responding to current treatment.
    The epigenetic lock may not be targeted.
    Bring to clinical team: is the current
    treatment addressing the lock mechanism
    identified in the epigenetic profile?

  RECONSTITUTING TRAJECTORY:
    Depth score falling then rising.
    Switch gene reactivation stalling
    or reversing.
    Interpretation: early geometric
    relapse signal.
    This may precede clinical detection
    by weeks to months.
    This is the highest-priority finding
    in the entire protocol.
    Communicate to clinical team
    immediately with the specific data.

  SHIFTING TRAJECTORY:
    The dominant gene correlates changing.
    A different set of genes driving
    depth at later time points.
    Interpretation: clonal evolution.
    A new cancer subpopulation with
    a different attractor geometry
    is emerging.
    The original reference geometry
    may no longer fully describe the
    dominant population.
    Flag for clinical discussion.
    Consider whether a new reference
    analysis is needed.

RECORD:
  Depth score at each time point.
  Velocity between each pair of points.
  Acceleration if three or more points.
  Trajectory classification.
  Time to signal if reconstituting
  trajectory detected.
```

---

## PART VIII — THE GEOMETRIC REPORT
### Format and Required Sections

---

```
The geometric report is the deliverable.
It is produced for every patient
who completes the analysis.
It is written for two audiences
simultaneously:
  The patient — in plain language.
  The clinical team — in precise language.

Both audiences receive the same document.
No separate versions.
One document. Complete. Honest.
```

### Report Structure

````markdown name=Individual_Protocol/GEOMETRIC_REPORT_TEMPLATE.md
# GEOMETRIC REPORT
## OrganismCore Attractor Framework
## Individual Patient Analysis

**Patient ID:** [Anonymous — assigned at intake]
**Cancer Type:** [lineage name]
**Biopsy Date:** [date of tissue collection]
**Analysis Date:** [date this report was produced]
**Reference Geometry Version:** [cancer type] — [document number and date]
**Analyst:** Eric Robert Lawson, OrganismCore
**ORCID:** https://orcid.org/0009-0002-0414-6544
**Framework Repository:** https://github.com/Eric-Robert-Lawson/attractor-oncology

---

## IMPORTANT — READ FIRST

```
This report is a geometric measurement
of your cancer attractor state.

It is not a medical diagnosis.
It is not a treatment recommendation.
It is not a clinical prognosis.

It is a measurement — derived from your
gene expression data using a mathematical
framework validated across 22+ cancer types —
of where your cancer cells sit in the
Waddington epigenetic landscape.

All clinical decisions remain entirely
with your licensed medical team.

The purpose of this report is to give you
and your clinical team additional geometric
information that is not available from
standard clinical instruments.

Read the findings. Bring them to your
oncologist. Ask the questions listed in
Section 7. The clinical team decides
what to do with the information.
```

---

## SECTION 1 — DEPTH SCORE

```
YOUR DEPTH SCORE: [value 0.00 to 1.00]

PERCENTILE: [value]th percentile in the
reference [cancer type] population.

This means your cancer cells are more
deeply trapped in their false attractor
than [percentile]% of patients in the
reference dataset.

INTERPRETATION: [Shallow / Moderate /
                 Deep / Very Deep]

PLAIN LANGUAGE:
[One paragraph in non-technical language
explaining what the depth score means
for this specific patient.
What are their cells trying to do?
How far have they traveled from normal?
What does this depth suggest about
the nature of their specific disease?]
```

---

## SECTION 2 — SWITCH GENE PROFILE

```
Switch genes are the genes that would
need to be active for your cancer cells
to complete their developmental journey
to their normal mature identity.
In your cancer cells, these genes are
suppressed — held off by the molecular
locks maintaining your cancer state.

[TABLE FORMAT:]

Gene    | Normal   | Your Value | Change  | Status
--------|----------|------------|---------|------------------
[gene1] | [value]  | [value]    | [%]     | [classification]
[gene2] | [value]  | [value]    | [%]     | [classification]
[gene3] | [value]  | [value]    | [%]     | [classification]
...

CLASSIFICATION KEY:
  DEEPLY SUPPRESSED: >2 standard deviations below normal
  SIGNIFICANTLY SUPPRESSED: 1-2 SD below normal
  MODERATELY SUPPRESSED: 0.3-1 SD below normal
  NEAR NORMAL: within 0.3 SD of normal
  UNEXPECTED: elevated when expected to be suppressed — flagged

SUMMARY:
[Which switch genes are most suppressed.
Which appear partially retained.
What this pattern suggests about the
specific geometry of this patient's
cancer state.]
```

---

## SECTION 3 — EPIGENETIC LOCK PROFILE

```
The epigenetic lock is what is maintaining
the suppression of your switch genes.
It is the molecular mechanism holding
your cancer cells in their trapped state.

DOMINANT LOCK TYPE: [PRC2 / CoREST / HDAC / Mixed / Atypical]

[TABLE FORMAT:]

Lock Gene | Your Value | Direction | vs Expected | Interpretation
----------|------------|-----------|-------------|---------------
EZH2      | [value]    | [UP/DOWN] | [MATCH/FLAG]| [interpretation]
KDM1A     | [value]    | [UP/DOWN] | [MATCH/FLAG]| [interpretation]
HDAC1     | [value]    | [UP/DOWN] | [MATCH/FLAG]| [interpretation]
HDAC2     | [value]    | [UP/DOWN] | [MATCH/FLAG]| [interpretation]
DNMT3A    | [value]    | [UP/DOWN] | [MATCH/FLAG]| [interpretation]

PLAIN LANGUAGE:
[What the dominant lock means in plain terms.
Which epigenetic agent class the geometry
suggests as relevant.
Any divergence from expected pattern noted
and explained honestly.]
```

---

## SECTION 4 — ATTRACTOR CLASSIFICATION

```
ATTRACTOR TYPE: [DIFFERENTIATION / SURVIVAL]

[DIFFERENTIATION ATTRACTOR — if applicable:]
Your cancer cells are blocked partway through
their developmental journey. They are stuck
before reaching their mature identity.
They are trying to become [normal cell type]
but cannot complete the journey.
The therapeutic logic for this geometry is:
dissolve the epigenetic block and restore
the developmental programme.

[SURVIVAL ATTRACTOR — if applicable:]
Your cancer cells have reached their mature
identity but cannot execute their natural
programmed death. They appear mature but
are artificially immortal.
The therapeutic logic for this geometry is:
restore the apoptotic exit pathway.
Remove the survival signal that is keeping
them alive past their natural endpoint.

PATIENT-SPECIFIC CONFIRMATION:
[Does this patient's data confirm the
expected attractor type for their cancer?
Any divergence noted and explained.]
```

---

## SECTION 5 — DRUG MAP POSITION

```
DEPTH STRATUM: [stratum name from drug map]

Based on your depth score and attractor
geometry, your cancer sits in the
[stratum] stratum of the [cancer type]
attractor landscape.

The geometry-derived drug map for this
stratum predicts:

PRIMARY PREDICTION:
  Target: [gene/pathway]
  Drug class: [class]
  Mechanism: [how this target dissolves
              the attractor at this depth]
  Literature status: [CONFIRMED BY CLINICAL
                      TRIAL / FDA APPROVED /
                      PHASE X TRIAL /
                      NOVEL PREDICTION]

SECONDARY PREDICTION:
  Target: [gene/pathway]
  Drug class: [class]
  Mechanism: [how this augments the primary]
  Literature status: [same classification]

EPIGENETIC AUGMENTATION:
  Based on your lock profile:
  [Lock-specific agent class]
  This is predicted to dissolve the
  specific epigenetic lock identified
  in Section 3.

GEOMETRIC CONTRAINDICATIONS:
  Based on your depth geometry, the
  following drug classes are predicted
  to be less effective at your specific
  attractor position:
  [List with geometric rationale]

CRITICAL STATEMENT:
  These are geometric predictions.
  They are not treatment recommendations.
  Whether any of these agents are
  appropriate for you specifically
  requires your full medical history,
  current treatment status, organ function,
  and clinical judgment.
  These predictions are provided to
  generate informed questions for your
  clinical team — not to replace them.
```

---

## SECTION 6 — TRAJECTORY (SERIAL DATA ONLY)

```
[Include only if multiple time points exist]

TRAJECTORY CLASSIFICATION:
[DISSOLVING / STABLE HIGH /
 RECONSTITUTING / SHIFTING]

TIME POINT DATA:

Date         | Depth Score | Velocity | Classification
-------------|-------------|----------|-----------------
[date 1]     | [score]     | —        | Baseline
[date 2]     | [score]     | [+/-]    | [interpretation]
[date 3]     | [score]     | [+/-]    | [interpretation]

TRAJECTORY INTERPRETATION:
[What the direction of change means
geometrically.
Is the attractor dissolving?
Is it reconstituting?
Is the rate of change accelerating
or decelerating?
Plain language explanation followed
by precise geometric statement.]

[IF RECONSTITUTING TRAJECTORY DETECTED:]
⚠️ PRIORITY SIGNAL:
The geometric trajectory shows the
attractor is reconstituting.
This is an early geometric signal
that may precede clinical detection
of progression.
Bring this finding to your oncologist
immediately with this report.
The specific data is in the table above.
This is the scenario where early
action has the most impact.
```

---

## SECTION 7 — QUESTIONS FOR YOUR CLINICAL TEAM

```
These questions are derived directly from
your geometric findings.
They are specific to your data.
They are designed to enable a more
informed conversation with your oncologist.

They are questions — not instructions.
What to do with the answers is your
clinical team's decision.

QUESTION 1: [Derived from switch gene profile]
[Specific question using the actual gene names
and values from this patient's analysis.]

QUESTION 2: [Derived from epigenetic lock profile]
[Specific question about the dominant lock
and whether it is being targeted.]

QUESTION 3: [Derived from depth score and drug map]
[Specific question about whether the depth stratum
changes the treatment selection rationale.]

QUESTION 4: [Derived from attractor classification]
[Specific question about differentiation vs
survival logic in current treatment.]

QUESTION 5: [Derived from any unexpected findings]
[If any unexpected findings were flagged —
a specific question about their significance.]

[IF SERIAL DATA — TRAJECTORY QUESTION:]
QUESTION 6: [Derived from trajectory]
[Specific question about what the geometric
trajectory means for the treatment plan
and whether it warrants any adjustment.]
```

---

## SECTION 8 — HONEST UNCERTAINTY

```
What the framework is confident about
in this analysis:

[List specific findings that are strongly
supported by the data and the reference
geometry. Be specific. Use the actual
values and gene names.]

What the framework is less certain about:

[List specific findings where the data
was ambiguous, where coverage was partial,
where the patient's profile diverged from
the expected pattern, or where the
prediction is novel rather than confirmed
by literature. Be equally specific.]

What this analysis cannot determine:

[Explicit list of things outside the
scope of the geometric analysis.
Clinical factors. Drug interactions.
Performance status. Comorbidities.
Anything that requires clinical judgment
beyond the expression data.]

Data quality notes:

[Any genes that were missing from the panel.
Any normalisation issues.
Any flags raised during data verification.
How these affect the confidence of the report.]
```

---

## SECTION 9 — FRAMEWORK STATUS AND CONSENT REMINDER

```
The OrganismCore attractor framework
has been validated retrospectively across
22+ cancer types using public gene expression
datasets. In every validated cancer type,
the framework correctly derived drug targets
confirmed by published pharmacology and
active clinical trials without prior
knowledge of the pharmacology.

The framework has not been validated in
a prospective clinical trial.
It has not been approved as a clinical
diagnostic instrument by any regulatory body.

This report is produced under the terms
of the informed consent document signed
on [consent date].

This report does not constitute medical
advice, diagnosis, or treatment recommendation.
All clinical decisions remain with your
licensed medical team.

Repository: https://github.com/Eric-Robert-Lawson/
            attractor-oncology

All analysis code, reference geometries,
and methodology are publicly available
and reproducible.

Contact: OrganismCore@proton.me
```

---

## DOCUMENT METADATA

```
Report version: 1.0
Framework version: [reference geometry version]
Analysis date: [date]
Analyst: Eric Robert Lawson
         OrganismCore
ORCID: https://orcid.org/0009-0002-0414-6544
```
````

---

## PART IX — CASE DOCUMENTATION
### What Is Preserved for Every Case

---

```
Every patient case produces a complete
documentation record.

The documentation serves four purposes:
  1. The patient can understand exactly
     how their report was produced.
  2. The analyst can reproduce the analysis
     at any future time point.
  3. The evidence base accumulates
     systematically.
  4. The framework improves as patterns
     across cases are observed.

CASE DIRECTORY STRUCTURE:

  CASE_[anonymous_ID]/
  │
  ├── intake_record.md
  │     Patient ID (anonymous)
  │     Cancer type
  │     Date of consent
  │     Date of data receipt
  │     Data format and quality tier
  │     Reference geometry version used
  │     Time points if serial
  │
  ├── consent/
  │     signed_consent.pdf
  │     (signed informed consent document)
  │
  ├── raw_data/
  │     [patient expression files —
  │      as received, completely unmodified]
  │     data_receipt_verification.md
  │     (checksums of received files)
  │
  ├── analysis/
  │     normalised_patient_data.csv
  │     depth_score_computation.csv
  │     switch_gene_profile.csv
  │     epigenetic_lock_profile.csv
  │     drug_map_position.csv
  │     trajectory_data.csv (if serial)
  │     analysis_log.txt
  │     (complete log of all computations —
  │      every step, every value, every
  │      decision recorded)
  │
  ├── report/
  │     GEOMETRIC_REPORT_[anonymous_ID].md
  │     (the final delivered report)
  │     report_delivery_record.md
  │     (date delivered, method, recipient)
  │
  └── followup/ (populated over time)
        [any clinical outcomes the patient
         chooses to share]
        [any corrections or updates to the
         report based on new information]
        [trajectory updates from new biopsies]

RETENTION:
  Case records are retained indefinitely
  unless the patient requests deletion.
  Deletion is immediate and complete
  upon request.
  No copies are retained after deletion.

ANONYMISATION:
  Patient ID codes are random.
  No personal identifiers in any
  file names or file contents.
  The mapping between patient ID code
  and patient identity exists only in
  the signed consent document.
  The consent document is stored
  separately from all analysis files.
```

---

## PART X — OUTCOME TRACKING
### Building the Evidence Base

---

```
The individual patient analyses are
not only service — they are evidence.

With patient consent, outcome tracking
builds the prospective evidence base
that will eventually validate the framework
at the clinical instrument level.

WHAT IS TRACKED (with consent):

  GEOMETRIC MEASUREMENTS:
    All depth scores at all time points.
    All trajectory classifications.
    All drug map positions.

  CLINICAL OUTCOMES (if patient shares):
    Treatment received.
    Response by clinical criteria.
    Whether geometric trajectory correlated
    with clinical response.
    Time to progression if it occurred.
    Survival at 6, 12, 24 months.

  CONCORDANCE ANALYSIS:
    Did the depth score stratum predict
    the correct drug class for this patient?
    Did the trajectory signal precede
    clinical detection of progression?
    Was the epigenetic lock target relevant
    to the treatment that worked?

  DISCORDANCE ANALYSIS:
    Where did the framework prediction
    diverge from clinical outcome?
    Why?
    What does the divergence teach?
    How should the reference geometry
    be updated?

OUTCOME TRACKING CONSENT:
  Outcome tracking requires additional
  explicit consent.
  It is separate from the analysis consent.
  It is always optional.
  The geometric analysis is provided
  regardless of whether the patient
  consents to outcome tracking.
  The analysis is never conditional
  on outcome tracking consent.

THE ACCUMULATION GOAL:
  Every case where geometric trajectory
  was measured AND clinical outcome is known
  is one data point toward prospective
  validation.

  50 cases with trajectory data and
  outcomes is a small dataset but
  a real prospective signal.

  100 cases begins to be statistically
  meaningful.

  This is how the framework moves from
  retrospective validation to prospective
  validation.

  One patient at a time.
  One case record at a time.
  One outcome documented at a time.

  The evidence base is built from
  the same patients who are helped.
  The help and the science are
  not separate activities.
  They are the same activity.
```

---

## PART XI — WHEN TO STOP
### Recognising the Limits

---

```
There are situations where the
individual patient protocol should
not proceed or should stop.

STOP BEFORE ANALYSIS:

  The reference geometry does not
  exist for the patient's cancer type.
  Action: establish reference geometry
  first using Workflow_Protocol.md.

  The data quality is insufficient
  (< 50% required gene coverage).
  Action: explain to patient what
  additional data is needed.
  Do not produce a report from
  insufficient data.

  The informed consent is incomplete.
  Action: do not accept data.
  Obtain complete consent first.

STOP DURING ANALYSIS:

  The patient's data shows no
  recognisable signal for their
  cancer type.
  Possible causes: wrong cancer type
  classification, sample mislabelling,
  significant data quality issue.
  Action: stop. Communicate with patient.
  Request clarification of diagnosis
  and sample identity before proceeding.

  The patient's expression profile is
  completely unlike any case in the
  reference population.
  Action: proceed with extreme caution.
  Produce a report that explicitly states
  the degree of divergence from reference.
  Do not force a classification that
  the data does not support.

STOP AFTER ANALYSIS:

  Never modify the report after delivery
  to make it more or less alarming than
  the data supports.
  If a correction is needed, issue a
  clearly labelled amended report
  with the reason for the amendment
  documented completely.

THE HARD LIMIT:

  If at any point in the process the
  analyst believes the patient may be
  using or intending to use the geometric
  report as a substitute for all clinical
  care rather than a supplement to it —

  Stop.
  Communicate directly with the patient.
  State clearly:
    "This report is designed to work
     alongside your clinical team.
     It is not designed to replace it.
     I am concerned that you may be
     planning to use this report without
     clinical oversight. I cannot in
     good conscience provide a report
     under those circumstances.
     Please engage with your medical team
     first. Then I can provide the report
     as an additional geometric perspective
     on what they are treating."

  This is not paternalism.
  It is the founding principle applied
  to its hardest case.

  "I do not want to take your money
   and promise you bullshit."

  Providing a report to a patient who
  intends to replace clinical care with it
  is promising them something the report
  cannot deliver.

  The report is a geometric instrument.
  It is powerful within its domain.
  It is not a clinical substitute.
  It was never designed to be.
```

---

## PART XII — THE FOUNDING PRINCIPLE
### Why This Protocol Exists in This Form

---

```
This protocol is as detailed as it is
because precision protects people.

Not precision as bureaucracy.
Precision as care.

Every section in this protocol exists
because of a specific harm that
imprecision would cause:

  The informed consent section exists
  because a person who is desperate
  and frightened deserves to know
  exactly what they are receiving
  before they receive it.
  Not after.

  The data quality verification exists
  because a report built on insufficient
  data is worse than no report.
  A wrong geometric picture is more
  dangerous than no picture.

  The honest uncertainty section exists
  because a report that presents all
  findings with equal confidence is lying.
  Some findings are strong.
  Some are uncertain.
  The patient deserves to know which
  is which.

  The stop conditions exist because
  there are situations where proceeding
  causes harm regardless of intent.
  Knowing when to stop is as important
  as knowing how to proceed.

  The outcome tracking section exists
  because helping one patient and
  building the evidence that helps
  the next patient are the same act.
  The evidence is not separate from
  the service.
  The service generates the evidence.

  The founding principle exists
  throughout every section because
  the founding principle is not
  a slogan — it is the structure
  of every decision this protocol makes.

  "I do not want to take people's money
   and promise them bullshit."

  That sentence is operationalised
  in every section of this protocol.

  The informed consent is that sentence
  made into a document.
  The data verification is that sentence
  made into a quality check.
  The honest uncertainty section is
  that sentence made into a required
  report component.
  The stop conditions are that sentence
  made into explicit decision rules.

  The protocol is the founding principle
  made precise enough to follow.
```

---

## PART XIII — THE RELATIONSHIP TO THE
##             POPULATION ANALYSES
### How Individual Cases Feed Back to the Framework

---

```
The population analyses establish
the reference geometry.

The individual patient analyses test it.

Every patient case is a prospective
data point that either confirms or
challenges the reference geometry.

WHEN A CASE CONFIRMS THE REFERENCE:
  The depth score places the patient
  in a stratum consistent with their
  clinical picture.
  The trajectory correlates with
  clinical response as predicted.
  Record the confirmation.
  The reference geometry is validated
  in a real prospective case.

WHEN A CASE CHALLENGES THE REFERENCE:
  The patient's profile diverges from
  the reference in a specific way.
  The trajectory does not match the
  prediction.
  Do not ignore this.
  Do not suppress it.
  Document it completely.
  Ask: what does this divergence teach?

  Possible causes of divergence:
    The patient's cancer is a subtype
    not well-represented in the reference
    population dataset.
    The reference geometry needs to be
    refined for this subtype.
    The framework prediction was wrong
    for this specific case.
    The clinical data that was compared
    is incomplete or ambiguous.

  In all cases: document the divergence
  exactly as it is. The wrong prediction
  protocol from Workflow_Protocol.md
  applies here too.
  Wrong predictions are information.
  They teach.
  They improve the reference geometry.
  They make the next patient's
  analysis more accurate.

THE FEEDBACK LOOP:

  Population analysis → reference geometry
  ↓
  Individual patient analysis
  ↓
  Outcome tracking (with consent)
  ↓
  Confirmation or divergence documented
  ↓
  Reference geometry updated if indicated
  ↓
  Next population analysis more accurate
  ↓
  Next individual patient analysis
  more accurate

  This is how the framework improves.
  Not by decree.
  Not by assumption.
  By following the data of real patients
  back into the reference geometry
  and updating it with what was learned.

  Each patient helps the next.
  Not through their suffering.
  Through the precision with which
  their data is handled and the
  honesty with which the results
  are recorded.
```

---

## DOCUMENT METADATA

```
document_id:    INDIVIDUAL_PATIENT_GEOMETRIC_
                ANALYSIS_PROTOCOL
folder:         Individual_Protocol/
series:         OrganismCore — Protocol Documents
author:         Eric Robert Lawson
date:           2026-03-04
version:        1.0
status:         COMPLETE — active protocol
                Updated as methodology evolves.
                Version number incremented with
                each substantive change.
                Change log maintained below.

repository:     https://github.com/Eric-Robert-Lawson/
                attractor-oncology

orcid:          https://orcid.org/0009-0002-0414-6544

contact:        OrganismCore@proton.me

related_documents:
  HOW_THIS_HELPS_YOU_TODAY.md
  Patient_Geometric_Sovereignty_Reasoning_Artifact.md
  Personalized_Attractor_Medicine_Reasoning_Artifact.md
  Cancer_Research/Workflow_Protocol.md
  Cancer_Research/OrganismCore_Cancer_Framework.md

founding_principle:
  "I do not want to take people's money
   and promise them bullshit."
   — Eric Robert Lawson, March 4, 2026

change_log:
  v1.0 — 2026-03-04
    Initial version.
    Complete protocol from data receipt
    to report delivery.
    All sections specified.
    Founding principle operationalised
    throughout.
```

---

*"The medical system will tell you what worked*
*in a population of people like you.*
*Your geometry will tell you what is happening*
*to your cells.*
*Both are necessary.*
*In the scenario that matters most,*
*yours is the measurement that saves your life."*

— Patient Geometric Sovereignty Reasoning Artifact
   Eric Robert Lawson, March 4, 2026
