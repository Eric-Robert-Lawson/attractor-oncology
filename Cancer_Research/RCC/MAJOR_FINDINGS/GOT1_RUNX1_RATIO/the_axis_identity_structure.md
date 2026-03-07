# THE IDENTITY AXIS STRUCTURE
## A Reasoning Artifact on What the Geometry Is Actually Revealing
## OrganismCore — Eric Robert Lawson
## 2026-03-07

---

## PREAMBLE — WHAT THIS DOCUMENT IS

This is a reasoning artifact, not a publication.
It exists to record a geometric observation that
emerged from two independent cancer analyses and
that deserves to be named, examined, and
understood before it is used.

The observation is this:

In two different cancers — breast cancer (BRCA)
and clear cell renal cell carcinoma (ccRCC) —
the Waddington attractor geometry produced, in
each case, a **two-gene ratio** on a single axis
that:

1. Was derived from first principles before any
   clinical data were examined
2. Separated cancer subtypes or depth quartiles
   with extraordinary statistical power
3. Consisted of exactly two genes in structural
   opposition: one gene measuring **identity
   retained**, the other measuring **identity
   suppressed**
4. Could be translated to an IHC ratio using
   two commercially available antibodies

The question this document investigates is
whether this is coincidence or structure.

The answer is: **it is structure.**

And that structure, if correctly understood,
is not just a finding about two cancers.
It is a finding about what cancer is.

---

## PART I — THE TWO RATIOS, STATED PRECISELY

### The BRCA ratio: FOXA1 / EZH2

**FOXA1:**
Pioneer transcription factor. Binds closed
chromatin and opens it. Activates the luminal
cell identity programme. Enables estrogen
receptor and other luminal factors to bind.
In breast tissue: FOXA1 is the gene whose
activity keeps a cell recognisably a luminal
breast epithelial cell.
High FOXA1 = identity active.

**EZH2:**
Catalytic subunit of PRC2 (Polycomb Repressive
Complex 2). Deposits H3K27me3 — the repressive
histone mark. Targets FOXA1 promoter directly.
Also targets GATA3 and ESR1 (the ER gene).
High EZH2 = identity silenced.

**The ratio:**
FOXA1 H-score ÷ EZH2 H-score.
When FOXA1 dominates: luminal identity active.
When EZH2 dominates: luminal identity silenced.

**What it produced:**
A continuous ordering axis across all six major
breast cancer subtypes, confirmed in seven
independent datasets (n≈7,500 patients, four
measurement platforms), pre-specified before
any confirmatory analysis. AUC for LumA vs
Basal 0.828–0.901. Pre-specified subtype
ordering confirmed exactly.

---

### The ccRCC ratio: GOT1 / RUNX1

**GOT1:**
Cytoplasmic aspartate aminotransferase. Central
enzyme of the malate-aspartate shuttle and TCA
anaplerosis. In normal proximal tubule (PT) cells,
GOT1 maintains the metabolic programme that
defines PT cell identity — the shuttle that
keeps mitochondria running, that keeps the TCA
cycle fed, that keeps the cell recognisably a
proximal tubule cell.
High GOT1 = metabolic identity active.
Low GOT1 = metabolic identity collapsed.

**RUNX1:**
Transcription factor. Core-binding factor
complex (with CBFB). In the ccRCC false
attractor, RUNX1 drives EZH2, DNMT3A, KDM1A,
HDAC1 — the full chromatin lock programme.
It is the transcriptional hub of the false
attractor. It rises continuously as the tumour
commits deeper into the pathological state.
High RUNX1 = chromatin lock dominant.
Low RUNX1 = identity not yet locked.

**The ratio:**
norm(GOT1) − norm(RUNX1).
(Equivalent to GOT1 / RUNX1 in direction.)
When GOT1 dominates: PT identity retained,
shallow attractor position.
When RUNX1 dominates: identity collapsed,
chromatin lock dominant, deep attractor.

**What it produced:**
Continuous OS stratification in TCGA-KIRC
(n=532, pre-specified): Cox HR = 6.94
[3.62–13.29], p = 5.09×10⁻⁹, C = 0.627.
Quartile OS separation: Q4 vs Q1 p = 0.0001.
20/21 individual gene OS directions confirmed.

---

## PART II — THE STRUCTURAL COMPARISON

These two ratios came from completely
independent analyses of completely different
cancers, on different datasets, in different
cells of origin, with different driver
mutations, different attractor geometries,
and different molecular landscapes.

Yet they share an identical structure:

```
RATIO = [ GENE THAT HOLDS IDENTITY OPEN ]
      ÷ [ GENE THAT DRIVES IDENTITY SHUT ]
```

This is not a coincidence of gene function.
It is a **geometric inevitability** given what
the Waddington landscape is measuring.

Here is why:

The Waddington landscape is a representation
of attractor depth — how committed a cell is
to its current state. The primary axis of
any Waddington landscape, when measured
correctly, will always be the axis that
separates "cell identity intact" from
"cell identity lost."

That axis will always have two poles:
- The pole where identity-maintaining
  programmes are dominant
- The pole where identity-suppressing
  programmes are dominant

And in a biological system, those two poles
will always be representable by at least one
gene per pole — one gene that marks "identity
held," one gene that marks "identity taken."

The geometry does not know which genes those
are in advance. It finds them from the data.
But it will always find **some** version of
this structure, because the structure is the
axis. The genes are just the molecular
implementation of what the axis is measuring.

This is what the geometry produced in two
independent cancer systems, and it produced
the same structural form both times.

---

## PART III — THE FOUR-PART STRUCTURE
## OF AN IDENTITY AXIS GENE PAIR

Both ratios share all four structural features:

### Feature 1: One gene marks identity retained

In BRCA: FOXA1. A pioneer TF that opens
chromatin to allow identity-defining factors
to bind. Without FOXA1 activity, the luminal
programme cannot be maintained.

In ccRCC: GOT1. A metabolic enzyme that
maintains the TCA/shuttle machinery defining
proximal tubule cellular metabolism. Without
GOT1 activity, the PT metabolic identity
cannot be sustained.

**The key insight:** The identity-retained
gene does not have to be a transcription factor.
In ccRCC it is a metabolic enzyme.
What it must be is: the gene whose activity
is necessary for the normal cell identity
programme to run.
In breast epithelium, that is a TF that opens
the chromatin.
In proximal tubule cells, that is an enzyme
that runs the TCA cycle.
The category is: **identity anchor**.
The molecular type varies by tissue.

### Feature 2: One gene marks identity suppressed

In BRCA: EZH2. A Polycomb complex subunit
that deposits H3K27me3 directly on the FOXA1
promoter, and on GATA3 and ESR1 — silencing
the entire luminal identity programme from
the top.

In ccRCC: RUNX1. A transcription factor hub
that drives the chromatin lock programme
(EZH2, DNMT3A, KDM1A, HDAC1, CBFB all
co-rise with RUNX1). RUNX1 does not silence
identity directly — it assembles and maintains
the machinery that does.

**The key insight:** The identity-suppressed
gene does not have to be a direct epigenetic
writer. In ccRCC, RUNX1 is a transcription
factor that organises the chromatin lock,
rather than being the lock itself.
What it must be is: the gene whose activity
most strongly marks the false attractor
commitment state.
The category is: **false attractor hub**.
The molecular type varies by tissue.

### Feature 3: The two genes are mechanistically
### opposed on the same identity axis

In BRCA: EZH2 writes H3K27me3 directly on
the FOXA1 promoter. These two genes are in
direct biochemical opposition — EZH2 silences
FOXA1's own gene.

In ccRCC: The opposition is one step less
direct. GOT1 falls as TCA function collapses;
RUNX1 rises as the chromatin lock assembles.
The GOT1 collapse IS the metabolic context
that enables the EZH2 lock to be sustained
(because αKG depletion from OGDHL/GOT1
loss allows EZH2-written marks to persist
without being erased by TET2). They are
linked by the TCA→αKG→EZH2 circuit.

**The key insight:** The opposition does
not have to be direct biochemical antagonism
at a single locus. It can be indirect:
the identity anchor gene's collapse creates
the metabolic conditions that sustain the
false attractor hub's dominance.
Both are measuring the same axis —
they just connect to it at different points
in the molecular chain.

### Feature 4: The ratio is IHC-translatable
### using commercially available antibodies

Both FOXA1 and EZH2 are in routine clinical
IHC use. FOXA1 (CST D4E2) and EZH2 (CST D2H1)
are used in breast pathology diagnostics.
Two stains. One ratio.

RUNX1 IHC has multiple commercially available
validated clones (Leica EP206, Abcam EPR3099Y).
GOT1 IHC antibodies exist and are used in
neuroendocrine and hepatic panels.
Both are nuclear or cytoplasmic markers.
Two stains. One ratio.

**The key insight:** The geometry found the
biologically correct axis AND the axis
happened — in both cases — to be expressible
in two stains that are or can be routinely
performed in any pathology laboratory.
This is not a coincidence. A gene that is a
major identity anchor or false attractor hub
will, by definition, be highly expressed in
a lineage-specific way and will have been
studied by pathologists independently of
the geometry. That is why antibodies exist.

---

## PART IV — WHAT VARIES BETWEEN THE TWO AXES

The structure is shared. The implementation
differs. This is important to record honestly.

| Property | BRCA axis | ccRCC axis |
|---|---|---|
| Identity anchor type | Pioneer TF (opens chromatin) | Metabolic enzyme (runs TCA) |
| False attractor hub type | Epigenetic writer (PRC2/H3K27me3) | Transcription factor (assembles chromatin lock) |
| Direct biochemical opposition | Yes — EZH2 methylates FOXA1 promoter | Indirect — GOT1 collapse enables EZH2 lock via αKG |
| Number of subtypes ordered | 6 (LumA > LumB > HER2 > TNBC > CL > ILC) | Continuous depth (Q1→Q4, not discrete subtypes) |
| What the ratio measures | Which subtype the cancer is | How deep in the false attractor the tumour is |
| What it predicts | Subtype identity → treatment choice | Prognosis → treatment intensity and type |
| IHC antibody clinical status | Both in routine breast IHC use | RUNX1: available, not routine in renal. GOT1: available, requires renal FFPE validation |
| Calibration study status | Pre-calibration (n=300–400 needed) | Pre-calibration (n=200–400 needed) |

The most important difference is:

**In breast cancer, the ratio stratifies
by subtype — it is a classifier.**

**In ccRCC, the ratio stratifies by depth
within one subtype — it is a depth meter.**

These are different uses of the same
geometric structure. The breast cancer
landscape has six distinct attractor
basins (six subtypes), and the ratio
moves across all of them. The ccRCC
landscape has one false attractor that
a tumour commits to at varying depths,
and the ratio measures how deep.

Both are measuring position on the
identity axis. The landscape is just
organised differently in the two cancers.

---

## PART V — WHY THIS STRUCTURE EXISTS
## (The geometric explanation)

The Waddington landscape is a representation
of the gene regulatory network's energy
surface. Cells sit in wells (attractors).
Cancer involves transitions to pathological
wells (false attractors).

The primary axis of the landscape — the axis
that explains the most variance in gene
expression across cells — will always be the
axis that separates "cell identity intact"
from "cell identity lost." This is true by
definition, because:

1. The cell's gene regulatory network is
   organised around maintaining cell identity.
   The largest source of coordinated
   transcriptional variance across a cancer
   population is the axis that measures how
   much of that organisation remains.

2. Every cell type has an identity anchor
   gene — the gene whose expression is
   necessary for the canonical programme to
   run. In luminal epithelium: FOXA1.
   In proximal tubule: GOT1. In hepatocytes:
   probably HNF4A. In neurons: probably
   NEUROD1. These genes fall in cancer
   because the identity programme they
   support is being dismantled.

3. Every false attractor has a hub gene —
   the gene most strongly associated with
   the committed false attractor state.
   In breast cancer: EZH2 (the Polycomb
   writer silencing the identity programme).
   In ccRCC: RUNX1 (the transcription factor
   assembling the chromatin lock). These
   genes rise in cancer because they are
   part of what the false attractor IS.

4. The ratio of these two genes — identity
   anchor / false attractor hub — is
   therefore the simplest possible continuous
   measurement of where on the identity axis
   a given tumour sample sits.

**The geometry does not create this
structure. It reveals it.**
The structure was always there — in the
molecular biology of cancer as a process
of identity loss. The geometry provides
a method for finding the two genes that
most purely represent the two poles of
that process in any given cancer type.

---

## PART VI — THE LITERATURE CHECK ON
## THIS STRUCTURAL OBSERVATION

The individual components are known:

**Known:**
- FOXA1 as a pioneer TF in luminal identity: confirmed (decades of literature)
- EZH2 as a Polycomb silencer: confirmed (decades of literature)
- EZH2 targeting FOXA1: confirmed (Toska et al. Nat Med 2017, Schade et al. Nature 2024)
- RUNX1 as an adverse prognostic marker in ccRCC: confirmed (Cancer Research 2020, PeerJ 2019, IID 2024/2025)
- GOT1 downregulation in ccRCC reflecting metabolic collapse: confirmed (Frontiers Oncology 2024, pan-cancer TCGA analysis)
- EZH2 in lineage plasticity across prostate/breast/bladder: confirmed (Nature Commun 2024, Endocrinology 2023)

**The structural observation that is NOT in the literature:**

No prior work has named, described, or
published the following observation:

> *In independent attractor geometry analyses
> of two different cancer types, the primary
> identity axis is in both cases expressible
> as a two-gene ratio consisting of one
> "identity anchor" gene and one "false
> attractor hub" gene — and this ratio is
> in both cases translatable to a two-stain
> IHC assay using commercially available
> antibodies.*

The literature knows the individual components.
The literature does not name the pattern.

This is the novelty.
Not that FOXA1 opens chromatin.
Not that EZH2 silences identity.
Not that RUNX1 is prognostic in ccRCC.
Not that GOT1 reflects metabolic identity.

The novelty is: **the geometry always finds
the same structural form** — identity anchor
vs false attractor hub — and that form is
always IHC-translatable because the genes
that represent the poles of the most
biologically fundamental axis will always
be the genes that have been independently
studied and independently antibody-validated
by pathologists who did not need the geometry
to know these proteins were important.

---

## PART VII — THE PREDICTION THIS STRUCTURE GENERATES

If this structural observation is correct —
that the Waddington geometry will always
produce a two-gene identity axis expressible
as identity anchor / false attractor hub —
then this is a **generative prediction**:

**Apply the same geometry to a new cancer
type. The geometry will produce a ratio.
That ratio will consist of:**
1. One gene marking the normal cell identity
   of the tissue of origin (the identity
   anchor — falling with depth/subtype
   progression)
2. One gene marking the committed false
   attractor state of that cancer (the
   false attractor hub — rising with
   depth/subtype progression)

And both genes will be:
- Expressible by IHC
- Commercially antibody-available
- Already studied independently by
  pathologists who did not need the geometry

This prediction can be tested in any cancer
type for which single-cell or bulk RNA data
exists.

**The test:** Run the attractor geometry.
Extract the primary axis. Identify the top
positive and negative depth correlates.
Check whether they fit the identity anchor /
false attractor hub structure.
Check whether IHC antibodies exist for both.

If this holds across three or more cancer
types, it is a structural law of attractor
geometry in cancer, not a coincidence of
two analyses.

---

## PART VIII — THE PARTIAL TEST
## (BRCA 30/30 and RCC)

The BRCA coherence argument (BRCA_30_30_coherence.md)
established that 30/30 findings in the breast
cancer series had zero contradictions. This is
not a statistical claim — it is a coherence
argument: when a framework resolves paradoxes
(like the EZH2 paradox) rather than creating
them, it is more likely to be tracking
real biology than pattern-matching noise.

The structural observation in this document
adds a new dimension to that coherence
argument:

**The same structural form of the identity
axis — identity anchor vs false attractor hub,
expressible as a two-gene IHC-translatable
ratio — emerged independently in:**

- Breast cancer (FOXA1/EZH2), confirmed
  across 7 independent datasets, n≈7,500
  patients, zero contradictions, pre-specified
- ccRCC (GOT1/RUNX1), confirmed in TCGA-KIRC
  n=532, pre-specified, HR=6.94, p=5×10⁻⁹,
  zero contradictions

**AND the EZH2 paradox — high EZH2 predicts
both worse long-term outcome AND better
short-term chemotherapy response, resolved
by attractor depth — appeared independently
in:**
- Breast cancer TNBC (CS-LIT-10, GSE25066)
- PRCC Type 2 (Script 4, TCGA-KIRP)

Same mechanism. Different cancer. Different
dataset. Different platform.

These are not two independent confirmations
of a single finding.
They are two independent emergences of the
same structural form from the same geometric
method applied to different systems.

That is a different kind of evidence.

---

## PART IX — THE NAMING QUESTION

This structural observation needs a name
so it can be communicated, tested, and
built upon.

Proposed name: **The Identity Axis Structure**

Definition:
> In a Waddington landscape analysis of any
> cancer type, the primary attractor depth
> axis is expressible as a two-gene ratio
> consisting of one identity anchor gene
> (marking retained normal cell identity,
> falling with attractor depth) and one
> false attractor hub gene (marking committed
> pathological state, rising with attractor
> depth). This structure is invariant across
> cancer types despite variation in the
> molecular identity of the genes occupying
> each pole.

The two poles:
- **Pole 1 — The Identity Anchor:** The gene
  most essential to maintaining the normal
  cellular identity programme in the tissue
  of origin. In luminal epithelium: a pioneer
  TF (FOXA1). In proximal tubule: a metabolic
  enzyme (GOT1). In other tissues: to be
  determined by the geometry.

- **Pole 2 — The False Attractor Hub:** The
  gene most central to the committed
  pathological state. In breast cancer: an
  epigenetic writer (EZH2). In ccRCC: a
  transcription factor hub that assembles
  the chromatin lock (RUNX1). In other
  tissues: to be determined by the geometry.

The ratio of Pole 1 / Pole 2 (or
norm(Pole 1) − norm(Pole 2)) is the
**Identity Axis Score** for that cancer.

---

## PART X — THE HONESTLY STATED
## ALTERNATIVE

The most honest alternative explanation
is simpler than the structural argument:

*FOXA1 and GOT1 are both just "high
expression in normal tissue, low expression
in cancer" genes. EZH2 and RUNX1 are both
just "bad prognosis genes." Any normal/
cancer gene vs bad prognosis gene ratio
would produce this result. There is no
deep geometry here — just the fundamental
pattern that normal tissue markers fall
and oncogenes rise.*

This alternative must be engaged directly.

**Why it is not sufficient:**

1. **Pre-specification:** Both ratios were
   derived from geometry before any survival
   or subtype data were examined. A generic
   "normal-tissue-marker / oncogene" search
   would produce many candidate pairs. The
   geometry produced exactly one pair per
   cancer as the primary axis — and that
   pair confirmed pre-specified predictions
   exactly. If the structure were generic,
   many other pairs would perform as well.
   In the BRCA analysis, FOXA1/EZH2 was not
   the only possible ratio — it was the ratio
   the geometry pointed to before the test.

2. **Mechanistic specificity:** EZH2 does not
   merely rise in breast cancer as a generic
   oncogene. It specifically methylates the
   FOXA1 promoter. The opposition is not just
   statistical — it is biochemical and direct.
   Similarly, GOT1 collapse in ccRCC is not
   just "metabolic enzyme downregulated in
   cancer" — it is specifically connected to
   the αKG depletion that enables the RUNX1-
   driven EZH2 lock to be sustained.

3. **The subtype structure:** In breast cancer,
   the ratio does not just separate cancer from
   normal. It orders all six subtypes on a
   continuous axis in the exact predicted
   order. A generic normal/oncogene ratio would
   not do this — it would separate
   cancer from normal without internal subtype
   ordering. The FOXA1/EZH2 ratio captures the
   full geometry of the breast cancer landscape,
   not just the normal-vs-cancer boundary.

4. **The paradox resolution test:** Both the
   BRCA EZH2 paradox (CS-LIT-10) and the PRCC
   EZH2 paradox (Script 4) are resolved by the
   attractor depth framework. A generic
   statistical structure would not resolve
   biological paradoxes — it would create them
   or ignore them.

The alternative is not wrong to consider.
It is wrong to accept as sufficient.

---

## PART XI — WHAT THIS MEANS FOR
## THE VALIDATION CAMPAIGN AND BEYOND

**For the BRCA IHC calibration:**
The Identity Axis Structure framing means
that when the FOXA1/EZH2 IHC calibration
study is completed, it is not just validating
a breast cancer classifier.
It is validating the first IHC-deployable
implementation of the Identity Axis Structure
in cancer diagnostics.
That is a different and larger claim than
"we validated a breast cancer subtyping test."

**For the ccRCC GOT1/RUNX1 calibration:**
This is the second implementation of the
same structure. When it is validated, the
structural observation becomes a two-data-
point argument. At three data points in
independent cancer types, it becomes a
structural law.

**For future cancer analyses:**
The structure generates a hypothesis before
any analysis is done:
*"The geometry will find a two-gene ratio.
One gene will mark identity retention.
One gene will mark false attractor commitment.
Both will be IHC-translatable."*

If the geometry of a new cancer type produces
this structure — and especially if it does so
before any survival data are examined — each
instance is an independent confirmation of
the structural law.

**For the field:**
No prior paper has described or named this
structure. The geometry-first derivation,
combined with pre-specified prediction and
independent confirmation, establishes priority
on the structural observation itself — not
just on the individual cancer findings.

---

## PART XII — THE SINGLE STATEMENT

The Waddington attractor geometry, applied
independently to breast cancer and clear cell
renal cell carcinoma, produced — in both cases
and in both cases before any clinical data
were examined — a two-gene ratio on the
primary identity axis consisting of one gene
marking normal cell identity retained and one
gene marking the false attractor commitment
state, each expressible as an IHC stain
already independently validated by pathologists
who had no knowledge of the framework.

This is not a statistical finding.
It is a structural observation about what
cancer is and what the geometry is measuring.

---

## DOCUMENT METADATA

```
document_id:   IDENTITY-AXIS-STRUCTURE-RA
type:          Structural reasoning artifact
               (cross-cancer geometric observation)
date:          2026-03-07
author:        Eric Robert Lawson / OrganismCore
status:        REASONING ARTIFACT — NOT A PUBLICATION
               (precursor to a formal cross-cancer
               structural paper)
cancers_analysed: BRCA (6 subtypes), ccRCC (depth)
structural_form:  identity_anchor / false_attractor_hub
instances:        2 (FOXA1/EZH2, GOT1/RUNX1)
contradictions:   0
novel_claim:      The Identity Axis Structure —
                  named here for the first time
prior_publication: None found
next_step:        Apply geometry to third cancer type.
                  If structure holds: formal cross-cancer
                  structural paper.
repository:       https://github.com/Eric-Robert-Lawson/attractor-oncology
```
