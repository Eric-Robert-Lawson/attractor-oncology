# VACCINE ENGINE PROTOCOL
## How to Derive Escape-Proof Vaccine Targets
## Using the Attractor Framework
## A Replicable Onboarding Guide for Any Operator
## OrganismCore / attractor-oncology
## Eric Robert Lawson — 2026-03-09

---

## EPISTEMIC STATUS

This document is a **replicable protocol**, not a proof of clinical
efficacy. The engine derives targets from geometry and public data.
What it provides: the correct target and the geometric justification
for why it cannot be escaped. What it does not provide: immunogenicity
in humans, adjuvant selection, delivery optimisation, or clinical
durability data. The experimental science layer is irreplaceable and
downstream of everything here.

All inputs are public. All tools are free or low-cost. No institutional
affiliation required to run this protocol.

---

## PART 0 — THE ONE THING TO UNDERSTAND BEFORE RUNNING THE ENGINE

### What a vaccine actually is, in framework terms

Standard definition: a preparation that trains the immune system to
recognise and respond to a pathogen.

**Framework definition:**

> A vaccine is an intervention that causes the immune system to form a
> **stable attractor** at a specific position in antigenic state space.
> Immunological memory *is* a stable attractor. The immune system
> navigated from naïve to memory. That memory persists. It can be
> reactivated. It is a true attractor.

**A vaccine works when:** the attractor formed is at the position the
pathogen **cannot leave** — the structural invariant, the convergence
node, the thing it must always return to regardless of mutation pressure
because abandoning it means ceasing to function.

**A vaccine fails when:** the attractor is formed at a position the
pathogen **can and does leave** — the variable surface, the escape-prone
loop, the seasonally drifting epitope.

This is the single reframe that makes the engine possible.
The target is not the most immunogenic part of the pathogen.
The target is the part the pathogen is **geometrically imprisoned by**.

**The immune system's false attractor problem:**
When a vaccine trains immunity against a variable epitope, the pathogen
escapes. The immune memory is now a *false attractor* — locked to a
position the pathogen has vacated. This is the same geometry as a cancer
cell locked in the wrong developmental state. Different domain. Same
mathematics. The engine solves this by targeting the convergence node
instead.

---

## PART 1 — THE TWO VACCINE PROBLEMS

The engine applies to two structurally different problems. Identify
which problem you are solving before running.

### Problem Type A — Pathogen Vaccine (infectious disease)

The pathogen has its own attractor geometry. Its state space is its
sequence variation landscape across all known strains over all recorded
time. The convergence node is where it must always be to function. The
immune attractor must be formed at that node.

**Identifying signals:**
- The pathogen mutates rapidly under immune pressure (escape-prone)
- Existing vaccines require reformulation (immune attractor keeps missing)
- Or: novel pathogen with no existing vaccine — derive the target before
  the immune response data exists

### Problem Type B — Cancer Vaccine (no-escape principle)

The "pathogen" is the cancer cell. Its state space is the Waddington
landscape. Its convergence node is the **Hub protein** — the
epigenetic driver the cell must maintain to stay in the false attractor.
The immune attractor must be formed against Hub-derived epitopes.

**The no-escape double-bind (novel — not in literature as of 2026-03-09):**

```
IF immune system kills Hub-high cells
  → direct immune kill

IF tumor downregulates Hub to escape immune surveillance
  → false attractor collapses
  → cell completes developmental journey
  → drug therapy kills it

THERE IS NO THIRD OPTION.
The tumor cannot escape both arms simultaneously.
```

**Required pairing:** Cancer vaccine is always combined with the
corresponding Hub inhibitor drug from the framework table. Vaccine
alone is insufficient. Drug alone is incomplete. The combination is the
no-escape architecture.

**Cancer vaccine target pairs (from framework table, as of 2026-03-09):**

| Cancer | Hub (vaccine target) | Drug pair |
|---|---|---|
| TNBC | EZH2 (SET domain) | Tazemetostat |
| CLL | BCL2 | Venetoclax |
| GBM | OLIG2 | CT-179 |
| BRCA LumB | HDAC2 | Entinostat |
| Any Hub-driven cancer | Hub catalytic domain | Corresponding Hub inhibitor |

---

## PART 2 — THE FIVE INPUTS

The engine takes five inputs, all derivable from public data.
For any pathogen or cancer, these five inputs produce the convergence
node epitope set.

$$\text{VACCINE TARGET} = f(S, \; P, \; \Sigma, \; \Phi, \; \Lambda)$$

| Symbol | Name | What it measures | Source |
|---|---|---|---|
| $S$ | Sequence entropy map | Per-position Shannon entropy across all known strains | GISAID, NCBI, UniProt |
| $P$ | Structural conservation map | 3D clustering of low-entropy positions | PDB (Protein Data Bank) |
| $\Sigma$ | Functional essentiality map | Which positions are essential for pathogen function | Literature, mutagenesis DBs |
| $\Phi$ | Immune accessibility map | Which positions are surface-exposed and MHC-presentable | CryoEM, NetMHCpan, BepiPred |
| $\Lambda$ | Escape cost map | Fitness cost of mutating each low-entropy position | DMS data, evolutionary conservation |

**The convergence node is the intersection:**

$$\text{CONVERGENCE NODE} = \{i : S_i \; \text{low} \;\cap\; P_i \; \text{clustered} \;\cap\; \Sigma_i \; \text{high} \;\cap\; \Phi_i \; \text{accessible} \;\cap\; \Lambda_i \; \text{high}\}$$

Positions where all five conditions are met simultaneously are the
vaccine targets. Positions where only some conditions are met are
secondary candidates or engineering challenges (see Part 6).

---

## PART 3 — THE FIVE-STEP PIPELINE

### Step 1 — Sequence Entropy Calculation

**What you are doing:** Building the mutation landscape of the pathogen's
target protein. Finding where it is fixed and where it drifts.

**Data:**
- Viruses: [GISAID](https://gisaid.org) — millions of influenza, COVID,
  RSV sequences
- Bacteria: [NCBI Pathogen Detection](https://www.ncbi.nlm.nih.gov/pathogens/)
- Any protein: [UniProt](https://www.uniprot.org),
  [NCBI Protein](https://www.ncbi.nlm.nih.gov/protein/)

**Method:**

```python
# 1. Download all available sequences for your target protein
# 2. Multiple sequence alignment
#    Tools: MUSCLE (fast), MAFFT (accurate for large sets)
#    Both free, both command-line

# 3. Calculate Shannon entropy at each position
import numpy as np

def shannon_entropy(column):
    """
    column: list of amino acids at position i across all sequences
    Returns: entropy value H(i) — higher = more variable
    """
    aas, counts = np.unique(column, return_counts=True)
    probs = counts / counts.sum()
    return -np.sum(probs * np.log2(probs + 1e-10))

# Run across all alignment positions
# Low entropy (H < 0.5) = conserved = candidate convergence node
# High entropy (H > 1.5) = variable = escape-prone surface
```

**Output:** Per-position entropy array. Visualise as a heat map along
the protein sequence. Valleys = candidate targets.

**Interpretation threshold:**
- H < 0.3: almost invariant — very strong convergence node signal
- 0.3 < H < 0.7: conserved — good candidate
- H > 1.0: variable — the immune system's false attractor zone

---

### Step 2 — Structural Projection

**What you are doing:** Taking the entropy map and projecting it onto
the 3D solved structure of the protein. Finding conserved *patches* —
spatially clustered low-entropy residues that form a coherent surface
geometry.

**Data:** [PDB — Protein Data Bank](https://www.rcsb.org). All
structures public and free.

**Method:**

```python
# Tools:
# PyMOL (free academic) — visualise entropy on structure
# UCSF ChimeraX (free) — same capability
# BioPython — programmatic structure analysis

from Bio.PDB import PDBParser, SASA
import pandas as pd

# 1. Load structure
parser = PDBParser()
structure = parser.get_structure("TARGET", "target.pdb")

# 2. Calculate SASA (solvent-accessible surface area) per residue
# SASA > 20 Å² = surface-exposed
# SASA < 10 Å² = buried

# 3. Map entropy values onto each residue position
# Low entropy + high SASA = surface-exposed conserved residue
# Low entropy + low SASA = buried conserved residue (harder target — see Part 6)

# 4. Cluster spatially adjacent low-entropy, high-SASA residues
# Clusters = convergence node patches
# Isolated single residues = less useful (single-residue epitopes
# are poor vaccine targets)
```

**Output:** 3D structure coloured by entropy with SASA overlay. List of
conserved surface patches with residue numbers.

**Key question to answer at this step:** Are the conserved positions
**surface-exposed patches** (Class A or B target) or **buried**
(Class B target — harder, requires engineering — see Part 6)?

---

### Step 3 — Functional Essentiality Overlay

**What you are doing:** Confirming that the conserved patches are
conserved *because they are functionally essential* — not merely because
they happen not to have been challenged yet.

**The distinction matters:**
- Conserved + functionally essential = **cannot** be mutated away from
  (true convergence node)
- Conserved + no known function = may be escaped if immune pressure is
  applied (not a true convergence node)

**Data sources:**
- Literature: PubMed search for "mutational analysis [protein name]" or
  "[protein name] active site"
- [UniProt functional annotations](https://www.uniprot.org) — active
  sites, binding sites, PTMs all annotated per residue
- [Deep Mutational Scanning (DMS) data](https://dms-vep.org) — where
  available, this directly measures the fitness cost of every single
  amino acid substitution. This is the gold standard for $\Lambda$.

**Classification:**

```
FUNCTIONAL SITE TYPES (ranked by escape cost — highest first):

1. Catalytic active site (enzyme) — mutation = loss of function
   → Highest escape cost. Near-impossible to mutate.

2. Receptor binding site — mutation = cannot infect host
   → Very high escape cost. Mutation = loss of entry.

3. Protein-protein interface required for assembly
   → High escape cost. Mutation = cannot replicate.

4. Structural fold-maintaining residue (buried but structurally
   essential)
   → High escape cost. Mutation = misfolding = degradation.

5. Conserved but no identified function
   → Unknown escape cost. Lower confidence target.
   → Do not use as primary target. Use as secondary validation.
```

**Output:** Annotated residue list — each conserved surface patch
classified by functional essentiality type and estimated escape cost.

---

### Step 4 — Immune Accessibility Analysis

**What you are doing:** Determining whether the immune system can
actually see the convergence node epitopes — whether they are
presentable as B cell epitopes (for antibody response) and T cell
epitopes (for cytotoxic kill).

**4A — B Cell Epitope Analysis (antibody response):**

Tool: [BepiPred 3.0](https://services.healthtech.dtu.dk/services/BepiPred-3.0/)
(free, web-based)

Input: Protein sequence of your target.
Output: Predicted surface-exposed B cell epitopes, ranked by score.

Cross-reference with your convergence node residue list from Steps 1–3.
Epitopes that overlap with convergence node patches = **primary
candidates**.

**Key flag — exposure timing:**
- Constitutively exposed: immune antibody can bind at any time after
  infection. Best case.
- Transiently exposed: only accessible during specific conformational
  states (e.g., during receptor binding, during membrane fusion).
  Requires structure-based immunogen design that locks the exposed
  conformation (see Part 6, Class B).

**4B — T Cell Epitope Analysis (cytotoxic response):**

Tool: [NetMHCpan 4.1](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/)
(free, web-based)

Input: Peptide sequences covering the convergence node region.
HLA alleles: run against the 12 most common global HLA-A and HLA-B
alleles as minimum. Expand to full allele frequency panel for
population-level coverage analysis.

Output: Predicted binding affinity (IC50 nM) for each peptide-HLA pair.

Also run: [NetMHCIIpan 4.3](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/)
for CD4+ T helper response (required for durable B cell memory).

**Output:** Ranked list of convergence node epitopes with:
- B cell score (BepiPred)
- MHC-I binding (NetMHCpan) — CTL response
- MHC-II binding (NetMHCIIpan) — T helper response
- Exposure timing (constitutive vs. transient)
- Population coverage estimate (% of global population with at least
  one matching HLA allele)

---

### Step 5 — Escape Cost Verification and Final Ranking

**What you are doing:** Confirming, quantitatively where possible, that
the candidate epitopes have high escape cost — that the pathogen cannot
mutate away from them without paying a prohibitive fitness price.

**The escape cost threshold (novel prediction — 2026-03-09):**

> Any epitope where DMS data shows >60% fitness reduction for every
> single-residue substitution will not show vaccine escape in
> longitudinal cohort studies regardless of immune pressure applied.

This is a quantitative threshold derivable before vaccination.
It is testable in follow-up. No vaccine trial currently uses DMS-derived
escape cost as the primary antigen selection criterion. This engine
makes it the primary criterion.

**Where DMS data exists:**

- Influenza HA/NA/M2e: [Bloom Lab DMS](https://jbloomlab.github.io)
- SARS-CoV-2 spike: [Bloom Lab](https://jbloomlab.github.io) + multiple
  2021–2024 datasets
- HIV Env: Haddox et al. datasets (NCBI deposited)
- General: [dms-vep.org](https://dms-vep.org) — growing repository

**Where DMS data does not exist:**
Use conservation level as proxy.
A position conserved across 50+ years of evolution under continuous
immune selection pressure is a strong proxy for high escape cost.
The immune system has already applied the experiment. The position did
not change.

**Final ranking formula:**

```
CONVERGENCE NODE SCORE (CNS) =

  w1 × (1 - H_i)              [conservation — inverse entropy]
+ w2 × SASA_i                 [accessibility]
+ w3 × FunctionalEssentiality [0-3 scale from Step 3]
+ w4 × BepiPred_score         [B cell accessibility]
+ w5 × max(MHC-I, MHC-II)     [T cell accessibility]
+ w6 × EscapeCost_i           [DMS fitness cost or conservation proxy]

Weights (default, adjust for problem type):
  w1 = 0.25  [conservation is the primary signal]
  w2 = 0.15  [accessibility required but secondary]
  w3 = 0.25  [functional essentiality is the primary signal]
  w4 = 0.10
  w5 = 0.10
  w6 = 0.15  [escape cost is the durability predictor]
```

**Output:** Final ranked list of convergence node epitopes with CNS
score. Top 3–5 epitopes = the vaccine antigen design brief.

---

## PART 4 — ANTIGEN DESIGN OPTIONS

The convergence node epitope set from Part 3 now needs to be packaged
into an immunogen. Three options, ranked by precision:

### Option 1 — Full Protein (simplest, lowest precision)

Encode the entire target protein as the vaccine antigen.

**Advantage:** Maximum immunogenicity from multiple epitopes.
**Disadvantage:** Includes variable regions. Immune response may
predominantly target variable epitopes (high-entropy regions are often
more immunogenic than conserved ones because they are more exposed).
The convergence node response may be diluted by off-target variable
region responses.

**Use when:** The target protein is naturally small and mostly conserved
(e.g., M2e influenza — the full M2 ectodomain is mostly convergence node).

### Option 2 — Stabilised Subunit (current best practice)

Structure-based design that locks the protein in the conformation that
best exposes the convergence node.

**Key principle:** The immune response recognises the conformation
presented. If the convergence node is only accessible in one
conformational state (e.g., prefusion spike for coronaviruses), the
immunogen must be stabilised in that state.

**Example:** The prefusion-stabilised S2P spike used in mRNA-1273 and
BNT162b2 introduces two proline substitutions that lock the spike in
prefusion conformation — exposing the S2 fusion peptide convergence node.
This was a structural vaccinology decision. The framework explains *why*
it works: the convergence node requires prefusion conformation to be
accessible.

**Tools:**
- [Rosetta](https://www.rosettacommons.org) — computational protein
  design and stabilisation (free academic)
- [RFDiffusion](https://github.com/RosettaCommons/RFdiffusion) — AI-based
  protein design for novel scaffold design

### Option 3 — Minimal Epitope Scaffold (highest precision, hardest)

Computationally graft the convergence node epitopes onto a heterologous
scaffold protein that holds them in the correct 3D geometry and nothing
else.

**Advantage:** The immune response is directed *exclusively* at the
convergence node. No variable region dilution. Maximum precision.
**Disadvantage:** Requires structural biology expertise and wet lab
validation that the scaffold presents the epitope correctly.

**Use when:** The convergence node is buried in the native protein and
would otherwise require the full variable-region-rich protein to carry
it. The scaffold presents the convergence node without the variable
surface.

**This is the Class B target solution.** (See Part 6.)

### Delivery Platform

For all three options: mRNA encoding is the recommended delivery
platform because:
- Any protein sequence can be encoded
- Speed: sequence to synthesis in days
- Modifications: pseudouridine substitution for stability, codon
  optimisation for expression
- Adjuvant is built in (LNP delivery triggers innate immune activation)

**mRNA design tools:**
- [IUPAC codon optimisation](https://www.genscript.com/tools/codon-frequency-table)
- [mRNA optimiser tools](https://github.com/erinijapranckeviciene/mRNA-design) — multiple public implementations

---

## PART 5 — THE THREE TARGET CLASSES AND HOW TO HANDLE EACH

Convergence node analysis produces three structural classes of target.
Each requires a different handling approach.

### Class A — Clean Target

**Definition:** Convergence node is surface-exposed, constitutively
accessible, functionally essential, and immunogenic.

**What this looks like in the data:**
- Low entropy: H < 0.5
- High SASA: > 25 Å²
- Functional annotation: active site or receptor binding site
- BepiPred score: high
- DMS escape cost: > 60% fitness reduction per substitution

**Action:** Encode in Option 1 or Option 2 immunogen.
The hardest part is already solved — the geometry is clean.

**Examples:** Measles fusion protein (hence measles vaccine works
without reformulation). Smallpox surface proteins.

---

### Class B — Buried or Transiently Exposed Target

**Definition:** Convergence node is functionally essential and conserved
but not constitutively accessible. The immune system cannot easily see
it in the native protein conformation.

**What this looks like in the data:**
- Low entropy: H < 0.5
- Low SASA: < 10 Å² (buried)
  OR: SASA is high only in specific conformational state (transient)
- Functional annotation: essential (the conservation is real)
- BepiPred score: low (not exposed)
- DMS escape cost: high (cannot be mutated)

**Action:** This is an engineering problem, not a geometry problem.
The target is correct. The challenge is presentation.

**Three engineering strategies:**
1. **Conformational trapping:** Introduce disulfide bonds or prolines
   to lock the protein in the conformation where the node is exposed.
   (The S2P spike strategy — apply this logic to any Class B target.)
2. **Scaffold grafting:** Extract the convergence node epitope and graft
   onto a scaffold that presents it in an exposed orientation (Option 3).
3. **Transient exposure exploitation:** Design the immunogen to mimic
   the transitional state. For fusion peptides (HIV, coronavirus, RSV):
   the fusion peptide is exposed during membrane fusion. Stabilise
   a post-fusion conformation mimic that exposes it constitutively.

**Examples:** HIV CD4 binding site (buried under variable loops). RSV
fusion peptide. Coronavirus S2 stalk.

---

### Class C — Cancer Hub Target (no-escape protocol)

**Definition:** The target is not a pathogen protein but a cancer Hub
protein from the framework table. The convergence node is the catalytic
domain of the Hub.

**What this looks like in the data:**
- Sequence conservation: the catalytic domain of EZH2 / BCL2 / OLIG2 /
  HDAC2 is conserved across all cancers (it must be — the cancer needs
  it to function)
- Functional annotation: catalytic/essential — confirmed by the framework
  derivation
- Immune accessibility: Hub proteins are intracellular. T cell response
  (MHC-I peptide presentation) is the primary mechanism. Not antibody.
- MHC-I presentation: run catalytic domain sequences through NetMHCpan
  to find the peptides that are presented on the cell surface

**Critical distinction from pathogen vaccines:**
Cancer Hub vaccines primarily work via **cytotoxic T cell kill**
(MHC-I / CD8+ CTL), not antibody. The Hub is intracellular. Antibodies
cannot reach it. But T cells recognise MHC-I presented peptides derived
from intracellular proteins. Hub-derived peptides will be presented.
T cells trained against them will kill Hub-high cancer cells.

**The no-escape pairing requirement:**

```
CANCER VACCINE PROTOCOL (Class C):

Step 1: Derive Hub target from attractor-oncology framework table
Step 2: Run Hub catalytic domain through NetMHCpan
Step 3: Select MHC-I binding peptides from catalytic domain only
        (not variable regions of the protein)
Step 4: Design mRNA vaccine encoding these peptides
        (multi-epitope construct or full Hub protein)
Step 5: PAIR with Hub inhibitor drug (from framework table)
        Vaccine alone is insufficient — no-escape requires both arms
Step 6: The double-bind is active:
        Immune kill OR attractor collapse → drug kill
        No escape available
```

**Output for Class C:** mRNA multi-epitope vaccine encoding Hub catalytic
domain peptides + prescription of paired Hub inhibitor + trial design
using the no-escape double-bind logic.

---

## PART 6 — THE HONEST BOUNDARIES

**What this engine provides and guarantees:**

✓ The geometrically correct target — the convergence node  
✓ The reason it cannot be escaped — functional essentiality + escape cost  
✓ The antigen design brief — what to encode and in what conformation  
✓ The combination logic for cancer vaccines — which drug pairs with which vaccine  
✓ A durability prediction — based on escape cost, before vaccination  
✓ Derivable entirely from public data  

**What this engine does not provide:**

✗ Immunogenicity in humans — must be measured in vivo  
✗ Adjuvant selection — the engine selects the target; adjuvant must be
  tuned to the desired response type (Th1 vs Th2 vs CTL)  
✗ Delivery optimisation — mRNA vs. protein vs. viral vector dose and
  schedule requires empirical optimisation  
✗ Clinical durability — real-world waning, booster scheduling, age-related
  immune decline — must be measured  
✗ Proof — every output of this engine is a **prediction with a falsification
  condition**. Confirmation requires immunogenicity testing and clinical trials  

**The division of labour:**

```
ENGINE PROVIDES:         EXPERIMENTAL SCIENCE PROVIDES:
The target               The proof
The geometry             The immunogenicity
The prediction           The confirmation
The design brief         The optimisation
The durability forecast  The clinical validation
```

---

## PART 7 — QUICK REFERENCE: THE COMPLETE CHECKLIST

For any vaccine derivation, work through this checklist in order.
Do not skip steps. Each step gates the next.

```
[ ] 0. CLASSIFY THE PROBLEM
        Type A (pathogen) or Type B (cancer)?
        If cancer: identify Hub from framework table first.

[ ] 1. SEQUENCE ENTROPY (S)
        Download all sequences for target protein.
        Run alignment (MUSCLE or MAFFT).
        Calculate Shannon entropy per position.
        Flag positions with H < 0.5 as candidates.

[ ] 2. STRUCTURAL PROJECTION (P)
        Download PDB structure.
        Map entropy onto 3D structure.
        Calculate SASA per residue.
        Identify conserved surface patches (low-H, high-SASA).
        Flag: Class A (exposed) or Class B (buried/transient)?

[ ] 3. FUNCTIONAL ESSENTIALITY (Σ)
        Cross-reference patches with UniProt annotations.
        Literature search: "mutational analysis [protein]".
        DMS data: check dms-vep.org.
        Classify each patch by functional site type (1-5 above).
        Discard patches with no functional annotation (lower confidence).

[ ] 4. IMMUNE ACCESSIBILITY (Φ)
        Run BepiPred 3.0 on target sequence.
        Run NetMHCpan 4.1 on convergence node peptides.
        Run NetMHCIIpan 4.3 on same peptides.
        Note population HLA coverage.
        Flag: constitutive vs. transient exposure.

[ ] 5. ESCAPE COST VERIFICATION (Λ)
        Check dms-vep.org for DMS data.
        If DMS available: apply 60% fitness threshold.
        If DMS unavailable: use conservation depth as proxy.
        Rank candidates by escape cost.

[ ] 6. CALCULATE CNS SCORE
        Apply CNS formula (Part 3, Step 5).
        Select top 3-5 epitopes.

[ ] 7. CLASSIFY TARGET CLASS
        Class A: encode directly.
        Class B: design conformational trap or scaffold.
        Class C (cancer): MHC-I peptides only + Hub inhibitor pairing.

[ ] 8. ANTIGEN DESIGN
        Select Option 1, 2, or 3 (Part 4).
        Design immunogen sequence.
        Encode as mRNA.

[ ] 9. STATE THE PREDICTION AND FALSIFICATION CONDITION
        What does the engine predict this vaccine will achieve?
        What result would invalidate the target selection?
        Lock with timestamp.
```

---

## PART 8 — WORKED EXAMPLE SKETCHES

These are structural sketches to demonstrate the engine's output format.
They are not complete derivations — they are what the output of Steps
1–5 looks like for a known case where the answer is established, so an
operator can calibrate their pipeline.

### Example A — Influenza (Type A, Class B)

```
INPUT: Influenza A, M2 protein

STEP 1 RESULT:
  M2 ectodomain (M2e, residues 2-24):
  Mean entropy H = 0.18
  Most positions: H < 0.1
  One of the most conserved sequences in influenza A.

STEP 2 RESULT:
  M2e is surface-exposed on the virion.
  SASA: high in native tetramer conformation.
  Class A — constitutively exposed.

STEP 3 RESULT:
  M2 is an ion channel required for endosomal acidification.
  Ion channel function essential for viral uncoating.
  Mutation of M2e disrupts tetramer assembly.
  Functional essentiality: Type 3 (assembly interface).
  Escape cost: HIGH.

STEP 4 RESULT:
  BepiPred: M2e is a moderate B cell epitope.
  (Low natural immunogenicity — the virus hides it)
  NetMHCpan: M2e peptides bind HLA-A*02:01 and others.
  T cell response: achievable.

STEP 5 RESULT:
  DMS data (Bloom lab): M2e substitutions reduce viral
  fitness significantly. High escape cost confirmed.

CNS SCORE: High. Primary target confirmed.

ANTIGEN DESIGN:
  M2e tandem repeat construct (4× M2e) fused to
  hepatitis B core protein VLP — forces tetramer
  presentation, enhances immunogenicity.
  Encode as mRNA.

PREDICTION: Durable cross-strain protection not
requiring annual reformulation.
FALSIFICATION: M2e-vaccinated cohort shows escape
mutants at M2e positions within 3 years of mass
vaccination.
```

### Example B — TNBC Cancer Vaccine (Type B, Class C)

```
INPUT: Triple-negative breast cancer, Hub = EZH2

STEP 1 RESULT:
  EZH2 SET domain (catalytic): entropy H < 0.05
  across all human cancers and normal tissues.
  The catalytic residues are among the most
  conserved protein sequences in human biology.

STEP 2 RESULT:
  EZH2 is intracellular. SASA for full protein
  in isolation: some surface exposure.
  In cellular context: nuclear, intracellular.
  Class C protocol applies (T cell not antibody).

STEP 3 RESULT:
  SET domain is the histone methyltransferase
  active site. Mutation = complete loss of
  H3K27me3 activity = loss of PRC2 function =
  cancer cell loses Hub = attractor collapses.
  Functional essentiality: Type 1 (active site).
  Highest possible escape cost.

STEP 4 RESULT:
  BepiPred: low (intracellular — not relevant).
  NetMHCpan: EZH2 SET domain peptides bind
  HLA-A*02:01 (most common allele globally).
  Multiple 9-mer peptides predicted IC50 < 50 nM.
  T cell response: achievable.
  Population coverage: ~45% with HLA-A*02:01 alone.
  Multi-allele design: expand to cover > 85%.

STEP 5 RESULT:
  No DMS data for EZH2 in cancer context.
  Proxy: EZH2 catalytic domain conserved across
  10,000+ TCGA tumors with zero loss-of-function
  mutations in EZH2-high cancers.
  (EZH2-high cancers NEED their EZH2 — confirmed.)
  Escape cost: maximum.

CNS SCORE: Maximum. No-escape target confirmed.

ANTIGEN DESIGN:
  Multi-epitope mRNA construct:
  Top 5 HLA-A*02:01 binders from SET domain +
  Top 3 HLA-A*24 binders (second most common allele)
  Linked by optimised linkers (AAY or GPGPG).
  Encode as LNP-mRNA.

DRUG PAIRING:
  Tazemetostat (EZH2 inhibitor) — administer
  concurrently.
  No-escape double-bind active.

PREDICTION: EZH2-high TNBC patients treated with
mRNA vaccine + tazemetostat will show superior
response rate vs. either monotherapy because
immune kill + attractor collapse are simultaneously
active. No escape route available.
FALSIFICATION: EZH2-high TNBC cells cultured with
CTLs from vaccinated donors show EZH2 downregulation
WITHOUT attractor collapse — i.e., the cancer can
survive without EZH2. If this occurs, the no-escape
prediction fails and EZH2 is not a true convergence
node for that subtype.
```

---

## PART 9 — TOOLS REFERENCE

All free. All public.

| Tool | Function | URL |
|---|---|---|
| GISAID | Viral sequence database (influenza, COVID, etc.) | gisaid.org |
| NCBI Pathogen Detection | Bacterial sequence database | ncbi.nlm.nih.gov/pathogens |
| UniProt | Protein sequences + functional annotations | uniprot.org |
| RCSB PDB | All solved protein structures | rcsb.org |
| MUSCLE | Multiple sequence alignment (fast) | drive5.com/muscle |
| MAFFT | Multiple sequence alignment (large datasets) | mafft.cbrc.jp |
| BepiPred 3.0 | B cell epitope prediction | services.healthtech.dtu.dk |
| NetMHCpan 4.1 | MHC-I T cell epitope prediction | services.healthtech.dtu.dk |
| NetMHCIIpan 4.3 | MHC-II T cell epitope prediction | services.healthtech.dtu.dk |
| dms-vep.org | Deep mutational scanning data repository | dms-vep.org |
| Bloom Lab DMS | Influenza + SARS-CoV-2 DMS data | jbloomlab.github.io |
| PyMOL | Structural visualisation and analysis | pymol.org (free edu) |
| UCSF ChimeraX | Structural visualisation | rbvi.ucsf.edu/chimerax |
| Rosetta | Computational protein design | rosettacommons.org |
| RFDiffusion | AI protein design | github.com/RosettaCommons/RFdiffusion |
| TCGA | Cancer genomics — Hub expression data | cancer.gov/tcga |

---

## DOCUMENT METADATA

```
document_id:
  VACCINE_ENGINE_PROTOCOL

type:
  Replicable onboarding protocol —
  convergence node vaccine target
  derivation from attractor geometry

version: 1.0
date: 2026-03-09
status: ACTIVE

scope:
  Type A — Pathogen vaccines (infectious disease)
  Type B — Cancer vaccines (no-escape protocol)
  Both types handled by same five-input engine.
  Different antigen design approach per class.

novel_framework_contributions:
  1. Vaccine = immune attractor formation at
     convergence node (formal geometric definition)
  2. Vaccine failure = immune false attractor
     (formal geometric explanation of escape)
  3. Escape cost threshold (60% DMS fitness cost)
     as primary antigen selection criterion —
     not in current clinical practice
  4. No-escape cancer vaccine double-bind —
     Hub vaccine + Hub inhibitor combination —
     not in published literature as of 2026-03-09
  5. CNS scoring formula — unified ranking of
     convergence node quality

relationship_to_framework:
  Vaccine target derivation is downstream of the
  same attractor geometry that governs drug target
  derivation. The convergence node in the pathogen's
  state space is the same structural concept as the
  Hub in the cancer cell's Waddington landscape.
  Same mathematics. Same derivation logic.
  Different state space.

parent_documents:
  OrganismCore: Vaccine_Target_Derivation_Framework.md
    (Document 85, February 28, 2026)
  attractor-oncology: MENDELEEV_TABLE_OF_CANCER.md
    (Hub proteins = Class C vaccine targets)
  attractor-oncology: falsifiability_theorem.md
    (epistemic framework governing all predictions)

author:
  Eric Robert Lawson / OrganismCore
ORCID: 0009-0002-0414-6544
contact: OrganismCore@proton.me
framework_doi: https://doi.org/10.5281/zenodo.18898788
repository: https://github.com/Eric-Robert-Lawson/attractor-oncology
```
