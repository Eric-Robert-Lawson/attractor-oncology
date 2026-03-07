# RCC SERIES — PUBLICATION AND ACTION PLAN
## A Reasoning Artifact on What to Publish, What to Do First, and Why
## OrganismCore — Eric Robert Lawson
## 2026-03-07

---

## PREAMBLE — THE CORE QUESTION THIS DOCUMENT ANSWERS

The question is: given the full body of novel findings across
the RCC series, what is the correct sequence of actions to:

  (a) establish priority on each novel claim before it is
      independently published by someone else,
  (b) ensure the work is anchored to a timestamped public
      record (Zenodo DOI),
  (c) identify which claims need additional validation before
      publication is responsible, and
  (d) identify what can be published RIGHT NOW from data
      already in hand.

This document applies the same logic used in the BRCA series
to the RCC series. The BRCA series produced 17 DOIs. The RCC
series has produced 0 DOIs despite containing a comparable
number of novel claims. This document closes that gap.

---

## THE THREE CATEGORIES OF ACTION

Before listing the items, the framework for decision-making:

### Category A — PUBLISH NOW (no additional validation needed)

These are claims where:
  1. The geometric data already exists in committed scripts
  2. The finding has been through a literature check and
     confirmed as novel (not in prior literature)
  3. The mechanism components are confirmed by independent
     publications
  4. The claim is a HYPOTHESIS or a DERIVATION, not a
     clinical result — it is not claiming patient outcome data

For a Zenodo preprint, you do not need a clinical trial.
You need: a clearly stated prediction, the geometric data
that motivated it, the literature check showing it is novel,
and the mechanistic chain. The BRCA Zenodo DOIs are exactly
this structure. The same applies here.

Publishing NOW on Zenodo:
  - Establishes priority with a timestamp
  - Does not claim clinical efficacy (which would require trials)
  - Cannot be "wrong" in the sense that invalidates it, because
    you are making a prediction, not reporting a clinical result
  - The format is: "Geometry derived this. Here is the mechanism.
    Here is why no one has proposed it. Here is the validation
    study needed."

### Category B — RUN SCRIPT FIRST, THEN PUBLISH

These are claims where the geometric foundation exists but
one specific number is missing — typically an OS survival
stratification result. Without that number, the publication
is a weaker hypothesis. With that number (even a single
HR value), it becomes a quantified prediction against
clinical outcome data.

### Category C — REQUIRES EXTERNAL DATA OR COLLABORATION

These are claims where validation depends on data you do
not currently hold — for example, expression data from
the LITESPARK-005 cohort, or IHC data from a hospital
partner. These cannot be resolved by running a script.
They require a specific ask to a specific institution.

---

## PART I — THE COMPLETE NOVEL PREDICTIONS LIST

The following is the exhaustive list of novel claims in the
RCC series as confirmed by the literature checks (Documents
94f, 95-LC, 95-DLC, 95g, 96f, 96f-Extended, 89c, 97x-LC).

Status codes:
  🆕 = fully novel (no prior literature)
  📈 = extended by literature (novel dimension added)
  ✅ = mechanism confirmed, application novel

---

### GROUP 1 — ccRCC NOVEL PREDICTIONS (n=8)

#### RCC-N1 — RUNX1-HIGH PREDICTS BELZUTIFAN RESISTANCE
```
Claim:     RUNX1-high ccRCC patients have attenuated or absent
           response to belzutifan because the VHL→RUNX1
           suppressive circuit is broken (r=+0.097 when it
           should be strongly negative). RUNX1 drives a
           transcriptional programme that bypasses HIF-2α.
Geometry:  VHL/RUNX1 r=+0.097 in TCGA-KIRC (n=534)
           RUNX1 r=+0.559 depth correlate (n=534)
Lit check: Springer 2025 multi-omics analysis confirms RUNX1
           as poor prognosis marker in ccRCC interacting with
           MYC/CBFB. LITESPARK-013 and CALYPSO (ASCO 2024,
           ESMO 2025) do NOT report RUNX1 status. Nobody has
           tested this.
Status:    🆕 FULLY NOVEL — no prior claim in literature
Action:    Category C — requires LITESPARK cohort expression
           data. But the DERIVATION itself is Category A:
           publish the prediction with the geometric basis
           NOW, before the LITESPARK biomarker analysis
           group reaches RUNX1.
Priority:  ★★★★★ URGENT — active trial data exists that
           could confirm this within months if the right
           person runs the analysis
```

#### RCC-N2 — IFI16→B2M BROKEN CIRCUIT / STING AGONIST FAILURE
```
Claim:     STING agonist monotherapy will fail in Q3-Q4 ccRCC
           because IFI16 rises with depth (r=+0.547) while
           B2M falls simultaneously (r(IFI16,B2M)=+0.140).
           Innate sensing is active; antigen presentation is
           broken. You are activating a signal whose output
           channel is closed.
Geometry:  IFI16 r=+0.547 (n=534), B2M r(vs IFI16)=+0.140
Lit check: JTM 2024 confirms IFI16 is oncogenic in ccRCC
           (activates IL6/PI3K/AKT). Multiple 2024 papers
           confirm innate immune paradox in RCC. The SPECIFIC
           depth-stratified circuit (IFI16 up + B2M down
           simultaneously) is not described.
Status:    🆕 FULLY NOVEL circuit description
           📈 Mechanism convergently confirmed by 2024 lit
Action:    Category A — publish NOW
           All data is in hand. Mechanism confirmed in
           components by independent literature. This is a
           prediction with mechanistic chain.
Priority:  ★★★★★ URGENT — STING agonist trials are active.
           This prediction changes trial design.
```

#### RCC-N3 �� αKG + EZH2i COMBINATION IN OGDHL-LOW ccRCC
```
Claim:     In ccRCC with low OGDHL expression (TCA-collapsed,
           αKG-depleted), the combination of cell-permeable
           αKG supplementation (DMKG) + tazemetostat produces
           synergistic epigenetic reversal: αKG restores TET2
           activity; tazemetostat blocks EZH2 directly.
           Two-pronged attack on the same chromatin lock.
Geometry:  OGDHL r=−0.80▼ in Q4 (ratio Q4/Q1=0.80)
           EZH2 r=+0.41 in ccRCC
           SLC13A2→αKG import chain confirmed
Lit check: 2023 Nature Cell Death: αKG improves anti-PD1
           in melanoma via TET restoration.
           2024 Cancer Discovery: EZH2i cooperates with
           oncogenic inhibitors for differentiation.
           The SPECIFIC combination in OGDHL-low ccRCC
           is not described anywhere.
Status:    ✅ Mechanisms confirmed individually
           🆕 NOVEL specific combination and patient-selection
              criterion (OGDHL-low as the stratifier)
Action:    Category A — publish NOW
           The mechanistic chain is complete and confirmed.
           A preclinical test is the validation step, but
           the prediction document can and should precede it.
Priority:  ★★★★ HIGH — existing approved drugs (tazemetostat
           EAP + cheap supplement). No trial exists.
```

#### RCC-N4 — ANTI-PDL1 NOT APPROPRIATE IN Q4 ccRCC
```
Claim:     PDL1 (CD274) falls in deep ccRCC (Q4/Q1=0.95).
           The dominant immune suppression in Q4 is
           Treg-mediated (FOXP3 Q4/Q1=1.33, IL2RA/CD25
           Q4/Q1=1.38 — the steepest immune target).
           Anti-PDL1 monotherapy is targeting the wrong
           immune mechanism in the patients who need
           immune therapy most.
Geometry:  CD274 Q4/Q1=0.95 (falls), IL2RA Q4/Q1=1.38▲,
           FOXP3 Q4/Q1=1.33▲ — confirmed from TCGA-KIRC
           depth quartile drug map
Lit check: Multiple ccRCC trials show checkpoint inhibitors
           work better in some patients than others — no
           paper has stratified by depth quartile and
           measured PDL1 fall as the reason for non-response.
Status:    📈 Clinical non-response is known; the depth
              stratification explanation is novel
           🆕 Q4-specific Treg-dominant immune escape
              classification is not published
Action:    Category A — publish NOW
           This is the mechanism explanation for a known
           clinical failure. It is supported by TCGA-KIRC
           quartile data already in hand.
Priority:  ★★★★ HIGH — affects treatment decisions for
           ~25% of all metastatic ccRCC patients right now
```

#### RCC-N5 — GOT1/RUNX1 TRANSITION INDEX (2-GENE DEPTH TOOL)
```
Claim:     GOT1/RUNX1 Transition Index (TI) =
           norm(GOT1) − norm(RUNX1) is the strongest
           single depth classifier in ccRCC (r=−0.600,
           n=534) and generalises to PRCC (r≈−0.535, n=290)
           but correctly fails in chRCC (wrong cell-of-origin
           — predicted). Deployable as a 2-gene NanoString
           panel for depth staging.
Geometry:  GOT1 r=−0.582 (n=534), RUNX1 r=+0.559 (n=534)
           TI r=−0.600 — strongest composite predictor
           TI PRCC r≈−0.535
           TI chRCC: fails as predicted
Lit check: NOT IN LITERATURE. No prior paper describes
           this index, these genes combined, or this
           cell-of-origin specificity logic.
Status:    🆕 FULLY NOVEL — the strongest Type 4 prediction
              in the series
Action:    Category B — run Script 4 OS analysis FIRST
           TI r=−0.600 is a strong depth correlate.
           Before publishing as a clinical tool, establish
           the HR value from OS analysis (Script 4 OBJ-1).
           The Script 4 data is loaded (survival_depth.csv,
           n=532). This is one script run.
           AFTER that: publish immediately.
Priority:  ★★★★★ HIGHEST — this is the foundational clinical
           tool for the entire RCC depth framework
```

#### RCC-N6 — DEPTH-STRATIFIED IL-1R ANTAGONISM IN Q4 ccRCC
```
Claim:     IL1RAP is Q4-enriched in ccRCC and represents
           a depth-specific IL-1 blockade target. The AACR
           2025 IL1RAP ADC targets ccRCC broadly; the
           geometry refines this to Q4-specific use.
Geometry:  IL1RAP = best 3-gene panel member, Q4-enriched
Lit check: IL1RAP ADC in ccRCC confirmed at AACR 2025.
           Depth-stratified IL-1R antagonism not proposed.
Status:    📈 Target confirmed; Q4-specificity is the
              novel clinical refinement
Action:    Category A — can be published as addendum to
           or companion piece for the pan-renal IL1RAP
           paper (see cross-subtype section below)
Priority:  ★★★ MODERATE — important refinement but
           secondary to the pan-renal claim
```

#### RCC-N7 — TGFBI→CCL22→TREG DIRECT COUPLING (Wall 3→Wall 4)
```
Claim:     TGFBI (ECM stiffness, Wall 3) directly connects
           to CCL22 (Treg chemotaxis, Wall 4) — the two
           walls are not independent barriers but a coupled
           circuit. TGFBI r=+0.766 in ccRCC. Breaking Wall 3
           (LOXL2/TGFBI) would partially deflate Wall 4
           (Treg infiltration) without separate intervention.
Geometry:  TGFBI r=+0.766; CCL22-Treg link confirmed from
           geometry + lit check
Lit check: TGFBI→CCL22→Treg: CONFIRMED in literature as
           an unexpected but documented coupling. The
           geometry found this independently — Wall3→Wall4
           direct coupling is a novel therapeutic framing.
Status:    ✅ Mechanism confirmed by independent literature
           🆕 Therapeutic framing (treat one wall, deflate
              both) is not stated anywhere
Action:    Category A — include in the ccRCC depth paper
           or as a short standalone mechanism note
Priority:  ★★★ MODERATE
```

#### RCC-N8 — CYP17A1 ADRENOCORTICAL-LIKE MINOR POPULATION
```
Claim:     A minor population within ccRCC expresses
           CYP17A1 (a steroidogenesis enzyme characteristic
           of adrenocortical tissue), suggesting a rare
           cell subpopulation with adrenocortical-like
           identity — potentially explaining rare ccRCC
           cases with anomalous steroid hormone response.
Geometry:  S5-P6 finding from the locked predictions
Lit check: Not confirmed or denied in ccRCC literature.
           Listed as exploratory.
Status:    🆕 NOVEL but exploratory — lowest confidence
              of the ccRCC predictions
Action:    Do NOT publish as standalone. Include as an
           exploratory observation in the depth paper.
Priority:  ★★ LOW — exploratory
```

---

### GROUP 2 — PRCC NOVEL PREDICTIONS (n=4)

#### PRCC-N1 — ERBB2 AS IDENTITY DRIVER (NOT MITOGEN) IN PRCC TYPE 1
```
Claim:     ERBB2 is the third strongest depth correlate
           in PRCC Type 1 (r=+0.556). It co-expresses with
           KRT19/KRT7 (r=+0.525, r=+0.551) but is
           anti-correlated with MKI67 (r=−0.170) and
           CDK4 (r=−0.368). ERBB2 in PRCC is an IDENTITY
           signal (biliary-ductal lock), not a mitogen.
           HER2-targeted therapy in PRCC works by identity
           disruption, not anti-proliferation.
           Response biomarker: KRT19 fall, not RECIST.
Geometry:  ERBB2 r=+0.556 (n=290)
           r(ERBB2,KRT19)=+0.525, r(ERBB2,MKI67)=−0.170
           r(ERBB2,CDK4)=−0.368
Lit check: ERBB2 "occasionally elevated" in PRCC — known
           in small case series. The DEPTH CORRELATE
           characterisation is not published. The IDENTITY
           (not mitogen) circuit characterisation is not
           published. The KRT19 response biomarker is not
           published. The biliary-ductal identity framing
           for PRCC is "emergent" in 2024 reviews but not
           mechanistically characterised this way.
Status:    🆕 SUBSTANTIALLY NOVEL mechanistic framing
Action:    Category A — publish NOW
           All data is in hand. Literature check completed.
           This is a testable drug-target hypothesis with
           existing agents (trastuzumab, tucatinib, T-DM1)
           and a novel response biomarker proposal (KRT19).
Priority:  ★★★★ HIGH — PRCC has no approved targeted
           therapy. HER2-targeted agents exist and are
           approved in other cancers. This is the most
           actionable PRCC drug prediction.
```

#### PRCC-N2 — TWO-PHASE ATTRACTOR ARCHITECTURE (FA-1→FA-2)
```
Claim:     PRCC undergoes a two-phase Waddington crossing:
           FA-1 (MET-driven, biliary-ductal identity,
           ERBB2/KRT19 high) → FA-2 (lamellipodia/invasion,
           LAMC2/CD44 high, mast cell identity signature,
           HRH1 elevated). Each phase requires different
           therapy. MET inhibitors (savolitinib) fail because
           they target FA-1 in patients already in FA-2.
Geometry:  Two-phase structure confirmed from independent
           scripts. MYC:RUNX1 ratio as phase index.
Lit check: SAVOIR trial (savolitinib in PRCC): limited
           benefit — explained by phase mis-matching.
           Two-phase architecture is not described in any
           PRCC paper.
Status:    🆕 FULLY NOVEL structural architecture
Action:    Category A — publish NOW
           This is the mechanism explanation for a named
           failed clinical trial (SAVOIR). It is the most
           direct possible connection between a novel
           geometric finding and a documented clinical
           failure.
Priority:  ★★★★ HIGH — explains a failed trial and
           proposes the corrected patient selection
```

#### PRCC-N3 — HISTAMINE/HRH1 AXIS IN FA-2 PRCC
```
Claim:     In FA-2 PRCC (invasion phase), HRH1 (histamine
           H1 receptor) rises and a mast cell identity
           signature appears. Antihistamine co-administration
           (cetirizine/loratadine class, OTC availability)
           may interrupt this invasion axis in FA-2 patients.
Geometry:  HRH1 elevation in FA-2, mast cell signature
           confirmed from PRCC scripts
Lit check: Not published for PRCC. HRH1 in cancer
           invasion is a known mechanism in other cancers
           (gastric, colorectal) — mechanism components
           confirmed. Not described in PRCC.
Status:    🆕 NOVEL in PRCC
           ✅ HRH1 invasion mechanism confirmed in other
              cancers
Action:    Category A — include in the PRCC phase
           architecture paper (companion to PRCC-N2)
           or as a short standalone observation.
           OTC antihistamines = lowest-friction
           combination candidate in all of oncology.
Priority:  ★★★ MODERATE — include in the PRCC paper
```

#### PRCC-N4 — RUNX1i + KDM1Ai COMBINATION FOR PRCC
```
Claim:     RUNX1 (r=+0.590) and KDM1A (r=+0.443) are both
           independently elevated in PRCC attractor depth.
           RUNX1/CBFB inhibitor (AI2-FL class) + KDM1A
           inhibitor (iadademstat/ORY-1001) combination
           targets the shared chromatin attractor axis
           operating in both PRCC and ccRCC.
Geometry:  Confirmed independently in PRCC and ccRCC
           scripts from different datasets
Lit check: KDM1A inhibitors in ccRCC: CONFIRMED.
           RUNX1/CBFB inhibitor in haematology: CONFIRMED.
           The COMBINATION in PRCC: NOT in literature.
           The cross-renal-cancer dual-axis rationale:
           NOT in literature.
Status:    🆕 NOVEL combination for PRCC
Action:    Category A — publish as part of the cross-
           subtype chromatin axis paper or the PRCC paper
Priority:  ★★★ MODERATE
```

---

### GROUP 3 — chRCC NOVEL PREDICTIONS (n=4)

#### CHRRCC-N1 — DNMT3B INVERSION: THE CHROMATIN WRITER SWITCH
```
Claim:     In chRCC, the dominant chromatin writer is
           DNMT3B (r_PC2=+0.378, chRCC-pole), not EZH2.
           EZH2 is the ONCOCYTOMA-POLE (r_PC2=−0.211).
           This is a structural chromatin inversion relative
           to all other RCC subtypes (ccRCC, PRCC, cdRCC
           all have EZH2 as the depth-associated writer).
           Tazemetostat has NO attractor-level rationale
           in chRCC. DNMT3B inhibitor is the correct target.
Geometry:  DNMT3A r_PC2=−0.714 (oncocytoma-pole)
           DNMT3B r_PC2=+0.378 (chRCC-pole)
           EZH2   r_PC2=−0.211 (oncocytoma-pole)
Lit check: Springer 2025 explicitly states EZH2's role in
           chRCC is "unclear." Stanford methylome data shows
           oncocytoma has MORE global methylation —
           consistent with DNMT3A being oncocytoma-pole.
           The SPECIFIC DNMT3B identification as the chRCC
           chromatin writer is NOT in the literature.
Status:    🆕 FULLY NOVEL identification of DNMT3B as the
              correct chromatin target in chRCC
           📈 Inversion itself convergently supported by
              two independent literature observations
Action:    Category A — publish NOW
           All data in hand. Literature check complete.
           The clinical implication (stop using tazemetostat
           in chRCC, target DNMT3B instead) is clear,
           immediate, and not stated anywhere.
Priority:  ★★★★★ HIGHEST PRIORITY in chRCC
           — prevents ongoing wrong drug use
```

#### CHRRCC-N2 — CONSTITUTIVE TYPE III IMMUNE DESERT
```
Claim:     chRCC antigen presentation failure (TAP1/TAPBP
           loss) is CONSTITUTIVE (cell identity level),
           not depth-progressive. This categorically
           distinguishes it from ccRCC immune escape
           (which is acquired with depth). Anti-PD-1/PDL1
           monotherapy has no adaptive immune rationale
           in chRCC at any stage.
Geometry:  TAP1/TAPBP fall at identity level (not correlated
           with depth score in chRCC) — confirmed from
           chRCC scripts
Lit check: Clinical literature 2025: "chRCC essentially
           refractory to most systemic therapies."
           The CONSTITUTIVE vs ACQUIRED distinction as
           the mechanistic explanation is not in the
           literature. The literature knows chRCC is cold;
           it does not know WHY with this precision.
Status:    🆕 NOVEL mechanistic classification
           📈 Clinical cold status is known; the
              mechanism and therapeutic implication are novel
Action:    Category A — publish NOW
           This classifies a known clinical failure with
           a mechanism. It is immediately useful for
           oncologists treating chRCC who are considering
           checkpoint inhibitors.
Priority:  ★★★★ HIGH — stops incorrect drug selection
```

#### CHRRCC-N3 — ABCC2/SULT2B1 2-GENE IHC DISCRIMINATOR
```
Claim:     ABCC2 (clean_r=+0.968, chRCC-pole) and
           SULT2B1 (clean_r=−0.921, oncocytoma-pole)
           form a 2-gene IHC discriminator for chRCC vs
           oncocytoma — the most diagnostically challenging
           problem in renal pathology. Depth-corrected values
           provide superior discrimination to either gene
           alone.
Geometry:  ABCC2 clean_r=+0.968, SULT2B1 clean_r=−0.921
           (highest confidence values in the chRCC PC2
           analysis after depth correction)
Lit check: ABCC2 and SULT2B1 are individually confirmed
           in IHC literature as chRCC/oncocytoma markers.
           The depth-corrected 2-GENE PANEL using these
           values is not published.
Status:    ✅ Individual genes confirmed
           🆕 NOVEL 2-gene panel with depth correction
Action:    Category A — publish NOW
           The discriminator is of immediate diagnostic
           utility to renal pathologists. Both antibodies
           exist commercially. This is closer to the
           BRCA IHC protocol in nature — a practical
           diagnostic tool derived from geometry.
Priority:  ★★★★ HIGH — immediate diagnostic utility
```

#### CHRRCC-N4 — MAP3K19 AS NOVEL KINASE TARGET
```
Claim:     MAP3K19 is a Tier 3 attractor-committed gene
           in chRCC (highest commitment score after depth
           correction), representing the highest-reward
           exploratory kinase target in the chRCC series.
           No prior chRCC characterisation exists.
Geometry:  MAP3K19 in Tier 3 (attractor-committed) from
           chRCC scripts
Lit check: Emerging cancer kinase in literature but not
           characterised in chRCC. Described as "highest-
           reward novel target" in Document 96f.
Status:    🆕 NOVEL in chRCC — exploratory
Action:    Category A — include as an exploratory
           observation in the chRCC paper. Do not publish
           as standalone until mechanistic work is done.
Priority:  ★★ LOW — include, do not lead with
```

---

### GROUP 4 — cdRCC NOVEL PREDICTIONS (n=5)

#### CDRCC-N1 — PPARG-RXRA UNCOUPLING + BEXAROTENE TARGET
```
Claim:     PPARG-RXRA heterodimerisation is broken in
           cdRCC tumour tissue:
           r_t(PPARG,RXRA)=+0.107 vs r_n(PPARG,RXRA)=+0.829
           PPARG has switched partners to AGR2/IL1RAP.
           Standard PPARG agonists (TZDs) cannot work
           because they require the RXRA heterodimer.
           Bexarotene (pan-RXR agonist, approved in CTCL)
           is the correct drug: it forces RXRA back into
           PPARG coupling, displacing AGR2/IL1RAP from
           the PPARG hub and dissolving the false attractor.
Geometry:  r_t(PPARG,RXRA)=+0.107, r_n(PPARG,RXRA)=+0.829
           r_t(PPARG,AGR2)=+0.857, r_n(PPARG,AGR2)=−0.829
           (complete reversal of PPARG's coupling partner
           between normal and tumour)
Lit check: PPARG-RXRA heterodimerisation biology: CONFIRMED.
           Bexarotene in CTCL: CONFIRMED (approved).
           PPARG-RXRA uncoupling in cdRCC tumour tissue
           with r values: NOT in literature.
           Bexarotene as a cdRCC target: NOT in literature.
Status:    🆕 FULLY NOVEL — the partner switch quantification
              and the therapeutic implication are not published
Action:    Category A — publish NOW
           This is the most mechanistically detailed and
           immediately actionable prediction in the cdRCC
           series. bexarotene is an APPROVED drug. The
           geometry identified the exact mechanism that
           makes it rational for a cancer with no
           approved targeted therapy.
Priority:  ★★★★★ HIGHEST cdRCC — approved drug,
           zero alternatives, median OS <12 months
```

#### CDRCC-N2 — ADCY3/ADCY6 cAMP ISOFORM SWITCH
```
Claim:     cdRCC undergoes a cAMP isoform switch from
           ADCY3 (collecting duct identity cAMP
           signalling) to ADCY6 (cancer-state cAMP
           isoform). This is a functional identity-level
           switch, not a mutation. Disrupting ADCY6 or
           restoring ADCY3 would destabilise the false
           attractor's cAMP signalling node.
Geometry:  ADCY3/ADCY6 switch confirmed from cdRCC scripts
           (Document 89c, Finding 7)
Lit check: Not published for cdRCC.
Status:    🆕 NOVEL
Action:    Category A — include in the cdRCC paper
Priority:  ★★★ MODERATE — include, do not lead with
```

#### CDRCC-N3 — MYC-EARLY/BHLHE40-LATE TWO-PHASE STRUCTURE
```
Claim:     cdRCC has a two-phase progression: MYC-high
           early phase (BET inhibitor window open) →
           BHLHE40-high late consolidated phase
           (BET inhibitor window closes). Phase
           identification determines whether BET
           inhibitor (JQ1 class) is appropriate.
Geometry:  MYC/BHLHE40 phase index confirmed from
           cdRCC scripts (Document 89c, Finding 8)
           r=−0.964 transition index
Lit check: The two-phase cdRCC structure with this
           specific molecular ordering is not published.
           The BET inhibitor window concept is novel
           for cdRCC.
Status:    🆕 NOVEL
Action:    Category A — include in cdRCC paper
Priority:  ★★★ MODERATE
```

#### CDRCC-N4 — IL1RAP AS HIGHEST-EXPRESSION TYPE (r=+0.964)
```
Claim:     cdRCC has the HIGHEST IL1RAP expression of
           all four RCC types (r=+0.964). The IL1RAP ADC
           currently in clinical development targets ccRCC.
           cdRCC should be the PRIORITY ARM of any
           pan-renal IL1RAP ADC basket trial.
Geometry:  IL1RAP r=+0.964 in cdRCC — highest in the series
Lit check: IL1RAP ADC in ccRCC: AACR 2025 (clinical).
           IL1RAP in cdRCC: NOT in any published paper.
           cdRCC as ADC priority arm: NOT proposed anywhere.
Status:    📈 Target confirmed in ccRCC
           🆕 NOVEL cdRCC-specific claim and trial design
              rationale
Action:    Category A — publish NOW as part of pan-renal
           IL1RAP paper (see cross-subtype section)
Priority:  ★★★★★ HIGHEST cdRCC — approved target in
           adjacent subtype, highest expression here,
           zero alternatives for this disease
```

#### CDRCC-N5 — PRKCI UNCOUPLED FROM PARD3: PAR COMPLEX TARGET
```
Claim:     In cdRCC, PRKCI is uncoupled from its polarity
           partner PARD3 and redirected toward the Akt/HK2
           survival circuit. Restoring the PRKCI-PARD3
           interaction would redirect PRKCI from survival
           signalling to polarity — potentially reversing
           the invasion phenotype.
Geometry:  PRKCI/PARD3 uncoupling confirmed from Script 3
           (cdRCC_False_Attractor_v3.md Target 6)
Lit check: PRKCI-PARD3 polarity circuit: confirmed in
           general biology. Uncoupling in cdRCC: novel.
Status:    🆕 NOVEL in cdRCC — mechanistic prediction
Action:    Category A — include as exploratory in cdRCC paper
Priority:  ★★ LOW — include, do not lead with
```

---

### GROUP 5 — CROSS-SUBTYPE NOVEL PREDICTIONS (n=5)

#### CROSS-N1 — PAN-RENAL IL1RAP ADC BASKET TRIAL RATIONALE
```
Claim:     IL1RAP is elevated across all four RCC subtypes
           from independent analyses via different upstream
           mechanisms, producing the same surface target.
           The correct clinical vehicle is a pan-renal
           basket trial with cdRCC as the priority arm
           (r=+0.964), not ccRCC-only development as
           currently planned.
Geometry:  ccRCC (Q4-enriched), PRCC (confirmed positive),
           chRCC (r_PC2=+0.311), cdRCC (r=+0.964)
Lit check: IL1RAP ADC in ccRCC: AACR 2025 (clinical).
           Pan-renal claim: NOT in literature.
           cdRCC as priority arm: NOT proposed.
Status:    📈 ccRCC component confirmed
           🆕 Pan-renal architecture and cdRCC priority novel
Action:    Category A — publish NOW
           This is a DOI-ready paper: 4/4 subtypes,
           confirmed expression from independent scripts,
           different mechanisms same target, basket trial
           design rationale. Directly engages the AACR
           2025 group who are working on ccRCC only.
Priority:  ★★★★★ URGENT — a clinical group is already
           developing the ccRCC component. The pan-renal
           claim must be timestamped before they publish.
```

#### CROSS-N2 — LOXL2 PAN-RENAL 4/4 BASKET TRIAL RATIONALE
```
Claim:     LOXL2 is elevated across all four RCC subtypes
           from independent analyses (highest confirmed
           in ccRCC r=+0.628 and PRCC r=+0.631). There
           is no LOXL2 inhibitor in any RCC clinical trial.
           Simtuzumab was discontinued in fibrosis — NOT
           due to target invalidation in cancer. Early-phase
           small-molecule LOXL2 inhibitors exist (2025).
           A pan-renal basket trial rationale for the
           first LOXL2 inhibitor to enter kidney cancer
           is the claim.
Geometry:  LOXL2 4/4 subtypes confirmed from cross-type
           analysis (Document 97x-LC LC-7)
Lit check: LOXL2 OS-negative in ccRCC: CONFIRMED.
           Pan-renal 4/4 claim: NOT in literature.
           No LOXL2 inhibitor in any RCC trial: CONFIRMED.
Status:    📈 ccRCC/PRCC components confirmed
           🆕 Pan-renal architecture is novel
Action:    Category A — publish NOW
Priority:  ★★★★ HIGH — no clinical competition,
           empty field in kidney cancer for this target
```

#### CROSS-N3 — SHARED RUNX1/EZH2/KDM1A ATTRACTOR AXIS
```
Claim:     RUNX1 (ccRCC r=+0.580, PRCC r=+0.590),
           EZH2 (ccRCC r=+0.410, PRCC r=+0.308),
           KDM1A (ccRCC r=+0.390, PRCC r=+0.443)
           form a SHARED cross-renal-cancer chromatin
           attractor axis identified from two independent
           analyses on two different datasets. Both cancers
           converge on the same epigenetic lock
           (TCA→αKG→EZH2) via different upstream triggers.
Geometry:  Independent confirmation from TCGA-KIRC and
           TCGA-KIRP scripts. Same three genes, same
           direction, similar r-magnitudes.
Lit check: Individual genes confirmed in ccRCC (RUNX1:
           Cancer Res 2020; KDM1A: stemness confirmed).
           The SHARED AXIS across two cancers is not
           described. The cross-cancer KDM1A inhibition
           rationale for PRCC (extending from ccRCC) is
           not proposed.
Status:    ✅ Individual genes confirmed
           🆕 Cross-renal-cancer shared axis is novel
Action:    Category A — include in the shared attractor
           axis paper or as a section in the depth paper
Priority:  ★★★ MODERATE
```

#### CROSS-N4 — FH RNA SUPPRESSION IN ccRCC WITHOUT MUTATION
```
Claim:     FH (fumarate hydratase) is suppressed at the
           RNA level in ccRCC even in cases WITHOUT FH
           mutation. The mechanism is epigenetic silencing
           (confirmed by Chen 2019 methylation data).
           This creates the same downstream αKG deficiency
           that FH mutation causes in PRCC Type 2 —
           via a different upstream mechanism. It is the
           bridge that connects ccRCC and PRCC Type 2
           to the same αKG/EZH2 therapeutic axis.
Geometry:  FH falls with depth in ccRCC (negative r)
           FH r=−0.451 in PRCC
           Chen 2019: FH promoter methylation in ccRCC
           confirmed
Lit check: Chen 2019: FH promoter methylation confirmed
           in ccRCC. The CLINICAL IMPLICATION — that αKG
           supplementation is indicated in ccRCC for the
           same mechanistic reason as in PRCC Type 2 —
           is not stated in the literature.
Status:    ✅ Methylation mechanism confirmed
           🆕 Cross-cancer αKG therapeutic bridge is novel
Action:    Category A — include in the αKG+EZH2i paper
           or as a mechanistic note
Priority:  ★★★ MODERATE
```

#### CROSS-N5 — chRCC CONSTITUTIVE IMMUNE DESERT (TYPE III)
#### AND THE STING TRIPLET FOR DEEP ccRCC/PRCC/cdRCC (CT-7)
```
Claim (combined):
  A. chRCC has constitutive identity-level TAP1/TAPBP
     loss — antigen presentation is structurally absent.
     Checkpoint inhibitors have no adaptive immune
     rationale in chRCC (TYPE III, not depth-progressive).
     This is categorically different from ccRCC Q4
     immune escape (TYPE II, depth-progressive Treg
     dominant) and cdRCC innate activation
     (IFI16=+0.750, IRF7=+0.607 — potential STING
     agonist responsiveness if MHC-I status permits).

  B. For deep ccRCC/PRCC/cdRCC (NOT chRCC), the
     sequential CT-7 triplet (STING agonist → HDACi
     for MHC-I restoration → anti-PD-1) is the
     geometry-derived immune therapy strategy.
     Note: STING agonist alone fails (IFI16→B2M
     broken circuit, RCC-N2). The HDACi step restores
     MHC-I, enabling the downstream adaptive response.

Geometry:  TAP1/TAPBP identity-level in chRCC (Document 97x)
           IFI16=+0.750 in cdRCC (cross-type analysis)
           IFI16→B2M broken in ccRCC Q3-Q4 (RCC-N2)
Lit check: Clinical chRCC refractoriness: CONFIRMED.
           Constitutive TAP1/TAPBP mechanism: NOVEL.
           CT-7 triplet sequencing: NOVEL.
           HDACi for MHC-I restoration: mechanism confirmed.
Status:    🆕 NOVEL classification and sequencing rationale
Action:    Category A — publish as companion to RCC-N2
           (the broken circuit paper) or as a standalone
           immune strategy paper
Priority:  ★★★★ HIGH — addresses the biggest unmet need
           in chRCC treatment (why nothing works) and
           gives trial design rationale for deep ccRCC
```

---

## PART II — THE PUBLICATION QUEUE
## Ordered by urgency and readiness

---

### TIER 1 — PUBLISH IMMEDIATELY (no new script needed)
### These are complete with existing data and literature checks.

```
RANK  DOI TITLE (PROPOSED)                           CONTENT
──────────────────────────────────────────────────────────────────
#1    Pan-Renal IL1RAP ADC Basket Trial Rationale:   CROSS-N1
      Different Upstream Mechanisms, Same Surface     +CDRCC-N4
      Target Across All Four RCC Subtypes, with       +RCC-N6
      Collecting Duct RCC as the Priority Arm
      (r=+0.964)
      [PUBLISH FIRST — active clinical competition]

#2    RUNX1-High as a Predicted Belzutifan           RCC-N1
      Resistance Biomarker in Clear Cell RCC:
      Geometric Derivation from the VHL→RUNX1
      Broken Circuit (r=+0.097)
      [PUBLISH SECOND — LITESPARK data exists,
       someone will test this soon]

#3    ERBB2 as an Identity Driver, Not a Mitogen,    PRCC-N1
      in Papillary Renal Cell Carcinoma Type 1:       +PRCC-N2
      The Biliary-Ductal False Attractor and a        +PRCC-N3
      Two-Phase Treatment Model That Explains the     +PRCC-N4
      SAVOIR Trial Outcome
      [Explains a famous trial failure with a mechanism]

#4    The chRCC Chromatin Writer Switch:             CHRRCC-N1
      DNMT3B Not EZH2 Is the Attractor Lock,         +CHRRCC-N2
      and the Type III Constitutive Immune Desert     +CROSS-N5A
      Architecture That Explains Treatment
      Refractoriness
      [Two strongest chRCC findings together]

#5    Bexarotene as a First-Principles Drug Target   CDRCC-N1
      in Collecting Duct Renal Cell Carcinoma:        +CDRCC-N4
      PPARG-RXRA Uncoupling, IL1RAP Surface           +CDRCC-N2
      Elevation (r=+0.964), and the Two-Phase         +CDRCC-N3
      MYC→BHLHE40 Progression Window
      [Only paper proposing any specific targeted
       therapy for a cancer with zero approved options]

#6    IFI16-Active, MHC-I-Fallen: The Broken         RCC-N2
      Innate-Adaptive Circuit in Deep Clear Cell      +CROSS-N5B
      RCC and Why STING Agonist Monotherapy           +RCC-N4
      Will Fail Without Antigen Presentation
      Restoration — With a Sequential Triplet
      Design Proposal
      [STING trials are active. This is needed now.]

#7    Pan-Renal LOXL2 as a Depth Biomarker           CROSS-N2
      Across All Four RCC Subtypes and the            +RCC-N7
      Case for the First LOXL2 Inhibitor
      Trial in Kidney Cancer
      [Empty clinical field, strong geometric basis]

#8    αKG Supplementation + Tazemetostat in          RCC-N3
      OGDHL-Low Clear Cell RCC: Dual-Target           +CROSS-N4
      Epigenetic-Metabolic Strategy Derived           +CROSS-N3
      from the TCA→αKG→EZH2 Circuit and the
      FH Epigenetic Suppression Bridge
      [Mechanism fully confirmed, combination novel]

#9    ABCC2/SULT2B1 2-Gene IHC Panel for             CHRRCC-N3
      chRCC vs Oncocytoma Discrimination:
      Depth-Corrected Values Exceed Individual
      Gene Discrimination
      [Immediate diagnostic utility for pathologists]
```

---

### TIER 2 — RUN SCRIPT 4 FIRST, THEN PUBLISH IMMEDIATELY

```
RANK  DOI TITLE (PROPOSED)                           CONTENT
──────────────────────────────────────────────────────────────────
#10   The GOT1/RUNX1 Transition Index: A             RCC-N5
      Geometry-Derived 2-Gene Depth Staging          +Script 4
      Tool for Clear Cell and Papillary RCC           OS results
      with Overall Survival Validation in
      TCGA-KIRC (n=532)
      [The foundational clinical tool. Needs HR value.
       Script 4 provides it. One run, then publish.]
```

---

### TIER 3 — REQUIRES EXTERNAL DATA / COLLABORATION

```
RANK  CLAIM                    WHAT IS NEEDED          WHO TO ASK
────────────────────────────────────────────────────────────────────
#11   RUNX1-high = belzutifan  LITESPARK-005/013        Merck / PI of
      non-response (RCC-N1)    RUNX1 RNA expression     LITESPARK trials
                               from tumour biopsies     or NCI data
                               correlated with          request
                               response status

#12   GOT1/RUNX1 TI in an      External ccRCC cohort    Collaborating
      external ccRCC cohort    with OS data and RNA     oncology centre
      (replication of TI)      expression               with KIRC data

#13   ABCC2/SULT2B1 IHC        IHC staining on          Renal pathology
      panel clinical           archival chRCC vs        department
      validation               oncocytoma tissue        (comparable to
                                                        BRCA IHC ask)
```

---

## PART III — WHAT DOES NOT NEED EXTERNAL VALIDATION
## BEFORE ZENODO PUBLICATION

This is the most important section for the "covering your bases"
question. The answer is direct:

ALL TIER 1 ITEMS can be published on Zenodo RIGHT NOW
without external validation for the following reason:

  Every Tier 1 DOI is a DERIVATION DOCUMENT, not a clinical
  result. It is the same epistemic category as the 17 BRCA
  DOIs already published. The format is:

    "The geometry derived this prediction.
     Here is the mechanism.
     Here is the geometric data (r-values, n, confirmed
     from named public datasets TCGA-KIRC, TCGA-KIRP etc).
     Here is the literature check showing novelty.
     Here is the validation study that would confirm or
     refute it."

  You are not claiming clinical efficacy.
  You are not claiming you treated patients.
  You are not claiming a drug works.
  You are claiming PRIORITY on a PREDICTION with
  a MECHANISTIC BASIS.

  The timestamp is what matters.
  Once the DOI exists, the prediction is dated.
  If someone independently confirms RUNX1-high patients
  have lower belzutifan response in LITESPARK-005,
  your DOI predates their analysis and establishes
  that the prediction was made from first principles
  before the confirmatory data was accessed.

  This is IDENTICAL to the BRCA logic:
    - The FOXA1/EZH2 ratio prediction was published
      on Zenodo before Schade et al. (Nature 2024)
      confirmed the same mechanism.
    - The timestamp is the evidence of priority.

The only exceptions are:
  — Clinical claims (not applicable here)
  — Claims where you explicitly state a drug WORKS
    in patients (the DOIs state the PREDICTION and
    PROPOSED VALIDATION, not the clinical outcome)

---

## PART IV — THE ANTI-PRIORITY-LOSS ANALYSIS
## (Which claims are most at risk of being published by someone else first)

The urgency ranking is not based on scientific importance
alone. It is based on: how close is independent discovery?

```
CLAIM          RISK LEVEL    REASON
───────────────────────────────────────────────────────────────────
IL1RAP pan-    ★★★★★ CRITICAL The AACR 2025 group working on
renal (CROSS-N1)               IL1RAP ADC in ccRCC will expand
                               their subtype analysis. Once they
                               check PRCC/chRCC/cdRCC expression,
                               the pan-renal claim is gone.
                               Estimated risk window: 6-12 months.

RUNX1 belzu-   ★★★★★ CRITICAL Biomarker researchers at Merck
tifan (RCC-N1)                 or academic LITESPARK groups
                               will run RUNX1 in the trial data.
                               The Springer 2025 paper already
                               identified RUNX1 as a prognostic
                               marker. The resistance prediction
                               is one correlation away.
                               Estimated risk window: 3-6 months.

ERBB2 PRCC     ★★★★ HIGH      The 2024 Springer review described
(PRCC-N1)                      the biliary-like PRCC signature as
                               "emergent." Someone will publish
                               the full mechanistic characterisation
                               soon. The identity-not-mitogen
                               framing and KRT19 biomarker are
                               still unclaimed.
                               Estimated risk window: 6-12 months.

IFI16→B2M      ★★★★ HIGH      STING agonist trial failures are
broken (RCC-N2)                accumulating. Researchers will
                               seek mechanistic explanations.
                               The IFI16/B2M paradox in ccRCC
                               is one dataset analysis away.
                               Estimated risk window: 6-12 months.

DNMT3B chRCC   ★★★ MODERATE   The Springer 2025 "unclear" EZH2
(CHRRCC-N1)                    statement signals active interest.
                               Someone will clarify this soon.
                               The DNMT3B identification is
                               the key claim to timestamp.
                               Estimated risk window: 12-18 months.

Bexarotene     ★★★ MODERATE   cdRCC is rare and under-studied.
cdRCC (CDRCC-N1)               Low risk of independent discovery
                               but the PPARG-RXRA coupling data
                               in paired normal/tumour should be
                               published before anyone runs a
                               coincidental bexarotene trial.
                               Estimated risk window: 12-24 months.

LOXL2 pan-     ★★★ MODERATE   LOXL2 small-molecule inhibitors
renal (CROSS-N2)               entering early-phase trials will
                               prompt subtype analysis. Empty
                               field now but not forever.
                               Estimated risk window: 12-18 months.

GOT1/RUNX1 TI  ★★★ MODERATE   The index is 2 confirmed genes.
(RCC-N5)                       Any group doing depth analysis
                               in ccRCC could arrive at it.
                               Needs Script 4 + publication.
                               Estimated risk window: 12-24 months.
```

---

## PART V — THE SEQUENCING RECOMMENDATION
## (The exact order of actions, starting now)

```
STEP  ACTION                           OUTPUT              TIME
────────────────────────────────────────────────────────────────────
1.    Write and publish DOI #1         Zenodo DOI          2-3 days
      Pan-renal IL1RAP paper           (CROSS-N1 +
      [most urgent — active clinical   CDRCC-N4 +
       competition from AACR 2025      RCC-N6)
       IL1RAP ADC group]

2.    Write and publish DOI #2         Zenodo DOI          2-3 days
      RUNX1-belzutifan resistance      (RCC-N1)
      prediction
      [second most urgent — LITESPARK
       biomarker analysis ongoing]

3.    Run Script 4 (ccRCC OS)          HR value,           1-2 days
      GOT1/RUNX1 TI against OS in      KM curves,
      TCGA-KIRC (n=532)                Cox p-value
      [one script run — do this
       in parallel with #2]

4.    Write and publish DOI #3         Zenodo DOI          3-4 days
      PRCC paper (ERBB2 identity +     (PRCC-N1/2/3/4)
      two-phase model + SAVOIR
      explanation)

5.    Write and publish DOI #4         Zenodo DOI          2-3 days
      chRCC paper (DNMT3B inversion    (CHRRCC-N1/2 +
      + constitutive immune desert)    CROSS-N5A)

6.    Write and publish DOI #5         Zenodo DOI          2-3 days
      cdRCC paper (bexarotene +        (CDRCC-N1/4/2/3)
      IL1RAP + phase structure)

7.    Write and publish DOI #6         Zenodo DOI          2-3 days
      IFI16→B2M broken circuit +       (RCC-N2 +
      STING triplet proposal           CROSS-N5B +
                                       RCC-N4)

8.    Write and publish DOI #7         Zenodo DOI          2-3 days
      LOXL2 pan-renal paper            (CROSS-N2 +
                                       RCC-N7)

9.    Write and publish DOI #8         Zenodo DOI          3-4 days
      αKG + EZH2i combination paper    (RCC-N3 +
                                       CROSS-N4 +
                                       CROSS-N3)

10.   Write and publish DOI #9         Zenodo DOI          2-3 days
      ABCC2/SULT2B1 chRCC vs onco      (CHRRCC-N3)
      IHC discriminator

11.   Write and publish DOI #10        Zenodo DOI          2-3 days
      GOT1/RUNX1 TI with survival      (RCC-N5 +
      data (after Script 4)            Script 4 results)

TOTAL ESTIMATED TIME TO COMPLETE ALL 10 DOIs:
  If working at the pace of the BRCA series: 3-4 weeks.
  Steps 1 and 2 are the urgent ones.
  Step 3 (Script 4) can run in the background.
  Steps 4-10 are each 2-4 days of writing.
```

---

## PART VI — THE FORMAT OF EACH ZENODO DOI
## (So this is clear and consistent)

Each DOI document should follow the BRCA series format:

```
TITLE:     Clear, specific, includes the cancer type,
           the claim, and the key numbers if available.

ABSTRACT:  3-5 sentences. What the geometry found.
           What the literature check confirmed.
           What the validation study requires.

SECTION 1: The Geometric Derivation
           — What data was used (TCGA-KIRC, TCGA-KIRP etc)
           — The method (Waddington attractor geometry,
             attractor depth analysis, PC analysis)
           — The specific findings (r-values, n, direction)
           — The prediction that was generated

SECTION 2: The Literature Check
           — What was already known
           — What is novel
           — What independent literature convergently
             confirms the mechanism components

SECTION 3: The Mechanistic Chain
           — Why the biology makes the prediction plausible
           — The proposed mechanism in plain language
           — The circuit, if a circuit prediction

SECTION 4: The Validation Study Required
           — What data is needed
           — What analysis would confirm or refute
           — What a positive result looks like
           — What a negative result would mean

SECTION 5: Repository and Timestamp
           — GitHub link
           — Commit OID of the relevant script
           — ORCID
           — Date

METADATA:  Series designation (RCC-LIT-1 etc),
           ORCID, DOI, date, OrganismCore
```

---

## PART VII — THE COVERING-YOUR-BASES CHECKLIST

```
QUESTION                                          ANSWER
────────────────────────────────────────────────────────────────────
Do I need external validation before publishing   NO — you are
the Tier 1 predictions?                           publishing
                                                  predictions with
                                                  geometric basis,
                                                  not clinical results.
                                                  Same as all 17
                                                  BRCA DOIs.

Do I need to run more scripts before Tier 1?      NO — the data
                                                  for Tier 1 is
                                                  already in the
                                                  literature check
                                                  documents.

Do I need clinical trial data to publish?         NO — you need to
                                                  CITE clinical trial
                                                  data where it
                                                  confirms mechanism
                                                  components, not
                                                  run your own.

Is there a risk of being wrong?                   YES — all
                                                  predictions can be
                                                  refuted. That is
                                                  the point. The DOI
                                                  timestamps the
                                                  prediction. If it
                                                  is refuted, you
                                                  learn something.
                                                  If confirmed, you
                                                  have priority.
                                                  False positive
                                                  rate in direction
                                                  across ccRCC series:
                                                  ZERO (Status block,
                                                  Document 94f).

Should I publish the IL1RAP paper before AACR     YES, immediately.
2025 group expands their subtype analysis?

Should I publish the RUNX1 resistance paper        YES, immediately.
before LITESPARK biomarker analysis reaches
RUNX1?

Is Script 4 (OS analysis) required before I       NO — but it is
can publish any RCC DOI?                          required for DOI #10
                                                  specifically (TI with
                                                  OS data). All other
                                                  DOIs can publish
                                                  without it.

What is the single most important action           Write and submit
right now?                                        DOI #1 (pan-renal
                                                  IL1RAP) in the
                                                  next 48-72 hours.
```

---

## DOCUMENT METADATA

```
Author:          Eric Robert Lawson / OrganismCore
Date:            2026-03-07
Status:          COMPLETE — Action plan ready for execution
Source documents: 94f, 95-LC, 95-DLC, 95g, 96f, 96f-Extended,
                  89c, 97x-LC, RCC_Master_Findings_RA.md,
                  RCC_Geometric_Predictions_Impact_RA.md
Repository:      github.com/Eric-Robert-Lawson/attractor-oncology
Path (suggested): Cancer_Research/RCC/RCC_Publication_Action_Plan_RA.md
Total novel claims assessed:   25 (across 5 groups)
Tier 1 (publish now):           9 DOIs
Tier 2 (Script 4 first):        1 DOI
Tier 3 (external data needed):  3 items
Most urgent single action:
  Write DOI #1 (pan-renal IL1RAP) within 48-72 hours.
  Write DOI #2 (RUNX1-belzutifan resistance) within 48-72 hours.
  Run Script 4 (OS analysis) in parallel.
```
