#### OV — Ovarian Cancer
```
Lineage:  Fallopian tube epithelial
          (STIC origin) or
          ovarian surface epithelial
Block:    HGSOC cells vs normal
          fallopian tube epithelium
Predicted switch genes:
  PAX8   — Müllerian epithelial TF
  WT1    — ovarian surface marker
  OVGP1  — oviductal glycoprotein
           (fallopian terminal)
Data:    GSE154600 or
         TCGA + scRNA-seq atlas
Doc: 90
```

#### UCEC — Endometrial Cancer
```
Lineage:  Endometrial epithelial
Block:    Endometrial cancer vs
          normal endometrium
Predicted switch genes:
  FOXA2  — confirmed LUAD
           ALSO expressed in
           endometrium —
           second test of FOXA2
           in a different tissue
  PGR    — progesterone receptor
           (endometrial terminal)
  HAND2  — endometrial stromal TF
Data:    GSE213216 or TCGA UCEC
         + scRNA atlas
Note:    FOXA2 overlap with LUAD
         would be second confirmation
         of FOXA2 as a switch gene
         in a different tissue context.
Doc: 91
```

---

### Session 4 — Rare and Aggressive

#### PDAC Subtypes
```
Classical vs Basal-like PDAC
Two false attractors within one cancer.
Same approach as BRCA LumA vs TNBC.
Doc: 92
```

#### Uveal Melanoma
```
Lineage:  Melanocyte/neural crest
Block:    Uveal melanoma vs normal
          uveal melanocytes
Predicted switch genes:
  MITF   — melanocyte master TF
  DCT    — dopachrome tautomerase
           (terminal melanin synthesis)
  TYRP1  — tyrosinase-related protein
Data:    GSE139829
Note:    MITF is the melanocyte
         equivalent of FOXA1 in breast.
         One master TF defining
         terminal melanocyte identity.
Doc: 93
```

#### Cutaneous Melanoma
```
Same logic as uveal but different
microenvironment and mutation spectrum.
MITF, DCT, TYRP1 as switch genes.
BRAF V600E drives the false attractor.
Data: GSE215120 or GSE72056
Doc: 94
```

#### Mesothelioma
```
Lineage:  Mesothelial cell
Block:    Mesothelioma vs normal
          mesothelium
Predicted switch genes:
  WT1    — mesothelial identity TF
  MSLN   — mesothelin (terminal)
  CALB2  — calretinin (terminal marker)
Data:    GSE195615
Doc: 95
```

#### Thyroid Cancer — PTC/FTC
```
Lineage:  Thyroid follicular cell
Block:    Thyroid cancer vs normal
          follicular cells
Predicted switch genes:
  NKX2-1 — confirmed LUAD partial
           ALSO the thyroid identity TF
           (TTF-1 is NKX2-1 —
           used in both lung and thyroid)
           Second test of NKX2-1 in
           its primary thyroid context
  PAX8   — thyroid specification TF
  TG     — thyroglobulin (terminal)
  TSHR   — TSH receptor (terminal)
Data:    GSE184362 or GSE213647
Note:    NKX2-1 was partial in LUAD
         (scaffold). In thyroid it may
         be a switch gene — its primary
         tissue. This would refine the
         scaffold/switch distinction:
         NKX2-1 is scaffold in lung
         but switch in thyroid.
Doc: 96
```

#### Neuroblastoma
```
Lineage:  Sympathetic neuron /
          chromaffin cell
Block:    Neuroblastoma cells vs
          normal sympathetic neurons
Predicted switch genes:
  PHOX2B — sympathetic neuron TF
  DBH    — dopamine beta-hydroxylase
           (terminal chromaffin)
  TH     — tyrosine hydroxylase
Data:    GSE137804 or GSE137804
Doc: 97
```

#### Medulloblastoma
```
Lineage:  Cerebellar granule neuron
Block:    MB cells vs normal
          granule neuron precursors
Predicted switch genes:
  ATOH1  — granule neuron TF
           (was in CRC panel, absent)
           Now test in correct tissue
  NEUROD1 — neuronal differentiation
  RBFOX3  — mature neuron marker
Data:    GSE119926 or GSE155446
Note:    ATOH1 was predicted in CRC
         but missing from panel.
         First test in its actual
         tissue of function.
Doc: 98
```

---

### Session 5 — Liquid Tumors and
### Rare Hematopoietic

#### DLBCL — Diffuse Large B Cell Lymphoma
```
Lineage:  Germinal center B cell
          → plasma cell
Block:    DLBCL cells vs normal GC
          B cells or plasma cells
Predicted switch genes:
  PRDM1  — plasma cell TF
  IRF4   — plasma cell identity
  BLIMP1 — terminal B cell
Data:    GSE181063 or GSE132509
Doc: 99
```

#### Follicular Lymphoma
```
Lineage:  Follicular B cell
Block:    FL cells vs normal
          follicular B cells
Predicted switch genes:
  BCL6   — germinal center master TF
           (loss of BCL6 is required
           for terminal B cell
           differentiation — this
           may INVERT the prediction:
           BCL6 ELEVATED in FL because
           it PREVENTS terminal
           differentiation)
  PRDM1  — downstream of BCL6
Data:    GSE181063
Note:    BCL6 is a repressor of
         terminal B cell
         differentiation. FL may be
         the first case where the
         false attractor is maintained
         by ACTIVE expression of a TF
         that represses the switch genes,
         rather than by suppression of
         the switch genes directly.
         Framework stress test.
Doc: 100
```

#### T Cell Lymphoma
```
Lineage:  T cell (various subtypes)
Block:    PTCL cells vs normal
          mature T cells
Predicted switch genes:
  TBX21  — Th1 identity TF
  GATA3  — Th2 identity TF
           (confirmed BRCA, ALL —
           third cancer test)
  RORC   — Th17 terminal TF
Data:    GSE188053
Doc: 101
```

#### Mast Cell Disease / MCL
```
Lineage:  Mast cell
Block:    Mastocytosis vs normal
          mast cells
Predicted switch genes:
  MITF   — mast cell identity TF
           (also melanocyte — second
           tissue test)
  TPSAB1 — tryptase (terminal)
Data:    GSE141560
Doc: 102
```

---

### Session 6 — Brain Tumors

#### LGG — Low Grade Glioma (IDH-mutant)
```
Lineage:  Oligodendrocyte or astrocyte
Block:    IDH-mutant glioma vs
          normal oligodendrocyte
          (same switch genes as GBM?)
Predicted switch genes:
  SOX10  — confirmed GBM 88.6%
  MBP    — confirmed GBM 89.6%
  PLP1   — confirmed GBM 83.4%
Data:    GSE131928 (Neftel — same
         dataset, IDH-mutant subset)
Note:    Direct comparison with GBM.
         IDH-mutant glioma is less
         aggressive than IDH-wt GBM.
         Does it show less switch gene
         suppression? This would be
         the first test of whether
         suppression magnitude
         correlates with clinical
         aggressiveness.
Doc: 103
```

#### Ependymoma
```
Lineage:  Ependymal cell
Block:    Ependymoma vs normal
          ependymal cells
Predicted switch genes:
  FOXJ1  — ependymal/ciliated TF
  CFAP126 — ciliogenesis (terminal)
Data:    GSE141383
Doc: 104
```

#### Oligodendroglioma
```
Lineage:  Oligodendrocyte (IDH+1p19q)
Block:    Oligodendroglioma cells vs
          normal oligodendrocytes
Same switch genes as GBM predicted.
Third test of SOX10/MBP/PLP1 axis.
Data:    GSE131928 subset
Doc: 105
```

---

### Session 7 — Cross-Tissue Validation

#### Matched Primary and Metastasis
```
Use LUAD dataset (GSE131907) which
contains primary tumor, brain
metastasis, and pleural effusion
from the same patients.

Question: Do metastatic cells show
MORE or LESS switch gene suppression
than primary tumor cells?

Prediction: Metastatic cells are
MORE suppressed — deeper in the
false attractor. The metastatic
state requires losing even more
differentiation identity.

This tests whether the false attractor
depth correlates with metastatic
potential.
Doc: 106
```

#### Chemotherapy-Resistant Subpopulations
```
Multiple datasets contain matched
pre- and post-treatment samples.

Question: Do resistant cells show
MORE switch gene suppression than
sensitive cells?

Prediction: Yes. Drug resistance
selects for deeper false attractor
states — cells that are more
dedifferentiated and therefore
more insensitive to
differentiation-based signals.

This would provide a mechanism
for acquired resistance:
resistance = deeper false attractor.

Data: GSE161533 (BRCA chemo)
      GSE150949 (LUAD EGFR inhibitor)
Doc: 107
```

#### Pediatric vs Adult Same Cancer
```
Compare pediatric GBM vs adult GBM
(already in Neftel dataset —
pediatric samples included).

Question: Same switch gene suppression
in pediatric GBM?

This tests whether the false attractor
is age-invariant.
Doc: 108
```

---

### Session 8 — Synthetic Lethal Tests

#### AML with DNMT3A mutation
```
Subset the AML dataset by mutation
status (DNMT3A, FLT3, NPM1).
Do different mutation backgrounds
show different depths of switch gene
suppression?

This connects mutation → epigenetic
state → false attractor depth.
Doc: 109
```

#### BRCA1/2 mutant vs sporadic BRCA
```
GSE176078 has BRCA1/2 mutation status.
Compare FOXA1/GATA3/ESR1 suppression
in BRCA1-mutant TNBC vs sporadic TNBC.

Prediction: BRCA1-mutant shows deeper
suppression — the germline mutation
predisposes to a deeper false attractor.
Doc: 110
```
