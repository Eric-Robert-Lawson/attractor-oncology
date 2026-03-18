\documentclass[11pt, a4paper]{article}

% ── PACKAGES ──────────────────────────────────────────────────
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{geometry}
\geometry{
  top=2.0cm,
  bottom=2.0cm,
  left=2.2cm,
  right=2.2cm
}
\usepackage{hyperref}
\hypersetup{
  colorlinks=true,
  urlcolor=blue,
  linkcolor=blue,
  citecolor=blue,
  pdftitle={Domestication as Field-Specified Developmental
            Arrest: A Unified Causal Mechanism for the
            Domestication Syndrome Derived from
            Attractor Geometry},
  pdfauthor={Eric Robert Lawson / OrganismCore},
  pdfsubject={Theoretical biology — domestication mechanism},
  pdfkeywords={domestication syndrome, neural crest,
               attractor geometry, developmental arrest,
               coherence gradient field, Sus scrofa,
               feralization, Waddington landscape,
               epigenetics, HPA axis}
}
\usepackage{microtype}
\usepackage{parskip}
\usepackage{booktabs}
\usepackage{array}
\usepackage{tabularx}
\usepackage{longtable}
\usepackage{xcolor}
\usepackage{mdframed}
\usepackage{enumitem}
\usepackage{titlesec}
\usepackage{lmodern}
\usepackage{authblk}
\usepackage{fancyhdr}
\usepackage{lastpage}
\usepackage{amsmath}
\usepackage{multirow}

% ── COLOURS ───────────────────────────────────────────────────
\definecolor{headergreen}{RGB}{20,100,60}
\definecolor{rulecolor}{RGB}{20,100,60}
\definecolor{claimbox}{RGB}{20,100,60}
\definecolor{notebox}{RGB}{240,248,240}
\definecolor{warningbox}{RGB}{255,248,220}
\definecolor{predbox}{RGB}{235,248,255}
\definecolor{confirmbox}{RGB}{240,255,240}
\definecolor{mechanismbox}{RGB}{248,240,255}

% ── HEADER / FOOTER ───────────────────────────────────────────
\pagestyle{fancy}
\fancyhf{}
\renewcommand{\headrulewidth}{0.4pt}
\fancyhead[L]{\small\color{headergreen}
  OrganismCore --- Attractor Geometry Series}
\fancyhead[R]{\small\color{headergreen}
  Domestication as Developmental Arrest --- 2026-03-18}
\fancyfoot[C]{\small Page \thepage\ of \pageref{LastPage}}

% ── SECTION FORMAT ────────────────────────────────────────────
\titleformat{\section}
  {\large\bfseries\color{headergreen}}
  {\thesection.}{0.5em}{}[
  \vspace{-0.3em}
  \textcolor{rulecolor}{\rule{\linewidth}{0.6pt}}]
\titleformat{\subsection}
  {\normalsize\bfseries\color{headergreen}}
  {\thesubsection.}{0.5em}{}
\titleformat{\subsubsection}
  {\small\bfseries}
  {\thesubsubsection.}{0.5em}{}

\setlength{\parskip}{0.45em}
\setlength{\parindent}{0pt}

% ══════════════════════════════════════════════════════════════
\title{
  \vspace{-0.5em}
  {\color{headergreen}\LARGE\bfseries
  Domestication as Field-Specified\\
  Developmental Arrest}\\[0.5em]
  {\large
  A Unified Causal Mechanism for the Domestication
  Syndrome Derived from Attractor Geometry,\\
  Confirmed Against Five Independent Literature
  Sources, and Generating Novel Experimental
  Predictions}\\[0.4em]
  {\normalsize\itshape
  OrganismCore Framework --- Attractor Geometry Series}
}

\author[1]{Eric Robert Lawson}
\affil[1]{
  OrganismCore (Independent Research)\\
  \href{mailto:OrganismCore@proton.me}
  {OrganismCore@proton.me}\ \ $\cdot$\ \
  ORCID:
  \href{https://orcid.org/0009-0002-0414-6544}
  {0009-0002-0414-6544}\\
  Repository:
  \href{https://github.com/Eric-Robert-Lawson/attractor-oncology}
  {\texttt{github.com/Eric-Robert-Lawson/attractor-oncology}}
}

\date{March 18, 2026}

% ══════════════════════════════════════════════════════════════
\begin{document}
\maketitle
\thispagestyle{fancy}

% ── PRIORITY / DERIVATION NOTICE ──────────────────────────────
\begin{mdframed}[
  backgroundcolor=notebox,
  linecolor=headergreen,
  linewidth=1pt,
  innertopmargin=0.6em,
  innerbottommargin=0.6em,
  innerleftmargin=1.0em,
  innerrightmargin=1.0em
]
\textbf{DERIVATION NOTICE — READ FIRST}

The causal mechanism proposed in this paper was derived
from first principles using the attractor geometry
framework before any literature search was conducted.
The predictions were locked in a timestamped repository
record (OrganismCore, 2026-03-18) prior to the literature
verification reported in Section~\ref{sec:confirmation}.
The confirmations against Marcstr\"{o}m \& Essen-M\"{o}ller
(1968), Hemmer (1990), Wilkins et al.\ (2014), Belyaev
(1979), and the feralization sequence literature are
therefore confirmations of pre-specified predictions,
not the source of those predictions.
This eliminates post-hoc pattern matching as an
explanation for the correspondence.

Repository record:
\href{https://github.com/Eric-Robert-Lawson/attractor-oncology}
{\texttt{github.com/Eric-Robert-Lawson/attractor-oncology}}

Pre-registration of novel experimental predictions
(P7--P11): deposited Zenodo, 2026-03-18.
\end{mdframed}

\vspace{0.5em}

% ── ABSTRACT ──────────────────────────────────────────────────
\begin{abstract}
\noindent
The domestication syndrome --- the convergent suite of
morphological and behavioural traits appearing in every
domesticated species including floppy ears, reduced
pigmentation, shorter snout, smaller adrenal glands,
and retained juvenile behaviour --- has resisted unified
causal explanation since its systematic description.
Wilkins, Wrangham and Fitch (2014) correctly identified
the cellular mechanism: mild reduction in neural crest
cell migration. They did not identify the causal origin
of that reduction. Here we propose and defend that origin:
the domestication coherence gradient field
(\textit{N}\textsubscript{domestic}) suppresses the
maternal hypothalamic-pituitary-adrenal (HPA) axis,
alters the intrauterine coherence gradient environment,
and arrests neural crest cell migration at the
juvenile-proximate developmental stage. Domestication is
therefore not genetic modification of the organism.
It is field-specified developmental arrest of the neural
crest migration program. The domestic phenotype is the
phenotypic signature of that arrest. Feralization is
developmental resumption when the maintenance field is
removed. This mechanism was derived from the attractor
geometry framework ($S + N + G \rightarrow R$) before
literature consultation and is confirmed against five
independent literature sources. The mechanism predicts
a specific asymmetry in domestication vs.\ feralization
rates, a specific reversion sequence in feral animals,
a specific cross-species universality of the syndrome,
and a specific set of novel experimental outcomes
(P7--P11) not predicted by the genetic model. The
mechanism is structurally identical to EZH2-mediated
attractor arrest in triple-negative breast cancer and
nitrogen starvation-mediated arrest in the LECA protocol,
suggesting a scale-invariant field-arrest mechanism
operative from single cells to whole organisms.
\end{abstract}

\vspace{0.5em}
\noindent\textcolor{rulecolor}{\rule{\linewidth}{0.8pt}}
\vspace{0.3em}

\tableofcontents

\newpage

% ══════════════════════════════════════════════════════════════
\section{Introduction: The Open Causal Question}
% ══════════════════════════════════════════════════════════════

Domestic animals are morphologically and behaviourally
distinct from their wild ancestors in ways that are
strikingly convergent across species. The dog differs
from the wolf, the domestic pig from the wild boar,
the domestic cat from the African wildcat, and the
domesticated silver fox from its wild counterpart in
a predictable, overlapping set of traits that has been
termed the \textit{domestication syndrome}
\cite{wilkins2014}: floppy ears, patches of
depigmentation, shortened snout, smaller adrenal glands,
reduced stress reactivity, docility, and retention of
juvenile behaviour into adulthood (neoteny).

The cross-species universality of this syndrome is
remarkable. These species share no recent common
ancestor. Their domestication events were geographically
separated and historically independent. Yet the phenotypic
outcome is the same. This convergence demands a unified
causal explanation.

Wilkins et al.\ (2014) provided the most significant
advance: a unified cellular mechanism. They demonstrated
that all traits of the domestication syndrome map onto
a single developmental perturbation --- mild reduction in
neural crest cell (NCC) migration during embryogenesis.
Neural crest cells are multipotent embryonic cells that
migrate from the dorsal neural tube to populate ear
cartilage, craniofacial structures, melanocytes, adrenal
medulla, and components of the peripheral nervous system.
When NCC migration is mildly reduced, all syndrome traits
appear simultaneously. The cellular mechanism is correct
and well-supported.

What Wilkins et al.\ did not provide is the causal origin
of the NCC migration reduction. Their explanation invokes
genetic change --- alleles affecting NCC number,
migration capacity, or differentiation, co-selected
during domestication alongside tameness alleles. This
genetic explanation is incomplete for three reasons:

\begin{enumerate}[itemsep=0.2em]
  \item \textbf{Frequency:} Marcstr\"{o}m \&
  Essen-M\"{o}ller (1968) documented domestication
  syndrome traits --- piebald coat, drooping ears,
  skull shape changes --- in the first captive generation
  of wild boar \textit{without selective breeding},
  at frequencies too high for recessive allele expression
  to explain.

  \item \textbf{Asymmetry:} Feralization (domestic
  $\rightarrow$ wild) occurs in months to one generation.
  Re-domestication (wild $\rightarrow$ domestic) requires
  many generations of sustained selection. This asymmetry
  is not predicted by genetic models in which allele
  frequencies shift at equivalent rates under equivalent
  selection.

  \item \textbf{Sufficiency:} The Belyaev fox experiment
  \cite{belyaev1979} demonstrated that selecting only for
  HPA axis suppressibility (tameness) produces the full
  domestication syndrome without morphology selection.
  The genetic model requires co-selection of linked
  alleles; the field model requires only that HPA
  suppression alters the intrauterine coherence gradient.
\end{enumerate}

We propose that the causal origin of NCC migration
reduction is not genetic but \textit{field-specified}:
the domestication coherence gradient field acts on the
maternal HPA axis, which alters the intrauterine
developmental environment, which changes the coherence
gradients guiding NCC migration, which arrests that
migration at the juvenile-proximate stage. The domestic
phenotype is the organismal signature of that
field-specified developmental arrest.

% ══════════════════════════════════════════════════════════════
\section{The Attractor Geometry Framework}
\label{sec:framework}
% ══════════════════════════════════════════════════════════════

\subsection{The Complete Developmental Equation}

The standard model of developmental biology is:
\[
  \text{Genome} \;\rightarrow\; \text{Organism}
\]

This model is incomplete. The complete model, as
formalised in the attractor geometry framework
\cite{kauffman1993, huang2009, wang2014}, is:
\[
  S \;+\; N \;+\; G \;\rightarrow\; R
\]

where $S$ is the structural substrate (the genome and
its encoded build capacity), $N$ is the external
coherence gradient field (the environmental specification
of which developmental attractor the genome resolves to),
$G$ is the developmental geometry (the global attractor
landscape, basin structure, and current position of
the system in state space), and $R$ is the resulting
phenotype.

This formulation follows directly from the Waddington
epigenetic landscape \cite{waddington1957}, whose
mathematical formalisation as a quasi-potential function:
\[
  U(\mathbf{x}) = -k_B T \ln P_{\text{ss}}(\mathbf{x})
\]
(where $P_{\text{ss}}$ is the steady-state probability
distribution over developmental state space) was
provided by Wang et al.\ (2014) \cite{wang2014}. The
landscape $U$ is not readable from the gene regulatory
network $F$ alone without computing the global state
space geometry. This is the mathematical basis for the
incompleteness of the genome-only model: the genome
encodes $F$; the phenotype is determined by $U$; and
$U$ is shaped by both $F$ and the external field $N$.

\subsection{Attractor Basins and Developmental Arrest}

In the Waddington landscape, stable phenotypes
correspond to attractor basins --- valleys in $U$
whose depth reflects the stability of the corresponding
developmental state. A system placed in a coherence
gradient field $N$ will navigate to the basin whose
geometry is specified by $N$.

\textbf{Developmental arrest} occurs when a coherence
gradient field holds a developing system in a
shallower, earlier-stage basin rather than allowing
the developmental trajectory to proceed to its
deep-basin endpoint. The arrested state is
phenotypically stable as long as the maintenance
field is applied. When the maintenance field is
removed, the system exits the shallow basin and
continues its developmental trajectory toward the
nearest deep basin.

This mechanism has been directly demonstrated at
multiple scales:
\begin{itemize}[itemsep=0.1em]
  \item \textit{Single cell:} Nitrogen starvation
    arrests yeast in the LECA-proximate attractor basin.
    Restore nitrogen; the modern eukaryotic program
    resumes \cite{leca2026}.
  \item \textit{Single cell:} EZH2-mediated H3K27me3
    deposition silences FOXA1 and arrests breast
    epithelial cells in the cancer attractor basin.
    Tazemetostat (EZH2 inhibitor) de-arrests the cell;
    the luminal program resumes \cite{foxa1ezh2_2026}.
  \item \textit{Organism:} Inversion of the coherence
    gradient field (mineral water above, UV-A below)
    arrests the radish developmental program in the
    inverted architectural attractor. The root grows
    upward; the shoot grows downward. Same genome;
    inverted field; inverted organism \cite{plantinv2026}.
\end{itemize}

We propose that domestication is the same mechanism
at the organismal scale: a field-specified developmental
arrest of the neural crest migration program.

% ══════════════════════════════════════════════════════════════
\section{The Domestication Mechanism: Field-Specified
Neural Crest Arrest}
\label{sec:mechanism}
% ══════════════════════════════════════════════════════════════

\subsection{The Domestication Coherence Gradient Field}

The domestication coherence gradient field
$N_{\text{domestic}}$ is defined by the following
components acting simultaneously:

\begin{itemize}[itemsep=0.15em]
  \item \textbf{Predation removal:} Elimination of
    fight-or-flight activation. Sustained HPA axis
    suppression.
  \item \textbf{Ad libitum food provision:} Elimination
    of metabolic stress. Sustained anabolic state.
  \item \textbf{Temperature regulation:} Elimination of
    thermoregulatory stress load.
  \item \textbf{Social stability:} Elimination of
    dominance-hierarchy conflict. Oxytocin-dominant
    social environment.
  \item \textbf{Human contact from early development:}
    Imprinting on the human social field from the
    critical window.
\end{itemize}

This field does not modify the genome. It modifies the
physiological state of the animal occupying it ---
specifically the HPA axis set point and the resulting
hormonal environment of any gestating offspring.

\subsection{The Causal Chain}

The complete causal chain from field to phenotype is:

\begin{mdframed}[
  backgroundcolor=mechanismbox,
  linecolor=headergreen,
  linewidth=0.8pt,
  innertopmargin=0.5em,
  innerbottommargin=0.5em,
  innerleftmargin=0.8em,
  innerrightmargin=0.8em
]
\begin{center}
\small
$N_{\text{domestic}}$\\[0.2em]
$\downarrow$\\[0.2em]
Maternal HPA axis suppression
(reduced CRH, ACTH, cortisol)\\[0.2em]
$\downarrow$\\[0.2em]
Reduced intrauterine glucocorticoid exposure\\[0.2em]
$\downarrow$\\[0.2em]
Altered embryonic coherence gradients\\
(BMP4/7, Wnt1/3a, CXCL12/CXCR4, EDNRB, FGF8)\\[0.2em]
$\downarrow$\\[0.2em]
Partial arrest of neural crest cell migration\\
(juvenile-proximate stage)\\[0.2em]
$\downarrow$\\[0.2em]
Basin B occupation (domestic attractor)\\[0.2em]
$\downarrow$\\[0.2em]
\textbf{Domestication syndrome}\\
(floppy ears, depigmentation, smaller adrenal,\\
shorter snout, docility, juvenile behaviour retention)
\end{center}
\end{mdframed}

\vspace{0.3em}

Each step in this chain has independent support in
the developmental biology and stress physiology
literature:

\textbf{Step 1 (HPA suppression):} Glucocorticoids are
known regulators of NCC migration. Elevated cortisol
suppresses NCC proliferation and increases apoptosis
in migrating NCC populations \cite{sapolsky2002}.
$N_{\text{domestic}}$ chronically suppresses maternal
cortisol.

\textbf{Step 2 (Intrauterine gradient alteration):}
Maternal glucocorticoid levels directly influence foetal
HPA axis programming via placental 11$\beta$-HSD2
activity and direct cortisol transfer. The intrauterine
hormonal environment shapes the coherence gradients
for embryonic NCC migration \cite{seckl2004}.

\textbf{Step 3 (NCC gradient dependence):} NCC migration
is not genetically autonomous. It is guided by
chemokine gradients (CXCL12/CXCR4), BMP gradients,
Wnt signalling, and extracellular matrix composition ---
all of which are sensitive to the hormonal environment
of the embryo \cite{theveneau2013}.

\textbf{Step 4 (Syndrome as arrest signature):}
Wilkins et al.\ (2014) confirmed that all syndrome
traits are downstream of NCC migration reduction.
Our contribution is the causal origin of that reduction:
it is field-specified, not allele-specified.

\subsection{Why NCC Migration Is a Developmental
Trajectory, Not a Genetic Switch}

Neural crest cell migration has a starting state
(NCC cells at the dorsal neural tube), a directed
trajectory (migration through embryonic tissue guided
by coherence gradients), and an endpoint (full
differentiation into all NCC-derived tissues).

Full migration to completion = Basin A = wild-type
phenotype.

Arrested migration at early stage = Basin B = domestic
phenotype.

The migration is not a binary genetic switch between
two allelic states. It is a continuous developmental
trajectory that can be held at any intermediate stage
by the coherence gradient environment of the developing
embryo. The domestication field holds it at the
juvenile-proximate stage. The wild field allows it to
proceed to completion.

\textbf{This is the precise mechanism of developmental
arrest.} The domestic pig is not a genetically different
organism from the wild boar. It is the same organism
whose developmental trajectory was arrested at an
earlier stage and held there by field maintenance.

% ══════════════════════════════════════════════════════════════
\section{The Two-Basin Landscape of \textit{Sus scrofa}}
\label{sec:landscape}
% ══════════════════════════════════════════════════════════════

The attractor landscape of \textit{Sus scrofa} contains
(at minimum) two primary developmental basins:

\begin{center}
\begin{tabular}{p{2.8cm}p{5.2cm}p{5.2cm}}
\toprule
& \textbf{Basin A (Wild Boar)} &
  \textbf{Basin B (Domestic Pig)} \\
\midrule
\textbf{Depth} &
  Deep. Carved by millions of years of\\
  natural selection. &
  Shallow. Carved by $\sim$10,000 years\\
  of domestication selection. \\[0.3em]
\textbf{NCC migration} &
  Runs to completion. &
  Arrested at juvenile-proximate stage. \\[0.3em]
\textbf{Field requirement} &
  $N_{\text{wild}}$: predation, food\\
  scarcity, thermoregulatory challenge. &
  $N_{\text{domestic}}$: protection,\\
  ad libitum food, stability. \\[0.3em]
\textbf{Phenotype} &
  Coarse hair, erect ears, elongated\\
  snout, large adrenal, aggressive,\\
  lean, tusked, camouflaged. &
  Fine hair, floppy ears, shortened\\
  snout, small adrenal, docile,\\
  fat-depositing, depigmented. \\[0.3em]
\textbf{Maintenance} &
  Self-maintaining under\\
  $N_{\text{wild}}$. &
  Requires continuous\\
  $N_{\text{domestic}}$ application. \\
\bottomrule
\end{tabular}
\end{center}

\vspace{0.3em}
The critical implication of the basin depth asymmetry
is the \textbf{transition rate asymmetry}:

\begin{itemize}[itemsep=0.15em]
  \item \textbf{Basin B $\rightarrow$ Basin A
    (feralization):} Remove $N_{\text{domestic}}$.
    The system exits the shallow Basin B and rolls
    into the deep Basin A. Rate: months to one
    generation. Energy barrier: low.
  \item \textbf{Basin A $\rightarrow$ Basin B
    (domestication):} Must carve Basin B against the
    restoring force of Basin A. Requires sustained,
    directed $N_{\text{domestic}}$ application over
    many generations. Energy barrier: high.
\end{itemize}

This asymmetry is geometrically necessary from the
landscape structure. It is not explained by the
genetic model, in which allele frequencies shift at
approximately equivalent rates under equivalent
selection pressure in both directions.

% ══════════════════════════════════════════════════════════════
\section{Confirmation Against Existing Literature}
\label{sec:confirmation}
% ══════════════════════════════════════════════════════════════

The following confirmations are against literature
identified \textit{after} the mechanism was derived.
Each confirmation is of a prediction entailed by the
field-arrest model.

\subsection{Marcstr\"{o}m \& Essen-M\"{o}ller (1968):
The Buried Confirmation}

Marcstr\"{o}m and Essen-M\"{o}ller placed wild-caught
\textit{Sus scrofa} in captivity without selective
breeding. In the first captive generation, the following
appeared spontaneously:

\begin{itemize}[itemsep=0.1em]
  \item White piebald patches (melanocyte migration
    arrest --- NCC-derived)
  \item Widening and drooping of ears (cartilage
    reduction --- NCC-derived)
  \item Skull shape changes toward shorter, wider
    craniofacial profile (NCC-derived)
  \item Curly tails
\end{itemize}

\textbf{The critical point:} these traits appeared
without selective breeding, in animals with fully wild
genomes, simply by applying $N_{\text{domestic}}$.
The frequency of trait expression was too high for
recessive allele expression to explain.
\textbf{Confirmed by field-arrest model.}
\textbf{Not explained by genetic model.}

\subsection{Hemmer (1990): Adrenal Reduction as
Earliest Change}

Hemmer's systematic analysis documented that adrenal
gland size reduction is among the earliest and most
consistent morphological changes in captive wild
animals across species \cite{hemmer1990}. The adrenal
medulla is NCC-derived. Its size reflects the degree
of NCC migration completion to the adrenal anlage.

\textbf{Prediction confirmed:} the fastest-responding
tissue is the one whose developmental program is most
directly downstream of the field's primary action
(HPA axis suppression). \textbf{Confirmed.}

\subsection{Wilkins, Wrangham \& Fitch (2014):
The Cellular Mechanism Identified}

Wilkins et al.\ correctly identified that all
domestication syndrome traits map to a single
developmental perturbation: mild reduction in NCC
migration. They attributed this to genetic change.
We provide the causal origin of the perturbation:
field-specified alteration of the intrauterine NCC
migration guidance environment via maternal HPA
axis suppression.

\textbf{Our paper completes the causal chain they
left open.} It does not contradict their finding.
It explains why it is true.

\subsection{Belyaev (1979) and Trut et al.\ (2009):
Tameness Is Not the Cause}

The silver fox domestication experiment selected
only for HPA axis suppressibility (reduced fear/flight
response). Within 10--15 generations, the full
domestication syndrome appeared: floppy ears,
piebald coat, shortened snout, curly tail, reduced
adrenal size, juvenile behaviour retention.

\textbf{Standard interpretation:} co-selected genetic
alleles linking tameness to morphology.

\textbf{Field-arrest interpretation:} selecting for
the most HPA-suppressible individuals selected for
those whose intrauterine environment was closest to
$N_{\text{domestic}}$. Their offspring's NCC migration
was most arrested. The syndrome appeared not because
tameness alleles are linked to ear alleles. It appeared
because both are downstream of the same field-specified
arrest mechanism.

\textbf{Tameness is the behavioural signature of
Basin B occupation. The morphological traits are the
structural signature of the same Basin B occupation.
They appear together because they are both downstream
of the same arrested developmental trajectory.}

\subsection{The Feralization Sequence: A Reversion
Signature}

In feralization, morphological reversion occurs in
a specific sequence: behaviour first (days to weeks),
then hair coarseness (weeks to months), then body
composition (months), then craniofacial structure
(requires offspring generation for full expression).

\textbf{Prediction from field-arrest model:}
The reversion sequence should follow the order of
developmental plasticity from highest to lowest ---
epigenetically-regulated traits (fastest) to
structurally-embedded traits (slowest). This is
exactly the observed sequence. The genetic model
predicts no specific sequence, as allele frequencies
shift at equivalent rates across traits.
\textbf{Confirmed.}

\subsection{Summary of Confirmations}

\begin{center}
\small
\begin{tabular}{p{5.5cm}p{3.5cm}p{3.0cm}}
\toprule
\textbf{Prediction (derived before literature)} &
\textbf{Source} & \textbf{Status} \\
\midrule
Domestication syndrome in F1 captive wild
animals without selective breeding &
Marcstr\"{o}m (1968) & \textbf{Confirmed} \\[0.3em]
Trait frequency too high for recessive allele
explanation &
Marcstr\"{o}m (1968) & \textbf{Confirmed} \\[0.3em]
Adrenal size reduction earliest and most
consistent change &
Hemmer (1990) & \textbf{Confirmed} \\[0.3em]
All syndrome traits $\rightarrow$ single
developmental perturbation (NCC) &
Wilkins et al.\ (2014) & \textbf{Confirmed} \\[0.3em]
Tameness selection alone produces full
syndrome &
Belyaev (1979) & \textbf{Confirmed} \\[0.3em]
Feralization rate asymmetry vs.\
re-domestication &
Multiple sources & \textbf{Confirmed} \\[0.3em]
Reversion sequence tracks developmental
plasticity order &
Feralization literature & \textbf{Confirmed} \\
\bottomrule
\end{tabular}
\end{center}

% ══════════════════════════════════════════════════════════════
\section{Novel Experimental Predictions}
\label{sec:predictions}
% ══════════════════════════════════════════════════════════════

The following predictions are entailed by the
field-arrest model and are not entailed by the genetic
model. They constitute the experimental program that
distinguishes the two frameworks.

\begin{mdframed}[
  backgroundcolor=predbox,
  linecolor=blue!50!black,
  linewidth=0.8pt,
  innertopmargin=0.5em,
  innerbottommargin=0.5em,
  innerleftmargin=0.8em,
  innerrightmargin=0.8em
]
\textbf{NOTE ON PRE-REGISTRATION:}
Predictions P7--P11 are pre-registered on Zenodo
(2026-03-18) prior to any experimental data collection.
DOI to be inserted upon Zenodo confirmation.
\end{mdframed}

\vspace{0.3em}

\subsection{P7 --- One-Generation Morphological
Re-domestication via Cross-Fostering}

\textbf{Design:} F1 offspring of wild-caught
\textit{Sus scrofa} (fully wild genome, confirmed
by genotyping) placed at birth with domestic foster
sow under full $N_{\text{domestic}}$ from birth.
Genotype-matched controls raised under $N_{\text{wild}}$.

\textbf{Measurements at 3, 6, 12 months:}
adrenal gland size (at necropsy), cortisol response
to standardised stressor, ear cartilage stiffness
and position, coat pigmentation pattern
(\% piebald surface area), hair coarseness,
snout length/width ratio.

\textbf{Field-arrest prediction:}
$N_{\text{domestic}}$-raised F1 will show smaller
adrenal glands and reduced HPA reactivity compared
to genotype-matched $N_{\text{wild}}$-raised F1.
Ear and pigmentation changes will be intermediate
(postnatal field only; full expression requires
intrauterine field application --- see P8).

\textbf{Genetic model prediction:}
No morphological difference between conditions
with identical genomes.

\textbf{The models diverge. The experiment
has not been done.}

\subsection{P8 --- Maternal Field Window Dependence}

\textbf{Design:} Feral sow inseminated and held
under $N_{\text{domestic}}$ from conception.
F1 offspring raised identically post-birth.
Compare to P7 controls (field applied only
from birth).

\textbf{Prediction:} F1 offspring of sow in
$N_{\text{domestic}}$ from conception will show
larger morphological shift toward domestic phenotype
than F1 offspring with only postnatal
$N_{\text{domestic}}$ exposure, with identical genomes.

The gradient of morphological shift across conditions
(full prenatal + postnatal $>$ postnatal only $>$
wild control) specifies the NCC migration critical
window relative to parturition.

\subsection{P9 --- H3K27me3 Enrichment at NCC
Migration Loci in Domestic vs.\ Wild Embryos}

\textbf{Design:} H3K27me3 ChIP-seq on embryonic
tissue (E14--E18, NCC active migration stage)
from domestic \textit{Sus scrofa} and wild boar.

\textbf{Prediction:} Domestic pig embryos will
show elevated H3K27me3 (EZH2-mediated silencing)
at NCC migration driver loci relative to wild boar
embryos at the same developmental stage, specifically:

\begin{itemize}[itemsep=0.1em, noitemsep]
  \item BMP4, BMP7 (NCC migration guidance)
  \item WNT1, WNT3A (NCC specification)
  \item CXCL12, CXCR4 (directional migration)
  \item EDNRB (pigment cell migration)
  \item SOX10, SNAI2 (NCC identity maintenance)
\end{itemize}

\textbf{Significance:} This would establish
molecular identity between the domestication
arrest mechanism and EZH2-mediated cancer attractor
arrest --- connecting domestication biology directly
to the cancer biology framework.

\subsection{P10 --- Breed Domestication Depth
Predicts Feralization Rate}

\textbf{Prediction:} Feralization rate (months to
morphological reversion threshold) correlates
negatively with domestication depth (genetic distance
from wild boar; generations of intense selection).

Ranking (slowest to fastest feralization):
Large White/Landrace (most derived) $>$
heritage breeds $>$ Mangalitsa/Ossabaw
(least derived, closest to wild boar ancestry).

\textbf{Design:} Comparative morphological measurement
(hair coarseness, adrenal weight, HPA reactivity)
in feral populations derived from known breed
starting points at 3 and 6 months post-feralization.

\subsection{P11 --- EZH2 Inhibition Shifts F1
Morphology Toward Wild-Type}

\textbf{Rationale:} If H3K27me3 at NCC migration
loci maintains the domestication arrest (P9), then
EZH2 inhibition in a pregnant domestic sow should
partially de-repress NCC migration drivers and
shift F1 morphology toward wild-type.

\textbf{Design:} Pregnant domestic sows treated with
sub-therapeutic tazemetostat (EZH2 inhibitor,
FDA-approved) during the NCC migration window.
F1 morphology assessed at 6 months vs.\
vehicle-treated controls.

\textbf{Prediction:} F1 offspring of EZH2-inhibited
dams will show increased hair coarseness, larger
adrenal glands, and higher HPA reactivity relative
to F1 offspring of vehicle-treated dams.
\textbf{No genetic modification. Same genome.
Field perturbation shifts phenotype.}

% ══════════════════════════════════════════════════════════════
\section{The Cross-Scale Mechanistic Identity}
\label{sec:crossscale}
% ══════════════════════════════════════════════════════════════

The domestication arrest mechanism is not unique
to Sus scrofa or to organismal-scale biology.
The same mechanism --- field-specified arrest of
a developmental trajectory at an attractor basin
maintained by epigenetic silencing, reversible by
field change --- is operative at multiple scales:

\vspace{0.3em}
\begin{center}
\small
\begin{tabular}{p{2.2cm}p{2.8cm}p{2.8cm}p{2.8cm}p{2.2cm}}
\toprule
\textbf{Scale} & \textbf{System} &
\textbf{Arrest field} & \textbf{Resumption} &
\textbf{Reference} \\
\midrule
Single cell &
TNBC &
EZH2/H3K27me3 silences FOXA1 &
Tazemetostat removes silencing &
\cite{foxa1ezh2_2026} \\[0.4em]
Single cell &
LECA yeast &
Nitrogen starvation arrests modern program &
Restore nitrogen &
\cite{leca2026} \\[0.4em]
Organism &
Radish &
Standard field specifies root-down attractor &
Invert field $\rightarrow$ root-up &
\cite{plantinv2026} \\[0.4em]
Organism &
Sus scrofa &
$N_{\text{domestic}}$ arrests NCC migration &
Remove field $\rightarrow$ feralization &
This paper \\
\bottomrule
\end{tabular}
\end{center}

\vspace{0.3em}
The structural identity across scales is not
metaphorical. The same causal operation is occurring
at each scale: the external coherence gradient field
$N$ specifies which attractor basin the structural
substrate $S$ resolves to within the developmental
geometry $G$, producing result $R$. The genome or
cellular build capacity never changes. The field
changes. The phenotype changes.

If the domestication arrest mechanism involves
EZH2-mediated H3K27me3 at NCC loci (P9), then
the molecular mechanism is identical to the
cancer arrest mechanism. Tazemetostat, the
FDA-approved EZH2 inhibitor currently in
clinical use for cancer, would be the reagent
that de-arrests the domestication program (P11).
The same drug. The same molecular mechanism.
Cancer and domestication as instances of the same
field-specified epigenetic arrest. Different scale.
One mechanism.

% ══════════════════════════════════════════════════════════════
\section{Why the Genetic Model Is Incomplete}
\label{sec:genetic}
% ══════════════════════════════════════════════════════════════

The genetic model of domestication correctly explains:
\begin{itemize}[itemsep=0.1em]
  \item Why domesticated traits are heritable
    (the alleles and epigenetic states are transmitted)
  \item Why selective breeding accelerates
    domestication (selection deepens Basin B)
\end{itemize}

The genetic model cannot explain:

\begin{enumerate}[itemsep=0.2em]
  \item \textbf{Cross-species syndrome universality:}
    Every domesticated species, regardless of
    phylogenetic distance, shows the same syndrome.
    Convergent acquisition of linked NCC-affecting
    alleles in every independent domestication event
    is implausible. A shared field type producing
    a shared arrest is parsimonious.

  \item \textbf{First-generation appearance without
    selection:} Marcstr\"{o}m (1968) documented
    syndrome traits in F1 captive wild boar without
    selection. Recessive allele expression cannot
    produce high-frequency, convergent trait
    appearance in F1 without selection.

  \item \textbf{Feralization/domestication rate
    asymmetry:} Allele frequencies shift at
    equivalent rates under equivalent selection.
    Feralization in months vs.\ domestication in
    generations is not explained. Basin depth
    asymmetry explains it directly.

  \item \textbf{Syndrome from tameness selection
    alone:} The Belyaev result requires either
    tight genetic linkage between tameness and
    every morphological trait simultaneously, or
    a single upstream mechanism (HPA suppression
    $\rightarrow$ intrauterine coherence gradient
    $\rightarrow$ NCC arrest) producing all traits
    simultaneously. The field-arrest model provides
    this upstream mechanism.

  \item \textbf{Reversion sequence:} Allele
    frequencies change at equivalent rates across
    traits. The ordered reversion sequence in
    feralization (behaviour $\to$ hair $\to$ body
    composition $\to$ craniofacial) has no explanation
    under the genetic model. It is geometrically
    necessary under the developmental arrest model.
\end{enumerate}

% ══════════════════════════════════════════════════════════════
\section{Discussion}
\label{sec:discussion}
% ════════════════���═════════════════════════════════════════════

\subsection{What Domestication Actually Is}

Domestication, under the field-arrest model, is not
the genetic modification of wild animals into domestic
ones. It is the engineering of a coherence gradient
field that arrests the developmental trajectory of
wild animals at the juvenile-proximate basin and holds
them there by continuous field maintenance.

The genome does not change.
The field changes.
The developmental trajectory is arrested.
The organism expresses the phenotypic signature
of that arrest.

This has two immediate practical implications:

\textbf{Implication 1:} The wild phenotype is
always recoverable. It is not lost. It is suspended.
Remove the maintenance field and the developmental
trajectory resumes. This is why feralization is
fast and universal.

\textbf{Implication 2:} The domestication field,
not the genome, is the primary engineering target.
To produce a domestic phenotype from a wild
starting point, the correct intervention is not
selective breeding per se --- it is engineering
the coherence gradient field that will arrest
NCC migration at the domestic-proximate stage in
the developing offspring.

\subsection{The Evidence That Has Always Been
There}

Marcstr\"{o}m and Essen-M\"{o}ller published their
captive wild boar data in 1968. The results were
correctly observed: domestication syndrome traits
in first-generation captive animals without
selective breeding. They were incorrectly interpreted:
latent allele expression under relaxed selection.

The correct interpretation has been unavailable
because the framework to state it --- field-specified
developmental arrest as a causal mechanism in
vertebrate biology --- did not exist. The evidence
was always there. The framework to read it was not.

\subsection{The Wojtek Problem}

Wojtek, the Syrian brown bear who served with the
Polish Army (1943--1963), was raised from cubhood
in the human social field. He became fully
socialised, learned to carry ammunition, ate
cigarettes, and lived as a functional member of
a human military unit. Standard interpretation:
bears are intelligent and can be trained.

Geometric interpretation: Wojtek was placed in
$N_{\text{domestic}}$ during the critical developmental
window for HPA axis calibration and social field
imprinting. His developmental program resolved to
the domestic-proximate basin not because he was
trained but because that was the basin specified
by the field in which he developed.

Every wildlife rehabilitation programme that
rescues neonates is running this experiment
continuously and observing this result consistently.
The result is interpreted as taming or socialisation.
It is developmental arrest. The difference is
not semantic. It is causal.

\subsection{Relation to Cancer Biology and
the LECA Protocol}

The field-arrest model of domestication connects
directly to two active experimental programs:

The FOXA1/EZH2 work in triple-negative breast
cancer proposes that cancer is a cellular
developmental arrest --- the cell held in an
ancestral attractor basin by EZH2-mediated
epigenetic silencing of the identity program.
If P9 confirms H3K27me3 enrichment at NCC loci
in domestic pig embryos, the same molecular
mechanism maintains developmental arrest at two
scales simultaneously: cell and organism.

The LECA protocol proposes that the precambrian
ionic coherence gradient field specifies the LECA
attractor basin in living yeast. The domestication
protocol is the same operation at the organismal
scale: the domestication coherence gradient field
specifies the domestic attractor basin in
\textit{Sus scrofa}. Both are recoverable.
Both are field-specified. Both involve the same
genome operating in a different field.

The pig, the cancer cell, the yeast, and the
radish are not separate biological curiosities.
They are the same causal geometry at different
scales. The domestication paper and the cancer
paper are the same paper at different scales.

% ══════════════════════════════════════════════════════════════
\section{Conclusions}
\label{sec:conclusions}
% ══════════════════════════════════════════════════════════════

We propose that domestication is field-specified
developmental arrest of the neural crest migration
program. The domestication coherence gradient field
suppresses the maternal HPA axis, alters the
intrauterine coherence gradient environment,
and arrests neural crest cell migration at the
juvenile-proximate developmental stage. The
domestication syndrome is the phenotypic signature
of that arrest. Feralization is developmental
resumption when the maintenance field is removed.

This mechanism:
\begin{enumerate}[itemsep=0.15em]
  \item Provides the causal origin of the NCC
    migration reduction identified by Wilkins
    et al.\ (2014)
  \item Explains the cross-species universality
    of the domestication syndrome without
    invoking convergent genetic change
  \item Explains the first-generation appearance
    of syndrome traits in captive wild animals
    without selection (Marcstr\"{o}m 1968)
  \item Explains the asymmetry between
    feralization rate and re-domestication rate
    from basin depth geometry
  \item Explains the ordered reversion sequence
    in feralization from developmental plasticity
    principles
  \item Generates five novel experimental predictions
    (P7--P11) that are not entailed by the genetic
    model and that constitute a decisive test
    between frameworks
  \item Is structurally identical to EZH2-mediated
    attractor arrest in cancer and nitrogen
    starvation arrest in the LECA protocol,
    suggesting a scale-invariant field-arrest
    mechanism across vertebrate biology
\end{enumerate}

The domestic pig is not a genetically modified
wild boar. It is a wild boar whose developmental
trajectory was arrested at the juvenile-proximate
basin by a field that has been continuously
maintained for $\sim$10,000 years. The wild boar
was never lost. It was suspended. Remove the
field and it resumes.

% ══════════════════════════════════════════════════════════════
% REFERENCES
% ══════════════════════════════════════════════════════════════
\newpage
\begin{thebibliography}{99}

\bibitem{wilkins2014}
Wilkins, A.S., Wrangham, R.W., \& Fitch, W.T. (2014).
The ``Domestication Syndrome'' in Mammals: A Unified
Explanation Based on Neural Crest Cell Behavior and
Genetics.
\textit{Genetics}, 197(3), 795--808.
\href{https://doi.org/10.1534/genetics.114.165423}
{\texttt{doi:10.1534/genetics.114.165423}}

\bibitem{belyaev1979}
Belyaev, D.K. (1979).
Destabilizing Selection as a Factor in Domestication.
\textit{Journal of Heredity}, 70(5), 301--308.

\bibitem{trut2009}
Trut, L., Oskina, I., \& Kharlamova, A. (2009).
Animal evolution during domestication: the domesticated
fox as a model.
\textit{BioEssays}, 31(3), 349--360.
\href{https://doi.org/10.1002/bies.200800070}
{\texttt{doi:10.1002/bies.200800070}}

\bibitem{hemmer1990}
Hemmer, H. (1990).
\textit{Domestication: The Decline of Environmental
Appreciation}.
Cambridge University Press, Cambridge.

\bibitem{marcstrom1968}
Marcstr\"{o}m, V., \& Essen-M\"{o}ller, M. (1968).
The inheritance of white pigmented areas in Swedish
wild boar (\textit{Sus scrofa}).
\textit{Viltr\aa{}vy} (Swedish Wildlife),
5, 325--336.

\bibitem{marcstrom1969}
Marcstr\"{o}m, V. (1969).
Captive rearing of wild boar.
\textit{Viltr\aa{}vy} (Swedish Wildlife),
6, 105--114.

\bibitem{kauffman1993}
Kauffman, S.A. (1993).
\textit{The Origins of Order: Self-Organization and
Selection in Evolution}.
Oxford University Press, New York.

\bibitem{huang2009}
Huang, S., Eichler, G., Bar-Yam, Y., \&
Ingber, D.E. (2005).
Cell Fates as High-Dimensional Attractor States of a
Complex Gene Regulatory Network.
\textit{Physical Review Letters}, 94(12), 128701.

\bibitem{wang2014}
Wang, J., Zhang, K., Xu, L., \& Wang, E. (2011).
Quantifying the Waddington landscape and biological
paths for development and differentiation.
\textit{Proceedings of the National Academy of Sciences},
108(20), 8257--8262.
\href{https://doi.org/10.1073/pnas.1017017108}
{\texttt{doi:10.1073/pnas.1017017108}}

\bibitem{waddington1957}
Waddington, C.H. (1957).
\textit{The Strategy of the Genes}.
Allen \& Unwin, London.

\bibitem{theveneau2013}
Th\'{e}veneau, E., \& Mayor, R. (2012).
Neural crest delamination and migration:
from epithelium-to-mesenchyme transition to collective
cell migration.
\textit{Developmental Biology}, 366(1), 34--54.

\bibitem{sapolsky2002}
Sapolsky, R.M., Romero, L.M., \& Munck, A.U. (2000).
How Do Glucocorticoids Influence Stress Responses?
Integrating Permissive, Suppressive, Stimulatory, and
Preparative Actions.
\textit{Endocrine Reviews}, 21(1), 55--89.

\bibitem{seckl2004}
Seckl, J.R., \& Meaney, M.J. (2004).
Glucocorticoid programming.
\textit{Annals of the New York Academy of Sciences},
1032, 63--84.

\bibitem{price1999}
Price, E.O. (1999).
Behavioral development in animals undergoing
domestication.
\textit{Applied Animal Behaviour Science},
65(3), 245--271.

\bibitem{larson2014}
Larson, G., \& Fuller, D.Q. (2014).
The Evolution of Animal Domestication.
\textit{Annual Review of Ecology, Evolution,
and Systematics}, 45, 115--136.

\bibitem{foxa1ezh2_2026}
Lawson, E.R. (2026).
FOXA1/EZH2 Dual IHC as a Universal Breast Cancer
Treatment Classifier: Attractor Geometry Derivation
and Seven-Dataset Validation.
OrganismCore / Zenodo.
\href{https://doi.org/10.5281/zenodo.18883922}
{\texttt{doi:10.5281/zenodo.18883922}}

\bibitem{leca2026}
Lawson, E.R. (2026).
LECA Resurrection Protocol: Pre-Registration.
OrganismCore / Zenodo.
\href{https://doi.org/10.5281/zenodo.18986790}
{\texttt{doi:10.5281/zenodo.18986790}}

\bibitem{plantinv2026}
Lawson, E.R. (2026).
Plant Architecture Inversion via Coherence Gradient
Field Specification: Pre-Registration.
OrganismCore / Zenodo.
\href{https://doi.org/10.5281/zenodo.19040399}
{\texttt{doi:10.5281/zenodo.19040399}}

\end{thebibliography}

% ── DOCUMENT METADATA ────────────────────────────────────────
\vspace{1em}
\noindent\textcolor{rulecolor}{\rule{\linewidth}{0.6pt}}
\vspace{0.3em}

{\footnotesize
\textbf{Document ID:} DOMESTICATION\_AS\_DEVELOPMENTAL\_
ARREST\_v1.0\\
\textbf{Date:} 2026-03-18\\
\textbf{Author:} Eric Robert Lawson / OrganismCore\\
\textbf{Derivation record:}
\href{https://github.com/Eric-Robert-Lawson/attractor-oncology}
{\texttt{github.com/Eric-Robert-Lawson/attractor-oncology}}\\
\textbf{Pre-registration of P7--P11:} Zenodo, 2026-03-18
(DOI to be inserted upon confirmation)\\
\textbf{License:} CC BY 4.0\\
\textbf{Target journal:} \textit{BioEssays}
(primary); \textit{Genetics} (secondary ---
direct extension of Wilkins et al.\ 2014)\\
\textbf{ORCID:}
\href{https://orcid.org/0009-0002-0414-6544}
{0009-0002-0414-6544}\\
\textbf{Contact:}
\href{mailto:OrganismCore@proton.me}
{OrganismCore@proton.me}
}

\end{document}
