# siRNA Design Studio

**Design, evaluation and laboratory preparation of small interfering RNA —
with the assumptions made explicit.**

A pipeline that takes a target sequence and returns ranked siRNA candidates,
their off-target risk assessed against a host transcriptome, their
conservation across pathogen isolates, and every oligonucleotide form needed
to order and test them at the bench.

[![DOI](https://zenodo.org/badge/1336237502.svg)](https://doi.org/10.5281/zenodo.22061890)
[![Tests](https://img.shields.io/badge/tests-110%20passing-brightgreen)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://sirna-design-studio.streamlit.app)


<!-- Add a screenshot:
![Interface](docs/screenshot.png)
-->

---

## Table of contents

- [The problem](#the-problem)
- [What the application does](#what-the-application-does)
- [What distinguishes it](#what-distinguishes-it)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Method](#method)
- [Modules](#modules)
- [Tests](#tests)
- [Comparison with existing tools](#comparison-with-existing-tools)
- [Limitations](#limitations)
- [References](#references)
- [Citation](#citation)

---

## The problem

RNA interference is a mature technique with an immature design process. A
21-nucleotide sequence complementary to a target will silence it — sometimes.
Whether it does depends on the thermodynamic asymmetry of the duplex, the
accessibility of the target site, the identity of the 5-terminal nucleotide,
the presence of unintended complementarity elsewhere in the transcriptome,
and, for a pathogen target, on whether the site is conserved across the
isolates one actually expects to encounter.

Existing design tools compute most of these. What they generally do not do is
say how the individual criteria were weighed against each other, whether the
resulting ranking survives a change in that weighting, or how a count of
off-target hits should be read against the number expected by chance in a
database of that size. A candidate list without those answers is difficult to
defend and difficult to reproduce.

This pipeline was written to make each of those decisions explicit and
inspectable.

---

## What the application does

### Conservation analysis

Given a set of isolate sequences, the application determines which regions
are present in the largest fraction of them. An siRNA directed at a conserved
region acts on the whole population and leaves the pathogen fewer routes of
escape by mutation.

Rather than computing Shannon entropy from a multiple sequence alignment, it
measures **exact k-mer coverage**: for each window in the reference, the
fraction of remaining isolates containing an identical string. The reasoning
is biological — an siRNA requires near-perfect complementarity, and a single
mismatch in the seed region is enough for the target not to be recognised.
What matters is therefore not conservation in the evolutionary sense but the
presence of one specific k-mer in a given isolate, which is measured directly
and needs no alignment software.

The coverage threshold is not assumed. The number of surviving windows plotted
against the threshold has a characteristic shape for a population with a
conserved core: a slow decline while variable windows are removed, then a
sharp drop. The point of maximum curvature marks the transition and is located
automatically. Where no such transition exists, the application says so rather
than returning an authoritative-looking number.

### Candidate design

Windows of 21, 22 and 24 nucleotides are generated and filtered according to
the selected host profile. The three lengths are not interchangeable: in
plants they enter different pathways — DCL4 for post-transcriptional silencing,
DCL2 for transitivity and systemic movement, DCL3 for RNA-directed DNA
methylation — and the application ranks them separately when asked to.

For each candidate it computes duplex free energy by the nearest-neighbour
method, the melting temperature, the folding energy of the guide strand alone,
the accessibility of the target site in its local mRNA context, and the
thermodynamic asymmetry that determines which strand the effector complex will
retain.

### Off-target assessment

The guide sequences are scanned against a host transcriptome in a single
streaming pass. Two classes of risk are distinguished: seed-region binding of
the kind that drives miRNA-like repression, and full-length matches capable of
directing cleavage. Mismatches are weighted by position, since a mismatch at
position 3 and one at position 19 are not biologically equivalent.

Hits in the intended target gene are recognised and excluded — without this,
every candidate would be rejected whenever the target is an endogenous gene.

### Construct assembly

Expression cassettes are assembled for a chosen promoter and cloning method.
The application knows that Pol III promoters constrain the first transcribed
nucleotide and require an internal poly-T terminator, that Golden Gate fails
if the insert contains the enzyme's own recognition site, and that Gateway
attB sites are longer than a short hairpin. Custom promoters and vectors can
be defined.

### Laboratory forms

A ranked sequence is not something that can be ordered. The application
generates the synthetic duplex with 3-prime overhangs, the shRNA hairpin,
annealing oligonucleotides for ligation, T7 templates for in vitro
transcription, stem-loop RT primers and a Northern probe for detection, and
four classes of specificity control.

### Reporting

An in silico synthesis report records the input with a checksum, every
parameter including those left at default, the metrics of each selected
candidate, the sequences to order and the limitations of the analysis. Output
in Markdown, printable HTML and Word format.

---

## What distinguishes it

Four things, each addressing a specific methodological gap rather than adding
a feature.

### 1. Positional weighting in off-target analysis

BLAST was built to find homology between long sequences. For a 21-nucleotide
query two of its properties work against the user: the default `word_size` of
11 means shorter matches are never seeded, so the **7–8 nucleotide seed region
is invisible by construction**; and the E-value scales with query length, so
the significance filter discards biologically meaningful hits.

BLAST also treats every position as equivalent. In an siRNA they are not.
This pipeline weights mismatches accordingly:

| guide position | weight | reason |
|---|---|---|
| 1 | 0.5 | does not pair with the target; housed in the MID pocket of AGO |
| 2–8 | 3.0 | seed region, determines recognition |
| 9–11 | 2.0 | cleavage site; a mismatch blocks endonucleolytic activity |
| 12+ | 1.0 | supports binding, does not determine recognition |

### 2. Correction for database size

An exact 8-nucleotide match occurs by chance once every 4⁸ = 65 536 positions
— roughly **1800 times in a genome of 1.2 × 10⁸ nucleotides**. A rule of the
form "reject on any seed hit" therefore rejects everything at genome scale.
So does a minimum cost taken over tens of thousands of random positions: the
minimum of many trials is low by construction.

Seed hits are reported as **enrichment relative to chance expectation**.
Rejection rests on a stricter criterion — at most three mismatches with an
intact seed — whose expected chance occurrence in a database of 10⁸
nucleotides is negligible, and which therefore needs no correction for
multiple testing.

### 3. Sensitivity analysis of the scoring weights

Weights are literature-informed but arbitrary. Nobody knows whether asymmetry
should carry 0.25 or 0.30.

The application varies each weight by ±25 % and recomputes the ranking. If
the top five is unchanged, the ordering reflects properties of the sequences
and this can be stated in the methods. If it changes, the ordering reflects
the choice of weights — and the user is told, rather than finding out from a
reviewer.

A Pareto front is computed alongside. A candidate scored highly but lying
outside the front implies another candidate is better in every criterion at
once; its position derives from the weighting, not from the sequence.

### 4. Host-specific rule profiles

Design rules are not universal across kingdoms, and treating them as such is a
common source of quiet error.

Plant AGO1 shows MID-pocket selectivity for 5-prime U, AGO2 for 5-prime A;
since the antiviral pathway runs mainly through these, the plant profile
requires one of them. Mammalian profiles filter the immunostimulatory motifs
recognised by TLR7 and TLR8 — a filter that is meaningless in plants, which
possess neither Toll-like receptors nor an interferon pathway. GC windows
differ. Relevant length classes differ.

Four profiles — plant, mammal, insect, neutral — make these differences
explicit. They are not cosmetic: on the same input, the mammalian profile
rejects candidates the plant profile retains, and vice versa.

---

## Installation

```bash
conda create -n sirna python=3.11 -y
conda activate sirna
pip install -r requirements.txt
```

`ViennaRNA` is optional. Without it the pipeline runs in reduced mode:
structural metrics are skipped, every other metric is unaffected, and a
warning is shown so the omission cannot pass unnoticed.

Detailed instructions for users unfamiliar with the command line are in
[`docs/INSTALL.md`](docs/INSTALL.md).

---

## Quick start

### Graphical interface

```bash
streamlit run app.py
```

Opens in a browser. Paste a sequence, press **Run analysis**.

### Command line

The worked example below uses the coding sequence of green fluorescent
protein, a reporter used across plants, animals and insects, and therefore a
convenient neutral test case for any host profile.

```bash
# design only
python main.py --target gfp.fasta --host plant --out results/

# with off-target analysis against a host transcriptome
python main.py --target gfp.fasta --host plant \
               --transcriptome host_cdna.fasta --out results/

# mammalian target, single length class
python main.py --target target_mrna.fasta --host mammal \
               --lengths 21 --out results/

# build expression cassettes as well
python main.py --target gfp.fasta --host plant \
               --promoter AtU6-1 --cloning goldengate_BsaI --out results/
```

Building the transcriptome index takes a few minutes; it is worth doing once
and reusing:

```bash
python main.py --build-index host_cdna.fasta --index-out host.idx
python main.py --target gfp.fasta --index host.idx --out results/
```

### Validation

```bash
python validate.py --data benchmark.csv --host mammal --out validation/
```

---

## Method

```
target sequence  ──┐
isolate set  ──────┤
host transcriptome ┘
        │
        ├─ 1  conservation
        │      exact k-mer coverage across isolates
        │      threshold from the knee of the coverage curve
        │
        ├─ 2  candidate generation
        │      sliding windows, 21 / 22 / 24 nt
        │
        ├─ 3  physicochemical filters
        │      GC range, homopolymers, 5-prime identity — per host profile
        │
        ├─ 4  thermodynamics
        │      nearest-neighbour model, Xia et al. 1998
        │      asymmetry as the primary metric
        │
        ├─ 5  secondary structure
        │      guide folding and target accessibility, ViennaRNA
        │
        ├─ 6  off-target
        │      streaming scan, positional weights,
        │      correction for database size,
        │      intended-target recognition
        │
        ├─ 7  scoring
        │      min-max normalisation → weights → ranking
        │      sensitivity analysis, Pareto front
        │
        ├─ 8  constructs
        │      promoter and cloning-method requirements
        │
        └─ 9  laboratory forms and report
```

---

## Modules

| module | function |
|---|---|
| `thermo.py` | nearest-neighbour thermodynamics, duplex asymmetry |
| `filters.py` | physicochemical filters, host miRNA collision check |
| `hosts.py` | organism profiles |
| `design.py` | candidate generation, secondary structure |
| `offtarget.py` | streaming off-target scan with positional weights |
| `conservation.py` | k-mer coverage, automated threshold optimisation |
| `scoring.py` | normalisation, ranking, sensitivity analysis, Pareto front |
| `constructs.py` | shRNA and expression cassette assembly |
| `oligos.py` | laboratory forms and specificity controls |
| `ncbi.py` | NCBI E-utilities retrieval with quality filtering |
| `report.py` | in silico synthesis report — Markdown, HTML, DOCX |
| `blast_criteria.py` | BLAST-based scoring with collinearity diagnostics |
| `validate.py` | validation against measured efficacy |
| `main.py` | command-line orchestrator |
| `app.py` | graphical interface |

---

## Tests

```bash
python -m pytest tests/ -v
```

110 tests across three files.

**Thermodynamics** — completeness of the parameter table, reverse-complement
symmetry of all sixteen dinucleotides, agreement with published values, a
hand-calculated duplex, independence of the result from which strand is
supplied. Asymmetry is tested separately: sign convention, antisymmetry under
strand exchange, and a value of zero for a self-complementary sequence.

**Off-target** — detection of planted matches, correct handling of empty
input, and **regression tests for three defects that occurred during
development**: the intended target counted as an off-target hit, a minimum
cost taken without correction for multiple testing, and a full transcriptome
index held in memory. Each now has a test preventing its return.

**Pipeline** — host profiles, filters, scoring (including that weights must
sum to one and that normalisation prevents scale domination), conservation,
controls, oligonucleotide forms, cassette assembly.

---

## Comparison with existing tools

| | this tool | pssRNAit | si-Fi21 |
|---|---|---|---|
| Trained efficacy model | no | yes, SVM | no |
| Reported correlation with measured efficacy | pending | R = 0.709 / 0.686 | not reported |
| Positional weighting of off-target mismatches | yes | no | no |
| Correction for database size | yes | no | no |
| Sensitivity analysis of weights | yes | no | no |
| Host-specific rule profiles | 4 | plant only | plant only |
| Conservation across isolates | yes | no | no |
| Laboratory forms and controls | yes | no | no |
| Reproducibility report | yes | no | no |
| Open source | yes | web only | yes |

The honest reading: **pssRNAit has the better efficacy model**, trained on
2431 measured siRNAs, and that is a real advantage this tool does not yet
match. What this tool offers instead is a documented and inspectable method —
every threshold justified, every weighting tested, every limitation stated —
together with the steps that follow design: conservation, controls, and the
sequences to order.

---

## Limitations

Stated plainly, since they bear on how results should be read.

**No experimental validation of plant predictions.** Every benchmark dataset
of measured siRNA efficacy is mammalian. A correlation obtained on such data
characterises the scoring function on mammalian sequences; it does not
validate the plant profile, whose distinguishing rules — the 5-prime U
preference of AGO1, the absence of TLR-mediated motif filtering, the narrower
GC window — are precisely those the mammalian data cannot test. This is a gap
in the field, not only in this tool.

**Manually set scoring weights.** Weights are literature-informed rather than
learned. The sensitivity analysis reports whether the ranking depends on them,
which is a mitigation rather than a solution.

**Pure-Python off-target scan.** About 0.5 million nucleotides per second; a
50 Mnt transcriptome takes roughly two minutes. Adequate for plant
transcriptomes. An FM-index would be appropriate at human genome scale.

**Approximated target accessibility.** Estimated from the free energy of a
local window rather than from unpaired probability derived from the Boltzmann
ensemble. The more accurate function is implemented but disabled by default
on grounds of speed.

**No ionic-strength correction.** The nearest-neighbour parameters apply to
1 M NaCl. Reported melting temperatures are reference values, not
physiological ones. Ranking is unaffected, since all candidates are treated
identically.

**Conservation detects substitutions only.** Insertions and deletions would
require a multiple sequence alignment. Shannon entropy is implemented for
input that has already been aligned.

**Off-target analysis addresses transcripts, not chromatin.** Effects mediated
by transcriptional silencing or by competition for the silencing machinery
are outside its scope.

---

## References

### Thermodynamics and structure

1. **Xia T., SantaLucia J., Burkard M.E., Kierzek R., Schroeder S.J., Jiao X.,
   Cox C., Turner D.H.** (1998) Thermodynamic parameters for an expanded
   nearest-neighbor model for formation of RNA duplexes with Watson-Crick base
   pairs. *Biochemistry* 37:14719–14735.
2. **Lorenz R., Bernhart S.H., Höner zu Siederdissen C., Tafer H., Flamm C.,
   Stadler P.F., Hofacker I.L.** (2011) ViennaRNA Package 2.0. *Algorithms for
   Molecular Biology* 6:26.
3. **Mathews D.H., Disney M.D., Childs J.L., Schroeder S.J., Zuker M.,
   Turner D.H.** (2004) Incorporating chemical modification constraints into a
   dynamic programming algorithm for prediction of RNA secondary structure.
   *PNAS* 101:7287–7292.

### Mechanism of RNA interference

4. **Fire A., Xu S., Montgomery M.K., Kostas S.A., Driver S.E., Mello C.C.**
   (1998) Potent and specific genetic interference by double-stranded RNA in
   *Caenorhabditis elegans*. *Nature* 391:806–811.
5. **Elbashir S.M., Harborth J., Lendeckel W., Yalcin A., Weber K.,
   Tuschl T.** (2001) Duplexes of 21-nucleotide RNAs mediate RNA interference
   in cultured mammalian cells. *Nature* 411:494–498.
6. **Khvorova A., Reynolds A., Jayasena S.D.** (2003) Functional siRNAs and
   miRNAs exhibit strand bias. *Cell* 115:209–216.
7. **Schwarz D.S., Hutvágner G., Du T., Xu Z., Aronin N., Zamore P.D.** (2003)
   Asymmetry in the assembly of the RNAi enzyme complex. *Cell* 115:199–208.
8. **Bartel D.P.** (2009) MicroRNAs: target recognition and regulatory
   functions. *Cell* 136:215–233.

### Design rules

9. **Reynolds A., Leake D., Boese Q., Scaringe S., Marshall W.S.,
   Khvorova A.** (2004) Rational siRNA design for RNA interference. *Nature
   Biotechnology* 22:326–330.
10. **Ui-Tei K., Naito Y., Takahashi F., Haraguchi T., Ohki-Hamazaki H.,
    Juni A., Ueda R., Saigo K.** (2004) Guidelines for the selection of highly
    effective siRNA sequences for mammalian and chick RNA interference.
    *Nucleic Acids Research* 32:936–948.
11. **Shabalina S.A., Spiridonov A.N., Ogurtsov A.Y.** (2006) Computational
    models with thermodynamic and composition features improve siRNA design.
    *BMC Bioinformatics* 7:65.
12. **Matveeva O., Nechipurenko Y., Rossi L., Moore B., Sætrom P.,
    Ogurtsov A.Y., Atkins J.F., Shabalina S.A.** (2007) Comparison of
    approaches for rational siRNA design leading to a new efficient and
    transparent method. *Nucleic Acids Research* 35:e63.
13. **Huesken D., Lange J., Mickanin C., Weiler J., Asselbergs F., Warner J.,
    Meloon B., Engel S., Rosenberg A., Cohen D., Labow M., Reinhardt M.,
    Natt F., Hall J.** (2005) Design of a genome-wide siRNA library using an
    artificial neural network. *Nature Biotechnology* 23:995–1001.
14. **Naito Y., Yoshimura J., Morishita S., Ui-Tei K.** (2009) siDirect 2.0:
    updated software for designing functional siRNA with reduced seed-dependent
    off-target effect. *BMC Bioinformatics* 10:392.

### Off-target effects and controls

15. **Jackson A.L., Bartz S.R., Schelter J., Kobayashi S.V., Burchard J.,
    Mao M., Li B., Cavet G., Linsley P.S.** (2003) Expression profiling reveals
    off-target gene regulation by RNAi. *Nature Biotechnology* 21:635–637.
16. **Birmingham A., Anderson E.M., Reynolds A., Ilsley-Tyree D., Leake D.,
    Fedorov Y., Baskerville S., Maksimova E., Robinson K., Karpilow J.,
    Marshall W.S., Khvorova A.** (2006) 3′ UTR seed matches, but not overall
    identity, are associated with RNAi off-targets. *Nature Methods*
    3:199–204.
17. **Buehler E., Chen Y.-C., Martin S.** (2012) C911: a bench-level control
    for sequence specific siRNA off-target effects. *PLoS ONE* 7(12):e51942.
18. **Judge A.D., Sood V., Shaw J.R., Fang D., McClintock K., MacLachlan I.**
    (2005) Sequence-dependent stimulation of the mammalian innate immune
    response by synthetic siRNA. *Nature Biotechnology* 23:457–462.
19. **Hornung V., Guenthner-Biller M., Bourquin C., Ablasser A., Schlee M.,
    Uematsu S., Noronha A., Manoharan M., Akira S., de Fougerolles A.,
    Endres S., Hartmann G.** (2005) Sequence-specific potent induction of
    IFN-α by short interfering RNA in plasmacytoid dendritic cells through
    TLR7. *Nature Medicine* 11:263–270.
20. **Altschul S.F., Gish W., Miller W., Myers E.W., Lipman D.J.** (1990)
    Basic local alignment search tool. *Journal of Molecular Biology*
    215:403–410.

### Plant RNA interference

21. **Mi S., Cai T., Hu Y., Chen Y., Hodges E., Ni F., Wu L., Li S., Zhou H.,
    Long C., Chen S., Hannon G.J., Qi Y.** (2008) Sorting of small RNAs into
    *Arabidopsis* argonaute complexes is directed by the 5′ terminal
    nucleotide. *Cell* 133:116–127.
22. **Takeda A., Iwasaki S., Watanabe T., Utsumi M., Watanabe Y.** (2008) The
    mechanism selecting the guide strand from small RNA duplexes is different
    among Argonaute proteins. *Plant and Cell Physiology* 49:493–500.
23. **Fahlgren N., Carrington J.C.** (2010) miRNA target prediction in plants.
    *Methods in Molecular Biology* 592:51–57.
24. **Carbonell A., López C., Daròs J.-A.** (2019) Fast-forward identification
    of highly effective artificial small RNAs against different tomato
    spotted wilt virus isolates. *Molecular Plant-Microbe Interactions*
    32:142–156.

### Antiviral strategy and multiplexing

25. **ter Brake O., Konstantinova P., Ceylan M., Berkhout B.** (2006)
    Silencing of HIV-1 with RNA interference: a multiple shRNA approach.
    *Molecular Therapy* 14:883–892.
26. **Taxman D.J., Livingstone L.R., Zhang J., Conti B.J., Iocca H.A.,
    Williams K.L., Lich J.D., Ting J.P.-Y., Reed W.** (2006) Criteria for
    effective design, construction, and gene knockdown by shRNA vectors.
    *BMC Biotechnology* 6:7.

### Methods

27. **Varkonyi-Gasic E., Wu R., Wood M., Walton E.F., Hellens R.P.** (2007)
    Protocol: a highly sensitive RT-PCR method for detection and
    quantification of microRNAs. *Plant Methods* 3:12.
28. **Satopää V., Albrecht J., Irwin D., Raghavan B.** (2011) Finding a
    "kneedle" in a haystack: detecting knee points in system behavior.
    *31st International Conference on Distributed Computing Systems
    Workshops*, pp. 166–171.

*Bibliographic details should be verified against the sources before citation.*

---

## Citation

```bibtex
@software{jarecka_sirna_design_studio_2026,
  author  = {Jarecka, Antonina},
  title   = {siRNA Design Studio: transparent design and evaluation
             of small interfering RNA},
  version = {1.0.0},
  year    = {2026},
  url     = {https://github.com/Tosia-J/sirna-design-studio}
}
```

Once a DOI is assigned, cite the **version DOI** rather than the concept DOI,
so that readers obtain exactly the code that was used.

---

## Contributing

Issues and pull requests are welcome. Particularly useful contributions would
be: a plant dataset of measured siRNA efficacy, which does not currently
exist at usable scale; an FM-index implementation of the off-target scan; and
additional host profiles with the supporting literature.

---

## Licence

MIT — see [LICENSE](LICENSE).
