# siRNA Design Studio

A pipeline for designing and evaluating small interfering RNA, with
host-specific rule profiles, position-weighted off-target analysis and
alignment-free conservation assessment across pathogen isolates.

<!-- Add once available:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![Tests](https://img.shields.io/badge/tests-110%20passing-brightgreen)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
-->

<!-- Add a screenshot:
![Application interface](docs/screenshot.png)
-->

---

## What this does differently

Several good siRNA design tools already exist. This one differs in four
respects, each addressing a specific methodological gap.

**Position-weighted off-target analysis.** BLAST treats every position as
equivalent and, with its default `word_size` of 11, does not seed matches to
the 7–8 nt seed region at all — so the most biologically important class of
off-target site is invisible by construction. This pipeline scores mismatches
with weights reflecting their actual contribution to target recognition:
3.0 for the seed (positions 2–8), 2.0 for the cleavage site (9–11), 1.0 for
the 3′ region.

**Correction for database size.** An exact 8-nucleotide match occurs by
chance roughly once every 65 536 positions, that is about 1800 times in a
genome of 1.2 × 10⁸ nucleotides. Rejecting a candidate on any seed hit, or on
a minimum cost taken over tens of thousands of random positions, rejects
everything at genome scale. Seed hits are therefore reported as enrichment
relative to chance expectation, and rejection rests on a criterion whose
expected chance occurrence is negligible.

**Explicit sensitivity analysis.** Scoring weights are literature-informed
but arbitrary. Each weight is varied by ±25 % and the ranking recomputed; the
stability of the top five is reported alongside the result. A ranking that
changes under this perturbation reflects the choice of weights rather than
properties of the sequence, and the user is told so.

**Host-specific profiles.** Design rules are not universal across kingdoms.
Plant AGO1 shows MID-pocket selectivity for 5′-U (Mi et al. 2008); plants
possess neither Toll-like receptors nor an interferon pathway, so the
mammalian immunostimulatory motif filter does not apply to them. Profiles for
plant, mammal, insect and a neutral default make these differences explicit
rather than implicit.

---

## Installation

```bash
conda create -n sirna python=3.11 -y
conda activate sirna
pip install -r requirements.txt
```

`ViennaRNA` is optional. Without it the pipeline runs in reduced mode:
structural metrics are skipped, everything else is unaffected, and a warning
is displayed.

---

## Quick start

**Graphical interface**

```bash
streamlit run app.py
```

**Command line**

```bash
python main.py --target gene.fasta --host plant --out results/
```

With off-target analysis:

```bash
python main.py --target gene.fasta --host plant \
               --transcriptome TAIR10_cdna.fasta --out results/
```

**Validation against measured efficacy**

```bash
python validate.py --data benchmark.csv --host mammal --out validation/
```

---

## Pipeline

```
target sequence
      │
      ├─ conservation      exact k-mer coverage across isolates,
      │                    threshold determined from the data
      │
      ├─ candidate generation      21 / 22 / 24 nt windows
      │
      ├─ physicochemical filters   host profile
      │
      ├─ thermodynamics            nearest-neighbour, Xia et al. 1998
      │                            asymmetry as the primary metric
      │
      ├─ secondary structure       ViennaRNA
      │
      ├─ off-target                streaming scan, positional weights
      │
      ├─ scoring                   normalisation → weights → ranking
      │                            + sensitivity analysis + Pareto front
      │
      ├─ constructs                shRNA, expression cassettes
      │
      └─ order sheet               duplex, hairpin, IVT templates,
                                   detection primers, C911 controls
```

---

## Modules

| module | function |
|---|---|
| `thermo.py` | nearest-neighbour thermodynamics, asymmetry |
| `filters.py` | physicochemical filters, host miRNA collision check |
| `hosts.py` | organism profiles |
| `design.py` | candidate generation, secondary structure |
| `offtarget.py` | streaming off-target scan with positional weights |
| `conservation.py` | k-mer coverage, automated threshold optimisation |
| `scoring.py` | normalisation, ranking, sensitivity analysis, Pareto front |
| `constructs.py` | shRNA and expression cassette assembly |
| `oligos.py` | laboratory forms and specificity controls |
| `ncbi.py` | NCBI E-utilities retrieval with quality filtering |
| `report.py` | in silico synthesis report (Markdown, HTML, DOCX) |
| `blast_criteria.py` | legacy BLAST scoring with collinearity diagnostics |
| `validate.py` | validation against measured efficacy |

---

## Tests

```bash
python -m pytest tests/ -v
```

110 tests. Coverage includes:

- integrity of the thermodynamic parameter table, including reverse-complement
  symmetry of all sixteen dinucleotides
- hand-calculated free energies checked against the published table
- antisymmetry of the asymmetry metric under strand exchange
- regression tests for three defects that occurred during development:
  the intended target counted as off-target, minimum cost without correction
  for multiple testing, and a full transcriptome index held in memory

---

## Limitations

Stated explicitly, since they bear on how results should be interpreted.

**No experimental validation of plant predictions.** All available benchmark
datasets of measured siRNA efficacy are mammalian. A correlation obtained on
them characterises the scoring function on mammalian data; it does not
validate the plant profile, whose distinguishing rules those data cannot test.
This is a limitation of the field rather than of the tool.

**Manually set scoring weights.** Weights are literature-informed rather than
learned. The sensitivity analysis reports whether the ranking depends on them,
but a model trained on measured efficacy data would be preferable.

**Pure-Python off-target scan.** Approximately 0.5 million nucleotides per
second; a transcriptome of 50 Mnt takes about two minutes. Adequate for plant
transcriptomes; an FM-index would be appropriate at human genome scale.

**Approximated target accessibility.** Estimated from window free energy
rather than from unpaired probability derived from the Boltzmann ensemble.
The more accurate function is implemented but disabled by default.

**No ionic-strength correction.** The nearest-neighbour parameters apply to
1 M NaCl. Reported melting temperatures are not physiological values.

**Conservation detects substitutions only.** Insertions and deletions would
require a multiple sequence alignment. Shannon entropy is implemented for
pre-aligned input.

---

## References

Core methodology:

1. Xia T. et al. (1998) *Biochemistry* 37:14719–14735 — nearest-neighbour
   thermodynamic parameters
2. Lorenz R. et al. (2011) *Algorithms Mol Biol* 6:26 — ViennaRNA
3. Khvorova A. et al. (2003) *Cell* 115:209–216 — strand bias
4. Schwarz D.S. et al. (2003) *Cell* 115:199–208 — asymmetry in RISC assembly
5. Bartel D.P. (2009) *Cell* 136:215–233 — seed site classification
6. Mi S. et al. (2008) *Cell* 133:116–127 — 5′ nucleotide directs AGO sorting
   in *Arabidopsis*
7. Buehler E. et al. (2012) *PLoS ONE* 7(12):e51942 — C911 specificity control

A complete list is given in the Methods tab of the application.

---

## Citation

```bibtex
@software{jarecka_sirna_design_studio,
  author  = {Jarecka, Antonina},
  title   = {siRNA Design Studio},
  version = {1.0.0},
  year    = {2026},
  url     = {https://github.com/YOUR-USERNAME/sirna-design-studio}
}
```

Once a DOI is assigned, cite the version DOI rather than the concept DOI, so
that readers obtain exactly the code that was used.

---

## Licence

MIT — see [LICENSE](LICENSE).
