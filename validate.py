#!/usr/bin/env python3
"""
=============================================================================
validate.py  --  VALIDATION AGAINST MEASURED SILENCING EFFICACY
=============================================================================

WHY THIS EXISTS

Every scoring scheme produces a ranking. The question that matters is whether
that ranking corresponds to anything measurable. Without a correlation against
experimentally determined efficacy, a scoring function is an opinion expressed
in numbers.

Published tools report this figure and it is the natural basis for comparison:

    pssRNAit (SVM)          R = 0.709 on the training set TR2431
                            R = 0.686 on the independent set TE249

A comparable number is what this script produces.

=============================================================================
BENCHMARK DATASETS
=============================================================================

    Huesken et al. 2005     2431 siRNAs, H1299 cells, branched-DNA readout.
                            The field standard. Nature Biotechnology 23:995.

    Reynolds et al. 2004    248 siRNAs. Nature Biotechnology 22:326.

    Vickers et al. 2003     80 siRNAs.

    Harborth et al. 2003    44 siRNAs.

    Ui-Tei et al. 2004      62 siRNAs. Nucleic Acids Research 32:936.

Composite sets combining several of these are also in circulation and are
preferable: training on Huesken alone and testing elsewhere gives markedly
weaker performance, which is attributed to heterogeneity of experimental
conditions between studies.

=============================================================================
A CAVEAT THAT MUST BE STATED, NOT BURIED
=============================================================================

Every one of these datasets is MAMMALIAN. There is no plant dataset of
comparable size. Consequences:

  - a correlation obtained here characterises the scoring function on
    mammalian data and licenses comparison with mammalian tools;
  - it does NOT validate the plant profile, whose distinguishing rules
    (5'-U preference of AGO1, absence of TLR-mediated motif filtering,
    narrower GC window) are precisely the ones the mammalian data cannot
    test;
  - reporting a mammalian correlation as though it validated plant
    predictions would be a misrepresentation.

The honest formulation is: the scoring function achieves R = x on mammalian
benchmark data; the plant profile remains unvalidated for want of a suitable
dataset, and this is a limitation of the field rather than of the tool.

=============================================================================
INPUT FORMAT
=============================================================================

A CSV or TSV with at least two columns:

    sequence   the sense (passenger) strand, or the guide strand if the
               --guide flag is used, 19-21 nt
    efficacy   measured knockdown; any monotonic scale, since Spearman
               correlation is used

Optional:

    target     accession or identifier of the target mRNA, used for
               accessibility if --targets is supplied

The Huesken data are distributed as supplementary material to the original
paper and through several later compilations; obtain them from the source
rather than from a redistribution of unknown provenance.

=============================================================================
USAGE
=============================================================================

    python validate.py --data huesken.csv --host mammal --out validation/
    python validate.py --data huesken.csv --guide --bootstrap 1000

Author: Antonina Jarecka
=============================================================================
"""

import argparse
import csv
import json
import math
import os
import random
import sys
from typing import Dict, List, Optional, Tuple

import thermo
import filters
import hosts
import scoring

try:
    import design
    VIENNA = design.VIENNA_DOSTEPNE
except ImportError:
    design = None
    VIENNA = False


# ============================================================================
# STATISTICS
# ============================================================================

def _rangi(wartosci: List[float]) -> List[float]:
    """Ranks with ties averaged, as required by Spearman's coefficient."""
    n = len(wartosci)
    idx = sorted(range(n), key=lambda i: wartosci[i])
    rangi = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and wartosci[idx[j + 1]] == wartosci[idx[i]]:
            j += 1
        srednia = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            rangi[idx[k]] = srednia
        i = j + 1
    return rangi


def pearson(x: List[float], y: List[float]) -> float:
    n = len(x)
    if n < 3:
        return 0.0
    sx, sy = sum(x) / n, sum(y) / n
    licz = sum((a - sx) * (b - sy) for a, b in zip(x, y))
    m1 = math.sqrt(sum((a - sx) ** 2 for a in x))
    m2 = math.sqrt(sum((b - sy) ** 2 for b in y))
    return licz / (m1 * m2) if m1 * m2 > 1e-12 else 0.0


def spearman(x: List[float], y: List[float]) -> float:
    """Rank correlation - insensitive to the scale on which efficacy is
    reported, which differs between studies."""
    return pearson(_rangi(x), _rangi(y))


def bootstrap_ci(x: List[float], y: List[float], n_iter: int = 1000,
                 alpha: float = 0.05, seed: int = 0
                 ) -> Tuple[float, float]:
    """Percentile bootstrap confidence interval for Spearman's rho."""
    rng = random.Random(seed)
    n = len(x)
    proby = []
    for _ in range(n_iter):
        idx = [rng.randrange(n) for _ in range(n)]
        proby.append(spearman([x[i] for i in idx], [y[i] for i in idx]))
    proby.sort()
    lo = proby[int(alpha / 2 * n_iter)]
    hi = proby[int((1 - alpha / 2) * n_iter) - 1]
    return lo, hi


def auc_top_bottom(wyniki: List[float], skutecznosc: List[float],
                   prog_frakcja: float = 0.25) -> float:
    """
    Area under the ROC curve for the practical question: does the score
    separate the most effective quartile from the least effective one?

    This is closer to how the tool is used than a correlation over the whole
    range - in practice one orders the top few candidates.
    """
    n = len(skutecznosc)
    posort = sorted(skutecznosc)
    prog_gora = posort[int((1 - prog_frakcja) * n)]
    prog_dol = posort[int(prog_frakcja * n)]

    dodatnie = [w for w, s in zip(wyniki, skutecznosc) if s >= prog_gora]
    ujemne = [w for w, s in zip(wyniki, skutecznosc) if s <= prog_dol]
    if not dodatnie or not ujemne:
        return 0.5

    lepsze = sum(1 for a in dodatnie for b in ujemne if a > b)
    remisy = sum(1 for a in dodatnie for b in ujemne if a == b)
    return (lepsze + 0.5 * remisy) / (len(dodatnie) * len(ujemne))


# ============================================================================
# SCORING OF BENCHMARK SEQUENCES
# ============================================================================

def ocen_sekwencje(sekwencje: List[str], profil: hosts.HostProfile,
                   jako_guide: bool = False,
                   verbose: bool = True) -> List[Dict]:
    """
    Computes the pipeline metrics for each benchmark sequence.

    Filters are NOT applied. Validation asks whether the score correlates
    with efficacy across the whole range; discarding low-scoring sequences
    first would remove exactly the cases that test the score.
    """
    out = []
    for i, s in enumerate(sekwencje):
        s = s.upper().replace('U', 'T')
        if set(s) - set('ACGT') or len(s) < 19:
            out.append(None)
            continue

        guide = s if jako_guide else thermo.revcomp_rna(s)
        guide = thermo.to_rna(guide)
        sense = thermo.to_rna(thermo.revcomp_rna(guide))

        dg_g, dg_p, asym = thermo.asymmetry(guide)
        g_pct = thermo.gc_content(sense)

        if VIENNA:
            _, mfe_g = design.fold_guide(guide)
        else:
            mfe_g = 0.0

        out.append({
            'nazwa': f'seq{i}',
            'guide_rna': guide,
            'asymetria': asym,
            'dostepnosc_celu': 0.0,     # no target context in benchmark data
            'mfe_guide': mfe_g,
            'gc_optymalnosc': 1.0 - abs(g_pct - 42.0) / 42.0,
            'koszt_offtarget': 20.0,    # no transcriptome in benchmark data
            'seed_czystosc': 0,
            'gc_proc': g_pct,
            'guide_5': guide[0],
        })
        if verbose and (i + 1) % 500 == 0:
            print(f'  scored {i + 1} / {len(sekwencje)}')
    return out


# ============================================================================
# INDIVIDUAL METRICS
# ============================================================================

def korelacje_pojedynczych(kandydaci: List[Dict],
                           skutecznosc: List[float]) -> Dict[str, Dict]:
    """
    Correlation of each metric taken alone.

    This matters more than the composite score. If a metric correlates with
    efficacy near zero, its weight in the composite is not justified by the
    data, whatever the literature says. If one metric carries almost all the
    signal, the others are decoration.
    """
    wynik = {}
    for m in ('asymetria', 'gc_optymalnosc', 'mfe_guide', 'gc_proc'):
        wart = [c[m] for c in kandydaci]
        wynik[m] = {
            'spearman': round(spearman(wart, skutecznosc), 4),
            'pearson': round(pearson(wart, skutecznosc), 4),
        }

    # 5'-terminal nucleotide, tested as a categorical variable
    wg_nt = {}
    for c, s in zip(kandydaci, skutecznosc):
        wg_nt.setdefault(c['guide_5'], []).append(s)
    wynik['guide_5_nucleotide'] = {
        nt: {'n': len(v), 'mean_efficacy': round(sum(v) / len(v), 4)}
        for nt, v in sorted(wg_nt.items()) if len(v) >= 10}

    return wynik


# ============================================================================
# MAIN
# ============================================================================

def wczytaj_dane(sciezka: str, kol_seq: str, kol_eff: str
                 ) -> Tuple[List[str], List[float]]:
    sep = '\t' if sciezka.endswith(('.tsv', '.txt')) else ','
    sekw, sk = [], []
    with open(sciezka, newline='') as fh:
        czyt = csv.DictReader(fh, delimiter=sep)
        if kol_seq not in czyt.fieldnames or kol_eff not in czyt.fieldnames:
            raise ValueError(
                f'Columns "{kol_seq}" and "{kol_eff}" required. '
                f'File contains: {czyt.fieldnames}')
        for w in czyt:
            try:
                e = float(w[kol_eff])
            except (TypeError, ValueError):
                continue
            sekw.append(w[kol_seq].strip())
            sk.append(e)
    return sekw, sk


def main() -> int:
    p = argparse.ArgumentParser(
        description='Validate the scoring function against measured efficacy',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    p.add_argument('--data', required=True, help='CSV/TSV benchmark file')
    p.add_argument('--seq-col', default='sequence')
    p.add_argument('--eff-col', default='efficacy')
    p.add_argument('--guide', action='store_true',
                   help='sequences are guide strands, not sense strands')
    p.add_argument('--host', default='mammal', choices=sorted(hosts.PROFILE))
    p.add_argument('--bootstrap', type=int, default=1000)
    p.add_argument('--out', default='validation')
    a = p.parse_args()

    os.makedirs(a.out, exist_ok=True)
    profil = hosts.get_profile(a.host)

    print('=' * 74)
    print('VALIDATION AGAINST MEASURED SILENCING EFFICACY')
    print('=' * 74)

    sekw, sk = wczytaj_dane(a.data, a.seq_col, a.eff_col)
    print(f'Loaded {len(sekw)} sequences from {a.data}')
    print(f'Host profile: {profil.nazwa}')
    if not VIENNA:
        print('WARNING: ViennaRNA unavailable - structural metrics excluded '
              'from the composite score.')

    kand = ocen_sekwencje(sekw, profil, jako_guide=a.guide)
    pary = [(c, e) for c, e in zip(kand, sk) if c is not None]
    kand = [c for c, _ in pary]
    sk = [e for _, e in pary]
    print(f'Scored {len(kand)} sequences '
          f'({len(sekw) - len(kand)} rejected as malformed)')

    if len(kand) < 30:
        print('ERROR: too few usable sequences for a meaningful correlation.',
              file=sys.stderr)
        return 1

    ranking = scoring.normalize_and_score(kand)
    wg_nazwy = {c['nazwa']: c for c in ranking}
    wyniki = [wg_nazwy[c['nazwa']]['wynik'] for c in kand]

    rho = spearman(wyniki, sk)
    r = pearson(wyniki, sk)
    lo, hi = bootstrap_ci(wyniki, sk, a.bootstrap)
    auc = auc_top_bottom(wyniki, sk)
    poj = korelacje_pojedynczych(kand, sk)

    print()
    print('COMPOSITE SCORE')
    print('-' * 74)
    print(f'  Spearman rho     {rho:+.4f}   '
          f'95% CI [{lo:+.4f}, {hi:+.4f}]  ({a.bootstrap} bootstrap samples)')
    print(f'  Pearson r        {r:+.4f}')
    print(f'  AUC top vs bottom quartile   {auc:.4f}')
    print()
    print('  For comparison, pssRNAit reports R = 0.709 on its training set '
          'and\n  R = 0.686 on an independent set.')

    print()
    print('INDIVIDUAL METRICS')
    print('-' * 74)
    for m, v in poj.items():
        if m == 'guide_5_nucleotide':
            continue
        print(f'  {m:<20} rho {v["spearman"]:+.4f}   r {v["pearson"]:+.4f}')
    print()
    print("  5'-terminal nucleotide of the guide strand:")
    for nt, v in poj['guide_5_nucleotide'].items():
        print(f'    {nt}  n = {v["n"]:>5}   mean efficacy '
              f'{v["mean_efficacy"]:+.4f}')

    print()
    print('INTERPRETATION')
    print('-' * 74)
    uwagi = []
    if abs(rho) < 0.15:
        uwagi.append('The composite score shows essentially no relationship '
                     'with measured efficacy on this dataset. The weighting '
                     'is not supported by the data.')
    elif abs(rho) < 0.35:
        uwagi.append('A weak but detectable relationship. Below what '
                     'published tools report; a trained model would be '
                     'expected to do better.')
    elif abs(rho) < 0.60:
        uwagi.append('A moderate relationship, in the range typical of '
                     'rule-based schemes without machine learning.')
    else:
        uwagi.append('A strong relationship, comparable with published '
                     'machine-learning tools.')

    najlepsza = max((m for m in poj if m != 'guide_5_nucleotide'),
                    key=lambda m: abs(poj[m]['spearman']))
    if abs(poj[najlepsza]['spearman']) > abs(rho):
        uwagi.append(
            f'The single metric "{najlepsza}" correlates more strongly '
            f'(rho {poj[najlepsza]["spearman"]:+.4f}) than the composite '
            f'score. This indicates the weighting is diluting rather than '
            f'combining the signal and should be revisited.')

    uwagi.append('All benchmark datasets are mammalian. This result '
                 'characterises the scoring function on mammalian data and '
                 'does not validate the plant profile, whose distinguishing '
                 'rules these data cannot test.')
    if not VIENNA:
        uwagi.append('ViennaRNA was unavailable, so structural metrics were '
                     'excluded. The result underestimates the full pipeline.')

    for u in uwagi:
        print(f'  - {u}')

    raport = {
        'dataset': a.data, 'n_sequences': len(kand), 'host_profile': a.host,
        'vienna_available': VIENNA,
        'composite': {'spearman': round(rho, 4), 'pearson': round(r, 4),
                      'ci_low': round(lo, 4), 'ci_high': round(hi, 4),
                      'auc_quartiles': round(auc, 4)},
        'individual_metrics': poj,
        'interpretation': uwagi,
        'reference_comparison': {'pssRNAit_training': 0.709,
                                 'pssRNAit_independent': 0.686},
    }
    sciezka = os.path.join(a.out, 'validation_report.json')
    with open(sciezka, 'w') as fh:
        json.dump(raport, fh, indent=2)

    tsv = os.path.join(a.out, 'scored_sequences.tsv')
    with open(tsv, 'w') as fh:
        fh.write('sequence\tguide\tscore\tmeasured_efficacy\tasymmetry\tgc\n')
        for c, w, e in zip(kand, wyniki, sk):
            fh.write(f'{thermo.revcomp_rna(c["guide_rna"])}\t'
                     f'{c["guide_rna"]}\t{w:.4f}\t{e}\t'
                     f'{c["asymetria"]:.3f}\t{c["gc_proc"]:.1f}\n')

    print()
    print(f'Written: {sciezka}\n         {tsv}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
