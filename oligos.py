"""
=============================================================================
oligos.py  --  LABORATORY-READY OLIGONUCLEOTIDE FORMS
=============================================================================

PURPOSE

A ranked guide-strand sequence is not something that can be ordered from a
supplier and used at the bench. Depending on the experimental route, several
distinct physical constructs are required, each with its own end chemistry,
overhangs and flanking sequence.

This module converts a selected candidate into every form needed downstream:

  SYNTHETIC DUPLEX ROUTE
    - guide strand with 2-nt 3' overhang
    - passenger strand with 2-nt 3' overhang
    - annealed duplex representation

  VECTOR-BASED ROUTE
    - shRNA hairpin
    - annealing oligos (top / bottom) with cloning overhangs

  IN VITRO TRANSCRIPTION ROUTE
    - DNA template with T7 promoter, sense and antisense

  DETECTION
    - stem-loop RT primer (Varkonyi-Gasic et al. 2007)
    - forward qPCR primer
    - Northern blot probe

  CONTROLS
    - C911 (positions 9-11 replaced by their complement)
    - seed-mutant
    - scrambled

=============================================================================
NOTE ON CONTROL CHOICE
=============================================================================

The scrambled control is the one most frequently used and the weakest.
Buehler et al. (2012) compared three control designs against a test set of
20 highly active siRNAs comprising 10 true and 10 false positives. Scrambled
controls lost activity in both groups, making them unable to distinguish
on-target from off-target effects. The C911 design - in which bases 9 to 11
are replaced by their complement - separated the two groups completely.

Reasoning: C911 preserves the seed region (positions 2-8) and therefore all
seed-driven off-target activity, while the central mismatches abolish
cleavage of the intended target. A difference between the parent siRNA and
its C911 control is therefore attributable to on-target activity.

This module generates both. C911 is recommended as the primary specificity
control; scrambled is provided because it remains conventional in the
literature.

Author: Antonina Jarecka
=============================================================================
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import random

COMP_DNA = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
COMP_RNA = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

T7_PROMOTER = 'TAATACGACTCACTATAGGG'

# Universal stem-loop scaffold (Varkonyi-Gasic et al. 2007). The six 3'
# nucleotides are replaced by the reverse complement of the small RNA's
# three terminal bases.
STEM_LOOP_SCAFFOLD = 'GTCGTATCCAGTGCAGGGTCCGAGGTATTCGCACTGGATACGAC'
UNIVERSAL_REVERSE = 'GTGCAGGGTCCGAGGT'


def rc_dna(s: str) -> str:
    return ''.join(COMP_DNA[c] for c in reversed(s.upper().replace('U', 'T')))


def rc_rna(s: str) -> str:
    return ''.join(COMP_RNA[c] for c in reversed(s.upper().replace('T', 'U')))


def to_dna(s: str) -> str:
    return s.upper().replace('U', 'T')


def to_rna(s: str) -> str:
    return s.upper().replace('T', 'U')


def gc(s: str) -> float:
    s = s.upper()
    return 100.0 * (s.count('G') + s.count('C')) / len(s) if s else 0.0


# ============================================================================
# CONTROLS
# ============================================================================

def c911_control(guide: str) -> str:
    """
    C911 control: bases 9, 10 and 11 replaced by their complement.

    Preserves the seed region and hence seed-mediated off-target activity,
    while central mismatches abolish target cleavage.

    Reference: Buehler et al. 2012, PLoS ONE 7(12):e51942.
    """
    g = list(to_rna(guide))
    if len(g) < 11:
        raise ValueError('Guide shorter than 11 nt - C911 not applicable.')
    for i in (8, 9, 10):          # 0-based positions 9, 10, 11
        g[i] = COMP_RNA[g[i]]
    return ''.join(g)


def c10_control(guide: str) -> str:
    """C10 control: base 10 only replaced by its complement (weaker variant)."""
    g = list(to_rna(guide))
    g[9] = COMP_RNA[g[9]]
    return ''.join(g)


def seed_mutant_control(guide: str) -> str:
    """
    Seed-mutant control: positions 2-8 replaced by their complement.

    Mirror image of C911 - abolishes seed-mediated off-target activity while
    leaving the 3' region intact. Used together with C911, it allows the two
    contributions to be separated.
    """
    g = list(to_rna(guide))
    for i in range(1, 8):
        g[i] = COMP_RNA[g[i]]
    return ''.join(g)


def scrambled_control(guide: str, target_mrna: Optional[str] = None,
                      seed: int = 42, max_attempts: int = 2000
                      ) -> Optional[str]:
    """
    Scrambled control: same nucleotide composition, shuffled order.

    Constraints applied:
      - no run of 4 identical nucleotides
      - GC content within 5 percentage points of the parent
      - no stretch of 8 nt complementary to the target mRNA, if supplied

    Returns None if no sequence satisfying the constraints is found.

    Note: see module docstring - this control is weaker than C911.
    """
    rng = random.Random(seed)
    g = to_rna(guide)
    bases = list(g)
    target = to_dna(target_mrna) if target_mrna else None

    for _ in range(max_attempts):
        rng.shuffle(bases)
        cand = ''.join(bases)
        if cand == g:
            continue
        if any(b * 4 in cand for b in 'ACGU'):
            continue
        if abs(gc(cand) - gc(g)) > 5.0:
            continue
        if target:
            d = to_dna(cand)
            if any(rc_dna(d[i:i + 8]) in target for i in range(len(d) - 7)):
                continue
        return cand
    return None


# ============================================================================
# DUPLEX FORMS
# ============================================================================

def synthetic_duplex(guide: str, overhang: str = 'dTdT') -> Dict[str, str]:
    """
    Synthetic siRNA duplex with 2-nt 3' overhangs.

    overhang:
      'dTdT'  - deoxythymidine, the conventional choice; cheaper and more
                nuclease-resistant (Elbashir et al. 2001)
      'UU'    - ribonucleotide
      'native'- retains the nucleotides present in the target sequence
    """
    g_core = to_rna(guide)
    p_core = rc_rna(g_core)

    if overhang == 'dTdT':
        oh = 'dTdT'
    elif overhang == 'UU':
        oh = 'UU'
    else:
        oh = ''

    return {
        'guide_5_3': f'{g_core}{oh}',
        'passenger_5_3': f'{p_core}{oh}',
        'guide_core': g_core,
        'passenger_core': p_core,
        'overhang': oh or 'none',
        'duplex_length': len(g_core),
    }


def duplex_diagram(guide: str, overhang: str = 'dTdT') -> str:
    """Text representation of the annealed duplex."""
    d = synthetic_duplex(guide, overhang)
    g, p = d['guide_core'], d['passenger_core']
    oh = d['overhang'] if d['overhang'] != 'none' else ''
    pad = ' ' * len(oh)
    return (f"passenger 5'-{pad}{p}{oh}-3'\n"
            f"              {'|' * len(g)}\n"
            f"guide     3'-{oh[::-1]}{g[::-1]}{pad}-5'")


# ============================================================================
# HAIRPIN AND CLONING OLIGOS
# ============================================================================

LOOPS = {
    'standard_9':  ('TTCAAGAGA',    'Classic 9-nt shRNA loop'),
    'standard_12': ('GTTCAAGAGAAC', '12-nt variant, self-complementary stem'),
    'mir30':       ('TAGTGAAGCCACAGATGTA', 'miR-30 derived, favours Drosha'),
    'low_struct':  ('CTGCTGCTGCT',  'Low intrinsic structure'),
}


def shrna(guide: str, loop_key: str = 'standard_12',
          guide_first: bool = True, terminator: str = '') -> Dict[str, str]:
    """
    shRNA hairpin as a DNA sequence.

    guide_first: places the guide strand at the 5' end, which favours its
    loading into the effector complex.

    terminator: for Pol III promoters, a poly-T tract appended to the insert.
    """
    g = to_dna(guide)
    p = rc_dna(g)
    loop, _ = LOOPS[loop_key]
    core = (g + loop + p) if guide_first else (p + loop + g)
    return {'shrna_dna': core + terminator,
            'loop': loop,
            'length': len(core + terminator),
            'arrangement': 'guide-loop-passenger' if guide_first
                           else 'passenger-loop-guide'}


def annealing_oligos(shrna_dna: str, overhang_5: str = 'CACC',
                     overhang_3: str = 'AAAC') -> Dict[str, str]:
    """
    Top and bottom oligonucleotides for annealing and ligation into a
    digested vector.

    The two oligos are annealed and ligated directly, without PCR. Default
    overhangs correspond to BsaI-digested vectors of the pENTR / pYLsgRNA
    type; replace them with the overhangs generated by your own vector.
    """
    top = overhang_5 + shrna_dna
    bottom = overhang_3 + rc_dna(shrna_dna)
    return {'top_oligo_5_3': top,
            'bottom_oligo_5_3': bottom,
            'top_length': len(top),
            'bottom_length': len(bottom),
            'note': 'Anneal at equimolar concentration: 95 C for 5 min, '
                    'then cool to room temperature over 45-60 min.'}


# ============================================================================
# IN VITRO TRANSCRIPTION
# ============================================================================

def ivt_template(sequence: str, both_strands: bool = True) -> Dict[str, str]:
    """
    DNA template for T7 in vitro transcription.

    T7 RNA polymerase requires the promoter followed by a G at the
    transcription start site; the promoter used here already supplies GGG.

    For dsRNA, both strands are transcribed separately and then annealed.
    """
    seq = to_dna(sequence)
    out = {'t7_promoter': T7_PROMOTER,
           'sense_template_5_3': T7_PROMOTER + seq,
           'note': 'Transcribe, treat with DNase, then purify.'}
    if both_strands:
        out['antisense_template_5_3'] = T7_PROMOTER + rc_dna(seq)
        out['note'] += (' For dsRNA: transcribe both strands separately, '
                        'mix in equimolar amounts, heat to 95 C and cool '
                        'slowly.')
    return out


# ============================================================================
# DETECTION
# ============================================================================

def stem_loop_rt_primer(small_rna: str) -> Dict[str, str]:
    """
    Stem-loop RT primer for detection of a small RNA by qPCR.

    Method: Varkonyi-Gasic et al. 2007, Plant Methods 3:12.

    The scaffold is constant; the six 3' nucleotides are the reverse
    complement of the six terminal nucleotides of the small RNA. The forward
    primer is the 5' portion of the small RNA, extended if necessary to reach
    an acceptable melting temperature.
    """
    s = to_dna(small_rna)
    rt = STEM_LOOP_SCAFFOLD + rc_dna(s[-6:])
    fwd_core = s[:min(len(s) - 6, 15)]
    # GC clamp added when the 5' portion is AT-rich
    fwd = ('GCGCG' + fwd_core) if gc(fwd_core) < 40 else ('GCG' + fwd_core)
    return {'rt_primer_5_3': rt,
            'forward_primer_5_3': fwd,
            'universal_reverse_5_3': UNIVERSAL_REVERSE,
            'rt_conditions': '16 C 30 min; then 60 cycles of '
                             '(30 C 30 s, 42 C 30 s, 50 C 1 s); 85 C 5 min',
            'reference': 'Varkonyi-Gasic et al. 2007, Plant Methods 3:12'}


def northern_probe(small_rna: str) -> Dict[str, str]:
    """DNA probe complementary to the small RNA, for Northern hybridisation."""
    s = to_dna(small_rna)
    p = rc_dna(s)
    return {'probe_5_3': p,
            'length': len(p),
            'gc_percent': round(gc(p), 1),
            'note': 'End-label with T4 polynucleotide kinase and gamma-32P '
                    'ATP, or order with a digoxigenin / biotin modification. '
                    'LNA-modified probes markedly improve sensitivity for '
                    'small RNAs.'}


# ============================================================================
# COMPLETE ORDER SHEET
# ============================================================================

@dataclass
class OrderSheet:
    name: str
    guide: str
    forms: Dict[str, Dict] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)


def build_order_sheet(name: str, guide: str,
                      target_mrna: Optional[str] = None,
                      loop_key: str = 'standard_12',
                      overhang: str = 'dTdT',
                      include: Optional[Tuple[str, ...]] = None
                      ) -> OrderSheet:
    """
    Generates every laboratory form for a single candidate.

    include: subset of
      'duplex', 'shrna', 'cloning', 'ivt', 'detection', 'controls'
      None means all of them.
    """
    wanted = include or ('duplex', 'shrna', 'cloning', 'ivt',
                         'detection', 'controls')
    sheet = OrderSheet(name=name, guide=to_rna(guide))

    if 'duplex' in wanted:
        sheet.forms['duplex'] = synthetic_duplex(guide, overhang)
        sheet.forms['duplex']['diagram'] = duplex_diagram(guide, overhang)

    if 'shrna' in wanted:
        sheet.forms['shrna'] = shrna(guide, loop_key)

    if 'cloning' in wanted:
        h = sheet.forms.get('shrna') or shrna(guide, loop_key)
        sheet.forms['cloning'] = annealing_oligos(h['shrna_dna'])

    if 'ivt' in wanted:
        sheet.forms['ivt'] = ivt_template(guide)

    if 'detection' in wanted:
        sheet.forms['detection'] = {
            **stem_loop_rt_primer(guide),
            **{f'northern_{k}': v for k, v in northern_probe(guide).items()},
        }

    if 'controls' in wanted:
        scr = scrambled_control(guide, target_mrna)
        ctrl = {
            'c911_5_3': c911_control(guide),
            'c911_note': 'RECOMMENDED specificity control. Retains the seed '
                         'region and therefore off-target activity, while '
                         'central mismatches abolish target cleavage '
                         '(Buehler et al. 2012).',
            'c10_5_3': c10_control(guide),
            'seed_mutant_5_3': seed_mutant_control(guide),
            'seed_mutant_note': ('Complement to C911 - abolishes '
                                 'seed-mediated off-target activity while '
                                 'leaving the 3-prime region intact.'),
            'scrambled_5_3': scr or 'not generated',
            'scrambled_note': 'Conventional but weak. Buehler et al. (2012) '
                              'showed that scrambled controls lose activity '
                              'in both true- and false-positive groups and '
                              'therefore cannot distinguish between them.',
        }
        sheet.forms['controls'] = ctrl
        if scr is None:
            sheet.warnings.append(
                'No scrambled sequence satisfying the constraints was found. '
                'Use the C911 control, which is preferable in any case.')

    if target_mrna is None:
        sheet.warnings.append(
            'Target mRNA not supplied - the scrambled control could not be '
            'checked for accidental complementarity to the target.')

    return sheet


def order_sheet_to_rows(sheet: OrderSheet) -> List[Dict[str, str]]:
    """Flattens an order sheet into rows suitable for a spreadsheet."""
    rows = []
    opis = {
        'guide_5_3': ('Guide strand (antisense)', 'RNA'),
        'passenger_5_3': ('Passenger strand (sense)', 'RNA'),
        'shrna_dna': ('shRNA hairpin', 'DNA'),
        'top_oligo_5_3': ('Cloning oligo, top', 'DNA'),
        'bottom_oligo_5_3': ('Cloning oligo, bottom', 'DNA'),
        'sense_template_5_3': ('IVT template, sense (T7)', 'DNA'),
        'antisense_template_5_3': ('IVT template, antisense (T7)', 'DNA'),
        'rt_primer_5_3': ('Stem-loop RT primer', 'DNA'),
        'forward_primer_5_3': ('qPCR forward primer', 'DNA'),
        'universal_reverse_5_3': ('qPCR universal reverse primer', 'DNA'),
        'northern_probe_5_3': ('Northern blot probe', 'DNA'),
        'c911_5_3': ('C911 control (recommended)', 'RNA'),
        'c10_5_3': ('C10 control', 'RNA'),
        'seed_mutant_5_3': ('Seed-mutant control', 'RNA'),
        'scrambled_5_3': ('Scrambled control (weak)', 'RNA'),
    }
    for grupa, zawartosc in sheet.forms.items():
        for klucz, wartosc in zawartosc.items():
            if klucz not in opis or not isinstance(wartosc, str):
                continue
            etykieta, typ = opis[klucz]
            rows.append({
                'candidate': sheet.name,
                'category': grupa,
                'form': etykieta,
                'chemistry': typ,
                'sequence_5_3': wartosc,
                'length_nt': len(wartosc.replace('dT', 'T')),
                'gc_percent': round(gc(wartosc.replace('dT', 'T')), 1),
            })
    return rows


if __name__ == '__main__':
    guide = 'UUGAAGUUCACCUUGAUGCCG'
    GFP = ('ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGAC'
           'GGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTAC'
           'GGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACC')

    sheet = build_order_sheet('A-21', guide, target_mrna=GFP)

    print('DUPLEX')
    print(sheet.forms['duplex']['diagram'])
    print()
    print('CONTROLS')
    print(f"  parent : {sheet.guide}")
    print(f"  C911   : {sheet.forms['controls']['c911_5_3']}")
    print(f"  seed   : {sheet.forms['controls']['seed_mutant_5_3']}")
    print(f"  scram  : {sheet.forms['controls']['scrambled_5_3']}")
    print()
    print('ORDER SHEET')
    for r in order_sheet_to_rows(sheet):
        print(f"  {r['form']:<32} {r['chemistry']:<4} "
              f"{r['length_nt']:>3} nt  {r['sequence_5_3'][:44]}")
    for w in sheet.warnings:
        print(f'  [!] {w}')
