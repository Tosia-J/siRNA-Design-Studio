"""
Tests for off-target analysis.

Several of these are regression tests for defects that were present in
earlier versions and are easy to reintroduce:

  - the intended target counted as an off-target hit, which rejected every
    candidate whenever the target gene was present in the reference
  - a minimum cost taken over a large database, which is low by chance and
    therefore rejected every candidate at genome scale
  - a full index of the transcriptome held in memory, which exhausted RAM
"""

import random
import pytest

import offtarget


COMP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}


def rc(s):
    return ''.join(COMP[c] for c in reversed(s.upper().replace('U', 'T')))


def zapisz_fasta(tmp_path, rekordy, nazwa='ref.fasta'):
    p = tmp_path / nazwa
    p.write_text(''.join(f'>{n}\n{s}\n' for n, s in rekordy))
    return str(p)


def losowa(n, seed=0):
    r = random.Random(seed)
    return ''.join(r.choice('ACGT') for _ in range(n))


# ============================================================================
# BASIC MECHANICS
# ============================================================================

class TestProbeConstruction:

    def test_probe_count_scales_with_queries_not_database(self):
        """
        The defining property of the streaming design: memory depends on the
        number of queries, not on the size of the reference. If this breaks,
        large transcriptomes will exhaust RAM again.
        """
        jeden = offtarget.OffTargetScanner({'a': losowa(21, 1)})
        dziesiec = offtarget.OffTargetScanner(
            {f'q{i}': losowa(21, i) for i in range(10)})
        assert dziesiec.liczba_sond < jeden.liczba_sond * 15
        assert jeden.liczba_sond < 100

    def test_positional_weights_ordering(self):
        """Seed positions must outweigh the 3' region."""
        assert offtarget.position_weight(3) > offtarget.position_weight(15)
        assert offtarget.position_weight(10) > offtarget.position_weight(15)
        assert offtarget.position_weight(1) < offtarget.position_weight(3)

    def test_seed_positions_carry_highest_weight(self):
        seed = [offtarget.position_weight(p) for p in range(2, 9)]
        reszta = [offtarget.position_weight(p) for p in range(12, 22)]
        assert min(seed) > max(reszta)


# ============================================================================
# DETECTION
# ============================================================================

class TestDetection:

    def test_perfect_match_detected(self, tmp_path):
        cel = losowa(200, 5)
        guide = rc(cel[50:71])
        f = zapisz_fasta(tmp_path, [('tx', cel)])
        r = offtarget.analyze_guides({'g': guide}, f)['g']
        assert r['n_trafien_pelnych'] >= 1
        assert r['min_koszt'] == pytest.approx(0.0)

    def test_absent_sequence_not_detected(self, tmp_path):
        f = zapisz_fasta(tmp_path, [('tx', 'ATGC' * 100)])
        r = offtarget.analyze_guides({'g': 'UUUUUCCCCCAAAAAGGGGGA'}, f)['g']
        assert r['n_trafien_pelnych'] == 0
        assert r['flaga'] == 'OK'

    def test_seed_match_counted(self, tmp_path):
        guide = 'UUGAAGUUCACCUUGAUGCCG'
        seed_2_8 = guide.replace('U', 'T')[1:8]
        kontekst = 'GGGGCCCC' + rc(seed_2_8) + 'A' + 'TTTTCCCC'
        f = zapisz_fasta(tmp_path, [('tx', kontekst)])
        r = offtarget.analyze_guides({'g': guide}, f)['g']
        assert r['n_8mer'] >= 1

    def test_empty_reference_handled(self, tmp_path):
        f = zapisz_fasta(tmp_path, [('tx', 'ACGT')])
        r = offtarget.analyze_guides({'g': losowa(21, 3)}, f)['g']
        assert r['n_trafien_pelnych'] == 0
        assert r['flaga'] == 'OK'


# ============================================================================
# REGRESSION: intended target must not count as off-target
# ============================================================================

class TestIntendedTargetExclusion:
    """
    Earlier versions rejected every candidate when the target gene was
    present in the reference, because a guide always matches its own target
    at zero cost.
    """

    def test_without_target_mrna_own_target_looks_like_offtarget(self,
                                                                 tmp_path):
        cel = losowa(300, 7)
        guide = rc(cel[100:121])
        f = zapisz_fasta(tmp_path, [('target_gene', cel)])
        r = offtarget.analyze_guides({'g': guide}, f)['g']
        assert r['flaga'].startswith('REJECT')

    def test_with_target_mrna_own_target_is_excluded(self, tmp_path):
        cel = losowa(300, 7)
        guide = rc(cel[100:121])
        f = zapisz_fasta(tmp_path, [('target_gene', cel)])
        r = offtarget.analyze_guides({'g': guide}, f, mrna_celu=cel)['g']
        assert r['flaga'] == 'OK'
        assert r['n_trafien_zamierzonych'] >= 1
        assert r['n_trafien_pelnych'] == 0

    def test_genuine_offtarget_still_detected_with_target_supplied(self,
                                                                   tmp_path):
        """
        Excluding the intended target must not suppress a real off-target
        hit in a different transcript.
        """
        cel = losowa(300, 11)
        inny = losowa(300, 12)
        guide = rc(cel[100:121])
        # plant the same site in an unrelated transcript
        inny = inny[:150] + cel[100:121] + inny[171:]
        f = zapisz_fasta(tmp_path, [('target_gene', cel), ('other', inny)])
        r = offtarget.analyze_guides({'g': guide}, f, mrna_celu=cel)['g']
        assert r['n_trafien_zamierzonych'] >= 1
        assert r['n_trafien_pelnych'] >= 1
        assert r['flaga'].startswith('REJECT')


# ============================================================================
# REGRESSION: multiple-testing correction
# ============================================================================

class TestChanceCorrection:
    """
    An exact 8-mer occurs by chance roughly once every 65 536 positions.
    Rejecting on any seed hit, or on a minimum cost taken over a large
    database, rejects everything at genome scale.
    """

    def test_chance_seed_hits_do_not_reject(self, tmp_path):
        baza = [(f'tx{i}', losowa(2000, 100 + i)) for i in range(60)]
        f = zapisz_fasta(tmp_path, baza)
        guide = losowa(21, 999).replace('T', 'U')
        r = offtarget.analyze_guides({'g': guide}, f)['g']
        assert r['n_8mer'] > 0, 'expected chance hits in 120 kb'
        assert not r['flaga'].startswith('REJECT')

    def test_enrichment_near_unity_for_random_sequence(self, tmp_path):
        baza = [(f'tx{i}', losowa(4000, 200 + i)) for i in range(40)]
        f = zapisz_fasta(tmp_path, baza)
        guide = losowa(21, 777).replace('T', 'U')
        r = offtarget.analyze_guides({'g': guide}, f)['g']
        assert r['wzbogacenie_8mer'] is not None
        assert 0.2 < r['wzbogacenie_8mer'] < 5.0

    def test_expected_count_matches_theory(self, tmp_path):
        """Expected chance hits should be database length divided by 4^8."""
        dl = 4000 * 40
        baza = [(f'tx{i}', losowa(4000, 300 + i)) for i in range(40)]
        f = zapisz_fasta(tmp_path, baza)
        r = offtarget.analyze_guides({'g': losowa(21, 5).replace('T', 'U')},
                                     f)['g']
        assert r['oczek_8mer_losowo'] == pytest.approx(dl / 65536, rel=0.02)

    def test_critical_hit_definition(self, tmp_path):
        """
        A critical hit requires at most three mismatches AND an intact seed.
        A hit with mismatches inside the seed must not be critical.
        """
        cel = losowa(300, 21)
        okno = cel[100:121]
        guide = rc(okno)
        # mutate a position inside the seed of the guide, i.e. near the 3'
        # end of the target window
        lista = list(cel)
        lista[100 + 17] = {'A': 'C', 'C': 'A', 'G': 'T', 'T': 'G'}[lista[100 + 17]]
        zmieniony = ''.join(lista)
        f = zapisz_fasta(tmp_path, [('tx', zmieniony)])
        r = offtarget.analyze_guides({'g': guide}, f)['g']
        assert r['n_krytycznych'] == 0


# ============================================================================
# DETERMINISM
# ============================================================================

class TestDeterminism:

    def test_same_input_gives_same_result(self, tmp_path):
        baza = [(f'tx{i}', losowa(1000, 400 + i)) for i in range(10)]
        f = zapisz_fasta(tmp_path, baza)
        g = {'g': losowa(21, 42).replace('T', 'U')}
        a = offtarget.analyze_guides(g, f)['g']
        b = offtarget.analyze_guides(g, f)['g']
        for k in ('n_8mer', 'n_trafien_pelnych', 'min_koszt', 'flaga'):
            assert a[k] == b[k]

    def test_record_order_does_not_change_counts(self, tmp_path):
        baza = [(f'tx{i}', losowa(1000, 500 + i)) for i in range(8)]
        f1 = zapisz_fasta(tmp_path, baza, 'a.fasta')
        f2 = zapisz_fasta(tmp_path, list(reversed(baza)), 'b.fasta')
        g = {'g': losowa(21, 43).replace('T', 'U')}
        a = offtarget.analyze_guides(g, f1)['g']
        b = offtarget.analyze_guides(g, f2)['g']
        assert a['n_8mer'] == b['n_8mer']
        assert a['n_trafien_pelnych'] == b['n_trafien_pelnych']
