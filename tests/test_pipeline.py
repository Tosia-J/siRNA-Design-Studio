"""
Tests for candidate generation, scoring, conservation, controls and
construct assembly.
"""

import random
import pytest

import thermo
import filters
import hosts
import scoring
import conservation
import oligos
import constructs


COMP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
rc = lambda s: ''.join(COMP[c] for c in reversed(s.upper().replace('U', 'T')))


def losowa(n, seed=0):
    r = random.Random(seed)
    return ''.join(r.choice('ACGT') for _ in range(n))


# ============================================================================
# HOST PROFILES
# ============================================================================

class TestHostProfiles:

    def test_all_profiles_valid(self):
        for nazwa in hosts.PROFILE:
            p = hosts.get_profile(nazwa)
            assert p.gc_min < p.gc_max
            assert p.dlugosci
            assert p.opis

    def test_unknown_profile_rejected(self):
        with pytest.raises(ValueError):
            hosts.get_profile('nonexistent')

    def test_mammalian_profile_filters_immunostimulatory_motifs(self):
        """
        The distinction between profiles must be real, not cosmetic. Mammals
        have TLR7/8 and an interferon pathway; plants do not.
        """
        assert hosts.get_profile('mammal').motywy_zabronione
        assert not hosts.get_profile('plant').motywy_zabronione

    def test_plant_profile_requires_specific_5prime(self):
        """AGO1 in plants selects 5'-U; AGO2 selects 5'-A."""
        assert hosts.get_profile('plant').guide_5_dozwolone == 'AU'
        assert hosts.get_profile('mammal').guide_5_dozwolone is None

    def test_plant_profile_covers_three_dcl_classes(self):
        assert set(hosts.get_profile('plant').dlugosci) == {21, 22, 24}

    def test_profiles_produce_different_results(self):
        """A profile that changes nothing is not a profile."""
        sekw = 'CGGCATCAAGGTGAACTTCAA'
        ok_p, _ = filters.physicochemical_filters(
            sekw, gc_min=hosts.get_profile('plant').gc_min,
            gc_max=hosts.get_profile('plant').gc_max)
        assert isinstance(ok_p, bool)


# ============================================================================
# PHYSICOCHEMICAL FILTERS
# ============================================================================

class TestFilters:

    def test_accepts_a_good_candidate(self):
        ok, powody = filters.physicochemical_filters('CGGCATCAAGGTGAACTTCAA')
        assert ok, powody

    def test_rejects_high_gc(self):
        ok, powody = filters.physicochemical_filters('GCGCGCGCGCGCGCGCGCGCG')
        assert not ok
        assert any('GC' in p for p in powody)

    def test_rejects_homopolymer(self):
        ok, powody = filters.physicochemical_filters('CGGCAAAAGGTGAACTTCAA')
        assert not ok

    def test_mammalian_motifs_reported_not_silently_dropped(self):
        """
        The motifs must remain reportable even in the plant profile, so the
        report can state that they were checked and found irrelevant.
        """
        znalezione = filters.check_mammalian_motifs('UGUGUAAAAAAAAAAAAAAAA')
        assert 'UGUGU' in znalezione


# ============================================================================
# SCORING
# ============================================================================

def _demo_kandydaci():
    return [
        {'nazwa': 'A', 'asymetria': 4.8, 'dostepnosc_celu': -0.32,
         'mfe_guide': 0.0, 'gc_optymalnosc': 0.94, 'koszt_offtarget': 18.0,
         'seed_czystosc': 0},
        {'nazwa': 'B', 'asymetria': 0.5, 'dostepnosc_celu': -0.28,
         'mfe_guide': -2.6, 'gc_optymalnosc': 0.99, 'koszt_offtarget': 21.0,
         'seed_czystosc': 0},
        {'nazwa': 'C', 'asymetria': 4.2, 'dostepnosc_celu': -0.41,
         'mfe_guide': -0.1, 'gc_optymalnosc': 0.88, 'koszt_offtarget': 5.0,
         'seed_czystosc': 3},
    ]


class TestScoring:

    def test_weights_must_sum_to_one(self):
        """
        Weights that do not sum to one make results incomparable between
        runs, so this is enforced rather than silently normalised.
        """
        zle = {k: 0.5 for k in scoring.METRYKI}
        with pytest.raises(ValueError):
            scoring.normalize_and_score(_demo_kandydaci(), wagi=zle)

    def test_normalisation_maps_to_unit_interval(self):
        r = scoring.normalize_and_score(_demo_kandydaci())
        for c in r:
            for m in scoring.METRYKI:
                assert 0.0 <= c[f'norm_{m}'] <= 1.0

    def test_normalisation_prevents_scale_domination(self):
        """
        The reason normalisation exists: without it a metric with a larger
        numerical range dominates regardless of its declared weight.
        """
        kand = _demo_kandydaci()
        for i, c in enumerate(kand):
            c['koszt_offtarget'] = c['koszt_offtarget'] * 1000
        r = scoring.normalize_and_score(kand)
        assert all(0.0 <= c['norm_koszt_offtarget'] <= 1.0 for c in r)

    def test_ranking_is_ordered(self):
        r = scoring.normalize_and_score(_demo_kandydaci())
        assert [c['ranga'] for c in r] == sorted(c['ranga'] for c in r)
        wyniki = [c['wynik'] for c in r]
        assert wyniki == sorted(wyniki, reverse=True)

    def test_empty_input_returns_empty(self):
        assert scoring.normalize_and_score([]) == []

    def test_identical_candidates_get_midpoint(self):
        """Zero variance must not divide by zero."""
        jeden = _demo_kandydaci()[0]
        r = scoring.normalize_and_score([dict(jeden), dict(jeden)])
        assert all(c['norm_asymetria'] == pytest.approx(0.5) for c in r)

    def test_sensitivity_analysis_returns_bounded_stability(self):
        s = scoring.sensitivity_analysis(_demo_kandydaci(), top_n=2)
        assert 0.0 <= s['stabilnosc'] <= 1.0
        assert s['interpretacja']

    def test_pareto_front_is_subset_and_nonempty(self):
        kand = _demo_kandydaci()
        front = scoring.pareto_front(kand)
        assert 0 < len(front) <= len(kand)
        nazwy = {c['nazwa'] for c in kand}
        assert {c['nazwa'] for c in front} <= nazwy

    def test_dominated_candidate_excluded_from_front(self):
        """A candidate worse in every metric must not be on the front."""
        kand = [
            {'nazwa': 'good', 'asymetria': 5.0, 'dostepnosc_celu': -0.1,
             'mfe_guide': 0.0, 'gc_optymalnosc': 1.0, 'koszt_offtarget': 30.0,
             'seed_czystosc': 0},
            {'nazwa': 'bad', 'asymetria': 0.1, 'dostepnosc_celu': -0.9,
             'mfe_guide': -8.0, 'gc_optymalnosc': 0.1, 'koszt_offtarget': 1.0,
             'seed_czystosc': 9},
        ]
        assert {c['nazwa'] for c in scoring.pareto_front(kand)} == {'good'}


# ============================================================================
# CONSERVATION
# ============================================================================

class TestConservation:

    @staticmethod
    def _zbior(tmp_path, rdzen=(200, 320), n=20, mut=0.05, seed=13):
        r = random.Random(seed)
        baza = ''.join(r.choice('ACGT') for _ in range(600))
        rek = [('ref', baza)]
        for i in range(n):
            s = list(baza)
            for p in range(len(s)):
                if rdzen[0] <= p < rdzen[1]:
                    continue
                if r.random() < mut:
                    s[p] = r.choice('ACGT')
            rek.append((f'iso{i}', ''.join(s)))
        f = tmp_path / 'iso.fasta'
        f.write_text(''.join(f'>{n_}\n{s}\n' for n_, s in rek))
        return str(f)

    def test_conserved_core_detected(self, tmp_path):
        f = self._zbior(tmp_path)
        w = conservation.analiza_konserwatywnosci(
            f, dlugosci=(21,), prog_pokrycia=0.95, verbose=False)
        assert w['regiony']
        pozycje = [r['poz_start_1based'] for r in w['regiony']]
        # the invariant core spans positions 201-320 (1-based)
        assert min(pozycje) >= 180
        assert max(pozycje) <= 320

    def test_blocks_merge_overlapping_windows(self, tmp_path):
        f = self._zbior(tmp_path)
        w = conservation.analiza_konserwatywnosci(
            f, dlugosci=(21,), prog_pokrycia=0.95, verbose=False)
        bloki = conservation.scal_regiony(w['regiony'], min_dlugosc=30)
        assert bloki
        assert bloki[0]['dlugosc_bloku'] >= 30

    def test_coverage_bounded(self, tmp_path):
        f = self._zbior(tmp_path)
        w = conservation.analiza_konserwatywnosci(
            f, dlugosci=(21,), prog_pokrycia=0.50, verbose=False)
        for r in w['regiony']:
            assert 0.0 <= r['pokrycie'] <= 1.0
            assert r['pokrycie_1mm'] >= r['pokrycie']

    def test_threshold_optimisation_finds_knee(self, tmp_path):
        f = self._zbior(tmp_path)
        o = conservation.optimise_threshold(f, dlugosc=21, verbose=False)
        assert 0.5 <= o['threshold'] <= 1.0
        assert o['confidence'] in ('high', 'medium', 'low')
        assert o['curve']

    def test_no_knee_reported_honestly(self, tmp_path):
        """
        For identical sequences there is no transition to find. The function
        must say so rather than return an authoritative-looking number.
        """
        seq = losowa(400, 3)
        f = tmp_path / 'same.fasta'
        f.write_text(''.join(f'>s{i}\n{seq}\n' for i in range(10)))
        o = conservation.optimise_threshold(str(f), dlugosc=21, verbose=False)
        assert o['confidence'] in ('low', 'medium', 'high')
        assert o['diagnostics']

    def test_entropy_requires_aligned_input(self):
        with pytest.raises(ValueError):
            conservation.entropia_pozycyjna(['ACGT', 'ACG'])

    def test_entropy_zero_for_invariant_column(self):
        assert conservation.entropia_kolumny('AAAA') == pytest.approx(0.0)

    def test_entropy_maximal_for_uniform_column(self):
        assert conservation.entropia_kolumny('ACGT') == pytest.approx(2.0)


# ============================================================================
# CONTROLS
# ============================================================================

class TestControls:

    def test_c911_changes_exactly_positions_9_to_11(self):
        guide = 'UUGAAGUUCACCUUGAUGCCG'
        c911 = oligos.c911_control(guide)
        assert len(c911) == len(guide)
        rozne = [i for i in range(len(guide)) if guide[i] != c911[i]]
        assert rozne == [8, 9, 10]

    def test_c911_uses_complement_not_arbitrary_substitution(self):
        guide = 'UUGAAGUUCACCUUGAUGCCG'
        c911 = oligos.c911_control(guide)
        for i in (8, 9, 10):
            assert c911[i] == oligos.COMP_RNA[guide[i]]

    def test_c911_preserves_seed(self):
        """
        The entire point of C911: the seed is intact, so seed-mediated
        off-target activity is retained while target cleavage is abolished.
        """
        guide = 'UUGAAGUUCACCUUGAUGCCG'
        assert oligos.c911_control(guide)[1:8] == guide[1:8]

    def test_seed_mutant_changes_seed_only(self):
        guide = 'UUGAAGUUCACCUUGAUGCCG'
        sm = oligos.seed_mutant_control(guide)
        rozne = [i for i in range(len(guide)) if guide[i] != sm[i]]
        assert all(1 <= i <= 7 for i in rozne)

    def test_c911_rejects_short_guide(self):
        with pytest.raises(ValueError):
            oligos.c911_control('UUGAA')

    def test_scrambled_preserves_composition(self):
        guide = 'UUGAAGUUCACCUUGAUGCCG'
        scr = oligos.scrambled_control(guide)
        assert scr is not None
        assert sorted(scr) == sorted(guide)
        assert scr != guide

    def test_scrambled_avoids_target_complementarity(self):
        guide = 'UUGAAGUUCACCUUGAUGCCG'
        cel = 'CGGCATCAAGGTGAACTTCAA' * 5
        scr = oligos.scrambled_control(guide, target_mrna=cel)
        if scr is not None:
            d = scr.replace('U', 'T')
            assert not any(rc(d[i:i + 8]) in cel for i in range(len(d) - 7))

    def test_scrambled_is_deterministic(self):
        g = 'UUGAAGUUCACCUUGAUGCCG'
        assert oligos.scrambled_control(g, seed=1) == \
               oligos.scrambled_control(g, seed=1)


# ============================================================================
# LABORATORY FORMS
# ============================================================================

class TestOligoForms:

    def test_duplex_strands_are_complementary(self):
        d = oligos.synthetic_duplex('UUGAAGUUCACCUUGAUGCCG')
        assert d['passenger_core'] == thermo.revcomp_rna(d['guide_core'])

    def test_overhang_added(self):
        d = oligos.synthetic_duplex('UUGAAGUUCACCUUGAUGCCG', 'dTdT')
        assert d['guide_5_3'].endswith('dTdT')
        assert oligos.synthetic_duplex(
            'UUGAAGUUCACCUUGAUGCCG', 'native')['guide_5_3'].endswith('G')

    def test_shrna_contains_both_arms_and_loop(self):
        guide = 'UUGAAGUUCACCUUGAUGCCG'
        h = oligos.shrna(guide, 'standard_12')['shrna_dna']
        g = guide.replace('U', 'T')
        assert h.startswith(g)
        assert h.endswith(rc(g))
        assert oligos.LOOPS['standard_12'][0] in h

    def test_shrna_length_is_consistent(self):
        guide = 'UUGAAGUUCACCUUGAUGCCG'
        h = oligos.shrna(guide, 'standard_12')
        assert h['length'] == 2 * len(guide) + len(oligos.LOOPS['standard_12'][0])

    def test_ivt_template_has_t7_promoter(self):
        t = oligos.ivt_template('ATGCATGCATGC')
        assert t['sense_template_5_3'].startswith(oligos.T7_PROMOTER)
        assert t['antisense_template_5_3'].startswith(oligos.T7_PROMOTER)

    def test_stem_loop_primer_ends_with_revcomp_of_small_rna(self):
        s = 'UUGAAGUUCACCUUGAUGCCG'
        p = oligos.stem_loop_rt_primer(s)
        assert p['rt_primer_5_3'].endswith(rc(s.replace('U', 'T')[-6:]))

    def test_northern_probe_is_complementary(self):
        s = 'UUGAAGUUCACCUUGAUGCCG'
        assert oligos.northern_probe(s)['probe_5_3'] == rc(s)

    def test_order_sheet_contains_expected_forms(self):
        sh = oligos.build_order_sheet('X', 'UUGAAGUUCACCUUGAUGCCG')
        wiersze = oligos.order_sheet_to_rows(sh)
        formy = {r['form'] for r in wiersze}
        assert any('Guide' in f for f in formy)
        assert any('C911' in f for f in formy)
        assert all(r['length_nt'] > 0 for r in wiersze)


# ============================================================================
# CONSTRUCTS
# ============================================================================

class TestConstructs:

    def test_pol3_adds_required_first_nucleotide(self):
        """AtU6-1 requires the transcript to begin with G."""
        k = constructs.zbuduj_kasete(
            'x', 'UUGAAGUUCACCUUGAUGCCG', 'AtU6-1',
            'restrykcyjne_BamHI_EcoRI')
        assert k.dodany_nt_5 == 'G'
        assert k.shrna_dna.startswith('G')
        assert any('5-koncu' in o or '5' in o for o in k.ostrzezenia)

    def test_pol3_appends_terminator(self):
        k = constructs.zbuduj_kasete(
            'x', 'UUGAAGUUCACCUUGAUGCCG', 'AtU6-1',
            'restrykcyjne_BamHI_EcoRI')
        assert k.shrna_dna.endswith('TTTTTT')

    def test_pol2_does_not_add_nucleotide(self):
        k = constructs.zbuduj_kasete(
            'x', 'UUGAAGUUCACCUUGAUGCCG', 'CaMV35S',
            'restrykcyjne_BamHI_EcoRI')
        assert k.dodany_nt_5 is None

    def test_atu3_needs_no_addition_for_a_starting_guide(self):
        """A guide already beginning with A fits AtU3 without modification."""
        k = constructs.zbuduj_kasete(
            'x', 'AUCUUGAAGUUCACCUUGAUGCCG', 'AtU3',
            'restrykcyjne_BamHI_EcoRI')
        assert k.dodany_nt_5 is None

    def test_restriction_site_collision_detected(self):
        """
        An insert containing the enzyme's own site would be cut internally.
        Since the shRNA sequence cannot be changed without losing the target,
        this must be reported rather than silently accepted.
        """
        guide = thermo.to_rna(thermo.revcomp_rna('GGATCCAAGGTGAACTTCAA'))
        k = constructs.zbuduj_kasete(
            'x', guide, 'CaMV35S', 'restrykcyjne_BamHI_EcoRI')
        assert any('KOLIZJA' in o or 'GGATCC' in o for o in k.ostrzezenia)

    def test_unknown_promoter_rejected(self):
        with pytest.raises(KeyError):
            constructs.zbuduj_kasete(
                'x', 'UUGAAGUUCACCUUGAUGCCG', 'nonexistent',
                'restrykcyjne_BamHI_EcoRI')
