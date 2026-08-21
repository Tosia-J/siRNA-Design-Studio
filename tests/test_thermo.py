"""
Tests for the thermodynamic module.

Reference values come from Xia et al. 1998, Biochemistry 37:14719-14735,
Table 4, and from the identities that any correct nearest-neighbour
implementation must satisfy.
"""

import math
import pytest

import thermo


# ============================================================================
# PARAMETER TABLE INTEGRITY
# ============================================================================

class TestParameterTable:
    """The tables must be complete and internally consistent."""

    def test_all_sixteen_dinucleotides_present(self):
        oczekiwane = {a + b for a in 'ACGU' for b in 'ACGU'}
        assert set(thermo.NN_DG37) == oczekiwane
        assert set(thermo.NN_DH) == oczekiwane

    @pytest.mark.parametrize('dinuk', sorted(
        {a + b for a in 'ACGU' for b in 'ACGU'}))
    def test_reverse_complement_symmetry(self, dinuk):
        """
        A nearest-neighbour stack read from either strand is the same
        physical stack, so its parameters must be identical. Violating this
        would make the free energy depend on which strand is supplied.
        """
        rc = thermo.revcomp_rna(dinuk)
        assert thermo.NN_DG37[dinuk] == pytest.approx(thermo.NN_DG37[rc])
        assert thermo.NN_DH[dinuk] == pytest.approx(thermo.NN_DH[rc])

    def test_all_stacks_are_stabilising(self):
        """Every Watson-Crick stack lowers free energy."""
        assert all(v < 0 for v in thermo.NN_DG37.values())

    def test_gc_stack_is_most_stable(self):
        """
        GC/GC is the strongest stack and AA/UU the weakest; this ordering is
        a basic property of the Turner model.
        """
        assert min(thermo.NN_DG37, key=thermo.NN_DG37.get) == 'GC'
        assert max(thermo.NN_DG37, key=thermo.NN_DG37.get) in ('AA', 'UU')

    def test_known_values_from_xia_table_4(self):
        """Spot checks against the published table."""
        assert thermo.NN_DG37['GC'] == pytest.approx(-3.42)
        assert thermo.NN_DG37['CG'] == pytest.approx(-2.36)
        assert thermo.NN_DG37['AA'] == pytest.approx(-0.93)
        assert thermo.NN_DG37['UA'] == pytest.approx(-1.33)
        assert thermo.INIT_DG == pytest.approx(4.09)
        assert thermo.TERM_AU_DG == pytest.approx(0.45)


# ============================================================================
# SEQUENCE HANDLING
# ============================================================================

class TestSequenceUtilities:

    def test_to_rna_converts_thymine(self):
        assert thermo.to_rna('ATGC') == 'AUGC'
        assert thermo.to_rna('atgc') == 'AUGC'

    def test_revcomp_is_an_involution(self):
        seq = 'AUGCAUGCAUGC'
        assert thermo.revcomp_rna(thermo.revcomp_rna(seq)) == seq

    def test_revcomp_known_case(self):
        assert thermo.revcomp_rna('AUGC') == 'GCAU'

    def test_gc_content(self):
        assert thermo.gc_content('GGCC') == pytest.approx(100.0)
        assert thermo.gc_content('AAUU') == pytest.approx(0.0)
        assert thermo.gc_content('AUGC') == pytest.approx(50.0)

    def test_homopolymer_detection(self):
        assert thermo.has_homopolymer('AAAAG', 4) is True
        assert thermo.has_homopolymer('AAAG', 4) is False
        assert thermo.has_homopolymer('GGGGA', 4) is True

    @pytest.mark.parametrize('bad', ['AUGX', 'AUG N', '123', 'AUGCZ'])
    def test_invalid_characters_rejected(self, bad):
        with pytest.raises(ValueError):
            thermo.duplex_dg(bad)

    def test_sequence_too_short_rejected(self):
        with pytest.raises(ValueError):
            thermo.duplex_dg('A')


# ============================================================================
# FREE ENERGY
# ============================================================================

class TestDuplexEnergy:

    def test_gc_rich_duplex_more_stable_than_au_rich(self):
        gc = thermo.duplex_dg('GCGCGCGCGCGCGCGCGCGCG')
        au = thermo.duplex_dg('AUAUAUAUAUAUAUAUAUAUA')
        assert gc < au

    def test_longer_duplex_more_stable(self):
        assert thermo.duplex_dg('GCGCGCGCGC') < thermo.duplex_dg('GCGCG')

    def test_energy_independent_of_strand_supplied(self):
        """
        The duplex has one free energy. Supplying either strand must give the
        same value - this follows from the symmetry tested above and is the
        property most likely to break if the table is edited carelessly.
        """
        seq = 'AUGCCGGAUUCGAUCGAUGCA'
        assert thermo.duplex_dg(seq) == pytest.approx(
            thermo.duplex_dg(thermo.revcomp_rna(seq)))

    def test_terminal_au_penalty_applied(self):
        """A duplex ending in A-U is penalised relative to one ending in G-C."""
        z_kara = thermo.duplex_dg('AGCGCGCGCU', include_init=True)
        bez_kary = thermo.duplex_dg('GGCGCGCGCC', include_init=True)
        stacki_a = sum(thermo.NN_DG37['AGCGCGCGCU'[i:i + 2]] for i in range(9))
        stacki_g = sum(thermo.NN_DG37['GGCGCGCGCC'[i:i + 2]] for i in range(9))
        assert (z_kara - stacki_a) - (bez_kary - stacki_g) == pytest.approx(
            2 * thermo.TERM_AU_DG)

    def test_manual_calculation_matches(self):
        """
        Full hand calculation for a short duplex, so that a change in the
        summation logic is caught even if the table is intact.

        5'-GCAU-3' : stacks GC, CA, AU
        """
        oczekiwane = (thermo.NN_DG37['GC'] + thermo.NN_DG37['CA']
                      + thermo.NN_DG37['AU'] + thermo.INIT_DG
                      + thermo.TERM_AU_DG)   # U at the 3' end only
        assert thermo.duplex_dg('GCAU') == pytest.approx(oczekiwane)


# ============================================================================
# ASYMMETRY - the metric the whole pipeline depends on
# ============================================================================

class TestAsymmetry:

    def test_definition_sign_convention(self):
        """
        asymmetry = dG5'(guide) - dG5'(passenger)

        A positive value means the guide 5' end is the weaker one, which is
        the configuration in which RISC loads the intended strand.
        """
        guide = 'UUGAAGUUCACCUUGAUGCCG'
        dg_g, dg_p, asym = thermo.asymmetry(guide)
        assert asym == pytest.approx(dg_g - dg_p)

    def test_au_rich_5prime_gives_positive_asymmetry(self):
        """A guide beginning A/U and ending G/C should score positively."""
        _, _, asym = thermo.asymmetry('AUAUAUGCGCGCGCGCGCGCG')
        assert asym > 0

    def test_gc_rich_5prime_gives_negative_asymmetry(self):
        """The reverse arrangement must give a negative value."""
        _, _, asym = thermo.asymmetry('GCGCGCGCGCGCGCGAUAUAU')
        assert asym < 0

    def test_asymmetry_is_antisymmetric(self):
        """
        Swapping guide and passenger must flip the sign exactly. If it does
        not, the two ends are being computed inconsistently.
        """
        guide = 'UUGAAGUUCACCUUGAUGCCG'
        passenger = thermo.revcomp_rna(guide)
        _, _, a1 = thermo.asymmetry(guide)
        _, _, a2 = thermo.asymmetry(passenger)
        assert a1 == pytest.approx(-a2)

    def test_palindrome_has_zero_asymmetry(self):
        """
        A self-complementary sequence has identical ends, so its asymmetry
        must vanish. This is the cleanest available null case.
        """
        _, _, asym = thermo.asymmetry('GCGCGCGCGCGCGC')
        assert asym == pytest.approx(0.0, abs=1e-9)

    def test_known_candidate_value(self):
        """Regression guard for the anti-GFP candidate used throughout."""
        _, _, asym = thermo.asymmetry('UUGAAGUUCACCUUGAUGCCG')
        assert asym == pytest.approx(4.83, abs=0.01)


# ============================================================================
# MELTING TEMPERATURE
# ============================================================================

class TestMeltingTemperature:

    def test_gc_rich_melts_higher(self):
        assert thermo.tm('GCGCGCGCGCGCGCGCGCGCG') > thermo.tm(
            'AUAUAUAUAUAUAUAUAUAUA')

    def test_longer_melts_higher(self):
        assert thermo.tm('GCGCGCGCGCGCGCGCGCGCG') > thermo.tm('GCGCGCGCGC')

    def test_concentration_dependence(self):
        """Tm rises with strand concentration for a bimolecular duplex."""
        seq = 'AUGCCGGAUUCGAUCGAUGCA'
        assert thermo.tm(seq, conc_M=1e-5) > thermo.tm(seq, conc_M=1e-8)

    def test_value_in_physically_plausible_range(self):
        """A 21-mer at 1 uM should melt somewhere between 20 and 100 C."""
        t = thermo.tm('UUGAAGUUCACCUUGAUGCCG')
        assert 20.0 < t < 100.0
