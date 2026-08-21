"""
=============================================================================
hosts.py  --  PROFILE GOSPODARZA
=============================================================================

PO CO TEN MODUL

Pierwsza wersja pipeline'u zakladala roslinnego gospodarza i na tej podstawie
usuwala filtr motywow immunostymulacyjnych. Bylo to poprawne dla roslin, ale
czyni algorytm nieuniwersalnym.

Reguly projektowania siRNA NIE SA uniwersalne miedzy krolestwami. Roznice
dotycza trzech obszarow:

1. ROZPOZNANIE IMMUNOLOGICZNE
   Ssaki: receptory TLR7/TLR8 rozpoznaja jednoniciowe RNA w endosomach,
   TLR3, PKR i RIG-I rozpoznaja dwuniciowe. Uruchamia to odpowiedz
   interferonowa. Zidentyfikowano motywy nasilajace to rozpoznanie:
   UGUGU (Judge i wsp. 2005), GUCCUUCAA (Hornung i wsp. 2005), UGGC.
   -> filtrowac przy celu ssaczym

   Rosliny: brak TLR, brak interferonow, brak PKR i RIG-I. Powyzsze motywy
   sa nieistotne. Wystepuje natomiast dsRNA-PTI, ale jest NIEZALEZNE od
   sekwencji, wiec niefiltrowalne.
   -> nie filtrowac, ale raportowac

   Owady: maja szlaki rozpoznawania dsRNA, ale bez interferonow. Motywy
   ssacze nie maja tu udokumentowanego znaczenia.
   -> nie filtrowac

2. PREFERENCJA 5'-KONCA NICI PROWADZACEJ
   Rosliny: AGO1 preferuje 5'-U, AGO2 preferuje 5'-A, AGO5 preferuje 5'-C.
   Szlak przeciwwirusowy idzie przez AGO1 i AGO2. -> wymagac U/A
   Ssaki: AGO2 nie wykazuje tak silnej selektywnosci MID; regula wynika
   glownie z asymetrii termodynamicznej. -> preferowac A/U, ale nie wymagac
   Owady: Ago2 w szlaku siRNA, preferencje slabiej opisane. -> jak ssaki

3. ZAKRES GC
   Rosliny: wezszy, 30-52%. AGO1 roslinne jest wrazliwsze na stabilnosc
   dupleksu.
   Ssaki: szerszy, 30-64% (Reynolds i wsp. 2004, Ui-Tei i wsp. 2004).

JAK UZYWAC

    from hosts import get_profile
    profil = get_profile('mammal')   # albo 'plant', 'insect', 'generic'

Profil przekazywany jest do filters.physicochemical_filters() i steruje
wlaczaniem poszczegolnych regul.

Autor: Antonina Jarecka
=============================================================================
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


# Motywy immunostymulacyjne u ssakow
MOTYWY_IMMUNO_SSAKI: Tuple[str, ...] = ('UGUGU', 'GUCCUUCAA', 'UGGC')


@dataclass
class HostProfile:
    """Zestaw regul projektowania zalezny od organizmu docelowego."""

    nazwa: str
    opis: str

    # zakres GC nici sensownej
    gc_min: float
    gc_max: float

    # 5'-koniec nici prowadzacej
    guide_5_dozwolone: Optional[str]      # None = brak wymogu
    guide_5_preferowane: str              # uzywane w scoringu, nie jako filtr

    # 5'-koniec nici pasazerskiej
    passenger_5_dozwolone: Optional[str]

    # motywy zabronione (filtr twardy)
    motywy_zabronione: Tuple[str, ...] = ()

    # motywy raportowane, ale niefiltrowane
    motywy_raportowane: Tuple[str, ...] = ()

    # progi termodynamiczne
    min_asymetria: float = 0.3
    max_mfe_guide: float = -3.0

    # dlugosci istotne biologicznie
    dlugosci: Tuple[int, ...] = (21,)
    opis_dlugosci: Dict[int, str] = field(default_factory=dict)

    # wykluczenia koncow ORF
    wyklucz_5prim_nt: int = 75
    wyklucz_3prim_nt: int = 50

    # uwagi do raportu
    uwagi: List[str] = field(default_factory=list)


PROFILE: Dict[str, HostProfile] = {

    'plant': HostProfile(
        nazwa='plant',
        opis='Angiosperms; antiviral pathway via AGO1/AGO2',
        gc_min=30.0, gc_max=52.0,
        guide_5_dozwolone='AU',
        guide_5_preferowane='U',
        passenger_5_dozwolone='GC',
        motywy_zabronione=(),
        motywy_raportowane=MOTYWY_IMMUNO_SSAKI,
        min_asymetria=0.3,
        max_mfe_guide=-3.0,
        dlugosci=(21, 22, 24),
        opis_dlugosci={
            21: 'DCL4 - PTGS, main antiviral pathway, mRNA cleavage',
            22: 'DCL2 - triggers transitivity (RDR6), secondary siRNA, systemic movement',
            24: 'DCL3 - RdDM, DNA methylation, transcriptional silencing',
        },
        uwagi=[
            'Mammalian immunostimulatory motifs are reported but not filtered '
            '- plants have neither TLR7/8 nor an interferon pathway.',
            'dsRNA-PTI is sequence-independent and therefore requires an '
            'experimental control (neutral dsRNA), not a filter.',
            'Checking the seed against host endogenous miRNAs is recommended '
            '(module filters.MirnaGuard).',
        ],
    ),

    'mammal': HostProfile(
        nazwa='mammal',
        opis='Mammals; AGO2, interferon pathway present',
        gc_min=30.0, gc_max=64.0,
        guide_5_dozwolone=None,
        guide_5_preferowane='AU',
        passenger_5_dozwolone=None,
        motywy_zabronione=MOTYWY_IMMUNO_SSAKI,
        motywy_raportowane=(),
        min_asymetria=0.3,
        max_mfe_guide=-3.0,
        dlugosci=(21,),
        opis_dlugosci={21: 'Canonical siRNA duplex with 2-nt 3-prime overhangs'},
        uwagi=[
            'Motifs UGUGU, GUCCUUCAA and UGGC are filtered as immunostimulatory '
            '(TLR7/8, interferon response).',
            'GC range wider than for plants, following Reynolds et al. 2004 '
            'and Ui-Tei et al. 2004.',
            'The 5-prime preference is treated as a scoring component rather '
            'than a hard filter - in mammals it follows mainly from '
            'thermodynamic asymmetry.',
        ],
    ),

    'insect': HostProfile(
        nazwa='insect',
        opis='Insects; Ago2 in the siRNA pathway, no interferons',
        gc_min=30.0, gc_max=60.0,
        guide_5_dozwolone=None,
        guide_5_preferowane='AU',
        passenger_5_dozwolone=None,
        motywy_zabronione=(),
        motywy_raportowane=MOTYWY_IMMUNO_SSAKI,
        min_asymetria=0.3,
        max_mfe_guide=-3.0,
        dlugosci=(21,),
        opis_dlugosci={21: 'Canonical siRNA duplex'},
        uwagi=[
            'Mammalian motifs reported, not filtered - no interferons.',
            'Induction of RNAi machinery gene expression by exogenous dsRNA and '
            'the effect of prior exposure have been described (Ye et al. 2019).',
        ],
    ),

    'generic': HostProfile(
        nazwa='generic',
        opis='Neutral profile - use when the host is undefined',
        gc_min=30.0, gc_max=55.0,
        guide_5_dozwolone=None,
        guide_5_preferowane='AU',
        passenger_5_dozwolone=None,
        motywy_zabronione=(),
        motywy_raportowane=MOTYWY_IMMUNO_SSAKI,
        min_asymetria=0.3,
        max_mfe_guide=-3.0,
        dlugosci=(21,),
        opis_dlugosci={21: 'Canonical siRNA duplex'},
        uwagi=[
            'The neutral profile applies shared rules only. Results should be '
            'treated as preliminary - use the appropriate profile for a '
            'specific host.',
        ],
    ),
}


def get_profile(nazwa: str) -> HostProfile:
    if nazwa not in PROFILE:
        raise ValueError(
            f'Nieznany profil: {nazwa}. Dostepne: {sorted(PROFILE)}')
    return PROFILE[nazwa]


def lista_profili() -> str:
    linie = ['AVAILABLE HOST PROFILES', '=' * 70]
    for n, p in PROFILE.items():
        linie.append(f'{n:<10} {p.opis}')
        linie.append(f'{"":10} GC {p.gc_min}-{p.gc_max}%, '
                     f'lengths {p.dlugosci}, '
                     f'motifs filtered: {len(p.motywy_zabronione)}')
    return '\n'.join(linie)


if __name__ == '__main__':
    print(lista_profili())
    print()
    p = get_profile('mammal')
    print(f'Profil {p.nazwa}: filtrowane motywy = {p.motywy_zabronione}')
    p = get_profile('plant')
    print(f'Profil {p.nazwa}: filtrowane motywy = {p.motywy_zabronione or "brak"}, '
          f'raportowane = {p.motywy_raportowane}')
