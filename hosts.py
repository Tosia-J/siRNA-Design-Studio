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

-----------------------------------------------------------------------------
CO SIE ZMIENILO W TEJ WERSJI (rozdzial profilu roslinnego)
-----------------------------------------------------------------------------

Dotychczasowy profil 'plant' mial jeden zestaw regul dla calej rosliny.
Jest to uproszczenie, ktore zaciera roznice miedzy dwoma odrebnymi szlakami.

Wang i wsp. 2023 (Genes Dev 37:103-118) pokazali w ukladzie rekonstytuowanym,
ze AGO4 - efektor szlaku RdDM - NIE MA ani wymagania dlugosci, ani preferencji
nukleotydu 5'. Nici prowadzace 21, 22, 23 i 24 nt programuja cieciezel z
porownywalna wydajnoscia, a A, U, C i G na koncu 5' dzialaja tak samo dobrze.
Przewaga 24 nt i bias 5'-A obserwowane in vivo pochodza od enzymow
wytwarzajacych siRNA - inicjacji Pol IV, RDR2 i kodu tnacego DCL3 - a nie od
samego AGO4.

Wynika stad, ze regula 'wymagaj 5'-U albo 5'-A' jest sluszna dla PTGS
(AGO1/AGO2, Mi i wsp. 2008), ale NIEUPRAWNIONA dla RdDM. Stosowanie jej do
kandydatow celujacych w chromatyne odrzuca poprawne sekwencje bez podstawy
biologicznej.

Stad dwa profile roslinne:

    plant_ptgs   AGO1/AGO2, cieciezel mRNA, wymog 5'-U/A
    plant_rddm   AGO4/AGO6, metylacja DNA, brak wymogu 5'

Klucz 'plant' pozostaje i wskazuje na plant_ptgs, zeby nie zepsuc istniejacych
wywolan.

-----------------------------------------------------------------------------
TRYB DOSTARCZANIA - DRUGA OS, NIEZALEZNA OD GOSPODARZA
-----------------------------------------------------------------------------

Dopuszczalne dlugosci zaleza nie tylko od organizmu, ale i od tego, JAK
czasteczka trafia do komorki. To dwie rozne sytuacje:

    'transgenic'  Roslina dostaje dluga spinke albo dsRNA i tnie ja sama.
                  Dlugosc produktu wybiera Dicer, nie projektant.
                  DCL4 -> 21, DCL2 -> 22, DCL3 -> 24.
                  23 nt NIE JEST osiagalne: DCL3 wytwarza 23-mer wylacznie
                  jako nic pasazerska dupleksu 24/23, a AGO4 tnie ja na
                  fragmenty 11 i 12 nt (Wang i wsp. 2023).

    'synthetic'   Kompleks dostaje gotowy dupleks - SIGS, nosnik, transfekcja.
                  Etap dicingu jest pominiety, wiec kod tnacy DCL3 nie
                  obowiazuje. AGO przyjmuje 21-24 nt.
                  23 nt STAJE SIE dopuszczalne i testuje dlugosc, ktora w
                  naturze nigdy nie wystepuje jako nic prowadzaca.

O tym, ktora nic zostanie zaladowana, decyduje asymetria nawisow 3'.
W dupleksie 24/23 nic 23-nt ma nawis 1 nt, a nic 24-nt nawis 2 nt; to
prawdopodobnie ten 2-nt nawis znakuje nic przeznaczona do zaladowania.
Przy projektowaniu dupleksu syntetycznego nalezy wiec dac nici prowadzacej
nawis 2 nt, a pasazerskiej 1 nt.

JAK UZYWAC

    from hosts import get_profile, dozwolone_dlugosci

    profil = get_profile('plant_rddm')
    dlugosci = dozwolone_dlugosci(profil, tryb='synthetic')   # (21, 22, 23, 24)
    dlugosci = dozwolone_dlugosci(profil, tryb='transgenic')  # (24,)

Profil przekazywany jest do filters.physicochemical_filters() i steruje
wlaczaniem poszczegolnych regul.

PISMIENNICTWO DODANE W TEJ WERSJI
    Wang F., Huang H.-Y., Huang J., Singh J., Pikaard C.S. 2023.
    Enzymatic reactions of AGO4 in RNA-directed DNA methylation.
    Genes & Development 37:103-118. doi:10.1101/gad.350240.122
    McCue A.D. i wsp. 2015. ARGONAUTE 6 bridges transposable element
    mRNA-derived siRNAs to the establishment of DNA methylation.
    EMBO J 34:20-35.

Autor: Antonina Jarecka
=============================================================================
"""

from dataclasses import dataclass, field, replace
from typing import Dict, List, Optional, Tuple


# Motywy immunostymulacyjne u ssakow
MOTYWY_IMMUNO_SSAKI: Tuple[str, ...] = ('UGUGU', 'GUCCUUCAA', 'UGGC')

# Tryby dostarczania czasteczki do komorki
TRYBY: Tuple[str, ...] = ('transgenic', 'synthetic')


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

    # --- pola dodane przy rozdziale profilu roslinnego -------------------
    # szlak wyciszania, ktorego dotycza reguly profilu
    sciezka: str = 'PTGS'

    # bialka Argonaute niosace ten szlak - do raportu i do uzasadnien
    argonaute: Tuple[str, ...] = ()

    # dlugosci osiagalne wylacznie przy dostarczeniu gotowego dupleksu,
    # bo etap dicingu jest wtedy pominiety
    dlugosci_tylko_syntetyczne: Tuple[int, ...] = ()

    # asymetria nawisow 3' decydujaca o wyborze nici (w nukleotydach)
    nawis_guide_3prim: int = 2
    nawis_passenger_3prim: int = 2

    # czy wymog 5'-konca ma podstawe w selektywnosci kieszeni MID
    wymog_5prim_uzasadniony: bool = True


PROFILE: Dict[str, HostProfile] = {

    'plant_ptgs': HostProfile(
        nazwa='plant_ptgs',
        opis='Angiosperms, post-transcriptional silencing; antiviral pathway '
             'via AGO1/AGO2',
        gc_min=30.0, gc_max=52.0,
        guide_5_dozwolone='AU',
        guide_5_preferowane='U',
        passenger_5_dozwolone='GC',
        motywy_zabronione=(),
        motywy_raportowane=MOTYWY_IMMUNO_SSAKI,
        min_asymetria=0.3,
        max_mfe_guide=-3.0,
        dlugosci=(21, 22),
        opis_dlugosci={
            21: 'DCL4 - PTGS, main antiviral pathway, mRNA cleavage',
            22: 'DCL2 - triggers transitivity (RDR6), secondary siRNA, '
                'systemic movement',
        },
        sciezka='PTGS',
        argonaute=('AGO1', 'AGO2'),
        dlugosci_tylko_syntetyczne=(),
        nawis_guide_3prim=2,
        nawis_passenger_3prim=2,
        wymog_5prim_uzasadniony=True,
        uwagi=[
            'The 5-prime requirement rests on MID-pocket selectivity of plant '
            'AGO1 (5-prime U) and AGO2 (5-prime A), Mi et al. 2008.',
            'Mammalian immunostimulatory motifs are reported but not filtered '
            '- plants have neither TLR7/8 nor an interferon pathway.',
            'dsRNA-PTI is sequence-independent and therefore requires an '
            'experimental control (neutral dsRNA), not a filter.',
            'Checking the seed against host endogenous miRNAs is recommended '
            '(module filters.MirnaGuard).',
        ],
    ),

    'plant_rddm': HostProfile(
        nazwa='plant_rddm',
        opis='Angiosperms, RNA-directed DNA methylation; AGO4/AGO6, '
             'transcriptional silencing',
        gc_min=30.0, gc_max=52.0,
        guide_5_dozwolone=None,          # AGO4 nie selekcjonuje po 5' koncu
        guide_5_preferowane='A',         # bias in vivo, ale tylko do scoringu
        passenger_5_dozwolone=None,
        motywy_zabronione=(),
        motywy_raportowane=MOTYWY_IMMUNO_SSAKI,
        min_asymetria=0.3,
        max_mfe_guide=-3.0,
        dlugosci=(24,),
        opis_dlugosci={
            21: 'AGO6, non-canonical RdDM from TE-derived siRNA '
                '(McCue et al. 2015) - synthetic delivery only',
            22: 'AGO6, non-canonical RdDM - synthetic delivery only',
            23: 'Never a natural guide: DCL3 produces 23-mers only as the '
                'passenger strand of a 24/23 duplex, sliced by AGO4 into '
                '11- and 12-nt fragments (Wang et al. 2023). Reachable only '
                'as a pre-formed synthetic duplex.',
            24: 'DCL3 - canonical RdDM, de novo cytosine methylation via DRM2',
        },
        sciezka='RdDM',
        argonaute=('AGO4', 'AGO6'),
        dlugosci_tylko_syntetyczne=(21, 22, 23),
        nawis_guide_3prim=2,
        nawis_passenger_3prim=1,
        wymog_5prim_uzasadniony=False,
        uwagi=[
            'No 5-prime requirement is applied. AGO4 loads guides beginning '
            'with A, U, C or G with comparable efficiency; the 5-prime A bias '
            'observed in vivo comes from Pol IV initiation, RDR2 and the DCL3 '
            'dicing code, not from AGO4 (Wang et al. 2023).',
            'No length requirement either: guides of 21 to 24 nt program '
            'slicing equally well in the reconstituted system.',
            'Strand choice is governed by 3-prime overhang asymmetry. Give '
            'the guide a 2-nt overhang and the passenger a 1-nt overhang when '
            'ordering a synthetic duplex.',
            'A 5-prime monophosphate raises loading and slicing five- to '
            'eightfold over a hydroxyl or triphosphate end.',
            'Guide 24 nt paired with a 12-nt passenger fragment still slices, '
            'so a partially occupied duplex is not necessarily inactive.',
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
        sciezka='PTGS',
        argonaute=('AGO2',),
        nawis_guide_3prim=2,
        nawis_passenger_3prim=2,
        wymog_5prim_uzasadniony=False,
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
        sciezka='PTGS',
        argonaute=('Ago2',),
        wymog_5prim_uzasadniony=False,
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
        sciezka='unspecified',
        argonaute=(),
        wymog_5prim_uzasadniony=False,
        uwagi=[
            'The neutral profile applies shared rules only. Results should be '
            'treated as preliminary - use the appropriate profile for a '
            'specific host.',
        ],
    ),
}


# Zgodnosc wsteczna: 'plant' nadal dziala i wskazuje na szlak przeciwwirusowy,
# bo taki byl sens tego klucza w poprzedniej wersji.
ALIASY: Dict[str, str] = {
    'plant': 'plant_ptgs',
    'plant_ago4': 'plant_rddm',
    'rddm': 'plant_rddm',
    'ptgs': 'plant_ptgs',
}


def get_profile(nazwa: str) -> HostProfile:
    klucz = ALIASY.get(nazwa, nazwa)
    if klucz not in PROFILE:
        dostepne = sorted(set(PROFILE) | set(ALIASY))
        raise ValueError(f'Nieznany profil: {nazwa}. Dostepne: {dostepne}')
    return PROFILE[klucz]


def dozwolone_dlugosci(profil: HostProfile, tryb: str = 'transgenic') -> Tuple[int, ...]:
    """Dlugosci osiagalne dla danego profilu przy danym sposobie dostarczenia.

    'transgenic' - roslina tnie prekursor sama, wiec dlugosc wybiera Dicer.
    'synthetic'  - kompleks dostaje gotowy dupleks, wiec etap dicingu odpada
                   i dopuszczalne staja sie takze dlugosci nieosiagalne
                   biogenetycznie (m.in. 23 nt w RdDM).
    """
    if tryb not in TRYBY:
        raise ValueError(f'Nieznany tryb: {tryb}. Dostepne: {list(TRYBY)}')
    if tryb == 'transgenic':
        return tuple(sorted(profil.dlugosci))
    return tuple(sorted(set(profil.dlugosci) | set(profil.dlugosci_tylko_syntetyczne)))


def profil_dla_trybu(profil: HostProfile, tryb: str = 'transgenic') -> HostProfile:
    """Kopia profilu z dlugosciami juz rozstrzygnietymi dla danego trybu.

    Wygodne, gdy design.py oczekuje jednego obiektu i nie ma sie dowiadywac
    o istnieniu trybow dostarczania.
    """
    return replace(profil, dlugosci=dozwolone_dlugosci(profil, tryb))


def opis_nawisow(profil: HostProfile) -> str:
    """Jednozdaniowa instrukcja do zamowienia dupleksu syntetycznego."""
    return (f'Guide strand: {profil.nawis_guide_3prim}-nt 3-prime overhang. '
            f'Passenger strand: {profil.nawis_passenger_3prim}-nt 3-prime '
            f'overhang. The asymmetry marks the strand intended for loading.')


def lista_profili() -> str:
    linie = ['AVAILABLE HOST PROFILES', '=' * 70]
    for n, p in PROFILE.items():
        linie.append(f'{n:<12} {p.opis}')
        linie.append(f'{"":12} pathway {p.sciezka}, '
                     f'{"/".join(p.argonaute) or "-"}, '
                     f'GC {p.gc_min}-{p.gc_max}%')
        linie.append(f'{"":12} lengths transgenic {dozwolone_dlugosci(p, "transgenic")}, '
                     f'synthetic {dozwolone_dlugosci(p, "synthetic")}, '
                     f'motifs filtered: {len(p.motywy_zabronione)}')
    if ALIASY:
        linie.append('-' * 70)
        linie.append('aliases: ' + ', '.join(f'{a} -> {t}' for a, t in ALIASY.items()))
    return '\n'.join(linie)


if __name__ == '__main__':
    print(lista_profili())
    print()

    p = get_profile('mammal')
    print(f'Profil {p.nazwa}: filtrowane motywy = {p.motywy_zabronione}')

    p = get_profile('plant')
    print(f'Profil {p.nazwa}: filtrowane motywy = {p.motywy_zabronione or "brak"}, '
          f'raportowane = {p.motywy_raportowane}')
    print(f'  wymog 5-prim: {p.guide_5_dozwolone} '
          f'(uzasadniony: {p.wymog_5prim_uzasadniony})')

    p = get_profile('plant_rddm')
    print(f'Profil {p.nazwa}: wymog 5-prim = {p.guide_5_dozwolone or "brak"} '
          f'(uzasadniony: {p.wymog_5prim_uzasadniony})')
    print(f'  dlugosci transgeniczne: {dozwolone_dlugosci(p, "transgenic")}')
    print(f'  dlugosci syntetyczne:   {dozwolone_dlugosci(p, "synthetic")}')
    print(f'  {opis_nawisow(p)}')
    print()
    print('23 nt:', p.opis_dlugosci[23])
