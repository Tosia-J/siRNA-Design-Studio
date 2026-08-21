"""
=============================================================================
thermo.py  --  TERMODYNAMIKA DUPLEKSU RNA/RNA
=============================================================================

STATUS: WLASNA IMPLEMENTACJA (zastepuje zewnetrzne kalkulatory)

Model: nearest-neighbour Watsona-Cricka dla RNA, parametry Xia i wsp. 1998
(Biochemistry 37:14719). To jest ten sam zestaw parametrow, ktorego uzywa
ViennaRNA do skladowej stackingowej -- roznica polega na tym, ze ViennaRNA
dodatkowo przeszukuje przestrzen struktur drugorzedowych algorytmem
dynamicznym (Zuker), a tu liczymy dupleks o ZNANEJ, pelnej komplementarnosci.

Dla dupleksu w pelni sparowanego przeszukiwanie jest zbedne, bo struktura
jest z gory zadana. Dlatego ten modul mozna napisac samodzielnie i wynik
bedzie identyczny co do parametrow -- w odroznieniu od faldowania pojedynczej
nici, gdzie wlasna implementacja bylaby gorsza (patrz fold.py).

CO LICZYMY I PO CO
------------------
1. dG37 calego dupleksu           -> ogolna sila wiazania siRNA z celem
2. dG 4 par zasad od konca 5' nici prowadzacej  -> ASYMETRIA
3. dG 4 par zasad od konca 5' nici pasazerskiej -> ASYMETRIA
4. Tm dupleksu

ASYMETRIA -- dlaczego to najwazniejsza liczba w calym module
------------------------------------------------------------
Kompleks RISC rozplata dupleks siRNA i zatrzymuje tylko jedna nic. Wybiera te,
ktorej koniec 5' jest SLABIEJ zwiazany -- bo od tego konca latwiej zaczac
rozplatanie. Jesli koniec 5' nici prowadzacej jest slabszy niz koniec 5' nici
pasazerskiej, RISC zaladuje wlasciwa nic. Jesli odwrotnie -- zaladuje nic
pasazerska, ktora nie jest komplementarna do celu, i siRNA nie zadziala, za to
moze wyciszyc cos przypadkowego.

Definiujemy:
    asymetria = dG5'(prowadzaca) - dG5'(pasazerska)

dG jest ujemne, wiec wartosc MNIEJ ujemna = slabsze wiazanie. Chcemy, zeby
koniec prowadzacej byl slabszy, czyli zeby asymetria byla DODATNIA.
Prog roboczy: > 0.3 kcal/mol. Powyzej 1.5 kcal/mol asymetria jest wyrazna.

Zrodlo reguly: Khvorova i wsp. 2003 (Cell 115:209), Schwarz i wsp. 2003
(Cell 115:199) -- dwie niezalezne prace, ktore odkryly to samo rownoczesnie.

Autor: Antonina Jarecka
=============================================================================
"""

from typing import Dict, Tuple

# --- Parametry Xia i wsp. 1998, dG37 w kcal/mol -----------------------------
# Klucz: dinukleotyd nici gornej czytany 5'->3'. Nic dolna jest domyslnie
# komplementarna. Np. 'GC' oznacza stack 5'-GC-3' / 3'-CG-5'.
# Pary rownowazne przez symetrie maja te sama wartosc (np. AA/UU == UU/AA).
NN_DG37: Dict[str, float] = {
    'AA': -0.93, 'UU': -0.93,
    'AU': -1.10,
    'UA': -1.33,
    'CU': -2.08, 'AG': -2.08,
    'CA': -2.11, 'UG': -2.11,
    'GU': -2.24, 'AC': -2.24,
    'GA': -2.35, 'UC': -2.35,
    'CG': -2.36,
    'GG': -3.26, 'CC': -3.26,
    'GC': -3.42,
}

# Entalpie (dH37, kcal/mol) -- potrzebne tylko do Tm
NN_DH: Dict[str, float] = {
    'AA': -6.82, 'UU': -6.82,
    'AU': -9.38,
    'UA': -7.69,
    'CU': -10.48, 'AG': -10.48,
    'CA': -10.44, 'UG': -10.44,
    'GU': -11.40, 'AC': -11.40,
    'GA': -12.44, 'UC': -12.44,
    'CG': -10.64,
    'GG': -13.39, 'CC': -13.39,
    'GC': -14.88,
}

INIT_DG = 4.09        # inicjacja helisy
INIT_DH = 3.61
TERM_AU_DG = 0.45     # kara za pare AU na koncu dupleksu
TERM_AU_DH = 3.72
R = 0.0019872         # stala gazowa, kcal/(mol*K)

COMP_RNA = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}


def to_rna(seq: str) -> str:
    return seq.upper().replace('T', 'U')


def revcomp_rna(seq: str) -> str:
    """Odwrotna komplementarnosc w alfabecie RNA."""
    s = to_rna(seq)
    return ''.join(COMP_RNA[c] for c in reversed(s))


def _validate(seq: str) -> str:
    s = to_rna(seq)
    bad = set(s) - set('ACGU')
    if bad:
        raise ValueError(f'Niedozwolone znaki w sekwencji: {sorted(bad)}')
    if len(s) < 2:
        raise ValueError('Sekwencja za krotka (min. 2 nt).')
    return s


def duplex_dg(seq: str, include_init: bool = True) -> float:
    """
    dG37 dupleksu w pelni sparowanego, gdzie `seq` to jedna z nici (5'->3').

    Sumujemy wklady wszystkich sasiadujacych par zasad, dodajemy inicjacje
    i kary za terminalne pary AU.
    """
    s = _validate(seq)
    dg = sum(NN_DG37[s[i:i + 2]] for i in range(len(s) - 1))
    if include_init:
        dg += INIT_DG
        if s[0] in 'AU':
            dg += TERM_AU_DG
        if s[-1] in 'AU':
            dg += TERM_AU_DG
    return dg


def duplex_dh(seq: str) -> float:
    """dH dupleksu -- pomocnicze do Tm."""
    s = _validate(seq)
    dh = sum(NN_DH[s[i:i + 2]] for i in range(len(s) - 1))
    dh += INIT_DH
    if s[0] in 'AU':
        dh += TERM_AU_DH
    if s[-1] in 'AU':
        dh += TERM_AU_DH
    return dh


def tm(seq: str, conc_M: float = 1e-6) -> float:
    """
    Temperatura topnienia dupleksu w stopniach Celsjusza.

    conc_M: stezenie calkowite nici (domyslnie 1 uM -- typowe dla oznaczen
    in vitro). Dla dupleksu niesamokomplementarnego uzywamy Ct/4.
    """
    s = _validate(seq)
    dh = duplex_dh(s)
    dg = duplex_dg(s)
    # dG = dH - T*dS  =>  dS = (dH - dG) / T  przy T = 310.15 K
    ds = (dh - dg) / 310.15
    t_kelvin = dh / (ds + R * _math_log(conc_M / 4.0))
    return t_kelvin - 273.15


def _math_log(x: float) -> float:
    import math
    return math.log(x)


def end_dg(seq_5to3: str, n_bp: int = 4) -> float:
    """
    dG pierwszych `n_bp` par zasad liczac od konca 5' podanej nici.

    UWAGA: celowo BEZ inicjacji i bez kar terminalnych. Porownujemy tu dwa
    konce tego samego dupleksu, wiec czlony stale i tak by sie skrocily,
    a ich pominiecie sprawia, ze liczba jest czysta suma stackingu.
    """
    s = _validate(seq_5to3)
    n = min(n_bp, len(s) - 1)
    return sum(NN_DG37[s[i:i + 2]] for i in range(n))


def asymmetry(guide_5to3: str, n_bp: int = 4) -> Tuple[float, float, float]:
    """
    Zwraca (dG5'_prowadzaca, dG5'_pasazerska, asymetria).

    asymetria > 0  ->  koniec 5' nici prowadzacej jest SLABSZY  ->  DOBRZE
    asymetria < 0  ->  RISC zaladuje nic pasazerska            ->  ZLE

    Argumentem jest nic PROWADZACA (antysensowna). Nic pasazerska
    wyliczamy jako jej odwrotna komplementarnosc.
    """
    guide = _validate(guide_5to3)
    passenger = revcomp_rna(guide)
    dg_g = end_dg(guide, n_bp)
    dg_p = end_dg(passenger, n_bp)
    return dg_g, dg_p, dg_g - dg_p


def gc_content(seq: str) -> float:
    s = _validate(seq)
    return 100.0 * (s.count('G') + s.count('C')) / len(s)


def has_homopolymer(seq: str, n: int = 4) -> bool:
    """Ciag n identycznych nukleotydow -- utrudnia synteze i sprzyja
    niespecyficznemu wiazaniu."""
    s = to_rna(seq)
    return any(b * n in s for b in 'ACGU')


if __name__ == '__main__':
    # Test na kandydacie A-21 z projektu anty-GFP
    g = 'UUGAAGUUCACCUUGAUGCCG'
    dg_g, dg_p, asym = asymmetry(g)
    print(f'nic prowadzaca : 5-{g}-3')
    print(f'nic pasazerska : 5-{revcomp_rna(g)}-3')
    print(f'GC             : {gc_content(g):.1f} %')
    print(f'dG dupleksu    : {duplex_dg(g):.2f} kcal/mol')
    print(f'Tm             : {tm(g):.1f} C')
    print(f'dG5 prowadzaca : {dg_g:.2f}')
    print(f'dG5 pasazerska : {dg_p:.2f}')
    print(f'ASYMETRIA      : {asym:+.2f}  ->  '
          f'{"OK, zaladuje sie wlasciwa nic" if asym > 0.3 else "ZLE"}')
