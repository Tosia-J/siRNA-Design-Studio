"""
=============================================================================
blast_criteria.py  --  KRYTERIA BLAST Z PUNKTACJA WAZONA
=============================================================================

PO CO TEN MODUL

Pierwotny pipeline oceny off-target opieral sie na siedmiu parametrach
zwracanych przez blastn, z recznie przypisanymi wagami:

    Mismatches         5/5   xnorm
    % Identity         5/5   xnorm  xinv
    E-value            4/5
    Alignment Length   4/5   xnorm
    Hit Count          4/5   xinv  xinv
    Bit Score          3/5
    Gap Opens          2/5   xinv

Zastapienie BLAST-a indeksem k-merowym (modul offtarget) usunelo te warstwe
z pipeline'u. Niniejszy modul ja przywraca, poniewaz:

  1. Jest to wynik juz uzyskany i opisany - usuniecie go bez zastapienia
     rownowaznikiem oznaczaloby utrate porownywalnosci z wczesniejszym etapem
     projektu.
  2. Parametry BLAST sa standardem w literaturze, wiec ich raportowanie
     ulatwia porownanie z innymi narzedziami.
  3. Rownolegle zastosowanie dwoch niezaleznych metod oceny off-target jest
     mocniejsze niz jedna.

Modul NIE zastepuje offtarget.py. Obie warstwy raportuje sie osobno.

=============================================================================
OCENA KRYTYCZNA ORYGINALNEGO ZESTAWU PARAMETROW
=============================================================================

Zestaw wymaga trzech korekt. Podaje je jawnie, poniewaz sa to zarzuty,
ktore recenzent postawi w pierwszej kolejnosci.

--- PROBLEM 1: WSPOLLINIOWOSC ALGEBRAICZNA ---

Trzy parametry sa ze soba powiazane rownaniem:

    % Identity = (Alignment Length - Mismatches - Gaps) / Alignment Length

Znajac dwa z nich, trzeci jest wyznaczony. Przypisanie im niezaleznych wag
(5/5, 5/5, 4/5) oznacza, ze ta sama informacja wchodzi do wyniku trzykrotnie,
z laczna waga 14/15 - czyli dominuje caly scoring.

--- PROBLEM 2: E-VALUE I BIT SCORE MIERZA TO SAMO ---

E-value wyprowadza sie z bit score wzorem:

    E = m * n * 2^(-S')

gdzie S' to bit score, m to dlugosc zapytania, n to rozmiar bazy. Dla
ustalonej bazy i ustalonej dlugosci zapytania E-value jest scisle malejaca
funkcja bit score. Sa to zatem dwa zapisy tej samej wielkosci, a nie dwa
niezalezne kryteria. Laczna waga 7/15 przypisana jednej informacji.

--- PROBLEM 3: GAP OPENS MA ZEROWA WARIANCJE ---

Przy zapytaniu dlugosci 21 nt i dopasowaniach o wysokiej identycznosci liczba
przerw wynosi praktycznie zawsze zero. Parametr o zerowej wariancji nie wnosi
informacji do rankingu, niezaleznie od przypisanej wagi.

--- PODSUMOWANIE ---

Z siedmiu parametrow niezaleznych wymiarow jest okolo trzech:

    (a) podobienstwo sekwencji   -- Mismatches / % Identity / Align. Length
    (b) istotnosc statystyczna   -- E-value / Bit Score
    (c) liczba trafien           -- Hit Count

Hit Count jest jedynym parametrem calkowicie niezaleznym od pozostalych
i prawdopodobnie najbardziej wartosciowym - mierzy, w ilu miejscach genomu
sekwencja moze zadzialac, a nie jak dobrze pasuje do jednego z nich.

--- OGRANICZENIE POWAZNIEJSZE NIZ WSPOLLINIOWOSC ---

Wszystkie siedem parametrow jest SLEPYCH NA POZYCJE niedopasowania.
Dopasowanie z jednym niedopasowaniem w pozycji 3 (region seed) i dopasowanie
z jednym niedopasowaniem w pozycji 19 otrzymuja identyczna punktacje, mimo ze
pierwsze jest biologicznie znacznie grozniejsze.

Dodatkowo blastn z domyslnym word_size = 11 w ogole nie zglasza dopasowan
regionu seed (7-8 nt), poniewaz sa ponizej progu zasiewu.

Z tego powodu modul offtarget.py pozostaje warstwa podstawowa, a niniejszy
modul warstwa uzupelniajaca i porownawcza.

=============================================================================

Autor: Antonina Jarecka
=============================================================================
"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import math
import csv


# ============================================================================
# DEFINICJA PARAMETROW
# ============================================================================

@dataclass
class Kryterium:
    nazwa: str
    waga_surowa: int          # w skali 1-5, jak w oryginalnym zestawieniu
    normalizuj: bool          # xnorm - normalizacja min-max do [0,1]
    odwroc: bool              # xinv - mniejsza wartosc = lepiej
    opis: str
    uwaga: str = ''


KRYTERIA: Dict[str, Kryterium] = {
    'mismatches': Kryterium(
        'Mismatches', 5, True, False,
        'Liczba niedopasowanych nukleotydow w dopasowaniu',
        'Wiecej niedopasowan = bezpieczniej, wiec BEZ odwrocenia.'),
    'pct_identity': Kryterium(
        '% Identity', 5, True, True,
        'Odsetek identycznych nukleotydow miedzy nicia antysensowna '
        'a genomem gospodarza',
        'Wspolliniowe z Mismatches i Alignment Length.'),
    'evalue': Kryterium(
        'E-value', 4, True, False,
        'Prawdopodobienstwo statystyczne niespecyficznego dopasowania',
        'Wyprowadzone z Bit Score - nie jest niezaleznym kryterium. '
        'Stosowac po transformacji logarytmicznej.'),
    'align_length': Kryterium(
        'Alignment Length', 4, True, True,
        'Dlugosc dopasowanego fragmentu',
        'Wspolliniowe z % Identity.'),
    'hit_count': Kryterium(
        'Hit Count', 4, True, True,
        'Laczna liczba trafien nici antysensownej w genomie gospodarza',
        'Jedyny parametr calkowicie niezalezny od pozostalych.'),
    'bit_score': Kryterium(
        'Bit Score', 3, True, True,
        'Wartosc liczbowa istotnosci dopasowania przed wazeniem',
        'Wspolliniowy z E-value.'),
    'gap_opens': Kryterium(
        'Gap Opens', 2, True, False,
        'Liczba przerw w dopasowaniu',
        'Przy zapytaniu 21 nt praktycznie zawsze zero - zerowa wariancja.'),
}

# Zestaw zredukowany - najmniejszy podzbior oryginalnej siodemki
# zachowujacy trzy niezalezne wymiary informacji
KRYTERIA_ZREDUKOWANE = ('mismatches', 'evalue', 'hit_count')


# ============================================================================
# ZESTAW PRZEPROJEKTOWANY
# ============================================================================
#
# Zestaw zredukowany usuwa redundancje, ale nadal opiera sie na parametrach
# BLAST, ktore powstaly do porownywania dlugich sekwencji. Ponizszy zestaw
# jest zaprojektowany od nowa pod konkretne pytanie biologiczne:
#
#     "W ilu miejscach genomu gospodarza ta sekwencja moze zadzialac
#      i jak grozne jest najgorsze z tych miejsc?"
#
# Sa to dwa niezalezne wymiary - liczba zagrozen i ich dotkliwosc - i zaden
# z nich nie da sie wyprowadzic z drugiego.
#
# ROZNICA WZGLEDEM ZESTAWU ORYGINALNEGO
#
#   oryginalny            : 7 parametrow, ~3 niezalezne wymiary,
#                           21 z 27 jednostek wagi na kryteria redundantne
#   przeprojektowany      : 3 parametry, 3 niezalezne wymiary,
#                           kazdy o jasnej interpretacji biologicznej
#
# UWAGA: zestaw ten nadal jest SLEPY NA POZYCJE niedopasowania. Rozroznienie
# miedzy niedopasowaniem w regionie seed a w regionie 3' wymaga wyjscia poza
# tabelaryczny format BLAST i realizowane jest w module offtarget.py.

KRYTERIA_PRZEPROJEKTOWANE: Dict[str, Kryterium] = {
    'n_trafien_wysokiej_identycznosci': Kryterium(
        'Trafienia >=85% identycznosci', 5, True, True,
        'Liczba dopasowan o identycznosci co najmniej 85% na odcinku '
        'co najmniej 16 nt',
        'Odpowiada na pytanie: ile miejsc w genomie moze zostac przeciete. '
        'Prog 85% na 16 nt odpowiada okolo 2-3 niedopasowaniom - powyzej '
        'tego ciecie przez RISC jest malo prawdopodobne.'),
    'identycznosc_najgorszego': Kryterium(
        'Identycznosc najgorszego trafienia', 5, True, True,
        'Najwyzsza identycznosc sposrod wszystkich trafien poza celem',
        'Odpowiada na pytanie: jak grozne jest najgorsze miejsce. '
        'Niezalezne od liczby trafien.'),
    'pokrycie_seed': Kryterium(
        'Trafienia z pelnym seedem', 4, True, True,
        'Liczba dopasowan obejmujacych pozycje 2-8 nici prowadzacej '
        'bez niedopasowan',
        'Jedyny parametr wrazliwy na pozycje mozliwy do wyliczenia '
        'z wyniku BLAST - wymaga kolumn qstart i qend.'),
}


def przelicz_kryteria_przeprojektowane(
        trafienia: List[Dict],
        dlugosc_zapytania: int = 21,
        prog_identycznosci: float = 85.0,
        prog_dlugosci: int = 16) -> Dict[str, float]:
    """
    Oblicza zestaw przeprojektowany z surowego wyniku BLAST.

    Wymaga kolumn: pident, length, qstart, qend (outfmt 6).
    """
    if not trafienia:
        return {'n_trafien_wysokiej_identycznosci': 0.0,
                'identycznosc_najgorszego': 0.0,
                'pokrycie_seed': 0.0}

    wysokie = [t for t in trafienia
               if t['pident'] >= prog_identycznosci
               and t['length'] >= prog_dlugosci]

    # region seed nici prowadzacej to pozycje 2-8 zapytania
    # dopasowanie musi je obejmowac w calosci
    z_seedem = [t for t in trafienia
                if t.get('qstart', 99) <= 2 and t.get('qend', 0) >= 8]

    return {
        'n_trafien_wysokiej_identycznosci': float(len(wysokie)),
        'identycznosc_najgorszego': max(t['pident'] for t in trafienia),
        'pokrycie_seed': float(len(z_seedem)),
    }


# ============================================================================
# WCZYTANIE WYNIKOW BLAST
# ============================================================================

# blastn -outfmt 6 zwraca kolumny:
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
KOLUMNY_OUTFMT6 = ('qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                   'gapopen', 'qstart', 'qend', 'sstart', 'send',
                   'evalue', 'bitscore')


def wczytaj_blast_outfmt6(sciezka: str) -> Dict[str, List[Dict]]:
    """
    Wczytuje wynik `blastn -outfmt 6` i grupuje trafienia wg identyfikatora
    zapytania.

    Zalecane wywolanie BLAST dla krotkich zapytan:

        blastn -task blastn-short -word_size 7 -evalue 1000 \\
               -query kandydaci.fasta -db transkryptom \\
               -outfmt 6 -out wynik.tsv

    Uwaga: parametry domyslne (word_size 11, evalue 10) pomijaja wiekszosc
    biologicznie istotnych dopasowan dla zapytan 21-nt.
    """
    wynik: Dict[str, List[Dict]] = {}
    with open(sciezka) as fh:
        for linia in fh:
            linia = linia.strip()
            if not linia or linia.startswith('#'):
                continue
            pola = linia.split('\t')
            if len(pola) < len(KOLUMNY_OUTFMT6):
                continue
            rek = dict(zip(KOLUMNY_OUTFMT6, pola))
            for k in ('pident', 'evalue', 'bitscore'):
                rek[k] = float(rek[k])
            for k in ('length', 'mismatch', 'gapopen'):
                rek[k] = int(rek[k])
            wynik.setdefault(rek['qseqid'], []).append(rek)
    return wynik


def agreguj_trafienia(trafienia: List[Dict]) -> Dict[str, float]:
    """
    Sprowadza liste trafien dla jednego kandydata do zestawu siedmiu
    parametrow.

    Zasada: dla parametrow opisujacych pojedyncze dopasowanie bierzemy
    NAJGORSZY przypadek, czyli trafienie o najnizszym E-value (najbardziej
    istotne). Hit Count zlicza wszystkie trafienia.
    """
    if not trafienia:
        # brak trafien - wartosci skrajnie korzystne
        return {'mismatches': 21.0, 'pct_identity': 0.0, 'evalue': 1000.0,
                'align_length': 0.0, 'hit_count': 0.0, 'bit_score': 0.0,
                'gap_opens': 0.0}

    najgorsze = min(trafienia, key=lambda t: t['evalue'])
    return {
        'mismatches': float(najgorsze['mismatch']),
        'pct_identity': float(najgorsze['pident']),
        'evalue': float(najgorsze['evalue']),
        'align_length': float(najgorsze['length']),
        'hit_count': float(len(trafienia)),
        'bit_score': float(najgorsze['bitscore']),
        'gap_opens': float(najgorsze['gapopen']),
    }


# ============================================================================
# PUNKTACJA
# ============================================================================

def _minmax(wartosci: List[float]) -> List[float]:
    lo, hi = min(wartosci), max(wartosci)
    if hi - lo < 1e-12:
        return [0.5] * len(wartosci)
    return [(w - lo) / (hi - lo) for w in wartosci]


def punktuj(kandydaci: List[Dict],
            zestaw: Optional[Tuple[str, ...]] = None,
            log_evalue: bool = True) -> List[Dict]:
    """
    Oblicza wynik wazony wedlug kryteriow BLAST.

    kandydaci : lista slownikow zawierajacych klucze z KRYTERIA
                (np. z agreguj_trafienia)
    zestaw    : ktore kryteria uwzglednic; None = wszystkie siedem
    log_evalue: czy zastosowac transformacje -log10 do E-value

    Zwraca kopie listy z polami norm_<kryterium> oraz wynik_blast.
    """
    if not kandydaci:
        return []

    # obsluga obu zestawow: oryginalnego i przeprojektowanego
    wszystkie = {**KRYTERIA, **KRYTERIA_PRZEPROJEKTOWANE}
    uzyte = zestaw if zestaw else tuple(KRYTERIA)
    suma_wag = sum(wszystkie[k].waga_surowa for k in uzyte)

    out = [dict(c) for c in kandydaci]

    for klucz in uzyte:
        kr = wszystkie[klucz]
        surowe = [float(c[klucz]) for c in out]

        if klucz == 'evalue' and log_evalue:
            # -log10(E) rosnie, gdy E maleje; dodajemy male epsilon
            surowe = [-math.log10(max(v, 1e-180)) for v in surowe]
            # po tej transformacji wieksza wartosc = bardziej istotne
            # trafienie = GORZEJ, wiec odwracamy
            surowe = [-v for v in surowe]
        elif kr.odwroc:
            surowe = [-v for v in surowe]

        znorm = _minmax(surowe) if kr.normalizuj else surowe
        for c, v in zip(out, znorm):
            c[f'norm_{klucz}'] = round(v, 4)

    for c in out:
        c['wynik_blast'] = round(
            sum(c[f'norm_{k}'] * wszystkie[k].waga_surowa / suma_wag
                for k in uzyte), 4)

    out.sort(key=lambda c: -c['wynik_blast'])
    for i, c in enumerate(out, 1):
        c['ranga_blast'] = i
    return out


# ============================================================================
# DIAGNOSTYKA WSPOLLINIOWOSCI
# ============================================================================

def korelacja_pearsona(x: List[float], y: List[float]) -> float:
    n = len(x)
    if n < 2:
        return 0.0
    sx, sy = sum(x) / n, sum(y) / n
    licz = sum((a - sx) * (b - sy) for a, b in zip(x, y))
    m1 = math.sqrt(sum((a - sx) ** 2 for a in x))
    m2 = math.sqrt(sum((b - sy) ** 2 for b in y))
    if m1 * m2 < 1e-12:
        return 0.0
    return licz / (m1 * m2)


def diagnostyka(kandydaci: List[Dict], prog: float = 0.85) -> Dict:
    """
    Sprawdza wspolliniowosc miedzy kryteriami na rzeczywistych danych.

    Zwraca pary o korelacji przekraczajacej prog oraz kryteria o zerowej
    wariancji. Wynik nalezy zamiescic w czesci metodycznej - jawne wykazanie
    wspolliniowosci i jej uwzglednienie jest mocniejsze niz przemilczenie.
    """
    klucze = [k for k in KRYTERIA if k in kandydaci[0]]
    dane = {k: [float(c[k]) for c in kandydaci] for k in klucze}

    zerowa_wariancja = [k for k in klucze
                        if max(dane[k]) - min(dane[k]) < 1e-12]

    pary = []
    for i, a in enumerate(klucze):
        for b in klucze[i + 1:]:
            r = korelacja_pearsona(dane[a], dane[b])
            if abs(r) >= prog:
                pary.append({'a': KRYTERIA[a].nazwa, 'b': KRYTERIA[b].nazwa,
                             'r': round(r, 3)})

    laczna_waga_redundantna = 0
    widziane = set()
    for p in pary:
        for nazwa in (p['a'], p['b']):
            klucz = next(k for k, v in KRYTERIA.items() if v.nazwa == nazwa)
            if klucz not in widziane:
                widziane.add(klucz)
                laczna_waga_redundantna += KRYTERIA[klucz].waga_surowa

    return {
        'pary_wspolliniowe': pary,
        'kryteria_zerowej_wariancji': [KRYTERIA[k].nazwa
                                       for k in zerowa_wariancja],
        'laczna_waga_kryteriow_redundantnych': laczna_waga_redundantna,
        'suma_wszystkich_wag': sum(KRYTERIA[k].waga_surowa for k in klucze),
        'zalecenie': (
            'Rozwaz zestaw zredukowany: '
            + ', '.join(KRYTERIA[k].nazwa for k in KRYTERIA_ZREDUKOWANE)
            if pary else 'Brak istotnej wspolliniowosci w tym zbiorze.'),
    }


def raport(kandydaci: List[Dict], diag: Dict, top_n: int = 10) -> str:
    L = ['PUNKTACJA WEDLUG KRYTERIOW BLAST', '=' * 78]
    L.append(f'{"#":<4}{"kandydat":<20}{"wynik":>8}   '
             + ' '.join(f'{KRYTERIA[k].nazwa[:8]:>9}' for k in KRYTERIA
                        if f'norm_{k}' in kandydaci[0]))
    L.append('-' * 78)
    for c in kandydaci[:top_n]:
        L.append(f'{c["ranga_blast"]:<4}{str(c.get("nazwa", "?"))[:20]:<20}'
                 f'{c["wynik_blast"]:>8.4f}   '
                 + ' '.join(f'{c[f"norm_{k}"]:>9.3f}' for k in KRYTERIA
                            if f'norm_{k}' in c))
    L.append('')
    L.append('DIAGNOSTYKA WSPOLLINIOWOSCI')
    L.append('-' * 78)
    if diag['pary_wspolliniowe']:
        for p in diag['pary_wspolliniowe']:
            L.append(f'  r = {p["r"]:+.3f}   {p["a"]} <-> {p["b"]}')
    else:
        L.append('  Nie wykryto par o korelacji powyzej progu.')
    if diag['kryteria_zerowej_wariancji']:
        L.append(f'  Zerowa wariancja: '
                 f'{", ".join(diag["kryteria_zerowej_wariancji"])}')
    L.append(f'  Laczna waga kryteriow redundantnych: '
             f'{diag["laczna_waga_kryteriow_redundantnych"]} '
             f'z {diag["suma_wszystkich_wag"]}')
    L.append(f'  {diag["zalecenie"]}')
    return '\n'.join(L)


if __name__ == '__main__':
    # Dane demonstracyjne odwzorowujace typowy wynik BLAST dla 21-nt zapytan
    demo = [
        {'nazwa': 'kand_A', 'mismatches': 4, 'pct_identity': 81.0,
         'evalue': 2.5, 'align_length': 21, 'hit_count': 2,
         'bit_score': 24.3, 'gap_opens': 0},
        {'nazwa': 'kand_B', 'mismatches': 2, 'pct_identity': 90.5,
         'evalue': 0.008, 'align_length': 21, 'hit_count': 7,
         'bit_score': 34.2, 'gap_opens': 0},
        {'nazwa': 'kand_C', 'mismatches': 6, 'pct_identity': 71.4,
         'evalue': 45.0, 'align_length': 21, 'hit_count': 1,
         'bit_score': 18.1, 'gap_opens': 0},
        {'nazwa': 'kand_D', 'mismatches': 3, 'pct_identity': 85.7,
         'evalue': 0.9, 'align_length': 21, 'hit_count': 4,
         'bit_score': 28.0, 'gap_opens': 0},
    ]

    print('--- ZESTAW PELNY (7 kryteriow) ---')
    r_pelny = punktuj(demo)
    d = diagnostyka(demo)
    print(raport(r_pelny, d))

    print()
    print('--- ZESTAW ZREDUKOWANY (3 kryteria) ---')
    r_zred = punktuj(demo, zestaw=KRYTERIA_ZREDUKOWANE)
    for c in r_zred:
        print(f'  {c["ranga_blast"]}. {c["nazwa"]:<10} '
              f'wynik {c["wynik_blast"]:.4f}')

    print()
    print('--- POROWNANIE KOLEJNOSCI ---')
    kol_pelny = [c['nazwa'] for c in r_pelny]
    kol_zred = [c['nazwa'] for c in r_zred]
    print(f'  pelny       : {kol_pelny}')
    print(f'  zredukowany : {kol_zred}')
    print(f'  {"IDENTYCZNA" if kol_pelny == kol_zred else "ROZNI SIE"}')
