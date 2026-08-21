"""
=============================================================================
conservation.py  --  ANALIZA KONSERWATYWNOSCI
=============================================================================

PO CO

Przy projektowaniu siRNA przeciw wirusowi cel jest ruchomy. Populacja
wirusowa to zbior izolatow roznych sekwencyjnie, a siRNA zaprojektowane na
jeden izolat moze nie dzialac na pozostale. Dodatkowo wirus moze uciec przez
mutacje w miejscu docelowym.

Rozwiazaniem jest kierowanie siRNA w regiony, ktore sa niemal identyczne we
wszystkich izolatach - bo tam wirus ma niewiele mozliwosci mutacji bez utraty
funkcji.

=============================================================================
DWIE METODY - I DLACZEGO DOMYSLNA JEST TA DRUGA
=============================================================================

METODA A: ENTROPIA SHANNONA Z DOPASOWANIA WIELOSEKWENCYJNEGO

Klasyczne podejscie. Wykonuje sie MSA (np. MAFFT), a nastepnie dla kazdej
kolumny liczy entropie:

    H(x) = -suma p_i * log2(p_i)

gdzie p_i to czestosc nukleotydu i w danej kolumnie. H = 0 oznacza kolumne
calkowicie niezmienna.

Zalety: standard w literaturze, dziala przy insercjach i delecjach.
Wady: wymaga zewnetrznego programu do MSA, a przy kilkuset genomach jest
kosztowna obliczeniowo. Ponadto mierzy konserwatywnosc POZYCJI, a nie
calego okna.

METODA B: POKRYCIE DOKLADNE K-MEROW  [DOMYSLNA]

Dla kazdego okna dlugosci L w sekwencji referencyjnej sprawdzamy, w jakim
odsetku pozostalych sekwencji wystepuje DOKLADNIE ten sam ciag.

Uzasadnienie wyboru: siRNA wymaga niemal doskonalej komplementarnosci, zeby
zadzialac. Pojedyncze niedopasowanie w regionie seed wystarcza, by cel nie
zostal rozpoznany. Interesuje nas wiec nie to, czy pozycje sa "konserwatywne"
w sensie ewolucyjnym, tylko czy KONKRETNY 21-mer wystepuje w danym izolacie.

To jest dokladnie ta wielkosc, ktora chcemy zmaksymalizowac - i mierzy sie ja
bezposrednio, bez dopasowania.

Dodatkowe zalety:
  - nie wymaga MSA ani zadnego programu zewnetrznego
  - zlozonosc liniowa wzgledem rozmiaru zbioru
  - wynik ma bezposrednia interpretacje: "to siRNA pasuje do 94% izolatow"

Ograniczenie: nie rozroznia, ktore konkretnie pozycje sie roznia. Dlatego
modul liczy dodatkowo pokrycie z tolerancja jednego niedopasowania oraz -
jesli wejscie jest dopasowane - entropie Shannona.

=============================================================================
JAK UZYWAC
=============================================================================

    from conservation import analiza_konserwatywnosci

    wynik = analiza_konserwatywnosci(
        sciezka_fasta='tobrfv_310_genomow.fasta',
        dlugosci=(21, 22, 24))

    # nastepnie w design.generuj_kandydatow przekazac wynik jako filtr

Autor: Antonina Jarecka
=============================================================================
"""

from collections import Counter, defaultdict
from typing import Dict, Iterator, List, Optional, Tuple
import math


# ============================================================================
# WEJSCIE
# ============================================================================

def czytaj_fasta(sciezka: str) -> Iterator[Tuple[str, str]]:
    nazwa, buf = None, []
    with open(sciezka) as fh:
        for linia in fh:
            linia = linia.rstrip()
            if linia.startswith('>'):
                if nazwa is not None:
                    yield nazwa, ''.join(buf).upper()
                nazwa, buf = linia[1:].split()[0], []
            else:
                buf.append(linia)
    if nazwa is not None:
        yield nazwa, ''.join(buf).upper()


def wczytaj_zbior(sciezka: str, maks: Optional[int] = None
                  ) -> List[Tuple[str, str]]:
    """
    Wczytuje zbior sekwencji. Normalizuje U -> T.

    maks: opcjonalne ograniczenie liczby sekwencji (przydatne przy testach
          na duzych zbiorach).
    """
    zbior = []
    for i, (n, s) in enumerate(czytaj_fasta(sciezka)):
        if maks and i >= maks:
            break
        zbior.append((n, s.replace('U', 'T')))
    if not zbior:
        raise ValueError(f'Brak sekwencji w pliku {sciezka}')
    return zbior


# ============================================================================
# METODA B - POKRYCIE K-MEROW
# ============================================================================

def zbuduj_zbiory_kmerow(sekwencje: List[str], k: int) -> List[set]:
    """
    Dla kazdej sekwencji buduje zbior wszystkich jej k-merow.

    Pamiec: len(sekwencje) * dlugosc * (rozmiar stringa k). Dla 310 genomow
    ToBRFV (~6400 nt) i k=21 to okolo 2 mln k-merow, czyli rzad 200-300 MB.
    Przy wiekszych zbiorach uzyc funkcji pokrycie_strumieniowe().
    """
    return [{s[i:i + k] for i in range(len(s) - k + 1)} for s in sekwencje]


def pokrycie_okna(kmer: str, zbiory: List[set]) -> float:
    """Odsetek sekwencji zawierajacych dokladnie ten k-mer."""
    if not zbiory:
        return 0.0
    return sum(1 for z in zbiory if kmer in z) / len(zbiory)


def pokrycie_strumieniowe(referencja: str, sekwencje: List[str],
                          k: int) -> Dict[int, float]:
    """
    Wariant oszczedny pamieciowo - przetwarza sekwencje po kolei, nie
    trzymajac wszystkich zbiorow k-merow naraz.

    Zwraca slownik: pozycja_0based -> odsetek pokrycia.
    """
    okna = [referencja[i:i + k] for i in range(len(referencja) - k + 1)]
    licznik = Counter()
    n = 0
    for s in sekwencje:
        n += 1
        obecne = {s[i:i + k] for i in range(len(s) - k + 1)}
        for idx, w in enumerate(okna):
            if w in obecne:
                licznik[idx] += 1
    return {i: licznik[i] / n for i in range(len(okna))} if n else {}


def _warianty_1_niedopasowanie(kmer: str) -> Iterator[str]:
    """Wszystkie warianty k-meru z dokladnie jednym podstawieniem."""
    for i, orig in enumerate(kmer):
        for b in 'ACGT':
            if b != orig:
                yield kmer[:i] + b + kmer[i + 1:]


def pokrycie_z_tolerancja(kmer: str, zbiory: List[set]) -> float:
    """
    Odsetek sekwencji zawierajacych k-mer dokladnie ALBO z jednym
    niedopasowaniem.

    Interpretacja: siRNA z jednym niedopasowaniem poza regionem seed czesto
    nadal dziala, choc slabiej. Ta wielkosc daje gorne oszacowanie zasiegu.
    """
    if not zbiory:
        return 0.0
    warianty = set(_warianty_1_niedopasowanie(kmer))
    warianty.add(kmer)
    trafione = sum(1 for z in zbiory if z & warianty)
    return trafione / len(zbiory)


# ============================================================================
# METODA A - ENTROPIA SHANNONA
# ============================================================================

def entropia_kolumny(kolumna: str, ignoruj_przerwy: bool = True) -> float:
    """
    Entropia Shannona jednej kolumny dopasowania, w bitach.

    H = 0      -> kolumna calkowicie niezmienna
    H = 2      -> wszystkie cztery nukleotydy rownie czeste (maksimum)
    """
    znaki = [z for z in kolumna
             if z in 'ACGT' or (not ignoruj_przerwy and z == '-')]
    if not znaki:
        return 0.0
    n = len(znaki)
    licz = Counter(znaki)
    return -sum((c / n) * math.log2(c / n) for c in licz.values())


def entropia_pozycyjna(dopasowanie: List[str]) -> List[float]:
    """
    Entropia dla kazdej kolumny dopasowania.

    UWAGA: wejscie musi byc DOPASOWANE (wszystkie sekwencje tej samej
    dlugosci). Funkcja nie wykonuje dopasowania - do tego uzyj MAFFT:

        mafft --auto wejscie.fasta > dopasowane.fasta
    """
    dlugosci = {len(s) for s in dopasowanie}
    if len(dlugosci) != 1:
        raise ValueError(
            f'Sekwencje maja rozne dlugosci ({sorted(dlugosci)[:5]}...). '
            f'Entropia wymaga dopasowania wielosekwencyjnego. '
            f'Uzyj MAFFT albo metody pokrycia k-merow, ktora dopasowania '
            f'nie wymaga.')
    L = dlugosci.pop()
    return [entropia_kolumny(''.join(s[i] for s in dopasowanie))
            for i in range(L)]


def srednia_entropia_okna(entropie: List[float], start: int,
                          dlugosc: int) -> float:
    okno = entropie[start:start + dlugosc]
    return sum(okno) / len(okno) if okno else 0.0


# ============================================================================
# ANALIZA ZBIORCZA
# ============================================================================

def analiza_konserwatywnosci(sciezka_fasta: str,
                             dlugosci: Tuple[int, ...] = (21,),
                             indeks_referencji: int = 0,
                             prog_pokrycia: float = 0.95,
                             maks_sekwencji: Optional[int] = None,
                             licz_entropie: bool = False,
                             verbose: bool = True) -> Dict:
    """
    Pelna analiza konserwatywnosci zbioru sekwencji.

    Zwraca slownik z:
      - 'referencja'  : (nazwa, sekwencja) uzyta jako odniesienie
      - 'n_sekwencji' : rozmiar zbioru
      - 'pokrycie'    : {dlugosc: {pozycja_1based: odsetek}}
      - 'regiony'     : lista regionow spelniajacych prog, posortowana
      - 'entropia'    : lista entropii pozycyjnych (jesli licz_entropie)
    """
    zbior = wczytaj_zbior(sciezka_fasta, maks_sekwencji)
    nazwy = [n for n, _ in zbior]
    seq = [s for _, s in zbior]

    nazwa_ref, ref = zbior[indeks_referencji]
    pozostale = [s for i, s in enumerate(seq) if i != indeks_referencji]

    if verbose:
        print(f'Wczytano {len(zbior)} sekwencji')
        print(f'Referencja: {nazwa_ref} ({len(ref)} nt)')
        dl = [len(s) for s in seq]
        print(f'Zakres dlugosci: {min(dl)}-{max(dl)} nt')

    wynik = {'referencja': (nazwa_ref, ref), 'n_sekwencji': len(zbior),
             'nazwy': nazwy, 'pokrycie': {}, 'regiony': [], 'entropia': None}

    # --- entropia, tylko jesli wejscie dopasowane ---
    if licz_entropie:
        try:
            wynik['entropia'] = entropia_pozycyjna(seq)
            if verbose:
                sr = sum(wynik['entropia']) / len(wynik['entropia'])
                print(f'Srednia entropia pozycyjna: {sr:.4f} bitow')
        except ValueError as e:
            if verbose:
                print(f'Entropia pominieta: {e}')

    # --- pokrycie k-merow ---
    for L in dlugosci:
        zbiory = zbuduj_zbiory_kmerow(pozostale, L)
        pokr = {}
        for i in range(len(ref) - L + 1):
            kmer = ref[i:i + L]
            if set(kmer) - set('ACGT'):
                continue
            pokr[i + 1] = pokrycie_okna(kmer, zbiory)
        wynik['pokrycie'][L] = pokr

        spelnia = [(p, v) for p, v in pokr.items() if v >= prog_pokrycia]
        if verbose:
            print(f'  dlugosc {L} nt: {len(spelnia)} / {len(pokr)} okien '
                  f'z pokryciem >= {prog_pokrycia:.0%}')

        for p, v in spelnia:
            wynik['regiony'].append({
                'dlugosc': L, 'poz_start_1based': p,
                'poz_end_1based': p + L - 1,
                'sekwencja': ref[p - 1:p - 1 + L],
                'pokrycie': round(v, 4),
                'pokrycie_1mm': round(
                    pokrycie_z_tolerancja(ref[p - 1:p - 1 + L], zbiory), 4),
            })

    wynik['regiony'].sort(key=lambda r: (-r['pokrycie'],
                                         r['poz_start_1based']))
    return wynik


def scal_regiony(regiony: List[Dict], min_dlugosc: int = 30) -> List[Dict]:
    """
    Scala nakladajace sie okna w ciagle bloki konserwatywne.

    Przydatne do wskazania "zlotych regionow" - dluzszych odcinkow, w ktorych
    kazde okno spelnia prog. Takie bloki sa lepszym celem niz pojedyncze okna,
    bo daja swobode w doborze dokladnej pozycji siRNA.
    """
    if not regiony:
        return []
    wg_dl = defaultdict(list)
    for r in regiony:
        wg_dl[r['dlugosc']].append(r)

    bloki = []
    for L, lista in wg_dl.items():
        lista = sorted(lista, key=lambda r: r['poz_start_1based'])
        akt_start = lista[0]['poz_start_1based']
        akt_koniec = lista[0]['poz_end_1based']
        pokrycia = [lista[0]['pokrycie']]
        for r in lista[1:]:
            if r['poz_start_1based'] <= akt_koniec + 1:
                akt_koniec = max(akt_koniec, r['poz_end_1based'])
                pokrycia.append(r['pokrycie'])
            else:
                if akt_koniec - akt_start + 1 >= min_dlugosc:
                    bloki.append({'dlugosc_okna': L,
                                  'start': akt_start, 'koniec': akt_koniec,
                                  'dlugosc_bloku': akt_koniec - akt_start + 1,
                                  'srednie_pokrycie': round(
                                      sum(pokrycia) / len(pokrycia), 4)})
                akt_start, akt_koniec = r['poz_start_1based'], r['poz_end_1based']
                pokrycia = [r['pokrycie']]
        if akt_koniec - akt_start + 1 >= min_dlugosc:
            bloki.append({'dlugosc_okna': L, 'start': akt_start,
                          'koniec': akt_koniec,
                          'dlugosc_bloku': akt_koniec - akt_start + 1,
                          'srednie_pokrycie': round(
                              sum(pokrycia) / len(pokrycia), 4)})
    bloki.sort(key=lambda b: -b['dlugosc_bloku'])
    return bloki


def filtruj_kandydatow(kandydaci: List[Dict], wynik_kons: Dict,
                       prog: float = 0.95) -> List[Dict]:
    """
    Dodaje kandydatom pole 'pokrycie_izolatow' i odfiltrowuje te ponizej
    progu.

    Uzycie w pipeline: wywolac po design.generuj_kandydatow, przed scoringiem.
    """
    out = []
    for k in kandydaci:
        L = k['dlugosc']
        poz = k['poz_start_1based']
        pokr = wynik_kons['pokrycie'].get(L, {}).get(poz)
        if pokr is None:
            continue
        k = dict(k)
        k['pokrycie_izolatow'] = round(pokr, 4)
        if pokr >= prog:
            out.append(k)
    return out


if __name__ == '__main__':
    # Demonstracja na sztucznym zbiorze "izolatow"
    import random
    random.seed(11)

    baza = ''.join(random.choice('ACGT') for _ in range(300))
    # region 100-160 pozostawiamy niezmienny, reszta mutuje
    izolaty = [('ref', baza)]
    for i in range(20):
        s = list(baza)
        for p in range(len(s)):
            if 100 <= p < 160:
                continue
            if random.random() < 0.04:
                s[p] = random.choice('ACGT')
        izolaty.append((f'izolat_{i}', ''.join(s)))

    with open('/tmp/izolaty.fasta', 'w') as fh:
        for n, s in izolaty:
            fh.write(f'>{n}\n{s}\n')

    w = analiza_konserwatywnosci('/tmp/izolaty.fasta', dlugosci=(21,),
                                 prog_pokrycia=0.95)
    print(f'\nZnaleziono {len(w["regiony"])} okien konserwatywnych')
    print('\nTOP 5:')
    for r in w['regiony'][:5]:
        print(f'  poz {r["poz_start_1based"]:>4}-{r["poz_end_1based"]:<4} '
              f'pokrycie {r["pokrycie"]:.1%}  '
              f'(z 1 niedopasowaniem {r["pokrycie_1mm"]:.1%})')

    bloki = scal_regiony(w['regiony'], min_dlugosc=30)
    print(f'\nBloki ciagle (>=30 nt): {len(bloki)}')
    for b in bloki:
        print(f'  {b["start"]}-{b["koniec"]} ({b["dlugosc_bloku"]} nt), '
              f'srednie pokrycie {b["srednie_pokrycie"]:.1%}')
    print('\n(Oczekiwane: blok w okolicy 100-160, tam wymuszono brak mutacji)')


# ============================================================================
# AUTOMATED THRESHOLD OPTIMISATION
# ============================================================================
#
# THE PROBLEM
#
# The coverage threshold is a free parameter and every choice is a trade-off:
#
#   threshold too high  -> few or no windows survive; the design has nothing
#                          to work with
#   threshold too low   -> conserved regions are not distinguished from
#                          variable ones and the filter does no work
#
# Setting it by hand at 0.95 is conventional but arbitrary, and the value that
# makes sense depends on how variable the particular pathogen population is.
# For a clonal population almost every window exceeds 0.95; for a diverse one
# almost none does.
#
# THE APPROACH
#
# The number of windows passing the threshold, plotted against the threshold,
# is a monotonically decreasing curve. For a sequence set containing genuinely
# conserved regions it has a characteristic shape: a slow decline while only
# variable windows are being removed, then a sharp drop once the conserved
# core is reached.
#
# The point of maximum curvature - the knee - marks the transition. Above it,
# raising the threshold discards conserved windows; below it, lowering it
# admits variable ones. It is therefore the natural operating point.
#
# The knee is located by the Kneedle method (Satopaa et al. 2011): the curve
# is normalised to the unit square, the straight line between its endpoints is
# subtracted, and the point of maximum residual is taken.
#
# WHEN THE METHOD FAILS
#
# If the curve has no knee - if it declines uniformly - the sequence set has
# no distinct conserved core and the returned value is not meaningful. The
# function reports this explicitly rather than returning a number that looks
# authoritative. Diagnostics are provided so the decision can be inspected.

def _knee_point(x: List[float], y: List[float]) -> Tuple[int, float]:
    """
    Locates the knee of a decreasing curve by the Kneedle method.

    Returns (index, normalised_residual). A residual below roughly 0.05
    indicates that no clear knee exists.
    """
    n = len(x)
    if n < 4:
        return 0, 0.0

    x0, x1 = x[0], x[-1]
    y0, y1 = y[0], y[-1]
    dx = (x1 - x0) or 1.0
    dy = (y1 - y0) or 1.0
    xn = [(v - x0) / dx for v in x]
    yn = [(v - y0) / dy for v in y]

    # residual against the chord joining the endpoints
    reszty = [yn[i] - xn[i] for i in range(n)]
    idx = max(range(n), key=lambda i: abs(reszty[i]))
    return idx, abs(reszty[idx])


def optimise_threshold(sciezka_fasta: str,
                       dlugosc: int = 21,
                       zakres: Tuple[float, float] = (0.50, 1.00),
                       krok: float = 0.01,
                       min_okien: int = 10,
                       maks_sekwencji: Optional[int] = None,
                       verbose: bool = True) -> Dict:
    """
    Determines a coverage threshold from the data instead of assuming one.

    min_okien: the smallest number of surviving windows considered usable.
               A threshold that leaves fewer is rejected regardless of the
               curve shape, since it cannot support a design.

    Returns a dictionary containing:
        'threshold'     recommended value
        'method'        how it was obtained
        'confidence'    'high' | 'medium' | 'low'
        'curve'         [(threshold, n_windows), ...] for plotting
        'diagnostics'   text explaining the decision
    """
    zbior = wczytaj_zbior(sciezka_fasta, maks_sekwencji)
    nazwa_ref, ref = zbior[0]
    pozostale = [s for _, s in zbior[1:]]

    if verbose:
        print(f'Optimising threshold on {len(zbior)} sequences, '
              f'window {dlugosc} nt')

    zbiory = zbuduj_zbiory_kmerow(pozostale, dlugosc)
    pokrycia = []
    for i in range(len(ref) - dlugosc + 1):
        kmer = ref[i:i + dlugosc]
        if set(kmer) - set('ACGT'):
            continue
        pokrycia.append(pokrycie_okna(kmer, zbiory))

    if not pokrycia:
        return {'threshold': 0.95, 'method': 'default',
                'confidence': 'low', 'curve': [],
                'diagnostics': 'No valid windows in the reference sequence.'}

    progi = []
    p = zakres[0]
    while p <= zakres[1] + 1e-9:
        progi.append(round(p, 4))
        p += krok

    krzywa = [(t, sum(1 for v in pokrycia if v >= t)) for t in progi]
    xs = [t for t, _ in krzywa]
    ys = [float(n) for _, n in krzywa]

    idx, reszta = _knee_point(xs, ys)
    prog_knee = xs[idx]
    n_knee = int(ys[idx])

    diag = []
    diag.append(f'Windows evaluated: {len(pokrycia)}')
    diag.append(f'Coverage range: {min(pokrycia):.1%} to {max(pokrycia):.1%}')
    mediana = sorted(pokrycia)[len(pokrycia) // 2]
    diag.append(f'Median coverage: {mediana:.1%}')
    diag.append(f'Knee residual: {reszta:.3f}')

    # decision
    if reszta < 0.05:
        prog = 0.95
        metoda = 'default (no knee detected)'
        pewnosc = 'low'
        diag.append(
            'The curve declines uniformly, so no distinct conserved core is '
            'present. The population may be uniformly variable, or uniformly '
            'clonal. The conventional value of 0.95 is returned; inspect the '
            'curve before relying on it.')
    elif n_knee < min_okien:
        # step back until enough windows survive
        kandydaci = [t for t, n in krzywa if n >= min_okien]
        prog = max(kandydaci) if kandydaci else zakres[0]
        metoda = f'knee lowered to retain at least {min_okien} windows'
        pewnosc = 'medium'
        diag.append(
            f'The knee at {prog_knee:.2f} leaves only {n_knee} windows, '
            f'too few for a design. The threshold was lowered to {prog:.2f}.')
    else:
        prog = prog_knee
        metoda = 'knee of the coverage curve (Kneedle)'
        pewnosc = 'high'
        diag.append(
            f'A clear knee at {prog:.2f} retains {n_knee} windows. Above this '
            f'point conserved windows begin to be discarded, below it '
            f'variable ones are admitted.')

    n_final = sum(1 for v in pokrycia if v >= prog)
    diag.append(f'Windows at the recommended threshold: {n_final}')

    if verbose:
        for d in diag:
            print(f'  {d}')

    return {'threshold': round(prog, 3),
            'method': metoda,
            'confidence': pewnosc,
            'n_windows': n_final,
            'knee_threshold': round(prog_knee, 3),
            'knee_residual': round(reszta, 4),
            'curve': krzywa,
            'coverage_values': pokrycia,
            'diagnostics': diag,
            'reference': nazwa_ref}
