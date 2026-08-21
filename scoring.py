"""
=============================================================================
scoring.py  --  PUNKTACJA WAGOWA I RANKING
=============================================================================

STATUS: WLASNA IMPLEMENTACJA

=============================================================================
BLAD, KTORY TEN MODUL NAPRAWIA
=============================================================================

W dotychczasowym podejsciu wagi byly przypisywane BEZPOSREDNIO surowym
metrykom. To nie dziala, i warto dokladnie zrozumiec dlaczego, bo to typowy
blad, ktory recenzent wychwytuje natychmiast.

Rozwazmy dwie metryki:
    asymetria dG    -- zakres realny mniej wiecej od -6 do +6 kcal/mol
    zawartosc GC    -- zakres 30 do 52 procent

Jesli obu nadamy wage 1.0 i zsumujemy, to GC zdominuje wynik po prostu
dlatego, ze jego liczby sa okolo dziesieciokrotnie wieksze. Waga "1.0" przy
GC oznacza w praktyce dziesieciokrotnie wiekszy wplyw niz waga "1.0" przy
asymetrii. Czyli deklarowane wagi NIE ODPOWIADAJA rzeczywistemu wplywowi
metryki na ranking.

Konsekwencja praktyczna: mozesz sadzic, ze zbudowalas system wazacy piec
kryteriow, a faktycznie zbudowalas system sortujacy po zawartosci GC z
niewielkim szumem.

ROZWIAZANIE
-----------
Normalizacja min-max KAZDEJ metryki do przedzialu [0, 1] W OBREBIE
analizowanego zbioru kandydatow, ZANIM przylozymy wagi. Dopiero wtedy waga
0.3 naprawde oznacza trzykrotnie mniejszy wplyw niz waga 0.9.

Dodatkowo: metryki, w ktorych "mniej znaczy lepiej" (np. liczba trafien
off-target), odwracamy, zeby wszystkie szly w te sama strone.

=============================================================================
DRUGA RZECZ: SKAD BRAC WAGI
=============================================================================

Wagi ponizej sa USTALONE RECZNIE na podstawie literatury i sa najslabszym
ogniwem calego pipeline'u. Uczciwie trzeba to napisac wprost w metodyce.

Droga do poprawy, w kolejnosci rosnacego nakladu:

  1. ANALIZA WRAZLIWOSCI (tanio, zrob od razu)
     Przelicz ranking dla wag zmienionych o +/-25 procent i sprawdz, czy
     pierwsza piatka sie zmienia. Jesli nie -- ranking jest odporny i mozesz
     to zaraportowac jako argument. Jesli tak -- masz problem i lepiej
     dowiedziec sie o tym teraz niz od recenzenta.
     Funkcja: sensitivity_analysis()

  2. RANKING PARETO zamiast sumy wazonej (srednio drogo)
     Zamiast zgniatac piec wymiarow do jednej liczby, wyznacz front Pareto:
     kandydatow, ktorych nie da sie poprawic w jednym kryterium bez
     pogorszenia innego. Znika arbitralnosc wag. Wada: front bywa liczny.
     Funkcja: pareto_front()

  3. MODEL UCZONY (drogo, kierunek docelowy)
     Random Forest albo gradient boosting trenowany na publicznych danych
     o zmierzonej skutecznosci siRNA (Huesken i wsp. 2005 -- 2431 zmierzonych
     siRNA, to standardowy zbior referencyjny w tej dziedzinie).
     UWAGA: te dane sa SSACZE. Model wytrenowany na nich przeniesiony na
     rosliny bedzie obarczony bledem, ktory trzeba jawnie omowic. Zbioru
     roslinnego o porownywalnej wielkosci po prostu nie ma -- i to samo w
     sobie jest ciekawym watkiem w dyskusji.

Autor: Antonina Jarecka
=============================================================================
"""

from typing import Dict, List, Tuple, Optional
import copy

# --- Definicja metryk -------------------------------------------------------
# kierunek: +1 = wieksza wartosc lepsza, -1 = mniejsza wartosc lepsza
METRYKI = {
    'asymetria':        {'kierunek': +1, 'waga': 0.25,
                         'opis': 'asymetria termodynamiczna koncow dupleksu'},
    'dostepnosc_celu':  {'kierunek': +1, 'waga': 0.20,
                         'opis': 'MFE okna mRNA na nt (mniej ujemne = dostepniejsze)'},
    'mfe_guide':        {'kierunek': +1, 'waga': 0.10,
                         'opis': 'MFE samej nici prowadzacej (blizej 0 = mniej sfaldowana)'},
    'gc_optymalnosc':   {'kierunek': +1, 'waga': 0.10,
                         'opis': 'bliskosc GC do optimum 42 procent'},
    'koszt_offtarget':  {'kierunek': +1, 'waga': 0.25,
                         'opis': 'minimalny koszt dopasowania off-target (wyzszy = bezpieczniej)'},
    'seed_czystosc':    {'kierunek': -1, 'waga': 0.10,
                         'opis': 'liczba trafien seed 8mer w transkryptomie'},
}


def _minmax(values: List[float]) -> List[float]:
    """Normalizacja do [0,1]. Przy zerowej wariancji zwraca same 0.5."""
    lo, hi = min(values), max(values)
    if hi - lo < 1e-12:
        return [0.5] * len(values)
    return [(v - lo) / (hi - lo) for v in values]


def normalize_and_score(kandydaci: List[Dict],
                        wagi: Optional[Dict[str, float]] = None
                        ) -> List[Dict]:
    """
    Normalizuje metryki i liczy wynik laczny.

    kandydaci: lista slownikow, kazdy MUSI zawierac klucze z METRYKI.
    wagi:      opcjonalne nadpisanie wag domyslnych.

    Zwraca kopie listy z dodanymi polami:
        norm_<metryka>  -- wartosc znormalizowana [0,1]
        wynik           -- suma wazona [0,1]
        ranga           -- pozycja w rankingu (1 = najlepszy)
    """
    if not kandydaci:
        return []

    w = {k: v['waga'] for k, v in METRYKI.items()}
    if wagi:
        w.update(wagi)

    # walidacja: czy wagi sumuja sie do 1
    suma_wag = sum(w[k] for k in METRYKI)
    if abs(suma_wag - 1.0) > 1e-6:
        raise ValueError(
            f'Wagi sumuja sie do {suma_wag:.4f}, a powinny do 1.0. '
            f'Inaczej wyniki nie beda porownywalne miedzy uruchomieniami.')

    out = [copy.deepcopy(c) for c in kandydaci]

    for nazwa, spec in METRYKI.items():
        surowe = [float(c[nazwa]) for c in out]
        # odwracamy metryki typu "mniej znaczy lepiej"
        if spec['kierunek'] == -1:
            surowe = [-v for v in surowe]
        znorm = _minmax(surowe)
        for c, v in zip(out, znorm):
            c[f'norm_{nazwa}'] = round(v, 4)

    for c in out:
        c['wynik'] = round(
            sum(c[f'norm_{n}'] * w[n] for n in METRYKI), 4)

    out.sort(key=lambda c: -c['wynik'])
    for i, c in enumerate(out, 1):
        c['ranga'] = i
    return out


def sensitivity_analysis(kandydaci: List[Dict],
                         perturbacja: float = 0.25,
                         top_n: int = 5) -> Dict[str, object]:
    """
    Sprawdza, czy pierwsza piatka rankingu jest odporna na zmiane wag.

    Dla kazdej metryki po kolei: zwieksz jej wage o `perturbacja`, przeskaluj
    pozostale tak, zeby suma dalej wynosila 1, przelicz ranking, porownaj
    zbior top_n z rankingiem bazowym.

    Zwraca stabilnosc: udzial kandydatow, ktorzy zostaja w top_n we WSZYSTKICH
    wariantach. Wartosc bliska 1.0 = ranking odporny.
    """
    baza = normalize_and_score(kandydaci)
    top_baza = {c.get('nazwa', c.get('id', str(i)))
                for i, c in enumerate(baza[:top_n])}

    wyniki_wariantow = {}
    przeciecie = set(top_baza)

    for metryka in METRYKI:
        w = {k: v['waga'] for k, v in METRYKI.items()}
        w[metryka] *= (1.0 + perturbacja)
        # przeskaluj wszystkie, zeby suma = 1
        s = sum(w.values())
        w = {k: v / s for k, v in w.items()}

        wariant = normalize_and_score(kandydaci, wagi=w)
        top_w = {c.get('nazwa', c.get('id', str(i)))
                 for i, c in enumerate(wariant[:top_n])}
        wyniki_wariantow[metryka] = sorted(top_w)
        przeciecie &= top_w

    return {
        'top_bazowy': sorted(top_baza),
        'warianty': wyniki_wariantow,
        'stabilny_rdzen': sorted(przeciecie),
        'stabilnosc': round(len(przeciecie) / max(len(top_baza), 1), 3),
        'interpretacja': (
            'ranking odporny na wagi' if len(przeciecie) / max(len(top_baza), 1) >= 0.8
            else 'ranking wrazliwy na wagi -- rozwaz front Pareto'),
    }


def pareto_front(kandydaci: List[Dict],
                 metryki: Optional[List[str]] = None) -> List[Dict]:
    """
    Wyznacza front Pareto -- kandydatow niezdominowanych.

    Kandydat A dominuje B, jesli A jest niegorszy we WSZYSTKICH kryteriach
    i lepszy w PRZYNAJMNIEJ JEDNYM.

    Zaleta wzgledem sumy wazonej: nie wymaga wag w ogole.
    """
    if metryki is None:
        metryki = list(METRYKI.keys())

    znorm = normalize_and_score(kandydaci)
    front = []
    for a in znorm:
        zdominowany = False
        for b in znorm:
            if a is b:
                continue
            niegorszy = all(b[f'norm_{m}'] >= a[f'norm_{m}'] for m in metryki)
            lepszy = any(b[f'norm_{m}'] > a[f'norm_{m}'] for m in metryki)
            if niegorszy and lepszy:
                zdominowany = True
                break
        if not zdominowany:
            front.append(a)
    return front


def raport_tekstowy(ranking: List[Dict], top_n: int = 10) -> str:
    linie = []
    linie.append('RANKING KANDYDATOW')
    linie.append('=' * 100)
    naglowek = (f"{'#':<3} {'nazwa':<14} {'wynik':>7} " +
                ' '.join(f'{m[:9]:>10}' for m in METRYKI))
    linie.append(naglowek)
    linie.append('-' * 100)
    for c in ranking[:top_n]:
        nazwa = str(c.get('nazwa', c.get('id', '?')))[:14]
        wiersz = f"{c['ranga']:<3} {nazwa:<14} {c['wynik']:>7.4f} "
        wiersz += ' '.join(f"{c[f'norm_{m}']:>10.3f}" for m in METRYKI)
        linie.append(wiersz)
    linie.append('')
    linie.append('Wartosci w kolumnach metryk sa ZNORMALIZOWANE do [0,1] '
                 'w obrebie tego zbioru.')
    linie.append('Wagi: ' + ', '.join(
        f'{m}={METRYKI[m]["waga"]}' for m in METRYKI))
    return '\n'.join(linie)


if __name__ == '__main__':
    # Demo na sztucznych danych pokazujace problem ze skala
    demo = [
        {'nazwa': 'K1', 'asymetria': 4.8, 'dostepnosc_celu': -0.32,
         'mfe_guide': 0.0, 'gc_optymalnosc': 0.94, 'koszt_offtarget': 18.0,
         'seed_czystosc': 0},
        {'nazwa': 'K2', 'asymetria': 0.5, 'dostepnosc_celu': -0.28,
         'mfe_guide': -2.6, 'gc_optymalnosc': 0.99, 'koszt_offtarget': 21.0,
         'seed_czystosc': 0},
        {'nazwa': 'K3', 'asymetria': 4.2, 'dostepnosc_celu': -0.41,
         'mfe_guide': -0.1, 'gc_optymalnosc': 0.88, 'koszt_offtarget': 5.0,
         'seed_czystosc': 3},
    ]
    r = normalize_and_score(demo)
    print(raport_tekstowy(r))
    print()
    s = sensitivity_analysis(demo, top_n=2)
    print('ANALIZA WRAZLIWOSCI')
    print(f"  top bazowy      : {s['top_bazowy']}")
    print(f"  stabilny rdzen  : {s['stabilny_rdzen']}")
    print(f"  stabilnosc      : {s['stabilnosc']}")
    print(f"  interpretacja   : {s['interpretacja']}")
    print()
    print('FRONT PARETO:', [c['nazwa'] for c in pareto_front(demo)])
